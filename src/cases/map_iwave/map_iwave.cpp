/* Script for the generation of internal waves over a ridge with a
 * linear density stratification */

/* ------------------ Top matter --------------------- */

// Required headers
#include "../../BaseCase.hpp"      // Support file containing default implementations of several functions
#include "../../Options.hpp"       // config-file parser
#include <random/normal.h>         // Blitz random number generator

using namespace ranlib;

// Tensor variables for indexing
blitz::firstIndex ii;
blitz::secondIndex jj;
blitz::thirdIndex kk;

/* ------------------ Define parameters --------------------- */

// Grid scales
double Lx, Ly, Lz;              // Grid lengths (m)
int    Nx, Ny, Nz;              // Number of points in x, y, z
double MinX, MinY, MinZ;        // Minimum x/y/z points (m)
// Mapped grid?
bool mapped;
// Grid types
DIMTYPE intype_x, intype_y, intype_z;
string grid_type[3];

// Physical parameters
double g, rot_f, rho_0;         // gravity accel (m/s^2), Coriolis frequency (s^-1), reference density (kg/m^3)
double visco;                   // viscosity (m^2/s)
double mu;                      // dynamic viscosity (kg/(m·s))
double kappa_rho;               // diffusivity of density (m^2/s)
// helpful constants
const int Num_tracers = 1;      // number of tracers (density and dyes)
const int RHO = 0;              // index for density

// Problem parameters
double tide_strength;           // maximum horizontal velocity forced by the tide (m/s)
double tide_period;             // tidal period (s)
double tide_M2;                 // angular tidal frequency (rad/s)
double buoyancy_freq;           // buoyancy frequency (s^-1)

// Topography parameters
double hill_height;             // height of hill (m)
double hill_centre;             // position of hill peak (m)
double hill_width;              // width of hill (m)

// Temporal parameters
double final_time;              // Final time (s)
double plot_interval;           // Time between field writes (s)
double dt_max;                  // maximum time step (s)

// Restarting options
bool restarting;                // are you restarting?
double initial_time;            // initial start time of simulation
int restart_sequence;           // output number to restart from

// Dump parameters
bool restart_from_dump;         // restarting from dump?
double compute_time;            // requested computation time
double avg_write_time;          // average time to write all output fields at one output
double real_start_time;         // real (clock) time when simulation begins
double compute_start_time;      // real (clock) time when computation begins (after initialization)

// other options
double perturb;                 // Initial velocity perturbation
bool compute_enstrophy;         // Compute enstrophy?
bool compute_dissipation;       // Compute dissipation?
bool compute_BPE;               // Compute background potential energy?
bool compute_internal_to_BPE;   // Compute BPE gained from internal energy?
bool compute_stresses_top;      // Compute top surface stresses?
bool compute_stresses_bottom;   // Compute bottom surface stresses?
bool write_pressure;            // Write the pressure field?
int iter = 0;                   // Iteration counter

/* ------------------ Adjust the class --------------------- */

class userControl : public BaseCase {
    public:
        // Grid and topography arrays
        DTArray *zgrid;                 // Full z-grid field
        Array<double,1> xx, yy, zz;     // 1D grid vectors
        Array<double,1> topo;           // topography vector
        DTArray *Hprime;                // derivative of topography vector

        // Arrays and operators for derivatives
        Grad * gradient_op;
        DTArray *temp1, *dxdydz;

        // Timing variables (for outputs and measuring time steps)
        int plot_number;        // plot output number
        double next_plot;       // time of next output write
        double comp_duration;   // clock time since computation began
        double clock_time;      // current clock time

        // Size of domain
        double length_x() const { return Lx; }
        double length_y() const { return Ly; }
        double length_z() const { return Lz; }

        // Resolution in x, y, and z
        int size_x() const { return Nx; }
        int size_y() const { return Ny; }
        int size_z() const { return Nz; }

        // Set expansions (FREE_SLIP, NO_SLIP (in vertical) or PERIODIC)
        DIMTYPE type_x() const { return intype_x; }
        DIMTYPE type_y() const { return intype_y; }
        DIMTYPE type_z() const { return intype_z; }

        // Record the gradient-taking object
        void set_grad(Grad * in_grad) { gradient_op = in_grad; }

        // Coriolis parameter, viscosity, and diffusivities
        double get_rot_f() const { return rot_f; }
        double get_visco() const { return visco; }
        double get_diffusivity(int t_num) const {
            return kappa_rho;
        }

        // Temporal parameters
        double init_time() const { return initial_time; }
        int get_restart_sequence() const { return restart_sequence; }
        double get_dt_max() const { return dt_max; }
        double get_next_plot() { return next_plot; }

        // Number of tracers (the first is density)
        int numtracers() const { return Num_tracers; }

        // Create mapped grid
        bool is_mapped() const { return mapped; }
        void do_mapping(DTArray & xg, DTArray & yg, DTArray & zg) {
            zgrid = alloc_array(Nx,Ny,Nz);

            // over-write zz to be between -1 and 1
            // (zz defined in automatic_grid below)
            zz = -cos(ii*M_PI/(Nz-1));      // Chebyshev in vertical

            // Define topography
            topo = hill_height*exp(-pow((xx(ii)-hill_centre)/hill_width,2));

            // make full grids
            xg = xx(ii) + 0*jj + 0*kk;
            yg = yy(jj) + 0*ii + 0*kk;
            zg = MinZ + 0.5*Lz*(1+zz(kk)) + 0.5*(1-zz(kk))*topo(ii);
            *zgrid = zg;

            // Write the arrays and matlab readers
            write_array(xg,"xgrid");
            write_reader(xg,"xgrid",false);
            if (Ny > 1) {
                write_array(yg,"ygrid");
                write_reader(yg,"ygrid",false);
            }
            write_array(zg,"zgrid");
            write_reader(zg,"zgrid",false);
        }

        /* Initialize velocities */
        void init_vels(DTArray & u, DTArray & v, DTArray & w) {
            if (master()) fprintf(stdout,"Initializing velocities\n");
            // if restarting
            if (restarting and !restart_from_dump) {
                init_vels_restart(u, v, w);
            } else if (restarting and restart_from_dump) {
                init_vels_dump(u, v, w);
            } else{
                // else have a near motionless field
                u = 0;
                v = -tide_strength * rot_f / tide_M2 * Lz / (Lz - topo(ii));
                w = 0;
                // Add a random perturbation to trigger any 3D instabilities
                int myrank;
                MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
                Normal<double> rnd(0,1);
                for (int i = u.lbound(firstDim); i <= u.ubound(firstDim); i++) {
                    rnd.seed(i);
                    for (int j = u.lbound(secondDim); j <= u.ubound(secondDim); j++) {
                        for (int k = u.lbound(thirdDim); k <= u.ubound(thirdDim); k++) {
                            u(i,j,k) += perturb*rnd.random();
                            w(i,j,k) += perturb*rnd.random();
                            if (Ny > 1 || rot_f != 0)
                                v(i,j,k) += perturb*rnd.random();
                        }
                    }
                }
                // Write the arrays
                write_array(u,"u",plot_number);
                write_array(w,"w",plot_number);
                if (Ny > 1 || rot_f != 0) {
                    write_array(v,"v",plot_number);
                }
            }
        }

        /* Initialize the tracers (density and dyes) */
        void init_tracers(vector<DTArray *> & tracers) {
            if (master()) fprintf(stdout,"Initializing tracers\n");
            DTArray & rho = *tracers[RHO];
            
            if (restarting and !restart_from_dump) {
                init_tracer_restart("rho",rho);
                rho -= pow(buoyancy_freq,2) / g * (*zgrid)(ii,jj,kk);
            } else if (restarting and restart_from_dump) {
                init_tracer_dump("rho",rho);
                rho -= pow(buoyancy_freq,2) / g * (*zgrid)(ii,jj,kk);
            } else {
                // Write the array of the perturbation density (no perturbation)
                rho = 0*ii + 0*jj + 0*kk;
                write_array(rho,"rho",plot_number);

                // density configuration (rho = rho'/rho_0)
                rho = -pow(buoyancy_freq,2) / g * (*zgrid)(ii,jj,kk);
                // Important! if mapped, and rho depends on z
                // then (*zgrid)(ii,jj,kk), must be used in stead of zz(kk)
            }
        }

        /* Forcing in the momentum equations */
        void forcing(double t, const DTArray & u, DTArray & u_f,
                const DTArray & v, DTArray & v_f, const DTArray & w, DTArray & w_f,
                vector<DTArray *> & tracers, vector<DTArray *> & tracers_f) {
            // southern hemisphere rotation
            u_f = -rot_f*v + tide_strength * cos(tide_M2 * t) *
                (tide_M2 - rot_f*rot_f/tide_M2 * Lz/(Lz-topo(ii)));
            v_f = +rot_f*u;
            w_f = -g*(*tracers[RHO]); // tracers[RHO] is rho'/rho_0
            *tracers_f[RHO] = 0;
        }

        /* Basic analysis: compute secondary variables, and save fields and diagnostics */
        void analysis(double time, DTArray & u, DTArray & v, DTArray & w,
                vector<DTArray *> & tracers, DTArray & pressure) {
            // Set-up
            if ( iter == 0 ) {
                if ( compute_enstrophy or compute_dissipation or
                        compute_stresses_top or compute_stresses_bottom ) {
                    temp1 = alloc_array(Nx,Ny,Nz);
                }
                if ( compute_stresses_top or compute_stresses_bottom ) {
                    // initialize the vector of the bottom slope (Hprime)
                    Hprime = alloc_array(Nx,Ny,1);
                    if (is_mapped()) {
                        bottom_slope(*Hprime, *zgrid, *temp1, gradient_op, grid_type, Nx, Ny, Nz);
                    } else {
                        topo = 0*ii;
                        *Hprime = 0*ii + 0*jj;
                    }
                }
                // Determine last plot if restarting from the dump file
                if (restart_from_dump) {
                    next_plot = (restart_sequence+1)*plot_interval;
                }
                // initialize the size of each voxel
                dxdydz = alloc_array(Nx,Ny,Nz);
                *dxdydz = (*get_quad_x())(ii)*(*get_quad_y())(jj)*(*get_quad_z())(kk);
                if (is_mapped()) {
                    *dxdydz = (*dxdydz)*(Lz-topo(ii))/Lz;
                }
            }
            // update clocks
            if (master()) {
                clock_time = MPI_Wtime();
                comp_duration = clock_time - compute_start_time;
            }

            /* Calculate and write out useful information */

            // Energy (PE assumes density is density anomaly)
            double ke_x = 0, ke_y = 0, ke_z = 0;
            if ( Nx > 1 ) {
                ke_x = pssum(sum(0.5*rho_0*(u*u)*(*dxdydz)));
            }
            if (Ny > 1 || rot_f != 0) {
                ke_y = pssum(sum(0.5*rho_0*(v*v)*(*dxdydz)));
            }
            if ( Nz > 1 ) {
                ke_z = pssum(sum(0.5*rho_0*(w*w)*(*dxdydz)));
            }
            double pe_tot;
            if (is_mapped()) {
                pe_tot = pssum(sum(rho_0*(1+(*tracers[RHO]))*g*((*zgrid)(ii,jj,kk) - MinZ)*(*dxdydz)));
            } else {
                pe_tot = pssum(sum(rho_0*(1+(*tracers[RHO]))*g*(zz(kk) - MinZ)*(*dxdydz)));
            }
            double BPE_tot = 0;
            if (compute_BPE) {
                compute_Background_PE(BPE_tot, *tracers[RHO], *dxdydz, Nx, Ny, Nz, Lx, Ly, Lz,
                        g, rho_0, iter, false, is_mapped(), topo);
            }
            // Conversion from internal energy to background potential energy
            double phi_i = 0;
            if (compute_internal_to_BPE) {
                compute_BPE_from_internal(phi_i, *tracers[RHO], kappa_rho, rho_0, g, Nz);
            }
            // viscous dissipation
            double diss_tot = 0;
            double max_diss = 0;
            if (compute_dissipation) {
                dissipation(*temp1, u, v, w, gradient_op, grid_type, Nx, Ny, Nz, mu);
                max_diss = psmax(max(*temp1));
                diss_tot = pssum(sum((*temp1)*(*dxdydz)));
            }
            // Vorticity / Enstrophy
            double max_vort_x = 0, enst_x_tot = 0;
            double max_vort_y = 0, enst_y_tot = 0;
            double max_vort_z = 0, enst_z_tot = 0;
            if (compute_enstrophy) {
                // x-vorticity
                if (Ny > 1 and Nz > 1) {
                    compute_vort_x(*temp1, v, w, gradient_op, grid_type);
                    max_vort_x = psmax(max(abs(*temp1)));
                    enst_x_tot = pssum(sum(0.5*pow(*temp1,2)*(*dxdydz)));
                }
                // y-vorticity
                if (Nx > 1 and Nz > 1) {
                    compute_vort_y(*temp1, u, w, gradient_op, grid_type);
                    max_vort_y = psmax(max(abs(*temp1)));
                    enst_y_tot = pssum(sum(0.5*pow(*temp1,2)*(*dxdydz)));
                }
                // z-vorticity
                if (Nx > 1 and Ny > 1) {
                    compute_vort_z(*temp1, u, v, gradient_op, grid_type);
                    max_vort_z = psmax(max(abs(*temp1)));
                    enst_z_tot = pssum(sum(0.5*pow(*temp1,2)*(*dxdydz)));
                }
            }
            // max of fields
            double max_u = psmax(max(abs(u)));
            double max_v = psmax(max(abs(v)));
            double max_w = psmax(max(abs(w)));
            double max_vel = psmax(max(sqrt(u*u + v*v + w*w)));
            double max_rho = psmax(max(*tracers[RHO]));
            double min_rho = psmin(min(*tracers[RHO]));
            // total mass (tracers[RHO] is non-dimensional density)
            double mass = pssum(sum(rho_0*(1+(*tracers[RHO]))*(*dxdydz)));
            // energy added/removed by tidal body force
            double E_tide = pssum(sum(
                        rho_0 * u * tide_strength * cos(tide_M2 * time) *
                        (tide_M2 - rot_f*rot_f/tide_M2 * Lz/(Lz-topo(ii))) * (*dxdydz)
                        ));

            if (master()) {
                // add diagnostics to buffers
                string header, line;
                add_diagnostic("Iter", iter,            header, line);
                add_diagnostic("Clock_time", comp_duration, header, line);
                add_diagnostic("Time", time,            header, line);
                add_diagnostic("Max_vel", max_vel,      header, line);
                add_diagnostic("Max_density", max_rho,  header, line);
                add_diagnostic("Min_density", min_rho,  header, line);
                add_diagnostic("Mass", mass,            header, line);
                add_diagnostic("PE_tot", pe_tot,        header, line);
                add_diagnostic("KE_from_forcing", E_tide,   header, line);
                if (compute_BPE) {
                    add_diagnostic("BPE_tot", BPE_tot,  header, line);
                }
                if (compute_internal_to_BPE) {
                    add_diagnostic("BPE_from_int", phi_i,   header, line);
                }
                if (compute_dissipation) {
                    add_diagnostic("Max_diss", max_diss,    header, line);
                    add_diagnostic("Diss_tot", diss_tot,    header, line);
                }
                if (Nx > 1) {
                    add_diagnostic("Max_u", max_u,  header, line);
                    add_diagnostic("KE_x", ke_x,    header, line);
                }
                if (Ny > 1 || rot_f != 0) {
                    add_diagnostic("Max_v", max_v,  header, line);
                    add_diagnostic("KE_y", ke_y,    header, line);
                }
                if (Nz > 1) {
                    add_diagnostic("Max_w", max_w,  header, line);
                    add_diagnostic("KE_z", ke_z,    header, line);
                }
                if (Ny > 1 && Nz > 1 && compute_enstrophy) {
                    add_diagnostic("Enst_x_tot", enst_x_tot, header, line);
                    add_diagnostic("Max_vort_x", max_vort_x, header, line);
                }
                if (Nx > 1 && Nz > 1 && compute_enstrophy) {
                    add_diagnostic("Enst_y_tot", enst_y_tot, header, line);
                    add_diagnostic("Max_vort_y", max_vort_y, header, line);
                }
                if (Nx > 1 && Ny > 1 && compute_enstrophy) {
                    add_diagnostic("Enst_z_tot", enst_z_tot, header, line);
                    add_diagnostic("Max_vort_z", max_vort_z, header, line);
                }

                // Write to file
                if (!(restarting and iter==0))
                    write_diagnostics(header, line, iter, restarting);
                // and to the log file
                fprintf(stdout,"[%d] (%.4g) %.4f: "
                        "%.4g %.4g %.4g %.4g\n",
                        iter,comp_duration,time,
                        max_u,max_v,max_w,max_rho);
            }

            // Top Surface Stresses
            if ( compute_stresses_top ) {
                stresses_top(u, v, w, *Hprime, *temp1, gradient_op, grid_type, mu, time, iter, restarting);
            }
            // Bottom Surface Stresses
            if ( compute_stresses_bottom ) {
                stresses_bottom(u, v, w, *Hprime, *temp1, gradient_op, grid_type, mu, time, iter, restarting);
            }

            /* Write to disk if at correct time */
            if ((time - next_plot) > -1e-6) {
                plot_number++;
                comp_duration = MPI_Wtime(); // time just before write (for dump)
                // Write the arrays
                write_array(u,"u",plot_number);
                write_array(w,"w",plot_number);
                if (Ny > 1 || rot_f != 0)
                    write_array(v,"v",plot_number);
                *temp1 = (*tracers[RHO]) + pow(buoyancy_freq,2) / g * (*zgrid)(ii,jj,kk);
                write_array(*temp1,"rho",plot_number);
                if (write_pressure)
                    write_array(pressure,"p",plot_number);
                // update next plot time
                next_plot = next_plot + plot_interval;

                // Find average time to write (for dump)
                clock_time = MPI_Wtime(); // time just after write
                avg_write_time = (avg_write_time*(plot_number-restart_sequence-1) 
                        + (clock_time - comp_duration))/(plot_number-restart_sequence);
                // Print information about plot outputs
                write_plot_times(time, clock_time, comp_duration, avg_write_time, plot_number, restarting);
            }

            // see if close to end of compute time and dump
            check_and_dump(clock_time, real_start_time, compute_time, time, avg_write_time,
                    plot_number, iter, u, v, w, tracers);
            // Change dump log file if successfully reached final time
            successful_dump(plot_number, final_time, plot_interval);
            // increase counter
            iter++;
        }

        // User specified variables to dump
        void write_variables(DTArray & u,DTArray & v, DTArray & w,
                vector<DTArray *> & tracers) {
            write_array(u,"u.dump");
            write_array(v,"v.dump");
            write_array(w,"w.dump");
            *temp1 = (*tracers[RHO]) + pow(buoyancy_freq,2) / g * (*zgrid)(ii,jj,kk);
            write_array(*temp1,"rho.dump");
        }

        // Constructor: Initialize local variables
        userControl():
            xx(split_range(Nx)), yy(Ny), zz(Nz),
            topo(split_range(Nx)), gradient_op(0),
            plot_number(restart_sequence),
            next_plot(initial_time + plot_interval)
    {   compute_quadweights(
            size_x(),   size_y(),   size_z(),
            length_x(), length_y(), length_z(),
            type_x(),   type_y(),   type_z());
        // Create one-dimensional arrays for the coordinates
        automatic_grid(MinX, MinY, MinZ, &xx, &yy, &zz);
    }
};

/* The ''main'' routine */
int main(int argc, char ** argv) {
    /* Initialize MPI.  This is required even for single-processor runs,
       since the inner routines assume some degree of parallelization,
       even if it is trivial. */
    MPI_Init(&argc, &argv);

    real_start_time = MPI_Wtime();     // start of simulation (for dump)
    /* ------------------ Define parameters from spins.conf --------------------- */
    options_init();

    option_category("Grid Options");
    add_option("Lx",&Lx,"Length of tank");
    add_option("Ly",&Ly,1.0,"Width of tank");
    add_option("Lz",&Lz,"Height of tank");
    add_option("Nx",&Nx,"Number of points in X");
    add_option("Ny",&Ny,1,"Number of points in Y");
    add_option("Nz",&Nz,"Number of points in Z");
    add_option("min_x",&MinX,0.0,"Minimum X-value");
    add_option("min_y",&MinY,0.0,"Minimum Y-value");
    add_option("min_z",&MinZ,0.0,"Minimum Z-value");

    option_category("Grid expansion options");
    string xgrid_type, ygrid_type, zgrid_type;
    add_option("type_x",&xgrid_type,
            "Grid type in X.  Valid values are:\n"
            "   FOURIER: Periodic\n"
            "   FREE_SLIP: Cosine expansion\n"
            "   NO_SLIP: Chebyhsev expansion");
    add_option("type_y",&ygrid_type,"FOURIER","Grid type in Y");
    add_option("type_z",&zgrid_type,"Grid type in Z");
    add_option("mapped_grid",&mapped,true,"Is the grid mapped?");

    option_category("Topography parameters");
    add_option("hill_height",&hill_height,"Height of hill");
    add_option("hill_centre",&hill_centre,"location of hill peak");
    add_option("hill_width",&hill_width,"Width of hill");

    option_category("Physical parameters");
    add_option("g",&g,9.81,"Gravitational acceleration");
    add_option("rot_f",&rot_f,0.0,"Coriolis parameter");
    add_option("rho_0",&rho_0,1000.0,"Reference density");
    add_option("visco",&visco,"Viscosity");
    add_option("kappa_rho",&kappa_rho,"Diffusivity of density");

    option_category("Problem parameters");
    add_option("tide_strength",&tide_strength,"maximum tidal velocity");
    add_option("tide_period",&tide_period,"tidal period");
    add_option("buoyancy_freq",&buoyancy_freq,"buoyancy frequency");

    option_category("Temporal options");
    add_option("final_time",&final_time,"Final time");
    add_option("plot_interval",&plot_interval,"Time between writes");
    add_option("dt_max",&dt_max,0.0,"Maximum time step. Zero value results in the default");

    option_category("Restart options");
    add_option("restart",&restarting,false,"Restart from prior output time.");
    add_option("restart_time",&initial_time,0.0,"Time to restart from");
    add_option("restart_sequence",&restart_sequence,-1,"Sequence number to restart from");

    option_category("Dumping options");
    add_option("restart_from_dump",&restart_from_dump,false,"If restart from dump");
    add_option("compute_time",&compute_time,-1.0,"Time permitted for computation");

    option_category("Other options");
    add_option("perturb",&perturb,"Initial perturbation in velocity");
    add_option("compute_enstrophy",&compute_enstrophy,true,"Calculate enstrophy?");
    add_option("compute_dissipation",&compute_dissipation,true,"Calculate dissipation?");
    add_option("compute_BPE",&compute_BPE,true,"Calculate BPE?");
    add_option("compute_internal_to_BPE",&compute_internal_to_BPE,true,
            "Calculate BPE gained from internal energy?");
    add_option("compute_stresses_top",&compute_stresses_top,false,"Calculate top surfaces stresses?");
    add_option("compute_stresses_bottom",&compute_stresses_bottom,false,"Calculate bottom surfaces stresses?");
    add_option("write_pressure",&write_pressure,false,"Write the pressure field?");

    option_category("Filter options");
    add_option("f_cutoff",&f_cutoff,0.6,"Filter cut-off frequency");
    add_option("f_order",&f_order,2.0,"Filter order");
    add_option("f_strength",&f_strength,20.0,"Filter strength");

    // Parse the options from the command line and config file
    options_parse(argc,argv);

    /* ------------------ Adjust and check parameters --------------------- */
    /* Now, make sense of the options received.  Many of these
     * can be directly used, but the ones of string-type need further procesing. */

    // adjust temporal values when restarting from dump
    if (restart_from_dump) {
        adjust_for_dump(restarting, initial_time, restart_sequence,
                final_time, compute_time, avg_write_time, Num_tracers, Nx, Ny, Nz);
    }

    // check restart sequence
    check_restart_sequence(restarting, restart_sequence, initial_time, plot_interval);

    // parse expansion types
    parse_boundary_conditions(xgrid_type, ygrid_type, zgrid_type, intype_x, intype_y, intype_z);
    // vector of expansion types
    grid_type[0] = xgrid_type;
    grid_type[1] = ygrid_type;
    grid_type[2] = zgrid_type;

    // adjust Ly for 2D
    if (Ny==1 and Ly!=1.0) {
        Ly = 1.0;
        if (master())
            fprintf(stdout,"Simulation is 2 dimensional, "
                    "Ly has been changed to 1.0 for normalization.\n");
    }

    /* ------------------ Derived parameters --------------------- */

    // Angular tidal frequency
    tide_M2 = 2*M_PI/tide_period;

    // Dynamic viscosity
    mu = visco*rho_0;
    // Maximum time step
    if (dt_max == 0.0) {
        // if dt_max not given in spins.conf, use the buoyancy frequency
        dt_max = 0.5/buoyancy_freq;
    }

    /* ------------------ Initialize --------------------- */

    // Create an instance of the above class
    userControl mycode;
    // Create a flow-evolver that takes its settings from the above class
    FluidEvolve<userControl> do_stuff(&mycode);
    // Initialize
    do_stuff.initialize();

    /* ------------------ Print some parameters --------------------- */

    compute_start_time = MPI_Wtime(); // beginning of simulation (after reading in data)
    double startup_time = compute_start_time - real_start_time;

    if (master()) {
        fprintf(stdout,"Internal wave generation problem\n");
        fprintf(stdout,"Using a %f x %f x %f grid of %d x %d x %d points\n",Lx,Ly,Lz,Nx,Ny,Nz);
        fprintf(stdout,"g = %f, rot_f = %f, rho_0 = %f\n",g,rot_f,rho_0);
        fprintf(stdout,"Time between plots: %g s\n",plot_interval);
        fprintf(stdout,"Initial velocity perturbation: %g\n",perturb);
        fprintf(stdout,"Filter cutoff = %f, order = %f, strength = %f\n",f_cutoff,f_order,f_strength);
        fprintf(stdout,"Approx. max buoyancy frequency squared: %g\n",pow(buoyancy_freq,2));
        fprintf(stdout,"Max time step: %g\n",dt_max);
        fprintf(stdout,"Start-up time: %.6g s.\n",startup_time);
    }

    /* ------------------ Run --------------------- */
    // Run until the end of time
    do_stuff.do_run(final_time);
    MPI_Finalize();
    return 0;
}
