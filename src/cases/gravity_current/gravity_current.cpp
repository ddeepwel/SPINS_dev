/* Script for the formation of a gravity current with zero initial velocity
 * and no topography */

/* ------------------ Top matter --------------------- */

// Required headers
#include "../../BaseCase.hpp"      // Support file containing default implementations of several functions
#include "../../Options.hpp"       // config-file parser
#include "../../Science.hpp"       // Science content
#include <random/normal.h>      // Blitz random number generator

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
// Grid types
DIMTYPE intype_x, intype_y, intype_z;

// Physical parameters
double g, rot_f, rho_0;         // gravity accel (m/s^2), Coriolis frequency (s^-1), reference density (kg/m^3)
double visco;                   // viscosity (m^2/s)
double mu;                      // dynamic viscosity (kg/(m·s))
double kappa_rho;               // diffusivity of density (m^2/s)
// helpful constants
const int Num_tracers = 1;      // number of tracers (density and dyes)
const int RHO = 0;              // index for rho

// Problem parameters
double delta_rho;               // density difference between different layers (% of reference density)
double delta_x;                 // horizontal transition length (m)
double Lmix;                    // Width of mixed region (m)

// Temporal parameters
double final_time;              // Final time (s)
double plot_interval;           // Time between field writes (s)

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
int iter = 0;                   // Iteration counter

// Maximum squared buoyancy frequency
double N2_max;

/* ------------------ Adjust the class --------------------- */

class userControl : public BaseCase {
    public:
        // Grid arrays
        Array<double,1> xx, yy, zz;

        // Timing variables (for outputs and measuring time steps)
        int plot_number;        // plot output number
        double next_plot;       // time of next output write
        double comp_duration;   // clock time since computation began
        double clock_time;      // current clock time

        /* Variables for Diagnostics */
        double max_u, max_v, max_w, max_vel, max_rho;
        double ke_x, ke_y, ke_z, ke_tot, pe_tot;

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

        // Coriolis parameter, viscosity, and diffusivities
        double get_visco() const { return visco; }
        double get_diffusivity(int t_num) const {
            return kappa_rho;
        }

        // Temporal parameters
        double init_time() const { return initial_time; }
        int get_restart_sequence() const { return restart_sequence; }

        // Number of tracers (the first is density)
        int numtracers() const { return Num_tracers; }

        /* Modify the timestep if necessary in order to land evenly on a plot time */
        double check_timestep (double intime, double now) {
            // Firstly, the buoyancy frequency provides a timescale that is not
            // accounted for with the velocity-based CFL condition.
            if (intime > 0.5/sqrt(N2_max)) {
                intime = 0.5/sqrt(N2_max);
            }
            // Now, calculate how many timesteps remain until the next writeout
            double until_plot = next_plot - now;
            int steps = ceil(until_plot / intime);
            // And calculate where we will actually be after (steps) timesteps
            // of the current size
            double true_fintime = steps*intime;

            // If that's close enough to the real writeout time, that's fine.
            if (fabs(until_plot - true_fintime) < 1e-6) {
                return intime;
            } else {
                // Otherwise, square up the timeteps.  This will always shrink the timestep.
                return (until_plot / steps);
            }
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
                v = 0;
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
                            if ( Ny > 1 )
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
            } else if (restarting and restart_from_dump) {
                init_tracer_dump("rho",rho);
            } else {
                // Density configuration
                rho = delta_rho*0.5*(1.0-tanh((xx(ii)-Lmix)/delta_x));
                // Write the array
                write_array(rho,"rho",plot_number);
            }
        }

        /* Forcing in the momentum equations */
        void forcing(double t, const DTArray & u, DTArray & u_f,
                const DTArray & v, DTArray & v_f, const DTArray & w, DTArray & w_f,
                vector<DTArray *> & tracers, vector<DTArray *> & tracers_f) {
            u_f = +rot_f*v;
            v_f = -rot_f*u;
            w_f = -g*(*tracers[RHO]);   // tracers[RHO] is defined as rho/rho_0
            *tracers_f[RHO] = 0;
        }

        void initialize_diagnostics_file() {
            if (master() and !restarting) {
                // create file for other diagnostics and write the column headers
                FILE * diagnos_file = fopen("diagnostics.txt","a");
                assert(diagnos_file);
                fprintf(diagnos_file,"Iter, Clock_time, Sim_time, "
                        "Max_U, Max_V, Max_W, Max_vel, "
                        "KE_x, KE_y, KE_z, Total_KE, Total_PE, "
                        "Max_density\n");
                fclose(diagnos_file);
            }
        }

        void write_diagnostics(double time) {
            if (master()) {
                /* add to the diagnostics file at each time step */
                FILE * diagnos_file = fopen("diagnostics.txt","a");
                assert(diagnos_file);
                fprintf(diagnos_file,"%d, %.12g, %.12f, "
                        "%.12g, %.12g, %.12g, %.12g, "
                        "%.12g, %.12g, %.12g, %.12g, %.12g, "
                        "%.12g\n",
                        iter,comp_duration,time,
                        max_u,max_v,max_w,max_vel,
                        ke_x,ke_y,ke_z,ke_tot,pe_tot,
                        max_rho);
                fclose(diagnos_file);
                /* and to the log file */
                fprintf(stdout,"[%d] (%.4g) %.4f: "
                        "%.4g %.4g %.4g "
                        "%.4g %.4g\n",
                        iter,comp_duration,time,
                        max_u,max_v,max_w,
                        ke_tot,max_rho);
            }
        }

        /* Basic analysis: compute secondary variables, and save fields and diagnostics */
        void analysis(double time, DTArray & u, DTArray & v, DTArray & w,
                vector<DTArray *> & tracers, DTArray & pressure) {
            // increase counter
            iter++;
            // Set-up
            if ( iter == 1 ) {
                // initialize the diagnostic files
                initialize_diagnostics_file();
                // Determine last plot if restarting from the dump case
                if (restart_from_dump) {
                    next_plot = (restart_sequence+1)*plot_interval;    
                }

            }
            // update clocks
            if (master()) {
                clock_time = MPI_Wtime();
                comp_duration = clock_time - compute_start_time;
            }

            /* Calculate and write out useful information */
            // Energy (PE assumes density is density anomaly)
            ke_x = pssum(sum(0.5*rho_0*(u*u)*
                        (*get_quad_x())(ii)*(*get_quad_y())(jj)*(*get_quad_z())(kk)));
            ke_y = pssum(sum(0.5*rho_0*(v*v)*
                        (*get_quad_x())(ii)*(*get_quad_y())(jj)*(*get_quad_z())(kk)));
            ke_z = pssum(sum(0.5*rho_0*(w*w)*
                        (*get_quad_x())(ii)*(*get_quad_y())(jj)*(*get_quad_z())(kk)));
            ke_tot = ke_x + ke_y + ke_z;
            pe_tot = pssum(sum(rho_0*(1+*tracers[RHO])*g*(zz(kk) - MinZ)*
                        (*get_quad_x())(ii)*(*get_quad_y())(jj)*(*get_quad_z())(kk)));
            // max of fields
            max_u = psmax(max(abs(u)));
            max_v = psmax(max(abs(v)));
            max_w = psmax(max(abs(w)));
            max_vel = psmax(max(sqrt(u*u + v*v + w*w)));
            max_rho = psmax(max(abs(*tracers[RHO])));

            // write to the diagnostic file
            write_diagnostics(time);

            /* Write to disk if at correct time */
            if ((time - next_plot) > -1e-6) {
                plot_number++;
                comp_duration = MPI_Wtime(); // time just before write (for dump)
                //Write the arrays
                write_array(u,"u",plot_number);
                write_array(w,"w",plot_number);
                if (Ny > 1 || rot_f != 0)
                    write_array(v,"v",plot_number);
                write_array(*tracers[RHO],"rho",plot_number);
                // update next plot time
                next_plot = next_plot + plot_interval;

                // Find average time to write (for dump)
                clock_time = MPI_Wtime(); // time just after write
                avg_write_time = (avg_write_time*(plot_number-restart_sequence-1) 
                        + (clock_time - comp_duration))/(plot_number-restart_sequence);
                // Print information about plot outputs
                if (master()) {
                    // in log file
                    fprintf(stdout,"*Write time: %.6g. Average write time: %.6g.\n",
                            clock_time - comp_duration, avg_write_time);
                    // track in a file
                    FILE * plottimes_file = fopen("plot_times.txt","a");
                    assert(plottimes_file);
                    if ( plot_number==restart_sequence+1 and !restarting )
                        fprintf(plottimes_file,"Output number, Simulation time (s), "
                                "Write time (s), Average write time (s)\n");
                    fprintf(plottimes_file,"%d, %.12f, %.12g, %.12g\n",
                                plot_number, time, clock_time - comp_duration, avg_write_time);
                    fclose(plottimes_file);
                }
            }

            // see if close to end of compute time and dump
            check_and_dump(clock_time, real_start_time, compute_time, time, avg_write_time,
                    plot_number, u, v, w, tracers);
            // Change dump log file if successfully reached final time
            successful_dump(plot_number, final_time, plot_interval);
        }

        // User specified variables to dump
        void write_variables(DTArray & u,DTArray & v, DTArray & w,
                vector<DTArray *> & tracers) {
            write_array(u,"u.dump");
            write_array(v,"v.dump");
            write_array(w,"w.dump");
            write_array(*tracers[RHO],"rho.dump");
        }

        // Constructor: Initialize local variables
        userControl():
            xx(split_range(Nx)), yy(Ny), zz(Nz),
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

    option_category("Physical parameters");
    add_option("g",&g,9.81,"Gravitational acceleration");
    add_option("rot_f",&rot_f,0.0,"Coriolis parameter");
    add_option("rho_0",&rho_0,1000.0,"Reference density");
    add_option("visco",&visco,"Viscosity");
    add_option("kappa_rho",&kappa_rho,"Diffusivity of density");

    option_category("Problem parameters");
    add_option("delta_rho",&delta_rho,"Density difference");
    add_option("delta_x",&delta_x,"Horizontal transition half-width");
    add_option("Lmix",&Lmix,"Width of collapse region");

    option_category("Temporal options");
    add_option("final_time",&final_time,"Final time");
    add_option("plot_interval",&plot_interval,"Time between writes");

    option_category("Restart options");
    add_option("restart",&restarting,false,"Restart from prior output time.");
    add_option("restart_time",&initial_time,0.0,"Time to restart from");
    add_option("restart_sequence",&restart_sequence,-1,"Sequence number to restart from");

    option_category("Dumping options");
    add_option("restart_from_dump",&restart_from_dump,false,"If restart from dump");
    add_option("compute_time",&compute_time,-1.0,"Time permitted for computation");

    option_category("Other options");
    add_option("perturb",&perturb,"Initial perturbation in velocity");

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

    // parse expansion types
    parse_boundary_conditions(xgrid_type, ygrid_type, zgrid_type, intype_x, intype_y, intype_z);

    // adjust Ly for 2D
    if (Ny==1 and Ly!=1.0) {
        Ly = 1.0;
        if (master())
            fprintf(stdout,"Simulation is 2 dimensional, "
                    "Ly has been changed to 1.0 for normalization.\n");
    }

    /* ------------------ Set correct initial time, and sequence --------------------- */

    if (restarting) {
        if (restart_sequence <= 0) {
            restart_sequence = int(initial_time/plot_interval);
        }
        if (master()) {
            fprintf(stdout,"Restart flags detected\n");
            fprintf(stdout,"Restarting from time %g, at sequence number %d\n",
                    initial_time,restart_sequence);
        }
    } else {
        // Not restarting, so set the initial sequence number
        // to the initial time / plot_interval
        restart_sequence = int(initial_time/plot_interval);
        if (fmod(initial_time,plot_interval) != 0.0) {
            if (master()) {
                fprintf(stdout,"Warning: the initial time (%g) does not appear "
                        "to be an even multiple of the plot interval (%g)\n",
                        initial_time,plot_interval);
            }
        }
    }

    /* ------------------ Derived parameters --------------------- */

    // Dynamic viscosity
    mu = visco*rho_0;
    // Maximum buoyancy frequency (squared) if the initial stratification was stable
    N2_max = g*delta_rho/(2*delta_x);

    /* ------------------ Print some parameters --------------------- */

    if (master()) {
        fprintf(stdout,"Dam break problem\n");
        fprintf(stdout,"Using a %f x %f x %f grid of %d x %d x %d points\n",Lx,Ly,Lz,Nx,Ny,Nz);
        fprintf(stdout,"g = %f, rot_f = %f, rho_0 = %f\n",g,rot_f,rho_0);
        fprintf(stdout,"Time between plots: %g s\n",plot_interval);
        fprintf(stdout,"Initial velocity perturbation: %g\n",perturb);
        fprintf(stdout,"Filter cutoff = %f, order = %f, strength = %f\n",f_cutoff,f_order,f_strength);
        fprintf(stdout,"Buoyancy frequency squared %g\n",N2_max);
    }

    /* ------------------ Do stuff --------------------- */

    // Create an instance of the above class
    userControl mycode;
    // Create a flow-evolver that takes its settings from the above class
    FluidEvolve<userControl> do_stuff(&mycode);
    // Initialize
    do_stuff.initialize();
    compute_start_time = MPI_Wtime(); // beginning of simulation (after reading in data)
    double startup_time = compute_start_time - real_start_time;
    if (master()) fprintf(stdout,"Start-up time: %.6g s.\n",startup_time);
    // Run until the end of time
    do_stuff.do_run(final_time);
    MPI_Finalize();
    return 0;
}
