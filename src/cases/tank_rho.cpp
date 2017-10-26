/* cases/tank_rho.cpp */
/* Generic script for two layer fluid
   with zero initial velocity
   and with topography */

// Required headers
#include "../TArray.hpp"        // Custom extensions to the library to support FFTs
#include "../NSIntegrator.hpp"  // Time-integrator for the Navier-Stokes equations
#include "../BaseCase.hpp"      // Support file that contains default implementations of several functions
#include "../Options.hpp"       // config-file parser
#include "../Science.hpp"       // Additional analysis routines
#include "../Par_util.hpp"
#include "../T_util.hpp"
#include <math.h>
#include <random/normal.h>      // Blitz random number generator
#include <blitz/array.h>        // Blitz++ array library
#include <mpi.h>                // MPI parallel library
#include <vector>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "../timing.hpp"
//#include <stdlib.h>

using namespace std;
using namespace TArrayn;
using namespace NSIntegrator;
using namespace ranlib;

// Tensor variables for indexing
blitz::firstIndex ii;
blitz::secondIndex jj;
blitz::thirdIndex kk;

/* ------------------ Parameters --------------------- */
// Grid scales
double Lx, Ly, Lz;              // (m)
int    Nx, Ny, Nz;              // Points in x, y (span), z directions
double MinX, MinY, MinZ;        // Minimum x/y/z points

// Physical constants
double g, rot_f;                // gravity accel (m/s^2), Coriolis frequency (s^-1)

// Stratification parameters
double rho_0;                   // reference density (kg/L)
double delta_rho;               // change in density (kg/L)
// pycnocline location parameters
double pyc_asym;                // % of depth to shift pycnocline above the mid-depth
double pyc_sep_perc;            // total separation of double pycnocline (as % of depth)
double h_perc;                  // pycnocline half-width (as % of depth)
double h_mix_perc;              // vertical half-width transition of mixed region (as % of depth)
// Horizontal stratification parameters
double delta_x;                 // horizontal transition length (m)
double Lmix;                    // Width of mixed region (m)
double Hmix;                    // Total height of mixed region (as % of depth)

// Hill parameters
double hill_height;             // height of hill (as % of depth)
double hill_centre;             // position of hill peak
double hill_width;              // width of hill

// Viscosity and diffusivity of density and tracers
static const int NUM_TRACER = 1;
static const int RHO = 0;
double VISCO;
double DIFFU_rho;
double DIFFU[ NUM_TRACER ];

// Temporal parameters
double plot_interval;           // Time between field writes (s)
double final_time;              // Final time (s)

// Vertical chain parameters
bool savechain;                 // (boolean) Flag to use chain or not
double chain1_start_time;       // time to start saving chain (s)
double chain1_end_time;         // time to stop saving chain (s)
double chain_plot_interval;     // time between chain writes (s)
// Chain locations (defined in main program)
int chain1_xindex, chain1_yindex;

// Initial velocity perturbation
double perturb;

// Dump parameters
bool restart_from_dump;
double real_start_time;
double compute_time;
double total_run_time;
double avg_write_time;

// Restarting options
bool restarting = false;
double restart_time;
double initial_time = 0;
int restart_sequence;

// Iteration counter
int itercount = 0;

/* ------------------ Derived parameters --------------------- */

// Pycnocline half-width
double h_halfwidth;
double h_mix_half;

// Flow speed
double c0;

// Squared maximum buoyancy frequency if the initial stratification was stable
double N2_max;

// Reynolds number
double Re;

/* ------------------ Initialize the class --------------------- */

class mapiw : public BaseCase {
    public:
        /* Grid arrays */
        DTArray * xgrid, * ygrid, * zgrid;
        Array<double,1> hill;

        /* Timing variables (for outputs and measuring time steps) */
        int plot_number;
        double last_plot, chain_last_plot;
        // variables for timing steps
        double t_step;
        double clock_time, step_start_time;

        /* Size of domain */
        double length_x() const { return Lx;}
        double length_y() const { return Ly;}
        double length_z() const { return Lz;}

        /* Resolution in X, Y, and Z */
        int size_x() const { return Nx; }
        int size_y() const { return Ny; }
        int size_z() const { return Nz; }

        /* Set expansions (FREE_SLIP, NO_SLIP (in vertical) or PERIODIC) */
        DIMTYPE type_x() const { return FREE_SLIP; }
        DIMTYPE type_y() const { return FREE_SLIP; }
        DIMTYPE type_z() const { return NO_SLIP; }
        DIMTYPE type_default() const { return PERIODIC; }

        /* Number of tracers */
        int numActive() const { return NUM_TRACER; }

        /* Viscosity and diffusivity */
        double get_visco() const { return VISCO; }
        double get_diffusivity(int t_num) const {
            switch (t_num) {
                case RHO:
                    return DIFFU[ RHO ];
                default:
                    abort();
            }
        }

        /* Create mapped grid */
        bool is_mapped() const {return true;}
        void do_mapping(DTArray & xg, DTArray & yg, DTArray & zg) {
            xgrid = alloc_array(Nx,Ny,Nz);
            ygrid = alloc_array(Nx,Ny,Nz);
            zgrid = alloc_array(Nx,Ny,Nz);

            Array<double,1> xx(split_range(Nx)), yy(Ny),zz(Nz);
            // Use periodic coordinates in horizontal
            xx = MinX + Lx*(ii+0.5)/Nx;     // x-coordinate
            yy = MinY + Ly*(ii+0.5)/Ny;     // y-coordinate
            zz = cos(ii*M_PI/(Nz-1));       // Chebyshev in vertical

            xg = xx(ii) + 0*jj + 0*kk;
            *xgrid = xg;

            yg = yy(jj) + 0*ii + 0*kk;
            *ygrid = yg;

            // a Gaussian hill
            hill = hill_height*Lz*exp(-pow((xx(ii)-hill_centre)/hill_width,2));
            zg = MinZ + 0.5*Lz*(1+zz(kk)) + 0.5*(1-zz(kk))*hill(ii);
            *zgrid = zg;

            //write_array(xg,"xgrid");
            //write_reader(xg,"xgrid",false);

            /*if (Ny > 1 || rot_f != 0) {
                write_array(yg,"ygrid");
                write_reader(yg,"ygrid",false);
            }

            write_array(zg,"zgrid");
            write_reader(zg,"zgrid",false);*/
        }

        /* Initial time */
        double init_time() const { return initial_time; }
        int get_restart_sequence() const { return restart_sequence; }

        /* Modify the timestep if necessary in order to land evenly on a plot time */
        double check_timestep (double intime, double now) {
            if (intime > 0.5/sqrt(N2_max)) {
                intime = 0.5/sqrt(N2_max);
            }
            // Now, calculate how many timesteps remain until the next writeout
            double until_plot = last_plot + plot_interval - now;
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

        /* Initialize velocities at the start of the run.  For this simple
           case, initialize all velocities to 0 */
        void init_vels(DTArray & u, DTArray & v, DTArray & w) {
            if (restarting and (!restart_from_dump)) {
                init_vels_restart(u, v, w);
            }
            else if (restarting and restart_from_dump) {
                init_vels_dump(u, v, w);
            }
            else{
                u = 0; // Use the Blitz++ syntax for simple initialization
                v = 0; // of an entire (2D or 3D) array with a single line
                w = 0; // of code.
                /* Add random initial perturbation */
                if (perturb > 0) {
                    int myrank;
                    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
                    Normal<double> rnd(0,1);
                    for (int i = u.lbound(firstDim); i<= u.ubound(firstDim); i++) {
                        rnd.seed(i);
                        for (int j = u.lbound(secondDim); j<= u.ubound(secondDim); j++) {
                            for (int k = u.lbound(thirdDim); k<= u.ubound(thirdDim); k++) {
                                u(i,j,k) += perturb*rnd.random();
                                v(i,j,k) += perturb*rnd.random();
                                w(i,j,k) += perturb*rnd.random();
                            }
                        }
                    }
                }

                // Also, write out the initial velocities and proper M-file readers
                /*write_reader(u,"u",true);
                write_reader(w,"w",true);
                write_array(u,"u",plot_number);
                write_array(w,"w",plot_number);
                if (Ny > 1 || rot_f != 0) {
                    write_reader(v,"v",true);
                    write_array(v,"v",plot_number);
                }*/
                return;
            }
        }

        /* Initialze the tracers (density) */
        void init_tracers(vector<DTArray *> & tracers) {
            DTArray & rho = *tracers[RHO];

            assert (tracers.size() == NUM_TRACER);
            if (restarting and (!restart_from_dump)) {
                init_tracer_restart("rho",rho);
            }
            else if (restarting and restart_from_dump) {
                init_tracer_dump("rho",rho);
            }
            else {
                // background stratification
                rho =  -0.25*delta_rho*tanh(((*zgrid)(ii,jj,kk)-(MinZ/Lz+0.5+pyc_asym-0.5*pyc_sep_perc)*Lz)/h_halfwidth);
                rho += -0.25*delta_rho*tanh(((*zgrid)(ii,jj,kk)-(MinZ/Lz+0.5+pyc_asym+0.5*pyc_sep_perc)*Lz)/h_halfwidth);
                rho = rho*0.5*(1.0+tanh(((*xgrid)(ii,jj,kk)-Lmix)/delta_x));
                // mixed region
                rho = rho + 0.5*(1.0-tanh(((*xgrid)(ii,jj,kk)-Lmix)/delta_x))
                    *(-0.25*delta_rho)*(
                            1.0+tanh(((*zgrid)(ii,jj,kk)-(MinZ/Lz+0.5+pyc_asym)*Lz+0.5*Hmix)/h_mix_half)
                           -1.0+tanh(((*zgrid)(ii,jj,kk)-(MinZ/Lz+0.5+pyc_asym)*Lz-0.5*Hmix)/h_mix_half));

                // write out arrays
                /*write_array(rho,"rho",0); 
                write_reader(rho,"rho",true);*/
            }
        }

        /* Forcing in the momentum equations */
        void forcing(double t, const DTArray & u, DTArray & u_f,
                const DTArray & v, DTArray & v_f, const DTArray & w,
                DTArray & w_f, vector<DTArray *> & tracers,
                vector<DTArray *> & tracers_f) {
            DTArray & rho = *tracers[RHO];
            u_f = +rot_f*v;
            v_f = -rot_f*u;
            w_f = -g*rho/rho_0;
            *tracers_f[RHO] = 0;
        }

        /* Basic analysis, to write out the field periodically */
        void analysis(double time, DTArray & u, DTArray & v, DTArray & w,
                vector<DTArray *> & tracer, DTArray & pressure) {
            /* If it is very close to the plot time, write data fields to disk */
            if ((time - last_plot - plot_interval) > -1e-6) {
                plot_number++;
                t_step = MPI_Wtime(); // time just before write (for dump)
                /*write_array(u,"u",plot_number);
                write_array(w,"w",plot_number);
                if (Ny > 1 || rot_f != 0) {
                    write_array(v,"v",plot_number);
                }
                write_array(*tracer[RHO],"rho",plot_number);*/
                last_plot = last_plot + plot_interval;

                // Find average time to write (for dump)
                clock_time = MPI_Wtime(); // time just afer write
                avg_write_time = (avg_write_time*(plot_number-restart_sequence-1) + (clock_time - t_step))/
                    (plot_number-restart_sequence);
                if (master()){
                    fprintf(stdout,"Last write time: %.6g. Average write time: %.6g.\n", clock_time - t_step, avg_write_time);
                }
                if (master()) fprintf(stdout,"*");
            }
            // increase counter and update clocks
            itercount++;
            if (itercount == 1){
                step_start_time = MPI_Wtime();
            }
            if (master()) {
                clock_time = MPI_Wtime();
                t_step = clock_time - step_start_time;
            }

            // Also, calculate and write out useful information: maximum u, v, w...
            double max_u = psmax(max(abs(u)));
            double max_v = psmax(max(abs(v)));
            double max_w = psmax(max(abs(w)));
            double max_ke = psmax(max(0.5*rho_0*(u*u + v*v + w*w)*Lx/Nx*Ly/Ny*Lz/Nz)); //only true for uniform grid
            double ke_tot = pssum(sum(0.5*rho_0*(u*u + v*v + w*w)*Lx/Nx*Ly/Ny*Lz/Nz));  //only true for uniform grid
            double max_rho = psmax(max(abs(*tracer[RHO])));
            if (master() && itercount == 1){
                double t_startup;
                t_startup = clock_time - real_start_time;
                fprintf(stdout,"Start-up time: %g s.\n",t_startup);
            }
            if (master() && itercount == 1) fprintf(stderr,"[Iter], (Clock time), Sim time:, Max U, Max V, Max W, Max KE, Total KE, Max Density\n");
            if (master()) fprintf(stderr,"[%d] (%.12g) %.6f: %.6g %.6g %.6g %.6g %.6g %.6g\n",
                    itercount,t_step,time,max_u,max_v,max_w,max_ke,ke_tot,max_rho);

            // Write out the chain if savechain is true
            if (savechain){
                int Iout,Jout,myrank,pointflag;
                MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
                if ((time >= chain1_start_time)
                        && ((time - chain_last_plot - chain_plot_interval) > -1e-6)
                        && (time < chain1_end_time)) {
                    Iout = chain1_xindex; Jout = chain1_yindex; pointflag = 0;
                    if ( Iout >= u.lbound(firstDim) && Iout <= u.ubound(firstDim) ) pointflag=1;
                    if (pointflag==1) {
                        /*write_chain("u_chain.txt",u, Iout, Jout, time);
                        write_chain("v_chain.txt",v, Iout, Jout, time);
                        write_chain("w_chain.txt",w, Iout, Jout, time);
                        write_chain("rho_chain.txt",*tracer[RHO], Iout, Jout, time);*/
                    }
                    chain_last_plot = chain_last_plot + chain_plot_interval;
                }
            }

            // Determine last plot if restarting from the dump case
            if (restart_from_dump and (itercount == 1)){
                last_plot = restart_sequence*plot_interval;    
            }
            // see if close to end of compute time and dump
            check_and_dump(clock_time, real_start_time, compute_time, time, avg_write_time,
                    plot_number, u, v, w, tracer);
            // Change dump log file if successfully reached final time
            // the dump time will be twice final time so that a restart won't actually start
            successful_dump(plot_number, final_time, plot_interval);
        }

        // User specified variables to dump
        void write_variables(DTArray & u,DTArray & v, DTArray & w,
                vector<DTArray *> & tracer) {
            /*write_array(u,"u.dump",-1);
            write_array(v,"v.dump",-1);
            write_array(w,"w.dump",-1);
            write_array(*tracer[RHO],"rho.dump",-1);*/
        }

        // Constructor
        mapiw(): // Initialization list for xx, yy and zz 1d grids
            hill(split_range(Nx))
    {   // Initialize the local variables
        plot_number = restart_sequence;
        last_plot = restart_time;
        chain_last_plot = chain1_start_time - chain_plot_interval;
    }
};


/* The ``main'' routine */
int main(int argc, char ** argv) {
    /* Initialize MPI.  This is required even for single-processor runs,
       since the inner routines assume some degree of parallelization,
       even if it is trivial. */

    f_order = 4; f_cutoff = 0.8; f_strength = -0.33;
    MPI_Init(&argc, &argv);

    real_start_time = MPI_Wtime();     // for dump
    /* ------------------ Define parameters from spins.conf --------------------- */
    options_init(); // Initialize options

    option_category("Restart options");
    add_switch("restart",&restarting,"Restart from prior output time.  OVERRIDES many other values.");
    add_option("restart_time",&restart_time,0.0,"Time to restart from");
    add_option("restart_sequence",&restart_sequence,-1,"Sequence number to restart from");

    option_category("Grid Options");
    add_option("Lx",&Lx,"Length of tank");
    add_option("Ly",&Ly,"Width of tank");
    add_option("Lz",&Lz,"Height of tank");
    add_option("Nx",&Nx,"Number of points in X");
    add_option("Ny",&Ny,1,"Number of points in Y");
    add_option("Nz",&Nz,"Number of points in Z");
    add_option("min_x",&MinX,0.0,"Minimum X-value");
    add_option("min_y",&MinY,0.0,"Minimum Y-value");
    add_option("min_z",&MinZ,0.0,"Minimum Z-value");

    option_category("Physical constants");
    add_option("g",&g,9.81,"Gravitational acceleration");
    add_option("rot_f",&rot_f,0.0,"Coriolis frequency");

    option_category("Stratification parameters");
    add_option("rho_0",&rho_0,"Reference density");
    add_option("delta_rho",&delta_rho,"density difference");
    add_option("pyc_asym",&pyc_asym,"percentage of depth to shift pycnocline");
    add_option("pyc_sep_perc",&pyc_sep_perc,"total separation of double pycnocline (as perc. of depth)");
    add_option("h_perc",&h_perc,"Pycnocline half-width as perc. of depth");
    add_option("h_mix_perc",&h_mix_perc,"Pycnocline half-width as perc. of depth");
    add_option("delta_x",&delta_x,"Horizontal transition half-width");
    add_option("Lmix",&Lmix,"Width of mixed region");
    add_option("Hmix",&Hmix,"Height of mixed region");

    option_category("Hill parameters");
    add_option("hill_height",&hill_height,"Height of hill (as percentage of depth)");
    add_option("hill_centre",&hill_centre,"location of hill peak");
    add_option("hill_width",&hill_width,"Width of hill");

    add_option("visco",&VISCO,"Viscosity");
    add_option("diffu_rho",&DIFFU_rho,"Diffusivity of density");	

    add_option("plot_interval",&plot_interval,"Time between writes");
    add_option("final_time",&final_time,"Final time");
    add_option("savechain",&savechain,false,"Flag to have save vertical chains or not");
    add_option("chain1_start_time",&chain1_start_time,"Time to start writing chain");
    add_option("chain1_end_time",&chain1_end_time,"Time to stop writing chain");
    add_option("chain_plot_interval",&chain_plot_interval,"Time between writes in chain");

    add_option("perturb",&perturb,"Initial perturbation in velocity");

    option_category("Dumping options");
    add_option("compute_time",&compute_time,-1.0,"Time permitted for computation");
    add_option("restart_from_dump",&restart_from_dump,false,"If restart from dump");

    options_parse(argc,argv);

    // Now, make sense of the options received.  Many of these values
    // can be directly used, but the ones of string-type need further
    // procesing.

    // Read dump_time.txt and check if past final time
    if (restart_from_dump){
        restarting = true;
        string dump_str;
        ifstream dump_file;
        dump_file.open ("dump_time.txt");

        getline (dump_file,dump_str); // ingnore 1st line

        getline (dump_file,dump_str);
        restart_time = atof(dump_str.c_str());

        getline (dump_file,dump_str); // ingore 3rd line

        getline (dump_file,dump_str);
        restart_sequence = atoi(dump_str.c_str());

        if (restart_time > final_time){
            // Die, ungracefully
            if (master()){
                fprintf(stderr,"Restart dump time (%.4g) is past final time (%.4g). Quitting now.\n",restart_time,final_time);
            }
            MPI_Finalize(); exit(1);
        }
    }
    if (compute_time > 0){
        avg_write_time = max(100.0*Nx*Ny*Nz/pow(512.0,3), 20.0);
    }

    /* ------------------ Derived parameters --------------------- */

    // Vertical chain locations
    chain1_xindex = Nx/2;
    chain1_yindex = Ny/2;

    // Pycnocline half-width
    h_halfwidth = h_perc*Lz;
    h_mix_half  = h_mix_perc*Lz;

    // Diffusivity information
    DIFFU[RHO] = DIFFU_rho;

    // Mode-2 wave speed
    c0 = 0.5*sqrt(g*h_halfwidth*delta_rho/rho_0);

    // Maximum buoyancy frequency (squared) if the initial stratification was stable
    N2_max = g/rho_0*delta_rho/(2*h_halfwidth);

    // Reynolds number
    Re = c0*h_halfwidth/VISCO;

    /* ------------------ Set correct initial time, and sequence --------------------- */
    if (restarting) {
        if (restart_sequence <= 0) {
            restart_sequence = int(restart_time/plot_interval);
        }
        if (master()) {
            fprintf(stdout,"Restart flags detected\n");
            fprintf(stdout,"Restarting from time %g, at sequence number %d\n",
                    restart_time,restart_sequence);
        }
        initial_time = restart_time;
    } else {
        // Not restarting, so set the initial sequence number
        // to the initial time / plot_interval
        restart_sequence = int(initial_time/plot_interval);
        if (fmod(initial_time,plot_interval) != 0.0) {
            if (master()) {
                fprintf(stdout,"Warning: the initial time (%g) does not appear to be an even multiple of the plot interval (%g)\n",
                        initial_time,plot_interval);
            }
        }
    }

    /* ------------------ Print some parameters --------------------- */
    if (master()) {
        fprintf(stdout,"Dam break problem\n");
        fprintf(stdout,"Using a %f x %f x %f grid of %d x %d x %d points\n",Lx,Ly,Lz,Nx,Ny,Nz);
        fprintf(stdout,"g = %f, rot_f = %f, rho_0 = %f, delta_rho = %f\n",g,rot_f,rho_0,delta_rho);
        fprintf(stdout,"Pycnocline half-width as %% of depth: h_perc = %g\n",h_perc);
        fprintf(stdout,"Pycnocline half-width: h = %g\n",h_halfwidth);
        fprintf(stdout,"Pycnocline vertical shift %%: zeta = %g\n",pyc_asym);
        fprintf(stdout,"Pycnocline separation as %% of depth: zeta_p = %g\n",pyc_sep_perc);
        fprintf(stdout,"Width of mixed region: L_mix = %g\n",Lmix);
        fprintf(stdout,"Height of mixed region: H_mix = %g\n",Hmix);
        fprintf(stdout,"Time between plots: %g s\n",plot_interval);
        fprintf(stdout,"Chain 1 indices: x_i = %d, y_i = %d\n",chain1_xindex,chain1_yindex);
        fprintf(stdout,"Initial velocity perturbation: %g\n",perturb);

        fprintf(stdout,"Stably-stratified phase speed %g\n",c0);
        fprintf(stdout,"Buoyancy frequency squared %g\n",N2_max);
        fprintf(stdout,"Reynolds number %g\n",Re);
    }

    /* ------------------ Do stuff --------------------- */
    mapiw mycode; // Create an instantiated object of the above class
    /// Create a flow-evolver that takes its settings from the above class
    FluidEvolve<mapiw> do_mapiw(&mycode);

    // Initialize
    double start_time = MPI_Wtime();
    do_mapiw.initialize();

    // Run until the end of time
    do_mapiw.do_run(final_time);

   double now = MPI_Wtime();
   if (master()) fprintf(stderr,"Total runtime complete in %gs\n",now-start_time);
   {
      int numprocs, myrank;
      MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
      MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

      for (int i = 0; i < numprocs; i++) {
         if (myrank == i) {
            MPI_Barrier(MPI_COMM_WORLD);
            fprintf(stderr,"Processor %d timing report\n",myrank);
            timing_stack_report();
            MPI_Barrier(MPI_COMM_WORLD);
         } else {
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
         }
      }
   }
    //if (master()) timing_stack_report();
    MPI_Finalize(); // Cleanly exit MPI
    return 0; // End the program
}
