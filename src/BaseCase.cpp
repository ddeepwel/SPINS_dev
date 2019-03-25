#include "BaseCase.hpp"
#include "Science.hpp"
#include "NSIntegrator.hpp"
#include "TArray.hpp"
#include <blitz/array.h>
#include <fstream>

//using namespace TArray;
using namespace NSIntegrator;
using blitz::Array;
using std::vector;

/* Call the source code writing function in the constructor */
extern "C" {
    void WriteCaseFileSource(void);
}
BaseCase::BaseCase(void)
{
    if (master()) WriteCaseFileSource();

    // Print version information
    if (master()) {
        fprintf(stdout,"-------------------------------------------\n");
        fprintf(stdout,"\nSPINS Version %d.%d.%d \n\n", MAJOR_VERSION, MINOR_VERSION, PATCH_VERSION);
        fprintf(stdout,"-------------------------------------------\n");
    }
}

/* Implementation of non-abstract methods in BaseCase */
int BaseCase::numActive() const { return 0; }
int BaseCase::numPassive() const { return 0; }
int BaseCase::numtracers() const { /* total number of tracers */
    return numActive() + numPassive();
}

int BaseCase::size_x() const {
    return size_cube();
}
int BaseCase::size_y() const {
    return size_cube();
}
int BaseCase::size_z() const {
    return size_cube();
}
double BaseCase::length_x() const {
    return length_cube();
}
double BaseCase::length_y() const {
    return length_cube();
}
double BaseCase::length_z() const {
    return length_cube();
}

DIMTYPE BaseCase::type_x() const {
    return type_default();
}
DIMTYPE BaseCase::type_y() const {
    return type_default();
}
DIMTYPE BaseCase::type_z() const {
    return type_default();
}
DIMTYPE BaseCase::type_default() const {
    return PERIODIC;
}

void BaseCase::tracer_bc_x(int t_num, double & dir, double & neu) const {
    if (!zero_tracer_boundary) {
        dir = 0; 
        neu = 1;
    }
    else {
        dir = 1;
        neu = 0;
    }
    return;
}
void BaseCase::tracer_bc_y(int t_num, double & dir, double & neu) const {
    if (!zero_tracer_boundary) {
        dir = 0; 
        neu = 1;
    }
    else {
        dir = 1;
        neu = 0;
    }
    return;
}
void BaseCase::tracer_bc_z(int t_num, double & dir, double & neu) const {
    if (!zero_tracer_boundary) {
        dir = 0; 
        neu = 1;
    }
    else {
        dir = 1;
        neu = 0;
    }
    return;
}
bool BaseCase::tracer_bc_forcing() const {
    return false;
}
bool BaseCase::is_mapped() const { // Whether this problem has mapped coordinates
    return false;
}
// Coordinate mapping proper, if is_mapped() returns true.  This features full,
// 3D arrays, but at least initially we're restricting ourselves to 2D (x,z)
// mappings
void BaseCase::do_mapping(DTArray & xgrid, DTArray & ygrid, DTArray & zgrid) {
    return;
}

/* Physical parameters */
double BaseCase::get_visco()            const { return 0; }
double BaseCase::get_diffusivity(int t) const { return 0; }
double BaseCase::get_rot_f()            const { return 0; }
int BaseCase::get_restart_sequence()    const { return 0; }
double BaseCase::get_next_plot()              { return 0; }

/* Initialization */
double BaseCase::init_time() const {
    return 0;
}
void BaseCase::init_tracers(vector<DTArray *> & tracers) {
    /* Initalize tracers one-by-one */
    if (numtracers() == 0) return; // No tracers, do nothing
    assert(numtracers() == int(tracers.size())); // Sanity check
    for (int i = 0; i < numtracers(); i++) {
        init_tracer(i, *(tracers[i]));
    }
}

double BaseCase::check_timestep(double t_step, double now) {
    // Check time step
    if (t_step < 1e-9) {
        // Timestep's too small, somehow stuff is blowing up
        if (master()) fprintf(stderr,"Tiny timestep (%e), aborting\n",t_step);
        return -1;
    } else if (t_step > get_dt_max()) {
        // Cap the maximum timestep size
        t_step = get_dt_max();
    }

    // Now, calculate how many timesteps remain until the next writeout
    double until_plot = get_next_plot() - now;
    int steps = ceil(until_plot / t_step);
    // Where will we be after (steps) timesteps of the current size?
    double real_until_plot = steps*t_step;
    // If that's close enough to the real writeout time, that's fine.
    if (fabs(until_plot - real_until_plot) < 1e-6) {
        return t_step;
    } else {
        // Otherwise, square up the timeteps.  This will always shrink the timestep.
    return (until_plot / steps);
    }
}

/* Forcing */
void BaseCase::forcing(double t,  DTArray & u, DTArray & u_f,
        DTArray & v, DTArray & v_f,  DTArray & w, DTArray & w_f,
        vector<DTArray *> & tracers, vector<DTArray *> & tracers_f) {
    /* First, if no active tracers then use the simpler form */
    if (numActive() == 0) {
        passive_forcing(t, u, u_f, v, v_f, w, w_f);
    } else {
        /* Look at split velocity-tracer forcing */
        vel_forcing(t, u_f, v_f, w_f, tracers);
        tracer_forcing(t, u, v, w, tracers_f);
    }
    /* And any/all passive tracers get 0 forcing */
    for (int i = 0; i < numPassive(); i++) {
        *(tracers_f[numActive() + i]) = 0;
    }
}
void BaseCase::passive_forcing(double t,  DTArray & u, DTArray & u_f,
        DTArray & v, DTArray & v_f,  DTArray &, DTArray & w_f) {
    /* Reduce to velocity-independent case */
    stationary_forcing(t, u_f, v_f, w_f);
}
void BaseCase::stationary_forcing(double t, DTArray & u_f, DTArray & v_f, 
        DTArray & w_f) {
    /* Default case, 0 forcing */
    u_f = 0;
    v_f = 0;
    w_f = 0;
}

/* Analysis */
void BaseCase::analysis(double t, DTArray & u, DTArray & v, DTArray & w,
        vector<DTArray *> tracer, DTArray & pres) {
    analysis(t,u,v,w,tracer);
}
void BaseCase::analysis(double t, DTArray & u, DTArray & v, DTArray & w,
        vector<DTArray *> tracer) {
    /* Do velocity and tracer analysis seperately */
    vel_analysis(t, u, v, w);
    for (int i = 0; i < numtracers(); i++) {
        tracer_analysis(t, i, *(tracer[i]));
    } 
}

void BaseCase::automatic_grid(double MinX, double MinY, double MinZ,
        Array<double,1> * xx, Array<double,1> * yy, Array<double,1> * zz){
    //Array<double,1> xx(split_range(size_x())), yy(size_y()), zz(size_z());
    bool xxa = false, yya = false, zza = false;
    if (!xx) {
        xxa = true; // Delete xx when returning
        xx = new Array<double,1>(split_range(size_x()));
    }
    if (!yy) {
        yya = true;
        yy = new Array<double,1>(size_y());
    }
    if (!zz) {
        zza = true;
        zz = new Array<double,1>(size_z());
    }
    Array<double,3> grid(alloc_lbound(size_x(),size_y(),size_z()),
            alloc_extent(size_x(),size_y(),size_z()),
            alloc_storage(size_x(),size_y(),size_z()));
    blitz::firstIndex ii;
    blitz::secondIndex jj;
    blitz::thirdIndex kk;

    // Generate 1D arrays
    if (type_x() == NO_SLIP) {
        *xx = MinX+length_x()*(0.5-0.5*cos(M_PI*ii/(size_x()-1)));
    } else {
        *xx = MinX + length_x()*(ii+0.5)/size_x();
    }
    *yy = MinY + length_y()*(ii+0.5)/size_y();
    if (type_z() == NO_SLIP) {
        *zz = MinZ+length_z()*(0.5-0.5*cos(M_PI*ii/(size_z()-1)));
    } else {
        *zz = MinZ + length_z()*(0.5+ii)/size_z();
    }

    // Write grid/reader
    grid = (*xx)(ii) + 0*jj + 0*kk;
    write_array(grid,"xgrid");
    write_reader(grid,"xgrid",false);

    if (size_y() > 1) {
        grid = 0*ii + (*yy)(jj) + 0*kk;
        write_array(grid,"ygrid");
        write_reader(grid,"ygrid",false);
    }

    grid = 0*ii + 0*jj + (*zz)(kk);
    write_array(grid,"zgrid");
    write_reader(grid,"zgrid",false);

    // Clean up
    if (xxa) delete xx;
    if (yya) delete yy;
    if (zza) delete zz;
}

/* Read velocities from regular output */
void BaseCase::init_vels_restart(DTArray & u, DTArray & v, DTArray & w){
    /* Restarting, so build the proper filenames and load the data into u, v, w */
    char filename[100];

    /* u */
    snprintf(filename,100,"u.%d",get_restart_sequence());
    if (master()) fprintf(stdout,"Reading u from %s\n",filename);
    read_array(u,filename,size_x(),size_y(),size_z());

    /* v, only necessary if this is an actual 3D run or if
       rotation is noNzero */
    if (size_y() > 1 || get_rot_f() != 0) {
        snprintf(filename,100,"v.%d",get_restart_sequence());
        if (master()) fprintf(stdout,"Reading v from %s\n",filename);
        read_array(v,filename,size_x(),size_y(),size_z());
    }

    /* w */
    snprintf(filename,100,"w.%d",get_restart_sequence());
    if (master()) fprintf(stdout,"Reading w from %s\n",filename);
    read_array(w,filename,size_x(),size_y(),size_z());
    return;
}

/* Read velocities from dump output */
void BaseCase::init_vels_dump(DTArray & u, DTArray & v, DTArray & w){
    /* Restarting, so build the proper filenames and load the data into u, v, w */

    /* u */
    if (master()) fprintf(stdout,"Reading u from u.dump\n");
    read_array(u,"u.dump",size_x(),size_y(),size_z());

    /* v, only necessary if this is an actual 3D run or if
       rotation is noNzero */
    if (size_y() > 1 || get_rot_f() != 0) {
        if (master()) fprintf(stdout,"Reading v from v.dump\n");
        read_array(v,"v.dump",size_x(),size_y(),size_z());
    }

    /* w */
    if (master()) fprintf(stdout,"Reading w from w.dump\n");
    read_array(w,"w.dump",size_x(),size_y(),size_z());
    return;
}


/* Read field from regular output */
void BaseCase::init_tracer_restart(const std::string & field, DTArray & the_tracer){
    /* Restarting, so build the proper filenames and load the data */
    char filename[100];

    snprintf(filename,100,"%s.%d",field.c_str(),get_restart_sequence());
    if (master()) fprintf(stdout,"Reading %s from %s\n",field.c_str(),filename);
    read_array(the_tracer,filename,size_x(),size_y(),size_z());
    return;
}

/* Read field from dump output */
void BaseCase::init_tracer_dump(const std::string & field, DTArray & the_tracer){
    /* Restarting, so build the proper filenames and load the data */
    char filename[100];

    snprintf(filename,100,"%s.dump",field.c_str());
    if (master()) fprintf(stdout,"Reading %s from %s\n",field.c_str(),filename);
    read_array(the_tracer,filename,size_x(),size_y(),size_z());
    return;
}

/* Read field from input data type */
void BaseCase::init_field(const std::string & field,
        const std::string & filename, DTArray & the_field, input_types input_data_type) {
    if (input_data_type == MATLAB) {
        // from matlab
        if (master()) fprintf(stdout,"Reading MATLAB-format %s (%d x %d) from %s\n",
                field.c_str(),size_x(),size_z(),filename.c_str());
        read_2d_slice(the_field,filename.c_str(),size_x(),size_z());
    } else if (input_data_type == CTYPE) {
        // from 2D spins output (generally, to go to 3D)
        if (master()) fprintf(stdout,"Reading CTYPE-format %s (%d x %d) from %s\n",
                field.c_str(),size_x(),size_z(),filename.c_str());
        read_2d_restart(the_field,filename.c_str(),size_x(),size_z());
    } else {
        // from full (2D or 3D) spins output
        // only suitable for initializing the grid
        if (master()) fprintf(stdout,"Reading %s from %s\n",field.c_str(),filename.c_str());
        read_array(the_field,filename.c_str(),size_x(),size_y(),size_z());
    }
}

/* write out vertical chain of data */
void BaseCase::write_chain(const char *filename, DTArray & val, int Iout, int Jout, double time) {
    FILE *fid=fopen(filename,"a");
    if (fid == NULL) {
        fprintf(stderr,"Unable to open %s for writing\n",filename);
        MPI_Finalize(); exit(1);
    }
    fprintf(fid,"%g",time);
    for (int ki=0; ki<size_z(); ki++) fprintf(fid," %g",val(Iout,Jout,ki));
    fprintf(fid,"\n");
    fclose(fid);
}

/* Read grid from regular output */
void BaseCase::init_grid_restart(const std::string & component,
        const std::string & filename, DTArray & grid) {
    if (master()) fprintf(stdout,"Reading %s from %s\n",component.c_str(),filename.c_str());
    read_array(grid,filename.c_str(),size_x(),size_y(),size_z());
    return;
}


/* Check and dump */
void BaseCase::check_and_dump(double clock_time, double real_start_time,
        double compute_time, double sim_time, double avg_write_time, int plot_number, int iter,
        DTArray & u, DTArray & v, DTArray & w, vector<DTArray *> & tracer){
    if (compute_time > 0) {
        // default is to not dump variables
        int do_dump = 0;

        // check if close to end of compute time
        if (master()) {
            double total_run_time = clock_time - real_start_time;
            double avg_clk_step = total_run_time/(iter+1); // +1 to accound for initial iter=0
            if (compute_time - total_run_time < 5*(avg_write_time + avg_clk_step)) {
                do_dump = 1; // true
            }
        }

        // Broadcast to all processors whether to dump or not
        MPI_Bcast(&do_dump, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // Dump variables if close to end of compute time
        if (do_dump == 1) {
            if (master()){
                fprintf(stdout,"Too close to final time, dumping!\n");
            }
            write_variables(u, v, w, tracer);

            // Write the dump time to a text file for restart purposes
            if (master()){
                FILE * dump_file;
                dump_file = fopen("dump_time.txt","w");
                assert(dump_file);
                fprintf(dump_file,"The dump time was:\n%.17g\n", sim_time);
                fprintf(dump_file,"The dump index was:\n%d\n", plot_number);
                fclose(dump_file);
            }

            // Die gracefully
            MPI_Finalize(); exit(0);
        }
    }
}

/* Change dump log file for successful completion*/
void BaseCase::successful_dump(int plot_number, double final_time, double plot_interval){
    if (master() && (plot_number == final_time/plot_interval)){
        // Write the dump time to a text file for restart purposes
        FILE * dump_file; 
        dump_file = fopen("dump_time.txt","w");
        assert(dump_file);
        fprintf(dump_file,"The dump 'time' was:\n%.12g\n", 2*final_time);
        fprintf(dump_file,"The dump index was:\n%d\n", plot_number);
        fclose(dump_file);
    }
}

// Explicitly instantiate common add_diagnostic templates, for integer and double
// types
template void BaseCase::add_diagnostic<int>(const string str, const int val,
        string & header, string & line);
template void BaseCase::add_diagnostic<double>(const string str, const double val,
        string & header, string & line);

// write the diagnostic file
void BaseCase::write_diagnostics(string header, string line,
        int iter, bool restarting) {
    // remove last two elements (comma and space)
    string clean_header = header.substr(0, header.size()-2);
    string clean_line   =   line.substr(0,   line.size()-2);
    // open file
    FILE * diagnos_file = fopen("diagnostics.txt","a");
    assert(diagnos_file);
    // print header
    if (iter == 0 and !restarting) {
        fprintf(diagnos_file,"%s\n",clean_header.c_str());
    }
    // print the line of values
    fprintf(diagnos_file, "%s\n", clean_line.c_str());
    // Close file
    fclose(diagnos_file);
}

// Write plot time information
void BaseCase::write_plot_times(double time, double clock_time, double comp_duration,
        double avg_write_time, int plot_number, bool restarting) {
    if (master()) {
        // in log file
        fprintf(stdout,"*Write time: %.6g. Average write time: %.6g.\n",
                clock_time - comp_duration, avg_write_time);
        // track in a file
        FILE * plottimes_file = fopen("plot_times.txt","a");
        assert(plottimes_file);
        if ( plot_number==get_restart_sequence()+1 and !restarting )
            fprintf(plottimes_file,"Output number, Simulation time (s), "
                    "Write time (s), Average write time (s)\n");
        fprintf(plottimes_file,"%d, %.17f, %.12g, %.12g\n",
                plot_number, time, clock_time - comp_duration, avg_write_time);
        fclose(plottimes_file);
    }
}

// parse file types
void parse_datatype(const string datatype, input_types & input_data_type) {
    if (datatype=="MATLAB") { input_data_type = MATLAB; }
    else if (datatype == "CTYPE") { input_data_type = CTYPE; }
    else if (datatype == "FULL") { input_data_type = FULL; }
    else {
        if (master())
            fprintf(stderr,"Invalid option %s received for file_type\n",datatype.c_str());
        MPI_Finalize(); exit(1);
    }
}

// parse expansion types 
void parse_boundary_conditions(const string xgrid_type, const string ygrid_type,
        const string zgrid_type, DIMTYPE & intype_x, DIMTYPE & intype_y, DIMTYPE & intype_z) {
    // x
    if      (xgrid_type == "FOURIER")   { intype_x = PERIODIC; }
    else if (xgrid_type == "FREE_SLIP") { intype_x = FREE_SLIP; }
    else if (xgrid_type == "NO_SLIP")   { intype_x = NO_SLIP; }
    else {
        if (master())
            fprintf(stderr,"Invalid option %s received for type_x\n",xgrid_type.c_str());
        MPI_Finalize(); exit(1);
    }
    // y
    if      (ygrid_type == "FOURIER")   { intype_y = PERIODIC; }
    else if (ygrid_type == "FREE_SLIP") { intype_y = FREE_SLIP; }
    else {
        if (master())
            fprintf(stderr,"Invalid option %s received for type_y\n",ygrid_type.c_str());
        MPI_Finalize(); exit(1);
    }
    // z
    if      (zgrid_type == "FOURIER")   { intype_z = PERIODIC; }
    else if (zgrid_type == "FREE_SLIP") { intype_z = FREE_SLIP; }
    else if (zgrid_type == "NO_SLIP")   { intype_z = NO_SLIP; }
    else {
        if (master())
            fprintf(stderr,"Invalid option %s received for type_z\n",zgrid_type.c_str());
        MPI_Finalize(); exit(1);
    }
}

// Top stresses
void BaseCase::stresses_top(TArrayn::DTArray & u, TArrayn::DTArray & v, TArrayn::DTArray & w,
        TArrayn::DTArray & Hprime, TArrayn::DTArray & temp, TArrayn::Grad * gradient_op,
        const string * grid_type, const double mu, double time, int itercount, bool restarting) {
    // set-up
    static DTArray *tx = alloc_array(size_x(),size_y(),1);
    static DTArray *ty = alloc_array(size_x(),size_y(),1);
    blitz::firstIndex ii;
    blitz::secondIndex jj;

    // top stress ( along channel - x )
    top_stress_x(*tx, u, temp, gradient_op, grid_type, size_z(), mu);
    // integral values (the first is the surface force)
    double fx_tot = pssum(sum(
                (*get_quad_x())(ii)*
                (*get_quad_y())(jj)*(*tx)));
    double fx_abs = pssum(sum(
                (*get_quad_x())(ii)*
                (*get_quad_y())(jj)*abs(*tx)));
    // extrema values
    double tx_max = psmax(max(*tx));
    double tx_min = psmin(min(*tx));

    // top stress ( across channel - y )
    top_stress_y(*ty, v, temp, gradient_op, grid_type, size_z(), mu);
    // integral values (the first is the surface force)
    double fy_tot = pssum(sum(
                (*get_quad_x())(ii)*
                (*get_quad_y())(jj)*(*ty)));
    double fy_abs = pssum(sum(
                (*get_quad_x())(ii)*
                (*get_quad_y())(jj)*abs(*ty)));
    // extrema values
    double ty_max = psmax(max(*ty));
    double ty_min = psmin(min(*ty));
    double ts_max = psmax(max(pow(pow(*tx,2)+pow(*ty,2),0.5)));

    // write to a stress diagnostic file
    if (master()) {
        FILE * stresses_file = fopen("stresses_top.txt","a");
        assert(stresses_file);
        if ( itercount==0 and !restarting )
            fprintf(stresses_file,"Time, "
                    "tx_max, tx_min, ty_max, ty_min, ts_max, "
                    "fx_tot, fx_abs, fy_tot, fy_abs\n");
        fprintf(stresses_file,"%.17f, "
                "%.17g, %.17g, %.17g, %.17g, %.17g, "
                "%.17g, %.17g, %.17g, %.17g\n",
                time,
                tx_max, tx_min, ty_max, ty_min, ts_max,
                fx_tot, fx_abs, fy_tot, fy_abs);
        fclose(stresses_file);
    }
}

// Bottom stresses
void BaseCase::stresses_bottom(TArrayn::DTArray & u, TArrayn::DTArray & v, TArrayn::DTArray & w,
        TArrayn::DTArray & Hprime, TArrayn::DTArray & temp, TArrayn::Grad * gradient_op,
        const string * grid_type, const double mu, double time, int itercount, bool restarting) {
    // set-up
    static DTArray *tx = alloc_array(size_x(),size_y(),1);
    static DTArray *ty = alloc_array(size_x(),size_y(),1);
    blitz::firstIndex ii;
    blitz::secondIndex jj;

    // bottom stress ( along channel - x )
    bottom_stress_x(*tx, Hprime, u, w, temp, gradient_op, grid_type, is_mapped(), mu);
    // integral values (the first is the surface force)
    double fx_tot = pssum(sum(
                (*get_quad_x())(ii)*pow(1+pow(Hprime,2),0.5)*
                (*get_quad_y())(jj)*(*tx)));
    double fx_abs = pssum(sum(
                (*get_quad_x())(ii)*pow(1+pow(Hprime,2),0.5)*
                (*get_quad_y())(jj)*abs(*tx)));
    // extrema values
    double tx_max = psmax(max(*tx));
    double tx_min = psmin(min(*tx));

    // bottom stress ( across channel - y )
    bottom_stress_y(*ty, Hprime, v, temp, gradient_op, grid_type, is_mapped(), mu);
    // integral values (the first is the surface force)
    double fy_tot = pssum(sum(
                (*get_quad_x())(ii)*pow(1+pow(Hprime,2),0.5)*
                (*get_quad_y())(jj)*(*ty)));
    double fy_abs = pssum(sum(
                (*get_quad_x())(ii)*pow(1+pow(Hprime,2),0.5)*
                (*get_quad_y())(jj)*abs(*ty)));
    // extrema values
    double ty_max = psmax(max(*ty));
    double ty_min = psmin(min(*ty));
    double ts_max = psmax(max(pow(pow(*tx,2)+pow(*ty,2),0.5)));

    // write to a stress diagnostic file
    if (master()) {
        FILE * stresses_file = fopen("stresses_bottom.txt","a");
        assert(stresses_file);
        if ( itercount==0 and !restarting )
            fprintf(stresses_file,"Time, "
                    "tx_max, tx_min, ty_max, ty_min, ts_max, "
                    "fx_tot, fx_abs, fy_tot, fy_abs\n");
        fprintf(stresses_file,"%.17f, "
                "%.17g, %.17g, %.17g, %.17g, %.17g, "
                "%.17g, %.17g, %.17g, %.17g\n",
                time,
                tx_max, tx_min, ty_max, ty_min, ts_max,
                fx_tot, fx_abs, fy_tot, fy_abs);
        fclose(stresses_file);
    }
}
