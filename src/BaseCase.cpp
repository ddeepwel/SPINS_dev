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

/* Read grid from regular output */
void BaseCase::init_grid_restart(const std::string & component,
        const std::string & filename, DTArray & grid) {
    if (master()) fprintf(stdout,"Reading %s from %s\n",component.c_str(),filename.c_str());
    read_array(grid,filename.c_str(),size_x(),size_y(),size_z());
    return;
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
