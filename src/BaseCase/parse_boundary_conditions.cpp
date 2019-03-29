#include "../BaseCase.hpp"

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
