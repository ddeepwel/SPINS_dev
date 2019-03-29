#include "../BaseCase.hpp"

void BaseCase::automatic_grid(double MinX, double MinY, double MinZ,
        Array<double,1> * xx, Array<double,1> * yy, Array<double,1> * zz) {
    // create arrays if needed
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
