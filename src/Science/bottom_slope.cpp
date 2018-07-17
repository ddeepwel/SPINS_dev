#include "../Science.hpp"
#include "math.h"
#include "../Par_util.hpp"
#include "stdio.h"
#include "../Split_reader.hpp"
#include "../T_util.hpp"
#include "../Parformer.hpp"
#include "../Sorter.hpp"
#include <numeric>

using blitz::Array;
using blitz::cos;
using namespace TArrayn;
using namespace NSIntegrator;
using namespace Transformer;


// Bottom slope
void bottom_slope(TArrayn::DTArray & Hprime, TArrayn::DTArray & zgrid,
        TArrayn::DTArray & temp, TArrayn::Grad * gradient_op,
        const string * grid_type, const int Nx, const int Ny, const int Nz) {
    // Set-up
    DTArray & z_x = *alloc_array(Nx,Ny,Nz);
    blitz::Range all = blitz::Range::all();
    blitz::firstIndex ii;
    blitz::secondIndex jj;
    blitz::thirdIndex kk;
    S_EXP expan[3];
    assert(gradient_op);

    // get bottom topography
    Array<double,1> xx(split_range(Nx));
    xx = zgrid(all,0,0);
    // put into temp array, and take derivative
    temp = xx(ii) + 0*jj + 0*kk;
    find_expansion(grid_type, expan, "zgrid");
    gradient_op->setup_array(&temp,expan[0],expan[1],expan[2]);
    gradient_op->get_dx(&z_x);
    // flatten to get 2D array
    Hprime(all,all,0) = z_x(all,all,0);
    delete &z_x, xx;
}
