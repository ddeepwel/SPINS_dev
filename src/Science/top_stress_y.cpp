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


// Top Stress (across channel - y)
void top_stress_y(TArrayn::DTArray & stress_y, TArrayn::DTArray & v,
        TArrayn::DTArray & temp, TArrayn::Grad * gradient_op,
        const string * grid_type, const int Nz, const double visco) {
    // Set-up
    blitz::Range all = blitz::Range::all();
    S_EXP expan[3];
    assert(gradient_op);
    assert((grid_type[2] == "NO_SLIP") && "Surface stress requires NO_SLIP vertical (z) boundary condition");

    // dv/dz
    find_expansion(grid_type, expan, "v");
    gradient_op->setup_array(&v,expan[0],expan[1],expan[2]);
    gradient_op->get_dz(&temp,false);
    // top stress
    stress_y(all,all,0) = -visco*temp(all,all,Nz-1);
}
