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

// Bottom Stress (across channel - y)
void bottom_stress_y(TArrayn::DTArray & stress_y, TArrayn::DTArray & Hprime,
        TArrayn::DTArray & v, TArrayn::DTArray & temp,
        TArrayn::Grad * gradient_op, const string * grid_type, const bool mapped,
        const double visco) {
    // Set-up
    blitz::Range all = blitz::Range::all();
    S_EXP expan[3];
    assert(gradient_op);
    assert((grid_type[2] == "NO_SLIP") && "Surface stress requires NO_SLIP vertical (z) boundary condition");
    // Across channel bottom stress can also be slightly inaccurate when x boundary condition is free slip
    // since v does not necessarily satisfy Neumann BCs and the derivative has Gibbs at domain edges
    // This is only true if there is velocity at left or right (in x) end wall corner, so is likely always small

    if (mapped) {
        // dv/dx
        find_expansion(grid_type, expan, "v");
        gradient_op->setup_array(&v,expan[0],expan[1],expan[2]);
        gradient_op->get_dx(&temp,false);
        // -v_x*H'
        stress_y(all,all,0) = -temp(all,all,0)*Hprime(all,all,0);
        // dv/dz
        gradient_op->setup_array(&v,expan[0],expan[1],expan[2]);
        gradient_op->get_dz(&temp,false);
        // add to -v_x*H'
        stress_y(all,all,0) = temp(all,all,0) + stress_y(all,all,0);
        // multiply by mu/sqrt(1+(H')^2)
        stress_y = visco/pow(1+pow(Hprime,2),0.5)*stress_y;
    } else {
        // dv/dz
        find_expansion(grid_type, expan, "v");
        gradient_op->setup_array(&v,expan[0],expan[1],expan[2]);
        gradient_op->get_dz(&temp,false);
        // bottom stress
        stress_y(all,all,0) = visco*temp(all,all,0);
    }
}
