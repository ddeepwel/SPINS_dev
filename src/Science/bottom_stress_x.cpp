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

// Bottom Stress (along channel - x)
void bottom_stress_x(TArrayn::DTArray & stress_x, TArrayn::DTArray & Hprime,
        TArrayn::DTArray & u, TArrayn::DTArray & w, TArrayn::DTArray & temp,
        TArrayn::Grad * gradient_op, const string * grid_type, const bool mapped,
        const double visco) {
    // Set-up
    blitz::Range all = blitz::Range::all();
    S_EXP expan[3];
    assert(gradient_op);
    assert((grid_type[2] == "NO_SLIP") && "Surface stress requires NO_SLIP vertical (z) boundary condition");
    // Along channel bottom stress can also be slightly inaccurate when x boundary condition is free slip
    // since u or w do not necessarily satisfy Neumann BCs and the derivative has Gibbs at domain edges
    // This is only true if there is velocity at left or right (in x) end wall corner, so is likely always small

    if (mapped) {
        // -du/dx
        find_expansion(grid_type, expan, "u");
        gradient_op->setup_array(&u,expan[0],expan[1],expan[2]);
        gradient_op->get_dx(&temp,false);
        temp = (-1)*temp;
        // dw/dz
        find_expansion(grid_type, expan, "w");
        gradient_op->setup_array(&w,expan[0],expan[1],expan[2]);
        gradient_op->get_dz(&temp,true);
        // 2H'*(w_z-u_x)
        stress_x(all,all,0) = 2*Hprime(all,all,0)*temp(all,all,0);

        // dw/dx
        find_expansion(grid_type, expan, "w");
        gradient_op->setup_array(&w,expan[0],expan[1],expan[2]);
        gradient_op->get_dx(&temp,false);
        // du/dz
        find_expansion(grid_type, expan, "u");
        gradient_op->setup_array(&u,expan[0],expan[1],expan[2]);
        gradient_op->get_dz(&temp,true);
        // (1-(H')^2)*(u_z+w_x)
        stress_x(all,all,0) += (1-pow(Hprime(all,all,0),2))*temp(all,all,0);
        // multiply by mu/(1+(H')^2)
        stress_x = visco/(1+pow(Hprime,2))*stress_x;
    } else {
        // du/dz
        find_expansion(grid_type, expan, "u");
        gradient_op->setup_array(&u,expan[0],expan[1],expan[2]);
        gradient_op->get_dz(&temp,false);
        // bottom stress
        stress_x(all,all,0) = visco*temp(all,all,0);
    }
}
