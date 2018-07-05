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


// Z-component of vorticity
void compute_vort_z(TArrayn::DTArray & vortz, TArrayn::DTArray & u, TArrayn::DTArray & v,
       TArrayn::Grad * gradient_op, const string * grid_type) {
    // Set-up
    S_EXP expan[3];
    assert(gradient_op);

    // Setup for du/dy
    find_expansion(grid_type, expan, "u");
    gradient_op->setup_array(&u,expan[0],expan[1],expan[2]);
    // get du/dy
    gradient_op->get_dy(&vortz,false);
    // Invert to get the negative
    vortz = (-1)*vortz;

    // Setup for dv/dx
    find_expansion(grid_type, expan, "v");
    gradient_op->setup_array(&v,expan[0],expan[1],expan[2]);
    // get dv/dx, and add to vortz
    gradient_op->get_dx(&vortz,true);
}
