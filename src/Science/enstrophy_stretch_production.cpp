#include "../Science.hpp"
/*#include "math.h"
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
*/
using namespace Transformer;

// Enstrophy production from vortex stetching/tilting
// Defined as   omega_i * omega_j * S_{ij}
//      or as   omega_i * omega_j * du_i/dx_j
// where index notation has been used (See Davidson)
// This is the vorticity dotted with the vorticity stretching term in the vorticity equation
void enstrophy_stretch_production(TArrayn::DTArray & enst_prod, TArrayn::DTArray & u,
        TArrayn::DTArray & v, TArrayn::DTArray & w, TArrayn::DTArray & temp1,
        TArrayn::DTArray & temp2, TArrayn::DTArray & temp3, TArrayn::Grad * gradient_op,
        const string * grid_type) {
    // Set-up
    assert(gradient_op);

    // x-dimension contribution
    vortex_stretch_x(temp3, u, v, w, temp1, temp2, gradient_op, grid_type);
    compute_vort_x(temp1, v, w, gradient_op, grid_type);
    // omega_x * vortex_stretch_x
    enst_prod = temp1 * temp3;

    // y-dimension contribution
    vortex_stretch_y(temp3, u, v, w, temp1, temp2, gradient_op, grid_type);
    compute_vort_y(temp1, u, w, gradient_op, grid_type);
    // omega_y * vortex_stretch_y
    enst_prod += temp1 * temp3;

    // z-dimension contribution
    vortex_stretch_z(temp3, u, v, w, temp1, temp2, gradient_op, grid_type);
    compute_vort_z(temp1, u, v, gradient_op, grid_type);
    // omega_z * vortex_stretch_z
    enst_prod += temp1 * temp3;
}
