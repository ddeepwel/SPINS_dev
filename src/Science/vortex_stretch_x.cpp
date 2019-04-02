#include "../Science.hpp"
using namespace Transformer;

// x component of vortex stretching (omega dot grad) u
void vortex_stretch_x(TArrayn::DTArray & vort_stretch, TArrayn::DTArray & u,
        TArrayn::DTArray & v, TArrayn::DTArray & w, TArrayn::DTArray & temp1,
        TArrayn::DTArray & temp2, TArrayn::Grad * gradient_op, const string * grid_type) {
    // Set-up
    S_EXP expan[3];
    assert(gradient_op);

    // x-vorticity
    compute_vort_x(temp1, v, w, gradient_op, grid_type);
    // get du/dx
    find_expansion(grid_type, expan, "u");
    gradient_op->setup_array(&u,expan[0],expan[1],expan[2]);
    gradient_op->get_dx(&temp2,false);
    // omega_x * du/dx
    vort_stretch = temp1 * temp2;

    // y-vorticity
    compute_vort_y(temp1, u, w, gradient_op, grid_type);
    // du/dy
    find_expansion(grid_type, expan, "u");
    gradient_op->setup_array(&u,expan[0],expan[1],expan[2]);
    gradient_op->get_dy(&temp2,false);
    // omega_y * du/dy
    vort_stretch += temp1 * temp2;

    // z-vorticity
    compute_vort_z(temp1, u, v, gradient_op, grid_type);
    // du/dz
    find_expansion(grid_type, expan, "u");
    gradient_op->setup_array(&u,expan[0],expan[1],expan[2]);
    gradient_op->get_dz(&temp2,false);
    // omega_z * du/dz
    vort_stretch += temp1 * temp2;
}
