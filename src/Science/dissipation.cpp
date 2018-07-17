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


// Viscous dissipation: 2*mu*e_ij*e_ij
void dissipation(TArrayn::DTArray & diss, TArrayn::DTArray & u, TArrayn::DTArray & v,
        TArrayn::DTArray & w, TArrayn::Grad * gradient_op, const string * grid_type,
        const int Nx, const int Ny, const int Nz, const double visco) {
    // Set-up
    static DTArray *temp = alloc_array(Nx,Ny,Nz);
    S_EXP expan[3];
    assert(gradient_op);

    // 1st term: e_11^2 = (du/dx)^2
    find_expansion(grid_type, expan, "u");
    gradient_op->setup_array(&u,expan[0],expan[1],expan[2]);
    gradient_op->get_dx(temp,false);
    diss = pow(*temp,2);
    // 2nd term: e_22^2 = (dv/dy)^2
    find_expansion(grid_type, expan, "v");
    gradient_op->setup_array(&v,expan[0],expan[1],expan[2]);
    gradient_op->get_dy(temp,false);
    diss += pow(*temp,2);
    // 3rd term: e_33^2 = (dw/dz)^2
    find_expansion(grid_type, expan, "w");
    gradient_op->setup_array(&w,expan[0],expan[1],expan[2]);
    gradient_op->get_dz(temp,false);
    diss += pow(*temp,2);
    // 4th term: 2e_12^2 = 2*(1/2*(u_y + v_x))^2
    // u_y
    find_expansion(grid_type, expan, "u");
    gradient_op->setup_array(&u,expan[0],expan[1],expan[2]);
    gradient_op->get_dy(temp,false);
    // v_x
    find_expansion(grid_type, expan, "v");
    gradient_op->setup_array(&v,expan[0],expan[1],expan[2]);
    gradient_op->get_dx(temp,true);
    diss += 2.0*pow(0.5*(*temp),2);
    // 5th term: 2e_13^2 = 2*(1/2*(u_z + w_x))^2
    // u_z
    find_expansion(grid_type, expan, "u");
    gradient_op->setup_array(&u,expan[0],expan[1],expan[2]);
    gradient_op->get_dz(temp,false);
    // w_x
    find_expansion(grid_type, expan, "w");
    gradient_op->setup_array(&w,expan[0],expan[1],expan[2]);
    gradient_op->get_dx(temp,true);
    diss += 2.0*pow(0.5*(*temp),2);
    // 6th term: 2e_23^2 = 2*(1/2*(v_z + w_y))^2
    // v_z
    find_expansion(grid_type, expan, "v");
    gradient_op->setup_array(&v,expan[0],expan[1],expan[2]);
    gradient_op->get_dz(temp,false);
    // w_y
    find_expansion(grid_type, expan, "w");
    gradient_op->setup_array(&w,expan[0],expan[1],expan[2]);
    gradient_op->get_dy(temp,true);
    diss += 2.0*pow(0.5*(*temp),2);
    // multiply by 2*mu
    diss *= 2.0*visco;
}
