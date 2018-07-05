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


// Enstrophy Density: 1/2*(vort_x^2 + vort_y^2 + vort_z^2)
void enstrophy_density(TArrayn::DTArray & enst, TArrayn::DTArray & u, TArrayn::DTArray & v,
        TArrayn::DTArray & w, TArrayn::Grad * gradient_op, const string * grid_type,
        const int Nx, const int Ny, const int Nz) {
    // initalize temporary array
    static DTArray *temp = alloc_array(Nx,Ny,Nz);

    // square vorticity components
    compute_vort_x(*temp, v, w, gradient_op, grid_type);
    enst = pow(*temp,2);
    compute_vort_y(*temp, u, w, gradient_op, grid_type);
    enst += pow(*temp,2);
    compute_vort_z(*temp, u, v, gradient_op, grid_type);
    enst += pow(*temp,2);
    enst = 0.5*enst;
}
