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


void compute_vorticity(TArrayn::DTArray & vortx, TArrayn::DTArray & vorty, TArrayn::DTArray & vortz,
        TArrayn::DTArray & u, TArrayn::DTArray & v, TArrayn::DTArray & w,
        TArrayn::Grad * gradient_op, const string * grid_type) {
    // compute each component
    compute_vort_x(vortx, v, w, gradient_op, grid_type);
    compute_vort_y(vorty, u, w, gradient_op, grid_type);
    compute_vort_z(vortz, u, v, gradient_op, grid_type);
}
