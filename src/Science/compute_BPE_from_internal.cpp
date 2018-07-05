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

// Internal energy converted to Background potential energy
// See Winters et al. 1995 (Available potential energy and mixing in density-stratified fluids)
// phi_i = - kappa * g * (integral of density at top boundary
//                        minus the integral of density at bottom boundary)
void compute_BPE_from_internal(double & phi_i, TArrayn::DTArray & rho,
        double kappa_rho, double rho_0, double g, int Nz, bool dimensional_rho) {
    // Tensor variables for indexing
    blitz::firstIndex ii;
    blitz::secondIndex jj;
    blitz::Range all = blitz::Range::all();

    if ( !dimensional_rho ) {
        phi_i = -kappa_rho * g * pssum(sum(
                    rho_0 * (rho(all,all,Nz-1) - rho(all,all,0))
                    * (*get_quad_x())(ii) * (*get_quad_y())(jj) ));
    } else {
        phi_i = -kappa_rho * g * pssum(sum(
                    (rho(all,all,Nz-1) - rho(all,all,0))
                    * (*get_quad_x())(ii) * (*get_quad_y())(jj) ));
    }
}
