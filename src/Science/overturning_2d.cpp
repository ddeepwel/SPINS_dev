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


/* Marek's overturning diagnotic.  For a given density, compute the total sum
   of z-levels (by vertical column) for which the density is statically 
   unstable.  That is, a completely stable stratification will return an Array
   full of 0s, and a completely unstable stratification will return an Array
   full of (zmax - zmin). */
Array<double,3> overturning_2d(
        Array<double,3> const & rho, // Density
        Array<double,1> const & zgrid, // Z-levels
        Dimension reduce // Dimension to consider the vertical
        ) {
    using namespace TArrayn;
    blitz::RectDomain<3> dims = rho.domain();
    // Now, for general behaviour over reduced dimensions, figure out min/max
    // for the output Array (and for the iteration inside)
    int szxmin, szymin, szzmin, szxmax, szymax, szzmax;
    szxmin = dims.lbound(firstDim);
    szymin = dims.lbound(secondDim);
    szzmin = dims.lbound(thirdDim);
    szxmax = dims.ubound(firstDim);
    szymax = dims.ubound(secondDim);
    szzmax = dims.ubound(thirdDim);

    // Assert that the currently-split dimension is fully available
    assert(dims.lbound(reduce) == 0);
    // Now, reset the "max" of the reduced dimension to the "min"
    switch(reduce) {
        case firstDim:
            szxmax = szxmin; break;
        case secondDim:
            szymax = szymin; break;
        case thirdDim:
            szzmax = szzmin; break;
    }
    // Define the output Array
    Array<double,3> output(blitz::Range(szxmin,szxmax), 
            blitz::Range(szymin,szymax),
            blitz::Range(szzmin,szzmax));

    // Now, loop over the output points and sum up the overturning water column
    double zdiff = zgrid(zgrid.ubound()) - zgrid(zgrid.lbound());
    // Calculate a threshold value -- otherwise simple rounding error can
    // cause erroneous overturning

    /* As an ad-hoc measure, set the threshold of "significant" overturning to
       the maximum of:
       1) 1e-8 * the maximum rho value, or
       2) an overturning of 1%, extended over the entire domain, that is
       2% * (max-min) * LZ / NZ */
    double maxrho = pvmax(rho);
    double minrho = pvmin(rho); 
    double thresh = fmax(1e-8*maxrho,fabs(zdiff) * (maxrho-minrho) * 1e-2
            / (zgrid.ubound(firstDim) - zgrid.lbound(firstDim)));
    for (int i = szxmin; i <= szxmax; i++) {
        for (int j = szymin; j <= szymax; j++) {
            for (int k = szzmin; k <= szzmax; k++) {
                /* Now, build a zplus/zminus pair of ranges for the density
                   and z-level differencing, and do the reduction.  Most of
                   the code duplication here arises because Blitz doesn't like
                   sum() reduction over anything but the last logical dimension */
                if (reduce == firstDim) {
                    blitz::Range zplus(rho.lbound(firstDim)+1,rho.ubound(firstDim));
                    blitz::Range zminus(rho.lbound(firstDim),rho.ubound(firstDim)-1);
                    output(i,j,k) = fabs(sum(
                                where(zdiff * (rho(zplus,j,k) - rho(zminus,j,k)) > thresh,
                                    zgrid(zplus) - zgrid(zminus), 0)));
                } else if (reduce == secondDim) {
                    blitz::Range zplus(rho.lbound(secondDim)+1,rho.ubound(secondDim));
                    blitz::Range zminus(rho.lbound(secondDim),rho.ubound(secondDim)-1);
                    output(i,j,k) = fabs(sum(
                                where(zdiff * (rho(i,zplus,k) - rho(i,zminus,k)) > thresh,
                                    zgrid(zplus) - zgrid(zminus), 0)));
                } else if (reduce == thirdDim) {
                    blitz::Range zplus(rho.lbound(thirdDim)+1,rho.ubound(thirdDim));
                    blitz::Range zminus(rho.lbound(thirdDim),rho.ubound(thirdDim)-1);
                    output(i,j,k) = fabs(sum(
                                where(zdiff * (rho(i,j,zplus) - rho(i,j,zminus)) > thresh,
                                    zgrid(zplus) - zgrid(zminus), 0)));
                }
            }
        }
    }

    return output;
}

