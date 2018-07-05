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

bool compare_pairs( pair<double, double> a, pair<double, double> b ) {
	return a.first < b.first;
}

// Compute Background Potential Energy (BPE)
void compute_Background_PE(double & BPE_tot, TArrayn::DTArray & rho, TArrayn::DTArray & quad3,
        int Nx, int Ny, int Nz, double Lx, double Ly, double Lz, double g,
        double rho_0, int iter, bool dimensional_rho, bool mapped, Array<double,1> hill) {
    // Tensor variables for indexing
    blitz::firstIndex ii;
    blitz::secondIndex jj;
    blitz::thirdIndex kk;
    // Set mpi parameters
    int myrank, numprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    // Arrays for sorting
    static Array<double,3> sort_rho, sort_quad;
    static vector<double> sort_hill(Nx), sort_dx(Nx);
    static vector<double> Lx_partsum(Nx), hillvol_partsum(Nx);
    double *local_vols =   new double[numprocs];
    double *local_vols_r = new double[numprocs];
    static double hill_max;
    static double hill_vol;
    static vector < pair<double, double> > height_width(Nx);

    // Stuff to do once at the beginning
    if (iter == 0) {
        // adjust if mapped
        if ( mapped ) {
            // information about the hill
            hill_vol = pssum(sum(hill*Ly*(*get_quad_x())));
            hill_max = psmax(max(hill));

            /* copy and sort the hill */
            double *hill_tmp  = new double[Nx];
            double *hill_tmp2 = new double[Nx];
            double *dx_tmp  =   new double[Nx];
            double *dx_tmp2 =   new double[Nx];
            // create temporary empty arrays for holding the grid and hill
            for (int II = 0; II < Nx; II++) {
                if ( (II >= hill.lbound(firstDim)) and (II <= hill.ubound(firstDim)) ) {
                    // populate these arrays with their processor's information
                    hill_tmp2[II] = hill(II);
                    dx_tmp2[II] = (*get_quad_x())(II);
                } else {
                    // else set to zero
                    hill_tmp2[II] = 0.;
                    dx_tmp2[II] = 0.;
                }
            }
            // share across processors
            MPI_Allreduce(hill_tmp2, hill_tmp, Nx, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(dx_tmp2,   dx_tmp,   Nx, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            // put into a double pair
            for (int II = 0; II < Nx; II++) {
                height_width[II].first  = hill_tmp[II];
                height_width[II].second = dx_tmp[II];
            }
            // sort the height_width pair according to the height
            sort(height_width.begin(), height_width.end(), compare_pairs);
            // put sorted pairs into separate vectors
            for (int II = 0; II < Nx; II++) {
                sort_hill[II] = height_width[II].first;
                sort_dx[II] = height_width[II].second;
            }
            // Compute cumulative sums of hill volume (reusing Lx_partsum)
            for (int II = 0; II < Nx; II++) {
                Lx_partsum[II] = sort_hill[II]*sort_dx[II]; 
            }
            std::partial_sum(Lx_partsum.begin(), Lx_partsum.end(), hillvol_partsum.begin());
            // Compute cumulative sums of dxs
            std::partial_sum(sort_dx.begin(), sort_dx.end(), Lx_partsum.begin());
            // clean-up
            delete[] dx_tmp;
            delete[] dx_tmp2;
            delete[] hill_tmp;
            delete[] hill_tmp2;
        }
    }

    // Compute sorted rho
    sortarray(rho, quad3, sort_rho, sort_quad);

    // Volume of memory-local portion of sorted array
    double local_vol = sum(sort_quad);

    // Share local volume information with other processors
    // ie. processor 0 has lightest fluid
    for (int II = 0; II < numprocs; II++) {
        local_vols[II] = (II <= myrank) ? local_vol : 0;
    }
    MPI_Allreduce(local_vols, local_vols_r, numprocs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // The total volume of fluid with greater or equal density
    // to the lightest element on the local domain
    double base_vol = local_vols_r[myrank];

    // Find the depth at which base_vol will fill to
    double tmpH, Area_star;
    if ( !mapped ) {
        tmpH = base_vol/Lx/Ly;
    } else {
        // check if base_vol will fill above the max hill height
        if ((base_vol + hill_vol)/Lx/Ly >= hill_max) {
            tmpH = (base_vol + hill_vol)/Lx/Ly;
        } else {
            // if it doesn't, find the depth where the fluid will fill to
            int II = 0;
            while ( (base_vol > Lx_partsum[II]*sort_hill[II] - hillvol_partsum[II]) and (II<Nx-1) ) {
                II++;
            }
            assert(II>0 && "Something bad happened leading up to the tmpH calculation.");
            // now subtract off the bit we went over by (draw yourself a picture, it works)
            tmpH = sort_hill[II] - (sort_hill[II]*Lx_partsum[II] - base_vol - hillvol_partsum[II])/Lx_partsum[II-1];
        }
    }

    // Now compute BPE from the sorted rho field
    double dH = 0;
    double BPE = 0;
    int LL = Nx-1;
    for (int II = sort_rho.lbound(firstDim); II <= sort_rho.ubound(firstDim); II++) {
        for (int KK = sort_rho.lbound(thirdDim); KK <= sort_rho.ubound(thirdDim); KK++) {
            for (int JJ = sort_rho.lbound(secondDim); JJ <= sort_rho.ubound(secondDim); JJ++) {
                // Find the surface area of the water at depth tmpH
                if ( mapped && (tmpH < hill_max) ) {
                    while ( (sort_hill[LL] > tmpH) && (LL > 0) ) {
                        LL--;
                    }
                    Area_star = Lx_partsum[LL]*Ly;
                } else {
                    Area_star = Lx*Ly;
                }

                // spread volume over domain and compute BPE for that cell
                dH = sort_quad(II,JJ,KK)/Area_star;
                if (dimensional_rho) {
                    // rho is dimensional density
                    BPE += g*sort_rho(II,JJ,KK)*(tmpH - 0.5*dH)*dH*Area_star;
                } else {
                    // rho is density anomaly
                    BPE += g*rho_0*(1+sort_rho(II,JJ,KK))*(tmpH - 0.5*dH)*dH*Area_star;
                }
                tmpH -= dH;
            }
        }
    }

    // Share BPE with other processors
    MPI_Allreduce(&BPE, &BPE_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // delete variables
    delete[] local_vols;
    delete[] local_vols_r;
}
