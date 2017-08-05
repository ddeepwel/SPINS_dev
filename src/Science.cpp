/* WARNING: Science Content!

   Implementation of various analysis routines */

#include "Science.hpp"
#include "math.h"
#include "Par_util.hpp"
#include "stdio.h"
#include "Split_reader.hpp"
#include "T_util.hpp"
#include "Parformer.hpp"
#include "Sorter.hpp"

// Marek's Overturning Diagnostic

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


// Read in a 2D-array from file and extend it to fill a full, 3D array in
// memory.  Unlike the following function, this uses the standard C storage
// order -- matlab uses the transpose of column-major ordering
void read_2d_restart(Array<double,3> & fillme, const char * filename, 
                  int Nx, int Ny) {

//   using blitz::ColumnMajorArray;
   using blitz::firstDim; using blitz::secondDim; using blitz::thirdDim;
   /* Get the local ranges we're interested in */
   blitz::Range xrange(fillme.lbound(firstDim),fillme.ubound(firstDim));
   blitz::Range zrange(fillme.lbound(thirdDim),fillme.ubound(thirdDim));
   
   /* Read the 2D slice from disk.  Matlab uses Column-Major array storage */

   blitz::GeneralArrayStorage<2> storage_order;
   blitz::Array<double,2> * sliced = 
      read_2d_slice<double>(filename,Nx,Ny,xrange,zrange,storage_order);

   /* Now, assign the slice to fill the 3D array */
   for(int y = fillme.lbound(secondDim); y <= fillme.ubound(secondDim); y++) {
      fillme(xrange,y,zrange) = (*sliced)(xrange,zrange);
   }
   delete sliced; 
}

// Read in a 2D file and interpret it as a 2D slice of a 3D array, for
// initialization with read-in-data from a program like MATLAB
void read_2d_slice(Array<double,3> & fillme, const char * filename, 
                  int Nx, int Ny) {

   using blitz::ColumnMajorArray;
   using blitz::firstDim; using blitz::secondDim; using blitz::thirdDim;
   /* Get the local ranges we're interested in */
   blitz::Range xrange(fillme.lbound(firstDim),fillme.ubound(firstDim));
   blitz::Range zrange(fillme.lbound(thirdDim),fillme.ubound(thirdDim));
   
   /* Read the 2D slice from disk.  Matlab uses Column-Major array storage */

   blitz::GeneralArrayStorage<2> storage_order = blitz::ColumnMajorArray<2>();
   blitz::Array<double,2> * sliced = 
      read_2d_slice<double>(filename,Nx,Ny,xrange,zrange,storage_order);

   /* Now, assign the slice to fill the 3D array */
   for(int y = fillme.lbound(secondDim); y <= fillme.ubound(secondDim); y++) {
      fillme(xrange,y,zrange) = (*sliced)(xrange,zrange);
   }
   delete sliced; 
}
   
// X-component of vorticity
void compute_vort_x(TArrayn::DTArray & vortx, TArrayn::DTArray & v, TArrayn::DTArray & w,
        TArrayn::Grad * gradient_op, const string * grid_type) {
    // Set-up
    S_EXP expan[3];
    assert(gradient_op);

    // Setup for dv/dz
    find_expansion(grid_type, expan, "v");
    gradient_op->setup_array(&v,expan[0],expan[1],expan[2]);
    // get dv/dz
    gradient_op->get_dz(&vortx,false);
    // Invert to get the negative
    vortx = (-1)*vortx;

    // Setup for dw/dy
    find_expansion(grid_type, expan, "w");
    gradient_op->setup_array(&w,expan[0],expan[1],expan[2]);
    // get dw/dy, and add to vortx
    gradient_op->get_dy(&vortx,true);
}

// Y-component of vorticity
void compute_vort_y(TArrayn::DTArray & vorty, TArrayn::DTArray & u, TArrayn::DTArray & w,
       TArrayn::Grad * gradient_op, const string * grid_type) {
    // Set-up
    S_EXP expan[3];
    assert(gradient_op);

    // Setup for dw/dx
    find_expansion(grid_type, expan, "w");
    gradient_op->setup_array(&w,expan[0],expan[1],expan[2]);
    // get dw/dx
    gradient_op->get_dx(&vorty,false);
    // Invert to get the negative
    vorty = (-1)*vorty;

    // Setup for du/dz
    find_expansion(grid_type, expan, "u");
    gradient_op->setup_array(&u,expan[0],expan[1],expan[2]);
    // get du/dz, and add to vorty
    gradient_op->get_dz(&vorty,true);
}

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

void compute_vorticity(TArrayn::DTArray & vortx, TArrayn::DTArray & vorty, TArrayn::DTArray & vortz,
        TArrayn::DTArray & u, TArrayn::DTArray & v, TArrayn::DTArray & w,
        TArrayn::Grad * gradient_op, const string * grid_type) {
    // compute each component
    compute_vort_x(vortx, v, w, gradient_op, grid_type);
    compute_vort_y(vorty, u, w, gradient_op, grid_type);
    compute_vort_z(vortz, u, v, gradient_op, grid_type);
}

bool compare_pairs( pair<double, double> a, pair<double, double> b ) {
	return a.first < b.first;
}

// Compute Background Potential Energy (BPE)
void compute_Background_PE(double & BPE_tot, TArrayn::DTArray & rho,
        int Nx, int Ny, int Nz, double Lx, double Ly, double Lz, double g,
        double rho_0, int iter, bool mapped, Array<double,1> hill) {
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
    static DTArray *quad3;     // quadrature weights
    static vector<double> sort_hill(Nx), sort_dx(Nx);
    double *local_vols =   new double[numprocs];
    double *local_vols_r = new double[numprocs];
    static double hill_max;
    static double hill_vol;
    static vector < pair<double, double> > height_width;

    // Stuff to do once at the beginning
    if (iter == 1) {
        // create array of voxels
        quad3 = alloc_array(Nx,Ny,Nz);
        *quad3 = (*get_quad_x())(ii)*(*get_quad_y())(jj)*(*get_quad_z())(kk);

        // adjust if mapped
        if ( mapped ) {
            *quad3 = (*quad3)*(Lz-hill(ii))/Lz;

            // information about the hill
            hill_vol = pssum(sum(hill*Ly*(*get_quad_x())));
            hill_max = psmax(max(hill));

            /* copy and sort the hill */
            double *hill_tmp =    new double[Nx];
            double *hill_tmptmp = new double[Nx];
            double *dx_tmp =    new double[Nx];
            double *dx_tmptmp = new double[Nx];
            // create temporary empty arrays for holding the grid and hill
            for (int II = 0; II < Nx; II++) {
                hill_tmptmp[II] = 0.;
                dx_tmptmp[II] = 0.;
            }
            // populate these arrays with their processor's information
            for (int II = hill.lbound(firstDim); II <= hill.ubound(firstDim); II++) {
                hill_tmptmp[II] = hill(II);
                dx_tmptmp[II] = (*get_quad_x())(II);
            }
            // share across processors
            MPI_Allreduce(hill_tmptmp, hill_tmp, Nx, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(dx_tmptmp,   dx_tmp,   Nx, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            // put into a double pair
            for (int II = 0; II < Nx; II++) {
                height_width.push_back(pair<double, double>(hill_tmp[II], dx_tmp[II]));
            }
            // sort the height_width pair according to the height
            sort(height_width.begin(), height_width.end(), compare_pairs);
            // put sorted pairs into separate vectors
            for (int II = 0; II < Nx; II++) {
                sort_hill[II] = height_width[II].first;
                sort_dx[II] = height_width[II].second;
            }
            // clean-up
            delete[] dx_tmp;
            delete[] dx_tmptmp;
            delete[] hill_tmp;
            delete[] hill_tmptmp;
        }
    }

    // Compute sorted rho
    sortarray(rho, *quad3, sort_rho, sort_quad);

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
        double tmp_vol = 0;
        tmpH = sort_hill[0];
        double dz_star;// = Lz/Nx/Ny/Nz;
        double Lx_star = 0;
        // check if base_vol will fill above the max hill height
        if ((base_vol + hill_vol)/Lx/Ly >= hill_max) {
            tmpH = (base_vol + hill_vol)/Lx/Ly;
        } else {
            // if it doesn't, find the depth where the fluid will fill to
            for (int II = 0; II < Nx-1; II++) {
                Lx_star += sort_dx[II];
                Area_star = Lx_star*Ly;
                dz_star = sort_hill[II+1] - sort_hill[II];
                tmp_vol += dz_star*Area_star;
                tmpH += dz_star;
                if (tmp_vol > base_vol) {
                    tmpH -= (tmp_vol - base_vol)/Area_star;
                    break;
                }
            }
        }
    }

    // Now compute BPE from the sorted rho field
    double dH = 0;
    double BPE = 0;
    for (int II = sort_rho.lbound(firstDim); II <= sort_rho.ubound(firstDim); II++) {
        for (int KK = sort_rho.lbound(thirdDim); KK <= sort_rho.ubound(thirdDim); KK++) {
            for (int JJ = sort_rho.lbound(secondDim); JJ <= sort_rho.ubound(secondDim); JJ++) {
                // Find the surface area of the water at depth tmpH
                if ( mapped && (tmpH < hill_max) ) {
                    double Lx_star = 0;
                    int LL = 0;
                    while ( ( sort_hill[LL] <= tmpH) && (LL < Nx) ) {
                        Lx_star += sort_dx[LL];
                        LL++;
                    }
                    Area_star = Lx_star*Ly;
                } else {
                    Area_star = Lx*Ly;
                }

                // spread volume over domain and compute BPE for that cell
                dH = sort_quad(II,JJ,KK)/Area_star;
                BPE += g*rho_0*(1+sort_rho(II,JJ,KK))*(tmpH - 0.5*dH)*dH*Area_star;
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

// Global arrays to store quadrature weights
Array<double,1> _quadw_x, _quadw_y, _quadw_z;

// Compute quadrature weights
void compute_quadweights(int szx, int szy, int szz, 
      double Lx, double Ly, double Lz,
      NSIntegrator::DIMTYPE DIM_X, NSIntegrator::DIMTYPE DIM_Y,
      NSIntegrator::DIMTYPE DIM_Z) {
   _quadw_x.resize(split_range(szx));
   _quadw_y.resize(szy); _quadw_z.resize(szz);
   if (DIM_X == NO_SLIP) {
      blitz::firstIndex ii;
      _quadw_x = 1;
      for (int k = 1; k <= (szx-2)/2; k++) {
         // From Trefethen, Spectral Methods in MATLAB
         // clenshaw-curtis quadrature weights
         _quadw_x -= 2*cos(2*k*M_PI*ii/(szx-1))/(4*k*k-1);
      }
      if ((szx%2))
         _quadw_x -= cos(M_PI*ii)/((szx-1)*(szx-1)-1);
      _quadw_x = 2*_quadw_x/(szx-1);
      if (_quadw_x.lbound(firstDim) == 0) {
         _quadw_x(0) = 1.0/(szx-1)/(szx-1);
      }
      if (_quadw_x.ubound(firstDim) == (szx-1)) {
         _quadw_x(szx-1) = 1.0/(szx-1)/(szx-1);
      }
      _quadw_x *= Lx/2;
   } else {
      // Trapezoid rule
      _quadw_x = Lx/szx;
   }
   if (DIM_Y == NO_SLIP) {
      blitz::firstIndex ii;
      _quadw_y = 1;
      for (int k = 1; k <= (szy-2)/2; k++) {
         // From Trefethen, Spectral Methods in MATLAB
         // clenshaw-curtis quadrature weights
         _quadw_y -= 2*cos(2*k*(M_PI*(ii)/(szy-1)))/(4*k*k-1);
      }
      if ((szy%2))
         _quadw_y -= cos(M_PI*ii)/((szy-1)*(szy-1)-1);
      _quadw_y = 2*_quadw_y/(szy-1);
      _quadw_y(0) = 1.0/(szy-1)/(szy-1);
      _quadw_y(szy-1) = 1.0/(szy-1)/(szy-1);
      _quadw_y *= Ly/2;
   } else {
      // Trapezoid rule
      _quadw_y = Ly/szy;
   }
   if (DIM_Z == NO_SLIP) {
      blitz::firstIndex ii;
      _quadw_z = 1;
      for (int k = 1; k <= (szz-2)/2; k++) {
         // From Trefethen, Spectral Methods in MATLAB
         // clenshaw-curtis quadrature weights
         _quadw_z -= 2*cos(2*k*(M_PI*(ii)/(szz-1)))/(4*k*k-1);
      }
      if ((szz%2))
         _quadw_z -= cos(M_PI*ii)/((szz-1)*(szz-1)-1);
      _quadw_z = 2*_quadw_z/(szz-1);
      _quadw_z(0) = 1.0/(szz-1)/(szz-1);
      _quadw_z(szz-1) = 1.0/(szz-1)/(szz-1);
      _quadw_z *= Lz/2;
   } else {
      // Trapezoid rule
      _quadw_z = Lz/szz;
   }
}
const blitz::Array<double,1> * get_quad_x() { 
   // Check whether the quad weight has been initialized
   if (_quadw_x.length(firstDim) <= 0) {
      assert(0 && "Error: quadrature weights were not initalized before use");
   }
   return &_quadw_x;
}
const blitz::Array<double,1> * get_quad_y() {
   if (_quadw_y.length(firstDim) <= 0) {
      assert(0 && "Error: quadrature weights were not initalized before use");
   }
   return &_quadw_y;
}
const blitz::Array<double,1> * get_quad_z() {
   if (_quadw_z.length(firstDim) <= 0) {
      assert(0 && "Error: quadrature weights were not initalized before use");
   }
   return &_quadw_z;
}

// function to parse the expansion types
void find_expansion(const string * grid_type, S_EXP * expan, string deriv_filename) {
    const int x_ind = 0;
    const int y_ind = 1;
    const int z_ind = 2;
    string prev_deriv, base_field;

    // check if field is a derivative field
    bool input_deriv = false;        // assume it's not a derivative
    int var_len = deriv_filename.length();
    if ( var_len > 2 ) {
        if ( deriv_filename.substr(var_len-2,1) == "_" ) {
            // if second last char is an underscore then its a derivative field
            input_deriv = true;
            prev_deriv = deriv_filename.substr(var_len-1,1);  // the completed derivative
            base_field = deriv_filename.substr(0,var_len-2);  // the differentiated field
        }
    }

    for ( int nn = 0; nn <= 2; nn++ ) {
        if      (grid_type[nn] == "FOURIER") { expan[nn] = FOURIER; }
        else if (grid_type[nn] == "NO_SLIP") { expan[nn] = CHEBY; }
        else if (grid_type[nn] == "FREE_SLIP") {
            // setup for a first derivative
            if ( deriv_filename == "u" or base_field == "u" ) {
                if      ( nn == x_ind ) { expan[nn] = SINE; }
                else if ( nn == y_ind ) { expan[nn] = COSINE; }
                else if ( nn == z_ind ) { expan[nn] = COSINE; }
            }
            else if ( deriv_filename == "v" or base_field == "v" ) {
                if      ( nn == x_ind ) { expan[nn] = COSINE; }
                else if ( nn == y_ind ) { expan[nn] = SINE; }
                else if ( nn == z_ind ) { expan[nn] = COSINE; }
            }
            else if ( deriv_filename == "w" or base_field == "w") {
                if      ( nn == x_ind ) { expan[nn] = COSINE; }
                else if ( nn == y_ind ) { expan[nn] = COSINE; }
                else if ( nn == z_ind ) { expan[nn] = SINE; }
            }
            else {
                if      ( nn == x_ind ) { expan[nn] = COSINE; }
                else if ( nn == y_ind ) { expan[nn] = COSINE; }
                else if ( nn == z_ind ) { expan[nn] = COSINE; }
            }
        }
    }

    // adjust if input field is a derivative field
    if ( input_deriv == true ) {
        if      ( prev_deriv == "x" ) { expan[x_ind] = swap_trig(expan[x_ind]); }
        else if ( prev_deriv == "y" ) { expan[y_ind] = swap_trig(expan[y_ind]); }
        else if ( prev_deriv == "z" ) { expan[z_ind] = swap_trig(expan[z_ind]); }
    }
}
// function to switch trig functions
S_EXP swap_trig( S_EXP the_exp ) {
    if ( the_exp == SINE ) {
        return COSINE; }
    else if ( the_exp == COSINE ) {
        return SINE; }
    else if ( the_exp == FOURIER ) {
        return FOURIER; }
    else if ( the_exp == CHEBY ) {
        return CHEBY; }
    else {
        MPI_Finalize(); exit(1); // stop
    }
}
