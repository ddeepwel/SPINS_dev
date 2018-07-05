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
