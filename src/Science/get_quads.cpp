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
