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
