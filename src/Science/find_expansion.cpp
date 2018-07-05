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
