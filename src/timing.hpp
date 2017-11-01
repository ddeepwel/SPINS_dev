/* timing.hpp -- header file for self-timing code */

#ifndef TIMING_HPP
#define TIMING_HPP 1  // Include this header file at most once

/* The timing infrastructure works with a stack-based model.

   At portions of the code to be instrumented, call timing_push("name")
   with some (short but) memortable name for the code segment; after the
   end of the section call timing_pop().  This underlying code will build
   the corresponding call-time-tree and measure (via MPI_Wtime) the time
   between the timing_push and the corresponding timing_pop.

   This model does require that timing be called in a predictable way, such
   that timing_pop is called before any possible exit from the code path.*/

#ifndef TIMING_ENABLE // If timing code is not enabled, define the respective functions as no-ops 

#define timing_push(x)
#define timing_pop()

#else
void timing_push(const char * name);
void timing_pop();
#endif

void timing_stack_report();



#endif
