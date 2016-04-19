// Development case for array sorting in parallel, via quickselect

#include <blitz/array.h>
#include "TArray.hpp"
#include <mpi.h>
#include <random/normal.h>
#include "Par_util.hpp"
#include <stdio.h>
#include <unistd.h>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <blitz/tinyvec-et.h>

#include "Sorter.hpp"

using std::vector;

using namespace ranlib;

namespace Sorter {

   parsorter::parsorter(int isize, MPI_Comm incomm):
      c(incomm), ASIZE(isize), tmparray(0) {
         MPI_Comm_rank(MPI_COMM_WORLD,&mygrank);
         MPI_Comm_rank(c,&myrank);
         MPI_Comm_size(c,&nproc);
         asizes.resize(nproc);
         ASIZE = isize;

         if (nproc > 1) {
//            double now = MPI_Wtime();
            MPI_Allgather(&ASIZE,1,MPI_INT,&asizes[0],1,MPI_INT,c);
            int color = (myrank >= (nproc/2));
            MPI_Comm_split(c,color,myrank,&nestedcomm);
//            mpilag += MPI_Wtime()-now;

            nested = new parsorter(isize,nestedcomm);
         } else {
            asizes[0] = ASIZE;
            nested = 0;
         }
   }

   parsorter::~parsorter() {
      if (nested) delete nested;
      if (tmparray) delete tmparray;
      //if (nproc > 1) MPI_Comm_free(&nestedcomm);
   }

   void parsorter::sortkvpairs(vector<kvpair> * inarray, vector<kvpair> * tarray) {
      vector<kvpair> &local = *inarray;

      if (nproc == 1) {
         //double now = MPI_Wtime();
         std::stable_sort(local.begin(),local.end());
         //sortlag += MPI_Wtime()-now;
         return;
      }

      if (tarray == 0) {
         if (tmparray == 0) tmparray = new vector<kvpair>(ASIZE);
         tarray = tmparray;
      } 
      vector<kvpair> &tmp = *tarray;

      // Perform a global search for the median;
      int gtoosmall = 0; // Number of global entries confirmed to be too small
      int gtoobig = 0; // Number of global entries confirmed to be too big
      typedef struct {int lt; int eq;} iterstats;
      iterstats giterstat = {0,0}, literstat = {0,0}; // Iteration statistics: number now less-than and equal-to

      int gtarget = 1; // Base-one index of the target, if we had a globally-sorted array
      for (int pp = 0; pp < (nproc/2); pp++) {
         // We want to put the median just after the halfway-processor
         gtarget += asizes[pp];
      }

      double gmedkey = 0; // Key of the median element

      int llbound = 0; 
      int lubound = ASIZE-1;

      double lmedkey = 0.0; // Key of the local median
      int lvalid = 0;

      vector<double> allmedians(nproc);
      double now;

      for(;;) {

         //gmedfindcount++;

//         now = MPI_Wtime(); 
         if (lubound >= llbound) {
            //kvpair med = median(local,llbound,lubound);
            kvpair med = pivotselect(local,llbound,lubound);
            lmedkey = med.key;
            lvalid = 1;
         } else {
            lmedkey = 999999;
            lvalid = 0;
         }
//         MPI_Barrier(c);
//         medlag += MPI_Wtime()-now;

         /*for (int pp = 0; pp < nproc; pp++) {
           if (myrank == pp) {
           fprintf(stdout,"%d: Local median (%d-%d) is %f\n",mygrank,llbound,lubound,lmedkey);
           usleep(100);
           MPI_Barrier(c);
           } else {
           MPI_Barrier(c);
           }
           }*/

         // Collect each of the candidate median values
//         now = MPI_Wtime();
         MPI_Allgather(&lmedkey,1,MPI_DOUBLE,&allmedians[0],1,MPI_DOUBLE,c);

         // Count how many processors actually had a valid sub-array to take a median on
         int gvalid = 0;
         MPI_Allreduce(&lvalid,&gvalid,1,MPI_INT,MPI_SUM,c);
//         mpilag += MPI_Wtime()-now;

         // Now, we know the array contains (allvalid) good entries, with the remainder set
         // to a very high value; that means we want to find the (allvalid)/2 element

//         now = MPI_Wtime();
         double est_pivot = median(allmedians,0,nproc-1,gvalid/2+1);
         /*if (myrank == 0) {
           fprintf(stdout,"Candidate median (of %d) is %f\n",gvalid,est_pivot);
           }*/

         if (lubound >= llbound) {
            // Pivot the array about the estimated median
            partition(local,llbound,lubound,kvpair(est_pivot),literstat.lt,literstat.eq);
         } else {
            literstat.lt = 0; literstat.eq = 0;
         }
//         partlag += MPI_Wtime()-now;
         /*for (int pp = 0; pp < nproc; pp++) {
           if (myrank == pp) {
           fprintf(stdout,"%d: %d below, %d equal %d above\n",pp,literstat.lt,literstat.eq,(lubound-llbound-literstat.lt-literstat.eq+1));
           usleep(100);
           MPI_Barrier(c);
           } else {
           MPI_Barrier(c);
           }
           }*/

         // Accumulate the local statistics to find out how many elements
         // of the global array were less than or equal to the candidate median
//         now = MPI_Wtime();
         MPI_Allreduce(&literstat,&giterstat,2,MPI_INT,MPI_SUM,c);
//         mpilag += MPI_Wtime()-now;

         /*if (myrank == 0) {
           fprintf(stdout,"Global: %d below, %d equal (+%d/%d)\n",giterstat.lt,giterstat.eq,gtoosmall,gtarget);
           }*/


         if (gtoosmall + giterstat.lt + giterstat.eq < gtarget) {
            // Candidate median was too small, so update bounds to exclude values
            // less than or equal to the candidate
            gtoosmall += giterstat.lt + giterstat.eq;
            llbound += literstat.lt + literstat.eq;

            /*if (myrank == 0) {
              fprintf(stdout,"Candidate median too small (%d)\n",gtoosmall);
              }*/
         } else if (gtoosmall + giterstat.lt >= gtarget) {
            // Candidate median was too large, so update bounds to exclude
            // values greater than or equal to the candidate
            lubound = llbound + literstat.lt - 1;
            /*if (myrank == 0) {
              fprintf(stdout,"Candidate median too large\n");
              }*/
         } else {
            gtoosmall += giterstat.lt;
            llbound += literstat.lt;
            /*if (myrank == 0) {
              fprintf(stdout,"Candidate median is just right\n");
              }*/
            break;
         }


      }

      // Now, accumulate per-processor stats on how many entries are lower than and equal to
      // the global median

      literstat.lt = llbound;
      iterstats procstats[nproc]; 

//      now = MPI_Wtime();
      MPI_Allgather(&literstat,2,MPI_INT,procstats,2,MPI_INT,c);

      // Number of at-median entries that are distributed with below-median entries
      // for even-splitting purposes
      int medbelow = gtarget - gtoosmall - 1;

      /*if (myrank == 0) {
        fprintf(stdout,"%d median entries to reclassify as less-than\n",medbelow);
        }*/
      for (int pp = 0; pp < nproc; pp++) {
         if (procstats[pp].eq >= medbelow) {
            procstats[pp].eq -= medbelow;
            procstats[pp].lt += medbelow;
            medbelow = 0;
            break;
         } else if (procstats[pp].eq > 0) {
            medbelow -= procstats[pp].eq;
            procstats[pp].lt += procstats[pp].eq;
            procstats[pp].eq = 0;
         }
      }



      int receivefrom[nproc];
      int sendto[nproc];
      for (int ii = 0; ii < nproc; ii++) {
         receivefrom[ii] = 0;
         sendto[ii] = 0;
      }
      // Bottom-half of processors 
      if (myrank < nproc/2) {

         // Calculate receive-counts from the top half
         int fillbefore = 0; // Number of spots that must be filled before this processor's allotment
         int fillme = ASIZE-procstats[myrank].lt; // number of spots that must be filled here
         //recfrom[myrank] = procstats[myrank].lt;
         //sendto[myrank] = procstats[myrank].lt;
         for (int pp = 0; pp < myrank; pp++) {
            fillbefore += (asizes[pp]-procstats[pp].lt);
         }
         for (int pp = nproc/2; pp < nproc; pp++) {
            int sendfromhere = procstats[pp].lt; // Number of entries available at this processor
            fillbefore -= sendfromhere;
            if (fillbefore < 0) {
               if (-fillbefore > fillme) {
                  receivefrom[pp] = fillme;
                  // If the array size does not change, sends and receives are symmetric
                  sendto[pp] = receivefrom[pp];
                  break;
               } else {
                  receivefrom[pp] = -fillbefore;
                  sendto[pp] = receivefrom[pp];
                  fillme -= (-fillbefore);
                  fillbefore = 0;
               }
            }
         }
      } else {
         int fillbefore = 0;
         int fillme = procstats[myrank].lt;
         for (int pp = nproc/2; pp < myrank; pp++) {
            fillbefore += procstats[pp].lt;
         }
         for (int pp = 0; pp < nproc/2; pp++) {
            int sendfromhere = asizes[pp]-procstats[pp].lt;
            fillbefore -= sendfromhere;
            if (fillbefore < 0) {
               if (-fillbefore > fillme) {
                  receivefrom[pp] = fillme;
                  sendto[pp] = receivefrom[pp];
                  break;
               } else {
                  receivefrom[pp] = -fillbefore;
                  sendto[pp] = receivefrom[pp];
                  fillme -= (-fillbefore);
                  fillbefore = 0;
               }
            }
         }
      }

      // Print array
      /*for (int pp = 0; pp < nproc; pp++) {
        if (myrank == pp) {
      //fprintf(stdout,"%d: %d+%d+%d\n",pp,literstat.lt,literstat.eq,ASIZE-(literstat.lt+literstat.eq));
      fprintf(stdout,"%3d: ",mygrank);
      for (int qq = 0; qq < nproc; qq++) {
      fprintf(stdout,"%3d ",receivefrom[qq]);
      }
      fprintf(stdout,"\n");
      usleep(100);
      MPI_Barrier(c);
      } else {
      MPI_Barrier(c);
      }
      }*/

      MPI_Request sendreqs[nproc], recreqs[nproc];
      int snum = 0, rnum = 0;

      // Copy the portion of our array that's kept from local to laprt
      if (myrank < nproc/2) {
         // Copy the first part
         memcpy(&tmp[0],&local[0],procstats[myrank].lt*sizeof(kvpair));
      } else {
         int keepnum = ASIZE-procstats[myrank].lt;
         memcpy(&tmp[procstats[myrank].lt],&local[procstats[myrank].lt],keepnum*sizeof(kvpair));
      }

      // Rolling offset used for send/receive caluclations
      int soffset = 0, roffset = 0;
      if (myrank < nproc/2) {
         soffset = procstats[myrank].lt;
         roffset = procstats[myrank].lt;
      }

      // Build non-blocking send/receives
      for (int pp = 0; pp < nproc; pp++) {
         if (sendto[pp] > 0) {
            MPI_Isend(&local[soffset],2*sendto[pp],MPI_DOUBLE,pp,0,c,&sendreqs[snum]);
            soffset += sendto[pp];
            snum++;
         }
         if (receivefrom[pp] > 0) {
            MPI_Irecv(&tmp[roffset],2*receivefrom[pp],MPI_DOUBLE,pp,0,c,&recreqs[rnum]);
            roffset += receivefrom[pp];
            rnum++;
         }
      }

      // Move the request handles to a single array
      MPI_Request areqs[2*nproc];
      MPI_Status astats[2*nproc];

      for (int kk = 0; kk < snum; kk++) areqs[kk] = sendreqs[kk];
      for (int kk = 0; kk < rnum; kk++) areqs[snum+kk] = recreqs[kk];

      MPI_Waitall(snum+rnum,areqs,astats);
//      mpilag += MPI_Wtime()-now;

//      now = MPI_Wtime();
      local = tmp;
//      othlag += MPI_Wtime()-now;

      nested->sortkvpairs(&local,&tmp);

   }
}// end namespace


// Outside of namespace

using TArrayn::DTArray;
using blitz::Array;
using namespace Sorter;

void sortarray(DTArray const &keys, DTArray const &vals,
               Array<double,3> & sortkeys, Array<double,3> & sortvals,
               MPI_Comm c) {
   /* Interface to sorting: take input 3D key/value arrays, sort them,
      and return the results in 1D output sortkeys/sortvalues; the output
      arrays are formally 3D for write_array purposes, but they will
      be of dimension <Nx,1,1> */

   // Static variables hold repeated state.  This will be double-checked
   // against the input arrays and communicator; in most use-cases we
   // should only call this function with a single unique array shape

   static parsorter * sortobj = 0;
   static MPI_Comm lastcomm = MPI_COMM_WORLD;
   static TinyVector<int,3> keybase;
   static int asize = 0, abase = 0;
   static bool warned = false;
   static vector<kvpair> pairs; // Array to hold key-value pairs for sorting
   bool initialize = false;

   if (!sortobj) {
      // First call -- initialize the above static arrays

      initialize = true;
   } else if (asize != keys.numElements() || lastcomm != c || 
              any(keybase != keys.base())) {
      // Also initialize, but complain about it

      delete sortobj;
      if (!warned) {
         fprintf(stderr,"WARNING: sortarray() called with multiple array shapes in a single program.\n"
                        "         this will result in poor performance.\n");
         warned = true;
      }

      initialize = true;
   }

   if (initialize) {
      asize = keys.numElements();
      pairs.resize(asize);
      lastcomm = c;
      keybase = keys.base();

      sortobj = new parsorter(asize,lastcomm);

      // Calculate the base index of our output array, since it is still
      // a single, global 1D array.  This is accomplished via MPI_Exscan,
      // which works like a cumulative sum exclusive of "this" element.

      int myrank;
      MPI_Comm_rank(lastcomm,&myrank);
      
      MPI_Exscan(&asize,&abase,1,MPI_INT,MPI_SUM,lastcomm);

      if (myrank == 0) abase = 0; // This is left undefined by Exscan

      // Debugging printout -- to be deleted
      fprintf(stdout,"Proc %d: local size %d elements, local base %d\n",myrank,asize,abase);
   } 

   for (int ii = 0; ii < asize; ii++) {
      // Copy keys and values from input array to pair-vector.  Since we're sorting this,
      // we don't need to care about the actual index-to-data mapping for the keys/values,
      // and we can copy in memory order
      pairs[ii].key = keys.data()[ii];
      pairs[ii].value = vals.data()[ii];
   }

   // Sort
   sortobj->sortkvpairs(&pairs);


   // If necessary, resize and rebase the sorted arrays
   if (sortkeys.numElements() != asize) sortkeys.resize(asize,1,1);
   if (sortvals.numElements() != asize) sortvals.resize(asize,1,1);

   // And rebase the sorted arrays if necessary
   if (sortkeys.lbound(0) != abase) {
      sortkeys.reindexSelf(TinyVector<int,3>(abase,0,0));
   }
   if (sortvals.lbound(0) != abase) {
      sortvals.reindexSelf(TinyVector<int,3>(abase,0,0));
   }

   // Copy the sorted results to sortkeys/sortvals
   for (int ii = 0; ii < asize; ii++) {
      sortkeys.data()[ii] = pairs[ii].key;
      sortvals.data()[ii] = pairs[ii].value;
   }


}




