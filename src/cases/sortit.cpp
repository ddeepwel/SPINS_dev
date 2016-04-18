// Development case for array sorting in parallel, via quickselect

#include <blitz/array.h>
#include "../TArray.hpp"
#include <mpi.h>
#include <random/normal.h>
#include "../Par_util.hpp"
#include <stdio.h>
#include <unistd.h>
#include <vector>
#include <algorithm>
#include <stdlib.h>

using std::vector;

using namespace ranlib;

class kvpair {
   public:
   double key;
   double value;

   kvpair(double ikey = 0, double ival = 0):
      key(ikey), value(ival) {
   }

   bool operator< (const kvpair & rhs) const {
      return (key < rhs.key);
   }
   bool operator== (const kvpair & rhs) const {
      return (key == rhs.key);
   }
   operator double() const { 
      return key; 
   };
};

template <typename T> 
   T pivotselect(const vector<T> &arr, int lbound, int ubound) {
   // Selects a candidate local pivot via median-of-3, using randomly-
   // selected entries
   int low = rand()%(1+ubound-lbound)+lbound;
   int mid = rand()%(1+ubound-lbound)+lbound;
   int high = rand()%(1+ubound-lbound)+lbound;
   //T lkey = arr[lbound], mkey = arr[mid], ukey = arr[ubound];
   T lkey = arr[low], mkey = arr[mid], ukey = arr[high];

   /*int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);*/

   T outkey;

   if (lkey <= mkey) {
      if (mkey <= ukey) outkey = mkey;
      else outkey = ukey;
   } else {
      if (lkey <= ukey) outkey = lkey;
      else outkey = ukey;
   }
/*   if (!myrank) {
      fprintf(stdout,"ps: [%d] %.2f [%d] %.2f [%d] %.2f = %f\n",lbound,double(lkey),
            mid,double(mkey),ubound,double(ukey),double(outkey));
   }*/
   return outkey;

}

int gmedfindcount = 0;
double medlag = 0;
double partlag = 0;
double mpilag = 0;
double othlag = 0;
double sortlag = 0;

template <typename T>
   void partition(vector<T> &arr, int lbound, int ubound, T pivotval, int &pivotlt, int &pivoteq) {
   // Partitions the array arr[] using the "Dutch National Flag" algorithm, such
   // that the array is split into three parts, containing:
   // *) Values strictly less than the pivot
   // *) Values equal to the pivot
   // *) Values strictly greater than the pivot

   // Initialize output values to 0
   pivotlt = 0; // Number of values in this array <= the pivot
   pivoteq = 0; // Number of values strictly equal to the pivot

   int lptr = lbound;
   int mptr = lbound;
   int hptr = ubound;

   T dummy;

   /*int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);*/

   while (mptr <= hptr) {
      // All values strictly below lptr are less than the pivot,
      // all values >= lptr and strictly below mptr are equal to the pivot,
      // and all values strictly greater than hptr are greater than the pivot.
      // The values at mptr and hptr are unknown, so we work by examining mptr:

      /*if (!myrank) {
         for (int ii = lbound; ii <= ubound; ii++) {
            char fchar;
            if (ii == lptr) fchar = '*';
            else if (ii == mptr) fchar = '^';
            else if (ii == hptr) fchar = '#';
            else fchar = ' ';
            fprintf(stdout,"%c%+4.2f ",fchar,double(arr[ii]));
         } fprintf(stdout,"\n");
      }*/
      if (arr[mptr] < pivotval) {
         // This element belongs in the lower bucket, which expands to arr[lptr].
         // This value was previously inside the middle bucket, so it is placed
         // at arr[mptr] and we can increment both low and middle pointers
         if (lptr != mptr) {
            dummy = arr[lptr];
            arr[lptr] = arr[mptr];
            arr[mptr] = dummy;
         }
         lptr += 1; mptr += 1;
      } else if (arr[mptr] == pivotval) {
         // This element belongs in the middle bucket, so it does not need to move
         mptr += 1;
      } else { // arr[mptr] > pivotval
         // This element belongs in the upper bucket, freeing up some room;
         // put arr[hptr] into its proper location
         if (pivotval > arr[hptr]) {
            // arr[hptr] belongs in the low bucket
            dummy = arr[hptr];
            arr[hptr] = arr[mptr]; // put middle-pointed value in top of array
            if (mptr != lptr) {
               arr[mptr] = arr[lptr]; // Move first equal-pivot value to top of middle bucket
            }
            arr[lptr] = dummy;
            // update stack pointers
            mptr++; lptr++; hptr--;
         } else if (pivotval == arr[hptr]) {
            // arr[hptr] belongs on the top of the middle bucket
            dummy = arr[hptr];
            arr[hptr] = arr[mptr];
            arr[mptr] = dummy;
            mptr++; hptr--;
         } else {
            // arr[hptr] is already properly in the top bucket -- don't move
            hptr--;
         }
      }
   }
   // Now, we can define pivotlt and pivotle in terms of the ultimate
   // pointer locations:
   pivotlt = lptr - lbound;
   pivoteq = mptr - lptr;
}

template <typename T>
   T median(vector<T> &array, int lbound, int ubound, int idx = -1) {
   // Finds the median of the input array via serial quickselect
   // Addendum: seeks the target index

   int ilbound = lbound, iubound = ubound;
   int partlt = 0, parteq = 0;

   // The objective is to partition the array such that this
   // value is in its sorted place.
   int targetidx;
   if (idx > 0) {
      targetidx = lbound + (idx-1);
   } else {
      targetidx = (ubound + lbound + 1)/2;
   }

   /*int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);*/

   /*if (!myrank) {
      fprintf(stdout,"Median: lbound %d, ubound %d, idx %d\n",lbound,ubound,targetidx);
      fprintf(stdout,"Called with %d %d %d\n",lbound,ubound,idx);
      fprintf(stdout,"(%d + %d + 1)/2 = %d\n",ubound,lbound,(ubound+lbound+1)/2);
   }*/

   while (ilbound <= targetidx && iubound > targetidx) {
      T pivot = pivotselect(array, ilbound, iubound);
      //if (!myrank) fprintf(stdout,"Pivot is %f [%d-%d = %d]\n",double(pivot),ilbound,iubound,1+iubound-ilbound);

      /*if (!myrank) {
         for (int ii = lbound; ii <= ubound; ii++) {
            char mychar=' ';
            if (ii == ilbound) mychar = '$';
            if (ii == iubound) mychar = '#';
            if (array[ii] == pivot) mychar = '!';
            //fprintf(stdout,"%d: %f\n",ii,double(array[ii]));
            fprintf(stdout,"%c%+4.2f ",mychar,double(array[ii]));
         } 
         fprintf(stdout,"\n");
      }*/
      partition(array,ilbound,iubound,pivot,partlt,parteq);
      /*if (!myrank) {
         for (int ii = lbound; ii <= ubound; ii++) {
            char mychar=' ';
            if (ii == ilbound) mychar = '$';
            if (ii == iubound) mychar = '#';
            if (array[ii] == pivot) mychar = '!';
            //fprintf(stdout,"%d: %f\n",ii,double(array[ii]));
            fprintf(stdout,"%c%+4.2f ",mychar,double(array[ii]));
         } 
         fprintf(stdout,"\n");
      }*/

      //fprintf(stdout,"\n%d <, %d =\n",partlt,parteq);*/

      if (ilbound + partlt >= targetidx) {
         iubound = ilbound + partlt;
      } else if (ilbound + partlt + parteq <= targetidx) {
         ilbound = ilbound + partlt + parteq;
      } else {
         break;
      }
   }
   return array[targetidx];
}


class parsorter {
   public:
      MPI_Comm c, nestedcomm;
      int nproc;
      int myrank, mygrank;
      vector<kvpair> * tmparray;
      vector<int> asizes;
      int ASIZE;

      parsorter * nested;

      parsorter(int isize, MPI_Comm incomm):
         c(incomm), ASIZE(isize), tmparray(0) {
         MPI_Comm_rank(MPI_COMM_WORLD,&mygrank);
         MPI_Comm_rank(c,&myrank);
         MPI_Comm_size(c,&nproc);
         asizes.resize(nproc);
         ASIZE = isize;

         if (nproc > 1) {
            double now = MPI_Wtime();
            MPI_Allgather(&ASIZE,1,MPI_INT,&asizes[0],1,MPI_INT,c);
            int color = (myrank >= (nproc/2));
            MPI_Comm_split(c,color,myrank,&nestedcomm);
            mpilag += MPI_Wtime()-now;

            nested = new parsorter(isize,nestedcomm);
         } else {
            asizes[0] = ASIZE;
            nested = 0;
         }
      }

      ~parsorter() {
         if (nested) delete nested;
         if (tmparray) delete tmparray;
         //if (nproc > 1) MPI_Comm_free(&nestedcomm);
      }

      void sortkvpairs(vector<kvpair> * inarray, vector<kvpair> * tarray = 0) {
         vector<kvpair> &local = *inarray;

         if (nproc == 1) {
      //   {
            double now = MPI_Wtime();
            std::stable_sort(local.begin(),local.end());
            sortlag += MPI_Wtime()-now;
            if (nproc == 1) return;
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

            gmedfindcount++;
           
            now = MPI_Wtime(); 
            if (lubound >= llbound) {
               //kvpair med = median(local,llbound,lubound);
               kvpair med = pivotselect(local,llbound,lubound);
               lmedkey = med.key;
               lvalid = 1;
            } else {
               lmedkey = 999999;
               lvalid = 0;
            }
            MPI_Barrier(c);
            medlag += MPI_Wtime()-now;

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
            now = MPI_Wtime();
            MPI_Allgather(&lmedkey,1,MPI_DOUBLE,&allmedians[0],1,MPI_DOUBLE,c);

            // Count how many processors actually had a valid sub-array to take a median on
            int gvalid = 0;
            MPI_Allreduce(&lvalid,&gvalid,1,MPI_INT,MPI_SUM,c);
            mpilag += MPI_Wtime()-now;

            // Now, we know the array contains (allvalid) good entries, with the remainder set
            // to a very high value; that means we want to find the (allvalid)/2 element

            now = MPI_Wtime();
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
            partlag += MPI_Wtime()-now;
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
            now = MPI_Wtime();
            MPI_Allreduce(&literstat,&giterstat,2,MPI_INT,MPI_SUM,c);
            mpilag += MPI_Wtime()-now;

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

         now = MPI_Wtime();
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
         mpilag += MPI_Wtime()-now;

         now = MPI_Wtime();
         local = tmp;
         othlag += MPI_Wtime()-now;

         nested->sortkvpairs(&local,&tmp);

      }
};









int main(int argc, char ** argv) {
   MPI_Init(&argc,&argv);
   Normal<double> rnd(0,1);

   // Get comm size, rank
   int myrank, nproc;
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
   MPI_Comm_size(MPI_COMM_WORLD,&nproc);
   
   // Array sizing
   int ASIZE = -1;
   if (argc > 1) ASIZE=atoi(argv[1]);
   if (ASIZE < nproc) ASIZE=nproc;
   ASIZE = (ASIZE/nproc) + (myrank < (ASIZE % nproc));
   vector<kvpair> local(ASIZE);

   vector<int> asizes(nproc);
   MPI_Allgather(&ASIZE,1,MPI_INT,&asizes[0],1,MPI_INT,MPI_COMM_WORLD);

   // Initialize random variable
   rnd.seed(myrank);
   srand(myrank);

   double overall;
   for (int loop = 1; loop < 10; loop++) {
   // Initialize array (keys only)
   for (int ii = 0; ii < ASIZE; ii++) {
      local[ii].key = rnd.random();
      local[ii].value = myrank+0.001*ii;
   }

   double now = MPI_Wtime();
   parsorter sortme(ASIZE,MPI_COMM_WORLD);

   sortme.sortkvpairs(&local);
   overall += MPI_Wtime()-now;
   MPI_Barrier(MPI_COMM_WORLD);
   }


   //std::sort(lpart.begin(),lpart.end());

   // Print array
   for (int pp = 0; pp < nproc; pp++) {
      if (myrank == pp) {
         fprintf(stdout,"--%-2d [%d, %d, %-8.3g + %-8.3g + %-8.3g + %-8.3g + %-8.3g = %-8.3g]--\n",myrank,ASIZE,gmedfindcount,medlag,partlag,mpilag,sortlag,othlag,overall);
         /*for (int kk = 0; kk < ASIZE; kk++) {
            fprintf(stdout,"%+.3f ",double(local[kk]));
         }*/
         /*for (int kk = 0; kk < ASIZE; kk++) {
            fprintf(stdout,"%+.3f ",lpart[kk].key);
         }
         fprintf(stdout,"\n");
         for (int kk = 0; kk < ASIZE; kk++) {
            fprintf(stdout,"%+.3f ",lpart[kk].value);
         }*/
         //fprintf(stdout,"\n");
         usleep(100);
         MPI_Barrier(MPI_COMM_WORLD);
      } else {
         MPI_Barrier(MPI_COMM_WORLD);
      }
   }

   MPI_Finalize();
   return 0;
}

