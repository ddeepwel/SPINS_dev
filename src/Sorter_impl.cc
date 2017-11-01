#ifndef SORTER_HPP
   #error "Sorter_impl.cc may only be #included in Sorter.hpp"
#endif

#ifndef SORTER_IMPL_CC
#define SORTER_IMPL_CC 1

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
   T median(vector<T> &array, int lbound, int ubound, int idx) {
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
#endif
