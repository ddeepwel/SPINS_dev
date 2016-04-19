#ifndef SORTER_HPP
#define SORTER_HPP 1

namespace Sorter {

   template <typename T>
      T pivotselect(const vector<T> &arr, int lbound, int ubound);

   template <typename T>
      void partition(vector <T> &arr, int lbound, int ubound, T pivotval, int &pivotlt, int &pivoteq);

   template <typename T>
      T median(vector <T> &array, int lbound, int ubound, int idx = -1);


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

   class parsorter {
      public: 
         MPI_Comm c, nestedcomm;
         int nproc, myrank, mygrank;
         vector<kvpair> * tmparray;
         vector<int> asizes;
         int ASIZE;

         parsorter * nested;

         parsorter(int isize, MPI_Comm incomm);
         ~parsorter();
         void sortkvpairs(vector<kvpair> * inarray, vector<kvpair> * tarray = 0);
   };


#include "Sorter_impl.cc"
};

// Outside the namespace

void sortarray(blitz::Array<double,3> const &keys, blitz::Array<double,3> const &values, 
               blitz::Array<double,3> & sortkeys,
               blitz::Array<double,3> & sortvals,
               MPI_Comm c = MPI_COMM_WORLD);

#endif
