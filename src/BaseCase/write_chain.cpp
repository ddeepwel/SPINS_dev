#include "../BaseCase.hpp"

/* write out vertical chain of data */
void BaseCase::write_chain(const char *filename, DTArray & val, int Iout, int Jout, double time) {
    FILE *fid=fopen(filename,"a");
    if (fid == NULL) {
        fprintf(stderr,"Unable to open %s for writing\n",filename);
        MPI_Finalize(); exit(1);
    }
    fprintf(fid,"%g",time);
    for (int ki=0; ki<size_z(); ki++) fprintf(fid," %g",val(Iout,Jout,ki));
    fprintf(fid,"\n");
    fclose(fid);
}
