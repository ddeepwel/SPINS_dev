#include <stdio.h>
#include <string.h>

/* Write the source code of the case to spinscase.cpp */
void WriteCaseFileSource(void)
{
    char* filename;
    if ( strcmp(casefilename, "cases/derivatives/derivatives.cpp") == 0 ) {
        filename = "derivatives.cpp";
    } else {
        filename = "spinscase.cpp";
    }
    FILE* fid;
    fid = fopen(filename,"w");
    fprintf(fid,"/* %s */\n%s",casefilename,casefilesource);
    fclose(fid);
}
