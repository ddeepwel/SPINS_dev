#include "../BaseCase.hpp"

/* Change dump log file for successful completion*/
void BaseCase::successful_dump(int plot_number, double final_time, double plot_interval) {
    if (master() && (plot_number == final_time/plot_interval)) {
        // Write the dump time to a text file for restart purposes
        FILE * dump_file;
        dump_file = fopen("dump_time.txt","w");
        assert(dump_file);
        fprintf(dump_file,"The dump 'time' was:\n%.12g\n", 2*final_time);
        fprintf(dump_file,"The dump index was:\n%d\n", plot_number);
        fclose(dump_file);
    }
}

