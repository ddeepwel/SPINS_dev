#include "../BaseCase.hpp"

// Write plot time information
void BaseCase::write_plot_times(double time, double clock_time, double comp_duration,
        double avg_write_time, int plot_number, bool restarting) {
    if (master()) {
        // in log file
        fprintf(stdout,"*Write time: %.6g. Average write time: %.6g.\n",
                clock_time - comp_duration, avg_write_time);
        // track in a file
        FILE * plottimes_file = fopen("plot_times.txt","a");
        assert(plottimes_file);
        if ( plot_number==get_restart_sequence()+1 and !restarting )
            fprintf(plottimes_file,"Output number, Simulation time (s), "
                    "Write time (s), Average write time (s)\n");
        fprintf(plottimes_file,"%d, %.17f, %.12g, %.12g\n",
                plot_number, time, clock_time - comp_duration, avg_write_time);
        fclose(plottimes_file);
    }
}
