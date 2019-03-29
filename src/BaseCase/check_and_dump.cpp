#include "../BaseCase.hpp"

/* Check and dump */
void BaseCase::check_and_dump(double clock_time, double real_start_time,
        double compute_time, double sim_time, double avg_write_time, int plot_number, int iter,
        DTArray & u, DTArray & v, DTArray & w, vector<DTArray *> & tracer) {
    if (compute_time > 0) {
        // default is to not dump variables
        int do_dump = 0;

        // check if close to end of compute time
        if (master()) {
            double total_run_time = clock_time - real_start_time;
            double avg_clk_step = total_run_time/(iter+1); // +1 to accound for initial iter=0
            if (compute_time - total_run_time < 5*(avg_write_time + avg_clk_step)) {
                do_dump = 1; // true
            }
        }

        // Broadcast to all processors whether to dump or not
        MPI_Bcast(&do_dump, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // Dump variables if close to end of compute time
        if (do_dump == 1) {
            if (master()){
                fprintf(stdout,"Too close to final time, dumping!\n");
            }
            write_variables(u, v, w, tracer);

            // Write the dump time to a text file for restart purposes
            if (master()){
                FILE * dump_file;
                dump_file = fopen("dump_time.txt","w");
                assert(dump_file);
                fprintf(dump_file,"The dump time was:\n%.17g\n", sim_time);
                fprintf(dump_file,"The dump index was:\n%d\n", plot_number);
                fclose(dump_file);
                fprintf(stdout,"Safety dump successful\n");
            }

            // Die gracefully
            MPI_Finalize(); exit(0);
        }
    }
}
