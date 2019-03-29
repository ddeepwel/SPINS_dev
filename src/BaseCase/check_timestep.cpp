#include "../BaseCase.hpp"

double BaseCase::check_timestep(double t_step, double now) {
    // Check time step
    if (t_step < 1e-9) {
        // Timestep's too small, somehow stuff is blowing up
        if (master()) fprintf(stderr,"Tiny timestep (%e), aborting\n",t_step);
        return -1;
    } else if (t_step > get_dt_max()) {
        // Cap the maximum timestep size
        t_step = get_dt_max();
    }

    // Now, calculate how many timesteps remain until the next writeout
    double until_plot = get_next_plot() - now;
    int steps = ceil(until_plot / t_step);
    // Where will we be after (steps) timesteps of the current size?
    double real_until_plot = steps*t_step;
    // If that's close enough to the real writeout time, that's fine.
    if (fabs(until_plot - real_until_plot) < 1e-6) {
        return t_step;
    } else {
        // Otherwise, square up the timeteps.  This will always shrink the timestep.
    return (until_plot / steps);
    }
}
