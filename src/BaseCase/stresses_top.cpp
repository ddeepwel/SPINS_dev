#include "../BaseCase.hpp"

// Top stresses
void BaseCase::stresses_top(TArrayn::DTArray & u, TArrayn::DTArray & v, TArrayn::DTArray & w,
        TArrayn::DTArray & Hprime, TArrayn::DTArray & temp, TArrayn::Grad * gradient_op,
        const string * grid_type, const double mu, double time, int itercount, bool restarting) {
    // set-up
    static DTArray *tx = alloc_array(size_x(),size_y(),1);
    static DTArray *ty = alloc_array(size_x(),size_y(),1);
    blitz::firstIndex ii;
    blitz::secondIndex jj;

    // top stress ( along channel - x )
    top_stress_x(*tx, u, temp, gradient_op, grid_type, size_z(), mu);
    // integral values (the first is the surface force)
    double fx_tot = pssum(sum(
                (*get_quad_x())(ii)*
                (*get_quad_y())(jj)*(*tx)));
    double fx_abs = pssum(sum(
                (*get_quad_x())(ii)*
                (*get_quad_y())(jj)*abs(*tx)));
    // extrema values
    double tx_max = psmax(max(*tx));
    double tx_min = psmin(min(*tx));

    // top stress ( across channel - y )
    top_stress_y(*ty, v, temp, gradient_op, grid_type, size_z(), mu);
    // integral values (the first is the surface force)
    double fy_tot = pssum(sum(
                (*get_quad_x())(ii)*
                (*get_quad_y())(jj)*(*ty)));
    double fy_abs = pssum(sum(
                (*get_quad_x())(ii)*
                (*get_quad_y())(jj)*abs(*ty)));
    // extrema values
    double ty_max = psmax(max(*ty));
    double ty_min = psmin(min(*ty));
    double ts_max = psmax(max(pow(pow(*tx,2)+pow(*ty,2),0.5)));

    // write to a stress diagnostic file
    if (master()) {
        FILE * stresses_file = fopen("stresses_top.txt","a");
        assert(stresses_file);
        if ( itercount==0 and !restarting )
            fprintf(stresses_file,"Time, "
                    "tx_max, tx_min, ty_max, ty_min, ts_max, "
                    "fx_tot, fx_abs, fy_tot, fy_abs\n");
        fprintf(stresses_file,"%.17f, "
                "%.17g, %.17g, %.17g, %.17g, %.17g, "
                "%.17g, %.17g, %.17g, %.17g\n",
                time,
                tx_max, tx_min, ty_max, ty_min, ts_max,
                fx_tot, fx_abs, fy_tot, fy_abs);
        fclose(stresses_file);
    }
}
