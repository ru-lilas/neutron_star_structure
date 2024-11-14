#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdarg.h>

/* math constants */
const double pi = 3.141592;

/* physical constants*/
const double c      = 2.9979e10; // speed of light in vacuum (cm s-1)
const double G      = 6.6743e-8; // Gravitational constant (erg cm g-2)

const char *filename_eachNS = "NS_profile_APR4.csv";
/* parameter */
const int ODE_order = 2;
double boundary_rest_mass_density[] = {
    0.0,
    1.512014e+14, // APR4
    5.011872e+14, // fixed 
    1.000000e+15 // fixed
    };
double kappa_coefficient[] = {
    3.5938855153e+13,
    4.6558609352e-08,
    4.2413098391e-17,
    1.2092051829e-15
    };
double adiabatic_index[] = {
    1.3569239500e+00,
    2.8300000000e+00,
    3.4450000000e+00,
    3.3480000000e+00
    };
double a_parameter[] = {
    0.0000000000e+00,
    8.3544443086e-03,
    9.9684658750e-03,
    8.0542854590e-03
    };

double minof(double x, double y) {
    if (x < y) {
        return x;
    }
    else {
        return y;
    }
}

void write_csv(int MODE, const char *filename, const char *format, ...) {
    FILE *csv_file;
    switch (MODE)
    {
        case 1:
            csv_file = fopen(("%s", filename), "w");
            break;
        case 0:
            csv_file = fopen(("%s", filename), "a");
            break;
        default:
            fprintf(stderr, "Error(weite_csv): Invalid mode(%d).\n", MODE);
            exit(EXIT_FAILURE);
    }

    if (csv_file == NULL) {
        fprintf(stderr, "Error: cannot opening csv file.\n");
        exit(EXIT_FAILURE);
    }

    char message[1024];
    va_list args;
    va_start(args, format);
        vfprintf(csv_file, format, args);
    va_end(args);
    fclose(csv_file);
}

double intpow(double base, int exponent) {
    double power = base;
    for (int i=0; i<(exponent-1); i++) {
        power *= base;
    }
    return power;
}

int region(double rest_mass_density) {

    int piecewise_i;

    if 
    ((rest_mass_density < boundary_rest_mass_density[1])) {
        piecewise_i = 0;
    }
    else if
    ((boundary_rest_mass_density[1] <= rest_mass_density) && (rest_mass_density < boundary_rest_mass_density[2])) {
        piecewise_i = 1;
    }
    else if
    ((boundary_rest_mass_density[2] <= rest_mass_density) && (rest_mass_density < boundary_rest_mass_density[3])) {
        piecewise_i = 2;
    }
    else if
    ((boundary_rest_mass_density[3] <= rest_mass_density)) {
        piecewise_i = 3;
    }
    else {
        fprintf(stderr, "Error(region): Can't determine the piecewise region for given rest_mass_density = %e. \n", rest_mass_density);
        exit(EXIT_FAILURE);
    }
    return piecewise_i;
}

double calc_kappa(int piecewise_i) {
/* kappa coefficient */
    double value;
    if ( 0 <= piecewise_i && piecewise_i < sizeof(kappa_coefficient))
    {
        value = kappa_coefficient[piecewise_i];
        return value;
    }
    else
    {
        fprintf(stderr, "Error(calc_kappa): Invalid value piecewise_i = %d\n", piecewise_i);
        exit(EXIT_FAILURE);
    }
}

double calc_gamma(int piecewise_i) {
    double value;
    if ( 0 <= piecewise_i && piecewise_i < sizeof(adiabatic_index))
    {
        value = adiabatic_index[piecewise_i];
        return value;
    }
    else
    {
        fprintf(stderr, "Error(calc_gamma): Invalid value piecewise_i = %d\n", piecewise_i);
        exit(EXIT_FAILURE);
    }
}

double calc_a_parameter(int piecewise_i) {
    double value;
    if ( 0 <= piecewise_i && piecewise_i < sizeof(a_parameter))
    {
        value = a_parameter[piecewise_i];
        return value;
    }
    else
    {
        fprintf(stderr, "Error(calc_a_parameter): Invalid value piecewise_i = %d\n", piecewise_i);
        exit(EXIT_FAILURE);
    }
}

double pressure(double rest_mass_density) {
    int i = region(rest_mass_density);
    double g    = calc_gamma(i);
    double k    = calc_kappa(i);
    double P = k*pow(rest_mass_density, g);
    return P;
}

double calc_energy_density(double rest_mass_density) {
    int piecewise_i = region(rest_mass_density);
    double g    = calc_gamma(piecewise_i);
    double k    = calc_kappa(piecewise_i);
    double a    = calc_a_parameter(piecewise_i);
    double P    = pressure(rest_mass_density);

    double e
    = (1+a)*rest_mass_density*intpow(c,2) + P / (g-1);

    return e;
}

double pressure_gradient(double radius, double rest_mass_density, double enclosed_mass) {
    double r    = radius;
    double rho  = rest_mass_density;
    double m    = enclosed_mass;
    int piecewise_i       = region(rho);
    double g    = calc_gamma(piecewise_i);
    double P    = pressure(rho);
    double e    = calc_energy_density(rho);

    double dPdr =
    - (G*m*e) / (intpow(r,2)*intpow(c,2))
    * (1 + P / e)
    * (1 + (4.0*pi*P*intpow(r,3)) / (m * intpow(c,2)))
    / (1 - (2.0*G*m) / (r * intpow(c,2)) );
}

void thermodynamic_value_writer(void) {
    double rest_mass_density_ini = 5.0e13;
    double rest_mass_density_fin = 3.0e15;
    double rest_mass_density     = rest_mass_density_ini;

    const char *filename = "[FOR_DEBUG]thermodynamic_value_propaty_APR4.csv";
    write_csv(
        1,
        filename,
        "rest_mass_density,kappa,gamma,energy_density,pressure\n"
        );
    for (int i = 1; i < 301; i++) {
        int piecewise_i = region(rest_mass_density);
        double k = calc_kappa(piecewise_i);
        double g = calc_gamma(piecewise_i);
        double e = calc_energy_density(rest_mass_density);
        double P = pressure(rest_mass_density);
        write_csv(
            0,
            filename,
            "%e,%e,%e,%e,%e\n",
            rest_mass_density,k,g,e,P
            );
        rest_mass_density = pow(10, log10(4+i)+13);
    }
}

void nancheck(double value, int err_line, int iteration) {
    if (isnan(value) != 0) {
        fprintf(stderr, "NaN detected. (iteration value %d)\n", iteration);
        return;
    }
}

void derivatives(int ODE_order, double t, double y[], double dydt[]) {
    double r    = t;
    double rho  = y[0];
    double m    = y[1];
    int piecewise_i       = region(rho);
    double g    = calc_gamma(piecewise_i);
    double P    = pressure(rho);
    double e    = calc_energy_density(rho);

    for (int i = 0; i < ODE_order; i++) {
        switch (i)
        {
        case 0:
            dydt[i] = 
            - (G*m*e*rho) / (g*P*intpow(r,2)*intpow(c,2))
            * (1 + P / e )
            * (1 + (4.0*pi*P*intpow(r,3)) / (m*intpow(c,2)) )
            / (1 - (2.0*G*m) / (r*intpow(c,2)) );
            break;

        case 1:
            dydt[i] = 
            4.0 * pi * intpow(r,2) * e / intpow(c,2);
            break;

        default:
            fprintf(stderr,"Error(derivatives): Invalid ODE_order = %d.\n", i);
            exit(EXIT_FAILURE);
            break;
        }
    }
    return;
}

void runge_kutta_csv(int MODE, int ODE_order, double time_value, double y[]) {
    const char *filename = "runge_kutta_output.csv";
    double r    = time_value;
    double rho  = y[0];
    double m    = y[1];
    double P    = pressure(rho);
    double e    = calc_energy_density(rho);

    switch (MODE)
    {
    case 0:
        write_csv(
            0,
            filename,
            "%e,%e,%e,%e,%e\n",
            r,rho,m,P,e
            );
        break;
    case 1:
        write_csv(
            1,
            filename,
            "radius,rest_mass_density,gravitational_mass,pressure,energy_density\n"
        );
        break;
    default:
        break;
    }
}

void runge_kutta_4(int ODE_order, double initial_t, double final_t, int iteration, double y[]) {
    
    /*determine timestep*/
    double timestep = (final_t - initial_t) / iteration;

    // reset step value
    double t = initial_t;
    // initial value of y[]
    double y_initial[ODE_order];
    for (int i = 0; i < ODE_order; i++) {
        y_initial[i] = y[i];
    }
    // crossing velocity
    double crossing_velocity[ODE_order];

    runge_kutta_csv(1,ODE_order,t,y);
    
/*-------------------------------------------------------------------*/
    double *k1      = (double *)malloc(ODE_order * sizeof(double));
    double *k2      = (double *)malloc(ODE_order * sizeof(double));
    double *k3      = (double *)malloc(ODE_order * sizeof(double));
    double *k4      = (double *)malloc(ODE_order * sizeof(double));
    double *ytmp    = (double *)malloc(ODE_order * sizeof(double));
/*-------------------------------------------------------------------*/

//---main loop-------------------------------------------------------------------------//
    for (int i = 0; i < iteration; i++) {

        // y[] in the previous step
        double y_pre[ODE_order];
        for (int j = 0; j < ODE_order; j++) {
            y_pre[j] = y[j];
        }

        // k1
        derivatives(ODE_order, t, y, k1);

        // k2
        for (int j = 0; j < ODE_order; j++) {
            ytmp[j] = y[j] + 0.5 * timestep * k1[j];
        }
        derivatives(ODE_order, t + 0.5 * timestep, ytmp, k2);

        // k3
        for (int j = 0; j < ODE_order; j++) {
            ytmp[j] = y[j] + 0.5 * timestep * k2[j];
        }
        derivatives(ODE_order, t + 0.5 * timestep, ytmp, k3);

        // k4
        for (int j = 0; j < ODE_order; j++) {
            ytmp[j] = y[j] + timestep * k3[j];
        }
        derivatives(ODE_order, t + timestep, ytmp, k4);

        /* update y */
        for (int j = 0; j < ODE_order; j++) {
            crossing_velocity[j] = (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]) / 6.0;
            double add_y =  crossing_velocity[j] * timestep;
            y[j] += add_y;
        }

    runge_kutta_csv(0,ODE_order,t,y);

    for (int j = 0; j < ODE_order; j++) {
        if (y[j] < 1.0){
            write_csv(0,filename_eachNS,"%e,%e,%e\n", y_initial[0],t,y[1]);
            return;
        }
    }

    timestep 
        = minof(
        -(y[0] / crossing_velocity[0])* 0.01,
        (y[1] / crossing_velocity[1])*0.01
        );
        
    t += timestep;
    }
    return;
}

/*=============================================================================*/
int main(void) {
    // thermodynamic_value_writer(); //FOR DEBUG
    double y[ODE_order];

    /* stepvalue */
    double  r_ini    = 1.0e2;
    double  r_fin    = 3.0e6;
    int     max_iteration = 10000;

    // density at the centre of NS
    double  central_density_ini = 1.0e14;
    double  central_density_fin = 2.6e15;
    // double  central_density = 2.0e15;
    double  central_density = central_density_ini;
    int N = 100;

    write_csv(1, filename_eachNS,"central_density,radius,mass\n");

    for (int i = 0; i < N; i++) {
        central_density = (1.0 + i*0.1)*intpow(10,14);
        // initial condition
        y[0] = central_density; // g cm-3
        y[1] = (4 * pi * intpow(r_ini, 3) * y[0]) / 3.0; // g

        runge_kutta_4(ODE_order, r_ini, r_fin, max_iteration, y);

    }
    for (int i = 0; i < N*3; i++) {
        central_density = (1.0 + i*0.01)*intpow(10,15);

        // initial condition
        y[0] = central_density; // g cm-3
        y[1] = (4 * pi * intpow(r_ini, 3) * y[0]) / 3.0; // g

        runge_kutta_4(ODE_order, r_ini, r_fin, max_iteration, y);

    }

return 0;
}