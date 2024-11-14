#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdarg.h>

// log level
#define LOG_DEBUG 0
#define LOG_INFO  1
#define LOG_WARN  2
#define LOG_ERROR 3
#define WRITING_DATA 0
#define WRITING_HEADER 1

/* math constants */
const double pi = 3.141592;

/* physical constants*/
const double c      = 2.9979e10; // speed of light in vacuum (cm s-1)
const double G      = 6.6743e-8; // Gravitational constant (erg cm g-2)
const double Msol   = 1.998e33; // solar mass(g)

/* parameter */
const int eq_order = 2;

//===Prototype Declaration=============================================================//
/* Firm subroutine (DO NOT CHANGE) */
double  intpow(double base, int exponent); // power for int value
void nancheck(double value, int err_line, int iteration); // NaN checker
void write_csv(int MODE, const char *filename, const char *format, ...); // output to csv
void write_log(int level, int err_line, const char *format, ...); // writing log file
double gamma_index(int piecewise_index) {
/* gamma index in EoS (DO NOT REWRITE VALUES) */
    double value;
    switch (piecewise_index)
    {
    case 0:
        value = 1.356924e+00;
        break;
    case 1:
        value = 2.830000e+00;
        break;
    case 2:
        value = 3.445000e+00;
        break;
    case 3:
        value = 3.348000e+00;
        break;
    default:
        return -1;
        break;
    }
    return value;
}
//-------------------------------------------------------------------------------------//

/* changeable */
double pressure(double rho); // Equation of State
int region(double rho); // determine value of piecewise subscript
void derivatives(int eq_order, double t, double y[], double dydt[]); // righthand side of ODEs
void runge_kutta_4(int eq_order, double t0, int N, double t_fin, double y[]); // rungekutta method (4th orders)
void runge_kutta_4_add_write(int MODE, double t, double y[]);
double kappa(int piecewise_index) {
/* kappa coefficient */
    double value;
    switch (piecewise_index)
    {
    case 0:
        value = 3.5948855153e+13;  // Hotokezaka et al.2013 eq.(3)
        break;
    case 1:
        value = 4.6558609352e-08;
        break;
    case 2:
        value = 4.2413098391e-17;
        break;
    case 3:
        value = 1.2092051829e-15;
        break;
    default:
        return -1;
        break;
    }
    return value;
}

//===MAIN========================================================================================//
int main(void)
{
    write_log(LOG_INFO, __LINE__, "Start to calculate.\n");

    /* dependent variable vector */
    double y[eq_order]; // y[0]:density y[1]:mass

    /* parameter setting to determine step */
    double  m_initial       = 1.0e22; // g
    double  m_final         = 3.0e33; // g
    int     max_iteration   = 10000;

    /* density at the centre of NS */
    double rho0_initial = 1.1e15; // g.cm-3

    /* output file reset */
    write_csv(-1, "runge_kutta_4_output.csv", "");

    // TODO: Insert a FOR sentence here when solving for various rho0
    double rho0 = rho0_initial;

    /* initial condition  */
    y[0] = rho0; // g cm-3
    y[1] = pow(((3.0*m_initial) / (4.0*pi*rho0)), 1.0/3.0); // cm

    /* runge kutta method */
    runge_kutta_4(eq_order, m_initial, max_iteration, m_final, y);

    write_log(LOG_INFO, __LINE__, "Finished calculating.\n");

return 0;
}
//===End of main====================================================================//

double pressure(double rho) {
    int i = region(rho);
    double P = kappa(i)*pow(rho, gamma_index(i));
    return P;
}

int region(double rho) {
    double boundary_rho[4] = { // g.cm-3
    0.0, // fixed
    1.512014e+14, // free parameter(APR4)
    5.011872e+14, // fixed 
    1.000000e+15, // fixed
    };

    int i;

    if 
    ((rho < boundary_rho[1])) {
        i = 0;
    }
    else if
    ((boundary_rho[1] <= rho) && (rho < boundary_rho[2])) {
        i = 1;
    }
    else if
    ((boundary_rho[2] <= rho) && (rho < boundary_rho[3])) {
        i = 2;
    }
    else if
    ((boundary_rho[3] <= rho)) {
        i = 3;
    }
    else {
        write_log(LOG_ERROR, __LINE__, "can't determine region for given rho = %e. \n", rho);
        exit(EXIT_FAILURE);
    }

    return i;
}

void derivatives(int eq_order, double t, double y[], double dydt[]){
    double m    = t; // g
    double rho  = y[0]; // g.cm-3
    double r    = y[1]; // cm
    int j = region(rho); // piecewise index
    double P    = pressure(rho);
    double TOV_part[4]; // product terms in TOV eq.
        // TOV_part[0] = -(G*m*pow(rho, 1.0-gamma_index(j)))/(kappa(j)*gamma_index(j)*intpow(r,2));
        TOV_part[0] = -(G*m*pow(rho, 2.0-gamma_index(j))) / (gamma_index(j) * kappa(j) * intpow(r,2));
        TOV_part[1] = (1.0 + (kappa(j) * pow(rho, gamma_index(j)-1)) / intpow(c,2));
        TOV_part[2] = (1.0 + (4 * pi * intpow(r,3) * kappa(j) * pow(rho, gamma_index(j))) / (m * intpow(c,2)));
        TOV_part[3] = (1.0 - ( (2*G*m) / (r * intpow(c,2)) ));
    // linear scale TOV eq.
    double drho_dr =
        TOV_part[0]*TOV_part[1]*TOV_part[2]/TOV_part[3];
    // mass conservation in linear scale
    double drdm =
        1.0/(4.0*pi*intpow(r,2)*rho);

    for (int i = 0; i < eq_order; i++){
        switch (i)
        {
        case 0: // TOV equation
            dydt[i] = rho/(gamma_index(j)*P) * drho_dr;
            break;

        case 1: // mass conservation
            dydt[i] = drdm;
            break;

        default:
            dydt[i] = 0;

            break;
        };
    };
return;
}

void runge_kutta_4(int eq_order, double t0, int N, double t_fin, double y[]) {

/* determine step */
// FIXME: change adaptive step 
double step = (t_fin - t0) / N;

double t = t0; // reset step value

double y_ini[eq_order]; // initial conditions

/* memorize initial conditions */
for (int i = 0; i < eq_order; i++) {
    y_ini[i] = y[i];
}

/*-------------------------------------------------------------------*/
double *k1      = (double *)malloc(eq_order * sizeof(double));
double *k2      = (double *)malloc(eq_order * sizeof(double));
double *k3      = (double *)malloc(eq_order * sizeof(double));
double *k4      = (double *)malloc(eq_order * sizeof(double));
double *ytmp    = (double *)malloc(eq_order * sizeof(double));
/*-------------------------------------------------------------------*/

//---writing csv-----------------------------------------------------------------------//
    char *filename = "runge_kutta_4_output.csv";
    write_csv(0, filename, "t");
    for (int i = 0; i< eq_order; i++) {
        write_csv(0, filename, ",y%d,k1[%d],k2[%d],k3[%d],k4[%d]", i,i,i,i,i);
    }
    runge_kutta_4_add_write(1,t,y);
    write_csv(0, filename, "\n");
//-------------------------------------------------------------------------------------//

//---main loop-------------------------------------------------------------------------//
for (int i = 0; i < N; i++) {

    double y_pre[eq_order]; // y in the previous step
    /* writing y_pre */
    for (int j = 0; j < eq_order; j++) {
        y_pre[j] = y[j];
    }

    // k1
    derivatives(eq_order, t, y, k1);

    // k2
    for (int j = 0; j < eq_order; j++) {
        ytmp[j] = y[j] + 0.5 * step * k1[j];
        nancheck(ytmp[j], __LINE__, i);
    }
    derivatives(eq_order, t + 0.5 * step, ytmp, k2);

    // k3
    for (int j = 0; j < eq_order; j++) {
        ytmp[j] = y[j] + 0.5 * step * k2[j];
        nancheck(ytmp[j], __LINE__, i);
    }
    derivatives(eq_order, t + 0.5 * step, ytmp, k3);

    // k4
    for (int j = 0; j < eq_order; j++) {
        ytmp[j] = y[j] + step * k3[j];
        nancheck(ytmp[j], __LINE__, i);
    }
    derivatives(eq_order, t + step, ytmp, k4);

    /*---writing csv-------------------------------------------------------------------*/
    write_csv(0, filename, "%e", t);
    for (int j = 0; j< eq_order; j++) {
        write_csv(0, filename, ",%e,%e,%e,%e,%e", y[j], k1[j], k2[j], k3[j], k4[j]);
    }
    runge_kutta_4_add_write(0,t,y);
    write_csv(0, filename, "\n");
    //---------------------------------------------------------------------------------//

    /* update y */
    for (int j = 0; j < eq_order; j++) {
        double add_y =  (step / 6.0) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
        y[j] += add_y;
    }

    if ((log10(pressure(pow(10,y[0]))))<1.0e-1) {
    // TODO: writing csv total mass and radius
    return;
    }

    /* next step */
    t += step;
}

    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(ytmp);

    write_log(LOG_WARN, __LINE__, "We reached max step N=%d\n",N);

return;
}

/*---runge_kutta subrutine-----------------------------------------------------------*/
// t,y以外に追加で何かしらのデータを書き込みたいときに使う
void runge_kutta_4_add_write(int MODE, double t, double y[]) {
    double r = y[1];
    double rho = y[0];
    double m = t;
    double P = pressure(rho);
    double dPdr = 
        -((G * m * rho) / intpow(r,2))
        *(1 + P / (rho * intpow(c,2)))
        *(1 + (4.0 * pi * intpow(r,3) * P) / (m * intpow(c,2)))
        /(1 - ( (2*G*m) / (r * intpow(c,2)) ))
        ;
    double H = -P / dPdr; // scale height

    switch (MODE)
    {
    case 0: // writing data
        write_csv(0, "runge_kutta_4_output.csv", ",%e,%e", P,H);
        break;
    case 1: // writing header
        write_csv(0, "runge_kutta_4_output.csv", ",P,H");
        break;
    default:
        write_log(LOG_ERROR, __LINE__, "Unknown writing MODE %d.\n",MODE);
        perror("Error was detected. See '.log' file\n");
        exit(EXIT_FAILURE);
        break;
    }
}
//===Firm subroutines (DO NOT CHANGE)============================================================//
void write_csv(int MODE, const char *filename, const char *format, ...) {
    FILE *csv_file;
    switch (MODE)
    {
        case -1:
            csv_file = fopen(("%s", filename), "w");
            break;
        case 0:
            csv_file = fopen(("%s", filename), "a");
    }

    if (csv_file == NULL) {
        perror("Error: cannot opening csv file.\n");
        exit(EXIT_FAILURE);
    }

    char message[1024];
    va_list args;
    va_start(args, format);
        vfprintf(csv_file, format, args);
    va_end(args);
    fclose(csv_file);
}

void write_log(int level, int err_line, const char *format, ...) {
    FILE *log_file = fopen(".log", "a");
    if (log_file == NULL) {
        printf("Error opening log file.\n");
        return;
    }

    // get current time
    time_t now;
    time(&now);
    struct tm *local = localtime(&now);

    // convert value of LOG_LEVEL into str
    const char *level_str;
    switch (level) {
        case LOG_DEBUG: level_str = "DEBUG"; break;
        case LOG_INFO:  level_str = "INFO"; break;
        case LOG_WARN:  level_str = "WARN"; break;
        case LOG_ERROR: level_str = "ERROR"; break;
        default: level_str = "UNKNOWN"; break;
    }

    // writing log file
    char message[100];
    fprintf(log_file, "%04d-%02d-%02d %02d:%02d:%02d [%s] [line %d] ",
            local->tm_year + 1900, local->tm_mon + 1, local->tm_mday,
            local->tm_hour, local->tm_min, local->tm_sec,
            level_str, err_line);

    va_list args;
    va_start(args, format);
    vfprintf(log_file, format, args);
    va_end(args);
    
    fclose(log_file);
}

double intpow(double base, int exponent) {
    double power = base;
    for (int i=0; i<(exponent-1); i++) {
        power *= base;
    }
    return power;
}

void nancheck(double value, int err_line, int iteration) {
    if (isnan(value) != 0) {
        write_log(LOG_WARN, err_line,
        "NaN detected. (iteration value %d)\n", iteration
        );
        exit(EXIT_FAILURE);
    }
}
