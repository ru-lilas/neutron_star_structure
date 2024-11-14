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
const double c      = 2.9979e10; // speed of light in vacuum (cm/s)
const double h      = 6.6260e-27; // Planck constant (erg s)
const double G      = 6.6743e-8; // Gravitational constant (erg cm g-2)
const double m_e    = 9.1093e-28; // electron mass (g)
const double m_p    = 1.672e-24; // proton mass(g)
const double Msol   = 1.998e33; // solar mass(g)
const double Rsol   = 6.957e10; // solar radius(cm)

/* parameter */
const int eq_order = 2;
// const double mu = 2.0;

//===Prototype Declaration=============================================================//
double  intpow(double base, int exponent); // power for int value
// double pressure(double rho);
int region(double rho); // determine value of piecewise subscript
void write_csv(int MODE, const char *filename, const char *format, ...); // output to csv
void write_log(int level, int err_line, const char *format, ...);
void nancheck(double value, int err_line, int iteration);
double pressure(double rho);
void derivatives(int eq_order, double t, double y[], double dydt[]);
void runge_kutta_4(int eq_order, double h0, double timestep, double step_final, double y[]);
void runge_kutta_4_add_write(int MODE, double t, double y[]);
//---gamma(fixed)----------------------------------------------------------------------//
double gamma_index(int piecewise_index) {
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
//---kappa-----------------------------------------------------------------------------//
double kappa(int piecewise_index) {
    double value;
    switch (piecewise_index)
    {
    case 0:
        value = 3.5938855153e+13;
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
    write_log(LOG_INFO, __LINE__, "Start calculation.\n");
    double r_initial = 1.0e0; // cm
    double r_interval = 1.0e3; // cm
    double r_final = 1.0e7; // cm

    // double rho0_initial = 1.0e13; // g cm-3
    double rho0_initial = 1.0e15; // g cm-3
    // double rho0_interval= 2.0e13; // g cm-3
    double rho0_interval= 1.0e14; // g cm-3
    // double rho0_final = 3.0e15; // g cm-3
    double rho0_final = 2.0e15; // g cm-3

    int rho0_step = (int)((rho0_final - rho0_initial)/rho0_interval);

    write_csv(-1, "MRrelation_APR4.csv", "rho,r,mass\n");
    write_csv(-1, "runge_kutta_4_output.csv", "");

    double rho0 = rho0_initial;

    for ( int j = 0; j < rho0_step; j++) {
        double y[eq_order]; // y[0]:density y[1]:mass

        // initialize y[0], y[1] for given rho0
        y[0] = rho0; // g cm-3
        y[1] = (4 * pi * intpow(r_initial, 3) * y[0]) / 3.0; // g

        // determine value of piecewise subscript
        int section = region(y[0]);

        runge_kutta_4(eq_order, r_initial, r_interval, r_final, y);

        rho0 += rho0_interval;
    }

return 0;
}
//===End of main====================================================================//

//===equation of state==============================================================//
double pressure(double rho) {
    int i = region(rho);
    double P = kappa(i)*pow(rho, gamma_index(i));
    return P;
}
//---region-------------------------------------------------------------------------//
int region(double rho) {
    double boundary_rho[4] = {
    0.0, // fixed
    1.512014e+14, // free parameter
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
        printf("Error in function 'region': can't determine region from given rho = %e.\n", rho);
        i = 0; // debug
        exit(EXIT_FAILURE);
    }

    return i;
}

//===derivatives====================================================================//
// FIXME: massが大きすぎる
void derivatives(int eq_order, double t, double y[], double dydt[]){
    double rho = y[0];
    double m = y[1];
    double P = pressure(rho);
    double r = t;
    int j = region(rho);

    for (int i = 0; i < eq_order; i++){
        switch (i)
        {
        case 0: // TOV equation
            dydt[i] =
            -((G * m * pow(rho, 2.0-gamma_index(j))) / ((kappa(j))*gamma_index(j)*intpow(r,2)))
            *(1 + (kappa(j) * pow(rho, gamma_index(j)-1)) / intpow(c,2))
            *(1 + (4 * pi * intpow(r,3) * kappa(j) * pow(rho, gamma_index(j))) / (m * intpow(c,2)))
            /(1 - ( (2*G*m) / (r * intpow(c,2)) ))
            ;
            break;
        case 1: // mass conservation
            dydt[i] = 4*pi*intpow(r,2)*rho;
            break;
        default:
            dydt[i] = 0;
            break;
        };
    };
return;
}

//===Runge Kutta method(4 orders)================================================================//
void runge_kutta_4(int eq_order, double h0, double timestep, double step_final, double y[])
{
int N = (int)( (step_final-h0) / timestep ); // total number of step
double t = h0; // step reset
double y0 = y[0];

// writing header
char *filename = "runge_kutta_4_output.csv";
write_csv(0, filename, "t");
for (int i = 0; i< eq_order; i++) {
    write_csv(0, filename, ",y%d", i);
}
runge_kutta_4_add_write(WRITING_HEADER,t,y);
write_csv(0, filename, "\n");

/*-------------------------------------------------------------------*/
double *k1      = (double *)malloc(eq_order * sizeof(double));
double *k2      = (double *)malloc(eq_order * sizeof(double));
double *k3      = (double *)malloc(eq_order * sizeof(double));
double *k4      = (double *)malloc(eq_order * sizeof(double));
double *ytmp    = (double *)malloc(eq_order * sizeof(double));
/*-------------------------------------------------------------------*/

for (int i = 0; i < N; i++) {

write_csv(0, filename, "%e", t);

if ((i % 1000)==0){
// printf("Runge Kutta i =%d\n", i);
}

for (int j = 0; j< eq_order; j++) {
    write_csv(0, filename, ",%e", y[j]);
}
runge_kutta_4_add_write(0,t,y);
write_csv(0, filename, "\n");

    // k1
    derivatives(eq_order, t, y, k1);

    // k2
    for (int j = 0; j < eq_order; j++) {
        ytmp[j] = y[j] + 0.5 * timestep * k1[j];
        // nancheck(ytmp[j], __LINE__, i);
    }
    // FIXME: ytmpが負の値をとるとNaNが返されてしまう
    derivatives(eq_order, t + 0.5 * timestep, ytmp, k2);

    // k3
    for (int j = 0; j < eq_order; j++) {
        ytmp[j] = y[j] + 0.5 * timestep * k2[j];
        // nancheck(ytmp[j], __LINE__, i);
    }
    derivatives(eq_order, t + 0.5 * timestep, ytmp, k3);

    // k4
    for (int j = 0; j < eq_order; j++) {
        ytmp[j] = y[j] + timestep * k3[j];
        // nancheck(ytmp[j], __LINE__, i);
    }
    derivatives(eq_order, t + timestep, ytmp, k4);

    // 古いt, y[0]を保存
    double pre_y0 = y[0];
    double pre_y1 = y[1];

    // yの更新
    for (int j = 0; j < eq_order; j++) {
        double add_y =  (timestep / 6.0) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
        y[j] += add_y;
    }

    if ( (isnan(y[0])) || (isnan(y[1])) || (pre_y0*y[0] < 0) || (pre_y1*y[1] < 0)) {
    // if ((pre_y0*y[0] < 0) || (pre_y1*y[1] < 0)) {
        
        double real_final_step = t+0.5*timestep;
        // ロールバック
        y[0] = pre_y0;
        y[1] = pre_y1;

    write_csv(0, "MRrelation_APR4.csv", "%e,%e,%e\n", y0, real_final_step, y[1]); // total mass
    write_csv(0, filename, "\n");

        return;
    }

    /* step variableを進める */
    t += timestep;
}

free(k1);
free(k2);
free(k3);
free(k4);
free(ytmp);

printf("We reached max N.\n");

return;
}
/*---runge_kutta subrutine-----------------------------------------------------------*/
// t,y以外に追加で何かしらのデータを書き込みたいときに使う
void runge_kutta_4_add_write(int MODE, double t, double y[]) {
    double r = t;
    double rho = y[0];
    double m = y[1];
    double P = pressure(y[0]);
    double dPdr = 
        -((G * m * rho) / intpow(r,2))
        *(1 + P / (rho * intpow(c,2)))
        *(1 + (4 * pi * intpow(r,3) * P) / (m * intpow(c,2)))
        /(1 - ( (2*G*m) / (r * intpow(c,2)) ))
        ;
    double H = -P / dPdr; // scale height

    switch (MODE)
    {
    case 0: // writing data
        write_csv(0, "runge_kutta_4_output.csv", ",%e,%e,%e", P,-dPdr,H);
        break;
    case 1: // writing header
        write_csv(0, "runge_kutta_4_output.csv", ",P,-dPdr,H");
        break;
    default:
        write_log(LOG_ERROR, __LINE__, "Unknown writing MODE %d.\n",MODE);
        perror("Error was detected. See '.log' file\n");
        exit(EXIT_FAILURE);
        break;
    }
}
//===変更不要の部品たち===============================================================//
/* csvへの出力 */
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
/* 整数冪を返す関数 */
double intpow(double base, int exponent) {
    double power = base;
    for (int i=0; i<(exponent-1); i++) {
        power *= base;
    }
    return power;
}
/* NaN check */
void nancheck(double value, int err_line, int iteration) {
    if (isnan(value) != 0) {
        write_log(LOG_WARN, err_line,
        "NaN detected. (iteration value %d)\n", iteration
        );
        exit(EXIT_FAILURE);
    }
}

/* log */
void write_log(int level, int err_line, const char *format, ...) {
    FILE *log_file = fopen(".log", "a");
    if (log_file == NULL) {
        printf("Error opening log file.\n");
        return;
    }

    // 現在の時刻を取得
    time_t now;
    time(&now);
    struct tm *local = localtime(&now);

    // ログレベルを文字列に変換
    const char *level_str;
    switch (level) {
        case LOG_DEBUG: level_str = "DEBUG"; break;
        case LOG_INFO:  level_str = "INFO"; break;
        case LOG_WARN:  level_str = "WARN"; break;
        case LOG_ERROR: level_str = "ERROR"; break;
        default: level_str = "UNKNOWN"; break;
    }

    // ログファイルに書き込み
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