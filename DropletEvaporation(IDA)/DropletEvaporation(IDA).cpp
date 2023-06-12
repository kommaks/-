
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <cmath>
using namespace std;

#include <ida/ida.h>                          /* prototypes for IDA fcts., consts.    */
#include <nvector/nvector_serial.h>           /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h>        /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h>        /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>          /* defs. of realtype, sunindextype      */
#include <sunnonlinsol/sunnonlinsol_newton.h> /* access to Newton SUNNonlinearSolver  */



 /* Problem Constants */

#define FTOL   RCONST(1.e-12) /* function tolerance */
#define STOL   RCONST(1.e-12) /* step tolerance     */

#define ZERO   RCONST(0.0)
#define PT25   RCONST(0.25)
#define PT5    RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)
#define TWO    RCONST(2.0)

#define PI     RCONST(3.1415926)
#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/* Macro to define dense matrix elements, indexed from 1. */

#define Ith(v,i)    NV_Ith_S(v,i-1)
#define IJth(A,i,j) SM_ELEMENT_D(A,i-1 ,j-1)

/* Prototypes of functions called by IDA */

int resrob(realtype tres, N_Vector yy, N_Vector yp,
    N_Vector resval, void* user_data);

/* Prototypes of private functions */
static void PrintHeader(realtype rtol, N_Vector avtol, N_Vector y);
static void PrintOutput(void* mem, realtype t, N_Vector y);
static void PrintRootInfo(int root_f1, int root_f2);
static int check_retval(void* returnvalue, const char* funcname, int opt);
static int check_ans(N_Vector y, realtype t, realtype rtol, N_Vector atol);

typedef struct UserDataStr{
    realtype* r;
    realtype* T;
    int N;
    int NEQ;
    int s;
    int kp;
    realtype Tr;
    realtype M;
    //My variables
    realtype* nevyaz;
    realtype* ACp;
    realtype* ADensity;
    realtype* ALambda;
    realtype a;
    realtype inter;
    realtype Ti;
    realtype p;
    realtype L_d;
    ofstream *outNevyaz;
    void* mykmem;
} UserData;

//Functions for IDA
void ExportToArray(double* T_vect, double& M, UserData *data, N_Vector yy, int N)
{
    for (int i = 0; i <= N + 1; i++)
    {
        T_vect[i] = Ith(yy, i);
        //cout << "T_vect  " << i << " =  " << Ith(yy, j - 1 + 1) << endl;
    }
    //cout << "data->Tr; = " << data->Tr << "\n";
    //T_vect[N] = data->Tr;
    //cout << "MyNx = " << myNx << "\n";
    //cout << "MyNm = " << myNm << "\n";
}
//------------------------------------------
// Functions from DropletWarming
//------------------------------------------

double Lagrange1(const vector<double>& x, const vector<double>& y, int n, double _x)
{
    double result = 0.0;
    for (int i = 0; i < n; i++)
    {
        double P = 1.0;
        for (int j = 0; j < n; j++)
        {
            if (j != i)
                P *= (_x - x[j]) / (x[i] - x[j]);
        }
        result += P * y[i];
    }
    return result;
}
//Linear interpolation
double LinIntWConductivity(double Tmprtr)
{
    /*
    const vector<double> LAGR_T = { 273.16, 278.16, 283.16, 288.16, 293.16, 298.16, 303.16, 308.16, 313.16, 318.16, 323.16, 328.16, 333.16,
        338.16, 343.16, 348.16, 353.16, 358.16, 363.16, 368.16, 373.16, 378.16, 383.16, 388.16, 393.16, 398.16, 403.16, 408.16, 413.16, 418.16,
        423.16, 428.16, 433.16, 438.16, 443.16, 448.16, 453.16, 458.16, 463.16, 468.16, 473.16, 478.16, 483.16, 488.16, 493.16, 498.16, 503.16,
        508.16, 513.16, 518.16, 523.16, 528.16, 533.16, 538.16, 543.16, 548.16, 553.16, 558.16, 563.16, 568.16, 573.16, 578.16, 583.16, 588.16,
        593.16, 598.16, 603.16, 608.16, 613.16, 618.16, 623.16, 628.16, 633.16, 638.16, 643.16 };
    const vector<double> LAGR_COND = { 0.56104, 0.57054, 0.58002, 0.58935, 0.59843, 0.60717, 0.61547, 0.6233, 0.6306, 0.63736, 0.64356, 0.64923,
        0.65436, 0.65897, 0.6631, 0.66676, 0.66999, 0.67282, 0.67526, 0.67734, 0.6791, 0.68054, 0.68169, 0.68257, 0.68319, 0.68356, 0.6837,
        0.68361, 0.6833, 0.68277, 0.68204, 0.6811, 0.67995, 0.6786, 0.67705, 0.67529, 0.67332, 0.67114, 0.66875, 0.66614, 0.66331, 0.66025,
        0.65696, 0.65343, 0.64964, 0.6456, 0.6413, 0.63671, 0.63184, 0.62667, 0.62118, 0.61537, 0.60923, 0.60274, 0.5959, 0.5887, 0.58113,
        0.57321, 0.56494, 0.55633, 0.54741, 0.53819, 0.52873, 0.51904, 0.50918, 0.49917, 0.48905, 0.47883, 0.46849, 0.458, 0.44735, 0.43652,
        0.4257, 0.41636, 0.42513 };
    */
    const vector<double> LAGR_T = { 273.16, 293.16, 313.16, 333.16, 353.16, 373.16, 393.16, 413.16, 433.16, 453.16,473.16, 493.16, 513.16,
        533.16, 553.16, 573.16, 593.16, 613.16, 633.16 };
    const vector<double> LAGR_COND = { 0.56104, 0.57054, 0.58002, 0.58935, 0.59843, 0.60717, 0.61547, 0.6233, 0.6306, 0.63736, 0.64356, 0.64923,
        0.65436, 0.65897, 0.6631, 0.66676, 0.66999, 0.67282, 0.67526, 0.67734, 0.6791, 0.68054, 0.68169, 0.68257, 0.68319, 0.68356, 0.6837,
        0.68361, 0.6833, 0.68277, 0.68204, 0.6811, 0.67995, 0.6786, 0.67705, 0.67529, 0.67332, 0.67114, 0.66875, 0.66614, 0.66331, 0.66025,
        0.65696, 0.65343, 0.64964, 0.6456, 0.6413, 0.63671, 0.63184, 0.62667, 0.62118, 0.61537, 0.60923, 0.60274, 0.5959, 0.5887, 0.58113,
        0.57321, 0.56494, 0.55633, 0.54741, 0.53819, 0.52873, 0.51904, 0.50918, 0.49917, 0.48905, 0.47883, 0.46849, 0.458, 0.44735, 0.43652,
        0.4257, 0.41636, 0.42513 };
    int nx = LAGR_T.size();
    int ny = LAGR_COND.size();
    vector <double> LAGR_COND_copy;
    int shag = 4;
    LAGR_COND_copy.resize(ny / shag + 1);
    int j = 0;
    for (int i = 0; i < ny; i += shag)
    {
        LAGR_COND_copy[j] = LAGR_COND[i];
        j += 1;
    }
    return Lagrange1(LAGR_T, LAGR_COND_copy, nx, Tmprtr);
}
//Formula setting for Heat Capacity and it's derivative, J/(kg*K)
double Cp(double Tmprtr, const double a)
{
    const double Cp_water = 4180.;
    const double Cp_steam = 2000.;
    //Molar mass of water
    const double W = 18 * pow(10, -3);
    const double T_boiling = 373;
    double AA, B, C, D, E;
    //Polinom for Tmprtr = (K) /1000
    double T_forCp = Tmprtr / 1000.;
    if (Tmprtr <= T_boiling)
    {
        AA = -203.6060;
        B = 1523.290;
        C = -3196.413;
        D = 2474.455;
        E = 3.855326;
    }
    else
    {
        AA = 30.09200;
        B = 6.832514;
        C = 6.793435;
        D = -2.534480;
        E = 0.082139;
    }
    //Get value
    if (a < pow(10, -8))
    {
        if (Tmprtr <= T_boiling)
            return Cp_water;
        else
            return Cp_steam;
    }
    else
        return (AA + B * T_forCp + C * pow(T_forCp, 2.) + D * pow(T_forCp, 3.) + E * pow(T_forCp, -2.)) / W;
}
//J/(kg*K^2)
double DfCp(double Tmprtr, const double a)
{
    //Molar mass of water
    const double W = 18 * pow(10, -3);
    const double T_boiling = 373;
    double AA, B, C, D, E;
    //Polinom for Tmprtr = (K) /1000
    double T_forCp = Tmprtr / 1000.;
    if (Tmprtr <= T_boiling)
    {
        AA = -203.6060;
        B = 1523.290;
        C = -3196.413;
        D = 2474.455;
        E = 3.855326;
    }
    else
    {
        AA = 30.09200;
        B = 6.832514;
        C = 6.793435;
        D = -2.534480;
        E = 0.082139;
    }
    //Get value
    if (a < pow(10, -8))
        return 0;
    else
        return (B + 2. * C * T_forCp + 3. * D * pow(T_forCp, 2.) - 2. * E * pow(T_forCp, -3.)) / 1000. / W;
}
//Formula setting for Density and it's derivative, kg/m^3
double Density(double Tmprtr, const double a)
{
    const double rho_water = 980.;
    const double rho_steam = 0.4;
    const double T_boiling = 373;
    const double R = 8.314462;
    const double p = pow(10, 5);
    const double W = 18 * pow(10, -3);
    if (a < pow(10, -8))
    {
        if (Tmprtr <= T_boiling)
            return rho_water;
        else
            return rho_steam;
    }
    else
        if (Tmprtr <= T_boiling)
            return rho_water;
        else
            return p * W / (R * Tmprtr);
}
//kg/(m^3*K)
double DfDensity(double Tmprtr, const double a)
{
    const double T_boiling = 373;
    const double R = 8.314462;
    const double p = pow(10, 5);
    const double W = 18 * pow(10, -3);
    if (a < pow(10, -8))
        return 0;
    else
    {
        if (Tmprtr <= T_boiling)
            return 0;
        else
            return -p * W / (R * pow(Tmprtr, 2.));
    }

}
//Formula setting for thermal conductivity and it's derivative, Watt/(m*K)
double Lambda(double Tmprtr, const double a, int flag_phase)
{
    const double lambda_water = 0.56;
    const double lambda_steam = 0.05;
    const double T_boiling = 373;
    if (a < pow(10, -8))
    {
        if (Tmprtr <= T_boiling)
            return lambda_water;
        else
            return lambda_steam;
    }
    else
    {
        if (Tmprtr <= T_boiling && flag_phase == 0)
            LinIntWConductivity(Tmprtr);
        else
        {
            //Introduce thermal conductivity INSTEAD of thermal diffusivity
            const double T_star = 647.096;
            double L0 = 2.443221 * pow(10, -3.);
            double L1 = 1.323095 * pow(10, -2.);
            double L2 = 6.770357 * pow(10, -3.);
            double L3 = -3.454586 * pow(10, -3.);
            double L4 = 4.096266 * pow(10, -4.);
            double lambda_star = pow(10, -3);
            double teta = Tmprtr / T_star;
            return lambda_star * (sqrt(teta) / (L0 + L1 / teta + L2 / pow(teta, 2.) + L3 / pow(teta, 3.) + L4 / pow(teta, 4.)));
        }
    }
}
//Watt/(m*K^2)
double DfLambda(double Tmprtr, const double a, int flag_phase)
{
    const double T_boiling = 373;
    if (a < pow(10, -8))
        return 0;
    else
    {
        if (Tmprtr <= T_boiling && flag_phase == 0)
            return 0;
        else
        {
            //Introduce thermal conductivity INSTEAD of thermal diffusivity
            const double T_star = 647.096;
            double L0 = 2.443221 * pow(10, -3.);
            double L1 = 1.323095 * pow(10, -2.);
            double L2 = 6.770357 * pow(10, -3.);
            double L3 = -3.454586 * pow(10, -3.);
            double L4 = 4.096266 * pow(10, -4.);
            double lambda_star = pow(10, -3);
            double teta = Tmprtr / T_star;
            double bracket_sum = L0 + L1 / teta + L2 / pow(teta, 2.) + L3 / pow(teta, 3.) + L4 / pow(teta, 4.);
            return lambda_star / (pow(bracket_sum, 2.0) * sqrt(teta)) * (0.5 * bracket_sum + L1 / teta + 2. * L2 / pow(teta, 2.)
                + 3. * L3 / pow(teta, 3.) + 4. * L4 / pow(teta, 4.)) / T_star;
        }
    }
}
//Setting arrays of parametres depending on Tmprtr:
void ArraysParameters(const double a, const int N, double *r, double *T_next,
    const int s, vector <double>& ACp, vector <double>& ADensity, vector <double>& ALambda, int n, double inter, double Ti, double p)
{
    //Define Parameters
    ACp.resize(N + 1);
    //ADfCp.resize(N + 1);
    ADensity.resize(N + 1);
    //ADfDensity.resize(N + 1);
    ALambda.resize(N + 1);
    //ADfLambda.resize(N + 1);

    for (int i = 0; i < s + 1; i++) {
        //cout << r[i] << endl;
        ACp[i] = Cp(T_next[i], a);
        //ADfCp[i] = DfCp(T_next[i], a);
        ADensity[i] = Density(T_next[i], a);
        //ADfDensity[i] = DfDensity(T_next[i], a);
        ALambda[i] = Lambda(T_next[i], a, 0);
        //ADfLambda[i] = DfLambda(T_next[i], a, 0);
    }
    for (int i = s + 1; i < N + 1; i++) {
        //cout << r[i] << endl;
        ACp[i] = Cp(T_next[i], a);
        //ADfCp[i] = DfCp(T_next[i], a);
        ADensity[i] = Density(T_next[i], a);
        //ADfDensity[i] = DfDensity(T_next[i], a);
        ALambda[i] = Lambda(T_next[i], a, 1);
        //ADfLambda[i] = DfLambda(T_next[i], a, 1);
    }

    //Plot graphics of Parameters as a function of Temperature
    ofstream Parameters;
    Parameters.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/Data_new/Parameters/Parameters_" + to_string(n) + ".dat");
    Parameters << "TITLE=\"" << "Graphics" << "\"" << endl;
    Parameters << R"(VARIABLES= "rj, m", "T, K", "Cp, J/kg*K", "DfCp, J/kg*K^2", "rho, kg/m^3", "Dfrho, kg/m^3*K", "Lambda, Watt/m*K", "DfLambda, Watt/m*K^2")" << "\n";
    int j = 0;
    for (j; j < s + 1; j++) {
        Parameters << r[j] << " " << T_next[j] << " " << ACp[j] << " " << ADensity[j] << " " << ALambda[j] << "\n";
    }
    //Interface
    Parameters << inter << " " << Ti << " " << ACp[s] << " " << ADensity[s] << " " << ALambda[s] << "\n";
    for (j + 1; j < N + 1; j++) {
        Parameters << r[j] << " " << T_next[j] << " " << ACp[j] << " " << ADensity[j] << " "<< ALambda[j] << "\n";
    }
    Parameters.close();
}
void GraphicsSystEqu (int n, double* r, double* T_next, vector<double>& ALambda, vector<double>& ADensity,
    const int N, const int s, double dM, double inter, double Ti)
{
    ofstream OutCurrentTemp;
    ofstream OutFlow;
    cout << "iout" << n << "\n";
    if (n % 10 == 0 || n < 40)
    {
        OutCurrentTemp.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/Data_new/Temp/Temp_" + to_string(n) + ".dat");
        OutCurrentTemp << "TITLE=\"" << "Graphics" << "\"" << "\n";
        OutCurrentTemp << R"(VARIABLES= "rj, m", "T, K", "Lambda, W/m*K" )" << "\n";
        int i = 0;
        for (i; i < s + 1; i++)
            OutCurrentTemp << r[i] << " " << T_next[i] << " " << ALambda[i] << "\n";
        //Interface
        OutCurrentTemp << inter << " " << Ti << " " << ALambda[s] << "\n";
        for (i + 1; i < N + 1; i++)
            OutCurrentTemp << r[i] << " " << T_next[i] << " " << ALambda[i] << "\n";
        OutCurrentTemp.close();
    }
    //Get flow movement
    OutFlow.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/Data_new/Flow/Flow_" + to_string(n) + ".dat");
    OutFlow << "TITLE=\"" << "Graphics" << "\"" << "\n";
    OutFlow << R"(VARIABLES= "rj, m", "q, W", "rho, kg/m^3", "u_r, m/s" )" << "\n";
    //Поток газа выводим
    double u_r;
    for (int j = s + 1; j < N; j++)
    {
        u_r = dM / (pow(r[j], 2.0) * ADensity[j]);
        OutFlow << r[j] << " " << ALambda[j] * (T_next[j] - T_next[j - 1]) / (r[j] - r[j - 1]) << " "
            << ADensity[j] << " " << u_r << "\n";
    }
    OutFlow.close();
}
void OpeningFiles(ofstream& OutSurfaceFlow, ofstream& OutLastStepNevyazka, ofstream& OutFirstStepNevyazka)
{
    //Collecting mass flow on interface
    OutSurfaceFlow.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/Data_new/SurfaceFlow.dat");
    OutSurfaceFlow << "TITLE=\"" << "Graphics" << "\"" << "\n";
    OutSurfaceFlow << R"(VARIABLES= "t, s", "T0, K", "Ts, K", "Ti, K", "Ts+1, K", "gradT_d, K/m", "gradT_g, K/m", "q_d, W", "q_g, W", "dM, kg/s",
 "rho, kg/m^3", "u_r, m/s", "p, m/(cell size) ", "inter, m" )" << "\n";
    //Collecting Last Step Nevyazka
    OutLastStepNevyazka.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/Data_new/LastStepNevyazka.dat");
    OutLastStepNevyazka << "TITLE=\"" << "Graphics" << "\"" << endl;
    OutLastStepNevyazka << R"(VARIABLES= "t, s", "F")" << endl;
    //Collecting First Step Nevyazka
    OutFirstStepNevyazka.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/Data_new/FirstStepNevyazka.dat");
    OutFirstStepNevyazka << "TITLE=\"" << "Graphics" << "\"" << endl;
    OutFirstStepNevyazka << R"(VARIABLES= "t, s", "F")" << endl;
}

void F_right(vector<double>& f, double *T_vect, double *r, vector <double>& ACp, vector <double>& ADensity, vector <double>& ALambda,
    const int N, const double a, double dM, int s, double inter, const double Ti, double p, const double L_d)
{
    double F;
    //Nevyazka on the left border
    int l_f_num = 0;
    F = 1. / (r[1] - r[0]) * pow((r[l_f_num + 1] + r[l_f_num]) / 2., 2.) * 0.5 * (ALambda[l_f_num + 1] + ALambda[l_f_num])
        * (T_vect[l_f_num + 1] - T_vect[l_f_num]) / (r[l_f_num + 1] - r[l_f_num]);
    f[l_f_num] = F / (ACp[l_f_num] * ADensity[l_f_num] * pow(r[l_f_num + 1], 2.) / 24.); 
    cout << "f[0]" << f[0] << "\n";

    //Nevyazka in Water
    for (int j = 1; j < s; j++)
    {
        F = 2. / (r[j + 1] - r[j - 1]) * (pow((r[j + 1] + r[j]) / 2., 2.) * 0.5 * (ALambda[j + 1] + ALambda[j])
            * (T_vect[j + 1] - T_vect[j]) / (r[j + 1] - r[j])
            - pow((r[j] + r[j - 1]) / 2., 2.) * 0.5 * (ALambda[j] + ALambda[j - 1])
            * (T_vect[j] - T_vect[j - 1]) / (r[j] - r[j - 1]));
        f[j] = F / (ACp[j] * ADensity[j] * pow(r[j], 2.));
        //cout << j << " f1_" << j << " " << f[j] << "\n";
    }
    //Nevyazka in Gaze
    for (int j = s + 2; j < N; j++)
    {
        F = 2. / (r[j + 1] - r[j - 1]) * (pow((r[j + 1] + r[j]) / 2., 2.) * 0.5 * (ALambda[j + 1] + ALambda[j])
            * (T_vect[j + 1] - T_vect[j]) / (r[j + 1] - r[j])
            - pow((r[j] + r[j - 1]) / 2., 2.) * 0.5 * (ALambda[j] + ALambda[j - 1])
            * (T_vect[j] - T_vect[j - 1]) / (r[j] - r[j - 1]))
            + ACp[j] * dM * (T_vect[j] - T_vect[j - 1]) / (r[j] - r[j - 1]);
        f[j] = F / (ACp[j] * ADensity[j] * pow(r[j], 2.));
    }
    //cout << "inter= " << inter << "\n";
    //Near interface Nevyazka
    double term1, term2, term3;
    double h = r[s + 1] - r[s];
    cout << "Ti" << Ti << "p" << p << "\n";

    //To the left of the interface
    double rs_avg = (r[s] + r[s - 1]) / 2.0;
    double rs_diff = inter - rs_avg;
    term1 = (1.0 / rs_diff) * (Lambda(Ti, a, 0) * pow(inter, 2.) * ((p / ((p + 1) * h)) * T_vect[s - 2] - ((p + 1) / (p * h)) * T_vect[s - 1]
        + ((2 * p + 1) / ((p + 1) * p * h)) * Ti) - pow(rs_avg, 2.) * 0.5 * (ALambda[s] + ALambda[s - 1]) * (T_vect[s] - T_vect[s - 1]) / (r[s] - r[s - 1]));
    term3 = ACp[s] * ADensity[s] * pow(r[s], 2.);
    f[s] = term1 / term3;

    //On the interface
    term1 = pow(inter, 2.) * Lambda(Ti, a, 1) * ((2.0 * p - 7.0) / ((p - 3.0) * (p - 4.0) * h) * Ti + (p - 4.0) / ((p - 3.0) * h) * T_vect[s + 2]
        + (p - 3.0) / ((4.0 - p) * h) * T_vect[s + 3]);
    term2 = pow(inter, 2.) * Lambda(Ti, a, 0) * (p / ((p + 1.0) * h) * T_vect[s - 2] - (p + 1.0) / (p * h) * T_vect[s - 1]
        + (2.0 * p + 1.0) / ((p + 1.0) * p * h) * Ti);
    f[N + 1] = term1 - term2 + L_d * dM;
    //To the right of the interface
    double rs_avg_plus = (r[s + 1] + r[s + 2]) / 2.0;
    double rs_diff_plus = rs_avg_plus - inter;
    term1 = (1.0 / rs_diff_plus) * (pow(rs_avg_plus, 2.) * 0.5 * (ALambda[s + 2] + ALambda[s + 1]) * (T_vect[s + 2] - T_vect[s + 1]) / (r[s + 2] - r[s + 1])
        - pow(inter, 2.) * Lambda(Ti, a, 1) * ((2 * p - 7) / ((p - 3) * (p - 4) * h) * Ti + (p - 4) / ((p - 3) * h) * T_vect[s + 2]
            + (p - 3) / ((4 - p) * h) * T_vect[s + 3]));
    term2 = -ACp[s + 1] * dM * ((T_vect[s + 1] - Ti) / (r[s + 1] - inter));
    term3 = ACp[s + 1] * ADensity[s + 1] * pow(r[s + 1], 2.);
    f[s + 1] = (term1 + term2) / term3;
}

void F_right_test(vector<double>& f, double* T_vect, double* r, vector <double>& ACp, vector <double>& ADensity, vector <double>& ALambda,
    const int N, const double a, double dM, int s, double inter, const double Ti, double p, const double L_d)
{
    double F;

    //Nevyazka in Water
    for (int j = 1; j < s + 1; j++)
    {
        F = 2. / (r[j + 1] - r[j - 1]) * (pow((r[j + 1] + r[j]) / 2., 2.) * 0.5 * (ALambda[j + 1] + ALambda[j])
            * (T_vect[j + 1] - T_vect[j]) / (r[j + 1] - r[j])
            - pow((r[j] + r[j - 1]) / 2., 2.) * 0.5 * (ALambda[j] + ALambda[j - 1])
            * (T_vect[j] - T_vect[j - 1]) / (r[j] - r[j - 1]));
        f[j] = F / (ACp[j] * ADensity[j] * pow(r[j], 2.));
        //cout << j << " f1_" << j << " " << f[j] << "\n";
    }
    //Nevyazka in Gaze
    for (int j = s + 1; j < N; j++)
    {
        F = 2. / (r[j + 1] - r[j - 1]) * (pow((r[j + 1] + r[j]) / 2., 2.) * 0.5 * (ALambda[j + 1] + ALambda[j])
            * (T_vect[j + 1] - T_vect[j]) / (r[j + 1] - r[j])
            - pow((r[j] + r[j - 1]) / 2., 2.) * 0.5 * (ALambda[j] + ALambda[j - 1])
            * (T_vect[j] - T_vect[j - 1]) / (r[j] - r[j - 1]))
            + ACp[j] * dM * (T_vect[j] - T_vect[j - 1]) / (r[j] - r[j - 1]);
        f[j] = F / (ACp[j] * ADensity[j] * pow(r[j], 2.));
    }
    //cout << "inter= " << inter << "\n";
    /*
    //Near interface Nevyazka
    double term1, term2, term3;
    double h = r[s + 1] - r[s];
    cout << "Ti" << Ti << "p" << p << "\n";

    //To the left of the interface
    double rs_avg = (r[s] + r[s - 1]) / 2.0;
    double rs_diff = inter - rs_avg;
    term1 = (1.0 / rs_diff) * (Lambda(Ti, a, 0) * pow(inter, 2.) * ((p / ((p + 1) * h)) * T_vect[s - 2] - ((p + 1) / (p * h)) * T_vect[s - 1]
        + ((2 * p + 1) / ((p + 1) * p * h)) * Ti) - pow(rs_avg, 2.) * 0.5 * (ALambda[s] + ALambda[s - 1]) * (T_vect[s] - T_vect[s - 1]) / (r[s] - r[s - 1]));
    term3 = ACp[s] * ADensity[s] * pow(r[s], 2.);
    f[s] = term1 / term3;

    //On the interface
    term1 = pow(inter, 2.) * Lambda(Ti, a, 1) * ((2.0 * p - 7.0) / ((p - 3.0) * (p - 4.0) * h) * Ti + (p - 4.0) / ((p - 3.0) * h) * T_vect[s + 2]
        + (p - 3.0) / ((4.0 - p) * h) * T_vect[s + 3]);
    term2 = pow(inter, 2.) * Lambda(Ti, a, 0) * (p / ((p + 1.0) * h) * T_vect[s - 2] - (p + 1.0) / (p * h) * T_vect[s - 1]
        + (2.0 * p + 1.0) / ((p + 1.0) * p * h) * Ti);
    f[N] = term1 - term2 + L_d * dM;

    //To the right of the interface
    double rs_avg_plus = (r[s + 1] + r[s + 2]) / 2.0;
    double rs_diff_plus = rs_avg_plus - inter;
    term1 = (1.0 / rs_diff_plus) * (pow(rs_avg_plus, 2.) * 0.5 * (ALambda[s + 2] + ALambda[s + 1]) * (T_vect[s + 2] - T_vect[s + 1]) / (r[s + 2] - r[s + 1])
        - pow(inter, 2.) * Lambda(Ti, a, 1) * ((2 * p - 7) / ((p - 3) * (p - 4) * h) * Ti + (p - 4) / ((p - 3) * h) * T_vect[s + 2]
            + (p - 3) / ((4 - p) * h) * T_vect[s + 3]));
    term2 = -ACp[s + 1] * dM * ((T_vect[s + 1] - Ti) / (r[s + 1] - inter));
    term3 = ACp[s + 1] * ADensity[s + 1] * pow(r[s + 1], 2.);
    f[s + 1] = (term1 + term2) / term3;
    */
}

double f(double q, const double N, const double h_min, const double x_r, const double x_l) //возвращает значение функции f(q) = ...

{
    return h_min * pow(q, N) - (x_r - x_l) * q - h_min + (x_r - x_l);
}

double df(double q, const double N, const double h_min, const double x_r, const double x_l) //возвращает значение производной

{
    return h_min * N * pow(q, N - 1) - (x_r - x_l);
}

double d2f(double q) // значение второй производной

{
    return 1.;
}

double DefineQ(const double N, const double h_min, const double x_r, const double x_l)
{
    int i = 0;//переменные для расчета итерации
    double x0, xn = 0;// вычисляемые приближения для корня
    double a, b;// границы отрезка, между которыми находится решение q
    a = 1;
    b = 1.2;
    cout << f(a, N, h_min, x_r, x_l) << " " << f(b, N, h_min, x_r, x_l) << endl;
    if (f(a, N, h_min, x_r, x_l) * f(b, N, h_min, x_r, x_l) > 0) // если знаки функции на краях отрезка одинаковые, то здесь нет корня
        cout << "\nError! No roots in this interval\n";
    else
    {
        if (f(a, N, h_min, x_r, x_l) * d2f(a) > 0) x0 = a; // для выбора начальной точки проверяем f(x0)*d2f(x0)>0 ?
        else x0 = b;
        xn = x0 - f(x0, N, h_min, x_r, x_l) / df(x0, N, h_min, x_r, x_l); // считаем первое приближение
        cout << ++i << "-th iteration = " << xn << "\n";
        while (fabs(x0 - xn) > pow(10, -8)) // пока не достигнем необходимой точности, будет продолжать вычислять
        {
            x0 = xn;
            xn = x0 - f(x0, N, h_min, x_r, x_l) / df(x0, N, h_min, x_r, x_l); // непосредственно формула Ньютона
            cout << ++i << "-th iteration = " << xn << "\n";
        }
        cout << "\nRoot = " << xn; // вывод вычисленного корня
    }
    std::cout << "\nHello World!\n";
    return xn;
}

void InitialGrid(int& N, const double x_l, const double x_r, const double T_l, const double T_r, vector<double>& r,
    vector<double>& T_cur, vector<double>& T_next, double h, const double q, const double h_min, const int big_number, const int const_params,
    const int Nd, const int N_uni, const int N_uni_near_center)
{
    T_cur.resize(N + 2);
    T_next.resize(N + 2);
    r.resize(N + 1);
    ofstream OutX;
    OutX.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/Data_new/X_grid.dat");
    OutX << "TITLE=\"" << "Graphics" << "\"" << endl;
    OutX << R"(VARIABLES= "j", "hj, m" )" << endl;
    //Non-uniform 1-dimensional grid 
    if (const_params == 0)
    {
        if (q < 1)
        {
            double h_cell = h;
            while (h_cell >= h_min)
            {
                h_cell = h_cell * q;
                cout << h_cell << endl;
                N = N + 1;
            }
            cout << N << endl;
        }

        if (q > 1 || q == 1)
        {
            h = h_min;
            r[0] = x_l;
            T_cur[0] = T_l;
            OutX << 0 << " " << h << " " << "\n";
            cout << "x0" << " " << r[0] << endl;
            for (int j = 1; j < N + 1; j++)
            {
                r[j] = r[j - 1] + h;
                h = h * q;
                if (r[j] <= x_l + 0.001)
                    T_cur[j] = T_l;
                else
                    T_cur[j] = T_r;
                cout << "T" << j << " " << T_cur[j] << endl;
                OutX << j << " " << h << " " << "\n";
                cout << "h" << h << " r" << j << " " << r[j] << endl;
            }
        }
        if (q < 1)
        {
            r[N] = x_r;
            T_cur[N] = T_r;
            cout << "N" << N << endl;
            for (int j = N - 1; j >= 0; j--)
            {
                cout << "T" << j + 1 << " " << T_cur[j + 1] << endl;
                OutX << j + 1 << " " << h << " " << "\n";
                cout << "h" << h << " r" << j + 1 << " " << r[j + 1] << endl;
                r[j] = r[j + 1] - h;
                h = h * q;
                T_cur[j] = T_r;
            }
            T_cur[0] = T_l;
            cout << "T" << 0 << " " << T_cur[0] << endl;
            OutX << 0 << " " << h << " " << "\n";
            cout << "h" << h << " r" << 0 << " " << r[0] << endl;
            cout << "Area size = " << r[N] - r[0] << endl;
        }
    }
    //Uniform for droplet and non-uniform for gas 1-dimensional grid 
    else if (const_params == 1)
    {
        if (q > 1 || q == 1)
        {
            //Setting the value on the boundaries and initial distribution
            r[0] = x_l;
            T_cur[0] = T_l;
            cout << "r0" << " " << r[0] << "\n";
            h = h_min;
            //r[1] = r[0] + h / 2.;
            //T_cur[1] = T_l;
            //cout << "r1" << " " << r[1] << "\n";
            int j = 1;
            //Uniform grid
            /*
            for (j; j < N_uni_near_center + 1; j++)
            {
                r[j] = r[j - 1] + h;
                if (j < Nd + 1)
                    T_cur[j] = T_l;
                else
                    T_cur[j] = T_r;
                //cout << "T" << j << " " << T_cur[j] << endl;
                OutX << j << " " << h << " " << "\n";
                cout << "h" << h << " r" << j << " " << r[j] << endl;
            }
            h = 2 * h;
            */
            for (j; j < N_uni + 1; j++)
            {
                r[j] = r[j - 1] + h;
                if (j < Nd + 1)
                    T_cur[j] = T_l;
                else
                    T_cur[j] = T_r;
                //cout << "T" << j << " " << T_cur[j] << endl;
                OutX << j << " " << h << " " << "\n";
                cout << "h" << h << " r" << j << " " << r[j] << " T= " << T_cur[j] << "\n";
            }
            //Non-uniform grid
            for (j + 1; j < N + 1; j++)
            {
                r[j] = r[j - 1] + h;
                h = h * q;
                T_cur[j] = T_r;
                //cout << "T" << j << " " << T_cur[j] << endl;
                OutX << j << " " << h << " " << "\n";
                cout << "h" << h << " r" << j << " " << r[j] << " T= " << T_cur[j] << "\n";
            }
        }
    }
    OutX.close();
}
/*
 *--------------------------------------------------------------------
 * Main Program
 *--------------------------------------------------------------------
 */
int Integrate_IDA(int N, vector<double>& r_vect, vector<double>& T_vect, vector <double>& ACp, vector <double>& ADensity, vector <double>& ALambda,
    vector<double>& nevyaz, double a, double& dM, int s, double inter, double Ti, double p, double L_d)
{
    void* mem;
    N_Vector yy, yp, avtol;
    realtype rtol, * yval, * ypval, * atval;
    realtype t0, tout1, tout, tret;
    int iout, retval, retvalr;
    SUNMatrix A;
    SUNLinearSolver LS;
    SUNNonlinearSolver NLS;
    SUNContext ctx;

    UserData* data = NULL;
    data = (UserData*)malloc(sizeof(UserData));
    if (data == NULL) {
        cout << "handle allocation failure" << "\n";
    }

    mem = NULL;
    yy = yp = avtol = NULL;
    yval = ypval = atval = NULL;
    A = NULL;
    LS = NULL;
    NLS = NULL;

    cout << "hey" << "\n";
    // Define the number of equations
    int NEQ = N + 2;
    int k = 0;

    data->outNevyaz = new ofstream;
    data->N = N;
    data->kp = k;
    data->r = new realtype[N + 2];
    data->T = new realtype[N + 2];
    data->NEQ = NEQ;
    data->Tr = T_vect[N];
    data->M = dM;
    data->s = s;
    data->a = a;
    data->inter = inter;
    data->p = p;
    data->L_d = L_d;

    for (int i = 0; i < N + 1; i++) {
        data->r[i] = r_vect[i];
        data->T[i] = T_vect[i];
        //cout << i << " = " << Ith(res_vect, i) << endl;
    }
    data->T[N + 1] = Ti;
    cout << "N= " << N << "p= " << p << "inter= " << inter << "L_d= " << L_d << "\n";

    /* Create SUNDIALS context */
    retval = SUNContext_Create(NULL, &ctx);
    if (check_retval(&retval, "SUNContext_Create", 1)) return(1);

    /* Allocate N-vectors. */
    yy = N_VNew_Serial(NEQ, ctx);
    if (check_retval((void*)yy, "N_VNew_Serial", 0)) return(1);
    yp = N_VClone(yy);
    if (check_retval((void*)yp, "N_VNew_Serial", 0)) return(1);
    avtol = N_VClone(yy);
    if (check_retval((void*)avtol, "N_VNew_Serial", 0)) return(1);

    /* Create and initialize  y, y', and absolute tolerance vectors. */
    yval = N_VGetArrayPointer(yy);
    ypval = N_VGetArrayPointer(yp);
    //Error rate
    rtol = RCONST(1.0e-3);
    atval = N_VGetArrayPointer(avtol);

    /* Integration limits by time */
    t0 = ZERO;
    tout1 = RCONST(0.01);

    //cout << NEQ << " = " << Ith(res_vect, NEQ + 1) << endl;
    for (int i = 0; i < N; i++) {
        Ith(avtol, i) = RCONST(1.0);
        Ith(yy, i) = T_vect[i];
        cout << "yy " << i << " = " << Ith(yy, i) << "\n";
    }
    Ith(avtol, N) = RCONST(1.0);
    Ith(yy, N) = data->T[N];
    Ith(avtol, N + 1) = RCONST(1.0);
    Ith(yy, N + 1) = data->T[N + 1];
    Ith(avtol, N + 2) = RCONST(1.0);
    //cout << "yy " << N << " = " << Ith(yy, N) << "\n";
    //Calculate Nevyazka's
    ArraysParameters(a, N, data->r, data->T, s, ACp, ADensity, ALambda, 0, inter, Ti, p);
    F_right_test(nevyaz, data->T, data->r, ACp, ADensity, ALambda, N, a, dM, s, inter, Ti, p, L_d);
    //Opening file for Nevyazka
    cout << "Отладочная информация: адрес переменной = " << &data << endl;
    // Открываем файловый поток и сохраняем его в структуре
    data->outNevyaz->open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/Data_new/Nevyaz.dat");

    // Проверяем, удалось ли открыть файл
    if (!data->outNevyaz->is_open()) {
        cout << "Ошибка открытия файла" << endl;
        return 0;
    }

    cout << "Отладочная информация: data=" << data << endl;
    *(data->outNevyaz) << "TITLE=\"" << "Graphics" << "\"" << endl;
    *(data->outNevyaz) << R"(VARIABLES= "i", "F")" << endl;
    //Calculate modul of Nevyazka
    double mod_nevyaz = 0;
    for (int i = 1; i < N; i++)
        mod_nevyaz += pow(nevyaz[i], 2.0);
    mod_nevyaz = pow(mod_nevyaz, 0.5);
    *(data->outNevyaz) << 0 << " " << mod_nevyaz << "\n";

    for (int i = 0; i < N; i++) {
        Ith(yp, i) = nevyaz[i];
        cout << "Ithypi = " << Ith(yp, i) << "\n";
    }
    Ith(yp, N) = nevyaz[N];
    Ith(yp, N + 1) = 0;
    cout << "Ithypi = " << Ith(yp, N) << "\n";

    PrintHeader(rtol, avtol, yy);

    /* Call IDACreate and IDAInit to initialize IDA memory */
    mem = IDACreate(ctx);
    if (check_retval((void*)mem, "IDACreate", 0)) return(1);

    retval = IDAInit(mem, resrob, t0, yy, yp);
    if (check_retval(&retval, "IDAInit", 1)) return(1);
    /* Call IDASVtolerances to set tolerances */

    retval = IDASVtolerances(mem, rtol, avtol);
    if (check_retval(&retval, "IDASVtolerances", 1)) return(1);

    retval = IDASetUserData(mem, data);
    if (check_retval(&retval, "IDASetUserData", 1)) return(1);
    retval = IDASetMaxNumSteps(mem, 20000);

    /* Create dense SUNMatrix for use in linear solves */
    A = SUNDenseMatrix(NEQ, NEQ, ctx);
    if (check_retval((void*)A, "SUNDenseMatrix", 0)) return(1);

    /* Create dense SUNLinearSolver object */
    LS = SUNLinSol_Dense(yy, A, ctx);
    if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) return(1);

    /* Attach the matrix and linear solver */
    retval = IDASetLinearSolver(mem, LS, A);
    if (check_retval(&retval, "IDASetLinearSolver", 1)) return(1);

    NLS = SUNNonlinSol_Newton(yy, ctx);
    if (check_retval((void*)NLS, "SUNNonlinSol_Newton", 0)) return(1);

    /* Attach the nonlinear solver */
    retval = IDASetNonlinearSolver(mem, NLS);
    if (check_retval(&retval, "IDASetNonlinearSolver", 1)) return(1);


    /* In loop, call IDASolve, print results, and test for error.
       Break out of loop when NOUT preset output times have been reached. */

    cout << "data->N" << data->N << "\n";

    iout = 0; tout = tout1;
    double tend = pow(10, 5);
    ofstream fout;
    while (iout < 10) {
        cout << "how is it going?" << "\n";
        retval = IDASolve(mem, tout, &tret, yy, yp, IDA_NORMAL);
        data->outNevyaz->close();
        cout << "how is it going2?" << "\n";
        ExportToArray(data->T, dM, data, yy, N);
        cout << "how is it going3?" << "\n";
        //PrintOutput(mem, tret, yy);
        cout << "t = " << tout << "\n";
        cout << "M = " << dM << "\n";

        //Ith(yp, N) = 0;
        fout.open("file" + to_string(tout * pow(10, 8)) + ".dat");
        fout << "TITLE=\"" << "Graphics" << "\"" << endl;
        fout << R"(VARIABLES= "x", "T", "lambda")" << endl;
        for (int i = 0; i < N; i++) {
            fout << r_vect[i] << "  " << T_vect[i] << "\n";
        }
        fout.close();
        //Collect information
        GraphicsSystEqu(iout, data->r, data->T, ALambda, ADensity, N, s, dM, inter, Ti);
        if (check_retval(&retval, "IDASolve", 1)) return(1);

        if (retval == IDA_SUCCESS) {
            iout++;
            tout += tout1;
        }
    }

    /* Print final statistics to the screen */
    printf("\nFinal Statistics:\n");
    retval = IDAPrintAllStats(mem, stdout, SUN_OUTPUTFORMAT_TABLE);

    /* Print final statistics to a file in CSV format */
   // FID = fopen("idaRoberts_dns_stats.csv", "w");
   // retval = IDAPrintAllStats(mem, FID, SUN_OUTPUTFORMAT_CSV);
   // fclose(FID);

    /* check the solution error */
    //retval = check_ans(yy, tret, rtol, avtol);

    /* Free memory */
    IDAFree(&mem);
    SUNNonlinSolFree(NLS);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
    N_VDestroy(avtol);
    N_VDestroy(yy);
    N_VDestroy(yp);
    SUNContext_Free(&ctx);

    delete(data->outNevyaz);
    delete[](data->r);
    delete[](data->T);
    free(data);

    return(retval);
}

void Solver(int N, vector<double>& r, vector<double>& T_cur, vector<double>& T_next,
    const double a, const double b, const double c, const int total_time, const int const_params,
    const double T_l, const double T_r, const int Nd, int m_max)
{
    double dM, q_g, q_d, Lambda_d, Lambda_g, Ti;
    const int N_minus = N - 1;
    int retval;
    vector <double> nevyaz(N + 2);
    vector <double> ACp, ADfCp, ADensity, ADfDensity, ALambda, ADfLambda;
    int flag_evap = 0;
    //Define initial surface position and mass flow on surface
    double p = 1.5;
    int s = Nd;
    double inter = r[s] + (p - 1) * (r[s + 1] - r[s]);
    //Define initial boiling tempereature for surface, K
    double T_cur_i = 300.;
    //double q = 20 * pow(10, 3);                  //W / m^2
    double L_d = 2258.2 * pow(10, 3.);              //J / kg 
    //const determination of dM
    //dM = 4 * PI * pow(inter, 2.0) * q / L_d;     // kg / s
    //dM equals zero because it is only warming time
    dM = 0.;                                       // kg / s
    //dM / (4 * PI);
    //Задаем начальную температуру следующего слоя как решение предыдущей
    Ti = T_cur_i;
    for (int j = 0; j < N + 1; j++)
        T_next[j] = T_cur[j];

    ofstream OutSurfaceFlow;
    ofstream OutLastStepNevyazka;
    ofstream OutFirstStepNevyazka;
    OpeningFiles(OutSurfaceFlow, OutLastStepNevyazka, OutFirstStepNevyazka);
    //TIME LOOP
    retval = Integrate_IDA(N, r, T_next, ACp, ADensity, ALambda, nevyaz, a, dM, s, inter, Ti, p, L_d);
}

int main(void)
{
    const int const_params = 1;          // 0 - Non-uniform grid, 1 - Mixed
    const double dt = 0.0001;
    const double dtau = 0.0001;
    int m_max = 1000;
    const int total_time = 200000;
    const int power_of = -3;              //const if a < pow(10, -8)
    const double a = pow(10, -(power_of + 3)), b = pow(10, -power_of), c = pow(10, -(power_of - 3));  //a = 1
    const double T_l = 300;
    const double T_r = 1500;
    double h = 0.2;
    //зададим минимальный возможный размер ячейки
    const double h_min = 0.0000243902439024;                                        //0.0005 * 2 / 41 = 0.0000243902439024= 24.3902439024 мкр;  
    //Number of cells
    int N = 100;
    const int Nd = 20;
    const int N_uni = 2 * Nd;
    const int N_uni_near_center = 5;
    const double x_l = 0.001;
    //const double R = 0.0005 м = 500 мкр
    const double x_r = 0.201; //зададим область = 20 см
    double x_uni = x_l + N_uni * h_min;
    double q = DefineQ(N - N_uni, h_min, x_r, x_uni);
    vector <double> r;
    vector <double> T_cur;
    vector <double> T_next;
    InitialGrid(N, x_l, x_r, T_l, T_r, r, T_cur, T_next, h, q, h_min, total_time, const_params, Nd, N_uni, N_uni_near_center);
    Solver(N, r, T_cur, T_next, a, b, c, total_time, const_params, T_l, T_r, Nd, m_max);
}


/*
 *--------------------------------------------------------------------
 * Functions called by IDA
 *--------------------------------------------------------------------
 */


 /*
  * Define the system residual function.
  */

static int resrob(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data)
{
    realtype* yval, * ypval, * rval;
    UserData *data;
    realtype* r_cells, * T_vect;

    //MyParams
    int s, k;
    double a, inter, p, dM, L_d;

    if (user_data == NULL) {
        cout << "handle null pointer error" << "\n";
    }
    data = (UserData*)user_data;

    T_vect = data->T;
    r_cells = data->r;
    s = data->s;
    a = data->a;
    inter = data->inter;
    p = data->p;
    dM = data->M;
    L_d = data->L_d;
    k = data->kp;

    int j;
    yval = N_VGetArrayPointer(yy);
    ypval = N_VGetArrayPointer(yp);
    rval = N_VGetArrayPointer(rr);
    //cout << "ypvalres0 = " << yval[0] << "\n";
    int myN = data->N;
    ExportToArray(T_vect, data->M, data, yy, data->N);
    vector <double> nevyaz(myN + 2);
    vector <double> ACp, ADensity, ALambda;
    ArraysParameters(a, myN, data->r, data->T, s, ACp, ADensity, ALambda, 0, inter, data->T[myN + 1], p);
    //Calculate Nevyazka's
    F_right_test(nevyaz, data->T, data->r, ACp, ADensity, ALambda, myN, a, dM, s, inter, data->T[myN + 1], p, L_d);
    for (j = 0; j < myN; j++) {
        rval[j] = nevyaz[j] - ypval[j];
        cout << "nevyaz[" << j << "]= " << nevyaz[j] << "\n";
        cout << "rval[" << j << "]= " << rval[j] << "\n";
    }
   
    cout << "nevyaz[0]= " << nevyaz[0] << "ypval[0]= " << ypval[0] << "\n";
    rval[myN] = 0;
    rval[myN + 1] = nevyaz[myN + 1];
    cout << "rval[myN + 1]= " << rval[myN + 1] << "\n";
    //Opening file for T on iteration
    k++;
    data->kp = k;
    //Calculate modul of Nevyazka
    double mod_nevyaz = 0;
    for (int i = 0; i < myN; i++)
        mod_nevyaz += pow(rval[i], 2.0);
    mod_nevyaz = pow(mod_nevyaz, 0.5);
    //
    *(data->outNevyaz) << k << " " << mod_nevyaz << "\n";
    //Define nevyazkas on each iteration
    ofstream OutIterNevyazka;
    OutIterNevyazka.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/Data_new/IterNevyazka/IterNevyazka_" + to_string(k) + ".dat");
    OutIterNevyazka << "TITLE=\"" << "Graphics" << "\"" << "\n";
    OutIterNevyazka << R"(VARIABLES= "rj, m", "F" )" << "\n";
    int i = 0;
    for (i; i < s + 1; i++)
        OutIterNevyazka << r_cells[i] << " " << rval[i] << "\n";
    //Interface
    //OutIterNevyazka << inter << " " << rval[myN] << "\n";
    for (i + 1; i < myN; i++)
        OutIterNevyazka << r_cells[i] << " " << rval[i] << "\n";
    OutIterNevyazka.close();
    //Define Temperature on each iteration
    ofstream OutIterCurrentTemp;
    OutIterCurrentTemp.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/Data_new/IterTemp/IterTemp_" + to_string(k) + ".dat");
    OutIterCurrentTemp << "TITLE=\"" << "Graphics" << "\"" << "\n";
    OutIterCurrentTemp << R"(VARIABLES= "rj, m", "T, K", "Lambda, W/m*K" )" << "\n";
    i = 0;
    for (i; i < s + 1; i++)
        OutIterCurrentTemp << r_cells[i] << " " << T_vect[i] << " " << ALambda[i] << "\n";
    //Interface
    OutIterCurrentTemp << inter << " " << T_vect[myN + 1] << " " << ALambda[s] << "\n";
    for (i + 1; i < myN + 1; i++)
        OutIterCurrentTemp << r_cells[i] << " " << T_vect[i] << " " << ALambda[i] << "\n";
    OutIterCurrentTemp.close();
    return(0);
}


/*
 *--------------------------------------------------------------------
 * Private functions
 *--------------------------------------------------------------------
 */

 /*
  * Print first lines of output (problem description)
  */

static void PrintHeader(realtype rtol, N_Vector avtol, N_Vector y)
{
    realtype* atval, * yval;

    atval = N_VGetArrayPointer(avtol);
    yval = N_VGetArrayPointer(y);

    printf("\nidaRoberts_dns: Robertson kinetics DAE serial example problem for IDA\n");
    printf("         Three equation chemical kinetics problem.\n\n");
    printf("Linear solver: DENSE, with user-supplied Jacobian.\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("Tolerance parameters:  rtol = %Lg   atol = %Lg %Lg %Lg \n",
        rtol, atval[0], atval[1], atval[2]);
    printf("Initial conditions y0 = (%Lg %Lg %Lg)\n",
        yval[0], yval[1], yval[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("Tolerance parameters:  rtol = %g   atol = %g %g %g \n",
        rtol, atval[0], atval[1], atval[2]);
    printf("Initial conditions y0 = (%g %g %g)\n",
        yval[0], yval[1], yval[2]);
#else
    printf("Tolerance parameters:  rtol = %g   atol = %g %g %g \n",
        rtol, atval[0], atval[1], atval[2]);
    printf("Initial conditions y0 = (%g %g %g)\n",
        yval[0], yval[1], yval[2]);
#endif
    printf("Constraints and id not used.\n\n");
    printf("-----------------------------------------------------------------------\n");
    printf("  t             y1           y2           y3");
    printf("      | nst  k      h\n");
    printf("-----------------------------------------------------------------------\n");
}

/*
 * Print Output
 */

static void PrintOutput(void* mem, realtype t, N_Vector y)
{
    realtype* yval;
    int retval, kused;
    long int nst;
    realtype hused;

    yval = N_VGetArrayPointer(y);

    retval = IDAGetLastOrder(mem, &kused);
    check_retval(&retval, "IDAGetLastOrder", 1);
    retval = IDAGetNumSteps(mem, &nst);
    check_retval(&retval, "IDAGetNumSteps", 1);
    retval = IDAGetLastStep(mem, &hused);
    check_retval(&retval, "IDAGetLastStep", 1);
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("%10.4Le %12.4Le %12.4Le %12.4Le | %3ld  %1d %12.4Le\n",
        t, yval[0], yval[1], yval[2], nst, kused, hused);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("%10.4e %12.4e %12.4e %12.4e | %3ld  %1d %12.4e\n",
        t, yval[0], yval[1], yval[2], nst, kused, hused);
#else
    printf("%10.4e %12.4e %12.4e %12.4e | %3ld  %1d %12.4e\n",
        t, yval[0], yval[1], yval[2], nst, kused, hused);
#endif
}

static void PrintRootInfo(int root_f1, int root_f2)
{
    printf("    rootsfound[] = %3d %3d\n", root_f1, root_f2);
    return;
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns an integer value so check if
 *            retval < 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
 */

static int check_retval(void* returnvalue, const char* funcname, int opt)
{
    int* retval;
    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && returnvalue == NULL) {
        fprintf(stderr,
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
        return(1);
    }
    else if (opt == 1) {
        /* Check if retval < 0 */
        retval = (int*)returnvalue;
        if (*retval < 0) {
            fprintf(stderr,
                "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
                funcname, *retval);
            return(1);
        }
    }
    else if (opt == 2 && returnvalue == NULL) {
        /* Check if function returned NULL pointer - no memory allocated */
        fprintf(stderr,
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
        return(1);
    }

    return(0);
}

/* compare the solution at the final time 4e10s to a reference solution computed
   using a relative tolerance of 1e-8 and absoltue tolerance of 1e-14 */
static int check_ans(N_Vector y, realtype t, realtype rtol, N_Vector atol)
{
    int      passfail = 0;        /* answer pass (0) or fail (1) retval */
    N_Vector ref;               /* reference solution vector        */
    N_Vector ewt;               /* error weight vector              */
    realtype err;               /* wrms error                       */

    /* create reference solution and error weight vectors */
    ref = N_VClone(y);
    ewt = N_VClone(y);

    /* set the reference solution data */
    //NV_Ith_S(ref, 0) = RCONST(5.2083474251394888e-08);
    //NV_Ith_S(ref, 1) = RCONST(2.0833390772616859e-13);
    //NV_Ith_S(ref, 2) = RCONST(9.9999994791631752e-01);

    /* compute the error weight vector, loosen atol */
    N_VAbs(ref, ewt);
    N_VLinearSum(rtol, ewt, RCONST(10.0), atol, ewt);
    if (N_VMin(ewt) <= ZERO) {
        fprintf(stderr, "\nSUNDIALS_ERROR: check_ans failed - ewt <= 0\n\n");
        return(-1);
    }
    N_VInv(ewt, ewt);

    /* compute the solution error */
    N_VLinearSum(ONE, y, -ONE, ref, ref);
    err = N_VWrmsNorm(ref, ewt);

    /* is the solution within the tolerances? */
    passfail = (err < ONE) ? 0 : 1;

    if (passfail) {
        //fprintf(stdout, "\nSUNDIALS_WARNING: check_ans error=%"GSYM"\n\n", err);
    }

    /* Free vectors */
    N_VDestroy(ref);
    N_VDestroy(ewt);

    return(passfail);
}
