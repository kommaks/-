
#include <iostream>
#include <fstream>
#include <C:\Users\user\source\Научная работа\DropletEvaporation(IDA)\consts.cpp>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <cmath>
using namespace std;

#include <ida/ida.h>   
#include <kinsol/kinsol.h>             /* access to KINSOL func., consts. */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector       */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix       */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype */
#include <sunnonlinsol/sunnonlinsol_newton.h> /* access to Newton SUNNonlinearSolver  */
/* Problem Constants */


//#define FTOL   RCONST(1.e-12) /* function tolerance */
//#define STOL   RCONST(1.e-12) /* step tolerance     */

#define ZERO   RCONST(0.0)
#define PT25   RCONST(0.25)
#define PT5    RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)
#define TWO    RCONST(2.0)

#define PI     RCONST(3.1415926)
#define E      RCONST(2.7182818)


typedef struct {
    realtype* x;
    realtype* T;
    realtype* Y_H2O;
    int Nx;
    int N_m;
    int NEQ;
    int Np_inter;
    int n_tout;
    realtype pp;
    realtype L_dp;
    realtype dMp;
    realtype ap;
    realtype Tl;
    realtype Tr;
    realtype M;
    realtype Tp_inter;
    ofstream* outNevyaz;
    int kp;
    void* mykmem;
} *UserData;

/* Accessor macro */
//#define Ith(v,i)    NV_Ith_S(v,i-1)
//#define IJth(A,i,j) SM_ELEMENT_D(A,i-1,j-1)

/* Accessor macro */
#define Ith(v,i)    NV_Ith_S(v,i)
#define IJth(A,i,j) SM_ELEMENT_D(A,i,j)

/* Functions Called by the IDA Solver */
static int resrob(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data);

void ExportToArray(vector<double>& T_vect, double& M, UserData data, N_Vector yy, int N_x)
{
    //cout << "data->Tl; = " << data->Tl << "\n";
    //T_vect[0] = data->Tl;
    for (int i = 0; i < N_x - 1; i++)
    {
        T_vect[i] = Ith(yy, i);
    }
    T_vect[N_x - 1] = data->Tr;
    //cout << "MyNx = " << myNx << "\n";
    //cout << "MyNm = " << myNm << "\n";
    /*
    M = Ith(yy, data->N_m + 1);
    j = 1;
    Y_vect[0] = 0.;
    Y_vect[0] = 0.;tout
    //cout << "M = " << M << "\n";
    for (int i = data->N_m + 1; i < data->NEQ; i++)
    {
        Y_vect[j] = Ith(yy, i + 1);
        //cout << "Y_vect  " << j << " =  " << Y_vect[j] << endl;
        j++;
    }
    Y_vect[j] = Y_vect[j - 1];
    */
}

void ExportToArray(double* T_vect, double& M, UserData data, N_Vector yy, int N_x)
{
    //cout << "data->Tl; = " << data->Tl << "\n";
    //T_vect[0] = data->Tl;
    for (int i = 0; i < N_x - 1; i++)
    {
        T_vect[i] = Ith(yy, i);
    }
    T_vect[N_x - 1] = data->Tr;
    //cout << "MyNx = " << myNx << "\n";
    //cout << "MyNm = " << myNm << "\n";
    //data->M = Ith(yy, data->N_m + 1);
    /*
    j = 1;
    Y_vect[0] = 0.;
    //cout << "M = " << M << "\n";
    for (int i = data->N_m + 1; i < data->NEQ; i++)
    {
        Y_vect[j] = Ith(yy, i + 1);
        //cout << "Y_vect  " << j << " =  " << Y_vect[j] << endl;
        j++;
    }
    Y_vect[j] = Y_vect[j - 1];
    */
}

/* Private Helper Functions */
static int check_retval(void* retvalvalue, const char* funcname, int opt);

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
    NV_Ith_S(ref, 0) = RCONST(5.2083474251394888e-08);
    NV_Ith_S(ref, 1) = RCONST(2.0833390772616859e-13);
    NV_Ith_S(ref, 2) = RCONST(9.9999994791631752e-01);

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

///////////////////////////////////////////////////////////////
//MyFunctions
///////////////////////////////////////////////////////////////
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
//Formula setting for Heat Capacity and it's derivative, J/(kg*K)= m^2/(s^2*K)
double Cp(double T, const double a)
{
    double Res;
    const double Cp_water = 4180.;
    const double Cp_steam = 2000.;
    const double T_boiling = 373;
    if (T <= T_boiling)
        Res = Cp_water;
    else
        Res = Cp_steam;
    /*
    //Molar mass of water
    const double W = 18 * pow(10, -3);
    double AA, B, C, D, EE;
    //Polinom for T = (K) /1000
    double T_forCp = T / 1000.;
    if (T <= T_boiling)
    {
        AA = -203.6060;
        B = 1523.290;
        C = -3196.413;
        D = 2474.455;
        EE = 3.855326;
    }
    else
    {
        AA = 30.09200;
        B = 6.832514;
        C = 6.793435;
        D = -2.534480;
        EE = 0.082139;
    }
    //Get value
    if (a < pow(10, -8))
    {

    }
    else
        Res = (AA + B * T_forCp + C * pow(T_forCp, 2.) + D * pow(T_forCp, 3.) + EE * pow(T_forCp, -2.)) / W;
    */
    //!cm
    return Res;
}
//J/(kg*K^2)
double DfCp(double T, const double a)
{
    /*
    //Molar mass of water
    const double W = 18 * pow(10, -3);
    const double T_boiling = 373;
    double AA, B, C, D, EE;
    //Polinom for T = (K) /1000
    double T_forCp = T / 1000.;
    if (T <= T_boiling)
    {
        AA = -203.6060;
        B = 1523.290;
        C = -3196.413;
        D = 2474.455;
        EE = 3.855326;
    }
    else
    {
        AA = 30.09200;
        B = 6.832514;
        C = 6.793435;
        D = -2.534480;
        EE = 0.082139;
    }
    //Get value
    if (a < pow(10, -8))
        return 0;
    else
        return (B + 2. * C * T_forCp + 3. * D * pow(T_forCp, 2.) - 2. * EE * pow(T_forCp, -3.)) / 1000. / W;
    */

    return 0;
}
//Formula setting for Density and it's derivative, kg/m^3
double Density(double T, const double a)
{

    double Res;
    const double rho_water = 980.;
    const double rho_steam = 0.4;
    const double T_boiling = 373;
    if (T <= T_boiling)
        Res = rho_water;
    else
        Res = rho_steam;
    /*
    const double R = 8.314462;
    const double pressure = pow(10, 5);
    const double W = 18 * pow(10, -3);
    if (a < pow(10, -8))
    {

    }

    else {
        if (T <= T_boiling)
            Res = rho_water;
        else
            Res = pressure * W / (R * T);
    }
    */
    //!cm
    return Res * pow(10, -6);
}
//kg/(m^3*K)
double DfDensity(double T, const double a)
{
    /*
    const double T_boiling = 373;
    const double R = 8.314462;
    const double p = pow(10, 5);
    const double W = 18 * pow(10, -3);
    if (a < pow(10, -8))
        return 0;
    else
    {
        if (T <= T_boiling)
            return 0;
        else
            return -p * W / (R * pow(T, 2.));
    }
    */
    return 0;
}
//Formula setting for thermal conductivity and it's derivative, Watt/(m*K)
double Lyambda(double T, const double a, int flag_phase)
{

    double Res;
    const double lambda_water = 0.56;
    //const double lambda_water = 0.1;
    const double lambda_steam = 0.05;
    const double T_boiling = 373.;
    /*
    if (T <= T_boiling && flag_phase == 0)
        Res = lambda_water;
    else
        Res = lambda_steam;
    */
    if (a < pow(10, -8))
    {
        if (T <= T_boiling)
            Res = lambda_water;
        else
            Res = lambda_steam;
    }
    else
    {
        if (T <= T_boiling && flag_phase == 0)
            Res = LinIntWConductivity(T);
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
            double teta = T / T_star;
            Res = lambda_star * (sqrt(teta) / (L0 + L1 / teta + L2 / pow(teta, 2.) + L3 / pow(teta, 3.) + L4 / pow(teta, 4.)));
        }
    }
    //!cm
    return Res * 0.01;
}
//Watt/(m*K^2)
double DfLambda(double T, const double a, int flag_phase)
{
    /*
    const double T_boiling = 373;
    if (a < pow(10, -8))
        return 0;
    else
    {
        if (T <= T_boiling && flag_phase == 0)
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
            double teta = T / T_star;
            double bracket_sum = L0 + L1 / teta + L2 / pow(teta, 2.) + L3 / pow(teta, 3.) + L4 / pow(teta, 4.);
            return lambda_star / (pow(bracket_sum, 2.0) * sqrt(teta)) * (0.5 * bracket_sum + L1 / teta + 2. * L2 / pow(teta, 2.)
                + 3. * L3 / pow(teta, 3.) + 4. * L4 / pow(teta, 4.)) / T_star;
        }
    }
    */
    return 0;
}
//Setting arrays of parametres depending on Tmprtr:
void ArraysParameters(const double a, const int N, realtype* r, realtype* T_next,
    const int s, vector <double>& ACp, vector <double>& ADensity, vector <double>& ALambda, int n, double inter, double p)
{
    int r_size = N;
    //Define Parameters
    ACp.resize(r_size);
    //ADfCp.resize(r.size());
    ADensity.resize(r_size);
    //ADfDensity.resize(r.size());
    ALambda.resize(r_size);
    //ADfLambda.resize(r.size());

    for (int i = 0; i < s; i++) {
        //cout << r[i] << endl;
        ACp[i] = Cp(T_next[i], a);
        //ADfCp[i] = DfCp(T_next[i], a);
        ADensity[i] = Density(T_next[i], a);
        //ADfDensity[i] = DfDensity(T_next[i], a);
        ALambda[i] = Lyambda(T_next[i], a, 0);
        //ADfLambda[i] = DfLambda(T_next[i], a, 0);
    }
    for (int i = s + 1; i < r_size; i++) {
        //cout << r[i] << endl;
        ACp[i] = Cp(T_next[i], a);
        //ADfCp[i] = DfCp(T_next[i], a);
        ADensity[i] = Density(T_next[i], a);
        //ADfDensity[i] = DfDensity(T_next[i], a);
        ALambda[i] = Lyambda(T_next[i], a, 1);
        //ADfLambda[i] = DfLambda(T_next[i], a, 1);
    }

    //Plot graphics of Parameters as a function of Temperature
    ofstream Parameters;
    Parameters.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/Data_new/Parameters/Parameters_" + to_string(n) + ".dat");
    Parameters << "TITLE=\"" << "Graphics" << "\"" << endl;
    Parameters << R"(VARIABLES= "rj, m", "T, K", "Cp, J/kg*K", "rho, kg/m^3", "Lambda, Watt/m*K")" << "\n";
    int j = 0;
    for (j; j < s; j++)
        Parameters << r[j] << " " << T_next[j] << " " << ACp[j] << " " << ADensity[j] << " " << ALambda[j] << "\n";
    //Interface
    //Left side
    Parameters << r[s] << " " << T_next[s] << " " << Cp(T_next[s], a) << " " << Density(T_next[s], a) << " "
        << Lyambda(T_next[s], a, 0) << "\n";
    //Right side
    Parameters << r[s] << " " << T_next[s] << " " << Cp(T_next[s], a) << " " << Density(T_next[s], a) << " "
        << Lyambda(T_next[s], a, 1) << "\n";
    for (j = s + 1; j < r_size; j++)
        Parameters << r[j] << " " << T_next[j] << " " << ACp[j] << " " << ADensity[j] << " " << ALambda[j] << "\n";
    Parameters.close();
}


void GraphicsSystEqu(int n, double* r, double* T_next, vector<double>& ALambda, vector<double>& ADensity,
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
        for (i; i < N + 1; i++)
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
double FAdiabatic(double T_center, double T_right, realtype* r, vector<double>& ALambda, vector<double>& ACp,
    vector<double>& ADensity)
{
    //
    //Nevyazka on the left border
    int l_f_num = 0;
    double rj_diff = r[1] - r[0];
    double rj_avg = pow((r[l_f_num + 1] + r[l_f_num]) / 2., 2.);
    double Lambda_j_half = 0.5 * (ALambda[l_f_num + 1] + ALambda[l_f_num]);
    double F = 1. / rj_diff * (rj_avg * Lambda_j_half
        * (T_right - T_center) / (r[l_f_num + 1] - r[l_f_num]));

    return F / (ACp[l_f_num] * ADensity[l_f_num] * pow(r[l_f_num + 1], 2.) / 24.);
    //cout << "f[1]" << f[1] << "\n";
}

double FDrop(double T_left, double T_center, double T_right, realtype* r, vector<double>& ALambda, vector<double>& ACp,
    vector<double>& ADensity, int j)
{
    //
    double rj_diff = r[j + 1] - r[j - 1];
    double rj_avg = (r[j + 1] + r[j]) / 2.;
    double rj_minus_avg = (r[j] + r[j - 1]) / 2.;
    double Lambda_j_half = 0.5 * (ALambda[j + 1] + ALambda[j]);
    double Lambda_j_minus_half = 0.5 * (ALambda[j] + ALambda[j - 1]);
    double F = 2. / rj_diff * (pow(rj_avg, 2.) * Lambda_j_half * (T_right - T_center) / (r[j + 1] - r[j])
        - pow(rj_minus_avg, 2.) * Lambda_j_minus_half * (T_center - T_left) / (r[j] - r[j - 1]));

    //double F = (2. / rj_diff) * (Lambda_j_half * (T_right - T_center) / (r[j + 1] - r[j])
      //  - Lambda_j_minus_half * (T_center - T_left) / (r[j] - r[j - 1]));
    return F / (ACp[j] * ADensity[j] * pow(r[j], 2.));
    //return F / (ACp[j] * ADensity[j]);
}

double FMinus(double T_bef_left, double T_left, double T_center, double T_right, realtype* r, vector<double>& ALambda,
    vector<double>& ACp, vector<double>& ADensity, double a, int s, double p, double& qs_g, double& qs_d)
{
    //
    double F;
    double h = r[s + 1] - r[s - 1];
    //To the left of the interface
    double rs_avg = (r[s - 1] + r[s - 2]) / 2.0;
    double rs_diff = r[s] - rs_avg;
    double Lambda_s_half = 0.5 * (ALambda[s - 1] + ALambda[s - 2]);
    double grad_i_left = (T_right - T_center) / (r[s] - r[s - 1]);
        //(((p / (p + 1.0)) * T_bef_left - ((p + 1.0) / p) * T_left + ((2.0 * p + 1.0) / ((p + 1.0) * p)) * T_right)) / h;
    double grad_s_l = (T_center - T_left) / (r[s - 1] - r[s - 2]);

    
    qs_g = grad_i_left * Lyambda(T_right, a, 0);
    qs_d = grad_s_l * Lambda_s_half;
    F = (1.0 / rs_diff) * (pow(r[s], 2.) * qs_g - pow(rs_avg, 2.) * qs_d);
    //F = (1.0 / rs_diff) * (Lyambda(T_right, a, 0) * grad_i_left - Lambda_s_half * grad_s_l);
    
    //cout << "T_right= " << T_right << "\n";
    cout << "grad_i_left= " << grad_i_left << "\n";
    //cout << "q-(s)= " << Lambda_s_half * grad_s_l << "\n";
    //cout << "q+(s)= " << Lyambda(T_right, a, 0) * grad_i_left << "\n";
    //cout << "Lambda_d= " << Lyambda(T_right, a, 0) << "\n";

    return F / (ACp[s - 1] * ADensity[s - 1] * pow(r[s - 1], 2.));
    //return F / (ACp[s - 1] * ADensity[s - 1]);
}

double FInterface(double T_2bef_left, double T_bef_left, double T_center, double T_aft_right, double T_2aft_right,
    realtype* r, vector<double>& ALambda, vector<double>& ACp, vector<double>& ADensity, double a, int s, double p,
    double& q_i_right, double& q_i_left)
{
    //double h = r[j + 1] - r[j];
    //
    double h = r[s + 1] - r[s - 1];
    double grad_i_right = (((2.0 * p - 7.0) / ((p - 3.0) * (p - 4.0))) * T_center
        + ((p - 4.0) / (p - 3.0)) * T_aft_right + ((p - 3.0) / (4.0 - p)) * T_2aft_right) / h;
    double grad_i_left = (((p / (p + 1.0)) * T_2bef_left - ((p + 1.0) / p) * T_bef_left + ((2.0 * p + 1.0) / ((p + 1.0) * p)) * T_center)) / h;
    q_i_right = Lyambda(T_center, a, 1) * grad_i_right;
    q_i_left = Lyambda(T_center, a, 0) * grad_i_left;
    /*
    cout << "h" << h << "\n";
    cout << "p" << p << "\n";
    cout << "Lyambda(T_center, a, 1)" << Lyambda(T_center, a, 1) << "\n";
    cout << "Lyambda(T_center, a, 0)" << Lyambda(T_center, a, 0) << "\n";
    cout << "T_bef_left" << T_bef_left << "\n";
    cout << "T_2bef_left" << T_2bef_left << "\n";
    cout << "T_center" << T_center << "\n";
    cout << "T_aft_right" << T_aft_right << "\n";
    cout << "T_2aft_right" << T_2aft_right << "\n";

    cout << "grad_i_right" << grad_i_right << "\n";
    cout << "grad_i_left" << grad_i_left << "\n";
    cout << "q_i_right" << q_i_right << "\n";
    cout << "q_i_left" << q_i_left << "\n";
    cout << "Ti" << T_center << "\n";
    */
    //return (q_i_right - q_i_left);
    return pow(r[s], 2.) * (q_i_right - q_i_left);
    //+ L_d * dM;
}

double FPlus(double T_left, double T_center, double T_right, double T_aft_right, realtype* r, vector<double>& ALambda,
    vector<double>& ACp, vector<double>& ADensity, double a, double dM, int s, double p)
{
    double F;
    double term1, term2;
    double h = r[s + 1] - r[s - 1];
    //To the right of the interface
    double rs_avg_plus = (r[s + 1] + r[s + 2]) / 2.0;
    double rs_diff_plus = rs_avg_plus - r[s];
    double Lambda_s_plus_half = 0.5 * (ALambda[s + 2] + ALambda[s + 1]);
    double grad_s_plus_r = (T_right - T_center) / (r[s + 2] - r[s + 1]);
    double grad_i_right = (((2.0 * p - 7.0) / ((p - 3.0) * (p - 4.0))) * T_left
        + ((p - 4.0) / (p - 3.0)) * T_right + ((p - 3.0) / (4.0 - p)) * T_aft_right) / h;

    //term1 = (1.0 / rs_diff_plus) * (Lambda_s_plus_half * grad_s_plus_r - Lyambda(T_left, a, 1) * grad_i_right);
    term1 = (1.0 / rs_diff_plus) * (pow(rs_avg_plus, 2.) * Lambda_s_plus_half * grad_s_plus_r
        - pow(r[s], 2.) * Lyambda(T_left, a, 1) * grad_i_right);
    //term2 = -dM * ((T_center - T_left) / (r[s + 1] - r[s]));
    F = term1;
    //+ term2;

    //cout << "Lambda_g= " << Lyambda(T_left, a, 1) << "\n";
    //cout << "T_left= " << T_left << "\n";
    //cout << "grad_s_plus_r= " << grad_s_plus_r << "\n";
    //cout << "q-(s+1)= " << Lyambda(T_left, a, 1) * grad_i_right << "\n";
    //cout << "q+(s+1)= " << Lambda_s_plus_half * grad_s_plus_r << "\n";
    //cout << "grad_i_right= " << grad_i_right << "\n";
    //cout << "grad_s_plus_r= " << grad_s_plus_r << "\n";
    //return F / (ACp[s + 1] * ADensity[s + 1]);
    return F / (ACp[s + 1] * ADensity[s + 1] * pow(r[s + 1], 2.));
}

double FSteam(double T_left, double T_center, double T_right, realtype* r, vector<double>& ALambda, vector<double>& ACp,
    vector<double>& ADensity, double dM, int j)
{
    double rj_diff = r[j + 1] - r[j - 1];
    double rj_avg = (r[j + 1] + r[j]) / 2.;
    double rj_minus_avg = (r[j] + r[j - 1]) / 2.;
    double Lambda_j_half = 0.5 * (ALambda[j + 1] + ALambda[j]);
    double Lambda_j_minus_half = 0.5 * (ALambda[j] + ALambda[j - 1]);

    //double F = 2. / rj_diff * (Lambda_j_half * (T_right - T_center) / (r[j + 1] - r[j])
       // - Lambda_j_minus_half * (T_center - T_left) / (r[j] - r[j - 1]));
        //- ACp[j] * dM * (T_center - T_left) / (r[j] - r[j - 1]);
    double F = (2. / rj_diff) * (pow(rj_avg, 2.) * Lambda_j_half * (T_right - T_center) / (r[j + 1] - r[j])
        - pow(rj_minus_avg, 2.) * Lambda_j_minus_half * (T_center - T_left) / (r[j] - r[j - 1]));
    // - ACp[j] * dM * (T_center - T_left) / (r[j] - r[j - 1]);

 //return F / (ACp[j] * ADensity[j]);
    return F / (ACp[j] * ADensity[j] * pow(r[j], 2.));

}
/*
void FRight(vector<double>& f, vector<double> T_cur, vector<double> T_next, vector<double> r, vector <double> ACp, vector <double> ADensity, vector <double> ALambda,
    const int N, const double a, const double b, const double c, const double d, const double dt, double dM, int s, double inter,
    const double Ti, double p, const double L_d, double& mod_nevyaz, int flag_evap)
{
    double F;

    //Nevyazka on the left border
    int l_f_num = 0;
    F = 1. / (r[1] - r[0]) * pow((r[l_f_num + 1] + r[l_f_num]) / 2., 2.) * 0.5 * (ALambda[l_f_num + 1] + ALambda[l_f_num])
        * (T_next[l_f_num + 1] - T_next[l_f_num]) / (r[l_f_num + 1] - r[l_f_num])
        - ACp[l_f_num] * ADensity[l_f_num] * pow(r[l_f_num + 1], 2.) / 24. * (T_next[l_f_num] - T_cur[l_f_num]) / dt;
    f[l_f_num] = -F;
    //cout << "f[1]" << f[1] << "\n";

    //Nevyazka in Water
    for (int j = 1; j < s + 1; j++)
    {
        F = FDrop(T_next[j - 1], T_next[j], T_next[j + 1], r, ALambda, j)
            - ACp[j] * ADensity[j] * pow(r[j], 2.) * (T_next[j] - T_cur[j]) / dt;
        f[j] = -F;
        //cout << j << " f1_" << j << " " << f[j] << "\n";
    }
    //Nevyazka in Gaze
    for (int j = s + 1; j < N; j++)
    {
        //cout << "(T_next[j] - T[j][n])" << (T_next[j] - T[j][n]) << endl;
        F = FSteam(T_next[j - 1], T_next[j], T_next[j + 1], r, ALambda, ACp, dM, j)
            - ACp[j] * ADensity[j] * pow(r[j], 2.) * (T_next[j] - T_cur[j]) / dt;
        f[j] = -F;
        //cout << j << " f1_" << j << " " << f[j] << "\n";
    }
    cout << "flag_evap= " << flag_evap << "\n";
    cout << "inter= " << inter << "\n";
    //Near interface Nevyazka
    double term1, term3;
    double h = r[s + 1] - r[s];
    cout << "Ti" << Ti << "p" << p << "\n";

    //To the left of the interface
    double rs_avg = (r[s] + r[s - 1]) / 2.0;
    double rs_diff = inter - rs_avg;
    double q_right = pow(inter, 2.) * Lyambda(Ti, a, 0) * (Ti - T_next[s]) / (inter - r[s]);
    //pow(inter, 2.) * Lyambda(Ti, a, 0) * ((p / ((p + 1) * h)) * T_next[s - 2] - ((p + 1) / (p * h)) * T_next[s - 1]
    //+ ((2 * p + 1) / ((p + 1) * p * h)) * Ti);
    double q_left = pow(rs_avg, 2.) * 0.5 * (ALambda[s] + ALambda[s - 1]) * (T_next[s] - T_next[s - 1]) / (r[s] - r[s - 1]);
    term1 = (1.0 / rs_diff) * (q_right - q_left);
    term3 = -ACp[s] * ADensity[s] * pow(r[s], 2.) * (T_next[s] - T_cur[s]) / dt;
    //cout << "term1= " << term1 << "\n";
    //cout << "term3= " << term3 << "\n";
    cout << "q-(s)= " << q_left << "\n";
    cout << "Lambda_d= " << Lyambda(Ti, a, 0) << "\n";
    cout << "q+(s)= " << q_right << "\n";
    F = term1 + term3;
    f[s] = -F;

    //On the interface
    q_right = pow(inter, 2.) * Lyambda(Ti, a, 1) * ((2.0 * p - 7.0) / ((p - 3.0) * (p - 4.0) * h) * Ti + (p - 4.0) / ((p - 3.0) * h) * T_next[s + 2]
        + (p - 3.0) / ((4.0 - p) * h) * T_next[s + 3]);
    q_left = pow(inter, 2.) * Lyambda(Ti, a, 0) * (p / ((p + 1.0) * h) * T_next[s - 2] - (p + 1.0) / (p * h) * T_next[s - 1]
        + (2.0 * p + 1.0) / ((p + 1.0) * p * h) * Ti);
    F = q_right - q_left + L_d * dM;
    //cout << "term1= " << term1 << "\n";
    //cout << "term2= " << term2 << "\n";
    //cout << "L_d * dM " << L_d * dM << "\n";
    f[N] = -F;

    //To the right of the interface
    double rs_avg_plus = (r[s + 1] + r[s + 2]) / 2.0;
    double rs_diff_plus = rs_avg_plus - inter;
    q_right = pow(rs_avg_plus, 2.) * 0.5 * (ALambda[s + 2] + ALambda[s + 1]) * (T_next[s + 2] - T_next[s + 1]) / (r[s + 2] - r[s + 1]);
    q_left = pow(inter, 2.) * Lyambda(Ti, a, 1) * (T_next[s + 1] - Ti) / (r[s + 1] - inter);
    //pow(inter, 2.) * Lyambda(Ti, a, 1) * ((2 * p - 7) / ((p - 3) * (p - 4) * h) * Ti + (p - 4) / ((p - 3) * h) * T_next[s + 2]
    //+ (p - 3) / ((4 - p) * h) * T_next[s + 3]);
    term1 = (1.0 / rs_diff_plus) * (q_right - q_left);
    term3 = -ACp[s + 1] * (ADensity[s + 1] * pow(r[s + 1], 2.) * (T_next[s + 1] - T_cur[s + 1]) / dt + dM * ((T_next[s + 1] - Ti) / (r[s + 1] - inter)));
    //cout << "term1= " << term1 << "\n";
    //cout << "term3= " << term3 << "\n";
    cout << "Lambda_g= " << Lyambda(Ti, a, 1) << "\n";
    cout << "q-(s+1)= " << q_left << "\n";
    cout << "q+(s+1)= " << q_right << "\n";
    F = term1 + term3;
    f[s + 1] = -F;
    //Calculate modul of Nevyazka
    for (int i = 1; i <= N; i++)
        mod_nevyaz += pow(f[i], 2.0);
    mod_nevyaz = pow(mod_nevyaz, 0.5);

}
*/
/*
void FRightTestCopyForArrays(vector <double>& f, realtype* T_vect, realtype* r, vector <double>& ACp, vector <double>& ADensity,
    vector <double>& ALambda, const int n, const int N, const double a, double dM, int s, double inter,
    const double Ti, double p, const double L_d, double& q_i_right, double& q_i_left)
{
    int N_minus = N - 1;
    double F;
    //Nevyazka on the left border
    int l_f_num = 0;
    double rj_diff = r[1] - r[0];
    double rj_avg = pow((r[l_f_num + 1] + r[l_f_num]) / 2., 2.);
    double Lambda_j_half = 0.5 * (ALambda[l_f_num + 1] + ALambda[l_f_num]);
    F = 1. / rj_diff * (rj_avg * Lambda_j_half
        * (T_vect[l_f_num + 1] - T_vect[l_f_num]) / (r[l_f_num + 1] - r[l_f_num]));
    f[l_f_num] = F / (ACp[l_f_num] * ADensity[l_f_num] * pow(r[l_f_num + 1], 2.) / 24.);
    //cout << "f[1]" << f[1] << "\n";
    //Nevyazka in Water
    for (int j = 1; j < s + 1; j++)
    {
        F = FDrop4Arrays(T_vect[j - 1], T_vect[j], T_vect[j + 1], r, ALambda, j);
        f[j] = F / (ACp[j] * ADensity[j] * pow(r[j], 2.));
        //cout << j << " f1_" << j << " " << f[j] << "\n";
    }
    //Nevyazka in Gaze
    for (int j = s + 1; j < N_minus; j++)
    {
        //cout << "(T_vect[j] - T[j][n])" << (T_vect[j] - T[j][n]) << endl;
        F = FSteam4Arrays(T_vect[j - 1], T_vect[j], T_vect[j + 1], r, ALambda, ACp, dM, j);
        f[j] = F / (ACp[j] * ADensity[j] * pow(r[j], 2.));
        //cout << j << " f1_" << j << " " << f[j] << "\n";
    }

    //cout << "flag_evap= " << flag_evap << "\n";
    //cout << "inter= " << inter << "\n";
    //Near interface Nevyazka
    double term1, term2, term3;
    double h = r[s + 1] - r[s];
    cout << "Ti" << Ti << "p" << p << "\n";
    //////////////////////////////////////////////////////////////////////////////////////////////////
    //Запишем градиенты потоков слева и справа
    double grad_i_right = (2.0 * p - 7.0) / ((p - 3.0) * (p - 4.0) * h) * Ti
        + (p - 4.0) / ((p - 3.0) * h) * T_vect[s + 2] + (p - 3.0) / ((4.0 - p) * h) * T_vect[s + 3];
    double grad_i_left = p / ((p + 1.0) * h) * T_vect[s - 2]
        - (p + 1.0) / (p * h) * T_vect[s - 1] + (2.0 * p + 1.0) / ((p + 1.0) * p * h) * Ti;
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    //To the left of the interface
    double rs_avg = (r[s] + r[s - 1]) / 2.0;
    double rs_diff = inter - rs_avg;
    double Lambda_s_half = 0.5 * (ALambda[s] + ALambda[s - 1]);
    //double grad_s_r = (p / ((p + 1) * h)) * T_vect[s - 2] - ((p + 1) / (p * h)) * T_vect[s - 1]
      //  + ((2 * p + 1) / ((p + 1) * p * h)) * Ti;
        //pow(inter, 2.) * Lyambda(Ti, a, 0) * (Ti - T_vect[s]) / (inter - r[s]);
    double grad_s_l = (T_vect[s] - T_vect[s - 1]) / (r[s] - r[s - 1]);

    term1 = (1.0 / rs_diff) * (pow(inter, 2.) * Lyambda(Ti, a, 0) * grad_i_left - pow(rs_avg, 2.) * Lambda_s_half * grad_s_l);
    term3 = ACp[s] * ADensity[s] * pow(r[s], 2.);
    //cout << "term1= " << term1 << "\n";
    //cout << "term3= " << term3 << "\n";
    cout << "q-(s)= " << Lambda_s_half * grad_s_l << "\n";
    cout << "q+(s)= " << Lyambda(Ti, a, 0) * grad_i_left << "\n";
    cout << "Lambda_d= " << Lyambda(Ti, a, 0) << "\n";

    f[s] = term1 / term3;

    //On the interface
    q_i_right = Lyambda(Ti, a, 1) * grad_i_right;
    q_i_left = Lyambda(Ti, a, 0) * grad_i_left;
    f[N_minus] = pow(inter, 2.) * (q_i_right - q_i_left) + L_d * dM;
    //f[N_minus] = q_i_right - q_i_left;
    //cout << "term1= " << term1 << "\n";
    //cout << "term2= " << term2 << "\n";
    //cout << "L_d * dM " << L_d * dM << "\n";

    //To the right of the interface
    double rs_avg_plus = (r[s + 1] + r[s + 2]) / 2.0;
    double rs_diff_plus = rs_avg_plus - inter;
    double Lambda_s_plus_half = 0.5 * (ALambda[s + 2] + ALambda[s + 1]);
    double grad_s_plus_r = (T_vect[s + 2] - T_vect[s + 1]) / (r[s + 2] - r[s + 1]);
    //double grad_s_plus_l = ((2.0 * p - 7.0) / ((p - 3.0) * (p - 4.0) * h) * Ti
      //  + (p - 4.0) / ((p - 3.0) * h) * T_vect[s + 2] + (p - 3.0) / ((4.0 - p) * h) * T_vect[s + 3]);
    //pow(inter, 2.) * Lyambda(Ti, a, 1) * ((2 * p - 7) / ((p - 3) * (p - 4) * h) * Ti + (p - 4) / ((p - 3) * h) * T_vect[s + 2]
    //+ (p - 3) / ((4 - p) * h) * T_vect[s + 3]);
    term1 = (1.0 / rs_diff_plus) * (pow(rs_avg_plus, 2.) * Lambda_s_plus_half * grad_s_plus_r
        - pow(inter, 2.) * Lyambda(Ti, a, 1) * grad_i_right);
    term2 = -dM * ((T_vect[s + 1] - Ti) / (r[s + 1] - inter));
    term3 = ACp[s + 1] * (ADensity[s + 1] * pow(r[s + 1], 2.));

    //cout << "term1= " << term1 << "\n";
    //cout << "term3= " << term3 << "\n";
    cout << "Lambda_g= " << Lyambda(Ti, a, 1) << "\n";
    cout << "q-(s+1)= " << Lyambda(Ti, a, 1) * grad_i_right << "\n";
    cout << "q+(s+1)= " << Lambda_s_plus_half * grad_s_plus_r << "\n";

    f[s + 1] = term1 / term3;
    //+ term2) / term3;

///////////////////////////////////
//Write Flows on interface
//Plot graphics of Flows as a function of Time
    ofstream Flows;
    Flows.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/Data_new/IntFlows/IntFlows_" + to_string(n) + ".dat");
    Flows << "TITLE=\"" << "Graphics" << "\"" << endl;
    Flows << R"(VARIABLES= "t, s", "q_g, W/m^2", "q_d, W/m^2",)" << "\n";
    int j = 0;
    for (j; j < s + 1; j++)
        Flows << r[j] << " " << q_i_right << " " << q_i_left << "\n";
    Flows.close();
    ///////////////////
}
*/

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

void InitialGrid(int N, const double x_l, const double T_l, const double T_r, vector<double>& r,
    vector<double>& T_cur, double T_cur_i, double h, const double q, const double h_min, const int const_params,
    const int N_inter, const int N_uni, const int N_uni_near_center)
{
    ofstream OutX;
    OutX.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/Data_new/X_grid.dat");
    OutX << "TITLE=\"" << "Graphics" << "\"" << endl;
    OutX << R"(VARIABLES= "j", "hj, m" )" << endl;
    /*
    //Setting the value on the boundaries and initial distribution
    r[0] = x_l;
    T_cur[0] = T_l;
    cout << "h" << h << " r" << 0 << " " << r[0] << " T= " << T_cur[0] << "\n";
    h = h_min;
    //cout << "r1" << " " << r[1] << "\n";
    double T_step = (T_r - T_l) / (N - 2);
    int j = 1;
    for (j; j < r.size(); j++)
    {
        //Define coordinates
        if (j == N_inter || j == N_inter + 1) {
            r[j] = r[j - 1] + h / 2.;
            T_cur[j] = T_cur[j - 1] + T_step / 2.;
        }
        else {
            r[j] = r[j - 1] + h;
            T_cur[j] = T_cur[j - 1] + T_step;
        }
        //cout << "T" << j << " " << T_cur[j] << endl;
        OutX << j << " " << h << " " << "\n";
        cout << "h" << h << " r" << j << " " << r[j] << " T= " << T_cur[j] << "\n";
    }
    OutX.close();
    */
    //Uniform for droplet and non-uniform for gas 1-dimensional grid 
    //Setting the value on the boundaries and initial distribution
    r[0] = x_l;
    T_cur[0] = T_l;
    cout << "r0" << " " << r[0] << "\n";
    h = h_min;
    double T_step = (T_r - T_l) / (N - 2 - N_inter + 0.5);
    //r[1] = r[0] + h / 2.;
    //T_cur[1] = T_l;
    //cout << "r1" << " " << r[1] << "\n";
    int j = 1;
    for (j; j < r.size(); j++)
    {
        //Define coordinates
        if (j == N_inter || j == N_inter + 1)
            r[j] = r[j - 1] + h / 2.;
        else
            r[j] = r[j - 1] + h;
        //Define Temperature
        if (j < N_inter)
            T_cur[j] = T_l;
        else if (j == N_inter)
            T_cur[j] = T_cur_i;
        //else if (j == N_inter + 1)
           // T_cur[j] = T_cur[j - 1] + T_step / 2.;
        //else if (j <= N_uni + 1)
            //T_cur[j] = T_cur[j - 1] + T_step;
        else
            T_cur[j] = T_r;
            //T_cur[j - 1];
        //cout << "T" << j << " " << T_cur[j] << endl;
        OutX << j << " " << h << " " << T_cur[j] << "\n";
        cout << "h" << h << " r" << j << " " << r[j] << " T= " << T_cur[j] << "\n";
    }
    /*
    //Non-uniform grid
    for (j; j < r.size(); j++)
    {
        r[j] = r[j - 1] + h;
        //h = h * q;
        T_cur[j] = T_r;
        //cout << "T" << j << " " << T_cur[j] << endl;
        OutX << j << " " << h << " " << T_cur[j] << "\n";
        cout << "h" << h << " r" << j << " " << r[j] << " T= " << T_cur[j] << "\n";
    }
    */
    OutX.close();
}

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

int Integrate_IDA(int N_x, vector<double>& x_vect, vector<double>& T_vect, vector<double>& Y_vect,
    double& M, int s, double tout1, int call, int number, int print_value, int cons_flag)
{
    void* mem;
    N_Vector yy, yp, avtol, cons;
    realtype rtol, abstol, * yval, * ypval, * atval;
    realtype t0, tout, tret;
    int iout, retval, retvalr;
    SUNMatrix A;
    SUNLinearSolver LS;
    SUNNonlinearSolver NLS;
    SUNContext ctx;
    UserData data;
    data = (UserData)malloc(sizeof * data);

    //int NEQ = NEQ_T + NEQ_Y;
    //MAXIM added
    int NEQ = N_x - 1;
    int k = 0;
    //Дальше не трогаем

    data->kp = k;
    data->outNevyaz = new ofstream;
    data->Nx = N_x;
    data->x = new realtype[N_x];
    data->T = new realtype[N_x];
    data->NEQ = NEQ;
    data->Tl = T_vect[0];
    data->Tr = T_vect[N_x - 1];
    //data->Tp_inter = T_vect[N_inter];
    data->Np_inter = s;
    data->M = M;

    int j = 0;
    for (int i = 0; i < N_x; i++) {
        data->x[i] = x_vect[i];
        data->T[i] = T_vect[i];
        //cout << i << " = " << Ith(res_vect, i + 1) << endl;
    }
    //Put Ti in T-massive as T[N_x] in struct data
    //data->T[N_x] = T_vect[N_x];
    cout << "check2" << "\n";

    mem = NULL;
    cons = yy = yp = avtol = NULL;
    yval = ypval = atval = NULL;
    A = NULL;
    LS = NULL;
    NLS = NULL;
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
    cons = N_VClone(yy);
    if (check_retval((void*)cons, "N_VNew_Serial", 0)) return(1);

    /* Create and initialize  y, y', and absolute tolerance vectors. */
    yval = N_VGetArrayPointer(yy);
    ypval = N_VGetArrayPointer(yp);
    rtol = RCONST(1.0e-3);
    abstol = RCONST(0.1);
    atval = N_VGetArrayPointer(avtol);
    /* Integration limits */
    t0 = ZERO;

    //заполняю идовский массив
    //MAXIM added
    double p = 1.5;
    data->x[s] = data->x[s - 1] + (p - 1) * (data->x[s + 1] - data->x[s - 1]);
    //double q = 20 * pow(10, 3);                  //W / m^2
    double L_d = 2258.2 * pow(10, 3.);              //J / kg 
    //const determination of dM
    //dM = 4 * PI * pow(inter, 2.0) * q / L_d;     // kg / s
    //dM equals zero because it is only warming time
    double dM = 0.;                                       // kg / s
    double a = 1.;
    //pow(10, -10);

    data->pp = p;
    data->L_dp = L_d;
    data->dMp = dM;
    data->ap = a;
    x_vect[s] = x_vect[s - 1] + (p - 1) * (x_vect[s + 1] - x_vect[s - 1]);
    vector <double> f(N_x);
    vector <double> ACp, ADfCp, ADensity, ADfDensity, ALambda, ADfLambda;
    ArraysParameters(a, N_x, data->x, data->T, s, ACp, ADensity, ALambda, 0, x_vect[s], p);
    //Collect dependence Ti от p
    double Ti_out;
    /*
    p = 1.001;
    while (p <= 2.)
    {
        Ti_out = (-Lyambda(T_vect[s], a, 1) * ((p - 4) / (p - 3) * T_vect[s + 2] + (p - 3) / (4 - p) * T_vect[s + 3])
            + Lyambda(T_vect[s], a, 0) * (p / (p + 1) * T_vect[s - 3] - (p + 1) / p * T_vect[s - 2]))
            / (Lyambda(T_vect[s], a, 1) * (2 * p - 7) / ((p - 3) * (p - 4)) - Lyambda(T_vect[s], a, 0) * (2 * p + 1) / ((p + 1) * p));
        cout << "p= " << p << "Ti_out" << Ti_out << "\n";
        p += 0.01;
    }
    p = 1.001;
    */

    double h = data->x[s + 1] - data->x[s - 1];
    double grad_i_right = (2.0 * p - 7.0) / ((p - 3.0) * (p - 4.0) * h) * data->T[s]
        + (p - 4.0) / ((p - 3.0) * h) * data->T[s + 2] + (p - 3.0) / ((4.0 - p) * h) * data->T[s + 3];
    double grad_i_left = p / ((p + 1.0) * h) * data->T[s - 3]
        - (p + 1.0) / (p * h) * data->T[s - 2] + (2.0 * p + 1.0) / ((p + 1.0) * p * h) * data->T[s];
    double qi_g = Lyambda(data->T[s], a, 1) * grad_i_right;
    double qi_d = Lyambda(data->T[s], a, 0) * grad_i_left;

    double qs_g = 0;
    double qs_d = 0;

    for (int i = 0; i < NEQ; i++) {
        atval[i] = RCONST(1.0e-6);
        Ith(yy, i) = T_vect[i];
    }
    //Ith(yy, NEQ_T) = T_vect[N_x];
    //Nevyazka
    for (j = 0; j < NEQ; j++)
    {
        if (j == 0)
            Ith(yp, j) = FAdiabatic(T_vect[0], T_vect[1], data->x, ALambda, ACp, ADensity);
        else if (j < s - 1)
            Ith(yp, j) = FDrop(T_vect[j - 1], T_vect[j], T_vect[j + 1], data->x, ALambda, ACp, ADensity, j);
        else if (j == s - 1) {
            Ith(yp, j) = FMinus(T_vect[s - 3], T_vect[s - 2], T_vect[s - 1], T_vect[s], data->x,
                ALambda, ACp, ADensity, a, s, p, qs_g, qs_d);
        } else if (j == s) {
            Ith(yp, j) = ZERO;
        }
        else if (j == s + 1) {
            Ith(yp, j) = FPlus(T_vect[s], T_vect[s + 1], T_vect[s + 2], T_vect[s + 3], data->x,
                ALambda, ACp, ADensity, a, dM, s, p);
        }
        else if (j < NEQ)
            Ith(yp, j) = FSteam(T_vect[j - 1], T_vect[j], T_vect[j + 1], data->x, ALambda, ACp, ADensity, dM, j);
    }
    //
    ofstream OutTolerances;
    OutTolerances.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/Data_new/AbsTolerances/AbsTol_"
        + to_string(k) + ".dat");
    OutTolerances << "TITLE=\"" << "Graphics" << "\"" << "\n";
    OutTolerances << R"(VARIABLES= "rj, m", "atval_j" )" << "\n";
    ofstream OutIterNevyazka;
    OutIterNevyazka.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/Data_new/IterNevyazka/IterNevyazka_"
        + to_string(k) + ".dat");
    OutIterNevyazka << "TITLE=\"" << "Graphics" << "\"" << "\n";
    OutIterNevyazka << R"(VARIABLES= "rj, m", "F", "dT/dt" )" << "\n";
    for (j = 0; j < NEQ; j++)
    {
        if (j == 0)
            OutIterNevyazka << data->x[j] << " " << abs(FAdiabatic(T_vect[0], T_vect[1], data->x, ALambda, ACp, ADensity)) << " " << abs(ypval[j]) << "\n";
        else if (j < s - 1) {
            OutIterNevyazka << data->x[j] << " " << abs(FDrop(T_vect[j - 1], T_vect[j], T_vect[j + 1], data->x, ALambda,
                ACp, ADensity, j)) << " " << abs(ypval[j]) << "\n";
        }
        else if (j == s - 1) {
            OutIterNevyazka << data->x[j] << " " << abs(FMinus(T_vect[s - 3], T_vect[s - 2], T_vect[s - 1], T_vect[s], data->x,
                ALambda, ACp, ADensity, a, s, p, qs_g, qs_d)) << " " << abs(ypval[j]) << "\n";
        }
        else if (j == s) {
            OutIterNevyazka << data->x[s] << " " << abs(FInterface(T_vect[s - 3], T_vect[s - 2], T_vect[s], T_vect[s + 2], T_vect[s + 3],
                data->x, ALambda, ACp, ADensity, a, s, p, qi_g, qi_d)) << " " << abs(ypval[j]) << "\n";
            cout << "NevyazInter" << abs(FInterface(T_vect[s - 3], T_vect[s - 2], T_vect[s], T_vect[s + 2], T_vect[s + 3],
                data->x, ALambda, ACp, ADensity, a, s, p, qi_g, qi_d)) << " " << abs(ypval[j]) << "\n";
        }
        else if (j == s + 1) {
            OutIterNevyazka << data->x[j] << " " << abs(FPlus(T_vect[s], T_vect[s + 1], T_vect[s + 2], T_vect[s + 3], data->x,
                ALambda, ACp, ADensity, a, dM, s, p)) << " " << abs(ypval[j]) << "\n";
        }
        else if (j > s + 1) {
            OutIterNevyazka << data->x[j] << " " << FSteam(T_vect[j - 1], T_vect[j], T_vect[j + 1], data->x,
                ALambda, ACp, ADensity, dM, j) << " " << abs(ypval[j]) << "\n";
        }
        OutTolerances << data->x[j] << " " << atval[j] << "\n";
    }
    OutTolerances.close();
    OutIterNevyazka.close();
    k++;
    data->kp = k;
    //Opening file for Nevyazka
    data->outNevyaz->open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/Data_new/Nevyaz.dat");
    // Проверяем, удалось ли открыть файл
    if (!data->outNevyaz->is_open()) {
        cout << "Ошибка открытия файла" << endl;
        return 0;
    }
    *(data->outNevyaz) << "TITLE=\"" << "Graphics" << "\"" << endl;
    *(data->outNevyaz) << R"(VARIABLES= "k", "F", "q_g, W/m^2",  "q_d, W/m^2" "qs_g, W/m^2", "qs_d, W/m^2",)" << endl;
    //Calculate modul of Nevyazka
    double mod_nevyaz = 0;
    for (int i = 0; i < N_x; i++)
        mod_nevyaz += pow(f[i], 2.0);
    mod_nevyaz = pow(mod_nevyaz, 0.5);
    *(data->outNevyaz) << k << " " << mod_nevyaz << " " << abs(qi_g) << " "
        << abs(qi_d) << " " << abs(qs_g) << " " << abs(qs_d) << " " << "\n";
    //Ith(yp, NEQ_T) = ZERO;
    cout << "Ith(yp, NEQ)" << Ith(yp, NEQ) << "\n";

    /* Call IDACreate and IDAInit to initialize IDA memory */
    mem = IDACreate(ctx);
    if (check_retval((void*)mem, "IDACreate", 0)) return(1);

    retval = IDAInit(mem, resrob, t0, yy, yp);
    if (check_retval(&retval, "IDAInit", 1)) return(1);
    /* Call IDASVtolerances to set tolerances */

    //retval = IDASVtolerances(mem, rtol, avtol);
    //if (check_retval(&retval, "IDASVtolerances", 1)) return(1);

    /* Call IDASStolerances to set tolerances */

    retval = IDASStolerances(mem, rtol, abstol);
    if (check_retval(&retval, "IDASStolerances", 1)) return(1);
    /*
    if (cons_flag == 1)
    {
        j = 1;
        for (int i = 1; i < NEQ_T; i++) {
            Ith(cons, i) = 0;

        }
        Ith(cons, NEQ_T) = 0;
        for (int i = NEQ_T + 1; i <= NEQ; i++) {
            Ith(cons, i) = 1;
        }
        retval = IDASetConstraints(mem, cons);

    }
    */
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

    /////////////////////////
    // Write Initial dT/dt, q_g, g_d
    ofstream TimeNevyaz;
    TimeNevyaz.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/TimeNevyaz/TimeNevyaz.dat");
    TimeNevyaz << "TITLE=\"" << "Graphics" << "\"" << endl;
    TimeNevyaz << R"(VARIABLES= "t, s", "T0, K", "Ts, K", "Ti, K", "Ts+1, K", "gradT_d, K/m", "gradT_g, K/m", "q_d, W", "q_g, W", "dM, kg/s",
 "rho, kg/m^3", "u_r, m/s", "p, m/(cell size) ", "inter, m" )" << endl;
    //dT/dt убрали
    TimeNevyaz << 0 << " " << T_vect[0] << " " << T_vect[s - 1] << " " << T_vect[s] << " " << T_vect[s + 1] << " " << 0 << " " << 0 << " "
        << abs(qi_d) * pow(10, 4.) << " " << abs(qi_g) * pow(10, 4.) << " " << data->M << " " << Density(T_vect[s], a) << " " << 0 << " "
        << data->pp << " " << data->x[s] << "\n";
    // Write Initial distribution of Temperatures
    ofstream fout;
    double mod_ypval;
    int i;
    fout.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/1file(Temp)/" +
        to_string(call) + "file" + to_string(0) + ".dat");
    fout << "TITLE=\"" << "Graphics" << "\"" << endl;
    fout << R"(VARIABLES= "rj, m", "T, K", "dT/dt", "Lambda, Watt/m*K", "Cp, J/kg*K", "rho, kg/m^3")" << endl;
    i = 0;
    for (i; i < s; i++) {
        fout << x_vect[i] << "  " << T_vect[i] << "  " << Ith(yp, i) << "  "
            << ALambda[i] << "  " << ACp[i] << "  " << ADensity[i] << "\n";
    }
    //Writing on interface
    //Right side
    fout << x_vect[s] << "  " << T_vect[s] << " " << Ith(yp, s) << "  "
        << Lyambda(T_vect[s], a, 0) << "  " << Cp(T_vect[s], a) << "  " << Density(T_vect[s], a) << "\n";
    //Left side
    fout << x_vect[s] << "  " << T_vect[s] << " " << Ith(yp, s) << "  "
        << Lyambda(T_vect[s], a, 1) << "  " << Cp(T_vect[s], a) << "  " << Density(T_vect[s], a) << "\n";
    for (i = s + 1; i < NEQ; i++) {
        fout << x_vect[i] << "  " << T_vect[i] << "  " << Ith(yp, i) << "  "
            << ALambda[i] << "  " << ACp[i] << "  " << ADensity[i] << "\n";
    }
    //For Right Boundary
    fout << x_vect[i] << "  " << T_vect[i] << "  " << ZERO << "  "
        << ALambda[i] << "  " << ACp[i] << "  " << ADensity[i] << "\n";
    fout.close();
    //Call Inital function for help
    //IDACalcIC(mem, IDA_Y_INIT, tout1);

    /* In loop, call IDASolve, print results, and test for error.
       Break out of loop when NOUT preset output times have been reached. */

    iout = 0; tout = tout1;
    data->n_tout = (tout / tout1);
    double tend = pow(10, 5);
    while (iout < print_value) {
        retval = IDASolve(mem, tout, &tret, yy, yp, IDA_NORMAL);
        ExportToArray(T_vect, M, data, yy, N_x);
        //PrintOutput(mem, tret, yy);
        //Print statistics on esach step on screen
        //printf("\nFinal Statistics on %d step:\n", iout + 1);
        //KINPrintAllStats(mem, stdout, SUN_OUTPUTFORMAT_TABLE);
        //IDAPrintAllStats(mem, stdout, SUN_OUTPUTFORMAT_TABLE);
        //Out dT/dt, q_d, g_g depending on timesteps
        /*
        //Calculate modul of Nevyazka
        mod_nevyaz = 0;
        for (int i = 0; i < N_x; i++)
            mod_nevyaz += pow(Ith(rr, i), 2.0);
        mod_nevyaz = pow(mod_nevyaz, 0.5);
        */
        h = data->x[s + 1] - data->x[s - 1];
        grad_i_right = (2.0 * p - 7.0) / ((p - 3.0) * (p - 4.0) * h) * data->T[s]
            + (p - 4.0) / ((p - 3.0) * h) * data->T[s + 2] + (p - 3.0) / ((4.0 - p) * h) * data->T[s + 3];
        grad_i_left = p / ((p + 1.0) * h) * data->T[s - 3]
            - (p + 1.0) / (p * h) * data->T[s - 2] + (2.0 * p + 1.0) / ((p + 1.0) * p * h) * data->T[s];
        qi_g = Lyambda(data->T[s], a, 1) * grad_i_right;
        qi_d = Lyambda(data->T[s], a, 0) * grad_i_left;
        TimeNevyaz << tout << " " << T_vect[0] << " " << T_vect[s - 1] << " " << T_vect[s] << " " << T_vect[s + 1] << " " << 0 << " " << 0 << " "
            << abs(qi_d) * pow(10, 4.) << " " << abs(qi_g) * pow(10, 4.) << " " << data->M << " " << Density(T_vect[s], a) << " " << 0 << " "
            << data->pp << " " << data->x[s] << "\n";
        //Out Files on timesteps
        if ((iout + 1) % number == 0 || iout < 16)
        {
            cout << "t = " << tout << "\n";
            //cout << "M = " << M << "\n";
            fout.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/1file(Temp)/" +
                to_string(call) + "file" + to_string(int(tout / tout1)) + ".dat");
            fout << "TITLE=\"" << "Graphics" << "\"" << endl;
            fout << R"(VARIABLES= "rj, m", "T, K", "dT/dt", "Lambda, Watt/m*K", "Cp, J/kg*K", "rho, kg/m^3")" << endl;
            /*
            mod_ypval = 0;
            for (int i = 0; i < N_x; i++)
                mod_ypval += pow(Ith(yp, i), 2.0);
            mod_ypval = pow(mod_ypval, 0.5);
            */
            i = 0;
            for (i; i < s; i++) {
                fout << x_vect[i] << "  " << T_vect[i] << "  "  << Ith(yp, i) << "  "
                    << ALambda[i] << "  " << ACp[i] << "  " << ADensity[i] << "\n";
            }
            //Writing on interface
            //Right side
            fout << x_vect[s] << "  " << T_vect[s] << " " << Ith(yp, s) << " "
                << Lyambda(T_vect[s], a, 0) << "  " << Cp(T_vect[s], a) << "  " << Density(T_vect[s], a) << "\n";
            //Left side
            fout << x_vect[s] << "  " << T_vect[s] << " " << Ith(yp, s) << " "
                << Lyambda(T_vect[s], a, 1) << "  " << Cp(T_vect[s], a) << "  " << Density(T_vect[s], a) << "\n";
            for (i = s + 1; i < NEQ; i++) {
                fout << x_vect[i] << "  " << T_vect[i] << "  " << Ith(yp, i) << " "
                    << ALambda[i] << "  " << ACp[i] << "  " << ADensity[i] << "\n";
            }
            //For Right Boundary
            fout << x_vect[i] << "  " << T_vect[i] << "  " << ZERO << "  "
                << ALambda[i] << "  " << ACp[i] << "  " << ADensity[i] << "\n";
            fout.close();
        }
        if (check_retval(&retval, "IDASolve", 1)) {
            return(1);
        }
        iout++;
        tout += tout1;
        data->n_tout = (tout / tout1);
    }

    TimeNevyaz.close();
    /* Print final statistics to the screen */
    printf("\nFinal Statistics:\n");
    retval = IDAPrintAllStats(mem, stdout, SUN_OUTPUTFORMAT_TABLE);



    /* Free memory */
    IDAFree(&mem);
    SUNNonlinSolFree(NLS);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
    N_VDestroy(avtol);
    N_VDestroy(yy);
    N_VDestroy(yp);
    N_VDestroy(cons);
    SUNContext_Free(&ctx);

    return(retval);
}


void MySolver(int& retval, vector<double>& Y_vect, double& M, int& N_center)
{
    //InTime
    cout << "INTEGRATE ALL\n\n\n";
    double tout1 = 0.0001;
    int number = 10;
    int print_value = pow(10, 5);
    cout << "print_value" << print_value << "\n";
    //Mymain
    const int const_params = 1;
    const double T_l = 300;
    const double T_r = 1500;
    //Если хочется отдельно задать Ti
    double T_cur_i = 337.637;
        //323.948 p = 1.1;
        //337.637 для p = 1.5 выдает невязку порядка 10^(-6);
        // p= 1.001 Ti_out= 321.232
    //зададим минимальный возможный размер ячейки
    const double h_min = 0.0025;        //cm                       //0.0005 * 2 / 41 = 0.0000243902439024= 24.3902439024 мкр;  
    //Number of cells
    int N = 102;                          //количество узлов. То есть ячеек N - 1, и еще раздвоили ячейку узлом интерфейса
    const int Nd = 20;
    const int N_uni = 2 * Nd;
    const int N_uni_near_center = 5;
    const double x_l = 0.;
    //const double R = 0.0005 м = 500 мкр
    const double x_r = 0.25;   //cm                   //для простой задачи зададим такую область
    double h = 0.2;
    double x_uni = x_l + N_uni * h_min;
    double q = DefineQ(N - N_uni, h_min, x_r, x_uni);
    vector <double> r;
    vector <double> T_cur;
    T_cur.resize(N);
    r.resize(N);
    cout << "check" << "\n";
    InitialGrid(N, x_l, T_l, T_r, r, T_cur, T_cur_i, h, q, h_min, const_params, Nd, N_uni, N_uni_near_center);
    retval = Integrate_IDA(N, r, T_cur, Y_vect, M, Nd, tout1, 1, number, print_value, 0);
}

int main()
{
    //init_consts(num_gas_species, num_react);
    //int N_x = 250;
    double b = 0.01;
    double M;
    double W, rho, Y_H2, Y_O2;
    int N_center;
    int retval;
    double w_dot;
    vector<double> x_vect;
    vector<double> T_vect;
    vector<double> Y_vect;
    double* my_x;
    //ofstream fout;
    //double tout1 = pow(10, -7);
    //int number = 500;
    //int print_value = 8000;
    double h, h_left, dTdx, maxdTdx = 0;
    double x_start, x_finish;
    int cons_flag = 0;
    {
        /*
        N_center = InitialData(N_x, x_vect, T_vect, Y_vect, M);
        Write_to_file("file_INITIAL", fout, x_vect,
            T_vect, Y_vect, M, N_x, 0);
        Add_elem(T_vect, Y_vect, x_vect, N_x, N_center, b);
        //Add_elem_spec(T_vect, Y_vect, x_vect, N_x, N_center, b);
        //Add_elem_spec(T_vect, Y_vect, x_vect, N_x, N_center, b);
        cout << "N_x = " << N_x << "\n";
        //retval = Integrate_Y(N_x, x_vect, T_vect, Y_vect, M, N_center, 0);
        Write_to_file("prefile1_Y", fout, x_vect,
            T_vect, Y_vect, M, N_x, 0);
        */
        //MAXIM added
        MySolver(retval, Y_vect, M, N_center);
        //Дальше не трогаем
        /*
        retval = Integrate_Y(N_x, x_vect, T_vect, Y_vect, M, N_center, 0);
        retval = Integrate(N_x, x_vect, T_vect, Y_vect, M, N_center);
        Write_to_file("file_part1", fout, x_vect,
            T_vect, Y_vect, M, N_x, 1);


        b = 0.05;
        Add_elem_spec(T_vect, Y_vect, x_vect, N_x, N_center, b);
        retval = Integrate_Y(N_x, x_vect, T_vect, Y_vect, M, N_center, 0);
        retval = Integrate(N_x, x_vect, T_vect, Y_vect, M, N_center);
        Write_to_file("file_part2", fout, x_vect,
            T_vect, Y_vect, M, N_x, 2);


        b = 0.01;
        Add_elem_spec(T_vect, Y_vect, x_vect, N_x, N_center, b);
        //Add_elem_spec(T_vect, Y_vect, x_vect, N_x, N_center, b);
        retval = Integrate_Y(N_x, x_vect, T_vect, Y_vect, M, N_center, 0);
        retval = Integrate(N_x, x_vect, T_vect, Y_vect, M, N_center);
        Write_to_file("file_part3", fout, x_vect,
            T_vect, Y_vect, M, N_x, 3);

        //b = 0.0012;
        //Add_elem_spec(T_vect, Y_vect, x_vect, N_x, N_center, b);
        ////Add_elem_spec(T_vect, Y_vect, x_vect, N_x, N_center, b);
        //retval = Integrate_Y(N_x, x_vect, T_vect, Y_vect, M, N_center, 0);
        //retval = Integrate(N_x, x_vect, T_vect, Y_vect, M, N_center);
        //Write_to_file("file_part3", fout, x_vect,
        //    T_vect, Y_vect, M, N_x, 4);
        */
    }
    return 0;
    //T_find();
}


static int resrob(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data)
{
    realtype* yval, * ypval, * rval;
    UserData data;
    realtype* x_cells, * T_vect, * Y_vect, * Tp_vect, * Yp_vect;
    double M;
    data = (UserData)user_data;
    T_vect = data->T;
    Y_vect = data->Y_H2O;
    x_cells = data->x;

    yval = N_VGetArrayPointer(yy);
    ypval = N_VGetArrayPointer(yp);
    rval = N_VGetArrayPointer(rr);
    //cout << "ypvalres0 = " << yval[0] << "\n";
    int myNx = data->Nx;
    ExportToArray(T_vect, data->M, data, yy, data->Nx);
    double Y_H2, Y_O2;
    double W;
    double rho;

    //MAXIM added
    int n_dt;
    int s = data->Np_inter;
    int k;
    int j;
    double p, L_d, dM, a;
    n_dt = data->n_tout;
    k = data->kp;
    p = data->pp;
    L_d = data->L_dp;
    dM = data->dMp;
    a = data->ap;
    //cout << "s= " << s << "\n";
    //cout << "T[s]= " << T_vect[s] << "\n";
    //cout << k << "ypval[s]= " << ypval[s] << "\n";
    x_cells[s] = x_cells[s - 1] + (p - 1) * (x_cells[s + 1] - x_cells[s - 1]);

    vector <double> f(myNx);
    vector <double> ACp, ADfCp, ADensity, ADfDensity, ALambda, ADfLambda;
    ArraysParameters(a, myNx, x_cells, T_vect, s, ACp, ADensity, ALambda, n_dt, data->x[s], p);
    double qi_g = 0;
    double qi_d = 0;
    double qs_g = 0;
    double qs_d = 0;
    //Nevyazka
    for (j = 0; j < myNx - 1; j++)
    {
        if (j == 0)
            rval[j] = FAdiabatic(T_vect[0], T_vect[1], data->x, ALambda, ACp, ADensity) - ypval[j];
        else if (j < s - 1)
            rval[j] = FDrop(T_vect[j - 1], T_vect[j], T_vect[j + 1], data->x, ALambda, ACp, ADensity, j) - ypval[j];
        else if (j == s - 1) {
            rval[j] = FMinus(T_vect[s - 3], T_vect[s - 2], T_vect[s - 1], T_vect[s], data->x,
                ALambda, ACp, ADensity, a, s, p, qs_g, qs_d) - ypval[j];
        }
        else if (j == s) {
            rval[j] = FInterface(T_vect[s - 3], T_vect[s - 2], T_vect[s], T_vect[s + 2], T_vect[s + 3],
                data->x, ALambda, ACp, ADensity, a, s, p, qi_g, qi_d);
            // / (Cp(T_vect[s], a) * Density(T_vect[s], a) * pow(data->x[s], 2.)) - ypval[j];
        }
        else if (j == s + 1) {
            rval[j] = FPlus(T_vect[s], T_vect[s + 1], T_vect[s + 2], T_vect[s + 3], data->x,
                ALambda, ACp, ADensity, a, dM, s, p) - ypval[j];
        }
        else if (j < myNx - 1) {
            rval[j] = FSteam(T_vect[j - 1], T_vect[j], T_vect[j + 1], data->x, ALambda, ACp, ADensity, dM, j) - ypval[j];
        }
    }

    /////////////////////////////////
    //Opening file for T on iteration
    //Define nevyazkas on each iteration
    if (k % (myNx + 1) == 0 && k > myNx)
    {
        //Calculate modul of Nevyazka
        double mod_nevyaz = 0;
        for (int i = 0; i < myNx; i++)
            mod_nevyaz += pow(rval[i], 2.0);
        mod_nevyaz = pow(mod_nevyaz, 0.5);
        //
        *(data->outNevyaz) << k << " " << mod_nevyaz << " " << abs(qi_g) << " "
            << abs(qi_d) << " " << abs(qs_g) << " " << abs(qs_d) << " " << "\n";
        ofstream OutIterNevyazka;
        OutIterNevyazka.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/Data_new/IterNevyazka/IterNevyazka_" + to_string(k) + ".dat");
        OutIterNevyazka << "TITLE=\"" << "Graphics" << "\"" << "\n";
        OutIterNevyazka << R"(VARIABLES= "rj, m", "F", "dT/dt" )" << "\n";
        int i = 0;
        for (i; i < s; i++)
            OutIterNevyazka << x_cells[i] << " " << abs(rval[i]) << " " << abs(ypval[i]) << "\n";
        //Interface
        OutIterNevyazka << x_cells[s] << " " << abs(rval[s]) << " " << abs(ypval[s]) << "\n";
        i++;
        for (i; i < myNx - 1; i++)
            OutIterNevyazka << x_cells[i] << " " << abs(rval[i]) << " " << abs(ypval[i]) << "\n";
        OutIterNevyazka.close();
        //Define Temperature on each iteration
        ofstream OutIterCurrentTemp;
        OutIterCurrentTemp.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/Data_new/IterTemp/IterTemp_" + to_string(k) + ".dat");
        OutIterCurrentTemp << "TITLE=\"" << "Graphics" << "\"" << "\n";
        OutIterCurrentTemp << R"(VARIABLES= "rj, m", "T, K", "Lambda, W/m*K" )" << "\n";
        i = 0;
        for (i; i < s; i++)
            OutIterCurrentTemp << x_cells[i] << " " << T_vect[i] << " " << ALambda[i] << "\n";
        //Interface
        OutIterCurrentTemp << x_cells[s] << " " << T_vect[s] << " " << ALambda[s] << "\n";
        i++;
        for (i; i < myNx; i++)
            OutIterCurrentTemp << x_cells[i] << " " << T_vect[i] << " " << ALambda[i] << "\n";
        OutIterCurrentTemp.close();
    }

    //Counter of the jacobian challenges k
    k++;
    data->kp = k;
    return(0);
}

static int check_retval(void* retvalvalue, const char* funcname, int opt)
{
    int* errretval;

    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && retvalvalue == NULL) {
        fprintf(stderr,
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
        return(1);
    }

    /* Check if retval < 0 */
    else if (opt == 1) {
        errretval = (int*)retvalvalue;
        if (*errretval < 0) {
            fprintf(stderr,
                "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
                funcname, *errretval);
            return(1);
        }
    }

    /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && retvalvalue == NULL) {
        fprintf(stderr,
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
        return(1);
    }

    return(0);
}
// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
