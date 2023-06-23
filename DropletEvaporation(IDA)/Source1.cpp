
#include <iostream>
#include <fstream>
//#include <C:\Users\user\source\Научная работа\DropletEvaporation(IDA)\consts.cpp>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <cmath>
using namespace std;

#include <ida/ida.h>   
//#include <kinsol/kinsol.h>             /* access to KINSOL func., consts. */
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

#define M_PI (3.141592653589793)
#define M_2PI (2.*M_PI)

#define T_BOILING RCONST(373.)


typedef struct {
    realtype* x;
    realtype* T;
    realtype* Y_H2O;
    int Nx;
    int N_m;
    int NEQ;
    int Np_inter;
    int n_tout;
    int flag_evapp;
    int Nd;
    realtype pp;
    realtype p_oldp;
    realtype L_dp;
    realtype dMp;
    realtype ap;
    realtype Tl;
    realtype Tr;
    realtype Tp_inter;
    realtype tout1p;
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
static int grob(realtype t, N_Vector yy, N_Vector yp,
    realtype* gout, void* user_data);

void ExportToArray(vector<double>& T_vect, double& dM, UserData data, N_Vector yy, int N_x)
{
    int s = data->Np_inter;
    int flag_evap = data->flag_evapp;
    //cout << "data->Tl; = " << data->Tl << "\n";
    //T_vect[0] = data->Tl;
    if (flag_evap == 0)
    {
        for (int i = 0; i < N_x - 1; i++)
        {
            T_vect[i] = Ith(yy, i);
        }
        T_vect[N_x - 1] = data->Tr;
        cout << "check\n";
    }
    if (flag_evap == 1)
    {
        for (int i = 0; i < N_x - 1; i++)
        {
            if (i == s)
                i++;
            T_vect[i] = Ith(yy, i);
        }
        T_vect[s] = data->T[s];
        dM = Ith(yy, s);
        cout << "Ith(yy,s)(in while)= " << Ith(yy, s) << "\n";
        T_vect[N_x - 1] = data->Tr;
    }
}

void ExportToArrayResrob(double* T_vect, double& dM, UserData data, N_Vector yy, int N_x)
{
    int s = data->Np_inter;
    int flag_evap = data->flag_evapp;
    //cout << "data->Tl; = " << data->Tl << "\n";
    //T_vect[0] = data->Tl;
    if (flag_evap == 0)
    {
        for (int i = 0; i < N_x - 1; i++)
        {
            T_vect[i] = Ith(yy, i);
        }
        T_vect[N_x - 1] = data->Tr;
    }
    if (flag_evap == 1)
    {
        for (int i = 0; i < N_x - 1; i++)
        {
            if (i == s)
                i++;
            T_vect[i] = Ith(yy, i);
        }
        T_vect[s] = data->T[s];
        cout << "Ith(yy,s)" << Ith(yy, s) << "\n";
        dM = Ith(yy, s);
        T_vect[N_x - 1] = data->Tr;
    }
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

/*
 * Root function routine. Compute functions g_i(t,y) for i = 0,1.
 */

static int grob(realtype t, N_Vector yy, N_Vector yp, realtype* gout,
    void* user_data)
{
    realtype* yval, Ti;
    UserData data;
    data = (UserData)user_data;

    int s = data->Np_inter;

    yval = N_VGetArrayPointer(yy);
    Ti = yval[s];
    gout[0] = Ti - T_BOILING;

    return(0);
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
//Formula setting for Heat Capacity and it's derivative, J/(kg*K)= cm^2/(s^2*K)
double Cp(double T, const double a)
{
    double Res;
    const double Cp_water = 4180.;
    const double Cp_steam = 2000.;
    const double T_boiling = 373;

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
        if (T <= T_boiling)
            Res = Cp_water;
        else
            Res = Cp_steam;
    }
    else
        Res = (AA + B * T_forCp + C * pow(T_forCp, 2.) + D * pow(T_forCp, 3.) + EE * pow(T_forCp, -2.)) / W;
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
//Formula setting for Density and it's derivative, kg/cm^3
double Density(double T, const double a)
{

    double Res;
    const double rho_water = 980.;
    const double rho_steam = 0.4;
    const double T_boiling = 373;

    const double R = 8.314462;
    const double pressure = pow(10, 5);
    const double W = 18 * pow(10, -3);
    if (a < pow(10, -8))
    {
        if (T <= T_boiling)
            Res = rho_water;
        else
            Res = rho_steam;
    }
    else {
        if (T <= T_boiling)
            Res = rho_water;
        else
            Res = pressure * W / (R * T);
    }
    //!cm
    return Res * pow(10, -6);
}
//kg/(cm^3*K)
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
//Formula setting for thermal conductivity and it's derivative, Watt/(cm*K)
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
//Watt/(cm*K^2)
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
    ACp[s] = Cp(T_next[s], a);
    ADensity[s] = Density(T_next[s], a);
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
    Parameters << R"(VARIABLES= "rj, cm", "T, K", "Cp, J/kg*K", "rho, kg/cm^3", "Lambda, Watt/cm*K")" << "\n";
    int j = 0;
    for (j; j < s; j++)
        Parameters << r[j] << " " << T_next[j] << " " << ACp[j] << " " << ADensity[j] << " " << ALambda[j] << "\n";
    //Interface
    //Left side
    Parameters << r[s] << " " << T_next[s] << " " << ACp[s] << " " << ADensity[s] << " "
        << Lyambda(T_next[s], a, 0) << "\n";
    //Right side
    Parameters << r[s] << " " << T_next[s] << " " << ACp[s] << " " << ADensity[s] << " "
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
        OutCurrentTemp << R"(VARIABLES= "rj, cm", "T, K", "Lambda, W/cm*K" )" << "\n";
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
    OutFlow << R"(VARIABLES= "rj, cm", "q, W", "rho, kg/cm^3", "u_r, cm/s" )" << "\n";
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
    OutSurfaceFlow << R"(VARIABLES= "t, s", "T0, K", "Ts, K", "Ti, K", "Ts+1, K", "gradT_d, K/cm", "gradT_g, K/cm", "q_d, W", "q_g, W", "dM, kg/s",
 "rho, kg/cm^3", "u_r, cm/s", "p, cm/(cell size) ", "inter, cm" )" << "\n";
    //Collecting Last Step Nevyazka
    OutLastStepNevyazka.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/Data_new/LastStepNevyazka.dat");
    OutLastStepNevyazka << "TITLE=\"" << "Graphics" << "\"" << endl;
    OutLastStepNevyazka << R"(VARIABLES= "t, s", "F")" << endl;
    //Collecting First Step Nevyazka
    OutFirstStepNevyazka.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/Data_new/FirstStepNevyazka.dat");
    OutFirstStepNevyazka << "TITLE=\"" << "Graphics" << "\"" << endl;
    OutFirstStepNevyazka << R"(VARIABLES= "t, s", "F")" << endl;
}

//Define roots of cubic equation
int Cubic(double* x, double a, double b, double c) {
    double q, r, r2, q3;
    q = (a * a - 3. * b) / 9.; r = (a * (2. * a * a - 9. * b) + 27. * c) / 54.;
    r2 = r * r; q3 = q * q * q;
    if (r2 < q3) {
        double t = acos(r / sqrt(q3));
        a /= 3.; q = -2. * sqrt(q);
        x[0] = q * cos(t / 3.) - a;
        x[1] = q * cos((t + M_2PI) / 3.) - a;
        x[2] = q * cos((t - M_2PI) / 3.) - a;
        return(3);
    }
    else {
        double aa, bb;
        if (r <= 0.) r = -r;
        aa = -pow(r + sqrt(r2 - q3), 1. / 3.);
        if (aa != 0.) bb = q / aa;
        else bb = 0.;
        a /= 3.; q = aa + bb; r = aa - bb;
        x[0] = q - a;
        x[1] = (-0.5) * q - a;
        x[2] = (sqrt(3.) * 0.5) * fabs(r);
        if (x[2] == 0.) return(2);
        return(1);
    }
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
    double grad_i_left = (((p / (p + 1.0)) * T_bef_left - ((p + 1.0) / p) * T_left + ((2.0 * p + 1.0) / ((p + 1.0) * p)) * T_right)) / h;
    double grad_s_l = (T_center - T_left) / (r[s - 1] - r[s - 2]);


    qs_g = grad_i_left * Lyambda(T_right, a, 0);
    qs_d = grad_s_l * Lambda_s_half;
    F = (1.0 / rs_diff) * (pow(r[s], 2.) * qs_g - pow(rs_avg, 2.) * qs_d);

    return F / (ACp[s - 1] * ADensity[s - 1] * pow(r[s - 1], 2.));
    //return F / (ACp[s - 1] * ADensity[s - 1]);
}

double FMinusM(double T_bef_left, double T_left, double T_center, double T_right, realtype* r, vector<double>& ALambda,
    vector<double>& ACp, vector<double>& ADensity, double a, double dM, int s, double p_old, double& p, double dt, double& qs_g, double& qs_d)
{
    //
    double F;
    double h = r[s + 1] - r[s - 1];
    //Redefine p
    double M_part = (dM * dt) / (4. * M_PI * Density(T_right, a) * pow(h, 3.));
    cout << "dt= " << dt << "\n";
    cout << "Density(T_right, a)= " << Density(T_right, a) << "\n";
    cout << "pow(h, 3.)= " << pow(h, 3.) << "\n";
    double p_array[3]{ 0 };
    int s_minus = s - 1; // assign a value to s_minus

    double c1 = (-1. / 4.);
    double c2 = (1 - p_old / 4. - s_minus);
    double c3 = (-1 + pow(p_old, 2) / 4. + 2 * s_minus - s_minus * s_minus);
    double c4 = (p_old - pow(p_old, 2) + pow(p_old, 3) / 4 - 2 * p_old * s_minus + p_old * pow(s_minus, 2) + s_minus * pow(p_old, 2)) - M_part;
    // the result is now stored in the variable 'result'
    double a_cubic = c2 / c1;
    double b = c3 / c1;
    double c = c4 / c1;

    Cubic(p_array, a_cubic, b, c);
    //cout << p_array[0] << "\n";
    //cout << p_array[1] << "\n";
    //cout << p_array[2] << "\n";

    if (M_part != 0 && dM < 1) {
        cout << "dM in interface" << dM << "\n";
        cout << "M_part" << M_part << "\n";
        cout << "p_before" << p << "\n";
        if (p_array[0] < 2. && p_array[0] > 1.)
            p = p_array[0];
        else if (p_array[1] < 2. && p_array[1] > 1.)
            p = p_array[1];                       //!!!!!!!!!maybe it wont be [1]
        else if (p_array[2] < 2. && p_array[2] > 1.)
            p = p_array[2];
        cout << "p" << p << "\n";
    }

    ///////////////////////////////////
    //To the left of the interface
    double rs_avg = (r[s - 1] + r[s - 2]) / 2.0;
    double rs_diff = r[s] - rs_avg;
    double Lambda_s_half = 0.5 * (ALambda[s - 1] + ALambda[s - 2]);
    double grad_i_left = (((p / (p + 1.0)) * T_bef_left - ((p + 1.0) / p) * T_left + ((2.0 * p + 1.0) / ((p + 1.0) * p)) * T_right)) / h;
    double grad_s_l = (T_center - T_left) / (r[s - 1] - r[s - 2]);


    qs_g = grad_i_left * Lyambda(T_right, a, 0);
    qs_d = grad_s_l * Lambda_s_half;
    F = (1.0 / rs_diff) * (pow(r[s], 2.) * qs_g - pow(rs_avg, 2.) * qs_d);

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

    return pow(r[s], 2.) * (q_i_right - q_i_left);
}

double FInterfaceM(double T_2bef_left, double T_bef_left, double T_center, double T_aft_right, double T_2aft_right, double dM,
    realtype* r, vector<double>& ALambda, vector<double>& ACp, vector<double>& ADensity, double L_d, double a, int s, double p_old, double& p, double dt,
    double& q_i_right, double& q_i_left)
{
    double h = r[s + 1] - r[s - 1];
    //Redefine p
    double M_part = (dM * dt) / (4. * M_PI * Density(T_center, a) * pow(h, 3.));
    double p_array[3]{ 0 };
    int s_minus = s - 1; // assign a value to s_minus

    double c1 = (-1. / 4.);
    double c2 = (1 - p_old / 4. - s_minus);
    double c3 = (-1 + pow(p_old, 2) / 4. + 2 * s_minus - s_minus * s_minus);
    double c4 = (p_old - pow(p_old, 2) + pow(p_old, 3) / 4 - 2 * p_old * s_minus + p_old * pow(s_minus, 2) + s_minus * pow(p_old, 2)) - M_part;
    // the result is now stored in the variable 'result'
    double a_cubic = c2 / c1;
    double b = c3 / c1;
    double c = c4 / c1;

    Cubic(p_array, a_cubic, b, c);
    //cout << p_array[0] << "\n";
    //cout << p_array[1] << "\n";
    //cout << p_array[2] << "\n";

    if (M_part != 0 && dM < 1) {
        cout << "dM in interface" << dM << "\n";
        cout << "M_part" << M_part << "\n";
        cout << "p_before" << p << "\n";
        if (p_array[0] < 2. && p_array[0] > 1.)
            p = p_array[0];
        else if (p_array[1] < 2. && p_array[1] > 1.)
            p = p_array[1];                       //!!!!!!!!!maybe it wont be [1]
        else if (p_array[2] < 2. && p_array[2] > 1.)
            p = p_array[2];
        cout << "p" << p << "\n";
    }


    ///////////////////////////////////
    double grad_i_right = (((2.0 * p - 7.0) / ((p - 3.0) * (p - 4.0))) * T_center
        + ((p - 4.0) / (p - 3.0)) * T_aft_right + ((p - 3.0) / (4.0 - p)) * T_2aft_right) / h;
    double grad_i_left = (((p / (p + 1.0)) * T_2bef_left - ((p + 1.0) / p) * T_bef_left + ((2.0 * p + 1.0) / ((p + 1.0) * p)) * T_center)) / h;
    q_i_right = Lyambda(T_center, a, 1) * grad_i_right;
    q_i_left = Lyambda(T_center, a, 0) * grad_i_left;

    //return pow(r[s], 2.) * (q_i_right - q_i_left);
    return pow(r[s], 2.) * (q_i_right - q_i_left) - L_d * dM;
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
    term2 = -ACp[s + 1] * dM * ((T_center - T_left) / (r[s + 1] - r[s]));
    F = term1 + term2;

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

double FPlusM(double T_left, double T_center, double T_right, double T_aft_right, realtype* r, vector<double>& ALambda,
    vector<double>& ACp, vector<double>& ADensity, double a, double dM, int s, double p_old, double& p, double dt)
{
    double F;
    double term1, term2;
    double h = r[s + 1] - r[s - 1];
    //Redefine p
    double M_part = (dM * dt) / (4. * M_PI * Density(T_left, a) * pow(h, 3.));
    double p_array[3]{ 0 };
    int s_minus = s - 1; // assign a value to s_minus

    double c1 = (-1. / 4.);
    double c2 = (1 - p_old / 4. - s_minus);
    double c3 = (-1 + pow(p_old, 2) / 4. + 2 * s_minus - s_minus * s_minus);
    double c4 = (p_old - pow(p_old, 2) + pow(p_old, 3) / 4 - 2 * p_old * s_minus + p_old * pow(s_minus, 2) + s_minus * pow(p_old, 2)) - M_part;
    // the result is now stored in the variable 'result'
    double a_cubic = c2 / c1;
    double b = c3 / c1;
    double c = c4 / c1;

    Cubic(p_array, a_cubic, b, c);
    //cout << p_array[0] << "\n";
    //cout << p_array[1] << "\n";
    //cout << p_array[2] << "\n";

    if (M_part != 0 && dM < 1) {
        cout << "dM in interface" << dM << "\n";
        cout << "M_part" << M_part << "\n";
        cout << "p_before" << p << "\n";
        if (p_array[0] < 2. && p_array[0] > 1.)
            p = p_array[0];
        else if (p_array[1] < 2. && p_array[1] > 1.)
            p = p_array[1];                       //!!!!!!!!!maybe it wont be [1]
        else if (p_array[2] < 2. && p_array[2] > 1.)
            p = p_array[2];
        cout << "p" << p << "\n";
    }


    /////////////////////////////
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
    term2 = -ACp[s + 1] * dM * ((T_center - T_left) / (r[s + 1] - r[s]));
    F = term1 + term2;

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
        - pow(rj_minus_avg, 2.) * Lambda_j_minus_half * (T_center - T_left) / (r[j] - r[j - 1]))
        - ACp[j] * dM * (T_center - T_left) / (r[j] - r[j - 1]);

    //return F / (ACp[j] * ADensity[j]);
    return F / (ACp[j] * ADensity[j] * pow(r[j], 2.));
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

double DefineP(const double N, const double h_min, const double x_r, const double x_l)
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
    OutX << R"(VARIABLES= "j", "hj, cm" )" << endl;
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
int ForReinitialisation(int N_x, vector<double>& x_vect, vector<double>& T_vect, vector<double>& Y_vect,
    double& dM, int& s, double tout1, int call, int number, int print_value, int cons_flag, double& p)
{
    return 0;
}

int Integrate_IDA(int N_x, vector<double>& x_vect, vector<double>& T_vect, vector<double>& Y_vect,
    double& dM, int& s, double tout1, int call, int number, int print_value, int cons_flag, double& p)
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
    int flag_evap = 0;

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
    data->dMp = dM;
    data->flag_evapp = flag_evap;
    data->tout1p = tout1;
    data->Nd = s;

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

    double p_old = p;
    data->x[s] = data->x[s - 1] + (p - 1) * (data->x[s + 1] - data->x[s - 1]);
    //double q = 20 * pow(10, 3);                  //W / cm^2
    double L_d = 2258.2 * pow(10, 3.);              //J / kg 
    double a = 1.;
    //pow(10, -10);

    data->pp = p;
    data->L_dp = L_d;
    data->ap = a;
    x_vect[s] = x_vect[s - 1] + (p - 1) * (x_vect[s + 1] - x_vect[s - 1]);
    vector <double> f(N_x);
    vector <double> ACp, ADfCp, ADensity, ADfDensity, ALambda, ADfLambda;
    ArraysParameters(a, N_x, data->x, data->T, s, ACp, ADensity, ALambda, 0, x_vect[s], p);
    //dM redefinition
    //double dM = (1. / 3.) * ADensity[s] * pow(x_vect[s], 3.);                      // kg / s
    //data->dMp = dM;
    /*
    //Collect dependence Ti от p
    double Ti_out;

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
        }
        else if (j == s) {
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
    OutTolerances << R"(VARIABLES= "rj, cm", "atval_j" )" << "\n";
    ofstream OutIterNevyazka;
    OutIterNevyazka.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/Data_new/IterNevyazka/IterNevyazka_"
        + to_string(k) + ".dat");
    OutIterNevyazka << "TITLE=\"" << "Graphics" << "\"" << "\n";
    OutIterNevyazka << R"(VARIABLES= "rj, cm", "F", "dT/dt" )" << "\n";
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
    *(data->outNevyaz) << R"(VARIABLES= "k", "F", "q_g, W/cm^2",  "q_d, W/cm^2" "qs_g, W/cm^2", "qs_d, W/cm^2",)" << endl;
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

    /* Call IDARootInit to specify the root function grob with 2 components */
    retval = IDARootInit(mem, 1, grob);
    if (check_retval(&retval, "IDARootInit", 1)) return(1);
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
    double u_r;
    ofstream TimeNevyaz;
    TimeNevyaz.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/TimeNevyaz/TimeNevyaz.dat");
    TimeNevyaz << "TITLE=\"" << "Graphics" << "\"" << endl;
    TimeNevyaz << R"(VARIABLES= "t, s", "T0, K", "Ts, K", "Ti, K", "Ts+1, K", "gradT_d, K/cm", "gradT_g, K/cm", "q_d, W", "q_g, W", "dM, kg/s",
 "rho, kg/cm^3", "u_r, cm/s", "p, cm/(cell size) ", "inter, cm" )" << endl;
    //dT/dt убрали
    u_r = dM / (pow(x_vect[s], 2.0) * Density(T_vect[s], a));
    TimeNevyaz << 0 << " " << T_vect[0] << " " << T_vect[s - 1] << " " << T_vect[s] << " " << T_vect[s + 1] << " " << 0 << " " << 0 << " "
        << abs(qi_d) << " " << abs(qi_g) << " " << data->dMp << " " << Density(T_vect[s], a) << " " << u_r << " "
        << data->pp << " " << data->x[s] << "\n";
    // Write Initial distribution of Temperatures
    ofstream fout;
    double mod_ypval;
    int i;
    fout.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/1file(Temp)/" +
        to_string(call) + "file" + to_string(0) + ".dat");
    fout << "TITLE=\"" << "Graphics" << "\"" << endl;
    fout << R"(VARIABLES= "rj, cm", "T, K", "dT/dt", "Lambda, Watt/cm*K", "Cp, J/kg*K", "rho, kg/cm^3")" << endl;
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


    char stats_filename[128];
    FILE* outStats;
    iout = 0; tout = tout1;
    data->n_tout = (tout / tout1);
    double tend = pow(10, 5);
    while (iout < print_value) {
        retval = IDASolve(mem, tout, &tret, yy, yp, IDA_NORMAL);

        //Print statistics on esach step on screen
        printf("\nFinal Statistics on %d step:\n", iout + 1);
        //KINPrintAllStats(mem, stdout, SUN_OUTPUTFORMAT_TABLE);
        //IDAPrintAllStats(mem, stdout, SUN_OUTPUTFORMAT_TABLE);
        /* Print final statistics to a file in CSV format */
        //if (tout)
        sprintf(stats_filename, "C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/Statistics/Stats_% f.csv", tout);  // generate fi
        outStats = fopen(stats_filename, "w");
        if (outStats != NULL) {  // check if file was opened successfully
            // write to the file
            fprintf(outStats, "This is file number %d\n", int(tout / tout1));
            IDAPrintAllStats(mem, outStats, SUN_OUTPUTFORMAT_TABLE);
            fclose(outStats);
        }
        else {
            printf("Error: unable to open file %s\n", stats_filename);
        }

        //break;
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
        //Define u_r
        u_r = dM / (pow(x_vect[s], 2.0) * Density(T_vect[s], a));
        TimeNevyaz << tout << " " << T_vect[0] << " " << T_vect[s - 1] << " " << T_vect[s] << " " << T_vect[s + 1] << " " << 0 << " " << 0 << " "
            << abs(qi_d) << " " << abs(qi_g) << " " << data->dMp << " " << Density(T_vect[s], a) << " " << u_r << " "
            << data->pp << " " << data->x[s] << "\n";
        //Out Files on timesteps
        if ((iout + 1) % number == 0 || iout < 160)
        {
            cout << "t = " << tout << "\n";
            //cout << "dM = " << dM << "\n";
            fout.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/1file(Temp)/" +
                to_string(call) + "file" + to_string(int(tout / tout1)) + ".dat");
            fout << "TITLE=\"" << "Graphics" << "\"" << endl;
            fout << R"(VARIABLES= "rj, cm", "T, K", "dT/dt", "Lambda, Watt/cm*K", "Cp, J/kg*K", "rho, kg/cm^3", "u_r, cm/s" )" << endl;
            /*
            mod_ypval = 0;
            for (int i = 0; i < N_x; i++)
                mod_ypval += pow(Ith(yp, i), 2.0);
            mod_ypval = pow(mod_ypval, 0.5);
            */
            i = 0;
            for (i; i < s; i++) {
                fout << data->x[i] << "  " << T_vect[i] << "  " << Ith(yp, i) << "  "
                    << Lyambda(T_vect[i], a, 0) << "  " << Cp(T_vect[i], a) << "  " << Density(T_vect[i], a) << " " << 0 << "\n";
            }
            //Writing on interface
            //Right side
            fout << data->x[s] << "  " << T_vect[s] << " " << Ith(yp, s) << " " << Lyambda(T_vect[s], a, 0) << "  "
                << Cp(T_vect[s], a) << "  " << Density(T_vect[s], a) << " " << dM / (pow(x_vect[s], 2.0) * Density(T_vect[s], a)) << "\n";
            //Left side
            fout << data->x[s] << "  " << T_vect[s] << " " << Ith(yp, s) << " " << Lyambda(T_vect[s], a, 1) << "  "
                << Cp(T_vect[s], a) << "  " << Density(T_vect[s], a) << " " << dM / (pow(x_vect[s], 2.0) * Density(T_vect[s], a)) << "\n";
            for (i = s + 1; i < NEQ; i++) {
                fout << data->x[i] << "  " << T_vect[i] << "  " << Ith(yp, i) << "  " << Lyambda(T_vect[i], a, 0) << "  "
                    << Cp(T_vect[i], a) << "  " << Density(T_vect[i], a) << " " << dM / (pow(x_vect[i], 2.0) * ADensity[i]) << "\n";
            }
            //For Right Boundary
            fout << data->x[i] << "  " << T_vect[i] << "  " << ZERO << "  " << Lyambda(T_vect[i], a, 0) << "  "
                << Cp(T_vect[i], a) << "  " << Density(T_vect[i], a) << " " << dM / (pow(x_vect[i], 2.0) * ADensity[i]) << "\n";
            fout.close();
        }
        if (check_retval(&retval, "IDASolve", 1)) {
            return(1);
        }
        //Catch when Ti = 373 and reinitialise IDA
        if (retval == IDA_ROOT_RETURN) {
            //Triggering the flag
            Ith(yy, s) = dM;
            flag_evap = 1;
            data->flag_evapp = flag_evap;
            //data->T[s] = T_vect[s];
            cout << "Turn on evap" << "\n";
            cout << "T_vect[s]= " << T_vect[s] << "data->T[s]= " << data->T[s] << "Ith(yy, s)= " << Ith(yy, s) << "\n";
            for (j = 0; j < NEQ; j++)
            {
                if (j == 0)
                    Ith(yp, j) = FAdiabatic(data->T[0], data->T[1], data->x, ALambda, ACp, ADensity);
                else if (j < data->Nd - 1)
                    Ith(yp, j) = FDrop(data->T[j - 1], data->T[j], data->T[j + 1], data->x, ALambda, ACp, ADensity, j);
                else if (j == data->Nd - 1) {
                    Ith(yp, j) = FMinus(data->T[s - 3], data->T[s - 2], data->T[s - 1], data->T[s], data->x,
                        ALambda, ACp, ADensity, a, s, p, qs_g, qs_d);
                }
                else if (j == data->Nd) {
                    Ith(yp, j) = ZERO;
                }
                else if (j == data->Nd + 1) {
                    Ith(yp, j) = FPlus(data->T[s], data->T[s + 1], data->T[s + 2], data->T[s + 3], data->x,
                        ALambda, ACp, ADensity, a, dM, s, p);
                }
                else if (j < NEQ)
                    Ith(yp, j) = FSteam(data->T[j - 1], data->T[j], data->T[j + 1], data->x, ALambda, ACp, ADensity, dM, j);
            }
            retval = IDAReInit(mem, tout, yy, yp);
            if (check_retval(&retval, "IDAReInit", 1)) return(1);
        }
        ExportToArray(T_vect, dM, data, yy, N_x);
        //PrintOutput(mem, tret, yy);
        //Renewing p
        if (data->pp <= 1.001) {
            data->pp = 1.999;
            data->x[s] = data->x[s - 1];
            data->Np_inter -= 1;
            s = data->Np_inter;
            data->T[s] = T_BOILING;

            //Reinitialization
            Ith(yy, s) = data->dMp;
            Ith(yy, s + 1) = data->T[s + 1];
            cout << "Going through a node" << "\n";
            cout << "T_vect[s]= " << T_vect[s] << "data->T[s]= " << data->T[s] << "Ith(yy, s)= " << Ith(yy, s) << "\n";
            for (j = 0; j < NEQ; j++)
            {
                if (j == 0)
                    Ith(yp, j) = FAdiabatic(data->T[0], data->T[1], data->x, ALambda, ACp, ADensity);
                else if (j < s - 1)
                    Ith(yp, j) = FDrop(data->T[j - 1], data->T[j], data->T[j + 1], data->x, ALambda, ACp, ADensity, j);
                else if (j == s - 1) {
                    Ith(yp, j) = FMinus(data->T[s - 3], data->T[s - 2], data->T[s - 1], data->T[s], data->x,
                        ALambda, ACp, ADensity, a, s, p, qs_g, qs_d);
                }
                else if (j == s) {
                    Ith(yp, j) = ZERO;
                }
                else if (j == s + 1) {
                    Ith(yp, j) = FPlus(data->T[s], data->T[s + 1], data->T[s + 2], data->T[s + 3], data->x,
                        ALambda, ACp, ADensity, a, dM, s, p);
                }
                else if (j < NEQ)
                    Ith(yp, j) = FSteam(data->T[j - 1], data->T[j], data->T[j + 1], data->x, ALambda, ACp, ADensity, dM, j);
            }
            retval = IDAReInit(mem, tout, yy, yp);
            if (check_retval(&retval, "IDAReInit", 1)) return(1);

        }
        p_old = data->pp;
        data->p_oldp = p_old;
        if (data->dMp > 300)
            break;
        //
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


void MySolver(int& retval, vector<double>& Y_vect, int& N_center)
{
    //InTime
    cout << "INTEGRATE ALL\n\n\n";
    double tout1 = 0.005;
    int number = 1;
    int print_value = pow(10, 3);
    cout << "print_value" << print_value << "\n";
    //Mymain
    const int const_params = 1;
    double dM = 0;
    double p = 1.1;
    const double T_l = 280;
    const double T_r = 1500;
    //Если хочется отдельно задать Ti
    double T_cur_i = 302.97;
    //316.777 p = 1.5 невязка 6.08293e-05 для T = 280
    //302.97 p = 1.1 невязка 6.0334e-05 для T = 280
    //323.948 p = 1.1; для T = 300
    //337.637 для p = 1.5 выдает невязку порядка 10^(-6) для T = 300;
    // p= 1.001 Ti_out= 321.232 для T = 300
//зададим минимальный возможный размер ячейки
    const double h_min = 0.0025;        //cm                       //0.0005 * 2 / 41 = 0.0000243902439024= 24.3902439024 мкр;  
    //Number of cells
    int N = 102;                          //количество узлов. То есть ячеек N - 1, и еще раздвоили ячейку узлом интерфейса
    int Nd = 20;
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
    retval = Integrate_IDA(N, r, T_cur, Y_vect, dM, Nd, tout1, 1, number, print_value, 0, p);
}

int main()
{
    //init_consts(num_gas_species, num_react);
    //int N_x = 250;
    double b = 0.01;
    double dM;
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
        //MAXIM added
        MySolver(retval, Y_vect, N_center);
        //Дальше не трогаем
    }
    return 0;
}

static int resrob(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data)
{
    realtype* yval, * ypval, * rval;
    UserData data;
    realtype* x_cells, * T_vect, * Y_vect, * Tp_vect, * Yp_vect;
    data = (UserData)user_data;
    T_vect = data->T;
    Y_vect = data->Y_H2O;
    x_cells = data->x;
    int NEQ = data->NEQ;

    yval = N_VGetArrayPointer(yy);
    ypval = N_VGetArrayPointer(yp);
    rval = N_VGetArrayPointer(rr);
    //cout << "ypvalres0 = " << yval[0] << "\n";
    int myNx = data->Nx;
    cout << "data->dMp_bef-export" << data->dMp << "\n";
    ExportToArrayResrob(T_vect, data->dMp, data, yy, data->Nx);
    double Y_H2, Y_O2;
    double W;
    double rho;

    //MAXIM added
    int n_dt;
    int s = data->Np_inter;
    int k;
    int j;
    int m;
    int Nd;
    double p, p_old, L_d, dM, a, flag_evap, tout1;
    n_dt = data->n_tout;
    k = data->kp;
    p = data->pp;
    p_old = data->p_oldp;
    L_d = data->L_dp;
    dM = data->dMp;
    a = data->ap;
    flag_evap = data->flag_evapp;
    tout1 = data->tout1p;
    Nd = data->Nd;

    //Nailing down Ti
    if (flag_evap == 1) {
        data->T[s] = T_BOILING;
        T_vect[s] = data->T[s];
    }
    x_cells[s] = x_cells[s - 1] + (p - 1) * (x_cells[s + 1] - x_cells[s - 1]);
    data->x[s] = x_cells[s];
    //cout << "s= " << s << "\n";
    cout << "T[s]= " << T_vect[s] << "\n";
    cout << "dM" << dM << "\n";

    vector <double> f(myNx);
    vector <double> ACp, ADfCp, ADensity, ADfDensity, ALambda, ADfLambda;
    ArraysParameters(a, myNx, x_cells, T_vect, s, ACp, ADensity, ALambda, n_dt, data->x[s], p);
    double qi_g = 0;
    double qi_d = 0;
    double qs_g = 0;
    double qs_d = 0;
    double Ti;
    cout << "flag_evap = " << flag_evap << "\n";
    //Nevyazka for Warming
    if (flag_evap == 0)
    {
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
    }
    else if (flag_evap == 1)
    {
        //Define number of nodes traversed 
        m = Nd - s;
        //Nevyazka for Evaporation
        cout << "p_at start of evaporation= " << p << "\n";
        cout << "data->dMp in evap cycle= " << data->dMp << "\n";
        //if (s == 20) {
        for (j = 0; j < myNx - 1; j++)
        {
            if (j == 0)
                rval[j] = FAdiabatic(T_vect[0], T_vect[1], data->x, ALambda, ACp, ADensity) - ypval[j];
            else if (j < s - 1)
                rval[j] = FDrop(T_vect[j - 1], T_vect[j], T_vect[j + 1], data->x, ALambda, ACp, ADensity, j) - ypval[j];
            else if (j == s - 1) {
                rval[s - 1] = FMinusM(T_vect[s - 3], T_vect[s - 2], T_vect[s - 1], T_BOILING, data->x,
                    ALambda, ACp, ADensity, a, data->dMp, s, p_old, p, tout1, qs_g, qs_d) - ypval[j];
                cout << "p_minus= " << p << "\n";
            }
            else if (j == s) {
                //Change p in Function for interface using changed dM 
                rval[s] = FInterfaceM(T_vect[s - 3], T_vect[s - 2], T_BOILING, T_vect[s + 2], T_vect[s + 3], data->dMp,
                    data->x, ALambda, ACp, ADensity, L_d, a, s, p_old, p, tout1, qi_g, qi_d);
                cout << "p_interface= " << p << "\n";
                // / (Cp(T_vect[s], a) * Density(T_vect[s], a) * pow(data->x[s], 2.)) - ypval[j];
            }
            else if (j == s + 1) {
                rval[s + 1] = FPlusM(T_BOILING, T_vect[s + 1], T_vect[s + 2], T_vect[s + 3], data->x,
                    ALambda, ACp, ADensity, a, data->dMp, s, p_old, p, tout1) - ypval[j];
                cout << "p_plus= " << p << "\n";
            }
            else if (j < myNx - 1) {
                rval[j] = FSteam(T_vect[j - 1], T_vect[j], T_vect[j + 1], data->x, ALambda, ACp, ADensity, data->dMp, j) - ypval[j];
            }
        }
        //}
            /*
        else {
            for (j = 0; j < myNx - 1; j++)
            {
                if (j == 0)
                    rval[j] = FAdiabatic(T_vect[0], T_vect[1], data->x, ALambda, ACp, ADensity) - ypval[j];
                else if (j < s - 1)
                    rval[j] = FDrop(T_vect[j - 1], T_vect[j], T_vect[j + 1], data->x, ALambda, ACp, ADensity, j) - ypval[j];
                else if (j == s - 1) {
                    rval[s - 1] = FMinusM(T_vect[s - 3], T_vect[s - 2], T_vect[s - 1], T_BOILING, data->x,
                        ALambda, ACp, ADensity, a, data->dMp, s, p_old, p, tout1, qs_g, qs_d) - ypval[j];
                    cout << "p_minus= " << p << "\n";
                }
                else if (j == s) {
                    rval[s + 1] = FPlusM(T_BOILING, T_vect[s + 1], T_vect[s + 2], T_vect[s + 3], data->x,
                        ALambda, ACp, ADensity, a, data->dMp, s, p_old, p, tout1) - ypval[j];
                    cout << "p_plus= " << p << "\n";
                }
                else if (j == Nd) {
                    //Change p in Function for interface using changed dM
                    rval[s] = FInterfaceM(T_vect[s - 3], T_vect[s - 2], T_BOILING, T_vect[s + 2], T_vect[s + 3], data->dMp,
                        data->x, ALambda, ACp, ADensity, L_d, a, s, p_old, p, tout1, qi_g, qi_d);
                    cout << "p_interface= " << p << "\n";
                }
                else if (j < myNx - 1) {
                    rval[j] = FSteam(T_vect[j - 1], T_vect[j], T_vect[j + 1], data->x, ALambda, ACp, ADensity, data->dMp, j) - ypval[j];
                }
            }
            */
            //}
        cout << "rval[s - 1]= " << rval[s - 1] << "\n";
        cout << "rval[s]= " << rval[s] << "\n";
        cout << "rval[s + 1]= " << rval[s + 1] << "\n";
        cout << "dM= " << dM << " data->dMp= " << data->dMp << "\n";
    }

    /////////////////////////////////
    //Opening file for T on iteration
    //Define nevyazkas on each iteration
    if (k % (NEQ + 2) == 0 && k > NEQ)
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
        OutIterNevyazka << R"(VARIABLES= "rj, cm", "F", "dT/dt" )" << "\n";
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
        OutIterCurrentTemp << R"(VARIABLES= "rj, cm", "T, K", "Lambda, W/cm*K" )" << "\n";
        i = 0;
        for (i; i < s; i++)
            OutIterCurrentTemp << x_cells[i] << " " << T_vect[i] << " " << ALambda[i] << "\n";
        //Interface
        OutIterCurrentTemp << x_cells[s] << " " << T_vect[s] << " " << Lyambda(T_vect[s], a, 0) << "\n";
        i++;
        for (i; i < myNx; i++)
            OutIterCurrentTemp << x_cells[i] << " " << T_vect[i] << " " << ALambda[i] << "\n";
        OutIterCurrentTemp.close();
    }

    //Counter of the jacobian challenges k
    k++;
    data->kp = k;
    data->pp = p;
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
