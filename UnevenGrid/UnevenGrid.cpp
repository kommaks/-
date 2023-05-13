// UnevenGrid.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

// VariableGrid.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

// NewtonLambdaConst.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//


#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

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
/* Now setting the initial value by myself
int DefineN(double h, const double q, const double x_l, const double x_r)
{
    if (q == 1)
    {
        return (x_r - x_l) / h;
    }
    if (q > 1)
    {
        return log(1 + (x_r - x_l) * (q - 1) / h) / log(q);
    }
    if (q < 1)
    {
        return log(1 - (x_r - x_l) * (1 - q) / h) / log(q);
    }
}
*/
//Linear interpolation
double LinIntWConductivity(double T)
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
        533.16, 553.16, 573.16, 593.16, 613.16, 633.16};
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
    for (int i = 0; i < ny; i+=shag)
    {
        LAGR_COND_copy[j] = LAGR_COND[i];
        j += 1;
    }
    return Lagrange1(LAGR_T, LAGR_COND_copy, nx, T);
}

//Formula setting for Heat Capacity and it's derivative, J / (kg * K)
double Cp(double T, const double a)
{
    const double Cp_water = 4180.;
    const double Cp_steam = 2000.;
    //Molar mass of water
    const double W = 18 * pow(10, -3);
    const double T_boiling = 373;
    double A, B, C, D, E;
    //Polinom for T = (K) /1000
    double T_forCp = T / 1000.;
    //For water
    if (T <= T_boiling)
    {
        A = -203.6060;
        B = 1523.290;
        C = -3196.413;
        D = 2474.455;
        E = 3.855326;
    }
    //For steam
    else
    {
        A = 30.09200;
        B = 6.832514;
        C = 6.793435;
        D = -2.534480;
        E = 0.082139;
    }
    //Return
    if (a == pow(10, -9))
        return 1;
    else if (a == pow(10, -8))
    {
        //For water
        if (T <= T_boiling)
            return Cp_water;
        //For steam
        else
            return Cp_steam;
    }
    else
        return (A + B * T_forCp + C * pow(T_forCp, 2.) + D * pow(T_forCp, 3.) + E * pow(T_forCp, -2.)) / W;
}
double DfCp(double T, const double a)
{
    //Molar mass of water
    const double W = 18 * pow(10, -3); // [kg / mol]
    const double T_boiling = 373;
    double A, B, C, D, E;
    //Polinom for T = (K) /1000
    double T_forCp = T / 1000.;
    //For water
    if (T <= T_boiling)
    {
        A = -203.6060;
        B = 1523.290;
        C = -3196.413;
        D = 2474.455;
        E = 3.855326;
    }
    //For steam
    else
    {
        A = 30.09200;
        B = 6.832514;
        C = 6.793435;
        D = -2.534480;
        E = 0.082139;
    }
    if (a == pow(10, -9))
        return 0.;
    else if (a == pow(10, -8))
    {
        //For water
        if (T <= T_boiling)
            return 0.;
        //For steam
        else
            return 0.;
    }
    else
        return ((B + 2. * C * T_forCp + 3. * D * pow(T_forCp, 2.) - 2. * E * pow(T_forCp, -3.)) / 1000.) / W;
}
//Formula setting for Density and it's derivative, kg/m^3
double Density(double T, const double a) 
{                                              
    const double rho_water = 980.;
    const double rho_steam = 0.4;
    const double T_boiling = 373;
    const double R = 8.314462;
    const double p = pow(10, 5);
    const double W = 18 * pow(10, -3);
    if (a == pow(10, -9))
        return 1;
    else if (a == pow(10, -8))
    {
        //For water
        if (T <= T_boiling)
            return rho_water;
        //For steam
        else
            return rho_steam;
    }
    else
        return p * W / (R * T);
}
double DfDensity(double T, const double a)
{
    const double T_boiling = 373;
    const double rho_steam = 0.4;
    const double R = 8.314462;
    const double p = pow(10, 5);
    const double W = 18 * pow(10, -3);
    if (a == pow(10, -9))
        return 0.;
    else if (a == pow(10, -8))
    {
        //For water
        if (T <= T_boiling)
            return 0.;
        //For steam
        else
            return 0.;
    }
    else
        return - p * W / (R * pow(T, 2.));
}
//Formula setting for thermal conductivity and it's derivative, Watt / (m*K)
double Lambda(double T, const double a, const double b, const double c)
{
    const double lambda_water = 0.56;
    const double lambda_steam = 0.05;
    const double T_boiling = 373;
    if (a == pow(10, -9))
        return a * pow(T, 2.0) + b * T + c;
    else if (a == pow(10, -8))
    {
        //For water
        if (T <= T_boiling)
            return lambda_water;
        //For steam
        else
            return lambda_steam;
    }
    else
    {
        if (T <= T_boiling)
            LinIntWConductivity(T);
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
            return lambda_star * (sqrt(teta) / (L0 + L1 / teta + L2 / pow(teta, 2.) + L3 / pow(teta, 3.) + L4 / pow(teta, 4.)));
        }
    }
}
double DfLambda(double T, const double a, const double b)
{
    const double T_boiling = 373;
    if (a == pow(10, -9))
        return 2 * a * T + b;
    else if (a == pow(10, -8))
    {
        //For water
        if (T <= T_boiling)
            return 0.;
        //For steam
        else
            return 0.;
    }
    else
    {
        if (T <= T_boiling)
            0.;
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
            double bracket_sum = pow(L0 + L1 / teta + L2 / pow(teta, 2.) + L3 / pow(teta, 3.) + L4 / pow(teta, 4.), 2.);
            return lambda_star * ((1. / (bracket_sum * T_star * sqrt(teta))) * (1. / 2. + (L1 / teta + 2. * L2 / pow(teta, 2.)
                + 3. * L3 / pow(teta, 3.) + 4. * L4 / pow(teta, 4.)) / bracket_sum));
        }
    }
}
//Plot Graphics of parameters
void GraphicsParameters(const double a, const double b, const double c, double T_l, const double T_r)
{
    vector <double> Tg;
    int dTg = 2;
    int Ng = (T_r - T_l) / dTg;
    Tg.resize(Ng + 1);
    Tg[0] = T_l;
    for (int i = 1; i < Ng + 1; i++) {
        Tg[i] = Tg[i - 1] + dTg;
    }
    ofstream Parameters;
    Parameters.open("C:/Users/user/source/Научная работа/UnevenGrid/Data/Parameters.dat");
    Parameters << "TITLE=\"" << "Graphics" << "\"" << endl;
    Parameters << R"(VARIABLES= "T", "Cp", "DfCp", "rho", "Dfrho", "Lambda", "DfLambda")" << endl;
    for (int i = 0; i < Ng + 1; i++) {
        //cout << r[i] << endl;
        Parameters << Tg[i] << " " << Cp(Tg[i], a) << " " << DfCp(Tg[i], a) << " " << Density(Tg[i], a) 
            << " " << DfDensity(Tg[i], a) << " " << Lambda(Tg[i], a, b, c) << " " << DfLambda(Tg[i], a, b) << "\n";
    }
    Parameters.close();
}
//Setting arrays of parametres depending on T:
void ArraysParameters(double T, const double a, const double b, const double c, const int N)
{
    //не закончено
    vector <double> Cp, DfCp, Density, DfDensity;
    Cp.resize(N + 1);
    DfCp.resize(N + 1);
    Density.resize(N + 1);
    DfDensity.resize(N + 1);
    return;
}
//Setting j-1 element of matrix
double DfLeft(vector<double>& r, vector<double>& T_next, const double a, const double b, const double c,
    int j)
{
    return (2. / (r[j + 1] - r[j - 1]) / pow(r[j], 2.)) * pow((r[j] + r[j - 1]) / 2., 2.)
        * (-DfLambda((T_next[j] + T_next[j - 1]) / 2., a, b) * (T_next[j] - T_next[j - 1]) / (r[j] - r[j - 1])
            + Lambda((T_next[j] + T_next[j - 1]) / 2., a, b, c) / (r[j] - r[j - 1]));
}
//Setting j element of matrix
double DfCenter(vector<double>& r, vector<double>& T_next, const double a, const double b, const double c,
    const double d, int j, const double dt, const double dtau)
{
    //if (a == pow(10, -9))
    //{
    //    return (2. / (r[j + 1] - r[j - 1]) / pow(r[j], 2.)) * (pow((r[j + 1] + r[j]) / 2., 2.)
    //       * (DfLambda((T_next[j] + T_next[j + 1]) / 2., a, b) * (T_next[j + 1] - T_next[j]) / (r[j + 1] - r[j])
    //            - Lambda((T_next[j] + T_next[j + 1]) / 2., a, b, c) / (r[j + 1] - r[j])) - pow((r[j] + r[j - 1]) / 2., 2.)
    //        * (DfLambda((T_next[j] + T_next[j - 1]) / 2., a, b) * (T_next[j] - T_next[j - 1]) / (r[j] - r[j - 1])
    //            + Lambda((T_next[j] + T_next[j - 1]) / 2., a, b, c) / (r[j] - r[j - 1]))) - 1. / dt;
    //}
    //else
    //{
    return (2. / (r[j + 1] - r[j - 1]) / pow(r[j], 2.)) * (pow((r[j + 1] + r[j]) / 2., 2.) * (DfLambda((T_next[j] + T_next[j + 1]) / 2., a, b) *
        (T_next[j + 1] - T_next[j]) / (r[j + 1] - r[j]) - Lambda((T_next[j] + T_next[j + 1]) / 2., a, b, c) / (r[j + 1] - r[j])) -
        pow((r[j] + r[j - 1]) / 2., 2.) * (DfLambda((T_next[j] + T_next[j - 1]) / 2., a, b) * (T_next[j] - T_next[j - 1]) / (r[j] - r[j - 1])
            + Lambda((T_next[j] + T_next[j - 1]) / 2., a, b, c) / (r[j] - r[j - 1]))) - (Cp(T_next[j], a) * Density(T_next[j], a)
                + (DfCp(T_next[j], a) * Density(T_next[j], a) + Cp(T_next[j], a) * DfDensity(T_next[j], a)) * T_next[j]) / dt;
        //- (Cp(T_next[j], a) * Density(T_next[j], a)) / dtau;
    //}
}
//Setting j+1 element of matrix
double DfRight(vector<double>& r, vector<double>& T_next, const double a, const double b, const double c,
    int j)
{
    return (2. / (r[j + 1] - r[j - 1]) / pow(r[j], 2.)) * pow((r[j + 1] + r[j]) / 2., 2.)
        * (DfLambda((T_next[j + 1] + T_next[j]) / 2., a, b) * (T_next[j + 1] - T_next[j]) / (r[j + 1] - r[j])
            + Lambda((T_next[j + 1] + T_next[j]) / 2., a, b, c) / (r[j + 1] - r[j]));
}

double AdiabaticBoundaryDfs(vector<double>& r, vector<double>& T_next, const double a, const double b, const double c,
    const int j, const double dt, const double dtau)
{
    return (2. / (r[2] - r[0]) / pow(r[1], 2.)) * (pow((r[2] + r[1]) / 2., 2.) * (DfLambda((T_next[1] + T_next[2]) / 2., a, b) *
        (T_next[2] - T_next[1]) / (r[2] - r[1]) - Lambda((T_next[1] + T_next[2]) / 2., a, b, c) / (r[2] - r[1]))) - (Cp(T_next[1], a)
            * Density(T_next[1], a) + (DfCp(T_next[1], a) * Density(T_next[1], a) + Cp(T_next[1], a) * DfDensity(T_next[1], a)) * T_next[1]) / dt
        - (Cp(T_next[1], a) * Density(T_next[1], a)) / dtau;
}
void Jacobian(vector<vector<double>>& J, vector<double>& r, vector<double>& T_next, const int N_minus, const double a,
    const double b, const double c, const double d, const double dt, const int const_params, const double dtau)
{
    if (const_params == 0)
    {
            J[1][1] = DfCenter(r, T_next, a, b, c, d, 1, dt, dtau);
            J[1][2] = DfRight(r, T_next, a, b, c, 1);
            //cout << "J[1][2]" << J[1][2] << endl;
            for (int j = 2; j < N_minus; j++) {
                J[j][j - 1] = DfLeft(r,T_next, a, b, c, j);
                J[j][j] = DfCenter(r, T_next, a, b, c, d, j, dt, dtau);
                J[j][j + 1] = DfRight(r, T_next, a, b, c, j);
            }
            J[N_minus][N_minus - 1] = DfLeft(r, T_next, a, b, c, N_minus);
            J[N_minus][N_minus] = DfCenter(r, T_next, a, b, c, d, N_minus, dt, dtau);
    }
    else if (const_params == 2 || const_params == 1)
    {
        //Jacobians for left boundary
        J[1][1] = AdiabaticBoundaryDfs(r, T_next, a, b, c, 1, dt, dtau);
        J[1][2] = DfRight(r, T_next, a, b, c, 1);
        for (int j = 2; j < N_minus; j++) {
            J[j][j - 1] = DfLeft(r, T_next, a, b, c, j);
            J[j][j] = DfCenter(r, T_next, a, b, c, d, j, dt, dtau);
            J[j][j + 1] = DfRight(r, T_next, a, b, c, j);
            cout << J[j][j - 1] << " " << J[j][j] << " " << J[j][j + 1] << "\n";
        }
        J[N_minus][N_minus - 1] = DfLeft(r, T_next, a, b, c, N_minus);
        J[N_minus][N_minus] = DfCenter(r, T_next, a, b, c, d, N_minus, dt, dtau);
    }
}

void F(vector<double>& f, vector<double>& T_cur, vector<double>& T_next, const int N_minus, vector<double>& r, const double a,
    const double b, const double c, const double d, const double dt, const int const_params)
{
    //cout << "Dens" << Density(T_next[0], a) << "\n";
    //cout << "Cp" << Cp(T_next[0], a) << "\n";
    //cout << "Lambda" << Lambda(T_next[0], a, b, c) << "\n";
    if (const_params == 0)
    {
        for (int j = 1; j < N_minus + 1; j++)
        {
            //cout << "(T_next[j] - T_cur[j])" << (T_next[j] - T_cur[j] << endl;
            f[j] = -((2. / (r[j + 1] - r[j - 1]) / pow(r[j], 2.)) * (pow((r[j + 1] + r[j]) / 2., 2.) * Lambda((T_next[j + 1] + T_next[j]) / 2., a, b, c)
                * (T_next[j + 1] - T_next[j]) / (r[j + 1] - r[j]) - pow((r[j] + r[j - 1]) / 2., 2.)
                * Lambda((T_next[j] + T_next[j - 1]) / 2., a, b, c) * (T_next[j] - T_next[j - 1])
                / (r[j] - r[j - 1])) - Cp(T_next[j], a) * Density(T_next[j], a) * (T_next[j] - T_cur[j]) / dt);
        }
    }
    else if (const_params == 2 || const_params == 1)
    {
        f[1] = -((2. / (r[2] - r[0]) / pow(r[1], 2.)) * (pow((r[2] + r[1]) / 2., 2.)
            * Lambda((T_next[2] + T_next[1]) / 2., a, b, c)
            * (T_next[2] - T_next[1]) / (r[2] - r[1])) - Cp(T_next[1], a) * Density(T_next[1], a) * (T_next[1] - T_cur[1]) / dt);
        //cout << "f[0]" << f[0] << "\n";
        for (int j = 2; j < N_minus + 1; j++)
        {
            //cout << "(T_next[j] - T[j][n])" << (T_next[j] - T[j][n]) << endl;
            f[j] = -((2. / (r[j + 1] - r[j - 1]) / pow(r[j], 2.)) * (pow((r[j + 1] + r[j]) / 2., 2.) * Lambda((T_next[j + 1] + T_next[j]) / 2., a, b, c)
                * (T_next[j + 1] - T_next[j]) / (r[j + 1] - r[j]) - pow((r[j] + r[j - 1]) / 2., 2.)
                * Lambda((T_next[j] + T_next[j - 1]) / 2., a, b, c) * (T_next[j] - T_next[j - 1]) / (r[j] - r[j - 1]))
                - Cp(T_next[j], a) * Density(T_next[j], a) * (T_next[j] - T_cur[j]) / dt);
            cout << j << " f" << j << " " << f[j] << endl;
            //cout << j << " " << T_next[j] << "\n";
            //cout << j << "Lambda " << Lambda(T_next[j], a, b, c) << "DfLambda " << DfLambda(T_next[j], a, b) << "\n";
        }
    }
}

void Progonka(vector<double>& T_next, vector<double>& alpha, vector<double>& beta, vector<double>& dT,
    vector<double>& nevyaz, vector < vector <double> >& J, const int N_minus, const int Nd)
{
    const double T_boiling = 373;
    //считаем коэффициенты слева
    //Define denominator for constant formules in method
    double zn;
    zn = J[1][1];
    //cout << "zn " << zn << endl;
    alpha[1] = -J[1][2] / zn;
    //cout << "a[1]" << alpha[1] << endl;
    beta[1] = nevyaz[1] / zn;
    //cout << "B[1]" << beta[1] << endl;
    //считаем коэффициенты для всех узлов
    for (int i = 2; i < N_minus; i++) {
        zn = J[i][i] + J[i][i - 1] * alpha[i - 1];
        //cout << "y" << 1 << " " << zn << endl;
        alpha[i] = -J[i][i + 1] / zn;
        //cout << "a" << i << " " << alpha[i] << endl;
        beta[i] = (nevyaz[i] - J[i][i - 1] * beta[i - 1]) / zn;
        //cout << "B" << i << " " << beta[i] << endl;
    }
    dT[N_minus] = (nevyaz[N_minus] - J[N_minus][N_minus - 1] * beta[N_minus - 1]) /
        (J[N_minus][N_minus] + J[N_minus][N_minus - 1] * alpha[N_minus - 1]);
    //cout << "T[N_minus] " << T[N_minus] << endl;
    T_next[N_minus] += dT[N_minus];
    //обновляем температуру на узлах, слева обновляем вручную, а dT остается = 0
    for (int i = N_minus - 1; i > 0; i--) {
        dT[i] = alpha[i] * dT[i + 1] + beta[i];
        //cout << "alpha" << i << " =" << alpha[i] << " beta" << i << " =" << beta[i] << "\n";
        T_next[i] += dT[i];
        if (i < Nd + 1 && T_next[i] > T_boiling)
        {
            dT[i] = 0;
            T_next[i] = T_boiling;
        }
        cout << "T_next" << i << " =" << T_next[i] << "\n";
    }
    //T_next[20] = T_next[21];
    T_next[0] = T_next[1];
    //cout << "T[0][k]" << T[0][k] << "\n";
}

void Solver(const int N, vector<double>& r, vector<double>& T_cur, vector<double>& T_next, const double d,
    const double a, const double b, const double c, const int total_time, const double dt, const int const_params, const double dtau,
    const int Nd)
{
    const int N_minus = N - 1;
    vector<double> alpha(N), beta(N), dT(N + 1), nevyaz(N);
    vector < vector <double> > J(N + 1, vector <double>(N + 1));
    ofstream OutNevyazka;
    ofstream OutCurrentTemp;
    ofstream OutStepCurrentTemp;
    //задаем начальную температуру следующего слоя как решение предыдущей
    for (int j = 0; j < N + 1; j++)
    {
        T_next[j] = T_cur[j];
    }    //Reading the initial distribution
    OutCurrentTemp.open("C:/Users/user/source/Научная работа/UnevenGrid/Data/Temp_" + to_string(0) + ".dat");
    OutCurrentTemp << "TITLE=\"" << "Graphics" << "\"" << endl;
    OutCurrentTemp << R"(VARIABLES= "rj", "T", "Lambda" )" << endl;
    for (int i = 0; i < N + 1; i++) {
        //cout << r[i] << endl;
        OutCurrentTemp << r[i] << " " << T_cur[i] << " " << Lambda(T_cur[i], a, b, c) << "\n";
    }
    OutCurrentTemp.close();
    for (int n = 0; n < total_time; n++)
    { 
        //метод Ньютона для n-ого временного слоя
        double mod_nevyaz = 10000;
        int s = 0;
        //Discrepancy for the current layer
        F(nevyaz, T_cur, T_next, N_minus, r, a, b, c, d, dt, const_params);
        //Collecting Nevyazka
        if (n % 10 == 0 || n < 20)
        {
            OutNevyazka.open("C:/Users/user/source/Научная работа/UnevenGrid/Data/Nevyazka_" + to_string(n) + ".dat");
            OutNevyazka << "TITLE=\"" << "Graphics" << "\"" << endl;
            OutNevyazka << R"(VARIABLES= "i", "F")" << endl;
        }
        //Getting Solution using cycle of Newton method
        while (mod_nevyaz > pow(10, -1))
        {
            mod_nevyaz = 0;
            Jacobian(J, r, T_next, N_minus, a, b, c, d, dt, const_params, dtau);
            cout << "s " << s << endl;
            //МЕТОД ПРОГОНКИ: T_next change from s to s+1
            Progonka(T_next, alpha, beta, dT, nevyaz, J, N_minus, Nd);
            //cout << "T0" << " =" << T[0] << endl;
            s += 1;
            F(nevyaz, T_cur, T_next, N_minus, r, a, b, c, d, dt, const_params);
            for (int i = 1; i < N_minus; i++) {
                //cout << "nevyaz[i]" << nevyaz[i] << endl;
                mod_nevyaz += pow(nevyaz[i], 2.0);
            }
            mod_nevyaz = pow(mod_nevyaz, 0.5);
            OutNevyazka << s << " " << mod_nevyaz << endl;
            cout << endl << "nevyazka = " << mod_nevyaz << endl;
            //Collecting Temperature and Lambda
            if (n % 10 == 0 || n < 20)
            {
                OutCurrentTemp.open("C:/Users/user/source/Научная работа/UnevenGrid/Data/Temp_" + to_string(n + 1) + ".dat");
                OutCurrentTemp << "TITLE=\"" << "Graphics" << "\"" << endl;
                OutCurrentTemp << R"(VARIABLES= "rj", "T", "Lambda" )" << endl;
                for (int i = 0; i < N + 1; i++) {
                    //cout << r[i] << endl;
                    OutCurrentTemp << r[i] << " " << T_next[i] << " " << Lambda(T_next[i], a, b, c) << "\n";
                }
                OutCurrentTemp.close();
            }
            //Collecting Step Temperature and Lambda
            if (n == 0)
            {
                OutStepCurrentTemp.open("C:/Users/user/source/Научная работа/UnevenGrid/Data1/StepTemp_" + to_string(s) + ".dat");
                OutStepCurrentTemp << "TITLE=\"" << "Graphics" << "\"" << endl;
                OutStepCurrentTemp << R"(VARIABLES= "rj", "T", "Lambda" )" << endl;
                for (int i = 0; i < N + 1; i++) {
                    //cout << r[i] << endl;
                    OutStepCurrentTemp << r[i] << " " << T_next[i] << " " << Lambda(T_next[i], a, b, c) << "\n";
                    //cout << r[i] << " " << T_next[i] << " " << Lambda(T_next[i], a, b, c) << "\n";
                }
                OutStepCurrentTemp.close();
            }
        }
        OutNevyazka.close();
        //задаем начальную температуру следующего слоя как решение предыдущей
        for (int j = 0; j < N + 1; j++)
        {
            T_cur[j] = T_next[j];
        }
    }
    OutCurrentTemp.open("C:/Users/user/source/Научная работа/UnevenGrid/Data/Temp_" + to_string(total_time) + ".dat");
    OutCurrentTemp << "TITLE=\"" << "Graphics" << "\"" << endl;
    OutCurrentTemp << R"(VARIABLES= "rj", "T", "Lambda" )" << endl;
    for (int i = 0; i < N + 1; i++) {
        //cout << r[i] << endl;
        OutCurrentTemp << r[i] << " " << T_next[i] << " " << Lambda(T_next[i], a, b, c) << "\n";
    }
    OutCurrentTemp.close();
    std::cout << "Done!\n";
}
    
/*
{
    int s;
    const int N_minus = N - 1;
    vector<double> alpha(N), beta(N), dT(N + 1), nevyaz(N);
    vector < vector <double> > J(N + 1, vector <double>(N + 1));

    ofstream OutNevyazka;

    for (int n = 0; n < total_time; n++)
    {
        //метод Ньютона для n-ого временного слоя
        F0(nevyaz, T, N_minus, r, d, a, b, c, n);
        double mod_nevyaz = 10000;
        s = 0;
        while (mod_nevyaz > pow(10, -6))
        {
            mod_nevyaz = 0;
            Jacobian(J, r, T, N_minus, d, a, b, c, n);
            cout << "s " << s << endl;
            //МЕТОД ПРОГОНКИ
            Progonka(T, alpha, beta, dT, nevyaz, J, N_minus, n);
            //КОНЕЦ МЕТОДА ПРОГОНКИ
            s++;
            //cout << "T0" << " =" << T[0] << endl;
            F(nevyaz, T, N_minus, r, d, a, b, c, n);
            for (int i = 1; i < N_minus; i++) {
                //cout << "nevyaz[i]" << nevyaz[i] << endl;
                mod_nevyaz += pow(nevyaz[i], 2.0);
            }
            mod_nevyaz = pow(mod_nevyaz, 0.5);
            OutNevyazka.open("Nevyazka_" + std::to_string(n) + ".dat");
            OutNevyazka << "TITLE=\"" << "Graphics" << "\"" << endl;
            OutNevyazka << R"(VARIABLES= "i", "F")" << endl;
            OutNevyazka << s << " " << mod_nevyaz << endl;
            cout << endl << "nevyazka = " << mod_nevyaz << endl;

            ofstream OutCurrentTemp;
            OutCurrentTemp.open("Temp_" + std::to_string(n) + ".dat");
            OutCurrentTemp << "TITLE=\"" << "Graphics" << "\"" << endl;
            OutCurrentTemp << R"(VARIABLES= "r/l", "T", "Lambda" )" << endl;
            for (int i = 0; i < N + 1; i++) {
                OutCurrentTemp << r[i] << " " << T[i][n] << " " << Lambda(T[i][n], a, b, c) << "\n";
            }
            OutCurrentTemp.close();
        }
    }   
    OutNevyazka.close();
    std::cout << "Done!\n";
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

void InitialGrid(int& N, const double x_l, const double x_r, const double T_l, const double T_r, vector<double>& r,
    vector<double>& T_cur, vector<double>& T_next, double h, const double q, const double h_min, const int big_number, const int const_params,
    const int Nd)
{
    T_cur.resize(N + 1);
    T_next.resize(N + 1);
    r.resize(N + 1);
    ofstream OutX;
    OutX.open("C:/Users/user/source/Научная работа/UnevenGrid/Data/X_grid.dat");
    OutX << "TITLE=\"" << "Graphics" << "\"" << endl;
    OutX << R"(VARIABLES= "j", "hj" )" << endl;
    //Non-uniform 1-dimensional grid 
    if (const_params == 0 || const_params == 2)
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
            cout << "dsfsd" << h << "\n";
            cout << "x0" << " " << r[0] << endl;
            h = h_min;
            int j = 1;
            //Grid for droplet
            for (j; j < Nd + 1; j++)
            {
                r[j] = r[j - 1] + h;
                T_cur[j] = T_l;

                cout << "T" << j << " " << T_cur[j] << endl;
                OutX << j << " " << h << " " << "\n";
                cout << "h" << h << " r" << j << " " << r[j] << endl;
            }
            for (j = Nd + 1; j < N + 1; j++)
            {
                r[j] = r[j - 1] + h;
                h = h * q;
                T_cur[j] = T_r;
                cout << "T" << j << " " << T_cur[j] << endl;
                OutX << j << " " << h << " " << "\n";
                cout << "h" << h << " r" << j << " " << r[j] << endl;
            }
        }
    }
    OutX.close();
}

int main()
{
    const int const_params = 1;          // 0 - r нестац, 1 - r нестац капля и газ, 2 - r нестац адиабатич слева  
    const double dt = 0.001;
    const double dtau = 100000000000000;
    const int total_time = 30;
    const int power_of = 5;              //=6, for a, b, c; =5 for drop - steam constant
    const double a = pow(10, -(power_of + 3)), b = pow(10, -power_of), c = pow(10, -(power_of - 3));
    const double d = 0.1;                          //коэффициент диффузии
    const double T_l = 300;
    const double T_r = 1500;
    double h = 0.2;
    //зададим минимальный возможный размер ячейки
    const double h_min = 0.000025;                 //0.0005 / 20 = 0.000025 = 25 мкр ; 
    //Number of cells
    int N = 100;
    const int Nd = 20;
    const double x_l = 0.0;
    //const double R = 0.05;                       //примерно 0.5 мм
    const double x_r = 1.0;                        //1 примерно равен 10 мм
    double R = x_l + Nd * h_min;                   //зададим радиус капли
    double q = DefineQ(N - Nd, h_min, x_r, R);
    vector <double> r;
    vector <double> T_cur;
    vector <double> T_next;
    GraphicsParameters(a, b, c, T_l, T_r);
    InitialGrid(N, x_l, x_r, T_l, T_r, r, T_cur, T_next, h, q, h_min, total_time, const_params, Nd);
    cout << "N" << N << endl;
    Solver(N, r, T_cur, T_next, d, a, b, c, total_time, dt, const_params, dtau, Nd);
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
