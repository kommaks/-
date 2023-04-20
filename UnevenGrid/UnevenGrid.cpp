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

//Formula setting for Heat Capacity
double Cp(double T, int a)
{
    const double T_boiling = 373;
    double A, B, C, D, E;
    //Polinom for T = (K) /1000
    double T_forCp = T / 1000.;
    if (T <= T_boiling)
    {
        A = -203.6060;
        B = 1523.290;
        C = -3196.413;
        D = 2474.455;
        E = 3.855326;
    }
    else
    {
        A = 30.09200;
        B = 6.832514;
        C = 6.793435;
        D = -2.534480;
        E = 0.082139;
    }
    //if (a == pow(10, -9))
    //    return 1;
    //else
        return A + B * T_forCp + C * pow(T_forCp, 2.) + D * pow(T_forCp, 3.) + E * pow(T_forCp, -2.);
}

double DfCp(double T, int a)
{
    const double T_boiling = 373;
    double A, B, C, D, E;
    //Polinom for T = (K) /1000
    double T_forCp = T / 1000.;
    if (T <= T_boiling)
    {
        A = -203.6060;
        B = 1523.290;
        C = -3196.413;
        D = 2474.455;
        E = 3.855326;
    }
    else
    {
        A = 30.09200;
        B = 6.832514;
        C = 6.793435;
        D = -2.534480;
        E = 0.082139;
    }
    //if (a == pow(10, -9))
    //    return 1.;
    //else
        return (B + 2. * C * T_forCp + 3. * D * pow(T_forCp, 2.) - 2. * E * pow(T_forCp, -3.)) / 1000.;
}


//Formula setting for Density
double Density(double T, int a) 
{                                              
    const double R = 8.314462;
    const double p = pow(10, 5);
    //if (a == pow(10, -9))
    //    return 1;
    //else {
        return p / (R * T);
    //}
}

double DfDensity(double T, int a)
{
    const double R = 8.314462;
    const double p = pow(10, 5);
    //if (a == pow(10, -9))
    //    return 1;
    //else
        return - p / (R * pow(T, 2.));
}

//Formula setting for thermal conductivity
double Lambda(double T, const double a, const double b, const double c)
{
    if (a == pow(10, -9))
        return a * pow(T, 2.0) + b * T + c;
    /*else
    //{
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
    */
}
double DfLambda(double T, const double a, const double b)
{
    if (a == pow(10, -9))
        return 2 * a * T + b;
    /*else
    //{
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
    */
}
//Setting j-1 element of matrix
double DfLeft(vector<double>& r, vector<vector<double>>& T, const double a, const double b, const double c,
    int j, int k)
{
    return (2. / (r[j + 1] - r[j - 1]) / pow(r[j], 2.)) * pow((r[j] + r[j - 1]) / 2., 2.)
        * (-DfLambda((T[j][k] + T[j - 1][k]) / 2., a, b) * (T[j][k] - T[j - 1][k]) / (r[j] - r[j - 1])
            + Lambda((T[j][k] + T[j - 1][k]) / 2., a, b, c) / (r[j] - r[j - 1]));
}
//Setting j element of matrix
double DfCenter(vector<double>& r, vector<vector<double>>& T, const double a, const double b, const double c, const double d,
    int j, int k)
{
    //Time step = 1
    double dt = 0.1;
    //if (a == pow(10, -9))
    //{
    //    return (2. / (r[j + 1] - r[j - 1]) / pow(r[j], 2.)) * (pow((r[j + 1] + r[j]) / 2., 2.)
    //       * (DfLambda((T[j][k] + T[j + 1][k]) / 2., a, b) * (T[j + 1][k] - T[j][k]) / (r[j + 1] - r[j])
    //            - Lambda((T[j][k] + T[j + 1][k]) / 2., a, b, c) / (r[j + 1] - r[j])) - pow((r[j] + r[j - 1]) / 2., 2.)
    //        * (DfLambda((T[j][k] + T[j - 1][k]) / 2., a, b) * (T[j][k] - T[j - 1][k]) / (r[j] - r[j - 1])
    //            + Lambda((T[j][k] + T[j - 1][k]) / 2., a, b, c) / (r[j] - r[j - 1]))) - 1. / dt;
    //}
    //else
    //{
    return (2. / (r[j + 1] - r[j - 1]) / pow(r[j], 2.)) * (pow((r[j + 1] + r[j]) / 2., 2.) * (DfLambda((T[j][k] + T[j + 1][k]) / 2., a, b) *
        (T[j + 1][k] - T[j][k]) / (r[j + 1] - r[j]) - Lambda((T[j][k] + T[j + 1][k]) / 2., a, b, c) / (r[j + 1] - r[j])) -
        pow((r[j] + r[j - 1]) / 2., 2.) * (DfLambda((T[j][k] + T[j - 1][k]) / 2., a, b) * (T[j][k] - T[j - 1][k]) / (r[j] - r[j - 1])
            + Lambda((T[j][k] + T[j - 1][k]) / 2., a, b, c) / (r[j] - r[j - 1]))) - 1 / dt;
    //}
}
//Setting j+1 element of matrix
double DfRight(vector<double>& r, vector<vector<double>>& T, const double a, const double b, const double c,
    int j, int k)
{
    return (2. / (r[j + 1] - r[j - 1]) / pow(r[j], 2.)) * pow((r[j + 1] + r[j]) / 2., 2.)
        * (DfLambda((T[j + 1][k] + T[j][k]) / 2., a, b) * (T[j + 1][k] - T[j][k]) / (r[j + 1] - r[j])
            + Lambda((T[j + 1][k] + T[j][k]) / 2., a, b, c) / (r[j + 1] - r[j]));
}

double AdiabaticBoundaryDfs(vector<double>& r, vector<vector<double>>& T, const double a, const double b, const double c,
    const int j, const int k)
{
    double dt = 0.1;
    return (2. / (r[2] - r[0]) / pow(r[1], 2.)) * (pow((r[2] + r[1]) / 2., 2.) * (DfLambda((T[1][k] + T[2][k]) / 2., a, b) *
        (T[2][k] - T[1][k]) / (r[2] - r[1]) - Lambda((T[1][k] + T[2][k]) / 2., a, b, c) / (r[2] - r[1]))) - 1 / dt;
}
void Jacobian(vector<vector<double>>& J, vector<double>& r, vector<vector<double>>& T, const int N_minus, const double a,
    const double b, const double c, const double d, int n, const int const_params)
{
    if (const_params == 0)
    {
            J[1][1] = DfCenter(r, T, a, b, c, d, 1, n + 1);
            /*
            cout << "J[1][1]" << J[1][1] << endl;
            cout << "DfL" << DfLambda((T[2] + T[1]) / 2., a, b) << endl;
            cout << "L" << Lambda((T[2] + T[1]) / 2., a, b, c) << endl;
            cout << "r[0]" << r[0] << endl;
            cout << "r[1]" << r[1] << endl;
            cout << "r[2]" << r[2] << endl;
            cout << "T[1]" << T[1] << endl;
            cout << "T[2]" << T[2] << endl;
            */
            J[1][2] = DfRight(r, T, a, b, c, 1, n + 1);
            //cout << "J[1][2]" << J[1][2] << endl;
            for (int j = 2; j < N_minus; j++) {
                J[j][j - 1] = DfLeft(r, T, a, b, c, j, n + 1);
                J[j][j] = DfCenter(r, T, a, b, c, d, j, n + 1);
                J[j][j + 1] = DfRight(r, T, a, b, c, j, n + 1);
            }
            J[N_minus][N_minus - 1] = DfLeft(r, T, a, b, c, N_minus, n + 1);
            J[N_minus][N_minus] = DfCenter(r, T, a, b, c, d, N_minus, n + 1);
    }
    else if (const_params == 2 || const_params == 1)
    {
        //Jacobians for left boundary
        J[1][1] = AdiabaticBoundaryDfs(r, T, a, b, c, 1, n + 1);
        J[1][2] = DfRight(r, T, a, b, c, 1, n + 1);
        for (int j = 2; j < N_minus; j++) {
            J[j][j - 1] = DfLeft(r, T, a, b, c, j, n + 1);
            J[j][j] = DfCenter(r, T, a, b, c, d, j, n + 1);
            J[j][j + 1] = DfRight(r, T, a, b, c, j, n + 1);
        }
        J[N_minus][N_minus - 1] = DfLeft(r, T, a, b, c, N_minus, n + 1);
        J[N_minus][N_minus] = DfCenter(r, T, a, b, c, d, N_minus, n + 1);
    }
}

void F(vector<double>& f, vector<vector<double>>& T, const int N_minus, vector<double>& r, const double a,
    const double b, const double c, const double d, const int n, const int const_params)
{
    //Time step = 1
    double dt = 0.1;
    if (const_params == 0)
    {
        cout << "Dens" << Density(T[0][n + 1], a) << "\n";
        cout << "Cp" << Cp(T[0][n + 1], a) << "\n";
        cout << "Lambda" << Lambda((T[1][n + 1] + T[1][n + 1]) / 2., a, b, c) << "\n";
        for (int j = 1; j < N_minus + 1; j++)
        {
            //cout << "(T[j][n + 1] - T[j][n])" << (T[j][n + 1] - T[j][n]) << endl;
            f[j] = -((2. / (r[j + 1] - r[j - 1]) / pow(r[j], 2.)) * (pow((r[j + 1] + r[j]) / 2., 2.) * Lambda((T[j + 1][n + 1] + T[j][n + 1]) / 2., a, b, c)
                * (T[j + 1][n + 1] - T[j][n + 1]) / (r[j + 1] - r[j]) - pow((r[j] + r[j - 1]) / 2., 2.)
                * Lambda((T[j][n + 1] + T[j - 1][n + 1]) / 2., a, b, c) * (T[j][n + 1] - T[j - 1][n + 1]) 
                / (r[j] - r[j - 1])) - (T[j][n + 1] - T[j][n]) / dt);
            
            //cout << "a =" << a << "\n";

            //cout << j << " f" << j << " " << f[j] << endl;
            /*
            cout << "T[" << j << "][n + 1]" << T[j][n + 1] << "\n";
            cout << "Cp" << Cp(T[j][n + 1]) << "\n";
            cout << "DfCp" << DfCp(T[j][n + 1]) << "\n";
            cout << "Density" << Density(T[j][n + 1]) << "\n";
            cout << "DfDensity" << DfDensity(T[j][n + 1]) << "\n";
            cout << "Lambda" << Lambda((T[j][n + 1] + T[j - 1][n + 1]) / 2., a, b, c) << "\n";
            cout << "DfLambda" << DfLambda((T[j][n + 1] + T[j - 1][n + 1]) / 2., a, b) << "\n";
            */
        }
    }
    else if (const_params == 2 || const_params == 1)
    {
        f[1] = -((2. / (r[2] - r[0]) / pow(r[1], 2.)) * (pow((r[2] + r[1]) / 2., 2.)
            * Lambda((T[2][n + 1] + T[1][n + 1]) / 2., a, b, c)
            * (T[2][n + 1] - T[1][n + 1]) / (r[2] - r[1])) - (T[1][n + 1] - T[1][n]) / dt);
        //cout << "f[0]" << f[0] << "\n";
        for (int j = 2; j < N_minus + 1; j++)
        {
            //cout << "(T[j][n + 1] - T[j][n])" << (T[j][n + 1] - T[j][n]) << endl;
            f[j] = -((2. / (r[j + 1] - r[j - 1]) / pow(r[j], 2.)) * (pow((r[j + 1] + r[j]) / 2., 2.) * Lambda((T[j + 1][n + 1] + T[j][n + 1]) / 2., a, b, c)
                * (T[j + 1][n + 1] - T[j][n + 1]) / (r[j + 1] - r[j]) - pow((r[j] + r[j - 1]) / 2., 2.)
                * Lambda((T[j][n + 1] + T[j - 1][n + 1]) / 2., a, b, c) * (T[j][n + 1] - T[j - 1][n + 1]) / (r[j] - r[j - 1]))
                - (T[j][n + 1] - T[j][n]) / dt);
            //cout << j << " f" << j << " " << f[j] << endl;
        }
    }
}

void Progonka(vector<vector<double>>& T, vector<double>& alpha, vector<double>& beta, vector<double>& dT, vector<double>& nevyaz, 
    vector < vector <double> >& J, const int N_minus, int k)
{
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
    T[N_minus][k] += dT[N_minus];
    //обновляем температуру на узлах
    cout << "T[N_minus][" << k << "] =" << T[N_minus][k] << endl;
    for (int i = N_minus - 1; i >= 1; i--) {
        dT[i] = alpha[i] * dT[i + 1] + beta[i];
        T[i][k] += dT[i];
        //cout << "T" << i  << " =" << T[i][k] << endl;
    }
    T[0][k] = T[1][k];
    //cout << "T[0][k]" << T[0][k] << "\n";
}

void Solver(const int N, vector<double>& r, vector<vector<double>>& T, const double d,
    const double a, const double b, const double c, const int total_time, const int const_params)
{
    const int N_minus = N - 1;
    vector<double> alpha(N), beta(N), dT(N + 1), nevyaz(N);
    vector < vector <double> > J(N + 1, vector <double>(N + 1));
    ofstream OutNevyazka;
    ofstream OutStepCurrentTemp;
    ofstream OutCurrentTemp;
    for (int n = 0; n < total_time; n++)
    { 
        //задаем начальную температуру следующего слоя как решение предыдущей
        for (int j = 0; j < N + 1; j++)
        {
            T[j][n + 1] = T[j][n];
        }
        //метод Ньютона для n-ого временного слоя
        double mod_nevyaz = 10000;
        int s = 0;
        F(nevyaz, T, N_minus, r, a, b, c, d, n, const_params);
        //Collecting Nevyazka
        if (n % 10 == 0 || n < 10)
        {
            OutNevyazka.open("C:/Users/user/source/Научная работа/UnevenGrid/Data/Nevyazka_" + to_string(n) + ".dat");
            OutNevyazka << "TITLE=\"" << "Graphics" << "\"" << endl;
            OutNevyazka << R"(VARIABLES= "i", "F")" << endl;
        }
        OutCurrentTemp.open("C:/Users/user/source/Научная работа/UnevenGrid/Data/Temp_" + to_string(0) + ".dat");
        OutCurrentTemp << "TITLE=\"" << "Graphics" << "\"" << endl;
        OutCurrentTemp << R"(VARIABLES= "rj", "T", "Lambda" )" << endl;
        for (int i = 0; i < N + 1; i++) {
            //cout << r[i] << endl;
            OutCurrentTemp << r[i] << " " << T[i][0] << " " << Lambda(T[i][0], a, b, c) << "\n";
        }
        OutCurrentTemp.close();
        //Getting Solution using cycle of Newton method
        while (mod_nevyaz > pow(10, -6))
        {
            //Collecting Nevyazka
            /*
            if (s % 10 == 0 || s < 10)
            {
                OutStepNevyazka.open("C:/Users/user/source/Научная работа/UnevenGrid/Data/Nevyazka_Step" + to_string(s) + ".dat");
                OutStepNevyazka << "TITLE=\"" << "Graphics" << "\"" << endl;
                OutStepNevyazka << R"(VARIABLES= "i", "F")" << endl;
            }
            */
            mod_nevyaz = 0;
            Jacobian(J, r, T, N_minus, a, b, c, d, n, const_params);
            cout << "s " << s << endl;
            //МЕТОД ПРОГОНКИ
            Progonka(T, alpha, beta, dT, nevyaz, J, N_minus, n + 1);
            //cout << "T0" << " =" << T[0] << endl;
            s += 1;
            F(nevyaz, T, N_minus, r, a, b, c, d, n, const_params);
            for (int i = 1; i < N_minus; i++) {
                //cout << "nevyaz[i]" << nevyaz[i] << endl;
                mod_nevyaz += pow(nevyaz[i], 2.0);
            }
            mod_nevyaz = pow(mod_nevyaz, 0.5);
            OutNevyazka << s << " " << mod_nevyaz << endl;
            cout << endl << "nevyazka = " << mod_nevyaz << endl;
            //Collecting Tempereature and Lambda
            if (n % 10 == 0 || n < 10)
            {
                OutCurrentTemp.open("C:/Users/user/source/Научная работа/UnevenGrid/Data/Temp_" + to_string(n + 1) + ".dat");
                OutCurrentTemp << "TITLE=\"" << "Graphics" << "\"" << endl;
                OutCurrentTemp << R"(VARIABLES= "rj", "T", "Lambda" )" << endl;
                for (int i = 0; i < N + 1; i++) {
                    //cout << r[i] << endl;
                    OutCurrentTemp << r[i] << " " << T[i][n + 1] << " " << Lambda(T[i][n + 1], a, b, c) << "\n";
                }
                OutCurrentTemp.close();
            }
            if (s % 10 == 0 || s < 10)
            {
                OutStepCurrentTemp.open("C:/Users/user/source/Научная работа/UnevenGrid/Data/Temp_Step" + to_string(s) + ".dat");
                OutStepCurrentTemp << "TITLE=\"" << "Graphics" << "\"" << endl;
                OutStepCurrentTemp << R"(VARIABLES= "rj", "T", "Lambda" )" << endl;
                for (int i = 0; i < N + 1; i++) {
                    //cout << r[i] << endl;
                    OutStepCurrentTemp << r[i] << " " << T[i][n + 1] << " " << Lambda(T[i][n + 1], a, b, c) << "\n";
                }
                OutStepCurrentTemp.close();
            }
        }
        OutNevyazka.close();
    }
    OutCurrentTemp.open("C:/Users/user/source/Научная работа/UnevenGrid/Data/Temp_" + to_string(total_time) + ".dat");
    OutCurrentTemp << "TITLE=\"" << "Graphics" << "\"" << endl;
    OutCurrentTemp << R"(VARIABLES= "rj", "T", "Lambda" )" << endl;
    for (int i = 0; i < N + 1; i++) {
        //cout << r[i] << endl;
        OutCurrentTemp << r[i] << " " << T[i][total_time] << " " << Lambda(T[i][total_time], a, b, c) << "\n";
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

float f(double q, const double N, const double h_min, const double x_r, const double x_l) //возвращает значение функции f(q) = ...

{
    return h_min * pow(q, N) - (x_r - x_l) * q - h_min + (x_r - x_l);
}

float df(float q, const double N, const double h_min, const double x_r, const double x_l) //возвращает значение производной

{
    return h_min * N * pow(q, N - 1) - (x_r - x_l);
}

float d2f(float q) // значение второй производной

{
    return 1.;
}

double DefineQ(const double N, const double h_min, const double x_r, const double x_l)
{
    int i = 0;//переменные для расчета итерации
    double x0, xn;// вычисляемые приближения для корня
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
    vector<vector<double>>& T, double h, const double q, const double h_min, const int big_number, const int const_params, const int Nd)
{
    T.resize(N + 1, vector<double>(big_number + 1));;
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
            T[0][0] = T_l;
            OutX << 0 << " " << h << " " << "\n";
            cout << "x0" << " " << r[0] << endl;
            for (int j = 1; j < N + 1; j++)
            {
                r[j] = r[j - 1] + h;
                h = h * q;
                if (r[j] <= x_l + 0.001)
                    T[j][0] = T_l;
                else
                    T[j][0] = T_r;
                cout << "T" << j << " " << T[j][0] << endl;
                OutX << j << " " << h << " " << "\n";
                cout << "h" << h << " r" << j << " " << r[j] << endl;
            }
        }
        if (q < 1)
        {
            r[N] = x_r;
            T[N][0] = T_r;
            cout << "N" << N << endl;
            for (int j = N - 1; j >= 0; j--)
            {
                cout << "T" << j + 1 << " " << T[j + 1][0] << endl;
                OutX << j + 1 << " " << h << " " << "\n";
                cout << "h" << h << " r" << j + 1 << " " << r[j + 1] << endl;
                r[j] = r[j + 1] - h;
                h = h * q;
                T[j][0] = T_r;
            }
            T[0][0] = T_l;
            cout << "T" << 0 << " " << T[0][0] << endl;
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
            T[0][0] = T_l;
            OutX << 0 << " " << h << " " << "\n";
            cout << "x0" << " " << r[0] << endl;
            h = h_min;
            int j = 1;
            //Grid for droplet
            for (j; j < Nd + 1; j++)
            {
                r[j] = r[j - 1] + h;
                T[j][0] = T_l;

                cout << "T" << j << " " << T[j][0] << endl;
                OutX << j << " " << h << " " << "\n";
                cout << "h" << h << " r" << j << " " << r[j] << endl;
            }
            for (j = Nd + 1; j < N + 1; j++)
            {
                r[j] = r[j - 1] + h;
                h = h * q;
                T[j][0] = T_r;
                cout << "T" << j << " " << T[j][0] << endl;
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
    const int total_time = 1000;
    const int power_of = 6;
    const double a = pow(10, -(power_of + 3)), b = pow(10, -power_of), c = pow(10, -(power_of - 3));
    const double d = 0.1;                          //коэффициент диффузии
    const double T_l = 400;
    const double T_r = 1500;
    double h = 0.2;
    //зададим минимальный возможный размер ячейки
    const double h_min = 0.0025;                   //0.05 / 20 = 0.0025 ; 
    int N = 100;
    const int Nd = 20;
    const double x_l = 0.0;
    //const double R = 0.05;                       //примерно 0.5 мм
    const double x_r = 1.0;                        //1 примерно равен 10 мм
    double R = x_l + Nd * h_min;                   //зададим радиус капли
    double q = DefineQ(N - Nd, h_min, x_r, R);
    vector<double> r;
    vector < vector <double> > T;
    InitialGrid(N, x_l, x_r, T_l, T_r, r, T, h, q, h_min, total_time, const_params, Nd);
    cout << "N" << N << endl;
    Solver(N, r, T, d, a, b, c, total_time, const_params);
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
