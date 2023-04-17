// NewtonLambdaConst.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//


#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

double Lambda(double T, const double a, const double b, const double c, const int const_params)
{
    if (const_params == 0)
    {
        return 1.0;
    }
    if (const_params == 1 || const_params == 2)
    {
        return a * pow(T, 2.0) + b * T + c;
    }
}
double DfLambda(double T, const double a, const double b)
{
    return 2 * a * T + b;
}
double DfLeft(vector<double>& x, vector<double>& T, double h, const double a, const double b, const double c,
    const int j, const int c_params)
{
    if (c_params == 1)
    {
        return (1.0 / pow(h, 2.0)) * (-DfLambda((T[j] + T[j - 1]) / 2.0, a, b) * (T[j] - T[j - 1])
            + Lambda((T[j] + T[j - 1]) / 2.0, a, b, c, c_params));
    }
    if (c_params == 2)
    {
        return (2. / (x[j + 1] - x[j - 1])) * (-DfLambda((T[j] + T[j - 1]) / 2., a, b) * (T[j] - T[j - 1]) / (x[j] - x[j - 1])
            + Lambda((T[j] + T[j - 1]) / 2., a, b, c, c_params) / (x[j] - x[j - 1]));
    }
}
double DfCenter(vector<double>& x, vector<double>& T, double h, const double a, const double b, const double c,
    const int j, const int c_params)
{
    if (c_params == 1)
    {
        return (1.0 / pow(h, 2.0)) * (DfLambda((T[j + 1] + T[j]) / 2.0, a, b) * (T[j + 1] - T[j])
            - Lambda((T[j + 1] + T[j]) / 2.0, a, b, c, c_params) - DfLambda((T[j] + T[j - 1]) / 2.0, a, b) * (T[j] - T[j - 1])
            - Lambda((T[j] + T[j - 1]) / 2.0, a, b, c, c_params));
    }
    if (c_params == 2)
    {
        return (2. / (x[j + 1] - x[j - 1])) * (DfLambda((T[j] + T[j + 1]) / 2., a, b) * (T[j + 1] - T[j]) / (x[j + 1] - x[j])
            - Lambda((T[j] + T[j + 1]) / 2., a, b, c, c_params) / (x[j + 1] - x[j]) - DfLambda((T[j] + T[j - 1]) / 2., a, b)
            * (T[j] - T[j - 1]) / (x[j] - x[j - 1]) - Lambda((T[j] + T[j - 1]) / 2., a, b, c, c_params) / (x[j] - x[j - 1]));
    }
}

double DfRight(vector<double>& x, vector<double>& T, double h, const double a, const double b, const double c,
    const int j, const int c_params)
{
    if (c_params == 1)
    {
        return (1.0 / pow(h, 2.0)) * (DfLambda((T[j + 1] + T[j]) / 2.0, a, b) * (T[j + 1] - T[j])
            + Lambda((T[j + 1] + T[j]) / 2.0, a, b, c, c_params));
    }
    if (c_params == 2)
    {
        return (2. / (x[j + 1] - x[j - 1])) * (DfLambda((T[j + 1] + T[j]) / 2., a, b) * (T[j + 1] - T[j]) / (x[j + 1] - x[j])
            + Lambda((T[j + 1] + T[j]) / 2., a, b, c, c_params) / (x[j + 1] - x[j]));
    }
}

void Jacobian(vector<vector<double>>& J, vector<double>& x, vector<double>& T, const int N_minus, double h,
    const double a, const double b, const double c, const int const_params)
{
    if (const_params == 0)
    {
        J[1][1] = -2.0 / pow(h, 2.0);
        J[1][2] = 1.0 / pow(h, 2.0);

        for (int i = 2; i < N_minus; i++) {
            J[i][i - 1] = 1.0 / pow(h, 2.0);
            J[i][i] = -2.0 / pow(h, 2.0);
            J[i][i + 1] = 1.0 / pow(h, 2.0);
        }
        J[N_minus][N_minus - 1] = 1.0 / pow(h, 2.0);
        J[N_minus][N_minus] = -2.0 / pow(h, 2.0);
    }
    if (const_params == 1 || const_params == 2)
    {
        J[1][1] = DfCenter(x, T, h, a, b, c, 1, const_params);
        J[1][2] = DfRight(x, T, h, a, b, c, 1, const_params);
        for (int j = 2; j < N_minus; j++) {
            J[j][j - 1] = DfLeft(x, T, h, a, b, c, j, const_params);
            J[j][j] = DfCenter(x, T, h, a, b, c, j, const_params);
            J[j][j + 1] = DfRight(x, T, h, a, b, c, j, const_params);
        }
        J[N_minus][N_minus - 1] = DfLeft(x, T, h, a, b, c, N_minus, const_params);
        J[N_minus][N_minus] = DfCenter(x, T, h, a, b, c, N_minus, const_params);
    }
}

void F(vector<double>& f, vector<double>& x, vector<double>& T, const int N_minus, double h,
    const double a, const double b, const double c, const int const_params)
{
    if (const_params == 0)
    {
        for (int j = 1; j < N_minus + 1; j++)
        {
            f[j] = (1.0 / pow(h, 2.0)) * (T[j - 1] - 2 * T[j] + T[j + 1]);
            //cout << j << " f" << j << " " << f[j] << endl;
        }

    }
    if (const_params == 1)
    {
        for (int j = 1; j < N_minus + 1; j++)
        {
            f[j] = - (1.0 / pow(h, 2.0)) * (Lambda((T[j + 1] + T[j]) / 2.0, a, b, c, const_params) * (T[j + 1] - T[j])
                - Lambda((T[j] + T[j - 1]) / 2.0, a, b, c, const_params) * (T[j] - T[j - 1]));
            //cout << j << " f" << j << " " << f[j] << endl;
        }
    }
    if (const_params == 2)
    {
        for (int j = 1; j < N_minus + 1; j++)
        {
            f[j] = - (2.0 / (x[j + 1] - x[j - 1])) * ((Lambda((T[j + 1] + T[j]) / 2.0, a, b, c, const_params) * (T[j + 1] - T[j])) /
                (x[j + 1] - x[j]) - (Lambda((T[j] + T[j - 1]) / 2.0, a, b, c, const_params) * (T[j] - T[j - 1]))/(x[j] - x[j - 1]));
            //cout << j << " f" << j << " " << f[j] << endl;
        }
    }
}

void Solver(const int N, vector<double>& x, vector <double>& T, const double a, const double b, const double c,
    double h, const int const_params)
{
    const int N_minus = N - 1;
    double zn;
    vector<double> alpha(N), beta(N), dT(N + 1), nevyaz(N);
    vector < vector <double> > J(N + 1, vector <double>(N + 1));
    double mod_nevyaz = 10000;
    int iter = 0;
    F(nevyaz, x, T, N_minus, h, a, b, c, const_params);
    ofstream OutNevyazka;
    OutNevyazka.open("Nevyazka.dat");
    OutNevyazka << "TITLE=\"" << "Graphics" << "\"" << endl;
    OutNevyazka << R"(VARIABLES= "i", "F")" << endl;
    //метод прогонки
    while (mod_nevyaz > pow(10, -6))
    {
        mod_nevyaz = 0;
        Jacobian(J, x, T, N_minus, h, a, b, c, const_params);
        cout << "iter " << iter << endl;
        //считаем коэффициенты слева
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
        cout << "T[N_minus] " << T[N_minus] << endl;
        T[N_minus] += dT[N_minus];
        //обновляем температуру на узлах
        cout << "T[N_minus]" << " =" << T[N_minus] << endl;
        for (int i = N_minus - 1; i >= 1; i--) {
            dT[i] = alpha[i] * dT[i + 1] + beta[i];
            T[i] += dT[i];
            //cout << "T" << i  << " =" << T[i] << endl;
        }
        //cout << "T0" << " =" << T[0] << endl;
        iter += 1;
        F(nevyaz, x, T, N_minus, h,  a, b, c, const_params);
        for (int i = 1; i < N_minus; i++) {
            //cout << "nevyaz[i]" << nevyaz[i] << endl;
            mod_nevyaz += pow(nevyaz[i], 2.0);
        }
        mod_nevyaz = pow(mod_nevyaz, 0.5);
        OutNevyazka << iter << " " << mod_nevyaz << endl;
        cout << "nevyazka = " << mod_nevyaz << endl;
        cout << endl;

        ofstream OutCurrentTemp;
        OutCurrentTemp.open("Temp.dat");
        OutCurrentTemp << "TITLE=\"" << "Graphics" << "\"" << endl;
        OutCurrentTemp << R"(VARIABLES= "x", "T", "Lambda" )" << endl;
        for (int i = 0; i < N + 1; i++) {
            OutCurrentTemp << x[i] << " " << T[i] << " " << Lambda(T[i], a, b, c, const_params) << "\n";
        }
        OutCurrentTemp.close();
    }
    OutNevyazka.close();
    std::cout << "Done!\n";
}
void InitialGrid(int& N, const double x_l, const double x_r, const double T_l, const double T_r,
    vector<double>& x, vector<double>& T, double h, const double q, const double h_min)
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
    T.resize(N + 1);
    x.resize(N + 1);
    ofstream OutX;
    OutX.open("X_grid.dat");
    OutX << "TITLE=\"" << "Graphics" << "\"" << endl;
    OutX << R"(VARIABLES= "j", "hj" )" << endl;
    if (q > 1 || q == 1)
    {
        x[0] = x_l;
        T[0] = T_l;
        OutX << 0 << " " << h << " " << "\n";
        cout << "x0" << " " << x[0] << endl;
        for (int j = 1; j < N + 1; j++)
        {
            x[j] = x[j - 1] + h;
            h = h * q;
            T[j] = T_r;
            cout << "T" << j << " " << T[j] << endl;
            OutX << j << " " << h << " " << "\n";
            cout << "h" << h << " x" << j << " " << x[j] << endl;
        }
    }
    if (q < 1)
    {
        x[N] = x_r;
        T[N] = T_r;
        cout << "N" << N << endl;
        for (int j = N - 1; j >= 0; j--)
        {
            cout << "T" << j + 1 << " " << T[j + 1] << endl;
            OutX << j + 1 << " " << h << " " << "\n";
            cout << "h" << h << " x" << j + 1 << " " << x[j + 1] << endl;
            x[j] = x[j + 1] - h;
            h = h * q;
            T[j] = T_r;
        }
        T[0] = T_l;
        cout << "T" << 0 << " " << T[0] << endl;
        OutX << 0 << " " << h << " " << "\n";
        cout << "h" << h << " x" << 0 << " " << x[0] << endl;
        cout << "b= " << x[N] - x[0] << endl;

    }
    OutX.close();
}

int main()
{
    const int const_params = 2;
    const int power_of = 6;
    const double a = pow(10, -(power_of + 3)), b = pow(10, -power_of), c = pow(10, -(power_of - 3));
    const double T_l = 300;
    const double T_r = 2000;
    //Вводим первый член геометрической прогрессии
    //для q < 1: h > (x_r - x_l) * (1 - q). Чтобы гарантированно: h >= (x_r - x_l)
    double h = 0.2;
    //зададим минимальный возможный размер ячейки
    const double h_min = 0.005;
    // использовать q >= 0.7
    const double q = 0.8;
    const double x_l = 0;
    const double x_r = 1;
    int N = 0;
    //= DefineN(h, q, x_l, x_r);
//cout << "N " << N << "log1" << endl;
    vector<double> x;
    vector<double> T;
    InitialGrid(N, x_l, x_r, T_l, T_r, x, T, h, q, h_min);
    cout << "N" << N << endl;
    Solver(N, x, T, a, b, c, h, const_params);
}

//ранняя версия для const_params = 0
/*
void Jacobian(vector<vector<double>>& J, const int N_minus, double h)
{
    J[1][1] = -2.0 / pow(h, 2.0);
    J[1][2] = 1.0 / pow(h, 2.0);

    for (int i = 2; i < N_minus; i++) {
        J[i][i - 1] = 1.0 / pow(h, 2.0);
        J[i][i] = -2.0 / pow(h, 2.0);
        J[i][i + 1] = 1.0 / pow(h, 2.0);
    }
    J[N_minus][N_minus - 1] = 1.0 / pow(h, 2.0);
    J[N_minus][N_minus] = -2.0 / pow(h, 2.0);
}

void F(vector<double>& f, vector<double>& T, const int N_minus, double h) {
    for (int j = 1; j < N_minus; j++)          // то есть до N
    {
        f[j] = (1.0 / pow(h, 2.0)) * (T[j + 1] - 2 * T[j] + T[j - 1]);
        //cout << j << "f[j]" << f[j] << endl;
    }
}

void Solver(const int N, vector<double>& x, vector <double>& T)
{
    const int N_minus = N - 1;
    double zn;
    double h = 1.0;
    vector<double> alpha(N), beta(N), dT(N + 1), nevyaz(N);
    vector < vector <double> > J(N + 1, vector <double>(N + 1));
    double mod_nevyaz = 10000.0;
    int iter = 0;
    F(nevyaz, T, N_minus, h);

    //метод прогонки
    while (mod_nevyaz > pow(10, -6))
    {
        mod_nevyaz = 0;
        Jacobian(J, N_minus, h);
        cout << "iter " << iter << endl;
        //считаем коэффициенты слева
        zn = J[1][1];
        alpha[1] = -J[1][2] / zn;
        beta[1] = nevyaz[1] / zn;
        //считаем коэффициенты для всех узлов
        for (int i = 2; i < N_minus; i++) {
            zn = J[i][i] + J[i][i - 1] * alpha[i - 1];
            alpha[i] = -J[i][i + 1] / zn;
            beta[i] = (nevyaz[i] - J[i][i - 1] * beta[i - 1]) / zn;
        }
        dT[N_minus] = (nevyaz[N_minus] - J[N_minus][N_minus - 1] * beta[N_minus - 1]) / (J[N_minus][N_minus] + J[N_minus][N_minus - 1] * alpha[N_minus - 1]);
        T[N_minus] -= dT[N_minus];
        //обновляем температуру на узлах
        cout << "T[N_minus]" << " =" << T[N_minus] << endl;
        for (int i = N_minus - 1; i > 0; i--) {
            dT[i] = alpha[i] * dT[i + 1] + beta[i];
            T[i] -= dT[i];
            cout << "T" << i << " =" << T[i] << endl;
        }
        cout << "T0" << " =" << T[0] << endl;
        iter += 1;
        F(nevyaz, T, N_minus, h);

        for (int i = 1; i < N_minus; i++) {
            //cout << "nevyaz[i]" << nevyaz[i] << endl;
            mod_nevyaz += pow(nevyaz[i], 2.0);
        }
        mod_nevyaz = pow(mod_nevyaz, 0.5);

        ofstream OutCurrentTemp;
        OutCurrentTemp.open("Temp.dat");
        OutCurrentTemp << "TITLE=\"" << "Graphics" << "\"" << endl;
        OutCurrentTemp << R"(VARIABLES= "x", "T" )" << endl;
        for (int i = 0; i < N + 1; i++) {
            OutCurrentTemp << x[i] << " " << T[i] << "\n";
        }
        OutCurrentTemp.close();

        cout << endl << "nevyazka = " << mod_nevyaz << " ";
        std::cout << std::endl << endl;

    }

    std::cout << "Done!\n";
}

void Init(const int N, const double T_l, const double T_r, vector<double>& x, vector<double>& T)
{
    T.resize(N + 1);
    x.resize(N + 1);
    x[0] = 0;
    double h = 1.0;
    T[0] = T_l;
    for (int i = 1; i <= N; i++)
    {
        T[i] = T_r;
        x[i] = x[i - 1] + h;
    }

}

int main()
{
    const int N = 50;
    const double T_l = 300;
    const double T_r = 2000;
    vector<double> x;
    vector<double> T;
    Init(N, T_l, T_r, x, T);
    Solver(N, x, T);

}
*/
// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
