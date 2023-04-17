// UnevenGrid.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

// VariableGrid.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

// NewtonLambdaConst.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//


#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

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
        //считается
    }
}

double Lambda(double T, const double a, const double b, const double c, const int flag_lambda)
{
    if (flag_lambda == 1)
    {
        return 1.;
    }
    else
    {
        return a * pow(T, 2.0) + b * T + c;
    }
}
double DfLambda(double T, const double a, const double b, const int flag_lambda)
{
    if (flag_lambda == 1)
    {
        return 0;
    }
    else
    {
        return 2 * a * T + b;
    }
}
double DfLeft(vector<double>& r, vector<double>& T, const double a, const double b, const double c,
    const int j, const int flag_grid, const int flag_lambda)
{
    if (flag_grid == 0)
    {
        return (2. / (r[j + 1] - r[j - 1])) * (-DfLambda((T[j] + T[j - 1]) / 2., a, b, flag_lambda) * (T[j] - T[j - 1]) / (r[j] - r[j - 1])
            + Lambda((T[j] + T[j - 1]) / 2., a, b, c, flag_lambda) / (r[j] - r[j - 1]));
    }
    else if (flag_grid == 1)
    {
        return (2. / (r[j] - r[j - 1])) * (1. / (r[j + 1] - r[j - 1]) - 1. / r[j]);
    }
    else
    {
        return (2. / (r[j + 1] - r[j - 1]) / pow(r[j], 2.)) * pow((r[j] + r[j - 1]) / 2., 2.)
            * (-DfLambda((T[j] + T[j - 1]) / 2., a, b, flag_lambda) * (T[j] - T[j - 1]) / (r[j] - r[j - 1])
            + Lambda((T[j] + T[j - 1]) / 2., a, b, c, flag_lambda) / (r[j] - r[j - 1]));
    }
}
double DfCenter(vector<double>& r, vector<double>& T, const double d, const double a, const double b, const double c,
    const int j, const int flag_grid, const int flag_lambda)
{
    //double dt = pow(r[j] - r[j - 1], 2.0) * d;
    if (flag_grid == 0)
    {
        return (2. / (r[j + 1] - r[j - 1])) * ((DfLambda((T[j] + T[j + 1]) / 2., a, b, flag_lambda) *(T[j + 1] - T[j]) / (r[j + 1] - r[j])
            - Lambda((T[j] + T[j + 1]) / 2., a, b, c, flag_lambda) / (r[j + 1] - r[j]))
            - (DfLambda((T[j] + T[j - 1]) / 2., a, b, flag_lambda) * (T[j] - T[j - 1]) / (r[j] - r[j - 1])
            + Lambda((T[j] + T[j - 1]) / 2., a, b, c, flag_lambda) / (r[j] - r[j - 1])));
    }
    else if (flag_grid == 1)
    {
        return -(2. / (r[j + 1] - r[j - 1])) * (1. / (r[j + 1] - r[j]) + 1. / (r[j] - r[j - 1])) + 2. / r[j] / (r[j] - r[j - 1]);
    }
    else
    {
        return (2. / (r[j + 1] - r[j - 1]) / pow(r[j], 2.)) * (pow((r[j + 1] + r[j]) / 2., 2.)
            * (DfLambda((T[j] + T[j + 1]) / 2., a, b, flag_lambda) * (T[j + 1] - T[j]) / (r[j + 1] - r[j])
                - Lambda((T[j] + T[j + 1]) / 2., a, b, c, flag_lambda) / (r[j + 1] - r[j])) -
            pow((r[j] + r[j - 1]) / 2., 2.) * (DfLambda((T[j] + T[j - 1]) / 2., a, b, flag_lambda) * (T[j] - T[j - 1]) / (r[j] - r[j - 1])
                + Lambda((T[j] + T[j - 1]) / 2., a, b, c, flag_lambda) / (r[j] - r[j - 1])));
    }
}

double DfRight(vector<double>& r, vector<double>& T, const double a, const double b, const double c,
    const int j, const int flag_grid, const int flag_lambda)
{
    if (flag_grid == 0)
    {
        return (2. / (r[j + 1] - r[j - 1])) * (DfLambda((T[j + 1] + T[j]) / 2., a, b, flag_lambda) * (T[j + 1] - T[j]) / (r[j + 1] - r[j])
            + Lambda((T[j + 1] + T[j]) / 2., a, b, c, flag_lambda) / (r[j + 1] - r[j]));
    }
    else if (flag_grid == 1)
    {
        return 2. / (r[j + 1] - r[j - 1]) / (r[j + 1] - r[j]);
    }
    else
    {
        return (2. / (r[j + 1] - r[j - 1]) / pow(r[j], 2.)) * pow((r[j + 1] + r[j]) / 2., 2.)
            * (DfLambda((T[j + 1] + T[j]) / 2., a, b, flag_lambda)
            * (T[j + 1] - T[j]) / (r[j + 1] - r[j]) + Lambda((T[j + 1] + T[j]) / 2., a, b, c, flag_lambda) / (r[j + 1] - r[j]));
    }
}

void Jacobian(vector<vector<double>>& J, vector<double>& r, vector<double>& T, const int N_minus, const double d,
    const double a, const double b, const double c, const int flag_grid, const int flag_lambda)
{

    J[1][1] = DfCenter(r, T, d, a, b, c, 1, flag_grid, flag_lambda);
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
    J[1][2] = DfRight(r, T, a, b, c, 1, flag_grid, flag_lambda);
    //cout << "J[1][2]" << J[1][2] << endl;
    for (int j = 2; j < N_minus; j++) {
        J[j][j - 1] = DfLeft(r, T, a, b, c, j, flag_grid, flag_lambda);
        J[j][j] = DfCenter(r, T, d, a, b, c, j, flag_grid, flag_lambda);
        J[j][j + 1] = DfRight(r, T, a, b, c, j, flag_grid, flag_lambda);
    }
    J[N_minus][N_minus - 1] = DfLeft(r, T, a, b, c, N_minus, flag_grid, flag_lambda);
    J[N_minus][N_minus] = DfCenter(r, T, d, a, b, c, N_minus, flag_grid, flag_lambda);
}

void F(vector<double>& f, vector<double>& T, const int N_minus, vector<double>& r,
    const double a, const double b, const double c, const int flag_grid, const int flag_lambda)
{
    if (flag_grid == 0)
    {
        for (int j = 1; j < N_minus + 1; j++)
        {
            f[j] = -(2. / (r[j + 1] - r[j - 1])) * (Lambda((T[j + 1] + T[j]) / 2., a, b, c, flag_lambda) * (T[j + 1] - T[j]) /
                (r[j + 1] - r[j]) - Lambda((T[j] + T[j - 1]) / 2., a, b, c, flag_lambda) * (T[j] - T[j - 1]) / (r[j] - r[j - 1]));
            //cout << j << " f" << j << " " << f[j] << endl;
        }
    }
    else if (flag_grid == 1)
    {
        for (int j = 1; j < N_minus + 1; j++)
        {
            f[j] = -((2. / (r[j + 1] - r[j - 1])) * ((T[j + 1] - T[j]) / (r[j + 1] - r[j]) - (T[j] - T[j - 1]) / (r[j] - r[j - 1]))
                + (2. / r[j]) * (T[j] - T[j - 1]) / (r[j] - r[j - 1]));
            cout << j << " f" << j << " " << f[j] << endl;
        }
    }
    else
    {
        for (int j = 1; j < N_minus + 1; j++)
        {
            f[j] = -(2. / (r[j + 1] - r[j - 1]) / pow(r[j], 2.))
                * (pow((r[j + 1] + r[j]) / 2., 2.) * Lambda((T[j + 1] + T[j]) / 2., a, b, c, flag_lambda) * (T[j + 1] - T[j]) / (r[j + 1] - r[j])
                - pow((r[j] + r[j - 1]) / 2., 2.) * Lambda((T[j] + T[j - 1]) / 2., a, b, c, flag_lambda) * (T[j] - T[j - 1]) / (r[j] - r[j - 1]));
            cout << j << " f" << j << " " << f[j] << endl;
        }
    }
}

void Progonka(vector<double>& T, vector<double>& alpha, vector<double>& beta, vector<double>& dT, vector<double>& nevyaz,
    vector < vector <double> >& J, const int N_minus)
{
    //считаем коэффициенты слева
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
    T[N_minus] += dT[N_minus];
    //обновляем температуру на узлах
    cout << "T[N_minus][" << " " << "] =" << T[N_minus] << endl;
    for (int i = N_minus - 1; i >= 1; i--) {
        dT[i] = alpha[i] * dT[i + 1] + beta[i];
        T[i] += dT[i];
        //cout << "T" << i  << " =" << T[i][n] << endl;
    }
}

void AnaliticSolver(const int N, vector<double>& r, vector <double>& T, const double T_l, const double T_r)
{
    //analytical solution for T|_x_l = T1 = 300K, T|_x_r = T2 = 2000K
    //lambda = const
    ofstream OutCurAnTemp;
    OutCurAnTemp.open("AnTemp.dat");
    OutCurAnTemp << "TITLE=\"" << "Graphics" << "\"" << endl;
    OutCurAnTemp << R"(VARIABLES= "rj", "T")" << endl;
    for (int i = 0; i < N + 1; i++)
    {
        T[i] = (T_r * r[N] - T_l * r[0]) / (r[N] - r[0]) + (T_l - T_r) * r[N] * r[0] / ((r[N] - r[0]) * r[i]);
        OutCurAnTemp << r[i] << " " << T[i] << "\n";
    }
    OutCurAnTemp.close();
}

void Solver(const int N, vector<double>& r, vector <double>& T, const double d, const double a, const double b, const double c,
    const int flag_grid, const int flag_lambda)
{
    const int N_minus = N - 1;
    double zn;
    vector<double> alpha(N), beta(N), dT(N + 1), nevyaz(N);
    vector < vector <double> > J(N + 1, vector <double>(N + 1));
    double mod_nevyaz = 10000;
    int iter = 0;
    F(nevyaz, T, N_minus, r, a, b, c, flag_grid, flag_lambda);
    ofstream OutNevyazka;
    OutNevyazka.open("Nevyazka.dat");
    OutNevyazka << "TITLE=\"" << "Graphics" << "\"" << endl;
    OutNevyazka << R"(VARIABLES= "i", "F")" << endl;
    //метод прогонки
    while (mod_nevyaz > pow(10, -6))
    {
        mod_nevyaz = 0;
        Jacobian(J, r, T, N_minus, d, a, b, c, flag_grid, flag_lambda);
        cout << "iter " << iter << endl;
        Progonka(T, alpha, beta, dT, nevyaz, J, N_minus);
        iter += 1;
        F(nevyaz, T, N_minus, r, a, b, c, flag_grid, flag_lambda);
        for (int i = 1; i < N_minus; i++) {
            //cout << "nevyaz[i]" << nevyaz[i] << endl;
            mod_nevyaz += pow(nevyaz[i], 2.0);
        }
        mod_nevyaz = pow(mod_nevyaz, 0.5);
        OutNevyazka << iter << " " << mod_nevyaz << endl;
        cout << endl << "nevyazka = " << mod_nevyaz << endl;

        ofstream OutCurrentTemp;
        OutCurrentTemp.open("Temp.dat");
        OutCurrentTemp << "TITLE=\"" << "Graphics" << "\"" << endl;
        OutCurrentTemp << R"(VARIABLES= "rj", "T", "Lambda" )" << endl;
        for (int i = 0; i < N + 1; i++) {
            //cout << r[i] << endl;
            OutCurrentTemp << r[i] << " " << T[i] << " " << Lambda(T[i], a, b, c, flag_lambda) << "\n";
        }
        OutCurrentTemp.close();
    }
    OutNevyazka.close();
    std::cout << "Done!\n";
}

//Search of q
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

void InitialGrid(int& N, const double x_l, const double x_r, const double T_l, const double T_r,
    vector<double>& r, vector<double>& T, double h, const double q, const double h_min)
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
    r.resize(N + 1);
    ofstream OutX;
    OutX.open("X_grid.dat");
    OutX << "TITLE=\"" << "Graphics" << "\"" << endl;
    OutX << R"(VARIABLES= "j", "hj" )" << endl;
    if (q > 1 || q == 1)
    {
        h = h_min;
        r[0] = x_l;
        T[0] = T_l;
        OutX << 0 << " " << h << " " << "\n";
        cout << "x0" << " " << r[0] << endl;
        for (int j = 1; j < N + 1; j++)
        {
            r[j] = r[j - 1] + h;
            h = h * q;
            if (r[j] <= x_l + 0.001)
                T[j] = T_l;
            else
                T[j] = T_r;
            cout << "T" << j << " " << T[j] << endl;
            OutX << j << " " << h << " " << "\n";
            cout << "h" << h << " r" << j << " " << r[j] << endl;
        }
    }
    if (q < 1)
    {
        r[N] = x_r;
        T[N] = T_r;
        cout << "N" << N << endl;
        for (int j = N - 1; j >= 0; j--)
        {
            cout << "T" << j + 1 << " " << T[j + 1] << endl;
            OutX << j + 1 << " " << h << " " << "\n";
            cout << "h" << h << " r" << j + 1 << " " << r[j + 1] << endl;
            r[j] = r[j + 1] - h;
            h = h * q;
            T[j] = T_r;
        }
        T[0] = T_l;
        cout << "T" << 0 << " " << T[0] << endl;
        OutX << 0 << " " << h << " " << "\n";
        cout << "h" << h << " r" << 0 << " " << r[0] << endl;
        cout << "Area size = " << r[N] - r[0] << endl;
    }
    OutX.close();
}

/*
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
r.resize(N + 1);
ofstream OutX;
OutX.open("X_grid.dat");
OutX << "TITLE=\"" << "Graphics" << "\"" << endl;
OutX << R"(VARIABLES= "j", "hj" )" << endl;
if (q > 1 || q == 1)
{
    r[0] = x_l;
    T[0] = T_l;
    OutX << 0 << " " << h << " " << "\n";
    cout << "x0" << " " << r[0] << endl;
    for (int j = 1; j < N + 1; j++)
    {
        r[j] = r[j - 1] + h;
        h = h * q;
        T[j] = T_r;
        cout << "T" << j << " " << T[j] << endl;
        OutX << j << " " << h << " " << "\n";
        cout << "h" << h << " r" << j << " " << r[j] << endl;
    }
}
if (q < 1)
{
    r[N] = x_r;
    T[N] = T_r;
    cout << "N" << N << endl;
    for (int j = N - 1; j >= 0; j--)
    {
        cout << "T" << j + 1 << " " << T[j + 1] << endl;
        OutX << j + 1 << " " << h << " " << "\n";
        cout << "h" << h << " r" << j + 1 << " " << r[j + 1] << endl;
        r[j] = r[j + 1] - h;
        h = h * q;
        T[j] = T_r;
    }
    T[0] = T_l;
    cout << "T" << 0 << " " << T[0] << endl;
    OutX << 0 << " " << h << " " << "\n";
    cout << "h" << h << " r" << 0 << " " << r[0] << endl;
    cout << "b= " << r[N] - r[0] << endl;

}
OutX.close();
*/

int main()
{
    const int flag_lambda = 0;     //0 - lambda(T), 1 - const 
    const int flag_grid = 2;       //0 - x,1 - r(1st scheme, lambda = const), else - r(2nd scheme)
    const int power_of = 6;
    const double a = pow(10, -(power_of + 3)), b = pow(10, -power_of), c = pow(10, -(power_of - 3));
    const double d = 0.1;     //коэффициент диффузии
    const double T_l = 300;
    const double T_r = 2000;
    //зададим минимальный возможный размер ячейки - первый член возрастающей геометрической прогрессии 
    const double h_min = 0.0001;
    int N = 50;
    double h = 0.1;
    const double x_l = 0.1;
    const double x_r = 1.1;
    double q = DefineQ(N, h_min, x_r, x_l);
    cout << "q" << q << endl;
        //= DefineN(h, q, x_l, x_r);
    //cout << "N " << N << "log1" << endl;
    vector<double> r;
    vector<double> T;
    InitialGrid(N, x_l, x_r, T_l, T_r, r, T, h, q, h_min);
    cout << "N" << N << endl;
    if (flag_grid != 0)
        AnaliticSolver(N, r, T, T_l, T_r);
    Solver(N, r, T, d, a, b, c, flag_grid, flag_lambda);
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

