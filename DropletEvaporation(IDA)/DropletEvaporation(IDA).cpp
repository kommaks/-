
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


#define FTOL   RCONST(1.e-12) /* function tolerance */
#define STOL   RCONST(1.e-12) /* step tolerance     */

#define ZERO   RCONST(0.0)
#define PT25   RCONST(0.25)
#define PT5    RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)
#define TWO    RCONST(2.0)

#define PI     RCONST(3.1415926)
#define E      RCONST(2.7182818)

static int num_gas_species = 9;
static int num_react = 22;
int cycle = 0;

double Y_N2 = 0.745187;
double Y_max = 1 - Y_N2;
double P = 0.101325;
double A = 6.85 * pow(10, 12);
double R = 8.314;
double Ea = 46.37 * 293.;
double koeff_l = 0.4;
double l = 0.5;
long int myiter = 0;
long int nniters;
double eps_x = pow(10, -8);
double T_start = 293.;
double T_finish = 7.289 * T_start;
double mol_weight = 21.0;
double gamma = 1.17;
double q = -43.28 * R * T_start / mol_weight;

typedef struct {
    realtype* x;
    realtype* T;
    realtype* Y_H2O;
    int Nx;
    int N_m;
    int NEQ;
    int N_centr;
    int n_tout;
    realtype Tl;
    realtype Tr;
    realtype M;
    realtype T_center;
    ofstream* outNevyaz;
    int kp;
    void* mykmem;
} *UserData;


/* Accessor macro */
#define Ith(v,i)    NV_Ith_S(v,i)
#define IJth(A,i,j) SM_ELEMENT_D(A,i,j)

/* Functions Called by the KINSOL Solver */
static int func(N_Vector u, N_Vector f, void* user_data);
static int func_Y(N_Vector u, N_Vector f, void* user_data);
static int resrob(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data);

void ExportToArray(vector<double>& T_vect, vector<double>& Y_vect, double& M, UserData data, N_Vector yy, int N_x)
{
    //cout << "data->Tl; = " << data->Tl << "\n";
    T_vect[0] = data->Tl;
    for (int i = 0; i < N_x - 1; i++)
    {
       T_vect[i] = Ith(yy, i);
    }
    T_vect[N_x - 1] = data->Tr;
    T_vect[N_x] = Ith(yy, N_x - 1);
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

void ExportToArray(double* T_vect, double* Y_vect, double& M, UserData data, N_Vector yy, int N_x)
{
    //cout << "data->Tl; = " << data->Tl << "\n";
    T_vect[0] = data->Tl;
    for (int i = 0; i < N_x - 1; i++)
    {
        T_vect[i] = Ith(yy, i);
    }
    T_vect[N_x - 1] = data->Tr;
    T_vect[N_x] = Ith(yy, N_x - 1);
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

void get_Y(double Y_H2O, double& Y_H2, double& Y_O2, double Y_N2)
{
    Y_H2 = (1 - Y_H2O - Y_N2) * 1. / 9.;
    Y_O2 = (1 - Y_H2O - Y_N2) * 8. / 9.;
}

double get_W(double Y_H2O, double Y_H2, double Y_O2, double Y_N2)
{
    return mol_weight;
}

double Cp_all(double T, double Y_H2O)
{

    //return R * gamma / (gamma - 1.) / mol_weight;
    return 2.737;
}

double Lambda(double T, double Y_H2O)
{
    return 2.9 * pow(10, -5) * pow(T, 0.7) * Cp_all(T, Y_H2O);
}

void Get_mole_fr(double& X_H2O, double& X_H2, double& X_O2, double& X_N2, double Y_H2O, double T)
{
    double Y_H2, Y_O2;
    get_Y(Y_H2O, Y_H2, Y_O2, Y_N2);
    double W = 0;

    W = get_W(Y_H2O, Y_H2, Y_O2, Y_N2);
    X_H2O = Y_H2O * W / phyc.mol_weight[6];
    X_H2 = Y_H2 * W / phyc.mol_weight[0];
    X_O2 = Y_O2 * W / phyc.mol_weight[2];
    X_N2 = Y_N2 * W / phyc.mol_weight[8];
}

double V(double Y_H2O, double Y_H2Or, double T, double Tr, realtype* x_cells, const int i)
{
    double h = x_cells[i + 1] - x_cells[i];
    double X_H2O, X_H2, X_O2, X_N2;
    Get_mole_fr(X_H2O, X_H2, X_O2, X_N2, Y_H2O, T);
    /*cout << "X_H2O = " << X_H2O << "\n";
    cout << "X_H2 = " << X_H2 << "\n";
    cout << "X_O2 = " << X_O2 << "\n";
    cout << "X_N2 = " << X_N2 << "\n";*/
    double X_H2Or, X_H2r, X_O2r, X_N2r;
    Get_mole_fr(X_H2Or, X_H2r, X_O2r, X_N2r, Y_H2Or, Tr);

    double Y_H2, Y_O2;
    get_Y((Y_H2O + Y_H2Or) / 2., Y_H2, Y_O2, Y_N2);
    double W = get_W((Y_H2O + Y_H2Or) / 2., Y_H2, Y_O2, Y_N2);
    double rho = P * W / R / T;
    double K = 1. - Y_N2;
    double Y = Y_H2O / K;
    double Yr = Y_H2Or / K;
    double D = Lambda((Tr + T) / 2., Y_H2O) / rho / Cp_all((T + Tr) / 2., Y_H2O);
    if ((Y_H2Or + Y_H2O) < eps_x) return 0;
    else return -D * (Yr - Y) / h / (Y + Yr) * 2.;
}

double rhoYkVk(double Y_H2O, double Y_H2Or, double T, double Tr, realtype* x_cells, const int i)
{
    double h = x_cells[i + 1] - x_cells[i];
    double X_H2O, X_H2, X_O2, X_N2;
    Get_mole_fr(X_H2O, X_H2, X_O2, X_N2, Y_H2O, T);
    double X_H2Or, X_H2r, X_O2r, X_N2r;
    Get_mole_fr(X_H2Or, X_H2r, X_O2r, X_N2r, Y_H2Or, Tr);
    double K = 1. - Y_N2;
    double Yr = Y_H2Or / K;
    double Y_H2, Y_O2;
    get_Y((Y_H2O + Y_H2Or) / 2., Y_H2, Y_O2, Y_N2);
    double W = get_W((Y_H2O + Y_H2Or) / 2., Y_H2, Y_O2, Y_N2);
    double rho = P * W / R / (T + Tr) * 2.;
    double D = Lambda((Tr + T) / 2., (Y_H2O + Y_H2Or) / 2.) / rho / Cp_all((T + Tr) / 2., (Y_H2O + Y_H2Or) / 2.);
    return -D * rho * (Y_H2Or - Y_H2O) / h;
}

int InitialData(int& Nx, vector<double>& x_vect, vector<double>& T_vect, vector<double>& Y_vect, double& M)
{
    double h = l / (Nx - 1);
    double x_start = koeff_l * l - l / 10.;
    double x_finish = koeff_l * l + l / 10.;
    int dN = (x_finish - x_start) / h;
    double j = 0;
    M = 1000 * 0.000871523;
    cout << "M = " << M << "\n";
    double W;

    for (int i = 0; i < Nx; i++) {
        x_vect[i] = h * i;
    }

    T_vect[0] = T_start;
    j = 0;
    for (int i = 0; i < Nx; i++) {
        if (x_vect[i] <= x_start)
        {
            Y_vect[i] = 0;
            T_vect[i] = T_start;
        }
        else if (x_vect[i] >= x_finish)
        {
            Y_vect[i] = (1. - Y_N2);
            T_vect[i] = T_finish;
        }
        else {
            Y_vect[i] = (1. - Y_N2) / dN * j;
            T_vect[i] = (T_finish - T_start) / dN * j + T_start;
            j++;
        }
    }
    /* for (int i = 0; i < Nx; i++) {
        T_vect[i] = 1111.31 * tanh((x_vect[i] - koeff_l * l) / 0.064) + 1111.31 + 293.;
    }*/

    for (int i = 0; i < Nx - 1; i++) {
        if (x_vect[i] <= koeff_l * l && x_vect[i + 1] > koeff_l * l)
            return i;
    }
}

void Add_elem(vector<double>& T, vector<double>& Y, vector<double>& x, int& N_x, int& N_center, double b)
{
    int j_t = 1;
    double T_max = 0, T_min = T[0];

    for (int i = 0; i < N_x; i++)
    {
        if (T[i] > T_max) T_max = T[i];
        if (T[i] < T_min) T_min = T[i];
    }

    while (j_t < N_x - 2)
    {
        if (fabs(T[j_t] - T[j_t - 1]) > b * (T_max - T_min))
        {
            T.insert(T.begin() + j_t, (T[j_t] + T[j_t - 1]) / 2.);
            Y.insert(Y.begin() + j_t, (Y[j_t] + Y[j_t - 1]) / 2.);
            x.insert(x.begin() + j_t, (x[j_t] + x[j_t - 1]) / 2.);
            N_x++;
            j_t++;
        }
        j_t++;
        //cout << "j_t = " << j_t << "\n";
    }
    T.insert(T.begin(), T[0]);
    Y.insert(Y.begin(), Y[0]);
    x.insert(x.begin(), -x[1]);
    N_x++;

    for (int k = 0; k < 3; k++)
    {
        T.insert(T.begin(), T[1]);
        Y.insert(Y.begin(), Y[1]);
        x.insert(x.begin(), 1.6 * x[0]);
        N_x++;
    }

    for (int i = 0; i < N_x - 1; i++) {
        if (x[i] <= koeff_l * l && x[i + 1] > koeff_l * l)
            N_center = i;
    }
}

double F_right(double T_left, double T_center, double T_right, double M, double Y_H2O, double Y_H2Or, realtype* x_cells, const int i)
{
    double h_left = x_cells[i] - x_cells[i - 1];
    double h = x_cells[i + 1] - x_cells[i];
    double Cp = Cp_all(T_center, Y_H2O);

    //cout << "Cp = " << Cp << "\n";
    //cout << "lambda = " << Lambda((T_right + T_center) / 2., Y_H2O) << "\n";
    double Y_H2, Y_O2;
    get_Y(Y_H2O, Y_H2, Y_O2, Y_N2);
    double W = get_W(Y_H2O, Y_H2, Y_O2, Y_N2);
    //cout << "W = " << W << "\n";
    double rho = P * W / R / T_center;
    //cout << "rho = " << rho << "\n";
    double K = 1. - Y_N2;
    double Y = 1. - Y_H2O / K;
    double dTdx = (h_left / h / (h + h_left) * T_right + (h - h_left) / h / h_left * T_center - h / h_left / (h + h_left) * T_left);
    //cout << "Hi = " << get_Hi(6, T_center) * pow(10, 3) << endl;
    double w_dot = A * rho * rho * Y * exp(-Ea / T_center);
    /*cout << "Lambda = " << -(2. / (h + h_left)) *
        (Lambda((T_right + T_center) / 2., Y_H2O) * (T_right - T_center) / h
            - Lambda((T_center + T_left) / 2., Y_H2O) * (T_center - T_left) / h_left) << "\n";
    cout << "Cp = " << Cp * M * dTdx << "\n";
    cout << "dT/dx = " << dTdx << "\n";
    cout << "M = " << M << "\n";*/
    //cout << " v = " << rho * Y_H2O / K * V(Y_H2O, Y_H2Or, T_center, T_right, x_cells, i) * Cp * dTdx << "\n";
    return -(2. / (h + h_left)) *
        (Lambda((T_right + T_center) / 2., Y_H2O) * (T_right - T_center) / h
            - Lambda((T_center + T_left) / 2., Y_H2O) * (T_center - T_left) / h_left)
        + Cp * M * dTdx
        + w_dot * q;
    //+ rho * Y_H2O / K  * V(Y_H2O, Y_H2Or, T_center, T_right, x_cells, i) * Cp * dTdx;
/*return -(2. / (h + h_left)) *
    (Lambda((T_right + T_center) / 2., Y_H2O) * (T_right - T_center) / h
        - Lambda((T_center + T_left) / 2., Y_H2O) * (T_center - T_left) / h_left)
    + Cp * M * dTdx
    + w_dot * pow(10, 3) * (-5.); */
}

double F_rightY(double T_left, double T_center, double T_right, double M, double Y_H2O_left, double Y_H2O_center, double Y_H2O_right, realtype* x_cells, const int i)
{
    double h_left = x_cells[i] - x_cells[i - 1];
    double h = x_cells[i + 1] - x_cells[i];
    double Y_H2, Y_O2;
    get_Y(Y_H2O_center, Y_H2, Y_O2, Y_N2);
    double W = get_W(Y_H2O_center, Y_H2, Y_O2, Y_N2);
    double rho = P * W / R / T_center;
    //cout << "rho = " << rho << "\n";
    double K = 1 - Y_N2;
    double Y = 1 - Y_H2O_center / K;
    double w_dot = K * A * rho * rho * Y * exp(-Ea / T_center);
    double h_right_half = (x_cells[i + 1] + x_cells[i]) / 2.;
    double h_left_half = (x_cells[i] + x_cells[i - 1]) / 2.;
    /*cout << "M = " << M * (h_left / h / (h + h_left) * Y_H2O_right + (h - h_left) / h / h_left * Y_H2O_center
       - h / h_left / (h + h_left) * Y_H2O_left) << "\n";
    cout << "wdot = " << w_dot << "\n";*/
    return M * (h_left / h / (h + h_left) * Y_H2O_right + (h - h_left) / h / h_left * Y_H2O_center - h / h_left / (h + h_left) * Y_H2O_left) - w_dot
        + (rhoYkVk(Y_H2O_center, Y_H2O_right, T_center, T_right, x_cells, i) - rhoYkVk(Y_H2O_left, Y_H2O_center, T_left, T_center, x_cells, i - 1))
        / (h_right_half - h_left_half);
}

int Integrate(int N_x, vector<double>& x_vect, vector<double>& T_vect, vector<double>& Y_vect, double& M, int N_center)
{
    SUNContext sunctx;
    UserData data;
    realtype fnormtol, scsteptol;
    N_Vector res_vect, s, s2, c;
    int glstr, mset, retval;
    void* kmem;
    SUNMatrix J;
    SUNLinearSolver LS;
    int NEQ_T = N_x - 2;
    int NEQ_Y = N_x - 2;
    int NEQ = NEQ_T + NEQ_Y;

    res_vect = NULL;
    s = c = s2 = NULL;
    kmem = NULL;
    J = NULL;
    LS = NULL;
    data = NULL;

    /* Create the SUNDIALS context that all SUNDIALS objects require */
    retval = SUNContext_Create(NULL, &sunctx);
    if (check_retval(&retval, "SUNContext_Create", 1)) return(1);

    /* User data */

    data = (UserData)malloc(sizeof * data);

    /* Create serial vectors of length NEQ */

    res_vect = N_VNew_Serial(NEQ, sunctx);
    if (check_retval((void*)res_vect, "N_VNew_Serial", 0)) return(1);

    s = N_VNew_Serial(NEQ, sunctx);
    if (check_retval((void*)s, "N_VNew_Serial", 0)) return(1);

    s2 = N_VNew_Serial(NEQ, sunctx);
    if (check_retval((void*)s2, "N_VNew_Serial", 0)) return(1);

    c = N_VNew_Serial(NEQ, sunctx);
    if (check_retval((void*)c, "N_VNew_Serial", 0)) return(1);

    //SetInitialGuess1(res_vect, data, NEQ);

    N_VConst(ONE, s); /* no scaling */
    N_VConst(ONE, s2);

    data->Nx = N_x;
    data->x = new realtype[N_x];

    data->Y_H2O = new realtype[N_x];
    data->T = new realtype[N_x];
    data->NEQ = NEQ;
    data->Tl = T_vect[0];
    data->T_center = T_vect[N_center];
    cout << "T_center = " << data->T_center << "\n";
    data->N_centr = N_center;
    int j = 0;
    for (int i = 0; i < N_x; i++) {
        data->x[i] = x_vect[i];
        data->T[i] = T_vect[i];
        data->Y_H2O[i] = Y_vect[i];
        //cout << i << " = " << Ith(res_vect, i + 1) << endl;
    }
    j = 1;
    //cout << NEQ << " = " << Ith(res_vect, NEQ + 1) << endl;
    for (int i = 1; i < NEQ_T; i++) {
        Ith(c, i) = 0.0;   /* no constraint on x1 */
        if (j == N_center)
        {
            j++;
        }
        Ith(res_vect, i) = T_vect[j];
        j++;
    }
    Ith(c, NEQ_T) = ONE;
    Ith(res_vect, NEQ_T) = M;
    data->N_m = NEQ_T - 1;
    for (int i = NEQ_T + 1; i <= NEQ; i++) {
        Ith(c, i) = 0.0;   /* constraint on x1 */
        Ith(s, i) = 2500. / 0.25;
        Ith(res_vect, i) = Y_vect[i - NEQ_T];
        //cout << "Yvect " << i - NEQ_T << " =  " << Y_vect[i - NEQ_T] - Y_max << endl;
    }
    fnormtol = FTOL; scsteptol = STOL;


    kmem = KINCreate(sunctx);
    if (check_retval((void*)kmem, "KINCreate", 0)) return(1);

    retval = KINSetUserData(kmem, data);
    if (check_retval(&retval, "KINSetUserData", 1)) return(1);
    retval = KINSetConstraints(kmem, c);
    if (check_retval(&retval, "KINSetConstraints", 1)) return(1);
    retval = KINSetFuncNormTol(kmem, fnormtol);
    if (check_retval(&retval, "KINSetFuncNormTol", 1)) return(1);
    retval = KINSetScaledStepTol(kmem, scsteptol);
    if (check_retval(&retval, "KINSetScaledStepTol", 1)) return(1);
    retval = KINSetMaxSetupCalls(kmem, 1);
    retval = KINSetNumMaxIters(kmem, 8000);
    retval = KINInit(kmem, func, res_vect);
    if (check_retval(&retval, "KINInit", 1)) return(1);


    /* Create dense SUNMatrix */
    J = SUNDenseMatrix(NEQ, NEQ, sunctx);
    if (check_retval((void*)J, "SUNDenseMatrix", 0)) return(1);

    /* Create dense SUNLinearSolver object */
    LS = SUNLinSol_Dense(res_vect, J, sunctx);
    if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) return(1);

    /* Attach the matrix and linear solver to KINSOL */
    retval = KINSetLinearSolver(kmem, LS, J);
    if (check_retval(&retval, "KINSetLinearSolver", 1)) return(1);

    glstr = 0;
    mset = 500;

    data->mykmem = kmem;
    retval = KINSol(kmem, res_vect, glstr, s, s2);
    if (check_retval(&retval, "KINSol", 1)) return(1);
    cout << "retval2 = " << retval << "\n";

    int myNx = data->Nx;
    int myNeq = data->NEQ;
    int myNm = data->N_m;
    T_vect[0] = data->Tl;
    ExportToArray(T_vect, Y_vect, M, data, res_vect, N_x);

  
    /* Free memory */
    printf("\nFinal statsistics:\n");
    retval = KINPrintAllStats(kmem, stdout, SUN_OUTPUTFORMAT_TABLE);
    N_VDestroy(res_vect);
    N_VDestroy(s);
    N_VDestroy(c);
    KINFree(&kmem);
    SUNLinSolFree(LS);
    SUNMatDestroy(J);
    free(data);
    SUNContext_Free(&sunctx);
    return 0;
}

int Integrate_Y(int N_x, vector<double>& x_vect, vector<double>& T_vect, vector<double>& Y_vect, double& M, int N_center, int diff)
{
    SUNContext sunctx;
    UserData data;
    realtype fnormtol, scsteptol;
    N_Vector res_vect, s, c;
    int glstr, mset, retval;
    void* kmem;
    SUNMatrix J;
    SUNLinearSolver LS;
    //int NEQ_T = N_x - 2;
    int NEQ_Y = N_x - 2;
    int NEQ = 1 + NEQ_Y;

    res_vect = NULL;
    s = c = NULL;
    kmem = NULL;
    J = NULL;
    LS = NULL;
    data = NULL;

    /* Create the SUNDIALS context that all SUNDIALS objects require */
    retval = SUNContext_Create(NULL, &sunctx);
    if (check_retval(&retval, "SUNContext_Create", 1)) return(1);

    /* User data */

    data = (UserData)malloc(sizeof * data);

    /* Create serial vectors of length NEQ */

    res_vect = N_VNew_Serial(NEQ, sunctx);
    if (check_retval((void*)res_vect, "N_VNew_Serial", 0)) return(1);

    s = N_VNew_Serial(NEQ, sunctx);
    if (check_retval((void*)s, "N_VNew_Serial", 0)) return(1);

    c = N_VNew_Serial(NEQ, sunctx);
    if (check_retval((void*)c, "N_VNew_Serial", 0)) return(1);

    N_VConst(ONE, s); /* no scaling */

    data->Nx = N_x;
    data->x = new realtype[N_x];

    data->Y_H2O = new realtype[N_x];
    data->T = new realtype[N_x];
    data->NEQ = NEQ;
    data->Tl = T_vect[0];
    data->T_center = T_vect[N_center];
    cout << "T_center = " << data->T_center << "\n";
    data->N_centr = N_center;
    int j = 0;
    for (int i = 0; i < N_x; i++) {
        data->x[i] = x_vect[i];
        data->T[i] = T_vect[i];
        data->Y_H2O[i] = Y_vect[i];
        //cout << i << " = " << Ith(res_vect, i + 1) << endl;
    }

    Ith(c, 1) = ONE;
    Ith(res_vect, 1) = M;
    data->N_m = 1 - 1;

    for (int i = 2; i <= NEQ; i++) {
        Ith(c, i) = 0.0;   /* no constraint on x1 */
        Ith(res_vect, i) = Y_vect[i - 1];
        //cout << "Yvect " << i - N_x + 1 << " =  " << Y_vect[i - N_x + 1] << endl;
    }
    fnormtol = FTOL; scsteptol = STOL;


    kmem = KINCreate(sunctx);
    if (check_retval((void*)kmem, "KINCreate", 0)) return(1);

    retval = KINSetUserData(kmem, data);
    if (check_retval(&retval, "KINSetUserData", 1)) return(1);
    retval = KINSetConstraints(kmem, c);
    if (check_retval(&retval, "KINSetConstraints", 1)) return(1);
    retval = KINSetFuncNormTol(kmem, fnormtol);
    if (check_retval(&retval, "KINSetFuncNormTol", 1)) return(1);
    retval = KINSetScaledStepTol(kmem, scsteptol);
    if (check_retval(&retval, "KINSetScaledStepTol", 1)) return(1);

    cout << "func_Y " << "\n";
    retval = KINInit(kmem, func_Y, res_vect);
    if (check_retval(&retval, "KINInit", 1)) return(1);

    /* Create dense SUNMatrix */
    J = SUNDenseMatrix(NEQ, NEQ, sunctx);
    if (check_retval((void*)J, "SUNDenseMatrix", 0)) return(1);

    /* Create dense SUNLinearSolver object */
    LS = SUNLinSol_Dense(res_vect, J, sunctx);
    if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) return(1);

    /* Attach the matrix and linear solver to KINSOL */
    retval = KINSetLinearSolver(kmem, LS, J);
    if (check_retval(&retval, "KINSetLinearSolver", 1)) return(1);

    glstr = 0;
    mset = 500;

    data->mykmem = kmem;
    retval = KINSol(kmem, res_vect, glstr, s, s);
    if (check_retval(&retval, "KINSol", 1)) return(1);
    //PrintFinalStats(kmem);
    cout << "retval = " << retval << "\n";

    int myNx = data->Nx;
    int myNeq = data->NEQ;
    int myNm = data->N_m;

    T_vect[0] = data->Tl;
    //cout << "MyNx = " << myNx << "\n";
    //cout << "MyNm = " << myNm << "\n";
    M = Ith(res_vect, myNm + 1);
    Y_vect[0] = 0.;
    cout << "M = " << M << "\n";
    for (int i = myNm + 1; i < N_x - 1; i++)
    {
        Y_vect[i] = Ith(res_vect, i + 1);
        //cout << "Y_vect  " << i << " =  " << Y_vect[i] << endl;
    }
    Y_vect[N_x - 1] = Y_vect[N_x - 2];
    /* Free memory */
    printf("\nFinal statsistics:\n");
    retval = KINPrintAllStats(kmem, stdout, SUN_OUTPUTFORMAT_TABLE);
    N_VDestroy(res_vect);
    N_VDestroy(s);
    N_VDestroy(c);
    KINFree(&kmem);
    SUNLinSolFree(LS);
    SUNMatDestroy(J);
    free(data);
    SUNContext_Free(&sunctx);
    return 0;
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
//Formula setting for Heat Capacity and it's derivative, J/(kg*K)
double Cp(double T, const double a)
{
    const double Cp_water = 4180.;
    const double Cp_steam = 2000.;
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
    {
        if (T <= T_boiling)
            return Cp_water;
        else
            return Cp_steam;
    }
    else
        return (AA + B * T_forCp + C * pow(T_forCp, 2.) + D * pow(T_forCp, 3.) + EE * pow(T_forCp, -2.)) / W;
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
    const double rho_water = 980.;
    const double rho_steam = 0.4;
    const double T_boiling = 373;
    const double R = 8.314462;
    const double p = pow(10, 5);
    const double W = 18 * pow(10, -3);
    if (a < pow(10, -8))
    {
        if (T <= T_boiling)
            return rho_water;
        else
            return rho_steam;
    }
    else
        if (T <= T_boiling)
            return rho_water;
        else
            return p * W / (R * T);
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
    /*
    const double lambda_water = 0.56;
    const double lambda_steam = 0.05;
    const double T_boiling = 373;
    if (a < pow(10, -8))
    {
        if (T <= T_boiling)
            return lambda_water;
        else
            return lambda_steam;
    }
    else
    {
        if (T <= T_boiling && flag_phase == 0)
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
    */
    return 0.05;
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
void ArraysParameters(const double a, const int N, vector <double>& r, vector <double>& T_next,
    const int s, vector <double>& ACp, vector <double>& ADensity, vector <double>& ALambda, int n,
    double inter, double Ti, double p)
{
    //Define Parameters
    ACp.resize(r.size());
    //ADfCp.resize(r.size());
    ADensity.resize(r.size());
    //ADfDensity.resize(r.size());
    ALambda.resize(r.size());
    //ADfLambda.resize(r.size());

    for (int i = 0; i < s + 1; i++) {
        //cout << r[i] << endl;
        ACp[i] = Cp(T_next[i], a);
        //ADfCp[i] = DfCp(T_next[i], a);
        ADensity[i] = Density(T_next[i], a);
        //ADfDensity[i] = DfDensity(T_next[i], a);
        ALambda[i] = Lyambda(T_next[i], a, 0);
        //ADfLambda[i] = DfLambda(T_next[i], a, 0);
    }
    for (int i = s + 1; i < r.size(); i++) {
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
    for (j; j < s + 1; j++) {
        Parameters << r[j] << " " << T_next[j] << " " << ACp[j] << " " << ADensity[j] << " " << ALambda[j] << "\n";
    }
    //Interface
    //Left side
    Parameters << inter << " " << Ti << " " << Cp(T_next[N], a) << " " << Density(T_next[N], a) << " "
        << Lyambda(T_next[N], a, 0) << "\n";
    //Right side
    Parameters << inter << " " << Ti << " " << Cp(T_next[N], a) << " " << Density(T_next[N], a) << " "
        << Lyambda(T_next[N], a, 1) << "\n";
    for (j; j < N; j++) {
        Parameters << r[j] << " " << T_next[j] << " " << ACp[j] << " " << ADensity[j] << " " << ALambda[j] << "\n";
    }
    Parameters.close();
}

void ArraysParametersCopyForArrays(const double a, const int N, realtype* r, realtype* T_next,
    const int s, vector <double>& ACp, vector <double>& ADensity, vector <double>& ALambda, int n,
    double inter, double Ti, double p)
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
        ALambda[i] = Lyambda(T_next[i], a, 0);
        //ADfLambda[i] = DfLambda(T_next[i], a, 0);
    }
    for (int i = s + 1; i < N + 1; i++) {
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
    for (j; j < s + 1; j++) {
        Parameters << r[j] << " " << T_next[j] << " " << ACp[j] << " " << ADensity[j] << " " << ALambda[j] << "\n";
    }
    //Interface
    //Left side
    Parameters << inter << " " << Ti << " " << Cp(T_next[N], a) << " " << Density(T_next[N], a) << " "
        << Lyambda(T_next[N], a, 0) << "\n";
    //Right side
    Parameters << inter << " " << Ti << " " << Cp(T_next[N], a) << " " << Density(T_next[N], a) << " "
        << Lyambda(T_next[N], a, 1) << "\n";
    for (j; j < N; j++) {
        Parameters << r[j] << " " << T_next[j] << " " << ACp[j] << " " << ADensity[j] << " " << ALambda[j] << "\n";
    }
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


double FDrop(double T_left, double T_center, double T_right, vector<double> r, vector<double> ALambda, int j)
{
    double rj_diff = r[j + 1] - r[j - 1];
    double rj_avg = pow((r[j + 1] + r[j]) / 2., 2.);
    double rj_minus_avg = pow((r[j] + r[j - 1]) / 2., 2.);
    double Lambda_j_half = 0.5 * (ALambda[j + 1] + ALambda[j]);
    double Lambda_j_minus_half = 0.5 * (ALambda[j] + ALambda[j - 1]);
    return 2. / rj_diff * (rj_avg * Lambda_j_half * (T_right - T_center) / (r[j + 1] - r[j])
        - rj_minus_avg * Lambda_j_minus_half * (T_center - T_left) / (r[j] - r[j - 1]));
}

double FSteam(double T_left, double T_center, double T_right, vector<double> r, vector<double> ALambda, vector<double>ACp,
    double dM, int j)
{
    double rj_diff = r[j + 1] - r[j - 1];
    double rj_avg = pow((r[j + 1] + r[j]) / 2., 2.);
    double rj_minus_avg = pow((r[j] + r[j - 1]) / 2., 2.);
    double Lambda_j_half = 0.5 * (ALambda[j + 1] + ALambda[j]);
    double Lambda_j_minus_half = 0.5 * (ALambda[j] + ALambda[j - 1]);
    return 2. / rj_diff * (rj_avg * Lambda_j_half * (T_right - T_center) / (r[j + 1] - r[j])
        - rj_minus_avg * Lambda_j_minus_half * (T_center - T_left) / (r[j] - r[j - 1]))
        - ACp[j] * dM * (T_center - T_left) / (r[j] - r[j - 1]);
}
double FDrop4Arrays(double T_left, double T_center, double T_right, realtype* r, vector<double> ALambda, int j)
{
    double rj_diff = r[j + 1] - r[j - 1];
    double rj_avg = pow((r[j + 1] + r[j]) / 2., 2.);
    double rj_minus_avg = pow((r[j] + r[j - 1]) / 2., 2.);
    double Lambda_j_half = 0.5 * (ALambda[j + 1] + ALambda[j]);
    double Lambda_j_minus_half = 0.5 * (ALambda[j] + ALambda[j - 1]);
    return 2. / rj_diff * (rj_avg * Lambda_j_half * (T_right - T_center) / (r[j + 1] - r[j])
        - rj_minus_avg * Lambda_j_minus_half * (T_center - T_left) / (r[j] - r[j - 1]));
}

double FSteam4Arrays(double T_left, double T_center, double T_right, realtype* r, vector<double> ALambda, vector<double>ACp,
    double dM, int j)
{
    double rj_diff = r[j + 1] - r[j - 1];
    double rj_avg = pow((r[j + 1] + r[j]) / 2., 2.);
    double rj_minus_avg = pow((r[j] + r[j - 1]) / 2., 2.);
    double Lambda_j_half = 0.5 * (ALambda[j + 1] + ALambda[j]);
    double Lambda_j_minus_half = 0.5 * (ALambda[j] + ALambda[j - 1]);
    return 2. / rj_diff * (rj_avg * Lambda_j_half * (T_right - T_center) / (r[j + 1] - r[j])
        - rj_minus_avg * Lambda_j_minus_half * (T_center - T_left) / (r[j] - r[j - 1]))
        - ACp[j] * dM * (T_center - T_left) / (r[j] - r[j - 1]);
}
void FRight(vector<double>& f, vector<double> T_cur, vector<double> T_next, vector<double> r, vector <double> ACp, vector <double> ADensity, vector <double> ALambda,
    const int N, const double a, const double b, const double c, const double d, const double dt, double dM, int s, double inter,
    const double Ti, double p, const double L_d, double& mod_nevyaz, int flag_evap)
{
    double F;
    /*
    //Nevyazka on the left border
    int l_f_num = 0;
    F = 1. / (r[1] - r[0]) * pow((r[l_f_num + 1] + r[l_f_num]) / 2., 2.) * 0.5 * (ALambda[l_f_num + 1] + ALambda[l_f_num])
        * (T_next[l_f_num + 1] - T_next[l_f_num]) / (r[l_f_num + 1] - r[l_f_num])
        - ACp[l_f_num] * ADensity[l_f_num] * pow(r[l_f_num + 1], 2.) / 24. * (T_next[l_f_num] - T_cur[l_f_num]) / dt;
    f[l_f_num] = -F;
    //cout << "f[1]" << f[1] << "\n";
    */
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

void FRightTest(vector <double>& f, vector <double>& T_vect, vector <double>& r, vector <double>& ACp, vector <double>& ADensity,
    vector <double>& ALambda, const int n, const int N, const double a, double dM, int s, double inter,
    const double Ti, double p, const double L_d, double& q_i_right, double& q_i_left)
{
    int N_minus = N - 1;
    double F;
    /*
    for (int i = 0; i <= N + 1; i++) {
        cout << "r[" << i << "]= " << r[i] << "\n";
        cout << "T[" << i << "]= " << T_vect[i] << "\n";
    }
    */

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
        F = FDrop(T_vect[j - 1], T_vect[j], T_vect[j + 1], r, ALambda, j);
        f[j] = F / (ACp[j] * ADensity[j] * pow(r[j], 2.));
        //cout << j << " f1_" << j << " " << f[j] << "\n";
    }
    //Nevyazka in Gaze
    for (int j = s + 1; j < N_minus; j++)
    {
        //cout << "(T_vect[j] - T[j][n])" << (T_vect[j] - T[j][n]) << endl;
        F = FSteam(T_vect[j - 1], T_vect[j], T_vect[j + 1], r, ALambda, ACp, dM, j);
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
void FRightTestCopyForArrays(vector <double>& f, realtype* T_vect, realtype* r, vector <double>& ACp, vector <double>& ADensity,
    vector <double>& ALambda, const int n, const int N, const double a, double dM, int s, double inter,
    const double Ti, double p, const double L_d, double& q_i_right, double& q_i_left)
{
    int N_minus = N - 1;
    double F;
    /*
    for (int i = 0; i <= N + 1; i++) {
        cout << "r[" << i << "]= " << r[i] << "\n";
        cout << "T[" << i << "]= " << T_vect[i] << "\n";
    }
    */

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
    vector<double>& T_cur, double h, const double q, const double h_min, const int const_params,
    const int Nd, const int N_uni, const int N_uni_near_center)
{
    ofstream OutX;
    OutX.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/Data_new/X_grid.dat");
    OutX << "TITLE=\"" << "Graphics" << "\"" << endl;
    OutX << R"(VARIABLES= "j", "hj, m" )" << endl;
    //Setting the value on the boundaries and initial distribution
    r[0] = x_l;
    T_cur[0] = T_l;
    cout << "h" << h << " r" << 0 << " " << r[0] << " T= " << T_cur[0] << "\n";
    h = h_min;
    //cout << "r1" << " " << r[1] << "\n";
    /*
    double T_step = 12;
    int j = 1;
    for (j; j < r.size(); j++)
    {
        r[j] = r[j - 1] + h;
        //h = h * q;
        T_cur[j] = T_cur[j - 1] + T_step;
        //cout << "T" << j << " " << T_cur[j] << endl;
        OutX << j << " " << h << " " << "\n";
        cout << "h" << h << " r" << j << " " << r[j] << " T= " << T_cur[j] << "\n";
    }
    OutX.close();
    */
    //Uniform for droplet and non-uniform for gas 1-dimensional grid 
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
        for (j; j < r.size(); j++)
        {
            r[j] = r[j - 1] + h;
            //h = h * q;
            T_cur[j] = T_r;
            //cout << "T" << j << " " << T_cur[j] << endl;
            OutX << j << " " << h << " " << "\n";
            cout << "h" << h << " r" << j << " " << r[j] << " T= " << T_cur[j] << "\n";
        }
    }
    OutX.close();
}


int Integrate_IDA(int N_x, vector<double>& x_vect, vector<double>& T_vect, vector<double>& Y_vect,
    double& M, int N_center, double tout1, int call, int number,int print_value, int cons_flag,
    int s)
{
    void* mem;
    N_Vector yy, yp, avtol, cons;
    realtype rtol, * yval, * ypval, * atval;
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
    int NEQ_T = N_x - 1;
    int NEQ = NEQ_T + 1;
    N_center = s;
    int k = 0;
    //Дальше не трогаем

    data->kp = k;
    data->outNevyaz = new ofstream;
    data->Nx = N_x;
    data->x = new realtype[N_x];
    data->T = new realtype[N_x + 1];
    data->NEQ = NEQ;
    data->Tl = T_vect[0];
    data->Tr = T_vect[N_x - 1];
    data->T_center = T_vect[N_center];
    data->N_centr = N_center;
    data->M = M;

    int j = 0;
    for (int i = 0; i < N_x; i++) {
        data->x[i] = x_vect[i];
        data->T[i] = T_vect[i];
        //cout << i << " = " << Ith(res_vect, i + 1) << endl;
    }
    //Put Ti in T-massive as T[N_x] in struct data
    data->T[N_x] = T_vect[N_x];
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
    rtol = RCONST(1.0e-6);
    atval = N_VGetArrayPointer(avtol);
    /* Integration limits */
    t0 = ZERO;

    //заполняю идовский массив
    //MAXIM added
    double p = 1.5;
    double inter = x_vect[s] + (p - 1) * (x_vect[s + 1] - x_vect[s]);
    //double q = 20 * pow(10, 3);                  //W / m^2
    double L_d = 2258.2 * pow(10, 3.);              //J / kg 
    //const determination of dM
    //dM = 4 * PI * pow(inter, 2.0) * q / L_d;     // kg / s
    //dM equals zero because it is only warming time
    double dM = 0.;                                       // kg / s
    double a = 1.;
        //pow(10, -10);
    vector <double> f(N_x);
    vector <double> ACp, ADfCp, ADensity, ADfDensity, ALambda, ADfLambda;
    ArraysParameters(a, N_x, x_vect, T_vect, s, ACp, ADensity, ALambda, 0, inter, T_vect[N_x], p);
    double q_i_g = 0;
    double q_i_d = 0;
    FRightTest(f, T_vect, x_vect, ACp, ADensity, ALambda, 0, N_x, a, dM, s, inter, T_vect[N_x], p, L_d, q_i_g, q_i_d);
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
    *(data->outNevyaz) << R"(VARIABLES= "i", "F",  "q_g, W/m^2",  "q_d, W/m^2")" << endl;
    //Calculate modul of Nevyazka
    double mod_nevyaz = 0;
    for (int i = 0; i < N_x; i++)
        mod_nevyaz += pow(f[i], 2.0);
    mod_nevyaz = pow(mod_nevyaz, 0.5);
    *(data->outNevyaz) << k << " " << mod_nevyaz << " " << abs(q_i_g / pow(inter, 2.)) << " "
        << abs(q_i_d / pow(inter, 2.)) << " " << "\n";

    for (int i = 0; i < NEQ_T; i++) {
        Ith(avtol, i) = RCONST(1.0e-8);
        Ith(yy, i) = T_vect[i];
        //cout << "yy " << i << " = " << Ith(yy, i) << "\n";
    }
    Ith(avtol, NEQ_T) = RCONST(1.0e-8);
    Ith(yy, NEQ_T) = T_vect[N_x];

    for (int i = 0; i < NEQ_T; i++) {
        Ith(yp, i) = f[i];
    }
    Ith(yp, NEQ_T) = ZERO;
    cout << "Ith(yp, NEQ_T)" << Ith(yp, NEQ_T) << "\n";
    //Дальше не трогаем
    /*
    //j = 1;
    //cout << NEQ << " = " << Ith(res_vect, NEQ + 1) << endl;
    for (int i = 1; i < NEQ_T; i++) {
        Ith(avtol, i) = RCONST(1.0e-8);
        if (j == N_center) j++;
        Ith(yy, i) = T_vect[j];
        //cout << "yy " << i << " = " << Ith(yy, i) << "\n";
        j++;
    }

    Ith(avtol, NEQ_T) = RCONST(1.0e-12);
    Ith(yy, NEQ_T) = M;
    data->N_m = NEQ_T - 1;
    //cout << "yy " << NEQ_T << " = " << Ith(yy, NEQ_T) << "\n";

    for (int i = NEQ_T + 1; i <= NEQ; i++) {
        Ith(avtol, i) = RCONST(1.0e-12);
        Ith(yy, i) = Y_vect[i - NEQ_T];
        //cout << "yy " << i << " = " << Ith(yy, i) << "\n";
    }

    j = 1;
    for (int i = 1; i < NEQ_T; i++) {
        if (j == N_center) j++;
        Ith(yp, i) = -F_right(T_vect[j - 1], T_vect[j], T_vect[j + 1], data->M, Y_vect[j], Y_vect[j + 1], data->x, j);
        //cout << "Ithypi = " << -F_right(T_vect[j - 1], T_vect[j], T_vect[j + 1], data->M, Y_vect[j], Y_vect[j + 1], data->x, j) << "\n";
        j++;

    }
    Ith(yp, NEQ_T) = 0;
    for (int i = NEQ_T + 1; i <= NEQ; i++) {
        j = i - NEQ_T;
        Ith(yp, i) = -F_rightY(T_vect[j - 1], T_vect[j], T_vect[j + 1], data->M, Y_vect[j - 1], Y_vect[j], Y_vect[j + 1], data->x, j);
    }
    */
    //
    /* Call IDACreate and IDAInit to initialize IDA memory */
    mem = IDACreate(ctx);
    if (check_retval((void*)mem, "IDACreate", 0)) return(1);

    retval = IDAInit(mem, resrob, t0, yy, yp);
    if (check_retval(&retval, "IDAInit", 1)) return(1);
    /* Call IDASVtolerances to set tolerances */

    retval = IDASVtolerances(mem, rtol, avtol);
    if (check_retval(&retval, "IDASVtolerances", 1)) return(1);
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
    
    //Write Initial distribution of Temperatures
    ofstream fout;
    int i;
    fout.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/1file(Temp)/" +
        to_string(call) + "file" + to_string(0) + ".dat");
    fout << "TITLE=\"" << "Graphics" << "\"" << endl;
    fout << R"(VARIABLES= "rj, m", "T, K", "Lambda, Watt/m*K", "Cp, J/kg*K", "rho, kg/m^3")" << endl;
    i = 0;
    for (i; i < s + 1; i++) {
        fout << x_vect[i] << "  " << T_vect[i] << "  "
            << ALambda[i] << "  " << ACp[i] << "  " << ADensity[i] << "\n";
    }
    //Writing on interface
    //Right side
    fout << inter << "  " << T_vect[N_x] << " "
        << Lyambda(T_vect[N_x], a, 0) << "  " << Cp(T_vect[N_x], a) << "  " << Density(T_vect[N_x], a) << "\n";
    //Left side
    fout << inter << "  " << T_vect[N_x] << " "
        << Lyambda(T_vect[N_x], a, 1) << "  " << Cp(T_vect[N_x], a) << "  " << Density(T_vect[N_x], a) << "\n";
    for (i; i < N_x; i++) {
        fout << x_vect[i] << "  " << T_vect[i] << "  "
            << ALambda[i] << "  " << ACp[i] << "  " << ADensity[i] << "\n";
    }
    fout.close();
    /* In loop, call IDASolve, print results, and test for error.
       Break out of loop when NOUT preset output times have been reached. */

    iout = 0; tout = tout1;
    data->n_tout = (tout / tout1);
    double tend = pow(10, 5);
    double Y_H2, Y_O2;
    double W, w_dot, rho;
    while (iout < print_value) {
        retval = IDASolve(mem, tout, &tret, yy, yp, IDA_NORMAL);
        ExportToArray(T_vect, Y_vect, M, data, yy, N_x);
        //PrintOutput(mem, tret, yy);
        if ((iout + 1) % number == 0 || (iout + 1) < 16)
        {
            cout << "t = " << tout << "\n";
            cout << "M = " << M << "\n";
            fout.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/1file(Temp)/" +
                to_string(call) + "file" + to_string(int(tout / tout1)) + ".dat");
            fout << "TITLE=\"" << "Graphics" << "\"" << endl;
            fout << R"(VARIABLES= "rj, m", "T, K", "Lambda, Watt/m*K", "Cp, J/kg*K", "rho, kg/m^3")" << endl;
            i = 0;
            for (i; i < s + 1; i++) {
                fout << x_vect[i] << "  " << T_vect[i] << "  "
                    << ALambda[i] << "  " << ACp[i] << "  " << ADensity[i] << "\n";
            }
            //Writing on interface
            //Right side
            fout << inter << "  " << T_vect[N_x] << " " 
                << Lyambda(T_vect[N_x], a, 0) << "  " << Cp(T_vect[N_x], a) << "  " << Density(T_vect[N_x], a) << "\n";
            //Left side
            fout << inter << "  " << T_vect[N_x] << " "
                << Lyambda(T_vect[N_x], a, 1) << "  " << Cp(T_vect[N_x], a) << "  " << Density(T_vect[N_x], a) << "\n";
            for (i; i < N_x; i++) {
                fout << x_vect[i] << "  " << T_vect[i] << "  "
                    << ALambda[i] << "  " << ACp[i] << "  " << ADensity[i] << "\n";
            }
            fout.close();
        }
        if (check_retval(&retval, "IDASolve", 1)) return(1);
        iout++;
        tout += tout1;
        data->n_tout = (tout / tout1);
    }

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

void Add_elem_spec(vector<double>& T, vector<double>& Y, vector<double>& x, int& N_x, int& N_center, double b)
{
    int j_t = 1;
    double T_max = 0, T_min = T[0];
    double Y_max = 0, Y_min = Y[N_x - 1];
    for (int i = 0; i < N_x; i++)
    {
        if (T[i] > T_max) T_max = T[i];
        if (T[i] < T_min) T_min = T[i];
    }

    while (j_t < N_x - 2)
    {
        if (fabs(T[j_t] - T[j_t - 1]) > b * (T_max - T_min))
        {
            T.insert(T.begin() + j_t, (T[j_t] + T[j_t - 1]) / 2.);
            Y.insert(Y.begin() + j_t, (Y[j_t] + Y[j_t - 1]) / 2.);
            x.insert(x.begin() + j_t, (x[j_t] + x[j_t - 1]) / 2.);
            N_x++;
            j_t++;
        }
        j_t++;
        //cout << "j_t = " << j_t << "\n";
    }

    j_t = 1;
    for (int i = 0; i < N_x; i++)
    {
        if (Y[i] > Y_max) Y_max = Y[i];
        if (Y[i] < Y_min) Y_min = Y[i];
    }

    while (j_t < N_x - 2)
    {
        if (fabs(Y[j_t] - Y[j_t - 1]) > b * (Y_max - Y_min))
        {
            T.insert(T.begin() + j_t, (T[j_t] + T[j_t - 1]) / 2.);
            Y.insert(Y.begin() + j_t, (Y[j_t] + Y[j_t - 1]) / 2.);
            x.insert(x.begin() + j_t, (x[j_t] + x[j_t - 1]) / 2.);
            N_x++;
            j_t++;
        }
        j_t++;
        //cout << "j_t = " << j_t << "\n";
    }
    for (int i = 0; i < N_x - 1; i++) {
        if (x[i] <= koeff_l * l && x[i + 1] > koeff_l * l)
            N_center = i;
    }
}

void Write_to_file(string str, ofstream& fout, vector<double>& x_vect,
    vector<double>& T_vect, vector<double>& Y_vect, double M, int N_x, int number) {
    double x_start, x_finish;
    fout.open(str + ".dat");
    fout << "TITLE=\"" << "Graphics" << "\"" << endl;
    fout << R"(VARIABLES= "x", "T", "Y fr", "rho", "v", "Cp", "lambda", "D")" << endl;
    double D, rho;
    double W;
    double Cp;
    for (int i = 0; i < N_x; i++) {
        W = get_W(0, 0, 0, 0);
        rho = P * W / R / T_vect[i];
        D = Lambda(T_vect[i], Y_vect[i]) / rho / Cp_all(0, 0);
        fout << x_vect[i] << "  " << T_vect[i] << " " << Y_vect[i] / (1. - Y_N2)
            << " " << P * mol_weight / R / T_vect[i]
            << " " << M / (P * mol_weight / R / T_vect[i])
            << " " << Cp_all(0, 0)
            << " " << Lambda(T_vect[i], Y_vect[i])
            << " " << D << endl;
    }
    fout.close();
    cout << "N_X " + to_string(number) + " = " << N_x << "\n";
    cout << "M_" + to_string(number) + "_block = " << M << endl;
    rho = P * W / R / T_vect[0];
    cout << "v = " << M / rho << "\n";
    for (int i = 1; i < N_x - 1; i++)
    {
        if (Y_vect[i] <= 0.1 * (1 - Y_N2) && Y_vect[i + 1] > 0.1 * (1 - Y_N2)) x_start = x_vect[i];
        if (Y_vect[i] <= 0.9 * (1 - Y_N2) && Y_vect[i + 1] > 0.9 * (1 - Y_N2)) x_finish = x_vect[i];
    }
    cout << "xfstat = " << x_finish - x_start << "\n";
    fout.open("h_" + to_string(number) + ".dat");
    fout << "TITLE=\"" << "Graphics" << "\"" << endl;
    fout << R"(VARIABLES= "x", "h_")" << endl;
    for (int i = 0; i < N_x - 1; i++) {
        fout << x_vect[i] << "  " << x_vect[i + 1] - x_vect[i] << endl;
    }
    fout.close();
}

void MySolver(int& retval, vector<double>& Y_vect, double& M, int& N_center)
{
    //InTime
    cout << "INTEGRATE ALL\n\n\n";
    double tout1 = pow(10, -4);
    int number = 100;
    int print_value = pow(10, 5);
    cout << "print_value" << print_value << "\n";
    //Mymain
    const int const_params = 1;
    const double T_l = 300;
    const double T_r = 1500;
    //Define initialtempereature for surface, K
    double T_cur_i = 300;
    //зададим минимальный возможный размер ячейки
    const double h_min = 0.000025;                                 //0.0005 * 2 / 41 = 0.0000243902439024= 24.3902439024 мкр;  
    //Number of cells
    int N = 101;                          //количество узлов. То есть ячеек N - 1
    const int Nd = 20;
    const int N_uni = 2 * Nd;
    const int N_uni_near_center = 5;
    const double x_l = 0.;
    //const double R = 0.0005 м = 500 мкр
    const double x_r = 0.0025; //для простой задачи зададим такую область
    double h = 0.2;
    double x_uni = x_l + N_uni * h_min;
    double q = DefineQ(N - N_uni, h_min, x_r, x_uni);
    vector <double> r;
    vector <double> T_cur;
    T_cur.resize(N + 1);
    r.resize(N);
    InitialGrid(N, x_l, T_l, T_r, r, T_cur, h, q, h_min, const_params, Nd, N_uni, N_uni_near_center);
    int s = Nd;
    T_cur[N] = T_cur_i;
    retval = Integrate_IDA(N, r, T_cur, Y_vect, M, N_center, tout1, 1, number, print_value, 0, s);
}

int main()
{
    init_consts(num_gas_species, num_react);
    int N_x = 250;
    double b = 0.01;
    double M;
    double W, rho, Y_H2, Y_O2;
    int N_center;
    int retval;
    double w_dot;
    vector<double> x_vect(N_x);
    vector<double> T_vect(N_x);
    vector<double> Y_vect(N_x);
    double* my_x;
    ofstream fout;
    //double tout1 = pow(10, -7);
    //int number = 500;
    //int print_value = 8000;
    double h, h_left, dTdx, maxdTdx = 0;
    double x_start, x_finish;
    int cons_flag = 0;
    {
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

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY KINSOL
 *--------------------------------------------------------------------
 */

 /*
  * System function for predator-prey system
  */


static int func(N_Vector u, N_Vector f, void* user_data)
{
    realtype* T, * fdata;
    realtype x1, l1, L1, x2, l2, L2;
    realtype* x_cells, * T_vect, * Y_vect;
    UserData data;
    double h_left, h;
    double tmp;
    double M;
    data = (UserData)user_data;
    x_cells = data->x;
    T_vect = data->T;
    Y_vect = data->Y_H2O;
    int myNx = data->Nx;
    int myNeq = data->NEQ;
    int myNm = data->N_m;

    T = N_VGetArrayPointer(u);
    fdata = N_VGetArrayPointer(f);
    double T_cer = data->T_center;
    int j = 0;
    //настоящий вектор T
    ExportToArray(T_vect, Y_vect, M, data, u, myNx);
    //cout << "Y_vect  " << j << " =  " << Y_vect[j] << endl;
    int fl = 0;
    //cout << "cycle = " << cycle << "\n";
    long int nniters;
    int retval = KINGetNumNonlinSolvIters(data->mykmem, &nniters);
    ofstream foutT;
    ofstream foutY;
    /*if (myiter < nniters)
    {
        foutT.open(to_string(nniters) + "T.dat");
        foutT << "TITLE=\"" << "Graphics" << "\"" << endl;
        foutT << R"(VARIABLES= "N", "FT")" << endl;
        cout << "nniters = " << nniters << "\n";
        foutY.open(to_string(nniters) + "Y.dat");
        foutY << "TITLE=\"" << "Graphics" << "\"" << endl;
        foutY << R"(VARIABLES= "N", "FY")" << endl;
        cout << "nniters = " << nniters << "\n";
    }*/

    for (j = 1; j < myNx - 1; j++) {
        fdata[j - 1] = F_right(T_vect[j - 1], T_vect[j], T_vect[j + 1], data->M, Y_vect[j], Y_vect[j + 1], x_cells, j);
        //cout << j - 1 << " " << fdata[j - 1] << "\n";
        //if (myiter < nniters)
        //{
        //    foutT << j - 1 << " " << fdata[j - 1] << "\n";
        //    //cout << j - 1 << " " << fdata[j - 1] << "\n";
        //}
    }
    //cout << "next\n";
    for (int i = 1; i < myNx - 1; i++) {
        fdata[j - 1] = F_rightY(T_vect[i - 1], T_vect[i], T_vect[i + 1], data->M, Y_vect[i - 1], Y_vect[i], Y_vect[i + 1], x_cells, i);
        //cout << "j - 1 = " << j - 1 << "\n";
        /*if (myiter < nniters)
        {
            foutY << i - 1 << " " << fdata[j - 1] << "\n";
        }*/
        j++;
    }
    cycle++;
    /* if (myiter < nniters)
     {
         myiter++;
         foutT.close();
         foutY.close();
     }*/
     //cout << "\n" << "\n" << "\n" << "\n" << "\n" << "\n";
    return(0);
}

static int func_Y(N_Vector u, N_Vector f, void* user_data)
{
    realtype* Y, * fdata;
    realtype x1, l1, L1, x2, l2, L2;
    realtype* x_cells, * T_vect, * Y_vect;
    UserData data;
    double h_left, h;
    double tmp;
    double M;
    data = (UserData)user_data;
    x_cells = data->x;
    T_vect = data->T;
    Y_vect = data->Y_H2O;
    int myNx = data->Nx;
    int myNeq = data->NEQ;
    int myNm = data->N_m;

    Y = N_VGetArrayPointer(u);
    fdata = N_VGetArrayPointer(f);
    double T_cer = data->T_center;

    //cout << "MyNx = " << myNx << "\n";
    //cout << "MyNm = " << myNm << "\n";
    data->M = Y[myNm];
    int j = 1;
    //cout << "M = " << data->M << "\n";
    for (int i = myNm + 1; i < myNx - 1; i++)
    {
        Y_vect[i] = Y[i];
        //cout << "Y_vect5555  " << i << " =  " << Y_vect[i] << endl;
    }
    Y_vect[0] = 0.;
    Y_vect[myNx - 1] = Y_vect[myNx - 2];
    int fl = 0;
    int Ncentr = data->N_centr;

    fdata[0] = F_right(T_vect[Ncentr - 1], T_vect[Ncentr], T_vect[Ncentr + 1], data->M, Y_vect[Ncentr], Y_vect[Ncentr + 1], x_cells, Ncentr);
    //cout << "fdata " << 0 << " = " << fdata[0] << "\n";
    //cout << "next\n";

    for (int i = 1; i < myNeq; i++) {
        fdata[i] = F_rightY(T_vect[i - 1], T_vect[i], T_vect[i + 1], data->M, Y_vect[i - 1], Y_vect[i], Y_vect[i + 1], x_cells, i);
        //cout << "j - 1 = " << j - 1 << "\n";
        //cout << "fdata " << i << " = " << fdata[i] << "\n";
    }
    //cycle++;
    //cout << "cycle = " << cycle << "\n";
    //cout << "\n" << "\n" << "\n" << "\n" << "\n" << "\n";
    return(0);
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
    int n_dt;
    n_dt = data->n_tout;
    int k;
    k = data->kp;
    int j;
    yval = N_VGetArrayPointer(yy);
    ypval = N_VGetArrayPointer(yp);
    rval = N_VGetArrayPointer(rr);
    //cout << "ypvalres0 = " << yval[0] << "\n";
    int myNx = data->Nx;
    ExportToArray(T_vect, Y_vect, data->M, data, yy, data->Nx);
    double Y_H2, Y_O2;
    double W;
    double rho;

    //MAXIM added
    double p = 1.5;
    int s = data->N_centr;
    cout << "s= " << s << "\n";
    cout << "T[s]= " << T_vect[s] << "\n";
    double inter = x_cells[s] + (p - 1) * (x_cells[s + 1] - x_cells[s]);
    //double q = 20 * pow(10, 3);                  //W / m^2
    double L_d = 2258.2 * pow(10, 3.);              //J / kg 
    //const determination of dM
    //dM = 4 * PI * pow(inter, 2.0) * q / L_d;     // kg / s
    //dM equals zero because it is only warming time
    double dM = 0.;                                       // kg / s
    double a = 1.;
        //pow(10., -10);
    vector <double> f(myNx);
    vector <double> ACp, ADfCp, ADensity, ADfDensity, ALambda, ADfLambda;
    ArraysParametersCopyForArrays(a, myNx, x_cells, T_vect, s, ACp, ADensity, ALambda, n_dt, inter, T_vect[myNx], p);
    double q_i_g = 0;
    double q_i_d = 0;
    FRightTestCopyForArrays(f, T_vect, x_cells, ACp, ADensity, ALambda, n_dt, myNx, a, dM, s, inter, T_vect[myNx], p, L_d,
        q_i_g, q_i_d);
    for (j = 0; j < myNx - 1; j++) {
            rval[j] = f[j] - ypval[j];
        //cout << "right  = " << j - 1 << "  k = " << k - 1 << "\n";
    }
    rval[myNx - 1] = f[myNx - 1];

    //Calculate modul of Nevyazka
    double mod_nevyaz = 0;
    for (int i = 0; i < myNx; i++)
        mod_nevyaz += pow(rval[i], 2.0);
    mod_nevyaz = pow(mod_nevyaz, 0.5);
    //
    *(data->outNevyaz) << k << " " << mod_nevyaz << " " << abs(q_i_g / pow(inter, 2.)) << " "
        << abs(q_i_d / pow(inter, 2.)) << " " << "\n";

    /////////////////////////////////
    //Opening file for T on iteration
    //Define nevyazkas on each iteration
    ofstream OutIterNevyazka;
    OutIterNevyazka.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/Data_new/IterNevyazka/IterNevyazka_" + to_string(k) + ".dat");
    OutIterNevyazka << "TITLE=\"" << "Graphics" << "\"" << "\n";
    OutIterNevyazka << R"(VARIABLES= "rj, m", "F" )" << "\n";
    int i = 0;
    for (i; i < s + 1; i++)
        OutIterNevyazka << x_cells[i] << " " << abs(rval[i]) << "\n";
    //Interface
    OutIterNevyazka << inter << " " << abs(rval[myNx - 1]) << "\n";
    for (i; i < myNx - 1; i++)
        OutIterNevyazka << x_cells[i] << " " << abs(rval[i]) << "\n";
    OutIterNevyazka.close();
    //Define Temperature on each iteration
    ofstream OutIterCurrentTemp;
    OutIterCurrentTemp.open("C:/Users/user/source/Научная работа/DropletEvaporation(IDA)/Data_new/IterTemp/IterTemp_" + to_string(k) + ".dat");
    OutIterCurrentTemp << "TITLE=\"" << "Graphics" << "\"" << "\n";
    OutIterCurrentTemp << R"(VARIABLES= "rj, m", "T, K", "Lambda, W/m*K" )" << "\n";
    i = 0;
    for (i; i < s + 1; i++)
        OutIterCurrentTemp << x_cells[i] << " " << T_vect[i] << " " << ALambda[i] << "\n";
    //Interface
    OutIterCurrentTemp << inter << " " << T_vect[myNx] << " " << ALambda[s] << "\n";
    for (i; i < myNx; i++)
        OutIterCurrentTemp << x_cells[i] << " " << T_vect[i] << " " << ALambda[i] << "\n";
    OutIterCurrentTemp.close();

    k++;
    data->kp = k;
    //////////////////////////////////////////////
    /*
    for (j = 1; j < myNx - 1; j++) {
        if (j != data->N_centr)
        {
            get_Y(Y_vect[j], Y_H2, Y_O2, Y_N2);
            W = get_W(Y_vect[j], Y_H2, Y_O2, Y_N2);
            //cout << "W = " << W << "\n";
            rho = P * W / R / T_vect[j];
            rval[j - 1] = F_right(T_vect[j - 1], T_vect[j], T_vect[j + 1], data->M, Y_vect[j], Y_vect[j + 1], x_cells, j) / rho + ypval[k];
            k++;
        }
        else
        {
            rval[j - 1] = F_right(T_vect[j - 1], T_vect[j], T_vect[j + 1], data->M, Y_vect[j], Y_vect[j + 1], x_cells, j);
        }
        //cout << "right  = " << j - 1 << "  k = " << k - 1 << "\n";
    }
    //cout << "next\n";
    //cout << "N_m = " << data->N_m << "\n";
    //cout << "M = " << yval[data->N_m] << "\n";
    for (int i = 1; i < myNx - 1; i++) {
        get_Y(Y_vect[i], Y_H2, Y_O2, Y_N2);
        W = get_W(Y_vect[i], Y_H2, Y_O2, Y_N2);
        //cout << "W = " << W << "\n";
        rho = P * W / R / T_vect[i];
        rval[j - 1] = F_rightY(T_vect[i - 1], T_vect[i], T_vect[i + 1], data->M, Y_vect[i - 1], Y_vect[i], Y_vect[i + 1], x_cells, i) / rho
            + ypval[j - 1];
        //cout << "right  = " << j - 1 << "  k = " << j - 1 << "\n";
        j++;
    }
    */
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
