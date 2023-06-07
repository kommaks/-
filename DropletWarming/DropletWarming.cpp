// DropletWarming.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;
#define PI 3.14159

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
	return Lagrange1(LAGR_T, LAGR_COND_copy, nx, T);
}
//Formula setting for Heat Capacity and it's derivative, J/(kg*K)
double Cp(double T, const double a)
{
	const double Cp_water = 4180.;
	const double Cp_steam = 2000.;
	//Molar mass of water
	const double W = 18 * pow(10, -3);
	const double T_boiling = 373;
	double AA, B, C, D, E;
	//Polinom for T = (K) /1000
	double T_forCp = T / 1000.;
	if (T <= T_boiling)
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
		if (T <= T_boiling)
			return Cp_water;
		else
			return Cp_steam;
	}
	else
		return (AA + B * T_forCp + C * pow(T_forCp, 2.) + D * pow(T_forCp, 3.) + E * pow(T_forCp, -2.)) / W;
}
//J/(kg*K^2)
double DfCp(double T, const double a)
{
	//Molar mass of water
	const double W = 18 * pow(10, -3);
	const double T_boiling = 373;
	double AA, B, C, D, E;
	//Polinom for T = (K) /1000
	double T_forCp = T / 1000.;
	if (T <= T_boiling)
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

}
//Formula setting for thermal conductivity and it's derivative, Watt/(m*K)
double Lambda(double T, const double a, int flag_phase)
{
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
}
//Watt/(m*K^2)
double DfLambda(double T, const double a, int flag_phase)
{
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
}
//Setting arrays of parametres depending on T:
void ArraysParameters(const double a, const double b, const double c, const int N, vector <double>& r, vector <double>& T_next,
	const int s, vector <double>& ACp, vector <double>& ADfCp, vector <double>& ADensity, vector <double>& ADfDensity,
	vector <double>& ALambda, vector <double>& ADfLambda, int n, double inter, double Ti, double p)
{
	//Define Parameters
	ACp.resize(N + 1);
	ADfCp.resize(N + 1);
	ADensity.resize(N + 1);
	ADfDensity.resize(N + 1);
	ALambda.resize(N + 1);
	ADfLambda.resize(N + 1);
	/*
	vector <double> Tg;
	int dTg = 1;
	int Ng = (T_r - T_l) / dTg;
	Tg.resize(Ng + 1);
	Tg[0] = T_l;
	//Setting Temperature array;
	for (int i = 0; i < N; i++) {
		Tg[i] = Tg[i - 1] + dTg;
	}
	*/
	for (int i = 0; i < s + 1; i++) {
		//cout << r[i] << endl;
		ACp[i] = Cp(T_next[i], a);
		ADfCp[i] = DfCp(T_next[i], a);
		ADensity[i] = Density(T_next[i], a);
		ADfDensity[i] = DfDensity(T_next[i], a);
		ALambda[i] = Lambda(T_next[i], a, 0);
		ADfLambda[i] = DfLambda(T_next[i], a, 0);
	}
	for (int i = s + 1; i < N + 1; i++) {
		//cout << r[i] << endl;
		ACp[i] = Cp(T_next[i], a);
		ADfCp[i] = DfCp(T_next[i], a);
		ADensity[i] = Density(T_next[i], a);
		ADfDensity[i] = DfDensity(T_next[i], a);
		ALambda[i] = Lambda(T_next[i], a, 1);
		ADfLambda[i] = DfLambda(T_next[i], a, 1);
	}

	//Plot graphics of Parameters as a function of Temperature
	ofstream Parameters;
	Parameters.open("C:/Users/user/source/Научная работа/DropletWarming/Data_new/Parameters/Parameters_" + to_string(n) + ".dat");
	Parameters << "TITLE=\"" << "Graphics" << "\"" << endl;
	Parameters << R"(VARIABLES= "rj, m", "T, K", "Cp, J/kg*K", "DfCp, J/kg*K^2", "rho, kg/m^3", "Dfrho, kg/m^3*K", "Lambda, Watt/m*K", "DfLambda, Watt/m*K^2")" << "\n";
	int j = 0;
	for (j; j < s + 1; j++) {
		Parameters << r[j] << " " << T_next[j] << " " << ACp[j] << " " << ADfCp[j] << " " << ADensity[j] << " " << ADfDensity[j] << " "
			<< ALambda[j] << " " << ADfLambda[j] << "\n";
	}
	//Interface
	Parameters << inter << " " << Ti << " " << ACp[s] << " " << ADfCp[s] << " " << ADensity[s] << " " << ADfDensity[s] << " "
		<< ALambda[s] << " " << ADfLambda[s] << "\n";
	for (j + 1; j < N + 1; j++) {
		Parameters << r[j] << " " << T_next[j] << " " << ACp[j] << " " << ADfCp[j] << " " << ADensity[j] << " " << ADfDensity[j] << " "
			<< ALambda[j] << " " << ADfLambda[j] << "\n";
	}
	Parameters.close();
}
void GraphicsSystEqu(int n, vector<double>& r, vector<double>& T_next, vector<double>& ALambda, vector<double>& ADensity,
	const int N, const int s, double dM, int flag_evap, double inter, double Ti)
{
	ofstream OutCurrentTemp;
	ofstream OutFlow;
	if (n % 10 == 0 || n < 40 || flag_evap == 1)
	{
		OutCurrentTemp.open("C:/Users/user/source/Научная работа/DropletWarming/Data_new/Temp/Temp_" + to_string(n) + ".dat");
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
	OutFlow.open("C:/Users/user/source/Научная работа/DropletWarming/Data_new/Flow/Flow_" + to_string(n) + ".dat");
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
	OutSurfaceFlow.open("C:/Users/user/source/Научная работа/DropletWarming/Data_new/SurfaceFlow.dat");
	OutSurfaceFlow << "TITLE=\"" << "Graphics" << "\"" << "\n";
	OutSurfaceFlow << R"(VARIABLES= "t, s", "T0, K", "Ts, K", "Ti, K", "Ts+1, K", "gradT_d, K/m", "gradT_g, K/m", "q_d, W", "q_g, W", "dM, kg/s",
 "rho, kg/m^3", "u_r, m/s", "p, m/(cell size) ", "inter, m" )" << "\n";
	//Collecting Last Step Nevyazka
	OutLastStepNevyazka.open("C:/Users/user/source/Научная работа/DropletWarming/Data_new/LastStepNevyazka.dat");
	OutLastStepNevyazka << "TITLE=\"" << "Graphics" << "\"" << endl;
	OutLastStepNevyazka << R"(VARIABLES= "t, s", "F")" << endl;
	//Collecting First Step Nevyazka
	OutFirstStepNevyazka.open("C:/Users/user/source/Научная работа/DropletWarming/Data_new/FirstStepNevyazka.dat");
	OutFirstStepNevyazka << "TITLE=\"" << "Graphics" << "\"" << endl;
	OutFirstStepNevyazka << R"(VARIABLES= "t, s", "F")" << endl;
}
//bad results wuth this function (was declined because of it)
double LambdaBetween(vector<double>& r, vector<double>& ALambda, int j, double r_between)
{
	return (r[j + 1] - r[j]) / ((r_between - r[j]) / ALambda[j] + (r[j + 1] - r_between) / ALambda[j + 1]);
}
//Setting j-1 element of matrix
void DfInterface(vector<double>& r, vector<double>& T_next, vector<vector<double>>& J, vector <double>& ACp, vector <double>& ADfCp, vector <double>& ADensity,
	vector <double>& ADfDensity, vector <double>& ALambda, vector <double>& ADfLambda, const int N, const double a, double dt, double dtau, double dM,
	double inter, double Ti, double p, int s)
{
	double h = (r[s + 1] + r[s]);
	double term1, term2, term3, term4;

	//To the left of the interface
	double rs_avg = (r[s] + r[s - 1]) / 2.0;
	double rs_diff = inter - rs_avg;
	J[s][s - 2] = (1.0 / rs_diff) * pow(inter, 2) * Lambda(Ti, a, 0) * (p / ((p + 1.0) * h));
	J[s][s - 1] = (1.0 / rs_diff) * (-pow(inter, 2) * Lambda(Ti, a, 0) * (p + 1.0) / (p * h) + pow(rs_avg, 2) * 0.5 * (ALambda[s] + ALambda[s - 1])
		/ (r[s] - r[s - 1]));
	J[s][s] = -(1.0 / rs_diff) * pow(rs_avg, 2) * 0.5 * (ALambda[s] + ALambda[s - 1]) / (r[s] - r[s - 1]) - ACp[s] * ADensity[s] * pow(r[s], 2) / dt;
	J[s][N] = (1.0 / rs_diff) * pow(inter, 2) * ((DfLambda(Ti, a, 0) * (p / ((p + 1.0) * h) * T_next[s - 2] - (p + 1.0) / (p * h) * T_next[s - 1] + (2.0 * p + 1.0)
		/ ((p + 1.0) * p * h) * Ti)) + Lambda(Ti, a, 0) * (2.0 * p + 1.0) / ((p + 1.0) * p * h));

	//On the interface(i = N)
	J[N][s - 2] = -pow(inter, 2.) / h * Lambda(Ti, a, 0) * (p / (p + 1.0));
	J[N][s - 1] = pow(inter, 2.) / h * Lambda(Ti, a, 0) * (p + 1.0) / p;
	J[N][N] = pow(inter, 2.) / h * (Lambda(Ti, a, 1) * (2.0 * p - 7.0) / ((p - 3.0) * (p - 4.0))
		- Lambda(Ti, a, 0) * (2.0 * p + 1.0) / ((p + 1.0) * p));
	J[N][s + 2] = pow(inter, 2.) / h * Lambda(Ti, a, 1) * (p - 4.0) / (p - 3.0);
	J[N][s + 3] = pow(inter, 2.) / h * Lambda(Ti, a, 1) * (p - 3.0) / (4.0 - p);

	//To the right of the interface
	double rs_plus_avg = (r[s + 1] + r[s + 2]) / 2.0;
	double rs_plus_diff = rs_plus_avg - inter;
	J[s + 1][N] = -(1.0 / rs_plus_diff) * pow(inter, 2.) * ((DfLambda(Ti, a, 1) * ((2.0 * p - 7.0) / ((p - 3.0) * (p - 4.0) * h) * Ti
		+ (p - 4.0) / ((p - 3.0) * h) * T_next[s + 2] + (p - 3.0) / ((4.0 - p) * h) * T_next[s + 3])) + Lambda(Ti, a, 1) * ((2.0 * p - 7.0) / ((p - 3.0) * (p - 4.0) * h)))
		+ ACp[s + 1] * dM * (1.0 / (r[s + 1] - inter));
	J[s + 1][s + 1] = -(1.0 / rs_plus_diff) * pow(rs_plus_avg, 2) * (0.5 * (ALambda[s + 2] + ALambda[s + 1]) / (r[s + 2] - r[s + 1]))
		- ACp[s + 1] * (ADensity[s + 1] * pow(r[s + 1], 2.) / dt + dM / (r[s + 1] - inter));
	J[s + 1][s + 2] = (1.0 / rs_plus_diff) * (pow(rs_plus_avg, 2) * (0.5 * (ALambda[s + 2] + ALambda[s + 1]) / (r[s + 2] - r[s + 1]))
		- pow(inter, 2.) * Lambda(Ti, a, 1) * (p - 4.0) / ((p - 3.0) * h));
	J[s + 1][s + 3] = -(1.0 / rs_plus_diff) * pow(inter, 2.) * Lambda(Ti, a, 1) * (p - 3.0) / ((4.0 - p) * h);
}
double DfLeft(vector<double>& r, vector<double>& T_next, vector<double>& ALambda, vector<double>& ADfLambda,
	const double a, const double b, const double c, int j)
{
	return 2. / (r[j + 1] - r[j - 1]) * pow((r[j] + r[j - 1]) / 2., 2.) * (-0.5 * (ADfLambda[j] + ADfLambda[j - 1])
		* (T_next[j] - T_next[j - 1]) / (r[j] - r[j - 1]) + 0.5 * (ALambda[j] + ALambda[j - 1]) / (r[j] - r[j - 1]));
}
//Setting j element of matrix
double DfCenter(vector<double>& r, vector<double>& T_next, vector <double>& ACp, vector <double>& ADfCp, vector <double>& ADensity,
	vector <double>& ADfDensity, vector <double>& ALambda, vector <double>& ADfLambda, const double a, const double b, const double c,
	const double d, int j, const double dt, const double dtau)
{
	return 2. / (r[j + 1] - r[j - 1]) * (pow((r[j + 1] + r[j]) / 2., 2.) * (0.5 * (ADfLambda[j + 1] + ADfLambda[j])
		* (T_next[j + 1] - T_next[j]) / (r[j + 1] - r[j]) - 0.5 * (ALambda[j + 1] + ALambda[j]) / (r[j + 1] - r[j]))
		- pow((r[j] + r[j - 1]) / 2., 2.) * (0.5 * (ADfLambda[j] + ADfLambda[j - 1]) * (T_next[j] - T_next[j - 1]) / (r[j] - r[j - 1])
			+ 0.5 * (ALambda[j] + ALambda[j - 1]) / (r[j] - r[j - 1])))
		- (ACp[j] * ADensity[j] + (ADfCp[j] * ADensity[j] + ACp[j] * ADfDensity[j]) * T_next[j]) * pow(r[j], 2.) / dt
		- ACp[j] * ADensity[j] * pow(r[j], 2.) / dtau;
}
//Setting j+1 element of matrix
double DfRight(vector<double>& r, vector<double>& T_next, vector<double>& ALambda, vector<double>& ADfLambda,
	const double a, const double b, const double c, int j)
{
	return 2. / (r[j + 1] - r[j - 1]) * pow((r[j + 1] + r[j]) / 2., 2.) * (0.5 * (ADfLambda[j + 1] + ADfLambda[j])
		* (T_next[j + 1] - T_next[j]) / (r[j + 1] - r[j]) + 0.5 * (ALambda[j + 1] + ALambda[j]) / (r[j + 1] - r[j]));
}
//Df for Adiabatic left boundary
double AdiabaticBoundaryDfCenter(vector<double>& r, vector<double>& T_next, vector <double>& ACp, vector <double>& ADfCp, vector <double>& ADensity,
	vector <double>& ADfDensity, vector <double>& ALambda, vector <double>& ADfLambda,
	const double a, const double b, const double c, const double dt, const double dtau)
{
	int l_f_num = 0;
	return 2. / (r[2] - r[0]) * (pow((r[l_f_num + 1] + r[l_f_num]) / 2., 2.) * (0.5 * (ADfLambda[l_f_num + 1] + ADfLambda[l_f_num])
		* (T_next[l_f_num + 1] - T_next[l_f_num]) / (r[l_f_num + 1] - r[l_f_num])
		- 0.5 * (ALambda[l_f_num + 1] + ALambda[l_f_num]) / (r[l_f_num + 1] - r[l_f_num])))
		- (ACp[l_f_num] * ADensity[l_f_num] + (ADfCp[l_f_num] * ADensity[l_f_num] + ACp[l_f_num] * ADfDensity[l_f_num]) * T_next[l_f_num]) * pow(r[l_f_num], 2.) / dt
		- ACp[l_f_num] * ADensity[l_f_num] * pow(r[l_f_num], 2.) / dtau;
}
double AdiabaticBoundaryDfRight(vector<double>& r, vector<double>& T_next, vector<double>& ALambda, vector<double>& ADfLambda,
	const double a, const double b, const double c, int j)
{
	return 2. / (r[2] - r[0]) * pow((r[j + 1] + r[j]) / 2., 2.) * (0.5 * (ADfLambda[j + 1] + ADfLambda[j])
		* (T_next[j + 1] - T_next[j]) / (r[j + 1] - r[j]) + 0.5 * (ALambda[j + 1] + ALambda[j]) / (r[j + 1] - r[j]));
}

void Jacobian(vector<vector<double>>& J, vector<double>& r, vector<double>& T_next, vector <double>& ACp, vector <double>& ADfCp,
	vector <double>& AD, vector <double>& ADfD, vector <double>& AL, vector <double>& ADfL, const int N_minus,
	const double a, const double b, const double c, const double d, const double dt, const double dtau, double dM, int s, double inter, double Ti, double p,
	int flag_evap)
{
	int l_f_num = 0;
	double r_temp, T_temp;
	//Jacobians for Water
	J[l_f_num][l_f_num] = AdiabaticBoundaryDfCenter(r, T_next, ACp, ADfCp, AD, ADfD, AL, ADfL, a, b, c, dt, dtau);
	J[l_f_num][l_f_num + 1] = AdiabaticBoundaryDfRight(r, T_next, AL, ADfL, a, b, c, l_f_num);
	for (int j = 1; j < s; j++) {
		J[j][j - 1] = DfLeft(r, T_next, AL, ADfL, a, b, c, j);
		J[j][j] = DfCenter(r, T_next, ACp, ADfCp, AD, ADfD, AL, ADfL, a, b, c, d, j, dt, dtau);
		J[j][j + 1] = DfRight(r, T_next, AL, ADfL, a, b, c, j);
	}
	//Jacobians for Gaze
	for (int j = s + 2; j < N_minus; j++) {
		J[j][j - 1] = DfLeft(r, T_next, AL, ADfL, a, b, c, j) + ACp[j] * dM / (r[j] - r[j - 1]);
		J[j][j] = DfCenter(r, T_next, ACp, ADfCp, AD, ADfD, AL, ADfL, a, b, c, d, j, dt, dtau) - ACp[j] * dM / (r[j] - r[j - 1]);
		J[j][j + 1] = DfRight(r, T_next, AL, ADfL, a, b, c, j);
	}
	J[N_minus][N_minus - 1] = DfLeft(r, T_next, AL, ADfL, a, b, c, N_minus) + ACp[N_minus] * dM / (r[N_minus] - r[N_minus - 1]);
	J[N_minus][N_minus] = DfCenter(r, T_next, ACp, ADfCp, AD, ADfD, AL, ADfL, a, b, c, d, N_minus, dt, dtau)
		- ACp[N_minus] * dM / (r[N_minus] - r[N_minus - 1]);
	//Define J[s], J[i], J[s + 1]
	DfInterface(r, T_next, J, ACp, ADfCp, AD, ADfD, AL, ADfL, N_minus + 1, a, dt, dtau, dM, inter, Ti, p, s);
}
void F(vector<double>& f, vector<double>& T_cur, vector<double>& T_next, vector<double>& r, vector <double>& ACp, vector <double>& ADensity, vector <double>& ALambda,
	const int N, const double a, const double b, const double c, const double d, const double dt, double dM, int s, double inter,
	const double Ti, double p, const double L_d, double& mod_nevyaz, int flag_evap)
{
	double F;
	//Nevyazka on the left border
	int l_f_num = 0;
	F = 2. / (r[2] - r[0]) * pow((r[l_f_num + 1] + r[l_f_num]) / 2., 2.) * 0.5 * (ALambda[l_f_num + 1] + ALambda[l_f_num])
		* (T_next[l_f_num + 1] - T_next[l_f_num]) / (r[l_f_num + 1] - r[l_f_num])
		- ACp[l_f_num] * ADensity[l_f_num] * pow(r[l_f_num], 2.) * (T_next[l_f_num] - T_cur[l_f_num]) / dt;
	f[l_f_num] = -F;
	//cout << "f[1]" << f[1] << "\n";
	//Nevyazka in Water
	for (int j = 1; j < s; j++)
	{
		F = 2. / (r[j + 1] - r[j - 1]) * (pow((r[j + 1] + r[j]) / 2., 2.) * 0.5 * (ALambda[j + 1] + ALambda[j])
			* (T_next[j + 1] - T_next[j]) / (r[j + 1] - r[j])
			- pow((r[j] + r[j - 1]) / 2., 2.) * 0.5 * (ALambda[j] + ALambda[j - 1])
			* (T_next[j] - T_next[j - 1]) / (r[j] - r[j - 1]))
			- ACp[j] * ADensity[j] * pow(r[j], 2.) * (T_next[j] - T_cur[j]) / dt;
		f[j] = -F;
		//cout << j << " f1_" << j << " " << f[j] << "\n";
	}
	//Nevyazka in Gaze
	for (int j = s + 2; j < N; j++)
	{
		//cout << "(T_next[j] - T[j][n])" << (T_next[j] - T[j][n]) << endl;
		F = 2. / (r[j + 1] - r[j - 1]) * (pow((r[j + 1] + r[j]) / 2., 2.) * 0.5 * (ALambda[j + 1] + ALambda[j])
			* (T_next[j + 1] - T_next[j]) / (r[j + 1] - r[j])
			- pow((r[j] + r[j - 1]) / 2., 2.) * 0.5 * (ALambda[j] + ALambda[j - 1])
			* (T_next[j] - T_next[j - 1]) / (r[j] - r[j - 1]))
			- ACp[j] * (ADensity[j] * pow(r[j], 2.) * (T_next[j] - T_cur[j]) / dt + dM * (T_next[j] - T_next[j - 1]) / (r[j] - r[j - 1]));
		f[j] = -F;
		//cout << j << " f1_" << j << " " << f[j] << "\n";
	}
	cout << "flag_evap= " << flag_evap << "\n";
	cout << "inter= " << inter << "\n";
	//Near interface Nevyazka
	double term1, term2, term3;
	double h = r[s + 1] - r[s];
	cout << "Ti" << Ti << "\n";
	//To the left of the interface
	double rs_avg = (r[s] + r[s - 1]) / 2.0;
	double rs_diff = inter - rs_avg;
	term1 = (1.0 / rs_diff) * (Lambda(Ti, a, 0) * pow(inter, 2.) * ((p / ((p + 1) * h)) * T_next[s - 2] - ((p + 1) / (p * h)) * T_next[s - 1]
		+ ((2 * p + 1) / ((p + 1) * p * h)) * Ti) - pow(rs_avg, 2) * 0.5 * (ALambda[s] + ALambda[s - 1]) * (T_next[s] - T_next[s - 1]) / (r[s] - r[s - 1]));
	term3 = -ACp[s] * ADensity[s] * pow(r[s], 2.) * (T_next[s] - T_cur[s]) / dt;
	cout << "term1= " << term1 << "\n";
	cout << "term3= " << term3 << "\n";
	cout << "q-(s)= " << -pow(rs_avg, 2) * 0.5 * (ALambda[s] + ALambda[s - 1]) * (T_next[s] - T_next[s - 1]) / (r[s] - r[s - 1]) << "\n";
	cout << "q+(s)= " << ALambda[s - 1] * ((p / ((p + 1) * h)) * T_next[s - 2] - ((p + 1) / (p * h)) * T_next[s - 1]
		+ ((2 * p + 1) / ((p + 1) * p * h)) * Ti) << "\n";
	F = term1 + term3;
	f[s] = -F;

	//On the interface
	term1 = pow(inter, 2.) * Lambda(Ti, a, 1) * ((2.0 * p - 7.0) / ((p - 3.0) * (p - 4.0) * h) * Ti + (p - 4.0) / ((p - 3.0) * h) * T_next[s + 2]
		+ (p - 3.0) / ((4.0 - p) * h) * T_next[s + 3]);
	term2 = pow(inter, 2.) * Lambda(Ti, a, 0) * (p / ((p + 1.0) * h) * T_next[s - 2] - (p + 1.0) / (p * h) * T_next[s - 1]
		+ (2.0 * p + 1.0) / ((p + 1.0) * p * h) * Ti);
	F = term1 - term2 + L_d * dM;
	cout << "term1= " << term1 << "\n";
	cout << "term2= " << term2 << "\n";
	cout << "L_d * dM " << L_d * dM << "\n";
	cout << "check1" << "\n";
	f[N] = -F;
	cout << "check2" << "\n";

	//To the right of the interface
	double rs_avg_plus = (r[s + 1] + r[s + 2]) / 2.0;
	double rs_diff_plus = rs_avg_plus - inter;
	term1 = (1.0 / rs_diff_plus) * (pow(rs_avg_plus, 2.) * 0.5 * (ALambda[s + 2] + ALambda[s + 1]) * (T_next[s + 2] - T_next[s + 1]) / (r[s + 2] - r[s + 1])
		-pow(inter, 2.) * Lambda(Ti, a, 1) * ((2 * p - 7) / ((p - 3) * (p - 4) * h) * Ti + (p - 4) / ((p - 3) * h) * T_next[s + 2]
			+ (p - 3) / ((4 - p) * h) * T_next[s + 3]));
	term3 = -ACp[s + 1] * (ADensity[s + 1] * pow(r[s + 1], 2) * (T_next[s + 1] - T_cur[s + 1]) / dt + dM * ((T_next[s + 1] - Ti) / (r[s + 1] - inter)));
	cout << "term1= " << term1 << "\n";
	cout << "term3= " << term3 << "\n";
	cout << "q-(s+1)= " << -ALambda[s + 2] * ((2 * p - 7) / ((p - 3) * (p - 4) * h) * Ti + (p - 4) / ((p - 3) * h) * T_next[s + 2]
		+ (p - 3) / ((4 - p) * h) * T_next[s + 3]) << "\n";
	cout << "q+(s+1)= " << pow(rs_avg_plus, 2.) * 0.5 * (ALambda[s + 2] + ALambda[s + 1]) * (T_next[s + 2] - T_next[s + 1]) / (r[s + 2] - r[s + 1]) << "\n";
	F = term1 + term3;
	f[s + 1] = -F;

	//Calculate modul of Nevyazka
	for (int i = 1; i <= N; i++)
		mod_nevyaz += pow(f[i], 2.0);
	mod_nevyaz = pow(mod_nevyaz, 0.5);
}
/*
void ThreeTies(vector<double>& r, vector<double>& T_cur, vector<double>& T_next, vector <double>& ACp, vector <double>& ADensity, vector <double>& ALambda,
	const double a, const double b, const double c, const double dt, const int N, double dM, int s, double inter, double& Ti, const double p, const double T_cur_i)
{
	double h = r[s + 1] - r[s];
	double ACp_i = ACp[N + 1];
	double ADensity_i = ACp[N + 1];
	// Чему будет равна лямбда от Ti? Надо задать в функции
	Ti = (pow(inter, 2) * ((Lambda(Ti, a, b, c, 1) * ((p - 4) / (p - 3) * T_next[s + 2] + (p - 3) / (4 - p) * T_next[s + 3])
		- Lambda(Ti, a, b, c, 0) * (p / (p + 1) * T_next[s - 2] - (p + 1) / p * T_next[s - 1])) / h + ((r[s + 1] - r[s - 1]) / 2) * ACp_i * (ADensity_i * T_cur_i / dt
			+ (dM / pow(inter, 2)) * (T_next[s] / (inter - r[s])))))
		/ (pow(inter, 2) * ((Lambda(Ti, a, b, c, 0) * (2 * p + 1) / ((p + 1) * p) - Lambda(Ti, a, b, c, 1) * (2 * p - 7) / ((p - 3) * (p - 4))) / h
			+ ((r[s + 1] - r[s - 1]) / 2) * ACp_i * (ADensity_i / dt + (dM / pow(inter, 2)) / (inter - r[s]))));
	//new
	Ti = ((pow((r[s + 1] + inter) / 2., 2.) * Lambda(Ti, a, b, c, 1) * ((p - 4) / (p - 3) * T_next[s + 2] + (p - 3) / (4 - p) * T_next[s + 3])
		- pow((r[s] + inter) / 2., 2.) * Lambda(Ti, a, b, c, 0) * (p / (p + 1) * T_next[s - 2] - (p + 1) / p * T_next[s - 1])) / h
		+ ((r[s + 1] - r[s]) / 2) * ACp_i * (dM * (T_next[s] / (inter - r[s]))))
		/ ((pow((r[s] + inter) / 2., 2.) * Lambda(Ti, a, b, c, 0) * (2 * p + 1) / ((p + 1) * p) - pow((r[s + 1] + inter) / 2., 2.)
			* Lambda(Ti, a, b, c, 1) * (2 * p - 7) / ((p - 3) * (p - 4))) / h
			+ ((r[s + 1] - r[s]) / 2) * ACp_i * (dM / (inter - r[s])));
	Ti = (-Lambda(Ti, a, b, c, 1) * ((p - 4) / (p - 3) * T_next[s + 2] + (p - 3) / (4 - p) * T_next[s + 3])
		+ Lambda(Ti, a, b, c, 0) * (p / (p + 1) * T_next[s - 2] - (p + 1) / p * T_next[s - 1]))
		/ (Lambda(Ti, a, b, c, 1) * (2 * p - 7) / ((p - 3) * (p - 4)) - Lambda(Ti, a, b, c, 0) * (2 * p + 1) / ((p + 1) * p));
	//
	T_next[s] = ((pow((r[s] + r[s - 1]) / 2, 2) * 0.5 * (ALambda[s] + ALambda[s - 1]) * (T_next[s - 1] / (r[s] - r[s - 1])))
		+ pow((r[s] + inter) / 2., 2.) * Lambda(Ti, a, b, c, 0) / h * ((p / (p + 1)) * T_next[s - 2] - ((p + 1) / p) * T_next[s - 1] + ((2 * p + 1) / ((p + 1) * p)) * Ti)
		+ ((inter - r[s - 1]) / 2) * ACp[s] * ADensity[s] * pow(r[s], 2) * T_cur[s] / dt)
		/ ((pow((r[s] + r[s - 1]) / 2, 2) * 0.5 * (ALambda[s] + ALambda[s - 1]) / (r[s] - r[s - 1]))
			+ (((inter - r[s - 1]) / 2) * ACp[s] * ADensity[s] * pow(r[s], 2) / dt));
	//
	T_next[s + 1] = ((pow((r[s + 1] + r[s + 2]) / 2., 2.)) * 0.5 * (ALambda[s + 2] + ALambda[s + 1]) * (T_next[s + 2] / (r[s + 2] - r[s + 1]))
		- pow((r[s + 1] + inter) / 2., 2.) * Lambda(Ti, a, b, c, 1) / h * ((2 * p - 7) / ((p - 3) * (p - 4)) * Ti + (p - 4) / (p - 3) * T_next[s + 2] + (p - 3) / (4 - p) * T_next[s + 3])
		+ (((r[s + 2] - inter) / 2) * ACp[s + 1] * (ADensity[s + 1] * pow(r[s + 1], 2.) * T_cur[s + 1] / dt + dM * (Ti / (r[s + 1] - inter))))) /
		(((pow((r[s + 1] + r[s + 2]) / 2., 2.)) * 0.5 * (ALambda[s + 2] + ALambda[s + 1]) / (r[s + 2] - r[s + 1]))
			+ (((r[s + 2] - inter) / 2) * ACp[s + 1] * (ADensity[s + 1] * pow(r[s + 1], 2) / dt + dM / (r[s + 1] - inter))));
	cout << "Lambda_g(T= " << Ti << ")" << Lambda(Ti, a, b, c, 1) << "\n";
	cout << "Lambda_d(T= " << Ti << ")" << Lambda(Ti, a, b, c, 0) << "\n";
}
*/

void Gauss(vector<double>& T_next, vector<double>& alpha, vector<double>& beta, vector<double>& dT,
	vector<double>& nevyaz, vector < vector <double> >& J, double& Ti, const int N_minus, const int s, int& flag_evap)
{
	const int N = N_minus + 1;
	const double T_boiling = 373.0;
	//Cчитаем коэффициенты для всех узлов
	int i, j, k;
	/*
	for (i = 0; i <= N; i++)
	{
		for (j = 0; j <= N; j++)
		{
			cout << "J[ " << i << "," << j << "]= " << J[i][j] << " ";
		}
		cout << "F,[ " << i << "]= " << nevyaz[i] << "\n" << "\n";
	}
	*/
	double d, sum, to_right;
	for (k = 0; k <= N; k++) // прямой ход
	{
		for (j = k + 1; j <= N; j++) //на каждой итерации 3 формулы, чтобы привести к трапециевидной матрице
		{
			d = J[j][k] / J[k][k]; // формула (1): обезразмерим элемент столбца на диагональный элемент в столбце
			for (i = k; i <= N; i++)
			{
				J[j][i] = J[j][i] - d * J[k][i]; // формула (2): 
			}
			nevyaz[j] = nevyaz[j] - d * nevyaz[k]; // формула (3)
		}
	}
	for (k = N; k >= 0; k--) // обратный ход
	{
		sum = 0;
		for (j = k + 1; j <= N; j++)
		{
			to_right = J[k][j] * dT[j]; // Произведение j-ого якобиана на значение j-ого вектора справа от диагонали
			sum = sum + to_right; // сумма произведений якобианов на значение j-ого вектора справа от диагонали
		}
		dT[k] = (nevyaz[k] - sum) / J[k][k]; // формула (4)
	}
	cout << "Roots of system: " << "\n";
	for (i = 0; i <= N_minus; i++)
	{
		T_next[i] += dT[i];
		if (i < s + 1 && T_next[i] > T_boiling)
			T_next[i] = T_boiling;
	}
	Ti += dT[N];
	if (Ti > T_boiling)
		Ti = T_boiling;
	//T_next[20] = T_next[21];
	//T_next[0] = T_next[1];
	//cout << "T[0][k]" << T[0][k] << "\n";
	cout << "s_Gauss= " << s << " T_next[s]" << T_next[s] << "\n";
}

void Progonka(vector<double>& T_next, vector<double>& alpha, vector<double>& beta, vector<double>& dT,
	vector<double>& nevyaz, vector < vector <double> >& J, const int N_minus, const int s, int& flag_evap)
{
	double l_f_num = 0;
	const double T_boiling = 373.0;
	//считаем коэффициенты слева
	//Define denominator for constant formules in method
	double zn;
	zn = J[l_f_num][l_f_num];
	//cout << "zn " << zn << endl;
	alpha[l_f_num] = -J[l_f_num][l_f_num + 1] / zn;
	//cout << "a[1]" << alpha[1] << endl;
	beta[l_f_num] = nevyaz[l_f_num] / zn;
	//cout << "B[1]" << beta[1] << endl;
	//считаем коэффициенты для всех узлов
	for (int i = 1; i < N_minus; i++) {
		//if (i == s + 2) {
			//zn = J[i][i] + J[i][i - 1] * alpha[s - 1];
			//alpha[i] = -J[i][i + 1] / zn;
			//beta[i] = (nevyaz[i] - J[i][i - 1] * beta[s - 1]) / zn;
		//}
		//else {
			zn = J[i][i] + J[i][i - 1] * alpha[i - 1];
			alpha[i] = -J[i][i + 1] / zn;
			beta[i] = (nevyaz[i] - J[i][i - 1] * beta[i - 1]) / zn;
		//}
		//Skip i = s & s + 1
		//if (i == s - 1)
			//i = s + 1;
	}
	dT[N_minus] = (nevyaz[N_minus] - J[N_minus][N_minus - 1] * beta[N_minus - 1]) /
		(J[N_minus][N_minus] + J[N_minus][N_minus - 1] * alpha[N_minus - 1]);
	//cout << "T[N_minus] " << T[N_minus] << endl;
	T_next[N_minus] += dT[N_minus];
	//обновляем температуру на узлах, слева обновляем вручную, а dT остается = 0
	for (int i = N_minus - 1; i >= 0; i--) {
		//if (i == s - 1)
			//dT[i] = alpha[i] * dT[s + 2] + beta[i];
		//else
		dT[i] = alpha[i] * dT[i + 1] + beta[i];
		T_next[i] += dT[i];
		if (i < s + 1 && T_next[i] > T_boiling)
		{
			//dT[i] = 0;
			T_next[i] = T_boiling;
		}
		//Skip i = s & s + 1
		//if (i == s + 2)
			//i = s;
	}
	//T_next[20] = T_next[21];
	//T_next[0] = T_next[1];
	//cout << "T[0][k]" << T[0][k] << "\n";
	cout << "s_Progonka= " << s << " T_next[s]" << T_next[s] << "\n";
}

void Solver(const int N, vector<double>& r, vector<double>& T_cur, vector<double>& T_next, const double d,
	const double a, const double b, const double c, const int total_time, const double dt, const int const_params,
	const double T_l, const double T_r, const double dtau, const int Nd, int m_max)
{
	int m;
	double dM, q_g, q_d, u_r, Lambda_d, Lambda_g, Ti, h;
	const int N_minus = N - 1;
	vector<double> alpha(N), beta(N), dT(N + 1), T_next_upd(N + 1), nevyaz(N + 1);
	vector < vector <double> > J(N + 1, vector <double>(N + 1));
	vector <double> ACp, ADfCp, ADensity, ADfDensity, ALambda, ADfLambda;
	int flag_evap = 0;
	//Define initial surface position and mass flow on surface
	double p = 1.5;
	int s = Nd;
	double inter = r[s] + (p - 1) * (r[s + 1] - r[s]);
	//Define initial boiling tempereature for surface, K
	double T_cur_i = 300.;
	//double q = 20 * pow(10, 3);                  //W / m^2
	double L_d = 2258.2 * pow(10, 3);              //J / kg 
	//const determination of dM
	//dM = 4 * PI * pow(inter, 2.0) * q / L_d;     // kg / s
	//dM equals zero because it is only warming time
	dM = 0.;                                       // kg / s
	//dM / (4 * PI);
	ArraysParameters(a, b, c, N, r, T_cur, s, ACp, ADfCp, ADensity, ADfDensity, ALambda, ADfLambda, 0, inter, T_cur_i, p);
	//Reading the initial distribution
	GraphicsSystEqu(0, r, T_cur, ALambda, ADensity, N, s, dM, flag_evap, inter, T_cur_i);
	//Задаем начальную температуру следующего слоя как решение предыдущей
	Ti = T_cur_i;
	for (int j = 0; j < N + 1; j++)
		T_next[j] = T_cur[j];

	ofstream OutSurfaceFlow;
	ofstream OutLastStepNevyazka;
	ofstream OutFirstStepNevyazka;
	OpeningFiles(OutSurfaceFlow, OutLastStepNevyazka, OutFirstStepNevyazka);
	//TIME LOOP
	ofstream OutNevyazka;
	int divide_count = 0;
	double mod_nevyaz_prev;
	double mod_nevyaz_cur = pow(10., 10.);
	for (int n = 1; n < total_time; n++)
	{
		//Check for flag_evap
		if (Ti >= 373.)
			flag_evap = 1;
		cout << "flag_evap= " << flag_evap << "\n";
		//Write Nevyazka on the last m-iteration
		OutLastStepNevyazka << (n - 1) * dt << " " << mod_nevyaz_cur << "\n";
		//Change Ti depending on the flow didderence and Writing flows, dM = dM/(4*PI)
		u_r = dM / (pow(inter, 2.0) * ADensity[s + 1]);
		h = r[s + 1] - r[s];
		OutSurfaceFlow << dt * (n - 1) << " " << T_next[0] << " " << T_next[s] << " " << Ti << " " << T_next[s + 1] << " "
			<< (p / (p + 1) * T_next[s - 2] - (p + 1) / p * T_next[s - 1] + (2 * p + 1) / ((p + 1) * p) * Ti) / (r[s + 1] - r[s]) << " "
			<< ((2 * p - 7) / ((p - 3) * (p - 4)) * Ti + (p - 4) / (p - 3) * T_next[s + 2] + (p - 3) / (4 - p) * T_next[s + 3]) / (r[s + 1] - r[s]) << " "
			<< Lambda(Ti, a, 0) / h * ((p / (p + 1)) * T_next[s - 2] - ((p + 1) / p) * T_next[s - 1] + ((2 * p + 1) / ((p + 1) * p)) * Ti) << " "
			<< Lambda(Ti, a, 1) / h * ((2 * p - 7) / ((p - 3) * (p - 4)) * Ti + (p - 4) / (p - 3) * T_next[s + 2] + (p - 3) / (4 - p) * T_next[s + 3]) << " "
			<< dM << " " << ADensity[s] << " " << u_r << " " << p << " " << inter << "\n";
		//Opening file for Nevyazka
		OutNevyazka.open("C:/Users/user/source/Научная работа/DropletWarming/Data_new/Nevyazka/Nevyazka_" + to_string(n) + ".dat");
		OutNevyazka << "TITLE=\"" << "Graphics" << "\"" << endl;
		OutNevyazka << R"(VARIABLES= "i", "F")" << endl;
		//Renew counter and Nevyazka
		mod_nevyaz_cur = pow(10., 10.);
		m = 0;
		//Getting Solution using cycle of Newton method, check if isnan() equals (mod == mod)
		while (mod_nevyaz_cur > pow(10., -1) && mod_nevyaz_cur == mod_nevyaz_cur && m < m_max)
		{
			//2nd stage(combination of two processes)
			if (flag_evap == 1)
			{
				//Incrementing the Reference Node for the Liquid / Gazeous Interface
				if (p >= 2.) {
					p = 1.001;
					s += 1;
				}
				//Decrementing the Reference Node for the Liquid / Gazeous Interface
				if (p <= 1.) {
					p = 1.999;
					s -= 1;
					//If the Entire Melt Layer Becomes Gazeous the Interace Boundary Condition is Not Applied
					if (s <= 5)
						cout << "Entire Melt Layer Becomes Gazeous, Time:", dt* n;
				}
				//Renew parameters
				ArraysParameters(a, b, c, N, r, T_next, s, ACp, ADfCp, ADensity, ADfDensity, ALambda, ADfLambda, n, inter, Ti, p);
				Lambda_g = (Lambda(Ti, a, 1) + ALambda[s + 2] + ALambda[s + 3]) / 3.;
				Lambda_d = (ALambda[s - 2] + ALambda[s - 1] + Lambda(Ti, a, 0)) / 3.;
				//Non-const determination of dM
				q_g = Lambda_g * ((2 * p - 7) / ((p - 3) * (p - 4)) * Ti + (p - 4) / (p - 3) * T_next[s + 2] + (p - 3) / (4 - p) * T_next[s + 3]);
				q_d = Lambda_d * (p / (p + 1) * T_next[s - 2] - (p + 1) / p * T_next[s - 1] + (2 * p + 1) / ((p + 1) * p) * Ti);
				dM = -pow(inter, 2.0) / L_d / (r[s + 1] - r[s]) * (q_g - q_d);
				p += dM * dt / (ADensity[s] * pow(inter, 2.0) * (r[s + 1] - r[s]));
				inter = r[s] + (p - 1) * (r[s + 1] - r[s]);
				cout << "p = " << p << "\n";
				cout << "inter = " << inter << "\n";
			}
			//Renew parameters
			ArraysParameters(a, b, c, N, r, T_next, s, ACp, ADfCp, ADensity, ADfDensity, ALambda, ADfLambda, n, inter, Ti, p);
			mod_nevyaz_prev = mod_nevyaz_cur;
			mod_nevyaz_cur = 0;
			F(nevyaz, T_cur, T_next, r, ACp, ADensity, ALambda, N, a, b, c, d, dt, dM, s, inter, Ti, p, L_d, mod_nevyaz_cur, flag_evap);
			cout << "nevyazka = " << mod_nevyaz_cur << "\n" << "\n";
			//Прежде чем пускать Невязку в прогонку, нужно сделать ее меньше, чем в прошлой итерации, но не более чем с 5 попыток
			//На 1ой итерации не зайдет сюда, потому что искусственно зададим огромное T_next_prev
			divide_count = 0;
			//Redefine previous iteration temperature
			while (mod_nevyaz_cur > mod_nevyaz_prev && divide_count < 5)
			{
				for (int i = 0; i < N; i++) {
					T_next[i] -= dT[i];
					dT[i] = dT[i] / 2.;
					//cout << "dT[i]" << dT[i] << "\n";
					T_next[i] += dT[i];
				}
				mod_nevyaz_cur = 0;
				F(nevyaz, T_cur, T_next, r, ACp, ADensity, ALambda, N, a, b, c, d, dt, dM, s, inter, Ti, p, L_d, mod_nevyaz_cur, flag_evap);
				cout << "nevyazka(modif) = " << mod_nevyaz_cur << "\n" << "\n";
				divide_count += 1;
			}
			//Collecting Nevyazka
			OutNevyazka << m << " " << mod_nevyaz_cur << "\n";
			if (m == 1)
				OutFirstStepNevyazka << n * dt << " " << mod_nevyaz_cur << "\n";
			Jacobian(J, r, T_next, ACp, ADfCp, ADensity, ADfDensity, ALambda, ADfLambda, N_minus, a, b, c, d, dt, dtau, dM, s, inter, Ti, p, flag_evap);
			cout << "m " << m << endl;
			//МЕТОД Гаусса: T_next change from m to m+1
			Gauss(T_next, alpha, beta, dT, nevyaz, J, Ti, N_minus, s, flag_evap);
			//Define 3 Temperatures near Interface
			//ThreeTies(r, T_cur, T_next, ACp, ADensity, ALambda, a, b, c, dt, N, dM, s, inter, Ti, p, T_cur_i);
			//Ti = (-ALambda[s + 2] * ((p - 4) / (p - 3) * T_next[s + 2] + (p - 3) / (4 - p) * T_next[s + 3])
			//    + ALambda[s - 1] * (p / (p + 1) * T_next[s - 2] - (p + 1) / p * T_next[s - 1]))
			//    / (ALambda[s + 2] * (2 * p - 7) / ((p - 3) * (p - 4)) - ALambda[s - 1] * (2 * p + 1) / ((p + 1) * p));
			cout << "Ti" << Ti << "\n";
			cout << flag_evap << "\n";
			//cout << "T0" << " =" << T[0] << endl;
			m += 1;
		}
		OutNevyazka.close();
		//Collecting Tempereature, Lambda, Flow and Nevyazka from Timestep
		GraphicsSystEqu(n, r, T_next, ALambda, ADensity, N, s, dM, flag_evap, inter, Ti);
		//задаем начальную температуру следующего слоя как решение предыдущей
		for (int j = 0; j < N + 1; j++)
			T_cur[j] = T_next[j];
		T_cur_i = Ti;
		cout << "s = " << s << "\n";
		if (mod_nevyaz_cur != mod_nevyaz_cur)
			break;
	}
	OutSurfaceFlow.close();
	OutLastStepNevyazka.close();
	cout << "Done!\n";
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
	T_cur.resize(N + 1);
	T_next.resize(N + 1);
	r.resize(N + 1);
	ofstream OutX;
	OutX.open("C:/Users/user/source/Научная работа/DropletWarming/Data_new/X_grid.dat");
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
				cout << "h" << h << " r" << j << " " << r[j] << endl;
			}
			//Non-uniform grid
			for (j; j < N + 1; j++)
			{
				r[j] = r[j - 1] + h;
				h = h * q;
				T_cur[j] = T_r;
				//cout << "T" << j << " " << T_cur[j] << endl;
				OutX << j << " " << h << " " << "\n";
				cout << "h" << h << " r" << j << " " << r[j] << endl;
			}
		}
	}
	OutX.close();
}

int main()
{
	const int const_params = 1;          // 0 - Non-uniform grid, 1 - Mixed
	const double dt = 0.00000001;
	const double dtau = 10000000000;
	int m_max = 100;
	const int total_time = 200000;
	const int power_of = 3;              //const if a < pow(10, -8)
	const double a = pow(10, -(power_of + 3)), b = pow(10, -power_of), c = pow(10, -(power_of - 3));  //a = 1
	const double d = 0.1;                          //коэффициент диффузии
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
	const double x_l = 0.0;
	//const double R = 0.0005 м = 500 мкр
	const double x_r = 0.2; //зададим область = 20 см
	double x_uni = x_l + N_uni * h_min;
	double q = DefineQ(N - N_uni, h_min, x_r, x_uni);
	vector <double> r;
	vector <double> T_cur;
	vector <double> T_next;
	cout << "kek " << "\n";
	InitialGrid(N, x_l, x_r, T_l, T_r, r, T_cur, T_next, h, q, h_min, total_time, const_params, Nd, N_uni, N_uni_near_center);
	cout << "lol " << "\n";
	cout << "N" << N << endl;
	Solver(N, r, T_cur, T_next, d, a, b, c, total_time, dt, const_params, T_l, T_r, dtau, Nd, m_max);
}

//helping parts from solid/liquid interface program
/*
!c .....Tracking the Solid / Liquid Interface using lagrangian Interpolation
if (r.GT.5.AND.r.LT.n_node - 5) then                                          !r > 5 и r < n_node - 5(проверка, находится ли не у границ)         !
	flag1 = 1                                                               !(на границах своя интерполяция будет с их учетом)
	!C .....Calculating the Amount of Solid formed in One Time Step
	10      call growth(T0, k, Latent, dx_s, dt_s, L_oref, L_iref, r, p, Tp, Ta)
	!c .....Mapping the Solid / Liquid Interface Location to the New Grid
	!m        solid_L = ((r - 2) * dx_s + p * dx_s) * L_oref                                !нам не нужно учитывать увеличение длины сетки
	!m        p = (solid_L - (r - 2) * dx_s * L_ref) / (dx_s * L_ref)
	!c .....Incrementing the Reference Node for the Solid / Liquid Interface
	if (p.GE.2.d0) then                                                     !когда p>2 переобозначаем на p = 1
		p = 1.001d0                                                         !так как i доходит до центра ячейки и p нужно ее догнать
		r = r + 1                                                             !то есть p движется от 1 до 2, увеличиваем опорный узел на 1
		!c .....If the Entire Melt Layer Becomes Solid the Interface Boundary Condition is Not Applied
		if (r.GE.n_node - 5) then
			flag1 = 0
			print*, "Entire Melt Layer Becomes Solid, Time:", time_s* Time_ref
			goto 100
			endif
			flag = 1
			endif
			!c .....Decrementing the Reference Node for the Solid / Liquid Interface
			if (p.LE.1.d0) then                                                     !тут p движется от 2 до 1 в обратную сторону
				p = 1.999d0
				r = r - 1
				!c .....If the Entire Melt Layer Becomes Liquid the Interace Boundary Condition is Not Applied
				if (r.LE. 5) then
					flag1 = 0
					print*, "Entire Melt Layer Becomes Liquid, Time:", time_s* Time_ref
					goto 100
					endif
					flag = 1
					end if
					!c .....Calculating the Amount of Solid Between the Reference Nodeand the Solid / Liquid Interface
					Iminus = (p - 1) * dx_s
					!c .....Calculating the Amount of Liquid Between the Reference Node and the Solid / Liquid Interface
					Iplus = (2 - p) * dx_s
					endif

	*/
	// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
	// Debug program: F5 or Debug > Start Debugging menu

	// Tips for Getting Started: 
	//   1. Use the Solution Explorer window to add/manage files
	//   2. Use the Team Explorer window to connect to source control
	//   3. Use the Output window to see build output and other messages
	//   4. Use the Error List window to view errors
	//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
	//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
