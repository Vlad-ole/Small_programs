#include "LRF.h"
#include <iostream>


using namespace std;

TRandom3 LRF::rnd;
StaticStuff LRF::staticStuff;
int LRF::N_instanse = 0;


LRF::LRF(double x, double y, double N_PE) : matrix_size(9), pitch(10)
{
	N_instanse++;

	lrf.resize(matrix_size);
	lrf_MC.resize(matrix_size);
	for (int j = 0; j < matrix_size; j++)
	{
		lrf[j].resize(matrix_size);
		lrf_MC[j].resize(matrix_size);
	}

	

	/* experimental matrix (x-ray, 2mm coll) 171123 20kV
	0.014	0.016	0.022	0.015	0.016
	0.022	0.029	0.053	0.044	0.019
	0.027	0.071	0.183	0.07	0.024
	0.025	0.054	0.081	0.05	0.027
	0.022	0.027	0.032	0.028	0.03
	*/

	/*lrf[0][0] = 0.014;
	lrf[0][1] = 0.016;
	lrf[0][2] = 0.022;
	lrf[0][3] = 0.015;
	lrf[0][4] = 0.016;

	lrf[1][0] = 0.022;
	lrf[1][1] = 0.029;
	lrf[1][2] = 0.053;
	lrf[1][3] = 0.044;
	lrf[1][4] = 0.019;

	lrf[2][0] = 0.027;
	lrf[2][1] = 0.071;
	lrf[2][2] = 0.183;
	lrf[2][3] = 0.07;
	lrf[2][4] = 0.024;

	lrf[3][0] = 0.025;
	lrf[3][1] = 0.054;
	lrf[3][2] = 0.081;
	lrf[3][3] = 0.05;
	lrf[3][4] = 0.027;

	lrf[4][0] = 0.022;
	lrf[4][1] = 0.027;
	lrf[4][2] = 0.032;
	lrf[4][3] = 0.028;
	lrf[4][4] = 0.03;*/

	//extrapolated matrix (x-ray, 2mm coll) 171123 20kV
	/*
	4.25E-03	4.84E-03	5.42E-03	6.01E-03	6.60E-03	6.01E-03	5.42E-03	4.84E-03	4.25E-03
	4.84E-03	8.50E-03	1.00E-02	1.15E-02	1.30E-02	1.15E-02	1.00E-02	8.50E-03	4.84E-03
	5.42E-03	1.00E-02	0.014	0.016	0.022	0.015	0.016	1.00E-02	5.42E-03
	6.01E-03	1.15E-02	0.022	0.029	0.053	0.044	0.019	1.15E-02	6.01E-03
	6.60E-03	1.30E-02	0.027	0.071	0.183	0.07	0.024	1.30E-02	6.60E-03
	6.01E-03	1.15E-02	0.025	0.054	0.081	0.05	0.027	1.15E-02	6.01E-03
	5.42E-03	1.00E-02	0.022	0.027	0.032	0.028	0.017	1.00E-02	5.42E-03
	4.84E-03	8.50E-03	1.00E-02	1.15E-02	1.30E-02	1.15E-02	1.00E-02	8.50E-03	4.84E-03
	4.25E-03	4.84E-03	5.42E-03	6.01E-03	6.60E-03	6.01E-03	5.42E-03	4.84E-03	4.25E-03
	*/
	
	lrf[0] = { 4.25E-03,	4.84E-03,	5.42E-03,	6.01E-03,	6.60E-03,	6.01E-03,	5.42E-03,	4.84E-03,	4.25E-03 };
	lrf[1] = { 4.84E-03,	8.50E-03,	1.00E-02,	1.15E-02,	1.30E-02,	1.15E-02,	1.00E-02,	8.50E-03,	4.84E-03 };
	lrf[2] = { 5.42E-03,	1.00E-02,	0.014,	0.016,	0.022,	0.015,	0.016,	1.00E-02,	5.42E-03 };
	lrf[3] = { 6.01E-03,	1.15E-02,	0.022,	0.029,	0.053,	0.044,	0.019,	1.15E-02,	6.01E-03 };
	lrf[4] = { 6.60E-03, 1.30E-02, 0.027, 0.071, 0.183, 0.07, 0.024, 1.30E-02, 6.60E-03 };
	lrf[5] = { 6.01E-03,	1.15E-02,	0.025,	0.054,	0.081,	0.05,	0.027,	1.15E-02,	6.01E-03 };
	lrf[6] = { 5.42E-03,	1.00E-02,	0.022,	0.027,	0.032,	0.028,	0.017,	1.00E-02,	5.42E-03 };
	lrf[7] = { 4.84E-03,	8.50E-03,	1.00E-02,	1.15E-02,	1.30E-02,	1.15E-02,	1.00E-02,	8.50E-03,	4.84E-03 };
	lrf[8] = { 4.25E-03,	4.84E-03,	5.42E-03,	6.01E-03,	6.60E-03,	6.01E-03,	5.42E-03,	4.84E-03,	4.25E-03 };
	
	//h2D = new TH2D("h2D", "h2D", 1000, -55, 55, 1000, -55, 55);
	//h2D = gr2D->GetHistogram();
	//h2D->Interpolate(x_tmp, y_tmp);

	//Print(0);
	SetTotalPE(N_PE);
	//Print(0);
	FillTGraph2D();
	//Print(0);
	SetXYCoordinates(x,y);
	//DrawTGraph2D();
	//Print(0);	
}

LRF::LRF(LRF A, LRF B) : matrix_size(9), pitch(10)
{
	N_instanse++;
	
	lrf.resize(matrix_size);
	lrf_MC.resize(matrix_size);
	for (int j = 0; j < matrix_size; j++)
	{
		lrf[j].resize(matrix_size);
		lrf_MC[j].resize(matrix_size);
	}
	
	for (int i = 0; i < matrix_size; i++)
	{
		for (int j = 0; j < matrix_size; j++)
		{
			lrf[i][j] = A.GetLRF()[i][j] + B.GetLRF()[i][j];
			lrf_MC[i][j] = A.GetLRF_MC()[i][j] + B.GetLRF_MC()[i][j];
		}
	}

	FillTGraph2D();
}

void LRF::Print(bool is_print_MC)
{
	cout << endl;
	if (is_print_MC) cout << "lrf_MC" << endl;
	else cout << "lrf" << endl;

	for (int i = 0; i < matrix_size; i++)
	{
		for (int j = 0; j < matrix_size; j++)
		{
			if (is_print_MC) cout << lrf_MC[i][j];
			else cout << lrf[i][j];
			
			if (j != matrix_size - 1) cout << "\t";
		}
		cout << endl;
	}
}

void LRF::SetTotalPE(double totalPE)
{
	for (int i = 0; i < matrix_size; i++)
	{
		for (int j = 0; j < matrix_size; j++)
		{
			lrf[i][j] *= totalPE;
		}
	}
}

void LRF::Generate()
{
	for (int i = 0; i < matrix_size; i++)
	{
		for (int j = 0; j < matrix_size; j++)
		{
			lrf_MC[i][j] = rnd.Poisson(lrf[i][j]);
		}
	}
}

double LRF::GetChi2()
{
	double chi = 0;
	double avr = 0;
	
	for (int i = 2; i <= 6; i++)
	{
		for (int j = 2; j <= 6; j++)
		{
			chi += pow(lrf_MC[i][j] - lrf[i][j], 2.0) / lrf[i][j];
		}
	}
	
	int ndf_estimated = 1;
	//int ndf_chi2 = 5*5 - ndf_estimated - 1;
	return chi;
}

std::vector<std::vector<double>> &LRF::GetLRF()
{
	return lrf;
}

std::vector<std::vector<double>> &LRF::GetLRF_MC()
{
	return lrf_MC;
}

void LRF::DrawTGraph2D()
{	
	TH2D *h = new TH2D("h", "", 1000, -50, 50, 1000, -50, 50);
	gr2D->SetHistogram(h);
	gr2D->Draw("TRI1");
}

void LRF::FillTGraph2D()
{
	vector<double> n_pe;
	//n_pe.resize(staticStuff.X_SiPM.size());
	for (int i = 0; i < matrix_size; i++)
	{
		for (int j = 0; j < matrix_size; j++)
		{
			n_pe.push_back(lrf[i][j]);
		}
	}

	double delta = 0;//https://root-forum.cern.ch/t/tgraph2d-interpolation-problem/17807
	for (int i = 0; i < staticStuff.X_SiPM.size(); i++)
	{
		staticStuff.X_SiPM[i] += (rnd.Uniform() - 0.5) * 2 * delta;
		//cout << staticStuff.X_SiPM[i] << "\t" << staticStuff.Y_SiPM[i] << "\t" << n_pe[i] << endl;
	}

	//cout << "&gr2D = " <<  &gr2D << endl;
	gr2D = new TGraph2D(staticStuff.X_SiPM.size(), &staticStuff.X_SiPM[0], &staticStuff.Y_SiPM[0], &n_pe[0]);
	ostringstream sst;
	sst << "gr2D_" << N_instanse;
	gr2D->SetName(sst.str().c_str());

	
	if (N_instanse > 1)
	{
		ostringstream sst_h2D_old;
		sst_h2D_old << "h2D_" << N_instanse - 1;
		delete gROOT->FindObject(sst_h2D_old.str().c_str());
	}
	ostringstream sst_h2D_new;
	sst_h2D_new << "h2D_" << N_instanse;
	h2D = new TH2D(sst_h2D_new.str().c_str(), sst_h2D_new.str().c_str(), 5000, -40, 40, 5000, -40, 40);
	
	//cout << "&gr2D = " << &gr2D << endl;


	h2D = gr2D->GetHistogram();
}

void LRF::SetXYCoordinates(double x, double y)
{
	//fill lrf[i][j] using h2D->Interpolate
	double delta = 0;//https://root-forum.cern.ch/t/tgraph2d-interpolation-problem/17807
	x_center = x;
	y_center = y;
	//cout << endl;
	int k = 0;
	for (int i = 0; i < matrix_size; i++)
	{
		for (int j = 0; j < matrix_size; j++)
		{			
			if (i >= 2 && i <= 6 && j >= 2 && j <= 6)
			{				
				lrf[i][j] = h2D->Interpolate(staticStuff.X_SiPM[k] - x_center + delta, staticStuff.Y_SiPM[k] - y_center + delta);
				//cout << staticStuff.X_SiPM[k] - x_center << "\t" << staticStuff.Y_SiPM[k] - y_center << "\t" << h2D->Interpolate(0, 0) << endl;
				//cout << gr2D->Interpolate(0,0) << endl;
			}
			k++;
		}
	}

	//update h2D, using lrf[i][j]
	vector<double> n_pe;
	for (int i = 0; i < matrix_size; i++)
	{
		for (int j = 0; j < matrix_size; j++)
		{
			n_pe.push_back(lrf[i][j]);
		}
	}
	if (gr2D != NULL)
	{
		delete gr2D;
		gr2D = new TGraph2D(staticStuff.X_SiPM.size(), &staticStuff.X_SiPM[0], &staticStuff.Y_SiPM[0], &n_pe[0]);
		ostringstream sst;
		sst << "gr2D_" << N_instanse << "SetXYCoordinates";
		gr2D->SetName(sst.str().c_str());
	}


	//cout << "gr2D->GetN() = " << gr2D->GetN() << endl;
	//for (int i = 0; i < gr2D->GetN(); i++)
	//{
	//	cout << i << "\t" << gr2D->GetPo << endl;
	//}

	//gr2D = new TGraph2D(staticStuff.X_SiPM.size(), &staticStuff.X_SiPM[0], &staticStuff.Y_SiPM[0], &n_pe[0]);
	h2D = gr2D->GetHistogram();
}



//int LRF::NumericalMinimization(const char * minName, const char *algoName, int randomSeed)
//{
//	// create minimizer giving a name and a name (optionally) for the specific
//	// algorithm
//	// possible choices are:
//	//     minName                  algoName
//	// Minuit /Minuit2             Migrad, Simplex,Combined,Scan  (default is Migrad)
//	//  Minuit2                     Fumili2
//	//  Fumili
//	//  GSLMultiMin                ConjugateFR, ConjugatePR, BFGS,
//	//                              BFGS2, SteepestDescent
//	//  GSLMultiFit
//	//   GSLSimAn
//	//   Genetic
//	ROOT::Math::Minimizer* minimum =
//		ROOT::Math::Factory::CreateMinimizer(minName, algoName);
//
//	// set tolerance , etc...
//	minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
//	minimum->SetMaxIterations(10000);  // for GSL
//	minimum->SetTolerance(0.001);
//	minimum->SetPrintLevel(1);
//
//	// create function wrapper for minimizer
//	// a IMultiGenFunction type
//	ROOT::Math::Functor f(&RosenBrock, 2);
//	double step[2] = { 0.01, 0.01 };
//	// starting point
//
//	double variable[2] = { -1., 1.2 };
//	if (randomSeed >= 0) {
//		TRandom2 r(randomSeed);
//		variable[0] = r.Uniform(-20, 20);
//		variable[1] = r.Uniform(-20, 20);
//	}
//
//	minimum->SetFunction(f);
//
//	// Set the free variables to be minimized !
//	minimum->SetVariable(0, "x", variable[0], step[0]);
//	minimum->SetVariable(1, "y", variable[1], step[1]);
//
//	// do the minimization
//	minimum->Minimize();
//
//	const double *xs = minimum->X();
//	std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): "
//		<< minimum->MinValue() << std::endl;
//
//	// expected minimum is 0
//	if (minimum->MinValue()  < 1.E-4  && f(xs) < 1.E-4)
//		std::cout << "Minimizer " << minName << " - " << algoName
//		<< "   converged to the right minimum" << std::endl;
//	else {
//		std::cout << "Minimizer " << minName << " - " << algoName
//			<< "   failed to converge !!!" << std::endl;
//		Error("NumericalMinimization", "fail to converge");
//	}
//
//	return 0;
//}
//


LRF::~LRF()
{
	//cout << "&gr2D = " << &gr2D << endl;
	
	// /*if (h2D != NULL)*/ delete h2D;
	// /*if (gr2D != NULL)*/ delete gr2D;
	
}

