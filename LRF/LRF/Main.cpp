#include <iostream>

#include "LRF.h"
#include "test.h"

//root cern
#include "TApplication.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TCut.h"
#include "TGraph2D.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include <TRandom3.h>
#include "TMath.h"

int Number_of_itterations = 0;

using namespace std;

//LRF *lrf_init_single;
//LRF *lrf_init_double;

//
//double RosenBrock(const double *xx)
//{
//	const Double_t x = xx[0];
//	const Double_t y = xx[1]; 
//	const Double_t cnst = xx[2];
//	const Double_t tmp1 = y - x*x;
//	const Double_t tmp2 = 1 - x;
//	return cnst * tmp1*tmp1 + tmp2*tmp2;
//}

double OneEvent(const double *xx)
{
	Number_of_itterations++;
	//if (Number_of_itterations % 10 == 0) 
	cerr << "Number_of_itterations = " << Number_of_itterations << endl;

	//xx[0] - xx[24] = lrf_MC

	const Double_t x1 = xx[25];
	const Double_t y1 = xx[26];
	const Double_t N1 = xx[27];

	LRF lrf1(x1, y1, N1);


	vector<double> E;//expected LRF values
	for (int i = 2; i <= 6; i++)
	{
		for (int j = 2; j <= 6; j++)
		{
			E.push_back(lrf1.GetLRF()[i][j]);
		}
	}

	double chi = 0;
	for (int i = 0; i < 25; i++)
	{
		chi += pow(xx[i] - E[i], 2.0) / E[i];
	}

	return chi;
}

double TwoEvents(const double *xx)
{
	Number_of_itterations++;
	//if (Number_of_itterations % 10 == 0) 
		cerr << "Number_of_itterations = " << Number_of_itterations << endl;
	
	//xx[0] - xx[24] = lrf_MC
	
	const Double_t x1 = xx[25];
	const Double_t y1 = xx[26];
	const Double_t N1 = xx[27];
	
	const Double_t x2 = xx[28];
	const Double_t y2 = xx[29];
	const Double_t N2 = xx[30];	

	LRF lrf1(x1, y1, N1);
	LRF lrf2(x2, y2, N2);
	LRF lrf_sum(lrf1, lrf2);

	vector<double> E;//expected LRF values
	for (int i = 2; i <= 6; i++)
	{
		for (int j = 2; j <= 6; j++)
		{
			E.push_back(lrf_sum.GetLRF()[i][j]);
		}
	}
	
	double chi = 0;
	for (int i = 0; i < 25; i++)
	{
		chi += pow(xx[i] - E[i], 2.0) / E[i];
	}
	
	
	//for (int i = 2; i <= 6; i++)
	//{
	//	for (int j = 2; j <= 6; j++)
	//	{
	//		chi += pow(lrf_MC[i][j] - lrf[i][j], 2.0) / lrf[i][j];
	//	}
	//}

	return chi;
}

double TwoEventsMLE(const double *xx)
{
	Number_of_itterations++;
	//if (Number_of_itterations % 10 == 0) 
	cerr << "Number_of_itterations = " << Number_of_itterations << endl;

	//xx[0] - xx[24] = lrf_MC

	const Double_t x1 = xx[25];
	const Double_t y1 = xx[26];
	const Double_t N1 = xx[27];

	const Double_t x2 = xx[28];
	const Double_t y2 = xx[29];
	const Double_t N2 = xx[30];

	LRF lrf1(x1, y1, N1);
	LRF lrf2(x2, y2, N2);
	LRF lrf_sum(lrf1, lrf2);

	vector<double> E;//expected LRF values
	for (int i = 2; i <= 6; i++)
	{
		for (int j = 2; j <= 6; j++)
		{
			E.push_back(lrf_sum.GetLRF()[i][j]);
		}
	}

	double chi = 0;
	for (int i = 0; i < 25; i++)
	{
		chi += pow(xx[i] - E[i], 2.0) / E[i];
	}

	return chi;
}


int NumericalMinimization(LRF *lrf, const char * minName = "Minuit2",
	const char *algoName = "",
	int randomSeed = -1)
{
	// create minimizer giving a name and a name (optionally) for the specific
	// algorithm
	// possible choices are:
	//     minName                  algoName
	// Minuit /Minuit2             Migrad, Simplex,Combined,Scan  (default is Migrad)
	//  Minuit2                     Fumili2
	//  Fumili
	//  GSLMultiMin                ConjugateFR, ConjugatePR, BFGS,
	//                              BFGS2, SteepestDescent
	//  GSLMultiFit
	//   GSLSimAn
	//   Genetic
	ROOT::Math::Minimizer* minimum =
		ROOT::Math::Factory::CreateMinimizer(minName, algoName);

	// set tolerance , etc...
	minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
	minimum->SetMaxIterations(10000);  // for GSL
	minimum->SetTolerance(0.001);
	minimum->SetPrintLevel(1);

	// create function wrapper for minimizer
	// a IMultiGenFunction type
	ROOT::Math::Functor f(&TwoEvents, 31);
	//double step[2] = { 0.01, 0.01 };
	// starting point

	//double variable[2] = { -1., 1.2 };
	//if (randomSeed >= 0) {
	//	TRandom2 r(randomSeed);
	//	variable[0] = r.Uniform(-20, 20);
	//	variable[1] = r.Uniform(-20, 20);
	//}

	minimum->SetFunction(f);

	vector<double> variable;
	variable.resize(31);

	vector<double> O;//observed values (LRF_MC)
	for (int i = 2; i <= 6; i++)
	{
		for (int j = 2; j <= 6; j++)
		{
			O.push_back(lrf->GetLRF_MC()[i][j]);
		}
	}
	for (int i = 0; i < 25; i++)
	{
		ostringstream sst;
		sst << i;
		minimum->SetFixedVariable(i, sst.str().c_str(), O[i]);
		variable[i] = O[i];
	}


	minimum->SetLimitedVariable(25, "x1", variable[25], 0.1, -19.9, 0);
	minimum->SetLimitedVariable(26, "y1", variable[26], 0.1, -19.9, 19.9);
	minimum->SetLimitedVariable(27, "N1", variable[27], 1, 1500 , 2500 );
	//minimum->SetFixedVariable(27, "N1", 2000);

	minimum->SetLimitedVariable(28, "x2", variable[28], 0.1, 0, 19.9);
	//minimum->SetFixedVariable(28, "x2", -15);
	
	minimum->SetLimitedVariable(29, "y2", variable[29], 0.1, -19.9, 19.9);
	//minimum->SetFixedVariable(29, "y2", -15);
	
	minimum->SetLimitedVariable(30, "N2", variable[30], 1, 1500 , 2500 );
	//minimum->SetFixedVariable(30, "N2", 2000);

	
	

	// do the minimization
	minimum->Minimize();

	const double *xs = minimum->X();
	const double *xs_err = minimum->Errors();
	cout << "N_itter = " << minimum->NIterations() << "; minName = " << minName << "; algoName = " << algoName << endl;
	std::cout << "Minimum: f(" << xs[25] << "," << xs[26] << "," << xs[27] << "; "
		<< xs[28] << "," << xs[29] << "," << xs[30] << "): " 
		<< minimum->MinValue() << std::endl;

	cout << "List of parameters:" << endl;
	double max_rel_error = 0;
	for (int i = 0; i < 31; i++)
	{
		double previous_max_rel_err = max_rel_error;
		double rel_err = 0;
		if (xs[i]) rel_err = fabs(xs_err[i] / xs[i]);
		cout << "xs[" << i << "] = " << xs[i] << " +- " << xs_err[i] << "; rel_err = " << rel_err << endl;

		if (rel_err > max_rel_error) max_rel_error = rel_err;
	}
	cout << "max_rel_error(%) = " << max_rel_error * 100 << endl;

	//cout << "List of parameters:" << endl;
	//for (int i = 0; i < 31; i++)
	//{
	//	double rel_err = 0;
	//	if (xs[i]) rel_err = xs_err[i] / xs[i];
	//	cout << "xs[" << i << "] = " << xs[i] << " +- " << xs_err[i] << "; rel_err = " << rel_err << endl;
	//}


	

	//// expected minimum is 0
	//if (minimum->MinValue()  < 1.E-4  && f(xs) < 1.E-4)
	//	std::cout << "Minimizer " << minName << " - " << algoName
	//	<< "   converged to the right minimum" << std::endl;
	//else {
	//	std::cout << "Minimizer " << minName << " - " << algoName
	//		<< "   failed to converge !!!" << std::endl;
	//	Error("NumericalMinimization", "fail to converge");
	//}

	return 0;
}

int NumericalMinimizationOneEvent(LRF *lrf, const char * minName = "Minuit2",
	const char *algoName = "",
	int randomSeed = -1)
{
	ROOT::Math::Minimizer* minimum1 =
		ROOT::Math::Factory::CreateMinimizer(minName, algoName);

	// set tolerance , etc...
	minimum1->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
	minimum1->SetMaxIterations(10000);  // for GSL
	minimum1->SetTolerance(0.001);
	minimum1->SetPrintLevel(1);

	// create function wrapper for minimizer
	// a IMultiGenFunction type
	ROOT::Math::Functor f1(&OneEvent, 28);

	minimum1->SetFunction(f1);

	vector<double> variable;
	variable.resize(28);

	vector<double> O;//observed values (LRF_MC)
	for (int i = 2; i <= 6; i++)
	{
		for (int j = 2; j <= 6; j++)
		{
			O.push_back(lrf->GetLRF_MC()[i][j]);
		}
	}
	for (int i = 0; i < 25; i++)
	{
		ostringstream sst;
		sst << i;
		minimum1->SetFixedVariable(i, sst.str().c_str(), O[i]);
		variable[i] = O[i];
	}
	
	minimum1->SetLimitedVariable(25, "x1", variable[25], 0.1, -19.9, 19.9);
	minimum1->SetLimitedVariable(26, "y1", variable[26], 0.1, -19.9, 19.9);
	minimum1->SetLimitedVariable(27, "N1", variable[27], 1, 1500, 2500);
	//minimum->SetFixedVariable(27, "N1", 2000);

	// do the minimization
	minimum1->Minimize();

	const double *xs = minimum1->X();
	const double *xs_err = minimum1->Errors();
	cout << "N_itter = " << minimum1->NIterations() << "; minName = " << minName << "; algoName = " << algoName << endl;
	std::cout << "Minimum: f(" << xs[25] << "," << xs[26] << "," << xs[27] << "): " << minimum1->MinValue() << std::endl;


	cout << "List of parameters:" << endl;
	double max_rel_error = 0;
	for (int i = 0; i < 28; i++)
	{
		double previous_max_rel_err = max_rel_error;
		double rel_err = 0;
		if (xs[i]) rel_err = fabs(xs_err[i] / xs[i]);
		cout << "xs[" << i << "] = " << xs[i] << " +- " << xs_err[i] << "; rel_err = " << rel_err << endl;

		if (rel_err > max_rel_error) max_rel_error = rel_err;
	}
	cout << "max_rel_error(%) = " << max_rel_error*100 << endl;

	return 0;
}

int main()
{
	TApplication theApp("theApp", 0, 0);

	double x1 = -15;
	double y1 = -15;
	double N1 = 2000;

	double x2 = 15;
	double y2 = 15;
	double N2 = 2000;

	double r = sqrt(pow((x1 - x2), 2.0) + pow((y1 - y2), 2.0));
	cout << "r = " << r << endl;

	cout << "**********************************************" << endl;
	cout << "lrf1" << endl;
	LRF lrf1(x1, y1, N1);
	lrf1.Print(0);
	lrf1.Generate();
	lrf1.Print(1);
	//lrf1.FillTGraph2D();
	//lrf1.DrawTGraph2D();
	cout << "**********************************************" << endl;

	cout << "**********************************************" << endl;
	cout << "lrf2" << endl;
	LRF lrf2(x2, y2, N2);
	lrf2.Print(0);
	lrf2.Generate();
	lrf2.Print(1);
	cout << "**********************************************" << endl;
	
	cout << "**********************************************" << endl;
	cout << "lrf_sum" << endl;
	LRF* lrf_sum = new LRF(lrf1, lrf2);
	lrf_sum->Print(0);
	lrf_sum->Print(1);
	cout << "**********************************************" << endl;


	LRF* lrf_1ev = new LRF(10, 10, 2000);
	lrf_1ev->Generate();

	// create minimizer giving a name and a name (optionally) for the specific
	// algorithm
	// possible choices are:
	//     minName                  algoName
	// Minuit /Minuit2             Migrad, Simplex,Combined,Scan  (default is Migrad)
	//  Minuit2                     Fumili2
	//  Fumili
	//  GSLMultiMin                ConjugateFR, ConjugatePR, BFGS,
	//                              BFGS2, SteepestDescent
	//  GSLMultiFit
	//   GSLSimAn
	//   Genetic

	NumericalMinimization(lrf_sum);
	//NumericalMinimizationOneEvent(lrf_1ev);
	
	//lrf_init_double = new LRF(0, 0, 2000);

	//LRF lrf1(0, 0, 2000);

	//NumericalMinimization();
	
	//LRF lrf_exp;
	//lrf_exp.SetTotalPE(1500);
	//lrf_exp.Print(0);
	
	//lrf_exp.Print(1);

	//lrf_exp.Generate();
	//lrf_exp.Print(1);
	
	//cout << "Chi2perNDF = " << lrf_exp.GetChi2perNDF() << endl;

	

	

	//cout << endl;
	//cout << "Chi2 = " << lrf_sum.GetChi2() << endl;
	//cout << "Chi2_critical (23 ndf) and at alpha=0.01 = " << 41.6 << endl;
	//
	
	//cout << "Prob = " << TMath::Prob(lrf_sum.GetChi2perNDF()*23, 23) << endl;

	//const int N_ch = 9*9;
	//TGraph2D* gr2D = new TGraph2D(N_ch, &x_position[0], &y_position[0], &n_pe_corrected[0]);
	
	//TH1F* h1 = new TH1F("Chi2perNDF", "Chi2perNDF", 100, 0, 5);
	//TH1F* h1_prob = new TH1F("h1_prob", "h1_prob", 100, 0, 1);
	//vector<double> Chi2perNDF;
	//vector<double> p;
	//
	/*for (int i = 0; i < 10000; i++)
	{
		lrf_sum.Generate();
		double chi_per_NDF = lrf_sum.GetChi2perNDF();
		double prob = TMath::Prob(chi_per_NDF * 23, 23);

		h1->Fill(chi_per_NDF);
		h1_prob->Fill(prob);

		Chi2perNDF.push_back(chi_per_NDF);
		p.push_back(prob);
	}*/

	/*TGraph* gr_Chi2NDF_vs_p = new TGraph(Chi2perNDF.size(), &Chi2perNDF[0], &p[0]);

	TCanvas *c1 = new TCanvas("c1", "");
	c1->Divide(2, 2, 0.01, 0.01);
	c1->cd(1);
	h1->Draw();
	h1->GetXaxis()->SetTitle("Chi2perNDF");
	
	c1->cd(2);
	h1_prob->Draw();

	c1->cd(3);
	gr_Chi2NDF_vs_p->Draw("AP");*/
	
	
	
	cout << endl << "All is ok" << endl;
	cerr << endl << "All is ok" << endl;
	theApp.Run();

	//system("pause");
	return 0;
}