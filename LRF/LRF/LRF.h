#pragma once
#include <vector>
#include <utility> 
#include "StaticStuff.h"

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
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"

class LRF
{
public:
	LRF(double x, double y, double N_PE);
	LRF(LRF A, LRF B);
	~LRF();

	void Generate();
	void Print(bool is_print_MC);
	double GetChi2();
	std::vector<std::vector<double>> &GetLRF();
	std::vector<std::vector<double>> &GetLRF_MC();
	static std::vector< double > X_SiPM;
	static std::vector< double > Y_SiPM;
	void DrawTGraph2D();
	//int NumericalMinimization(const char * minName = "Minuit2", const char *algoName = "", int randomSeed = -1);
	//static double const RosenBrock(const double *xx);

private:
	static int N_instanse;
	
	std::vector<std::vector<double>> lrf;
	std::vector<std::vector<double>> lrf_MC;

	void SetTotalPE(double totalPE);
	void SetXYCoordinates(double x, double y);
	void FillTGraph2D();
	void Set_XY_SiPM();
	

	static TRandom3 rnd;
	TGraph2D* gr2D;
	TH2D *h2D;
	const int matrix_size;
	const double pitch;
	static StaticStuff staticStuff; // constructor runs once, single instance
	double x_center;
	double y_center;
	//double totalPE;
};

