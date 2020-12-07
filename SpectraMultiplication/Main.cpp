#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include "Math/Polynomial.h"
#include "Math/Interpolator.h"

#include <algorithm>

//Root cern
#include "TApplication.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TTree.h"
#include "TThread.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TColor.h"
#include "TStyle.h"
#include "TLegend.h"


using namespace std;

int main(int argc, char *argv[])
{
	TApplication theApp("theApp", &argc, argv);//let's add some magic! https://root.cern.ch/phpBB3/viewtopic.php?f=3&t=22972
	
	bool is_PMT = 0;

	double x, y;
	string common_path = "D:\\git_repositories\\Small_programs\\SpectraMultiplication\\data\\";
	ifstream file_data_LAr_scint(common_path + "LAr_scint_Heindl_2011.txt");
	vector<double> nm_LAr_scint;
	vector<double> amp_LAr_scint;	
	while (file_data_LAr_scint.good())
	{
		file_data_LAr_scint >> x >> y;
		nm_LAr_scint.push_back(x);
		amp_LAr_scint.push_back(y);
	}
	TGraph *gr_LAr_scint = new TGraph(nm_LAr_scint.size(), &nm_LAr_scint[0], &amp_LAr_scint[0]);
	gr_LAr_scint->Draw();
	gr_LAr_scint->GetYaxis()->SetRangeUser(0, 1.2);
	gr_LAr_scint->SetTitle("");
	gr_LAr_scint->GetXaxis()->SetTitle("Wavelength [nm]");
	gr_LAr_scint->SetName("gr_LAr_scint");



	ifstream file_data_PDE_48V(common_path + "PDE_48V.txt");
	vector<double> nm_PDE_48V;
	vector<double> amp_PDE_48V;
	while (file_data_PDE_48V.good())
	{
		file_data_PDE_48V >> x >> y;
		nm_PDE_48V.push_back(x);
		amp_PDE_48V.push_back(y);
	}
	TGraph *gr_PDE_48V = new TGraph(nm_PDE_48V.size(), &nm_PDE_48V[0], &amp_PDE_48V[0]);
	gr_PDE_48V->Draw("same");
	gr_PDE_48V->SetLineColor(kGreen);
	gr_PDE_48V->SetName("gr_PDE_48V");


	ifstream file_data_PMT(common_path + "PMT.txt");
	vector<double> nm_PMT;
	vector<double> amp_PMT;
	while (file_data_PMT.good())
	{
		file_data_PMT >> x >> y;
		nm_PMT.push_back(x);
		amp_PMT.push_back(y);
	}
	TGraph *gr_PMT = new TGraph(nm_PMT.size(), &nm_PMT[0], &amp_PMT[0]);
	gr_PMT->Draw("same");
	gr_PMT->SetLineColor(kGreen+2);
	gr_PMT->SetName("gr_PMT");



	ifstream file_data_NBrS_4_6Td(common_path + "NBrS_4_6Td.txt");
	vector<double> nm_NBrS_4_6Td;
	vector<double> amp_NBrS_4_6Td;
	while (file_data_NBrS_4_6Td.good())
	{
		file_data_NBrS_4_6Td >> x >> y;
		nm_NBrS_4_6Td.push_back(x);
		amp_NBrS_4_6Td.push_back(y);
	}
	TGraph *gr_NBrS_4_6Td = new TGraph(nm_NBrS_4_6Td.size(), &nm_NBrS_4_6Td[0], &amp_NBrS_4_6Td[0]);
	gr_NBrS_4_6Td->Draw("same");
	gr_NBrS_4_6Td->SetLineColor(kBlue);
	gr_NBrS_4_6Td->SetName("gr_NBrS_4_6Td");


	ifstream file_data_T_acrylic_1_4mm(common_path + "T_acrylic_1_4mm.txt");
	vector<double> nm_T_acrylic_1_4mm;
	vector<double> amp_T_acrylic_1_4mm;
	while (file_data_T_acrylic_1_4mm.good())
	{
		file_data_T_acrylic_1_4mm >> x >> y;
		nm_T_acrylic_1_4mm.push_back(x);
		amp_T_acrylic_1_4mm.push_back(y);
	}
	TGraph *gr_T_acrylic_1_4mm = new TGraph(nm_T_acrylic_1_4mm.size(), &nm_T_acrylic_1_4mm[0], &amp_T_acrylic_1_4mm[0]);
	gr_T_acrylic_1_4mm->Draw("same");
	gr_T_acrylic_1_4mm->SetLineColor(kRed);
	gr_T_acrylic_1_4mm->SetName("gr_T_acrylic_1_4mm");


	vector<double> nm_mult;
	vector<double> amp_mult;	
	for (int i = 300; i < 1000; i++)
	{	

		if (is_PMT)
		{
			y = gr_PMT->Eval(i) * gr_T_acrylic_1_4mm->Eval(i) * gr_LAr_scint->Eval(i);
		}
		else
		{
			y = gr_PDE_48V->Eval(i) * gr_T_acrylic_1_4mm->Eval(i) * gr_LAr_scint->Eval(i);
		}

		nm_mult.push_back(i);
		amp_mult.push_back(y);
	}
	TGraph *gr_mult = new TGraph(nm_mult.size(), &nm_mult[0], &amp_mult[0]);
	gr_mult->Draw("same");
	gr_mult->SetLineColor(kMagenta);
	gr_mult->SetName("gr_mult");
	cout << "Integral = " << gr_mult->Integral() /*sum*/ << endl;

	auto legend = new TLegend(0.1, 0.7, 0.48, 0.9);
	//legend->SetHeader("The Legend Title", "C"); // option "C" allows to center the header
	legend->AddEntry("gr_LAr_scint", "S1 spectrum (Heindl et al.)", "l");
	legend->AddEntry("gr_PDE_48V", "SiPM PDE at 48V", "l");
	legend->AddEntry("gr_PMT", "PMT QE (T corrected)", "l");
	legend->AddEntry("gr_NBrS_4_6Td", "S2 NBrS spectrum at 4.6 Td", "l");
	legend->AddEntry("gr_T_acrylic_1_4mm", "Acrylic plate (1.5 mm) transmittance", "l");

	if (is_PMT) legend->AddEntry("gr_mult", "PMT_QE * Acrylic_Trans * Emiss_spectrum_Heindl", "l");
	if (!is_PMT) legend->AddEntry("gr_mult", "SiPM_PDE * Acrylic_Trans * Emiss_spectrum_Heindl", "l");
	legend->Draw();

	//x-ray
	//const double Y_PE_per_e = 35.6/7.9E6 /*1.04E-6*/;
	//const double dOmega = 2.15E-3 /*3.30E-4*/; 
	//const double W = 23.6; // eV/e-
	//double Y_max = Y_PE_per_e / (dOmega * gr_mult->Integral());
	//cout << "Y_max [ph/(e*nm)] = " << Y_max << endl;
	//Y_max *= (1E6 / W);

	//alpha
	double Y_PE_per_MeV;
	double dOmega;
	double Y_max;

	if (is_PMT)
	{
		Y_PE_per_MeV = 0.66 / 5.5;//4PMT
		dOmega = 4 * 0.5E-2;
		Y_max = Y_PE_per_MeV / (dOmega * gr_mult->Integral());
	}
	else
	{
		Y_PE_per_MeV = 0.415 / 5.5; //23_SiPM
		dOmega = 2.21E-3;
		Y_max = Y_PE_per_MeV / (dOmega * gr_mult->Integral() * 0.50 / 0.55); //for SiPMs at 46V
	}

	cout << "Y_max [ph/(MeV*nm)] = " << Y_max << endl;


	double sum = 0;	
	for (int i = 400; i < 1000; i++)
	{
		sum += Y_max * /*gr_NBrS_4_6Td->Eval(i)*/ gr_LAr_scint->Eval(i);
	}
	cout << "Absolute yield [ph/MeV] = " << sum << endl;



	cout << endl;
	cout << "all is ok" << endl;
	theApp.Run();
	system("pause");
	return 0;
}