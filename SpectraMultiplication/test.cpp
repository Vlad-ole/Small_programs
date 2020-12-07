#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include "Math/Polynomial.h"
#include "Math/Interpolator.h"

#include <algorithm>


//Root cern
#include "TApplication.h"
#include "TH1F.h"
//#include "TROOT.h"
//#include "TStopwatch.h"
//#include "TFile.h"
//#include "TTree.h"
//#include "TThread.h"
//#include "TSystem.h"
//#include "TH1F.h"
//#include "TH2F.h"
//#include "TCanvas.h"
//#include "TGraph.h"
//#include "TColor.h"
//#include "TStyle.h"
#include <TRandom.h>
#include <TRandom2.h>
#include "TROOT.h"
#include "TStyle.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"


using namespace std;

double rand(int N, TRandom2 &rnd)
{
	double res = 0;
	
	for (int i = 0; i < N; i++)
	{
		res += rnd.Exp(1);
	}

	return res;
}

int main(int argc, char *argv[])
{
	TApplication theApp("theApp", 0, 0);//let's add some magic! https://root.cern.ch/phpBB3/viewtopic.php?f=3&t=22972

	int x_max = 35;
	TH1F *h1 = new TH1F("h1", "h1", x_max, 0, x_max);
	TH1F *h2 = new TH1F("h2", "h2", x_max, 0, x_max);
	TH1F *h3 = new TH1F("h3", "h3", x_max, 0, x_max);
	TH1F *h5 = new TH1F("h5", "h5", x_max, 0, x_max);
	TH1F *h8 = new TH1F("h8", "h8", x_max, 0, x_max);
	TH1F *h10 = new TH1F("h10", "h10", x_max, 0, x_max);
	TH1F *h15 = new TH1F("h15", "h15", x_max, 0, x_max);


	TRandom2 rnd(1);
	TRandom2 rnd2(2);

	double Ex1 = 0;
	double Ex2 = 0;
	int N_events = 10000;
	for (int i = 0; i < 10000; i++)
	{
		double val_1 = rand(1, rnd);	
		double val_2 = rand(2, rnd);
		double val_10 = rand(10, rnd);
		
		h1->Fill(val_1);
		h2->Fill(val_2);
		h3->Fill( rand(3, rnd));
		h5->Fill( rand(5, rnd));
		h8->Fill( rand(8, rnd));
		h10->Fill(val_10);
		h15->Fill(rand(15, rnd));

		Ex1 += val_10;
		Ex2 += val_10*val_10;
	}
	Ex1 /= N_events;
	Ex2 /= N_events;
	double Var = Ex2 - Ex1*Ex1;
	cout << "Var =" << Var << endl;

	h1->Draw();
	gPad->SetLogy();
	h1->Scale(1 / h1->Integral());
	h1->GetXaxis()->SetTitle("Ne");
	h1->GetYaxis()->SetTitle("Counts");
	h1->GetXaxis()->SetRangeUser(0,35);
	h1->SetStats(0);

	h2->Draw("same");
	h2->Scale(1 / h2->Integral());
	h2->SetLineColor(kRed);

	h3->Draw("same");
	h3->Scale(1 / h3->Integral());
	h3->SetLineColor(kBlack);

	h5->Draw("same");
	h5->Scale(1 / h5->Integral());
	h5->SetLineColor(kOrange);

	h8->Draw("same");
	h8->Scale(1 / h8->Integral());
	h8->SetLineColor(kGreen);

	h10->Draw("same");
	h10->Scale(1 / h10->Integral());
	h10->SetLineColor(kMagenta);

	h15->Draw("same");
	h15->Scale(1 / h15->Integral());
	h15->SetLineColor(kOrange+3);

	auto legend = new TLegend(0.1, 0.7, 0.2, 0.9);
	//legend->SetHeader("The Legend Title", "C"); // option "C" allows to center the header
	legend->AddEntry(h1, "1e", "l");
	legend->AddEntry(h2, "2e", "l");
	legend->AddEntry(h3, "3e", "l");
	legend->AddEntry(h5, "5e", "l");
	legend->AddEntry(h8, "8e", "l");
	legend->AddEntry(h10, "10e", "l");
	legend->AddEntry(h15, "15e", "l");
	legend->Draw();

	TLine *line = new TLine(5, 0, 5, 0.1);
	line->SetLineColor(kRed);
	line->SetLineWidth(5);
	line->Draw();

	TLine *line2 = new TLine(10, 0, 10, 0.1);
	line2->SetLineColor(kRed);
	line2->SetLineWidth(5);
	line2->Draw();

	
	theApp.Run();
	system("pause");
}