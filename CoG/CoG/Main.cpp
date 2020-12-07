#include <vector>
#include <iostream>

#include <TRandom.h>
#include <TRandom2.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TGraph.h>
#include "TApplication.h"
#include "TROOT.h"
#include "TStyle.h"

using namespace std;

int main(int argc, char *argv[])
{
	TApplication theApp("theApp", &argc, argv);//let's add some magic! https://root.cern.ch/phpBB3/viewtopic.php?f=3&t=22972

	double x_min = -23;
	double x_max = 23;
	int n_bins = 9;
	double bin_step = (x_max - x_min) / n_bins;

	TH1F *h_SiPM_mu_0 = new TH1F("h_SiPM_mu_0", "h_SiPM_mu_0", n_bins, x_min, x_max);
	TH1F *h_SiPM_mu_non0 = new TH1F("h_SiPM_mu_non0", "h_SiPM_mu_non0", n_bins, x_min, x_max);
	TH1F *h_original_mu_non0 = new TH1F("h_original_mu_non0", "h_original_mu_non0", 100, -25, 25);
	TH1F *h_original_mu_0 = new TH1F("h_original_mu_0", "h_original_mu_0", 100, -25, 25);
	TH1F *h_cog_raw_mu0 = new TH1F("h_cog_raw_mu0", "h_cog_raw_mu0", 500, x_min, x_max);

	TRandom2 rnd;
	TRandom2 rnd_u(100);

	cout << "SiPM width = " << h_SiPM_mu_0->GetBinWidth(1) << endl;
	cout << "Channel pitch = " << h_SiPM_mu_0->GetBinWidth(1) << endl;
	for (int i = 1; i <= n_bins; i++)
	{
		cout << i << "\t" << h_SiPM_mu_0->GetBinLowEdge(i) << "\t" << (h_SiPM_mu_0->GetBinLowEdge(i) + h_SiPM_mu_0->GetBinWidth(i)) << endl;
	}

	int N_events = 10;
	vector< vector<int> > PE_in_SiPM_matrix;
	vector<double> x_real;
	vector<double> x_cog_raw;
	double width = 26;
	double gauss_sigma = 5;
	double wide_gauss_sigma = 25;
	double wide_gauss_amp = 1;
	int Npe = 50;

	for (int i = 0; i < N_events; i++)
	{
		//if (i % 10000 == 0)
			//cout << "event = " << i << endl;
		

		for (double val_u = -20; val_u <= 20; val_u++)
		{
			x_real.push_back(val_u);

			double x_cog_tmp = 0;
			int k_repeats = 50;
			for (int k = 0; k < k_repeats; k++)
			{

				//double val_u = rnd_u.Uniform(-width / 2.0, width / 2.0);
				vector<int> PE_in_SiPM_matrix_tmp(n_bins);

				//int Npe_in_event = 10;
				int Npe_in_event_counter = 0;
				double val = rnd.Gaus(val_u, gauss_sigma) + wide_gauss_amp*rnd.Gaus(val_u, wide_gauss_sigma);
				h_original_mu_non0->Fill(val);
				
				while (Npe_in_event_counter < Npe)
				{
					val = rnd.Gaus(val_u, gauss_sigma) + wide_gauss_amp*rnd.Gaus(val_u, wide_gauss_sigma);
					
					for (int j = 1; j <= n_bins; j++)
					{
						if (j % 2 == 1)
						{
							if ((val >(x_min + bin_step*(j - 1))) && (val < x_min + bin_step*j))
							{
								h_SiPM_mu_non0->Fill(val);
								PE_in_SiPM_matrix_tmp[j - 1]++;
								Npe_in_event_counter++;
							}

						}
					}
				}
				
				double x_cog = PE_in_SiPM_matrix_tmp[0] * -20 + PE_in_SiPM_matrix_tmp[2] * -10 +
					PE_in_SiPM_matrix_tmp[4] * 0 + PE_in_SiPM_matrix_tmp[6] * 10 + PE_in_SiPM_matrix_tmp[8] * 20;

				x_cog /= (PE_in_SiPM_matrix_tmp[0] + PE_in_SiPM_matrix_tmp[2] + PE_in_SiPM_matrix_tmp[4]
					+ PE_in_SiPM_matrix_tmp[6] + PE_in_SiPM_matrix_tmp[8]);

				x_cog_tmp += x_cog;

				

			}//k

			x_cog_tmp /= k_repeats;
			x_cog_raw.push_back(x_cog_tmp);

			//PE_in_SiPM_matrix.push_back(PE_in_SiPM_matrix_tmp);
		}//val_u
	}//events



	//mu = 0;
	for (int i = 0; i < 100000; i++)
	{
		//if (i % 10000 == 0)
			//cout << "event = " << i << endl;
		
		
		
		double val_mu_0 = rnd.Gaus(0, gauss_sigma) + wide_gauss_amp*rnd.Gaus(0, wide_gauss_sigma);
		h_original_mu_0->Fill(val_mu_0);
		vector<int> PE_in_SiPM_matrix_tmp(n_bins);
		
		//int Npe_in_event = 10;
		int Npe_in_event_counter = 0;
		while (Npe_in_event_counter < Npe)
		{
			Npe_in_event_counter++;
			val_mu_0 = rnd.Gaus(0, gauss_sigma) + wide_gauss_amp*rnd.Gaus(0, wide_gauss_sigma);
			for (int j = 1; j <= n_bins; j++)
			{
				if (j % 2 == 1)
				{
					if ((val_mu_0 >(x_min + bin_step*(j - 1))) && (val_mu_0 < x_min + bin_step*j))
					{
						h_SiPM_mu_0->Fill(val_mu_0);
						PE_in_SiPM_matrix_tmp[j - 1]++;
						Npe_in_event_counter++;
					}

				}
			}
		}

		double x_cog = PE_in_SiPM_matrix_tmp[0] * -20 + PE_in_SiPM_matrix_tmp[2] * -10 +
			PE_in_SiPM_matrix_tmp[4] * 0 + PE_in_SiPM_matrix_tmp[6] * 10 + PE_in_SiPM_matrix_tmp[8] * 20;

		x_cog /= (PE_in_SiPM_matrix_tmp[0] + PE_in_SiPM_matrix_tmp[2] + PE_in_SiPM_matrix_tmp[4]
			+ PE_in_SiPM_matrix_tmp[6] + PE_in_SiPM_matrix_tmp[8]);

		h_cog_raw_mu0->Fill(x_cog);
	}




	TCanvas *c1 = new TCanvas("c1", "c1");
	c1->Divide(2, 2, 0.01, 0.01);

	c1->cd(1);
	h_SiPM_mu_0->Draw();
	h_SiPM_mu_0->Scale(0.1/h_SiPM_mu_0->Integral());
	h_original_mu_0->Draw("same");
	h_original_mu_0->Scale(1/h_original_mu_0->Integral());

	c1->cd(2);
	h_SiPM_mu_non0->Draw();
	h_SiPM_mu_non0->Scale(0.1 / h_SiPM_mu_non0->Integral());
	h_original_mu_non0->Draw("same");
	h_original_mu_non0->Scale(1 / h_original_mu_non0->Integral());

	c1->cd(3);
	gStyle->SetOptFit(1112);
	h_cog_raw_mu0->SetStats(1);
	h_cog_raw_mu0->Draw();
	h_cog_raw_mu0->Fit("gaus");

	c1->cd(4);
	TGraph *gr = new TGraph(x_cog_raw.size(), &x_cog_raw[0], &x_real[0]);
	gr->Draw("AP");
	gr->GetXaxis()->SetTitle("x_sim");
	gr->GetYaxis()->SetTitle("x_real");
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(1);

	cout << endl;
	cout << "all is ok" << endl;
	theApp.Run();
	system("pause");
	return 0;
}



//
//int cog()
//{
//	TRandom2 rnd;
//	TRandom2 rnd_u(100);
//
//	TH1F *h_original_mu0 = new TH1F("h_original_mu0", "h_original_mu0", 100, -25, 25);
//	TH1F *h_reco = new TH1F("h_reco", "h_reco", 100, -25, 25);
//	TH1F *h_original_mu_non0 = new TH1F("h_original_mu_non0", "h_original_mu_non0", 100, -25, 25);
//
//	double x_min = -23;
//	double x_max = 23;
//	int n_bins = 9;
//	double bin_step = (x_max - x_min) / n_bins;
//	TH1F *h_SiPM_mu0 = new TH1F("h_SiPM_mu0", "h_SiPM_mu0", n_bins, x_min, x_max);
//
//	double smearing = 10.0 / 2;
//
//	int N_events = 1;
//	vector< vector<int> > PE_in_SiPM_matrix;
//
//
//	for (int i = 0; i < N_events; i++)
//	{
//		vector<int> PE_in_SiPM_matrix_tmp;
//		//PE_in_SiPM_matrix.resize(n_bins);
//
//		int Npe_in_event = 10;
//		int Npe_in_event_counter = 0;
//		//while (Npe_in_event_counter < Npe_in_event);
//		{
//			cout << "start: Npe_in_event_counter" << Npe_in_event_counter << endl;
//			double val = rnd.Gaus(0, 6.5);
//			h_original_mu0->Fill(val);
//			//h_SiPM->Fill(val);
//			for (int j = 1; j <= n_bins; j++)
//			{
//				if (j % 2 == 1)
//				{
//					if ((val > (x_min + bin_step*(j - 1))) && (val < x_min + bin_step*j))
//					{
//						h_SiPM_mu0->Fill(val);
//						PE_in_SiPM_matrix_tmp.push_back(val);
//						Npe_in_event_counter++;
//					}
//					else
//					{
//						PE_in_SiPM_matrix_tmp.push_back(0);
//					}
//
//				}
//			}
//			cout << "end: Npe_in_event_counter" << Npe_in_event_counter << endl;
//		}
//
//		h_original_mu_non0->Fill(rnd.Gaus(rnd_u.Uniform(-smearing, smearing), 6.5));
//		//PE_in_SiPM_matrix.push_back(PE_in_SiPM_matrix_tmp);
//
//		for (int k = 0; k < PE_in_SiPM_matrix_tmp.size(); k++)
//		{
//			cout << k << "\t" << PE_in_SiPM_matrix_tmp[k] << endl;
//		}
//
//
//	}
//
//	TCanvas *c1 = new TCanvas("c1", "c1");
//	c1->Divide(2, 2, 0.01, 0.01);
//
//	c1->cd(1);
//	h_original_mu0->Draw();
//
//	c1->cd(2);
//	h_SiPM_mu0->Draw();
//
//	c1->cd(3);
//	h_original_mu_non0->Draw();
//
//	for (int j = 1; j <= n_bins; j++)
//	{
//		cout << j << "\t" << h_SiPM_mu0->GetBinContent(j) << endl;
//	}
//
//	cout << "cog; rnd = " << rnd.Gaus(0, 1) << "; n_bins = " << n_bins << endl;
//
//	for (int j = 0; j < n_bins; j++)
//	{
//		//cout << j << "\t" << PE_in_SiPM_matrix[0][j] << endl; 	
//	}
//
//	return 0;
//}
