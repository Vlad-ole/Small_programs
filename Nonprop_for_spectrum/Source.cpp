#include <iostream>
#include <vector>
#include <fstream>

#include "Math/Polynomial.h"
#include "Math/Interpolator.h"

#include <algorithm>

using namespace std;



int main()
{
	double x, y;
	
	ifstream file_data("D:\\git_repositories\\Small_programs\\Nonprop_for_spectrum\\YAP_Ce_rel_662keV.dat");
	vector<double> Ev;
	vector<double> Nonpropv;
	while (file_data.good())
	{
		file_data >> x >> y;
		Ev.push_back(x);
		Nonpropv.push_back(y);
	}
	ROOT::Math::Interpolator inter(Ev.size(), ROOT::Math::Interpolation::kLINEAR);
	inter.SetData(Ev, Nonpropv);


	ifstream file_in("D:\\git_repositories\\Small_programs\\Nonprop_for_spectrum\\input.dat");
	vector<double> E2v;
	vector<double> Countsv;
	while (file_in.good())
	{
		file_in >> x >> y;
		E2v.push_back(x);
		Countsv.push_back(y);
	}
	

	ofstream file_out("D:\\git_repositories\\Small_programs\\Nonprop_for_spectrum\\output.dat");
	for (int i = 0; i < E2v.size(); i++)
	{
		file_out << inter.Eval(E2v[i])*E2v[i] << endl;
	}
	


	system("pause");
	return 0;
}