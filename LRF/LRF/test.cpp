#include "test.h"

using namespace std;

StaticStuff test::staticStuff;
int test::N_instanse = 0;

test::test()
{
	N_instanse++;
	/*vector<double> x;
	vector<double> y;
	for (int i = 0; i < 10; i++)
	{
		x.push_back(i);
		y.push_back(i);
	}
	gr = new TGraph(x.size(), &x[0], &y[0]);*/

	FillTGraph2D();
}

void test::FillTGraph2D()
{
	vector<double> n_pe;
	for (int i = 0; i < 10; i++)
	{
		n_pe.push_back(i);
	}

	//cout << "&gr2D = " << &gr2D << endl;	
	gr2D = new TGraph2D(staticStuff.X_SiPM.size(), &staticStuff.X_SiPM[0], &staticStuff.Y_SiPM[0], &n_pe[0]);
	ostringstream sst;
	sst << "gr2D_" << N_instanse;
	gr2D->SetName(sst.str().c_str());
	//cout << "&gr2D = " << &gr2D << endl;


	//h2D = gr2D->GetHistogram();
}


test::~test()
{
	delete gr2D;
}
