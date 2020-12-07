#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

int main()
{
	FILE *f = fopen("D:\\git_repositories\\PhD_dark_matter\\Double_channel\\Double_channel\\in.dat", "r");
	ofstream file_out("D:\\git_repositories\\PhD_dark_matter\\Double_channel\\Double_channel\\out.dat");

	if (f == NULL)
	{
		cout << "Can't open this file: " << endl;
		system("pause");
		return 0;
	}

	double x, y;
	vector<double> xv;
	vector<double> yv;

	while (!feof(f))
	{
		fscanf(f, "%lf %lf\n", &x, &y);
		xv.push_back(x);
		yv.push_back(y);
	}

	for (int i = 0; i < xv.size() / 2; i++)
	{
		int j = i * 2;
		file_out << (xv[j] + xv[j+1]) / 2.0 << "\t" << yv[j] + yv[j+1] << endl;
	}



	system("pause");
	return 0;
}