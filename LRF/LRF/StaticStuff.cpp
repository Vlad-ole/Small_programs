#include "StaticStuff.h"
#include <iostream>
using namespace std;

StaticStuff::StaticStuff()
{
	int matrix_size = 9;
	double pitch = 10;
	for (int i = 0; i < matrix_size; i++)
	{
		for (double j = -40; j <= 40; j += pitch)
		{
			X_SiPM.push_back(j);
		}
	}

	for (double j = 40; j >= -40; j -= pitch)
	{
		for (int i = 0; i < matrix_size; i++)
		{
			Y_SiPM.push_back(j);
		}
	}


	//test
	/*for (int i = 0; i < X_SiPM.size(); i++)
	{
		cout << X_SiPM[i] << "\t" << Y_SiPM[i] << endl;
	}*/
}


StaticStuff::~StaticStuff()
{
}
