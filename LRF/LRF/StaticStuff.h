#pragma once

#include <vector>

class StaticStuff
{
public:
	StaticStuff();
	~StaticStuff();

public:
	std::vector< double > X_SiPM;
	std::vector< double > Y_SiPM;
};

