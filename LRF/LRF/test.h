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

class test
{
public:
	test();
	~test();
private:
	//TGraph* gr;
	void FillTGraph2D();
	static StaticStuff staticStuff; // constructor runs once, single instance
	TGraph2D* gr2D;
	static int N_instanse;
};

