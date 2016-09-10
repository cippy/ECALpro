#ifndef testsInEB_h
#define testsInEB_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <TTree.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TPaletteAxis.h>

#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>

#include "calibAnaEcalEB.h"

class testsInEB : public calibAnaEcalEB {
 public:

  testsInEB(TTree *tree); 
  virtual ~testsInEB() { std::cout<<"~testsInEB() called"<<std::endl; }

  ///////////////////////////////////////
  // public member functions
  virtual void setHistograms();
  virtual void draw2Dmap(TH2D*);
  virtual void set2DmapMaxZaxisVector();
  virtual void setVerticalRangeInHisto();
  virtual void Loop();
  virtual Double_t getDVinCrystal(const Int_t&, const Int_t&); // use the file given by Francesca to get DeltaV between the APDs in gain 50 in each crystal. 

  TH2D* sigma_noGaps = NULL;
  TH2D* resoInCrystalWithLowDV = NULL;
  TH2D* resoInCrystalWithMediumDV = NULL;
  TH2D* resoInCrystalWithHighDV = NULL;

  ////////////////////////////////////////////////////
  // variables used in member functions, such as Loop()


  ///////////////////////////////////////////////////
  // member functions to access protected data member

  ///////////////////////////////////////////////////
  // private or protected data members


};


#endif
