#include <TROOT.h>
#include <TAttFill.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TColor.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TStyle.h>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h                                                                                                                
#include <cstdio>
#include <cmath>
#include <sstream> //to use ostringstream to convert numbers to string in c++                                       

using namespace std;

// use as:
// root -b -q 'convert_txt_toTH2.C+("<filename>.txt","<filename>.root")'

void convert_txt_toTH2(const string txtfileName = "", const string rootfileName = "") {

  Double_t dV = -1.0;
  
  ifstream inputFile(txtfileName.c_str());

  Int_t NbinsX_2Dmap = 360;
  Double_t lowerX_2Dmap = 0.5;
  Double_t upperX_2Dmap = 360.5;
  Int_t NbinsY_2Dmap = 171;
  Double_t lowerY_2Dmap = -85.5;
  Double_t upperY_2Dmap = 85.5;

  // open file to store histogram
  TFile *rootFile = new TFile((rootfileName).c_str(),"RECREATE");
  if (!rootFile || !rootFile->IsOpen()) {
    cout << "Error: file \"" << rootfileName << "\" was not opened." << endl;
    exit(EXIT_FAILURE);
  }

  TH2F *hEB = new TH2F("hEB","cystals map in EB",NbinsX_2Dmap,lowerX_2Dmap,upperX_2Dmap,NbinsY_2Dmap,lowerY_2Dmap,upperY_2Dmap);
  TH2F *hEB_noGaps = new TH2F("hEB_noGaps","cystals map in EB without gaps",NbinsX_2Dmap,lowerX_2Dmap,upperX_2Dmap,NbinsY_2Dmap,lowerY_2Dmap,upperY_2Dmap);
  
  // file format is --> a b c , that is ieta, iphi dV                 
  Int_t a,b;
  Double_t c;
  
  if (inputFile.is_open()) {

    while ((dV < 0.0) && (inputFile >> a >> b >> c )) {

      hEB->Fill((Double_t)b,(Double_t)a,(Double_t)c);  // b is iphi
      Int_t abs_a = fabs(a);
      if (abs_a != 1 && abs_a != 25 && abs_a != 26 && abs_a != 45 && abs_a != 46 && abs_a != 65 && abs_a != 66 && abs_a != 85) {
	if ( (b%20 != 0) && (b%20 != 1) ) {  // gaps are in 1,20,21,40,41, ecc...
	  hEB_noGaps->Fill((Double_t)b,(Double_t)a,(Double_t)c); 
	}

      }

    }

  } else {
    std::cout << "Error: could not open file " << txtfileName << std::endl;
    exit(EXIT_FAILURE);
  }

  inputFile.close();

  rootFile->Write();
  rootFile->Close();
  delete rootFile;

}
