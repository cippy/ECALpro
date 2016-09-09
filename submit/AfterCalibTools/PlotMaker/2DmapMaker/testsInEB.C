#define testsInEB_cxx

#include "testsInEB.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TPaletteAxis.h>

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

#ifdef testsInEB_cxx

testsInEB::testsInEB(TTree *tree) : calibAnaEcalEB(tree) {

  ////////////////////////////                                                                                                                                         
  //initializing data members                                                                                                                                          
  ///////////////////////////                                                                                                                                          
  
  //////////////////////////                                                                                                                                           
  // protected data members                                                                                                                                            

  ////////////////////////////////                                                                                                                                     
  // public data members                    

  Init(tree);

}

//===============================================                                 

void testsInEB::setHistograms() {

  // hSignal = new TH2D("hSignal",Form("Signal in %s",EBorEE.c_str()),NbinsX_2Dmap,lowerX_2Dmap,upperX_2Dmap,NbinsY_2Dmap,lowerY_2Dmap,upperY_2Dmap);
  // hBackground = new TH2D("hBackground",Form("Background in %s",EBorEE.c_str()),NbinsX_2Dmap,lowerX_2Dmap,upperX_2Dmap,NbinsY_2Dmap,lowerY_2Dmap,upperY_2Dmap);
  // SoverB = new TH2D("SoverB",Form("S/B in %s",EBorEE.c_str()),NbinsX_2Dmap,lowerX_2Dmap,upperX_2Dmap,NbinsY_2Dmap,lowerY_2Dmap,upperY_2Dmap);     
  // SoverSqrtSplusB = new TH2D("SoverSqrtSplusB",Form("S/sqrt(S+B) in %s",EBorEE.c_str()),NbinsX_2Dmap,lowerX_2Dmap,upperX_2Dmap,NbinsY_2Dmap,lowerY_2Dmap,upperY_2Dmap);
  // SigmaMeanOverMean = new TH2D("SigmaMeanOverMean",Form("sigma(fit_mean)/fit_mean in %s",EBorEE.c_str()),NbinsX_2Dmap,lowerX_2Dmap,upperX_2Dmap,NbinsY_2Dmap,lowerY_2Dmap,upperY_2Dmap); //mean is the fit mean (should be the pi0 peak mass)  
  // mean = new TH2D("mean",Form("fit_mean in %s",EBorEE.c_str()),NbinsX_2Dmap,lowerX_2Dmap,upperX_2Dmap,NbinsY_2Dmap,lowerY_2Dmap,upperY_2Dmap);
  sigma = new TH2D("sigma",Form("fit_sigma in %s",EBorEE.c_str()),NbinsX_2Dmap,lowerX_2Dmap,upperX_2Dmap,NbinsY_2Dmap,lowerY_2Dmap,upperY_2Dmap);

  // plots with resolution in crystal (taken as sigma)
  resoInCrystalWithLowDV = new TH2D("resoInCrystalWithLowDV",Form("fit_sigma in %s crystals with low DV",EBorEE.c_str()),
				    NbinsX_2Dmap,lowerX_2Dmap,upperX_2Dmap,NbinsY_2Dmap,lowerY_2Dmap,upperY_2Dmap);
  resoInCrystalWithMediumDV = new TH2D("resoInCrystalWithMediumDV",Form("fit_sigma in %s crystals with medium DV",EBorEE.c_str()),
				       NbinsX_2Dmap,lowerX_2Dmap,upperX_2Dmap,NbinsY_2Dmap,lowerY_2Dmap,upperY_2Dmap);
  resoInCrystalWithHighDV = new TH2D("resoInCrystalWithHighDV",Form("fit_sigma in %s crystals with high dV",EBorEE.c_str()),
				     NbinsX_2Dmap,lowerX_2Dmap,upperX_2Dmap,NbinsY_2Dmap,lowerY_2Dmap,upperY_2Dmap);

  // th2dVector.push_back(hSignal);
  // th2dVector.push_back(hBackground);
  // th2dVector.push_back(SoverB);
  // th2dVector.push_back(SoverSqrtSplusB);
  // th2dVector.push_back(SigmaMeanOverMean);
  // th2dVector.push_back(mean);
  th2dVector.push_back(sigma);
  th2dVector.push_back(resoInCrystalWithLowDV);
  th2dVector.push_back(resoInCrystalWithMediumDV);
  th2dVector.push_back(resoInCrystalWithHighDV);

  th2dMinZaxisVector.push_back(0.005);
  th2dMinZaxisVector.push_back(0.005);  
  th2dMinZaxisVector.push_back(0.005);
  th2dMinZaxisVector.push_back(0.005);


}  

//=============================================== 

void testsInEB::set2DmapMaxZaxisVector() {

  // th2dMaxZaxisVector.push_back(hSignal->GetBinContent(hSignal->GetMaximumBin()));
  // th2dMaxZaxisVector.push_back(hBackground->GetBinContent(hBackground->GetMaximumBin()));
  // th2dMaxZaxisVector.push_back(10e9);
  // th2dMaxZaxisVector.push_back(10e9); // when this value is very large (bigger than the default) use the default to plot axis                  
  // th2dMaxZaxisVector.push_back(0.0125);//0.02                                                                                   
  // th2dMaxZaxisVector.push_back(0.140);
  th2dMaxZaxisVector.push_back(0.015);
  th2dMaxZaxisVector.push_back(0.015);
  th2dMaxZaxisVector.push_back(0.015);
  th2dMaxZaxisVector.push_back(0.015);

}

//===============================================                                                                                                                      

void testsInEB::setVerticalRangeInHisto() {
  
  for (UInt_t i = 0; i < th2dVector.size(); i++) {
  
    if (th2dMinZaxisVector[i] < 0.0) th2dMinZaxisVector[i] = 0.0;
    th2dVector[i]->SetMinimum(th2dMinZaxisVector[i]);

    // if the maximum choosen by user is bigger than default, don't do anything, otherwise set the user value as the maximum                                           
    if (th2dMaxZaxisVector[i] < th2dVector[i]->GetBinContent(th2dVector[i]->GetMaximumBin())) {
  
      th2dVector[i]->SetMaximum(th2dMaxZaxisVector[i]);

    }

  }


}

//=========================================================

void testsInEB::draw2Dmap(TH2D* hist2d) {

  string canvasName(hist2d->GetName());
  canvasName = "c_" + canvasName;
  TCanvas *c = new TCanvas(canvasName.c_str(),canvasName.c_str());
  string name = wwwPath + hist2d->GetName() + "_" + EBorEE;  // name  (with path) of file to save canvas: EBorEE can be "EB" or "EEp" or "EEm" 

  hist2d->Draw("COLZ");
  if (EBorEE == "EB") {
    hist2d->GetXaxis()->SetTitle("i #phi");
    hist2d->GetYaxis()->SetTitle("i #eta");
  } else {
    hist2d->GetXaxis()->SetTitle("iX");
    hist2d->GetYaxis()->SetTitle("iY");
  } 
  hist2d->GetXaxis()->SetTitleSize(0.06);
  hist2d->GetXaxis()->SetTitleOffset(0.7);
  hist2d->GetYaxis()->SetTitleSize(0.06);
  hist2d->GetYaxis()->SetTitleOffset(0.8);
  hist2d->SetStats(0);
  hist2d->Draw("COLZ");
  // after drawing, fix the palette                                                                                                                                  
  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*)hist2d->GetListOfFunctions()->FindObject("palette");
  if (!palette || palette == NULL) {
    cout << "Error in function draw2Dmap(): palette not found. ABORT" << endl;
    exit(EXIT_FAILURE);
  }
  // the following lines move the palette. Choose the values you need for the position.                                                                              
  palette->SetX1NDC(0.91);
  palette->SetX2NDC(0.94);
  //palette->SetY1NDC(0.2);                                                                                                                                          
  //palette->SetY2NDC(0.8);                                                                                                                                          
  gPad->Modified();
  gPad->Update();
  // end of palette fixes                                                                                                                                             
  c->SaveAs((name + ".pdf").c_str());
  c->SaveAs((name + ".png").c_str());

 
}

//===============================================                                                                                                                      

Double_t testsInEB::getDVinCrystal(const Int_t& i_phi, const Int_t& i_eta) {

  Double_t dV = -1.0;

  string fileName = "crystal_list.txt";
  ifstream inputFile(fileName.c_str());
  
  // file format is --> a b c , that is ieta, iphi dV                                              
  Int_t a,b;
  Double_t c;

  if (inputFile.is_open()) {

    while ((dV < 0.0) && (inputFile >> a >> b >> c )) {
 
      if (a == i_eta && b == i_phi) dV = c;

    }

  } else {
    std::cout << "Error in testsInEB::getDVinCrystal(const Int_t&, const Int_t&): could not open file " << fileName << std::endl;
    exit(EXIT_FAILURE);
  }

  inputFile.close();

  return dV;

}

//===============================================                                                                     

void testsInEB::Loop()
{  

  if (fChain == 0) return;

  setHistograms();

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;                                                                                                                              
      
    if (jentry % 100000 == 0) cout << jentry << endl;

    if ((abs(Backgr) > 0.00001) && (abs(Signal) > 0.00001)) { // avoid empty crystals due to masked towers or whatever                                                 

      normalizedS = Signal * fit_Snorm;
      normalizedB = Backgr * fit_Bnorm;

      // to avoid that in 2D maps points below lower threshold in z axis are drawn white (as if they are empty), fill with the maximum between threshold and value     
      // hSignal->Fill((Double_t)iphi,(Double_t)ieta,max(th2dMinZaxisVector[0],(Double_t)normalizedS));
      // hBackground->Fill((Double_t)iphi,(Double_t)ieta,max(th2dMinZaxisVector[1],(Double_t)normalizedB));
      // SoverB->Fill((Double_t)iphi,(Double_t)ieta,max(th2dMinZaxisVector[2],(Double_t)normalizedS/normalizedB));
      // SoverSqrtSplusB->Fill((Double_t)iphi,(Double_t)ieta, max(th2dMinZaxisVector[3],(Double_t)normalizedS/sqrt(normalizedS + normalizedB)));
      // SigmaMeanOverMean->Fill((Double_t)iphi,(Double_t)ieta, max(th2dMinZaxisVector[4],(Double_t)fit_mean_err/fit_mean));
      // mean->Fill((Double_t)iphi,(Double_t)ieta, max(th2dMinZaxisVector[5],(Double_t)fit_mean));
      sigma->Fill((Double_t)iphi,(Double_t)ieta, max(th2dMinZaxisVector[0],(Double_t)fit_sigma));

      Double_t deltaV = getDVinCrystal(iphi, ieta);
      if ( deltaV < 1.0 ) resoInCrystalWithLowDV->Fill((Double_t)iphi,(Double_t)ieta, max(th2dMinZaxisVector[1],(Double_t)fit_sigma));
      else if ( deltaV < 4.0 ) resoInCrystalWithMediumDV->Fill((Double_t)iphi,(Double_t)ieta, max(th2dMinZaxisVector[2],(Double_t)fit_sigma));
      else if ( deltaV >= 4.0 ) resoInCrystalWithHighDV->Fill((Double_t)iphi,(Double_t)ieta, max(th2dMinZaxisVector[3],(Double_t)fit_sigma));

    }

  }

  // set preference for max value in the vertical scale
  set2DmapMaxZaxisVector();  
  // now set the vertical axis maximum value based on user input (will use the least between the default and the user input). 
  setVerticalRangeInHisto();

  for ( UInt_t i = 0; i < th2dVector.size(); i++ ) {

    draw2Dmap(th2dVector[i]);

  }

}

#endif
