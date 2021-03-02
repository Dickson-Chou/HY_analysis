// draw the graph after the selection
// example code to run 2016 NSSM MC X->Y+H samples
// .L xAna_nano_nssm.C++
// xAna_nano_nssm("test.root") or xAna_nano_nssm("input.txt")
// example root file is at /afs/cern.ch/work/s/syu/public/forTiKai/nssm_nano.root

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TH1D.h>
#include <TFile.h>
#include "untuplizer.h"
#include <TLorentzVector.h>
#include "setNCUStyle.C"
#include "tdrstyle.C"
#include "CMS_lumi.C"
#include "CMS_lumi.h"
const double L2016 = 35.9*1000; // unit pb
// xs unit pb & number 

using namespace std;

//define punzi_eq
  	double punzi(double sigeff, double bg){
  		return sigeff/(1+TMath::Power(bg,0.5));
  	}


void opt(std::string outputFileName = "opt.root"){

//open file
  TFile* signal = TFile::Open("signal.root");
  TFile* back_500_700 = TFile::Open("back_500_700.root");
  TFile* back_700_1000 = TFile::Open("back_700_1000.root");
  TFile* back_1000_1500 = TFile::Open("back_1000_1500.root");
  TFile* back_1500_2000 = TFile::Open("back_1500_2000.root");
  TFile* back_2000_inf = TFile::Open("back_2000_inf.root");
  TFile* ttbar = TFile::Open("ttbar.root");

  TH1F* HBB = (TH1F*) signal->Get("signal_d_ 300");
  

  TH1F* HBB_5 = (TH1F*) back_500_700->Get("YPN");
  TH1F* HBB_7 = (TH1F*) back_700_1000->Get("YPN");
  TH1F* HBB_10 = (TH1F*) back_1000_1500->Get("YPN");
  TH1F* HBB_15 = (TH1F*) back_1500_2000->Get("YPN");
  TH1F* HBB_20 = (TH1F*) back_2000_inf->Get("YPN");
  TH1F* HBB_tt = (TH1F*) ttbar->Get("YPN");

  TH1F* heve_s = (TH1F*) signal->Get("heve");
  TH1F* heve_5 = (TH1F*) back_500_700->Get("heve");
  TH1F* heve_7 = (TH1F*) back_700_1000->Get("heve");
  TH1F* heve_10 = (TH1F*) back_1000_1500->Get("heve");
  TH1F* heve_15 = (TH1F*) back_1500_2000->Get("heve");
  TH1F* heve_20 = (TH1F*) back_2000_inf->Get("heve");

// merge back_HBB
  TH1F* bHBB = new TH1F("bHBB","",100,0,1.0);

  TH1F* punzi_Net = new TH1F("punzi_d_300","",100,0,1.0);
  

  TH1F* effs_b = new TH1F("S vs B","",100,0,1);
  	std::vector<double> xs;
 		xs.push_back(29980);
 		xs.push_back(6334);
 		xs.push_back(1088);
 		xs.push_back(99.11);
 		xs.push_back(20.23);

  	std::vector<double> number;
  		number.push_back(heve_5->GetBinContent(1));
  		number.push_back(heve_7->GetBinContent(1));
  		number.push_back(heve_10->GetBinContent(1));
  		number.push_back(heve_15->GetBinContent(1));
  		number.push_back(heve_20->GetBinContent(1));


  	std::vector<double> weight;
  	double W_B = 0;
  	for (int i = 0 ; i < 5 ; i++){
  		W_B = xs[i]*35900/number[i];
  		weight.push_back(W_B);
  		}
  	for(int i = 0 ; i < 5 ; i++){
  		cout << "weight[i]" << weight[i] << endl;
  		}
  
  	bHBB->Add(HBB_5,weight[0]);
  	bHBB->Add(HBB_7,weight[1]);
  	bHBB->Add(HBB_10,weight[2]);
  	bHBB->Add(HBB_15,weight[3]);
  	bHBB->Add(HBB_20,weight[4]);

 //punzi
  	double nSigEventTotal = heve_s->GetBinContent(2);
  	double nSigEvent = HBB->Integral();
  	double nBgEvent = bHBB->Integral();
  	int nBin = HBB->GetSize();
  	cout << nBin << endl;
  	vector <vector<float>> effs(2);//2-D vector
  	float Teffs[100], Teffb[100];
  	vector <vector<double>> effs_t(2);//with total
  	vector <vector<float>> effb(2);//bg
  	vector <vector<double>> effb_t(2);//with total
  	vector <vector<double>> punziList(2);
  	double event[2][2] = {0};

  	cout << nSigEvent << "," << nBin << endl; 
  	for(int i = 0 ; i < nBin ; i++){ 
  	  	event[0][0] += HBB->GetBinContent(i);//from 0 to end for sig
  	  	event[1][1] += HBB->GetBinContent(nBin-i-1);//from end to 0 for sig
  	  	event[0][1] += bHBB->GetBinContent(i);
  	  	event[1][0] += bHBB->GetBinContent(nBin-i-1);
  	  	effs[0].push_back(event[0][0]/nSigEvent);//(sig/other)
  	  	effs[1].push_back(event[1][1]/nSigEvent);
  	  	effs_t[0].push_back(event[0][0]/nSigEventTotal);//(sig/total)
  	  	effs_t[1].push_back(event[1][1]/nSigEventTotal);
  	  	effb[0].push_back(event[0][1]/nBgEvent);
  	  	effb[1].push_back(event[1][0]/nBgEvent);
  	  	punziList[0].push_back(punzi(event[0][0]/nSigEvent,event[0][1]));
  	  	punziList[1].push_back(punzi(event[1][1]/nSigEvent,event[1][0]));
  	  	if(effb[1][i] > 1){
  	  		effb[1][i] = 1;
  	  	}
  	  }
  	  	reverse(effs[1].begin(),effs[1].end());
  		reverse(effb[1].begin(),effb[1].end());
  		reverse(punziList[1].begin(),punziList[1].end());
  		
  	
  	for(int i = 1 ; i < nBin-1 ; i++){
  		punzi_Net->Fill(i/double(nBin-1),punziList[1][i]);
  		Teffs[i] = effs[1][i];
  		Teffb[i] = 1 - effb[1][i];
  		cout << Teffb[i] << " || " << Teffs[i] << " || " << punziList[1][i] << endl; 
  	}
  	TGraph *sigEff_vs_bkgEff = new TGraph(101, Teffs, Teffb);

  	sigEff_vs_bkgEff->SetMaximum(1.2);
  	sigEff_vs_bkgEff->Draw("ap");
  	sigEff_vs_bkgEff->SetTitle("");
  	sigEff_vs_bkgEff->GetXaxis()->SetTitle("Sig Efficiency");
  	sigEff_vs_bkgEff->GetXaxis()->SetLimits(0 , 1.0);
  	sigEff_vs_bkgEff->GetXaxis()->SetTickSize(0.03);
  	sigEff_vs_bkgEff->GetXaxis()->SetTitleOffset(1.2); 
  	sigEff_vs_bkgEff->GetXaxis()->SetLabelOffset(0.015); 
  	sigEff_vs_bkgEff->GetYaxis()->SetTitle("Bkg rejection Efficiency");
  	sigEff_vs_bkgEff->GetYaxis()->SetTitleOffset(1.3);
  	//sigEff_vs_bkgEff->GetYaxis()->SetNdivisions(505);
 	sigEff_vs_bkgEff->GetYaxis()->SetTickSize(0.03);
  	sigEff_vs_bkgEff->GetYaxis()->SetLabelOffset(0.005);
  	sigEff_vs_bkgEff->GetHistogram()->SetMaximum(1.2); // along
  	sigEff_vs_bkgEff->GetHistogram()->SetMinimum(0.);

  	sigEff_vs_bkgEff->SetLineColor(kRed);
  	sigEff_vs_bkgEff->SetLineWidth(3);
  	sigEff_vs_bkgEff->SetMarkerColor(kBlue);
  	sigEff_vs_bkgEff->SetMarkerSize(1);
  	sigEff_vs_bkgEff->SetMarkerStyle(20);
  	sigEff_vs_bkgEff->Draw("p");

  	


  	THStack *hs = new THStack("hs","histograms");
  	bHBB->SetFillColor(kBlue);
  	hs->Add(bHBB);
  	HBB_tt->SetFillColor(kGreen);
  	hs->Add(HBB_tt);
//  	HBB->SetLineColor(kRed);
//  	hs->Add(HBB);

  	auto c1 = new TCanvas();
  		punzi_Net->GetXaxis()->SetTitle("Particle Net");
  		punzi_Net->GetYaxis()->SetTitle("Punzi_significance");
  		punzi_Net->GetYaxis()->SetTitleOffset(1.4);
		punzi_Net->SetLineColor(kRed);
  		punzi_Net->SetLineWidth(2);
  		punzi_Net->Draw("hist");
/*  	
  	auto c2 = new TCanvas();
  		hs->Draw();
*/
//  	auto c3 = new TCanvas();
//  		sigEff_vs_bkgEff->Draw("l");

}

