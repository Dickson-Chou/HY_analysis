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
const double L2016 = 35.9*1000; // unit pb
// xs unit pb & number 


 const double xs_500_700 = 29980;
 const double xs_700_1000 = 6334;
 const double xs_1000_1500 = 1088;
 const double xs_1500_2000 = 99.11;
 const double xs_2000_inf = 20.23;

 

using namespace std;

void draw(std::string outputFileName="result.root"){ 
/*
void draw(std::string filename="*.root", std::string outputFileName="result.root"){
           
  std::vector<std::string> inputFiles;

  // check first if this is a root file
  if(filename.find(".root")!=string::npos)
    {
      cout << "This is a single input root file" << endl;
      inputFiles.push_back(filename);
    }
  else // assume this is a text file
    {
      cout << "This is a text file " << endl;
      ifstream fin;
      fin.open(filename.data());
      string temp;
      fin >> temp;
      while(!fin.eof())
	{
	  inputFiles.push_back(temp);
	  fin >> temp;
	}
      cout << "There are " << inputFiles.size() << " files" << endl;
      for(unsigned int ifile=0; ifile < inputFiles.size(); ifile++)
	cout << "Input file " << ifile << " is " << inputFiles[ifile] << endl;
    }
  cout << "Output file name is " << outputFileName << endl;

  setNCUStyle(true);
  //get TTree from file ...  
  TreeReader data(inputFiles,"variable");
  TTree* thisTree = data.GetTree();
*/
 
//  Float_t*  FatJet_ParticleNetMD_ProbXbb = data.GetPtrFloat("FatJet_ParticleNetMD_ProbXbb");
//  Float_t*  I_FatJet_ParticleNetMD_ProbXbb = data.GetPtrFloat("I_FatJet_ParticleNetMD_ProbXbb");


  //get TH1F data from file ...
  TFile* signal = TFile::Open("signal.root");
  TFile* back_500_700 = TFile::Open("back_500_700.root");
  TFile* back_700_1000 = TFile::Open("back_700_1000.root");
  TFile* back_1000_1500 = TFile::Open("back_1000_1500.root");
  TFile* back_1500_2000 = TFile::Open("back_1500_2000.root");
  TFile* back_2000_inf = TFile::Open("back_2000_inf.root");
  TFile* ttbar = TFile::Open("ttbar.root");

  TTree *tSig = (TTree *)signal->Get("variable");
  TTree *tback_500_700 = (TTree *)back_500_700->Get("variable");
  TTree *tback_1500_2000 = (TTree *)back_1500_2000->Get("variable");
  TTree *tback_2000_inf = (TTree *)back_2000_inf->Get("variable");
  

  TH1F* Ydistri_initial = (TH1F*) signal->Get("Ydistri_initial");
  TH1F* Ydistri_P1 = (TH1F*) signal->Get("Ydistri_P1");
  TH1F* Ydistri_P2 = (TH1F*) signal->Get("Ydistri_P2");
  TH1F* Ydistri_P3 = (TH1F*) signal->Get("Ydistri_P3");
  TH1F* Ydistri_P4 = (TH1F*) signal->Get("Ydistri_P4");
  TH1F* Ydistri_P5 = (TH1F*) signal->Get("Ydistri_P5");
  TH1F* Ydistri_P6 = (TH1F*) signal->Get("Ydistri_P6");
  TH1F* Ydistri_P7 = (TH1F*) signal->Get("Ydistri_P7");
  TH1F* Ydistri_P8 = (TH1F*) signal->Get("Ydistri_P8");

  TH1F* Unequal_Ydistri_initial = (TH1F*) signal->Get("Unequal_Ydistri_initial");
  TH1F* Unequal_Ydistri_P1 = (TH1F*) signal->Get("Unequal_Ydistri_P1");
  TH1F* Unequal_Ydistri_P2 = (TH1F*) signal->Get("Unequal_Ydistri_P2");
  TH1F* Unequal_Ydistri_P3 = (TH1F*) signal->Get("Unequal_Ydistri_P3");
  TH1F* Unequal_Ydistri_P4 = (TH1F*) signal->Get("Unequal_Ydistri_P4");
  TH1F* Unequal_Ydistri_P5 = (TH1F*) signal->Get("Unequal_Ydistri_P5");
  TH1F* Unequal_Ydistri_P6 = (TH1F*) signal->Get("Unequal_Ydistri_P6");
  TH1F* Unequal_Ydistri_P7 = (TH1F*) signal->Get("Unequal_Ydistri_P7");
  TH1F* Unequal_Ydistri_P8 = (TH1F*) signal->Get("Unequal_Ydistri_P8");

  TH1F* HBB = (TH1F*) signal->Get("HBB_D");
  TH1F* hbb = (TH1F*) signal->Get("hbb");

  TH1F* HBB_5 = (TH1F*) back_500_700->Get("HBB_D");
  TH1F* HBB_7 = (TH1F*) back_700_1000->Get("HBB_D");
  TH1F* HBB_10 = (TH1F*) back_1000_1500->Get("HBB_D");
  TH1F* HBB_15 = (TH1F*) back_1500_2000->Get("HBB_D");
  TH1F* HBB_20 = (TH1F*) back_2000_inf->Get("HBB_D");
  TH1F* HBB_tt = (TH1F*) ttbar->Get("HBB_D");
 

  TH1F* HDeep = (TH1F*) signal->Get("HDeep");
  TH1F* HDeep_5 = (TH1F*) back_500_700->Get("HDeep");
  TH1F* HDeep_7 = (TH1F*) back_700_1000->Get("HDeep");
  TH1F* HDeep_10 = (TH1F*) back_1000_1500->Get("HDeep");
  TH1F* HDeep_15 = (TH1F*) back_1500_2000->Get("HDeep");
  TH1F* HDeep_20 = (TH1F*) back_2000_inf->Get("HDeep");

  TH1F* Hdeep = (TH1F*) signal->Get("Hdeep");

  TH1F* HT_5 = (TH1F*) back_500_700->Get("HT");
  TH1F* HT_7 = (TH1F*) back_700_1000->Get("HT");
  TH1F* HT_10 = (TH1F*) back_1000_1500->Get("HT");
  TH1F* HT_15 = (TH1F*) back_1500_2000->Get("HT");
  TH1F* HT_20 = (TH1F*) back_2000_inf->Get("HT");


 // TH1F* HBB = new TH1F("HBB"," ",100,0,1);


//  TCanvas *c1 = new TCanvas("c1","c1",3);
  
  TH1F* bHBB = new TH1F("bHBB","",100,0,1);
  TH1F* bHDeep = new TH1F("bHDeep","",100,0,1);
  TH1F* THT = new TH1F("THT","",100,0,2500);
  TH1F* Total = (TH1F*) signal->Get("heve");

  TH1F* heve_5 = (TH1F*) back_500_700->Get("heve");
  TH1F* heve_7 = (TH1F*) back_700_1000->Get("heve");
  TH1F* heve_10 = (TH1F*) back_1000_1500->Get("heve");
  TH1F* heve_15 = (TH1F*) back_1500_2000->Get("heve");
  TH1F* heve_20 = (TH1F*) back_2000_inf->Get("heve");

  TH1F* eff_final = new TH1F("Eff_final","",16,0,16);
  TH1F* eff1 = new TH1F("EFF1","",16,0,16);
  TH1F* eff2 = new TH1F("EFF2","",16,0,16);
  TH1F* eff3 = new TH1F("EFF3","",16,0,16);
  TH1F* eff4 = new TH1F("EFF4","",16,0,16);
  TH1F* eff5 = new TH1F("EFF5","",16,0,16);
  TH1F* eff6 = new TH1F("EFF6","",16,0,16);
  TH1F* eff7 = new TH1F("EFF7","",16,0,16);
  TH1F* eff8 = new TH1F("EFF8","",16,0,16);

  TH1F* eff_final_UN = new TH1F("Eff_final_UN","",16,0,16);
  TH1F* eff1_UN = new TH1F("EFF1_UN","",16,0,16);
  TH1F* eff2_UN = new TH1F("EFF2_UN","",16,0,16);
  TH1F* eff3_UN = new TH1F("EFF3_UN","",16,0,16);
  TH1F* eff4_UN = new TH1F("EFF4_UN","",16,0,16);
  TH1F* eff5_UN = new TH1F("EFF5_UN","",16,0,16);
  TH1F* eff6_UN = new TH1F("EFF6_UN","",16,0,16);
  TH1F* eff7_UN = new TH1F("EFF7_UN","",16,0,16);
  TH1F* eff8_UN = new TH1F("EFF8_UN","",16,0,16);


  
  		//definition
  Float_t YM[16]={90,100,125,150,200,250,300,400,500,600,700,800,900,1000,1200,1400};//vector of ymass points
  double ym[17]={70,95,120,145,195,245,295,395,495,595,695,795,895,995,1195,1395,1595};
  std::vector<float> YDistri_initial;
  std::vector<float> YDistri;
  std::vector<float> total;
  std::vector<float> Ryy(16,0);
  std::vector<float> Ryt(16,0);
  

  

  for(int i = 1 ; i < 17 ; i++ ){
//  	YDistri_initial.push_back(Ydistri_initial->GetBinContent(i));
  	YDistri.push_back(Ydistri_initial->GetBinContent(i));
  	total.push_back(Total->GetBinContent(3));
  }
 
  //th1f distri draw
 
/*
  for(int i = 0 ; i < 16 ; i++){

  	Ryy[i] = YDistri[i]/YDistri_initial[i];
  	Ryt[i] = YDistri[i]/total[i];
  	effyy->Fill(i,Ryy[i]);
  	effyy->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
  	effyt->Fill(i,Ryt[i]);
  	effyt->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
  }
 */ 
//  cout << "initial = " << YDistri_initial[2] << endl;
//  cout << "Ydistri = " << YDistri[2] << endl;
//  cout << "npass[1] = " << total[2] << endl;

//divide::
  auto c1 = new TCanvas();
/*
  eff1_UN->Divide(Unequal_Ydistri_P1,Unequal_Ydistri_initial,1,1,"B");
  eff2_UN->Divide(Unequal_Ydistri_P2,Unequal_Ydistri_initial,1,1,"B");
  eff3_UN->Divide(Unequal_Ydistri_P3,Unequal_Ydistri_initial,1,1,"B");
  eff4_UN->Divide(Unequal_Ydistri_P4,Unequal_Ydistri_initial,1,1,"B");
  eff5_UN->Divide(Unequal_Ydistri_P5,Unequal_Ydistri_initial,1,1,"B");
  eff6_UN->Divide(Unequal_Ydistri_P6,Unequal_Ydistri_initial,1,1,"B");
  eff7_UN->Divide(Unequal_Ydistri_P7,Unequal_Ydistri_initial,1,1,"B");
  eff8_UN->Divide(Unequal_Ydistri_P8,Unequal_Ydistri_initial,1,1,"B");

 */

  eff1->Divide(Ydistri_P1,Ydistri_initial,1,1,"B");
  eff2->Divide(Ydistri_P2,Ydistri_initial,1,1,"B");
  eff3->Divide(Ydistri_P3,Ydistri_initial,1,1,"B");
  eff4->Divide(Ydistri_P4,Ydistri_initial,1,1,"B");
  eff5->Divide(Ydistri_P5,Ydistri_initial,1,1,"B");
  eff6->Divide(Ydistri_P6,Ydistri_initial,1,1,"B");
  eff7->Divide(Ydistri_P7,Ydistri_initial,1,1,"B");
//  eff8->Divide(Ydistri_P8,Ydistri_initial,1,1,"B");



  for(int i = 0 ; i < 16 ; i++){
  	
  	eff1->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
  	eff2->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
 	eff3->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
 	eff4->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i])); 
  	eff5->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
  	eff6->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
  	eff7->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
//  	eff8->GetXaxis()->SetBinLabel(i+1,Form("%4.2f",YM[i]));
  	
//  	eff3_UN->Fill(YM[i]);
//  	eff3_UN->SetBins(16,ym);
  }

  gStyle->SetOptTitle(kFALSE);
  gStyle->SetOptStat("");//cancel the label
  	
 // eff->GetMaximum()
 // subleading cut eff
  eff1->SetMinimum(0);
  eff1->SetLineColor(kRed);
  eff1->SetMarkerColor(kRed);
  eff1->SetMarkerStyle(kFullCircle);
  eff2->SetLineColor(kBlue+2);
  eff2->SetMarkerColor(kBlue+2);
  eff2->SetMarkerStyle(kFullSquare);
  eff3->SetLineColor(kGreen+3);
  eff3->SetMarkerColor(kGreen+3);
  eff3->SetMarkerStyle(kFullTriangleUp);
  eff4->SetLineColor(kRed+3);
  eff4->SetMarkerColor(kRed+3);
  eff4->SetMarkerStyle(kFullTriangleDown);

  eff1->Draw();
  eff2->Draw("SAME");
  eff3->Draw("SAME");
  eff4->Draw("SAME");

  auto legend = new TLegend(0.6,0.25,0.85,0.5);//(x1,y1,x2,y2)
//  legend->SetHeader("cut2 to cut 5","C"); // option "C" allows to center the header
  legend->SetNColumns(2);
  legend->AddEntry(eff1,"subleading","lep");
  legend->AddEntry(eff2,"+leading","lep");
  legend->AddEntry(eff3,"nGoodpair","lep");
  legend->AddEntry(eff4,"masspass","lep");
  legend->Draw();

  c1->Modified();
  c1->Print("Divide_eff_sub.jpg");

  auto c2 = new TCanvas();

  eff5->SetMinimum(0);
  eff5->SetLineColor(kRed);
  eff5->SetMarkerColor(kRed);
  eff5->SetMarkerStyle(kFullCircle);
  eff6->SetLineColor(kBlue+2);
  eff6->SetMarkerColor(kBlue+2);
  eff6->SetMarkerStyle(kFullSquare);
  eff7->SetLineColor(kGreen+3);
  eff7->SetMarkerColor(kGreen+3);
  eff7->SetMarkerStyle(kFullTriangleUp);
//  eff8->SetLineColor(kRed+3);
//  eff8->SetMarkerColor(kRed+3);
//  eff8->SetMarkerStyle(kFullTriangleDown);


  eff6->Draw();
  eff7->Draw("SAME");
  eff5->Draw("SAME");
//  eff8->Draw("SAME");

  auto legend_1 = new TLegend(0.65,0.65,0.9,0.9);//(x1,y1,x2,y2)
//  legend_1->SetHeader("cut6 to cut 10","C"); // option "C" allows to center the header
  legend_1->SetNColumns(2);
  legend_1->AddEntry(eff5,"deepTag","lep");
  legend_1->AddEntry(eff6,"particle net","lep");
  legend_1->AddEntry(eff7,"deep+particle","lep");
//  legend_1->AddEntry(eff8,"cut9","lep");

  legend_1->Draw();

  c2->Modified();
  c2->Print("Divide_eff_lead.jpg");

  auto c3 = new TCanvas();

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
  	cout << "weight[" << i << "]" << weight[i] << endl;
  }
  
  bHBB->Add(HBB_5,weight[0]);
  bHBB->Add(HBB_7,weight[1]);
  bHBB->Add(HBB_10,weight[2]);
  bHBB->Add(HBB_15,weight[3]);
  bHBB->Add(HBB_20,weight[4]);
  //h->Integral(begin,end,option);
 // bHBB->Integral(1,100)
  bHBB->SetLineColor(kRed);
  bHBB->SetLineWidth(2);
  HBB->SetLineColor(kBlue);
  HBB->SetLineWidth(2);
  HBB_tt->SetLineColor(kGreen-7);
  HBB_tt->SetLineWidth(2);


  HBB->DrawNormalized("hist");
  bHBB->DrawNormalized("histsame");
  HBB_tt->DrawNormalized("histsame");

  auto legend_particle = new TLegend(0.6,0.25,0.85,0.5);//(x1,y1,x2,y2)
//  legend->SetHeader("cut2 to cut 5","C"); // option "C" allows to center the header
  legend_particle->SetNColumns(2);
  legend_particle->AddEntry(bHBB,"QCD_PNet","f");
  legend_particle->AddEntry(HBB,"signal_PNet","f");
  legend_particle->AddEntry(HBB_tt,"ttbar","f");
  legend_particle->Draw();

  auto c4 = new TCanvas();

 /*
  bHDeep->Add(HDeep_5,weight[0]);
  bHDeep->Add(HDeep_7,weight[1]);
  bHDeep->Add(HDeep_10,weight[2]);
  bHDeep->Add(HDeep_15,weight[3]);
  bHDeep->Add(HDeep_20,weight[4]);
  bHDeep->SetLineColor(kRed);
  bHDeep->SetFillColor(kRed);

  bHDeep->Draw("hist");
  HDeep->Draw("hist SAME");

  auto legend_deep = new TLegend(0.6,0.25,0.85,0.5);//(x1,y1,x2,y2)
//  legend->SetHeader("cut2 to cut 5","C"); // option "C" allows to center the header
  legend_deep->SetNColumns(2);
  legend_deep->AddEntry(bHDeep,"QCD_Deep","f");
  legend_deep->AddEntry(HDeep,"signal_Deep","f");
  legend_deep->Draw();
  */
 
  THT->Add(HT_5,weight[0]);
  THT->Add(HT_7,weight[1]);
  THT->Add(HT_10,weight[2]);
  THT->Add(HT_15,weight[3]);
  THT->Add(HT_20,weight[4]);
  c4->SetLogy();
  THT->Draw("hist");
  /*
  c2->Divide(2,1);
//  auto h1 = (TH1F*)Ydistri_initial->Clone();
//  h1->Reset();
  c2->cd(1);
  Ydistri_P3->SetLineColor(kBlue+3);
  Ydistri_P4->SetLineColor(kOrange+9);

  
  Ydistri_P3->SetFillColor(kBlue+3);
  Ydistri_P4->SetFillColor(kOrange+9);
//  Ydistri_P3->SetTitle("Y Distribution of N-1");
  Ydistri_P3->GetXaxis()->SetTitle("Y Mass points");
  Ydistri_P3->GetYaxis()->SetTitle("Number of Y");
  Ydistri_P3->Draw("hist");
  Ydistri_P4->Draw("hist SAME");

  auto legend_2 = new TLegend(0.7,0.7,0.9,0.9);//(x1,y1,x2,y2)
//  legend_2->SetHeader("Distri of ","C"); // option "C" allows to center the header
//  legend_2->SetNColumns(2);
  legend_2->AddEntry(Ydistri_P3," no particle net","f");
  legend_2->AddEntry(Ydistri_P4,"particle net","f");
  legend_2->Draw();

  c2->cd(2);
  eff3->SetLineColor(kBlue+3);
  eff3->SetMarkerColor(kBlue+3);
  eff3->SetMarkerStyle(kFullCircle);
  eff4->SetLineColor(kOrange+9);
  eff4->SetMarkerColor(kOrange+9);
  eff4->SetMarkerStyle(kFullCircle);
//  eff3->GetXaxis()->SetTitle("Y Mass points");
//  eff3->GetYaxis()->SetTitle("efficiency");
//  eff3->GetXaxis()->SetTitleSize(1);
//  eff3->GetYaxis()->SetTitleSize(1);

  eff3->Draw();
  eff4->Draw("SAME");

  auto legend_3 = new TLegend(0.7,0.7,0.9,0.9);//(x1,y1,x2,y2)
//  legend_3->SetHeader("cut6 to cut 10","C"); // option "C" allows to center the header
//  legend_1->SetNColumns(2);
  legend_3->AddEntry(eff3,"no particle net","lep");
  legend_3->AddEntry(eff4,"particle net","lep");
  legend_3->Draw();

  c2->Modified();
  c2->Print("N-1_Net.pdf");

  auto c3 = new TCanvas();
  Unequal_Ydistri_P3->SetLineColor(kBlue+3);
  Unequal_Ydistri_P4->SetLineColor(kOrange+9);

  Unequal_Ydistri_P3->SetFillColor(kBlue+3);
  Unequal_Ydistri_P4->SetFillColor(kOrange+9);

  Unequal_Ydistri_P3->Draw("hist");
  Unequal_Ydistri_P4->Draw("hist SAME");

  auto legend_4 = new TLegend(0.7,0.7,0.9,0.9);//(x1,y1,x2,y2)
//  legend_2->SetHeader("Distri of ","C"); // option "C" allows to center the header
//  legend_2->SetNColumns(2);
  legend_4->AddEntry(Unequal_Ydistri_P3," no particle net","f");
  legend_4->AddEntry(Unequal_Ydistri_P4,"particle net","f");
  legend_4->Draw();

  c3->Modified();
  c3->Print("N-1_Net_un.pdf");
  
*/
/*
  HBB->SetLineColor(kBlue+3);
  HBB->SetFillColor(kBlue+3);
  hbb->SetLineColor(kOrange);
  hbb->SetFillColor(kOrange);
  HBB->Draw();
  hbb->Draw("SAME");

  auto l = new TLegend(0.15,0.7,0.35,0.9);
  l->AddEntry(HBB,"No particle net","f");
  l->AddEntry(hbb,"paritcle net","f");
  l->Draw();
  c1->Draw();

  auto c2 = new TCanvas();

  HDeep->SetLineColor(kGreen+3);
  HDeep->SetFillColor(kGreen+3);
  Hdeep->SetLineColor(kBlue-2);
  Hdeep->SetFillColor(kBlue-2);
  HDeep->Draw();
  Hdeep->Draw("SAME");

  auto l2 = new TLegend(0.15,0.7,0.35,0.9);
  l2->AddEntry(HDeep,"No DeepAK8","f");
  l2->AddEntry(Hdeep,"DeepAK8","f");
  l2->Draw();
  c2->Draw();

 */
//TGraphAsymmErrors::

//  c1->Update();

//  auto gr = new TGraphAsymmErrors(Ydistri , Ydistri_initial , "cl = 0.683 b(1,1) mode");

//  Ydistri_P1->Draw("PLC PMC");
//  Ydistri_P2->Draw("SAME PLC PMC");
//  c1->Modified();
//  c1->Print("TGraphAsymmErrors_eff.pdf");
  
  //store 
  TFile* outFile = new TFile(outputFileName.data(),"recreate");
//  Ydistri_initial->Write();
//  eff1->Write();
//  eff2->Write();
//  eff3->Write();
//  eff3_UN->Draw();
//  gr->Write();
  outFile->Close();
}






