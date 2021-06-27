
//======================================================
//update:2021/4/20 different delta_R for particle net
//======================================================

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TH1D.h>
#include <TFile.h>
#include "untuplizer.h"
#include <TLorentzVector.h>

using namespace std;
void pn_xAna(std::string filename="nano_117.root", std::string outputFileName="pn.root"){
           
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
  
  //get TTree from file ...  
  TreeReader data(inputFiles,"Events");
  TTree* thisTree = data.GetTree();

  // check if this tree has the required branches for trigger
  std::vector<std::string> trigNames;
  trigNames.push_back("HLT_PFHT800");
  trigNames.push_back("HLT_PFHT900");
  trigNames.push_back("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5");
  trigNames.push_back("HLT_AK8PFJet360_TrimMass30");
  trigNames.push_back("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20");
  trigNames.push_back("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50");
  trigNames.push_back("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50");

  std::vector<std::string> tempNames = trigNames;

  // first check if the tirgger-branch exists or not
  for(unsigned int itrig=0; itrig < tempNames.size(); itrig++)
    {
      TBranch* thisBranch = thisTree->FindBranch(tempNames[itrig].data());
      if(thisBranch==NULL){
	cerr << "Branch: " << tempNames[itrig] << " is not present in the tree! " << endl;
	trigNames.erase(trigNames.begin()+itrig);
      }
    }

  unsigned int nTrigs= trigNames.size();
  cout << "Available number of trigger paths = " << nTrigs << endl;
  


  
  Long64_t nTotal=0;
  Long64_t nPass[20]={0};
  double ym[17]={70,95,120,145,195,245,295,395,495,595,695,795,895,995,1195,1395,1595};
  Float_t YM[16]={90,100,125,150,200,250,300,400,500,600,700,800,900,1000,1200,1400};//vector of ymass points
  float_t deltaR[24]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4};//VALUE OF DELTA_R
  

  std::vector<float> M_H;
  std::vector<float> M_Y;

  

  const unsigned int nLabels=10;
//  TCanvas *c1 = new TCanvas("c1","c1",3);
  
  TH1F* heve=new TH1F("heve","",nLabels,-0.5,20.5);
  TH1F* m_H = new TH1F("m_H"," ",80,100,140);
  TH1F* m_Y = new TH1F("m_Y"," ",90,140,320);
  

 
  vector<TH1F*> YPN200(24);
  vector<TH1F*> YPN250(24);
  vector<TH1F*> YPN200_d(24);
  vector<TH1F*> YPN250_d(24);
  YPN200[0] = new TH1F("YPN","",100,0,1);

  
  for(int i = 0 ; i < 24 ; i++){
    YPN200[i] = (TH1F*)YPN200[0]->Clone(Form("YPN200_%4.1f",deltaR[i]));
    YPN250[i] = (TH1F*)YPN200[0]->Clone(Form("YPN250_%4.1f",deltaR[i]));
    YPN200_d[i] = (TH1F*)YPN200[0]->Clone(Form("YPN200_d_%4.1f",deltaR[i]));
    YPN250_d[i] = (TH1F*)YPN200[0]->Clone(Form("YPN250_d_%4.1f",deltaR[i]));
  }




  TH1F* HT = new TH1F("HT","",100,0,2500);

  heve->SetYTitle("Number of Events");
  heve->LabelsOption("v");
//  YMass->SetYTitle("Number of YMass");
//  YMass->LabelsOption("v");
  


  //store the variable into the tree called data
  TTree* variable = new TTree("variable", "variable");
  
  Float_t FatJet_Pt[2],FatJet_Eta[2],FatJet_Mass[2],FatJet_Msoftdrop[2],FatJet_BtagHbb[2],FatJet_ParticleNetMD_ProbXbb,FatJet_DeepTagMD_ZHbbvsQCD;
  Float_t FatJet_Phi[2],I_FatJet_ParticleNetMD_ProbXbb,I_FatJet_DeepTagMD_ZHbbvsQCD;
  
  
  variable->Branch("FatJet_Pt",&FatJet_Pt,"FatJet_Pt[2]/F");
  variable->Branch("FatJet_Eta",&FatJet_Eta,"FatJet_Eta[2]/F");
  variable->Branch("FatJet_Phi",&FatJet_Phi,"FatJet_Phi[2]/F");
  variable->Branch("FatJet_Mass",&FatJet_Mass,"FatJet_Mass[2]/F");
  variable->Branch("FatJet_Msoftdrop",&FatJet_Msoftdrop,"FatJet_Msoftdrop[2]/F");
  variable->Branch("FatJet_BtagHbb",&FatJet_BtagHbb,"FatJet_BtagHbb[2]/F");
  variable->Branch("FatJet_ParticleNetMD_ProbXbb",&FatJet_ParticleNetMD_ProbXbb,"FatJet_ParticleNetMD_ProbXbb/F");
  variable->Branch("I_FatJet_ParticleNetMD_ProbXbb",&I_FatJet_ParticleNetMD_ProbXbb,"I_FatJet_ParticleNetMD_ProbXbb/F");
  variable->Branch("I_FatJet_DeepTagMD_ZHbbvsQCD",&I_FatJet_DeepTagMD_ZHbbvsQCD,"I_FatJet_DeepTagMD_ZHbbvsQCD/F");
  variable->Branch("FatJet_DeepTagMD_ZHbbvsQCD",&FatJet_DeepTagMD_ZHbbvsQCD,"FatJet_DeepTagMD_ZHbbvsQCD/F");
  





  const char *label[nLabels];
  label[0]="Total";
  for(unsigned int i=1; i< nLabels; i++){
    label[i] = Form("Cut %d",i);
  }

  // start looping over events
  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){

    if (jEntry % 5000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());

    data.GetEntry(jEntry);
    nTotal++;
    heve->Fill(label[0],1.);


 
    //1. trigger 
    
    
    Bool_t passTrigger=false;
      
    for(unsigned int itrig=0; itrig< nTrigs; itrig++)
      {
	if(data.GetBool(trigNames[itrig].data())==true)
	  {
	    passTrigger=true;
	    break;
	  }
      } // end of loop over required trigger paths

    if(!passTrigger)continue;
    nPass[0]++;
    heve->Fill(label[1],1.);

    // loop over fatjets
    Float_t*  FatJet_pt = data.GetPtrFloat("FatJet_pt");
    Float_t*  FatJet_eta = data.GetPtrFloat("FatJet_eta");
    Float_t*  FatJet_phi = data.GetPtrFloat("FatJet_phi");
    Float_t*  FatJet_tau1 = data.GetPtrFloat("FatJet_tau1");
    Float_t*  FatJet_tau2 = data.GetPtrFloat("FatJet_tau2");
    Float_t*  FatJet_mass = data.GetPtrFloat("FatJet_mass");
    Float_t*  FatJet_msoftdrop = data.GetPtrFloat("FatJet_msoftdrop");
    Float_t  SoftActivityJetHT = data.GetFloat("SoftActivityJetHT");
    Float_t*  FatJet_btagHbb = data.GetPtrFloat("FatJet_btagHbb");
    Float_t*  FatJet_deepTagMD_ZHbbvsQCD = data.GetPtrFloat("FatJet_deepTagMD_ZHbbvsQCD");
    Float_t*  FatJet_ParticleNetMD_probXbb = data.GetPtrFloat("FatJet_ParticleNetMD_probXbb");


    HT->Fill(SoftActivityJetHT);
    UInt_t nFatJet = data.GetInt("nFatJet");
    UInt_t nGoodPair = 0;
    UInt_t nGoodPairJet = 0;
    Float_t Max_pt = 0;
    Float_t sub_Max_pt = 0;
//    Int_t leadingID = 0;
//    Int_t subleadingID = 0;
    Int_t leadingID[24][1] = {0};
    Int_t subleadingID[24][1] = {0};
    Int_t MID = 0;
    
  //  Float_t m_H = 0;

    Bool_t masspass = false;
    Bool_t pt_etapass_sub = false;
    Bool_t pt_etapass = false;
    Bool_t deep = false;
    Bool_t deep_sub = false;
    Bool_t net = false;



// ==========================================================================
// selection
// ==========================================================================
for(UInt_t i = 0 ; i < 24 ; i++){
    for(UInt_t ij=0; ij < nFatJet; ij++){

      if(FatJet_pt[ij] < 200)continue;
      if(fabs(FatJet_eta[ij]) > 2.4)continue;
        pt_etapass_sub = true;

      for(UInt_t jj=0; jj < ij; jj++){

	     if(FatJet_pt[jj] < 300)continue;
	     if(fabs(FatJet_eta[jj]) > 2.4)continue;
            pt_etapass = true;

       

	// eta difference
	     if(fabs(FatJet_eta[ij]-FatJet_eta[jj]) > deltaR[i])continue;
       
         if(FatJet_pt[jj] > Max_pt){
            Max_pt = FatJet_pt[jj];
            leadingID[i][0] = jj;
        }

	     nGoodPair++;
	
      } // end of inner jet loop

      if(FatJet_pt[ij] > sub_Max_pt && sub_Max_pt < Max_pt){
        sub_Max_pt = FatJet_pt[ij] ;
        subleadingID[i][0] = ij;
      
      }
         
    } // end of outer jet loop
}// detalR loop

    if(pt_etapass_sub){
      nPass[1]++;
      heve->Fill(label[2],1.);
      
  }  
    if(pt_etapass){
      nPass[2]++;
      heve->Fill(label[3],1.);
      
  }
    if(nGoodPair<1)continue;
    nPass[3]++;
    heve->Fill(label[4],1.);
    
for(UInt_t i = 0 ; i < 24 ; i++){
    if(100 < FatJet_msoftdrop[leadingID[i][0]] && FatJet_msoftdrop[leadingID[i][0]] < 140 && 160 < FatJet_msoftdrop[subleadingID[i][0]]){
      M_H.push_back(FatJet_msoftdrop[leadingID[i][0]]);
      M_Y.push_back(FatJet_msoftdrop[subleadingID[i][0]]);
      nGoodPairJet++;
      masspass = true;
      if(FatJet_pt[leadingID[i][0]] > 300 && FatJet_pt[subleadingID[i][0]] > 200){
//        I_FatJet_ParticleNetMD_ProbXbb = FatJet_ParticleNetMD_probXbb[subleadingID] ; 
        if(160 < FatJet_msoftdrop[subleadingID[i][0]] && FatJet_msoftdrop[subleadingID[i][0]] < 250){
        YPN200[i]->Fill(FatJet_ParticleNetMD_probXbb[subleadingID[i][0]]); 
        }
        if(200 < FatJet_msoftdrop[subleadingID[i][0]] && FatJet_msoftdrop[subleadingID[i][0]] < 300){
        YPN250[i]->Fill(FatJet_ParticleNetMD_probXbb[subleadingID[i][0]]); 
        }
        
        //FOR DEEPTAG selection 
/*        if(FatJet_deepTagMD_ZHbbvsQCD[subleadingID] > 0.8){
          for(int iymass=0; iymass < nYmass; iymass++){
            if(data.GetBool(YmassNames[iymass].data())==true){
              YPN_d[iymass]->Fill(FatJet_ParticleNetMD_probXbb[subleadingID]);
            }
          }
        }  //end of deeptag then see the distribution of particle net
*/
      }//end of pt selection
    }//end of masspass
     
    

    if(100 < FatJet_msoftdrop[subleadingID[i][0]] && FatJet_msoftdrop[subleadingID[i][0]] < 140 && 160 < FatJet_msoftdrop[leadingID[i][0]]){
      M_H.push_back(FatJet_msoftdrop[subleadingID[i][0]]);
      M_Y.push_back(FatJet_msoftdrop[leadingID[i][0]]);
      nGoodPairJet++;
      masspass = true;
      if(FatJet_pt[leadingID[i][0]] > 200 && FatJet_pt[subleadingID[i][0]] > 300){
//        I_FatJet_ParticleNetMD_ProbXbb = FatJet_ParticleNetMD_probXbb[leadingID] ; 
        if(160 < FatJet_msoftdrop[leadingID[i][0]] && FatJet_msoftdrop[leadingID[i][0]] < 250){
        YPN200[i]->Fill(FatJet_ParticleNetMD_probXbb[leadingID[i][0]]);
        }
        if(200 < FatJet_msoftdrop[leadingID[i][0]] && FatJet_msoftdrop[leadingID[i][0]] < 300){
        YPN250[i]->Fill(FatJet_ParticleNetMD_probXbb[leadingID[i][0]]);
        }
        
//FOR DEEPTAG selection
/*
        if(FatJet_deepTagMD_ZHbbvsQCD[leadingID] > 0.8){
          for(int iymass=0; iymass < nYmass; iymass++){
            if(data.GetBool(YmassNames[iymass].data())==true){
              YPN_d[iymass]->Fill(FatJet_ParticleNetMD_probXbb[leadingID]);
            }
          }
        }//end of DEEPTAG SELECTION
*/   
      }//END OF PT
   }//END OF MASSPASS
}//deltaR loop    



// store the leading & subleading information into the vector(tree)
/*
    FatJet_Pt[0] = FatJet_pt[leadingID];
    FatJet_Pt[1] = FatJet_pt[subleadingID];
    FatJet_Eta[0] = FatJet_eta[leadingID];
    FatJet_Eta[1] = FatJet_eta[subleadingID];
    FatJet_Phi[0]=FatJet_phi[leadingID];
    FatJet_Phi[1]=FatJet_phi[subleadingID];  
    
    FatJet_Mass[0] = FatJet_mass[leadingID];
    FatJet_Mass[1] = FatJet_mass[subleadingID];

    FatJet_Msoftdrop[0] = FatJet_msoftdrop[leadingID];
    FatJet_Msoftdrop[1] = FatJet_msoftdrop[subleadingID];
    FatJet_BtagHbb[0] = FatJet_btagHbb[leadingID];
    FatJet_BtagHbb[1] = FatJet_btagHbb[subleadingID];
    
 */   
    if(masspass == true){
      nPass[4]++;
      heve->Fill(label[5],1.);
  }


    variable->Fill();


    

    

  } // end of loop over events

  

  for(int i=0;i<20;i++)
    if(nPass[i]>0)cout << "nPass["<<i<<"]= " << nPass[i] << endl;



  for(int i = 0 ; i < nPass[4] ; i++ ){
    m_H->Fill( M_H[i] );
    m_Y->Fill( M_Y[i] );

  }

  TFile* outFile = new TFile(outputFileName.data(),"recreate");

  
  heve->Write();
  m_H->Write();
  m_Y->Write();
  for(int i = 0 ; i< 24 ; i++){
  YPN200[i]->Write();
  YPN250[i]->Write();
  }

  
  
  

  HT->Write();
  variable->Write();
  outFile->Close();
  
  
  
  
 } // end of macro
 