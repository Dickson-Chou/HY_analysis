// .L xAna_nano_nssm.C++
// xAna_nano_nssm("test.root") or xAna_nano_nssm("input.txt")

//======================================================
//update:2021/6/22
//focusing on specific mass range
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
void old_xAna(std::string filename="devdatta_nanoAOD_nssm.root", std::string outputFileName="histo.root"){
           
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
  
  // second check if the Ymass-brance exists or not
  std::vector<std::string> YmassNames;
  YmassNames.push_back("GenModel_YMass_90");
  YmassNames.push_back("GenModel_YMass_100");
  YmassNames.push_back("GenModel_YMass_125");
  YmassNames.push_back("GenModel_YMass_150");
  YmassNames.push_back("GenModel_YMass_200");
  YmassNames.push_back("GenModel_YMass_250");
  YmassNames.push_back("GenModel_YMass_300");
  YmassNames.push_back("GenModel_YMass_400");
  YmassNames.push_back("GenModel_YMass_500");
  YmassNames.push_back("GenModel_YMass_600");
  YmassNames.push_back("GenModel_YMass_700");
  YmassNames.push_back("GenModel_YMass_800");
  YmassNames.push_back("GenModel_YMass_900");
  YmassNames.push_back("GenModel_YMass_1000");
  YmassNames.push_back("GenModel_YMass_1200");
  YmassNames.push_back("GenModel_YMass_1400");


  std::vector<std::string> TYmassNames = YmassNames;

  for(unsigned int iymass=0; iymass < TYmassNames.size(); iymass++)
  {
    TBranch* massBranch = thisTree->FindBranch(TYmassNames[iymass].data());
    if(massBranch==NULL){
  cerr << "Branch: " << TYmassNames[iymass] << " is not present in the tree! " << endl;
  YmassNames.erase(YmassNames.begin()+iymass);
    }
  }

  unsigned int nYmass = YmassNames.size();
  cout << "Available number of Ymass value = " << nYmass << endl;



  
  Long64_t nTotal=0;
  Long64_t nPass[20]={0};
  double ym[17]={70,95,120,145,195,245,295,395,495,595,695,795,895,995,1195,1395,1595};
  Float_t YM[16]={90,100,125,150,200,250,300,400,500,600,700,800,900,1000,1200,1400};//vector of ymass points
  Float_t eta[24]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4};//delta eta value
  std::vector<float> npassym(16,0);//vector of the number of ym in each mass.
  std::vector<float> npassym_initial(16,0);//vector of the number of ym without selection
  std::vector<float> npassym_P1(16,0);
  std::vector<float> npassym_P2(16,0);
  std::vector<float> npassym_P3(16.0);
  std::vector<float> npassym_P4(16,0);
  std::vector<float> npassym_P5(16,0);
  std::vector<float> npassym_P6(16,0);
  std::vector<float> npassym_P7(16,0);
  std::vector<float> npassym_P8(16,0);
  std::vector<float> npassym_P9(16,0);
  std::vector<float> npassym_P10(16,0);
  std::vector<float> npassym_P11(16,0);
  std::vector<float> npassym_P12(16,0);
  

  std::vector<float> M_H;
  std::vector<float> M_Y;
  

  

  const unsigned int nLabels=10;
//  TCanvas *c1 = new TCanvas("c1","c1",3);
  
  TH1F* heve=new TH1F("heve","",nLabels,-0.5,20.5);
  TH1F* m_H = new TH1F("m_H"," ",80,100,140);
  TH1F* m_Y = new TH1F("m_Y"," ",110,140,360);
 // TH1F* H_FatJet_ParticleNetMD_probXbb = new TH1F("H_FatJet_ParticleNetMD_probXbb"," ",100,0,1);

//  TH1F* m_Y = new TH1F("m_Y"," ",)
//  TH1F* YMass = new TH1F("YMass","",2,0,2);
 
  
  
//==============================================
//histogram of particle net with different y mass range
//==============================================


  vector<TH1F*> YPN_M90(24);
  vector<TH1F*> YPN_M100(24);
  vector<TH1F*> YPN_M125(24);
  vector<TH1F*> YPN_M150(24);
  vector<TH1F*> YPN_M200(24);
  vector<TH1F*> YPN_M250(24);
  vector<TH1F*> YPN_M300(24);
  vector<TH1F*> YPN_M400(24);
  vector<TH1F*> YPN_M500(24);
  vector<TH1F*> YPN_M600(24);
  vector<TH1F*> YPN_M700(24);
  vector<TH1F*> YPN_M800(24);
  vector<TH1F*> YPN_M900(24);
  vector<TH1F*> YPN_M1000(24);
  vector<TH1F*> YPN_M1200(24);
  vector<TH1F*> YPN_M1400(24);

  vector<TH1F*> YM200(24);
  vector<TH1F*> YM250(24);
  vector<TH1F*> YM300(24);
  vector<TH1F*> YM200_l(24);
  vector<TH1F*> YM250_l(24);
  vector<TH1F*> YM300_l(24);
  vector<TH1F*> YM200_s(24);
  vector<TH1F*> YM250_s(24);
  vector<TH1F*> YM300_s(24);
  vector<TH1F*> YM200_p(24);
  vector<TH1F*> YM250_p(24);
  vector<TH1F*> YM300_p(24);
  

  YPN_M90[0] = new TH1F("YPN","",100,0,1);
  YM200[0] = new TH1F("ym","",110,140,360);
  
  for(int i = 0 ; i < 24 ; i++){
    YPN_M90[i] = (TH1F*)YPN_M90[0]->Clone(Form("signal90_%4.1f",eta[i]));
    YPN_M100[i] = (TH1F*)YPN_M90[0]->Clone(Form("signal100_%4.1f",eta[i]));
    YPN_M125[i] = (TH1F*)YPN_M90[0]->Clone(Form("signal125_%4.1f",eta[i]));
    YPN_M150[i] = (TH1F*)YPN_M90[0]->Clone(Form("signal150_%4.1f",eta[i]));
    YPN_M200[i] = (TH1F*)YPN_M90[0]->Clone(Form("signal200_%4.1f",eta[i]));
    YPN_M250[i] = (TH1F*)YPN_M90[0]->Clone(Form("signal250_%4.1f",eta[i]));
    YPN_M300[i] = (TH1F*)YPN_M90[0]->Clone(Form("signal300_%4.1f",eta[i]));
    YPN_M400[i] = (TH1F*)YPN_M90[0]->Clone(Form("signal400_%4.1f",eta[i]));
    YPN_M500[i] = (TH1F*)YPN_M90[0]->Clone(Form("signal500_%4.1f",eta[i]));
    YPN_M600[i] = (TH1F*)YPN_M90[0]->Clone(Form("signal600_%4.1f",eta[i]));
    YPN_M700[i] = (TH1F*)YPN_M90[0]->Clone(Form("signal700_%4.1f",eta[i]));
    YPN_M800[i] = (TH1F*)YPN_M90[0]->Clone(Form("signal800_%4.1f",eta[i]));
    YPN_M900[i] = (TH1F*)YPN_M90[0]->Clone(Form("signal900_%4.1f",eta[i]));
    YPN_M1000[i] = (TH1F*)YPN_M90[0]->Clone(Form("signal1000_%4.1f",eta[i]));
    YPN_M1200[i] = (TH1F*)YPN_M90[0]->Clone(Form("signal1200_%4.1f",eta[i]));
    YPN_M1400[i] = (TH1F*)YPN_M90[0]->Clone(Form("signal1400_%4.1f",eta[i]));
    YM200[i] = (TH1F*)YM200[0]->Clone(Form("YM200_%4.1f",eta[i]));
    YM250[i] = (TH1F*)YM200[0]->Clone(Form("YM250_%4.1f",eta[i]));
    YM300[i] = (TH1F*)YM200[0]->Clone(Form("YM300_%4.1f",eta[i]));
    YM200_l[i] = (TH1F*)YM200[0]->Clone(Form("YM200_l_%4.1f",eta[i]));
    YM250_l[i] = (TH1F*)YM200[0]->Clone(Form("YM250_l_%4.1f",eta[i]));
    YM300_l[i] = (TH1F*)YM200[0]->Clone(Form("YM300_l_%4.1f",eta[i]));
    YM200_s[i] = (TH1F*)YM200[0]->Clone(Form("YM200_s_%4.1f",eta[i]));
    YM250_s[i] = (TH1F*)YM200[0]->Clone(Form("YM250_s_%4.1f",eta[i]));
    YM300_s[i] = (TH1F*)YM200[0]->Clone(Form("YM300_s_%4.1f",eta[i]));
    YM200_p[i] = (TH1F*)YM200[0]->Clone(Form("YM200_p_%4.1f",eta[i]));
    YM250_p[i] = (TH1F*)YM200[0]->Clone(Form("YM250_p_%4.1f",eta[i]));
    YM300_p[i] = (TH1F*)YM200[0]->Clone(Form("YM300_p_%4.1f",eta[i]));
  }

  TH1F* DDE = new TH1F("DDE","",100,0,10);


  TH1F* HT = new TH1F("HT","",100,0,2500);

  heve->SetYTitle("Number of Events");
  heve->LabelsOption("v");
//  YMass->SetYTitle("Number of YMass");
//  YMass->LabelsOption("v");
  


  //store the variable into the tree called data
  TTree* variable = new TTree("variable", "variable");
  
  Float_t FatJet_Pt[2],FatJet_Eta[2],FatJet_Mass[2],FatJet_Msoftdrop[2],FatJet_BtagHbb[2],FatJet_ParticleNetMD_ProbXbb,FatJet_DeepTagMD_ZHbbvsQCD;
  Float_t FatJet_Phi[2],I_FatJet_ParticleNetMD_ProbXbb,I_FatJet_DeepTagMD_ZHbbvsQCD;
  Bool_t GenModel_YMAss_150,GenModel_YMAss_90,GenModel_YMAss_100,GenModel_YMAss_125,GenModel_YMAss_200,GenModel_YMAss_250,GenModel_YMAss_300,GenModel_YMAss_400;
  Bool_t GenModel_YMAss_500,GenModel_YMAss_600,GenModel_YMAss_700,GenModel_YMAss_800,GenModel_YMAss_900,GenModel_YMAss_1000,GenModel_YMAss_1200,GenModel_YMAss_1400;
  
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
  variable->Branch("GenModel_YMAss_90",&GenModel_YMAss_90,"GenModel_YMAss_90/O");
  variable->Branch("GenModel_YMAss_100",&GenModel_YMAss_100,"GenModel_YMAss_100/O");
  variable->Branch("GenModel_YMAss_125",&GenModel_YMAss_125,"GenModel_YMAss_125/O");
  variable->Branch("GenModel_YMAss_150",&GenModel_YMAss_150,"GenModel_YMAss_150/O");
  variable->Branch("GenModel_YMAss_200",&GenModel_YMAss_200,"GenModel_YMAss_200/O");
  variable->Branch("GenModel_YMAss_250",&GenModel_YMAss_250,"GenModel_YMAss_250/O");
  variable->Branch("GenModel_YMAss_300",&GenModel_YMAss_300,"GenModel_YMAss_300/O"); 
  variable->Branch("GenModel_YMAss_400",&GenModel_YMAss_400,"GenModel_YMAss_400/O");
  variable->Branch("GenModel_YMAss_500",&GenModel_YMAss_500,"GenModel_YMAss_500/O");
  variable->Branch("GenModel_YMAss_600",&GenModel_YMAss_600,"GenModel_YMAss_600/O");
  variable->Branch("GenModel_YMAss_700",&GenModel_YMAss_700,"GenModel_YMAss_700/O");
  variable->Branch("GenModel_YMAss_800",&GenModel_YMAss_800,"GenModel_YMAss_800/O");
  variable->Branch("GenModel_YMAss_900",&GenModel_YMAss_900,"GenModel_YMAss_900/O");
  variable->Branch("GenModel_YMAss_1000",&GenModel_YMAss_1000,"GenModel_YMAss_1000/O");
  variable->Branch("GenModel_YMAss_1200",&GenModel_YMAss_1200,"GenModel_YMAss_1200/O");
  variable->Branch("GenModel_YMAss_1400",&GenModel_YMAss_1400,"GenModel_YMAss_1400/O");





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


  // store the initial number of the y into the vector 
  for(unsigned int i = 0 ; i < nYmass ; i++ ){   
    if(data.GetBool(YmassNames[i].data())==true){
      npassym_initial[i] += 1;
     }
    }
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
    Int_t* GenPart_genPartIdxMother = data.GetPtrInt("GenPart_genPartIdxMother");
    Int_t* GenPart_pdgId = data.GetPtrInt("GenPart_pdgId");
    Float_t* GenPart_pt = data.GetPtrFloat("GenPart_pt");
    Float_t* GenPart_eta = data.GetPtrFloat("GenPart_eta");
    Float_t* GenPart_phi = data.GetPtrFloat("GenPart_phi");
    Float_t* GenPart_mass = data.GetPtrFloat("GenPart_mass");

    Bool_t GenModel_YMass_200 = data.GetBool("GenModel_YMass_200");
    Bool_t GenModel_YMass_250 = data.GetBool("GenModel_YMass_250");
    Bool_t GenModel_YMass_300 = data.GetBool("GenModel_YMass_300");
/*
    Bool_t GenModel_YMass_90 = data.GetBool("GenModel_YMass_90");
    Bool_t GenModel_YMass_100 = data.GetBool("GenModel_YMass_100");
    Bool_t GenModel_YMass_125 = data.GetBool("GenModel_YMass_125");
    Bool_t GenModel_YMass_150 = data.GetBool("GenModel_YMass_150");
    Bool_t GenModel_YMass_200 = data.GetBool("GenModel_YMass_200");
    Bool_t GenModel_YMass_250 = data.GetBool("GenModel_YMass_250");
    Bool_t GenModel_YMass_300 = data.GetBool("GenModel_YMass_300");
    Bool_t GenModel_YMass_400 = data.GetBool("GenModel_YMass_400");
    Bool_t GenModel_YMass_500 = data.GetBool("GenModel_YMass_500");
    Bool_t GenModel_YMass_600 = data.GetBool("GenModel_YMass_600");
    Bool_t GenModel_YMass_700 = data.GetBool("GenModel_YMass_700");
    Bool_t GenModel_YMass_800 = data.GetBool("GenModel_YMass_800");
    Bool_t GenModel_YMass_900 = data.GetBool("GenModel_YMass_900");
    Bool_t GenModel_YMass_1000 = data.GetBool("GenModel_YMass_1000");
    Bool_t GenModel_YMass_1200 = data.GetBool("GenModel_YMass_1200");
    Bool_t GenModel_YMass_1400 = data.GetBool("GenModel_YMass_1400");
    Bool_t GenModel_YMass_123 = data.GetBool("GenModel_YMass_123"); */
/*
    std::vector<bool> v;
    for(int i=0; i<YmassNames.size();i++){
      
      Bool_t tmp = data.GetBool(YmassNames[i].data());
      v.push_back(tmp);
    
  }
    for(int i = 0 ; i < YmassNames.size();i++){
      cout << v[i] << " " << YmassNames[i] << endl;
    }
    // bool_t .. = data.GetBool("");

*/
    HT->Fill(SoftActivityJetHT);
    UInt_t nFatJet = data.GetInt("nFatJet");
    UInt_t nGenPart = data.GetInt("nGenPart");
    UInt_t nGoodPair = 0;
    UInt_t nGoodPairJet = 0;
    Float_t Max_pt = 0;
    Float_t sub_Max_pt = 0;
    Int_t leadingID[24][1] = {0};
    Int_t subleadingID[24][1] = {0};
  //  Int_t leadingID = 0;
  //  Int_t subleadingID = 0;
    Int_t MID = 0;
    
    
  //  Float_t m_H = 0;

    Bool_t masspass = false;
    Bool_t pt_etapass_sub = false;
    Bool_t pt_etapass = false;
    Bool_t deep = false;
    Bool_t deep_sub = false;
    Bool_t net = false;
    Bool_t tau21pass_sub = false;
    Bool_t tau21pass = false;

     //TLorentz
    vector<vector<TLorentzVector>> v(2);
    std::vector<float> deltaR;
    TLorentzVector v_H;
    TLorentzVector v_Y;

// =========================================================================
// matching
// =========================================================================

    for(UInt_t ij=0; ij < nGenPart ; ij++){
      
      if(abs(GenPart_pdgId[ij]) == 5){
        MID = GenPart_genPartIdxMother[ij];
        
        if(GenPart_pdgId[MID] == 25){
          v_H.SetPtEtaPhiM(GenPart_pt[ij],GenPart_eta[ij],GenPart_phi[ij],GenPart_mass[ij]);
          v[0].push_back(v_H);
        }
        else if(GenPart_pdgId[MID] == 35){
          v_Y.SetPtEtaPhiM(GenPart_pt[ij],GenPart_eta[ij],GenPart_phi[ij],GenPart_mass[ij]);
          v[1].push_back(v_Y);
        }

      }

    }
    deltaR.push_back(v[0][0].DeltaR(v[0][1]));//Hbb
    deltaR.push_back(v[1][0].DeltaR(v[1][1]));//Ybb
    //======HbYb
    deltaR.push_back(v[0][0].DeltaR(v[1][0]));
    deltaR.push_back(v[0][0].DeltaR(v[1][1]));
    deltaR.push_back(v[0][1].DeltaR(v[1][0]));
    deltaR.push_back(v[0][1].DeltaR(v[1][1]));
/*
    cout << v[0][0].Eta()-v[0][1].Eta() << endl;
    cout << v[0][0].Phi()-v[0][1].Phi() << endl;

    for(int i = 0 ; i < 6 ; i++){
      cout << "deltaR" << deltaR[i] << endl;
//      cout << v[0][0].Pt() << endl;
    }
*/

    if(deltaR[0] > 0.8 && deltaR[1] > 0.8)continue;

// ==========================================================================
// selection
// ==========================================================================
  for(UInt_t i = 0 ; i < 24 ; i++){
    for(UInt_t ij=0; ij < nFatJet; ij++){

      if(FatJet_pt[ij] < 200)continue;
      if(fabs(FatJet_eta[ij]) > 2.4)continue;
        pt_etapass_sub = true;
      
//      if(!(FatJet_deepTagMD_ZHbbvsQCD[ij] > 0.8 && FatJet_ParticleNetMD_probXbb[ij] > 0.85))continue;
//      if(FatJet_ParticleNetMD_probXbb[ij] > 0.85){
//      net_sub = true;
//    }

      for(UInt_t jj=0; jj < ij; jj++){

	     if(FatJet_pt[jj] < 300)continue;
	     if(fabs(FatJet_eta[jj]) > 2.4)continue;
       pt_etapass = true;
	     
//       if(!(FatJet_deepTagMD_ZHbbvsQCD[jj] > 0.8 && FatJet_ParticleNetMD_probXbb[jj] > 0.85))continue;
//       if(FatJet_ParticleNetMD_probXbb[jj] > 0.85){
//       net = true;
//     }
       

// eta difference
// when doing N-1 for delta eta  
       DDE->Fill(fabs(FatJet_eta[ij]-FatJet_eta[jj]));
	     if(fabs(FatJet_eta[ij]-FatJet_eta[jj]) > eta[i])continue;
       
       if(FatJet_pt[jj] > Max_pt){
        Max_pt = FatJet_pt[jj];
        leadingID[i][0] = jj;
      }
/*      
       if(jj == leadingID[i][0]){
        FatJet_tau21_j = tau21_j;
      }
*/	
	     nGoodPair++;
	
      } // end of inner jet loop

      if(FatJet_pt[ij] > sub_Max_pt && sub_Max_pt < Max_pt){
        sub_Max_pt = FatJet_pt[ij] ;
        subleadingID[i][0] = ij;
      
      }
         
    } // end of outer jet loop
  }//end of different eta

    if(pt_etapass_sub){
      nPass[1]++;
      heve->Fill(label[2],1.);
      for(unsigned int iymass=0 ; iymass < nYmass; iymass++){
        if(data.GetBool(YmassNames[iymass].data())==true){
        npassym_P1[iymass]+=1;
      }
    }
  }  
    if(pt_etapass){
      nPass[2]++;
      heve->Fill(label[3],1.);
      for(unsigned int iymass=0 ; iymass < nYmass; iymass++){
        if(data.GetBool(YmassNames[iymass].data())==true){
        npassym_P2[iymass]+=1;
      }
    }
  }
    if(nGoodPair<1)continue;
    nPass[3]++;
    heve->Fill(label[4],1.);
    for(unsigned int iymass=0 ; iymass < nYmass; iymass++){
        if(data.GetBool(YmassNames[iymass].data())==true){
        npassym_P3[iymass]+=1;
      }
    }
  
  for(UInt_t i = 0 ; i < 24 ; i++){

     for(int iymass=0; iymass < nYmass; iymass++){
      if(data.GetBool(YmassNames[4].data())==true){
            YM200_l[i]->Fill(FatJet_msoftdrop[leadingID[i][0]]);
            YM200_s[i]->Fill(FatJet_msoftdrop[subleadingID[i][0]]);
          }
      if(data.GetBool(YmassNames[5].data())==true){
            YM250_l[i]->Fill(FatJet_msoftdrop[leadingID[i][0]]);
            YM250_s[i]->Fill(FatJet_msoftdrop[subleadingID[i][0]]);
          } 
      if(data.GetBool(YmassNames[6].data())==true){
            YM300_l[i]->Fill(FatJet_msoftdrop[leadingID[i][0]]);
            YM300_s[i]->Fill(FatJet_msoftdrop[subleadingID[i][0]]);
          }  
    }      
    if(100 < FatJet_msoftdrop[leadingID[i][0]] && FatJet_msoftdrop[leadingID[i][0]] < 140 &&  160 < FatJet_msoftdrop[subleadingID[i][0]] && FatJet_msoftdrop[subleadingID[i][0]] < 250){
      nGoodPairJet++;
      masspass = true;
      if(FatJet_pt[leadingID[i][0]] > 300 && FatJet_pt[subleadingID[i][0]] > 200){
        for(int iymass=0; iymass < nYmass; iymass++){
          if(data.GetBool(YmassNames[4].data())==true){
            YPN_M200[i]->Fill(FatJet_ParticleNetMD_probXbb[subleadingID[i][0]]);
            YM200[i]->Fill(FatJet_msoftdrop[subleadingID[i][0]]);
            if(FatJet_ParticleNetMD_probXbb[subleadingID[i][0]] > 0.8){
              YM200_p[i]->Fill(FatJet_msoftdrop[subleadingID[i][0]]);
            }
          }
        }
      }
    }
    if(100 < FatJet_msoftdrop[subleadingID[i][0]] && FatJet_msoftdrop[subleadingID[i][0]] < 140 &&  160 < FatJet_msoftdrop[leadingID[i][0]] && FatJet_msoftdrop[leadingID[i][0]] < 250){
      nGoodPairJet++;
      masspass = true;
      if(FatJet_pt[subleadingID[i][0]] > 300 && FatJet_pt[leadingID[i][0]] > 200){
        for(int iymass=0; iymass < nYmass; iymass++){
          if(data.GetBool(YmassNames[4].data())==true){
            YPN_M200[i]->Fill(FatJet_ParticleNetMD_probXbb[leadingID[i][0]]);
            YM200[i]->Fill(FatJet_msoftdrop[leadingID[i][0]]);
            if(FatJet_ParticleNetMD_probXbb[leadingID[i][0]] > 0.8){
              YM200_p[i]->Fill(FatJet_msoftdrop[leadingID[i][0]]);
            }
          }
        }
      }
    }
    if(100 < FatJet_msoftdrop[leadingID[i][0]] && FatJet_msoftdrop[leadingID[i][0]] < 140 &&  200 < FatJet_msoftdrop[subleadingID[i][0]] && FatJet_msoftdrop[subleadingID[i][0]] < 300){
      nGoodPairJet++;
      masspass = true;
      if(FatJet_pt[leadingID[i][0]] > 300 && FatJet_pt[subleadingID[i][0]] > 200){
        for(int iymass=0; iymass < nYmass; iymass++){
          if(data.GetBool(YmassNames[5].data())==true){
            YPN_M250[i]->Fill(FatJet_ParticleNetMD_probXbb[subleadingID[i][0]]);
            YM250[i]->Fill(FatJet_msoftdrop[subleadingID[i][0]]);
            if(FatJet_ParticleNetMD_probXbb[subleadingID[i][0]] > 0.8){
              YM250_p[i]->Fill(FatJet_msoftdrop[subleadingID[i][0]]);
            }
          }
        }
      }
    }
    if(100 < FatJet_msoftdrop[subleadingID[i][0]] && FatJet_msoftdrop[subleadingID[i][0]] < 140 &&  200 < FatJet_msoftdrop[leadingID[i][0]] && FatJet_msoftdrop[leadingID[i][0]] < 300){
      nGoodPairJet++;
      masspass = true;
      if(FatJet_pt[subleadingID[i][0]] > 300 && FatJet_pt[leadingID[i][0]] > 200){
        for(int iymass=0; iymass < nYmass; iymass++){
          if(data.GetBool(YmassNames[5].data())==true){
            YPN_M250[i]->Fill(FatJet_ParticleNetMD_probXbb[leadingID[i][0]]);
            YM250[i]->Fill(FatJet_msoftdrop[leadingID[i][0]]);
            if(FatJet_ParticleNetMD_probXbb[leadingID[i][0]] > 0.8){
              YM250_p[i]->Fill(FatJet_msoftdrop[leadingID[i][0]]);
            }
          }
        }
      }
    }

    
/*    if(100 < FatJet_msoftdrop[leadingID[i][0]] && FatJet_msoftdrop[leadingID[i][0]] < 140 &&  160 < FatJet_msoftdrop[subleadingID[i][0]] ){
      M_H.push_back(FatJet_msoftdrop[leadingID[i][0]]);
      M_Y.push_back(FatJet_msoftdrop[subleadingID[i][0]]);
      nGoodPairJet++;
      masspass = true;
      if(FatJet_pt[leadingID[i][0]] > 300 && FatJet_pt[subleadingID[i][0]] > 200){
        for(int iymass=0; iymass < nYmass; iymass++){
          if(data.GetBool(YmassNames[0].data())==true){
            YPN_M90[i]->Fill(FatJet_ParticleNetMD_probXbb[subleadingID[i][0]]);
          }
          if(data.GetBool(YmassNames[1].data())==true){
            YPN_M100[i]->Fill(FatJet_ParticleNetMD_probXbb[subleadingID[i][0]]);
          }
          if(data.GetBool(YmassNames[2].data())==true){
            YPN_M125[i]->Fill(FatJet_ParticleNetMD_probXbb[subleadingID[i][0]]);
          }
          if(data.GetBool(YmassNames[3].data())==true){
            YPN_M150[i]->Fill(FatJet_ParticleNetMD_probXbb[subleadingID[i][0]]);
          }  
          if(data.GetBool(YmassNames[4].data())==true){
            if( 160 < FatJet_msoftdrop[subleadingID[i][0]] < 250 ){
              YPN_M200[i]->Fill(FatJet_ParticleNetMD_probXbb[subleadingID[i][0]]);
              YM200[i]->Fill(FatJet_msoftdrop[subleadingID[i][0]]);
              if(FatJet_ParticleNetMD_probXbb[subleadingID[i][0]] > 0.8){
                YM200_p[i]->Fill(FatJet_msoftdrop[subleadingID[i][0]]);  
              }
            }
          }
          if(data.GetBool(YmassNames[5].data())==true){
            if( 200 < FatJet_msoftdrop[subleadingID[i][0]] < 300){
              YPN_M250[i]->Fill(FatJet_ParticleNetMD_probXbb[subleadingID[i][0]]);
              YM250[i]->Fill(FatJet_msoftdrop[subleadingID[i][0]]);
              if(FatJet_ParticleNetMD_probXbb[subleadingID[i][0]] > 0.8){
                YM250_p[i]->Fill(FatJet_msoftdrop[subleadingID[i][0]]);  
              }
            }
          }
          if(data.GetBool(YmassNames[6].data())==true){
            YPN_M300[i]->Fill(FatJet_ParticleNetMD_probXbb[subleadingID[i][0]]);
            YM300[i]->Fill(FatJet_msoftdrop[subleadingID[i][0]]);
            if(FatJet_ParticleNetMD_probXbb[subleadingID[i][0]] > 0.8){
              YM300_p[i]->Fill(FatJet_msoftdrop[subleadingID[i][0]]);  
            }
          }
          if(data.GetBool(YmassNames[7].data())==true){
            YPN_M400[i]->Fill(FatJet_ParticleNetMD_probXbb[subleadingID[i][0]]);
          }
          if(data.GetBool(YmassNames[8].data())==true){
            YPN_M500[i]->Fill(FatJet_ParticleNetMD_probXbb[subleadingID[i][0]]);
          }
          if(data.GetBool(YmassNames[9].data())==true){
            YPN_M600[i]->Fill(FatJet_ParticleNetMD_probXbb[subleadingID[i][0]]);
          }
          if(data.GetBool(YmassNames[10].data())==true){
            YPN_M700[i]->Fill(FatJet_ParticleNetMD_probXbb[subleadingID[i][0]]);
          }
          if(data.GetBool(YmassNames[11].data())==true){
            YPN_M800[i]->Fill(FatJet_ParticleNetMD_probXbb[subleadingID[i][0]]);
          }
          if(data.GetBool(YmassNames[12].data())==true){
            YPN_M900[i]->Fill(FatJet_ParticleNetMD_probXbb[subleadingID[i][0]]);
          }
          if(data.GetBool(YmassNames[13].data())==true){
            YPN_M1000[i]->Fill(FatJet_ParticleNetMD_probXbb[subleadingID[i][0]]);
          }
          if(data.GetBool(YmassNames[14].data())==true){
            YPN_M1200[i]->Fill(FatJet_ParticleNetMD_probXbb[subleadingID[i][0]]);
          }
          if(data.GetBool(YmassNames[15].data())==true){
            YPN_M1400[i]->Fill(FatJet_ParticleNetMD_probXbb[subleadingID[i][0]]);
          }
         
        }//FOR OF Y DIFFERENT MASS RANGE

      }//end of pt selection
    }//end of masspass
*/

   
    
/*    if(100 < FatJet_msoftdrop[subleadingID[i][0]] && FatJet_msoftdrop[subleadingID[i][0]] < 140 &&  160 < FatJet_msoftdrop[leadingID[i][0]]  ){
      M_H.push_back(FatJet_msoftdrop[subleadingID[i][0]]);
      M_Y.push_back(FatJet_msoftdrop[leadingID[i][0]]);
      nGoodPairJet++;
      masspass = true;
      if(FatJet_pt[leadingID[i][0]] > 200 && FatJet_pt[subleadingID[i][0]] > 300){
        for(int iymass=0; iymass < nYmass; iymass++){
          if(data.GetBool(YmassNames[0].data())==true){
            YPN_M90[i]->Fill(FatJet_ParticleNetMD_probXbb[leadingID[i][0]]);
          }
          if(data.GetBool(YmassNames[1].data())==true){
            YPN_M100[i]->Fill(FatJet_ParticleNetMD_probXbb[leadingID[i][0]]);
          }
          if(data.GetBool(YmassNames[2].data())==true){
            YPN_M125[i]->Fill(FatJet_ParticleNetMD_probXbb[leadingID[i][0]]);
          }
          if(data.GetBool(YmassNames[3].data())==true){
            YPN_M150[i]->Fill(FatJet_ParticleNetMD_probXbb[leadingID[i][0]]);
          }  
          if(data.GetBool(YmassNames[4].data())==true){
            if( 160 < FatJet_msoftdrop[leadingID[i][0]] < 250){
              YPN_M200[i]->Fill(FatJet_ParticleNetMD_probXbb[leadingID[i][0]]);
              YM200[i]->Fill(FatJet_msoftdrop[leadingID[i][0]]);
              if(FatJet_ParticleNetMD_probXbb[leadingID[i][0]] > 0.8){
                YM200_p[i]->Fill(FatJet_msoftdrop[leadingID[i][0]]);  
              }
            }
          }
          if(data.GetBool(YmassNames[5].data())==true){
            if( 200 < FatJet_msoftdrop[leadingID[i][0]] < 300){
              YPN_M250[i]->Fill(FatJet_ParticleNetMD_probXbb[leadingID[i][0]]);
              YM250[i]->Fill(FatJet_msoftdrop[leadingID[i][0]]);
              if(FatJet_ParticleNetMD_probXbb[leadingID[i][0]] > 0.8){
                YM250_p[i]->Fill(FatJet_msoftdrop[leadingID[i][0]]);  
              }
            }
          }
          if(data.GetBool(YmassNames[6].data())==true){
            YPN_M300[i]->Fill(FatJet_ParticleNetMD_probXbb[leadingID[i][0]]);
            YM300[i]->Fill(FatJet_msoftdrop[leadingID[i][0]]);
            if(FatJet_ParticleNetMD_probXbb[leadingID[i][0]] > 0.8){
              YM300_p[i]->Fill(FatJet_msoftdrop[leadingID[i][0]]);  
            }
          }
          if(data.GetBool(YmassNames[7].data())==true){
            YPN_M400[i]->Fill(FatJet_ParticleNetMD_probXbb[leadingID[i][0]]);
          }
          if(data.GetBool(YmassNames[8].data())==true){
            YPN_M500[i]->Fill(FatJet_ParticleNetMD_probXbb[leadingID[i][0]]);
          }
          if(data.GetBool(YmassNames[9].data())==true){
            YPN_M600[i]->Fill(FatJet_ParticleNetMD_probXbb[leadingID[i][0]]);
          }
          if(data.GetBool(YmassNames[10].data())==true){
            YPN_M700[i]->Fill(FatJet_ParticleNetMD_probXbb[leadingID[i][0]]);
          }
          if(data.GetBool(YmassNames[11].data())==true){
            YPN_M800[i]->Fill(FatJet_ParticleNetMD_probXbb[leadingID[i][0]]);
          }
          if(data.GetBool(YmassNames[12].data())==true){
            YPN_M900[i]->Fill(FatJet_ParticleNetMD_probXbb[leadingID[i][0]]);
          }
          if(data.GetBool(YmassNames[13].data())==true){
            YPN_M1000[i]->Fill(FatJet_ParticleNetMD_probXbb[leadingID[i][0]]);
          }
          if(data.GetBool(YmassNames[14].data())==true){
            YPN_M1200[i]->Fill(FatJet_ParticleNetMD_probXbb[leadingID[i][0]]);
          }
          if(data.GetBool(YmassNames[15].data())==true){
            YPN_M1400[i]->Fill(FatJet_ParticleNetMD_probXbb[leadingID[i][0]]);
          }

        }//FOR OF Y DIFFERENT MASS RANGE   
      }//END OF PT
   }//END OF MASSPASS
 */  
  }//different eta value
      
 
    

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
      for(unsigned int iymass=0 ; iymass < nYmass; iymass++){
        if(data.GetBool(YmassNames[iymass].data())==true){
        npassym_P4[iymass]+=1;
        
      }
    }
  }

    


    variable->Fill();


    

    

  } // end of loop over events

  
 
  
  
  Float_t nytotal = 0 ;
  for(int i = 0 ; i < 16 ; i++){
    nytotal += npassym_P7[i];
  }
  cout << "nytotal = " << nytotal << endl;
  cout << "Number of total events = " << nTotal << endl;

  for(int i=0;i<20;i++)
    if(nPass[i]>0)cout << "nPass["<<i<<"]= " << nPass[i] << endl;

  
  

  





 

  
  // writing example output file
  TFile* outFile = new TFile(outputFileName.data(),"recreate");
//  gr->Write();
  
  
  heve->Write();
  m_H->Write();
  m_Y->Write();
  
  for(int i = 0 ; i < 24 ; i++){
  //  YPN_M90[i]->Write();
  //  YPN_M100[i]->Write();
  //  YPN_M125[i]->Write();
  //  YPN_M150[i]->Write();
    YPN_M200[i]->Write();
    YPN_M250[i]->Write();
  //  YPN_M300[i]->Write();
  //  YPN_M400[i]->Write();
  //  YPN_M500[i]->Write();
  //  YPN_M600[i]->Write();
  //  YPN_M700[i]->Write();
  //  YPN_M700[i]->Write();
  //  YPN_M800[i]->Write();
  //  YPN_M900[i]->Write();
  //  YPN_M1000[i]->Write();
  //  YPN_M1200[i]->Write();
  //  YPN_M1400[i]->Write();
    YM200[i]->Write();
    YM250[i]->Write();
    YM300[i]->Write();
    YM200_l[i]->Write();
    YM250_l[i]->Write();
    YM300_l[i]->Write();
    YM200_s[i]->Write();
    YM250_s[i]->Write();
    YM300_s[i]->Write();
    YM200_p[i]->Write();
    YM250_p[i]->Write();
    YM300_p[i]->Write();
  }
  DDE->Write();
  

  HT->Write();
//  YMass->Write();
  variable->Write();
  outFile->Close();
  
  
  
  
 } // end of macro
 




