#include <iostream>
#include <vector>
#include <string>

#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TRandom.h>
#include <TVector.h>

// RooFit stuff
#include <RooBukinPdf.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooPlot.h>

// ATLAS Style!
#include "/global/homes/j/jpasner/atlasrootstyle/AtlasLabels.C"
#include "/global/homes/j/jpasner/atlasrootstyle/AtlasStyle.C"

using namespace std;

// ***************************************************************************************
// ***********                          Main Function                          ***********
// ***************************************************************************************
void ttbar_study() {
  //SetAtlasStyle();
  //Delphes pythia8 ttbar sample from Andrea
  //TFile *myFile = TFile::Open("/global/projecta/projectdirs/atlas/jpasner/andrea_Augmented_ntuple/user.asciandr.17912905.total.evttree.root");

  //Marco's Hbb samples
  TFile *myFile = TFile::Open("/afs/cern.ch/work/b/battagl/public/HbbISR/tuples/evttree-mc16_13TeV.309450.PowhegPy8EG_NNLOPS_nnlo_30_ggH125_bb_kt200.deriv.DAOD_EXOT8.e6281_e5984_s3126_r10201_r10210_p3529.v4.root");

  //TFile *myFile = TFile::Open("/global/projecta/projectdirs/atlas/jpasner/andrea_Augmented_ntuple/user.asciandr.mc16d_13TeV.410471.ntuple_AUGMENTED_evttree.root/user.asciandr.17912905._000001.evttree.root");
  //TFile *myFile = TFile::Open("/global/projecta/projectdirs/atlas/jpasner/storage/marcos_ntuples/evttree-mc16_13TeV.410471.PhPy8EG_A14_ttbar_hdamp258p75_allhad.deriv.DAOD_EXOT8.e6337_e5984_s3126_r10201_r10210_p3529.v4.root");
  TFile output_file("ggH_output_file.root","RECREATE"); // Store all output in 1 file
  TTree *tree = (TTree*) myFile->Get("evttree");

  //Triggers
  bool HLT_ht1000_L1J100 = 0;
  bool HLT_j420_a10_lcw_L1J100 = 0;
  bool HLT_j420_a10r_L1J100 = 0;
  bool HLT_j380 = 0;
  bool HLT_4j100 = 0;
  tree->SetBranchAddress("HLT_ht1000_L1J100", &HLT_ht1000_L1J100);
  tree->SetBranchAddress("HLT_j420_a10_lcw_L1J100", &HLT_j420_a10_lcw_L1J100);
  tree->SetBranchAddress("HLT_j420_a10r_L1J100", &HLT_j420_a10r_L1J100);
  tree->SetBranchAddress("HLT_j380", &HLT_j380);
  tree->SetBranchAddress("HLT_4j100", &HLT_4j100);

  // Partons
  int nPartons = 0;
  float partonEta = 0;
  float partonPhi = 0;
  vector<float> *partonStatus = 0;
  vector<float> *partonPdgId = 0;
  vector<float> *partonBarcode = 0;
  vector<float> *partonParentPdgId = 0;
  vector<float> *partonParentBarcode = 0;
  vector<float> *partonPt = 0;
  vector<float> *partonPx = 0;
  vector<float> *partonPy = 0;
  vector<float> *partonPz = 0;
  vector<float> *partonE = 0;
  vector<float> *partonM = 0;
  TVector3 parton3vector;
  vector<TVector3> b_parton3vectors;
  vector<TLorentzVector> b_parton4vectors;
  TLorentzVector parton4vector;
  tree->SetBranchAddress("partonStatus",&partonStatus);
  tree->SetBranchAddress("partonParentPdgId",&partonParentPdgId);
  tree->SetBranchAddress("partonParentBarcode",&partonParentBarcode);
  tree->SetBranchAddress("partonBarcode",&partonBarcode);
  tree->SetBranchAddress("partonPdgId",&partonPdgId);
  tree->SetBranchAddress("partonPt",&partonPt);
  tree->SetBranchAddress("partonPx",&partonPx);
  tree->SetBranchAddress("partonPy",&partonPy);
  tree->SetBranchAddress("partonPz",&partonPz);
  tree->SetBranchAddress("partonE",&partonE);
  tree->SetBranchAddress("partonM",&partonM);

  //Bosons
  float bosonEta = 0;
  float bosonPhi = 0;
  vector<float> *bosonStatus = 0;
  vector<float> *bosonPdgId = 0;
  vector<float> *bosonPt = 0;
  vector<float> *bosonPx = 0;
  vector<float> *bosonPy = 0;
  vector<float> *bosonPz = 0;
  vector<float> *bosonE = 0;
  vector<float> *bosonM = 0;
  TVector3 boson3vector;
  TLorentzVector boson4vector;
  tree->SetBranchAddress("bosonStatus",&bosonStatus);
  tree->SetBranchAddress("bosonPdgId",&bosonPdgId);
  tree->SetBranchAddress("bosonPt",&bosonPt);
  tree->SetBranchAddress("bosonPx",&bosonPx);
  tree->SetBranchAddress("bosonPy",&bosonPy);
  tree->SetBranchAddress("bosonPz",&bosonPz);
  tree->SetBranchAddress("bosonE",&bosonE);
  tree->SetBranchAddress("bosonM",&bosonM);

  //Fatjets
  int nFatJets = 0;
  float fatJetEta = 0;
  float fatJetPhi = 0;
  vector<float> *fatJetPdgId = 0;
  vector<float> *fatJetPt = 0;
  vector<float> *fatJetPx = 0;
  vector<float> *fatJetPy = 0;
  vector<float> *fatJetPz = 0;
  vector<float> *fatJetE = 0;
  vector<float> *fatJetM = 0;
  vector<float> *fatJetNGhostH = 0;
  vector<float> *fatJetNVRJets = 0;
  TVector3 fatJet3vector;
  TLorentzVector fatJet4vector;
  tree->SetBranchAddress("fatJetPt",&fatJetPt);
  tree->SetBranchAddress("fatJetPx",&fatJetPx);
  tree->SetBranchAddress("fatJetPy",&fatJetPy);
  tree->SetBranchAddress("fatJetPz",&fatJetPz);
  tree->SetBranchAddress("fatJetE",&fatJetE);
  tree->SetBranchAddress("fatJetM",&fatJetM);
  tree->SetBranchAddress("fatJetNGhostH",&fatJetNGhostH);
  tree->SetBranchAddress("fatJetNVRJets",&fatJetNVRJets);

  //trackJetBHadron
  int ntrackJetBHadron = 0;
  float trackJetBHadronEta = 0;
  float trackJetBHadronPhi = 0;
  vector<float> bHad_NLepSec;
  vector<float> bHad_NLepTer;
  vector<float> *trackJetBHadronNLepSec = 0;
  vector<float> *trackJetBHadronNLepTer = 0;
  vector<float> *trackJetBHadronPdgId = 0;
  vector<float> *trackJetBHadronPt = 0;
  vector<float> *trackJetBHadronPx = 0;
  vector<float> *trackJetBHadronPy = 0;
  vector<float> *trackJetBHadronPz = 0;
  vector<float> *trackJetBHadronE = 0;
  vector<float> *trackJetBHadronM = 0;
  TVector3 trackJetBHadron3vector;
  TVector3 bHadron3vector_paired_to_bParton0;
  TVector3 bHadron3vector_paired_to_bParton1;
  TLorentzVector bHadron4vector_paired_to_bParton0;
  TLorentzVector bHadron4vector_paired_to_bParton1;
  TLorentzVector trackJetBHadron4vector;
  tree->SetBranchAddress("trackJetBHadronNLepSec",&trackJetBHadronNLepSec);
  tree->SetBranchAddress("trackJetBHadronNLepTer",&trackJetBHadronNLepTer);
  tree->SetBranchAddress("trackJetBHadronPt",&trackJetBHadronPt);
  tree->SetBranchAddress("trackJetBHadronPx",&trackJetBHadronPx);
  tree->SetBranchAddress("trackJetBHadronPy",&trackJetBHadronPy);
  tree->SetBranchAddress("trackJetBHadronPz",&trackJetBHadronPz);
  tree->SetBranchAddress("trackJetBHadronE",&trackJetBHadronE);
  tree->SetBranchAddress("trackJetBHadronM",&trackJetBHadronM);

  //trackJet
  float trackJetEta = 0;
  float trackJetPhi = 0;
  vector<float> *trackJetIdFatJet = 0;
  vector<float> *trackJetPt = 0;
  vector<float> *trackJetPx = 0;
  vector<float> *trackJetPy = 0;
  vector<float> *trackJetPz = 0;
  vector<float> *trackJetE = 0;
  vector<float> *trackJetM = 0;
  TVector3 trackJet3vector;
  TLorentzVector trackJet4vector;
  tree->SetBranchAddress("trackJetIdFatJet",&trackJetIdFatJet);
  tree->SetBranchAddress("trackJetPt",&trackJetPt);
  tree->SetBranchAddress("trackJetPx",&trackJetPx);
  tree->SetBranchAddress("trackJetPy",&trackJetPy);
  tree->SetBranchAddress("trackJetPz",&trackJetPz);
  tree->SetBranchAddress("trackJetE",&trackJetE);
  tree->SetBranchAddress("trackJetM",&trackJetM);

  //vrJet
  float vrJetEta = 0;
  float vrJetPhi = 0;
  vector<float> *vrJetIdFatJet = 0;
  vector<float> *vrJetHybBEff_77 = 0;
  vector<float> *vrJetPt = 0;
  vector<float> *vrJetPx = 0;
  vector<float> *vrJetPy = 0;
  vector<float> *vrJetPz = 0;
  vector<float> *vrJetE = 0;
  vector<float> *vrJetM = 0;
  TVector3 vrJet3vector;
  TLorentzVector vrJet4vector;
  tree->SetBranchAddress("vrJetIdFatJet",&vrJetIdFatJet);
  tree->SetBranchAddress("vrJetHybBEff_77",&vrJetHybBEff_77);
  tree->SetBranchAddress("vrJetPt",&vrJetPt);
  tree->SetBranchAddress("vrJetPx",&vrJetPx);
  tree->SetBranchAddress("vrJetPy",&vrJetPy);
  tree->SetBranchAddress("vrJetPz",&vrJetPz);
  tree->SetBranchAddress("vrJetE",&vrJetE);
  tree->SetBranchAddress("vrJetM",&vrJetM);


  // Histograms!
  TH1F *h_dR_missingParton_oppFatJet_250 = new TH1F("dR_missingParton_oppFatJet_250","deltaR between oppFatJet and missing parton, 250 < pT < 450",90,0,3);
  TH1F *h_dR_missingParton_oppFatJet_450 = new TH1F("dR_missingParton_oppFatJet_450","deltaR between oppFatJet and missing parton, 450 < pT < 650",90,0,3);
  TH1F *h_dR_missingParton_oppFatJet_650 = new TH1F("dR_missingParton_oppFatJet_650","deltaR between oppFatJet and missing parton, 650 < pT",90,0,3);

  TH2D *containment_250 = new TH2D("containment_250","250 < Away Fatjet pT < 450",2,0,2,3,0,3);
  TH2D *containment_450 = new TH2D("containment_450","450 < Away Fatjet pT < 650",2,0,2,3,0,3);
  TH2D *containment_650 = new TH2D("containment_650","650 < Away Fatjet pT",2,0,2,3,0,3);

  TH1F *h_0B_2W_mass = new TH1F("contains_0B_2W","conatins 0B 2W opposite Fatjet mass",30,0,300);
  TH1F *h_0B_1W_mass = new TH1F("contains_0B_1W","conatins 0B 1W opposite Fatjet mass",30,0,300);
  TH1F *h_0B_0W_mass = new TH1F("contains_0B_0W","conatins 0B 0W opposite Fatjet mass",30,0,300);
  TH1F *h_1B_2W_mass = new TH1F("contains_1B_2W","conatins 1B 2W opposite Fatjet mass",30,0,300);
  TH1F *h_1B_1W_mass = new TH1F("contains_1B_1W","conatins 1B 1W opposite Fatjet mass",30,0,300);
  TH1F *h_1B_0W_mass = new TH1F("contains_1B_0W","conatins 1B 0W opposite Fatjet mass",30,0,300);
  TH1F* h_containment_mass_array[2][3];
  h_containment_mass_array[0][0] = h_0B_0W_mass;
  h_containment_mass_array[0][1] = h_0B_1W_mass;
  h_containment_mass_array[0][2] = h_0B_2W_mass;
  h_containment_mass_array[1][0] = h_1B_0W_mass;
  h_containment_mass_array[1][1] = h_1B_1W_mass;
  h_containment_mass_array[1][2] = h_1B_2W_mass;

  TH1F *h_contained_fatJetM = new TH1F("contained_fatJet_mass","Contained vs. Uncontained oppFatJetMass",60,0,300);
  TH1F *h_uncontained_fatJetM = new TH1F("uncontained_fatJet_mass","uncontained_fatJet_mass",60,0,300);
  TH1F *h_contained_fatJetPt = new TH1F("contained_fatJet_pT","Contained vs. Uncontained oppFatJetPt",60,250,1000);
  TH1F *h_uncontained_fatJetPt = new TH1F("uncontained_fatJet_pT","uncontained_fatJet_pT",60,250,1000);

  TH1F *h_nVRJets_oppFatJet = new TH1F("nVRJets_in_oppFatJet","Number of VRJets in oppFatJet that pass pT cut",7,0,7);
  TH1F *h_nVRJets_nearby_oppFatJet = new TH1F("nVRJets_nearby_oppFatJet","VRJets near oppFatJet that pass pT cut (1.0 < dR < 2.0)",7,0,7);

  TH1F *h_dR_VRjet_oppFatJet = new TH1F("VRjet_oppFatJet","dR between VRjets and oppFatJet",100,0,5);
  TH1F *h_dR_VRjet_pTCut_oppFatJet = new TH1F("VRjet_pTCut_oppFatJet","dR between VRjets and oppFatJet that pass pT cut",100,0,5);


  TH1F *h_3vr_fatJetM = new TH1F("3vr_fatJet_mass","Contained (3 vrJets) vs. Uncontained (2 vrJets) oppFatJetMass",60,0,300);
  TH1F *h_2vr_fatJetM = new TH1F("2vr_fatJet_mass","2vr_fatJet_mass",60,0,300);

  TH1F *h_3vr_chosen_fatJetM = new TH1F("3vr_chosen_fatJet_mass","Contained (3 vr_chosenJets) vs. Uncontained (2 vr_chosenJets) oppFatJetMass",60,0,300);
  TH1F *h_2vr_chosen_fatJetM = new TH1F("2vr_chosen_fatJet_mass","2vr_chosen_fatJet_mass",60,0,300);
  TH1F *h_2vr_corrected_fatJetM = new TH1F("2vr_corrected_fatJetM","2vr_corrected_fatJetM",60,0,300);

  TH2F *dR_W_fatJet_0 = new TH2F("dR_W_fatJet_0","250 < Away Fatjet pT < 450",100,0,1.5,100,0,1.5);
  TH2F *dR_W_fatJet_1 = new TH2F("dR_W_fatJet_1","450 < Away Fatjet pT < 650",100,0,1.5,100,0,1.5);
  TH2F *dR_W_fatJet_2 = new TH2F("dR_W_fatJet_2","650 < Away Fatjet pT",100,0,1.5,100,0,1.5);
  TH1F *h_nFatJets = new TH1F("nFatJets","Histogram of number of fatjets in event",7,0,7);
  TH1F *h_nBosons = new TH1F("nBosons","Histogram of number of W bosons in event",10,0,10);
  TH1F *h_dR_wBosons = new TH1F("dR_wBosons","angle between w bosons",100,0,7);
  TH1F *h_dPhi_bPartons = new TH1F("dPhi_bPartons","phi angle between b partons from tops",100,0,3.14);
  TH1F *h_wSpread = new TH1F("wSpread","Histogram of phi between opposite fatJet and all w partons",100,0,3.14);

  // Needed variables
  int pass_preSelection = 0;
  int pass_statusCodes = 0;
  int above_140 = 0;
  int counter = 0;
  int w_in_opp = 0;
  int b_in_opp = 0;
  TLorentzVector closest_vrJet;
  vector<TLorentzVector> escaped_parton4vectors;
  vector<TLorentzVector> bottom4vectors;
  vector<TLorentzVector> fatJet4vectors;
  vector<TLorentzVector> in_vrJet4vectors;
  vector<TLorentzVector> out_vrJet4vectors;
  vector<TLorentzVector> wBoson4vectors;
  vector<TLorentzVector> all_wParton4vectors;
  vector<TLorentzVector> opp_wParton4vectors;
  vector<TLorentzVector> opp_bParton4vectors;
  vector<int> parentVector;
  float tt_mass = 0;

  int nVRjets = 0;
  int nVR_B_jets = 0;

  int w_parton_in_fatJet_250 = 0;
  int b_parton_in_fatJet_250 = 0;
  int w_parton_in_fatJet_450 = 0;
  int b_parton_in_fatJet_450 = 0;
  int w_parton_in_fatJet_650 = 0;
  int b_parton_in_fatJet_650 = 0;
  int w_parton_in_fatjet = 0;
  int b_parton_in_fatjet = 0;

  //********************************************************************************************
  int truth_study = 0;
  //********************************************************************************************

  int nEvent = tree->GetEntries();
  cout << "***** Starting Event Loop *****" << endl;
  cout << "***** nEvent = " << nEvent << " *****" << endl;
  for(int event_itr = 0; event_itr < nEvent; event_itr++) {
    //if(event_itr > 10000) {
    //  break;
    //}
    tree->GetEntry(event_itr);
    if((HLT_ht1000_L1J100 || HLT_j420_a10_lcw_L1J100 || HLT_j420_a10r_L1J100 || HLT_j380 || HLT_4j100) &&fatJetPt->size() > 1) {
      if(fatJetPt->at(0) > 400 && fatJetPt->at(1) > 250 && fatJetM->at(0) < 145.0 && fatJetM->at(0) > 105.0) {
        // Require signal fatjet pT > 400GeV and opposite fatjet pT > 250GeV.  Also require mass of signal fatjet is higgs-like
        pass_preSelection++;
        escaped_parton4vectors.clear();
        parentVector.clear();
        bottom4vectors.clear();
        fatJet4vectors.clear();
        in_vrJet4vectors.clear();
        out_vrJet4vectors.clear();
        wBoson4vectors.clear();
        opp_bParton4vectors.clear();
        all_wParton4vectors.clear();
        opp_wParton4vectors.clear();
        b_in_opp = 0;
        w_in_opp = 0;
        tt_mass = 0;
        w_parton_in_fatJet_250 = 0;
        b_parton_in_fatJet_250 = 0;
        w_parton_in_fatJet_450 = 0;
        b_parton_in_fatJet_450 = 0;
        w_parton_in_fatJet_650 = 0;
        b_parton_in_fatJet_650 = 0;
        b_parton_in_fatjet = 0;
        w_parton_in_fatjet = 0;
        nVRjets = 0;
        nVR_B_jets = 0;
        closest_vrJet.SetPxPyPzE(0,0,0,0);

        h_nFatJets->Fill(fatJetPt->size());


        for(int fatJet_itr = 0; fatJet_itr < fatJetPt->size(); fatJet_itr++) {
          fatJet3vector.SetXYZ(fatJetPx->at(fatJet_itr),fatJetPy->at(fatJet_itr),fatJetPz->at(fatJet_itr));
          fatJet4vector.SetPxPyPzE(fatJetPx->at(fatJet_itr),fatJetPy->at(fatJet_itr),fatJetPz->at(fatJet_itr),fatJetE->at(fatJet_itr));
          fatJet4vectors.push_back(fatJet4vector);

          if(fatJet_itr == 1 && fatJetNVRJets->at(fatJet_itr) >= 3) {
            h_3vr_fatJetM->Fill(fatJet4vector.M());
          }
          else if(fatJet_itr == 1 && fatJetNVRJets->at(fatJet_itr) == 2) {
            h_2vr_fatJetM->Fill(fatJet4vector.M());
          }
        } // End of fatJet loop

        for(int vrJet_itr = 0; vrJet_itr < vrJetPt->size(); vrJet_itr++) {
          vrJet3vector.SetXYZ(vrJetPx->at(vrJet_itr),vrJetPy->at(vrJet_itr),vrJetPz->at(vrJet_itr));
          vrJet4vector.SetPxPyPzE(vrJetPx->at(vrJet_itr),vrJetPy->at(vrJet_itr),vrJetPz->at(vrJet_itr),vrJetE->at(vrJet_itr));

          h_dR_VRjet_oppFatJet->Fill(vrJet4vector.DeltaR(fatJet4vectors[1]));
          if(vrJet4vector.Pt() > 5) {
            h_dR_VRjet_pTCut_oppFatJet->Fill(vrJet4vector.DeltaR(fatJet4vectors[1]));
          }

          if(vrJet4vector.Pt() > 5 && vrJet4vector.DeltaR(fatJet4vectors[1]) < 1.0) {
            in_vrJet4vectors.push_back(vrJet4vector);
            nVRjets++;
          }
          else if(vrJet4vector.Pt() > 5 && vrJet4vector.DeltaR(fatJet4vectors[1]) > 1.0 && vrJet4vector.DeltaR(fatJet4vectors[1]) < 2.0) {
            out_vrJet4vectors.push_back(vrJet4vector);
          }
          if(vrJet4vector.Pt() > 5 && vrJetHybBEff_77->at(vrJet_itr) == 1 && vrJet4vector.DeltaR(fatJet4vectors[1]) < 1.0) {
            nVR_B_jets++;
          }
        } // End of vrJet loop

        h_nVRJets_oppFatJet->Fill(in_vrJet4vectors.size());
        if(nVRjets == 2) {
          h_nVRJets_nearby_oppFatJet->Fill(out_vrJet4vectors.size());
        }

        if(nVRjets == 2 && out_vrJet4vectors.size() > 0) {
          closest_vrJet = out_vrJet4vectors[0];
          for(int n = 0; n < out_vrJet4vectors.size(); n++) {
            if(closest_vrJet.DeltaR(fatJet4vectors[1]) > out_vrJet4vectors[n].DeltaR(fatJet4vectors[1])) {
              closest_vrJet = out_vrJet4vectors[n]; 
            }
          }
          h_2vr_corrected_fatJetM->Fill((fatJet4vectors[1] + closest_vrJet).M());
        }

        if(nVRjets >= 3 && nVR_B_jets > 0) {
          h_3vr_chosen_fatJetM->Fill(fatJet4vectors[1].M());
        }
        if(nVRjets == 2) {
          h_2vr_chosen_fatJetM->Fill(fatJet4vectors[1].M());
        }

        /*

        for(int boson_itr = 0; boson_itr < bosonPt->size(); boson_itr++) {
          boson3vector.SetXYZ(bosonPx->at(boson_itr),bosonPy->at(boson_itr),bosonPz->at(boson_itr));
          boson4vector.SetPxPyPzE(bosonPx->at(boson_itr),bosonPy->at(boson_itr),bosonPz->at(boson_itr),bosonE->at(boson_itr));
          bosonPhi = boson4vector.Phi();
          bosonEta = boson4vector.Eta();
          if(abs(bosonPdgId->at(boson_itr)) == 24 && bosonStatus->at(boson_itr) == 22) {
            wBoson4vectors.push_back(boson4vector);
          }
        }

        if(partonPt->size() > 1) {
          for(int parton_itr = 0; parton_itr < partonPt->size(); parton_itr++) {
            parton4vector.SetPxPyPzE(partonPx->at(parton_itr),partonPy->at(parton_itr),partonPz->at(parton_itr),partonE->at(parton_itr));
            if(abs(partonParentPdgId->at(parton_itr)) == 6 && partonStatus->at(parton_itr) == 23 && abs(partonPdgId->at(parton_itr)) == 5) {
              bottom4vectors.push_back(parton4vector);
              if(abs(parton4vector.DeltaPhi(fatJet4vectors[1])) < (M_PI/2.0)) {
                opp_bParton4vectors.push_back(parton4vector);
              }
            } // End of top loop if statement
            if(abs(partonParentPdgId->at(parton_itr)) == 24) {
              // Store all partons from a W
              all_wParton4vectors.push_back(parton4vector);
              if(abs(parton4vector.DeltaPhi(fatJet4vectors[1])) < (M_PI/2.0)) {
                opp_wParton4vectors.push_back(parton4vector);
                parentVector.push_back(partonParentPdgId->at(parton_itr));
              } // If parton within deltaR 1.0 of opposite fatJet
            } // If parent was a W
          } // End of parton loop
        } // More than 1 parton in event

        for(int num = 0; num < all_wParton4vectors.size(); num++) {
          // Plot deltaR between opposite fatJet and all partons from a W boson
          h_wSpread->Fill(abs(all_wParton4vectors[num].DeltaPhi(fatJet4vectors[1])));
        }

        if(bottom4vectors.size() > 2) {
          cout << "WARNING: More than 3 bottom quarks satisfy top association requirement" << endl;
        }
        if(opp_bParton4vectors.size() != 1) {
          continue;
        }
        if(bottom4vectors.size() < 2) {
          continue;
        }
        if(opp_wParton4vectors.size() != 2) {
          // Require two partons from a W within DeltaR < 1.0 of the opposite fatJet
          continue;
        }
        if(parentVector[0] != parentVector[1]) {
          // The two wPartons came from different parents!
          continue;
        }

        pass_statusCodes++;

        for(int boson_itr = 0; boson_itr < bosonPt->size(); boson_itr++) {
          boson3vector.SetXYZ(bosonPx->at(boson_itr),bosonPy->at(boson_itr),bosonPz->at(boson_itr));
          boson4vector.SetPxPyPzE(bosonPx->at(boson_itr),bosonPy->at(boson_itr),bosonPz->at(boson_itr),bosonE->at(boson_itr));
          bosonPhi = boson4vector.Phi();
          bosonEta = boson4vector.Eta();
          if(abs(bosonPdgId->at(boson_itr)) == 24 && bosonStatus->at(boson_itr) == 22) {
            if(abs(fatJet4vectors[1].DeltaPhi(boson4vector)) < 1.5) {
              if(fatJetPt->at(1)  > 250 && fatJetPt->at(1) < 450) {
                dR_W_fatJet_0->Fill(abs(boson4vector.DeltaR(fatJet4vectors[1])), abs(opp_bParton4vectors[0].DeltaR(fatJet4vectors[1])));
              }
              if(fatJetPt->at(1)  > 450 && fatJetPt->at(1) < 650) {
                dR_W_fatJet_1->Fill(abs(boson4vector.DeltaR(fatJet4vectors[1])), abs(opp_bParton4vectors[0].DeltaR(fatJet4vectors[1])));
              }
              else if(fatJetPt->at(1)  > 650) {
                dR_W_fatJet_2->Fill(abs(boson4vector.DeltaR(fatJet4vectors[1])), abs(opp_bParton4vectors[0].DeltaR(fatJet4vectors[1])));
              }
            } // If W boson opposite of signal candidate
          } // If W boson
        } // End boson loop

        if(fatJetPt->at(1)  > 250 && fatJetPt->at(1) < 450) {
          if(abs(opp_wParton4vectors[0].DeltaR(fatJet4vectors[1])) < 1.0) {
            w_parton_in_fatJet_250++;
          }
          else escaped_parton4vectors.push_back(opp_wParton4vectors[0]);
          if(abs(opp_wParton4vectors[1].DeltaR(fatJet4vectors[1])) < 1.0) {
            w_parton_in_fatJet_250++;
          }
          else escaped_parton4vectors.push_back(opp_wParton4vectors[1]);
          if(abs(opp_bParton4vectors[0].DeltaR(fatJet4vectors[1])) < 1.0) {
            b_parton_in_fatJet_250++;
          }
          else escaped_parton4vectors.push_back(opp_bParton4vectors[0]);
          containment_250->Fill(b_parton_in_fatJet_250,w_parton_in_fatJet_250);
        }
        else if(fatJetPt->at(1)  > 450 && fatJetPt->at(1) < 650) {
          if(abs(opp_wParton4vectors[0].DeltaR(fatJet4vectors[1])) < 1.0) {
            w_parton_in_fatJet_450++;
          }
          else escaped_parton4vectors.push_back(opp_wParton4vectors[0]);
          if(abs(opp_wParton4vectors[1].DeltaR(fatJet4vectors[1])) < 1.0) {
            w_parton_in_fatJet_450++;
          }
          else escaped_parton4vectors.push_back(opp_wParton4vectors[1]);
          if(abs(opp_bParton4vectors[0].DeltaR(fatJet4vectors[1])) < 1.0) { 
            b_parton_in_fatJet_450++;
          }
          else escaped_parton4vectors.push_back(opp_bParton4vectors[0]);
          containment_450->Fill(b_parton_in_fatJet_450,w_parton_in_fatJet_450);
        }
        else if(fatJetPt->at(1)  > 650) {
          if(abs(opp_wParton4vectors[0].DeltaR(fatJet4vectors[1])) < 1.0) {
            w_parton_in_fatJet_650++;
          }
          else escaped_parton4vectors.push_back(opp_wParton4vectors[0]);
          if(abs(opp_wParton4vectors[1].DeltaR(fatJet4vectors[1])) < 1.0) {
            w_parton_in_fatJet_650++;
          }
          else escaped_parton4vectors.push_back(opp_wParton4vectors[1]);
          if(abs(opp_bParton4vectors[0].DeltaR(fatJet4vectors[1])) < 1.0) {
            b_parton_in_fatJet_650++;
          }
          else escaped_parton4vectors.push_back(opp_bParton4vectors[0]);
          containment_650->Fill(b_parton_in_fatJet_650,w_parton_in_fatJet_650);
        }

        if((b_parton_in_fatJet_250 + w_parton_in_fatJet_250) == 2) {
          h_dR_missingParton_oppFatJet_250->Fill(fatJet4vectors[1].DeltaR(escaped_parton4vectors[0]));
        }
        else if((b_parton_in_fatJet_450 + w_parton_in_fatJet_450) == 2) {
          h_dR_missingParton_oppFatJet_450->Fill(fatJet4vectors[1].DeltaR(escaped_parton4vectors[0]));
        }
        else if((b_parton_in_fatJet_650 + w_parton_in_fatJet_650) == 2) {
          h_dR_missingParton_oppFatJet_650->Fill(fatJet4vectors[1].DeltaR(escaped_parton4vectors[0]));
        }

        b_parton_in_fatjet = b_parton_in_fatJet_250 + b_parton_in_fatJet_450 + b_parton_in_fatJet_650;
        w_parton_in_fatjet = w_parton_in_fatJet_250 + w_parton_in_fatJet_450 + w_parton_in_fatJet_650;

        h_containment_mass_array[b_parton_in_fatjet][w_parton_in_fatjet]->Fill(fatJet4vectors[1].M());

        if((b_parton_in_fatjet + w_parton_in_fatjet) == 3 ) {
          h_contained_fatJetPt->Fill(fatJet4vectors[1].Pt());
          h_contained_fatJetM->Fill(fatJet4vectors[1].M());
          counter++;
        }
        else{
          h_uncontained_fatJetPt->Fill(fatJet4vectors[1].Pt());
          h_uncontained_fatJetM->Fill(fatJet4vectors[1].M());
        }
      
        h_dPhi_bPartons->Fill(abs(bottom4vectors[0].DeltaPhi(bottom4vectors[1])));
        h_dR_wBosons->Fill(wBoson4vectors[0].DeltaR(wBoson4vectors[1]));
        h_nBosons->Fill(wBoson4vectors.size());
        */
      } // If Signal Candidate fatjet looks like higgs
    } // End of trigger if statement
  } // End of tree loop

  TCanvas *c_3vr_fatJetM = new TCanvas("3vr_fatJetM","3vr_fatJetM",800,800);
  TLegend *leg_3vr_fatJetM = new TLegend(0.6,0.6,0.9,0.7);
  leg_3vr_fatJetM->AddEntry(h_3vr_fatJetM,"Contains 3 vr jets");
  leg_3vr_fatJetM->AddEntry(h_2vr_fatJetM,"Conatins 2 vr jets");
  h_3vr_fatJetM->GetXaxis()->SetTitle("opposite fatJet Mass [GeV]");
  h_3vr_fatJetM->SetLineColor(2);
  h_2vr_fatJetM->SetLineColor(4);
  h_3vr_fatJetM->Draw();
  h_2vr_fatJetM->Draw("same");
  leg_3vr_fatJetM->Draw("same");
  c_3vr_fatJetM->SetRightMargin(0.18);
  c_3vr_fatJetM->SetLeftMargin(0.18);
  c_3vr_fatJetM->SaveAs("vrContainment_fatJetM.pdf");

  TCanvas *c_3vr_chosen_fatJetM = new TCanvas("3vr_chosen_fatJetM","3vr_chosen_fatJetM",800,800);
  TLegend *leg_3vr_chosen_fatJetM = new TLegend(0.6,0.6,0.9,0.7);
  leg_3vr_chosen_fatJetM->AddEntry(h_3vr_chosen_fatJetM,"3 deltaR associated vr jets && b-tag > 0");
  leg_3vr_chosen_fatJetM->AddEntry(h_2vr_chosen_fatJetM,"2 deltaR associated vr jets");
  //leg_3vr_chosen_fatJetM->AddEntry(h_2vr_corrected_fatJetM,"Corrected");
  h_3vr_chosen_fatJetM->GetXaxis()->SetTitle("opposite fatJet Mass [GeV]");
  h_3vr_chosen_fatJetM->SetLineColor(2);
  h_2vr_chosen_fatJetM->SetLineColor(4);
  //h_2vr_corrected_fatJetM->SetLineColor(6);
  h_3vr_chosen_fatJetM->Draw();
  h_2vr_chosen_fatJetM->Draw("same");
  //h_2vr_corrected_fatJetM->Draw("same");
  leg_3vr_chosen_fatJetM->Draw("same");
  c_3vr_chosen_fatJetM->SetRightMargin(0.18);
  c_3vr_chosen_fatJetM->SetLeftMargin(0.18);
  c_3vr_chosen_fatJetM->SaveAs("vr_chosen_Containment_fatJetM.pdf");

  TCanvas *c_dR_VRjet_oppFatJet = new TCanvas("dR_VRjet_oppFatJet","dR_VRjet_oppFatJet",800,800);
  TLegend *leg_dR_VRjet_oppFatJet = new TLegend(0.6,0.6,0.9,0.7);
  leg_dR_VRjet_oppFatJet->AddEntry(h_dR_VRjet_oppFatJet, "dR_VRjet_oppFatJet");
  leg_dR_VRjet_oppFatJet->AddEntry(h_dR_VRjet_pTCut_oppFatJet, "dR_VRjet_pTCut_oppFatJet");
  h_dR_VRjet_oppFatJet->GetXaxis()->SetTitle("dR VRJet and oppFatJet");
  h_dR_VRjet_oppFatJet->SetLineColor(2);
  h_dR_VRjet_pTCut_oppFatJet->SetLineColor(4);
  h_dR_VRjet_oppFatJet->Draw();
  h_dR_VRjet_pTCut_oppFatJet->Draw("same");
  leg_dR_VRjet_oppFatJet->Draw("same");
  c_dR_VRjet_oppFatJet->SetRightMargin(0.18);
  c_dR_VRjet_oppFatJet->SetLeftMargin(0.18);
  c_dR_VRjet_oppFatJet->SaveAs("dR_VRjet_oppFatJet.pdf");

  TCanvas *c_nVRJets_oppFatJet = new TCanvas("nVRJets_oppFatJet","nVRJets_oppFatJet",800,800);
  h_nVRJets_oppFatJet->GetXaxis()->SetTitle("Number of VRJets");
  h_nVRJets_oppFatJet->Scale(1.0/h_nVRJets_oppFatJet->Integral());
  h_nVRJets_oppFatJet->Draw();
  c_nVRJets_oppFatJet->SetRightMargin(0.18);
  c_nVRJets_oppFatJet->SetLeftMargin(0.18);
  c_nVRJets_oppFatJet->SaveAs("nVRJets_oppFatJet.pdf");

  TCanvas *c_nVRJets_nearby_oppFatJet = new TCanvas("nVRJets_nearby_oppFatJet","nVRJets_nearby_oppFatJet",800,800);
  h_nVRJets_nearby_oppFatJet->GetXaxis()->SetTitle("Number of VRJets near oppFatJet");
  h_nVRJets_nearby_oppFatJet->Scale(1.0/h_nVRJets_nearby_oppFatJet->Integral());
  h_nVRJets_nearby_oppFatJet->Draw();
  c_nVRJets_nearby_oppFatJet->SetRightMargin(0.18);
  c_nVRJets_nearby_oppFatJet->SetLeftMargin(0.18);
  c_nVRJets_nearby_oppFatJet->SaveAs("nVRJets_nearby_oppFatJet.pdf");

  cout << "Pass Pre-Selections: " << pass_preSelection << endl;

  TAxis *axis = h_3vr_chosen_fatJetM->GetXaxis();
  int xmin = 140;
  int xmax = 300;
  int bmin = axis->FindBin(xmin); 
  int bmax = axis->FindBin(xmax); 
  double integral = h_3vr_chosen_fatJetM->Integral(bmin,bmax);
  integral -= h_3vr_chosen_fatJetM->GetBinContent(bmin)*(xmin-axis->GetBinLowEdge(bmin))/
              axis->GetBinWidth(bmin);
  integral -= h_3vr_chosen_fatJetM->GetBinContent(bmax)*(axis->GetBinUpEdge(bmax)-xmax)/
              axis->GetBinWidth(bmax);
  cout << "n-events above m = 140GeV: " << integral << endl;

  /*
  cout << "Pass StatusCode and Pdg: " << pass_statusCodes << endl;
  cout << "counter = " << counter << endl;

  TCanvas *c_dR_missingParton_oppFatJet_250 = new TCanvas("dR_missingParton_oppFatJet_250","dR_missingParton_oppFatJet_250",800,800);
  h_dR_missingParton_oppFatJet_250->GetXaxis()->SetTitle("dR oppFatJet and missing parton");
  h_dR_missingParton_oppFatJet_250->Draw();
  c_dR_missingParton_oppFatJet_250->SetRightMargin(0.18);
  c_dR_missingParton_oppFatJet_250->SetLeftMargin(0.18);
  c_dR_missingParton_oppFatJet_250->SaveAs("dR_missingParton_oppFatJet_250.pdf");

  TCanvas *c_dR_missingParton_oppFatJet_450 = new TCanvas("dR_missingParton_oppFatJet_450","dR_missingParton_oppFatJet_450",800,800);
  h_dR_missingParton_oppFatJet_450->GetXaxis()->SetTitle("dR oppFatJet and missing parton");
  h_dR_missingParton_oppFatJet_450->Draw();
  c_dR_missingParton_oppFatJet_450->SetRightMargin(0.18);
  c_dR_missingParton_oppFatJet_450->SetLeftMargin(0.18);
  c_dR_missingParton_oppFatJet_450->SaveAs("dR_missingParton_oppFatJet_450.pdf");

  TCanvas *c_dR_missingParton_oppFatJet_650 = new TCanvas("dR_missingParton_oppFatJet_650","dR_missingParton_oppFatJet_650",800,800);
  h_dR_missingParton_oppFatJet_650->GetXaxis()->SetTitle("dR oppFatJet and missing parton");
  h_dR_missingParton_oppFatJet_650->Draw();
  c_dR_missingParton_oppFatJet_650->SetRightMargin(0.18);
  c_dR_missingParton_oppFatJet_650->SetLeftMargin(0.18);
  c_dR_missingParton_oppFatJet_650->SaveAs("dR_missingParton_oppFatJet_650.pdf");

  TCanvas *c_containment_mass_array = new TCanvas("containment","Mass of opposite faJet missing (x=B,y=W) partons",800,800);
  c_containment_mass_array->Divide(2,3);
  c_containment_mass_array->cd(1);
  h_containment_mass_array[0][2]->Draw();
  c_containment_mass_array->cd(2);
  h_containment_mass_array[1][2]->Draw();
  c_containment_mass_array->cd(3);
  h_containment_mass_array[0][1]->Draw();
  c_containment_mass_array->cd(4);
  h_containment_mass_array[1][1]->Draw();
  c_containment_mass_array->cd(5);
  h_containment_mass_array[0][0]->Draw();
  c_containment_mass_array->cd(6);
  h_containment_mass_array[1][0]->Draw();
  c_containment_mass_array->SaveAs("containment_mass_array.pdf");

  TCanvas *c_containment_250 = new TCanvas("containment_250","containment_250",800,800);
  containment_250->GetXaxis()->SetTitle("number of b partons contained in oppFatJet");
  containment_250->GetYaxis()->SetTitle("number of W parton contained in oppFatJet");
  containment_250->SetStats(0);
  containment_250->Scale(1/(containment_250->Integral()));
  containment_250->Draw("COLZ TEXT"); // Draw with colors and bin value in text
  c_containment_250->SetRightMargin(0.18);
  c_containment_250->SetLeftMargin(0.18);
  c_containment_250->SaveAs("containment_250.pdf");

  TCanvas *c_containment_450 = new TCanvas("containment_450","containment_450",800,800);
  containment_450->GetXaxis()->SetTitle("number of b partons contained in oppFatJet");
  containment_450->GetYaxis()->SetTitle("number of W parton contained in oppFatJet");
  containment_450->SetStats(0);
  containment_450->Scale(1/(containment_450->Integral()));
  containment_450->Draw("COLZ TEXT"); // Draw with colors and bin value in text
  c_containment_450->SetRightMargin(0.18);
  c_containment_450->SetLeftMargin(0.18);
  c_containment_450->SaveAs("containment_450.pdf");

 TCanvas *c_containment_650 = new TCanvas("containment_650","containment_650",800,800);
  containment_650->GetXaxis()->SetTitle("number of b partons contained in oppFatJet");
  containment_650->GetYaxis()->SetTitle("number of W parton contained in oppFatJet");
  containment_650->SetStats(0);
  containment_650->Scale(1/(containment_650->Integral()));
  containment_650->Draw("COLZ TEXT"); // Draw with colors and bin value in text
  c_containment_650->SetRightMargin(0.18);
  c_containment_650->SetLeftMargin(0.18);
  c_containment_650->SaveAs("containment_650.pdf");

  TCanvas *c_contained_fatJetM = new TCanvas("contained_fatJetM","contained_fatJetM",800,800);
  TLegend *leg_contained_fatJetM = new TLegend(0.2,0.7,0.4,0.8);
  leg_contained_fatJetM->AddEntry(h_contained_fatJetM,"All partons contained");
  leg_contained_fatJetM->AddEntry(h_uncontained_fatJetM,"Some partons not contained");
  h_contained_fatJetM->GetXaxis()->SetTitle("opposite fatJet Mass [GeV]");
  h_contained_fatJetM->SetLineColor(2);
  h_uncontained_fatJetM->SetLineColor(4);
  h_contained_fatJetM->Draw();
  h_uncontained_fatJetM->Draw("same");
  leg_contained_fatJetM->Draw("same");
  c_contained_fatJetM->SetRightMargin(0.18);
  c_contained_fatJetM->SetLeftMargin(0.18);
  c_contained_fatJetM->SaveAs("containment_fatJetM.pdf");

  TCanvas *c_contained_fatJetPt = new TCanvas("contained_fatJetPt","contained_fatJetPt",800,800);
  TLegend *leg_contained_fatJetPt = new TLegend(0.4,0.7,0.7,0.8);
  leg_contained_fatJetPt->AddEntry(h_contained_fatJetPt,"All partons contained");
  leg_contained_fatJetPt->AddEntry(h_uncontained_fatJetPt,"Some partons not contained");
  h_contained_fatJetPt->GetXaxis()->SetTitle("opposite fatJet pT [GeV]");
  h_contained_fatJetPt->SetLineColor(2);
  h_uncontained_fatJetPt->SetLineColor(4);
  h_contained_fatJetPt->Draw();
  h_uncontained_fatJetPt->Draw("same");
  leg_contained_fatJetPt->Draw("same");
  c_contained_fatJetPt->SetRightMargin(0.18);
  c_contained_fatJetPt->SetLeftMargin(0.18);
  c_contained_fatJetPt->SaveAs("containment_fatJetPt.pdf");

  TCanvas *c_dR_wBosons = new TCanvas("dR_wBosons","dR_wBosons",800,800);
  h_dR_wBosons->GetXaxis()->SetTitle("dR between W bosons");
  h_dR_wBosons->Draw();
  c_dR_wBosons->SetRightMargin(0.18);
  c_dR_wBosons->SetLeftMargin(0.18);
  c_dR_wBosons->SaveAs("dR_wBosons.pdf");

  TCanvas *c_dPhi_bPartons = new TCanvas("dPhi_bPartons","dPhi_bPartons",800,800);
  h_dPhi_bPartons->GetXaxis()->SetTitle("dPhi between b partons from top quarks");
  h_dPhi_bPartons->Draw();
  c_dPhi_bPartons->SetRightMargin(0.18);
  c_dPhi_bPartons->SetLeftMargin(0.18);
  c_dPhi_bPartons->SaveAs("dPhi_bPartons.pdf");

  TCanvas *c_nBosons = new TCanvas("nBosons","nBosons",800,800);
  h_nBosons->GetXaxis()->SetTitle("Number of W bosons in event");
  h_nBosons->Draw();
  c_nBosons->SetRightMargin(0.18);
  c_nBosons->SetLeftMargin(0.18);
  c_nBosons->SaveAs("nBosons.pdf");

  TCanvas *c_nFatJets = new TCanvas("nFatJets","nFatJets",800,800);
  h_nFatJets->GetXaxis()->SetTitle("Number of fatjets in event");
  h_nFatJets->Draw();
  c_nFatJets->SetRightMargin(0.18);
  c_nFatJets->SetLeftMargin(0.18);
  c_nFatJets->SaveAs("nFatJets.pdf");

  TCanvas *c_wSpread = new TCanvas("wSpread","wSpread",800,800);
  h_wSpread->GetXaxis()->SetTitle("Delta Phi of opp. fatJet and partons from a W");
  h_wSpread->Draw();
  c_wSpread->SetRightMargin(0.18);
  c_wSpread->SetLeftMargin(0.18);
  c_wSpread->SaveAs("wSpread.pdf");

  TCanvas *c_dR_W_fatJet_0 = new TCanvas("dR_W_fatJet_0","dR_W_fatJet",800,800);
  dR_W_fatJet_0->GetXaxis()->SetTitle("dR of W boson and opposite fatjet");
  dR_W_fatJet_0->GetYaxis()->SetTitle("dR of b parton and opposite fatJet");
  dR_W_fatJet_0->Draw("COLZ");
  c_dR_W_fatJet_0->SetRightMargin(0.18);
  c_dR_W_fatJet_0->SetLeftMargin(0.18);
  c_dR_W_fatJet_0->SaveAs("dR_W_fatJet_0.pdf");

  TCanvas *c_dR_W_fatJet_1 = new TCanvas("dR_W_fatJet_1","dR_W_fatJet",800,800);
  dR_W_fatJet_1->GetXaxis()->SetTitle("dR of W boson and opposite fatjet");
  dR_W_fatJet_1->GetYaxis()->SetTitle("dR of b parton and opposite fatJet");
  dR_W_fatJet_1->Draw("COLZ");
  c_dR_W_fatJet_1->SetRightMargin(0.18);
  c_dR_W_fatJet_1->SetLeftMargin(0.18);
  c_dR_W_fatJet_1->SaveAs("dR_W_fatJet_1.pdf");

  TCanvas *c_dR_W_fatJet_2 = new TCanvas("dR_W_fatJet_2","dR_W_fatJet",800,800);
  dR_W_fatJet_2->GetXaxis()->SetTitle("dR of W boson and opposite fatjet");
  dR_W_fatJet_2->GetYaxis()->SetTitle("dR of b parton and opposite fatJet");
  dR_W_fatJet_2->Draw("COLZ");
  c_dR_W_fatJet_2->SetRightMargin(0.18);
  c_dR_W_fatJet_2->SetLeftMargin(0.18);
  c_dR_W_fatJet_2->SaveAs("dR_W_fatJet_2.pdf");
  */

  output_file.Write();
}
