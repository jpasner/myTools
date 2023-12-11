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
void resolved_jets_study() {
  //SetAtlasStyle();
  // Higgs sample
  TFile *myFile = TFile::Open("/global/projecta/projectdirs/atlas/jpasner/storage/marcos_ntuples/evttree-mc16_13TeV.410471.PhPy8EG_A14_ttbar_hdamp258p75_allhad.deriv.DAOD_EXOT8.e6337_e5984_s3126_r10201_r10210_p3529.v4.root");
  //TFile output_file("output_file.root","RECREATE"); // Store all output in 1 file
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
  vector<float> *partonPdgId = 0;
  vector<float> *partonBarcode = 0;
  vector<float> *partonParentPdgId = 0;
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
  tree->SetBranchAddress("partonParentPdgId",&partonParentPdgId);
  tree->SetBranchAddress("partonBarcode",&partonParentPdgId);
  tree->SetBranchAddress("partonPdgId",&partonPdgId);
  tree->SetBranchAddress("partonPt",&partonPt);
  tree->SetBranchAddress("partonPx",&partonPx);
  tree->SetBranchAddress("partonPy",&partonPy);
  tree->SetBranchAddress("partonPz",&partonPz);
  tree->SetBranchAddress("partonE",&partonE);
  tree->SetBranchAddress("partonM",&partonM);

  //Bosons
  int nBosons = 0;
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
  TVector3 fatJet3vector;
  TLorentzVector fatJet4vector;
  tree->SetBranchAddress("fatJetPt",&fatJetPt);
  tree->SetBranchAddress("fatJetPx",&fatJetPx);
  tree->SetBranchAddress("fatJetPy",&fatJetPy);
  tree->SetBranchAddress("fatJetPz",&fatJetPz);
  tree->SetBranchAddress("fatJetE",&fatJetE);
  tree->SetBranchAddress("fatJetM",&fatJetM);
  tree->SetBranchAddress("fatJetNGhostH",&fatJetNGhostH);

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

  // Histograms!
  TH2F *dR_hadrons_fatJet = new TH2F("dR_partons_fatJet","dR between b hadron jets and higgs associated fatJet",30,0,1.5,30,0,1.5);
  TH1F *dR_bb = new TH1F("dR_bb","dR between bb pair associated with higgs",30,0,5);
  TH1F *m_higgs = new TH1F("m_higgs","higgs mass",30,0,300);
  TH1F *pT_higgs = new TH1F("pT_higgs","higgs pT",100,0,1500);
  TH1F *pT_b1 = new TH1F("pT_b1","pT of higher pT b-hadron",20,0,700);
  TH1F *pT_b2 = new TH1F("pT_b2","pT of lower pT b-hadron",20,0,700);

  // Needed variables
  float bParton_energy_sum = 0;
  float bHadron_bParton_0_angle = 0;
  float bHadron_bParton_1_angle = 0;
  float b_parton_opening_angle = 0;
  float deltaE_bPartons = 0;
  float bHadron_fatJet_angle = 0;
  float angle = 0;
  float dR = 0;
  float dR_parton0_higgs_container = 0;
  float dR_parton1_higgs_container = 0;
  float dR_hadron0_fatJet_container = 0;
  float dR_hadron1_fatJet_container = 0;
  float deltaPt_b_partons = 0;
  int bHadron_paired_to_bParton0 = -1; // Negative 1 => parton hasn't been paired with a b-hadron
  int bHadron_paired_to_bParton1 = -1;
  int bParton_matched_fatJetId = -1; //fatJetId that the bParton in question has been matched to
  int HiggsAssociationFlag = 0; // Flag indicating the fatJet has been associated to a higgs (using specified mechanism)
  int SemiLepBhadDecay_flag = 0; // 1 -> includes semiLeptonic decay, 0 -> noLep decay
  int counter = 0;

  //Variables for b-hadrons falling outside of fatJet
  int n_bins = 5;
  int bin_iterator = 0; 
  float Higgs_pt_Bins[] = {300,365,407,449,512,800};
  vector<float> n_higgs_ptBins {0,0,0,0,0};
  vector<float> n_hadron_outside_ptBins {0,0,0,0,0};
  vector<float> n_oneMatch_ptBins {0,0,0,0,0};
  vector<float> n_noMatch_in_ptBin {0,0,0,0,0};

  //*******************************     Controlling flags     **********************************
  //********************************************************************************************

  int higgs_fatJet_association_flag = 1; // 1 -> My Association, 0 -> fatJetNGhostH association

  //********************************************************************************************
  //********************************************************************************************

  int nEvent = tree->GetEntries();
  cout << "***** Starting Event Loop *****" << endl;
  for(int event_itr = 0; event_itr < nEvent; event_itr++) {
    //if(event_itr > 50) {
    //  break;
    //}
    tree->GetEntry(event_itr);
    if( HLT_ht1000_L1J100 || HLT_j420_a10_lcw_L1J100 || HLT_j420_a10r_L1J100 || HLT_j380 || HLT_4j100) {

      for(int boson_itr = 0; boson_itr < bosonPt->size(); boson_itr++){
        // Reset all flags for higgs event
        SemiLepBhadDecay_flag = 0;
        bHad_NLepSec.clear();
        bHad_NLepTer.clear();
        bParton_matched_fatJetId = -1; //fatJetId that the bParton in question has been matched to
        bHadron_paired_to_bParton0 = -1; // Reset to negative 1 => parton hasn't been paired with a b-hadron
        bHadron_paired_to_bParton1 = -1;
        b_parton4vectors.clear();
        b_parton3vectors.clear();
        deltaE_bPartons = 0;
        b_parton_opening_angle = 0;

        if(bosonPdgId->at(boson_itr) == 25 && bosonStatus->at(boson_itr) == 62) {
          // Parton Loop
          for(int parton_itr = 0; parton_itr < partonPt->size(); parton_itr++) {
            parton3vector.SetXYZ(partonPx->at(parton_itr),partonPy->at(parton_itr),partonPz->at(parton_itr));
            parton4vector.SetPxPyPzE(partonPx->at(parton_itr),partonPy->at(parton_itr),partonPz->at(parton_itr),partonE->at(parton_itr));
            partonPhi = parton4vector.Phi();
            partonEta = parton4vector.Eta();

            if(abs(partonPdgId->at(parton_itr)) == 5 && partonParentPdgId->at(parton_itr) == 25) {
              // Parton is a b-quark and its parent was a higgs
              b_parton3vectors.push_back(parton3vector);
              b_parton4vectors.push_back(parton4vector);
            }
          } // End of Parton loop

          if(b_parton3vectors.size() > 2) {
            cout << "WARNING: More than 2 b partons have been associated with the higgs.  This is unexpected behavior in this plotting script!!!" << endl;
            cout << "***** Exiting *****" << endl;
            exit(EXIT_FAILURE);
          }

          // Opening angle between b-partons and smearing resolution study
          //if(b_parton3vectors.size() > 1 && bosonPt->at(boson_itr) > 450 && bosonPdgId->at(boson_itr) == 25) {
          if(b_parton3vectors.size() > 1) {
            dR_bb->Fill(b_parton4vectors.at(0).DeltaR(b_parton4vectors.at(1)));
            pT_higgs->Fill(bosonPt->at(boson_itr));
            if(b_parton4vectors.at(0).DeltaR(b_parton4vectors.at(1)) > 1.0) {
              counter++;
            }
            // Found two b-partons paired with the Higgs 
            b_parton_opening_angle = b_parton3vectors.at(0).Angle(b_parton3vectors.at(1));
            deltaPt_b_partons = abs(b_parton4vectors.at(0).Pt() - b_parton4vectors.at(1).Pt()) / bosonPt->at(boson_itr);
          } // End of if checking we have 2 b-partons paired to the Higgs

          // traackJetBhadron Loop
          for(int trackJetBHadron_itr = 0; trackJetBHadron_itr < trackJetBHadronPt->size(); trackJetBHadron_itr++) {
            trackJetBHadron3vector.SetXYZ(trackJetBHadronPx->at(trackJetBHadron_itr),trackJetBHadronPy->at(trackJetBHadron_itr),trackJetBHadronPz->at(trackJetBHadron_itr));
            trackJetBHadron4vector.SetPxPyPzE(trackJetBHadronPx->at(trackJetBHadron_itr),trackJetBHadronPy->at(trackJetBHadron_itr),trackJetBHadronPz->at(trackJetBHadron_itr),trackJetBHadronE->at(trackJetBHadron_itr));
            trackJetBHadronPhi = trackJetBHadron4vector.Phi();
            trackJetBHadronEta = trackJetBHadron4vector.Eta();

            // This block should assign each of the 2 higgs matched b-partons to one (and only one)
            // b-hadron in the event and then saves the iterator number, 3 vector and 4 vector for those hadrons 
            if(trackJetBHadronPt->at(trackJetBHadron_itr) > 0) { // Bad trackJetBHadronPt's are given negative values
              bHadron_bParton_0_angle = trackJetBHadron3vector.Angle(b_parton3vectors[0]);
              bHadron_bParton_1_angle = trackJetBHadron3vector.Angle(b_parton3vectors[1]);
              if(bHadron_bParton_0_angle < .1 || bHadron_bParton_1_angle < .1) { // Make sure the bHadron is within .01 angle of either of the b-partons
                if( bHadron_bParton_0_angle < bHadron_bParton_1_angle) { // Assign the hadron to the parton it is closest to
                  bHadron3vector_paired_to_bParton0 = trackJetBHadron3vector;
                  bHadron4vector_paired_to_bParton0 = trackJetBHadron4vector;
                  bHadron_paired_to_bParton0 = trackJetBHadron_itr;
                }
                else {
                  bHadron3vector_paired_to_bParton1 = trackJetBHadron3vector;
                  bHadron4vector_paired_to_bParton1 = trackJetBHadron4vector;
                  bHadron_paired_to_bParton1 = trackJetBHadron_itr;
                }
              }
            }
          } // End of trackJetBHadron loop

          // In this block you check that both partons have been paired with a b hadron
          if(bHadron_paired_to_bParton0 > -1 && bHadron_paired_to_bParton1 > -1) {
            if(trackJetIdFatJet->at(bHadron_paired_to_bParton0) == trackJetIdFatJet->at(bHadron_paired_to_bParton1)) {
              // Both partons tagged to the same fatJet and that fatJet is ghost associated to a Higgs
              bParton_matched_fatJetId = trackJetIdFatJet->at(bHadron_paired_to_bParton0);
            }
          }


          // FatJet Loop
          for(int fatJet_itr = 0; fatJet_itr < fatJetPt->size(); fatJet_itr++) {
            // Reject events where you don't get both b-trackjets (hadrons)
            if( b_parton3vectors.size() < 2) {
              continue;
            }
            // higgs_fatJet_association_flag == 0 means use the fatJetNGhostH association
            if(fatJetNGhostH->at(fatJet_itr) == 1 && higgs_fatJet_association_flag == 0) {
              HiggsAssociationFlag = 1;
            }
            // higgs_fatJet_association_flag == 1 means use my fatJet to Higgs association
            else if(fatJet_itr == bParton_matched_fatJetId && higgs_fatJet_association_flag == 1) {
              HiggsAssociationFlag = 1;
            }
            // No higgs-fatJet association
            else {
              HiggsAssociationFlag = 0;
            }

            if(fatJetPt->at(fatJet_itr) > 450 && HiggsAssociationFlag == 1) {
              fatJet3vector.SetXYZ(fatJetPx->at(fatJet_itr),fatJetPy->at(fatJet_itr),fatJetPz->at(fatJet_itr));
              fatJet4vector.SetPxPyPzE(fatJetPx->at(fatJet_itr),fatJetPy->at(fatJet_itr),fatJetPz->at(fatJet_itr),fatJetE->at(fatJet_itr));
              fatJetPhi = fatJet4vector.Phi();
              fatJetEta = fatJet4vector.Eta();

              //counter++;

              if(bHadron4vector_paired_to_bParton0.DeltaR(bHadron4vector_paired_to_bParton1) > 1.0) {
                //counter++;
              }

              // dR between fatJet and b Hadron jets
              dR_hadron0_fatJet_container = bHadron4vector_paired_to_bParton0.DeltaR(fatJet4vector);
              dR_hadron1_fatJet_container = bHadron4vector_paired_to_bParton1.DeltaR(fatJet4vector);

              // Angle between Fatjet and Higgs
              angle = fatJet3vector.Angle(boson3vector);
              // Mass of ghost associated fatJets

              // DeltaR for ghost associated fatJets and higgs
              dR = fatJet4vector.DeltaR(boson4vector);
              // Eta for ghost associated fatJet and higgs

              if(fatJetM->at(fatJet_itr) < 105 || fatJetM->at(fatJet_itr) > 145 ) {
                // Events with fatjet mass inside range specified on last line
                //counter++;
                m_higgs->Fill(fatJetM->at(fatJet_itr));
                //dR_bb->Fill(bHadron4vector_paired_to_bParton0.DeltaR(bHadron4vector_paired_to_bParton1));
                if(bHadron4vector_paired_to_bParton0.Pt() > bHadron4vector_paired_to_bParton1.Pt()) {
                  pT_b1->Fill(bHadron4vector_paired_to_bParton0.Pt());
                  pT_b2->Fill(bHadron4vector_paired_to_bParton1.Pt());
                  dR_hadrons_fatJet->Fill(dR_hadron0_fatJet_container,dR_hadron1_fatJet_container);
                }
                else{
                  pT_b1->Fill(bHadron4vector_paired_to_bParton1.Pt());
                  pT_b2->Fill(bHadron4vector_paired_to_bParton0.Pt());
                  dR_hadrons_fatJet->Fill(dR_hadron1_fatJet_container,dR_hadron0_fatJet_container);
                }
              }
            }
          } // End of Fatjet loop
        } // End of if statement checking for higgs PdgId and Status Code = 62
      } // End of boson loop
    } // End of trigger if statement
  } // End of tree loop

  TCanvas *c_dR_hadrons_fatJet = new TCanvas("dR_hadrons_fatJet","dR_hadrons_fatJet",800,800);
  dR_hadrons_fatJet->GetXaxis()->SetTitle("dR of higher Pt b hadron jet");
  dR_hadrons_fatJet->GetYaxis()->SetTitle("dR of lower Pt b hadron jet");
  dR_hadrons_fatJet->Draw("COLZ");
  c_dR_hadrons_fatJet->SetRightMargin(0.18);
  c_dR_hadrons_fatJet->SetLeftMargin(0.18);
  c_dR_hadrons_fatJet->SaveAs("dR_hadrons_fatJet.pdf");

  TCanvas *c_dR_bb = new TCanvas("dR_bb","dR_bb",800,800);
  dR_bb->GetXaxis()->SetTitle("dR of higgs b-partons");
  dR_bb->Draw();
  c_dR_bb->SetRightMargin(0.18);
  c_dR_bb->SetLeftMargin(0.18);
  c_dR_bb->SaveAs("dR_bb.pdf");

  TCanvas *c_m_higgs = new TCanvas("m_higgs","m_higgs",800,800);
  m_higgs->GetXaxis()->SetTitle("m_higgs");
  m_higgs->Draw();
  c_m_higgs->SetRightMargin(0.18);
  c_m_higgs->SetLeftMargin(0.18);
  c_m_higgs->SaveAs("m_higgs.pdf"); 

  TCanvas *c_pT_higgs = new TCanvas("pT_higgs","pT_higgs",800,800);
  pT_higgs->GetXaxis()->SetTitle("pT_higgs");
  pT_higgs->Draw();
  c_pT_higgs->SetRightMargin(0.18);
  c_pT_higgs->SetLeftMargin(0.18);
  c_pT_higgs->SaveAs("pT_higgs.pdf"); 


  TCanvas *c_pT_bHadrons = new TCanvas("pT_bHadrons","pT_bHadrons",800,800);
  TLegend *leg_pT_bHadrons = new TLegend(0.4,0.7,0.7,0.8);
  leg_pT_bHadrons->AddEntry(pT_b1, "Higher pT b hadron from Higgs");
  leg_pT_bHadrons->AddEntry(pT_b2, "Lower pT b hadron from Higgs");
  pT_b2->GetXaxis()->SetTitle("pT_bHadrons");
  pT_b1->SetLineColor(2);
  pT_b2->SetLineColor(4);
  pT_b2->Draw();
  pT_b1->Draw("same");
  leg_pT_bHadrons->Draw("same");
  c_pT_bHadrons->SaveAs("pT_bHadrons.pdf","recreate"); 

  //output_file.Write();

  cout << counter << endl;
}
