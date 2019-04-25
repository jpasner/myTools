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
  TH2F *dR_W_fatJet_0 = new TH2F("dR_W_fatJet_0","400 < Signal Fatjet pT < 600",100,0,1.5,100,0,1.5);
  TH2F *dR_W_fatJet_1 = new TH2F("dR_W_fatJet_1","600 < Signal Fatjet pT < 800",100,0,1.5,100,0,1.5);
  TH2F *dR_W_fatJet_2 = new TH2F("dR_W_fatJet_2","800 < Signal Fatjet pT",100,0,1.5,100,0,1.5);
  TH1F *h_nFatJets = new TH1F("nFatJets","Histogram of number of fatjets in event",7,0,7);
  TH1F *h_nBosons = new TH1F("nBosons","Histogram of number of W bosons in event",5,0,5);
  TH1F *h_dR_wBosons = new TH1F("dR_wBosons","angle between w bosons",100,0,7);
  TH1F *h_dR_bPartons = new TH1F("dR_bPartons","angle between b partons from tops",100,0,7);
  TH1F *h_wSpread = new TH1F("wSpread","Histogram of distance from opposite fatJet to all w partons",100,0,5);
  TH1F *h_wTrackSpread = new TH1F("wTrackSpread","Histogram of distance from opposite fatJet to w trackJets",100,0,2);
  TH1F *h_bTrackSpread = new TH1F("bTrackSpread","Histogram of distance from opposite fatJet to b trackJet",100,0,2);
  

  // Needed variables
  int counter = 0;
  int w_in_opp = 0;
  int b_in_opp = 0;
  int trackJet_found = 0;
  int nBoson = 0;
  vector<TLorentzVector> bottom4vectors;
  vector<TLorentzVector> fatJet4vectors;
  vector<TLorentzVector> wBoson4vectors;
  vector<TLorentzVector> all_wParton4vectors;
  vector<TLorentzVector> opp_wParton4vectors;
  vector<TLorentzVector> oppFatJet_wPartons_trackJet4vectors;
  vector<TLorentzVector> oppFatJet_bPartons_trackJet4vectors;
  vector<int> parentVector;
  float tt_mass = 0;

  //********************************************************************************************
  //********************************************************************************************

  int nEvent = tree->GetEntries();
  cout << "***** Starting Event Loop *****" << endl;
  for(int event_itr = 0; event_itr < nEvent; event_itr++) {
    //if(event_itr > 500) {
    //  break;
    //}
    tree->GetEntry(event_itr);
    if((HLT_ht1000_L1J100 || HLT_j420_a10_lcw_L1J100 || HLT_j420_a10r_L1J100 || HLT_j380 || HLT_4j100) &&fatJetPt->size() > 1) {
      if(fatJetPt->at(0) > 400 && fatJetPt->at(1) > 250 && fatJetM->at(0) < 145.0 && fatJetM->at(0) > 105.0) {
        // Require signal fatjet pT > 400GeV and opposite fatjet pT > 250GeV.  Also require mass of signal fatjet is higgs-like
        parentVector.clear();
        bottom4vectors.clear();
        fatJet4vectors.clear();
        wBoson4vectors.clear();
        all_wParton4vectors.clear();
        opp_wParton4vectors.clear();
        oppFatJet_bPartons_trackJet4vectors.clear();
        oppFatJet_wPartons_trackJet4vectors.clear();
        b_in_opp = 0;
        w_in_opp = 0;
        tt_mass = 0;
        h_nFatJets->Fill(fatJetPt->size());

        for(int fatJet_itr = 0; fatJet_itr < fatJetPt->size(); fatJet_itr++) {
          fatJet3vector.SetXYZ(fatJetPx->at(fatJet_itr),fatJetPy->at(fatJet_itr),fatJetPz->at(fatJet_itr));
          fatJet4vector.SetPxPyPzE(fatJetPx->at(fatJet_itr),fatJetPy->at(fatJet_itr),fatJetPz->at(fatJet_itr),fatJetE->at(fatJet_itr));
          fatJet4vectors.push_back(fatJet4vector);
        } // End of fatJet loop


        if(partonPt->size() > 1) {
          if(abs(partonPdgId->at(0))==6&&abs(partonPdgId->at(0))==abs(partonPdgId->at(1))&&partonPdgId->at(0)!=partonPdgId->at(1)) {
            //cout << "******" << endl;
            for(int parton_itr = 0; parton_itr < partonPt->size(); parton_itr++) {
              parton4vector.SetPxPyPzE(partonPx->at(parton_itr),partonPy->at(parton_itr),partonPz->at(parton_itr),partonE->at(parton_itr));
              if(abs(partonParentPdgId->at(parton_itr)) == 6 && partonStatus->at(parton_itr) == 23 && abs(partonPdgId->at(parton_itr)) == 5) {
                bottom4vectors.push_back(parton4vector);
              } // End of top loop if statement
              if(abs(partonParentPdgId->at(parton_itr)) == 24) {
                // Store all partons from a W
                all_wParton4vectors.push_back(parton4vector);
                if(parton4vector.DeltaR(fatJet4vectors[1]) < 1.0) {
                  //cout << partonParentPdgId->at(parton_itr) << " " << parton4vector.DeltaR(fatJet4vectors[1]) << endl;
                  opp_wParton4vectors.push_back(parton4vector);
                  parentVector.push_back(partonParentPdgId->at(parton_itr));
                } // If parton within deltaR 1.0 of opposite fatJet
              } // If parent was a W
            } // End of parton loop
          } // First two partons are top quarks with opposite sign
        } // More than 1 parton in event

        if(bottom4vectors.size() > 2) {
          cout << "WARNING: More than 3 bottom quarks satisfy top association requirement" << endl;
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

        counter++;

        for(int trackJet_itr = 0; trackJet_itr < trackJetPt->size(); trackJet_itr++) {
          trackJet4vector.SetPxPyPzE(trackJetPx->at(trackJet_itr),trackJetPy->at(trackJet_itr),trackJetPz->at(trackJet_itr),trackJetE->at(trackJet_itr));
          if(trackJet4vector.DeltaR(bottom4vectors[0]) < .1 && bottom4vectors[0].DeltaR(fatJet4vectors[1]) < 1.0) {
            oppFatJet_bPartons_trackJet4vectors.push_back(trackJet4vector);
            if(trackJetIdFatJet->at(trackJet_itr) == 1) {
              b_in_opp++;
            }
          }
          else if(trackJet4vector.DeltaR(bottom4vectors[1]) < .1 && bottom4vectors[1].DeltaR(fatJet4vectors[1]) < 1.0) {
            oppFatJet_bPartons_trackJet4vectors.push_back(trackJet4vector);
            if(trackJetIdFatJet->at(trackJet_itr) == 1) {
              b_in_opp++;
            }
          }
          else if(trackJet4vector.DeltaR(opp_wParton4vectors[0]) < .1) {
            oppFatJet_wPartons_trackJet4vectors.push_back(trackJet4vector);
            if(trackJetIdFatJet->at(trackJet_itr) == 1) {
              w_in_opp++;
            }
          }
          else if(trackJet4vector.DeltaR(opp_wParton4vectors[1]) < .1) {
            oppFatJet_wPartons_trackJet4vectors.push_back(trackJet4vector);
            if(trackJetIdFatJet->at(trackJet_itr) == 1) {
              w_in_opp++;
            }
          }
        }
        if(oppFatJet_wPartons_trackJet4vectors.size() < 2) {
          continue;
        }

        if(oppFatJet_wPartons_trackJet4vectors[0].E() == 0 || oppFatJet_wPartons_trackJet4vectors[1].E() == 0) {
          continue;
        }

        if(oppFatJet_bPartons_trackJet4vectors.size() != 1) {
          continue;
        }

        if(oppFatJet_bPartons_trackJet4vectors[0].E() == 0) {
          continue;
        }

        if(!(b_in_opp == 1 && w_in_opp == 2)) {
         cout << "B : " << b_in_opp << endl;
        cout << "W : " << w_in_opp << endl;
        cout << "Wooops" << endl;
        }

      

        trackJet_found++;

        nBoson = 0;
        for(int boson_itr = 0; boson_itr < bosonPt->size(); boson_itr++) {
          boson3vector.SetXYZ(bosonPx->at(boson_itr),bosonPy->at(boson_itr),bosonPz->at(boson_itr));
          boson4vector.SetPxPyPzE(bosonPx->at(boson_itr),bosonPy->at(boson_itr),bosonPz->at(boson_itr),bosonE->at(boson_itr));
          bosonPhi = boson4vector.Phi();
          bosonEta = boson4vector.Eta();

          if(abs(bosonPdgId->at(boson_itr)) == 24) {
            nBoson++;
            wBoson4vectors.push_back(boson4vector);
            if(fatJet4vectors[1].DeltaR(boson4vector) < 1.5) {
              if(fatJetPt->at(0)  > 400 && fatJetPt->at(0) < 600) {
                dR_W_fatJet_0->Fill(boson4vector.DeltaR(fatJet4vectors[1]), oppFatJet_bPartons_trackJet4vectors[0].DeltaR(fatJet4vectors[1]));
              }
              if(fatJetPt->at(0)  > 600 && fatJetPt->at(0) < 800) {
                dR_W_fatJet_1->Fill(boson4vector.DeltaR(fatJet4vectors[1]), oppFatJet_bPartons_trackJet4vectors[0].DeltaR(fatJet4vectors[1]));
              }
              else if(fatJetPt->at(0)  > 800) {
                dR_W_fatJet_2->Fill(boson4vector.DeltaR(fatJet4vectors[1]), oppFatJet_bPartons_trackJet4vectors[0].DeltaR(fatJet4vectors[1]));
              }


            } // If W boson opposite of signal candidate
          } // If W boson
        } // End boson loop

        for(int num = 0; num < all_wParton4vectors.size(); num++) {
          // Plot deltaR between opposite fatJet and all partons from a W boson
          h_wSpread->Fill(all_wParton4vectors[num].DeltaR(fatJet4vectors[1]));
        }
        h_wTrackSpread->Fill(oppFatJet_wPartons_trackJet4vectors[0].DeltaR(fatJet4vectors[1]));
        h_wTrackSpread->Fill(oppFatJet_wPartons_trackJet4vectors[1].DeltaR(fatJet4vectors[1]));
        h_bTrackSpread->Fill(oppFatJet_bPartons_trackJet4vectors[0].DeltaR(fatJet4vectors[1]));

        h_dR_bPartons->Fill(bottom4vectors[0].DeltaR(bottom4vectors[1]));
        h_dR_wBosons->Fill(wBoson4vectors[0].DeltaR(wBoson4vectors[1]));
        h_nBosons->Fill(wBoson4vectors.size());
      } // If Signal Candidate fatjet looks like higgs
    } // End of trigger if statement
  } // End of tree loop

  TCanvas *c_wTrackSpread = new TCanvas("wTrackSpread","wTrackSpread",800,800);
  h_wTrackSpread->GetXaxis()->SetTitle("dR opposite fatJet to w trackJets");
  h_wTrackSpread->Draw();
  c_wTrackSpread->SetRightMargin(0.18);
  c_wTrackSpread->SetLeftMargin(0.18);
  c_wTrackSpread->SaveAs("wTrackSpread.pdf");

  TCanvas *c_bTrackSpread = new TCanvas("bTrackSpread","bTrackSpread",800,800);
  h_bTrackSpread->GetXaxis()->SetTitle("dR opposite fatJet to b trackJet");
  h_bTrackSpread->Draw();
  c_bTrackSpread->SetRightMargin(0.18);
  c_bTrackSpread->SetLeftMargin(0.18);
  c_bTrackSpread->SaveAs("bTrackSpread.pdf");

  TCanvas *c_dR_wBosons = new TCanvas("dR_wBosons","dR_wBosons",800,800);
  h_dR_wBosons->GetXaxis()->SetTitle("dR between W bosons");
  h_dR_wBosons->Draw();
  c_dR_wBosons->SetRightMargin(0.18);
  c_dR_wBosons->SetLeftMargin(0.18);
  c_dR_wBosons->SaveAs("dR_wBosons.pdf");

  TCanvas *c_dR_bPartons = new TCanvas("dR_bPartons","dR_bPartons",800,800);
  h_dR_bPartons->GetXaxis()->SetTitle("dR between b partons from top quarks");
  h_dR_bPartons->Draw();
  c_dR_bPartons->SetRightMargin(0.18);
  c_dR_bPartons->SetLeftMargin(0.18);
  c_dR_bPartons->SaveAs("dR_bPartons.pdf");

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
  h_wSpread->GetXaxis()->SetTitle("Delta R of opp. fatJet and partons from a W");
  h_wSpread->Draw();
  c_wSpread->SetRightMargin(0.18);
  c_wSpread->SetLeftMargin(0.18);
  c_wSpread->SaveAs("wSpread.pdf");

  TCanvas *c_dR_W_fatJet_0 = new TCanvas("dR_W_fatJet_0","dR_W_fatJet",800,800);
  dR_W_fatJet_0->GetXaxis()->SetTitle("dR of W and opposite fatjet");
  dR_W_fatJet_0->GetYaxis()->SetTitle("dR of b trackJet and opposite fatJet");
  dR_W_fatJet_0->Draw("COLZ");
  c_dR_W_fatJet_0->SetRightMargin(0.18);
  c_dR_W_fatJet_0->SetLeftMargin(0.18);
  c_dR_W_fatJet_0->SaveAs("dR_W_fatJet_0.pdf");

  TCanvas *c_dR_W_fatJet_1 = new TCanvas("dR_W_fatJet_1","dR_W_fatJet",800,800);
  dR_W_fatJet_1->GetXaxis()->SetTitle("dR of W and opposite fatjet");
  dR_W_fatJet_1->GetYaxis()->SetTitle("dR of b trackJet and opposite fatJet");
  dR_W_fatJet_1->Draw("COLZ");
  c_dR_W_fatJet_1->SetRightMargin(0.18);
  c_dR_W_fatJet_1->SetLeftMargin(0.18);
  c_dR_W_fatJet_1->SaveAs("dR_W_fatJet_1.pdf");

  TCanvas *c_dR_W_fatJet_2 = new TCanvas("dR_W_fatJet_2","dR_W_fatJet",800,800);
  dR_W_fatJet_2->GetXaxis()->SetTitle("dR of W and opposite fatjet");
  dR_W_fatJet_2->GetYaxis()->SetTitle("dR of b trackJet and opposite fatJet");
  dR_W_fatJet_2->Draw("COLZ");
  c_dR_W_fatJet_2->SetRightMargin(0.18);
  c_dR_W_fatJet_2->SetLeftMargin(0.18);
  c_dR_W_fatJet_2->SaveAs("dR_W_fatJet_2.pdf");

  cout << "Events containing 2 b partons from \"good\" top quarks:" << counter << endl;
  cout << "Events containing 2 tops and we find a b-trackjet: " << trackJet_found << endl;

  //output_file.Write();
}
