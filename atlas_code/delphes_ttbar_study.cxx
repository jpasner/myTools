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
void delphes_ttbar_study() {
  //SetAtlasStyle();
  // Higgs sample
  //TFile *myFile = TFile::Open("/global/projecta/projectdirs/atlas/jpasner/storage/marcos_ntuples/evttree-mc16_13TeV.410471.PhPy8EG_A14_ttbar_hdamp258p75_allhad.deriv.DAOD_EXOT8.e6337_e5984_s3126_r10201_r10210_p3529.v4.root");
  TFile *myFile = TFile::Open("/global/homes/j/jpasner/2_evttree-ttbar-pythia8-NNPDF30-DelphesATLAS.root");
  //TFile output_file("output_file.root","RECREATE"); // Store all output in 1 file
  TTree *tree = (TTree*) myFile->Get("evttree");

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
  TVector3 fatJet3vector;
  TLorentzVector fatJet4vector;
  tree->SetBranchAddress("fatJetPt",&fatJetPt);
  tree->SetBranchAddress("fatJetPx",&fatJetPx);
  tree->SetBranchAddress("fatJetPy",&fatJetPy);
  tree->SetBranchAddress("fatJetPz",&fatJetPz);
  tree->SetBranchAddress("fatJetE",&fatJetE);
  tree->SetBranchAddress("fatJetM",&fatJetM);

  //trackJet
  float trackJetEta = 0;
  float trackJetPhi = 0;
  vector<float> *trackJetPt = 0;
  vector<float> *trackJetPx = 0;
  vector<float> *trackJetPy = 0;
  vector<float> *trackJetPz = 0;
  vector<float> *trackJetE = 0;
  vector<float> *trackJetM = 0;
  TVector3 trackJet3vector;
  TLorentzVector trackJet4vector;
  tree->SetBranchAddress("trackJetPt",&trackJetPt);
  tree->SetBranchAddress("trackJetPx",&trackJetPx);
  tree->SetBranchAddress("trackJetPy",&trackJetPy);
  tree->SetBranchAddress("trackJetPz",&trackJetPz);
  tree->SetBranchAddress("trackJetE",&trackJetE);
  tree->SetBranchAddress("trackJetM",&trackJetM);

  // Histograms!
  TH2D *containment_400 = new TH2D("containment_400","400 < Signal Fatjet pT < 600",2,0,2,3,0,3);
  TH2D *containment_600 = new TH2D("containment_600","600 < Signal Fatjet pT < 800",2,0,2,3,0,3);
  TH2D *containment_800 = new TH2D("containment_800","800 < Signal Fatjet pT",2,0,2,3,0,3);
  TH1F *h_contained_fatJetM = new TH1F("contained_fatJet_mass","Contained vs. Uncontained oppFatJetMass",30,0,300);
  TH1F *h_uncontained_fatJetM = new TH1F("uncontained_fatJet_mass","uncontained_fatJet_mass",30,0,300);
  TH1F *h_contained_fatJetPt = new TH1F("contained_fatJet_pT","Contained vs. Uncontained oppFatJetPt",30,250,1000);
  TH1F *h_uncontained_fatJetPt = new TH1F("uncontained_fatJet_pT","uncontained_fatJet_pT",30,250,1000);

  TH2F *dPhi_W_fatJet_0 = new TH2F("dPhi_W_fatJet_0","400 < Signal Fatjet pT < 600",100,0,1.5,100,0,1.5);
  TH2F *dPhi_W_fatJet_1 = new TH2F("dPhi_W_fatJet_1","600 < Signal Fatjet pT < 800",100,0,1.5,100,0,1.5);
  TH2F *dPhi_W_fatJet_2 = new TH2F("dPhi_W_fatJet_2","800 < Signal Fatjet pT",100,0,1.5,100,0,1.5);
  TH1F *h_nFatJets = new TH1F("nFatJets","Histogram of number of fatjets in event",7,0,7);
  TH1F *h_nBosons = new TH1F("nBosons","Histogram of number of W bosons in event",5,0,5);
  TH1F *h_dPhi_wBosons = new TH1F("dPhi_wBosons","angle between w bosons",100,0,7);
  TH1F *h_dPhi_bPartons = new TH1F("dPhi_bPartons","angle between b partons from tops",100,0,7);
  TH1F *h_dPhi_tPartons = new TH1F("dPhi_tPartons","angle between t partons",100,0,7);
  TH1F *h_wSpread = new TH1F("wSpread","Histogram of distance from opposite fatJet to all w partons",100,0,5);
  TH1F *h_wTrackSpread = new TH1F("wTrackSpread","Histogram of distance from opposite fatJet to w trackJets",100,0,2);
  TH1F *h_bTrackSpread = new TH1F("bTrackSpread","Histogram of distance from opposite fatJet to b trackJet",100,0,2);
  

  // Needed variables
  int bottom = 0;
  int anti_bottom = 0;
  int pass_preSelection = 0;
  int pass_statusCodes = 0;
  int counter = 0;
  int w_in_opp = 0;
  int b_in_opp = 0;
  int nBoson = 0;
  vector<TLorentzVector> bottom4vectors;
  vector<TLorentzVector> fatJet4vectors;
  vector<TLorentzVector> wBoson4vectors;
  vector<TLorentzVector> all_wParton4vectors;
  vector<TLorentzVector> opp_wParton4vectors;
  vector<TLorentzVector> opp_bParton4vectors;
  vector<int> parentVector;
  float tt_mass = 0;

  int w_parton_in_fatJet_400 = 0;
  int b_parton_in_fatJet_400 = 0;
  int w_parton_in_fatJet_600 = 0;
  int b_parton_in_fatJet_600 = 0;
  int w_parton_in_fatJet_800 = 0;
  int b_parton_in_fatJet_800 = 0;

  //********************************************************************************************
  //********************************************************************************************

  int nEvent = tree->GetEntries();
  cout << "***** Starting Event Loop *****" << endl;
  for(int event_itr = 0; event_itr < nEvent; event_itr++) {
    //if(event_itr > 500) {
    //  break;
    //}
    tree->GetEntry(event_itr);
    if(fatJetPt->at(0) > 400 && fatJetPt->at(1) > 250 && fatJetM->at(0) < 145.0 && fatJetM->at(0) > 105.0) {
      // Require signal fatjet pT > 400GeV and opposite fatjet pT > 250GeV.  Also require mass of signal fatjet is higgs-like
      pass_preSelection++;
      parentVector.clear();
      bottom4vectors.clear();
      fatJet4vectors.clear();
      wBoson4vectors.clear();
      opp_bParton4vectors.clear();
      all_wParton4vectors.clear();
      opp_wParton4vectors.clear();
      b_in_opp = 0;
      w_in_opp = 0;
      tt_mass = 0;
      w_parton_in_fatJet_400 = 0;
      b_parton_in_fatJet_400 = 0;
      w_parton_in_fatJet_600 = 0;
      b_parton_in_fatJet_600 = 0;
      w_parton_in_fatJet_800 = 0;

      h_nFatJets->Fill(fatJetPt->size());

      for(int fatJet_itr = 0; fatJet_itr < fatJetPt->size(); fatJet_itr++) {
        fatJet3vector.SetXYZ(fatJetPx->at(fatJet_itr),fatJetPy->at(fatJet_itr),fatJetPz->at(fatJet_itr));
        fatJet4vector.SetPxPyPzE(fatJetPx->at(fatJet_itr),fatJetPy->at(fatJet_itr),fatJetPz->at(fatJet_itr),fatJetE->at(fatJet_itr));
        fatJet4vectors.push_back(fatJet4vector);
      } // End of fatJet loop

      for(int parton_itr2 = 0; parton_itr2 < partonPt->size(); parton_itr2++) {
          if(partonStatus->at(parton_itr2) == 23) {
            if(partonPdgId->at(parton_itr2) == 5 && partonParentPdgId->at(parton_itr2) == 6) {
              bottom++;
            }
            else if(partonPdgId->at(parton_itr2) == -5 && partonParentPdgId->at(parton_itr2) == -6) {
              anti_bottom++;
            }
          } 
       }

      if(partonPt->size() > 1) {
        //if(abs(partonPdgId->at(0))==6&&abs(partonPdgId->at(0))==abs(partonPdgId->at(1))&&partonPdgId->at(0)!=partonPdgId->at(1)) {
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
                //cout << partonParentPdgId->at(parton_itr) << " " << parton4vector.DeltaPhi(fatJet4vectors[1]) << endl;
                opp_wParton4vectors.push_back(parton4vector);
                parentVector.push_back(partonParentPdgId->at(parton_itr));
              } // If parton within deltaR 1.0 of opposite fatJet
            } // If parent was a W
          } // End of parton loop
        //} // First two partons are top quarks with opposite sign
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
        // Require two partons from a W within DeltaPhi < 1.0 of the opposite fatJet
        continue;
      }
      if(parentVector[0] != parentVector[1]) {
        // The two wPartons came from different parents!
        continue;
      }
      pass_statusCodes++;

      nBoson = 0;
      for(int boson_itr = 0; boson_itr < bosonPt->size(); boson_itr++) {
        boson3vector.SetXYZ(bosonPx->at(boson_itr),bosonPy->at(boson_itr),bosonPz->at(boson_itr));
        boson4vector.SetPxPyPzE(bosonPx->at(boson_itr),bosonPy->at(boson_itr),bosonPz->at(boson_itr),bosonE->at(boson_itr));
        bosonPhi = boson4vector.Phi();
        bosonEta = boson4vector.Eta();

        if(abs(bosonPdgId->at(boson_itr)) == 24) {
          nBoson++;
          wBoson4vectors.push_back(boson4vector);
          if(abs(fatJet4vectors[1].DeltaPhi(boson4vector)) < 1.5) {
            if(fatJetPt->at(0)  > 400 && fatJetPt->at(0) < 600) {
              dPhi_W_fatJet_0->Fill(abs(boson4vector.DeltaPhi(fatJet4vectors[1])), abs(bottom4vectors[0].DeltaPhi(fatJet4vectors[1])));
            }
            if(fatJetPt->at(0)  > 600 && fatJetPt->at(0) < 800) {
              dPhi_W_fatJet_1->Fill(abs(boson4vector.DeltaPhi(fatJet4vectors[1])), abs(bottom4vectors[0].DeltaPhi(fatJet4vectors[1])));
            }
            else if(fatJetPt->at(0)  > 800) {
              dPhi_W_fatJet_2->Fill(abs(boson4vector.DeltaPhi(fatJet4vectors[1])), abs(bottom4vectors[0].DeltaPhi(fatJet4vectors[1])));
            }
          } // If W boson opposite of signal candidate
        } // If W boson
      } // End boson loop

      if(fatJetPt->at(0)  > 400 && fatJetPt->at(0) < 600) {
        if(abs(opp_wParton4vectors[0].DeltaR(fatJet4vectors[1])) < 1.0) w_parton_in_fatJet_400++;
        if(abs(opp_wParton4vectors[1].DeltaR(fatJet4vectors[1])) < 1.0) w_parton_in_fatJet_400++;
        if(abs(opp_bParton4vectors[0].DeltaR(fatJet4vectors[1])) < 1.0) b_parton_in_fatJet_400++;
        containment_400->Fill(b_parton_in_fatJet_400,w_parton_in_fatJet_400);
      }
      if(fatJetPt->at(0)  > 600 && fatJetPt->at(0) < 800) {
        if(abs(opp_wParton4vectors[0].DeltaR(fatJet4vectors[1])) < 1.0) w_parton_in_fatJet_600++;
        if(abs(opp_wParton4vectors[1].DeltaR(fatJet4vectors[1])) < 1.0) w_parton_in_fatJet_600++;
        if(abs(opp_bParton4vectors[0].DeltaR(fatJet4vectors[1])) < 1.0) b_parton_in_fatJet_600++;
        containment_600->Fill(b_parton_in_fatJet_600,w_parton_in_fatJet_600);
      }
      else if(fatJetPt->at(0)  > 800) {
        if(abs(opp_wParton4vectors[0].DeltaR(fatJet4vectors[1])) < 1.0) w_parton_in_fatJet_800++;
        if(abs(opp_wParton4vectors[1].DeltaR(fatJet4vectors[1])) < 1.0) w_parton_in_fatJet_800++;
        if(abs(opp_bParton4vectors[0].DeltaR(fatJet4vectors[1])) < 1.0) b_parton_in_fatJet_800++;
        containment_800->Fill(b_parton_in_fatJet_800,w_parton_in_fatJet_800);
      }

      if((b_parton_in_fatJet_400 + w_parton_in_fatJet_400 == 3) || (b_parton_in_fatJet_600 + w_parton_in_fatJet_600 == 3) || (b_parton_in_fatJet_800 + w_parton_in_fatJet_800 == 3)) {
        h_contained_fatJetPt->Fill(fatJet4vectors[1].Pt()); 
        h_contained_fatJetM->Fill(fatJet4vectors[1].M());
        counter++;
      }
      else{
        h_uncontained_fatJetPt->Fill(fatJet4vectors[1].Pt()); 
        h_uncontained_fatJetM->Fill(fatJet4vectors[1].M());
      }

      h_dPhi_bPartons->Fill(abs(bottom4vectors[0].DeltaPhi(bottom4vectors[1])));
      h_dPhi_wBosons->Fill(abs(wBoson4vectors[0].DeltaPhi(wBoson4vectors[1])));
      h_nBosons->Fill(wBoson4vectors.size());
    } // If Signal Candidate fatjet looks like higgs
  } // End of tree loop

  cout << "Pass Pre-Selections: " << pass_preSelection << endl;
  cout << "Pass StatusCode and Pdg: " << pass_statusCodes << endl;
  cout << "counter = " << counter << endl;

  cout << "nBottom = " << bottom << endl;
  cout << "nAntiBottom = " << anti_bottom << endl;

 TCanvas *c_containment_400 = new TCanvas("containment_400","containment_400",800,800);
  containment_400->GetXaxis()->SetTitle("number of b partons contained in oppFatJet");
  containment_400->GetYaxis()->SetTitle("number of W parton contained in oppFatJet");
  containment_400->SetStats(0);
  containment_400->Scale(1/(containment_400->Integral()));
  containment_400->Draw("COLZ TEXT"); // Draw with colors and bin value in text
  c_containment_400->SetRightMargin(0.18);
  c_containment_400->SetLeftMargin(0.18);
  c_containment_400->SaveAs("containment_400.pdf");

 TCanvas *c_containment_600 = new TCanvas("containment_600","containment_600",800,800);
  containment_600->GetXaxis()->SetTitle("number of b partons contained in oppFatJet");
  containment_600->GetYaxis()->SetTitle("number of W parton contained in oppFatJet");
  containment_600->SetStats(0);
  containment_600->Scale(1/(containment_600->Integral()));
  containment_600->Draw("COLZ TEXT"); // Draw with colors and bin value in text
  c_containment_600->SetRightMargin(0.18);
  c_containment_600->SetLeftMargin(0.18);
  c_containment_600->SaveAs("containment_600.pdf");

 TCanvas *c_containment_800 = new TCanvas("containment_800","containment_800",800,800);
  containment_800->GetXaxis()->SetTitle("number of b partons contained in oppFatJet");
  containment_800->GetYaxis()->SetTitle("number of W parton contained in oppFatJet");
  containment_800->SetStats(0);
  containment_800->Scale(1/(containment_800->Integral()));
  containment_800->Draw("COLZ TEXT"); // Draw with colors and bin value in text
  c_containment_800->SetRightMargin(0.18);
  c_containment_800->SetLeftMargin(0.18);
  c_containment_800->SaveAs("containment_800.pdf");

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

  TCanvas *c_dPhi_wBosons = new TCanvas("dPhi_wBosons","dPhi_wBosons",800,800);
  h_dPhi_wBosons->GetXaxis()->SetTitle("dPhi between W bosons");
  h_dPhi_wBosons->Draw();
  c_dPhi_wBosons->SetRightMargin(0.18);
  c_dPhi_wBosons->SetLeftMargin(0.18);
  c_dPhi_wBosons->SaveAs("dPhi_wBosons.pdf");

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
  h_wSpread->GetXaxis()->SetTitle("DeltaPhi of opp. fatJet and partons from a W");
  h_wSpread->Draw();
  c_wSpread->SetRightMargin(0.18);
  c_wSpread->SetLeftMargin(0.18);
  c_wSpread->SaveAs("wSpread.pdf");

  TCanvas *c_dPhi_W_fatJet_0 = new TCanvas("dPhi_W_fatJet_0","dPhi_W_fatJet",800,800);
  dPhi_W_fatJet_0->GetXaxis()->SetTitle("dPhi of W and opposite fatjet");
  dPhi_W_fatJet_0->GetYaxis()->SetTitle("dPhi of b parton and opposite fatJet");
  dPhi_W_fatJet_0->Draw("COLZ");
  c_dPhi_W_fatJet_0->SetRightMargin(0.18);
  c_dPhi_W_fatJet_0->SetLeftMargin(0.18);
  c_dPhi_W_fatJet_0->SaveAs("dPhi_W_fatJet_0.pdf");

  TCanvas *c_dPhi_W_fatJet_1 = new TCanvas("dPhi_W_fatJet_1","dPhi_W_fatJet",800,800);
  dPhi_W_fatJet_1->GetXaxis()->SetTitle("dPhi of W and opposite fatjet");
  dPhi_W_fatJet_1->GetYaxis()->SetTitle("dPhi of b parton and opposite fatJet");
  dPhi_W_fatJet_1->Draw("COLZ");
  c_dPhi_W_fatJet_1->SetRightMargin(0.18);
  c_dPhi_W_fatJet_1->SetLeftMargin(0.18);
  c_dPhi_W_fatJet_1->SaveAs("dPhi_W_fatJet_1.pdf");

  TCanvas *c_dPhi_W_fatJet_2 = new TCanvas("dPhi_W_fatJet_2","dPhi_W_fatJet",800,800);
  dPhi_W_fatJet_2->GetXaxis()->SetTitle("dPhi of W and opposite fatjet");
  dPhi_W_fatJet_2->GetYaxis()->SetTitle("dPhi of b parton and opposite fatJet");
  dPhi_W_fatJet_2->Draw("COLZ");
  c_dPhi_W_fatJet_2->SetRightMargin(0.18);
  c_dPhi_W_fatJet_2->SetLeftMargin(0.18);
  c_dPhi_W_fatJet_2->SaveAs("dPhi_W_fatJet_2.pdf");
  //output_file.Write();

}
