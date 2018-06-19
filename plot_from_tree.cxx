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
#include "/Users/Tibbins/work/atlasstyle/AtlasLabels.C"
#include "/Users/Tibbins/work/atlasstyle/AtlasStyle.C"

using namespace std;

void smear_E(TLorentzVector *v, double A) {
  // Smear the energy of this TlorentzVector by the value A, makeing sure to keep the mass invariant
  // To keep the mass invariant you must scale the momentum by the calcualted value B (see function)
  double E = v->E();
  double P = v->P();
  double Px = v->Px();
  double Py = v->Py();
  double Pz = v->Pz();
  double B = 0;

  // Add 1 to A -> E * (1 + A) = Smeared_E  
  A = 1 + A;
  // Scale factor (B) for ijk components of momentum such that mass doesn't change when E is scaled by A
  B = sqrt((E*E*(A*A - 1) + P*P) / (P*P)); 

  v->SetPxPyPzE(B*Px, B*Py, B*Pz, A*E);

  return;
}

void plot_TH1F(TH1F *hist, string output_filename, string title, string x_axis, string y_axis, string options, int canvas_width, int canvas_height){
  // Generic plotting function for TH1F plots
  TCanvas *can = new TCanvas(output_filename.c_str(),output_filename.c_str(),canvas_width,canvas_height);
  hist->SetTitle(title.c_str());
  hist->GetXaxis()->SetTitle(x_axis.c_str());
  hist->GetYaxis()->SetTitle(y_axis.c_str());
  hist->SetStats(true);
  //gStyle->SetOptStat(1);
  hist->Draw(options.c_str());
  can->SaveAs(output_filename.c_str());
  return;
}

void make_bin_plots(vector<TH1F*> hist_vector, string output_filename) {
  // This function is used for plotting the mass, pT and Energy bins of the fatJet Mass resolution process
  TCanvas *can = new TCanvas(output_filename.c_str(),output_filename.c_str(),800,600);
  can->Divide(3,2);
  for(int i = 0; i < hist_vector.size(); i++) {
    can->cd(i+1);// The +1 is because the TCanvas' are labeled starting with 1 instead of 0
    hist_vector[i]->Draw();
  }
  can->SaveAs(output_filename.c_str());
  return;
}

double gaus_width_of_FWHM(TH1F *hist) {
  // Fit FWHM area of peak with gaus and return width
  int peak_bin = hist->GetMaximumBin();
  double peak_value = hist->GetBinContent(peak_bin);
  double peak_location = hist->GetBinCenter(peak_bin);
  int low_FWHM_xbin = 0;
  hist->GetBinWithContent(peak_value / 2.0, low_FWHM_xbin, 0, peak_bin, 5000); // Fill low_FWHM_xbin with bin number (below peak_bin) that contains y = (y_peak / 2)
  int high_FWHM_xbin = 0;
  hist->GetBinWithContent(peak_value / 2.0, high_FWHM_xbin, peak_bin, hist->GetNbinsX(), 5000); // Fill high_FWHM_xbin with bin number (above peak_bin) that contains y = (y_peak / 2)
  double xLow = hist->GetBinCenter(low_FWHM_xbin); // Set low x boundary for fit to the lower boundry of the FWHM
  double xHigh = hist->GetBinCenter(high_FWHM_xbin); // Set high x boundary for fit to the high boundry of the FWHM

  hist->SetLineWidth(1);
  hist->Fit("gaus","","",xLow,xHigh); // Gaussian fit inside of FWHM region
  hist->GetFunction("gaus")->SetLineColor(4);
  hist->GetFunction("gaus")->SetLineWidth(1);
  return hist->GetFunction("gaus")->GetParameter(2); // Retrun width result from fit
}

double bukin_fit_width(TH1F *input_hist) {
  // Fit Smeared mass distribution using a bukin and return its width
  RooRealVar smeared_mass("smeared_mass","smeared_mass",input_hist->GetXaxis()->GetXmin() ,input_hist->GetXaxis()->GetXmax());

  // Declare histogram to be fit along and a variable (smeared_Mass) to be filled with fit results
  RooDataHist hist("hist","histogram of smeared Mass",smeared_mass,input_hist);

  // Declare bukin variables and function. Don't forget to give variables ranges or they will be considered constant in the fit
  RooRealVar peak_position("peak_position","location of peak on x axis",125,100,150);
  RooRealVar sigma("sigma","resolution of Mass defined as: FWHM divided by 2*sqrt(2*log(2))=2.35",20,0,60);
  RooRealVar peak_asymmetry("peak_asymmetry","peak asymmetry parameter",0,-0.5,0.5);
  RooRealVar left_tail("left_tail","parameter of the left tail",0,-0.5,0.5);
  RooRealVar right_tail("right_tail","parameter of the right tail",0,-0.5,0.5);
  RooBukinPdf bukin("my bukin name","my bukin title",smeared_mass,peak_position,sigma,peak_asymmetry,left_tail,right_tail);

  // Create frame for RooFit plot result to be plotted on
  RooPlot* frame = smeared_mass.frame();
  hist.plotOn(frame,RooFit::MarkerSize(0.1),RooFit::XErrorSize(0));

  // Fit the PDF to the data and plot on the frame we just created
  bukin.fitTo(hist);
  bukin.plotOn(frame,RooFit::LineWidth(1));
  frame->SetTitle(input_hist->GetTitle());
  frame->Draw();

  // Store results from fit
  return(sigma.getVal()); // return mass resolution result from fit
}

double bukin_fit_mean(TH1F *input_hist) {
  // Fit Smeared mass distribution using a bukin and return its width
  RooRealVar smeared_mass("smeared_mass","smeared_mass",50 ,200);

  // Declare histogram to be fit along and a variable (smeared_Mass) to be filled with fit results
  RooDataHist hist("hist","histogram of smeared Mass",smeared_mass,input_hist);

  // Declare bukin variables and function. Don't forget to give variables ranges or they will be considered constant in the fit
  RooRealVar peak_position("peak_position","location of peak on x axis",125,100,150);
  RooRealVar sigma("sigma","resolution of Mass defined as: FWHM divided by 2*sqrt(2*log(2))=2.35",20,0,60);
  RooRealVar peak_asymmetry("peak_asymmetry","peak asymmetry parameter",0,-0.5,0.5);
  RooRealVar left_tail("left_tail","parameter of the left tail",0,-0.5,0.5);
  RooRealVar right_tail("right_tail","parameter of the right tail",0,-0.5,0.5);
  RooBukinPdf bukin("my bukin name","my bukin title",smeared_mass,peak_position,sigma,peak_asymmetry,left_tail,right_tail);

  // Create frame for RooFit plot result to be plotted on
  RooPlot* frame = smeared_mass.frame();
  hist.plotOn(frame,RooFit::MarkerSize(0.1),RooFit::XErrorSize(0));

  // Fit the PDF to the data and plot on the frame we just created
  bukin.fitTo(hist);
  bukin.plotOn(frame,RooFit::LineWidth(1));
  frame->SetTitle(input_hist->GetTitle());
  TCanvas *can = new TCanvas("fit_canvas","fit_canvas",800,600);
  frame->Draw();
  string fit_output_filename = input_hist->GetName();
  fit_output_filename.append("_bukin_fit.pdf");

  TLegend *leg_bukin_fit_width = new TLegend(0.7,0.9,0.95,0.95);
  leg_bukin_fit_width->SetHeader(Form("Bukin width: (%f)", sigma.getVal()),"C");
  leg_bukin_fit_width->Draw();

  can->SaveAs(fit_output_filename.c_str());

  // Store results from fit
  return(peak_position.getVal()); // return mass resolution result from fit
}

void fit_histo(vector<TH1F*> hist_vector, vector<TH1F*> bins_hist_vector, vector<double> &var_mean_vector, vector<double> &mass_mean_vector, vector<double> &mass_sigma_vector, vector<double> &mass_sigma_error_vector, vector<double> &mass_error_low_vector, vector<double> &mass_error_high_vector, vector<double> &bin_low_edge_vector, vector<double> &bin_high_edge_vector, double xLow, double xHigh, string output_filename) {
  // Specialty function for extracting mass resolution of higgs associated fatJet in bins of E or pT of the Higgs or fatJet.
  // Implement and use bukin pdf for fit

  TCanvas *can = new TCanvas(output_filename.c_str(),output_filename.c_str(),800,600);
  can->Divide(3,2);

  for( int i = 0; i < hist_vector.size(); i++) {
    can->cd(i+1);// The +1 is because the TCanvas' are labeled starting with 1 instead of 0
    // Observable name and its range to be fit
    RooRealVar fatJet_mass("fatJet_mass","fatJet_mass",xLow,xHigh);

    // Declare histogram to be fit along and a variable (fatJet_Mass) to be filled with fit results
    RooDataHist hist("hist","histogram of fatJet Mass",fatJet_mass,hist_vector[i]); 

    // Declare bukin variables and function. Don't forget to give variables ranges or they will be considered constant in the fit
    RooRealVar peak_position("peak_position","location of peak on x axis",100,100,150);
    RooRealVar sigma("sigma","resolution of fatJet Mass defined as: FWHM divided by 2*sqrt(2*log(2))=2.35",10,0,50);
    RooRealVar peak_asymmetry("peak_asymmetry","peak asymmetry parameter",0,-1,1);
    RooRealVar left_tail("left_tail","parameter of the left tail",0,-1,1);
    RooRealVar right_tail("right_tail","parameter of the right tail",0,-1,1);
    RooBukinPdf bukin("my bukin name","my bukin title",fatJet_mass,peak_position,sigma,peak_asymmetry,left_tail,right_tail);

    // Create frame for RooFit plot result to be plotted on
    RooPlot* frame = fatJet_mass.frame();
    hist.plotOn(frame);

    // Fit the PDF to the data and plot on the frame we just created
    bukin.fitTo(hist);
    bukin.plotOn(frame);
    frame->SetTitle(hist_vector[i]->GetTitle()); 
    frame->Draw();

    // Store results from fit
    mass_mean_vector.push_back(peak_position.getVal());  // Store mean value from fit, should be ~125GeV for higgs
    mass_sigma_vector.push_back(sigma.getVal()); // Store mass resolution result from fit
    mass_sigma_error_vector.push_back(sigma.getError()); // Store mass resolution error from fit
    var_mean_vector.push_back(bins_hist_vector[i]->GetMean()); // Store mean of the E or pT for this E or pT bin

    // This is just to get the x-axis error bars for the fatJet mass resolution plot
    // Such that they show the choosen bin edges in pT or E of the higgs or fatJet
    mass_error_low_vector.push_back(bins_hist_vector[i]->GetMean() - bin_low_edge_vector[i]);
    mass_error_high_vector.push_back(bin_high_edge_vector[i] - bins_hist_vector[i]->GetMean());
  }
  can->SaveAs(output_filename.c_str());
  return;
}

// ***************************************************************************************
// ***********                          Main Function                          ***********
// ***************************************************************************************
void plot_from_tree() {
  SetAtlasStyle();
  // Higgs sample
  TFile *myFile = TFile::Open("/Users/Tibbins/work/hbbISR/samples/old/rel_21_signal/evttree-mc16_13TeV.345342.PowhegPy8EG_NNLOPS_nnlo_30_ggH125_bb.deriv.DAOD_EXOT8.e5739_e5984_s3126_r9364_r9315_p3309.root");
  //TFile *myFile = TFile::Open("/Users/Tibbins/work/hbbISR/samples/february_19_2018_samples/evttree-mc16_13TeV.309450.PowhegPy8EG_NNLOPS_nnlo_30_ggH125_bb_kt200.deriv.DAOD_EXOT8.e6281_e5984_s3126_r9364_r9315_p3309.v3.root");
  //TFile *myFile = TFile::Open("/Users/Tibbins/work/hbbISR/samples/february_19_2018_samples/evttree-mc16_13TeV.304708.Sherpa_CT10_Zqq_Pt500_1000.deriv.DAOD_EXOT8.e4692_s3126_r9364_r9315_p3309.v3.root");
  TFile output_file("output_file.root","RECREATE"); // Store all output in 1 file
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

  // ***** Start of Resolution Plot variable declaration ****

  // Mass Resolution Histogram in fatJet Pt bins
  vector<int> fatJet_pT_binning_vector = {20,20,20,20,20};
  vector<double> fatJet_pT_bin_low_edge = {300,365,407,449,512};
  vector<double> fatJet_pT_bin_high_edge = {365,407,449,512,800};
  vector<double> fatJet_pT_mean;
  vector<double> fatJet_pT_error_low;
  vector<double> fatJet_pT_error_high;

  vector<double> fatJet_bin_fatJet_mass_mean_pT;
  vector<double> fatJet_bin_fatJet_mass_sigma_pT;
  vector<double> fatJet_bin_fatJet_mass_sigma_error_pT;

  TH1F *fatJet_bin_M0 =  new TH1F("fatJet_bin_M0","Fatjet mass for fatJet pT Range 0",fatJet_pT_binning_vector[0],65,190);
  TH1F *fatJet_bin_M1 =  new TH1F("fatJet_bin_M1","Fatjet mass for fatJet pT Range 1",fatJet_pT_binning_vector[1],65,190);
  TH1F *fatJet_bin_M2 =  new TH1F("fatJet_bin_M2","Fatjet mass for fatJet pT Range 2",fatJet_pT_binning_vector[2],65,190);
  TH1F *fatJet_bin_M3 =  new TH1F("fatJet_bin_M3","Fatjet mass for fatJet pT Range 3",fatJet_pT_binning_vector[3],65,190);
  TH1F *fatJet_bin_M4 =  new TH1F("fatJet_bin_M4","Fatjet mass for fatJet pT Range 4",fatJet_pT_binning_vector[4],65,190);
  vector<TH1F*> fatJet_pT_mass_res_hist_vector = {fatJet_bin_M0,fatJet_bin_M1,fatJet_bin_M2,fatJet_bin_M3,fatJet_bin_M4};
  TH1F *fatJet_bin_FJ0 =  new TH1F("fatJet_bin_FJ0","fatJet pT Range 0",fatJet_pT_binning_vector[0],fatJet_pT_bin_low_edge[0],fatJet_pT_bin_high_edge[0]);
  TH1F *fatJet_bin_FJ1 =  new TH1F("fatJet_bin_FJ1","fatJet pT Range 1",fatJet_pT_binning_vector[1],fatJet_pT_bin_low_edge[1],fatJet_pT_bin_high_edge[1]);
  TH1F *fatJet_bin_FJ2 =  new TH1F("fatJet_bin_FJ2","fatJet pT Range 2",fatJet_pT_binning_vector[2],fatJet_pT_bin_low_edge[2],fatJet_pT_bin_high_edge[2]);
  TH1F *fatJet_bin_FJ3 =  new TH1F("fatJet_bin_FJ3","fatJet pT Range 3",fatJet_pT_binning_vector[3],fatJet_pT_bin_low_edge[3],fatJet_pT_bin_high_edge[3]);
  TH1F *fatJet_bin_FJ4 =  new TH1F("fatJet_bin_FJ4","fatJet pT Range 4",fatJet_pT_binning_vector[4],fatJet_pT_bin_low_edge[4],fatJet_pT_bin_high_edge[4]);
  vector<TH1F*> fatJet_pT_bins_hist_vector = {fatJet_bin_FJ0,fatJet_bin_FJ1,fatJet_bin_FJ2,fatJet_bin_FJ3,fatJet_bin_FJ4};

  // Mass Resolution Histogram stuff in fatJet Energy bins
  vector<int> fatJet_E_binning_vector = {20,20,20,20,20};
  vector<double> fatJet_E_bin_low_edge = {100,430,508,605,800};
  vector<double> fatJet_E_bin_high_edge = {430,508,605,800,2000};
  vector<double> fatJet_E_mean;
  vector<double> fatJet_E_error_low;
  vector<double> fatJet_E_error_high;

  vector<double> fatJet_bin_fatJet_mass_mean_E;
  vector<double> fatJet_bin_fatJet_mass_sigma_E;
  vector<double> fatJet_bin_fatJet_mass_sigma_error_E;

  TH1F *fatJet_bin_ME0 =  new TH1F("fatJet_bin_ME0","Fatjet mass for fatJet E Range 0",fatJet_E_binning_vector[0],65,190);
  TH1F *fatJet_bin_ME1 =  new TH1F("fatJet_bin_ME1","Fatjet mass for fatJet E Range 1",fatJet_E_binning_vector[1],65,190);
  TH1F *fatJet_bin_ME2 =  new TH1F("fatJet_bin_ME2","Fatjet mass for fatJet E Range 2",fatJet_E_binning_vector[2],65,190);
  TH1F *fatJet_bin_ME3 =  new TH1F("fatJet_bin_ME3","Fatjet mass for fatJet E Range 3",fatJet_E_binning_vector[3],65,190);
  TH1F *fatJet_bin_ME4 =  new TH1F("fatJet_bin_ME4","Fatjet mass for fatJet E Range 4",fatJet_E_binning_vector[4],65,190);
  vector<TH1F*> fatJet_E_mass_res_hist_vector = {fatJet_bin_ME0,fatJet_bin_ME1,fatJet_bin_ME2,fatJet_bin_ME3,fatJet_bin_ME4};
  TH1F *fatJet_bin_FJE0 =  new TH1F("fatJet_bin_FJE0","fatJet E Range 0",fatJet_E_binning_vector[0],fatJet_E_bin_low_edge[0],fatJet_E_bin_high_edge[0]);
  TH1F *fatJet_bin_FJE1 =  new TH1F("fatJet_bin_FJE1","fatJet E Range 1",fatJet_E_binning_vector[1],fatJet_E_bin_low_edge[1],fatJet_E_bin_high_edge[1]);
  TH1F *fatJet_bin_FJE2 =  new TH1F("fatJet_bin_FJE2","fatJet E Range 2",fatJet_E_binning_vector[2],fatJet_E_bin_low_edge[2],fatJet_E_bin_high_edge[2]);
  TH1F *fatJet_bin_FJE3 =  new TH1F("fatJet_bin_FJE3","fatJet E Range 3",fatJet_E_binning_vector[3],fatJet_E_bin_low_edge[3],fatJet_E_bin_high_edge[3]);
  TH1F *fatJet_bin_FJE4 =  new TH1F("fatJet_bin_FJE4","fatJet E Range 4",fatJet_E_binning_vector[4],fatJet_E_bin_low_edge[4],fatJet_E_bin_high_edge[4]);
  vector<TH1F*> fatJet_E_bins_hist_vector = {fatJet_bin_FJE0,fatJet_bin_FJE1,fatJet_bin_FJE2,fatJet_bin_FJE3,fatJet_bin_FJE4};

  // Mass Resolution Histogram in Higgs Pt bins
  vector<int> higgs_pT_binning_vector = {20,20,20,20,20};
  vector<double> higgs_pT_bin_low_edge = {300,365,407,449,512};
  vector<double> higgs_pT_bin_high_edge = {365,407,449,512,800};
  vector<double> higgs_pT_mean;
  vector<double> higgs_pT_error_low;
  vector<double> higgs_pT_error_high;

  vector<double> higgs_bin_fatJet_mass_mean_pT;
  vector<double> higgs_bin_fatJet_mass_sigma_pT;
  vector<double> higgs_bin_fatJet_mass_sigma_error_pT;

  TH1F *higgs_bin_M0 =  new TH1F("higgs_bin_M0","Fatjet mass for Higgs pT Range 0",higgs_pT_binning_vector[0],65,190);
  TH1F *higgs_bin_M1 =  new TH1F("higgs_bin_M1","Fatjet mass for Higgs pT Range 1",higgs_pT_binning_vector[1],65,190);
  TH1F *higgs_bin_M2 =  new TH1F("higgs_bin_M2","Fatjet mass for Higgs pT Range 2",higgs_pT_binning_vector[2],65,190);
  TH1F *higgs_bin_M3 =  new TH1F("higgs_bin_M3","Fatjet mass for Higgs pT Range 3",higgs_pT_binning_vector[3],65,190);
  TH1F *higgs_bin_M4 =  new TH1F("higgs_bin_M4","Fatjet mass for Higgs pT Range 4",higgs_pT_binning_vector[4],65,190);
  vector<TH1F*> higgs_pT_mass_res_hist_vector = {higgs_bin_M0,higgs_bin_M1,higgs_bin_M2,higgs_bin_M3,higgs_bin_M4};
  TH1F *higgs_bin_H0 =  new TH1F("higgs_bin_H0","Higgs pT Range 0",higgs_pT_binning_vector[0],higgs_pT_bin_low_edge[0],higgs_pT_bin_high_edge[0]);
  TH1F *higgs_bin_H1 =  new TH1F("higgs_bin_H1","Higgs pT Range 1",higgs_pT_binning_vector[1],higgs_pT_bin_low_edge[1],higgs_pT_bin_high_edge[1]);
  TH1F *higgs_bin_H2 =  new TH1F("higgs_bin_H2","Higgs pT Range 2",higgs_pT_binning_vector[2],higgs_pT_bin_low_edge[2],higgs_pT_bin_high_edge[2]);
  TH1F *higgs_bin_H3 =  new TH1F("higgs_bin_H3","Higgs pT Range 3",higgs_pT_binning_vector[3],higgs_pT_bin_low_edge[3],higgs_pT_bin_high_edge[3]);
  TH1F *higgs_bin_H4 =  new TH1F("higgs_bin_H4","Higgs pT Range 4",higgs_pT_binning_vector[4],higgs_pT_bin_low_edge[4],higgs_pT_bin_high_edge[4]);
  vector<TH1F*> higgs_pT_bins_hist_vector = {higgs_bin_H0,higgs_bin_H1,higgs_bin_H2,higgs_bin_H3,higgs_bin_H4};

  // Mass Resolution Histogram stuff in Higgs Energy bins
  vector<int> higgs_E_binning_vector = {20,20,20,20,20};
  vector<double> higgs_E_bin_low_edge = {100,430,508,605,800};
  vector<double> higgs_E_bin_high_edge = {430,508,605,800,2000};
  vector<double> higgs_E_mean;
  vector<double> higgs_E_error_low;
  vector<double> higgs_E_error_high;

  vector<double> higgs_bin_fatJet_mass_mean_E;
  vector<double> higgs_bin_fatJet_mass_sigma_E;
  vector<double> higgs_bin_fatJet_mass_sigma_error_E;

  TH1F *higgs_bin_ME0 =  new TH1F("higgs_bin_ME0","Fatjet mass for Higgs E Range 0",higgs_E_binning_vector[0],65,190);
  TH1F *higgs_bin_ME1 =  new TH1F("higgs_bin_ME1","Fatjet mass for Higgs E Range 1",higgs_E_binning_vector[1],65,190);
  TH1F *higgs_bin_ME2 =  new TH1F("higgs_bin_ME2","Fatjet mass for Higgs E Range 2",higgs_E_binning_vector[2],65,190);
  TH1F *higgs_bin_ME3 =  new TH1F("higgs_bin_ME3","Fatjet mass for Higgs E Range 3",higgs_E_binning_vector[3],65,190);
  TH1F *higgs_bin_ME4 =  new TH1F("higgs_bin_ME4","Fatjet mass for Higgs E Range 4",higgs_E_binning_vector[4],65,190);
  vector<TH1F*> higgs_E_mass_res_hist_vector = {higgs_bin_ME0,higgs_bin_ME1,higgs_bin_ME2,higgs_bin_ME3,higgs_bin_ME4};
  TH1F *higgs_bin_HE0 =  new TH1F("higgs_bin_HE0","Higgs E Range 0",higgs_E_binning_vector[0],higgs_E_bin_low_edge[0],higgs_E_bin_high_edge[0]);
  TH1F *higgs_bin_HE1 =  new TH1F("higgs_bin_HE1","Higgs E Range 1",higgs_E_binning_vector[1],higgs_E_bin_low_edge[1],higgs_E_bin_high_edge[1]);
  TH1F *higgs_bin_HE2 =  new TH1F("higgs_bin_HE2","Higgs E Range 2",higgs_E_binning_vector[2],higgs_E_bin_low_edge[2],higgs_E_bin_high_edge[2]);
  TH1F *higgs_bin_HE3 =  new TH1F("higgs_bin_HE3","Higgs E Range 3",higgs_E_binning_vector[3],higgs_E_bin_low_edge[3],higgs_E_bin_high_edge[3]);
  TH1F *higgs_bin_HE4 =  new TH1F("higgs_bin_HE4","Higgs E Range 4",higgs_E_binning_vector[4],higgs_E_bin_low_edge[4],higgs_E_bin_high_edge[4]);
  vector<TH1F*> higgs_E_bins_hist_vector = {higgs_bin_HE0,higgs_bin_HE1,higgs_bin_HE2,higgs_bin_HE3,higgs_bin_HE4};

  // Mass resolution in bins of opening angle between b-partons
  vector<int> partonAngle_binning_vector = {20,20,20,20,20};
  vector<double> partonAngle_bin_low_edge = {0,0.33,0.44,0.535,0.645};
  vector<double> partonAngle_bin_high_edge = {0.33,0.44,0.535,0.645,3.2};
  vector<double> partonAngle_mean;
  vector<double> partonAngle_error_low;
  vector<double> partonAngle_error_high;

  vector<double> partonAngle_bin_fatJet_mass_mean;
  vector<double> partonAngle_bin_fatJet_mass_sigma;
  vector<double> partonAngle_bin_fatJet_mass_sigma_error;

  TH1F *partonAngle_bin_M0 =  new TH1F("partonAngle_bin_M0","Fatjet mass for parton opening angle Range 0",partonAngle_binning_vector[0],65,190);
  TH1F *partonAngle_bin_M1 =  new TH1F("partonAngle_bin_M1","Fatjet mass for parton opening angle Range 1",partonAngle_binning_vector[1],65,190);
  TH1F *partonAngle_bin_M2 =  new TH1F("partonAngle_bin_M2","Fatjet mass for parton opening angle Range 2",partonAngle_binning_vector[2],65,190);
  TH1F *partonAngle_bin_M3 =  new TH1F("partonAngle_bin_M3","Fatjet mass for parton opening angle Range 3",partonAngle_binning_vector[3],65,190);
  TH1F *partonAngle_bin_M4 =  new TH1F("partonAngle_bin_M4","Fatjet mass for parton opening angle Range 4",partonAngle_binning_vector[4],65,190);
  vector<TH1F*> partonAngle_mass_res_hist_vector = {partonAngle_bin_M0,partonAngle_bin_M1,partonAngle_bin_M2,partonAngle_bin_M3,partonAngle_bin_M4};
  TH1F *partonAngle_bin_0 =  new TH1F("partonAngle_bin_0","parton opening angle Range 0",partonAngle_binning_vector[0],partonAngle_bin_low_edge[0],partonAngle_bin_high_edge[0]);
  TH1F *partonAngle_bin_1 =  new TH1F("partonAngle_bin_1","parton opening angle Range 1",partonAngle_binning_vector[1],partonAngle_bin_low_edge[1],partonAngle_bin_high_edge[1]);
  TH1F *partonAngle_bin_2 =  new TH1F("partonAngle_bin_2","parton opening angle Range 2",partonAngle_binning_vector[2],partonAngle_bin_low_edge[2],partonAngle_bin_high_edge[2]);
  TH1F *partonAngle_bin_3 =  new TH1F("partonAngle_bin_3","parton opening angle Range 3",partonAngle_binning_vector[3],partonAngle_bin_low_edge[3],partonAngle_bin_high_edge[3]);
  TH1F *partonAngle_bin_4 =  new TH1F("partonAngle_bin_4","parton opening angle Range 4",partonAngle_binning_vector[4],partonAngle_bin_low_edge[4],partonAngle_bin_high_edge[4]);
  vector<TH1F*> partonAngle_bins_hist_vector = {partonAngle_bin_0,partonAngle_bin_1,partonAngle_bin_2,partonAngle_bin_3,partonAngle_bin_4};

  // Mass Resolution Histogram in delta E of partons bins
  vector<int> deltaE_bPartons_binning_vector = {20,20,20,20,20};
  vector<double> deltaE_bPartons_bin_low_edge = {0,65,136,225,340};
  vector<double> deltaE_bPartons_bin_high_edge = {65,136,225,340,1000};
  vector<double> deltaE_bPartons_mean;
  vector<double> deltaE_bPartons_error_low;
  vector<double> deltaE_bPartons_error_high;

  vector<double> deltaE_bPartons_bin_fatJet_mass_mean;
  vector<double> deltaE_bPartons_bin_fatJet_mass_sigma;
  vector<double> deltaE_bPartons_bin_fatJet_mass_sigma_error;

  TH1F *deltaE_bPartons_bin_M0 =  new TH1F("deltaE_bPartons_bin_M0","Fatjet mass for deltaE_bPartons Range 0",deltaE_bPartons_binning_vector[0],65,190);
  TH1F *deltaE_bPartons_bin_M1 =  new TH1F("deltaE_bPartons_bin_M1","Fatjet mass for deltaE_bPartons Range 1",deltaE_bPartons_binning_vector[1],65,190);
  TH1F *deltaE_bPartons_bin_M2 =  new TH1F("deltaE_bPartons_bin_M2","Fatjet mass for deltaE_bPartons Range 2",deltaE_bPartons_binning_vector[2],65,190);
  TH1F *deltaE_bPartons_bin_M3 =  new TH1F("deltaE_bPartons_bin_M3","Fatjet mass for deltaE_bPartons Range 3",deltaE_bPartons_binning_vector[3],65,190);
  TH1F *deltaE_bPartons_bin_M4 =  new TH1F("deltaE_bPartons_bin_M4","Fatjet mass for deltaE_bPartons Range 4",deltaE_bPartons_binning_vector[4],65,190);
  vector<TH1F*> deltaE_bPartons_mass_res_hist_vector = {deltaE_bPartons_bin_M0,deltaE_bPartons_bin_M1,deltaE_bPartons_bin_M2,deltaE_bPartons_bin_M3,deltaE_bPartons_bin_M4};
  TH1F *deltaE_bPartons_bin_0 =  new TH1F("deltaE_bPartons_bin_0","deltaE_bPartons Range 0",deltaE_bPartons_binning_vector[0],deltaE_bPartons_bin_low_edge[0],deltaE_bPartons_bin_high_edge[0]);
  TH1F *deltaE_bPartons_bin_1 =  new TH1F("deltaE_bPartons_bin_1","deltaE_bPartons Range 1",deltaE_bPartons_binning_vector[1],deltaE_bPartons_bin_low_edge[1],deltaE_bPartons_bin_high_edge[1]);
  TH1F *deltaE_bPartons_bin_2 =  new TH1F("deltaE_bPartons_bin_2","deltaE_bPartons Range 2",deltaE_bPartons_binning_vector[2],deltaE_bPartons_bin_low_edge[2],deltaE_bPartons_bin_high_edge[2]);
  TH1F *deltaE_bPartons_bin_3 =  new TH1F("deltaE_bPartons_bin_3","deltaE_bPartons Range 3",deltaE_bPartons_binning_vector[3],deltaE_bPartons_bin_low_edge[3],deltaE_bPartons_bin_high_edge[3]);
  TH1F *deltaE_bPartons_bin_4 =  new TH1F("deltaE_bPartons_bin_4","deltaE_bPartons Range 4",deltaE_bPartons_binning_vector[4],deltaE_bPartons_bin_low_edge[4],deltaE_bPartons_bin_high_edge[4]);
  vector<TH1F*> deltaE_bPartons_bins_hist_vector = {deltaE_bPartons_bin_0,deltaE_bPartons_bin_1,deltaE_bPartons_bin_2,deltaE_bPartons_bin_3,deltaE_bPartons_bin_4};

  // ***** End of Resolution Plot variable declarations *****

  // 1D histograms
  TH1F *higgs_pT = new TH1F("Higgs Pt","Higgs Pt",60,0,1000);
  TH1F *higgs_pT_with_higgsFatJet = new TH1F("Higgs Pt FatJet Turn On","Higgs Pt FatJet Turn On",60,0,1000);
  TH1F *higgs_pT_with_higgsFatJet_250GeVfatJet = new TH1F("Higgs Pt FatJet Turn On fatJet pT > 250GeV","Higgs Pt FatJet Turn On fatJet pT > 250GeV",60,0,1000);
  TH1F *higgs_fatJet_angle = new TH1F("higgs_fatJet_angle","Angle between higgs and fatJet",100,0,3.2);
  TH1F *higgs_fatJet_eta = new TH1F("Eta_higgs_fatJet","Eta between higgs and Associated FatJet (fatJetEta - bosonEta)",100,-1,1);
  TH1F *b_parton_opening_angle_hist = new TH1F("b_parton_opening_angle_hist","Opening angle between b's",100,0,3.2);
  TH1F *M_fatJet = new TH1F("M_fatJet","M of Fatjet with info about semiLeptonic decay of bHad associated with Higgs Fatjet",150,0,400);
  TH1F *M_fatJet_noLep = new TH1F("M_fatJet_noLep","M of Fatjet no semi Leptonic decays",150,0,400);
  TH1F *M_fatJet_semiLep = new TH1F("M_fatJet_semiLep","M of Fatjet semi Leptonic decays",150,0,400);
  TH1F *higgs_fatJet_dR = new TH1F("higgs_fatJet_dR","delta R of Boson vs Fatjet",150,0,6.4);
  TH1F *bHadron_bParton_0_angle_hist = new TH1F("bHadron_bParton_0_angle","Angle between bHadron and bParton #0 angle",100,0,.1);
  TH1F *bHadron_bParton_1_angle_hist = new TH1F("bHadron_bParton_1_angle","Angle between bHadron and bParton #1 angle",100,0,.1);
  TH1F *b_parton_energy = new TH1F("b Parton Energy","Truth matched b parton energy summed",60,0,3000);
  TH1F *fatJet_energy = new TH1F("fatJet Energy","Higgs Associated fatJet energy",60,0,3000);
  TH1F *rand0 = new TH1F("rand0","x = 14 , y =11",100,-0.2,0.2);
  TH1F *energy_resolution = new TH1F("energy_resolution","energy_resolution",150,-1,1);
  TH1F *signal_mass = new TH1F("signal_mass","Mass of signal",150,0,300);
  TH1F *bkg_mass = new TH1F("bkg_mass","Mass of background",100,0,300);

  // 2D histograms
  TH2F *dR_partons_higgs = new TH2F("dR_partons_higgs","dR between b partons and higgs",100,0,6.4,100,0,6.4);
  TH2F *dR_hadrons_fatJet = new TH2F("dR_partons_fatJet","dR between b hadron jets and higgs associated fatJet",100,0,6.4,100,0,6.4);
  TH2F *pT_angle_partons = new TH2F("pT_angle_partons","delta pT vs opening angle of two b-partons",60,0,1.0,60,0,3.2);

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

  //*******************************     Smearing Control flags     **********************************
  //vector<float> x_axis_randN_sigma_vector {0.010,.07};
  //vector<float> y_axis_randN_sigma_vector {0.005,.05};
  //vector<float> x_axis_randN_sigma_vector {0.000001};
  //vector<float> y_axis_randN_sigma_vector {0.000001};
  //vector<float> x_axis_randN_sigma_vector {0.000001,0.001,0.03,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.3};
  //vector<float> y_axis_randN_sigma_vector {0.000001,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11};
  //vector<float> x_axis_randN_sigma_vector {0.000001,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1};
  //vector<float> y_axis_randN_sigma_vector {0.000001,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05};
  //vector<float> x_axis_randN_sigma_vector {0.000001,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.06,0.065,0.07};
  //vector<float> y_axis_randN_sigma_vector {0.000001,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05};

  vector<float> x_axis_randN_sigma_vector {0.000001,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.06,0.065,0.07};
  vector<float> y_axis_randN_sigma_vector {0.000001,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2};

  int nSmears = 88520; // Number of events that fulfil my parton-higgs association
  //int nSmears = 356900; // Number of events that fulfil my parton-higgs association
  gRandom = new TRandom3(0);
  //*************************************************************************************************

  // Stuff for making smeared delta E vs smeared opening angle of two b-partons plot
  TLorentzVector b1_4vector;
  TLorentzVector b2_4vector;
  int rand_itr = 0;
  float randN_temp = 0;
  float smeared_bPair_energy_change = 0;
  float relative_bPair_energy_change = 0;
  float smeared_bPair_angular_change = 0;
  float smeared_bPair_mass = 0;
  float bPair_mass = 0;
  float bPair_mass_change = 0;
  float smeared_bPair_E = 0;
  float bPair_E = 0;
  float sigma = 0;
  int smear_x_bins = x_axis_randN_sigma_vector.size();
  int smear_y_bins = y_axis_randN_sigma_vector.size();
  vector < vector < vector < vector<float> > > > randN_vector; // [x bins] [y bins] [b1_energy,b2_energy,b1_theta,b2_theta,b1_phi,b2_phi]
  vector < vector <TH1F*> > smeared_angle_hist_vector; // [x bins] [y bins]
  vector < vector <TH1F*> > smeared_energy_hist_vector; // [x bins] [y bins]
  vector < vector <TH1F*> > smeared_mass_hist_vector; // [x bins] [y bins]
  vector < vector <double> > angular_change_vector (9); // Vector of vectors which will contain the mass resolution for selected populations
  vector < vector <double> > energy_change_vector (9); // Vector of vectors which will contain the mass resolution for selected populations
  TH2F *smeared_bPartons = new TH2F("smeared_bPartons","fatJet Mass resolution for smeared b-partons",smear_x_bins,0,1,smear_y_bins,0,1);


  // Generate Random Numbers for b parton angle and energy smearing plot
  for(int x_itr = 0; x_itr < smear_x_bins; x_itr++){
    smeared_angle_hist_vector.push_back(vector <TH1F*>());
    smeared_energy_hist_vector.push_back(vector <TH1F*>());
    smeared_mass_hist_vector.push_back(vector <TH1F*>());
    randN_vector.push_back(vector < vector < vector<float> > >());
    for(int y_itr = 0; y_itr < smear_y_bins; y_itr++) {
      ostringstream mass_histogramNameStream;
      mass_histogramNameStream << "Smeared mass bin: x =" << x_itr+1 << ", y = " << y_itr+1;
      smeared_mass_hist_vector[x_itr].push_back(new TH1F(mass_histogramNameStream.str().c_str(),mass_histogramNameStream.str().c_str(),100,0,500));

      ostringstream angle_histogramNameStream;
      angle_histogramNameStream << "Smeared angle bin: x =" << x_itr+1 << ", y = " << y_itr+1;
      smeared_angle_hist_vector[x_itr].push_back(new TH1F(angle_histogramNameStream.str().c_str(),angle_histogramNameStream.str().c_str(),100,-1,1));

      ostringstream energy_histogramNameStream;
      energy_histogramNameStream << "Smeared energy bin: x =" << x_itr+1 << ", y = " << y_itr+1;
      smeared_energy_hist_vector[x_itr].push_back(new TH1F(energy_histogramNameStream.str().c_str(),energy_histogramNameStream.str().c_str(),100,-0.2,0.2));

      randN_vector[x_itr].push_back(vector < vector<float> >());
      for(int rand_var_itr = 0; rand_var_itr < 6; rand_var_itr++){
        if(rand_var_itr < 2) {
        // Smearing energy corresponds to the y-axis
          sigma = y_axis_randN_sigma_vector[y_itr];
        }
        else {
        // rand_var_itr 2 or greater corresponds to smearing angle which is on the x-axis
          sigma = x_axis_randN_sigma_vector[x_itr];
        }

        randN_vector[x_itr][y_itr].push_back(vector<float>());
        for(int rand_itr = 0; rand_itr < nSmears; rand_itr++) {
          // [x bins] [y bins] [b1_energy,b2_energy,b1_theta,b2_theta,b1_phi,b2_phi]
          randN_temp = gRandom->Gaus(0,sigma);
          randN_vector[x_itr][y_itr][rand_var_itr].push_back(randN_temp);
        }
      }
    }
  } // End of aandom number generation for smearing

  //*******************************     Controlling flags     **********************************
  int higgs_fatJet_association_flag = 1; // 1 -> My Association, 0 -> fatJetNGhostH association
  int doResolutionPlots = 0; // Turn on or off fatJet mass resolution plots
  int doSmearPlot = 0; // Turn on or off the bParton smearing for simulated mass resolution plots
  int doSensitivityStudy = 0; // Turn on or off the sensitivity studyin the choosen mass peak region
  //********************************************************************************************
  if(doSensitivityStudy == 1) {doSmearPlot = 1;} // The sensitivity study uses the mass resolutions found in the smeared_bParton mass matrix

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
        smeared_bPair_mass = 0;
        bPair_mass = 0;
        bPair_mass_change = 0;
        smeared_bPair_E = 0;
        bPair_E = 0;
        smeared_bPair_energy_change = 0;
        relative_bPair_energy_change = 0;
        smeared_bPair_angular_change = 0;
        b1_4vector.SetPxPyPzE(0,0,0,0);
        b2_4vector.SetPxPyPzE(0,0,0,0);


        if(bosonPdgId->at(boson_itr) == 25 && bosonStatus->at(boson_itr) == 62) {
          higgs_pT->Fill(bosonPt->at(boson_itr)); // Used for making higgs-fatJet associated turn on curve

          boson3vector.SetXYZ(bosonPx->at(boson_itr),bosonPy->at(boson_itr),bosonPz->at(boson_itr));
          boson4vector.SetPxPyPzE(bosonPx->at(boson_itr),bosonPy->at(boson_itr),bosonPz->at(boson_itr),bosonE->at(boson_itr));
          bosonPhi = boson4vector.Phi();
          bosonEta = boson4vector.Eta();

          // Parton Loop
          for(int parton_itr = 0; parton_itr < partonPt->size(); parton_itr++) {
            parton3vector.SetXYZ(partonPx->at(parton_itr),partonPy->at(parton_itr),partonPz->at(parton_itr));
            parton4vector.SetPxPyPzE(partonPx->at(parton_itr),partonPy->at(parton_itr),partonPz->at(parton_itr),partonE->at(parton_itr));
            partonPhi = parton4vector.Phi();
            partonEta = parton4vector.Eta();

            if(abs(partonPdgId->at(parton_itr)) == 5 && partonParentPdgId->at(parton_itr) == 25) {
              bParton_energy_sum = bParton_energy_sum + partonE->at(parton_itr);
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
          if(b_parton3vectors.size() > 1) {
            b_parton_energy->Fill(bParton_energy_sum);
            bParton_energy_sum = 0;
            dR_parton0_higgs_container = b_parton4vectors.at(0).DeltaR(boson4vector);
            dR_parton1_higgs_container = b_parton4vectors.at(1).DeltaR(boson4vector);

            if(b_parton4vectors.at(0).Pt() > b_parton4vectors.at(1).Pt()) {
              dR_partons_higgs->Fill(dR_parton0_higgs_container,dR_parton1_higgs_container);
            }
            else{
              dR_partons_higgs->Fill(dR_parton1_higgs_container,dR_parton0_higgs_container);
            }

            deltaE_bPartons = abs(b_parton4vectors.at(0).E() - b_parton4vectors.at(1).E());
            b_parton_opening_angle = b_parton3vectors.at(0).Angle(b_parton3vectors.at(1));
            b_parton_opening_angle_hist->Fill(b_parton_opening_angle);
            deltaPt_b_partons = abs(b_parton4vectors.at(0).Pt() - b_parton4vectors.at(1).Pt()) / bosonPt->at(boson_itr);
            pT_angle_partons->Fill(deltaPt_b_partons,b_parton_opening_angle);

            // Smearing Histogram stuff
            bPair_mass = (b_parton4vectors.at(0) + b_parton4vectors.at(1)).M();
            bPair_E = b_parton4vectors.at(0).E() + b_parton4vectors.at(1).E();
            if(doSmearPlot == 1) {
              for(int smear_itr = 0; smear_itr < 10; smear_itr++) {
                // This loop does a repeated smearing of the same event to increase statistics in the result
                // Makes sure to never re-use a random number by iterating rand_itr at the end
                for(int x_itr = 0; x_itr < smear_x_bins; x_itr++){
                  for(int y_itr = 0; y_itr < smear_y_bins; y_itr++) {
                    // Smear the b-pair's 6 numbers by the number stored at rand_itr in the randN_vector storage vector
                    b1_4vector = b_parton4vectors.at(0);
                    b2_4vector = b_parton4vectors.at(1);

                    // Smear each energy and produced smeared delta Energy
                    //randN_vector [x bins] [y bins] [b1_energy,b2_energy,b1_theta,b2_theta,b1_phi,b2_phi] 
                    //b1_4vector = b1_4vector * (1 + randN_vector[x_itr][y_itr][0].at(rand_itr));
                    //b2_4vector = b2_4vector * (1 + randN_vector[x_itr][y_itr][1].at(rand_itr));
                    smear_E(&b1_4vector, randN_vector[x_itr][y_itr][0].at(rand_itr));
                    smear_E(&b2_4vector, randN_vector[x_itr][y_itr][1].at(rand_itr));

        
                    // Smear both angles and produce smeared opening angle
                    b1_4vector.SetTheta(b1_4vector.Theta() + randN_vector[x_itr][y_itr][2].at(rand_itr));
                    b2_4vector.SetTheta(b2_4vector.Theta() + randN_vector[x_itr][y_itr][3].at(rand_itr));
                    b1_4vector.SetPhi(b1_4vector.Phi() + randN_vector[x_itr][y_itr][4].at(rand_itr));
                    b2_4vector.SetPhi(b2_4vector.Phi() + randN_vector[x_itr][y_itr][5].at(rand_itr));
                    
                    if(b1_4vector.DeltaR(b2_4vector) < 1.0) {
                      if(x_itr == 13 && y_itr == 10){
                        rand0->Fill(randN_vector[x_itr][y_itr][0].at(rand_itr));
                        rand0->SetTitle(Form("RandN Gaussian for energy smearing with width: (%f)",y_axis_randN_sigma_vector[y_itr]));
                      }

                      smeared_bPair_mass = (b1_4vector + b2_4vector).M();
                      smeared_mass_hist_vector[x_itr][y_itr]->Fill(smeared_bPair_mass);

                      smeared_bPair_angular_change = b1_4vector.Angle(b2_4vector.Vect()) - b_parton_opening_angle;
                      smeared_bPair_energy_change = (b1_4vector.E() + b2_4vector.E()) - bPair_E; 
                      relative_bPair_energy_change = ((b1_4vector.E() + b2_4vector.E()) - bPair_E) / bPair_E;
                      smeared_angle_hist_vector[x_itr][y_itr]->Fill(smeared_bPair_angular_change);
                      smeared_energy_hist_vector[x_itr][y_itr]->Fill(relative_bPair_energy_change);

                      // Stuff for making TGraph plot of pairs colored how much their mass changed from the truth
                      // This will show whether it was the angle or energy smearing that contributed most to their change in mass
                      bPair_mass_change = abs(smeared_bPair_mass - bPair_mass);
                      //smeared_bPair_energy_change = abs(smeared_bPair_energy_change); // Make positive for TGraph plot
                      //smeared_bPair_angular_change = abs(smeared_bPair_angular_change); // Make positive for TGraph plot

                      if(bPair_mass_change < 5) {
                        angular_change_vector[0].push_back(smeared_bPair_angular_change);
                        energy_change_vector[0].push_back(relative_bPair_energy_change);
                      }
                      else if(bPair_mass_change < 10) {
                        angular_change_vector[1].push_back(smeared_bPair_angular_change);
                        energy_change_vector[1].push_back(relative_bPair_energy_change);
                      }
                      else if(bPair_mass_change < 20) {
                        angular_change_vector[2].push_back(smeared_bPair_angular_change);
                        energy_change_vector[2].push_back(relative_bPair_energy_change);
                      }
                      else if(bPair_mass_change < 30) {
                        angular_change_vector[3].push_back(smeared_bPair_angular_change);
                        energy_change_vector[3].push_back(relative_bPair_energy_change);
                      }
                      else if(bPair_mass_change < 40) {
                        angular_change_vector[4].push_back(smeared_bPair_angular_change);
                        energy_change_vector[4].push_back(relative_bPair_energy_change);
                      }
                      else if(bPair_mass_change < 50) {
                        angular_change_vector[5].push_back(smeared_bPair_angular_change);
                        energy_change_vector[5].push_back(relative_bPair_energy_change);
                      }
                        else if(bPair_mass_change < 60) {
                        angular_change_vector[6].push_back(smeared_bPair_angular_change);
                        energy_change_vector[6].push_back(relative_bPair_energy_change);
                      }
                      else if(bPair_mass_change < 70) {
                        angular_change_vector[7].push_back(smeared_bPair_angular_change);
                        energy_change_vector[7].push_back(relative_bPair_energy_change);
                      }
                      else if(bPair_mass_change < 80) {
                        angular_change_vector[8].push_back(smeared_bPair_angular_change);
                        energy_change_vector[8].push_back(relative_bPair_energy_change);
                      }
                    }
                  } // End of y_axis for loop
                } // End of x_axis for loop          
                rand_itr++;
              }
            } // End of doSmearPlot if statement
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
              bHadron_bParton_0_angle_hist->Fill(b_parton3vectors[0].Angle(bHadron3vector_paired_to_bParton0));
              bHadron_bParton_1_angle_hist->Fill(b_parton3vectors[1].Angle(bHadron3vector_paired_to_bParton1));

              // Store information about which bHadron's decay semi-leptonically so we can display info on fatJet mass plot
              bHad_NLepSec.push_back(trackJetBHadronNLepSec->at(bHadron_paired_to_bParton0));
              bHad_NLepSec.push_back(trackJetBHadronNLepSec->at(bHadron_paired_to_bParton1));
              bHad_NLepTer.push_back(trackJetBHadronNLepTer->at(bHadron_paired_to_bParton0));
              bHad_NLepTer.push_back(trackJetBHadronNLepTer->at(bHadron_paired_to_bParton1));

              if((bHad_NLepSec.at(0) + bHad_NLepTer.at(0) + bHad_NLepSec.at(1) + bHad_NLepTer.at(1)) > 0) {
                SemiLepBhadDecay_flag = 1; // Includes semi leptonic decay of bHadron tagged to Higgs
              }
              else {
                SemiLepBhadDecay_flag = 0; // Non-leptonic decay of bHadron tagged to higgs
              }
            }
          }

          // FatJet Loop
          for(int fatJet_itr = 0; fatJet_itr < fatJetPt->size(); fatJet_itr++) {
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
            // For fatJet turn on curve without fatJetPt > 250GeV cut
            if(HiggsAssociationFlag == 1) {
              higgs_pT_with_higgsFatJet->Fill(bosonPt->at(boson_itr)); 
            }

            if(fatJetPt->at(fatJet_itr) > 250 && HiggsAssociationFlag == 1) {
              energy_resolution->Fill((fatJetE->at(fatJet_itr) - bPair_E) / bPair_E);
              fatJet_energy->Fill(fatJetE->at(fatJet_itr)); // For comparing fatJet Energy to energy of 2 higgs associated trackJetBHadrons
              higgs_pT_with_higgsFatJet_250GeVfatJet->Fill(bosonPt->at(boson_itr)); // For fatJet turn on curve including the fatJetPt > 250GeV cut

              fatJet3vector.SetXYZ(fatJetPx->at(fatJet_itr),fatJetPy->at(fatJet_itr),fatJetPz->at(fatJet_itr));
              fatJet4vector.SetPxPyPzE(fatJetPx->at(fatJet_itr),fatJetPy->at(fatJet_itr),fatJetPz->at(fatJet_itr),fatJetE->at(fatJet_itr));
              fatJetPhi = fatJet4vector.Phi();
              fatJetEta = fatJet4vector.Eta();

              // dR between fatJet and b Hadron jets
              dR_hadron0_fatJet_container = bHadron4vector_paired_to_bParton0.DeltaR(fatJet4vector);
              dR_hadron1_fatJet_container = bHadron4vector_paired_to_bParton1.DeltaR(fatJet4vector);
              if(bHadron4vector_paired_to_bParton0.Pt() > bHadron4vector_paired_to_bParton1.Pt()) {
                dR_hadrons_fatJet->Fill(dR_hadron0_fatJet_container,dR_hadron1_fatJet_container);
              }
              else{
                dR_hadrons_fatJet->Fill(dR_hadron1_fatJet_container,dR_hadron0_fatJet_container);
              }
              // Angle between Fatjet and Higgs
              angle = fatJet3vector.Angle(boson3vector);
              higgs_fatJet_angle->Fill(angle);
              // Mass of ghost associated fatJets
              M_fatJet->Fill(fatJetM->at(fatJet_itr));
              if(SemiLepBhadDecay_flag == 1) { // One or both bHadron(s) associated with higgs fatJet has a semi-leptonic decay
                M_fatJet_semiLep->Fill(fatJetM->at(fatJet_itr));
              }
              else {
                M_fatJet_noLep->Fill(fatJetM->at(fatJet_itr));
              }
              // DeltaR for ghost associated fatJets and higgs
              dR = fatJet4vector.DeltaR(boson4vector);
              higgs_fatJet_dR->Fill(dR);
              // Eta for ghost associated fatJet and higgs
              higgs_fatJet_eta->Fill(fatJetEta-bosonEta);

              if(fatJetM->at(fatJet_itr) > 65 && fatJetM->at(fatJet_itr) < 190 ) {
                // Mass Resolution of Higgs associated fatJets in delta E of partons bins (with mass window)
                if(deltaE_bPartons > deltaE_bPartons_bin_low_edge[0] && deltaE_bPartons < deltaE_bPartons_bin_high_edge[0]) { 
                  deltaE_bPartons_bin_M0->Fill(fatJetM->at(fatJet_itr));
                  deltaE_bPartons_bin_0->Fill(deltaE_bPartons);
                }
                else if(deltaE_bPartons > deltaE_bPartons_bin_low_edge[1] && deltaE_bPartons < deltaE_bPartons_bin_high_edge[1]) { 
                  deltaE_bPartons_bin_M1->Fill(fatJetM->at(fatJet_itr));
                  deltaE_bPartons_bin_1->Fill(deltaE_bPartons);
                }
                else if(deltaE_bPartons > deltaE_bPartons_bin_low_edge[2] && deltaE_bPartons < deltaE_bPartons_bin_high_edge[2]) { 
                  deltaE_bPartons_bin_M2->Fill(fatJetM->at(fatJet_itr));
                  deltaE_bPartons_bin_2->Fill(deltaE_bPartons);
                }
                else if(deltaE_bPartons > deltaE_bPartons_bin_low_edge[3] && deltaE_bPartons < deltaE_bPartons_bin_high_edge[3]) { 
                  deltaE_bPartons_bin_M3->Fill(fatJetM->at(fatJet_itr));
                  deltaE_bPartons_bin_3->Fill(deltaE_bPartons);
                }
                else if(deltaE_bPartons > deltaE_bPartons_bin_low_edge[4] && deltaE_bPartons <  deltaE_bPartons_bin_high_edge[4]) { 
                  deltaE_bPartons_bin_M4->Fill(fatJetM->at(fatJet_itr));
                  deltaE_bPartons_bin_4->Fill(deltaE_bPartons);
                }

                // Mass resolution of Higgs associated fatJets in fatJet parton opening angle bins (with mass window)
                if(b_parton_opening_angle > partonAngle_bin_low_edge[0] && b_parton_opening_angle < partonAngle_bin_high_edge[0]) { 
                  partonAngle_bin_M0->Fill(fatJetM->at(fatJet_itr));
                  partonAngle_bin_0->Fill(b_parton_opening_angle);
                }
                else if(b_parton_opening_angle > partonAngle_bin_low_edge[1] && b_parton_opening_angle < partonAngle_bin_high_edge[1]) { 
                  partonAngle_bin_M1->Fill(fatJetM->at(fatJet_itr));
                  partonAngle_bin_1->Fill(b_parton_opening_angle);
                }
                else if(b_parton_opening_angle > partonAngle_bin_low_edge[2] && b_parton_opening_angle < partonAngle_bin_high_edge[2]) { 
                  partonAngle_bin_M2->Fill(fatJetM->at(fatJet_itr));
                  partonAngle_bin_2->Fill(b_parton_opening_angle);
                }
                else if(b_parton_opening_angle > partonAngle_bin_low_edge[3] && b_parton_opening_angle < partonAngle_bin_high_edge[3]) { 
                  partonAngle_bin_M3->Fill(fatJetM->at(fatJet_itr));
                  partonAngle_bin_3->Fill(b_parton_opening_angle);
                }
                else if(b_parton_opening_angle > partonAngle_bin_low_edge[4] && b_parton_opening_angle <  partonAngle_bin_high_edge[4]) { 
                  partonAngle_bin_M4->Fill(fatJetM->at(fatJet_itr));
                  partonAngle_bin_4->Fill(b_parton_opening_angle);
                }

                // Mass resolution of Higgs associated fatJets in fatJet pT bins (with mass window)
                if(fatJetPt->at(fatJet_itr) > fatJet_pT_bin_low_edge[0] && fatJetPt->at(fatJet_itr) < fatJet_pT_bin_high_edge[0]) { 
                  fatJet_bin_M0->Fill(fatJetM->at(fatJet_itr));
                  fatJet_bin_FJ0->Fill(fatJetPt->at(fatJet_itr));
                }
                else if(fatJetPt->at(fatJet_itr) > fatJet_pT_bin_low_edge[1] && fatJetPt->at(fatJet_itr) < fatJet_pT_bin_high_edge[1]) { 
                  fatJet_bin_M1->Fill(fatJetM->at(fatJet_itr));
                  fatJet_bin_FJ1->Fill(fatJetPt->at(fatJet_itr));
                }
                else if(fatJetPt->at(fatJet_itr) > fatJet_pT_bin_low_edge[2] && fatJetPt->at(fatJet_itr) < fatJet_pT_bin_high_edge[2]) { 
                  fatJet_bin_M2->Fill(fatJetM->at(fatJet_itr));
                  fatJet_bin_FJ2->Fill(fatJetPt->at(fatJet_itr));
                }
                else if(fatJetPt->at(fatJet_itr) > fatJet_pT_bin_low_edge[3] && fatJetPt->at(fatJet_itr) < fatJet_pT_bin_high_edge[3]) { 
                  fatJet_bin_M3->Fill(fatJetM->at(fatJet_itr));
                  fatJet_bin_FJ3->Fill(fatJetPt->at(fatJet_itr));
                }
                else if(fatJetPt->at(fatJet_itr) > fatJet_pT_bin_low_edge[4] && fatJetPt->at(fatJet_itr) <  fatJet_pT_bin_high_edge[4]) { 
                  fatJet_bin_M4->Fill(fatJetM->at(fatJet_itr));
                  fatJet_bin_FJ4->Fill(fatJetPt->at(fatJet_itr));
                }

                // Mass resolution of Higgs associated fatJets in fatJet Energy bins (with mass window)
                if(fatJetE->at(fatJet_itr) > fatJet_E_bin_low_edge[0] && fatJetE->at(fatJet_itr) < fatJet_E_bin_high_edge[0]) { 
                  fatJet_bin_ME0->Fill(fatJetM->at(fatJet_itr));
                  fatJet_bin_FJE0->Fill(fatJetE->at(fatJet_itr));
                }
                else if(fatJetE->at(fatJet_itr) > fatJet_E_bin_low_edge[1] && fatJetE->at(fatJet_itr) < fatJet_E_bin_high_edge[1]) { 
                  fatJet_bin_ME1->Fill(fatJetM->at(fatJet_itr));
                  fatJet_bin_FJE1->Fill(fatJetE->at(fatJet_itr));
                }
                else if(fatJetE->at(fatJet_itr) > fatJet_E_bin_low_edge[2] && fatJetE->at(fatJet_itr) < fatJet_E_bin_high_edge[2]) { 
                  fatJet_bin_ME2->Fill(fatJetM->at(fatJet_itr));
                  fatJet_bin_FJE2->Fill(fatJetE->at(fatJet_itr));
                }
                else if(fatJetE->at(fatJet_itr) > fatJet_E_bin_low_edge[3] && fatJetE->at(fatJet_itr) < fatJet_E_bin_high_edge[3]) { 
                  fatJet_bin_ME3->Fill(fatJetM->at(fatJet_itr));
                  fatJet_bin_FJE3->Fill(fatJetE->at(fatJet_itr));
                }
                else if(fatJetE->at(fatJet_itr) > fatJet_E_bin_low_edge[4] && fatJetE->at(fatJet_itr) <  fatJet_E_bin_high_edge[4]) { 
                  fatJet_bin_ME4->Fill(fatJetM->at(fatJet_itr));
                  fatJet_bin_FJE4->Fill(fatJetE->at(fatJet_itr));
                }

                // Mass resolution of Higgs associated fatJets in Higgs pT bins (with mass window)
                if(bosonPt->at(boson_itr) > higgs_pT_bin_low_edge[0] && bosonPt->at(boson_itr) < higgs_pT_bin_high_edge[0]) { 
                  higgs_bin_M0->Fill(fatJetM->at(fatJet_itr));
                  higgs_bin_H0->Fill(bosonPt->at(boson_itr));
               }
                else if(bosonPt->at(boson_itr) > higgs_pT_bin_low_edge[1] && bosonPt->at(boson_itr) < higgs_pT_bin_high_edge[1]) { 
                  higgs_bin_M1->Fill(fatJetM->at(fatJet_itr));
                  higgs_bin_H1->Fill(bosonPt->at(boson_itr));
                }
                else if(bosonPt->at(boson_itr) > higgs_pT_bin_low_edge[2] && bosonPt->at(boson_itr) < higgs_pT_bin_high_edge[2]) { 
                  higgs_bin_M2->Fill(fatJetM->at(fatJet_itr));
                  higgs_bin_H2->Fill(bosonPt->at(boson_itr));
                }
                else if(bosonPt->at(boson_itr) > higgs_pT_bin_low_edge[3] && bosonPt->at(boson_itr) < higgs_pT_bin_high_edge[3]) { 
                  higgs_bin_M3->Fill(fatJetM->at(fatJet_itr));
                  higgs_bin_H3->Fill(bosonPt->at(boson_itr));
                }
                else if(bosonPt->at(boson_itr) > higgs_pT_bin_low_edge[4] && bosonPt->at(boson_itr) <  higgs_pT_bin_high_edge[4]) { 
                  higgs_bin_M4->Fill(fatJetM->at(fatJet_itr));
                  higgs_bin_H4->Fill(bosonPt->at(boson_itr));
                }

                // Mass resolution of Higgs associated fatJets in Higgs Energy bins (with mass window)
                if(bosonE->at(boson_itr) > higgs_E_bin_low_edge[0] && bosonE->at(boson_itr) < higgs_E_bin_high_edge[0]) { 
                  higgs_bin_ME0->Fill(fatJetM->at(fatJet_itr));
                  higgs_bin_HE0->Fill(bosonE->at(boson_itr));
                }
                else if(bosonE->at(boson_itr) > higgs_E_bin_low_edge[1] && bosonE->at(boson_itr) < higgs_E_bin_high_edge[1]) { 
                  higgs_bin_ME1->Fill(fatJetM->at(fatJet_itr));
                  higgs_bin_HE1->Fill(bosonE->at(boson_itr));
                }
                else if(bosonE->at(boson_itr) > higgs_E_bin_low_edge[2] && bosonE->at(boson_itr) < higgs_E_bin_high_edge[2]) { 
                  higgs_bin_ME2->Fill(fatJetM->at(fatJet_itr));
                  higgs_bin_HE2->Fill(bosonE->at(boson_itr));
                }
                else if(bosonE->at(boson_itr) > higgs_E_bin_low_edge[3] && bosonE->at(boson_itr) < higgs_E_bin_high_edge[3]) { 
                  higgs_bin_ME3->Fill(fatJetM->at(fatJet_itr));
                  higgs_bin_HE3->Fill(bosonE->at(boson_itr));
                }
                else if(bosonE->at(boson_itr) > higgs_E_bin_low_edge[4] && bosonE->at(boson_itr) <  higgs_E_bin_high_edge[4]) { 
                  higgs_bin_ME4->Fill(fatJetM->at(fatJet_itr));
                  higgs_bin_HE4->Fill(bosonE->at(boson_itr));
                }
              }
            }
          } // End of Fatjet loop

          // Study of rate of b-Hadrons that fall outside of the fatJet in bins of higgs pT
          if(b_parton3vectors.size() == 2) {
            if(bosonPt->at(boson_itr) > higgs_pT_bin_low_edge[0] && bosonPt->at(boson_itr) < higgs_pT_bin_high_edge[0]) { 
              bin_iterator = 0;
            }
            else if(bosonPt->at(boson_itr) > higgs_pT_bin_low_edge[1] && bosonPt->at(boson_itr) < higgs_pT_bin_high_edge[1]) { 
              bin_iterator = 1;
            }
            else if(bosonPt->at(boson_itr) > higgs_pT_bin_low_edge[2] && bosonPt->at(boson_itr) < higgs_pT_bin_high_edge[2]) { 
              bin_iterator = 2;
            }
            else if(bosonPt->at(boson_itr) > higgs_pT_bin_low_edge[3] && bosonPt->at(boson_itr) < higgs_pT_bin_high_edge[3]) { 
              bin_iterator = 3;
            }
            else if(bosonPt->at(boson_itr) > higgs_pT_bin_low_edge[4] && bosonPt->at(boson_itr) < higgs_pT_bin_high_edge[4]) { 
              bin_iterator = 4;
            }
            // Fill appropriate bin (bin_iterator) of vectors
            if(bHadron_paired_to_bParton0 != -1 && bHadron_paired_to_bParton1 != -1) {
              n_higgs_ptBins[bin_iterator]++; // Number of events where we pair both b's from the higgs to a bHadron
            }
            if( bHadron_paired_to_bParton0 == -1 && bHadron_paired_to_bParton1 == -1) {
              n_noMatch_in_ptBin[bin_iterator]++; // Number of events where neither b from the higgs is paired to a bHadron
            }
            else if( bHadron_paired_to_bParton0 == -1 || bHadron_paired_to_bParton1 == -1) {
              n_oneMatch_ptBins[bin_iterator]++; // Number of events where one b from the higgs is paired to a bHadron
            }
            else if(trackJetIdFatJet->at(bHadron_paired_to_bParton0) != trackJetIdFatJet->at(bHadron_paired_to_bParton1)) {
              counter++;
              n_hadron_outside_ptBins[bin_iterator]++; // Number of events where both b's are paired to a bHadron but are associated with different fatJets
            }
          }
        } // End of if statement checking for higgs PdgId and Status Code = 62
      } // End of boson loop
    } // End of trigger if statement
  } // End of tree loop

  // ******************************************************
  // *****           Start of Plotting Area           *****
  // ******************************************************
  
  TCanvas *c_pT_angle_partons = new TCanvas("pT_angle_partons","delta pT vs angle of b-partons",800,800);
  pT_angle_partons->GetXaxis()->SetTitle("(higgs b-parton delta pT) / (higgs pT)");
  pT_angle_partons->GetYaxis()->SetTitle("angle between higgs b-partons");
  pT_angle_partons->Draw("COLZ");
  c_pT_angle_partons->SetRightMargin(0.18);
  c_pT_angle_partons->SetLeftMargin(0.18);
  c_pT_angle_partons->SaveAs("pT_angle_partons.pdf");

  TCanvas *c_dR_partons_higgs = new TCanvas("dR_partons_higgs","dR_partons_higgs",800,800);
  dR_partons_higgs->GetXaxis()->SetTitle("dR of higher Pt b parton");
  dR_partons_higgs->GetYaxis()->SetTitle("dR of lower Pt b parton");
  dR_partons_higgs->Draw("COLZ");
  c_dR_partons_higgs->SetRightMargin(0.18);
  c_dR_partons_higgs->SetLeftMargin(0.18);
  c_dR_partons_higgs->SaveAs("dR_partons_higgs.pdf");

  TCanvas *c_dR_hadrons_fatJet = new TCanvas("dR_hadrons_fatJet","dR_hadrons_fatJet",800,800);
  dR_hadrons_fatJet->GetXaxis()->SetTitle("dR of higher Pt b hadron jet");
  dR_hadrons_fatJet->GetYaxis()->SetTitle("dR of lower Pt b hadron jet");
  dR_hadrons_fatJet->Draw("COLZ");
  c_dR_hadrons_fatJet->SetRightMargin(0.18);
  c_dR_hadrons_fatJet->SetLeftMargin(0.18);
  c_dR_hadrons_fatJet->SaveAs("dR_hadrons_fatJet.pdf");

  TCanvas *c_compare_fatJetE_sumBpartonE = new TCanvas("E [GeV]","compare fatJet E to the sum of truth Higgs b Parton E",800,600);
  c_compare_fatJetE_sumBpartonE->Divide(2,1);
  c_compare_fatJetE_sumBpartonE->cd(1);
  b_parton_energy->SetTitle("#splitline{Comparison of fatJet (Red) and summed}{higgs b-partons (Blue) Energy}");
  b_parton_energy->SetLineColor(4);
  b_parton_energy->Draw();
  fatJet_energy->SetLineColor(2);
  fatJet_energy->Draw("same");
  c_compare_fatJetE_sumBpartonE->cd(2);
  TH1F *normalized_b_parton_energy = (TH1F*)b_parton_energy->Clone();
  TH1F *normalized_fatJet_energy = (TH1F*)fatJet_energy->Clone();
  float b_parton_integral = 0;
  float fatJet_integral = 0;
  b_parton_integral = normalized_b_parton_energy->GetEntries();
  fatJet_integral = normalized_fatJet_energy->GetEntries();
  normalized_b_parton_energy->Scale(1/b_parton_integral);
  normalized_fatJet_energy->Scale(1/fatJet_integral);
  normalized_fatJet_energy->SetTitle("Same as left but both nomralized to 1");
  normalized_fatJet_energy->Draw();
  normalized_b_parton_energy->Draw("same");
  c_compare_fatJetE_sumBpartonE->SaveAs("fatJet_bPartons_energy_comparison.pdf");

  TCanvas *c_m_fatJet = new TCanvas("M_FatJet","Mass of FatJet plot",800,600);
  TLegend *leg_m_fatJet = new TLegend(0.51,0.6,0.95,0.75);
  leg_m_fatJet->AddEntry(M_fatJet,"Total");
  leg_m_fatJet->AddEntry(M_fatJet_semiLep,"SemiLeptonic bHad Decays");
  leg_m_fatJet->AddEntry(M_fatJet_noLep,"No Leptonic bHad Decays");
  M_fatJet_semiLep->SetLineColor(4);
  M_fatJet_noLep->SetLineColor(2);
  M_fatJet->SetLineColor(1);
  M_fatJet->GetXaxis()->SetTitle("fatJet mass [GeV]");
  M_fatJet->Draw();
  M_fatJet_semiLep->Draw("same");

  M_fatJet_noLep->Draw("same");
  leg_m_fatJet->Draw("same");
  c_m_fatJet->SaveAs("M_FatJet.pdf","recreate");

  // Function for plotting TH1F
  // void plot_TH1F(char* output_filename, TH1F *hist, char* title, char* x_axis, char* y_axis, char* options, int canvas_width, int canvas_height){
  gStyle->SetOptStat(1);
  plot_TH1F(bHadron_bParton_0_angle_hist,"bHadron_bParton_0_angle_hist.pdf","Angle Vetween bHadron and bParton #0","Angle [rad]","#","",800,600);
  plot_TH1F(bHadron_bParton_1_angle_hist,"bHadron_bParton_1_angle_hist.pdf","Angle Vetween bHadron and bParton #1","Angle [rad]","#","",800,600);
  plot_TH1F(higgs_fatJet_eta,"higgs_fatJet_eta.pdf","Eta between truth higgs and associated fatJet","eta","#","",800,600);
  plot_TH1F(higgs_fatJet_angle,"higgs_fatJet_angle.pdf","Angle between truth higgs and associated fatJet","Angle [rad]","#","",800,600);
  plot_TH1F(higgs_fatJet_dR,"higgs_fatJet_dR.pdf","dR between truth higgs and associated fatJet","dR","#","",800,600);
  plot_TH1F(b_parton_opening_angle_hist,"b_parton_opening_angle.pdf","Opening angle between b's","Angle [rad]","#","",800,600);
  plot_TH1F(rand0,"rand0.pdf","rand0","rand N","#","",800,600);
  gStyle->SetOptStat(0);

  //Energy_resolution from reconstruction
  TCanvas *c_energy_resolution = new TCanvas("energy_resolution","energy_resolution",800,600);

  // *** Fit fatJet Energy Resolution distribution using a bukin function***
  //RooRealVar energy_resolution_var("energy_resolution","energy_resolution",energy_resolution->GetXaxis()->GetXmin() ,energy_resolution->GetXaxis()->GetXmax());
  RooRealVar energy_resolution_var("energy_resolution","energy_resolution", -0.5 , 0.5);

  // Declare histogram to be fit along and a variable (smeared_Mass) to be filled with fit results
  RooDataHist er_hist("er_hist","histogram of smeared Mass",energy_resolution_var,energy_resolution);

  // Declare bukin variables and function. Don't forget to give variables ranges or they will be considered constant in the fit
  RooRealVar er_peak_position("peak_position","location of peak on x axis",0,-1,1);
  RooRealVar er_sigma("sigma","resolution of Mass defined as: FWHM divided by 2*sqrt(2*log(2))=2.35",10,0,30);
  RooRealVar er_peak_asymmetry("peak_asymmetry","peak asymmetry parameter",0,-0.5,0.5);
  RooRealVar er_left_tail("left_tail","parameter of the left tail",0,-0.5,0.5);
  RooRealVar er_right_tail("right_tail","parameter of the right tail",0,-0.5,0.5);
  RooBukinPdf er_bukin("my bukin name","my bukin title",energy_resolution_var,er_peak_position,er_sigma,er_peak_asymmetry,er_left_tail,er_right_tail);

  // Create frame for RooFit plot result to be plotted on
  RooPlot* er_frame = energy_resolution_var.frame();
  er_hist.plotOn(er_frame,RooFit::MarkerSize(0.1),RooFit::XErrorSize(0));

  // Fit the PDF to the data and plot on the frame we just created
  er_bukin.fitTo(er_hist);
  er_bukin.plotOn(er_frame,RooFit::LineWidth(1));
  er_frame->SetTitle(energy_resolution->GetTitle());
  er_frame->GetXaxis()->SetTitle("(fatJetE - bPairE) / bPairE");
  er_frame->GetYaxis()->SetTitle("#");
  er_frame->Draw();

  TLegend *leg_energy_resolution = new TLegend(0.7,0.9,0.95,0.95);
  leg_energy_resolution->SetHeader(Form("Bukin width: (%f)", er_sigma.getVal()),"C");
  leg_energy_resolution->Draw();
  c_energy_resolution->SaveAs("energy_resolution.pdf");

  // Turn On Curves
  higgs_pT_with_higgsFatJet->Divide(higgs_pT);
  plot_TH1F(higgs_pT_with_higgsFatJet,"higgs_fatJet_turnOn_Pt.pdf","FatJet TurnOn in higgs pT","pT [GeV]","%","",800,600);
  higgs_pT_with_higgsFatJet_250GeVfatJet->Divide(higgs_pT);
  plot_TH1F(higgs_pT_with_higgsFatJet_250GeVfatJet,"higgs_fatJet_turnOn_Pt_250GeVfatJet.pdf","FatJet > 250GeV TurnOn in higgs pT","pT [GeV]","%","",800,600);
  

  // Study of b-Hadrons falling outside of fatJet
  TCanvas *c_b_outside_fatJet = new TCanvas("b_outside_fatJet","b_outside_fatJet",800,600);
  TH1F *hist_hadron_outside_ptBins = new TH1F("hist_hadron_outside_ptBins","hist_hadron_outside_ptBins",n_bins, Higgs_pt_Bins);
  TH1F *hist_oneMatch_ptBins = new TH1F("hist_oneMatch_ptBins","hist_oneMatch_ptBins",n_bins, Higgs_pt_Bins);
  TH1F *hist_noMatch_in_ptBin = new TH1F("hist_noMatch_in_ptBin","hist_noMatch_in_ptBin",n_bins, Higgs_pt_Bins);
  for(int iterator = 0; iterator < 5; iterator++) {
    hist_hadron_outside_ptBins->SetBinContent(iterator+1,n_hadron_outside_ptBins[iterator]/n_higgs_ptBins[iterator]);
    //hist_hadron_outside_ptBins->SetBinContent(iterator+1,n_hadron_outside_ptBins[iterator]);
    hist_oneMatch_ptBins->SetBinContent(iterator+1,n_oneMatch_ptBins[iterator]/n_higgs_ptBins[iterator]);
    hist_noMatch_in_ptBin->SetBinContent(iterator+1,n_noMatch_in_ptBin[iterator]/n_higgs_ptBins[iterator]);
  }

  hist_oneMatch_ptBins->SetLineColor(1);
  hist_hadron_outside_ptBins->SetLineColor(2);
  hist_noMatch_in_ptBin->SetLineColor(4);

  hist_oneMatch_ptBins->SetMaximum(0.35);
  hist_oneMatch_ptBins->SetMinimum(0);

  //hist_oneMatch_ptBins->Draw("same");
  hist_hadron_outside_ptBins->GetXaxis()->SetTitle("Higgs pT [GeV]");
  hist_hadron_outside_ptBins->GetYaxis()->SetTitle("b outside / total");
  hist_hadron_outside_ptBins->GetYaxis()->SetTitleOffset(1.6);
  
  hist_hadron_outside_ptBins->Draw("same");
  //hist_noMatch_in_ptBin->Draw("same");i

  ATLASLabel(0.5,0.87,"Simulation Internal");

  c_b_outside_fatJet->SaveAs("b_outside_fatJet.pdf");

  // ***********************************************************
  // ******      b-parton smearing 2D resolution plot      *****
  // ***********************************************************
  if(doSmearPlot == 1) {
    int canvas_itr = 1;
    double mass_width = 0;
    
    TCanvas *c_smeared_mass_histograms = new TCanvas("smeared mass histograms","mass histograms",1000,1000);
    TCanvas *c_smeared_angle_histograms = new TCanvas("smeared angle histograms","angle histograms",1000,1000);
    TCanvas *c_smeared_energy_histograms = new TCanvas("smeared energy histograms","energy histograms",1000,1000);
    c_smeared_mass_histograms->Divide(smear_x_bins,smear_y_bins);
    c_smeared_angle_histograms->Divide(smear_x_bins,smear_y_bins);
    c_smeared_energy_histograms->Divide(smear_x_bins,smear_y_bins);
    for(int y_itr = 0; y_itr < smear_y_bins; y_itr++) {
      smeared_bPartons->GetYaxis()->SetBinLabel(y_itr+1,to_string(y_axis_randN_sigma_vector[y_itr]).c_str());
      for(int x_itr = 0; x_itr < smear_x_bins; x_itr++){
        smeared_bPartons->GetXaxis()->SetBinLabel(x_itr+1,to_string(x_axis_randN_sigma_vector[x_itr]).c_str());
        c_smeared_mass_histograms->cd(canvas_itr);
        //smeared_mass_hist_vector[x_itr][y_itr]->Fit("gaus","","",85,165);
        smeared_mass_hist_vector[x_itr][y_itr]->Draw();
        mass_width = (int)(gaus_width_of_FWHM(smeared_mass_hist_vector[x_itr][y_itr]) / 0.01) * 0.01; // Round to nearest one-hundreth
        //mass_width = (int)(bukin_fit_width(smeared_mass_hist_vector[x_itr][y_itr]) / 0.01) * 0.01; // Round to nearest one-hundreth
        //mass_width = (int)(smeared_mass_hist_vector[x_itr][y_itr]->GetRMS() / 0.01) * 0.01; // Round to nearest one-hundreth
        //mass_width = smeared_mass_hist_vector[x_itr][y_itr]->GetFunction("gaus")->GetParameter(2);

        smeared_bPartons->SetBinContent(smeared_bPartons->GetBin(x_itr+1,y_itr+1), mass_width);

        c_smeared_angle_histograms->cd(canvas_itr);
        smeared_angle_hist_vector[x_itr][y_itr]->Draw();
        //smeared_angle_hist_vector[x_itr][y_itr]->Fit("gaus","","",-0.1,0.1);

        c_smeared_energy_histograms->cd(canvas_itr);
        //gStyle->SetStatW(0.4);
        //smeared_energy_hist_vector[x_itr][y_itr]->SetStats(true);
        smeared_energy_hist_vector[x_itr][y_itr]->Draw();
        TLegend *leg_smeared_energy_hist = new TLegend(0.7,0.9,0.95,0.95);
        smeared_energy_hist_vector[x_itr][y_itr]->Fit("gaus","","",smeared_energy_hist_vector[x_itr][y_itr]->GetXaxis()->GetXmin(),smeared_energy_hist_vector[x_itr][y_itr]->GetXaxis()->GetXmax());
        smeared_energy_hist_vector[x_itr][y_itr]->GetFunction("gaus")->SetLineColor(4);
        smeared_energy_hist_vector[x_itr][y_itr]->GetFunction("gaus")->SetLineWidth(1);
        leg_smeared_energy_hist->SetHeader(Form("Gaus Width: (%f)", gaus_width_of_FWHM(smeared_energy_hist_vector[x_itr][y_itr])),"C");
        //leg_smeared_energy_hist->SetHeader(Form("Gaus Mean: (%f)", (int)(smeared_energy_hist_vector[x_itr][y_itr]->GetFunction("gaus")->GetParameter(1) / 0.01) * 0.01),"C");
        leg_smeared_energy_hist->Draw();

        canvas_itr++;
      }
    }
    c_smeared_mass_histograms->SaveAs("smeared_mass_histograms.pdf");
    c_smeared_angle_histograms->SaveAs("smeared_angle_histograms.pdf");
    c_smeared_energy_histograms->SaveAs("smeared_energy_histograms.pdf");

    //gStyle->SetOptStat(0);
    TCanvas *c_smeared_bPartons = new TCanvas("smeared_bPartons","smeared bPartons",800,600);
    smeared_bPartons->GetXaxis()->SetTitle("phi/theta angular resolution [rad]");
    smeared_bPartons->GetYaxis()->SetTitle("energy resolution [original X value]");
    smeared_bPartons->GetZaxis()->SetTitle("mass resolution [GeV]");
    smeared_bPartons->SetStats(false);
    c_smeared_bPartons->SetRightMargin(0.18);
    c_smeared_bPartons->SetLeftMargin(0.18);
    smeared_bPartons->SetMarkerSize(1);
    smeared_bPartons->Draw("COLZ TEXT");
    c_smeared_bPartons->SaveAs("smeared_bPartons.pdf");

    // TGraph Smearing plot
    TCanvas *c_bPair_resolutions = new TCanvas("bPair_resolution","bPair Resolutions",800,600);
    TGraph *mass_res_5 = new TGraph(angular_change_vector[0].size(), &angular_change_vector[0][0], &energy_change_vector[0][0]);
    TGraph *mass_res_10 = new TGraph(angular_change_vector[1].size(), &angular_change_vector[1][0], &energy_change_vector[1][0]);
    TGraph *mass_res_20 = new TGraph(angular_change_vector[2].size(), &angular_change_vector[2][0], &energy_change_vector[2][0]);
    TGraph *mass_res_30 = new TGraph(angular_change_vector[3].size(), &angular_change_vector[3][0], &energy_change_vector[3][0]);
    TGraph *mass_res_40 = new TGraph(angular_change_vector[4].size(), &angular_change_vector[4][0], &energy_change_vector[4][0]);
    TGraph *mass_res_50 = new TGraph(angular_change_vector[5].size(), &angular_change_vector[5][0], &energy_change_vector[5][0]);
    TGraph *mass_res_60 = new TGraph(angular_change_vector[6].size(), &angular_change_vector[6][0], &energy_change_vector[6][0]);
    TGraph *mass_res_70 = new TGraph(angular_change_vector[7].size(), &angular_change_vector[7][0], &energy_change_vector[7][0]);
    TGraph *mass_res_80 = new TGraph(angular_change_vector[8].size(), &angular_change_vector[8][0], &energy_change_vector[8][0]);

    mass_res_80->SetTitle("smeared bPairs colored by [smeared mass - truth mass]");
    mass_res_80->GetXaxis()->SetTitle("change in opening angle [rad]");
    mass_res_80->GetYaxis()->SetTitle("relative change in bPair energy [GeV]");

    mass_res_80->SetMarkerColor(1);mass_res_80->SetMarkerStyle(20);mass_res_80->Draw("AP same");
    mass_res_70->SetMarkerColor(2);mass_res_70->SetMarkerStyle(20);mass_res_70->Draw("P same");
    mass_res_60->SetMarkerColor(3);mass_res_60->SetMarkerStyle(20);mass_res_60->Draw("P same");
    mass_res_50->SetMarkerColor(4);mass_res_50->SetMarkerStyle(20);mass_res_50->Draw("P same");
    mass_res_40->SetMarkerColor(5);mass_res_40->SetMarkerStyle(20);mass_res_40->Draw("P same");
    mass_res_30->SetMarkerColor(6);mass_res_30->SetMarkerStyle(20);mass_res_30->Draw("P same");
    mass_res_20->SetMarkerColor(7);mass_res_20->SetMarkerStyle(20);mass_res_20->Draw("P same");
    mass_res_10->SetMarkerColor(8);mass_res_10->SetMarkerStyle(20);mass_res_10->Draw("P same");
    mass_res_5->SetMarkerColor(9);mass_res_5->SetMarkerStyle(20);mass_res_5->Draw("P same");

    TLegend *leg_bPair_res = new TLegend(0.8,0.6,0.95,0.9);
    leg_bPair_res->AddEntry(mass_res_5,"0 < dM < 5 GeV");
    leg_bPair_res->AddEntry(mass_res_10,"5 < dM< 10 GeV");
    leg_bPair_res->AddEntry(mass_res_20,"10 < dM< 20 GeV");
    leg_bPair_res->AddEntry(mass_res_30,"20 < dM< 30 GeV");
    leg_bPair_res->AddEntry(mass_res_40,"30 < dM< 40 GeV");
    leg_bPair_res->AddEntry(mass_res_50,"40 < dM < 50 GeV");
    leg_bPair_res->AddEntry(mass_res_60,"50 < dM< 60 GeV");
    leg_bPair_res->AddEntry(mass_res_70,"60 < dM < 70 GeV");
    leg_bPair_res->AddEntry(mass_res_80,"70 < dM < 80 GeV");
    leg_bPair_res->Draw("same");

    c_bPair_resolutions->SetRightMargin(0.15);
    c_bPair_resolutions->SetLeftMargin(0.15);

    c_bPair_resolutions->SaveAs("bPair_resolutions.png");

    // ****************************************
    // ******      Sensitivity study      *****
    // ****************************************
    if(doSensitivityStudy == 1) {
      double integral_result = 0;
      double sigma = 0;
      double fatJet_M_bukin_mean = bukin_fit_mean(M_fatJet);
      int lower_mass_bin = 0;;
      int higher_mass_bin = 0;
      int low_x_bin = 4;
      int high_x_bin = 7;
      int low_y_bin = 4;
      int high_y_bin = 7;

      TCanvas *c_mass_resolution_matrix = new TCanvas("mass_resolution_matrix","mass_resolution_matrix",800,600);
      TH2F *mass_resolution_matrix = new TH2F("mass_resolution_matrix","Fatjet Mass Resolution Matrix",3,0,1,3,0,1);

      for(int y_bin_itr = low_y_bin; y_bin_itr < high_y_bin; y_bin_itr++) {
        for(int x_bin_itr = low_x_bin; x_bin_itr < high_x_bin; x_bin_itr++) {
          sigma = smeared_bPartons->GetBinContent(x_bin_itr,y_bin_itr);
          lower_mass_bin = M_fatJet->GetXaxis()->FindBin(fatJet_M_bukin_mean - sigma);
          higher_mass_bin = M_fatJet->GetXaxis()->FindBin(fatJet_M_bukin_mean + sigma);

          mass_resolution_matrix->GetYaxis()->SetBinLabel(y_bin_itr - low_y_bin + 1,to_string(y_axis_randN_sigma_vector[y_bin_itr]).c_str());
          mass_resolution_matrix->GetXaxis()->SetBinLabel(x_bin_itr - low_x_bin + 1,to_string(x_axis_randN_sigma_vector[x_bin_itr]).c_str());

          integral_result = M_fatJet->Integral(lower_mass_bin,higher_mass_bin);
          mass_resolution_matrix->SetBinContent(mass_resolution_matrix->GetBin(x_bin_itr - low_x_bin + 1, y_bin_itr - low_y_bin + 1),integral_result);
        }
      }

      mass_resolution_matrix->SetStats(false);
      mass_resolution_matrix->SetMarkerSize(1);
      mass_resolution_matrix->Draw("COLZ TEXT");
      c_mass_resolution_matrix->SetRightMargin(0.18);
      c_mass_resolution_matrix->SetLeftMargin(0.18);
      c_mass_resolution_matrix->SaveAs("mass_resolution_matrix.pdf");

    } // End of doSensitivityStudy if statement
  } // End of doSmearPlot if statement

  // ************************************************
  // ******      Resolution Plotting Stuff      *****
  // ************************************************

  if(doResolutionPlots == 1) {
    //  Fit function declaration:  void fit_histo(vector<TH1F*> hist_vector, vector<TH1F*> bins_hist_vector, vector<double> &var_mean_vector, vector<double> &mass_mean_vector, vector<double> &mass_sigma_vector, vector<double> &mass_sigma_error_vector, vector<double> &mass_error_low_vector, vector<double> &mass_error_high_vector, vector<double> &bin_low_edge_vector, vector<double> &bin_high_edge_vector, double xLow, double xHigh) 
    // Generate Graph of Resolution for fatJet mass in parton opening Angle pT bins
    cout << "***********************************************************************************" << endl;
    cout << "**** Generate Graph of Resolution for fatJet mass in parton opening angle bins ****" << endl;
    cout << "***********************************************************************************" << endl;
    fit_histo(deltaE_bPartons_mass_res_hist_vector, deltaE_bPartons_bins_hist_vector, deltaE_bPartons_mean, deltaE_bPartons_bin_fatJet_mass_mean, deltaE_bPartons_bin_fatJet_mass_sigma, deltaE_bPartons_bin_fatJet_mass_sigma_error, deltaE_bPartons_error_low, deltaE_bPartons_error_high, deltaE_bPartons_bin_low_edge, deltaE_bPartons_bin_high_edge, 65, 190, "deltaE_bPartons_fatJet_Mass_plots.pdf");

    // Generate Graph of Resolution for fatJet mass in parton opening Angle pT bins
    cout << "***********************************************************************************" << endl;
    cout << "**** Generate Graph of Resolution for fatJet mass in parton opening angle bins ****" << endl;
    cout << "***********************************************************************************" << endl;
    fit_histo(partonAngle_mass_res_hist_vector, partonAngle_bins_hist_vector, partonAngle_mean, partonAngle_bin_fatJet_mass_mean, partonAngle_bin_fatJet_mass_sigma, partonAngle_bin_fatJet_mass_sigma_error, partonAngle_error_low, partonAngle_error_high, partonAngle_bin_low_edge, partonAngle_bin_high_edge, 65, 190, "partonAngle_fatJet_Mass_plots.pdf");
    
    // Generate Graph of Resolution for fatJet mass in fatJet pT bins
    cout << "************************************************************************" << endl;
    cout << "**** Generate Graph of Resolution for fatJet mass in fatJet pT bins ****" << endl;
    cout << "************************************************************************" << endl;
    fit_histo(fatJet_pT_mass_res_hist_vector, fatJet_pT_bins_hist_vector, fatJet_pT_mean, fatJet_bin_fatJet_mass_mean_pT, fatJet_bin_fatJet_mass_sigma_pT, fatJet_bin_fatJet_mass_sigma_error_pT, fatJet_pT_error_low, fatJet_pT_error_high, fatJet_pT_bin_low_edge, fatJet_pT_bin_high_edge, 65, 190, "fatJet_pt_fatJet_Mass_plots.pdf");

    // Generate Graph of Resolution for fatJet mass in fatJet E bins
    cout << "************************************************************************" << endl;
    cout << "**** Generate Graph of Resolution for fatJet mass in fatJet E bins ****" << endl;
    cout << "************************************************************************" << endl;
    fit_histo(fatJet_E_mass_res_hist_vector, fatJet_E_bins_hist_vector, fatJet_E_mean, fatJet_bin_fatJet_mass_mean_E, fatJet_bin_fatJet_mass_sigma_E, fatJet_bin_fatJet_mass_sigma_error_E, fatJet_E_error_low, fatJet_E_error_high, fatJet_E_bin_low_edge, fatJet_E_bin_high_edge, 65, 190, "fatJet_energy_fatJet_Mass_plots.pdf");

    // Generate Graph of Resolution for fatJet mass in Higgs pT bins
    cout << "***********************************************************************" << endl;
    cout << "**** Generate Graph of Resolution for fatJet mass in Higgs pT bins ****" << endl;
    cout << "***********************************************************************" << endl;
    fit_histo(higgs_pT_mass_res_hist_vector, higgs_pT_bins_hist_vector, higgs_pT_mean, higgs_bin_fatJet_mass_mean_pT, higgs_bin_fatJet_mass_sigma_pT, higgs_bin_fatJet_mass_sigma_error_pT, higgs_pT_error_low, higgs_pT_error_high, higgs_pT_bin_low_edge, higgs_pT_bin_high_edge, 65, 190, "higgs_pT_fatJet_Mass_plots.pdf");

    // Generate Graph of Resolution for fatJet mass in Higgs E bins
    cout << "***********************************************************************" << endl;
    cout << "**** Generate Graph of Resolution for fatJet mass in Higgs E bins ****" << endl;
    cout << "***********************************************************************" << endl;
    fit_histo(higgs_E_mass_res_hist_vector, higgs_E_bins_hist_vector, higgs_E_mean, higgs_bin_fatJet_mass_mean_E, higgs_bin_fatJet_mass_sigma_E, higgs_bin_fatJet_mass_sigma_error_E, higgs_E_error_low, higgs_E_error_high, higgs_E_bin_low_edge, higgs_E_bin_high_edge, 65, 190, "higgs_energy_fatJet_Mass_plots.pdf");
    
    // *****  Resolution Plotting mess *****
    // Resolution plots in parton opening angle bins

    TCanvas *c_mass_res_deltaE_bPartons = new TCanvas("fatJet Mass Resolution in bins of deltaE of the b Partons","fatJet Mass Resolution in bins of deltaE of the b Partons",800,600);
    TGraphAsymmErrors *fatJet_mass_res_deltaE_bPartons_bins = new TGraphAsymmErrors(deltaE_bPartons_mass_res_hist_vector.size(), &deltaE_bPartons_mean[0], &deltaE_bPartons_bin_fatJet_mass_sigma[0], &deltaE_bPartons_error_low[0], &deltaE_bPartons_error_high[0], &deltaE_bPartons_bin_fatJet_mass_sigma_error[0], &deltaE_bPartons_bin_fatJet_mass_sigma_error[0]);
    fatJet_mass_res_deltaE_bPartons_bins->SetTitle("Higgs Associated Fatjet Mass Resolution in exclusive deltaE of the b partons bins 0 < angle < 500");
    fatJet_mass_res_deltaE_bPartons_bins->SetMinimum(0);
    fatJet_mass_res_deltaE_bPartons_bins->SetMaximum(20);
    fatJet_mass_res_deltaE_bPartons_bins->Draw();
    c_mass_res_deltaE_bPartons->SaveAs("fatJet_massResolution_deltaE_bPartons_bins.pdf");

    TCanvas *c_mass_res_partonAngle = new TCanvas("fatJet Mass Resolution in parton opening angle bins","fatJet Mass Resoltuion in parton opening angle bins",800,600);
    TGraphAsymmErrors *fatJet_mass_res_partonAngle_bins = new TGraphAsymmErrors(partonAngle_mass_res_hist_vector.size(), &partonAngle_mean[0], &partonAngle_bin_fatJet_mass_sigma[0], &partonAngle_error_low[0], &partonAngle_error_high[0], &partonAngle_bin_fatJet_mass_sigma_error[0], &partonAngle_bin_fatJet_mass_sigma_error[0]);
    fatJet_mass_res_partonAngle_bins->SetTitle("Higgs Associated Fatjet Mass Resolution in exclusive parton opening angle bins 0 < angle < 3.2");
    fatJet_mass_res_partonAngle_bins->SetMinimum(0);
    fatJet_mass_res_partonAngle_bins->SetMaximum(20);
    fatJet_mass_res_partonAngle_bins->Draw();
    c_mass_res_partonAngle->SaveAs("fatJet_massResolution_partonAngle_bins.pdf");
    
    // Resolution plots in fatJet bins
    TCanvas *c_mass_res_fatJet_pT = new TCanvas("fatJet Mass Resolution in fatJet pT bins","fatJet Mass Resoltuion in fatJet pT bins",800,600);
    TGraphAsymmErrors *fatJet_mass_res_pT_bins = new TGraphAsymmErrors(fatJet_pT_mass_res_hist_vector.size(), &fatJet_pT_mean[0], &fatJet_bin_fatJet_mass_sigma_pT[0], &fatJet_pT_error_low[0], &fatJet_pT_error_high[0], &fatJet_bin_fatJet_mass_sigma_error_pT[0], &fatJet_bin_fatJet_mass_sigma_error_pT[0]);
    fatJet_mass_res_pT_bins->SetTitle("Higgs Associated Fatjet Mass Resolution in exclusive fatJet pT bins 300 < pT < 1500");
    fatJet_mass_res_pT_bins->SetMinimum(0);
    fatJet_mass_res_pT_bins->SetMaximum(20);
    fatJet_mass_res_pT_bins->Draw();
    c_mass_res_fatJet_pT->SaveAs("fatJet_massResolution_fatJet_pT_bins.pdf");

    TCanvas *c_mass_res_fatJet_E = new TCanvas("fatJet Mass Resolution in fatJet E bins","fatJet Mass Resoltuion in fatJet E bins",800,600);
    TGraphAsymmErrors *fatJet_mass_res_E_bins = new TGraphAsymmErrors(fatJet_E_mass_res_hist_vector.size(), &fatJet_E_mean[0], &fatJet_bin_fatJet_mass_sigma_E[0], &fatJet_E_error_low[0], &fatJet_E_error_high[0], &fatJet_bin_fatJet_mass_sigma_error_E[0], &fatJet_bin_fatJet_mass_sigma_error_E[0]);
    fatJet_mass_res_E_bins->SetTitle("Higgs Associated Fatjet Mass Resolution in exclusive fatJet E bins 100 < E < 2000");
    fatJet_mass_res_E_bins->SetMinimum(0);
    fatJet_mass_res_E_bins->SetMaximum(20);
    fatJet_mass_res_E_bins->Draw();
    c_mass_res_fatJet_E->SaveAs("fatJet_massResolution_fatJet_E_bins.pdf");

    // Resolution plots in higgs bins
    TCanvas *c_mass_res_higgs_pT = new TCanvas("fatJet Mass Resolution in Higgs pT bins","fatJet Mass Resoltuion in Higgs pT bins",800,600);
    TGraphAsymmErrors *higgs_mass_res_pT_bins = new TGraphAsymmErrors(higgs_pT_mass_res_hist_vector.size(), &higgs_pT_mean[0], &higgs_bin_fatJet_mass_sigma_pT[0], &higgs_pT_error_low[0], &higgs_pT_error_high[0], &higgs_bin_fatJet_mass_sigma_error_pT[0], &higgs_bin_fatJet_mass_sigma_error_pT[0]);
    higgs_mass_res_pT_bins->SetTitle("Higgs Associated Fatjet Mass Resolution in exclusive Truth Higgs pT bins 300 < pT < 1500");
    higgs_mass_res_pT_bins->SetMinimum(0);
    higgs_mass_res_pT_bins->SetMaximum(20);
    higgs_mass_res_pT_bins->Draw();
    c_mass_res_higgs_pT->SaveAs("Higgs_massResolution_Higgs_pT_bins.pdf");

    TCanvas *c_mass_res_higgs_E = new TCanvas("fatJet Mass Resolution in Higgs E bins","fatJet Mass Resoltuion in Higgs E bins",800,600);
    TGraphAsymmErrors *higgs_mass_res_E_bins = new TGraphAsymmErrors(higgs_E_mass_res_hist_vector.size(), &higgs_E_mean[0], &higgs_bin_fatJet_mass_sigma_E[0], &higgs_E_error_low[0], &higgs_E_error_high[0], &higgs_bin_fatJet_mass_sigma_error_E[0], &higgs_bin_fatJet_mass_sigma_error_E[0]);
    higgs_mass_res_E_bins->SetTitle("Higgs Associated Fatjet Mass Resolution in exclusive Truth Higgs E bins 100 < E < 2000");
    higgs_mass_res_E_bins->SetMinimum(0);
    higgs_mass_res_E_bins->SetMaximum(20);
    higgs_mass_res_E_bins->Draw();
    c_mass_res_higgs_E->SaveAs("Higgs_massResolution_Higgs_E_bins.pdf");

    // Function declaration: void make_bin_plots(vector<TH1F*> hist_vector, string output_filename) 
    // Resolution Plots in parton opening angle bins
    //make_bin_plots(deltaE_bPartons_mass_res_hist_vector,"deltaE_bPartons_fatJet_Mass_plots.pdf");
    make_bin_plots(deltaE_bPartons_bins_hist_vector,"deltaE_bPartons_bin_plots.pdf");

    //make_bin_plots(partonAngle_mass_res_hist_vector,"partonAngle_fatJet_Mass_plots.pdf");
    make_bin_plots(partonAngle_bins_hist_vector,"partonAngle_bin_plots.pdf");

    // Resolution Plots in FatJet Pt bins
    //make_bin_plots(fatJet_pT_mass_res_hist_vector,"fatJet_pt_fatJet_Mass_plots.pdf");
    make_bin_plots(fatJet_pT_bins_hist_vector,"fatJet_pT_bin_plots.pdf");

    // Resolution Plots in FatJet E bins
    //make_bin_plots(fatJet_E_mass_res_hist_vector,"fatJet_energy_fatJet_Mass_plots.pdf");
    make_bin_plots(fatJet_E_bins_hist_vector,"fatJet_energy_bin_plots.pdf");

    // Resolution Plots in Higgs pT bins
    //make_bin_plots(higgs_pT_mass_res_hist_vector,"higgs_pT_fatJet_Mass_plots.pdf");
    make_bin_plots(higgs_pT_bins_hist_vector,"higgs_pT_bin_plots.pdf");

    // Resolution Plots in Higgs E bins
    //make_bin_plots(higgs_E_mass_res_hist_vector,"higgs_energy_fatJet_Mass_plots.pdf");
    make_bin_plots(higgs_E_bins_hist_vector,"higgs_energy_bin_plots.pdf");
  }
  cout << counter << endl;
  output_file.Write();
}
