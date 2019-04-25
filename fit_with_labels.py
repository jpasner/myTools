from ROOT import *
import ROOT
import os
import re
import sys
import math
import copy
# Run in Batch mode to suppress plots being displayed to screen (which is slow for long distance servers)
ROOT.gROOT.SetBatch(ROOT.kTRUE)

# Imports ATLAS style for plotting
ROOT.gROOT.LoadMacro("/global/homes/j/jpasner/atlasstyle/AtlasUtils.C")
ROOT.gROOT.LoadMacro("/global/homes/j/jpasner/atlasstyle/AtlasStyle.C")
ROOT.gROOT.LoadMacro("/global/homes/j/jpasner/atlasstyle/AtlasLabels.C")
SetAtlasStyle()

input_file_list = [sys.argv] # Get files from command line
my_input_file_list = copy.deepcopy(input_file_list[0]) # Make our own copy, leave original alone
my_input_file_list.pop(0) # Get rid of first object which should be the script itself

hist_handle = ROOT.TH1F()

if len(sys.argv) < 2:
  print "Please incidcate physics process files containing histograms to be plotted"
  print "Format: python fit_with_labels.py file1.root file2.root (etc...)"
  # Generally this looks like: python fit_with_labels.py folder/*.root
  exit()

hist_handle

for i,iFile in enumerate(my_input_file_list):
  my_input_file = ROOT.TFile(iFile,"read")
  print "Found likelihoodFunction"
  hist_handle = my_input_file.Get("likelihoodFunction")
  likelihoodFunction_hist = copy.deepcopy(hist_handle)

  hist_handle = my_input_file.Get("h_sr")
  h_sr = copy.deepcopy(hist_handle)

canvas = ROOT.TCanvas("Canvas for Histogram Plotting","Canvas for Histogram Plotting",800,600)

# ********** likelihoodFunction bukin fit code**********

# Fit distribution using a bukin and return its width
SIGNAL_normalisation = RooRealVar("SIGNAL_normalisation","SIGNAL_normalisation",likelihoodFunction_hist.GetXaxis().GetXmin() ,likelihoodFunction_hist.GetXaxis().GetXmax());

# Have to convert RooRealVar into RooArgList since python doesn't allow implicit converstions
argList = RooArgList(SIGNAL_normalisation)

# Declare histogram to be fit along and a variable (smeared_Mass) to be filled with fit results
hist = RooDataHist("hist","histogram of smeared Mass",argList,likelihoodFunction_hist);

# Declare bukin variables and function. Don't forget to give variables ranges or they will be considered constant in the fit
# 2Tag
peak_position = RooRealVar("peak_position","location of peak on x axis",130,100,200);
sigma = RooRealVar("sigma","resolution of Mass defined as: FWHM divided by 2*sqrt(2*log(2))=2.35",50,0,100);
peak_asymmetry = RooRealVar("peak_asymmetry","peak asymmetry parameter",.14,-1,1);
left_tail = RooRealVar("left_tail","parameter of the left tail",-.5,-10,10);
right_tail = RooRealVar("right_tail","parameter of the right tail",-.2,-10,10);
bukin = RooBukinPdf("my bukin name","my bukin title",SIGNAL_normalisation,peak_position,sigma,peak_asymmetry,left_tail,right_tail);


# 1Tag
#peak_position = RooRealVar("peak_position","location of peak on x axis",1100,900,1200);
#sigma = RooRealVar("sigma","resolution of Mass defined as: FWHM divided by 2*sqrt(2*log(2))=2.35",130,50,250);
#peak_asymmetry = RooRealVar("peak_asymmetry","peak asymmetry parameter",.14,-0.3,0.3);
#left_tail = RooRealVar("left_tail","parameter of the left tail",-.5,-1.5,0);
#right_tail = RooRealVar("right_tail","parameter of the right tail",-.2,-1,1);
#bukin = RooBukinPdf("my bukin name","my bukin title",SIGNAL_normalisation,peak_position,sigma,peak_asymmetry,left_tail,right_tail);

# 0Tag
#peak_position = RooRealVar("peak_position","location of peak on x axis",1900,1500,3000);
#sigma = RooRealVar("sigma","resolution of Mass defined as: FWHM divided by 2*sqrt(2*log(2))=2.35",200,150,350);
#peak_asymmetry = RooRealVar("peak_asymmetry","peak asymmetry parameter",.14,-1,1);
#left_tail = RooRealVar("left_tail","parameter of the left tail",-.5,-2,2);
#right_tail = RooRealVar("right_tail","parameter of the right tail",-.2,-2,2);
#bukin = RooBukinPdf("my bukin name","my bukin title",SIGNAL_normalisation,peak_position,sigma,peak_asymmetry,left_tail,right_tail);

# Create frame for RooFit plot result to be plotted on
frame = SIGNAL_normalisation.frame();
hist.plotOn(frame,ROOT.RooFit.MarkerSize(0.1),ROOT.RooFit.XErrorSize(0));

# Fit the PDF to the data and plot on the frame we just created
bukin.fitTo(hist);
bukin.plotOn(frame,ROOT.RooFit.LineWidth(1));
frame.SetTitle(likelihoodFunction_hist.GetTitle());
likelihoodFunction_hist.Draw()
frame.Draw("same");

# Store results from fit
#return(sigma.getVal()); # return mass resolution result from fit


# ATLAS Label
text=[]
suffix = 'Internal'
lumi = 80.7
xlabel = .55
ylabel = .8

text.append('#font[72]{ATLAS} %s'%suffix)
text.append("Simulation, #sqrt{s} = 13 TeV")
text.append("ttbar 2tag Region, %s fb^{-1}"%lumi)

latext=None
for i in range(len(text)):
        if latext==None: latext=text[i]
        else: latext='#splitline{%s}{%s}'%(latext,text[i])

Tl = ROOT.TLatex()
Tl.SetNDC()

k_factor = peak_position.getVal()/h_sr.Integral()
k_factor_sigma = sigma.getVal()/h_sr.Integral()

leg = ROOT.TLegend(.6,.6,.9,.7)


leg.AddEntry(likelihoodFunction_hist, "likelihoodFunction")
leg.AddEntry(frame, "Bukin Fit:  " + str(round(peak_position.getVal(),2)) + "#pm" + str(round(sigma.getVal(),2)))
leg.AddEntry(0, "nEvents Pre-fit in this range: " + str(round(h_sr.Integral(),2)), "")
leg.AddEntry(0, "k factor: " + str(round(k_factor,2)) + "#pm" + str(round(k_factor_sigma,2)), "")  
leg.SetBorderSize(0)

print h_sr.Integral()

Tl.DrawLatex(xlabel, ylabel, latext)
leg.Draw()

canvas.SaveAs("likelihoodFunction_bukinFit_result.pdf")
# *********************************
