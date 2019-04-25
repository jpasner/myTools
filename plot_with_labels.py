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
  print "Format: python plot.py file1.root file2.root (etc...)"
  # Generally this looks like: python plot.py folder/*.root
  exit()

hist_handle

for i,iFile in enumerate(my_input_file_list):
  my_input_file = ROOT.TFile(iFile,"read")
  if "ttbar" in iFile:
    print "Found ttbar"
    my_input_file.GetObject("hbbisr_GhostVR30Rmax4Rmin02TrackJet_fj0pt480_fj1pt250_tjet2_ptsort0_hpt480_nbtag1in2ext_MV2c10_FixedCutBEff_77/muon0_dPhi",hist_handle)
    ttbar_dPhi_histogram = copy.deepcopy(hist_handle)
    my_input_file.GetObject("hbbisr_GhostVR30Rmax4Rmin02TrackJet_fj0pt480_fj1pt250_tjet2_ptsort0_hpt480_nbtag2in2ext_MV2c10_FixedCutBEff_77/muon0_pt",hist_handle)
    ttbar_pT_histogram = copy.deepcopy(hist_handle)
  elif "qcd" in iFile:
    print "Found qcd"
    my_input_file.GetObject("hbbisr_GhostVR30Rmax4Rmin02TrackJet_fj0pt480_fj1pt250_tjet2_ptsort0_hpt480_nbtag1in2ext_MV2c10_FixedCutBEff_77/muon0_dPhi",hist_handle)
    qcd_dPhi_histogram = copy.deepcopy(hist_handle)
    my_input_file.GetObject("hbbisr_GhostVR30Rmax4Rmin02TrackJet_fj0pt480_fj1pt250_tjet2_ptsort0_hpt480_nbtag2in2ext_MV2c10_FixedCutBEff_77/muon0_pt",hist_handle)
    qcd_pT_histogram = copy.deepcopy(hist_handle)

canvas = ROOT.TCanvas("Canvas for Histogram Plotting","Canvas for Histogram Plotting",800,600)



# ********** dPhi plots **********
#ttbar_dPhi_hist ogram.Scale(1.0/ttbar_histogram.Integral())
ttbar_dPhi_histogram.Rebin(3)
ttbar_dPhi_histogram.GetXaxis().SetTitle("dPhi(leading muon - Hcand)")

#qcd_dPhi_histogram.Scale(1.0/qcd_histogram.Integral())
qcd_dPhi_histogram.SetLineColor(2)
qcd_dPhi_histogram.Rebin(3)

ttbar_dPhi_histogram.Draw()
qcd_dPhi_histogram.Draw("same")

# ATLAS Label
text=[]
suffix = 'Internal'
lumi = 80.7
xlabel = .4
ylabel = .8

text.append('#font[72]{ATLAS} %s'%suffix)
text.append("Simulation, #sqrt{s} = 13 TeV")
text.append("1 tag Region, %s fb^{-1}"%lumi)

latext=None
for i in range(len(text)):
        if latext==None: latext=text[i]
        else: latext='#splitline{%s}{%s}'%(latext,text[i])

Tl = ROOT.TLatex()
Tl.SetNDC()

leg = ROOT.TLegend(.4,.6,.5,.7)
leg.AddEntry(ttbar_dPhi_histogram, "t#bart")
leg.AddEntry(qcd_dPhi_histogram, "QCD")
leg.SetBorderSize(0)

Tl.DrawLatex(xlabel, ylabel, latext)
leg.Draw()

canvas.SaveAs("dPhi.pdf")


# ********** pT plots **********
#ttbar_pT_hist ogram.Scale(1.0/ttbar_histogram.Integral())

#qcd_pT_histogram.Scale(1.0/qcd_histogram.Integral())
qcd_pT_histogram.SetLineColor(2)
qcd_pT_histogram.GetYaxis().SetTitle("Events [/ 5 GeV]")

qcd_pT_histogram.Draw()
ttbar_pT_histogram.Draw("same")

# ATLAS Label
text=[]
suffix = 'Internal'
lumi = 80.7
xlabel = .4
ylabel = .8

text.append('#font[72]{ATLAS} %s'%suffix)
text.append("Simulation, #sqrt{s} = 13 TeV")
text.append("Signal Region, %s fb^{-1}"%lumi)

latext=None
for i in range(len(text)):
        if latext==None: latext=text[i]
        else: latext='#splitline{%s}{%s}'%(latext,text[i])

Tl = ROOT.TLatex()
Tl.SetNDC()

leg = ROOT.TLegend(.4,.6,.5,.7)
leg.AddEntry(ttbar_dPhi_histogram, "t#bar{t}")
leg.AddEntry(qcd_pT_histogram, "QCD")
leg.SetBorderSize(0)

Tl.DrawLatex(xlabel, ylabel, latext)
leg.Draw()

canvas.SaveAs("pT.pdf")
# *********************************
