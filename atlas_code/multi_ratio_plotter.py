import sys
import ROOT
import copy
import re
from array import *

#   INFO:
#     Edit list_of_distributions to decide which processes to creat the ratio for
#     Script expects you to provide an electron and muon input file with the format python ratio_plotter.py file1.root file2.root

ROOT.gROOT.LoadMacro("/afs/cern.ch/user/j/jpasner/atlasstyle/AtlasStyle.C")

def ratioplot( h_baseline, h_array, color_array, distribution):
  "Takes care of creating canvas with both histograms on top and ratio on bottom"
  # Define the Canvas
  c =  ROOT.TCanvas("c", "canvas", 800, 800)

  # Upper plot will be in pad1
  pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
  pad1.SetBottomMargin(0) # Upper and lower plot are joined
  pad1.SetGridx()         # Vertical grid
  pad1.Draw()             # Draw the upper pad: pad1
  pad1.cd()               # pad1 becomes the current pad
  h_baseline.SetStats(0)          # No statistics on upper plot
  h_baseline.Draw()               # Draw h_baseline
  for i,histo in enumerate(h_array):
    histo.Draw("same")         # Draw h2 on top of h_baseline

  # Do not draw the Y axis label on the upper plot and redraw a small
  # axis instead, in order to avoid the first label (0) to be clipped.
  #h_baseline.GetYaxis().SetLabelSize(0.)
  #axis = ROOT.TGaxis( -5, 20, -5, 220, 20,220,510,"")
  #axis.SetLabelFont(43) # Absolute font size in pixel (precision 3)
  #axis.SetLabelSize(15)
  #axis.Draw()

  # Chiara reccommendation for h_baseline.GetYaxis().SetRange() ??

  # lower plot will be in pad
  c.cd()          # Go back to the main canvas before defining pad2
  pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
  pad2.SetTopMargin(0)
  pad2.SetBottomMargin(0.2)
  pad2.SetGridx() # vertical grid
  pad2.Draw()
  pad2.cd()       # pad2 becomes the current pad

  h_ratio = []

  # Define the ratio plot
  for i,histo in enumerate(h_array):
    h_ratio.append(h_baseline.Clone("h_ratio"))
    h_ratio[i].SetMinimum(-2.0)  # Define Y ..
    h_ratio[i].SetMaximum(10.0) # .. range
    h_ratio[i].Sumw2()
    h_ratio[i].SetStats(0)      # No statistics on lower plot
    h_ratio[i].SetLineColor(color_array[i])
    h_ratio[i].Divide(histo)
    #h_ratio[i].SetMarkerStyle(21)
    h_ratio[i].SetTitle("") # Remove the ratio title
    h_ratio[i].Draw("same HIST")       # Draw the ratio plot
    h_array[i].SetLineColor(color_array[i])
    h_array[i].SetLineWidth(2)

  # h_baseline settings
  h_baseline.SetLineColor(1)
  h_baseline.SetLineWidth(2)

  # Y axis h_baseline plot settings
  h_baseline.GetYaxis().SetTitleSize(20)
  h_baseline.GetYaxis().SetTitleFont(43)
  h_baseline.GetYaxis().SetTitleOffset(1.55)

  # h2 settings
  #h2.SetLineColor(ROOT.kRed)
  #h2.SetLineWidth(2)

  # Ratio plot (h_ratio) settings

  # Y axis ratio plot settings
  ratio_name = h_baseline.GetName() + '/' + 'other'
  
  for i,histo in enumerate(h_array):
    h_ratio[i].GetYaxis().SetTitle(ratio_name)
    h_ratio[i].GetYaxis().SetNdivisions(505)
    h_ratio[i].GetYaxis().SetTitleSize(20)
    h_ratio[i].GetYaxis().SetTitleFont(43)
    h_ratio[i].GetYaxis().SetTitleOffset(1.55)
    h_ratio[i].GetYaxis().SetLabelFont(43) # Absolute font size in pixel (precision 3)
    h_ratio[i].GetYaxis().SetLabelSize(15)

    # X axis ratio plot settings
    h_ratio[i].GetXaxis().SetTitleSize(20)
    h_ratio[i].GetXaxis().SetTitleFont(43)
    h_ratio[i].GetXaxis().SetTitleOffset(4.)
    h_ratio[i].GetXaxis().SetLabelFont(43) # Absolute font size in pixel (precision 3)
    h_ratio[i].GetXaxis().SetLabelSize(15)

  c.SaveAs("sigFatJet_ratio.pdf")


################################################################################################

input_file_list = [sys.argv] # Get files from command line

if len(sys.argv) is not 2:
  print "Please provide the electron and muon input files you want the raio of using the format: python ratio_plotter.py file1.root file2.root"
  exit()

my_input_file_list = copy.deepcopy(input_file_list[0]) # Personal copy
my_input_file_list.pop(0) # fixes bug where script would think the first file was the script itself 
print my_input_file_list[0]
my_input_file = ROOT.TFile(my_input_file_list[0],"read")

list_of_distributions = ['corrected_passCut_sigFatJet_mass','corrected_noPassCut_sigFatJet_mass','passCut_sigFatJet_mass','noPassCut_sigFatJet_mass']
color_array = [2,4,6,8]

histogram_array = []

baseline_histogram = copy.deepcopy(my_input_file.Get("sigFatJet_mass"))
baseline_histogram.SetName('sigFatJet_mass')

for j,iDist in enumerate(list_of_distributions):
  print "Generating ratio plot for " + iDist

  histogram_array.append(copy.deepcopy(my_input_file.Get(iDist)))
  print histogram_array[j]
  histogram_array[j].SetName(iDist)

ratioplot(baseline_histogram, histogram_array, color_array, iDist)
