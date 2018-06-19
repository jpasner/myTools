import sys
import ROOT
import copy
import re
from array import *

#   INFO:
#     Edit list_of_distributions to decide which processes to creat the ratio for
#     Script expects you to provide an electron and muon input file with the format python ratio_plotter.py file1.root file2.root

ROOT.gROOT.LoadMacro("/afs/cern.ch/user/j/jpasner/atlasstyle/AtlasStyle.C")

def ratioplot(h1 ,h2 ,distribution):
  "Takes care of creating canvas with both histograms on top and ratio on bottom"
  # Define the Canvas
  c =  ROOT.TCanvas("c", "canvas", 800, 800)

  # Upper plot will be in pad1
  pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
  pad1.SetBottomMargin(0) # Upper and lower plot are joined
  pad1.SetGridx()         # Vertical grid
  pad1.Draw()             # Draw the upper pad: pad1
  pad1.cd()               # pad1 becomes the current pad
  h1.SetStats(0)          # No statistics on upper plot
  h1.Draw()               # Draw h1
  h2.Draw("same")         # Draw h2 on top of h1

  # Do not draw the Y axis label on the upper plot and redraw a small
  # axis instead, in order to avoid the first label (0) to be clipped.
  #h1.GetYaxis().SetLabelSize(0.)
  #axis = ROOT.TGaxis( -5, 20, -5, 220, 20,220,510,"")
  #axis.SetLabelFont(43) # Absolute font size in pixel (precision 3)
  #axis.SetLabelSize(15)
  #axis.Draw()

  # Chiara reccommendation for h1.GetYaxis().SetRange() ??

  # lower plot will be in pad
  c.cd()          # Go back to the main canvas before defining pad2
  pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
  pad2.SetTopMargin(0)
  pad2.SetBottomMargin(0.2)
  pad2.SetGridx() # vertical grid
  pad2.Draw()
  pad2.cd()       # pad2 becomes the current pad

  # Define the ratio plot
  h3 = h1.Clone("h3")
  h3.SetLineColor(ROOT.kBlack)
  h3.SetMinimum(0)  # Define Y ..
  h3.SetMaximum(1.5) # .. range
  h3.Sumw2()
  h3.SetStats(0)      # No statistics on lower plot
  h3.Divide(h2)
  h3.SetMarkerStyle(21)
  h3.Draw("ep")       # Draw the ratio plot

  # h1 settings
  h1.SetLineColor(ROOT.kBlue+1)
  h1.SetLineWidth(2)

  # Y axis h1 plot settings
  h1.GetYaxis().SetTitleSize(20)
  h1.GetYaxis().SetTitleFont(43)
  h1.GetYaxis().SetTitleOffset(1.55)

  # h2 settings
  h2.SetLineColor(ROOT.kRed)
  h2.SetLineWidth(2)

  # Ratio plot (h3) settings
  h3.SetTitle("") # Remove the ratio title

  # Y axis ratio plot settings
  ratio_name = h1.GetName() + '/' + h2.GetName()
  h3.GetYaxis().SetTitle(ratio_name)
  h3.GetYaxis().SetNdivisions(505)
  h3.GetYaxis().SetTitleSize(20)
  h3.GetYaxis().SetTitleFont(43)
  h3.GetYaxis().SetTitleOffset(1.55)
  h3.GetYaxis().SetLabelFont(43) # Absolute font size in pixel (precision 3)
  h3.GetYaxis().SetLabelSize(15)

  # X axis ratio plot settings
  h3.GetXaxis().SetTitleSize(20)
  h3.GetXaxis().SetTitleFont(43)
  h3.GetXaxis().SetTitleOffset(4.)
  h3.GetXaxis().SetLabelFont(43) # Absolute font size in pixel (precision 3)
  h3.GetXaxis().SetLabelSize(15)

  c.SaveAs(distribution + ".pdf")


################################################################################################

input_file_list = [sys.argv] # Get files from command line

if len(sys.argv) is not 3:
  print "Please provide the electron and muon input files you want the raio of using the format: python ratio_plotter.py file1.root file2.root"
  exit()

my_input_file_list = copy.deepcopy(input_file_list[0]) # Personal copy
my_input_file_list.pop(0) # fixes bug where script would think the first file was the script itself 

#list_of_distributions = ['data','Zbb']
#list_of_distributions = ['data','Z1l','Z1c','Z1b']
list_of_distributions = ['data','W1l','W1c','W1b']

for j,iDist in enumerate(list_of_distributions):
  print "Generating ratio plot for " + iDist

  for i,iFile in enumerate(my_input_file_list):
    my_input_file = ROOT.TFile(iFile,"read")

    # Check that file exists, exit if it doesn't
    if 'El' in iFile and my_input_file.IsOpen():
      electron_histogram = copy.deepcopy(my_input_file.Get(iDist))
      electron_histogram.SetName('Electron')
    elif 'Mu' in iFile and my_input_file.IsOpen():
      muon_histogram = copy.deepcopy(my_input_file.Get(iDist))
      muon_histogram.SetName('Muon')
    else:
      print "Couldn't open this file: "
      my_input_file.Print()
      break

  ratioplot(muon_histogram, electron_histogram, iDist)
