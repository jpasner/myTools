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

#
# Script to stack HbbISR histograms for different physics processes stored in different files
#

def choose_color(sample_name ,histo):
  "Give correct colors to samples to match Francesco's colors"
  if "hist-ttbar" in sample_name:
    histo.SetFillColor(ROOT.kGreen+2) 
  elif "hist-qcd" in sample_name:
    #histo.SetFillColor(ROOT.UCYellow) This is this color Francesco used but I can't get it to work
    histo.SetFillColor(ROOT.kOrange+1) 
  elif "hist-sig" in sample_name:
    histo.SetFillColor(ROOT.kBlack) 
  elif "hist-singletop" in sample_name:
    histo.SetFillColor(ROOT.kPink+2) 
  elif "hist-w" in sample_name:
    histo.SetFillColor(ROOT.kBlue) 
  elif "hist-z" in sample_name:
    histo.SetFillColor(ROOT.kRed) 
  else:
    print "********************************************************************************"
    print "Input files didn't match expected file name formatting for choosing histo colors"
    print "********************************************************************************"

def ratio_plot(data_histo ,mc_sum_histo ,mc_stack_histo ,legend ,distribution):
  "Takes care of creating canvas with both histograms on top and ratio on bottom"
  # Define the Canvas
  c =  ROOT.TCanvas("c", "canvas", 800, 800)

  # Upper plot will be in pad1
  pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
  pad1.SetBottomMargin(0) # Upper and lower plot are joined
  pad1.SetGridx()         # Vertical grid
  pad1.Draw()             # Draw the upper pad: pad1
  pad1.cd()               # pad1 becomes the current pad
  data_histo.SetStats(0)          # No statistics on upper plot
  
  #data_histo.Scale(mc_sum_histo.Integral()/data_histo.Integral());
  mc_stack_histo.Draw("hist e")         # Draw mc_histo on top of data_histo
  data_histo.Draw("same")               # Draw data_histo
  ATLASLabel(0.45,0.85,"Preliminary")
  legend.Draw()

  # Do not draw the Y axis label on the upper plot and redraw a small
  # axis instead, in order to avoid the first label (0) to be clipped.
  #data_histo.GetYaxis().SetLabelSize(0.)
  #axis = ROOT.TGaxis( -5, 20, -5, 220, 20,220,510,"")
  #axis.SetLabelFont(43) # Absolute font size in pixel (precision 3)
  #axis.SetLabelSize(15)
  #axis.Draw()

  # Chiara reccommendation for data_histo.GetYaxis().SetRange() ??

  # lower plot will be in pad
  c.cd()          # Go back to the main canvas before defining pad2
  pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
  pad2.SetTopMargin(0)
  pad2.SetBottomMargin(0.2)
  pad2.SetGridx() # vertical grid
  pad2.Draw()
  pad2.cd()       # pad2 becomes the current pad

  # Define the ratio plot
  h3 = data_histo.Clone("h3")
  h3.SetLineColor(ROOT.kBlack)
  h3.SetMinimum(0)  # Define Y ..
  h3.SetMaximum(2) # .. range
  h3.Sumw2()
  h3.SetStats(0)      # No statistics on lower plot


  h3.Divide(mc_sum_histo)
  h3.SetMarkerStyle(21)
  h3.Draw("ep")       # Draw the ratio plot

  # data_histo settings
  data_histo.SetLineColor(ROOT.kBlue+1)
  data_histo.SetLineWidth(2)

  # Y axis data_histo plot settings
  data_histo.GetYaxis().SetTitleSize(20)
  data_histo.GetYaxis().SetTitleFont(43)
  data_histo.GetYaxis().SetTitleOffset(1.55)

  # mc_sum_histo settings
  #mc_sum_histo.SetLineColor(ROOT.kRed)
  #mc_sum_histo.SetLineWidth(2)

  # Ratio plot (h3) settings
  h3.SetTitle("") # Remove the ratio title

  # Y axis ratio plot settings
  ratio_name = 'data / mc_stack'
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

  c.SaveAs(distribution)


################################################################################################


# Analysis defined b-tagging regions

regions_list = {
                "hbbisr_GhostVR30Rmax4Rmin02TrackJet_fj0pt480_fj1pt250_tjet2_ptsort0_hpt480_nbtag2in2ext_MV2c10_FixedCutBEff_77",
                "hbbisr_GhostVR30Rmax4Rmin02TrackJet_fj0pt480_fj1pt250_tjet2_ptsort0_hpt480_nbtag0in2ext_MV2c10_FixedCutBEff_85",
                "hbbisr_GhostVR30Rmax4Rmin02TrackJet_fj0pt480_fj1pt250_tjet2_ptsort0_hpt480_nbtag0in2ext_MV2c10_FixedCutBEff_77",
                "hbbisr_GhostVR30Rmax4Rmin02TrackJet_fj0pt480_fj1pt250_tjet2_ptsort0_hpt480_nlnotbtag2in2ext_MV2c10_FixedCutBEff_77",
                "hbbisr_GhostVR30Rmax4Rmin02TrackJet_fj0pt480_fj1pt250_tjet2_ptsort0_hpt480_nlnotbtag1in2_MV2c10_FixedCutBEff_77",
                "hbbisr_GhostVR30Rmax4Rmin02TrackJet_fj0pt480_fj1pt250_tjet2_ptsort0_hpt480_nbtag1in2ext_MV2c10_FixedCutBEff_77",
               }

# Analysis defined mass regions for HCand / ZPrimeCand
mass_ranges_list = {#"M50to100",
#                    "M75to125",
#                    "M100to150",
#                    "M125to175",
#                    "M150to200",
#                    "M175to225",
#                    "M200to250",
#                    "M225to275",
#                    "M250to300",
                   }

# List all the variables you want stack plots made of
variables_to_plot = {"Hcand_m",
                    }

# Make TStack plot for each variable in all mass_range + region combinations
stack_hist_dict = {}

input_file_list = [sys.argv] # Get files from the command line

if len(sys.argv) < 2:
  print "Please incidcate physics process files containing histograms to be stacked"
  print "Format: python stack_plot.py file1.root file2.root (etc...)"
  # Generally this looks like: python stack_plot.py folder/*.root
  exit()

hist_handle = ROOT.TH1F() # Storage variable to dump histograms into using TDirectory::GetObject()
my_input_file_list = copy.deepcopy(input_file_list[0]) # Make our own copy, leave original alone
my_input_file_list.pop(0) # Get rid of first object which should be the script itself

stack_legend = ROOT.TLegend(0.8,0.75,0.95,0.95)

luminosity = 54.11

for i,iFile in enumerate(my_input_file_list):
  my_input_file = ROOT.TFile(iFile,"read")
  for j,iRegion in enumerate(regions_list):
    for k,iVar in enumerate(variables_to_plot): # Iterate through list of variables you wish to plot for each region
      # First we get the mass region inclusive histograms
      my_input_file.GetObject(iRegion+"/"+iVar,hist_handle)
      if "data" not in iFile:
        # Scale monte carlo by luminosity of data we're comparing to
        hist_handle.Scale(1000 * luminosity) 

      hist_handle.SetName(iFile)
      # Generate name for inclusive mass region histograms
      histogram_name = iRegion + "/" + iVar
      if i == 0:
        # Create dictionary indexed by histogram_name.  Each entry has 0) A THStack 1) A list of copies of the included histograms 2) A data histogram
        stack_hist_dict[histogram_name] = [ROOT.THStack(histogram_name,histogram_name),[copy.deepcopy(hist_handle)],ROOT.TH1F()]
      else:
        stack_hist_dict[histogram_name][1].append(copy.deepcopy(hist_handle))

      if "data" in iFile:
        # Don't add data to stack and doesn't need to choose a color
        stack_hist_dict[histogram_name][2] = copy.deepcopy(hist_handle)
      else:
        choose_color(iFile, stack_hist_dict[histogram_name][1][i])
        stack_hist_dict[histogram_name][0].Add(stack_hist_dict[histogram_name][1][i])

      if j == 0:
        # Use python's sub command and a regex expression to reduce histogram name for legend
        histogram_legend_name = re.sub(r'.*hist-', '', iFile)
        histogram_legend_name = re.sub(r'.root', '', histogram_legend_name)
        stack_legend.AddEntry(stack_hist_dict[histogram_name][1][i],histogram_legend_name)

# Not currently looking at individual mass ranges
###################################################################################################
#      for l,iMassRegion in enumerate(mass_ranges_list):
#        # Next we go through the individual mass region histograms
#        my_input_file.GetObject(iRegion+"/"+iMassRegion+"/"+iVar,hist_handle)
#        hist_handle.SetName(iFile)
#        # Generate name for exclusive mass region histograms
#        histogram_name = iRegion + "/" + iMassRegion + "/" + iVar
#        if i == 0:
#          # Create dictionary indexed by histogram_name.  Each entry has a THStack and a list of copies of the included histograms
#          stack_hist_dict[histogram_name] = [ROOT.THStack(histogram_name,histogram_name),[copy.deepcopy(hist_handle)]]
#        else:
#          stack_hist_dict[histogram_name][1].append(copy.deepcopy(hist_handle))
#        stack_hist_dict[histogram_name][1][i].SetFillColor(i+1)
#        stack_hist_dict[histogram_name][0].Add(stack_hist_dict[histogram_name][1][i])
###################################################################################################

print "\n ********** BEGIN PLOTTING ********** \n"
      
stack_hist_dict_key_list = stack_hist_dict.keys()
canvas = ROOT.TCanvas("Canvas for Stack Histogram Plotting","Canvas for Stack Histogram Plotting",800,600)
#canvas.SetLogy()

# Draw 
for m,iStack in enumerate(stack_hist_dict_key_list):
  print "\n Running: " + iStack
  stack_hist_dict[iStack][0].Draw("hist e") # Must specify "hist" for filled histograms, "e" gives errors error bars
  stack_hist_dict[iStack][0].GetXaxis().SetRange(stack_hist_dict[iStack][0].GetXaxis().FindBin(1.0) - 1,stack_hist_dict[iStack][0].GetXaxis().FindBin(2.0) - 1)
  stack_hist_dict[iStack][0].GetXaxis().SetTitle("Z' candidate large-R jet mass [GeV]")
  stack_hist_dict[iStack][0].Draw("hist e") # Must specify "hist" for filled histograms, "e" gives errors error bars
  stack_hist_dict[iStack][2].Draw("same")
  string = iStack + ".pdf"
  string = string.replace("/",".")
  stack_legend.Draw()
  
  # Draw only the stack plot
  #canvas.SaveAs("no_ratio_" + string)

  #ratio_plot(data_histo ,mc_sum_histo ,mc_stack_histo ,legend ,distribution)
  ratio_plot(stack_hist_dict[iStack][2] ,stack_hist_dict[iStack][0].GetStack().Last() ,stack_hist_dict[iStack][0] , stack_legend ,string)

