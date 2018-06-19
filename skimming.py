import os
import sys
import ROOT
import copy
import re
from array import *

##ROOT.gROOT.LoadMacro("/afs/cern.ch/user/j/jpasner/atlasstyle/AtlasStyle.C")

# INFO:
#   Binning for 2B2J or 1B1J is hardcoded in the remove_whitespace function
#   Merging last 2 or 3 bins of 1B1J variable is hardcoded in the merge_bins function
#   To turn on or off merging you have to comment two lines in the code.  One for nominal samples, one for the systematics
#   Edit keep_list to indicate which process to keep


# This implementation is specifically for the Eq2bTagSel_diBJets_btagSumFixBin variable
def remove_whitespace(histogram):
  "Creates new copy of histogram with smaller x-axis"

  ##### FOR REFERENCE: Federico's original bin edges #####
    # 2B2J Reference:
    #bin_edges = {-2, -0.8241525, -0.354075, -0.1755727, -0.065094, -0.0230671, 0.351695, 
    #           0.8217725, 1.0002748, 1.1107535, 1.1527804, 1.29185, 1.4703523, 1.580831, 
    #           1.6228579, 1.6488546, 1.7593333, 1.8013602, 1.869812, 1.9118389, 1.9538658, 2}
    # 1B1J Reference:
    # btagFixBin bin edges to remove empty bins in data so that WSMaker works (hopefully!)
    #bin_edges = [-1., 0.1758475, 0.645925, 0.8244273, 0.934906, 0.9769329, 0.9977155, 1.0]
  ########################################################

  # 2B2J bin edges and bin array offset
#  bin_offset = 11
#  bin_edges = [1.29185, 1.4703523, 1.580831, 1.6228579, 1.6488546, 1.7593333, 1.8013602, 1.869812, 1.9118389, 1.9538658, 2]

  # 1B1J bin edges and bin array offset
  bin_offset = 2
  #Chiara: this is true for old Z inputs  bin_edges = [0.645925, 0.8244273, 0.934906, 0.9769329, 0.9977155, 1.0]
  bin_edges = [0.645925, 0.8244273, 0.934906, 0.9769329, 1.0] #Chiara: true for new inputs
  
  bin_array = array('d',bin_edges)
  n_bins = len(bin_array) - 1

  #new_histogram = ROOT.TH1F(histogram.GetName(),histogram.GetTitle(),n_bins,bin_array)
  new_histogram = ROOT.TH1F(histogram.GetName(),histogram.GetTitle(),n_bins,bin_edges[0],bin_edges[-1])

  for i in range(n_bins+1):
    new_histogram.SetBinContent(i ,histogram.GetBinContent(i+bin_offset))
    new_histogram.SetBinError(i ,histogram.GetBinError(i+bin_offset))
  
  new_histogram.SetEntries(histogram.GetEntries())

  return new_histogram

################################################################################################
def merge_bins(histogram):
  "Allows for merging of histogram bins"


  # Merge last 2 bins of 4 bin 1B1J variable
  #new_histogram = ROOT.TH1F(histogram.GetName(),histogram.GetTitle(),3,histogram.GetXaxis().GetXmin(),histogram.GetXaxis().GetXmax())
  #new_histogram.SetBinContent(1, histogram.GetBinContent(1))
  #new_histogram.SetBinContent(2, histogram.GetBinContent(2))
  #new_histogram.SetBinContent(3, histogram.GetBinContent(3) + histogram.GetBinContent(4))

  #new_histogram.SetBinError(1, histogram.GetBinError(1))
  #new_histogram.SetBinError(2, histogram.GetBinError(2))
  #new_histogram.SetBinError(3, histogram.GetBinError(3) + histogram.GetBinError(4))

  #new_histogram.SetEntries(histogram.GetEntries())

  # Merge last 2 bins of 5 bin 1B1J variable
  new_histogram = ROOT.TH1F(histogram.GetName(),histogram.GetTitle(),4,histogram.GetXaxis().GetXmin(),histogram.GetXaxis().GetXmax())
  new_histogram.SetBinContent(1, histogram.GetBinContent(1))
  new_histogram.SetBinContent(2, histogram.GetBinContent(2))
  new_histogram.SetBinContent(3, histogram.GetBinContent(3))
  new_histogram.SetBinContent(4, histogram.GetBinContent(4) + histogram.GetBinContent(5))

  new_histogram.SetBinError(1, histogram.GetBinError(1))
  new_histogram.SetBinError(2, histogram.GetBinError(2))
  new_histogram.SetBinError(3, histogram.GetBinError(3))
  new_histogram.SetBinError(4, histogram.GetBinError(4) + histogram.GetBinError(5))

  new_histogram.SetEntries(histogram.GetEntries())

  # Merge last 3 bins of 5 bin 1B1J variable
  
  #new_histogram = ROOT.TH1F(histogram.GetName(),histogram.GetTitle(),3)
  #new_histogram.SetBinContent(0, histogram.GetBinContent(0))
  #new_histogram.SetBinContent(1, histogram.GetBinContent(1))
  #new_histogram.SetBinContent(2, histogram.GetBinContent(2) + histogram.GetBinContent(3))

  #new_histogram.SetBinError(0, histogram.GetBinError(0))
  #new_histogram.SetBinError(1, histogram.GetBinError(1))
  #new_histogram.SetBinError(2, histogram.GetBinError(2) + histogram.GetBinError(3))

  #new_histogram.SetEntries(histogram.GetEntries())

  return new_histogram

################################################################################################
def import_list_of_systematics():
  "Returns the number of lines (systematics) read in and the names of the systematics themselves"

  print "--- Opening systematics_to_extract.txt ---"

  # Check that systematics_to_extract.txt is open, otherwise exit
  try: 
    systematicsListFile = open("systematics_to_extract.txt")
  except IOError:
    print "No systematicsListFile found. Exiting . . ."
    exit()

  lines = systematicsListFile.read().splitlines()
  print "List of Systematics found inside systematicsListFile:\n"
  print lines
  print "\nThese " + str(len(lines)) + " systmeatics will be extracted for each distribution if they exist."

  print "------------------------------------------\n"

  return len(lines), lines 

################################################################################################

input_file_list = [sys.argv] # Get files from command line

if len(sys.argv) < 2:
  print "Please indicate which hist files you would like to skim with the format: python skimming.py file1.root file2.root (etc...)" 
  #generrally *.root
  exit()

my_input_file_list = copy.deepcopy(input_file_list[0]) # Personal copy
my_input_file_list.pop(0) # fixes bug where script would think the first file was the script itself 

# Read systematics_to_extract.txt to find list of systematics to extract for each distribution
number_of_systematics,list_of_systematics = import_list_of_systematics()

# Variables you wish to remove (order by most exclusionary first)
#remove_list = ['El']
remove_list = []

# Variables you wish to keep from files
#keep_list = ['_Mu_Eq2bTagSel_diBJets_btagSumFixBin']
#keep_list = ['Mu_2bTagSel_diBJets_btagSumFixBin'] #Camilla's version doesn't have the 'Eq'?
#keep_list = ['Mu_1bTagSel_leadBJet_btagFixBin'] #Region for single b-tag tests
#keep_list = ['El_Eq1bTagSel_leadBJet_btagFixBin'] #For semen's 31 January preliminary systematics sample

# New Naming convention variable to keep from file
# 2B2J format
#keep_list = ['El_2B2J_diBJets_bTagWeightSumOf2Jets_inCalibBin'] #Z 2B2J el
#keep_list = ['Mu_2B2J_diBJets_bTagWeightSumOf2Jets_inCalibBin']  #Z 2B2J mu
# 1B1J format
#keep_list = ['El_1B1J_leadBJet_bTagWeight_inCalibBin'] #Z case - 1B1J el
#keep_list = ['Mu_1B1J_leadBJet_bTagWeight_inCalibBin'] #Z case - 1B1J mu

##keep_list = ['Mu_1B1JHighMTW_leadJet_bTagWeight_inCalibBin'] #W case - 1B1J - muon
#keep_list = ['El_1B1JHighMTW_leadJet_bTagWeight_inCalibBin'] #W case - 1B1J - electron
#keep_list = ['El_2B2JHighMTW_diBJets_bTagWeightSumOf2Jets_inCalibBin'] #W case - 2B2J - electron
keep_list = ['Mu_2B2JHighMTW_diBJets_bTagWeightSumOf2Jets_inCalibBin'] #W case - 2B2J - muon


# Print objects to keep and objects to remove:
print "Keeping histograms containing: " + str(keep_list)
print "Ignoring histograms containing: " + str(remove_list) + "\n"

# Make output file for skimming results if it doesn't already exist
file_path = os.getcwd() + "/output/"
directory = os.path.dirname(file_path)

try:
  os.stat(directory)
except:
  os.mkdir(directory)  

# Loop over all files specified in command line input
for i,iFile in enumerate(my_input_file_list):
  my_input_file = ROOT.TFile(iFile,"read")

  # Check that file exists, exit if it doesn't
  if not my_input_file.IsOpen():
    print "Couldn't open this file: "
    my_input_file.Print()
    break

  # Get list of all objects (keys) in file
  input_object_list = my_input_file.GetListOfKeys()

  # Check that file contains object keys, exit if it doesn't
  if not input_object_list:
    print "No object keys found in file: "
    my_input_file.Print()
    break

  # Output file creation
  input_file_name = my_input_file.GetName()
  print "--- Opening root file: " + input_file_name + " ---"
  output_file_name = 'output/' + re.sub(r"^.*hist-","",input_file_name)
  output_file = ROOT.TFile(output_file_name,"recreate")

  # Check if input file contains Systematics subdirectory
  if my_input_file.GetDirectory("Systematics"):
    print "  *** " + input_file_name + " contains Systematics subdirectory ***"
    output_file.mkdir("Systematics/")
    systematics_flag = True
  else:
    systematics_flag = False

  # Loop over all object keys in file
  for j in range(input_object_list.GetSize()):

    # Remove all histograms that we don't want (specified above)
    checker = False
    for k,iUnwanted in enumerate(remove_list):

      # Check to see if histogram name includes a variable we don't want to keep
      if iUnwanted in input_object_list.At(j).GetName():
        checker = True

    # If histogram contains variable we don't want to keep, we move to the next histogram
    if checker:
      continue

    # Search for all desired variables (specified above)
    for l,iWanted in enumerate(keep_list):

      # Look for any object keys that contain the desired variable name, Write to file
      if iWanted in input_object_list.At(j).GetName():        
        histogram = my_input_file.Get(input_object_list.At(j).GetName()) # extract histogram from file

        # Name and axis formatting stuff goes here
        histogram = remove_whitespace(histogram) # Re-bins the inputs into equal width bins and removes unfilled bins from btagweight distributions  

##Chiara: activate only if the merging is needed!
##        histogram = merge_bins(histogram) # Meges togeather bins in histogram.  Currently set up for merging last 2 or last 3 bins

        new_histogram_name = histogram.GetName().replace(iWanted,'') # Removes name of variable from histogram name

        # If statmentes to remove montecarlo labales to make sure final variable names are correct
        if '_improved' in new_histogram_name:
          new_histogram_name = new_histogram_name.replace('_improved','')

        if 'v22' in new_histogram_name:
          new_histogram_name = new_histogram_name.replace('v22','_v22')

        if 'MadZee' in new_histogram_name:
          new_histogram_name = new_histogram_name.replace('MadZee','Z')

        if 'MadWtaunu' in new_histogram_name:
          new_histogram_name = new_histogram_name.replace('MadWtaunu','W')

        if 'MadWenu' in new_histogram_name:
          new_histogram_name = new_histogram_name.replace('MadWenu','W')

        if 'MadZmumu' in new_histogram_name:
          new_histogram_name = new_histogram_name.replace('MadZmumu','Z')

        if 'MadWmunu' in new_histogram_name:
          new_histogram_name = new_histogram_name.replace('MadWmunu','W')

        temp_string_list = new_histogram_name.split('_')

        # Fix formatting for flavours so that WSMaker postfit plots work properly (parser doesn't like '_')
        if '_L' in new_histogram_name:
          new_histogram_name = temp_string_list[0] + '1l'
        elif '_1C' in new_histogram_name:
          new_histogram_name = temp_string_list[0] + '1c'
        elif '_1B' in new_histogram_name:
          new_histogram_name = temp_string_list[0] + '1b'
        elif '_NC' in new_histogram_name:
          new_histogram_name = temp_string_list[0] + 'cc'
        elif '_NB' in new_histogram_name:
          new_histogram_name = temp_string_list[0] + 'bb'
        elif temp_string_list[0] == 'W':
          new_histogram_name = temp_string_list[0] + 'jets'
        elif temp_string_list[0] == 'Z':
          new_histogram_name = temp_string_list[0] + 'jets'
        else:
          new_histogram_name = temp_string_list[0]

# Note that there are different methods for extracting systematics
# Exclusive Option: List systematics you want to extract for each distribution in systematics_to_extract.txt
#   Only listed systematics will be extracted (if they exist!)
# Inclusive Option: Iterate accross all systematics and extract every one that exists for the given distributions
#   This will give you all of the systematics that were generated, but can take a LONG TIME
################################################################################################

        # Extract systematics listed inside of systematics_to_extract.txt
        # From above: number_of_systematics and list_of_systematics were extracted by import_list_of_systematics() from systematics_to_extract.txt
        if systematics_flag:
          sys_directory = my_input_file.GetDirectory("Systematics")
          
          for n in range(number_of_systematics):
            for variation in ['__1up','__1down']: # Need to extract both up and down systematic histograms
              sys_histogram_name = histogram.GetTitle() + "_" + list_of_systematics[n] + variation # Construct name of systematic histogram you intend to extract
              sys_histogram = sys_directory.Get(sys_histogram_name) # extract SPECIFIC systematic histogram from Systematics subdirectory

              if sys_histogram:
                sys_histogram = remove_whitespace(sys_histogram) # Re-bins the inputs into equal width bins and removes unfilled bins from btagweight distributions 
                sys_histogram = merge_bins(sys_histogram) # Meges togeather bins in histogram.  Currently set up for merging last 2 or last 3 bins
                new_sys_histogram_name = sys_histogram.GetName().replace(histogram.GetTitle(),'') # Get full systematic histogram name and isolate Systematic description by removing sample/region name
                new_sys_histogram_name = new_histogram_name + new_sys_histogram_name # Prepend name of sample onto Systematic Histogram Name

                print '  + Saving Systematics histogram: ' + sys_histogram.GetTitle() + ' as ' + new_sys_histogram_name
                output_file.cd("Systematics/")
                sys_histogram.SetName(new_sys_histogram_name)        
                sys_histogram.Write()
                sys_histogram.Delete()
                output_file.cd()

              else: 
                print "  - Systematic Histogram doesn't exist: " + sys_histogram_name
 
################################################################################################

        # Extract ALL sytematics (this takes much longer!!!)
#        if systematics_flag:
#          sys_directory = my_input_file.GetDirectory("Systematics")
#          sys_object_list = sys_directory.GetListOfKeys()
#
#          for n in range(sys_object_list.GetSize()):
#            print "Number " + str(n) + "out of " + str(sys_object_list.GetSize()) + " " + sys_object_list.At(n).GetName() 
#            if histogram.GetTitle() in sys_object_list.At(n).GetTitle():
#              sys_histogram = sys_directory.Get(sys_object_list.At(n).GetName()) # extract systematic histogram from Systematics subdirectory
#              sys_histogram = remove_whitespace(sys_histogram)
#              new_sys_histogram_name = sys_histogram.GetName().replace(histogram.GetTitle(),'') # Get full systematic histogram name and isolate Systematic description by removing sample/region name
#              new_sys_histogram_name = new_histogram_name + new_sys_histogram_name # Prepend name of sample onto Systematic Histogram Name
#
#              print '  Saving Systematics histogram: ' + sys_object_list.At(n).GetName() + ' as ' + new_sys_histogram_name
#              output_file.cd("Systematics/")
#              sys_histogram.SetName(new_sys_histogram_name)        
#              sys_histogram.Write()
#              sys_histogram.Delete()
#              output_file.cd()
        
################################################################################################

        # Save sample histogram inside of base file (not in Systematics subdireectory)
        print 'Saving histogram: ' + input_object_list.At(j).GetName() + ' as ' + new_histogram_name
        histogram.SetName(new_histogram_name)        
        histogram.Write()
        histogram.Delete()

  # Write file and close
  print "+++ Writing root file: " + output_file_name + " +++ \n"
  output_file.Write()
  output_file.Close()
  my_input_file.Close()

print "*** Done ***"

