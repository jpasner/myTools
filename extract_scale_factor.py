import sys
import copy

#
# This tool was developed by Jake Pasner to extract the Scale Factors, Errors and Relative Erros
# from the output files of the fit cross checking in the VHF WSMaker.
# example file: fccs/FitCrossChecks_Z_unmerged_fullSys_MG.FloatBCLcc_VHF_13TeV_2B2J_Z_Combination_diBJets_bTagWeightSumOf2jets_inCalibBin_combined/output_0.log
#

output_file = open("output.csv", "w") # Comma separated variable files are the easiest format to bring into keynote or excel 
input_file_list = [sys.argv] # Get fcc Scale Factor (SF) files from command line

if len(sys.argv) < 2:
  print "Please indicate which fcc output.log files you would like to skim with the format: python extract_scale_factor.py output_0.log (etc...)"
  sys.exit()

my_input_file_list = copy.deepcopy(input_file_list[0]) # Personal copy
my_input_file_list.pop(0) # fixes bug where script would think the first file was the script itself 

for i,input_file in enumerate(my_input_file_list): # Loop through all fcc files we want to extract SF from
  output_file.write("Input file: " + my_input_file_list[i] + "\n") #Write name of fcc file before extracting data

  output_file.write("Distribution, Final Value, Error, Relative Error \n") # Setup of comma separated variable columns
  my_input_file = open(input_file,"r")
  flag = 0
  for line in my_input_file:
    if flag == 1: # Flag that indicates begining of section containing SF results from fit
      output_vector = line.split()
      if "alpha" in line: # The nomralisation SF always appear before the results for alpha systematics so we leave the file when we see "alpha"
        break
      print line
      print output_vector
      if len(output_vector) == 4:
        output_file.write("{}, {}, {}, {}\n".format(output_vector[0], float(output_vector[1]), float(output_vector[3]), float(output_vector[3])/float(output_vector[1]))) #Formatting of data
      else:
        print "Issue with formatting of line: " + line # A little bit of bug finding code
        print "Exiting program"
        sys.exit()
    if "  --------------------  --------------------------" in line: # This dashed line indicates the begining of the results section of the fcc output file
      flag = 1
  my_input_file.close()
  print "\n"

output_file.close()
