import pandas as pd
import os
from pandas_profiling import ProfileReport
 

# os.chdir(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
os.chdir('/home/adh/_projects/eWRIMS_parties')
cwd = os.getcwd()
 

# Full_list_InFile = os.path.join(cwd, "water-rights-parties.csv")
Full_list_InFile = os.path.join(cwd, "EWRIMS_PARTY_DUP.xlsx") # Don't directly link to specific folders!!!!
 
##################  reading in full non-QAed data
full_data = pd.read_excel(Full_list_InFile)
# full_data = pd.read_csv(Full_list_InFile)
 
# profile = ProfileReport(full_data, title="Pandas Profiling Report")
profile = ProfileReport(full_data, title="Pandas Profiling Report", explorative=True)
# profile.to_widgets()
# profile.to_notebook_iframe()
profile.to_file("EWRIMS_PARTY_DUP-report.html")
