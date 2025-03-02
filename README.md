# NCBI2GO
Convert GB files in NCBI to structured data

# Usage Instructions

# 1. Import the Script into PyCharm
Open PyCharm.
Go to File > Open, then select the folder containing "NCBI2GO.py" and click Open.
Alternatively, you can drag and drop the script directly into your PyCharm project.

# 2. Install Required Packages
Ensure you have the following dependencies installed. If not, use the following commands to install them:
    pip install biopython         # Install Bio package  
    pip install pandas            # Install pandas package  
    pip install taxonomy-ranks    # Install taxonomy_ranks package

# 3. Set File Paths
Modify the script to specify the correct file paths based on your directory structure. Example:
  GenBank file path
  gb_file = "DataCleaning/Aves_RefSeq.gb"
  Excel file path
  excel_file = "DataCleaning/Aves_RefSeq.xlsx"
  Read the normalization dictionary, or use an empty dictionary if unavailable
  dict_file = "DataCleaning/EmptyDictionary.xlsx"

# 4. Run the Script
After updating the file paths, save the script.
In PyCharm, click the green "Run" button to execute the script, or right-click the file and select Run 'NCBI2GO.py'.

