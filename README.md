# EIMS-Peptide-Fragment-Predictor
Makes predictions of fragments generating when subjecting peptides to electron ionization mass spectrometry
The program is run from the command line. To run, navigate to the folder containing the python script (predictFragments2023.py). The monoisotopic element mass definition file (monoisotopic_element_masses.csv), element index file (element_index.csv), fragment definition text file (AA_data.csv), and license file (license.txt) should also be present in the same directory as the python script file.

Run the script using Python 3, for example:
  python3 predictFragments2023.py
  
The main menu will be displayed once the program runs. Entering the integer associated with the menu option and pressing Enter will execute the indicated action(s) in the program. A csv file will be saved with results in the source code folder after predictions are generated.

Note that to conduct multiple fragmentation predictions (second program option), you will need to have prepared a text file containing the sequences you want to conduct the analysis on listed one sequence per line in one letter code.
