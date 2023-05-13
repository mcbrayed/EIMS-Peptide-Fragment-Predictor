#predictFragments is a program that generates predicted fragment m/z values for peptides subjected to electroionization mass spectrometry
#    Copyright (C) 2023  Dominic McBrayer

    #This program is free software: you can redistribute it and/or modify
   # it under the terms of the GNU General Public License as published by
    #the Free Software Foundation, version 3 of the License
    
    #This program is distributed in the hope that it will be useful,
  # but WITHOUT ANY WARRANTY; without even the implied warranty of
   # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   # GNU General Public License for more details.

   # You should have received a copy of the GNU General Public License
   # along with this program.  If not, see <https://www.gnu.org/licenses/>.


#program code prepared by: Dominic McBrayer and Christina Signoretti 2020
                          #Dominic McBrayer 2021
                          #Dominic McBrayer 2022
                          #Dominic McBrayer 2023



#****************IMPORTS********************************
import tkinter #can be used to build a graphical user interface
from tkinter import filedialog #used to access the system's graphical file manager for locating files to load (more user friendly for users not comfortable with command line prompts)


#*****************CLASS DEFINITIONS**********************************
class AminoAcid():
    def __init__(self,name,formula,sidechain_class,frags):
        self.name = name
        self.formula = formula #the formula for the amino acid (gets folded into the peptide's overall formula)
        self.sidechain_class = sidechain_class #e.g. alcohol, acidic
        self.fragmentations = frags #list of different loss formulas
        self.frag_label = ""

    def set_formula(self,formula):
        self.formula = formula

    def get_frag_label(self):
        return self.frag_label

    def __str__(self):
        return self.name

    def add_formula(self,formula):
        final = []
        for i in range(len(self.formula)):
            final.append(self.formula[i]+formula[i])
        self.formula = final

    def get_name(self):
        "'Returns the name designation for the amino acid'"
        return self.name

    def sub_formula(self,formula):
        from operator import sub
        sub = list(map(sub, self.formula, formula))
        self.formula = sub

    def get_formula(self):
        return self.formula

    def fragment(self): #will focus on side chain fragmentation
        aa_frags = [self.duplicate()]
        for frag in self.fragmentations:
            if sum(frag[1]) == 0: #frag[1] references the formula list not the string label
                pass
            else:
                temp_aa = self.duplicate()
                if "charged" in frag[0].lower(): #if the fragmentation description contains "charged" keyword, indicating a nonrelative frag event, will just generate a fragment of the defined formula rather than subtracting the formula
                    temp_aa.set_formula(frag[1])
                else:
                    temp_aa.sub_formula(frag[1])
                temp_aa.frag_label += "("+frag[0]+")" #update amino acid's name with the fragmentation type
                aa_frags.append(temp_aa)
        return aa_frags

    def duplicate(self): #reinitializes an AminoAcid object with the same values as the original and returns the duplicate
        #ALWAYS update this function after modifying the class
        copy_AA = AminoAcid(self.name,self.formula,self.sidechain_class,self.fragmentations)
        return copy_AA

class PepFrag():
    def __init__(self,amino_acids,fragmentation_record = "Initial",form_peptide_bonds = False):
        self.amino_acids = amino_acids #in a sequence form
        self.fragmentation_record = fragmentation_record
        if form_peptide_bonds == True:
            self.form_peptide_bonds()
        self.update_formula()

    def __len__(self): #generates and returns the length of the peptide fragment
        return len(self.amino_acids)

    def get_seq_str(self):
        s = ""
        for aa in self.amino_acids:
            s += str(aa)
        return s
    
    def get_frag_record(self):
        return self.fragmentation_record

    def get_csv_str(self,element_conver_dict,index_conver_dict): #generates tab-delimited output for the fragment
        s = ""
        s += self.get_seq_str()
        s += "\t"
        s += str(self.get_mz_value(element_conver_dict,index_conver_dict))+"\t"
        frag_record_list = self.fragmentation_record.split("\n")
        log = ""
        for entry in frag_record_list:
            log += entry+"|"
        s += log+"\n"
        return s

    def __str__(self): #generates a text-based representation of the PepFrag object
        s = ""
        s += self.get_seq_str()
        frag_record_list = self.fragmentation_record.split("\n")
        for entry in frag_record_list:
            s += "\n\t"+entry
        return s

    def duplicate(self):
        #ALWAYS update this function after modifying the class
        copy_of_amino_acids = []
        for aa in self.amino_acids:
            copy_of_amino_acids.append(aa.duplicate())
        copy_of_pep_frag = PepFrag(copy_of_amino_acids, self.fragmentation_record)
        return copy_of_pep_frag

    def update_formula(self): #resets the formula attribute and then recalculates based off status of current amino acids
        self.formula = [0,0,0,0,0,0]
        for aa in self.amino_acids:
            self.add_formula(aa.get_formula())

    def form_peptide_bonds(self):
        for i in range(len(self.amino_acids)):
            self.amino_acids[i].name += str(i+1)
            if (i+1) >= len(self.amino_acids):
                break
            self.amino_acids[i].sub_formula([0,1,0,1,0,0]) #subtracts OH from the formula for the first amino acid
            self.amino_acids[i+1].sub_formula([0,1,0,0,0,0]) #subtracts H from the formula of the next amino acid
        self.update_formula() # will go through amino acids and add up their formulas once corrected

    def add_formula(self,formula): #takes an atom count formula and adds its values to those of the current formula
        final = []
        for i in range(len(self.formula)):
            final.append(self.formula[i]+formula[i])
        self.formula = final

    def get_formula(self):
        return self.formula

    def sub_formula(self,formula): #takes an atom count formula and subtracts its values from those of the current formula
        from operator import sub
        sub = list(map(sub, self.formula, formula))
        self.formula = sub

    def fragment(self): # Connects everything together to step systematically through all fragmentation types
        print("////////////Conducting Initial fragmentation\\\\\\\\\\\\\\\\")
        total_fragments = []
        #total_fragments = [self.duplicate()]
        type_A_frags = self.type_A_frag() #do not add to total fragments (will be included as part of side chain analysis)
        print("Found a total of "+str(len(type_A_frags))+" type A fragments (not counting side chains)")
        type_A_side_chain_frags = []
        for frag in type_A_frags:
            type_A_side_chain_frags.extend(frag.side_chain_frag(False)) #since type A loses C-terminus, will not include C-terminus in analysis
        total_fragments.extend(type_A_side_chain_frags)
        print("Found a total of "+str(len(type_A_side_chain_frags))+" fragments after accounting for side chains")
        type_B_frags = self.type_B_frag()
        print("Found a total of "+str(len(type_B_frags))+" type B fragments (not counting side chains)")
        type_B_side_chain_frags = []
        for frag in type_B_frags:
            type_B_side_chain_frags.extend(frag.side_chain_frag())
        print("Found a total of "+str(len(type_B_side_chain_frags))+" fragments after accounting for side chains")
        total_fragments.extend(type_B_side_chain_frags)
        type_C_frags = self.type_C_frag()
        print("Found a total of "+str(len(type_C_frags))+" type C fragments (not counting side chains)")
        type_C_side_chain_frags = []
        for frag in type_C_frags:
            type_C_side_chain_frags.extend(frag.side_chain_frag())
        print("Found a total of "+str(len(type_C_side_chain_frags))+" fragments after accounting for side chains")
        total_fragments.extend(type_C_side_chain_frags)
        print("**Conducting molecular ion side chain fragmentation")
        side_chain_frags = self.side_chain_frag() #this is for molecular ions that might fragment at side chains but not along the backbone
        
        print("Found a total of "+str(len(side_chain_frags))+" side chain ('molecular') fragments")

        total_fragments.extend(side_chain_frags)
        total_fragments,remove_count = remove_duplicate_fragments(total_fragments) #remove any duplicates
        print("Initial grand total of "+str(len(total_fragments))+" possible peptide fragments")
        print("Removed "+str(remove_count)+" duplicates from grand total")
        print("Found a GRAND total of "+str(len(total_fragments))+" possible peptide fragments")
        
        return total_fragments

    def update_log(self, new_data):
        self.fragmentation_record += "\n"+new_data

    def type_A_frag(self):
        print("***Type A Fragmentation Function")
        type_A_backbone_cleavages_list = []
        for i in range(len(self.amino_acids)):
            #make a copy of Self
            pos_frag = self.duplicate()
            #line below applies to INTERNAL amino acids (need if/else style condition to handle C-terminal case)
            if i == len(self.amino_acids)-1:#do c-terminus stuff [1, 1, 0, 2, 0, -1]
                pos_frag.amino_acids[i].sub_formula([1, 1, 0, 2, 0, -1])
            else:
                pos_frag.amino_acids[i].sub_formula([1, 0, 0, 1, 0, -1]) #subtracts 1C and 1O, and adds +1 charge subtracts: [1, 0, 0, 1, 0, -1]
                pos_frag.amino_acids = pos_frag.amino_acids[:i+1] #keep the amino acids in the N-terminal (postively charged) fragment
            pos_frag.update_log(pos_frag.amino_acids[i].get_name()+"(type_a)")
            pos_frag.update_formula()
            #generated fragment (and store in a list) --> fragments contain amino acids with updated formulas corresponding to fragmentation
            type_A_backbone_cleavages_list.append(pos_frag)
        return type_A_backbone_cleavages_list

    def type_B_frag(self):
        print("***Type B Fragmentation Function")
        type_B_backbone_cleavages_list = []
        for i in range(len(self.amino_acids)):
            positive_B_frag = self.duplicate()
            if i == len(self.amino_acids)-1:
                pass # This cleaves peptide bond. Last residue does not have another residue to cleave with it
            else:
                positive_B_frag.update_log(positive_B_frag.amino_acids[i].get_name()+"(type_b)") #double check that correct amino acid is being referenced in log
                positive_B_frag.amino_acids[i+1].sub_formula([0, 0, 0, 0, 0, -1])
                positive_B_frag.amino_acids = positive_B_frag.amino_acids[i+1:] #slice selects for the amino acids in the positively charged fragment
                # This is subtracting the backbone of neg type b frag without the side chain C2H3NO-
            positive_B_frag.update_formula()
            type_B_backbone_cleavages_list.append(positive_B_frag)
        return type_B_backbone_cleavages_list

    def type_C_frag(self): # occurs for aromatic residues after the N-H on backbone. ONLY for aromatic residues. F, Y, W residues.
        print("***Type C Fragmentation Function")
        type_C_backbone_cleavages_list = []
        for i in range(len(self.amino_acids)):
            positive_C_frag = self.duplicate()
            positive_C_frag.update_log(positive_C_frag.amino_acids[i].get_name()+"(type_c)") # this gets the AA name first. We need this for if else statements
            if positive_C_frag.amino_acids[i].get_name()[0] in "FYW": #F, Y, and W are the three aromatic residues that undergo type C cleavage
                ND_C_frag = self.duplicate()
                ND_C_frag.update_log(positive_C_frag.amino_acids[i].get_name()+"(type_c-nd)") #nd stands for "non-dominant" will use to indicate lower confidence in a fragment
                if i == 0:
                    ND_C_frag.amino_acids[i].sub_formula([0,2,1,0,0,-1])
                    positive_C_frag.amino_acids[i].sub_formula([0, 3, 1, 0, 0, -1]) # subtracts the 1N and 1H off backbone and 1H for H rearrangement
                else:
                    ND_C_frag.amino_acids[i].sub_formula([0,1,1,0,0,-1])
                    positive_C_frag.amino_acids[i].sub_formula([0,2,1,0,0,-1])
                positive_C_frag.amino_acids[i].fragmentations = [("ND",[0,0,0,0,0,0])] #treat C-cleavage associated side chains as "pre-fragmented" since they are lost as part of C-cleavage, so remove frag definitions
                ND_C_frag.amino_acids[i].fragmentations = [("ND",[0,0,0,0,0,0])]
                positive_C_frag.amino_acids = positive_C_frag.amino_acids[i:] # want to keep aa in C-terminal (pos charged) fragment
                ND_C_frag.amino_acids = ND_C_frag.amino_acids[i:]
                ND_C_frag.update_formula()
                type_C_backbone_cleavages_list.append(ND_C_frag)
                positive_C_frag.update_formula()
                type_C_backbone_cleavages_list.append(positive_C_frag)
        return type_C_backbone_cleavages_list

    def side_chain_frag(self,correct_c_terminus = True):
        print("***Side chain frag function")
        frag_index_dict = {}
        all_fragments = []
        frag_site_counter = 0
        peptide_fragments = []
        c_terminus_data = [("(c_term intact)",[0,0,0,0,0,0]),("(c_term water)",[0,2,0,1,0,0]),("(c_term hydroxyl)",[0,1,0,1,0,0]),("(c_term carboxylic)",[1,1,0,2,0,0])]
        c_terminus_name = self.amino_acids[-1].get_name()
        new_frag_site = False
        for i in range(len(self.amino_acids)):
            sub_frags = self.amino_acids[i].fragment() #get a list of amino acids labeled with their frag type and with their individual frag formulas updated
            if len(sub_frags) == 0: #ignore sites that don't have defined fragments
                pass
            elif len(self.amino_acids) == 1 and correct_c_terminus == False: #looking at a single peptide and not factoring in C-terminus. So all initial possibilities are the total number of possibilities
                # print("single amino acid and not factoring c-terminus degradation in")
                for aa in sub_frags:
                    pep_frag = self.duplicate()
                    pep_frag.amino_acids[int(aa.get_name()[1:])-1] = aa
                    pep_frag.update_formula()
                    if aa.get_frag_label() != "":
                        pep_frag.update_log(aa.get_name()+aa.get_frag_label())
                    peptide_fragments.append(pep_frag)
                
            else:
                for aa in sub_frags:
                    # print("aa frag type = ",aa.get_frag_label())
                    if "charged" in aa.get_frag_label(): #charged frags should only be counted once per site - not mixed and matched. Add them to final export list now, but do not include in the combination fragment analysis
                        # print("Exclusive case identified, current frag type not undergoing mix and match process")
                        pep_frag = self.duplicate()
                        pep_frag.amino_acids = [aa.duplicate()]
                        pep_frag.update_formula()
                        pep_frag.update_log(aa.get_name()+aa.get_frag_label())
                        peptide_fragments.append(pep_frag)
                    else:
                        # print("adding amino acid frag to mix and match process")
                        all_fragments.append(aa)
                        frag_index_dict[self.amino_acids[i].get_name()] = i
                        new_frag_site = True
            if new_frag_site:
                frag_site_counter += 1 
        if correct_c_terminus == True:
            if self.amino_acids[-1].get_name() not in frag_index_dict:
                print("C-terminus not available to fragment")
                pass #indicates that C-terminus is not available to correct (e.g. the loss of side chains that retain the charge)
            else:
                print("Adding C-terminal fragments")
                for correction in c_terminus_data: #adds c-terminus fragments to C-terminal amino acid
                    temp_aa = self.amino_acids[-1].duplicate()
                    temp_aa.sub_formula(correction[1])
                    temp_aa.frag_label = correction[0]
                    all_fragments.append(temp_aa)
                frag_site_counter += 1
        import itertools
        combinations = itertools.combinations(all_fragments,frag_site_counter)
        combo_list = []
        filter_list = []
        for entry in combinations:
            temp_list = []
            s = "" #used to populate the filter list
            for aa in entry:
                s += aa.get_name()+aa.get_frag_label()
                temp_list.append(aa)
            filter_list.append(s)
            combo_list.append(temp_list)
        print("found "+str(len(combo_list))+" unfiltered combination(s)")
        # print("filter_list = ",filter_list)
        filtered_combos = []
        for i in range(len(filter_list)):
            keep = True
            # print("*******considering combination: "+filter_list[i])
            for frag_index in frag_index_dict: #use keys in the frag index dict to find aa names
                if "c_term" in filter_list[i]:
                    if filter_list[i].count("c_term") > 1:
                        keep = False
                        # print("rejecting due to multiple combined c_term degradations")
                    
                elif filter_list[i].count(frag_index) > 1:
                    # print("rejecting invalid combination = "+filter_list[i])
                    keep = False
                        
                if filter_list[i].count("side_chain") > 1:
                    # print("rejecting due to multiple combined side chain frags at one site")
                    keep = False
            if keep == True:
                filtered_combos.append(combo_list[i])
                # print("keeping the combination!")
        print("found "+str(len(filtered_combos))+" filtered combination(s)")
        #now ready to actually generate the peptide fragment_lists
        test_counter = 0
        # print("assembling fragments from amino acid fragmentation events")
        for combo in filtered_combos:
            pep_frag = self.duplicate() #make a copy of the current unfragmented peptide
            for aa in combo:
                # print("updating amino acids from combinations")
                for frag in aa.fragmentations:
                    # print("considering possible frag: ",frag)
                    if "("+frag[0]+")" == aa.get_frag_label() or frag[0] == aa.get_frag_label():
                        # print("frag matches current modification")
                        pep_frag.amino_acids[frag_index_dict[aa.get_name()]].sub_formula(frag[1])
                # print("considering c-terminus")
                for frag in c_terminus_data:
                    # print("considering possible frag: ",frag)
                    if frag[0] == aa.get_frag_label():
                        # print("frag matches current modification")
                        pep_frag.amino_acids[frag_index_dict[aa.get_name()]].sub_formula(frag[1])
                if aa.get_frag_label() == "" or aa.get_frag_label() == "(c_term intact)": #update peptide records to reflect modifications
                    # print("skipping log update, intact peptide")
                    pass
                else:
                    # print("label is not empty and is not an intact c-term")
                    pep_frag.update_log(aa.get_name()+aa.get_frag_label())
                if pep_frag.get_formula()[-1] == 0:
                    pep_frag.add_formula([0,0,0,0,0,1]) #if uncharged, add a positive charge (need to update if adding negative ion support or if formula format changes)
                pep_frag.update_formula()
                # print("***adding pep_frag***\n",pep_frag)
                peptide_fragments.append(pep_frag)
                test_counter += 1

        return peptide_fragments

    def get_monoisotopic_mass(self, element_mass_conversion_dict,index_element_conversion_dict):
        peptide_formula = self.get_formula()
        uncorrected_mass = 0
        for i in range(len(peptide_formula)-1):
            uncorrected_mass += peptide_formula[i]*element_mass_conversion_dict[index_element_conversion_dict[i]]
        corrected_mass = self.get_mass_correction(peptide_formula[-1])
        final_mono_mass = uncorrected_mass + corrected_mass
        return final_mono_mass

    def get_mass_correction(self, charge):
        if charge > 0:
            return(-0.00054858 * charge)
        elif charge < 0:
            return(0.00054858 * -charge)
        else: # accounts for when charge = 0
            return(charge)

    def get_mz_value(self,element_mass_conversion_dict,index_element_conversion_dict):
        mono_mass = self.get_monoisotopic_mass(element_mass_conversion_dict,index_element_conversion_dict)
        if self.formula[-1] == 0:
            return mono_mass
        else:
            return abs(mono_mass/self.formula[-1]) #divides the mass by the peptide's charge

#****************************FUNCTION DEFINITIONS********************************
def remove_duplicate_fragments(fragments): #can be used to remove duplicates (as determined by the string representation of fragments)
    descriptions = []
    filtered_frags = []
    removed_count = 0
    for frag in fragments:
        if str(frag) not in descriptions: #this fragment has not been included, yet, so add to the filtered list
            descriptions.append(str(frag))
            filtered_frags.append(frag)
        else: #duplicate fragment already included, skip
            removed_count += 1
    return filtered_frags,removed_count

def get_user_input():
    sequence = input("Please input peptide sequence in one-letter code format: ")
    sequence.strip()
    return sequence

def get_AA_object(sequence):
    Amino_Acid_Data = load_aa_data()
    Amino_Acids = []
    for aa in sequence:
        aa_obj = Amino_Acid_Data[aa].duplicate()
        Amino_Acids.append(aa_obj)
    return Amino_Acids #(name, formula, sidechain_class, frags)

def load_element_masses(): #requires definition file monoisotopic_element_masses.csv to be present in the root directory
    with open("monoisotopic_element_masses.csv", "r") as f:
        index_mass_conversion_dict = {}
        for line in f:
            if "Element" not in line:
                line_data = line.split("\n")
                extracted_data = line_data[0]
                comma_data = extracted_data.split(",")
                index_mass_conversion_dict[comma_data[0]] = float(comma_data[1]) #convert the mass values into float (so they are treated as numbers later)
        return index_mass_conversion_dict

def load_element_index(): #requires definition file element_index.csv to be present in the root directory
    index_element_conversion_dict = {}
    with open("element_index.csv","r") as f:
        for line in f:
            data = line.split("\n")[0].split(",") #remove the newline characters from the data string and split on the comma (should give a two entry list)
            index_element_conversion_dict[int(data[0])] = data[1] #use the first entry as an integer key and set the element's string symbol (second entry) as the value
    return index_element_conversion_dict

def choose_file_directory(type = "load"): #find the path to a file
    root = tkinter.Tk() #allows use of OS graphical interface to find the file path
    root.withdraw() #hide the tk window
    initial_dir = load_default_folder_path()
    if type == "folder":
        directory = filedialog.askdirectory(initialdir = initial_dir)
        if directory != "":
            save_default_folder_path(directory)
        return directory
    if type == "load":
        directory = filedialog.askopenfilename(initialdir = initial_dir)
    if type == "save":
        directory = filedialog.asksaveasfilename(initialdir = initial_dir,filetypes = [('CSV','.csv')])
        if ".csv" not in directory: #always export files as .csv files
            directory += ".csv"
    if directory != "":
        path = directory.split("/")
        path.pop(-1) #remove the actual file name from the path
        new_default_folder_path = ""
        for x in range(len(path)):
            new_default_folder_path += path[x]+"/"
        save_default_folder_path(new_default_folder_path)
    return directory

def load_default_folder_path(): #opens a text file that stored the last directory accessed when loading a file
    try:
        with open("lastdir.txt","r") as f:
            default_folder = next(f)
        return default_folder
    except FileNotFoundError:
        return "."
    except StopIteration:
        return "."

def save_default_folder_path(path_string):
    with open("lastdir.txt","w") as f:
        f.write(path_string)

def load_seq_list(directory):
    ls = []
    with open(directory,"r") as f:
        for line in f:
            ls.append(line.split("\n")[0])
    return ls

def run_multiple_predictions(seq_list,e_mass_convert_dict,index_convert_dict,output_directory = ".//MultiFragResults.csv"):
    print("Running multiple predictions using provided sequence list")
    output = "Frag Seq\tPredicted m/z\tLog\n"
    for seq in seq_list:
        output += "NEW SEQUENCE ("+seq+")\n"
        AA_obj = get_AA_object(seq)
        peptide = PepFrag(AA_obj,form_peptide_bonds = True)
        fragments = peptide.fragment()
        fragments.append(peptide)
        fragments = remove_negative_fragments(fragments,e_mass_convert_dict,index_convert_dict)
        fragments,remove_count = remove_duplicate_fragments(fragments) #remove any duplicates
        print("A total of "+str(remove_count)+" duplicate fragments were removed from the initial grand total")
        for frag in fragments:
            output += frag.get_csv_str(e_mass_convert_dict,index_convert_dict)

    print("Saving Fragment Results")
    try:
        with open(output_directory,"w") as f:
            f.write(output)
    except FileNotFoundError:
        print("Could not write file: "+output_directory)
    return output
         
def save_predictions(fragments,e_mass_convert_dict,index_convert_dict,output_directory = ".//FragResults.csv"):
    print("\nSaving Fragment Prediction Results")
    output = "Frag Seq\tPredicted m/z\tLog\n"
    for frag in fragments:
        output += frag.get_csv_str(e_mass_convert_dict,index_convert_dict)
    try:
        with open(output_directory,"w") as f:
            f.write(output)
    except FileNotFoundError:
        print("Could not write file: "+output_directory)   

def remove_negative_fragments(total_fragments,element_mass_conversion_dict,index_element_conversion_dict):
    final_fragments = []
    remove_count = 0
    for f in total_fragments: #filter out any negative mass fragments since those are not possible
        if f.get_monoisotopic_mass(element_mass_conversion_dict,index_element_conversion_dict) > 0:
            final_fragments.append(f)
        else:
            remove_count += 1
            
    print("Removed "+str(remove_count)+" impossible fragment(s) (negative mass)")
    return final_fragments
    
def load_aa_data():
    Amino_Acid_Data = {}
    with open("AA_data.csv","r") as f:
        for line in f:
            if "name,formula" in line.lower():
                pass
            else:                                          #example data line from file = 'A,3|7|1|2|0|0,nonpolar,0|0|0|0|0|0~2|3|1|2|0|0'
                comma_data = line.split("\n")[0].split(",") #using example above should give list like this ['A','3|7|1|2|0|0','nonpolar','0|0|0|0|0|0~2|3|1|2|0|0']
                name = comma_data[0]
                #print("name = ",name)
                formula_data = comma_data[1].split("|") #using example should give list like this: ['3','7','1','2','0','0']
                formula = []
                for value in formula_data:
                    formula.append(int(value)) #using example after loop list should be: [3,7,1,2,0,0]
                #print("formula = ",formula)
                class_def = comma_data[2]
                #print("class_def = ",class_def)
                fragment_lists = comma_data[3].split("~") #using example should give list like this: ['0|0|0|0|0|0','2|3|1|2|0|0']
                fragments = []
                for ls in fragment_lists:
                    temp_list = ls.split("|") #using example first iteration should give: ['0','0','0','0','0','0']
                    final_list = []            #using example second iteration should give: ['2','3','1','2','0','0']
                    for j in range(1,len(temp_list)):
                        final_list.append(int(temp_list[j])) #using example first iteration should give: [0,0,0,0,0,0]
                    fragments.append(("side_chain: "+temp_list[0],final_list))        #using example second iteration should give: [2,3,1,2,0,0]
                                        #final fragments list should look this: [[0,0,0,0,0,0],[2,3,1,2,0,0]]
               
                Amino_Acid_Data[name] = AminoAcid(name,formula,class_def,fragments)
    return Amino_Acid_Data
   
def display_warranty():
    print("\n15. Disclaimer of Warranty.\n"+
        "THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY\n"+
        "APPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT\n"+
        "HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM 'AS IS' WITHOUT WARRANTY\n"+
        "OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO,\n"+
        "THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR\n"+
        "PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM\n"+
        "IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF\n"+
        "ALL NECESSARY SERVICING, REPAIR OR CORRECTION.\n\n"+
        "16. Limitation of Liability.\n"
        "IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING\n"+
        "WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES AND/OR CONVEYS\n"+
        "THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY\n"+
        "GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE\n"+
        "USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF\n"+
        "DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD\n"+
        "PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS),\n"+
        "EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF\n"
        "SUCH DAMAGES.\n\n"
    )

def display_license():
    with open("license.txt","r") as f:
        for line in f:
            print(line)
   
def main():
    print("\npredictFragments  Copyright (C) 2023  Dominic McBrayer\n"+
        "This program comes with ABSOLUTELY NO WARRANTY; for details type 'show w'.\n"+
        "This is free software, and you are welcome to redistribute it\n"+
        "under certain conditions; type 'show c' for full details."
    )
    
    
    element_mass_conversion_dict = load_element_masses()
    index_element_conversion_dict = load_element_index()

    fragments = []
    output_directory = "./Output/"
    while True:
        print("\nChoose an option:")
        print("1) Generate Fragments from Sequence")
        print("2) Run Multi-Fragment Predictions from Sequence List")
        print("3) Exit")
        choice = input()
        if choice == "show w":
            display_warranty()
        
        elif choice == "show c":
            display_license()
        
        elif choice == "1":
            sequence = get_user_input()
            AA_obj = get_AA_object(sequence)
            peptide = PepFrag(AA_obj,form_peptide_bonds = True)
            
            fragments = peptide.fragment()
            fragments = remove_negative_fragments(fragments,element_mass_conversion_dict,index_element_conversion_dict)
            if len(peptide) == 1: #indicates single amino acid "peptide"
                fragments,remove_count = remove_duplicate_fragments(fragments) #remove any duplicates
                print("A total of "+str(remove_count)+" duplicate fragments were removed from the initial grand total")
            print("*"*15+"FRAGMENTATION LIST"+"*"*15)
            print("ORIGINAL PEPTIDE")
            print(peptide)
            print("\t\tMonomass = "+str(peptide.get_monoisotopic_mass(element_mass_conversion_dict,index_element_conversion_dict))+" "+str(peptide.get_formula()))
            print("\t\tM/z value = "+str(peptide.get_mz_value(element_mass_conversion_dict,index_element_conversion_dict)))
            for frag in fragments:
                print("*New Fragment*\n",frag)
                print("\t\tMonomass = "+str(frag.get_monoisotopic_mass(element_mass_conversion_dict,index_element_conversion_dict))+" "+str(frag.get_formula()))
                print("\t\tM/z value = "+str(frag.get_mz_value(element_mass_conversion_dict,index_element_conversion_dict)))        
        elif choice == "2":
            print("Choose seq list")
            seq_file = choose_file_directory()
            if seq_file != "":
                seq_list = load_seq_list(seq_file)
                run_multiple_predictions(seq_list,element_mass_conversion_dict,index_element_conversion_dict)
            else:
                print("File not recognized, aborting")
       
        elif choice == "3":
            return
        else:
            print("Please select from options 1-3")

if __name__ == "__main__":
    main()
