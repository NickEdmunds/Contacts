#!/usr/bin/env python2

# Author:  Nick Edmunds
# Program: Sort_contacts_v1_1.py
# Function: Calculate the contacts from a model pdb file (dimer) and them compare them to those from a Voronota-contacts run.
# Use Euclidean distance of 8A as specified in CASP15 literature.
# Score the contacts: 
# Short side-chains: A distance of 5.7 scores 0.8. Each 0.1 either side of that score (0.8 - 0.026). So 3.4 should be 0.2 and less than 3.3 will be set 0.1.
# The idea is that 5.7 is close to the perfect contact distance and much closer than this and it could just be an accidental clash.
# Same for the distances above 5.7 (-0.026 each 0.1) until 0.80 = 0.2. There will be no conatcs greater than 0.80. The idea is that these are getting less perfect so score less well.
# Long side-chains: A distance of 8.0 scores 0.8. Each 0.1 less than that distance will go down by 0.01. So score =distance/10. With long side-chains the likelihood of a close contact is less.
# Hydrophobic side-chains: A distance of 4.0 scores 0.8 as these tend to stack. Each 0.1 above that distance (0.8 - 0.015), so 8.0 should be 0.2.
# For distances below 4.0 (-0.03 each 0.1) until 2.0 = 0.2. Less than 2.0 scores 0.1.
# For mixed chain length/hydrophobicity - treat as short chains.
# For the Voro scores - literally if it's in my list and the Voro list then add 0.4 to the score with a maximum of 1.0.
# (CDA score - Contact Distance Agreement)

# Call functions and create empty lists to store Chain, Residue and Atom co-ordinates.
import numpy
import pandas as pd
from scipy.spatial import distance
Sort_contacts=[]
Sort_Voro=[]
Sort_Final=[]
Last_list1=[]
Last_list2=[]

# Sub-routine to select and choose the Target file
file_choice=raw_input('Please choose a Target for the _Contacts.txt and _Voro_CB_contacts.txt files: ')

# Open the output contacts file and the T1078_DeepHomo_contacts.txt file and order them both by residue number.
with open(file_choice+"_Contacts.txt","r+") as Cont, open (file_choice+"_Voro_CB_contacts.txt", "r+") as Voro:
    # List through the Atom_coords1 list
    for line in Cont:
        # Split the line and turn the residue numbers to integers
        sp=line.split()
        Res1=int(sp[0])
        ResA=sp[1]
        Res2=int(sp[3])
        ResB=sp[4]
        Dist=float(sp[6])
# For short side-chains
        if (ResA=='CYS' or ResA=='SER' or ResA=='GLY' or ResA=='THR' or ResA=='HIS' or ResA=='ASN' or ResA=='ASP') and (ResB=='CYS' or ResB=='SER' or ResB=='GLY' or ResB=='THR' or  ResB=='HIS' or ResB=='ASN' or ResB=='ASP') and Dist >= 5.7:
	    Score= round(0.8-((Dist-5.7)*0.26),2)
    	elif (ResA=='CYS' or ResA=='SER' or ResA=='GLY' or ResA=='THR' or ResA=='HIS' or ResA=='ASN' or ResA=='ASP') and (ResB=='CYS' or ResB=='SER' or ResB=='GLY' or ResB=='THR' or  ResB=='HIS' or ResB=='ASN' or ResB=='ASP') and Dist < 5.7 and Dist >=3.4:
	    Score= round(0.8-((5.7-Dist)*0.26),2)
    	elif (ResA=='CYS' or ResA=='SER' or ResA=='GLY' or ResA=='THR' or ResA=='HIS' or ResA=='ASN' or ResA=='ASP') and (ResB=='CYS' or ResB=='SER' or ResB=='GLY' or ResB=='THR' or  ResB=='HIS' or ResB=='ASN' or ResB=='ASP') and Dist <3.4:
	    Score= 0.1
# For short long-chains	    
	elif (ResA=='GLU' or ResA=='GLN' or ResA=='TRP' or ResA=='LYS' or ResA=='ARG') and (ResB=='GLU' or ResB=='GLN' or ResB=='TRP' or ResB=='LYS' or  ResB=='ARG'):
	    Score= round((Dist/10),2)
# For Hydrophobic side-chains	    
	elif (ResA=='PHE' or ResA=='TYR' or ResA=='VAL' or ResA=='ILE' or ResA=='ALA' or ResA=='PRO' or ResA=='MET' or ResA=='LEU') and (ResB=='PHE' or ResB=='TYR' or ResB=='VAL' or ResB=='ILE' or ResB=='ALA' or ResB=='PRO' or ResB=='MET' or ResB=='LEU') and Dist >= 4.0:
	    Score= round(0.8-((Dist-4.0)*0.15),2)
	elif (ResA=='PHE' or ResA=='TYR' or ResA=='VAL' or ResA=='ILE' or ResA=='ALA' or ResA=='PRO' or ResA=='MET' or ResA=='LEU') and (ResB=='PHE' or ResB=='TYR' or ResB=='VAL' or ResB=='ILE' or ResB=='ALA' or ResB=='PRO' or ResB=='MET' or ResB=='LEU') and Dist < 4.0 and Dist >=2.0:
	    Score= round(0.8-((4.0-Dist)*0.3),2)
	elif (ResA=='PHE' or ResA=='TYR' or ResA=='VAL' or ResA=='ILE' or ResA=='ALA' or ResA=='PRO' or ResA=='MET' or ResA=='LEU') and (ResB=='PHE' or ResB=='TYR' or ResB=='VAL' or ResB=='ILE' or ResB=='ALA' or ResB=='PRO' or ResB=='MET' or ResB=='LEU') and Dist < 2.0:
	    Score= 0.1
# For Mixed side-chains
	else:
	    if Dist >= 5.7:
	        Score= round(0.8-((Dist-5.7)*0.26),2)
            if Dist < 5.7 and Dist >=3.4:
	        Score= round(0.8-((5.7-Dist)*0.26),2)
	    if Dist <3.4:
	        Score= 0.1
 # Add the new lines to the list. Res1 and 2 are numerical forms of sp[0] and sp[2] so we're adding residue number and name from chain A and B and then the distance and score.
	Sort_contacts.append((Res1, ResA, Res2, ResB, Dist, Score))
    for line in Voro:
        # Split the line and turn the residue numbers to integers and the score to a float  A 6 42 PRO CB B 120 3536 PRO CB 1.8078 4.16518
        sp=line.split()
        ResN1=int(sp[1])
        ResN2=int(sp[6])
        DistN=round(float(sp[11]),1)
        ScorN=0.4
        # Add the new lines to the list. Res 1 and 2 are again numerical so add Res number, name, res number, name and distance. Set the score to 0.4
        Sort_Voro.append((ResN1,  sp[3], ResN2, sp[8], DistN, ScorN))
# Order the residues in the list by residue 1 and then residue 2
Sort_contacts.sort(key=lambda tup: (tup[0], tup[2]))
# Order the residues in the list by residue 1 and then residue 2
Sort_Voro.sort(key=lambda tup: (tup[0], tup[2]))
print "See the contents of the list: Sort_contacts - contacts from T1078o.pdb_Contacts.txt"
print Sort_contacts
print len(Sort_contacts)
print "See the contents of the list: Sort_Voro - contacts from Voro_CB_contacts"
print Sort_Voro
print len(Sort_Voro)
# Now compare the two lists and update Sort_contacts with the data in Sort_Voro
from collections import OrderedDict
# Initialise two ordered Dictionaries from Sort_contacts and Sort_Voro lists
SC= OrderedDict((v[:-1], v[-1]) for v in Sort_contacts)
SV= OrderedDict((v[:-1], v[-1]) for v in Sort_Voro)
# Update SF dictionary as SC with the values from SV summed where matched
SF={k: SC.get(k, 0) + SV.get(k, 0) for k in set(SC) | set(SV)}
# Turn the results back into a new final list called Sort_Final
Sort_Final=list(SF.items())
print Sort_Final
# Print the contents of Sort_Final to the output file.
with open(file_choice+"_Final_Contacts.txt","w") as Final:
#    for element in Sort_Final:
    Final.write('\n'.join('%s %s' % x for x in Sort_Final))
Final.close()
# Read in the contents of Final_Contacts.txt so that it can be formatted correctly and any scores > 1.0 are reset to 1.0
with open(file_choice+"_Final_Contacts.txt","r+") as Form:
    # List through the list
    for line in Form:
        # Split the line and turn the residue numbers to integers
        sp=line.split()
        ResA=sp[0]
        score=float(sp[5])
        if score >1.0:
            score=1.0
        else: 
            score=score
	Last_list1.append(("A"+ResA[1:len(ResA)], score))
	Last_list2.append(("B"+sp[2], score))
# Print the contents of Final_list to the FINAL output file.
with open(file_choice+"_CASP_Contacts_A.txt","w") as CASPA, open (file_choice+"_CASP_Contacts_B.txt", "w") as CASPB:
    CASPA.write('\n'.join('%s %s' % x for x in Last_list1))
    CASPB.write('\n'.join('%s %s' % x for x in Last_list2))
CASPA.close()
CASPB.close()
#End
