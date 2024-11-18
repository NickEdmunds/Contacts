#!/usr/bin/env python2

# Author:  Nick Edmunds
# Program: Clean_up_contacts_v1_2.py
# Function: Import the 2 contacts A and B files output by Sort_contacts.py, make the lists unique and re-order and generally put the list into the required CASP format. As below:
# PFRMAT QA
# TARGET T1109
# AUTHOR ModFOLDdock
# METHOD ModFOLDdock: multiple clustering and single model methods are used for scoring global accuracy of multimers. Interface accuracy is scored using ModFOLDIA.
# MODEL 1
# QMODE 2
# ranked_21.pdb 0.9615 0.9616 A22:0.9725 A23:0.9832 A24:0.9914 A26:0.9674 A27:0.9932 A28:0.9804    | <- marke 100 chars.
# A30:0.9752 A31:0.9915 A34:0.9943 A47:0.9696 A49:0.9654 A56:0.9716 A60:0.9689 A61:0.9946 A62:0.9682 <- new line start with new record.
# A212:0.9856 A213:0.9786 A215:0.9740 A217:0.9723 A218:0.9766 B22:0.9981 B23:0.9985 B24:0.9976       <- chain B carry straight on from A
# 29-May-22 - Modified to output only unique residues from each chain. It does this by turning the list into an np.array and then keeping the residue entry with the highest score vale.
#	      It then converst the np.array back to a list so that it can be prined and formatted correctly.

# Call functions and create empty lists to store Chain, Residue and Atom co-ordinates.
import numpy as np
import pandas as pd
Sort_A =[]
int_A  =[]
Final_A=[]
Sort_B =[]
int_B  =[]
Final_B=[]

# Sub-routine to select and choose the Target file
file_choice=raw_input('Add the Target name: ')
top_model=raw_input('Add the name for the best model for this target: ')

# Open the output contacts A file .._CASP_contacts_A.txt.
with open(file_choice+"_CASP_Contacts_A.txt","r+") as contA:
    # List through the entries in contA
    for line in contA:
        # Split the line into Residue (Res) and Score (Scr) and turn the residue numbers to integers
        line=line.strip()
        sp=line.split(",")
        Res=sp[0]
        Scr=sp[1]
        # Add the records to a list along with the number part of Res.
        Sort_A.append((int(Res[1:4]), float(Scr.lstrip())))
        # Order the residues in the list by residue 1 and then residue 2
        Sort_A.sort()
Arr1 = np.array(Sort_A)
s=np.array(list((x, max(Arr1[Arr1[:,0]==x, 1])) for x in np.unique(Arr1[:,0])))
int_A = s.tolist()
Final_A=["A"+str(int(i[0]))+":"+str(i[1]) for i in int_A]

# Open the output contacts B file .._CASP_contacts_B.txt.
with open(file_choice+"_CASP_Contacts_B.txt","r+") as contB:
    # List through the entries in contB
    for line in contB:
        # Split the line into Residue (Res) and Score (Scr) and turn the residue numbers to integers
        line=line.strip()
        sp=line.split(",")
        Res=sp[0]
        Scr=sp[1]
        # Add the records to a list along with the number part of Res.
        Sort_B.append((int(Res[1:4]), float(Scr.lstrip())))
        # Order the residues in the list by residue 1 and then residue 2
        Sort_B.sort()
Arr2 = np.array(Sort_B)
t=np.array(list((x, max(Arr2[Arr2[:,0]==x, 1])) for x in np.unique(Arr2[:,0])))
int_B = t.tolist()
Final_B=["B"+str(int(i[0]))+":"+str(i[1]) for i in int_B]

# Print the contents of Final_A to the output file.
with open(file_choice+"_final_contact_list.txt","w") as Final:
     Final.write("PFRMAT QA" '\n')
     Final.write("TARGET "+file_choice+ '\n')
     Final.write("AUTHOR NE" '\n')
     Final.write("METHOD ModFOLDdock (regular, R and S) for multimer global (Score) and interface (QScore). Contacts.py for per residue Interface score." '\n')
     Final.write("MODEL 1" '\n')
     Final.write("QMODE 2" '\n')
     Final.write(""+top_model+" Add Score and QScore here" '\n')
     Final.write(" ".join(Final_A))
     Final.write(" ".join(Final_B))
Final.close()
#End
