#!/usr/bin/env python2

# Author:  Nick Edmunds
# Program: Contacts_v1_1.py
# Function: Calculate the contacts from a model pdb file (dimer) and then compare them to those from a contact predictor.
# Use Euclidean distance - to start with use 8A as that was the distance cutoff used in ModFold for the CDA score.
# (CDA score - Contact Distance Agreement). The output produces a file file_choice+"_Contacts.txt where file_choice is the pdb input file name.
# This lists the residue number, residue name, atom from chain A then chain B plus the distance, e.g. 5 THR CB 119 SER CB 7.1
# The output is now limited to CB-CB < 8A.

# Call functions and create empty lists to store Chain, Residue and Atom co-ordinates.
import numpy
import pandas as pd
from scipy.spatial import distance
Atom_coords1 = []
Atom_coords2 = []
Res_coords1 = []
Res_coords2 = []
Full_coords=[]
Unique_coords=[]
Unique_res=[]
Sorted_contacts=[]
Sorted_Deep=[]

# Sub-routine to select and choose the pdb file
file_choice=raw_input('Please enter a full pdb filename as input: ')
fname, fextn = file_choice.split('.')
fname = fname[0:5]
chain1=raw_input('Please enter the chain idetifier for the first chain (A,C etc): ')
chain2=raw_input('Please enter the chain idetifier for the second chain (B,D etc): ')

# Open the pdb file in read mode and get chain1 coordinates. The variable indicating the file is pdb.
with open(file_choice,"r") as pdb:
    # Read each line of the pdb file
    for line in pdb:
        # Split the lines into individual parts
        sp=line.split()
        # Check whether each line is for ATOM and chain 1
        if sp[0]=='ATOM' and sp[4]==chain1 and (sp[2]=='CB' or (sp[3]=='GLY' and sp[2]=='CA')):
        # Append Atom#, Atom name, Residue, Chain ID, Residue#, and x,y,z coordinates to the Atom_coords1 list
            Atom_coords1.append((sp[1], sp[2], sp[3], sp[4], sp[5], float(sp[6]), float(sp[7]),float(sp[8])))
        # Append Residue, Chain ID, Residue# to the Res_coords1 list
            Res_coords1.append((sp[3], sp[4], sp[5]))
            Unique_coords1=list(set(Res_coords1))
            
# Open the pdb file again and get chain2 coordinates. The variable indicating the file is pdb2.
with open(file_choice,"r") as pdb2:
    # Read each line of the pdb file
    for line in pdb2:
        # Split the lines into individual parts
        sp=line.split()
        # Check whether each line is for ATOM and chain 2
        if sp[0]=='ATOM' and sp[4]==chain2 and (sp[2]=='CB' or (sp[3]=='GLY' and sp[2]=='CA')):
        # Append Atom#, Atom name, Residue, Chain ID, Residue#, and x,y,z coordinates to the Atom_coords2 list
            Atom_coords2.append((sp[1], sp[2], sp[3], sp[4], sp[5], float(sp[6]), float(sp[7]),float(sp[8])))
        # Append Residue, Chain ID, Residue# to the Res_coords2 list
            Res_coords2.append((sp[3], sp[4], sp[5]))
            Unique_coords2=set(Res_coords2)

# Print out the number of records found for Chain 1 and 2
print "There are " + str(len(Unique_coords1)) + " residues in Chain " + chain1
print "There are " + str(len(Unique_coords2)) + " residues in Chain " + chain2

# Open the output file and compare the coordinates of the two files
with open(fname+"_Contacts.txt","w") as Cont:
    # List through the Atom_coords1 list
    for i in range(len(Atom_coords1)):
        s=Atom_coords1[i]
        # Create a nested for loop to list through the Atom_coords2 list
        for j in range(len(Atom_coords2)):
            p=Atom_coords2[j]
            # Calculate the distance between the coordinates in the Atom_coords1 and Atom_coords2 lists
            dist = round(distance.euclidean([s[5], s[6], s[7]], [p[5], p[6], p[7]]),4)
            # Check whether the distance is less than 8A
            if dist <=8:
                Full_coords.append((s[4], s[2], s[1], p[4], p[2], p[1], round(dist,1)))
# Print the contents of Full_coords to the output file.
    Cont.write('\n'.join('%s %s %s %s %s %s %s' % x for x in Full_coords))
# Print out the total number of contacts found between Chain 1 and 2
print "There are " + str(len(Full_coords)) + " CB-CB atom contacts between chain " + chain1 + " and " + chain2 + ". These are listed in " + fname+"_Contacts.txt"
Cont.close()
# End
