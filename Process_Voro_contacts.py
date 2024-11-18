#!/usr/bin/env python2

# Author:  Nick Edmunds
# Program: Process_Voro_contacts_v1_1.py
# Function: To stript all the extraneous information out of the file Voro_expanded_contacts.txt which is made from a run of
# voronota-contacts -i T1078o.pdb --contacts-query '--no-same-chain --no-solvent' --sum-at-end > Voro_contacts.txt followed by:
# cat Voro_contacts.txt | voronota expand-descriptors | column -t > Voro_expanded_contacts.txt

# Call functions and create empty lists to store Chain, Residue and Atom co-ordinates.
import numpy
import pandas as pd
from scipy.spatial import distance
Voro_contacts=[]

# Sub-routine to select and choose the Target file
file_choice=raw_input('Please choose a Target for the _expanded_contacts.txt file: ')

# Open the output contacts file from the two Voro processes Voro_expanded_contacts.txt.
with open(file_choice+"_expanded_contacts.txt","r+") as Voro, open (file_choice+"_Voro_CB_contacts.txt", "w") as Voro_out:
    # List through the contacts list
    for line in Voro:
        # Split the line and select just: A:chain(0), res(1), atom(3), res_nm(5), atm_nm(6), B:chain(7), res(8), atom(10), res_nm(12), atom_nm(13), SA(14), dist(15)
        sp=line.split()
        # Add the lines to the list only where both chain A atom and chain B atom is CB.
	if sp[6]=='CB' and sp[13]=='CB':
	    Voro_contacts.append((sp[0], sp[1], sp[3], sp[5], sp[6], sp[7], sp[8], sp[10], sp[12], sp[13], sp[14], sp[15]))
# Print the contents of Voro_contacts to the output file.
    Voro_out.write('\n'.join('%s %s %s %s %s %s %s %s %s %s %s %s' % x for x in Voro_contacts))
print "See the contents of the list: Voro_contacts"
print Voro_contacts
print len(Voro_contacts)
Voro_out.close()
