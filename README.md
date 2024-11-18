Run Contacts.py
---------------
To identify the first list of contacts from the top model identified by ModFOLDdock simply by calculating those interchain CB atoms < 8A apart.
Outputs a lis of A-B contacts with their distance.
$./Contacts.py
Please enter a full pdb filename as input: T1109_top_model.pdb
Please enter the chain idetifier for the first chain (A,C etc): A
Please enter the chain idetifier for the second chain (B,D etc): B
Output: <target>_Contacts.txt

Run VoroMQA contacts query * See below
--------------------------
To get the second list of contacts from VoroMQA.
$ voronota-contacts -i T1109_top_model.pdb --contacts-query '--no-same-chain --no-solvent' > T1109_voro_contacts.txt
$ cat T1109_voro_contacts.txt | voronota expand-descriptors | column -t > T1109_expanded_contacts.txt
--
Output: <target>_expanded_contacts.txt - This file is already supplied for H1106 and voroMQA may not be installed on your system.

Run Process_Voro_contacts.py
----------------------------
To stript all the extraneous information out of the file Voro_expanded_contacts.txt
$./Process_Voro_contacts.py
Please choose a Target for the _expanded_contacts.txt file: T1109 (automatically attaches "_expanded_contacts.txt")
Output: <target>_Voro_CB_contacts.txt

Run Sort_contacts.py
--------------------
Calculates a contacts score using <model>_Contacts.txt and then compares them to those in <target>_Voro_CB_contacts.txt
Score the contacts - see program header for calculations. Then add 0.4 to the score (maximum of 1.0) if contact is also in <target>_Voro_CB_contacts.txt.
$./Sort_contacts.py
Please choose a Target for the _Contacts.txt and _Voro_CB_contacts.txt files: T1109
Output: <target>_Final_Contacts.txt as intermediate file, then <target>_CASP_Contacts_A.txt and <target>_CASP_Contacts_B.txt

Run Clean_up_contacts.py
------------------------
To format the two A and B target lists into that required for submission to CASP. Need to create a unique A and B list from the input files, taking the highest score for each residue.
$ ./Clean_up_contacts.py
Add the Target name: H1106
Add the name for the best model for this target: H1106TS367_3.pdb
Input: <target>_CASP_Contacts_A.txt and <target>_CASP_Contacts_B.txt
Output: H1106_final_contact_list.txt (this is the formatted submission file)
