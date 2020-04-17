# INdigestion
What this program does:
"I Need a digestion" prints expected gel band sizes from restriction digests. The user can specify which enzymes they have on-hand and what band sizes they prefer.

&nbsp;&nbsp;Dependencies: [biopython](https://biopython.org/)
![biopython-logo](https://raw.githubusercontent.com/amcrabtree/synteny-mapper/master/images/biopython_logo_white.png)

<b>Files Required for Program to Run:</b>
1. Annotated plasmid file(s). Must be in genbank format. PLASMID MUST BE CIRCULARIZED. GENE INSERT MUST BE ANNOTATED AS "GENE". These requirements are easy to accomplish with the program Ape (a plasmid editor), which will save your output to the correct format for this program. Either that or manually adjust the text in the genbank file. Note that this program selects only enzymes or enzyme pairs that cuts at least once within the specified gene. 
2. Reference files:
* my_enzymes.csv - <i>a csv containing the user's enzymes on-hand (to be modified by user)</i>
* enzyme_dictionary.csv - <i>a csv containing reference information for enzyme names and cutting patterns.</i> Note that not all enzymes are on here yet. Add more enzymes as necessary. NEB has a nice list of recognition sites <a href="https://www.neb.com/tools-and-resources/selection-charts/alphabetized-list-of-recognition-specificities">here</a>. 
<p>&nbsp;</p>

<b>Adjusting Script & Reference Files Before Running: </b>
- Band size parameters can be changed within the python script:
![indigestion_script.jpeg](https://raw.githubusercontent.com/amcrabtree/INdigestion/master/images/indigestion_script.jpeg)
- Reference file "my_enzymes.csv" needs to be updated with your lab's list of enzymes. Make sure they are capitalized correctly too. They must match the names as printed in "enzyme_dictionary.csv" (case-sensitive). 
- Reference file "enzyme_dictionary.csv" must include all enzymes you have in your lab, though it can contain more. Ideally this file would have information for all restriction enzymes available, but I didn't want to add them all myself. There are hundreds. 
<p>&nbsp;</p>

<b>Program Input:</b>
- Full filename(s) of plasmid(s). When entering a list of files, separate by commas and no spaces between filenames. Include extension as part of the name. 
<p>&nbsp;</p>

<b>Program Output:</b>

![indigestion_input_output.jpeg](https://raw.githubusercontent.com/amcrabtree/INdigestion/master/images/indigestion_input_output.jpeg)
