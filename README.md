# INdigestion
-------------
"I Need a digestion" prints expected gel band sizes from restriction digests. The user can specify which enzymes they have on-hand and what band sizes they prefer.

When entering a list of files, separate by commas and no spaces between filenames. 

This program works on circular plasmids only. 

This program assumes that you are digesting a vector with a gene insert. It selects only enzymes or enzyme pairs that cuts in this insert. You must have the gene insert annotated as "gene" in Ape (a plasmid editor) or other similar program that saves to genbank format. 

![indigestion_input_output.jpeg](https://raw.githubusercontent.com/amcrabtree/INdigestion/master/images/indigestion_input_output.jpeg)

Band size parameters can be changed within the python script:
![indigestion_script.jpeg](https://raw.githubusercontent.com/amcrabtree/INdigestion/master/images/indigestion_script.jpeg)

Requirements:
- Biopython
- Reference files
   * my_enzymes.csv - <i>a csv containing the user's enzymes on-hand (to be modified by user)</i>
   * enzyme_dictionary.csv - <i>a csv containing reference information for enzyme names and cutting patterns</i>
