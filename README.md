# INdigestion
What this program does:
"<u>I</u> <u>N</u>eed a <u>digestion</u>" ("INdigestion.py") is a python command-line app which prints expected gel band sizes from restriction digests. The user can specify which enzymes they have in their lab and what band sizes they prefer. For individuals new to command-line interface, there is also a GUI ("INdigestion_zenity.sh") which facilitates proper input of arguments and runs the command line app via a series of prompts. 
<p>&nbsp;</p>

<b>Dependencies</b>
* [python 3](https://www.python.org/downloads/), [biopython](https://biopython.org/), [zenity](https://linuxconfig.org/how-to-use-graphical-widgets-in-bash-scripts-with-zenity)

Downloading zenity
&nbsp;&nbsp;&nbsp;&nbsp;MacOS using [Homebrew](https://formulae.brew.sh/formula/zenity):
```
brew install zenity
```
&nbsp;&nbsp;&nbsp;&nbsp;Linux (or windows users using [Ubuntu](https://zoomadmin.com/HowToInstall/UbuntuPackage/zenity) app):
```
sudo apt-get update -y
sudo apt-get install -y zenity
```
<p>&nbsp;</p>

<b>CLI Usage</b>
```
python3 INdigestion.py -p my_plasmid.gb
```
<p>&nbsp;</p>

<b>Options</b>
flag | description
------------ | -------------
-p	[ARG]	| (required) plasmid file: <i>.gb .gbk .genbank .ape .fasta .fa </i>
-e	[ARG]	| file with list of your enzymes ("my_enzymes.csv"), if not in script folder
-n	[ARG]	| minimum number of desired bands
-x	[ARG]	| maximum number of desired bands
-s	[ARG]	| smallest size of desired bands
-b	[ARG]	| biggest size of desired bands
-g	[ARG]	| minimum gap size between bands
-i		| no insert in plasmid
-f		| plasmid in fasta format (no insert)
-t		| tests program with .gb test file
-h		| help (print options)
<p>&nbsp;</p>

<b>Files Required for Program to Run:</b>
1. Annotated plasmid file(s). Must be in genbank format. PLASMID MUST BE CIRCULARIZED. GENE INSERT MUST BE ANNOTATED AS "GENE". These requirements are easy to accomplish with the program Ape (a plasmid editor), which will save your output to the correct format for this program. Either that or manually adjust the text in the genbank file. Note that this program selects only enzymes or enzyme pairs that cuts at least once within the specified gene. 
2. Reference files:
   * my_enzymes.csv - <i>a csv containing the user's enzymes on-hand (to be modified by user)</i>
   * enzyme_dictionary.csv - <i>a csv containing reference information for enzyme names and cutting patterns.</i> Note that not all enzymes are on here yet. Add more enzymes as necessary. NEB has a nice list of recognition sites [here](https://www.neb.com/tools-and-resources/selection-charts/alphabetized-list-of-recognition-specificities). 
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
