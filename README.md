# INdigestion
What this program does:
"<u>I</u> <u>N</u>eed a <u>digestion</u>" ("INdigestion.py") is a python command-line app which prints expected gel band sizes from restriction digests. The user can specify which enzymes they have in their lab and what band sizes they prefer. For individuals preferring a graphic user interface, I'm working on a web app version of this using Plotly Dash. Stay tuned. 

<b>Dependencies</b>
* [python 3](https://www.python.org/downloads/), [biopython](https://biopython.org/)

<b>CLI Usage</b>
```
python3 INdigestion.py -p my_plasmid.gb
```

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
-i		| plasmid has a gene insert requiring an internal cut
-t		| tests program with associated test file
-h		| help (print options)
<p>&nbsp;</p>

<b>Files Required for Program to Run:</b>
1. <b>Annotated plasmid file.</b> Must be in genbank, ape, or fasta format. Program will only process circular plasmids. GENE INSERT MUST BE ANNOTATED AS "GENE". This can be done using [Ape (a plasmid editor)](https://jorgensen.biology.utah.edu/wayned/ape/) or a standard text editor. When the -i flag is selected, the program will select enzymes or enzyme pairs that cut at least once within the annotated gene. 
2. <b>CSV with a list of your enzymes.</b> If you store the file as "my_enzymes.csv" in the same folder as the script, you won't need to use the -e flag to specify where your enzyme file is located. 

<p>&nbsp;</p>

<b>Adjusting Script & Reference Files Before Running: </b>
- Band size parameters can be changed within the python script:
![indigestion_script.jpeg](https://raw.githubusercontent.com/amcrabtree/INdigestion/master/images/indigestion_script.jpeg)
- Reference file "my_enzymes.csv" needs to be updated with your lab's list of enzymes. 
<p>&nbsp;</p>

<b>Program Input</b>
- Plasmid file – Genbank or fasta file containing plasmid sequence
- Enzyme CSV file (“my_enzymes.csv”) – File contains the names of the enzymes you have on hand; update to match your stocks; doesn't need to be included in arguments unless you move it outside the file with the python script
<p>&nbsp;</p>

<b>Program Output</b>
- Sizes of expected bands print on command line screen. 
![indigestion_input_output.jpeg](https://raw.githubusercontent.com/amcrabtree/INdigestion/master/images/indigestion_input_output.jpeg)
<p>&nbsp;</p>

<b>Examples</b>

view all the options (“--help" = “-h”)
```
python INdigestion.py --help
```

makes sure program is working right (“--test" = “-t”)
```
python INdigestion.py --test	
```

minimum input
```
python INdigestion.py -p test_plasmid.gb	
```

cut within an insert required
```
python INdigestion.py -p test_plasmid.gb -i	
```

change minimum number of bands, cut within insert
```
python INdigestion.py -p test_plasmid.gb -i -n 4
```
