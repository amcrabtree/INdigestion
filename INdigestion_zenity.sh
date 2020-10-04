#!/usr/bin/env sh

## alignment pipeline
## Authors: Angela Crabtree

################ CHECK PYTHON ##################
PY=$(python3 --version)
REFNAME=$(echo "${PY}" | cut -d" " -f2)
if [[ ${REFNAME} == "command" ]]
then 
	zenity --error \
		--width=400 --height=100 \
		--text "You don't have python3 installed\!\n\nInstall python3 (https://www.python.org/downloads/). "
	exit 0
fi

BIOPY=$(python3 -c "import Bio"; echo $?)
if [[ ${BIOPY} == "1" ]]
then
	zenity --error \
		--width=400 --height=100 \
		--text "You don't have biopython installed\!\n\nInstall biopython (https://biopython.org/wiki/Download). "
	exit 0
fi

################ START ASSEMBLING COMMAND LINE INPUT ##################
S_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
SCRIPT="${S_DIR}/INdigestion.py"
CLI=""
## choose plasmid file
FILE=$(zenity --file-selection --text "Choose your plasmid file")
EXTN=$(basename "${FILE}" | cut -d"." -f2)   # store name of extension
if [[ "${FILE}" == "" ]] 
then 
	exit 0 
elif [[ ${EXTN} != "fasta" ]] && [[ ${EXTN} != "fa" ]] && [[ ${EXTN} != "gb" ]] && [[ ${EXTN} != "gbk" ]] && [[ ${EXTN} != "genbank" ]] && [[ ${EXTN} != "ape" ]]
then 
	zenity --error --width=400 --height=200 \
		--text "Please select a plasmid file with one of the following extensions:\n\n.fasta \n.fa \n.gb \n.gbk \n.genbank \n.ape"
	exit 0
fi

if [[ ${EXTN} == "fasta" ]] || [[ ${EXTN} == "fa" ]]
then
	CLI="${CLI} -f"
else
	ANS_i=$(zenity --question --text "Does your file have an annotated insert?"; echo $?)
		if [[ ${ANS_i} == "1" ]]
		then
			CLI="${CLI} -i"
		fi
fi

## save a file-info-specific CLI variable so I can loop with different parameters
PRECLI="${CLI}"
LOOP="0" # initialize loop
while [ ${LOOP} == "0" ]
do
	################ CHANGING BAND SIZES ##################
	ANS=$(zenity --question --text "Do you want to use default parameters?"; echo $?)
	if [[ ${ANS} == "1" ]]
	then
		# minimum number of desired bands
		ANS_n=$(zenity --scale --text "What is the MINIMUM number of bands you want to see in your gel?" --min-value 1 --max-value 10 --value 2 --step 1)
		CLI="${PRECLI} -n ${ANS_n}"
	
		# maximum number of desired bands
		ANS_x=$(zenity --scale --text "What is the MAXIMUM number of bands you want to see in your gel?" --min-value 2 --max-value 10 --value 6 --step 1)
		CLI="${CLI} -x ${ANS_x}"
	
		# smallest size of desired bands
		ANS_s=$(zenity --scale --text "What is the SMALLEST band size you want to see in your gel?" --min-value 100 --max-value 1000 --value 500 --step 50)
		CLI="${CLI} -s ${ANS_s}"
	
		# biggest size of desired bands
		ANS_b=$(zenity --scale --text "What is the LARGEST band size you want to see in your gel?" --min-value 1000 --max-value 8000 --value 4000 --step 100)
		CLI="${CLI} -b ${ANS_b}"
	
		# minimum gap size between bands
		ANS_g=$(zenity --scale --text "What is the MINIMUM gap size you want to see in your gel?" --min-value 50 --max-value 500 --value 100 --step 50)
		CLI="${CLI} -g ${ANS_g}"
	fi

	################ RUN FULL COMMAND ##################
	
	echo "\npython3 \"${SCRIPT}\" -p \"${FILE}\" ${CLI}\n"
	python3 "${SCRIPT}" -p "${FILE}" ${CLI}

	################# OPTIONAL LOOP ##################
	LOOP=$(zenity --question --text "Do you want to try different parameters?"; echo $?)
done

exit 0