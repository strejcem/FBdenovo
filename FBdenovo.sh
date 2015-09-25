#!/bin/bash

## v1.01
## FrameBot de novo bash script
## May 2015
## by Michal Strejcek (strejcem@vscht.cz)


# modify this to the full path to the FrameBot.jar
RDPframebot='java -jar /home/lovec/opt/RDPTools/FrameBot.jar'

fastastring() {
	#transform fasta file from block format to single-line format	
	while read line; do
        	if [ "${line:0:1}" = ">" ]; then
			if [ ! -z "$seq" ]; then
				echo "$seq"
			fi		
			unset seq
			echo "$line"
		else		
			seq="${seq}${line}"
		fi
	done
	echo "$seq"
}

FBdenovo() {

	echo "using $in as input file"
	echo "FB identity treshold set to $id"
	echo "FB min AA length is set to $len"
	# removing results-files
	rm -f "$in"_denovo_prot.fasta
	rm -f "$in"_denovo_nucl.fasta
	rm -f "$in"_ref.STOP	# .STOP file contains unprocessed sequences that contain STOP codon - often non-target sequences

	mkdir FrameBot_$in
	cd FrameBot_$in
	
	# 0th cycle 
	i=0
	fastastring < ../$in > "$in".string
	cp "$in".string STOP_check	# preparation for STOP codon check
	head -2 < STOP_check > temp	# translation of the most abundant sequence = "ref"
	$RDPframebot translate 11 temp ref 1
	rm -f "$i"_ref.STOP
	while [[ $(tail -2 ref | grep "*" -c) > 0 ]]; do	#checks for STOP codon
			cat temp >> "$i"_ref.STOP
				tail -n+3 < STOP_check > STOP_check1	# skips the sequence
			mv STOP_check1 STOP_check
			head -2 < STOP_check > temp	# selects new reference
			$RDPframebot translate 11 temp ref 1
	done

	# running global (seqs usually trimmed to the same length) FrameBot with the "ref" as a reference
	$RDPframebot framebot -f -15 -N -a global -l "$len" -o "$i" -i "$id" ref "$in".string
	cat 0_corr_prot.fasta >> ../"$in"_denovo_prot.fasta
	cat 0_corr_nucl.fasta >> ../"$in"_denovo_nucl.fasta
	
	# 1+ cycles
	while [ -s "$i"_corr_nucl.fasta ] || [ -s ref ]; do
		cp "$i"_failed_nucl.fasta STOP_check	
		head -2 < STOP_check > temp	
		$RDPframebot translate 11 temp ref 1
		rm -f "$i"_ref.STOP
		while [[ $(tail -2 ref | grep "*" -c) > 0 ]]; do 
			cat temp >> "$i"_ref.STOP
			tail -n+3 < STOP_check > STOP_check1
			mv STOP_check1 STOP_check
			head -2 < STOP_check > temp
			$RDPframebot translate 11 temp ref 1
		done
	
		$RDPframebot framebot -f -15 -N -a global -l "$len" -o "$[$i+1]" -i "$id" ref "$i"_failed_nucl.fasta
		cat "$[$i+1]"_corr_prot.fasta >> ../"$in"_denovo_prot.fasta
		cat "$[$i+1]"_corr_nucl.fasta >> ../"$in"_denovo_nucl.fasta
		if [ -f "$i"_ref.STOP ]; then
			cat "$i"_ref.STOP > ../"$in_"ref.STOP	# print unprocessed sequences that cotain STOP codon
		fi
		echo "iteration "$i""
		i=$[$i+1]
	done
	
	rm -f temp
	rm -f ref

	echo "FrameBot de novo - finished"
	echo "Sequences are NOT dereplicated"
}

usage="\n$(basename "$0") -i <fasta_file> [-t <0..1>] [-l <integer>]
Where:
	-h  show this message
	-i  input file (sorted fasta file; most = first)
	-t  FrameBot internal parameter for identity threshold (default=0.4)
	-l  FrameBot internal parameter for minimal AA length (default=100)

- requires installed FrameBot https://github.com/rdpstaff/Framebot and specify the full path to it in the beginning of this script
- the input must be sorted by sequence abundance (first is the most)
- the input sequences should be trimmed to the same lenght (e.g. 400bp)" 

in=
id=0.4
len=100
while getopts "hi:tl" OPTION; do
	case $OPTION in
	h)
		echo -e "$usage"
		exit
		;;
	i)
		in=$OPTARG
		;;
	t)
		id=$OPTARG
		;;
	l)
		len=$OPTARG     
		;;
	*)
		echo -e "$usage" >&2
		exit 1
		;;
	esac
done
if [ -z "$in" ]; then
	echo -e "\nillegal operation"
	echo -e "$usage" >&2
	exit 1
else
	if [ ! -f "$in" ];then
		echo -e "\nThe input file does not exist!"
		echo -e "$usage"
		exit 1
fi
fi

FBdenovo
