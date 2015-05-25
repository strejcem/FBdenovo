#!/bin/bash
# usage: refbot.sh sorted.fasta

# designed for sequences trimmed to 400bp

if [ -z "$1" ]
then
	echo "Usage: `basename $0` [fasta file]"
	echo "Input DNA sequences must be sorted by abundance (most = first)"
	exit
fi

FrameBot='java -jar /home/lovec/opt/FrameBot1.0/dist/FrameBot.jar'	#FrameBot.jar path
id=0.4	#FB identity cutoff
len=100	#FB Aminoacid length cutoff

	
fastastring()
{
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

FBdenovo()
{
	# removing results-files
	rm -f "$1"_denovo_prot.fasta
	rm -f "$1"_denovo_nucl.fasta
	rm -f "$1"_ref.STOP	# .STOP file contains unprocessed sequences that contain STOP codon - often non-target sequences


	mkdir FrameBot_$1
	cd FrameBot_$1
	
	# 0th cycle 
	i=0
	fastastring < ../$1 > "$1".string
	cp "$1".string STOP_check	# preparation for STOP codon check
	head -2 < STOP_check > temp	# translation of the most abundant sequence = "ref"
	$FrameBot translate 11 temp ref 1
	rm -f "$i"_ref.STOP
	while [[ $(tail -2 ref | grep "*" -c) > 0 ]]; do	#checks for STOP codon
			cat temp >> "$i"_ref.STOP
				tail -n+3 < STOP_check > STOP_check1	# skips the sequence
			mv STOP_check1 STOP_check
			head -2 < STOP_check > temp	# selects new reference
			$FrameBot translate 11 temp ref 1
	done

	# running global (seqs usually trimmed to the same length) FrameBot with the "ref" as a reference
	$FrameBot framebot -f -15 -N -a global -l "$len" -o "$i" -i "$id" ref "$1".string
	cat 0_corr_prot.fasta >> ../"$1"_denovo_prot.fasta
	cat 0_corr_nucl.fasta >> ../"$1"_denovo_nucl.fasta
	
	# 1+ cycles
	while [ -s "$i"_corr_nucl.fasta ] || [ -s ref ]; do
		cp "$i"_failed_nucl.fasta STOP_check	
		head -2 < STOP_check > temp	
		$FrameBot translate 11 temp ref 1
		rm -f "$i"_ref.STOP
		while [[ $(tail -2 ref | grep "*" -c) > 0 ]]; do 
			cat temp >> "$i"_ref.STOP
			tail -n+3 < STOP_check > STOP_check1
			mv STOP_check1 STOP_check
			head -2 < STOP_check > temp
			$FrameBot translate 11 temp ref 1
		done
	
		$FrameBot framebot -f -15 -N -a global -l "$len" -o "$[$i+1]" -i "$id" ref "$i"_failed_nucl.fasta
		cat "$[$i+1]"_corr_prot.fasta >> ../"$1"_denovo_prot.fasta
		cat "$[$i+1]"_corr_nucl.fasta >> ../"$1"_denovo_nucl.fasta
		cat "$i"_ref.STOP > ../"$1_"ref.STOP	# print unprocessed sequences that cotain STOP codon
		i=$[$i+1]
	done
echo "FrameBot de novo - finished"
echo "Sequences are NOT dereplicated"
}

FBdenovo $1
