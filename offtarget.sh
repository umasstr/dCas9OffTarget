#!/bin/bash
o=output
i=null
m=8
#

if [ "$1" == "-h" ] || [ "$1" == "--help" ] || [ "$1" == "" ]; then
  echo "
		-b	<BED file>
		-o	<output dir ; default = output>
		-r	<reference genome used for alignment>
		-i	<input/background/IgG BED>
		-m	<# mismatches permitted; default = 8>
		-s	<sgRNA sequence ; use N for any nucleodide, does not count against mismatches>
		-B  <BAM>
		-I <Input BAM>"
  exit 0
fi

while getopts b:o:r:i:m:s:B: options; do
	case $options in 
		b) b=$OPTARG;; #bedfile
		o) o=$OPTARG;; #output dir
		r) r=$OPTARG;; #reference genome
		i) i=$OPTARG;; #input/background/IgG
		m) m=$OPTARG;; #mismatches
		s) s=$OPTARG;; #sgrna sequence to fasta
		B) B=$OPTARG;; #BAM

	esac
done

mkdir $o

echo ">sgRNA" > $o/sgRNA.fa && echo $s >> $o/sgRNA.fa

if [ "$B" == "null" ]; then
	echo "no BAM detected"
	exit 0
fi

if [ "$I" == "null" ]; then
	echo "no Input BAM detected"
	if [ "$i" == "null" ]; then
		echo -n "No background input entered. "
		sleep 1
		if [ "$b" == "null" ]; then
			macs2 callpeak -t $B -f BAM -g hs -o ./ -n ${B%.bam}
			mv "${B%.bam}"_peaks.narrowPeak "${B%.bam}".bed
			b="${B%.bam}".bed
		fi
		echo -n "Processing"
		for x in {1..5}; do echo -n "." && sleep 1;done
		echo ""
		bedtools getfasta -bed $b -fi $r > $o/${b%.bed*}.fa
	else
		if [ "$b" == "null" ]; then
			macs2 callpeak -t $B -f BAM -g hs -o ./ -n ${B%.bam}
			mv "${B%.bam}"_peaks.narrowPeak "${B%.bam}".bed
			b="${B%.bam}".bed
		fi
		echo -n "Subtracting background"
		for x in {1..5}; do echo -n "." && sleep 1;done
		echo ""
		bedtools intersect -v -a $b -b $i > $o/${b%.bed}.bg.bed
		echo -n "Processing"
		for x in {1..5}; do echo -n "." && sleep 1;done
		echo ""
		bedtools getfasta -bed $o/${b%.bed}.bg.bed -fi $r > $o/${b%.bed}.fa
		
	fi
else
	macs2 callpeak -t $B -c $I -f BAM -g hs -o ./ -n ${B%.bam}
	mv "${B%.bam}"_peaks.narrowPeak "${B%.bam}".bed
	b="${B%.bam}".bed
	bedtools getfasta -bed $b -fi $r > $o/${b%.bed*}.fa
fi
bedtools getfasta -bed $b -fi $r > $o/${b%.bed}.fa

patman -e $m -D $o/${b%.bed}.fa -P $o/sgRNA.fa -a -o $o/${b%.bed}.pat

paste -d'\t'  <(sed 's/:/\t/' $o/${b%.bed}.pat | sed 's/-/\t/' | sed 's/\t\t/\t/g' | awk '{x=$2+$4-1;y=$2+$5;print$1"\t"x"\t"y"\t"$1":"$2"-"$3"\t"$8"\t"$7}') <(cut -f5 $o/${b%.bed}.pat) | sed 's/\t\t/\t/g' | grep -v chrUn | grep -v chrM | grep -v random |sort -u -k4 | sort -k5,5n > $o/${b%.bed}.pat.bed 

bedtools getfasta -s -fi $r -bed $o/${b%.*}.pat.bed > $o/${b%.bed}.out.fasta

paste -d'\t' <(cut -f1-3 $o/${b%.*}.pat.bed) <( cut -f4-7 <(bedtools multicov -bams $B -bed <(paste -d'\t' <(cut -f4 $o/${b%.*}.pat.bed |sed 's/:/\t/' | sed 's/-/\t/' | sed 's/\t\t/\t/g') <(cut -f4-6 $o/${b%.*}.pat.bed) ) ) ) <(grep -vF '>' $o/${b%.bed}.out.fasta) | sort -u -k4,4 | awk '{ if ($7 < 100) {++removed;next} else {print $0} }' | sort -nr -k7,7 > $o/${b%.*}.final.bed

