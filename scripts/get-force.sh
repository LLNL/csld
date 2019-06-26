#!/bin/bash

if (( $# < 1 )) ; then
    echo -e "Usage: (single OUTCAR mode) $0 OUTCAR_file [tag (total,polaron)] [step]"
    echo -e "Usage: (batch mode) $0 -d dir1 dir2 ..."
    exit 0
fi

tag="TOTAL-FORCE"
dline=1 # one additional line of "-----"
Start=4
STEP="-1"

grep_occurence() {
    #    echo "Cannot obtain the $2 -th occurence"
    if [[ "$2" < 0 ]]; then # e.g. last line by default
	ngrep=$2
	ngrep=$((-ngrep))
	grep "$1" -n|tac|awk "NR==$ngrep"|sed 's/:.*//'
    elif [[ "$2" > 0 ]]; then
	grep "$1" -n| awk "NR==$2"|sed 's/:.*//'
    else
	exit -1
    fi
}


if [[ "$1" == "-d" ]]; then
    # batch mode
    batch=1
    shift
    outfs=( "$@" )
else
    batch=0
    outfs=($1)
fi

for f in "${outfs[@]}"; do
    if [[ "$batch" == 0 ]]; then
	outf=$f
	if [[ "$2" == "polaron" ]]; then
	    tag="Polaron Force"
	    dline=0
	    Start=2
	fi
	if [[ ! -z "$3" ]]; then
	    STEP="$3"
	fi
    else
	if [ ! -d "$f" ]; then
	    echo "WARNING: directory $f not found"
	    continue
	fi
	if [ -s "$f/OUTCAR" ]; then
	    outf="$f/OUTCAR"
	elif [ -s "%f/OUTCAR.bz2" ]; then
	    outf="$f/OUTCAR.bz2"
	elif [ -s "%f/OUTCAR.gz" ]; then
	    outf="$f/OUTCAR.gz"
	else
	    echo "WARNING: OUTCAR in directory $f not found"
	    continue
	fi
    fi

    if [[ "$outf" == *.bz2 ]]; then
	CAT=bzcat
    elif [[ "$outf" == *.gz ]]; then
	CAT=zcat
    else
	CAT=cat
    fi

    if [[ ! "$STEP" =~ ^-?[0-9]+$ ]]; then
	exit -1
    fi

    NIONS=`$CAT $outf|grep 'NIONS =' |head -n1|awk '{print $NF}'`
    if [ -z "$NIONS" ]; then
	exit -2
    fi

    line=`$CAT $outf | grep_occurence "$tag" $STEP`
    if [ -z "$line" ]; then
	exit -1
    fi
    line=$((line+dline))

    if [ $batch == 0 ]; then
	$CAT $outf| awk -v s=$Start "(NR>$line)&&(NR<=$line+$NIONS)"'{printf("%-16s %-16s %-16s\n", $(s), $(s+1), $(s+2))}'
    else
	$CAT $outf| awk -v s=$Start "(NR>$line)&&(NR<=$line+$NIONS)"'{printf("%-16s %-16s %-16s\n", $(s), $(s+1), $(s+2))}' > $f/force.txt
    fi
done
#$CAT $1| grep 'total drift' -B $NLINE |tail -n $NLINE|awk 'BEGIN{flag=0} {if(flag) { if (match($0, "total drift")) {exit} print $4, $5, $6;}  if (match($0, "TOTAL-FORCE")){ flag=1}}'  | grep '[0-9]'
#    $CAT $1| grep 'total drift' -B $NLINE |tail -n $NLINE|awk 'BEGIN{flag=1} {if(flag) { if (match($0, "total drift")) {exit} print $4, $5, $6;}  if (match($0, "TOTAL-FORCE")){ flag=1}}'  | grep '[0-9]'
