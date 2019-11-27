#! /bin/bash 

###
# ADD CODES TO LOAD MODULES HERE
#
# Broad institute: 
# source /broad/software/scripts/useuse
# use .plink2-2.00a
# use OpenblasR
#
# erisone.partners.org
# module add plink/2.0a
# module add R-mkl/3.3.2
###


moddir=$1
valdir=$2
valtag=$3

if [ $# -lt 3 ]; then
    echo "Usage: nps_score.plink2.sh workdir valdir valdatasetID"
    exit
fi

# work dir
if [ ! -d "$moddir" ]; then
    echo "ERROR: NPS data directory does not exist: $moddir"
    exit 1
fi


if [ $# -gt 3 ]; then
    winshifts=( $@ )
    winshifts=("${winshifts[@]:3}")
else
    
    # auto-detect window shifts
    echo -n "Detecting window shifts..."

    numwinshifts=`find $moddir -name "win_*.*.*.Q.RDS" -exec basename '{}' \; | grep -o "^win_[0-9]*" | sort -u | sed 's/win_//' | wc -l`

    if [ $numwinshifts -eq 0 ]; then
	echo " ERROR: autodetect failed"
	exit 1
    else
	echo -n ": $numwinshifts shifts detected"
    fi

    winshifts=`find $moddir -name "win_*.*.*.Q.RDS" -exec basename '{}' \; | grep -o "^win_[0-9]*" | sort -u | sed 's/win_//'`

    echo -n " ("
    echo -n $winshifts
    echo ")"
fi

traindir=`Rscript -e "args <- readRDS(\"$moddir/args.RDS\"); cat(\":NPS:\\t\", args[[\"traindir\"]], sep='');" | grep -F ':NPS:' | cut -f2 `

traintag=`Rscript -e "args <- readRDS(\"$moddir/args.RDS\"); cat(\":NPS:\\t\", args[[\"traintag\"]], sep='');" | grep -F ':NPS:' | cut -f2 `

for winshift in $winshifts
do
    echo "----- Shifted by $winshift -----"

    for chrom in `seq 1 22`
    do

	modtag="$traintag.win_${winshift}"
	prefix="$valdir/chrom${chrom}.${valtag}"

	tail -n +2 $traindir/chrom$chrom.$traintag.snpinfo | paste - $moddir/$modtag.adjbetahat.chrom$chrom.txt | cut -f 3,6,7 > $moddir/$modtag.adjbetahat_plink2.chrom$chrom.txt

	if [ -f $prefix.bgen ]; then

	    if [ -f $prefix.sample ]; then

		echo "Calculating PRS for $prefix.bgen and $prefix.sample..."

		input="--bgen $prefix.bgen --sample $prefix.sample"
	
	    else
		echo "ERROR: $prefix.sample missing for $prefix.bgen"
		exit 1
	    fi

	elif [ -f $prefix.bed ]; then

	    echo "Calculating PRS for $prefix.bed..."
	
	    input="--bfile $prefix"

	elif [ -f $prefix.ped ]; then

	    echo "Calculating PRS for $prefix.ped..."
	
	    input="--file $prefix"

	else
	    echo "ERROR: cannot find $prefix/.bgen/.bed/.ped"
	    exit 1
	fi

	plink2 $input --score $moddir/$modtag.adjbetahat_plink2.chrom$chrom.txt 1 2 3 no-mean-imputation --out $moddir/$modtag.predY.$valtag.chrom$chrom --threads 1 --memory 4000
	
    done
done

echo "Done"
