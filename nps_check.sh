#! /bin/bash 

if [ $# -lt 1 ]; then
    echo "Usage:"
    echo "    nps_check.sh workdir"
    echo ""
    echo "To manually control nps_check, you can use the following options:"
    echo "    nps_check.sh stdgt traindir traintag"    
    echo "    nps_check.sh init workdir"
    echo "    nps_check.sh gwassig workdir"
    echo "    nps_check.sh decor workdir winshift1 [winshift2 ...]"
    echo "    nps_check.sh prune workdir winshift1 [winshift2 ...]"
    echo "    nps_check.sh prep_part workdir winshift1 [winshift2 ...]"
    echo "    nps_check.sh part workdir winshift1 [winshift2 ...]"
    echo "    nps_check.sh reweight workdir winshift1 [winshift2 ...]"
    echo "    nps_check.sh score workdir valtag winshift1 [winshift2 ...]"

    exit 1
fi

step=$1
status=0

# Check Rscript and R version
Rver=`Rscript -e 'cat(":NPS:\t", version$major, "\n", sep='')' | grep -F ':NPS:' | cut -f2`
Rver_string=`Rscript -e 'cat(":NPS:\t", version$version.string, "\n", sep='')' | grep -F ':NPS:' | cut -f2`

if [ $? != 0 ]; then 
   echo "ERROR: cannot run Rscript"
   exit 2
fi

if [ $Rver -lt 3 ]; then 
   echo "ERROR: R-3.0 or later is required: $Rver_string"
   exit 2
fi 

# Automatic mode
if [ $# -eq 1 ]; then
    if [ $step != "stdgt" ] && [ $step != "init" ] && [ $step != "decor" ] && [ $step != "prune" ] && [ $step != "prep_part" ] && [ $step != "part" ] && [ $step != "reweight" ] && [ $step != "score" ]; then

	workdir=$1

	echo "NPS data directory: $workdir"

	# work dir
	if [ ! -d "$workdir" ]; then
	    echo "ERROR: NPS data directory does not exist: $workdir"
	    exit 1
	fi

	# init
	auto=`ls -t $workdir/args.RDS $workdir/tail_betahat.*.table $workdir/*.Q.RDS $workdir/*.pruned.table $workdir/win_*.part.RDS $workdir/win_*.trPT.*.RDS $workdir/win_*.PTwt.RDS $workdir/*.win_*.predY.chrom*.txt 2> /dev/null | head -n 1 | grep -F args.RDS | wc -l `
	
	if [ $auto != 0 ]; then
	    ./nps_check.sh init $workdir
	    exit $?
	fi
	
	# auto-detect window shifts
	echo -n "Detecting window shifts..."

	numwinshifts=`find $workdir -name "win_*.*.*.Q.RDS" -exec basename '{}' \; | grep -o "^win_[0-9]*" | sort -u | sed 's/win_//' | wc -l`

	if [ $numwinshifts -eq 0 ]; then
	    echo " ERROR: autodetect failed"
	    exit 1
	else
	    echo -n ": $numwinshifts shifts detected"
	fi

	winshifts=`find $workdir -name "win_*.*.*.Q.RDS" -exec basename '{}' \; | grep -o "^win_[0-9]*" | sort -u | sed 's/win_//'`

	echo -n " ("
	echo -n $winshifts
	echo ")"


	# score
	auto=`ls -t $workdir/args.RDS $workdir/tail_betahat.*.table $workdir/*.Q.RDS $workdir/*.pruned.table $workdir/win_*.part.RDS $workdir/win_*.trPT.*.RDS $workdir/win_*.PTwt.RDS $workdir/*.win_*.predY.chrom*.txt 2> /dev/null | head -n 1 | grep -F .predY.chrom | wc -l `

	if [ $auto != 0 ]; then

	    # auto-detect valtag
	    echo -n "Detecting validation dataset tag..."
	    
	    scorefp=`ls -t $workdir/*.predY.chrom*.txt | head -n 1`
	    valtag=`basename $scorefp | sed 's/\.win_[0-9]*\.predY\.chrom[0-9]*\.txt//'`
	    if [ -z "$valtag" ]; then
		echo " ERROR: autodetect failed"
		exit 1
	    else 
		echo $valtag
	    fi

	    ./nps_check.sh score $workdir $valtag $winshifts
	    exit $?
	fi
    
	# reweight
	auto=`ls -t $workdir/args.RDS $workdir/tail_betahat.*.table $workdir/*.Q.RDS $workdir/*.pruned.table $workdir/win_*.part.RDS $workdir/win_*.trPT.*.RDS $workdir/win_*.PTwt.RDS $workdir/*.win_*.predY.chrom*.txt 2> /dev/null | head -n 1 | grep -F .PTwt. | wc -l `
    
	if [ $auto != 0 ]; then
	    ./nps_check.sh reweight $workdir $winshifts
	    exit $?
	fi

	# part
	auto=`ls -t $workdir/args.RDS $workdir/tail_betahat.*.table $workdir/*.Q.RDS $workdir/*.pruned.table $workdir/win_*.part.RDS $workdir/win_*.trPT.*.RDS $workdir/win_*.PTwt.RDS $workdir/*.win_*.predY.chrom*.txt 2> /dev/null | head -n 1 | grep -F .trPT. | wc -l `
    
	if [ $auto != 0 ]; then
	    ./nps_check.sh part $workdir $winshifts
	    exit $?
	fi

	# prep_part
	auto=`ls -t $workdir/args.RDS $workdir/tail_betahat.*.table $workdir/*.Q.RDS $workdir/*.pruned.table $workdir/win_*.part.RDS $workdir/win_*.trPT.*.RDS $workdir/win_*.PTwt.RDS $workdir/*.win_*.predY.chrom*.txt 2> /dev/null | head -n 1 | grep -F .part.RDS | wc -l `
	
	if [ $auto != 0 ]; then
	    ./nps_check.sh prep_part $workdir $winshifts
	    exit $?
	fi

	# prune
	auto=`ls -t $workdir/args.RDS $workdir/tail_betahat.*.table $workdir/*.Q.RDS $workdir/*.pruned.table $workdir/win_*.part.RDS $workdir/win_*.trPT.*.RDS $workdir/win_*.PTwt.RDS $workdir/*.win_*.predY.chrom*.txt 2> /dev/null | head -n 1 | grep -F .pruned.table | wc -l `
    
	if [ $auto != 0 ]; then
	    ./nps_check.sh decor $workdir $winshifts
	    
	    # always check decor as well
	    
	    if [ $? != 0 ]; then 
		exit $?
	    fi

	    ./nps_check.sh prune $workdir $winshifts

	    if [ $? != 0 ]; then 
		exit $?
	    fi
	
	    exit $?
	fi

	# decor
	auto=`ls -t $workdir/args.RDS $workdir/tail_betahat.*.table $workdir/*.Q.RDS $workdir/*.pruned.table $workdir/win_*.part.RDS $workdir/win_*.trPT.*.RDS $workdir/win_*.PTwt.RDS $workdir/*.win_*.predY.chrom*.txt 2> /dev/null | head -n 1 | grep -F .Q.RDS | wc -l `

	if [ $auto != 0 ]; then
	    ./nps_check.sh decor $workdir $winshifts
	    exit $?
	fi
	
	# gwassig
	auto=`ls -t $workdir/args.RDS $workdir/tail_betahat.*.table $workdir/*.Q.RDS $workdir/*.pruned.table $workdir/win_*.part.RDS $workdir/win_*.trPT.*.RDS $workdir/win_*.PTwt.RDS $workdir/*.win_*.predY.chrom*.txt 2> /dev/null | head -n 1 | grep -F tail_betahat. | wc -l `
    
	if [ $auto != 0 ]; then
	    ./nps_check.sh gwassig $workdir
	    exit $?
	fi
	
    fi

    echo "ERROR: cannot automatically figure out the previous step"
    exit 1
fi


# Manual mode 
if [ $step == "stdgt" ]; then
    echo "Verifying nps_stdgt:"
    
    if [ $# -ne 3 ]; then
	echo "Usage: nps_check.sh stdgt traindir traintag"
	exit 1
    fi
    
    traindir=$2
    traintag=$3
    
    for chrom in `seq 1 22`
    do 
	filepre="$traindir/chrom$chrom.$traintag"
	echo -n "Checking $filepre ..."

	if [ ! -s $filepre.meandos ]; then
	    echo "FAIL: .meandos missing or empty"
	    status=1
	    continue
	fi

	if [ ! -s $filepre.snpinfo ]; then
	    echo "FAIL: .snpinfo missing or empty"
	    status=1
	    continue
	fi

	if [ ! -s $filepre.stdgt.gz ]; then
	    echo "FAIL: .stdgt.gz missing or empty"
	    status=1
	    continue
	fi

	gzip -t $filepre.stdgt.gz 
	
	if [ $? != 0 ]; then 
	    echo "FAIL: .stdgt.gz broken"
	    status=1
	    continue
	fi

	echo "OK"
    done

    if [ $status != 0 ]; then 
	echo "FAILED"
    fi
    
    exit $status

elif [ $step == "init" ]; then
    echo "Verifying nps_$step:"

    if [ $# -ne 2 ]; then
	echo "Usage: nps_check.sh init workdir"
	exit 1
    fi

    workdir=$2

    echo -n "Checking $workdir/args.RDS ..."

    if [ ! -f $workdir/args.RDS ]; then
	echo "FAIL"
	exit 1
    fi
    
    ver=`Rscript -e "args <- readRDS(\"$workdir/args.RDS\"); cat(\":NPS:\\t\", args[[\"VERSION\"]], sep='');" | grep -F ':NPS:' | cut -f2 `
    
    echo "OK (version $ver)"

    echo -n "Checking $workdir/log ..."

    if [ ! -d $workdir/log ]; then
	echo "FAIL"
	exit 1
    fi

    echo "OK"

    # check stdgt
    traintag=`Rscript -e "args <- readRDS(\"$workdir/args.RDS\"); cat(\":NPS:\\t\", args[[\"traintag\"]], sep='');" | grep -F ':NPS:' | cut -f2 `
    traindir=`Rscript -e "args <- readRDS(\"$workdir/args.RDS\"); cat(\":NPS:\\t\", args[[\"traindir\"]], sep='');" | grep -F ':NPS:' | cut -f2 `

    ./nps_check.sh stdgt $traindir $traintag
    exit $?
fi

if [ $step == "gwassig" ]; then

    echo "Verifying nps_$step:"

    if [ $# -lt 2 ]; then
	echo "Usage: nps_check.sh $step workdir"
	exit 1
    fi

    workdir=$2

    for chrom in `seq 1 22`
    do 
	logfile="$workdir/log/nps_$step.Rout.$chrom"
	echo -n "Checking $logfile ..."
	
	if [ ! -f $logfile ]; then
	    echo "FAIL (missing)"
	    status=1
	    continue
	fi
	
	last=`grep -w Done $logfile | tail -n 1`

	if [ "$last" != "Done" ]; then
	    echo "FAIL (incomplete)"
	    status=1
	    continue
	fi

	echo "OK"
    done

    if [ $status != 0 ]; then
	echo "FAILED"
	exit $status
    fi

    traintag=`Rscript -e "args <- readRDS(\"$workdir/args.RDS\"); cat(\":NPS:\\t\", args[[\"traintag\"]], sep='');" | grep -F ':NPS:' | cut -f2 `
    traindir=`Rscript -e "args <- readRDS(\"$workdir/args.RDS\"); cat(\":NPS:\\t\", args[[\"traindir\"]], sep='');" | grep -F ':NPS:' | cut -f2 `

    # check tail_betahat files
    for chrom in `seq 1 22`
    do 
	tailbetahatfile="$workdir/tail_betahat.$chrom.table"
	
	echo -n "Checking $tailbetahatfile ..."
	
	if [ ! -s $tailbetahatfile ]; then
	    echo "FAIL (missing or empty)"
	    status=1
	    continue
	fi
	
	M1=`tail -n +2 $traindir/chrom$chrom.$traintag.snpinfo | wc -l`
	M2=`cat $tailbetahatfile | wc -l`
	
	if [ $M1 != $M2 ]; then
	    echo "FAIL (incomplete)"
	    status=1
	    continue
	fi
	
	echo "OK"
    done
		
    if [ $status != 0 ]; then
	echo "FAILED"
	exit $status
    fi

    # check summstat files
    for chrom in `seq 1 22`
    do 
	summstatfile="$workdir/harmonized.summstats.txt.$chrom"
	
	echo -n "Checking $summstatfile ..."
	
	if [ ! -s $summstatfile ]; then
	    echo "FAIL (missing or empty)"
	    status=1
	    continue
	fi
	
	M1=`tail -n +2 $traindir/chrom$chrom.$traintag.snpinfo | wc -l`
	M2=`tail -n +2 $summstatfile | wc -l`
	
	if [ $M1 != $M2 ]; then
	    echo "FAIL (incomplete)"
	    status=1
	    continue
	fi
	
	echo "OK"
    done
		
    if [ $status != 0 ]; then
	echo "FAILED"
	exit $status
    fi
    
    # check timestamp
    echo -n "Checking timestamp..."
    
    for chrom in `seq 1 22`
    do 
	outdated=`find $workdir/ -name "harmonized.summstats.txt.$chrom" ! -newer "$workdir/args.RDS" | wc -l`
	
	if [ $outdated != 0 ]; then
	    echo "FAIL (outdated gwassig data)"
	    exit 1
	fi
    done
    
    echo "OK"
    
elif [ $step == "decor" ] || [ $step == "prune" ] || [ $step == "part" ]; then
    echo "Verifying nps_$step:"

    if [ $# -lt 3 ]; then
	echo "Usage: nps_check.sh $step workdir winshift1 wishift2 ..."
	exit 1
    fi

    workdir=$2
    cmdargs=( $@ )
    argslen=${#cmdargs[@]}

    for (( k=2; k<argslen; k++ ))
    do
	winshift=${cmdargs[$k]}

	echo "----- Shifted by $winshift -----"

	for chrom in `seq 1 22`
	do 
	    logfile="$workdir/log/nps_$step.Rout.$winshift.$chrom"
	    echo -n "Checking $logfile ..."

	    if [ ! -f $logfile ]; then
		echo "FAIL (missing)"
		status=1
		continue
	    fi
	
	    last=`grep -w Done $logfile | tail -n 1`

	    if [ "$last" != "Done" ]; then
		echo "FAIL (incomplete)"
		status=1
		continue
	    fi

	    echo "OK"
	done

	if [ $status != 0 ]; then
	    echo "FAILED"
	    exit $status
	fi

	if [ $step == "prune" ]; then
	
	    echo -n "Checking window count..." 

	    win1=`ls -l $workdir/win_$winshift.*.Q.RDS | wc -l`
	    win2=`ls -l $workdir/win_$winshift.*.pruned.table | wc -l`

	    if [ $win1 != $win2 ]; then
		echo "FAIL ($win1 != $win2)"
		exit 1
	    else
		echo "OK ($win1 windows)"
	    fi

	    echo -n "Checking timestamp..."
	    
	    for chrom in `seq 1 22`
	    do 

		decorfile=`ls -t $workdir/win_$winshift.$chrom.*.Q.RDS | head -n 1`
		outdated=`find $workdir/ -name "win_$winshift.$chrom.*.pruned.table" ! -newer "$decorfile" | wc -l`

		if [ $outdated != 0 ]; then
		    echo "FAIL (outdated pruning data for chr$chrom)"
		    exit 1
		fi
	    done
	    
	    echo "OK"

	elif [ $step == "part" ]; then 

	    for chrom in `seq 1 22`
	    do 
		
		trPT="$workdir/win_$winshift.trPT.$chrom.RDS"
		
		echo -n "Checking $trPT ..."

		if [ ! -s $trPT ]; then
		    echo "FAIL (missing or empty)"
		    status=1
		    continue
		fi

		dim=`Rscript -e "trPT <- readRDS(\"$trPT\"); cat(\":NPS:\\t\", paste(dim(trPT), collapse=' x '), sep='');" | grep -F ':NPS:' | cut -f2 `

		echo "OK ($dim)"
	    done

	    if [ $status != 0 ]; then
		echo "FAILED"
		exit $status
	    fi
	    
	    echo -n "Checking timestamp ..."
	    
	    outdated=`find $workdir/ -name "win_$winshift.trPT.*.RDS" ! -newer "$workdir/win_$winshift.part.RDS" | grep -v tail.RDS | wc -l`

	    if [ $outdated != 0 ]; then
		echo "FAIL (outdated trPT data)"
		exit 1
	    fi

	    echo "OK"
	fi
    done

elif [ $step == "prep_part" ]; then
    echo "Verifying nps_$step:"

    if [ $# -lt 3 ]; then
	echo "Usage: nps_check.sh $step workdir winshift1 wishift2 ..."
	exit 1
    fi

    workdir=$2
    cmdargs=( $@ )
    argslen=${#cmdargs[@]}

    for (( k=2; k<argslen; k++ ))
    do
	winshift=${cmdargs[$k]}

	echo "----- Shifted by $winshift -----"

	partfile="win_$winshift.part.RDS"
	prevfile=`ls -t $workdir/win_$winshift.*.pruned.table | head -n 1`
	    
	echo -n "Checking $workdir/$partfile ..."

	if [ ! -s $workdir/$partfile ]; then
	    echo "FAIL (missing or empty)"
	    exit 1
	fi

	echo "OK"

	echo -n "Checking timestamp ..."

	outdated=`find $workdir/ -name "$partfile" ! -newer "$prevfile" | wc -l`    

	if [ $outdated != 0 ]; then
	    echo "FAIL (outdated partition files)"
	    exit 1
	fi

	echo "OK"
    done

elif [ $step == "reweight" ]; then
    echo "Verifying nps_$step:"

    if [ $# -lt 3 ]; then
	echo "Usage: nps_check.sh $step workdir winshift1 wishift2 ..."
	exit 1
    fi

    workdir=$2

    echo -n "Checking S0 weight ..."
    
    if [ ! -s "$workdir/PTwt.tail.RDS" ]; then
	echo "FAIL (missing or empty $workdir/PTwt.tail.RDS)"
	exit 1
    fi

    echo "OK"
    
    cmdargs=( $@ )
    argslen=${#cmdargs[@]}

    for (( k=2; k<argslen; k++ ))
    do
	winshift=${cmdargs[$k]}

	echo "----- Shifted by $winshift -----"

	ptwtfile="$workdir/win_$winshift.PTwt.RDS"
	
	echo -n "Checking partition weights ..."
	if [ ! -s "$ptwtfile" ]; then
	    echo "FAIL (missing or empty)"
	    exit 1
	fi
	
	dim=`Rscript -e "PTwt <- readRDS(\"$ptwtfile\"); cat(\":NPS:\\t\", paste(dim(PTwt), collapse=' x '), sep='')" | grep -F ':NPS:' | cut -f2 `

	echo "OK ($dim)"

	echo -n "Checking timestamp ..."

	prevfile=`ls -t $workdir/win_$winshift.trPT.*.RDS | head -n 1`
	outdated=`find $workdir/ -name "win_$winshift.PTwt.RDS" ! -newer "$prevfile" | wc -l`

	if [ $outdated != 0 ]; then
	    echo "FAIL (outdated PTwt data)"
	    exit 1
	fi

	echo "OK"

	# back2snpeff
	traintag=`Rscript -e "args <- readRDS(\"$workdir/args.RDS\"); cat(\":NPS:\\t\", args[[\"traintag\"]], sep='');" | grep -F ':NPS:' | cut -f2 `
	traindir=`Rscript -e "args <- readRDS(\"$workdir/args.RDS\"); cat(\":NPS:\\t\", args[[\"traindir\"]], sep='');" | grep -F ':NPS:' | cut -f2 `

	for chrom in `seq 1 22`
	do 
		
	    snpeff="$workdir/$traintag.win_$winshift.adjbetahat.chrom$chrom.txt"
	    
	    echo -n "Checking $snpeff ..."

	    if [ ! -s $snpeff ]; then
		echo "FAIL (missing or empty)"
		status=1
		continue
	    fi

	    M1=`tail -n +2 $traindir/chrom$chrom.$traintag.snpinfo | wc -l`
	    M2=`cat $snpeff | wc -l`

	    if [ $M1 != $M2 ]; then
		echo "FAIL (marker count mismatch: $M1 != $M2)"
		status=1
		continue
	    fi

	    echo "OK"
	done

	if [ $status != 0 ]; then 
	    echo "FAILED"
	    exit $status
	fi
	    
	echo -n "Checking timestamp ..."

	outdated=`find $workdir/ -name "$traintag.win_$winshift.adjbetahat.chrom*.txt" ! -newer "$workdir/win_$winshift.PTwt.RDS" | wc -l`


	if [ $outdated != 0 ]; then
	    echo "FAIL (outdated snpeff data)"
	    exit 1
	fi

	echo "OK"
	
    done

elif [ $step == "score" ]; then
    echo "Verifying nps_$step:"

    if [ $# -lt 4 ]; then
	echo "Usage: nps_check.sh $step workdir valtag winshift1 wishift2 ..."
	exit 1
    fi

    workdir=$2
    valtag=$3

    cmdargs=( $@ )
    argslen=${#cmdargs[@]}

    for (( k=3; k<argslen; k++ ))
    do
	winshift=${cmdargs[$k]}

	echo "----- Shifted by $winshift -----"

	traintag=`Rscript -e "args <- readRDS(\"$workdir/args.RDS\"); cat(\":NPS:\\t\", args[[\"traintag\"]], sep='');" | grep -F ':NPS:' | cut -f2 `

	modtag="$traintag.win_${winshift}"

	for chrom in `seq 1 22`
	do 

	    scorefile="$workdir/$modtag.predY.chrom$chrom.txt"

	    echo -n "Checking $scorefile ..."

	    if [ ! -s $scorefile ]; then
		echo "FAIL (missing or empty)"
		status=1
		continue
	    fi
	    
	    # check line number
	    if [ $chrom -ne "1" ]; then

		N0=`cat $scorefile | wc -l`

		N=`cat $workdir/$modtag.predY.chrom1.txt | wc -l`

		if [ $N != $N0 ]; then
		    echo "FAIL (incomplete)"
		    status=1
		    continue
		fi
	    fi

	    echo "OK (N=$N)"
	done

	if [ $status != 0 ]; then 
	    echo "FAILED"
	    exit $status
	fi
	
	echo -n "Checking timestamp ..."

	prevfile=`ls -t $workdir/$modtag.adjbetahat.chrom*.txt | head -n 1`
	outdated=`find $workdir/ -name "$modtag.predY.chrom*.txt" ! -newer "$prevfile" | wc -l`

	if [ $outdated != 0 ]; then
	    echo "FAIL (outdated score data)"
	    exit 1
	fi

	echo "OK"
    done

else 
    echo "ERROR: unknown NPS step: $step"
    exit 1
fi

exit $status

