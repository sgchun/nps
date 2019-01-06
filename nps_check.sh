#! /bin/bash 

if [ $# -eq 0 ]; then
    echo "Usage: nps_check.sh nps_command ..." 
    echo "nps_command:"
    echo "    stdgt traindir traintag"
    echo "    init workdir"
    echo "    decor workdir winshift1 winshift2 ..."
    echo "    prune workdir winshift1 winshift2 ..."
    echo "    gwassig workdir winshift1 winshift2 ..."
    echo "    prep_part workdir winshift1 winshift2 ..."
    echo "    part workdir winshift1 winshift2 ..."
    echo "    weight workdir winshift1 winshift2 ..."
    echo "    back2snpeff workdir winshift1 winshift2 ..."
    echo "    score workdir valdir valtag winshift1 winshift2 ..."

    exit 1
fi

step=$1
status=0

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
	    echo "FAIL: .meandos missing"
	    status=1
	    continue
	fi

	if [ ! -s $filepre.snpinfo ]; then
	    echo "FAIL: .snpinfo missing"
	    status=1
	    continue
	fi

	if [ ! -s $filepre.stdgt.gz ]; then
	    echo "FAIL: .stdgt.gz missing"
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
    
    exit $status

elif [ $step == "init" ]; then

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
    
    ver=`Rscript -e "args <- readRDS(\"$workdir/args.RDS\"); cat(args[[\"VERSION\"]]);" | tail -n 1`
    
    echo "OK (version $ver)"

    echo -n "Checking $workdir/log ..."

    if [ ! -d $workdir/log ]; then
	echo "FAIL"
	exit 1
    fi

    echo "OK"

    # older log files
    outdated=`find $workdir/log/ -name "*.Rout.*" ! -newer "$workdir/args.RDS" | wc -l`
    if [ $outdated -gt 0 ]; then
	echo "WARNING: Potentially outdated log files in $workdir/log/"
	echo "$outdated Rout files are found older than $workdir/args.RDS"
    fi

elif [ $step == "decor" ] || [ $step == "prune" ] || [ $step == "gwassig" ] || [ $step == "part" ] || [ $step == "back2snpeff" ]; then
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
	
	    last=`tail -n 1 $logfile`

	    if [ "$last" != "Done" ]; then
		echo "FAIL (incomplete)"
		status=1
		continue
	    fi

	    echo "OK"
	done

	if [ $status != 0 ]; then
	    exit $status
	fi

	if [ $step == "prune" ]; then
	
	    echo -n "Checking window count..." 

	    if [ $winshift == 0 ]; then
		win1=`ls -l $workdir/win.*.Q.RDS | wc -l`
		win2=`ls -l $workdir/win.*.pruned.table | wc -l`

		if [ $win1 != $win2 ]; then
		    echo "FAIL ($win1 != $win2)"
		    exit 1
		else
		    echo "OK ($win1 windows)"
		fi

		echo -n "Checking timestamp..."
		
		decorfile=`ls -t $workdir/win.*.Q.RDS | head -n 1`
		outdated=`find $workdir/ -name "win.*.pruned.table" ! -newer "$decorfile" | wc -l`

		if [ $outdated != 0 ]; then
		    echo "FAIL (outdated pruning data)"
		    exit 1
		fi
		
		echo "OK"

	    else
		win1=`ls -l $workdir/win_$winshift.*.Q.RDS | wc -l`
		win2=`ls -l $workdir/win_$winshift.*.pruned.table | wc -l`

		if [ $win1 != $win2 ]; then
		    echo "FAIL ($win1 != $win2)"
		    exit 1
		else
		    echo "OK ($win1 windows)"
		fi

		echo -n "Checking timestamp..."

		decorfile=`ls -t $workdir/win_$winshift.*.Q.RDS | head -n 1`
		outdated=`find $workdir/ -name "win_$winshift.*.pruned.table" ! -newer "$decorfile" | wc -l`

		if [ $outdated != 0 ]; then
		    echo "FAIL (outdated pruning data)"
		    exit 1
		fi
		
		echo "OK"

	    fi

	elif [ $step == "gwassig" ]; then

	    # check timestamp
	    echo -n "Checking timestamp..."

	    if [ $winshift == 0 ]; then

		prunefile=`ls -t $workdir/win.*.pruned.table | head -n 1`
		outdated=`find $workdir/ -name "win.*.pruned.tailfix.table" ! -newer "$prunefile" | wc -l`

		if [ $outdated != 0 ]; then
		    echo "FAIL (outdated gwassig data)"
		    exit 1
		fi

		outdated=`find $workdir/ -name "tail_betahat.*.table" ! -newer "$prunefile" | wc -l`

		if [ $outdated != 0 ]; then
		    echo "FAIL (outdated gwassig data)"
		    exit 1
		fi

		outdated=`find $workdir/ -name "trPT.*.tail.RDS" ! -newer "$prunefile" | wc -l`

		if [ $outdated != 0 ]; then
		    echo "FAIL (outdated gwassig data)"
		    exit 1
		fi

	    else

		prunefile=`ls -t $workdir/win_$winshift.*.pruned.table | head -n 1`
		outdated=`find $workdir/ -name "win_$winshift.*.pruned.tailfix.table" ! -newer "$prunefile" | wc -l`

		if [ $outdated != 0 ]; then
		    echo "FAIL (outdated gwassig data)"
		    exit 1
		fi
	    fi

	    echo "OK"

	elif [ $step == "part" ]; then 

	    for chrom in `seq 1 22`
	    do 
		
		if [ $winshift == 0 ]; then 
		    trPT="$workdir/trPT.$chrom.RDS"
		else 
		    trPT="$workdir/win_$winshift.trPT.$chrom.RDS"
		fi
		
		echo -n "Checking $trPT ..."

		if [ ! -s $trPT ]; then
		    echo "FAIL (missing)"
		    status=1
		    continue
		fi

		dim=`Rscript -e "trPT <- readRDS(\"$trPT\"); cat(dim(trPT));" | tail -n 1`
		dim=`echo $dim | sed 's/ / x /g'`

		echo "OK ($dim)"
	    done

	    if [ $status != 0 ]; then
		exit $status
	    fi
	    
	    echo -n "Checking timestamp ..."
	    if [ $winshift == 0 ]; then 
		outdated=`find $workdir/ -name "trPT.*.RDS" ! -newer "$workdir/part.RDS" | grep -v tail.RDS | wc -l`
	    else 
		outdated=`find $workdir/ -name "win_$winshift.trPT.*.RDS" ! -newer "$workdir/win_$winshift.part.RDS" | grep -v tail.RDS | wc -l`
	    fi

	    if [ $outdated != 0 ]; then
		echo "FAIL (outdated trPT data)"
		exit 1
	    fi

	    echo "OK"

	elif [ $step == "back2snpeff" ]; then 

	    traintag=`Rscript -e "args <- readRDS(\"$workdir/args.RDS\"); cat(args[[\"traintag\"]]);" | tail -n 1`
	    traindir=`Rscript -e "args <- readRDS(\"$workdir/args.RDS\"); cat(args[[\"traindir\"]]);" | tail -n 1`

	    for chrom in `seq 1 22`
	    do 
		
		if [ $winshift == 0 ]; then 
		    snpeff="$workdir/$traintag.adjbetahat.chrom$chrom.txt"
		else 
		    snpeff="$workdir/$traintag.win_$winshift.adjbetahat.chrom$chrom.txt"
		fi
		
		echo -n "Checking $snpeff ..."

		if [ ! -s $snpeff ]; then
		    echo "FAIL (missing)"
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
		exit $status
	    fi
	    
	    echo -n "Checking timestamp ..."

	    if [ $winshift == 0 ]; then 

		outdated=`find $workdir/ -name "$traintag.adjbetahat.chrom*.txt" ! -newer "$workdir/PTwt.RDS" | wc -l`

	    else 

		outdated=`find $workdir/ -name "$traintag.win_$winshift.adjbetahat.chrom*.txt" ! -newer "$workdir/win_$winshift.PTwt.RDS" | wc -l`

	    fi

	    if [ $outdated != 0 ]; then
		echo "FAIL (outdated snpeff data)"
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

	if [ $winshift == 0 ]; then 
	    partfile="part.RDS"
	    prevfile=`ls -t $workdir/win.*.pruned.table $workdir/win.*.pruned.tailfix.table | head -n 1`

	else
	    partfile="win_$winshift.part.RDS"
	    prevfile=`ls -t $workdir/win_$winshift.*.pruned.table $workdir/win_$winshift.*.pruned.tailfix.table | head -n 1`
	    
	fi

	echo -n "Checking $workdir/$partfile ..."

	if [ ! -s $workdir/$partfile ]; then
	    echo "FAIL (missing)"
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

elif [ $step == "weight" ]; then
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

	if [ $winshift == 0 ]; then 

	    echo -n "Checking S0 weight ..."

	    if [ ! -s "$workdir/PTwt.tail.RDS" ]; then
		echo "FAIL (missing $workdir/PTwt.tail.RDS)"
		exit 1
	    fi

	    echo "OK"

	    ptwtfile="$workdir/PTwt.RDS"
	else
	    ptwtfile="$workdir/win_$winshift.PTwt.RDS"
	fi

	echo -n "Checking partition weights ..."
	if [ ! -s "$ptwtfile" ]; then
	    echo "FAIL (missing)"
	    exit 1
	fi
	
	dim=`Rscript -e "PTwt <- readRDS(\"$ptwtfile\"); cat(dim(PTwt));" | tail -n 1`
	dim=`echo $dim | sed 's/ / x /g'`

	echo "OK ($dim)"

	echo -n "Checking timestamp ..."

	if [ $winshift == 0 ]; then 
	    prevfile=`ls -t $workdir/trPT.*.RDS | head -n 1`
	    outdated=`find $workdir/ -name "PTwt*.RDS" ! -newer "$prevfile" | wc -l`
	else 
	    prevfile=`ls -t $workdir/win_$winshift.trPT.*.RDS | head -n 1`
	    outdated=`find $workdir/ -name "win_$winshift.PTwt.RDS" ! -newer "$prevfile" | wc -l`
	fi

	if [ $outdated != 0 ]; then
	    echo "FAIL (outdated PTwt data)"
	    exit 1
	fi

	echo "OK"
    done

elif [ $step == "score" ]; then
    echo "Verifying nps_$step:"


    if [ $# -lt 5 ]; then
	echo "Usage: nps_check.sh $step workdir valdir valtag winshift1 wishift2 ..."
	exit 1
    fi

    workdir=$2
    valdir=$3
    valtag=$4

    cmdargs=( $@ )
    argslen=${#cmdargs[@]}

    for (( k=4; k<argslen; k++ ))
    do
	winshift=${cmdargs[$k]}

	echo "----- Shifted by $winshift -----"


	traintag=`Rscript -e "args <- readRDS(\"$workdir/args.RDS\"); cat(args[[\"traintag\"]]);" | tail -n 1`

	if [ $winshift == 0 ]; then
	    modtag=$traintag
	else
	    modtag="$traintag.win_${winshift}"
	fi

	for chrom in `seq 1 22`
	do 

	    scorefile="$valdir/$modtag.predY.chrom$chrom.txt"

	    echo -n "Checking $scorefile ..."

	    if [ ! -s $scorefile ]; then
		echo "FAIL (missing)"
		status=1
		continue
	    fi
	    
	    # check line number
	    N=`zcat $valdir/chrom${chrom}.${valtag}.dosage.gz | head -n 1 | tr " " "\n" | tail -n +7 | wc -l`

	    N0=`cat $scorefile | wc -l`

	    if [ $N != $N0 ]; then
		echo "FAIL (incomplete)"
		status=1
		continue
	    fi

	    echo "OK (N=$N)"
	done

	if [ $status != 0 ]; then 
	    exit $status
	fi
	
	echo -n "Checking timestamp ..."

	prevfile=`ls -t $workdir/$modtag.adjbetahat.chrom*.txt | head -n 1`
	outdated=`find $valdir/ -name "$modtag.predY.chrom*.txt" ! -newer "$prevfile" | wc -l`

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

