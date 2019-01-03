#! /bin/bash

#####
# Load qctool module here

# (Broad institute)
source /broad/software/scripts/useuse
use .qctool-2.0-rc2

QCTOOL=qctool
#####

BGEN_PATH_BASE=$1
OUTPUT_DIR=$2

for CHROM in `seq 1 22`; 
do 
    echo "chrom$CHROM:"

    BGEN_PATH=`echo $BGEN_PATH_BASE | sed "s/\#/$CHROM/"`
    QCTOOL_FILE="$OUTPUT_DIR/chrom$CHROM.qctool_snpstats.txt"
    MFI_FILE="$OUTPUT_DIR/chrom$CHROM.mfi.txt"
    
    echo "    bgen input: $BGEN_PATH"

    echo "    Running qctool..."
    
    $QCTOOL -g $BGEN_PATH -osnp $QCTOOL_FILE -snp-stats

    echo "    Capturing the following columns:"
    cat $QCTOOL_FILE |grep -v "^#"  | head -n 1 | cut -d' ' -f 1-2,4-6,14,15,17 | sed 's/ /\t/g'

    echo "    SNP info output: $MFI_FILE"
    cat $QCTOOL_FILE |grep -v "^#"  | tail -n +2 | cut -d' ' -f 1-2,4-6,14,15,17 | sed 's/ /\t/g' > $MFI_FILE


    echo "    OK"

done

echo "Done"
