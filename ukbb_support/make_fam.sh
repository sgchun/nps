#! /bin/bash 

WORK_DIR=$1
COHORT_NAME=$2

echo "Making $WORK_DIR/$COHORT_NAME.fam"

zcat $WORK_DIR/chrom1.$COHORT_NAME.QC2.dosage.gz | head -n 1 | tr ' ' "\n" | tail -n +7 | awk '{ print( $1 " " $1 " 0 0 0 -9" ) }' > $WORK_DIR/$COHORT_NAME.fam
