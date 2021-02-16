DIR=""
FASTQ_DIR=$DIR"fastq/"

COLLAPSE_SCRIPT="/path/to/collapse_identical_merged_reads_min_counts.py"

SEQUENCE_TO_TRIM_AT_END_OF_READS="AGATCGGAAGAGCGTCGTGTAGGGAAAGA"
BP_TO_TRIM_FROM_BEGINNING_OF_READS=20
OVERWRITE_FILES_IF_SAME_NAME=false
MIN_READS_FOR_COLLAPSING=2

TRIMMED_READS_DIR=$DIR"trimmed_fastq/"
mkdir $TRIMMED_READS_DIR

MERGED_READS_DIR=$DIR"merged_fastq/"
mkdir $MERGED_READS_DIR

COLLAPSED_READS_DIR=$DIR"collapsed_fastq/"
mkdir $COLLAPSED_READS_DIR

ALIGNMENTS_DIR=$DIR"alignments/"
mkdir $ALIGNMENTS_DIR

PROCESSED_COUNTS_DIR=$DIR"processed_counts/"
mkdir $PROCESSED_COUNTS_DIR

## trim reads for quality and TruSeq adapters
cd $FASTQ_DIR
for sample in *_1.fq.gz; 
do prefix=${sample%_1.fq.gz};
echo "trimming reads for sample:" $prefix;
READ1_FASTQ=$FASTQ_DIR$prefix"_1.fq.gz"
READ2_FASTQ=$FASTQ_DIR$prefix"_2.fq.gz"
echo $prefix; 
bbduk.sh in1=$READ1_FASTQ in2=$READ2_FASTQ out1=$TRIMMED_READ1_FASTQ out2=$TRIMMED_READ2_FASTQ literal=$SEQUENCE_TO_TRIM overwrite=$OVERWRITE_FILES_IF_SAME_NAME ftl=$BP_TO_TRIM_FROM_BEGINNING_OF_READS ktrim=r k=23 mink=11 hdist=1 tpe tbo; 
done

## merge reads
cd $TRIMMED_READS_DIR
for sample in *cleaned_R1.fastq.gz; do prefix=${sample%cleaned_R1.fastq.gz}; 
echo "merging reads for sample:" $prefix;
TRIMMED_READ1_FASTQ=$TRIMMED_READS_DIR$prefix"_trimmed_R1.fastq.gz"
TRIMMED_READ2_FASTQ=$TRIMMED_READS_DIR$prefix"_trimmed_R2.fastq.gz"
MERGED_FASTQ=$MERGED_READS_DIR$prefix"_cleaned_merged.fastq"
bbmerge.sh in1=$TRIMMED_READ1_FASTQ in2=$TRIMMED_READ2_FASTQ out=$MERGED_FASTQ overwrite=$OVERWRITE_FILES_IF_SAME_NAME; 
done

## collapse identical reads and change read name to read counts
cd $MERGED_READS_DIR
for sample in *merged.fastq; 
do prefix=${sample%cleaned_merged.fastq}; 
echo "collapsing reads for sample:" $prefix;
MERGED_FASTQ=$MERGED_READS_DIR$prefix"_cleaned_merged.fastq"
COLLAPSED_FASTQ=$COLLAPSED_READS_DIR$prefix"_collapsed.fastq"
python $COLLAPSE_SCRIPT $MERGED_FASTQ $COLLAPSED_FASTQ $MIN_READS_FOR_COLLAPSING
done
