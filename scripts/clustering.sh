#!/bin/bash
# DESCRIPTION
#    Generation of consensus sequences to obtain
#    lowest possible error rate. Part of the longread-UMI-HIV-pipeline.
#    
# IMPLEMENTATION
#    author   Laurens Lambrechts     
#    license  GNU General Public License
# TODO


### Description ----------------------------------------------------------------

USAGE="
-- longread_umi clustering: Cluster related UMI consensus sequences.
   Clustering for UMI consensus sequences and creates a consensus sequence with >=3x coverage. 
   Reads are initially grouped by read clustering at 99.5% identity and a centroid sequence is picked.
   The centroid sequence is used as a mapping reference for all reads in the cluster.
   
usage: $(basename "$0" .sh) [-h -b] (-c file -F string -R string -m value -M value -p value -o dir -i dir -t value ) 

where:
    -h  Show this help text.
    -c  UMI consensus file.
    -F  Forward tagging primer sequence.
    -R  Reverse tagging primer sequence.
    -m  Minimum read length.
    -M  Maximum read length.
    -p  Percentage identity for usearch [Default = 0.995].
    -o  Output directory.
    -i  Input directory.
    -t  Number of threads to use. [Default = 1]
    -b  Debug flag. Keep temp files. [Default = NO]
"

### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hzc:F:R:m:M:p:i:o:t:b' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    c) CONSENSUS_FILE=$OPTARG;;
    F) FW2=$OPTARG;;
    R) RV2=$OPTARG;;
    m) MIN_LENGTH=$OPTARG;;
    M) MAX_LENGTH=$OPTARG;;
    p) IDEN=$OPTARG;;  
    o) OUT_DIR=$OPTARG;;
    i) IN_DIR=$OPTARG;;
    t) THREADS=$OPTARG;;
    b) DEBUG="YES";;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${CONSENSUS_FILE+x} ]; then echo "-c $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${OUT_DIR+x} ]; then echo "-o $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${IN_DIR+x} ]; then echo "-i $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${MIN_LENGTH+x} ]; then echo "-m $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${MAX_LENGTH+x} ]; then echo "-M $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${IDEN+x} ]; then echo "-p is missing. Defaulting to 99.5 identity."; IDEN=0.995; fi;
if [ -z ${FW2+x} ]; then echo "-F $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${RV2+x} ]; then echo "-R $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${THREADS+x} ]; then echo "-t is missing. Defaulting to 1 thread."; THREADS=1; fi;
if [ -z ${DEBUG+x} ]; then DEBUG="NO"; fi;

### Source commands and subscripts -------------------------------------
. $LONGREAD_UMI_PATH/scripts/dependencies.sh # Path to dependencies script

# Program paths

export SEQTK=seqtk
export GNUPARALLEL=parallel
export RACON=racon
export MINIMAP2=minimap2
export GAWK=gawk
export SAMTOOLS=samtools
export BCFTOOLS=bcftools
export CUTADAPT=cutadapt
export PORECHOP_UMI=porechop
export FILTLONG=filtlong
export BWA=bwa
export USEARCH=usearch

### Custom functions ----------------------------------------------------------
cluster_consensus () {
  # Input
  local IN=$(cat) # Single cluster sequence
  local VARIANT_OUT=$1 # Output dir for cluster
  local CLUSTER_READ_DIR=$2 # Path to cluster reads
  local CONSENSUS_FILE=$3 # Path to consensus reads
  local CLUSTER_THREADS=$4 # Number of threads allocated

  # Name format
  local CLUSTER_NAME=$(echo "$IN" |\
    $GAWK 'NR==1{gsub(">|;.*","",$0); print $0}')

  # Prepare output folder
  local CLUSTER_OUT=$VARIANT_OUT/${CLUSTER_NAME}
  mkdir $CLUSTER_OUT
  echo "$IN" > $CLUSTER_OUT/${CLUSTER_NAME}_con.fa

  # Map cluster reads to cluster sequences
  $MINIMAP2 \
    -t $CLUSTER_THREADS \
    -ax map-ont\
    $CLUSTER_OUT/${CLUSTER_NAME}_con.fa \
    $CLUSTER_READ_DIR/${CLUSTER_NAME} |\
    $SAMTOOLS view -b - \
    > $CLUSTER_OUT/temp.bam

  $SAMTOOLS sort \
    -o $CLUSTER_OUT/${CLUSTER_NAME}.bam \
    $CLUSTER_OUT/temp.bam
  $SAMTOOLS index $CLUSTER_OUT/${CLUSTER_NAME}.bam
  rm $CLUSTER_OUT/temp.bam   

    $BCFTOOLS mpileup \
    -Ov \
    -d 1000000 \
    -L 1000000 \
    -a "FORMAT/AD,FORMAT/DP" \
    -f $CLUSTER_OUT/${CLUSTER_NAME}_con.fa \
    $CLUSTER_OUT/${CLUSTER_NAME}.bam \
    > $CLUSTER_OUT/${CLUSTER_NAME}_pileup.vcf

 $BCFTOOLS norm \
    -f $CLUSTER_OUT/${CLUSTER_NAME}_con.fa \
    $CLUSTER_OUT/${CLUSTER_NAME}_pileup.vcf \
    -Ov \
    -o $CLUSTER_OUT/${CLUSTER_NAME}_norm.vcf

  $BCFTOOLS view -i 'AD[0:1]/FORMAT/DP>0.5' \
    -Oz \
    $CLUSTER_OUT/${CLUSTER_NAME}_norm.vcf \
    > $CLUSTER_OUT/${CLUSTER_NAME}_polish.vcf
  $BCFTOOLS index $CLUSTER_OUT/${CLUSTER_NAME}_polish.vcf

  $BCFTOOLS consensus $CLUSTER_OUT/${CLUSTER_NAME}_polish.vcf \
    -f $CLUSTER_OUT/${CLUSTER_NAME}_con.fa |\
    $SEQTK seq -l0 - |\
      $GAWK -v name="$CLUSTER_NAME" -v out="$CLUSTER_OUT" '
        /^>/{
          ++n;
          if (n == 1){
            split($0,size,";")
            gsub("size=", ";size=", size[2])
            print ">" name size[2] > out "/" name "_consensus.fa"
            getline
            print > out "/" name "_consensus.fa"
          }
        } END {
          if (n > 1){msg = "warning"} else {msg="ok"}
          print name, n, msg >> out "/bincon_log.txt"
        }'
}
export -f cluster_consensus

### Megabin Consensus pipeline -------------------------------------------------

# Threads handling
CLUSTER_JOBS=10
CLUSTER_THREADS=$(( $THREADS/$CLUSTER_JOBS ))
[ "$CLUSTER_THREADS" -gt 1 ] || CLUSTER_THREADS=1

# Format names
CONSENSUS_NAME=${CONSENSUS_FILE##*/}
CONSENSUS_NAME=${CONSENSUS_NAME%.*}

# Prepare output folders
if [ -d "$OUT_DIR" ]; then
  echo "Output folder exists. Exiting..."
  exit 0
fi

LOG_DIR=$OUT_DIR/logs

mkdir $OUT_DIR
mkdir $OUT_DIR/clusters


# Trim UMI consensus data
longread_umi trim_amplicon \
  -d $IN_DIR         `# Path to consensus data`\
  -p "$CONSENSUS_NAME.fa"    `# Consensus file pattern. Regex must be flanked by '"..."'`\
  -o $OUT_DIR          `# Output folder`\
  -F $FW2                 `# Forward primer sequence`\
  -R $RV2                 `# Reverse primer sequence`\
  -m $MIN_LENGTH          `# Min read length`\
  -M $MAX_LENGTH          `# Max read length` \
  -t $THREADS             `# Number of threads` \
  -l $OUT_DIR/logs 

mv $OUT_DIR/"$CONSENSUS_NAME.fa" $OUT_DIR/"$CONSENSUS_NAME".trimmed.fa

# Hard mask homopolymers
$SEQTK seq -l0 $OUT_DIR/"$CONSENSUS_NAME".trimmed.fa |\
$GAWK '
  /^>/{print}
  !/^>/{
    gsub(/A{3,}/, "AA", $0)
    gsub(/T{3,}/, "TT", $0)
    gsub(/C{3,}/, "CC", $0)
    gsub(/G{3,}/, "GG", $0)
    print $0
  }' > $OUT_DIR/m_temp.fa

# Usearch find uniques
$USEARCH \
  -fastx_uniques $OUT_DIR/m_temp.fa \
  -strand both \
  -fastaout $OUT_DIR/u_temp.fa  \
  -uc $OUT_DIR/u_temp.uc \
  -sizeout

# Usearch cluster @ 99.5 identity

$USEARCH  \
    -cluster_fast $OUT_DIR/u_temp.fa  \
    -id $IDEN \
    -strand both \
    -centroids $OUT_DIR/u_temp_995.fa \
    -uc $OUT_DIR/u_temp_995_masked.uc \
    -maxrejects 0 \
    -sizein \
    -sizeout \
    -sort length 


# Bin reads according to clusters
$GAWK \
  -v R_FA="$OUT_DIR/"$CONSENSUS_NAME".trimmed.fa" \
  -v C_UC="$OUT_DIR/u_temp_995_masked.uc" \
  -v OUT_DIR="$OUT_DIR" \
  -v DEBUG="$DEBUG" \
  '
  # Bin reads based on clustering results
  (FILENAME != R_FA && $1 ~ /S|H/) {
    sub(";size.*", "", $9)
    sub(";size.*", "", $10)
    if ($10 ~ /\*/){ 
      READS[$9]=$9 # Assign read to own cluster
    } else if ($10 !~ /\*/){
      READS[$9]=$10 # Assign read to new cluster
      for (READ in READS){
        if (READS[READ] == $9){
          READS[READ]=$10 # Update reads in same cluster to new cluster
        }
      }
    }
    # Temp output
      print $9, READS[$9], FILENAME > OUT_DIR "/clusters/temp_cluster_assign.txt" 
  }
(FILENAME == C_UC){
    sub(";size.*", "", $9)
    # Remove cluster < 2 reads and assign new cluster names
    if ($1 == "C" && $3+0 > 2){
      CLUSTER_C[$9]="Cluster" $2 ";size=" $3 ";"
      for (READ in READS){
        if (READS[READ] in CLUSTER_C){
          READS_CLUSTERS[READ]=CLUSTER_C[READS[READ]]
        }
      }
    }
  }
  # extract single 
(FILENAME == C_UC) {
  sub(";size.*", "", $9)
  sub(";size.*", "", $10)
  if ($1 == "C" && $3+0 == 1){
    print $9 > "temp_single.lst"
  }
}
# get duo
(FILENAME == C_UC){
    sub(";size.*", "", $9)
    # Remove cluster < 2 reads and assign new cluster names
    if ($1 == "C" && $3+0 == 2){
      DUO_C[$9]="Duo" $2 ";size=" $3 ";"
      for (READ in READS){
        if (READS[READ] in DUO_C){
          READS_DUO[READ]=DUO_C[READS[READ]]
        }
      }
    }
  }
  # Output filtered clusters
  (FILENAME == R_FA && FNR%2==1){
    READ_HEADER=$0
    sub("^>", "", READ_HEADER)
    if (READ_HEADER in CLUSTER_C){
      # Format cluster name
      CLUSTER_NAME=READS_CLUSTERS[READ_HEADER]
      sub(";.*", "", CLUSTER_NAME)
      # Output header to centroid file
      print ">" CLUSTER_C[READ_HEADER] > OUT_DIR "/centroids.fa"
      # Output header to read cluster file
      print ">" READ_HEADER > OUT_DIR "/clusters/" CLUSTER_NAME
      getline
      # Output sequence to centroid file
      print $0 > OUT_DIR "/centroids.fa"
      # Output sequence to read cluster file
      print $0 > OUT_DIR "/clusters/" CLUSTER_NAME
    } else if (READ_HEADER in READS_CLUSTERS){
      # Format cluster name
      CLUSTER_NAME=READS_CLUSTERS[READ_HEADER]
      sub(";.*", "", CLUSTER_NAME)
      # Output header to read cluster file
      print ">" READ_HEADER > OUT_DIR "/clusters/" CLUSTER_NAME
      getline
      # Output sequence to read cluster file
      print $0 > OUT_DIR "/clusters/" CLUSTER_NAME
    }
  }
  # Output filtered duos
  (FILENAME == R_FA && FNR%2==1){
    READ_HEADER=$0
    sub("^>", "", READ_HEADER)
    if (READ_HEADER in DUO_C){
      # Format cluster name
      CLUSTER_NAME=READS_DUO[READ_HEADER]
      sub(";.*", "", CLUSTER_NAME)
      # Output header to read cluster file
      print ">" READ_HEADER > OUT_DIR "/clusters/" CLUSTER_NAME".fa"
      getline
      # Output sequence to read cluster file
      print $0 > OUT_DIR "/clusters/" CLUSTER_NAME".fa"
    } else if (READ_HEADER in READS_DUO){
      # Format cluster name
      CLUSTER_NAME=READS_DUO[READ_HEADER]
      sub(";.*", "", CLUSTER_NAME)
      # Output header to read cluster file
      print ">" READ_HEADER > OUT_DIR "/clusters/" CLUSTER_NAME".fa"
      getline
      # Output sequence to read cluster file
      print $0 > OUT_DIR "/clusters/" CLUSTER_NAME".fa"
    }
    next
  }
  ' \
  $OUT_DIR/u_temp.uc \
  $OUT_DIR/u_temp_995_masked.uc \
  $OUT_DIR/"$CONSENSUS_NAME".trimmed.fa

# Generating consensus from detected clusters
VARIANT_OUT=$OUT_DIR/clustering_consensus
mkdir -p $VARIANT_OUT


cat $OUT_DIR/centroids.fa | $SEQTK seq -l0 - |\
  $GNUPARALLEL \
    --env cluster_consensus \
    --progress \
    -j $CLUSTER_JOBS \
    --recstart ">" \
    -N 1 \
    --pipe \
    "cat | cluster_consensus $VARIANT_OUT $OUT_DIR/clusters \
    $CONSENSUS_FILE $CLUSTER_THREADS"

cat $VARIANT_OUT/*/*consensus.fa > $OUT_DIR/"$CONSENSUS_NAME"_consensus_temp.fa

$USEARCH \
  -fastx_uniques $OUT_DIR/"$CONSENSUS_NAME"_consensus_temp.fa \
  -strand both \
  -fastaout $OUT_DIR/"$CONSENSUS_NAME"_consensus.fa  \
  -uc $OUT_DIR/"$CONSENSUS_NAME"_consensus.uc \
  -sizein \
  -sizeout 

# Extracting non clustered bins (=1) & failed trimming bins

diff <( grep '^>' $IN_DIR/"$CONSENSUS_NAME.fa" | sort ) <( grep '^>' $OUT_DIR/"$CONSENSUS_NAME".trimmed.fa | sort ) | awk 'sub(/^< >/, "")' > $OUT_DIR/"temp_trim.fail.lst"

$SEQTK subseq $OUT_DIR/"$CONSENSUS_NAME".trimmed.fa temp_single.lst > $OUT_DIR/"$CONSENSUS_NAME".single.trimmed.fa
$SEQTK subseq $IN_DIR/"$CONSENSUS_NAME.fa" $OUT_DIR/"temp_trim.fail.lst" > $OUT_DIR/"$CONSENSUS_NAME".FAIL.trimmed.fa

#
for i in $(ls $OUT_DIR/clusters/Duo* | sed 's/.fa//g'); 
do grep ">" "$i".fa  | cut -c2- > $OUT_DIR/clusters/temp_names;

D1=$(cat $OUT_DIR/clusters/temp_names | awk -F ';' '{print $2}' | sed 's/ubs=//g' | head -1)
D2=$(cat $OUT_DIR/clusters/temp_names | awk -F ';' '{print $2}' | sed 's/ubs=//g' | tail -1)
N=$(echo "$i" | awk -F '/' '{print $11}')

if [ $D1 -gt $D2 ];
then 
  head -1 $OUT_DIR/clusters/temp_names | seqtk subseq "$i".fa  - | sed "s/>.*/&_$N/" > "$i".duo.trimmed.fa ;
elif [ $D1 -lt $D2 ]; 
  then
  tail -1 $OUT_DIR/clusters/temp_names | seqtk subseq "$i".fa  - | sed "s/>.*/&_$N/" > "$i".duo.trimmed.fa ;
else
  head -1 $OUT_DIR/clusters/temp_names | seqtk subseq "$i".fa  - | sed "s/>.*/&_$N/" > "$i".duo.trimmed.fa 
fi
done; 

cat $OUT_DIR/clusters/*.duo.trimmed.fa > $OUT_DIR/"$CONSENSUS_NAME".duo.trimmed.fa

if [ "$DEBUG" = "NO" ]; then
  rm $OUT_DIR/*temp* ;
  rm $OUT_DIR/centroids.fa
  rm $OUT_DIR/clusters/*temp*
  rm $OUT_DIR/clusters/*.duo.trimmed.fa
fi

### Testing
exit 0


