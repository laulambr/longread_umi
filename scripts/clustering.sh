#!/bin/bash
# DESCRIPTION
#    Phasing and variant calling of consensus sequences to obtain
#    lowest possible error rate. Part of the longread-UMI-pipeline.
#    
# IMPLEMENTATION
#    author   Søren Karst (sorenkarst@gmail.com)
#             Ryan Ziels (ziels@mail.ubc.ca)
#    license  GNU General Public License
# TODO


### Description ----------------------------------------------------------------

USAGE="
-- longread_umi clustering: Cluster related UMI consensus sequences.
   Clustering for UMI consensus sequences and creates a consensus sequence with >=3x coverage. 
   Reads are initially grouped by read clustering at 99.5% identity and a centroid sequence is picked.
   The centroid sequence is used as a mapping reference for all reads in the cluster.
   
usage: $(basename "$0" .sh) [-h -b] (-c file -F string -R string -m value -M value -o dir -i dir -t value ) 

where:
    -h  Show this help text.
    -c  UMI consensus file.
    -F  Forward tagging primer sequence.
    -R  Reverse tagging primer sequence.
    -m  Minimum read length.
    -M  Maximum read length.
    -o  Output directory.
    -i  Input directory.
    -t  Number of threads to use. [Default = 1]
    -b  Debug flag. Keep temp files. [Default = NO]
"

### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hzc:F:R:m:M:i:o:t:b' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    c) CONSENSUS_FILE=$OPTARG;;
    F) FW2=$OPTARG;;
    R) RV2=$OPTARG;;
    m) MIN_LENGTH=$OPTARG;;
    M) MAX_LENGTH=$OPTARG;;  
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
bam_read_split() {
  ### Description
  # Input is a bam file containing reads mapped to a reference.
  # The function splits the bam file per read. The output is per
  # read mapping information in sam format.
 
  # Input 
  IN_BAM=$1

  # Split
  $SAMTOOLS view -h $IN_BAM |\
    $GAWK 'NR<=3{h[++l]=$0; next}
      NR>3{s[++n]=$0; next}
      END {
        for (i in s){
          print h[1]"\n"h[2]"\n"h[3]"\n"s[i]
      }
    }'
}
export -f bam_read_split

extract_vars(){
  ### Description
  # Input is a sam file of a single read mapped to a reference.
  # The function extracts base information at diagnostic
  # coordinates based on a reference positions. The output is a 
  # read name and diagnostic bases, which are used as a bin id.
   
  # Input
  local IN=$(cat)
  local REG=$1
  local REF=$2

  # Format name
  local UMI_NAME=$(echo "$IN" |\
    $GAWK 'NR==4{print $1}')

  # mpileup
  if [ -z "$REG" ]; then
    echo "${UMI_NAME} all"
  elif [ ! -z "$REG" ]; then
    echo "$IN" |\
      $SAMTOOLS view -b - |\
      $BCFTOOLS mpileup \
        -Ov \
        -a "FORMAT/AD,FORMAT/DP" \
        -f $REF \
        -t $REG \
        - |\
      $GAWK -v un="$UMI_NAME" '!/^#/{
        gsub(",|<\\*>", "", $5)
        df=length($4)-length($5)
        if ($5 == ""){p = p $2 $4}
        else if (df > 0){
          for(c=1;c<=df;c++){d=d"-"}
          p = p $2 $5 d
        }
        else if (df <= 0){p = p $2 $5}
      } END {
        print un, p
      }'
  fi
  # Testing
  # samtools view -h Cluster25.bam | head -n 4 | extract_vars "$REG" Cluster25_con.fa > test.txt
}
export -f extract_vars

var_bin() {
  ### Description
  # Input is file of reads in fasta format and a list of read names
  # and diagnostic bases (bin id). The extracts reads and divides them
  # into thier assigned bins. 

  # Input
  local IN=$1
  local READS=$2
  local CN=$3

  # Binning
  cat $IN |\
  $GAWK -v out="$CN" 'NR==FNR{
    seq[$1] = $2
    if (!($2 in cl)){vn[$2]=++j}
    cl[$2]++
    next
  }
  FNR%2==1{
    # Read consensus read name
    sn = substr($0, 2);
    # Check read has been binned
    if (sn in seq){
      # Write read to varians bin if bin is big enough
      if (cl[seq[sn]] >= 3){
        print ">" sn > out"_var"vn[seq[sn]]"_bin.fa"
        getline
        print > out"_var"vn[seq[sn]]"_bin.fa"
        print "var"vn[seq[sn]], seq[sn] > out"variant_names.txt"
      }
      # Write read to none bin if bin is too small
      else if (cl[seq[i]] < 3){
        print ">" sn > out"_none.fa"
        getline
        print > out"_none.fa"
      }
    sn = ""
    }    
  }' - <($SEQTK seq -l0 $READS)
  
  # Testing
  #samtools view -h Cluster0.bam | head -n 4 |\
  # extract_vars "$REG" Cluster0_con.fa | var_bin \
  #../../../final/consensus_sracon_medaka_medaka_b30.fa test
  # seqtk seq -l0 ../../clusters/Cluster0 > test.fa
}
export -f var_bin

phased_consensus(){
  ### Description
  # Input is bin file of reads in fasta format and an output
  # name. The binned reads are used to create a raw consensus
  # using usearch, which is polished using bcftools.

  # Input
  local PB=$1
  local PB_NAME=$2

  # Format
  local OUT=${PB%/*}

  # Usearch raw consensus
  $USEARCH \
    -cluster_fast $PB \
    -id 0.99 \
    -strand both \
    -consout $OUT/${PB_NAME}_bincon.fa \
    -sizeout \
    -relabel ${PB_NAME}

  # Polish with $BCFTOOLS
  $MINIMAP2 \
    -t 1 \
    -ax map-ont \
    $OUT/${PB_NAME}_bincon.fa \
    $PB |\
    $SAMTOOLS view -b - \
    > $OUT/${PB_NAME}_temp.bam

  $SAMTOOLS sort \
    -o $OUT/${PB_NAME}.bam \
    $OUT/${PB_NAME}_temp.bam
  $SAMTOOLS index $OUT/${PB_NAME}.bam
  rm $OUT/${PB_NAME}_temp.bam   

  $BCFTOOLS mpileup \
    -Ov \
    -d 1000000 \
    -L 1000000 \
    -a "FORMAT/AD,FORMAT/DP" \
    -f $OUT/${PB_NAME}_bincon.fa \
    $OUT/${PB_NAME}.bam \
    > $OUT/${PB_NAME}_pileup.vcf

  $BCFTOOLS norm \
    -f $OUT/${PB_NAME}_bincon.fa \
    $OUT/${PB_NAME}_pileup.vcf \
    -Ov \
    -o $OUT/${PB_NAME}_norm.vcf

  $BCFTOOLS view -i 'AD[0:1]/FORMAT/DP>0.5' \
    -Oz \
    $OUT/${PB_NAME}_norm.vcf \
    > $OUT/${PB_NAME}_polish.vcf.gz
  $BCFTOOLS index $OUT/${PB_NAME}_polish.vcf.gz
    
  $BCFTOOLS consensus $OUT/${PB_NAME}_polish.vcf.gz \
    -f $OUT/${PB_NAME}_bincon.fa |\
    $SEQTK seq -l0 - |\
      $GAWK -v name="$PB_NAME" -v out="$OUT" '
        /^>/{
          ++n;
          if (n == 1){
            split($0,size,";")
            gsub("size=", ";depth=", size[2])
            print ">" name size[2] > out "/" name "_variant.fa"
            getline
            print > out "/" name "_variant.fa"
          }
        } END {
          if (n > 1){msg = "warning"} else {msg="ok"}
          print name, n, msg >> out "/bincon_log.txt"
        }'
}
export -f phased_consensus

phasing () {
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
    -ax asm5\
    $CLUSTER_OUT/${CLUSTER_NAME}_con.fa \
    $CLUSTER_READ_DIR/${CLUSTER_NAME} |\
    $SAMTOOLS view -b - \
    > $CLUSTER_OUT/temp.bam

  $SAMTOOLS sort \
    -o $CLUSTER_OUT/${CLUSTER_NAME}.bam \
    $CLUSTER_OUT/temp.bam
  $SAMTOOLS index $CLUSTER_OUT/${CLUSTER_NAME}.bam
  rm $CLUSTER_OUT/temp.bam   

  # Call variants
  $BCFTOOLS mpileup \
    -Ov \
    -d 1000000 \
    -L 1000000 \
    -a "FORMAT/AD,FORMAT/DP" \
    -f $CLUSTER_OUT/${CLUSTER_NAME}_con.fa \
    $CLUSTER_OUT/${CLUSTER_NAME}.bam \
    > $CLUSTER_OUT/variants.vcf
  $BCFTOOLS view -i 'AD[0:1-]>2' \
    $CLUSTER_OUT/variants.vcf \
    > $CLUSTER_OUT/calls.vcf

  # Convert variants to list
  REG=$(
    $GAWK '
      # Skip comment lines
      !/^#/ {
        reg[$1":"$2]++
       } END {
         for (i in reg){
           if(++n == 1) regout=i
           else if(++n > 1) regout=regout","i
         }
         print regout
       }' $CLUSTER_OUT/calls.vcf

      # Homopolymers effects should have been removed.
      # If not test: && /tolower($5) !~ /a{3,}|t{3,}|c{3,}|g{3,}/
      # but remember this code does not take SNPs close to hps in to account
      # Exampel: GAAA -> AAAA
  )

  # Phasing

  ## Export extract_vars and dependencies
  export -f extract_vars
  export SAMTOOLS=$SAMTOOLS
  export BCFTOOLS=$BCFTOOLS
  export GAWK=$GAWK
  export SEQTK=$SEQTK 

  bam_read_split $CLUSTER_OUT/${CLUSTER_NAME}.bam |\
    $GNUPARALLEL \
      --env extract_vars \
      -j $CLUSTER_THREADS \
      -L 4 \
	  -N 1 \
	  --pipe \
      "cat | extract_vars \"$REG\" \"$CLUSTER_OUT/${CLUSTER_NAME}_con.fa\"" \
      > $CLUSTER_OUT/read_variants.txt

  var_bin \
    $CLUSTER_OUT/read_variants.txt \
    $CONSENSUS_FILE \
    $CLUSTER_OUT/${CLUSTER_NAME}

  # Consensus from phased reads

  ## Export phased_consensus and dependencies
  export -f phased_consensus
  export SAMTOOLS=$SAMTOOLS
  export BCFTOOLS=$BCFTOOLS
  export GAWK=$GAWK
  export MINIMAP2=$MINIMAP2
  export USEARCH=$USEARCH

  find $CLUSTER_OUT -type f -name "*_bin.fa" |\
    $GNUPARALLEL \
	  --env phased_consensus \
	  -j $CLUSTER_THREADS \
      --rpl '{name} s:.*/::; s/_bin.fa$//' \
      "phased_consensus {} {name}"
}
export -f phasing

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

### Consensus phasing pipeline -------------------------------------------------

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
    -id 0.995 \
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

if [ "$DEBUG" = "NO" ]; then
  rm $OUT_DIR/*temp* ;
  rm $OUT_DIR/centroids.fa
fi

### Testing
exit 0

