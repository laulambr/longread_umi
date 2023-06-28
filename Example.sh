#!/bin/bash

# Install

wget https://raw.githubusercontent.com/laulambr/longread_umi_hiv/master/scripts/install_conda.sh

#wget https://raw.githubusercontent.com/SorenKarst/longread_umi/develop/scripts/install_conda.sh


bash ./install_conda.sh

/data/thesisguest1/Tools/usearch

# medaka fix
pip install protobuf==3.20.*



# STEP1

conda activate longread_umi_HIV || source activate longread_umi_HIV


############# SETTING VARIABLES #############

#primers:
#HIV primer for tagging (inner round Pinzone et al):

hiv_f=AAGTAGTGTGTGCCCGTCTGTTGTGTGAC
hiv_r=GGAAAGTCCCCAGCGGAAAGTCCCTTGTAG

ont_f=CAAGCAGAAGACGGCATACGAGAT
ont_r=AATGATACGGCGACCACCGAGATC

#These should be fixed as standard
#set working directory
wrk=/data/thesisguest1/work/UMI

run_ID=RUN26_bis

prj=$wrk/data_"$run_ID"

sample_id="$run_ID"_guppy_v6.4.2_sup_BC


cd $prj


############# FILTER FILES for HIV #############
#https://www.samformat.info/sam-format-flag

for bc in {07..12}; do 
mkdir $prj/$sample_id"$bc"
#minimap2 -t 30 -ax map-ont $CONDA_PREFIX/longread_umi/scripts/HXB2.fasta $prj/$sample_id"$bc".fastq | samtools fastq -n -f 0x4 - > $prj/$sample_id"$bc"/$sample_id"$bc".nonHIV.fastq 
minimap2 -t 30 -ax map-ont $CONDA_PREFIX/longread_umi/scripts/HXB2.fasta $prj/$sample_id"$bc".fastq | samtools fastq -n -F 0x4 - > $prj/$sample_id"$bc"/$sample_id"$bc".HIV.fastq

done





############# run UMI bin recovery and polishing pipeline #############


#Reads between 100-10000 bp
#medaka model: r103_min_high_g360
#primer pos: 150 (see primer_pos: mean 111)

for bc in {07..12}; do 
longread_umi nanopore_pipeline -d $prj/$sample_id"$bc"/$sample_id"$bc".HIV.fastq -o $prj/$sample_id"$bc"/umi_out_all_length_$sample_id"$bc"/ -v 15 -n 10 -s 200 -e 200 -m 100 -M 10000 -f $ont_f -F $hiv_f  -r $ont_r -R $hiv_r -c 3 -p 1 -q r103_hac_g507 -U 'r103_min_high_g360' -t 30 -T 8 > $prj/$sample_id"$bc"/pipeline_sample_all_fast_nohup.log
done




############# demultiplex UMI bins of each ONT fastq file by PCR barcode #############



mkdir $prj/demultiplex-participants_replicates

### Demultiplex the UMI for each replicate in different participant files

for bc in {07..12}; do 
longread_umi demultiplex -r $prj/$sample_id"$bc"/$sample_id"$bc".HIV.fastq -c $prj/$sample_id"$bc"/umi_out_all_length_$sample_id"$bc"/raconx3_medakax1/consensus_raconx3_medakax1.fa -n '1-8' -o $prj/demultiplex-participants_replicates/$sample_id"$bc" -t 20 -b $CONDA_PREFIX/longread_umi/scripts/barcodes_HIV.tsv -u $prj/$sample_id"$bc"/umi_out_all_length_$sample_id"$bc"/umi_binning/read_binning/umi_bin_map.txt
done

### Mark demultiplexed participant samples from each replicate with a replicate tag



replicate=1

for bc in {07..12}; do  
  for pcr_barcode in {1..8}; do
    awk -v pat="$replicate" '/^>/ {$1=$1";rep0"pat""} 1' $prj/demultiplex-participants_replicates/$sample_id"$bc"/barcode"$pcr_barcode".fa > $prj/demultiplex-participants_replicates/barcode"$pcr_barcode".rep0"$replicate".fa
    awk -v pat="$replicate" '/^>/ {$1=$1";rep0"pat""} 1' $prj/demultiplex-participants_replicates/$sample_id"$bc"/undetermined.fa > $prj/demultiplex-participants_replicates/undetermined.rep0"$replicate".fa
  done
  replicate=$((replicate+1))
  find $prj/demultiplex-participants_replicates/ -type f -empty -delete
done


### Combine UMIs found in different replicate  each participant sample

for pcr_barcode in {1..8}; do
  cat $prj/demultiplex-participants_replicates/barcode"$pcr_barcode".rep*.fa >  $prj/demultiplex-participants_replicates/"$run_ID".barcode"$pcr_barcode".replicates.fa
  cat $prj/demultiplex-participants_replicates/undetermined.rep*.fa >  $prj/demultiplex-participants_replicates/"$run_ID".undetermined.replicates.fa
  find $prj/demultiplex-participants_replicates/ -type f -empty -delete
done 



today=$(date +"%Y_%m_%d")

 rm $prj/demultiplex-participants_replicates/umi_bin_replicate_stats_$today.txt

 echo -e "Run_ID\tBarcode\tID1\tID2\tID3\tID4\tID5\tID6\tID7\tID8\tundertermined" >> $prj/demultiplex-participants_replicates/umi_bin_replicate_stats_$today.txt
 

for bc in {07..12}; do 
    cd $prj/demultiplex-participants_replicates/$sample_id"$bc";
    ID1=$(seqkit stats barcode1.fa -T | awk 'NR!=1 {print $4}');
    ID2=$(seqkit stats barcode2.fa -T | awk 'NR!=1 {print $4}');
    ID3=$(seqkit stats barcode3.fa -T | awk 'NR!=1 {print $4}');
    ID4=$(seqkit stats barcode4.fa -T | awk 'NR!=1 {print $4}');
    ID5=$(seqkit stats barcode5.fa -T | awk 'NR!=1 {print $4}');
    ID6=$(seqkit stats barcode6.fa -T | awk 'NR!=1 {print $4}');
    ID7=$(seqkit stats barcode7.fa -T | awk 'NR!=1 {print $4}');
    ID8=$(seqkit stats barcode8.fa -T | awk 'NR!=1 {print $4}');
    undet=$(seqkit stats undetermined.fa -T | awk 'NR!=1 {print $4}');
    echo -e $run_ID'\t'BC"$bc"'\t'"$ID1"'\t'$ID2'\t'$ID3'\t'$ID4'\t'$ID5'\t'$ID6'\t'$ID7'\t'$ID8'\t'$undet >> $prj/demultiplex-participants_replicates/umi_bin_replicate_stats_$today.txt;  

done

****

### Check if UMI crossover occured, 



mkdir $prj/demultiplex-participants_replicates/check_crossover

replicate=1

for bc in {07..12}; do 
#Names of included umi bins
awk -v pat="$replicate" -F ';' '/^>/ {print substr($1,2,length($1)-5)";rep0"pat""} ' $prj/$sample_id"$bc"/umi_out_all_length_$sample_id"$bc"/raconx3_medakax1/consensus_raconx3_medakax1.fa > $prj/demultiplex-participants_replicates/$sample_id"$bc"/umi_$sample_id"$bc".list
#Obtain & relabel umi bins for each replicate
awk -v pat="$replicate" -F ';' '/^>/ {$0=$1";rep0"pat""} 1' $prj/$sample_id"$bc"/umi_out_all_length_$sample_id"$bc"/umi_binning/umi_ref/umi_ref.fa > $prj/demultiplex-participants_replicates/$sample_id"$bc"/umi_ref_renamed.fa
#extract the used umi bins
seqtk subseq  $prj/demultiplex-participants_replicates/$sample_id"$bc"/umi_ref_renamed.fa $prj/demultiplex-participants_replicates/$sample_id"$bc"/umi_$sample_id"$bc".list >> $prj/demultiplex-participants_replicates/check_crossover/umi.included.fa
#Names of included umi bins linked to real name
awk -v pat="$replicate" -F ';' '/^>/ {print substr($1,2,length($1)-5)";rep0"pat"\t"substr($1,2,length($1)-1)";"$2";rep0"pat} ' $prj/$sample_id"$bc"/umi_out_all_length_$sample_id"$bc"/raconx3_medakax1/consensus_raconx3_medakax1.fa >> $prj/demultiplex-participants_replicates/check_crossover/umi_$sample_id.binsizecol.list
replicate=$((replicate+1))
done

# rename umi ref names to real names
seqkit replace -p '(.+)$' -r '{kv}' -k $prj/demultiplex-participants_replicates/check_crossover/umi_$sample_id.binsizecol.list $prj/demultiplex-participants_replicates/check_crossover/umi.included.fa > $prj/demultiplex-participants_replicates/check_crossover/umi.included.renamed.fa


# Cluster UMI pairs and detect doubles == crosssover
usearch \
  -fastx_uniques \
   $prj/demultiplex-participants_replicates/check_crossover/umi.included.renamed.fa \
  -fastaout $prj/demultiplex-participants_replicates/check_crossover/umi_crossover.fa \
  -minuniquesize 2 \
  -sizeout \
  -uc $prj/demultiplex-participants_replicates/check_crossover/umi_crossover.txt \
  -strand both

#awk '$1=="H" {print $0}' $prj/demultiplex-participants_replicates/check/umi_crossover.txt
#awk '$1=="C" && $3>1 {print $0}' $prj/demultiplex-participants_replicates/check/umi_crossover.txt


### 
# list of non crossover bins
awk '$1=="C" && $3<2 {print $9}' $prj/demultiplex-participants_replicates/check_crossover/umi_crossover.txt >$prj/demultiplex-participants_replicates/check_crossover/bins.noncrossover.ok.list

# export file to R
awk '$1=="H" {print $9 "_" $10}' $prj/demultiplex-participants_replicates/check_crossover/umi_crossover.txt > $prj/demultiplex-participants_replicates/check_crossover/R_input.txt

### filter bins

# Declare an associative array to store the maximum values and corresponding column names
declare -A max_values
declare -A max_columns

# Read the file line by line
while IFS= read -r line; do
  # Split the line into individual entries
  entries=(${line//_/ })

  # Extract the column values
  column1=${entries[0]}
  column2=${entries[1]}

  # Extract the numbers from column 1 and column 2
  number1=$(echo "$column1" | awk -F 'ubs=|;rep' '{print $2}')
  number2=$(echo "$column2" | awk -F 'ubs=|;rep' '{print $2}')

  number1=$(grep -oP 'ubs=\K\d+' <<< "$column1")
  number2=$(grep -oP 'ubs=\K\d+' <<< "$column2")

# Compare the numbers and update the maximum value and corresponding column name in the arrays
if [[ -z ${max_values[$column2]} || $number1 -gt ${max_values[$column2]} ]]; then
  if [[ $number2 -gt $number1 ]]; then
    max_values[$column2]=$number2
    max_columns[$column2]=$column2
  else
    max_values[$column2]=$number1
    max_columns[$column2]=$column1
  fi
fi

done < $prj/demultiplex-participants_replicates/check_crossover/R_input.txt

# Print the maximum values and corresponding column names to the output file
for key in "${!max_values[@]}"; do
  echo ${max_columns[$key]} >> $prj/demultiplex-participants_replicates/check_crossover/bins.crossover.keep.list
done



# list of umi_bins to keep
cat $prj/demultiplex-participants_replicates/check_crossover/bins.crossover.keep.list $prj/demultiplex-participants_replicates/check_crossover/bins.noncrossover.ok.list > $prj/demultiplex-participants_replicates/check_crossover/final.list

# select bins that remain
seqtk subseq $prj/demultiplex-participants_replicates/check_crossover/umi.included.renamed.fa $prj/demultiplex-participants_replicates/check_crossover/final.list >> $prj/demultiplex-participants_replicates/check_crossover/umi.included.fa

# final
for pcr_barcode in {1..8}; do
  seqtk subseq $prj/demultiplex-participants_replicates/"$run_ID".barcode"$pcr_barcode".replicates.fa $prj/demultiplex-participants_replicates/check_crossover/final.list >> $prj/demultiplex-participants_replicates/"$run_ID".barcode"$pcr_barcode".replicates.crossoverremoved.fa
done 
  seqtk subseq $prj/demultiplex-participants_replicates/"$run_ID".undetermined.replicates.fa $prj/demultiplex-participants_replicates/check_crossover/final.list >> $prj/demultiplex-participants_replicates/"$run_ID".undetermined.replicates.crossoverremoved.fa





### Tidy

mv $prj/demultiplex-participants_replicates/check_crossover/umi_crossover.* /$prj/demultiplex-participants_replicates
rm -rf $prj/demultiplex-participants_replicates/barcode*.fa $prj/demultiplex-participants_replicates/check_crossover $prj/demultiplex-participants_replicates/undetermined*.fa $prj/demultiplex-participants_replicates/$sample_id*
find $prj/demultiplex-participants_replicates/ -type f -empty -delete



/data/thesisguest1/work/UMI/data_RUN26_bis/RUN26_bis_guppy_v6.4.2_sup_BC07
conda remove --name longread_umi_HIV --all