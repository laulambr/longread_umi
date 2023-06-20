#
source activate longread_umi_custom

#set working directory
wrk=/data/thesisguest1/work/UMI

#set the following parameters
#what is the name of your project folder where everything will run?

run_ID=RUN07
prj=$wrk/data_"$run_ID"

sample_id="$run_ID"_guppy_v5.0.17_sup_BC

raw_data_loc=$prj/raw_fastq_guppy_v5.0.17_sup


mkdir $prj/demultiplex-participants_replicates

### Demultiplex the UMI for each replicate in different participant files

for bc in {13..28}; do 
longread_umi demultiplex -r $prj/$sample_id"$bc"/$sample_id"$bc".HIV.fastq -c $prj/$sample_id"$bc"/umi_out_all_$sample_id"$bc"/raconx3_medakax1/consensus_raconx3_medakax1.fa -n '1-8' -o $prj/demultiplex-participants_replicates/$sample_id"$bc" -t 20 -b $CONDA_PREFIX/longread_umi/scripts/barcodes_HIV.tsv -u $prj/$sample_id"$bc"/umi_out_all_$sample_id"$bc"/umi_binning/read_binning/umi_bin_map.txt
done

### Mark demultiplexed participant samples from each replicate with a replicate tag


replicate=1

for bc in {13..28}; do 
  for barcode in {1..8}; do
    awk -v pat="$replicate" '/^>/ {$1=$1";rep0"pat""} 1' $prj/demultiplex-participants_replicates/$sample_id"$bc"/barcode"$barcode".fa > $prj/demultiplex-participants_replicates/barcode"$barcode".rep0"$replicate".fa
    awk -v pat="$replicate" '/^>/ {$1=$1";rep0"pat""} 1' $prj/demultiplex-participants_replicates/$sample_id"$bc"/undetermined.fa > $prj/demultiplex-participants_replicates/undetermined.rep0"$replicate".fa
  done
  replicate=$((replicate+1))
done




### Check if UMI crossover occured, 
mkdir $prj/demultiplex-participants_replicates/check

replicate=1

for bc in {13..28}; do 
#Names of included umi bins
awk -v pat="$replicate" -F ';' '/^>/ {print substr($1,2,length($1)-5)";rep0"pat""} ' $prj/$sample_id"$bc"/umi_out_all_$sample_id"$bc"/raconx3_medakax1/consensus_raconx3_medakax1.fa > $prj/demultiplex-participants_replicates/$sample_id"$bc"/umi_$sample_id"$bc".list
#Obtain & relabel umi bins for each replicate
awk -v pat="$replicate" -F ';' '/^>/ {$0=$1";rep0"pat""} 1' $prj/$sample_id"$bc"/umi_out_all_$sample_id"$bc"/umi_binning/umi_ref/umi_ref.fa > $prj/demultiplex-participants_replicates/$sample_id"$bc"/umi_ref_renamed.fa
#extract the used umi bins
seqtk subseq  $prj/demultiplex-participants_replicates/$sample_id"$bc"/umi_ref_renamed.fa $prj/demultiplex-participants_replicates/$sample_id"$bc"/umi_$sample_id"$bc".list >> $prj/demultiplex-participants_replicates/check/umi.included.fa
replicate=$((replicate+1))
done


# Cluster UMI pairs and detect doubles == crosssover
usearch \
  -fastx_uniques \
   $prj/demultiplex-participants_replicates/check/umi.included.fa \
  -fastaout $prj/demultiplex-participants_replicates/check/umi_crossover.fa \
  -minuniquesize 2 \
  -sizeout \
  -uc $prj/demultiplex-participants_replicates/check/umi_crossover.txt \
  -strand both

#awk '$1=="H" {print $0}' $prj/demultiplex-participants_replicates/check/umi_crossover.txt
#awk '$1=="C" && $3>1 {print $0}' $prj/demultiplex-participants_replicates/check/umi_crossover.txt


cd $prj/demultiplex-participants_replicates/



gawk \
'
($3 != "1"&&$1=="C") {
  print $9, $9, $3 >> "umi_crossover_list.txt"
  print $9 >> "cross.ids.list"
}
($1 ~ /S|H/ && $10 !~ /\*/) { print $9, $10, "part" >> "umi_crossover_list.txt"} 
' \
umi_crossover.txt


for barcode in {1..8}; do
while read i; do
umi=$(awk -F ';' -v var="$i" '$0==var {print $1}' cross.ids.list);
rep=$(awk -F ';' -v var="$i" '$0==var {print $2}' cross.ids.list);
if cat "$run_ID".barcode"$barcode".replicates.fa | grep $umi"bin" | grep -q $rep; 
  then echo $i ',' "$run_ID".barcode"$barcode" >> cross_patient.txt
fi
done < ./cross.ids.list
done
rm ./cross.ids.list















### Combine replicate tagged for each participant samples 

for barcode in {1..8}; do
  cat $prj/demultiplex-participants_replicates/barcode"$barcode".rep*.fa >  $prj/demultiplex-participants_replicates/"$run_ID".barcode"$barcode".replicates.fa
  cat $prj/demultiplex-participants_replicates/undetermined.rep*.fa >  $prj/demultiplex-participants_replicates/"$run_ID".undetermined.replicates.fa
done 

mv $prj/demultiplex-participants_replicates/check/umi_crossover.* /$prj/demultiplex-participants_replicates

rm -rf $prj/demultiplex-participants_replicates/barcode*.fa $prj/demultiplex-participants_replicates/check $prj/demultiplex-participants_replicates/undetermined*.fa $prj/demultiplex-participants_replicates/$sample_id*








### in own directory
#for bc in {13..18}; do 
#longread_umi demultiplex -r $prj/$sample_id"$bc"/$sample_id"$bc".HIV.fastq -c $prj/$sample_id"$bc"/umi_out_all_$sample_id"$bc"/raconx3_medakax1/consensus_raconx3_medakax1.fa -n '1-8' -o $prj/$sample_id"$bc"/umi_out_all_$sample_id"$bc"/demultiplex -t 20 -b $CONDA_PREFIX/longread_umi/scripts/barcodes_HIV.tsv -u $prj/$sample_id"$bc"/umi_out_all_$sample_id"$bc"/umi_binning/read_binning/umi_bin_map.txt
#done


