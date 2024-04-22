# longread_umi_hiv

Tool set for analyzing HIV-PULSE data. This repository is a modified branch of the [`longread_umi`](https://github.com/SorenKarst/longread_umi) pipeline developed by [Karst et al., 2021, Nature Methods](https://doi.org/10.1038/s41592-020-01041-y). 

**Table of contents**
- [Installation](#installation)
- [Usage](#usage)

**Citation**  
Lambrechts et al. (2023). HIV-PULSE: A long-read sequencing assay for high-throughput near full-length HIV-1 proviral genome characterization [NAR](https://academic.oup.com/nar/article/51/20/e102/7306674).

## Installation

### 1. Conda

1. Download and install [usearch](https://drive5.com/usearch/download.html) >=10 according to instructions on the website
   
   
2. Download installer script from terminal  
   ```
   wget https://raw.githubusercontent.com/laulambr/longread_umi_hiv/master/scripts/install_conda.sh
   ```
   
3. Run installation script from terminal and follow instructions 
   ```
   bash ./install_conda.sh
   
4. If miniconda was installed along with the pipeline, initiate conda and refresh terminal before using pipeline.  
   ```
   conda init; source ~/.bashrc
   
5. Activate and deactivate conda environment
   
   ```
   conda activate longread_umi_HIV
   ...
   conda deactivate
 
6. Test if installation was succesfull by running following command when environment is activated
   ```
   longread_umi -h

### 2. Install medaka
1. Create directory for medaka installation
    ```
    mkdir /path/to/medaka
    ```

2. Install medaka
    ```
    cd /path/to/medaka
    python3 -m venv medaka --prompt "medaka"
    source medaka/bin/activate
    pip install --upgrade pip
    pip install medaka
    medaka tools download_models
    deactivate
    ```   
2. Edit medaka location in scripts/dependencies.sh
   ```
   export MEDAKA_ENV_START="source /path/to/medaka/bin/activate"
   ```

## Usage

### Main tools


`longread_umi nanopore`: Generate UMI bin sequences from Nanopore data.

`longread_umi demultiplex`: Dual barcode demultiplexing.

`longread_umi cluster`: Cluster related UMI bin sequences. 

```   
usage: 

longread_umi nanopore_pipeline [-h] [ -k flag] (-d file -v value -o dir -s value) 
(-e value -m value -M value -f string -F string -r string -R string )
( -c value -p value -n value -u dir -U string -t value -T value ) 

where:
    -h  Show this help text.
    -d  Single file containing raw Nanopore data in fastq format.
    -v  Minimum read coverage for using UMI consensus sequences for 
        variant calling.
    -o  Output directory.
    -s  Check start of read up to s bp for UMIs.
    -e  Check end of read up to f bp for UMIs.
    -m  Minimum read length.
    -M  Maximum read length.
    -f  Forward adaptor sequence. 
    -F  Forward primer sequence.
    -r  Reverse adaptor sequence.
    -R  Reverse primer sequence.
    -c  Number of iterative rounds of consensus calling with Racon.
    -p  Number of iterative rounds of consensus calling with Medaka.
    -q  Medaka model used for polishing. r941_min_high_g360, r103_min_high_g360 etc.
    -k  Flag for keeping failed bins in output.
    -n  Process n number of bins. If not defined all bins are processed.
        Pratical for testing large datasets.
    -u  Directory with UMI binned reads.
    -U  UMI filter settings. Define settings for:
        - UMI match error mean (UMEM): Mean match error between reads in a bin
          and the UMI reference.
        - UMI match error SD (UMESD): Standard deviation for match error between
          reads in a bin and the UMI reference.
        - Bin cluster ratio (BCR): Ratio between UMI bin size and UMI cluster size.
        - Read orientation ratio (ROR): n(+ strand reads)/n(all reads). '0' is the
          means disabled.
        Settings can be provided as a string: 'UMEM/UMESD/BCR/ROR'
        Or as a preset:
        - 'r941_min_high_g360' == '3;2;6;0.3'
        - 'r103_min_high_g360' == '3;2.5;12;0.3'
    -t  Number of threads to use.
    -T  Number of medaka jobs to start. Threads pr. job is threads/jobs.
        [Default = 1].		
```

```   
usage: 

longread_umi demultiplex [-h] (-c file -r file -u file -o dir -b file)
(-p string -n range -t value) 

where:
    -h  Show this help text.
    -c  UMI consensus sequences that need demultiplexing.
    -r  Raw read sequences that were used to generate
        the consensus sequences.
    -u  List of raw read names and UMI bin assignments.
    -o  Output directory.
    -b  File containing barcodes. 
        [Default = "$BARCODES"].
    -p  Barcode name prefix [Default = 'barcode'].
    -n  Barcode numbers used. [Default  = '1-120'].
    -t  Number of threads used.
```

```   
usage: 

longread_umi clustering [-h -b] (-c file -F string -R string -m value -M value -p value -o dir -i dir -t value ) 

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
```
### Example script

See [following file](https://github.com/laulambr/longread_umi_hiv/blob/master/Example.sh) for example commands to perform HIV-PULSE data analysis.
  
## License
[GNU General Public License, version 3](LICENSE)