#!/bin/bash
# DESCRIPTION
#    Install longread_umi_HIV as conda environment.
#
# IMPLEMENTATION
#    author   Søren Karst (sorenkarst@gmail.com)
#             Ryan Ziels (ziels@mail.ubc.ca)
#             Laurens Lambrechts (laurens.lambrechts@ugent.be)
#    license  GNU General Public License

# Terminal input
BRANCH=${1:-master} # Default to master branch

# Check conda installation ----------------------------------------------------
if [[ -z $(which conda) ]]; then
  # Ask to install
  read -t 1 -n 10000 discard # Clears stdin before read
  read \
    -n 1 \
    -p "Conda not found. Install miniconda3 (y/n)? " \
    ASK_CONDA_INSTALL    
  
  if [ "$ASK_CONDA_INSTALL" == "y" ]; then
    # Install conda
    [ -f Miniconda3-py38_4.8.3-Linux-x86_64.sh ] ||\
      wget "https://repo.anaconda.com/miniconda/Miniconda3-py38_4.8.3-Linux-x86_64.sh"
    bash ./Miniconda3-py38_4.8.3-Linux-x86_64.sh
    echo ""
    echo "#-----------------------------------------------------------"
    echo "Miniconda installed"
    echo ""
    echo "Re-run install_conda.sh to continue longread_umi installation."
    echo ""
    exec bash
  else
    echo ""
	echo "Installation aborted..."
    echo ""
    exit 1 
  fi
else
  echo ""
  echo "Conda found"
  echo "version: $(conda -V)"
  echo ""
fi

# Install longread-UMI conda env ----------------------------------------------
echo ""
echo "Installing longread_umi_HIV conda environment.."
echo ""

# Define conda env yml
echo "name: dev_longread_umi_HIV
channels:
- conda-forge
- bioconda
- defaults
dependencies:
- seqtk=1.3
- python =3.7
- parallel=20191122
- racon=1.4.20
- minimap2=2.17
- pip
- gawk=4.1.3
- cutadapt=2.7
- filtlong=0.2.0
- bwa=0.7.17
- samtools=1.11
- bcftools=1.11
- git
- seqkit=2.4.0
- pip:
  - medaka==1.4.4
" > ./dev_longread_umi_HIV.yml

# Install conda env
conda env create -f ./dev_longread_umi_HIV.yml

eval "$(conda shell.bash hook)"
conda activate dev_longread_umi_HIV || source activate dev_longread_umi_HIV

# Install porechop
$CONDA_PREFIX/bin/pip install \
  git+https://github.com/rrwick/Porechop.git  

# Download longread-UMI from git
git clone \
  --branch "$BRANCH" \
  https://github.com/laulambr/longread_umi_hiv.git \
  $CONDA_PREFIX/longread_umi

# Modify adapters.py
cp \
  $CONDA_PREFIX/longread_umi/scripts/adapters.py \
  $CONDA_PREFIX/lib/python3.7/site-packages/porechop/adapters.py

# Create links to pipeline
find \
  $CONDA_PREFIX/longread_umi/ \
  -name "*.sh" \
  -exec chmod +x {} \;
  
ln -s \
  $CONDA_PREFIX/longread_umi/longread_umi.sh \
  $CONDA_PREFIX/bin/longread_umi
  
  
# Create link to usearch installation
read -t 1 -n 10000 discard
read \
  -p "Type path to usearch excutable and press enter:  " \
  USEARCH_PATH

USEARCH_PATH_F=$(sed -e 's/^"//' -e 's/"$//' <<< "$USEARCH_PATH")
unset USEARCH_PATH

if [[ ! -x "$USEARCH_PATH_F" ]]; then
  echo "File '$USEARCH_PATH_F' is not executable or found."
  read -t 1 -n 10000 discard
  read \
    -n 1 \
    -p "Attempt to make '$USEARCH_PATH_F' excutable (y/n)? " \
    ASK_USEARCH_X    
    echo ""
  if [ "$ASK_USEARCH_X" == "y" ]; then
    chmod +x "$USEARCH_PATH_F"
  else
    echo ""
    echo "Installation aborted ..."
    echo ""
    exit 1 
  fi
fi


ln -s \
  "$USEARCH_PATH_F" \
  $CONDA_PREFIX/bin/usearch  
  
# Check installation
if [[ -z $(which longread_umi) ]]; then
  echo ""
  echo "Can't locate longread_umi_HIV"
  echo "longread_umi_HIV installation failed..."
  echo ""
else
  echo ""
  echo "longread_umi_HIV installation success..."
  echo ""
  echo "Path to conda environment: $CONDA_PREFIX"
  echo "Path to pipeline files: $CONDA_PREFIX/longread_umi"
  echo ""
  echo ""
fi

conda deactivate

# Cleanup
read -t 1 -n 10000 discard
read \
  -n 1 \
  -p "Cleanup install files (y/n)? " \
  CLEAN_INSTALL
  echo ""
  
if [ "$CLEAN_INSTALL" == "y" ]; then
  if [ -f Miniconda3-latest-Linux-x86_64.sh  ]; then 
    rm -f ./Miniconda3-latest-Linux-x86_64.sh
  fi
  if [ -f install_conda.sh  ]; then 
    rm -f ./install_conda.sh
  fi
  if [ -f longread_umi.yml  ]; then 
    rm -f ./longread_umi_HIV.yml
  fi
fi
