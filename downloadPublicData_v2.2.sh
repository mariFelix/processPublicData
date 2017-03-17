#!/bin/bash

# author: Marianne S. Felix
# marianne.sabourin-felix.1@ulaval.ca
# Version : 2.1
# 2016-02-12
# 2016-02-23
#
# downloadPublicData.sh

# Custom display in shell
red=`tput setaf 1`
bold=`tput bold`
normal=`tput sgr0`
ERROR="$red[ERROR]$normal"

# Help
menu="
$bold--------------------------------------------------------
 downloadPublicData.sh version 2.1
--------------------------------------------------------$normal
 This script dowload public data (SRA format) 
 and do the Quality Check.

   Usage   : ./downloadPublicData.sh

$bold [Dependecies] $normal: - SRA Toolkit
                 - FastQC

$red CAUTION $normal: This script will only process files with a valid accession name.

           If some file(s) can't be found on NCBI SRA database, you will
           need to restart this script with the correct accession name
           of the datasets that didn't work the first time.
           
           This script assumes that every dataset dowloaded at the same
           time are all single-end or all paired-end.
$bold--------------------------------------------------------$normal
"

# Display help
if [[ $1 == -h ]] || [[ $1 == --help ]]
then
    echo "$menu"
    exit
fi


###########################
#                         #
#   Input of files        #
#                         #
###########################

echo "$menu"

# User's input
read -p "Enter your output folder name, followed by [ENTER] : " outputFolder

# Validation of output directory (must not exists)
while [ -d $outputFolder ]
do
    echo " $ERROR : Directory $outputFolder exists !"
    read -p "Do you want to replace existing directory ? [yes|no] : " choice
    
    until [[ $choice == yes ]] || [[ $choice == no ]]
    do
        echo " $ERROR : Please answer yes or no !"
        read -p "Do you want to replace existing directory ? [yes|no] : " choice
    done
    
    case $choice in

        yes)
            echo "Directory $outputFolder will be overwritten !"
            rm -r $outputFolder
            ;;
            
        no)
            echo "Please, choose an other output name."
            echo "Exiting program..."
            exit
            ;;
        *)
            echo " $red[ERROR101]$normal : An internal error occured, please report it to marianne.sabourin-felix.1@ulaval.ca"
            exit
            ;;
    esac
done

# Sigle end (SE) or Paired end (PE) sequencing
read -p "Does the sequencing of your data is single end (SE) or paired end (PE) ? [SE|PE] : " end
until [[ $end == SE ]] || [[ $end == PE ]]
do
    echo " $ERROR : Please answer SE or PE !"
    read -p "Does the sequencing of your dataset is single end (SE) or paired end (PE) ? [SE|PE] : " end
done

### Input of file(s) ###
# TODO Accept fastq.gz files (from ENCODE) 
read -p "Type the number of SRA dataset(s) you want to download, followed by [ENTER] : " numFiles

# Validation of user's input
until [[ "$numFiles" =~ ^[1-9]+ ]]
do
    echo " $ERROR : The number of SRA dataset(s) must be an integer greater than zero!"
    read -p "Type the number of SRA dataset(s) you want to download, followed by [ENTER] : " numFiles
done

# User's file input
for (( fileNo=1; fileNo<=$numFiles; fileNo++ ))
do
    read -p "Enter your dataset name #$fileNo, followed by [ENTER] : " inputFile
    files[fileNo]=$inputFile
done

echo "--------------------------------------------------------"


###########################
#                         #
#  Download and analysis  #
#                         #
###########################

mkdir $outputFolder
cd $outputFolder

# Index the number of file not found
error=0

for file in ${files[@]}
do
    ### Dowload from SRA database ###
    echo "Dowloading $file to fastq.gz format..."
    if [[ $end == SE ]]
    then
        # Single-end
        fastq-dump --gzip $file 2> /dev/null
    else
        # Paired-end
        fastq-dump --split-files --gzip $file 2> /dev/null
    fi
    
    # Do not continue if the file can't be found
    if [[ ! -f $file.fastq.gz ]] && [[ ! -f $file"_1.fastq.gz" ]] && [[ ! -f $file"_2.fastq.gz" ]]
    then
        echo " $ERROR : Dataset $file not found on the NCBI SRA database !"
        error=$((error+1))
        continue
    fi

    echo "Dataset $file downloaded !"
    
    ### Perform Quality Check ###
    if [[ $end == SE ]]
    then
        # Single-end
        echo "Quality check for file $file.fastq.gz..."
        fastqc -q $file.fastq.gz
        echo "Quality check done for file $file.fastq.gz !"
        
        # Remove redundant files
        rm $file"_fastqc.zip"
    else
        # Paired-end (forward)
        echo "Quality check for files" $file"_1.fastq.gz..."
        fastqc -q $file"_1.fastq.gz"
        echo "Quality check done for file" $file"_1.fastq.gz !"
        # Remove redundant files
        rm $file"_1_fastqc.zip"
        
        # Skip file_2 (if PE was choose but the dataset was single-end)
        [[ ! -f $file"_2.fastq.gz" ]] && continue
        
        # Paired-end (reverse)
        echo "Quality check for files" $file"_2.fastq.gz..."
        fastqc -q $file"_2.fastq.gz"
        echo "Quality check done for file" $file"_2.fastq.gz !"
        # Remove redundant files
        rm $file"_2_fastqc.zip"
    fi    
done

#if [[ $end == PE ]] && [[ ! `find . -maxdepth 1 -type f -name "*_2.fastq.gz" 2>/dev/null` ]]
#then
    
    
#fi

#if [[ `ls $outputFolder/*.fastq.gz 2>&1 /dev/null` ]]
if [[ $error < $numFiles ]]
then
    if [[ $error != 0 ]]
    then
        echo "--------------------------------------------------------"
        echo " $ERROR : A problem occured while dowloading the files."
        echo "           Please check the $error input names that didn't download and try again !"
    fi
    
    echo "--------------------------------------------------------"
    echo "Conversion to *.fastq.gz and FastQC analysis done !"
    echo "After the quality check analysis you will be able to choose your parameters for the trimming step ! :)"
else
    echo "--------------------------------------------------------"
    echo " $ERROR : A problem occured while dowloading the files."
    echo "           Please check all your input names and try again !"  
fi

cd ../
