#!/bin/bash

# author: Marianne S. Felix
# marianne.sabourin-felix.1@ulaval.ca
# Version : 3.0
# 2017-01-12
# 2017-06-29 : move mkdir methExtract in the loop in order
#              to analyse many datasets at the same time.
# 2018-04-25 : Allow single-end alignment and analysis
#
# processWGBSdata_V3.sh

# Custom display in shell
red=`tput setaf 1`
bold=`tput bold`
normal=`tput sgr0`
ERROR="$red[ERROR]$normal"

# Help
menu="
$bold--------------------------------------------------------
 processWGBSdata.sh version 3.0
--------------------------------------------------------$normal
 This script do the processing of WGBS data.
 
 This should be used after the trimmingPublicData.sh script.

   Usage   : ./processWGBSdata.sh

$bold [Dependecies] $normal: - Bismark
                 - Bowtie2
                 - Samtools
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
#   Input of files and    #
#   parameters            #
#                         #
###########################

echo "$menu"

# Input directory (must exists and contain *.fastq.gz files)
read -ep "Enter your input folder name, followed by [ENTER] : " inputFolder
until [ -d $inputFolder ] && [[ `find $inputFolder/. -maxdepth 1 -type f -name "*.fastq.gz" 2>/dev/null` ]]
do
    if [ -d $inputFolder ]
    then
        echo " $ERROR : The input folder must contain *.fastq.gz files !"
    else
        echo "$ERROR : Folder $inputFolder not found !"
    fi
    read -ep "Please enter your input folder name, followed by [ENTER] : " inputFolder
done


# Output directory (must not exists)
read -ep "Enter your output folder name, followed by [ENTER] : " outputFolder
while [ -d $outputFolder ]
do
    echo " $ERROR : Directory $outputFolder exists !"
    read -ep "Do you want to replace existing directory ? [yes|no] : " choice
    
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
            echo "$red[ERROR101]$normal : An internal error occured, please report it to marianne.sabourin-felix.1@ulaval.ca"
            exit
            ;;
    esac
done


# Reference genome
read -ep "Enter path to the folder where the Bisulfite_Genome folder is, \
followed by [ENTER] : " refGen

until [ -d $refGen/"Bisulfite_Genome" ] && [ -d $refGen ] && 
([[ `find $refGen/. -maxdepth 1 -type f -name "*.fa" 2>/dev/null` ]] || 
 [[ `find $refGen/. -maxdepth 1 -type l -name "*.fa" 2>/dev/null` ]])
do
    if [ -d $refGen ]
    then
        if [ -d $refGen/"Bisulfite_Genome" ]
        then
            echo "Folder $refGen must contain the genome.fa file !"
        else
            echo "Folder $refGen must contain Bisulfite_Genome folder !"
        fi
        
    else
        echo "Folder $refGen not found !"
    fi
    
    read -ep "Enter path to the folder where the Bisulfite_Genome folder is, \
followed by [ENTER] : " refGen

done


# Path to Bowtie2
read -ep "Enter path to the folder where Bowtie2 is (ex: /usr/local/software/bowtie2-2.2.9), followed by [ENTER] : " bowtie

until [ -d $bowtie ] && [ -f $bowtie/"bowtie2" ]
do
    if [ -d $bowtie ]
    then
        if ! [ -f $bowtie/"bowtie2" ]
        then
            echo "Folder $bowtie must contain the bowtie2 programm !"
        fi
        
    else
        echo "Folder $refGen not found !"
    fi

    read -ep "Enter path to the folder where Bowtie2 is, followed by [ENTER] : " bowtie
done


# Sigle end (SE) or Paired end (PE) sequencing
read -p "Does the sequencing of your data is single end (SE) or paired end (PE) ? [SE|PE] : " end
until [[ $end == SE ]] || [[ $end == PE ]]
do
    echo " $ERROR : Please answer SE or PE !"
    read -p "Does the sequencing of your dataset is single end (SE) or paired end (PE) ? [SE|PE] : " end
done


# Number of threads
read -p "Type the number of thread(s) you want to use for Bowtie2 (x4), followed by [ENTER] : " threads
until [[ "$threads" =~ ^[1-9]+ ]]
do
    echo " $ERROR : The number of thread(s) must be an integer greater than 0 !"
    read -p "Type the number of thread(s) you want to use for Bowtie2 (x4), followed by [ENTER] : " threads
done


#########################
#                       #
####     M A I N     ####
#                       #
#########################

StartTime=$(date +%s)

mkdir $outputFolder

echo "-----------------------"
echo "   BISMARK ALIGNMENT   "
echo "-----------------------"


### Bismark alignment ###
# Paired-end
if [[ $end == PE ]]
then
    for forward in $inputFolder/*_1_paired_trim.fastq.gz;
    do
        prefix=`basename $forward _1_paired_trim.fastq.gz`
    
        reverse=$prefix"_2_paired_trim.fastq.gz"
    
        bismark --bowtie2 --non_directional --path_to_bowtie $bowtie -p 4 $refGen \
        -1 $forward -2 $inputFolder/$reverse -o $outputFolder --basename $prefix

    done
fi


# Single-end
if [[ $end == SE ]]
then
    for file in $inputFolder/*.fastq.gz
    do
        prefix=`basename $file _trim.fastq.gz`
    
        bismark --bowtie2 --non_directional --path_to_bowtie $bowtie -p 4 $refGen \
        $file -o $outputFolder --basename $prefix
        
    done
fi



echo "-----------------------------------"
echo "   BISMARK METHYLATION EXTRACTOR   "
echo "-----------------------------------"

### Bismark methylation extractor ###

# Single-end
if [[ $end == SE ]]
then
    for file in $outputFolder/*.bam
    do

        prefix=`basename $file .bam`
        methExtract=$outputFolder/methExtract_$prefix
        mkdir $methExtract
    
        bismark_methylation_extractor -o $methExtract \
        --ignore 3 --ignore_3prime 3 --multicore $threads \
        --bedGraph --zero_based --cutoff 2 --cytosine_report \
        --genome_folder $refGen $file
    done
fi


# Paired-end
if [[ $end == PE ]]
then
    for file in $outputFolder/*.bam
    do

        prefix=`basename $file _pe.bam`
        methExtract=$outputFolder/methExtract_$prefix
        mkdir $methExtract
    
        bismark_methylation_extractor -p -o $methExtract \
        --ignore 3 --ignore_3prime 3 --ignore_r2 3 --ignore_3prime_r2 3 \
        --multicore $threads --bedGraph --zero_based --cutoff 2 --cytosine_report \
        --genome_folder $refGen $file
    done
fi


echo "-----------------"
echo "   COMPRESSING   "
echo "-----------------"

### Compress unused files ###
mv $methExtract/*.bedGraph $outputFolder
tar -zcvf $methExtract.tar.gz $methExtract
#rm -r $methExtract


echo "--------------------------"
echo "   SORTING AND INDEXING   "
echo "--------------------------"

### Sorting and indexing ###
for file in $outputFolder/*.bam
do
    out=`basename $file .bam`_sorted
    samtools sort --threads $threads -o $outputFolder/$out.bam \
    -T $outputFolder/$out $file
    samtools index $outputFolder/$out.bam
    
    #rm $file
done



EndTime=$(date +%s)
ElapsedTime=$(($EndTime - $StartTime))

echo "--------------------------------------------------------"
echo "Bismark alignment, sorting, indexing and methExtract done in \
$(($ElapsedTime / 3600 )) hours $((($ElapsedTime % 3600) / 60)) minutes \
$(($ElapsedTime % 60)) seconds !"

