#!/bin/bash

# author: Marianne S. Felix
# marianne.sabourin-felix.1@ulaval.ca
# Version : 1.0
# 2018-01-07
#
# ChIP-seq_Alignment-Bowtie2.sh

# Custom display in shell
red=`tput setaf 1`
bold=`tput bold`
normal=`tput sgr0`
ERROR="$red[ERROR]$normal"

# Help
menu="
$bold--------------------------------------------------------
 ChIP-seq_Alignment-Bowtie2.sh version 1.0
--------------------------------------------------------$normal
 This script do the alignment of public data.
 
 This should be used after the trimmingPublicData.sh script.

   Usage   : ./ChIP-seq_Alignment-Bowtie2.sh

$bold [Dependecies] $normal: - Bowtie2
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

# Sigle end (SE) or Paired end (PE) sequencing
read -p "Does the sequencing of your data is single end (SE) or paired end (PE) ? [SE|PE] : " end
until [[ $end == SE ]] || [[ $end == PE ]]
do
    echo " $ERROR : Please answer SE or PE !"
    read -p "Does the sequencing of your dataset is single end (SE) or paired end (PE) ? [SE|PE] : " end
done


# Reference genome basename
# refGenomes/HhrDNA_OSplus49/index/HhrDNA_OSplus49
read -ep "Enter path to the indexed reference genome basename \
(i.e. for sacCer3.1.bt2 etc.: /path/to/index/sacCer3), followed by [ENTER] : " refGen

until [ -f $refGen.1.bt2 ] && [ -f $refGen.2.bt2 ] && \
      [ -f $refGen.3.bt2 ] && [ -f $refGen.4.bt2 ] && \
      [ -f $refGen.rev.1.bt2 ] && [ -f $refGen.rev.2.bt2 ] && [ -n $refGen ]
do
    if [ -z $refGen ] || [ -n "$refGen"* ]
    then
        echo " $ERROR : reference genome basename $refGen not found !"
    elif [ -d $refGen ]
    then
        echo " $ERROR : reference genome basename $refGen is a directory !"
    else
        echo " $ERROR : Some index files are missing for reference genome basename $refGen..."
        echo "Required files are : name.1.bt2, name.2.bt2, name.3.bt2, name.4.bt2, name.rev.1.bt2 and name.rev.2.bt2 "
    fi
    read -ep "Please enter path to the indexed reference genome basename, followed by [ENTER] : " refGen

done

# Number of threads
read -p "Type the number of thread(s) you want to use for Bowtie2, followed by [ENTER] : " threads
until [[ "$threads" =~ ^[1-9]+ ]]
do
    echo " $ERROR : The number of thread(s) must be an integer greater than 0 !"
    read -p "Type the number of thread(s) you want to use for Bowtie2, followed by [ENTER] : " threads
done

echo "--------------------------------------------------------"


###########################
#                         #
#  Alignment, sorting     #
#  and indexing           #
#                         #
###########################

StartTime=$(date +%s)

mkdir $outputFolder


if [[ $end == PE ]]
then
    for forward in $inputFolder/*_1_paired_trim.fastq.gz;
    do
        prefix=`basename $forward _1_paired_trim.fastq.gz`
    
        reverse=$prefix"_2_paired_trim.fastq.gz"
        
        # Alignment
        echo "Running Bowtie2 alignment for file" `basename $forward`"..."
        bowtie2 -p $threads -x $refGen -1 $forward -2 $inputFolder/$reverse \
        | samtools view -F 4 -bS - > $outputFolder/$prefix"_aln.bam"
        
        # Sorting
        echo "Sorting file $prefix""_aln.bam by chromosomal coordinates..."
        samtools sort --threads $threads -o $outputFolder/$prefix"_sorted.bam" \
        -T $outputFolder/$prefix"_sorted" $outputFolder/$prefix"_aln.bam" 
        echo "Alignment file $prefix""_sorted.bam created !"
        
        # Indexing
        echo "Creating index for file $prefix""_sorted.bam..."
        samtools index $outputFolder/$prefix"_sorted.bam"
        echo "Index file $prefix""_sorted.bam.bai created !"
        
        # Remove alignment file
        rm $outputFolder/$prefix"_aln.bam"

    done
fi


if [[ $end == SE ]]
then
    for file in $inputFolder/*.fastq.gz
    do
        prefix=`basename $file .fastq.gz`
        
        # Alignment
        echo "Running Bowtie2 alignment for file" `basename $file`"..."
        bowtie2 -p $threads -x $refGen -U $file | samtools view -F 4 -bS - \
        > $outputFolder/$prefix"_aln.bam"
        echo "Alignment file $prefix""_aln.bam created !"
        
        # Sorting
        echo "Sorting file $prefix""_aln.bam by chromosomal coordinates..."
        samtools sort --threads $threads -o $outputFolder/$prefix"_sorted.bam" \
        -T $outputFolder/$prefix"_sorted" $outputFolder/$prefix"_aln.bam"
        echo "Alignment file $prefix""_sorted.bam created !"
        
        # Indexing
        echo "Creating index for file $prefix""_sorted.bam..."
        samtools index $outputFolder/$prefix"_sorted.bam"
        echo "Index file $prefix""_sorted.bam.bai created !"
        
        # Remove alignment file
        rm $outputFolder/$prefix"_aln.bam"
        
    done
fi


EndTime=$(date +%s)
ElapsedTime=$(($EndTime - $StartTime))

echo "--------------------------------------------------------"
echo "Bowtie2 alignment, sorting and indexing done in $(($ElapsedTime / 3600 )) \
hours $((($ElapsedTime % 3600) / 60)) minutes $(($ElapsedTime % 60)) seconds !"

echo "You can now visually inspect your data or do the peakcalling with MACS2 ! :)"


