#!/bin/bash

# author: Marianne S. Felix
# marianne.sabourin-felix.1@ulaval.ca
# Version : 2.0
# 2016-02-12
# 2016-08-25 : PE option added
#
# trimmingPublicData.sh

# Custom display in shell
red=`tput setaf 1`
bold=`tput bold`
normal=`tput sgr0`
ERROR="$red[ERROR]$normal"

# Help
menu="
$bold--------------------------------------------------------
 trimmingPublicData.sh version 2.0
--------------------------------------------------------$normal
 This script do the trimming of public ChIP-Seq data
 and do the Quality Check.
 
 This should be used after the downloadPublicData.sh script.

   Usage   : ./trimmingPublicData.sh

$bold [Dependecies] $normal: - Trimmomatic-0.33
                 - FastQC
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

# Number of threads
read -p "Type the number of thread(s) you want to use for trimmomatic, followed by [ENTER] : " threads
until [[ "$threads" =~ ^[1-9]+ ]]
do
    echo " $ERROR : The number of thread(s) must be an integer greater than 0 !"
    read -p "Type the number of thread(s) you want to use for trimmomatic, followed by [ENTER] : " threads
done

# Phred score encoding (must be 33 or 64)
read -p "Does your dataset encoding is phred33 (Sanger, Illumina >= 1.8) or phred64 (Illumina < 1.8) ? [33|64] : " phred
until [[ $phred == 33 ]] || [[ $phred == 64 ]]
do
    echo " $ERROR : Please answer 33 or 64 !"
    read -p "Does your dataset encoding is phred33 (Sanger, Illumina >= 1.8) or phred64 (Illumina < 1.8) ? [33|64] : " phred
done

# LEADING (must be between 0 and 41)
read -p "Type the minimum score you want to use for trimmomatic LEADING option, followed by [ENTER] : " leading
while (( $leading < 0 )) || (( $leading > 41 )) || [ -z $leading ] || [[ ${leading:0:1} == - ]]
do
    echo " $ERROR : The minimum score must be an integer between 0 and 41 !"
    read -p "Type the minimum score you want to use for trimmomatic LEADING option, followed by [ENTER] : " leading
done

# TRAILING (must be between 0 and 41)
read -p "Type the minimum score you want to use for trimmomatic TRAILING option, followed by [ENTER] : " trailing
while (( $trailing < 0 )) || (( $trailing > 41 )) || [ -z $trailing ] || [[ ${trailing:0:1} == - ]]
do
    echo " $ERROR : The minimum score must be an integer between 0 and 41 !"
    read -p "Type the minimum score you want to use for trimmomatic TRAILING option, followed by [ENTER] : " trailing
done

# MINLEN (must be greater than 0)
read -p "Type the minimum length you want to use for trimmomatic MINLEN option, followed by [ENTER] : " minlen
until [[ "$minlen" =~ ^[1-9]+ ]]
do
    echo " $ERROR : The minimum length must be an integer between 0 and 100 !"
    read -p "Type the minimum length you want to use for trimmomatic MINLEN option, followed by [ENTER] : " minlen
done

# Adapters
read -p "Do you want to remove adapters from dataset ? [yes|no] : " adaptChoice
until [[ $adaptChoice == yes ]] || [[ $adaptChoice == no ]]
do
    echo " $ERROR : Please answer yes or no !"
    read -p "Do you want to remove adapters from dataset ? [yes|no] : " adaptChoice
done

# Adapter file
# TODO Accept *.fasta format
if [[ $adaptChoice == yes ]]
then
    read -ep "Type the path to the adapter file ex: /home/App/Trimmomatic-0.33/adapters/TruSeq3-SE.fa, followed by [ENTER] : " adapterFile
    until [ -f $adapterFile -o -h $adapterFile ] && [[ ${adapterFile##*.} == fa ]]
    do
        if [ -f $adapterFile -o -h $adapterFile -o -z $adapterFile ]
        then
            echo " $ERROR : Adapter file must be in *.fa format !"
        else
            echo " $ERROR : Adapter file $adapterFile not found !"
        fi
    
        read -ep "Type the path to the adapter file ex: /home/App/Trimmomatic-0.33/adapters/TruSeq3-SE.fa, followed by [ENTER] : " adapterFile
    done
    
    # Set the ILLUMINACLIP default settings (if the user say no to clipChoice)
    seedMismatch=2
    palindromeClipThreshold=30
    simpleClipThreshold=10
    
    read -p "Do you want to change the ILLUMINACLIP settings ? (Default : ILLUMINACLIP:adapterFile:2:30:10) [yes|no] : " clipChoice
    until [[ $clipChoice == yes ]] || [[ $clipChoice == no ]]
    do
        echo " $ERROR : Please answer yes or no !"
        read -p read -p "Do you want to change the ILLUMINACLIP settings ? (Default : ILLUMINACLIP:adaptFile:2:30:10) [yes|no] : " clipChoice
    done
    
    # TODO Change the range of possible answer...
    if [[ $clipChoice == yes ]]
    then
        read -p "Type the maximal seed mismatch you want to allow for trimmomatic ILLUMINACLIP option, followed by [ENTER] : " seedMismatch
        until [[ "$seedMismatch" =~ ^[0-9]+ ]]
        do
            echo " $ERROR : The maximal seed mismatch must be an integer equal or greater than 0 !"
            read -p "Type the maximal seed mismatch you want to allow for trimmomatic ILLUMINACLIP option (Default=2), followed by [ENTER] : " seedMismatch
        done
        
        read -p "Type the minimum score for palindrome match you want to allow for trimmomatic ILLUMINACLIP option (Defaut=30), followed by [ENTER] : " palindromeClipThreshold
        until [[ "$palindromeClipThreshold" =~ ^[0-9]+ ]]
        do
            echo " $ERROR : The maximal seed mismatch must be an integer equal or greater than 0 !"
            read -p "Type the minimum score for palindrome match you want to allow for trimmomatic ILLUMINACLIP option (Defaut=30), followed by [ENTER] : " palindromeClipThreshold
        done
        
        read -p "Type the minimum score for simple match you want to allow for trimmomatic ILLUMINACLIP option (Defaut=10), followed by [ENTER] : " simpleClipThreshold
        until [[ "$simpleClipThreshold" =~ ^[0-9]+ ]]
        do
            echo " $ERROR : The maximal seed mismatch must be an integer equal or greater than 0 !"
            read -p "Type the minimum score for simple match you want to allow for trimmomatic ILLUMINACLIP option (Defaut=10), followed by [ENTER] : " simpleClipThreshold
        done
    fi
        
fi

echo "--------------------------------------------------------"


###########################
#                         #
#  Trimming and analysis  #
#                         #
###########################

StartTime=$(date +%s)

mkdir $outputFolder

# Find Trimmomatic path
pathTrimmomatic=`locate Trimmomatic-0.33/trimmomatic-0.33.jar`

#/home/bickj/software/Trimmomatic-0.33/adapters/TruSeq2-PE.fa

# TODO Handle PE datasets
for file in $inputFolder/*.fastq.gz
do
    prefix=`basename $file .fastq.gz`
    
    # Skip file_2 (only process file_1)
    if [[ $prefix == *_2 ]]
    then
        continue
    
    elif [[ $prefix == *_1 ]]
    then
        
        echo "Running Trimmomatic trimming for file" `basename $file`"..."
        
        # If only file_1.fastq.gz if found it's probably because it hase been downoloaded as PE when it was SE
        baseFile=`basename $prefix _1`
        if [[ ! `find $inputFolder/. -maxdepth 1 -type f -name $baseFile"_2.fastq.gz" 2>/dev/null` ]]
        then
            end=SE
        
        else
            # Parameters assignation
            forward=$file
            reverse=$inputFolder/$baseFile"_2.fastq.gz"
            
            outForPaired=$outputFolder/$baseFile"_1_paired_trim.fastq.gz"
            outForUnpaired=$outputFolder/$baseFile"_1_unpaired_trim.fastq.gz"
            
            outRevPaired=$outputFolder/$baseFile"_2_paired_trim.fastq.gz"
            outRevUnpaired=$outputFolder/$baseFile"_2_unpaired_trim.fastq.gz"
            
            #echo "forward : $forward"
            #echo "reverse : $reverse"
            #echo "outForPaired : $outForPaired"
            #echo "outForUnpaired : $outForUnpaired"
            #echo "outRevPaired : $outRevPaired"
            #echo "outRevUnpaired : $outRevUnpaired"
            
        fi
    fi
    
    
    if [[ $end == SE ]] && [[ $adaptChoice == no ]]
    then
        java -jar $pathTrimmomatic $end -threads $threads -phred$phred $file \
        $outputFolder/$prefix"_trim.fastq.gz" LEADING:$leading TRAILING:$trailing \
        MINLEN:$minlen
    
    elif [[ $end == SE ]] && [[ $adaptChoice == yes ]]
    then
        java -jar $pathTrimmomatic $end -threads $threads -phred$phred $file \
        $outputFolder/$prefix"_trim.fastq.gz" LEADING:$leading TRAILING:$trailing \
        ILLUMINACLIP:$adapterFile:$seedMismatch:$palindromeClipThreshold:$simpleClipThreshold \
        MINLEN:$minlen
    
    elif [[ $end == PE ]] && [[ $adaptChoice == no ]]
    then
        java -jar $pathTrimmomatic $end -threads $threads -phred$phred $forward \
        $reverse $outForPaired $outForUnpaired $outRevPaired $outRevUnpaired \
        LEADING:$leading TRAILING:$trailing MINLEN:$minlen
    
    elif [[ $end == PE ]] && [[ $adaptChoice == yes ]]
    then
        java -jar $pathTrimmomatic $end -threads $threads -phred$phred $forward \
        $reverse $outForPaired $outForUnpaired $outRevPaired $outRevUnpaired \
        ILLUMINACLIP:$adapterFile:$seedMismatch:$palindromeClipThreshold:$simpleClipThreshold \
        LEADING:$leading TRAILING:$trailing MINLEN:$minlen
    
    else
        echo "$red[ERROR102]$normal : An internal error occured, please report it to marianne.sabourin-felix.1@ulaval.ca"
        exit
    fi
    
    
    # Quality Check
    if [[ $end == SE ]]
    then
        # Single-end
        echo "Quality check for file $file.fastq.gz..."
        fastqc -q $outputFolder/$prefix"_trim.fastq.gz"
        echo "Quality check done for file $file.fastq.gz !"
        
        # Remove redundant files
        rm $outputFolder/$prefix"_trim.fastc.zip"
    
    else
        
        for file in $outForPaired $outForUnpaired $outRevPaired $outRevUnpaired
        do
            echo "Quality check for file $file..."
            fastqc -q $file
            echo "Quality check done for file $file !"
            # Remove redundant files
            fastqcZip=`basename $file .fastq.gz`_fastqc.zip
            rm $outputFolder/$fastqcZip
        done
    fi  

done


EndTime=$(date +%s)
ElapsedTime=$(($EndTime - $StartTime))

echo "--------------------------------------------------------"
echo "Trimmomatic trimming and FastQC analysis done in $(($ElapsedTime / 3600 )) \
hours $((($ElapsedTime % 3600) / 60)) minutes $(($ElapsedTime % 60)) seconds !"

echo "If the quality check analysis is not right, you can re-run this script with other parameters ! :)"
echo "--------------------------------------------------------"
echo " $bold[TIP]$normal : You will may want to change the minimum score for LEADING or TRAILING option."
if [[ $adaptChoice == yes ]]
then
    echo " $bold[TIP]$normal : You will may want to change the Adapter file."
    echo " $bold[TIP]$normal : You will may want to change the ILLUMINACLIP settings."
    echo "         - min score for palindrome match should be in the range of 30 (up to 50)"
    echo "         - min score for simple match should be between 7 and 15"
fi


