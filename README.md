# processPublicData

2017-03-17

Marianne S. Felix

marianne.sabourin-felix@hotmail.com
-----------------------------------

## INTRODUCTION

Processing public datasets.

## WARNINGS

These scripts are provided "as is", without warranty of any kind. You can redistribute and/or modify it under the terms of the GNU General Public License v3.

## SOFTWARE DEPENDENCIES

downloadPublicData_v.2.2 :

* SRA Toolkit
* FastQC

trimmingPublicData_v2 :

* Trimmomatic-0.33
* FastQC

alignPublicData :

* Bowtie2
* Samtools

processWGBSdata :

* Bismark
* Bowtie2
* Samtools


## HOW TO USE IT

###### downloadPublicData_v2.2

To run the scripts do :

```
./downloadPublicData_v2.2.sh 
```

Then follow the instructions :

```
Enter your output folder name, followed by [ENTER] :
```

```
Does the sequencing of your data is single end (SE) or paired end (PE) ? [SE|PE] :
```

```
Type the number of SRA dataset(s) you want to download, followed by [ENTER] :
```

```
Enter your dataset name #X, followed by [ENTER] :
...
```
When the file(s) is/are downloaded the following message will appear :

```
Dowloading SRRXXXXXXX to fastq.gz format...
Dataset SRRXXXXXXX downloaded !"
...
```
    
```
Quality check for file SRRXXXXXXX.fastq.gz...
Quality check done for file SRRXXXXXXX.fastq.gz
...
```

###### trimmingPublicData_v2

To run the scripts do :

```
./trimmingPublicData_v2.sh
```

Then follow the instructions :

```
Enter your input folder name, followed by [ENTER] :
```

```
Enter your output folder name, followed by [ENTER] :
```

```
Does the sequencing of your data is single end (SE) or paired end (PE) ? [SE|PE] :
```

```
Type the number of thread(s) you want to use for trimmomatic, followed by [ENTER] :
```

```
Does your dataset encoding is phred33 (Sanger, Illumina >= 1.8) or phred64 (Illumina < 1.8) ? [33|64] :
```

```
Type the minimum score you want to use for trimmomatic LEADING option, followed by [ENTER] :
```

```
Type the minimum score you want to use for trimmomatic TRAILING option, followed by [ENTER] :
```

```
Type the minimum length you want to use for trimmomatic MINLEN option, followed by [ENTER] :
```

```
Do you want to remove adapters from dataset ? [yes|no] :
```

```
Type the path to the adapter file ex: /home/App/Trimmomatic-0.33/adapters/TruSeq3-SE.fa, followed by [ENTER] :
```

```
Do you want to change the ILLUMINACLIP settings ? (Default : ILLUMINACLIP:adapterFile:2:30:10) [yes|no] :
```

```
Running Trimmomatic trimming for file SRRXXXXXXX.fastq.gz...
Quality check for file SRRXXXXXXX.fastq.gz...
Quality check done for file SRRXXXXXXX.fastq.gz !
```

When the file(s) is/are processed the following message will appear :

```
Trimmomatic trimming and FastQC analysis done in X hours X minutes X seconds !

If the quality check analysis is not right, you can re-run this script with other parameters ! :)
--------------------------------------------------------
[TIP] : You will may want to change the minimum score for LEADING or TRAILING option.
[TIP] : You will may want to change the Adapter file.
[TIP] : You will may want to change the ILLUMINACLIP settings.
         - min score for palindrome match should be in the range of 30 (up to 50)
         - min score for simple match should be between 7 and 15
```

###### alignPublicData

```
./alignPublicData.sh
```

```

```

```

```

```

```

```

```

```

```

```

```


###### processWGBSdata



## OUTPUT FILES



## HOW IT WORKS


