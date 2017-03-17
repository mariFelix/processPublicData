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

To run the scripts do :

```
./downloadPublicData_v2.2.sh 
```

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
```

```
Dowloading SRRXXXXXXX to fastq.gz format...
```

```
```







## OUTPUT FILES



## HOW IT WORKS


