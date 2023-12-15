#!/bin/bash  -l 

#file_names=()
#files_dir="/rhome/zli529/lab/RNA_machine_learning"

#cd $files_dir

#for names in $(cat sample_names.txt);do
#    file_names+=("$names")
#done


public_dir="/rhome/zli529/lab/RNA_machine_learning/public_data/pairEnd_data"
#cd $public_dir

singleEnd_data="/rhome/zli529/lab/RNA_machine_learning/singleEnd_data/pairEnd_data2"
ENA_single="/bigdata/lerochlab/zli529/RNA_machine_learning/ENA_RNA/single"
#cd $singleEnd_data
single_list="ENAsingle_list.txt"
for names in $ENA_single/*;do
    echo "$names" >> "$single_list"
done


#mkdir $singleEnd_data

#for names in $(cat sample_names.txt); do
    
#    if [ -e "${names}_2.part_.fastq" ] ; then
#        echo "${names}_2.fastq exists."
#    else 
#        mv "${names}_1.fastq" "$singleEnd_data"
#    fi
#done

#mkdir "fastq1"
#mkdir "fastq2"

#fastq1_list="ENA1_list.txt"
#fastq2_list="ENA2_list.txt"
#for names in /bigdata/lerochlab/zli529/RNA_machine_learning/ENA_RNA/pairEnd_data/*"1.part"*;do
#    echo "$names"  >> "$fastq1_list" 
#done

#for names in /bigdata/lerochlab/zli529/RNA_machine_learning/ENA_RNA/pairEnd_data/*"2.part"*;do
#    echo "$names"  >> "$fastq2_list" 
#done

