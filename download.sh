#!/bin/bash  -l


rna_names_2=()
input_dir="/rhome/zli529/lab/RNA_machine_learning"

cd $input_dir
for names in $(cat sample_names_3.txt); do
    
    rna_names_2+=("$names")
done


echo ${rna_names_2[@]}



for file_names in ${rna_names_2[@]};do 
    cd ~/lab/RNA_machine_learning/singleEnd_data
    echo "downloading $file_names"
    ~/lab/RNA_machine_learning/sratoolkit.3.0.7-ubuntu64/bin/fastq-dump --split-files $file_names

done



