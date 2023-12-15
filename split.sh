#!/bin/bash  -l



read1=()
read2=()


ENA_RNA="/bigdata/lerochlab/zli529/RNA_machine_learning/ENA_RNA"
#input_dir="/bigdata/lerochlab/zli529/RNA_machine_learning/public_data"
for sample in "$ENA_RNA"/*_1.fastq.gz;do
    samplename=$(basename "$sample")
    read1+=("$samplename")
done

#echo "R1 files are ${read1[@]}"

for sample in "$ENA_RNA"/*_2.fastq.gz;do
    samplename=$(basename "$sample")
    read2+=("$samplename")
done

#echo "R2 files are ${read2[@]}"

module load seqkit
#public_dir="/rhome/zli529/lab/RNA_machine_learning/public_data"
singleEnd_data="/rhome/zli529/lab/RNA_machine_learning/singleEnd_data"

#cd $public_dir
cd $ENA_RNA
for i in {1..15};do

    seqkit split2 -1 ${read1[i]} -2 ${read2[i]} -p 10 -O pairEnd_data 
done     


#for sample in "$ENA_RNA"/*_1.fastq;do
#    samplename=$(basename "$sample")
#    read1+=("$samplename")
#done
#echo "R2 files are ${read1[@]}"

#for i in {1..15};do
#    seqkit split2  ${read1[i]}  -p 10 -O pairEnd_data
#done


