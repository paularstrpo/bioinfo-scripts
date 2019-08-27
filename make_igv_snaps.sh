#!/bin/bash


bam_path='/path/to/your/bam_files/'
igv_snapshotter='/path/to/IGV-snapshot-automator/make_IGV_snapshots.py'
wdir='/path/to/your/working_directory'
ref_genome="hg38"

# only output the scripts
for i in $(ls raw_data/readcounts/dna/*.bed); do
    tumor=${i##*/}
    tumor=${tumor%%_vs*}
    normal=${i##*vs_}
    normal=${normal%_*bed}
    tag=${tumor}_v_${normal}

    mkdir -p ${wdir}/${tag}

    python $igv_snapshotter -o ${wdir}/${tag} -mem 4000 -r $i -nf4 -ht 1000 -nosnap -g ${ref_genome} ${bam_path}/${tumor}/${tumor}.sorted.markdup.recal.bam  ${bam_path}/${normal}/${normal}.sorted.markdup.recal.bam

    mv ${wdir}/${tag}/IGV_snapshots.bat ${wdir}/${tag}_IGV_snapshots.bat 
    rm -r ${wdir}/${tag}/
done
