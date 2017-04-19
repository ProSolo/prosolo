#!/usr/bin/env bash

SEQS=( ref alt )
COVS=( 10 20 30 40 50 60 70 80 90 100 0 )

READ_LN=10
INS_SIZE=12

REF="ref.fa"

bwa index ${REF}
samtools faidx ${REF}

for seq in ${SEQS[@]}
do
    for cov in ${COVS[@]}
    do
        OUT_PREFIX="${seq}.${cov}X"
        TMP_PREFIX="${OUT_PREFIX}.100"
        if [ ${cov} == 0 ]
        then
            samtools view -O bam -H -o ${OUT_PREFIX}.bam ${seq}.10X.bam
            samtools index ${OUT_PREFIX}.bam
        else
            # sequences have different starting letters, making for individual seed factors
            #SEED_FACTOR=`printf "%d" "'${seq:0:1}"`
            SEED_FACTOR=1
            art_illumina -ss HS25 --noALN --rndSeed $((${cov}*${SEED_FACTOR})) -i ${seq}.fa --paired -l $READ_LN --sdev 1 -m $INS_SIZE -f 100 --out ${TMP_PREFIX}.
            bwa aln ${REF} ${TMP_PREFIX}.1.fq >${TMP_PREFIX}.1.sai
            bwa aln ${REF} ${TMP_PREFIX}.2.fq >${TMP_PREFIX}.2.sai
            bwa sampe ${REF} ${TMP_PREFIX}.1.sai ${TMP_PREFIX}.2.sai ${TMP_PREFIX}.1.fq ${TMP_PREFIX}.2.fq | samtools sort -O bam -o ${TMP_PREFIX}.bam
            rm ${TMP_PREFIX}.1.sai ${TMP_PREFIX}.2.sai ${TMP_PREFIX}.1.fq ${TMP_PREFIX}.2.fq
            samtools index ${TMP_PREFIX}.bam
            if [ $cov == 100 ]
            then
                cp ${TMP_PREFIX}.bam ${OUT_PREFIX}.bam
            else
                ZERO_PAD=`printf "%02d" ${cov}`
                samtools view -O bam -s ${cov}.${ZERO_PAD} -o ${OUT_PREFIX}.bam ${TMP_PREFIX}.bam
            fi
            samtools index ${OUT_PREFIX}.bam
            rm ${TMP_PREFIX}.bam ${TMP_PREFIX}.bam.bai
        fi
    done
    for x in 2 4 6
    do
        ln -s ${seq}.10X.bam ${seq}.${x}X.bam
        samtools index ${seq}.${x}X.bam
    done
    for x in 14 16
    do
        ln -s ${seq}.20X.bam ${seq}.${x}X.bam
        samtools index ${seq}.${x}X.bam
    done
    for x in 24 26 28
    do
        ln -s ${seq}.30X.bam ${seq}.${x}X.bam
        samtools index ${seq}.${x}X.bam
    done
done

while read -a S
do
    if [ ${S[0]:0:1} == "#" ]
    then
        continue
    fi
    # position at reference
    P=${S[1]}
    POS=${P}-${P}
    # reference coverage bulk
    RB=${S[4]}
    # alt coverage bulk
    AB=${S[5]}
    # reference coverage single cell
    RS=${S[6]}
    # alt coverage single cell
    AS=${S[7]}
    samtools merge -c -p --reference ${REF} -R ref:${POS} -O bam bulk.pos${P}.bam ref.${RB}X.bam alt.${AB}X.bam
    if [ $((${RS}+${AS})) == 10 ]
    then
        samtools merge -c -p --reference ${REF} -R ref:${POS} -O bam tmp10.bam ref.$((${RS}*10))X.bam alt.$((${AS}*10))X.bam
        samtools view -s $((${RB}+${AS})).1 -O bam -o single-cell.pos${P}.bam tmp10.bam
        rm tmp10.bam
    else
        samtools merge -c -p --reference ${REF} -R ref:${POS} -O bam tmp50.bam ref.$((${RS}*2))X.bam alt.$((${AS}*2))X.bam
        samtools view -s $((${AB}+${RS})).5 -O bam -o single-cell.pos${P}.bam tmp50.bam
        rm tmp50.bam
    fi
done < simulation-input.tsv
rm *X.bam *X.bam.bai

samtools merge -c -p --reference ${REF} -O bam bulk.bam bulk.pos*.bam
samtools index bulk.bam
rm bulk.pos*.bam
samtools merge -c -p --reference ${REF} -O bam single-cell.bam single-cell.pos*.bam
samtools index single-cell.bam
rm single-cell.pos*.bam

samtools mpileup -A -Q 0 -f ref.fa -l pos.tsv single-cell.bam | awk 'BEGIN { OFS="\t" } { print $1,$2,$5 }' | sed -r -e 's/\$|\^[^ ]//g' | awk 'BEGIN { OFS="\t" } { print $1,$2,split($3,ref,/[,.]/)-1,split($3,alt,/[^,.]/)-1 }' >single-cell.alt-ref-cov.tsv
samtools mpileup -A -Q 0 -l pos.tsv single-cell.bam | awk 'BEGIN { OFS="\t" } { print $1,$2,$5 }' | sed -r -e 's/\$|\^[^ ]//g' | awk 'BEGIN { OFS="\t" } { print $1,$2,split($3,A,/[Aa]/)-1,split($3,C,/[Cc]/)-1,split($3,G,/[Gg]/)-1,split($3,T,/[Tt]/)-1 }' >single-cell.nt-cov.tsv
samtools mpileup -A -Q 0 -f ref.fa -l pos.tsv bulk.bam | awk 'BEGIN { OFS="\t" } { print $1,$2,$5 }' | sed -r -e 's/\$|\^[^ ]//g' | awk 'BEGIN { OFS="\t" } { print $1,$2,split($3,ref,/[,.]/)-1,split($3,alt,/[^,.]/)-1 }' >bulk.alt-ref-cov.tsv
samtools mpileup -A -Q 0 -l pos.tsv bulk.bam | awk 'BEGIN { OFS="\t" } { print $1,$2,$5 }' | sed -r -e 's/\$|\^[^ ]//g' | awk 'BEGIN { OFS="\t" } { print $1,$2,split($3,A,/[Aa]/)-1,split($3,C,/[Cc]/)-1,split($3,G,/[Gg]/)-1,split($3,T,/[Tt]/)-1 }' >bulk.nt-cov.tsv

samtools mpileup -f ref.fa -l pos.tsv -Q 0 --BCF -o candidates.bcf bulk.bam single-cell.bam
