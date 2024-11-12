file='emission_frequencies.txt'

minimap2 -a -x map-ont hybrid_ref.fasta pure_reads/EM11.fastq.gz > alignments/EM11.sam
samtools sort alignments/EM11.sam > alignments/EM11.bam
samtools index alignments/EM11.bam

minimap2 -a -x map-ont hybrid_ref.fasta pure_reads/EM60.fastq.gz > alignments/EM60.sam
samtools sort alignments/EM60.sam > alignments/EM60.bam
samtools index alignments/EM60.bam

python write_emission_frequencies.py \
            --bam alignments/EM11.bam \
            --msa_refs msa_refs.fasta \
            --name 'EM11' \
            --out $file \

python write_emission_frequencies.py \
            --bam alignments/EM60.bam \
            --msa_refs msa_refs.fasta \
            --name 'EM60' \
            --out $file \