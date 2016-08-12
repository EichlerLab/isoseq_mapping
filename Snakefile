import operator
import os

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

shell.prefix("source config.sh;")
configfile: "config.json"

SPECIES = config["species"].split(",")

def _get_isoseq_file_by_name(wildcards):
    return config[wildcards.species]['isoseq']

rule collect_isoseq_alignments:
    input:
        expand("{species}.gmap.filtered.sort.split.novel_exons.tbl", species=SPECIES),
        expand("{species}.gmap.filtered.bed", species=SPECIES)

rule generate_summary_table:
    input: "{file}.bed"
    output: "{file}.tbl"
    params: sge_opts="-l mfree=4G"
    shell:
        """cat {input} | \
        awk '$4>10' | \
        bedtools coverage -a - -b /net/eichler/vol2/eee_shared/assemblies/hg38/wgac/wgac_merged.bed | \
        bedtools coverage -a - -b /net/eichler/vol2/eee_shared/assemblies/hg38/repeats_and_trf.bed | \
        cut -f 1,2,3,4,5,9,13 | \
        bedtools sort | \
        bedtools intersect -a - -b /net/eichler/vol2/eee_shared/assemblies/hg38/genes/refGene.bed12 -wao | \
        cut -f 1-7,11 | \
        bedtools groupby -i - -g 1,2,3,4,5,6,7 -c 8 -o distinct | \
        bedtools sort | \
        awk 'BEGIN {{print "#chr\tstart\tend\tisoseq_hits\tlength\tsegdup_frac\trepeat_frac\tgene"}} {{print}}' > {output}"""

rule find_novel_exons:
    input: "{file}.bed", reference_exons=config[config["reference"]]["exons"]
    output: "{file}.novel_exons.bed",
    params: sge_opts="-l mfree=16G", min_exon_length="10"
    shell:
        "bedtools intersect -v -a {input[0]} -b {input.reference_exons} | bedtools merge -c 1 -o count | sort -nrk4 | awk -v OFS=\"\t\" '{{if ($3-$2 >= {params.min_exon_length}) print $0, $3-$2}}' | bedtools sort > {output}"

rule split_and_sort_bam_to_bed:
    input: "{file}.bam"
    output: "{file}.split.bed"
    params: sge_opts="-l mfree=16G"
    shell:
        "bedtools bamtobed -split -i {input} | bedtools sort > {output}"

rule sort_bam_to_bed:
    input: "{file}.bam"
    output: "{file}.bed"
    params: sge_opts="-l mfree=16G"
    shell:
        "bedtools bamtobed -i {input} -bed12 -color | bedtools sort > {output}"

rule index_bam:
    input: "{file}.bam"
    output: "{file}.bam.bai"
    params: sge_opts="-l mfree=16G"
    shell:
        "samtools index {wildcards.file}.bam"

rule sort_bam:
    input: "{file}.bam"
    output: "{file}.sort.bam"
    params: sge_opts="-l mfree=16G"
    shell:
        "samtools sort {input} {wildcards.file}.sort"

rule sam_to_bam:
    input: "{file}.sam"
    output: "{file}.bam"
    params: sge_opts="-l mfree=32G"
    shell:
        "samtools view -bS {input} > {wildcards.file}.bam;"

rule filter_isoseq_alignments:
    input: "{alignment}.sam"
    output: "{alignment}.filtered.sam"
    params: sge_opts="-l mfree=32G"
    shell:
        """samtools view -h {input} | awk '{{if ($0 ~ /^@/) {{print}} else {{for (i=1; i<=NF; i+=1) {{if ($i ~ "NM:i:") {{split($i, a, ":"); if (a[3]/length($10) < 0.2) print}}}}}};  }}' | samtools view -hS - > {output} """

rule align_isoseq_with_gmap:
    input: isoseq=_get_isoseq_file_by_name
    output: "{species}.gmap.sam"
    params: sge_opts="-l mfree=4G -pe serial 12 -l disk_free=5G", index=config[config["reference"]]["gmap_db"], index_name=config[config["reference"]]["gmap_db_name"]
    shell:
        "rsync -r --bwlimit=50000 {params.index} /var/tmp/{params.index_name};"
        "/net/eichler/vol5/home/mchaisso/software/bin/gmap -D /var/tmp/{params.index_name}/ -d {params.index_name} -f samse -n 0 -t 12 {input.isoseq} > /var/tmp/isoseq.gmap.sam;"
        "mv /var/tmp/isoseq.gmap.sam {output};"

rule align_isoseq_with_STAR:
    input: isoseq=_get_isoseq_file_by_name
    output: "{species}.star.sam"
    params: sge_opts="-l mfree=4G -pe serial 4", index=config[config["reference"]]["star_index"]
    shell: """~chrismh/src/STAR/bin/Linux_x86_64/STARlong \
        --runMode alignReads \
        --outSAMattributes NH HI NM MD \
        --readNameSeparator space \
        --outFilterMultimapScoreRange 1 \
        --outFilterMismatchNmax 2000 \
        --scoreGapNoncan -20 \
        --scoreGapGCAG -4 \
        --scoreGapATAC -8 \
        --scoreDelOpen -1 \
        --scoreDelBase -1 \
        --scoreInsOpen -1 \
        --scoreInsBase -1 \
        --alignEndsType Local \
        --seedSearchStartLmax 50 \
        --seedPerReadNmax 100000 \
        --seedPerWindowNmax 1000 \
        --alignTranscriptsPerReadNmax 100000 \
        --alignTranscriptsPerWindowNmax 10000 \
        --runThreadN 15 \
        --genomeDir {params.index} --readFilesIn {input.isoseq} > {output}"""
