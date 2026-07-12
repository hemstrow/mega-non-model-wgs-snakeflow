
rule trim_reads_pe:
    input:
        unpack(get_fastq),
    output:
        r1=temp("results/bqsr-round-{bqsr_round}/trimmed/{sample}---{unit}.1.fastq.gz"),
        r2=temp("results/bqsr-round-{bqsr_round}/trimmed/{sample}---{unit}.2.fastq.gz"),
        #r1_unpaired=temp("results/bqsr-round-{bqsr_round}/trimmed/{sample}---{unit}.1.unpaired.fastq.gz"),
        #r2_unpaired=temp("results/bqsr-round-{bqsr_round}/trimmed/{sample}---{unit}.2.unpaired.fastq.gz"),
        html="results/bqsr-round-{bqsr_round}/qc/fastp/{sample}---{unit}.html",
        json="results/bqsr-round-{bqsr_round}/qc/fastp/{sample}---{unit}.json"
    conda:
        "../envs/fastp.yaml"
    log:
        out="results/bqsr-round-{bqsr_round}/logs/trim_reads_pe/{sample}---{unit}.log",
        err="results/bqsr-round-{bqsr_round}/logs/trim_reads_pe/{sample}---{unit}.err"
    benchmark:
        "results/bqsr-round-{bqsr_round}/benchmarks/trim_reads_pe/{sample}---{unit}.bmk"
    params:
        trim_settings=config["params"]["fastp"]["pe"]["trimmer"],
    shell:
        " fastp -i {input.r1} -I {input.r2} "
        "       -o {output.r1} -O {output.r2} "
        "       -h {output.html} -j {output.json} "
        "  {params.trim_settings} > {log.out} 2> {log.err} "


# eca modified this.  The idea is to give 4 threads to bwa.
# and it will get 4 cores and also take all the memory you'd
# expect for those cores.  Sedna's machines are almost all
# 20 core units, so this should fill them up OK.
rule map_reads:
    input:
        reads = [
            "results/bqsr-round-{bqsr_round}/trimmed/{sample}---{unit}.1.fastq.gz",
            "results/bqsr-round-{bqsr_round}/trimmed/{sample}---{unit}.2.fastq.gz"
        ],
        idx=rules.bwa_index.output,
    output:
        temp("results/bqsr-round-{bqsr_round}/mapped/{sample}---{unit}.sorted.bam"),
    log:
        "results/bqsr-round-{bqsr_round}/logs/map_reads/{sample}---{unit}.log",
    benchmark:
        "results/bqsr-round-{bqsr_round}/benchmarks/map_reads/{sample}---{unit}.bmk"
    params:
        extra=get_read_group,
        sorting="samtools",
        sort_order="coordinate",
        sort_extra=""
    threads: 4
    resources:
        mem_mb=19200,
        time="23:59:59"
    wrapper:
        "v1.23.3/bio/bwa/mem"



rule filter_bams:
    input:
        get_all_bams_of_common_sample
    output:
        bam="results/bqsr-round-{bqsr_round}/mkdup/{sample}.bam",
        bai="results/bqsr-round-{bqsr_round}/mkdup/{sample}.bai",
        metrics="results/bqsr-round-{bqsr_round}/qc/mkdup/{sample}.metrics.txt",
        merged=temp("results/bqsr-round-{bqsr_round}/mkdup/{sample}.merge.bam"),
        nsort=temp("results/bqsr-round-{bqsr_round}/mkdup/{sample}.nsort.bam"),
        markdup=temp("results/bqsr-round-{bqsr_round}/mkdup/{sample}.markdup.bam"),
        psort=temp("results/bqsr-round-{bqsr_round}/mkdup/{sample}.psort.bam"),
        q1=temp("results/bqsr-round-{bqsr_round}/mkdup/{sample}.q1.bam"),
        fixmate=temp("results/bqsr-round-{bqsr_round}/mkdup/{sample}.fixmate.bam"),
        namesort=temp("results/bqsr-round-{bqsr_round}/mkdup/{sample}.namesort.bam"),
        flt=temp("results/bqsr-round-{bqsr_round}/mkdup/{sample}.flt.bam"),
    log:
        "results/bqsr-round-{bqsr_round}/logs/filter_bams/{sample}.log",
    benchmark:
        "results/bqsr-round-{bqsr_round}/benchmarks/filter_bams/{sample}.bmk"
    params:
        mapQ=config["filtering"]["bams"]["initial_mapQ"]
    conda:
        "../envs/samtools.yaml"
    shell:
        '''
        # merge and report
	echo "Beginning merge." > {log}
        samtools merge {output.merged} {input} 2>> {log}
        echo "Initial merged bam:" > {output.metrics}
	samtools flagstat {output.merged} >> {output.metrics}

	# fixmate, sort, and remove dups
	echo "Beginning Duplicate Removal." >> {log}
	echo "Namesort 1." >> {log}
	samtools sort -n -o {output.nsort} {output.merged} # sort by name
        echo "Fixmate." >> {log}
	samtools fixmate -r -m {output.nsort} {output.fixmate} 2>> {log} # fixmate
	echo "psort." >> {log}
        samtools sort -o {output.psort} {output.fixmate} 2>> {log} # sort by position
	echo "markdup." >> {log}
        samtools markdup -r {output.psort} {output.markdup} 2>> {log} # remove dups
	echo "Duplicates Removed:" >> {output.metrics}
        samtools flagstat {output.markdup} >> {output.metrics}

	# remove poorly mapped
	echo "Beginning mapping filter." >> {log}
        samtools view -q {params.mapQ} -b {output.markdup} > {output.q1} 2>> {log} # remove poorly mapped
	echo "Namesort 2." >> {log}
        samtools sort -n -o {output.namesort} {output.q1} 2>> {log} # sort by name again
	echo "Poorly Mapped Removed:" >> {output.metrics}
        samtools flagstat {output.namesort} >> {output.metrics}

	# remove improper pairs and finalize
	echo "Beginning improper pair removal." >> {log}
        samtools fixmate -m {output.namesort} {output.fixmate} 2>> {log} # filter bad mates again
	echo "Imp. Pair Filter." >> {log}
        samtools view -f 0x2 -b {output.fixmate} > {output.flt} 2>> {log} # remove improper pairs
	echo "Psort 2." >> {log}
        samtools sort -o {output.bam} {output.flt} 2>> {log} # sort by position again
        samtools index {output.bam} 2>> {log} # index
	echo "Improper Pairs Removed (final bam):" >> {output.metrics}
        samtools flagstat {output.bam} >> {output.metrics}
	echo "Done." >> {log}
        '''


