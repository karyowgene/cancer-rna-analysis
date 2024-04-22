# Define source directory
source_dir = "/user/HDD/rna_seq"

rule all:
    input:
        expand("{source_dir}/multiqc_report.html", source_dir=source_dir),
        expand("kallisto_quant_{sample}/abundance.tsv", sample=[sample.split('_R1')[0] for sample in glob_wildcards(os.path.join(source_dir, "{sample}_R1_paired.fastq.gz")).sample])

rule rename_files:
    input:
        r1 = os.path.join(source_dir, "{sample_id}_R1_001.fastq.gz"),
        r2 = os.path.join(source_dir, "{sample_id}_R2_001.fastq.gz")
    output:
        r1_out = os.path.join(source_dir, "{sample_id}_R1.fastq.gz"),
        r2_out = os.path.join(source_dir, "{sample_id}_R2.fastq.gz")
    shell:
        """
        mv {input.r1} {output.r1_out}
        mv {input.r2} {output.r2_out}
        """

rule fastqc:
    input:
        fastq = os.path.join(source_dir, "{sample}_R1.fastq.gz")
    output:
        html = os.path.join(source_dir, "fastqc_{sample}_R1.html")
    shell:
        "fastqc {input.fastq} -o {source_dir}"

rule trimmomatic:
    input:
        r1 = os.path.join(source_dir, "{sample}_R1.fastq.gz"),
        r2 = os.path.join(source_dir, "{sample}_R2.fastq.gz")
    output:
        r1_paired = os.path.join(source_dir, "{sample}_R1_paired.fastq.gz"),
        r1_unpaired = os.path.join(source_dir, "{sample}_R1_unpaired.fastq.gz"),
        r2_paired = os.path.join(source_dir, "{sample}_R2_paired.fastq.gz"),
        r2_unpaired = os.path.join(source_dir, "{sample}_R2_unpaired.fastq.gz")
    params:
        adapters = "TruSeq3-PE-2.fa"
    shell:
        """
        trimmomatic PE -threads 4 {input.r1} {input.r2} \
            {output.r1_paired} {output.r1_unpaired} {output.r2_paired} {output.r2_unpaired} \
            ILLUMINACLIP:{params.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:151
        """

rule fastqc_trimmed:
    input:
        fastq = os.path.join(source_dir, "{sample}_R1_paired.fastq.gz")
    output:
        html = os.path.join(source_dir, "fastqc_{sample}_R1_paired.html")
    shell:
        "fastqc {input.fastq} -o {source_dir}"

rule multiqc:
    input:
        expand(os.path.join(source_dir, "fastqc_{sample}_R1_paired.html"), sample=glob_wildcards(os.path.join(source_dir, "{sample}_R1_paired.fastq.gz")).sample)
    output:
        html = os.path.join(source_dir, "multiqc_report.html")
    shell:
        "multiqc {source_dir} -o {source_dir}"

rule kallisto:
    input:
        r1 = os.path.join(source_dir, "{sample}_R1_paired.fastq.gz"),
        r2 = os.path.join(source_dir, "{sample}_R2_paired.fastq.gz")
    output:
        directory("kallisto_quant_{sample}")
    params:
        index = "/user/HDD2/RNA/rnaref/Homo_sapiens.GRCh38.cdna.all.release-110.idx"
    shell:
        """
        kallisto quant -i {params.index} -o kallisto_quant_{wildcards.sample} {input.r1} {input.r2}
        """
