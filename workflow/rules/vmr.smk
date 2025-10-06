# Unzip inputs
rule unzip_assembly:
    input:
        fa=lambda wildcards: next(
            f for f in [
                os.path.join(ASSEMBLY_DIR, f"{wildcards.sample}.fa.gz"),
                os.path.join(ASSEMBLY_DIR, f"{wildcards.sample}.fna.gz"),
                os.path.join(ASSEMBLY_DIR, f"{wildcards.sample}.fasta.gz")
            ] if os.path.exists(f)
        )
    output: fa=os.path.join(ASSEMBLY_DIR, "{sample}.fna")
    shell:
        "gunzip -c {input.fa} > {output.fa}"

rule normalize_extension_assembly:
    input:
        lambda wildcards: next(
            f for f in [
                os.path.join(ASSEMBLY_DIR, f"{wildcards.sample}.fa"),
                os.path.join(ASSEMBLY_DIR, f"{wildcards.sample}.fasta")
            ] if os.path.exists(f)
        )
    output: os.path.join(ASSEMBLY_DIR, "{sample}.fna")
    shell:
        "cp {input} {output}"

rule unzip_reads:
    input: os.path.join(READS_DIR, "{sample}_{number}.fastq.gz")
    output: os.path.join(READS_DIR, "{sample}_{number}.fastq")
    shell:
        "gunzip -c {input} > {output}"

# Protein prediction
rule pyrodigal:
    input:
        fa=os.path.join(ASSEMBLY_DIR, "{sample}.fna")
    output:
        proteins=os.path.join(RESULTS_DIR, "proteins/proteins_{sample}.faa"),
        genes=os.path.join(RESULTS_DIR, "genes/genes_{sample}.fna")
    conda: os.path.join(ENV_DIR, "annotation.yaml")
    threads: config['pyrodigal']['threads']
    shell:
        "pyrodigal -j {threads} -i {input.fa} -a {output.proteins} -d {output.genes} -p meta"

# HMM search for USCGs and MCPs
rule hmmsearch_uscg:
    input:
        proteins=os.path.join(RESULTS_DIR, "proteins/proteins_{sample}.faa"),
        uscg="../data/all_USCG.hmm"
    output:
        tbl=os.path.join(RESULTS_DIR, "hmm/uscg_hits_{sample}.tbl"),
        txt=os.path.join(RESULTS_DIR, "hmm/uscg_hmmsearch_{sample}.txt")
    conda: os.path.join(ENV_DIR, "annotation.yaml")
    threads: config['hmmsearch']['threads']
    shell:
        "hmmsearch --cut_ga --cpu {threads} --tblout {output.tbl} {input.uscg} {input.proteins} > {output.txt}"

rule hmmsearch_mcp:
    input:
        proteins=os.path.join(RESULTS_DIR, "proteins/proteins_{sample}.faa"),
        mcp="../data/MCPs_MK_Selection.hmm"
    output:
        tbl=os.path.join(RESULTS_DIR, "hmm/mcp_hits_{sample}.tbl"),
        txt=os.path.join(RESULTS_DIR, "hmm/mcp_hmmsearch_{sample}.txt")
    conda: os.path.join(ENV_DIR, "annotation.yaml")
    threads: config['hmmsearch']['threads']
    shell:
        "hmmsearch -E 1e-11 --cpu {threads} --tblout {output.tbl} {input.mcp} {input.proteins} > {output.txt}"

# Extracting sequences of identified marker genes (keeping only the best hit for MCPs)
rule extract_ids_uscg:
    input:
        tbl=os.path.join(RESULTS_DIR, "hmm/uscg_hits_{sample}.tbl")
    output:
        ids=os.path.join(RESULTS_DIR, "hmm/uscg_hits_{sample}_ids.txt")
    shell:
        "grep -v '^#' {input.tbl} | awk '{{print $1}}' > {output.ids}"

rule filter_best_mcp:
    input:
        tbl=os.path.join(RESULTS_DIR, "hmm/mcp_hits_{sample}.tbl")
    output:
        all_ids=temporary(os.path.join(RESULTS_DIR, "hmm/all_mcp_hits_{sample}_ids.txt")),
        ids=os.path.join(RESULTS_DIR, "hmm/mcp_hits_{sample}_ids.txt")
    shell:
        r"""
        grep -v '^#' {input.tbl} | awk '{{print $1}}' > {output.all_ids} && \
        awk '
        {{
            n = split($0, parts, "_");
            suffix = parts[n];
            base_id = "";
            for (i = 1; i < n; i++) {{
                base_id = base_id (i == 1 ? "" : "_") parts[i];
            }}
            if (!(base_id in seen) || suffix < seen[base_id]) {{
                seen[base_id] = suffix;
                full_id[base_id] = $0;
            }}
        }}
        END {{
            for (id in full_id)
                print full_id[id];
        }}
        ' {output.all_ids} > {output.ids}
        """

rule extract_sequences_marker:
    input:
        genes=os.path.join(RESULTS_DIR, "genes/genes_{sample}.fna"),
        ids=os.path.join(RESULTS_DIR, "hmm/{marker}_hits_{sample}_ids.txt")
    output:
        fna=os.path.join(RESULTS_DIR, "hmm/{marker}_genes_{sample}.fna")
    conda: os.path.join(ENV_DIR, "annotation.yaml")
    shell:
        "seqtk subseq {input.genes} {input.ids} | seqkit rmdup -n - > {output.fna}"

# Find viral contigs within the assembly
rule filter_length_2000:
    output: os.path.join(RESULTS_DIR, "genomad", "filtered2000_{sample}.fna")
    input: os.path.join(ASSEMBLY_DIR, "{sample}.fna")
    conda: os.path.join(ENV_DIR, "annotation.yaml")
    threads: config['seqkit']['threads']
    shell:
        """(date && seqkit seq -j {threads} --remove-gaps -o {output} -m 2000 {input} && date) &> {log}"""

rule db_genomad:
    output: os.path.join(RESULTS_DIR, "dbs", "genomad_db", "genomad_marker_metadata.tsv")
    conda: os.path.join(ENV_DIR, "genomad.yaml")
    shell:
        """
        (date && cd $(dirname $(dirname {output})) &&
        wget -nc https://zenodo.org/records/14886553/files/genomad_db_v1.9.tar.gz && 
        tar --skip-old-files -zxvf genomad_db_v1.9.tar.gz && date) &> {log}
        """

rule genomad:
    output: os.path.join(RESULTS_DIR, "genomad", "geNomad_{sample}.filtered2000.renamed.contigs", "{sample}.filtered2000.renamed.contigs_summary", "{sample}.filtered2000.renamed.contigs_virus.fna")
    input: 
        assembly = os.path.join(RESULTS_DIR, "genomad", "filtered2000_{sample}.fna"),
        db = os.path.join(RESULTS_DIR, "dbs", "genomad_db", "genomad_marker_metadata.tsv")
    conda: os.path.join(ENV_DIR, "genomad.yaml")
    threads: config['genomad']['threads']
    shell:
        "(date && genomad end-to-end --threads {threads} --enable-score-calibration --max-fdr 0.05 {input.assembly} $(dirname $(dirname {output})) $(dirname {input.db}) && date) &> {log}"

rule extract_viruses:
    input:
        genomad = os.path.join(RESULTS_DIR, "viruses", "genomad", "geNomad_{sample}.filtered2000.renamed.contigs", "{sample}.filtered2000.renamed.contigs_summary", "{sample}.filtered2000.renamed.contigs_virus.fna")
    output:
        fna=os.path.join(RESULTS_DIR, "sequences/viruses_genes_{sample}.fna")
    conda: os.path.join(ENV_DIR, "annotation.yaml")
    shell:
        "ln -s {input.genomad} > {output.fna}"

# Subsample reads and map to extracted genes
rule subsample_reads:
    input:
        fq=os.path.join(READS_DIR, "{sample}_{number}.fastq")
    output:
        fq_sub=os.path.join(RESULTS_DIR, "reads_subsampled/subsample_{sample}_{number}.fastq")
    conda: os.path.join(ENV_DIR, "annotation.yaml")
    threads: config['seqkit']['threads']
    shell:
        "seqkit sample -j {threads} -p 0.1 -s 42 -o {output.fq_sub} {input.fq}"

rule minimap2_marker_mapping:
    input:
        ref=os.path.join(RESULTS_DIR, "sequences/{marker}_genes_{sample}.fna"),
        r1=os.path.join(RESULTS_DIR, "reads_subsampled/subsample_{sample}_R1.fastq"),
        r2=os.path.join(RESULTS_DIR, "reads_subsampled/subsample_{sample}_R2.fastq")
    output:
        bam=os.path.join(RESULTS_DIR, "mapping/minimap2_{sample}_mapped_on_{marker}.bam")
    conda: os.path.join(ENV_DIR, "mapping.yaml")
    threads: config['minimap2']['threads']
    shell:
        "minimap2 -ax sr -I 16G -t {threads} {input.ref} {input.r1} {input.r2} | "
        "samtools view -@ {threads} -bS -F 4 | "
        "samtools sort -@ {threads} -o {output.bam}"

# Calculate coverage of marker genes
rule coverm_filter_marker:
    input:
        bam=os.path.join(RESULTS_DIR, "mapping/minimap2_{sample}_mapped_on_{marker}.bam")
    output:
        out=os.path.join(RESULTS_DIR, "mapping/minimap2_{sample}_mapped_on_{marker}_filtered.bam")
    params:
        min_read_identity=90,
        min_read_coverage=75
    conda: os.path.join(ENV_DIR, "mapping.yaml")
    threads: config['coverm']['threads']
    shell:
        "coverm filter -t {threads} -b {input.bam} --min-read-percent-identity {params.min_read_identity} --min-read-aligned-percent {params.min_read_coverage} -o {output.out}"

rule coverm_mean_marker:
    input:
        bam=os.path.join(RESULTS_DIR, "mapping/minimap2_{sample}_mapped_on_{marker}_filtered.bam")
    output:
        out=os.path.join(RESULTS_DIR, "coverm/{marker}_coverm_mean_{sample}.out")
    conda: os.path.join(ENV_DIR, "mapping.yaml")
    threads: config['coverm']['threads']
    shell:
        "coverm contig -t {threads} -b {input.bam} -m 'mean' -o {output.out}"

rule coverm_length_marker:
    input:
        bam=os.path.join(RESULTS_DIR, "mapping/minimap2_{sample}_mapped_on_{marker}_filtered.bam")
    output:
        out=os.path.join(RESULTS_DIR, "coverm/{marker}_coverm_length_{sample}.out")
    conda: os.path.join(ENV_DIR, "mapping.yaml")
    threads: config['coverm']['threads']
    shell:
        "coverm contig -t {threads} -b {input.bam} -m 'length' -o {output.out}"

rule coverm_count_marker:
    input:
        bam=os.path.join(RESULTS_DIR, "mapping/minimap2_{sample}_mapped_on_{marker}_filtered.bam")
    output:
        out=os.path.join(RESULTS_DIR, "coverm/{marker}_coverm_count_{sample}.out")
    conda: os.path.join(ENV_DIR, "mapping.yaml")
    threads: config['coverm']['threads']
    shell:
        "coverm contig -t {threads} -b {input.bam} -m 'count' -o {output.out}"

rule calculate_rpkm_marker:
    input:
        length_file=os.path.join(RESULTS_DIR, "coverm/{marker}_coverm_length_{sample}.out"),
        count_file=os.path.join(RESULTS_DIR, "coverm/{marker}_coverm_count_{sample}.out"),
        r1_file=os.path.join(RESULTS_DIR, "reads_subsampled/subsample_{sample}_R1.fastq"),
        r2_file=os.path.join(RESULTS_DIR, "reads_subsampled/subsample_{sample}_R2.fastq")
    output: os.path.join(RESULTS_DIR, "coverm/{marker}_rpkm_{sample}.tsv")
    shell:
        r"""
        total_reads=$(( $(wc -l < {input.r1_file})/4 + $(wc -l < {input.r2_file})/4 ))
        paste {input.length_file} {input.count_file} | awk -v total_reads=$total_reads 'BEGIN{{OFS="\t"}} {{
            if($1==$3) {{ rpkm = ($4 * 1e9)/($2 * total_reads); print $1,$2,$4,rpkm }}
        }}' > {output}
        """

rule coverm_rpkm_marker:
    input: os.path.join(RESULTS_DIR, "coverm/{marker}_rpkm_{sample}.tsv")
    output: os.path.join(RESULTS_DIR, "coverm/{marker}_coverm_rpkm_{sample}.out")
    shell:
        r"""
        cut -f 1,4 {input} > {output}
        """

# Calculate counts per sample
rule counts_per_category_uscg:
    input:
        tbl=os.path.join(RESULTS_DIR, "hmm/uscg_hits_{sample}.tbl"),
        coverm=os.path.join(RESULTS_DIR, "coverm/uscg_coverm_{method}_{sample}.out")
    output:
        hits=temporary(os.path.join(RESULTS_DIR, "headerless_uscg_hits_{sample}_with_{method}.tsv")),
        merged=temporary(os.path.join(RESULTS_DIR, "merged_{sample}_with_{method}.tsv")),
        counts=os.path.join(RESULTS_DIR, "counts/counts_per_category_uscg_{sample}_with_{method}.tsv")
    shell:
        r"""
        grep -v '^#' {input.tbl} | awk '{{print $1, $3}}' > {output.hits} && \
        awk 'NR==FNR{{a[$1]=$2; next}} NR==1{{print $0, "Category"; next}} {{print $0, a[$1]}}' {output.hits} {input.coverm} > {output.merged} && \
        awk '
        BEGIN {{
            OFS="\t"
            map["Ribosomal_L24e"]="ribosomal_L24"
            map["Ribosomal_S17e"]="Ribosomal_S17"
            map["Ribosomal_S19e"]="Ribosomal_S19"
            map["Ribosomal_S3Ae"]="Ribosomal_S3_C"
            map["Ribosomal_S8e"]="Ribosomal_S8"
        }}
        NR==1 {{ next }}
        {{
            cat=$3
            if (cat in map) cat=map[cat]
            sum[cat]+=$2
        }}
        END {{
            print "Category","Sum"
            for (c in sum) print c, sum[c]
        }}' {output.merged} > {output.counts}
        """

rule counts_virus:
    input:
        coverm=os.path.join(RESULTS_DIR, "coverm/{marker}_coverm_{method}_{sample}.out")
    output:
        counts=os.path.join(RESULTS_DIR, "counts/counts_{marker}_{sample}_with_{method}.tsv")
    shell:
        r"""
        awk 'NR>1 {{sum+=$2}} END {{print "{wildcards.sample}\t" sum}}' {input.coverm} > {output.counts}
        """

# Summarize statistics across samples
rule summary_stats_uscg:
    input:
        counts=expand(os.path.join(RESULTS_DIR, "counts/counts_per_category_uscg_{sample}_with_{{method}}.tsv"), sample=SAMPLES)
    output:
        summary=os.path.join(RESULTS_DIR, "summary_stats_hosts_with_uscg_and_{method}.tsv")
    shell:
        r"""
        echo -e "Sample\tMaxUSCG\tMeanUSCG\tMedianUSCG" > {output.summary}
        for sample in {input.counts}; do
            name=$(basename ${{sample}} | sed 's/counts_per_category_uscg_//g' | sed 's/.tsv//g')
            vals=($(awk 'NR>1 {{print $2}}' ${{sample}} | sort -n))
            n=${{#vals[@]}}
            if [ $n -gt 0 ]; then
                sum=0
                for v in "${{vals[@]}}"; do
                    sum=$(echo "$sum + $v" | bc)
                done
                mean=$(echo "scale=4; $sum / $n" | bc)
                maxv=${{vals[$((n-1))]}}
                if (( n % 2 == 1 )); then
                    median=${{vals[$((n/2))]}}
                else
                    median=$(echo "scale=4; (${{vals[$((n/2 - 1))]}} + ${{vals[$((n/2))]}}) / 2" | bc)
                fi
            else
                mean=0
                maxv=0
                median=0
            fi
            echo -e "${{name}}\t${{maxv}}\t${{mean}}\t${{median}}" >> {output.summary}
        done
        """

rule summary_stats_mcp:
    input:
        counts=expand(os.path.join(RESULTS_DIR, "counts/counts_{{marker}}_{sample}_with_{{method}}.tsv"), sample=SAMPLES)
    output:
        summary=os.path.join(RESULTS_DIR, "summary_stats_virus_with_{marker}_and_{method}.tsv")
    shell:
        r"""
        echo -e "Sample\tSumMCP" > {output.summary}
        cat {input.counts} >> {output.summary}
        """

rule vmr_calculation_with_mcp:
    input:
        uscg_summary=os.path.join(RESULTS_DIR, "summary_stats_hosts_with_uscg_and_{method}.tsv"),
        mcp_summary=os.path.join(RESULTS_DIR, "summary_stats_virus_with_mcp_and_{method}.tsv")
    output:
        vmr=os.path.join(RESULTS_DIR, "vmr_results_with_mcp_and_{method}.tsv")
    shell:
        r"""
        echo -e "Sample\tSumMCP\tMaxUSCG\tMeanUSCG\tMedianUSCG\tVMR_MaxUSCG\tVMR_MeanUSCG\tVMR_MedianUSCG" > {output.vmr}
        join -1 1 -2 1 {input.mcp_summary} {input.uscg_summary} | awk '
            BEGIN {{ OFS="\t" }}
            NR>1 {{
                sample=$1
                sum_mcp=$2
                max_uscg=$3
                mean_uscg=$4
                median_uscg=$5
                if (max_uscg > 0) {{ vmr_max = sum_mcp / max_uscg }} else {{ vmr_max = "NA" }}
                if (mean_uscg > 0) {{ vmr_mean = sum_mcp / mean_uscg }} else {{ vmr_mean = "NA" }}
                if (median_uscg > 0) {{ vmr_median = sum_mcp / median_uscg }} else {{ vmr_median = "NA" }}
                print sample, sum_mcp, max_uscg, mean_uscg, median_uscg, vmr_max, vmr_mean, vmr_median
            }}' >> {output.vmr}
        """

rule vmr_calculation_with_virus:
    input:
        uscg_summary=os.path.join(RESULTS_DIR, "summary_stats_uscg_with_{method}.tsv"),
        mcp_summary=os.path.join(RESULTS_DIR, "summary_stats_virus_with_{method}.tsv")
    output:
        vmr=os.path.join(RESULTS_DIR, "vmr_results_with_virus_and_{method}.tsv")
    shell:
        r"""
        echo -e "Sample\tSumMCP\tMaxUSCG\tMeanUSCG\tMedianUSCG\tVMR_MaxUSCG\tVMR_MeanUSCG\tVMR_MedianUSCG" > {output.vmr}
        join -1 1 -2 1 {input.mcp_summary} {input.uscg_summary} | awk '
            BEGIN {{ OFS="\t" }}
            NR>1 {{
                sample=$1
                sum_mcp=$2
                max_uscg=$3
                mean_uscg=$4
                median_uscg=$5
                if (max_uscg > 0) {{ vmr_max = sum_mcp / max_uscg }} else {{ vmr_max = "NA" }}
                if (mean_uscg > 0) {{ vmr_mean = sum_mcp / mean_uscg }} else {{ vmr_mean = "NA" }}
                if (median_uscg > 0) {{ vmr_median = sum_mcp / median_uscg }} else {{ vmr_median = "NA" }}
                print sample, sum_mcp, max_uscg, mean_uscg, median_uscg, vmr_max, vmr_mean, vmr_median
            }}' >> {output.vmr}
        """

# Generate the VMR results for the different methods
rule vmr_calculation_with_mcp_and_mean:
    input: os.path.join(RESULTS_DIR, "vmr_results_with_mcp_and_mean.tsv")

rule vmr_calculation_with_mcp_and_rpkm:
    input: os.path.join(RESULTS_DIR, "vmr_results_with_mcp_and_rpkm.tsv")

rule vmr_calculation_with_virus_and_mean:
    input: os.path.join(RESULTS_DIR, "vmr_results_with_virus_and_mean.tsv")

rule vmr_calculation_with_virus_and_rpkm:
    input: os.path.join(RESULTS_DIR, "vmr_results_with_virus_and_rpkm.tsv")
