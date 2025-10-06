# 0. Unzip inputs
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

# 1. Protein prediction
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

# 2. HMM search for USCGs and MCPs
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

# 3. Extracting sequences of identified marker genes (keeping only the best hit for MCPs)
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
        fna=os.path.join(RESULTS_DIR, "hmm_sequences/{marker}_genes_{sample}.fna")
    conda: os.path.join(ENV_DIR, "annotation.yaml")
    shell:
        "seqtk subseq {input.genes} {input.ids} | seqkit rmdup -n - > {output.fna}"

# 4. Subsample reads and map to extracted genes
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
        ref=os.path.join(RESULTS_DIR, "hmm_sequences/{marker}_genes_{sample}.fna"),
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

# 5. Calculate coverage of marker genes
rule coverm_marker:
    input:
        bam=os.path.join(RESULTS_DIR, "mapping/minimap2_{sample}_mapped_on_uscg.bam")
    output:
        out=os.path.join(RESULTS_DIR, "coverm/uscg_coverm_mean_{sample}.out")
    params:
        min_read_identity=95,
        min_read_coverage=75
    conda: os.path.join(ENV_DIR, "mapping.yaml")
    threads: config['coverm']['threads']
    shell:
        "coverm contig -t {threads} -b {input.bam} -m 'mean' --min-read-percent-identity {params.min_read_identity} --min-read-aligned-percent {params.min_read_coverage} -o {output.out}"

# 6. Calculate counts per sample
rule counts_per_category_uscg:
    input:
        tbl=os.path.join(RESULTS_DIR, "hmm/uscg_hits_{sample}.tbl"),
        coverm=os.path.join(RESULTS_DIR, "coverm/uscg_coverm_mean_{sample}.out")
    output:
        hits=temporary(os.path.join(RESULTS_DIR, "headerless_uscg_hits_{sample}.tsv")),
        merged=temporary(os.path.join(RESULTS_DIR, "merged_{sample}.tsv")),
        counts=os.path.join(RESULTS_DIR, "counts/counts_per_category_uscg_{sample}.tsv")
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

rule counts_mcp:
    input:
        coverm=os.path.join(RESULTS_DIR, "coverm/mcp_coverm_mean_{sample}.out")
    output:
        counts=os.path.join(RESULTS_DIR, "counts/counts_mcp_{sample}.tsv")
    shell:
        r"""
        awk 'NR>1 {{sum+=$2}} END {{print "{wildcards.sample}\t" sum}}' {input.coverm} > {output.counts}
        """

# 7. Summarize statistics across samples
rule summary_stats_uscg:
    input:
        counts=expand(os.path.join(RESULTS_DIR, "counts/counts_per_category_uscg_{sample}.tsv"), sample=SAMPLES)
    output:
        summary=os.path.join(RESULTS_DIR, "summary_stats_uscg.tsv")
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
        counts=expand(os.path.join(RESULTS_DIR, "counts/counts_mcp_{sample}.tsv"), sample=SAMPLES)
    output:
        summary=os.path.join(RESULTS_DIR, "summary_stats_mcp.tsv")
    shell:
        r"""
        echo -e "Sample\tSumMCP" > {output.summary}
        cat {input.counts} >> {output.summary}
        """

rule vmr_calculation:
    input:
        uscg_summary=os.path.join(RESULTS_DIR, "summary_stats_uscg.tsv"),
        mcp_summary=os.path.join(RESULTS_DIR, "summary_stats_mcp.tsv")
    output:
        vmr=os.path.join(RESULTS_DIR, "vmr_results.tsv")
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
