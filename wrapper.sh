#!/bin/bash
# vmr: wrapper to run Snakemake from anywhere with custom directories

# Get the directory of the wrapper script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Workflow and config relative to the wrapper
WORKFLOW="$SCRIPT_DIR/workflow"
CONFIG_FILE="$SCRIPT_DIR/config/config.yaml"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --assembly_dir) assembly_dir="$2"; shift 2 ;;
        --reads_dir) reads_dir="$2"; shift 2 ;;
        --result_dir) results_dir="$2"; shift 2 ;;
        --help) echo "Usage: vmr --assembly_dir DIR --reads_dir DIR --result_dir DIR"; exit 0 ;;
        *) echo "Unknown option $1"; exit 1 ;;
    esac
done

# Update config.yaml
sed -i "s|^reads_dir: .*|reads_dir: \"$reads_dir\"|" "$CONFIG_FILE"
sed -i "s|^assembly_dir: .*|assembly_dir: \"$assembly_dir\"|" "$CONFIG_FILE"
sed -i "s|^results_dir: .*|results_dir: \"$results_dir\"|" "$CONFIG_FILE"

# Run Snakemake from anywhere
echo "$SNAKEFILE"
SMK_CLUSTER="sbatch -q {cluster.qos} -N {cluster.nodes} -n {cluster.n} -c {threads} -t {cluster.time}"
snakemake \
    --directory "$WORKFLOW" \
    --configfile "$CONFIG_FILE" \
    --use-conda \
    --jobs 50 \
    --cluster-config "config/slurm.yaml" \
    --cluster "$SMK_CLUSTER" \
    --keep-going \
    -p vmr_calculation
