import os

# Load the config file
configfile: "config.yaml"

# Extract sample names and their corresponding fastq files
samples = list(config["samples"].keys())
indices = config["indices"]
nodes = [node for node in config["node_indices"].keys()]
node_to_indices = config["node_indices"]


# Define a dictionary to map indices to nodes
index_to_node = {}
for node, indices in node_to_indices.items():
    for index in indices:
        index_to_node[index] = node

# Define all the final outputs, including node-specific jobs and merged BAM file
rule all:
    input:
        expand(
            "results/{sample}_{index}.bam",
            sample=samples,
            index=indices,
        ),
        "results/merged.bam",


# Function to create dynamic rules for mapping reads
def create_map_reads_rule(sample, node, index):
    print(f"Creating rule for {sample} on {node} with index {index}")
    slurm_node = index_to_node[index]
    rule:
        input:
            fastq=f"data/{sample}.fastq",
            index=lambda wildcards: config["indices"][index],
        output:
            bam=f"results/{sample}_{index}.bam",
        resources:
            mem_mb=32000,
            cpus_per_task=4,
            slurm_extra=f"--nodelist={slurm_node}",
        shell:
            """
            bowtie2 -x {input.index} -U {input.fastq} -S {output.bam}
            """


# Dynamically create mapping rules for each node and index
for sample in samples:
    for index, node in index_to_node.items():
            create_map_reads_rule(sample, node, index)


# Rule to merge BAM files using samtools
rule merge_bam:
    input:
        expand(
            "results/{sample}_{index}.bam",
            sample=samples,
            index=indices,
        ),
    output:
        merged_bam="results/merged.bam",
    resources:
        mem_mb=32000,
        cpus_per_task=4,
    shell:
        """
        samtools merge -f {output.merged_bam} {input}
        """
