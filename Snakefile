def input_for_do_stuff(wildcards):
    try:
        return next(Path("data").glob("*.tar"))
    except StopIteration as err:
        raise Exception("No tar found") from err

rule unpack_gisaid_tar:
    input: input_for_do_stuff 
    output:
        sequences = "data/gisaid.fasta",
        metadata = "data/metadata_omicron.tsv",
    shell:
        """
        tar -xvf {input} -C data/unpacked
        mv data/unpacked/*metadata.tsv {output.metadata}
        mv data/unpacked/*.fasta {output.sequences}
        """

rule download_nextclade_dataset:
    output: directory("data/nextclade_dataset")
    shell: "nextclade dataset get --name='sars-cov-2' --output-dir={output}"

genes = [
    "E",
    "M",
    "N",
    "ORF1a",
    "ORF1b",
    "ORF3a",
    "ORF6",
    "ORF7a",
    "ORF7b",
    "ORF8",
    "ORF9b",
    "S",
]

rule exclude_outliers:
    input:
        sequences = "data/gisaid.fasta",
        metadata = "data/metadata_omicron.tsv",
        exclude = "data/exclude.txt",
    output:
        sampled_sequences = "data/filtered.fasta",
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --output {output.sampled_sequences}
        """

rule join_ref_fasta:
    input:
        reference = "data/reference_seq.fasta",
        omicron = "data/filtered.fasta",
    output:
        "data/omicron.fasta",
    shell:
        "cat {input.reference} {input.omicron} > {output}"

rule join_ref_meta:
    input:
        reference = "data/root_meta.tsv",
        omicron = "data/metadata_omicron.tsv",
    output:
        "data/metadata.tsv",
    shell:
        "cat {input.omicron} {input.reference} > {output}"

rule nextclade:
    input: 
        fasta = "data/omicron.fasta",
        dataset = rules.download_nextclade_dataset.output
    output:
        alignment = "pre-processed/nextclade/premask.aligned.fasta",
    shell:
        """
        nextclade run \
            --input-fasta {input.fasta} \
            --input-dataset {input.dataset} \
            --output-dir pre-processed/nextclade \
            --jobs 0 \
            -n premask
        """

# no mask
# rule mask:
#     input:
#         alignment = "pre-processed/nextclade/omicron.aligned.fasta",
#     output:
#         alignment = "builds/masked.fasta",
#     shell: "mv {input.alignment} {output.alignment}"


rule mask:
    input:
        alignment = "pre-processed/nextclade/premask.aligned.fasta",
    output:
        alignment = "builds/masked.fasta",
    shell:
        """
        python3 scripts/mask-alignment.py \
            --alignment {input.alignment} \
            --mask-from-beginning 300 \
            --mask-from-end 300 \
            --mask-site-file data/exclude_sites_light.txt \
            --output {output.alignment}
        """

rule nextclade_after_mask:
    input: 
        fasta = "builds/masked.fasta",
        dataset = rules.download_nextclade_dataset.output
    output:
        alignment = "pre-processed/nextclade/omicron.aligned.fasta",
        translations = expand("pre-processed/nextclade/omicron.gene.{gene}.fasta", gene=genes),
    shell:
        """
        nextclade run \
            --input-fasta {input.fasta} \
            --input-dataset {input.dataset} \
            --output-dir pre-processed/nextclade \
            --jobs 0 \
            -n omicron
        """

rule tree:
    input:
        alignment = rules.nextclade_after_mask.output.alignment,
    output:
        tree = "builds/tree_raw.nwk"
    params:
        args = "'-ninit 10 -n 4 -czb'",
        # exclude_sites = "data/exclude_sites.txt"
        exclude_sites = "data/exclude_none.txt"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --tree-builder-args {params.args} \
            --output {output.tree} \
            --exclude-sites {params.exclude_sites} \
            --nthreads 8
        """

rule refine:
    input:
        tree = rules.tree.output.tree,
        alignment = rules.mask.output.alignment,
        metadata = "data/metadata.tsv",
    output:
        tree = "builds/tree.nwk",
        node_data = "builds/branch_lengths.json"
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --root MN908947 \
            --alignment {input.alignment} \
            --divergence-units mutations \
            --output-tree {output.tree} \
            --metadata {input.metadata} \
            --output-node-data {output.node_data} \
        """

rule ancestral:
    input:
        tree = rules.refine.output.tree,
        alignment = rules.mask.output.alignment
    output:
        node_data = "builds/nt_muts.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference} \
            --infer-ambiguous 2>&1 | tee {log}
        """


rule translate:
    input:
        tree = rules.refine.output.tree,
        translations = expand("pre-processed/nextclade/omicron.gene.{gene}.fasta", gene=genes),
        reference = "data/reference_seq.fasta",
        genemap = "data/annotation.gff",
    output:
        node_data = "builds/aa_muts.json",
    params:
        genes = genes
    shell:
        """
        python3 scripts/explicit_translation.py \
            --tree {input.tree} \
            --annotation {input.genemap} \
            --reference {input.reference} \
            --translations {input.translations:q} \
            --genes {params.genes} \
            --output {output.node_data} 2>&1 | tee {log}
        """

def _get_node_data_by_wildcards(wildcards):
    """Return a list of node data files to include for a given build's wildcards.
    """
    # Define inputs shared by all builds.
    wildcards_dict = dict(wildcards)
    inputs = [
        rules.refine.output.node_data,
        rules.ancestral.output.node_data,
        rules.translate.output.node_data,
    ]
    inputs = [input_file.format(**wildcards_dict) for input_file in inputs]
    return inputs

rule export:
    input:
        tree = rules.refine.output.tree,
        node_data = _get_node_data_by_wildcards,
        auspice_config = "data/auspice_config.json",
        metadata = "data/metadata.tsv",
    output:
        auspice_json = "auspice/auspice.json",
    shell:
        """
        export AUGUR_RECURSION_LIMIT=10000;
        augur export v2 \
            --tree {input.tree} \
            --node-data {input.node_data} \
            --include-root-sequence \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json} \
            --metadata {input.metadata} \
            --description data/description.md \
            # --title "Phylogenetic analysis of the ORF1ab genes of the Omicron variant"
        """
