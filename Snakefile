rule build:
    input:
        "auspice/ncov_21K-diversity.json",
        "auspice/ncov_21K-diversity_unmasked.json"

rule build_all:
    input:
        "auspice/ncov_21K-all-diversity.json",

rule deploy:
    input:
        "deploy/21K/latest.json",
        "deploy/21K/latest_unmasked.json",

rule deploy_all:
    input:
        "deploy/21K-all/latest.json",

def input_for_do_stuff(wildcards):
    try:
        return next(Path("data").glob("*.tar"))
    except StopIteration as err:
        raise Exception("No tar found") from err

rule unpack_gisaid_tar:
    input: input_for_do_stuff 
    output:
        sequences = "builds/gisaid.fasta",
        metadata = "builds/metadata_omicron.tsv",
    shell:
        """
        tar -xvf {input} -C data/unpacked
        tsv-uniq -H -f strain data/unpacked/*metadata.tsv >{output.metadata}
        seqkit rmdup -i data/unpacked/*.fasta >{output.sequences}
        """

rule download_nextclade_dataset:
    output: directory("builds/nextclade_dataset")
    shell: "nextclade dataset get --name='sars-cov-2' --output-dir={output}"

rule download_lat_longs:
    output: "builds/lat_longs.tsv"
    params:
        url = "https://raw.githubusercontent.com/nextstrain/ncov/master/defaults/lat_longs.tsv"
    shell: 
        """
        curl {params.url} | \
        sed "s/North Rhine Westphalia/North Rhine-Westphalia/g" | \
        sed "s/Baden-Wuerttemberg/Baden-Wurttemberg/g" \
        > {output}
        """

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

rule subsample_fasta:
    input: "builds/gisaid.fasta"
    output: "builds/subsampled.fasta"
    shell: "seqkit sample -n 1500 {input} >{output}"

rule exclude_outliers:
    input:
        sequences = "builds/subsampled.fasta",
        metadata = "builds/metadata_omicron.tsv",
        exclude = "data/{build}/exclude.txt",
    output:
        sampled_sequences = "builds/{build}/filtered.fasta",
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
        omicron = "builds/{build}/filtered.fasta",
    output:
        "builds/{build}/omicron.fasta",
    shell:
        "cat {input.reference} {input.omicron} > {output}"

rule join_ref_meta:
    input:
        reference = "data/root_meta.tsv",
        omicron = "builds/metadata_omicron.tsv",
    output:
        "builds/metadata_raw.tsv",
    shell:
        "cat {input.omicron} {input.reference} > {output}"

rule remove_false_meta_linebreaks:
    input:
        "builds/metadata_raw.tsv",
    output:
        "builds/metadata.tsv",
    run:
        import re
        with open(input[0], "r") as f:
            string_in = f.read()
        string_out = re.sub(r"\n(?!(hCoV-19/|MN908947))"," ",string_in)
        with open(output[0], "w") as f: 
            f.write(str(string_out))

rule nextclade:
    input: 
        fasta = "builds/{build}/omicron.fasta",
        dataset = rules.download_nextclade_dataset.output
    output:
        alignment = "pre-processed/{build}/nextclade/premask.aligned.fasta",
        translations = expand("pre-processed/{{build}}/nextclade/premask.gene.{gene}.fasta", gene=genes),
    shell:
        """
        nextalign \
            --sequences {input.fasta} \
            --reference "builds/nextclade_dataset/reference.fasta" \
            --genemap "builds/nextclade_dataset/genemap.gff" \
            --genes E,M,N,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S \
            --output-dir pre-processed/{wildcards.build}/nextclade \
            --output-basename premask \
            --jobs 0
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
        alignment = "pre-processed/{build}/nextclade/premask.aligned.fasta",
    output:
        alignment = "builds/{build}/masked.fasta",
    shell:
        """
        python3 scripts/mask-alignment.py \
            --alignment {input.alignment} \
            --mask-from-beginning 100 \
            --mask-from-end 200 \
            --mask-terminal-gaps \
            --mask-site-file data/exclude_sites_light.txt \
            --output {output.alignment}
        """

rule nextclade_after_mask:
    input: 
        fasta = "builds/{build}/masked.fasta",
        dataset = rules.download_nextclade_dataset.output
    output:
        alignment = "pre-processed/{build}/nextclade/omicron.aligned.fasta",
        translations = expand("pre-processed/{{build}}/nextclade/omicron.gene.{gene}.fasta", gene=genes),
    shell:
        """
        nextalign \
            --sequences {input.fasta} \
            --reference "builds/nextclade_dataset/reference.fasta" \
            --genemap "builds/nextclade_dataset/genemap.gff" \
            --genes E,M,N,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S \
            --output-dir pre-processed/{wildcards.build}/nextclade \
            --output-basename omicron \
            --jobs 0
        """

rule tree:
    input:
        alignment = rules.nextclade_after_mask.output.alignment,
    output:
        tree = "builds/{build}/tree_raw.nwk"
    params:
        args = "'-ninit 10 -n 4 -czb -T AUTO -ntmax 10'",
        # exclude_sites = "data/exclude_sites.txt"
        exclude_sites = "data/exclude_none.txt"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --tree-builder-args {params.args} \
            --output {output.tree} \
            --exclude-sites {params.exclude_sites}
        """

rule refine:
    input:
        tree = rules.tree.output.tree,
        alignment = rules.mask.output.alignment,
        metadata = "builds/metadata.tsv",
    output:
        tree = "builds/{build}/tree.nwk",
        node_data = "builds/{build}/branch_lengths.json"
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
            # --timetree \
            --keep-polytomies \
            --clock-rate 0.0006 \
            --clock-std-dev 0.003 \
            --coalescent skyline \
            --date-inference marginal \
            --date-confidence \
            --no-covariance
        """

rule ancestral:
    input:
        tree = rules.refine.output.tree,
        alignment = rules.mask.output.alignment
    output:
        node_data = "builds/{build}/nt_muts.json"
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

rule ancestral_unmasked:
    input:
        tree = rules.refine.output.tree,
        alignment = rules.mask.input.alignment
    output:
        node_data = "builds/{build}/nt_muts_unmasked.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference} \
            --keep-ambiguous \
            --keep-overhangs \
            2>&1 | tee {log}
        """


rule translate:
    input:
        tree = rules.refine.output.tree,
        translations = expand("pre-processed/{{build}}/nextclade/omicron.gene.{gene}.fasta", gene=genes),
        reference = "data/reference_seq.fasta",
        genemap = "data/annotation.gff",
    output:
        node_data = "builds/{build}/aa_muts.json",
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

rule translate_unmasked:
    input:
        tree = rules.refine.output.tree,
        translations = expand("pre-processed/{{build}}/nextclade/premask.gene.{gene}.fasta", gene=genes),
        reference = "data/reference_seq.fasta",
        genemap = "data/annotation.gff",
    output:
        node_data = "builds/{build}/aa_muts_unmasked.json",
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

rule recency:
    input:
        metadata= "builds/metadata.tsv",
    output:
        node_data = "builds/{build}/recency.json",
    shell:
        """
        python3 scripts/construct_recency.py \
            --metadata {input.metadata} \
            --output {output} 2>&1 | tee {log}
        """

def _get_node_data_by_wildcards(wildcards):
    """Return a list of node data files to include for a given build's wildcards.
    """
    # Define inputs shared by all builds.
    wildcards_dict = dict(wildcards)
    inputs = [
        rules.refine.output.node_data,
        rules.recency.output.node_data,
    ]
    inputs = [input_file.format(**wildcards_dict) for input_file in inputs]
    return inputs

rule export:
    input:
        tree = rules.refine.output.tree,
        node_data = _get_node_data_by_wildcards,
        ancestral = "builds/{build}/nt_muts{masking}",
        translate = "builds/{build}/aa_muts{masking}",
        auspice_config = "data/auspice_config.json",
        lat_longs = rules.download_lat_longs.output,
        metadata = "builds/metadata.tsv",
    output:
        auspice_json = "auspice/ncov_{build}-diversity{masking}",
    shell:
        """
        export AUGUR_RECURSION_LIMIT=10000;
        augur export v2 \
            --tree {input.tree} \
            --node-data {input.node_data} {input.ancestral} {input.translate} \
            --include-root-sequence \
            --auspice-config {input.auspice_config} \
            --lat-longs {input.lat_longs} \
            --output {output.auspice_json} \
            --metadata {input.metadata} \
            --description data/description.md \
            # --title "Phylogenetic analysis of the ORF1ab genes of the Omicron variant"
        """

rule deploy_single:
    input: "auspice/ncov_{build}-diversity{masking}",
    output: 'deploy/{build}/latest{masking}',
    shell: 
        """
        nextstrain deploy s3://nextstrain-neherlab {input} 2>&1 && touch {output}
        """