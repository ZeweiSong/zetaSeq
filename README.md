<img src="zetaSeq.png" alt="zetaSeq" width="200"/>

`zetaSeq` is a Python library and program collections I wrote for myself when dealing with amplion data and later shotgun data of microbial communities. It is not a full pipeline, but contains various tools I found useful. Most of them were created because I cannot find the same function on the Internet. The order of sections however is organized in a workflow way of a metagenomic analysis.

[1. Installation](https://gitlab.genomics.cn/meer/zetaSeq#installation)

[2. Basic sequence input and output](https://gitlab.genomics.cn/meer/zetaSeq#basic-sequence-input-and-output)

[3. Qulity control](https://gitlab.genomics.cn/meer/zetaSeq#qulity-control)

[4. Sample preprocess](https://gitlab.genomics.cn/meer/zetaSeq#sample-preprocess)

[5. Alignment](https://gitlab.genomics.cn/meer/zetaSeq#alignment)

[6. A Bayesian profiler](https://gitlab.genomics.cn/meer/zetaSeq#a-bayesian-profiler)

[7. Biom format](https://gitlab.genomics.cn/meer/zetaSeq#biom-format)

[8. Contigs](https://gitlab.genomics.cn/meer/zetaSeq#contigs)

[9. Ribosomal RNA](https://gitlab.genomics.cn/meer/zetaSeq#ribosomal-rna)

[10. NCBI Taxdump](https://gitlab.genomics.cn/meer/zetaSeq#ncbi-taxdump)

[11. DereCo](https://gitlab.genomics.cn/meer/zetaSeq#dereco)

# Installation

The only thing you need is a Python 3 environment. Anaconda is a good choice (https://www.anaconda.com/products/distribution).

    python –-version

Clone or download zetaSeq repository to your path.

    git clone https://gitlab.genomics.cn/meer/zetaSeq
    cd zetaSeq
    chmod +x biom_concat contister contister_gather keep_pair_aln lnster rrnaster solve_profile truncaster
    export PATH="your/path/to/zetaSeq:$PATH"
    export PYTHONPATH="your/path/to/zetaSeq:$PYTHONPATH"

# Basic sequence input and output
`zetaSeq.io` is a module for handing FASTA/Q data. Sequences are returned as a list of tuple (label, sequence) for FASTA, or (label, sequence, +, quality) for FASTQ.

    #!/usr/bin/env python
    from zetaSeq import io as seqIO
    fasta_file = 'my_fasta.fna'
    seqs = [i for i in seqIO.sequence(fasta_file)]
    print(len(seq))
    print(seqs[0])
    print(seqs[0][0])
    print(seqs[0][1])
    
    seqs[0][0] = 'This is a reverse compliment of me'
    seqs[0][1] = seqIO.revcomp(seqs[0][1])
    count = seqIO.write_seqs(seqs, 'my_output.fasta', fastx='a', gz=False)

`zetaSeq.io.sequence(file_path)` Read FASTA/Q sequences as a list of tuple. FASTA sequences **has to be** in one line.

`zetaSeq.io.sequence_multiline_fasta(file_path)` Read multiline FASTA sequences as a list of tuple. You will need this one when reading downloaded genomes, which usually wraps sequences by 80 bases.

`zetaSeq.io.sequence_twin(file1_path, file2_path)` Read paired FASTA/Q sequences, return two list of tuples. This is useful when dealing with raw paired-end data.

`zetaSeq.io.sequence_bytes(file_path)` As same as zetaSeq.io.sequence(), but this one read in a trunk of sequences by bytes, which is 40% faster. But usually reading sequences is not the bottleneck of our analysis.
2

# Qulity control

`truncaster` Trim sequences or pair of sequences until they reach the **expect error rate** threshold (default = 0.01). Sequences shorter than minimal length (default = 100) will be discarded. For pair sequences, both in one pair has to pass the length threshold.

**Expected error** is the sum of all P in a given sequence, in which P = 10^(-Q/10). It can be seemed as the expected number of bases to be wrong in this sequence.
**Expected error rate** is expected error per 100 bp.

    truncaster -r1 my_sequences.1.fq.gz -r2 my_sequeces.2.fq.gz -o my_qced_sequences -r 0.01 -l 100 # trim pair end sequences
    truncaster -r1 my_sequences.fq.gz -o my_qced_sequences.fq.gz # Trim single end sequences, -r and -l has default value.
    truncaster -r1 my_sequences.1.fq.gz -r2 my_sequeces.2.fq.gz -o my_qced_seque -t 1000000 # set a larger trunk

# Sample preprocess

`lnster` Create softlinks and rename raw FASTQ files.

    lnster path_list.tsv

You've got your sample sequenced, but the raw data is scatteered in the path of the data releaser. Before launch your awsome pipelines, we would like to softlink all files into a single folder, and renamed them with meaning names. `lnster` need a tab-delimited file as input. Put the/path/to/your/raw/data in the fist column, and sample name in the second column, for example:

    /path/to/the/releaser/data/lanexxx1  SampleA
    /path/to/the/releaser/data/lanexxx3  SampleB
    /path/to/the/releaser/data/lanexxx5  SampleC

`lnster` will check each path and look for FASTQ files. All softlinks will be put under folder data/samples/. So for the tsv file above, you will get:

    ls -l data/samples/
    SampleA.1.fq.gz
    SampleA.2.fq.gz
    SampleB.1.fq.gz
    SampleB.2.fq.gz
    SampleC.1.fq.gz
    SampleC.2.fq.gz

This is my favorite way to prepare working directory for my **snakemake** pipelines.

# Alignment

`keep_paired_aln` Given a pair-end sequence sets. If R1 align to target2 with Aln_1, and R2 align to targets with Aln_2, keeps the paired part in Aln_1 and Aln_2 as a single alignment sets.

    keep_paired_aln -a my_aln_1.b6 my_aln_2.b6 -o my_aln.b6

# A Bayesian profiler

`solve_profile` Bayesian Inference Of Multi-Alignments on genome SetS (BIOMASS). Given an alignment and associated taxonomy of the targets, BIOMASS perform a Bayesian inference on the read count of all presenting target. The inferred count can be feeded back for a permutation.

The alignment can be from amplicon data or shotgun data, and has to be in a aligning mode that report all best ties, or all ties above the threshold.

    solve_profile -a my_alignment -t path/to/the/reference/taxonomy.tsv -l species -n my_sample -p 10 -o my_profile.tsv

# Biom format

`biom_concat` Concatenate multiple biom format profiles into a single biom file. You may need to install biom format `pip install biom-format` first.

    biom convert -i my_profile_1.tsv -o my_profile_1.biom --to-json # convert tsv to biom
    biom convert -i my_profile_2.tsv -o my_profile_2.biom --to-json
    biom convert -i my_profile_3.tsv -o my_profile_3.biom --to-json
    biom_concat -i my_profile_1.tsv my_profile_2.tsv my_-profile_3.tsv -o my_profile.biom # Combine the three profiles
    biom convert -i my_profle.biom -o my_profile.tsv --to-tsv # Convert the biom back to tsv

# Contigs

Another path for metagenomics is to assembly randomly sequenced reads into contigs and then genome bins. `contister` gives you a quick brief on the contigs you got, `contister_gather` put all report together for a table view.

    contister -i my_contig.fa -o my_contig.json -m 2000 # `contister` report in `JSON` format. Use `-m` to set a length cutoff.
    contister_gather -i $(ls -d *.json) -o contigs_report.tsv # Gather all JSON report under current directory into a tab-delimited report.

# Ribosomal RNA

Now we got a set of contigs, one straight forward thing is to find all ribosomal RNA sequences within them. Ribosomal RNA (16S, 23S, 18S, 28S, 5S, 5.8S, 12S) is the most popular taxonomic marker so far. Mining them out may give us some hint on the present species.

`rrnaster` Predict ribosomal RNA sequences from a contig. You need to download Rfam covaraince models as references.

    mkdir -p rfam/cm/
    wget http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.tar.gz;mv Rfam.tar.gz rfam/cm/;tar xzvf Rfam.tar.gz
    rrnaster -g my_contig.fa -d /path/to/rfam/cm/ -r all -o rrna_report.tsv -c rrna_concensus.fa -t 4 # rrnaster can search models for all Kingdoms at the same time, predicted region overlapped will be clustered later based on hit E-value.

# NCBI Taxdump

NCBI is nasty for using a string of number known as `taxid` to organism their genomes. While `taxid` solve the issue of redundancy, it makes us difficult to track the taxonomy of NCBI genomes.

`taxdumper` parse scientifc names from names.dmp, taxonomic hierarchy from nodes.dmp, and put them together in a graph. It is not easy for us to search taxid for its scientifc name, or a name for its taxid, as well as doing uptracing and downtracing.

    taxdumper download
    taxdumper search # search is interactive, so jsutt follow the guide.

`search_taxdump_for_assmebly_summary` use the same method, we can add a real taxonomy to all genomes available on NCBI.

    wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt # This is the summary of all genomes on NCBI with taxid only.
    taxdumper.py download
    search_taxdump_for_assmebly_summary
    head ncbi_genbank_genomes.txt # You can find the nicely parsed taxonomy under this file.

# DereCo

Remove redundant (duplicated) sequences from a group of contigs generated from a co-assembly approach (samples assembled invididually, not combined, usually due to large sample number).

    # Preparea working folder
    ./dereco prep -i data/contigs.fa -o derep_out/ -c 1000000

    # We set the length cutoff to 1M,
    # sequences longer than 1000000 will be saved in one file,
    # and sequences shorter than 1M will be saved in another file.
    # First we will dereplicate the large file by aligning it to itself.
    # Then we dereplicate the small file by aligning it to the new large file.
    # Finaly, we dereplicate the small file further by aligning it to itself.
    ./dereco cacl -i derep_out/ -t 8
