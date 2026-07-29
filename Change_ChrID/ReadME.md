Genome sequences downloaded from NCBI may feature chromosome identifiers different from the naming convention used in our submitted assembly. Two steps are required to convert chromosome IDs and strand orientations to reproduce the results presented in this study.

Step 1: After downloading the genome assembly, check the chromosome identifiers (formatted like CM148818.1). We use the H. limbatus assembly ASM5614945v1 as an example to demonstrate the conversion procedure. First, refer to the Chromosomes table on the NCBI genome page, which contains columns for Chromosome and GenBank accession. Create a mapping file named hli_rn formatted as shown below:

CM148818.1	chr01	for
CM148819.1	chr02	for
CM148820.1	chr03	for
CM148821.1	chr04	for
CM148822.1	chr05	for
CM148823.1	chr06	for
CM148824.1	chr07	for
CM148825.1	chr08	for
CM148826.1	chr09	for
CM148827.1	chr10	for
CM148828.1	chr11	for
CM148829.1	chr12	for
CM148830.1	chr13	for
CM148831.1	chr14	for
CM148832.1	chr15	for
CM148833.1	chr16	for
CM148834.1	chr17	for
CM148835.1	chr18	for
CM148836.1	chr19	for
CM148837.1	chr20	for

Next, run the provided script:
./rename_genome_fa.sh ncbi_hli.fa hli_rn hli_temp.fa

This command generates the original H. limbatus genome assembly submitted to NCBI. To facilitate interspecific synteny comparisons and result reproducibility, chromosome IDs must be further converted according to homologous relationships with H. sajori.

Step 2: Execute the same script with the supplied hli_rn_file:
./rename_genome_fa.sh hli_temp.fa hli_rn_file hli.fa

The output file hli.fa contains genome sequences with chromosome IDs consistent with those used for all analyses in this manuscript. Identical conversion workflows apply to the remaining four assemblies. The file prefixes correspond to species/assemblies as follows:
hge: Northern Hyporhamphus gernaerti
hli: Hyporhamphus limbatus
shge: Southern Hyporhamphus gernaerti
xhge: X haplotype of male Northern H. gernaerti
yhge: Y haplotype of male Northern H. gernaerti
