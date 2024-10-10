## Background

Extension to the Crop Haplotypes website that "provides an interactive graphic visualisation of the shared haplotypes between the wheat genome assemblies generated as part of the 10+ Wheat Genome Project (Brinton et al, 2020; Walkowiak et al, 2020).", that I produced as part of my Master thesis of the MSc. in Agrobiotechnology by the Justus Liebig University of Giessen.

## Aim

These R scripts enable extracting useful quality properties from the genome assemblies, such as the alignment density or alignment continuity. These can be useful for decision-making on the haplotype information shown on the website, since this only allows a maximum resolution of 1 Mbp and combined chromosome level assemblies, which are highly dense in alignments, with scaffold-level assemblies, that are prone to regions with low alignment density.

## Code

- 'functions_hbpa.r' is required by every script in this repository

- 'template_hbpa.r' can be used to personalize parameters and source the script. All the information will come up on the console and in plots

- 'demonstration_rdma_hbpa.r' uses parameters from the haploblock RDMa, subject topic for the author's master thesis and can serve as an example of use

- 'Haplotype-based Pangenome Analysis in Wheat.pdf' presents an overview from the last script's outcome

- 'table_assembly_chr_length.txt' contains the chromosome length for all chromosomes in the 15 wheat assemblies

- 'alignment_props_supplementary_hbpa.r' was used for the master thesis to generate data on alignment properties from all rds files obtained by Brinton and colleagues and to apply a statistical analysis to prove the effects of using scaffolds as query in pairwise comparisons

## References
Brinton, J., Ramirez-Gonzalez, R. H., Simmonds, J., Wingen, L., Orford, S., Griffiths, S., 10 Wheat Genome Project, Haberer, G., Spannagl, M., Walkowiak, S., Pozniak, C., & Uauy, C. (2020). A haplotype-led approach to increase the precision of wheat breeding. Communications biology, 3(1), 712. https://doi.org/10.1038/s42003-020-01413-2
