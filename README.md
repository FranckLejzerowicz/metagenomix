# metagenomix

## Description

Running the metagenomics data analysis from raw fastqs to read-mapped MAGs 
creation and scrutiny for multiple taxonomic, gene and genome targets.
This also eases importing issues and notably run this on a High-Performance 
Computer (HPC), either [Torque](http://docs.adaptivecomputing.com/torque/4-0-2/help.htm) 
or [Slurm](https://slurm.schedmd.com/documentation.html) scheduler, and it 
heavily relies on a you having installed and configured properly multiple conda
environments (explanations and environment.yml files provided in Wiki). 
This tool wrapping other tools could be done using 
[snakemake](https://snakemake.readthedocs.io/en/stable/) but this is meant to 
only create the scripts to run on your samples, not execute these scripts.

## Installation
```
pip install --upgrade git+https://github.com/FranckLejzerowicz/metagenomix.git
```

*_Note that python and pip should be python3_

### Depencency for jobs preparation

[Xpbs](https://github.com/FranckLejzerowicz/Xpbs) must be installed as it 
allows automatic preparation of HPC scripts from the basic bash scripts
written here. For Xpbs to work, it is necessary that the user provide edit the 
config.txt file of this tool (simply adding the email address for job completion
[as explained here](https://github.com/FranckLejzerowicz/Xpbs#requisite)).

### Depencencies for jobs to run

You must install all the databases and softwares (lists below) that the pipeline
will allow you to run. I recommend using conda environments for each tool and 
then specify the name of the environment for this tool in the run parameters
configuration file. I am working on an installer of conda environments to 
alleviate this time-consuming step, but since some tools may require you to edit
configurations for you system, this will not be fully integrated.  

#### Databases

* [PFAM](http://pfam.xfam.org/): allows searching HMM profiles per keywords from the Pfam-A catalogue,
(see [here](https://doi.org/10.1093/nar/gkp985) or 
[here](https://academic.oup.com/nar/article/38/suppl_1/D211/3112325?searchresult=1#64944278)) 
which you need to download: Pfam-A files "Pfam-A.hmm" and "Pfam-A.hmm.dat.gz"
[here](http://ftp.ebi.ac.uk/pub/databases/Pfam/releases)

* [dbCAN](https://bcb.unl.edu/dbCAN/): allows searching HMM profiles per CAZy subfamily.

* [MAR](https://mmp2.sfb.uit.no/databases/): Allow querying the MarDB and MarRef databases for marine microbes. 

* [GTDB](https://gtdb.ecogenomic.org/)


#### Tools

* [Antismash](https://antismash.secondarymetabolites.org/#!/start) Search genomes for secondary metabolite biosynthetic gene clusters([paper](https://academic.oup.com/nar/article/45/W1/W36/3778252))
* [Anvi'o](https://anvio.org/) ([paper](https://www.nature.com/articles/s41564-020-00834-3) Integrated multi-omics at scale. [paper](https://peerj.com/articles/1319/?testing))
* [atropos](https://github.com/jdidion/atropos) Trimming. ([paper](http://journal.embnet.org/index.php/embnetjournal/article/view/200))
* [CheckM](https://github.com/Ecogenomics/CheckM) Assess the quality of microbial genomes recovered from isolates, single cells, and metagenomes. ([paper](https://genome.cshlp.org/content/25/7/1043))
* [cutadapt](https://cutadapt.readthedocs.io/en/stable/) Quality control. ([paper](http://journal.embnet.org/index.php/embnetjournal/article/view/200))
* [Deeplasmid](https://github.com/wandreopoulos/deeplasmid) Separates plasmids from chromosomal sequences (ML). ([paper](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkab1115/6454267))
* [DeepVirFinder](https://github.com/jessieren/DeepVirFinder) Identifying viruses from metagenomic data by (deep learning). ([paper](https://link.springer.com/article/10.1007/s40484-019-0187-4))
* [DIAMOND](https://github.com/bbuchfink/diamond) Accelerated BLAST compatible local sequence aligner. ([paper](https://www.nature.com/articles/s41592-021-01101-x))
* [dRep](https://github.com/MrOlm/drep) Rapid comparison and dereplication of genomes. ([paper](https://www.nature.com/articles/ismej2017126))
* [FastQC](https://github.com/s-andrews/FastQC) A quality control analysis tool for high throughput sequencing data.
* [FLASh](http://www.cbcb.umd.edu/software/flash) Overlapping reads pairs merging. ([paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3198573/))
* [graphlan](https://github.com/biobakery/graphlan) High-quality circular representations of taxonomic and phylogenetic trees.
* [GRASP2](https://sourceforge.net/projects/grasp2/) gene-centric homolog search tool ([paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2818-1))
* [GroopM](https://github.com/Ecogenomics/GroopM) Metagenomic binning suite. ([paper](https://peerj.com/articles/603/))
* [GTDB-tk](https://ecogenomics.github.io/GTDBTk/index.html) Taxonomic classification of bacterial and archaeal genomes based on [GTDB](https://gtdb.ecogenomic.org/). ([paper](https://academic.oup.com/bioinformatics/article/36/6/1925/5626182))
* [HUMAnN 3](https://huttenhower.sph.harvard.edu/humann/) Abundance profiling of microbial taxa and metabolic pathways. ([paper](https://elifesciences.org/articles/65088))
* [I-TASSER](https://zhanggroup.org/I-TASSER/) Protein structure prediction and structure-based function annotation. ([paper](https://www.nature.com/articles/nmeth.3213))
* [Integron_Finder](https://github.com/gem-pasteur/Integron_Finder) Find integrons in bacterial genomes. ([paper](https://academic.oup.com/nar/article/44/10/4539/2516972?login=false))
* [inStrain](https://github.com/MrOlm/inStrain) Stain-level analyses for co-occurring genome populations. ([paper](https://www.nature.com/articles/s41587-020-00797-0))
* [IonCom](https://zhanggroup.org/IonCom/) Ligand-specific method for small ligand (including metal and acid radical ions) binding site prediction. ([paper](https://academic.oup.com/bioinformatics/article/32/21/3260/2415108?login=true))
* [Kneadata](https://github.com/biobakery/kneaddata) Quality control.
* [Kraken2](https://github.com/DerrickWood/kraken2) ([paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0))
* [MacSyFinder](https://github.com/gem-pasteur/macsyfinder) Detection of macromolecular systems in protein datasets using systems modelling and similarity search. ([paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0110726))
* [mapDamage2](https://ginolhac.github.io/mapDamage/) Tracking and quantifying damage patterns in ancient DNA sequences. ([paper](https://academic.oup.com/bioinformatics/article/29/13/1682/184965?login=true))
* [metaclade2](http://gitlab.lcqb.upmc.fr/vicedomini/metaclade2) Multi-source domain annotation. ([paper](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0532-2))
* [Metamarker](https://bitbucket.org/mkoohim/metamarker) De novo pipeline to discover novel metagenomic biomarkers. ([paper](https://academic.oup.com/bioinformatics/article/35/19/3812/5368527))
* [metaphlan](https://huttenhower.sph.harvard.edu/metaphlan/) Taxon abundance profiler based on clade-specific marker genes from ~17,000 reference genomes. ([paper](https://elifesciences.org/articles/65088))
* [MetaWRAP](https://github.com/bxlab/metaWRAP) Flexible pipeline for genome-resolved metagenomic data analysis. ([paper](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0541-1))
* [Metaxa2](https://microbiology.se/software/metaxa2/) Identification and taxonomic classification of SSU and LSU rRNA in metagenomes. ([paper](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12399))
* [MIDAS](https://github.com/snayfach/MIDAS) Strain-level genomic variation ([paper](https://genome.cshlp.org/content/26/11/1612.short))
* [MOCAT2](https://mocat.embl.de/) Taxonomic and functional abundance profiling, reads assembler and gene prediction. ([paper](https://academic.oup.com/bioinformatics/article/32/16/2520/1743334))
* [MultiQC](https://github.com/ewels/MultiQC) Quality control. ([paper](https://academic.oup.com/bioinformatics/article/32/19/3047/2196507))
* [panphlan](https://github.com/segatalab/panphlan) Strain-level profiler for gene composition of individual strains. ([paper](https://elifesciences.org/articles/65088))
* [PhyloFlash](http://hrgv.github.io/phyloFlash/) Reconstruct the SSU rRNAs and explore phylogenetic composition of an illumina (meta)genomic dataset. ([paper](https://journals.asm.org/doi/10.1128/mSystems.00920-20))
* [phylophlan](https://huttenhower.sph.harvard.edu/phylophlan/) ([paper](https://www.nature.com/articles/s41467-020-16366-7))
* [PLASS](https://github.com/soedinglab/plass) Protein-Level ASSembler ([paper](https://www.nature.com/articles/s41592-019-0437-4))
* [Prokka](https://github.com/tseemann/prokka) Rapid prokaryotic genome annotation ([paper](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517))
* [Prodigal](https://github.com/hyattpd/Prodigal) Gene Prediction ([paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-119))
* [QUAST](http://cab.spbu.ru/software/quast/) Genome assembly evaluation tool. ([paper](http://quast.sourceforge.net/))
* [run_dbcan](https://github.com/linnabrown/run_dbcan) Search for CAZymes. ([paper](https://doi.org/10.1093/nar/gkx894))
* [seqtk](https://github.com/lh3/seqtk) Toolkit for processing sequences in FASTA/Q formats. ([paper](https://docs.csc.fi/apps/seqtk/))
* [SHOGUN](https://github.com/knights-lab/SHOGUN) Taxonomic and functional profiling. ([paper](https://academic.oup.com/bioinformatics/article/36/13/4088/5828930?login=true))
* [SignalP](https://services.healthtech.dtu.dk/service.php?SignalP-6.0) Predicts presence of signal peptides and location of their cleavage sites in proteins. ([paper](https://link.springer.com/article/10.1007/s10930-019-09838-3))
* [Simka](https://gatb.inria.fr/software/simka/) De novo metagenomes comparison using k-mer spectra and classical ecological distances. ([paper](https://peerj.com/articles/cs-94/))
* [SPAdes](http://cab.spbu.ru/files/release3.15.3/manual.html) Assembler. ([paper](https://currentprotocols.onlinelibrary.wiley.com/doi/abs/10.1002/cpbi.102))
* [SqueezeMeta](https://github.com/jtamames/SqueezeMeta) Complete pipeline for metagenomic analysis. ([paper](https://www.frontiersin.org/articles/10.3389/fmicb.2018.03349/full))
* [strainphlan](http://segatalab.cibio.unitn.it/tools/strainphlan/) Strain-level population genomics tool. ([paper](https://elifesciences.org/articles/65088))
* [strainsifter](https://github.com/bhattlab/StrainSifter) Detection of a bacterial strain in one or more metagenome(s). ([paper](https://www.nature.com/articles/s41591-018-0202-8?sf200127381=1))
* [StrainPro](https://github.com/hsinnan75/StrainPro#download) Strain-level profiling. ([paper](https://www.biorxiv.org/content/10.1101/807149v1.article-metrics))
* [TMHMM](https://services.healthtech.dtu.dk/service.php?TMHMM-2.0) Prediction of transmembrane helices in proteins. ([paper](https://www.sciencedirect.com/science/article/abs/pii/S0022283600943158?via%3Dihub))
* [vamb](https://github.com/RasmussenLab/) Variational autoencoder for metagenomic binning. ([paper](https://www.nature.com/articles/s41587-020-00777-4))
* [VirStrain](https://github.com/liaoherui/VirStrain) RNA virus strain-level identification for short reads. ([paper](https://www.biorxiv.org/content/10.1101/2020.12.21.423722v2.abstract))
* [Woltka](https://github.com/qiyunzhu/woltka) Taxonomic and functional profiling. [preprint](https://www.biorxiv.org/content/10.1101/2021.04.04.438427v1.abstract)
* [WIsH](https://github.com/soedinglab/WIsH) Predict prokaryotic hosts from metagenomic phage contigs ([paper](https://academic.oup.com/bioinformatics/article/33/19/3113/3964377))
* [YAMB](https://github.com/laxeye/YAMB) (Yet Another Metagenome Binner) - semi-automatic pipeline for metagenomic contigs binning. ([paper](https://www.biorxiv.org/content/10.1101/521286v1))
* Will keep growing - that's the point of this pipeline

## Input

IN CONSTRUCTION

## Outputs

IN CONSTRUCTION
