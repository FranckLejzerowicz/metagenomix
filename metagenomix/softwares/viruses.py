# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
from os.path import basename
from metagenomix._io_utils import io_update, to_do, status_update


def viralverify_cmd(
        self,
        contigs: str,
        out: str
) -> str:
    """Collect ViralVerify command.

    Parameters
    ----------
    self : Commands class instance
        .databases : dict
            Path to the reference databases
        .outputs : dict
            All outputs
    contigs : str
        Path to the input contigs fasta file
    out : str
        Path to the output folder

    Returns
    -------
    cmd : str
        ViralVerify command
    """
    cmd, cmd_rm = '', ''
    if self.soft.params['path']:
        cmd += 'export PATH=$PATH:%s\n' % self.soft.params['path']

    if contigs.endswith('.fa.gz') or contigs.endswith('.fasta.gz'):
        cmd += 'gunzip -c %s > %s\n' % (contigs, contigs.rstrip('.gz'))
        cmd_rm += 'rm %s\n' % contigs.rstrip('.gz')

    cmd += 'viralverify'
    cmd += ' -f %s' % contigs.rstrip('.gz')
    cmd += ' -o %s' % out
    if self.soft.params['db']:
        cmd += ' --db'
    if self.soft.params['p']:
        cmd += ' --p'
    cmd += ' -t %s' % self.soft.params['cpus']
    cmd += ' --thr %s' % self.soft.params['thr']
    cmd += ' --hmm %s\n' % '%s/Pfam-A.hmm' % self.databases.paths.get('pfam')
    cmd += 'for i in %s/*; do if [ ! -d $i ]; then gzip -q $i; fi; done\n' % out
    cmd += 'gzip %s/Prediction_results_fasta/*fasta\n' % out
    cmd += cmd_rm
    return cmd


def viralverify(self) -> None:
    """viralVerify classifies contigs (output of metaviralSPAdes or other
    assemblers) as viral, non-viral or uncertain, based on gene content. Also
    for non-viral contigs it can optionally provide plasmid/non-plasmid
    classification.

    viralVerify predicts genes in the contig using Prodigal in the
    metagenomic mode, runs hmmsearch on the predicted proteins and classifies
    the contig as vrial or non-viral by applying the Naive Bayes classifier
    (NBC). For the set of predicted HMMs, viralVerify uses trained NBC to
    classify this set to be viral or chromosomal.

    To improve results in the case of metagenomes with possible host
    contamination, we recommend users to filter out reads that align to the
    host genome prior to assembly. Since viralVerify is based on gene
    classification, it can be used on contigs on any length, and short
    viruses can be detected as long as they contain a recognizable
    virus-specific gene. To help analyze the rapidly growing amount of novel
    data, we have added a script that allows users to construct their own
    training database from a set of viral, chromosomal and plasmid contigs,
    as well as custom HMM database.

    References
    ----------
    Antipov, Dmitry, et al. "Metaviral SPAdes: assembly of viruses from
    metagenomic data." Bioinformatics 36.14 (2020): 4126-4129.

    Notes
    -----
    GitHub  : https://github.com/ablab/viralVerify
    Paper   : https://doi.org/10.1093/bioinformatics/btaa490

    Parameters
    ----------
    self : Commands class instance
        .soft.prev : str
            Previous software in the pipeline
        .dir : str
            Path to pipeline output folder for CoCoNet
        .soft.sam_pool : str
            Sample or co-assembly group name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    if self.soft.prev not in self.config.tools['assembling']:
        sys.exit('[viralVerify] can only be run after assembly (on contigs)')

    pfam = '%s/Pfam-A.hmm' % self.databases.paths.get('pfam')
    if not self.config.dev and to_do(pfam):
        sys.exit('[viralVerify] Needs the Pfam .hmm database in database yaml')

    for (tech, group), inputs in self.inputs[self.sam_pool].items():

        contigs = inputs[0]
        to_dos = status_update(self, tech, [contigs], group=group)
        out = '/'.join([self.dir, tech, self.sam_pool, group])
        self.outputs['dirs'].append(out)

        base = basename(contigs).rsplit('.', 2)[0]
        out_fp = '%s/%s_result_table.csv.gz' % (out, base)
        self.outputs['outs'][(tech, group)] = out_fp

        if self.config.force or to_do(out_fp):
            key = (tech, group)
            cmd = viralverify_cmd(self, contigs, out)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=contigs, o_d=out, key=key)
            self.soft.add_status(tech, self.sam_pool, 1, group=group)
        else:
            self.soft.add_status(tech, self.sam_pool, 0, group=group)


def coconet_cmd(
        self,
        fp: str,
        bam: list,
        out_dir: str
) -> str:
    """Collect CoCoNet command.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    fp : str
        Path to the contigs file
    bam : list
        Path(s) to spades read mapping BAM files
    out_dir : str
        Path to the output folder

    Returns
    -------
    cmd : str
        CoCoNet command
    """
    cmd = 'coconet run'
    cmd += ' --fasta %s' % fp
    cmd += ' --threads %s' % self.soft.params['cpus']
    for boolean in ['quiet', 'no_rc', 'silent', 'continue',
                    'patience', 'recruit_small_contigs']:
        if self.soft.params[boolean]:
            cmd += ' --%s' % boolean
    for p in [
        'debug', 'flag', 'min_ctg_len', 'min_prevalence', 'min_dtr_size',
        'min_mapping_quality', 'min_aln_coverage', 'fragment_step',
        'n_train', 'n_test', 'batch_size', 'test_batch', 'load_batch',
        'cover_filters', 'cover_kernel', 'cover_stride', 'merge_neurons',
        'kmer', 'wsize', 'wstep', 'n_frags', 'max_neighbors', 'n_clusters',
        'vote_threshold', 'gamma1', 'gamma2', 'fragment_length', 'theta',
        'test_ratio', 'learning_rate'
    ]:
        cmd += ' --%s %s' % (p.replace('_', '-'), self.soft.params[p])
    for p in ['tlen_range', 'compo_neurons', 'cover_neurons']:
        if self.soft.params[p]:
            cmd += ' --%s %s' % (p.replace('_', '-'),
                                 ' '.join(map(str, self.soft.params[p])))
    if bam:
        for b in bam:
            print(kjsrbf)
        cmd += ' --bam %s' % ' '.join(bam)
    cmd += ' --output %s\n' % out_dir
    return cmd


def coconet(self) -> None:
    """CoCoNet (Composition and Coverage Network) is a binning method for viral
    metagenomes. It leverages deep learning to abstract the modeling of the
    k-mer composition and the coverage for binning contigs assembled form
    viral metagenomic data. Specifically, our method uses a neural network to
    learn from the metagenomic data a flexible function for predicting the
    probability that any pair of contigs originated from the same genome.
    These probabilities are subsequently combined to infer bins, or clusters
    representing the species present in the sequenced samples. Our approach
    was specifically designed for diverse viral metagenomes, such as those
    found in environmental samples (e.g., oceans, soil, etc.).

    References
    ----------
    Arisdakessian, Cédric G., et al. "CoCoNet: an efficient deep learning
    tool for viral metagenome binning." Bioinformatics 37.18 (2021): 2803-2810.

    Notes
    -----
    GitHub  : https://github.com/Puumanamana/CoCoNet
    Docs    : https://coconet.readthedocs.io/
    Paper   : https://doi.org/10.1093/bioinformatics/btab213

    Parameters
    ----------
    self : Commands class instance
        .soft.prev : str
            Previous software in the pipeline
        .dir : str
            Path to pipeline output folder for CoCoNet
        .soft.sam_pool : str
            Sample or co-assembly group name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    if self.soft.prev not in self.config.tools['assembling']:
        sys.exit('[coconet] can only be run after assembly')

    bams = {}
    if 'mapping_spades' in self.softs:
        bams = self.softs['mapping_spades'].outputs

    for (tech, group), inputs in self.inputs[self.sam_pool].items():

        to_dos = status_update(self, tech, [inputs[0]], group=group)

        bam = bams.get((tech, group), [])

        out_dir = '/'.join([self.dir, tech, self.sam_pool, group])
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'][group] = out_dir

        log_fp = '%s/coconet.log' % out_dir
        if self.config.force or to_do(log_fp):
            contigs = inputs[0]
            cmd = coconet_cmd(self, contigs, bam, out_dir)
            key = (tech, group)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=contigs, o_d=out_dir, key=key)
            self.soft.add_status(tech, self.sam_pool, 1, group=group)
        else:
            self.soft.add_status(tech, self.sam_pool, 0, group=group)


def threecac(self) -> None:
    """3CAC is a three-class classifier designed to classify contigs in mixed
    metagenome assemblies as phages, plasmids, chromosomes, or uncertain.

    References
    ----------
    Pu, Lianrong, and Ron Shamir. "3CAC: improving the classification of
    phages and plasmids from metagenomic assemblies using assembly graphs."
    bioRxiv (2021).

    Notes
    -----
    GitHub  : https://github.com/Shamir-Lab/3CAC
    Paper   : https://doi.org/10.1101/2021.11.05.467408

    Parameters
    ----------
    self : Commands class instance
        .soft.prev : str
            Previous software in the pipeline
        .dir : str
            Path to pipeline output folder for CoCoNet
        .soft.sam_pool : str
            Sample or co-assembly group name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    pass


def deepvirfinder(self):
    """DeepVirFinder predicts viral sequences using deep learning method. The
    method has good prediction accuracy for short viral sequences, so it can
    be used to predict sequences from the metagenomic data.

    References
    ----------
    Ren, J., Song, K., Deng, C., Ahlgren, N.A., Fuhrman, J.A., Li, Y., Xie,
    X., Poplin, R. and Sun, F., 2020. Identifying viruses from metagenomic
    data using deep learning. Quantitative Biology, 8(1), pp.64-77.

    Notes
    -----
    GitHub  : https://github.com/jessieren/DeepVirFinder
    Paper   : https://doi.org/10.1007/s40484-019-0187-4

    Parameters
    ----------
    self : Commands class instance
        .soft.prev : str
            Previous software in the pipeline
        .dir : str
            Path to pipeline output folder for DeepVirFinder
        .soft.sam_pool : str
            Sample or co-assembly group name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    pass


def wish(self):
    """WIsH can identify bacterial hosts from metagenomic data, keeping good
    accuracy even on smaller contigs.

    References
    ----------
    Galiez, C., Siebert, M., Enault, F., Vincent, J. and Söding, J.,
    2017. WIsH: who is the host? Predicting prokaryotic hosts from
    metagenomic phage contigs. Bioinformatics, 33(19), pp.3113-3114.

    Notes
    -----
    GitHub  : https://github.com/soedinglab/WIsH
    Paper   : https://doi.org/10.1093/bioinformatics/btx383

    Parameters
    ----------
    self : Commands class instance
        .soft.prev : str
            Previous software in the pipeline
        .dir : str
            Path to pipeline output folder for WiSH
        .soft.sam_pool : str
            Sample or co-assembly group name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    pass


def virstrain(self):
    """An RNA virus strain-level identification tool for short reads.

    References
    ----------
    Liao, H., Cai, D. and Sun, Y., 2022. VirStrain: a strain identification
    tool for RNA viruses. Genome Biology, 23(1), p.38.

    Notes
    -----
    GitHub  : https://github.com/liaoherui/VirStrain
    Paper   : https://doi.org/10.1186/s13059-022-02609-x

    Parameters
    ----------
    self : Commands class instance
        .soft.prev : str
            Previous software in the pipeline
        .dir : str
            Path to pipeline output folder for VirStrain
        .soft.sam_pool : str
            Sample or co-assembly group name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    pass
