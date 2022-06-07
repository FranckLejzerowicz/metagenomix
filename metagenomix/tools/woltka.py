# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources
from os.path import basename, isfile, splitext

RESOURCES = pkg_resources.resource_filename('metagenomix', 'resources')


def get_aligments(inputs: dict, outputs: dict) -> dict:
    """Get the alignment paths per sample.

    Parameters
    ----------
    inputs : dict
        Input files
    outputs : dict
        All outputs

    Returns
    -------
    alis : dict
        Alignments
    """
    alis = {}
    for sam, sam_inputs in inputs.items():
        for ali_tax_db, files in sam_inputs.items():
            if ali_tax_db == ('bowtie2', 'tax', 'wol'):
                alis[sam] = [x for x in files if x.endswith('.sam')][0]
    outputs['io']['I']['f'].update([alis[x] for x in alis])
    return alis


def write_woltka_map(alis: dict, out_dir: str, outputs: dict) -> str:
    """Write the mapping file that servers as input to Woltka.

    Parameters
    ----------
    alis : dict
        Alignments.
    out_dir : str
        Path to the output folder.
    outputs: dict
        All outputs

    Returns
    -------
    woltka_map : str
        Path to the output woltka samples file.
    """
    woltka_map = '%s/samples.map' % out_dir
    for idx, sam in enumerate(alis.keys()):
        echo = 'echo -e "%s\\t%s"' % (sam, alis[sam])
        if idx:
            outputs['cmds'].append('%s >> %s' % (echo, woltka_map))
        else:
            outputs['cmds'].append('%s > %s' % (echo, woltka_map))
    return woltka_map


def get_gotus(out_dir: str, woltka_map: str, outputs: dict) -> str:
    """Get the Woltka command for the gotu classification.

    Parameters
    ----------
    out_dir : str
        Path to pipeline output folder for SHOGUN.
    woltka_map : str
        Path to the output woltka samples file.
    outputs: dict
        All outputs

    Returns
    -------
    genomes_out : str
        Path to the output genomes per sample table.
    """
    genomes_out = '%s/genomes.tsv' % out_dir
    if not isfile(genomes_out):
        cmd = 'woltka gotu'
        cmd += ' -i %s' % woltka_map
        cmd += ' -o %s' % genomes_out
        outputs['cmds'].append(cmd)
        outputs['outs'].append(genomes_out)
        outputs['io']['O']['f'].add(genomes_out)
    return genomes_out


def get_tax_cmd(woltka_map: str, out_dir: str, path: str,
                outputs: dict) -> str:
    """Get the taxonomic classification outputs and prepare the
    Woltka commands for this classification.

    Parameters
    ----------
    woltka_map : str
        Path to the Woltka input file.
    out_dir : str
        Path to pipeline output folder for SHOGUN.
    path : str
        Path to the woltka database.
    outputs: dict
        All outputs

    Returns
    -------
    tax_outmap : str
        Path to the folder containing the taxonomic maps.
    """
    tax_out, tax_outmap = '%s/taxa' % out_dir, '%s/taxmap' % out_dir
    tax_outputs = ['phylum', 'genus', 'species', 'none']
    tax_todo = []
    for tdx, tax_output in enumerate(tax_outputs):
        cur_tax_output = '%s/%s.tsv' % (tax_out, tax_output)
        outputs['cmds'].append(cur_tax_output)
        if not isfile(cur_tax_output):
            tax_todo.append(tax_output)
            outputs['io']['O']['f'].add(cur_tax_output)

    taxid = '%s/taxonomy/taxid.map' % path
    nodes = '%s/taxonomy/nodes.dmp' % path
    names = '%s/taxonomy/names.dmp' % path
    outputs['io']['I']['f'].update([taxid, nodes, names])
    if len(tax_todo):
        cur_cmd = '\n# taxonomic\n'
        cur_cmd += 'woltka classify'
        cur_cmd += ' -i %s' % woltka_map
        cur_cmd += ' --map %s' % taxid
        cur_cmd += ' --nodes %s' % nodes
        cur_cmd += ' --names %s' % names
        cur_cmd += ' --rank %s' % ','.join(tax_outputs)
        cur_cmd += ' --add-rank'
        cur_cmd += ' --add-lineage'
        cur_cmd += ' --name-as-id'
        cur_cmd += ' --to-tsv'
        cur_cmd += ' --outmap %s' % tax_outmap
        cur_cmd += ' -o %s' % tax_out
        outputs['cmds'].append(cur_cmd)
    return tax_outmap


def classif_go(woltka_map: str, out_dir: str, path: str, tax_outmap: str,
               outputs: dict):
    """Get the taxonomic classification outputs and prepare the Woltka
    commands for this classification.

    Parameters
    ----------
    woltka_map : str
        Path to the Woltka input file.
    out_dir : str
        Path to pipeline output folder for SHOGUN.
    path : str
        Path to the woltka database.
    tax_outmap : str
        Path to the folder containing the taxonomic maps.
    outputs: dict
        All outputs
    """
    coords = '%s/release/function/coords.txt.xz' % path
    uniref_map = '%s/release/function/uniref/uniref.map.xz' % path
    uniref_names = '%s/release/function/uniref/uniref.names.xz' % path
    outputs['io']['I']['f'].update([coords, uniref_map, uniref_names])
    go_rt = '%s/release/function/go' % path
    gos = ['process', 'function', 'component']
    woltka_fun_out = '%s/go' % out_dir
    for go in gos:
        cmd = '\n# %s [no stratification]\n' % go
        cmd += 'woltka classify'
        cmd += ' -i %s' % woltka_map
        cmd += ' --coords %s' % coords
        cmd += ' --map-as-rank'
        cmd += ' --rank %s' % go
        cmd += ' --map %s' % uniref_map
        cmd += ' --to-tsv'
        cur_map = '%s/%s.map.xz' % (go_rt, go)
        cmd += ' --map %s' % cur_map
        outputs['io']['I']['f'].add(cur_map)
        cur_out = '%s/%s.tsv' % (woltka_fun_out, go)
        cmd += ' -o %s' % cur_out
        if not isfile(cur_out):
            outputs['cmds'].append(cmd)
            outputs['io']['O']['f'].add(cur_out)
        outputs['outs'].append(cur_out)

    stratifications = ['phylum', 'family', 'genus', 'species']
    for stratification in stratifications:
        woltka_fun_out = '%s/go_%s' % (out_dir, stratification)
        for go in gos:
            cmd = '\n# %s [%s]\n' % (go, stratification)
            cmd += 'woltka classify'
            cmd += ' -i %s' % woltka_map
            cmd += ' --coords %s' % coords
            cmd += ' --map-as-rank'
            cmd += ' --rank %s' % go
            cmd += ' --stratify %s/%s' % (tax_outmap, stratification)
            cmd += ' --map %s' % uniref_map
            cmd += ' --to-tsv'
            cur_map = '%s/%s.map.xz' % (go_rt, go)
            cmd += ' --map %s' % cur_map
            go_out = '%s/%s.tsv' % (woltka_fun_out, go)
            cmd += ' -o %s' % go_out
            if not isfile(go_out):
                outputs['cmds'].append(cmd)
                outputs['io']['O']['f'].add(go_out)
            outputs['outs'].append(go_out)


def classif_genes(woltka_map: str, out_dir: str, path: str, outputs: dict):
    """Get the Woltka commands for the gene-level classification.

    Parameters
    ----------
    woltka_map : str
        Path to the Woltka input file.
    out_dir : str
        Path to pipeline output folder for SHOGUN.
    path : str
        Path to the woltka database.
    outputs : dict
        All outputs

    Returns
    -------
    genes : str
        Path to the genes classification output.
    """
    coords = '%s/release/function/coords.txt.xz' % path
    genes = '%s/wol_genes.biom' % out_dir
    if not isfile(genes):
        cmd = '\n# per gene\n'
        cmd += 'woltka classify'
        cmd += ' -i %s' % woltka_map
        cmd += ' --coords %s' % coords
        cmd += ' -o %s' % genes
        outputs['cmds'].append(cmd)
        outputs['io']['O']['f'].add(genes)
    else:
        outputs['io']['I']['f'].add(genes)
    outputs['outs'].append(genes)
    return genes


def classif_uniref(genes: str, out_dir: str, path: str, outputs: dict):
    """Get the Woltka commands for the uniref-level classification.

    Parameters
    ----------
    genes : str
        Path to the genes classification output.
    out_dir : str
        Path to pipeline output folder for SHOGUN.
    path : str
        Path to the woltka database.
    outputs : dict
        All outputs

    Returns
    -------
    uniref : str
        Path to the uniref classification.
    """
    uniref_map = '%s/release/function/uniref/uniref.map.xz' % path
    uniref_names = '%s/release/function/uniref/uniref.names.xz' % path
    uniref = '%s/wol_uniref.biom' % out_dir
    if not isfile(uniref):
        cmd = 'woltka tools collapse'
        cmd += ' --input %s' % genes
        cmd += ' --map %s' % uniref_map
        cmd += ' --names %s' % uniref_names
        cmd += ' --output %s' % uniref
        outputs['cmds'].append(cmd)
        outputs['io']['O']['f'].add(uniref)
    else:
        outputs['io']['I']['f'].add(uniref)
    outputs['outs'].append(uniref)
    return uniref


def classif_eggnog(uniref: str, out_dir: str, outputs: dict):
    """Get the Woltka commands for the eggnog-level classification.

    Parameters
    ----------
    uniref : str
        Path to the uniref classification.
    out_dir : str
        Path to pipeline output folder for SHOGUN.
    outputs : dict
        All outputs
    """
    biom = '%s/eggnog/wol_eggnog.biom' % out_dir
    if not isfile(biom):
        cmd = 'woltka tools collapse'
        cmd += '--input %s' % uniref
        cmd += ' --map /projects/wol/release/annotation/eggnog/eggnog.map.xz'
        cmd += ' --output %s\n\n' % biom
        outputs['cmds'].append(cmd)
        outputs['io']['O']['f'].add(biom)
    else:
        outputs['io']['I']['f'].add(biom)
    tsv = '%s.tsv' % splitext(biom)[0]
    if not isfile(tsv):
        cmd = 'biom convert -i %s -o %s.tmp --to-tsv\n' % (biom, tsv)
        cmd += 'tail -n +2 %s.tmp > %s\n' % (tsv, tsv)
        cmd += 'rm %s.tmp\n' % tsv
        outputs['cmds'].append(cmd)
        outputs['io']['O']['f'].add(tsv)
    outputs['outs'].extend([tsv, biom])


def classif_cazy(genes: str, out_dir: str, path: str, outputs: dict):
    """Get the Woltka commands for the cazy-level classification.

    Parameters
    ----------
    genes : str
        Path to the genes classification.
    out_dir : str
        Path to pipeline output folder for SHOGUN.
    path : str
        Path to the woltka database.
    outputs : dict
        All outputs
    """
    cazy_map = '%s/cazymes/cazy/mapping_3tools.txt' % path
    biom = '%s/cazy/wol_cazy.biom' % out_dir
    if not isfile(biom):
        cmd = 'woltka tools collapse'
        cmd += '--input %s' % genes
        cmd += ' --map %s' % cazy_map
        cmd += ' --output %s\n\n' % biom
        outputs['cmds'].append(cmd)
        outputs['io']['O']['f'].add(biom)
    else:
        outputs['io']['I']['f'].add(biom)
    tsv = '%s.tsv' % splitext(biom)[0]
    if not isfile(tsv):
        cmd = 'biom convert -i %s -o %s.tmp --to-tsv\n' % (biom, tsv)
        cmd += 'tail -n +2 %s.tmp > %s\n' % (tsv, tsv)
        cmd += 'rm %s.tmp\n' % tsv
        outputs['cmds'].append(cmd)
        outputs['io']['O']['f'].add(tsv)
    outputs['outs'].extend([biom, tsv])


def classif_metacyc(genes: str, out_dir: str, path: str, outputs: dict):
    """Get the Woltka commands for the metacyc-level classification.

    Parameters
    ----------
    genes : str
        Path to the genes classification.
    out_dir : str
        Path to pipeline output folder for SHOGUN.
    path : str
        Path to the woltka database.
    outputs : dict
        All outputs
    """
    metacyc_dir = '%s/release/function/metacyc' % path
    metacyc = [('protein', 'protein_name.txt', 'protein.map.xz'),
               ('enzrxn', 'enzrxn_name.txt', 'protein-to-enzrxn.txt'),
               ('reaction', 'reaction_name.txt', 'enzrxn-to-reaction.txt'),
               ('pathway', 'pathway_name.txt', 'reaction-to-pathway.txt'),
               ('super_pathway', '', 'pathway-to-super_pathway.txt'),
               ('regulation', '', 'enzrxn-to-regulation.txt'),
               ('regulator', '', 'regulation-to-regulator.txt'),
               ('ec', '', 'reaction-to-ec.txt')]
    files = [genes]
    woltka_fun_out = '%s/metacyc' % out_dir
    cmd = ''
    for idx, (level, names, maps) in enumerate(metacyc):
        input_fp = files[idx]
        if not input_fp.endswith('.biom'):
            input_biom = '%s.biom' % splitext(input_fp)[0]
            if not isfile(input_biom):
                cmd += 'biom convert'
                cmd += ' -i %s' % input_fp
                cmd += ' -o %s' % input_biom
                cmd += ' --to-hdf5'
                cmd += ' --table-type="OTU table"\n'
                outputs['io']['O']['f'].add(input_biom)
            else:
                outputs['io']['I']['f'].add(input_biom)
        else:
            input_biom = input_fp
        cmd += '\n# %s [no stratification]\n' % level
        biom = '%s/metacyc-%s_pre.biom' % (woltka_fun_out, level)
        if not isfile(biom):
            cmd += 'woltka tools collapse'
            if names:
                cmd += ' --input %s' % input_biom
                cmd += ' --names %s/%s' % (metacyc_dir, names)
            else:
                if level == 'super_pathway':
                    cmd += ' --input %s' % files[4]
                elif level == 'regulation':
                    cmd += ' --input %s' % files[2]
                elif level == 'regulator':
                    cmd += ' --input %s' % files[6]
                elif level == 'ec':
                    cmd += ' --input %s' % files[3]
                else:
                    cmd += ' --input %s' % files[4]
            cmd += ' --map %s/%s' % (metacyc_dir, maps)
            cmd += ' --output %s\n' % biom
            outputs['io']['O']['f'].add(biom)
        else:
            outputs['io']['I']['f'].add(biom)
        tsv = '%s/metacyc-%s_pre.tsv' % (woltka_fun_out, level)
        if not isfile(tsv):
            cmd += 'biom convert'
            cmd += ' -i %s' % biom
            cmd += ' -o %s.tmp' % tsv
            cmd += ' --to-tsv\n'
            cmd += 'tail -n +2 %s.tmp > %s\n' % (tsv, tsv)
            cmd += 'rm %s.tmp\n' % tsv
            outputs['io']['O']['f'].add(tsv)
        else:
            outputs['io']['I']['f'].add(tsv)
        files.extend([biom, tsv])
    if cmd:
        outputs['cmds'].append(cmd)
    outputs['outs'].extend(files)


def classif_kegg(uniref: str, out_dir: str, path: str, outputs: dict):
    """Get the Woltka commands for the kegg-level classification.

    Parameters
    ----------
    uniref : str
        Path to the uniref classification.
    out_dir : str
        Path to pipeline output folder for SHOGUN.
    path : str
        Path to the woltka database.
    outputs : dict
        All outputs
    """
    ko_names_maps = [
        ('ko', 'ko.name', 'ko.map.xz', ''),
        ('query', '', '', ''),
        ('ko-cog', '', 'ko-to-cog.txt', 'ko'),
        ('ko-disease', 'disease_name.txt', 'ko-to-disease.txt', 'ko'),
        ('ko-ec', '', 'ko-to-ec.txt', 'ko'),
        ('ko-go', '', 'ko-to-go.txt', 'ko'),
        ('ko-module', 'module_name.txt', 'ko-to-module.txt', 'ko'),
        ('ko-pathway', 'pathway_name.txt', 'ko-to-pathway.txt', 'ko'),
        ('ko-reaction', 'reaction_name.txt', 'ko-to-reaction.txt', 'ko'),
        ('ko-module-class', '', 'module-to-class.txt', 'ko-module'),
        ('ko-module-compound', 'compound_name.txt',
         'module-to-compound.txt', 'ko-module'),
        ('ko-module-ko', 'ko_name.txt', 'module-to-ko.txt', 'ko-module'),
        ('ko-module-pathway', 'pathway_name.txt',
         'module-to-pathway.txt', 'ko-module'),
        ('ko-module-reaction', 'reaction_name.txt',
         'module-to-reaction.txt', 'ko-module'),
        ('ko-pathway-class', '', 'pathway-to-class.txt', 'ko-pathway'),
        ('ko-pathway-compound', 'compound_name.txt',
         'pathway-to-compound.txt', 'ko-pathway'),
        ('ko-pathway-disease', 'disease_name.txt',
         'pathway-to-disease.txt', 'ko-pathway'),
        ('ko-pathway-ko', 'ko_name.txt', 'pathway-to-ko.txt', 'ko-pathway'),
        ('ko-pathway-module', 'module_name.txt',
         'pathway-to-module.txt', 'ko-pathway'),
        ('ko-reaction-compound', 'compound_name.txt',
         'reaction-to-compound.txt', 'ko-reaction'),
        ('ko-reaction-enzyme', '', 'reaction-to-enzyme.txt', 'ko-reaction'),
        ('ko-reaction-ko', 'ko_name.txt', 'reaction-to-ko.txt', 'ko-reaction'),
        ('ko-reaction-left_compound', 'compound_name.txt',
         'reaction-to-left_compound.txt', 'ko-reaction'),
        ('ko-reaction-module', 'module_name.txt',
         'reaction-to-module.txt', 'ko-reaction'),
        ('ko-reaction-pathway', 'pathway_name.txt',
         'reaction-to-pathway.txt', 'ko-reaction'),
        ('ko-reaction-rclass', 'rclass_name.txt',
         'reaction-to-rclass.txt', 'ko-reaction'),
        ('ko-reaction-right_compound', 'compound_name.txt',
         'reaction-to-right_compound.txt', 'ko-reaction')]
    cmd = ''
    files = []
    kegg_maps = '%s/kegg_queried' % out_dir
    for (level, name, maps, prev) in ko_names_maps:
        if maps:
            biom = '%s/kegg/kegg-%s_pre.biom' % (out_dir, level)
            tsv = '%s/kegg/kegg-%s_pre.tsv' % (out_dir, level)
            if not prev:
                if not isfile(tsv):
                    cmd += 'woltka tools collapse'
                    cmd += ' --input %s' % uniref
                    cmd += ' --names %s/release/function/kegg/%s' % (path, name)
                    cmd += ' --map %s/release/function/kegg/%s' % (path, maps)
                    cmd += ' --output %s\n' % biom

                    cmd += ' biom convert -i %s' % biom
                    cmd += ' -o %s.tmp --to-tsv\n' % tsv
                    cmd += ' tail -n +2 %s.tmp\n' % tsv
                    cmd += ' > %s\n' % tsv
                    cmd += ' rm %s.tmp\n' % tsv
                    outputs['io']['O']['f'].update([biom, tsv])
                else:
                    outputs['io']['I']['f'].add(biom)
            else:
                input_fp = '%s/kegg/kegg-%s.biom' % (out_dir, level)
                if not isfile(tsv):
                    cmd += 'woltka tools collapse'
                    cmd += ' --input %s' % input_fp
                    if name:
                        cmd += ' --names %s/%s' % (kegg_maps, name)
                    cmd += ' --map %s/%s' % (kegg_maps, maps)
                    cmd += ' --output %s\n' % biom
                    cmd += 'biom convert -i %s\n' % biom
                    cmd += ' -o %s.tmp --to-tsv\n' % tsv
                    cmd += 'tail -n +2 %s.tmp > %s\n' % (tsv, tsv)
                    cmd += 'rm %s.tmp\n\n' % tsv
                    outputs['io']['O']['f'].update([biom, tsv])
                else:
                    outputs['io']['I']['f'].add(biom)
        else:
            if not isfile('%s/kegg_info.txt' % kegg_maps):
                cmd += 'cd %s\n' % kegg_maps
                cmd += 'cp %s %s/%s\n' % (tsv, kegg_maps, basename(tsv))
                kegg_query = '%s/wol/kegg_query.py' % RESOURCES
                cmd += 'python3 %s %s\n' % (kegg_query, basename(tsv))
                outputs['io']['O']['d'].add(kegg_maps)
            else:
                outputs['io']['I']['d'].add(kegg_maps)
    if cmd:
        outputs['cmds'].append(cmd)
    outputs['outs'].extend(files)


def woltka(prev: str, out_dir: str, inputs: dict, path: str) -> dict:
    """

    Parameters
    ----------
    prev : str
        Previous software
    out_dir : str
        Path to pipeline output folder
    inputs : dict
        Input files
    path : str
        Path to Web of Life database

    Returns
    -------
    outputs : dict
        All outputs
    """
    outputs = {
        'io': {'I': {'d': set(), 'f': set()}, 'O': {'d': set(), 'f': set()}},
        'cmds': [], 'dirs': [], 'outs': []}
    alis = get_aligments(inputs, outputs)
    woltka_map = write_woltka_map(alis, out_dir, outputs)
    get_gotus(out_dir, woltka_map, outputs)
    tax_outmap = get_tax_cmd(woltka_map, out_dir, path, outputs)
    classif_go(woltka_map, out_dir, path, tax_outmap, outputs)
    genes = classif_genes(woltka_map, out_dir, path, outputs)
    uniref = classif_uniref(genes, out_dir, path, outputs)
    classif_eggnog(uniref, out_dir, outputs)
    classif_cazy(genes, out_dir, path, outputs)
    classif_metacyc(genes, out_dir, path, outputs)
    classif_kegg(uniref, out_dir, path, outputs)
    return outputs
