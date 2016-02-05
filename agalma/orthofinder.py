#! /usr/bin/env python
__author__="Andres Perez-Figueroa"
__date__ ="$28-jan-2015 08:00:00$"
__version_="0.0.1"

# Script to copy and prepare the assemblies files for their use in OrthoFinder

import os
import subprocess
from glob import glob
from collections import namedtuple
from agalma import config
from agalma import database
import biolite.database
from biolite import diagnostics
from biolite import report
from biolite import utils
from biolite import wrappers
from biolite import workflows
from biolite.pipeline import Pipeline
from biolite.pipeline import BasePipeline

pipe = Pipeline('orthofinder', "Pipeline for executing orthofinder")

pipe.add_arg('load_ids', nargs='+', metavar='LOAD_ID', help="""
    Specifies which datasets to analyze; A space delimited list of run_id's
    corresponding to the load runs that populated the database""")

pipe.add_arg('--seq_type', '-t', metavar='TYPE', default='protein', help="""
    Choose whether comparisons are based on 'nucleotide', 'protein', or
    'masked_protein' sequences.""")

pipe.add_arg(
    '--genome_type', '-g', nargs="+", metavar='TYPE',
    default=['nuclear'], help="""
        Select sequences with genome type: nuclear, mitochondrial, and/or
        plastid.""")

pipe.add_arg(
    '--molecule_type', '-m', nargs="+", metavar='TYPE',
    default=['protein-coding'], help="""
        Select sequences with molecule type: protein-coding, ribosomal-small,
        ribosomal-large, ribosomal, and/or unknown.""")

SpeciesData = namedtuple('SpeciesData', "name ncbi_id itis_id catalog_id")

### STAGES ###

@pipe.stage
def lookup_species(load_ids): #Taken directly from homologize.py
    """Lookup the species data for each run"""

    species = {}

    # Given a run_id, returns a named tuple with the following elements:
    #
    # name        The species names
    # ncbi_id        The NCBI taxon id
    # itis_id        The ITIS taxon id
    # catalog_id    The agalma catalog id

    for load_id in load_ids:
        row = biolite.database.execute("""
            SELECT catalog.species, catalog.ncbi_id, catalog.itis_id, catalog.id
            FROM catalog, runs
            WHERE runs.run_id=? AND catalog.id=runs.id;""",
            (load_id,)).fetchone()
        if not row:
            utils.die("Couldn't find species data for run ID %s" % load_id)
        if row[0] is None:
            utils.die("Species name is empty for catalog ID '%s'" % row[3])
        species[load_id] = SpeciesData(*row)
        diagnostics.log(str(load_id), row)

    diagnostics.log("species", species)
    ingest('species')


@pipe.stage
def prepare_fasta(id, _run_id,load_ids, species, seq_type,
        molecule_type, genome_type, outdir):
    """ Write fasta file adapted for OrthoFinder """


    fasta_dir = utils.safe_mkdir('fastas_%s' % (_run_id))


    nloads = 0
    nseqs = 0

    for load_id in load_ids:
        nloads +=1
        taxon = species[load_id].name
        fasta = os.path.join(fasta_dir, taxon + '.fa')
        with open(fasta, 'w') as ffasta:
            for record in database.load_seqs(
                    load_id, taxon, seq_type,                                molecule_type, genome_type):
                newid = taxon + "@" + str(record.id)
                nseqs += 1
                utils.write_fasta(ffasta, record.seq, newid)
    if not nseqs:
        utils.die("no sequences were written to the FASTA file")
    diagnostics.log('nloads', nloads)
    diagnostics.log('nseqs', nseqs)


    ingest('fasta_dir', 'fasta')

#@pipe.stage
#def orthofinder(id, _run_id,load_ids, species, outdir):
#    """ Run OrthoFinder """

if __name__ == "__main__":
    # Run the pipeline.
    pipe.run()
    # Push the local diagnostics to the global database.
    diagnostics.merge()
