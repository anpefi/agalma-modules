# Modification of the homologize.py from Agalma to include filtering by SRV
# (Score Ratio Values) instead of raw Bit-scores. It also allows for change the
# Inflation parameter of MCL. nucleotide secuences are compared by means of
# blastn.
# Andrés Pérez-Figueroa. 2016.

# Agalma - Tools for processing gene sequence data and automating workflows
# Copyright (c) 2012-2014 Brown University. All rights reserved.
#
# This file is part of Agalma.
#
# Agalma is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Agalma is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Agalma.  If not, see <http://www.gnu.org/licenses/>.

"""
Identifies homologous sequences across datasets. Takes assembly or group of
assemblies and prepares a set of comprehensive comparisons between them. First
an all by all BLAST is run with a stringent threshold, and the hits which match
above a given score are used as edges between two transcripts which then form a
graph. This graph is the basis of a series of comparative scores.
"""

import os
import numpy
import networkx as nx
from collections import namedtuple
from glob import glob
from agalma import config
from agalma import database
import biolite.database
from biolite import diagnostics
from biolite import report
from biolite import utils
from biolite import wrappers
from biolite import workflows
from biolite.pipeline import Pipeline

pipe = Pipeline('homologize', __doc__)

### ARGUMENTS ###

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

pipe.add_arg('--min_bitscore', type=float, default=0.0, help="""
	Filter out BLAST hit edges with a SRV (Blast Score Ratio Values) below min_bitscore.""")

pipe.add_arg('--min_overlap', type=float, default=0.5, help="""
	Filter out BLAST hit edges with a max HSP length less than this fraction
	of target length.""")

pipe.add_arg('--min_nodes', type=int, default=4, help="""
	Minimum number of nodes to be retained as a homolgous gene cluster.""")

pipe.add_path_arg('--genesets', nargs='+', default=None, help="""
	Concatenate the specified amino acid FASTAs and use as the reference for
	a blastx or blastp comparison, instead of performing an all-to-all
	comparison. FASTA headers must start with a non-digit.""")

pipe.add_path_arg('--blast_hits', default=None, help="""
	Use a pre-computed blast alignment when starting from the parse_edges
	stage.""")

pipe.add_arg('--mcl_inflation', type=float, default=2.1, help="""
 	Sets the main inflation value for the mcl algorithm. This value is the
	main handle for affecting cluster granularity. It is usually chosen
	somewhere in the range [1.2-5.0].""")

SpeciesData = namedtuple('SpeciesData', "name ncbi_id itis_id catalog_id")

### STAGES ###

@pipe.stage
def lookup_species(load_ids):
	"""Lookup the species data for each run"""

	species = {}

	# Given a run_id, returns a named tuple with the following elements:
	#
	# name		The species names
	# ncbi_id		The NCBI taxon id
	# itis_id		The ITIS taxon id
	# catalog_id	The agalma catalog id

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
def write_fasta(
		id, _run_id, load_ids, species, seq_type,
		molecule_type, genome_type, outdir):
	"""Write sequences from the Agalma database to a FASTA file"""

	blast_dir = utils.safe_mkdir('allvall_blast_%s_%s' % (id, _run_id))
	fasta = os.path.join(blast_dir, 'all.fa')

	# The nodes file contains a header, which describes the attributes, as well
	# as a line for every node (transcript) in the analysis.
	nodes = os.path.join(outdir, 'nodes.txt')

	nloads = 0
	nseqs = 0
	nbases = 0

	with open(nodes, 'w') as fnodes, open(fasta, 'w') as ffasta:
		print >>fnodes, "label\tid\tassembly\tassembly_number"
		for load_id in load_ids:
			nloads += 1
			taxon = species[load_id].name
			for record in database.load_seqs(
									load_id, taxon, seq_type,
									molecule_type, genome_type):
				nseqs += 1
				nbases += len(record.seq)
				# The id of the node (the second column) is a unique identifier
				# and because we already have one in the table, use that here.
				print >>fnodes, "%s\t%d\t%s\t%d" % (
				            record.header, record.id, load_id, nloads)
				utils.write_fasta(ffasta, record.seq, record.id)

	if not nseqs:
		utils.die("no sequences were written to the FASTA file")

	diagnostics.log_path(nodes, 'nodes')

	diagnostics.log('nloads', nloads)
	diagnostics.log('nseqs', nseqs)
	diagnostics.log('nbases', nbases)

	ingest('blast_dir', 'fasta', 'nodes')


@pipe.stage
def prepare_blast(blast_dir, fasta, seq_type, genesets):
	"""Prepare all-by-all BLAST database and command list"""

	db = os.path.join(blast_dir, 'db')

	if genesets:
		if seq_type == 'nucleotide':
			program = 'blastx'
		else:
			program = 'blastp'
		for geneset in genesets:
			utils.cat_to_file(geneset, db+'.fa')
			# Also add genesets to the query file, to find overlapping genes
			# between sets.
			utils.cat_to_file(geneset, fasta)
		wrappers.MakeBlastDB(db+'.fa', db, 'prot')
	else:
		if (seq_type == 'masked_protein') or (seq_type == 'protein'):
			dbtype = 'prot'
			program = 'blastp'
		elif seq_type == 'nucleotide':
			dbtype = 'nucl'
			program = 'blastn'
		else:
			utils.die("unrecognized sequence type '%s'" % seq_type)
		wrappers.MakeBlastDB(fasta, db, dbtype)

	command = program + " -db db -evalue 1e-20" \
	                  + " -outfmt \"6 qseqid sseqid bitscore qlen length\""
	commands = workflows.blast.split_query(fasta, command, 100000, blast_dir)

	ingest('commands')


@pipe.stage
def run_blast(blast_dir, commands, outdir):
	"""Run all-by-all BLAST"""

	blast_hits = os.path.join(blast_dir, 'hits.tab')
	wrappers.Parallel(
		commands, "--joblog blast.log --resume-failed --halt 1",
		stdout=blast_hits, cwd=blast_dir)
	ingest('blast_hits')


@pipe.stage
def parse_edges(_run_id, seq_type, blast_hits, min_overlap, min_bitscore, min_nodes):
	"""Parse BLAST hits into edges weighted by bitscore"""

	graph = nx.Graph()
	edge_file = 'allvall_edges_%s.abc' % _run_id

	nseqs = 0
	nedges = {
		'all': 0,
		'non-self': 0,
		'passed-overlap': 0,
		'passed-bitscore': 0,
		'passed-bitscore-unique': 0}
	max_bitscore = None
	last_query = "" # For identification of a new query in the table


	for f in glob(blast_hits):
		for line in open(f):
			# fields:
			# qseqid sseqid bitscore qlen length
			id_from, id_to, bitscore, qlen, length = line.rstrip().split()
			# Sometimes blast outputs a bad query id if the query is longer
			# than 10Kb
			if id_from.startswith('Query_'):
				utils.info("discarding bad hit with query id '%s'" % id_from)
				continue
			bitscore = float(bitscore)
			length = float(length) / float(qlen)
			# Correct for nucleotide vs. amino acid length  ### WHY? lenght here is the overlap proportion
			#if seq_type == 'nucleotide'
			#	length *= 3.0
			# Filter out self hits, low scoring hits, and short hits
			nedges['all'] += 1
			if id_from != last_query:
				max_bitscore = bitscore ## The self-hit is not always at the top as the file is sorted by e-values and ties happens.
				last_query = id_from	## however the self-hit always have the same bitscore as the first hit.
			if id_from != id_to:
				nedges['non-self'] += 1
				if length > min_overlap:
					nedges['passed-overlap'] += 1
					srv = float(bitscore) / float(max_bitscore)  #SRV: BLAST Score Ratio Values. Blom et al. BMC Bioinformatics 2009, 10:154 doi:10.1186/1471-2105-10-154
					if srv > min_bitscore:
						nedges['passed-bitscore'] += 1
						# If an edge already exists between the nodes, update
						# its score to be the max
						if graph.has_edge(id_from, id_to):
							e = graph.edge[id_from][id_to]
							e['score'] = max(e['score'], bitscore)
						else:
							nedges['passed-bitscore-unique'] += 1
							graph.add_node(id_from)
							graph.add_node(id_to)
							graph.add_edge(id_from, id_to, score=bitscore)

	diagnostics.prefix.append('nedges')
	diagnostics.log_dict(nedges)
	diagnostics.prefix.pop()

	with open(edge_file, 'w') as f:
		for subgraph in nx.connected_component_subgraphs(graph):
			nnodes = subgraph.number_of_nodes()
			if nnodes >= min_nodes:
				nseqs += nnodes
				for id_from, id_to in subgraph.edges_iter():
					print >>f, "%s\t%s\t%f" % (id_from, id_to, subgraph.edge[id_from][id_to]['score'])

	if not nseqs:
		utils.die("no sequences were written to the FASTA file")

	diagnostics.log('nseqs', nseqs)

	ingest('edge_file')


@pipe.stage
def mcl_cluster(edge_file, mcl_inflation, outdir):
	"""Run mcl on all-by-all graph to form gene clusters"""
	cluster_file = os.path.join(outdir, 'mcl_gene_clusters.abc')
	wrappers.Mcl(edge_file, "--abc -I ", str(mcl_inflation), " -o", cluster_file)
	ingest('cluster_file')


@pipe.stage
def load_mcl_cluster(_run_id, cluster_file):
	"""Load cluster file from mcl into homology database"""

	nseqs = 0
	hist = {}

	database.execute("BEGIN")
	database.execute("DELETE FROM homology WHERE run_id=?;", (_run_id,))
	with open(cluster_file, 'r') as f:
		cluster_id = 0
		for line in f:
			cluster = filter(lambda s : s[0].isdigit(), line.rstrip().split())
			n = len(cluster)
			hist[n] = hist.get(n, 0) + 1
			if n >= 4:
				nseqs += n
				for seq_id in cluster:
					database.execute("""
						INSERT INTO homology (run_id, component_id, sequence_id)
						VALUES (?,?,?);""",
						(_run_id, cluster_id, seq_id))
				cluster_id += 1
	database.execute("COMMIT")

	utils.info(
		"histogram of gene cluster sizes:\n",
		'\n '.join("%d\t:\t%d" % (k, hist[k]) for k in sorted(hist)))

	diagnostics.log('nseqs', nseqs)
	diagnostics.log('histogram', hist)


### RUN ###

if __name__ == "__main__":
	# Run the pipeline.
	pipe.run()
	# Push the local diagnostics to the global database.
	diagnostics.merge()

### REPORT ###

class Report(report.BaseReport):
	def init(self):
		self.name = pipe.name
		# Lookups
		self.lookup('load_ids', diagnostics.INIT, 'load_ids')
		self.str2list('load_ids')
		self.lookup('species', 'homologize.lookup_species')
		self.lookup('nseqs', 'homologize.write_fasta', 'nseqs')
		# Generators
		self.generator(self.species_table)
		self.generator(self.component_histogram)

	def species_table(self):
		"""
		Summary of all species processed.
		"""
		html = []
		if self.check('nseqs'):
			html.append("<p>Total sequences: %s</p>" % self.data.nseqs)
		if self.check('load_ids', 'species'):
			headers = ("Load ID", "Catalog ID", "Species", "NCBI ID", "ITIS ID")
			rows = []
			for load_id in self.data.load_ids:
				row = self.data.species.get(str(load_id), '()')
				row = diagnostics.str2list(row)
				if len(row) == 4:
					rows.append((load_id, row[3], row[0], row[1], row[2]))
				else:
					rows.append((load_id, '-', '-', '-', '-'))
			html += self.table(rows, headers)
		return html

	def component_histogram(self):
		"""
		Distribution of the number of nodes in each cluster.
		"""
		hist = {}
		sql = """
			SELECT component_id, COUNT(*)
			FROM homology
			WHERE run_id=?
			GROUP BY component_id;"""
		for _, count in database.execute(sql, (self.run_id,)):
			hist[count] = hist.get(count, 0) + 1
		if hist:
			hist = numpy.array([(k,hist[k]) for k in sorted(hist.iterkeys())])
			imgname = "%d.component.hist.png" % self.run_id
			props = {
				'title': "Distribution of Cluster Sizes",
				'xlabel': "# Nodes in Cluster",
				'ylabel': "Frequency"}
			return [self.histogram_overlay(imgname, [hist], props=props)]

# vim: noexpandtab ts=4 sw=4
