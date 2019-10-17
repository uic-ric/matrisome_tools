#!/usr/bin/env python
################################################################################
# Script : ecm_process.py
# Author : George Chlipala
# Created: Feb 13, 2019
# -- Description ----------------------------------------
# Script to process peptide and protein results and then annotate proteins using ECM database
# -- Requirements ---------------------------------------
# Python
# Numpy
################################################################################

import sys, argparse, os, traceback, numpy, re
from collections import defaultdict, Counter

class Protein:
    
    def __init__(self):
        # counters for PSMs (spectra)
        self.unique_psms = 0
        self.nonunique_psms = 0
        self.interspec_psms = 0
        # counters for EIC (integrated area of extracted ion chromatograms)
        self.unique_eic = 0.0
        self.nonunique_eic = 0.0
        self.interspec_eic = 0.0
        # counters for peptides (from PSM data)
        self.unique_psm_peptides = 0
        self.nonunique_psm_peptides = 0
        self.interspec_psm_peptides = 0
        # counters for peptides (from EIC data)
        self.unique_eic_peptides = 0
        self.nonunique_eic_peptides = 0
        self.interspec_eic_peptides = 0

    def add_unique_psms(self, psms):
        self.unique_psms += psms
        self.unique_psm_peptides += 1

    def add_nonunique_psms(self, psms):
        self.nonunique_psms += psms
        self.nonunique_psm_peptides += 1

    def add_interspec_psms(self, psms):
        self.interspec_psms += psms
        self.interspec_psm_peptides += 1

    def add_unique_eic(self, eic):
        self.unique_eic += eic
        self.unique_eic_peptides += 1

    def add_nonunique_eic(self, eic):
        self.nonunique_eic += eic
        self.nonunique_eic_peptides += 1

    def add_interspec_eic(self, eic):
        self.interspec_eic += eic
        self.interspec_eic_peptides += 1

##
# Class to store a list of parsed proteins from the peptide list report
class ProteinList:

    # regex for determining the reference index for a protein match
    db_index = re.compile('^([0-9])::(.+)$')
    
    @staticmethod
    def parse_protein_ids(protein_list):
        for protein in protein_list.strip('"').split("/"):
            a_match = ProteinList.db_index.match(protein)
            if a_match:
                yield (int(a_match.group(1)), a_match.group(2))
            else:
                yield (1, protein)

    def __init__(self):
        self.proteins = defaultdict(Protein)

    def iteritems(self):
        return self.proteins.iteritems()

    # Method to parse a peptide PSM table from OpenMS
    def parse_psm_table(self, table_file):
        psm_table = numpy.genfromtxt(table_file, delimiter="\t", comments="#", 
            skip_header=2, usecols=(1,2,4),
            dtype=('U256', numpy.uint32, numpy.uint32), names=True,
            converters={0: lambda s: s.strip('\"') })

        for row in psm_table:
            protein_list = list(ProteinList.parse_protein_ids(row[0]))
            if row[1] == 1:
                protein_id = protein_list[0][1]
                self.proteins[protein_id].add_unique_psms(row[2])
            else:
                is_interspec = len(set([x[0] for x in protein_list])) > 1
                for protein in protein_list:
                    if is_interspec:
                        self.proteins[protein[1]].add_interspec_psms(row[2])
                    else: 
                        self.proteins[protein[1]].add_nonunique_psms(row[2])

    # Method to parse a peptide EIC table from OpenMS
    def parse_eic_table(self, table_file):
        eic_table = numpy.genfromtxt(table_file, delimiter="\t", comments="#", 
            skip_header=2, usecols=(1,2,4), names=True,
            dtype=('U256', numpy.uint32, numpy.float64), 
            converters={0: lambda s: s.strip('\"')})

        for row in eic_table:
            protein_list = list(ProteinList.parse_protein_ids(row[0]))
            if row[1] == 1:
                protein_id = protein_list[0][1]
                self.proteins[protein_id].add_unique_eic(row[2])
            else:
                is_interspec = len(set([x[0] for x in protein_list])) > 1
                for protein in protein_list:
                    if is_interspec:
                        self.proteins[protein[1]].add_interspec_eic(row[2])
                    else: 
                        self.proteins[protein[1]].add_nonunique_eic(row[2])

##
# class to contain information about a gene record from the matrisome DB
class Gene:

    def __init__(self, gene, gene_name, division, category, species='NA'):
        self.gene = gene
        self.gene_name = gene_name
        self.division = division
        self.category = category
        self.species = species

##
# Function to safely retrieve an index for a item in a list.
# Will return False if the item does not exist in the list
def safe_index(a_list, column):
    try:
        return a_list.index(column)
    except ValueError:
        return False 

class MatrisomeDB:

    # regex patterns for uniprot and refseq IDs
    uniprot_pattern = re.compile('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}')
    uniprot_label_pattern = re.compile('^(sp|tr)\|([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})\|.+')
    refseq_pattern = re.compile('\w{2}_\d{1,}\.\d{1,}')

    empty_gene = Gene('', '', 'Other', 'Other')

    ## 
    # Parse the database file, columns are...
    # Division        Category        Gene Symbol     Gene Name       Species   Synonyms        MGI_IDs MGI_IDs Links   UniProt_IDs     Refseq_IDs      Notes
    def __init__(self, db_file):

        is_io = hasattr(db_file, 'read')
        if not is_io:
            db_file = open(db_file, 'r')

        header = next(db_file).strip().split("\t")
        uniprot_col = safe_index(header, "UniProt_IDs")
        refseq_col = safe_index(header, "Refseq_IDs")
        spec_col = safe_index(header, "Species")

        if spec_col:
            self.mixed_species = True
        else:
            self.mixed_species = False

        self.gene_name = dict()
        self.mgi = dict()
        self.uniprot = dict()
        self.refseq = dict()
        for line in db_file:
            parts = line.strip().split("\t")
            gene = Gene(parts[2], parts[3], parts[0], parts[1])
            if spec_col:
                gene.species = parts[spec_col].strip()
            self.gene_name[parts[2]] = gene 
            if uniprot_col and len(parts) > uniprot_col:
                for an_id in parts[uniprot_col].strip().split(':'):
                    if an_id != '':
                        self.uniprot[an_id] = gene
            if refseq_col and len(parts) > refseq_col:
                for an_id in parts[refseq_col].strip().split(':'):
                    if an_id != '':
                        self.refseq[an_id] = gene

        if not is_io:
            db_file.close()

    ##
    # Get a gene record for a protein ID.
    # Accepts Uniprot, Refseq and gene names as IDs
    def get(self, anid):
        if MatrisomeDB.uniprot_pattern.match(anid):
            return self.uniprot.get(anid, MatrisomeDB.empty_gene)
        elif MatrisomeDB.refseq_pattern.match(anid):
            return self.refseq.get(anid, MatrisomeDB.empty_gene)
        else:
            match = MatrisomeDB.uniprot_label_pattern.match(anid)
            if match:
                return self.uniprot.get(match.group(2), MatrisomeDB.empty_gene)
            else:
                return self.gene_name.get(anid, MatrisomeDB.empty_gene)
            
# Main routine
if __name__ == "__main__":
    try:
        # Parse arguments from command line
        parser = argparse.ArgumentParser(description="Combine quantitated protein results then annotate using ECM database.")
        parser.add_argument('-p', '--peptide_psm', help="Peptide PSM file")
        parser.add_argument('-e', '--peptide_eic', help="Peptide EIC file")
        parser.add_argument('-o', '--output', help="Output table. Default is STDOUT", default='-')
        parser.add_argument('-d', '--db', help="ECM master file", required=True)
        parser.add_argument('-u', '--union', help="Output union of protein files.  Default is intersection", action='store_true', default=False)
        parser.add_argument('-n', '--no_header', help="Do not include column headers in output", action='store_true', default=False)
        parser.add_argument('-s', '--sample_name', help="Sample name to include in output table.")
        opts = parser.parse_args()

        if opts.peptide_psm is None and opts.peptide_eic is None:
            sys.stderr.write("Must supply either a PSM or EIC table")
            exit(1)

        protein_list = ProteinList()

        # Get the PSM table (count of spectra)
        if opts.peptide_psm is not None:
            protein_list.parse_psm_table(opts.peptide_psm)
         
        # Get the EIC table (area under the curve)
        if opts.peptide_eic is not None:
            protein_list.parse_eic_table(opts.peptide_eic)

        # Load database
        db = MatrisomeDB(opts.db) 

        # Setup the output format
        table_format = ["{protein_id}", "{gene.gene}", "{gene.gene_name}", "{gene.division}", "{gene.category}"]
        if db.mixed_species:
            table_format.append("{gene.species}")

        if opts.peptide_psm is not None and opts.peptide_eic is not None:
            table_format.extend(["{protein.unique_psm_peptides}", "{protein.nonunique_psm_peptides}", "{protein.interspec_psm_peptides}",
                "{protein.unique_psms}", "{protein.nonunique_psms}", "{protein.interspec_psms}",
                "{protein.unique_eic}", "{protein.nonunique_eic}", "{protein.interspec_eic}"])
            header = ( "Unique peptides", "Non-unique peptides", "Interspecies peptides", 
                "Unique spectra", "Non-unique spectra", "Interspecies spectra",
                "Unique EIC", "Non-unique EIC", "Interspecies EIC")
        elif opts.peptide_psm is not None:
            table_format.extend(["{protein.unique_psm_peptides}", "{protein.nonunique_psm_peptides}", "{protein.interspec_psm_peptides}",
                "{protein.unique_psms}", "{protein.nonunique_psms}", "{protein.interspec_psms}"])
            header = ( "Unique peptides", "Non-unique peptides", "Interspecies peptides", 
                "Unique spectra", "Non-unique spectra", "Interspecies spectra")
        else:
            table_format.extend(["{protein.unique_eic_peptides}", "{protein.nonunique_eic_peptides}", "{protein.interspec_eic_peptides}",
                "{protein.unique_eic}", "{protein.nonunique_eic}", "{protein.interspec_eic}"])
            header = ( "Unique peptides", "Non-unique peptides", "Interspecies peptides", 
                "Unique EIC", "Non-unique EIC", "Interspecies EIC")

        table_format = "\t".join(table_format)

        # Write the output
        output = sys.stdout if opts.output == "-" else open(opts.output, 'w')
        if not opts.no_header:
            # Print the header
            if opts.sample_name is not None:
                output.write("## Sample: " + opts.sample_name + "\n")
            output.write("\t".join(("ProteinID", "Gene symbol", "Gene name", "Division", "Category")))
            if db.mixed_species:
                output.write("\tSpecies")
            output.write("\t" + "\t".join(header) + "\n")

        for protein_id, protein in protein_list.iteritems():
            gene = db.get(protein_id)
            output.write(table_format.format(protein_id=protein_id, gene=gene, protein=protein) + "\n")
             
        if opts.output != '-':
            output.close()

    except Exception as e:
        traceback.print_exc(file=sys.stderr)
