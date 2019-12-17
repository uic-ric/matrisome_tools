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
import logging
import urllib3
from collections import defaultdict, Counter
from uniprot import UniprotAPI

##
# Function to safely retrieve an index for a item in a list.
# Will return False if the item does not exist in the list
def safe_index(a_list, column):
    try:
        return a_list.index(column)
    except ValueError:
        return False 

class Protein:
    
    def __init__(self, gene=None):
        self.gene = gene
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
# class to store peptide entries
class Peptide:

    def __init__(self, peptide_seq=""):
        self.peptide_seq = peptide_seq
        self.protein_ids = list()
        self.psms = 0
        self.eic_int = 0

##
# Class to store a list of parsed peptides
class PeptideList:

    # regex for determining the reference index for a protein match
    db_index = re.compile('^([0-9]+)::(.+)$')
    
    # Method to parse the proteins IDs from an OpenMS peptide output
    @staticmethod
    def parse_protein_ids(protein_list):
        protein_ids = protein_list.strip("\"").split("/")
        for protein in protein_ids:
            if "::" in protein:
                yield protein.split("::")[1]
            else:
                yield protein

    def __init__(self):
        self.peptides = defaultdict(Peptide)
        self.protein_ids = set()

    def iteritems(self):
        return self.peptides.iteritems()

    # Method to parse a peptide PSM table from OpenMS
    def parse_psm_table(self, table_file):
        logging.debug("Parsing PSM table in file {}".format(table_file))
        psm_table = numpy.genfromtxt(table_file, delimiter="\t", comments="#", 
            skip_header=2, usecols=(0,1,2,4),
            dtype=("U4096", "U4096", numpy.uint32, numpy.uint32), names=True,
            converters={0: lambda s: s.strip("\"") })

        for row in psm_table:
            if len(self.peptides[row[0]].protein_ids) == 0:
                these_protein_ids = list(PeptideList.parse_protein_ids(row[1]))
                self.peptides[row[0]].protein_ids = these_protein_ids
                self.protein_ids = self.protein_ids.union(these_protein_ids)
            self.peptides[row[0]].psms = row[3]

    # Method to parse a peptide EIC table from OpenMS
    def parse_eic_table(self, table_file):
        logging.debug("Parsing EIC table in file {}".format(table_file))
        eic_table = numpy.genfromtxt(table_file, delimiter="\t", comments="#", 
            skip_header=2, usecols=(0,1,2,4), names=True,
            dtype=("U4096", "U4096", numpy.uint32, numpy.float64), 
            converters={0: lambda s: s.strip("\"")})

        for row in eic_table:
            if len(self.peptides[row[0]].protein_ids) == 0:
                these_protein_ids = list(PeptideList.parse_protein_ids(row[1]))
                self.peptides[row[0]].protein_ids = these_protein_ids
                self.protein_ids = self.protein_ids.union(these_protein_ids)
            self.peptides[row[0]].eic_int = row[3]


    # Method to generate a protein list, dictionary of Protein objects index by protein ID, 
    # from the parsed peptide lists and annotate using the Matrisome DB
    def gen_protein_list(self, db):
        logging.debug("Generating protein list for {} peptides.".format(len(self.peptides)))
        proteins = defaultdict(Protein)
        for peptide in self.peptides.values():
            gene_records = [ db.get(protein_id) for protein_id in peptide.protein_ids ] 
            if len(gene_records) == 1:
                if peptide.psms != 0:
                    proteins[peptide.protein_ids[0]].add_unique_psms(peptide.psms)
                if peptide.eic_int != 0:
                    proteins[peptide.protein_ids[0]].add_unique_eic(peptide.eic_int)
                if proteins[peptide.protein_ids[0]].gene is None:
                    proteins[peptide.protein_ids[0]].gene = gene_records[0]
            else:
                organisms = set([ gene.organism for gene in gene_records ])
                is_interspec = len(organisms) > 1
                for p, protein_id in enumerate(peptide.protein_ids):
                    if peptide.psms != 0:
                        if is_interspec:
                            proteins[protein_id].add_interspec_psms(peptide.psms)
                        else:
                            proteins[protein_id].add_nonunique_psms(peptide.psms)
                    if peptide.eic_int != 0:
                        if is_interspec:
                            proteins[protein_id].add_interspec_eic(peptide.eic_int)
                        else:
                            proteins[protein_id].add_nonunique_eic(peptide.eic_int)
                    if proteins[protein_id].gene is None:
                        proteins[protein_id].gene = gene_records[p]
        logging.debug("Found {} proteins.".format(len(proteins)))
        return proteins

##
# class to contain information about a gene record from the matrisome DB
class Gene:

    def __init__(self, gene, division="Other", category="Other", description=None, organism=None):
        self.gene = gene
        self.division = division
        self.category = category
        self.description = description
        self.organism = organism

class MatrisomeDB:

    # regex patterns for uniprot and refseq IDs
    UNIPROT_PATTERN = re.compile('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}')
    UNIPROT_LABEL_PATTERN = re.compile('^(sp|tr)\|([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})\|.+')
    REFSEQ_PATTERN = re.compile('\w{2}_\d{1,}\.\d{1,}')

    ## 
    # Parse the database file, columns are...
    # Gene Symbol   Division        Category
    def __init__(self, db_file, protein_ids):
        logging.debug("Loading Matrisome DB from {}".format(db_file))
        # Check if the db_file is an IO object or a filename
        is_io = hasattr(db_file, "read")
        if not is_io:
            db_file = open(db_file, "r")

        header = next(db_file).strip().split("\t")

        symbol_index = safe_index(header, "Gene Symbol")
        division_index = safe_index(header, "Division")
        category_index = safe_index(header, "Category")

        self.gene_db = dict()
        for line in db_file:
            parts = line.strip().split("\t")
            self.gene_db[parts[symbol_index]] = Gene(parts[symbol_index], parts[division_index], parts[category_index])

        if not is_io:
            db_file.close()

        self.load_ids(protein_ids)

    ##
    # Load a set of IDs into the database
    def load_ids(self, protein_ids, retries=3):
        uniprot_ids = list()
        refseq_ids = list()
        for protein_id in protein_ids:
            if MatrisomeDB.UNIPROT_PATTERN.match(protein_id):
                uniprot_ids.append(protein_id)
            elif MatrisomeDB.REFSEQ_PATTERN.match(protein_id):
                refseq_ids.append(protein_id)
            else:
                match = MatrisomeDB.UNIPROT_LABEL_PATTERN.match(protein_id)
                if match:
                    uniprot_ids.append(match.group(2)) 
        if len(refseq_ids) > 0:
            attempt = 0
            while attempt < retries:
                try:
                    self.refseq_map = UniprotAPI.convert("P_REFSEQ_AC", refseq_ids)
                    uniprot_ids.extend(self.refseq_map.values())
                    break
                except urllib3.exceptions.ProtocolError:
                    attempt += 1
                    if attempt < retries:
                        logging.debug("Error mapping Uniprot IDs.  Will try again.")
                    else:
                        logging.debug("Failed to map Uniprot IDs.")
                        raise 
                except:
                    raise
        uniprot_ids = filter(lambda x: x != "", uniprot_ids)  
        self.uniprot_data = dict()
        attempt = 0
        while attempt < retries:
            try:
                for entry in UniprotAPI.get(uniprot_ids):
                    for acc_id in entry.accession_ids:
                        self.uniprot_data[acc_id] = entry
                break
            except urllib3.exceptions.ProtocolError:
                attempt += 1
                if attempt < retries:
                    logging.debug("Error mapping Uniprot IDs.  Will try again.")
                else:
                    logging.debug("Failed to map Uniprot IDs.")
                    raise 
            except:
                raise

    ##
    # Get a gene record for a protein ID.
    # Accepts Uniprot, Refseq and gene names as IDs
    def get(self, anid):
        # Check if it is a refseq ID.  If so, map to Uniprot ID
        if MatrisomeDB.REFSEQ_PATTERN.match(anid):
            anid = self.refseq_map.get(anid, None)
        else:
            # Otherwise, check if it a Swissprot/Trembl style and strip out the Uniprot ID
            match = MatrisomeDB.UNIPROT_LABEL_PATTERN.match(anid)
            if match:
                anid = match.group(2)

        if anid is not None:
            uniprot_entry = self.uniprot_data.get(anid, None)
            if uniprot_entry is not None:
                gene_name = str(uniprot_entry.get_genename().name)
                gene_record = self.gene_db.get(gene_name, Gene(gene_name))
                if gene_record.description is None:
                    gene_record.description = uniprot_entry.get_recname()
                if gene_record.organism is None:
                    gene_record.organism = uniprot_entry.get_organism()
                return gene_record

        return Gene("")
                    

# Main routine
if __name__ == "__main__":
    try:
        # Parse arguments from command line
        parser = argparse.ArgumentParser(description="Combine quantitated protein results then annotate using ECM database.")
        parser.add_argument("-p", "--peptide_psm", help="Peptide PSM file")
        parser.add_argument("-e", "--peptide_eic", help="Peptide EIC file")
        parser.add_argument("-o", "--output", help="Output table. Default is STDOUT", default="-")
        parser.add_argument("-d", "--db", help="ECM master file", required=True)
        parser.add_argument("-u", "--union", help="Output union of protein files.  Default is intersection", action="store_true", default=False)
        parser.add_argument("-n", "--no_header", help="Do not include column headers in output", action="store_true", default=False)
        parser.add_argument("-s", "--sample_name", help="Sample name to include in output table.")
        parser.add_argument("-v", "--verbose", help="Generate verbose output.", action="store_true", default=False)
        opts = parser.parse_args()

        if opts.peptide_psm is None and opts.peptide_eic is None:
            sys.stderr.write("Must supply either a PSM or EIC table")
            exit(1)

        if opts.verbose:
            logging.basicConfig(level=logging.DEBUG)

        peptide_list = PeptideList()

        # Get the PSM table (count of spectra)
        if opts.peptide_psm is not None:
            peptide_list.parse_psm_table(opts.peptide_psm)
         
        # Get the EIC table (area under the curve)
        if opts.peptide_eic is not None:
            peptide_list.parse_eic_table(opts.peptide_eic)

        # Load matrisome database for the parsed protein IDs
        db = MatrisomeDB(opts.db, peptide_list.protein_ids) 

        # Generate a protein list using the database
        # Protein list will be a dict of Protein objects indexed by protein ID (Uniprot accession ID)
        protein_list = peptide_list.gen_protein_list(db)

        # Setup the output format
        table_format = ["{protein_id}", "{gene.gene}", "{gene.description}", "{gene.division}", "{gene.category}", "{gene.organism}" ]

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
        output = sys.stdout if opts.output == "-" else open(opts.output, "w")
        if not opts.no_header:
            # Print the header
            if opts.sample_name is not None:
                output.write("## Sample: " + opts.sample_name + "\n")
            output.write("\t".join(("ProteinID", "Gene symbol", "Gene name", "Division", "Category", "Species")))
            output.write("\t" + "\t".join(header) + "\n")

        # Write out the protein data
        for protein_id, protein in protein_list.iteritems():
            output.write(table_format.format(protein_id=protein_id, gene=protein.gene, protein=protein) + "\n")
             
        if opts.output != "-":
            output.close()

    except Exception as e:
        traceback.print_exc(file=sys.stderr)
