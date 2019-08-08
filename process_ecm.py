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
from collections import defaultdict

##
# function to merge two tables by row IDs
def merge_tables(table_a, table_b, union=False):
    merged_dict = defaultdict(MergedRow)
    for row in table_a:
        merged_dict[row[0]].a = row.tolist()[1:]

    for row in table_b:
        merged_dict[row[0]].b = row.tolist()[1:]
    
    len_a = len(table_a[0]) - 1
    len_b = len(table_b[0]) - 1
    
    merged_array = list()
    
    for row_id, row in merged_dict.iteritems():
        if union or row.has_both():
            merged_array.append((row_id, ) + row.merge(len_a, len_b))

    merged_dtype = list()
    for i in range(len_a + 1):
        merged_dtype.append( (table_a.dtype.names[i], table_a.dtype[i]) )
    for i in range(1, len_b + 1):
        merged_dtype.append( (table_b.dtype.names[i], table_b.dtype[i]) )

    return numpy.array(merged_array, dtype=merged_dtype)  

##
# Class to contain a merged row.  Used by merge_tables function
class MergedRow:
    
    def __init__(self):
        self.a = None
        self.b = None

    def merge(self, len_a=0, len_b=0):
        if self.a is None:
            self.a = ( numpy.nan, ) * len_a  
        if self.b is None:
            self.b = ( numpy.nan, ) * len_b  
        return self.a + self.b

    def has_both(self):
        return self.a is not None and self.b is not None

class Gene:

    def __init__(self, gene, gene_name, division, category):
        self.gene = gene
        self.gene_name = gene_name
        self.division = division
        self.category = category

class MatrisomeDB:

    # regex patterns for uniprot and refseq IDs
    uniprot_pattern = re.compile('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}')
    refseq_pattern = re.compile('\w{2}_\d{1,}\.\d{1,}')

    empty_gene = Gene('', '', 'Other', 'Other')

    ## 
    # Parse the database file, columns are...
    # Division        Category        Gene Symbol     Gene Name       Synonyms        MGI_IDs MGI_IDs Links   UniProt_IDs     Refseq_IDs      Notes
    def __init__(self, db_file):

        is_io = hasattr(db_file, 'read')
        if not is_io:
            db_file = open(db_file, 'r')

        header = next(db_file).strip().split("\t")
        uniprot_col = header.index("UniProt_IDs")
        refseq_col = header.index("Refseq_IDs")

        self.gene_name = dict()
        self.mgi = dict()
        self.uniprot = dict()
        self.refseq = dict()
        for line in db_file:
            parts = line.strip().split("\t")
            gene = Gene(parts[2], parts[3], parts[0], parts[1])
            self.gene_name[parts[2]] = gene 
            if uniprot_col != -1 and len(parts) > uniprot_col:
                for an_id in parts[uniprot_col].strip().split(':'):
                    if an_id != '':
                        self.uniprot[an_id] = gene
            if refseq_col != -1 and len(parts) > refseq_col:
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
            return self.gene_name.get(anid, MatrisomeDB.empty_gene)
            
# Main routine
if __name__ == "__main__":
    try:
        # Parse arguments from command line
        parser = argparse.ArgumentParser(description="Combine quantitated protein results then annotate using ECM database.")
        parser.add_argument('-p', '--protein_psm', help="Protein PSM file")
        parser.add_argument('-e', '--protein_eic', help="Protein EIC file")
        parser.add_argument('-o', '--output', help="Output table. Default is STDOUT", default='-')
        parser.add_argument('-d', '--db', help="ECM master file", required=True)
        parser.add_argument('-u', '--union', help="Output union of protein files.  Default is intersection", action='store_true', default=False)
        parser.add_argument('-n', '--no_header', help="Do not include column headers in output", action='store_true', default=False)
        parser.add_argument('-s', '--sample_name', help="Sample name to include in output table.")
        opts = parser.parse_args()

        if opts.protein_eic is None and opts.protein_psm is None:
            sys.stderr.write("Must supply either a PSM or EIC table")
            exit(1)

        # Get the PSM table (count of spectra)
        psm_table = None
        if opts.protein_psm is not None:
            psm_table = numpy.genfromtxt(opts.protein_psm, delimiter="\t", comments="#", skip_header=3, usecols=(0,3,4),
                dtype=[('Protein_ID', 'U256'), ('Peptides_spectra', numpy.int8), ('PSM_count', numpy.int8)],
                converters={0: lambda s: s.strip('\"') })
         
        # Get the EIC table (area under the curve)
        eic_table = None
        if opts.protein_eic is not None:
            eic_table = numpy.loadtxt(opts.protein_eic, delimiter="\t", comments="#", skiprows=3, usecols=(0,3,4),
                dtype=[('Protein_ID', 'U256'), ('Peptides_EIC', numpy.int8), ('EIC_intensity', numpy.float64)], 
                converters={0: lambda s: s.strip('\"')})

        if psm_table is not None and eic_table is not None:
            # Need to merge the tables
            table = merge_tables(psm_table, eic_table, opts.union)
            table_format = ["{row[1]:d}", "{row[3]:d}", "{row[2]:d}", "{row[4]:.02f}"]
            header = ( table.dtype.names[1], table.dtype.names[3], table.dtype.names[2], table.dtype.names[4] )
        elif psm_table is not None:
            table = psm_table
            table_format = ["{row[1]:d}", "{row[2]:d}"]
            header = table.dtype.names[1:2]
        else:
            table = eic_table
            table_format = ["{row[1]:d}", "{row[2]:.02f}"]
            header = table.dtype.names[1:2]

        # Setup the output format
        table_format = "\t".join(["{row[0]}", "{gene.gene_name}", "{gene.division}", "{gene.category}"] + table_format)

        # Load database
        db = MatrisomeDB(opts.db) 

        # Write the output
        output = sys.stdout if opts.output == "-" else open(opts.output, 'w')
        if not opts.no_header:
            # Print the header
            if opts.sample_name is not None:
                output.write("## Sample: " + opts.sample_name + "\n")
            output.write("\t".join((table.dtype.names[0], "Gene name", "Division", "Category") + header) + "\n")

        for row in table:
            gene = db.get(row[0])
            output.write(table_format.format(row=row, gene=gene) + "\n")
             
        if opts.output != '-':
            output.close()

    except Exception as e:
        traceback.print_exc(file=sys.stderr)
