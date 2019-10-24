#!/usr/bin/env python
################################################################################
# Script : uniprot_lib.py
# Author : George Chlipala
# Created: Oct 22, 2019
# -- Description ----------------------------------------
# Custom python library for parsing Uniprot entries
# -- Requirements ---------------------------------------
# Python 2.7.6
# Requests library
################################################################################

__author__ = "George Chlipala"
__credits__ = ["George Chlipala"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "George Chlipala"
__email__ = "gchlip2@uic.edu"
__status__ = "Production"

import re
import urllib3
import cStringIO
import logging
from collections import defaultdict

# Supress SNI error messages.
urllib3.disable_warnings()

class UniprotAPI:

    # Base URL for Uniprot
    UNIPROT_API_BASE = "https://www.uniprot.org"
    # Number of entries to process at a single time.
    chunk_size = 5000
    logger = logging.getLogger('uniprot.UniprotAPI')

    ##
    # Function to get a set of Uniprot entries
    @staticmethod
    def get(*entries):
        all_entries = list()
        for e in entries:
            if isinstance(e, basestring):
                all_entries.append(e)
            else:
                all_entries.extend(e)
        # Check if the entry list is greater than the chunk size (5000)
        # If so, perform the request a chunk at a time.
        response_data = ""
        e_len = len(all_entries)
        UniprotAPI.logger.debug("Retrieving Uniprot information for {} entries.".format(e_len))
        if e_len > UniprotAPI.chunk_size:
            UniprotAPI.logger.debug("Entry list greater than chunk size of {}. Will retrieve in batches.".format(UniprotAPI.chunk_size))
            for n in range(0, e_len, UniprotAPI.chunk_size):
                max_n = min(n + UniprotAPI.chunk_size, e_len)
                response_data += UniprotAPI.__get_entries__(all_entries[n:max_n])
        else:
            response_data += UniprotAPI.__get_entries__(all_entries)
        
        return UniprotEntry.parse(cStringIO.StringIO(response_data))

    ##
    # Hidden method to actually perform the get request.
    @staticmethod
    def __get_entries__(entries, from_type="ACC+ID", output_format="txt"):
        http = urllib3.PoolManager(headers={'User-Agent': 'urllib3-python'})
        req = http.request('POST', UniprotAPI.UNIPROT_API_BASE + "/uploadlists/", redirect=False, 
            fields={ "file": ("entries.txt", "\n".join(entries) + "\n", "text/plain"),
                "format": output_format, "from": from_type, "to": "ACC" })
        UniprotAPI.logger.debug("Performing HTTP request to Uniprot with format={}, to={}, from={} and {} entries.".format(output_format, from_type, "ACC", len(entries)))
        # The following is a work around for the redirect.
        # The orginal request is performed as a POST, however Uniprot then sends a 302 (redirect)
        # to a URL that should be accessed via GET. 
        # urllib3 tries the redirect as a POST and Uniprot generates a 500 (server error)
        UniprotAPI.logger.debug("HTTP status: {}".format(req.status))
        if req.status == 302:
            location = req.get_redirect_location()
            if not location.startswith('http'):
                location = UniprotAPI.UNIPROT_API_BASE + location
            req = http.request('GET', location)
        return req.data

    ##
    # Function to convert IDs
    @staticmethod
    def convert(id_type, *ids):
        all_ids = list()
        for i in ids:
            if isinstance(i, basestring):
                all_ids.append(i)
            else:
                all_ids.extend(i)

        id_map = dict()
        e_len = len(all_ids)
        UniprotAPI.logger.debug("Converting {} IDs from {} to Uniprot accession IDs.".format(e_len, id_type))
        if e_len > UniprotAPI.chunk_size:
            UniprotAPI.logger.debug("Entry list greater than chunk size of {}. Will retrieve in batches.".format(UniprotAPI.chunk_size))
            for n in range(0, e_len, UniprotAPI.chunk_size):
                max_n = min(n + UniprotAPI.chunk_size, e_len)
                id_map.update(UniprotAPI.convert(id_type, all_ids[n:max_n]))
        else:
            data = UniprotAPI.__get_entries__(all_ids, id_type, "tab") 
            # http = urllib3.PoolManager(headers={'User-Agent': 'urllib3-python'})
            # req = http.request('POST', UniprotAPI.UNIPROT_API_BASE + "/uploadlists/", encode_multipart=False,
            #    fields={ "format": "tab", "from": id_type, "to": "ACC", "query" : " ".join(all_ids) })
            data = cStringIO.StringIO(data)
            headers = next(data)
            for line in data:
                parts = line.strip().split("\t")
                id_map[parts[0]] = parts[1]
        return id_map
    
##
# Class for Uniprot entries
class UniprotEntry(object):

    # regexs need for parsing
    DE_TYPE = re.compile("^ *(Rec|Alt|Sub)Name: *([A-Za-z]+)=(.+);$")
    SEMICOLON_SEP = re.compile('; +')
    COMMA_SEP = re.compile(', *')

    ##
    # Method to create a new entry object.  Typically only called by parsing method
    def __init__(self):
        self.entry_id = ""
        self.status = ""
        self.accession_ids = list()
        self.record_names = { "Rec": dict(), "Alt": dict(), "Sub": dict() }
        self.contains_names = list()
        self.includes_names = list()
        self.name_flags = set()
        self.gene_names = list()
        self.organism = None
        self.sequence = ""
        self.references = list()
        self.comments = dict()
        self.db_refs = defaultdict(list)
        self.evidence = (0, '')
        self.keywords = list()

    ##
    # static method to parse entries from an I/O handle.
    # This method acts as a generator
    @staticmethod
    def parse(handle):
    
        record = UniprotEntry()
        ref_num = 0
        comment_topic = None
        subname = 0 
        gn_data = ""
        
        for line in handle: 
            # If the line is "//", then it is the end of a record.
            if line[:2] == "//": 
                # Add the current GN data 
                record.gene_names.append(UniprotGeneName(gn_data))
                # Sort all references by ref_number
                record.references.sort(key=lambda ref: ref.ref_number)
                yield record 
                record = UniprotEntry() 
                ref_num = 0
                comment_topic = None
                subname = 0
                gn_data = ""
                continue

            if line.startswith(" "):
                record.sequence += line[6:].strip().replace(" ", "")
                continue

            line_code = line[:2]
            data = line[5:].strip() 

            try:
                if line_code == "ID":
                    # Identification line:  ID   EntryName Status; SequenceLength.
                    parts = data.split(";",1) 
                    seq_len = parts[1].rstrip(".")
                    parts = parts[0].split(" ")
                    record.entry_id = parts[0]
                    record.status = parts[1].strip()
                elif line_code == "AC":
                    # Accession number line: AC   AC_number_1;[ AC_number_2;]...[ AC_number_N;]
                    record.accession_ids.extend(UniprotEntry.SEMICOLON_SEP.split(data.strip(";")))
                elif line_code == "DE":
                    # Description lines
                    match = UniprotEntry.DE_TYPE.match(data)
                    if match:
                        de_name_data = UniprotString(match.group(3))
                        current_name_type = match.group(1)
                        if subname == 1:
                            record.contains_names[-1][current_name_type][match.group(2)] = de_name_data
                        elif subname == 2:
                            record.includes_names[-1][current_name_type][match.group(2)] = de_name_data
                        else:
                            record.record_names[current_name_type][match.group(2)] = de_name_data
                    elif data.strip() == "Contains:":
                        subname = 1
                        record.contains_names.append({ "Rec": dict(), "Alt": dict(), "Sub": dict() })
                    elif data.strip() == "Includes:":
                        subname = 2
                        record.includes_names.append({ "Rec": dict(), "Alt": dict(), "Sub": dict() })
                    elif data.lstrip().startswith("Flags: "):
                        current_name_type = None
                        subname = 3
                        parts = data.split(":")
                        record.name_flags.add(parts[1].strip().rstrip(";"))
                    else:
                        parts = data.split("=",1)
                        de_name_data = UniprotString(parts[1].strip(";"))
                        if subname == 1:
                            record.contains_names[-1][current_name_type][parts[0].strip()] = de_name_data
                        elif subname == 2:
                            record.includes_names[-1][current_name_type][parts[0].strip()] = de_name_data
                        else:
                            record.record_names[current_name_type][parts[0].strip()] = de_name_data
                elif line_code == "GN":
                    # Gene names...
                    # GN   Name=<name>; Synonyms=<name1>[, <name2>...]; OrderedLocusNames=<name1>[, <name2>...];
                    # GN   ORFNames=<name1>[, <name2>...];
                    if data.strip() == "and":
                        record.gene_names.append(UniprotGeneName(gn_data))
                        gn_data = ""
                    else:
                        gn_data += data + " "
                elif line_code == "OS":
                    # Organism species line
                    record.organism = data.rstrip(".")
                elif line_code == "RN":
                    # Reference number line.  Start of a reference record
                    record.references.append(UniprotReference(data))
                    ref_num = record.references[-1].ref_number
                elif line_code == "RP" and ref_num > 0: 
                    # Reference position
                    record.references[-1].position += data.strip()
                elif line_code == "RC" and ref_num > 0: 
                    # Reference comment
                    for comment in UniprotEntry.SEMICOLON_SEP.split(data.rstrip(';')):
                        if "=" in comment:
                            comment_token, comment_data = comment.split("=",1)
                            record.references[-1].comment[comment_token] = comment_data
                        else:
                            record.references[-1].comment[comment_token] += comment_data
                elif line_code == "RX" and ref_num > 0: 
                    # Reference cross reference
                    for item in UniprotEntry.SEMICOLON_SEP.split(data.rstrip(';')):
                        parts = item.split("=",1)
                        record.references[-1].cross_refs[parts[0]] = parts[1]
                elif line_code == "RG" and ref_num > 0: 
                    # Reference group (author affiliation)
                    record.references[-1].group += data.strip().rstrip(';')
                elif line_code == "RA" and ref_num > 0: 
                    # Reference authors
                    for item in UniprotEntry.COMMA_SEP.split(data.rstrip(",").rstrip(';')):
                        record.references[-1].authors.append(item)
                elif line_code == "RT" and ref_num > 0: 
                    # Reference title
                    record.references[-1].title += data.lstrip('"').rstrip(";").rstrip('"')
                elif line_code == "RL" and ref_num > 0: 
                    # Reference location, e.g. journal information
                    record.references[-1].location += data.strip()
                elif line_code == "CC":
                    # Entry comments
                    if data.startswith('-!-'):
                        data = data[3:].strip()
                        (comment_topic, data) = data.split(':',1)
                        record.comments[comment_topic] = data.strip()
                    elif comment_topic is not None:
                        record.comments[comment_topic] += data.strip()
                elif line_code == "DR": 
                    # Database cross reference 
                    # TODO: should better organize the DB refs
                    parts = UniprotEntry.SEMICOLON_SEP.split(data.strip())
                    record.db_refs[parts[0]].append(parts[1:])  
                elif line_code == "PE":
                    # Protein evidence line
                    parts = data.rstrip(";").split(":")
                    record.evidence = ( int(parts[0].strip()), parts[1] )
                elif line_code == "KW":
                    record.keywords.extend(UniprotEntry.SEMICOLON_SEP.split(data.rstrip(";")))
                # elif line_code == "SQ":
                    # Sequence header TODO
            except:
                logging.warning("Unrecognized line format: {}".format(line))
                raise

    ##
    # Method to get a RecName
    def get_recname(self, name_type="Full"):
        return self.record_names["Rec"].get(name_type, None)

    ##
    # Method to get an AltName
    def get_altname(self, name_type="Full"):
        return self.record_names["Alt"].get(name_type, None)

    ##
    # Method to get a SubName
    def get_subname(self, name_type="Full"):
        return self.record_names["Sub"].get(name_type, None)

    ##
    # Method to get a Gene Name (symbol)
    def get_genename(self, index=0):
        return self.gene_names[index]

    ##
    # Method to get the organism name
    def get_organism(self):
        return self.organism

    ##
    # Method to get the entry ID
    def get_entry_id(self):
        return self.entry_id
    
    ##
    # Method to get the entry status
    def get_status(self):
        return self.status

    ##
    # Method to get the accession IDs
    def get_accession_ids(self):
        return self.accession_ids


##
# Class to store reference (R*) information from a Uniprot record
class UniprotReference:

    COMPLEX_RN_PATTERN = re.compile('\[([0-9]+)\] {(.+)}$')
    
    def __init__(self, rn_data):
                
        match = UniprotReference.COMPLEX_RN_PATTERN.match(rn_data)
        if match:
            self.ref_number = int(match.group(1))
            self.ref_tags = [ x.strip() for x in match.group(2).split(",") ] 
        else:
            self.ref_number = int(rn_data.strip().replace('[','').replace(']',''))
            self.ref_tags = list()
        self.position = ""
        self.comment = dict()
        self.cross_refs = dict()
        self.group = ""
        self.authors = list()
        self.title = ""
        self.location = ""

    def get_authors(self):
        return self.authors

    def get_cross_ref(self, db):
        return self.cross_refs.get(db, None)
        
    def get_comment(self, token):
        return self.comment.get(token, None)
    
    def __str__(self):
        return ", ".join(self.authors) + ' "' + self.title + '" ' + self.location        

##
# Class to store gene names (GN) information from Uniprot
class UniprotGeneName:

    def __init__(self, gn_data):
        self.name = None
        self.synonyms = list()
        self.ordered_locus_names = list()
        self.orf_names = list()

        for name_data in UniprotEntry.SEMICOLON_SEP.split(gn_data):
            name_parts = name_data.split("=",1)
            if name_parts[0] == "Name":
                self.name = UniprotString(name_parts[1])
            elif name_parts == "Synonyms":
                self.synonyms = [ UniprotString(x) for x in UniprotEntry.COMMA_SEP.split(name_parts[1].strip()) ]
            elif name_parts == "OrderedLocusNames":
                self.ordered_locus_names = [ UniprotString(x) for x in UniprotEntry.COMMA_SEP.split(name_parts[1].strip()) ]
            elif name_parts == "ORFNames":
                self.orf_names = [ UniprotString(x) for x in UniprotEntry.COMMA_SEP.split(name_parts[1].strip()) ]

class UniprotString:

    EVIDENCE_PATTERN = re.compile('^ *(.+) {(.+)} *$')
    
    def __init__(self, content):
        match = UniprotString.EVIDENCE_PATTERN.match(content)
        if match:
            self.content = match.group(1)
            self.evidence = [ UniprotEvidence(x.strip()) for x in match.group(2).split(",") ]
        else:
            self.content = content
            self.evidence = list()

    def __str__(self):
        return self.content

class UniprotEvidence:

    def __init__(self, evidence_string):
        if "|" in evidence_string:
            parts = evidence_string.split("|", 1)
            self.source = parts[1].split(":",1)
        else:
            self.source = None

        self.eco = evidence_string[4:]
