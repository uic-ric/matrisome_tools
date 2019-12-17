################################################################################
# Script : uniprot.py
# Author : George Chlipala
# Created: Oct 22, 2019
# -- Description ----------------------------------------
# Custom python library for parsing Uniprot entries
# -- Requirements ---------------------------------------
# Python 2.7.6
# urllib3 library
################################################################################
"""
.. module:: uniprot
   :platform: Unix, Windows
   :synopsis: Module to retrieve information from Uniprot (http://www.uniprot.org)

.. moduleauthor:: George Chlipala <gchlip2@uic.edu>

"""

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
    """API class for Uniprot.  This class has static methods/functions for retreiving Uniprot records and converting protein/gene IDs
    """

    # Base URL for Uniprot
    UNIPROT_API_BASE = "https://www.uniprot.org"
    # Number of entries to process at a single time.
    chunk_size = 5000
    logger = logging.getLogger("uniprot.UniprotAPI")

    ##
    # Function to get a set of Uniprot entries
    @staticmethod
    def get(*entries, **kwargs):
        """
        Function to retrieve records from Uniprot
        
        Args:
            *args: Uniprot accession IDs to fetch.

        Keyword Args:
            id_type: ID type.  Default is "ACC+ID", a.k.a UniprotKB AC/ID. Refer to https://www.uniprot.org/help/api_idmapping for a list of valid ID types

        Returns:
            A generator of UniprotEntry objects
        """
        all_entries = list()
        for e in entries:
            if isinstance(e, basestring):
                all_entries.append(e)
            else:
                all_entries.extend(e)
        id_type = kwargs.get("id_type", "ACC+ID")
        # Check if the entry list is greater than the chunk size (5000)
        # If so, perform the request a chunk at a time.
        response_data = ""
        e_len = len(all_entries)
        UniprotAPI.logger.debug("Retrieving Uniprot information for {} entries.".format(e_len))
        if e_len > UniprotAPI.chunk_size:
            UniprotAPI.logger.debug("Entry list greater than chunk size of {}. Will retrieve in batches.".format(UniprotAPI.chunk_size))
            for n in range(0, e_len, UniprotAPI.chunk_size):
                max_n = min(n + UniprotAPI.chunk_size, e_len)
                response_data += UniprotAPI.__get_entries__(all_entries[n:max_n], id_type)
        else:
            response_data += UniprotAPI.__get_entries__(all_entries, id_type)
        
        return UniprotEntry.parse(cStringIO.StringIO(response_data))

    ##
    # Hidden method to actually perform the get request.
    @staticmethod
    def __get_entries__(entries, from_type="ACC+ID", output_format="txt", to_type="ACC"):
        """
        Function to perform API calls to Uniprot.  This function is called by the get() and convert() functions and should not be called directly.
        
        Args:
            entries: a list of entries to fetch
            from_type: ID type. Default is "ACC+ID", a.k.a UniprotKB AC/ID. Refer to https://www.uniprot.org/help/api_idmapping for a list of valid ID types
            output_format: Format of output.  Default is "txt".  Other allowed values are "tab", "xml", "rdf", "fasta" and "gff". 
            to_type: ID type. Default is "ACC", a.k.a UniprotKB AC. Refer to https://www.uniprot.org/help/api_idmapping for a list of valid ID types. Only recognized if output_format is "tab"
    
        Returns:
            A dict mapping input IDs with converted IDs.
        """
        # Set destination ID type to "ACC" unless the output_format is "ACC"
        if output_format != "tab":
            to_type = "ACC"
        UniprotAPI.logger.debug("Performing HTTP request to Uniprot with format={}, from={}, to={} and {} entries.".format(output_format, from_type, to_type, len(entries)))
        http = urllib3.PoolManager(headers={"User-Agent": "urllib3-python"})
        req = http.request("POST", UniprotAPI.UNIPROT_API_BASE + "/uploadlists/", redirect=False, 
            fields={ "file": ("entries.txt", "\n".join(entries) + "\n", "text/plain"),
                "format": output_format, "from": from_type, "to": to_type })
        # The following is a work around for the redirect.
        # The orginal request is performed as a POST, however Uniprot then sends a 302 (redirect)
        # to a URL that should be accessed via GET. 
        # urllib3 tries the redirect as a POST and Uniprot generates a 500 (server error)
        UniprotAPI.logger.debug("HTTP status: {}".format(req.status))
        if req.status == 302:
            location = req.get_redirect_location()
            if not location.startswith("http"):
                location = UniprotAPI.UNIPROT_API_BASE + location
            req = http.request("GET", location)
        return req.data

    ##
    # Function to convert IDs
    @staticmethod
    def convert(id_type, *ids, **kwargs):
        """
        Function to convert protein IDs.
        
        Args:
            id_type: source ID type. Refer to https://www.uniprot.org/help/api_idmapping for a list of valid ID types
            *args: IDs to convert
    
        Keyword Args:
            to_type: Destination ID type. Default is "ACC", a.k.a. Uniprot accession ID.  Refer to https://www.uniprot.org/help/api_idmapping for a list of valid ID types

        Returns:
            A dict mapping input IDs with converted IDs.
        """
        all_ids = list()
        for i in ids:
            if isinstance(i, basestring):
                all_ids.append(i)
            else:
                all_ids.extend(i)

        to_type = kwargs.get("to_type", "ACC")
        id_map = dict()
        e_len = len(all_ids)
        UniprotAPI.logger.debug("Converting {} IDs from {} to Uniprot accession IDs.".format(e_len, id_type))
        if e_len > UniprotAPI.chunk_size:
            UniprotAPI.logger.debug("Entry list greater than chunk size of {}. Will retrieve in batches.".format(UniprotAPI.chunk_size))
            for n in range(0, e_len, UniprotAPI.chunk_size):
                max_n = min(n + UniprotAPI.chunk_size, e_len)
                id_map.update(UniprotAPI.convert(id_type, all_ids[n:max_n], to_type=to_type))
        else:
            data = UniprotAPI.__get_entries__(all_ids, id_type, "tab", to_type)
            # http = urllib3.PoolManager(headers={"User-Agent": "urllib3-python"})
            # req = http.request("POST", UniprotAPI.UNIPROT_API_BASE + "/uploadlists/", encode_multipart=False,
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
    """Class for Uniprot entries (records).

    :ivar references: References/citations for the entry.
    :ivar evidence_level: Evidence level for the existence of the protein. 
        Possible values are 1 (evidence at protein level), 2 (evidence at transcript level), 3 (inferred from homology), 4 (predicted), and 5 (uncertain)
    :ivar keywords: Keywords for the entry.
    """

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
        self.evidence_level = 0
        self.evidence_text = ""
        self.keywords = list()

    ##
    # static method to parse entries from an I/O handle.
    # This method acts as a generator
    @staticmethod
    def parse(handle):
        """Function to parse UniProt entries/records from an I/O handle.  
        Data should be in UniProt KB format, a.k.a UniProt text format.
        Details about this format can be found at https://web.expasy.org/docs/userman.html

        Args:
            handle: I/O handle.  Can be from file or HTTP request.

        Returns:
            A generator of UniprotEntry objects.
        """
    
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
                    for comment in UniprotEntry.SEMICOLON_SEP.split(data.rstrip(";")):
                        if "=" in comment:
                            comment_token, comment_data = comment.split("=",1)
                            record.references[-1].comment[comment_token] = comment_data
                        else:
                            record.references[-1].comment[comment_token] += comment_data
                elif line_code == "RX" and ref_num > 0: 
                    # Reference cross reference
                    for item in UniprotEntry.SEMICOLON_SEP.split(data.rstrip(";")):
                        parts = item.split("=",1)
                        record.references[-1].cross_refs[parts[0]] = parts[1]
                elif line_code == "RG" and ref_num > 0: 
                    # Reference group (author affiliation)
                    record.references[-1].group += data.strip().rstrip(";")
                elif line_code == "RA" and ref_num > 0: 
                    # Reference authors
                    for item in UniprotEntry.COMMA_SEP.split(data.rstrip(",").rstrip(";")):
                        record.references[-1].authors.append(item)
                elif line_code == "RT" and ref_num > 0: 
                    # Reference title
                    record.references[-1].title += data.lstrip(""").rstrip(";").rstrip(""")
                elif line_code == "RL" and ref_num > 0: 
                    # Reference location, e.g. journal information
                    record.references[-1].location += data.strip()
                elif line_code == "CC":
                    # Entry comments
                    if data.startswith("-!-"):
                        data = data[3:].strip()
                        (comment_topic, data) = data.split(":",1)
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
                    record.evidence_level = int(parts[0].strip())
                    record.evidence_text = parts[1]
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
        """Get the recommended name (RecName) of the entry.

        Args:
            name_type: Name type to return. Valid name types are "Full", "Short", or "EC". Default is "Full"
        
        Returns:
            Name as an UniprotString
        """
        return self.record_names["Rec"].get(name_type, None)

    ##
    # Method to get an AltName
    def get_altname(self, name_type="Full"):
        """Get alternate names (AltName) of the entry.

        Args:
            name_type: Name type to return. Valid name types are "Full", "Short", or "EC". Default is "Full"
        
        Returns:
            Name as an UniprotString
        """
        return self.record_names["Alt"].get(name_type, None)

    ##
    # Method to get a SubName
    def get_subname(self, name_type="Full"):
        """Get submitted name (SubName) of the entry.

        Args:
            name_type: Name type to return. Valid name types are "Full", "Short", or "EC". Default is "Full"
        
        Returns:
            Name as an UniprotString
        """
        return self.record_names["Sub"].get(name_type, None)

    ##
    # Method to get a Gene Name (symbol)
    def get_genename(self, index=0):
        """Get gene name of the entry.

        Args:
            index: Index of the name to return. Default is 0, i.e. the first name in the entry.
        
        Returns:
            Name as an UniprotGeneName object.
        """
        return self.gene_names[index]

    ##
    # Method to get the organism name
    def get_organism(self):
        """Get the source organisms of the entry.

        Returns:
            Organism as a string.
        """
        return self.organism

    ##
    # Method to get the entry ID
    def get_entry_id(self):
        """Get the entry ID (UniProtKB ID)

        Returns:
            Entry ID as a string
        """
        return self.entry_id
    
    ##
    # Method to get the entry status
    def get_status(self):
        """Get the entry status.

        Returns:
            Entry status, i.e. "Reviewed" (SwissProt) or "Unreviewed" (TrEMBL)
        """
        return self.status

    ##
    # Method to get the accession IDs
    def get_accession_ids(self):
        """Get the UniProt accession IDs associated with the entry.

        Returns:
            list of Accession IDs
        """
        return self.accession_ids

    ##
    # Method to get comments
    def get_comment(self, topic):
        """Get a comments for the entry.

        Args:
            topic: UniProt topic. Refer to https://web.expasy.org/docs/userman.html#CC_line for a list of valid topics

        Returns:
            Comment for the selected topic. If no comment exists for the topic will return None.
        """
        return self.comments.get(topic, None)

##
# Class to store reference (R*) information from a Uniprot record
class UniprotReference:
    """Class to contain UniProt reference/citation information
    
    :ivar position: This describes the extent of the work relevant to the entry carried out by the authors.
    :ivar group: The consortium name associated with the reference/citation
    :ivar title: Title of the paper or other work.
    :ivar location: Publication location of the reference/citation, e.g. journal name with volume and page information.
    """

    COMPLEX_RN_PATTERN = re.compile('\[([0-9]+)\] {(.+)}$')
    
    def __init__(self, rn_data):
                
        match = UniprotReference.COMPLEX_RN_PATTERN.match(rn_data)
        if match:
            self.ref_number = int(match.group(1))
            self.ref_tags = [ x.strip() for x in match.group(2).split(",") ] 
        else:
            self.ref_number = int(rn_data.strip().replace("[","").replace("]",""))
            self.ref_tags = list()
        self.position = ""
        self.comment = dict()
        self.cross_refs = dict()
        self.group = ""
        self.authors = list()
        self.title = ""
        self.location = ""

    def get_authors(self):
        """Get the authors for the reference.

        Returns:
            list of authors
        """
        return self.authors

    def get_cross_ref(self, db):
        """Get the reference IDs for the reference.

        Args:
            db: Reference database to select, e.g. PubMed, DOI.

        Returns:
            list of IDs for the selected database.
        """
        return self.cross_refs.get(db, None)
        
    def get_comment(self, token):
        """Get comments for the reference

        Args:
            token: Valid tokens are STRAIN, PLASMID, TRANSPOSON, and TISSUE 
        
        Returns:
            Comments for the specified token
        """
        return self.comment.get(token, None)
    
    def __str__(self):
        return ", ".join(self.authors) + " \"" + self.title + "\" " + self.location        

##
# Class to store gene names (GN) information from Uniprot
class UniprotGeneName:
    """Class to contain UniProt gene names. All names are stored as UniprotString objects.
    
    :ivar name: Primary gene name.
    :ivar synonyms: Synonyms
    :ivar ordered_locus_names: Ordered locus names, a.k.a. OLN, ORF numbers, CDS numbers or Gene numbers.
    :ivar orf_names: Open reading frame (ORF) names, a.k.a. sequencing names or contig names or temporary ORF names
    """

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


# Class to store Uniprot strings with evidence information
class UniprotString(str):
    """Class to store a string from UniProt that may contain additional evidence information.

    :ivar evidence: List of evidence annotations
    """

    EVIDENCE_PATTERN = re.compile('^ *(.+) {(.+)} *$')
    
    def __new__(cls, content):
        match = UniprotString.EVIDENCE_PATTERN.match(content)
        if match:
            content = match.group(1)
            evidence = [ UniprotEvidence(x.strip()) for x in match.group(2).split(",") ]
        else:
            evidence = list()
        self = super(UniprotString, cls).__new__(cls, content)
        self.evidence = evidence
        return self

    def with_evidence(self):
        """Method to return the string with evidence annotations

        Returns:
            String with evidence annotations
        """
        return str(self) + " {" + ", ".join([str(e) for e in self.evidence]) + "}"

##
# Class to store Uniprot evidence information
class UniprotEvidence:
    """Class to store UniProt evidence attributions/annotations.

    :ivar eco: Evidence Codes Ontology ID for the attribution/annotation.
    :ivar source_database: Source database of the evidence attribution/annotation.
    :ivar source_id: Source ID of the evidence attribution/annotation.
    """

    def __init__(self, evidence_string):
        self.source_database = None
        if "|" in evidence_string:
            parts = evidence_string.split("|", 1)
            if ":" in parts[1]:
                self.source_database, self.source_id = parts[1].split(":",1)
            else:
                self.source_id = parts[1]
            evidence_string = parts[0]
        else:
            self.source_id = None

        self.eco = evidence_string[4:]

    def __str__(self):
        if self.source_database is None:
            return "ECO:" + self.eco
        else:
            return "ECO:" + self.eco + "|" + self.source_database + ":" + self.source_id
