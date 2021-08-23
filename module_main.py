from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
from Bio import Align
from Bio.Blast import Record as blastRecord
from Bio.Blast import NCBIXML
from dataclasses import dataclass
from math import floor
from Bio import SeqFeature
from Bio import SeqRecord
from module_applications import hmmscan
import re
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from os import listdir
import urllib
import os
import pickle


"""
import importlib
import sys
"""

# Display Settings ------------------------------------------------------------------------

Entrez.email = 'ahatoum@ua.edu'
pd.set_option('display.expand_frame_repr', False)  # Prints all columns
pd.set_option("display.max_rows", 6000, "display.max_columns", 100, "max_colwidth", 40)
genome_path = 'genomes/gbk_format/'
neighb_path = "/Users/baslan/PycharmProjects/NHI_proximity_To_Defences/neighbours/" # THIS HAS TO BE FULL PATH! (to be passed to hmmer)

# Utility Functions -----------------------------------------------------------------------

def read_file(file):
    text_file = open(file, "r")
    my_text = text_file.read()
    text_file.close()
    return my_text

def write_file(path, text):
    with open(path, 'w') as outfile:
        outfile.write(text)
        outfile.close()

def parse_neighbourhood_filename(filename):

    m = re.search(r'^(.+)_[nopb]__pfams_(.+).txt$', filename)
    genome_accession = m.group(1)
    locus_tag = m.group(2)
    return genome_accession, locus_tag

# This was used once to correct the genbank file extensions.
def rename_files():
    import os

    path = "genomes/gbk_format/"

    for filename in listdir(path):
        if filename.endswith('.txt'):
            filename_new = filename.split('.txt')[0] + '.gbk'
            print()
            print(path+filename)
            print(path+filename_new)

            # os.rename(path+filename, path+filename_new)

def read_fasta():
    # todo: make this useful.
    from Bio import SeqIO
    for seq_record in SeqIO.parse("proteins/WP_033072565.1.txt", "fasta"):

        print(seq_record.id)
        print(repr(seq_record.seq))
        print(seq_record.seq)
        print(len(seq_record))

def get_nc_nz_number(accession):

    """
    If you provide a non NZ/NC acession number, you get the NC/NZ equivalent that are in the same assembly report.

    Usage:

    ids = get_nc_nz_number('CP013020.1')
    print(ids)

    """

    handle = Entrez.esearch(db="assembly", term=accession)  # finds the assebmly with that accesion in it.
    record = Entrez.read(handle)
    print()
    print(record)
    refseq_ids = list()

    for id in record['IdList']:

        esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
        esummary_record = Entrez.read(esummary_handle)

        accession_id = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['Synonym']['Genbank']  # ['AssemblyAccession']

        print()
        print(accession_id)
        print()
        print(esummary_record)

        rec = Entrez.read(Entrez.esearch(db="nucleotide", term=accession_id, idtype="acc"))  # This is the unexpected magic line
        refseq_id_list = rec['IdList']
        print()
        print(rec)

        for rfid in refseq_id_list:
            refseq_ids.append(rfid)

    return refseq_ids

# code source: https://www.biostars.org/p/141581/
def entrez_efetch_experiment():

    handle = Entrez.efetch(db="assembly", id='33028', rettype='gb', retmode='text')
    text = handle.read()
    print(text)

def entrez_esearch_experiment():

    # {'Count': '1', 'RetMax': '1', 'RetStart': '0', 'IdList': ['33028'], 'TranslationSet': [], 'TranslationStack': [{'Term': 'CP000029.1[All Fields]', 'Field': 'All Fields', 'Count': '1', 'Explode': 'N'}, 'GROUP'], 'QueryTranslation': 'CP000029.1[All Fields]'}

    handle = Entrez.esearch(db="assembly", term="CP000029.1", retmax=814, idtype="acc")
    record = Entrez.read(handle)
    handle.close()
    print(record['IdList'][0])

def entrez_summary_experiment():

    # esummary_handle = Entrez.esummary(db="nucleotide", id='559341', report="full")
    # esummary_record = Entrez.read(esummary_handle)
    # print(esummary_record)

    esummary_handle = Entrez.esummary(db="assembly", id='33028', report="full")
    esummary_record = Entrez.read(esummary_handle)
    esummary_handle.close()
    print(esummary_record)

# Explore
def get_assembly_summary(id):
    """
    Get esummary for an entrez id
    code source: https://dmnfarrell.github.io/bioinformatics/assemblies-genbank-python
    """
    from Bio import Entrez
    esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
    esummary_record = Entrez.read(esummary_handle)
    return esummary_record

# Eexplore
def get_assemblies(term, download=True, path='assemblies'):
    """Download genbank assemblies for a given search term.
    Args:
        term: search term, usually organism name
        download: whether to download the results
        path: folder to save to

    source: https://dmnfarrell.github.io/bioinformatics/assemblies-genbank-python
    """

    from Bio import Entrez
    #provide your own mail here
    Entrez.email = "A.N.Other@example.com"
    handle = Entrez.esearch(db="assembly", term=term, retmax='200')
    record = Entrez.read(handle)
    ids = record['IdList']
    print (f'found {len(ids)} ids')
    links = []
    for id in ids:
        #get summary
        summary = get_assembly_summary(id)
        #get ftp link
        url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
        if url == '':
            continue
        label = os.path.basename(url)
        #get the fasta link - change this to get other formats
        link = os.path.join(url,label+'_genomic.fna.gz')
        print (link)
        links.append(link)
        if download == True:
            #download link
            urllib.request.urlretrieve(link, f'{label}.fna.gz')
    return links

# Downloaders -------------------------------------------------------------------------

def download_protein_from_entrez(fasta_path, id):

    try:
        with Entrez.efetch(db="protein", rettype="fasta", retmode="text", id=id) as handle:
            seq_record = SeqIO.read(handle, "fasta")
            SeqIO.write(seq_record, fasta_path + id + '.txt', "fasta")
            print(id, ' successful...')
            return seq_record

    except:

        print(id, ' encountered a problem.')
        return 0

def download_genome(path, accession_numbers, rettype, retmode):

    """
    For genbank full, use rettype = 'gbwithparts'
    For regular, use rettype = 'gb'
    retmode = 'xml' or 'text'

    Pay attention to _gb file splitter if you change extentions.

    """

    filenames = [filename.split('_gb')[0] for filename in listdir(path)]

    n=0
    L = len(accession_numbers)

    for ID in accession_numbers:

        n+=1

        if ID in filenames:
            print('{}/{}: {} already downloaded...'.format(n, L, ID))
            continue

        time.sleep(0.2)
        t0 = time.perf_counter()

        try:
            handle = Entrez.efetch(db="nuccore", id=ID, rettype=rettype, retmode=retmode)
            text = handle.read()

        except:
            text = "None"
        t1 = time.perf_counter()

        if retmode == 'text':
            filename = '{}_{}.gbk'.format(ID, rettype)
            fopen_par =  "w"
        elif retmode == 'xml':
            filename = '{}_{}.xml'.format(ID, rettype)
            fopen_par =  "wb"
        else:
            raise ValueError('retmode not recognized')

        print('{}/{}: {} -- downloaded in {} seconds'.format(n, L, filename, round(t1 - t0, 2)))
        file = open(path + filename, fopen_par)
        file.write(text)
        file.close()

def download_genomes_from_blast_results(hit_list_file):

    PTBR = ReadBlastOutputXML(hit_list_file)  # 'all_features'
    accession_numbers = [p.accession for p in PTBR.get_blast_hit_object_list()]
    download_genome(genome_path, accession_numbers, rettype='gbwithparts', retmode='text')

# TOOLKIT ------------------------------------------------------------------------------

# This is a collection of tools from o'th attempt AND Mostly obsolete!
# It was designed to download and stitch together protein sequences from a given list of protein accession numbers.
class UniqueHomologs:

    input_wp_text = 'input files/wp_numbers.txt'
    output_wp_csv = 'output files/wp_numbers_counted.csv'
    output_stitched = 'output files/fastas_stitched.txt'
    fasta_path = 'output files/Downloaded Fasta Files/'

    def download_and_stitch_fasta_files(self):

        ids = self.get_ids_from_counted(self.output_wp_csv)

        for id in ids:
            download_protein_from_entrez(self.fasta_path, id)
            time.sleep(0.2)

        self.stitch_fasta_files(ids, self.fasta_path, self.output_stitched)

    def count_wps(self):

        lines = read_file(self.input_wp_text)
        lines = lines.splitlines()

        wp_dict = dict()
        for line in lines:

            wp_number = (re.search(r'\[protein_id=(WP_[0-9]+\.[0-9]) *\]', line)).group(1)

            if wp_number not in wp_dict.keys():
                wp_dict[wp_number] = 1

            else:
                wp_dict[wp_number] = wp_dict[wp_number] + 1

        rows = []
        for wp_number in wp_dict.keys():
            row = [wp_number, wp_dict[wp_number]]
            rows.append(row)

        df = pd.DataFrame(rows, columns=['WP_Number', 'Count'])
        df.sort_values(by=['Count'], inplace=True, ascending=False)
        df.set_index('WP_Number', inplace=True)
        print(df)
        df.to_csv(self.output_wp_csv, sep=',')

    @staticmethod
    def get_ids_from_counted(path):

        df = pd.read_csv(path)
        df.set_index('WP_Number', inplace=True)
        ids = list(df.index)
        return ids

    @staticmethod
    def stitch_fasta_files(ids, fasta_path, output_path):

        """Reads individual fasta files and combines them into a single one. """

        long_string = ''

        for idnum in ids:

            text = read_file(fasta_path + idnum + '.txt')
            text = re.sub(r'\n+', '', text)
            text = re.sub(r']', ']\n', text)
            long_string += text + '\n\n'

        print(long_string)

        write_file(output_path, long_string)

class FractionAnalyzer:

    """

    This reads in a multiple sequence alignment file in 'clustal' format (e.g., MUSCLE)
    and compares letter by letter each alignment to the reference as it is defined.
    Comparison means it increments the corresponding match by one if the subject alignment
    letter is the same as the reference.
    Then, at the end, it chucks all counters that correspond to a '-' location in the reference
    sequence.

    set0 = ["sample_treefiles/04_26_21_muscle-I20210427-040201-0008-84212899-p1m.clw", 'WP_002489608.1']
    set1 = ["sample_treefiles/test_alignment.txt", 'WP_095669550.1']
    fa = mba.FractionAnalyzer(set0[0], reference_id=set0[1])
    fa.run_fraction_analysis()
    """

    def __init__(self, path, frmt = 'clustal', reference_id = 'WP_002489608.1'):

        # reference_id is NHI accession number.
        alignment: Align.MultipleSeqAlignment = AlignIO.read(path, frmt)
        a: SeqRecord.SeqRecord
        self.alignment_dict = dict()
        self.reference_id = reference_id

        for a in alignment:
            if a.id != reference_id:
                self.alignment_dict[a.id] = a.seq
            else:
                self.reference_seq = a.seq

        n1 = len(alignment)
        n2 = len(self.alignment_dict.keys())
        self.nhomologs = n2

        if  n1 - 1 != n2:

            for a in alignment: print(a.id)
            print()
            print()
            for key in self.alignment_dict.keys(): print(key)
            raise ValueError('Alignment set is not unique')

        print(self.reference_seq)
        seq = str(self.reference_seq).replace('-','')
        print(seq)

    def run_fraction_analysis(self):

        def compare(ref, s, c):

            for string in s:

                for i, ref_aa in enumerate(ref):
                    if string[i] == ref_aa:
                        c[i] += 1

            return c

        # CALCULATE ----------------------------------

        s1 = self.reference_seq
        s = [self.alignment_dict[key] for key in self.alignment_dict.keys()]

        c = np.zeros(len(s1))
        c1 = compare(self.reference_seq, s, c)

        s1_list = [l for l in s1]
        rows = list()
        for i, l in enumerate(s1_list):
            if l != '-':
                rows.append([l, c1[i]])

        position = list()
        aa = list()
        count = list()

        pandas_rows = list()
        for j, row in enumerate(rows):

            position.append(j+1)
            aa.append(row[0])
            count.append(round(100 * row[1]/self.nhomologs, 2))
            print(j+1, row[0], row[1])
            pandas_rows.append([j+1, row[0], row[1], round(100 * row[1]/self.nhomologs, 2)])

        df = pd.DataFrame(pandas_rows, columns = ['Position', 'AA', 'Count', 'Percent'])
        df.set_index('Position', inplace=True)
        df.to_csv('fraction_analysis.csv', sep =',')

        # PLOT ----------------------------------

        plt.bar(position, count, width=0.3)
        plt.xlabel('Position')
        plt.ylabel('Count (%)')
        plt.title('NHI Homologs Fractions')
        plt.show()

@dataclass
class BlastHit:

    ridnum: str = ''
    accession: str = ''
    description: str = ''
    hit_range: str = ''
    e_value: str = ''
    total_score: str = ''
    query_cover: str = ''
    per_ident: str = ''
    query: str = ''
    match: str = ''
    sbjct: str = ''

    def __str__(self):

        s = ['RID#: {}'.format(self.ridnum),
             'Accession: {}'.format(self.accession),
             'Description: {}'.format(self.description),
             'Hit Range: {}'.format(self.hit_range),
             'Expect: {}'.format(self.e_value),
             'Total Score: {}'.format(self.total_score),
             'Query Coverage: {}'.format(self.query_cover),
             'Per Ident: {}'.format(self.per_ident)]

        return '\n' + '\n'.join(s)

@dataclass
class gbkProtein:

    parse_case: str = ''
    location: str = ''
    locus_tag: str = ''
    refseq: str = ''
    product: str = ''
    protein_id: str = ''
    translation: str = ''
    note: str = ''
    pseudo: str = ''

    def __str__(self):

        s = ['CASE: {}'.format(self.parse_case),
             'locn: {}'.format(str(self.location)),
             'ltag: {}'.format(self.locus_tag),
             'rseq: {}'.format(self.refseq),
             'prdt: {}'.format(self.product),
             'prid: {}'.format(self.protein_id),
             'trns: {}...'.format(self.translation[0:20]),
             'note: {}'.format(self.note)]

        return '\n' + '\n'.join(s)

class gbkGenome:

    def __init__(self, accession_name):

        global genome_path
        file = genome_path + '{}_gbwithparts.gbk'.format(accession_name)

        records = SeqIO.parse(file, "gb")  # it is an iterator
        self.record = next(records)
        self.accession_name = accession_name

    def __str__(self):

        """
        Data Structure Info:

        type(self.record) is -> <class 'Bio.SeqRecord.SeqRecord'>
        self.record.annotations.keys() are:
        dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'keywords','source', 'organism', 'taxonomy', 'references', 'comment', 'structured_comment'])

        self.record.id  # self.record.annotations['accessions'][0]

        """

        a = [
             'Accession: {}'.format(self.record.id),
             'Description: {}'.format(self.record.description),
             'Taxonomy: {}'.format('/'.join(self.record.annotations['taxonomy'])),
             'Organism: {}'.format(self.record.annotations['organism']),
             'Number of Features: {}'.format(len(self.record.features))
             ]

        return '\n'.join(a)

    @staticmethod
    def convert_feature_to_custom_protein_class(feature):

        Protein = gbkProtein()

        if feature is None:
            return Protein

        Protein.location = str(feature.location)
        quals = feature.qualifiers
        qual_keys = quals.keys()

        # keys(['locus_tag', 'inference', 'note', 'codon_start', 'transl_table', 'product', 'protein_id', 'translation'])

        if 'locus_tag' in qual_keys:
            Protein.locus_tag = quals['locus_tag'][0]

        if 'inference' in qual_keys:

            m = re.search(r'RefSeq:([A-Z]+_[0-9]+\.[0-9])', quals['inference'][0])
            if m is not None:
                Protein.refseq = m.group(1)
            else:

                Protein.refseq = quals['inference'][0]

        if 'product' in qual_keys:
            Protein.product = quals['product'][0]

        if 'protein_id' in qual_keys:
            Protein.protein_id = quals['protein_id'][0]

        if 'translation' in qual_keys:
            Protein.translation = quals['translation'][0]

        if 'note' in qual_keys:
            Protein.note = quals['note'][0]

        if 'pseudo' in qual_keys:
            Protein.pseudo = '<True>'

        return Protein

    def get_coincident_features_from_hit_range(self, range_string, mode='CDS'):

        """
        range_string looks something like c23456-23490, or 1700-19050
        Left leg (ll) is the smaller boundary of this range, and right leg (rl) is the other.
        It retrieves from genome features both where the left and right locations are contained within.
        You can only search CDS or you can search all features using the mode.

        If the search cannot find at least one feature containing one feature, then it shrinks the hit range by 5% from left and right
        and tries one more time. (if both x0, y0 fell on both sides of a feature in the intergenic regions.)

        --------

        # Problematic One, should try the shrinkage.
        # Should find: locus_tag = H7R81_RS18575, protein_id = WP_182969753.1
        # mba.cmnd_single_hit_extraction('NZ_AP022290.1', '4051988-4053541')

        # Should find: locus_tag = CKU_RS00150, protein_id = WP_000632676.1
        # mba.cmnd_single_hit_extraction('NZ_CP014791.1', '38333-39748')


        # One Legged test_case
        # Should find: locus_tag = N506_RS02030, protein_id = WP_095669550.1
        # mba.cmnd_single_hit_extraction('NZ_CP006803.1', '405152-406207')

        # Guaranteed No find case
        # mba.cmnd_single_hit_extraction('NZ_CP014791.1', '1e9-2e9')

        """

        def parse_range_string(range_string):

            m = re.search(r'(c)?([0-9]+)-([0-9]+)', range_string)

            if m is None:
                raise ValueError('PROBLEM PARSING GIVEN HIT RANGE!')
            else:
                is_comp = True if m.group(1) == 'c' else False  # print(is_comp, x0, y0)
                n1 = int(m.group(2))
                n2 = int(m.group(3))
                x0 = min(n1, n2) + 1
                y0 = max(n1, n2) - 1


            return x0, y0

        def calculate_left_percent_coverage(x, feature):
            coverage = feature.location.end - x
            full = feature.location.end - feature.location.start
            percent_coverage = round(100 * coverage / full, 2)
            return percent_coverage

        def calculate_right_percent_coverage(y, feature):
            coverage = y - feature.location.start
            full = feature.location.end - feature.location.start
            percent_coverage = round(100 * coverage / full, 2)
            return percent_coverage

        def find_legs(x0, y0, features):

            ll = None
            rl = None
            lc = None
            rc = None

            for feature in features:
                if x0 in feature:
                    ll = feature
                    lc = calculate_left_percent_coverage(x0, feature)
                if y0 in feature:
                    rl = feature
                    rc = calculate_right_percent_coverage(y0, feature)

            return ll, rl, lc, rc

        x0, y0 = parse_range_string(range_string)

        if mode == 'CDS':

            ll, rl, lc, rc = find_legs(x0, y0, self.get_cds_list())

        elif mode == 'all_features':

            ll, rl, lc, rc = find_legs(x0, y0, self.record.features)

        else:
            raise ValueError('Mode is invalid!')

        if ll is None and rl is None:

            print('Trying in range shrink mode...')

            hit_length = y0-x0
            x0 += int(hit_length * 0.05)
            y0 -= int(hit_length * 0.05)

            ll, rl, lc, rc  = find_legs(x0, y0, self.record.features)

        if ll is not None and ll.type == 'source':
            ll, lc = None, None
        if rl is not None and rl.type == 'source':
            rl, rc = None, None

        return ll, rl, lc, rc

    def find_protein_by_id(self, protein_id):

        features = list()
        for feature in self.get_cds_list():
            quals = feature.qualifiers
            qual_keys = quals.keys()
            if 'protein_id' in qual_keys:
                if protein_id == quals['protein_id'][0]:
                    features.append(feature)

        return features

    def get_neighbours_using_locus_tag(self, path, locus_tag, window_size=10):

        """
        Given a locus tag, this function gets the neighbours within the windowsize that are (i) cds, (ii) have a translation.
        Then it exports a fasta file of these proteins preserving the order and runs an hmmscan on it.

        If circle_indicator is _n_, it indicates in the filename that some neigbours were acquired with negative indices.
        If circle_indicator is _p_, it indicates in the filename that some neigbours were acquired with indices greater than len.
        If circle_indicator is _o_, it indicates in the filename that all neighbours were acquired with in-range indices.

        The indices referred to above are the indices of the coding_list.
        """

        def get_neighbour_indices(a, i, window_size):

            indices = list()
            imax = len(a) - 1

            lower_bound = i - window_size
            upper_bound = i + window_size

            if lower_bound < 0 and upper_bound > imax:
                circle_indicator = '_b_'

            elif lower_bound < 0:
                circle_indicator = '_n_'

            elif upper_bound > imax:
                circle_indicator = '_p_'

            else:
                circle_indicator = '_o_'

            if 2 * window_size + 1 > imax + 1:
                return list(range(0, imax + 1)), circle_indicator

            for m in range(lower_bound, upper_bound + 1):

                if m > imax:
                    indices.append(m - imax - 1)

                else:
                    indices.append(m)

            return indices, circle_indicator

        origin_index = None
        string = ''

        # Get coding_list features -------------------------------------------

        coding_list = list()
        for feature in self.record.features:
            if feature.type == "CDS" and 'translation' in feature.qualifiers.keys():
                coding_list.append(feature)

        # Determine the coding_list index of the subject locus_tag -----------

        for feature in coding_list:
            if 'locus_tag' in feature.qualifiers.keys():
                if locus_tag == feature.qualifiers['locus_tag'][0]:
                    origin_index = coding_list.index(feature)

        # Determine neighbour_indices ----------------------------------------

        neighbour_indices, circle_indicator = get_neighbour_indices(coding_list, origin_index, window_size)

        # Get Neighbours -----------------------------------------------------

        for m in neighbour_indices:

            ltg = coding_list[m].qualifiers['locus_tag'][0]
            pid = coding_list[m].qualifiers['protein_id'][0]
            prd = coding_list[m].qualifiers['product'][0]
            trl = coding_list[m].qualifiers['translation'][0]

            s = '>' + ltg + ' ' + pid + ' ' + prd + '\n' + trl + '\n'
            string += s

        input_name  = self.accession_name + circle_indicator + '_ngbrs_' + locus_tag + '.fa'
        output_name = self.accession_name + circle_indicator + '_pfams_' + locus_tag + '.txt'

        write_file(path + input_name, string)

        return input_name, output_name

    def get_cds_list(self):

        return [feature for feature in self.record.features if feature.type == "CDS"]

    def print_cds_features(self):

        cds_feature_list = self.get_cds_list()
        if any(cds_feature_list):
            print('\n'.join([str(cds) for cds in cds_feature_list]))
        else:
            print('No CDS features found!')

    def print_features(self):

        print()
        print('FEATURES:')
        print('------------------------------------------------')
        for feature in self.record.features:
            print(feature)
            print('\n\n')

class ReadBlastOutputXML:

    """
    The toolkit below is complete:

    The main tool is 'BlastResultsProcessor'.
    It uses the parsed hit-list text from tblastn and extracts the exact protein hits from downloaded genomes (in the genomes folder).
    Then, it creates a dataframe summarizing the results, as well as, creating a fasta-format text file (stiched) for constructing the
    phylogenetic tree.

    Use 'parse_hit_list' method to parse a tblastn hits text file.
    IT DOES NEEED the BlastHit dataclass at the top (which is different than the xmlProtein and gbkProtein).
    Use 'print_hit_list' to print what 'parse_hit_list' returns.

    The method 'download_all_genomes_from_list' uses 'download_genome' to iteratively
    download the genomes extracted by the 'BlastResultsProcessor' from the hitlist output textfile.

    ---

    # main Genome parser from text (gbk)
    # print(a.convert_feature_to_custom_protein_class(a.cds_list[2]))

    """

    all_rec_attributes = ['alignments',
                          'application',
                          'blast_cutoff',
                          'database',
                          'database_length',
                          'database_letters',
                          'database_name',
                          'database_sequences',
                          'date',
                          'descriptions',
                          'dropoff_1st_pass',
                          'effective_database_length',
                          'effective_hsp_length',
                          'effective_query_length',
                          'effective_search_space',
                          'effective_search_space_used',
                          'expect',
                          'filter',
                          'frameshift',
                          'gap_penalties',
                          'gap_trigger',
                          'gap_x_dropoff',
                          'gap_x_dropoff_final',
                          'gapped',
                          'hsps_gapped',
                          'hsps_no_gap',
                          'hsps_prelim_gapped',
                          'hsps_prelim_gapped_attemped',
                          'ka_params',
                          'ka_params_gap',
                          'matrix',
                          'multiple_alignment',
                          'num_good_extends',
                          'num_hits',
                          'num_letters_in_database',
                          'num_seqs_better_e',
                          'num_sequences',
                          'num_sequences_in_database',
                          'posted_date',
                          'query',
                          'query_id',
                          'query_length',
                          'query_letters',
                          'reference',
                          'sc_match',
                          'sc_mismatch',
                          'threshold',
                          'version',
                          'window_size']

    def __init__(self, path):

        self.blast_record: blastRecord = NCBIXML.read(open(path))  # use .parse if multiple results in it to iterate.

        m = re.search(r'/([A-Z0-9]+)-Alignment.xml', path)
        self.ridnum = m.group(1)

    def get_blast_hit_object_list(self):

        blast_hit_object_list = list()
        for alignment in self.blast_record.alignments:

            split_alignment_title = alignment.title.split('|')
            accession_name = split_alignment_title[3]
            description = split_alignment_title[4]

            for hsp in alignment.hsps:

                if hsp.frame[1] < 0:
                    hit_range = 'c{}-{}'.format(hsp.sbjct_end, hsp.sbjct_start)
                else:
                    hit_range = '{}-{}'.format(hsp.sbjct_start, hsp.sbjct_end)

                score = hsp.bits
                query_cover = self.calculate_query_coverage(hsp.query_start,hsp.query_end, self.blast_record.query_length)
                per_ident = self.calculate_percent_identity(hsp.query_start, hsp.query_end, hsp.identities)

                bh = BlastHit()
                bh.ridnum = self.ridnum
                bh.accession = accession_name
                bh.hit_range = hit_range
                bh.description = description
                bh.total_score = score
                bh.query_cover = query_cover
                bh.e_value = hsp.expect
                bh.per_ident = per_ident
                bh.query = hsp.query
                bh.match = hsp.match
                bh.sbjct = hsp.sbjct

                blast_hit_object_list.append(bh)

        return blast_hit_object_list

    def detailed_report(self):

        hsp: blastRecord.HSP
        alignment: blastRecord.Alignment

        print('\n######################################################################################\n')
        print('Query ID: {}'.format(self.blast_record.query_id))
        print('Application: {}'.format(self.blast_record.application))
        print('Query Effective Len: {}'.format(self.blast_record.effective_query_length))
        print('Query Len: {}'.format(self.blast_record.query_length))
        print('Query: {}'.format(self.blast_record.query))
        print('Num Hits: {}'.format(self.blast_record.num_hits))
        print('Matrix: {}'.format(self.blast_record.matrix))
        print('# of Sequences in Database: {}'.format(self.blast_record.num_sequences_in_database))
        print('Gap Penalties: {}'.format(self.blast_record.gap_penalties))
        print('Filter: {}'.format(self.blast_record.filter))
        print('KA Parameters: {}'.format(self.blast_record.ka_params))
        print('Version: {}'.format(self.blast_record.version))
        print('\n######################################################################################\n')

        print(dir(self.blast_record))

        for alignment in self.blast_record.alignments:  # record.alignments is a list

            print('\n')
            print("Alignment --------------------------------------------------------------------")
            print("hit_def: ", alignment.hit_def)
            print("hit_id: ", alignment.hit_id)
            print("sequence(title): ", alignment.title)
            print("length: ", alignment.length)

            for i, hsp in enumerate(alignment.hsps):
                print()
                print("\tHSP #{}".format(i))
                print("\te value: ", hsp.expect)
                print("\tscore: ", hsp.score)
                print("\tidentities :", hsp.identities)
                print("\tpositives :", hsp.positives)
                print("\tnum_alignments :", hsp.num_alignments)
                print("\talign_length :", hsp.align_length)
                print("\tgaps :", hsp.gaps)
                print("\tframe :",
                      hsp.frame)  # todo: If y in frame = (x,y) tuple is negative it is complement range (vice versa)
                print("\tbits :", hsp.bits)
                print("\tQ_start:", hsp.query_start)
                print("\tQ_end:", hsp.query_end)
                print("\tS_start:",
                      hsp.sbjct_start)  # start is always < end. complement is indicated in the second number in frame (negative is complement)
                print("\tS_end:", hsp.sbjct_end)
                print("\tstrand:", hsp.strand)
                print("\tQ:", hsp.query[0:75] + "...")
                print("\tM:", hsp.match[0:75] + "...")
                print("\tS:", hsp.sbjct[0:75] + "...")
                print('\t### DERIVED ###')
                print('\tQuery Coverage: ', self.calculate_query_coverage(hsp.query_start,hsp.query_end, self.blast_record.query_length))
                print('\tPer Ident: ', self.calculate_percent_identity(hsp.query_start, hsp.query_end, hsp.identities))

    def xml_all_header_info(self):
        attributes = [a for a in dir(self.blast_record) if not a.startswith('__')]
        attributes.remove('alignments')
        attributes.remove('descriptions')
        for att in attributes:
            print(att, ': ', getattr(self.blast_record, att))

        print('\n\n')
        for description in self.blast_record.descriptions:
            print(description)

    @staticmethod
    def calculate_query_coverage(query_start, query_end, query_length):

        return '{}%'.format(floor(100 * (query_end-query_start+1)/query_length))

    @staticmethod
    def calculate_percent_identity(query_start, query_end, identities):

        return '{}%'.format(round(100 * identities/(query_end-query_start+1), 3))

class NeighbourhoodFileAnalyzer:

    """

    USAGE:
    a = NeighbourhoodFileAnalyzer(neighb_path, "CP033163.1_o__pfams_EA459_08735.txt")
    print(a.count_pf("PF03060"))
    print(a.report_origin_stats())
    print(a.df_neigbr.iloc[1]['accession1'])
    """

    def __init__(self, path, filename):

        def tokenize_line(line):

            tokens = line.split()
            tokenized_line = tokens[0:22]
            tokenized_line.append(' '.join(tokens[22:]))

            return tokenized_line

        # Parse Filename ------------------------------------------------------------------

        self.genome_accession, self.locus_tag = parse_neighbourhood_filename(filename)

        columns = ['Target Name',
                   'accession1',
                   'tlen',
                   'query name',
                   'accession2',
                   'qlen',
                   'seq E-Value',
                   'seq score',
                   'seq bias',
                   'dom #',
                   'dom of',
                   'dom c-Evalue',
                   'dom i-Evalue',
                   'dom score',
                   'dom bias',
                   'hmm from',
                   'hmm to',
                   'ali from',
                   'ali to',
                   'env from',
                   'env to',
                   'acc',
                   'description of target']

        # Parse File ------------------------------------------------------------------

        raw_text = read_file(path + filename)
        lines = raw_text.splitlines()
        lines = [line for line in lines if not line.startswith('#')]
        self.tokenized_neigbr_lines = list()
        self.tokenized_origin_lines = list()

        for line in lines:

            tokenized_line = tokenize_line(line)

            if self.locus_tag in line:
                self.tokenized_origin_lines.append(tokenized_line)
            else:
                self.tokenized_neigbr_lines.append(tokenized_line)

        self.df_origin = pd.DataFrame(self.tokenized_origin_lines, columns = columns)
        self.df_neigbr = pd.DataFrame(self.tokenized_neigbr_lines, columns = columns)

        self.df_origin['accession1'] = self.df_origin['accession1'].apply(lambda s: s.split('.')[0])
        self.df_neigbr['accession1'] = self.df_neigbr['accession1'].apply(lambda s: s.split('.')[0])

    def count_pf(self, pf_number):
        pf_list = self.df_neigbr['accession1'].to_list()
        m = pf_list.count(pf_number)
        return m

    def report_origin_stats(self):

        # orig_locus_tags = list(set(self.df_origin['query name'].to_list()))
        origin_pfam_list = self.df_origin['accession1'].to_list()

        return origin_pfam_list

def mapped_from_blh_list(blh_list, mode='all_features', debug_mode=False):

    def debug_report(hit, ll, rl, lc, rc):

        print('\n\n')
        print('------------------------------------------------------------------------')
        print('Target Genome', hit.accession)
        print('Target Range', hit.hit_range)
        print('------------------------------------------------------------------------')
        print()
        print('GENOME FEATUREs:')
        print()
        print('LEFT: {}'.format(lc))
        print(ll)
        print('RIGHT: {}'.format(rc))
        print(rl)

    N = len(blh_list)
    rows = list()

    hit: BlastHit
    for i, hit in enumerate(blh_list):

        genome_hit = gbkGenome(hit.accession)

        ll, rl, lc, rc = genome_hit.get_coincident_features_from_hit_range(hit.hit_range, mode=mode)

        print('{}/{} - Processing: {}'.format(i + 1, N, hit.accession))

        if ll == rl:
            feature = ll
            feature_ambiguity = 0
        else:

            print(hit.accession, lc, rc, hit.hit_range)
            if (lc is None) and (rc is not None):
                feature = rl
                feature_ambiguity = 100-rc
            elif (lc is not None) and (rc is None):
                feature = ll
                feature_ambiguity = 100-lc
            elif (lc is None) and (rc is None):
                feature = None
                feature_ambiguity = 'N/A'
            elif lc > rc:
                feature = ll
                feature_ambiguity = rc
            elif lc < rc:
                feature = rl
                feature_ambiguity = lc
            else:
                raise ValueError('Leg coverage is equal!')

        if debug_mode: debug_report(hit, ll, rl, lc, rc)

        # Prepare the row for dataframe -----------------------------

        genome_annotated = len(genome_hit.record.features)
        genome_accession = genome_hit.record.id
        genome_description = genome_hit.record.description
        genome_organism = genome_hit.record.annotations['organism']
        genome_taxonomy = '/'.join(genome_hit.record.annotations['taxonomy'])
        genome_context = 'plasmid' if 'plasmid' in genome_description else 'chromosome'

        gbk_protein = genome_hit.convert_feature_to_custom_protein_class(feature)  # feature = None is handled in the method.

        feature_location = gbk_protein.location
        feature_locus_tag = gbk_protein.locus_tag
        feature_protein_id = gbk_protein.protein_id
        feature_product = gbk_protein.product
        feature_refseq = gbk_protein.refseq
        feature_translation = gbk_protein.translation
        feature_pseudo = gbk_protein.pseudo

        row = [hit.ridnum,
               genome_annotated,
               genome_accession.strip(),
               hit.total_score,
               hit.query_cover,
               hit.e_value,
               hit.per_ident,
               hit.hit_range,
               feature_location,
               feature_ambiguity,
               genome_description.strip(),
               genome_context,
               genome_organism.strip(),
               genome_taxonomy.strip(),
               feature_locus_tag.strip(),
               feature_protein_id.strip(),
               feature_product.strip(),
               feature_refseq.strip(),
               feature_pseudo,
               feature_translation]

        rows.append(row)

    columns = ['RID Job #',
               '# of Annotations',
               'Genome Accession',
               'Total Score',
               'Query Cover',
               'E Value',
               'Per. Identity.',
               'Hit Range',
               'Feature Location',
               'Loc. Ambiguity',
               'Description',
               'Context',
               'Organism',
               'Taxonomy',
               'Locus Tag',
               'Protein ID',
               'Product',
               'Ref Seq',
               'Pseudo',
               'Translation']

    df = pd.DataFrame(rows, columns = columns)
    df.set_index('Genome Accession', inplace=True)
    df.to_csv('mapped_blast_results.csv', sep=',')
    print(df)

def stitched_from_mapped(mapped_csv_path):

    # seq_record = SeqIO.read('proteins/{}.txt'.format(pid), "fasta")
    # seq_record = mba.download_protein_from_entrez('proteins/', pid)

    """

    This creates a stitched file from a mapped_csv (output from mapped_from_blh_list)
    It collapses all the entries with the same protein id into a single one and describe it with the shortest "organism" name (and count)
    Since the stitched file was not reliable, I went ahead and redownloaded them from entrez for each protein id.

    # stitched_from_mapped("2021-05-02 - manually_curated/8KMXX35H01R_and_8SUNF67Z013_mapped_blast_results_AH_Final.csv")

    """

    df = pd.read_csv(mapped_csv_path)

    cols = ['RID Job #',
            'Genome Accession',
            'Description',
            'Organism',
            'Taxonomy',
            'Locus Tag',
            'Protein ID',
            'Product',
            'Translation']

    df = df[cols]
    protein_ids = set(df['Protein ID'].to_list())

    rows = list()
    for pid in protein_ids:

        df_set = df[df['Protein ID'] == pid]

        plasmid = ''
        descriptions = df_set['Description'].tolist()
        n_descriptions = len(descriptions)
        for description in descriptions:
            if 'plasmid' in description: plasmid = 'plasmid'

        ridnum_list = list(set(df_set['RID Job #'].tolist()))
        ridnums = ', '.join(ridnum_list)

        organisms = list(set(df_set['Organism'].tolist()))
        organisms = sorted(organisms, key=len, reverse=False)
        organism = '; '.join(organisms)  # organisms[0] -> you can use this if you want only one organism name in there with the shortest string.

        tax = set(df_set['Taxonomy'].tolist())
        tax = list(set([t.split('/')[1] for t in tax]))[0]

        product = list(set(df_set['Product'].tolist()))[0]

        seq = list(set(df_set['Translation'].tolist()))

        if len(seq) > 1: print('Protein ID = {} is associated with more than one translation'.format(pid))

        row = '>{} {} [{}:{}, {}][{}] {{}}'.format(pid, product, tax, organism, n_descriptions, plasmid, ridnums)

        rows.append(row)
        rows.append(seq[0])

        print('--------------------------------------------------------------')
        print(row)
        print()
        print(seq[0])
        if len(seq) > 1: print('MULTIPLE SEQUENCES!!!!!')
        # time.sleep(0.5)

    text = '\n'.join(rows)
    write_file('stitched_fasta.txt', text)

def neighbourhood_files_from_mapped(mapped_csv_path, window_size=20):

    global neighb_path
    processed_accessions_list = list()

    for file in listdir(neighb_path):
        if '_pfams_' in file:
            genome_accession, locus_tag = parse_neighbourhood_filename(file)
            processed_accessions_list.append(genome_accession + '@' + locus_tag)

    df = pd.read_csv(mapped_csv_path)
    cols = ['RID Job #',
            'Genome Accession',
            'Description',
            'Organism',
            'Taxonomy',
            'Locus Tag',
            'Protein ID',
            'Product',
            'Translation']
    df = df[cols]
    index_list = df.index.to_list()
    L = len(index_list)

    ga_list = df['Genome Accession'].to_list()
    ga_set = set(ga_list)

    print('Redundant Genomes:', len(ga_list)-len(ga_set))
    print('Processed genome/locus pairs in the target folder: ', len(processed_accessions_list))

    pairs_in_csv = list()

    for n, index in enumerate(index_list):

        acc = df.loc[index]['Genome Accession']
        ltg = df.loc[index]['Locus Tag']

        pair = '{}@{}'.format(acc, ltg)
        pairs_in_csv.append(pair)

        if  pd.isna(ltg):

            print('{}/{}: {} has no locus tags'.format(n+1, L, acc))

        elif pair in processed_accessions_list:

            print('{}/{}: {}@{} already processed.'.format(n+1, L, acc, ltg))

        else:

            print(acc, ltg)
            t0 = time.perf_counter()

            genome = gbkGenome(acc)  # "NZ_CP059973.1", "GD631_RS01820"
            input_name, output_name = genome.get_neighbours_using_locus_tag(neighb_path, ltg, window_size=window_size)
            time.sleep(0.15)

            # Here is where the hmmer executable is ran to get pfams.
            hmmscan(neighb_path + input_name, neighb_path + output_name)

            t1 = time.perf_counter()
            print('{}/{}: {} processed in {} seconds'.format(n+1, L, input_name, round(t1 - t0, 2)))

    # The code below is just to inform if there are redundant genome accession + locus pairs in the mapped csv.
    # If you do, the number of csv entries won't match the number of neighbourhood analysis files generated.

    pair_dict = dict()
    for pair in pairs_in_csv:
        if pair not in pair_dict.keys():
            pair_dict[pair] = 1
        else:
            pair_dict[pair] = pair_dict[pair] + 1

    for pair in pair_dict.keys():
        if pair_dict[pair] > 1:
            print(pair, pair_dict[pair])

def neighbourhood_analysis(pfam_filename, neighb_path, mapped_csv_path):

    # Initialize

    df_pfams = pd.read_csv(pfam_filename)
    pfams = df_pfams['Pfam'].to_list()
    pfam_dict = dict()
    for pfam in pfams: pfam_dict[pfam] = [0, []]

    neighb_file_list = list()
    rows = list()

    # Code below just constructs a handy lookup table to map genome accessions to their descriptions and context.

    df_lookup = dict()
    df_mapped_csv = pd.read_csv(mapped_csv_path)
    for index in df_mapped_csv.index:
        accession   = df_mapped_csv.loc[index]['Genome Accession']
        description = df_mapped_csv.loc[index]['Description']
        context = df_mapped_csv.loc[index]['Context']
        df_lookup[accession] = [description, context]

    # Read into memory all the pfam files from neighbours.

    for file in listdir(neighb_path):
        if '_pfams_' in file: neighb_file_list.append(file)

    # CREATE results_main.csv ----------------------------------------------------------------------------------

    for file in neighb_file_list:

        detected_pfams = ''
        total_number_of_pfams = 0
        genome_accession, locus_tag = parse_neighbourhood_filename(file)
        obj = NeighbourhoodFileAnalyzer(neighb_path, file)
        origin_pfam_list = obj.report_origin_stats()

        origin_descriptor = '{}: [{}]'.format(locus_tag, ', '.join(origin_pfam_list))
        definition = df_lookup[genome_accession][0]
        location = df_lookup[genome_accession][1]

        for pfam in pfams:

            n = obj.count_pf(pfam)

            if n > 0:

                pfam_dict[pfam][1].append([genome_accession, n])
                detected_pfams += '{}:{}, '.format(pfam, n)
                total_number_of_pfams += n

        row = [genome_accession, definition, location, origin_descriptor, detected_pfams, total_number_of_pfams]
        rows.append(row)

    df = pd.DataFrame(rows, columns=['File Seq ID', 'Definition', 'Location', 'Origin Description', 'Detected pfams', 'Total'])
    df.sort_values(by=['Total'], inplace=True, ascending=False)
    df.set_index('File Seq ID', inplace=True)
    df.to_csv('results_main.csv', sep=',')
    print(df)

    # CREATE results_pfam.csv ----------------------------------------------------------------------------------

    rows = list()
    df_pfams.set_index('Pfam', inplace=True)
    for pfam in pfam_dict.keys():
        file_list = pfam_dict[pfam][1]
        N = len(file_list)
        name = str(df_pfams.loc[pfam]['Name'])
        rows.append([pfam, name, N, file_list])

    df2 = pd.DataFrame(rows, columns=['Pfam', 'Name', 'Number of Files', 'File List & Occurrences'])
    df2.sort_values(by=['Number of Files'], inplace=True, ascending=False)
    df2['Number of Files'] = df2['Number of Files'].map(str)
    df2.set_index('Pfam', inplace=True)
    df2.to_csv('results_pfam.csv', sep=',')
    print()
    print(df2)

def process_neighbourhood_analysis_files(file_paths):

    n_vec = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    rows = list()
    for key in file_paths.keys():

        results = pd.read_csv(file_paths[key])
        N = len(results)
        p_vec = list()

        for n in n_vec:
            df = results[results['Total'] >= n]
            ratio = round(100 * (len(df) / N), 2)
            p_vec.append(ratio)

        rows.append(p_vec)
        print(key)
        print(n_vec)
        print(p_vec)
        print()

    rows.append(n_vec)
    df_neighbourhood = pd.DataFrame(rows)
    df_t = df_neighbourhood.T

    df_t.columns = ['10AB', '20AB', '30AB', '40AB', 'x']
    df_t.set_index('x', inplace=True)
    df_t.to_csv('compiled_neighbourhood.csv', sep=',')
    print(df_t)

def plot_neighbourhood(path):

    df = pd.read_csv(path)

    x = df['x']
    x_major_ticks = np.arange(1, 11, 1)
    y_major_ticks = np.arange(0, 91, 10)

    plt.scatter(x=x, y=df['40AB'], c='blue', label='40AB', marker='o')
    plt.plot(x, df['40AB'], c='black')

    plt.scatter(x=x, y=df['30AB'], c='blue', label='30AB', marker='s')
    plt.plot(x, df['30AB'], c='black')

    plt.scatter(x=x, y=df['20AB'], c='blue', label='20AB', marker='*')
    plt.plot(x, df['20AB'], c='black')

    plt.scatter(x=x, y=df['10AB'], c='blue', label='10AB', marker='P')
    plt.plot(x, df['10AB'], c='black')
    plt.legend()
    plt.title('Neighbourhood Analysis')
    plt.xlabel('defense related pfams >')
    plt.ylabel('Fraction of Genomes (%)')
    plt.xticks(x_major_ticks)
    plt.yticks(y_major_ticks)
    plt.grid(color='black', linestyle='--', linewidth=0.5, which='major')
    plt.grid(color='green', linestyle='--', linewidth=0.5, which='minor')

    plt.show()

# Assembly Downloader

def read_summary_file(file):
    text_file = open(file, "r")
    my_text = text_file.read()
    text_file.close()
    return my_text

def read_large_text_linebyline(file, term):

    with open(file) as f:

        p = 1
        for line in f:
            if term in line: print(line)
            # if p < 3: print(line)
            p += 1

        print(p)


# COMMANDS -----------------------------------------------------------------------------

class Commands:

    @staticmethod
    def cmnd_single_hit_extraction(acc_name, hit_range):

        # acc_name = 'NZ_LR962337.1'  # 'HE980450.1'
        # hit_range = '50241-51803' # '56785-56830'

        genome_hit = gbkGenome(acc_name)  #
        ll, rl, lc, rc = genome_hit.get_coincident_features_from_hit_range(hit_range, mode='CDS')

        print(acc_name)
        print('Hit Range: ' + hit_range)
        print()
        print('L Overlap: {}'.format(lc))
        print(ll)
        print()
        print('R Overlap: {}'.format(rc))
        print(rl)

    # onetime
    @staticmethod
    def cmnd_find_genomes_in_green_file_not_in_tblastn():

        # This is to compare tekirs results to what our tblastn results gave us. What's in there and what's not.

        def get_df_green(filename):

            df_green = pd.read_csv(filename)
            df_green['Genome Accession'] = ''
            df_green['ProteinID'] = ''

            for index in df_green.index:
                string = df_green.loc[index]['Protein_ID']
                m = re.search(r'lcl\|(.+)_cds_(.+\.[0-9])_[0-9]+', string)
                df_green.loc[index]['Genome Accession'] = m.group(1)
                df_green.loc[index]['ProteinID'] = m.group(2)

            df_green = df_green.drop(columns=['RefSeq_ID', 'Protein_ID', 'Locus_Tag', 'Assembly_ID', 'Status', 'Location'])
            df_green['Genome in Mapped'] = False
            df_green['Protein ID in Mapped'] = False
            df_green.set_index('assembly_accession', inplace=True)
            df_green.sort_values(by='ProteinID', inplace=True)

            return df_green

        pd.set_option("max_colwidth", 40)

        df_green = get_df_green('104_original_unique.csv')  # df_green = pd.read_csv('from_green_file.csv')
        df = pd.read_csv('mapped_blast_results.csv')

        df_genomes = df['Genome Accession'].to_list()
        df_protein_id = df['Protein ID'].to_list()

        wp_dic = dict()
        for index in df_green.index:

            genome_accession = df_green.loc[index]['Genome Accession']
            protein_id = df_green.loc[index]['ProteinID']

            if genome_accession in df_genomes:
                df_green.loc[index, 'Genome in Mapped'] = True

            if protein_id in df_protein_id:
                df_green.loc[index, 'Protein ID in Mapped'] = True
                wp_dic[protein_id] = True
            else:
                wp_dic[protein_id] = False


        df_result = df_green[df_green['Genome in Mapped'] == False]
        df_result.to_csv('df_result.csv', sep = ',')

        print('\n\n--------------------------------')
        N =  len(set(df_green['ProteinID'].to_list()))# len(wp_dic.keys())
        n = 0
        for key in wp_dic.keys():
            if wp_dic[key]: n +=1
        print('{}/{} is in there'.format(n, N))

    @staticmethod
    def cmnd_parse_back_mapped():

        """
        This was intented to be used one time to read back a 'mapped_blast_results.csv' after Asma's manual curation
        and collapse all the entries with the same protein id into a single one and describe it with the shortest "organism" name (and count)
        Since the stitched file was not reliable, I went ahead and redownloaded them from entrez for each protein id.
        :return:
        """
        df = pd.read_csv('manually_curated_04_25_2021/04_24_21_mapped_blast_results_AH2.csv')

        cols = ['Genome Accession',
                'Description',
                'Organism',
                'Taxonomy',
                'Locus Tag',
                'Protein ID',
                'Product']

        df = df[cols]
        protein_ids = set(df['Protein ID'].to_list())

        rows = list()
        for pid in protein_ids:

            df_set = df[df['Protein ID'] == pid]
            n_descriptions = len(df_set['Description'].tolist())
            organisms = list(set(df_set['Organism'].tolist()))
            organisms = sorted(organisms, key=len, reverse=False)
            organism = organisms[0]
            tax = set(df_set['Taxonomy'].tolist())
            tax = list(set([t.split('/')[1] for t in tax]))[0]
            product = list(set(df_set['Product'].tolist()))[0]

            row = '>{} {} [{}:{}, {}]'.format(pid, product, tax, organism, n_descriptions)
            seq_record = SeqIO.read('proteins/{}.txt'.format(pid), "fasta")  # seq_record = mba.download_protein_from_entrez('proteins/', pid)
            rows.append(row)
            rows.append(str(seq_record.seq))

            print('--------------------------------------------------------------')
            print(row)
            print()
            print(seq_record)
            # time.sleep(0.5)

        text = '\n'.join(rows)
        write_file('new_stitched.txt', text)

    @staticmethod
    def cmnd_verify_protein_in_mapped():

        df = pd.read_csv('manually_curated_04_25_2021/04_24_21_mapped_blast_results_AH2.csv')
        for index in df.index:
            genome_accession = df.loc[index]['Genome Accession']
            protein_id = df.loc[index]['Protein ID']
            print()
            print()
            print('>CASE: ', genome_accession, protein_id)
            try:
                g0 = gbkGenome(genome_accession)
                fs = g0.find_protein_by_id(protein_id)
                if len(fs) == 0: print('The protein is not in here!')
                for f in fs: print(f)
            except:
                print('-----------------------------------------------')
                print('An Error Occurred!')
                print('-----------------------------------------------')

    @staticmethod
    def cmnd_parse_back_stiched():
        text = read_file('manually_curated_04_25_2021/04_24_21_stitched_fasta copy_AH2_ONELINE.txt')
        lines = text.splitlines()

        for line in lines:
            line_sections = line.split('@')
            string = line_sections[0]
            m = re.search(r'^>.+\.[0-9]_(.+) Ref', string)
            if m is None: print(string)

    @staticmethod
    def cmnd_temp():

        pairs = [('NZ_CP014791.1', 'WP_000632676.1'), ('NZ_AP022290.1', 'WP_000632676.1')]
        for pair in pairs:
            print('PAIR: ', pair)
            g0 = gbkGenome(pair[0])
            fs = g0.find_protein_by_id(pair[1])
            for f in fs: print(f)

    @staticmethod
    def main_RID_8SUNF67Z013():

        """This filters contigs and WGSs from 8SUNF67Z013, as well any any hits after that whose acession is in 8KMXX35H01R"""

        # xml = "2021-04-22-RID-83NB9BH901R/04_22_21_83NB9BH901R-Alignment_53567_all_proks_250_hits.xml"
        # a = ReadBlastOutputXML(xml)
        # a.process('CDS', False)
        # a.detailed_report()
        # for blh in a.get_blast_hit_object_list(): print(blh)
        # ReadBlastOutputXML.stitched_from_mapped('mapped_blast_results.csv')
        # a = ReadBlastOutputXML('BLAST_OUTPUTS/2021-04-29 - 8NVFJVEU013 - tblastn/8NVFJVEU013-Alignment.xml')
        # a = ReadBlastOutputXML("BLAST_OUTPUTS/2021-04-28 - 8KMXX35H01R - tblastn/8KMXX35H01R-Alignment.xml")
        # blh_list = a.get_blast_hit_object_list()
        # for blh in blh_list: print(blh)
        # mapped_from_blh_list(blh_list[0:5], mode ='CDS', debug_mode=False)
        # stitched_from_mapped('mapped_blast_results.csv')
        # accession_numbers = [hit.accession for hit in blist_filtered]
        # download_genome(genome_path, accession_numbers, rettype='gbwithparts', retmode='text')

        a = ReadBlastOutputXML("BLAST_OUTPUTS/2021-04-28 - 8KMXX35H01R - tblastn/8KMXX35H01R-Alignment.xml")
        alist = a.get_blast_hit_object_list()
        alist_acessions = [h.accession for h in alist]

        b = ReadBlastOutputXML("BLAST_OUTPUTS/2021-04-30 - 8SUNF67Z013 - tblastn/8SUNF67Z013-Alignment.xml")
        blist = b.get_blast_hit_object_list()

        blist_filtered = list()
        for hit in blist:
            if ('contig' not in hit.description) and ('shotgun' not in hit.description) and (
                    hit.accession not in alist_acessions):
                blist_filtered.append(hit)

        blist_filtered = sorted(blist_filtered, key=lambda x: x.description, reverse=False)

        for hit in blist_filtered: print(hit.description)

        mapped_from_blh_list(blist_filtered, mode='CDS', debug_mode=False)

    @staticmethod
    def main_RID_8KMXX35H01R():

        a = ReadBlastOutputXML("BLAST_OUTPUTS/2021-04-28 - 8KMXX35H01R - tblastn/8KMXX35H01R-Alignment.xml")
        alist = a.get_blast_hit_object_list()
        mapped_from_blh_list(alist, mode='CDS', debug_mode=False)

    @staticmethod
    def main_RID_8TFSMP75013():

        final_set = ['LN864705.1',
                     'CP031055.1',
                     'CP013020.1',
                     'CP072355.1',
                     'AP024464.1',
                     'AE011190.1',
                     'CP033953.1',
                     'CP009999.1',
                     'CP019727.1',
                     'CP023305.1',
                     'FP929043.1',
                     'FP929042.1',
                     'CP032160.1',
                     'KT186230.1',
                     'CP033163.1',
                     'FP929053.1',
                     'JN172911.1',
                     'CP073280.1',
                     'FP929062.1',
                     'AB450918.1',
                     'KT716758.1',
                     'AP017936.1',
                     'CP009909.1']

        b = ReadBlastOutputXML("BLAST_OUTPUTS/2021-04-30 - 8TFSMP75013 - tblastn/8TFSMP75013-Alignment.xml")
        blist = b.get_blast_hit_object_list()
        blist_filtered = [blh for blh in blist if blh.accession in final_set]
        mapped_from_blh_list(blist_filtered, mode='CDS', debug_mode=False)

    @staticmethod
    def compare_and_filter_blast_results_using_description_word_sets():
        def get_word_set(string: str):

            """
            this is for comparing genome description lines. The goal is to get rid of all the fluff words and
            reduce the description to a set of crutial words for comparing to another description from another
            blast result file
            """
            string = re.sub(r'\s*,\s*', ' ', string)
            string = re.sub(r'\s*:\s*', ' ', string)
            string = re.sub(r'unnamed_?[0-9]+', '', string)
            string = string.replace('chromosome', '')
            string = string.replace('plasmid', '')
            string = string.replace('genome', '')
            string = string.replace('complete', '')
            string = string.replace('sequence', '')
            string = string.replace('assembly', '')
            string = string.replace('DNA', 'chromosome')
            string = string.replace('strain', '')
            string = string.replace('str.', '')
            string = string.replace('sp.', '')
            string = string.replace('subsp.', '')
            string = string.replace('isolate', '')

            s = string.split()
            s = [word for word in s if len(word) > 1]
            s = set(s)
            return s

        """
        This reads two blast jobs and filters one of them by comparing to other based on chopping descriptions to unique identifying words.
        so that only the new ones that are not in the other RID job is processed.
        """

        a = ReadBlastOutputXML("BLAST_OUTPUTS/2021-04-28 - 8KMXX35H01R - tblastn/8KMXX35H01R-Alignment.xml")
        alist = a.get_blast_hit_object_list()
        a_descriptions = [get_word_set(hit.description) for hit in alist]

        b = ReadBlastOutputXML("BLAST_OUTPUTS/2021-04-30 - 8TFSMP75013 - tblastn/8TFSMP75013-Alignment.xml")
        blist = b.get_blast_hit_object_list()
        b_descriptions = [(get_word_set(hit.description), '{} {}'.format(hit.accession, hit.description)) for hit in
                          blist]

        for b_desc in b_descriptions:

            # print(b_desc)

            if b_desc[0] not in a_descriptions: print(b_desc[1])

    @staticmethod
    def check_mapped_csv_for_unique_translations():

        # This checks a mapped_csv (which is already filtered for unique protein ids)
        # if there are different protein ids with identical sequences.

        df = pd.read_csv("2021-05-02 - manually_curated/final mapped csv/8KMXX35H01R_and_8SUNF67Z013_mapped_blast_results_AH_Final.csv")
        translations = df['Translation'].to_list()
        translations_set = set(translations)
        print(len(translations))
        print(len(translations_set))

    @staticmethod
    def cmnd_check_if_proteins_are_all_unique_sequences():
        path = "manually_curated_04_25_2021/unique proteins/"
        files = [filename for filename in listdir(path) if '.txt' in filename]

        seqs = list()
        for file in files:
            seq_record = SeqIO.read(path + file, "fasta")
            seqs.append(seq_record.seq)

        seqs_set = set(seqs)

        print(len(seqs))
        print(len(seqs_set))

# SCRIPT -------------------------------------------------------------------------------

"""
Proteins we missed from the fuckers list.
{'WP_186731239.1', 'WP_147000324.1', 'WP_147001163.1', 'WP_194957901.1', 'WP_194853477.1', 'WP_156765833.1', 'WP_014566650.1', 'WP_193155371.1', 'WP_080276158.1', 'WP_000035333.1', 'WP_189695520.1', 'WP_165964726.1', 'WP_011679784.1', 'WP_033715838.1'}

df0 = pd.read_csv("code_development_files/FUCKER - Previous Analysis/ANALYSIS - NHI Unique Homolog /completed analysis/All Bacteria/104_unique_Nhi_homologs.csv")
fucker_ids = set(df0['WP_Number'].to_list())
df = pd.read_csv("2021-05-02 - manually_curated/8KMXX35H01R_and_8SUNF67Z013_mapped_blast_results_AH_Final.csv")
protein_ids = set(df['Protein ID'].to_list())
print(fucker_ids-protein_ids)

"BLAST_OUTPUTS/2021-04-29 - 8NVFJVEU013 - tblastn/8NVFJVEU013-Alignment.xml"
"BLAST_OUTPUTS/2021-04-28 - 8KMXX35H01R - tblastn/8KMXX35H01R-Alignment.xml"
"2021-05-02 - manually_curated/muscle-I20210503-054218-0867-84510672-p1m.clw"


"""

u = 5

if u == 1:
    ws = 40
    neighbourhood_files_from_mapped('2021-05-02 - manually_curated/step 1 - final mapped csv/8KMXX35H01R_and_8SUNF67Z013_and_8TFSMP75013_mapped_blast_results_AH_Final.csv', window_size=ws)

if u == 2:
    mapped_csv_path = '2021-05-02 - manually_curated/step 1 - final mapped csv/8KMXX35H01R_and_8SUNF67Z013_and_8TFSMP75013_mapped_blast_results_AH_Final.csv'
    neighb_path = "/Users/baslan/PycharmProjects/NHI_proximity_To_Defences/neighbours_window_size 40AB/"
    neighbourhood_analysis('PFAMS_defense_related.csv', neighb_path, mapped_csv_path)

if u == 3:

    file_paths = dict()
    file_paths['10AB'] = "2021-05-02 - manually_curated/step 4 - neighbourhood/Analysis Results/defense related_10AB/results_main_10AB.csv"
    file_paths['20AB'] = "2021-05-02 - manually_curated/step 4 - neighbourhood/Analysis Results/defense related_20AB/results_main_20AB.csv"
    file_paths['30AB'] = "2021-05-02 - manually_curated/step 4 - neighbourhood/Analysis Results/defense related_30AB/results_main_30AB.csv"
    file_paths['40AB'] = "2021-05-02 - manually_curated/step 4 - neighbourhood/Analysis Results/defense related_40AB/results_main_40AB.csv"
    process_neighbourhood_analysis_files(file_paths)

if u == 4:

    path = 'compiled_neighbourhood.csv'

    plot_neighbourhood(path)

if u == 5:

    fname = 'assemblies/assembly_summary/assembly_summary.txt'

    read_large_text_linebyline(fname, 'GCF_009662395')

