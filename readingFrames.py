import re
import sys

# use s = re.sub('[^acgtACGT]+', '', seq).upper() (or lower?) may or may not need the + sign
# for regex to remove all non alpha

# dictionary found here: https://pythonforbiologists.com/dictionaries/
gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'
}


def read_dna(filename):
    """
    :param filename:
    :return:
    """
    # input will  be seq.fasta.txt
    with open(filename, 'r') as infile:

        seq = infile.read()
        # loop through the file, check against the dictionary
        s = re.sub('[^acgtACGT]', '', seq).upper()
        # parse the data. use the regex
        for key in gencode.keys():
            s = s.replace(key, gencode[key])
        print(s)
        # read from filename and create dictionary

        # create the output file
        output = open("output.txt", "w")
        output.write(s)
        pass


def complement(seq, comp):
    """
    :param seq:
    :param comp:
    :return:
    /* input the exception here */
    :raises:
    """
    pass



def find_reading_frame(seq, pos, codon_table):
    pass


def write_frames(gene, rev_complement, filename):
    pass


def main():
    if len(sys.argv) == 1:
        print("Please specify an input file")
        return
    read_dna(sys.argv[1])
    pass


if __name__ == 'main':
    main()

main()
