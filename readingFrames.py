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
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
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
        stringLength = len(s)
        # parse the data. use the regex
        output = "5'3' Frame 1 \r\n"
        output += process_frames(0, s, stringLength)
        output += "\r\n\r\n5'3' Frame 2 \r\n"
        output += process_frames(1, s, stringLength)
        output += "\r\n\r\n5'3' Frame 3 \r\n"
        output += process_frames(2, s, stringLength)

        reversed = s[::-1]
        print(reversed)
        output += "\r\n\r\n3'5' Frame 1 \r\n"
        output += process_frames(0, reversed, stringLength)
        output += "\r\n\r\n3'5' Frame 2 \r\n"
        output += process_frames(1, reversed, stringLength)
        output += "\r\n\r\n3'5' Frame 3 \r\n"
        output += process_frames(2, reversed, stringLength)



        # create the output file
        output_file = open("output.txt", "w")
        output_file.write(output)
        pass



def process_frames(startingIndex, input, stringLength):
    i = startingIndex
    returnString = ""
    while i < stringLength:
        check = input[i]
        for j in range(i + 1, i + 3):
            if j >= stringLength:
                continue
            check += input[j]
            # get the next 2
        if check not in gencode:
            break
        encoded = gencode[check]
        returnString += encoded
        i = i + 3
    return returnString


def main():
    if len(sys.argv) == 1:
        print("Please specify an input file")
        return
    read_dna(sys.argv[1])
    pass

main()
