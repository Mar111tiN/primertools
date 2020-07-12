import re

def file2str(file):
    '''
    returns a string from a text file
    '''

    with open(file, 'r') as file:
        return file.read().upper().replace('\n', '')


def get_chrom(chrom, chroms_folder='chroms'):
    '''
    convenience function returning the chromosome sequence without
    the header sequences
    a folder containing the split chromosomes has to be provided
    '''

    if 'CHR' in str(chrom).upper():
        chrom = chrom.upper().replace('CHR', '')
    chrom_file = f"{chroms_folder}/chr{chrom}.fa"
    chrom_seq = file2str(chrom_file)
    chrom_seq = re.sub(r'^>CHR[0-9XY]+', '', chrom_seq)
    chrom = {
        "name": f"chr{chrom}",
        "sequence": chrom_seq,
        "length": len(chrom_seq)
    }
    return chrom 