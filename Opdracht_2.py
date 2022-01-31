#!/usr/bin/python3

"""
docstring
"""

__author__ = "Joshua Tolhuis"
__version__ = 0.1

import sys
import argparse
from Bio import AlignIO
import math


CODON_DICT = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'}


def argument_parser():
    """
    creates parser
    :return parser:
    """
    # create parser
    parser = argparse.ArgumentParser(
        description="This script determines whether entering a specified snp on a given postion\n"
                    "produces a severe mutation in amino acid sequence.",
        usage="python3 snp_annotation.py -n <SNP nucleotide>  -p <positie in het gen> -m <MSA>"
    )

    parser.add_argument('-n', '--SNP', type=str, metavar='',  required=True,
                        help='A single character representing the SNP nucleotide')
    parser.add_argument('-p', '--position', type=int, metavar='', required=True,
                        help='A number indicating the position of the SNP in the nucleotide sequence')
    parser.add_argument('-m', '--MSA', type=str, metavar='', required=True,
                        help='A clustal file of the MSA')
    return parser.parse_args()


def read_nucleotide_sequence(filename):
    """
    This function reads in a nucleotide sequence from a given file and returns
    the sequence as string.
        :return: The nucleotide sequence
    """
    with open(filename, 'r') as file:
        return ''.join(line.strip() for line in file if not line.startswith('>'))


def get_new_triplets(snp_pos, snp, sequence):
    """
    This function finds the old and new triplets of the snp, it assumes that the sequence provided
    is already in the correct reading frame, i.e the first triplet is the start codon at indices
    0, 1 and 2
    :param snp_pos: the position of the snp in the sequence (starting from 0) as int
    :param snp: the new nucleotide as string of length 1
    :param sequence: A nucleotide sequence as string
    :return triplets: The old and new (after snp) triplets at the snp position
    """
    # get the position of the snp within its triplet
    pos_in_triplet = snp_pos % 3
    # use this to select the triplet from the sequence as list (list aren't immutable)
    triplet = list(sequence[snp_pos - pos_in_triplet:snp_pos + (3 - pos_in_triplet)])
    # save original triplet
    old_triplet = ''.join(triplet)
    # insert snp and rejoin the new triplet
    triplet[pos_in_triplet] = snp
    triplet = ''.join(triplet)

    return old_triplet, triplet


def detect_change(triplets):
    """
    detects whether or not the snp creates any change in amino acid sequence
    :param triplets: tuple of two triplet strings of length 3 each
    :return boolean: boolean indicating if there's a change or not
    """
    return CODON_DICT[triplets[0]] == CODON_DICT[triplets[1]]


def read_msa(filename):
    """
    This functions parses the score lines from the msa file to a single string that will be
    the length of the longest sequence in the alignment.
    :param filename: name of msa file
    :return msa_scores: found scores
    """
    # start with empty string
    msa_scores = ''
    # open file
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith(' '):
                # each line that starts with a space contains scores after charcter 20
                msa_scores += line[20:len(line)].replace('\n', '')
    return msa_scores


def get_amino_acid_position(position, sequence, msa_scores):
    """
    changes given nucleotide position to position of mutation in msa
    :param position: snp nucleotide position
    :param sequence: nucleotide sequence (needed for length)
    :param msa_scores: msa scores found in the msa
    :return: position of mutation which can be used to index in msa_scores
    """
    position = math.ceil(position/3) - 1
    # if the msa is longer than the amount of amino acids in the sequence:
    if len(msa_scores) > len(sequence)/3:
        # let position start from difference
        position += len(msa_scores) - len(sequence)/3
    return int(position)


def main(args):
    """
    Main function
    :param args: arguments passed on from the commandline
    """
    # get values from arg parser
    snp = args.SNP.upper()
    position = args.position
    # check postion
    if not len(snp) == 1 or snp not in 'ATCG':
        print("error: SNP must be a single letter and may only be A, T, C or G")
        sys.exit()

    msa = args.MSA

    # get sequence from file
    sequence = read_nucleotide_sequence('./data/nucleotide/creatinekinase_nuc_hm.fasta')
    if position > len(sequence):
        print("error: position can't be higher than length of sequence")
        sys.exit()
    # get old and new triplets, inserting snp
    triplets = get_new_triplets(position, snp, sequence)
    # check for change in amino acid and get new amino acid
    change = detect_change(triplets)

    # if there is no change, there is no problem, exit program
    if not change:
        print('SNP did not create any changes in amino acid sequence: 4/4')
        sys.exit()
    # else, check severity of change
    # get msa scores from the msa files
    msa_scores = read_msa(msa)
    # get position of amino acid mutation in msa
    msa_pos = get_amino_acid_position(position, sequence, msa_scores)
    if msa_scores[msa_pos] == ' ':
        print("insignificant, mutation in non conserved domain: 4/4")
    if msa_scores[msa_pos] == '.':
        print("barely severe, mutation in mildly conserved domain: 3/4")
    if msa_scores[msa_pos] == ':':
        print("mildly severe, mutation in moderately conserved domain: 2/4")
    if msa_scores[msa_pos] == '*':
        print("severe, mutation in highly conserved domain: 1/4")

    return 0


if __name__ == '__main__':
    exit_code = main(argument_parser())
    sys.exit(exit_code)
