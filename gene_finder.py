# -*- coding: utf-8 -*-
"""
gene_finder testing 1

@author: Meaghen Sausville

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    if (nucleotide == 'A'):
        return 'T'
    elif (nucleotide == 'T'):
        return 'A'
    elif (nucleotide == 'G'):
        return 'C'
    elif (nucleotide == 'C'):
        return 'G'
    else:
        return "Not a valid nucleotide"


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    reversed1 = dna[::-1]
    output = ""
    for dnareversed in reversed1:
        output = output + get_complement(dnareversed)
    return output


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    i = 0  # index (where in the "line" the gene falls)
    output = ""  # gene from start codon without end codon

    while i < len(dna):
        codon = dna[i:i+3]
        if codon == 'TAG' or codon == 'TAA' or codon == 'TGA':
            return output
        output = output + codon
        i = i+3

    return output


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    i = 0
    ALLORF = []

    while i < len(dna):
        codon = dna[i:i+3]
        if codon == 'ATG':
            ORF = rest_of_ORF(dna[i:])
            ALLORF.append(ORF) #.append adds a string (ORF) to a list of strings (ALLORF)
            i += len(ORF)
        else:
            i = i+3
    return ALLORF


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    i = 0
    output = []
    while i < 3:
        output += find_all_ORFs_oneframe(dna[i:])  # += means take the right side of equation and add that to output, basically "output = output plus this."
        i = i+1
    return output


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    TotalORFs = []

    TotalORFs += find_all_ORFs(dna)
    TotalORFs += find_all_ORFs(get_reverse_complement(dna))

    return TotalORFs


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    ORFs = find_all_ORFs_both_strands(dna)
    longestORF = max(ORFs, key = len)

    return longestORF


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """

        i = 0
        For i < num_trials:
            dna = shufflestring(dna)[1]
            i += 1
        result = longest_ORF(dna)
        return result

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    i = 0
    aminostring = ''
    while i < len(dna):
        codon = dna[i:i+3]
        i += 3
        if len(codon) == 3:
            amino_acid = aa_table[codon]
            aminostring = aminostring + amino_acid
    return aminostring


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    output = coding_strand_to_AA(longest_ORF(dna))
    return output


if __name__ == "__main__":
    import doctest
    doctest.testmod()
