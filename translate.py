#! /usr/bin/env python3

import sys

def translate_sequence(rna_sequence, genetic_code):
    amino = ""
    seq = rna_sequence.upper()
    for i in range(0,len(seq)-len(seq)%3,3):
        if len(seq)%3 != 0:
            r = len(seq)%3
            seq = rna_sequence[:len(rna_sequence)-r]
            if genetic_code[seq[i:i+3]] == "*":
                break
            else:
                amino += genetic_code[seq[i:i+3]]
        elif len(seq)%3 == 0:
            if genetic_code[seq[i:i+3]] == "*":
                break
            else:
                amino += genetic_code[seq[i:i+3]]
    return amino
    """Translates a sequence of RNA into a sequence of amino acids.
    Translates `rna_sequence` into string of amino acids, according to the
    `genetic_code` given as a dict. Translation begins at the first position of
    the `rna_sequence` and continues until the first stop codon is encountered
    or the end of `rna_sequence` is reached.

    If `rna_sequence` is less than 3 bases long, or starts with a stop codon,
    an empty string is returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    str
        A string of the translated amino acids.
    """

def get_all_translations(rna_sequence, genetic_code):
    amino_acids = []
    pep = []
    seq = rna_sequence.upper()
    frames = [seq[i:i+3] for i in range(0, len(seq))]
    for i in range(0, len(seq)):
        frames = []
        frames.append([seq[i:i+3] for i in range(0,len(seq),3)])
        frames.append([seq[i:i+3] for i in range(1,len(seq),3)])
        frames.append([seq[i:i+3] for i in range(2,len(seq),3)])
    print(frames)
    for i in range(0,len(frames)):
        start = 0 #changing this number removes the number of peptides we can find
        while start < len(frames[i]):
            if frames[i][start] == "AUG":
                #amino_acids.append(frames[i][start:]) 
                for stop in range(start,len(frames[i])):
                    if frames[i][stop]=="UAA" or  frames[i][stop]=="UGA" or frames[i][stop]=="UAG":
                        amino_acids.append(frames[i][start:stop])
                        start=stop+1
                        break
                                 
            start+=1
    for i in range(0,len(amino_acids),1):
       final_seq = ''.join(amino_acids[i])
       pep.append(translate_sequence(final_seq,genetic_code))
    return pep
    
# You are getting close on this! If you look at the error messages, you will see that it is printing one protein sequences in some cases where it expects 2. 
# It looks like your code is only grabbing a start codon if the protein sequence has a stop codon with the given start codon.
# It should be set to start coding at a start codon, even if there is not a stop codon in that reading frame.
# See if you can change your for stop in range() for loop to use .append on your amino_acids even if you don't run into a stop codon.


    """Get a list of all amino acid sequences encoded by an RNA sequence.

    All three reading frames of `rna_sequence` are scanned from 'left' to
    'right', and the generation of a sequence of amino acids is started
    whenever the start codon 'AUG' is found. The `rna_sequence` is assumed to
    be in the correct orientation (i.e., no reverse and/or complement of the
    sequence is explored).

    The function returns a list of all possible amino acid sequences that
    are encoded by `rna_sequence`.

    If no amino acids can be translated from `rna_sequence`, an empty list is
    returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    list
        A list of strings; each string is an sequence of amino acids encoded by
        `rna_sequence`.
    """

def get_reverse(sequence):
    rev = sequence[::-1]
    rev = rev.upper()
    if sequence is None:
    	rev = ""
    return rev

    """Reverse orientation of `sequence`.

    Returns a string with `sequence` in the reverse order.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> get_reverse('AUGC')
    'CGUA'
    """

def get_complement(sequence):
    seq = sequence.upper()
    comp = seq.replace('A', 'u').replace('U', 'a').replace('C', 'g').replace('G', 'c')
    comp = comp.upper()
    return comp
    """Get the complement of a `sequence` of nucleotides.

    Returns a string with the complementary sequence of `sequence`.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> get_complement('AUGC')
    'UACG'
    """

def reverse_and_complement(sequence):
    rev = sequence[::-1]
    rev = rev.upper()
    comp = rev.replace('A', 'u').replace('U', 'a').replace('C', 'g').replace('G', 'c')
    comp_rev = comp.upper()
    return comp_rev
    """Get the reversed and complemented form of a `sequence` of nucleotides.

    Returns a string that is the reversed and complemented sequence
    of `sequence`.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> reverse_and_complement('AUGC')
    'GCAU'
    """
    

def get_longest_peptide(rna_sequence, genetic_code):
    """Get the longest peptide encoded by an RNA sequence.

    Explore six reading frames of `rna_sequence` (the three reading frames of
    `rna_sequence`, and the three reading frames of the reverse and complement
    of `rna_sequence`) and return (as a string) the longest sequence of amino
    acids that it encodes, according to the `genetic_code`.

    If no amino acids can be translated from `rna_sequence` nor its reverse and
    complement, an empty string is returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    str
        A string of the longest sequence of amino acids encoded by
        `rna_sequence`.
    """    


if __name__ == '__main__':
    genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'}
    rna_seq = ("AUG"
            "UAC"
            "UGG"
            "CAC"
            "GCU"
            "ACU"
            "GCU"
            "CCA"
            "UAU"
            "ACU"
            "CAC"
            "CAG"
            "AAU"
            "AUC"
            "AGU"
            "ACA"
            "GCG")
    longest_peptide = get_longest_peptide(rna_sequence = rna_seq,
            genetic_code = genetic_code)
    assert isinstance(longest_peptide, str), "Oops: the longest peptide is {0}, not a string".format(longest_peptide)
    message = "The longest peptide encoded by\n\t'{0}'\nis\n\t'{1}'\n".format(
            rna_seq,
            longest_peptide)
    sys.stdout.write(message)
    if longest_peptide == "MYWHATAPYTHQNISTA":
        sys.stdout.write("Indeed.\n")
