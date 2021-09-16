# Author: Jack Scacco
# Date: 9/15/21
# File: solutions.py

import sys

def extract_seqs(filename):
    """Given a fastq file, return a list of the DNA sequences in that file."""

    all_seqs = []
    
    with open(filename, 'r') as f:
        
        this_line = f.readline()
        prev_line = ""

        # this_line will be empty (i.e. False) once we reach EOF
        while this_line:

            # If the line we read in before this one starts with "@cluster",
            # then we know this line contains a sequence.
            # (we chop off the last char bc it is \n)
            if prev_line[:8] == "@cluster":
                all_seqs.append(this_line[:-1])

            prev_line = this_line
            this_line = f.readline()

        return all_seqs
        

def count_k_mers(seqs, k):
    """Given a list of sequences, return a frequency table of each k-mer."""

    k_mer_freqs = {}
    
    for s in seqs:
        # Prevent indexing errors by making sure k is not larger than
        # the length of the seq
        if k <= len(s):
            # i represents the starting index of each k-mer in the seq
            for i in range(len(s) - k):
                k_mer = s[i:i+k]
                if k_mer in k_mer_freqs:
                    k_mer_freqs[k_mer] += 1
                else:
                    k_mer_freqs[k_mer] = 1

    return k_mer_freqs
                

def problem1a():
    k = 21
    seqs = extract_seqs("SP1.fastq")
    k_mer_freqs = count_k_mers(seqs, k)

    # View k_mer_freqs
    # for k in k_mer_freqs:
    #    print("\nk-mer:", k, "\nfrequency:", k_mer_freqs[k])

    # Make sure we aren't missing any. (this works bc all seqs are same length)
    # Compute expected number of k_mers
    num_seqs = len(seqs)
    seq_len = len(seqs[0])
    k_mers_per_seq = seq_len - k
    total_k_mers_expected = num_seqs * k_mers_per_seq

    # Compute acutal number of k_mers
    total_k_mers_actual = 0
    for k in k_mer_freqs:
        total_k_mers_actual += k_mer_freqs[k]

    assert(total_k_mers_expected == total_k_mers_actual)

    
def main():
    problem1a()

if __name__ == "__main__":
    main()
