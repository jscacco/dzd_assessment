# Author: Jack Scacco
# Date: 9/15/21
# File: solutions.py

import sqlite3
from sqlite3 import Error
import random

BASES = ["A", "C", "G", "T"]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                        Problem 1a                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

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
            for i in range(len(s) - k + 1):
                k_mer = s[i:i+k]
                if k_mer in k_mer_freqs:
                    k_mer_freqs[k_mer] += 1
                else:
                    k_mer_freqs[k_mer] = 1

    return k_mer_freqs
                

def problem1a():
    """Synthesizes problem 1a"""
    
    k = 21
    seqs = extract_seqs("SP1.fastq")
    k_mer_freqs = count_k_mers(seqs, k)

    # View k_mer_freqs
    for k in k_mer_freqs:
       print("\nk-mer:", k, "\nfrequency:", k_mer_freqs[k])

    # Make sure we aren't missing any. (this works bc all seqs are same length)
    # Compute expected number of k_mers
    num_seqs = len(seqs)
    seq_len = len(seqs[0])
    k_mers_per_seq = seq_len - k + 1
    total_k_mers_expected = num_seqs * k_mers_per_seq

    # Compute acutal number of k_mers
    total_k_mers_actual = 0
    for k in k_mer_freqs:
        total_k_mers_actual += k_mer_freqs[k]

    assert(total_k_mers_expected == total_k_mers_actual)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                        Problem 1b                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def create_connection(db_file):
    """Creates a new database file. Copied from 
       https://www.sqlitetutorial.net/sqlite-python/creating-database/"""

    conn = None
    try:
        # remember to close this connection!
        conn = sqlite3.connect(db_file)
        print(sqlite3.version)
        
    except Error as e:
        print(e)

    return conn


def create_table(conn, create_table_sql):
    """Creates a table from the create_table_sql statement. Copied from tutorial."""
    try:
        c = conn.cursor()
        c.execute(create_table_sql)
    except Error as e:
        print(e)


def create_kmer_table():
    """Creates the kmer table (doesn't need to be called after the first time)."""
    
    sql_create_kmer_table = \
        """ CREATE TABLE IF NOT EXISTS kmer (
        k_mer text PRIMARY KEY,
        count integer NOT NULL
        ); """

    # Create connection
    database = "data.db"
    conn = create_connection(database)

    # Create table
    if conn is not None:
        create_table(conn, sql_create_kmer_table)
    else:
        print("Error! Cannot create the database connection.")

    conn.close()


def enter_kmer(conn, kmer_info):
    """Given a Connection object and a tuple representing a kmer and its frequency,
       add that info into the kmer table. Returns the rowid of that info."""

    sql = \
        """ INSERT INTO kmer(k_mer,count)
        VALUES(?,?) """

    cur = conn.cursor()
    cur.execute(sql, kmer_info)
    conn.commit()

    
def fill_kmer_table(kmer_freqs, database):
    """Given a dict of kmer freqs and a database file, fill the kmer table
       in that file with the info stored in the dict."""

    conn = create_connection(database)
    
    with conn:
        for kmer in kmer_freqs:
            kmer_info = (kmer, kmer_freqs[kmer])
            enter_kmer(conn, kmer_info)

    conn.close()

    
def problem1b():
    """Synthesizes problem 1b."""

    # No need to call this after the first time
    # create_kmer_table()
    
    k = 21
    seqs = extract_seqs("SP1.fastq")
    k_mer_freqs = count_k_mers(seqs, k)
    database = "data.db"

    # This code was tested using the command-line sqlite3 commands described in the tutorial:
    # >sqlite3 data.sb
    # sqlite> .header on
    # sqlite> .mode column
    # sqlite> SELECT * from kmer;

    print("Filling table (this may take a while)...")
    fill_kmer_table(k_mer_freqs, database)
    print("Done filling table.")

    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                        Problem 2                                          #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def num_letters_diff(first, second):
    """Given two strings of the same length (first and second), return the number of 
       letters by which they differ."""

    assert(len(first) == len(second))

    num_diff = 0
    
    for i in range(len(first)):
        if first[i] != second[i]:
            num_diff += 1

    return num_diff


def match(kseq, seq):
    """Given a k-mer and a DNA sequency longer than k, return a set containing
       all k-mers in seq that match kseq with at most two letters different."""

    matches = []
    k = len(kseq)

    # Look at all possible k-mers in seq and see how many letters they differ from kseq
    for i in range(len(seq) - k + 1):
        k_mer = seq[i:i+k]
        if k_mer not in matches and num_letters_diff(kseq, k_mer) <= 2:
            matches.append(k_mer)

    # Convert the list to a set and return
    return set(matches)


def problem2():
    """Synthesizes problem 2."""

    # kseq = "ACGT"
    # seq = "ACACACGT"

    kseq = ""
    seq = ""

    # Modify these to your heart's content for random test cases
    seq_len = 2000
    kseq_len = 5
    
    # Test on arbitrarily large kseq and seq (less likely to have matches the longer
    # kseq is, though)
    for i in range(seq_len):
        if i < k_seq_len:
            kseq += random.choice(BASES)
        seq += random.choice(BASES)

    print("kseq:", kseq)
    print("seq:", seq)
    
    matches = match(kseq, seq)
    print("matches:")
    for m in matches:
        print(m)

        
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                          Main                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def main():
    # problem1a()
    # problem1b()
    problem2()
    # print(match("ACGT", "ACACACGT"))
    
    
if __name__ == "__main__":
    main()
