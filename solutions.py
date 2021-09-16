# Author: Jack Scacco
# Date: 9/15/21
# File: solutions.py

import sqlite3
from sqlite3 import Error

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
            for i in range(len(s) - k):
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
    return cur.lastrowid
    
    
def fill_kmer_table(kmer_freqs, database):
    """Given a dict of kmer freqs and a database file, fill the kmer table
       in that file with the info stored in the dict."""

    conn = create_connection(database)
    
    with conn:
        for kmer in kmer_freqs:
            kmer_info = (kmer, kmer_freqs[kmer])
            # (Not sure if we need this line, but I'm gonna stick with the tutorial
            kmer_id = enter_kmer(conn, kmer_info)

    conn.close()
    
def problem1b():
    """Synthesizes problem 1b."""

    # create_kmer_table()
    
    k = 21
    seqs = extract_seqs("SP1.fastq")
    k_mer_freqs = count_k_mers(seqs, k)
    database = "data.db"

    print("Filling table (this may take a while)...")
    fill_kmer_table(k_mer_freqs, database)
    
    
def main():
        
    problem1b()

if __name__ == "__main__":
    main()
