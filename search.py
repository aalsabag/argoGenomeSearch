import os
import argparse
import math
import time
import glob
import shutil

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_protein

def regular_search(seq, file_name):
    i = 0
    total_count = 0
    result_file = os.path.basename(file_name)
    f= open(f"result_{result_file}","w+")
    for seq_record in SeqIO.parse(file_name, "fasta"):
        total_count += 1
        index = seq_record.seq.find(seq)
        if (index != -1):
            i += 1
            print(f"Found Something!")
            print(f"ID: {seq_record.id}")
            f.write(f"ID: {seq_record.id}\n")
    f.close()
    print(total_count)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', help="Specify file for which to search. It should be in the container or be a downloadable link", type = str, default = "influenza_xsmall.faa", required = False)
    parser.add_argument('-s', '--sequence', help="Specify sequence of amino acids. ie. NDVTSL", type = str)
    parser.add_argument('-d', '--directory', help="Optional. Copy results to certain directory", type = str)
    arguments = parser.parse_args()

    start_time = time.time() #Timing execution
    regular_search(file_name = arguments.file, seq = arguments.sequence)
    if arguments.directory is not None:
        for file in glob.glob("result*"):
            shutil.copy(file, arguments.directory)
    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()
