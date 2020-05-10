import os
import argparse
import math
import time

def split_file(file_name, split):
    num_lines = sum(1 for line in open(file_name))
    lines_per_file = math.ceil(num_lines / split) #split into 10 files by defualt
    smallfile = None
    keep_going = False
    file_number = 0
    with open(file_name) as bigfile:
        for lineno, line in enumerate(bigfile):
            if lineno % lines_per_file == 0 or keep_going == True:
                if line.startswith(">"): #make sure to split on a new entry
                    if smallfile:
                        smallfile.close()
                    file_number += 1
                    small_filename = 'small_file_{}_{}'.format(file_name,file_number)
                    smallfile = open(small_filename, "w")
                    keep_going = False
                else:
                    keep_going = True
            smallfile.write(line)
        if smallfile:
            smallfile.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', help="Specify file for which to search. It should be in the container or be a downloadable link", type = str, default = "influenza_xsmall.faa", required = False)
    parser.add_argument('-s', '--split', help="Specify how many files we should split this into. Default is 10", type = int, default = 10, required = False)
    arguments = parser.parse_args()

    start_time = time.time() #Timing execution
    split_file(file_name = arguments.file, split = arguments.split)
    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()
