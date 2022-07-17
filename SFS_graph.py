#!/usr/bin/env python3

import sys, time, argparse, subprocess, os, signal, shutil, copy
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

Description = """
Rimerge
"""

#------------------------------------------------------------
def remove_file(file_path):
    os.remove(file_path)

def remove_dir(dir_path):
    shutil.rmtree(dir_path)

def move_dir_content(dir_source, dir_dest):
    files = os.listdir(dir_source)
    for f in files:
        shutil.move(dir_source + "/" + f, dir_dest)

#------------------------------------------------------------
# execute command: return True if everything OK, False otherwise
def execute_command(command, seconds=10000000):
    try:
        process = subprocess.Popen(command.split(), preexec_fn=os.setsid)
        process.wait(timeout=seconds)
    except subprocess.CalledProcessError:
        print("Error executing command line:")
        print("\t"+ command)
        return False
    except subprocess.TimeoutExpired:
        os.killpg(os.getpgid(process.pid), signal.SIGTERM)
        print("Command exceeded timeout:")
        print("\t"+ command)
        return False
    return True

#------------------------------------------------------------
# Return a list of super sample specific strings from a list of sample specific
def super_specific_strings(sample_specific_strings):
    out_super_sample_specifcs = list()
    super_sample_specific = copy.deepcopy(sample_specific_strings[0])
    i = 1
    while i < len(sample_specific_strings) - 1:
        if (sample_specific_strings[i][1] < super_sample_specific[1] + len(super_sample_specific[0])):
            super_sample_specific[0] += sample_specific_strings[i][0][int(sample_specific_strings[i][1]) - int(super_sample_specific[1]) : ]
        else:
            out_super_sample_specifcs.append(super_sample_specific)
            super_sample_specific = copy.deepcopy(sample_specific_strings[i])
        i += 1

    out_super_sample_specifcs.append(super_sample_specific)

    return out_super_sample_specifcs

def extract_mems(lengths_file_path, pointers_file_path):
    lenghts_lines = list()
    pointers_lines = list()

    with open(lengths_file_path) as lengths_file:
        lenght_lines = lengths_file.readlines()

    with open(pointers_file_path) as pointers_file:
        pointers_lines = pointers_file.readlines()

    mems = dict()
    for record_id in range(int(len(lenght_lines) / 2)):
        read_name = lenght_lines[record_id * 2].strip()
        lengths   = [int(l) for l in lenght_lines[record_id * 2 + 1].strip().split(" ")]
        pointers  = [int(p) for p in pointers_lines[record_id * 2 + 1].strip().split(" ")]

        mems[read_name] = list()
        mems[read_name].append([pointers[0], pointers[0] + lengths[0], lengths[0] , 0])

        i = 1
        while i < (len(lengths) - 2):
            if (lengths[i] < lengths[i + 1]):
                mems[read_name].append([pointers[i + 1], pointers[i + 1] + lengths[i + 1],  lengths[i + 1] , i + 1])
            i += 1
    return mems

def read_sss_file(sfs_file_path):
    sample_specific_strings = dict()
    max_mem_per_read = dict()
    with open(sfs_file_path, 'rb') as sss_file:
        sss_read_name_size = int.from_bytes(sss_file.read(8), 'little')
        while sss_read_name_size != 0:
            sss_list = list()

            sss_read_name = sss_file.read(sss_read_name_size).decode('ascii')
            sss_max_mem_pos = int.from_bytes(sss_file.read(8), 'little')
            sss_max_mem_idx = int.from_bytes(sss_file.read(8), 'little')
            sss_max_mem_len = int.from_bytes(sss_file.read(8), 'little')

            max_mem_per_read[sss_read_name] = [sss_max_mem_pos, sss_max_mem_idx, sss_max_mem_len]

            ss_list_size = int.from_bytes(sss_file.read(8), 'little')
            for ss_i in range(ss_list_size):
                sss_string_size = int.from_bytes(sss_file.read(8), 'little')
                sss_string      = sss_file.read(sss_string_size).decode('ascii')
                sss_read_pos    = int.from_bytes(sss_file.read(8), 'little')
                sss_ref_pos     = int.from_bytes(sss_file.read(8), 'little')

                if (sss_string_size >= 3):
                    sss_list.append([sss_string, sss_read_pos, sss_ref_pos])

            sample_specific_strings[sss_read_name] = sss_list

            # keep looping
            sss_read_name_size = int.from_bytes(sss_file.read(8), 'little')

    return sample_specific_strings, max_mem_per_read

# Very basic workflow
# - build moni index on reference, assuming fasta
# - align reads computing SFSs
# - build themisto on the SFSs
def main():
    parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-r', '--reference', help='reference file name', type=str, required=True, dest='reference')
    parser.add_argument('-p', '--pattern', help='the input query', type=str, required=True, dest='pattern')
    args = parser.parse_args()

    # Build index
    build_index = "moni build -r {} -f".format(args.reference)
    execute_command(build_index)

    # Compute SFSs
    compute_sfs = "moni sample-specific -i {} -p {}".format(args.reference, args.pattern)
    execute_command(compute_sfs)

    # Build Graph
    pointers_file_name = args.pattern + "_" + args.reference + ".pointers"
    lengths_file_name  = args.pattern + "_" + args.reference + ".lengths"
    sfs_file_name      = args.pattern + "_" + args.reference + "_0.ss.tmp.out"
    out_sequences      = args.pattern + ".sequences.fasta"
    out_colors         = args.pattern + ".colors.txt"
    out_graph_prefix   = args.pattern + "_" + args.reference + "_DBG"

    colors_list = []

    with open(out_sequences, 'w') as out_sequences_handle:
        string_reference = str(SeqIO.read(args.reference, "fasta").seq)
        color = 2 # 0,1 reserved for mems and SFSs
        with open(args.pattern) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                colors_list.append(color)
                color += 1
                SeqIO.write(record, out_sequences_handle, "fasta")

        # get mems
        mems_per_read = extract_mems(lengths_file_name, pointers_file_name)
        for read, mems_list in mems_per_read.items():
            for mem in mems_list:
                colors_list.append(0)
                seq = Seq(string_reference[mem[0]:mem[1]])
                record = SeqRecord(seq, "MEM From:\t{}".format(read), "", "")
                SeqIO.write(record, out_sequences_handle, "fasta")

        # get sss
        sss_per_read, _ = read_sss_file(sfs_file_name)
        for read, sss_list in sss_per_read.items():
            for super_sss in super_specific_strings(sss_list):
                colors_list.append(1)
                seq = Seq(super_sss[0])
                record = SeqRecord(seq, "SSS From:\t{}".format(read), "", "")
                SeqIO.write(record, out_sequences_handle, "fasta")

    with open(out_colors, 'wt') as out_colors_handle:
        for color in colors_list:
            out_colors_handle.write(str(color) + '\n')

    # Create graph
    create_graph = "themisto build --node-length 30 -i {} -c {} --index-prefix {} --temp-dir tmp".format(out_sequences, out_colors, out_graph_prefix)
    execute_command(create_graph)

if __name__ == '__main__':
    main()