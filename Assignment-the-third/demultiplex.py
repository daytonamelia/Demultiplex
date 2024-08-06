#!/usr/bin/env python

import bioinfo, argparse, gzip, re
import matplotlib.pyplot as plt

def get_args():
    parser = argparse.ArgumentParser(description="Given 4 files from an illumina lane, checks for matching, hopped, unknown indexes.")
    parser.add_argument("-r1", "--read1", help="Path to R1 file", type=str, required=True)
    parser.add_argument("-r2", "--read2", help="Path to R2 file", type=str, required=True)
    parser.add_argument("-r3", "--read3", help="Path to R3 file", type=str, required=True)
    parser.add_argument("-r4", "--read4", help="Path to R4 file", type=str, required=True)
    parser.add_argument("-i", "--index", help="Path to index file for parsing.", type=str, required=True)
    parser.add_argument("-o", "--output", help="Output directory path.", type=str, default="demultiplexed/")
    parser.add_argument("-c", "--cutoff", help="Phred cutoff for reads, meaning anything BELOW will be cutoff.", default=30)
    parser.add_argument("-p", "--phredscore", help="Phred score scheme.", type=int, default=33)
    return parser.parse_args()

def valid_indexes(file: str) -> dict:
    '''Takes an input file and returns a dictionary of valid illumina barcodes.'''
    index_dict = {}
    with open(file, "r") as fh:
        for line in fh:
            line = line.strip("\n")
            line_list = re.split("\t", line)
            if line_list[0].startswith('s'):
                continue
            index_dict[line_list[3]] = line_list[4]
    return index_dict

def lowest_qscore(qscore: str, phred: int) -> str:
    '''Takes a string of qscore and returns the lowest qscore across the string.'''
    lowestq = "z"
    for score in qscore:
        score_ord = ord(score) - phred
        if score_ord < (ord(lowestq) - phred):
            lowestq = score
    return lowestq

print("Parsing index file and making/opening output files...", flush=True)

args = get_args()
index_list = valid_indexes(args.index)  # valid indexes from index file

# make matched index output files
matched_indexes = {}
output_files = []
for key, val in index_list.items():
    matched_indexes[val] = 0
    output_files.append(f"{args.output}demux_matched_{val}_R1.fq")
    output_files.append(f"{args.output}demux_matched_{val}_R2.fq")

# make hopped index output files
hopped_indexes = {}
output_files.append(f"{args.output}demux_hopped_R1.fq")
output_files.append(f"{args.output}demux_hopped_R2.fq")

# make unknown index output files
unknown_indexes = {}
output_files.append(f"{args.output}demux_unknown_R1.fq")
output_files.append(f"{args.output}demux_unknown_R2.fq")

# open all output files
# https://stackoverflow.com/questions/29550290/how-to-open-a-list-of-files-in-python
openoutput = {filename: open(filename, 'w') for filename in output_files}

print("Reading input files...", flush=True)

# open and analyze all input files
with gzip.open(args.read1, "rt") as r1, gzip.open(args.read2, "rt") as r2, gzip.open(args.read3, "rt") as r3, gzip.open(args.read4, "rt") as r4:
    line_counter = 0
    while True:
        # get each line
        line1 = r1.readline()
        line1 = line1.strip('\n')
        line2 = r2.readline()
        line2 = line2.strip('\n')
        line3 = r3.readline()
        line3 = line3.strip('\n')
        line4 = r4.readline()
        line4 = line4.strip('\n')
        if line1.startswith('@'):
            headerstate = True
        # meat of the code in here...
        if headerstate and line_counter != 0:
            # check if indexes are valid indexes
            index1 = record2[1]
            index1_qscore = record2[3]
            index2 = bioinfo.revcomp(record3[1])
            index2_qscore = record2[3]
            # index 1 is unknown
            if index1 not in index_list.values() or bioinfo.convert_phred(lowest_qscore(index1_qscore, args.phredscore)) < args.cutoff:
                # write to R1 file
                openoutput[f"{args.output}demux_unknown_R1.fq"].write(f"{record1[0]} {index1}-{index2}\n")
                openoutput[f"{args.output}demux_unknown_R1.fq"].write(f"{record1[1]}\n")
                openoutput[f"{args.output}demux_unknown_R1.fq"].write(f"{record1[2]}\n")
                openoutput[f"{args.output}demux_unknown_R1.fq"].write(f"{record1[3]}\n")
                # write to R2 file
                openoutput[f"{args.output}demux_unknown_R2.fq"].write(f"{record4[0]} {index1}-{index2}\n")
                openoutput[f"{args.output}demux_unknown_R2.fq"].write(f"{record4[1]}\n")
                openoutput[f"{args.output}demux_unknown_R2.fq"].write(f"{record4[2]}\n")
                openoutput[f"{args.output}demux_unknown_R2.fq"].write(f"{record4[3]}\n")
                # make stats for unknown_indexes
                if f"{index1}-{index2}" not in unknown_indexes:
                    unknown_indexes[f"{index1}-{index2}"] = 1
                else:
                    unknown_indexes[f"{index1}-{index2}"] +=1
            # index 1 is OK but index 2 is unknown
            elif index2 not in index_list.values() or bioinfo.convert_phred(lowest_qscore(index2_qscore, args.phredscore)) < args.cutoff:
                # write to R1 file
                openoutput[f"{args.output}demux_unknown_R1.fq"].write(f"{record1[0]} {index1}-{index2}\n")
                openoutput[f"{args.output}demux_unknown_R1.fq"].write(f"{record1[1]}\n")
                openoutput[f"{args.output}demux_unknown_R1.fq"].write(f"{record1[2]}\n")
                openoutput[f"{args.output}demux_unknown_R1.fq"].write(f"{record1[3]}\n")
                # write to R2 file
                openoutput[f"{args.output}demux_unknown_R2.fq"].write(f"{record4[0]} {index1}-{index2}\n")
                openoutput[f"{args.output}demux_unknown_R2.fq"].write(f"{record4[1]}\n")
                openoutput[f"{args.output}demux_unknown_R2.fq"].write(f"{record4[2]}\n")
                openoutput[f"{args.output}demux_unknown_R2.fq"].write(f"{record4[3]}\n")
                # make stats for unknown_indexes
                if f"{index1}-{index2}" not in unknown_indexes:
                    unknown_indexes[f"{index1}-{index2}"] = 1
                else:
                    unknown_indexes[f"{index1}-{index2}"] +=1
            else: # indexes are not unknown. are they hopped or do they match?
                if index1 == index2: # they match!
                    # write to R1 file
                    openoutput[f"{args.output}demux_matched_{index1}_R1.fq"].write(f"{record1[0]} {index1}-{index2}\n")
                    openoutput[f"{args.output}demux_matched_{index1}_R1.fq"].write(f"{record1[1]}\n")
                    openoutput[f"{args.output}demux_matched_{index1}_R1.fq"].write(f"{record1[2]}\n")
                    openoutput[f"{args.output}demux_matched_{index1}_R1.fq"].write(f"{record1[3]}\n")
                    # write to R2 file
                    openoutput[f"{args.output}demux_matched_{index2}_R2.fq"].write(f"{record4[0]} {index1}-{index2}\n")
                    openoutput[f"{args.output}demux_matched_{index2}_R2.fq"].write(f"{record4[1]}\n")
                    openoutput[f"{args.output}demux_matched_{index2}_R2.fq"].write(f"{record4[2]}\n")
                    openoutput[f"{args.output}demux_matched_{index2}_R2.fq"].write(f"{record4[3]}\n")
                    # make stats for matched_indexes - this is the one case where the dictionary is already initialized
                    matched_indexes[f"{index1}"] +=1
                else: # they hopped!
                    openoutput[f"{args.output}demux_hopped_R1.fq"].write(f"{record1[0]} {index1}-{index2}\n")
                    openoutput[f"{args.output}demux_hopped_R1.fq"].write(f"{record1[1]}\n")
                    openoutput[f"{args.output}demux_hopped_R1.fq"].write(f"{record1[2]}\n")
                    openoutput[f"{args.output}demux_hopped_R1.fq"].write(f"{record1[3]}\n")
                    # write to R2 file
                    openoutput[f"{args.output}demux_hopped_R2.fq"].write(f"{record4[0]} {index1}-{index2}\n")
                    openoutput[f"{args.output}demux_hopped_R2.fq"].write(f"{record4[1]}\n")
                    openoutput[f"{args.output}demux_hopped_R2.fq"].write(f"{record4[2]}\n")
                    openoutput[f"{args.output}demux_hopped_R2.fq"].write(f"{record4[3]}\n")
                    # make stats for hopped_indexes
                    if f"{index1}-{index2}" not in hopped_indexes:
                        hopped_indexes[f"{index1}-{index2}"] = 1
                    else:
                        hopped_indexes[f"{index1}-{index2}"] +=1
        # we need to break out of the loop while making sure the last record gets written to a file!
        if line1 == "" and headerstate: # we have no more lines AND we output the last line
            break
        if line1 == "": # no more lines but we need to output the last line
            headerstate = True
            continue
        # we're done analyzing the record so we make a new record in this if-else
        if headerstate:
            record1 = [line1]
            record2 = [line2]
            record3 = [line3]
            record4 = [line4]
        else: # add to records
            record1.append(line1)
            record2.append(line2)
            record3.append(line3)
            record4.append(line4)
        headerstate = False
        line_counter += 1

# close all write files
# https://stackoverflow.com/questions/29550290/how-to-open-a-list-of-files-in-python
for file in openoutput.values():
    file.close()

print("Printing stats to statsfile...", flush=True)

# make stats file
with open (f"demux_stats.txt", "w") as statsfile:
    hopped_sum = 0
    unknown_sum = 0
    matched_sum = 0
    statsfile.write(f"Hopped Indexes:\n")
    for key,value in hopped_indexes.items():
        hopped_sum += value
        statsfile.write(f"\t{key}\t{value}\n")
        statsfile.write(f"{value/(line_counter/4)*100}% of all reads\n")
    statsfile.write(f"\n{hopped_sum} hopped reads\n")
    statsfile.write(f"{hopped_sum/(line_counter/4)*100}% hopped reads\n")

    statsfile.write(f"\nMatched Indexes:\n")
    for key, value in matched_indexes.items():
        matched_sum += value
        statsfile.write(f"\t{key}\t{value}\n")
        statsfile.write(f"{value/(line_counter/4)*100}% of all reads\n")
    statsfile.write(f"\n{matched_sum} matched reads\n")
    statsfile.write(f"{matched_sum/(line_counter/4)*100}% matched reads\n")

    statsfile.write(f"\nUnknown Indexes:\n")
    for key,value in unknown_indexes.items():
        unknown_sum += value
        statsfile.write(f"\t{key}\t{value}\n")
        statsfile.write(f"{value/(line_counter/4)*100}% of all reads\n")
    statsfile.write(f"\n{unknown_sum} unknown reads\n")
    statsfile.write(f"{unknown_sum/(line_counter/4)*100}% unknown reads\n")
    
print(f"Making histograms...", flush=True)

# histogram of matched indexes
fig, ax = plt.subplots()
ax.bar(matched_indexes.keys(), matched_indexes.values(), edgecolor='k')
ax.set_title('Matched Indexes Count')
ax.minorticks_on()
ax.set_xlabel('Index')
ax.set_ylabel('Count')
fig.autofmt_xdate(rotation=45)
plt.savefig(f"matched_hist.png")

print("Done successfully!", flush=True)
