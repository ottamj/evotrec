##################################################################################
##  evotrec.py                                                                  ##
##  (c) 2024 Andreas Ott                                                        ##
##################################################################################

import os as os
from Bio.SeqIO.FastaIO import SimpleFastaParser
from datetime import datetime
import hammingdist

def retrieve_metadata(
    timeseries_flag,
    refseq_id,
    input_afasta
):
    print(f"Preparing {args.afasta}...\n")
    filename = input_afasta.split('.')[0]
    number_of_sequences = 0
    refseq = ""
    dates = []
    if timeseries_flag != False:
        try:
            start_date = refseq_id.split('|')[-2]
            datetime.strptime(start_date, '%Y-%m-%d')
        except:
            print(f"Invalid date in sequence {refseq_id}.")
            quit()
        end_date = start_date
        with open(input_afasta, 'r') as fasta_file:
            for title, seq in SimpleFastaParser(fasta_file):
                number_of_sequences += 1
                if title == refseq_id:
                    refseq = seq
                try:
                    date = title.split('|')[-2]
                    datetime.strptime(date, '%Y-%m-%d')
                except:
                    print(f"Invalid date in sequence {title}.")
                    quit()
                if date < start_date:
                    print(f"Date before reference sequence date ({start_date}) in sequence {title}.")
                    quit()
                elif date > end_date:
                    end_date = date
                dates.append(date)
        timerange = (datetime.fromisoformat(end_date) - datetime.fromisoformat(start_date)).days + 1
        print(f"Start date: {start_date}")
        print(f"End date: {end_date}")
        print(f"Time range: {timerange} days\n")
    else:
        with open(input_afasta, 'r') as fasta_file:
            for title, seq in SimpleFastaParser(fasta_file):
                number_of_sequences += 1
                if title == refseq_id:
                    refseq = seq
        start_date = None
        end_date = None
        timerange = 1
    if refseq == "":
        print(f"Reference sequence not found.")
        quit()
    print(f"Number of sequences: {number_of_sequences}\n")

    return filename, refseq, dates, start_date, end_date, timerange


def murit(
    in_dist,
    out_timedist,
    start_date,
    timerange,
    dates
):
    with open(in_dist, "r") as dist_file, open(out_timedist, "w") as timedist_file:
        for line in dist_file:
            i, j, dist = [int(x) for x in line.strip().split(" ")]
            if dist == 1:
                days_from_refseq = [(datetime.fromisoformat(dates[k]) - datetime.fromisoformat(start_date)).days + 1 for k in [i,j]]
                timedist = max(days_from_refseq)
            elif dist == 2:
                timedist = timerange + 1
            elif dist == 3:
                timedist = timerange + 2
            timedist_file.write(f"{i} {j} {timedist}\n")


def retrieve_snv_cycles(
    in_ripser
):
    print('Analyzing cycles...')
    batch = []
    cycle = []
    go_flag = False
    with open(in_ripser, "r") as file:
        for line in file:
            stripped_line = line.replace("\n", " ").strip()
            if 'persistent homology intervals in dim 1:' in stripped_line:
                go_flag = True
            if go_flag == True:
                if stripped_line[0] == "[":
                    batch.append(' '.join(cycle))
                    cycle = []
                cycle.append(stripped_line)
        batch.append(' '.join(cycle))
        batch.pop(0)
    snv_cycles = []
    snv_indices = set()
    for raw_cycle in batch:
        zero_edge_flag = False
        interval_string = raw_cycle.split(':')[0]
        cycle = []
        edges = []
        for piece in raw_cycle.split('{')[-1].split('[')[1:]: #select exhaustive cycles (2nd row in Ripser output)
            edge = sorted([int(piece.split(',')[0]), int(piece.split(',')[1].split(']')[0])])
            length = int(piece.split(')')[0].split('(')[-1])
            cycle.append([edge, length])
            edges.append(edge)
            if length == 0:
                zero_edge_flag = True
                break
        if zero_edge_flag is False:
            snv_cycles.append([[[pair[0][0], pair[0][1]], pair[1]] for pair in cycle])
            for edge in edges:
                snv_indices.update(edge)
    snv_indices = sorted(list(snv_indices))

    return snv_cycles, snv_indices


def retrieve_sequences_in_cycles(
    snv_indices,
    afasta_path
):
    sequences = {}
    with open(afasta_path, 'r') as file:
        for i, (title, seq) in enumerate(SimpleFastaParser(file)):
            try:
                if i == snv_indices[0]:
                    sequences[i]=seq
                    snv_indices.pop(0)
            except:
                break

    return sequences


def retrieve_mutations_in_cycles(
        snv_cycles,
        sequences_in_snv_cycles,
        refseq
):
    mutations_in_snv_cycles = []
    for cycle in snv_cycles:
        mutations_per_cycle = []
        for edge, length in cycle:
            sequence_pair = []
            for vertex in edge:
                sequence_pair.append(sequences_in_snv_cycles[vertex])
            for site, pair in enumerate(zip(*sequence_pair), 1): #count of sites starts at 1 (ie count is 1-indexed)
                if (pair[0] != pair[1]) and ('-' not in pair):
                    if pair[0] == refseq[site-1]:
                        mutations_per_cycle.append([(site, pair[0], pair[1]), edge, length])
                    elif pair[1] == refseq[site-1]:
                        mutations_per_cycle.append([(site, pair[1], pair[0]), edge, length])
        lengths = []
        for mutation, edge, length in mutations_per_cycle:
            lengths.append(length)
        if len(lengths) > 0:
            max_length = max(lengths)
        else:
            max_length = 0
        mutations_per_cycle_max_length = [[mutation, edge, max_length] for mutation, edge, length in mutations_per_cycle]
        mutations_in_snv_cycles.append(mutations_per_cycle_max_length)

    return mutations_in_snv_cycles


def expand_timeseries(mutation, count, timerange):
    tri_timeseries = [str(mutation[0]), mutation[1], mutation[2]]
    tri = 0
    for day in range(timerange):
        try: tri = tri + count[day+1]
        except: pass
        tri_timeseries.append(str(tri))

    return ','.join(tri_timeseries)


def tri_analysis(
    filename,
    timeseries_flag,
    timerange,
    mutations_in_snv_cycles
):
    print('Computing tRI...\n')
    tri_count = {}
    cycles_count = {}
    noted_edges = []
    for mutations in mutations_in_snv_cycles:
        noted_snvs_in_cycle = []
        for mutation, edge, day in mutations: #edge length is time difference to start date + 1
            if (mutation not in noted_snvs_in_cycle) and (edge not in noted_edges):
                tri_count.setdefault(mutation, {})
                tri_count[mutation].setdefault(day, 0)
                tri_count[mutation][day] += 1
                noted_snvs_in_cycle.append(mutation)
                noted_edges.append(edge)
    with open(f"{filename}.csv", 'w') as out_csv:
        if timeseries_flag == True:
            out_csv.write(','.join(['POS', 'REF', 'ALT'] + [str(i) for i in range(1, timerange + 1)])+'\n')
        else:
            out_csv.write('POS,REF,ALT,TRI\n')
        for mutation, count in tri_count.items():
            out_csv.write(expand_timeseries(mutation, count, timerange)+"\n")
    if len(tri_count) == 0:
        print('No tRI signal found.')
    else:
        print(f"Results written to {filename}.csv.")


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(exit_on_error=False)
    parser.add_argument("afasta")
    parser.add_argument("refseq_id")
    parser.add_argument("--timeseries", nargs='?', const=True, default=False)
    args, unknown_args = parser.parse_known_args()
    if len(unknown_args) > 0:
        print("Invalid arguments.")
        quit()
    try:
        str(args.afasta)
        str(args.refseq_id)
        str(args.timeseries)
    except ValueError:
        print("Invalid arguments.")
        quit()

    print('''
===============================================================
EVOtRec -- Topological recurrence analysis of genome alignments
(c) 2024 Andreas Ott
===============================================================
''')

    filename, refseq, dates, start_date, end_date, timerange = retrieve_metadata(
        timeseries_flag=args.timeseries,
        refseq_id=args.refseq_id,
        input_afasta=args.afasta
    )

    print('Hammingdist...')
    data = hammingdist.from_fasta(args.afasta, max_distance=4)
    data.dump_sparse(filename+'.dist', threshold=3)

    if args.timeseries != False:
        print('MuRiT...')
        murit(in_dist=filename+'.dist', out_timedist=filename+'.timedist', start_date=start_date, timerange=timerange, dates=dates)
        print('Ripser...')
        os.system(f"ripser --dim 1 --format sparse --threshold {timerange+1} {filename + '.timedist'} > {filename + '.ripser'}")
    else:
        print('Ripser...')
        os.system(f"ripser --dim 1 --format sparse --threshold 2 {filename + '.dist'} > {filename + '.ripser'}")

    snv_cycles, snv_indices = retrieve_snv_cycles(
        in_ripser=filename+'.ripser'
    )

    sequences_in_cycles = retrieve_sequences_in_cycles(
        snv_indices=snv_indices,
        afasta_path=args.afasta
    )

    mutations_in_snv_cycles = retrieve_mutations_in_cycles(
        snv_cycles=snv_cycles,
        sequences_in_snv_cycles=sequences_in_cycles,
        refseq=refseq
    )

    tri_analysis(
        filename=filename,
        timeseries_flag=args.timeseries,
        timerange=timerange,
        mutations_in_snv_cycles=mutations_in_snv_cycles
    )
