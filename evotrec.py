##################################################################################
##  evotrec.py                                                                  ##
##  (c) 2024 Andreas Ott                                                        ##
##################################################################################
"""
evotrec.py

This script performs topological recurrence analysis of genome alignments using persistent homology.

Functions:
    retrieve_metadata(timeseries_flag, refseq_id, input_afasta):
        Retrieves metadata from the input FASTA file, including the reference sequence, dates, and time range.

    murit(in_dist, out_timedist, start_date, timerange, dates):
        Converts Hamming distances to time distances and writes the results to a file.

    retrieve_snv_cycles(in_ripser):
        Analyzes cycles from Ripser output and retrieves single nucleotide variant (SNV) cycles.

    retrieve_sequences_in_cycles(snv_indices, afasta_path):
        Retrieves sequences corresponding to the SNV indices from the input FASTA file.

    retrieve_mutations_in_cycles(snv_cycles, sequences_in_snv_cycles, refseq):
        Identifies mutations in the SNV cycles and associates them with the reference sequence.

    expand_timeseries(mutation, count, timerange):
        Expands the mutation count into a time series format.

    tri_analysis(filename, timeseries_flag, timerange, mutations_in_snv_cycles):
        Computes the topological recurrence index (tRI) and writes the results to a CSV file.

Main Execution:
    Parses command-line arguments and performs the topological recurrence analysis workflow, including:
        - Retrieving metadata
        - Calculating Hamming distances
        - Converting distances to time distances (if timeseries flag is set)
        - Running Ripser for persistent homology analysis
        - Analyzing SNV cycles
        - Retrieving sequences and mutations in cycles
        - Performing tRI analysis and writing results to a CSV file
"""

import os as os
from Bio.SeqIO.FastaIO import SimpleFastaParser
from datetime import datetime
import hammingdist


def retrieve_metadata(timeseries_flag, refseq_id, input_afasta):
    """
    Retrieve metadata from a given FASTA file.

    Parameters:
    timeseries_flag (bool): Flag indicating if the sequences are part of a time series.
    refseq_id (str): The reference sequence identifier.
    input_afasta (str): Path to the input FASTA file.

    Returns:
    tuple: A tuple containing:
        - filename (str): The base name of the input FASTA file.
        - refseq (str): The reference sequence.
        - dates (list): List of dates extracted from the sequence titles.
        - start_date (str or None): The start date of the time series, if applicable.
        - end_date (str or None): The end date of the time series, if applicable.
        - timerange (int): The number of days in the time range, or 1 if not a time series.

    Raises:
    SystemExit: If an invalid date is found in the sequence titles or if the reference sequence is not found.
    """
    # Print the preparation message
    print(f"Preparing {args.afasta}...\n")

    # Initialize variables
    filename = input_afasta.split(".")[0]
    number_of_sequences = 0
    refseq = ""
    dates = []

    # If timeseries_flag is set, process the sequences as part of a time series
    if timeseries_flag:
        try:
            # Extract and validate the start date from the reference sequence ID
            start_date = refseq_id.split("|")[-2]
            datetime.strptime(start_date, "%Y-%m-%d")
        except ValueError:
            print(f"Invalid date in sequence {refseq_id}.")
            quit()

        end_date = start_date
        with open(input_afasta, "r") as fasta_file:
            for title, seq in SimpleFastaParser(fasta_file):
                number_of_sequences += 1
                if title == refseq_id:
                    refseq = seq
                try:
                    # Extract and validate the date from the sequence title
                    date = title.split("|")[-2]
                    datetime.strptime(date, "%Y-%m-%d")
                except ValueError:
                    print(f"Invalid date in sequence {title}.")
                    quit()
                # Check if the date is before the reference sequence date
                if date < start_date:
                    print(
                        f"Date before reference sequence date ({start_date}) in sequence {title}."
                    )
                    quit()
                elif date > end_date:
                    end_date = date
                dates.append(date)

        # Calculate the time range in days
        timerange = (
            datetime.fromisoformat(end_date) - datetime.fromisoformat(start_date)
        ).days + 1
        print(f"Start date: {start_date}")
        print(f"End date: {end_date}")
        print(f"Time range: {timerange} days\n")
    else:
        # If timeseries_flag is not set, just count the sequences and find the reference sequence
        with open(input_afasta, "r") as fasta_file:
            for title, seq in SimpleFastaParser(fasta_file):
                number_of_sequences += 1
                if title == refseq_id:
                    refseq = seq
        start_date = None
        end_date = None
        timerange = 1

    # Check if the reference sequence was found
    if not refseq:
        print(f"Reference sequence not found.")
        quit()

    print(f"Number of sequences: {number_of_sequences}\n")

    return filename, refseq, dates, start_date, end_date, timerange


def murit(in_dist, out_timedist, start_date, timerange, dates):
    """
    Processes a distance file and writes a time distance file based on given parameters.

    Args:
        in_dist (str): Path to the input distance file. The file should contain lines with three integers separated by spaces.
        out_timedist (str): Path to the output time distance file.
        start_date (str): The reference start date in ISO format (YYYY-MM-DD).
        timerange (int): A base value used to calculate time distances for certain conditions.
        dates (list of str): A list of dates in ISO format (YYYY-MM-DD) corresponding to the indices in the distance file.

    The function reads the input distance file line by line, processes each line to calculate a time distance based on the given parameters,
    and writes the results to the output time distance file. The time distance is calculated as follows:
        - If the distance is 1, the time distance is the maximum number of days from the reference start date for the given indices.
        - If the distance is 2, the time distance is `timerange + 1`.
        - If the distance is 3, the time distance is `timerange + 2`.
    """
    with open(in_dist, "r") as dist_file, open(out_timedist, "w") as timedist_file:
        for line in dist_file:
            # Parse the line to get the indices and distance
            i, j, dist = [int(x) for x in line.strip().split(" ")]

            # Calculate the time distance based on the given distance
            if dist == 1:
                # Calculate the number of days from the reference start date for the given indices
                days_from_refseq = [
                    (
                        datetime.fromisoformat(dates[k])
                        - datetime.fromisoformat(start_date)
                    ).days
                    + 1
                    for k in [i, j]
                ]
                # The time distance is the maximum number of days
                timedist = max(days_from_refseq)
            elif dist == 2:
                # For distance 2, the time distance is timerange + 1
                timedist = timerange + 1
            elif dist == 3:
                # For distance 3, the time distance is timerange + 2
                timedist = timerange + 2

            # Write the calculated time distance to the output file
            timedist_file.write(f"{i} {j} {timedist}\n")


def retrieve_snv_cycles(in_ripser):
    """
    Analyzes the output from Ripser to retrieve single nucleotide variant (SNV) cycles.

    This function reads a Ripser output file, identifies cycles in dimension 1, and extracts
    the cycles and their corresponding edges. It returns the cycles and the indices of the
    SNVs involved in these cycles.

    Args:
        in_ripser (str): The file path to the Ripser output file.

    Returns:
        tuple: A tuple containing:
            - snv_cycles (list): A list of cycles, where each cycle is represented as a list of edges.
              Each edge is a list containing two vertices and the length of the edge.
            - snv_indices (list): A sorted list of unique SNV indices involved in the cycles.
    """
    # Print a message indicating the start of cycle analysis
    print("Analyzing cycles...")

    # Initialize variables to store batch of cycles and individual cycle
    batch = []
    cycle = []
    go_flag = False

    # Open the Ripser output file for reading
    with open(in_ripser, "r") as file:
        for line in file:
            # Remove newline characters and strip leading/trailing whitespace
            stripped_line = line.replace("\n", " ").strip()

            # Check if the line indicates the start of persistent homology intervals in dimension 1
            if "persistent homology intervals in dim 1:" in stripped_line:
                go_flag = True

            # If the flag is set, process the lines to extract cycles
            if go_flag == True:
                if stripped_line[0] == "[":
                    # If the line starts with '[', it indicates the start of a new cycle
                    batch.append(" ".join(cycle))
                    cycle = []
                cycle.append(stripped_line)

        # Append the last cycle to the batch and remove the first empty entry
        batch.append(" ".join(cycle))
        batch.pop(0)
    # Initialize lists to store SNV cycles and a set to store unique SNV indices
    snv_cycles = []
    snv_indices = set()

    # Iterate over each raw cycle in the batch
    for raw_cycle in batch:
        zero_edge_flag = False  # Flag to check for zero-length edges
        interval_string = raw_cycle.split(":")[0]  # Extract interval string
        cycle = []  # List to store edges and their lengths for the current cycle
        edges = []  # List to store edges for the current cycle

        # Process each piece of the raw cycle to extract edges and their lengths
        for piece in raw_cycle.split("{")[-1].split("[")[1:]:
            edge = sorted(
                [int(piece.split(",")[0]), int(piece.split(",")[1].split("]")[0])]
            )
            length = int(piece.split(")")[0].split("(")[-1])
            cycle.append([edge, length])
            edges.append(edge)
            if length == 0:
                zero_edge_flag = True  # Set flag if a zero-length edge is found
                break

        # If no zero-length edge is found, add the cycle and update SNV indices
        if zero_edge_flag is False:
            snv_cycles.append([[[pair[0][0], pair[0][1]], pair[1]] for pair in cycle])
            for edge in edges:
                snv_indices.update(edge)

    # Convert the set of SNV indices to a sorted list
    snv_indices = sorted(list(snv_indices))

    return snv_cycles, snv_indices


def retrieve_sequences_in_cycles(snv_indices, afasta_path):
    """
    Retrieve sequences from a FASTA file at specified indices.

    This function reads a FASTA file and retrieves sequences at the positions
    specified in the `snv_indices` list. The sequences are returned in a
    dictionary where the keys are the indices and the values are the sequences.

    Args:
        snv_indices (list of int): A list of indices specifying which sequences
                                   to retrieve from the FASTA file.
        afasta_path (str): The path to the FASTA file.

    Returns:
        dict: A dictionary where the keys are the indices from `snv_indices`
              and the values are the corresponding sequences from the FASTA file.
    """
    sequences = {}
    with open(afasta_path, "r") as file:
        for i, (_, seq) in enumerate(SimpleFastaParser(file)):
            try:
                # Check if the current index matches the first index in snv_indices
                if i == snv_indices[0]:
                    sequences[i] = seq  # Store the sequence in the dictionary
                    snv_indices.pop(0)  # Remove the first index from snv_indices
            except:
                break

    return sequences


def retrieve_mutations_in_cycles(snv_cycles, sequences_in_snv_cycles, refseq):
    """
    Retrieve mutations in cycles of single nucleotide variants (SNVs).

    This function identifies mutations in cycles of SNVs by comparing sequences
    in each cycle to a reference sequence. It returns a list of mutations for
    each cycle, including the site of mutation, the original and mutated bases,
    the edge in the cycle, and the maximum length of the edges in the cycle.

    Args:
        snv_cycles (list): A list of cycles, where each cycle is a list of tuples.
                           Each tuple contains an edge (a pair of vertices) and a length.
        sequences_in_snv_cycles (dict): A dictionary mapping vertices to their corresponding sequences.
        refseq (str): The reference sequence to compare against.

    Returns:
        list: A list of lists, where each inner list contains mutations for a cycle.
              Each mutation is represented as a list containing the site of mutation,
              the original base, the mutated base, the edge, and the maximum length of the edges in the cycle.
    """
    mutations_in_snv_cycles = []

    for cycle in snv_cycles:
        mutations_per_cycle = []

        for edge, length in cycle:
            sequence_pair = []

            # Collect sequences for the vertices in the edge
            for vertex in edge:
                sequence_pair.append(sequences_in_snv_cycles[vertex])

            # Compare sequences site by site
            for site, pair in enumerate(
                zip(*sequence_pair), 1
            ):  # count of sites is 1-indexed
                if (pair[0] != pair[1]) and (
                    "-" not in pair
                ):  # Check for mutations and ignore gaps
                    if pair[0] == refseq[site - 1]:
                        mutations_per_cycle.append(
                            [(site, pair[0], pair[1]), edge, length]
                        )
                    elif pair[1] == refseq[site - 1]:
                        mutations_per_cycle.append(
                            [(site, pair[1], pair[0]), edge, length]
                        )

        lengths = [length for mutation, edge, length in mutations_per_cycle]

        # Determine the maximum length in the cycle
        max_length = max(lengths) if lengths else 0

        # Update mutations with the maximum length
        mutations_per_cycle_max_length = [
            [mutation, edge, max_length]
            for mutation, edge, length in mutations_per_cycle
        ]

        # Append the processed mutations for the current cycle
        mutations_in_snv_cycles.append(mutations_per_cycle_max_length)

    return mutations_in_snv_cycles


def expand_timeseries(mutation, count, timerange):
    """
    Expands a timeseries based on mutation data and count over a specified time range.

    Args:
        mutation (list): A list containing mutation data. The first element is expected to be a string,
                         and the second and third elements are expected to be integers.
        count (list): A list of integers representing counts for each day.
        timerange (int): The number of days to expand the timeseries over.

    Returns:
        str: A comma-separated string representing the expanded timeseries.
    """
    # Initialize the timeseries list with mutation details
    tri_timeseries = [str(mutation[0]), mutation[1], mutation[2]]
    tri = 0

    # Iterate over the range of days
    for day in range(timerange):
        try:
            tri = tri + count[day + 1]
        except:
            pass

        # Append the current value of tri to the timeseries
        tri_timeseries.append(str(tri))

    return ",".join(tri_timeseries)


def tri_analysis(filename, timeseries_flag, timerange, mutations_in_snv_cycles):
    """
    Analyzes temporal recurrence index (tRI) from mutation data and writes the results to a CSV file.

    Parameters:
    filename (str): The base name of the output CSV file.
    timeseries_flag (bool): If True, the output CSV will include a timeseries header.
    timerange (int): The range of time points to consider in the analysis.
    mutations_in_snv_cycles (list): A list of mutation cycles, where each cycle is a list of tuples.
                                     Each tuple contains (mutation, edge, day).
    """
    print("Computing tRI...\n")
    tri_count = {}
    noted_edges = []

    # Iterate through each set of mutations in the cycles
    for mutations in mutations_in_snv_cycles:
        noted_snvs_in_cycle = []

        # Iterate through each mutation, edge, and day in the mutations
        for (
            mutation,
            edge,
            day,
        ) in mutations:  # edge length is time difference to start date + 1
            if (mutation not in noted_snvs_in_cycle) and (edge not in noted_edges):
                tri_count.setdefault(mutation, {})
                tri_count[mutation].setdefault(day, 0)
                tri_count[mutation][
                    day
                ] += 1  # Increment the count for the mutation on the given day
                noted_snvs_in_cycle.append(mutation)
                noted_edges.append(edge)

    # Write the results to a CSV file
    with open(f"{filename}.csv", "w") as out_csv:
        if timeseries_flag:
            # Write the header for timeseries data
            out_csv.write(
                ",".join(
                    ["POS", "REF", "ALT"] + [str(i) for i in range(1, timerange + 1)]
                )
                + "\n"
            )
        else:
            # Write the header for non-timeseries data
            out_csv.write("POS,REF,ALT,TRI\n")

        # Write the mutation data to the CSV file
        for mutation, count in tri_count.items():
            out_csv.write(expand_timeseries(mutation, count, timerange) + "\n")

    # Print a message based on whether any tRI signal was found
    if len(tri_count) == 0:
        print("No tRI signal found.")
    else:
        print(f"Results written to {filename}.csv.")


if __name__ == "__main__":

    import argparse
    import textwrap

    # Set up argument parser
    parser = argparse.ArgumentParser(exit_on_error=False)
    parser.add_argument("afasta")
    parser.add_argument("refseq_id")
    parser.add_argument("--timeseries", nargs="?", const=True, default=False)
    args, unknown_args = parser.parse_known_args()

    # Check for unknown arguments
    if len(unknown_args) > 0:
        print("Invalid arguments.")
        quit()

    # Validate arguments
    try:
        str(args.afasta)
        str(args.refseq_id)
        str(args.timeseries)
    except ValueError:
        print("Invalid arguments.")
        quit()

    # Print program header    
    print(textwrap.dedent("""
        ===============================================================
        EVOtRec -- Topological recurrence analysis of genome alignments
        (c) 2024 Andreas Ott
        ===============================================================
        """)
    )

    # Retrieve metadata from input files
    filename, refseq, dates, start_date, end_date, timerange = retrieve_metadata(
        timeseries_flag=args.timeseries,
        refseq_id=args.refseq_id,
        input_afasta=args.afasta,
    )

    # Perform Hamming distance calculation
    print("Hammingdist...")
    data = hammingdist.from_fasta(args.afasta, max_distance=4)
    data.dump_sparse(filename + ".dist", threshold=3)

    # Perform MuRiT and Ripser analysis based on timeseries flag
    if args.timeseries != False:
        print("MuRiT...")
        murit(
            in_dist=filename + ".dist",
            out_timedist=filename + ".timedist",
            start_date=start_date,
            timerange=timerange,
            dates=dates,
        )
        print("Ripser...")
        os.system(
            f"ripser --dim 1 --format sparse --threshold {timerange+1} {filename + '.timedist'} > {filename + '.ripser'}"
        )
    else:
        print("Ripser...")
        os.system(
            f"ripser --dim 1 --format sparse --threshold 2 {filename + '.dist'} > {filename + '.ripser'}"
        )

    # Retrieve SNV cycles and sequences
    snv_cycles, snv_indices = retrieve_snv_cycles(in_ripser=filename + ".ripser")

    sequences_in_cycles = retrieve_sequences_in_cycles(
        snv_indices=snv_indices, afasta_path=args.afasta
    )

    # Retrieve mutations in cycles
    mutations_in_snv_cycles = retrieve_mutations_in_cycles(
        snv_cycles=snv_cycles,
        sequences_in_snv_cycles=sequences_in_cycles,
        refseq=refseq,
    )

    # Perform tRI analysis and write results to a CSV file
    tri_analysis(
        filename=filename,
        timeseries_flag=args.timeseries,
        timerange=timerange,
        mutations_in_snv_cycles=mutations_in_snv_cycles,
    )
