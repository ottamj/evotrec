##################################################################################
##  evotrec.py                                                                  ##
##  (c) 2024 Andreas Ott                                                        ##
##################################################################################
"""
evotrec.py

This script performs topological recurrence analysis of genome alignments using
persistent homology.

References:
    - http://arxiv.org/abs/2106.07292
    - https://arxiv.org/abs/2207.03394


Functions:
    retrieve_metadata(timeseries_flag, refseq_id, input_afasta):
        Retrieves metadata from the input FASTA file, including the reference sequence,
        dates, and time range.

    murit(in_dist, out_timedist, start_date, timerange, dates):
        Converts Hamming distances to Rips transformed distances and writes the results
        to a file.

    retrieve_snv_cycles(in_ripser):
        Analyzes cycles from Ripser output and retrieves single nucleotide variant (SNV)
        cycles.

    retrieve_sequences_in_cycles(snv_indices, afasta_path):
        Retrieves sequences corresponding to the SNV indices from the input FASTA file.

    retrieve_mutations_in_cycles(snv_cycles, sequences_in_snv_cycles, refseq):
        Identifies mutations in the SNV cycles and associates them with the reference
        sequence.

    expand_timeseries(mutation, count, timerange):
        Expands the mutation count into a time series format.

    tri_analysis(filename, timeseries_flag, timerange, mutations_in_snv_cycles):
        Computes the topological recurrence index (tRI) and writes the results to a CSV
        file.

Main Execution:
    Parses command-line arguments and performs the topological recurrence analysis
    workflow:
        - Retrieving metadata
        - Calculating Hamming distances
        - Converting distances to Rips transformed distances (if timeseries flag is set)
        - Running Ripser (with representatives output) for persistent homology analysis
        - Analyzing SNV cycles
        - Retrieving sequences and mutations in cycles
        - Performing tRI analysis and writing results to a CSV file
"""
# Max 2025-02-10
# 1. Meine neuen Analysen zeigen Site-Patterns >2 (also nicht nur Z/2Z sondern auch Z/3Z, Z/4Z)
#   -> Entsprechend sollte man denke ich in tri_analysis(...) einen Codeblock ergänzen,
#   der diese Fälle auschließt.
#   -> Ich habe einen Vorschlag eingefügt (auskommentiert)

# 2. Man sieht außerdem Site-Patterns =1
#   -> kann man lösen indem man SNV-Zykeln betrachtet bei denen alle Sequenzen im Zykel an den im
#   Zykel beteiligten Stellen keinen Dash haben
#   -> entsprechende Code-Blöcke habe ich in retrieve_mutations_in_cycles(...) ergänzt.

# 3. translation to refseq positions (siehe Ende von main())

# 4. Bei Anpassung von tRI-Analyse muss man eventuell auch expand_timeseries (siehe comments in tri_analysis())

# 5. Übergang zu Tupeln für edges: [[a,b],dist(a,b)] -> ((a,b),dist(a,b)) besseres Handling mit set() und weniger Speicherplatz

# 6. Keine statische Analyse mehr: Stattdessen Datum SLICE_DATE = 'YYYY-MM-DD'<= end_date -> tRI_ouptut:
#   tRI_analysis_{ref_seq_date}_to_end_date
#   tRI_analysis_slice_date = slice von tRI_analysis_slice_date

# 7. Um Lesbarkeit von Code zu verbessern, evtl. aliases einfügen. Code-Beispiel:
#   Edge = Tuple[int, int]               # Represents an edge between two nodes (i, j)
#   EdgeWithLength = Tuple[Edge, int]    # Represents an edge with its associated length ((i, j), length)
#   Cycle = List[EdgeWithLength]         # Represents a cycle as a list of edges with lengths
#   Triangle = List[EdgeWithLength]      # Represents a triangle as a list of three edges with lengths
#   Funktionsbeispiel:
#   def identify_and_remove_noted_triangles(
#       cycles: List[Cycle],
#       triangles: List[Triangle]
#   ) -> List[Cycle]:
#   """
#   DOC_STRING
#   """
#   ...CODE...
#   return ...
# 
# Michael 2025-08-06
# 5. erledigt
# 7. erledigt

import os
from datetime import datetime
from typing import List, Tuple, Dict, Optional

from Bio.SeqIO.FastaIO import SimpleFastaParser
import hammingdist


# Type aliases for improved code readability

# Represents an edge between two nodes (i, j)
Edge = Tuple[int, int]

# Represents an edge with its associated length ((i, j), length)
EdgeWithLength = Tuple[Edge, int]

# Represents a cycle as a list of edges with lengths
Cycle = List[EdgeWithLength]

# Represents a triangle as a list of three edges with lengths
Triangle = List[EdgeWithLength]

# Represents a mutation (position, reference_base, alternate_base)
Mutation = Tuple[int, str, str]

# Mutation with edge and length information
MutationEntry = List[Tuple[Mutation, Edge, int]]


def validate_date(date_str: str, context: str) -> None:
    """
    Validates a date string in the format 'YYYY-MM-DD'
    and raises a ValueError with context if invalid.
    """
    try:
        datetime.strptime(date_str, "%Y-%m-%d")
    except ValueError as e:
        raise ValueError(f"Invalid date in {context}: {e}")


def retrieve_metadata(
    input_afasta: str, refseq_id: str, timeseries_flag: bool
) -> Tuple[str, str, List[str], Optional[str], Optional[str], int]:
    """
    Retrieve metadata from a given FASTA file.

    Parameters:
    input_afasta (str): Path to the input FASTA file.
    refseq_id (str): The reference sequence identifier.
    timeseries_flag (bool): Flag indicating if the sequences are part of a time series.

    Returns:
    Tuple[str, str, List[str], Optional[str], Optional[str], int]: A tuple containing:
        - filename (str): The base name of the input FASTA file.
        - refseq (str): The reference sequence.
        - dates (List[str]): List of dates extracted from the sequence titles.
        - start_date (Optional[str]): The start date of the time series, if applicable.
        - end_date (Optional[str]): The end date of the time series, if applicable.
        - timerange (int): The number of days in the time range or 1 if not time series.

    Raises:
        ValueError: If an invalid date is found in the sequence titles
        or if the reference sequence is not found.
    """
    print(f"Preparing {input_afasta}...\n")

    filename = input_afasta.split(".")[0]
    number_of_sequences = 0
    refseq = ""
    dates = []

    if timeseries_flag:
        start_date = refseq_id.split("|")[-2]
        validate_date(start_date, f"sequence {refseq_id}")

        end_date = start_date
        with open(input_afasta, "r") as fasta_file:
            for title, seq in SimpleFastaParser(fasta_file):
                number_of_sequences += 1
                if title == refseq_id:
                    refseq = seq

                date = title.split("|")[-2]
                validate_date(date, f"sequence {title}")

                if date < start_date:
                    raise ValueError(
                        f"Date before reference sequence date ({start_date}) "
                        f"in sequence {title}."
                    )
                elif date > end_date:
                    end_date = date

                dates.append(date)

        timerange = (
            datetime.fromisoformat(end_date) - datetime.fromisoformat(start_date)
        ).days + 1
        print(f"Start date: {start_date}")
        print(f"End date: {end_date}")
        print(f"Time range: {timerange} days\n")
    else:
        with open(input_afasta, "r") as fasta_file:
            for title, seq in SimpleFastaParser(fasta_file):
                number_of_sequences += 1
                if title == refseq_id:
                    refseq = seq
        start_date = None
        end_date = None
        timerange = 1

    if not refseq:
        raise ValueError("Reference sequence not found.")

    print(f"Number of sequences: {number_of_sequences}\n")

    return filename, refseq, dates, start_date, end_date, timerange


def murit(
    input_distance_file: str,
    output_distance_file: str,
    start_date: str,
    timerange: int,
    dates: List[str],
) -> None:
    """
    Implements a simplified version of the Rips transformation, outlined in the
    referenced article, for sequence alignments with temporal information.

    The function reads the input distance file line by line, processes each line to
    calculate a Rips transformed distance based on the given parameters, and writes the
    results to the output file. The transformed distance is calculated as follows:
        - If the distance is 1, the transformed distance is the maximum number of days
            from the reference start date for the given indices.
        - If the distance is 2, the transformed distance is `timerange + 1`.
        - If the distance is 3, the transformed distance is `timerange + 2`.

    References:
        - https://arxiv.org/abs/2207.03394

    Args:
        input_distance_file (str): Path to the input distance file.
            The file should contain lines with three integers separated by spaces,
            representing two indices and a distance (1, 2, or 3).
        output_distance_file (str): Path to the output file.
        start_date (str): The reference start date in ISO format (YYYY-MM-DD)
            used to compute the temporal distance in days.
        timerange (int): The base value added to distances 2 or larger.
        dates (List[str]): A list of dates in ISO format (YYYY-MM-DD), where
             each date corresponds to an index in the distance file.

    Output:
        The result is written to `output_distance_file`
        Each line contains three integers: two indices and the transformed distance.
    """
    with open(input_distance_file, "r") as dist_file, open(
        output_distance_file, "w"
    ) as timedist_file:
        for line in dist_file:
            i, j, dist = [int(x) for x in line.strip().split(" ")]

            # Calculate the transformed distance based on the given distance
            if dist == 1:
                # The transformed distance is the maximum number of days
                days_from_refseq = [
                    (
                        datetime.fromisoformat(dates[k])
                        - datetime.fromisoformat(start_date)
                    ).days
                    + 1
                    for k in [i, j]
                ]
                timedist = max(days_from_refseq)
            else:
                # For distance 2 or larger, add the timerange to the distance
                timedist = timerange + dist - 1

            # Write the calculated distance to the output file
            timedist_file.write(f"{i} {j} {timedist}\n")


def retrieve_snv_cycles(in_ripser: str) -> Tuple[List[Cycle], List[int]]:
    """
    Analyzes the output from Ripser to retrieve single nucleotide variant (SNV) cycles.

    This function reads a Ripser output file, identifies cycles in dimension 1, and
    extracts the exhaustive cycle representatives. It returns the cycle representatives
    and the indices of the sequences involved in these cycles.

    Args:
        in_ripser (str): The file path to the Ripser output file.

    Returns:
        Tuple[List[Cycle], List[int]]: A tuple containing:
            - snv_cycles (List[Cycle]): A list of cycles, where each cycle is represented
                as a list of edges with their lengths.
            - snv_indices (List[int]): A sorted list of sequence indices involved in cycles.

    Example:
        The input file should have the following structure:

        ```
        ... (other lines)
        persistent homology intervals in dim 0:
        [0,1):  {[2], [3]}
        ... (other lines)
        [0, ):  {[3]:1}
        persistent homology intervals in dim 1:
        [1,2):
        {[0,1] (1), [0,2] (1), [1,3] (1), [2,3] (1)}
        {[0,1] (1), [0,2] (1), [1,3] (1), [2,3] (1)}
        [1,2):
        TODO (edit Max)
        [1,2):
        {[4101,7675] (1), [6638,7675] (1), [4101,7803] (1), [6638,7803] (1)}
        {[4101,7675] (1), [6638,7675] (1), [4101,7803] (1), [6638,7803] (1)}
        [1,2):
        {[3871,7691] (1), [4993,7691] (1), [3483,7800] (1), [4993,7800] (1), [3483,7801] (1), [3871,7801] (1)}
        {[3871,7691] (1), [4993,7691] (1), [3483,7800] (1), [4993,7800] (1), [3483,7801] (1), [3871,7801] (1)}
        ... (other lines)
        ```

        Each interval [b,d) is a persistence interval, each {...} its representative.
        TODO (edit Max): a persistence interval with two representatives, where the first/second one is the exhaustive representative.
        The lines under "persistence intervals in dim 1:" are the ones that are
        processed to extract the cycle representatives.
    """
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

            # Set flag to start processing the lines for 1d cycles
            if "persistent homology intervals in dim 1:" in stripped_line:
                go_flag = True

            # If the flag is set, process the lines to extract cycles
            if go_flag is True:
                if stripped_line[0] == "[":
                    # If the line starts with '[' it's the start of a new cycle
                    batch.append(" ".join(cycle))
                    cycle = []
                cycle.append(stripped_line)

        # Append the last cycle to the batch
        batch.append(" ".join(cycle))
        # remove the first entry: 'persistent homology intervals in dim 1:'
        batch.pop(0)

    # Initialize lists to store SNV cycles and a set to store unique SNV indices
    snv_cycles = []
    snv_indices = set()

    # Iterate over each raw cycle in the batch
    for raw_cycle in batch:
        zero_edge_flag = False  # Flag to check for zero-length edges
        cycle = []  # List to store edges and their lengths for the current cycle
        edges = []  # List to store edges for the current cycle

        # Process each piece of the raw cycle to extract edges and their lengths
        for piece in raw_cycle.split("{")[-1].split("[")[1:]:
            edge = tuple(sorted(
                [int(piece.split(",")[0]), int(piece.split(",")[1].split("]")[0])]
            ))
            length = int(piece.split(")")[0].split("(")[-1])
            cycle.append((edge, length))
            edges.append(edge)
            if length == 0:
                zero_edge_flag = True  # Set flag if a zero-length edge is found
                break
        # TODO (edit Max): maybe remove -> we have no duplicates
        # If no zero-length edge is found, add the cycle and update SNV indices
        if zero_edge_flag is False:
            snv_cycles.append([((pair[0][0], pair[0][1]), pair[1]) for pair in cycle])
            for edge in edges:
                snv_indices.update(edge)

    # Convert the set of SNV indices to a sorted list
    snv_indices = sorted(list(snv_indices))

    return snv_cycles, snv_indices


def retrieve_sequences_in_cycles(
    snv_indices: List[int], afasta_path: str
) -> Dict[int, str]:
    """
    Retrieve sequences from a FASTA file at specified indices.

    This function reads a FASTA file and retrieves sequences at the positions
    specified in the `snv_indices` list. The sequences are returned in a
    dictionary where the keys are the indices and the values are the sequences.

    Args:
        snv_indices (List[int]): An increasing list of indices specifying which
                sequences to retrieve from the FASTA file.
        afasta_path (str): The path to the FASTA file.

    Returns:
        Dict[int, str]: A dictionary where the keys are the indices from `snv_indices`
              and the values are the corresponding sequences from the FASTA file.
    """
    sequences = {}
    with open(afasta_path, "r") as file:
        for i, (_, seq) in enumerate(SimpleFastaParser(file)):
            if not snv_indices:
                break  # Exit the loop if snv_indices is empty
            if i == snv_indices[0]:
                sequences[i] = seq  # Store the sequence in the dictionary
                snv_indices.pop(0)  # Remove the first index from snv_indices

    return sequences


def retrieve_mutations_in_cycles(
    snv_cycles: List[Cycle], sequences_in_snv_cycles: Dict[int, str], refseq: str
) -> List[MutationEntry]:
    """
    Retrieve mutations in SNV cycles.

    This function identifies mutations in SNV cycles by comparing sequences
    in each cycle to a reference sequence. It returns a list of mutations for
    each cycle, including the site of mutation, the original and mutated bases,
    the edge in the cycle, and the maximum length of the edges in the cycle.

    Args:
        snv_cycles (List[Cycle]): A list of cycles, where each cycle is a list of edges with lengths.
        sequences_in_snv_cycles (Dict[int, str]): A dictionary mapping vertices to their
            corresponding sequences.
        refseq (str): The reference sequence to compare against.

    Returns:
        List[MutationEntry]: A list of lists, where each inner list collects the mutations that
            appear in a cycle. Each mutation is represented as a list
            of the form [(POS, REF, ALT), edge, max_length].
    """
    mutations_in_snv_cycles = []

    for cycle in snv_cycles:
        mutations_per_cycle = []
        # Comment Max: Handling of dash sites
        # involved_sites = []
        # involved_vertices = []

        for edge, length in cycle:
            sequence_pair = []

            # Collect sequences for the vertices in the edge
            for vertex in edge:
                sequence_pair.append(sequences_in_snv_cycles[vertex])
                # Comment Max: Handling of dash sites
                # involved_vertices.append(vertex)

            # Compare sequences site by site
            for site, pair in enumerate(
                zip(*sequence_pair), 1  # count of sites is 1-indexed
            ):
                if (pair[0] != pair[1]) and ("-" not in pair):
                    # Check for mutations and ignore gaps
                    if pair[0] == refseq[site - 1]:
                        mutations_per_cycle.append(
                            ((site, pair[0], pair[1]), edge, length)
                        )
                        # Comment Max: Handling of dash sites
                        # involved_sites.append(site)
                    elif pair[1] == refseq[site - 1]:
                        mutations_per_cycle.append(
                            ((site, pair[1], pair[0]), edge, length)
                        )
                        # Comment Max: Handling of dash sites
                        # involved_sites.append(site)

            # Comment Max: Handling of dash sites
            # involved_sites = list(set(involved_sites))
            # involved_vertices = list(set(involved_vertices))

        # dash_flag = True
        # for vertex in involved_vertices:
        #     for site in involved_sites:
        #         if sequences_in_snv_cycles[vertex][site - 1] == "-":
        #             sequences_in_snv_cycles[vertex][site - 1]
        #             dash_flag = False
        # if dash_flag == False:
        #    continue

        lengths = [length for _, _, length in mutations_per_cycle]

        # Determine the maximum length in the cycle
        max_length = max(lengths) if lengths else 0

        # Update mutations with the maximum length
        mutations_per_cycle_max_length = [
            (mutation, edge, max_length) for mutation, edge, _ in mutations_per_cycle
        ]

        # Append the processed mutations for the current cycle
        mutations_in_snv_cycles.append(mutations_per_cycle_max_length)

    return mutations_in_snv_cycles


def expand_timeseries(mutation: Mutation, count: Dict[int, int], timerange: int) -> str:
    """
    Helper function to generate a string representation of tRI timeseries data, used
    in writing a CSV file.

    Args:
        mutation (Mutation): A tuple representing a mutation in the form (POS, REF, ALT)
        count (Dict[int, int]): A dictionary mapping days to tRI counts for each day.
        timerange (int): The number of days to expand the timeseries over.

    Returns:
        str: A comma-separated string representing the expanded timeseries.
    """
    # Initialize the timeseries list with mutation details
    tri_timeseries = [str(mutation[0]), mutation[1], mutation[2]]
    tri_value = 0

    # Iterate over the range of days
    for day in range(timerange):
        tri_value += count.get(day + 1, 0)

        # Append the current value of tri to the timeseries
        tri_timeseries.append(str(tri_value))

    return ",".join(tri_timeseries)


def tri_analysis(
    mutations_in_snv_cycles: List[MutationEntry],
    timerange: int,
    output_filename: str,
    timeseries_flag: bool,
) -> None:
    """
    Performs tRI analysis from mutation data and writes the results to a CSV file.

    Parameters:
    mutations_in_snv_cycles (List[MutationEntry]): A list of mutation cycles, where each cycle is
        a list of the form [(POS, REF, ALT), edge, max_length].
    timerange (int): The range of time points to consider in the analysis.
    output_filename (str): The base name of the output CSV file (.csv extension will be added)
    timeseries_flag (bool): If True, output CSV will capture the tRI timeseries data.
    """
    print("Computing tRI...\n")
    tri_count = {}
    noted_edges = []

    # Iterate through each set of mutations in the cycles
    for mutations in mutations_in_snv_cycles:
        noted_snvs_in_cycle = []
        # Comment Max: Code-block für Site-Patterns
        # # Identify site patterns !=2.
        # # Note that in other cycles these SNVs can be realized in a 2-pattern.
        # # Therefore, the corresponding SNVs are ignored in the present cycle only.
        # for current_entry in mutations:
        #     number_of_sites = len([entry[0][0] for entry in mutations if entry[0][0] == current_entry[0][0]])
        #     if number_of_sites != 2:
        #         noted_snvs_in_cycle.append(current_entry[0])

        # Iterate through each mutation, edge, and day in the mutations

        # Comment Max: Code-Block für Korrektur von Bedinungen für noted_edge und noted_snvs_in_cycle
        # for mutation, edge, day in mutations: #edge length is time difference to start date + 1
        #     if edge not in noted_edges:
        #         noted_edges.append(edge)
        #         if mutation not in noted_snvs_in_cycle:
        #             tri_count.setdefault(mutation, {})
        #             tri_count[mutation].setdefault(day, 0)
        #             tri_count[mutation][day] += 1
        #             noted_snvs_in_cycle.append(mutation)
        # tri_table = []
        # for mutation, count in tri_count.items():
        #     tri_table.append(expand_timeseries(mutation, count, timerange))

        # if len(tri_count) == 0:
        #     print('No tRI signal found.')

        # tri_table_df = pd.DataFrame(tri_table, columns=['POS', 'REF', 'ALT'] + [str(i) for i in range(1, timerange + 1)])

        # return tri_table_df

        # Eventuell muss man dann expand_timeseries wie in evotrec2 nehmen:
        # def expand_timeseries(mutation, count, timerange):
        #     tri_timeseries = [mutation[0], mutation[1], mutation[2]]
        #     tri = 0
        #     for day in range(timerange):
        #         try: tri = tri + count[day+1]
        #         except: pass
        #         tri_timeseries.append(tri)

        #     return tri_timeseries

        for mutation, edge, day in mutations:
            # TODO: (Andreas) copy in correct version of the code
            if (mutation not in noted_snvs_in_cycle) and (edge not in noted_edges):
                tri_count.setdefault(mutation, {})
                tri_count[mutation].setdefault(day, 0)
                tri_count[mutation][
                    day
                ] += 1  # Increment the count for the mutation on the given day
                noted_snvs_in_cycle.append(mutation)
                noted_edges.append(edge)

    # Write the results to a CSV file
    with open(f"{output_filename}.csv", "w") as out_csv:
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
        print(f"Results written to {output_filename}.csv.")


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
    print(
        textwrap.dedent(
            """
        ===============================================================
        EVOtRec -- Topological recurrence analysis of genome alignments
        (c) 2024 Andreas Ott
        ===============================================================
        """
        )
    )

    # Retrieve metadata from input files
    filename, refseq, dates, start_date, end_date, timerange = retrieve_metadata(
        input_afasta=args.afasta,
        refseq_id=args.refseq_id,
        timeseries_flag=args.timeseries,
    )

    # Perform Hamming distance calculation
    print("Hammingdist...")
    distances = hammingdist.from_fasta(args.afasta, max_distance=4)
    distances.dump_sparse(filename + ".dist", threshold=3)

    # Perform MuRiT and Ripser analysis based on timeseries flag
    if args.timeseries is not False:
        print("MuRiT...")
        # Type assertion: start_date and dates are guaranteed to be non-None when timeseries_flag is True
        assert (
            start_date is not None
        ), "start_date should not be None when timeseries_flag is True"
        assert (
            dates is not None
        ), "dates should not be None when timeseries_flag is True"
        murit(
            input_distance_file=filename + ".dist",
            output_distance_file=filename + ".timedist",
            start_date=start_date,
            timerange=timerange,
            dates=dates,
        )
        print("Ripser...")
        os.system(
            f"ripser --dim 1 --format sparse --threshold {timerange+1} "
            f"{filename + '.timedist'} > {filename + '.ripser'}"
        )
    else:
        print("Ripser...")
        os.system(
            f"ripser --dim 1 --format sparse --threshold 2 "
            f"{filename + '.dist'} > {filename + '.ripser'}"
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
        mutations_in_snv_cycles=mutations_in_snv_cycles,
        timerange=timerange,
        output_filename=filename,
        timeseries_flag=args.timeseries,
    )

# Comment Max: code aus evotrec2:
# tri_table_df = tri_analysis(
#         filename=filename,
#         timeseries_flag=args.timeseries,
#         timerange=timerange,
#         mutations_in_snv_cycles=mutations_in_snv_cycles
#     )

#     def get_reference_site(row, **kwargs):
#         return kwargs['ref_positions'][row.loc['POS']-1]

#     def get_nt_mutation(row):
#         return row.loc['REF']+str(row.loc['REF_POS'])+row.loc['ALT']

# #    part_of_genome = ["test", 1, 10]
#     part_of_genome = ["spike", 21563, 25384]
# #    part_of_genome = ["whole", 266, 29674]

#     ref_positions = []
#     i = part_of_genome[1]
#     for x in refseq:
#         if x != '-':
#             ref_pos = i
#             i += 1
#         else:
#             ref_pos = -1
#         ref_positions.append(ref_pos)

#     tri_table_df['REF_POS'] = tri_table_df.apply(get_reference_site, axis=1, ref_positions=ref_positions)
#     tri_table_df['nt_mutation'] = tri_table_df.apply(get_nt_mutation, axis=1)
#     tri_table_df = tri_table_df.drop(columns=['POS', 'REF_POS', 'REF', 'ALT'])
#     tri_table_df = tri_table_df.set_index('nt_mutation')

#     print(tri_table_df)

#     tri_table_df.to_csv(f"{filename}_tri.csv")
