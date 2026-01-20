import os
import sys
import pytest
import subprocess
import filecmp

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from evotrec import (  # noqa: E402
    retrieve_metadata,
    murit,
    retrieve_snv_cycles,
    retrieve_sequences_in_cycles,
    retrieve_mutations_in_cycles,
    expand_timeseries,
    tri_analysis,
)


@pytest.fixture
def setup_data():
    return {
        # for retrieve_metadata
        "input_fasta": "test.fasta",
        "refseq_id": "seq0|2020-01-01",
        "timeseries_flag": True,
        # for murit
        "in_dist": "test.dist",
        "out_timedist": "test.timedist",
        "start_date": "2020-01-01",
        "timerange": 4,
        "dates": ["2020-01-01", "2020-01-02", "2020-01-03", "2020-01-04"],
        # for retrieve_snv_cycles
        "in_ripser": "test.ripser",
        # for retrieve_sequences_in_cycles
        "sequences_in_snv_cycles": {0: "AA", 1: "AC", 2: "CA", 3: "CC"},
        "snv_cycles": [[((0, 1), 1), ((0, 2), 1), ((1, 3), 1), ((2, 3), 1)]],
        "snv_indices": [0, 1, 2, 3],
        # for retrieve_mutations_in_cycles
        "refseq": "AA",
        # for expand_timeseries
        "mutation": [0, "A", "C"],
        "count": {0: 0, 1: 1, 2: 2, 3: 1},
        # for tri_analysis
        "output_filename": "test",
        "mutations_in_snv_cycles": [[(((2, "A", "C"), (0, 1), 1))]],
    }


def test_retrieve_metadata(mocker, setup_data):
    # Mock the file reading operation for the input fasta file
    mocked_open = mocker.mock_open(
        read_data=">seq0|2020-01-01\nAA\n>seq1|2020-01-02\nAC\n>seq2|2020-01-03\nCA\n>seq3|2020-01-04\nCC\n"
    )
    mocker.patch("builtins.open", mocked_open)

    # Get the fixture data
    data = setup_data

    # Call the function under test
    result = retrieve_metadata(
        input_afasta=data["input_fasta"],
        refseq_id=data["refseq_id"],
        timeseries_flag=data["timeseries_flag"],
    )

    # Verify the output
    assert result == (
        "test",
        "AA",
        ["2020-01-01", "2020-01-02", "2020-01-03", "2020-01-04"],
        "2020-01-01",
        "2020-01-04",
        4,
    )


def test_murit(mocker, setup_data):
    # Mock the file reading operation for the input distance file
    mocked_open = mocker.mock_open(
        read_data="1 0 1\n2 0 1\n2 1 2\n3 0 2\n3 1 1\n3 2 1\n"
    )
    mocker.patch("builtins.open", mocked_open)

    # Get the fixture data
    data = setup_data

    # Call the function under test
    murit(
        input_distance_file=data["in_dist"],
        output_distance_file=data["out_timedist"],
        start_date=data["start_date"],  # Corrected parameter name
        timerange=data["timerange"],
        dates=data["dates"],
    )

    # Verify the output file contains the expected content
    handle = mocked_open()
    handle.write.assert_any_call("1 0 2\n")
    handle.write.assert_any_call("2 0 3\n")
    handle.write.assert_any_call("2 1 5\n")
    handle.write.assert_any_call("3 0 5\n")
    handle.write.assert_any_call("3 1 4\n")
    handle.write.assert_any_call("3 2 4\n")


def test_retrieve_snv_cycles(mocker, setup_data):
    # Mock the file reading operation for the input distance file
    mocked_open = mocker.mock_open(
        read_data="persistent homology intervals in dim 1:\
            \n [1,2):\n\
            {[0,1] (1), [0,2] (1), [1,3] (1), [2,3] (1)}\n\
            {[0,1] (1), [0,2] (1), [1,3] (1), [2,3] (1)}\n"
    )
    mocker.patch("builtins.open", mocked_open)

    # Get the fixture data
    data = setup_data

    # Call the function under test
    result = retrieve_snv_cycles(data["in_ripser"])

    # Verify the output
    assert result == (
        [[((0, 1), 1), ((0, 2), 1), ((1, 3), 1), ((2, 3), 1)]],
        [0, 1, 2, 3],
    )


def test_retrieve_sequences_in_cycles(mocker, setup_data):
    # Mock the file reading operation for the input distance file
    mocked_open = mocker.mock_open(
        read_data=">seq0|2020-01-01|\nAA\n>seq1|2020-01-02|\nAC\n>seq2|2020-01-03|\nCA\n>seq3|2020-01-04|\nCC\n"
    )
    mocker.patch("builtins.open", mocked_open)

    # Get the fixture data
    data = setup_data

    # Call the function under test
    result = retrieve_sequences_in_cycles(data["snv_indices"], data["input_fasta"])

    # Verify the output
    assert result == {0: "AA", 1: "AC", 2: "CA", 3: "CC"}


def test_retrieve_mutations_in_cycles(setup_data):
    # Get the fixture data
    data = setup_data

    # Call the function under test
    result = retrieve_mutations_in_cycles(
        data["snv_cycles"], data["sequences_in_snv_cycles"], data["refseq"]
    )

    # Verify the output
    assert result == [
        [
            ((2, "A", "C"), (0, 1), 1),
            ((1, "A", "C"), (0, 2), 1),
            ((1, "A", "C"), (1, 3), 1),
            ((2, "A", "C"), (2, 3), 1),
        ]
    ]


def test_expand_timeseries(setup_data):
    # Get the fixture data
    data = setup_data

    # Call the function under test
    result = expand_timeseries(data["mutation"], data["count"], data["timerange"])

    # Verify the output
    assert result == "0,A,C,1,3,4,4"


def test_tri_analysis(mocker, setup_data):
    # Mock the file reading operation for the input distance file
    mock_file = mocker.patch("builtins.open", mocker.mock_open())

    # Get the fixture data
    data = setup_data

    # Call the function under test
    tri_analysis(
        mutations_in_snv_cycles=data["mutations_in_snv_cycles"],
        output_filename=data["output_filename"],
        timerange=data["timerange"],
        timeseries_flag=data["timeseries_flag"],
    )

    # Verify the output file contains the expected content
    mock_file().write.assert_called_with("2,A,C,1,1,1,1\n")


def test_evotrec_pipeline():
    # Define the command to run the evotrec.py script
    input_file = "tests/test.fasta"
    ref_seq = "NC_045512.2|Severeacuterespiratorysyndromecoronavirus2isolateWuhan-Hu-1,completegenome|China|2019-12-30"

    command = ["python", "evotrec.py", input_file, ref_seq]

    # Run the command
    result = subprocess.run(command, capture_output=True, text=True)

    # Check if the script ran successfully
    assert result.returncode == 0, f"Script failed with error: {result.stderr}"

    # Compare the produced ripser output file with the expected ripser output file
    produced_ripser_output = input_file.replace(".fasta", ".ripser")
    expected_ripser_output = "tests/expected.ripser"
    assert filecmp.cmp(
        produced_ripser_output, expected_ripser_output
    ), "The produced ripser output does not match the expected output."

    # Compare the produced tri output file with the expected tri output file
    produced_tri_output = input_file.replace(".fasta", ".csv")
    expected_tri_output = "tests/expected.csv"
    assert filecmp.cmp(
        produced_tri_output, expected_tri_output
    ), "The produced tri output does not match the expected output."


def test_evotrec_pipeline_timeseries():
    # Define the command to run the evotrec.py script with --timeseries flag
    input_file = "tests/test.fasta"
    ref_seq = "NC_045512.2|Severeacuterespiratorysyndromecoronavirus2isolateWuhan-Hu-1,completegenome|China|2019-12-30"

    command = ["python", "evotrec.py", input_file, ref_seq, "--timeseries"]

    # Run the command
    result = subprocess.run(command, capture_output=True, text=True)

    # Check if the script ran successfully
    assert result.returncode == 0, f"Script failed with error: {result.stderr}"

    # Compare the produced ripser output file with the expected ripser output file
    produced_ripser_output = input_file.replace(".fasta", ".ripser")
    expected_ripser_output = "tests/expected_timeseries.ripser"
    assert filecmp.cmp(
        produced_ripser_output, expected_ripser_output
    ), "The produced ripser output does not match the expected output."

    # Compare the produced tri output file with the expected tri output file
    produced_tri_output = input_file.replace(".fasta", ".csv")
    expected_tri_output = "tests/expected_timeseries.csv"
    assert filecmp.cmp(
        produced_tri_output, expected_tri_output
    ), "The produced tri output does not match the expected output."
