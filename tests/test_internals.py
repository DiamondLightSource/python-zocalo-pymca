import shutil
from pathlib import Path

import numpy as np
import pymca_zocalo.internals as internals
from pytest import fixture, raises


@fixture(scope="session")
def mock_visit_directory(tmp_path_factory):
    base_dir = tmp_path_factory.mktemp("base_dir")
    # Note, this works because base_dir has four levels, putting the visit at the 5th level (matching /dls structure)
    custom_dir = base_dir / "nt37183-3"
    custom_dir.mkdir(parents=True)
    return custom_dir


# @fixture
# def mock_cfg_file(monkeypatch, tmp_path):

#     temp_file = tmp_path / "pymca_new.cfg"
#     shutil.copyfile(Path("tests/pymca_new.cfg"), temp_file)

#     def mock_get_file_path():
#         return str(temp_file)
#     monkeypatch.setattr('/dls_sw/i03/software/pymca/pymca_new.cfg', mock_get_file_path)
#     return temp_file


def test_GIVEN_examplar_data_in_h5file_and_cfg_file_WHEN_run_auto_pymca_called_THEN_produces_expected_files(
    mock_visit_directory,
):
    test_data_file = Path("tests/sample_data.h5")
    temp_input_file = mock_visit_directory / test_data_file.name
    test_config_file = Path("tests/pymca_new.cfg")
    temp_config_file = mock_visit_directory / test_config_file.name
    shutil.copyfile(test_data_file, temp_input_file)
    shutil.copyfile(test_config_file, temp_config_file)
    internals.run_auto_pymca(
        str(temp_input_file),
        "0.0",
        "0.8",
        "(0, 0, 0)",
        "1.0",
        12700.0,
    )
    tmp_file_basename = Path(temp_input_file).stem
    output_mca = mock_visit_directory / f"{tmp_file_basename}.mca"
    output_html = mock_visit_directory / f"{tmp_file_basename}.html"
    output_png = mock_visit_directory / f"{tmp_file_basename}.png"
    pymca_dir = mock_visit_directory / "processed" / "pymca"
    output_fit = pymca_dir / "out" / f"{test_data_file.stem}.h5"
    output_results = pymca_dir / f"{test_data_file.stem}.results.dat"

    for file in [output_mca, output_html, output_png, output_fit, output_results]:
        assert file.is_file(), f"File '{file.name}' not written"


# def test_GIVEN_examplar_data_in_h5file_WHEN_run_auto_pymca_called_THEN_produces_expected_files(mock_visit_directory, mock_cfg_file):
#     test_data_file = Path("tests/sample_data.h5")
#     temp_input_file = mock_visit_directory / test_data_file.name
#     test_config_file =Path("tests/pymca_new.cfg")
#     temp_config_file = mock_visit_directory / test_config_file.name
#     shutil.copyfile(test_data_file, temp_input_file)
#     shutil.copyfile(test_config_file, temp_config_file)
#     internals.run_auto_pymca(
#         str(temp_input_file),
#         "0.0",
#         "0.8",
#         "(0, 0, 0)",
#         "1.0",
#         18000.0,
#     )
#     tmp_file_basename = Path(temp_input_file).stem
#     output_mca = mock_visit_directory / f"{tmp_file_basename}.mca"
#     output_html = mock_visit_directory / f"{tmp_file_basename}.html"
#     output_png = mock_visit_directory / f"{tmp_file_basename}.png"
#     pymca_dir = mock_visit_directory / "processed" / "pymca"
#     output_fit = pymca_dir / "out" / f"{test_data_file.stem}.h5"
#     output_results = pymca_dir / f"{test_data_file.stem}.results.dat"

#     for file in [output_mca, output_html, output_png, output_fit, output_results]:
#         assert file.is_file(), f"File '{file.name}' not written"


def test_GIVEN_fitted_peaks_WHEN_call_parse_raw_fluoro_THEN_get_expected_output():
    channel_energy = np.arange(0, 6, dtype=float)
    channel_counts = np.array([250.0, 1.1, 250.0, 250.0, 250.0])
    offset = 1.0
    beam_energy = 5.0
    peaks = [
        ("Se-K", 500.0, 40.0),
        ("Fe-K", 250.0, 20.0),
        ("Pb-L", 50.0, 10.0),
        ("Au-L", 40.0, 15.0),
        ("Pt-L", 30.0, 14.0),
        ("Hg-L", 10.0, 8.0),
    ]

    output_txt = internals.parse_raw_fluoro(
        channel_energy, channel_counts, beam_energy, peaks, offset
    )

    assert (
        "Element\tCounts\t%age\tExpected Emission Energies" in output_txt
    ), "Missing column labels"
    assert "Pt-L" in output_txt, "Expected peak missing from output"
    assert "Hg-L" not in output_txt, "Unexpected peak in output"
    assert "Counts (total): 1001.1" in output_txt, "Incorrect total counts"
    assert "(background): 9.1" in output_txt, "Incorrect background counts"


def test_GIVEN_low_counts_WHEN_call_parse_raw_fluoro_THEN_output_warning():
    channel_energy = np.arange(0, 6, dtype=float)
    channel_counts = np.array([4.0, 1.1, 4.0, 4.0, 4.0])
    offset = 1.0
    beam_energy = 5.0
    peaks = [
        ("Se-K", 500.0, 40.0),
        ("Fe-K", 250.0, 20.0),
        ("Pb-L", 50.0, 10.0),
        ("Au-L", 40.0, 15.0),
        ("Pt-L", 30.0, 14.0),
        ("Hg-L", 10.0, 8.0),
    ]

    output_txt = internals.parse_raw_fluoro(
        channel_energy, channel_counts, beam_energy, peaks, offset
    )
    assert (
        "No fluorescence peaks detected" in output_txt
    ), "No peaks warning message not written"
    assert "Se-K" not in output_txt, "Peak written to output when not supposed to"
    assert "Counts (total): 17.1" in output_txt, "Incorrect total counts"
    assert "(background): 9.1" in output_txt, "Incorrect background counts"


def test_GIVEN_no_peaks_WHEN_call_parse_raw_fluoro_THEN_output_warning():
    channel_energy = np.arange(0, 6, dtype=float)
    channel_counts = np.array([250.0, 1.1, 250.0, 250.0, 250.0])
    offset = 1.0
    beam_energy = 5.0
    peaks = []

    output_txt = internals.parse_raw_fluoro(
        channel_energy, channel_counts, beam_energy, peaks, offset
    )
    assert (
        "Counts found but no peaks of select elements fitted" in output_txt
    ), "Warning message not written"
    assert "Counts (total): 1001.1" in output_txt, "Incorrect total counts"
    assert "(background): 9.1" in output_txt, "Incorrect background counts"


def test_GIVEN_unsuitable_cutoff_WHEN_call_parse_raw_fluoro_THEN_raise_ValueError():
    channel_energy = np.arange(0, 6, dtype=float)
    channel_counts = np.array([250.0, 1.1, 250.0, 250.0, 250.0])
    offset = 6.0
    beam_energy = 5.0
    peaks = [
        ("Se-K", 500.0, 40.0),
        ("Fe-K", 250.0, 20.0),
        ("Pb-L", 50.0, 10.0),
        ("Au-L", 40.0, 15.0),
        ("Pt-L", 30.0, 14.0),
        ("Hg-L", 10.0, 8.0),
    ]

    with raises(ValueError):
        internals.parse_raw_fluoro(
            channel_energy, channel_counts, beam_energy, peaks, offset
        )


# @patch("pymca_zocalo.internals.open")
# def test_GIVEN_input_file_not_found_WHEN_find_cut_off_energy_called_THEN_raises_exception(
#     mock_open,
# ):
#     mock_open.side_effect = FileNotFoundError()
#     with raises(FileNotFoundError):
#         find_cut_off_energy("test_file_name", 0)


# @patch("pymca_zocalo.internals.open", new_callable=mock_open, read_data="1 0")
# def test_GIVEN_energy_of_1_and_cutoff_of_0_WHEN_find_cut_off_energy_called_THEN_returns_1(
#     mock_file,
# ):
#     line_of_greater_energy = find_cut_off_energy("test_file_name", 0)
#     assert line_of_greater_energy == 1


# @patch("pymca_zocalo.internals.open")
# def test_different_energies(mock_file):
#     test_data = "0 1\n5 10\n"
#     mock_open(mock_file, read_data=test_data)
#     line_of_greater_energy = find_cut_off_energy("test_file_name", 2)
#     assert line_of_greater_energy == 1


# @patch("pymca_zocalo.internals.open")
# def test_file_has_non_floats_different_energies(mock_file):
#     test_data = "hyadasd 1\n5 10\n"
#     mock_open(mock_file, read_data=test_data)
#     line_of_greater_energy = find_cut_off_energy("test_file_name", 2)
#     assert line_of_greater_energy == 1


# @patch("pymca_zocalo.internals.open")
# def test_file_nothing_in_then_gives_none(mock_file):
#     test_data = ""
#     mock_open(mock_file, read_data=test_data)
#     line_of_greater_energy = find_cut_off_energy("test_file_name", 2)
#     assert line_of_greater_energy is None
