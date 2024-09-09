import shutil
from pathlib import Path

import numpy as np
import pymca_zocalo.internals as internals
from pytest import fixture, raises


@fixture
def mock_visit_directory(tmp_path):
    # Note, this works because tmp_path has four levels, putting the visit at the 5th level (matching /dls structure)
    custom_dir = tmp_path / "nt37183-3"
    custom_dir.mkdir(parents=True)
    yield custom_dir
    shutil.rmtree(custom_dir)


def test_GIVEN_examplar_data_in_h5file_and_cfg_file_WHEN_run_auto_pymca_called_THEN_produces_expected_files(
    mock_visit_directory,
):
    test_data_file = Path("tests/sample_data.h5")
    temp_input_file = mock_visit_directory / test_data_file.name
    shutil.copyfile(test_data_file, temp_input_file)
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


def test_GIVEN_examplar_data_in_mca_format_WHEN_run_auto_pymca_called_THEN_produces_expected_files(
    mock_visit_directory,
):
    test_mca_file = Path("tests/sample_data.mca")
    test_dat_file = Path("tests/sample_data.dat")
    test_cfg_file = Path("tests/pymca_new.cfg")
    temp_mca_file = mock_visit_directory / test_mca_file.name
    temp_dat_file = mock_visit_directory / test_dat_file.name
    temp_cfg_file = mock_visit_directory / test_cfg_file.name
    shutil.copyfile(test_mca_file, temp_mca_file)
    shutil.copyfile(test_dat_file, temp_dat_file)
    shutil.copyfile(test_cfg_file, temp_cfg_file)
    internals.run_auto_pymca(
        str(temp_mca_file),
        "0.0",
        "0.4",
        "(-1309, -227.0, 230.0)",
        "1.0",
        18000.0,
        src_cfg_file=temp_cfg_file,
    )
    tmp_file_basename = Path(temp_mca_file).stem
    output_mca = mock_visit_directory / f"{tmp_file_basename}.mca"
    output_html = mock_visit_directory / f"{tmp_file_basename}.html"
    output_png = mock_visit_directory / f"{tmp_file_basename}.png"
    pymca_dir = mock_visit_directory / "processed" / "pymca"
    output_fit = pymca_dir / "out" / f"{tmp_file_basename}.h5"
    output_results = pymca_dir / f"{tmp_file_basename}.results.dat"

    for file in [output_mca, output_html, output_png, output_fit, output_results]:
        assert file.is_file(), f"File '{file.name}' not written"


def test_GIVEN_examplar_data_in_dat_format_WHEN_run_auto_pymca_called_THEN_produces_expected_files(
    mock_visit_directory,
):
    test_mca_file = Path("tests/sample_data.mca")
    test_dat_file = Path("tests/sample_data.dat")
    test_cfg_file = Path("tests/pymca_new.cfg")
    temp_mca_file = mock_visit_directory / test_mca_file.name
    temp_dat_file = mock_visit_directory / test_dat_file.name
    temp_cfg_file = mock_visit_directory / test_cfg_file.name
    shutil.copyfile(test_mca_file, temp_mca_file)
    shutil.copyfile(test_dat_file, temp_dat_file)
    shutil.copyfile(test_cfg_file, temp_cfg_file)
    internals.run_auto_pymca(
        str(temp_dat_file),
        "0.0",
        "0.4",
        "(-1309, -227.0, 230.0)",
        "1.0",
        18000.0,
        src_cfg_file=temp_cfg_file,
    )
    tmp_file_basename = Path(temp_mca_file).stem
    output_mca = mock_visit_directory / f"{tmp_file_basename}.mca"
    output_html = mock_visit_directory / f"{tmp_file_basename}.html"
    output_png = mock_visit_directory / f"{tmp_file_basename}.png"
    pymca_dir = mock_visit_directory / "processed" / "pymca"
    output_fit = pymca_dir / "out" / f"{tmp_file_basename}.h5"
    output_results = pymca_dir / f"{tmp_file_basename}.results.dat"

    for file in [output_mca, output_html, output_png, output_fit, output_results]:
        assert file.is_file(), f"File '{file.name}' not written"


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
