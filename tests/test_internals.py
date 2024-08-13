from unittest.mock import mock_open, patch

from pymca_zocalo.internals import find_cut_off_energy
from pytest import raises


@patch("pymca_zocalo.internals.open")
def test_GIVEN_input_file_not_found_WHEN_find_cut_off_energy_called_THEN_raises_exception(
    mock_open,
):
    mock_open.side_effect = FileNotFoundError()
    with raises(FileNotFoundError):
        find_cut_off_energy("test_file_name", 0)


@patch("pymca_zocalo.internals.open", new_callable=mock_open, read_data="1 0")
def test_GIVEN_energy_of_1_and_cutoff_of_0_WHEN_find_cut_off_energy_called_THEN_returns_1(
    mock_file,
):
    line_of_greater_energy = find_cut_off_energy("test_file_name", 0)
    assert line_of_greater_energy == 1


@patch("pymca_zocalo.internals.open")
def test_different_energies(mock_file):
    test_data = "0 1\n5 10\n"
    mock_open(mock_file, read_data=test_data)
    line_of_greater_energy = find_cut_off_energy("test_file_name", 2)
    assert line_of_greater_energy == 1


@patch("pymca_zocalo.internals.open")
def test_file_has_non_floats_different_energies(mock_file):
    test_data = "hyadasd 1\n5 10\n"
    mock_open(mock_file, read_data=test_data)
    line_of_greater_energy = find_cut_off_energy("test_file_name", 2)
    assert line_of_greater_energy == 1


@patch("pymca_zocalo.internals.open")
def test_file_nothing_in_then_gives_none(mock_file):
    test_data = ""
    mock_open(mock_file, read_data=test_data)
    line_of_greater_energy = find_cut_off_energy("test_file_name", 2)
    assert line_of_greater_energy is None
