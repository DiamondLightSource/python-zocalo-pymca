import os
import re
import shutil
from datetime import datetime
from importlib import resources
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import xraylib as xrl
from PyMca5.PyMca import McaAdvancedFitBatch
from PyMca5.PyMcaIO import ConfigDict


def parse_raw_fluoro(channel_energy, channel_counts, peaks, cutoff_channel):
    """Calculate total and background counts from raw spectrum, then
    return results as formatted text for top 5 fitted peaks if there
    is sufficient signal above the background"""
    output_txt = []
    edge_mapper = {
        "K": xrl.K_SHELL,
        "L": xrl.L3_SHELL,
    }

    line_mapper = {
        "K": xrl.KL3_LINE,
        "L": xrl.L3M5_LINE,
    }
    # Get total counts and background counts up to a cutoff energy
    if not cutoff_channel:
        cutoff_channel = len(channel_energy)
    total_count = np.sum(channel_counts[0:cutoff_channel])
    # Maximum of 2 counts per channel contribute to background
    background_count = np.sum(np.minimum(2, channel_counts[0:cutoff_channel]))

    if (total_count - background_count) > 100:
        if len(peaks):
            output_txt.append("Element\tCounts\t%age\tExpected Emission Energies")
            for peak_num, (emis_line, counts, _sigma) in enumerate(peaks, start=1):
                symbol, edge = emis_line.split("-")
                z = xrl.SymbolToAtomicNumber(symbol)
                edgesEnergies = "<b>{:g}</b>,{:g}".format(
                    xrl.LineEnergy(z, line_mapper[edge]) * 1000.0,
                    xrl.EdgeEnergy(z, edge_mapper[edge]) * 1000.0,
                )
                if counts >= 100000:
                    pk_counts = int(counts)
                else:
                    pk_counts = round(counts, 1)
                output_txt.append(
                    "{}\t{:g}\t{:g}\t{}".format(
                        emis_line,
                        pk_counts,
                        round(100 * counts / total_count, 1),
                        edgesEnergies,
                    )
                )
                if peak_num == 5:
                    break
        else:
            output_txt.append(
                "Counts found but no peaks of select elements fitted, check the spectrum manually"
            )
    else:
        output_txt.append("No fluorescence peaks detected, try a higher transmission")

    output_txt.append(
        "\nCounts (total): {:g} (background): {:g}".format(
            total_count, background_count
        )
    )
    return "\n".join(output_txt)


# Parses a file and sorts spectral peaks with largest area first
# and prints them to stdout
def parse_spec_fit(name):
    try:
        f = h5py.File(name, "r")
        parameters = f[list(f.keys())[0] + "/xrf_fit/results/parameters"]
        all_fit = filter(
            lambda name: not name.startswith("Scatter")
            and not name.endswith("_errors"),
            parameters,
        )

        def mapper(name):
            return [
                name.replace("_", "-"),
                parameters[name][0, 0],
                parameters[name + "_errors"][0, 0],
            ]

        mapped = map(mapper, all_fit)
        peaks = sorted(mapped, key=lambda p: p[1], reverse=True)

    except KeyError:
        raise KeyError(
            f"PyMCA fit parameters not found at '/xrf_fit/results/parameters' in '{name}'"
        )

    return peaks


def configure_cfg(cfg_file, peaks, energy_kev):
    with open(cfg_file, "r") as file:
        lines = file.readlines()
        peaks_section = False
        peaks_section_end = False
        peaks_found = False
        peak_line_num = len(lines)

        for line_num, line in enumerate(lines):
            # Replace energy in file with actual photon energy
            if match := re.match(r"energy = (\d+(\.\d+)?),", line):
                lines[line_num] = line.replace(match.group(1), str(energy_kev))
            # Look for peaks
            if line.strip() == "[peaks]":
                peaks_section = True
            elif peaks_section and not peaks_section_end:
                if "[" in line:
                    peaks_section_end = True
                    peak_line_num = line_num
                if "=" in line:
                    peaks_found = True
        # Add peaks to file if not found
        if not peaks_found:
            if not peaks_section:
                lines.append("\n[peaks]\n")
                peak_line_num += 1
            lines.insert(peak_line_num, peaks)

    with open(cfg_file, "w") as file:
        file.writelines(lines)


def parse_elements(cutoff_energy, first_channel_energy, beamline):
    """
    Filters a dictionary of element/edge combinations likely to be encountered in an
     XRF experiment and returns those that are within the fittable region of the spectrum.
    Args:
        cutoff_energy (float): The energy cutoff defining the scattering region of the
         spectrum. Used as an upper limit for element filtering. Expressed in electron
         volts (eV).
        first_channel_energy (float): The energy level of the first channel. Used as a
          lower limit for element filtering. Expressed in electron volts (eV).
        beamline (str): The beamline that the MCA spectrum was recorded on.

    Returns:
        dict[str]: A dictionary of the filtered elements and edges. The element and
        edge form a key-value pair.
    """
    elements = {
        "Ti": "K",
        "V": "K",
        "Cr": "K",
        "Mn": "K",
        "Fe": "K",
        "Co": "K",
        "Ni": "K",
        "Cu": "K",
        "Zn": "K",
        "As": "K",
        "Se": "K",
        "Br": "K",
        "Sr": "K",
        "Mo": "K",
        "I": "L",
        "Xe": "L",
        "Gd": "L",
        "W": "L",
        "Os": "L",
        "Ir": "L",
        "Pt": "L",
        "Au": "L",
        "Hg": "L",
        "Pb": "L",
    }
    # Additional elements for i23
    i23_elements = {
        "Na": "K",
        "K": "K",
        "Mg": "K",
        "Mo": "L",
        "P": "K",
        "S": "K",
        "Ag": "L",
        "Cl": "K",
        "Ca": "K",
        "Cd": "L",
    }

    edge_mapper = {
        "K": xrl.K_SHELL,
        "L": xrl.L3_SHELL,
    }

    def _check_edge(symbol, edge):
        Z = xrl.SymbolToAtomicNumber(symbol)
        shell = edge_mapper[edge]
        return (
            float(cutoff_energy) > xrl.LineEnergy(Z, shell) * 1000.0
            and float(first_channel_energy) < xrl.LineEnergy(Z, shell) * 1000.0
        )

    if beamline == "i23":
        elements.update(i23_elements)

    valid_peaks = {}
    for symbol, edge in elements.items():
        if _check_edge(symbol, edge):
            valid_peaks[symbol] = edge

    return valid_peaks


def spectrum_to_mca(channel_counts, calib, output_file, input_file, selection):
    # Converts a spectrum as a 1D numpy array of counts to an MCA file
    timestamp = datetime.now().strftime("%a %b %d %H:%M:%S %Y")
    header = [
        f"#F {output_file}",
        f"#D {timestamp}",
        "",
        f"#S 1 {input_file} {selection['entry']}",
        f"#D {timestamp}",
        "#@MCA %16C",
        f"#@CHANN {len(channel_counts)} 0 {len(channel_counts)-1} 1",
        f"#@CALIB {calib['zero']} {calib['gain']} 0.0",
    ]
    header = "\n".join(header)

    data_out = ["@A"]
    for _i, counts in enumerate(channel_counts, start=1):
        # Create string with removed trailing 0s
        counts_str = str(counts).rstrip("0").rstrip(".")
        data_out.append(counts_str)
        if not _i % 16 and _i != len(channel_counts):
            data_out.append("\\\n")

    data_out = " ".join(data_out)

    with open(output_file, "w") as f:
        f.write(header + "\n" + data_out)


def plot_fluorescence_spectrum(
    out_file, input_file, channel_energy, channel_counts, beam_energy, cutoff_energy
):
    i_cutoff = np.argmax(channel_energy > cutoff_energy)

    # Create plot of spectrum, split by the cutoff
    plt.figure(figsize=(8, 6))
    plt.title(f"Fluorescence Spectrum {input_file}", fontsize=12)
    plt.xlabel("Energy (eV)", fontsize=12)
    plt.ylabel("Number of Counts", fontsize=12)
    plt.xlim(channel_energy[0], beam_energy)
    plt.minorticks_on()
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.plot(
        channel_energy[:i_cutoff],
        channel_counts[:i_cutoff],
        label="Data 1",
        linestyle="-",
        linewidth=1,
        color="darkorange",
    )
    plt.plot(
        channel_energy[i_cutoff:],
        channel_counts[i_cutoff:],
        label="Data 2",
        linestyle="-",
        linewidth=0.5,
        color="deepskyblue",
    )
    if i_cutoff:
        plt.axvline(x=cutoff_energy, color="black", linestyle="--", linewidth=1)
    plt.savefig(out_file, format="png")


def read_h5file(src_data_file, h5path):
    h5path_parts = Path(h5path).parts
    if h5path_parts[0] == os.sep:
        root_group = h5path_parts[1]
        y_data_path = "/" + "/".join(h5path_parts[2:])
    else:
        root_group = h5path_parts[0]
        y_data_path = "/" + "/".join(h5path_parts[1:])
    selection = {"entry": root_group, "y": y_data_path}

    # Read in data from h5 file
    calibration_path = Path(h5path).parent / "calibration"
    with h5py.File(src_data_file, "r") as hf:
        cal_list = hf[str(calibration_path)]
        calibration = {"zero": cal_list[0], "gain": cal_list[1]}
        # Squeeze to collapse 1-length dimensions
        channel_counts = np.squeeze(hf[str(h5path)])
        # Calculate channel energies (in eV) from calibration
    channel_energy = np.array(
        [
            (calibration["zero"] + calibration["gain"] * _i) * 1000
            for _i in range(len(channel_counts))
        ]
    )

    return selection, calibration, channel_counts, channel_energy


def run_auto_pymca(
    src_data_file,
    omega,
    transmission,
    samplexyz,
    acq_time,
    beam_energy,
    src_cfg_file=None,
    h5path=None,
):
    # Set paths
    src_data_dir = Path(src_data_file).parent
    filename_stem = Path(src_data_file).stem
    filename = Path(src_data_file).name
    filepath_parts = Path(src_data_file).parts
    visit_dir = Path(*filepath_parts[:6])
    rel_dir_path = Path(*filepath_parts[6:-1])
    mca_path = src_data_dir / f"{filename_stem}.mca"
    beamline = filepath_parts[2]

    config_changes = {}
    energy_kev = float(beam_energy) / 1000.0

    if h5py.is_hdf5(src_data_file):
        src_cfg_file = resources.files("pymca_zocalo.data") / "pymca_new.cfg"
        if not h5path:
            h5path = "entry/instrument/detector/data"
        selection, calibration, channel_counts, channel_energy = read_h5file(
            src_data_file, h5path
        )
        # Create a .mca file (not used by this code but a more user friendly file format for PyMCA GUI)
        spectrum_to_mca(channel_counts, calibration, mca_path, src_data_file, selection)

        config_changes["detector"] = {
            "zero": calibration["zero"],
            "gain": calibration["gain"],
        }

    elif src_data_file.endswith((".dat", ".mca")):
        src_data_file = mca_path
        dat_path = src_data_dir / f"{filename_stem}.dat"
        selection = {}
        spectrum_data = np.loadtxt(dat_path, delimiter="\t", skiprows=3, usecols=(0, 1))
        channel_energy = spectrum_data[:, 0]
        channel_counts = spectrum_data[:, 1]

        if src_cfg_file is None:
            src_cfg_file = Path("/dls_sw") / beamline / "software/pymca/pymca_new.cfg"
        else:
            src_cfg_file = Path(src_cfg_file)

    else:
        raise ValueError(
            "Invalid input file - Input file must be in .dat, .mca, or hdf5 format"
        )

    if not src_cfg_file.exists():
        raise FileNotFoundError(f"Config file '{src_cfg_file}' does not exist")

    output_dir = visit_dir / "processed" / "pymca" / rel_dir_path
    data_dir = output_dir / "data"
    results_dir = output_dir / "out"
    data_file = data_dir / filename
    cfg_file = output_dir / f"{filename_stem}.cfg"

    data_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    shutil.copyfile(src_data_file, data_file)
    # Set a cutoff region 500eV below Compton edge to exclude scattering peaks from fit
    compton_edge = beam_energy / (1 + 2 * beam_energy / 5.110e5)
    cutoff_energy = compton_edge - 500

    config_dict = ConfigDict.ConfigDict()
    config_dict.read(src_cfg_file)

    first_channel_energy = config_dict["detector"]["zero"] * 1000.0

    peaks = parse_elements(cutoff_energy, first_channel_energy, beamline)

    if not data_file.exists():
        raise FileNotFoundError(f"File {data_file} does not exist")

    # Edit and write new config to file
    for main_dict, sub_dict in config_changes.items():
        config_dict[main_dict].update(sub_dict)
    if not config_dict.get("peaks"):
        config_dict["peaks"] = peaks
    config_dict["fit"]["energy"][0] = energy_kev
    fit_xmax = config_dict["fit"]["xmax"]
    cutoff_energy_kev = cutoff_energy / 1000.0
    cutoff_channel = int(
        (cutoff_energy_kev - config_dict["detector"]["zero"])
        / config_dict["detector"]["gain"]
    )
    if cutoff_channel < fit_xmax:
        config_dict["fit"]["xmax"] = cutoff_channel
    config_dict["fit"]["scatterflag"] = 0
    config_dict.write(cfg_file)

    if not cfg_file.exists():
        raise FileNotFoundError(f"Failed to write config file to {cfg_file}")

    # Create a batch fit object
    pymca_batch_fit_obj = McaAdvancedFitBatch.McaAdvancedFitBatch(
        str(cfg_file), [str(data_file)], str(results_dir), 0, 250.0, selection=selection
    )

    # ProcessList method fits data and writes results to .h5 file
    pymca_batch_fit_obj.processList()
    fit_data_file = results_dir / f"{filename_stem}.h5"

    if not fit_data_file.exists():
        raise FileNotFoundError(f"Results file {fit_data_file} could not be opened")

    fitted_peaks = parse_spec_fit(fit_data_file)

    results_file = output_dir / f"{filename_stem}.results.dat"
    results_txt = [f"{peak[0]} {peak[1]} {peak[2]}" for peak in fitted_peaks]
    results_txt = "\n".join(results_txt)
    with open(results_file, "w") as f:
        f.write(results_txt)

    pymca_output = parse_raw_fluoro(
        channel_energy, channel_counts, fitted_peaks, cutoff_channel
    )

    plot_output_file = src_data_dir / f"{filename_stem}.png"

    plot_fluorescence_spectrum(
        plot_output_file,
        src_data_file,
        channel_energy,
        channel_counts,
        beam_energy,
        cutoff_energy,
    )

    snapshot_dir = src_data_dir / "jpegs" / rel_dir_path

    # Find crystal snapshot if it exists
    try:
        snapshot_files = sorted(
            snapshot_dir.glob(f"{filename_stem}*[0-9].png"),
            key=lambda p: p.stat().st_mtime,
        )
        snapshot_filename = snapshot_files[0].name
    except Exception:
        snapshot_filename = ""

    html_file = src_data_dir / f"{filename_stem}.html"

    timestamp = datetime.now().strftime("%a %b %-d, %Y - %T")

    html_file_contents = f"""<html><head><title>{filename}</title></head><body>
<style type="text/css">
table td {{ padding: 5px;}}
table.table2 {{ border-collapse: collapse;}}
table.table2 td {{ border-style: none; font-size: 80%; padding: 2px;}}
</style>
<table border width=640><tr><td align=center>{timestamp}<br/><a href="{rel_dir_path}">{filename}</a></td></tr>
<tr><td><table class="table2" width=100%>
<tr><th colspan=4 align=center>Fluorescence Spectrum</td></tr>
<tr><td>Beamline Energy:</td><td>{beam_energy}eV</td><td>Omega:</td><td>{omega}&deg;</td></tr>
<tr><td>Acq Time:</td><td>{acq_time}s</td><td>Trans:</td><td>{transmission}%</td></tr>
<tr><td>Sample Position:</td><td colspan=3>{samplexyz}</td></tr>
</table></td></tr>
<tr><td>Automated PyMca results<br /><pre>{pymca_output}</pre><a href="http://www.diamond.ac.uk/dms/MX/Common/Interpreting-AutoPyMCA/Interpreting%20AutoPyMCA.pdf">Guide to AutoPyMCA</a> (pdf)</td></tr>
<tr><td><img src="{plot_output_file}" /></td></tr>
<tr><td><img src="{snapshot_dir/snapshot_filename}" alt="Snapshot not taken" width=640 /></td></tr></table>
"""

    with open(html_file, "w") as f:
        f.write(html_file_contents)
