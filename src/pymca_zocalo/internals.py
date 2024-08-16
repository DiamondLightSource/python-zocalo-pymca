import os
import re
import shutil
from datetime import datetime
from glob import glob
from pathlib import Path, PurePath

import h5py
import matplotlib.pyplot as plt
import numpy as np
import xraylib as xrl
from PyMca5.PyMca import McaAdvancedFitBatch


def parse_raw_fluoro(channel_energy, channel_counts, beam_energy, peaks, offset):
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
    cutoff_energy = beam_energy - offset
    cutoff_idx = np.argmax(channel_energy > cutoff_energy)
    if not cutoff_idx:
        raise ValueError(
            "No cutoff index found - all recorded channels are too close in energy to the beam"
        )
    total_count = np.sum(channel_counts[0:cutoff_idx])
    # Maximum of 2 counts per channel contribute to background
    background_count = np.sum(np.minimum(2, channel_counts[0:cutoff_idx]))

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


def parse_elements(energy):
    """
    Filters a dictionary of element/edge combinations likely to be encountered in an
     XRF experiment and returns those that are below the given photon energy.
    Args:
        energy (float): The energy level to filter the elements by, expressed in
         electron volts (eV).

    Returns:
        str: A string representation of the filtered elements. Each element is on
         a new line and is formatted as 'symbol = edge'.
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

    edge_mapper = {
        "K": xrl.K_SHELL,
        "L": xrl.L3_SHELL,
    }

    def _check_edge(element):
        symbol, edge = element
        Z = xrl.SymbolToAtomicNumber(symbol)
        shell = edge_mapper[edge]
        return float(energy) > xrl.EdgeEnergy(Z, shell) * 1000.0

    return "\n".join(
        (
            "{} = {}".format(*element)
            for element in filter(_check_edge, elements.items())
        )
    )


def find_cut_off_energy(inputFile, cutoffenergy):
    with open(inputFile, "r") as f:
        for i, line in enumerate(f):
            try:
                energy = float(line.split()[0])
                if energy > float(cutoffenergy):
                    return i
            except (ValueError, IndexError):
                pass


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
    outFile, inputFile, channel_energy, channel_counts, beam_energy, offset=1200
):
    cutoff_energy = beam_energy - offset
    i_cutoff = np.argmax(channel_energy > cutoff_energy)

    # Create plot of spectrum, split by the cutoff
    plt.figure(figsize=(8, 6))
    plt.title(f"Fluorescence Spectrum {inputFile}", fontsize=12)
    plt.xlabel("Energy (eV)", fontsize=12)
    plt.ylabel("Number of Counts", fontsize=12)
    plt.xlim(0, beam_energy)
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
    plt.axvline(x=cutoff_energy, color="black", linestyle="--", linewidth=1)
    plt.savefig(outFile, format="png")


def run_auto_pymca(
    inputFile,
    omega,
    transmission,
    samplexyz,
    acqTime,
    beam_energy,
    CFGFile=None,
    peaksFile=None,
    h5path=None,
):
    if not inputFile.endswith((".dat", ".mca", ".h5")):
        raise ValueError("inputFile must end with .dat, .mca or .h5")
    # For historical purposes, .dat given as input by GDA but .mca file is used
    if inputFile.endswith(".dat"):
        inputFile = os.path.splitext(inputFile)[0] + ".mca"

    selection = {}
    if h5py.is_hdf5(inputFile):
        if not h5path:
            h5path = "entry/data/data"
        h5path_parts = Path(h5path).parts
        if h5path_parts[0] == os.sep:
            root_group = h5path_parts[1]
            y_data_path = "/" + "/".join(h5path_parts[2:])
        else:
            root_group = h5path_parts[0]
            y_data_path = "/" + "/".join(h5path_parts[1:])
        selection = {"entry": root_group, "y": y_data_path}

    FilePrefix = os.path.splitext(os.path.basename(inputFile))[0]
    file_name = os.path.basename(inputFile)
    pure_path = PurePath(inputFile).parts
    VisitDir = pure_path[:6]
    RestOfFilename = os.path.join(*pure_path[6:])
    RestOfDirs = os.path.dirname(RestOfFilename)
    BEAMLINE = pure_path[2]

    peaks = None
    cutoff_offset = 1000

    if CFGFile is None:
        CFGFile = os.path.join("/dls_sw", BEAMLINE, "software/pymca/pymca_new.cfg")

    if not os.path.isfile(CFGFile):
        raise FileNotFoundError(f"Config file '{CFGFile}' does not exist")

    energy_keV = float(beam_energy) / 1000.0

    OutputDir = os.path.join(*VisitDir, "processed/pymca", RestOfDirs)
    DataDir = os.path.join(OutputDir, "data")
    ResultsDir = os.path.join(OutputDir, "out")
    cfg_path = os.path.join(OutputDir, FilePrefix + ".cfg")

    Path(OutputDir).mkdir(parents=True, exist_ok=True)
    Path(DataDir).mkdir(parents=True, exist_ok=True)
    Path(ResultsDir).mkdir(parents=True, exist_ok=True)

    os.chdir(OutputDir)
    shutil.copy2(inputFile, DataDir)
    shutil.copyfile(CFGFile, cfg_path)

    if peaksFile is not None:
        with open(peaksFile) as f:
            peaks = f.read()
    else:
        peaks = parse_elements(beam_energy)

    configure_cfg(cfg_path, peaks, energy_keV)

    file_path = os.path.join(DataDir, file_name)

    if not os.path.isfile(cfg_path):
        raise FileNotFoundError(f"File {cfg_path} does not exist")

    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"File {file_path} does not exist")

    b = McaAdvancedFitBatch.McaAdvancedFitBatch(
        cfg_path, [file_path], "out", 0, 250.0, selection=selection
    )
    # ProcessList method fits data and writes results to .h5 file
    b.processList()
    DatOut = os.path.join(ResultsDir, FilePrefix + ".h5")

    if not os.path.isfile(DatOut):
        raise FileNotFoundError(f"Results file {DatOut} could not be opened")

    peaks = parse_spec_fit(DatOut)

    ResultsFile = os.path.join(OutputDir, FilePrefix) + ".results.dat"
    results_txt = [f"{peak[0]} {peak[1]} {peak[2]}" for peak in peaks]
    results_txt = "\n".join(results_txt)
    with open(ResultsFile, "w") as f:
        f.write(results_txt)

    # Read and plot the spectrum data
    if h5py.is_hdf5(inputFile):
        # Read the calibration from config file
        calib = {}
        with open(cfg_path) as f:
            for line in f:
                # Uses regex to extract the zero and gain calibration values
                if match := re.match(r"^(zero|gain) = (\d+(?:\.\d+)?)$", line.strip()):
                    param, value = match.groups()
                    calib[param] = float(value)
        # Read in spectrum data counts
        with h5py.File(inputFile, "r") as hf:
            # Squeeze to collapse 1-length dimensions
            channel_counts = np.squeeze(hf[h5path])
        # Calculate channel energies (in eV) from calibration
        channel_energy = np.array(
            [
                (calib["zero"] + calib["gain"] * _i) * 1000
                for _i in range(len(channel_counts))
            ]
        )
        mca_out = os.path.splitext(inputFile)[0] + ".mca"
        # Create a .mca file (not used by this code but a more user friendly file format for PyMCA GUI)
        spectrum_to_mca(channel_counts, calib, mca_out, inputFile, selection)

    else:
        rawDatFile = os.path.splitext(inputFile)[0] + ".dat"
        spectrum_data = np.loadtxt(
            rawDatFile, delimiter="\t", skiprows=3, usecols=(0, 1)
        )
        channel_energy = spectrum_data[:, 0]
        channel_counts = spectrum_data[:, 1]

    pymca_output = parse_raw_fluoro(
        channel_energy, channel_counts, beam_energy, peaks, cutoff_offset
    )

    outFile = os.path.splitext(inputFile)[0] + ".png"

    plot_fluorescence_spectrum(
        outFile,
        inputFile,
        channel_energy,
        channel_counts,
        beam_energy,
        offset=cutoff_offset,
    )

    pure_path = PurePath(inputFile).parts
    currentvisit = pure_path[5]
    ScanName = os.path.basename(inputFile)
    RelNameHtml = os.path.join(*pure_path[6:])
    OutputDir = os.path.dirname(inputFile).replace(
        currentvisit, os.path.join(currentvisit, "jpegs")
    )
    ScanNumber = os.path.splitext(ScanName)[0]

    # Find crystal snapshot if it exists
    try:
        RelPNGfile = sorted(
            glob(os.path.join(OutputDir, ScanNumber) + "*[0-9].png"),
            key=os.path.getmtime,
        )
        pure_path = PurePath(RelPNGfile[0]).parts
        RelPNGfile = os.path.join(*pure_path[7:])
    except Exception:
        RelPNGfile = ""

    Path(OutputDir).mkdir(parents=True, exist_ok=True)

    HtmlFile = os.path.splitext(inputFile)[0] + ".html"

    timestamp = datetime.now().strftime("%a %b %-d, %Y - %T")

    HtmlFileContents = f"""<html><head><title>{ScanName}</title></head><body>
<style type="text/css">
table td {{ padding: 5px;}}
table.table2 {{ border-collapse: collapse;}}
table.table2 td {{ border-style: none; font-size: 80%; padding: 2px;}}
</style>
<table border width=640><tr><td align=center>{timestamp}<br/><a href="{RelNameHtml}">{ScanName}</a></td></tr>
<tr><td><table class="table2" width=100%>
<tr><th colspan=4 align=center>Fluorescence Spectrum</td></tr>
<tr><td>Beamline Energy:</td><td>{beam_energy}eV</td><td>Omega:</td><td>{omega}&deg;</td></tr>
<tr><td>Acq Time:</td><td>{acqTime}s</td><td>Trans:</td><td>{transmission}%</td></tr>
<tr><td>Sample Position:</td><td colspan=3>{samplexyz}</td></tr>
</table></td></tr>
<tr><td>Automated PyMca results<br /><pre>{pymca_output}</pre><a href="http://www.diamond.ac.uk/dms/MX/Common/Interpreting-AutoPyMCA/Interpreting%20AutoPyMCA.pdf">Guide to AutoPyMCA</a> (pdf)</td></tr>
<tr><td><img src="{outFile}" /></td></tr>
<tr><td><img src="jpegs/{RelPNGfile}" alt="Snapshot not taken" width=640 /></td></tr></table>
"""

    with open(HtmlFile, "w") as f:
        f.write(HtmlFileContents)
