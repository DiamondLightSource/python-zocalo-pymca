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

EDGE_MAPPER = {
    "K": xrl.K_SHELL,
    "L": xrl.L3_SHELL,
}

LINE_MAPPER = {
    "K": xrl.KL3_LINE,
    "L": xrl.L3M5_LINE,
}


# edge is here the beam energy...
def parse_raw_fluoro(counts, channel_energy, beam_energy, results_file):
    rv = []
    # Get total counts and background counts up to a cutoff energy
    i_cutoff = np.argmax(channel_energy > beam_energy - 1250.0)
    total_count = np.sum(counts[0:i_cutoff])
    # Maximum of 2 counts per channel contribute to background
    background_count = np.sum(np.minimum(2, counts[0:i_cutoff]))

    if (total_count - background_count) > 100:
        # Results file contains lines like
        # Cu-K 7929.020000 90.714300
        with open(results_file) as f:
            rv.append("Element\tCounts\t%age\tExpected Emission Energies")
            for _i, line in enumerate(f):
                el, pk, conf = line.strip().split()
                symbol, edge = el.split("-")
                Z = xrl.SymbolToAtomicNumber(symbol)
                edgesEnergies = "<b>{:g}</b>,{:g}".format(
                    xrl.LineEnergy(Z, LINE_MAPPER[edge]) * 1000.0,
                    xrl.EdgeEnergy(Z, EDGE_MAPPER[edge]) * 1000.0,
                )
                if float(pk) >= 100000:
                    pk_counts = int(float(pk))
                else:
                    pk_counts = round(float(pk), 1)
                rv.append(
                    "{}\t{:g}\t{:g}\t{}".format(
                        el,
                        pk_counts,
                        round(100 * float(pk) / total_count, 1),
                        edgesEnergies,
                    )
                )
                if _i == 5:
                    break
    else:
        rv.append("No fluorescence peaks detected, try a higher transmission")

    rv.append(
        "\nCounts (total): {:g} (background): {:g}".format(
            total_count, background_count
        )
    )
    return rv


# Parses a file and sorts spectral peaks with largest area first
# and prints them to stdout
def parse_spec_fit(name):
    f = h5py.File(name, "r")
    parameters = f[list(f.keys())[0] + "/xrf_fit/results/parameters"]
    all_fit = filter(
        lambda name: not name.startswith("Scatter") and not name.endswith("_errors"),
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

    return peaks


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

    def _check_edge(element):
        symbol, edge = element
        Z = xrl.SymbolToAtomicNumber(symbol)
        shell = EDGE_MAPPER[edge]
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


def plot_fluorescence_spectrum(
    outFile,
    inputFile,
    channel_energy,
    channel_counts,
    beam_energy,
):
    cutoff_energy = beam_energy - 1000
    i_cutoff = np.argmax(channel_energy > cutoff_energy)

    # Create plot of spectrum split by the cutoff
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
        root_group = "entry"
        y_data_path = "/data/data"
        selection = {"entry": root_group, "y": y_data_path}

    FilePrefix = os.path.splitext(os.path.basename(inputFile))[0]
    file_name = os.path.basename(inputFile)
    pure_path = PurePath(inputFile).parts
    VisitDir = pure_path[:6]
    RestOfFilename = os.path.join(*pure_path[6:])
    RestOfDirs = os.path.dirname(RestOfFilename)
    BEAMLINE = pure_path[2]

    peaks = None

    if CFGFile is None:
        CFGFile = os.path.join("/dls_sw", BEAMLINE, "software/pymca/pymca_new.cfg")
        peaks = parse_elements(beam_energy)
    elif peaksFile is not None:
        with open(peaksFile) as f:
            peaks = f.read()

    if not os.path.isfile(CFGFile):
        raise FileNotFoundError(f"Config file '{CFGFile}' does not exist")

    energy_keV = float(beam_energy) / 1000.0

    OutputDir = os.path.join(*VisitDir, "processed/pymca", RestOfDirs)
    DataDir = os.path.join(OutputDir, "data")
    ResultsDir = os.path.join(OutputDir, "out")

    Path(OutputDir).mkdir(parents=True, exist_ok=True)
    Path(DataDir).mkdir(parents=True, exist_ok=True)
    Path(ResultsDir).mkdir(parents=True, exist_ok=True)

    os.chdir(OutputDir)
    shutil.copy2(inputFile, DataDir)

    with open(CFGFile, "r") as f1, open(
        os.path.join(OutputDir, FilePrefix + ".cfg"), "w"
    ) as f2:
        for line in f1:
            # Replace the default energy with the beam energy
            if match := re.match(r"energy = (\d+(\.\d+)?),", line):
                line = line.replace(match.group(1), str(energy_keV))
            print(line, file=f2)

        if peaks is not None:
            print(peaks, file=f2)

    file_path = os.path.join(DataDir, file_name)
    cfg_path = os.path.join(OutputDir, FilePrefix + ".cfg")

    if not os.path.isfile(cfg_path):
        raise FileNotFoundError(f"File {cfg_path} does not exist")

    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"File {file_path} does not exist")

    b = McaAdvancedFitBatch.McaAdvancedFitBatch(
        cfg_path, file_path, "out", 0, 250.0, selection=selection
    )
    # ProcessList method fits data and writes results to .h5 file
    b.processList()
    DatOut = os.path.join(ResultsDir, FilePrefix + ".h5")

    if not os.path.isfile(DatOut):
        raise FileNotFoundError(f"Results file {DatOut} could not be opened")

    peaks = parse_spec_fit(DatOut)

    ResultsFile = os.path.join(OutputDir, FilePrefix) + ".results.dat"

    with open(ResultsFile, "w") as f:
        for p in peaks:
            print("{} {} {}".format(p[0], p[1], p[2]), file=f)

    if inputFile.endswith(".h5"):
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
            channel_counts = np.squeeze(hf["entry/data/data"])
        # Calculate channel energies (in eV) from calibration
        channel_energy = np.array(
            [
                calib["zero"] + calib["gain"] * _i * 1000
                for _i in range(len(channel_counts))
            ]
        )

    else:
        rawDatFile = os.path.splitext(inputFile)[0] + ".dat"
        spectrum_data = np.loadtext(rawDatFile, delimiter=" ", skiprows=3)
        channel_energy = spectrum_data[:, 0]
        channel_counts = spectrum_data[:, 1]

    pymca_output = parse_raw_fluoro(
        channel_energy, channel_counts, beam_energy, ResultsFile
    )
    pymca_output = "\n".join(pymca_output)

    outFile = os.path.splitext(inputFile)[0] + ".png"

    plot_fluorescence_spectrum(
        outFile,
        inputFile,
        channel_energy,
        channel_counts,
        beam_energy,
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
