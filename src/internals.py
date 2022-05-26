import os
from pathlib import PurePath, Path
import xraylib as xrl
import shutil
import getpass
import smtplib
from email.message import EmailMessage

from PyMca5.PyMca import McaAdvancedFitBatch

import h5py

from glob import glob

import subprocess

from datetime import datetime

EDGE_MAPPER = {
    'K': xrl.K_SHELL,
    'L': xrl.L3_SHELL,
}

LINE_MAPPER = {
    'K': xrl.KL3_LINE,
    'L': xrl.L3M5_LINE,
}

# edge is here the beam energy...
def parse_raw_fluoro(fn, edge, trans, acqTime):
    rv = []

    with open(fn) as f:
        total_count = 0
        background_count = 0
        for i,l in enumerate(f):
            if i > 2:
                e, cts = l.strip().split()
                total_count += float(cts)
                background_count += min(float(cts), 2)

                # crude but works I guess...
                if float(e) > (float(edge) - 1250.0):
                    break

    if (total_count - background_count) > 100:
        FilePrefix = os.path.splitext(os.path.basename(fn))[0]
        pure_path = PurePath(fn).parts
        VisitDir = pure_path[:6]
        RestOfFilename = os.path.join(*pure_path[6:])
        RestOfDirs = os.path.dirname(RestOfFilename)
        OutputDir = os.path.join(*VisitDir, 'processed/pymca', RestOfDirs)
        elf = os.path.join(OutputDir, FilePrefix) + '.results.dat'

        els = []
        # contains lines like
        # Cu-K 7929.020000 90.714300
        with open(elf) as f:
            rv.append("Element\tCounts\t%age\tExpected Emission Energies")
            for i,l in enumerate(f):
                el, pk, conf = l.strip().split()
                symbol, edge = el.split('-')
                Z = xrl.SymbolToAtomicNumber(symbol)
                edgesEnergies = "<b>{:g}</b>,{:g}".format(xrl.LineEnergy(Z, LINE_MAPPER[edge]) * 1000.0, xrl.EdgeEnergy(Z, EDGE_MAPPER[edge]) * 1000.0)
                if float(pk) >= 100000:
                    counts = int(float(pk))
                else:
                    counts = round(float(pk),1)
                rv.append("{}\t{:g}\t{:g}\t{}".format(el, counts, round(100*float(pk)/total_count, 1), edgesEnergies))
                if i == 5:
                    break
    else:
        rv.append("No fluorescence peaks detected, try a higher transmission")

    rv.append("\nCounts (total): {:g} (background): {:g}".format(total_count, background_count))
    return rv

# Parses a file and sorts spectral peaks with largest area first
# and prints them to stdout
def parse_spec_fit(name):
    f = h5py.File(name, 'r')
    parameters = f[list(f.keys())[0] + '/xrf_fit/results/parameters']
    all_fit = filter(lambda name: not name.startswith('Scatter') and not name.endswith('_errors'), parameters)

    def mapper(name):
        return [name.replace('_', '-'), parameters[name][0, 0], parameters[name + '_errors'][0, 0]]

    mapped = map(mapper, all_fit)

    peaks = sorted(mapped, key = lambda p: p[1], reverse=True)

    return peaks

def parse_elements(energy): 
    """ energy: expressed in eV!"""
    elements = {
        'Ti' : 'K', \
        'V' : 'K', \
        'Cr' : 'K', \
        'Mn' : 'K', \
        'Fe' : 'K', \
        'Co' : 'K', \
        'Ni' : 'K', \
        'Cu' : 'K', \
        'Zn' : 'K', \
        'As' : 'K', \
        'Se' : 'K', \
        'Br' : 'K', \
        'Sr' : 'K', \
        'Mo' : 'K', \
        'I' : 'L', \
        'Xe' : 'L', \
        'Gd' : 'L', \
        'W' : 'L', \
        'Os' : 'L', \
        'Ir' : 'L', \
        'Pt' : 'L', \
        'Au' : 'L', \
        'Hg' : 'L', \
        'Pb' : 'L'}

    def _check_edge(element):
        symbol, edge = element
        Z = xrl.SymbolToAtomicNumber(symbol)
        shell = EDGE_MAPPER[edge]
        return float(energy) > xrl.EdgeEnergy(Z, shell) * 1000.0

    return '\n'.join(map(lambda element : '{} = {}'.format(*element), filter(_check_edge, elements.items())))

def find_cut_off_energy(inputFile, cutoffenergy):
    with open(inputFile, 'r') as f:
        for i, line in enumerate(f):
            try:
                energy = float(line.split()[0])
                if energy > float(cutoffenergy):
                    return i
            except:
                pass

def run_auto_pymca(SpectrumFile, energy, transmission, acqTime, CFGFile=None, peaksFile=None):
    FilePrefix = os.path.splitext(os.path.basename(SpectrumFile))[0]
    pure_path = PurePath(SpectrumFile).parts
    VisitDir = pure_path[:6]
    RestOfFilename = os.path.join(*pure_path[6:])
    RestOfDirs = os.path.dirname(RestOfFilename)
    BEAMLINE = pure_path[2]

    peaks = None

    if CFGFile is None:
        CFGFile = os.path.join('/dls_sw', BEAMLINE, 'software/pymca/pymca_new.cfg')
        peaks = parse_elements(energy)
    elif peaksFile is not None:
        with open(peaksFile) as f:
            peaks = f.read()


    if not os.path.isfile(CFGFile):
        CFGFile = '/dls_sw/i03/software/pymca/pymca.cfg'

    energy_keV = float(energy) / 1000.0

    OutputDir = os.path.join(*VisitDir, 'processed/pymca', RestOfDirs)
    DataDir = os.path.join(OutputDir, 'data')
    ResultsDir = os.path.join(OutputDir, 'out')

    Path(OutputDir).mkdir(parents=True, exist_ok=True)
    Path(DataDir).mkdir(parents=True, exist_ok=True)
    Path(ResultsDir).mkdir(parents=True, exist_ok=True)

    os.chdir(OutputDir)
    shutil.copy2(SpectrumFile, DataDir)

    with open(CFGFile, 'r') as f1, open(os.path.join(OutputDir, FilePrefix + '.cfg'), 'w') as f2:
        for line in f1:
            line.replace('20000.0', str(energy_keV))
            print(line, file=f2)

        if peaks is not None:
            print(peaks, file=f2)

    if not os.path.isfile(os.path.join(OutputDir, FilePrefix + '.cfg')):
        raise FileNotFoundError(os.path.join(OutputDir, FilePrefix + '.cfg'))

    if not os.path.isfile(os.path.join(DataDir, FilePrefix + '.mca')):
        raise FileNotFoundError(os.path.join(DataDir, FilePrefix + '.mca'))

    b = McaAdvancedFitBatch.McaAdvancedFitBatch(os.path.join(OutputDir, FilePrefix + '.cfg'), os.path.join(DataDir, FilePrefix + '.mca'), 'out', 0, 250.0)
    b.processList()

    # the results are written to an HDF5 file
    DatOut = os.path.join(ResultsDir, FilePrefix + '.h5')

    if not os.path.isfile(DatOut):
        raise FileNotFoundError('{} could not be opened'.format(DatOut))

    peaks = parse_spec_fit(DatOut)

    ResultsFile = os.path.join(OutputDir, FilePrefix) + '.results.dat'

    with open(ResultsFile, 'w') as f:
        for p in peaks:
            print('{} {} {}'.format(p[0], p[1], p[2]), file=f)

    rawDatFile = os.path.splitext(SpectrumFile)[0] + '.dat'

    html = parse_raw_fluoro(rawDatFile, energy, transmission, acqTime)

    if getpass.getuser() == 'gda2' or getpass.getuser() == 'awf63395':
        with open(ResultsFile) as f:
            msg = EmailMessage()
            msg.set_content(f.read())

        msg['Subject'] = "pymca zocalo: {}".format(SpectrumFile)
        msg['From'] = "donotreply@pymca-zocalo.diamond.ac.uk"
        if getpass.getuser() == 'gda2':
            msg['To'] = 'Dave.Hall@diamond.ac.uk'
        else:
            msg['To'] = 'Tom.Schoonjans@diamond.ac.uk'

        with smtplib.SMTP('localhost') as smtp:
            smtp.send_message(msg)

    return html


def plot_fluorescence_spectrum(inputFile, omega, transmission, samplexyz, acqTime, energy, xfeFluorescenceSpectrumID, CFGFile=None, peaksFile=None):
    if not inputFile.endswith('.dat'):
        raise ValueError('inputFile must end with .dat')
    MCAFile = os.path.splitext(inputFile)[0] + '.mca'

    outputPymca = run_auto_pymca(MCAFile, energy, transmission, acqTime, CFGFile=CFGFile, peaksFile=peaksFile)
    outputPymca = '\n'.join(outputPymca)

    if getpass.getuser() == 'gda2' or getpass.getuser() == 'awf63395':
        msg = EmailMessage()
        msg.set_content(outputPymca)

        msg['Subject'] = f"pymca zocalo: {inputFile} trans: {transmission} time {acqTime} energy {energy}"
        msg['From'] = "donotreply@pymca-zocalo.diamond.ac.uk"
        if getpass.getuser() == 'gda2':
            msg['To'] = 'juan.sanchez-weatherby@diamond.ac.uk, mark.williams@diamond.ac.uk'
        else:
            msg['To'] = 'Tom.Schoonjans@diamond.ac.uk'

        with smtplib.SMTP('localhost') as smtp:
            smtp.send_message(msg)

    outFile = os.path.splitext(inputFile)[0] + '.png'
    cutoffenergy = float(energy) - 1000.0

    cutoffline = find_cut_off_energy(inputFile, cutoffenergy)

    gnuplot_input = \
f"""set term png size 800,600
set title "Fluorescence Spectrum {inputFile}"
set xlabel "Energy (eV)"
set ylabel "Number of Counts"
set xrange [0:{energy}]
set mxtics 4
set mytics 4
set grid xtics ytics mxtics mytics
set style arrow 1 nohead lt 4 lw 5
set style arrow 2 nohead lt 39 lw 1
set arrow from {cutoffenergy},graph(0) to {cutoffenergy},graph(1) as 2
set out '{outFile}'
set nokey
set style line 1 lt 1 
set style line 2 lt 39
plot '{inputFile}' every ::0::{cutoffline} using 1:2 with lines ls 1, '{inputFile}' every ::{cutoffline} using 1:2 with lines ls 2
"""

    process_rv = subprocess.run(['/usr/bin/gnuplot'], input=gnuplot_input, universal_newlines=True, encoding='utf8')

    pure_path = PurePath(inputFile).parts
    currentvisit = pure_path[5]
    ScanName = os.path.basename(inputFile)
    RelNameHtml = os.path.join(*pure_path[6:])
    OutputDir = os.path.dirname(inputFile).replace(currentvisit, os.path.join(currentvisit, 'jpegs'))
    ScanNumber = os.path.splitext(ScanName)[0]

    # RelPNGfile is allowed to be empty apparently...
    try:
        RelPNGfile = sorted(glob(os.path.join(OutputDir, ScanNumber) + '*[0-9].png'), key=os.path.getmtime)
        pure_path = PurePath(RelPNGfile[0]).parts
        RelPNGfile = os.path.join(*pure_path[7:])
    except Exception as e:
        # print('RelPNGfile exception: ' + str(e))
        RelPNGfile = ''

    Path(OutputDir).mkdir(parents=True, exist_ok=True)

    HtmlFile = os.path.splitext(inputFile)[0] + '.html'

    timestamp = datetime.now().strftime('%a %b %-d, %Y - %T')

    HtmlFileContents = \
f'''<html><head><title>{ScanName}</title></head><body>
<style type="text/css">
table td {{ padding: 5px;}}
table.table2 {{ border-collapse: collapse;}}
table.table2 td {{ border-style: none; font-size: 80%; padding: 2px;}}
</style>
<table border width=640><tr><td align=center>{timestamp}<br/><a href="{RelNameHtml}">{ScanName}</a></td></tr>
<tr><td><table class="table2" width=100%>
<tr><th colspan=4 align=center>Fluorescence Spectrum</td></tr>
<tr><td>Beamline Energy:</td><td>{energy}eV</td><td>Omega:</td><td>{omega}&deg;</td></tr>
<tr><td>Acq Time:</td><td>{acqTime}s</td><td>Trans:</td><td>{transmission}%</td></tr>
<tr><td>Sample Position:</td><td colspan=3>{samplexyz}</td></tr>
</table></td></tr>
<tr><td>Automated PyMca results<br /><pre>{outputPymca}</pre><a href="http://www.diamond.ac.uk/dms/MX/Common/Interpreting-AutoPyMCA/Interpreting%20AutoPyMCA.pdf">Guide to AutoPyMCA</a> (pdf)</td></tr>
<tr><td><img src="{outFile}" /></td></tr>
<tr><td><img src="jpegs/{RelPNGfile}" alt="Snapshot not taken" width=640 /></td></tr></table>
'''

    with open(HtmlFile, 'w') as f:
        f.write(HtmlFileContents)
