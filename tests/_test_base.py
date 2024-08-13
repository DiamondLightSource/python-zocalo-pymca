import unittest
import os
import tempfile
import pathlib
import shutil
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pymca_zocalo.internals


TEST_DATA = {
    'dat': '20190513_11_26_46_Cu.dat',
    'mca': '20190513_11_26_46_Cu.mca',
    'cfg': 'pymca_new.cfg',
    }

class Base(object):
    def setUp(self):
        self.visitDir = self._get_visit_dir()
        self.visitDir.mkdir(parents=True)

        for f in TEST_DATA.values():
            shutil.copy2(os.path.join(os.path.dirname(os.path.abspath(__file__)), f), self.visitDir)

        self.SpectrumFile = os.path.join(self.visitDir, os.path.basename(TEST_DATA['mca']))
        self.energy = '14000.0' # eV
        self.transmission = '100.0'
        self.acqTime = '4.0'
        self.omega = '1.0'
        self.samplexyz = '1.0'
        self.xfeFluorescenceSpectrumID = '25'
        self.CFGFile = os.path.join(self.visitDir, os.path.basename(TEST_DATA['cfg']))
        self.DatFile = os.path.join(self.visitDir, os.path.basename(TEST_DATA['dat']))

        # append peaks to the cfg file
        peaks = pymca_zocalo.internals.parse_elements(self.energy)
        with open(self.CFGFile, 'a') as f:
            print(peaks, file = f)


    def test_run_auto_pymca(self):
        pymca_zocalo.internals.run_auto_pymca(self.SpectrumFile, self.energy, self.transmission, self.acqTime, self.CFGFile)

    def test_plot_fluorescence_spectrum(self):
        pymca_zocalo.internals.plot_fluorescence_spectrum(self.DatFile, self.omega, self.transmission, self.samplexyz, self.acqTime, self.energy, self.xfeFluorescenceSpectrumID, CFGFile=self.CFGFile)

        outFile = os.path.splitext(self.DatFile)[0] + '.png'
        self.assertTrue(os.path.exists(outFile))
        HtmlFile = os.path.splitext(self.DatFile)[0] + '.html'

class TestSimple(Base, unittest.TestCase):
    def _get_visit_dir(self):
        tempdir = tempfile.TemporaryDirectory(dir='/tmp', prefix='beamline-')
        return pathlib.Path(tempdir.name, 'data', 'year', 'visit')

class TestDepthOne(Base, unittest.TestCase):
    def _get_visit_dir(self):
        tempdir = tempfile.TemporaryDirectory(dir='/tmp', prefix='beamline-')
        return pathlib.Path(tempdir.name, 'data', 'year', 'visit', 'test')

class TestDepthThree(Base, unittest.TestCase):
    def _get_visit_dir(self):
        tempdir = tempfile.TemporaryDirectory(dir='/tmp', prefix='beamline-')
        return pathlib.Path(tempdir.name, 'data', 'year', 'visit', 'test1', 'test2', 'test3')
