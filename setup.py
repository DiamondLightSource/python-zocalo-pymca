from setuptools import setup

# these lines allow the version to be specified in Makefile.private
import os

version = os.environ.get("MODULEVER", "0.0")

setup(
    name="pymca-zocalo",
    version=version,
    description="Module",
    author="Tom Schoonjans",
    author_email="Tom.Schoonjans@diamond.ac.uk",
    packages=["pymca_zocalo"],
    install_requires=["workflows>=1.7", "zocalo", "procrunner", "setuptools", "PyMca5"],
    entry_points={
        "workflows.services": [
            "DLSPyMcaFitter = pymca_zocalo:DLSPyMcaFitter",
        ],
    },
    zip_safe=False,
)
