from setuptools import setup

setup(
    name="pymca-zocalo",
    version="0.1.0",
    description="PyMca components for automated data processing with Zocalo at Diamond Light Source",
    author="Tom Schoonjans",
    author_email="Tom.Schoonjans@diamond.ac.uk",
    packages=["pymca_zocalo"],
    install_requires=["workflows>=1.7", "zocalo", "setuptools", "PyMca5"],
    entry_points={
        "workflows.services": [
            "DLSPyMcaFitter = pymca_zocalo:DLSPyMcaFitter",
        ],
    },
    zip_safe=False,
    license="BSD license",
    test_suite="tests",
    url="https://github.com/DiamondLightSource/python-zocalo-pymca",
)
