from setuptools import setup

setup(
    name="pymca-zocalo",
    version="0.2.1",
    description="PyMca components for automated data processing with Zocalo at Diamond Light Source",
    author="Tom Schoonjans",
    author_email="DataAnalysis@diamond.ac.uk",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    packages=["pymca_zocalo"],
    install_requires=[
        "pymca",
        "setuptools",
        "workflows>=1.7",
        "xraylib",
        "zocalo",
    ],
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
