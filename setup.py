from setuptools import find_packages, setup

setup(
    name="pymca-zocalo",
    version="1.0.0",
    description="PyMca components for automated data processing with Zocalo at Diamond Light Source",
    author="Tom Schoonjans",
    author_email="DataAnalysis@diamond.ac.uk",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3.11",
    ],
    packages=find_packages("src"),
    package_dir={"": "src"},
    include_package_data=True,
    package_data={"": ["*.cfg"]},
    install_requires=[
        # pymca and xraylib required but listing here causes error when zocalo resolves environment
        # "pymca",
        "workflows>=1.7",
        # "xraylib",
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
