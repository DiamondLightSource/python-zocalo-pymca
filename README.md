## PyMca-Zocalo: PyMca components for automated data processing with Zocalo at Diamond Light Source
![Build Status](https://github.com/DiamondLightSource/python-zocalo-pymca/actions/workflows/tests.yml/badge.svg?branch=main)

Repository containing the pymca fitter zocalo service, that automates the fitting of XRF spectra to identify elements present in a sample, along with their relative abundance.

The service is deployed in the zocalo namespace kubernetes deployment that is managed in [container-tools](https://gitlab.diamond.ac.uk/scisoft/container-tools). The container image is created via a GitHub action in this repository, which automatically creates and and pushes an image with a unique tag when commits are made to the main branch.

Note to non-members of the Diamond Light Source organisation: While this repository is public and a dockerfile exists, this service is specific to the Diamond Light Source infrastructure and is not a plug and play service. It requires an active installation of rabbitmq and Zocalo and it is tailored to the structure of the DLS file-system.


