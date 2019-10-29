from __future__ import absolute_import, division, print_function

import errno
import os
import re
import sys
import traceback

import workflows.recipe
from workflows.services.common_service import CommonService
from .internals import plot_fluorescence_spectrum

class DLSPyMcaFitter(CommonService):
    """A service that takes an XRF dataset and sends it to PyMca for fitting"""

    _service_name = "DLS PyMca Fitter"

    _logger_name = "dlstbx.services.pymca_fitter"

    def initializing(self):
        """Subscribe to a queue. Received messages must be acknowledged."""
        self.log.info("PyMca fitter service starting")
        workflows.recipe.wrap_subscribe(
            self._transport,
            "pymca.fitter",
            self.pymca_fitter_call,
            acknowledgement=True,
            log_extender=self.extend_log
        )

    def pymca_fitter_call(self, rw, header, message):
        """Call dispatcher"""
        args = list()
        args.append(rw.recipe_step.get("parameters", {}).get("inputFile"))
        args.append(rw.recipe_step.get("parameters", {}).get("omega"))
        args.append(rw.recipe_step.get("parameters", {}).get("transmission"))
        args.append(rw.recipe_step.get("parameters", {}).get("samplexyz"))
        args.append(rw.recipe_step.get("parameters", {}).get("acqTime"))
        args.append(rw.recipe_step.get("parameters", {}).get("energy"))
        args.append(rw.recipe_step.get("parameters", {}).get("xfeFluorescenceSpectrumID"))

        self.log.debug("Commands: %s", ' '.join(args))
        try:
            plot_fluorescence_spectrum(*args)
        except Exception as e:
            exc_info = sys.exc_info()
            self.log.error("Error running PyMca: {}\n{}".format(str(e), traceback.format_exception(*exc_info)))

            rw.transport.nack(header)
            return
        self.log.info("%s was successfully processed", rw.recipe_step.get("parameters", {}).get("inputFile"))
        rw.transport.ack(header)

