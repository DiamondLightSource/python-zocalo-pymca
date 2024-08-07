from __future__ import absolute_import, division, print_function

import workflows.recipe
from workflows.services.common_service import CommonService

from .internals import run_auto_pymca


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
            log_extender=self.extend_log,
        )

    def pymca_fitter_call(self, rw, header, message):
        """Call dispatcher"""
        params = rw.recipe_step.get("parameters", {}).copy()

        # Deal with un-populated placeholders with None
        req_placeholders = [
            "{inputFile}",
            "{energy}",
        ]

        opt_placeholders = [
            "{omega}",
            "{transmission}",
            "{samplexyz}",
            "{acqTime}",
            "{cfgFile}",
            "{h5path}",
        ]

        for key, value in params.items():
            if value in req_placeholders:
                self.log.error(
                    f"Error running PyMca - Missing rquired parameter: '{key}'"
                )
                rw.transport.ack(header)
                return
            if value in opt_placeholders:
                params[key] = None

        self.log.debug(f"pymca_fitter running with params {params}")
        try:
            run_auto_pymca(
                params["inputFile"],
                params["omega"],
                params["transmission"],
                params["samplexyz"],
                params["acqTime"],
                float(params["energy"]),
                cfgFile=params.get("cfgFile"),
                peaksFile=params.get("peaksFile"),
                h5path=params.get("h5path"),
            )
        except Exception as e:
            self.log.warning(f"Error running PyMca: {e}", exc_info=True)
            rw.transport.ack(header)
            return
        self.log.info(
            "%s was successfully processed",
            rw.recipe_step.get("parameters", {}).get("inputFile"),
        )
        rw.transport.ack(header)
