# perturbation-status-ps

Independent Omnibenchmark wrapper for the frozen, config-driven PS
perturbation-status estimator. The module rebuilds its own input from the
canonical H5AD and does not depend on the Mixscape/Mixscale shared signature.
It accepts the stage-level shared-artifact flag for Omnibenchmark contract
compatibility but deliberately does not read that artifact.
