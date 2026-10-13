# perturbation-status-ps

Independent Omnibenchmark wrapper for the frozen, config-driven PS
perturbation-status estimator. The module rebuilds its own input from the
canonical H5AD and does not depend on the Mixscape/Mixscale shared signature.
It accepts the stage-level shared-artifact flag for Omnibenchmark contract
compatibility but deliberately does not read that artifact.

The reviewed PS implementation and dataset configuration are versioned under
`implementation/`, while `environment/ps.def` fixes the R/Python environment
and upstream scMAGeCK and presto revisions. The module does not access the
original benchmark development directory at runtime.
