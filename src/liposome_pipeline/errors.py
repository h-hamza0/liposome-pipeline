class PipelineError(Exception):
    """Base class for all errors raised by this package."""


class ConfigError(PipelineError):
    """The pipeline configuration file is malformed or incomplete."""


class GromacsError(PipelineError):
    """A GROMACS invocation failed."""


class SlurmError(PipelineError):
    """A SLURM submission failed or timed out."""
