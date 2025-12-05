"""Settings module for plotting."""

from pydantic import Field

from koopmans.pydantic import BaseModel


class PlottingConfig(BaseModel):
    """Model for the configuration of plotting settings."""

    degauss: float = Field(default=0.05, description="Gaussian broadening (in eV) for the DOS interpolation, as in QE.")
    nstep: int = Field(default=1000, description="Number of steps for the plot of the interpolated DOS.")
    Emin: float | None = Field(default=None, description="Minimum energy for the plot of the interpolated DOS.")
    Emax: float | None = Field(default=None, description="Maximum energy for the plot of the interpolated DOS.")
