from pydantic import Field, field_validator, model_validator
from typing_extensions import Self

from koopmans.pydantic import BaseModel


class MLConfig(BaseModel):
    train: bool = Field(False, description="train a machine learning model to predict the screening parameters")
    test: bool = Field(False, description="test the machine learning model")
    predict: bool = Field(False, description="use a machine learning model to predict the screening parameters")
    model_file: str | None = Field(
        None, description="JSON file containing the ML model information (generated from a prior training calculation)")
    n_max: int = Field(4, gt=0, description="The maximum expansion coefficient n for radial basis functions")
    l_max: int = Field(4, gt=0, description="The maximum angular expansion coefficient")
    r_min: float = Field(0.5, gt=0.0, description="The width of the narrowest radial basis function")
    r_max: float = Field(4.0, gt=0.0, description="The width of the broadest radial basis function")
    alphas_from_file: bool = Field(
        False, description="If True, read the screening coefficients from file instead of calculating them ab initio")
    train_on_the_fly: bool = Field(
        False, description="If True, the ML-model gets trained after the calculation of each orbital. If False, the ML-model gets trained at the end of each snapshot")
    occ_and_emp_together: bool = Field(
        True, description="If True, there will be one ML model for both occupied and empty states. If False, there will be one ML Model for occupied states and one for empty states")
    estimator: str = Field('ridge_regression', description="What to use as the estimator for the ML model")
    descriptor: str = Field('orbital_density', description="What to use as the descriptor for the ML model")

    @model_validator(mode='after')
    def mutually_exclusive_modes(self) -> Self:
        if [self.predict, self.train, self.test].count(True) > 1:
            raise ValueError(
                'Training, testing, and using the ML model are mutually exclusive; change `ml:predict` '
                '/`ml:train`/`ml:test` so that at most one is `True`')
        return self

    @field_validator('r_min', mode='after')
    @classmethod
    def check_small_rmin(cls, v: float) -> float:
        if v < 0.5:
            warn(
                "Small values of `r_min` (<0.5) can lead to problems in the construction of the radial basis. "
                f"The provided value is `r_min = {v}`.")
        return v

    @model_validator(mode='after')
    def check_rmin_less_than_rmax(self) -> Self:
        if not self.r_min < self.r_max:
            raise ValueError(f"`r_min` is larger or equal to `r_max = {self.r_max}`.")
