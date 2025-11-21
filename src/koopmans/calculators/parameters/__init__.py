from koopmans.pydantic import BaseModel
from koopmans import settings
from pydantic import Field, model_validator
from typing_extensions import Self
from typing import Any
from .pw import PWInputParameters
from .wannier90 import Wannier90InputParameters

class CalculatorsInputParameters(BaseModel):
    """A collection of calculator configuration dictionaries."""

    kcp: settings.KoopmansCPSettingsDict = Field(default_factory=settings.KoopmansCPSettingsDict)
    kcw_ham: settings.KoopmansHamSettingsDict = Field(default_factory=settings.KoopmansHamSettingsDict)
    kcw_screen: settings.KoopmansScreenSettingsDict = Field(default_factory=settings.KoopmansScreenSettingsDict)
    kcw_wannier: settings.Wann2KCSettingsDict = Field(default_factory=settings.Wann2KCSettingsDict)
    ph: settings.PhSettingsDict = Field(default_factory=settings.PhSettingsDict)
    projwfc: settings.ProjwfcSettingsDict = Field(default_factory=settings.ProjwfcSettingsDict)
    pw: PWInputParameters = Field(default_factory=PWInputParameters)
    pw2wannier: settings.PW2WannierSettingsDict = Field(default_factory=settings.PW2WannierSettingsDict)
    wann2kcp: settings.Wann2KCPSettingsDict = Field(default_factory=settings.Wann2KCPSettingsDict)
    ui: settings.UnfoldAndInterpolateSettingsDict = Field(default_factory=settings.UnfoldAndInterpolateSettingsDict)
    ui_occ: settings.UnfoldAndInterpolateSettingsDict = Field(default_factory=settings.UnfoldAndInterpolateSettingsDict)
    ui_emp: settings.UnfoldAndInterpolateSettingsDict = Field(default_factory=settings.UnfoldAndInterpolateSettingsDict)
    wjl: settings.WannierJLSettingsDict = Field(default_factory=settings.WannierJLSettingsDict)
    w90: Wannier90InputParameters | None = None
    w90_up: Wannier90InputParameters | None = None
    w90_down: Wannier90InputParameters | None = None

    @model_validator(mode='after')
    def w90_exclusivity(self) -> Self:
        if (self.w90_up is None) != (self.w90_down is None):
            raise ValueError("Both 'w90_up' and 'w90_down' must be specified together.")
        elif (self.w90 is None) != (self.w90_up is None):
            raise ValueError("Either specify 'w90' for non-spin-polarized calculations or both 'w90_up' and 'w90_down' for spin-polarized calculations.")
        return self
    
    @model_validator(mode='before')
    @classmethod
    def flatten_ui_occ_and_emp_subdictionaries(cls, data: dict[str, Any]) -> dict[str, Any]:
        if data['ui']:
            if 'occ' not in data['ui'] and 'emp' not in data['ui']:
                return data
            universal_settings = {k: v for k, v in data['ui'].items() if k not in ['occ', 'emp']}
            for filling in ['occ', 'emp']:
                # Add any keywords in the spin subdictionary
                data[f'ui_{filling}'] = data['ui'].pop(filling, {})
                # Add any keywords in the main dictionary
                data[f'ui_{filling}'].update(universal_settings)
        return data

    @model_validator(mode='before')
    @classmethod
    def flatten_w90_up_down_subdictionaries(cls, data: dict[str, Any]) -> dict[str, Any]:
        """Allow w90_up and w90_down to be specified as subdictionaries of w90.

        e.g.
        ```
        "w90": {
            "num_iter": 100,
            "up": {
                "projections": ...
            },
            "down": {
                "projections": ...
            }
        ```
        
        """
        if data['w90']:
            if 'up' not in data['w90'] and 'down' not in data['w90']:
                return data
            universal_settings = {k: v for k, v in data['w90'].items() if k not in ['up', 'down']}
            for spin in ['up', 'down']:
                # Add any keywords in the spin subdictionary
                data[f'w90_{spin}'] = data['w90'].pop(spin, {})
                # Add any keywords in the main dictionary
                data[f'w90_{spin}'].update(universal_settings)
            data['w90'] = None
        return data