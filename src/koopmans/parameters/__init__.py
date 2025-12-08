from typing import Any

from pydantic import Field, model_validator
from typing_extensions import Self

from koopmans import settings
from koopmans.parameters.ui import UnfoldAndInterpolateConfig
from koopmans.pydantic import BaseModel

from .kcp import KCPInputParameters
from .pw import PWInputParameters
from .pw2wannier90 import PW2Wannier90InputParameters
from .wannier90 import Wannier90InputParameters


class CalculatorsInputParameters(BaseModel):
    """A collection of calculator configuration dictionaries."""

    kcp: KCPInputParameters = Field(default_factory=lambda: KCPInputParameters())
    kcw_ham: settings.KoopmansHamSettingsDict = Field(default_factory=settings.KoopmansHamSettingsDict)
    kcw_screen: settings.KoopmansScreenSettingsDict = Field(default_factory=settings.KoopmansScreenSettingsDict)
    kcw_wannier: settings.Wann2KCSettingsDict = Field(default_factory=settings.Wann2KCSettingsDict)
    ph: settings.PhSettingsDict = Field(default_factory=settings.PhSettingsDict)
    projwfc: settings.ProjwfcSettingsDict = Field(default_factory=settings.ProjwfcSettingsDict)
    pw: PWInputParameters = Field(default_factory=lambda: PWInputParameters())
    pw2wannier90: PW2Wannier90InputParameters = Field(default_factory=lambda: PW2Wannier90InputParameters())
    wann2kcp: settings.Wann2KCPSettingsDict = Field(default_factory=settings.Wann2KCPSettingsDict)
    ui: UnfoldAndInterpolateConfig = Field(default_factory=lambda: UnfoldAndInterpolateConfig())
    ui_occ: UnfoldAndInterpolateConfig = Field(default_factory=lambda: UnfoldAndInterpolateConfig())
    ui_emp: UnfoldAndInterpolateConfig = Field(default_factory=lambda: UnfoldAndInterpolateConfig())
    wjl: settings.WannierJLSettingsDict = Field(default_factory=settings.WannierJLSettingsDict)
    wannier90: Wannier90InputParameters | None = None
    wannier90_up: Wannier90InputParameters | None = None
    wannier90_down: Wannier90InputParameters | None = None

    @model_validator(mode='after')
    def wannier90_exclusivity(self) -> Self:
        if (self.wannier90_up is None) != (self.wannier90_down is None):
            raise ValueError("Both 'wannier90_up' and 'wannier90_down' must be specified together.")
        elif (self.wannier90 is None) != (self.wannier90_up is None):
            raise ValueError(
                "Either specify 'wannier90' for non-spin-polarized calculations or both 'wannier90_up' and 'wannier90_down' for spin-polarized calculations.")
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
    def flatten_wannier90_up_down_subdictionaries(cls, data: dict[str, Any]) -> dict[str, Any]:
        """Allow wannier90_up and wannier90_down to be specified as subdictionaries of wannier90.

        e.g.
        ```
        "wannier90": {
            "num_iter": 100,
            "up": {
                "projections": ...
            },
            "down": {
                "projections": ...
            }
        ```

        """
        if data['wannier90']:
            if 'up' not in data['wannier90'] and 'down' not in data['wannier90']:
                return data
            universal_settings = {k: v for k, v in data['wannier90'].items() if k not in ['up', 'down']}
            for spin in ['up', 'down']:
                # Add any keywords in the spin subdictionary
                data[f'wannier90_{spin}'] = data['wannier90'].pop(spin, {})
                # Add any keywords in the main dictionary
                data[f'wannier90_{spin}'].update(universal_settings)
            data['wannier90'] = None
        return data
