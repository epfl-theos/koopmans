from typing import Any

from ._utils import SettingsDict


class WannierJLSettingsDict(SettingsDict):
    def __init__(self, **kwargs) -> None:

        super().__init__(valid=['task', 'nval', 'outdir-val', 'outdir-cond', 'run-disentangle', 'run-optrot',
                                'run-maxloc', 'rotate-unk', 'binary'],
                         defaults={'task': 'splitvc', 'outdir-val': 'val', 'outdir-cond': 'cond'},
                         are_paths=['outdir-val', 'outdir-cond'],
                         **kwargs)

    @property
    def _other_valid_keywords(self):
        return []

    def __setitem__(self, key: str, value: Any):
        if key == 'task':
            if value != 'splitvc':
                raise NotImplementedError('task != splitvc is not yet supported')
        return super().__setitem__(key, value)
