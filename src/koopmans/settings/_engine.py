from ._utils import Setting, SettingsDictWithChecks


class EngineSettingsDict(SettingsDictWithChecks):

    def __init__(self, **kwargs) -> None:
        settings = [
            Setting('npool',
                    'Number of pools for parallelizing over kpoints (should be commensurate with the k-point grid)',
                    int, None, None),
            Setting('from_scratch',
                    'if True, will delete any preexisting workflow and start again; '
                    'if False, will resume a workflow from where it was last up to',
                    bool, True, (True, False)),
            Setting('keep_tmpdirs',
                    'If False, delete all of the temporary directories at the end of the calculation',
                    bool, True, (True, False))]

        super().__init__(settings=settings, physicals=[], **kwargs)

    @property
    def _other_valid_keywords(self):
        return []
