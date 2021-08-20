import itertools
import numpy as np
import pandas as pd


class Band(object):
    def __init__(self, index=None, filled=True, group=None, alpha=None, error=None, dct={}):
        self.index = index
        self.filled = filled
        self.group = group
        self.alpha_history = []
        self.error_history = []
        self.alpha = alpha
        self.error = error
        if dct:
            self.fromdct(dct)

    def fromdct(self, dct):
        for k, v in dct.items():
            assert hasattr(self, k)
            if v is not None:
                setattr(self, k, v)

    def todict(self):
        dct = self.__dict__
        dct['__koopmans_name__'] = self.__class__.__name__
        dct['__koopmans_module__'] = self.__class__.__module__
        return dct

    @property
    def alpha(self):
        assert len(self.alpha_history) > 0, 'Band does not have screening parameters'
        return self.alpha_history[-1]

    @alpha.setter
    def alpha(self, value):
        if value:
            self.alpha_history.append(value)

    @property
    def error(self):
        assert len(self.error_history) > 0, 'Band does not have error data'
        return self.error_history[-1]

    @error.setter
    def error(self, value):
        if value:
            self.error_history.append(value)


class Bands(object):
    def __init__(self, n_bands=None, bands=None, dct={}, **kwargs):
        if bands is None and n_bands:
            self._bands = [Band(i + 1) for i in range(n_bands)]
        elif bands and n_bands is None:
            self._bands = bands
        elif not dct:
            raise ValueError('The arguments "n_bands" and "bands" are mutually exclusive')
        for k, v in kwargs.items():
            assert hasattr(self, k)
            if v:
                setattr(self, k, v)
        if dct:
            self.fromdct(dct)

    def __iter__(self):
        for b in self._bands:
            yield b

    def fromdct(self, dct):
        for k, v in dct.items():
            if v is not None:
                setattr(self, k, v)

    def todict(self):
        dct = self.__dict__
        dct['__koopmans_name__'] = self.__class__.__name__
        dct['__koopmans_module__'] = self.__class__.__module__
        return dct

    def get(self, filled=None, group=None, to_solve=None):
        if to_solve is None:
            band_subset = self._bands
        elif to_solve:
            band_subset = self.bands_to_solve
        return Bands(bands=[b for b in band_subset if (getattr(b, 'filled') == filled or filled is None)
                            and (getattr(b, 'group') == group or group is None)])

    def __getitem__(self, key):
        return self._bands[key]

    def num(self, filled=None):
        if filled is None:
            return len(self._bands)
        elif filled is True:
            return len([b for b in self._bands if b.filled])
        elif filled is False:
            return len([b for b in self._bands if not b.filled])
        else:
            raise ValueError(f'Invalid choice "filled" = {filled}')

    @property
    def filling(self):
        return [b.filled for b in self._bands]

    @filling.setter
    def filling(self, value):
        for b, v in zip(self._bands, value):
            b.filled = v

    @property
    def indices(self):
        return [b.index for b in self._bands]

    @property
    def groups(self):
        return [b.group for b in self._bands]

    @groups.setter
    def groups(self, value):
        assert len(value) == len(
            self._bands), f'You tried to set the orbital groups with a list of length {len(value)} != {self.num()}'
        for i, v in enumerate(value):
            self._bands[i].group = v

    @property
    def to_solve(self):
        # Update which bands to solve explicitly

        # If groups have not been assigned, solve all bands
        if None in self.groups:
            return self._bands

        # If not, work out which bands to solve explicitly
        groups_found = set([])
        to_solve = []

        for band in [b for b in self._bands[::-1] if b.filled] + [b for b in self._bands if not b.filled]:
            # Looping through the filled bands from highest to lowest, then empty bands from
            # lowest to highest
            if band.group not in groups_found:
                groups_found.add(band.group)
                to_solve.append(band)

        if groups_found != set(self.groups):
            raise ValueError('Splitting of orbitals into groups failed')

        return sorted(to_solve, key=lambda x: x.index)

    @property
    def alphas(self):
        # This returns the alpha values for the iteration number where we have alpha for all bands
        i = min([len(b.alpha_history) for b in self._bands]) - 1
        if i == -1:
            raise AttributeError()
        return [b.alpha_history[i] for b in self._bands]

    @alphas.setter
    def alphas(self, value):
        self.update_alphas(value)

    def update_alphas(self, value, group=None):
        '''
        Sets the band's screening parameters to the value provided
         - "value" can be a scalar, a list, or a pandas DataFrame of the alpha_history
         - if "group" is provided then it applies this value to the orbitals belonging to this group only
        '''

        if isinstance(value, pd.DataFrame):
            assert group is None, 'Cannot update only one group via a pandas DataFrame'
            for b, alpha_history in zip(self._bands, np.transpose(value.values.tolist())):
                # Make sure to exclude NaNs
                b.alpha_history = [a for a in alpha_history.tolist() if not np.isnan(a)]
            return

        if isinstance(value, float):
            value = [value for _ in range(self.num())]
        assert len(value) == len(
            self._bands), f'You tried to set the orbital alphas with a list of length {len(value)} != {self.num()}'
        for i, v in enumerate(value):
            if group:
                if self._bands[i].group != group:
                    continue
            self._bands[i].alpha = v

    @property
    def errors(self):
        return [b.error for b in self._bands]

    @errors.setter
    def errors(self, value):
        self.update_errors(value)

    def update_errors(self, value, group=None):
        '''
        Sets the band's residual error to the value provided
         - "value" can be a scalar, a list, or a pandas DataFrame
         - if "group" is provided then it applies this value to the orbitals belonging to this group only
        '''

        if isinstance(value, pd.DataFrame):
            assert group is None, 'Cannot update only one group via a pandas DataFrame'
            if not value.empty:
                for b, error_history in zip(self._bands, np.transpose(value.values.tolist())):
                    b.error_history = [e for e in error_history.tolist() if not np.isnan(e)]
            return

        if isinstance(value, float):
            value = [value for _ in range(self.num())]
        assert len(value) == len(
            self._bands), f'You tried to set the orbital errors with a list of length {len(value)} != {self.num()}'
        for i, v in enumerate(value):
            if group:
                if self._bands[i].group != group:
                    continue
            self._bands[i].error = v

    @property
    def alpha_history(self):
        # Create an array of values padded with NaNs
        data = np.array(list(itertools.zip_longest(*[b.alpha_history for b in self._bands], fillvalue=np.nan)))
        return pd.DataFrame(data, columns=self.indices)

    @property
    def error_history(self):
        # Create an array of values padded with NaNs
        data = np.array(list(itertools.zip_longest(*[b.error_history for b in self._bands], fillvalue=np.nan)))
        if data.size == 0:
            data = None
        return pd.DataFrame(data, columns=self.indices)

    def print_history(self):
        # Printing out a progress summary
        print('\nalpha')
        print(self.alpha_history)
        if not self.error_history.empty:
            print('\nDelta E_i - epsilon_i (eV)')
            print(self.error_history)
        print('')
