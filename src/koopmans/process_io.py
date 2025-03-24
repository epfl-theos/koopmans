import numpy as np
from pydantic import BaseModel


class IOModel(BaseModel):
    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        for key in self.dict().keys():
            cond = (getattr(self, key) == getattr(other, key))
            if isinstance(cond, np.ndarray):
                cond = cond.all()
            if not cond:
                return False
        return True

    def todict(self):
        dct = self.model_dump()
        dct['__koopmans_name__'] = self.__class__.__name__
        dct['__koopmans_module__'] = self.__module__
        return dct

    @classmethod
    def fromdict(cls, dct):
        return cls(**dct)
