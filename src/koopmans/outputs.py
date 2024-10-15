from pydantic import BaseModel


class OutputModel(BaseModel):

    def todict(self):
        dct = self.model_dump()
        dct['__koopmans_name__'] = self.__class__.__name__
        dct['__koopmans_module__'] = self.__module__
        return dct

    @classmethod
    def fromdict(cls, dct):
        return cls(**dct)
