from typing import TextIO

from koopmans.ml import AbstractMLModel
from koopmans.utils import serialization


def read_ml_model(fd: TextIO):
    return serialization.decode(fd.read())


def write_ml_model(obj: AbstractMLModel, fd: TextIO):
    fd.write(serialization.encode(obj))
