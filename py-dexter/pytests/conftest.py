import numpy
import pytest
import dexter
import matplotlib


@pytest.fixture(autouse=True)
def add_imports(doctest_namespace):
    matplotlib.use("agg")  # Disable interactive plots
    doctest_namespace["dex"] = dexter
    doctest_namespace["np"] = numpy
