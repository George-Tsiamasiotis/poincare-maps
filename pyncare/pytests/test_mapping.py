import pytest
from pyncare import Mapping


def test_mapping_fields(mapping: Mapping):
    """Tests that the fields derived from the wrapped Rust object are the correct type."""
    assert isinstance(mapping.section, str)
    assert isinstance(mapping.alpha, float)
    assert isinstance(mapping.intersections, int)


def test_immutability(mapping: Mapping):
    """Tests that current fields are immutable."""
    with pytest.raises(AttributeError):
        mapping.section = ""
        mapping.alpha += 1
        mapping.intersections += 1


def test_repr(mapping: Mapping):
    str(mapping)
