"""Common helper types to be used across all submodules."""

from typing import Any


class _ReprStrImpl:
    """To be inherited by any wrapper types to avoid repetition."""

    _r: Any

    def __repr__(self) -> str:
        return self._r.__repr__()

    def __str__(self) -> str:
        return self._r.__repr__()
