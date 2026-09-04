"""Common helper types to be used across all submodules."""

from typing import Any, TypeAlias

_PyAny: TypeAlias = Any
"""The exported type."""


class _ReprStrImpl:
    """Exports the wrapped type's `__str__` and `__repr__`."""

    _r: _PyAny

    def __repr__(self) -> str:
        return self._r.__repr__()

    def __str__(self) -> str:
        return self._r.__repr__()


class _RustTypeWrapper:
    """Wrappers around the exported Rust types."""

    _r: _PyAny

    @classmethod
    def _wrap(cls, _r: _PyAny) -> Any:
        """Creates a wrapped type from the corresponding exported type."""
        obj = cls.__new__(cls)
        obj._r = _r
        return obj
