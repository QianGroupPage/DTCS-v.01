"""
Color Map: a way for all species to have the same color within runs of the
the program. Once you restart, it'll re-pick the colors, but that's fine.
"""
from typing import Union
import collections.abc
import itertools

import matplotlib.colors
import sympy as sym

from dtcs.common import util

_COLORS = itertools.cycle(['red', 'green', 'orange', 'blue', 'purple', 'pink',
                           'yellow', 'gray', 'cyan'])

class _ColorManagerSingleton(collections.abc.Mapping):

    singleton = None

    def __new__(cls, *args, **kwargs):
        """Make this class a singleton, so only one exists."""
        if not cls.singleton:
            cls.singleton = super().__new__(cls, *args, **kwargs)
            return cls.singleton
        else:
            return cls.singleton

    def __init__(self):
        super().__init__()
        self._map = {}

    def hex(self, name: Union[str, sym.Symbol]) -> str:
        """Return the color as hex (#RRGGBB)"""
        return matplotlib.colors.to_hex(self.__getitem__(name))

    def rgb(self, name: Union[str, sym.Symbol]) -> tuple:
        """Return the color as rgb, a tuple of (R, G, B) in [0, 1]"""
        return matplotlib.colors.to_rgb(self.__getitem__(name))

    def rgb256(self, name: Union[str, sym.Symbol]) -> tuple:
        """Return the color as rgb, a tuple of (R, G, B) in [0, 255]"""
        return tuple(round(value * 255) for value in self.rgb(name))

    def rgba(self, name: Union[str, sym.Symbol]) -> tuple:
        """Return the color as rgba, a tuple of (R, G, B, A) in [0, 1]"""
        return matplotlib.colors.to_rgba(self.__getitem__(name))

    def rgba256(self, name: Union[str, sym.Symbol]) -> tuple:
        """Return the color as rgba, a tuple of (R, G, B, A) in [0, 255]"""
        return tuple(round(value * 255) for value in self.rgba(name))

    def _get_color(self, name: Union[str, sym.Symbol]) -> str:
        """
        Get the color corresponding to the name from _map. If there is none, it
        will set one.
        """
        name = util.symbol_to_name(name)
        if name in self._map:
            return self._map[name]
        else:
            return self._default(name)

    def _set_color(self, name: Union[str, sym.Symbol], color: str):
        """
        Set the color corresponding to the given name.
         - This will attempt to convert all colors into hex.
         - It will also pick a default if you provide it with nothing.
        """
        name = util.symbol_to_name(name)
        if not color:
            self._default(name)
            return
        color = matplotlib.colors.to_hex(color, keep_alpha=False)
        self._map[name] = color

    def _default(self, name: Union[str, sym.Symbol]) -> str:
        """Pick a color for that name and save it."""
        # TODO(Andrew) Add a warning/info/debug here.
        self._set_color(name, next(_COLORS))
        return self._get_color(name)

    def __getitem__(self, key: Union[str, sym.Symbol]) -> str:
        """Get the color globally corresponding to the given name."""
        return self._get_color(key)

    def __setitem__(self, key: Union[str, sym.Symbol], value: str):
        """Set the color globally corresponding to the given name."""
        return self._set_color(key, value)

    def __len__(self):
        return len(self._map)

    def __iter__(self):
        return iter(self._map)

    def __repr__(self):
        return dict.__repr__(self._map)

    def __str__(self):
        return dict.__str__(self._map)


color_map = _ColorManagerSingleton()
