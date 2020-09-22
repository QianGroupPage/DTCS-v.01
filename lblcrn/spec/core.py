"""TODO

I had a lot of fun with this one.
- Andrew
"""

import abc
import collections
import copy
import warnings
from typing import Iterable, List

import monty.json  # TODO: Fix style?

import lblcrn


__author__ = 'Andrew Bogdan'
__email__ = 'andrewbogdan@lbl.gov'


class SpecMeta(abc.ABCMeta):

    _schema = set()
    _default = {}

    def __new__(mcs, name, bases, attrs):
        """Make Specs inherit their schema and default."""
        schema = set()
        default = {}

        # For all SpecMeta being based, inherit schema and default.
        for base in bases:
            if isinstance(base, mcs):
                schema.update(base._schema)
                default.update(base._default)

        # Update from attrs last to give the subbest class highest priority
        schema.update(attrs.get('_schema', set()))
        default.update(attrs.get('_default', {}))
        # Ensure that schema has all the keys in default.
        schema.update(default.keys())

        assert [hash(value) for value in default.values()] or True, \
            'Defaults should not be mutable.'

        attrs['_schema'] = schema
        attrs['_default'] = default

        return super().__new__(mcs, name, bases, attrs)


class Spec(collections.abc.MutableMapping,
           monty.json.MSONable,
           metaclass=SpecMeta):
    """TODO

    It's basically a wrapper around a dictionary (.spec) through which it
    keeps track of all its instance variables.
    """

    spec = None

    def __init__(self, **kwargs):
        self.spec = {}
        self.update(kwargs)

    def as_dict(self, sanitize=True) -> dict:
        """Return a MSON-serializable dict representation."""
        d = {
            '@module': self.__class__.__module__,
            '@class': self.__class__.__name__,
            '@version': lblcrn.__version__,  # TODO: Better way to do this?
        }
        d.update(self.spec)
        if sanitize:
            d = monty.json.jsanitize(d, strict=True)
        return d

    @classmethod
    def from_dict(cls, d: dict):
        """Load from a dict representation."""
        valid_spec = {}
        for key, value in d.items():
            # Add all the key/value pairs except monty reserved ones.
            if '@' in key:
                continue
            valid_spec[key] = value

        return cls(**valid_spec)

    def __copy__(self):
        return self.__class__(**self.spec)

    def __deepcopy__(self, memodict={}):
        spec_deepcopy = {}
        for key, value in self.spec.items():
            spec_deepcopy[key] = copy.deepcopy(value)
        return self.__class__(**spec_deepcopy)

    def __dir__(self) -> List[str]:
        attrs = list(super().__dir__())
        attrs.extend(key for key in self.spec if hasattr(self, key))
        return attrs

    def __len__(self) -> int:
        return len(self.spec)

    def __getitem__(self, key):
        """Forward to self.spec.

        Returns:
            It follows the search order:
            (1) If the name is in self.spec, return self.spec[name].
            Note: everything after this point is handled by __missing__
            (2) If the name is a key in self._default, give a default value:
                - If the default value is callable (e.g. list or dict), it
                    will call it and save it to spec.
                - If the default value isn't callable, it won't be saved to
                    the spec, but it's guaranteed to not be mutable.
            (3) If the name is in self._schema, return None.
            (4) Otherwise, raise an KeyError.
        """
        if key in self.spec:
            return self.spec[key]
        else:
            return self.__missing__(key)

    def __setitem__(self, key, value):
        """Forwards to the spec.

        Warns you if you try to set a key which isn't in self._schema.
        """
        if key not in self._schema:
            warnings.warn(f'\'{key}\' is not the schema for spec '
                          f'{self.__class__.__name__}.')
        self.spec[key] = value

    def __delitem__(self, key):
        """Forwards to the spec."""
        if key in self.spec:
            del self.spec[key]
        elif key in self._schema:
            # It's not in the spec, so there's nothing to delete.
            pass
        else:
            raise KeyError(f'\'{key}\'')

    def __iter__(self) -> Iterable:
        return iter(self.spec)

    def __contains__(self, key) -> bool:
        """Returns if the key is in the spec.

        It will not recognize default values, unless they are mutable ones
        and were already generated.
        """
        return key in self.spec

    def __missing__(self, key):
        """Uses self._schema and self.default to deal with missing keys.

        Returns:
            (1) If the name is a key in self._default, give a default value:
                - If the default value is callable (e.g. list or dict), it
                    will call it and save it to spec.
                - If the default value isn't callable, it won't be saved to
                    the spec, but it's guaranteed to not be mutable.
            (2) If the name is in self._schema, return None.
            (3) Otherwise, raise an KeyError.
        """
        if key in self._default:
            if callable(self._default[key]):
                self.spec[key] = self._default[key]()
                return self.spec[key]
            else:
                return self._default[key]
        elif key in self._schema:
            return None
        else:
            raise KeyError(f'\'{key}\'')

    def __getattr__(self, name):
        """Forwards to __getitem__."""
        try:
            return self[name]
        except IndexError:
            raise AttributeError(f'\'{self.__class__.__name__}\' object'
                                 f'has no attribute \'{name}\'')

    def __setattr__(self, name, value):
        """Forwards names not used by the class to __setitem__."""
        if hasattr(type(self), name):
            super().__setattr__(name, value)
        else:
            self[name] = value

    def __delattr__(self, name):
        """Delete the attributes in a way consistent with __setattr__."""
        if hasattr(type(self), name):
            super().__delattr__(name)
        else:
            try:
                del self[name]
            except KeyError:
                raise AttributeError(f'\'{self.__class__.__name__}\' object'
                                     f'has no attribute \'{name}\'')

    def __str__(self):
        if self.spec:
            spec_str = ''
            for key, value in self.spec.items():
                spec_str += f'{key}:\t{str(value)}\n'
            spec_str = spec_str[:-1]

            # Indent by one tab
            s = f'{self.__class__.__name__}:\n'
            for line in spec_str.split('\n'):
                s += '\t' + line + '\n'

            return s[:-1]
        else:
            return f'{self.__class__.__name__}()'

    def __repr__(self):
        if self.spec:
            s = f'{self.__class__.__name__}('
            for key, value in self.spec.items():
                s += f'{key}={repr(value)}, '
            s = s[:-2]
            s += ')'
            return s
        else:
            return f'{self.__class__.__name__}()'

