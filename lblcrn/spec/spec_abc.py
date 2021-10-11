"""TODO(Andrew)

I'm so sorry for using a metaclass
- Andrew Bogdan

Credit: Ye Wang, Andrew Bogdan
"""

from typing import Any, List, DefaultDict, Dict, Optional, Iterable, Union

import collections
import copy
import warnings

import monty.json

import lblcrn


class _SpecView(collections.abc.MutableMapping):
    def __init__(self, dic):
        self._dict = dic

    def __getitem__(self, item):
        if not isinstance(item, str):
            raise TypeError('Specs only accept string keys.')
        if item.startswith('_'):
            raise KeyError(f'{item}')
        return self._dict.__getitem__(item)

    def __setitem__(self, key, value):
        if not isinstance(key, str):
            raise TypeError('Specs only accept string keys.')
        if key.startswith('_'):
            raise ValueError(f'Cannot set protected keys.')
        return self._dict.__setitem__(key, value)

    def __delitem__(self, item):
        if not isinstance(item, str):
            raise TypeError('Specs only accept string keys.')
        if item.startswith('_'):
            raise ValueError(f'Cannot delete protected keys.')
        return self._dict.__delitem__(item)

    def __iter__(self):
        for key in self._dict.keys():
            if not key.startswith('_'):
                yield key

    def __len__(self):
        return len(list(iter(self)))

    def __copy__(self):
        return _SpecView(dict(self))

    def __deepcopy__(self):
        return _SpecView(copy.deepcopy(dict(self)))

    def __str__(self):
        return str(dict(self))

    def __repr__(self):
        return repr(dict(self))


class SpecABC(monty.json.MSONable):
    """TODO

    TODO(Andrew): I should make it hard to tell that it's special. Devs and
     users shouldn't be able to tell it exists, just that serializing is
     super easy. Don't make it act like a dict, make a sublcass act like a
     dict. This should just be object++.
    """

    def __init__(self, **kwargs):
        # TODO(Andrew) Check and throw warnings if I get unknown kwargs
        self.spec.update(kwargs)

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
        decode = monty.json.MontyDecoder().process_decoded
        spec = {}
        for key, value in d.items():
            # Add all the key/value pairs except monty reserved ones.
            if '@' in key:
                continue
            spec[key] = decode(value)

        return cls(**spec)

    @property
    def spec(self):
        return _SpecView(self.__dict__)

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

class Spec(SpecABC):
    """TODO(Andrew)"""
    def __init__(self,
                 name: str = '',
                 description: str = '',
                 **kwargs):
        self.name = name
        self.description = description
        super().__init__(**kwargs)


class SpecCollection(Spec,
                     collections.abc.MutableSequence):
    """A class used to store and manage a collection of input elements.

    Attributes:
        elements: list of Specs in the collection.
            Each element should be a subclass of Spec.

    Methods:
        TODO(list methods, incl. ones inherited from MutableSequence
    """
    def __init__(self,
                 *elements: List[Spec],
                 name: str = '',
                 description: str = '',):
        super().__init__(name=name, description=description)

        self.elements = []
        self._names = {}
        elements = elements or []
        for ele in elements:
            if not isinstance(ele, Spec):
                raise TypeError('Elements must be a subclass of Spec')
            if ele.name in self._names:
                raise ValueError('Element names must be unique.')
            self._names[ele.name] = len(self.elements)
            self.elements.append(ele)

    @classmethod
    def from_dict(cls, d: dict):
        """Load from a dict representation."""
        decode = monty.json.MontyDecoder().process_decoded
        spec = {}
        for key, value in d.items():
            # Add all the key/value pairs except monty reserved ones.
            if '@' in key:
                continue
            spec[key] = decode(value)

        return cls(*spec.pop('elements'), **spec)

    # --- Helpful Accessor Methods --------------------------------------------
    # TODO(Andrew): Ye originally made this where it was cached, so that you
    #  could access it a bunch quickly. I'll re-create this later, but I
    #  removed it for ease of development. I commented out the cache lines.
    #  You'd need to make a _cache variable that didn't auto-serialize.

    def by_type(self) -> DefaultDict:
        """Make a dictionary from type of an element to a list of elements with the type.

        Returns:
             a DefaultDict which defaults to an empty list, element type -> list of elements.
        """
        res = collections.defaultdict(list)
        for ele in self.elements:
            res[type(ele)].append(ele)
        return res

    def by_subclass(self) -> DefaultDict:
        """Make a dictionary from class of an element to a list of elements
        with that class as a superclass. TODO(Andrew)

        Returns:
             a DefaultDict which defaults to an empty list, element class ->
                 list of elements.
        """
        res = collections.defaultdict(list)
        for ele in self.elements:
            for cls in ele.__class__.__mro__:
                res[cls].append(ele)
        return res

    def by_name(self) -> Dict:
        """Make a dictionary from the name of an element to the element.

        Returns:
            a dictionary mapping of element name -> element.
        """
        return {ele.name: ele for ele in self.elements}

    # --- Sequence Implementation ---------------------------------------------
    def insert(self, index, value: Spec):
        # TODO(Andrew) only supports integer indexing
        if not isinstance(value, Spec):
            raise TypeError('Elements must be a subclass of Spec.')
        if value.name in self._names:
            raise ValueError('Element names must be unique.')
        self.elements.insert(index, value)
        self._names = {self.elements[index].name: index for index in
                       range(len(self.elements))}

    def __getitem__(self, item: Union[int, str]):
        # TODO(Andrew) Supports integer and name-indexing
        if isinstance(item, int):
            return self.elements.__getitem__(item)
        elif isinstance(item, str):
            return self.elements.__getitem__(self._names[item])
        else:
            raise TypeError('Index must be int or str.')

    def __setitem__(self, key: Union[int, str], value: Spec):
        # TODO(Andrew) Supports integer and name-indexing
        if not isinstance(value, Spec):
            raise TypeError('Elements must be a subclass of Spec.')
        if isinstance(key, int):
            return self.elements.__setitem__(key, value)
        elif isinstance(key, str):
            if key != value.name:
                raise ValueError('Key must be equal to spec.name')
            self._names[key] = len(self.elements)
            return self.elements.__setitem__(self._names[key], value)
        else:
            raise TypeError('Index must be int or str.')

    def __delitem__(self, item: Union[int, str]):
        # TODO(Andrew) Supports integer and name-indexing
        if isinstance(item, int):
            return self.elements.__delitem__(item)
        elif isinstance(item, str):
            return self.elements.__delitem__(self._names[item])
        else:
            raise TypeError('Index must be int or str.')

    def __len__(self):
        return len(self.elements)

    def __iter__(self):
        return iter(self.elements)

    # --- Additional Magic ----------------------------------------------------
    def __contains__(self, item: str):
        if not isinstance(item, str):
            raise TypeError('Checking containment accepts names only.')
        return item in self._names

    # TODO(Andrew) implement iteritems ?

    # @property
    # def values_by_name(self) -> Dict:
    #     """Make a dictionary from the name of an element to the element.
    #
    #     Returns:
    #         a dictionary mapping of element name -> element value.
    #     """
    #     return {ele.name: ele.value for ele in self.elements}