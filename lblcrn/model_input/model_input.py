from typing import Any, List, DefaultDict, Dict, Optional
from collections.abc import Iterable
from collections import defaultdict
import monty.json


class InputElement(monty.json.MSONable):
    """A class used to store and manage input elements.

    """
    def __init__(self,
                 name: str = "",
                 description: str = "",
                 value: Optional[Any] = None,
                 ):
        self.description = description
        self.name = name
        self.value = value


# isobaric_group
#
# isothermal_group
#
#
# # A file to Keep uploading rules
# dose 1 torr of water gas to silver 111 surface, keep temperature at Room Temperature.
# action:
# data file that save the rules
# data file, expanding
#
# data.rules.hydration
# data.rules.co2_reduction
# data.rules.haber_bosch


class InputCollection(InputElement):
    """A class used to store and manage a collection of input elements.

    Terms:
        this collection: an object of this class;

    Attributes:
        elements (list): list of all input elements in the collection.
                         each element may or may not be of type InputElement or its
                         descendant class.
        elements_by_name (dictionary): a dictionary from input elements.
    """
    def __init__(self,
                 name: str = "",
                 description: str = "",
                 elements: List = []):
        super().__init__(name, description)
        self.elements = elements
        self.elements_by_type = self.make_elements_by_type()
        self.elements_by_name = self.make_elements_by_name()
        self.values_by_name = self.make_values_by_name()

    def make_elements_by_type(self) -> DefaultDict:
        """Make a dictionary from type of an element to a list of elements with the type.

        Returns:
             a DefaultDict which defaults to an empty list, element type -> list of elements.
        """
        res = defaultdict(list)
        for ele in self.elements:
            res[type(ele)].append(ele)
        return res

    def make_elements_by_name(self) -> Dict:
        """Make a dictionary from the name of an element to the element.

        Returns:
            a dictionary mapping of element name -> element.
        """
        return {ele.name: ele for ele in self.elements}

    def make_values_by_name(self) -> Dict:
        """Make a dictionary from the name of an element to the element.

        Returns:
            a dictionary mapping of element name -> element value.
        """
        return {ele.name: ele.value for ele in self.elements}

    def _append_and_update(self, element: Any):
        """Internal method to append an element to the collection and update it for each
        dictionary where the element will appear.
        """
        self.elements.append(element)
        self.elements_by_type = self.make_elements_by_type()
        self.elements_by_name = self.make_elements_by_name()
        self.values_by_name = self.make_values_by_name()
        self._update_elements(self.elements)

    def _update_elements(self, elements):
        """
        Prototypical method to update any dictionary formats of an element.

        :param elements: a list of elements;
        :return: None
        """
        pass

    def append(self, element: Any):
        """Append an element to this collection.
        """
        if element.name in self.elements_by_name:
            print(f"Failed to add {element} to InputCollection {self}\n" +
                  f"{self.elements_by_name[element.name]} is already in the InputCollection with name" +
                  f" {element.name}.")
            return
        self._append_and_update(element)

    def extend(self, obj: Iterable):
        """Extend a collection of elements to this collection.
        """
        if isinstance(obj, Iterable):
            for ele in obj:
                self.append(ele)
        else:
            print(f"{obj} is of type {type(obj)}. {self}.extend() requires an Iterable input.")

    def __add__(self, ele: InputElement):
        """Add an input element to this collection.

        Supports self += input element
        """
        self.append(ele)

    def __iter__(self):
        return iter(self.elements)

    def __len__(self):
        return len(self)
