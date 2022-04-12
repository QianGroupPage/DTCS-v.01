from typing import Any, List, DefaultDict, Dict, Optional
from src.model_input.model_input import InputElement, InputCollection


class ConditionElement:
    """A class used to store and manage input elements.

    Future work: Study the inheritance relationship from InputElement
    """
    def __init__(self,
                 name: str = "",
                 description: str = "",
                 value: Optional[Any] = None,
                 unit: str = ""
                 ):
        super().__init__(name, description, value)
        self.unit = unit


class ConditionCollection(InputCollection):
    """A class used to store and manage a collection of input elements.

    Terms:
        this collection: an object of this class;


    Attributes:
        elements (list): list of all input elements in the collection.
                         each element may or may not be of type InputElement or its
                         descendant class.
        elements_by_name (dictionary): a dictionary from input elements
    """
    def __init__(self,
                 name: str = "",
                 description: str = "",
                 elements: List = []):
        super().__init__(name=name, description=description, elements=elements)

