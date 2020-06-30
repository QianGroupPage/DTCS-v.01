"""Test data file access.

It's easier to debug resources like this.
"""

import unittest

from lblcrn import *


class TestExamples(unittest.TestCase):
    """Tests examples."""

    def test_examples(self):
        for example in lblcrn_examples.list_all():
            lblcrn_examples.load(example)


if __name__ == '__main__':
    unittest.main()
