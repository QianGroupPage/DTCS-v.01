"""Test data file access.

It's easier to debug resources like this.
"""

import unittest

from dtcs import *


class TestExamples(unittest.TestCase):
    """Tests examples."""

    def test_examples(self):
        for example in dtcs_examples.list_all():
            dtcs_examples.load(example)


if __name__ == '__main__':
    unittest.main()
