"""Test resource access.

It's easier to debug resources like this.
"""

import unittest

from lblcrn import *


class TestResources(unittest.TestCase):
    """Tests resource access."""

    def test_examples(self):
        cts = lblcrn_examples.load('nickel-cts')
        self.assertEqual(type(cts), CRNTimeSeries)


if __name__ == '__main__':
    unittest.main()
