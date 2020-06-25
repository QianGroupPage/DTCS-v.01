"""
setup.py
--------

For package lblcrn.
"""

import re
from setuptools import setup, find_packages

# Collection version information
version_info = open('lblcrn/__init__.py').read()

pattern = r'__version__\s?=\s?[\'\"](.*)[\'"]'
match = re.search(pattern, version_info)
if match:
    __version__ = match.group(1)

# Collect readme file
with open('README.md', 'r') as readme_file:
    readme = readme_file.read()

if __name__ == '__main__':
    # Setup
    setup(name='lblcrn',
          version=__version__,
          package_data={'lblcrn': ['resources/*.csv']},
          packages=find_packages(),
          install_requires=[
              'matplotlib',
              'monty',
              'numpy',
              'pandas',
              'scipy',
              'seaborn',
              'sklearn',
              'sympy'],
          requires=[
              'jupyterlab',
          ],
          long_description=readme,
          )
