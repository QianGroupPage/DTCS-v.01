"""Setup package lblcrn.

This:
- Gets version information (__version__)
- Reads README.md into long_description
- Calls setuptools.setup()
"""

import glob
import re

import setuptools

# --- Collect version information --------------------------------------------
version_info = open('lblcrn/__init__.py').read()

pattern = r'__version__\s?=\s?[\'\"](.*)[\'"]'
match = re.search(pattern, version_info)
if match:
    __version__ = match.group(1)
else:
    __version__ = None

# --- Collect README ---------------------------------------------------------
with open('README.md', 'r') as readme_file:
    readme = readme_file.read()

# --- Get data files ---------------------------------------------------------
DATA_FP = 'lblcrn/'
DATA_GLOBS = [('examples/objects', ['examples/objects/*.json'])]

data_files = []
for path, patterns in DATA_GLOBS:
    files = []
    for pattern in patterns:
        files.extend([file.replace('\\', '/') for file in glob.glob(pattern)])
    data_files.append((DATA_FP + path, files))

# --- Setup ------------------------------------------------------------------
if __name__ == '__main__':
    setuptools.setup(name='lblcrn',
                     version=__version__,
                     packages=setuptools.find_packages(),
                     package_data={'lblcrn': ['resources/*.csv']},
                     data_files=data_files,
                     install_requires=[
                         'matplotlib',
                         'monty',
                         'numpy',
                         'pandas',
                         'pygame',
                         'scipy',
                         'seaborn',
                         'setuptools',
                         'sklearn',
                         'sympy'],
                     requires=[
                         'jupyterlab',
                     ],
                     long_description=readme,
                     )
