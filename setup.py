"""Setup package lblcrn.

This:
- Gets version information (__version__)
- Reads README.md into long_description
- Calls setuptools.setup()
"""

import collections
import fnmatch
import os
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
# The directory under venv/ to store datafiles.
DATA_FP = 'lblcrn/'

# The files to include. These are triples of the form
# (destination path, path to search, [file patterns])
# The search is recursive.
INCLUDE = [
    ('examples/objects', 'examples/objects', ['*.json']),
    ('docs', 'docs/build', ['*']),
]

mapping = collections.defaultdict(list)
for new_root, search_root, patterns in INCLUDE:
    for dir_path, subdirs, files in os.walk(search_root):
        clean_path = dir_path.replace(os.sep, '/')
        new_path = DATA_FP + clean_path.replace(search_root, new_root)
        for file in files:
            if any([fnmatch.fnmatch(file, pattern) for pattern in patterns]):
                mapping[new_path].append(clean_path + '/' + file)

data_files = [(dir_name, files) for dir_name, files in mapping.items()]
#import pprint; pprint.pprint(data_files); input('>')

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
                         'sympy',
                         'ffmpeg',
                         'signal'],
                     requires=[
                         'jupyterlab',
                     ],
                     long_description=readme,
                     )
