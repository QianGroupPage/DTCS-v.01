from setuptools import setup, find_packages

with open('README.md', 'r') as readme_file:
	readme=readme_file.read()

setup(name='crnsim',
	version='0.1',
	package_data={'lblcrn': ['resources/*.csv']},
	packages=find_packages(),
	requires=[
		'jupyterlab',
		'matplotlib',
		'numpy',
		'scipy'
		'sympy',
		'sklearn',],
	long_description=readme,
)