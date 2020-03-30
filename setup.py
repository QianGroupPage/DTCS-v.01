from setuptools import setup, find_packages

with open('README.md', 'r') as readme_file:
	readme=readme_file.read()

setup(name='crnsim',
	version='0.1',
	packages=['homogenous_crn'],
	package_dir={'': 'src'},
	requires=[
		'jupyterlab',
		'matplotlib',
		'numpy',
		'scipy'
		'sympy',
		'sklearn',],
	long_description=readme,
)