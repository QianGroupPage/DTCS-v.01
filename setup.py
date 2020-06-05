from setuptools import setup, find_packages

with open('README.md', 'r') as readme_file:
	readme=readme_file.read()

setup(name='lblcrn',
	version='0.1',
	package_data={'lblcrn': ['resources/*.csv']},
	packages=find_packages(),
	install_requires=[
		'jupyterlab',
		'matplotlib',
		'numpy',
		'pandas',
		'scipy',
		'sympy',
		'sklearn',],
	long_description=readme,
)