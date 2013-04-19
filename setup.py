
from distutils.core import setup

setup(
    name='PyGChem',
    version='0.1.0',
    author='B.Bovy',
    author_email='bbovy@ulg.ac.be',
    packages=['pygchem',],
    scripts=[],
    url='http://pypi.python.org/pypi/PyGChem/',
    license='LICENSE.txt',
    description='Useful tools for GEOS-Chem simulations',
    long_description=open('README.txt').read(),
    install_requires=[
        "matplotlib >= 1.1.1",
        "numpy >= 1.5.0",
    ],
)
