
import os
from distutils.core import setup

from pygchem import __version__, __project_url__

NAME = 'PyGChem'
LIBNAME = 'pygchem'

def get_package_data(name, extlist):
    """Return data files for package *name* with extensions in *extlist*"""
    flist = []
    # Workaround to replace os.path.relpath (not available until Python 2.6):
    offset = len(name)+len(os.pathsep)
    for dirpath, _dirnames, filenames in os.walk(name):
        for fname in filenames:
            if (not fname.startswith('.')
                and os.path.splitext(fname)[1] in extlist):
                flist.append(os.path.join(dirpath, fname)[offset:])
    print flist
    return flist

def get_subpackages(name):
    """Return subpackages of package *name*"""
    splist = []
    for dirpath, _dirnames, _filenames in os.walk(name):
        if os.path.isfile(os.path.join(dirpath, '__init__.py')):
            splist.append(".".join(dirpath.split(os.sep)))
    print splist
    return splist

def get_packages():
    """Return package list"""
    packages = get_subpackages(LIBNAME)  # + get_subpackages('pygchemplugins')
    return packages

setup(name=NAME,
      version=__version__,
      description='Python interface to GEOS-Chem 3D Chemistry Transport Model',
      long_description=open('README.md').read(),
      keywords='Interface GEOS-Chem Chemistry Transport Model CTM Atmospheric',
      author='Benoit Bovy',
      author_email='bbovy@ulg.ac.be',
      license='GPLv3',
      url=__project_url__,
      download_url='%s/archive/%s.tar.gz' % (__project_url__, __version__),
      platforms=['any'],
      packages=get_packages(),
      package_data={LIBNAME: get_package_data(LIBNAME, ('.dat'))},
      scripts=[],
      requires=["matplotlib (>= 1.1.0)",
                "numpy (>= 1.5.0)"])
