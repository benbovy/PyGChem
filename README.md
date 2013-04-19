PyGChem
=======

'PyGCHem' is a high-level, user-friendly Python interface to the 
GEOS-Chem Model, a global 3-D chemical transport model (CTM)
for atmospheric composition.

The project is currently at an early stage of development. 
This version is already usable, although many features still have to be
implemented. The first release will provide:
-   An object oriented library and a collection of routines (Python package)
    for model setup and parameterization (e.g., the global chemistry 
    mechanism), for I/O operations, and for processing/plotting outputs
    from the model ;
-   A collection of scripts, based on the library, to achieve common tasks
    related to GEOS-Chem simulations ;
-   Several GUI elements (Qt widgets) for common visualizations and
    animations ;
-   A customized GEOS-Chem Shell (using IPython or bash) to setup/run 
    simulations and perform post-processing/plotting tasks.


Dependencies
------------

Requirements:
-   Python 2.7 or higher (not 3.x)
-   Numpy 1.5 or higher
-   Matplotlib 1.0 or higher
-   Basemap extension for matplotlib

Optional modules:
-   IPython
-   PyQt4 4.4+
-   NetworkX

We recommend the installation of a scientific python distribution:
-   Python XY (Windows): <http://code.google.com/p/pythonxy/>
-   Enthought Canopy or EPD (Multi-platform, free academic license): 
    <https://www.enthought.com/>


Installation
------------

-   Download the code from github as zip archive [1][2], unzip and cd into
    'PyGChem-x-y-z'
-   Run 'python setup.py install' to install 'PyGChem'.
  

[1] <https://github.com/benbovy/PyGChem>
[2] You can also use the 'git' software with the following command:
    ```git clone https://github.com/benbovy/PyGChem.git```


Documentation
-------------

Full documentation with examples and illustrations is not available yet,
but Python classes and routines are well documented in the code (docstrings).


License
-------

Copyright (C) 2013 Benoît Bovy

Parts of the code are taken from the 'gchem' project 
<https://github.com/gkuhl/gchem/> Copyright (C) 2012 Gerrit Kuhlmann.

Licensed under the terms of the GPL License v3.
See LICENSE.txt for more information.


Contact
-------

Benoît Bovy <bbovy@ulg.ac.be>

Affiliation:
GIRPAS - Department of Astrophysics, Geophysics and Oceanography
University of Liège (Belgium)
<http://girpas.astro.ulg.ac.be/>


