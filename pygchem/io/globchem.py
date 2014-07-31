# -*- coding: utf-8 -*-

# module pygchem.io.globchem
# pygchem: Python interface for GEOS-Chem Chemistry Transport Model
#
# Copyright (C) 2012 Benoit Bovy
# see license.txt for more details
#
# Last modification: 03/2013

"""
Module for reading and writing GEOS-Chem simulation files related to the
global chemistry mechanism

- This module reads/writes GEOS-Chem simulation files "globchem.dat",
"jv_spec.dat", "jv_spec_aod.dat", "jv_atms.dat" and "ratj.d".
- This module also implements exporting of the global chemistry mechanism
to an sqlite database. It allows a deep exploration of the chemistry
mechanism trough complex queries.
- Functions in this module can be called directly or trough the classes
defined in the :module:`pygchem.glochem` module.

External dependencies: fortranformat, sqlite3
"""

import os
import collections

from pygchem.utils.fortranformat import FortranRecordReader
from pygchem.utils.fortranformat import FortranRecordWriter


def read_globchem_dat(filename):
    """
    Read and parse the "globchem.dat" globchem_file (SMVGEARII syntax).

    Parameters
    ----------
    filename : string
        Filename or path to the global chemistry globchem_file (globchem.dat)

    Returns
    -------
    header : string
        header (commented text) of the file
    mb_groups : dict
        mass balance groups (id, status): used for mass balance capability of
        chemical solvers (SMVGEARII or KPP)
    species : dict
        chemical species
    reac_kn : dict
        list of chemical kinectic reactions
    reac_ph : dict
        list of photolysis reactions
    species_links : list of tuples
        links between species made by the kinectic reactions

    """

    species = dict()
    reac_kn = dict()
    reac_ph = dict()
    species_links = list()

    # keys of 'species' and 'reac' dictionnaries
    species_keys = ['ord_n', 'fline', 'status', 'id', 'absorption',
                    'atomic_mass', 'init_concentration', 'mb_groups',
                    'formula', 'name']
    reac_keys = ['ord_n', 'fline', 'status', 'flag', 'comment',
                 'rate', 'reactants', 'products', 'coefs']

    def sanatize_str_inlist(l):
        """
        clean (strip...) string items in a list
        """
        return([x.strip() if type(x) is str else x for x in l])

    def go_2_begin(globchem_file, fline):
        """
        jump lines in the text file until reading the line
        with the word "BEGIN"
        """
        jumptext = ""
        while True:
            text = globchem_file.readline().strip()
            fline += 1
            if(text == "BEGIN"):
                break
            else:
                jumptext += text + '\n'
        return jumptext

    def convert_species_item(species_item, mb_groups):
        """
        Make the conversion (type and/or format), for each species, between
        the blocks of information encoded in the file globchem.dat and the
        dictionary items used in pygchem.
        """
        #species_item['absorption'] == bool(species_item['absorption'])
        species_item['absorption'] = False  # not used in GEOS-Chem
        species_item['atomic_mass'] = float(species_item['atomic_mass'])
        species_item['init_concentration'] = [float(i) for i in
                                            species_item['init_concentration']]
        mb_vals = [int(i) for i in species_item['mb_groups']]
        species_item['mb_groups'] = dict(zip(mb_groups.keys(), mb_vals))

    def trim_rlist(reac_list):
        """
        Take a 'reaction list' (e.g. reactants, products or stoichiometric
        coefficients) and return a list with removed empty items
        """
        # filter below permits to take only items that are not '' (if string)
        # or 0 (if int or float) because bool('') and bool(0) both return False
        return filter(bool, reac_list)

    def read_reac(globchem_file, fline):
        """
        Reads a 'chemical reaction' bloc of lines in
        globchem.dat

        Returns
        -------
        - reac_item: dict
            an item of the reaction dictionnary (see ...)
        """

        reac_r1_fmt = FortranRecordReader('(A1,1X,I4,1X,ES8.2,1X,ES8.1,1X,'
                                          'I6,1X,I1,1X,A2,F6.2,1X,2(F6.0,1X),'
                                          'A20)')
        reac_r2_fmt = FortranRecordReader('(7X,ES8.2,1X,ES8.1,1X,I6,1X,I1,3X,'
                                          'F6.2,1X,2(F6.0,1X))')
        reac_eq_fmt = FortranRecordReader('(4(A1,0PF5.3,A14)/'
                                          '4(A1,0PF5.3,A14)/'
                                          '4(A1,0PF5.3,A14)/'
                                          '4(A1,0PF5.3,A14)/'
                                          '4(A1,0PF5.3,A14))')

        # re-arrange reaction rate parameters:
        # parameter names for values (one-line) read in globchem.dat
        # and correspondence between structure in globchem.dat and structure
        # of pygchem.globchem module, depending on the reaction flag.
        # note that reac_rates_keys are modified for additional lines of
        # parameters (e.g., for line 1, 'AR' becomes 'AR1')
        # TODO: move this information to a separate function to make it 
        # available to export functions
        reac_rates_keys = ['AR', 'BR', 'CR', 'FCV', 'FCT1', 'FCT2']
        mod_ratep = {'' : {'arr': ['AR', 'BR', 'CR']},
                     'P': {'arrlow': ['AR', 'BR', 'CR'],
                           'arrhigh': ['AR1', 'BR1', 'CR1'],
                           'falloff': ['FCV', 'FCT1', 'FCT2']},
                     'E': {'arr': ['AR', 'BR', 'CR']},
                     'X': {'arrK0': ['AR', 'BR', 'CR'],
                           'arrK2': ['AR1', 'BR1', 'CR1'],
                           'arrK3': ['AR2', 'BR2', 'CR2']},
                     'Y': {'arrK0': ['AR', 'BR', 'CR']},
                     'Z': {'arrK1': ['AR', 'BR', 'CR'],
                           'arrK2': ['AR1', 'BR1', 'CR1']},
                     'A': {'arrK1': ['AR', 'BR', 'CR'],
                           'xcarbon': 'AR1'},
                     'B': {'arrK1': ['AR', 'BR', 'CR'],
                           'xcarbon': 'AR1'},
                     'C': {'arrK1': ['AR', 'BR', 'CR']},
                     'D': {'arrK1': ['AR', 'BR', 'CR']},
                     'V': {'arrK1': ['AR', 'BR', 'CR'],
                           'arrK2': ['AR1', 'BR1', 'CR1']},
                     'G': {'arrK1': ['AR', 'BR', 'CR'],
                           'arrK2': ['AR1', 'BR1', 'CR1']},
                    }

        # check if end of reaction list
        l_globchem = globchem_file.readline().strip()
        if(l_globchem[0:4] == "9999"):
            fline += 1
            return "end"

        # read reaction info
        reac_row1 = sanatize_str_inlist(reac_r1_fmt.read(l_globchem)) + [fline]
        reac_item = [reac_row1[i] for i in [1, -1, 0, 6, -2]]

        # read reaction rate parameters
        n_addrates = reac_row1.pop(5)
        reac_rates = dict(zip(reac_rates_keys,
                              reac_row1[2:5] + reac_row1[6:9]))
        if(n_addrates > 0):
            for i in range(0, n_addrates):
                new_keys = [s + str(i + 1) for s in reac_rates_keys]
                l_globchem = globchem_file.readline()
                fline += 1
                reac_row2 = reac_r2_fmt.read(l_globchem)
                del reac_row2[3]  # remove slot for number of rate coef lines
                reac_rates.update(dict(zip(new_keys, reac_row2)))

        # re-arrange rate parameters using corr_ratep
        reac_flag = reac_item[3].upper()
        mod_reac_rates = dict()
        if reac_flag in mod_ratep:
            for k, v in mod_ratep[reac_flag].items():
                if (isinstance(v, collections.Iterable)
                    and not isinstance(v, basestring)):
                    mod_reac_rates[k] = [reac_rates[n] for n in v]
                else:
                    mod_reac_rates[k] = reac_rates[v]
        else:
           mod_reac_rates.update(reac_rates)

        reac_item.append(mod_reac_rates)

        # read reaction equation
        l_globchem = ""
        for i in range(0, 5):
            l_globchem += globchem_file.readline()
            fline += 1
        reac_eq = sanatize_str_inlist(reac_eq_fmt.read(l_globchem))

        del reac_eq[0]
        if(reac_eq[-2] is None):
            reac_eq[-2] = 0.0

        reac_item.append(trim_rlist([reac_eq[i] for i in xrange(1, 13, 3)]))
        reac_item.append(trim_rlist([reac_eq[i] for i in xrange(13, 61, 3)]))
        reac_item.append(trim_rlist([reac_eq[i] for i in range(12, 60, 3)]))

        return reac_item

    # path to additional info files
    _dir_path = os.path.dirname(os.path.abspath(__file__))
    finfo_path = os.path.join(_dir_path, "..", "data")

    # read species.dat (get more info about species)
    species_info = {}
    s_fmt = FortranRecordReader('(A11,A32,A40)')
    f_species_dat = open(os.path.join(finfo_path, "species.dat"), "r")
    f_species_dat.readline()
    for line in f_species_dat:
        row = s_fmt.read(line)
        row = [s.strip() for s in row]
        species_info[row[0]] = row[1:]
    f_species_dat.close()
    species_info_blank = ['' for s in xrange(0, len(species_info.keys()[0]))]

    # globchem.dat file
    globchem_file = open(filename, "r")

    # get header/comments of globchem.dat
    fline = 1
    species_subheader = _format_subheader("Species List")
    header = go_2_begin(globchem_file, fline).replace(species_subheader, "")

    # read in species
    n_species = 0

    #    mass balance groups (id and status)
    species_mb_id_fmt = FortranRecordReader('20X,A3,8(1X,A3)')
    l_globchem = globchem_file.readline()
    mb_id = sanatize_str_inlist(species_mb_id_fmt.read(l_globchem))

    species_mb_st_fmt = FortranRecordReader('21X,A1,8(3X,A1)')
    l_globchem = globchem_file.readline()
    mb_st = sanatize_str_inlist(species_mb_st_fmt.read(l_globchem))

    mb_groups = dict(zip(mb_id, mb_st))
    fline += 2

    species_r1_fmt = FortranRecordReader('(A1,1X,A14,A2,1X,0PF6.2,4(1PE10.3))')
    while True:
        l_globchem = globchem_file.readline().strip()
        if(l_globchem == "END"):
            break

        species_r1 = sanatize_str_inlist(species_r1_fmt.read(l_globchem))
        species_this_key = species_r1[1]
        species_vals = ([n_species, fline] +
                        species_r1[0:4] +
                        [species_r1[4:]] +
                        [[int(i) for i in globchem_file.readline().split()]])
        try:
            species_vals += species_info[species_this_key]
        except KeyError:
            species_vals += species_info_blank

        species[species_this_key] = dict(zip(species_keys, species_vals))
        convert_species_item(species[species_this_key], mb_groups)
        del species_vals
        n_species += 1
        fline += 2

    # read in reactions
    go_2_begin(globchem_file, fline)
    mem_p_ratep = dict()
    n_reac = 0
    while True:
        reac_vals = read_reac(globchem_file, fline)
        if(reac_vals == "end"):
            break
        else:
            # case for reaction with E special flag (reverse equilibrium)
            # --> add rate params of previous reaction with P flag
            reac_flag = reac_vals[3]
            if reac_flag == 'P':
                mem_p_ratep = reac_vals[5]
            elif reac_flag == 'E':
                reac_vals[5].update(mem_p_ratep)

            reac_kn[n_reac] = dict(zip(reac_keys, reac_vals))
            reac_kn[n_reac]['id'] = n_reac
            n_reac += 1
        del reac_vals

    go_2_begin(globchem_file, fline)
    n_reac = 0
    while True:
        reac_vals = read_reac(globchem_file, fline)
        if(reac_vals == "end"):
            break
        else:
            reac_ph[n_reac] = dict(zip(reac_keys, reac_vals))
            reac_ph[n_reac]['id'] = n_reac
            n_reac += 1
        del reac_vals

    # identify links between species (through reactions): parse reac_kn
    link_id = 0
    for reac_key, reac_vals in reac_kn.iteritems():
        for reactant in reac_vals['reactants']:
            if(reactant != '' and reactant != 'M'):
                for product in reac_vals['products']:
                    if(product != ''):
                        species_links.append((link_id, reac_key,
                                              reactant, product))
                        link_id += 1

    # close globchem_file
    globchem_file.close()

    return header, mb_groups, species, reac_kn, reac_ph

def write_globchem_dat(filename, header, mb_groups, species, reac_kn, reac_ph,
                       overwrite=False):
    """
    Write global chemistry mechanism to a "globchem.dat" globchem_file
    (as used by GeosCHEM, following the smvgearII syntax)

    Parameters
    ----------
    filename : string
        filename or path to the global chemistry globchem_file (globchem.dat)
    header : string
        header (commented text) to append to the file
    mb_groups : dict
        mass balance groups (id, status): used for mass balance capability of
        chemical solvers (SMVGEARII or KPP)
    species : dict
        chemical species
    reac_kn : dict
        list of chemical kinectic reactions
    reac_ph : dict
        list of photolysis reactions
    overwrite : bool
        overwrite file if exists
    """

    # output globchem.dat file
    if os.path.exists(filename) and not overwrite:
        raise IOError("file %s already exists" % filename)
    globchem_file = open(filename, "w")

    globchem_file.write(header)
    globchem_file.write(_format_subheader("Species List"))

    globchem_file.write(_format_subheader("Kinetic reactions"))

    globchem_file.write(_format_subheader("Photolysis reactions"))


def _format_subheader(headtxt):
    """
    Return a formatted sub-header title for text files as follows
    
    #========================================================================
    # headtxt
    #======================================================================== 
    """
    subheader = "#=========================================================="\
                "===================\n# %s\n#==============================="\
                "==============================================" % headtxt
    return subheader
