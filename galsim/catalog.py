# Copyright (c) 2012-2017 by the GalSim developers team on GitHub
# https://github.com/GalSim-developers
#
# This file is part of GalSim: The modular galaxy image simulation toolkit.
# https://github.com/GalSim-developers/GalSim
#
# GalSim is free software: redistribution and use in source and binary forms,
# with or without modification, are permitted provided that the following
# conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions, and the disclaimer given in the accompanying LICENSE
#    file.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions, and the disclaimer given in the documentation
#    and/or other materials provided with the distribution.
#
"""@file catalog.py
Routines for controlling catalog input/output with GalSim.
"""

from six import iteritems, iterkeys, itervalues
from builtins import zip
import galsim
import os
import numpy as np

class Catalog(object):
    """A class storing the data from an input catalog.

    Each row corresponds to a different object to be built, and each column stores some item of
    information about that object (e.g. flux or half_light_radius).

    Initialization
    --------------

    @param file_name    Filename of the input catalog. (Required)
    @param dir          Optionally a directory name can be provided if `file_name` does not
                        already include it.
    @param file_type    Either 'ASCII' or 'FITS'.  If None, infer from `file_name` ending.
                        [default: None]
    @param comments     The character used to indicate the start of a comment in an
                        ASCII catalog.  [default: '#']
    @param hdu          Which hdu to use for FITS files.  [default: 1]

    Attributes
    ----------

    After construction, the following attributes are available:

        nobjects   The number of objects in the catalog.
        ncols      The number of columns in the catalog.
        isfits     Whether the catalog is a fits catalog.
        names      For a fits catalog, the valid column names.

    """
    _req_params = { 'file_name' : str }
    _opt_params = { 'dir' : str , 'file_type' : str , 'comments' : str , 'hdu' : int }
    _single_params = []
    _takes_rng = False

    # _nobjects_only is an intentionally undocumented kwarg that should be used only by
    # the config structure.  It indicates that all we care about is the nobjects parameter.
    # So skip any other calculations that might normally be necessary on construction.
    def __init__(self, file_name, dir=None, file_type=None, comments='#', hdu=1,
                 _nobjects_only=False):

        # First build full file_name
        self.file_name = file_name.strip()
        if dir is not None:
            import os
            self.file_name = os.path.join(dir,self.file_name)

        if file_type is None:
            import os
            name, ext = os.path.splitext(file_name)
            if ext.lower().startswith('.fit'):
                file_type = 'FITS'
            else:
                file_type = 'ASCII'
        file_type = file_type.upper()
        if file_type not in ['FITS', 'ASCII']:
            raise ValueError("file_type must be either FITS or ASCII if specified.")
        self.file_type = file_type
        self.comments = comments
        self.hdu = hdu

        if file_type == 'FITS':
            self.readFits(hdu, _nobjects_only)
        elif file_type == 'ASCII':
            self.readAscii(comments, _nobjects_only)
        else:
            raise ValueError("Invalid file_type %s"%file_type)

    # When we make a proxy of this class (cf. galsim/config/stamp.py), the attributes
    # don't get proxied.  Only callable methods are.  So make method versions of these.
    def getNObjects(self) : return self.nobjects
    def isFits(self) : return self.isfits

    def readAscii(self, comments, _nobjects_only=False):
        """Read in an input catalog from an ASCII file.
        """
        # If all we care about is nobjects, this is quicker:
        if _nobjects_only:
            # See the script devel/testlinecounting.py that tests several possibilities.
            # An even faster version using buffering is possible although it requires some care
            # around edge cases, so we use this one instead, which is "correct by inspection".
            with open(self.file_name) as f:
                if (len(comments) == 1):
                    c = comments[0]
                    self.nobjects = sum(1 for line in f if line[0] != c)
                else:
                    self.nobjects = sum(1 for line in f if not line.startswith(comments))
            return

        # Read in the data using the numpy convenience function
        # Note: we leave the data as str, rather than convert to float, so that if
        # we have any str fields, they don't give an error here.  They'll only give an
        # error if one tries to convert them to float at some point.
        self.data = np.loadtxt(self.file_name, comments=comments, dtype=bytes)
        # Convert the bytes to str.  For Py2, this is a no op.
        self.data = self.data.astype(str)

        # If only one row, then the shape comes in as one-d.
        if len(self.data.shape) == 1:
            self.data = self.data.reshape(1, -1)
        if len(self.data.shape) != 2:
            raise IOError('Unable to parse the input catalog as a 2-d array')

        self.nobjects = self.data.shape[0]
        self.ncols = self.data.shape[1]
        self.isfits = False

    def readFits(self, hdu, _nobjects_only=False):
        """Read in an input catalog from a FITS file.
        """
        from galsim._pyfits import pyfits, pyfits_version
        with pyfits.open(self.file_name) as fits:
            raw_data = fits[hdu].data
        if pyfits_version > '3.0':
            self.names = raw_data.columns.names
        else:
            self.names = raw_data.dtype.names
        self.nobjects = len(raw_data.field(self.names[0]))
        if (_nobjects_only): return
        # The pyfits raw_data is a FITS_rec object, which isn't picklable, so we need to
        # copy the fields into a new structure to make sure our Catalog is picklable.
        # The simplest is probably a dict keyed by the field names, which we save as self.data.
        self.data = {}
        for name in self.names:
            self.data[name] = raw_data.field(name)
        self.ncols = len(self.names)
        self.isfits = True

    def get(self, index, col):
        """Return the data for the given `index` and `col` in its native type.

        For ASCII catalogs, `col` is the column number.
        For FITS catalogs, `col` is a string giving the name of the column in the FITS table.

        Also, for ASCII catalogs, the "native type" is always str.  For FITS catalogs, it is
        whatever type is specified for each field in the binary table.
        """
        if self.isfits:
            if col not in self.names:
                raise KeyError("Column %s is invalid for catalog %s"%(col,self.file_name))
            if index < 0 or index >= self.nobjects:
                raise IndexError("Object %d is invalid for catalog %s"%(index,self.file_name))
            if index >= len(self.data[col]):
                raise IndexError("Object %d is invalid for column %s"%(index,col))
            return self.data[col][index]
        else:
            icol = int(col)
            if icol < 0 or icol >= self.ncols:
                raise IndexError("Column %d is invalid for catalog %s"%(icol,self.file_name))
            if index < 0 or index >= self.nobjects:
                raise IndexError("Object %d is invalid for catalog %s"%(index,self.file_name))
            return self.data[index, icol]

    def getFloat(self, index, col):
        """Return the data for the given `index` and `col` as a float if possible
        """
        return float(self.get(index,col))

    def getInt(self, index, col):
        """Return the data for the given `index` and `col` as an int if possible
        """
        return int(self.get(index,col))

    def __repr__(self):
        s = "galsim.Catalog(file_name=%r, file_type=%r"%(self.file_name, self.file_type)
        if self.comments != '#':
            s += ', comments=%r'%self.comments
        if self.hdu != 1:
            s += ', hdu=%r'%self.hdu
        s += ')'
        return s

    def __str__(self): return "galsim.Catalog(file_name=%r)"%self.file_name

    def __eq__(self, other): return repr(self) == repr(other)
    def __ne__(self, other): return not self.__eq__(other)
    def __hash__(self): return hash(repr(self))


class Dict(object):
    """A class that reads a python dict from a file.

    After construction, it behaves like a regular python dict, with one exception.
    In order to facilitate getting values in a hierarchy of fields, we allow the '.'
    character to chain keys together for the get() method.  So,

        >>> d.get('noise.properties.variance')

    is expanded into

        >>> d['noise']['properties']['variance']

    Furthermore, if a "key" is really an integer, then it is used as such, which accesses
    the corresponding element in a list.  e.g.

        >>> d.get('noise_models.2.variance')

    is equivalent to

        >>> d['noise_models'][2]['variance']

    This makes it much easier to access arbitrary elements within parameter files.

    Caveat: The above prescription means that an element whose key really has a '.' in it
    won't be accessed correctly.  This is probably a rare occurrence, but the workaround is
    to set `key_split` to a different character or string and use that to chain the keys.


    @param file_name    Filename storing the dict.
    @param dir          Optionally a directory name can be provided if `file_name` does not
                        already include it. [default: None]
    @param file_type    Options are 'Pickle', 'YAML', or 'JSON' or None.  If None, infer from
                        `file_name` extension ('.p*', '.y*', '.j*' respectively).
                        [default: None]
    @param key_split    The character (or string) to use to split chained keys.  (cf. the
                        description of this feature above.)  [default: '.']
    """
    _req_params = { 'file_name' : str }
    _opt_params = { 'dir' : str , 'file_type' : str, 'key_split' : str }
    _single_params = []
    _takes_rng = False

    def __init__(self, file_name, dir=None, file_type=None, key_split='.'):

        # First build full file_name
        self.file_name = file_name.strip()
        if dir is not None:
            import os
            self.file_name = os.path.join(dir,self.file_name)

        if file_type is None:
            import os
            name, ext = os.path.splitext(self.file_name)
            if ext.lower().startswith('.p'):
                file_type = 'PICKLE'
            elif ext.lower().startswith('.y'):
                file_type = 'YAML'
            elif ext.lower().startswith('.j'):
                file_type = 'JSON'
            else:
                raise ValueError('Unable to determine file_type from file_name ending')

        file_type = file_type.upper()
        if file_type not in ['PICKLE','YAML','JSON']:
            raise ValueError("file_type must be one of Pickle, YAML, or JSON if specified.")
        self.file_type = file_type

        self.key_split = key_split

        if file_type == 'PICKLE':
            try:
                import cPickle as pickle
            except ImportError:
                import pickle
            with open(self.file_name, 'rb') as f:
                self.dict = pickle.load(f)
        elif file_type == 'YAML':
            import yaml
            with open(self.file_name, 'r') as f:
                self.dict = yaml.load(f)
        elif file_type == 'JSON':
            import json
            with open(self.file_name, 'r') as f:
                self.dict = json.load(f)
        else:
            raise ValueError("Invalid file_type %s"%file_type)

    def get(self, key, default=None):
        # Make a list of keys according to our key_split parameter
        chain = key.split(self.key_split)
        d = self.dict
        while len(chain):
            k = chain.pop(0)

            # Try to convert to an integer:
            try: k = int(k)
            except ValueError: pass

            # If there are more keys, just set d to the next in the chanin.
            if chain: d = d[k]
            # Otherwise, return the result.
            else:
                if k not in d and default is None:
                    raise ValueError("key=%s not found in dictionary"%key)
                return d.get(k,default)

        raise ValueError("Invalid key=%s given to Dict.get()"%key)

    # The rest of the functions are typical non-mutating functions for a dict, for which we just
    # pass the request along to self.dict.
    def __len__(self):
        return len(self.dict)

    def __getitem__(self, key):
        return self.dict[key]

    def __contains__(self, key):
        return key in self.dict

    def __iter__(self):
        return self.dict.__iter__

    def keys(self):
        return self.dict.keys()

    def values(self):
        return self.dict.values()

    def items(self):
        return self.dict.items()

    def iterkeys(self):
        return iterkeys(self.dict)

    def itervalues(self):
        return itervalues(self.dict)

    def iteritems(self):
        return iteritems(self.dict)

    def __repr__(self):
        s = "galsim.Dict(file_name=%r, file_type=%r"%(self.file_name, self.file_type)
        if self.key_split != '.':
            s += ', key_split=%r'%self.key_split
        s += ')'
        return s

    def __str__(self): return "galsim.Dict(file_name=%r)"%self.file_name

    def __eq__(self, other): return repr(self) == repr(other)
    def __ne__(self, other): return not self.__eq__(other)
    def __hash__(self): return hash(repr(self))



class OutputCatalog(object):
    """A class for building up a catalog for output, typically storing truth information
    about a simulation.

    Each row corresponds to a different object, and each column stores some item of
    information about that object (e.g. flux or half_light_radius).

    Note: no type checking is done when the data are added in addRow().  It is up to
    the user to make sure that the values added for each row are compatible with the
    types given here in the `types` parameter.

    Initialization
    --------------

    @param names    A list of names for the output columns.
    @param types    A list of types for the output columns. [default: None, which assumes all
                    columns are float]

    Attributes
    ----------

    After construction, the following attributes are available:

        nobjects    The number of objects so far in the catalog.
        ncols       The number of columns in the catalog.
        names       The names of the columns.
        types       The types of the columns.
        rows        The rows of data that have been accumulated so far.

    """
    # Watch out for this "Gotcha".  Using _rows=[] as the default argument doesn't work!
    # https://pythonconquerstheuniverse.wordpress.com/2012/02/15/mutable-default-arguments/
    def __init__(self, names, types=None, _rows=(), _sort_keys=()):
        self.names = names
        if types is None:
            self.types = [ float for i in names ]
        else:
            self.types = types
        self.rows = list(_rows)
        self.sort_keys = list(_sort_keys)

    @property
    def nobjects(self): return len(self.rows)
    @property
    def ncols(self): return len(self.names)

    # Again, when we use this through a proxy, we need getters for the attributes.
    def getNames(self): return self.names
    def getTypes(self): return self.types
    def setTypes(self, types): self.types = types
    def getNObjects(self): return self.nobjects
    def getNCols(self): return self.ncols

    def addRow(self, row, sort_key=None):
        """Add a row of data to the catalog.

        Warning: no type checking is done at this point.  If the values in the row do not
        match the column types, you may get an error when writing, or you may lose precision,
        depending on the nature of the mismatch.

        @param row          A list with one item per column in the same order as the names list.
        @param sort_key     If the rows may be added out of order, you can provide a sort_key,
                            which will be used at the end to re-sort the rows.
        """
        if len(row) != self.ncols:
            raise ValueError("Length of row does not match the number of columns")
        self.rows.append(tuple(row))
        if sort_key is None:
            self.sort_keys.append(self.nobjects)
        else:
            self.sort_keys.append(sort_key)

    def write(self, file_name, dir=None, file_type=None, prec=8):
        """Write the catalog to a file.

        @param file_name    The name of the file to write to.
        @param dir          Optionally a directory name can be provided if `file_name` does not
                            already include it. [default: None]
        @param file_type    Which kind of file to write to. [default: determine from the file_name
                            extension]
        @param prec         Output precision for ASCII. [default: 8]
        """
        if dir is not None:
            import os
            file_name = os.path.join(dir,file_name)

        # Figure out which file type the catalog is
        if file_type is None:
            import os
            name, ext = os.path.splitext(file_name)
            if ext.lower().startswith('.fit'):
                file_type = 'FITS'
            else:
                file_type = 'ASCII'
        file_type = file_type.upper()
        if file_type not in ['FITS', 'ASCII']:
            raise ValueError("file_type must be either FITS or ASCII if specified.")

        if file_type == 'FITS':
            self.writeFits(file_name)
        elif file_type == 'ASCII':
            self.writeAscii(file_name, prec)
        else:
            raise ValueError("Invalid file_type %s"%file_type)

    def makeData(self):
        """Returns a numpy array of the data as it should be written to an output file.
        """
        cols = zip(*self.rows)

        dtypes = []
        new_cols = []
        for col, name, t in zip(cols, self.names, self.types):
            name = str(name)  # numpy will barf if the name is a unicode string
            dt = np.dtype(t) # just used to categorize the type into int, float, str
            if dt.kind in np.typecodes['AllInteger']:
                dtypes.append( (name, int) )
                new_cols.append(col)
            elif dt.kind in np.typecodes['AllFloat']:
                dtypes.append( (name, float) )
                new_cols.append(col)
            elif t == galsim.Angle:
                dtypes.append( (name + ".rad", float) )
                new_cols.append( [ val.rad for val in col ] )
            elif t == galsim.PositionI:
                dtypes.append( (name + ".x", int) )
                dtypes.append( (name + ".y", int) )
                new_cols.append( [ val.x for val in col ] )
                new_cols.append( [ val.y for val in col ] )
            elif t == galsim.PositionD:
                dtypes.append( (name + ".x", float) )
                dtypes.append( (name + ".y", float) )
                new_cols.append( [ val.x for val in col ] )
                new_cols.append( [ val.y for val in col ] )
            elif t == galsim.Shear:
                dtypes.append( (name + ".g1", float) )
                dtypes.append( (name + ".g2", float) )
                new_cols.append( [ val.g1 for val in col ] )
                new_cols.append( [ val.g2 for val in col ] )
            else:
                col = [ str(s).encode() for s in col ]
                maxlen = np.max([ len(s) for s in col ])
                dtypes.append( (name, str, maxlen) )
                new_cols.append(col)

        data = np.array(list(zip(*new_cols)), dtype=dtypes)

        sort_index = np.argsort(self.sort_keys)
        data = data[sort_index]

        return data

    def writeAscii(self, file_name, prec=8):
        """Write catalog to an ASCII file.

        @param file_name    The name of the file to write to.
        @param prec         Output precision for floats. [default: 8]
        """
        data = self.makeData()

        width = prec+8
        header_form = ""
        for i in range(len(data.dtype.names)):
            header_form += "{%d:^%d} "%(i,width)
        header = header_form.format(*data.dtype.names)

        fmt = []
        for name in data.dtype.names:
            dt = data.dtype[name]
            if dt.kind in np.typecodes['AllInteger']:
                fmt.append('%%%dd'%(width))
            elif dt.kind in np.typecodes['AllFloat']:
                fmt.append('%%%d.%de'%(width,prec))
            else:
                fmt.append('%%%ds'%(width))

        try:
            np.savetxt(file_name, data, fmt=fmt, header=header)
        except (AttributeError, TypeError):  # pragma: no cover
            # header was added with version 1.7, so do it by hand if not available.
            with open(file_name, 'w') as fid:
                fid.write('#' + header + '\n')
                np.savetxt(fid, data, fmt=fmt)

    def writeFits(self, file_name):
        """Write catalog to a FITS file.

        @param file_name    The name of the file to write to.
        """
        tbhdu = self.writeFitsHdu()
        galsim.fits.writeFile(file_name, tbhdu)

    def writeFitsHdu(self):
        """Write catalog to a FITS hdu.

        @returns an HDU with the FITS binary table of the catalog.
        """
        # Note to developers: Because of problems with pickling in older pyfits versions, this
        # code is duplicated in galsim/config/extra_truth.py, BuildTruthHDU.  If you change
        # this function, you should update BuildTruthHDU as well.
        from galsim._pyfits import pyfits

        data = self.makeData()

        cols = []
        for name in data.dtype.names:
            dt = data.dtype[name]
            if dt.kind in np.typecodes['AllInteger']:
                cols.append(pyfits.Column(name=name, format='J', array=data[name]))
            elif dt.kind in np.typecodes['AllFloat']:
                cols.append(pyfits.Column(name=name, format='D', array=data[name]))
            else:
                cols.append(pyfits.Column(name=name, format='%dA'%dt.itemsize, array=data[name]))

        cols = pyfits.ColDefs(cols)

        # Depending on the version of pyfits, one of these should work:
        try:
            tbhdu = pyfits.BinTableHDU.from_columns(cols)
        except AttributeError:  # pragma: no cover
            tbhdu = pyfits.new_table(cols)
        return tbhdu

    def __repr__(self):
        def make_type_str(t):
            s = repr(t)
            if s[1:5] == 'type': return s[7:-2]
            elif s[1:6] == 'class': return s[8:-2]
            else: return s
        type_str = "( " + ", ".join([ make_type_str(t) for t in self.types ]) + " )"
        return "galsim.OutputCatalog(names=%r, types=%s, _rows=%r)"%(
                self.names, type_str, self.rows)

    def __str__(self):
        return "galsim.OutputCatalog(name=%r)"%self.names

    def __eq__(self, other): return repr(self) == repr(other)
    def __ne__(self, other): return not self.__eq__(other)
    def __hash__(self): return hash(repr(self))

