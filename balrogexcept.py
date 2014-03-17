# 101 --imagein does not exist
# 102 --weight in does not exist
# 103 --psfin does not exist
# 104 --psfin not given
#
# 111 --imagein is not a FITS file
# 112 --weightin is not a FITS file
# 113 --psfin is not a FITS file
#
# 121 --imageext does not exist
# 122 --weightext does not exist

# 131 image missing NAXIS1
# 132 image missing NAXIS2
# 133 weight missing NAXIS1
# 134 weight missing NAXIS2
# 135 image and weight have different dimensions
# 136 xmin larger than xmax
# 137 ymin larger than ymax

# 140 sextractor path does not exist

# 201 Cannot create --outdir
# 202 Cannot create output subdirectory


#### config.py syntax errors

# 305 asked for attribute rules.%s_sersic.%s
# 405 asked for attribute sampled.%s_sersic.%s

# 401 asked for attribute of sampled that doesn't exist
# 301 asked for attribute of rules that doesn't exist

# 402 sampled itself has no indexing, and even if it did it wouldn't be reassignable
# 302 rules itself has no indexing

# 403 sampled reassignment
# 303 rule not understood

# 404 sampled index out of range
# 304 rules index out of range

# 306 tried to assign one of the sersic components of rules to something other than an array

#  -1 ngal direct change attempt
#  -2 nProfiles direct change attempt

# 501 Function return something other than array of the right length
# 502 Funciton arg error
# 503 Catalog arg error
# 504 Catalog file doesn't exist
# 505 Catalog ext doesn't exist
# 506 Catalog col doesn't exist


class BaseException(Exception):
    def __init__(self, *args):
        self.code = args[0]
        arguments = args[1:]
        self.init(*arguments)

    def __str__(self):
        return repr(self.msg)


class OutdirWriteError(BaseException):
    def init(self, dir, fdir):
        self.msg = "ERROR code: %i. Attempting to create --outdir %s failed. Could not create %s" %(self.code,dir,fdir)


class SubdirWriteError(BaseException):
    def init(self, dir):
        self.msg = "ERROR code: %i. Could not create output subdirectory %s" %(self.code,dir)


class ImageInputError(BaseException):
    def init(self, label, arg, val):
        self.msg = "ERROR code: %i. Given input %s file does not exist: --%s %s" %(self.code, label, arg, val)


class PsfInputError(BaseException):
    def init(self, image):
        self.msg = "ERROR code: %i. Input image file other than default was given: --imagein %s, but no --psfin was given. No sensible default exists for assuming a --psfin." %(self.code, image)


class FitsFileError(BaseException):
    def init(self, label, arg, val):
        self.msg = "ERROR code: %i. Given input %s file was not recognized as a FITS file by pyfits and could not be opened: --%s %s" %(self.code, label, arg, val)


class FitsExtError(BaseException):
    def init(self, label, arg, val, iarg,ival):
        self.msg = "ERROR code: %i. Given input %s extension does not exist: --%s %s, --%s %s" %(self.code, label, arg,val, iarg,ival)


class FitsHeaderError(BaseException):
    def init(self, label,keyword, iarg,ival, earg,eval):
        self.msg = "ERROR code: %i. Given input %s extension was not recognized as image because it is missing keyword %s: --%s %s, --%s %s" %(self.code,label,keyword, iarg,ival, earg,eval)


class SizeMismatchError(BaseException):
    def init(self, ix,iy, wx,wy):
        self.msg = "ERROR code: %i. Input image and input weight have different dimensions (col, row): image = (%i, %i), weight = (%i, %i)" %(self.code, ix,iy,wx,wy)


class SizeError(BaseException):
    def init(self, label, minimum,maximum):
        self.msg = "ERROR code: %i. --%smin = %i is greater than --%smax = %i" %(self.code, label,minimum, label,maximum)


class SampledAttributeError(BaseException):
    def init(self, name):
        self.msg = 'ERROR code: %i. Asked for an attribute of sampled which does not exist: sampled.%s. Only sampled.{x,y,g1,g2,magnification,sersicindex,halflightradius,magnitude,axisratio,beta} are valid.' %(self.code, name)

class RulesAttributeError(BaseException):
    def init(self, name):
        self.msg = 'ERROR code: %i. Asked for an attribute of rules which does not exist: rules.%s. Only rules.{x,y,g1,g2,magnification,sersicindex,halflightradius,magnitude,axisratio,beta} are valid.' %(self.code, name)

class RulesComponentAttributeError(BaseException):
    def init(self):
        self.msg = 'ERROR code: %i. Sersic components of rules do not have any attributes. So you cannot get or set any' %(self.code)

class SampledComponentAttributeError(BaseException):
    def init(self):
        self.msg = 'ERROR code: %i. Sersic components of sampled do not have any attributes. So you cannot get or set any' %(self.code)


class SampledIndexingError(BaseException):
    def init(self):
        self.msg = 'ERROR code: %i. The object sampled itself has no indexing. Only sampled.{sersicindex,halflightradius,magnitude,axisratio,beta} are indexed.' %(self.code)

class RulesIndexingError(BaseException):
    def init(self):
        self.msg = 'ERROR code: %i. The object rules itself has no indexing. Only rules.{sersicindex,halflightradius,magnitude,axisratio,beta} are indexed.' %(self.code)


class SampledAssignmentError(BaseException):
    def init(self, name):
        self.msg = 'ERROR code: %i. Attempted to (re)assign sampled.%s. You cannot (re)assign any attributes of the sampling results. This is only to be used to access the results. Assign rules to change results.' %(self.code, name)

class RulesAssignmentError(BaseException):
    def init(self, label):
        self.msg = 'ERROR code: %i. The rule you gave for %s was not understood. Rules can be a single value, an array of length ngal, an attribute of sampled, a Catalog() statement or a Function() statement.' %(self.code, label)


class SampledIndexOutOfRange(BaseException):
    def init(self, name, size):
        self.msg = "ERROR code: %i. Attempted to index sampled.%s beyond it's size of %i." %(self.code, name, size)

class RulesIndexOutOfRange(BaseException):
    def init(self, name, size):
        self.msg = "ERROR code: %i. Attempted to index rules.%s beyond it's size of %i." %(self.code, name, size)


class RulesAssignmentNoArrayError(BaseException):
    def init(self):
        self.msg = 'ERROR code: %i. Attempted an illegal reassingment of sersic component rules where the number of components is greater than 1. These must be arrays of the same length as the number of profiles' %(self.code)


class RulesnProfilesError(BaseException):
    def init(self, name):
        self.msg = "ERROR code: %i. You've deduced that rules has an attribute called rules.nProfiles. However, you're not allowed to change it directly. Use InitializeSersic()" %(self.code)


class RulesNgalError(BaseException):
    def init(self):
        self.msg = "ERROR code: %i. You've deduced that rules has an attribute called rules.ngal. However, you're not allowed to change it." %(self.code)


class FunctionReturnError(BaseException):
    def init(self, name):
        self.msg = "ERROR code: %i. Illegal return value for the function %s. All Function rules must return an array of length ngal." %(self.code, name)

class FunctionArgError(BaseException):
    def init(self, label):
        self.msg = "ERROR code: %i. Must specify %s with sampling type Function" %(self.code, label)

class CatalogArgError(BaseException):
    def init(self, label):
        self.msg = "ERROR code: %i. Must specify %s with sampling type Catalog" %(self.code, label)

class CatalogFileError(BaseException):
    def init(self, file):
        self.msg = "ERROR code: %i. Catalog file specified does not exist: %s." %(self.code, file)

class CatalogExtError(BaseException):
    def init(self, file, ext):
        self.msg = "ERROR code: %i. Catalog extenstion specified does not exist: %i, %s." %(self.code, ext, file)

class CatalogColError(BaseException):
    def init(self, file, ext, col):
        self.msg = "ERROR code: %i. Catalog column specified does not exist: %s, %s[%i]." %(self.code, col, file, ext)


class SextractorPathError(BaseException):
    def init(self, path):
        self.msg = "ERROR code: %i. Path to sextractor does not exist: %s." %(self.code, path)
