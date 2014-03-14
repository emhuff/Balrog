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

# 201 Cannot create --outdir
# 202 Cannot create output subdirectory


class OutdirWriteError(Exception):

    def __init__(self, code, dir, fdir):
        self.code = code
        self.msg = "ERROR code: %i. Attempting to create --outdir %s failed. Could not create %s" %(self.code,dir,fdir)

    def __str__(self):
        return repr(self.msg)


class SubdirWriteError(Exception):

    def __init__(self, code, dir):
        self.code = code
        self.msg = "ERROR code: %i. Could not create output subdirectory %s" %(self.code,dir)

    def __str__(self):
        return repr(self.msg)



class ImageInputError(Exception):

    def __init__(self, code, label, arg, val):
        self.code = code
        self.msg = "ERROR code: %i. Given input %s file does not exist: --%s %s" %(self.code, label, arg, val)

    def __str__(self):
        return repr(self.msg)


class PsfInputError(Exception):

    def __init__(self, code, image):
        self.code = code
        self.msg = "ERROR code: %i. Input image file other than default was given: --imagein %s, but no --psfin was given. No sensible default exists for assuming a --psfin." %(self.code, image)

    def __str__(self):
        return repr(self.msg)


class FitsFileError(Exception):
    
    def __init__(self, code, label, arg, val):
        self.code = code
        self.msg = "ERROR code: %i. Given input %s file was not recognized as a FITS file by pyfits and could not be opened: --%s %s" %(self.code, label, arg, val)

    def __str__(self):
        return repr(self.msg)


class FitsExtError(Exception):
    
    def __init__(self, code, label, arg, val, iarg,ival):
        self.code = code
        self.msg = "ERROR code: %i. Given input %s extension does not exist: --%s %s, --%s %s" %(self.code, label, arg,val, iarg,ival)

    def __str__(self):
        return repr(self.msg)


class FitsHeaderError(Exception):
    
    def __init__(self, code,label,keyword, iarg,ival, earg,eval):
        self.code = code
        self.msg = "ERROR code: %i. Given input %s extension was not recognized as image because it is missing keyword %s: --%s %s, --%s %s" %(self.code,label,keyword, iarg,ival, earg,eval)

    def __str__(self):
        return repr(self.msg)


class SizeMismatchError(Exception):
    
    def __init__(self, code, ix,iy, wx,wy):
        self.code = code
        self.msg = "ERROR code: %i. Input image and input weight have different dimensions (col, row): image = (%i, %i), weight = (%i, %i)" %(code, ix,iy,wx,wy)

    def __str__(self):
        return repr(self.msg)


class SizeError(Exception):
    
    def __init__(self, code, label, minimum,maximum):
        self.code = code
        self.msg = "ERROR code: %i. --%smin = %i is greater than --%smax = %i" %(code, label,minimum, label,maximum)

    def __str__(self):
        return repr(self.msg)
