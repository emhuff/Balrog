"""Some methods to get metadata for multi-epoch run on a DES tile
We need:
i) List of SE image paths
For each of these we need:
a) psfEx file path
b) SCAMP astrometry file
c) gain
d) zeropoint"""
import desdb
import numpy as np
from string import Template
import subprocess

COADD_TABLE="Y1A1_COADD"
query_template=Template("select $columns from $tables where $conditions")
balrog_bin="/home/maccrann/code/Balrog/balrog.py"
funpack_bin="funpack"

def funpack(fz_file):
    subprocess.call([funpack_bin,fz_file])

class Coadd(desdb.files.Coadd):
    """The desdb.files.Coadd class already pretty much does what we need"""
    """Here I've just supplemented the _load_srclist method so it also gets the gain(s!)"""
    def load_psf_path(self,psf_path=None):
        #First check we've loaded the image path...we need this
        if psf_path is None:
            try:
                image_path=self['image_url']
            except KeyError:
                print "image path not set, need this to set default psf path"
                raise KeyError
            self['psf_path']=image_path.replace(".fits.fz","_psfcat.psf")
        else:
            self['psf_path']=psf_path            

    def load(self, srclist=False):
        super(Coadd,self).load(srclist=srclist)
        self['det_image_url']=self['image_url'].replace(self['band']+".fits","det.fits")

    def _load_srclist(self):
        """Same as original, but adding gaina and gainb, and whatever else turns out 
        to be imporant"""

        """
        This new query works because the new system is more sane
            - guaranteed number of levels between coadd and SE image
            - single coadd run for all bands
        
        Thanks to Bob Armstrong for the new query
        Note for non-psf homogenized coadds, will are using one less level.
        See Bob's email.
        """

        desdata=desdb.files.get_des_rootdir()
        query="""
        SELECT
            magzp,
            coadd.band as band,
            d.id,
            d.gaina,
            d.gainb,
            loc.run,
            loc.exposurename as expname,
            loc.ccd
        FROM
            coadd_src,coadd,image c,image d, location loc
        WHERE
            coadd.band='{band}'
            and coadd_src.coadd_imageid=coadd.id
            and coadd.run='{coadd_run}'
            and c.id=coadd_src.src_imageid
            and c.parentid=d.id
            and loc.id = d.id\n"""


        query=query.format(band=self['band'],
                           coadd_run=self['coadd_run'])

        res = self.conn.quick(query, show=self.verbose)

        df=desdb.DESFiles(fs=self.fs)
        srclist=[]
        for r in res:
            for type in ['image','bkg','seg','cat']:
                ftype='red_%s' % type
                url=df.url(ftype,
                           run=r['run'],
                           expname=r['expname'],
                           ccd=r['ccd'])
                r[ftype] = url

            r['astro_refine'] = df.url('astro_refine',
                                       coadd_run=self['coadd_run'],
                                       expname=r['expname'],
                                       ccd=r['ccd'])


            srclist.append(r)

        self.srclist=srclist
        return

    def funpack(self):
        subprocess.call([funpack_bin,self['image_url']])
        self['image_url_funpacked']=(self['image_url']).rsplit('.',1)[0]        

    def call_balrog(self,**kwargs):
        #need to funpack coadd to run sextractor on it
        subprocess.call([funpack_bin,self['image_url']])
        self['image_url_funpacked']=(self['image_url']).rsplit('.',1)[0]
        print 'funpacked image path:',self['image_url_funpacked']

        balrog_args=['-ie','0','-we','1','-ft','-n','100']
        #the essentials
        balrog_args+=['-i',self['image_url_funpacked']]
        balrog_args+=['-w',self['image_url_funpacked']]
        balrog_args+=['-p',self['psf_path']]
        for key,value in kwargs.iteritems():
            if type(value)==bool:
                balrog_args.append("-%s"%key)
            else:
                balrog_args+=["-%s"%key,value]
        subprocess_args=[balrog_bin]+balrog_args
        print subprocess_args
        subprocess.call(subprocess_args)

def tilename_to_id(conn,tilename,band):
    query=query_template.substitute(columns="id",tables=COADD_TABLE,conditions="tilename='%s' and band='%s'"%(tilename,band))
    d=conn.quick(query)
    return d[0]['id']

def sync_coadd(run):
    """Use desdb des-sync-coadd to sync coadd stuff
    This actually syncs all the bands and all the astro_refine headers which is overkill really...
    But can be more selective when optimizing all this."""
    subprocess.call(["./des-sync-coadd-nm",run]) #Gets the coadd images, catalogs and psf catalogs
    subprocess.call(["./des-sync-coadd-nm","-a",run]) #Gets the astro_refine headers too
    subprocess.call(["./des-sync-coadd-nm","-q",run]) #Gets the segmaps

    
def setup_tile(tilename,band='i',sync=False):
    
    conn=desdb.Connection()
    id=tilename_to_id(conn,tilename,band)

    #Make a coadd object
    c = Coadd(id=id)

    c.load(srclist=True)
    #Sync coadd if necessary
    if sync:
        sync_coadd(c['coadd_run'])

    #Set psf path
    c.load_psf_path()
    return c

def balrog_SE(srclist_entry,**kwargs):
    balrog_args=['-ie','0','-we','1','-ft','-n','5']
    subprocess.call([funpack_bin,srclist_entry['red_image']])
    srclist_entry['red_image_funpacked']=(srclist_entry['red_image']).rsplit('.',1)[0]
    balrog_args+=['-i',srclist_entry['red_image_funpacked']]
    balrog_args+=['-w',srclist_entry['red_image_funpacked']]
    balrog_args+=['-p',srclist_entry['red_image'].replace(".fits.fz","_psfcat.psf")]
    balrog_args+=['-sh',srclist_entry['astro_refine']]
    
    for key,value in kwargs.iteritems():
        if type(value)==bool:
            balrog_args.append("-%s"%key)
        else:
            balrog_args+=["-%s"%key,value]
    subprocess_args=[balrog_bin]+balrog_args
    print subprocess_args
    subprocess.call(subprocess_args)

def main():

    tilename="DES0356-5331"   #One of the Y1A1 testbed tiles!
    band="i"
    c=setup_tile(tilename,band='i')
    #Run balrog on coadd...
    #c.call_balrog(o="coadd",pc="config_coadd.py")
    balrog_SE(c.srclist[0],o="se",pc="config.py")

if __name__=="__main__":
    main()
    
