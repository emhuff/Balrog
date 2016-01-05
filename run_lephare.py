# To run this routine, you need to install Le Phare package that you can find at:
# http://www.cfht.hawaii.edu/~arnouts/LEPHARE/download.html
# You need to define LEPHAREDIR and LEPHAREWORK paths as explained in the INSTALL file of Le Phare package
#
# Add the DESjul13 and vista repertories in ~/wherever_lephare_is/lephare_dev/filt/
#
# setenv LEPHAREWORK '/share/splinter/sjouvel/SIMUL/LePhare/WORK/'
# LEPHAREWORK will be used automatically by Le Phare to put logs and data Le Phare needs while running
# 

import pdb
import pyfits as pf
import numpy as np
import os
import subprocess

########### PATHS ##############
balrog_dir=os.environ['BALROGDIR']
lephare_dir=os.environ['LEPHAREDIR']
catalogs_dir="default_example/output/balrog_cat/"
input_dir=os.path.join(balrog_dir,catalogs_dir)
lephare_output=balrog_dir+'/lephare_output/'

########### INPUT/OUTPUT FILES ##############
lephare_config = balrog_dir+"/astro_config/zphot.para"
lp_output_para=balrog_dir+"astro_config/zphot_output_PDZ.para"
# temporarily reading cosmos.fits
in_file=os.path.join(balrog_dir+"cosmos.fits")
ascii_file=os.path.join(in_file+".ascii")
out_file=os.path.join(in_file+".ascii.lephare")

########### MAGNITUDES AND PRIORS TO RUN LE PHARE ##############
mag_names = ['MAPP_G_SUBARU','MAPP_R_SUBARU','MAPP_I_SUBARU','MAPP_Z_SUBARU','MAPP_HSC_Y']
emag_names = ['ERR_G_SUBARU','ERR_R_SUBARU','ERR_I_SUBARU','ERR_Z_SUBARU','ERR_HSC_Y']
nz_prior = '3,1,3'
adapt = 'YES'

def make_input_catalogue():
# formatting catalogues for Le Phare input and save it in ascii
	print "=====> Creating input catalogue for Le Phare <====="
	out = []
	in_catalog= pf.open(in_file)
	# object id
	obj_id = in_catalog[1].columns[0].name
	out.extend([in_catalog[1].data[obj_id]])
	# magnitudes, errors
	for F in range(len(mag_names)):
	        out.extend([in_catalog[1].data[mag_names[F]],in_catalog[1].data[emag_names[F]]])
	# context parameter => see Le Phare documentation
	context = [0] * len(in_catalog[1].data[obj_id]) + np.power(2,len(mag_names)) - 1 
	out.extend([context])
	# specz if available
	specz = in_catalog[1].data['Z']
	out.extend([specz])
	# SAVE TABLE IN ASCII
	out_tmp = np.vstack(out).T
	fmt = '%d '+'%f '*2*len(mag_names)+'%d %f'
	np.savetxt(ascii_file, out_tmp, fmt=fmt)

def lephare():
##### LE PHARE RUN
	print "=====> Running Le Phare <====="
	# create a repertory to put Le Phare libraries
	if (os.path.isdir(lephare_output) == False): subprocess.call('mkdir '+lephare_output,shell='True')
	# running Le Phare
	#cmd = lephare_dir+"/source/sedtolib -t G -c "+lephare_config 
	#subprocess.call(cmd,shell='True') 
	#cmd = lephare_dir+"/source/sedtolib -t S -c "+lephare_config 
	#subprocess.call(cmd,shell='True') 
	#cmd = lephare_dir+"/source/sedtolib -t Q -c "+lephare_config 
	#subprocess.call(cmd,shell='True') 
	#cmd = lephare_dir+"/source/filter -t G -c "+lephare_config 
	#subprocess.call(cmd,shell='True') 
	#cmd = lephare_dir+"/source/mag_gal -t G -c "+lephare_config 
	#subprocess.call(cmd,shell='True') 
	#cmd = lephare_dir+"/source/mag_star -c "+lephare_config 
	#subprocess.call(cmd,shell='True') 
	#cmd = lephare_dir+"/source/mag_gal -t Q -c "+lephare_config
	#subprocess.call(cmd,shell='True') 
	cmd = lephare_dir+"/source/zphota -t G -c "+lephare_config+" -CAT_IN "+ascii_file+" -CAT_OUT "+out_file+" -PARA_OUT "+lp_output_para+ " -AUTO_ADAPT "+adapt+" -NZ_PRIOR "+nz_prior 
	subprocess.call(cmd,shell='True') 
	# moving diretories to lephare_output directory
	cmd = 'mv *.dat *lephare* '+lephare_output
	subprocess.call(cmd,shell='True')

def main():
    #make_input_catalogue()
    lephare()

if __name__ == '__main__':
    main()

