import nibabel as nib
import numpy as np
from sys import argv
#import ipdb

def andb(arrs):
	result = np.zeros(arrs[0].shape)
	for aa in arrs: result+=np.array(aa,dtype=np.int)
	return result

def niwrite(data,affine, name , header=None):
	#stdout.write(" + Writing file: %s ...." % name) 
	
	thishead = header
	if thishead == None:
		thishead = head.copy()
		thishead.set_data_shape(list(data.shape))

	outni = nib.Nifti1Image(data,affine,header=thishead)
	outni.to_filename(name)
	print 'done.'

#syntax is:
#leftsurfvolts rightsurfvolts affvolts subcortmask templatesegmentation nativesegmentation prefix

lhvol = nib.load(argv[1])
aff = lhvol.get_affine()
head = lhvol.get_header()
lhd = np.squeeze(lhvol.get_data())
rhd = np.squeeze(nib.load(argv[2]).get_data())
vold = nib.load(argv[3]).get_data()
subm = nib.load(argv[4]).get_data()
tsegm = nib.load(argv[5]).get_data()
nsegm = nib.load(argv[6]).get_data()

#ipdb.set_trace()

out = np.zeros(lhd.shape)
wmmask = andb([andb([tsegm==1,nsegm==1])==2,andb([tsegm==21,nsegm==21])==2])==1
out[wmmask,:] = vold[wmmask,:]
lhdmask = lhd.mean(-1)!=0
rhdmask = rhd.mean(-1)!=0
out[lhdmask!=0,:] = lhd[lhdmask!=0,:]
out[rhdmask!=0,:] = rhd[rhdmask!=0,:]
out[andb([lhdmask!=0,rhdmask!=0])>1] = 0
out[subm!=0,:]=vold[subm!=0,:]

niwrite(out,aff,argv[7])

