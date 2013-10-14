import sys
import commands
from re import split as resplit
import re
from os import system,getcwd,mkdir,chdir,popen
import os.path
from string import rstrip,split
from optparse import OptionParser,OptionGroup
import cPickle as pickle

execdir = os.path.dirname(sys.argv[0])

#Filename parser for NIFTI and AFNI files
def dsprefix(idn):
	def prefix(datasetname):
		return split(datasetname,'+')[0]
	if len(split(idn,'.'))!=0:
		if split(idn,'.')[-1]=='HEAD' or split(idn,'.')[-1]=='BRIK' or split(idn,'.')[-2:]==['BRIK','gz']:
			return prefix(idn)
		elif split(idn,'.')[-1]=='nii' and not split(idn,'.')[-1]=='nii.gz':
			return '.'.join(split(idn,'.')[:-1])
		elif split(idn,'.')[-2:]==['nii','gz']:
			return '.'.join(split(idn,'.')[:-2])
		else:
			return prefix(idn)
	else:
		return prefix(idn)

def dssuffix(idna):
	suffix = idna.split(dsprefix(idna))[-1]
	spl_suffix=suffix.split('.')
	if len(spl_suffix[0])!=0 and spl_suffix[0][0] == '+': return spl_suffix[0]
	else: return suffix

def parse_mecode(dstr):
	dsinputs=dsprefix(dstr)
	prefix=resplit(r'[\[\],]',dsinputs)[0]
	datasets=resplit(r'[\[\],]',dsinputs)[1:-1]
	trailing=resplit(r'[\]+]',dsinputs)[-1]
	isf= dssuffix(dstr)
	setname = prefix+''.join(datasets)+trailing
	return sorted([prefix+echo+trailing+isf for echo in datasets]),setname

def gohome(go=1): 
	subdir = '%s/%s' % (options.basepath,options.sid)
	if go: cmds.append("cd %s"  % (subdir))
	return subdir

def gomeicadir(go=1):
	homedir = gohome(0)
	meicadir="%s/meica.%s" % (homedir,parse_mecode(options.func)[1])
	if go: cmds.append("cd %s" % (meicadir))
	return meicadir

def gosumadir(go=1):
	homedir = gohome(0)
	sumadir = "%s/%s_SUMA" % (homedir,options.sid)
	if go: cmds.append("cd %s" % sumadir)
	return sumadir

def gofeatsdir(go=1):
	homedir = gohome(0)
	featsdir = "%s/%s_feats.%s" % (homedir,options.sid,parse_mecode(options.func)[1])
	if go: cmds.append("cd %s" % featsdir)
	return featsdir

def gomfeatsdir(go=1):
	featsdir = gofeatsdir(0)
	mfeatsdir = "%s/meica_feats"  % (featsdir)
	if go: cmds.append("cd %s" % mfeatsdir)
	return mfeatsdir

def gosfeatsdir(go=1):
	featsdir = gofeatsdir(0)
	sfeatsdir = "%s/spca_feats"  % (featsdir)
	if go: cmds.append("cd %s" % sfeatsdir)
	return sfeatsdir

def getfsid():
	return options.fs.rstrip('/').split('/')[-1]

def getsessfn():
	return 'clust_%s_%s.pkl' % (options.sid,parse_mecode(options.func)[1])

def getscriptfn():
	return '.clust_%s_%s.sh' % (options.sid,parse_mecode(options.func)[1])

#A grand tour of modern mathematics...
#run_steps_all = ['setup','suma','get_func','meica_pp','tedanaX','fnirt_anat','fnirt_func','surfsmooth','surfcomb','prep_feats','save_space','meica_feats','spca','spca_feats','done']
#run_steps_all = ['setup','suma','get_func','meica_pp','tedanaX','fnirt_anat','fnirt_func','surfsmooth','surfcomb','prep_feats','save_space','meica_feats','done']
run_steps_all = ['setup','suma','get_func','fnirt_anat','meica_pp','z','fnirt_tedanaX','done']
run_steps_str = ','.join(run_steps_all)

from optparse import OptionParser
parser=OptionParser()
parser.add_option('',"--sid",dest='sid',help="Subject ID",default='test')
parser.add_option('',"--sess",dest='sess',help="Load options from session file",default='')
parser.add_option('',"--sourcepath",dest='sourcepath',help="Source directory for functionals",default=os.getcwd())
parser.add_option('',"--basepath",dest='basepath',help="Base path",default=os.getcwd())
parser.add_option('',"--func",dest='func',help="ME functional data code",default='')
parser.add_option('',"--fs",dest='fs',help="FreeSurfer directory",default='')
parser.add_option('',"--TEs",dest='TEs',help="Echo times (comma separated)",default='')
parser.add_option('',"--fwhm",dest='fwhm',help="FWHM smoothing kernel, default 0mm",default='0mm')
parser.add_option('',"--tpattern",dest='tpattern',help="Default read from data.",default='')
parser.add_option('',"--steps",dest='steps',help="Do procedures, default all: %s " % run_steps_str,default=run_steps_str)
parser.add_option('',"--stepin",dest='stepin',help="Step to start from (for resuming, for example)",default='')
parser.add_option('',"--noexec",dest='noexec',action='store_true',help="Just write script, dont run", default=False)
parser.add_option('',"--jobs",dest='jobs',help="Number of CPUs to use", default=2)
parser.add_option('',"--vol",dest='vol',action='store_true',help="Do steps for volumetric analysis only",default=False)
parser.add_option('',"--anat",dest='anat',help="Raw anatomical for freesurfer skull-strip (vol mode)",default='')
(options,args) = parser.parse_args()


#Load existing session if found
if options.sess!='': 
	ifh = open(options.sess,'rb')
	options_in = pickle.load(ifh)
	ifh.close()
	options_in.steps = options.steps
	options_in.stepin = options.stepin
	options_in.noexec = options.noexec
	options_in.jobs = options.jobs
	options = options_in

#Initialize
cmds=["echo \"\nRunning clustering preprocessing script\n\"","export OMP_NUM_THREADS=%s" % options.jobs,"export MKL_NUM_THREADS=%s" % options.jobs]
if options.stepin!='':
	startstep=run_steps_all.index(options.stepin)
	run_steps=run_steps_all[startstep:]
elif options.steps=='all' or options.steps==run_steps_str:
	run_steps=run_steps_all
else:
	run_steps=options.steps.split(',')
if options.vol:
	if 'suma' in run_steps: run_steps[run_steps.index('suma')] = 'fs_skullstrip'
	for ss in ['surfsmooth','surfcomb','suma']: 
		if ss in run_steps: run_steps.remove(ss)

print run_steps

#Do setup
if 'setup' in run_steps:
	subdir = gohome(0)
	cmds.append('mkdir %s' % subdir)
	gohome()
	featsdir = gofeatsdir(0)
	cmds.append("mkdir %s" % featsdir) 
else: gohome()

if 'fs_skullstrip' in run_steps:
	homedir = gohome()
	if options.vol and options.anat!='':
		cmds.append('recon-all -autorecon1 -subjid %s -i %s' % (options.fs.split('/')[-1],options.anat) ) 
	cmds.append('mri_convert %s/mri/brainmask.mgz ./anat_do.nii.gz ; mri_convert %s/mri/T1.mgz ./anat.nii.gz' % (options.fs,options.fs) )
	cmds.append("3dWarp -overwrite -deoblique -prefix anat_do.nii.gz anat_do.nii.gz; 3dAutobox -overwrite -prefix anat_do.nii.gz anat_do.nii.gz" )
	cmds.append('3dWarp -overwrite -deoblique -prefix anat.nii.gz anat.nii.gz')
	cmds.append("ln -s anat_do.nii.gz anat_ns.nii.gz " )

#Make SUMA Spec
if 'suma' in run_steps and not options.vol:
	fsid = getfsid()
	if not os.path.exists("%s/SUMA" % (options.fs)):
		cmds.append("cd %s" % options.fs)
		cmds.append("@SUMA_Make_Spec_FS -sid %s" % (fsid) )
	sumadir = gosumadir(0)
	homedir=gohome()
	cmds.append("if [ ! -e %s ]; then ln -s %s/SUMA %s; fi" % (sumadir,options.fs,sumadir))
	cmds.append("if [ ! -e %s/anat_do.nii.gz ]; then 3dWarp -overwrite -deoblique -prefix %s/anat_do.nii.gz %s/brainmask.nii; 3dAutobox -overwrite -prefix %s/anat_do.nii.gz %s/anat_do.nii.gz; fi" % (homedir,homedir,sumadir, homedir, homedir))
	cmds.append("if [ ! -e %s/anat_ns.nii.gz ]; then ln -s %s/anat_do.nii.gz %s/anat_ns.nii.gz; fi " % (homedir,homedir,homedir))
	cmds.append("if [ ! -e %s/anat.nii.gz ]; then 3dWarp -overwrite -deoblique -prefix %s/anat.nii.gz %s/T1.nii; fi" % (homedir,homedir,sumadir))

#After preparing directories and grabbing anatomicals, set error exit mode
cmds.append('set -e')

#Transfer functional data into home
if 'get_func' in run_steps:
	gohome()
	for ds in parse_mecode(options.func)[0]: 
		cmds.append('ln -s %s/%s* .' % (options.sourcepath,ds))
	
#Run ME-ICA preprocessing
if 'meica_pp' in run_steps:
	homedir = gohome()
	if options.tpattern!='': tps='--tpattern=%s' % options.tpattern
	else: tps=''
	options.fwhm='0mm'
	cmds.append("python27 /home/kundup/code/me-ica/meica.py -d \"%s\" -e %s -f %s -b 0 --t2salign -a anat_ns.nii.gz --keep_int %s --no_skullstrip --pp_only --OVERWRITE" % (options.func,options.TEs,options.fwhm,tps))
	meicadir = gomeicadir(0)
	featsdir = gofeatsdir(0)
	cmds.append("cp %s/motion.1D %s" % (meicadir,homedir) )
	cmds.append("cp %s/zcat_ffd.nii.gz %s/" % (meicadir,featsdir) )

#Run tedanaX
if 'tedanaX' in run_steps:
	gomeicadir()
	cmds.append("python27 /home/kundup/code/me-ica/meica.libs/tedanaX.py -e %s  -d zcat_ffd.nii.gz --sourceTEs=0 --kdaw=10 --rdaw=1 --initcost=pow3 --finalcost=tanh --OC --conv=2.5e-5" % options.TEs)

#Compute non-linear warp on anatomical
if 'fnirt_anat' in run_steps:
	gohome()
	fnirt_anat_script="""
if [ ! -e anat_nlwarp.nii.gz ]; then
3dresample -overwrite -master anat_do.nii.gz -prefix anat_nsf.nii.gz -inset anat_ns.nii.gz -rmode NN
3daxialize -overwrite -orient RPI -prefix anat_nsf.nii.gz anat_nsf.nii.gz
3daxialize -overwrite -orient RPI -prefix anat_do.nii.gz -overwrite anat_do.nii.gz
flirt -ref /usr/local/fsl/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz -in anat_nsf.nii.gz -omat anat_MNI.flirt.mat -out anat_lMNI.nii
fnirt --ref=/usr/local/fsl/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz --in=anat_do.nii.gz --aff=anat_MNI.flirt.mat --iout=anat_nlMNI --cout=anat_nlwarp
fi
"""
	cmds.append(fnirt_anat_script)

#Run ME-ICA FNIRT preprocessing
if 'fnirt_meica_pp' in run_steps:
	homedir = gohome()
	meicadir = gomeicadir()
	nte = len(options.TEs.split(','))
	mnitmpl = "/usr/local/fsl/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz"
	cmds.append("bash %s/meica_fnirt.sh %i %s/anat_do.nii.gz %s %s/anat_nlwarp.nii.gz   " % (execdir,nte,homedir,mnitmpl,homedir))

#Do tedanaX on fnirted me-ica
if 'fnirt_tedanaX' in run_steps:
	gomeicadir()
	cmds.append('cd fnirt')
	cmds.append("python27 /home/kundup/code/me-ica/meica.libs/tedanaX.py -e %s  -d zcat_ffd.nii.gz --sourceTEs=-1 --kdaw=20 --rdaw=2 --initcost=tanh --finalcost=tanh --OC --conv=2.5e-5" % options.TEs)

#Do tedanaX on fnirted me-ica
if 'feats_tedanaX_redo' in run_steps:
	gofeatsdir()
	cmds.append('if [ -d TED ]; then rm -rf TED; fi ')
	cmds.append('mkdir TED; cp meica_mix.1D TED')
	cmds.append("python27 /home/kundup/code/me-ica/meica.libs/tedanaX.py -e %s  -d zcat_ffd.nii.gz --sourceTEs=-1 --kdaw=20 --rdaw=2 --initcost=tanh --finalcost=tanh --OC --conv=2.5e-5 --mix=meica_mix.1D" % options.TEs)
	cmds.append('mv TED/betas_OC.nii TED/betas_hik_OC.nii TED/hik_ts_OC.nii TED/accepted.txt .')
	cmds.append('gzip -f *nii')

#Apply non-linear warp on functional
if 'fnirt_func' in run_steps:
	homedir = gohome()
	meicadir = gomeicadir(0)
	featsdir = gofeatsdir(0)
	cmds.append("3daxialize -prefix %s/ts_OC_rs.nii.gz -orient RPI -overwrite %s/TED/ts_OC.nii " % (meicadir,meicadir))
	cmds.append("3dresample -dxyz 2 2 2 -overwrite -master %s/anat_do.nii.gz -inset %s/ts_OC_rs.nii.gz -prefix %s/ts_OC_rs.nii.gz" % (homedir,meicadir,meicadir))
	cmds.append("fslmaths %s/ts_OC_rs.nii.gz -add 0. %s/ts_OC_rs.nii.gz " % (meicadir,meicadir))
	cmds.append("applywarp --ref=/usr/local/fsl/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz --in=%s/ts_OC_rs.nii.gz --out=%s/ts_OC_mni.nii.gz --warp=%s/anat_nlwarp.nii.gz" % (meicadir,featsdir,homedir))

#Do surface smoothing
if 'surfsmooth' in run_steps and not options.vol:
	sumadir = gosumadir(0)
	meicadir = gomeicadir()
	homedir=gohome(0)
	cmds.append("bash %s/surf_code/surfacesmooth-debug6.sh %s %s %i %s/anat_ns.nii.gz" % (execdir,sumadir,getfsid(),len(options.TEs.split(',') ),homedir))
	cmds.append("cp %s/lh_nct %s/rh_nct %s" % (meicadir,meicadir,homedir))

#Optimally combine surface smoothed data
if 'surfcomb' in run_steps and not options.vol:
	gomeicadir()
	featsdir = gofeatsdir(0)
	cmds.append("python27 /home/kundup/code/me-ica/meica.libs/optcom_surf.py -e %s -p %s/LH_ts.1D e*.lh.1D.dset" % (options.TEs,featsdir))
	cmds.append("python27 /home/kundup/code/me-ica/meica.libs/optcom_surf.py -e %s -p %s/RH_ts.1D e*.rh.1D.dset" % (options.TEs,featsdir))

if 'prep_feats' in run_steps:
	meicadir = gomeicadir(0)	
	featsdir = gofeatsdir(0)
	cmds.append("cp %s/TED/meica_mix.1D %s/" % (meicadir,featsdir) )
	cmds.append("cp %s/TED/accepted.txt %s/" % (meicadir,featsdir) )
	cmds.append("cp %s/TED/comp_table.txt %s/" % (meicadir,featsdir) )
	cmds.append("cp %s/TED/feats_OC2.nii* %s/" % (meicadir,featsdir) )
	cmds.append("cp %s/TED/betas_OC.nii* %s/" % (meicadir,featsdir) )
	cmds.append("cp %s/TED/ts_OC.nii %s/" % (meicadir,featsdir) )
	cmds.append("cp %s/TED/hik_ts_OC.nii* %s/" % (meicadir,featsdir) )
	gofeatsdir()
	cmds.append("gzip -f *nii")

if 'prep_fnirt_feats' in run_steps:
	meicadir = gomeicadir(0)	
	featsdir = gofeatsdir(0)
	cmds.append("cp %s/fnirt/zcat_ffd.nii.gz %s/" % (meicadir,featsdir) )
	cmds.append("cp %s/fnirt/TED/meica_mix.1D %s/" % (meicadir,featsdir) )
	cmds.append("cp %s/fnirt/TED/accepted.txt %s/" % (meicadir,featsdir) )
	cmds.append("cp %s/fnirt/TED/comp_table.txt %s/" % (meicadir,featsdir) )
	cmds.append("cp %s/fnirt/TED/feats_OC2.nii* %s/" % (meicadir,featsdir) )
	cmds.append("cp %s/fnirt/TED/betas_OC.nii* %s/" % (meicadir,featsdir) )
	cmds.append("cp %s/fnirt/TED/ts_OC.nii %s/" % (meicadir,featsdir) )
	cmds.append("cp %s/fnirt/TED/hik_ts_OC.nii* %s/" % (meicadir,featsdir) )
	gofeatsdir()
	cmds.append("gzip -f *nii")

if 'save_space' in run_steps:
	homedir = gohome(0)
	meicadir = gomeicadir(0)
	cmds.append("cp %s/motion.1D %s" % (meicadir,homedir) )
	cmds.append("rm -rf %s" % meicadir)

#Prepare features based on high dimensional ME-ICA
if 'meica_feats' in run_steps:
	mfeats = gomfeatsdir(0)
	featsdir = gofeatsdir(0)
	cmds.append("if [ ! -e %s ]; then mkdir %s; fi" % (mfeats,mfeats)) 
	cmds.append("cd %s" % (mfeats))
	cmds.append("3dresample -overwrite -inset %s/ts_OC_mni.nii.gz -prefix %s/ts_OC_mni.nii.gz -master %s/deepstruct_nlmni.nii.gz" % (featsdir,featsdir,execdir)) #should be moved to fnirt_func
	if options.vol:
		cmds.append("python27  %s/surf_code/svfeats.py --vol=%s/ts_OC_mni.nii.gz --vol1Dmask %s/deepstruct_nlmni.nii.gz --mix=%s/meica_mix.1D  --mixacc=%s/accepted.txt  --noempty=0" % (execdir,featsdir,execdir,featsdir,featsdir) )
	else:
		cmds.append("python27  %s/surf_code/svfeats.py --vol=%s/ts_OC_mni.nii.gz --vol1Dmask %s/deepstruct_nlmni.nii.gz --lsurf=%s/LH_ts.1D --rsurf=%s/RH_ts.1D --mix=%s/meica_mix.1D  --mixacc=%s/accepted.txt  --noempty=0" % (execdir,featsdir,execdir,featsdir,featsdir,featsdir,featsdir))
		cmds.append("mv ts_OC_mni_masked_2c.1D.dset SC_feats_norm.1D")
		cmds.append("mv LH_ts_2c.1D.dset LH_feats_norm.1D")
		cmds.append("mv RH_ts_2c.1D.dset RH_feats_norm.1D")
	cmds.append("gzip -f *nii")

if 'ts_reg' in run_steps:
	homedir = gohome(0)
	featsdir = gofeatsdir()
	cmds.append("1d_tool.py -overwrite -derivative -infile %s/motion.1D -write %s/motion_deriv.1D" % (homedir,featsdir))
	cmds.append("1dcat %s/motion.1D motion_deriv.1D > motion_regr.1D" % (homedir))
	cmds.append("3dBandpass -overwrite -despike -prefix ts_OC_mni_bp.nii.gz -ort motion_regr.1D 0.01 0.1  ts_OC_mni.nii.gz ")

#Do sparse PCA
if 'spca' in run_steps:
	sfeats = gosfeatsdir(0)
	cmds.append("if [ ! -e %s ]; then mkdir %s; fi" % (sfeats,sfeats))
	cmds.append("cd %s" % (sfeats))
	#cmds.append("""nc=`python27 -c "import numpy; ct=numpy.loadtxt('../comp_table.txt'); print len(ct[ct[:,1]>55.5,0])"`""")
	cmds.append("""nc=`python27 -c "print len(open('../accepted.txt').readlines()[0].split(','))"`""")
	cmds.append("python27 %s/surf_code/spca.py --data=../hik_ts_OC.nii.gz --means=../ts_OC.nii.gz --nc=${nc} --detrend A --jobs %s" % (execdir,options.jobs))

#Generate sparse PCA features
if 'spca_feats' in run_steps:
	sfeats = gosfeatsdir(0)
	featsdir = gofeatsdir(0)
	cmds.append("if [ ! -e %s ]; then mkdir %s; fi" % (sfeats,sfeats))
	sfeats = gosfeatsdir()
	if options.vol:
		cmds.append("python27 %s/surf_code/svfeats.py --vol=%s/ts_OC_mni.nii.gz --vol1Dmask %s/deepstruct_nlmni.nii.gz --mix=%s/meica_mix.1D  --mixacc=%s/accepted.txt --dmix=hik_pca_sV_cz.1D --noempty=0" % (execdir,featsdir,execdir,featsdir,featsdir))
	else:	
		cmds.append("python27 %s/surf_code/svfeats.py --vol=%s/ts_OC_mni.nii.gz --vol1Dmask %s/deepstruct_nlmni.nii.gz --lsurf=%s/LH_ts.1D --rsurf=%s/RH_ts.1D --mix=%s/meica_mix.1D  --mixacc=%s/accepted.txt --dmix=hik_pca_sV_cz.1D --noempty=0" % (execdir,featsdir,execdir,featsdir,featsdir,featsdir,featsdir))
		cmds.append("mv ts_OC_mni_masked_2c.1D.dset SC_feats_norm.1D")
		cmds.append("mv LH_ts_2c.1D.dset LH_feats_norm.1D")
		cmds.append("mv RH_ts_2c.1D.dset RH_feats_norm.1D")
	cmds.append("gzip -f *nii")

if 'done' in run_steps:
	featsdir=gofeatsdir()
	cmds.append('touch done')


#Save and run script
cmdout = '\n'.join(cmds)
sessfn = getsessfn()
ofh = open(sessfn,'wb')
pickle.dump(options,ofh)
ofh.close()
scriptfn= getscriptfn()
ofh = open(scriptfn,'w')
ofh.write(cmdout)
ofh.close()
if not options.noexec: os.system('bash %s' % (scriptfn))


