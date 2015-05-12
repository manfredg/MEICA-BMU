

import sys
import commands
from re import split as resplit
import re
from os import system,getcwd,mkdir,chdir,popen
import os.path
from string import rstrip,split
from optparse import OptionParser,OptionGroup
import cPickle as pickle
#import ipdb

#Initiate session
class Sess:
	def __init__(self):
		return None
S = Sess

run_steps_all = ['setup','mepre','suma','sumavol','svmerge','svtedana','export']
run_steps_str = ','.join(run_steps_all)

from optparse import OptionParser
parser=OptionParser()
parser.add_option('-e',"",dest='tes',help="ex: -e 14.5,38.5,62.5  Echo times (in ms)",default='')
parser.add_option('-d',"",dest='dsinputs',help="ex: -d RESTe1.nii.gz,RESTe2.nii.gz,RESTe3.nii.gz",default='')
parser.add_option('-a',"",dest='anat',help="ex: -a mprage.nii.gz  Anatomical dataset (optional)",default='')
parser.add_option('-b',"",dest='basetime',help="ex: -b 10s OR -b 10v  Time to steady-state equilibration in seconds(s) or volumes(v). Default 0. ",default='0')
parser.add_option('-f',"",dest='fsid',help="ex: -f sub_12345   FreeSurfer subject ID ",default=None)
parser.add_option('-p',"",dest='prefix',help="ex: -p sub_12345   Prefix for directory of output from this pipeline ",default=None)
parser.add_option('-t',"",dest='template',help="MNI template space with surface and volume",default=False)
extopts=OptionGroup(parser,"Extra options")
extopts.add_option('',"--fwhm",dest='fwhm',help="Desired volumetric spatial FWHM smoothing, converted to surface equivalent ",default=None)
extopts.add_option('',"--ld",dest='ld',help="Linear mesh density, default 60 ",default=None)
extopts.add_option('',"--mex",dest='mex',help="Additional options sent to ME-ICA",default='')
extopts.add_option('',"--align_args",dest='align_args',help="Additional arguments to anatomical-functional co-registration routine",default='')
extopts.add_option('',"--ted_args",dest='ted_args',help="Additional arguments to TE-dependence analysis routine",default='')
parser.add_option_group(extopts)
runopts=OptionGroup(parser,"Run options")
runopts.add_option('',"--sourcepath",dest='sourcepath',help="Source directory for functionals",default=os.getcwd())
runopts.add_option('',"--jobs",dest='jobs',help="Number of CPUs to use", default=2)
runopts.add_option('',"--steps",dest='steps',help="Do procedures, default all: %s " % run_steps_str,default=run_steps_str)
runopts.add_option('',"--noexec",dest='noexec',action='store_true',help="Just write script, dont run", default=False)
parser.add_option_group(runopts)
(options,args) = parser.parse_args()

run_steps = options.steps.split(',')

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


def logcomment(comment,level=3): 
	majmark='\n'
	leading='--------'
	if level==3: majmark=''
	if level==1: 
		leading="+* "
		cmds.append("""echo "\n++++++++++++++++++++++++" """)
		majmark=''
	cmds.append("""%secho %s"%s" """ % (majmark,leading,comment))

def dirmakeenter(dir,enter=1):
	cmds.append("if [ ! -d %s ]; then mkdir %s; fi" % (dir,dir))
	if enter: cmds.append("cd %s" % (dir))

def gohome(go=1): 
	subdir = '%s' % (S.bdir)
	if go: cmds.append("cd %s"  % (subdir))
	return subdir

def gomeicadir(go=1):
	meicadir="%s/meica.%s" % (gomeicaPdir(0),S.setname)
	if go: cmds.append("cd %s" % (meicadir))
	return meicadir

def gomeicaPdir(go=1):
	S.meicaPdir = os.path.join(gohome(0),'meicaP.%s' % options.prefix)
	if go==1: cmds.append("cd %s" % S.meicaPdir)
	elif go==3: dirmakeenter(S.meicaPdir)
	return S.meicaPdir

def gosumadir(go=1):
	homedir = gomeicaPdir(0)
	sumadir = "%s/%s_SUMA" % (homedir,options.fsid)
	if go: cmds.append("cd %s" % sumadir)
	return sumadir

def gofsdir(go=1):
	fsdir = os.path.join(os.environ['SUBJECTS_DIR'],options.fsid)
	if go: cmds.append("cd %s" % fsdir)
	return fsdir

def gotmpldir(go=1):
	tmpldir = options.template
	if tmpldir[0] != '/': tmpldir = os.path.join(os.path.abspath(os.path.curdir),tmpldir)
	if go: cmds.append("cd %s" % tmpldir)
	return tmpldir

def gospmdir(go=1):
	spmdir= '%s/%s_inspace' % (gomeicaPdir(0),S.setname)
	if go: cmds.append('cd %s' % spmdir)
	return spmdir

def gontmdir(go=1):
	ntmdir= '%s/%s_innative' % (gomeicaPdir(0),S.setname)
	if go: cmds.append('cd %s' % ntmdir)
	return ntmdir

def godecodir(go=1):
	decodir = '%s/%s_deco' % (gomeicaPdir(0),S.setname)
	if go: cmds.append('cd %s' % decodir)
	return decodir

def parseade():
	#Parse locations
	S.bdir = rstrip(popen('pwd').readlines()[0])
	S.codedir=os.path.dirname(os.path.abspath(os.path.expanduser(sys.argv[0])))
	if os.path.exists('%s/meica.libs' % (S.codedir)): 
		S.tedanapath = 'meica.libs/tedana.py'
	S.fspath = gofsdir(0)
	#Parse TEs
	S.tes=split(options.tes,',')
	#Parse functional dataset specification
	S.outprefix=options.prefix
	S.isf = dssuffix(options.dsinputs)
	if '[' in options.dsinputs:
		dsinputs=dsprefix(options.dsinputs)
		S.prefix=resplit(r'[\[\],]',dsinputs)[0]
		S.dsids=sorted(resplit(r'[\[\],]',dsinputs)[1:-1])
		S.trailing=resplit(r'[\]+]',dsinputs)[-1]
		S.setname=os.path.basename(S.prefix+''.join(S.dsids)+S.trailing)
		S.datasets = [S.prefix + dd + S.trailing+S.isf for dd in S.dsids]
	else:
		datasets_in = options.dsinputs.split(',')
		S.dsids = [str(ii) for ii in sorted(range(1,len(datasets_in)+1))]
		S.datasets = datasets_in
		S.prefix = dsprefix(datasets_in[0])
		S.trailing=''
		S.setname=S.prefix
	#Parse anatomical name
	S.anatsource = os.path.dirname(options.anat)
	S.anatfn = os.path.basename(options.anat)
	if S.anatsource=='': S.anatsource = os.path.abspath(os.path.curdir)
	S.anatprefix=dsprefix(os.path.basename(options.anat))	
	S.sourceanatpath = os.path.join(S.anatsource,S.anatprefix+dssuffix(options.anat))
	S.nsmprage="%s_ns.nii.gz" % dsprefix(S.anatprefix)
	#Parse functional files
	S.funcsource = os.path.dirname(options.dsinputs)
	S.funcfn = os.path.basename(options.dsinputs)
	if S.funcsource == '': S.funcsource = os.path.abspath(os.path.curdir)
	S.funcpaths = [os.path.join(S.funcsource,os.path.basename(ff)) for ff in S.datasets]
	#Parse smoothing parameters
	S.volsmopt = '--fwhm=%s' % (options.fwhm) if options.fwhm!=None else ''
	S.surfsmopt = '--FWHM %f' %  (float(S.volsmo) * sqrt(2./3)) if options.fwhm!=None else ''
	#Parse template name
	S.tmplname = os.path.basename(options.template.rstrip('/')).split('.mpt')[0]
	
	print S.anatsource, S.funcsource, S.anatfn, S.funcfn

parseade()

if __name__=='__main__':

	cmds = []
	cmds.append('#!/bin/bash' )
	cmds.append('#%s' % ' '.join(sys.argv) )
	cmds.append('set -e')

	if 'setup' in run_steps:
		gomeicaPdir(3)
		cmds.append("ln -sf %s* %s/" % ( S.sourceanatpath,gomeicaPdir(0)))
		for ff in S.funcpaths:
			cmds.append("ln -sf %s* %s/" % (ff,gomeicaPdir(0)))

	
	#Gradwarp and fieldmap steps here

	if 'mepre' in run_steps:
		logcomment("Performing ME-ICA preprocessing steps basic functional corrections and anatomical-functional coregistration",level=1)
		gomeicaPdir(1)
		infargs='--OVERWRITE --qwarp --fres 2.5 --space=%s/stdvol+tlrc --keep_int --pp_only %s' % (options.template,S.volsmopt) #!! --qwarp should be here but ommitted to save time 
		cmds.append('%s %s/meica.py -e %s -d \"%s\" -a %s %s %s %s %s' % (sys.executable,S.codedir,options.tes,S.funcfn,S.anatfn,infargs,options.align_args,options.ted_args,options.mex))
		dirmakeenter('%s_inspace' % S.setname,0)		
		cmds.append('cp %s/*_in.nii.gz %s_inspace/' % (gomeicadir(0),S.setname))
		cmds.append('cp %s/eBvrmask.nii.gz %s_inspace/' % (gomeicadir(0),S.setname))
		infargs='--fres 2.5 --RESUME --pp_only %s' % S.volsmopt
		cmds.append('%s %s/meica.py -e %s -d \"%s\" -a %s %s %s %s %s' % (sys.executable,S.codedir,options.tes,S.funcfn,S.anatfn,infargs,options.align_args,options.ted_args,options.mex))
		dirmakeenter('%s_innative' % S.setname,0)
		cmds.append('cp %s/*_in.nii.gz %s_innative/' % (gomeicadir(0),S.setname))
		cmds.append('cp %s/eBvrmask.nii.gz %s_innative/' % (gomeicadir(0),S.setname))

	if 'suma' in run_steps:
		logcomment("Doing SUMA steps",level=1)
		fsdir = gofsdir(1)
		#Create SUMA session if it doesn't exist
		cmds.append("if [ ! -d %s/SUMA ]; then cd %s; @SUMA_Make_Spec_FS -sid %s; fi" % (S.fspath, S.fspath, options.fsid) )
		cmds.append("if [ ! -d %s ]; then ln -s %s/SUMA %s; fi" % (gosumadir(0),S.fspath,gosumadir(0)))
		#Warp to me-ica's skullstripped anat
		gomeicaPdir()
		cmds.append("for ff in `ls %s/*rank*nii`; do 3daxialize -overwrite -prefix $ff $ff; done" % (gosumadir(0))  )
		cmds.append("@SUMA_AlignToExperiment -exp_anat %s/%s -surf_anat %s/%s_SurfVol+orig -atlas_followers -followers_interp NN -strip_skull surf_anat -overwrite_resp O -surf_anat_followers %s/aseg_rank.nii %s/aparc.a2009s+aseg_rank.nii" % (gomeicaPdir(0),S.nsmprage,gosumadir(0),options.fsid,gosumadir(0),gosumadir(0)) )		
		gosumadir(1)
		#MapIcosahedron for subject???
		cmds.append("hemis=\"lh rh\"; for hemi in $hemis; do MapIcosahedron -overwrite -spec std.141.%s_${hemi}.spec -ld 200 -prefix ld200_ -morph std.141.${hemi}.sphere.asc; done" % (options.fsid))
		cmds.append('cd %s/%s_innative' % (gomeicaPdir(0),S.setname))
		cmds.append("hemis=\"lh rh\"; for hemi in $hemis; do")
		for dii in S.dsids:
			cmds.append('cd %s/%s_innative' % (gomeicaPdir(0),S.setname))
			cmds.append("3dVol2Surf -overwrite -spec %s/std.141.%s_${hemi}.spec -surf_A smoothwm -surf_B pial -sv %s/%s_SurfVol_ns_Alnd_Exp+orig -grid_parent e%s_in.nii.gz -map_func ave -f_steps 15 -f_index nodes -outcols_NSD_format -out_niml ./e%s_in_${hemi}.std_141.niml.dset" % (gosumadir(0),options.fsid,gomeicaPdir(0),options.fsid,dii,dii) )
			cmds.append("ConvertDset -o_niml_bi -input ./e%s_in_${hemi}.std_141.niml.dset -pad_to_node ld141 -overwrite -prefix ./e%s_in_${hemi}.std_141.niml.dset" % (dii,dii) )	
			#if S.surfsmopt!='': cmds.append("SurfSmooth -overwrite  -input e${echon}_vr_bp.${hemi}.niml.dset -target_fwhm 12mm -spec ${sumadir}/ld60_${sub}_${hemi}.spec -surf_A pial -surf_B smoothwm -met HEAT_07 -add_index -output ./e${echon}_vr_bp_sm.${hemi}.1D.dset ")
		cmds.append('done')

	if 'sumavol' in run_steps:
		logcomment("Doing SUMA surface to volume mapping",level=1)
		dset_list = [ '%s/%s_innative/e%s_in_${hemi}.std_141.niml.dset' % (gomeicaPdir(0),S.setname,dii) for dii in S.dsids ]
		S.us_dset_list = [ '%s/%s_innative/ld200_e%s_in_${hemi}.std_141.niml.dset' % (gomeicaPdir(0),S.setname,dii) for dii in S.dsids ]	
		S.dsetmap_list = '-NN_dset_map ' + ' -NN_dset_map '.join(dset_list)
		gosumadir(1)
		cmds.append("hemis=\"lh rh\"; for hemi in $hemis; do MapIcosahedron -overwrite -spec std.141.%s_${hemi}.spec -ld 200 -prefix ld200_ -morph std.141.${hemi}.sphere.asc %s; done" % (options.fsid,S.dsetmap_list))
		cmds.append('cd %s/%s_innative' % (gomeicaPdir(0),S.setname))
		cmds.append("hemis=\"lh rh\"; for hemi in $hemis; do ln -fs %s %s/%s; done" % (' '.join(S.us_dset_list),gomeicaPdir(0),'%s_inspace' % S.setname)   )
		S.us_dset_list = [ '%s/%s_innative/ld200_e%s_in_${hemi}.std_141.niml.dset' % (gomeicaPdir(0),S.setname,dii) for dii in S.dsids ]	
		cmds.append('cd %s/%s_inspace' % (gomeicaPdir(0),S.setname))
		cmds.append('3dresample -overwrite -inset %s/abtemplate.nii.gz -dxyz 2.5 2.5 2.5 -prefix abtemplate_epi.nii.gz -rmode Li ' % gomeicadir(0))
		for dset in S.us_dset_list:		
			cmds.append("hemis=\"lh rh\"; for hemi in $hemis; do 3dSurf2Vol -spec %s/ld200_std.141.%s_${hemi}.spec -surf_A pial -surf_B smoothwm -grid_parent abtemplate_epi.nii.gz -sv %s/%s_SurfVol_Alnd_Exp+orig.HEAD -datum float -f_steps 20 -map_func mode -sdata %s -overwrite -prefix %s.nii.gz; done" % (gotmpldir(0),S.tmplname,gotmpldir(0),S.tmplname,dset,dsprefix(os.path.basename(dset)) ))

	if 'svmerge' in run_steps:
		logcomment("Merging surface and volume warped datasets",level=1)
		gospmdir(1)
		cmds.append("3dNwarpApply -overwrite -nwarp '%s/%s_ns_atnl_WARP.nii.gz %s/%s_ns2at.aff12.1D' -master %s/abtemplate_epi.nii.gz -dxyz 2.5 -source %s/aseg_rank_Alnd_Exp+orig -interp NN -prefix %s/sub_aseg_rank_epi.nii.gz" % (gomeicaPdir(0),S.anatprefix,gomeicaPdir(0),S.anatprefix,gospmdir(0),gomeicaPdir(0),gospmdir(0)))
		cmds.append("3dresample -overwrite -inset  %s/aseg_rank_Alnd_Exp+orig -master %s/abtemplate_epi.nii.gz -rmode NN -prefix %s/tmpl_aseg_rank_epi.nii.gz"% (gotmpldir(0),gospmdir(0),gospmdir(0)) )
		cmds.append("3dresample -overwrite -inset  %s/deepstruct.nii.gz -master %s/abtemplate_epi.nii.gz -rmode NN -prefix %s/deepstruct_epi.nii.gz"% (gotmpldir(0),gospmdir(0),gospmdir(0)) )
		merges_zcat = []		
		for dsid in S.dsids:
			cmds.append("%s %s/meica.libs/surfvolmerge.py ld200_e%s_in_lh.std_141.niml.dset.nii.gz ld200_e%s_in_rh.std_141.niml.dset.nii.gz e%s_in.nii.gz deepstruct_epi.nii.gz tmpl_aseg_rank_epi.nii.gz sub_aseg_rank_epi.nii.gz e%s_svmerge.nii.gz" % (sys.executable,S.codedir,dsid,dsid,dsid,dsid))
			merges_zcat.append('e%s_svmerge.nii.gz' % dsid)
		cmds.append("3dZcat -prefix zcat_ffd.nii.gz %s" % ' '.join(merges_zcat))
		dirmakeenter(godecodir(0),1)
		cmds.append("mv %s/zcat_ffd.nii.gz %s" % (gospmdir(0),godecodir(0)) )

	if 'svtedana' in run_steps:
		logcomment("Executing decomposition steps",level=1)
		godecodir(1)
		cmds.append("%s %s/meica.libs/tedana.py -e %s  -d %s/zcat_ffd.nii.gz --sourceTEs=-1 --kdaw=10 --rdaw=1 --initcost=tanh --finalcost=tanh --conv=2.5e-5" % (sys.executable,S.codedir,options.tes,godecodir(0)))
		cmds.append("3dTstat -prefix ./mean.nii.gz TED/dn_ts_OC.nii")
		cmds.append("3dcalc -a mean.nii.gz -b TED/hik_ts_OC.nii -expr 'a+b' -prefix ./hik_ts_OCm.nii.gz")
		cmds.append("%s %s/meica.libs/godec.py -p 3 -t 3 -r 1 -d hik_ts_OCm.nii.gz" % (sys.executable,S.codedir))

	if 'export' in run_steps:
		logcomment("Copying data to start directory",level=1)
		godecodir(1)
		cmds.append("cp %s/noise_nvnr1k2p3t3.nii %s/%s_hiks.nii" %  (godecodir(0),S.bdir,options.prefix) )
		cmds.append("cp %s/TED/betas_hik_OC.nii %s/%s_mefc.nii" %  (godecodir(0),S.bdir,options.prefix) )
		cmds.append("cp %s/TED/ts_OC.nii %s/%s_tsoc.nii" %  (godecodir(0),S.bdir,options.prefix) )
		cmds.append("nifti_tool -mod_hdr -mod_field sform_code 2 -mod_field qform_code 2 -infiles %s/%s_tsoc.nii -overwrite" % (S.bdir,options.prefix) )
		cmds.append("nifti_tool -mod_hdr -mod_field sform_code 2 -mod_field qform_code 2 -infiles %s/%s_hiks.nii -overwrite" % (S.bdir,options.prefix) )
		cmds.append("nifti_tool -mod_hdr -mod_field sform_code 2 -mod_field qform_code 2 -infiles %s/%s_mefc.nii -overwrite" % (S.bdir,options.prefix) )
		cmds.append("gzip -f %s/%s_hiks.nii" %  (S.bdir,options.prefix) )
		cmds.append("gzip -f %s/%s_mefc.nii" %  (S.bdir,options.prefix) )
		cmds.append("gzip -f %s/%s_tsoc.nii" %  (S.bdir,options.prefix) )
		cmds.append("cp %s/motion.1D %s/%s_motion.1D" %  (gomeicadir(0),S.bdir,options.prefix) )
		cmds.append("1d_tool.py -overwrite -derivative -infile %s/motion.1D -write %s/motion_dt.1D" % (gomeicadir(0),gomeicadir(0)))
		cmds.append("1dcat %s/motion.1D %s/motion_dt.1D > %s/%s_moreg.1D" % (gomeicadir(0),gomeicadir(0),S.bdir,options.prefix))
		
	
	outscr='\n'.join(cmds)
	scrfile = '_meicaP_%s.sh' % options.prefix
	ofh = open(scrfile,'w')
	ofh.write(outscr)
	ofh.close()
	
	if not options.noexec:
		os.system("bash %s/%s"  % (os.getcwd(),scrfile) )








