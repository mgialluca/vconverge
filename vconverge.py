import os
import sys
from astropy.io import ascii
#import corner
import numpy as np
import re
import shutil
import bigplanet as bp
import json
import copy
from scipy.stats import ks_2samp  
from scipy.special import kl_div
#import matplotlib.pyplot as plt

# Parameter sweep until user defined params converge

def extract_info_vcnv(vcnvFile):
	vcnv = open(vcnvFile, 'r')
	lines = vcnv.readlines()
	params_to_conv = []
	hold = []
	sConvergenceMethod = 'KL'
	for i in range(len(lines)):
		if lines[i].split() == []:
			pass
		elif lines[i].split()[0] == 'sVspaceFile':
			vspFi = lines[i].split()[1]
		elif lines[i].split()[0] == 'iStepSize':
			StepSize = int(lines[i].split()[1])
		elif lines[i].split()[0] == 'iMaxSteps':
			MaxSteps = int(lines[i].split()[1])
		elif lines[i].split()[0] == 'sConvergenceMethod':
			ConvMethod = lines[i].split()[1]
		elif lines[i].split()[0] == 'fConvergenceCondition':
			ConvCondit = float(lines[i].split()[1])
		elif lines[i].split()[0] == 'iNumberOfConvergences':
			ConvNum = int(lines[i].split()[1])
		elif lines[i].split()[0] == 'sObjectFile':
			if hold != [] and len(hold) > 1:
				params_to_conv.append(hold)
			hold = []
			hold.append(lines[i].split()[1])
		elif lines[i].split()[0] == 'saConverge':
			if len(lines[i].split()) != 3:
				raise IOError('Regarding line: %s Please be sure to only specify final/initial after saConverge, then the name of the param to converge, vconverge\'s saConverge option cannot accept any more or less than these inputs' % lines[i])
			if lines[i].split()[1] != 'final' and lines[i].split()[1] != 'initial':
				raise IOError('Formatting of line: %s Is incorrect, please use: saConverge [initial/final] [name of param to converge]' % lines[i])
			if hold == []:
				raise IOError('Please specify the body you would like %s to converge for' % lines[i].split()[2])
			param = lines[i].split()[2]+','+lines[i].split()[1]
			hold.append(param)
	if hold != [] and len(hold) > 1:
		params_to_conv.append(hold)

	return vspFi, StepSize, MaxSteps, ConvMethod, ConvCondit, ConvNum, params_to_conv

def extract_info_vsp(vspFile): # Extracts relevant info from vspace.in file for vconverge including srcfolder, destfolder, and predefined prior file info (if applicable)
	vspog = open(vspFile, 'r')
	linesog = vspog.readlines()
	PrimeFi = 'vpl.in'
	for i in range(len(linesog)):
		if linesog[i].split() == []:
			pass
		#elif linesog[i].split()[0] == 'srcfolder': # Old Syntax, before backwards compatibility break
		elif linesog[i].split()[0] == 'sSrcFolder': # Code will want to know source folder
			if linesog[i].split()[1] == '.': # Current working dir is srcfolder, get full path for vspacetmp srcfolder
				source_fold = os.getcwd()
			else:
				source_fold = linesog[i].split()[1]
		#elif linesog[i].split()[0] == 'destfolder': # Old Syntax, before backwards compatibility break
		elif linesog[i].split()[0] == 'sDestFolder': # Code will want to know dest folder
			dest_fold = linesog[i].split()[1]
		#elif linesog[i].split()[0] == 'samplemode': # Old Syntax, before backwards compatibility break
		elif linesog[i].split()[0] == 'sSampleMode': # If sample mode isn't random lets stop now
			modename = linesog[i].split()[1]
			if modename.startswith('r') or modename.startswith('R'):
				mode = 'random'
			else:
				raise IOError('vconverge is currently only compatible with random mode')
		#elif linesog[i].split()[0] == 'randsize': # Old Syntax, before backwards compatibility break
		elif linesog[i].split()[0] == 'iRandSize': 
			initial_sim_size = int(linesog[i].split()[1])
		#elif linesog[i].split()[0] == 'trialname': # Old Syntax, before backwards compatibility break
		elif linesog[i].split()[0] == 'sTrialName': 
			triname = linesog[i].split()[1]
		#elif linesog[i].split()[0] == 'file' or linesog[i].split()[0] == 'sPrimaryFile': # Old Syntax, before backwards compatibility break
		elif linesog[i].split()[0] == 'sBodyFile' or linesog[i].split()[0] == 'sPrimaryFile': # Set the current file in case we need to append var names for predefined prior mode
			curr_fi = linesog[i].split()[1]
			curr_fi = curr_fi.split('.')[0]
			if linesog[i].split()[0] == 'sPrimaryFile':
				PrimeFi = linesog[i].split()[1]

	return source_fold, dest_fold, triname, PrimeFi, initial_sim_size#, prior_files, prior_vars, prior_vars_cols
			

def create_tmp_vspin(vspFile, RunIndex, stepsize): # Creates a temporary vspace.in file to run for steps subsequent to original run
	vspog = open(vspFile, 'r') # Read in original (og) vspace file provided as input by user
	linesog = vspog.readlines()
	vspog.close()
	vsptmp = open('vconverge_tmp/vspace_tmp.in', 'w') # Create new temp vspace file in vconverge tmp directory
	for i in range(len(linesog)):
		if linesog[i].split() == []:
			pass
		#elif linesog[i].split()[0] == 'destfolder': # Old Syntax, before backwards compatibility break
		elif linesog[i].split()[0] == 'sDestFolder': # Destination folder needs to be in vconverge_tmp, unique runs will be copied into OG destfolder
			vsptmp.write('destfolder vconverge_tmp/Step_'+str(RunIndex)+'\n')
		#elif linesog[i].split()[0] == 'trialname': # Old Syntax, before backwards compatibility break
		elif linesog[i].split()[0] == 'sTrialName': # Trial name must be unique to not overwrite any data during copy
			vsptmp.write('trialname Step'+str(RunIndex)+'_'+linesog[i].split()[1]+'\n')
		#elif linesog[i].split()[0] == 'randsize': # Old Syntax, before backwards compatibility break
		elif linesog[i].split()[0] == 'iRandSize': 
			vsptmp.write('randsize '+str(stepsize)+'\n')
		elif re.search("\[", linesog[i]) != None: # look for predefined prior files, change to tmp name
			spl = re.split("[\[\]]", linesog[i])
			var = spl[0].strip() # Current variable
			values = spl[1].split(",")
			for j in range(len(values)):
				values[j] = values[j].strip()
			if values[2][0] == 'p':
				priorhold = values[0].split('/')
				tmpprior = priorhold[len(priorhold)-1]
				tmpprior = 'tmp_'+tmpprior
				vsptmp.write(var+' [vconverge_tmp/'+tmpprior+', '+values[1]+', '+values[2]+', '+values[3]+'] '+spl[2].strip()+'\n')
			else:
				vsptmp.write(linesog[i]+'\n')
		else:
			vsptmp.write(linesog[i]+'\n')
	vsptmp.close()


def create_tmp_prior_files(RunIndex, og_triname, dst_fold): # Create new prior files deleting previously used priors
	if RunIndex == 1:
		prev = open(os.path.join(dst_fold, og_triname+'PriorIndicies.json'), 'r')
	else:
		prev = open(os.path.join(dst_fold, 'Step_'+str(RunIndex-1)+'/Step'+str(RunIndex-1)+'_'+og_triname+'PriorIndicies.json'), 'r')
	prev = json.load(prev)
	prev = json.loads(prev)
	for i in prev:
		finame = i.split('/')[len(i.split('/')) - 1]
		extension = finame[len(finame) - 3]+finame[len(finame) - 2]+finame[len(finame) - 1]
		if extension == 'npy':
			priorfi = np.load(i)
			newpriorfi = np.delete(priorfi, prev[i], axis=0)
			if RunIndex == 1:
				np.save('vconverge_tmp/tmp_'+finame, newpriorfi)
			else:
				np.save('vconverge_tmp/'+finame, newpriorfi)
		elif extension == 'txt' or extension == 'dat':
			priorfi = ascii.read(i)
			del(priorfi[np.array(prev[i])])
			if RunIndex == 1:
				ascii.write(priorfi, 'vconverge_tmp/tmp_'+finame, format='fixed_width', delimiter=' ', overwrite=True)
			else:
				ascii.write(priorfi, 'vconverge_tmp/'+finame, format='fixed_width', delimiter=' ', overwrite=True)

def vconverge(vcnvFile):
	print('in the function')
	# make the temporary directory
	if os.path.exists('vconverge_tmp'):
		shutil.rmtree('vconverge_tmp')
	os.mkdir('vconverge_tmp')

	# extract required info from the vconverge.in file and the vspace.in file respectively
	vspFile, StepSize, MaxSteps, ConvMethod, ConvCondit, ConvNum, params_to_conv = extract_info_vcnv(vcnvFile)
#	src_fold, dst_fold, og_triname, primeFi, pfiles, pvars, pvarscols = extract_info_vsp(vspFile)
	src_fold, dst_fold, og_triname, primeFi, initialsims = extract_info_vsp(vspFile)

	# Check that vconverge can find the files associated with the bodies that have requested converging parameters
	# These files should be among the files copied by vspace and should be found in the srcfolder defined in the vspace.in file
	for i in params_to_conv:
		if not os.path.isfile(os.path.join(src_fold, i[0])):
			raise IOError('%s does not exist. Please specify body/primary file source directories in the vspace.in file. Give only the name of the file to option sObjectFile for vconverge.' % os.path.join(src_fold, i[0]))

	# Extract the sName of the bodies for each requested converging parameter then delete the list entry containing the name of the file
	# append sName (body name that will be used by VPLanet) to each parameter associated with it to build up a string of bodyname_parametername_finalorinitial
	for i in params_to_conv:
		curr_fi = open(os.path.join(src_fold, i[0]), 'r')
		currlines = curr_fi.readlines()
		curr_fi.close()
		for k in range(len(currlines)):
			if currlines[k].split() == []:
				pass
			if currlines[k].split()[0] == 'sName':
				body = currlines[k].split()[1]
				break
		for k in range(len(i)):
			if k != 0:
				i[k] = body+','+i[k]
		i.pop(0)
	
	# turn list of lists into one single list for simplicity's sake
	hold = []
	for i in params_to_conv:
		for k in range(len(i)):
			hold.append(i[k])
	params_to_conv = hold
	del(hold)

	# Get the future name of the log files from the primary file (sSystemName + .log) 
	if not os.path.isfile(os.path.join(src_fold, primeFi)):
		raise IOError('vconverge cannot find %s. Please specify the primary file and source directory in the vspace.in file as vspace would expect it, this is also how vconverge will expect it.' % os.path.join(src_fold, primeFi))
	primefi_open = open(os.path.join(src_fold, primeFi))
	primelines = primefi_open.readlines()
	primefi_open.close()
	for i in range(len(primelines)):
		if primelines[i].split() == []:
			pass
		elif primelines[i].split()[0] == 'sSystemName':
			vplanet_logfile = primelines[i].split()[1]+'.log'
			break

	# Create a dictionary to hold the values to converge from every simulation
	converge_dict = {}
	for i in params_to_conv:
		converge_dict[i] = []

	#Run Vspace on OG
	os.system('python -m vspace '+str(vspFile))
	#os.system('multi-planet '+str(vspFile))
	os.system('python -m multiplanet '+str(vspFile))
	#Run Multi-planet on OG
	RunIndex = 1
	create_tmp_vspin(vspFile, RunIndex, StepSize) # Make the temporary vspace file
	create_tmp_prior_files(RunIndex, og_triname, dst_fold) # Make the temporary prior files

	# go through initial set and extract the values of the converging parameters for all sims
	body = []
	variable = []
	finit = []
	for i in params_to_conv: # extract info on the bodies, variables to converge, and whether they're final or initial
		body.append(i.split(',')[0])
		variable.append('('+i.split(',')[1]+')')
		finit.append(i.split(',')[2])

	body = np.array(body)
	variable = np.array(variable)
	finit = np.array(finit)

	for subdir, dirs, files in os.walk(dst_fold): # Loop through all the sims and extract info from each log
		if subdir != dst_fold:
			curr_file = open(os.path.join(subdir, vplanet_logfile), 'r')
			curr_lines = curr_file.readlines()
			curr_file.close()
			finithold = None
			currbody = None
			for i in range(len(curr_lines)):
				if len(curr_lines[i].split()) > 2:
					if curr_lines[i].split()[1] == 'FINAL':
						finithold = 'final'
#						print(curr_lines[i].split()[1])
					elif curr_lines[i].split()[1] == 'INITIAL':
						finithold = 'initial'
#						print(curr_lines[i].split()[1])
					elif curr_lines[i].split()[1] == 'BODY:':
						currbody = curr_lines[i].split()[2]
#						print(currbody)
					elif curr_lines[i].split()[0] in variable:
#						index = variable.index(curr_lines[i].split()[0])
						index = np.where(variable == curr_lines[i].split()[0])[0]
						if currbody in body[index] and finithold in finit[index]:
							bodyindex = np.where(body[index] == currbody)[0]
							finitindex = np.where(finit[index] == finithold)[0]
							trueindex = None
							for g in bodyindex:
#								print('g: ',g)
								for h in finitindex:
#									print('h: ',h)
									if g == h:
										trueindex = g
							if trueindex is not None:
#								print('body '+str(body[index][trueindex])+', variable '+str(variable[index][trueindex])+'_'+str(finit[index][trueindex])+' AKA '+str(curr_lines[i].split()[0]))
#								print('In params_to_conv: '+str(np.array(params_to_conv)[index][trueindex]))
								converge_dict[np.array(params_to_conv)[index][trueindex]].append(float(curr_lines[i].split()[len(curr_lines[i].split()) - 1]))

	# Check to make sure all requested parameters are real
	for i in params_to_conv:
		if converge_dict[i] == []:
			raise IOError('%s is not being simulated by VPLanet, or is not a VPLanet recognizable parameter. Please check spelling, modules used, etc.' % i)

	# --------------------------------------- LOOP START --------------------------------------------------
	# At this point, the initial training set has been processed and recorded. Now the loop will begin systematically adding steps of simulations of size StepSize
	# After every step, the code will check for convergence. At this point KL convergence is the only usable method, more may be added in future.
	Num_Convs = 0 # Keep track of the number of times the code reports consecutive convergence
	Totnumconvs = 0 # Keep track of total number of times convergence is reached
	quants = list(np.linspace(0.01, 0.99, 99)) # array from 0.01 to 0.99 by 0.01 for the quantile calculations
	f = open('vconverge_results.txt', 'w')
	
	while Num_Convs < ConvNum and RunIndex <= MaxSteps: # until the number of convergences has been met or the max number of allowed steps has been taken

		if ConvMethod == 'KL_Quantiles': # If convergence methods use quantiles, build up previous quantiles
			prev_quant = {} # get quantiles before next step is taken
			for i in params_to_conv:
				prev_quant[i] = np.quantile(converge_dict[i], quants)
		elif ConvMethod == 'KS_pval'or ConvMethod == 'KS_statistic': # If convergence method uses true values, save previous true values before next step
			prev_dict = copy.deepcopy(converge_dict) # deep copy so the dictionaries don't refer to the same object

		# Run Vspace
		# Run Multi-planet
		os.system('python -m vspace vconverge_tmp/vspace_tmp.in')
		#os.system('multi-planet vconverge_tmp/vspace_tmp.in')
		os.system('python -m multiplanet vconverge_tmp/vspace_tmp.in')

		# Go through sims on this step and append to the running list of values (converge_dict)
		for subdir, dirs, files in os.walk('vconverge_tmp/Step_'+str(RunIndex)): # Loop through all the sims and extract info from each log
			if subdir != 'vconverge_tmp/Step_'+str(RunIndex):
				shutil.copytree(subdir, os.path.join(dst_fold, subdir.split('/')[len(subdir.split('/'))-1])) # Copy sim into user defined dest folder so they have the data after tmp files are deleted
				curr_file = open(os.path.join(subdir, vplanet_logfile), 'r')
				curr_lines = curr_file.readlines()
				curr_file.close()
				finithold = None
				currbody = None
				for i in range(len(curr_lines)):
					if len(curr_lines[i].split()) > 2:
						if curr_lines[i].split()[1] == 'FINAL':
							finithold = 'final'
#							print(curr_lines[i].split()[1])
						elif curr_lines[i].split()[1] == 'INITIAL':
							finithold = 'initial'
#							print(curr_lines[i].split()[1])
						elif curr_lines[i].split()[1] == 'BODY:':
							currbody = curr_lines[i].split()[2]
#							print(currbody)
						elif curr_lines[i].split()[0] in variable:
							index = np.where(variable == curr_lines[i].split()[0])[0]
							if currbody in body[index] and finithold in finit[index]:
								bodyindex = np.where(body[index] == currbody)[0]
								finitindex = np.where(finit[index] == finithold)[0]
								trueindex = None
								for g in bodyindex:
#									print('g: ',g)
									for h in finitindex:
#										print('h: ',h)
										if g == h:
											trueindex = g
								if trueindex is not None:
#									print('body '+str(body[index][trueindex])+', variable '+str(variable[index][trueindex])+'_'+str(finit[index][trueindex])+' AKA '+str(curr_lines[i].split()[0]))
#									print('In params_to_conv: '+str(np.array(params_to_conv)[index][trueindex]))
									converge_dict[np.array(params_to_conv)[index][trueindex]].append(float(curr_lines[i].split()[len(curr_lines[i].split()) - 1]))
		if ConvMethod == 'KL_Quantiles':
			curr_quant = {} # get current quantiles, after step has been taken
			for i in params_to_conv:
				curr_quant[i] = np.quantile(converge_dict[i], quants)

		f.write('\n')
		f.write('----- Step '+str(RunIndex)+' -----\n\n')
		if ConvMethod == 'KS_statistic':
			f.write('KS Statistics:\n')
		elif ConvMethod == 'KS_pval':
			f.write('KS P-Values:\n')
		elif ConvMethod == 'KL_Quantiles':
			f.write('KL Divergence:\n')

		# For the KS test, check the p-value or the ks statistic (based on user input)
		if ConvMethod == 'KS_pval' or ConvMethod == 'KS_statistic':
			ks = {}
			converged = []
			for i in params_to_conv:
				ksstat, pval = ks_2samp(prev_dict[i], converge_dict[i])
				if ConvMethod == 'KS_pval':
					ks[i] = pval
					f.write(i+' ----- '+str(pval)+'\n')
				elif ConvMethod == 'KS_statistic':
					ks[i] = ksstat
					print(i+' KS Statistic ---- '+str(ksstat))
					f.write(i+' ----- '+str(ksstat)+'\n')
				if ConvMethod == 'KS_pval' and pval >= ConvCondit:
					converged.append(True)
				elif ConvMethod == 'KS_statistic' and ksstat <= ConvCondit:
					converged.append(True)
				else:
					converged.append(False) # If even 1 parameter has not converged, it's False and break the loop

		f.write('\n')

		if False not in converged: # If all params have converged, add 1 to the Num_Convs (number of times the model has converged consecutively)
			Num_Convs = Num_Convs + 1
			Totnumconvs = Totnumconvs + 1
			f.write('Converged = True\n')
			f.write('Number of Consecutive Convergences: '+str(Num_Convs)+'\n')
		elif False in converged: # If even one param did not converged, make sure Num_Convs (number of times the model has converged consecutively) is set to 0
			Num_Convs = 0
			f.write('Converged = False\n')

		RunIndex = RunIndex + 1

		# Create new vspace.in and prior files
		create_tmp_vspin(vspFile, RunIndex, StepSize) # Make the temporary vspace file
		create_tmp_prior_files(RunIndex, og_triname, 'vconverge_tmp') # Make the temporary prior files

	if RunIndex >= MaxSteps:
		print('MaxSteps reached')
		print('Number of Convergences: '+str(Num_Convs))


	f.write('\n')
	f.write('-------------- VCONVERGE STATS ----------------\n\n')
	f.write('Number of Steps Taken: '+str(RunIndex-1)+'\n')
	totalsims = initialsims + ((RunIndex-1)*StepSize)
	f.write('Total Number of Simulations Ran: '+str(totalsims)+'\n')
	if RunIndex < MaxSteps:
		f.write(str(ConvNum)+' Consecutive Convergences Achieved, vconverge run sucessful\n')
	else:
		f.write('Max Steps Reached. '+str(Num_Convs)+' Consecutive Convergences Achieved\n')
	f.write('Total Number of Convergences during vconverge run: '+str(Totnumconvs)+'\n')
#	if False not in converged:
#		f.write('Converged - True\n')
#	else:
#		f.write('Converged - False\n')
#	f.write('Number of Consecutive Convergences: '+str(Num_Convs)+'\n\n')
#	if ConvMethod == 'KS_pval':
#		f.write('p-values from KS Test:\n')
#	elif ConvMethod == 'KS_statistic':
#		f.write('KS Statistics from KS Test:\n')
#	if ConvMethod == 'KS_pval' or ConvMethod == 'KS_statistic':
#		for i in params_to_conv:
#			f.write(i+' ----- '+str(ks[i])+'\n')
	f.close()

	# Save dictionary of converging params
	dicthold = json.dumps(converge_dict)
	cparams = open(os.path.join(dst_fold, 'Converged_Param_Dictionary.json'), 'w')
	json.dump(dicthold, cparams)
	cparams.close()

	for i in converge_dict:
		for k in range(len(converge_dict[i])-1):
			if converge_dict[i][k] == converge_dict[i][k+1]:
				same = True
			else:
				same = False
				break
		if same == True:
			converge_dict[i][0] = converge_dict[i][0]+1e-10

        ######## UNCOMMENT BELOW FOR CORNERS PLOTS ##############

	# Group converging parameters by body for corners plotting
	#hold = []
	#for i in range(len(body)):
	#	if body[i] not in hold:
	#		hold.append(body[i])
	#groups = []
	#for i in hold:
	#	hold2 = []
	#	for k in range(len(body)):
	#		if i == body[k]:
	#			hold2.append(params_to_conv[k])
	#	groups.append(hold2)

	# Make Corners plots
	#for i in groups:
	#	currbody = i[0].split(',')[0]
	#	cornerdat = []
	#	for k in range(len(converge_dict[i[0]])):
	#		hold3 = []
	#		for j in range(len(i)):
	#			hold3.append(converge_dict[i[j]][k])
	#		cornerdat.append(hold3)
	#	fig = corner.corner(cornerdat, quantiles=[0.16, 0.5, 0.84], color='g', show_titles=True, title_fmt='.4f', labels=i)
	#	plt.savefig(os.path.join(dst_fold, currbody+'.png'))
	#	plt.close(fig)
	

	return converge_dict, converged


cdic, conved = vconverge(sys.argv[1])
