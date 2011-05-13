################################################
#                                              #
#    PHASE B: KINGFISH Spectroscopic Pipeline  #
#       BETA version tested in HIPE 7.0.815    #
#         person to blame: Kevin Croxall       #
#            (aside from the NHSC)             #
#                 May 13, 2011                 #
################################################

# A list of the Phase A pools that are to be processed
# need to be saved and referenced here at the begining of the Phase A pipeline.  

phasealist = simpleAsciiTableReader(file = "/Users/kcroxall/phaseA.dat") #UPDATE to the correct file location
ndim = phasealist[0].data.dimensions[0]

from herschel.pacs.signal import MaskViewer
verbose = 1
n=1
for n in range(0,ndim):
	# ------------------------------------------------------------------------------
	#         Processing      Level 1 -> Level 2
	# ------------------------------------------------------------------------------
	# GET THE DATA
	# If you instead begin here from a slicedCubes you saved to a pool with saveSlicedCopy, 
	# then (using the poolName you set when saving):
	name=str(phasealist[0].data[n]) 
	slicedCubes = readSliced(name)
	galname = slicedCubes.meta["object"].string
	calTree=getCalTree()
	#then do a for loop over the length of the array for the full processing...
	linelist = Double1d()
	linelist.append(slicedCubes.refs[0].product["BlockTable"]['LineId'].data[0])
	for i in range (1,slicedCubes.getRefs().size()):
		if (slicedCubes.refs[i].product["BlockTable"]['LineId'].data[0] != slicedCubes.refs[i-1].product["BlockTable"]['LineId'].data[0]): 
			linelist.append(slicedCubes.refs[i].product["BlockTable"]['LineId'].data[0])
	
	linewave = Double1d()
	for qq in range (0,linelist.dimensions[0]):
		count=0
		sum=0
		for i in range (0,slicedCubes.getRefs().size()):
			if (str(slicedCubes.refs[i].product["BlockTable"]['LineId'].data[0]) == str(int(linelist[qq]))): 
				count = count+1
				sum = sum+slicedCubes.refs[i].product["wave"].data[0]
		linewave.append(int(sum/count))
	
	for lineloop in range (0,linelist.dimensions[0]):
		if (linewave[lineloop] > 190): 
			linename = "NII205"
			linecent = 205.17823*(1+z)
		if (linewave[lineloop] > 140)&(linewave[lineloop] < 190): 
			linename = "CII"
			linecent = 157.7409*(1+z)
		if (linewave[lineloop] > 100)&(linewave[lineloop] < 140): 
			linename = "NII122"
			linecent = 121.89757*(1+z)
		if (linewave[lineloop] > 75)&(linewave[lineloop] < 100): 
			linename = "OIII"
			linecent = 88.356*(1+z)
		if (linewave[lineloop] > 60)&(linewave[lineloop] < 75): 
			linename = "OI"
			linecent = 63.183705
		if verbose: slicedSummary(slicedCubes)
#PHIL NOTED THAT CONTEXT LAYERS CHANGED IN 7.0.1786 MAKE SURE THINGS ARE READING AND REFERENCEING CORRECTLY
		lineId      = [int(linelist[lineloop])]
		wavelength  = []
		rasterLine  = []
		rasterCol   = []
		nodCycle    = ""
		scical      = ""
		sliceNumber = []
		onOff="ON"
		sCubesOn = selectSlices(slicedCubes, lineId=lineId, wavelength=wavelength, rasterLine=rasterLine,\
         rasterCol=rasterCol, onOff=onOff, nodCycle=nodCycle, scical=scical, sliceNumber=sliceNumber, verbose=verbose)
		onOff="OFF"
		sCubesOff = selectSlices(slicedCubes, lineId=lineId, wavelength=wavelength, rasterLine=rasterLine,\
         rasterCol=rasterCol, onOff=onOff, nodCycle=nodCycle, scical=scical, sliceNumber=sliceNumber, verbose=verbose)
###################### BETA VERSION of using the BASELINE AS AN OFF ######################################
###################### Could be supplemented by the transient correction?
###################### Needs a better model and to treats sharp jumps...
		Dotcloudmod = getSlicedCopy(sCubesOn)
		for rasti in range (0,Dotcloudmod.getRefs().size()):
			for q in range (0,5):
				for qq in range (0,5):
					for j in range (0,16):
						wave = Double1d()
						flux = Double1d()
						reset = Double1d()
						for i in range (0,Dotcloudmod.refs[rasti].product["wave"].data.dimensions[0]/16):
							wave.append(Dotcloudmod.refs[rasti].product["wave"].data[j+i*16,qq,q])
							reset.append(Dotcloudmod.refs[rasti].product["Status"]["RESETINDEX"].data[j+i*16])
							flux.append(Dotcloudmod.refs[rasti].product["flux"].data[j+i*16,qq,q])
						ind = flux.where((wave<linecent-0.6).or(wave>linecent+0.6))
						myModel = PolynomialModel(15)	
						myFitter = Fitter(reset[ind], myModel)
						ind2 = flux.where(IS_NAN(flux))
						flux2=flux
						flux2[ind2] = 180
						flux2[ind2] = SUM(flux2)/flux2.dimensions[0]
						fitresults = myFitter.fit(flux2[ind])
						nsize = flux.dimensions[0]
						u = MIN(reset) + Double1d.range(nsize) * ((MAX(reset) - MIN(reset)) / nsize)
						umodel = myModel(u)
						newflux = flux - umodel
						for i in range (0,Dotcloudmod.refs[rasti].product["wave"].data.dimensions[0]/16):
							Dotcloudmod.refs[rasti].product["flux"].data[j+i*16,qq,q] = newflux[i]	
		sCubesOn = getSlicedCopy(Dotcloudmod)
		Dotcloudmod = getSlicedCopy(sCubesOff)
		for rasti in range (0,Dotcloudmod.getRefs().size()):
			for q in range (0,5):
				for qq in range (0,5):
					for j in range (0,16):
						wave = Double1d()
						flux = Double1d()
						reset = Double1d()
						for i in range (0,Dotcloudmod.refs[rasti].product["wave"].data.dimensions[0]/16):
							wave.append(Dotcloudmod.refs[rasti].product["wave"].data[j+i*16,qq,q])
							reset.append(Dotcloudmod.refs[rasti].product["Status"]["RESETINDEX"].data[j+i*16])
							flux.append(Dotcloudmod.refs[rasti].product["flux"].data[j+i*16,qq,q])
						ind = flux.where((wave<157.5).or(wave>158.2))
						myModel = PolynomialModel(15)	
						myFitter = Fitter(reset[ind], myModel)
						ind2 = flux.where(IS_NAN(flux))
						flux2=flux
						flux2[ind2] = 180
						flux2[ind2] = SUM(flux2)/flux2.dimensions[0]
						fitresults = myFitter.fit(flux2[ind])
						nsize = flux.dimensions[0]
						u = MIN(reset) + Double1d.range(nsize) * ((MAX(reset) - MIN(reset)) / nsize)
						umodel = myModel(u)
						newflux = flux - umodel
						for i in range (0,Dotcloudmod.refs[rasti].product["wave"].data.dimensions[0]/16):
							Dotcloudmod.refs[rasti].product["flux"].data[j+i*16,qq,q] = newflux[i]	
		sCubesOff = getSlicedCopy(Dotcloudmod)
############################################################
############################################################
############################################################
if verbose: 
		slicedSummary(sCubesOn)
		slicedSummary(sCubesOff)
		x,y = 2,2
		p6 = plotCubes(sCubesOn,[],x=x,y=y,masks=sCubesOn.refs[0].product.maskTypes)
		p6 = plotCubes(sCubesOff,p6,x=x,y=y,masks=sCubesOff.refs[0].product.maskTypes)
		plotMasks = String1d(["GLITCH"])
		p6b = plotCubes(sCubesOn,[],x=x,y=y,masks=sCubesOn.refs[0].product.maskTypes)
		p6b = plotCubesMasks(sCubesOn,p6b,x=x,y=y,masks=plotMasks)
		plotname = galname + "_" + linename + "pretrim.eps"
		p6.saveAsEPS(plotname)     #saves in the kcroxall directory
		
		# Note that for the final cube rebinning it is recommended that: 
		#    - For single scan SED or Nyquist sampled range scans, it is recommended to perform the final rebinning with 
		#        oversample=2, upsample>=1, which corresponds to the native resolution of the instrument
		#    - For line scan, deep range scans or Nyquist sampled range scans with repetition factor > 1, 
		#        oversampling > 2 is made possible by the high degree of redundancy provided by the observation
		waveGrid=wavelengthGrid(sCubesOn.refs[0].product, oversample=2, upsample=1, calTree = calTree)
		sCubesOn  = activateMasks(sCubesOn, String1d(["GLITCH","UNCLEANCHOP","RAWSATURATION","SATURATION","GRATMOVE", "BADPIXELS"]), exclusive = True)
		sCubesOff = activateMasks(sCubesOff, String1d(["GLITCH","UNCLEANCHOP","RAWSATURATION","SATURATION","GRATMOVE", "BADPIXELS"]), exclusive = True)
		sCubesOn  = specFlagOutliers(sCubesOn, waveGrid, nSigma=3, nIter=4,saveStats=1)
		sCubesOff = specFlagOutliers(sCubesOff, waveGrid, nSigma=2, nIter=4,saveStats=1)
		if verbose:
			x,y = 2,2
			p7 = plotCubes(sCubesOn,[],x=x,y=y,masks=sCubesOn.refs[0].product.maskTypes)
			plotMasks = String1d(["OUTLIERS"])
			p7 = plotCubesMasks(sCubesOn,p7,x=x,y=y,masks=plotMasks)
			p8 = plotCubes(sCubesOn,[],x=x,y=y,masks=sCubesOn.refs[0].product.maskTypes)
			plotname = galname + "_" + linename + "posttrim.eps"
			p8.saveAsEPS(plotname)     #saves in the kcroxall directory
			plotname = galname + "_" + linename + "posttrim_flag.eps"
			p7.saveAsEPS(plotname)     #saves in the kcroxall directory
		
		sCubesOn  = activateMasks(sCubesOn, String1d(["GLITCH","UNCLEANCHOP","RAWSATURATION","SATURATION","GRATMOVE", "BADPIXELS", "OUTLIERS"]), exclusive = True)
		sCubesOff = activateMasks(sCubesOff, String1d(["GLITCH","UNCLEANCHOP","RAWSATURATION","SATURATION","GRATMOVE", "BADPIXELS", "OUTLIERS"]), exclusive = True)
		slicedRebinnedCubesOn = specWaveRebin(sCubesOn, waveGrid)
		slicedRebinnedCubesOff = specWaveRebin(sCubesOff, waveGrid)
		if verbose:
			slicedSummary(slicedRebinnedCubesOn)
			slicedSummary(slicedRebinnedCubesOff)
			# Sky footprint: the second is overplotted on the first, p8 is first created and then replotted on
			p8 = plotCubesRaDec(slicedRebinnedCubesOn)
			p8 = plotCubesRaDec(slicedRebinnedCubesOff,p8)
			plotname = galname + "_" + linename + "_radec.eps"
			p8.saveAsEPS(plotname)     #saves in the kcroxall directory
		
		slicedRebinnedCubesOn  = specAverageCubes(slicedRebinnedCubesOn)
		slicedRebinnedCubesOff = specAverageCubes(slicedRebinnedCubesOff)
		if verbose: p7 = plotCubes(slicedRebinnedCubesOn,[],x=x,y=y)
		if verbose:
			x,y = 2,2
			# here we are adding to a previously created plot, but if you deleted plot "p6" from above, then first
			p6 = plotCubes(sCubesOn,[],x=x,y=y,masks=sCubesOn.refs[0].product.maskTypes)
			p6 = plotCubes(sCubesOff,p6,x=x,y=y,masks=sCubesOff.refs[0].product.maskTypes)
			# and then
			p6.addLayer(LayerXY(slicedRebinnedCubesOn.refs[0].product["waveGrid"].data,slicedRebinnedCubesOn.refs[0].product["image"].data[:,x,y]))
			p6.addLayer(LayerXY(slicedRebinnedCubesOff.refs[0].product["waveGrid"].data,slicedRebinnedCubesOff.refs[0].product["image"].data[:,x,y]))
			plotname = galname + "_" + linename + "withspec.eps"
			p6.saveAsEPS(plotname)     #saves in the kcroxall directory
		
		slicedRebinnedCubesOnERR = getSlicedCopy(slicedRebinnedCubesOn)
		for slice in range(len(slicedRebinnedCubesOnERR.refs)):
			slicedRebinnedCubesOnERR.refs[slice].product["image"].data = slicedRebinnedCubesOnERR.refs[slice].product["stddev"].data / SQRT(slicedRebinnedCubesOnERR.refs[slice].product["exposure"].data)
		
		# Subtract the off.  
		# Rebin and subtract. It should be that all off positions were 
		# averaged in specAverageCubes (=> slicedRebinnedCubesOff has only one slice)
		# First create a slicedCubes to hold the difference slicedCubes, and this is the easiest way to do this:
		slicedDiffCubes = getSlicedCopy(slicedRebinnedCubesOn)
		# Then fill that. As there should be only one slice in slicedRebinnedCubesOff, we subtract that single slice from 
		# all the On slices.
		for slice in range(len(slicedRebinnedCubesOn.refs)):
			slicedDiffCubes.refs[slice].product.setFlux(slicedRebinnedCubesOn.refs[slice].product.image - slicedRebinnedCubesOff.refs[1].product.image)
		
#deal with more than 1 off cube in the data...
# ADD IN A TEST TO DETERMINE THE CLOSEST OFF...
		slicedModOffCube = getSlicedCopy(slicedRebinnedCubesOff)
		for q in range (0,5):
			for qq in range (0,5):
				y = slicedModOffCube.refs[1].product["image"].data[15:slicedModOffCube.refs[1].product["image"].data.dimensions[0]-15,q,qq]
				x = slicedModOffCube.refs[1].product["waveGrid"].data[15:slicedModOffCube.refs[1].product["image"].data.dimensions[0]-15]
				p=PlotXY(x,y)
				#knots = Double1d.range( 17 ) * 10      # make knots from 0 to 160
				#myModel = CubicSplinesModel(knots)
				myModel = PolynomialModel(3)	
				myFitter = Fitter(x, myModel)
				ind = y.where(IS_NAN(y))
				y2=y
				y2[ind] = 180
				print ind
				print SUM(y)/y.dimensions[0]
				y2[ind] = SUM(y2)/y2.dimensions[0]
				yweights = slicedModOffCube.refs[1].product["stddev"].data[15:slicedModOffCube.refs[1].product["image"].data.dimensions[0]-15,q,qq] / SQRT(slicedModOffCube.refs[1].product["exposure"].data[15:slicedModOffCube.refs[1].product["image"].data.dimensions[0]-15,q,qq])
				yweights[ind]=0.0001
				fitresults = myFitter.fit(y2)
#				fitresults = myFitter.fit(y2,yweights)
				nsize = y.dimensions[0]+30
				u = MIN(x) + Double1d.range(nsize) * ((MAX(x) - MIN(x)) / nsize)
				# Apply the model
				umodel = myModel(u)
				slicedModOffCube.refs[1].product["image"].data[:,q,qq] = umodel
				nsize = y.dimensions[0]
				u = MIN(x) + Double1d.range(nsize) * ((MAX(x) - MIN(x)) / nsize)
				# Apply the model
				umodel = myModel(u)
				p[1] = LayerXY(x, umodel, name = "Fit", color = java.awt.Color.green)
		
		slicedDiffCubesMod = getSlicedCopy(slicedRebinnedCubesOn)
		for slice in range(len(slicedRebinnedCubesOn.refs)):
			slicedDiffCubesMod.refs[slice].product.setFlux(slicedRebinnedCubesOn.refs[slice].product.image - slicedModOffCube.refs[1].product.image)
		if verbose:
			slice = 0
			p10 = plotCube5x5(slicedModOffCube.refs[1].product)
			p11 = plotCube5x5(slicedRebinnedCubesOff.refs[1].product)
			p9 = plotCube5x5(slicedDiffCubesMod.refs[slice].product)
			p9 = plotCube5x5(slicedDiffCubes.refs[slice].product)
			p10 = plotCube5x5(slicedRebinnedCubesOn.refs[slice].product)
			plotname = galname + "_" + linename + "_5x5_slice0_nooff.eps"
			p10.saveAsEPS(plotname)     #saves in the kcroxall directory
			plotname = galname + "_" + linename + "_5x5_slice0_offsub.eps"
			p9.saveAsEPS(plotname)     #saves in the kcroxall directory
		name=galname + "_"+linename+"_suboff_pipe_PhaseB"
		saveSlicedCopy(slicedDiffCubes,name)
		name=galname + "_"+linename+"_submod_pipe_PhaseB"
		saveSlicedCopy(slicedDiffCubesMod,name)
		name=galname + "_"+linename+"_nooff_pipe_PhaseB"
		saveSlicedCopy(slicedRebinnedCubesOn,name)
		name=galname + "_"+linename+"_err_pipe_PhaseB"
		saveSlicedCopy(slicedRebinnedCubesOnERR,name)


# WHICH ONE TO SAVE??????
# DO WE WANT ONLY MAPS OF EVERYTHING TOGETHER?
# OR ALSO THE MAPS OF THE INDIVIDUAL AORS???


# End Phase B

#projectedCube_difmod = specProject(slicedDiffCubesMod, outputPixelsize=2.85)
#projectedCube_diff = specProject(slicedDiffCubes, outputPixelsize=2.85)
#projectedCube_nooff = specProject(slicedRebinnedCubesOn, outputPixelsize=2.85)

#projectedCube_nooff = specProject(slicedRebinnedCubesOn, outputPixelsize=1)

print "you have arrived"

print "CONGRATULATIONS! Phase B complete!"

################################################
#                                              #
#    PHASE B: KINGFISH Spectroscopic Pipeline  #
#       BETA version tested in HIPE 7.0.815    #
#         person to blame: Kevin Croxall       #
#            (aside from the NHSC)             #
#                 May 13, 2011                 #
################################################
