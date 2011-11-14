################################################
#                                              #
#    PHASE A: KINGFISH Spectroscopic Pipeline  #
#      BETA version tested in HIPE 8.0.3067    #
#         person to blame: Kevin Croxall       #
#            (aside from the NHSC)             #
#                 Nov 14, 2011                 #
################################################

# All data must be imported into HIPE data pools.  The pipeline will not unpack the HSA tarballs
# and import them into HIPE on its own.  A list of the pools and OBSIDs that are to be processed
# need to be saved and referenced here at the begining of the Phase A pipeline.  Data are saved 
# into new pools that are sorted by AOR and line.  The frames can be saved chopped up by raster
# or as one large data stream.
def getMedian(numericValues):
  theValues = SORT(numericValues)
  if len(theValues) % 2 == 1:
    return theValues[(len(theValues)+1)/2-1]
  else:
    lower = theValues[len(theValues)/2-1]
    upper = theValues[len(theValues)/2]
    return (float(lower + upper)) / 2 

def ktrans_pos(frame):
	table = TableDataset()
	ascii = AsciiTableTool()
	raarray = Double1d()
	decarray = Double1d()
	rastarray = Double1d()
	spatialarray = Double1d()
	obsid = str(frame.meta["obsid"].long)
	galname = frame.meta["object"].string
	for rasti in range (0,frame.getRefs().size()):
		for spatialpix in range (0,25):
			#print "Spatial Pixel ",spatialpix
			RApix = getMedian(frame.refs[rasti].product["Ra"].data[8,spatialpix,:])
			DECpix = getMedian(frame.refs[rasti].product["Dec"].data[8,spatialpix,:])
			raarray.append(RApix)
			decarray.append(DECpix)
			rastarray.append(rasti)
			spatialarray.append(spatialpix)
	table["Raster"] = Column(rastarray)
	table["Spatialpix"] = Column(spatialarray)
	table["RA"] = Column(raarray)
	table["DEC"] = Column(decarray)
	ascii.formatter=FixedWidthFormatter(sizes=[10,10,25,25])
	tabname = "/home/kcroxall/" + galname + "_" + obsid + "_postest.tab"   # UPDATE PATH
	ascii.save(tabname, table)

def ktrans_frame(frame):
	z=frame.meta["redshiftValue"].double/300000
	for rasti in range (1,frame.getRefs().size()):
		print "Raster", rasti
		wavecent = ((frame.refs[rasti].product["BlockTable"]["MinWave"].data[0] + frame.refs[rasti].product["BlockTable"]["MaxWave"].data[0])/2)
		chanwidth = 25
		linewidth = 0.03
		linecent = wavecent
		if (wavecent > 190): 
			chanwidth = 55
			linewidth = 0.18
			linecent = 205.17823*(1+z)
		if (wavecent > 150)&(wavecent < 170): 
			linewidth = 0.18
			chanwidth = 45
			linecent = 157.7409*(1+z)
		if (wavecent > 110)&(wavecent < 130): 
			chanwidth = 45
			linewidth = 0.18
			linecent = 121.89757*(1+z)
		if (wavecent > 75)&(wavecent < 95): 
			chanwidth = 45
			linewidth = 0.065
			linecent = 88.356*(1+z)
		if (wavecent > 60)&(wavecent < 70): 
			chanwidth = 75
			linewidth = 0.065
			linecent = 63.183705*(1+z)
		print "Line at ", linecent
		for spatialpix in range (0,25):
			print "Spatial Pixel ",spatialpix
			for specpix in range (1,17):
				#print "Spectral Pixel ", specpix
				wave = frame.refs[rasti].product["Wave"].data[specpix,spatialpix,:]
				flux = frame.refs[rasti].product["Signal"].data[specpix,spatialpix,:]
				fluxnew = frame.refs[rasti].product["Signal"].data[specpix,spatialpix,:]
				reset = frame.refs[rasti].product["Status"]["RESETINDEX"].data[:]
				if (flux[33] == flux[33])&(flux[53] == flux[53]):
					for pix in range (0,reset.dimensions[0]-1):
						#print pix
						dist = ABS(reset[:] - reset[pix])
						nearby = dist.where((dist<chanwidth).and((wave<linecent-linewidth).or(wave>linecent+linewidth)))
						subval = getMedian(flux[nearby])
						flux2 = flux[nearby]
						dist2 = ABS(flux2 - subval)
						mask = dist2.where(IS_FINITE)
						dist2c = dist2[mask]
						meandist2 = SUM(dist2c)/dist2c.dimensions[0]
						nearby2 = dist2.where(dist2<meandist2)
						subval2 = SUM(flux2[nearby2])/flux2[nearby2].dimensions[0]
						fluxnew[pix] = flux[pix]-subval2
				frame.refs[rasti].product["Signal"].data[specpix,spatialpix,:] = fluxnew[:]
			System.gc()
	del(wave,flux,fluxnew,reset,spatialpix,specpix,rasti,dist,nearby,subval,flux2,dist2,meandist2,nearby2,subval2)
	return (frame)


Pversion = "PhaseA_v2.0"
Hversion = "8.0.3067"
Cversion = "32"

poollist = simpleAsciiTableReader(file = "/home/kcroxall/ngc3077test.lst")     #UPDATE to the correct file location
obsidlist = simpleAsciiTableReader(file = "/home/kcroxall/ngc3077test.lst")   #UPDATE to the correct file location
ndim = poollist[0].data.dimensions[0]

from java.awt import Color

for n in range(0,ndim):
	# ------------------------------------------------------------------------------
	#        Preparation
	# ------------------------------------------------------------------------------
	print poollist[0].data[n]
	print obsidlist[0].data[n]
	print n,poollist[0].data.dimensions[0]
	obs    = getObservation(obsidlist[0].data[n], verbose=True, useHsa=0, poolLocation=None, poolName=str(poollist[0].data[n]))
	# verbose: 0 - silent, execute the pipeline only
	#	   1 - will trigger diagnostic output on the screen, plots, and displays
	verbose = 0
	if verbose: obsSummary(obs)
	cameraR = "red"
	cameraB = "blue"
	calTree = getCalTree(obs=obs)
	#calTreeb = obs.calibration	# to use the cal files packaged with the obs rather than the latest on-machine cal products
	if verbose: 
		print calTree
		print calTree.common
		print calTree.spectrometer
	
	obs = checkForAnomaly70(obs)
	pacsPropagateMetaKeywords(obs,'0', obs.level0)
	level0 = PacsContext(obs.level0)
	
	slicedFramesR  = SlicedFrames(level0.fitted.getCamera(cameraR).product)		#	scientific data
	slicedRawRampR = level0.raw.getCamera(cameraR).product				#	raw data for one pixel
	slicedDmcHeadR = level0.dmc.getCamera(cameraR).product				#	the mechanisms' status information
	slicedFramesB  = SlicedFrames(level0.fitted.getCamera(cameraB).product)
	slicedRawRampB = level0.raw.getCamera(cameraB).product
	slicedDmcHeadB = level0.dmc.getCamera(cameraB).product    
	galname = slicedFramesR.meta["object"].string
	slicedFramesR.meta["Pversion"]=StringParameter(Pversion)
	slicedFramesR.meta["Hversion"]=StringParameter(Hversion)
	slicedFramesR.meta["Cversion"]=StringParameter(Cversion)
	slicedFramesB.meta["Pversion"]=StringParameter(Pversion)
	slicedFramesB.meta["Hversion"]=StringParameter(Hversion)
	slicedFramesB.meta["Cversion"]=StringParameter(Cversion)
	
	if verbose:
		slicedSummary(slicedFramesR)			# Overview of the basic structure of the data
		slicedSummary(slicedFramesB)
		p0 = slicedSummaryPlot(slicedFramesR,signal=1)	# Grating position & raw signal of central pixel
# ------------------------------------------------------------------------------
#        Processing      Level 0 -> Level 0.5
# ------------------------------------------------------------------------------
	print "Level 0 -> Level 0.5"
	slicedFramesR = specFlagSaturationFrames(slicedFramesR, rawRamp = slicedRawRampR, calTree=calTree, copy=1)
	slicedFramesR = specConvDigit2VoltsPerSecFrames(slicedFramesR, calTree=calTree)
	slicedFramesR = detectCalibrationBlock(slicedFramesR)
	slicedFramesR = addUtc(slicedFramesR, obs.auxiliary.timeCorrelation)
	slicedFramesR = specAddInstantPointing(slicedFramesR, obs.auxiliary.pointing, calTree = calTree, orbitEphem = obs.auxiliary.orbitEphemeris, horizonsProduct = obs.auxiliary.horizons)    
	slicedFramesR = specExtendStatus(slicedFramesR, calTree=calTree)
	slicedFramesR = convertChopper2Angle(slicedFramesR, calTree=calTree)
	slicedFramesR = specAssignRaDec(slicedFramesR, calTree=calTree)
	print "Now B"
	slicedFramesB = specFlagSaturationFrames(slicedFramesB, rawRamp = slicedRawRampB, calTree=calTree, copy=1)
	slicedFramesB = specConvDigit2VoltsPerSecFrames(slicedFramesB, calTree=calTree)
	slicedFramesB = detectCalibrationBlock(slicedFramesB)
	slicedFramesB = addUtc(slicedFramesB, obs.auxiliary.timeCorrelation)
	slicedFramesB = specAddInstantPointing(slicedFramesB, obs.auxiliary.pointing, calTree = calTree, orbitEphem = obs.auxiliary.orbitEphemeris, horizonsProduct = obs.auxiliary.horizons)    
	slicedFramesB = specExtendStatus(slicedFramesB, calTree=calTree)
	slicedFramesB = convertChopper2Angle(slicedFramesB, calTree=calTree)
	slicedFramesB = specAssignRaDec(slicedFramesB, calTree=calTree)
	if verbose: 						# show footprints for selected raster position(s)
		slicedPlotPointingOnOff(slicedFramesR)
	print "WaveCalc R"
	slicedFramesR = waveCalc(slicedFramesR, calTree=calTree)
	slicedFramesR = specCorrectHerschelVelocity(slicedFramesR, obs.auxiliary.orbitEphemeris, obs.auxiliary.pointing, obs.auxiliary.timeCorrelation)
	slicedFramesR = findBlocks(slicedFramesR, calTree = calTree)
	slicedFramesR = specFlagBadPixelsFrames(slicedFramesR, calTree=calTree)
	slicedFramesR = flagChopMoveFrames(slicedFramesR, dmcHead=slicedDmcHeadR, calTree=calTree)
	slicedFramesR = flagGratMoveFrames(slicedFramesR, dmcHead=slicedDmcHeadR, calTree=calTree)
	print "Now B"
	slicedFramesB = waveCalc(slicedFramesB, calTree=calTree)
	slicedFramesB = specCorrectHerschelVelocity(slicedFramesB, obs.auxiliary.orbitEphemeris, obs.auxiliary.pointing, obs.auxiliary.timeCorrelation)
	slicedFramesB = findBlocks(slicedFramesB, calTree = calTree)
	slicedFramesB = specFlagBadPixelsFrames(slicedFramesB, calTree=calTree)
	slicedFramesB = flagChopMoveFrames(slicedFramesB, dmcHead=slicedDmcHeadB, calTree=calTree)
	slicedFramesB = flagGratMoveFrames(slicedFramesB, dmcHead=slicedDmcHeadB, calTree=calTree)
	print "Delete superfluous products to ease RAM usage"
	del(slicedRawRampR,slicedRawRampB,slicedDmcHeadR,slicedDmcHeadB,level0)
	System.gc()
	if verbose:
		slicedSummary(slicedFramesR)						# an overview of the slicedFrames contents
		slicedSummary(slicedFramesB)
		maskSummary(slicedFramesR)
		maskSummary(slicedFramesB)
		p1 = slicedSummaryPlot(slicedFramesR,signal=0)		# Show the basic data structure, without the signal
		p1 = slicedSummaryPlot(slicedFramesB,signal=0)		
	print "Rules!"
	rules = [SlicingRule("LineId",1),SlicingRule("RasterLineNum",1),SlicingRule("RasterColumnNum",1),SlicingRule("NoddingPosition",1),SlicingRule("NodCycleNum",1),SlicingRule("IsOutOfField",1),SlicingRule("Band",1)]
	slicedFramesR = pacsSliceContext(slicedFramesR, slicingRules = rules, removeUndefined=1)
	slicedFramesR = specAddGratingCycleStatus(slicedFramesR)
	slicedFramesB= pacsSliceContext(slicedFramesB, slicingRules = rules, removeUndefined=1)
	slicedFramesB = specAddGratingCycleStatus(slicedFramesB)
	if verbose:
		slicedSummary(slicedFramesR)
		slicedSummary(slicedFramesB)
		p2 = slicedSummaryPlot(slicedFramesR,signal=0)
		p2 = slicedSummaryPlot(slicedFramesB,signal=0)
# ------------------------------------------------------------------------------
#         Processing      Level 0.5 -> Level 1
# ------------------------------------------------------------------------------
	print "Level 0.5 -> Level 1.0"
	slicedFramesR = activateMasks(slicedFramesR, String1d([" "]), exclusive = True)
	slicedFramesR = specFlagGlitchFramesQTest(slicedFramesR, copy=1)
	slicedFramesB = activateMasks(slicedFramesB, String1d([" "]), exclusive = True)
	slicedFramesB = specFlagGlitchFramesQTest(slicedFramesB, copy=1)
	if verbose:
		slicedSummary(slicedFramesR)
		slicedSummary(slicedFramesB)
		p3 = slicedSummaryPlot(slicedFramesR,signal=0)
		p3 = slicedSummaryPlot(slicedFramesB,signal=0)
		slice = 1
		p4 = plotSignalBasic(slicedFramesR, slice=slice)	# Single pixel, signal vs wavelength: Only unmasked datapoints are plotted
		p4 = plotSignalBasic(slicedFramesB, slice=slice)
		MaskViewer(slicedFramesR.get(slice))			# Inspect timeline of signals and masked signals
		MaskViewer(slicedFramesB.get(slice))
	print "ACTIVATE QUALITY INFO!"
	slicedFramesR = activateMasks(slicedFramesR, slicedFramesR.get(0).getMaskTypes())
	slicedFramesR = convertSignal2StandardCap(slicedFramesR, calTree=calTree)
	slicedFramesB = activateMasks(slicedFramesB, slicedFramesB.get(0).getMaskTypes())
	slicedFramesB = convertSignal2StandardCap(slicedFramesB, calTree=calTree)
	
	calBlockR = selectSlices(slicedFramesR,scical="cal").get(0)
	csResponseAndDarkR = specDiffCs(calBlockR, calTree = calTree)
	calBlockB = selectSlices(slicedFramesB,scical="cal").get(0)
	csResponseAndDarkB = specDiffCs(calBlockB, calTree = calTree)
	slicedFramesR = specSubtractDark(slicedFramesR, calTree=calTree)
	slicedFramesB = specSubtractDark(slicedFramesB, calTree=calTree)
	print "rsrfCal"
	slicedFramesR = rsrfCal(slicedFramesR, calTree=calTree)
	slicedFramesB = rsrfCal(slicedFramesB, calTree=calTree)
	print "specRespCal"
	slicedFramesR = specRespCal(slicedFramesR, csResponseAndDark=csResponseAndDarkR)
	slicedFramesB = specRespCal(slicedFramesB, csResponseAndDark=csResponseAndDarkB)
#	slicedFramesR = specRespCal(slicedFramesR, calTree=calTree)
#	slicedFramesB = specRespCal(slicedFramesB, calTree=calTree)	
	print "NHSC transient correction"
	slicedFramesR2 = getSlicedCopy(slicedFramesR)
	if verbose:
		slicedSummary(slicedFramesR2)
		slice = 3
		p5 = plotSignalBasic(slicedFramesR2, slice=slice)
		slice = 2
		p5off = plotSignalBasic(slicedFramesR2, slice=slice)
		module = 12
		p5 = plotTransient(slicedFramesR2, slice=slice, module=module, color=java.awt.Color.black, title="Slice "+str(slice)+" Module = "+str(module))
	slicedFramesR2 = specLongTermTransient(slicedFramesR2)
	if verbose:
		slicedSummary(slicedFramesR2)
		slice = 3
		p5ton = plotSignalBasic(slicedFramesR2, slice=slice)
		slice = 2
		p5toff = plotSignalBasic(slicedFramesR2, slice=slice)
		p5 = plotTransient(slicedFramesR2, slice=slice, module=module, color=java.awt.Color.black, title="Slice "+str(slice)+" Module = "+str(module))
	slicedFramesB2 = getSlicedCopy(slicedFramesB)
	if verbose:
		slicedSummary(slicedFramesB2)
		slice = 2
		p5 = plotSignalBasic(slicedFramesB2, slice=slice)
		slice = 3
		p5off = plotSignalBasic(slicedFramesB2, slice=slice)
		p5 = plotTransient(slicedFramesB2, slice=slice, module=module, color=java.awt.Color.black, title="Slice "+str(slice)+" Module = "+str(module))
	slicedFramesB2 = specLongTermTransient(slicedFramesB2)
	if verbose:
		slicedSummary(slicedFramesB2)
		slice = 2
		p5ton = plotSignalBasic(slicedFramesB2, slice=slice)
		slice = 3
		p5toff = plotSignalBasic(slicedFramesB2, slice=slice)
		p5 = plotTransient(slicedFramesB2, slice=slice, module=module, color=java.awt.Color.black, title="Slice "+str(slice)+" Module = "+str(module))
	p=PlotXY(slicedFramesR.refs[3].product["Status"]["RESETINDEX"].data,slicedFramesR.refs[8].product["Signal"].data[2,12,:])
	p[2]=LayerXY(slicedFramesR2.refs[3].product["Status"]["RESETINDEX"].data,slicedFramesR2.refs[8].product["Signal"].data[2,12,:])
	print "END NHSC transient correction"
	slicedFramesR3 = getSlicedCopy(slicedFramesR2)
	ktrans_pos(slicedFramesR3)
	slicedFramesR3 = ktrans_frame(slicedFramesR3)
	p[3]=LayerXY(slicedFramesR3.refs[3].product["Status"]["RESETINDEX"].data,slicedFramesR3.refs[8].product["Signal"].data[2,12,:])
	slicedFramesB3 = getSlicedCopy(slicedFramesB2)
	slicedFramesB3 = ktrans_frame(slicedFramesB3)
	
#saved frames here for test purposes "ngc3077_phaseAtesting_R" etc
#name="ngc3077_phaseAtesting_R2"
#saveSlicedCopy(slicedFramesR2,name)
	slicedFramesR = selectSlices(slicedFramesR,scical="sci")	#no trans
	slicedFramesR2 = selectSlices(slicedFramesR2,scical="sci")	#Dtrans
	slicedFramesR3 = selectSlices(slicedFramesR3,scical="sci")	#ktrans
	slicedFramesB = selectSlices(slicedFramesB,scical="sci")
	slicedFramesB2 = selectSlices(slicedFramesB2,scical="sci")
	slicedFramesB3 = selectSlices(slicedFramesB3,scical="sci")
	#code from Kevin to trim unwanted data from red and blue frames based on gpr
	nslice = slicedFramesB2["MasterBlockTable"]["FramesNo"].data[slicedFramesB2["MasterBlockTable"]["FramesNo"].data.dimensions[0]-1]
	qc = 0
	for qq in range(0,nslice+1):
		gpr = slicedFramesB.refs[qq].product["Status"]["GPR"].data[0]
		if (gpr > 660000)&(gpr < 690000) | (gpr > 930000)&(gpr < 950000) | (gpr > 200000)&(gpr < 260000): del(slicedFramesB2.refs[qc]) 
		qc = qc+1
		if (gpr > 660000)&(gpr < 690000) | (gpr > 930000)&(gpr < 950000) | (gpr > 200000)&(gpr < 260000): qc=qc-1
	nslice = slicedFramesB3["MasterBlockTable"]["FramesNo"].data[slicedFramesB3["MasterBlockTable"]["FramesNo"].data.dimensions[0]-1]
	qc = 0 
	for qq in range(0,nslice+1):
		gpr = slicedFramesB.refs[qq].product["Status"]["GPR"].data[0]
		if (gpr > 660000)&(gpr < 690000) | (gpr > 930000)&(gpr < 950000) | (gpr > 200000)&(gpr < 260000): del(slicedFramesB3.refs[qc]) 
		qc = qc+1
		if (gpr > 660000)&(gpr < 690000) | (gpr > 930000)&(gpr < 950000) | (gpr > 200000)&(gpr < 260000): qc=qc-1
	nslice = slicedFramesR2["MasterBlockTable"]["FramesNo"].data[slicedFramesR2["MasterBlockTable"]["FramesNo"].data.dimensions[0]-1]
	qc = 0
	for qq in range(0,nslice+1):
		gpr = slicedFramesR.refs[qq].product["Status"]["GPR"].data[0]
		if (gpr > 380000)&(gpr < 420000) | (gpr > 500000)&(gpr < 540000): del(slicedFramesR2.refs[qc]) 
		qc = qc+1
		if (gpr > 380000)&(gpr < 420000) | (gpr > 500000)&(gpr < 540000): qc=qc-1
	nslice = slicedFramesR3["MasterBlockTable"]["FramesNo"].data[slicedFramesR3["MasterBlockTable"]["FramesNo"].data.dimensions[0]-1]
	qc = 0
	for qq in range(0,nslice+1):
		gpr = slicedFramesR.refs[qq].product["Status"]["GPR"].data[0]
		if (gpr > 380000)&(gpr < 420000) | (gpr > 500000)&(gpr < 540000): del(slicedFramesR3.refs[qc]) 
		qc = qc+1
		if (gpr > 380000)&(gpr < 420000) | (gpr > 500000)&(gpr < 540000): qc=qc-1
	if verbose:slicedSummary(slicedFramesR)
	if verbose:slicedSummary(slicedFramesR2)
	if verbose:slicedSummary(slicedFramesR3)
	if verbose:slicedSummary(slicedFramesB)
	if verbose:slicedSummary(slicedFramesB2)
	if verbose:slicedSummary(slicedFramesB3)
	del(slicedFramesR,slicedFramesB)
	#end Kevin section	
	slicedCubesRK = specFrames2PacsCube(slicedFramesR3)
	slicedCubesRD = specFrames2PacsCube(slicedFramesR2)
	slicedCubesBK = specFrames2PacsCube(slicedFramesB3)
	slicedCubesBD = specFrames2PacsCube(slicedFramesB2)
	if verbose:slicedSummary(slicedCubesRK)
	if verbose:slicedSummary(slicedCubesBK)
	if verbose:slicedSummary(slicedCubesRD)
	if verbose:slicedSummary(slicedCubesBD)
	slicedCubesRK_NOFF = getSlicedCopy(slicedCubesRK)
	slicedCubesBK_NOFF = getSlicedCopy(slicedCubesBK)
	slicedCubesRD_NOFF = getSlicedCopy(slicedCubesRD)
	slicedCubesBD_NOFF = getSlicedCopy(slicedCubesBD)
	print "Spectral Flat Fielding"
	# 1. Flag outliers and rebin
	waveGridR=wavelengthGrid(slicedCubesRK, oversample=2, upsample=1, calTree=calTree)
	slicedCubesRK = activateMasks(slicedCubesRK, String1d(["GLITCH","UNCLEANCHOP","NOISYPIXELS","RAWSATURATION","SATURATION","GRATMOVE", "BADPIXELS"]), exclusive = True)
	slicedCubesRK = specFlagOutliers(slicedCubesRK, waveGridRK, nSigma=5, nIter=2)
	slicedCubesRK = activateMasks(slicedCubesRK, String1d(["GLITCH","UNCLEANCHOP","NOISYPIXELS","RAWSATURATION","SATURATION","GRATMOVE", "OUTLIERS", "BADPIXELS"]), exclusive = True)
	slicedRebinnedCubesR = specWaveRebin(slicedCubesRK, waveGridRK)
	waveGridRD=wavelengthGrid(slicedCubesRD, oversample=2, upsample=1, calTree=calTree)
	slicedCubesRD = activateMasks(slicedCubesRD, String1d(["GLITCH","UNCLEANCHOP","NOISYPIXELS","RAWSATURATION","SATURATION","GRATMOVE", "BADPIXELS"]), exclusive = True)
	slicedCubesRD = specFlagOutliers(slicedCubesRD, waveGridRD, nSigma=5, nIter=2)
	slicedCubesRD = activateMasks(slicedCubesRD, String1d(["GLITCH","UNCLEANCHOP","NOISYPIXELS","RAWSATURATION","SATURATION","GRATMOVE", "OUTLIERS", "BADPIXELS"]), exclusive = True)
	slicedRebinnedCubesRD = specWaveRebin(slicedCubesRD, waveGridRD)
	# 2. Mask the spectral lines
	widthDetect =  2.5
	threshold   = 10.
	print "find lines..."
	slicedCubesMaskRK = slicedMaskLines(slicedCubesRK,slicedRebinnedCubesRK, widthDetect=widthDetect, widthMask=2.5, threshold=threshold, copy=1, verbose=verbose, maskType="INLINE")
	slicedCubesMaskRD = slicedMaskLines(slicedCubesRD,slicedRebinnedCubesRD, widthDetect=widthDetect, widthMask=2.5, threshold=threshold, copy=1, verbose=verbose, maskType="INLINE")
	#print "set lines..."
	#slicedCubesMaskRK = slicedMaskLines(slicedCubesRK,slicedRebinnedCubesRK, lineList=[157.7409,121.89757,205.17823,88.356,63.183705], widthDetect=widthDetect, widthMask=2.5, copy=1, verbose=verbose, maskType="INLINE")
	#slicedCubesMaskRD = slicedMaskLines(slicedCubesRD,slicedRebinnedCubesRD, lineList=[157.7409,121.89757,205.17823,88.356,63.183705], widthDetect=widthDetect, widthMask=2.5, copy=1, verbose=verbose, maskType="INLINE")
	# 3. Actual spectral flatfielding
	slopeInContinuum = 1
	slicedCubesRK = slicedSpecFlatFieldLine(slicedCubesMaskRK, scaling=1, copy=1, maxrange=[50.,230.], slopeInContinuum=slopeInContinuum)
	slicedCubesRD = slicedSpecFlatFieldLine(slicedCubesMaskRD, scaling=1, copy=1, maxrange=[50.,230.], slopeInContinuum=slopeInContinuum)
	# 4. Rename mask OUTLIERS to OUTLIERS_1 (specFlagOutliers would refuse to overwrite OUTLIERS) & deactivate mask INLINE
	slicedCubesRK.renameMask("OUTLIERS", "OUTLIERS_1")
	slicedCubesRK = deactivateMasks(slicedCubesRK, String1d(["INLINE", "OUTLIERS_1"]))
	if verbose: maskSummary(slicedCubesRK, slice=0)
	slicedCubesRD.renameMask("OUTLIERS", "OUTLIERS_1")
	slicedCubesRD = deactivateMasks(slicedCubesRD, String1d(["INLINE", "OUTLIERS_1"]))
	if verbose: maskSummary(slicedCubesRD, slice=0)
	# 5. Remove intermediate results
	del waveGridRK, slicedRebinnedCubesRK, slicedCubesMaskRK, waveGridRD, slicedRebinnedCubesRD, slicedCubesMaskRD
	# --- End of Spectral Flat Fielding
	# ------------------------------------------------------------------------------
	#         Processing Level 1.0
	# ------------------------------------------------------------------------------
	
	if verbose:slicedSummary(slicedCubesRD)
	if verbose:slicedSummary(slicedCubesRK)
	if verbose:slicedSummary(slicedCubesRK_NOFF)
	if verbose:slicedSummary(slicedCubesRD_NOFF)

	nameRK=galname+"_"+str(obsidlist[0].data[n])+"_"+cameraR+"_"+Hversion+"_"+Pversion+"_Ktrans"
	nameBK=galname+"_"+str(obsidlist[0].data[n])+"_"+cameraB+"_"+Hversion+"_"+Pversion+"_Ktrans"
	nameRD=galname+"_"+str(obsidlist[0].data[n])+"_"+cameraR+"_"+Hversion+"_"+Pversion+"_Dtrans"
	nameBD=galname+"_"+str(obsidlist[0].data[n])+"_"+cameraB+"_"+Hversion+"_"+Pversion+"_Dtrans"
	if (slicedCubesRK.getRefs().size() > 1):saveSlicedCopy(slicedCubesRK,nameRK)
	if (slicedCubesBK.getRefs().size() > 1):saveSlicedCopy(slicedCubesBK,nameBK)
	if (slicedCubesRD.getRefs().size() > 1):saveSlicedCopy(slicedCubesRD,nameRD)
	if (slicedCubesBD.getRefs().size() > 1):saveSlicedCopy(slicedCubesBD,nameBD)
	#delete products before cycling to the next galaxy
	print "finished with " + str(poollist[0].data[n])
	System.gc()
	del(slicedFramesR,slicedFramesB,slicedFramesR2,slicedFramesB2,slicedCubesR,slicedCubesB,slicedCubesR2,slicedCubesB2,nameR,nameB)
		
	# End Phase A



print "you have arrived"

print "CONGRATULATIONS! Phase A complete!"

################################################
#                                              #
#    PHASE A: KINGFISH Spectroscopic Pipeline  #
#      BETA version tested in HIPE 8.0.3067    #
#         person to blame: Kevin Croxall       #
#            (aside from the NHSC)             #
#                 Nov 14, 2011                 #
################################################

