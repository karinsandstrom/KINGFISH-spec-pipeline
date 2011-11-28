################################################
#                                              #
#    PHASE A: KINGFISH Spectroscopic Pipeline  #
#      BETA version tested in HIPE 8.0.3215    #
#         person to blame: Kevin Croxall       #
#            (aside from the NHSC)             #
#                 Nov 28, 2011                 #
################################################

# All data must be imported into HIPE data pools.  The pipeline will not unpack the HSA tarballs
# and import them into HIPE on its own.  A list of the pools and OBSIDs that are to be processed
# need to be saved and referenced here at the begining of the Phase A pipeline.  Data are saved 
# into new pools that are sorted by AOR and line.  The frames can be saved chopped up by raster
# or as one large data stream.

Pversion = "PhaseA_v3.0"
Hversion = "8.0.3215"
Cversion = "32"
shortver = "v30"
poollist = simpleAsciiTableReader(file = "/Users/kcroxall/poolfile3077.dat")     #UPDATE to the correct file location
obsidlist = simpleAsciiTableReader(file = "/Users/kcroxall/obsidfile3077.dat")   #UPDATE to the correct file location

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
	verbose = 1
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
	slicedDmcHeadR = level0.dmc.getCamera(cameraR).product    			#	the mechanisms' status information
	slicedFramesB  = SlicedFrames(level0.fitted.getCamera(cameraB).product)
	slicedRawRampB = level0.raw.getCamera(cameraB).product
	slicedDmcHeadB = level0.dmc.getCamera(cameraB).product    
	slicedFramesR.meta["Pversion"]=StringParameter(Pversion)
	slicedFramesR.meta["Hversion"]=StringParameter(Hversion)
	slicedFramesR.meta["Cversion"]=StringParameter(Cversion)
	slicedFramesB.meta["Pversion"]=StringParameter(Pversion)
	slicedFramesB.meta["Hversion"]=StringParameter(Hversion)
	slicedFramesB.meta["Cversion"]=StringParameter(Cversion)
	galname = slicedFramesR.meta["object"].string
	
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
	print "save PhaseA products"
	nameR=galname+"_"+str(obsidlist[0].data[n])+"_"+cameraR+"_"+shortver
	nameB=galname+"_"+str(obsidlist[0].data[n])+"_"+cameraB+"_"+shortver
	if (slicedFramesR.getRefs().size() > 1):saveSlicedCopy(slicedFramesR,nameR)
	if (slicedFramesB.getRefs().size() > 1):saveSlicedCopy(slicedFramesB,nameB)
	#delete products before cycling to the next galaxy
	print "finished with " + str(poollist[0].data[n])
	System.gc()
	del(slicedFramesR,slicedFramesB,nameR,nameB)
	# End Phase A



print "you have arrived"

print "CONGRATULATIONS! Phase A complete!"

################################################
#                                              #
#    PHASE A: KINGFISH Spectroscopic Pipeline  #
#      BETA version tested in HIPE 8.0.3215    #
#         person to blame: Kevin Croxall       #
#            (aside from the NHSC)             #
#                 Nov 28, 2011                 #
################################################

