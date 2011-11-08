################################################
#                                              #
#    PHASE A: KINGFISH Spectroscopic Pipeline  #
#      BETA version tested in HIPE 8.0.3067    #
#         person to blame: Kevin Croxall       #
#            (aside from the NHSC)             #
#                 Nov  8, 2011                 #
################################################

# All data must be imported into HIPE data pools.  The pipeline will not unpack the HSA tarballs
# and import them into HIPE on its own.  A list of the pools and OBSIDs that are to be processed
# need to be saved and referenced here at the begining of the Phase A pipeline.  Data are saved 
# into new pools that are sorted by AOR and line.  The frames can be saved chopped up by raster
# or as one large data stream.

Pversion = "PhaseA_v1.6"
Hversion = "8.0.3067"
Cversion = "32"

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
	print "tmp"
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
	
	slicedFramesB2 = specLongTermTransient(slicedFramesB2)
	
	if verbose:
		slicedSummary(slicedFramesB2)
		slice = 2
		p5ton = plotSignalBasic(slicedFramesB2, slice=slice)
		slice = 3
		p5toff = plotSignalBasic(slicedFramesB2, slice=slice)
	
	p=PlotXY(slicedFramesR.refs[3].product["Status"]["RESETINDEX"].data,slicedFramesR.refs[3].product["Signal"].data[2,12,:])
	p[2]=LayerXY(slicedFramesR2.refs[3].product["Status"]["RESETINDEX"].data,slicedFramesR2.refs[3].product["Signal"].data[2,12,:])
	print "Level 0.5 -> Level 1"
	print "END NHSC transient correction"
	slicedFramesR = selectSlices(slicedFramesR,scical="sci")
	slicedFramesB = selectSlices(slicedFramesB,scical="sci")
	slicedFramesR2 = selectSlices(slicedFramesR2,scical="sci")
	slicedFramesB2 = selectSlices(slicedFramesB2,scical="sci")
	
	#code from Kevin to trim unwanted data from red and blue frames based on gpr
	slicedFramesR3 = getSlicedCopy(slicedFramesR)
	slicedFramesB3 = getSlicedCopy(slicedFramesB)
	nslice = slicedFramesB2["MasterBlockTable"]["FramesNo"].data[slicedFramesB2["MasterBlockTable"]["FramesNo"].data.dimensions[0]-1]
	qc =1 
	for qq in range(1,nslice+1):
		gpr = slicedFramesB.refs[qq].product["Status"]["GPR"].data[0]
		if (gpr > 660000)&(gpr < 690000) | (gpr > 930000)&(gpr < 950000) | (gpr > 200000)&(gpr < 260000): del(slicedFramesB2.refs[qc]) 
		qc = qc+1
		if (gpr > 660000)&(gpr < 690000) | (gpr > 930000)&(gpr < 950000) | (gpr > 200000)&(gpr < 260000): qc=qc-1
	
	nslice = slicedFramesB3["MasterBlockTable"]["FramesNo"].data[slicedFramesB3["MasterBlockTable"]["FramesNo"].data.dimensions[0]-1]
	qc =1 
	for qq in range(1,nslice+1):
		gpr = slicedFramesB.refs[qq].product["Status"]["GPR"].data[0]
		if (gpr > 660000)&(gpr < 690000) | (gpr > 930000)&(gpr < 950000) | (gpr > 200000)&(gpr < 260000): del(slicedFramesB3.refs[qc]) 
		qc = qc+1
		if (gpr > 660000)&(gpr < 690000) | (gpr > 930000)&(gpr < 950000) | (gpr > 200000)&(gpr < 260000): qc=qc-1
	
	nslice = slicedFramesR2["MasterBlockTable"]["FramesNo"].data[slicedFramesR2["MasterBlockTable"]["FramesNo"].data.dimensions[0]-1]
	qc =1 
	for qq in range(1,nslice+1):
		gpr = slicedFramesR.refs[qq].product["Status"]["GPR"].data[0]
		if (gpr > 380000)&(gpr < 420000) | (gpr > 500000)&(gpr < 540000): del(slicedFramesR2.refs[qc]) 
		qc = qc+1
		if (gpr > 380000)&(gpr < 420000) | (gpr > 500000)&(gpr < 540000): qc=qc-1
	
	nslice = slicedFramesR3["MasterBlockTable"]["FramesNo"].data[slicedFramesR3["MasterBlockTable"]["FramesNo"].data.dimensions[0]-1]
	qc =1 
	for qq in range(1,nslice+1):
		gpr = slicedFramesR.refs[qq].product["Status"]["GPR"].data[0]
		if (gpr > 380000)&(gpr < 420000) | (gpr > 500000)&(gpr < 540000): del(slicedFramesR3.refs[qc]) 
		qc = qc+1
		if (gpr > 380000)&(gpr < 420000) | (gpr > 500000)&(gpr < 540000): qc=qc-1
	
	if verbose:slicedSummary(slicedFramesR2)
	if verbose:slicedSummary(slicedFramesB2)
	if verbose:slicedSummary(slicedFramesR3)
	if verbose:slicedSummary(slicedFramesB3)
	del(slicedFramesR,slicedFramesB)
	#end Kevin section	
	
	print "Refine the spectral flatfield"
	slicedFramesR = specFlatFieldRange(slicedFramesR3,polyOrder=5, minWaveRangeForPoly=4., verbose=1)
	slicedFramesB = specFlatFieldRange(slicedFramesB3,polyOrder=5, minWaveRangeForPoly=4., verbose=1)
	slicedFramesR2 = specFlatFieldRange(slicedFramesR2,polyOrder=5, minWaveRangeForPoly=4., verbose=1)
	slicedFramesB2 = specFlatFieldRange(slicedFramesB2,polyOrder=5, minWaveRangeForPoly=4., verbose=1)
	
#####################
#   We could save here before projecting... Apply Off subtraction here?
#   Pros: simpler extraction of every pixel
#####################
	
	print "Frames -> Cubes"
	slicedCubesR = specFrames2PacsCube(slicedFramesR)
	slicedCubesB = specFrames2PacsCube(slicedFramesB)
	slicedCubesR2 = specFrames2PacsCube(slicedFramesR2)
	slicedCubesB2 = specFrames2PacsCube(slicedFramesB2)
	
	if verbose:slicedSummary(slicedCubesR)
	if verbose:slicedSummary(slicedCubesB)
	if verbose:slicedSummary(slicedCubesR2)
	if verbose:slicedSummary(slicedCubesB2)
	
	slicedCubesR.meta["Pversion"]=StringParameter(Pversion)
	slicedCubesR.meta["Hversion"]=StringParameter(Hversion)
	slicedCubesR.meta["Cversion"]=StringParameter(Cversion)
	slicedCubesB.meta["Pversion"]=StringParameter(Pversion)
	slicedCubesB.meta["Hversion"]=StringParameter(Hversion)
	slicedCubesB.meta["Cversion"]=StringParameter(Cversion)
	slicedCubesR2.meta["Pversion"]=StringParameter(Pversion)
	slicedCubesR2.meta["Hversion"]=StringParameter(Hversion)
	slicedCubesR2.meta["Cversion"]=StringParameter(Cversion)
	slicedCubesB2.meta["Pversion"]=StringParameter(Pversion)
	slicedCubesB2.meta["Hversion"]=StringParameter(Hversion)
	slicedCubesB2.meta["Cversion"]=StringParameter(Cversion)
	nameR=galname+"_"+str(obsidlist[0].data[n])+"_"+cameraR+"_"+Hversion+"_"+Pversion
	nameB=galname+"_"+str(obsidlist[0].data[n])+"_"+cameraB+"_"+Hversion+"_"+Pversion
	nameR2=galname+"_"+str(obsidlist[0].data[n])+"_"+cameraR+"_"+Hversion+"_"+Pversion+"_Dtrans"
	nameB2=galname+"_"+str(obsidlist[0].data[n])+"_"+cameraB+"_"+Hversion+"_"+Pversion+"_Dtrans"
	if (slicedCubesR.getRefs().size() > 1):saveSlicedCopy(slicedCubesR,nameR)
	if (slicedCubesB.getRefs().size() > 1):saveSlicedCopy(slicedCubesB,nameB)
	if (slicedCubesR2.getRefs().size() > 1):saveSlicedCopy(slicedCubesR2,nameR2)
	if (slicedCubesB2.getRefs().size() > 1):saveSlicedCopy(slicedCubesB2,nameB2)
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
#                 Nov  8, 2011                 #
################################################

