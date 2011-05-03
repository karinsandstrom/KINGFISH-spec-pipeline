################################################
#                                              #
#    PHASE A: KINGFISH Spectroscopic Pipeline  #
#       BETA version tested in HIPE 7.0.815    #
#         person to blame: Kevin Croxall       #
#            (aside from the NHSC)             #
#                 Mar 28, 2011                 #
################################################

# All data must be imported into HIPE data pools.  The pipeline will not unpack the HSA tarballs
# and import them into HIPE on its own.  A list of the pools and OBSIDs that are to be processed
# need to be saved and referenced here at the begining of the Phase A pipeline.  Data are saved 
# into new pools that are sorted by AOR and line.  The frames can be saved chopped up by raster
# or as one large data stream.

Pversion = "PhaseA_v4"
Hversion = "7.0.815"
Cversion = "16"

poollist = simpleAsciiTableReader(file = "/home/kcroxall/poolfile4.dat")         #UPDATE to the correct file location
obsidlist = simpleAsciiTableReader(file = "/home/kcroxall/obsidfile4.dat")   #UPDATE to the correct file location
ndim = poollist[0].data.dimensions[0]

from herschel.pacs.signal import MaskViewer

for n in range(0,ndim):
	# ------------------------------------------------------------------------------
	#        Preparation
	# ------------------------------------------------------------------------------
	print poollist[0].data[n]
	print obsidlist[0].data[n]
	print n,poollist[0].data.dimensions[0]
	useHsa = 0   #change to 1 to download....
	obs    = getObservation(obsidlist[0].data[n], verbose=True, useHsa=useHsa, poolLocation=None, poolName=str(poollist[0].data[n]))
	if useHsa: obs = saveObservation(obs, poolLocation=None, poolName=poollist[0].data[n])
# verbose: 0 - silent, execute the pipeline only
#	       1 - will trigger diagnostic output on the screen, plots, and displays
	verbose = 0
# updateObservationContext
#      0 - do not update the observation context
#	   1 - update the observation context
	updateObservationContext = 0
	if verbose: obsSummary(obs)
	cameraR = "red"
	cameraB = "blue"
	calTree = getCalTree(obs=obs)
#	calTreeb = obs.calibration 			# to use the cal files packages with the obs rather than the on machine latest cal
	if verbose: 
		print calTree
		print calTree.common
		print calTree.spectrometer
	
	pacsPropagateMetaKeywords(obs,'0', obs.level0)
	level0 = PacsContext(obs.level0)
	slicedFramesR  = SlicedFrames(level0.fitted.getCamera(cameraR).product)
	slicedRawRampR = level0.raw.getCamera(cameraR).product
	slicedDmcHeadR = level0.dmc.getCamera(cameraR).product    
	slicedFramesB  = SlicedFrames(level0.fitted.getCamera(cameraB).product)
	slicedRawRampB = level0.raw.getCamera(cameraB).product
	slicedDmcHeadB = level0.dmc.getCamera(cameraB).product    
	
	if verbose:
		slicedSummary(slicedFramesR)
		p0 = slicedSummaryPlot(slicedFramesR,signal=1)
	
# ------------------------------------------------------------------------------
#        Processing      Level 0 -> Level 0.5
# ------------------------------------------------------------------------------
	print "Level 0 -> Level 0.5"
	slicedFramesR = specFlagSaturationFrames(slicedFramesR, rawRamp = slicedRawRampR, calTree=calTree, copy=1)
	slicedFramesR = specConvDigit2VoltsPerSecFrames(slicedFramesR, calTree=calTree)
	slicedFramesR = detectCalibrationBlock(slicedFramesR)
	slicedFramesR = addUtc(slicedFramesR, obs.auxiliary.timeCorrelation)
	slicedFramesR = specAddInstantPointing(slicedFramesR, obs.auxiliary.pointing, calTree = calTree, orbitEphem = obs.auxiliary.orbitEphemeris, horizonsProduct = None)    
	slicedFramesR = specExtendStatus(slicedFramesR, calTree=calTree)
	slicedFramesR = convertChopper2Angle(slicedFramesR, calTree=calTree)
	slicedFramesR = specAssignRaDec(slicedFramesR, calTree=calTree)
	slicedFramesB = specFlagSaturationFrames(slicedFramesB, rawRamp = slicedRawRampB, calTree=calTree, copy=1)
	slicedFramesB = specConvDigit2VoltsPerSecFrames(slicedFramesB, calTree=calTree)
	slicedFramesB = detectCalibrationBlock(slicedFramesB)
	slicedFramesB = addUtc(slicedFramesB, obs.auxiliary.timeCorrelation)
	slicedFramesB = specAddInstantPointing(slicedFramesB, obs.auxiliary.pointing, calTree = calTree, orbitEphem = obs.auxiliary.orbitEphemeris, horizonsProduct = None)    
	slicedFramesB = specExtendStatus(slicedFramesB, calTree=calTree)
	slicedFramesB = convertChopper2Angle(slicedFramesB, calTree=calTree)
	slicedFramesB = specAssignRaDec(slicedFramesB, calTree=calTree)
	
	if verbose: 
		slicedPlotPointingOnOff(slicedFramesR)
	
	slicedFramesR = waveCalc(slicedFramesR, calTree=calTree)
	slicedFramesR = specCorrectHerschelVelocity(slicedFramesR, obs.auxiliary.orbitEphemeris, obs.auxiliary.pointing, obs.auxiliary.timeCorrelation)
	slicedFramesR = findBlocks(slicedFramesR, calTree = calTree)
	slicedFramesR = specFlagBadPixelsFrames(slicedFramesR, calTree=calTree)
	slicedFramesR = flagChopMoveFrames(slicedFramesR, dmcHead=slicedDmcHeadR, calTree=calTree)
	slicedFramesR = flagGratMoveFrames(slicedFramesR, dmcHead=slicedDmcHeadR, calTree=calTree)
	slicedFramesB = waveCalc(slicedFramesB, calTree=calTree)
	slicedFramesB = specCorrectHerschelVelocity(slicedFramesB, obs.auxiliary.orbitEphemeris, obs.auxiliary.pointing, obs.auxiliary.timeCorrelation)
	slicedFramesB = findBlocks(slicedFramesB, calTree = calTree)
	slicedFramesB = specFlagBadPixelsFrames(slicedFramesB, calTree=calTree)
	slicedFramesB = flagChopMoveFrames(slicedFramesB, dmcHead=slicedDmcHeadB, calTree=calTree)
	slicedFramesB = flagGratMoveFrames(slicedFramesB, dmcHead=slicedDmcHeadB, calTree=calTree)
	
	if verbose:
		slicedSummary(slicedFramesR)
		maskSummary(slicedFramesR)
		p1 = slicedSummaryPlot(slicedFramesR,signal=0)
	
	rules = [SlicingRule("LineId",1),SlicingRule("RasterLineNum",1),SlicingRule("RasterColumnNum",1),SlicingRule("NoddingPosition",1),SlicingRule("NodCycleNum",1),SlicingRule("IsOutOfField",1),SlicingRule("Band",1)]
	slicedFramesR = pacsSliceContext(slicedFramesR, slicingRules = rules, removeUndefined=1)
	slicedFramesR = addIndexInCycle(slicedFramesR)
	slicedFramesB= pacsSliceContext(slicedFramesB, slicingRules = rules, removeUndefined=1)
	slicedFramesB = addIndexInCycle(slicedFramesB)
	
	if verbose:
		slicedSummary(slicedFramesR)
		p2 = slicedSummaryPlot(slicedFramesR,signal=0)
	
	addSliceMetaData(slicedFramesR)
	addSliceMetaData(slicedFramesB)
	if updateObservationContext:
	  obs = updateObservation(obs, cameraR, "0.5", slicedFrames=slicedFramesR)
	  obs = updateObservation(obs, cameraB, "0.5", slicedFrames=slicedFramesB)
	print "Delete superfluous products to ease RAM usage"
	del(slicedRawRampR,slicedRawRampB,useHsa,slicedDmcHeadR,slicedDmcHeadB)
	System.gc()
# ------------------------------------------------------------------------------
#         Processing      Level 0.5 -> Level 1
# ------------------------------------------------------------------------------
	print "Level 0.5 -> Level 1"
	#code from Kevin to trim unwanted data from red and blue frames based on gpr
	slicedFramesR2 = getSlicedCopy(slicedFramesR)
	slicedFramesB2 = getSlicedCopy(slicedFramesB)
	
	nslice = slicedFramesB2["MasterBlockTable"]["FramesNo"].data[slicedFramesB2["MasterBlockTable"]["FramesNo"].data.dimensions[0]-1]
	qc =1 
	for qq in range(1,nslice+1):
		gpr = slicedFramesB.refs[qq].product["Status"]["GPR"].data[0]
		if (gpr > 660000)&(gpr < 690000) | (gpr > 930000)&(gpr < 950000) | (gpr > 200000)&(gpr < 260000): del(slicedFramesB2.refs[qc]) 
		qc = qc+1
		if (gpr > 660000)&(gpr < 690000) | (gpr > 930000)&(gpr < 950000) | (gpr > 200000)&(gpr < 260000): qc=qc-1
	
	
	if verbose:slicedSummary(slicedFramesB2)
	
	nslice = slicedFramesR2["MasterBlockTable"]["FramesNo"].data[slicedFramesR2["MasterBlockTable"]["FramesNo"].data.dimensions[0]-1]
	qc =1 
	for qq in range(1,nslice+1):
		gpr = slicedFramesR.refs[qq].product["Status"]["GPR"].data[0]
		if (gpr > 380000)&(gpr < 420000) | (gpr > 500000)&(gpr < 540000): del(slicedFramesR2.refs[qc]) 
		qc = qc+1
		if (gpr > 380000)&(gpr < 420000) | (gpr > 500000)&(gpr < 540000): qc=qc-1
	
	if verbose:slicedSummary(slicedFramesR2)
	del(slicedFramesR,slicedFramesB)
	System.gc()
	#end Kevin section	
	if verbose:	maskSummary(slicedFramesR2)
	slicedFramesR2 = activateMasks(slicedFramesR2, String1d([" "]), exclusive = True)
	slicedFramesB2 = activateMasks(slicedFramesB2, String1d([" "]), exclusive = True)
	if verbose: maskSummary(slicedFramesR2,slice=0)
	slicedFramesR2 = specFlagGlitchFramesQTest(slicedFramesR2, copy=1)
	slicedFramesB2 = specFlagGlitchFramesQTest(slicedFramesB2, copy=1)
	if verbose:
		slicedSummary(slicedFramesR2)
		p3 = slicedSummaryPlot(slicedFramesR2,signal=0)
		slice = 1
		p4 = plotSignalBasic(slicedFramesR2, slice=slice)
		MaskViewer(slicedFramesR2.get(slice))
	
	slicedFramesR2 = activateMasks(slicedFramesR2, slicedFramesR2.get(0).getMaskTypes())
	slicedFramesR2 = convertSignal2StandardCap(slicedFramesR2, calTree=calTree)
	slicedFramesR2 = selectSlices(slicedFramesR2,scical="sci")
	slicedFramesR2 = specSubtractDark(slicedFramesR2, calTree=calTree)
	slicedFramesR2 = rsrfCal(slicedFramesR2, calTree=calTree)
	slicedFramesR2 = specRespCal(slicedFramesR2, calTree=calTree)
	
	slicedFramesB2 = activateMasks(slicedFramesB2, slicedFramesB2.get(0).getMaskTypes())
	slicedFramesB2 = convertSignal2StandardCap(slicedFramesB2, calTree=calTree)
	slicedFramesB2 = selectSlices(slicedFramesB2,scical="sci")
	slicedFramesB2 = specSubtractDark(slicedFramesB2, calTree=calTree)
	slicedFramesB2 = rsrfCal(slicedFramesB2, calTree=calTree)
	slicedFramesB2 = specRespCal(slicedFramesB2, calTree=calTree)
	if verbose:
		slicedSummary(slicedFramesR2)
		slice = 0
		p5 = plotSignalBasic(slicedFramesR2, slice=slice)
		slice = 1
		p5off = plotSignalBasic(slicedFramesR2, slice=slice)
	
	slicedFramesR2 = specFlatFieldRange(slicedFramesR2,polyOrder=5, minWaveRangeForPoly=4., verbose=1)
	slicedFramesB2 = specFlatFieldRange(slicedFramesB2,polyOrder=5, minWaveRangeForPoly=4., verbose=1)
	nslicesR = len(slicedFramesR2.refs)
	raR = Double1d(nslicesR)
	decR = Double1d(nslicesR)
	
	nslicesB = len(slicedFramesB2.refs)
	raB = Double1d(nslicesB)
	decB = Double1d(nslicesB)
	
	for i in range(nslicesR):
	    frame = slicedFramesR2.refs[i].product
	    raR[i]=MEDIAN(frame.ra[8,12,:])
	    decR[i]=MEDIAN(frame.dec[8,12,:])
	
	for i in range(nslicesB):
	    frame = slicedFramesB2.refs[i].product
	    raB[i]=MEDIAN(frame.ra[8,12,:])
	    decB[i]=MEDIAN(frame.dec[8,12,:])
	
	
	print "Frames -> Cubes"
	if verbose: PlotXY(raR)
	if verbose: PlotXY(decR)
	slicedCubesR = specFrames2PacsCube(slicedFramesR2)
	slicedCubesB = specFrames2PacsCube(slicedFramesB2)
	slicedCubesR.meta["Pversion"]=StringParameter(Pversion)
	slicedCubesR.meta["Hversion"]=StringParameter(Hversion)
	slicedCubesR.meta["Cversion"]=StringParameter(Cversion)
	slicedCubesB.meta["Pversion"]=StringParameter(Pversion)
	slicedCubesB.meta["Hversion"]=StringParameter(Hversion)
	slicedCubesB.meta["Cversion"]=StringParameter(Cversion)
	nameR=str(poollist[0].data[n])+"_"+cameraR+"_pipe_PhaseA" 
	nameB=str(poollist[0].data[n])+"_"+cameraB+"_pipe_PhaseA" 
	
	saveSlicedCopy(slicedCubesR,nameR)
	saveSlicedCopy(slicedCubesB,nameB)
	#delete products before cycling to the next galaxy
	print "finished with " + str(poollist[0].data[n])
	System.gc()
	del(decB,decR,gpr,frame,level0,nslice,nslicesB,nslicesR,qc,qq,raB,raR,slicedFramesR2,slicedFramesB2)
		
	# End Phase A



print "you have arrived"

print "CONGRATULATIONS! Phase A complete!"

################################################
#                                              #
#    PHASE A: KINGFISH Spectroscopic Pipeline  #
#       BETA version tested in HIPE 7.0.815    #
#         person to blame: Kevin Croxall       #
#            (aside from the NHSC)             #
#                 Mar 28, 2011                 #
################################################

