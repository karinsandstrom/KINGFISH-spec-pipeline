################################################
#                                              #
#    PHASE A: KINGFISH Spectroscopic Pipeline  #
#            Tested in HIPE 7.0.1931           #
#         person to blame: Kevin Croxall       #
#            (aside from the NHSC)             #
#                 Jun 3, 2011                  #
################################################

# All data must be imported into HIPE data pools.  The pipeline will not unpack the HSA tarballs
# and import them into HIPE on its own.  A list of the pools and OBSIDs that are to be processed
# need to be saved and referenced here at the begining of the Phase A pipeline.  Data are saved 
# into new pools that are sorted by AOR and line.  The frames can be saved chopped up by raster
# or as one large data stream.


Pversion = "PhaseA_v0.2"
Hversion = "7.0.1931"
Cversion = "26"

poollist = simpleAsciiTableReader(file = "/Users/kcroxall/poolfile3077.dat")     #UPDATE to the correct file location
obsidlist = simpleAsciiTableReader(file = "/Users/kcroxall/obsidfile3077.dat")   #UPDATE to the correct file location
ndim = poollist[0].data.dimensions[0]

from herschel.pacs.signal import MaskViewer
from herschel.pacs.signal import SlicedFrames
from herschel.pacs.signal.context import *
from herschel.pacs.spg.common import *
from herschel.pacs.spg.spec import *
from herschel.pacs.cal import *
from herschel.ia.numeric import *
from herschel.ia.jconsole import *
from herschel.pacs.spg.pipeline import *  

for n in range(0,ndim):
	# ------------------------------------------------------------------------------
	#        Preparation
	# ------------------------------------------------------------------------------
	print poollist[0].data[n]
	print obsidlist[0].data[n]
	print n,poollist[0].data.dimensions[0]
	useHsa = 0   #change to 1 to download....
#	obs    = getObservation('1342214219', verbose=True, useHsa=useHsa, poolLocation=None, poolName='1342214219')
	obs    = getObservation(obsidlist[0].data[n], verbose=True, useHsa=useHsa, poolLocation=None, poolName=str(poollist[0].data[n]))
	if useHsa: obs = saveObservation(obs, poolLocation=None, poolName=poollist[0].data[n])
# verbose: 0 - silent, execute the pipeline only
#	       1 - will trigger diagnostic output on the screen, plots, and displays
	verbose = 1
# updateObservationContext
#      0 - do not update the observation context
#	   1 - update the observation context
	updateObservationContext = 0
	if verbose: obsSummary(obs)
	cameraR = "red"
	cameraB = "blue"
	obs = checkForAnomaly70(obs)		# Check for decMec Anomaly 70 and attach the QualityInformation WHY? WHAT DOES THIS DO?
	pacsPropagateMetaKeywords(obs,'0', obs.level0)
	level0 = PacsContext(obs.level0)
	calTree = getCalTree(obs=obs)
#	calTreeb = obs.calibration	# to use the cal files packaged with the obs rather than the latest on-machine cal products
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
#	slicedFramesR = specAddInstantPointing(slicedFramesR, obs.auxiliary.pointing, calTree = calTree, orbitEphem = obs.auxiliary.orbitEphemeris, horizonsProduct = None)		#horizons never works for me....
	slicedFramesR = specExtendStatus(slicedFramesR, calTree=calTree)
	slicedFramesR = convertChopper2Angle(slicedFramesR, calTree=calTree)
	slicedFramesR = specAssignRaDec(slicedFramesR, calTree=calTree)
	slicedFramesB = specFlagSaturationFrames(slicedFramesB, rawRamp = slicedRawRampB, calTree=calTree, copy=1)
	slicedFramesB = specConvDigit2VoltsPerSecFrames(slicedFramesB, calTree=calTree)
	slicedFramesB = detectCalibrationBlock(slicedFramesB)
	slicedFramesB = addUtc(slicedFramesB, obs.auxiliary.timeCorrelation)
	slicedFramesB = specAddInstantPointing(slicedFramesB, obs.auxiliary.pointing, calTree = calTree, orbitEphem = obs.auxiliary.orbitEphemeris, horizonsProduct = obs.auxiliary.horizons)    
#	slicedFramesB = specAddInstantPointing(slicedFramesB, obs.auxiliary.pointing, calTree = calTree, orbitEphem = obs.auxiliary.orbitEphemeris, horizonsProduct = None)		#horizons never works for me....
	slicedFramesB = specExtendStatus(slicedFramesB, calTree=calTree)
	slicedFramesB = convertChopper2Angle(slicedFramesB, calTree=calTree)
	slicedFramesB = specAssignRaDec(slicedFramesB, calTree=calTree)
	
	if verbose: 					# show footprints for selected raster position(s)
		slicedPlotPointingOnOff(slicedFramesR)
		slicedPlotPointingOnOff(slicedFramesB)
	
	slicedFramesR = waveCalc(slicedFramesR, calTree=calTree)
	slicedFramesR = specCorrectHerschelVelocity(slicedFramesR, obs.auxiliary.orbitEphemeris, obs.auxiliary.pointing, obs.auxiliary.timeCorrelation, horizonsProduct = None)
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
	
	print "Delete superfluous products to ease RAM usage"
	del(slicedRawRampR,slicedRawRampB,useHsa,slicedDmcHeadR,slicedDmcHeadB,level0)
	System.gc()
	
	if verbose:
		slicedSummary(slicedFramesR)						# an overview of the slicedFrames contents
		slicedSummary(slicedFramesB)
		maskSummary(slicedFramesR)
		maskSummary(slicedFramesB)
		p1 = slicedSummaryPlot(slicedFramesR,signal=0)		# Show the basic data structure, without the signal
		p1 = slicedSummaryPlot(slicedFramesB,signal=0)		
	
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
	addSliceMetaData(slicedFramesR)
	addSliceMetaData(slicedFramesB)
#	slicedDmcHeadR = pacsSliceContextByReference(slicedDmcHeadR,slicedFramesR)     ??
#	slicedDmcHeadB = pacsSliceContextByReference(slicedDmcHeadB,slicedFramesB)     ??
	if updateObservationContext:
	  obs = updateObservation(obs, cameraR, "0.5", slicedFrames=slicedFramesR, calTree = calTree)
	  obs = updateObservation(obs, cameraB, "0.5", slicedFrames=slicedFramesB, calTree = calTree)
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
	
	nslice = slicedFramesR2["MasterBlockTable"]["FramesNo"].data[slicedFramesR2["MasterBlockTable"]["FramesNo"].data.dimensions[0]-1]
	qc =1 
	for qq in range(1,nslice+1):
		gpr = slicedFramesR.refs[qq].product["Status"]["GPR"].data[0]
		if (gpr > 380000)&(gpr < 420000) | (gpr > 500000)&(gpr < 540000): del(slicedFramesR2.refs[qc]) 
		qc = qc+1
		if (gpr > 380000)&(gpr < 420000) | (gpr > 500000)&(gpr < 540000): qc=qc-1

	
	if verbose:slicedSummary(slicedFramesR2)
	if verbose:slicedSummary(slicedFramesB2)
	del(slicedFramesR,slicedFramesB)
	System.gc()
	#end Kevin section	
	slicedFramesR2 = activateMasks(slicedFramesR2, String1d([" "]), exclusive = True)
	slicedFramesR2 = specFlagGlitchFramesQTest(slicedFramesR2, copy=1)
	slicedFramesB2 = activateMasks(slicedFramesB2, String1d([" "]), exclusive = True)
	slicedFramesB2 = specFlagGlitchFramesQTest(slicedFramesB2, copy=1)
	if verbose:
		slicedSummary(slicedFramesR2)
		slicedSummary(slicedFramesB2)
		p3 = slicedSummaryPlot(slicedFramesR2,signal=0)
		p3 = slicedSummaryPlot(slicedFramesB2,signal=0)
		slice = 1
		p4 = plotSignalBasic(slicedFramesR2, slice=slice)	# Single pixel, signal vs wavelength: Only unmasked datapoints are plotted
		p4 = plotSignalBasic(slicedFramesB2, slice=slice)
		MaskViewer(slicedFramesR2.get(slice))			# Inspect timeline of signals and masked signals
		MaskViewer(slicedFramesB2.get(slice))
	
	slicedFramesR2 = activateMasks(slicedFramesR2, slicedFramesR2.get(0).getMaskTypes())
	slicedFramesR2 = convertSignal2StandardCap(slicedFramesR2, calTree=calTree)
	slicedFramesB2 = activateMasks(slicedFramesB2, slicedFramesB2.get(0).getMaskTypes())
	slicedFramesB2 = convertSignal2StandardCap(slicedFramesB2, calTree=calTree)
	
	calFrameR2 = activateMasks(slicedFramesR2.getCal(0), slicedFramesR2.getCal(0).getMaskTypes())
	csResponseAndDarkR2 = specDiffCs(calFrameR2, calTree = calTree)
	slicedFramesR2 = specSubtractDark(slicedFramesR2, calTree=calTree)
	calFrameB2 = activateMasks(slicedFramesB2.getCal(0), slicedFramesB2.getCal(0).getMaskTypes())
	csResponseAndDarkB2 = specDiffCs(calFrameB2, calTree = calTree)
	slicedFramesB2 = specSubtractDark(slicedFramesB2, calTree=calTree)
	
	slicedFramesR2 = rsrfCal(slicedFramesR2, calTree=calTree)
	slicedFramesB2 = rsrfCal(slicedFramesB2, calTree=calTree)
	slicedFramesR2 = specRespCal(slicedFramesR2, calTree=calTree)
	slicedFramesB2 = specRespCal(slicedFramesB2, calTree=calTree)
	
	print "transient corr"
	if verbose:
	      slice = 1
	      module = 12
	      title = "Slice "+str(slice)+" - Module "+str(module)
	      plot = plotTransient(slicedFramesR2, slice=slice, module=module, color=Color.BLUE, title=title)
	
	slicedFramesR2 = specLongTermTransient(slicedFramesR2)
	if verbose: plot = plotTransient(slicedFramesR2, p=plot, module=module, slice=slice,color=Color.RED)
	if verbose:
	      title = "Slice "+str(slice)+" - Module "+str(module)
	      plot = plotTransient(slicedFramesB2, slice=slice, module=module, color=Color.BLUE, title=title)
	
	slicedFramesB2 = specLongTermTransient(slicedFramesB2)
	if verbose: plot = plotTransient(slicedFramesB2, p=plot, module=module, slice=slice,color=Color.RED)
	
	slicedFramesR2 = selectSlices(slicedFramesR2,scical="sci")
	slicedFramesB2 = selectSlices(slicedFramesB2,scical="sci")
	if verbose:
		slicedSummary(slicedFramesR2)
		slicedSummary(slicedFramesB2)
		slice = 0					### CHECK SLICES.... commented for CN not UnC... slice 1,0 not 0,1 ???
		p5 = plotSignalBasic(slicedFramesR2, slice=slice)		# Single pixel, signal vs wavelength
		slice = 1
		p5off = plotSignalBasic(slicedFramesR2, slice=slice)
		p6 = plotSignalBasic(slicedFramesB2, slice=slice)
		p6off = plotSignalBasic(slicedFramesB2, slice=slice)
	print "Refine the spectral flatfield"
	slicedFramesR2 = specFlatFieldRange(slicedFramesR2,polyOrder=5, minWaveRangeForPoly=4., verbose=1)
	slicedFramesB2 = specFlatFieldRange(slicedFramesB2,polyOrder=5, minWaveRangeForPoly=4., verbose=1)
	
#####################
#   We could save here before projecting... Apply Off subtraction here?
#   Pros: simpler extraction of every pixel
#####################
	print "Frames -> Cubes"
	slicedCubesR = specFrames2PacsCube(slicedFramesR2)
	slicedCubesB = specFrames2PacsCube(slicedFramesB2)
	
	slicedCubesR.meta["Pversion"]=StringParameter(Pversion)
	slicedCubesR.meta["Hversion"]=StringParameter(Hversion)
	slicedCubesR.meta["Cversion"]=StringParameter(Cversion)
	slicedCubesB.meta["Pversion"]=StringParameter(Pversion)
	slicedCubesB.meta["Hversion"]=StringParameter(Hversion)
	slicedCubesB.meta["Cversion"]=StringParameter(Cversion)
	nameR=galname+"_"+str(obsidlist[0].data[n])+"_"+cameraR+"_"+Pversion
	nameB=galname+"_"+str(obsidlist[0].data[n])+"_"+cameraB+"_"+Pversion
	saveSlicedCopy(slicedCubesR,nameR)
	saveSlicedCopy(slicedCubesB,nameB)
#	saveSlicedCopy(slicedFramesR3,nameR)
#	saveSlicedCopy(slicedFramesB3,nameB)
	#delete products before cycling to the next galaxy
	print "finished with " + str(poollist[0].data[n])
	System.gc()
	del(gpr,nslice,qc,qq,slicedFramesR2,slicedFramesB2,slicedCubesR,slicedCubesB,nameR,nameB,reset,slicedSpecFlagOutliers,slicedSpecWaveRebin,slicedWavelengthGrid,calFrameB2,cakFrameR2,csResponseAndDarkB2,csResponseAndDarkR2)
		
	# End Phase A



print "you have arrived"

print "CONGRATULATIONS! Phase A complete!"

################################################
#                                              #
#    PHASE A: KINGFISH Spectroscopic Pipeline  #
#            Tested in HIPE 7.0.1931           #
#         person to blame: Kevin Croxall       #
#            (aside from the NHSC)             #
#                 Jun 3, 2011                  #
################################################

