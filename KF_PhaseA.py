################################################
#                                              #
#    PHASE A: KINGFISH Spectroscopic Pipeline  #
#            Tested in HIPE 7.0.1786           #
#         person to blame: Kevin Croxall       #
#            (aside from the NHSC)             #
#                 May 13, 2011                 #
################################################

# All data must be imported into HIPE data pools.  The pipeline will not unpack the HSA tarballs
# and import them into HIPE on its own.  A list of the pools and OBSIDs that are to be processed
# need to be saved and referenced here at the begining of the Phase A pipeline.  Data are saved 
# into new pools that are sorted by AOR and line.  The frames can be saved chopped up by raster
# or as one large data stream.

Pversion = "PhaseA_v6"
Hversion = "7.0.1786"
Cversion = "16"

poollist = simpleAsciiTableReader(file = "/home/kcroxall/poolfile4.dat")     #UPDATE to the correct file location
obsidlist = simpleAsciiTableReader(file = "/home/kcroxall/obsidfile4.dat")   #UPDATE to the correct file location
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
	verbose = 0
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
#	calTreeb = obs.calibration 			# to use the cal files packages with the obs rather than the on machine latest cal
	if verbose: 
		print calTree
		print calTree.common
		print calTree.spectrometer
	
	slicedFramesR  = SlicedFrames(level0.fitted.getCamera(cameraR).product)
	slicedRawRampR = level0.raw.getCamera(cameraR).product
	slicedDmcHeadR = level0.dmc.getCamera(cameraR).product    
	slicedFramesB  = SlicedFrames(level0.fitted.getCamera(cameraB).product)
	slicedRawRampB = level0.raw.getCamera(cameraB).product
	slicedDmcHeadB = level0.dmc.getCamera(cameraB).product    
	
	if verbose:
		slicedSummary(slicedFramesR)
		p0 = slicedSummaryPlot(slicedFramesR,signal=1)
#try :
#    hp = obs.auxiliary.refs["HorizonsProduct"].product
#except :
#    print "WARNING : No Horizons found !"
#    hp = None	

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
	
	if verbose:
		slicedSummary(slicedFramesR)
		maskSummary(slicedFramesR)
		p1 = slicedSummaryPlot(slicedFramesR,signal=0)
	
	rules = [SlicingRule("LineId",1),SlicingRule("RasterLineNum",1),SlicingRule("RasterColumnNum",1),SlicingRule("NoddingPosition",1),SlicingRule("NodCycleNum",1),SlicingRule("IsOutOfField",1),SlicingRule("Band",1)]
	slicedFramesR = pacsSliceContext(slicedFramesR, slicingRules = rules, removeUndefined=1)
	slicedFramesR = specAddGratingCycleStatus(slicedFramesR, copy = True)
	slicedFramesB= pacsSliceContext(slicedFramesB, slicingRules = rules, removeUndefined=1)
	slicedFramesB = specAddGratingCycleStatus(slicedFramesB, copy = True)
#	slicedFramesR = addIndexInCycle(slicedFramesR) old tasks...
#	slicedFramesB = addIndexInCycle(slicedFramesB)
	if verbose:
		slicedSummary(slicedFramesR)
		p2 = slicedSummaryPlot(slicedFramesR,signal=0)
		p2 = slicedSummaryPlot(slicedFramesB,signal=0)
	addSliceMetaData(slicedFramesR)
	addSliceMetaData(slicedFramesB)
#	slicedDmcHeadR = pacsSliceContextByReference(slicedDmcHeadR,slicedFramesR)     ??
#	slicedDmcHeadB = pacsSliceContextByReference(slicedDmcHeadB,slicedFramesB)     ??
	if updateObservationContext:
	  obs = updateObservation(obs, cameraR, "0.5", slicedFrames=slicedFramesR)
	  obs = updateObservation(obs, cameraB, "0.5", slicedFrames=slicedFramesB)
#	slicedFramesR = cleanupSlicedFrames(slicedFramesR)
#	slicedFramesB = cleanupSlicedFrames(slicedFramesB)
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

	
##############  psuedo BLM
#OI
#for qq in range (1,slicedFramesB2.getRefs().size()):
#	for i in range (0,18):
#		for j in range (0,25):
#			reset = slicedFramesB2.refs[qq].product["Status"]["RESETINDEX"].data
#			gpr = slicedFramesB2.refs[qq].product["Status"]["GPR"].data
#			flux = slicedFramesB2.refs[qq].product["Signal"].data
#			ind = gpr.where(((gpr < 400000).and(gpr > 397740)).or((gpr > 387200).and(gpr < 389180)))
#			flux[i,j,ind] = float('NaN')
#	slicedFramesB2.refs[qq].product["Signal"].data = flux
#
#p1=PlotXY(slicedFramesB2.refs[3].product["Wave"].data[8,8,:],slicedFramesB2.refs[qq].product["Signal"].data[8,8,:])
##CII,NII122
#p1=PlotXY(slicedFramesR2.refs[2].product["Wave"].data[8,8,:],slicedFramesR2.refs[2].product["Signal"].data[8,8,:],line=0)
#
#for qq in range (1,slicedFramesR2.getRefs().size()):
#	for i in range (0,18):
#		for j in range (0,25):
#			reset = slicedFramesR2.refs[qq].product["Status"]["RESETINDEX"].data
#			gpr = slicedFramesR2.refs[qq].product["Status"]["GPR"].data
#			flux = slicedFramesR2.refs[qq].product["Signal"].data
#			ind = gpr.where(((gpr < 937400).and(gpr > 934600)).or((gpr > 919400).and(gpr < 922200)).or((gpr < 679000).and(gpr > 676000)).or((gpr < 663800).and(gpr > 661000)))
#			flux[i,j,ind] = float('NaN')
#	slicedFramesR2.refs[qq].product["Signal"].data = flux
#
#p1[1]=LayerXY(slicedFramesR2.refs[2].product["Wave"].data[8,8,:],slicedFramesR2.refs[2].product["Signal"].data[8,8,:],line=0)
##############  END psuedo BLM
	
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
	slicedFramesB2 = activateMasks(slicedFramesB2, slicedFramesB2.get(0).getMaskTypes())
	slicedFramesR2 = addQualityInformation(slicedFramesR2)
	slicedFramesB2 = addQualityInformation(slicedFramesB2)
	
	slicedFramesR2 = activateMasks(slicedFramesR2, slicedFramesR2.get(0).getMaskTypes())
	slicedFramesR2 = convertSignal2StandardCap(slicedFramesR2, calTree=calTree)
	slicedFramesB2 = activateMasks(slicedFramesB2, slicedFramesB2.get(0).getMaskTypes())
	slicedFramesB2 = convertSignal2StandardCap(slicedFramesB2, calTree=calTree)
	
	calFrameR2 = activateMasks(slicedFramesR2.getCal(0), slicedFramesR2.getCal(0).getMaskTypes())
	csResponseAndDarkR2 = specDiffCs(calFrameR2, calTree = calTree)
	slicedFramesR2 = specSubtractDark(slicedFramesR2, calTree = calTree, scical = "sci", keepall = False)
	calFrameB2 = activateMasks(slicedFramesB2.getCal(0), slicedFramesB2.getCal(0).getMaskTypes())
	csResponseAndDarkB2 = specDiffCs(calFrameB2, calTree = calTree)
	slicedFramesB2 = specSubtractDark(slicedFramesB2, calTree = calTree, scical = "sci", keepall = False)
	
#	slicedFramesR2 = selectSlices(slicedFramesR2,scical="sci")
#	slicedFramesR2 = specSubtractDark(slicedFramesR2, calTree=calTree)
#	slicedFramesB2 = selectSlices(slicedFramesB2,scical="sci")
#	slicedFramesB2 = specSubtractDark(slicedFramesB2, calTree=calTree)
	
	slicedFramesR2 = rsrfCal(slicedFramesR2, calTree=calTree)
	slicedFramesB2 = rsrfCal(slicedFramesB2, calTree=calTree)
	slicedFramesR2 = specRespCal(slicedFramesR2, calTree=calTree)
	slicedFramesB2 = specRespCal(slicedFramesB2, calTree=calTree)
	if verbose:
		slicedSummary(slicedFramesR2)
		slice = 0
		p5 = plotSignalBasic(slicedFramesR2, slice=slice)
		slice = 1
		p5off = plotSignalBasic(slicedFramesR2, slice=slice)
	#p1=PlotXY(slicedFramesR2.refs[1].product["Status"]["RESETINDEX"].data,slicedFramesR2.refs[1].product["Signal"].data[8,15,:],line=0, symbol = 14, symbolSize = 1)
	slicedFramesR3 = getSlicedCopy(slicedFramesR2)
#	slicedFramesR3 = specLongTermTransient(slicedFramesR3)
	#p1=PlotXY(slicedFramesR3.refs[1].product["Status"]["RESETINDEX"].data,slicedFramesR3.refs[1].product["Signal"].data[8,15,:],line=0, symbol = 14, symbolSize = 1)
	slicedFramesB3 = getSlicedCopy(slicedFramesB2)
#	slicedFramesB3 = specLongTermTransient(slicedFramesB3)
	
#	slicedFramesR3 = specFlatFieldRange(slicedFramesR3,polyOrder=5, minWaveRangeForPoly=4., verbose=1)
#	slicedFramesB3 = specFlatFieldRange(slicedFramesB3,polyOrder=5, minWaveRangeForPoly=4., verbose=1)
#####################
#   We could save here before projecting... Apply Off subtraction here?
#   Pros: simpler extraction of every pixel
#####################
	print "Frames -> Cubes"
	slicedCubesR = specFrames2PacsCube(slicedFramesR3)
	slicedCubesB = specFrames2PacsCube(slicedFramesB3)
	slicedCubesR.meta["Pversion"]=StringParameter(Pversion)
	slicedCubesR.meta["Hversion"]=StringParameter(Hversion)
	slicedCubesR.meta["Cversion"]=StringParameter(Cversion)
	slicedCubesB.meta["Pversion"]=StringParameter(Pversion)
	slicedCubesB.meta["Hversion"]=StringParameter(Hversion)
	slicedCubesB.meta["Cversion"]=StringParameter(Cversion)
	nameR='n1266_noBLMsim_red_phaseA'
	nameB='n1266_noBLMsim_blue_phaseA'
#	nameR=str(obsidlist[0].data[n])+"_"+cameraR+"_pipe6_PhaseA" 
#	nameB=str(obsidlist[0].data[n])+"_"+cameraB+"_pipe6_PhaseA" 
	saveSlicedCopy(slicedCubesR,nameR)
	saveSlicedCopy(slicedCubesB,nameB)
#	saveSlicedCopy(slicedFramesR3,nameR)
#	saveSlicedCopy(slicedFramesB3,nameB)
	#delete products before cycling to the next galaxy
	print "finished with " + str(poollist[0].data[n])
	System.gc()
	del(gpr,level0,nslice,qc,qq,slicedFramesR2,slicedFramesB2,slicedFramesR3,slicedFramesB3,slicedCubesR,slicedCubesB,nameR,nameB,flux,i,j,reset,slicedSpecFlagOutliers,slicedSpecWaveRebin,slicedWavelengthGrid,calFrameB2,cakFrameR2,csResponseAndDarkB2,csResponseAndDarkR2)
		
	# End Phase A



print "you have arrived"

print "CONGRATULATIONS! Phase A complete!"

################################################
#                                              #
#    PHASE A: KINGFISH Spectroscopic Pipeline  #
#            Tested in HIPE 7.0.1786           #
#         person to blame: Kevin Croxall       #
#            (aside from the NHSC)             #
#                 May 13, 2011                 #
################################################

