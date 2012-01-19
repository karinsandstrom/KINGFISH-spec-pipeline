################################################
#                                              #
#    PHASE B: KINGFISH Spectroscopic Pipeline  #
#      BETA version tested in HIPE 8.0.3215    #
#         person to blame: Kevin Croxall       #
#            (aside from the NHSC)             #
#                Jan 19, 2012                  #
################################################
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
	camera = slicedFrames.meta["camera"].string
	for rasti in range (0,frame.getRefs().size()):
		for spatialpix in range (0,25):
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
	tabname = homename + galname + "_" + obsid + "_" + camera + "_postest.tab"   # UPDATE PATH
	ascii.save(tabname, table)

def ktrans_frame(frame):
	z=frame.meta["redshiftValue"].double/300000
	for rasti in range (0,frame.getRefs().size()):
		print "Raster", rasti, " of",frame.getRefs().size()-1
		wavecent = ((frame.refs[rasti].product["BlockTable"]["MinWave"].data[0] + frame.refs[rasti].product["BlockTable"]["MaxWave"].data[0])/2)
		chanwidth = 75
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
		for spatialpix in range (0,25):
			for specpix in range (1,17):
				wave = frame.refs[rasti].product["Wave"].data[specpix,spatialpix,:]
				flux = frame.refs[rasti].product["Signal"].data[specpix,spatialpix,:]
				fluxnew = frame.refs[rasti].product["Signal"].data[specpix,spatialpix,:]
				reset = frame.refs[rasti].product["Status"]["RESETINDEX"].data[:]
				continuum = getMedian(flux)
				if (flux[33] == flux[33])&(flux[53] == flux[53]):
					for pix in range (0,reset.dimensions[0]-1):
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
						fluxnew[pix] = flux[pix]-subval2#+continuum
				frame.refs[rasti].product["Signal"].data[specpix,spatialpix,:] = fluxnew[:]
		System.gc()
	del(wave,flux,fluxnew,reset,spatialpix,specpix,rasti,dist,nearby,subval,flux2,dist2,meandist2,nearby2,subval2)
	return (frame)

def ktrans_frame_posvel_perpix(frame,table):
	counter = 0
	zgal=frame.meta["redshiftValue"].double/300000
	for rasti in range (0,frame.getRefs().size()):
		print "Raster", rasti, " of",frame.getRefs().size()-1
		wavecent = ((frame.refs[rasti].product["BlockTable"]["MinWave"].data[0] + frame.refs[rasti].product["BlockTable"]["MaxWave"].data[0])/2)
		chanwidth = 75
		linewidth = .00001
		linecent = wavecent
		stepsize = 0.006241358
		rastline = .00001
		if (wavecent > 190): 
			chanwidtht = 60
			rastline = 205.17823
			linewidtht = 0.115  #inst - 0.105
			stepsize = 0.00406
		if (wavecent > 150)&(wavecent < 170): 
			chanwidtht = 42 #bol -42
			rastline = 157.7409
			linewidtht = 0.15 #bol - 0.14 # inst - 0.126
			stepsize = 0.0062
		if (wavecent > 110)&(wavecent < 130): 
			chanwidtht = 45
			rastline = 121.89757
			linewidtht = 0.14  #inst - 0.118
			stepsize = 0.0072
		if (wavecent > 75)&(wavecent < 95): 
			chanwidtht = 45
			rastline = 88.356
			linewidtht =0.042 #inst - 0.037
			stepsize = 0.002
		if (wavecent > 60)&(wavecent < 70): 
			chanwidtht = 55
			rastline =  63.183705
			linewidtht = 0.025 #inst - 0.018
			stepsize = 0.001
		print "Line at ", rastline
		for spatialpix in range (0,25):
			medvel=float(table[4].data[counter])
			minvel=float(table[5].data[counter])
			maxvel=float(table[6].data[counter])
			counter += 1
			z = medvel/300000
			if (z == 0): z=zgal		#gaurd against vel=0/NaN values in the vel files.... hopefully will go away...
			zmin = minvel/300000
			zmax = maxvel/300000
			linecent = rastline*(1+z)
			linemax = 157.7409*(1+zmax)
			linemin = 157.7409*(1+zmin)
			if (0.15 < (linemax-linemin)/2):
				ratio = ((linemax-linemin)/2) / 0.15
				linewidth = linewidtht*ratio
				chanwidth = 2*linewidth/stepsize
			else: 
				linewidth=linewidtht
				chanwidth=chanwidtht
			if (linecent != linecent):linecent=wavecent
			print "Spatial Pixel ",spatialpix," LineCenter ",linecent," FWHM ",linewidth," BolattoFWHM_CII ",(linemax-linemin)/2," ChanWidth ",chanwidth," CalcChanWidth ",2*linewidth/stepsize
			for specpix in range (1,17):
				#print "Spectral Pixel ", specpix
				wave = frame.refs[rasti].product["Wave"].data[specpix,spatialpix,:]
				flux = frame.refs[rasti].product["Signal"].data[specpix,spatialpix,:]
				fluxnew = frame.refs[rasti].product["Signal"].data[specpix,spatialpix,:]
				reset = frame.refs[rasti].product["Status"]["RESETINDEX"].data[:]
				continuum = getMedian(flux)
				#print "Spatial Pix ", spatialpix, " Spectral Pixel ", specpix, " Continuum ", continuum, " raster ", rasti, " linewidth = ",linewidth 
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
						#meandist2 = getMedian(dist2c)
						nearby2 = dist2c.where(dist2<meandist2+2)
						subval2 = SUM(flux2[nearby2])/flux2[nearby2].dimensions[0]
						fluxnew[pix] = flux[pix]-subval2#+continuum
				frame.refs[rasti].product["Signal"].data[specpix,spatialpix,:] = fluxnew[:]
	return (frame)


# A list of the Phase A pools that are to be processed
# need to be saved and referenced here at the begining of the Phase A pipeline.  

phasealist = simpleAsciiTableReader(file = "/home/kcroxall/phaseA_dec22.lst") #UPDATE to the correct file location
homename = "/home/kcroxall/"							#UPDATE to the correct file location
ndim = phasealist[0].data.dimensions[0]
verbose = 0
postrans = 1
for n in range(0,ndim):
	name=str(phasealist[0].data[n]) 
	slicedFrames = readSliced(name)
	galname = slicedFrames.meta["object"].string
	obsid = str(slicedFrames.meta["obsid"].long)
	camera = slicedFrames.meta["camera"].string
	print "transient correction for ", n 
	slicedFrames2 = getSlicedCopy(slicedFrames)
	#code from Kevin to trim unwanted data from red and blue frames based on gpr
	if (camera == "SPECBLUE"):
		nslice = slicedFrames2["MasterBlockTable"]["FramesNo"].data[slicedFrames2["MasterBlockTable"]["FramesNo"].data.dimensions[0]-1]
		qc = 0
		for qq in range(0,nslice+1):
			gpr = slicedFrames.refs[qq].product["Status"]["GPR"].data[0]
			if (gpr > 660000)&(gpr < 690000) | (gpr > 930000)&(gpr < 950000) | (gpr > 200000)&(gpr < 260000): del(slicedFrames2.refs[qc]) 
			qc = qc+1
			if (gpr > 660000)&(gpr < 690000) | (gpr > 930000)&(gpr < 950000) | (gpr > 200000)&(gpr < 260000): qc=qc-1
	if (camera == "SPECRED"):
		nslice = slicedFrames2["MasterBlockTable"]["FramesNo"].data[slicedFrames2["MasterBlockTable"]["FramesNo"].data.dimensions[0]-1]
		qc = 0
		for qq in range(0,nslice+1):
			gpr = slicedFrames.refs[qq].product["Status"]["GPR"].data[0]
			if (gpr > 380000)&(gpr < 420000) | (gpr > 500000)&(gpr < 540000): del(slicedFrames2.refs[qc]) 
			qc = qc+1
			if (gpr > 380000)&(gpr < 420000) | (gpr > 500000)&(gpr < 540000): qc=qc-1
	if (slicedFrames2.getRefs().size() < 2): 
		print "I need valid data to make this work"
		continue
	slicedFrames = selectSlices(slicedFrames,scical="sci")
	slicedFrames2 = selectSlices(slicedFrames2,scical="sci")
	if verbose:slicedSummary(slicedFrames)
	if verbose:slicedSummary(slicedFrames2)
	ktrans_pos(slicedFrames2)
	slicedFrames = getSlicedCopy(slicedFrames2)
	if verbose:
		slicedSummary(slicedFrames2)
		slice = 2
		p5 = plotSignalBasic(slicedFrames2, slice=slice)
		slice = 1
		p5off = plotSignalBasic(slicedFrames2, slice=slice)
		module = 12
		p5 = plotTransient(slicedFrames2, slice=slice, module=module, color=java.awt.Color.black, title="Slice "+str(slice)+" Module = "+str(module))
	if postrans:
		from herschel.ia.io.ascii import AsciiParser
		atrt = AsciiTableReaderTask()
		arrname = homename + "tables/" + galname+"_"+obsid+"_"+camera+"_veltab.txt"
		posarr = atrt(file=arrname, parserDelim='\t', parserGuess=AsciiParser.GUESS_TRY,parserSkip=1)
		slicedFrames2 = ktrans_frame_posvel_perpix(slicedFrames2,posarr)
	else:slicedFrames2 = ktrans_frame(slicedFrames2)
	if verbose:
		slicedSummary(slicedFrames2)
		slice = 2
		p5 = plotSignalBasic(slicedFrames2, slice=slice)
		slice = 1
		p5off = plotSignalBasic(slicedFrames2, slice=slice)
		module = 12
		p5 = plotTransient(slicedFrames2, slice=slice, module=module, color=java.awt.Color.black, title="Slice "+str(slice)+" Module = "+str(module))
	p=PlotXY(slicedFrames.refs[2].product["Status"]["RESETINDEX"].data,slicedFrames.refs[2].product["Signal"].data[2,12,:],line=0,symbol = 5, symbolSize = 2)
	p[2]=LayerXY(slicedFrames2.refs[2].product["Status"]["RESETINDEX"].data,slicedFrames2.refs[2].product["Signal"].data[2,12,:],line=0,symbol = 5, symbolSize = 2)
	p.xaxis.title.text="Time [Reset Counter]"
	p.yaxis.title.text="Signal [milibarts]"
	p.title.text="Raster 2, Spaxel 12, Pixel 2"
	plotname = galname + "_" + camera + "_" + obsid  + "_transex.eps"
	p.saveAsEPS(plotname)     #saves in the home directory
	plotname = galname + "_" + camera + "_" + obsid  + "_transex.png"
	p.saveAsPNG(plotname)     #saves in the home directory
	#end Kevin section	
	print "Spectral Flat Fielding"
	slicedFrames2 = specFlatFieldRange(slicedFrames2,polyOrder=5, minWaveRangeForPoly=4., verbose=1)	#old version tested
	slicedCubes = specFrames2PacsCube(slicedFrames2)
	del(slicedFrames2,slicedFrames,qq,qc)
	if verbose:slicedSummary(slicedCubes)
	#slicedCubes_NOFF = getSlicedCopy(slicedCubes)								#New Version to test
	## 1. Flag outliers and rebin
	#waveGrid=wavelengthGrid(slicedCubes, oversample=2, upsample=1, calTree=calTree)
	#slicedCubes = activateMasks(slicedCubes, String1d(["GLITCH","UNCLEANCHOP","NOISYPIXELS","RAWSATURATION","SATURATION","GRATMOVE", "BADPIXELS"]), exclusive = True)
	#slicedCubes = specFlagOutliers(slicedCubes, waveGrid, nSigma=5, nIter=2)
	#slicedCubes = activateMasks(slicedCubes, String1d(["GLITCH","UNCLEANCHOP","NOISYPIXELS","RAWSATURATION","SATURATION","GRATMOVE", "OUTLIERS", "BADPIXELS"]), exclusive = True)
	#slicedRebinnedCubesR = specWaveRebin(slicedCubes, waveGrid)
	## 2. Mask the spectral lines
	#widthDetect =  2.5
	#threshold   = 10.
	#print "find lines..."
	#slicedCubesMask = slicedMaskLines(slicedCubes,slicedRebinnedCubes, widthDetect=widthDetect, widthMask=2.5, threshold=threshold, copy=1, verbose=verbose, maskType="INLINE")
	##print "set lines..."
	##slicedCubesMask = slicedMaskLines(slicedCubes,slicedRebinnedCubes, lineList=[157.7409,121.89757,205.17823,88.356,63.183705], widthDetect=widthDetect, widthMask=2.5, copy=1, verbose=verbose, maskType="INLINE")
	## 3. Actual spectral flatfielding
	#slopeInContinuum = 1
	#slicedCubes = slicedSpecFlatFieldLine(slicedCubesMask, scaling=1, copy=1, maxrange=[50.,230.], slopeInContinuum=slopeInContinuum)
	## 4. Rename mask OUTLIERS to OUTLIERS_1 (specFlagOutliers would refuse to overwrite OUTLIERS) & deactivate mask INLINE
	#slicedCubes.renameMask("OUTLIERS", "OUTLIERS_1")
	#slicedCubes = deactivateMasks(slicedCubes, String1d(["INLINE", "OUTLIERS_1"]))
	#if verbose: maskSummary(slicedCubes, slice=0)
	## 5. Remove intermediate results
	#del waveGrid, slicedRebinnedCubes, slicedCubesMask
	## --- End of Spectral Flat Fielding
	#
	# ------------------------------------------------------------------------------
	#         Processing      Level 1 -> Level 2
	# ------------------------------------------------------------------------------
	calTree=getCalTree()
	linelist = Double1d()
	linelist.append(slicedCubes.refs[1].product["BlockTable"]['LineId'].data[0])
	for i in range (2,slicedCubes.getRefs().size()):
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
	z=slicedCubes.meta["redshiftValue"].double/300000
	for lineloop in range (0,linelist.dimensions[0]):
		if (linewave[lineloop] > 190): 
			linename = "NII205"
		if (linewave[lineloop] > 140)&(linewave[lineloop] < 190): 
			linename = "CII"
		if (linewave[lineloop] > 100)&(linewave[lineloop] < 140): 
			linename = "NII122"
		if (linewave[lineloop] > 75)&(linewave[lineloop] < 100): 
			linename = "OIII"
		if (linewave[lineloop] > 60)&(linewave[lineloop] < 75): 
			linename = "OI"
		lineId      = [int(linelist[lineloop])]
		wavelength  = []
		rasterLine  = []
		rasterCol   = []
		nodCycle    = ""
		scical      = ""
		sliceNumber = []
		onOff="ON"
		sCubesOn = selectSlices(slicedCubes, lineId=lineId, wavelength=wavelength, rasterLine=rasterLine,rasterCol=rasterCol, onOff=onOff, nodCycle=nodCycle, scical=scical, sliceNumber=sliceNumber, verbose=verbose)
		onOff="OFF"
		sCubesOff = selectSlices(slicedCubes, lineId=lineId, wavelength=wavelength, rasterLine=rasterLine,rasterCol=rasterCol, onOff=onOff, nodCycle=nodCycle, scical=scical, sliceNumber=sliceNumber, verbose=verbose)
		if verbose: 
			slicedSummary(sCubesOn)
			slicedSummary(sCubesOff)
			x,y = 2,2
			p6 = plotCubes(sCubesOn,[],x=x,y=y,masks=sCubesOn.refs[0].product.maskTypes, symbol=5, symbolSize=2)
			p6 = plotCubes(sCubesOff,p6,x=x,y=y,masks=sCubesOff.refs[0].product.maskTypes, symbol=5, symbolSize=2)
			plotMasks = String1d(["GLITCH"])
			p6b = plotCubes(sCubesOn,[],x=x,y=y,masks=sCubesOn.refs[0].product.maskTypes)
			p6b = plotCubesMasks(sCubesOn,p6b,x=x,y=y,masks=plotMasks)
		wave = Double1d()
		flux = Double1d()
		reset = Double1d()
		for i in range (0,sCubesOn.refs[1].product["wave"].data.dimensions[0]/16):
			wave.append(sCubesOn.refs[1].product["wave"].data[0+i*16,2,2])
			reset.append(sCubesOn.refs[1].product["Status"]["RESETINDEX"].data[1+i*16])
			flux.append(sCubesOn.refs[1].product["flux"].data[0+i*16,2,2])
		med=getMedian(flux)
		p10=PlotXY(wave,flux,line=0,symbol = 5, symbolSize = 1,yrange=[med-30,med+75])
		for j in range (0,16):
			wave = Double1d()
			flux = Double1d()
			reset = Double1d()
			for i in range (0,sCubesOn.refs[0].product["wave"].data.dimensions[0]/16):
				wave.append(sCubesOn.refs[0].product["wave"].data[j+i*16,2,2])
				reset.append(sCubesOn.refs[0].product["Status"]["RESETINDEX"].data[j+i*16])
				flux.append(sCubesOn.refs[0].product["flux"].data[j+i*16,2,2])
			p10[j+2]=LayerXY(wave,flux,line=0,symbol = 5, symbolSize = 1)
		p10.xaxis.title.text="Wavelength [microns]"
		p10.yaxis.title.text="Signal [milibarts]"
		p10.title.text="Raster 1, Spaxel 2-2, All Pixels"
		plotname = galname + "_" + linename + "_" + obsid  + "_centspax.eps"
		p10.saveAsEPS(plotname)     #saves in the home directory
		plotname = galname + "_" + linename + "_" + obsid  + "_centspax.png"
		p10.saveAsPNG(plotname)     #saves in the home directory
		waveGrid=wavelengthGrid(sCubesOn.refs[1].product, oversample=2, upsample=1, calTree = calTree)
		sCubesOn  = activateMasks(sCubesOn, String1d(["GLITCH","UNCLEANCHOP","RAWSATURATION","SATURATION","GRATMOVE", "BADPIXELS"]), exclusive = True)
		sCubesOn  = specFlagOutliers(sCubesOn, waveGrid, nSigma=3, nIter=4,saveStats=1)
		sCubesOff = activateMasks(sCubesOff, String1d(["GLITCH","UNCLEANCHOP","RAWSATURATION","SATURATION","GRATMOVE", "BADPIXELS"]), exclusive = True)
		sCubesOff = specFlagOutliers(sCubesOff, waveGrid, nSigma=2, nIter=4,saveStats=1)
		if verbose:
			x,y = 2,2
			p7 = plotCubes(sCubesOn,[],x=x,y=y,masks=sCubesOn.refs[0].product.maskTypes, symbol = 14, symbolSize = 1)
			plotMasks = String1d(["OUTLIERS"])
			p7 = plotCubesMasks(sCubesOn,p7,x=x,y=y,masks=plotMasks)
		p8 = plotCubes(sCubesOn,[],x=2,y=2,masks=sCubesOn.refs[0].product.maskTypes, symbol = 14, symbolSize = 1)
		p8.xaxis.title.text="Wavelength [microns]"
		p8.yaxis.title.text="Signal [milibarts]"
		p8.title.text="All Rasters, Spaxel 2-2, Outliers Clipped"
		plotname = galname + "_" + linename + "_" + obsid  + "_rebinFLG.eps"
		p8.saveAsEPS(plotname)     #saves in the home directory
		plotname = galname + "_" + linename + "_" + obsid  + "_rebinFLG.png"
		p8.saveAsPNG(plotname)     #saves in the home directory
		sCubesOn  = activateMasks(sCubesOn, String1d(["GLITCH","UNCLEANCHOP","RAWSATURATION","SATURATION","GRATMOVE", "BADPIXELS", "OUTLIERS"]), exclusive = True)
		slicedRebinnedCubesOn = specWaveRebin(sCubesOn, waveGrid)
		sCubesOff = activateMasks(sCubesOff, String1d(["GLITCH","UNCLEANCHOP","RAWSATURATION","SATURATION","GRATMOVE", "BADPIXELS", "OUTLIERS"]), exclusive = True)
		slicedRebinnedCubesOff = specWaveRebin(sCubesOff, waveGrid)
		if verbose:
			slicedSummary(slicedRebinnedCubesOn)
			slicedSummary(slicedRebinnedCubesOff)
			p8 = plotCubesRaDec(slicedRebinnedCubesOn)
			p8 = plotCubesRaDec(slicedRebinnedCubesOff,p8)
		slicedRebinnedCubesOn  = specAverageCubes(slicedRebinnedCubesOn)
		slicedRebinnedCubesOff = specAverageCubes(slicedRebinnedCubesOff)
		p9 = plotCubes(slicedRebinnedCubesOn,[],x=2,y=2,title="All Rasters, central position",symbol = 14, symbolSize = 3)
		p9.xaxis.title.text="Wavelength [microns]"
		p9.yaxis.title.text="Signal [milibarts]"
		for i in range (0,slicedRebinnedCubesOn.getRefs().size()):
			p9[i].line=2
		exp = Double1d(RESHAPE(slicedRebinnedCubesOn.refs[0].product["exposure"].data[:,2,2]))
		tmp = MAX(NAN_FILTER(exp))
		exp /= tmp
		min_exp = MIN(NAN_FILTER(exp))
		max_exp = MAX(NAN_FILTER(exp))
		margin = 0.05 * (max_exp - min_exp)
		exp_start = min_exp - margin
		exp_end = max_exp + margin
		p9[100] = LayerXY(slicedRebinnedCubesOn.refs[0].product["waveGrid"].data,slicedRebinnedCubesOn.refs[0].product["exposure"].data[:,2,2]/tmp,color=java.awt.Color.black,yrange = [exp_start,exp_end])
		p9[100].yaxis.title.text="EXPOSURE DEPTH"
		plotname = galname + "_" + linename + "_" + obsid  + "_allRastWave.eps"
		p9.saveAsEPS(plotname)     #saves in the home directory
		plotname = galname + "_" + linename + "_" + obsid  + "_allRastWave.png"
		p9.saveAsPNG(plotname)     #saves in the home directory
		if verbose:
			x,y = 2,2
			p6 = plotCubes(sCubesOn,[],x=x,y=y,masks=sCubesOn.refs[0].product.maskTypes)
			p6 = plotCubes(sCubesOff,p6,x=x,y=y,masks=sCubesOff.refs[0].product.maskTypes)
			p6.addLayer(LayerXY(slicedRebinnedCubesOn.refs[0].product["waveGrid"].data,slicedRebinnedCubesOn.refs[0].product["image"].data[:,x,y]))
			p6.addLayer(LayerXY(slicedRebinnedCubesOff.refs[0].product["waveGrid"].data,slicedRebinnedCubesOff.refs[0].product["image"].data[:,x,y]))
		slicedRebinnedCubesOnERR = getSlicedCopy(slicedRebinnedCubesOn)
		for slice in range(len(slicedRebinnedCubesOnERR.refs)):
			slicedRebinnedCubesOnERR.refs[slice].product["image"].data = slicedRebinnedCubesOnERR.refs[slice].product["stddev"].data / SQRT(slicedRebinnedCubesOnERR.refs[slice].product["exposure"].data)
		print 'Subtract the off'
		slicedDiffCubes = getSlicedCopy(slicedRebinnedCubesOn)
		if (len(slicedRebinnedCubesOff.refs) > 2):
			image = slicedRebinnedCubesOff.refs[1].product.image
			for i in range(2,len(slicedRebinnedCubesOff.refs)):
				image2 = slicedRebinnedCubesOff.refs[i].product.image
				image = image + image2
			image = image/(len(slicedRebinnedCubesOff.refs)-1)
			slicedRebinnedCubesOff.refs[0].product.image = image
		if (len(slicedRebinnedCubesOff.refs) == 2):slicedRebinnedCubesOff.refs[0].product.image = slicedRebinnedCubesOff.refs[1].product.image
		for slice in range(len(slicedRebinnedCubesOn.refs)):
			slicedDiffCubes.refs[slice].product.setFlux(slicedRebinnedCubesOn.refs[slice].product.image - slicedRebinnedCubesOff.refs[0].product.image)
			# ADD IN A TEST TO DETERMINE THE CLOSEST OFF...  Rather than averaging them???
		Spectrum = extractSpaxelSpectrum(slicedDiffCubes, slice=slice, spaxelX=2, spaxelY=2)
		sumflux = Spectrum.flux*0
		plot = PlotXY()
		for slice in range(len(slicedRebinnedCubesOn.refs)):
			for i in range (0,5):
				for j in range (0,5):
					Spectrum = extractSpaxelSpectrum(slicedDiffCubes, slice=slice, spaxelX=i, spaxelY=j)
					layer = LayerXY(Spectrum.wave,Spectrum.flux)
					plot.addLayer(layer)
					sumflux += Spectrum.flux
		plot.xaxis.title.text="Wavelength [microns]"
		plot.yaxis.title.text="Signal [milibarts]"
		plot.title.text="All Pixels"
		plotav = PlotXY(Spectrum.wave,sumflux/25/len(slicedRebinnedCubesOn.refs))
		plotav.title.text="Averagre Spectrum"
		plotav.xaxis.title.text="Wavelength [microns]"
		plotav.yaxis.title.text="Signal [milibarts]"
		plotname = galname + "_" + linename + "_" + obsid  + "_allSpec.eps"
		plot.saveAsEPS(plotname)     #saves in the home directory
		plotname = galname + "_" + linename + "_" + obsid  + "_allSpec.png"
		plot.saveAsPNG(plotname)     #saves in the home directory
		plotname = galname + "_" + linename + "_" + obsid  + "_avSpec.eps"
		plotav.saveAsEPS(plotname)     #saves in the home directory
		plotname = galname + "_" + linename + "_" + obsid  + "_avSpec.png"
		plotav.saveAsPNG(plotname)     #saves in the home directory
		name=galname + "_"+linename + "_" + obsid +"_suboff_pipe_PhaseB_veltun15"
		saveSlicedCopy(slicedDiffCubes,name)
		name=galname + "_"+linename + "_" + obsid +"_err_pipe_PhaseB_veltun15"
		saveSlicedCopy(slicedRebinnedCubesOnERR,name)
	del(waveGrid,wavelength,tmp,sumflux,sum,Spectrum,slicedRebinnedCubesOn,exp,exp_end,exp_start,gpr,i,j,layer,lineloop,max_exp,med,min_exp,nodCycle,p10,p8,p9,plot,plotav,qq,rasterCol,rasterLine,sCubesOff,sCubesOn,slice,slicedDiffCubes,slicedRebinnedCubesOff,scical,slicedRebinnedCubesOnERR,wave,sliceNumber,reset,plotname,onOff,flux,count)


print "you have arrived"

print "CONGRATULATIONS! Phase B complete!"
#projectedCube_1 = specProject(slicedDiffCubes.refs[0].product, outputPixelsize=2.85)
#projectedCube_2 = specProject(slicedDiffCubes.refs[1].product, outputPixelsize=2.85)
#projectedCube_3 = specProject(slicedDiffCubes.refs[2].product, outputPixelsize=2.85)
#projectedCube_4 = specProject(slicedDiffCubes.refs[3].product, outputPixelsize=2.85)
#projectedCube_diff = specProject(slicedDiffCubes, outputPixelsize=2.85)

################################################
#                                              #
#    PHASE B: KINGFISH Spectroscopic Pipeline  #
#      BETA version tested in HIPE 8.0.3215    #
#         person to blame: Kevin Croxall       #
#            (aside from the NHSC)             #
#                Jan 19, 2012                  #
################################################
