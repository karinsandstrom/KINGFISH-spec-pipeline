################################################
#                                              #
#    PHASE B: KINGFISH Spectroscopic Pipeline  #
#      BETA version tested in HIPE 8.0.2050    #
#         person to blame: Kevin Croxall       #
#            (aside from the NHSC)             #
#                Oct 26, 2011                  #
################################################
def getra(cube,raster):
	ra = cube.refs[raster].product["ra"].data[:,2,2]
	medra=getMedian(ra)
	return (medra)

def getdec(cube,raster):
	dec = cube.refs[raster].product["dec"].data[:,2,2]
	meddec=getMedian(dec)
	return (meddec)

def getpositions(cube):
	raarray = Double1d()
	decarray = Double1d()
	rastarray = Double1d()
	for rasti in range (0,cube.getRefs().size()):
		deci = getdec(cube,rasti)
		rai = getra(cube,rasti)
		raarray.append(rai)
		decarray.append(deci)
		rastarray.append(rasti)
	posarray = [raarray,decarray,rastarray]
	return (posarray)

def getMedian(numericValues):
  theValues = SORT(numericValues)
  if len(theValues) % 2 == 1:
    return theValues[(len(theValues)+1)/2-1]
  else:
    lower = theValues[len(theValues)/2-1]
    upper = theValues[len(theValues)/2]
    return (float(lower + upper)) / 2 

# A list of the Phase A pools that are to be processed
# need to be saved and referenced here at the begining of the Phase A pipeline.  

phasealist = simpleAsciiTableReader(file = "/Users/kcroxall/dtrans.lst") #UPDATE to the correct file location
ndim = phasealist[0].data.dimensions[0]
verbose = 0
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
	obsid = str(slicedCubes.meta["obsid"].long)
	calTree=getCalTree()
	#then do a for loop over the length of the array for the full processing...
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
			chanwidth = 45
			linewidth = 0.18
			linecent = 205.17823*(1+z)
		if (linewave[lineloop] > 140)&(linewave[lineloop] < 190): 
			linename = "CII"
			linewidth = 0.18
			chanwidth = 45
			linecent = 157.7409*(1+z)
		if (linewave[lineloop] > 100)&(linewave[lineloop] < 140): 
			linename = "NII122"
			chanwidth = 45
			linewidth = 0.18
			linecent = 121.89757*(1+z)
		if (linewave[lineloop] > 75)&(linewave[lineloop] < 100): 
			linename = "OIII"
			chanwidth = 45
			linewidth = 0.065
			linecent = 88.356*(1+z)
		if (linewave[lineloop] > 60)&(linewave[lineloop] < 75): 
			linename = "OI"
			chanwidth = 75
			linewidth = 0.065
			linecent = 63.183705*(1+z)
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
		sCubesOn = selectSlices(slicedCubes, lineId=lineId, wavelength=wavelength, rasterLine=rasterLine,rasterCol=rasterCol, onOff=onOff, nodCycle=nodCycle, scical=scical, sliceNumber=sliceNumber, verbose=verbose)
		onOff="OFF"
		sCubesOff = selectSlices(slicedCubes, lineId=lineId, wavelength=wavelength, rasterLine=rasterLine,rasterCol=rasterCol, onOff=onOff, nodCycle=nodCycle, scical=scical, sliceNumber=sliceNumber, verbose=verbose)
#pre-Ktrans plot
		pos=getpositions(sCubesOn)
		table = TableDataset()
		ascii = AsciiTableTool()
		table["Raster"] = Column(pos[2])
		table["RA"] = Column(pos[0])
		table["DEC"] = Column(pos[1])
		ascii.formatter=FixedWidthFormatter(sizes=[10,25,25])
		tabname = "/Users/kcroxall/" + galname + "_" + linename + "_pos.tab"   # UPDATE PATH
		ascii.save(tabname, table)
		wave = Double1d()
		flux = Double1d()
		reset = Double1d()
		gpr = Double1d()
		x=2
		y=2
		for i in range (0,sCubesOn.refs[0].product["wave"].data.dimensions[0]/16):
			wave.append(sCubesOn.refs[0].product["wave"].data[1+i*16,x,y])
			reset.append(sCubesOn.refs[0].product["Status"]["RESETINDEX"].data[1+i*16])
			flux.append(sCubesOn.refs[0].product["flux"].data[1+i*16,x,y])
			gpr.append(sCubesOn.refs[0].product["Status"]["GPR"].data[1+i*16])
		ind = flux.where((wave<linecent-linewidth).or(wave>linecent+linewidth))
		med=getMedian(flux)
		p1=PlotXY(reset,flux,line=0,symbol = 5, symbolSize = 1,yrange=[med-30,med+75])
		p1.xaxis.title.text="Time (Readout Counter)"
		p1.yaxis.title.text="Signal [milibarts]"
		p1.title.text="Raster 1, Spaxel 2-2, Pixel 1"
		p1[100] = LayerXY(reset[ind], flux[ind],line=0, symbol = 14, symbolSize = 1)
		y = Double1d(RESHAPE(gpr))
		min_y = MIN(NAN_FILTER(y))
		max_y = MAX(NAN_FILTER(y))
		margin = 0.05 * (max_y - min_y)
		y_start = min_y - margin - 75000
		y_end = max_y + margin
		p1[50] = LayerXY(reset, gpr, line=2, stroke=0.2, symbol=14, symbolSize=4,color=java.awt.Color.green,yrange = [y_start,y_end])
		p1[50].yaxis.title.text="GPR"
		plotname = galname + "_" + linename + "_" + obsid   + "_" + "preKtrans.eps"
		p1.saveAsEPS(plotname)     #saves in the kcroxall directory
		plotname = galname + "_" + linename + "_" + obsid   + "_" + "preKtrans.png"
		p1.saveAsPNG(plotname)     #saves in the kcroxall directory
###################### BETA VERSION of using the BASELINE AS AN OFF ######################################
###################### Could be supplemented by the transient correction?
###################### Needs a better model and to treats sharp jumps...
		Dotcloudmod = getSlicedCopy(sCubesOn)
		for rasti in range (0,Dotcloudmod.getRefs().size()):
			print "Raster", rasti
			for q in range (0,5):
				print "X = ",q
				for qq in range (0,5):
#					print "Y = ", qq
					for j in range (0,16):
						wave = Double1d()
						flux = Double1d()
						fluxnew = Double1d()
						reset = Double1d()
						for i in range (0,Dotcloudmod.refs[rasti].product["wave"].data.dimensions[0]/16):
							wave.append(Dotcloudmod.refs[rasti].product["wave"].data[j+i*16,qq,q])
							reset.append(Dotcloudmod.refs[rasti].product["Status"]["RESETINDEX"].data[j+i*16])
							flux.append(Dotcloudmod.refs[rasti].product["flux"].data[j+i*16,qq,q])
							fluxnew.append(Dotcloudmod.refs[rasti].product["flux"].data[j+i*16,qq,q])
						for pix in range (0,reset.dimensions[0]):
							dist = ABS(reset[:] - reset[pix])
							nearby = dist.where((dist<chanwidth).and((wave<linecent-linewidth).or(wave>linecent+linewidth)))
							subval = getMedian(flux[nearby])
							flux2 = flux[nearby]
							dist2 = ABS(flux2 - subval)
							meandist2 = SUM(dist2)/dist2.dimensions[0]
							nearby2 = dist2.where(dist2<meandist2)
							subval2b = SUM(flux2[nearby2])/flux2[nearby2].dimensions[0]
							fluxnew[pix] = flux[pix]-subval2b
						for i in range (0,Dotcloudmod.refs[rasti].product["wave"].data.dimensions[0]/16):
							Dotcloudmod.refs[rasti].product["flux"].data[j+i*16,qq,q] = fluxnew[i]	
			print "cleaning garbage"
			System.gc()
		sCubesOn = getSlicedCopy(Dotcloudmod)
		Dotcloudmod = getSlicedCopy(sCubesOff)
		for rasti in range (0,Dotcloudmod.getRefs().size()):
			print "Raster", rasti
			for q in range (0,5):
				print "X = ",q
				for qq in range (0,5):
#					print "Y = ", qq
					for j in range (0,16):
						wave = Double1d()
						flux = Double1d()
						fluxnew = Double1d()
						reset = Double1d()
						for i in range (0,Dotcloudmod.refs[rasti].product["wave"].data.dimensions[0]/16):
							wave.append(Dotcloudmod.refs[rasti].product["wave"].data[j+i*16,qq,q])
							reset.append(Dotcloudmod.refs[rasti].product["Status"]["RESETINDEX"].data[j+i*16])
							flux.append(Dotcloudmod.refs[rasti].product["flux"].data[j+i*16,qq,q])
							fluxnew.append(Dotcloudmod.refs[rasti].product["flux"].data[j+i*16,qq,q])
						for pix in range (0,reset.dimensions[0]):
							dist = ABS(reset[:] - reset[pix])
							nearby = dist.where((dist<chanwidth).and((wave<linecent-linewidth).or(wave>linecent+linewidth)))
							subval = getMedian(flux[nearby])
							flux2 = flux[nearby]
							dist2 = ABS(flux2 - subval)
							meandist2 = SUM(dist2)/dist2.dimensions[0]
							nearby2 = dist2.where(dist2<meandist2)
							subval2b = SUM(flux2[nearby2])/flux2[nearby2].dimensions[0]
							fluxnew[pix] = flux[pix]-subval2b
						for i in range (0,Dotcloudmod.refs[rasti].product["wave"].data.dimensions[0]/16):
							Dotcloudmod.refs[rasti].product["flux"].data[j+i*16,qq,q] = fluxnew[i]	
			print "cleaning garbage"
			System.gc()
		sCubesOff = getSlicedCopy(Dotcloudmod)
############################################################
############################################################
############################################################
		wave = Double1d()
		flux = Double1d()
		reset = Double1d()
		gpr = Double1d()
		x=2
		y=2
		for i in range (0,sCubesOn.refs[0].product["wave"].data.dimensions[0]/16):
			wave.append(sCubesOn.refs[0].product["wave"].data[1+i*16,x,y])
			reset.append(sCubesOn.refs[0].product["Status"]["RESETINDEX"].data[1+i*16])
			flux.append(sCubesOn.refs[0].product["flux"].data[1+i*16,x,y])
			gpr.append(sCubesOn.refs[0].product["Status"]["GPR"].data[1+i*16])
		ind = flux.where((wave<linecent-linewidth).or(wave>linecent+linewidth))
		med=getMedian(flux)
		p1=PlotXY(reset,flux,line=0,symbol = 5, symbolSize = 1,yrange=[med-30,med+75])
		p1.xaxis.title.text="Time (Readout Counter)"
		p1.yaxis.title.text="Signal [milibarts]"
		p1.title.text="Raster 1, Spaxel 2-2, Pixel 1"
		p1[100] = LayerXY(reset[ind], flux[ind],line=0, symbol = 14, symbolSize = 1)
		y = Double1d(RESHAPE(gpr))
		min_y = MIN(NAN_FILTER(y))
		max_y = MAX(NAN_FILTER(y))
		margin = 0.05 * (max_y - min_y)
		y_start = min_y - margin - 75000
		y_end = max_y + margin
		p1[50] = LayerXY(reset, gpr, line=2, stroke=0.2, symbol=14, symbolSize=4,color=java.awt.Color.green,yrange = [y_start,y_end])
		p1[50].yaxis.title.text="GPR"
		plotname = galname + "_" + linename + "_" + obsid + "_" + "postKtrans.eps"
		p1.saveAsEPS(plotname)     #saves in the kcroxall directory
		plotname = galname + "_" + linename + "_" + obsid + "_" + "postKtrans.png"
		p1.saveAsPNG(plotname)     #saves in the kcroxall directory
		if verbose: 
			slicedSummary(sCubesOn)
			slicedSummary(sCubesOff)
			x,y = 2,2
			p6 = plotCubes(sCubesOn,[],x=x,y=y,masks=sCubesOn.refs[0].product.maskTypes, symbol=5, symbolSize=2)
			p6 = plotCubes(sCubesOff,p6,x=x,y=y,masks=sCubesOff.refs[0].product.maskTypes, symbol=5, symbolSize=2)
			plotMasks = String1d(["GLITCH"])
			p6b = plotCubes(sCubesOn,[],x=x,y=y,masks=sCubesOn.refs[0].product.maskTypes)
			p6b = plotCubesMasks(sCubesOn,p6b,x=x,y=y,masks=plotMasks)
			plotname = galname + "_" + linename + "_" + obsid  + "pretrim.eps"
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
		p10.saveAsEPS(plotname)     #saves in the kcroxall directory
		plotname = galname + "_" + linename + "_" + obsid  + "_centspax.png"
		p10.saveAsPNG(plotname)     #saves in the kcroxall directory
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
#			p7.saveAsEPS(plotname)     #saves in the kcroxall directory
		p8 = plotCubes(sCubesOn,[],x=2,y=2,masks=sCubesOn.refs[0].product.maskTypes, symbol = 14, symbolSize = 1)
		p8.xaxis.title.text="Wavelength [microns]"
		p8.yaxis.title.text="Signal [milibarts]"
		p8.title.text="All Rasters, Spaxel 2-2, Outliers Clipped"
		plotname = galname + "_" + linename + "_" + obsid  + "_rebinFLG.eps"
		p8.saveAsEPS(plotname)     #saves in the kcroxall directory
		plotname = galname + "_" + linename + "_" + obsid  + "_rebinFLG.png"
		p8.saveAsPNG(plotname)     #saves in the kcroxall directory
		sCubesOn  = activateMasks(sCubesOn, String1d(["GLITCH","UNCLEANCHOP","RAWSATURATION","SATURATION","GRATMOVE", "BADPIXELS", "OUTLIERS"]), exclusive = True)
		slicedRebinnedCubesOn = specWaveRebin(sCubesOn, waveGrid)
		sCubesOff = activateMasks(sCubesOff, String1d(["GLITCH","UNCLEANCHOP","RAWSATURATION","SATURATION","GRATMOVE", "BADPIXELS", "OUTLIERS"]), exclusive = True)
		slicedRebinnedCubesOff = specWaveRebin(sCubesOff, waveGrid)
		if verbose:
			slicedSummary(slicedRebinnedCubesOn)
			slicedSummary(slicedRebinnedCubesOff)
			# Sky footprint: the second is overplotted on the first, p8 is first created and then replotted on
			p8 = plotCubesRaDec(slicedRebinnedCubesOn)
			p8 = plotCubesRaDec(slicedRebinnedCubesOff,p8)
			plotname = galname + "_" + linename + "_" + obsid  + "_radec.eps"
#			p8.saveAsEPS(plotname)     #saves in the kcroxall directory
		slicedRebinnedCubesOn  = specAverageCubes(slicedRebinnedCubesOn)
		slicedRebinnedCubesOff = specAverageCubes(slicedRebinnedCubesOff)
		p9 = plotCubes(slicedRebinnedCubesOn,[],x=2,y=2,title="All Rasters, central position",symbol = 14, symbolSize = 3)
		p9.xaxis.title.text="Wavelength [microns]"
		p9.yaxis.title.text="Signal [milibarts]"
		lowv = Double1d(2)
		hiwv = Double1d(2)
		lowv[0] = linecent - linewidth
		lowv[1] = linecent - linewidth
		hiwv[0] = linecent + linewidth
		hiwv[1] = linecent + linewidth
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
		bar = Double1d(2)
		bar[0]=0
		bar[1]=1
		p9[200] = LayerXY(lowv,bar,color=java.awt.Color.black,line=3,yrange = [exp_start,exp_end])
		p9[201] = LayerXY(hiwv,bar,color=java.awt.Color.black,line=3,yrange = [exp_start,exp_end])
		p9[100].yaxis.title.text="EXPOSURE DEPTH"
		plotname = galname + "_" + linename + "_" + obsid  + "_allRastWave.eps"
		p9.saveAsEPS(plotname)     #saves in the kcroxall directory
		plotname = galname + "_" + linename + "_" + obsid  + "_allRastWave.png"
		p9.saveAsPNG(plotname)     #saves in the kcroxall directory
		if verbose:
			x,y = 2,2
			# here we are adding to a previously created plot, but if you deleted plot "p6" from above, then first
			p6 = plotCubes(sCubesOn,[],x=x,y=y,masks=sCubesOn.refs[0].product.maskTypes)
			p6 = plotCubes(sCubesOff,p6,x=x,y=y,masks=sCubesOff.refs[0].product.maskTypes)
			# and then
			p6.addLayer(LayerXY(slicedRebinnedCubesOn.refs[0].product["waveGrid"].data,slicedRebinnedCubesOn.refs[0].product["image"].data[:,x,y]))
			p6.addLayer(LayerXY(slicedRebinnedCubesOff.refs[0].product["waveGrid"].data,slicedRebinnedCubesOff.refs[0].product["image"].data[:,x,y]))
		slicedRebinnedCubesOnERR = getSlicedCopy(slicedRebinnedCubesOn)
		for slice in range(len(slicedRebinnedCubesOnERR.refs)):
			slicedRebinnedCubesOnERR.refs[slice].product["image"].data = slicedRebinnedCubesOnERR.refs[slice].product["stddev"].data / SQRT(slicedRebinnedCubesOnERR.refs[slice].product["exposure"].data)
# Subtract the off.  
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
#deal with more than 1 off cube in the data...
# ADD IN A TEST TO DETERMINE THE CLOSEST OFF...
		Spectrum = extractSpaxelSpectrum(slicedDiffCubes, slice=slice, spaxelX=2, spaxelY=2)
		sumflux = Spectrum.flux*0
		plot = PlotXY()
		for slice in range(len(slicedRebinnedCubesOn.refs)):
			for i in range (0,5):
				for j in range (0,5):
#					spaxelX = i
#					spaxelY = j
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
		plot.saveAsEPS(plotname)     #saves in the kcroxall directory
		plotname = galname + "_" + linename + "_" + obsid  + "_allSpec.png"
		plot.saveAsPNG(plotname)     #saves in the kcroxall directory
		plotname = galname + "_" + linename + "_" + obsid  + "_avSpec.eps"
		plotav.saveAsEPS(plotname)     #saves in the kcroxall directory
		plotname = galname + "_" + linename + "_" + obsid  + "_avSpec.png"
		plotav.saveAsPNG(plotname)     #saves in the kcroxall directory
		name=galname + "_"+linename + "_" + obsid +"_suboff_pipe_PhaseB_ktrans"
		saveSlicedCopy(slicedDiffCubes,name)
		name=galname + "_"+linename + "_" + obsid +"_err_pipe_PhaseB_ktrans"
		saveSlicedCopy(slicedRebinnedCubesOnERR,name)
	del(x,y,y_end,y_start,waveGrid,wavelength,tmp,sumflux,sum,Spectrum,slicedRebinnedCubesOn,bar,exp,exp_end,exp_start,gpr,hiwv,ind,i,j,layer,linecent,lineId,lineloop,linewidth,lowv,margin,max_exp,max_y,med,min_exp,min_y,nodCycle,p1,p10,p8,p9,plot,plotav,qq,rasterCol,rasterLine)
	del(sCubesOff,sCubesOn,slice,slicedDiffCubes,slicedRebinnedCubesOff,scical,slicedRebinnedCubesOnERR,wave,sliceNumber,reset,plotname,onOff,flux,count)

#projectedCube_1 = specProject(slicedDiffCubes.refs[0].product, outputPixelsize=2.85)
#projectedCube_2 = specProject(slicedDiffCubes.refs[1].product, outputPixelsize=2.85)
#projectedCube_3 = specProject(slicedDiffCubes.refs[2].product, outputPixelsize=2.85)
#projectedCube_4 = specProject(slicedDiffCubes.refs[3].product, outputPixelsize=2.85)
#projectedCube_diff = specProject(slicedDiffCubes, outputPixelsize=2.85)

# End Phase B

print "you have arrived"

print "CONGRATULATIONS! Phase B complete!"

################################################
#                                              #
#    PHASE B: KINGFISH Spectroscopic Pipeline  #
#      BETA version tested in HIPE 8.0.2050    #
#         person to blame: Kevin Croxall       #
#            (aside from the NHSC)             #
#                Oct 26, 2011                  #
################################################
