################################################
#                                              #
#    PHASE B: KINGFISH Spectroscopic Pipeline  #
#      BETA version tested in HIPE 7.0.1931    #
#         person to blame: Kevin Croxall       #
#            (aside from the NHSC)             #
#                 Jun 6, 2011                  #
################################################
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
	obsid = str(slicedCubes.meta["obsid"].long)
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

	z=slicedCubes.meta["redshiftValue"].double/300000
	
	for lineloop in range (0,linelist.dimensions[0]):
		if (linewave[lineloop] > 190): 
			linename = "NII205"
			linewidth = 0.14
			linecent = 205.17823*(1+z)
		if (linewave[lineloop] > 140)&(linewave[lineloop] < 190): 
			linename = "CII"
			linewidth = 0.13
			linecent = 157.7409*(1+z)
		if (linewave[lineloop] > 100)&(linewave[lineloop] < 140): 
			linename = "NII122"
			linewidth = 0.12
			linecent = 121.89757*(1+z)
		if (linewave[lineloop] > 75)&(linewave[lineloop] < 100): 
			linename = "OIII"
			linewidth = 0.04
			linecent = 88.356*(1+z)
		if (linewave[lineloop] > 60)&(linewave[lineloop] < 75): 
			linename = "OI"
			linewidth = 0.02
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
		sCubesOn = selectSlices(slicedCubes, lineId=lineId, wavelength=wavelength, rasterLine=rasterLine,\
         rasterCol=rasterCol, onOff=onOff, nodCycle=nodCycle, scical=scical, sliceNumber=sliceNumber, verbose=verbose)
		onOff="OFF"
		sCubesOff = selectSlices(slicedCubes, lineId=lineId, wavelength=wavelength, rasterLine=rasterLine,\
         rasterCol=rasterCol, onOff=onOff, nodCycle=nodCycle, scical=scical, sliceNumber=sliceNumber, verbose=verbose)
#pre-Ktrans plot
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
					print "Y = ", qq
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
							nearby = dist.where((dist<45).and((wave<linecent-linewidth).or(wave>linecent+linewidth)))
							subval = getMedian(flux[nearby])
							flux2 = flux[nearby]
							dist2 = ABS(flux2 - subval)
							meandist2 = SUM(dist2)/dist2.dimensions[0]
							nearby2 = dist2.where(dist2<meandist2)
							subval2 = getMedian(flux2[nearby2])
							subval2b = SUM(flux2[nearby2]))/flux2[nearby2].dimensions[0]
							fluxnew[pix] = flux[pix]-subval2
						for i in range (0,Dotcloudmod.refs[rasti].product["wave"].data.dimensions[0]/16):
							Dotcloudmod.refs[rasti].product["flux"].data[j+i*16,qq,q] = fluxnew[i]	
					System.gc()

		sCubesOn = getSlicedCopy(Dotcloudmod)
		Dotcloudmod = getSlicedCopy(sCubesOff)
		for rasti in range (0,Dotcloudmod.getRefs().size()):
			print "Raster", rasti
			for q in range (0,5):
				print "X = ",q
				for qq in range (0,5):
					print "Y = ", qq
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
							nearby = dist.where((dist<45).and((wave<linecent-linewidth).or(wave>linecent+linewidth)))
							subval = getMedian(flux[nearby])
							flux2 = flux[nearby]
							dist2 = ABS(flux2 - subval)
							meandist2 = SUM(dist2)/dist2.dimensions[0]
							nearby2 = dist2.where(dist2<meandist2)
							subval2 = getMedian(flux2[nearby2])
							subval2b = SUM(flux2[nearby2]))/flux2[nearby2].dimensions[0]
							fluxnew[pix] = flux[pix]-subval2
						for i in range (0,Dotcloudmod.refs[rasti].product["wave"].data.dimensions[0]/16):
							Dotcloudmod.refs[rasti].product["flux"].data[j+i*16,qq,q] = fluxnew[i]	
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
		plotname = galname + "_" + linename + "_" + obsid  + "_" + "postKtrans.eps"
		p1.saveAsEPS(plotname)     #saves in the kcroxall directory
		plotname = galname + "_" + linename + "_" + obsid  + "_" + "postKtrans.png"
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
			p6.saveAsEPS(plotname)     #saves in the kcroxall directory

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

# Note that for the final cube rebinning it is recommended that: 
#    - For single scan SED or Nyquist sampled range scans, it is recommended to perform the final rebinning with 
#        oversample=2, upsample>=1, which corresponds to the native resolution of the instrument
#    - For line scan, deep range scans or Nyquist sampled range scans with repetition factor > 1, 
#        oversampling > 2 is made possible by the high degree of redundancy provided by the observation
#
# The oversample factor is used to increase the number of wavelength 
# bins by the formula bins*oversample, where the number of bins is based 
# on the theoretical resolution of your observation. The upsample factor 
# specifies how many shifts per wavelength bin to make while rebinning. 
# Standard products are generated with oversample=2 and upsample=3 values. 
#
		waveGrid=wavelengthGrid(sCubesOn.refs[0].product, oversample=2, upsample=1, calTree = calTree)
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

		p8 = plotCubes(sCubesOn,[],x=x,y=y,masks=sCubesOn.refs[0].product.maskTypes, symbol = 14, symbolSize = 1)
		p8.xaxis.title.text="Wavelength [microns]"
		p8.yaxis.title.text="Signal [milibarts]"
		p8.title.text="All Rasters, Spaxel 2-2, Outliers Clipped"
		plotname = galname + "_" + linename + "_" + obsid  + "_rebinFLG.eps"
		p8.saveAsEPS(plotname)     #saves in the kcroxall directory
		plotname = galname + "_" + linename + "_" + obsid  + "_rebinFLG.png"
		p8.saveAsPNG(plotname)     #saves in the kcroxall directory
		plotname = galname + "_" + linename + "_" + obsid  + "_posttrim_flag.eps"


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
		p9 = plotCubes(slicedRebinnedCubesOn,[],x=x,y=y,title="All Rasters, central position",symbol = 14, symbolSize = 3)
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
			plotname = galname + "_" + linename + "withspec.eps"
#			p6.saveAsEPS(plotname)     #saves in the kcroxall directory

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

		Spectrum = extractSpaxelSpectrum(slicedDiffCubes, slice=slice, spaxelX=2, spaxelY=2)
		sumflux = Spectrum.flux*0
		plot = PlotXY()
		for slice in range(len(slicedRebinnedCubesOn.refs)):
			for i in range (0,5):
				for j in range (0,5):
					spaxelX = i
					spaxelY = j
					Spectrum = extractSpaxelSpectrum(slicedDiffCubes, slice=slice, spaxelX=spaxelX, spaxelY=spaxelY)
					layer = LayerXY(Spectrum.wave,Spectrum.flux)
					plot.addLayer(layer)
					sumflux += Spectrum.flux

		plot.xaxis.title.text="Wavelength [microns]"
		plot.yaxis.title.text="Signal [milibarts]"
		plotav = PlotXY(Spectrum.wave,sumflux/25/4)
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
		


		if verbose:
			slice = 0
			p10 = plotCube5x5(slicedModOffCube.refs[1].product)
			p11 = plotCube5x5(slicedRebinnedCubesOff.refs[1].product)
			p9 = plotCube5x5(slicedDiffCubesMod.refs[slice].product)
			p9 = plotCube5x5(slicedDiffCubes.refs[slice].product)
			p10 = plotCube5x5(slicedRebinnedCubesOn.refs[slice].product)
			plotname = galname + "_" + linename + "_" + obsid  + "_5x5_slice0_nooff.eps"
			p10.saveAsEPS(plotname)     #saves in the kcroxall directory
			plotname = galname + "_" + linename + "_" + obsid  + "_5x5_slice0_offsub.eps"
			p9.saveAsEPS(plotname)     #saves in the kcroxall directory
		name=galname + "_"+linename + "_" + obsid +"_suboff_pipe_PhaseB_med"
		saveSlicedCopy(slicedDiffCubes,name)
		name=galname + "_"+linename + "_" + obsid +"_submod_pipe_PhaseB_med"
		saveSlicedCopy(slicedDiffCubesMod,name)
		name=galname + "_"+linename + "_" + obsid +"_nooff_pipe_PhaseB_med"
		saveSlicedCopy(slicedRebinnedCubesOn,name)
		name=galname + "_"+linename + "_" + obsid +"_err_pipe_PhaseB_med"
		saveSlicedCopy(slicedRebinnedCubesOnERR,name)


# WHICH ONE TO SAVE??????
# DO WE WANT ONLY MAPS OF EVERYTHING TOGETHER?
# OR ALSO THE MAPS OF THE INDIVIDUAL AORS???
projectedCube_1 = specProject(slicedDiffCubes.refs[0].product, outputPixelsize=2.85)
projectedCube_2 = specProject(slicedDiffCubes.refs[1].product, outputPixelsize=2.85)
projectedCube_3 = specProject(slicedDiffCubes.refs[2].product, outputPixelsize=2.85)
projectedCube_4 = specProject(slicedDiffCubes.refs[3].product, outputPixelsize=2.85)

# End Phase B

projectedCube_difmod = specProject(slicedDiffCubesMod, outputPixelsize=2.85)
projectedCube_diff = specProject(slicedDiffCubes, outputPixelsize=2.85)
projectedCube_nooff = specProject(slicedRebinnedCubesOn, outputPixelsize=2.85)
projectedCube_nooff = specProject(slicedRebinnedCubesOn, outputPixelsize=1)

print "you have arrived"

print "CONGRATULATIONS! Phase B complete!"

################################################
#                                              #
#    PHASE B: KINGFISH Spectroscopic Pipeline  #
#      BETA version tested in HIPE 7.0.1931    #
#         person to blame: Kevin Croxall       #
#            (aside from the NHSC)             #
#                 Jun 6, 2011                  #
################################################
