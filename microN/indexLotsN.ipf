#pragma rtGlobals=1		// Use modern global access method.
#pragma ModuleName=indexLots
#pragma version = 2.37
#include "ArrayOf3dOrientsN", version>=2.58
#include "DepthResolvedQueryN", version>=1.52
#include "IndexingN", version>=4.45
#include "xmlMultiIndex", version>=1.66

Static StrConstant strainConstraints="110111"


//#include  "microGeometry", version>=2.3
//#include  "LatticeSym", version>=3.1


//  June 25, 2007   version 1.01
// Keep un-indexed image information (so I can still use the intensity information)
//
//	Aug 1, 2007, with version 2, the lines in indexingLIst start with image name and not two indicies. each line is really a list now
//
//	Apr 15, 2009, created to the new geometry version
//
//	Apr 17, 2009, added Tetragonal SymmetryOps
//
//	Aug 31, 2009, added the include of xmlMultiIndex
//
//	Dec 1,  2009, re-did SymmetryOps, to use symmetry from Space Group, should work for everything
//
//	Apr 22, 2014, changed IndexLots to be more modern, also now can write an xml file like the cluster.
//
//	Dec 27, 2016, Add2XmlFile() now returns the full path to the xml file, moved the xmlns=... out of <step> into <AllSteps>.


Menu "Rotations"
	"-"
	SubMenu("Index Lots")
		"Index lots of reconstructed Images",IndexLots("","","","",NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN)
		MenuItemIfWaveClassExists("Write 3d Orientations file","indexationResultList","TEXT:1"),Write3dOrientsFile($"","")
		MenuItemsWaveClassOnGraph("    Find Closest point to cursor","OrientationSliceWave*;RGBfullOrientationMap",""),findClosestPoint(NaN,NaN)
		MenuItemsWaveClassOnGraph("    Put Markers at Reference Points","OrientationSliceWave*;RGBfullOrientationMap",""),putMarkersAtRefPoints("FF")
		SubMenu("Twins")
			"Enable BCC Twins",EnableBCC_Twins()
			"Disable BCC Twins",DisableBCC_Twins()
		End
	End
	SubMenu("Trim Top Surface & Bottom")
		"Trim Both Top & Bottom",TrimTopAndBottomSurfaces(NaN)
		"   Trim Top",TrimTopSurface3d()
		"   Trim Bottom",TrimBotSurface3d(NaN)
	End
End



//if (str2num(StringByKey("IGORVERS",IgorInfo(0)))>=6.0)
//	DoAlert 0, "You are running Igor 6, you cannot Interpolate with version 6!  Should work now"
//	print "You are running Igor 6, you cannot Interpolate with version 6!  Should work now,  see line with enoise(1e-4)"
//endif

// IndexLots("reconPath","WhWi_","1-1080","19-67",0,0,2,30,18,0.5)

//	IndexLots("7597-14982","50-80")
//	IndexLots("1-39036","50-80")					// index everything!

//	IntensityVsDepth("1-7596","0-95")			// X1
//	IntensityVsDepth("7597-15192","0-95")		// X2
//	IntensityVsDepth("15193-22788","0-95")		// X3
//	IntensityVsDepth("22789-30384","0-95")		// X4
//	IntensityVsDepth("30385-37980","0-95")		// X5
//	IntensityVsDepth("37981-45576","0-95")		// X6


//	IndexLots("1-8979","20-75")

//	root:wireA:wireA_3724_55,root:wireA:wireA_3724,root:wireA:wireA_3724_30



//$refMat	0.0469741,-0.0324687,-0.0004100,0.0323781,0.0468902,-0.0037312,0.0024582,0.0028368,0.0569813 // reference matrix identifies zero angle, same order as data
////
////     g11 g12 g13        h        x
//// g = g21 g22 g23,  and  k = g x  y
////     g31 g32 g33        l        z
////
//// where (x,y,z) are in the X,H,F coordinate system, g is recip matrix
////
////	1st 3 columns X,Y,Z locations of voxel in beamline coordinate system (microns)
////	next 9 columns orientation matrix g11 g12 g13 g21 g22 g23 g31 g32 g33
////
//RotationMatrix
//0.0460287 0.00615766 -0.0332325		 -0.0103891 0.0560083 -0.00401165		 0.0321618 0.00927953 0.0462652

//	¥IndexLots("reconPath","Inconel-718_","1-2","13-15", 140, 0,0,2,72,17,0.1,0, minPeakWidth=1.13,maxPeakWidth=18,keVmaxTest=30,FitPeaksFunc=FitPeaksWithExternal,saveXML=1,printIt=1)




//Function IndexLots(pathName,filePrefix,positions,depths, threshAboveAvg, h,k,l,cone,keVmaxCalc,angleTolerance,doStrain,[minPeakWidth,maxPeakWidth,minSep,keVmaxTest,maxNum,FitPeaksFunc])
Function IndexLots(pathName,filePrefix,positions,depths, threshAboveAvg, h,k,l,cone,keVmaxCalc,angleTolerance,doStrain,[minPeakWidth,maxPeakWidth,minSep,keVmaxTest,maxNum,FitPeaksFunc,saveXML,printIt])
	String pathName					// name of path to use
	String filePrefix					// first part of file name, e.g. "WH_"
	String positions, depths		// string ranges (non-existant file names are jsut skipped), if positions=="", then just use one index
	Variable threshAboveAvg		// threshold above average value in FitPeaksStepWise (try 10 or 20 are good numbers)
	Variable h,k,l						// hkl near center (approx)
	Variable cone						// cone angle to the central hkl
	Variable keVmaxCalc				// 14, maximum energy to calculate (keV)
	Variable angleTolerance		// 0.25, angular tolerance (deg)
	Variable doStrain					// flag, true=also compute the strain (=1 use BL, =2 useXHF, =3 use Sample coordinates)

	Variable minPeakWidth			// min fwhm for a peak (in either x or y) =1.2
	Variable maxPeakWidth			// max fwhm for a peak (in either x or y) =35
	Variable minSep					// min separation between two peaks =50
	Variable keVmaxTest				// max energy to check for the final pass in the indexation, =50
	Variable maxNum					// maximum number of peaks to find, normally goes to completion
	FUNCREF FitPeaksProtoE FitPeaksFunc
	Variable saveXML					// flag to save the xml file, default is True
	Variable printIt

	minPeakWidth = ParamIsDefault(minPeakWidth) ? 1.2 : minPeakWidth
	maxPeakWidth = ParamIsDefault(maxPeakWidth) ? 35 : maxPeakWidth
	minSep = ParamIsDefault(minSep) ? 50 : minSep
	keVmaxTest = ParamIsDefault(keVmaxTest) ? 50 : keVmaxTest
	maxNum = ParamIsDefault(maxNum) ? Inf : maxNum
	saveXML = ParamIsDefault(saveXML) ? NaN : saveXML
	saveXML = numtype(saveXML) ? 1 : !(!saveXML)
	printIt = ParamIsDefault(printIt) ? NaN : printIt
	printIt = numtype(printIt) ? 0 : !(!printIt)
	if (ParamIsDefault(FitPeaksFunc))
		FUNCREF FitPeaksProtoE FitPeaksFunc=FitPeaksWithExternal
	endif
	String funcName = FitPeaksFunc($"",NaN,NaN,NaN,NaN,whoami=1)		// get name of FitPeaksFunc
	if (strlen(funcName)<1)
		Abort "Illegal FitPeaksFunc passed to IndexLots"
	endif

	Variable minNpeaks=5					// minimum # of peaks successfully indexed

	String list = depthResolve#GetFileRootAndDoubleRange(pathName,filePrefix,positions,depths)
	if (strlen(list)<1)
		return 1
	endif
	pathName = StringByKey("pathName",list,"=")
	filePrefix = StringByKey("filePrefix",list,"=")
	positions = StringByKey("positions",list,"=")
	depths = StringByKey("depths",list,"=")
	String fileRoot = StringByKey("fileRoot",list,"=")

	if (!(threshAboveAvg>0) || numtype(h+k+l) || !(cone>0) || !(keVmaxCalc>0) || !(angleTolerance>0) || !(doStrain>=0) || !(doStrain<=3))
		threshAboveAvg = !(threshAboveAvg>0) ? 20 : threshAboveAvg
		h = numtype(h) ? 0 : h
		k = numtype(k) ? 0 : k
		l = numtype(l) ? 1 : l
		cone = cone>0 ? cone : 30
		keVmaxCalc = keVmaxCalc>0 ? keVmaxCalc : 17
		angleTolerance = angleTolerance>0 ? angleTolerance : 0.5
		doStrain = (!(doStrain>=0) || !(doStrain<=3)) ? 0 : limit(doStrain,0,3)
		String hklStr
		sprintf hklStr,"%g %g %g",h,k,l
		Prompt threshAboveAvg, " threshold above average value in peak fitting"
		Prompt hklStr,"central (hkl)"
		Prompt cone,"cone angle from the hkl (degree)"
		Prompt keVmaxCalc,"maximum energy for calculating indexing (keV)"
		Prompt angleTolerance,"angular tolerance (degree)"
		Prompt doStrain, "also compute the strain",popup,"NO Strain;Beam LIne coords;* XHF coords;Sample (outward normal)"
		doStrain += 1
		String funcList = FunctionList("Fit*", ";", "KIND:2;NPARAMS:7;VALTYPE:4")
//		funcList = RemoveFromList("FitPeaksProto",funcList)
		Prompt funcName,"peak fitting function",popup,funcList
		DoPrompt "(hkl) for indxing",hklStr,cone,keVmaxCalc,angleTolerance,doStrain,threshAboveAvg,funcName
		if (V_flag)
			return 1
		endif
		doStrain -= 1
		hklStr= ReplaceString(",",hklStr, " ")
		sscanf hklStr, "%g %g %g",h,k,l
		if (V_flag!=3)
			h = NaN ; k=NaN ; l=NaN
		endif
		FUNCREF FitPeaksProtoE FitPeaksFunc=$funcName
		funcName = FitPeaksFunc($"",NaN,NaN,NaN,NaN,whoami=1)	// get name of FitPeaksFunc
		printIt = 1
	endif

	if (NumberByKey("printIt",list,"=") || printIt)
		printf " ¥IndexLots(\"%s\",\"%s\",\"%s\",\"%s\", %g, %g,%g,%g, %g, %g, %g, %d",pathName,filePrefix,positions,depths,threshAboveAvg,h,k,l,cone,keVmaxCalc,angleTolerance,doStrain
		if (!ParamIsDefault(minPeakWidth))
			printf ", minPeakWidth=%g" minPeakWidth
		endif
		if (!ParamIsDefault(maxPeakWidth))
			printf ", maxPeakWidth=%g" maxPeakWidth
		endif
		if (!ParamIsDefault(minSep))
			printf ", minSep=%g" minSep
		endif
		if (!ParamIsDefault(keVmaxTest))
			printf ", keVmaxTest=%g" keVmaxTest
		endif
		if (!ParamIsDefault(maxNum))
			printf ", maxNum=%g" maxNum
		endif
		if (!ParamIsDefault(FitPeaksFunc))
			printf ", FitPeaksFunc=%s" funcName
		endif
		if (!ParamIsDefault(saveXML))
			printf ", saveXML=%g" saveXML
		endif
		if (!ParamIsDefault(printIt))
			printf ", printIt=%g" printIt
		endif
		printf ")\r"
	endif	
	if (!(threshAboveAvg>0) || numtype(h+k+l) || !(cone>0) || !(keVmaxCalc>0) || !(angleTolerance>0) || !(doStrain>=0) || !(doStrain<=3))
		DoAlert 0,"invalid values"
		return 1
	endif
	printf "fitting peaks using peak widths between [%g, %g],  min peak separation of %g,  and a threshold above average of %g\r",minPeakWidth,maxPeakWidth,minSep,threshAboveAvg
	// FitPeaksStepWise(image,1.2,35,50,10)	// FitPeaksStepWise(image,minPeakWidth,maxPeakWidth,minSep,threshAboveAvg)
	printf "equiv \t\t%s(image,%g,%g,%g,%g)\r",funcName,minPeakWidth,maxPeakWidth,minSep,threshAboveAvg
	printf "equiv \t\tIndexAndDisplay(FullPeakList,%g,%g,%g, %g,%g,%g,%g)\r",keVmaxCalc,keVmaxTest,angleTolerance, h,k,l,cone
	Variable epoch = DateTime, seconds
	printf "running in data folder  '%s'\r",GetDataFolder(1)

	String progressWin = ProgressPanelStart("",stop=1,showTime=1,status="Starting")	// display a progress bar
	Make/N=50/T/O IndexingList=""										// holds results of each indexing
	String str, indexNote =  ReplaceStringByKey("waveClass","","indexationResultList","=")
	PathInfo pathName
	if (V_flag)
		indexNote = ReplaceStringByKey("imageFilePath",indexNote,S_path,"=")
	endif
	indexNote = ReplaceStringByKey("filePrefix",indexNote,filePrefix,"=")
	if (strlen(positions))
		indexNote = ReplaceStringByKey("positions",indexNote,positions,"=")
	endif
	indexNote = ReplaceStringByKey("depths",indexNote,depths,"=")
	sprintf str,"{%g,%g,%g}",h,k,l
	indexNote = ReplaceStringByKey("hkl",indexNote,str,"=")
	indexNote = ReplaceNumberByKey("cone",indexNote,cone,"=")
	indexNote = ReplaceNumberByKey("keVmaxCalc",indexNote,keVmaxCalc,"=")
	indexNote = ReplaceNumberByKey("angleTolerance",indexNote,angleTolerance,"=")
	indexNote = ReplaceNumberByKey("minPeakWidth",indexNote,minPeakWidth,"=")
	indexNote = ReplaceNumberByKey("maxPeakWidth",indexNote,maxPeakWidth,"=")
	indexNote = ReplaceNumberByKey("minSep",indexNote,minSep,"=")
	indexNote = ReplaceNumberByKey("threshAboveAvg",indexNote,threshAboveAvg,"=")
	indexNote = ReplaceNumberByKey("minPeakWidth",indexNote,minPeakWidth,"=")
	indexNote = ReplaceNumberByKey("maxPeakWidth",indexNote,maxPeakWidth,"=")
	indexNote = ReplaceNumberByKey("minSep",indexNote,minSep,"=")
	indexNote = ReplaceNumberByKey("keVmaxTest",indexNote,keVmaxTest,"=")
	if (doStrain)
		indexNote = ReplaceStringByKey("strainCoords",indexNote,SelectString(doStrain-2,"BeamLine","XHF","Sample"),"=")
		printf "saving strain in %s coordinates\r",SelectString(doStrain-2,"BeamLine","XHF","Sample")
	endif
	Note/K IndexingList, indexNote

	SVAR imageExtension=root:Packages:imageDisplay:imageExtension
	Variable Npoints=0
	Variable a,b,c,alpha,bet,gam
	Variable depth,NpatternsFound
	Variable X1,Y1,Z1
	Variable lo=Inf,hi=-Inf
	String wnote, fname,imageName, xmlStr="", xmlFileName=""
	Variable Nimages = max(1,ItemsInRange(positions))*max(1,ItemsInRange(depths))
	Variable i=1, j, N=0, indexed=0, once, stopped=0
	for (i=str2num(positions),once=1; (!numtype(i) || once) && !stopped; i=NextInRange(positions,i))
		once = 0
		if (numtype(i))
			sprintf fname, "%s%d%s",fileRoot,str2num(depths),imageExtension
		else
			sprintf fname, "%s%d_%d%s",fileRoot,i,str2num(depths),imageExtension
		endif
		GetFileFolderInfo/P=$pathName/Q/Z fname
		if (V_flag)																	// file not found, skip this position
			continue
		endif

		for (j=str2num(depths); !numtype(j); j=NextInRange(depths,j))
			if (numtype(i))
				sprintf fname, "%s%d%s",fileRoot,j,imageExtension
			else
				sprintf fname, "%s%d_%d%s",fileRoot,i,j,imageExtension
			endif
			imageName = LoadGenericImageFile(fname)
			if (exists(imageName)!=1)
				continue
			endif
			Wave image = $imageName
			imageName = NameOfWave(image)
			seconds = DateTime-epoch
			sprintf str,"%s    %.2f sec/image",imageName,(seconds/N)
			if (ProgressPanelUpdate(progressWin,100*N/Nimages,status=str))// update progress bar
				stopped = 1
				break																	//   and break out of loop
			endif
			N += 1

			wnote = note(image)
	 		depth = NumberByKey("depth", wnote,"=")						// values read from image file
			if (numtype(depth))
	 			depth = NumberByKey("depthSi", wnote,"=")				// old style name
			endif
			X1 = NumberByKey("X1", wnote,"=")								// positioner positions
			Y1 = NumberByKey("Y1", wnote,"=")
			Z1 = NumberByKey("Z1", wnote,"=")
			if (abs(sum(image))<2)												// image is empty (try to process when negative)
				KillWaves/Z image
				continue
			endif

			Wave FullPeakList = $FitPeaksFunc(image,minPeakWidth,maxPeakWidth,minSep,threshAboveAvg,maxNum=maxNum)
			if (!WaveExists(FullPeakList))
				KillWaves/Z image
				continue
			endif
			KillWaves/Z image
			wnote = note(FullPeakList)
			str = ""
			str = ReplaceStringByKey("image",str,imageName,"=")
			if (numtype(depth)==0)
				str = ReplaceNumberByKey("depth",str,depth,"=")
			endif
			str = ReplaceNumberByKey("X1",str,X1,"=")
			str = ReplaceNumberByKey("Y1",str,Y1,"=")
			str = ReplaceNumberByKey("Z1",str,Z1,"=")
			str = ReplaceNumberByKey("totalPeakIntensity",str,NumberByKey("totalPeakIntensity",wnote,"="),"=")
			str = ReplaceNumberByKey("totalIntensity",str,NumberByKey("totalIntensity",wnote,"="),"=")

			if (DimSize(FullPeakList,0)>=minNpeaks)						// enough peaks found, try to index
				Wave FullPeakIndexed=IndexAndDisplay(FullPeakList,keVmaxCalc,keVmaxTest,angleTolerance, h,k,l,cone,printIt=0)
				if (WaveExists(FullPeakIndexed))
					wnote = note(FullPeakIndexed)
					NpatternsFound = NumberByKey("NpatternsFound",wnote,"=")
			 		if (NpatternsFound>0)										// a valid indexation, save the information
						indexed += 1												// counts number of patterns indexed
						str = ReplaceNumberByKey("NpatternsFound",str,NpatternsFound,"=")
						str = ReplaceNumberByKey("Nindexed",str,NumberByKey("Nindexed",wnote,"="),"=")
						str = ReplaceNumberByKey("NiData",str,NumberByKey("NiData",wnote,"="),"=")
						str = ReplaceNumberByKey("goodness",str,NumberByKey("goodness0",wnote,"="),"=")
						str = ReplaceNumberByKey("rms_error",str,NumberByKey("rms_error0",wnote,"="),"=")
						str = ReplaceStringByKey("rlattice",str,StringBykey("recip_lattice0",wnote,"="),"=")
						str = ReplaceStringByKey("EulerAngles",str,StringBykey("EulerAngles0",wnote,"="),"=")
						if (doStrain)
							Wave epsilon = $(DeviatoricStrainRefine(0,strainConstraints,coords=doStrain))
						endif
						if (WaveExists(epsilon))
							sscanf StringByKey("lattice_constants_refined",note(epsilon),"="),"{%g %g %g %g %g %g}",a,b,c,alpha,bet,gam
							if (V_flag==6)
								str = ReplaceNumberByKey("a",str,a,"=")
								str = ReplaceNumberByKey("b",str,b,"=")
								str = ReplaceNumberByKey("c",str,c,"=")
								str = ReplaceNumberByKey("alpha",str,alpha,"=")
								str = ReplaceNumberByKey("bet",str,bet,"=")
								str = ReplaceNumberByKey("gam",str,gam,"=")
							endif
							str = ReplaceNumberByKey("exx",str,epsilon[0][0],"=")
							str = ReplaceNumberByKey("eyy",str,epsilon[1][1],"=")
							str = ReplaceNumberByKey("ezz",str,epsilon[2][2],"=")
							str = ReplaceNumberByKey("exy",str,epsilon[0][1],"=")
							str = ReplaceNumberByKey("exz",str,epsilon[0][2],"=")
							str = ReplaceNumberByKey("eyz",str,epsilon[1][2],"=")
						endif
					endif
				endif
			endif
			if (doStrain && Npoints==0)										// make sure something is there for the first one
				if (!keyInList("exx",str,"=",";"))
					str = ReplaceNumberByKey("exx",str,NaN,"=")		// somethiing must be here for first one, needed in Write3dOrientsFile()
					str = ReplaceNumberByKey("a",str,NaN,"=")
				endif
			endif

			if (saveXML)
				xmlStr = CreateXmlStepStr(FullPeakList,FullPeakIndexed)
				xmlFileName = Add2XmlFile("reconPath",xmlFileName,xmlStr)
			endif

			if (DimSize(IndexingList,0)<=Npoints)							// wave is not long enough
				Save/T/O/P=home indexinglist as "indexinglist.itx"	// in case things die
				Redimension/N=(DimSize(IndexingList,0)+50) IndexingList	// increase size
			endif
			lo = min(j,lo)															// just save max range of depth indicies
			hi = max(j,hi)
			IndexingList[Npoints] = str
			Npoints += 1
		endfor
	endfor
	Redimension/N=(Npoints) IndexingList									// set to final size
	if (stopped)
		printf "\r***   STOPPED by User   ***\r\r"
	endif
	printf "Successfully indexed %d of the %d images processed.\r",indexed,N

	seconds = DateTime-epoch
	printf "maximum depth range from indexinglist is [%d, %d]\r",lo,hi
	print "done, total execution time is  ",Secs2Time(seconds,5,0)
	if (strlen(xmlFileName))
		printf "saved the xml file to  \"%s\"\r",xmlFileName
	else
		printf "NO xml file saved"
	endif
	DoWindow/K $progressWin
	if (seconds>10*60)															// if execution took more than 10 min, automatically save
		print "This took more than 10min, so save the experiment"
		SaveExperiment
	endif
End
//
Function/S FitPeaksProtoE(image,minPeakWidth,boxSize,maxRfactor,threshAboveAvg,[mask,whoami,peakShape,smoothing,thresholdRatio,maxNum,printIt])
	Wave image						// the image from the file
	Variable minPeakWidth		// minimum width for an acceptable peak (pixels)
	Variable boxSize				// maximum width for an acceptable peak (pixels)
	Variable maxRfactor			// max R-factor
	Variable threshAboveAvg	// threshold above average value,  use 25
	Wave mask						// starting mask for the fit (this mask is unchanged by this routine)
	Variable whoami				// if passed, then just return the name of this routine
	String peakShape				// Lorentzian (default), or Gaussian
	Variable smoothing			// if TRUE, do a smoothing operation before peak search
	Variable thresholdRatio	// set threshold to (ratio*[std dev] + avg) if threshold passed as NaN (optional)
	Variable maxNum				// maximum number of peaks to fit (defaults to unlimited)
	Variable printIt

	return ""
End


////		 	NewFitPeaks(image,110,40)
////		 	NewFitPeaks(image,500,40)
////			NewFitPeaks(image,500,60)
////			ImageStats/M=1 image
////			NewFitPeaks(image,V_avg*105+90,40)
////			FitPeaksStepWise(image,2,20,145,90)
////			FitPeaksStepWise(image,2,20,145,130)
////			FitPeaksStepWise(image,2,20,50,40)
////			FitPeaksStepWise(image,2,35,50,40)
////			FitPeaksStepWise(image,1.2,35,50,10)							// FitPeaksStepWise(image,minPeakWidth,maxPeakWidth,minSep,threshAboveAvg)
//			FitPeaksStepWise(image,minPeakWidth,maxPeakWidth,minSep,threshAboveAvg)
//			Wave FullPeakList=FullPeakList
//			if (DimSize(FullPeakList,0)<6)										// not enough peaks found ot index
//				KillWaves image
//				continue
//			endif
//			KillWaves/Z FullPeakIndexed
//
////			IndexAndDisplay(FullPeakList,14,26,0.25, 0,2,2, 50)
////			IndexAndDisplay(FullPeakList,14,26,0.25, 5,1,5, 22)
////			IndexAndDisplay(FullPeakList,14,39,0.25, 1,2,3,72)
////			IndexAndDisplay(FullPeakList,14,35,0.25, 0,2,2,45)
////			IndexAndDisplay(FullPeakList,14,35,0.25, 0,2,2,45)
//			IndexAndDisplay(FullPeakList,keVmaxCalc,50,angleTolerance, h,k,l,cone)	// IndexAndDisplay(FullPeakList,keVmaxCalc,keVmaxTest,angleTolerance,hp,kp,lp,cone)
//			KillWaves image






// Trim off bad points above the surface and too deep
Function TrimTopAndBottomSurfaces(bottomLevel)
	Variable bottomLevel
	TrimTopsurface3d()
	TrimBotSurface3d(bottomLevel)
End
//
Function TrimBotSurface3d(level)				// trim points that are too deep
	Variable level									// trim all points with TotalIntensity3d less than level.  If level in (0,1) then use level of (V_max/level)

	Wave R3dX=R3dX, R3dH=R3dH, R3dF=R3dF, Intens3D=TotalIntensity3d
	if (!WaveExists(R3dX) || !WaveExists(R3dX) || !WaveExists(R3dX) || !WaveExists(Intens3D))
		DoAlert 0, "the waves R3dX, R3dH, R3dF, or TotalIntensity3d do not exists in the current folder"
	endif
	Variable printIt=level<1
	WaveStats/M=1/Q Intens3D
	if (!(level>0))
		level = 1/12
		Prompt level,"level to trim bottom (fraction of max, or intensity if > 1)"
		DoPrompt "level",level
		if (V_flag)
			return 1
		endif
		printIt = stringmatch(GetRTStackInfo(0),"*PopMenuProc;TrimBotSurface3d;") ? 2 : 1
	endif
	level = (0<level && level<1)? V_max*level : level
	if (!(level>0))
		DoAlert 0, "cannot trim bottom to level = "+num2str(level)
	endif
	if (printIt)
		Variable d, f = round(360*level/V_max)		// now answer  = f/100
		d= gcd(f,360)
		String fstr=""
		if (abs(level/V_max - f/360)<1e-5)
			sprintf fstr,"%d/%d",f/d,360/d
		else
			fstr = num2str(level/V_max)
		endif
		printf "%sTrimBotSurface3d(%g)		//  or  TrimBotSurface3d(%s)\r",SelectString(printIt==2,"","¥"),level,fstr
	endif

	Variable NX=DimSize(R3dX,0), NH=DimSize(R3dX,1), NF=DimSize(R3dX,2)
	Make/N=(NF)/O/D top_temp_
	Wave wtop = top_temp_
	SetScale/P x,0,1,"",wtop

	Make/N=(NX,NH)/O/D bottomHeight=NaN
	SetScale/P x,DimOffset(R3dX,0),DimDelta(R3dX,0),WaveUnits(R3dX,0), bottomHeight
	SetScale/P y,DimOffset(R3dX,1),DimDelta(R3dX,1),WaveUnits(R3dX,1), bottomHeight
	SetScale d 0,0,WaveUnits(R3dX,2), bottomHeight

//	Variable i,j, m=ceil((Fmax-DimOffset(R3dX,2)) / DimDelta(R3dX,2)+.01)
	Variable i,j, m
	for (j=0;j<NH;j+=1)
		for (i=0;i<NX;i+=1)
			wtop = Intens3D[i][j][p]
			WaveStats/Q/M=1 wtop
			if (!(V_min<level))				// never gets low enough
				continue
			endif
			for (m=V_maxloc;m<NF;m+=1)		// search forward from peak to find level
				if (!(wtop[m]>level))
					break
				endif
			endfor
			if (m>=NF)							// never got low enough
				continue
			endif
			bottomHeight[i][j] = DimOffset(R3dX,2) + m*DimDelta(R3dX,2)
			R3dX[i][j][m,NF-1] = NaN
			R3dH[i][j][m,NF-1] = NaN
			R3dF[i][j][m,NF-1] = NaN
		endfor
	endfor
	KillWaves/Z top_temp_
	return 0
End
//
Function TrimTopSurface3d()					// find the top surface, and trim all points above it
	Wave R3dX=R3dX, R3dH=R3dH, R3dF=R3dF, Intens3D=TotalIntensity3d
	if (!WaveExists(R3dX) || !WaveExists(R3dX) || !WaveExists(R3dX) || !WaveExists(Intens3D))
		DoAlert 0, "the waves R3dX, R3dH, R3dF, or TotalIntensity3d do not exists in the current folder"
	endif
	WaveStats/M=1/Q Intens3D
	Variable level = V_max*0.5
	if (stringmatch(GetRTStackInfo(0),"*PopMenuProc;TrimTopSurface3d;"))
		printf "¥TrimTopSurface3d()				//  set surface at intensity = %g\r",level
	else
		printf "set surface at intensity = %g\r",level
	endif

	Variable NX=DimSize(R3dX,0), NH=DimSize(R3dX,1), NF=DimSize(R3dX,2)
	Make/N=(NF)/O/D top_temp_
	Wave wtop = top_temp_
	SetScale/P x,0,1,"",wtop


	Make/N=(NX,NH)/O/D surfaceHeight=NaN
	SetScale/P x,DimOffset(R3dX,0),DimDelta(R3dX,0),WaveUnits(R3dX,0), surfaceHeight
	SetScale/P y,DimOffset(R3dX,1),DimDelta(R3dX,1),WaveUnits(R3dX,1), surfaceHeight
	SetScale d 0,0,WaveUnits(R3dX,2), surfaceHeight
	Variable i,j, m
	for (j=0;j<NH;j+=1)
		for (i=0;i<NX;i+=1)
			wtop = Intens3D[i][j][p]
			WaveStats/Q/M=1 wtop
			if (!(V_max > level))				// never gets high enough
				continue
			endif
			for (m=V_maxloc;m>=0;m-=1)		// search backwards from peak to find level
				if (!(wtop[m]>level))
					break
				endif
			endfor
			if (m<0)							// never got low enough
				continue
			endif
			surfaceHeight[i][j] = DimOffset(R3dX,2) + m*DimDelta(R3dX,2)
			R3dX[i][j][0,m] = NaN
			R3dH[i][j][0,m] = NaN
			R3dF[i][j][0,m] = NaN
		endfor
	endfor
	KillWaves/Z top_temp_
	return 0
End
//
//Function topSurface()
//	Wave XX=XX, HH=HH, FF=FF, totalIntensity=totalIntensity
//
//	Variable Xo=XX[0]
//	Variable i0,i1
//	Variable m, i, N=numpnts(XX)
//	Variable level
//
//
//	for (i1=0; XX[i1]==Xo && i1<N; i1+=1)	//find i1, end point
//	endfor
//	N = i1
//
//	Make/N=(100)/O Hsurface=NaN,Fsurface=NaN
//	SetScale d 0,0,"µm", Fsurface,Hsurface
//	i0 = 0
//	m = 0
//	do
//		for (i1=i0; FF[i1+1]>FF[i1] && i1<N; i1+=1)	//find i1, end point
//		endfor
//
//		WaveStats/M=1/Q/R=[i0,i1]/Z totalIntensity
////		if (i0==0)
//			level = (totalIntensity[i0]+V_max)/2		// find level between start and V_max in [0,V_maxloc]
////			printf "looking for %s = %g\r",NameOfWave(totalIntensity),level
////		endif
//
//		FindLevel/P/Q/R=[i0,V_maxloc] totalIntensity, level
//		if (!V_flag)
//			if (m>=numpnts(Hsurface))
//				Redimension/N=(m+50) Hsurface,Fsurface
//			endif
//			Hsurface[m] = HH(V_LevelX)
//			Fsurface[m] = FF(V_LevelX)
////			printf "in [%d, %g] at (H,F) = (%g, %G)\r",i0,i1,HH(V_LevelX),FF(V_LevelX)
//			m += 1
//		else
//			printf "FindLevel failed for range [%d, %d], max is only %g\r",i0,i1,V_max
//		endif
//		i0 = i1+1
//	while(XX[i0]==Xo)
//	Redimension/N=(m) Hsurface,Fsurface
//End














Function Write3dOrientsFile(IndexingList,irefList)
	Wave/T IndexingList						// list of indexing information
	String irefList								// list of iref values to use for reference matricies (usually only one value is needed)
	if (!WaveExists(IndexingList) || !(str2num(irefList)>=0))
		String IndexingListName = SelectString(WaveExists(IndexingList),"",NameOfWave(IndexingList))
		Prompt IndexingListName,"Text Wave with Indexations",popup,WaveListClass("indexationResultList","*","TEXT:1")
		irefList = SelectString(!(str2num(irefList)>=0),irefList,"1")
		Prompt irefList,"list of indecis into IndexingList for the reference directions (usually only one value)"
		DoPrompt "",IndexingListName,irefList
		if (V_flag)
			return 1
		endif
		printf "¥Write3dOrientsFile(%s,\"%s\")\t\t\tdata folder  =  %s\r",IndexingListName,irefList,GetDataFolder(1)
		Wave/T IndexingList=$IndexingListName
	endif
	irefList = ReplaceString(",",irefList,";")					// change all commas -> semi-colons, so irefList can be separated by comma, semi-colon, or spaces
	irefList = ReplaceString("  ",irefList," ")					// change all double spaces to single
	irefList = ReplaceString(" ",irefList,";")					// change all spaces -> semi-colons
	if (!WaveExists(IndexingList) || !(str2num(irefList)>=0))
		DoAlert 0,"Bad inputs, do nothing"
		return 1
	endif
	Variable N=DimSize(IndexingList,0)
	if (N<1)
		DoAlert 0, NameOfWave(IndexingList)+"is empty, do nothing"
		return 1
	endif
	Variable i, iref, NrefMat=ItemsInList(irefList)						// index into IndexingList[], value from irefList
	for (i=0;i<NrefMat;i+=1)
		iref = str2num(StringFromList(i,irefList))
		if (!(NumberByKey("NpatternsFound",IndexingList[iref],"=")>0))
			DoAlert 0, NameOfWave(IndexingList)+"["+num2str(iref)+"] was not indexed, try again"
			return 1
		endif
	endfor

	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))									//fill the lattice structure with test values
		DoAlert 0, "no crystal structure found"
		return 1
	endif
	initSymmetryOperations()											// initialize all symmetry operations
	Variable useSymmetry												// flag, know how to reduce this symmetry
	useSymmetry = strlen(LatticeSym#MakeSymmetryOps(xtal))>0			// make a wave with the symmetry operation

//	if (xtal.SpaceGroup >=195)											// cubic material
//		String/G root:Packages:Lattices:Rotations:SymmetryOpsPath = "root:Packages:Lattices:Rotations:CubicSymmetryOps"
//		useSymmetry = 1
//	elseif (xtal.SpaceGroup ==limit(xtal.SpaceGroup,75,142))			// tetragonal material
//		String/G root:Packages:Lattices:Rotations:SymmetryOpsPath = "root:Packages:Lattices:Rotations:TetragonalSymmetryOps"
//		useSymmetry = 1
//	elseif (xtal.SpaceGroup ==limit(xtal.SpaceGroup,168,194))		// hexagonal material
//		String/G root:Packages:Lattices:Rotations:SymmetryOpsPath = "root:Packages:Lattices:Rotations:HexagonalSymmetryOps"
//		useSymmetry = 1
//	else																	// triclinic material
//		String/G root:Packages:Lattices:Rotations:SymmetryOpsPath = "root:Packages:Lattices:Rotations:TriclinicSymmetryOps"
//		useSymmetry = 1
//	endif

	Variable X0=0, Y0=0, Z0=0
	Y0 = 23.335
	Z0 = -26.163

	Wave rho = $MakeUnique3x3Mat($"")								// make a new 3x3 matrix to hold rotation from rl to g0 = rho  x rl
	String wName = UniqueName("mat3x3_",1,0)
	Make/N=(3,3,NrefMat)/D $wName
	Wave rhos=$wName
	Make/O/D/N=(3,3,NrefMat) g0
	g0 = (p==q)
	Wave g0Temp = $MakeUnique3x3Mat($"")							// make a new 3x3 matrix to hold rotation from rl to g0 = rho  x rl
	MatrixOp/O gm = Identity(3)*1.0
	Make/N=(3,3)/O/D rot33
	rot33 = (p==q)
	Wave rl = $MakeUnique3x3Mat($"")								// make a new 3x3 matrix to hold the standard reciprocal lattice
	rl[0][0] = xtal.as0 ;	rl[0][1] = xtal.bs0 ;	rl[0][2] = xtal.cs0	// reciprocal lattice in standard orientation (needed for sym reduction)
	rl[1][0] = xtal.as1 ;	rl[1][1] = xtal.bs1 ;	rl[1][2] = xtal.cs1
	rl[2][0] = xtal.as2 ;	rl[2][1] = xtal.bs2 ;	rl[2][2] = xtal.cs2
	for (i=0;i<NrefMat;i+=1)
		iref = str2num(StringFromList(i,irefList))					// use iref to get the refrence lattice
		Variable ax,ay,az,bx,by,bz,cx,cy,cz
		sscanf  StringByKey("rlattice",IndexingList[iref],"="),"{{ %g, %g, %g}{ %g,%g,%g}{%g, %g,%g}}",ax,ay,az,bx,by,bz,cx,cy,cz
		g0[0][0][i]=ax;		g0[1][0][i]=ay;		g0[2][0][i]=az
		g0[0][1][i]=bx;		g0[1][1][i]=by;		g0[2][1][i]=bz
		g0[0][2][i]=cx;		g0[1][2][i]=cy;		g0[2][2][i]=cz
		g0Temp = g0[p][q][i]
		MatrixOp/O rho = g0Temp x Inv(rl)
		rhos[][][i] = rho[p][q]
	endfor

	Variable fid
	Open/C="R*ch"/M="file to save result of indexing (reciprocal lattices)"/P=home/T="TEXT" fid
	if (strlen(S_fileName)<1)
		DoAlert 0, "Cannot open output file"
		KillWaves/Z gm, g0, g0Temp, rot33, rho, rhos, rl
		return 1
	endif
	Variable fTotalPeakIntensity, fTotalIntensity,fepsilon				// flags
	fTotalPeakIntensity = keyInList("totalPeakIntensity",IndexingList[0],"=",";")
	fTotalIntensity = keyInList("totalIntensity",IndexingList[0],"=",";")
	fepsilon = keyInList("exx",IndexingList[0],"=",";")
	String wnote = note(IndexingList)
	printf "writing result list of orientations to file =  '%s'\r",S_fileName
	fprintf fid,"$filetype	ArrayOf3dOrientsFile\r"
	fprintf fid,"$X0	%g				// X(micron) of center for cylindrical symmetry, beam line coordinates\r",X0
	fprintf fid,"$Y0	%g				// Y(micron) of center for cylindrical symmetry\r",Y0
	fprintf fid,"$Z0	%g				// Z(micron) of center for cylindrical symmetry\r",Z0
	fprintf fid,"$dX	0.0				// X component of unit vector in direction of cylindrical axis, beam line coordinates\r"
	fprintf fid,"$dY	0.707107		// Y component of unit vector in direction of cylindrical axis\r"
	fprintf fid,"$dZ	-0.707107		// Z component of unit vector in direction of cylindrical axis\r"
	fprintf fid,"$NrefMat	%d		// number of reference matricies to use\r",NrefMat
	fprintf fid,"$irefList	%s		// iref numbers for each ref matrix\r",ReplaceString(";",irefList,",")
	for (i=0;i<NrefMat;i+=1)
		fprintf fid,"$%s	%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f // reference matrix identifies zero angle, same order as data\r","refMat"+SelectString(i,"",num2istr(i)),g0[0][0][i],g0[0][1][i],g0[0][2][i],g0[1][0][i],g0[1][1][i],g0[1][2][i],g0[2][0][i],g0[2][1][i],g0[2][2][i]
	endfor
	fprintf fid,"//\r"
	fprintf fid,"//     g11 g12 g13        h        x\r"
	fprintf fid,"// g = g21 g22 g23,  and  k = g x  y\r"
	fprintf fid,"//     g31 g32 g33        l        z\r"
	fprintf fid,"//\r"
	if (strlen(StringByKey("imageFilePath",wnote,"=")))
		fprintf fid,"$imageFilePath	%s\r",StringByKey("imageFilePath",wnote,"=")
	endif
	fprintf fid,"$filePrefix	%s\r",StringByKey("filePrefix",wnote,"=")
	fprintf fid,"$positions	%s\r",StringByKey("positions",wnote,"=")
	fprintf fid,"$depths	%s\r",StringByKey("depths",wnote,"=")
	if (strlen(StringByKey("hkl",wnote,"=")))
		fprintf fid,"$hkl	%s\r",StringByKey("hkl",wnote,"=")
		fprintf fid,"$cone	%s\r",StringByKey("cone",wnote,"=")
		fprintf fid,"$keVmaxCalc	%s\r",StringByKey("keVmaxCalc",wnote,"=")
		fprintf fid,"$angleTolerance	%s\r",StringByKey("angleTolerance",wnote,"=")
		fprintf fid,"$minPeakWidth	%s\r",StringByKey("minPeakWidth",wnote,"=")
		fprintf fid,"$maxPeakWidth	%s\r",StringByKey("maxPeakWidth",wnote,"=")
		fprintf fid,"$minSep	%s\r",StringByKey("minSep",wnote,"=")
		fprintf fid,"$threshAboveAvg	%s\r",StringByKey("threshAboveAvg",wnote,"=")
	endif
	if (strlen(StringByKey("strainCoords",wnote,"=")) && fepsilon)
		fprintf fid,"$strainCoords	%s		// coordinate system for the strains",StringByKey("strainCoords",wnote,"=")
	endif
	fprintf fid,"//\r"
	fprintf fid,"// where (x,y,z) are in the X,H,F coordinate system, g is recip matrix\r"
	fprintf fid,"//\r"
	fprintf fid,"//	1st 3 columns X,Y,Z locations of voxel in beamline coordinate system (microns)\r"
	fprintf fid,"//	next 9 columns orientation matrix g11 g12 g13 g21 g22 g23 g31 g32 g33\r"
	fprintf fid,"//\r"
	fprintf fid,"$XYZ_matrix 	X  Y  Z  ReciprocalLattice"
	if (NrefMat>1)
		fprintf fid,"  irefMat"
	endif
	if (fTotalPeakIntensity)
		fprintf fid,"  totalPeakIntensity"
	endif
	if (fTotalIntensity)
		fprintf fid,"  totalIntensity"
	endif
	if (fepsilon)
		fprintf fid,"  exx  eyy  ezz  exy  exz  eyz  a_nm  b_nm  c_nm  alpha  bet  gam"
	endif
	fprintf fid,"\r"

	Variable xx,yy,zz, xo, yo, zo, depth
	Variable totalPeakIntensity, totalIntensity
	Variable exx,eyy,ezz,exy,exz,eyz
	Variable a,b,c,alpha,bet,gam
	String list=IndexingList[0]
	xo = NumberByKey("X1",list,"=")									// position of sample positioner
	yo = NumberByKey("Y1",list,"=")
	zo = NumberByKey("Z1",list,"=")
	Variable angle,minA=Inf,maxA=-Inf, angleMin,irefMin
	for (i=0;i<N;i+=1)
		list = IndexingList[i]
		totalPeakIntensity = NumberByKey("totalPeakIntensity",list,"=")	
		totalIntensity = NumberByKey("totalIntensity",list,"=")	
		depth = NumberByKey("depth",list,"=")							// depth measured long the beam (z-direction)
		if (numtype(depth))
			depth = NumberByKey("depthSi",list,"=")					// try old style tag
		endif
		depth = numtype(depth) ? 0 : depth								// if depth not present, then set it to zero
		xx = -(NumberByKey("X1",list,"=")-xo)
		yy = -(NumberByKey("Y1",list,"=")-yo)
		zz = -(NumberByKey("Z1",list,"=")-zo) + depth
		gm = NaN
		if (NumberByKey("NpatternsFound",list,"=")>0)
			sscanf  StringByKey("rlattice",list,"="),"{{ %g, %g, %g}{ %g,%g,%g}{%g, %g,%g}}",ax,ay,az,bx,by,bz,cx,cy,cz	// three reciprocal lattice vectors
			// this is the reciprocal lattice that was just read in, invert to get matrix to write out
			gm[0][0]=ax;		gm[1][0]=ay;		gm[2][0]=az
			gm[0][1]=bx;		gm[1][1]=by;		gm[2][1]=bz
			gm[0][2]=cx;		gm[1][2]=cy;		gm[2][2]=cz
			for (iref=0,angleMin=Inf; iref<NrefMat; iref+=1)			// find irefMin
				rho[][] = rhos[p][q][iref]
				if (useSymmetry)										// symmetry reduce this structure
					angle = symReducedRecipLattice(rl,rho,gm,rot33)	// angle between go and gm
				else
					angle = rotationAngleOfMat(rho)					// do not yet know how to symmetry reduce other systems
				endif
				if (angle < angleMin)
					irefMin = iref
					angleMin = angle
				endif
			endfor														// after this loop, know irefMin
			rho = rhos[p][q][irefMin]									// set rho to the best value
			if (useSymmetry)											// symmetry reduce this structure
				angle = symReducedRecipLattice(rl,rho,gm,rot33)		// angle between go and gm
			else
				MatrixOp/O rot33 = gm x Inv(rho x rl)	// r = gm x rl^-1 x rho^-1
				angle = rotationAngleOfMat(rot33)						// do not yet know how to symmetry reduce other systems
			endif
			maxA = max(maxA,angle)
			minA = min(minA,angle)
			g0Temp = g0[p][q][irefMin]
			MatrixOp/O gm = rot33 x g0Temp
			exx = NumberByKey("exx",list,"=")	
			eyy = NumberByKey("eyy",list,"=")	
			ezz = NumberByKey("ezz",list,"=")	
			exy = NumberByKey("exy",list,"=")	
			exz = NumberByKey("exz",list,"=")	
			eyz = NumberByKey("eyz",list,"=")	
			a = NumberByKey("a",list,"=")	
			b = NumberByKey("b",list,"=")	
			c = NumberByKey("c",list,"=")	
			alpha = NumberByKey("alpha",list,"=")	
			bet = NumberByKey("bet",list,"=")	
			gam = NumberByKey("gam",list,"=")	
		else
			exx = NaN;		eyy = NaN;		ezz = NaN
			exy = NaN;		exz = NaN;		eyz = NaN
			a = NaN;		b = NaN;		c = NaN
			alpha = NaN;	bet = NaN;		gam = NaN
		endif
		fprintf fid,"%g %g %g %g %g %g %g %g %g %g %g %g"xx,yy,zz,gm[0][0],gm[0][1],gm[0][2],gm[1][0],gm[1][1],gm[1][2],gm[2][0],gm[2][1],gm[2][2]
		if (NrefMat>1)
			fprintf fid,"  %d",irefMin
		endif
		if (fTotalPeakIntensity)
			fprintf fid,"  %g",totalPeakIntensity
		endif
		if (fTotalIntensity)
			fprintf fid,"  %g",totalIntensity
		endif
		if (fepsilon)
			fprintf fid,"  %g %g %g %g %g %g  %.8g %.8g %.8g %.8g %.8g %.8g",exx,eyy,ezz,exy,exz,eyz,a,b,c,alpha,bet,gam
		endif
		fprintf fid,"\r"
	endfor
	Close fid
	printf "rotation angles range over [%g¡, %g¡]\r",minA,maxA
	KillWaves/Z gm, g0, g0Temp, rot33, rho, rhos, rl
	return 0
End	


// returns rot, rotation matrix from rho to B, which is symmetry reduced to smallest total angle
// added the possibility of twin  operations (which should be nested inside the symmetry operations)
Function symReducedRecipLattice(rl,rho,B,rot)
	Wave rl								// standard recip lattice
	Wave rho								// A = rho x rl, rotation from standard (= identity) to reference orientation (=A)
	Wave B								// final recip lattice
	Wave rot								// resulting symmetry reduced rotation matrix

	String wName
	Wave roti = root:Packages:Lattices:Rotations:rotMat
	Wave symOp = root:Packages:Lattices:Rotations:SymMatrix
	Wave symOps=$StrVarOrDefault("root:Packages:Lattices:SymOps:SymmetryOpsPath","")	// wave with symmetry ops
	Variable N = NumberByKey("Nproper",note(symOps),"=")	// number of symmetyry operations to check
	N = numtype(N) ? DimSize(symOps,0) : N
	Variable i, trace, traceMax = -4
	for (i=0;i<N;i+=1)
		symOp = symOps[i][p][q]									// symmetry operation
		MatrixOp/O roti = B x Inv(rho x symOp x rl)				// r = gm x rl^-1 x sym-1 x rho^-1
		trace = MatrixTrace(roti)
		if (trace > traceMax)											// save value for optimum rotation matrix
			traceMax = trace
			rot = roti
		endif
	endfor

	// now process the twins if needed
	//
	// I am assuming that the order of multiplication "twinOp x symOp" or "symOp x twinOp" does not matter
	String twinOpsPathList = StrVarOrDefault("root:Packages:Lattices:Rotations:twinOpsPathList","")	// list of path to wavew with twin operations
	Variable NtwinTypes = ItemsInList(twinOpsPathList)
	if (NtwinTypes>0)
		Wave twinOp=root:Packages:Lattices:Rotations:twinMatrix
		Variable j,m,Nt												// number of twin operations in one twinOps
		for (m=0;m<NtwinTypes;m+=1)							// loop over list of operations
			Wave twinOps=$StringFromList(m,twinOpsPathList)
			Nt = WaveExists(twinOps) ? DimSize(twinOps,0) : 0
			for (j=0;j<Nt;j+=1)
				twinOp = twinOps[j][p][q]							// twin operation
				for (i=0;i<N;i+=1)
					symOp = symOps[i][p][q]						// symmetry operation
					MatrixOp/O roti = B x Inv(rho x twinOp x symOp x rl)	// r = gm x rl^-1 x sym-1 x rho^-1
					trace = MatrixTrace(roti)
					if (trace > traceMax)								// save value for optimum rotation matrix
						traceMax = trace
						rot = roti
					endif
				endfor
			endfor
		endfor
	endif

	Variable cosine = (traceMax-1)/2
	cosine = max(min(cosine,1),-1)
	Variable angle = acos(cosine)*180/PI
	return angle
End
//
//// returns rot, rotation matrix from rho to B, which is symmetry reduced to smallest total angle
//Function symReducedRecipLattice(rl,rho,B,rot)
//	Wave rl								// standard recip lattice
//	Wave rho								// A = rl x rho, rotation from standard to reference
//	Wave B								// final recip lattice
//	Wave rot								// resulting symmetry reduced rotation matrix
//
//	String wName
//	Wave roti = root:Packages:Lattices:Rotations:rotMat
//	Wave symOp = root:Packages:Lattices:Rotations:SymMatrix
//	Wave symOps=$StrVarOrDefault("root:Packages:Lattices:Rotations:SymmetryOpsPath","")	// wave with symmetry ops
//	Variable N = NumberByKey("Nproper",note(symOps),"=")	// number of symmetyry operations to check
//	N = numtype(N) ? DimSize(symOps,0) : N
//	Variable i, trace, traceMax = -4
//	for (i=0;i<N;i+=1)
//		symOp = symOps[i][p][q]									// symmetry operation
//
////		MatrixOp/O roti = B x Inv(rl) x Inv(symOp) x Inv(rho)	// r = gm x rl^-1 x sym-1 x rho^-1
//		MatrixOp/O roti = B x Inv(rho x symOp x rl)	// r = gm x rl^-1 x sym-1 x rho^-1
//		trace = MatrixTrace(roti)
//		if (trace > traceMax)											// save value for optimum rotation matrix
//			traceMax = trace
//			rot = roti
////print "\r***"
//		endif
////printf"\t% d  % d  % d	trace =%g  -> %g¡\r",symOp[0][0],symOp[0][1],symOp[0][2],trace, rotationAngleOfMat(roti)
////printf"\t% d  % d  % d     angle=%g\r",symOp[1][0],symOp[1][1],symOp[1][2],rotationAngleOfMat(symOp)
////printf"\t% d  % d  % d\r",symOp[2][0],symOp[2][1],symOp[2][2]
//	endfor
//	Variable cosine = (traceMax-1)/2
//	cosine = max(min(cosine,1),-1)
//	Variable angle = acos(cosine)*180/PI
//	return angle
//End


//
// returns rot, rotation matrix from A to B, which is symmetry reduced to smallest total angle
Static Function symReducedRotation(A,B,rot)
	Wave A,B									// two orientation matricies
	Wave rot									// resulting symmetry reduced rotation matrix
	String wName

	Wave BASi = root:Packages:Lattices:Rotations:rotMat
	Wave symOp = root:Packages:Lattices:Rotations:SymMatrix
	Wave symOps=$StrVarOrDefault("root:Packages:Lattices:SymOps:SymmetryOpsPath","")	// wave with symmetry ops
	Variable N = NumberByKey("Nproper",note(symOps),"=")	// number of symmetyry operations to check
	N = numtype(N) ? DimSize(symOps,0) : N
	Variable i, trace, traceMax = -4
	for (i=0;i<N;i+=1)
		symOp = symOps[i][p][q]									// symmetry operation
		MatrixOp/O BASi = B x Inv(A) x symOp						// B(A^-1)(S^i)
		trace = MatrixTrace(BASi)
		if (trace > traceMax)											// save value for optimum rotation matrix
			traceMax = trace
			rot = BASi
		endif
	endfor
	Variable cosine = (traceMax-1)/2
	cosine = max(min(cosine,1),-1)
	Variable angle = acos(cosine)*180/PI
	return angle
End


Function initSymmetryOperations()
	NewDataFolder/O root:Packages
	NewDataFolder/O root:Packages:Lattices:Rotations
	Make/N=(3,3)/D/O root:Packages:Lattices:Rotations:SymMatrix,root:Packages:Lattices:Rotations:rotMat
	Make/N=(3,3)/D/O root:Packages:Lattices:Rotations:TwinMatrix
End
//Function initSymmetryOperations()
//	Variable doit = 0
//	doit = (exists("root:Packages:Lattices:Rotations:SymmetryOpsPath")!=2) ? 1 : doit
//	if (!doit)
//		return 0
//	endif
//	NewDataFolder/O root:Packages
//	NewDataFolder/O root:Packages:Lattices:Rotations
////	MakeCubicSymmetryOps("root:Packages:Lattices:Rotations:")
////	MakeTetragonalSymmetryOps("root:Packages:Lattices:Rotations:")
////	MakeHexagonalSymmetryOps("root:Packages:Lattices:Rotations:")
////	MakeTriclinicSymmetryOps("root:Packages:Lattices:Rotations:")
//	Make/N=(3,3)/D/O root:Packages:Lattices:Rotations:SymMatrix,root:Packages:Lattices:Rotations:rotMat
//	Make/N=(3,3)/D/O root:Packages:Lattices:Rotations:TwinMatrix
//End


Function EnableBCC_Twins()
	String/G root:Packages:Lattices:Rotations:twinOpsPathList
	SVAR twinOpsPathList=root:Packages:Lattices:Rotations:twinOpsPathList
	String str = "root:Packages:Lattices:Rotations:BCCTwinOps"
	if (WhichListItem(str,twinOpsPathList)<0)
		twinOpsPathList = AddListItem(str,twinOpsPathList)
	endif
	MakeBCCtwinOps("root:Packages:Lattices:Rotations:")
	printf "\ralso will check the symmetry operations for BCC twins\r"
	printf "		checking for 70¡ rotations about 011 types\r\r"
End
//
Function DisableBCC_Twins()
	SVAR twinOpsPathList=root:Packages:Lattices:Rotations:twinOpsPathList
	if (SVAR_Exists(twinOpsPathList))
		twinOpsPathList = RemoveFromList("root:Packages:Lattices:Rotations:BCCTwinOps",twinOpsPathList)
	endif
	printf "\rno longer checking for BCC twins\r\r"
End

Function/S MakeBCCtwinOps(fldr)				// make a wave with the BCC twin operations, rotations about (011), ±70¡
	String fldr										// folder used to put the symmetryOperations

	Make/N=(3,3)/D/O root:Packages:Lattices:Rotations:TwinMatrix
	String wName = fldr+"BCCTwinOps"
	Make/N=(1+3,3,3)/D/O $wName=0	// there are 3 unique (011) types, and each gets 1 rotation + identity = 4
	Wave ops = $wName
	wName = UniqueName("mat",1,0)
	Make/N=(3,3)/D/O $wName
	Wave mat=$wName
	wName = UniqueName("axis",1,0)
	Make/N=3/D/O $wName
	Wave axis=$wName

	Variable angle=acos(1/3)*180/PI			// gives 70.529¡

	// the identity, 1 operation
	ops[0][][] = (q==r)

	// the 1 rotation about each of the (110), adds 1 operations
	axis = {1,1,0}
	rotationMatAboutAxis(axis,angle,mat)		// 70.529¡ about x-axis
	ops[1][][] = mat[q][r]

	// the 1 rotation about each of the (110), adds 1 operation
	axis = {1,0,1}
	rotationMatAboutAxis(axis,angle,mat)		// 70.529¡ about x-axis
	ops[2][][] = mat[q][r]

	// the 1 rotation about each of the (110), adds 1 operation
	axis = {0,1,1}
	rotationMatAboutAxis(axis,angle,mat)		// 70.529¡ about x-axis
	ops[3][][] = mat[q][r]

	KillWaves/Z mat,axis
	Note/K ops
	Note ops, "Nproper="+num2istr(DimSize(ops,0))+";"	// alll of the operations are proper
	return GetWavesDataFolder(ops,2)
End


//	WH_2890_39   is at		Y1 = -3691.05;	Z1 = -5095.45
//	from point 0,  				yo = -3682.05;	zo = -5086.5


// find the entry in indexingList[] that is the closest contributor to this (H,F) point
Function findClosestPoint(H,F)
	Variable H,F
	Variable printIt = (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(GetRTStackInfo(0),"*ButtonProc;findClosestPoint;"))
	String wName=""
	if (numtype(H+F) && printIt)
		String c,csrList=""
		Variable k
		for (k=char2num("A");k<=char2num("J");k+=1)
			c = num2char(k)
			if (strlen(csrInfo($c))>0)
				csrList += c+";"
			endif
		endfor
		if (strlen(csrList)<1)
			DoAlert 0, "no inputs and no cursors on image"
			ShowInfo
			return NaN
		endif
		c = StringFromLIst(0,csrList)
		if (ItemsInList(csrList)>1)		// more than one cursor on graph, so ask
			Prompt c,"cursor to use",popup,csrList
			DoPrompt "cursor",c
			if (V_flag)
				return NaN
			endif
		endif
		H = hcsr($c)
		F = vcsr($c)
		wName = CsrWave($c,"",1)
	endif
	if (numtype(H+F))
		if (printIt)
			DoAlert 0, "no inputs and no cursors on image"
			ShowInfo
		endif
		return NaN
	endif

	Variable X=0
	String fldr= GetDataFolder(1)
	if (strlen(wName))
		fldr = GetWavesDataFolder(ImageNameToWaveRef("",wName),1)
	endif
	if (!stringmatch(StrVarOrDefault(fldr+"SliceNormal","X"),"X"))
		if (printIt)
			DoAlert 0, "not a cut at constant X"
		endif
		return NaN
	endif
	X = NumVarOrDefault(fldr+"SliceValue",0)

	Wave XX=XX, HH=HH, FF=FF, indicies=indicies
	if (!WaveExists(XX) || !WaveExists(HH) || !WaveExists(FF))
		DoAlert 0,"Could not find XX, HH or FF"
		return NaN
	endif
	Wave/T indexingList=::IndexingList
	if (!WaveExists(indexingList))
		Wave/T indexingList=:::IndexingList
	endif

	Variable i, imin=-1, minDist, dist, N=numpnts(FF)
	for (i=0,minDist=Inf; i<N; i+=1)
		dist = (X-XX[i])^2 + (H-HH[i])^2 + (F-FF[i])^2
		if (dist<minDist)
			imin = i
			minDist = dist
		endif
	endfor
	minDist = sqrt(minDist)
	if (imin>=N || imin<0)
		if (printIt)
			print "ERROR, could not find closest point"
		endif
		return NaN
	endif
	imin = WaveExists(indicies) ? indicies[imin] : imin

	String imageName=""
	if (WaveExists(indexingList))
		ImageName = StringByKey("image",IndexingList[imin],"=")
	elseif (!WaveExists(indexingList) && exists("imageNames")==1)
		Wave/T imageNames=imageNames
		imageName = imageNames[imin]
	endif

	if (printIt)
		printf "for point (%g, %g, %gµm)%s,    closest at iref=%d,   dist=%g   for  '%s'\r",X,H,F,SelectString(strlen(imageName),""," on "+wName),imin,minDist,imageName
		if (WaveExists(indexingList))
			print (IndexingList[imin])[0,140]+"  . . ."
		endif
	endif
	return imin
End




// put a marker on the graph at the ref points in list, if list==FF, then use irefList from wave note
Function putMarkersAtRefPoints(list)
	String list

	Wave HH=HH, FF=FF
	if (!WaveExists(HH) || !WaveEXists(FF))
		DoAlert 0,"cannot find HH or FF, probably in wrong data folder"
		return 1
	endif
	if (stringmatch(list,"FF"))					// for FF, get irefList from the wave note
		list = StringByKey("irefList", note(FF),"=")
	endif
	list = ReplaceString("  ",list," ")			// convert space or comma separated lists to semi-colon
	list = ReplaceString(" ",list,";")
	list = ReplaceString(",",list,";")

	String win = StringFromList(0,WinList("*",";", "WIN:1" ))
	if (strlen(win)<1)
		return 1
	endif
	Variable i,iref
	for (i=0;i<ItemsInList(list);i+=1)
		iref = (str2num(StringFromList(i,list)))
		printf "[iref==%d]       \tmarker at (H,F) = (%.3g,  %.3g)\r",iref,HH[iref],FF[iref]
		DrawMarker(HH[iref],FF[iref],1.5,1.5,"x",color="black",thick=2,win=win)
	endfor
	return 0
End





//
//	Variable depthSi = NumberByKey("depthSi",list,"=")	// depth measured long the beam (z-direction)
//	Variable xo,yo, zo
//	xo = NumberByKey("X1",IndexingList[0],"=")			// position of first point
//	yo = NumberByKey("Y1",IndexingList[0],"=")
//	zo = NumberByKey("Z1",IndexingList[0],"=")
//	Variable xv = -(NumberByKey("X1",list,"=")-xo)
//	Variable yy = -(NumberByKey("Y1",list,"=")-yo)
//	Variable zz = -(NumberByKey("Z1",list,"=")-zo) + depthSi
//
//	Variable Y0=20.5, Z0=23.33
//
//	Variable X0=0
//	Variable H0=(Y0+Z0)/sqrt(2), F0=(-Y0+Z0)/sqrt(2)
//	Variable H1, F1, X1
//	X1 = xv - X0
//	H1 = (yy+zz)/sqrt(2) - H0
//	F1 = (-yy+zz)/sqrt(2) - F0


//// now using the one from depthResolve.ipf
//	Static Function/T GetFileRootAndDoubleRange(pathName,filePrefix,positions,depths)
//		String pathName					// name of path to use
//		String filePrefix					// first part of file name, e.g. "WH_"
//		String positions, depths			// string ranges
//	
//		String fileRoot=""
//		Variable printIt=0
//		if (strlen(pathName)<1 || strlen(positions)<1 || strlen(depths)<1)
//			pathName = SelectString(strlen(pathName),"reconPath",pathName)
//			Prompt pathName,"name of path that points to the reconstructed images"
//			Prompt positions,"range of positions, the first index (non-existant ones are skipped)"
//			Prompt depths, "range of depths, the second index"
//			DoPrompt "index ranges",pathName,positions,depths
//			if (V_flag)
//				return ""
//			endif
//			printIt = 1
//		endif
//		pathName = SelectString(stringmatch(pathName,CleanupName(pathName,0)),"",pathName)
//		if (!strlen(pathName) || ItemsInRange(positions)<1 || ItemsInRange(depths)<1)
//			DoAlert 0, "nothing to do, positions='"+positions+"',   depths='"+depths+"'  path='"+pathName+"'"
//			return ""
//		endif
//	
//		PathInfo $pathName
//		String path = S_path
//		if (strlen(filePrefix)<1 || !V_flag || strlen(S_path)<1)	// un-assigned path, or no file prefix
//			Variable refNum
//			Open/D/M="pick any one of the images"/R/T=".spe" refNum
//			if (strlen(S_fileName)<1)
//				DoAlert 0,"no file specified"
//				return ""
//			endif
//			filePrefix = ParseFilePath(3,S_fileName, ":", 1, 0)
//			path = ParseFilePath(1,S_fileName,":",1,0)
//			NewPath/O/Z $pathName,path
//			filePrefix = filePrefix[0,strsearch(filePrefix,"_",Inf,1)-1]
//			filePrefix = filePrefix[0,strsearch(filePrefix,"_",Inf,1)]
//			if (!strlen(filePrefix))
//				DoAlert 0,"unable to identify file prefix"
//				return ""
//			endif
//			printIt = 1
//		endif
//		String list = ReplaceStringByKey("fileRoot","",path+filePrefix,"=")
//		list = ReplaceStringByKey("positions",list,positions,"=")
//		list = ReplaceStringByKey("depths",list,depths,"=")
//		list = ReplaceStringByKey("pathName",list,pathName,"=")
//		list = ReplaceStringByKey("filePrefix",list,filePrefix,"=")
//		list = ReplaceNumberByKey("printIt",list,printIt,"=")
//		return list
//	End






Function/T Add2XmlFile(path,fname,xmlStr)
	String path
	String fname
	String xmlStr

	if (strlen(xmlStr)<1)
		return ""
	endif

	String buf
	Variable f
	GetFileFolderInfo/P=$path/Q/Z=1 fname
	if (!V_isfile || V_isFolder)
		Open/D/M="New XML file for output"/P=$path/T=".xml"/Z=2 f as fname
		if (V_flag)
			return ""
		endif
		fname = S_fileName
	endif

	Open/A/P=$path/Z=1 f as fname
	if (strlen(S_fileName)<1)
		return ""
	endif

	FStatus f
	if (V_logEOF<10)
		buf = "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n\n<AllSteps xmlns=\"http://sector34.xray.aps.anl.gov/34ide:indexResult\">\n\n"
		FSetPos f, 0
		FBinWrite f, buf
	else
		FSetPos f, V_logEOF-12
	endif
	FBinWrite f, xmlStr
	buf = "\n</AllSteps>\n"
	FBinWrite f, buf
	Close f
	return fname
End
//
Static Function/T CreateXmlStepStr(FullPeakList,FullPeakIndexed)
	Wave FullPeakList,FullPeakIndexed
	if (!WaveExists(FullPeakList))
		return ""
	endif
	Variable Npeaks=DimSize(FullPeakList,0), val
	if (Npeaks<1)
		return ""
	endif
	String pnote=note(FullPeakList), str

	String xml="<step>\n"
	xml += addXMLstrFromTagVals("title",pnote,"",1)
	xml += addXMLstrFromTagVals("sampleName",pnote,"",1)
	xml += addXMLstrFromTagVals("userName",pnote,"",1)
	xml += addXMLstrFromTagVals("beamLine",pnote,"beamline",1)
	xml += addXMLstrFromTagVals("scanNum",pnote,"",1)
	xml += addXMLstrFromTagVals("file_time",pnote,"date",1)
	xml += addXMLstrFromTagVals("BeamBad",pnote,"beamBad",1)
	xml += addXMLstrFromTagVals("CCDshutter",pnote,"",1)
	xml += addXMLstrFromTagVals("LightOn",pnote,"lightOn",1)
	xml += addXMLstrFromTagVals("MonoMode",pnote,"monoMode",1)

	Variable X1=NumberByKey("X1",pnote,"="), Y1=NumberByKey("Y1",pnote,"="), Z1=NumberByKey("Z1",pnote,"=")
	if (numtype(X1+Y1+Z1)==0)
		xml += "\t<Xsample>"+num2str(X1)+"</Xsample>\n"
		xml += "\t<Ysample>"+num2str(Y1)+"</Ysample>\n"
		xml += "\t<Zsample>"+num2str(Z1)+"</Zsample>\n"
	endif

	xml += addXMLnumFromTagVals("depth",pnote,"",1,9)
	xml += addXMLnumFromTagVals("keV",pnote,"energy",1,7,unit="keV")
	xml += addXMLnumFromTagVals("HutchTemperature",pnote,"hutchTemperature",1,4)
	xml += addXMLnumFromTagVals("sampleDistance",pnote,"",1,9)

	xml += "\t<detector>\n"			// do the detector
	xml += addXMLstrFromTagVals("imageFileName",pnote,"inputImage",2)
	xml += addXMLstrFromTagVals("detectorID",pnote,"",2)
	xml += addXMLnumFromTagVals("exposure",pnote,"",2,9,unit="sec")
	xml += addXMLnumFromTagVals("xDimDet",pnote,"Nx",2,9)
	xml += addXMLnumFromTagVals("yDimDet",pnote,"Ny",2,9)

	xml += addXMLnumFromTagVals("totalIntensity",pnote,"totalSum",2,9)
	xml += addXMLnumFromTagVals("totalPeakIntensity",pnote,"sumAboveThreshold",2,9)
	xml += addXMLnumFromTagVals("NumPixelsAboveThreshold",pnote,"numAboveThreshold",2,9)

	Variable startx=NumberByKey("startx",pnote,"="), endx=NumberByKey("endx",pnote,"="), groupx=NumberByKey("groupx",pnote,"=")
	Variable starty=NumberByKey("starty",pnote,"="), endy=NumberByKey("endy",pnote,"="), groupy=NumberByKey("groupy",pnote,"=")
	if (numtype(startx+endx+groupx + starty+endy+groupy)==0)
		sprintf str "\t\t<ROI startx=\"%d\" endx=\"%d\" groupx=\"%d\" starty=\"%d\" endy=\"%d\" groupy=\"%d\"></ROI>",startx,endx,groupx , starty,endy,groupy
		xml += str + "\n"
	endif

	str = "<peaksXY"
	val = NumberByKey("minPeakWidth",pnote,"=")
	str += SelectString(numtype(val), " minwidth=\""+num2str(val)+"\"", "")
	val = NumberByKey("threshAboveAvg",pnote,"=")
	str += SelectString(numtype(val), " threshold=\""+num2str(val)+"\"", "")
	val = NumberByKey("maxRfactor",pnote,"=")
	str += SelectString(numtype(val), " maxRfactor=\""+num2str(val)+"\"", "")
	val = NumberByKey("maxPeakWidth",pnote,"=")
	str += SelectString(numtype(val), " maxwidth=\""+num2str(val)+"\"", "")
	val = NumberByKey("boxsize",pnote,"=")
	str += SelectString(numtype(val), " boxsize=\""+num2str(val)+"\"", "")
	str += " Npeaks=\""+num2str(Npeaks)+"\">"
	xml += "\t\t"+str+"\n"

	Make/N=(Npeaks)/FREE ww
	ww = FullPeakList[p][0]
	sprintf str,"<Xpixel>%s</Xpixel>", vec2str(ww,maxPrint=1000,bare=1,sep=" ")
	xml += "\t\t\t"+str+"\n"

	ww = FullPeakList[p][1]
	sprintf str,"<Ypixel>%s</Ypixel>", vec2str(ww,maxPrint=1000,bare=1,sep=" ")
	xml += "\t\t\t"+str+"\n"

	ww = FullPeakList[p][10]
	sprintf str,"<Integral>%s</Integral>", vec2str(ww,maxPrint=1000,bare=1,sep=" ")
	xml += "\t\t\t"+str+"\n"

	if (DimSize(FullPeakList,1)>11)
		ww = FullPeakList[p][11]
		sprintf str,"<Intens>%s</Intens>", vec2str(ww,maxPrint=1000,bare=1,sep=" ")
		xml += "\t\t\t"+str+"\n"
	endif

	ww = FullPeakList[p][4]
	ww /= 2										// change FWHM to HWHM
	sprintf str,"<hwhmX unit=\"pixel\">%s</hwhmX>", vec2str(ww,maxPrint=1000,bare=1,sep=" ")
	xml += "\t\t\t"+str+"\n"

	ww = FullPeakList[p][5]
	ww /= 2										// change FWHM to HWHM
	sprintf str,"<hwhmY unit=\"pixel\">%s</hwhmY>", vec2str(ww,maxPrint=1000,bare=1,sep=" ")
	xml += "\t\t\t"+str+"\n"

	ww = FullPeakList[p][8]
	sprintf str,"<tilt unit=\"degree\">%s</tilt>", vec2str(ww,maxPrint=1000,bare=1,sep=" ")
	xml += "\t\t\t"+str+"\n"

	ww = FullPeakList[p][9]
	sprintf str,"<chisq>%s</chisq>", vec2str(ww,maxPrint=1000,bare=1,sep=" ")
	xml += "\t\t\t"+str+"\n"

	xml += "\t\t</peaksXY>\n"
	xml += "\t</detector>\n"

	// do the indexed result
	if (WaveExists(FullPeakIndexed))
		String inote=note(FullPeakIndexed), hklPrefer, mStr, unit, runit=""
		Variable Npattern,m
		unit = StringByKey("lengthUnit",inote,"=")
		if (strlen(unit))
			runit = " unit=\"1/"+unit+"\""
			unit = " unit=\""+unit+"\""
		endif

		str = "\t<indexing"
		val = NumberByKey("Nindexed",inote,"=")
		str += SelectString(numtype(val), " Nindexed=\""+num2str(val)+"\"", "")
		str += " Npeaks=\""+num2str(Npeaks)+"\""

		Npattern = NumberByKey("NpatternsFound",inote,"=")
		str += SelectString(numtype(Npattern), " Npatterns=\""+num2str(Npattern)+"\"", "")

		val = NumberByKey("keVmaxCalc",inote,"=")
		str += SelectString(numtype(val), " keVmaxCalc=\""+num2str(val)+"\"", "")

		val = NumberByKey("keVmaxTest",inote,"=")
		str += SelectString(numtype(val), " keVmaxTest=\""+num2str(val)+"\"", "")

		val = NumberByKey("angleTolerance",inote,"=")
		str += SelectString(numtype(val), " angleTolerance=\""+num2str(val)+"\"", "")

		val = NumberByKey("cone",inote,"=")
		str += SelectString(numtype(val), " cone=\""+num2str(val)+"\"", "")

		hklPrefer = StringByKey("hklPrefer",inote,"=")
		hklPrefer = ReplaceString(",",hklPrefer," ")
		hklPrefer = ReplaceString("  ",hklPrefer," ")
		hklPrefer = ReplaceString("}",hklPrefer,"")
		hklPrefer = ReplaceString("{",hklPrefer,"")
		hklPrefer = ReplaceString("\"",hklPrefer,"")
		hklPrefer = ReplaceString("'",hklPrefer,"")
		str += SelectString(strlen(hklPrefer), "", " hklPrefer=\""+hklPrefer+"\"")
		str += ">\n"
		str += "\t\t<!-- Result of indexing. -->\n"
		xml += str

		Variable NindexedMax = DimSize(FullPeakIndexed,0)
		Make/N=(NindexedMax)/D/FREE wi
		Make/N=3/D/FREE vec
		Variable Nindexed
		for (m=0;m<Npattern;m+=1)
			Redimension/N=(NindexedMax) wi
			wi = FullPeakIndexed[p][0][m]
			WaveStats/Q/M=1 wi
			Nindexed = V_npnts
			Redimension/N=(Nindexed) wi
			mStr = num2istr(m)
			str = "\t\t<pattern num=\""+mStr+"\""
			val = NumberByKey("rms_error"+mStr,inote,"=")
			str += SelectString(numtype(val), " rms_error=\""+num2str(val)+"\"", "")
			val = NumberByKey("goodness"+mStr,inote,"=")
			str += SelectString(numtype(val), " goodness=\""+num2str(val)+"\"", "")
			val = NumberByKey("goodness"+mStr,inote,"=")
			str += SelectString(numtype(val), " Nindexed=\""+num2str(Nindexed)+"\"", "")
			xml += str + ">\n"

			str = StringByKey("recip_lattice"+mStr,inote,"=")
			Wave w33 = matString2mat33(str)
			xml += "\t\t\t<recip_lattice"+runit+">\n"
			vec = w33[p][0]
			xml += "\t\t\t\t<astar>"+vec2str(vec,places=9,bare=1,sep=" ")+"</astar>\n"
			vec = w33[p][1]
			xml += "\t\t\t\t<bstar>"+vec2str(vec,places=9,bare=1,sep=" ")+"</bstar>\n"
			vec = w33[p][2]
			xml += "\t\t\t\t<cstar>"+vec2str(vec,places=9,bare=1,sep=" ")+"</cstar>\n"
			xml += "\t\t\t</recip_lattice>\n"

			xml += "\t\t\t<hkl_s>\n"
			wi = FullPeakIndexed[p][3][m]
			xml += "\t\t\t\t<h>"+vec2str(wi,places=9,bare=1,sep=" ")+"</h>\n"
			wi = FullPeakIndexed[p][4][m]
			xml += "\t\t\t\t<k>"+vec2str(wi,places=9,bare=1,sep=" ")+"</k>\n"
			wi = FullPeakIndexed[p][5][m]
			xml += "\t\t\t\t<l>"+vec2str(wi,places=9,bare=1,sep=" ")+"</l>\n"
			xml += "\t\t\t</hkl_s>\n"
			xml += "\t\t</pattern>\n"
		endfor

		xml += "\t\t<xtl>\n"
		xml += addXMLstrFromTagVals("structureDesc",inote,"",3)

		STRUCT crystalStructure xtal
		FillCrystalStructDefault(xtal)
		str = xtal.sourceFile
		if (strlen(str))
			xml += "\t\t\t<xtlFile>" + str + "</xtlFile>\n"
		endif
	
		xml += addXMLnumFromTagVals("SpaceGroup",inote,"",3,3)
		Wave wlat = str2vec(StringByKey("latticeParameters",inote,"="))
		xml += "\t\t\t<latticeParameters"+unit+">"+vec2str(wlat,places=9,bare=1,sep=" ")+"</latticeParameters>\n"
		xml += "\t\t</xtl>\n"

	xml += "\t</indexing>\n"
	endif

	xml += "</step>\n"
	//	print ReplaceString("\n",xml,"\r")
	return xml
End
//
Static Function/T addXMLnumFromTagVals(key,wnote,XMLname,Ntabs,places,[unit])
	String key, wnote
	String XMLname
	Variable Ntabs		// number of leading tabs
	Variable places
	String unit
	unit = SelectString(ParamIsDefault(unit),unit,"")

	Ntabs = numtype(Ntabs) ? 0 : Ntabs
	Ntabs = limit(round(Ntabs),0,20)
	places = numtype(places) || places<1 ? 5 : places
	places = limit(round(places),1,16)
	XMLname = SelectString(strlen(XMLname),key,XMLname)
	Variable val = NUmberByKey(key,wnote,"=")
	if (numtype(val)==0)
		String str, fmt
		if (strlen(unit))
			fmt = "<%s unit=\"%s\">%."+num2istr(places)+"g</%s>\n"
			sprintf str,fmt,XMLname,unit,val,XMLname
		else
			fmt = "<%s>%."+num2istr(places)+"g</%s>\n"
			sprintf str,fmt,XMLname,val,XMLname
		endif
		return PadString("",Ntabs,0x09)+str
	endif
	return ""
End
//
Static Function/T addXMLstrFromTagVals(key,wnote,XMLname,Ntabs,[unit])
	String key, wnote
	String XMLname
	Variable Ntabs		// number of leading tabs
	String unit
	unit = SelectString(ParamIsDefault(unit),unit,"")

	Ntabs = numtype(Ntabs) ? 0 : Ntabs
	Ntabs = limit(round(Ntabs),0,20)
	XMLname = SelectString(strlen(XMLname),key,XMLname)
	String str = StringByKey(key,wnote,"=")
	if (strlen(str))
		if (strlen(unit))
			return PadString("",Ntabs,0x09)+"<"+XMLname+" unit=\""+unit+"\">"+str+"</"+XMLname+">\n"
		else
			return PadString("",Ntabs,0x09)+"<"+XMLname+">"+str+"</"+XMLname+">\n"
		endif
	endif
	return ""
End
