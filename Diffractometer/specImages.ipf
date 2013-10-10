#pragma rtGlobals=1		// Use modern global access method.
#pragma ModuleName=specImage
#pragma version = 0.45
#pragma IgorVersion = 6.2
#include "spec", version>=2.25
#include "Diffractometer", version >=0.26
#include "ImageDisplayScaling", version >=1.86
#include "QspaceVolumesView", version >=1.02
#include "GizmoZoomTranslate", version>=1.35
#include "Masking", version>=1.01

Static Constant hc_keVnm = 1.239841856			// h*c (keV-nm)

// #define USE_EXTRA_PILATUS_SPEC_FILE			// put this in your procedure file to use the funky file extra with extra Pilatus info

Static Constant MAX_Q_SPACE_DIM = 100				// for different value use:    OverRide Constant MAX_Q_SPACE_DIM = 300


Menu "Qspace"
	"Make Qspace Volume from spec Images...",MakeQspaceVolumeSpec("")
	help = {"Make and Fill a 3D Qspace array from a spec scan with a 2D detector."}
	MenuItemIfWaveClassExists("Make a Mask from an image...","speImage*;rawImage*","DIMS:2"), MakeMask($"")
	help = {"Make a mask from an image, the mask can then be used when making a q-space volume."}
	"-"
End

Menu "Data"
	Submenu "spec"
		"Load Pilatuf spec Image...",LoadSpecImage(NaN,NaN)
	End
End

Menu "Load Waves"
	"Load Pilatuf spec Image...",LoadSpecImage(NaN,NaN)
End

Menu "Diffractometer"
	SubMenu "spec Images"
		"Load spec Image...",LoadSpecImage(NaN,NaN)
		help = {"Load one image from a spec scan with 2D detector."}
		"Re Plot spec Image...", rePlotDiffractometerImage(NaN,NaN)
		help = {"Re-Plot one image from a diffractometer scan."}
		"Set spec names for normalizing...",specImage#SetGlobalSpecialSpecNames()
	End
	"-"
End


Static Function AfterFileOpenHook(refNum,file,pathName,type,creator,kind)
	Variable refNum, kind
	String file,pathName,type,creator
	if ((kind==1) || (kind==2))		// an experiment (packed or unpacked)
		init_specImage("")
	endif
	return 0
End
Static Function IgorStartOrNewHook(IgorApplicationNameStr)
	String IgorApplicationNameStr
	init_specImage("")
	return 0
End


//  ============================================================================  //
//  ============================= Start of Qspace Histogram ============================  //

// for single thread operation, add the following line to your main procedure window
//		#define USE_SINGLE_THREAD

#ifdef USE_SINGLE_THREAD

//Function/WAVE MakeQspaceVolumeSpec(scanRange,[mask])
//	String scanRange						// combine all images from every point in each scan from scanRange
//	Wave mask								// used pixels that are 1, ignore pixels that are 0
//
//	String motors = StrVarOrDefault("root:Packages:Diffractometer:DiffractometerAxisNames","")
//	Variable Nm=ItemsInList(motors)
//	STRUCT detectorGeometry d
//	if (FillDetectorStructDefault(d))				//fill the detector structure with current values
//		DoAlert 0, "No Detector Structure, please set one"
//		return $""
//	elseif (Nm<1)
//		DoAlert 0, "Unable to find list of motor names, set Diffractometer Type"
//		return $""
//	endif
//	if (ItemsInRange(scanRange)<1)				// bad or empty range, prompt user
//		Variable last = NumVarOrDefault("root:Packages:spec:lastScan",NaN)
//		scanRange = SelectString(numtype(last),num2str(last),"")
//		String maskName=NameOfWave(mask), maskList=WaveListClass("imageMask","*","DIMS:2")
//		Prompt maskName,"mask to use with all images",popup, "_none_;"+maskList
//		Prompt scanRange,"range of spec scans to process together, e.g. \"1-20\""
//		if (ItemsInList(maskList)>0)
//			DoPrompt "scans",scanRange,maskName
//		else
//			DoPrompt "scans",scanRange
//		endif
//		if (V_flag)
//			return $""
//		endif
//		maskName = SelectString(stringmatch(maskName,"_none_"),maskName,"")
//		Wave mask = $maskName
//		printf "•MakeQspaceVolumeSpec(\"%s\"", scanRange
//		if (strlen(maskName))
//			printf ", mask=%s",maskName
//		endif
//		printf ")\r"
//	endif
//	if (ItemsInRange(scanRange)<1)				// number of scans
//		return $""
//	endif
//
//	if (!WaveExists(mask))						// try again to enter mask
//		String imageFldr, fldrSav
//		maskList = WaveListClass("imageMask","*","DIMS:2")
//		for (i=0;i<ItemsInRange(scanRange);i+=1)
//			sprintf imageFldr,"%sspec%d:images:",StrVarOrDefault("root:Packages:spec:specDataFolder",""),scanNum
//			if (DataFolderExists(imageFldr))
//				SetDataFolder imageFldr
//				maskList += WaveListClass("imageMask","*","DIMS:2")
//				SetDataFolder fldrSav
//			endif
//		endfor
//		fldrSav= GetDataFolder(1)
//		print maskList
//	endif
//
//
//return $""
//
//	print "Using a Single Thread"
//	STRUCT boundingVolume v
//	initBoundingVolumeStruct(v)
//
//	Make/N=(Nm)/D/FREE A
//	Make/N=3/D/FREE q0,q1
//	String wnote, str, mname, imageName
//	Variable m, keV, val
//	Variable scanNum,N								// N counts number of images to process
//	Variable filterTrans,I0, AvgI0=0,NAvgI0=0, scale
//	Variable ip,Np									// for looping over points in one scan
//
//	Variable NimagesExpected=0
//	for (scanNum=str2num(scanRange); numtype(scanNum)==0; scanNum=NextInRange(scanRange,scanNum))
//		NimagesExpected += specInfo(scanNum,"Npoints")
//	endfor
//
//	Wave pxpy=$""
//	String progressWin = ProgressPanelStart("",stop=1,showTime=1,status="part 1 of 2")
//	Variable stop=0, preExists=0
//	for (N=0, scanNum=str2num(scanRange); numtype(scanNum)==0 && stop==0; scanNum=NextInRange(scanRange,scanNum))
//		for (m=0;m<Nm;m+=1)
//			A[m] = specInfo(scanNum,StringFromList(m,motors))
//		endfor
//		keV = specInfo(scanNum,"DCM_energy")
//		keV = numtype(keV) ? 10*hc_keVnm/specInfo(scanNum,"wavelength") : keV
//		Np = specInfo(scanNum,"Npoints")			// number of points in this scan
//		for (ip=0;ip<Np;ip+=1)					// for each point in this scan
//			if (mod(N,50)==0)
//				if (ProgressPanelUpdate(progressWin,N/NimagesExpected*100))
//					stop = 1
//					break
//				endif
//			endif
//			imageName = LoadSpecImage(scanNum,ip,extras="quiet:1")
//			preExists=(strsearch(imageName,"preExists:",0)==0)
//			if (preExists)
//				Wave image = $StringByKey("preExists", imageName)
//			else
//				Wave image = $imageName
//			endif
//			if (!WaveExists(image))
//				sprintf str,"Could not load scan %g,  point %g",scanNum,ip
//				print str
//				DoAlert 0,str
//				return $""
//			endif
//			if (!WaveExists(pxpy))
//				Wave pxpy = diffractometer#MakePixelListFromImage(image)		// first time through, make pixel list (full-chip, unbinned)
//			endif
//			wnote = note(image)
//			if (!preExists)
//				KillWaves/Z image					// only kill those waves that I have just loaded
//			endif
//
//			for (m=0;m<Nm;m+=1)
//				mname = CleanupName(StringFromList(m,motors),0)
//				val = NumberByKey(mname,wnote,"=")
//				A[m] = numtype(val) ? A[m] : val
//			endfor
//			val = NumberByKey("keV",wnote,"=")	;	keV = numtype(val) ? keV : val
//			diffractometer#QofPixel(d,keV,d.Nx/2,d.Ny/2,A,q1)
//			if (N==0)
//				q0 = q1
//			endif
//			extendQvolumeFromImage(d,keV,A,v)
//
//			I0 = NumberByKey("I0",wnote,"=")
//			if (I0>0)
//				AvgI0 += I0	
//				NAvgI0 += 1
//			endif
//			N+=1
//		endfor
//	endfor
//	AvgI0 = (numtype(AvgI0)==0 || AvgI0>0 || NAvgI0>0) ? AvgI0/NAvgI0 : 1
//	Variable seconds1 = SecondsInProgressPanel(progressWin)
//	printf "  the first part took  %s\r",Secs2Time(seconds1,5,0)
//	printf "Found %g images to process,  with <I0> = %g\r",N,AvgI0
//	print "Measured Q-Volume is: ",boundingVolumeStruct2str(v)
//
//	if (WaveExists(mask))							// make a 1D version of mask (matches with pxpy)
//		Duplicate/FREE mask, mask1D
//		Redimension/N=(DimSize(pxpy,0)) mask1D
//	else
//		Wave mask1D = $""
//	endif
//
//	String QspaceNameEnd=num2str(str2num(scanRange))
//	QspaceNameEnd += SelectString(ItemsInRange(scanRange)>1,"","_"+num2str(lastInRange(scanRange)))
//	String QspaceName="Qspace"+QspaceNameEnd
//	Make/N=3/D/FREE q01 = q1-q0				// find dQ
//	Variable dQ = norm(q01)/(N-1)
//print "dQ start = ",dQ
//	dQ = max(dQ,max(max(v.dx,v.dy),v.dz)/MAX_Q_SPACE_DIM)
//print "dQ after = ",dQ
//	Variable Nx=ceil(v.dx / dQ)+1, Ny=ceil(v.dy / dQ)+1, Nz=ceil(v.dz / dQ)+1
//
//	Variable Nqtot = Nx*Ny*Nz
//	Make/N=(Nx,Ny,Nz)/D/O $QspaceName=0
//	Wave Qspace=$QspaceName
//	Make/N=(Nx,Ny,Nz)/I/FREE QspaceNorm=0
//	SetScale/P x v.xlo,dQ,"nm\S-1\M", Qspace, QspaceNorm
//	SetScale/P y v.ylo,dQ,"nm\S-1\M", Qspace, QspaceNorm
//	SetScale/P z v.zlo,dQ,"nm\S-1\M", Qspace, QspaceNorm
//	v.xhi = (v.xlo) + (Nx-1)*dQ
//	v.yhi = (v.ylo) + (Ny-1)*dQ
//	v.zhi = (v.zlo) + (Nz-1)*dQ
//	updateBoundingVolumeStruct(v)
//	print "Histogram Q-Volume is: ",boundingVolumeStruct2str(v)
//	printf "dQ = %g 1/nm,   N = [%g, %g, %g],   %g points\r",dQ,Nx,Ny,Nz,Nqtot
//	Variable unWrapped=0							// flags that unwrapping occured
//
//	ProgressPanelUpdate(progressWin,0,status="part 2 of 2",resetClock=1)
//	Variable i
//Variable findQsec=0, findQtimer, histQsec=0, histQtimer
//	for (i=0, scanNum=str2num(scanRange); numtype(scanNum)==0 && stop==0; scanNum=NextInRange(scanRange,scanNum))
//		for (m=0;m<Nm;m+=1)
//			A[m] = specInfo(scanNum,StringFromList(m,motors))
//		endfor
//		keV = specInfo(scanNum,"DCM_energy")
//		keV = numtype(keV) ? 10*hc_keVnm/specInfo(scanNum,"wavelength") : keV
//		Np = specInfo(scanNum,"Npoints")			// number of points in this scan
//		for (ip=0;ip<Np;ip+=1)					// for each point in this scan
//			if (mod(i,5)==0)
//				if (ProgressPanelUpdate(progressWin,i/N*100))
//					stop = 1
//					break
//				endif
//			endif
//			imageName = LoadSpecImage(scanNum,ip,extras="quiet:1")
//			preExists=(strsearch(imageName,"preExists:",0)==0)
//			if (preExists)
//				Wave image = $StringByKey("preExists", imageName)
//			else
//				Wave image = $imageName
//			endif
//			if (!WaveExists(image))
//				sprintf str,"Could not load scan %g,  point %g",scanNum,ip
//				print str
//				DoAlert 0,str
//				return $""
//			endif
//			wnote = note(image)
//			for (m=0;m<Nm;m+=1)
//				mname = CleanupName(StringFromList(m,motors),0)
//				val = NumberByKey(mname,wnote,"=")
//				A[m] = numtype(val) ? A[m] : val
//			endfor
//			val = NumberByKey("keV",wnote,"=")	;	keV = numtype(val) ? keV : val
//			filterTrans = NumberByKey("filterTrans",wnote,"=")
//			filterTrans = filterTrans>0 ? filterTrans : 1
//			I0 = NumberByKey("I0",wnote,"=")
//			I0 = I0 > 0 ? I0 : 1
//			scale = AvgI0 / (filterTrans*I0)		// correction for filters and I0
//
//			#ifdef EnablePilatusUnWrapping		// to enable the unwrappng functions
//				unWrapped = UnWrapImage(image) || unWrapped
//			#endif
//
//			Redimension/N=(numpnts(image)) image
//findQtimer=startMSTimer
//			Wave qvecs =  diffractometer#QofPixelVEC(d,keV,pxpy,A)	// ALL q-vectors for pixels in image
//findQsec += stopMSTimer(findQtimer)
//histQtimer=startMSTimer
//			HistogramOneImage(qvecs,image,Qspace, QspaceNorm,scale,mask1D)
//histQsec += stopMSTimer(histQtimer)
//			WaveClear qvecs
//			if (!preExists)
//				KillWaves/Z image					// only kill those waves that I have just loaded
//			endif
//			i += 1
//		endfor
//	endfor
//	Variable seconds2=SecondsInProgressPanel(progressWin)
//	printf "  the 2nd part took  %s,   the whole process took  %s\r",Secs2Time(seconds2,5,0),Secs2Time(seconds1+seconds2,5,0)
//	DoWindow/K $progressWin
//printf "finding qvecs took %g sec,   histograming took %g sec\r",findQsec*1e-6,histQsec*1e-6
//
//	MatrixOP/FREE eq0 = equal(QspaceNorm,0)		// set to 1 if point not used
//	Variable zeros=sum(eq0)						// number of zero points in Qspace
//
//	Duplicate/FREE QspaceNorm,QspaceNorm1
//	QspaceNorm1 = max(QspaceNorm,1)			// avoid divide by zero
//	Qspace /= QspaceNorm1							// normalize by number of times each point sampled
//	wnote = ReplaceStringByKey("waveClass",wnote,"Qspace3D","=")
//	for (m=0;m<Nm;m+=1)
//		wnote = RemoveByKey(StringFromList(m,motors),wnote,"=")
//	endfor
//	wnote = ReplaceNumberByKey("zeros", wnote,zeros,"=")
//	Wave qvec = BraggPeak(Qspace)
//	wnote = ReplaceStringByKey("qPeak", wnote,vec2str(qvec,bare=1),"=")
//
//	String command = specInfoT(scanNum,"specCommand")
//	String file_time = specInfoT(scanNum,"timeWritten")
//	String title = specInfoT(scanNum,"comment")
//	wnote = ReplaceStringByKey("scanRange", wnote,scanRange,"=")
//	if (strlen(command))
//		wnote = ReplaceStringByKey("command", wnote,command,"=")
//	endif
//	if (strlen(file_time))
//		wnote = ReplaceStringByKey("file_time", wnote,file_time,"=")
//	endif
//	if (strlen(title))
//		wnote = ReplaceStringByKey("title", wnote,title,"=")
//	endif
//	Note/K Qspace, wnote
//	Note/K QspaceNorm, ReplaceStringByKey("waveClass", wnote,"Qspace3DNorm","=")
//
//	WaveStats/Q/M=1 Qspace
//	printf "V_npnts = %g, V_numNaNs = %g, V_numINFs = %g\r",V_npnts, V_numNaNs, V_numINFs
//	printf "V_avg = %g,  V_Sum = %g,  V_sdev = %g\r",V_avg, V_Sum, V_sdev
//	printf "there are %g zeros, and %g non-zero points  (total = %g)\r",zeros,Nqtot-zeros,Nqtot
//	printf "min[%g, %g, %g] = %g\r",V_minRowLoc,V_minColLoc,V_minLayerLoc,V_min
//	printf "max[%g, %g, %g] = %g\r",V_maxRowLoc,V_maxColLoc,V_maxLayerLoc,V_max
//	#ifdef EnablePilatusUnWrapping
//		printf "%s\r", SelectString(unWrapped,"NO Unwrapping needed","Images were UNwrapped")
//	#endif
//	MakeSampledPoints(QspaceNorm,name="SampledPoints"+QspaceNameEnd+"XYZ")
//	MakeGizmocubeCorners(Qspace)
//	Qspace = Qspace==0 ? NaN : Qspace
//
//	return Qspace
//End
////
//Static Function HistogramOneImage(qvecs,image1D,Qspace, QspaceNorm,scale,mask1D)
//	Wave qvecs									// q vector for each pixel in detector
//	Wave image1D
//	Wave Qspace, QspaceNorm
//	Variable scale
//	Wave mask1D
//	scale = numtype(scale) || scale==0 ? 1 : scale
//
//	Variable N=numpnts(image1D), i
//	Make/N=3/D/FREE Qxyz0=DimOffset(Qspace,p), idxyz=(1/DimDelta(Qspace,p))
//	MatrixOP/FREE mxyz = (qvecs - rowRepeat(Qxyz0,N)) * rowRepeat(idxyz,N)
//	if (WaveExists(mask1D))
//		for (i=0;i<N;i+=1)
//			if (mask1D[i])
//				Qspace[mxyz[i][0]][mxyz[i][1]][mxyz[i][2]] += image1D[i] * scale
//				QspaceNorm[mxyz[i][0]][mxyz[i][1]][mxyz[i][2]] += 1
//			endif
//		endfor
//	else
//		for (i=0;i<N;i+=1)
//			Qspace[mxyz[i][0]][mxyz[i][1]][mxyz[i][2]] += image1D[i] * scale
//			QspaceNorm[mxyz[i][0]][mxyz[i][1]][mxyz[i][2]] += 1
//		endfor
//	endif
//End

#else			// MULTI THREADED

Function/WAVE MakeQspaceVolumeSpec(scanRange,[mask,Nthreads,doConvex])
	String scanRange						// combine all images from every point in each scan from scanRange
	Wave mask								// used pixels that are 1, ignore pixels that are 0
	Variable Nthreads						// number of threads to use
	Variable doConvex						// if True, make the convex hull
	doConvex = ParamIsDefault(doConvex) ? 0 : doConvex
	doConvex = numtype(doConvex) ? NaN : !(!doConvex)

	String motors = StrVarOrDefault("root:Packages:Diffractometer:DiffractometerAxisNames","")
	Variable Nm=ItemsInList(motors)
	if (Nm<1)
		DoAlert 0, "Unable to find list of motor names, set Diffractometer Type"
		return $""
	endif

	Variable printIt=0
	if (ItemsInRange(scanRange)<1 || numtype(doConvex))	// bad or empty range, prompt user
		if (ItemsInRange(scanRange)<1)
			Variable last = NumVarOrDefault("root:Packages:spec:lastScan",NaN)
			scanRange = SelectString(numtype(last),num2str(last),"")
		endif
		Prompt scanRange,"range of spec scans to process together, e.g. \"1-20\""
		Prompt doConvex,"Also compute Convex Hull (can take a while)",popup,"No Convex Hull;Compute Convex Hull"
		doConvex = numtype(doConvex) ? 1 : doConvex
		doConvex += 1
		DoPrompt "scans",scanRange, doConvex
		if (V_flag)
			return $""
		endif
		doConvex -= 1
		printIt = 1
	endif
	if (ItemsInRange(scanRange)<1)				// number of scans
		return $""
	endif
	Variable i
	String str, maskName,maskList=""
	if (!WaveExists(mask))						// if mask not passed, try to make maskList
		String imageFldr
		maskList = WaveListClass("imageMask","*","DIMS:2")	// check for masks in current folder
		for (i=0;i<ItemsInRange(scanRange);i+=1)				// and also check for masks in images folder for each scan
			sprintf imageFldr,"%sspec%d:images:",StrVarOrDefault("root:Packages:spec:specDataFolder",""),PointInRange(i,scanRange)
			maskList += WaveListClass("imageMask","*","DIMS:2",fldr=imageFldr)
		endfor
	endif
	if (strlen(maskList))							// need to ask user if they want to use a mask
		Prompt maskName,"mask to use with all images",popup, "_none_;"+maskList
		DoPrompt "Mask",maskName
		if (V_flag)
			return $""
		endif
		maskName = SelectString(stringmatch(maskName,"_none_"),maskName,"")
		printIt = printIt || (strlen(maskName)>0)
		Wave mask = $maskName
	elseif(WaveExists(mask))
		maskName = GetWavesDataFolder(mask,3)
	else
		maskName = ""
	endif
	if (printIt)
		printf "•MakeQspaceVolumeSpec(\"%s\"", scanRange
		if (strlen(maskName))
			printf ", mask=%s",maskName
		endif
		if (!ParamIsDefault(doConvex) || doConvex)
			printf ", doConvex=%g",doConvex
		endif
		if (!ParamIsDefault(Nthreads))
			printf ", Nthreads=%g",Nthreads
		endif
		printf ")\r"
	endif
	if (ParamIsDefault(Nthreads) || numtype(Nthreads) || Nthreads<=0)
		Nthreads = ThreadProcessorCount
		Nthreads -= (Nthreads>=4) ? 2 : 0			// when Nthreads is 4 or greater, lower by two
	endif
	printf "Using MultiThreading with %g threads\r",Nthreads

	STRUCT boundingVolume v
	initBoundingVolumeStruct(v)
	Make/N=(Nm)/D/FREE A
	Make/N=3/D/FREE q0,q1
	String mname, wnote, imageName
	Variable m, keV, val
	Variable scanNum,N								// N counts number of images to process
	Variable filterTrans,I0, AvgI0=0,NAvgI0=0, scale
	Variable ip,Np									// for looping over points in one scan

	String fileName=StrVarOrDefault("root:Packages:spec:specDefaultFile","")
	String path=StrVarOrDefault("root:Packages:spec:specDefaultPath","raw")
	LoadRangeOfSpecScans(scanRange,fileName,path)// ensure that the spec scans have already been loaded
	Variable NimagesExpected=0
	for (scanNum=str2num(scanRange); numtype(scanNum)==0; scanNum=NextInRange(scanRange,scanNum))
		NimagesExpected += specInfo(scanNum,"Npoints")
	endfor
	if (!(NimagesExpected>0))
		str = "Not All images in range '"+scanRange+"' were loaded"
		print str
		DoAlert 0, str
		return $""
	endif

	wnote = LoadSpecImageInfo(str2num(scanRange),0)
	STRUCT detectorGeometry d
	if (FillDetectorStructDefault(d,StringByKey("detectorID",wnote,"=")))	//fill the detector structure with current values
		DoAlert 0, "No Detector Structure, please set one"
		return $""
	endif

	Variable Nalloc=NimagesExpected
	Make/N=(Nalloc,Nm)/D/FREE A_W=NaN
	Make/N=(Nalloc)/D/FREE keV_W=NaN, scanNum_W=NaN, ip_W=NaN

	String sPart = SelectString(doConvex,"2","3")
	String progressWin = ProgressPanelStart("",stop=1,showTime=1,status="part 1 of "+sPart)
	Variable stop=0, preExists=0, Npixels=0
	for (N=0, scanNum=str2num(scanRange); numtype(scanNum)==0 && stop==0; scanNum=NextInRange(scanRange,scanNum))
		for (m=0;m<Nm;m+=1)
			A[m] = specInfo(scanNum,StringFromList(m,motors))
		endfor
		keV = specInfo(scanNum,"DCM_energy")
		keV = numtype(keV) ? 10*hc_keVnm/specInfo(scanNum,"wavelength") : keV
		Np = specInfo(scanNum,"Npoints")			// number of points in this scan
		for (ip=0;ip<Np;ip+=1)					// for each point in this scan
			if (mod(N,50)==0)
				if (ProgressPanelUpdate(progressWin,N/NimagesExpected*100))
					stop = 1
					break
				endif
			endif

			wnote = LoadSpecImageInfo(scanNum,ip)
			Npixels = max(Npixels,NumberByKey("xdim",wnote,"=")*NumberByKey("ydim",wnote,"="))
			if (!(Npixels>0))
				sprintf str,"Could not load scan %g,  point %g",scanNum,ip
				print str
				DoAlert 0,str
				return $""
			endif

			for (m=0;m<Nm;m+=1)
				mname = CleanupName(StringFromList(m,motors),0)
				val = NumberByKey(mname,wnote,"=")
				A[m] = numtype(val) ? A[m] : val
			endfor
			val = NumberByKey("keV",wnote,"=")	;	keV = numtype(val) ? keV : val
			diffractometer#QofPixel(d,keV,d.Nx/2,d.Ny/2,A,q1)	// q0 & q1 only used to find dq, not volume
			if (N==0)
				q0 = q1
			endif
			extendQvolumeFromImage(d,keV,A,v)	// uses the four corners to extend volume in Q-space

			I0 = NumberByKey("I0",wnote,"=")
			if (I0>0)
				AvgI0 += I0	
				NAvgI0 += 1
			endif

			if (N>=Nalloc)
				Nalloc += 500						// need to extend arrays
				Redimension/N=(Nalloc) keV_W, scanNum_W, ip_W
				Redimension/N=(Nalloc,-1) A_W
			endif

			// save these values to make next part easier to multithread
			scanNum_W[N] = scanNum				// save scan number
			ip_W[N] = ip							// save point in scan
			A_W[N][] = A[q]						// save the motor positions
			keV_W[N] = keV						// save the energy
			N+=1
		endfor
	endfor
	Variable seconds1 = SecondsInProgressPanel(progressWin)
	if (stop)
		DoWindow/K $progressWin
		print "Stopped in first part by the user after",Secs2Time(seconds1,5,0)
		return $""
	endif
	printf "  the first part took  %s\r",Secs2Time(seconds1,5,0)

	AvgI0 = (numtype(AvgI0)==0 || AvgI0>0 || NAvgI0>0) ? AvgI0/NAvgI0 : 1
	Redimension/N=(N) keV_W, scanNum_W, ip_W
	Redimension/N=(N,-1) A_W
	printf "Found %g images to process,  with <I0> = %g\r",N,AvgI0
	if (Npixels<1 || N<1)
		print "found number of pixels=%g,   and   number of images = %g,   nothing to do.\r",Npixels,N
		return $""
	endif
	print "Measured Q-Volume is: ",boundingVolumeStruct2str(v)

	if (Npixels==numpnts(mask))					// make a 1D version of mask (same number of pixels as image)
		Duplicate/FREE mask, mask1D
		Redimension/N=(Npixels) mask1D
	else
		Wave mask1D = $""							// no mask
	endif

	String QspaceNameEnd=num2str(str2num(scanRange))
	QspaceNameEnd += SelectString(ItemsInRange(scanRange)>1,"","_"+num2str(lastInRange(scanRange)))
	String QspaceName="Qspace"+QspaceNameEnd
	Make/N=3/D/FREE q01 = q1-q0				// find dQ
	Variable dQ = norm(q01)/(N-1)
print "dQ start = ",dQ
	dQ = max(dQ,max(max(v.dx,v.dy),v.dz)/MAX_Q_SPACE_DIM)
print "dQ after = ",dQ
	Variable Nx=ceil(v.dx / dQ)+1, Ny=ceil(v.dy / dQ)+1, Nz=ceil(v.dz / dQ)+1
	Variable Nqtot = Nx*Ny*Nz

	Variable bad = 0								// check for reasonable sized histogram
	bad += numtype(Nx+Ny+Nz+ v.xlo +v.ylo + v.zlo + dQ)!=0
	bad += (dQ<=0)
	bad += (Nx<1 || Ny<1 || Nz<1)					// not too small
	bad += (Nx>2000 || Ny>2000 || Nz>2000)		//  & not too big
	bad += (Nqtot >5e7)							//  & not too big
	if (bad)
		sprintf str,"Unreasonable sized histogram (%g  x  %g  x  %g)\r %s",Nx,Ny,Nz,boundingVolumeStruct2str(v)
		DoAlert 0, str
		print str
		return $""
	endif

	Make/N=(Nx,Ny,Nz)/D/O $QspaceName=0
	Wave Qspace=$QspaceName
	Make/N=(Nx,Ny,Nz)/I/FREE QspaceNorm=0
	SetScale/P x v.xlo,dQ,"nm\S-1\M", Qspace, QspaceNorm
	SetScale/P y v.ylo,dQ,"nm\S-1\M", Qspace, QspaceNorm
	SetScale/P z v.zlo,dQ,"nm\S-1\M", Qspace, QspaceNorm
	v.xhi = (v.xlo) + (Nx-1)*dQ
	v.yhi = (v.ylo) + (Ny-1)*dQ
	v.zhi = (v.zlo) + (Nz-1)*dQ
	updateBoundingVolumeStruct(v)
	print "Histogram Q-Volume is: ",boundingVolumeStruct2str(v)
	printf "dQ = %g 1/nm,   N = [%g, %g, %g],   %g points\r",dQ,Nx,Ny,Nz,Nqtot
	Variable unWrapped=0							// flags that unwrapping occured

	ProgressPanelUpdate(progressWin,0,status="part 2 of "+sPart,resetClock=1)
	String strStruct
	StructPut/S/B=2 d, strStruct
	Variable iThreads
	Variable/G threadGroupID = ThreadGroupCreate(Nthreads)
	for(i=0;i<Nthreads;i+=1)
		ThreadStart threadGroupID,i,HistogramOneImageMT()
	endfor

	Variable iDone, dfrStatus, j, imod=max(Nthreads,5), iNext=-1
	for (iThreads=0,iDone=0, j=0; j<N; j+=1)		// This loop actually does the histogramming
		if (iDone>iNext)
			iNext = iDone + imod
			if (ProgressPanelUpdate(progressWin,iDone/N*100))
				break
			endif
		endif

		imageName = LoadSpecImage(scanNum_W[j],ip_W[j],extras="quiet:1")
		preExists=(strsearch(imageName,"preExists:",0)==0)
		if (preExists)
			Wave image = $StringByKey("preExists", imageName)
		else
			Wave image = $imageName
		endif
		if (!WaveExists(image))
			sprintf str,"Could not load scan %g,  point %g",scanNum_W[j],ip_W[j]
			DoWindow/K $progressWin
			print str
			DoAlert 0,str
			return $""
		endif
		wnote = note(image)
		filterTrans = NumberByKey("filterTrans",wnote,"=")
		filterTrans = filterTrans>0 ? filterTrans : 1
		I0 = NumberByKey("I0",wnote,"=")
		I0 = I0 > 0 ? I0 : 1
		scale = AvgI0 / (filterTrans*I0)			// correction for filters and I0

		#ifdef EnablePilatusUnWrapping			// to enable the unwrappng functions
			unWrapped = UnWrapImage(image) || unWrapped
		#endif

		NewDataFolder/S forThread					// set up data folder for calculation
		Variable/G :keV=keV_W[j]
		Variable/G :Nx=Nx, :Ny=Ny, :Nz=Nz
		Variable/G :xlo = v.xlo, :ylo = v.ylo, :zlo = v.zlo
		Variable/G :dQ=dQ, :scale=scale
		String/G :strStruct=strStruct

		Wave pxpy_MT = diffractometer#MakePixelListFromImage(image,pxpyName="pxpy_MT")	// re-make pixel list (full-chip, unbinned pixel units)
		Duplicate A, A_MT
		A_MT = A_W[j][p]
		Duplicate image, image1D
		if (!preExists)
			KillWaves/Z image						// only kill those images that I have just loaded, not pre-existing images
		endif
		Redimension/N=(numpnts(image1D)) image1D
		if (WaveExists(mask1D))
			Duplicate mask1D, mask1D_MT
		endif
		WaveClear pxpy_MT, A_MT, image1D, mask1D_MT
		ThreadGroupPutDF threadGroupID,:			// Send current data folder to input queue
		iThreads += 1								// number of threads currently executing

		DFREF dfr= ThreadGroupGetDFR(threadGroupID,0)	// Get results from a finished thread
		do											// process any finished threads
			dfrStatus = AccumulateOneHistogramThread(dfr,Qspace,QspaceNorm)
			if (dfrStatus)							// found finished thread, process it
				iThreads -= 1						// number of threads currently executing
				iDone += 1
			endif

			if (dfrStatus)
				DFREF dfr= ThreadGroupGetDFR(threadGroupID,0)	// find a finished thread
				dfrStatus = DatafolderRefStatus(dfr)
			elseif (iThreads>=(Nthreads+1))
				sleep/C=-1/T 1
				DFREF dfr= ThreadGroupGetDFR(threadGroupID,0)	// find a finished thread
				dfrStatus = 1
			endif
		while(dfrStatus)
	endfor

	do												// process all remaining threads
		DFREF dfr= ThreadGroupGetDFR(threadGroupID,0)	// Get results from a finished thread
		if (AccumulateOneHistogramThread(dfr,Qspace,QspaceNorm))
			iThreads -= 1							// number of threads currently executing
		elseif (iThreads>0)							// still some threads left to wait for
			sleep/C=-1/T 1
		endif
	while(iThreads>0)

	Variable seconds2=SecondsInProgressPanel(progressWin), seconds3=0
	printf "  the 2nd part took  %s\r",Secs2Time(seconds2,5,0)

	if(ThreadGroupRelease(threadGroupID) == -2)	// Terminate the HistogramOneImageMT by setting an abort flag
		str = "Thread would not quit normally, had to force kill it. Restart Igor."
		DoAlert 0, str
		Print str
	else
	KillVariables threadGroupID
	endif

	MatrixOP/FREE eq0 = equal(QspaceNorm,0)		// set to 1 if point not used
	Variable zeros=sum(eq0)						// number of zero points in Qspace

	Duplicate/FREE QspaceNorm,QspaceNorm1
	QspaceNorm1 = max(QspaceNorm,1)			// avoid divide by zero
	Qspace /= QspaceNorm1							// normalize by number of times each point sampled
	wnote = ReplaceStringByKey("waveClass",wnote,"Qspace3D","=")
	for (m=0;m<Nm;m+=1)
		wnote = RemoveByKey(StringFromList(m,motors),wnote,"=")
	endfor
	wnote = ReplaceNumberByKey("zeros", wnote,zeros,"=")
	Wave qvec = BraggPeak(Qspace)
	wnote = ReplaceStringByKey("qPeak", wnote,vec2str(qvec,bare=1),"=")
	wnote = ReplaceNumberByKey("Nthreads", wnote,Nthreads,"=")
	wnote = ReplaceNumberByKey("executionTime", wnote,seconds1+seconds2,"=")

	scanNum = str2num(scanRange)				// first scan number
	String command = specInfoT(scanNum,"specCommand")
	String file_time = specInfoT(scanNum,"timeWritten")
	String title = specInfoT(scanNum,"comment")
	wnote = ReplaceStringByKey("scanRange", wnote,scanRange,"=")
	if (strlen(command))
		wnote = ReplaceStringByKey("command", wnote,command,"=")
	endif
	if (strlen(file_time))
		str = specTime2Igor(file_time,2)			// try to convert ugly spec time to nicer time string
		str = SelectString(strlen(str),file_time,str)
		wnote = ReplaceStringByKey("file_time", wnote,str,"=")
	endif
	if (strlen(title))
		wnote = ReplaceStringByKey("title", wnote,title,"=")
	endif
	Note/K Qspace, wnote
	Note/K QspaceNorm, ReplaceStringByKey("waveClass", wnote,"Qspace3DNorm","=")

	WaveStats/Q/M=1 Qspace
	printf "V_npnts = %g, V_numNaNs = %g, V_numINFs = %g\r",V_npnts, V_numNaNs, V_numINFs
	printf "V_avg = %g,  V_Sum = %g,  V_sdev = %g\r",V_avg, V_Sum, V_sdev
	printf "there are %g zeros, and %g non-zero points  (total = %g)\r",zeros,Nqtot-zeros,Nqtot
	printf "min[%g, %g, %g] = %g\r",V_minRowLoc,V_minColLoc,V_minLayerLoc,V_min
	printf "max[%g, %g, %g] = %g\r",V_maxRowLoc,V_maxColLoc,V_maxLayerLoc,V_max
	#ifdef EnablePilatusUnWrapping
		printf "%s\r", SelectString(unWrapped,"NO Unwrapping needed","Images were UNwrapped")
	#endif
	Wave sampled=MakeSampledPoints(QspaceNorm,name="SampledPoints"+QspaceNameEnd+"XYZ")
	wnote = note(sampled)
	if (strlen(StringByKey("sourceWaveFullPath",wnote,"="))<1)
		wnote = ReplaceStringByKey("sourceWave",wnote,NameOfWave(Qspace),"=")
		wnote = ReplaceStringByKey("sourceWaveFullPath",wnote,GetWavesDataFolder(Qspace,2),"=")
		Note/K sampled, wnote
	endif

	MakeGizmocubeCorners(Qspace)
	Qspace = Qspace==0 ? NaN : Qspace

	if (doConvex)
		ProgressPanelUpdate(progressWin,0,status="part 3 of 3, computing Convex Bounding Volume",resetClock=1)
		QspaceVolumesView#MakeConvexHullFrom3D(Qspace,name="ConvexHullTrianlges"+QspaceNameEnd)
		seconds3 = SecondsInProgressPanel(progressWin)
		printf "  the 3rd part took  %s,   the whole process took  %s\r",Secs2Time(seconds3,5,0),Secs2Time(seconds1+seconds2+seconds3,5,0)
	endif
	DoWindow/K $progressWin
	printf "  the whole process took  %s\r",Secs2Time(seconds1+seconds2+seconds3,5,0)

	if (strlen(GetRTStackInfo(2))<1)
		MakeGizmoQspace3D(Qspace)
	endif
	return Qspace
End
//
Static Function AccumulateOneHistogramThread(dfr,Qspace,QspaceNorm)	// returns dfrStatus (0== non-existant)
	DFREF dfr										// a data folder reference, holds result of one histogram
	Wave Qspace
	Wave QspaceNorm

	Variable dfrStatus=DatafolderRefStatus(dfr)
	if (dfrStatus==0)								// an invalid folder reference
		return 0
	endif

	Wave QspaceMT=dfr:Qspace,  QspaceNormMT=dfr:QspaceNorm
	if (WaveExists(QspaceMT) && WaveExists(QspaceNormMT))
		MatrixOP/FREE QMT = Qspace + QspaceMT	// accumulate in final arrays
		MatrixOP/FREE QnormMT = QspaceNorm + QspaceNormMT
		// Qspace += QspaceMT
		// QspaceNorm += QspaceNormMT
	endif
	WaveClear QspaceMT, QspaceNormMT
	KillDataFolder dfr								// Redundant because dfr refers to a free data folder

	if (WaveExists(QMT) && WaveExists(QnormMT))
		MatrixOP/O/S Qspace = QMT				// transfer to in final
		MatrixOP/O/S QspaceNorm = QnormMT
		WaveClear QMT, QnormMT
	endif
	return dfrStatus
End
//
//Static Function AccumulateOneHistogramThread(dfr,Qspace,QspaceNorm)	// returns dfrStatus (0== non-existant)
//	DFREF dfr										// a data folder reference, holds result of one histogram
//	Wave Qspace
//	Wave QspaceNorm
//
//	Variable dfrStatus=DatafolderRefStatus(dfr)
//	if (dfrStatus==0)								// an invalid folder reference
//		return 0
//	endif
//
//	Wave QspaceMT=dfr:Qspace
//	Wave QspaceNormMT=dfr:QspaceNorm
//	if (WaveExists(QspaceMT) && WaveExists(QspaceNormMT))
//		Qspace += QspaceMT						// accumulate in final arrays
//		QspaceNorm += QspaceNormMT
//	endif
//	WaveClear QspaceMT, QspaceNormMT
//	KillDataFolder dfr								// Redundant because dfr refers to a free data folder
//	return dfrStatus
//End
////
ThreadSafe Static Function HistogramOneImageMT()
	do
		do	// wait for a free thread
			DFREF dfr = ThreadGroupGetDFR(0,inf)	// Get free data folder from input queue
			if (DataFolderRefStatus(dfr) == 0)
				if( GetRTError(2) )				// New in 6.2 to allow this distinction:
					Print "worker closing down due to group release"
				endif
			else
				break
			endif
		while(1)

		WAVE pxpy=dfr:pxpy_MT, A=dfr:A_MT
		SVAR strStruct=dfr:strStruct
		NVAR keV=dfr:keV
		STRUCT detectorGeometry d					// one detector geometry for whole image
		StructGet/S/B=2 d, strStruct		
		NewDataFolder/S outDF						// set to new data folder
		Wave qvecs = diffractometer#QofPixelVEC(d,keV,pxpy,A)	// ALL q-vectors for pixels in image
		WaveClear pxpy, A

		NVAR Nx=dfr:Nx, Ny=dfr:Ny, Nz=dfr:Nz
		NVAR xlo=dfr:xlo, ylo=dfr:ylo, zlo=dfr:zlo, dQ=dfr:dQ
		Variable bad = !WaveExists(qvecs)
		bad += numtype(Nx+Ny+Nz+xlo+ylo+zlo+dQ)!=0
		bad += (dQ<=0)
		bad += (Nx<1 || Ny<1 || Nz<1)				// not too small
		bad += (Nx>2000 || Ny>2000 || Nz>2000)	//  & not too big
		bad += (Nx*Ny*Nz > 5e7)					//  & not too big
		if (bad)										// ERROR, skip this one
			printf "found bad=%g,   WaveExists(qvecs) =%g\r",bad,WaveExists(qvecs)
			printf "Nx=%g,  Ny=%g,  Nz=%g,  dQ = %g\r",Nx,Ny,Nz,dQ
			printf "xlo = %g, ylo = %g, zlo = %g\r",xlo,ylo,zlo

		else
			WAVE image1D=dfr:image1D			// 1D, unwrapped, version of image to process
			WAVE mask1D=dfr:mask1D_MT			// this may be missing
			NVAR scale=dfr:scale
			scale = numtype(scale) || scale==0 ? 1 : scale

			Make/N=(Nx,Ny,Nz)/D/O Qspace=0		// these two get made in outDF
			Make/N=(Nx,Ny,Nz)/I QspaceNorm=0	// note, Qspace & QspaceNorm get combined with full versions
			SetScale/P x xlo,dQ,"nm\S-1\M", Qspace, QspaceNorm
			SetScale/P y ylo,dQ,"nm\S-1\M", Qspace, QspaceNorm
			SetScale/P z zlo,dQ,"nm\S-1\M", Qspace, QspaceNorm

			Variable N=numpnts(image1D), i
			Make/N=3/D/FREE Qxyz0=DimOffset(Qspace,p), idxyz=(1/DimDelta(Qspace,p))
			MatrixOP/FREE mxyz = (qvecs - rowRepeat(Qxyz0,N)) * rowRepeat(idxyz,N)
			if (WaveExists(mask1D))
				for (i=0;i<N;i+=1)
					if (mask1D[i])
						Qspace[mxyz[i][0]][mxyz[i][1]][mxyz[i][2]] += image1D[i] * scale
						QspaceNorm[mxyz[i][0]][mxyz[i][1]][mxyz[i][2]] += 1
					endif
				endfor
			else
				for (i=0;i<N;i+=1)
					Qspace[mxyz[i][0]][mxyz[i][1]][mxyz[i][2]] += image1D[i] * scale
					QspaceNorm[mxyz[i][0]][mxyz[i][1]][mxyz[i][2]] += 1
				endfor
			endif
		endif
		WaveClear image1D, mask1D
		WaveClear Qspace, QspaceNorm
		WaveClear qvecs, Qxyz0, idxyz, mxyz
		ThreadGroupPutDF 0,:						// Put current data folder in output queue
		KillDataFolder dfr							// We are done with the input data folder
	while(1)

	return 0
End
#endif

Static Function extendQvolumeFromImage(d,keV,A,v)
	STRUCT detectorGeometry &d
	Variable keV
	Wave A
	STRUCT boundingVolume &v

	Variable Nx = d.Nx, Ny = d.Ny
	Make/N=3/D/FREE qvec

	diffractometer#QofPixel(d,keV,0,0,A,qvec)			// (0,0) corner of detector
	extendBoundingVolumeStruct(v,qvec)

	diffractometer#QofPixel(d,keV,Nx-1,0,A,qvec)		// (Xmax,0) corner of detector
	extendBoundingVolumeStruct(v,qvec)

	diffractometer#QofPixel(d,keV,0,Ny-1,A,qvec)		// (0,Ymax) corner of detector
	extendBoundingVolumeStruct(v,qvec)

	diffractometer#QofPixel(d,keV,Nx-1,Ny-1,A,qvec)	// (Xmax,Ymax) corner of detector
	extendBoundingVolumeStruct(v,qvec)
End


//Function/S Rescale3DbyQ4(Qspace)
//	Wave Qspace
//
//	Variable qx,qy,qz
//	sscanf StringByKey("qpeak", note(Qspace),"="), "%g, %g, %g", qx,qy,qz
//	if (V_flag!=3)
//		return ""
//	endif
//
//	Variable dx=DimDelta(Qspace,0), dy=DimDelta(Qspace,1), dz=DimDelta(Qspace,2)
//	Variable x0=DimOffset(Qspace,0), y0=DimOffset(Qspace,1), z0=DimOffset(Qspace,2)
//	Duplicate/FREE Qspace, Q4
//	Q4 = Qspace[p][q][r] * ( (p*dx+x0-qx)^2 + (q*dy+y0-qy)^2 + (r*dz+z0-qz)^2) ^2
//	String wName = NameOfWave(Qspace)+"Q4"
//	Duplicate/O Q4, $wName
//	Note/K $wName,ReplaceNumberByKey("Qpower",note(Qspace),4,"=")
//	printf "rescaled '%s' by q^4, using a Bragg peak at {%g, %g, %g}(1/nm), result in '%s'\r",NameOfWave(Qspace),qx,qy,qz,wName
//	return wName
//End


Static Function/WAVE MakeSampledPoints(Qspace,[name])
	Wave Qspace
	String name
	name = SelectString(ParamIsDefault(name),name,"SampledPointsXYZ")
	name = SelectString(strlen(name),"SampledPointsXYZ",name)

	Variable Nx=DimSize(Qspace,0), Ny=DimSize(Qspace,1), Nz=DimSize(Qspace,2)
	Variable N=Nx*Ny*Nz

	Variable dx=DimDelta(Qspace,0), dy=DimDelta(Qspace,1), dz=DimDelta(Qspace,2)
	Variable x0=DimOffset(Qspace,0), y0=DimOffset(Qspace,1), z0=DimOffset(Qspace,2)
	Variable xx,yy,zz

	Make/N=(N,3)/O $name=NaN
	Wave sampledXYZ = $name

	Variable i,j,k, m
	for (k=0,m=0;k<Nz;k+=1)
		zz = k*dz + z0
		for (j=0;j<Ny;j+=1)
			yy = j*dy + y0
			for (i=0;i<Nx;i+=1)
				xx = i*dx + x0
				if (Qspace[i][j][k]>0)
					sampledXYZ[m][0] = xx
					sampledXYZ[m][1] = yy
					sampledXYZ[m][2] = zz
					m += 1
				endif
			endfor
		endfor
	endfor
	Redimension/N=(m,-1) sampledXYZ

	String wnote="waveClass=sampledPoints;"
	wnote = ReplaceStringByKey("sourceWave",wnote,NameOfWave(Qspace),"=")
	wnote = ReplaceStringByKey("sourceWaveFullPath",wnote,GetWavesDataFolder(Qspace,2),"=")
	Note/K sampledXYZ, wnote
	return sampledXYZ
End

//  ============================= End of Qspace Histogram =============================  //
//  ============================================================================  //



//  ============================================================================  //
//  ============================ Start of spec Image Plotting ============================  //

Function/T rePlotDiffractometerImage(scanNum,pindex,[extras])
	Variable scanNum											// scanNumber in spec file
	Variable pindex												// point in scanNum
	String extras												// example		"quiet:1;multiple:0;"
	extras = SelectString(ParamIsDefault(extras),extras,"")

	if (!(scanNum>0) || !(pindex>=0) || numtype(scanNum+pindex))
		scanNum = scanNum>0 ? scanNum : NumVarOrDefault("root:Packages:spec:lastScan", 1)
		pindex = pindex>=0 ? pindex : 0
		Prompt scanNum,"spec Scan Number"
		Prompt pindex,"point in spec scan"
		DoPrompt "spec Scan & Point", scanNum,pindex
		if (V_flag)
			return ""
		endif
		if (!(scanNum>0) || numtype(scanNum))
			printf "scan number %g is invalid\r",scanNum
		endif
		if (!(pindex>=0) || numtype(pindex))
			printf "point number %g is invalid\r",pindex
		endif
		printf "rePlotDiffractometerImage(%g,%g",scanNum,pindex
		if (!ParamIsDefault(extras))
			printf ", extras=\"%s\"",extras
		endif
		printf ")\r"
	endif
	String 	specFldr
	sprintf specFldr,"%sspec%d:",StrVarOrDefault("root:Packages:spec:specDataFolder",""),scanNum
	if (!DataFolderExists(specFldr))
		printf "spec scan %g has not been loaded, you must load it first\r",scanNum
		return ""
	endif
	DFREF dfr=$specFldr										// Reference to folder with spec scan
	Wave ww=$(specFldr+GetIndexedObjNameDFR(dfr,1,0))
	if (!WaveExists(ww))
		return ""
	elseif (pindex<0 || pindex>=DimSize(ww,0) || numtype(pindex))
		printf "spec scan %g has been loaded, BUT it has only %d points, you asked for point=%g\r",scanNum,DimSize(ww,0),pindex
		return ""
	endif

	String imageFldr,wName									// name of wave to hold image
	String specPath = StrVarOrDefault("root:Packages:spec:specDefaultFile","")
	String name=ParseFilePath(3, specPath,":", 0,0)
	sprintf wName,"%s_S%d_%d",name,scanNum,pindex
	wName = CleanupName(wName,0)
	sprintf specFldr,"%sspec%d:",StrVarOrDefault("root:Packages:spec:specDataFolder",""),scanNum
	sprintf imageFldr,"%sspec%d:images:",StrVarOrDefault("root:Packages:spec:specDataFolder",""),scanNum
	Wave image = $(imageFldr+wName)
	if (WaveExists(image))
		GraphDiffractometerImage(image)
	else
		String str
		sprintf str, "In scan %g, image #%g is not loaded",scanNum,pindex
		printf str
		DoAlert 0, str
	endif
	return ""
End
//
Static Function/T GraphDiffractometerImage(image)
	Wave image
	if (!WaveExists(image))
		return ""
	endif
	String win=StringFromList(0,FindGraphsWithWave(image))
	if (strlen(win))												// image is already on graph, just bring graph to front
		DoWindow/F $win
		return win
	endif
	Variable faint=55000											// or try 4500

	String wnote=note(image), title, str
	Variable scanNum=NumberByKey("scanNum",wnote,"="), point=NumberByKey("scanPoint",wnote,"=")
	Variable keV=specInfo(scanNum,"DCM_energy")
	keV = numtype(keV) ? specInfo(scanNum,"keV") : keV
	Variable filters=NumberByKey("filters",wnote,"="), filterTrans=NumberByKey("filterTrans",wnote,"=")
	Variable Nx=DimSize(image,0), Ny=DimSize(image,1)
	Variable dMarker=round(max(Nx,Ny)/10)
	dMarker = min(dMarker,min(Nx,Ny))
	String comment = specInfoT(scanNum,"comment")
	String timeWritten = specInfoT(scanNum,"timeWritten")
	String command = specInfoT(scanNum,"specCommand")
	String fileName=StringByKey("file_name",wnote,"=")
	sprintf title,"scan%d - #%d,   %s",scanNum,point,comment
	if (strlen(timeWritten))
		title += "\r"+specTime2Igor(timeWritten,2)
	endif
	if (numtype(filters+filterTrans)==0 && filters>0)
		sprintf str, "filter[%s] = %.5g",specProc#filterNumber2Name(filters),filterTrans
		title += "\r"+str
	endif

	String DiffractometerAxisNames = StrVarOrDefault("root:Packages:Diffractometer:DiffractometerAxisNames","")
	Variable Naxes = ItemsInList(DiffractometerAxisNames)
	if (Naxes>0)													// add the diffractometer angles
		String aName, str1
		Variable i, angle
		str = ""
		for (i=0;i<Naxes;i+=1)
			aName = StringFromList(i,DiffractometerAxisNames)
			angle = NumberByKey(aName,wnote,"=")				// first check wavenote for angle
			angle = numtype(angle) ? specInfo(scanNum,aName) : angle	// angle is not in wavenote, so check specInfo()
			sprintf str1, "%s = %.3f°",aName,angle
			str += SelectString(i>0,"",", ")+str1					// accumulate the "name=angle, "
		endfor
		title += "\r\\Zr070"+str+"\\M"
	endif
	if (strlen(command))
		title += "\r\\Zr075"+command+"\\M"
	endif
	if (strlen(fileName))
		title += "\r\\Zr075"+fileName+"\\M"
	endif

	Display /W=(265,50,1001,350)								// build a new graph
	AppendImage image
	ModifyImage $NameOfWave(image) ctab= {*,*,Terrain256,1}
	ModifyGraph grid=1, tick=2, mirror=1, minor=1, lowTrip=0.001, axOffset(left)=-1.1,axOffset(bottom)=-1
	ModifyGraph gridRGB=(faint,faint,65535), gridStyle=1
	Label left "y  (\\U)"
	Label bottom "x  (\\U)"
	SetAxis/A/R left
	TextBox/C/N=textTitle/F=0/B=1/A=LB/X=2/Y=4 title
	Variable px=(Nx-1)/2, py=(Ny-1)/2
	STRUCT detectorGeometry d
	String id=StringByKey("detectorID",wnote,"=")
	if (!FillDetectorStructDefault(d,id))							//fill the detector structure with current values
		px = d.px0
		py = d.py0
	endif
	DrawMarker(px,py,dMarker,dMarker,"cross gap",dash=2,layer="UserAxes")
	SetWindow kwTopWin,hook(peakInfo)=specImage#getQofImageHook
	DoUpdate
	win = StringFromList(0,WinList("*",";","WIN:"))
	SetAspectToSquarePixels(win)
	return win
End
//
//		SetWindow kwTopWin,hook(peakInfo)=getQofImageHook
Static Function getQofImageHook(s)									// Shift-mouseDown=follow mouse and show Qvector
	STRUCT WMWinHookStruct &s
	String win = (s.winName)

	if (s.eventCode==5)												// this lets you leave the window up temporaily
		return 0
	elseif (s.eventMod!=3)
		Tag/K/N=indexedPeakInfo/W=$win
		s.doSetCursor= 0
		s.cursorCode = 0											// for nothing
		return 0
	endif

	String imageNameStr=StringFromList(0,ImageNameList(win,";"))
	Wave image = ImageNameToWaveRef(win,imageNameStr)
	if (!WaveExists(image))
		return 1
	endif
	String wnote=note(image)

	String id=StringByKey("detectorID",wnote,"=")
	STRUCT detectorGeometry d
	if (FillDetectorStructDefault(d,id))									//fill the geometry structure with current values
		return 1
	endif

	Variable Naxes=NumVarOrDefault("root:Packages:Diffractometer:Naxes",NaN)
	String DiffractometerAxisNames = StrVarOrDefault("root:Packages:Diffractometer:DiffractometerAxisNames","")
	if (numtype(Naxes)||strlen(DiffractometerAxisNames)<1 || Naxes<1)
		return 1
	endif
	Variable scanNum=NumberByKey("scanNum",wnote,"=")
	Variable keV=specInfo(scanNum,"DCM_energy")
	keV = numtype(keV) ? specInfo(scanNum,"keV") : keV
	keV = numtype(kev) ? NumberByKey("keV",wnote,"=") : keV
	Make/N=(Naxes)/D/FREE A=NumberByKey(StringFromLIst(p,DiffractometerAxisNames),wnote,"=")
	A = numtype(A[p]) ? specInfo(scanNum,StringFromLIst(p,DiffractometerAxisNames)) : A[p]	// angle is not in wavenote, so check specInfo()
	if (numtype(scanNum+keV+sum(A)))
		return 1
	endif

	GetWindow $win psize
	Variable vert = limit((s.mouseLoc.v-V_top)/(V_bottom-V_top),0,1)	// fractional position on graph
	Variable horiz = limit((s.mouseLoc.h-V_left)/(V_right-V_left),0,1)
	Variable theta, px,py											// pixel position at the mouse-click
	GetAxis/W=$win/Q bottom
	if (V_flag)
		return 1
	endif
	px = (V_max-V_min)*horiz + V_min
	GetAxis/W=$win/Q left
	if (V_flag)
		return 1
	endif
	py = (V_min-V_max)*vert + V_max

	Make/N=3/D/FREE qvec=NaN
	theta = diffractometer#QofPixel(d,keV,px,py,A,qvec)			// returns Q of pixel in SAMPLE (not xtal or beam line) frame

	STRUCT sampleStructure sa	
	String str = ParseFilePath(1,GetWavesDataFolder(image,1),":",1,0)+"sampleStructStr"
	String strStruct=StrVarOrDefault(str,"")						// fill the sample structure with values in spec scan directory
	if (strlen(strStruct)>0)										// need this to get hkl
		StructGet/S/B=2 sa, strStruct								// found structure information, load into s
		Wave hkl = HKLofPixel(sa,d,keV,px,py,A)					// returns hkl of pixel in beam CRYSTAL frame, (not sample frame, not BL coords)
	else
		Wave hkl=$""
	endif

	String tagStr=""
	sprintf str,"\f01q\f00\\BBL\\M = {%.3f, %.3f, %.3f} nm\\S-1\\M |q|= %.3f",qvec[0],qvec[1],qvec[2],norm(qvec)
	tagStr += str
	if (WaveExists(hkl))
		sprintf str,"\r\f01hkl\f00 = [%.3f, %.3f, %.3f]",hkl[0],hkl[1],hkl[2]
		tagStr += str
	endif
	sprintf str,"\r\\[0pixel [%.2f, %.2f],  \\F'Symbol'\\Zr1002q\\]0 = %.4f\\F'Symbol'∞\\]0",px,py,2*theta*180/PI
	tagStr += str

	px = limit(px,0,DimSize(image,0)-1)							// needed in case (px,py) is outside the image
	py = limit(py,0,DimSize(image,1)-1)
	Variable index = round(px) + DimSize(image,0)*round(py)		// convert px,py back to index into image, for positioning text box
	GetAxis/W=$win/Q bottom
	horiz = (px-V_min)/ (V_max-V_min)							// determine anchor for text box
	GetAxis/W=$win/Q left
	vert = (py-V_max)/(V_min-V_max)
	Variable x0 = horiz<0.5 ? 5 : -5,  y0 = vert<0.5 ? -5 : 5
	String anchor = SelectString(horiz<0.5,"R","L")+ SelectString(vert<0.5,"B","T")
	Tag/C/N=indexedPeakInfo/W=$win/A=$anchor/F=2/L=2/X=(x0)/Y=(y0)/P=1 $imageNameStr,index,tagStr
	s.doSetCursor= 1
	s.cursorCode = 21												// point marker
	DoUpdate
	return 1														// 1 means do not send back to Igor for more processing
End

//  ============================ End of spec Image Plotting =============================  //
//  ============================================================================  //



//  ============================================================================  //
//  ============================ Start of spec Image Loading ============================  //

Function/T LoadSpecImage(scanNum,pindex[,extras])
	Variable scanNum											// scanNumber in spec file
	Variable pindex												// point in scanNum
	String extras												// example		"quiet:1;multiple:0;"
	extras = SelectString(ParamIsDefault(extras),extras,"")
	Variable quiet = NumberByKey("quiet",extras)				// for quiet=1, if any problems just return
	Variable multiple = NumberByKey("multiple",extras)		// load multiple images, NOT YET WORKING
	quiet = numtype(quiet) ? 0 : !(!quiet)
	multiple = numtype(multiple) ? 0 : !(!multiple)

	if (!(scanNum>0) || !(pindex>=0) || numtype(scanNum+pindex))
		scanNum = scanNum>0 ? scanNum : NumVarOrDefault("root:Packages:spec:lastScan", 1)
		pindex = pindex>=0 ? pindex : 0
		Prompt scanNum,"spec Scan Number"
		Prompt pindex,"point in spec scan"
		DoPrompt "spec Scan & Point", scanNum,pindex
		if (V_flag)
			return ""
		endif

		if (!(scanNum>0) || numtype(scanNum))
			printf "scan number %g is invalid\r",scanNum
		endif
		if (!(pindex>=0) || numtype(pindex))
			printf "point number %g is invalid\r",pindex
		endif
	endif

	String 	specFldr
	sprintf specFldr,"%sspec%d:",StrVarOrDefault("root:Packages:spec:specDataFolder",""),scanNum
	if (!DataFolderExists(specFldr))
		printf "spec scan %g has not been loaded, you must load it first\r",scanNum
		return ""
	endif

	DFREF dfr=$specFldr										// Reference to folder with spec scan
	Wave ww=$(specFldr+GetIndexedObjNameDFR(dfr,1,0))
	if (!WaveExists(ww))
		return ""
	elseif (pindex<0 || pindex>=numpnts(ww) || numtype(pindex))
		printf "spec scan %g has been loaded, BUT it has only %d points, you asked for point=%g\r",scanNum,numpnts(ww),pindex
		return ""
	endif
	String fullName = specImageFileName(scanNum,pindex)		// full path+name to file
	GetFileFolderInfo/Q/Z=1 fullName							// check if file exists
	PathInfo home
	if (!V_isFile && V_flag)									// look for a path relative to "home" path, if home exists
//		String str = ":"+ParseFilePath(0,S_path,":",1,0)+":"
//		Variable ii = strsearch(fullName,str,0)
//		ii += strlen(str)
//		fullName = S_path+fullName[ii,Inf]
		fullName = S_path+fullName
		GetFileFolderInfo/P=home/Q/Z=1 fullName				// again, check if file exists
		if (!V_isFile)											// file does not exist, give up looking
			if (!quiet)
				printf "Could not find file  '%s'\r",fullName
			endif
			return ""
		endif
	endif

	String wnote, fldrSav= GetDataFolder(1)
	SetDataFolder $specFldr
	NewDataFolder/O/S images
	String imageName											// simple name of image in image folder
	if (ParamIsDefault(extras))
		wnote = LoadSpecImageInfo(scanNum,pindex)
		imageName = LoadGenericImageFile(fullName)
	else
		wnote = LoadSpecImageInfo(scanNum,pindex,extras=extras)
		imageName = LoadGenericImageFile(fullName,extras=extras)
	endif
	Variable preExists=(strsearch(imageName,"preExists:",0)==0)
	if (preExists)
		Wave image = $StringByKey("preExists", imageName)
	else
		Wave image = $imageName
	endif

	SetDataFolder fldrSav
	if (!quiet)
		printf "LoadSpecImage(%g,%g",scanNum,pindex
		if (!ParamIsDefault(extras))
			printf ", quiet=%g",quiet
		endif
		if (!ParamIsDefault(extras))
			printf ", multiple=%g",multiple
		endif
		if (WaveExists(image))
			printf ") // wave = '%s'\r",GetWavesDataFolder(image,2)+SelectString(preExists,"",",    a preExisting wave")
		else
			printf ") // NO wave loaded\r"
		endif
	endif
	if (!WaveExists(image))
		return ""
	endif

	Note/K image, wnote
	if (!quiet)
		GraphDiffractometerImage(image)
	endif
	return SelectString(preExists,"","preExists:")+GetWavesDataFolder(image,2)
End


Function/T LoadSpecImageROI(scanNum,pindex,i0,i1,j0,j1,[extras])
	Variable scanNum				// scanNumber in spec file
	Variable pindex					// point in scanNum
	Variable i0,i1,j0,j1			// pixel range of ROI (if i1 or j1<0 then use whole image)
	String extras					// example		"quiet:1;multiple:0;"
	extras = SelectString(ParamIsDefault(extras),extras,"")

	String imageName
	if (ParamIsDefault(extras))
		imageName = LoadSpecImage(scanNum,pindex)
	else
		imageName = LoadSpecImage(scanNum,pindex,extras=extras)
	endif
	Variable preExists=(strsearch(imageName,"preExists:",0)==0)
	if (preExists)
		Wave image = $StringByKey("preExists", imageName)
	else
		Wave image = $imageName
	endif
	if (!WaveExists(image))
		return ""
	endif
	String wnote = note(image)

	Variable ydim, xdim									// x,y size of whole array
	xdim = DimSize(image,0)
	ydim = DimSize(image,1)
	i1 = (i1<1) ? xdim-1 : i1							// -1 flags use whole range
	j1 = (j1<1) ? ydim-1 : j1

	i0 = max(round(i0),0)
	i1 = max(round(i1),0)
	j0 = max(round(j0),0)
	j1 = max(round(j1),0)

	Variable nx = i1-i0+1
	Variable ny = j1-j0+1
	if (nx<1 || ny<1)									// nothing to read
		KillWaves/Z image
		return ""
	endif

	Duplicate/R=[i0,i1][j0,j1]/FREE image, ROIwave	// save the ROI I want
	Redimension/N=(nx,ny,-1) image
	image = ROIwave									// re-set image to ROI values
	WAVEClear ROIwave

	wnote = ReplaceNumberByKey("startx",wnote,i0,"=")
	wnote = ReplaceNumberByKey("endx",wnote,i1,"=")
	wnote = ReplaceNumberByKey("starty",wnote,j0,"=")
	wnote = ReplaceNumberByKey("endy",wnote,j1,"=")	
	Note/K image,wnote
	return SelectString(preExists,"","preExists:")+GetWavesDataFolder(image,2)
End

Function/T LoadSpecImageInfo(scanNum,pindex[,extras])
	Variable scanNum											// scanNumber in spec file
	Variable pindex												// point in scanNum
	String extras												// example		"quiet:1;multiple:0;"
	extras = SelectString(ParamIsDefault(extras),extras,"")
	Variable multiple = NumberByKey("multiple",extras)		// load multiple images, NOT YET WORKING
	multiple = numtype(multiple) ? 0 : !(!multiple)
	Variable quiet = NumberByKey("quiet",extras)				// deault is quiet=1,  for quiet=1, if any problems just return ""
	quiet = numtype(quiet) ? 1 : !(!quiet)

	if (!(scanNum>0) || !(pindex>=0) || numtype(scanNum+pindex) && !quiet)
		scanNum = scanNum>0 ? scanNum : NumVarOrDefault("root:Packages:spec:lastScan", 1)
		pindex = pindex>=0 ? pindex : 0
		Prompt scanNum,"spec Scan Number"
		Prompt pindex,"point in spec scan"
		DoPrompt "spec Scan & Point", scanNum,pindex
		if (V_flag)
			return ""
		endif

		if (!(scanNum>0) || numtype(scanNum))
			printf "scan number %g is invalid\r",scanNum
		endif
		if (!(pindex>=0) || numtype(pindex))
			printf "point number %g is invalid\r",pindex
		endif
		quiet = 0
	endif
	String 	specFldr
	sprintf specFldr,"%sspec%d:",StrVarOrDefault("root:Packages:spec:specDataFolder",""),scanNum
	if (!DataFolderExists(specFldr))
		if (!quiet)
			printf "spec scan %g has not been loaded, you must load it first\r",scanNum
		endif
		return ""
	endif

	DFREF dfr=$specFldr										// Reference to folder with spec scan
	Wave ww=$(specFldr+GetIndexedObjNameDFR(dfr,1,0))
	if (!WaveExists(ww))
		return ""
	elseif (pindex<0 || pindex>=numpnts(ww) || numtype(pindex))
		printf "spec scan %g has been loaded, BUT it has only %d points, you asked for point=%g\r",scanNum,numpnts(ww),pindex
		return ""
	endif
	String fullName = specImageFileName(scanNum,pindex)		// full path+name to file
	GetFileFolderInfo/Q/Z=1 fullName							// check if file exists
	PathInfo home
	if (!V_isFile && V_flag)									// look for a path relative to "home" path, if home exists
//		String str = ":"+ParseFilePath(0,S_path,":",1,0)+":"
//		Variable ii = strsearch(fullName,str,0)
//		ii += strlen(str)
//		fullName = S_path+fullName[ii,Inf]
		fullName = S_path+fullName
		GetFileFolderInfo/P=home/Q/Z=1 fullName				// again, check if file exists
		if (!V_isFile)											// file does not exist, give up looking
			if (!quiet)
				printf "Could not find file  '%s'\r",fullName
			endif
			return ""
		endif
	endif

	String list = ReadGenericHeader(fullName)
//	String list
//	if (strlen(extras))
//		list = ReadGenericHeader(fullName,extras=extras)
//	else
//		list = ReadGenericHeader(fullName)
//	endif

	if (!quiet)
		printf "LoadSpecImageInfo(%g,%g",scanNum,pindex
		if (!ParamIsDefault(extras))
			printf ", quiet=%g",quiet
		endif
		if (!ParamIsDefault(extras))
			printf ", multiple=%g",multiple
		endif
		printf ") // wave = '%s'\r",fullName
	endif
	if (strlen(list)==0)				// probably could not find the image, so nothing read
		if (!quiet)
			print "ERROR  -- Nothing read, probably could not find the image"
		endif
		return ""
	endif

	String filterTransWaveName = GetSpecialWaveName_spec(scanNum,"filterTrans")
	String I0WaveName = GetSpecialWaveName_spec(scanNum,"I0WaveName")
	String EnergyWaveName = GetSpecialWaveName_spec(scanNum,"EnergyWaveName")

	list = ReplaceNumberByKey("scanNum",list,scanNum,"=")
	list = ReplaceNumberByKey("scanPoint",list,pindex,"=")

	String motors = StrVarOrDefault("root:Packages:Diffractometer:DiffractometerAxisNames",""), mot
	Variable i,val, Nm=ItemsInList(motors)
	for (i=0;i<Nm;i+=1)
		mot = StringFromList(i,motors)
		val = specInfo(scanNum,mot)
		if (numtype(val)==0)
			list = ReplaceNumberByKey(mot,list,val,"=")
		endif
	endfor

	String str
	Variable N=CountObjects(specFldr,1)
	if (N>0)
		for (i=0;i<N;i+=1)
			str = GetIndexedObjName(specFldr,1,i)
			if (StringMatch(str,"Epoch"))					// stop when I get to Epoch
				break
			endif
			Wave ww=$(specFldr+str)						// since this is before Epoch, it must be a positioner
			if (WaveExists(ww))							// save positions in wavenote
				list = ReplaceNumberByKey(str,list,ww[pindex],"=")
			endif
		endfor
	endif

	if (exists(specFldr+"filters"))							// special for filters
		Wave ww=$(specFldr+"filters")
		Variable filters = ww[pindex]
		filters = numtype(filters) ? 0 : filters
		if (filters>=0)
			list = ReplaceNumberByKey("filters",list,filters,"=")
		endif
	endif
	if (exists(specFldr+EnergyWaveName)==1)			// special for energy
		Wave ww=$(specFldr+EnergyWaveName)
		Variable keV = ww[pindex]
		if (numtype(keV)==0)
			list = RemoveByKey("DCM_energy",list,"=")
			list = ReplaceNumberByKey("keV",list,keV,"=")
		endif
	endif
	Variable filterTrans=NaN
	if (exists(specFldr+filterTransWaveName)==1)		// special for computed tranmission from filters
		Wave ww=$(specFldr+filterTransWaveName)
		filterTrans = ww[pindex]
	endif
	filterTrans = numtype(filterTrans) ? specInfo(scanNum,"PF4_trans") : filterTrans
	if (filterTrans>0 && filterTrans<=1)
		list = ReplaceNumberByKey("filterTrans",list,filterTrans,"=")
	endif
	if (exists(specFldr+I0WaveName)==1)					// special for I0
		Wave ww=$(specFldr+I0WaveName)
		Variable I0 = ww[pindex]
		if (I0>0 && numtype(I0)==0)
			list = ReplaceNumberByKey("I0",list,I0,"=")
		endif
	endif

	STRUCT sampleStructure sa	
//	str = ParseFilePath(1,GetWavesDataFolder(image,1),":",1,0)+"sampleStructStr"
//	String strStruct=StrVarOrDefault(str,"")				// fill the sample structure with values in spec scan directory
	String strStruct=StrVarOrDefault(specFldr+"sampleStructStr","")	// fill the sample structure with values in spec scan directory
	Make/N=9/D/FREE UBs=NaN
	if (strlen(strStruct)>1)
		StructGet/S/B=2 sa, strStruct						// found structure information, load into s
		UBs = sa.UB[p]
	else														// could not get recip from sampleStructure, use specInfo
		str = specInfoT(scanNum,"UB")
		str = TrimFrontBackWhiteSpace(str)
		str = ReplaceString("  ",str," ")					// change all double spaces to single spaces
		if (ItemsInList(str," ")==9)
			Make/N=9/D/FREE UBs
			UBs = str2num(StringFromList(p,str," "))
//print str
			Redimension/N=(3,3) UBs
			String DiffractometerName = StrVarOrDefault("root:Packages:Diffractometer:DiffractometerName","XYZ")
			String name = "diffractometer#"+DiffractometerName+"Spec2BL"
			FUNCREF protoSpec2BL getSpec2BL = $name
			Wave spec2BL = getSpec2BL(0)
			MatrixOP/FREE UBs = spec2BL x UBs * 10		// spec uses Angstrom, I use nm
			Redimension/N=(9) UBs
		endif
	endif
	if (numtype(sum(UBs))==0)							// save the reciprocal lattice in image wave note
		sprintf str "{{%g,%g,%g}{%g,%g,%g}{%g,%g,%g}}",UBs[0],UBs[1],UBs[2], UBs[3],UBs[4],UBs[5], UBs[6],UBs[7],UBs[8]
		list = ReplaceStringByKey("recip_lattice0",list,str,"=")
	endif

#ifdef USE_EXTRA_PILATUS_SPEC_FILE
	String moreSpec = fetchPilatusInfoExtraSpec(fullName,scanNum=scanNum,pindex=pindex)	// only gets information from the log file
	if (strlen(moreSpec))
		list = MergeKeywordLists(list,moreSpec,0,"=","")	// give priority to list
	endif
#endif
	return list
End


Function/T specImageFileName(scanNum,pindex,[specPath])
	Variable scanNum
	Variable pindex
	String specPath
	specPath = SelectString(ParamIsDefault(specPath),specPath,StrVarOrDefault("root:Packages:spec:specDefaultFile",""))
	String fullName=""
	if (exists("PilatusTiffReadHeader")==6)
		String path = ParseFilePath(1, specPath,":", 1,0)
		String name=ParseFilePath(3, specPath,":", 0,0)
		sprintf fullName, "%simages:%s:S%03d:%s_S%03d_%05d.tif", path,name,scanNum,name,scanNum,pindex
		// sprintf fullName, "%simages:%s:S%03d:%s_S%03d_%05d.tif", ParseFilePath(1, specPath,":", 1,0),name,scanNum,name,scanNum,pindex
	elseif (exists("ReadHDF5header")==6)
		fullName=""						// not yet implemented
	elseif (exists("TiffReadHeader")==6)
		fullName=""						// not yet implemented
	elseif (exists("WinViewReadHeader")==6)
		fullName=""						// not yet implemented
	endif
	return fullName
End



Static Function SetGlobalSpecialSpecNames()			// set global values used find things like Io, filter transmission, ...
	String I0WaveName=StrVarOrDefault("root:Packages:spec:I0WaveName","I0")						// name of wave in spec that holds I0 to be used for normalizations
	String filterTransWaveName=StrVarOrDefault("root:Packages:spec:filterTransWaveName","trans")	// name of wave in spec that holds filter transmission
	String EnergyWaveName=StrVarOrDefault("root:Packages:spec:EnergyWaveName","Energy")	// name of wave in spec that holds mono energy

	Prompt I0WaveName,"spec wave containing the incident beam intensity, used for normalizing"
	Prompt filterTransWaveName,"spec wave containing the filter tramsmission, used for normalizing"
	Prompt EnergyWaveName,"spec wave containing the mono energy"
	DoPrompt "Special spec Wave Names", I0WaveName,filterTransWaveName, EnergyWaveName
	if (V_flag)
		return 1
	endif
	String/G root:Packages:spec:I0WaveName = I0WaveName
	String/G root:Packages:spec:filterTransWaveName = filterTransWaveName
	String/G root:Packages:spec:EnergyWaveName = EnergyWaveName
	return 0
End


// This function is needed because it will in future also have code that looks in the spec file to get the result from each scan in a spec file
ThreadSafe Function/T GetSpecialWaveName_spec(scanNum,id)
	Variable  scanNum
	String id				// can be "filterTrans" or "I0" or "Io"

	String name = specInfoT(scanNum,id)					// try to get the name from this spec scan
	if (strlen(name)<1)
		name = StrVarOrDefault("root:Packages:spec:"+id,"")
	endif
//	String name = specInfoT(scanNum,id+"Name")			// try to get the name from this spec scan
//	if (strlen(name)<1)
//		name = StrVarOrDefault("root:Packages:spec:"+id+"WaveName","")
//	endif
	return name
End


ThreadSafe Function/WAVE getRLfrom3DWave(scatter,point)		// find reciprocal lattice, this is called by external routines
	Wave scatter
	Variable point
	if (!WaveExists(scatter))
		return $""
	endif
	if (WaveDims(scatter)==3)
		String recip_lattice0=StringByKey("recip_lattice0",note(scatter),"=")
		Wave RL = matString2mat33(recip_lattice0)
	else
		Wave RL=$""
	endif
	return RL
End
//
ThreadSafe Static Function/WAVE matString2mat33(str)			// returns a FREE wave with (3x3) matrix
	String str
	Variable a0,a1,a2, b0,b1,b2, c0,c1,c2
	sscanf str, "{{%g,%g,%g}{%g,%g,%g}{%g,%g,%g}}", a0,a1,a2, b0,b1,b2, c0,c1,c2
	if (V_flag==9)
		Make/N=(3,3)/D/FREE mat
		mat[0][0]= {a0,a1,a2}					// the reciprocal lattice
		mat[0][1]= {b0,b1,b2}
		mat[0][2]= {c0,c1,c2}
		return mat
	else
		return $""
	endif
End

//  ============================ End of spec Image Loading =============================  //
//  ============================================================================  //



//  ============================================================================  //
//  ============================= Start of spec Image Init ==============================  //

Function init_specImage(imageType)
	String imageType					// image type to include if none avaiable
	InitGizmoZoomTranslate()
	specInitPackage()
	if (imageLoadersAvailable())
		return 0
	endif
//	if (exists("root:Packages:spec:I0WaveName")!=2 || exists("root:Packages:spec:filterTransWaveName")!=2 || exists("root:Packages:spec:EnergyWaveName")!=2)
//		SetGlobalSpecialSpecNames()
//	endif
	if (exists("root:Packages:Diffractometer:DiffractometerName")!=2 || exists("root:Packages:Diffractometer:DiffractometerAxisNames")!=2 || exists("root:Packages:Diffractometer:Naxes")!=2)
		diffractometer#SelectDiffractometer("")
		diffractometer#SetDetectorParameters(NaN,NaN,NaN,NaN,"","","")
	endif

	String str,ipf,types="",vers="", typeVers="PilatusTiffs:2.10;HDF5:0.26;WinView:2.10;Tiff simple:1.14"
	Variable i, ver
	for (i=0;i<ItemsInList(typeVers);i+=1)
		str = StringFromList(i,typeVers)
		types += StringFromList(0,str,":")+";"
		vers += StringFromList(1,str,":")+";"
	endfor

	if (WhichListItem(imageType,types)<0)
		Prompt imageType,"Image Type",popup,types
		DoPrompt "Select an image type",imageType
		if (V_flag)
			DoAlert 0, "You Do NOT have any image loaders available"
			return 1
		endif
	endif
	ver = NumberByKey(imageType,typeVers)
	if (WhichListItem(imageType,types)<0)
		DoAlert 0, "You Do NOT have any image loaders available"
		return 1
	endif
	ipf = SelectString(StringMatch(imageType,"HDF5"),imageType,"HDF5images")
	ipf = SelectString(StringMatch(ipf,"Tiff simple"),ipf,"Tiff")

	if (ver>0)
		sprintf str, "INSERTINCLUDE  \"%s\", version>=%g",ipf,ver
	else
		sprintf str, "INSERTINCLUDE  \"%s\"",ipf
	endif
	Execute/P str
	Execute/P "COMPILEPROCEDURES "
	print "Included support for ", imageType

	if (exists("root:Packages:spec:I0WaveName")!=2 || exists("root:Packages:spec:filterTransWaveName")!=2 || exists("root:Packages:spec:EnergyWaveName")!=2)
		SetGlobalSpecialSpecNames()
	endif
	return 0
End

//  ============================== End of spec Image Init ==============================  //
//  ============================================================================  //
