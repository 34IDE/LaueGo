#pragma rtGlobals=3		// Use modern global access method.
#pragma version = 0.06
#pragma IgorVersion = 6.3
#pragma ModuleName=Stats3D



Function/S PeakMax3D(vol, [HW,printIt])	// finds max of data in a 3D volume
	Wave vol
	Variable HW										// HW used for the fit
	Variable printIt
	HW = ParamIsDefault(HW) || numtype(HW) || HW<=0 ? NaN : HW
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRtStackInfo(2))==0 : printIt
	if (!(WaveDims(vol)==3 && numpnts(vol)>0))						// 3D volume array
		return ""
	endif

	WaveStats/M=1/Q vol
	Variable maxVal = V_max
	Make/D/FREE xyzMax={V_maxRowLoc,V_maxColLoc, V_maxLayerLoc}
	Make/N=3/I/FREE ijk = round( (xyzMax[p] - DimOffset(vol,p)) / DimDelta(vol,p) )
	Variable npnt = ijk[2]*DimSize(vol,1)*DimSize(vol,0) + ijk[1]*DimSize(vol,0) + ijk[0]
	if (numtype(HW) || HW<=0)
		Make/N=(3)/D/FREE sizes = (DimSize(vol,p)-1) * DimDelta(vol,p)
		HW = WaveMin(sizes)/2
	endif

	if (printIt)
		String name=SelectString(WaveType(vol,2)==1, "", NameOfWave(vol))
		printf "max = %g   @ %s,  %s[%d][%d][%d],  pnt#=%d\r",maxVal, vec2str(xyzMax,sep=", "), NameOfWave(vol),ijk[0],ijk[1],ijk[2],npnt
	endif

	String list=FitPeakIn3D(vol,xyzMax,HW,printIt=printIt)
	list = ReplaceNumberByKey("iStart",list,ijk[0],"=")
	list = ReplaceNumberByKey("jStart",list,ijk[1],"=")
	list = ReplaceNumberByKey("kStart",list,ijk[2],"=")
	return list
End


Function IntegralOfVolume(subVol,[bkg])
	// returns the Integral of a volume, Integral{ f(x,y,z) dx dy dz }
	// this can optionally remove a bkg value from all the valid points
	Wave subVol
	Variable bkg
	bkg = ParamIsDefault(bkg) || numtype(bkg) ? 0 : bkg	// default bkg is zero, bkg=0 integrates everything
	//	Sometimes, a large part of subVol may be zero because no pixels correspond to that voxel, 
	//	they should not be included in this integral. 
	//	So to subtract the bkg, I need to know the number of voxels actually used.
	//	If bkg==0, then none is subtracted, so it all becomes simple again.

	Variable NusedVoxels=NumberByKey("NusedVoxels",note(subVol),"=")
	NusedVoxels = NusedVoxels<=0 ? NaN : NusedVoxels


	if (bkg && numtype(NusedVoxels))							// a bkg was passed, but cannot find NusedVoxels
		// did not provide NusedVoxels, but have bkg, ignore all voxels that are < (1/5 bkg)
		Duplicate/FREE subVol, measuredVol
		measuredVol = measuredVol<=(bkg/5) ? NaN : measuredVol
		WaveStats/M=1/Q measuredVol
		NusedVoxels = V_npnts
		WaveClear measuredVol
	else																	// bkg is zero or NusedVoxels was passed
		WaveStats/M=1/Q subVol
		NusedVoxels = numtype(NusedVoxels) ? V_npnts : NusedVoxels
	endif
	Variable total = V_Sum - (NusedVoxels*bkg)				// sum of desired voxels
	total *= DimDelta(subVol,0)*DimDelta(subVol,1)*DimDelta(subVol,2)	// convert sum to integral
	return total
End


Function/S FitPeakIn3D(space3D,startXYZ,HWx,[HWy,HWz, stdDev, printIt])
	Wave space3D
	Wave startXYZ						// starting xyz for the fit
	Variable HWx,HWy,HWz			// starting half widths dx, dy, dz for the fit
	Wave stdDev							// errors in space3D, standard deviation of each value in space3D
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt)? strlen(GetRTStackInfo(2))==0 : !(!printIt)
	HWy = ParamIsDefault(HWy) || numtype(HWy) || HWy<=0 ? HWx : HWy	// HWy & HWz default to HWx
	HWz = ParamIsDefault(HWz) || numtype(HWz) || HWz<=0 ? HWx : HWz

	if (!WaveExists(space3D) || !WaveExists(startXYZ))
		return ""
	elseif (WaveDims(space3D)!=3 || numpnts(startXYZ)!=3 || numtype(sum(startXYZ)))
		return ""
	elseif (numtype(HWx+HWy+HWz) || HWx<=0 || HWy<=0 || HWz<=0)
		return ""
	endif

	Make/N=3/D/FREE HW={HWx,HWy,HWz}			// half widths in x, y, & z for start of fit
	if (WaveExists(stdDev))
		Wave errWave = stdDev
	else
		Duplicate/FREE space3D, errWave		// this assumes that one count in detector == 1 photon
		errWave = sqrt(space3D)				// this is shot noise weighting
	endif

	// set the starting point for the fit
	Variable maxVal=WaveMax(space3D),minVal=WaveMin(space3D)
	Make/N=8/D/O W_coef = {NaN,NaN,startXYZ[0],HW[0]/4,startXYZ[1],HW[1]/4,startXYZ[2],HW[2]/4}
	if (maxVal>-minVal)			// a positive peak
		W_coef[1] = maxVal
		W_coef[0] = minVal==0 ? 1 : minVal
	else								// a negative peak
		W_coef[1] = minVal
		W_coef[0] = maxVal==0 ? 1 : maxVal
	endif

	Variable V_FitOptions = (printIt ? 0 : 4)
	Variable V_FitError=0, V_FitQuitReason=0
	FuncFitMD/Q Stats3D#Gaussian3DFitFunc, W_coef, space3D/W=errWave/I=1
	Variable chisq = V_chisq
	Wave W_sigma=W_sigma
	Make/N=3/T/FREE units=WaveUnits(space3D,p)
	Make/N=3/D/FREE xyzFit=W_coef[2*p + 2]

	Variable err=V_FitError, swap
	err = err || !isPointInVolume(space3D,xyzFit)	// verify if point is in a volume, works for triplets and for a 3D array

	if (err)
		print "ERROR -- Fit failed, peak not inside of fitting volume"
		print FitErrorString(V_FitError,V_FitQuitReason)
		return ""
	endif

	Wave hkl=$""
#if exists("diffractometer#sample2crystal")
	STRUCT sampleStructure sa	
	String str = ParseFilePath(1,GetWavesDataFolder(space3D,1),":",1,0)+"sampleStructStr"
	String strStruct=StrVarOrDefault(str,"")				// fill the sample structure with values in spec scan directory
	StructGet/S/B=2 sa, strStruct								// found structure information, load into s
	Wave hkl = diffractometer#sample2crystal(sa,xyzFit)	// rotate qvec from sample-location frame into crystal based frame, the hkl
#endif
	if (!WaveExists(hkl))
		FUNCREF getRLfrom3DWaveProto getRL=$"getRLfrom3DWave"
		Wave RL = getRl(space3D,NaN)
		if (WaveExists(RL))
			MatrixOP/FREE hkl = Inv(RL) x xyzFit
		endif
	endif

	if (printIt)
		units = SelectString(stringmatch(units[p],"nm\\S-1*"),units[p],"1/nm")		// two different ways of writing 1/nm
		units = SelectString(strlen(units[p]),""," ("+units[p]+")")						// add () when something present
		Make/N=3/T/FREE Qstr=SelectString(strsearch(units[p],"1/nm",0)>=0,"","Q")	// is it Q?
		printf "Gaussian Fit has offset = %s,  amp = %s,  chisq = %g\r",ValErrStr(W_coef[0],W_sigma[0],sp=1),ValErrStr(W_coef[1],W_sigma[1],sp=1),chisq
		printf "  <%sx> = %s%s,   FWHMx = %s%s\r",Qstr[0],ValErrStr(xyzFit[0],W_sigma[2],sp=1),units[0],ValErrStr(W_coef[3],W_sigma[3],sp=1),units[0]
		printf "  <%sy> = %s%s,   FWHMy = %s%s\r",Qstr[1],ValErrStr(xyzFit[1],W_sigma[4],sp=1),units[1],ValErrStr(W_coef[5],W_sigma[5],sp=1),units[1]
		printf "  <%sz> = %s%s,   FWHMz = %s%s\r",Qstr[2],ValErrStr(xyzFit[2],W_sigma[6],sp=1),units[1],ValErrStr(W_coef[7],W_sigma[7],sp=1),units[2]
		if (stringmatch(units[0],units[1]) && stringmatch(units[0],units[2]))
			Variable xyzMag=norm(xyzFit), dxyzMag=sqrt(W_sigma[2]^2 + W_sigma[4]^2 + W_sigma[6]^2)
			printf "  <|%s|> = %s%s\r",Qstr[0],ValErrStr(xyzMag,dxyzMag,sp=1),units[0]
		endif

		Make/N=3/I/FREE ijk, Nxyz=DimSize(space3D,p)
		ijk = round( (xyzFit[p]-DimOffset(space3D,p)) / DimDelta(space3D,p) )
		ijk = limit(round(ijk[p]),0,Nxyz[p]-1)
		Variable i=ijk[0], j=ijk[1], k=ijk[2]
		Variable N = k*Nxyz[0]*Nxyz[1] + j*Nxyz[0] + i
		printf "  Closest point to peak is the %s[%g, %g, %g] or [%d] = %g\r",NameOfWave(space3D),i,j,k,N,space3D[i][j][k]

		if (WaveExists(hkl))
			printf "hkl Peak Center = %s\r",vec2str(hkl)
		endif
	endif

	String keyVals=""
	keyVals = ReplaceNumberByKey("offset",keyVals,W_coef[0],"=")
	keyVals = ReplaceNumberByKey("offsetErr",keyVals,W_sigma[0],"=")
	keyVals = ReplaceNumberByKey("amp",keyVals,W_coef[1],"=")
	keyVals = ReplaceNumberByKey("ampErr",keyVals,W_sigma[1],"=")
	keyVals = ReplaceNumberByKey("Xc",keyVals,W_coef[2],"=")
	keyVals = ReplaceNumberByKey("XcErr",keyVals,W_sigma[2],"=")
	keyVals = ReplaceNumberByKey("Yc",keyVals,W_coef[4],"=")
	keyVals = ReplaceNumberByKey("YcErr",keyVals,W_sigma[4],"=")
	keyVals = ReplaceNumberByKey("Zc",keyVals,W_coef[6],"=")
	keyVals = ReplaceNumberByKey("ZcErr",keyVals,W_sigma[6],"=")
	keyVals = ReplaceNumberByKey("FWX",keyVals,W_coef[3],"=")
	keyVals = ReplaceNumberByKey("FWXErr",keyVals,W_sigma[3],"=")
	keyVals = ReplaceNumberByKey("FWY",keyVals,W_coef[5],"=")
	keyVals = ReplaceNumberByKey("FWYErr",keyVals,W_sigma[5],"=")
	keyVals = ReplaceNumberByKey("FWZ",keyVals,W_coef[7],"=")
	keyVals = ReplaceNumberByKey("FWZErr",keyVals,W_sigma[7],"=")
	keyVals = ReplaceNumberByKey("chisq",keyVals,V_chisq,"=")
	if (strlen(units[0]))
		keyVals = ReplaceStringByKey("Xunit",keyVals,units[0],"=")
	endif
	if (strlen(units[1]))
		keyVals = ReplaceStringByKey("Yunit",keyVals,units[1],"=")
	endif
	if (strlen(units[2]))
		keyVals = ReplaceStringByKey("Zunit",keyVals,units[2],"=")
	endif
	if (WaveExists(hkl))
		keyVals = ReplaceStringByKey("hklPeakCenter",keyVals,vec2str(hkl,sep=","),"=")
	endif
	return keyVals
End



Function/S XX_FitPeakAt3Dmarker(space3D,Qc,QxHW,[QyHW,QzHW,printIt])
	Wave space3D
	Wave Qc								// center of sub-volume to fit
	Variable QxHW,QyHW,QzHW		// half widths dQz, dQy, dQz for the sub volume
	Variable printIt

	printIt = ParamIsDefault(printIt) ? 0 : !(!printIt)
	QyHW = ParamIsDefault(QyHW) ? NaN : QyHW	// QyHW & QzHW default to QwhX
	QzHW = ParamIsDefault(QzHW) ? NaN : QzHW
	if (WaveExists(space3D))
		if (WaveDims(space3D)!=3)
			return ""					// quit if really bad input passed
		endif
	endif
	if (WaveExists(Qc))
		if (numpnts(Qc)!=3)
			return ""					// quit if really bad input passed
		endif
	endif
	String spaceList = WaveListClass("GizmoXYZ;Qspace3D*","*","DIMS:3")
	if (!WaveExists(space3D))					// no space3D passed
		if (ItemsInList(spaceList)==1)		// only 1 choice, use it
			Wave space3D = $StringFromList(0,spaceList)
		elseif (ItemsInList(spaceList)<1)	// no acceptable choices, quit
			return ""
		endif
	endif
	String QcList=WaveList("*",";","DIMS:2,MAXROWS:1,MINROWS:1,MAXCOLS:3,MINCOLS:3" )+WaveList("*",";","DIMS:1,MAXROWS:3,MINROWS:3" )
	if (!WaveExists(Qc))						// no Qc passed
		if (ItemsInList(QcList)==1)			// only 1 choice, use it
			Wave Qc = $StringFromList(0,QcList)
		elseif (ItemsInList(QcList)<1)		// no acceptable choices, quit
			return ""
		endif
	endif
	if (numpnts(Qc)==3 && WaveType(Qc,2)==2)	// a free wave was passed, deal with it
		QcList = "_free_;"+QcList
	endif

	if (!WaveExists(space3D) || !WaveExists(Qc) || !(QxHW>0))
		printIt = 1
		String QcName=NameOfWave(Qc), spaceName = NameOfWave(space3D)
		QcName = SelectString(strlen(QcName),"gizmoScatterMarker",QcName)
		QxHW = QxHW>0 ? QxHW : 1
		Prompt spaceName,"3D space",popup,spaceList
		Prompt QcName,"Center of 3D space",popup,QcList
		Prompt QxHW,"X-HW in Volume to Use"
		Prompt QyHW,"Y-HW in Volume to Use (NaN defaults to X-HW)"
		Prompt QzHW,"Z-HW in Volume to Use (NaN defaults to X-HW)"
		DoPrompt "Fit 3D Peak",spaceName,QcName,QxHW,QyHW,QzHW
		if (V_flag)
			return ""
		endif
		printf "FitPeakAt3Dmarker(%s, %s, %g",spaceName,QcName,QxHW
		if (QyHW>0)
			printf ", QyHW=%g",QyHW
		endif
		if (QzHW>0)
			printf ", QzHW=%g",QzHW
		endif
		printf ")\r"
		Wave space3D=$spaceName
		if (!stringmatch(QcName,"_free_"))
			Wave Qc=$QcName						// this is for a free wave which does not have a name
		endif
	endif
	QyHW = QyHW>0 ? QyHW : QxHW
	QzHW = QzHW>0 ? QzHW : QxHW
	if (!WaveExists(space3D) || !WaveExists(Qc) || !(QxHW>0) || !(QyHW>0) || !(QzHW>0))
		return ""
	endif

	Make/N=3/D/FREE Qhw={QxHW,QyHW,QzHW}	// half widths dQz, dQy, dQz for the sub volume
	if (!WaveExists(space3D) || !WaveExists(Qc))
		return ""
	elseif (numtype(sum(Qhw)+sum(Qc)))
		return ""
	endif

	Make/N=3/D/FREE Qlo, Qhi, iLo,iHi,Np
	Qlo = Qc[p]-Qhw[p]
	Qhi = Qc[p]+Qhw[p]
	Np = DimSize(space3D,p)
	iLo = (Qlo[p]-DimOffSet(space3D,p))/DimDelta(space3D,p)
	iHi = (QHi[p]-DimOffSet(space3D,p))/DimDelta(space3D,p)
	iLo = limit(floor(iLo),0,Np-1)
	iHi = limit(ceil(iHi),0,Np-1)
	Np = iHi[p]-iLo[p]+1
	if (0 && printIt)
		printWave(Qlo,Name="Qlo")
		printWave(Qhi,Name="Qhi")
		printWave(iLo,Name="iLo")
		printWave(iHi,Name="iHi")
		printWave(Np,Name="Np")
	endif
	if (WaveMin(Np)<=0)
		return ""
	endif
	Make/N=(Np[0],Np[1],Np[2])/FREE/D sub3D, stdDev
	sub3D = space3D[iLo[0]+p][iLo[1]+q][iLo[2]+r]
	SetScale/P x, iLo[0]*DimDelta(space3D,0)+DimOffset(space3D,0), DimDelta(space3D,0),"",sub3D
	SetScale/P y, iLo[1]*DimDelta(space3D,1)+DimOffset(space3D,1), DimDelta(space3D,1),"",sub3D
	SetScale/P z, iLo[2]*DimDelta(space3D,2)+DimOffset(space3D,2), DimDelta(space3D,2),"",sub3D
	stdDev = sub3D					// this assumes that one count in detector == 1 photon

	// set the starting point for the fit
	Variable maxVal=WaveMax(sub3D),minVal=WaveMin(sub3D)
	Make/N=8/D/O W_coef = {NaN,NaN,Qc[0],Qhw[0]/4,Qc[1],Qhw[1]/4,Qc[2],Qhw[2]/4}
	if (maxVal>-minVal)			// a positive peak
		W_coef[1] = maxVal
		W_coef[0] = minVal==0 ? 1 : minVal
	else								// a negative peak
		W_coef[1] = minVal
		W_coef[0] = maxVal==0 ? 1 : maxVal
	endif

//	Variable V_FitOptions=2, V_FitError=0, V_FitQuitReason=0		// V_FitOptions=2 means robust fitting
	Variable V_FitError=0, V_FitQuitReason=0
	FuncFitMD/Q Stats3D#Gaussian3DFitFunc, W_coef, sub3D/W=stdDev/I=1
	Variable chisq = V_chisq
	Wave W_sigma=W_sigma
	Make/N=3/T/FREE units=WaveUnits(space3D,p)
	Make/N=3/D/FREE Qo=W_coef[2*p + 2]

	Variable err = V_FitError
	err = err || !( abs(Qc[0]-Qo[0]) < Qhw[0] )
	err = err || !( abs(Qc[1]-Qo[1]) < Qhw[0] )
	err = err || !( abs(Qc[2]-Qo[2]) < Qhw[0] )
	if (err)
		print "ERROR -- Fit failed, peak not inside of fitting volume"
		print FitErrorString(V_FitError,V_FitQuitReason)
		return ""
	endif

	Wave hkl=$""
#if exists("diffractometer#sample2crystal")
	STRUCT sampleStructure sa	
	String str = ParseFilePath(1,GetWavesDataFolder(space3D,1),":",1,0)+"sampleStructStr"
	String strStruct=StrVarOrDefault(str,"")				// fill the sample structure with values in spec scan directory
	StructGet/S/B=2 sa, strStruct								// found structure information, load into s
	Wave hkl = diffractometer#sample2crystal(sa,Qo)		// rotate qvec from sample-location frame into crystal based frame, the hkl
#endif
	if (!WaveExists(hkl))
		FUNCREF getRLfrom3DWaveProto getRL=$"getRLfrom3DWave"
		Wave RL = getRl(space3D,NaN)
		if (WaveExists(RL))
			MatrixOP/FREE hkl = Inv(RL) x Qo
		endif
	endif

	if (printIt)
		units = SelectString(stringmatch(units[p],"nm\\S-1*"),units[p],"1/nm")		// two different ways of writing 1/nm
		units = SelectString(strlen(units[p]),""," ("+units[p]+")")						// add () when something present
		Make/N=3/T/FREE Qstr=SelectString(strsearch(units[p],"1/nm",0)>=0,"","Q")	// is it Q?
		printf "Gaussian Fit has offset = %s,  amp = %s,  chisq = %g\r",ValErrStr(W_coef[0],W_sigma[0],sp=1),ValErrStr(W_coef[1],W_sigma[1],sp=1),chisq
		printf "  <%sx> = %s%s,   FWHMx = %s%s\r",Qstr[0],ValErrStr(Qo[0],W_sigma[2],sp=1),units[0],ValErrStr(W_coef[3],W_sigma[3],sp=1),units[0]
		printf "  <%sy> = %s%s,   FWHMy = %s%s\r",Qstr[1],ValErrStr(Qo[1],W_sigma[4],sp=1),units[1],ValErrStr(W_coef[5],W_sigma[5],sp=1),units[1]
		printf "  <%sz> = %s%s,   FWHMz = %s%s\r",Qstr[2],ValErrStr(Qo[2],W_sigma[6],sp=1),units[1],ValErrStr(W_coef[7],W_sigma[7],sp=1),units[2]
		if (stringmatch(units[0],units[1]) && stringmatch(units[0],units[2]))
			Variable Qmag=norm(Qo), dQmag=sqrt(W_sigma[2]^2 + W_sigma[4]^2 + W_sigma[6]^2)
			printf "  <|%s|> = %s%s\r",Qstr[0],ValErrStr(Qmag,dQmag,sp=1),units[0]
		endif
		if (WaveDims(space3D)==3)
			Np = DimSize(space3D,p)
			Make/N=3/D/FREE ijk
			ijk = (Qo[p]-DimOffset(space3D,p)) / DimDelta(space3D,p)
			ijk = limit(round(ijk[p]),0,Np[p]-1)
			Variable i=ijk[0], j=ijk[1], k=ijk[2]
			Variable N = k*Np[0]*Np[1] + j*Np[0] + i
			printf "  Closest point to peak is the %s[%g, %g, %g] or [%d] = %g\r",NameOfWave(space3D),i,j,k,N,space3D[i][j][k]
		endif
		if (WaveExists(hkl))
			printf "hkl Peak Center = %s\r",vec2str(hkl)
		endif
	endif

	String keyVals=""
	keyVals = ReplaceNumberByKey("offset",keyVals,W_coef[0],"=")
	keyVals = ReplaceNumberByKey("offsetErr",keyVals,W_sigma[0],"=")
	keyVals = ReplaceNumberByKey("amp",keyVals,W_coef[1],"=")
	keyVals = ReplaceNumberByKey("ampErr",keyVals,W_sigma[1],"=")
	keyVals = ReplaceNumberByKey("Xc",keyVals,W_coef[2],"=")
	keyVals = ReplaceNumberByKey("XcErr",keyVals,W_sigma[2],"=")
	keyVals = ReplaceNumberByKey("Yc",keyVals,W_coef[4],"=")
	keyVals = ReplaceNumberByKey("YcErr",keyVals,W_sigma[4],"=")
	keyVals = ReplaceNumberByKey("Zc",keyVals,W_coef[6],"=")
	keyVals = ReplaceNumberByKey("ZcErr",keyVals,W_sigma[6],"=")
	keyVals = ReplaceNumberByKey("FWX",keyVals,W_coef[3],"=")
	keyVals = ReplaceNumberByKey("FWXErr",keyVals,W_sigma[3],"=")
	keyVals = ReplaceNumberByKey("FWY",keyVals,W_coef[5],"=")
	keyVals = ReplaceNumberByKey("FWYErr",keyVals,W_sigma[5],"=")
	keyVals = ReplaceNumberByKey("FWZ",keyVals,W_coef[7],"=")
	keyVals = ReplaceNumberByKey("FWZErr",keyVals,W_sigma[7],"=")
	keyVals = ReplaceNumberByKey("chisq",keyVals,V_chisq,"=")
	if (strlen(units[0]))
		keyVals = ReplaceStringByKey("Xunit",keyVals,units[0],"=")
	endif
	if (strlen(units[1]))
		keyVals = ReplaceStringByKey("Yunit",keyVals,units[1],"=")
	endif
	if (strlen(units[2]))
		keyVals = ReplaceStringByKey("Zunit",keyVals,units[2],"=")
	endif
	if (WaveExists(hkl))
		keyVals = ReplaceStringByKey("hklPeakCenter",keyVals,vec2str(hkl,sep=","),"=")
	endif
	return keyVals
End
//
Static Function Gaussian3DFitFunc(w,xx,yy,zz) : FitFunc
	Wave w
	Variable xx
	Variable yy
	Variable zz

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ Variable ss = ((xx-x0)/xFWHM)^2 + ((yy-y0)/yFWHM)^2 + ((zz-z0)/zFWHM)^2
	//CurveFitDialog/ ss *= 4*ln(2)
	//CurveFitDialog/ ss = max(-500,ss)
	//CurveFitDialog/ f(xx,yy,zz) = offset + amp*exp(-ss)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 3
	//CurveFitDialog/ xx
	//CurveFitDialog/ yy
	//CurveFitDialog/ zz
	//CurveFitDialog/ Coefficients 8
	//CurveFitDialog/ w[0] = offset
	//CurveFitDialog/ w[1] = amp
	//CurveFitDialog/ w[2] = x0
	//CurveFitDialog/ w[3] = xFWHM
	//CurveFitDialog/ w[4] = y0
	//CurveFitDialog/ w[5] = yFWHM
	//CurveFitDialog/ w[6] = z0
	//CurveFitDialog/ w[7] = zFWHM

	Variable ss = ((xx-w[2])/w[3])^2 + ((yy-w[4])/w[5])^2 + ((zz-w[6])/w[7])^2
	ss *= 4*ln(2)
	ss = max(-500,ss)
	return w[0] + w[1]*exp(-ss)
End


Static Function GaussianCross3DFitFunc(w,xx,yy,zz) : FitFunc
	Wave w
	Variable xx
	Variable yy
	Variable zz

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ x0 -= w[2]  ;  y0 -= w[3]  ;  z0 -= w[4]
	//CurveFitDialog/ Variable ss = (x0/xFWHM)^2 + (y0/yFWHM)^2 + (z0/zFWHM)^2 + x0*y0/Cxy + x0*z0/Cxz + y0*z0/Cyz
	//CurveFitDialog/ f(x0,y0,z0) = bkg + amp * exp(-ss)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 3
	//CurveFitDialog/ x0
	//CurveFitDialog/ y0
	//CurveFitDialog/ z0
	//CurveFitDialog/ Coefficients 11
	//CurveFitDialog/ w[0] = bkg
	//CurveFitDialog/ w[1] = amp
	//CurveFitDialog/ w[2] = x0
	//CurveFitDialog/ w[3] = xFWHM
	//CurveFitDialog/ w[4] = y0
	//CurveFitDialog/ w[5] = yFWHM
	//CurveFitDialog/ w[6] = z0
	//CurveFitDialog/ w[7] = zFWHM
	//CurveFitDialog/ w[8] = Cxy
	//CurveFitDialog/ w[9] = Cxz
	//CurveFitDialog/ w[10] = Cyz

	xx -= w[2]
	yy -= w[3]
	zz -= w[4]
	Variable ss = (xx/w[3])^2 + (yy/w[5])^2 + (zz/w[7])^2 + xx*yy/w[8] + xx*zz/w[9] + yy*zz/w[10]
	ss *= 4*ln(2)
	ss = max(-500,ss)
	return w[0] + w[1]*exp(-ss)
End


Function/WAVE centerOf3Ddata(ww3D)	// finds center of data, works for triplets and for a 3D array
	Wave ww3D
	if (!WaveExists(ww3D))
		return $""
	endif

	if (WaveDims(ww3D)==2 && DimSize(ww3D,1)==3 && DimSize(ww3D,0)>0)	// triplets
		MatrixOP/FREE center = sumCols(ww3D)/numRows(ww3D)
		Redimension/N=3 center
	elseif (WaveDims(ww3D)==3 && numpnts(ww3D)>0)						// 3D array
		Make/N=3/D/FREE ccc
		ccc = DimOffset(ww3D,p) + DimDelta(ww3D,p)*DimSize(ww3D,p)/2
		Wave center = ccc
	else
		WAVE center = $""
	endif
	return center
End


Function/WAVE CenterOfMass3D(w3D)	// returns a 3-vector with the center of mass
	Wave w3D									// a 3D wave (with scaling)

	WaveStats/M=1/Q w3D					// using WaveStats instead of WaveSum() avoids problems with NaN's
	Variable mass = V_sum
	Make/N=3/D/FREE com

	Duplicate/FREE w3D, oneComp
	Redimension/D oneComp

	oneComp *= x
	WaveStats/M=1/Q oneComp
	com[0] = V_sum / mass				// x-component of com

	oneComp = w3D
	oneComp *= y
	WaveStats/M=1/Q oneComp
	com[1] = V_sum / mass				// y-component of com

	oneComp = w3D
	oneComp *= z
	WaveStats/M=1/Q oneComp
	com[2] = V_sum / mass				// z-component of com

	return com
End


Function/WAVE ExtractSubVolume(volumeAll,Vc,HWx,[HWy,HWz])
	Wave volumeAll
	Wave Vc						// point at center of sub-volume (3 vector, scaled coordinates)
	Variable HWx				// Half Width of sub-volume in units of volumeAll (scaled coordinates)
	Variable HWy,HWz
	HWy = ParamIsDefault(HWy) || numtype(HWy) || HWy<=0 ? HWx : HWy
	HWz = ParamIsDefault(HWz) || numtype(HWz) || HWz<=0 ? HWx : HWz

	Make/N=3/D/FREE HW = {HWx,HWy,HWz}
	if (WaveDims(volumeAll)!=3 || numpnts(volumeAll)<2 || numtype(sum(HW)+sum(Vc)) || numpnts(Vc)!=3)
		return $""
	endif
	Make/N=3/D/FREE dxyz=DimDelta(volumeAll,p)
	Make/N=3/I/FREE ijk = round( (Vc[p] - DimOffset(volumeAll,p)) / dxyz[p] )
	Make/N=3/I/FREE Nxyz=DimSize(volumeAll,p)
	Make/N=3/I/FREE ijkLo0, ijkHi0, ijkLo,ijkHi
	ijkLo0 = limit(round(ijk - HW[p]/dxyz[p]), 0, Nxyz[p]-1)		// ijkLo0 goes with first x
	ijkHi0 = limit(round(ijk + HW[p]/dxyz[p]), 0, Nxyz[p]-1)		// ijkHi0 goes with last x
	Make/N=3/D/FREE xyz0 = ijkLo0[p]*dxyz[p] + DimOffset(volumeAll,p)
	ijkLo = min(ijkLo0[p], ijkHi0[p])			// changes when dx is negative
	ijkHi = max(ijkLo0[p], ijkHi0[p])			// now have ijkLo[p] < ijkHi[p]
	Make/N=3/I/FREE Nsub = ijkHi[p] - ijkLo[p] + 1

	if (WaveMin(Nsub)<1)
		return $""
	endif
	Make/N=(Nsub[0],Nsub[1],Nsub[2])/FREE/D subV=0
	SetScale/P x xyz0[0],dxyz[0],WaveUnits(volumeAll,0), subV
	SetScale/P y xyz0[1],dxyz[1],WaveUnits(volumeAll,0), subV
	SetScale/P z xyz0[2],dxyz[2],WaveUnits(volumeAll,0), subV
	subV = volumeAll[p+ijkLo[0]][q+ijkLo[1]][r+ijkLo[2]]

	String str, wnote=note(volumeAll)
	sprintf str, "%d,%d,%d,%d,%d,%d", ijkLo[0],ijkHi[0],ijkLo[1],ijkHi[1],ijkLo[2],ijkHi[2]
	wnote = ReplaceStringByKey("roi3D",wnote,str,"=")
	Note/K subV, wnote
	return subV
End


Function isPointInVolume(vol,vec)	// verify if point is in a volume, works for triplets and for a 3D array
	Wave vol								// represents the volume
	Wave vec								// must be a 3-vec

	if (numpnts(vec)!=3)
		return 0
	endif

	if (WaveDims(vol)==2 && DimSize(vol,1)==3 && DimSize(vol,0)>0)	// triplets
		MatrixOP/FREE lo = minVal(col(vol,0))
		MatrixOP/FREE hi = maxVal(col(vol,0))
		if (vec[0]<lo[0] || vec[0]>hi[0])
			return 0
		endif

		MatrixOP/FREE lo = minVal(col(vol,1))
		MatrixOP/FREE hi = maxVal(col(vol,1))
		if (vec[1]<lo[0] || vec[1]>hi[0])
			return 0
		endif

		MatrixOP/FREE lo = minVal(col(vol,2))
		MatrixOP/FREE hi = maxVal(col(vol,2))
		if (vec[2]<lo[0] || vec[2]>hi[0])
			return 0
		endif

	elseif (WaveDims(vol)==3 && numpnts(vol)>0)						// 3D array
		Make/N=3/D/FREE ijk = round( (vec[p]-DimOffset(vol,p))/DimDelta(vol,p) )
		if (WaveMin(ijk)<0)
			return 0
		endif
		Make/N=3/I/FREE tooHigh = ijk[p] > (DimSize(vol,p)-1)
		if (WaveMax(tooHigh))
			return 0
		endif
	endif

	return 1
End



//  ======================================================================================  //
//  ========================= Start of Bounding Volume Structure =========================  //

Structure boundingVolume			// a generic bounding volume, probably in k-space
	double	xlo, xhi					// range of x, these 6 form the corners of a box
	double	ylo, yhi					// range of y
	double	zlo, zhi					// range of z
	double	xW, yW, zW				// box size is xW x yW x zW
	int16		Nx, Ny, Nz				// dimensions of the array
	double	dx, dy, dz				// step size along box
	double	vol						// volume of box = xW*yW*zW
EndStructure
//
//ThreadSafe Function initBoundingVolumeStruct(v)
Function initBoundingVolumeStruct(v)
	STRUCT boundingVolume &v
	v.xlo = Inf  ;	v.ylo = Inf ;	v.zlo = Inf
	v.xhi = -Inf ;	v.yhi = -Inf ;	v.zhi = -Inf
	v.Nx = 0     ;	v.Ny = 0 ;		v.Nz = 0			// default dimension
	updateBoundingVolumeStruct(v)
End
//
//ThreadSafe Function updateBoundingVolumeStruct(v)
Function updateBoundingVolumeStruct(v)
	STRUCT boundingVolume &v
	v.xW = abs(v.xhi - v.xlo)
	v.yW = abs(v.yhi - v.ylo)
	v.zW = abs(v.zhi - v.zlo)
	v.vol = v.xW * v.yW * v.zW

	if (v.Nx > 1)					// prefer to use Nx
		v.dx = v.xW / (v.Nx - 1)
	elseif (v.dx > 0)				// N not good, try to use dx
		v.Nx = ceil(v.xW / v.dx) + 1
	endif

	if (v.Ny > 1)					// prefer to use Ny
		v.dy = v.yW / (v.Ny - 1)
	elseif (v.dy > 0)				// N not good, try to use dy
		v.Ny = ceil(v.yW / v.dy) + 1
	endif

	if (v.Nz > 1)					// prefer to use Nz
		v.dz = v.zW / (v.Nz - 1)
	elseif (v.dz > 0)				// N not good, try to use dz
		v.Nz = ceil(v.zW / v.dz) + 1
	endif

	v.dx = v.Nx>1 ? v.xW / (v.Nx - 1) : 0		// re-update the dx
	v.dy = v.Ny>1 ? v.yW / (v.Ny - 1) : 0
	v.dz = v.Nz>1 ? v.zW / (v.Nz - 1) : 0
End


//ThreadSafe Function/S boundingVolumeStruct2str(v)
Function/S boundingVolumeStruct2str(v)
	STRUCT boundingVolume &v
	String str
//	sprintf str,"Vol = X=[%g, %g] ÆX=%g,  Y=[%g, %g] ÆY=%g,  Z=[%g, %g] ÆZ=%g,  Vol=%g\r",v.xlo,v.xhi,v.xW, v.ylo,v.yhi,v.yW, v.zlo,v.zhi,v.zW,v.vol
	sprintf str,"Vol = X=[%g, %g] Nx=%g,  Y=[%g, %g] Ny=%g,  Z=[%g, %g] Nz=%g,  Vol=%g\r",v.xlo,v.xhi,v.Nx, v.ylo,v.yhi,v.Ny, v.zlo,v.zhi,v.Nz,v.vol
	return str
End


//ThreadSafe Function/S boundingVolumeStructEncode(v)
Function/S boundingVolumeStructEncode(v)		// convert v --> str
	STRUCT boundingVolume &v

	String str
	String fmt = "xlo:%.15g,xhi:%.15g,ylo:%.15g,yhi:%.15g,zlo:%.15g,zhi:%.15g,xW:%.15g,yW:%.15g,zW:%.15g,Nx:%.15g,Ny:%.15g,Nz:%.15g,dx:%.15g,dy:%.15g,dz:%.15g,vol:%.15g"
	sprintf str,fmt, v.xlo,v.xhi,v.ylo,v.yhi,v.zlo,v.zhi,v.xW,v.yW,v.zW,v.Nx,v.Ny,v.Nz,v.dx,v.dy,v.dz,v.vol
	return str
End


//ThreadSafe Function boundingVolumeStructDecode(str,v)
Function boundingVolumeStructDecode(str,v)		// convert str --> v
	STRUCT boundingVolume &v
	String str

	String fmt = "xlo:%g,xhi:%g,ylo:%g,yhi:%g,zlo:%g,zhi:%g,xW:%g,yW:%g,zW:%g,Nx:%g,Ny:%g,Nz:%g,dx:%g,dy:%g,dz:%g,vol:%g"
	Variable xlo,xhi,ylo,yhi,zlo,zhi,xW,yW,zW,Nx,Ny,Nz,dx,dy,dz,vol
	sscanf str, fmt, xlo,xhi,ylo,yhi,zlo,zhi,xW,yW,zW,Nx,Ny,Nz,dx,dy,dz,vol
	Variable i = V_flag
	if (i!=16)
		initBoundingVolumeStruct(v)
		return 1
	endif
	v.xlo = xlo	;	v.xhi = xhi	;	v.xW = xW	;	v.Nx = Nx	;	v.dx = dx
	v.ylo = ylo	;	v.yhi = yhi	;	v.yW = yW	;	v.Ny = Ny	;	v.dy = dy
	v.zlo = zlo	;	v.zhi = zhi	;	v.zW = zW	;	v.Nz = Nz	;	v.dz = dz
	v.vol = vol
	return 0
End


//ThreadSafe Function copyBoundingVolumeStructs(vin,vout)		// note, vall may be v1 or v2
Function copyBoundingVolumeStructs(vin,vout)		// note, vall may be v1 or v2
	STRUCT boundingVolume &vin, &vout
	vout.xlo = vin.xlo ;	vout.ylo = vin.ylo ;	vout.zlo = vin.zlo
	vout.xhi = vin.xhi ;	vout.yhi = vin.yhi ;	vout.zhi = vin.zhi
	vout.xW = vin.xW   ;	vout.yW = vin.yW   ;	vout.zW = vin.zW
	vout.Nx = vin.Nx   ;	vout.Ny = vin.Ny   ;	vout.Nz = vin.Nz
	vout.dx = vin.dx   ;	vout.dy = vin.dy   ;	vout.dz = vin.dz
	vout.vol = vin.vol
End
//
//ThreadSafe Function extendBoundingVolumeStruct(v,vec)
Function extendBoundingVolumeStruct(v,vec)
	// extend v to include the vector vec, maintain the dx,dy,dz when doing this
	STRUCT boundingVolume &v
	Wave vec

	updateBoundingVolumeStruct(v)									// ensure that dx,dy,dz are calculated
	Variable xlo = min(xlo,vec[0]), xhi = max(xhi,vec[0])	// new ranges
	Variable ylo = min(ylo,vec[1]), yhi = max(yhi,vec[1])
	Variable zlo = min(zlo,vec[2]), zhi = max(zhi,vec[2])
	Variable Nx,Ny,Nz
	Nx = ceil(abs(xhi-xlo) / v.dx) + 1							// new Nx,Ny,Nz
	Ny = ceil(abs(yhi-ylo) / v.dy) + 1
	Nz = ceil(abs(zhi-zlo) / v.dz) + 1
	v.Nx = Nx>0 && numtype(Nx)==0 ? Nx : 0
	v.Ny = Ny>0 && numtype(Ny)==0 ? Ny : 0
	v.Nz = Nz>0 && numtype(Nz)==0 ? Nz : 0

	v.xlo = xlo ;		v.xhi = xhi
	v.ylo = ylo ;		v.yhi = yhi
	v.zlo = zlo ;		v.zhi = zhi
	updateBoundingVolumeStruct(v)
End


//ThreadSafe Function insideBoundingVolumeStruct(v,vec)		// returns TRUE if vec is inside v
Function insideBoundingVolumeStruct(v,vec)		// returns TRUE if vec is inside v
	STRUCT boundingVolume &v
	Wave vec		// either a simple vec[3], or vec[N][3], (actually it only uses the first 3 values, so 3 could be 4)

	Variable dim=WaveDims(vec), inside=0
	if (dim==1)
		if (DimSize(vec,0)<3)
			return 0
		endif
		inside = v.xlo <= vec[0] && vec[0] <= v.xhi
		inside = inside && (v.ylo <= vec[1] && vec[1] <= v.yhi)
		inside = inside && (v.zlo <= vec[2] && vec[2] <= v.zhi)
	elseif (dim==2)						// list of vectors
		if (DimSize(vec,1)<3)
			return 0
		endif
		Variable lo = v.xlo, hi = v.xhi
		MatrixOP/FREE outside0 = maxVal( greater(lo,col(vec,0)) || greater(col(vec,0),hi) )

		lo = v.ylo	;	hi = v.yhi
		MatrixOP/FREE outside1 = maxVal( greater(lo,col(vec,1)) || greater(col(vec,1),hi ) )

		lo = v.zlo	;	hi = v.zhi
		MatrixOP/FREE outside2 = maxVal( greater(lo,col(vec,2)) || greater(col(vec,2),hi) )
		inside = !(outside0[0] || outside1[0] || outside2[0])
	else
		inside = 0							// vec cannot be 3D
	endif
	return inside
End


//ThreadSafe Function combineBoundingVolumeStructs(v1,v2,vall)	// note, vall may be v1 or v2
Function combineBoundingVolumeStructs(v1,v2,vall)	// note, vall may be v1 or v2
	STRUCT boundingVolume &v1, &v2, &vall

	updateBoundingVolumeStruct(v1)
	updateBoundingVolumeStruct(v2)

	Variable xlo=min(v1.xlo, v2.xlo), xhi=max(v1.xhi, v2.xhi)
	Variable ylo=min(v1.ylo, v2.ylo), yhi=max(v1.yhi, v2.yhi)
	Variable zlo=min(v1.zlo, v2.zlo), zhi=max(v1.zhi, v2.zhi)
	Variable dx=max(v1.dx, v2.dx), dy=max(v1.dy, v2.dy), dz=max(v1.dz, v2.dz)
	Variable Nx,Ny,Nz
	Nx = ceil(abs(xhi-xlo) / dx) + 1						// new Nx,Ny,Nz
	Ny = ceil(abs(yhi-ylo) / dy) + 1
	Nz = ceil(abs(zhi-zlo) / dz) + 1
	vall.Nx = Nx>0 && numtype(Nx)==0 ? Nx : 0
	vall.Ny = Ny>0 && numtype(Ny)==0 ? Ny : 0
	vall.Nz = Nz>0 && numtype(Nz)==0 ? Nz : 0
	vall.xlo = xlo ;		vall.xhi = xhi
	vall.ylo = ylo ;		vall.yhi = yhi
	vall.zlo = zlo ;		vall.zhi = zhi
	updateBoundingVolumeStruct(vall)
End


// set a boundingVolume structure from the dimensions of a 3D wave
Function Fill_boundingVolume_from_3Dwave(v,w3D)
	STRUCT boundingVolume &v
	Wave w3D

	Variable lo,hi
	lo = DimOffset(w3D,0)
	hi = lo + DimDelta(w3D,0)*(DimSize(w3D,0)-1)
	OrderValues(lo,hi)	
	v.xlo = lo
	v.xhi = hi

	lo = DimOffset(w3D,1)
	hi = lo + DimDelta(w3D,1)*(DimSize(w3D,1)-1)
	OrderValues(lo,hi)	
	v.ylo = lo
	v.yhi = hi

	lo = DimOffset(w3D,2)
	hi = lo + DimDelta(w3D,2)*(DimSize(w3D,2)-1)
	OrderValues(lo,hi)	
	v.zlo = lo
	v.zhi = hi

	v.Nx = DimSize(w3D,0)
	v.Ny = DimSize(w3D,1)
	v.Nz = DimSize(w3D,2)
	updateBoundingVolumeStruct(v)
End

//  ========================== End of Bounding Volume Structure ==========================  //
//  ======================================================================================  //



//  ======================================================================================  //
//  =================================== Start of Init ====================================  //

//Static Function InitStats3D()
//End

//  ==================================== End of Init =====================================  //
//  ======================================================================================  //



