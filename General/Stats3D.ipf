#pragma rtGlobals=3		// Use modern global access method.
#pragma version = 0.01
#pragma IgorVersion = 6.3
#pragma ModuleName=Stats3D


// need to work on ExtractSubVolume, the whole swap thing is not right


Static Function/WAVE ExtractSubVolume(volumeAll,Vc,HWx,[HWy,HWz])
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
	Make/N=3/I/FREE ijk = round( (Vc[p] - DimOffset(volumeAll,p)) / DimDelta(volumeAll,p) )
	Make/N=3/I/FREE Nxyz=DimSize(volumeAll,p)
	Make/N=3/I/FREE ijkLo, ijkHi
	ijkLo = limit(round(ijk - HW[p]/DimDelta(volumeAll,p)), 0, Nxyz[p]-1)
	ijkHi = limit(round(ijk + HW[p]/DimDelta(volumeAll,p)), 0, Nxyz[p]-1)
	Variable swap
//	if (ijkLo[0]>ijkHi[0])			// ensure proper order
//		swap = ijkLo[0]
//		ijkLo[0] = ijkHi[0]
//		ijkHi[0] = swap
//	endif
//	if (ijkLo[1]>ijkHi[1])
//		swap = ijkLo[1]
//		ijkLo[1] = ijkHi[1]
//		ijkHi[1] = swap
//	endif
//	if (ijkLo[2]>ijkHi[2])
//		swap = ijkLo[2]
//		ijkLo[2] = ijkHi[2]
//		ijkHi[2] = swap
//	endif
//	ijkLo = max(ijkLo[p],0)
//	ijkHi = min(ijkHi[p], Nxyz[p]-1)
	Make/N=3/D/FREE xyz0 = ijkLo[p]*DimDelta(volumeAll,p) + DimOffset(volumeAll,p)
	Make/N=3/I/FREE Nsub = abs(ijkHi[p] - ijkLo[p] + 1)

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

Function/T PeakMax3D(vol, [printIt])	// finds max of data in a 3D volume
	Wave vol
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRtStackInfo(2))==0 : printIt
	if (!(WaveDims(vol)==3 && numpnts(vol)>0))						// 3D volume array
		return ""
	endif

	WaveStats/M=1/Q vol
	Variable maxVal = V_max
	Make/D/FREE xyzMax={V_maxRowLoc,V_maxColLoc, V_maxLayerLoc}
	Make/N=(3)/D/FREE sizes = (DimSize(vol,p)-1) * DimDelta(vol,p)
	Make/N=3/I/FREE ijk = round( (xyzMax[p] - DimOffset(vol,p)) / DimDelta(vol,p) )
	Variable npnt = ijk[2]*DimSize(vol,1)*DimSize(vol,0) + ijk[1]*DimSize(vol,0) + ijk[0]
	Variable HW = WaveMin(sizes)/2
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


Function/T FitPeakIn3D(space3D,startXYZ,HWx,[HWy,HWz,printIt])
	Wave space3D
	Wave startXYZ						// starting xyz for the fit
	Variable HWx,HWy,HWz			// half widths dx, dy, dz for the sub volume
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

	Make/N=3/D/FREE HW={HWx,HWy,HWz}			// half widths in x, y, & z for the sub volume
	Duplicate/FREE space3D, stdDev	// this assumes that one count in detector == 1 photon

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
	FuncFitMD/Q GMarkers#Gaussian3DFitFunc, W_coef, space3D/W=stdDev/I=1
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



Function/T XX_FitPeakAt3Dmarker(space3D,Qc,QxHW,[QyHW,QzHW,printIt])
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
	FuncFitMD/Q GMarkers#Gaussian3DFitFunc, W_coef, sub3D/W=stdDev/I=1
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
//  =================================== Start of Init ====================================  //

//Static Function InitStats3D()
//End

//  ==================================== End of Init =====================================  //
//  ======================================================================================  //



