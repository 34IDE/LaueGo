#pragma rtGlobals=3		// Use modern global access method.
#pragma version = 0.09
#pragma IgorVersion = 6.3
#pragma ModuleName=Stats3D




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


Function FitPeakIn3D(space3D,GP, HWx,[HWy,HWz, startXYZ, stdDev, func3D, printIt])
	Wave space3D
	STRUCT Generic3DPeakStructure &GP	// holds result of fitting
	Variable HWx,HWy,HWz			// starting half widths dx, dy, dz for the fit
	Wave startXYZ						// starting xyz for the fit (defaults to position of peak)
	Wave stdDev							// errors in space3D, standard deviation of each value in space3D
	String func3D						// name of 3D peak fitting function, defaults to "Gaussian3DFitFunc"
	Variable printIt
	func3D = SelectString(strlen(func3D), "Gaussian3DFitFunc", func3D)
	printIt = ParamIsDefault(printIt) || numtype(printIt)? strlen(GetRTStackInfo(2))==0 : !(!printIt)
	HWy = ParamIsDefault(HWy) || numtype(HWy) || HWy<=0 ? HWx : HWy	// HWy & HWz default to HWx
	HWz = ParamIsDefault(HWz) || numtype(HWz) || HWz<=0 ? HWx : HWz

	initGeneric3DPeakStructure(GP)
	FUNCREF Peak3DFitFuncProto FitFunc = $func3D
	func3D = StringByKey("NAME",FuncRefInfo(FitFunc))
	func3D = SelectString(strlen(func3D), "Gaussian3DFitFunc", func3D)

	if (!WaveExists(space3D))		// is space3D valid?
		return 1
	elseif (WaveDims(space3D)!=3)
		return 1
	elseif (strlen(func3D)<1)
		return 1
	endif

	if (!WaveExists(startXYZ))
		WaveStats/M=1/Q space3D	// set starting point to the position of max value
		Make/D/FREE xyz0 = {V_maxRowLoc,V_maxColLoc, V_maxLayerLoc}
	else
		Wave xyz0 = startXYZ
	endif
	if (!WaveExists(xyz0))
		return 1
	elseif (numpnts(xyz0)!=3 || numtype(sum(xyz0)))
		return 1
	elseif (numtype(HWx+HWy+HWz) || HWx<=0 || HWy<=0 || HWz<=0)
		return 1
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
	Make/N=8/D/O W_coef = {NaN,NaN,xyz0[0],HW[0]/4,xyz0[1],HW[1]/4,xyz0[2],HW[2]/4}
	if (maxVal>-minVal)			// a positive peak
		W_coef[1] = maxVal
		W_coef[0] = minVal==0 ? 1 : minVal
	else								// a negative peak
		W_coef[1] = minVal
		W_coef[0] = maxVal==0 ? 1 : maxVal
	endif

	Variable V_FitOptions = (printIt ? 0 : 4)
	Variable V_FitError=0, V_FitQuitReason=0
	FuncFitMD/Q FitFunc, W_coef, space3D/W=errWave/I=1
	Variable chisq = V_chisq
	Wave W_sigma=W_sigma
	W_coef[3] = abs(W_coef[3])		// FWx, FWy, FWz, must be >=0
	W_coef[5] = abs(W_coef[5])
	W_coef[7] = abs(W_coef[7])
	Make/N=3/T/FREE units=WaveUnits(space3D,p)
	Make/N=3/D/FREE xyzFit=W_coef[2*p + 2]

	Variable err=V_FitError, swap
	err = err || !isPointInVolume(space3D,xyzFit)	// verify if point is in a volume, works for triplets and for a 3D array
	if (err)
		print "ERROR -- Fit failed, peak not inside of fitting volume"
		print FitErrorString(V_FitError,V_FitQuitReason)
		return 1
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

	Make/N=3/I/FREE ijk, Nxyz=DimSize(space3D,p)
	ijk = round( (xyzFit[p]-DimOffset(space3D,p)) / DimDelta(space3D,p) )
	ijk = limit(round(ijk[p]),0,Nxyz[p]-1)

	GP.OK = 1
	GP.bkg = W_coef[0]		;		GP.bkgErr = W_sigma[0]
	GP.amp = W_coef[1]		;		GP.ampErr = W_sigma[1]
	GP.x = W_coef[2]			;		GP.xErr = W_sigma[2]
	GP.y = W_coef[4]			;		GP.yErr = W_sigma[4]
	GP.z = W_coef[6]			;		GP.zErr = W_sigma[6]
	GP.FWx = W_coef[3]		;		GP.FWxErr = W_sigma[3]
	GP.FWy = W_coef[5]		;		GP.FWyErr = W_sigma[5]
	GP.FWz = W_coef[7]		;		GP.FWzErr = W_sigma[7]
	GP.chisq = V_chisq

	GP.ix = ijk[0]				;		GP.iy = ijk[1]				;	GP.iz = ijk[2]
	GP.maxValue = space3D[ijk[0]][ijk[1]][ijk[2]]

	GP.Xunit = SelectString(strlen(units[0]),"", units[0])
	GP.Yunit = SelectString(strlen(units[1]),"", units[1])
	GP.Zunit = SelectString(strlen(units[2]),"", units[2])
	GP.wname = SelectString(WaveType(space3D,2)==1, "", NameOfWave(space3D))
	GP.funcName = func3D

	if (numpnts(W_coef)>10)
		GP.Cxy = W_coef[8]		;		GP.CxyErr = W_sigma[8]
		GP.Cxz = W_coef[9]		;		GP.CxzErr = W_sigma[9]
		GP.Cyz = W_coef[10]		;		GP.CyzErr = W_sigma[10]
	endif

	if (WaveExists(hkl))
		GP.hkl[0] = hkl[0]
		GP.hkl[1] = hkl[1]
		GP.hkl[2] = hkl[2]
	endif

	if (printIt)
		printGeneric3DPeakStructure(GP)
	endif

	return 0				// no error
End


Function Peak3DFitFuncProto(w,xx,yy,zz) : FitFunc
	Wave w
	Variable xx
	Variable yy
	Variable zz
	return Gaussian3DFitFunc(w,xx,yy,zz)
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
//  ========================== Start of 3D Peak Shape Structure ==========================  //

Static Constant PeakShapeUnitMaxLen = 50
Structure Generic3DPeakStructure	// describes the Generic fit to a 3D peak
	int32		OK								// values are OK, 1==OK, 0==BAD
	double	bkg, bkgErr					// background and error
	double	amp, ampErr					// amplitude and error
	double	x, xErr						// x-center and error
	double	y, yErr						// y-center and error
	double	z, zErr						// z-center and error
	double	FWx, FWxErr					// FWHM in x-direction and error
	double	FWy, FWyErr					// FWHM in y-direction and error
	double	FWz, FWzErr					// FWHM in z-direction and error
	double	chisq							// chisq for the fit
	double	Cxy, CxyErr					// OPTIONAL, xy cross term and error, use NaN if cross-term not fit
	double	Cxz, CxzErr					// OPTIONAL, xz cross term and error
	double	Cyz, CyzErr					// OPTIONAL, yz cross term and error
	char		Xunit[PeakShapeUnitMaxLen]	// OPTIONAL, name of x-units
	char		Yunit[PeakShapeUnitMaxLen]	// OPTIONAL, name of y-units
	char		Zunit[PeakShapeUnitMaxLen]	// OPTIONAL, name of z-units
	double	hkl[3]						// OPTIONAL, contains possible hkl
	int32		ix,iy,iz						// closest points in 3d volume to {x,y,z}
	double	maxValue						// value at [ix,iy,iz]
	char		wname[50]					// OPTIONAL, name of wave that was fit
	char		funcName[50]				// OPTIONAL, name of function used to fit
EndStructure


//ThreadSafe Function updateGeneric3DPeakStructure(GP)
//	STRUCT Generic3DPeakStructure &GP
//	Variable bad = numtype(GP.bkg + GP.amp)
//	bad = bad || numtype(GP.x + GP.y + GP.z)
//	bad = bad || numtype(GP.FWx + GP.FWy + GP.FWz)
//	GP.OK = bad ? 0 : 1
//End


ThreadSafe Function initGeneric3DPeakStructure(GP)
	STRUCT Generic3DPeakStructure &GP
	GP.OK = 0							// BAD values
	GP.bkg = NaN	;	GP.bkgErr = NaN
	GP.amp = NaN	;	GP.ampErr = NaN
	GP.x = NaN		;	GP.xErr = NaN
	GP.y = NaN		;	GP.yErr = NaN
	GP.z = NaN		;	GP.zErr = NaN
	GP.FWx = NaN	;	GP.FWxErr = NaN
	GP.FWy = NaN	;	GP.FWyErr = NaN
	GP.FWz = NaN	;	GP.FWzErr = NaN
	GP.chisq = NaN
	GP.Cxy = NaN	;	GP.CxyErr = NaN
	GP.Cxz = NaN	;	GP.CxzErr = NaN
	GP.Cyz = NaN	;	GP.CyzErr = NaN
	GP.ix = 0		;	GP.iy = 0		; GP.iz = 0
	GP.maxValue = NaN
	GP.wname = ""
	GP.funcName = ""
	GP.Xunit = ""	;	GP.Yunit = ""	;	GP.Zunit = ""
	GP.hkl[0] = NaN;	GP.hkl[1] = NaN;	GP.hkl[2] = NaN
End


ThreadSafe Function copyGeneric3DPeakStructure(in,dest)
	STRUCT Generic3DPeakStructure &in, &dest

	dest.OK		= in.OK
	dest.bkg		= in.bkg		;	dest.bkgErr	= in.bkgErr
	dest.amp		= in.amp		;	dest.ampErr	= in.ampErr
	dest.x		= in.x		;	dest.xErr	= in.xErr
	dest.y		= in.y		;	dest.yErr	= in.yErr
	dest.z		= in.z		;	dest.zErr	= in.zErr
	dest.FWx		= in.FWx		;	dest.FWxErr	= in.FWxErr
	dest.FWy		= in.FWy		;	dest.FWyErr	= in.FWyErr
	dest.FWz		= in.FWz		;	dest.FWzErr	= in.FWzErr
	dest.chisq	= in.chisq
	dest.Cxy		= in.Cxy		;	dest.CxyErr	= in.CxyErr
	dest.Cxz		= in.Cxz		;	dest.CxzErr	= in.CxzErr
	dest.Cyz		= in.Cyz		;	dest.CyzErr	= in.CyzErr
	dest.ix		= in.ix		;	dest.iy		= in.iy	;	dest.iz = in.iz
	dest.maxValue = in.maxValue
	dest.funcName	= in.funcName
	dest.wname	= in.wname

	dest.Xunit = in.Xunit
	dest.Yunit = in.Yunit
	dest.Zunit = in.Zunit
	dest.hkl[0] = in.hkl[0];	dest.hkl[1] = in.hkl[1]	;	dest.hkl[2] = in.hkl[2]
End


Function printGeneric3DPeakStructure(GP)
	STRUCT Generic3DPeakStructure &GP

	if (!(GP.OK))
		print "Generic3DPeakStructure is INVALID"
		return 1
	endif

	Make/N=3/T/FREE units = {GP.Xunit, GP.Yunit, GP.Zunit}
	units = SelectString(stringmatch(units[p],"nm\\S-1*"),units[p],"1/nm")		// two different ways of writing 1/nm
	units = SelectString(strlen(units[p]),""," ("+units[p]+")")						// add () when something present
	Make/N=3/T/FREE Qstr=SelectString(strsearch(units[p],"1/nm",0)>=0,"","Q")	// is it Q?

	String funcName = SelectString(strlen(GP.funcName), "Generic", GP.funcName)
	printf "%s Fit has offset = %s,  amp = %s,  chisq = %g\r",funcName,ValErrStr(GP.bkg,GP.bkgErr,sp=1),ValErrStr(GP.amp,GP.ampErr,sp=1), GP.chisq
	printf "  <%sx> = %s%s,   FWHMx = %s%s\r",Qstr[0],ValErrStr(GP.x,GP.xErr,sp=1),units[0],ValErrStr(GP.FWx,GP.FWxErr,sp=1),units[0]
	printf "  <%sy> = %s%s,   FWHMy = %s%s\r",Qstr[0],ValErrStr(GP.y,GP.yErr,sp=1),units[0],ValErrStr(GP.FWy,GP.FWyErr,sp=1),units[0]
	printf "  <%sz> = %s%s,   FWHMz = %s%s\r",Qstr[0],ValErrStr(GP.z,GP.zErr,sp=1),units[0],ValErrStr(GP.FWz,GP.FWzErr,sp=1),units[0]

	if (stringmatch(units[0],units[1]) && stringmatch(units[0],units[2]))
		Variable xyzMag=sqrt((GP.x)^2 + (GP.y)^2 + (GP.z)^2)
		Variable dxyzMag=sqrt((GP.xErr)^2 + (GP.yErr)^2 + (GP.zErr)^2)
		printf "  <|%s|> = %s%s\r",Qstr[0],ValErrStr(xyzMag,dxyzMag,sp=1),units[0]
	endif
	printf "  Closest point to peak is the %s[%g, %g, %g] = %g\r",GP.wname,GP.ix,GP.iy,GP.iz,GP.maxValue

	Make/N=3/D/FREE hkl = GP.hkl[p]
	if (numtype(sum(hkl))==0)
		printf "hkl Peak Center = %s\r",vec2str(hkl)
	endif
	return 1
End

//  =========================== End of 3D Peak Shape Structure ===========================  //
//  ======================================================================================  //



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


//	*************************************************************************************
//	***** DEPRECATED ********* DEPRECATED ********* DEPRECATED ********* DEPRECATED *****
//	****************************** PeakMax3D is DEPRECATED ******************************
//		PeakMax3D is deprecated, use FitPeakIn3D() directly instead.
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

	STRUCT Generic3DPeakStructure GP	// holds result of fitting
	Variable err = FitPeakIn3D(vol,GP,HW,startXYZ=xyzMax,printIt=printIt)
	if (err)
		return ""
	endif

	String keyVals=""
	keyVals = ReplaceNumberByKey("offset",keyVals,GP.bkg,"=")
	keyVals = ReplaceNumberByKey("offsetErr",keyVals,GP.bkgErr,"=")
	keyVals = ReplaceNumberByKey("amp",keyVals,GP.amp,"=")
	keyVals = ReplaceNumberByKey("ampErr",keyVals,GP.ampErr,"=")
	keyVals = ReplaceNumberByKey("Xc",keyVals,GP.x,"=")
	keyVals = ReplaceNumberByKey("XcErr",keyVals,GP.xErr,"=")
	keyVals = ReplaceNumberByKey("Yc",keyVals,GP.y,"=")
	keyVals = ReplaceNumberByKey("YcErr",keyVals,GP.yErr,"=")
	keyVals = ReplaceNumberByKey("Zc",keyVals,GP.z,"=")
	keyVals = ReplaceNumberByKey("ZcErr",keyVals,GP.zErr,"=")


	keyVals = ReplaceNumberByKey("FWX",keyVals,GP.FWx,"=")
	keyVals = ReplaceNumberByKey("FWXErr",keyVals,GP.FWxErr,"=")
	keyVals = ReplaceNumberByKey("FWY",keyVals,GP.FWy,"=")
	keyVals = ReplaceNumberByKey("FWYErr",keyVals,GP.FWyErr,"=")
	keyVals = ReplaceNumberByKey("FWZ",keyVals,GP.FWz,"=")
	keyVals = ReplaceNumberByKey("FWZErr",keyVals,GP.FWzErr,"=")
	keyVals = ReplaceNumberByKey("chisq",keyVals,GP.chisq,"=")

	if (numtype(GP.Cxy + GP.Cxz + GP.Cxz)==0)
		keyVals = ReplaceNumberByKey("Cxy",keyVals,GP.Cxy,"=")
		keyVals = ReplaceNumberByKey("CxyErr",keyVals,GP.CxyErr,"=")
		keyVals = ReplaceNumberByKey("Cxz",keyVals,GP.Cxz,"=")
		keyVals = ReplaceNumberByKey("CxzErr",keyVals,GP.CxzErr,"=")
		keyVals = ReplaceNumberByKey("Cyz",keyVals,GP.Cyz,"=")
		keyVals = ReplaceNumberByKey("CyzErr",keyVals,GP.CyzErr,"=")
	endif

	if (strlen(GP.Xunit))
		keyVals = ReplaceStringByKey("Xunit",keyVals,GP.Xunit,"=")
	endif
	if (strlen(GP.Yunit))
		keyVals = ReplaceStringByKey("Yunit",keyVals,GP.Yunit,"=")
	endif
	if (strlen(GP.Zunit))
		keyVals = ReplaceStringByKey("Zunit",keyVals,GP.Zunit,"=")
	endif

	if ( numtype(GP.hkl[0] + GP.hkl[1] + GP.hkl[2])==0 )
		Make/N=3/D/FREE hkl = GP.hkl[p]
		keyVals = ReplaceStringByKey("hklPeakCenter",keyVals,vec2str(hkl,sep=","),"=")
	endif

	keyVals = ReplaceNumberByKey("iStart",keyVals,ijk[0],"=")
	keyVals = ReplaceNumberByKey("jStart",keyVals,ijk[1],"=")
	keyVals = ReplaceNumberByKey("kStart",keyVals,ijk[2],"=")
	return keyVals
End
//	****************************** PeakMax3D is DEPRECATED ******************************
//	***** DEPRECATED ********* DEPRECATED ********* DEPRECATED ********* DEPRECATED *****
//	*************************************************************************************





//  ======================================================================================  //
//  =================================== Start of Init ====================================  //

//Static Function InitStats3D()
//End

//  ==================================== End of Init =====================================  //
//  ======================================================================================  //



