#pragma rtGlobals=1		// Use modern global access method.
#pragma IgorVersion = 6.11
#pragma version = 2.05
#pragma ModuleName=fwhm

// with v 2.0, major revision, started using structures
//		The old functions are still there, but with new innards

Static Constant MAX_PEAK_COEFS = 6


Menu "Analysis"
	"-"
	MenuItemIfWavesExists("Simple Peak Parameters of ...","*","DIMS:1,MINROWS:3"), FindSimplePeakParameters($"",$"")
	MenuItemIfWavesExists("Peak Parameters [W_coef]","W_coef","DIMS:1,MINROWS:4"), FWHM_Coefs()
	MenuItemIfFitPeakOnGraph("Fit Displayed Peak to a Resolution","peakWave;Resolution*","DIMS:1"), FitPeakOnGraph("")
End



//  ======================================================================================  //
//  ============================= Start of Fit a Peak Wave  ==============================  //

Function/S MenuItemIfFitPeakOnGraph(item,classes,optionsStr)
	String item						// string that you want to appear in menu
	String classes					// semi-colon separated list of wave classes, can use "*" in classes
	String optionsStr				// options that are found in WaveList function optionsStr
	Variable invisible			// controls menu item when conditions not met: true -> menu item absent, false or absent -> menu item grayed out
	Variable all					// when all is TRUE, then all of the classes in waveClassList must be present, not just one
	String list = WaveListClass(classes,"*",optionsStr,all=all)
	if (strlen(list) && strlen(WinList("*",";","WIN:1")))
		return item
	endif
	return "("+item
End

// This routine is used when fitting to data on a graph
Function FitPeakOnGraph(gName)
	String gName						// graph name, use "" for top graph
	String pkList=WaveListClass("peakWave;Resolution*","*","DIMS:1")
	if (ItemsInLIst(pkList)<1)
		print "ERROR -- FitPeakOnGraph(), cannot find any peak waves"
		return 1
	endif

	Variable fitWidth=1
	String wListAll=TraceNameList(gName,";",1)
	String yname, pname, wList=""
	Variable Nall=ItemsInList(wListAll), N,i
	for (i=0,N=0;i<Nall;i+=1)
		yname = StringFromList(i,wListAll)
		Wave ww = $yname
		if (StringMatch(yname,"fit_*") || strsearch(note(ww),yname,0)==0)	// reject fit results
			continue
		endif
		wList += yname+";"
	endfor

	Wave yPeak=$"", yData=$"", xw=$""
	if (ItemsInLIst(pkList)==1)
		pname = StringFromList(0,pkList)
		Wave yPeak = $pname
	endif

	if (ItemsInList(wList)<1)				// no data
		return 1
	elseif (ItemsInList(wList)==1 && WaveExists(yPeak))	// nothing to ask
		yname = StringFromList(0,wList)
		Wave yData = $yname
	else
		Prompt fitWidth, "Fit width of reolution peak",popup,"Fixed Resolution Width;Resolution Width can Change"
		Prompt pname,"Peak Resolution Wave", popup, pkList
		Prompt yname,"Data Wave", popup, wList
		fitWidth = fitWidth ? 2 : 1
		DoPrompt "Data", yname,pname,fitWidth
		if (V_flag)
			return 1
		endif
		Wave yData = $yname
		Wave yPeak = $pname
		fitWidth = fitWidth==2
	endif
	if (!WaveExists(yPeak))
		print "ERROR -- FitPeakOnGraph(), cannot find the peak wave"
		return 1
	elseif (!WaveExists(yData))
		print "ERROR -- FitPeakOnGraph(), cannot find the data wave"
		return 1
	endif
	Wave xData = XWaveRefFromTrace(gName,yname)
	Wave xPeak = XWaveRefFromTrace(gName,pname)
	Variable/C pz = GetCursorRangeFromGraph(gName,0)

	if (numtype(pz)==0)
		FitPeakWave(yData,xData,yPeak,xPeak, fitW=fitWidth, printIt=1,pLo=real(pz),pHi=imag(pz))
	else
		FitPeakWave(yData,xData,yPeak,xPeak, fitW=fitWidth, printIt=1)
	endif

	return 0
End



// This routine is used when fitting programactically
Function/WAVE FitPeakWave(yw,xw,yp,xp,[fitW,pLo,pHi,printIt,peakStruct])
	Wave yw, xw				// y and x data to be fit, (xw is optional, otherwise y-data assumed scaled)
	Wave yp, xp				// y and x values of the reference peak wave (xp is optional, but xp MUST be monotonic)
	Variable fitW			// flag, True=fit width, False=width is fixed
	Variable pLo,pHi		// optional point range to fit
	Variable printIt
	STRUCT PeakShapeStructure &peakStruct
	fitW = ParamIsDefault(fitW) || numtype(fitW) ? 0 : fitW
	pLo = ParamIsDefault(pLo) || numtype(pLo) ? -Inf : pLo
	pHi = ParamIsDefault(pHi) || numtype(pHi) ? Inf : pHi
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : !(!printIt)

	if (!WaveExists(yw))
		print "ERROR -- Input data wave does not exist"
		return $""
	elseif (!WaveExists(yp))
		print "ERROR -- Input peak wave does not exist"
		return $""
	endif

	STRUCT PeakWaveFitStruct fs
	Make/N=(fitW ? 4 : 3)/D/FREE cw
	Wave sw = CalcSimplePeakParameters(yw,xw,1)	// {bkg, amplitude, x0, FWHM, net, COM}
	Wave sp = CalcSimplePeakParameters(yp,xp,1)	// {bkg, amplitude, x0, FWHM, net, COM}

	cw[0] = sw[0]							// background
	cw[0] = abs(cw[0]) > 0 ? cw[0] : sw[1]/100
	cw[1] = sw[1] / sp[1]				// starting amplitude
	cw[2] = sw[2]							// starting xOffset
	cw[2] = abs(cw[2]) > 0 ? cw[2] : sw[3]/10
	if (abs(cw[1])==0)
		print "ERROR -- Input data has no peak, it is flat"
		return $""
	endif

	Wave fs.yPeak = yp
	Wave fs.xPeak = xp
	fs.X0 = NumberByKey("X0",note(yp),"=")
	fs.X0 = numtype(fs.X0) ? sp[2] : fs.X0
	if (fitW)
		fs.FWHM0 = NumberByKey("FWHM",note(yp),"=")
		fs.FWHM0 = numtype(fs.FWHM0) ? sp[3] : fs.FWHM0
		if (numtype(fs.FWHM0))
			print "ERROR -- Unable to find FWHM of the peak wave"
			return $""
		endif
		cw[3] = numtype(sw[3]) ? fs.FWHM0 : sw[3]
	else
		fs.FWHM0 = NaN
	endif
	Duplicate/FREE cw, cwStart

	if (WaveExists(xp))
		WaveStats/M=1/Q xp
		fs.xLo = V_min
		fs.xHi = V_max
		fs.iLo = V_minRowLoc
		fs.iHi = V_maxRowLoc
	else
		fs.iLo = 0
		fs.iHi = numpnts(yp)-1
		if (DimDelta(yp,0)<0)
			Variable swap = fs.iLo
			fs.iLo = fs.xHi
			fs.iHi = swap
		endif
		fs.xLo = pnt2x(yp,fs.iLo)
		fs.xHi = pnt2x(yp,fs.iHi)
	endif

	String S_Info=""
	Variable V_FitError=0, V_FitQuitReason=0
	if (ItemsInLIst(FindGraphsWithWave(yw)))
		FuncFit fwhm#fitUsingPeakWaveFunc, cw, yw[pLo,pHi]/X=xw/NWOK/D /STRC=fs
	else
		FuncFit fwhm#fitUsingPeakWaveFunc, cw, yw[pLo,pHi]/X=xw/NWOK /STRC=fs
	endif
	Wave W_sigma

	String wnote="waveClass:fitCoefsPeakWave;"
	wnote = ReplaceNumberByKey("V_chisq",wnote,V_chisq,"=")
	wnote = ReplaceNumberByKey("V_FitError",wnote,V_FitError,"=")
	wnote = ReplaceNumberByKey("V_FitQuitReason",wnote,V_FitQuitReason,"=")
	if (WaveExists(W_sigma))
		wnote = ReplaceStringByKey("W_sigma",wnote,vec2str(W_sigma,sep=","),"=")
	endif
	wnote = MergeKeywordLists(wnote,S_Info,0,"=",";")
	Note/K cw, wnote

	STRUCT PeakShapeStructure pkLocal
	initPeakShapeStructure(pkLocal)
	pkLocal.N = numpnts(cw)
	if (V_FitError==0)
		Duplicate/O cw, W_coef
		Variable i
		for (i=0; i<pkLocal.N; i+=1)
			pkLocal.coef[i] = cw[i]
			pkLocal.sigma[i] = W_sigma[i]
		endfor
		pkLocal.x0 = cw[2]
		pkLocal.dx0 = W_sigma[2]
		if (fitW)
			pkLocal.FWHM = cw[3]
			pkLocal.dFWHM = W_sigma[3]
		else
			pkLocal.FWHM = sp[3]
		endif
		pkLocal.amp = cw[1]
		pkLocal.damp = W_sigma[1]
		pkLocal.bkg = cw[0]
		pkLocal.dbkg = W_sigma[0]
		pkLocal.type = "PeakWave"
		pkLocal.COM = computeCOM(yw,xw)
		pkLocal.net = sp[4]*cw[1]
		pkLocal.dnet = sp[4]*W_sigma[1]
		pkLocal.yunits = WaveUnits(yw,-1)
		if (WaveExists(xw))
			pkLocal.xunits = WaveUnits(xw,-1)
		else
			pkLocal.xunits = WaveUnits(yw,0)
		endif
	endif
	if (!ParamIsDefault(peakStruct))			// copy to peakStruct if one was passed
		copyPeakShapeStructure(peakStruct,pkLocal)
	endif

	if (printIt)
		String inData, peakData
		if (WaveExists(xw))
			sprintf inData, "%s vs %s",prettyWaveName(yw),prettyWaveName(xw)
		else
			inData = prettyWaveName(yw)
		endif
		if (WaveExists(xp))
			sprintf peakData, "%s vs %s",prettyWaveName(yp),prettyWaveName(xp)
		else
			peakData = prettyWaveName(yp)
		endif
		printf "Fit of data={%s}  to  peak={%s}\r",inData,peakData

		printf "  Started with coefs = %s\r",vec2str(cwStart)
		if (V_FitError)
			printf "***** Fit Failed with V_FitError = %g,   and  V_FitQuitReason = %g\r",V_FitError,V_FitQuitReason
		else
			printf "  After fit (chiSq=%g),  W_coef = %s,  W_sigma = %s\r",V_chisq,vec2str(W_coef),vec2str(W_sigma)
			print PeakShapeStruct2str(pkLocal)
		endif
	endif
	return cw
End
//
Static Structure PeakWaveFitStruct
	Wave cw							// {bkg, amp, xOffset}  or  {bkg, amp, xOffset, fwhm}
	Variable x
	STRUCT WMFitInfoStruct fi	// optional WMFitInfoStruct
	Wave yPeak
	Wave xPeak						// optional, if not present, assume yPeak is scaled
	Variable X0						// center of {yPeak,xPeak}, usually X0=0
	Variable FWHM0					// FWHM of yPeak, as given (only used when fitting fwhm)
	Variable xLo,iLo				// lowest x value and index into yp where it occurs
	Variable xHi,iHi				// highest x value and index into yp where it occurs
EndStructure
//
Static Function fitUsingPeakWaveFunc(s) : FitFunc		// the fitting function
	Struct PeakWaveFitStruct &s
	// cw = {bkg, amp, xOffset}  or  cw = {bkg, amp, xOffset, fwhm}
	Variable yval, xval=s.x - s.cw[2] + s.X0
	xval *= numpnts(s.cw)>3 ? s.FWHM0/s.cw[3] : 1
	if (xval<=s.xLo)
		yval = s.yPeak[s.iLo]
	elseif (xval>=s.xHi)
		yval = s.yPeak[s.iHi]
	else
		yval = WaveExists(s.xPeak) ? s.yPeak[BinarySearchInterp(s.xPeak,xval)] : s.yPeak(xval)
	endif
	return s.cw[0] + s.cw[1] * yval
end
//
Static Function/S prettyWaveName(ww)
	Wave ww

	if (!WaveExists(ww))
		return ""
	elseif (WaveType(ww,2)==2)
		return "FREE"
	endif

	String strFull= GetWavesDataFolder(ww,2)
	String strRoot=ReplaceString("root:",strFull,"",0,1)
	String strLocal=ReplaceString(GetDataFolder(1),strFull,"",0,1)

	if (strsearch(strRoot,":",0)<0)
		return strRoot
	elseif (strsearch(strLocal,":",0)<0)
		return strLocal
	endif
	return strFull
End
//Static Function/S prettyWaveName(ww)
//	Wave ww
//
//	if (!WaveExists(ww))
//		return ""
//	elseif (WaveType(ww,2)==2)
//		return "FREE"
//	endif
//
//	String str0= GetWavesDataFolder(ww,2)
//	String str=ReplaceString("root:",str0,"",0,1)
//	return SelectString(strsearch(str,":",0)<0, str0, str)
//End

//Function Make_yPeakTestData()
//	Make/N=51/D/O yPeak
//	SetScale/I x -2,2,"", yPeak
//	yPeak = exp(-(x/0.5)^2)
//	Variable ar = area(ypeak)
//	ypeak /= ar
//	String wnote="waveClass=peakWave;"
//	wnote = ReplaceNumberByKey("FWHM",wnote,0.832554611157698,"=")
//	wnote = ReplaceNumberByKey("X0",wnote,0,"=")
//	Note/K yPeak, wnote
//
//	Make/N=500/O yData, xData
//	SetScale/I x -20,20,"", xData
//	xData = x
//	SetScale/P x 0,1,"", yData, xData
//	yData = 20*exp(-((xData[p]-4)/2)^2) + gnoise(2)
//End

//  ============================== End of Fit a Peak Wave  ===============================  //
//  ======================================================================================  //




//  ======================================================================================  //
//  ================================ Start of basic FWHM  ================================  //

Function FWHM_Coefs([type,printIt])
	// uses W_coef to return FWHM, and probably print other statistics
	String type						// user can pass type, or it is read from wave note
	Variable printIt
	type = SelectString(ParamIsDefault(type), type, "")
	printIt = ParamIsDefault(printIt) ? NaN : printIt
	printIt = numtype(printIt) ? 1 : !(!printIt)

	Wave W_coef=W_coef, W_sigma=W_sigma
	if (!WaveExists(W_coef))
		return NaN
	endif
	Variable Nc=DimSize(W_coef,0)
	if (DimSize(W_sigma,0)!=Nc)
		Wave W_sigma=$""
	endif

	String FitList=""
	if (Nc==4)
		FitList="Lorentzian;Gaussian;"
	elseif (Nc==6)
		FitList="Simple;"
	elseif (Nc==5)
		FitList="Voigt;PearsonVII;"
	else
		return NaN
	endif

	type = SelectString(strlen(type), StringByKey("type",note(W_coef),"="), type)
	if (strlen(type)<1)
		if (ItemsInList(FitList)<1)
			return NaN
		elseif (ItemsInList(FitList)==1)
			type = StringFromList(0,FitList)
		elseif (ItemsInList(FitList)>1)
			Prompt type,"type of W_coef", popup, FitList
			DoPrompt "W_coef type",type
			if (V_flag)
				return NaN
			endif
		endif
		printIt = 1
	endif

	STRUCT PeakShapeStructure pk
	if (GenericPeakShapeStruct(pk,W_coef,W_sigma,type=type))
		if (printIt)
			printf "ERROR -- Unable to use W_coef of type \"%s\"\r",type
		endif
		return NaN
	endif
	if (printIt)
		print PeakShapeStruct2str(pk)
	endif
	return pk.FWHM
End


Function FindSimplePeakParameters(yw,xw,[W_coef,useBkg,printIt])
	// a user command to get Simple peak stats from yw (or yw vs xw), and probably print them
	// it also creates the W_coef & W_sigma waves containing {bkg, amplitude, x0, FWHM, net, COM}
	// returns 0 if all went OK, 1 on an error
	Wave yw,xw			// yw vs xw, or just yw using scaling
	Wave/D W_coef		// wave to hold the results, one will be created if not present
	Variable useBkg											// if true, find a background, if false NO Background
	Variable printIt
	useBkg = ParamIsDefault(useBkg) ? NaN : useBkg
	useBkg = numtype(useBkg) ? 1 : !(!useBkg)		// generic default is to use the background
	printIt = ParamIsDefault(printIt) ? NaN : printIt
	printIt = numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : !(!printIt)

	if (!WaveExists(yw))
		String yName,xName
		Prompt yName, "input wave", popup, WaveList("*",";","DIMS:1")
		Prompt xName, "input Xwave", popup, "-calculated-;"+WaveList("*",";","DIMS:1")
		Prompt useBkg, "background", popup, "find the background;no background"
		DoPrompt "wave(s) defining peak", yName,xName,useBkg
		if (V_flag)
			return 1
		endif
		useBkg = useBkg==1
		Wave yw=$yName, xw=$xName
		printIt = 1
	endif
	if (printIt)
		yName = SelectString(WaveExists(yw),"$\"\"",NameOfWave(yw))
		xName = SelectString(WaveExists(xw),"$\"\"",NameOfWave(xw))
		printf "FindSimplePeakParameters(%s, %s",yName,xName
		if (!ParamIsDefault(printIt))
			printf ", printIt=%g",printIt
		endif
		printf ")\r"
	endif
	if (!WaveExists(yw))
		return 1
	endif

	Wave Wc = CalcSimplePeakParameters(yw,xw,useBkg)
	Variable Np=DimSize(Wc,0)
	if (WaveExists(W_coef))
		Redimension/N=(Np)/D W_coef
	else
		Make/N=(Np)/O/D W_coef=NaN
	endif
	W_coef = Wc
	CopyScales Wc, W_coef
	Note/K W_coef, note(Wc)
	Duplicate/O W_coef, W_sigma
	W_sigma = NaN					// no known errors yet

	if (printIt)
		STRUCT PeakShapeStructure pk
		Simple2PeakShapeStruct(pk,W_coef,$"")
		print PeakShapeStruct2str(pk)
	endif
	return 0
End
//
ThreadSafe Function/WAVE CalcSimplePeakParameters(yw,xw,useBkg)	// look at peak in yw vs xw, put results into W_coef
	// calculate the peak stats for peak in yw, or yw vs xw, and returns result as FREEE wave
	// Wc = {bkg, amplitude, x0, FWHM, net, COM}   Simple peak stats
	Wave yw,xw
	Variable useBkg											// if true, find a background, if false NO Background
	useBkg = numtype(useBkg) ? 1 : !(!useBkg)		// generic default is to use the background
	if (!WaveExists(yw))
		return $""
	endif

	String xunits, yunits=WaveUnits(yw,-1)
	if (WaveExists(xw))
		xunits = WaveUnits(xw,-1)
	else
		xunits = WaveUnits(yw,0)
	endif

	Variable level, pStart, p1,p2						// points values of the half widths
	Variable dx=DimDelta(yw,0), xoff=DimOffset(yw,0), N=DimSize(yw,0)
	Variable bkg, amp, center, FWHM, net, COM		// peak parameters to find
	COM = computeCOM(yw,xw)								// computes center of mass
	bkg = useBkg ? ( yw[0] + yw[N-1] ) / 2 : 0
	WaveStats/Q yw
	level = ((V_max+V_min)/2)								// level of the half-width
	Variable dip=abs(V_max-bkg) < abs(V_min-bkg)	//a negative peak, a dip
	if (dip)														//a negative peak
		amp = V_min-bkg
		pStart = round((V_minloc-xoff)/dx)				// pStart in points (lowest point)
	else															// a positive peak
		amp = V_max-bkg
		pStart = round((V_maxloc-xoff)/dx)				// pStart in points (highest point)
	endif
	p1 = FindLevelFromPoint(yw,level,pStart,-1,dip)	// find left side
	p2 = FindLevelFromPoint(yw,level,pStart,1,dip)	// find right side

	if (numtype(p1+p2))
		return $""
	elseif (WaveExists(xw))
		FWHM = abs(xw[p2]-xw[p1])
		center=(xw[p2]+xw[p1])/2
		net = areaXY(xw,yw) - (xw[N-1]-xw[0])*bkg
	else
		FWHM = abs((p2-p1)*dx)
		center = dx*(p2+p1)/2 + xoff
		net = area(yw) - bkg*(N-1)*dx
	endif

	Make/N=6/D/FREE Wc={bkg, amp, center, FWHM, net, COM}
	SetScale/P x 0,1, xunits, Wc
	SetScale d 0,0, yunits, Wc

	String wnote="waveClass=peakFitCoefs;"
	wnote = ReplaceStringByKey("type",wnote,"Simple","=")
	wnote = ReplaceStringByKey("yWave",wnote,GetWavesDataFolder(yw,2),"=")
	if (WaveExists(xw))
		wnote = ReplaceStringByKey("xWave",wnote,GetWavesDataFolder(xw,2),"=")
	endif
	if (strlen(xunits))
		wnote = ReplaceStringByKey("xunits",wnote,xunits,"=")
	endif
	if (strlen(yunits))
		wnote = ReplaceStringByKey("yunits",wnote,yunits,"=")
	endif
	Note/K Wc, wnote
	return Wc
End
//
ThreadSafe Static Function FindLevelFromPoint(wy,level,p0,direction,dip)
	// returns point where wy[point] crosses level
	// it searches for the crossing from p0 going in the indicated direction from p0
	// if wy[p0] > level,  then it assumes a positive peak
	// if wy[p0] < level,  then it assumes a negative peak (a dip)
	// if wy[p0] == level, then it returns p0
	Wave wy					// wave with values, assumed to be 1D
	Variable level			// value in wy[] to find
	Variable p0				// point at start of search
	Variable direction	// (+ is >0), (- is <=0) direction
	Variable dip			// true for a dip, a positive peak is false
								// if numtype(dip) is true, then auto-detect dips

	if (!WaveExists(wy) || numtype(level+p0+direction) || p0<0)
		return NaN
	endif
	Variable i,N=DimSize(wy,0), point=NaN
	p0 = round(p0)
	if (N<p0)								// starting point out of range
		return NaN
	elseif (wy[p0]==level)				// we are there already, no need to hunt
		return p0
	endif
	dip = numtype(dip) ? wy[p0]<level : dip	// set dip

	WaveStats/Q/M=1 wy
	if (!(level==limit(level,V_min,V_max)))
		return NaN							// level is not in range of wy
	endif

	if (dip)									// negative peak
		if (direction>0)					// + direction
			for (i=p0; i<N && wy[i]<level; i+=1)
			endfor
		else									// - direction
			for (i=p0; i>=0 && wy[i]<level; i-=1)
			endfor
		endif
	else										// positive peak
		if (direction>0)					// + direction
			for (i=p0; i<N && wy[i]>level; i+=1)
			endfor
		else									// - direction
			for (i=p0; i>=0 && wy[i]>level; i-=1)
			endfor
		endif
	endif

	if (!(i==limit(i,0,N-1)))			// check if i is a valid point in wave
		point = NaN
	elseif (i==p0 || wy[i]==level)	// no need to interpolate, this is the answer
		point = i
	else										// wy[i]<level, level<wy[i+1]
		i -= direction>0 ? 1 : 0
		point = (level-wy[i])/(wy[i+1]-wy[i]) + i
	endif
	return point
End

//  ================================= End of basic FWHM  =================================  //
//  ======================================================================================  //




//  ============================ Start of PeakShapeStructure  ============================  //
//  ======================================================================================  //

// Generic description of a peak
Structure PeakShapeStructure
	int32 N					// number of values in coef[] & sigma[]
	double coef[MAX_PEAK_COEFS]	// original coefs from W_coef, only up to MAX_PEAK_COEFS values
	double sigma[MAX_PEAK_COEFS]	// original sigmas from W_sigma (optional)
	double x0, dx0			// center of peak & the error
	double FWHM	, dFWHM	// Full Width Half Max & the error
	double amp, damp		// amplitude & the error
	double bkg, dbkg		// background value & the error
	double net, dnet		// net area of peak (background subtracted) & the error
	double shape, dshape// optional shape parameter (usually only for Voigt or PearsonVII) & the error
	double shape1, dshape1// optional extra shape parameter, porbably not used & the error
	double COM, dCOM		// center of mass & the error
	char xunits[40]		// optional x-units
	char yunits[40]		// optional y-units
	char type[40]			// type of shape, e.g. "Simple;Lorentzian;Gaussian;Voigt;PearsonVII"
EndStructure
//
ThreadSafe Function copyPeakShapeStructure(f,i)
	// copy PeakShapeStructure i into f
	STRUCT PeakShapeStructure &f, &i
	f.N = i.N
	f.x0 = i.x0 ;			f.dx0 = i.dx0
	f.FWHM = i.FWHM ;		f.dFWHM = i.dFWHM
	f.amp = i.amp ;		f.damp = i.damp
	f.bkg = i.bkg ;		f.dbkg = i.dbkg
	f.net = i.net ;		f.dnet = i.dnet
	f.shape = i.shape ;	f.dshape = i.dshape
	f.shape1 = i.shape1;f.dshape1 = i.dshape1
	f.COM = i.COM ;		f.dCOM = i.dCOM
	f.xunits = i.xunits;f.yunits = i.yunits
	f.type = i.type
	Variable m
	for (m=0;m<MAX_PEAK_COEFS;m+=1)
		f.coef[m] = i.coef[m]
		f.sigma[m] = i.sigma[m]
	endfor
End
//
ThreadSafe Function initPeakShapeStructure(pk)
	// initialze PeakShapeStructure to default/empty values
	STRUCT PeakShapeStructure &pk
	pk.N = 0
	pk.x0 = NaN ;		pk.dx0 = NaN
	pk.FWHM = NaN ;	pk.dFWHM = NaN
	pk.amp = NaN ;		pk.damp = NaN
	pk.bkg = NaN ;		pk.dbkg = NaN
	pk.net = NaN ;		pk.dnet = NaN
	pk.shape = NaN ;	pk.dshape = NaN
	pk.shape1 = NaN;	pk.dshape1 = NaN
	pk.COM = NaN ;		pk.dCOM = NaN
	pk.xunits = "" ;	pk.yunits = ""
	pk.type = ""
	Variable m
	for (m=0;m<MAX_PEAK_COEFS;m+=1)
		pk.coef[m] = NaN
		pk.sigma[m] = NaN
	endfor
End
//
ThreadSafe Function/T PeakShapeStruct2str(pk,[short])
	// return a string suitable for printing or display of a PeakShapeStructure
	STRUCT PeakShapeStructure &pk
	Variable short
	short = ParamIsDefault(short) || numtype(short) ? 1 : !(!short)

	String out,str, type=pk.type
	type = SelectString(strlen(type),"Unknown",type)
	if (WhichListItem(type,"Simple;Lorentzian;Gaussian;")>0)	// no shape parameters for this peak
		pk.shape = NaN
		pk.shape1 = NaN
	endif
	if (WhichListItem(type,"Voigt;PearsonVII;")>0)		// no shape1 parameter for this peak
		pk.shape = NaN
	endif
	String xunits = SelectString(strlen(pk.xunits), "", "("+pk.xunits+")")
	String yunits = SelectString(strlen(pk.yunits), "", "("+pk.yunits+")")
	String xyunits = ""
	if (strlen(xunits+yunits))
		xyunits = pk.yunits
		xyunits += SelectString(strlen(xyunits),""," x ")
		xyunits += pk.xunits
		xyunits = "("+xyunits+")"
	endif
	String sFWHM=ValErrStr(pk.FWHM,pk.dFWHM)+xunits, sCen=ValErrStr(pk.x0,pk.dx0)+xunits
	String sAmp=ValErrStr(pk.amp,pk.damp)+yunits, sBkg=ValErrStr(pk.bkg,pk.dbkg)+yunits
	String sNet=ValErrStr(pk.net,pk.dnet)+xyunits

	if (short)
		sprintf out, "%s:  FWHM=%s, Center=%s, Amp=%s, Bkg=%s Integral=%s",type,sFWHM,sCen,sAmp,sBkg,sNet
	else
		sprintf out, "%s:  FWHM=%s, Center=%s\r    Amp=%s, Bkg=%s Integral=%s",type,sFWHM,sCen,sAmp,sBkg,sNet
	endif
	if (numtype(pk.COM)==0)
		sprintf str,",  COM = %s", ValErrStr(pk.COM,pk.dCOM)
		out += str
	endif
	if (numtype(pk.shape)==0)
		sprintf str,"%sshape = %s", SelectString(short, "\r",",  "),ValErrStr(pk.shape,pk.dshape)
		out += str
	endif
	if (numtype(pk.shape1)==0)
		sprintf str,",  shape1 = %s", ValErrStr(pk.shape1,pk.dshape1)
		out += str
	endif
	return out
End
//
ThreadSafe Function GenericPeakShapeStruct(pk,wcoef,wsigma,[type])
	// given a W_coef (and maybe W_sigma or type), fill the PeakShapeStructure
	STRUCT PeakShapeStructure &pk
	Wave wcoef, wsigma
	String type											// user can pass type, or it is read from wave note
	type = SelectString(ParamIsDefault(type), type, "")
	if (WaveType(wcoef,1)!=1)						// wcoef exists, and it is numeric
		return 1
	endif

	if (strlen(type)<1)
		type=StringByKey("type",note(wcoef),"=")
	endif
	if (stringmatch(type"Simple*"))
		return Simple2PeakShapeStruct(pk,wcoef,wsigma)
	elseif (stringmatch(type"Gauss*"))
		return Gaussian2PeakShapeStruct(pk,wcoef,wsigma)
	elseif (stringmatch(type"Lorentz*"))
		return Lorentzian2PeakShapeStruct(pk,wcoef,wsigma)
	endif
//	printf "ERROR -- in GenericPeakShapeStruct(...), Unable to identify type of Peak Fit from \"%s\"\r",NameOfWave(wcoef)
	return 1
End
//
ThreadSafe Static Function Simple2PeakShapeStruct(pk,wcoef,wsigma)
	// for a Simple peak result in wcoef (and optionally wsigma), fill in PeakShapeStructure
	// y = Simple peak,  K0=bkg, K1=amplitude, K2=x0, K3=FWHM
	STRUCT PeakShapeStructure &pk
	Wave wcoef, wsigma
	if (WaveType(wcoef,1)!=1)						// wcoef exists, and it is numeric
		return 1
	elseif (DimSize(wcoef,0)<6)
		return 1
	endif
	initPeakShapeStructure(pk)
	Variable i,N=min(DimSize(wcoef,0),6)
	pk.N = N
	for (i=0;i<N;i+=1)
		pk.coef[i] = wcoef[i]
		if (WaveExists(wsigma) && DimSize(wsigma,0)>i)
			pk.sigma[i] = wsigma[i]
		endif
	endfor
	String wnote=note(wcoef)
	String type=StringByKey("type",wnote,"=")
	if (strlen(type)>1 && !StringMatch(type, "Simple*"))
		printf "ERROR -- Simple2PeakShapeStruct(...), class should be \"Simple\", but it is \"%s\"\r",type
	endif

	pk.type = "Simple"
	pk.amp  = wcoef[1]
	pk.x0   = wcoef[2]
	pk.bkg  = wcoef[0]
	pk.FWHM = wcoef[3]
	pk.net  = wcoef[4]
	pk.COM  = wcoef[5]
	if (DimSize(wsigma,0)>=5)
		pk.damp  = wsigma[1]
		pk.dx0   = wsigma[2]
		pk.dbkg  = wsigma[0]
		pk.dFWHM = wsigma[3]
		pk.dnet  = wsigma[4]
		pk.dCOM  = wsigma[5]
	endif

	String xunits=WaveUnits(wcoef,0), yunits=WaveUnits(wcoef,-1)
	if (strlen(xunits)<1)
		xunits = StringByKey("xunits",wnote,"=")
	endif
	if (strlen(yunits)<1)
		yunits = StringByKey("yunits",wnote,"=")
	endif
	pk.xunits = xunits
	pk.yunits = yunits
	return 0
End
//
ThreadSafe Static Function Gaussian2PeakShapeStruct(pk,wcoef,wsigma)
	// for a Gaussian peak fit (given by wcoef & wsigma), fill in PeakShapeStructure
	// y = K0+K1*exp(-((x-K2)/K3)^2)
	STRUCT PeakShapeStructure &pk
	Wave wcoef, wsigma
	if (WaveType(wcoef,1)!=1)						// wcoef exists, and it is numeric
		return 1
	elseif (DimSize(wcoef,0)<4)
		return 1
	endif
	initPeakShapeStructure(pk)
	Variable i,N=min(DimSize(wcoef,0),4)
	pk.N = N
	for (i=0;i<N;i+=1)
		pk.coef[i] = wcoef[i]
		if (WaveExists(wsigma))
			pk.sigma[i] = wsigma[i]
		endif
	endfor
	pk.type = "Gaussian"
	pk.amp = wcoef[1]
	pk.x0 = wcoef[2]
	pk.bkg = wcoef[0]
	pk.FWHM = 2*wcoef[3]*sqrt(ln(2))
	pk.net = wcoef[1] * abs(wcoef[3]) * sqrt(PI)
	pk.xunits = WaveUnits(wcoef,0)
	pk.yunits = WaveUnits(wcoef,-1)
	if (DimSize(wsigma,0)>=4)
		pk.damp = wsigma[1]
		pk.dx0 = wsigma[2]
		pk.dbkg = wsigma[0]
		pk.dFWHM = 2*wsigma[3]*sqrt(ln(2)) 
		pk.dnet = pk.net * sqrt((wsigma[1]/wcoef[1])^2 + (wsigma[3]/wcoef[3])^2)
	endif
	return 0
End
//
ThreadSafe Static Function Lorentzian2PeakShapeStruct(pk,wcoef,wsigma)
	// for a Lorentzian peak fit (given by wcoef & wsigma), fill in PeakShapeStructure
	// y = K0+K1/((x-K2)^2+K3)
	STRUCT PeakShapeStructure &pk
	Wave wcoef, wsigma
	if (WaveType(wcoef,1)!=1)						// wcoef exists, and it is numeric
		return 1
	elseif (DimSize(wcoef,0)<4)
		return 1
	endif
	initPeakShapeStructure(pk)
	Variable i,N=min(DimSize(wcoef,0),4)
	pk.N = N
	for (i=0;i<N;i+=1)
		pk.coef[i] = wcoef[i]
		if (WaveExists(wsigma))
			pk.sigma[i] = wsigma[i]
		endif
	endfor
	pk.type = "Lorentzian"
	pk.amp = wcoef[1]/wcoef[3]
	pk.x0 = wcoef[2]
	pk.bkg = wcoef[0]
	pk.FWHM = 2*sqrt(wcoef[3])
	pk.net = PI*wcoef[1]/sqrt(wcoef[3])
	pk.xunits = WaveUnits(wcoef,0)
	pk.yunits = WaveUnits(wcoef,-1)
	if (DimSize(wsigma,0)>=4)
		pk.damp = pk.amp * sqrt((wsigma[1]/wcoef[1])^2 + (wsigma[3]/wcoef[3])^2)
		pk.dx0 = wsigma[2]
		pk.dbkg = wsigma[0]
		pk.dFWHM = wsigma[3]
		pk.dnet = pk.net * sqrt((wsigma[1]/wcoef[1])^2 + (wsigma[3]/wcoef[3]/2)^2)
	endif
	return 0
End

//  ======================================================================================  //
//  ============================= End of PeakShapeStructure  =============================  //




// ==============================================================================================
// ================ DEPRECATED  DEPRECATED  DEPRECATED  DEPRECATED  DEPRECATED  =================
// ================== The Folowing Functions are DERECATED, use the ones above ==================
// ==============================================================================================


Function FWHM_Gaussian()
	Wave W_coef=W_coef
	STRUCT PeakShapeStructure pk
	Variable/G V_FWHM=NaN
	if (Gaussian2PeakShapeStruct(pk,W_coef,$"")==0)
		V_FWHM = pk.FWHM
	endif
	return V_FWHM
End
//Function FWHM_Gaussian()
//	Variable/G V_FWHM
//	Wave W_coef=W_coef
//	if (!WaveExists(W_coef))
//		V_FWHM = NaN
//		return NaN
//	endif
//	V_FWHM = 2*W_coef[3]*sqrt(ln(2))
//	return V_FWHM
//End

Function FWHM_Lorentz()
	Wave W_coef=W_coef
	STRUCT PeakShapeStructure pk
	Variable/G V_FWHM=NaN
	if (Lorentzian2PeakShapeStruct(pk,W_coef,$"")==0)
		V_FWHM = pk.FWHM
	endif
	return V_FWHM
End
//Function FWHM_Lorentz()
//	Variable/G V_FWHM
//	Wave W_coef=W_coef
//	if (!WaveExists(W_coef))
//		V_FWHM = NaN
//		return NaN
//	endif
//	V_FWHM = 2*sqrt(W_coef[3])
//	return V_FWHM
//End


Function FWHM_peak(ww)
	Wave ww
	Wave Wc = CalcSimplePeakParameters(ww,$"",0)	// look at pea in yw vs xw, put results into W_coef
	if (!WaveExists(Wc))
		return NaN
	endif
	Variable/G V_center=Wc[2]
	return Wc[3]
End
//Function FWHM_peak(ww)
//	Wave ww
//	WaveStats/Q ww
//	Make/N=2/FREE/D ww_temp__
//	FindLevels /D=ww_temp__ /N=6 /Q ww, ((V_max+V_min)/2)
//	if (V_LevelsFound<2)
//		return NaN					// no peak
//	endif
//
//	Variable p1,p2				// points values of the half widths
//	Variable x1,x2				// x values of halfwidths
//	p1 = BinarySearch(ww_temp__, V_maxloc)
//	x2 = ww_temp__[p1+1]
//	x1 = ww_temp__[p1]
//	Variable/G V_center=(x2+x1)/2
//	return abs(x2-x1)
//End


Function FWHM_peakXY(wy,wx)
	Wave wy,wx
	if (!WaveExists(wx))
		return NaN
	endif
	Wave Wc = CalcSimplePeakParameters(wy,wx,0)	// look at pea in yw vs xw, put results into W_coef
	if (!WaveExists(Wc))
		return NaN
	endif
	Variable/G V_center=Wc[2]
	return Wc[3]
End
//Function FWHM_peakXY(wy,wx)
//	Wave wy,wx
//
//	WaveStats/Q wy
//	Make/N=2/FREE/D ww_temp__
//	FindLevels /D=ww_temp__ /N=6 /P/Q wy, ((V_max+V_min)/2)
//	if (V_LevelsFound<2)
//		return NaN					// no peak
//	endif
//
//	Variable pm,p1,p2, dx			// points of mid, HW point, and width
//	pm = x2pnt(wy,V_maxloc)		// point of the max
//	p1 = BinarySearch(ww_temp__, pm)
//	p2 = ww_temp__[p1+1]
//	p1 = ww_temp__[p1]
//	dx = abs(wx[p2]-wx[p1])
//	Variable/G V_center=(wx[p2]+wx[p1])/2
//	return dx
//End


Function PlainFWHMofWave(yName,xName)
	String yName, xName

	Make/N=5/D/FREE wC
	if (FindSimplePeakParameters($yName,$xName,W_coef=wC))
		return NaN
	endif
	Variable/G V_FWHM=wC[3], V_center=wC[2]
	return wC[3]
End
//Function PlainFWHMofWave(yw,xw)
//	String yw, xw
//	if (stringmatch(xw,"none"))
//		xw = "_calculated_"
//	endif
//	if (exists(yw)!=1 || (exists(xw)!=1 && !stringmatch(xw,"_calculated_") ))
//		Prompt yw, "input wave", popup, WaveList("*",";","DIMS:1")
//		Prompt xw, "input Xwave", popup, "_calculated_;"+WaveList("*",";","DIMS:1")
//		DoPrompt "wave(s) defining peak", yw,xw
//	endif
//	if (V_flag)
//		Abort
//	endif
//
//	Variable/G V_FWHM
//	if (cmpstr("_calculated_",xw)==0)
//		V_FWHM=FWHM_peak($yw)
//	else
//		V_FWHM=FWHM_peakXY($yw,$xw)
//	endif
//	NVAR V_center=V_center
//	if (ItemsInList(GetRTStackInfo(0))<2)
//		printf "PlainFWHMofWave(\"%s\",\"%s\") = %g, centered at %g\r",yw,xw,V_FWHM,V_center
//	endif
//	return V_FWHM
//End


