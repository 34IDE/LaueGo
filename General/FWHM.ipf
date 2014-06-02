#pragma rtGlobals=1		// Use modern global access method.
#pragma IgorVersion = 6.11
#pragma version = 2.04
#pragma ModuleName=fwhm

// with v 2.0, major revision, started using structures
//		The old functions are still there, but with new innards

Static Constant MAX_PEAK_COEFS = 6


Menu "Analysis"
	"-"
	MenuItemIfWavesExists("Simple Peak Parameters of ...","*","DIMS:1,MINROWS:3"), FindSimplePeakParameters($"",$"")
	MenuItemIfWavesExists("Peak Parameters [W_coef]","W_coef","DIMS:1,MINROWS:4"), FWHM_Coefs()
End



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
	short = ParamIsDefault(short) ? NaN : short
	short = numtype(short) ? 0 : !(!short)

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
		sprintf out, "%s:  FWHM=%s, Center=%s\rAmp=%s, Bkg=%s Integral=%s",type,sFWHM,sCen,sAmp,sBkg,sNet
	else
		sprintf out, "%s:  FWHM=%s, Center=%s, Amp=%s, Bkg=%s Integral=%s",type,sFWHM,sCen,sAmp,sBkg,sNet
	endif
	if (numtype(pk.COM)==0)
		sprintf str,",  COM = %s", ValErrStr(pk.COM,pk.dCOM)
		out += str
	endif
	if (numtype(pk.shape)==0)
		sprintf str,"%sshape = %s", SelectString(short, ",  ","\r"),ValErrStr(pk.shape,pk.dshape)
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


