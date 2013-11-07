#pragma rtGlobals=1		// Use modern global access method.
#include "Utility_JZT"
#pragma ModuleName=JZTanalysis
#pragma version = 0.03


Menu "Analysis"
	"-"
	MenuItemIfWavesExists("Scale Factor Between Waves...","*","DIMS:1"), ScaleFactorBetweenWaves($"",$"",$"",$"",ask=1)
	MenuItemIfWavesExists("re-bin wave","*","DIMS:1"), rebin_wave("",NaN,NaN)
	MenuItemIfWavesExists("Combine Nearby Ponts","*","DIMS:1"), CombineNearbyPoints("","","",NaN,NaN)
End


Function rebin_wave(inwave,nchan,renormalize)
	String inwave					// name of wave to bin up
	Variable nchan					// number of points to bin together
	Variable renormalize		// set to 0 or 1
	if (!(nchan>0) || numtype(renormalize) || exists(inwave)!=1)
		nchan = (!(nchan>0)) ? 10 : nchan
		renormalize = numtype(renormalize) ? 0 : renormalize
		inwave = SelectString(exists(inwave)==1,"",inwave)
		Prompt inwave,"input wave for rebinning",popup, WaveList("*", ";", "")
		Prompt nchan,"number of points in input wave to sum together"
		Prompt renormalize,"renormalize result",popup,"true;false"
		renormalize += 1
		DoPrompt "re-bin wave",inwave,nchan,renormalize
		if (V_flag || nchan<=1)
			return 1
		endif
		renormalize -= 1
		printf "rebin_wave(\"%s\", %g, %g)\r",inwave,nchan,renormalize
	endif

	Wave in = $inwave
	String outwave = inwave+"_"+num2istr(nchan)
	Variable nin=numpnts(in)
	Variable nout = ceil(nin/nchan)			//number of points in outwave wave
	Variable lo_edge						//used for scaling the $outwave
	String W_xunits = StringByKey("XUNITS",WaveInfo($inwave, 0))

	Duplicate/O $inwave $outwave
	SetScale/P x 0,1,"", $outwave
	Wave out = $outwave
	out = sum(out,p*nchan,p*nchan+nchan-1)
	Redimension/N=(nout) $outwave
	lo_edge = (nchan-1)/2 * deltax(in) + leftx(in)	//find left edge for scaling
	SetScale/P x, (lo_edge), (deltax(in)*nchan), W_xunits, $outwave
	if (renormalize)
		out /= nchan						//renormalize the binned wave
	endif
End



Function ScaleFactorBetweenWaves(cc0,xx0,cc1,xx1,[x0,x1,p0,p1,ask,printIt])
	Wave cc0,xx0			// first set of (x,y) pairs
	Wave cc1,xx1			// second set of (x,y) pairs
	Variable x0, x1		// range in x values of (cc0,xx0) to use, cannot specify both x0 & p0
	Variable p0, p1		// range in points of (cc0,xx0) to use, defaults to whole range if none given
	Variable ask
	Variable printIt
	ask = ParamIsDefault(ask) ? 0 : ask
	ask = numtype(ask) ? 0 : ask
	Variable printItDefault = strlen(GetRTStackInfo(2))<1 || ask
	printIt = ParamIsDefault(printIt) ? printItDefault : printIt
	printIt = numtype(printIt) ? 0 : printIt
	if (!ParamIsDefault(p0) && !ParamIsDefault(x0))
		print "ERROR -- ScaleFactorBetweenWaves(), both p0 & x0 specified"
		return NaN
	elseif (!ParamIsDefault(p1) && !ParamIsDefault(x1))
		print "ERROR -- ScaleFactorBetweenWaves(), both p1 & x1 specified"
		return NaN
	endif
	x0 = ParamIsDefault(x0) ? -Inf : x0
	x0 = numtype(x0) ? -Inf : x0
	x1 = ParamIsDefault(x1) ? Inf : x1
	x1 = numtype(x1) ? Inf : x1
	p0 = ParamIsDefault(p0) ? NaN : p0
	p0 = numtype(p0) ? NaN : p0
	p1 = ParamIsDefault(p1) ? NaN : p1
	p1 = numtype(p1) ? NaN : p1
	String c0name=NameOfWave(cc0), c1name=NameOfWave(cc1), x0name=NameOfWave(xx0),x1name=NameOfWave(xx1)
	String wList = WaveList("*",";","DIMS:1")
	c0name = SelectString(strlen(c0name),StringFromList(0,wList),c0Name)
	c1name = SelectString(strlen(c1name),StringFromList(1,wList),c1name)
	Prompt c0name "Wave1 Y1 values", popup, wList
	Prompt c1name "Wave1 Y2 values", popup, wList
	Prompt x0name "Wave1 X1 values", popup, "_none_;"+wList
	Prompt x1name "Wave1 X2 values", popup, "_none_;"+wList
	Prompt x0, "low X value on first wave"
	Prompt x1, "high X value on first wave"
	Prompt p0, "low point value on first wave"
	Prompt p1, "high point value on first wave"
	if (ask)
		DoPrompt "Select Waves",c0name,c1name,x0name,x1name,x0,x1,p0,p1
		if (V_flag)
			return NaN
		endif
		Wave cc0 = $c0name
		Wave xx0 = $x0name
		Wave cc1 = $c1name
		Wave xx1 = $x1name
		printIt = 1
	elseif (!WaveExists(cc0) || !WaveExists(cc1))
		if (ItemsInList(wList)<2)
			return NaN
		endif
		DoPrompt "Select Waves",c0name,c1name,x0name,x1name
		if (V_flag)
			return NaN
		endif
		Wave cc0 = $c0name
		Wave xx0 = $x0name
		Wave cc1 = $c1name
		Wave xx1 = $x1name
		printIt = 1
	elseif (!WaveExists(cc0))
		if (ItemsInList(wList)<1)
			return NaN
		endif
		DoPrompt "Select Waves",c0name,x0name
		if (V_flag)
			return NaN
		endif
		Wave cc0 = $c0name
		Wave xx0 = $x0name
		printIt = 1
	elseif (!WaveExists(cc1))
		if (ItemsInList(wList)<1)
			return NaN
		endif
		DoPrompt "Select Waves",c1name,x1name
		if (V_flag)
			return NaN
		endif
		Wave cc1 = $c1name
		Wave xx1 = $x1name
		printIt = 1
	endif

	if (printIt)
		c0name = NameOfWave(cc0)
		c1name = NameOfWave(cc1)
		x0name = NameOfWave(xx0)
		x1name = NameOfWave(xx1)
		c0name = SelectString(strlen(c0name),"$\"\"",c0name)
		c1name = SelectString(strlen(c1name),"$\"\"",c1name)
		x0name = SelectString(strlen(x0name),"$\"\"",x0name)
		x1name = SelectString(strlen(x1name),"$\"\"",x1name)
		printf "ScaleFactorBetweenWaves(%s, %s, %s, %s",c0name,x0name,c1name,x1name
		if (numtype(p0)==0)
			printf ", p0=%g",p0
		endif
		if (numtype(p1)==0)
			printf ", p1=%g",p1
		endif
		if (numtype(x0)==0)
			printf ", x0=%g",x0
		endif
		if (numtype(x1)==0)
			printf ", x1=%g",x1
		endif
		printf ")\r"
	endif
	if (!WaveExists(cc0) || !WaveExists(cc1))
		print "ERROR -- ScaleFactorBetweenWaves(), one of the data waves cc0 or cc1 do not exist"
		return NaN
	elseif (WaveExists(xx0) && DimSize(xx0,0)<DimSize(cc0,0))
			print "ERROR -- ScaleFactorBetweenWaves(), xx0 is shorter than cc0"
			return NaN
	elseif (WaveExists(xx1) && DimSize(xx1,0)<DimSize(cc1,0))
			print "ERROR -- ScaleFactorBetweenWaves(), xx1 is shorter than cc1"
			return NaN
	endif

	Variable pMax=DimSize(cc0,0)-1
	if (numtype(p0)==0)
		x0 = WaveExists(xx0) ? xx0[p0] : pnt2x(cc0,p0)
	else
		p0 = numtype(x0) ? 0 : point_from_xw(xx0,cc0,x0,pMax)
		x0 = WaveExists(xx0) ? xx0[p0] : pnt2x(cc0,p0)
	endif
	if (numtype(p1)==0)
		x1 = WaveExists(xx0) ? xx0[p1] : pnt2x(cc0,p1)
	else
		p1 = numtype(x1) ? pMax : point_from_xw(xx0,cc0,x1,pMax)
		x1 = WaveExists(xx0) ? xx0[p1] : pnt2x(cc0,p1)
	endif
	if (p1<p0)
		printf "ERROR -- ScaleFactorBetweenWaves(), p1<p0  (%g < %g)\r",p1,p0
		return NaN
	endif

	// possibly reduce range in case cc1 is shorter
	Variable imax = DimSize(cc1,0)-1
	if (WaveExists(xx1))
		imax = min(imax,DimSize(xx1,0)-1)
		x0 = max(x0,xx1[0])
		x1 = min(x1,xx1[imax])
	else
		x0 = max(x0,leftx(cc1))
		x1 = min(x1,pnt2x(cc1,imax))
	endif

	// recalc p0 & p1, using new (possibly limited x0,x1)
	if (WaveExists(xx0))
		p0 = point_from_xw(xx0,cc0,x0,0)
		p1 = point_from_xw(xx0,cc0,x1,1)
	else
		p0 = limit(x2pnt(cc0,x0),0,pMax)
		p1 = limit(x2pnt(cc0,x1),0,pMax)
	endif
	if (numtype(p0+p1) || p1<p0)
		printf "ERROR -- ScaleFactorBetweenWaves(), p1<p0  (%g < %g)\r",p1,p0
		return NaN
	endif

	// recalc the (x0, x1) that go with the (p0, p1)
	if (WaveExists(xx0))
		x0 = xx0[p0]
		x1 = xx0[p1]
	else
		x0 = pnt2x(cc0,p0)
		x1 = pnt2x(cc0,p1)
	endif

	Make/N=(p1-p0+1)/FREE/D wy0
	wy0 = cc0[p+p0]
	if (WaveExists(xx0))
		Make/N=(p1-p0+1)/FREE/D wx0
		wx0 = xx0[p+p0]
	else
		SetScale/P x, pnt2x(cc0,p0), DimDelta(cc0,0),"", wy0
	endif
	// no longer need p0,p1, process all of wy0, and range of wy0 is [x0,x1]

	Duplicate/FREE wy0, wy1					// create wy1 to match range of wy0
	if (WaveExists(wx0) && WaveExists(xx1))
		wy1 = interp(wx0[p],xx1,cc1)
	elseif (WaveExists(xx1))					// have an xx1, but no wx0, interpolate
		wy1 = interp(x,xx1,cc1)
	elseif (WaveExists(wx0))					// have an wx0, but no xx1, interpolate
		wy1 = cc1(wx0[p])
	else												// have neither wx0 nor xx1, maybe interpolate
		wy1 = cc1(x)
	endif
	Make/N=(numpnts(wy0))/B/U/FREE flags	// use flags to get rid of NaN's or Inf's
	flags = numtype(wy0+wy1) ? 1 : 0
	wy0 = flags ? 0 : wy0
	wy1 = flags ? 0 : wy1

	// I now have a wy1 whose x-values all line up with the wy0
	//
	//	This routine based on min of   SUM{ [Yi - a*f(Xi)]^2 } = V
	// with Yi = wy0[i] and Zi=wy1[i]
	//	Since I have already interpolated, SUM{ [Yi - a*Zi]^2 } = V
	//	dV/da = SUM{ 2*[Yi - a*Zi] * -Zi) }
	//	dV/da = -2 * SUM{ Yi*Zi - a*Zi^2 }
	//	for dV/da = 0 implies
	//	0 = -2 * SUM{ Yi*Zi - a*Zi^2 }
	//	0 = SUM{ Yi*Zi } - a*SUM{ Zi^2 }
	//	SUM{ Yi*Zi } = a*SUM{ Zi^2 }
	//	a = SUM{ Yi*Zi } / SUM{ Zi^2 }
	// so calculate SUM{ Yi*Zi } and SUM{ Zi^2 }
	MatrixOp/FREE yz_zz = sum(wy0*wy1) / sum(wy1*wy1)	// SUM{ Yi*Zi } / SUM{ Zi^2 }

	Variable scale = yz_zz[0]
	if (printIt)
		printf "  the scaling in range [%g, %g] is:       \"%s\" = %g * \"%s\"\r",x0,x1,NameOfWave(cc0),scale,NameOfWave(cc1)
	endif
	return scale
End
//
Static Function point_from_xw(wx,wy,x0,pMax)
	Wave wx										// optional x wave
	Wave wy										// for y-wave scaling (no wx)
	Variable x0
	Variable pMax								// max allowed value of p0

	Variable p0
	if (WaveExists(wx))
		p0 = BinarySearch(wx,x0)				// convert x0 to p0
		p0 = p0==-1 ? 0 : p0
		p0 = p0==-2 ? Inf : p0
	else
		p0 = x2pnt(wy,x0)
	endif
	p0 = limit(p0,0,pMax)
	return p0
End






// added proper handling of NaN's in ywave
Function CombineNearbyPoints(yName,xName,yErr,dx,firstX)
	String yName			// ="Det"
	String xName			// ="dE"
	String yErr				// ="Det_err"
	Variable dx				// max difference between two points considered the same (=0.05)
	Variable firstX			// =-1e30

	Prompt yName,"Name of Y wave",popup, WaveList("*", ";", "")
	Prompt xName,"Name of X wave",popup, WaveList("*", ";", "")
	Prompt yErr,"Name of Y error wave",popup, "_none_;_root_N_;"+WaveList("*", ";", "")
	Prompt dx, "Max difference between two points considered the same"
	Prompt firstX, "X value of where to start Averaging"

	dx = abs(dx)
	if (exists(yName)!=1)
		yName=""
	endif
	if (exists(xName)!=1)
		xName=""
	endif
	if (strlen(yName)<1 || strlen(xName)<1 || dx<=0 || numtype(firstX)==2)
		DoPrompt "Combine nearby points", yName,xName,yErr,dx,firstX
		if (V_flag)
			return 1
		endif
		printf "CombineNearbyPoints(\"%s\", \"%s\", \"%s\", %g, %g)\r",yName,xName,yErr,dx,firstX
	endif
	PauseUpdate

	Variable err = (cmpstr("_none_",yErr)) && (strlen(yErr)>0)		// true when yErr wave given
	Variable root_N=(!cmpstr("_root_N_",yErr))						// used sqrt(N) for error
	if (root_N)
		err = 0
	endif

	if (exists(yName)!=1 || exists(xName)!=1)
		Abort "trying to combine waves that do not exist!"
	endif
	Wave yw=$yName
	Wave xw=$xName
	if (numpnts(yw)!=numpnts(xw))
		Abort "'"+yName+"' and '"+xName+"'  must be the same length"
	endif
	if (err)
		Wave yew=$yErr
		if (numpnts(yw)!=numpnts(yew))
			Abort "'"+yName+"' and '"+yErr+"'  must be the same length"
		endif
	endif
	if (!cmpstr(yName,xName))
		Abort "xwave and ywave must be different"
	endif

	if (err)
		sort xw,xw,yw,yew
	else
		sort xw,xw,yw
	endif

	// remove all NaN's in either xw, yw, or yew
	Variable i = numpnts(xw)-1								// start checking from the end of the waves
	if (err)
		do
			if (numtype(xw[i]+yw[i]+yew[i])==2)// found a NaN
				DeletePoints i, 1, yw,xw,yew
			endif
			i -= 1
		while(i>=0)
	else
		do
			if (numtype(xw[i]+yw[i])==2)		// found a NaN
				DeletePoints i, 1, yw,xw
			endif
			i -= 1
		while(i>=0)
	endif

	Variable sumY,sumX,n,j
	Variable sumE, sigma2
	i=max(0,BinarySearch(xw,firstX))
	do
		if( abs(xw[i+1]-xw[i])<dx )
			n = 0
			do					// make n the number of points to combine it should be at least 2
				n += 1
			while( (abs(xw[i+n]-xw[i])<dx) && ((i+n)<numpnts(xw)))
			sumY = 0
			sumX = 0
			sumE = 0
			sigma2 = 1
			j = 0
			do					// sum up n points and reduce length of xName, yName, and yErr
				if (err)
					sigma2 = (yew[j+i])^2
				endif
				if (root_N)
					sigma2 = abs(yw[j+i])			// sigma = sqrt(y[i])
				endif
				sumX += xw[j+i]
				sumY += yw[j+i] / sigma2
				sumE += 1/sigma2
				j += 1
			while (j<n)

			xw[i] = sumX/n
			xw[(i+1),] = xw[p+n-1]
			yw[i] = sumY / sumE
			yw[(i+1),] = yw[p+n-1]
			if (err)
				yew[i] = sqrt(1/sumE)
				yew[(i+1),] = yew[p+n-1]
				Redimension/N=(numpnts(xw)-n+1) yw,xw,yew
			else
				Redimension/N=(numpnts(xw)-n+1) yw,xw
			endif
		endif

		i += 1
	while (i<(numpnts(xw)-1))
	return 0
End
//Function CombineNearbyPointsOLD(yName,xName,yErr,dx,firstX)
//	String yName			// ="Det"
//	String xName			// ="dE"
//	String yErr				// ="Det_err"
//	Variable dx				// max difference between two points considered the same (=0.05)
//	Variable firstX			// =-1e30
//
//	Prompt yName,"Name of Y wave",popup, WaveList("*", ";", "")
//	Prompt xName,"Name of X wave",popup, WaveList("*", ";", "")
//	Prompt yErr,"Name of Y error wave",popup, "_none_;_root_N_;"+WaveList("*", ";", "")
//	Prompt dx, "Max difference between two points considered the same"
//	Prompt firstX, "X value of where to start Averaging"
//
//	dx = abs(dx)
//	if (exists(yName)!=1)
//		yName=""
//	endif
//	if (exists(xName)!=1)
//		xName=""
//	endif
//	if (strlen(yName)<1 || strlen(xName)<1 || dx<=0 || numtype(firstX))
//		DoPrompt "Combine nearby points", yName,xName,yErr,dx,firstX
//	endif
//	PauseUpdate
//
//	Variable err = (cmpstr("_none_",yErr)) && (strlen(yErr)>0)		// true when yErr wave given
//	Variable root_N=(!cmpstr("_root_N_",yErr))						// used sqrt(N) for error
//	if (root_N)
//		err = 0
//	endif
//
//	if (exists(yName)!=1 || exists(xName)!=1)
//		Abort "trying to combine waves that do not exist!"
//	endif
//	Wave yw=$yName
//	Wave xw=$xName
//	if (numpnts(yw)!=numpnts(xw))
//		Abort "'"+yName+"' and '"+xName+"'  must be the same length"
//	endif
//	if (err)
//		Wave yew=$yErr
//		if (numpnts(yw)!=numpnts(yew))
//			Abort "'"+yName+"' and '"+yErr+"'  must be the same length"
//		endif
//	endif
//	if (!cmpstr(yName,xName))
//		Abort "xwave and ywave must be different"
//	endif
//
//	if (err)
//		sort xw,xw,yw,yew
//	else
//		sort xw,xw,yw
//	endif
//
//	Variable Ndups=0			// number of duplicate sets found
//	Variable sumY,sumX,n,j,i=BinarySearch(xw,firstX)
//	Variable sumE, sigma2
//	i = max(0,i)
////	print "starting at point",i
//	do
//
//		if( abs(xw[i+1]-xw[i])<dx )
//			Ndups += 1
//			n = 0
//			do					// make n the number of points to combine it should be at least 2
//				n += 1
//			while( (abs(xw[i+n]-xw[i])<dx) && ((i+n)<numpnts(xw)))
//			sumY = 0
//			sumX = 0
//			sumE = 0
//			sigma2 = 1
//			j = 0
//			do					// sum up n points and reduce length of xName, yName, and yErr
//				if (err)
//					sigma2 = (yew[j+i])^2
//				endif
//				if (root_N)
//					sigma2 = abs(yw[j+i])			// sigma = sqrt(y[i])
//				endif
//				sumX += xw[j+i]
//				sumY += yw[j+i] / sigma2
//				sumE += 1/sigma2
//				j += 1
//			while (j<n)
//
//			xw[i] = sumX/n
//			xw[(i+1),] = xw[p+n-1]
//			yw[i] = sumY / sumE
//			yw[(i+1),] = yw[p+n-1]
//			if (err)
//				yew[i] = sqrt(1/sumE)
//				yew[(i+1),] = yew[p+n-1]
//				Redimension/N=(numpnts(xw)-n+1) yw,xw,yew
//			else
//				Redimension/N=(numpnts(xw)-n+1) yw,xw
//			endif
//		endif
//
//		i += 1
//	while (i<(numpnts(xw)-1))
////	print "    found ",Ndups," duplicate sets"
//End






// this revision of gcf uses the the "PrimeFactors" command available with Igor 5
// is also assumes that 0 can be divided by anything, so gcf of {0,4,8} is 4
Function gcf(factors)				// find greatest common factor of a wave whose values are all integers
	String factors					// either a list of factors or the name of a wave containing the factors

	String tempName=""
	Variable N,i
	if (exists(factors)==1 && strsearch(factors,";",0)<0)		// passed a wave name
		Wave wav=$factors
	else
		N = ItemsInList(factors)								// a list was passed
		if (N<1)
			return NaN
		endif
		tempName = UniqueName("tempPrimes",1,0)
		Make/N=(N) $tempName
		Wave wav=$tempName
		for (i=0;i<N;i+=1)
			wav[i] = str2num(StringFromList(i,factors))
		endfor
	endif
	N=numpnts(wav)
	if (N<1)
		return NaN
	endif

	String primeListsName = UniqueName("ListOfPrimes",1,0)
	Make/N=(N)/T $primeListsName
	Wave/T primeLists=$primeListsName
	primeLists = ""
	String allPrimes="", item
	Variable j,gcf
	for (i=0;i<N;i+=1)
		PrimeFactors/Q  wav[i]
		Wave W_PrimeFactors=W_PrimeFactors
		for (j=0;j<numpnts(W_PrimeFactors);j+=1)
			primeLists[i]=AddListItem(num2istr(W_PrimeFactors[j]), primeLists[i])
		endfor
		allPrimes+= primeLists[i]
		if (wav[i]==0)									// ensure "0" is there for a 0
			primeLists[i]=AddListItem("0", primeLists[i])
		endif
	endfor

	Variable m,found,jN = ItemsInList(allPrimes)
	j = 0
	gcf=1
	for (j=0;j<jN;j+=1)								// loop over all of the primes
		item = StringFromList(j,allPrimes)
		found = 1
		for (i=0;i<N;i+=1)
			m = WhichListItem(item,primeLists[i])
			if (m>=0)
				primeLists[i] = RemoveListItem(m,primeLists[i])
			elseif (cmpstr(primeLists[i],"0;"))
				found=0
			endif
		endfor
		if (found)
			gcf *= str2num(item)
		endif
	endfor
	KillWaves/Z primeLists,W_PrimeFactors, $tempName
	return gcf
End

