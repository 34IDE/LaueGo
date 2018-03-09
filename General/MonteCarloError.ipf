#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma version = 1.00
#pragma ModuleName=MCerror

Static Constant MCerrorMaxXvals=20,  MCerrorMaxYvals=10

Static Constant tol = 2e-4				// fractional tolerance change in result for stopping


Menu "Macros"
	"Propogate Errors thru Function...", ErrorThruFunction("", $"", $"")
End

//	August, 2017, 0.01, start
//
//	March 8, 2018, 1.00, changed structure and everything else to allow funcName() to return either a scalar or a vector of values
//		Now you can have multiple output values in addition to multiple input values.
//		There are TWO styles of functions:
//		genericRadius(xVals) or  MonteCarloFuncitonProto(xVals)  returns a real scalar
//		genericRadiusW(xVals) or MonteCarloFuncitonProtoW(xVals) returns a vector of reals



//			Error Propogation thru a Function, using a Monte Carlo method
//
// Determine the error in the output of a function based on the errors (standard deviations) of the input values.
// funcName (or func) is a function of Nx variables with each variable having its own standard deviation.
// You can use either ErrorThruFunction() or ErrorFromMonteCarlo() to calculate the standard deviation of func(x-variables)
// arising from the standard deviation of each of the x inputs.
//
// func must have a single argument that is a wave with the x-values, an example is:
//
//	Function genericRadius(xVal)					// computes the distance from origin of xVal in N-dimensional space
//		Wave xVal										// the N input X values
//		return norm(xVal)
//	End
//
//
//	Function/WAVE genericRadiusW(xVal)			// computes the distance from origin of xVal in N-dimensional space, AND angle
//		Wave xVal										// the N input X values
//		Make/N=2/FREE out
//		out[0] = norm(xVal)							// return radius
//		out[1] = atan2(xVal[1],xVal[0]) * 180/PI	// returns angle (degree)
//		return out										// returns a vector
//	End
//
//
//
// two examples of its use is:
//
//	Function test1()									// test using ErrorThruFunction()
//		Make/N=3/D/FREE xWave={3,1,0.2}, xErr={0.01,1,0.1}			// separate wave for x and x-error
//		ErrorThruFunction("genericRadius", xWave, xErr, printIt=1)
//
//		print " "
//		Make/N=(3,2)/D/FREE xValErr = {{3,1,0.2}, {0.01,1,0.1}}	// x and x-error as two columns of one wave
//		ErrorThruFunction("genericRadius", xValErr, $"", printIt=1)
//
//		print " "
//		ErrorThruFunction("genericRadiusW", xValErr, $"", printIt=1)
//	End
//
//
//	Function test2()									// test using ErrorFromMonteCarlo()
//		STRUCT MCerrorStructure MC				// describes the error from a Monte Carlo error analysis
//		MCerror#Init_MCerrorStructure(MC)
//		FUNCREF MonteCarloFuncitonProto func = genericRadius
//		MC.Nx = 3
//		MC.xVal[0] = 3		;	MC.xErr[0] = 0.01
//		MC.xVal[1] = 1		;	MC.xErr[1] = 1
//		MC.xVal[2] = 0.2	;	MC.xErr[2] = 0.1
//		MC.funcName = StringByKey("NAME",FuncRefInfo(func))
//		if (MCerror#ErrorFromMonteCarlo(MC))				// returns value and error of func(xVal,xErr)
//			print "ERROR -- failure calling ErrorFromMonteCarlo()"
//		else
//			print MCerror#MCerrorStructure2str(MC)
//		endif
//	End
//
//
//
//	or you can simply run:
//				ErrorThruFunction("", $"", $"")
//	and you will be prompted for the function name and Xwaves
//
//
//	The MC is an MCerrorStructure that contains all of the details (both inputs and outputs).
//	See the definition of MCerrorStructure below for details.
//
//	The first 4 items in MC are the inputs to the problem {funcName, Nx, xVal[], xErr[]}
//	And the last 6 {N,Ntotal,value,err,avg,sdev} are the result of the Monte Carlo analysis.
//
//	NOTE, avg is likely to be shifted a bit from value if func(...) is not flat at xVal.
//	e.g. consider when func(x) = exp(x), then symmetric errors in x tend to favor "+" excursions over "-" excursions.






Function/S ErrorThruFunction(funcName, xWave, errWave, [MC, Nmax, printIt])
	// returns result of error propogation as a "key=value;" list
	String funcName						// function(vector), returns either a scalar or a vector
	Wave xWave								// x values (or 2 columns with x & err)
	Wave errWave							// errors, if xWave has 2 columns, then errWave can be NULL
	STRUCT MCerrorStructure &MC		// OPTIONAL, will be filled by this function (only an output)
	Variable Nmax							// OPTIONAL, max number of evaluations to use. If not passed, use 1e6, it usually stops first
	Variable printIt
	Nmax = ParamIsDefault(Nmax) || numtype(Nmax) || Nmax<2 ? 1e6 : Nmax
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRtStackInfo(2))==0 : printIt

	if (!WaveExists(xWave) || strlen(FunctionInfo(funcName))<1)
		String wlist1=WaveList("*",";","DIMS:1,CMPLX:0,")
		String wlist2=WaveList("*",";","DIMS:2,MINCols:2,CMPLX:0,")
		String flist=FunctionList("*",";","KIND:2,NPARAMS:1,VALTYPE:9,WIN:Procedure")
		String xName=NameOfWave(xWave), errName=NameOfWave(errWave)
		Prompt funcName, "Function Name", popup, flist
		Prompt xName, "X Values", popup, wlist1+wlist2
		Prompt errName, "Errors of X Values", popup, wlist1+"--none--"
		DoPrompt "Choose for Errors", funcName, xName, errName
		if (V_flag)
			return ""
		endif
		errName = ReplaceString("--none--",errName,"")
		Wave xWave=$xName
		Wave errWave=$errName
		printIt = 1
	endif
	if (printIt)
		xName = SelectString(WaveExists(xWave), "$\"\"", NameOfWave(xWave))
		errName = SelectString(WaveExists(errWave), "$\"\"", NameOfWave(errWave))
		printf "ErrorThruFunction(\"%s\", %s, %s)\r", funcName,xName,errName
	endif
	Variable NX=DimSize(xWave,0)
	if (!(NX>0))
		return ""
	endif

	if (WaveExists(errWave))
		Make/N=(DimSize(errWave,0))/D/FREE errs=NaN
		errs = errWave[p][0]
	elseif (DimSize(xWave,1)>1)
		Make/N=(NX)/D/FREE errs=NaN
		errs = xWave[p][1]
	endif
	Variable Ne=DimSize(errs,0)
	if (Ne>NX)
		Redimension/N=(NX) errs
	elseif (Ne<NX)
		return ""
	endif

	Variable iReturn=NumberByKey("RETURNTYPE",FunctionInfo(funcName)), isProto=1, isVec
	if (iReturn==4)									// Scalar return
		FUNCREF MonteCarloFuncitonProto func = $funcName
		isProto = NumberByKey("ISPROTO",FuncRefInfo(func))
		isVec = 0
	elseif (iReturn==16384)						// Wave return
		FUNCREF MonteCarloFuncitonProtoW funcW = $funcName
		isProto = NumberByKey("ISPROTO",FuncRefInfo(funcW))
		isVec = 1
	endif
	if (isProto)
		return ""
	endif

	STRUCT MCerrorStructure MClocal					// describes the error from a Monte Carlo error analysis
	Init_MCerrorStructure(MClocal)
	MClocal.funcName = StringByKey("NAME",FuncRefInfo(func))+StringByKey("NAME",FuncRefInfo(funcW))

	MClocal.Nmax = Nmax
	MClocal.NX = NX
	MClocal.NY = 1											// don't yet know length of vector returned, but must have at least 1
	Variable i
	for (i=0;i<NX;i+=1)
		MClocal.xVal[i] = xWave[i][0]
		MClocal.xErr[i] = errs[i]
	endfor

	Variable err = ErrorFromMonteCarlo(MClocal)// get errors for func(xVal,xErr)
	if (!ParamIsDefault(MC))							// MC was passed, fill it
		MC = MClocal
	endif

	if (printIt && err)
		print "ERROR -- failure calling ErrorFromMonteCarlo()"
	elseif (printIt)
			print MCerrorStructure2str(MClocal)
	endif

	return MCerror#MCerrorStructure2KeyVals(MClocal)
End


Static Function ErrorFromMonteCarlo(MC)		// fills MC structure with info about errors of func(xVal), loops until avg is stable
	STRUCT MCerrorStructure &MC					// describes the error from a Monte Carlo error analysis

	Variable tick=stopMSTimer(-2)
	Make/N=(MC.NX)/D/FREE xVal = MC.xVal[p]
	Make/N=(MC.NX)/D/FREE xErr = MC.xErr[p]
	Variable iReturn=NumberByKey("RETURNTYPE",FunctionInfo(MC.funcName)), isProto=1, isVec
	if (iReturn==4)									// Scalar return
		FUNCREF MonteCarloFuncitonProto func = $MC.funcName
		isProto = NumberByKey("ISPROTO",FuncRefInfo(func))
		isVec = 0
	elseif (iReturn==16384)						// Wave return
		FUNCREF MonteCarloFuncitonProtoW funcW = $MC.funcName
		isProto = NumberByKey("ISPROTO",FuncRefInfo(funcW))
		isVec = 1
	endif
	if (isProto)
		return 1
	endif
	Variable vTemp

	Duplicate/FREE xVal, xs
	xs = xVal
	if (isVec)
		Wave value = funcW(xs)						// set the target values, which is funcW(xVal)
		MC.NY = numpnts(value)
		if (MC.NY > MCerrorMaxYvals)
			printf "ERROR -- ErrorFromMonteCarlo(),  %s() returns %d values, can only accept %d (change MCerrorMaxYvals)\r", MC.funcName, MC.NY, MCerrorMaxYvals
			return 1
		endif
	else
		vTemp  = func(xs)								// set the target values, which is func(xVal)
		Make/N=1/D/FREE ys1, value1
		Wave value=value1, ys=ys1
		value = vTemp
		MC.NY = 1
	endif
	if (numtype(sum(value)))
		return 1
	endif

	Make/N=(MC.NY)/D/FREE sumY2=0, sumY=0, avg=NaN, sdev, diff, xbar, changing
	Variable i, N, Ntotal
	for (i=0,N=0,Ntotal=0; i< MC.Nmax; i+=1)// loop until errors are good enough, or reach Nmax evaluations
		xs = xVal[p] + gnoise(xErr[p])			// next set of x values with errors added
		if (isVec)
			Wave ys = funcW(xs)						// value of funcW with this set of x's
		else
			vTemp = func(xs)							// value of func with this set of x's
			ys = vTemp
		endif
		Ntotal += 1										// total number of function calls (regardless of wheter they are valid)
		if (numtype(sum(ys))==0)
			sumY2 += ys*ys								// accumumlate valid values for statistics
			sumY += ys
			N += 1										// this is a valid result (ys is not Inf or NaN)
		endif

		if (mod(i,50)==1)								// check avg and sdev every 50 evaluations
			xbar = sumY/N								// Note, xbar, diff, avg, sdev are VECTOR equations (implicit loops)
			diff = xbar - avg							// difference between current_avg and last_avg
			avg = xbar									// update current values of avg & sdev
			sdev = sqrt( (sumY2 - 2*xbar*sumY + N*xbar*xbar)/(N-1) )
			changing = abs(diff/sdev)
			if (N>50 && WaveMax(changing)<tol)	// done when avg is changing by less than tol*sdev
				break
			endif
		endif
	endfor

	// assign results
	MC.N = N												// number of VALID evaluations of function
	MC.Ntotal = Ntotal								// total number of times function was called
	for (i=0;i<MC.NY;i+=1)
		MC.value[i] = value[i]						// set the target value and standard deviation about target value
		MC.err[i] = sqrt( (sumY2[i] - 2*value[i]*sumY[i] + N*value[i]*value[i])/(N-1) )
		MC.avg[i] = avg[i]							// average value of function (using only valid values), NOT value
		MC.sdev[i] = sdev[i]						// standard deviation using N evaluations
	endfor
	MC.seconds = (stopMSTimer(-2)-tick)*1e-6

	return (N<3)
End
//
Function MonteCarloFuncitonProto(xVal)		// function(vector) --> scalar
	Wave xVal											// the input X values
	return NaN
End
//
Function/WAVE MonteCarloFuncitonProtoW(xVal)	// function(vector) --> vector (a wave with Y values)
	Wave xVal											// the input X values
	return $""											// should contain wave with Y values
End


Structure MCerrorStructure			// describes the error from a Monte Carlo error analysis
	// These first 4 are filled in by the user
	char		funcName[50]				// name of function used to evaluate, can return either scalar or wave
	int32		NX								// number of X values (never more than MCerrorMaxXvals)
	int32		NY								// number of Y values (never more than MCerrorMaxYvals)
	double	xVal[MCerrorMaxXvals]	// X values for input to func
	double	xErr[MCerrorMaxXvals]	// errors for each of the X values (standard deviation)
	double	Nmax							// max number of evaluations to try. It usually stops first
	//
	// The following are computed by ErrorFromMonteCarlo()
	int32		N								// number of VALID evaluations of function
	int32		Ntotal						// total number of evaluations of function
	double	value[MCerrorMaxYvals]	// value of func(xVal)
	double	err[MCerrorMaxYvals]	// standard deviation of N evaluations about value (NOT average)
	double	avg[MCerrorMaxYvals]	// average value of function (this may be different than value)
	double	sdev[MCerrorMaxYvals]	// standard deviation of N evaluations about average
	double	seconds						// execution time (second)
	// NOTE, avg is likely to be shifted a bit from value if func(...) is not flat at xVal.
	//		Consider the case of func(x) = exp(x), then symmetric errors in x tend to favor + excursions over - ones.
EndStructure


Static Function/S Init_MCerrorStructure(MC)
	STRUCT MCerrorStructure &MC
	MC.N = 0				;		MC.Ntotal = 0
	MC.funcName = ""
	MC.NX = 0			;		MC.NY = 0
	MC.Nmax = 0			;		MC.seconds = 0
	Variable i
	for (i=0;i<MCerrorMaxXvals;i+=1)
		MC.xVal[i] = NaN		;		MC.xErr[i] = NaN
	endfor
	for (i=0;i<MCerrorMaxYvals;i+=1)
		MC.value[i] = NaN		;		MC.err[i] = NaN
		MC.avg[i] = NaN		;		MC.sdev[i] = NaN
	endfor
End
//
Static Function/S MCerrorStructure2str(MC)
	STRUCT MCerrorStructure &MC					// describes the error from a Monte Carlo error analysis
	String str, valStr="", out=MC.funcName+"("
	Variable i, spx=(MC.NX > 1 ? 0 : 1), spy=(MC.NY > 1 ? 0 : 1)
	for (i=0;i<MC.NX;i+=1)
		out += ValErrStr(MC.xVal[i],MC.xErr[i],sp=spx) + ", "
	endfor
	out += ")"
	out = ReplaceString(", )",out,")")
	for (i=0;i<MC.NY;i+=1)
		valStr += ValErrStr(MC.value[i], MC.err[i], sp=spy)+", "
	endfor
	valStr = TrimEnd(valStr, chars=", ")
	valStr = SelectString(MC.NY>1, valStr, "{"+valStr+"}")
	String maxString = SelectString(MC.Ntotal < MC.Nmax, ", reached Nmax","")
	if (MC.N == MC.Ntotal)
		sprintf str, " = %s,    [from %d evaluations%s]", valStr, MC.N, maxString
	else
		sprintf str, " = %s,    [from %d valid evaluations out of %d calls%s]", valStr, MC.N, MC.Ntotal, maxString
	endif

	out += str

	Make/N=(MC.NY)/D/FREE maxDiff = abs(MC.value[p] - MC.avg[p]) > (MC.err[i])/10
	if (sum(maxDiff))									// one of the .avg and .value don't match, so show them too
		for (i=0;i<MC.NY;i+=1)
			sprintf str, ",  <y> = %s",ValErrStr(MC.avg[i], MC.sdev[i], sp=spy)
			out += str
		endfor
	endif
	sprintf str, "\r\t\ttook %.3g sec", MC.seconds
	out += str
	return out
End
//
Static Function/S MCerrorStructure2KeyVals(MC)
	STRUCT MCerrorStructure &MC					// describes the error from a Monte Carlo error analysis
	String keyVal=""
	if (strlen(MC.funcName))
		keyVal = ReplaceStringByKey("funcName",keyVal,MC.funcName,"=")
	endif
	Variable NX=limit(MC.NX,0,MCerrorMaxXvals), NY=limit(MC.NY,0,MCerrorMaxYvals)
	keyVal = ReplaceNumberByKey("NX",keyVal,NX,"=")
	keyVal = ReplaceNumberByKey("NY",keyVal,NY,"=")
	keyVal = ReplaceNumberByKey("Nmax",keyVal,MC.Nmax,"=")
	keyVal = ReplaceNumberByKey("seconds",keyVal,MC.seconds,"=")
	if (NX>0)
		Make/N=(NX)/D/FREE xVal=MC.Xval[p], xErr=MC.xErr[p]
		keyVal = ReplaceStringByKey("xVal",keyVal,vec2str(xVal,places=15,bare=1,sep=","),"=")
		keyVal = ReplaceStringByKey("xErr",keyVal,vec2str(xErr,places=15,bare=1,sep=","),"=")
		keyVal = ReplaceNumberByKey("N",keyVal,MC.N,"=")
		keyVal = ReplaceNumberByKey("Ntotal",keyVal,MC.Ntotal,"=")
	endif
	if (NY>0)
		Make/N=(NY)/D/FREE value=MC.value[p], err=MC.err[p], avg=MC.avg[p], sdev=MC.sdev[p]
		keyVal = ReplaceStringByKey("value",keyVal,vec2str(value,places=15,bare=1,sep=","),"=")
		keyVal = ReplaceStringByKey("err",keyVal,vec2str(err,places=15,bare=1,sep=","),"=")
		keyVal = ReplaceStringByKey("avg",keyVal,vec2str(avg,places=15,bare=1,sep=","),"=")
		keyVal = ReplaceStringByKey("sdev",keyVal,vec2str(sdev,places=15,bare=1,sep=","),"=")
		keyVal = ReplaceStringByKey("value",keyVal,vec2str(value,places=15,bare=1,sep=","),"=")
	endif
	return keyVal
End
