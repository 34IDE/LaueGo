#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma version = 0.10
#pragma ModuleName=MCerror

Static Constant MCerrorMaxXvals = 20

Menu "Macros"
	"Propogate Errors thru Function...", ErrorThruFunction("", $"", $"")
End

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
//
// two examples of its use is:
//
//	Function test1()									// test using ErrorThruFunction()
//		Make/N=3/D/FREE xWave={3,1,0.2}, xErr={0.01,1,0.1}			// separate wave for x and x-error
//		ErrorThruFunction("genericRadius", xWave, xErr, printIt=1)
//
//		Make/N=(3,2)/D/FREE xValErr = {{3,1,0.2}, {0.01,1,0.1}}	// x and x-error as two columns of one wave
//		ErrorThruFunction("genericRadius", xValErr, $"", printIt=1)
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
//	e.g. consider when func(x) = exp(x), then symmetric errors in x tend to favor + excursions over - ones.






Function/S ErrorThruFunction(funcName, xWave, errWave, [MC, printIt])
	// returns result of error propogation as a "key=value;" list
	String funcName
	Wave xWave								// x values (or 2 columns with x & err)
	Wave errWave							// errors, if xWave has 2 columns, then errWave can be NULL
	STRUCT MCerrorStructure &MC		// OPTIONAL, will be filled by this function (only an output)
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRtStackInfo(2))==0 : printIt

	if (!WaveExists(xWave) || strlen(FunctionInfo(funcName))<1)
		String wlist1=WaveList("*",";","DIMS:1,CMPLX:0,")
		String wlist2=WaveList("*",";","DIMS:2,MINCols:2,CMPLX:0,")
		String flist=FunctionList("*",";","KIND:2,NPARAMS:1,VALTYPE:1,WIN:Procedure")
		flist += FunctionList("*",";","KIND:2,NPARAMS:1,VALTYPE:1,WIN:MonteCarloError")
		flist = RemoveFromList("MonteCarloFuncitonProto",flist)
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
	Variable Nx=DimSize(xWave,0)
	if (!(Nx>0))
		return ""
	endif

	if (WaveExists(errWave))
		Make/N=(DimSize(errWave,0))/D/FREE errs=NaN
		errs = errWave[p][0]
	elseif (DimSize(xWave,1)>1)
		Make/N=(Nx)/D/FREE errs=NaN
		errs = xWave[p][1]
	endif
	Variable Ne=DimSize(errs,0)
	if (Ne>Nx)
		Redimension/N=(Nx) errs
	elseif (Ne<Nx)
		return ""
	endif

	FUNCREF MonteCarloFuncitonProto func = $funcName
	if (NumberByKey("ISPROTO",FuncRefInfo(func)))
		return ""
	endif
	STRUCT MCerrorStructure MClocal					// describes the error from a Monte Carlo error analysis
	Init_MCerrorStructure(MClocal)
	MClocal.funcName = StringByKey("NAME",FuncRefInfo(func))

	MClocal.Nx = Nx
	Variable i
	for (i=0;i<Nx;i+=1)
		MClocal.xVal[i] = xWave[i][0]
		MClocal.xErr[i] = errs[i]
	endfor

	String out = "ERROR -- failure calling ErrorFromMonteCarlo()"
	if (!ErrorFromMonteCarlo(MClocal))				// returns value and error of func(xVal,xErr)
		out = MCerrorStructure2str(MClocal)
	endif
	if (!ParamIsDefault(MC))
		MC = MClocal
	endif

	if (printIt)
		print out
	endif
	return MCerror#MCerrorStructure2KeyVals(MClocal)
End


Static Function ErrorFromMonteCarlo(MC,[Nmax])// fills MC structure with info about errors of func(xVal), loops until avg is stable
	STRUCT MCerrorStructure &MC					// describes the error from a Monte Carlo error analysis
	Variable Nmax										// max number of evaluations to use. If not passed, use 1e6, it usually stops first
	Nmax = ParamIsDefault(Nmax) || numtype(Nmax) || Nmax<1 ? 1e6 : Nmax

	Make/N=(MC.Nx)/D/FREE xVal = MC.xVal[p]
	Make/N=(MC.Nx)/D/FREE xErr = MC.xErr[p]
	FUNCREF MonteCarloFuncitonProto func = $MC.funcName
	if (NumberByKey("ISPROTO",FuncRefInfo(func)))
		return 1
	endif

	Duplicate/FREE xVal, xs
	xs = xVal
	Variable value = func(xs)						// set the target value, which is func(xVal)
	if (numtype(value))
		return 1
	endif

	Variable sumY=0, sumY2=0, N=0, Ntotal=0
	Variable i, avg, sdev, diff, ys, xbar
	for (i=0,avg=NaN; i<Nmax; i+=1)
		xs = xVal[p] + gnoise(xErr[p])			// next set of x values with errors added
		ys = func(xs)									// value of func with this set of x's
		Ntotal += 1										// total number of function calls (regardless of wheter they are valid)
		if (numtype(ys)==0)
			sumY2 += ys*ys								// accumumlate for statistics
			sumY += ys
			N += 1										// this is a valid result (ys is not Inf or NaN)
		endif

		if (mod(i,50)==1)								// check avg and sdev every 50 evaluations
			xbar = sumY/N
			diff = abs(avg - xbar)
			avg = xbar									// update current values of avg & sdev
			sdev = sqrt( (sumY2 - 2*xbar*sumY + N*xbar*xbar)/(N-1) )
			if (N>50 && diff < abs(sdev*1e-5))	// done if avg is changing by less than sdev/1e5
				break
			endif
		endif
	endfor

	// assign results
	MC.value = value									// set the target value and standard deviation about target value
	MC.err = sqrt( (sumY2 - 2*value*sumY + N*value*value)/(N-1) )
	MC.N = N												// number of VALID evaluations of function
	MC.Ntotal = Ntotal								// total number of times function was called
	MC.avg = avg										// average value of function (using only valid values), NOT value
	MC.sdev = sdev										// standard deviation using N evaluations
	return (N<3)
End
//
Function MonteCarloFuncitonProto(xVal)
	Wave xVal											// the input X values
	return NaN
End


Structure MCerrorStructure			// describes the error from a Monte Carlo error analysis
	// These first 4 are filled in by the user
	char		funcName[50]				// name of function used to evaluate
	int32		Nx								// number of X values (never more than MCerrorMaxXvals)
	double	xVal[MCerrorMaxXvals]	// X values for input to func
	double	xErr[MCerrorMaxXvals]	// errors for each of the X values (standard deviation)
	//
	// The following are computed by ErrorFromMonteCarlo()
	int32		N								// number of VALID evaluations of function
	int32		Ntotal						// total number of evaluations of function
	double	value							// value of func(xVal)
	double	err							// standard deviation of N evaluations about value (NOT average)
	double	avg							// average value of function (this may be different than value)
	double	sdev							// standard deviation of N evaluations about average
	// NOTE, avg is likely to be shifted a bit from value if func(...) is not flat at xVal.
	//		Consider the case of func(x) = exp(x), then symmetric errors in x tend to favor + excursions over - ones.
EndStructure
//

Static Function/S Init_MCerrorStructure(MC)
	STRUCT MCerrorStructure &MC					// describes the error from a Monte Carlo error analysis
	MC.N = 0				;		MC.Ntotal = 0
	MC.value = NaN		;		MC.err = NaN
	MC.avg = NaN		;		MC.sdev = NaN
	MC.funcName = ""
	MC.Nx = 0
	Variable i
	for (i=0;i<MCerrorMaxXvals;i+=1)
		MC.xVal[i] = NaN
		MC.xErr[i] = NaN
	endfor
End
//
Static Function/S MCerrorStructure2str(MC)
	STRUCT MCerrorStructure &MC					// describes the error from a Monte Carlo error analysis
	String str, out = MC.funcName + "("
	Variable i
	for (i=0;i<MC.Nx;i+=1)
		out += ValErrStr(MC.xVal[i],MC.xErr[i]) + ", "
	endfor
	out += ")"
	out = ReplaceString(", )",out,")")
	if (MC.N == MC.Ntotal)
		sprintf str, " = %s,    [from %d valid evaluations of %s()]", ValErrStr(MC.value,MC.err,sp=1), MC.N, MC.funcName
	else
		sprintf str, " = %s,    [from %d valid evaluations out of %d calls to %s()]", ValErrStr(MC.value,MC.err,sp=1), MC.N, MC.Ntotal, MC.funcName
	endif
	out += str
	if (abs(MC.value - MC.avg)> (MC.err)/10)	// .avg and .value don't match, so show them too
		sprintf str, ",  <y> = %s",ValErrStr(MC.avg,MC.sdev)
		out += str
	endif
	return out
End
//
Static Function/S MCerrorStructure2KeyVals(MC)
	STRUCT MCerrorStructure &MC					// describes the error from a Monte Carlo error analysis
	String keyVal=""
	if (strlen(MC.funcName))
		keyVal = ReplaceStringByKey("funcName",keyVal,MC.funcName,"=")
	endif
	Variable Nx=limit(MC.Nx,0,MCerrorMaxXvals)
	keyVal = ReplaceNumberByKey("Nx",keyVal,MC.Nx,"=")
	if (Nx>0)
		Make/N=(NX)/D/FREE xVal=MC.Xval[p], xErr=MC.xErr[p]
		keyVal = ReplaceStringByKey("xVal",keyVal,vec2str(xVal,places=15,bare=1,sep=","),"=")
		keyVal = ReplaceStringByKey("xErr",keyVal,vec2str(xErr,places=15,bare=1,sep=","),"=")
		keyVal = ReplaceNumberByKey("N",keyVal,MC.N,"=")
		keyVal = ReplaceNumberByKey("Ntotal",keyVal,MC.Ntotal,"=")
		keyVal = ReplaceNumberByKey("value",keyVal,MC.value,"=")
		keyVal = ReplaceNumberByKey("err",keyVal,MC.err,"=")
		keyVal = ReplaceNumberByKey("avg",keyVal,MC.avg,"=")
		keyVal = ReplaceNumberByKey("sdev",keyVal,MC.sdev,"=")
	endif
	return keyVal
End
