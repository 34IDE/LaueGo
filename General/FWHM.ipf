#pragma rtGlobals=1		// Use modern global access method.
#pragma version = 1.4


Menu "Analysis"
	"-"
	"FWHM of wave", PlainFWHMofWave("","")
	"FWHM Gaussian",  print FWHM_Gaussian()
	"FWHM Lorentzian",  print FWHM_Lorentz()
End


Function FWHM_Gaussian()
	Variable/G V_FWHM
	Wave W_coef=W_coef
	if (!WaveExists(W_coef))
		V_FWHM = NaN
		return NaN
	endif
	V_FWHM = 2*W_coef[3]*sqrt(ln(2))
	return V_FWHM
End

Function FWHM_Lorentz()
	Variable/G V_FWHM
	Wave W_coef=W_coef
	if (!WaveExists(W_coef))
		V_FWHM = NaN
		return NaN
	endif
	V_FWHM = 2*sqrt(W_coef[3])
	return V_FWHM
End


Function FWHM_peak(ww)
	Wave ww
	WaveStats/Q ww
	Make/N=2/FREE/D ww_temp__
	FindLevels /D=ww_temp__ /N=6 /Q ww, ((V_max+V_min)/2)
	if (V_LevelsFound<2)
		return NaN					// no peak
	endif

	Variable p1,p2				// points values of the half widths
	Variable x1,x2				// x values of halfwidths
	p1 = BinarySearch(ww_temp__, V_maxloc)
	x2 = ww_temp__[p1+1]
	x1 = ww_temp__[p1]
	Variable/G V_center=(x2+x1)/2
	return abs(x2-x1)
End


Function FWHM_peakXY(wy,wx)
	Wave wy,wx

	WaveStats/Q wy
	Make/N=2/FREE/D ww_temp__
	FindLevels /D=ww_temp__ /N=6 /P/Q wy, ((V_max+V_min)/2)
	if (V_LevelsFound<2)
		return NaN					// no peak
	endif

	Variable pm,p1,p2, dx			// points of mid, HW point, and width
	pm = x2pnt(wy,V_maxloc)		// point of the max
	p1 = BinarySearch(ww_temp__, pm)
	p2 = ww_temp__[p1+1]
	p1 = ww_temp__[p1]
	dx = abs(wx[p2]-wx[p1])
	Variable/G V_center=(wx[p2]+wx[p1])/2
	return dx
End


Function PlainFWHMofWave(yw,xw)
	String yw, xw
	if (stringmatch(xw,"none"))
		xw = "_calculated_"
	endif
	if (exists(yw)!=1 || (exists(xw)!=1 && !stringmatch(xw,"_calculated_") ))
		Prompt yw, "input wave", popup, WaveList("*",";","DIMS:1")
		Prompt xw, "input Xwave", popup, "_calculated_;"+WaveList("*",";","DIMS:1")
		DoPrompt "wave(s) defining peak", yw,xw
	endif
	if (V_flag)
		Abort
	endif

	Variable/G V_FWHM
	if (cmpstr("_calculated_",xw)==0)
		V_FWHM=FWHM_peak($yw)
	else
		V_FWHM=FWHM_peakXY($yw,$xw)
	endif
	NVAR V_center=V_center
	if (ItemsInList(GetRTStackInfo(0))<2)
		printf "PlainFWHMofWave(\"%s\",\"%s\") = %g, centered at %g\r",yw,xw,V_FWHM,V_center
	endif
	return V_FWHM
End


