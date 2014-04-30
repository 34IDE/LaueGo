#pragma rtGlobals=1		// Use modern global access method.
#pragma version = 1.41
#pragma ModuleName=Refinement
#include "Indexing", version>=1.3
#include "Triangulate", version>=2.1

Static Constant hc = 1.239841857			// keV-nm

Menu "micro"
	SubMenu "Refinement"
		"Fill Refinement Data from Indexing", FillGeoRefineDataFromIndexing($"",$"")
		"Edit Data for Refinement...", EditGeoRefinementData(NaN)
		"Refine Geometry...",MakeRefinePanel()
		"Refine Lattice --> get strain",doLatticeRefinement($"")
	End
	SubMenu("Tables")
		"Refinement Data", EditGeoRefinementData(NaN)
		help={"Show/edit the data used for Refinement"}
	End
End




Function MeasuredToIndexedField(scale)				// draw vector distortion arrows from measured to indexed spots
	Variable scale									// scale factor for arrow length (about 200)
	scale = numtype(scale) ? 200 : scale
	Wave RefineData=RefineData
	if (!WaveExists(RefineData))
		Abort "cannot find 'RefineData'"
	endif
	STRUCT microGeometry geo
	if (FillGeometryStructDefault(geo))
		Abort "need to read the geometry"
	endif
	STRUCT crystalLattice xtal
	if (FillLatticeStructDefault(xtal))
		Abort "need to read an xtal"
	endif

	printf "Draw vectors from measured spots to calculated (i.e. indexed) spots, exagerated by x %g\r",scale
	Wave qvec = $microGeo#MakeUnique3Vector($"")
	Variable/C pz
	Variable px,py,h,k,l
	Variable dx, dy
	Variable i, N=DimSize(RefineData,0)
	printf "  hkl\r"
	for (i=0;i<N;i+=1)
		if (hklStr2hkl(GetDimLabel(RefineData,0,i),h,k,l))
			continue							// skip this, no hkl entered for this point
		endif
		px = RefineData[i][3]					// actual measured spot position
		py = RefineData[i][4]
		pixel2q(geo,px,py,$"")
		qvec = RefineData[i][p]
		pz = q2pixel(geo,qvec)
		dx = (real(pz) - px) * scale
		dy = (imag(pz) - py) * scale
		printf "(%s)		(%.3f,  \t%.3f)		(%.3f,  \t%.3f)		Æ=(%.3f,  %.3f)\r",hkl2str(h,k,l),px,py,real(pz),imag(pz),dx/scale,dy/scale
		// draw a line from measured spot to calculated (i.e. indexed) spot
		DrawArrow("",px,py,px+dx,py+dy)
	endfor
	KillWaves/Z qvec,M_Lower,M_Upper,W_LUPermutation,M_x
End
Static Function DrawArrow(win,x1,y1,x2,y2)
	String win
	Variable x1,y1,x2,y2
	SetDrawLayer/W=$win UserFront
	SetDrawEnv/W=$win xcoord=bottom,ycoord=left,arrow=1
	DrawLine/W=$win x1,y1, x2,y2
End





// This actually optimizes the reciprocal lattice, from which the lattice is extracted
Function doLatticeRefinement(RefineData,[useLast,addNoise])
	Wave RefineData
	Variable useLast
	Variable addNoise
	useLast = numtype(useLast) || !useLast ? 0 : 1	// ensure valid logical value
	addNoise = numtype(addNoise) || !addNoise ? 0 : 1	// ensure valid logical value
	if (ParamIsDefault(useLast))
		useLast = 0
	endif
	if (ParamIsDefault(addNoise))
		addNoise = 0
	endif
	if (!WaveExists(RefineData))
		if (stringmatch(WaveListClass("RefinementData","*",""),"RefineData;"))
			Wave RefineData = :RefineData		// RefineData exists, and it is the only one
		else
			Variable startPoint					// 1=usual, 2=last, 3=usual+noise
			startPoint = addNoise ? 3 : 1
			startPoint = useLast ? 2 : startPoint
			String RefineDataStr = SelectString(WaveExists(RefineData),"",NameOfWave(RefineData))
			Prompt RefineDataStr,"Indexed Peak List",popup,WaveListClass("RefinementData","*","")
			Prompt startPoint, "starting point of refinement",popup,"Start Fresh;Start with Last Values;Start Fresh with Noise"
			DoPrompt/Help="3D-Xray Diffraction[Indexing]" "Refinement data",RefineDataStr,startPoint
			if (V_flag)
				return 1
			endif
			Wave RefineData=$RefineDataStr
			addNoise = startPoint==3
			useLast = startPoint==2
		endif
		if (!WaveExists(RefineData))
			return 1
		endif
		printf "¥doLatticeRefinement(%s,useLast=%d,addNoise=%d)\r",NameOfWave(RefineData),useLast,addNoise
		printf "starting refinement using data from '%s'\r",NameOfWave(RefineData)
	endif

	String wnote = note(RefineData)
	Variable as0,as1,as2						// componenets of a*
	Variable bs0,bs1,bs2						// componenets of b*
	Variable cs0,cs1,cs2						// componenets of c*
	String recip_lattice = StringByKey("recip_lattice",wnote,"=")	// the starting point
	sscanf recip_lattice, "{{%g,%g,%g}{%g,%g,%g}{%g,%g,%g}}",as0,as1,as2,bs0,bs1,bs2,cs0,cs1,cs2
	if (V_flag!=9)
		DoAlert 0, "Unable to read recip_lattice"
		return 1
	endif

	// change recip lattice vectors from ideal system into Wenge system
	STRUCT microGeometry geo
	if (FillGeometryStructDefault(geo))
		Doalert 0, "cannot find a geometry structure"
		return 1
	endif
	Make/N=3/O/D root:Packages:geometry:totalError_ki,root:Packages:geometry:totalError_qhat
	Make/N=3/O/D root:Packages:geometry:totalError_qvec
	Wave qvec=root:Packages:geometry:totalError_qvec		// a Q vector, not normalized
	qvec = {as0, as1, as2}
	IdealBeamLine2wenge(geo,qvec,qvec)
	as0=qvec[0]  ;	as1=qvec[1]  ;	as2=qvec[2]
	qvec = {bs0, bs1, bs2}
	IdealBeamLine2wenge(geo,qvec,qvec)
	bs0=qvec[0]  ;	bs1=qvec[1]  ;	bs2=qvec[2]
	qvec = {cs0, cs1, cs2}
	IdealBeamLine2wenge(geo,qvec,qvec)
	cs0=qvec[0]  ;	cs1=qvec[1]  ;	cs2=qvec[2]
	STRUCT crystalLattice xtal0				// starting, or undeformed lattice
	LatticeFromRecip(xtal0, as0,as1,as2,  bs0,bs1,bs2,  cs0,cs1,cs2)	// find lattice constants from recip
	printf "before optimizing, undeformed lattice: \r"
	printf "		a = %.7f		b = %.7f		c = %.7f (nm)\r",xtal0.a,xtal0.b,xtal0.c
	printf "		alpha = %.5f	beta = %.5f	gam = %.5f¡\r",xtal0.alpha,xtal0.beta,xtal0.gam
	Variable i, N, haveEnergy=0
	for (i=0;i<DimSize(RefineData,0);i+=1)// see if there is any energy data entered
		if (numtype(RefineData[i][7]+RefineData[i][6])==0)
			haveEnergy = 1
			break
		endif
	endfor
	String funcName 						// name of function to fit
	if (haveEnergy)
		N = 9								// number of parameters to optimize
		funcName = "totalErrorLattice"
	else
		N = 8								// number of parameters to optimize (only deviatoric strain)
		funcName = "totalErrorLatticeDev"
		Variable Vsc = (as1*bs2-as2*bs1)*cs0 + (as2*bs0-as0*bs2)*cs1 + (as0*bs1-as1*bs0)*cs2
		wnote = ReplaceNumberByKey("Vsc",wnote,Vsc,"=")
		Note/K RefineData, wnote			// need to pass length of c* when only doing deviatoric strain
	endif
	if (exists(funcName)!=6)
		DoAlert 0,"bad choice of optimizing parameters, no such function '"+funcName+"'"
		return 1
	endif
	if (strlen(FindTableWithWave(RefineData)))
		DoAlert 0, "You must close the table with Refinement Data before you can refine"
		return 1
	endif

	if (useLast)								// either add some noise, or use last starting point
		String lastValues = StrVarOrDefault(":LastOptimizeValues","")
		if (strlen(lastValues))
			sscanf lastValues, "%g,%g,%g,%g,%g,%g,%g,%g,%g",as0,as1,as2,bs0,bs1,bs2,cs0,cs1,cs2
		endif
	elseif(addNoise)
		Variable len = sqrt(as0^2 + as1^2 + as2^2) * 5e-3		// add extra strain to starting point
		as0 += enoise(len);		as1 += enoise(len);		as2 += enoise(len)
		bs0 += enoise(len);	bs1 += enoise(len);	bs2 += enoise(len)
		cs0 += enoise(len);		cs1 += enoise(len);		cs2 += enoise(len)
	endif
	Make/N=(N)/O/D RefineLat_typXWave, RefineLat_xWave
	Wave typXWave=RefineLat_typXWave, xWave=RefineLat_xWave
	if (haveEnergy)
		xWave = {as0,as1,as2,bs0,bs1,bs2,cs0,cs1,cs2}
	else
		xWave = {as0,as1,as2,bs0,bs1,bs2,cs0,cs1}
	endif
	WaveStats/M=1/Q xWave
	if (V_numNans+V_numINFs)			// bad data points
		KillWaves/Z RefineLat_typXWave, RefineLat_xWave
		return 1
	endif
	typXWave = xWave

	Variable ystart = totalErrorLattice(RefineData, as0,as1,as2, bs0,bs1,bs2, cs0,cs1,cs2)
	printf "starting optimization with xWave = %s,  with an error of %g\r", vec2str(xWave,places=7), ystart
	if (numtype(ystart))
		DoAlert 0, "totalErrorLattice() returns NaN at start"
		return 1
	endif

	Variable timer = startMSTimer
	if (haveEnergy)
		Optimize/M={2,1}/R=typXWave/X=xWave/Y=(ystart) $funcName,RefineData
	else
		Optimize/M={0,1}/R=typXWave/X=xWave/Y=(ystart) $funcName,RefineData
	endif
	Variable seconds = stopMSTimer(timer)/1e6
	as0 = xWave[0]  ;	as1 = xWave[1]  ;	as2 = xWave[2]		// optimized values
	bs0 = xWave[3]  ;	bs1 = xWave[4]  ;	bs2 = xWave[5]
	cs0 = xWave[6]  ;	cs1 = xWave[7]
	if (haveEnergy)
		cs2 = xWave[8]
	else
	//		Vsc = ((a*)x(b*))¥c
	//		(a*) x (b*) = {(as1*bs2-as2*bs1), (as2*bs0-as0*bs2), (as0*bs1-as1*bs0)}
	//		Vsc = (as1*bs2-as2*bs1)*cs0 + (as2*bs0-as0*bs2)*cs1 + (as0*bs1-as1*bs0)*cs2
		cs2 = (Vsc - (as1*bs2-as2*bs1)*cs0 - (as2*bs0-as0*bs2)*cs1) / (as0*bs1-as1*bs0)
	endif
	String/G :LastOptimizeValues
	SVAR LastOptimizeValues=:LastOptimizeValues
	sprintf LastOptimizeValues,"%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g",as0,as1,as2,bs0,bs1,bs2,cs0,cs1,cs2
	if (exists("W_OptGradient")==1)
		Wave W_OptGradient=W_OptGradient
		printf "gradients = %s\r",vec2str(W_OptGradient,places=3)
	endif
	printf "Optimization took %.2f sec\r",seconds
	printf "		a* = {%+.4f,  %+.4f,  %+.4f}\r",as0,as1,as2
	printf "		b* = {%+.4f,  %+.4f,  %+.4f}\r",bs0,bs1,bs2
	printf "		c* = {%+.4f,  %+.4f,  %+.4f}\r",cs0,cs1,cs2
	KillWaves/Z W_OptGradient
	KillWaves/Z M_Lower,M_Upper,W_LUPermutation,M_x
	KillWaves/Z RefineLat_typXWave, RefineLat_xWave

	STRUCT crystalLattice xtal1					// final, or deformed lattice
	LatticeFromRecip(xtal1, as0,as1,as2,  bs0,bs1,bs2,  cs0,cs1,cs2)	// find lattice constants from recip
	printf "after optimizing, deformed lattice:\r"
	printf "		a = %.7f		b = %.7f		c = %.7f (nm)\r",xtal1.a,xtal1.b,xtal1.c
	printf "		alpha = %.5f	beta = %.5f	gam = %.5f¡\r",xtal1.alpha,xtal1.beta,xtal1.gam

	Variable maxStrain = StrainFromLat2Lat(xtal0,xtal1,!haveEnergy)			// only deviatoric if no energies
	Wave strainTensor=strainTensor
	Variable dilation = MatrixTrace(strainTensor)/3
	String dilationStr = " \tdilation = 0"
	if (abs(dilation)>1e-12)
		sprintf dilationStr, " \tdilation = %.1e",MatrixTrace(strainTensor)/3	// dilation
	endif
	printf "StrainTensor=\t%+.6f \t%+.6f \t%+.6f \tmax strain = %.1e\r",strainTensor[0][0],strainTensor[0][1],strainTensor[0][2],maxStrain
	printf "\t\t\t\t%+.6f \t%+.6f \t%+.6f%s\r",strainTensor[1][0],strainTensor[1][1],strainTensor[1][2],dilationStr
	printf "\t\t\t\t%+.6f \t%+.6f \t%+.6f\r",strainTensor[2][0],strainTensor[2][1],strainTensor[2][2]
	if (V_flag)
		printf "Optimize ended with error code = %d\r",V_flag
	endif
	return V_flag
End
//
//
Function totalErrorLatticeDev(w,as0,as1,as2, bs0,bs1,bs2, cs0,cs1)		// just deviatoric strain
	Wave w
	Variable as0,as1,as2						// componenets of a*
	Variable bs0,bs1,bs2						// componenets of b*
	Variable cs0,cs1							// componenets of c* (except for last one)

	Variable cs2, Vsc = NumberByKey("Vsc",note(w),"=")
	//		Vsc = ((a*)x(b*))¥(c*)
	cs2 = (Vsc - (as1*bs2-as2*bs1)*cs0 - (as2*bs0-as0*bs2)*cs1) / (as0*bs1-as1*bs0)
	return totalErrorLattice(w, as0,as1,as2, bs0,bs1,bs2, cs0,cs1,cs2)
End
//
//	Optimization function (we want to minimize this)
//	for each measured data point, the following information may be available
//		if a piece of info was not measured, then set it to NaN, and it will be ignored
//		e.g. if you know only one energy, then most of the w[i][6] will be NaN
//	w[i][0]			h	qhat is direction of a spot computed strictly from the recip x hkl
//	w[i][1]			k
//	w[i][2]			l
//	w[i][3]			measured x-pixel position (full chip & unbinned)
//	w[i][4]			measured y-pixel position (full chip & unbinned)
//	w[i][5]			weighting for use with x-y part
//	w[i][6]			measured energy (keV)
//	w[i][7]			weighting for this reflection (only used with energy)
//
Function totalErrorLattice(w, as0,as1,as2, bs0,bs1,bs2, cs0,cs1,cs2)
	Wave w										// parameters, including measured data
	Variable as0,as1,as2						// componenets of a*
	Variable bs0,bs1,bs2						// componenets of b*
	Variable cs0,cs1,cs2						// componenets of c*

	Wave qvec=root:Packages:geometry:totalError_qvec		// a Q vector, not normalized
	STRUCT microGeometry geo
	FillGeometryStructDefault(geo)
	Variable theta2pixel = (2*geo.dd/(geo.dpsy/geo.NyCCD))	// scales theta to pixels

	Variable d									// d-spacing (nm)
	Variable i, N=DimSize(w,0)				// number of hkl points measured
	Variable/C pz
	Variable err=1e-322						// really small, but not exactly zero
	Variable weight, thetaSpot, dSpot
	Variable h,k,l
	for (i=0;i<N;i+=1)
		if (hklStr2hkl(GetDimLabel(w,0,i),h,k,l)) // get the h,k,l for this spot
			continue							// skip this, no hkl entered for this point
		endif
		qvec[0] = h*as0 + k*bs0 + l*cs0		// calculate q-vector
		qvec[1] = h*as1 + k*bs1 + l*cs1
		qvec[2] = h*as2 + k*bs2 + l*cs2
		d = 2*PI/normalize(qvec)				// get d-spacing from the |qvec|
		if (numtype(d))
			continue							// cannot compute d, skip ths point
		endif

		if (!numtype(w[i][3]+w[i][4]))		// have valid pixel for this point
			pz = q2pixel(geo,qvec)
			weight = w[i][5]
			weight = numtype(weight) ? 1 : weight
			err += (real(pz)-w[i][3])^2 * weight	// the pixel errors
			err += (imag(pz)-w[i][4])^2 * weight
		endif
		if (!numtype(w[i][6]))				// have energy for this point
			thetaSpot = pixel2q(geo,w[i][3],w[i][4],qvec)	// Bragg angle (rad) of actual measured spot
			dspot = hc/w[i][6]/2/sin(thetaSpot)	// d-spacing from measured E and measured spot
			weight = w[i][7]
			weight = (numtype(weight) ? 1 : weight)*(theta2pixel/d)// makes err in d scale like pixels
			err += ((dspot-d)*weight)^2		// error in d from energy , with its own weighting
		endif
	endfor
	err = (err < 1e-321) ? NaN : err
	return err
End

Static Function LatticeFromRecip(xtal, as0,as1,as2,  bs0,bs1,bs2,  cs0,cs1,cs2)	// find lattice constants from recip
	STRUCT crystalLattice &xtal								// resultant xtal
	Variable as0,as1,as2,  bs0,bs1,bs2,  cs0,cs1,cs2		// the reicprocal lattice vectors (input)

	Variable Vsc											// volume of the reciprocal lattice, triple product
	Vsc = as0*(bs1*cs2-bs2*cs1) + as1*(bs2*cs0-bs0*cs2) + as2*(bs0*cs1-bs1*cs0)	// a¥(b x c)
	Variable pv = 2*PI/Vsc								// scale factor

	// compute direct lattice from reciprocal lattice
	Variable a0,a1,a2,  b0,b1,b2,  c0,c1,c2				//components of the direct lattice vectors
	a0=(bs1*cs2-bs2*cs1)*pv ;	 a1=(bs2*cs0-bs0*cs2)*pv ;	a2=(bs0*cs1-bs1*cs0)*pv	// (b x c)*2¹/Vc
	b0=(cs1*as2-cs2*as1)*pv ;	b1=(cs2*as0-cs0*as2)*pv ; 	b2=(cs0*as1-cs1*as0)*pv	// (c x a)*2¹/Vc
	c0=(as1*bs2-as2*bs1)*pv ;	c1=(as2*bs0-as0*bs2)*pv ;	c2=(as0*bs1-as1*bs0)*pv	// (a x b)*2¹/Vc

	xtal.a = sqrt(a0^2 + a1^2 + a2^2)
	xtal.b = sqrt(b0^2 + b1^2 + b2^2)
	xtal.c = sqrt(c0^2 + c1^2 + c2^2)
	xtal.alpha =	acos((b0*c0 + b1*c1 + b2*c2)/(xtal.b * xtal.c)) * 180/PI
	xtal.beta =	acos((a0*c0 + a1*c1 + a2*c2)/(xtal.a * xtal.c)) * 180/PI
	xtal.gam =	acos((a0*b0 + a1*b1 + a2*b2)/(xtal.a * xtal.b)) * 180/PI

	xtal.a0 = a0;	xtal.a1 = a1;	xtal.a2 = a2
	xtal.b0 = b0;	xtal.b1 = b1;	xtal.b2 = b2
	xtal.c0 = c0;	xtal.c1 = c1;	xtal.c2 = c2

	xtal.as0 = as0;	xtal.as1 = as1;	xtal.as2 = as2
	xtal.bs0 = bs0;	xtal.bs1 = bs1;	xtal.bs2 = bs2
	xtal.cs0 = cs0;	xtal.cs1 = cs1;	xtal.cs2 = cs2
//	LatticeSym#setDirectRecip(xtal)
	return 0
End
//Function testLatticeFromRecip()
//	STRUCT crystalLattice xtal0								// start xtal
//	FillLatticeStructDefault(xtal0)
//	xtal0.a = 0.5
//	xtal0.b = 0.6
//	xtal0.c = 0.7
//	xtal0.alpha = 95
//	xtal0.beta = 85
//	xtal0.gam = 110
//	LatticeSym#setDirectRecip(xtal0)
//
//	STRUCT crystalLattice xtal1								// resultant xtal
//	LatticeFromRecip(xtal1, xtal0.as0,xtal0.as1,xtal0.as2,  xtal0.bs0,xtal0.bs1,xtal0.bs2,  xtal0.cs0,xtal0.cs1,xtal0.cs2)
//	print_crystalLattice(xtal0)
//	print "\r=======================================================\r\r"
//	print_crystalLattice(xtal1)
//	print/d xtal0.Vc
//	print/d xtal1.Vc,"   ",(xtal0.Vc - xtal1.Vc) / xtal0.Vc
//	print  "delta lengths = ",(xtal1.a-xtal0.a),"   ",(xtal1.b-xtal0.b),"   ",(xtal1.c-xtal0.c)
//	print "delta angles = ",(xtal1.alpha-xtal0.alpha),"   ",(xtal1.beta-xtal0.beta),"   ",(xtal1.gam-xtal0.gam)
//End






Function doGeoRefinement(RefineData,hstr)
	Wave RefineData
	String hstr				// the hold string, a 1 means hold that parameter constant {alphad, betad, xbet, xgam, dd, xc, yc}
	Variable Np=strlen(hstr)						// total number of parameters that can be fit
	if (Np!=9)
		DoAlert 0, "if you do not know how to set hstr, use the Refinement Panel"
		return 1
	endif
	Variable i, Nopt									// number of parameters to optimize
	for (i=0,Nopt=0; i<Np; i+=1)
		Nopt += (str2num(hstr[i])==0)
	endfor
	String fitList="", allParams = "xalfd;xbetd;xbet;xgam;dd;xc;yc;dpsx;dpsy"
	for (i=0;i<Np;i+=1)
		fitList += SelectString(str2num(hstr[i]),StringFromList(i,allParams)+",","")
	endfor
	fitList = fitList[0,strsearch(fitList,",",Inf,1)-1]// remove trailing ","

	if (WhichListItem("dd",fitList,",")>=0 && WhichListItem("dpsx",fitList,",")>=0 && WhichListItem("dpsy",fitList,",")>=0)
		DoAlert 0, "Invalid combination, cannot fit height (dd) and the ccd size (dpsx & dpsy)"
		return 1
	endif

	if (!WaveExists(RefineData))
		if (stringmatch(WaveListClass("RefinementData","*",""),"RefineData;"))
			Wave RefineData = :RefineData		// RefineData exists, and it is the only one
		else
			String RefineDataStr = SelectString(WaveExists(RefineData),"",NameOfWave(RefineData))
			Prompt RefineDataStr,"Indexed Peak List",popup,WaveListClass("RefinementData","*","")
			DoPrompt/Help="3D-Xray Diffraction[Indexing]" "Refinement data",RefineDataStr
			if (V_flag)
				return 1
			endif
			Wave RefineData=$RefineDataStr
		endif
		if (!WaveExists(RefineData))
			return 1
		endif
		printf "¥doGeoRefinement(%s,\"%s\")\r",NameOfWave(RefineData),hstr
		printf "starting refinement using data from '%s'\r",NameOfWave(RefineData)
	endif

	String funcName = "totalErrorGeo"+num2istr(Nopt)// form function name from number of fitted parameters
	if (exists(funcName)!=6)
		DoAlert 0,"bad choice of optimizing parameters, no such function '"+funcName+"'"
		return 1
	endif
	if (strlen(FindTableWithWave(RefineData)))
		DoAlert 0, "You must close the table with Refinement Data before you can refine"
		return 1
	endif

	String wnote = note(RefineData), valList
	Variable xalfd=NumberByKey("xalfd",wnote,"="), xbetd=NumberByKey("xbetd",wnote,"=")
	Variable xbet=NumberByKey("xbet",wnote,"="), xgam=NumberByKey("xgam",wnote,"=")
	Variable dd=NumberByKey("dd",wnote,"=")
	Variable xc=NumberByKey("xcent",wnote,"="), yc=NumberByKey("ycent",wnote,"=")
	Variable dpsx=NumberByKey("dpsx",wnote,"="), dpsy=NumberByKey("dpsy",wnote,"=")
	sprintf valList,"xalfd=%g;xbetd=%g;xbet=%g;xgam=%g;dd=%g;xc=%g;yc=%g;dpsx=%g;dpsy=%g",xalfd,xbetd, xbet,xgam, dd, xc,yc,dpsx,dpsy

	Make/N=3/O/D root:Packages:geometry:totalError_ki,root:Packages:geometry:totalError_qhat
	Make/N=(Nopt)/O/D RefineGeo_typXWave, RefineGeo_xWave
	Wave typXWave=RefineGeo_typXWave, xWave=RefineGeo_xWave
	for(i=0;i<Nopt;i+=1)
		xWave[i] = NumberByKey(StringFromList(i,fitList,","),valList,"=")// fill starting array with values
	endfor
	typXWave = xWave
	wnote = ReplaceStringByKey("fitParameters",wnote,fitList,"=")
	Note/K RefineData, wnote
	printf "Optimizing '%s'\r",ReplaceString(",",fitList,", ")
	Variable ystart = totalErrorGeo(RefineData,xalfd,xbetd, xbet,xgam, dd, xc,yc,dpsx,dpsy)
	printf "starting point is xWave = %s,  with an error of %g\r", vec2str(typXWave,places=7), ystart
	Variable timer = startMSTimer
	Optimize/M={2,1}/R=typXWave/X=xWave/Y=(ystart) $funcName,RefineData
	printf "Optimization took %.2f sec\r",stopMSTimer(timer)/1e6
	if (V_flag)
		printf "Optimization failed with V_flag =",V_flag
		DoAlert 0, "Optimization failed, V_flag = "+num2str(V_flag)
		KillWaves/Z W_OptGradient
		KillWaves/Z M_Lower,M_Upper,W_LUPermutation,M_x
		KillWaves/Z RefineGeo_typXWave, RefineGeo_xWave
		return V_flag
	endif
	for(i=0;i<Nopt;i+=1)							// replace values in valList with optimized values
		valList = ReplaceNumberByKey(StringFromList(i,fitList,","),valList,xWave[i],"=")
	endfor
	STRUCT microGeometry geo
	Variable foundGeo = (!FillGeometryStructDefault(geo))
	geo.xalfd=NumberByKey("xalfd",valList,"=");	geo.xbetd=NumberByKey("xbetd",valList,"=")
	geo.xbet=NumberByKey("xbet",valList,"=");	geo.xgam=NumberByKey("xgam",valList,"=")
	geo.dd=NumberByKey("dd",valList,"=")
	geo.xcent=NumberByKey("xc",valList,"=");		geo.ycent=NumberByKey("yc",valList,"=")
	geo.dpsx=NumberByKey("dpsx",valList,"=");	geo.dpsy=NumberByKey("dpsy",valList,"=")

	printf "after optimizing, \r"
	printf "		xalfd = %.5f	xbetd = %.5f\r",geo.xalfd,geo.xbetd
	printf "		xbeta = %.5f	xgamma = %.5f\r",geo.xbet,geo.xgam
	printf "		xc = %.3f		yc = %.3f\r",geo.xcent,geo.ycent
	printf "		dd = %.4f\r",geo.dd
	printf "		dpsx = %.4f 	dpsy = %.4f\r",geo.dpsx,geo.dpsy
	if (WhichListItem("dpsx",fitList,",")>=0 || WhichListItem("dpsy",fitList,",")>=0)
		dd = geo.dd
		dpsx = geo.dpsx
		dpsy = geo.dpsy
		print "  or"
		Variable scale = round((dpsx+dpsy)/2)/((dpsx+dpsy)/2)
		printf "		dd = %.4f\r",dd*scale
		printf "		dpsx = %.4f 	dpsy = %.4f\r",dpsx*scale,dpsy*scale
		print "  or"
		scale = round(50)/dpsx
		printf "		dd = %.4f\r",dd*scale
		printf "		dpsx = %.4f 	dpsy = %.4f\r",dpsx*scale,dpsy*scale
	endif

	if (foundGeo)
		String strStruct								// need a string version to pass to MakeGeometryParametersPanel()
		StructPut/S/B=2 geo, strStruct					// put back into string so I can pass it
		Variable/G root:Packages:geometry:PanelValues:dirty=0
		DoWindow/K GeometrySet						// ensure a new windwo
		MakeGeometryParametersPanel(strStruct)
	endif
	KillWaves/Z W_OptGradient
	KillWaves/Z M_Lower,M_Upper,W_LUPermutation,M_x
	KillWaves/Z RefineGeo_typXWave, RefineGeo_xWave
	return 0
End
//
Function totalErrorGeo1(w,x1)
	Wave w
	Variable x1
	String wnote = note(w)
	Variable xalfd=NumberByKey("xalfd",wnote,"="), xbetd=NumberByKey("xbetd",wnote,"=")
	Variable xbet=NumberByKey("xbet",wnote,"="), xgam=NumberByKey("xgam",wnote,"=")
	Variable dd=NumberByKey("dd",wnote,"=")
	Variable xc=NumberByKey("xcent",wnote,"="), yc=NumberByKey("ycent",wnote,"=")
	Variable dpsx=NumberByKey("dpsx",wnote,"="), dpsy=NumberByKey("dpsy",wnote,"=")
	String item,list = StringByKey("fitParameters",wnote,"=")
	item = StringFromList(0,list,",")
	setParameterValue(item,x1,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	return totalErrorGeo(w,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
End
//
Function totalErrorGeo2(w,x1,x2)
	Wave w
	Variable x1,x2
	String wnote = note(w)
	Variable xalfd=NumberByKey("xalfd",wnote,"="), xbetd=NumberByKey("xbetd",wnote,"=")
	Variable xbet=NumberByKey("xbet",wnote,"="), xgam=NumberByKey("xgam",wnote,"=")
	Variable dd=NumberByKey("dd",wnote,"=")
	Variable xc=NumberByKey("xcent",wnote,"="), yc=NumberByKey("ycent",wnote,"=")
	Variable dpsx=NumberByKey("dpsx",wnote,"="), dpsy=NumberByKey("dpsy",wnote,"=")
	String list = StringByKey("fitParameters",wnote,"=")
	setParameterValue(StringFromList(0,list,","),x1,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(1,list,","),x2,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	return totalErrorGeo(w,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
End
//
Function totalErrorGeo3(w,x1,x2,x3)
	Wave w
	Variable x1,x2,x3
	String wnote = note(w)
	Variable xalfd=NumberByKey("xalfd",wnote,"="), xbetd=NumberByKey("xbetd",wnote,"=")
	Variable xbet=NumberByKey("xbet",wnote,"="), xgam=NumberByKey("xgam",wnote,"=")
	Variable dd=NumberByKey("dd",wnote,"=")
	Variable xc=NumberByKey("xcent",wnote,"="), yc=NumberByKey("ycent",wnote,"=")
	Variable dpsx=NumberByKey("dpsx",wnote,"="), dpsy=NumberByKey("dpsy",wnote,"=")
	String list = StringByKey("fitParameters",wnote,"=" )
	setParameterValue(StringFromList(0,list,","),x1,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(1,list,","),x2,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(2,list,","),x3,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	return totalErrorGeo(w,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
End
//
Function totalErrorGeo4(w,x1,x2,x3,x4)
	Wave w
	Variable x1,x2,x3,x4
	String wnote = note(w)
	Variable xalfd=NumberByKey("xalfd",wnote,"="), xbetd=NumberByKey("xbetd",wnote,"=")
	Variable xbet=NumberByKey("xbet",wnote,"="), xgam=NumberByKey("xgam",wnote,"=")
	Variable dd=NumberByKey("dd",wnote,"=")
	Variable xc=NumberByKey("xcent",wnote,"="), yc=NumberByKey("ycent",wnote,"=")
	Variable dpsx=NumberByKey("dpsx",wnote,"="), dpsy=NumberByKey("dpsy",wnote,"=")
	String list = StringByKey("fitParameters",wnote,"=" )
	setParameterValue(StringFromList(0,list,","),x1,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(1,list,","),x2,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(2,list,","),x3,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(3,list,","),x4,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	return totalErrorGeo(w,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
End
//
Function totalErrorGeo5(w,x1,x2,x3,x4,x5)
	Wave w
	Variable x1,x2,x3,x4,x5
	String wnote = note(w)
	Variable xalfd=NumberByKey("xalfd",wnote,"="), xbetd=NumberByKey("xbetd",wnote,"=")
	Variable xbet=NumberByKey("xbet",wnote,"="), xgam=NumberByKey("xgam",wnote,"=")
	Variable dd=NumberByKey("dd",wnote,"=")
	Variable xc=NumberByKey("xcent",wnote,"="), yc=NumberByKey("ycent",wnote,"=")
	Variable dpsx=NumberByKey("dpsx",wnote,"="), dpsy=NumberByKey("dpsy",wnote,"=")
	String list = StringByKey("fitParameters",wnote,"=" )
	setParameterValue(StringFromList(0,list,","),x1,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(1,list,","),x2,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(2,list,","),x3,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(3,list,","),x4,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(4,list,","),x5,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	return totalErrorGeo(w,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
End
//
Function totalErrorGeo6(w,x1,x2,x3,x4,x5,x6)
	Wave w
	Variable x1,x2,x3,x4,x5,x6
	String wnote = note(w)
	Variable xalfd=NumberByKey("xalfd",wnote,"="), xbetd=NumberByKey("xbetd",wnote,"=")
	Variable xbet=NumberByKey("xbet",wnote,"="), xgam=NumberByKey("xgam",wnote,"=")
	Variable dd=NumberByKey("dd",wnote,"=")
	Variable xc=NumberByKey("xcent",wnote,"="), yc=NumberByKey("ycent",wnote,"=")
	Variable dpsx=NumberByKey("dpsx",wnote,"="), dpsy=NumberByKey("dpsy",wnote,"=")
	String list = StringByKey("fitParameters",wnote,"=" )
	setParameterValue(StringFromList(0,list,","),x1,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(1,list,","),x2,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(2,list,","),x3,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(3,list,","),x4,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(4,list,","),x5,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(5,list,","),x6,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	return totalErrorGeo(w,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
End
//
Function totalErrorGeo7(w,x1,x2,x3,x4,x5,x6,x7)
	Wave w
	Variable x1,x2,x3,x4,x5,x6,x7
	String wnote = note(w)
	Variable xalfd=NumberByKey("xalfd",wnote,"="), xbetd=NumberByKey("xbetd",wnote,"=")
	Variable xbet=NumberByKey("xbet",wnote,"="), xgam=NumberByKey("xgam",wnote,"=")
	Variable dd=NumberByKey("dd",wnote,"=")
	Variable xc=NumberByKey("xcent",wnote,"="), yc=NumberByKey("ycent",wnote,"=")
	Variable dpsx=NumberByKey("dpsx",wnote,"="), dpsy=NumberByKey("dpsy",wnote,"=")
	String list = StringByKey("fitParameters",wnote,"=" )
	setParameterValue(StringFromList(0,list,","),x1,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(1,list,","),x2,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(2,list,","),x3,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(3,list,","),x4,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(4,list,","),x5,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(5,list,","),x6,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(6,list,","),x7,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	return totalErrorGeo(w,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
End
//
Function totalErrorGeo8(w,x1,x2,x3,x4,x5,x6,x7,x8)
	Wave w
	Variable x1,x2,x3,x4,x5,x6,x7,x8
	String wnote = note(w)
	Variable xalfd=NumberByKey("xalfd",wnote,"="), xbetd=NumberByKey("xbetd",wnote,"=")
	Variable xbet=NumberByKey("xbet",wnote,"="), xgam=NumberByKey("xgam",wnote,"=")
	Variable dd=NumberByKey("dd",wnote,"=")
	Variable xc=NumberByKey("xcent",wnote,"="), yc=NumberByKey("ycent",wnote,"=")
	Variable dpsx=NumberByKey("dpsx",wnote,"="), dpsy=NumberByKey("dpsy",wnote,"=")
	String list = StringByKey("fitParameters",wnote,"=" )
	setParameterValue(StringFromList(0,list,","),x1,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(1,list,","),x2,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(2,list,","),x3,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(3,list,","),x4,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(4,list,","),x5,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(5,list,","),x6,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(6,list,","),x7,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(7,list,","),x8,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	return totalErrorGeo(w,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
End
//
Function totalErrorGeo9(w,x1,x2,x3,x4,x5,x6,x7,x8,x9)
	Wave w
	Variable x1,x2,x3,x4,x5,x6,x7,x8,x9
	String wnote = note(w)
	Variable xalfd=NumberByKey("xalfd",wnote,"="), xbetd=NumberByKey("xbetd",wnote,"=")
	Variable xbet=NumberByKey("xbet",wnote,"="), xgam=NumberByKey("xgam",wnote,"=")
	Variable dd=NumberByKey("dd",wnote,"=")
	Variable xc=NumberByKey("xcent",wnote,"="), yc=NumberByKey("ycent",wnote,"=")
	Variable dpsx=NumberByKey("dpsx",wnote,"="), dpsy=NumberByKey("dpsy",wnote,"=")
	String list = StringByKey("fitParameters",wnote,"=" )
	setParameterValue(StringFromList(0,list,","),x1,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(1,list,","),x2,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(2,list,","),x3,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(3,list,","),x4,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(4,list,","),x5,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(5,list,","),x6,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(6,list,","),x7,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(7,list,","),x8,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	setParameterValue(StringFromList(8,list,","),x9,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	return totalErrorGeo(w,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
End
//
Static Function setParameterValue(item,xx,xalfd,xbetd,xbet,xgam,dd,xc,yc,dpsx,dpsy)
	String item				// name of a variable to pass
	Variable xx
	Variable &xalfd,&xbetd,&xbet,&xgam,&dd,&xc,&yc,&dpsx,&dpsy	// variable values
	strswitch(item)
		case "xalfd":
			xalfd = xx
			break
		case "xbetd":
			xbetd = xx
			break
		case "xbet":
			xbet = xx
			break
		case "xgam":
			xgam = xx
			break
		case "dd":
			dd = xx
			break
		case "xc":
			xc = xx
			break
		case "yc":
			yc = xx
			break
		case "dpsx":
			dpsx = xx
			break
		case "dpsy":
			dpsy = xx
			break
		default:
		Abort "unknown fit parameter '"+item+"'"
	endswitch
	return 0
End
//
//Static Function/S Wave2str(w,m)
//	Wave w						// wave
//	Variable m					// number of places
//	if (!WaveExists(w))
//		return "{}"
//	endif
//	String num,fmt,str="{"
//	sprintf fmt,"%%.%dg",m
//	Variable i,N=min(numpnts(w),15)	// at most 15 terms
//	for (i=0;i<N-1;i+=1)
//		sprintf num,fmt,w[i]
//		str += num+","
//	endfor
//	sprintf num,fmt,w[i]
//	str += num
//	str += SelectString(numpnts(w)>15,"",",...")
//	return str+"}"
//End



//	Optimization function (we want to minimize this)
//	for each measured data point, the following information may be available
//		if a piece of info was not measured, then set it to NaN, and it will be ignored
//		e.g. if you know only one energy, then most of the w[i][6] will be NaN
//	w[i][0]			qhat, x	qhat is direction of a spot computed strictly from the recip x hkl
//	w[i][1]			qhat, y		NOTE, qhat is assumed to be normalized
//	w[i][2]			qhat, z
//	w[i][3]			measured x-pixel position (full chip & unbinned)
//	w[i][4]			measured y-pixel position (full chip & unbinned)
//	w[i][5]			weighting for use with x-y part
//	w[i][6]			measured energy (keV)
//	w[i][7]			weighting for this reflection (only used with energy)
//	w[i][8]			theta pre-computed from measured energy and known d-spacing (radian)
//
Function totalErrorGeo(w,xalfd, xbetd,xbet, xgam,dd,xc, yc, dpsx, dpsy)
	Wave w										// parameters, including measured data
	Variable xalfd, xbetd						// variables to fit
	Variable xbet, xgam
	Variable dd
	Variable xc, yc
	Variable dpsx, dpsy

	Wave qhat=root:Packages:geometry:totalError_qhat	// holds q-vector for a reflection
	Wave ki=root:Packages:geometry:totalError_ki		// incident beam direction

	String strStruct=StrVarOrDefault(":geoStructStr","")	// set to values in current directory
	if (strlen(strStruct)<1)
		strStruct=StrVarOrDefault("root:Packages:geometry:geoStructStr","")	// try the default values
	endif
	STRUCT microGeometry geo
	StructGet/S/B=2 geo, strStruct				// load geometry
	geo.xalfd = xalfd								// set local structure to contain the input values being fitted
	geo.xbetd = xbetd
	geo.xbet = xbet
	geo.xgam = xgam
	geo.dd = dd
	geo.xcent = xc
	geo.ycent = yc
	geo.dpsx = dpsx
	geo.dpsy = dpsy

	GeometryUpdateCalc(geo)
	Variable validXtal=0
	STRUCT crystalLattice xtal
	strStruct = StrVarOrDefault(":xtalStructStr","")
	if (strlen(strStruct))
		StructGet/S/B=2 xtal, strStruct
		validXtal = 1
	endif
	Variable theta2pixel = (2*geo.dd/(geo.dpsy/geo.NyCCD))	// scales theta to pixels

	Variable i, N=DimSize(w,0)				// number of hkl points measured
	Variable thetaIndex							// Bragg angle of a point from indexation (radian)
	Variable thetaEne							// Bragg angle from energy & known d (radian)
	Variable/C pz
	Variable err=1e-322						// really small, but not exactly zero
	Variable weight
	for (i=0;i<N;i+=1)
		if (!numtype(w[i][0]+w[i][1]+w[i][2]+w[i][3]+w[i][4]))	// have Q-data for this point
			qhat = w[i][p]
			pz = q2pixel(geo,qhat)
			weight = w[i][5]
			weight = numtype(weight) ? 1 : weight
			err += (real(pz)-w[i][3])^2 * weight	// the pixel errors
			err += (imag(pz)-w[i][4])^2 * weight
		endif
		if (!numtype(w[i][8]))				// have energy (actually theta from the energy) for this point
			ki = geo.ki[p]
			thetaIndex = asin(-MatrixDot(ki,qhat))	// computed Bragg angle (radian)
			weight = w[i][7]
			weight = (numtype(weight) ? 1 : weight)*theta2pixel	// makes err in theta scale like pixels
			err += ((thetaIndex-w[i][8])*weight)^2 	// error in theta from energy , with its own weighting
		endif
	endfor
	err = (err < 1e-321) ? NaN : err
	return err
End






// Make (if it does not exist) and view the data needed for refining geometry parameters
Function FillGeoRefineDataFromIndexing(FullPeakList,FullPeakIndexed)
	Wave FullPeakList,FullPeakIndexed
	if (!WaveExists(FullPeakList) || !WaveExists(FullPeakIndexed))
		String peakListStr = SelectString(WaveExists(FullPeakList),"",NameOfWave(FullPeakList))
		String indexPeakList = SelectString(WaveExists(FullPeakIndexed),"",NameOfWave(FullPeakIndexed))
		Prompt peakListStr,"wave with fitted peaks",popup,WaveListClass("FittedPeakList","*","DIMS:2,MAXCOLS:11,MINCOLS:11")
		Prompt indexPeakList,"Indexed Peak List",popup,WaveListClass("IndexedPeakList","*","")
		DoPrompt/Help="3D-Xray Diffraction[Indexing]" "list of peaks and indexing",peakListStr,indexPeakList
		if (V_flag)
			return 1
		endif
		Wave FullPeakList=$peakListStr
		Wave FullPeakIndexed = $indexPeakList
		if (!WaveExists(FullPeakList) || !WaveExists(FullPeakIndexed))
			return 1
		endif
		printf "for the list of fitted pixels = '%s',   and for the indexed peaks using '%s'\r",NameOfWave(FullPeakList),NameOfWave(FullPeakIndexed)
	endif

	String wnote = note(FullPeakList)
	Variable startx,groupx, starty,groupy			// ROI of the original image
	startx = NumberByKey("startx",wnote,"=")
	groupx = NumberByKey("groupx",wnote,"=")
	starty = NumberByKey("starty",wnote,"=")
	groupy = NumberByKey("groupy",wnote,"=")
	startx = numtype(startx) ? 1 : startx
	groupx = numtype(groupx) ? 1 : groupx
	starty = numtype(starty) ? 1 : starty
	groupy = numtype(groupy) ? 1 : groupy

	STRUCT microGeometry geo
	String strStruct=StrVarOrDefault(":geoStructStr","")	// set to values in current directory
	if (strlen(strStruct)<1)
		strStruct=StrVarOrDefault("root:Packages:geometry:geoStructStr","")	// try the default values
	endif
	if (strlen(strStruct)<1)								// geo structure not found
		DoAlert 0, "no geometry structure found, did you forget to set it?"
		return 1
	endif
	StructGet/S/B=2 geo, strStruct							// load geometry
	wnote = ""
	wnote = ReplaceNumberByKey("dpsx",wnote,geo.dpsx,"=")
	wnote = ReplaceNumberByKey("dpsy",wnote,geo.dpsy,"=")
	wnote = ReplaceNumberByKey("NxCCD",wnote,geo.NxCCD,"=")
	wnote = ReplaceNumberByKey("NyCCD",wnote,geo.NyCCD,"=")
	wnote = ReplaceNumberByKey("xcent",wnote,geo.xcent,"=")
	wnote = ReplaceNumberByKey("ycent",wnote,geo.ycent,"=")
	wnote = ReplaceNumberByKey("ddOffset",wnote,geo.ddOffset,"=")
	wnote = ReplaceNumberByKey("dd",wnote,geo.dd,"=")
	wnote = ReplaceNumberByKey("xbet",wnote,geo.xbet,"=")
	wnote = ReplaceNumberByKey("xgam",wnote,geo.xgam,"=")
	wnote = ReplaceNumberByKey("xalfd",wnote,geo.xalfd,"=")
	wnote = ReplaceNumberByKey("xbetd",wnote,geo.xbetd,"=")
	wnote = ReplaceNumberByKey("wireF",wnote,geo.wire.F,"=")
	wnote = ReplaceNumberByKey("wireH0",wnote,geo.wire.H0,"=")
	wnote = ReplaceNumberByKey("wireHyc",wnote,geo.wire.Hyc,"=")
	wnote = ReplaceNumberByKey("wireX",wnote,geo.wire.X,"=")
	wnote = ReplaceNumberByKey("wiredia",wnote,geo.wire.dia,"=")
	strStruct = StrVarOrDefault(":xtalStructStr","")
	if (strlen(strStruct))
		STRUCT crystalLattice xtal
		StructGet/S/B=2 xtal, strStruct					// load xtal structure
		wnote = ReplaceNumberByKey("xtal.a",wnote,xtal.a,"=")
		wnote = ReplaceNumberByKey("xtal.b",wnote,xtal.b,"=")
		wnote = ReplaceNumberByKey("xtal.c",wnote,xtal.c,"=")
		wnote = ReplaceNumberByKey("xtal.alpha",wnote,xtal.alpha,"=")
		wnote = ReplaceNumberByKey("xtal.beta",wnote,xtal.beta,"=")
		wnote = ReplaceNumberByKey("xtal.gam",wnote,xtal.gam,"=")
		wnote = ReplaceNumberByKey("xtal.structNum",wnote,xtal.structNum,"=")
	endif
	wnote = ReplaceStringByKey("recip_lattice",wnote,StringByKey("recip_lattice0",note(FullPeakIndexed),"="),"=")
	wnote = ReplaceStringByKey("FullPeakList",wnote,GetWavesDataFolder(FullPeakList,2),"=" )
	wnote = ReplaceStringByKey("FullPeakIndexed",wnote,GetWavesDataFolder(FullPeakIndexed,2),"=" )
	wnote = ReplaceStringByKey("waveClass",wnote,"RefinementData","=" )

	Variable i,N = DimSize(FullPeakIndexed,0)	// number of measured points for the refinement	
	Variable Nmeas = DimSize(FullPeakList,0)	// number of measured points
	Make/N=(N,11)/O/D RefineData
	RefineData = NaN
	Note/K RefineData,wnote
	SetDimLabel 1,0,'q^ x',RefineData			//	RefineData[i][0]		qhat, x	qhat is direction of spot calculated only from recip x hkl
	SetDimLabel 1,1,'q^ y',RefineData			//	RefineData[i][1]		qhat, y		NOTE, qhat is assumed to be normalized
	SetDimLabel 1,2,'q^ z',RefineData			//	RefineData[i][2]		qhat, z
	SetDimLabel 1,3,'x pixel',RefineData		//	RefineData[i][3]		measured x-pixel position (full chip & unbinned)
	SetDimLabel 1,4,'y pixel',RefineData		//	RefineData[i][4]		measured y-pixel position (full chip & unbinned)
	SetDimLabel 1,5,'pixel weight',RefineData	//	RefineData[i][5]		weighting for use with x-y part
	SetDimLabel 1,6,'Energy (keV)',RefineData	//	RefineData[i][6]		measured energy (keV)
	SetDimLabel 1,7,'Energy weight',RefineData	//	RefineData[i][7]		weighting for this reflection (only used with energy)
	SetDimLabel 1,8,'theta (rad)',RefineData	//	RefineData[i][8]		theta auto-computed from RefineData[i][6] and known d-spacing
	SetDimLabel 1,9,'X linearized pixel',RefineData	//	RefineData[i][9]	same as RefineData[i][3], but after distortion correction
	SetDimLabel 1,10,'Y linearized pixel',RefineData//	RefineData[i][10]
	SetDimLabel 0,-1,'(hkl)',RefineData
	String str
	Variable dist, distMin, j,jmin, px,py
	for (i=0;i<N;i+=1)
		RefineData[i][0,2] = FullPeakIndexed[i][q]	// qhat
		px = FullPeakIndexed[i][9]					// computed peak position from orientation matrix
		py = FullPeakIndexed[i][10]
		// find measured peak that is closest
		for (j=0,distMin=Inf;j<Nmeas;j+=1)
			dist = (px-FullPeakList[j][0])^2 + (py-FullPeakList[j][1])^2
			if (dist < distMin)
				distMin = dist
				jmin = j
			endif
		endfor
//		RefineData[i][3,4] = FullPeakList[jmin][q-3]
		px = FullPeakList[jmin][0]
		py = FullPeakList[jmin][1]
		px = (startx-1) + px*groupx + (groupy-1)/2	// change to unbinned full chip pixels
		py = (starty-1) + py*groupy + (groupy-1)/2
		RefineData[i][3] = px
		RefineData[i][4] = py
		sprintf str,"(%d, %d, %d)",FullPeakIndexed[i][3],FullPeakIndexed[i][4],FullPeakIndexed[i][5]
		SetDimLabel 0,i,$str,RefineData
		RefineData[i][5] = 1
	endfor
	if (exists("setE")==6)
		Execute "setE()"
	endif
	Make/N=3/O/D root:Packages:geometry:totalError_ki,root:Packages:geometry:totalError_qhat
	Make/N=3/O/D root:Packages:geometry:totalError_qvec
	return 0
End
//
//
//
// Make (if it does not exist) and view the data needed for refining geometry parameters
Function EditGeoRefinementData(N)
	Variable N					// number of measured points for the refinement
	if (!(N>0) && exists("RefineData")!=1)	// if wave exits, no need to make one
		N = 10
		Prompt N,  "number of data points for refinement"
		DoPrompt "number of points",N
		if (V_flag)
			return 1
		endif
		if (!(N>0))				// invalid size, must be > 0
			return 1
		endif
	endif
	if (exists("RefineData")!=1)
		Make/N=(N,11)/O/D RefineData
		RefineData = NaN
		Note/K RefineData,"waveClass=RefinementData"
	else
		N = 0
	endif

	SetDimLabel 1,0,'q^ x',RefineData			//	RefineData[i][0]		qhat, x	qhat is direction of spot calculated only from recip x hkl
	SetDimLabel 1,1,'q^ y',RefineData			//	RefineData[i][1]		qhat, y		NOTE, qhat is assumed to be normalized
	SetDimLabel 1,2,'q^ z',RefineData			//	RefineData[i][2]		qhat, z
	SetDimLabel 1,3,'x pixel',RefineData		//	RefineData[i][3]		measured x-pixel position (full chip & unbinned)
	SetDimLabel 1,4,'y pixel',RefineData		//	RefineData[i][4]		measured y-pixel position (full chip & unbinned)
	SetDimLabel 1,5,'pixel weight',RefineData	//	RefineData[i][5]		weighting for use with x-y part
	SetDimLabel 1,6,'Energy (keV)',RefineData	//	RefineData[i][6]		measured energy (keV)
	SetDimLabel 1,7,'Energy weight',RefineData	//	RefineData[i][7]		weighting for this reflection (only used with energy)
	SetDimLabel 1,8,'theta (rad)',RefineData	//	RefineData[i][8]		theta auto-computed from RefineData[i][6] and known d-spacing
	SetDimLabel 1,9,'X linearized pixel',RefineData	//	RefineData[i][9]	same as RefineData[i][3], but after distortion correction
	SetDimLabel 1,10,'Y linearized pixel',RefineData//	RefineData[i][10]
	Variable i
	SetDimLabel 0,-1,'(hkl)',RefineData
	for (i=0;i<N;i+=1)
		SetDimLabel 0,i,'h k l',RefineData
	endfor
	String table = FindTableWithWave(RefineData)
	if (strlen(table))
		DoWindow/F $table
	else
		Edit/K=1/W=(5,44,757,513) RefineData.ld
		ModifyTable width(Point)=30,width(RefineData.d)=78,width(RefineData.l)=80
		SetWindow kwTopWin,hook(onClose)=RefinementDataTableCloseHook 
		SetWindow kwTopWin userdata(onKillWave)=GetWavesDataFolder(RefineData,2)
	endif
	// this is needed for totalErrorGeo(), do it here just in case
	Make/N=3/O/D root:Packages:geometry:totalError_ki,root:Packages:geometry:totalError_qhat
	Make/N=3/O/D root:Packages:geometry:totalError_qvec
	return 0
End
//
// finds the name of an existing table containing the wave w0
// if the wave cannot be found in existing tables, it returns an empty string.
Static Function/S FindTableWithWave(w0)	// find the table window which contains the specified wave
	Wave w0
	if (!WaveExists(w0))
		return ""
	endif
	String name0=GetWavesDataFolder(w0,2)		// full path name of desired wave
	String table,tlist = WinList("*",";","WIN:2")
	Variable i, m,Nm=ItemsInList(tlist)
	for (m=0;m<Nm;m+=1)						// loop over all displayed tables
		table = StringFromList(m,tlist)				// table name
		i = 0
		do											// check each wave in this table
			Wave wav = WaveRefIndexed(table,i,1)	// 1 means table data
			if (!WaveExists(wav))					// no more waves in this table
				break
			endif
			if (stringmatch(GetWavesDataFolder(wav,2),name0))
				return table						// found our wave, return table name
			endif
			i += 1
		while(1)
	endfor
	return ""										// no exising displayed table found
End
//
Function RefinementDataTableCloseHook(s)
	STRUCT WMWinHookStruct &s
	if (s.eventCode!=2)								//  only process kill events
		return 0
	endif

	String wName = GetUserData("","","onKillWave")
	Wave wav = WaveRefIndexed(s.winName,0,1)	// 1 means table data
	if (!WaveExists(wav))
		return 0
	elseif (!stringmatch(GetWavesDataFolder(wav,2),wName))
		return 0
	elseif (!stringmatch(StringByKey("waveClass",note(wav),"="),"RefinementData"))
		return 0
	elseif(DimSize(wav,1)<11)					// not enough columns (probably an older version)
		DoAlert 0, "RefineData[][] does not have enough columns, remake it"
		return 0
	endif

	Variable validXtal=0
	STRUCT crystalLattice xtal
	String strStruct = StrVarOrDefault(":xtalStructStr","")
	if (strlen(strStruct))
		StructGet/S/B=2 xtal, strStruct
		validXtal = 1
		String wnote = note(wav)
		wnote = ReplaceNumberByKey("xtal.a",wnote,xtal.a,"=")
		wnote = ReplaceNumberByKey("xtal.b",wnote,xtal.b,"=")
		wnote = ReplaceNumberByKey("xtal.c",wnote,xtal.c,"=")
		wnote = ReplaceNumberByKey("xtal.alpha",wnote,xtal.alpha,"=")
		wnote = ReplaceNumberByKey("xtal.beta",wnote,xtal.beta,"=")
		wnote = ReplaceNumberByKey("xtal.gam",wnote,xtal.gam,"=")
		wnote = ReplaceNumberByKey("xtal.structNum",wnote,xtal.structNum,"=")
	endif

	Make/N=3/O/D RefinementDataTaCloseHook_qhat
	Wave qhat = RefinementDataTaCloseHook_qhat
	Variable i,N=DimSize(wav,0)
	Variable d, h,k,l
	wav[][8] = NaN
	Variable px,py
	Wave xymap=root:Packages:geometry:xymap
	for (i=0;i<N;i+=1)
		if (!numtype(wav[i][0]+wav[i][1]+wav[i][2]+wav[i][3]+wav[i][4]))		// make sure that all qhat are normalized
			qhat = wav[i][p]
			normalize(qhat)
			wav[i][0,2] = qhat[q]
			wav[i][5] = numtype(wav[i][5]) ? 1 : wav[i][5]	// ensure valid weighting
			px = wav[i][3]
			py = wav[i][4]
			peakcorrection2(xymap,px,py)						// makes px,py distortion corrected (used to accelerate optimization)
			wav[i][9] = px										// store distortion corrected pixels too, this way OptimizeAllTilts()
			wav[i][10] = py									// does not have to do distorion corrections.
		else
			wav[i][5] = NaN
		endif
		if (!numtype(wav[i][6]) && validXtal)					// for valid energy input & xtal data
			wav[i][7] = numtype(wav[i][7]) ? 1 : wav[i][7]	// ensure valid weighting
			if (hklStr2hkl(GetDimLabel(wav,0,i),h,k,l))
				continue
			endif
 			d = 2*PI/QvecFrom_HKL(xtal,h,k,l,$"")
			wav[i][8] = asin(hc/(2*d*wav[i][6]))			// alywas fill in computed theta
		else
			wav[i][7] = NaN
		endif
	endfor
	KillWaves/Z RefinementDataTaCloseHook_qhat
	return 1
End


Function QvecFrom_HKL(xtal,h,k,l,Qvec)	// returns length of qvector, and optionally Qvec too.
	STRUCT crystalLattice &xtal
	Variable h,k,l
	Wave Qvec					// optional wave to recieve the result of:  (recip) x (hkl)

	Variable qx,qy,qz
	qx = h*xtal.as0 + k*xtal.bs0 + l*xtal.cs0
	qy = h*xtal.as1 + k*xtal.bs1 + l*xtal.cs1
	qz = h*xtal.as2 + k*xtal.bs2 + l*xtal.cs2
	if (WaveExists(Qvec))
		Qvec[0] = qx
		Qvec[1] = qy
		Qvec[2] = qz
	endif
	return sqrt(qx*qx+qy*qy+qz*qz)
End

Function hklStr2hkl(hklStr,h,k,l)
	String hklStr
	Variable &h,&k,&l
	hklStr = ReplaceString(",",hklStr," ")
	hklStr = ReplaceString("(",hklStr," ")
	hklStr = ReplaceString(")",hklStr," ")
	sscanf hklStr,"%d %d %d",h,k,l
	if (V_flag!=3)
		h = NaN
		k = NaN
		l = NaN
		return 1
	endif
	return 0
End



Static Function StrainFromLat2Lat(xtal0,xtal1,deviatoric)	// find strain matrix that goes from xtal0 to xtal1
	STRUCT crystalLattice &xtal0
	STRUCT crystalLattice &xtal1
	Variable deviatoric									// if true, only report deviatoric strain
	//	A x (abc) = (a'b'c')
	//	A = (a'b'c') x (abc)^-1
	//	A = lat1 x lat0^-1
	Wave lat0 = $(microGeo#MakeUnique3x3Mat($""))
	Wave lat1 = $(microGeo#MakeUnique3x3Mat($""))
	lat0[0,2][0] = {xtal0.a0,xtal0.a1,xtal0.a2}
	lat0[0,2][1] = {xtal0.b0,xtal0.b1,xtal0.b2}
	lat0[0,2][2] = {xtal0.c0,xtal0.c1,xtal0.c2}
	lat1[0,2][0] = {xtal1.a0,xtal1.a1,xtal1.a2}
	lat1[0,2][1] = {xtal1.b0,xtal1.b1,xtal1.b2}
	lat1[0,2][2] = {xtal1.c0,xtal1.c1,xtal1.c2}
	MatrixOp/O strainTensor = lat1 x Inv(lat0)
	MatrixOp/O strainTensor = (strainTensor + strainTensor^t)/2 - Identity(3)
	if (deviatoric)										// make the trace zero for deviatoric strain
		Variable Delta = MatrixTrace(strainTensor)/3
		MatrixOp/O strainTensor = strainTensor - Delta*Identity(3)
	endif
	WaveStats/Q/M=1 strainTensor
	KillWaves/Z lat0,lat1
	return max(abs(V_max),abs(V_min))
End




Function MakeRefinePanel()
	if (strlen(WinList("RefinementPanel","","WIN:64")))
		DoWindow/F RefinementPanel
		return 0
	endif
	NewPanel/K=1/W=(728,259,907,425) as "Refinement"
	DoWindow/C RefinementPanel
	SetDrawLayer UserBack
	SetDrawEnv fname= "Lucida Grande"
	DrawText 11,22,"Refine Parameters"
	CheckBox alphad,pos={25,30},size={56,15},title="alphad"
		CheckBox alphad,help={"detector tilt (roll)"},fSize=12,value= 1
	CheckBox betad,pos={100,30},size={52,15},title="betad"
		CheckBox betad,help={"detector tilt (pitch, forward/back)"},fSize=12,value= 1
	CheckBox xbet,pos={25,50},size={45,15},title="xbet"
		CheckBox xbet,help={"beam tilt (pitch, forward/back)"},fSize=12,value= 1
	CheckBox xgam,pos={100,50},size={50,15},title="xgam",help={"beam tilt (yaw)"}
		CheckBox xgam,fSize=12,value= 1
	CheckBox dd,pos={25,70},size={32,15},title="dd"
		CheckBox dd,help={"detector tilt (height of detector center)"},fSize=12,value= 0
	CheckBox xc,pos={25,90},size={32,15},title="xc"
		CheckBox xc,help={"detector tilt (detector center X)"},fSize=12,value= 0
	CheckBox yc,pos={100,90},size={32,15},title="yc"
		CheckBox yc,help={"detector tilt (detector center Y)"},fSize=12,value= 0
	CheckBox dpsx,pos={25,110},size={46,15},title="dpsx"
		CheckBox dpsx,help={"size of detector X direction"},fSize=12,value= 0
	CheckBox dpsy,pos={100,110},size={46,15},title="dpsy"
		CheckBox dpsy,help={"size of detector Y direction"},fSize=12,value= 0
	Button refine,pos={33,138},size={100,20},proc=RefineGeoButtonProc,title="Refine"
End
//
Function RefineGeoButtonProc(ctrlName) : ButtonControl
	String ctrlName
	String hstr=""
	ControlInfo/W=RefinementPanel alphad
	hstr += num2istr(!V_Value)
	ControlInfo/W=RefinementPanel betad
	hstr += num2istr(!V_Value)
	ControlInfo/W=RefinementPanel xbet
	hstr += num2istr(!V_Value)
	ControlInfo/W=RefinementPanel xgam
	hstr += num2istr(!V_Value)
	ControlInfo/W=RefinementPanel dd
	hstr += num2istr(!V_Value)
	ControlInfo/W=RefinementPanel xc
	hstr += num2istr(!V_Value)
	ControlInfo/W=RefinementPanel yc
	hstr += num2istr(!V_Value)
	ControlInfo/W=RefinementPanel dpsx
	hstr += num2istr(!V_Value)
	ControlInfo/W=RefinementPanel dpsy
	hstr += num2istr(!V_Value)
	print "¥|"
	if (	stringmatch(StringFromLIst(0,WinList("RefinementPanel",";","WIN:64")),"RefinementPanel"))
		DoWindow/K RefinementPanel
	endif
	doGeoRefinement($"",hstr)
End

