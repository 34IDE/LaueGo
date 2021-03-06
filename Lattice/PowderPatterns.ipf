#pragma TextEncoding = "UTF-8"		// For details execute DisplayHelpTopic "The TextEncoding Pragma"
#pragma rtGlobals=3		// Use modern global access method.
#pragma version = 1.03
#pragma IgorVersion = 6.3
#pragma ModuleName=powder
#requiredPackages "LatticeSym;"
#initFunctionName "Init_PowderPatternLattice()"
#include "LatticeSym" version>=7.00


Menu "Analysis"
	SubMenu "Lattice"
		"-"
		"Calculate Powder Lines ...",CalcPowderLines(NaN)
		help={"Calculate the Powder Lines for a crystal structure"}
		MenuItemIfWaveClassExists("Make Powder Pattern from Lines ...","PowderLines*",""),PowderPatternFromLines($"",NaN)
		help={"Calculate an actual Powder Pattern from a set of Powder Lines"}
		MenuItemIfWaveClassExists("  Graph of Powder Pattern","PowderPattern*",""),GraphPowderPattern($"")
		help={"Show graph of Powder Pattern (with its Powder Lines)"}
		MenuItemIfWaveClassExists("  Graph of Powder Lines","PowderLines*",""),GraphPowderLines($"")
		help={"Show graph of Powder Lines"}
		MenuItemIfWaveClassExists("  Table of Powder Lines","PowderLines*",""),TablePowderLines($"")
		help={"Show table of Powder Lines"}
	End
End

Static Constant hc_keVnm = 1.239841856			// h*c (keV-nm)
Static Constant N_POWDER_LINES_COLUMNS = 8		// number of coumns in a powder lines wave


Static Function AfterFileOpenHook(refNum,file,pathName,type,creator,kind)
	Variable refNum, kind
	String file,pathName,type,creator
	if ((kind==1) || (kind==2) || (kind==12))		// an experiment (packed or unpacked), or an ipf file
		Init_PowderPatternLattice()
	endif
	return 0
End
Static Function IgorStartOrNewHook(IgorApplicationNameStr)
	String IgorApplicationNameStr
	Init_PowderPatternLattice()
	return 0
End


Static Function PowderPatternPopMenuProc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa
	if (pa.eventCode != 2)
		return 0
	endif

	if (strsearch(pa.popStr,"Calculate Powder Lines",0,2)>=0)
		printf "%sCalcPowderLines(NaN)\r", BULLET
		CalcPowderLines(NaN)
	elseif (strsearch(pa.popStr,"Make Powder Pattern from Lines",0,2)>=0)
		printf "%sPowderPatternFromLines($\"\",NaN)\r", BULLET
		PowderPatternFromLines($"",NaN)
	elseif (strsearch(pa.popStr,"Graph of Powder Pattern",0,2)>=0)
		printf "%sGraphPowderPattern($\"\")\r", BULLET
		print GraphPowderPattern($"")
	elseif (strsearch(pa.popStr,"Graph of Powder Lines",0,2)>=0)
		printf "%sGraphPowderLines($\"\")\r", BULLET
		print GraphPowderLines($"")
	elseif (strsearch(pa.popStr,"Table of Powder Lines",0,2)>=0)
		printf "%sTablePowderLines($\"\")\r", BULLET
		print TablePowderLines($"")
	endif
	return 0
End


// Create a powder pattern from a list of powder lines
Function/WAVE PowderPatternFromLines(lines,fwhmQ,[theta])
	Wave lines
	Variable fwhmQ
	Variable theta										// a flag, if true, then use theta for the x-axis
	fwhmQ = numtype(fwhmQ) || fwhmQ<=0 ? NaN : fwhmQ
	theta = ParamIsDefault(theta) || numtype(theta) ? 0 : !(!theta)

	Variable printIt = strlen(GetRTStackInfo(2))<1
	if (!WaveExists(lines) || numtype(fwhmQ))
		printIt = 1
		String options
		sprintf options "DIMS:2,MINCOLS:%d,MAXCOLS:%d", N_POWDER_LINES_COLUMNS,N_POWDER_LINES_COLUMNS
		String list = reverseList(WaveListClass("PowderLines*","*",options))
		Variable i = ItemsInLIst(list)
		i = WaveExists(lines) ? -1 : i
		if (i==0)
			DoAlert 0, "No waves of this type in current folder"
			return $""
		elseif (i==1)
			Wave lines = $StringFromList(0,list)
		endif

		if (!WaveExists(lines) || numtype(fwhmQ))
			fwhmQ = fwhmQ>0 ? fwhmQ : NumVarOrDefault("root:Packages:Lattices:Powder:defaultPeakFWHMQ",0.5)
			String wName = StringFromList(0,list)
			Prompt wName, "Wave with Powder Lines",popup,list
			Prompt fwhmQ,"FWHM of lines (1/nm)"
			DoPrompt "Powder Lines",wName,fwhmQ
			if (V_flag)
				return $""
			endif
			Wave lines = $wName
		endif
		if (!WaveExists(lines))
			DoAlert 0, "No Powder Lines"
			return $""
		endif
	endif
	fwhmQ = numtype(fwhmQ) || fwhmQ<=0 ? NaN : fwhmQ

	Variable keV=NaN
	if (WaveExists(lines))
		keV = NumberByKey("keV",note(lines),"=")
	endif

	theta = numtype(keV+lines[0][4]) ? 0 : theta	// may need to override input
	if (WaveExists(lines) && fwhmQ>0 && numtype(keV+lines[0][4])==0 && ParamIsDefault(theta))
		Prompt theta, "Q or Theta",popup,"Q on x-axis;Theta on x-axis"
		theta += 1
		DoPrompt "Choose X-axis",theta
		theta -= 1
		if (V_flag)
			return $""
		endif			
	endif

	if (printIt)
		printf "PowderPatternFromLines(%s,%g",NameOfWave(lines),fwhmQ
		if (!ParamIsDefault(theta) || theta)
			printf ", theta=%g", theta
		endif
		printf ")\r"
	endif
	if (!WaveExists(lines) || numtype(fwhmQ))
		return $""
	endif
	Variable/G root:Packages:Lattices:Powder:defaultPeakFWHMQ=fwhmQ		// save as default

	Variable Nlines=DimSize(lines,0)
	Variable Qwidth = lines[Nlines-1][0] + 2*fwhmQ	// width in Q(1/nm)
	Variable N = ceil(Qwidth/fwhmQ * 5)		// number of points in curve
	N = theta ? 2*N : N							// sample at higher resolution because we will interpolate back down

	wName = AddEndingToWaveName(NameOfWave(lines),SelectString(theta,"_Q","_Theta"))
	Make/N=(N)/O $wName/WAVE=intens = 0	// Wave to receive intensity
	Variable useF2 = numtype(lines[0][6])!=0
	String leftLabel = SelectString(useF2,"Intensity","|F|\\S2\\M")
	//	SetScale d 0,0,leftLabel, intens
	String class = "PowderPattern" + SelectString(theta,"Qnm","Theta")
	String wnote = ReplaceStringByKey("waveClass",note(lines),class,"=")
	wnote = ReplaceNumberByKey("FWHMQ",wnote,fwhmQ,"=")
	wnote = ReplaceNumberByKey("useF2",wnote,useF2,"=")
	wnote = ReplaceStringByKey("PowderLines",wnote,GetWavesDataFolder(lines,2),"=")
	wnote = ReplaceStringByKey("leftLabel",wnote,leftLabel,"=")
	Variable dQ = Qwidth/(N-1)
	String xUnits = SelectString(theta,"1/nm",DEGREESIGN)

	// first calculate in Q whether we want Q or theta
	Variable m,in, a = 4*ln(2)/(fwhmQ*fwhmQ)	// gaussian = exp(-a*x^2) in Q (not theta)
	for (m=0;m<Nlines;m+=1)
		SetScale/P x -lines[m][0],dQ,"", intens	// rescale so line[m] is zero at central Q
		in = useF2 ? lines[m][5] : lines[m][6]
		intens += in * exp(-a * x^2)
	endfor
	SetScale/I x 0,Qwidth,xUnits, intens			// re-set x to correct Q scaling

	if (theta)
		Duplicate/FREE intens, intensQ
		SetScale/I x 0,Qwidth,"", intensQ			// set to Q scaling
		Variable sine = limit(Qwidth * (hc_keVnm/keV) / (4*PI),0,1)
		Variable width = asin(sine) * 180/PI		// width in degrees
		N /= 2											// we used an N that was too big before, so reduce it here
		Redimension/N=(N) intens
		Make/N=(N)/FREE  Qs
		SetScale/I x 0,width,xUnits, intens, Qs	// set to Theta scaling
		Qs = 4*PI*sin(x*PI/180)*keV/hc_keVnm		// Q at each angle
		intens = intensQ(Qs[p])						// fill array at each theta
	endif

	Note/K intens, wnote
	if (printIt)
		printf "Created Powder Pattern in '%s',  using  %d  lines,  with FWHM = %g (1/nm)",NameOfWave(intens),Nlines,fwhmQ
		if (theta)
			printf ",  Using Theta scaling"
		endif
		printf "\r"
	endif
	GraphPowderPattern(intens)						// create  the graph or bring it to the front
	return intens
End

Function/WAVE CalcPowderLines(Qmax,[keV,Polarization,scaling])
	Variable Qmax					// max Q (1/nm)
	Variable keV					// x-ray energy (keV)
	Variable Polarization		// 1=sigma, 0=unpolarized, -1=pi, or anything in the range [-1,1]
	Variable scaling				// 0=traditional (100 max),  1=un-scaled
	keV = ParamIsDefault(keV) ? NaN : keV
	Polarization = ParamIsDefault(Polarization) ? 0 : Polarization		// default to Un-Polarized
	Polarization = numtype(Polarization) ? 0 : Polarization
	Polarization = limit(Polarization,-1,1)
	scaling = ParamIsDefault(scaling) || numtype(scaling) ? NumVarOrDefault("root:Packages:Lattices:Powder:intensityScaling",0) : scaling

	Variable printIt = (Qmax<=0 || numtype(Qmax)) || strlen(GetRTStackInfo(2))<1
	if (Qmax<=0 || numtype(Qmax))
		Prompt Qmax,"Maximum Q in powder pattern (1/nm)"
		Prompt keV,"Energy of Beam (keV), OPTIONAL"
		Prompt Polarization,"Polarization",popup,"Sigma (Synchrotrons);Un-Polarized (Tube Source);Pi;User Supplied..."
		Prompt scaling, "Intensity Scaling", popup, "Traditional (max=100);UN-scaled"
		Qmax = Qmax>0 ? Qmax : NumVarOrDefault("root:Packages:Lattices:Powder:defaultQmax",100)
		keV = keV>0 ? keV : NumVarOrDefault("root:Packages:Lattices:Powder:defaultKeV",NaN)
		if (Polarization==-1)
			Polarization = 3
		elseif (Polarization==0)
			Polarization = 2
		elseif (Polarization==1)
			Polarization = 1
		else
			Polarization = 4
		endif	
		scaling += 1
		DoPrompt "Qmax",Qmax,keV,Polarization,scaling
		if (V_flag)
			return $""
		endif
		if (Polarization==1)
			Polarization = 1
		elseif (Polarization==2)
			Polarization = 0
		elseif (Polarization==3)
			Polarization = -1
		elseif (Polarization==4)
			Polarization = NaN
		endif	
		if (!(Polarization>=-1 && Polarization<=1))
			Prompt Polarization,"Enter Polarization in range [-1,1]"
			DoPrompt "Polarization", Polarization
			if (V_flag)
				return $""
			endif
		endif
		scaling -= 1
	endif
	scaling = scaling ? 1 : 0
	if (printIt)
		printf "CalcPowderLines(%g",Qmax
		if (!ParamIsDefault(keV) || numtype(keV)!=2)
			printf ", keV=%g",keV
		endif
		printf ", Polarization=%g",Polarization
		if (!ParamIsDefault(scaling) || scaling!=0)
			printf ", scaling=%g",scaling
		endif
		printf ")\r"
	endif
	if (Qmax<=0 || numtype(Qmax) || !(Polarization>=-1 && Polarization<=1))
		return $""
	endif
	keV = keV<=0 || numtype(keV) ? NaN : keV
	Variable lambda = hc_keVnm/keV
	if (lambda>0)
		Qmax = min(Qmax,4*PI/lambda)				// with keV, Q is limited to theta=90 degree
	endif
	Variable/G root:Packages:Lattices:Powder:defaultQmax=Qmax
	Variable/G root:Packages:Lattices:Powder:defaultKeV=keV
	Variable/G root:Packages:Lattices:Powder:intensityScaling=scaling

	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))					//fill the lattice structure with 'current' values
		DoAlert 0, "ERROR  -- CalcPowderLines()\rNo lattice structure found"
		return $""
	endif
	Variable dim = xtal.dim
	dim = dim==2 ? 2 : 3
	String desc = xtal.desc
	desc = ReplaceString(";",desc," ")
	Variable astar, bstar, cstar=NaN
	if (dim==2)
		astar = sqrt((xtal.as0)^2 + (xtal.as1)^2)
		bstar = sqrt((xtal.bs0)^2 + (xtal.bs1)^2)
	else
		astar = sqrt((xtal.as0)^2 + (xtal.as1)^2 + (xtal.as2)^2)
		bstar = sqrt((xtal.bs0)^2 + (xtal.bs1)^2 + (xtal.bs2)^2)
		cstar = sqrt((xtal.cs0)^2 + (xtal.cs1)^2 + (xtal.cs2)^2)
	endif
	Variable hmax=floor(Qmax/astar), kmax=floor(Qmax/bstar), lmax=floor(Qmax/cstar)
	lmax = dim==2 ? 0 : lmax

	String wNote = "waveClass=PowderLines;"
	wNote = ReplaceNumberByKey("Qmax",wNote,Qmax,"=")
	if (keV>0)
		wNote = ReplaceNumberByKey("keV",wNote,keV,"=")
	endif
	String xtalTitle = desc + SelectString(strlen(xtal.SpaceGroupID), "", "  "+xtal.SpaceGroupID)
	xtalTitle += SelectString(numtype(xtal.SpaceGroupIDnum), "", "  "+getHMsym2(xtal.SpaceGroupIDnum), "")
	wNote = ReplaceStringByKey("xtalTitle",wNote,xtalTitle,"=")

	Variable NlinesMax=200
	Make/N=(NlinesMax,7)/D/FREE lines=NaN

	Variable F2min = magSqr( FstructMax(xtal,Qmax,keV=keV) )/ 50
	F2min = numtype(F2min) || F2min<=0 ? 0.001 : F2min
	Variable qhkl,Pol=NaN,Lorentz=NaN,theta=NaN,F2
	Variable h,k,l, m
	for (m=0,l=0; l<=lmax; l=LatticeSym#incrementIndex(l))
		for (k=0; k<=kmax; k=LatticeSym#incrementIndex(k))
			for (h=0; h<=hmax; h=LatticeSym#incrementIndex(h))
				qhkl = 2*PI/dSpacing(xtal,h,k,l)
				if (qhkl>Qmax || qhkl<=0)
					continue
				endif
				if (keV > 0)
					theta = asin((qhkl*lambda/(4*PI)))
					Lorentz = 2/ (sin(theta)*sin(2*theta))
					Pol = (1+Polarization)/2 + (1-Polarization)/2 * cos(2*theta)^2
					// LP = (1+cos(2*theta)^2) / (sin(theta)*sin(2*theta))
					F2 = magsqr(Fstruct(xtal,h,k,l,keV=keV))
				else
					F2 = magsqr(Fstruct(xtal,h,k,l))
				endif
				if (F2 < F2min)			// was F2 < 0.001
					continue
				endif
				if (m>=NlinesMax)
					NlinesMax += 200
					Redimension/N=(NlinesMax,-1) lines
				endif
				lines[m][0] = qhkl
				lines[m][1] = h
				lines[m][2] = k
				lines[m][3] = dim==2 ? NaN : l
				lines[m][4] = theta*180/PI			// save angle in degrees
				lines[m][5] = F2
				lines[m][6] = F2 * Lorentz * Pol	// was F2*LP, have yet to include mult
				m += 1
			endfor
		endfor
	endfor
	NlinesMax = m											// the largest multiplicity
	Redimension/N=(NlinesMax,-1) lines				// trim off extras

	Make/N=(m)/FREE/D Qs, QsIndex, temp
	Qs = lines[p][0]										// Q values for every possible hkl
	MakeIndex Qs, QsIndex								//  to sort with increasing Q
	Variable i
	for (i=0;i<DimSize(lines,1);i+=1)
		temp = lines[p][i]								// for each column in lines
		IndexSort QsIndex, temp						// sort this column
		lines[][i] = temp[p]							// replace column in lines with sorted column
	endfor
	WaveClear temp
	Redimension/N=(NlinesMax+1,-1) lines			// add point at end to stop next loop
	lines[NlinesMax][0] = NaN							// ensures that following do-while loop will end

	String str = "_"+num2istr(Qmax)+"nm"
	String wname = CleanupName(desc,0)
	wname = ReplaceString("__",wname,"_")
	wname = ReplaceString("__",wname,"_")
	wname = ReplaceString("__",wname,"_")
	wname = wname[0,31-strlen(str)-1] + str
	Make/N=(NlinesMax,N_POWDER_LINES_COLUMNS)/D/O $wname/WAVE=lineShort = NaN
	// columns are:
	//	0	Q (1/nm)
	//	1	h
	//	2	k
	//	3	l							not used for 2D
	//	4	theta (degree)
	//	5	| F| ^2
	//	6	Intensity
	//	7	multiplicity
	SetDimLabel 1,0,Q_nm,lineShort ;	SetDimLabel 1,1,h,lineShort ;		SetDimLabel 1,2,k,lineShort ;			SetDimLabel 1,3,l,lineShort
	SetDimLabel 1,4,theta,lineShort ;	SetDimLabel 1,5,F2,lineShort ;		SetDimLabel 1,6,Intensity,lineShort ;	SetDimLabel 1,7,mult,lineShort

	str = wname[0,31-10-1]+"_hklStr"
	Make/N=(NlinesMax)/T/O $str/WAVE=hklStr = ""
	Note/K hklStr, "waveClass=HKLlabelsPowderLines;LinesWave="+GetWavesDataFolder(lineShort,2)+";"
	wNote = ReplaceStringByKey("hklStrWave",wNote,str,"=")

	str += "All"
	Make/N=(NlinesMax)/T/O $str/WAVE=hklStrAll = ""
	Note/K hklStrAll, "waveClass=HKLlabelsPowderLinesAll;LinesWave="+GetWavesDataFolder(lineShort,2)+";"
	wNote = ReplaceStringByKey("hklStrAllWave",wNote,str,"=")
	wNote = ReplaceNumberByKey("dim",wNote,dim,"=")

	Make/N=(NlinesMax,dim)/FREE hkls					// temp list of hkls for line, used to pick the "nicest" hkl
	Variable mult, Nlines=0
	for (Nlines=0,i=0; i<NlinesMax;)
		qhkl = lines[i][0]									// Q of this line
		mult = 0
		hkls = NaN
		do
			lineShort[Nlines][0,N_POWDER_LINES_COLUMNS-2] = lines[i][q]
			hkls[mult][] = lines[i][q+1]					// accumulate list of all hkls in this line
			hklStrAll[Nlines] += "("+hkl2str(lines[i][1],lines[i][2],lines[i][3],dim=dim)+")+"
			mult += 1
			i += 1
		while(abs(lines[i][0]-qhkl)<1e-10)
		hklStrAll[Nlines] = ReplaceString(", ",hklStrAll[Nlines],",")
		hklStrAll[Nlines] = TrimBoth(hklStrAll[Nlines],chars="+")
		lineShort[Nlines][7] = mult
		Wave niceHKL = nicestHKL(hkls)
		lineShort[Nlines][1,dim] = niceHKL[q-1]		// re-set hkl to the "nicest" one
		hklStr[Nlines] = nicestHKLstr(hkls)			// better for labeling, it includes multipliticites
		Nlines += 1
	endfor
	Redimension/N=(Nlines,-1) lineShort, hklStr, hklStrAll		// trim to final size
	lineShort[][6] = lineShort[p][6] * lineShort[p][7]	// increas intens by multiplicity

	if (scaling==0)											// traditional scaling, a max of 100
		MatrixOP/FREE mmm = 100 / maxVal(col(lineShort,6))
		Variable scale100 = mmm[0]
		if (numtype(scale100)==0 && scale100>0)
			lineShort[][6] *= scale100
		endif
	endif
	wNote = ReplaceStringByKey("scaling",wNote,SelectString(scaling,"Traditional(max=100)","UN-scaled"),"=")
	Note/K lineShort, wNote

	if (printIt)
		printf "Created Powder Lines in '%s',  with %d lines,  range of hkl is [%s%d, %s%d, %s%d]\r",wname,Nlines,PLUSMINUS,hmax,PLUSMINUS,kmax,PLUSMINUS,lmax
	endif
	return lineShort
End
//
Static Function/WAVE nicestHKL(hklsIn)		// returns a free wave that is the nicest one from the list of hkl's
	// The nicest hkl has L>K>H, as much as possible
	Wave hklsIn
	Variable dim=DimSize(hklsIn,1)
	if (DimSize(hklsIn,0)<1 || (dim!=2 && dim!=3))
		return $""
	endif
	WaveStats/Q/M=1 hklsIn
	Variable N=V_npnts/dim
	if (N<1)
		return $""
	endif

	Make/N=(N,dim)/FREE hkls=hklsIn[p][q]
	if (dim!=2)
		ImageStats/M=1/G={0,N-1,2,2} hkls
		hkls = hkls[p][2]==V_max ? hkls[p][q] : NaN// remove all hkl with L<V_max
	endif
	ImageStats/M=1/G={0,N-1,1,1} hkls
	hkls = hkls[p][1]==V_max ? hkls[p][q] : NaN	// remove all hkl with K<V_max
	ImageStats/M=1/G={0,N-1,0,0} hkls
	hkls = hkls[p][0]==V_max ? hkls[p][q] : NaN	// remove all hkl with H<V_max

	ImageStats/M=1/G={0,N-1,0,dim-1} hkls				// find the best (probably only) one left
	Make/N=(dim)/FREE hkl=hkls[V_maxRowLoc][p]
	return hkl
End
//Function test_nicestHKL()
//	Make/N=(5,3)/FREE hkls
//	hkls[0][0]= {2,0,1,0,1}
//	hkls[0][1]= {1,-2,2,1,0}
//	hkls[0][2]= {0,1,0,2,2}
//	Wave hkl = nicestHKL(hkls)
//	print "hkl = ",vec2str(hkl)
//End


Static Function/T nicestHKLstr(hklsIn)		// returns a free wave that is the nicest one from the list of hkl's
	// The nicest hkl has L>K>H, as much as possible
	Wave hklsIn
	Variable dim=DimSize(hklsIn,1)
	if (DimSize(hklsIn,0)<1 || (dim!=2 && dim!=3))
		return ""
	endif
	WaveStats/Q/M=1 hklsIn
	Variable N=V_npnts/dim
	if (N<1)
		return ""
	endif
	Make/N=(N,dim)/FREE hkls=hklsIn[p][q]

	if (dim>2)
		ImageStats/M=1/G={0,N-1,2,2} hkls
		hkls = hkls[p][2]==V_max ? hkls[p][q] : NaN	// remove all hkl with L<V_max
	endif
	ImageStats/M=1/G={0,N-1,1,1} hkls
	hkls = hkls[p][1]==V_max ? hkls[p][q] : NaN		// remove all hkl with K<V_max
	ImageStats/M=1/G={0,N-1,0,0} hkls
	hkls = hkls[p][0]==V_max ? hkls[p][q] : NaN		// remove all hkl with H<V_max

	ImageStats/M=1/G={0,N-1,0,dim-1} hkls					// find the best (probably only) one left
	Make/N=(dim)/FREE hkl=hkls[V_maxRowLoc][p]
	Variable h=hkl[0], k=hkl[1],   l = (dim==2) ? NaN : hkl[2]
	String hklStr = "("+hkl2str(h,k,l,dim=dim)+")"

	// remove this hkl from hkls and do again
	hkls = hklsIn[p][q]

	MatrixOP/FREE checks = sumRows(abs(hkls))
	Redimension/N=(N) checks					// I am assuming that if the summs differ than the hkl is of different type
	hkl = abs(hkl)
	Variable i,m, check0=sum(hkl)
	checks = checks[p]!=check0				// keep all items with checks[]==1
	for (i=0,m=0; i<N; i+=1)
		if (checks[i])
			hkls[m][] = hkls[i][q]
			m += 1
		endif
	endfor
	N = m
	if (N>0)
		Redimension/N=(N,-1) hkls
		String str = nicestHKLstr(hkls)
		hklStr += SelectString(strlen(str),""," + "+str)
	endif
	return hklStr
End
//	Function test_nicestHKLstr()
//	//	Make/N=(5,3)/FREE hkls
//	//	hkls[0][0]= {2,0,1,0,1}
//	//	hkls[0][1]= {1,-2,2,1,0}
//	//	hkls[0][2]= {0,1,0,2,2}
//	
//		Make/N=(32,3)/FREE hkls
//		hkls[0][0]= {-1,1,3,-5,-3,-3,1,-1,1,1,-1,-3,-1,5,-1,5,5,-5,1,5,-5,1,-1,-3,3,-5,1,3,3,-1,-1,1}
//		hkls[0][1]= {5,-1,3,1,-3,3,1,1,5,-1,-1,3,5,1,-5,-1,-1,-1,-5,1,-1,1,-1,-3,-3,1,5,3,-3,1,-5,-5}
//		hkls[0][2]= {1,5,3,-1,3,-3,-5,5,-1,-5,-5,3,-1,1,1,-1,1,-1,-1,-1,1,5,5,-3,-3,1,1,-3,3,-5,-1,1}
//	
//	//	Make/N=(8,3)/FREE hkls
//	//	hkls[0][0]= {1,-1,-1,-1,1,1,-1,1}
//	//	hkls[0][1]= {-1,-1,1,-1,1,-1,1,1}
//	//	hkls[0][2]= {1,1,-1,-1,1,-1,1,-1}
//	
//		printf "hkl = '%s'\r",powder#nicestHKLstr(hkls)
//	End



Function GraphPowderPattern(ww)
	Wave ww

	if (!WaveExists(ww))
		String list = reverseList(WaveListClass("PowderPattern*","*",""))
		Variable i = ItemsInLIst(list)
		if (i<1)
			DoAlert 0, "No waves of this type in current folder"
			return 1
		elseif (i==1)
			Wave ww = $StringFromList(0,list)
		else
			String wName = StringFromList(0,list)
			Prompt wName, "Wave with Powder Pattern",popup,list
			DoPrompt "Powder Pattern to Graph",wName
			if (V_flag)
				return 1
			endif			
			Wave ww = $wName
		endif
		if (!WaveExists(ww))
			DoAlert 0, "No waves of name "+wName
			return 1
		endif
	endif

	String graph = StringFromList(0,WindowsWithWave(ww,1))
	if (strlen(graph)>0)
		DoWindow/F $graph
		return 0
	endif

	String wNote = note(ww), bottomLabel
	Wave lines = $StringByKey("PowderLines",wnote,"=")
	Variable useF2 = NumberByKey("useF2", wNote,"=")
	String leftLabel = StringByKey("leftLabel",wNote,"=")
	Variable iy = useF2 ? 5 : 6
	Variable ix = 0
	if (StringMatch(StringByKey("waveClass",note(ww),"="),"PowderPatternQnm"))
		ix = 0
	elseif (StringMatch(StringByKey("waveClass",note(ww),"="),"PowderPatternTheta"))
		ix = 4
	elseif (StringMatch(WaveUnits(ww,0),"1/nm"))
		ix = 0
	elseif (StringMatch(WaveUnits(ww,0),DEGREESIGN))
		ix = 4
	endif
#if (IgorVersion()<7)
	if (strlen(WaveUnits(ww,-1))==0)
		bottomLabel = SelectString(ix,"Q  (\\Enm\\S-1\\M)", "\\Zr150\\F'Symbol'q\260\\F]0 \\E\\M")	// "\260" = option-5, degree sign in Symbol font
	else
		bottomLabel = SelectString(ix,"Q  (\\U)", "\\Zr150\\F'Symbol'q\\F]0\\U\\M")
	endif
#else
	if (strlen(WaveUnits(ww,-1))==0)
		bottomLabel = SelectString(ix,"Q  (\\Enm\\S-1\\M)", LetterName2Unicode("theta")+DEGREESIGN+" \\E\\M")
	else
		bottomLabel = SelectString(ix,"Q  (\\U)", LetterName2Unicode("theta")+"\\U\\M")
	endif
#endif

	Display /W=(92,274,1135,578) ww
	if (WaveExists(lines))
		AppendToGraph	lines[*][iy] vs lines[*][ix]
		String lname = NameOfWave(lines)
		ModifyGraph mode($lname)=1, rgb($lname)=(1600,1600,65535)
	endif
	ModifyGraph tick=2, mirror=1, minor=1, lowTrip=0.001

	Label left leftLabel
	Label bottom bottomLabel
	SetAxis/A/N=1/E=1

	Variable keV = NumberByKey("keV",wnote,"=")
	Variable FWHMQ = NumberByKey("FWHMQ",wnote,"=")
	String title = "\\JR\\Zr125"
	title += StringByKey("xtalTitle", wnote,"=")
	title += SelectString(keV>0 || FWHMQ>0,"","\\Zr075")
	title += SelectString(keV>0,"","\r"+num2str(keV)+" keV")
	title += SelectString(FWHMQ>0,"","\rFWHM = "+num2str(FWHMQ)+"(nm\\S-1\\M)")
	TextBox/C/N=textTitle/F=0/B=1 title
	SetWindow kwTopWin,hook(powderLineHook)=powder#ShowPowderLinesWindowHook
End


Function GraphPowderLines(ww,[theta])
	Wave ww
	Variable theta					// flag, if true, show theta rather than Q
	theta = ParamIsDefault(theta) || numtype(theta) ? 0 : !(!theta)

	if (!WaveExists(ww))
		String list = reverseList(WaveListClass("PowderLines*","*",""))
		Variable i = ItemsInLIst(list)
		if (i<1)
			DoAlert 0, "No waves of this type in current folder"
			return 1
		elseif (i==1)
			Wave ww = $StringFromList(0,list)
		else
			String wName = StringFromList(i-1,list)
			Prompt wName, "Wave with Powder Lines",popup,list
			DoPrompt "Powder Lines",wName
			if (V_flag)
				return 1
			endif			
			Wave ww = $wName
		endif
		if (!WaveExists(ww))
			DoAlert 0, "No waves of name "+wName
			return 1
		endif
	endif

	String graph = StringFromList(0,WindowsWithWave(ww,1))
	if (strlen(graph)>0)
		DoWindow/F $graph
		return 0
	endif

	if (ParamIsDefault(theta) && numtype(ww[0][4])==0)
		Prompt theta, "Q or Theta",popup,"Q on x-axis;Theta on x-axis"
		theta += 1
		DoPrompt "Choose X-axis",theta
		theta -= 1
		if (V_flag)
			return 1
		endif			
	endif
	Variable iy=numtype(ww[0][6]) ? 5 : 6
	Variable ix = (theta && numtype(ww[0][4])==0) ? 4 : 0
	String leftLabel = SelectString(iy==6,"|F|\\S2\\M","Intensity")
#if (IgorVersion()<7)
	String bottomLabel = SelectString(ix,"Q  (nm\\S-1\\M)", "\\Zr150\\F'Symbol'q\260\\F]0\\M")	// "\260" = option-5, degree sign in Symbol font
#else
	String bottomLabel = SelectString(ix,"Q  (nm\\S-1\\M)", LetterName2Unicode("theta")+DEGREESIGN+"\\M")
#endif

	String wnote = note(ww)
	Variable keV = NumberByKey("keV",wnote,"=")
	String title = "\\JR\\Zr125"
	title += StringByKey("xtalTitle", wnote,"=")
	title += SelectString(keV>0,"","\r"+num2str(keV)+"\\M keV")

	Display /W=(273,348,1151,732) ww[*][iy] vs ww[*][ix]
	ModifyGraph mode=1, tick=2, mirror=1, minor=1, lowTrip=0.001, rgb=(1600,1600,65535)
	Label left leftLabel
	Label bottom bottomLabel
	TextBox/C/N=textTitle/B=1/F=0 title
	SetAxis/A/N=1/E=1 
	SetWindow kwTopWin,hook(powderLineHook)=powder#ShowPowderLinesWindowHook
	return 0
End


Static Function ShowPowderLinesWindowHook(s)
	STRUCT WMWinHookStruct &s

//	if ((s.eventMod==3) && (s.eventCode ==3  || s.eventCode == 4))	// mouse down or moved with Shift key held down
	if (s.eventMod==3)	// mouse down and Shift key held down
		String win = s.winName
		String wList=TraceNameList(win,";",1)
		Wave lines=$""
		Variable i
		for (i=0;i<ItemsInList(wList);i+=1)
			Wave lines=TraceNameToWaveRef(win,StringFromList(i,wList))
			if (WaveInClass(lines,"PowderLines"))
				break
			endif
			Wave lines=$""
		endfor
		if (!WaveExists(lines))
			return 0
		endif
		String wNote=note(lines)
		Wave/T hklStr = $(GetWavesDataFolder(lines,1)+StringByKey("hklStrWave",wNote,"="))
		Wave/T hklStrAll = $(GetWavesDataFolder(lines,1)+StringByKey("hklStrWave",wNote,"=")+"All")

		Variable dim=NumberByKey("dim",wNote,"=")
		dim = dim==2 ? 2 : 3
		Variable Xval = AxisValFromPixel(win,"bottom",s.mouseLoc.h)
		String xrange=StringByKey("XRANGE",TraceInfo(win,NameOfWave(lines),0))
		Variable isTheta=StringMatch(xrange,"*[4]")		// [0] is Q, [4] is theta
		Variable keV=NumberByKey("keV",wNote,"=")
		isTheta = isTheta && keV>0

		Variable iaxis = isTheta ? 4 : 0
		Make/N=(DimSize(lines,0))/FREE delta = abs(lines[p][iaxis] - Xval)
		WaveStats/Q/M=1 delta
		WaveClear delta
		Variable iLine = V_minloc
		Xval = lines[iLine][iaxis]							// actual X value of the line, not the mouse
  		Variable Yval = lines[iLine][6]

		String tagStr, anchor, str
		Variable Qval=lines[iLine][0], theta=lines[iLine][4], mult=lines[iLine][7]
		if (WaveExists(hklStr))
			tagStr = hklStr[iLine]
		else
			tagStr = "("+hkl2str(lines[iLine][1],lines[iLine][2],lines[iLine][3], dim=dim)+")"
		endif
//		tagStr += "\\Zr075\r"
		tagStr += "\\Zr085\r"
		if (WaveExists(hklStrAll))
			if (ItemsInList(hklStrAll[iLine],"+") <= 4)
				tagStr += hklStrAll[iLine]+"\r"
			endif
		endif
//		sprintf str, "Q = %g nm\\S-1\M\\Zr075",Qval
		sprintf str, "Q = %g nm\\S-1\M\\Zr085",Qval
		tagStr += str
		if (theta>0)
			sprintf str, "\rtheta = %g"+DEGREESIGN,theta
			tagStr += str
		endif
		sprintf str, "\rmult = %g",mult
		tagStr += str

		GetWindow/Z $win  psize
		Variable X0 = Xval>AxisValFromPixel(win,"bottom",(V_left+V_right)*0.5) ? -3 : 3
//		Variable Y0 = Yval>AxisValFromPixel(win,"left",(V_top+V_bottom)*0.5) ? -7 : 7
		Variable Y0 = Yval>AxisValFromPixel(win,"left",(V_top+V_bottom)*0.3) ? -7 : 7
		anchor = SelectString(X0>0,"R","L")
		anchor += SelectString(Y0>0,"T","B")
#if (IgorVersion()<7)
		Tag/C/N=lineHKL/W=$win/A=$anchor/F=0/L=2/X=(X0)/Y=(Y0)/P=1 $NameOfWave(lines),iLine,tagStr
#else
		Tag/C/N=lineHKL/W=$win/A=$anchor/F=0/L=2/X=(X0)/Y=(Y0)/G=(0,0,0,45875)/B=1/P=1 $NameOfWave(lines),iLine,tagStr
#endif
		DoUpdate
		return 1
	else											// clean up
		Tag/K/N=lineHKL/W=$(s.winName)
	endif
	return 0										// 0 if nothing done, else 1
End



Function TablePowderLines(w)
	Wave w
	if (!WaveExists(w))
		String list = reverseList(WaveListClass("PowderLines*","*",""))
		Variable i = ItemsInLIst(list)
		if (i<1)
			DoAlert 0, "No waves of this type in current folder"
			return 1
		elseif (i==1)
			Wave w = $StringFromList(0,list)
		else
			String wName = StringFromList(i-1,list)
			Prompt wName, "Wave with Powder Lines",popup,list
			DoPrompt "Powder Lines",wName
			if (V_flag)
				return 1
			endif			
			Wave w = $wName
		endif
		if (!WaveExists(w))
			DoAlert 0, "No waves of name "+wName
			return 1
		endif
	endif

	String table = StringFromList(0,WindowsWithWave(w,2))
	if (strlen(table)>0)
		DoWindow/F $table
		return 0
	endif
	String fontName=GetDefaultFont("")
	Variable fontSize = 12
	//	Variable right = WaveExists(hklStr) ? 709 : 588
	Variable right = 588

	Wave/T hklStr = $(GetWavesDataFolder(w,1)+StringByKey("hklStrWave",note(w),"="))
	Variable hklsWidth = 0
	if (WaveExists(hklStr))
		Make/N=(DimSize(hklStr,0))/I/FREE lengths = FontSizeStringWidth(fontName,fontSize,0,hklStr[p])
		hklsWidth = limit(4+WaveMax(lengths), 50, 200)
		right = right + hklsWidth
	endif

	Variable top=44
	Variable height=TableHeight(w,top,fontName,fontSize)
	Edit/K=1/W=(5,top,right,top+height) w.ld
	ModifyTable width(Point)=1,title(w.d)="Powder Lines",width(w.l)=20, format(w.d)=3,width(w.d)=68
	if (WaveExists(hklStr) && hklsWidth)
		AppendToTable hklStr.y
		ModifyTable width(hklStr.y)=hklsWidth, title(hklStr.d)="hkl's"
	endif
	ModifyTable font=fontName, size=fontSize
End

Static Function TableHeight(ww,top,fontName,fontSize)
	Wave ww
	Variable top
	String fontName
	Variable fontSize
	if (strlen(fontName)<1)
		fontName = GetDefaultFont("")
	endif
	fontSize = fontSize<6 || fontSize>64 || numtype(fontSize) ? 12 : fontSize

	String screen1=StringByKey("SCREEN1",IgorInfo(0))
	Variable i = strsearch(screen1,"RECT=",0)
	String rect = StringByKey("RECT",screen1[i,Inf],"=")
	Variable scrHeight=str2num(StringFromList(3,rect,","))
	scrHeight = numtype(scrHeight) || scrHeight<600 ? 600 : scrHeight
	Variable Nr=DimSize(ww,0), Nc=DimSize(ww,1)
	Variable height = 86 + (FontSizeHeight(fontName,fontSize,0))*Nr

	Variable hasColLabels
	for (i=0,hasColLabels=0; i<Nc; i+=1)
		hasColLabels = hasColLabels || strlen(GetDimLabel(ww,1,i))
	endfor
	if (hasColLabels)
		height +=  Nr*4								// leave room for column labels
	endif

	height =  min(scrHeight-top-5, height)	// make it fit on the screen
	return height
End



Function Init_PowderPatternLattice()
	if (!DataFolderExists("root:Packages:Lattices:Powder") && exists("Get_f")!=6)		// do not have Cromer-Liberman, warn the user
		DoAlert/T="Cromer-Liberman not Loaded" 1, "Cromer-Liberman NOT Loaded.\rFor accurate line intensities include Cromer-Liberman?"
		if (V_flag==1)
			Execute/P "INSERTINCLUDE  \"CromerLiberman\", version>=1.7";Execute/P "COMPILEPROCEDURES "
		endif
	endif
	NewDataFolder/O root:Packages
	NewDataFolder/O root:Packages:Lattices
	NewDataFolder/O root:Packages:Lattices:Powder
	if (exists("root:Packages:Lattices:Powder:defaultPeakFWHMQ")!=2)
		Variable/G root:Packages:Lattices:Powder:defaultPeakFWHMQ=0.5
	endif
	if (exists("root:Packages:Lattices:Powder:defaultQmax")!=2)
		Variable/G root:Packages:Lattices:Powder:defaultQmax=100
	endif
	if (exists("root:Packages:Lattices:Powder:defaultKeV")!=2)
		Variable/G root:Packages:Lattices:Powder:defaultKeV=NaN
	endif
	if (exists("root:Packages:Lattices:Powder:intensityScaling")!=2)
		Variable/G root:Packages:Lattices:Powder:intensityScaling=0		// 0=traditional(100) ,1=un-scaled
	endif
End
