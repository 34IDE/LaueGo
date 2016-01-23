#pragma rtGlobals=3		// Use modern global access method.
#pragma version = 0.16
#pragma IgorVersion = 6.3
#pragma ModuleName=powder
#requiredPackages "LatticeSym;"
#initFunctionName "Init_PowderPatternLattice()"


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

Static Constant hc_keVnm = 1.239841856				// h*c (keV-nm)
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
		printf "¥CalcPowderLines(NaN)\r"
		CalcPowderLines(NaN)
	elseif (strsearch(pa.popStr,"Make Powder Pattern from Lines",0,2)>=0)
		printf "¥PowderPatternFromLines($\"\",NaN)\r"
		PowderPatternFromLines($"",NaN)
	elseif (strsearch(pa.popStr,"Graph of Powder Pattern",0,2)>=0)
		printf "¥GraphPowderPattern($\"\")\r"
		print GraphPowderPattern($"")
	elseif (strsearch(pa.popStr,"Graph of Powder Lines",0,2)>=0)
		printf "¥GraphPowderLines($\"\")\r"
		print GraphPowderLines($"")
	elseif (strsearch(pa.popStr,"Table of Powder Lines",0,2)>=0)
		printf "¥TablePowderLines($\"\")\r"
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
	String leftLabel = SelectString(useF2,"Intensity","| F |\\S2\\M")
	//	SetScale d 0,0,leftLabel, intens
	String class = "PowderPattern" + SelectString(theta,"Qnm","Theta")
	String wnote = ReplaceStringByKey("waveClass",note(lines),class,"=")
	wnote = ReplaceNumberByKey("FWHMQ",wnote,fwhmQ,"=")
	wnote = ReplaceNumberByKey("useF2",wnote,useF2,"=")
	wnote = ReplaceStringByKey("PowderLines",wnote,GetWavesDataFolder(lines,2),"=")
	wnote = ReplaceStringByKey("leftLabel",wnote,leftLabel,"=")
	Variable dQ = Qwidth/(N-1)
	String xUnits = SelectString(theta,"1/nm","¡")

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
		SetScale/I x 0,Qwidth,"", intensQ		// set to Q scaling
		Variable width = asin(Qwidth * (hc_keVnm/keV) / (4*PI)) * 180/PI	// width in degrees
		N /= 2									// we used an N that was too big before, so reduce it here
		Redimension/N=(N) intens
		Make/N=(N)/FREE  Qs
		SetScale/I x 0,width,xUnits, intens, Qs	// set to Theta scaling
		Qs = 4*PI*sin(x*PI/180)*keV/hc_keVnm	// Q at each angle
		intens = intensQ(Qs[p])					// fill array at each theta
	endif

	Note/K intens, wnote
	if (printIt)
		printf "Created Powder Pattern in '%s',  using  %d  lines,  with FWHM = %g (1/nm)",NameOfWave(intens),Nlines,fwhmQ
		if (theta)
			printf ",  Using Theta scaling"
		endif
		printf "\r"
	endif
	GraphPowderPattern(intens)				// create  the graph or bring it to the front
	return intens
End
//
Static Function/T AddEndingToWaveName(wName,waveNameEnd)
	// create a new wave name from wav but with a new ending, needed due to maxIgorNameLen
	String wName
	String waveNameEnd

	String path=""								// a possible path preceeding the name
	Variable i=strsearch(wName,":",Inf,1)
	if (i>=0)
		path = wName[0,i]						// strip off the path and save it
		wName = wName[i+1,Inf]				// wName, now only the name part (no path)
	endif
	wName = wName[0,maxIgorNameLen-strlen(waveNameEnd)-1]
	wName = path + wName + waveNameEnd	// reassemble the full name
	return wName
End

Function/WAVE CalcPowderLines(Qmax,[keV,Polarization])
	Variable Qmax				// max Q (1/nm)
	Variable keV				// x-ray energy (keV)
	Variable Polarization		// 1=sigma, 0=unpolarized, -1=pi, or anything in the range [-1,1]
	keV = ParamIsDefault(keV) ? NaN : keV
	Polarization = ParamIsDefault(Polarization) ? 0 : Polarization		// default to Un-Polarized
	Polarization = numtype(Polarization) ? 0 : Polarization
	Polarization = limit(Polarization,-1,1)

	Variable printIt = (Qmax<=0 || numtype(Qmax)) || strlen(GetRTStackInfo(2))<1
	if (Qmax<=0 || numtype(Qmax))
		Prompt Qmax,"Maximum Q in powder pattern (1/nm)"
		Prompt keV,"Energy of Beam (keV), OPTIONAL"
		Prompt Polarization,"Polarization",popup,"Sigma (Synchrotrons);Un-Polarized (Tube Source);Pi;User Supplied..."
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
		DoPrompt "Qmax",Qmax,keV,Polarization
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
	endif
	if (printIt)
		printf "CalcPowderLines(%g",Qmax
		if (!ParamIsDefault(keV) || numtype(keV)!=2)
			printf ", keV=%g",keV
		endif
		printf ", Polarization=%g",Polarization
		printf ")\r"
	endif
	if (Qmax<=0 || numtype(Qmax) || !(Polarization>=-1 && Polarization<=1))
		return $""
	endif
	keV = keV<=0 || numtype(keV) ? NaN : keV
	Variable lambda = hc_keVnm/keV
	if (lambda>0)
		Qmax = min(Qmax,4*PI/lambda)				// with keV, Q is limited to theta=90¡
	endif
	Variable/G root:Packages:Lattices:Powder:defaultQmax=Qmax
	Variable/G root:Packages:Lattices:Powder:defaultKeV=keV

	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))					//fill the lattice structure with 'current' values
		DoAlert 0, "ERROR  -- CalcPowderLines()\rNo lattice structure found"
		return $""
	endif
	String desc = xtal.desc
	desc = ReplaceString(";",desc," ")
	Variable astar, bstar, cstar
	astar = sqrt((xtal.as0)^2 + (xtal.as1)^2 + (xtal.as2)^2)
	bstar = sqrt((xtal.bs0)^2 + (xtal.bs1)^2 + (xtal.bs2)^2)
	cstar = sqrt((xtal.cs0)^2 + (xtal.cs1)^2 + (xtal.cs2)^2)
	Variable hmax=floor(Qmax/astar), kmax=floor(Qmax/bstar), lmax=floor(Qmax/astar)

	String wNote = "waveClass=PowderLines;"
	wNote = ReplaceNumberByKey("Qmax",wNote,Qmax,"=")
	if (keV>0)
		wNote = ReplaceNumberByKey("keV",wNote,keV,"=")
	endif
	wNote = ReplaceStringByKey("xtalName",wNote,desc,"=")

	Variable NlinesMax=200
	Make/N=(NlinesMax,7)/D/FREE lines=NaN

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
				if (F2 < 0.001)
					continue
				endif
				if (m>=NlinesMax)
					NlinesMax += 200
					Redimension/N=(NlinesMax,-1) lines
				endif
				lines[m][0] = qhkl
				lines[m][1] = h
				lines[m][2] = k
				lines[m][3] = l
				lines[m][4] = theta*180/PI			// save angle in degrees
				lines[m][5] = F2
				lines[m][6] = F2 * Lorentz * Pol		// was F2*LP
				m += 1
			endfor
		endfor
	endfor
	NlinesMax = m
	Redimension/N=(NlinesMax,-1) lines				// trim off extras

	Make/N=(m)/FREE/D Qs, QsIndex, temp
	Qs = lines[p][0]									// Q values for every possible hkl
	MakeIndex Qs, QsIndex								//  to sort with increasing Q
	Variable i
	for (i=0;i<DimSize(lines,1);i+=1)
		temp = lines[p][i]								// for each column in lines
		IndexSort QsIndex, temp							// sort this column
		lines[][i] = temp[p]							// replace column in lines with sorted column
	endfor
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
	//	3	l
	//	4	theta (degree)
	//	5	| F| ^2
	//	6	Intensity
	//	7	multiplicity
	SetDimLabel 1,0,Q_nm,lineShort ;	SetDimLabel 1,1,h,lineShort ;		SetDimLabel 1,2,k,lineShort ;			SetDimLabel 1,3,l,lineShort
	SetDimLabel 1,4,theta,lineShort ;	SetDimLabel 1,5,F2,lineShort ;		SetDimLabel 1,6,Intensity,lineShort ;	SetDimLabel 1,7,mult,lineShort

	str = wname[0,31-7-1]+"_hklStr"
	Make/N=(NlinesMax)/T/O $str/WAVE=hklStr = ""
	Note/K hklStr, "waveClass=HKLlabelsPowderLines;LinesWave="+GetWavesDataFolder(lineShort,2)+";"
	wNote = ReplaceStringByKey("hklStrWave",wNote,str,"=")

	Make/N=(NlinesMax,3)/FREE hkls					// temp list of hkls for line, used to pick the "nicest" hkl
	Variable mult, Nlines=0
	for (Nlines=0,i=0; i<NlinesMax;)
		qhkl = lines[i][0]								// Q of this line
		mult = 0
		hkls = NaN
		do
			lineShort[Nlines][0,N_POWDER_LINES_COLUMNS-2] = lines[i][q]
			hkls[mult][] = lines[i][q+1]				// accumulate list of all hkls in this line
			mult += 1
			i += 1
		while(abs(lines[i][0]-qhkl)<1e-10)
		lineShort[Nlines][7] = mult
		Wave niceHKL = nicestHKL(hkls)
		lineShort[Nlines][1,3] = niceHKL[q-1]			// re-set hkl to the "nicest" one
		hklStr[Nlines] = nicestHKLstr(hkls)			// better for labeling, it includes multipliticites
		Nlines += 1
	endfor
	Redimension/N=(Nlines,-1) lineShort, hklStr		// trim to final size

	Note/K lineShort, wNote
	if (printIt)
		String plusMinus = SelectString(stringmatch(IgorInfo(2),"Macintosh"),"",num2char(-79))	// for Windows use ALT keycode 0177
		printf "Created Powder Lines in '%s',  with %d lines,  range of hkl is [%s%d, %s%d, %s%d]\r",wname,Nlines,plusMinus,hmax,plusMinus,kmax,plusMinus,lmax
	endif
	return lineShort
End
//
Static Function/WAVE nicestHKL(hklsIn)		// returns a free wave that is the nicest one from the list of hkl's
	// The nicest hkl has L>K>H, as much as possible
	Wave hklsIn
	if (DimSize(hklsIn,0)<1 || DimSize(hklsIn,1)!=3)
		return $""
	endif
	WaveStats/Q/M=1 hklsIn
	Variable N=V_npnts/3
	if (N<1)
		return $""
	endif
	Make/N=(N,3)/FREE hkls=hklsIn[p][q]

	ImageStats/M=1/G={0,N-1,2,2} hkls
	hkls = hkls[p][2]==V_max ? hkls[p][q] : NaN	// remove all hkl with L<V_max
	ImageStats/M=1/G={0,N-1,1,1} hkls
	hkls = hkls[p][1]==V_max ? hkls[p][q] : NaN	// remove all hkl with K<V_max
	ImageStats/M=1/G={0,N-1,0,0} hkls
	hkls = hkls[p][0]==V_max ? hkls[p][q] : NaN	// remove all hkl with H<V_max

	ImageStats/M=1/G={0,N-1,0,2} hkls		// find the best (probably only) one left
	Make/N=3/FREE hkl=hkls[V_maxRowLoc][p]
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
	if (DimSize(hklsIn,0)<1 || DimSize(hklsIn,1)!=3)
		return ""
	endif
	WaveStats/Q/M=1 hklsIn
	Variable N=V_npnts/3
	if (N<1)
		return ""
	endif
	Make/N=(N,3)/FREE hkls=hklsIn[p][q]

	ImageStats/M=1/G={0,N-1,2,2} hkls
	hkls = hkls[p][2]==V_max ? hkls[p][q] : NaN	// remove all hkl with L<V_max
	ImageStats/M=1/G={0,N-1,1,1} hkls
	hkls = hkls[p][1]==V_max ? hkls[p][q] : NaN	// remove all hkl with K<V_max
	ImageStats/M=1/G={0,N-1,0,0} hkls
	hkls = hkls[p][0]==V_max ? hkls[p][q] : NaN	// remove all hkl with H<V_max

	ImageStats/M=1/G={0,N-1,0,2} hkls		// find the best (probably only) one left
	Make/N=3/FREE hkl=hkls[V_maxRowLoc][p]
	String hklStr = "("+hkl2str(hkl[0],hkl[1],hkl[2])+")"

	// remove this hkl from hkls and do again
	hkls = hklsIn[p][q]

	MatrixOP/FREE checks = sumRows(abs(hkls))
	Redimension/N=(N) checks					// I am assuming that if the summs differ than the hkl is of different type
	Variable i,m, check0 = abs(hkl[0])+abs(hkl[1])+abs(hkl[2])
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

	String graph = StringFromList(0,FindGraphsWithWave(ww))
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
	elseif (StringMatch(WaveUnits(ww,0),"¡"))
		ix = 4
	endif

	if (strlen(WaveUnits(ww,-1))==0)
		bottomLabel = SelectString(ix,"Q  (\\Enm\\S-1\\M)","\\Zr150\\F'Symbol'q°\\F]0 \\E\\M")
	else
		bottomLabel = SelectString(ix,"Q  (\\U)","\\Zr150\\F'Symbol'q\\F]0\\U\\M")
	endif

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
	title += StringByKey("xtalName", wnote,"=")
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

	String graph = StringFromList(0,FindGraphsWithWave(ww))
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
	String leftLabel = SelectString(iy==6,"| F |\\S2\\M","Intensity")
	String bottomLabel = SelectString(ix,"Q  (nm\\S-1\\M)","\\Zr150\\F'Symbol'q°\\F]0\\M")

	String wnote = note(ww)
	Variable keV = NumberByKey("keV",wnote,"=")
	String title = "\\JR\\Zr125"
	title += StringByKey("xtalName", wnote,"=")
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

	if ((s.eventMod==3) && (s.eventCode ==3  || s.eventCode == 4))	// mouse down or moved with Shift key held down
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
			tagStr = "("+hkl2str(lines[iLine][1],lines[iLine][2],lines[iLine][3])+")"
		endif
		sprintf str, "\\Zr075\rQ = %g nm\\S-1\M\\Zr075",Qval
		tagStr += str
		if (theta>0)
			sprintf str, "\rtheta = %g¡",theta
			tagStr += str
		endif
		sprintf str, "\rmult = %g",mult
		tagStr += str

		GetWindow/Z $win  psize
		Variable X0 = Xval>AxisValFromPixel(win,"bottom",(V_left+V_right)/2) ? -3 : 3
		Variable Y0 = Yval>AxisValFromPixel(win,"left",(V_top+V_bottom)/2) ? -7 : 7
		anchor = SelectString(X0>0,"R","L")
		anchor += SelectString(Y0>0,"T","B")
		Tag/C/N=lineHKL/W=$win/A=$anchor/F=0/L=2/X=(X0)/Y=(Y0)/P=1 $NameOfWave(lines),iLine,tagStr
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

	String table = FindTableWithWave(w)
	if (strlen(table)>0)
		DoWindow/F $table
		return 0
	endif

	Wave/T hklStr = $(GetWavesDataFolder(w,1)+StringByKey("hklStrWave",note(w),"="))
	Variable right = WaveExists(hklStr) ? 709 : 588

	Edit/K=1/W=(5,44,right,434) w.ld
	ModifyTable width(Point)=1,title(w.d)="Powder Lines",width(w.l)=20, format(w.d)=3,width(w.d)=68
	if (WaveExists(hklStr))
		AppendToTable hklStr.y
		ModifyTable width(hklStr.y)=122, title(hklStr.d)="hkl's"
	endif
End


// MOVE this to JZT_Utility ???
//
Static Function/S FindTableWithWave(w)	// find the graph window which contains the specified wave
	Wave w
	if (!WaveExists(w))
		return ""
	endif
	String name0=GetWavesDataFolder(w,2), out=""
	String win,wlist = WinList("*",";","WIN:2"), clist, cwin
	Variable i,m,Nm=ItemsInList(wlist)
	for (m=0;m<Nm;m+=1)
		win = StringFromList(m,wlist)
		CheckDisplayed/W=$win w
		if (V_flag>0)
			out += win+";"
		else
			clist = ChildWindowList(win)
			for (i=0;i<ItemsInLIst(clist);i+=1)
				cwin = StringFromList(i,clist)
				CheckDisplayed/W=$(win+"#"+cwin) w
				if (V_flag>0)
					out += win+"#"+cwin+";"
					break
				endif
			endfor
		endif
	endfor
	return out
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
End
