#pragma rtGlobals=1		// Use modern global access method.
#pragma ModuleName=LaueSimulation
#pragma version = 1.28
#pragma IgorVersion = 6.11

#include  "microGeometryN", version>=1.85
#include  "LatticeSym", version>=5.14

Static Constant hc_keVnm = 1.2398419739			// h*c (keV-nm)
Static Constant Nmax_DEFAULT = 150
Static Constant Elo_DEFAULT = 5
Static Constant Ehi_DEFAULT = 30


Menu LaueGoMainMenuName
	SubMenu "Laue Simulation"
		"Make Simulated Laue Pattern", MakeSimulatedLauePattern(NaN,NaN)
		"  re-Plot a Laue Pattern",DisplaySimulatedLauePattern($"")
		"Make wave with kf's...", Make_kf_Sim($"", printIt=1)
		"Make new Recip Lattice rotated about an existing...", MakeRotatedRecipLattice($"","",NaN, printIt=1)
	End
End
Menu "Analysis"
	SubMenu "Laue Simulation"
		"Make Simulated Laue Pattern", MakeSimulatedLauePattern(NaN,NaN)
		"  re-Plot a Laue Pattern",DisplaySimulatedLauePattern($"")
		"Make wave with kf's...", Make_kf_Sim($"", printIt=1)
		"Make new Recip Lattice rotated about an existing...", MakeRotatedRecipLattice($"","",NaN, printIt=1)
	End
End


// ==================================================================================================
// =============================== Start of Laue Simulation Make Sim  ===============================

Function/WAVE MakeSimulatedLauePattern(Elo,Ehi,[h,k,l,recipSource,Nmax,detector,printIt])
	Variable Elo,Ehi					// energy range (keV)
	Variable h,k,l						// central hkl
	Wave recipSource					// OPTIONAL, either a (3,3) mat with recip lattice, OR a waveClass=IndexedPeakList* with a recip_lattice0
	Variable Nmax						// maximum number of reflection to generate
	Variable detector					// detector number [0,MAX_Ndetectors-1]
	Variable printIt
	Nmax = ParamIsDefault(Nmax) || Nmax<=1 ? Nmax_DEFAULT : min(Nmax,2000)
	detector = ParamIsDefault(detector) ? 0 : detector
	detector = detector==round(limit(detector,0,MAX_Ndetectors-1)) ? detector : NaN
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt
	Variable h0=h, k0=k, l0=l							// rename central hkl

	STRUCT crystalStructure xtal						// crystal structure
	if (FillCrystalStructDefault(xtal))				//fill the geometry structure with current default values
		DoAlert 0,"Unable to load Crystal Structure"
		return $""
	endif
	STRUCT microGeometry geo
	if (FillGeometryStructDefault(geo,alert=1))	//fill the geometry structure with current default values
		return $""
	endif

	String dList = "", color
	Variable Ndetectors, i, idet
	for (i=0,Ndetectors=0; i<MAX_Ndetectors; i+=1)
		if (geo.d[i].used)
			Ndetectors += 1
			detector = numtype(detector) ? i : detector
			color = detectorID2color(geo.d[i].detectorID)
			dList += num2istr(i)+SelectString(strlen(color),";"," ("+color+");")
		endif
	endfor
	if (Ndetectors<1)												// nothing there
		return $""
	elseif (Ndetectors>1 && ParamIsDefault(detector))	// force a choice of which detector to use
		detector = NaN
	endif
	if (WaveExists(recipSource))
		Wave recip = recipSource
	endif

	Variable itemp = (detector>=0 && detector<MAX_Ndetectors) ? detector : 0
	Variable startx=0, endx=(geo.d[itemp].Nx)-1, groupx=1	// used to define the ROI
	Variable starty=0, endy=(geo.d[itemp].Ny)-1, groupy=1

	String recipName=NameOfWave(recipSource)
	Variable err = (numtype(Elo+Ehi) || Elo<0 || Ehi<=Elo || numtype(detector))
	err = err || ( !WaveExists(recipSource) && ( numtype(h0+k0+l0) || (h0==0 && k0==0 && l0==0) ) )
	if (err)
		if (numtype(h0+k0+l0) || (h0==0 && k0==0 && l0==0))
			h0=0;	k0=0;	l0=1										// default to the (001) reflection
		endif
		Elo = (numtype(Elo) || Elo<0) ? Elo_DEFAULT : Elo
		Ehi = (numtype(Ehi) || Ehi<=Elo) ? max(Ehi_DEFAULT,15+Elo) : Ehi
		Prompt Elo,"low energy cutoff (keV)"
		Prompt Ehi,"high energy cutoff (keV)"
		Prompt h0,"h"
		Prompt k0,"k"
		Prompt l0,"l"
		Prompt Nmax,"maximum number of reflections to accept"
		String recipList = "{h0,k0,l0} at center;"+WaveListClass("IndexedPeakList*","*",""), str
		recipList += WaveList("*",";","DIMS:2,MINROWS:3,MAXROWS:3,MINCOLS:3,MAXCOLS:3")
		recipList = RemoveDuplicatesFromList(recipList,matchCase=0)	// remove any duplicate names
		Prompt recipName,"Orientation Source",popup,recipList
		DoPrompt "hkl & energy",h0,Elo,k0,Ehi,l0,recipName,Nmax
		if (V_flag)
			return $""
		endif
		Prompt startx,"X start of ROI [0, "+num2istr(geo.d[itemp].Nx -1)+"]"
		Prompt starty,"X start of ROI [0, "+num2istr(geo.d[itemp].Ny -1)+"]"
		Prompt endx,"X end of ROI [0, "+num2istr(geo.d[itemp].Nx -1)+"]"
		Prompt endy,"Y end of ROI [0, "+num2istr(geo.d[itemp].Ny -1)+"]"
		Prompt groupx,"binning in X direction"
		Prompt groupy,"binning in Y direction"
		String detectorName=StringFromList(itemp,dList)
		Prompt detectorName,"Detector",popup,dList
		DoPrompt "ROI (full chip un-binned pixels)",startx,starty,endx,endy,groupx,groupy,detectorName
		if (V_flag)
			return $""
		endif
		detector = str2num(detectorName)
		endx = min(geo.d[detector].Nx-1,endx)
		endy = min(geo.d[detector].Ny-1,endy)

		Wave recipSource = $recipName
		if (WaveExists(recipSource))
			if (DimSize(recipSource,0)==3 && DimSize(recipSource,1)==3)
				Wave recip = recipSource					// recipSource is a 3x3 mat, all done
				Wave RL0 = recipFrom_xtal(xtal)		// returns a FREE wave with reciprocal lattice
				MatrixOP/FREE rho = recip x Inv(RL0)	// should be a rotaion matrix if recip is same as current xtal
				MatrixOP/FREE err0 = sum(abs(rho x rho^t - identity(3)))
				Variable ok=(err0[0]<1e-3 && MatrixDet(rho)>0.95)
				if (!ok)
					sprintf str, "Sum( (R x R^t) - Ident) = %g,  and |R| = %g\rContinue?\r",err0[0], MatrixDet(rho)
					DoAlert 1, str
					ok = V_flag==1
				endif
				if (!ok)
					sprintf str, "ERROR -- The given reicprocal lattice from 3x3 wave:\r    '%s'\ris not a proper rotation of the current Xtal", recipName
					print str
					DoAlert 0, str
					return $"" 
				endif
			else												// try to get recip from wave note
				Wave recip = decodeMatFromStr(StringByKey("recip_lattice0", note(recipSource),"="))	
			endif
		else
		endif
	endif
	recipName = SelectString(WaveExists(recipSource),"",NameOfWave(recipSource))
	if (printIt)
		printf "%sMakeSimulatedLauePattern(%g, %g",BULLET,Elo,Ehi
		if (WaveExists(recipSource))
			printf ", recipSource=%s",recipName
		else
			printf ", h=%g,k=%g,l=%g", h0,k0,l0
		endif
		if (Nmax != Nmax_DEFAULT)
			printf ", Nmax=%g", Nmax
		endif
		if (!ParamIsDefault(detector) || detector!=0)
			printf ", detector=%g", detector
		endif		
		printf ")\r"
	endif
	err = (numtype(Elo+Ehi) || Elo<0 || Ehi<=Elo || numtype(detector))
	err = err || ( !WaveExists(recipSource) && ( numtype(h0+k0+l0) || (h0==0 && k0==0 && l0==0) ) )
	if (err)
		return $""
	endif

	Make/N=3/D/FREE hkl={h0,k0,l0}, qcenter
	pixel2q(geo.d[detector],(geo.d[detector].Nx -1)/2,(geo.d[detector].Ny -1)/2, qcenter)
	if (WaveExists(recip))								// calculate hkl at detector center
		MatrixOP/FREE hkl = Inv(recip) x qcenter
		MatrixOP/FREE fff = maxVal(abs(hkl))
		hkl = round(hkl[p] * 24 / fff[0])
		h0 = hkl[0]
		k0 = hkl[1]
		l0 = hkl[2]
		lowestOrderHKL(h0,k0,l0)
		hkl = {h0, k0, l0}
	elseif (norm(hkl)>1e-5)								// make recip from {h0,k0,l0}
		// qcenter is Q vector that diffracts to the center of the detector, we want (h0,k0,l0) to be parallel to qcenter
		if (printIt)
			printf "vector from sample to detector center = %s\r",vec2str(qcenter,sep=", ")
		endif
		Wave recip = recipFrom_xtal(xtal)				// returns a FREE wave with reciprocal lattice
		// compute rotation that puts (h0,k0,l0) along qcenter (i.e. puts reference reflection on detector)
		MatrixOp/O/FREE qhat = recip x hkl
		normalize(qhat)
		Cross qhat, qcenter								// W_Cross is the new rotation axis
		Wave W_Cross=W_Cross
		Variable theta = atan2(norm(W_Cross),MatrixDot(qcenter,qhat))*180/PI
		Wave rho = rotationMatFromAxis(W_Cross,theta)
		KillWaves/Z W_Cross
		MatrixOp/O/FREE recip = rho x recip
	endif
	if (!WaveExists(recip))
		return $""
	endif

	String FullPeakIndexedName
	sprintf FullPeakIndexedName "SimulatedPeaks%d%d%d%s",abs(h0),abs(k0),abs(l0),detectorID2color(geo.d[detector].detectorID)
	Make/N=(Nmax,13,1)/O $FullPeakIndexedName
	Wave PeakIndexed = $FullPeakIndexedName
	SetScale d 0,0,"pixel", PeakIndexed

	//	find highest 2theta --> thetaMax
	Variable thetaMax
	thetaMax = pixel2q(geo.d[detector],0,0,$"")	// check all four corners of the detecgtor to find the maximum theta
	thetaMax = max(thetaMax,pixel2q(geo.d[detector],0,geo.d[detector].Ny -1,$""))
	thetaMax = max(thetaMax,pixel2q(geo.d[detector],geo.d[detector].Nx -1,0,$""))
	thetaMax = max(thetaMax,pixel2q(geo.d[detector],geo.d[detector].Nx -1,geo.d[detector].Ny -1,$""))
	Make/N=3/D/FREE ki={0,0,1}							//	this is a convention

	Wave hklRange = Find_hkl_range(xtal,geo.d[detector],recip, Elo,Ehi, startx,endx,starty,endy)
	Variable Nh,Nk,Nl
	Nh = (hklRange[1] - hklRange[0] + 1)
	Make/N=(Nh)/I/FREE hRange
	hRange = hklRange[0] + p
	Make/N=(Nh)/FREE size = hRange
	size += size<0 ? 0.2 : 0
	size = abs(size)
	Sort size, hRange

	Nk = (hklRange[3] - hklRange[2] + 1)
	Make/N=(Nk)/I/FREE kRange
	kRange = hklRange[2] + p
	Make/N=(Nk)/FREE size = kRange
	size += size<0 ? 0.2 : 0
	size = abs(size)
	Sort size, kRange

	Nl = (hklRange[5] - hklRange[4] + 1)
	Make/N=(Nl)/I/FREE lRange
	lRange = hklRange[4] + p
	Make/N=(Nl)/FREE size = lRange
	size += size<0 ? 0.2 : 0
	size = abs(size)
	Sort size, lRange
	WaveClear size

	Variable/C pz
	Variable px,py, keV, Qlen, sintheta, Nspots
	String progressWin = ProgressPanelStart("",stop=1,showTime=1)	// display a progress bar
	Variable ih,ik,il, breakAll=0
	for (il=0,Nspots=0; il<Nl && !breakAll; il+=1)
		l = lRange[il]
		for (ik=0; ik<Nk && !breakAll; ik+=1)
			k = kRange[ik]
			for (ih=0; ih<Nh && !breakAll; ih+=1)
				h = hRange[ih]
				hkl = {h,k,l}
				if (parallel_hkl_exists(hkl,Nspots,PeakIndexed))			// already got this reflection
					continue
				elseif (h==0 && k==0 && l==0)
					continue
				endif

				MatrixOp/O/FREE qhat = recip x hkl		// make the qhat to test
				Qlen = normalize(qhat)
				pz = q2pixel(geo.d[detector],qhat)
				px = real(pz)
				py = imag(pz)
				if (!(px>=startx && px<=endx && py>=starty && py<=endy))
					continue
				endif
				sintheta = -MatrixDot(ki,qhat)
				keV = Qlen*hc_keVnm/(4*PI*sintheta)
				if (keV!=limit(keV,Elo,Ehi) && keV<Inf)
					continue
				endif
				if (!allowedHKL(h,k,l,xtal))
					continue
				endif

				// passed all the tests, a good reflection, add it
				PeakIndexed[Nspots][0][0] = qhat[0]
				PeakIndexed[Nspots][1][0] = qhat[1]
				PeakIndexed[Nspots][2][0] = qhat[2]
				PeakIndexed[Nspots][3][0] = h
				PeakIndexed[Nspots][4][0] = k
				PeakIndexed[Nspots][5][0] = l
				PeakIndexed[Nspots][6][0] = intensityOfPeak(qhat,hkl,Elo,Ehi)
				PeakIndexed[Nspots][7][0] = keV
				PeakIndexed[Nspots][9][0] = (px-(startx-FIRST_PIXEL)-(groupx-1)/2)/groupx	// change to binned pixels
				PeakIndexed[Nspots][10][0] = (py-(starty-FIRST_PIXEL)-(groupy-1)/2)/groupy	// pixels are still zero based
				PeakIndexed[Nspots][11][0] = detector		// detector number
				Nspots += 1
				if (Nspots+1 >= Nmax)
					breakAll = 1			//h=Inf;	k=Inf;	l=inf						// force loops to end
				endif
			endfor
		endfor
		if (ProgressPanelUpdate(progressWin,il/Nl*100	))	// update progress bar
			break														//   and break out of loop
		endif
	endfor
	Variable executionTime = SecondsInProgressPanel(progressWin)
	DoWindow/K $progressWin
	Redimension/N=(Nspots,-1,-1) PeakIndexed		// trim to exact size
	PeakIndexed[][8][0] = 0								// error is always zero for a calculated spot
	PeakIndexed[][12][0] = (PeakIndexed[p][6][0])^0.3	// last column is just intensity^0,3, used for plotting
	if (printIt)
		printf "hkl ranges:  h=[%g,%g],  k=[%g,%g],  l=[%g,%g]\r", hklRange[0],hklRange[1], hklRange[2],hklRange[3], hklRange[4],hklRange[5]
		printf "calculated %d simulated spots into the wave '%s',   took %s\r",Nspots,FullPeakIndexedName,ElapsedTime2Str(executionTime)
	endif

	String wnote=ReplaceStringByKey("waveClass","","IndexedPeakListSimulate","=")
	if (strlen(recipName)==0)
		sprintf str,"%g,%g,%g",h0,k0,l0
		wnote = ReplaceStringByKey("hkl",wnote,str,"=")
	else
		wnote = ReplaceStringByKey("recipSource",wnote,GetWavesDataFolder(recipSource,2),"=")
	endif
	wnote = ReplaceNumberByKey("Elo",wnote,Elo,"=")
	wnote = ReplaceNumberByKey("Ehi",wnote,Ehi,"=")
	wnote = ReplaceNumberByKey("startx",wnote,startx,"=")
	wnote = ReplaceNumberByKey("endx",wnote,endx,"=")
	wnote = ReplaceNumberByKey("groupx",wnote,groupx,"=")
	wnote = ReplaceNumberByKey("starty",wnote,starty,"=")
	wnote = ReplaceNumberByKey("endy",wnote,endy,"=")
	wnote = ReplaceNumberByKey("groupy",wnote,groupy,"=")
	wnote = ReplaceNumberByKey("NpatternsFound",wnote,Nspots>0 ? 1 : 0,"=")
	wnote = ReplaceNumberByKey("Nindexed",wnote,Nspots,"=")
	wnote = ReplaceNumberByKey("thetaMax",wnote,thetaMax*180/PI,"=")
	sprintf str,"%g,%g",hklRange[0],hklRange[1]
	wnote = ReplaceStringByKey("hRange",wnote,str,"=")
	sprintf str,"%g,%g",hklRange[2],hklRange[3]
	wnote = ReplaceStringByKey("kRange",wnote,str,"=")
	sprintf str,"%g,%g",hklRange[4],hklRange[5]
	wnote = ReplaceStringByKey("lRange",wnote,str,"=")
	wnote = ReplaceNumberByKey("executionTime",wnote,round(executionTime*1000)/1000,"=")
	sprintf str,"{%g, %g, %g, %g, %g, %g}",xtal.a,xtal.b,xtal.c,xtal.alpha,xtal.beta,xtal.gam
	wnote = ReplaceStringByKey("latticeParameters",wnote,str,"=")
	wnote = ReplaceStringByKey("lengthUnit",wnote,"nm","=")
	wnote = ReplaceNumberByKey("SpaceGroup",wnote,xtal.SpaceGroup,"=")
	wnote = ReplaceStringByKey("SpaceGroupID",wnote,xtal.SpaceGroupID,"=")
	for (i=0;i<xtal.N;i+=1)
		sprintf str, "{%s %g %g %g %g}",xtal.atom[i].name,xtal.atom[i].x,xtal.atom[i].y,xtal.atom[i].z,xtal.atom[i].occ
		wnote = ReplaceStringByKey("AtomDesctiption"+num2istr(i+1),wnote,str,"=")
	endfor
	wnote = ReplaceStringByKey("structureDesc",wnote,xtal.desc,"=")
	sprintf str,"{{%g,%g,%g}{%g,%g,%g}{%g,%g,%g}}",recip[0][0],recip[1][0],recip[2][0],recip[0][1],recip[1][1],recip[2][1],recip[0][2],recip[1][2],recip[2][2]
	wnote = ReplaceStringByKey("recip_lattice0",wnote,str,"=")
	if (WaveExists(rho))
		sprintf str,"{{%g,%g,%g}{%g,%g,%g}{%g,%g,%g}}",rho[0][0],rho[1][0],rho[2][0],rho[0][1],rho[1][1],rho[2][1],rho[0][2],rho[1][2],rho[2][2]
		wnote = ReplaceStringByKey("rotation_matrix0",wnote,str,"=")
	endif

	XYZ2pixel(geo.d[detector],ki,px,py,onDetector=1)			// check for incident beam hitting detector after sample
	if (numtype(px+py))
		qhat = -ki
		XYZ2pixel(geo.d[detector],qhat,px,py, onDetector=1)	// check for incident beam hitting detector before sample
	endif
	if (numtype(px+py))
		pz = perpPixelOnDetector(geo.d[detector])				// pixel at perpendicualr direction to incident beam, and azimuth is multiple of 45¡
		px = real(pz)
		py = imag(pz)
	endif
	if (numtype(px+py)==0)											// incident beam hits detector, draw a cross at that point
		sprintf str,"{%g,%g}",px,py
		wnote = ReplaceStringByKey("pixel000",wnote,str,"=")
	endif
	Note/K PeakIndexed, wnote

	SetDimLabel 1,0,Qx,PeakIndexed				;	SetDimLabel 1,1,Qy,PeakIndexed
	SetDimLabel 1,2,Qz,PeakIndexed				;	SetDimLabel 1,3,h,PeakIndexed
	SetDimLabel 1,4,k,PeakIndexed				;	SetDimLabel 1,5,l,PeakIndexed
	SetDimLabel 1,6,Intensity,PeakIndexed	;	SetDimLabel 1,7,keV,PeakIndexed
	SetDimLabel 1,8,angleErr,PeakIndexed		;	SetDimLabel 1,9,pixelX,PeakIndexed
	SetDimLabel 1,10,pixelY,PeakIndexed		;	SetDimLabel 1,11,detNum,PeakIndexed	
	SetDimLabel 1,12,Intensity_03,PeakIndexed

	if (printIt && Nspots>0)
		String gName = StringFromLIst(0,WindowsWithWave(PeakIndexed,1))
		if (strlen(gName))
			DoWindow/F $gName
		else
			DisplaySimulatedLauePattern(PeakIndexed)
		endif
	endif
	return PeakIndexed
End
//
Static Function parallel_hkl_exists(hkl,N,pw)
	Wave hkl
	Variable N				// number of rows in pw to check
	Wave pw

	Variable tol=(1-1e-5), maxDot
	N = min(N,DimSize(pw,0))
	if (N<1)
		maxDot = -2
	elseif (N==1)
		Make/N=3/D/FREE hklsHave=pw[0][p+3]
		maxDot = MatrixDot(hkl,hklsHave) / (norm(hkl)*norm(hklsHave))
	else
		Make/N=(N,3)/D/FREE hklsHave=pw[p][q+3]
		MatrixOp/FREE maxDot0 = maxVal(sumRows(NormalizeRows(hklsHave) * NormalizeRows(RowRepeat(hkl,N))))
		maxDot = maxDot0[0]
	endif
	return (maxDot > tol)
End
//
//	Static Function parallel_hkl_exists(h,k,l,N,pw)
//		Variable h,k,l
//		Variable N				// number of rows in pw to check
//		Wave pw
//	
//		N = min(N,DimSize(pw,0))
//		Variable dot,i,tol=(1-1e-5)
//		for (i=0;i<N;i+=1)
//			dot = h*pw[i][3][0] + k*pw[i][4][0] + l*pw[i][5][0]
//			dot /= sqrt(h*h + k*k + l*l)
//			dot /= sqrt(pw[i][3][0]^2 + pw[i][4][0]^2 + pw[i][5][0]^2)
//			if (dot>tol)
//				return 1
//			endif
//		endfor
//		return 0
//	End
//
Static Function/WAVE Find_hkl_range(xtal,d,recip, Elo,Ehi, startx,endx,starty,endy)
	STRUCT crystalStructure &xtal					// crystal structure
	STRUCT detectorGeometry &d						// structure definition for a detector
	Wave recip
	Variable Elo,Ehi
	Variable startx,endx,starty,endy

	Variable midx=(startx+endx)/2, midy=(starty+endy)/2, Nc=9
	Make/N=(Nc,2)/D/FREE corner		// points on detector center, 4 corners, and midpoint of each edge
	corner[0][0]= {midx,  startx,midx  ,endx  ,endx,  endx,  midx,  startx,startx}
	corner[0][1]= {midy,  starty,starty,starty,midy,  endy,  endy,  endy  ,midy}

	Make/N=3/D/FREE qvec, qhat, ki={0,0,1}
	Make/N=3/D/FREE hklHi=-Inf, hklLo=Inf
	Variable i, Qlen, sintheta
	for (i=0;i<Nc;i+=1)
		pixel2q(d,corner[i][0],corner[i][1],qhat)
		sintheta = -MatrixDot(ki,qhat)

		Qlen = Elo/hc_keVnm * (4*PI*sintheta)
		qvec = Qlen * qhat
		MatrixOp/FREE hkl = Inv(recip) x qvec
		hklLo = min(hklLo,hkl)	;	hklHi = max(hklHi,hkl)

		Qlen = Ehi/hc_keVnm * (4*PI*sintheta)
		qvec = Qlen * qhat
		MatrixOp/FREE hkl = Inv(recip) x qvec
		hklLo = min(hklLo,hkl)	;	hklHi = max(hklHi,hkl)
	endfor
	hklLo = floor(hklLo)	;	hklHi = ceil(hklHi)

	Make/N=8/D/FREE hklRange={hklLo[0],hklHi[0], hklLo[1],hklHi[1], hklLo[2],hklHi[2]}
	return hklRange
End
//
Static Function intensityOfPeak(qhat,hklIN,Elo,Ehi)
	Wave qhat
	Wave hklIN
	Variable Elo, Ehi

	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))			//fill the lattice structure with test values
		DoAlert 0, "ERROR findClosestHKL()\rno lattice structure found"
		return 1
	endif

	Make/N=3/D/FREE hkl0=hklIN
	hkl0 = abs(hklIN)
	Variable hklMax=WaveMax(hkl0), iDiv
	Make/N=(hklMax,3)/D/FREE hkls=hklIN[q]/(p+1)
	MatrixOp/FREE hkls = abs(hkls - round(hkls))
	hkls = hkls > 1e-5 ? NaN : p
	MatrixOp/FREE flags = sumRows(hkls)
	flags = numtype(flags) ? NaN : p+1
	iDiv = WaveMax(flags)
	hkl0 = hklIN/iDiv										// hkl0 has all common divisors removed

	FUNCREF spectrumProto spectrumFunc = $"spectrum"
	Make/N=3/D/FREE kf, ki={0,0,1}
	kf = ki - 2*MatrixDot(ki,qhat)*qhat
	Variable keV, keV0
	keV0 = hc_keVnm / ( 2 * dSpacing(xtal,hkl0[0],hkl0[1],hkl0[2]) * sin(acos(MatrixDot(ki,kf))/2) )

	Variable Fhkl, mu, intens=0
	Variable istart = max(floor(Elo/keV0),1), i
	for (i=istart,keV=i*keV0; keV<Ehi; i+=1,keV+=keV0)	// for each harmonic
		if (keV>Elo)
			mu = LatticeSym#muOfXtal(xtal, keV)
			mu = numtype(mu) ? 1 : mu								// in case Cromer not loaded
			Fhkl = magSqr(Fstruct(xtal,i*hkl0[0],i*hkl0[1],i*hkl0[2], keV=keV))
			if (Fhkl > 1e-5)
				intens += Fhkl * spectrumFunc(keV) / mu
			endif
		endif
	endfor
	return 1e-4 * intens
End

Function spectrumProto(keV)
	Variable keV
	if (numtype(kev) || keV<=0 || keV>30)
		return 0
	endif
	Variable eV = keV * 1000
	return 1
End

// ================================ End of Laue Simulation Make Sim  ================================
// ==================================================================================================



// ==================================================================================================
// ============================= Start of Laue Simulation Display Sim  ==============================

Function DisplaySimulatedLauePattern(FullPeakIndexed)
	Wave FullPeakIndexed

	if (!WaveExists(FullPeakIndexed))
		String wName="", list=reverseList(WaveListClass("IndexedPeakList*","*","MINCOLS:11"))
		if (ItemsInlist(list)==1)
			wName = StringFromList(0,list)
		elseif (ItemsInlist(list)>1)
			Prompt wName,"List of Indexed Peaks", popup, list
			DoPrompt "indexed peaks",wName
			if (V_flag)
				return 1
			endif
		endif
		Wave FullPeakIndexed=$wName
	endif
	if (!WaveExists(FullPeakIndexed))
		return 1
	endif
	String win = StringFromList(0,WindowsWithWave(FullPeakIndexed,1))
	if (strlen(win))
		DoWindow/F $win
		return 0
	endif

	Variable dY=NumberByKey("endy",note(FullPeakIndexed),"=")*0.4
	Variable dX=NumberByKey("endx",note(FullPeakIndexed),"=")*0.4
	dY = numtype(dY) ? 586 : min(dY,586)				// set graph size depending on detector
	dX = numtype(dX) ? 586 : min(dX,586)
	Variable left=132, top=156, right=left+dX, bottom=top+dY	// was /W=(132,156,785,742)
	Display/K=1/W=(left,top,right,bottom) FullPeakIndexed[*][10][0] vs FullPeakIndexed[*][9][0]
	AppendToGraph FullPeakIndexed[*][10][0] vs FullPeakIndexed[*][9][0]
	GraphSimulateLaueStyle()
	if (exists("getSimulatedPeakInfoHook")==6)
		SetWindow kwTopWin,hook(peakInfo)=getSimulatedPeakInfoHook
	endif
	return 0
End
//
Function GraphSimulateLaueStyle()
	Wave ww = TraceNameToWaveRef("",StringFromList(0,TraceNameList("",";",1)))
	if (!WaveExists(ww))
		Abort "Unable to find any trace on top graph"
	endif

	String wnote = note(ww)
	Variable startx=NumberByKey("startx",wnote,"=")
	Variable starty=NumberByKey("starty",wnote,"=")
	Variable endx=NumberByKey("endx",wnote,"=")
	Variable endy=NumberByKey("endy",wnote,"=")
	Variable groupx=NumberByKey("groupx",wnote,"=")
	Variable groupy=NumberByKey("groupy",wnote,"=")
	Variable xlo = startx/groupx, xhi = endx/groupx
	Variable ylo = starty/groupx, yhi = endy/groupx
	String wname = StringFromList(0,TraceNameList("",";",1))

	ModifyGraph/Z width={Aspect,(xhi-xlo)/(yhi-ylo)}
	ModifyGraph/Z grid=1, tick=2, mirror=1, minor=1, gridStyle=1
	ModifyGraph/Z lowTrip=0.001, standoff=0, axOffset(bottom)=-1
	ModifyGraph mode=3, marker[0]=19, marker[1]=8, rgb=(0,0,0), msize=2
	Variable iColumnMarker = DimSize($wname,1)>12 ? 12 : 6
	ModifyGraph zmrkSize={$wname[*][iColumnMarker][0],*,*,2,12}
	ModifyGraph zColor[0]={$wname[*][7][0],*,*,Rainbow256}
	ModifyGraph grid=1,gridRGB=(45000,45000,65535)
	ModifyGraph axOffset(left) = (stringmatch(IgorInfo(2),"Macintosh" )) ? -1.1 : -2.1 		// mac, or pc

	Variable h,k,l,Elo,Ehi
	String title, str = StringByKey("hkl",wnote,"=")
	h = str2num(StringFromList(0,str,","))
	k = str2num(StringFromList(1,str,","))
	l = str2num(StringFromList(2,str,","))
	Elo=NumberByKey("Elo",wnote,"=")
	Ehi=NumberByKey("Ehi",wnote,"=")
	if (numtype(h+k+l)==0)
		title = StringByKey("structureDesc",wnote,"=")
		title += " ("+hkl2str(h,k,l)+")"
		title += "\rE=["+num2str(Elo)+","+num2str(Ehi)+"]keV"
		TextBox/C/N=textTitle/F=0/B=1/X=3/Y=2 title
	else
		title = StringByKey("structureDesc",wnote,"=")
		str = getHMsym2(SpaceGroupID2num(StringByKey("SpaceGroup",wnote,"=")))
		title += SelectString(strlen(str),"", "  "+str)
		title += "\r\\Zr077recip source = \""+StringByKey("recipSource",wnote,"=")+"\""
		sprintf str, "\r\\Zr120E = [%g, %g] keV", Elo,Ehi
		title += str

		Wave recipWave = $StringByKey("recipSource",wnote,"=")
		String hklStr = StringBykey("hkl",note(recipWave),"=")
		Variable angle = NumberBykey("angle",note(recipWave),"=")
		if (strlen(hklStr) && numtype(angle)==0 && DimSize(recipWave,0)==3 && DimSize(recipWave,1)==3)
			sprintf str, "\ra %g¡ rotation about the (%s)",angle,hklStr
			title += str
		endif
		if (IgorVersion()>=7)
			TextBox/C/N=title/F=0/G=(0,0,0,39321)/B=(65535,65535,65535,32768)/A=LT/X=4/Y=3 title
		else
			TextBox/C/N=title/F=0/A=LT/X=4/Y=3 title
		endif
	endif

	SetAxis bottom xlo-0.5,xhi+0.5
	SetAxis left yhi+0.5, ylo-0.5
	Label/Z left "y (\\U)"
	Label/Z bottom "x (\\U)"
	ShowInfo

	str = StringByKey("pixel000",wnote,"=")			// draw a cross at direct beam (if it intercepts the detector)
	if (strlen(str))
		Variable x0,y0
		sscanf str,"{%g,%g}",x0,y0
		if (V_flag==2 && numtype(x0+y0)==0)
			Variable dx = 0.08*abs(xhi-xlo), dy=0.08*(yhi-ylo)
			SetDrawLayer/K UserFront
			SetDrawEnv xcoord=bottom, ycoord=left, linethick=0.5 ,dash=1
			DrawLine (x0),(y0-dy),x0,(y0+dy)
			SetDrawEnv xcoord=bottom,ycoord=left, linethick=0.5 ,dash=1
			DrawLine (x0-dx),y0,(x0+dx),y0
		endif
	endif
End
//
Function getSimulatedPeakInfoHook(s)	// Command=fitted peak,  Shift=Indexed peak,  Command-Shift=strained peak,  &  Command-mouseDown=follow mouse
	STRUCT WMWinHookStruct &s

	// use shift to look at the indexing result, use command to look at fitted peaks
	String win = (s.winName)
	if (!(s.eventMod&2) || s.keycode)						// eventMod==2 is shift key, 8 is command, 1 is mouse, and no regular key held down
		Tag/K/N=indexedPeakInfo/W=$win
		return 0											// shift key not down, ignore
	endif

	Wave wTrace = TraceNameToWaveRef("",StringFromList(0,TraceNameList("",";",1)))
	if (!WaveExists(wTrace))
		return 1
	endif

	GetWindow $win psize
	Variable vert = limit((s.mouseLoc.v-V_top)/(V_bottom-V_top),0,1)	// fractional position on graph
	Variable horiz = limit((s.mouseLoc.h-V_left)/(V_right-V_left),0,1)

	Variable mx,my											// pixel position at the mouse
	GetAxis/W=$win/Q bottom
	if (V_flag)
		return 1
	endif
	Variable xlo=V_min, xhi=V_max
	mx = (xhi-xlo)*horiz + xlo

	GetAxis/W=$win/Q left
	if (V_flag)
		return 1
	endif
	Variable ylo=V_min, yhi=V_max
	my = (ylo-yhi)*vert + yhi

	String tagStr="", str, wnote=note(wTrace), SpaceGroupID
	Variable h,k,l,keV
	Variable dist2, m
	Variable px,py,i,N

	dist2=Inf
	m=-1													// find the indexed peak closest to the mouse-click
	N=DimSize(wTrace,0)
	for (i=0;i<N;i+=1)										// search through all indexed peaks
		px = wTrace[i][9][0]								// binned pixels
		py = wTrace[i][10][0]
		if ((px-mx)^2+(py-my)^2 < dist2)
			m = i
			dist2 = (px-mx)^2+(py-my)^2
		endif
	endfor
	if (m<0)
		return 1											// no indexed peak found
	endif
	h=wTrace[m][3][0]
	k=wTrace[m][4][0]
	l=wTrace[m][5][0]

	px = wTrace[m][9][0]
	py = wTrace[m][10][0]
	keV=wTrace[m][7][0]
	SpaceGroupID = StringByKey("SpaceGroupID",wnote,"=")
	if (strlen(SpaceGroupID)<1)
		SpaceGroupID = StringByKey("SpaceGroup",wnote,"=")
	endif
	sprintf str,"hkl=(%d %d %d),   %.4f keV\rpixel(%.2f, %.2f),   #%d",h,k,l,keV,px,py,m
	tagStr = "\\Zr090Calculated peak position\r" + str
	if (latticeSym#isValidSpaceGroupID(SpaceGroupID,3))
		sprintf str, "\r%s    Space Group %s", getHMsym2(SpaceGroupID2num(SpaceGroupID)),SpaceGroupID
		tagStr += str
	endif

	Variable theta = NaN									// find theta for this spot
	STRUCT microGeometry geo
	if (!FillGeometryStructDefault(geo))					//fill the geometry structure with current default values
		Make/N=3/O/D PeakInfoHook_qhatBL
		Wave qBL=PeakInfoHook_qhatBL
		Variable startx,groupx, starty,groupy				// ROI of the original image
		startx = NumberByKey("startx",wnote,"=")
		groupx = NumberByKey("groupx",wnote,"=")
		starty = NumberByKey("starty",wnote,"=")
		groupy = NumberByKey("groupy",wnote,"=")
		startx = numtype(startx) ? FIRST_PIXEL : startx
		groupx = numtype(groupx) ? 1 : groupx
		starty = numtype(starty) ? FIRST_PIXEL : starty
		groupy = numtype(groupy) ? 1 : groupy
		Variable pxUnb = (startx-FIRST_PIXEL) + groupx*px + (groupx-1)/2	// change to un-binned pixels
		Variable pyUnb = (starty-FIRST_PIXEL) + groupy*py + (groupy-1)/2	// pixels are still zero based
		if (numtype(pxUnb+pyUnb)==0)
			Variable dNum = max(detectorNumFromID(geo, StringByKey("detectorID", wnote,"=")),0)
			theta = pixel2q(geo.d[dNum],pxUnb,pyUnb,qBL)*180/PI				// get theta, and q^ in Beam Line system
		endif
		KillWaves/Z PeakInfoHook_qhatBL
	endif
	if (numtype(theta)==0)
		sprintf str "\r\F'Symbol'q\F]0 = %.3f¡",theta
		tagStr += str
	endif

	px = limit(px,xlo,xhi)									// needed in case (px,py) is outside the image
	py = limit(py,yhi,ylo)
	horiz = (px-xlo)/ (xhi-xlo)
	vert = (py-yhi)/(ylo-yhi)
	Variable x0 = horiz<0.5 ? 5 : -5,  y0 = vert<0.5 ? -5 : 5
	String anchor = SelectString(horiz<0.5,"R","L")+ SelectString(vert<0.5,"B","T")
	Tag/C/N=indexedPeakInfo/W=$win/A=$anchor/F=2/L=2/X=(x0)/Y=(y0)/P=1 $NameOfWave(wTrace),m,tagStr
	DoUpdate
	return 1												// 1 means do not send back to Igor for more processing
End
//
Static Function/C perpPixelOnDetector(d)
	// returns pixel on detector where a ray perpendicular to z-axis and with azimuth of n*45¡ hits detector
	// if no perpendicular ray hits the detector, return cmplx(NaN,NaN)
	STRUCT detectorGeometry &d				// geometry parameters for the detector

	Make/N=(8,3)/D/FREE xyz=0				// all 8 kf directions (every 45¡ and perp to z-axis)
	xyz[][0] = cos(p*PI/4)
	xyz[][1] = sin(p*PI/4)
	Wave pxpy = XYZ2pixelVEC(d,xyz, onDetector=1)	// returns pixel position (px,py) where each xyz vector will intercept detector
	Make/D/FREE pzC = {d.Nx/2, d.Ny/2}				// center of detector

	MatrixOP/FREE dist2 = sumRows(magSqr(pxpy - rowRepeat(pzC,8)))
	WaveStats/M=1/Q dist2									// find [px,py] that lies cloest to detector center
	Variable px=pxpy[V_minloc][0], py=pxpy[V_minloc][1]
	// print "closest azimuth is %g¡\r",V_minloc*45,"   ",px,"  ",py
	return cmplx(px,py)						// failed, closest pixel does not lie on detector
End

// ============================== End of Laue Simulation Display Sim  ===============================
// ==================================================================================================



// ==================================================================================================
// ================================ Start of Laue Simulation Utility ================================

Function/WAVE Make_kf_Sim(FullPeakIndexed, [printIt])
	// make wave with kf vectors from Qhats in FullPeakIndexed
	Wave FullPeakIndexed
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt
	if (!WaveExists(FullPeakIndexed))
		String wName="", list=reverseList(WaveListClass("IndexedPeakList*","*",""))
		if (ItemsInlist(list)==1)
			wName = StringFromList(0,list)
		elseif (ItemsInlist(list)>1)
			Prompt wName,"List of Indexed Peaks", popup, list
			DoPrompt "indexed peaks",wName
			if (V_flag)
				return $""
			endif
		endif
		Wave FullPeakIndexed=$wName
		printIt = 1
	endif
	if (!WaveExists(FullPeakIndexed))
		return $""
	endif
	Variable N=DimSize(FullPeakIndexed,0)
	if (N<1)
		return $""
	endif
	if (printIt)
		printf "¥Make_kf_Sim(%s)\r",NameOfWave(FullPeakIndexed)
	endif

	Make/N=3/D/FREE ki={0,0,1}
	Make/N=(N,3)/D/FREE qhat=FullPeakIndexed[p][q]
	MatrixOp/FREE qhat = NormalizeRows(qhat)

	// kf^ = ki^ - 2*(ki^ . q^)*q^		note: (ki^ . q^) must be NEGATIVE (both Bragg & Laue)
	String name = NameOfWave(FullPeakIndexed)+"_kf"
	name = ReplaceString("__",CleanupName(name,0),"_")
	Make/N=(N,3)/D/O $name/WAVE=kf
	MatrixOP/FREE dots = sumRows(RowRepeat(ki,N)*qhat)
	kf = ki[q] - 2*dots[p]*qhat[p][q]

	String wNote=note(FullPeakIndexed)
	wNote = ReplaceStringByKey("waveClass",wNote,"kfSimulate","=")
	wNote = ReplaceStringByKey("source",wNote,NameOfWave(FullPeakIndexed),"=")
	Note/K kf, wNote
	if (printIt)
		printf "\t\tSaving kf wave in: \"%s\"\r",name
	endif
	return kf
End


Function/WAVE MakeRotatedRecipLattice(recipSource,hklStr,angle,[printIt])
	Wave recipSource		// either a class="IndexedPeakList" or a 3x3 matrix
	String hklStr
	Variable angle
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetrtStackInfo(2))==0 : printIt

	String recipName=NameOfWave(recipSource)
	Wave hkl = minStr2Vec(hklStr,3)
	Variable dim=3

	if (!WaveExists(recipSource) || numpnts(hkl)<3 || numtype(angle) || abs(angle)>360)
		angle = numtype(angle) ? 60 : angle
		hklStr = SelectString(strlen(hklStr), "111", hklStr)
		String recipList = WaveListClass("IndexedPeakList*,StrainPeakList","*","")
		recipList += WaveList("*",";","DIMS:2,MINROWS:3,MAXROWS:3,MINCOLS:3,MAXCOLS:3")
		recipList += WaveListClass("recipLattice*","*","DIMS:2,MINROWS:2,MAXROWS:3,MINCOLS:2,MAXCOLS:3")
		// recipList += WaveListClass("recipLattice*","*","DIMS:2,MINROWS:3,MAXROWS:3,MINCOLS:3,MAXCOLS:3")
		recipList += ";-standard orientation-"
		recipList = RemoveFromList("epsilon;epsilonBL;epsilonSample;epsilonXHF", recipList)
		Prompt recipName,"Orientation Source",popup,recipList
		Prompt hklStr "hkl of axis"
		Prompt angle, "rotation angle about hkl (degree)"
		DoPrompt "hkl & energy",recipName, hklStr, angle
		if (V_flag)
			return $""
		endif
		if (CmpStr(recipName, "-standard orientation-")==0)
			STRUCT crystalStructure xtal			// temporary crystal structure
			FillCrystalStructDefault(xtal)			//fill the lattice structure with default values
			dim = xtal.dim
			Wave recip0 = recipFrom_xtal(xtal)	// returns a FREE wave with reciprocal lattice
			Wave recipSource = recip0
		else
			Wave recipSource = $recipName
		endif
		Wave hkl = minStr2Vec(hklStr,3)
	endif
	if (printIt)
		if (CmpStr(recipName, "-standard orientation-")==0)
			printf "%sMakeRotatedRecipLattice($\"standard orientation\", \"%s\", %g)\r", BULLET,hklStr,angle
		else
			printf "%sMakeRotatedRecipLattice(%s, \"%s\", %g)\r", BULLET,recipName,hklStr,angle
		endif
	endif

	if (!WaveExists(recipSource) || numpnts(hkl)<3 || numtype(angle) || abs(angle)>360)
		return $""
	endif


	String wnote=note(recipSource), recipStr="", key
	recipStr = StringByKey("recip_lattice_refined",wnote,"=")
	key = "recip_lattice_refined"
	if (strlen(recipStr)==0)
		recipStr = StringByKey("recip_lattice",wnote,"=")
		key = "recip_lattice"
	endif
	if (strlen(recipStr)==0)
		recipStr = StringByKey("recip_lattice0",wnote,"=")
		key = "recip_lattice0"
	endif

	Make/N=(3,3)/D/FREE rl=NaN
	if (strlen(recipStr))						// set rl from a recipStr
		Wave rl = str2recip(recipStr)
	endif
	if (strlen(recipStr)==0 && WaveDims(recipSource)==2 && DimSize(recipSource,0)==3 && DimSize(recipSource,1)==3)
		rl = recipSource[p][q]				// set rl from a 3x3 matrix
	endif

	// now rotate rl.  Need to rotate about the q(hkl), not hkl itself
	MatrixOP/FREE axis = normalize(rl x hkl)
	Wave rho = rotationMatFromAxis(axis,angle)
	MatrixOP/FREE rl = rho x rl

	String wname = SelectString(CmpStr(recipName, "-standard orientation-"),"standard", recipName)
	wname = ReplaceString("orientation",wname,"")
	wname = ReplaceString(" ",wname,"")
	String ending
	sprintf ending, "_%s_%d", hklStr,angle
	ending = ReplaceString(" ",ending,"")
	ending = ReplaceString("__",ending,"_")
	wname = AddEndingToWaveName(wname,ending)
	Make/N=(3,3)/D/O $wname/WAVE=recip = rl
	SetDimLabel 1,0,$"a*",recip
	SetDimLabel 1,1,$"b*",recip
	SetDimLabel 1,2,$"c*",recip
	wnote = "waveClass=recipLattice"
	wnote = ReplaceStringByKey("hkl",wnote,hklStr,"=")
	wnote = ReplaceNumberByKey("angle",wnote,angle,"=")
	Note/K recip, wnote

	if (printIt)
		if (strlen(recipStr))
			printf "Rotated key='%s' from waveNote of '%s' around the (%s) by %g¡\r",key,recipName,hklStr,angle
		else
			printf "Rotated reciprocal lattice in '%s' around the (%s) by %g¡\r",recipName,hklStr,angle
		endif
		String str
		sprintf str, "New rotated reciprocal lattice:  \"%s\"\r    a*          b*          c*   ",wname
		printWave(rl,name=str,brief=1,zeroThresh=1e-8)
	endif
	return recip
End

// ==================================================================================================
// ================================= End of Laue Simulation Utility =================================



// ==================================================================================================
// ================================= Start of Laue Simulation Panel =================================

Function/T FillLaueSimParametersPanel(strStruct,hostWin,left,top)
	String strStruct									// optional passed value of xtal structure, this is used if passed
	String hostWin										// name of home window
	Variable left, top									// offsets from the left and top

	SetWindow kwTopWin,userdata(LaueSimPanelName)=hostWin+"#LaueSimPanel"
	NewPanel/K=1/W=(left,top,left+221,top+445+30)/HOST=$hostWin
	ModifyPanel frameStyle=0, frameInset=0
	RenameWindow #,LaueSimPanel

	Button buttonMakeLaueSim,pos={29,5},size={160,20},proc=LaueSimulation#LaueSimButtonProc,title="Make a Laue Simulation"
	Button buttonMakeLaueSim,help={"Make a Laue Simulation"}

	Button buttonLaueSimRePlot,pos={29,40},size={160,20},proc=LaueSimulation#LaueSimButtonProc,title="Re-Plot a Simulation"
	Button buttonLaueSimRePlot,help={"Re-Plot a Laue Simulation"}

	Button buttonLaueSim_kf,pos={29,75},size={160,20},proc=LaueSimulation#LaueSimButtonProc,title="Make kf Wave"
	Button buttonLaueSim_kf,help={"Make wave with kf vectors from Qhats in FullPeakIndexed"}

	Button buttonLaueSim_Rotate pos={29,110},size={160,20},proc=LaueSimulation#LaueSimButtonProc,title="Make New Rotated R.L."
	Button buttonLaueSim_Rotate,help={"Make a new Reciprocal Lattice rotated about an hkl of an existing recip lattice"}

	PopupMenu popupLaueSimTables,pos={29,145},size={66,20},proc=LaueSimulation#TablesLaueSimPopMenuProc,title="Display Tables..."
	PopupMenu popupLaueSimTables,help={"Show table of indexed peak positions and any reicprocal lattices"}
	PopupMenu popupLaueSimTables,fSize=14,mode=0,value= #"\"Simulated Peaks;Indexed Peaks;Strain Peaks;Reciprocal Lattices\""

	EnableDisableLaueSimControls(hostWin+"#LaueSimPanel")
	return "#LaueSimPanel"
End
//
Static Function EnableDisableLaueSimControls(win)				// here to enable/disable
	String win												// window (or sub window) to use
	Variable d

	Button buttonMakeLaueSim,win=$win,disable=0
	d = strlen(WaveListClass("IndexedPeakListSimulate*","*",""))<1 ? 2 : 0
	Button buttonLaueSimRePlot,win=$win,disable=d
	Button buttonLaueSim_kf,win=$win,disable=d
	d = strlen(WaveListClass("FittedPeakList*;IndexedPeakListSimulate;IndexedPeakList,StrainPeakList*,recipLattice*","*","")) ? 0 : 2
	PopupMenu popupLaueSimTables,win=$win,disable=0
End
//
Static Function LaueSimButtonProc(B_Struct) : ButtonControl
	STRUCT WMButtonAction &B_Struct
	if (B_Struct.eventCode != 2)
		return 0
	endif
	String ctrlName=B_Struct.ctrlName

	if (stringmatch(ctrlName,"buttonMakeLaueSim"))
		MakeSimulatedLauePattern(NaN,NaN,printIt=1)
	elseif (stringmatch(ctrlName,"buttonLaueSimRePlot"))
		DisplaySimulatedLauePattern($"")
	elseif (stringmatch(ctrlName,"buttonLaueSim_kf"))
		Make_kf_Sim($"", printIt=1)
	elseif (stringmatch(ctrlName,"buttonLaueSim_Rotate"))
		MakeRotatedRecipLattice($"","",NaN, printIt=1)
	endif
	EnableDisableLaueSimControls(GetUserData("microPanel","","LaueSimPanelName"))
End

Static Function TablesLaueSimPopMenuProc(ctrlName,popNum,popStr) : PopupMenuControl
	String ctrlName
	Variable popNum
	String popStr
	if (stringmatch(popStr,"Simulated Peaks"))
		DisplayTableOfWave($"",classes="IndexedPeakListSimulate",promptStr="Wave with list of Simulated peaks")
	elseif (stringmatch(popStr,"Indexed Peaks"))
		DisplayTableOfWave($"",classes="IndexedPeakList",promptStr="Wave with list of indexed peaks")
	elseif (stringmatch(popStr,"Strain Peaks"))
		DisplayTableOfWave($"",classes="StrainPeakList*",promptStr="Wave with result of Strain Refinement")
	elseif (stringmatch(popStr,"Reciprocal Lattices"))
		DisplayTableOfWave($"",classes="recipLattice*",promptStr="Wave with a reciprocal lattice",options="DIMS:2,MINCOLS:2,MAXCOLS:3,MINROWS:2,MAXROWS:3")
	endif
End

// ================================== End of Laue Simulation Panel ==================================
// ==================================================================================================
