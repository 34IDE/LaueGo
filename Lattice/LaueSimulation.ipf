#pragma rtGlobals=1		// Use modern global access method.
#pragma ModuleName=LaueSimulation
#pragma version = 1.16
#pragma IgorVersion = 6.11

#include  "microGeometryN", version>=1.85
#include  "LatticeSym", version>=5.14


Menu LaueGoMainMenuName
	SubMenu "Laue Simulation"
		"Make Simulated Laue Pattern", MakeSimulatedLauePattern(NaN,NaN)
		"  re-Plot a Laue Pattern",DisplaySimulatedLauePattern($"")
	End
End
Menu "Analysis"
	SubMenu "Laue Simulation"
		"Make Simulated Laue Pattern", MakeSimulatedLauePattern(NaN,NaN)
		"  re-Plot a Laue Pattern",DisplaySimulatedLauePattern($"")
	End
End



Static Constant hc = 1.239841857		// keV-nm

Function/WAVE MakeSimulatedLauePattern(Elo,Ehi,[h,k,l,recipSource,Nmax,detector,printIt])
	Variable Elo,Ehi					// energy range (keV)
	Variable h,k,l						// central hkl
	Wave recipSource					// OPTIONAL, either a (3,3) mat with recip lattice, OR a waveClass=IndexedPeakList* with a recip_lattice0
	Variable Nmax						// maximum number of reflection to generate
	Variable detector					// detector number [0,MAX_Ndetectors-1]
	Variable printIt
	Nmax = ParamIsDefault(Nmax) ? 100 : Nmax
	Nmax = (Nmax>0 && Nmax<=2000) ? Nmax : 100
	detector = ParamIsDefault(detector) ? 0 : detector
	detector = detector==round(limit(detector,0,MAX_Ndetectors-1)) ? detector : NaN
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt
	Variable h0=h, k0=k, l0=l								// rename central hkl

	STRUCT crystalStructure xtal							// crystal structure
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
	if (Ndetectors<1)											// nothing there
		return $""
	elseif (Ndetectors>1 && ParamIsDefault(detector))		// force a choice of which detector to use
		detector = NaN
	endif

	Variable startx=0, endx=(geo.d[0].Nx)-1, groupx=1	// used to define the ROI
	Variable starty=0, endy=(geo.d[0].Ny)-1, groupy=1

	String recipName=NameOfWave(recipSource)
	Variable err = (numtype(Elo+Ehi) || Elo<0 || Ehi<=Elo || numtype(detector))
	err = err || ( !WaveExists(recipSource) && ( numtype(h0+k0+l0) || (h0==0 && k0==0 && l0==0) ) )
	if (err)
		if (numtype(h0+k0+l0) || (h0==0 && k0==0 && l0==0))
			h0=0;	k0=0;	l0=1									// default to the (001) reflection
		endif
		Elo = (numtype(Elo) || Elo<0) ? 6 : Elo
		Ehi = (numtype(Ehi) || Ehi<=Elo) ? max(25,15+Elo) : Ehi
		Prompt Elo,"low energy cutoff (keV)"
		Prompt Ehi,"high energy cutoff (keV)"
		Prompt h0,"h"
		Prompt k0,"k"
		Prompt l0,"l"
		Prompt Nmax,"maximum number of reflections to accept"
//		String recipList = "{h0,k0,l0} at center;"+WaveListClass("IndexedPeakList*;IndexedPeakListSimulate*","*",""), str
		String recipList = "{h0,k0,l0} at center;"+WaveListClass("IndexedPeakList*","*",""), str
		recipList += WaveList("*",";","DIMS:2,MINROWS:3,MAXROWS:3,MINCOLS:3,MAXCOLS:3")
		Prompt recipName,"Orientation Source",popup,recipList
		DoPrompt "hkl & energy",h0,Elo,k0,Ehi,l0,recipName,Nmax
		if (V_flag)
			return $""
		endif
		Prompt startx,"X start of ROI [0, "+num2istr(geo.d[0].Nx -1)+"]"
		Prompt starty,"X start of ROI [0, "+num2istr(geo.d[0].Ny -1)+"]"
		Prompt endx,"X end of ROI [0, "+num2istr(geo.d[0].Nx -1)+"]"
		Prompt endy,"Y end of ROI [0, "+num2istr(geo.d[0].Ny -1)+"]"
		Prompt groupx,"binning in X direction"
		Prompt groupy,"binning in Y direction"
		String detectorName=StringFromList(detector,dList)
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
			else													// try to get recip from wave note
				Wave recip = decodeMatFromStr(StringByKey("recip_lattice0", note(recipSource),"="))	
			endif
		else
		endif
	endif
	recipName = SelectString(WaveExists(recipSource),"",NameOfWave(recipSource))
	if (printIt)
		printf "MakeSimulatedLauePattern(%g, %g",Elo,Ehi
		if (WaveExists(recipSource))
			printf ", recipSource=%s",recipName
		else
			printf ", h=%g,k=%g,l=%g", h0,k0,l0
		endif
		if (!ParamIsDefault(Nmax) && Nmax!=100)
			printf ", Nmax=%g", Nmax
		endif
		if (!ParamIsDefault(detector) && detector!=0)
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
	if (WaveExists(recip))									// calculate hkl at detector center
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
		Cross qhat, qcenter									// W_Cross is the new rotation axis
		Wave W_Cross=W_Cross
		Variable theta = atan2(norm(W_Cross),MatrixDot(qcenter,qhat))*180/PI
		Make/N=(3,3)/O/D/FREE rho
		microGeo#rotationMatFromAxis(W_Cross,theta,rho)
		KillWaves/Z W_Cross
		MatrixOp/O/FREE recip = rho x recip
	endif
	if (!WaveExists(recip))
		return $""
	endif

	String FullPeakIndexedName
	sprintf FullPeakIndexedName "SimulatedPeaks%d%d%d%s",abs(h0),abs(k0),abs(l0),detectorID2color(geo.d[detector].detectorID)
	Make/N=(Nmax,12,1)/O $FullPeakIndexedName
	Wave PeakIndexed = $FullPeakIndexedName

	//	find highest 2theta --> thetaMax
	Variable thetaMax
	thetaMax = pixel2q(geo.d[detector],0,0,$"")	// check all four corners of the detecgtor to find the maximum theta
	thetaMax = max(thetaMax,pixel2q(geo.d[detector],0,geo.d[detector].Ny -1,$""))
	thetaMax = max(thetaMax,pixel2q(geo.d[detector],geo.d[detector].Nx -1,0,$""))
	thetaMax = max(thetaMax,pixel2q(geo.d[detector],geo.d[detector].Nx -1,geo.d[detector].Ny -1,$""))
	Make/N=3/D/FREE ki={0,0,1}							//	this is a convention

	Variable dmin = hc/(2*Ehi*sin(thetaMax))
	Variable hmax = ceil(xtal.a / dmin)
	Variable kmax = ceil(xtal.b / dmin)
	Variable lmax = ceil(xtal.c / dmin)

	Variable/C pz
	Variable px,py, keV, Qlen, sintheta, Nspots, icnt=0
	String progressWin = ProgressPanelStart("",stop=1,showTime=1)	// display a progress bar
	for (l=0,Nspots=0; abs(l)<=lmax; l = l<0 ? -l : -(l+1))
		for (k=0; abs(k)<=kmax; k = k<0 ? -k : -(k+1))					// for kmax=4, k={0,-1,1,-2,2,-3,3,-4,4}
			for (h=0; abs(h)<=hmax; h = h<0 ? -h : -(h+1))
				hkl = {h,k,l}
				if (parallel_hkl_exists(h,k,l,Nspots,PeakIndexed))		// already got this reflection
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
				keV = Qlen*hc/(4*PI*sintheta)
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
				PeakIndexed[Nspots][7][0] = keV
				PeakIndexed[Nspots][9][0] = (px-(startx-FIRST_PIXEL)-(groupx-1)/2)/groupx		// change to binned pixels
				PeakIndexed[Nspots][10][0] = (py-(starty-FIRST_PIXEL)-(groupy-1)/2)/groupy	// pixels are still zero based
				PeakIndexed[Nspots][11][0] = detector	// detector number
				Nspots += 1
				if (Nspots+1 >= Nmax)
					h=Inf;	k=Inf;	l=inf					// force loops to end
				endif
			endfor
		endfor
		icnt += 1
		if (ProgressPanelUpdate(progressWin,icnt/(2*lmax+1)*100	))	// update progress bar
			break													//   and break out of loop
		endif
	endfor
	Variable executionTime = SecondsInProgressPanel(progressWin)
	DoWindow/K $progressWin
	Redimension/N=(Nspots,-1,-1) PeakIndexed		// trim to exact size
	PeakIndexed[][8][0] = 0								// error is always zero for a calculated spot
	PeakIndexed[][6][0] = 1								// set all intensities to 1
	if (printIt)
		printf "calculated %d simulated spots into the wave '%s',   took %s\r",Nspots,FullPeakIndexedName,Secs2Time(executionTime,5,0)
	endif

	String wnote=ReplaceStringByKey("waveClass","","IndexedPeakListSimulate","=")
	if (strlen(recipName)==0)
		sprintf str,"%g,%g,%g",h0,k0,l0
		wnote = ReplaceStringByKey("hkl",wnote,str,"=")
	else
		wnote = ReplaceStringByKey("recipSource",wnote,recipName,"=")
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
	sprintf str,"{%g,%g,%g}",hmax,kmax,lmax
	wnote = ReplaceStringByKey("hklmax",wnote,str,"=")
	wnote = ReplaceNumberByKey("executionTime",wnote,round(executionTime*1000)/1000,"=")
	sprintf str,"{%g, %g, %g, %g, %g, %g}",xtal.a,xtal.b,xtal.c,xtal.alpha,xtal.beta,xtal.gam
	wnote = ReplaceStringByKey("latticeParameters",wnote,str,"=")
	wnote = ReplaceStringByKey("lengthUnit",wnote,"nm","=")
	wnote = ReplaceNumberByKey("SpaceGroup",wnote,xtal.SpaceGroup,"=")
	wnote = ReplaceNumberByKey("SpaceGroupID",wnote,xtal.SpaceGroupID,"=")
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

	XYZ2pixel(geo.d[detector],ki,px,py)				// check for incident beam hitting detector after sample
	if (numtype(px+py))
		qhat = -ki
		XYZ2pixel(geo.d[detector],qhat,px,py)			// check for incident beam hitting detector before sample
	endif
	if ((px>=startx && px<=endx && py>=starty && py<=endy))	// incident beam hits detector, draw a cross at that point
		sprintf str,"{%g,%g}",px,py
		wnote = ReplaceStringByKey("pixel000",wnote,str,"=")
	endif
	Note/K PeakIndexed, wnote

	SetDimLabel 1,0,Qx,PeakIndexed		;	SetDimLabel 1,1,Qy,PeakIndexed
	SetDimLabel 1,2,Qz,PeakIndexed		;	SetDimLabel 1,3,h,PeakIndexed
	SetDimLabel 1,4,k,PeakIndexed			;	SetDimLabel 1,5,l,PeakIndexed
	SetDimLabel 1,6,Intensity,PeakIndexed	;	SetDimLabel 1,7,keV,PeakIndexed
	SetDimLabel 1,8,angleErr,PeakIndexed	;	SetDimLabel 1,9,pixelX,PeakIndexed
	SetDimLabel 1,10,pixelY,PeakIndexed	;	SetDimLabel 1,11,detNum,PeakIndexed	

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
Static Function parallel_hkl_exists(h,k,l,N,pw)
	Variable h,k,l
	Variable N				// number of rows in pw to check
	Wave pw

	N = min(N,DimSize(pw,0))
	Variable dot,i,tol=(1-1e-5)
	for (i=0;i<N;i+=1)
		dot = h*pw[i][3][0] + k*pw[i][4][0] + l*pw[i][5][0]
		dot /= sqrt(h*h + k*k + l*l)
		dot /= sqrt(pw[i][3][0]^2 + pw[i][4][0]^2 + pw[i][5][0]^2)
		if (dot>tol)
			return 1
		endif
	endfor
	return 0
End


Function DisplaySimulatedLauePattern(FullPeakIndexed)
	Wave FullPeakIndexed

	if (!WaveExists(FullPeakIndexed))
		String wName
		Prompt wName,"List of Indexed Peaks",popup,reverseList(WaveListClass("IndexedPeakList*","*",""))
		DoPrompt "indexed peaks",wName
		if (V_flag)
			return 1
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

	Display /W=(132,156,785,742) FullPeakIndexed[*][10][0] vs FullPeakIndexed[*][9][0]
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

	ModifyGraph/Z width={Aspect,(xhi-xlo)/(yhi-ylo)}
	ModifyGraph/Z grid=1, tick=2, mirror=1, minor=1, gridStyle=1
	ModifyGraph/Z lowTrip=0.001, standoff=0, axOffset(bottom)=-1
	ModifyGraph mode=3,marker=19,rgb=(1,4,52428),msize=2
	ModifyGraph grid=1,gridRGB=(45000,45000,65535)
	ModifyGraph axOffset(left) = (stringmatch(IgorInfo(2),"Macintosh" )) ? -1.1 : -2.1 		// mac, or pc

	Variable h,k,l,Elo,Ehi
	String str = StringByKey("hkl",wnote,"=")
	h = str2num(StringFromList(0,str,","))
	k = str2num(StringFromList(1,str,","))
	l = str2num(StringFromList(2,str,","))
	Elo=NumberByKey("Elo",wnote,"=")
	Ehi=NumberByKey("Ehi",wnote,"=")
	if (numtype(h+k+l)==0)
		str = StringByKey("structureDesc",wnote,"=")
		str += " ("+hkl2str(h,k,l)+")"
		str += "\rE=["+num2str(Elo)+", "+num2str(Ehi)+"]keV"
		TextBox/C/N=textTitle/F=0/B=1 str
	endif

	SetAxis bottom xlo-0.5,xhi+0.5
	SetAxis left yhi+0.5, ylo-0.5
	Label/Z left "y  \\U"
	Label/Z bottom "x  \\U"
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
	tagStr = "\\Zr090Indexed peak position\r" + str
	if (latticeSym#isValidSpaceGroupID(SpaceGroupID))
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



//=======================================================================================
// ================================= Start of Laue Simulation set panel =================================

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
	endif
	EnableDisableLaueSimControls(GetUserData("microPanel","","LaueSimPanelName"))
End

// ================================= End of Laue Simulation set panel ==================================
//=======================================================================================

