#pragma rtGlobals=1		// Use modern global access method.
#pragma ModuleName=StereoGraphicProjection
#pragma version = 2.85
#include "LatticeSym", version>=4.13


Menu "Analysis"
	"Stereographic Projection", MakeStereo(NaN,NaN,NaN,NaN,NaN,NaN)
End


Function MakeStereo(Hz,Kz,Lz,hklmax,fntSize,phi,[Qmax,hklPerp, WulffStepIn,WulffPhiIn,markers,outlines])
	Variable Hz,Kz,Lz			// hkl of the pole
	Variable hklmax
	Variable fntSize
	Variable phi				// aximuthal angle (degree)
	Variable Qmax				// maximum Q (1/nm)
	Wave hklPerp				// optional hkl giving the perpendicular direction long=0
	Variable WulffStepIn, WulffPhiIn	// for Wulff Net, use WulffStep=0 for no Wulff Net
	Variable markers			// flag, show markers for each hkl
	String outlines			// string with outlines of detectors to display
	Qmax = ParamIsDefault(Qmax) ? Inf : Qmax
	Qmax = Qmax > 0 ? Qmax : Inf
	markers = ParamIsDefault(markers) ? 0 : markers
	markers = numtype(markers) ? 0 : !(!markers)	// if no Qmax, then Q are limited by only hkl
	outlines = SelectString(ParamIsDefault(outlines),outlines,"")
	// format of outlines is detctorOutline0;detctorOutline1;detctorOutline2;
	// detctorOutline0 = "r=65535:g=43688:b=32768:hkl0=-0.2 0.5 -0.8:hkl1=-0.3 0.7 -0.5:hkl2=0.3 0.8 -0.5:hkl3=0.2 0.5 -0.8:;
	// the r,g,b give the color to use, and each of the hkli are a normalized vector in hkl units

	if (!(fntSize>0) || !(hklmax>0) || numtype(Hz+Kz+Lz) || numtype(phi))
		Hz = numtype(Hz) ? 3 : Hz
		Kz = numtype(Kz) ? -1 : Kz
		Lz = numtype(Lz) ? 7 : Lz
		fntSize = fntSize>0 ? fntSize : 9
		hklmax = hklmax>0 ? hklmax : 6
		phi = numtype(phi) ? 0 : phi
		String hklStr
		sprintf hklStr,"%g %g %g",Hz,Kz,Lz
		Prompt hklStr,"(hkl) of the pole"
		Prompt fntSize,"size of font for hkl"
		Prompt hklmax,"maximum hkl to accept"
		Prompt Qmax,"maximum Q to show (a redundant limit) (1/nm)"
		Prompt phi,"rotation angle counter-clockwise about pole (degree)"
		Prompt markers,"Also show 'x' at hkl's",popup,"No Markers;Show Markers"
		markers = markers+1
		DoPrompt "hkl max & font size",hklStr,hklmax,Qmax,fntSize,phi, markers
		if (V_flag)
			return 1
		endif
		markers = markers==2
		sscanf hklStr, "%g %g %g",Hz,Kz,Lz
		printf "MakeStereo(%g,%g,%g,  %g,%g, %g",Hz,Kz,Lz,hklmax,fntSize,phi
		if (!(numtype(Qmax)==1 && Qmax>0))
			printf ", Qmax=%g",Qmax
		endif
		if (!ParamIsDefault(hklPerp))
			printf ", hklPerp={%g,%g,%g}",hklPerp[0],hklPerp[1],hklPerp[2]
		endif
		if (markers)
			printf ", markers=%g",markers
		endif
		if (!ParamIsDefault(WulffStepIn))
			printf ", WulffStepIn=%g",WulffStepIn
		endif
		if (!ParamIsDefault(WulffPhiIn))
			printf ", WulffPhiIn=%g",WulffPhiIn
		endif
		printf ")"
		if (strlen(outlines))
			printf "\t\t// also given data to display detector outlines"
		endif
		printf "\r"
	endif
	if (!(fntSize>0) || !(hklmax>0) || numtype(Hz+Kz+Lz))
		return 1
	endif

	if (!StringMatch(GetRTStackInfo(2),"StereoOfIndexedPattern"))
		if (setLattice())
			return 1
		endif
	endif
	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))		//fill the lattice structure with current values
		DoAlert 0, "No Lattice, please set one"
		MakeLatticeParametersPanel("")
		return 1
	endif
	if (LatticeSym#LatticeBad(xtal))
		DoAlert 0, "Lattice is INvalid, re-set it"
		MakeLatticeParametersPanel("")
		return 1
	endif

	Variable maxWaveLen=9000
	Make/N=(maxWaveLen)/FREE longs=NaN, lats=NaN
	Make/N=(maxWaveLen,3)/W/FREE hkls
	Make/N=3/D/FREE qhat
	Variable qmag										// magnitude of qvector

	Variable hmax, kmax, lmax						// get separate limits on h,k, & l based on hklmax and Qmax
	qmag = sqrt(xtal.as0^2 + xtal.as1^2 + xtal.as2^2)
	hmax = min(ceil(Qmax/qmag),hklmax)
	qmag = sqrt(xtal.bs0^2 + xtal.bs1^2 + xtal.bs2^2)
	kmax = min(ceil(Qmax/qmag),hklmax)
	qmag = sqrt(xtal.cs0^2 + xtal.cs1^2 + xtal.cs2^2)
	lmax = min(ceil(Qmax/qmag),hklmax)
	Variable H,K,L
	Variable N=0, lat, long
	String str
	for (L=0;L<=lmax;L=signedAdvance(L))
		for (K=0;K<=kmax;K=signedAdvance(K))
			for (H=0;H<=hmax;H=signedAdvance(H))
				if (H==0 && K==0 && L==0)			// skip the (000)
					continue
				endif
				if (magsqr(Fstruct(xtal,h,k,l)) < 0.001)
					continue
				endif
				qhat[0] = H*xtal.as0+K*xtal.bs0+L*xtal.cs0		// qmag*qhat =  recip x hkl
				qhat[1] = H*xtal.as1+K*xtal.bs1+L*xtal.cs1
				qhat[2] = H*xtal.as2+K*xtal.bs2+L*xtal.cs2
				qmag = normalize(qhat)
				if (qmag > Qmax)
					continue
				endif
				if (ParamIsDefault(hklPerp))
					hkl2LongLat(Hz,Kz,Lz,qhat,long,lat,xtal)
				else
					hkl2LongLat(Hz,Kz,Lz,qhat,long,lat,xtal,hklPerp=hklPerp)
				endif
				if (lat>=-0.01 && N<maxWaveLen)
					longs[N] = long
					lats[N] = lat
					hkls[N][0]=H  ;  hkls[N][1]=K  ;  hkls[N][2]=L
					N += 1
				endif
			endfor
		endfor
	endfor
	Redimension/N=(N) longs,lats
	lats = abs(lats[p]) < 1e-12 ? 0 : lats[p]	// This is needed because of a bug in Project
	longs += phi										// add on the azimuthal rotation
	Project/M=1 longs, lats
	Wave W_XProjection=W_XProjection, W_YProjection=W_YProjection
	Duplicate/FREE W_XProjection, Xspots
	Duplicate/FREE W_YProjection, Yspots
	KillWaves/Z W_XProjection, W_YProjection

	// remove duplicates
	Variable i, j
	for (i=N-1;i>=0;i-=1)
		for (j=0;j<i;j+=1)
			if (abs(Xspots[i]-Xspots[j])+abs(Yspots[i]-Yspots[j]) < 0.01)
				Xspots[i] = NaN
				Yspots[i] = NaN
				break
			endif
		endfor
	endfor

	Display/W=(358.5,41,848.25,518.75)/K=1	// draw the graph
	NewFreeAxis/O/L emptyLeft
	NewFreeAxis/O/B emptyBottom
	ModifyGraph gfMult=60,width={Aspect,1},tick=3,noLabel=2,standoff=0,axThick=0,freePos=0,freePos=0
	SetAxis emptyLeft -2,2
	SetAxis emptyBottom -2,2
//	sprintf str, "\Zr142#%d   %s\r\\Zr141(%d %d %d)\r90\\F'Symbol'°\\F]0",xtal.SpaceGroup,getSymString(xtal.SpaceGroup),Hz,Kz,Lz
//	sprintf str, "\Zr142#%d   %s\r\\Zr141(%s)\r90\\F'Symbol'°\\F]0",xtal.SpaceGroup,getSymString(xtal.SpaceGroup),hkl2str(Hz,Kz,Lz)
	if (strlen(xtal.desc))
		sprintf str, "\\JR\\Zr200%s\r\\Zr071#%d   %s\r\\Zr141(%s)\r90\\F'Symbol'°\\F]0",xtal.desc,xtal.SpaceGroup,getSymString(xtal.SpaceGroup),hkl2str(Hz,Kz,Lz)
	else
		sprintf str, "\\JR\\Zr142#%d   %s\r\\Zr141(%s)\r90\\F'Symbol'°\\F]0",xtal.SpaceGroup,getSymString(xtal.SpaceGroup),hkl2str(Hz,Kz,Lz)
	endif
	TextBox/N=pole/F=0/A=RT/X=0/Y=0/B=1 str
	Button zoomIn,pos={5,5},size={30,20},proc=ZoomStereoButtonProc,title="+"
	Button zoomOut,pos={5,25},size={30,20},proc=ZoomStereoButtonProc,title="Ñ"
	Variable/G WulffNetStep=NumVarOrDefault("WulffNetStep",0), WulffNetPhi=NumVarOrDefault("WulffNetPhi",0)
	WulffNetStep = ParamIsDefault(WulffStepIn) ? WulffNetStep : WulffStepIn	// use passed values if given
	WulffNetPhi = ParamIsDefault(WulffPhiIn) ? WulffNetPhi : WulffPhiIn
	SetVariable setvarWulffStep,pos={47,7},size={76,18},proc=StereoGraphicProjection#UpdateWulffNetProc,title="step"
	SetVariable setvarWulffStep,fSize=12,limits={0,90,5},value= WulffNetStep
	SetVariable setvarWulffPhi,pos={64,30},size={57,17},proc=StereoGraphicProjection#UpdateWulffNetProc,title="f"
	SetVariable setvarWulffPhi,font="Symbol",fSize=14
	SetVariable setvarWulffPhi,limits={-90,90,1},value= WulffNetPhi
	SetVariable setvarWulffPhi,disable=(WulffNetStep ? 0 : 1)

	String list = ReplaceNumberByKey("stereoRadius","",90,"=")
	list = ReplaceNumberByKey("Hz",list,Hz,"=")
	list = ReplaceNumberByKey("Kz",list,Kz,"=")
	list = ReplaceNumberByKey("Lz",list,Lz,"=")
	SetWindow kwTopWin, userdata(stereoRadius)=list

	SetDrawLayer/K UserFront
	Wave Dpolys = MakeDetectorPolyWaves(outlines,Hz,Kz,Lz,phi,hklPerp,xtal)
	if (WaveExists(Dpolys))
		Variable Nd=DimSize(Dpolys,0), Np=DimSize(Dpolys,1)
		Make/N=(Np)/D/FREE xy
		for (i=0;i<Nd;i+=1)							// draw each detector polygon
			xy = Dpolys[i][p]
			if (numtype(sum(xy))==0)				// valid polygon vertices
				Wave rgb = str2vec(StringByKey("rgb"+num2istr(i),note(Dpolys),"="))
				SetDrawEnv xcoord=emptyBottom,ycoord=emptyLeft,linefgc=(rgb[0],rgb[1],rgb[2]),linethick=2,fillpat=0
				DrawPoly/ABS 0,0, 1,1, {xy[0],xy[1],xy[2],xy[3],xy[4],xy[5],xy[6],xy[7],xy[8],xy[9]}
			endif
		endfor
	endif

	if (1)												// Draw round marker at center
		SetDrawEnv xcoord= emptyBottom,ycoord= emptyLeft,fillpat= 0
		SetDrawEnv linefgc= (43520,43520,43520)
		DrawOval -0.01,-0.01,0.01,0.01
	else
		SetDrawEnv xcoord= emptyBottom,ycoord= prel
		DrawLine 0,0.49,0,0.51						// draw a small cross at the center
		SetDrawEnv xcoord= prel,ycoord= emptyLeft
		DrawLine 0.49,0,0.51,0
	endif

	SetDrawEnv fillpat= 0
	DrawOval 0,0,1,1									// draw big circle
	Variable low
	Variable wid = markers ? 0.05 : 0
	for (i=0;i<N;i+=1)								// draw all of the hkl on the graph
		if (numtype(Xspots[i]+Yspots[i]))
			continue
		endif
//		hklstr = hkl2str(hkls[i][0],hkls[i][1],hkls[i][2])
//		hklstr = ReplaceString(",",hklstr,"")
		hklstr = hkl2IgorBarStr(hkls[i][0],hkls[i][1],hkls[i][2])	// changes negatives to a bar over the number, only for displaying, not printing
		low = lowOrder(hkls[i][0],hkls[i][1],hkls[i][2])	// flag, a low order reflection, show in larger font & bold
		SetDrawEnv xcoord=emptyBottom,ycoord=emptyLeft,textxjust=1,textyjust=1,fstyle=low,fsize=(fntsize+4*low)
		DrawText Xspots[i],Yspots[i]+2*wid,hklstr
		if (wid>0)
			DrawMarker(Xspots[i],Yspots[i],wid,wid,"X",color="0,0,0")
		endif
	endfor

	STRUCT WMSetVariableAction sva				// update the Wulff Net
	sva.eventCode=3
	sva.dval = WulffNetStep
	sva.ctrlName = "setvarWulffStep"
	sva.win = StringFromList(0,winlist("*",";","WIN:1"))
	UpdateWulffNetProc(sva)

	KillWaves/Z W_cross
End
//
Static Function signedAdvance(i)
	Variable i
	i = round(i)
	return (i<0) ? abs(i) : -(i+1)
End
//
Function ZoomStereoButtonProc(B_Struct) : ButtonControl
	STRUCT WMButtonAction &B_Struct
	if (B_Struct.eventCode != 2)
		return 0
	endif

	String list = GetUserData("","","stereoRadius")
	Variable Hz = NumberByKey("Hz",list,"=")
	Variable Kz = NumberByKey("Kz",list,"=")
	Variable Lz = NumberByKey("Lz",list,"=")

	Variable step = 5									// default step size
	step = (B_Struct.eventMod == 8) ? 1 : step		// command-key, 1¡ steps
	step = (B_Struct.eventMod == 4) ? 0.1 : step	// option-key, divide by 0.1¡ steps
	Variable stereoRadius = NumberByKey("stereoRadius",list,"=")
	stereoRadius = numtype(stereoRadius) ? 90 : stereoRadius
	if (stringmatch(B_Struct.ctrlName,"zoomIn"))
		stereoRadius = (B_Struct.eventMod == 2) ? 5 : stereoRadius-step	// zoom full way for shift-key
	elseif (stringmatch(B_Struct.ctrlName,"zoomOut"))
		stereoRadius = (B_Struct.eventMod == 2) ? 90 : stereoRadius+step
	else
		return 1
	endif
	stereoRadius = limit(stereoRadius,5,90)
	Variable d = 2*sin(stereoRadius*PI/180) / (1+cos(stereoRadius*PI/180))
	SetAxis emptyLeft -d,d
	SetAxis emptyBottom -d,d
	String str
	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))		//fill the lattice structure with current values
		sprintf str, "\\Zr200(%s)\r%g\\F'Symbol'°\\F]0",hkl2str(Hz,Kz,Lz),stereoRadius
	else
		if (strlen(xtal.desc))
			sprintf str, "\\JR\\Zr200%s\r\\Zr071#%d   %s\r\\Zr141(%s)\r%g\\F'Symbol'°\\F]0",xtal.desc,xtal.SpaceGroup,getSymString(xtal.SpaceGroup),hkl2str(Hz,Kz,Lz),stereoRadius
		else
			sprintf str, "\\JR\\Zr142#%d   %s\r\\Zr141(%s)\r%g\\F'Symbol'°\\F]0",xtal.SpaceGroup,getSymString(xtal.SpaceGroup),hkl2str(Hz,Kz,Lz),stereoRadius
		endif
	endif
	TextBox/C/N=pole str
	list = ReplaceNumberByKey("stereoRadius",list,stereoRadius,"=")
	SetWindow kwTopWin, userdata(stereoRadius)=list
	return 0
End
//
Static Function UpdateWulffNetProc(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	if (!(sva.eventCode & 3))						// only process on mouse up or enter key
		return 0
	endif
	String str
	Variable dval = sva.dval, step=NumVarOrDefault("WulffNetStep",10)
	if (stringmatch(sva.ctrlName,"setvarWulffStep"))
		str = MakeWulffNet(dval,NumVarOrDefault("WulffNetPhi",0))
	elseif (stringmatch(sva.ctrlName,"setvarWulffPhi"))
		str = MakeWulffNet(step,dval)
	endif
	SetVariable setvarWulffPhi,disable=(step ? 0 : 1)

	String win = sva.win
	Variable takeOff = !step || strlen(str)<1		// flag, take Wulff net off of graph
	Wave wy=$StringFromList(0,str), wx=$StringFromList(1,str), wc=$StringFromList(2,str)
	takeOff = (!WaveExists(wy) || !WaveExists(wx)) || !WaveExists(wx) || takeOff	// if either wave does not exist then takeOff = 1
	if (!WaveExists(wy) || !WaveExists(wx) || !WaveExists(wc))
		Wave wx = WulffNetX, wy=WulffNetY, wc=WulffNetC
	endif
	Variable onGraph = WhichListItem(NameOfWave(wy),TraceNameList(win,";",1))>=0
	if (!takeOff && !onGraph)						// WulffNetY not on graph, put it on
		AppendToGraph/L=emptyLeft/B=emptyBottom/W=$win WulffNetY vs WulffNetX
		ModifyGraph/Z/W=$win zColor(WulffNetY)={WulffNetC,0,1,Grays,0}
	elseif(takeOff && onGraph)					// remove from graph and kill waves
		RemoveFromGraph/Z/W=$win $NameOfWave(wy)
		KillWaves/Z wx,wy,wc
	endif
	return 0
End
//
Static Function/T MakeWulffNet(step,phi,[da])
	Variable step						// angle step size for lines of longitude and lattitude (degree)
	Variable phi						// azimuthal orientation (degree)
	Variable da							// resolution (degree)
	da = ParamIsDefault(da) ? 1 : da			// default resolution is 1¡ (make smaller only if you want to print a very pretty picture)
	phi = !(phi>-180 && phi<180) ? 0 : phi	// set to 0 if phi not in (-180,180)
	if (!(step>0 && step<=90))
		return ""
	endif

	Variable Npnts = ceil(2*(181/da)*(181/step))
	Make/N=(Npnts)/O WulffNetY,WulffNetX, WulffNetC=.91	// initially, X is long and Y is lat
	Variable lat, long, n
	for(n=0, lat=-90+step; lat<(90-da); lat+=step)		// loops to make lines of lattitude
		for (long=-90; long<=90; long+=da, n+=1)
			WulffNetY[n] = long
			WulffNetX[n] = lat
			if (abs(mod(lat,20))<da)
				WulffNetC[n] = 0
			elseif (abs(mod(lat,10))<da)
				WulffNetC[n] = .7
			endif
		endfor
		WulffNetY[n] = NaN
		WulffNetX[n] = NaN
		n += 1
	endfor
	for (long=-90+step; long<=(90-step); long+=step)
		for (lat=-90; lat<=90; lat+=da, n+=1)
			WulffNetY[n] = long
			WulffNetX[n] = lat
			if (abs(mod(long,20))<da)
				WulffNetC[n] = 0
			elseif (abs(mod(long,10))<da)
				WulffNetC[n] = .7
			endif
		endfor
		WulffNetY[n] = NaN
		WulffNetX[n] = NaN
		n += 1
	endfor
	n -= 1												// trim off the last point (it is NaN)
	Redimension/N=(N) WulffNetY,WulffNetX
	if (step>=20)
		WulffNetC = 0									// just a few lines, make them all dark
	endif
	WulffNetX = abs(WulffNetX[p]) < 1e-12 ? 0 : WulffNetX[p]	// This is needed because of a bug in Project
	Project/M=1/C={0,0} WulffNetY, WulffNetX			// Project/M=1 longs, lats
	Wave W_XProjection=W_XProjection, W_YProjection=W_YProjection
	Variable c=cos(phi*PI/180), s=sin(phi*PI/180)		// apply the rotation
	WulffNetX = c*W_YProjection[p] - s*W_XProjection[p]
	WulffNetY = s*W_YProjection[p] + c*W_XProjection[p]
	KillWaves/Z W_XProjection, W_YProjection
	return GetWavesDataFolder(WulffNetY,2)+";"+GetWavesDataFolder(WulffNetX,2)+";"+GetWavesDataFolder(WulffNetC,2)
End
//
Static Function hkl2LongLat(Hz,Kz,Lz,qhat,long,lat,xtal,[hklPerp])
	Wave qhat					// unit vector in direction of qvector
	Variable &long, &lat
	Variable Hz,Kz,Lz			// hkl of the pole
	STRUCT crystalStructure &xtal
	Wave hklPerp				// optional hkl giving the perpendicular direction long=0

	Make/N=3/O/D/FREE xhat, yhat, zhat
	zhat[0] = Hz*xtal.as0+Kz*xtal.bs0+Lz*xtal.cs0	// zhat =  recip x {Hz,Kz,Lz}
	zhat[1] = Hz*xtal.as1+Kz*xtal.bs1+Lz*xtal.cs1
	zhat[2] = Hz*xtal.as2+Kz*xtal.bs2+Lz*xtal.cs2
	if (ParamIsDefault(hklPerp))
		FindPerpVector(zhat,xhat)
	else
		xhat[0] = hklPerp[0]*xtal.as0+hklPerp[1]*xtal.bs0+hklPerp[2]*xtal.cs0	// zhat =  recip x hklPerp
		xhat[1] = hklPerp[0]*xtal.as1+hklPerp[1]*xtal.bs1+hklPerp[2]*xtal.cs1
		xhat[2] = hklPerp[0]*xtal.as2+hklPerp[1]*xtal.bs2+hklPerp[2]*xtal.cs2
	endif

	Cross zhat, xhat
	Wave W_Cross=W_Cross
	yhat = W_Cross
	normalize(xhat)
	normalize(zhat)
	normalize(yhat)

	lat = 90 - acos(min(MatrixDot(zhat,qhat),1))*180/PI
	long = atan2(MatrixDot(yhat,qhat),MatrixDot(xhat,qhat))*180/PI
//	long += ParamIsDefault(hklPerp) ? 0 : 90	// makes hklPerp point to the right
	long += ParamIsDefault(hklPerp) ? 0 : -90	// makes hklPerp point to the left
	return 0
End
//
Static Function FindPerpVector(ref,perp)	// find a vector perpendicular to ref
	Wave ref								// reference vector, find a vec perp to this
	Wave perp							// resulting perpendicular vector, not normalized
	Make/N=3/O/D/FREE hat

	hat = abs(ref)
	WaveStats/M=1/Q hat
	hat = 0
	hat[V_maxLoc] = 1
	Variable dot = MatrixDot(ref,hat)
	perp = ref - dot*hat
	if (norm(perp)<1e-9)				// ref is exacly along x, y, or z
		perp = 0
		perp[V_minloc]=1
	endif
	Cross ref, perp
	Wave W_cross=W_cross
	perp = W_cross
//	KillWaves/Z W_cross
	return 0
End

Static Function lowOrder(h,k,l)		// returns true if (hkl) is a low order reflection
	Variable h,k,l
	h = abs(h)
	k = abs(k)
	l = abs(l)

	Variable zeros = (!h + !k + !l)
	if (zeros >= 2)
		return 1								// (100) type 
	elseif (zeros==1 && (h==k || h==l || k==l))
		return 1								// (110) type
	elseif (h==k && h==l)
		return 1								// (111) type
	endif
	return 0									// not a low order reflection
End


Static Function/WAVE MakeDetectorPolyWaves(outlines,Hz,Kz,Lz,phi,hklPerp,xtal)
	String outlines
	Variable Hz,Kz,Lz
	Variable phi
	Wave hklPerp
	STRUCT crystalStructure &xtal

	Variable Nd=ItemsInList(outlines)
	if (Nd<1)
		return $""
	endif

	Make/N=(Nd,10)/FREE Dpolys=NaN							// detector polygons, sequence of 5 (x,y) pairs
	String wnote=""
	Make/N=3/D/FREE qhat
	Make/N=5/D/FREE Dlongs, Dlats
	Make/N=(4,3)/FREE hklCorner
	String str
	Variable m, i, long,lat
	for (i=0;i<Nd;i+=1)											// loop over each detector polygon
		str = StringFromList(i,outlines)
		if (strlen(str)<1)
			continue
		endif
		Wave rgb = str2vec(StringByKey("rgb",str,"=",":"))
		if (!WaveExists(rgb) || numpnts(rgb)!=3)		// set bad rgbs to Gray
			Make/FREE/N=3 rgb=40000
		elseif (WaveMin(rgb)<0 || WaveMax(rgb)>65535 || numtype(sum(rgb)))
			rgb = 40000	
		endif
		Redimension/W/U rgb
		Wave hkl0 = str2vec(StringByKey("hkl0",str,"=",":"))	// the four corners
		Wave hkl1 = str2vec(StringByKey("hkl1",str,"=",":"))
		Wave hkl2 = str2vec(StringByKey("hkl2",str,"=",":"))
		Wave hkl3 = str2vec(StringByKey("hkl3",str,"=",":"))
		if (!(norm(hkl0)>0 && norm(hkl1)>0 && norm(hkl2)>0 && norm(hkl3)>0))
			continue
		endif
		hklCorner[0][] = hkl0[q]
		hklCorner[1][] = hkl1[q]
		hklCorner[2][] = hkl2[q]
		hklCorner[3][] = hkl3[q]
		Dlongs = NaN
		Dlats = NaN
		Variable H,K,L
		for (m=0;m<4;m+=1)
			qhat = hklCorner[m][p]
			H = hklCorner[m][0]
			K = hklCorner[m][1]
			L = hklCorner[m][2]
			qhat[0] = H*xtal.as0+K*xtal.bs0+L*xtal.cs0	// qmag*qhat =  recip x hkl
			qhat[1] = H*xtal.as1+K*xtal.bs1+L*xtal.cs1
			qhat[2] = H*xtal.as2+K*xtal.bs2+L*xtal.cs2
			normalize(qhat)
			if (!WaveExists(hklPerp))
				hkl2LongLat(Hz,Kz,Lz,qhat,long,lat,xtal)
			else
				hkl2LongLat(Hz,Kz,Lz,qhat,long,lat,xtal,hklPerp=hklPerp)
			endif
			if (lat>=-0.01)
				Dlongs[m] = long
				Dlats[m] = lat
			endif
		endfor
		Dlongs[4] = Dlongs[0]
		Dlats[4] = Dlats[0]
		Dlats = abs(Dlats[p]) < 1e-12 ? 0 : Dlats[p]	// This is needed because of a bug in Project
		Dlongs += phi												// add on the azimuthal rotation
		if (numtype(sum(Dlongs)+sum(Dlats))==0)			// detector appears valid
			Project/M=1 Dlongs, Dlats
			Wave W_XProjection=W_XProjection, W_YProjection=W_YProjection
			Dpolys[i][] = mod(q,2)==0 ? W_XProjection[floor(q/2)] : W_YProjection[floor(q/2)]
			wnote = ReplaceStringByKey("rgb"+num2istr(i),wnote,vec2str(rgb,bare=1,sep=" "),"=")
		endif
	endfor
	KillWaves/Z W_XProjection, W_YProjection
	Note/K Dpolys, wnote
	return Dpolys
End





Function InitStereoGraphicPackage()
	InitLatticeSymPackage()							// used to initialize this package
End