#pragma rtGlobals=1		// Use modern global access method.
#pragma ModuleName=microGeo
#pragma version = 2.65
#define MICRO_VERSION_ORIGINAL
#include  "LatticeSym", version>=3.54
Constant FIRST_PIXEL=1							// use for one-based pixels (old Winview way)

//Static StrConstant GeoWebServer = "www.uni.aps.anl.gov/34ide"
//Static StrConstant GeoWebServer = "top.uni.aps.anl.gov/34ide"
//Static StrConstant GeoWebServer = "www.aps.anl.gov/Sectors/33_34/34ide"
Static StrConstant GeoWebServer = "sector34.xray.aps.anl.gov/34ide"


// version 2.4, changed pixel2XYZ() to return position in microns (not mm)
// version 2.45, added ability to set geometry from web server, and the wire axis
// version 2.51, transfered grain math items from ArryaOf3dOrients.ipf to here
// version 2.55, modified distortion correction.  added offset to correction so that correction at (xc,yc) is zero
// version 2.61, added the FLIP_QX

// with 2.1, I can  reverse the sense of pixel X using the constant CCD_X_REVERSE
Constant CCD_X_REVERSE = 0								// 0 is no reverse, 1 is reverse
//Constant USE_DISTORTION_DEFAULT = 0					// default is to NOT USE distortion
Constant USE_DISTORTION_DEFAULT = 1					// default is TO USE distortion
//Constant cosThetaWire=0.707106781186547, sinThetaWire=0.707106781186548		//  angle wire moves, usually 45 degrees
Strconstant EPICS_PREFIX="34ide:geometry:"				// prefix for the geometry PVs
//
//		put the following line in your procedure file to activate
//	#define FLIP_QX


Menu "Help"
	"Micro Diffraction",DisplayHelpTopic/K=1 "3D-Xray Diffraction"
End
Menu "micro", dynamic
	"-"
	"ÐÐ> micro Panel", /Q, MakeMicroPanel(-1)
	"-"
	"Set Geometry Parameters...",MakeGeometryParametersPanel("")
	help={"Show panel used to set (or view) all of the geometry parameters"}
	Submenu DistortionCorrectionMenuItem(-1)
		DistortionCorrectionMenuItem(0) ,/Q, setUseDistortion(0)		// 0 is OFF
		DistortionCorrectionMenuItem(1) ,/Q, setUseDistortion(1)		// 1 is ON
	End
End
//
Function/S DistortionCorrectionMenuItem(m)				// note, "!\022(" causes a menu item to be checked and disabled
	Variable m
	Init_microGeo()
	NVAR useDistortion=root:Packages:geometry:useDistortion
	if (m==-1)											// -1 is special for the menu header
		return "Distortion Correction is "+SelectString(useDistortion,"OFF","ON")
	elseif (m %^ useDistortion)								// only one is true
		return SelectString(useDistortion,"Start","Stop")+" Using Distortion Correction"
	else
		return "!\022("+SelectString(useDistortion,"Not ","")+"Using Distortion Correction"
	endif
End

Function MICRO_VERSION_ORIGINAL_Func()
	return 0
End

Function MakeMicroPanel(tab)						// makes the main microPanel
	Variable tab
	if (WinType("micropanel")==7)
		DoWindow/F microPanel
		if (!(tab>=0 && tab<=4))
			ControlInfo/W=microPanel tabMicro
			tab = V_Value
		endif
	else
		NewPanel /W=(672,60,991,683)/K=1/N=microPanel
		TabControl tabMicro,pos={0,0},size={317,18},proc=microGeneralTabProc
		TabControl tabMicro,help={"The main micro-beam functions"},userdata(tabNum)="3"
//		TabControl tabMicro,tabLabel(0)="Index",tabLabel(1)="E-W scan"
		TabControl tabMicro,tabLabel(0)="Index",tabLabel(1)="E scans"
		TabControl tabMicro,tabLabel(2)="3D-arrays",tabLabel(3)="Geo"
		TabControl tabMicro,tabLabel(4)="Xtal"
		MoveWindowToCorner("microPanel","tr")	// re-position to the top right corner
	endif
	tab = (tab>=0 && tab<=4) ? tab : 3				// default to geometry
	TabControl tabMicro,value=tab					// set control to tab
	STRUCT WMTabControlAction  tca				// open panel to the geometry tab
	tca.ctrlName="tabMicro"
	tca.win="microPanel"
	tca.eventCode=2
	tca.tab = tab
	tca.ctrlRect.bottom = 18
	microGeneralTabProc(tca)
End
//
Function microGeneralTabProc(tca) : TabControl		// changes the panel for each tab when a tab is selected
	STRUCT WMTabControlAction &tca
	switch( tca.eventCode )
		case -1:									// control killed
			return 1
		case 2:										// mouse up
			Variable tab = tca.tab
			break
	endswitch
	if (NumberByKey("tabNum",GetUserData(tca.win,"tabMicro","tabNum"),"=")==tab)
		return 0
	endif
	TabControl $tca.ctrlName, win=$tca.win,userdata(tabNum)=num2istr(tab)
	String win = GetUserData(tca.win,"tabMicro","panelName")
	if (WinType(win)==7 && strlen(win))
		KillWindow $win
	endif

	String fillFunctionList = "Index;EWscan;Arrays3d;Geometry;Lattice"
	String funcName = "Fill"+StringFromList(tab,fillFunctionList)+"ParametersPanel"
	if (exists(funcName)!=6)
		funcName = "Load"+funcName[4,Inf]
	endif
	if (exists(funcName)==6)
		FUNCREF FillGeometryParametersPanel FillParametersPanel = $funcName
		win = tca.win+FillParametersPanel("","microPanel",45,2+tca.ctrlRect.bottom)
	else
		win = ""
	endif
	TabControl $tca.ctrlName, win=$tca.win,userdata(panelName)=win
	SetActiveSubwindow microPanel
	return 0
End
//
Function/T LoadLatticeParametersPanel(strStruct,hostWin,left,top)
	String strStruct									// optional passed value of xtal structure, this is used if passed
	String hostWin									// name of home window
	Variable left, top									// offsets from the left and top

	SetWindow kwTopWin,userdata(LatticePanelName)=hostWin+"#LatticePanel"
	NewPanel/K=1/W=(left,top,left+221,top+365)/HOST=$hostWin
	ModifyPanel frameStyle=0, frameInset=0
	RenameWindow #,LatticePanel
	Button buttonInitLattice,pos={33,31},size={150,50},title="Init Lattice\rpackage",proc=LoadPackageButtonProc
	Button buttonInitLattice,help={"Load lattice symmetry procedures"}
	return "#LatticePanel"
End
//
Function/T LoadIndexParametersPanel(strStruct,hostWin,left,top)
	String strStruct									// optional passed value of xtal structure, this is used if passed
	String hostWin									// name of home window
	Variable left, top									// offsets from the left and top

	SetWindow kwTopWin,userdata(IndexPanelName)=hostWin+"#IndexPanel"
	NewPanel/K=1/W=(left,top,left+221,top+365)/HOST=$hostWin
	ModifyPanel frameStyle=0, frameInset=0
	RenameWindow #,IndexPanel
	Button buttonInitIndex,pos={33,31},size={150,50},title="Init Indexing\rpackage",proc=LoadPackageButtonProc
	Button buttonInitIndex,help={"Load Package for Fitting Peaks and Indexing patterns"}
	return "#IndexPanel"
End
//
Function/T LoadEWscanParametersPanel(strStruct,hostWin,left,top)
	String strStruct									// optional passed value of xtal structure, this is used if passed
	String hostWin									// name of home window
	Variable left, top									// offsets from the left and top

	SetWindow kwTopWin,userdata(EWscanPanelName)=hostWin+"#EWscanPanel"
	NewPanel/K=1/W=(left,top,left+221,top+365)/HOST=$hostWin
	ModifyPanel frameStyle=0, frameInset=0
	RenameWindow #,EwscanPanel
	Button buttonInitEWscan,pos={33,31},size={150,50},title="Init E-W scan\rpackage",proc=LoadPackageButtonProc
	Button buttonInitEWscan,help={"Load Package for processing E-W scans"}
	return "#EWscanPanel"
End
//Function/T LoadEWscanParametersPanel(strStruct,hostWin,left,top)
//	String strStruct									// optional passed value of xtal structure, this is used if passed
//	String hostWin									// name of home window
//	Variable left, top									// offsets from the left and top
//
//	SetWindow kwTopWin,userdata(EWscanPanelName)=hostWin+"#EWscanPanel"
//	NewPanel/K=1/W=(left,top,left+221,top+365)/HOST=$hostWin
//	ModifyPanel frameStyle=0, frameInset=0
//	RenameWindow #,EwscanPanel
//	if (strlen(WinList("EnergyWireScans.*", "", "WIN:128")))
//		TitleBox titleLoadedEnergyWire,pos={18,37},size={180,40},title="\\JCmacros already loaded,\rsee the menubar"
//		TitleBox titleLoadedEnergyWire,fSize=16,frame=0
//	else
//		Button buttonInitEWscan,pos={33,31},size={150,50},title="Init E-W scan\rpackage",proc=LoadPackageButtonProc
//		Button buttonInitEWscan,help={"Load Package for processing E-W scans"}
//	endif
//	return "#EWscanPanel"
//End
//
Function/T LoadArrays3dParametersPanel(strStruct,hostWin,left,top)
	String strStruct									// optional passed value of xtal structure, this is used if passed
	String hostWin									// name of home window
	Variable left, top									// offsets from the left and top

	SetWindow kwTopWin,userdata(Arrays3dPanelName)=hostWin+"#Arrays3dPanel"
	NewPanel/K=1/W=(left,top,left+221,top+365)/HOST=$hostWin
	ModifyPanel frameStyle=0, frameInset=0
	RenameWindow #,Arrays3dPanel
	if (strlen(WinList("ArrayOf3dOrients.*", "", "WIN:128")))
		TitleBox titleLoadedArray3d,pos={18,37},size={180,40},title="\\JCmacros already loaded,\rsee the menubar"
		TitleBox titleLoadedArray3d,fSize=16,frame=0
	else
		Button buttonInitArrays3d,pos={33,31},size={150,50},title="Init 3d Arrays\rpackage",proc=LoadPackageButtonProc
		Button buttonInitArrays3d,help={"Processes 3D array of matricies into viewable of rotataion angles, dislocation tensors, GND, etc."}
	endif
	return "#Arrays3dPanel"
End
//
Function LoadPackageButtonProc(ctrlName) : ButtonControl
	String ctrlName

	if (stringmatch(ctrlName,"buttonInitLattice"))
		Execute/P "INSERTINCLUDE  \"LatticeSym\", version>=3.32";Execute/P "COMPILEPROCEDURES ";Execute/P "InitLatticeSymPackage()"
	elseif (stringmatch(ctrlName,"buttonInitIndex"))
		Execute/P "INSERTINCLUDE  \"Indexing\", version>=2.59";Execute/P "COMPILEPROCEDURES ";Execute/P "InitIndexingPackage()"
	elseif (stringmatch(ctrlName,"buttonInitIndexLots"))
		Execute/P "DELETEINCLUDE  \"Indexing\", version>=2.59"
		Execute/P "INSERTINCLUDE  \"IndexLots\", version>=2.17"
		Execute/P "COMPILEPROCEDURES ";Execute/P "initSymmetryOperations()"
//		Execute/P "INSERTINCLUDE  \"IndexLots\", version>=2.17";Execute/P "COMPILEPROCEDURES ";Execute/P "initSymmetryOperations()"
	elseif (stringmatch(ctrlName,"buttonInitEWscan"))
		Execute/P "INSERTINCLUDE  \"EnergyWireScans\", version>=0.981";Execute/P "COMPILEPROCEDURES "
	elseif (stringmatch(ctrlName,"buttonInitArrays3d"))
		Execute/P "INSERTINCLUDE  \"ArrayOf3dOrients\", version>=1.8";Execute/P "COMPILEPROCEDURES "
	else
		return 1
	endif

	TabControl tabMicro,value=3				// set control to tab
	STRUCT WMTabControlAction  tca			// open panel to the geometry tab
	tca.ctrlName="tabMicro"
	tca.win="microPanel"
	tca.eventCode=2
	tca.ctrlRect.bottom = 18
	tca.tab = 3
	microGeneralTabProc(tca)
	return 0
End
//
Static Function MoveWindowToCorner(win,corner)
	String win
	String corner			//  one of "TL", "TR", "BL", "BR"
	corner = UpperStr(corner)
	win = SelectString(strlen(win)>0,StringFromList(0,WinList("*",";","WIN:71")),win)
	if (strlen(win)<1 || WhichListItem(corner,"TL;TR;BL;BR")<0)
		return 1
	endif

	Variable leftScreen,topScreen,rightScreen,bottomScreen,  top, left
	if (stringmatch(IgorInfo(2),"Windows"))					// use the inside of the Igor Application Window on Windows
		GetWindow kwFrameInner,wsize
		leftScreen = V_left ;	rightScreen = V_right
		topScreen = V_top ;		bottomScreen = V_bottom
	else	
		String list = StringByKey("SCREEN1",IgorInfo(0))		// for a Mac, use the whole screen (except menubar)
		Variable i=strsearch(list,"RECT=",0)
		sscanf list[i+5,Inf], "%g,%g,%g,%g", leftScreen,topScreen,rightScreen,bottomScreen
		topScreen += 44										// leave room for menubar on Mac
	endif

	GetWindow $win,wsize										// size of window to move
	Variable width=V_right-V_left, height=V_bottom-V_top		// get size of the window
	top = stringmatch(corner[0],"T") ? topScreen : (bottomScreen-height)	// desired top left position of window
	left = stringmatch(corner[1],"L") ? leftScreen : (rightScreen-width)
	MoveWindow/W=$win left,top,left+width,top+height
	return 0
End
//Function test()
//	Display/N=testWindow
//	Variable delay=0.5
//	microGeo#MoveWindowToCorner("testWindow","tl")
//	Sleep/B/S delay
//	microGeo#MoveWindowToCorner("testWindow","bl")
//	Sleep/B/S delay
//	microGeo#MoveWindowToCorner("testWindow","br")
//	Sleep/B/S delay
//	microGeo#MoveWindowToCorner("testWindow","tr")
//	Sleep/B/S delay
//	DoWindow/K testWindow
//End













// to convert from a ROI pixel (zero based) to full unbinned chip use:
//
//	px = (startx-1) + px*groupx + (groupy-1)/2	// pixel is zero based here
//	py = (starty-1) + py*groupy + (groupy-1)/2




//  Note about the origin.	March 22, 2006
//
// Igor uses (0,0) as the first point in the image.   Wenge, uses (1,1) as the first point
// So:
// 
// changed peakcorrection2() to correct for this, nothing else.  Be Careful of meaning of geo.xcent and geo.ycent
// Wenge's program computes them with origin at (1,1), I use the origin (0,0).  So I should take his values minus 1.



// Beam line corodinate system (X=out door, Y=up, Z=downstream)
// Wenge corodinate system (X=out door, Y=upstream, Z=up)
//
// Xbeam = Xwenge
// Ybeam = Zwenge
// Zbeam = -Ywenge
//
// Xwenge = Xbeam
// Ywenge = -Zbeam
// Zwenge = Ybeam
//
Function YZ2F(y,z)		// F = -Y*cos(angle) + Z*sin(angle)
	Variable Y,Z
	Variable cosTheta = NumVarOrDefault("root:Packages:geometry:cosThetaWire",cos(PI/4))
	Variable sinTheta = NumVarOrDefault("root:Packages:geometry:sinThetaWire",sin(PI/4))
	return -Y*cosTheta + Z*sinTheta
End
Function YZ2H(Y,Z)		// H =  Y*sin(angle) + Z*cos(angle)
	Variable Y,Z
	Variable cosTheta = NumVarOrDefault("root:Packages:geometry:cosThetaWire",cos(PI/4))
	Variable sinTheta = NumVarOrDefault("root:Packages:geometry:sinThetaWire",sin(PI/4))
	return Y*sinTheta + Z*cosTheta
End
Function HF2Y(H,F)	// Y =  H*sin(angle) - F*cos(angle)
	Variable H,F
	Variable cosTheta = NumVarOrDefault("root:Packages:geometry:cosThetaWire",cos(PI/4))
	Variable sinTheta = NumVarOrDefault("root:Packages:geometry:sinThetaWire",sin(PI/4))
	return H*sinTheta - F*cosTheta
End
Function HF2Z(H,F)	// Z =  H*cos(angle) + F*sin(angle)
	Variable H,F
	Variable cosTheta = NumVarOrDefault("root:Packages:geometry:cosThetaWire",cos(PI/4))
	Variable sinTheta = NumVarOrDefault("root:Packages:geometry:sinThetaWire",sin(PI/4))
	return H*cosTheta + F*sinTheta
End

//	This is only called by GenericWaveNoteInfo() in Utility_JZT.ipf
Function/T MoreWaveNoteInfo(ww,list)		// additions to the list, only called by GenericWaveNoteInfo() in Utility_JZT.ipf
	Wave ww
	String list
	if (!WaveExists(ww))
		return list
	endif

	Variable Y,Z,H,F
	Y = NumberByKey("Y1",list,"=")
	Z = NumberByKey("Z1",list,"=")
	H = NumberByKey("H1",list,"=")			// get H & F to test for their existance
	F = NumberByKey("F1",list,"=")
	if (!numtype(Y+Z) && numtype(H+F))		// if Y1 & Z1 exist, and H1 & F1 do not
		H = YZ2H(Y,Z)		// H =  Y*sin(angle) + Z*cos(angle)
		F = YZ2F(Y,Z)		// F = -Y*cos(angle) + Z*sin(angle)
		list = ReplaceNumberByKey("H1",list,H,"=")
		list = ReplaceNumberByKey("F1",list,F,"=")
	endif

	Y = NumberByKey("Y2",list,"=")
	Z = NumberByKey("Z2",list,"=")
	H = NumberByKey("H2",list,"=")
	F = NumberByKey("F2",list,"=")
	if (!numtype(Y+Z) && numtype(H+F))		// if Y2 & Z2 exist, and H2 & F2 do not
		H = YZ2H(Y,Z)
		F = YZ2F(Y,Z)
		list = ReplaceNumberByKey("H2",list,H,"=")
		list = ReplaceNumberByKey("F2",list,F,"=")
	endif

	return list
End




Structure microGeometry	// structure definition
	int32 NxCCD, NyCCD		// # of un-binned pixels in CCD, also called xdim, ydim
	double dpsx, dpsy			// size of CCD (mm)
	double xcent, ycent		// center of CCD in un-binned 1 based pixels directly over the Si calibration std.
	double ddOffset			// dd = CCDy + ddOffset, CCDy is CCD height reported by Heidenhain (µm)
	double dd					// height of CCD (µm)
	double xbet, xgam			// incident beam tilt (degrees)
	double xalfd, xbetd		// CCD tilt (degrees)
	double ki[3]				// unit vector in inpt direction
	double rded00, rded01, rded02	// rotation matrix for detector
	double rded10, rded11, rded12
	double rded20, rded21, rded22
	double rde00, rde01, rde02		// rotatation matrix from Wenge coord to Ideal BeamLine system
	double rde10, rde11, rde12
	double rde20, rde21, rde22
	STRUCT wireGeometry wire
EndStructure
//
Structure wireGeometry		// structure definition
	double F					// F of the wire for the Si calibration wire scan (µm)
	double H0					// H of wire centered on the direct beam (µm)
	double Hyc					// H of wire centered on the vertical over the Si (µm)
	double X					// X of wire for Si (µm)
	double dia					// diameter of wire (µm)
	double xyzSi[3]			// wire beam line coordinates that would put wire center at the Si position, computed (µm)
	double xyz0[3]			// wire beam line coordinates when wire is on the incident beam & Fwire is F, computed (µm)
	double axis[3]				// unit vector in direction of wire axis
EndStructure
//
Function printGeometry(geo)
	STRUCT microGeometry &geo
	printf "using geomery parameters:\r"
	printf "	NxCCD=%d, NyCCD=%d		// number of un-binned pixels in CCD, also called xdim, ydim\r",geo.NxCCD,geo.NyCCD
	printf "	dpsx=%g, dpsy=%g			// size of CCD (mm)\r",geo.dpsx,geo.dpsy
	printf "	xcent=%.3f, ycent=%.3f	// center of CCD in un-binned 1 based pixels\r",geo.xcent,geo.ycent
	printf "	dd=%g						// height of CCD (µm)\r",geo.dd
	printf "	ddOffset=%g					// dd = CCDy + ddOffset (µm)\r",geo.ddOffset
	printf "	xbet=%.5f, xgam=%.5f		// incident beam tilt (degrees)\r",geo.xbet,geo.xgam
	printf "	xalfd=%.5f, xbetd=%.5f	// CCD tilt (degrees)\r",geo.xalfd,geo.xbetd
	printf "	ki = {%.4f, %.4f, %.4f}	// incident beam direction (Wenge coordinates)\r",geo.ki[0],geo.ki[1],geo.ki[2]
	if (NumVarOrDefault("root:Packages:geometry:printVerbose",0))
		printf "			{%+.6f, %+.6f, %+.6f}	// rotation matrix for detector (Wenge coordinates)\r",geo.rded00, geo.rded01, geo.rded02
		printf "	rded =	{%+.6f, %+.6f, %+.6f}\r",geo.rded10, geo.rded11, geo.rded12
		printf "			{%+.6f, %+.6f, %+.6f}\r",geo.rded20, geo.rded21, geo.rded22
		printf "			{%+.6f, %+.6f, %+.6f}	// rotation matrix for incident beam (Wenge coordinates)\r",geo.rde00, geo.rde01, geo.rde02
		printf "	rde =	{%+.6f, %+.6f, %+.6f}\r",geo.rde10, geo.rde11, geo.rde12
		printf "			{%+.6f, %+.6f, %+.6f}\r",geo.rde20, geo.rde21, geo.rde22
//		printf "			{%+.4f, %+.4f, %+.4f}	// rotation matrix for detector (Wenge coordinates)\r",geo.rded00, geo.rded01, geo.rded02
//		printf "	rded =	{%+.4f, %+.4f, %+.4f}\r",geo.rded10, geo.rded11, geo.rded12
//		printf "			{%+.4f, %+.4f, %+.4f}\r",geo.rded20, geo.rded21, geo.rded22
//		printf "			{%+.4f, %+.4f, %+.4f}	// rotation matrix for incident beam (Wenge coordinates)\r",geo.rde00, geo.rde01, geo.rde02
//		printf "	rde =	{%+.4f, %+.4f, %+.4f}\r",geo.rde10, geo.rde11, geo.rde12
//		printf "			{%+.4f, %+.4f, %+.4f}\r",geo.rde20, geo.rde21, geo.rde22
	endif
	printf "		for the wire:\r"
	printf "	F=%.2f							// F of the wire for the Si calibration wire scan (µm)\r",geo.wire.F
	printf "	H0=%.2f						// H of wire centered on the direct beam (µm)\r",geo.wire.H0
	printf "	Hyc=%.2f						// H of wire centered on the vertical over the Si (µm)\r",geo.wire.Hyc
	printf "	X=%g								// X of wire for Si (µm)\r",geo.wire.X
	printf "	diameter=%.2f						// diameter of wire (µm)\r",geo.wire.dia
	printf "	wire axis direction = {%.6f, %.6f, %.6f}	// wire (PM500 wire coordinates) on the incident beam (µm)\r",geo.wire.axis[0],geo.wire.axis[1],geo.wire.axis[2]
	printf "	Si = {%.2f, %.2f, %.2f}	// wire (PM500 wire coordinates) of the Si position (µm)\r",geo.wire.xyzSi[0],geo.wire.xyzSi[1],geo.wire.xyzSi[2]
	printf "	@beam = {%.2f, %.2f, %.2f}	// wire (PM500 wire coordinates) on the incident beam (µm)\r",geo.wire.xyz0[0],geo.wire.xyz0[1],geo.wire.xyz0[2]
End


// These two routines convert betweeen Wenge and Ideal Beam Line coordinate systems.  They are used to get to/from Indexing program (Euler)
Function wenge2IdealBeamline(geo,qW,qBL)
	STRUCT microGeometry &geo
	Wave qW, qBL
	Variable q0,q1,q2
	q0 = qW[0]*geo.rde00 + qW[1]*geo.rde01 + qW[2]*geo.rde02	// rotate qW from Wenge to Ideal Beam-line
	q1 = qW[0]*geo.rde10 + qW[1]*geo.rde11 + qW[2]*geo.rde12	//   Ideal Beam LIne has incident beam along {0,0,1}
	q2 = qW[0]*geo.rde20 + qW[1]*geo.rde21 + qW[2]*geo.rde22	// 	Wenge has incident near {0,-1,0}, but tilted a bit
	qBL = {q0,q1,q2}
End
//
Function IdealBeamLine2wenge(geo,qBL,qW)
	STRUCT microGeometry &geo
	Wave qBL, qW
	Variable q0,q1,q2
	q0 = qBL[0]*geo.rde00 + qBL[1]*geo.rde10 + qBL[2]*geo.rde20	// rotate qBL from Ideal Beam-line to Wenge
	q1 = qBL[0]*geo.rde01 + qBL[1]*geo.rde11 + qBL[2]*geo.rde21	//   Ideal Beam LIne has incident beam along {0,0,1}
	q2 = qBL[0]*geo.rde02 + qBL[1]*geo.rde12 + qBL[2]*geo.rde22	// 	Wenge has incident near {0,-1,0}, but tilted a bit
	qW = {q0,q1,q2}
End

Function EulerMatrix(alpha,bet,gam)	// Make the rotation matrix M_Euler from Euler angles (degree)
	Variable alpha,bet,gam				// Euler angles, passed as degrees
	alpha *= PI/180					// need radians internally
	bet *= PI/180
	gam *= PI/180
	Make/N=(3,3)/O/D M_Euler
	if (alpha<=-2*PI || alpha>=2*PI || bet<=-2*PI || bet>=2*PI || gam<=-2*PI || gam>=2*PI)
		M_Euler = NaN
		return 1
	endif
	Variable ca=cos(alpha), cb=cos(bet), sa=sin(alpha), sb=sin(bet)
	Variable cg = cos(gam), sg = sin(gam)
	M_Euler[0][0] = cg*cb*ca - sg*sa
	M_Euler[1][0] = -sg*cb*ca - cg*sa
	M_Euler[2][0] = sb*ca
	M_Euler[0][1] = cg*cb*sa + sg*ca
	M_Euler[1][1] = -sg*cb*sa + cg*ca
	M_Euler[2][1] = sb*sa
	M_Euler[0][2] = -cg*sb
	M_Euler[1][2] = sg*sb
	M_Euler[2][2] = cb
	return 0
End

Function MovingSiOrigin(geo)			// this calculates Wenge's moving origin
	// see the procedure "calctotalrange" in "ThreeDimX_RayMicroscopy.pro" for the origin of this calculation
	STRUCT microGeometry &geo
	Variable SiHstart = -7030							// this one parameter is not in the geometry strucutre

	Variable cosTheta = NumVarOrDefault("root:Packages:geometry:cosThetaWire",cos(PI/4))
	Variable sinTheta = NumVarOrDefault("root:Packages:geometry:sinThetaWire",sin(PI/4))
	Variable ddc = geo.dd*1e3							// 40220
	Variable wirepitchy = geo.dpsy/geo.NxCCD * 1000
	Variable yoff = 0									// full chip image
	Variable CCDypitch=24.0, CCDysize=24*2084		// un-binned images, yes, these are hard coded

	Variable startwirevertical, startwirehorizontal, sourcemin
	startwirevertical=(SiHstart-geo.wire.H0)*cosTheta
	startwirehorizontal=(SiHstart-geo.wire.Hyc)*sinTheta
	sourcemin=CCDYsize-(1+geo.ycent-yoff)*wirepitchy-yoff*CCDypitch-ddc*(CCDYsize-(1+geo.ycent-yoff)*wirepitchy-yoff*CCDypitch-startwirehorizontal)/(ddc-startwirevertical)
	return -sourcemin
End
//Function TestMovingSiOrigin()
//	STRUCT microGeometry geo
//	FillGeometryStructDefault(geo)
//	geo.dd = 40.22
//	geo.ycent = 923.22
//	geo.wire.H0 = -7277
//	geo.wire.Hyc = -6858.27
//	Variable depth = MovingSiOrigin(geo)
//	printf "Si is %g µm from the moving origin\r",depth
//	// using these number, the results are:     Si is 243.067 µm from the moving origin
//End



//	Function test()
//		STRUCT microGeometry geo
//		SVAR geoStr=geoStr
//		StructGet/S/B=2 geo, geoStr		// retrieve the information from geoStr (also gets all of geo.wire.xxx)
//	
//		print SelectString(USE_DISTORTION_DEFAULT,"no distortion","using distortion")
//	//	Variable px,py
//		Variable depth, edge = 1
//		Variable i, N=1000
//	
//		Variable ccdy,ccdz
//		Variable ccdy0, dccdy=27					// ccdy ranges from ccdy0±dccdy
//		Variable ccdz0,dccdz						// ccdz ranges from ccdz0±dccdz
//		Variable ccdz1,ccdz2						// ends of z range on ccd
//		Wave ccd=root:Packages:geometry:pixel2depth_ccd
//		pixel2XYZ(geo,1000,1000,ccd)			// returns xyz of pixel on CCD relative to the Si position (Wenge coordinates)
//		ccdy0 = ccd[2]								// change from Wenge to beamline
//		pixel2XYZ(geo,0,0,ccd)
//		ccdz1 = -ccd[1]								// change from wenge coord to beam line
//		pixel2XYZ(geo,2000,2000,ccd)
//		ccdz2 = -ccd[1]
//		ccdz0 = (ccdz1+ccdz2)/2					// center of allowed ccdz values
//		dccdz = abs(ccdz1-ccdz2)/2				// range of allowed ccdz values
//	
//		Wave Xw = $MakeUnique3Vector($"")
//		Xw = {.5, -6119, -800}
//	
//	 	Variable timer = startMSTimer
//		for (i=0;i<N;i+=1)
//	//		px = round(enoise(1000))+1000
//	//		py = round(enoise(1000))+1000
//	//		depth = pixel2depth(geo,px,py,Xw,edge)
//			ccdy = ccdy0+enoise(dccdy)
//			ccdz = ccdz0+enoise(dccdz)
//			depth = pixelXYZ2depth(geo,ccdy,ccdz,Xw,edge)
//		endfor
//		Variable seconds = stopMSTimer(timer)*1e-6
//		Variable total = (2084^2)/N * seconds
//		printf "execution took %g sec,     a whole image will take %s\r",seconds,Secs2Time(total,5,1)
//		KillWaves Xw
//	End


//Function xxxTemp()
//	STRUCT microGeometry geo
//	SVAR geoStr=geoStr
//	StructGet/S/B=2 geo, geoStr		// retrieve the information from geoStr (also gets all of geo.wire.xxx)
//
//	Wave ki = $MakeUnique3Vector($"")
//	ki = geo.ki[p]
//	printf "ki = {%g, %g, %g},   |ki| = %g\r",ki[0],ki[1],ki[2],norm(ki)
//	KillWaves ki
//End



// Given the three points, P=pixel, S=source, W=wire:
// The next three routines do calculations for the unknown point when the other two are known
//
Function calcWireH(geo,depth,py,pz,edge)	// calc H of wire that is tangent to line from depth on beam to pixel on ccd
	STRUCT microGeometry &geo
	Variable depth						// current Z coordinate in sample relative to the Si position, not distance from Si (µm)
	Variable py,pz						// y and z component of pixel (xyz) position in beam line coords (origin at Si)
	Variable edge						// 1=leading edge (usual),  0=trailing edge

	Variable R=(geo.wire.dia)/2, Fo=geo.wire.F		// wire radius, Fo of the wire
	Variable kiy=geo.ki[2], kiz=-geo.ki[1]			// convert from Wenge coords to beam line
	Variable sy,sz									// y,z coordinates of sample point beamline coords (origin at Si)
	sy = kiy/kiz*depth								// project along ki a distance depth,  sz = kiz*depth/kiz,  sx = kix*depth/kiz
	sz = depth

	Variable Xoy, Xoz								// point where wire intersects incident beam relative to Si position (µm)
	Xoy = geo.wire.xyz0[1] - geo.wire.xyzSi[1]		// origin at Si position
	Xoz = geo.wire.xyz0[2] - geo.wire.xyzSi[2]

	// eqn of line that wire travels along is:   y = z*tan(theta) + b
	Variable tanTheta,b		
	Variable cosTheta = NumVarOrDefault("root:Packages:geometry:cosThetaWire",cos(PI/4))
	Variable sinTheta = NumVarOrDefault("root:Packages:geometry:sinThetaWire",sin(PI/4))
	tanTheta = sinTheta/cosTheta					// = tan(theta), the slope
	b = Xoy - Xoz*tanTheta							// b comes from point that wire passes through on the incident beam

	// need eqn of line parallel to line from S to P, but moved over by distance r (=wire radius):     z = m*y + beta1
	Variable sinphi = (py-sy) / sqrt((py-sy)^2 + (pz-sz)^2)	// sin(phi), phi is angle of ray to pixel measured from z toward y
	Variable m = (pz-sz)/(py-sy)					// slope of line (it is translated, but parallel to P,S line
	Variable beta0 = sz-m*sy						// intercept for line through center of parallel to line from S to P
	Variable beta1 = (edge ? beta0-r/sinphi : beta0+r/sinphi)	// intercept of line from S to P tangent to leading (or trailing) edge

	Variable z = (m*b+beta1) / (1-m*tanTheta)	// z of wire at intersection of two lines
	z += geo.wire.xyzSi[2]							// convert from Si at origin to PM500 origin

	Variable H = (z-Fo*sinTheta)/cosTheta			// from:   Z =  H*cos(angle) + F*sin(angle)
	return H										// H of the wire, the returned value (wire origin)
End
//
//
// for a ray coming from depth, and tangent to wire, where will it hit detector, returns CCD(y,z) position as complex value
Function/C depthWire2pixelXYZ(geo,depth,Xw,xyzMin,xyzMax,edge)
	STRUCT microGeometry &geo
	Variable depth									// depth along incident beam, relative to the Si position (µm)
	Wave Xw										// wire (x,y,z) beam line coordinates in wire units (µm)
	Wave xyzMin,xyzMax							// two points on the detector (assumed to be ends) relative to Si
	Variable edge									// 1=leading edge (usual),  0=trailing edge
	Variable testing = (WhichListItem("TestdepthWire2pixelXYZ", GetRTStackInfo(0))>=0)

	Variable R = (geo.wire.dia)/2					// radius of wire (µm)
	Variable kiy=geo.ki[2], kiz=-geo.ki[1]			// components of ki[] in beam line coords (not Wenge coords)
	Variable sz = depth, sy = kiy/kiz * depth		// source position, relative to the Si position (beam line coords)

	Variable vy,vz,v								// (vy,vz) vector from source (at depth) to center of wire
	vy = (Xw[1] - geo.wire.xyzSi[1]) - sy
	vz = (Xw[2] - geo.wire.xyzSi[2]) - sz
	v = sqrt(vy*vy + vz*vz)						// length of v

	Variable phi0 = atan2(vy,vz)					// angle from positive z to tangent ray measured from z toward y
	Variable dphi =asin(R/v)						// angle between v and ray at tangent point (always positive)
	Variable phi = phi0 + (edge ? -dphi : dphi)
	Variable cotPhi = 1/tan(phi)

	// line from P to S is z = y*cot(phi) + b
	Variable b = sz - sy*cotPhi

	// line along surface of detector, line between xyzMin[] & xyzMax[],   y = m*z + beta1
	Variable m = (xyzMax[1]-xyzMin[1])/(xyzMax[2]-xyzMin[2])
	Variable beta1 = xyzMin[1] - m*xyzMin[2]

	Variable ccdy, ccdz								// y, z parts of the pixel location relative to Si (beam line coords)
	ccdy = (m*b + beta1) / (1-m*cotPhi)			// y value of intersection
	ccdz = ccdy*cotPhi + b							// from line P to S,   z = y*cot(phi) + b

	if (testing)
		Variable cosTheta = NumVarOrDefault("root:Packages:geometry:cosThetaWire",cos(PI/4))
		Variable sinTheta = NumVarOrDefault("root:Packages:geometry:sinThetaWire",sin(PI/4))
		printf "Hwire=%g, Fwire=%g,   Xw=(%g, %g, %g)\r",Xw[1]*sinTheta + Xw[2]*cosTheta, -Xw[1]*cosTheta + Xw[2]*sinTheta,Xw[0], Xw[1], Xw[2]
		printf "s = (%g, %g),   v = (%g, %g),   R=%gµm,  |wire-sample| = %gµm\r",sy,sz,vy,vz,R,v
		printf "line from source to tangent is:   z = %g*y + %g\r",cotPhi,b
		printf "line between xyzMin & xyzMax is:   y = %g*z + %g\r",m,beta1
		// check that ccdy,ccdz satisfy both equations
		printf " is z-y*cot(phi) = %g == %g=b?, Æ=%g\r",ccdz-ccdy*cotPhi,b,ccdz-ccdy*cotPhi-b
		printf " is y-m*z = %g == %g=beta?,  Æ=%g\r",ccdy-m*ccdz,beta1,ccdy-m*ccdz-beta1
	endif
	return cmplx(ccdy,ccdz)
End
//
//
// returns depth (starting point) of ray ending at pixel (px,py) that is tangent to leading edge of the wire
// the returned depth is relative to the Si position (µm)
Function pixel2depth(geo,px,py,Xw,edge)
	STRUCT microGeometry &geo
	Variable px,py						// uncorrected unbinned pixel position (pixels) (zero based)
	Wave Xw							// wire (x,y,z) beam line coordinates in wire units (µm)
	Variable edge						// 1=leading edge (usual),  0=trailing edge
	Variable depth=NaN					// computed depth relative to the Si position (µm)

	Wave ccd=root:Packages:geometry:pixel2depth_ccd

	pixel2XYZ(geo,px,py,ccd)			// returns xyz of pixel on CCD (µm) relative to the Si position (Wenge coordinates)
	Variable swap
	swap = ccd[1]							// convert xyz from Wenge to beam line coordinates
	ccd[1] = ccd[2]
	ccd[2] = -swap

	depth = pixelXYZ2depth(geo,ccd[1],ccd[2],Xw,edge)
	return depth						// relative to the Si position (µm)
End
//
// Returns depth (starting point) of ray ending at point xyz (probably a pixel location) that is tangent
// to leading edge of the wire.  The returned depth is relative to the Si position (µm)
// depth is the change in z value from the Si position, not the distance
Function pixelXYZ2depth(geo,ccdy,ccdz,Xw,edge)
	STRUCT microGeometry &geo
	Variable ccdy, ccdz					// y, z parts of the xyz location of pixel (beam line coords)(µm), origin at Si
	Wave Xw								// wire (x,y,z) beam line coordinates in wire PM500 units (µm)
	Variable edge							// 1=leading edge (usual),  0=trailing edge
	Variable depth							// computed depth relative to the Si position (µm)

	Variable kiy,kiz						// components of ki[] in beam line coords (not Wenge coords)
	kiy = geo.ki[2]
	kiz = -geo.ki[1]

	Variable Xwy,Xwz						// wire center relative to the Si position (using beam line coords)
	Xwy = Xw[1] - geo.wire.xyzSi[1]
	Xwz = Xw[2] - geo.wire.xyzSi[2]

	Variable Xoy, Xoz						// point where wire intersects incident beam relative to Si position (µm)
	Xoy = geo.wire.xyz0[1] - geo.wire.xyzSi[1]
	Xoz = geo.wire.xyz0[2] - geo.wire.xyzSi[2]

	Variable vy,vz, v
	vy = Xwy-ccdy						// (vy,vz) vector from pixel on ccd toward the wire center
	vz = Xwz-ccdz
	v = sqrt(vy*vy + vz*vz)				// length of vector (0,vy,vz) = Xw[]-P[]

	Variable R=(geo.wire.dia)/2			// wire radius
	Variable dphi = asin(R/v)			// angle between wire center and tangent point of wire (from ccd pixel)
	Variable phi0 = atan2(vz,vy)		// angle from yhat to V (V is vector from pixel to wire center)
	Variable tanphi = tan((edge ? phi0-dphi : phi0+dphi))

	// line from pixel to tangent point is:   z = y*tan(phio±dphi) + b
	Variable b = ccdz - ccdy*tanphi

	// line of the beam goes through Xo and with direction given by ki[]
	// line is:  y = z * (kiy/kiz) + beta
	Variable beta1 = Xoy - Xoz*(kiy/kiz)

	depth = (beta1*tanphi+b) / (1-tanphi*kiy/kiz)	// this is the z value of the intercept
	return depth						// this is depth (along z) relative to the Si
End
//
//
// Returns depth (starting point) of ray ending at point xyz (probably a pixel location) that is tangent
// to leading edge of the wire.  The returned depth is relative to the Si position (µm)
//	Function pixelXYZ2depth(geo,ccdy,ccdz,Xw,edge)
//		STRUCT microGeometry &geo
//		Variable ccdy, ccdz					// y, z parts of the xyz location of pixel (beam line coords)
//		Wave Xw							// wire (x,y,z) beam line coordinates in wire units (µm)
//		Variable edge						// 1=leading edge (usual),  0=trailing edge
//		Variable depth=NaN					// computed depth relative to the Si position (µm)
//	
//		Variable Xwy,Xwz					// wire center relative to the Si position (using beam line coords)
//		Xwy = Xw[1] - geo.wire.xyzSi[1]
//		Xwz = Xw[2] - geo.wire.xyzSi[2]
//	
//		ccdy -= geo.wire.xyzSi[1]			// change point on ccd to be relative to the Si position
//		ccdz -= geo.wire.xyzSi[2]
//	
//		Variable kiy,kiz						// components of ki[] in beam line coords (not Wenge coords)
//		kiy = geo.ki[2]
//		kiz = -geo.ki[1]
//	
//		Variable Xoy, Xoz					// point where wire intersects incident beam relative to Si position (µm)
//		Xoy = geo.wire.xyz0[1] - geo.wire.xyzSi[1]
//		Xoz = geo.wire.xyz0[2] - geo.wire.xyzSi[2]
//	
//		Variable r = (geo.wire.dia)/2		// radius of wire (µm)
//	
//		// want to find intersction of two lines
//		//	1) line going through Xo[3] along the ki[3] direction
//		//	2) line going through ccd[3] and tangent to wire at Xw[3] postion
//	
//		Variable d							// distance from ccd to wire center in yz plane
//		d = sqrt((Xwy-ccdy)^2 + (Xwz-ccdz)^2)	// distance from ccd to wire center in yz plane
//		Variable phi = asin(r/d)			// angle between wire center and tangent point (from ccd pixel)
//	
//		Variable vy,vz
//		vy = Xwy-ccdy						// (vy,vz) vector from pixel on ccd toward the wire center
//		vz = Xwz-ccdz
//	
//		// calculate theta, the angle of the correction vector, v[] (whose length is r)
//		// angle measured from  y^ (toward z^) and then points from ccd[] toward the tangent point
//		// atan2(vz,vy) is angle from ccd to the wire center, then correct to reach tangent point
//		// leading edge atan2()-(phi+90¡),  trailing atan2+(phi+90¡)
//		Variable theta =  atan2(vz,vy) + (edge ? -phi-PI/2 : phi+PI/2)
//		vy += r*cos(theta)					// and add the correction
//		vz += r*sin(theta)
//		// v now points from ccd[] toward the leading tangent point
//	
//		// so solve the following vector equation, line from ccd intersects incident beam.
//		//  the intersection is ccd+a*v == Xo+b*ki,   where  c[] = Xo[]-ccd[]
//		//
//		// to solve for a & b we need two equations, so use the y-z plane, giving two equations in a and b
//		//			a*vy -b*kiy = cy		// the y and z components of above vector equation
//		//			a*vz - b*kiz = cz
//		//
//		// multiply second eqn by vy/vz and subtract to eliminate the a term, and solve for b
//		//			b * [ kiz*vy/vz - kiy] = cz*vy/vz - py
//		//
//		Variable cy=Xoy-ccdy, cz=Xoz-ccdz, b
//		b = (cz*vy/vz-cy) / (kiz*vy/vz-kiy)
//		depth = Xoz+b*kiz					// from equation for the incident beam:  xyz = Xo[]+b*ki[]
//		return depth						// relative to the Si position (µm)
//	End









//	Function test_convert()
//		STRUCT microGeometry geo
//		// first set the geometry structure
//		geo.NxCCD=2084	; geo.NyCCD=2084		// # of un-binned pixels in CCD, also called xdim, ydim
//		geo.dpsx=50010	; geo.dpsy=50168		// size of CCD (µm)
//		geo.xcent=1091.53	; geo.ycent=939.941	// center of CCD in un-binned pixels
//		geo.dd=30017.3							// height of CCD (µm)
//		geo.xbet=0.31037	; geo.xgam=0.75037	// incident beam tilt (degrees)
//		geo.xalfd=0.10525	; geo.xbetd=0.171		// CCD tilt (degrees)
//		geo.wire.F = 3830
//		geo.wire.H0 = -4823
//		geo.wire.Hyc = -5000
//		geo.wire.X = 0
//		geo.wire.dia = 51
//		GeometryUpdateCalc(geo)					// this sets rest of geometry, formerly called detectortilt()
//		String/G geoStr
//		StructPut/S/B=2 geo, geoStr				// store the information into geoStr (this stores all of geo.wire too)
//		//
//		//	you can read the struct with:
//		//		SVAR geoStr=geoStr
//		//		StructGet/S/B=2 geo, geoStr		// retrieve the information from geoStr (also gets all of geo.wire.xxx)
//		//
//		printf "using parameters:\r"
//		printf "	NxCCD=%d, NyCCD=%d		// # of un-binned pixels in CCD, also called xdim, ydim\r",geo.NxCCD,geo.NyCCD
//		printf "	dpsx=%g, dpsy=%g			// size of CCD (µm)\r",geo.dpsx,geo.dpsy
//		printf "	xcent=%.2f, ycent=%g	// center of CCD in un-binned pixels\r",geo.xcent,geo.ycent
//		printf "	dd=%g						// height of CCD (µm)\r",geo.dd
//		printf "	xbet=%g, xgam=%g		// incident beam tilt (degrees)\r",geo.xbet,geo.xgam
//		printf "	xalfd=%g, xbetd=%g		// CCD tilt (degrees)\r",geo.xalfd,geo.xbetd
//		printf "		for the wire:\r"
//		printf "	F=%g							// F of the wire for the Si calibration wire scan\r",geo.wire.F
//		printf "	H0=%g							// H of wire centered on the direct beam\r",geo.wire.H0
//		printf "	Hyc=%g						// H of wire centered on the vertical over the Si\r",geo.wire.Hyc
//		printf "	X=%g								// X of wire for Si\r",geo.wire.X
//		printf "	diameter=%g						// diameter of wire \r",geo.wire.dia
//		printf "	Si = {%.2f, %.2f, %.2f}	// wire coordinates of the Si position (µm)\r",geo.wire.xyzSi[0],geo.wire.xyzSi[1],geo.wire.xyzSi[2]
//		printf "	@beam = {%.2f, %.2f, %.2f}	// wire coordinates on the incident beam (µm)\r\r",geo.wire.xyz0[0],geo.wire.xyz0[1],geo.wire.xyz0[2]
//	
//		Make/O/N=(2,2) testXY					// test values
//		testXY= {{1024,1098.4,100,200,300,1042,0},{1024,917.88,1500,1600,350,2040,2040}}
//		Make/O/D testTheta = {42.8988911219471,45.3535835733482,34.9344263723845,32.976490326739,55.5618162886343,24.0805,27.8326}
//	
//		Make/N=3/O/D qhat
//		String errStr
//		Variable/C zpix
//		Variable i,theta, px,py
//		for (i=0;i<DimSize(testXY,0);i+=1)
//			px = testXY[i][0]
//			py = testXY[i][1]
//			theta = pixel2q(geo,px,py,qhat)*180/PI
//			sprintf errStr, "\t\ttheta ERROR = %g",theta-testTheta[i]
//			errStr = SelectString(abs(theta-testTheta[i])>1e-13,"",errStr)
//			printf "theta(%g, %g) = %g¡,		q = {%.5f, %.5f, %.5f}%s\r",px,py,theta,qhat[0],qhat[1],qhat[2],errStr
//			zpix = q2pixel(geo,qhat)				// returns pixel position as a complex number cmplx(px,py)
//			if (abs(	px-real(zpix))+abs(py-imag(zpix)) > 1e-6)		// cannot use 1e-12 with the distortion inversion
//				printf "ERROR, recomputed pixel = (%g, %g)\r",real(zpix),imag(zpix) ; beep
//			endif
//		endfor
//		KillWaves/Z testXY, qhat, testTheta
//	End



//	Function tt1()		// a test for peakcorrection2() and its inverse peakUncorrect()
//	//	Wave xymap = $loadCCDCorrectionTable("")
//	//	printf "loaded distortion table into '%s'\r", GetWavesDataFolder(xymap,2)
//	
//		Wave xymap=root:Packages:geometry:xymap
//	//	Wave xymap = xymap
//	
//		if (!USE_DISTORTION_DEFAULT)
//			Abort "not using distortion, see 'Constant distort =' at top"
//		endif
//	
//	
//		Variable px0,py0,px1,py1
//		px0=200  ;  py0=200
//		px1=px0  ;  py1=py0
//		peakcorrection2(xymap,px1,py1)
//		printf "x = %g  -->  %g,    y = %g  -->  %g,     Æx = %g,  Æy = %g",px0,px1,py0,py1,px1-px0,py1-py0
//		peakUncorrect(xymap,px1,py1)
//		printf "     and invert to get (%g, %g),    Æ=(%g, %g)\r",px1,py1,px0-px1,py0-py1
//	
//		px0=942.054  ;  py0=1673.23
//		px1=px0  ;  py1=py0
//		peakcorrection2(xymap,px1,py1)
//		printf "x = %g  -->  %g,    y = %g  -->  %g,     Æx = %g,  Æy = %g",px0,px1,py0,py1,px1-px0,py1-py0
//		peakUncorrect(xymap,px1,py1)
//		printf "     and invert to get (%g, %g),    Æ=(%g, %g)\r",px1,py1,px0-px1,py0-py1
//	
//		px0=500  ;  py0=500
//		px1=px0  ;  py1=py0
//		peakcorrection2(xymap,px1,py1)
//		printf "x = %g  -->  %g,    y = %g  -->  %g,     Æx = %g,  Æy = %g",px0,px1,py0,py1,px1-px0,py1-py0
//		peakUncorrect(xymap,px1,py1)
//		printf "     and invert to get (%g, %g),    Æ=(%g, %g)\r",px1,py1,px0-px1,py0-py1
//	
//		px0=1042  ;  py0=1042
//		px1=px0  ;  py1=py0
//		peakcorrection2(xymap,px1,py1)
//		printf "x = %g  -->  %g,    y = %g  -->  %g,     Æx = %g,  Æy = %g",px0,px1,py0,py1,px1-px0,py1-py0
//		peakUncorrect(xymap,px1,py1)
//		printf "     and invert to get (%g, %g),    Æ=(%g, %g)\r",px1,py1,px0-px1,py0-py1
//	
//	
//		px0=93.09*2  ;  py0=143.96*2
//		px1=px0  ;  py1=py0
//		peakcorrection2(xymap,px1,py1)
//		printf "x = %g  -->  %g,    y = %g  -->  %g,     Æx = %g,  Æy = %g",px0,px1,py0,py1,px1-px0,py1-py0
//		peakUncorrect(xymap,px1,py1)
//		printf "     and invert to get (%g, %g),    Æ=(%g, %g)\r",px1,py1,px0-px1,py0-py1
//	End


// convert qvector to pixel position
Function/C q2pixel(geo,qvec)						// returns pixel position as a complex number cmplx(px,py)
	STRUCT microGeometry &geo
	Wave qvec										// qvec need not be normalized
	Variable px,py									// final pixel position, full chip unbinned 0 based pixels

	Wave qhat=root:Packages:geometry:q2pixel_qhat, ki=root:Packages:geometry:q2pixel_ki
	Wave kout=root:Packages:geometry:q2pixel_kout, vecB=root:Packages:geometry:q2pixel_vecB
	Wave mat=root:Packages:geometry:q2pixel_mat
	qhat = qvec
	#ifdef FLIP_QX
	qhat[0] = -qhat[0]
	#endif
	normalize(qhat)

	ki = geo.ki[p]
	normalize(ki)
	Variable qLen = -2*MatrixDot(qhat,ki)			// length of qhat, note (q^ dot -ki) always positive
	kout = qhat*qLen + ki							// ko - ki = q
	normalize(kout)

	//	k k^ = a px^ + b py^ + dd {001}				// need to solve for a & b, (know all unit vectors and dd)
	//  -dd {001} = a px^ + b py^ - k k^			// rewrite this way to solve for a, b, k, (I know dd)
	// vecB = mat x {a,b,k}
	vecB = {0,0,-geo.dd}
	mat[0][0] = geo.rded00
	mat[1][0] = geo.rded01
	mat[2][0] = geo.rded02
	mat[0][1] = geo.rded10
	mat[1][1] = geo.rded11
	mat[2][1] = geo.rded12
	mat[][2] = -kout[p]
//	MatrixOp/O M_x = Inv(mat) x vecB
	MatrixLUD mat
	MatrixLUBkSub M_Lower, M_Upper, W_LUPermutation, vecB
	Wave M_x = M_x
	M_x[0] = CCD_X_REVERSE ? -M_x[0] : M_x[0]
	px = geo.xcent + M_x[0]*geo.NxCCD/geo.dpsx	// this formula is for 1 based pixels
	py = geo.ycent - M_x[1]*geo.NyCCD/geo.dpsy
	px -= 1											// convert from 1 based to 0 based pixels
	py -= 1
	Wave xymap=root:Packages:geometry:xymap
	peakUncorrect(xymap,px,py)	// go from true peak position to the measured peak position (put distortion back in), just invert peakcorrection2()

//	KillWaves/Z M_Lower,M_Upper,W_LUPermutation,M_x, kout,ki,qhat,vecB,mat
	return cmplx(px,py)
End



// convert xp,yp positions on screen into absolute coordinates x,y,z (Wenge coords)
Function pixel2q(geo,px,py,qhat)
	STRUCT microGeometry &geo
	Variable px,py						// pixel position, 0 based, first pixel is (0,0), NOT (1,1)
	Wave qhat							// q-vector in Wenge coordinates

	Wave ki=root:Packages:geometry:pixel2q_ki, kout=root:Packages:geometry:pixel2q_kout
	ki = geo.ki[p]									// incident beam direction
	pixel2XYZ(geo,px,py,kout)						// kout is in direction of pixel in xyz coords
	normalize(kout)

	Variable theta = acos(MatrixDot(kout,ki))/2	// Bragg angle (radians)
	if (WaveExists(qhat))
		qhat = kout-ki								// qhat bisects kout and -ki
		normalize(qhat)
		#ifdef FLIP_QX
		qhat[0] = -qhat[0]
		#endif
	endif
	return theta
End
//
// convert px,py positions on detector into x,y,z (Wenge coords) with origin at Si position (µm)
Function pixel2XYZ(geo,px,py,xyz)
	STRUCT microGeometry &geo
	Variable px,py								// pixel is zero based on input, first pixel is (0,0), NOT (1,1)
	Wave xyz										// absolute position of pixel in Wenge coords (a 3-vector) (µm)

	Wave xymap=root:Packages:geometry:xymap
	peakcorrection2(xymap,px,py)				// convert pixel on CCD to a true distance in pixels (takes zero based pixels)

	Variable xxst,yyst
	xxst = (1+px-geo.xcent)*geo.dpsx/geo.NxCCD	// the "1+" is to make 1 based pixels
	xxst = CCD_X_REVERSE ? -xxst : xxst
	yyst = -(1+py-geo.ycent)*geo.dpsy/geo.NyCCD
	Variable xcc,ycc,ddc
	xcc = geo.rded00*xxst + geo.rded10*yyst
	ycc = geo.rded01*xxst + geo.rded11*yyst
	ddc = geo.rded02*xxst + geo.rded12*yyst	// dd is up
	xyz = {xcc,ycc,ddc+geo.dd}					// position in mm
	xyz *= 1000									// convert to µm
	return 0
End


// DEPRECATED		This routine should probably not be used.   And it is probably not used by anyone else either.
// convert pixel location xyz[3](µm) to pixel values
Function/C xyz2pixel(geo,xyz)						// returns pixel values as a complex number cmplx(px,py)
	STRUCT microGeometry &geo
	Wave xyz										// absolute position of pixel in Wenge coords (a 3-vector) (µm), origin at Si
	Variable px,py									// final pixel position, returned (zero based pixels)

	Wave vecB=root:Packages:geometry:q2pixel_vecB
	Wave mat=root:Packages:geometry:q2pixel_mat

	Variable xcc,ycc,ddc								//	xyz = {xcc,ycc,ddc+geo.dd}
	xcc = xyz[0]									// coordinates with origin at center of ccd, (xc,yc)
	ycc = xyz[1]
	ddc = xyz[2] - geo.dd*1000
	vecB = {xcc,ycc,ddc}
	vecB /= 1000									// convert from µm to mm

	mat[0][0] = geo.rded00
	mat[1][0] = geo.rded01
	mat[2][0] = geo.rded02
	mat[0][1] = geo.rded10
	mat[1][1] = geo.rded11
	mat[2][1] = geo.rded12
	mat[0][2] = geo.rded20
	mat[1][2] = geo.rded21
	mat[2][2] = geo.rded22

	MatrixOp/O vecInCCD = Inv(mat) x vecB
	vecInCCD[0] = CCD_X_REVERSE ? -vecInCCD[0] : vecInCCD[0]
	px = geo.xcent + vecInCCD[0]*geo.NxCCD/geo.dpsx	// here pixels are 1 based
	py = geo.ycent - vecInCCD[1]*geo.NyCCD/geo.dpsy
	KillWaves/Z vecInCCD
	px -= 1											// convert from 1 to 0 based pixels
	py -= 1

	Wave xymap=root:Packages:geometry:xymap
	peakUncorrect(xymap,px,py)	// go from true peak position to the measured peak position (put distortion back in), just invert peakcorrection2()
	return cmplx(px,py)							// return 0 based pixels
End



//correct the peak positon from measured to actual, Wenge Yang 5/19/2003
// This assumes a zero based coordinate on input, origin is (0,0), not (1,1)
Function/C peakcorrection2(xymap,px,py)		// returns cmplx(dx,dy)
	Wave xymap								// distortion map
	Variable &px,&py							// un-binned full frame pixel location on input, changed to distorted value at end

	Variable useDistortion = NumVarOrDefault("root:Packages:geometry:useDistortion",USE_DISTORTION_DEFAULT)
	if (!useDistortion || !WaveExists(xymap))	// do not distort
		return cmplx(0,0)
	endif
	Variable dx, dy								// the calculated corrections
	// here we start the calculation of distortion correction just like Wenge
	px += 1									// convert origin from (0,0) to (1,1)
	py += 1


	// This section added to speed up the distortion correction when correcting every pixel, not just a few fitted peaks
	if (exists("root:Packages:geometry:tempCachedDistortionMap")==1)	// this section is when using precomputed distortion calculations
		Wave distortionMap = root:Packages:geometry:tempCachedDistortionMap
		if (NumberByKey("use",note(distortionMap),"="))
			if (WaveDims(distortionMap)!=3 || DimSize(distortionMap,2)!=2)
				Abort "peakcorrection2(), The wave root:Packages:geometry:tempCachedDistortionMap, has an illegal size"
			endif
			Variable i,j, Ni = DimSize(distortionMap,0), Nj = DimSize(distortionMap,1)
			i = (px - DimOffset(distortionMap,0)) / DimDelta(distortionMap,0)
			j = (py - DimOffset(distortionMap,1)) / DimDelta(distortionMap,1)
			if (i<0 || j<0 || i>=Ni || j>=Nj)
				Abort "peakcorrection2(), 'root:Packages:geometry:tempCachedDistortionMap' has a scaling that does not fit input pixels"
			endif
			dx = distortionMap[i][j][0]
			dy = distortionMap[i][j][1]
			px += dx - 1							// add correction and convert origin from (1,1) to (0,0)
			py += dy - 1
			return cmplx(dx,dy)
		endif
	endif


	Variable nx=DimSize(xymap,0)-1, ny=DimSize(xymap,1)-1
	String noteStr = note(xymap)
	if (!stringmatch(StringByKey("CCDname",noteStr,"="),"White2084x2084"))
		px -= 1									// convert pixels back to 0 based
		py -= 1
		return cmplx(0,0)						// this routine is only for (2084 x 2084) CCD
	endif
	Variable cornerX0,cornerY0,cornerX1,cornerY1
	cornerX0 = NumberByKey("cornerX0",noteStr,"=")-1
	cornerY0 = NumberByKey("cornerY0",noteStr,"=")-1
	cornerX1 = NumberByKey("cornerX1",noteStr,"=")-1
	cornerY1 = NumberByKey("cornerY1",noteStr,"=")-1
//	printf "cornerX0=%d,  cornerY0=%d,  cornerX1=%d,  cornerY1=%d\r",cornerX0,cornerY0,cornerX1,cornerY1

	Make/O/D/N=4 cornerxy
	Make/N=4/O tempx,tempy
	tempx = {0,nx,0,nx}
	tempy = {0,0,ny,ny}
	cornerxy[0] = xymap[tempx[cornerX0]][tempy[cornerX0]][0]
	cornerxy[1] = xymap[tempx[cornerY0]][tempy[cornerY0]][1]
	cornerxy[2] = xymap[tempx[cornerX1]][tempy[cornerX1]][0]
	cornerxy[3] = xymap[tempx[cornerY1]][tempy[cornerY1]][1]
	KillWaves/Z tempx,tempy
//	printWave(cornerxy)

	Variable x0=xymap[0][0][0], y0=xymap[0][0][1]
	Variable xxstep, xystep, yxstep, yystep
	xxstep = (xymap[nx][0][0]-xymap[0][0][0]) / nx
	xystep = (xymap[nx][0][1]-xymap[0][0][1]) / nx
	yxstep = (xymap[0][ny][0]-xymap[0][0][0]) / ny
	yystep = (xymap[0][ny][1]-xymap[0][0][1]) / ny

	Variable xcent, ycent
//	xcent = (round(((px-x0)*yystep-(py-y0)*yxstep)/(xxstep*yystep-xystep*yxstep))<nx)>0
//	ycent = (round(((px-x0)*xystep-(py-y0)*xxstep)/(xystep*yxstep-yystep*xxstep))<ny)>0
	xcent = round(((px-x0)*yystep-(py-y0)*yxstep)/(xxstep*yystep-xystep*yxstep))
	xcent = limit(xcent,0,nx)
	ycent = round(((px-x0)*xystep-(py-y0)*xxstep)/(xystep*yxstep-yystep*xxstep))
	ycent = limit(ycent,0,ny)

//	printf "px=%g,  py=%g,  xymap[xcent=%d][ycent=%d][0]=%g,  xymap[%d][%d][2]=%g\r",px, py,xcent,ycent,xymap[xcent][ycent][0],xcent,ycent,xymap[xcent][ycent][1]

	dx=0  ;  dy=0								// initialize the correction
	if (px <= cornerxy[0] || px >= cornerxy[2] || py < cornerxy[1] || py >= cornerxy[3])
		dx = xymap[xcent][ycent][2]			// outside measured range, use default
		dy = xymap[xcent][ycent][3]
	else
		Variable xcent2, ycent2
		Variable r1,r2,r3,r4
		if (px >= xymap[xcent][ycent][0] && py >= xymap[xcent][ycent][1])
//			xcent2=(xcent+1)<nx>0
//			ycent2=(ycent+1)<ny>0
			xcent2 = limit(xcent+1,0,nx)
			ycent2 = limit(ycent+1,0,ny)
			r1=sqrt((px-xymap[xcent][ycent][0])^2+(py-xymap[xcent][ycent][1])^2)
			r2=sqrt((px-xymap[xcent][ycent2][0])^2+(py-xymap[xcent][ycent2][1])^2)
			r3=sqrt((px-xymap[xcent2][ycent][0])^2+(py-xymap[xcent2][ycent][1])^2)
			r4=sqrt((px-xymap[xcent2][ycent2][0])^2+(py-xymap[xcent2][ycent2][1])^2)
			dx=(xymap[xcent][ycent][2]/r1+xymap[xcent][ycent2][2]/r2+xymap[xcent2][ycent][2]/r3+xymap[xcent2][ycent2][2]/r4)
			dx=dx/(1/r1+1/r2+1/r3+1/r4)
			dy=(xymap[xcent][ycent][3]/r1+xymap[xcent][ycent2][3]/r2+xymap[xcent2][ycent][3]/r3+xymap[xcent2][ycent2][3]/r4)
			dy=dy/(1/r1+1/r2+1/r3+1/r4)
		endif
		if (px >= xymap[xcent][ycent][0] && py < xymap[xcent][ycent][1])
//			xcent2=(xcent+1)<nx>0
//			ycent2=(ycent-1)<ny>0
			xcent2 = limit(xcent+1,0,nx)
			ycent2 = limit(ycent-1,0,ny)
			r1=sqrt((px-xymap[xcent][ycent][0])^2+(py-xymap[xcent][ycent][1])^2)
			r2=sqrt((px-xymap[xcent][ycent2][0])^2+(py-xymap[xcent][ycent2][1])^2)
			r3=sqrt((px-xymap[xcent2][ycent][0])^2+(py-xymap[xcent2][ycent][1])^2)
			r4=sqrt((px-xymap[xcent2][ycent2][0])^2+(py-xymap[xcent2][ycent2][1])^2)
			dx=(xymap[xcent][ycent][2]/r1+xymap[xcent][ycent2][2]/r2+xymap[xcent2][ycent][2]/r3+xymap[xcent2][ycent2][2]/r4)
			dx=dx/(1/r1+1/r2+1/r3+1/r4)
			dy=(xymap[xcent][ycent][3]/r1+xymap[xcent][ycent2][3]/r2+xymap[xcent2][ycent][3]/r3+xymap[xcent2][ycent2][3]/r4)
			dy=dy/(1/r1+1/r2+1/r3+1/r4)
		endif
		if (px < xymap[xcent][ycent][0] && py >= xymap[xcent][ycent][1])
//			xcent2=(xcent-1)<nx>0
//			ycent2=(ycent+1)<ny>0
			xcent2 = limit(xcent-1,0,nx)
			ycent2 = limit(ycent+1,0,ny)
			r1=sqrt((px-xymap[xcent][ycent][0])^2+(py-xymap[xcent][ycent][1])^2)
			r2=sqrt((px-xymap[xcent][ycent2][0])^2+(py-xymap[xcent][ycent2][1])^2)
			r3=sqrt((px-xymap[xcent2][ycent][0])^2+(py-xymap[xcent2][ycent][1])^2)
			r4=sqrt((px-xymap[xcent2][ycent2][0])^2+(py-xymap[xcent2][ycent2][1])^2)
			dx=(xymap[xcent][ycent][2]/r1+xymap[xcent][ycent2][2]/r2+xymap[xcent2][ycent][2]/r3+xymap[xcent2][ycent2][2]/r4)
			dx=dx/(1/r1+1/r2+1/r3+1/r4)
			dy=(xymap[xcent][ycent][3]/r1+xymap[xcent][ycent2][3]/r2+xymap[xcent2][ycent][3]/r3+xymap[xcent2][ycent2][3]/r4)
			dy=dy/(1/r1+1/r2+1/r3+1/r4)
		endif
		if (px < xymap[xcent][ycent][0] && py < xymap[xcent][ycent][1])
//			xcent2=(xcent-1)<nx>0
//			ycent2=(ycent-1)<ny>0
			xcent2 = limit(xcent-1,0,nx)
			ycent2 = limit(ycent-1,00,ny)
			r1=sqrt((px-xymap[xcent][ycent][0])^2+(py-xymap[xcent][ycent][1])^2)
			r2=sqrt((px-xymap[xcent][ycent2][0])^2+(py-xymap[xcent][ycent2][1])^2)
			r3=sqrt((px-xymap[xcent2][ycent][0])^2+(py-xymap[xcent2][ycent][1])^2)
			r4=sqrt((px-xymap[xcent2][ycent2][0])^2+(py-xymap[xcent2][ycent2][1])^2)
			dx=(xymap[xcent][ycent][2]/r1+xymap[xcent][ycent2][2]/r2+xymap[xcent2][ycent][2]/r3+xymap[xcent2][ycent2][2]/r4)
			dx=dx/(1/r1+1/r2+1/r3+1/r4)
			dy=(xymap[xcent][ycent][3]/r1+xymap[xcent][ycent2][3]/r2+xymap[xcent2][ycent][3]/r3+xymap[xcent2][ycent2][3]/r4)
			dy=dy/(1/r1+1/r2+1/r3+1/r4)
		endif
	endif

	Variable dx0=NumberByKey("dxCenter",noteStr,"="), dy0=NumberByKey("dyCenter",noteStr,"=")
	dx -= numtype(dx0) ? 0 : dx0										// offset to make (dx,dy) zero at (xc,yc)
	dy -= numtype(dy0) ? 0 : dy0
	px += dx - 1															// add correction and convert origin from (1,1) to (0,0)
	py += dy - 1

	KillWaves/Z cornerxy
	return cmplx(dx,dy)
End
//
//	//correct the peak positon from measured to actual, Wenge Yang 5/19/2003
//	// This assumes a zero based coordinate on input, origin is (0,0), not (1,1)
//	Function/C peakcorrection2(xymap,px,py)		// returns cmplx(dx,dy)
//		Wave xymap								// distortion map
//		Variable &px,&py							// un-binned full frame pixel location on input, changed to distorted value at end
//	
//		Variable useDistortion = NumVarOrDefault("root:Packages:geometry:useDistortion",USE_DISTORTION_DEFAULT)
//		if (!useDistortion || !WaveExists(xymap))	// do not distort
//			return cmplx(0,0)
//		endif
//		px += 1									// convert origin from (0,0) to (1,1)
//		py += 1
//	
//		Variable nx=DimSize(xymap,0)-1, ny=DimSize(xymap,1)-1
//		String noteStr = note(xymap)
//		if (!stringmatch(StringByKey("CCDname",noteStr,"="),"White2084x2084"))
//			return cmplx(0,0)						// this routine is only for (2084 x 2084) CCD
//		endif
//		Variable cornerX0,cornerY0,cornerX1,cornerY1
//		cornerX0 = NumberByKey("cornerX0",noteStr,"=")-1
//		cornerY0 = NumberByKey("cornerY0",noteStr,"=")-1
//		cornerX1 = NumberByKey("cornerX1",noteStr,"=")-1
//		cornerY1 = NumberByKey("cornerY1",noteStr,"=")-1
//	//	printf "cornerX0=%d,  cornerY0=%d,  cornerX1=%d,  cornerY1=%d\r",cornerX0,cornerY0,cornerX1,cornerY1
//	
//		Make/O/D/N=4 cornerxy
//		Make/N=4/O tempx,tempy
//		tempx = {0,nx,0,nx}
//		tempy = {0,0,ny,ny}
//		cornerxy[0] = xymap[tempx[cornerX0]][tempy[cornerX0]][0]
//		cornerxy[1] = xymap[tempx[cornerY0]][tempy[cornerY0]][1]
//		cornerxy[2] = xymap[tempx[cornerX1]][tempy[cornerX1]][0]
//		cornerxy[3] = xymap[tempx[cornerY1]][tempy[cornerY1]][1]
//		KillWaves/Z tempx,tempy
//	//	printWave(cornerxy)
//	
//		Variable x0=xymap[0][0][0], y0=xymap[0][0][1]
//		Variable xxstep, xystep, yxstep, yystep
//		xxstep = (xymap[nx][0][0]-xymap[0][0][0]) / nx
//		xystep = (xymap[nx][0][1]-xymap[0][0][1]) / nx
//		yxstep = (xymap[0][ny][0]-xymap[0][0][0]) / ny
//		yystep = (xymap[0][ny][1]-xymap[0][0][1]) / ny
//	
//		Variable xcent, ycent
//	//	xcent = (round(((px-x0)*yystep-(py-y0)*yxstep)/(xxstep*yystep-xystep*yxstep))<nx)>0
//	//	ycent = (round(((px-x0)*xystep-(py-y0)*xxstep)/(xystep*yxstep-yystep*xxstep))<ny)>0
//		xcent = round(((px-x0)*yystep-(py-y0)*yxstep)/(xxstep*yystep-xystep*yxstep))
//		xcent = limit(xcent,0,nx)
//		ycent = round(((px-x0)*xystep-(py-y0)*xxstep)/(xystep*yxstep-yystep*xxstep))
//		ycent = limit(ycent,0,ny)
//	
//	//	printf "px=%g,  py=%g,  xymap[xcent=%d][ycent=%d][0]=%g,  xymap[%d][%d][2]=%g\r",px, py,xcent,ycent,xymap[xcent][ycent][0],xcent,ycent,xymap[xcent][ycent][1]
//	
//		Variable dx, dy								// the corrections
//		if (px <= cornerxy[0] || px >= cornerxy[2] || py < cornerxy[1] || py >= cornerxy[3])
//			dx = xymap[xcent][ycent][2]			// outside measured range, use default
//			dy = xymap[xcent][ycent][3]
//		else
//			Variable xcent2, ycent2
//			Variable r1,r2,r3,r4
//			if (px >= xymap[xcent][ycent][0] && py >= xymap[xcent][ycent][1])
//	//			xcent2=(xcent+1)<nx>0
//	//			ycent2=(ycent+1)<ny>0
//				xcent2 = limit(xcent+1,0,nx)
//				ycent2 = limit(ycent+1,0,ny)
//				r1=sqrt((px-xymap[xcent][ycent][0])^2+(py-xymap[xcent][ycent][1])^2)
//				r2=sqrt((px-xymap[xcent][ycent2][0])^2+(py-xymap[xcent][ycent2][1])^2)
//				r3=sqrt((px-xymap[xcent2][ycent][0])^2+(py-xymap[xcent2][ycent][1])^2)
//				r4=sqrt((px-xymap[xcent2][ycent2][0])^2+(py-xymap[xcent2][ycent2][1])^2)
//				dx=(xymap[xcent][ycent][2]/r1+xymap[xcent][ycent2][2]/r2+xymap[xcent2][ycent][2]/r3+xymap[xcent2][ycent2][2]/r4)
//				dx=dx/(1/r1+1/r2+1/r3+1/r4)
//				dy=(xymap[xcent][ycent][3]/r1+xymap[xcent][ycent2][3]/r2+xymap[xcent2][ycent][3]/r3+xymap[xcent2][ycent2][3]/r4)
//				dy=dy/(1/r1+1/r2+1/r3+1/r4)
//			endif
//			if (px >= xymap[xcent][ycent][0] && py < xymap[xcent][ycent][1])
//	//			xcent2=(xcent+1)<nx>0
//	//			ycent2=(ycent-1)<ny>0
//				xcent2 = limit(xcent+1,0,nx)
//				ycent2 = limit(ycent-1,0,ny)
//				r1=sqrt((px-xymap[xcent][ycent][0])^2+(py-xymap[xcent][ycent][1])^2)
//				r2=sqrt((px-xymap[xcent][ycent2][0])^2+(py-xymap[xcent][ycent2][1])^2)
//				r3=sqrt((px-xymap[xcent2][ycent][0])^2+(py-xymap[xcent2][ycent][1])^2)
//				r4=sqrt((px-xymap[xcent2][ycent2][0])^2+(py-xymap[xcent2][ycent2][1])^2)
//				dx=(xymap[xcent][ycent][2]/r1+xymap[xcent][ycent2][2]/r2+xymap[xcent2][ycent][2]/r3+xymap[xcent2][ycent2][2]/r4)
//				dx=dx/(1/r1+1/r2+1/r3+1/r4)
//				dy=(xymap[xcent][ycent][3]/r1+xymap[xcent][ycent2][3]/r2+xymap[xcent2][ycent][3]/r3+xymap[xcent2][ycent2][3]/r4)
//				dy=dy/(1/r1+1/r2+1/r3+1/r4)
//			endif
//			if (px < xymap[xcent][ycent][0] && py >= xymap[xcent][ycent][1])
//	//			xcent2=(xcent-1)<nx>0
//	//			ycent2=(ycent+1)<ny>0
//				xcent2 = limit(xcent-1,0,nx)
//				ycent2 = limit(ycent+1,0,ny)
//				r1=sqrt((px-xymap[xcent][ycent][0])^2+(py-xymap[xcent][ycent][1])^2)
//				r2=sqrt((px-xymap[xcent][ycent2][0])^2+(py-xymap[xcent][ycent2][1])^2)
//				r3=sqrt((px-xymap[xcent2][ycent][0])^2+(py-xymap[xcent2][ycent][1])^2)
//				r4=sqrt((px-xymap[xcent2][ycent2][0])^2+(py-xymap[xcent2][ycent2][1])^2)
//				dx=(xymap[xcent][ycent][2]/r1+xymap[xcent][ycent2][2]/r2+xymap[xcent2][ycent][2]/r3+xymap[xcent2][ycent2][2]/r4)
//				dx=dx/(1/r1+1/r2+1/r3+1/r4)
//				dy=(xymap[xcent][ycent][3]/r1+xymap[xcent][ycent2][3]/r2+xymap[xcent2][ycent][3]/r3+xymap[xcent2][ycent2][3]/r4)
//				dy=dy/(1/r1+1/r2+1/r3+1/r4)
//			endif
//			if (px < xymap[xcent][ycent][0] && py < xymap[xcent][ycent][1])
//	//			xcent2=(xcent-1)<nx>0
//	//			ycent2=(ycent-1)<ny>0
//				xcent2 = limit(xcent-1,0,nx)
//				ycent2 = limit(ycent-1,00,ny)
//				r1=sqrt((px-xymap[xcent][ycent][0])^2+(py-xymap[xcent][ycent][1])^2)
//				r2=sqrt((px-xymap[xcent][ycent2][0])^2+(py-xymap[xcent][ycent2][1])^2)
//				r3=sqrt((px-xymap[xcent2][ycent][0])^2+(py-xymap[xcent2][ycent][1])^2)
//				r4=sqrt((px-xymap[xcent2][ycent2][0])^2+(py-xymap[xcent2][ycent2][1])^2)
//				dx=(xymap[xcent][ycent][2]/r1+xymap[xcent][ycent2][2]/r2+xymap[xcent2][ycent][2]/r3+xymap[xcent2][ycent2][2]/r4)
//				dx=dx/(1/r1+1/r2+1/r3+1/r4)
//				dy=(xymap[xcent][ycent][3]/r1+xymap[xcent][ycent2][3]/r2+xymap[xcent2][ycent][3]/r3+xymap[xcent2][ycent2][3]/r4)
//				dy=dy/(1/r1+1/r2+1/r3+1/r4)
//			endif
//		endif
//		px += dx - 1								// add correction and convert origin from (1,1) to (0,0)
//		py += dy - 1
//		KillWaves/Z cornerxy
//		return cmplx(dx,dy)
//	End
//
Function/C peakUncorrect(xymap,px,py)		// go from true peak to the measured peak position (put distortion back in), just invert peakcorrection2()
	Wave xymap								// distortion map
	Variable &px,&py							// un-binned full frame 0 based pixel with distortion correction, changed to un-distorted value at end

	Variable useDistortion = NumVarOrDefault("root:Packages:geometry:useDistortion",USE_DISTORTION_DEFAULT)
	if (!useDistortion || !WaveExists(xymap))	// do not distort
		return cmplx(px,py)
	endif

	Variable maxIter=6, tol=1e-4				// at most 6 iterations, and a tolerance of 1e-4

	Variable pxd=px, pyd=py					// save the starting point (the required distorted position)
	Variable/C dz
	Variable i=0,err = Inf
	for (i=0; i<maxIter && err>tol; i+=1)		// end when it changes by less than tol or after maxIter iterations
		dz = peakcorrection2(xymap,px,py)		// returns changed px,py
		err = abs(px-pxd)+abs(py-pyd)
		px = pxd - real(dz)
		py = pyd - imag(dz)
	endfor
	return cmplx(px,py)
End



//=======================================================================================
//=======================================================================================



//Constant NxCCD=2084, NyCCD=2084		// # of un-binned pixels in CCD, also called xdim, ydim
//Constant dpsx=50010, dpsy=50168			// size of CCD (µm)
//Constant xcent=1091.53, ycent=939.941	// center of CCD in un-binned pixels
//Constant dd=30017.3						// height of CCD (µm)
//Constant xbet=0.31037, xgam=0.75037	// incident beam tilt (degrees)
//Constant xalfd=0.10525, xbetd=0.171		// CCD tilt (degrees)

// output from test(), use these to check the calculation
// test()
//	using parameters:
//		NxCCD=2084, NyCCD=2084		// # of un-binned pixels in CCD, also called xdim, ydim
//		dpsx=50010, dpsy=50168			// size of CCD (µm)
//		xcent=1091.53, ycent=939.941	// center of CCD in un-binned pixels
//		dd=30017.3						// height of CCD (µm)
//		xbet=0.31037, xgam=0.75037		// incident beam tilt (degrees)
//		xalfd=0.10525, xbetd=0.171		// CCD tilt (degrees)
//
//	theta(1024, 1024) = 42.8989¡
//	theta(1098.4, 917.88) = 45.3536¡
//	theta(100, 1500) = 34.9344¡
//	theta(200, 1600) = 32.9765¡
//	theta(300, 350) = 55.5618¡








// calculate rotation matrix from detector tilts angles
// note: in IDL matrices are defined as transpose of the usual matrices
// 
// usese geo to set ki[3], rded[3][3], and rde[3][3]
Function GeometryUpdateCalc(geo)			// formerly callled detectortilt()
	STRUCT microGeometry &geo

	Init_microGeo()
	// first convert all angles from degrees --> radians
	Variable xbr = -geo.xbet*PI/180, xcr=-geo.xgam*PI/180
	Variable xard=geo.xalfd*PI/180, xbrd=geo.xbetd*PI/180

	// rotation matrix for incident beam, use it to set ki[]
	//
	// rot(xbr) =	1		0			0
	//				0	  cos(xbr)	sin(xbr)
	//				0	-sin(xbr)	cos(xbr)
	// rot(xbr) rotates a vector about x negative. So +xbetd move the beam upward as it should.
	//
	// rot(xcr) =	  cos(xcr)	sin(xcr)	0
	//				-sin(xcr)	cos(xcr)	0
	//				  	0			0		1
	// rot(xcr) rotates a vector about z (which points upward) in a negative direction.  So +xgam is backward.
	//
	// rde = rot(xbr) x root(xcr)
	//
	Wave rde = root:Packages:geometry:GeoUpdate_rde
	rde[0][0] = cos(xcr)							// rotation matrix for incident beam
	rde[1][0] = -sin(xcr)*cos(xbr)
	rde[2][0] = sin(xbr)*sin(xcr)
	rde[0][1] = sin(xcr)
	rde[1][1] = cos(xbr)*cos(xcr)
	rde[2][1] = -cos(xcr)*sin(xbr)
	rde[0][2] = 0
	rde[1][2] = sin(xbr)
	rde[2][2] = cos(xbr)
	Wave ki = root:Packages:geometry:GeoUpdate_ki
//	ki = -rde[1][p]									// ki = (0,-1,0) x rde
	ki = {0,-1,0}
	MatrixMultiply ki/T, rde
	Wave M_product=M_product
	ki = M_product
	KillWaves M_product
	normalize(ki)
	geo.ki[0] = ki[0]								// in Wenge coordinates (X=out door, Y=upstream, Z=up)
	geo.ki[1] = ki[1]								// note, cannot automatically iterate a vector assignment when it is in a structure
	geo.ki[2] = ki[2]

	geo.rde00 =   rde[0][0]							// rotation matrix that goes from Wenge system to IDEAL Beam-Line
	geo.rde10 =   rde[2][0]							// this is like rde[][], but includes the swap of Y=-z and Z=y
	geo.rde20 = -rde[1][0]							// NOTE, this is a rotation matrix, so its Transpose is its inverse
	geo.rde01 =   rde[0][1]
	geo.rde11 =   rde[2][1]
	geo.rde21 = -rde[1][1]
	geo.rde02 =   rde[0][2]
	geo.rde12 =   rde[2][2]
	geo.rde22 = -rde[1][2]
	//				// an example of how to use geo.rdexx
	//			qBL[0] = qW[0]*geo.rde00 + qW[1]*geo.rde01 + qW[2]*geo.rde02	// transform Wenge to Ideal Beam-line
	//			qBL[1] = qW[0]*geo.rde10 + qW[1]*geo.rde11 + qW[2]*geo.rde12	//   with incident beam along {0,0,1}
	//			qBL[2] = qW[0]*geo.rde20 + qW[1]*geo.rde21 + qW[2]*geo.rde22

	// rotation matrix for detector
	//
	// rot(a) =  cos(a)		0	  sin(a)			// rotation about y axis
	//			  	0			1		1
	//			-sin(a)		0	  cos(a)
	//
	// rot(b) =	1		0			0				// rotation about x axis
	//				0	  cos(xbr)	sin(xbr)
	//				0	-sin(xbr)	cos(xbr)
	//
	// rded = rot(xard) x root(xbrd)
	//
	geo.rded00 = cos(xard)
	geo.rded10 = 0
	geo.rded20 = -sin(xard)
	geo.rded01 = -sin(xard)*sin(xbrd)
	geo.rded11 = cos(xbrd)
	geo.rded21 = -sin(xbrd)*cos(xard)
	geo.rded02 = sin(xard)*cos(xbrd)
	geo.rded12 = sin(xbrd)
	geo.rded22 = cos(xard)*cos(xbrd)

	Variable len = sqrt(geo.wire.axis[0]*geo.wire.axis[0] + geo.wire.axis[1]*geo.wire.axis[1] + geo.wire.axis[2]*geo.wire.axis[2])
	if (numtype(len)==0)
		geo.wire.axis[0] /= len  ;  geo.wire.axis[1] /= len  ;  geo.wire.axis[2] /= len
	else
		geo.wire.axis[0] = 1  ;  geo.wire.axis[1] = 0  ;  geo.wire.axis[2] = 0
	endif

	// compute wire coordinates that would put wire center at the Si position (µm)
	Variable ir2=1/sqrt(2), dZ
	Variable H0=geo.wire.H0, Hyc=geo.wire.Hyc, F=geo.wire.F

	geo.wire.xyz0[0] = geo.wire.X						// wire coordinates when wire is on the incident beam, computed (µm)
	geo.wire.xyz0[1] = (H0-F)*ir2
	geo.wire.xyz0[2] = (H0+F)*ir2

	geo.wire.xyzSi[2] = (Hyc+F)*ir2					// Z of the Si position in wire units & beam line coordinates (µm)
	dZ = (Hyc-H0)*ir2 								// dist (along z) from where wire intersects beam to Si
	geo.wire.xyzSi[1] = (H0-F)*ir2 - dZ*ki[2]/ki[1]	// Y of the Si in wire units & beam line coordinates (µm)
	geo.wire.xyzSi[0] = geo.wire.X - dZ*ki[0]/ki[1]	// X of the Si in wire units & beam line coordinates (µm)

	Wave xymap=root:Packages:geometry:xymap
	resetCenterDistortionMap(xymap,geo.xcent,geo.ycent)	// update distortion map so that correction at (xc,yc) goes to zero
End
//
Static Function resetCenterDistortionMap(xymap,xc,yc)	// set distortion correction of xc & yc into note of xymap, used so correction to (xc,yc) is zero
	Wave xymap												// distortion map
	Variable xc,yc												// values that should get zero distortion correction, should be geo.xcent,geo.ycent
	if (!WaveExists(xymap))								// wave does not exist, do nothing
		return 1
	endif

	String wnote = note(xymap)
	wnote = ReplaceNumberByKey("dxCenter",wnote,0,"=")			// set correction to zero
	wnote = ReplaceNumberByKey("dyCenter",wnote,0,"=")
	Note/K xymap, wnote													// re-set value in wave note
	if (numtype(xc+yc))
		return 1															// invalid (xc,yc) passed so leave corrections at zero
	endif

	Variable px=xc-1, py=yc-1
	Variable/C dp = peakcorrection2(xymap,px,py)					// returns cmplx(dx,dy), px,py are changed
	wnote = ReplaceNumberByKey("dxCenter",wnote,real(dp),"=")	// add this to corrected to get proper value
	wnote = ReplaceNumberByKey("dyCenter",wnote,imag(dp),"=")
	Note/K xymap, wnote													// store values in wave note
End



//Function MakeGeometryParametersPanel(strStruct)
//	String strStruct									// optional passed value of geo structure, this is used if passed
//	if (strlen(WinList("GeometrySet",";","WIN:64")))	// window alreay exits, bring it to front
//		DoWindow/F GeometrySet
//		return 0
//	endif
//
//	Init_microGeo()
//	NewDataFolder/O root:Packages:geometry:PanelValues// ensure that the needed data folders exist
//	// create globals for the panel
//	Variable/G root:Packages:geometry:PanelValues:xalfd
//	Variable/G root:Packages:geometry:PanelValues:xbetd
//	Variable/G root:Packages:geometry:PanelValues:dd
//	Variable/G root:Packages:geometry:PanelValues:ddOffset
//	Variable/G root:Packages:geometry:PanelValues:NxCCD
//	Variable/G root:Packages:geometry:PanelValues:NyCCD
//	Variable/G root:Packages:geometry:PanelValues:dpsx
//	Variable/G root:Packages:geometry:PanelValues:dpsy
//	Variable/G root:Packages:geometry:PanelValues:xcent
//	Variable/G root:Packages:geometry:PanelValues:ycent
//	Variable/G root:Packages:geometry:PanelValues:xbet
//	Variable/G root:Packages:geometry:PanelValues:xgam
//	Variable/G root:Packages:geometry:PanelValues:dia
//	Variable/G root:Packages:geometry:PanelValues:H0
//	Variable/G root:Packages:geometry:PanelValues:Hyc
//	Variable/G root:Packages:geometry:PanelValues:F
//	NVAR xalfd = root:Packages:geometry:PanelValues:xalfd
//	NVAR xbetd = root:Packages:geometry:PanelValues:xbetd
//	NVAR dd = root:Packages:geometry:PanelValues:dd
//	NVAR ddOffset = root:Packages:geometry:PanelValues:ddOffset
//	NVAR NxCCD = root:Packages:geometry:PanelValues:NxCCD
//	NVAR NyCCD = root:Packages:geometry:PanelValues:NyCCD
//	NVAR dpsx = root:Packages:geometry:PanelValues:dpsx
//	NVAR dpsy = root:Packages:geometry:PanelValues:dpsy
//	NVAR xcent = root:Packages:geometry:PanelValues:xcent
//	NVAR ycent = root:Packages:geometry:PanelValues:ycent
//	NVAR xbet = root:Packages:geometry:PanelValues:xbet
//	NVAR xgam = root:Packages:geometry:PanelValues:xgam
//	NVAR dia = root:Packages:geometry:PanelValues:dia
//	NVAR H0 = root:Packages:geometry:PanelValues:H0
//	NVAR Hyc = root:Packages:geometry:PanelValues:Hyc
//	NVAR F = root:Packages:geometry:PanelValues:F
//
//	STRUCT microGeometry geo
//	geo.NxCCD = 0
//	if (strlen(strStruct)>0)
//		StructGet/S/B=2 geo, strStruct					// found structure information, load into geo
//		Variable/G root:Packages:geometry:PanelValues:dirty=1	// for passed structure
//	else
//		FillGeometryStructDefault(geo)					//fill the geometry structure with current values
//		Variable/G root:Packages:geometry:PanelValues:dirty=0
//	endif
//	xalfd = geo.xalfd										// set the panel globals from the struct
//	xbetd = geo.xbetd
//	dd = geo.dd
//	ddOffset = geo.ddOffset
//	NxCCD = geo.NxCCD
//	NyCCD = geo.NyCCD
//	dpsx = geo.dpsx
//	dpsy = geo.dpsy
//	xcent = geo.xcent
//	ycent = geo.ycent
//	xbet = geo.xbet
//	xgam = geo.xgam
//	dia = geo.wire.dia
//	H0 = geo.wire.H0
//	Hyc = geo.wire.Hyc
//	F = geo.wire.F
//
//	NewPanel /K=1 /W=(675,60,895,653)
//	DoWindow/C GeometrySet
//	SetDrawLayer UserBack
//	SetDrawEnv fname= "Lucida Grande",fsize= 16
//	DrawText 30,18,"Detector Position / tilt"
//	SetDrawEnv fname= "Lucida Grande",fsize= 16
//	DrawText 30,162,"CCD size / resolution"
//	SetDrawEnv fname= "Lucida Grande",fsize= 16
//	DrawText 30,263,"Incident beam direction"
//	SetDrawEnv fname= "Lucida Grande",fsize= 16
//	DrawText 30,323,"Scanning wire"
//	SetVariable xalfd,pos={78,20},size={92,15},proc=GeoPanelVarChangedProc,title="xalfd¡"
//	SetVariable xalfd,help={"alpha tilt of detector (degree)"}
//	SetVariable xalfd,limits={-1,1,0},value= root:Packages:geometry:PanelValues:xalfd,bodyWidth= 60
//	SetVariable xbetd,pos={71,40},size={99,15},proc=GeoPanelVarChangedProc,title="xbetad¡"
//	SetVariable xbetd,help={"beta tilt of detector (degree)"}
//	SetVariable xbetd,limits={-1,1,0},value= root:Packages:geometry:PanelValues:xbetd,bodyWidth= 60
//	SetVariable dd,pos={69,60},size={101,15},proc=GeoPanelVarChangedProc,title="dd (mm)"
//	SetVariable dd,help={"height of detector (mm)"}
//	SetVariable dd,limits={0,2000,0},value= root:Packages:geometry:PanelValues:dd,bodyWidth= 60
//	SetVariable ddOffset,pos={38,80},size={132,15},proc=GeoPanelVarChangedProc,title="dd Offset (mm)"
//	SetVariable ddOffset,help={"offset between motor and dd,  dd = CCDy+ddOffset"}
//	SetVariable ddOffset,limits={-100,100,0},value= root:Packages:geometry:PanelValues:ddOffset,bodyWidth= 60
//	SetVariable xcent,pos={48,100},size={122,15},proc=GeoPanelVarChangedProc,title="xcent (pixel)",format="%.3f"
//	SetVariable xcent,help={"central pixel in CCD X, 1 based & unbinned"}
//	SetVariable xcent,limits={0,5000,0},value= root:Packages:geometry:PanelValues:xcent,bodyWidth= 60
//	SetVariable ycent,pos={48,120},size={122,15},proc=GeoPanelVarChangedProc,title="ycent (pixel)",format="%.3f"
//	SetVariable ycent,help={"central pixel in CCD y, 1 based & unbinned"}
//	SetVariable ycent,limits={0,5000,0},value= root:Packages:geometry:PanelValues:ycent,bodyWidth= 60
//	SetVariable NxCCD,pos={76,165},size={94,15},proc=GeoPanelVarChangedProc,title="NxCCD"
//	SetVariable NxCCD,help={"x size of CCD (pixels)"},format="%d"
//	SetVariable NxCCD,limits={1,5000,0},value= root:Packages:geometry:PanelValues:NxCCD,bodyWidth= 60
//	SetVariable NyCCD,pos={76,185},size={94,15},proc=GeoPanelVarChangedProc,title="NyCCD"
//	SetVariable NyCCD,help={"y size of CCD (pixels)"},format="%d"
//	SetVariable NyCCD,limits={1,5000,0},value= root:Packages:geometry:PanelValues:NyCCD,bodyWidth= 60
//	SetVariable dpsx,pos={58,204},size={112,15},proc=GeoPanelVarChangedProc,title="dpsx (mm)"
//	SetVariable dpsx,help={"x size of CCD (mm)"}
//	SetVariable dpsx,limits={0,500,0},value= root:Packages:geometry:PanelValues:dpsx,bodyWidth= 60
//	SetVariable dpsy,pos={58,222},size={112,15},proc=GeoPanelVarChangedProc,title="dpsy (mm)"
//	SetVariable dpsy,help={"y size of CCD (mm)"}
//	SetVariable dpsy,limits={0,500,0},value= root:Packages:geometry:PanelValues:dpsy,bodyWidth= 60
//	SetVariable xbet,pos={76,265},size={94,15},proc=GeoPanelVarChangedProc,title="xbeta¡"
//	SetVariable xbet,help={"beta tilt of incident beam (degree)"}
//	SetVariable xbet,limits={-1,1,0},value= root:Packages:geometry:PanelValues:xbet,bodyWidth= 60
//	SetVariable xgam,pos={77,285},size={93,15},proc=GeoPanelVarChangedProc,title="xgam¡"
//	SetVariable xgam,help={"gamma tilt of incident beam (degree)"}
//	SetVariable xgam,limits={-12,12,0},value= root:Packages:geometry:PanelValues:xgam,bodyWidth= 60
//	SetVariable dia,pos={44,325},size={126,15},proc=GeoPanelVarChangedProc,title="wire dia (µm)"
//	SetVariable dia,help={"diameter of wire (micron)"}
//	SetVariable dia,limits={1,500,0},value= root:Packages:geometry:PanelValues:dia,bodyWidth= 60
//	SetVariable H0,pos={45,345},size={125,15},proc=GeoPanelVarChangedProc,title="wire H0 (µm)"
//	SetVariable H0,help={"H of wire where it cuts incident beam"}
//	SetVariable H0,limits={-50000,50000,0},value= root:Packages:geometry:PanelValues:H0,bodyWidth= 60
//	SetVariable Hyc,pos={40,365},size={130,15},proc=GeoPanelVarChangedProc,title="wire Hyc (µm)"
//	SetVariable Hyc,help={"H of wire where it cuts vetical ray to (xc,yc)"}
//	SetVariable Hyc,limits={-50000,50000,0},value= root:Packages:geometry:PanelValues:Hyc,bodyWidth= 60
//	SetVariable wireF,pos={52,385},size={118,15},proc=GeoPanelVarChangedProc,title="wire F (µm)"
//	SetVariable wireF,help={"F of wire at H0 & Hyc"}
//	SetVariable wireF,limits={-50000,50000,0},value= root:Packages:geometry:PanelValues:F,bodyWidth= 60
//	Button buttonSaveDefault,pos={35,415},size={150,20},proc=GeometryPanelButtonProc,title="Save as Default"
//	Button buttonSaveDefault,help={"saves these values as default in the current data folder"}
//	Button buttonSaveLocal,pos={35,440},size={150,20},proc=GeometryPanelButtonProc,title="Save in Current Fldr"
//	Button buttonSaveLocal,help={"saves these values as default in the current data folder"}
//	Button buttonEPICS,pos={35,465},size={150,20},proc=GeometryPanelButtonProc,title="Load from EPICS"
//	Button buttonEPICS,help={"Fill the values in this panel from EPICS"}
//	Button buttonFillFromFile,pos={35,490},size={150,20},proc=GeometryPanelButtonProc,title="Load from a file"
//	Button buttonFillFromFile,help={"Fill the values in this panel from a file"}
//	Button buttonWriteToFile,pos={35,515},size={150,20},proc=GeometryPanelButtonProc,title="Write to a file"
//	Button buttonWriteToFile,help={"Write values in this panel to a file"}
//	Button buttonPrintGeo,pos={35,540},size={150,20},proc=GeometryPanelButtonProc,title="Print Geo to History"
//	Button buttonPrintGeo,help={"print geometry values to the history"}
//	Button buttonCancel,pos={35,565},size={150,20},proc=GeometryPanelButtonProc,title="Cancel"
//	Button buttonCancel,help={"Close this panel, do not save any changes"}
//	Button buttonSetdd,pos={171,57},size={45,18},proc=GeometryPanelButtonProc,title="image"
//	Button buttonSetdd,help={"set dd using CCDy from top image"}
//	SetWindow kwTopWin,hook(close)=SetGeoPanelHook
//End
Function MakeGeometryParametersPanel(strStruct)
	String strStruct									// optional passed value of geo structure, this is used if passed
	if (strlen(WinList("GeometrySet",";","WIN:64")))	// window alreay exits, bring it to front
		DoWindow/F GeometrySet
		return 0
	endif
	NewPanel /K=1 /W=(675,60,895,653)
	DoWindow/C GeometrySet
	FillGeometryParametersPanel(strStruct,"GeometrySet",0,0)
	SetWindow kwTopWin,hook(close)=SetGeoPanelHook
End
//
Function/T FillGeometryParametersPanel(strStruct,hostWin,left,top)
	String strStruct									// optional passed value of geo structure, this is used if passed
	String hostWin										// name of home window
	Variable left, top									// offsets from the left and top

	Init_microGeo()
	NewDataFolder/O root:Packages:geometry:PanelValues// ensure that the needed data folders exist
	// create globals for the panel
	Variable/G root:Packages:geometry:PanelValues:xalfd
	Variable/G root:Packages:geometry:PanelValues:xbetd
	Variable/G root:Packages:geometry:PanelValues:dd
	Variable/G root:Packages:geometry:PanelValues:ddOffset
	Variable/G root:Packages:geometry:PanelValues:NxCCD
	Variable/G root:Packages:geometry:PanelValues:NyCCD
	Variable/G root:Packages:geometry:PanelValues:dpsx
	Variable/G root:Packages:geometry:PanelValues:dpsy
	Variable/G root:Packages:geometry:PanelValues:xcent
	Variable/G root:Packages:geometry:PanelValues:ycent
	Variable/G root:Packages:geometry:PanelValues:xbet
	Variable/G root:Packages:geometry:PanelValues:xgam
	Variable/G root:Packages:geometry:PanelValues:dia
	Variable/G root:Packages:geometry:PanelValues:H0
	Variable/G root:Packages:geometry:PanelValues:Hyc
	Variable/G root:Packages:geometry:PanelValues:F
	Variable/G root:Packages:geometry:PanelValues:axisX
	Variable/G root:Packages:geometry:PanelValues:axisY
	Variable/G root:Packages:geometry:PanelValues:axisZ

	NVAR xalfd = root:Packages:geometry:PanelValues:xalfd
	NVAR xbetd = root:Packages:geometry:PanelValues:xbetd
	NVAR dd = root:Packages:geometry:PanelValues:dd
	NVAR ddOffset = root:Packages:geometry:PanelValues:ddOffset
	NVAR NxCCD = root:Packages:geometry:PanelValues:NxCCD
	NVAR NyCCD = root:Packages:geometry:PanelValues:NyCCD
	NVAR dpsx = root:Packages:geometry:PanelValues:dpsx
	NVAR dpsy = root:Packages:geometry:PanelValues:dpsy
	NVAR xcent = root:Packages:geometry:PanelValues:xcent
	NVAR ycent = root:Packages:geometry:PanelValues:ycent
	NVAR xbet = root:Packages:geometry:PanelValues:xbet
	NVAR xgam = root:Packages:geometry:PanelValues:xgam
	NVAR dia = root:Packages:geometry:PanelValues:dia
	NVAR H0 = root:Packages:geometry:PanelValues:H0
	NVAR Hyc = root:Packages:geometry:PanelValues:Hyc
	NVAR F = root:Packages:geometry:PanelValues:F
	NVAR axisX = 	root:Packages:geometry:PanelValues:axisX
	NVAR axisY = 	root:Packages:geometry:PanelValues:axisY
	NVAR axisZ = 	root:Packages:geometry:PanelValues:axisZ

	STRUCT microGeometry geo
	geo.NxCCD = 0
	if (strlen(strStruct)>0)
		StructGet/S/B=2 geo, strStruct					// found structure information, load into geo
		Variable/G root:Packages:geometry:PanelValues:dirty=1	// for passed structure
	else
		FillGeometryStructDefault(geo)					//fill the geometry structure with current values
		Variable/G root:Packages:geometry:PanelValues:dirty=0
	endif
	xalfd = geo.xalfd										// set the panel globals from the struct
	xbetd = geo.xbetd
	dd = geo.dd
	ddOffset = geo.ddOffset
	NxCCD = geo.NxCCD
	NyCCD = geo.NyCCD
	dpsx = geo.dpsx
	dpsy = geo.dpsy
	xcent = geo.xcent
	ycent = geo.ycent
	xbet = geo.xbet
	xgam = geo.xgam
	dia = geo.wire.dia
	H0 = geo.wire.H0
	Hyc = geo.wire.Hyc
	F = geo.wire.F
	axisX = geo.wire.axis[0]
	axisY = geo.wire.axis[1]
	axisZ = geo.wire.axis[2]

	SetWindow kwTopWin,userdata(GeoPanelName)=hostWin+"#geoPanel"
	NewPanel/K=1/W=(left,top,left+221,top+603)/HOST=$hostWin
	ModifyPanel frameStyle=0, frameInset=0
	RenameWindow #,geoPanel
	SetDrawLayer UserBack
	SetDrawEnv fname= "Lucida Grande",fsize= 16
	DrawText 30,25,"Detector Position / tilt"
	SetDrawEnv fname= "Lucida Grande",fsize= 16
	DrawText 30,169,"CCD size / resolution"
	SetDrawEnv fname= "Lucida Grande",fsize= 16
	DrawText 30,270,"Incident beam direction"
	SetDrawEnv fname= "Lucida Grande",fsize= 16
	DrawText 30,330,"Scanning wire"
	SetVariable xalfd,pos={78,27},size={92,15},proc=GeoPanelVarChangedProc,title="xalfd¡"
	SetVariable xalfd,help={"alpha tilt of detector (degree)"}
	SetVariable xalfd,limits={-1,1,0},value= root:Packages:geometry:PanelValues:xalfd,bodyWidth= 60
	SetVariable xbetd,pos={71,47},size={99,15},proc=GeoPanelVarChangedProc,title="xbetad¡"
	SetVariable xbetd,help={"beta tilt of detector (degree)"}
	SetVariable xbetd,limits={-1,1,0},value= root:Packages:geometry:PanelValues:xbetd,bodyWidth= 60
	SetVariable dd,pos={69,67},size={101,15},proc=GeoPanelVarChangedProc,title="dd (mm)"
	SetVariable dd,help={"height of detector (mm)"}
	SetVariable dd,limits={0,2000,0},value= root:Packages:geometry:PanelValues:dd,bodyWidth= 60
	SetVariable ddOffset,pos={38,87},size={132,15},proc=GeoPanelVarChangedProc,title="dd Offset (mm)"
	SetVariable ddOffset,help={"offset between motor and dd,  dd = CCDy+ddOffset"}
	SetVariable ddOffset,limits={-100,100,0},value= root:Packages:geometry:PanelValues:ddOffset,bodyWidth= 60
	SetVariable xcent,pos={48,107},size={122,15},proc=GeoPanelVarChangedProc,title="xcent (pixel)",format="%.3f"
	SetVariable xcent,help={"central pixel in CCD X, 1 based & unbinned"}
	SetVariable xcent,limits={0,5000,0},value= root:Packages:geometry:PanelValues:xcent,bodyWidth= 60
	SetVariable ycent,pos={48,127},size={122,15},proc=GeoPanelVarChangedProc,title="ycent (pixel)",format="%.3f"
	SetVariable ycent,help={"central pixel in CCD y, 1 based & unbinned"}
	SetVariable ycent,limits={0,5000,0},value= root:Packages:geometry:PanelValues:ycent,bodyWidth= 60
	SetVariable NxCCD,pos={76,172},size={94,15},proc=GeoPanelVarChangedProc,title="NxCCD"
	SetVariable NxCCD,help={"x size of CCD (pixels)"},format="%d"
	SetVariable NxCCD,limits={1,5000,0},value= root:Packages:geometry:PanelValues:NxCCD,bodyWidth= 60
	SetVariable NyCCD,pos={76,192},size={94,15},proc=GeoPanelVarChangedProc,title="NyCCD"
	SetVariable NyCCD,help={"y size of CCD (pixels)"},format="%d"
	SetVariable NyCCD,limits={1,5000,0},value= root:Packages:geometry:PanelValues:NyCCD,bodyWidth= 60
	SetVariable dpsx,pos={58,211},size={112,15},proc=GeoPanelVarChangedProc,title="dpsx (mm)"
	SetVariable dpsx,help={"x size of CCD (mm)"}
	SetVariable dpsx,limits={0,500,0},value= root:Packages:geometry:PanelValues:dpsx,bodyWidth= 60
	SetVariable dpsy,pos={58,229},size={112,15},proc=GeoPanelVarChangedProc,title="dpsy (mm)"
	SetVariable dpsy,help={"y size of CCD (mm)"}
	SetVariable dpsy,limits={0,500,0},value= root:Packages:geometry:PanelValues:dpsy,bodyWidth= 60
	SetVariable xbet,pos={76,272},size={94,15},proc=GeoPanelVarChangedProc,title="xbeta¡"
	SetVariable xbet,help={"beta tilt of incident beam (degree)"}
	SetVariable xbet,limits={-1,1,0},value= root:Packages:geometry:PanelValues:xbet,bodyWidth= 60
	SetVariable xgam,pos={77,292},size={93,15},proc=GeoPanelVarChangedProc,title="xgam¡"
	SetVariable xgam,help={"gamma tilt of incident beam (degree)"}
	SetVariable xgam,limits={-12,12,0},value= root:Packages:geometry:PanelValues:xgam,bodyWidth= 60
	SetVariable dia,pos={44,332},size={126,15},proc=GeoPanelVarChangedProc,title="wire dia (µm)"
	SetVariable dia,help={"diameter of wire (micron)"}
	SetVariable dia,limits={1,500,0},value= root:Packages:geometry:PanelValues:dia,bodyWidth= 60
	SetVariable H0,pos={45,352},size={125,15},proc=GeoPanelVarChangedProc,title="wire H0 (µm)"
	SetVariable H0,help={"H of wire where it cuts incident beam"}
	SetVariable H0,limits={-50000,50000,0},value= root:Packages:geometry:PanelValues:H0,bodyWidth= 60
	SetVariable Hyc,pos={40,372},size={130,15},proc=GeoPanelVarChangedProc,title="wire Hyc (µm)"
	SetVariable Hyc,help={"H of wire where it cuts vetical ray to (xc,yc)"}
	SetVariable Hyc,limits={-50000,50000,0},value= root:Packages:geometry:PanelValues:Hyc,bodyWidth= 60
	SetVariable wireF,pos={52,392},size={118,15},proc=GeoPanelVarChangedProc,title="wire F (µm)"
	SetVariable wireF,help={"F of wire at H0 & Hyc"}
	SetVariable wireF,limits={-50000,50000,0},value= root:Packages:geometry:PanelValues:F,bodyWidth= 60

//	SetVariable wireAxisX,pos={52,392},size={30,15},proc=GeoPanelVarChangedProc,title="wire F (µm)"
//	SetVariable wireAxisX,help={"direction of wire axis, beam line coords"}
//	SetVariable wireAxisX,limits={1,1,0},value= root:Packages:geometry:PanelValues:axisX,bodyWidth= 60
//	SetVariable wireAxisY,pos={52,392},size={30,15},proc=GeoPanelVarChangedProc
//	SetVariable wireAxisY,help={"direction of wire axis, beam line coords"}
//	SetVariable wireAxisY,limits={1,1,0},value= root:Packages:geometry:PanelValues:axisY,bodyWidth= 60
//	SetVariable wireAxisZ,pos={52,392},size={30,15},proc=GeoPanelVarChangedProc
//	SetVariable wireAxisZ,help={"direction of wire axis, beam line coords"}
//	SetVariable wireAxisZ,limits={1,1,0},value= root:Packages:geometry:PanelValues:axisZ,bodyWidth= 60

	SetVariable wireAxisX,pos={9,407},size={68,15},proc=GeoPanelVarChangedProc,title="axis"
	SetVariable wireAxisX,help={"direction of wire axis, beam line coords"},format="%.6f"
	SetVariable wireAxisX,limits={1,1,0},value= root:Packages:geometry:PanelValues:axisX,bodyWidth= 45
	SetVariable wireAxisY,pos={84,407},size={45,15},proc=GeoPanelVarChangedProc,title=" "
	SetVariable wireAxisY,help={"direction of wire axis, beam line coords"},format="%.6f"
	SetVariable wireAxisY,limits={1,1,0},value= root:Packages:geometry:PanelValues:axisY,bodyWidth= 45
	SetVariable wireAxisZ,pos={138,407},size={45,15},proc=GeoPanelVarChangedProc,title=" "
	SetVariable wireAxisZ,help={"direction of wire axis, beam line coords"},format="%.6f"
	SetVariable wireAxisZ,limits={1,1,0},value= root:Packages:geometry:PanelValues:axisZ,bodyWidth= 45

	Button buttonSaveDefault,pos={35,422},size={150,20},proc=GeometryPanelButtonProc,title="Save as Default"
	Button buttonSaveDefault,help={"saves these values as default in the current data folder"}
	Button buttonSaveLocal,pos={35,447},size={150,20},proc=GeometryPanelButtonProc,title="Save in Current Fldr"
	Button buttonSaveLocal,help={"saves these values as default in the current data folder"}
	Button buttonEPICS,pos={35,472},size={150,20},proc=GeometryPanelButtonProc,title="Load from EPICS"
	Button buttonEPICS,help={"Fill the values in this panel from EPICS"}
//	Button buttonFillFromFile,pos={35,497},size={150,20},proc=GeometryPanelButtonProc,title="Load from a file"
//	Button buttonFillFromFile,help={"Fill the values in this panel from a file"}
	Button buttonFillFromFile,pos={35,497},size={105,20},proc=GeometryPanelButtonProc,title="Load from file"
	Button buttonFillFromFile,help={"Fill the values in this panel from a file"}
	Button buttonFillFromWeb,pos={145,497},size={40,20},proc=GeometryPanelButtonProc,title="web"
	Button buttonFillFromWeb,help={"Fill the values in this panel from the web"}
	Button buttonWriteToFile,pos={35,522},size={150,20},proc=GeometryPanelButtonProc,title="Write to a file"
	Button buttonWriteToFile,help={"Write values in this panel to a file"}
	Button buttonPrintGeo,pos={35,547},size={150,20},proc=GeometryPanelButtonProc,title="Print Geo to History"
	Button buttonPrintGeo,help={"print geometry values to the history"}
	Button buttonCancel,pos={35,572},size={150,20},proc=GeometryPanelButtonProc,title="Cancel"
	Button buttonCancel,help={"Close this panel, do not save any changes"}
	Button buttonSetdd,pos={171,64},size={45,18},proc=GeometryPanelButtonProc,title="image"
	Button buttonSetdd,help={"set dd using CCDy from top image"}
	return "#geoPanel"
End
//
// this is called by all of the varables, and it just sets the dirty flag for the panel
Function GeoPanelVarChangedProc(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName
	NVAR dirty = root:Packages:geometry:PanelValues:dirty
	dirty = 1
End
//
Function SetGeoPanelHook(s)
	STRUCT WMWinHookStruct &s
	if (s.eventCode !=2)							// only process kill events
		return 0
	endif
	Variable dirty = NumVarOrDefault("root:Packages:geometry:PanelValues:dirty",1)
	if (!dirty)
		return 0								// nothing was changed, do not ask any questions
	endif

	// values were changed in the panel and not saved, ask the user what to do.
	Variable action=1
	Prompt action, "choose an action, do not press cancel",Popup,"Save in default data folder;Save in current data folder;just close window, do not save these values"
	do
		DoPrompt/Help="3D-Xray Diffraction[geometry]" "changes were made, but not saved",action
	while(V_flag)
	// action = V_flag ? 3 : action				// if cancel was pushed, do not save, treat it like a "just quit"
	GeometryPanelButtonProc(StringFromList(action-1,"buttonSaveDefault;buttonSaveLocal;buttonCancel"))
End
//
Function GeometryPanelButtonProc(ctrlName) : ButtonControl
	String ctrlName

	STRUCT microGeometry geo
	NVAR xalfd = root:Packages:geometry:PanelValues:xalfd		// connect to the globals
	NVAR xbetd = root:Packages:geometry:PanelValues:xbetd
	NVAR dd = root:Packages:geometry:PanelValues:dd
	NVAR ddOffset = root:Packages:geometry:PanelValues:ddOffset
	NVAR NxCCD = root:Packages:geometry:PanelValues:NxCCD
	NVAR NyCCD = root:Packages:geometry:PanelValues:NyCCD
	NVAR dpsx = root:Packages:geometry:PanelValues:dpsx
	NVAR dpsy = root:Packages:geometry:PanelValues:dpsy
	NVAR xcent = root:Packages:geometry:PanelValues:xcent
	NVAR ycent = root:Packages:geometry:PanelValues:ycent
	NVAR xbet = root:Packages:geometry:PanelValues:xbet
	NVAR xgam = root:Packages:geometry:PanelValues:xgam
	NVAR dia = root:Packages:geometry:PanelValues:dia
	NVAR H0 = root:Packages:geometry:PanelValues:H0
	NVAR Hyc = root:Packages:geometry:PanelValues:Hyc
	NVAR F = root:Packages:geometry:PanelValues:F
	NVAR axisX = root:Packages:geometry:PanelValues:axisX
	NVAR axisY = root:Packages:geometry:PanelValues:axisY
	NVAR axisZ = root:Packages:geometry:PanelValues:axisZ
	NVAR dirty = root:Packages:geometry:PanelValues:dirty

	if (stringmatch(ctrlName,"buttonCancel"))				// close window, do not save numbers
		String win = StringFromLIst(0,WinList("*",";","WIN:64"))
		SetWindow $win,hook(close)= $""
		DoWindow/K $win
		dirty = 0
		return 0
	elseif (stringmatch(ctrlName,"buttonSaveDefault"))
		String/G root:Packages:geometry:geoStructStr=""	// use default location
		SVAR geoStructStr = root:Packages:geometry:geoStructStr
	elseif (stringmatch(ctrlName,"buttonSaveLocal"))
		String/G geoStructStr=""							// save it in local folder
		SVAR geoStructStr = geoStructStr
	elseif (stringmatch(ctrlName,"buttonEPICS"))
		if (FillGeometryStructFromEpics(geo))
			return 1
		endif
		xalfd = geo.xalfd										// set the panel globals from the struct
		xbetd = geo.xbetd
		dd = geo.dd/1000
		ddOffset = geo.ddOffset/1000
		NxCCD = geo.NxCCD
		NyCCD = geo.NyCCD
		dpsx = geo.dpsx
		dpsy = geo.dpsy
		xcent = geo.xcent
		ycent = geo.ycent
		xbet = geo.xbet
		xgam = geo.xgam
		dia = geo.wire.dia
		H0 = geo.wire.H0
		Hyc = geo.wire.Hyc
		F = geo.wire.F
		axisX = geo.wire.axis[0]
		axisY = geo.wire.axis[1]
		axisZ = geo.wire.axis[2]
		dirty = 1
		return 0
	elseif (stringmatch(ctrlName,"buttonFillFromWeb"))
		Variable fnum
		Open/D/M="Select raw, NOT reconstructed SPE image"/P=home/R/T=".SPE" fnum
		if (strlen(S_fileName)<1)
			return 1
		endif
		GetFileFolderInfo/Q/Z S_fileName
		if (V_Flag)
			return 1
		endif
		if (GeoFromWeb(V_modificationDate,geo))
			return 1
		endif
		xalfd = geo.xalfd										// set the panel globals from the struct
		xbetd = geo.xbetd
		dd = geo.dd
		ddOffset = geo.ddOffset
		NxCCD = geo.NxCCD
		NyCCD = geo.NyCCD
		dpsx = geo.dpsx
		dpsy = geo.dpsy
		xcent = geo.xcent
		ycent = geo.ycent
		xbet = geo.xbet
		xgam = geo.xgam
		dia = geo.wire.dia
		H0 = geo.wire.H0
		Hyc = geo.wire.Hyc
		F = geo.wire.F
		axisX = geo.wire.axis[0]
		axisY = geo.wire.axis[1]
		axisZ = geo.wire.axis[2]
		dirty = 1
		return 0
	elseif (stringmatch(ctrlName,"buttonFillFromFile"))
		if (ReadGeoFromKeyFile("","home",geo))
			return 1
		endif
		xalfd = geo.xalfd										// set the panel globals from the struct
		xbetd = geo.xbetd
		dd = geo.dd
		ddOffset = geo.ddOffset
		NxCCD = geo.NxCCD
		NyCCD = geo.NyCCD
		dpsx = geo.dpsx
		dpsy = geo.dpsy
		xcent = geo.xcent
		ycent = geo.ycent
		xbet = geo.xbet
		xgam = geo.xgam
		dia = geo.wire.dia
		H0 = geo.wire.H0
		Hyc = geo.wire.Hyc
		F = geo.wire.F
		axisX = geo.wire.axis[0]
		axisY = geo.wire.axis[1]
		axisZ = geo.wire.axis[2]
		dirty = 1
		return 0
	elseif (stringmatch(ctrlName,"buttonWriteToFile"))
		geo.xbet = xbet										// set all the geo values from global values
		geo.xgam = xgam
		geo.xalfd = xalfd
		geo.xbetd = xbetd
		geo.dd = dd
		geo.ddOffset = ddOffset
		geo.NxCCD = NxCCD
		geo.NyCCD = NyCCD
		geo.dpsx = dpsx
		geo.dpsy = dpsy
		geo.xcent = xcent
		geo.ycent = ycent
		geo.wire.dia = dia
		geo.wire.F = F
		geo.wire.H0 = H0
		geo.wire.Hyc = Hyc
		geo.wire.axis[0] = axisX
		geo.wire.axis[1] = axisY
		geo.wire.axis[2] = axisZ
		WriteGeoToKeyFile("","home",geo,"")
		return 0
	elseif (stringmatch(ctrlName,"buttonPrintGeo"))		// print shown geometry parameters to the history
		geo.xbet = xbet										// set all the geo values from global values
		geo.xgam = xgam
		geo.xalfd = xalfd
		geo.xbetd = xbetd
		geo.dd = dd
		geo.ddOffset = ddOffset
		geo.NxCCD = NxCCD
		geo.NyCCD = NyCCD
		geo.dpsx = dpsx
		geo.dpsy = dpsy
		geo.xcent = xcent
		geo.ycent = ycent
		geo.wire.dia = dia
		geo.wire.F = F
		geo.wire.H0 = H0
		geo.wire.Hyc = Hyc
		geo.wire.axis[0] = axisX
		geo.wire.axis[1] = axisY
		geo.wire.axis[2] = axisZ
		GeometryUpdateCalc(geo)							// calculate other values
		printGeometry(geo)									// print them all out
		return 0
	elseif (stringmatch(ctrlName,"buttonSetdd"))			// set dd using values from top image
		Wave image=ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";")))
		if (!WaveExists(image))
			DoAlert 0,"image cannot be found for top graph"
			return 1
		endif
		Variable CCDy = NumberByKey("CCDy",note(image),"=")
		if (numtype(CCDy))
			DoAlert 0,"unable to read CCD height from image '"+NameOfWave(image)+"'"
			return 1
		endif
		dd = ddOffset+CCDy
		dirty = 1
		return 0
	else
		Abort "''GeometryPanelButtonProc' called from unsupported control"
	endif

	// store the global values from panel into the 'geoStructStr''
	geo.xalfd = xalfd
	geo.xbetd = xbetd
	geo.dd = dd
	geo.ddOffset = ddOffset
	geo.NxCCD = NxCCD
	geo.NyCCD = NyCCD
	geo.dpsx = dpsx
	geo.dpsy = dpsy
	geo.xcent = xcent
	geo.ycent = ycent
	geo.xbet = xbet
	geo.xgam = xgam
	geo.wire.dia = dia
	geo.wire.H0 = H0
	geo.wire.Hyc = Hyc
	geo.wire.F = F
	geo.wire.axis[0] = axisX
	geo.wire.axis[1] = axisY
	geo.wire.axis[2] = axisZ
	GeometryUpdateCalc(geo)
	StructPut/S/B=2 geo, geoStructStr
	dirty = 0
	SavePackagePreferences/FLSH=1 "microGeo","microGeoPrefs",0,geo
End
//
Function FillGeometryStructFromEpics(geo)		//fill the geometry structure from EPICS (uses caget)
	STRUCT microGeometry &geo

	String geoList = "xalfd;xbetd;dd;ddOffset;NxCCD;NyCCD;dpsx;dpsy;xcent;ycent;xbet;xgam;wireDia;wireH0;wireH1;wireF0"
	String item,pvList=""
	Variable i

	// make list of pv names for each item in geoList (items in geo that have EPICS pv's)
	for (i=0,item=StringFromList(0,geoList); strlen(item); i+=1,item=StringFromList(i,geoList))
		pvList += EPICS_PREFIX+item+";"
	endfor

	// use FUNCREF incase this computer does not have the function get_mult_PV()
	String funcSpec = SelectString(exists("get_mult_PV")==6 , "protoEPICS_geo", "get_mult_PV")
	FUNCREF protoEPICS_geo  func=$funcSpec
	String epicsOut = func(pvList)				// get values from EPICS
	if (strlen(epicsOut)<2)
		return 1
	endif
	// turn result of get_mult_PV() into a keyword=value list
	Variable n=strlen(EPICS_PREFIX)
	pvList = ""
	for (i=0,item=StringFromList(0,epicsOut); strlen(item); i+=1,item=StringFromList(i,epicsOut))
		if (strsearch(item,EPICS_PREFIX,0)==0)
			item = item[n,Inf]
		endif
		pvList += item+";"
	endfor
	geo.xbet = NumberByKey("xbet",pvList,"=")				// set all the geo values from keword list
	geo.xgam = NumberByKey("xgam",pvList,"=")
	geo.xalfd = NumberByKey("xalfd",pvList,"=")
	geo.xbetd = NumberByKey("xbetd",pvList,"=")
	geo.dd = NumberByKey("dd",pvList,"=")*1000				// convert mm to µm
	geo.ddOffset = NumberByKey("ddOffset",pvList,"=")*1000	// convert mm to µm
	geo.NxCCD = NumberByKey("NxCCD",pvList,"=")
	geo.NyCCD = NumberByKey("NyCCD",pvList,"=")
	geo.dpsx = NumberByKey("dpsx",pvList,"=")				// convert mm to µm
	geo.dpsy = NumberByKey("dpsy",pvList,"=")				// convert mm to µm
	geo.xcent = NumberByKey("xcent",pvList,"=")
	geo.ycent = NumberByKey("ycent",pvList,"=")

	geo.wire.dia = NumberByKey("wireDia",pvList,"=")
	geo.wire.F = NumberByKey("wireF0",pvList,"=")
	geo.wire.H0 = NumberByKey("wireH0",pvList,"=")
	geo.wire.Hyc = NumberByKey("wireH1",pvList,"=")
	geo.wire.axis[0] = NumberByKey("wireaxisX",pvList,"=")
	geo.wire.axis[1] = NumberByKey("wireaxisY",pvList,"=")
	geo.wire.axis[2] = NumberByKey("wireaxisZ",pvList,"=")
	GeometryUpdateCalc(geo)
	return 0
End
//
Function/S protoEPICS_geo(str)			// proto function for get_mult_PV(), also contains default values for testing
	String str
	DoAlert 0, "EPICS not available, goemetry unchanged"
	return ""
//	if (!stringmatch(StringFromList(0,GetRTStackInfo(0)),"GeometryPanelButtonProc"))
//		print "EPICS not available, Using default geo values"
//	endif
//	DoAlert 1, "EPICS not available, Use default geo values?"
//	if (V_flag!=1)
//		return ""
//	endif
//	str = "yum:geometry:xalfd=0.15605;yum:geometry:xbetd=0.19635;yum:geometry:dd=40.0143;yum:geometry:ddOffset=0.0143;"
//	str += "yum:geometry:NxCCD=2084;yum:geometry:NyCCD=2084;yum:geometry:dpsx=50.01;yum:geometry:dpsy=50.168;"
//	str += "yum:geometry:xcent=1089.7;yum:geometry:ycent=942.901;yum:geometry:xbet=0.30943;yum:geometry:xgam=0.75192;"
//	str += "yum:geometry:wireDia=52;yum:geometry:wireH0=-5305.9;yum:geometry:wireH1=-4889.03;yum:geometry:wireF0=3830; "
//	return str
End
//
//Function testFillGeometryStructFromEpics()
//	STRUCT microGeometry geo
//	FillGeometryStructFromEpics(geo)
//	printGeometry(geo)
//	String/G geoSructStr
//	StructPut/S/B=2 geo,geoSructStr
//End





// write the geometry structure to a keyword file
Function WriteGeoToKeyFile(fileName,path,geo,fileNote)
	String fileName									// full path name to the file
	String path										// name of an Igor path to use
	STRUCT microGeometry &geo					// the structure to fill from the file
	String fileNote									// a note about this file
	GeometryUpdateCalc(geo)						// calculate other values
	Variable f										// file id
	Open/C="R*ch"/M="new geometry parameters file"/P=$path/T="TEXT"/Z=2 f as fileName
	if (V_flag)
		DoAlert 0, "nothing written to file"
		return 1
	endif
	if( (strlen(fileName)+strlen(fileNote))==0 )
		Prompt fileNote, "add a descriptive note to this file?"
		DoPrompt "descriptive note",fileNote
		if (V_flag)
			fileNote=""
		endif
	endif
	fileName = S_fileName
	fprintf f,"$filetype		geometryFile\n"
//	fprintf f,"$geometryFile\n"
	fprintf f, "$dateWritten	%s\n", date()
	fprintf f, "$timeWritten	%s (%g)\n", Secs2Time(DateTime,3,1),date2secs(-1,-1,-1)/3600
	fprintf f, "$EPOCH			%.0f		// seconds from midnight January 1, 1904\n",DateTime
	if (strlen(fileNote))
		fprintf f,"$fileNote		%s\n",ReplaceString("\r",fileNote,"\\r")
	endif
	fprintf f,"$xalfd			%.5f			// CCD tilt (degrees)\n",geo.xalfd
	fprintf f,"$xbetd			%.5f\n",geo.xbetd
	fprintf f,"$dd				%.4f			// height of CCD (µm)\n",geo.dd
	fprintf f,"$ddOffset		%.4f			// dd = CCDy + ddOffset (µm)\n",geo.ddOffset
	fprintf f,"$NxCCD			%d			// number of un-binned pixels in CCD, also called xdim, ydim\n",geo.NxCCD
	fprintf f,"$NyCCD			%d\n",geo.NyCCD
	fprintf f,"$dpsx			%.3f			// size of CCD (mm)\n",geo.dpsx
	fprintf f,"$dpsy			%.3f\n",geo.dpsy
	fprintf f,"$xcent			%.3f		// center of CCD in un-binned 1 based pixels\n",geo.xcent
	fprintf f,"$ycent			%.3f\n",geo.ycent
	fprintf f,"$xbet			%.5f			// incident beam tilt (degrees)\n",geo.xbet
	fprintf f,"$xgam			%.5f\n",geo.xgam
	fprintf f,"$dia			%.2f			// diameter of wire (µm)\n",geo.wire.dia
	fprintf f,"$H0				%.2f		// H of wire centered on the direct beam (µm)\n",geo.wire.H0
	fprintf f,"$Hyc			%.2f		// H of wire centered on the vertical over the Si (µm)\n",geo.wire.Hyc
	fprintf f,"$F				%.2f			// F of the wire for the Si calibration wire scan (µm)\n",geo.wire.F
	fprintf f,"$wireAxis		{%.6f,%.6f,%.6f} // unit vector along wire axis\n",geo.wire.axis[0],geo.wire.axis[1],geo.wire.axis[2]
	Close f
	String str
	sprintf str, "finished wrting geometry to file '%s'",fileName
	DoAlert 0,str
	return 0
End


// read in the geometry structure from a keyword file
Function ReadGeoFromKeyFile(fileName,path,geo)
	String fileName									// full path name to the file
	String path										// name of an Igor path to use
	STRUCT microGeometry &geo					// the structure to fill from the file

	String list = keyStrFromFile(fileName,"geometryFile",path)// read in all of the tagged values into a keyword list
	if (strlen(list)<1)
		return 1
	endif

	Variable xalfd, xbetd, dd, ddOffset, NxCCD, NyCCD, dpsx, dpsy, xcent, ycent, xbet, xgam, dia, H0, Hyc, F,Ax,Ay,Az
	xalfd = NumberByKey("xalfd",list,"=")
	xbetd = NumberByKey("xbetd",list,"=")
	dd = NumberByKey("dd",list,"=")
	ddOffset = NumberByKey("ddOffset",list,"=")
	NxCCD = NumberByKey("NxCCD",list,"=")
	NyCCD = NumberByKey("NyCCD",list,"=")
	dpsx = NumberByKey("dpsx",list,"=")
	dpsy = NumberByKey("dpsy",list,"=")
	xcent = NumberByKey("xcent",list,"=")
	ycent = NumberByKey("ycent",list,"=")
	xbet = NumberByKey("xbet",list,"=")
	xgam = NumberByKey("xgam",list,"=")
	dia = NumberByKey("dia",list,"=")
	H0 = NumberByKey("H0",list,"=")
	Hyc = NumberByKey("Hyc",list,"=")
	F = NumberByKey("F",list,"=")
	sscanf StringByKey("wireAxis",list,"="),"{%g,%g,%g}",Ax,Ay,Az
	if (V_flag!=3)
		Ax=NaN  ;  Ay=NaN  ;  Az=NaN 
	endif

	geo.xbet = numtype(xbet) ? geo.xbet : xbet		// set all the geo if they are present
	geo.xgam = numtype(xgam) ? geo.xgam : xgam
	geo.xalfd = numtype(xalfd) ? geo.xalfd : xalfd
	geo.xbetd = numtype(xbetd) ? geo.xbetd : xbetd
	geo.dd = numtype(dd) ? geo.dd : dd
	geo.ddOffset = numtype(ddOffset) ? geo.ddOffset : ddOffset
	geo.NxCCD = numtype(NxCCD) ? geo.NxCCD : NxCCD
	geo.NyCCD = numtype(NyCCD) ? geo.NyCCD : NyCCD
	geo.dpsx = numtype(dpsx) ? geo.dpsx : dpsx
	geo.dpsy = numtype(dpsy) ? geo.dpsy : dpsy
	geo.xcent = numtype(xcent) ? geo.xcent : xcent
	geo.ycent = numtype(ycent) ? geo.ycent : ycent
	geo.wire.dia = numtype(dia) ? geo.wire.dia : dia
	geo.wire.F = numtype(F) ? geo.wire.F : F
	geo.wire.H0 = numtype(H0) ? geo.wire.H0 : H0
	geo.wire.Hyc = numtype(Hyc) ? geo.wire.Hyc : Hyc
	if (numtype(Ax+Ay+Az)==0)
		geo.wire.axis[0] = Ax
		geo.wire.axis[1] = Ay
		geo.wire.axis[2] = Az
	endif
	GeometryUpdateCalc(geo)						// calculate other values
	printf "Loaded geometry information from '%s'\r",StringByKey("keyStrFileName",list,"=")
	return 0
End
//
//Function testReadGeoFromKeyFile()
//	STRUCT microGeometry geo
//	Variable i = ReadGeoFromKeyFile("geometry.txt",geo)
//	printGeometry(geo)
//End
//
Static Function/S keyStrFromFile(fname,ftype,path)// read in all of the tagged values from a file into a keyword list
	String fname										// full path name to file with tagged geometry values
	String ftype										// the required file identifier, included as a tag (ftype is optional)
	String path										// name of Igor path

	Variable refNum
	Open/M="file containing tagged values"/P=$path/R/Z=2 refNum as fname
	if (strlen(S_fileName)<1 || !refNum)
		return ""
	endif

	Variable i
	String line
	if (strlen(ftype))								// an ftype is present, so check the file
		Variable OK=0
		FReadLine refNum, line
		line[strlen(line)-1,Inf] = " "				// change the trailing term char to a space
		if (strsearch(line,"$filetype",0)==0)		// starts with "$filetype"
			OK = char2num(line[9])<=32			// $filetype ends with whitespace
			i = strsearch(line,"//",0)				// search for comment identifier
			i = (i<0) ? strlen(line) : i				// points to start of comment id (or end of string)
			line = line[0,i-1]						// trim off comment
			line = ReplaceString("\t",line,";")		// change all tabs, commas, and spaces to semi-colons
			line = ReplaceString(",",line,";")
			line = ReplaceString(" ",line,";")
			OK = OK && WhichListItem(ftype,line)>=0
		else
			ftype = "$"+ftype
			OK = (strsearch(line,ftype,0)==0 && char2num(line[strlen(ftype)])<=32)
		endif
		if (!OK)			// if (strsearch(line,ftype,0) || char2num(line[strlen(ftype)])>32)
			printf "The file must start with '%s', but the first line is '%s'\r",ftype,line
			Close refNum
			return ""
		endif
	endif
//	print "loading geometry parameters from ",S_fileName
	String tagName,value,list = ReplaceStringByKey("keyStrFileName","",S_fileName,"=")
	Variable dollar = char2num("$")
	do
		FReadLine refNum, line
		if (char2num(line)==dollar)				// if it starts with a $, probable tag so process it
			i = strsearch(line,"//",0)				// strip off comments
			if (i>=0)
				line = line[0,i-1]
			endif
			for (i=0;char2num(line[i+1])>32;i+=1)// find end of tag, it ends with a space or lower
			endfor
			tagName = line[1,i]
			if (strlen(tagName)<1)
				continue
			endif

			// check if tag is already in list, only take the first instance of a tag, not the last
			if (keyInList(tagName,list,"=",""))			// check if this key is already in list
				continue
			endif

			for (i=i+1;char2num(line[i])<=32;i+=1)	// find first non-white space
			endfor
			value = line[i,Inf]							// value associated with tagName

			for (i=strlen(value)-1;char2num(value[i])<=32 && i>0;i-=1)	// strip off trailing whitespace
			endfor
			value = ReplaceString(";",value[0,i],":")
			list = ReplaceStringByKey(tagName,list,value,"=")
		endif
	while (strlen(line))							// go until EOF
	Close refNum
	return list
End
//
//	$geometryFile
//	$NxCCD		2084					// number of un-binned pixels in CCD, also called xdim, ydim
//	$NyCCD		2084
//	$dpsx		50.01						// size of CCD (mm)
//	$dpsy		50.168
//	$xcent		1089.70					// center of CCD in un-binned pixels
//	$ycent		942.901
//	$dd			40.0143					// height of CCD (µm)
//	$ddOffset	0.0143					// dd = CCDy + ddOffset (µm)
//	$xbet		0.30943					// incident beam tilt (degrees)
//	$xgam		0.75192
//	$xalfd		0.15605					// CCD tilt (degrees)
//	$xbetd		0.19635
//	//	for the wire:
//	$F			3830						// F of the wire for the Si calibration wire scan (µm)
//	$H0		-5305.9					// H of wire centered on the direct beam (µm)
//	$Hyc		-4889.03					// H of wire centered on the vertical over the Si (µm)
//	$X			0							// X of wire for Si (µm)
//	$dia		52							// diameter of wire (µm)
//	$wireAxis	{1, 0, 0.008727}		// direction of wire axis, should be close to {1,0,0}


// with an Igor epoch, retrieve the geo data form the web server
Function GeoFromWeb(epoch,geo)
	Variable epoch					// Igor epoch, also value of V_modificationDate from GetFileFolderInfo
	STRUCT microGeometry &geo

	String ddate, ttime				// ddate="2007-03-16",ttime="14:32:55"
	Variable im,id,iy
	sscanf Secs2Date(epoch,-1),"%d/%d/%d",id,im,iy
	sprintf ddate "%04d-%02d-%02d",iy,im,id
	ttime = Secs2Time(epoch,3)

	String cmd
//	sprintf cmd "do shell script \"curl -s -f -m 20 http://www.uni.aps.anl.gov/34ide/index.php\\\?tag=geo\\\&date=%sT%s\"",ddate,ttime
	sprintf cmd "do shell script \"curl -s -f -m 20 http://%s/index.php\\\?tag=geo\\\&date=%sT%s\"",GeoWebServer,ddate,ttime
	//	printf "cmd = %s\r\r",cmd

	ExecuteScriptText cmd
	// printf "%s", S_value  ;  print " \r\r"
	if (strsearch(S_value,"ERROR",0)==1)
		return 1
	endif

	Variable N=strlen(S_value), dquote=char2num("\"")			// ASCII value of a double-quote
	if (char2num(S_value[0])==dquote)							// remove leading double-quote from Apple Script
		S_value = S_value[1,N-1]
		N -= 1
	endif
	if (char2num(S_value[N-1])==dquote)						// remove trailing double-quote from Apple Script
		S_value = S_value[0,N-2]
		N -= 1
	endif
	String list = keyStrFromBuffer(S_value)

	Variable xalfd, xbetd, dd, ddOffset, NxCCD, NyCCD, dpsx, dpsy, xcent, ycent, xbet, xgam, dia, H0, Hyc, F,Ax,Ay,Az
	xalfd = NumberByKey("xalfd",list,"=")
	xbetd = NumberByKey("xbetd",list,"=")
	dd = NumberByKey("dd",list,"=")
	ddOffset = NumberByKey("ddOffset",list,"=")
	NxCCD = NumberByKey("NxCCD",list,"=")
	NyCCD = NumberByKey("NyCCD",list,"=")
	dpsx = NumberByKey("dpsx",list,"=")
	dpsy = NumberByKey("dpsy",list,"=")
	xcent = NumberByKey("xcent",list,"=")
	ycent = NumberByKey("ycent",list,"=")
	xbet = NumberByKey("xbet",list,"=")
	xgam = NumberByKey("xgam",list,"=")
	dia = NumberByKey("dia",list,"=")
	H0 = NumberByKey("H0",list,"=")
	Hyc = NumberByKey("Hyc",list,"=")
	F = NumberByKey("F",list,"=")
	sscanf StringByKey("wireAxis",list,"="),"{%g,%g,%g}",Ax,Ay,Az
	if (V_flag!=3)
		Ax=NaN  ;  Ay=NaN  ;  Az=NaN 
	endif

	geo.xbet = numtype(xbet) ? geo.xbet : xbet		// set all the geo if they are present
	geo.xgam = numtype(xgam) ? geo.xgam : xgam
	geo.xalfd = numtype(xalfd) ? geo.xalfd : xalfd
	geo.xbetd = numtype(xbetd) ? geo.xbetd : xbetd
	geo.dd = numtype(dd) ? geo.dd : dd
	geo.ddOffset = numtype(ddOffset) ? geo.ddOffset : ddOffset
	geo.NxCCD = numtype(NxCCD) ? geo.NxCCD : NxCCD
	geo.NyCCD = numtype(NyCCD) ? geo.NyCCD : NyCCD
	geo.dpsx = numtype(dpsx) ? geo.dpsx : dpsx
	geo.dpsy = numtype(dpsy) ? geo.dpsy : dpsy
	geo.xcent = numtype(xcent) ? geo.xcent : xcent
	geo.ycent = numtype(ycent) ? geo.ycent : ycent
	geo.wire.dia = numtype(dia) ? geo.wire.dia : dia
	geo.wire.F = numtype(F) ? geo.wire.F : F
	geo.wire.H0 = numtype(H0) ? geo.wire.H0 : H0
	geo.wire.Hyc = numtype(Hyc) ? geo.wire.Hyc : Hyc
	if (numtype(Ax+Ay+Az)==0)
		geo.wire.axis[0] = Ax
		geo.wire.axis[1] = Ay
		geo.wire.axis[2] = Az
	endif
	GeometryUpdateCalc(geo)						// calculate other values
	if (ItemsInList(GetRTStackInfo(0))<3)		// print everything if run from command line
		printf "Loaded geometry information from web on  %s,  %s",Secs2Date(epoch,2),Secs2Time(epoch,1)
		printf "     this geometry created on   %s  %s\r",StringByKey("dateWritten",list,"="),StringByKey("timeWritten",list,"=")
	endif
	return 0
End
//
Static Function/S keyStrFromBuffer(buf)// read in all of the tagged values from a file into a keyword list
	String buf

	buf = ReplaceString("\r",buf,"\n")
	buf = ReplaceString("\t",buf," ")
	buf += "\n"

	Variable i, i0=0,i1
	String line
	String tagName,value,list = ""
	Variable dollar = char2num("$")
	do
		i1 = strsearch(buf, "\n",i0)
		line = buf[i0,i1-1]
		if (char2num(line)==dollar)				// if it starts with a $, probable tag so process it
			i = strsearch(line,"//",0)				// strip off comments
			if (i>=0)
				line = line[0,i-1]
			endif
			for (i=0;char2num(line[i+1])>32;i+=1)// find end of tag, it ends with a space or lower
			endfor
			tagName = line[1,i]
			if (strlen(tagName)<1)
				continue
			endif

			// check if tag is already in list, only take the first instance of a tag, not the last
			if (keyInList(tagName,list,"=",""))			// check if this key is already in list
				continue
			endif

			for (i=i+1;char2num(line[i])<=32;i+=1)	// find first non-white space
			endfor
			value = line[i,Inf]							// value associated with tagName

			for (i=strlen(value)-1;char2num(value[i])<=32 && i>0;i-=1)	// strip off trailing whitespace
			endfor
			value = ReplaceString(";",value[0,i],":")
			list = ReplaceStringByKey(tagName,list,value,"=")
		endif
		i0 = i1 + 1
	while (strlen(line))
	return list
End



Function FillGeometryStructDefault(geo)					//fill the geometry structure with test values
	STRUCT microGeometry &geo								// returns 0 if something set, 0 is nothing done

	String strStruct=StrVarOrDefault(":geoStructStr","")// set to values in current directory
	if (strlen(strStruct)<1)
		strStruct=StrVarOrDefault("root:Packages:geometry:geoStructStr","")	// try the default values
	endif
	if (strlen(strStruct)>1)
		StructGet/S/B=2 geo, strStruct						// found structure information, load into geo
	else
		LoadPackagePreferences/MIS=1 "microGeo","microGeoPrefs",0,geo
		if (V_flag)
			return 1											// did nothing, nothing found
		endif
	endif
	GeometryUpdateCalc(geo)
	return 0
End

//=======================================================================================
//=======================================================================================



Function setUseDistortion(use)
	Variable use
	Init_microGeo()
	NVAR useDistortion=root:Packages:geometry:useDistortion
	Variable before = useDistortion							// value on entry
	if (numtype(use)==0)									//  valid input, set the global
		useDistortion = use
	else														// need to prompt user
		Prompt use, "Use Distortion Correction?", popup, "NO Distortion Correction;Use Distortion Correction"
		use = useDistortion+1								// default value
		DoPrompt/Help="3D-Xray Diffraction[Distortion Correction]" "Distortion Correction", use
		if (!V_flag)
			useDistortion = use - 1
		endif
	endif
	if (useDistortion && exists("root:Packages:geometry:xymap")!=1)
		LoadStandardDistortion()
	endif
	if (ItemsInList(GetRTStackInfo(0))<2)
		if (before == useDistortion)
			print SelectString(useDistortion,"Distortion Correction remains OFF","Continuing to Use Distortion Correction")
		else
			printf "¥¥Distortion Correction changed from %s to %s\r",SelectString(before,"OFF","ON"),SelectString(useDistortion,"OFF","ON")
		endif
	endif
	return useDistortion
End


//	xymap,Nx,Ny,cornerX0,cornerY0,cornerX1,cornerY1
Function/S loadCCDCorrectionTable(mName)
	String mName						// name for xymap file
	String CCDxyMapfile = "CCD_distorMay03_corr.dat"
	mName = SelectString(strlen(mName)>0,"xymap",CleanupName(mName,0))

	if (ItemsInList(GetRTStackInfo(0))<2)
		print "Loading the correction table ..."
	endif
	Variable refNum
	Open/M="distortion table"/P=home/R/T="TEXT"/Z=2 refNum as CCDxyMapfile
	if (V_flag)
		return ""
	endif
	String line
	FReadLine refNum, line

	Variable Nx,Ny,cornerX0,cornerY0,cornerX1,cornerY1
	sscanf line, "%d %d %d %d %d %d",Nx,Ny,cornerX0,cornerY0,cornerX1,cornerY1
	if (V_flag!=6)
		Close refNum
		return ""
	endif

	String noteStr=""
	noteStr = ReplaceNumberByKey("cornerX0",noteStr,cornerX0,"=")
	noteStr = ReplaceNumberByKey("cornerY0",noteStr,cornerY0,"=")
	noteStr = ReplaceNumberByKey("cornerX1",noteStr,cornerX1,"=")
	noteStr = ReplaceNumberByKey("cornerY1",noteStr,cornerY1,"=")
	noteStr = ReplaceStringByKey("CCDname",noteStr,"White2084x2084","=")

	Make/N=(Nx,Ny,4)/O $mName
	Wave xymap = $mName
	Note/K xymap
	Note xymap, noteStr

	Variable i,j,x1,y1,x2,y2
	for (j=0;j<Ny;j+=1)
		for (i=0;i<Ny;i+=1)
			FReadLine refNum, line
			sscanf line, "%g %g %g %g",x1,y1,x2,y2
			if (V_flag!=4)
				Close refNum
				return ""
			endif
			xymap[i][j][0] = x1
			xymap[i][j][1] = y1
			xymap[i][j][2] = x2
			xymap[i][j][3] = y2
		endfor
	endfor
	Close refNum
	return GetWavesDataFolder(xymap,2)
End
//
Function LoadStandardDistortion()
	String str= FunctionPath("LoadStandardDistortion")
	if (strlen(str)<2)
		return 1
	endif
	String fileName = ParseFilePath(1, str, ":", 1, 0)+"xymap.itx"
	LoadWave/O/T fileName
	if (V_flag!=1)
		DoAlert 0, "Unable to load xymap.itx"
		return 1
	elseif (!stringmatch(StringFromList(0,S_waveNames),"xymap"))
		DoAlert 0,"Loaded something other than xyamp, see history"
		return 1
	endif

	// move to root:Packages:geometry if not already there
	if (stringmatch(GetDataFolder(1),"root:Packages:geometry:"))
		return 0
	endif
	Duplicate/O xymap root:Packages:geometry:xymap
	KillWaves xymap
	printf "¥¥loaded standard distortion from '%s' into file root:Packages:geometry:xymap\r",fileName
	return 0
End


Function MakeDistortionMap()
	Variable N=20
	Make/N=(N,N)/O totalDistortion
	Make/N=(N,N)/O/C xyDistortion
	SetScale/I x 0,2083,"pixels", totalDistortion
	SetScale/I y 0,2083,"pixels", totalDistortion
	SetScale d 0,0,"pixels", totalDistortion
	Variable x0=DimOffset(totalDistortion,0), dx=DimDelta(totalDistortion,0)
	Variable y0=DimOffset(totalDistortion,1), dy=DimDelta(totalDistortion,1)

	Wave xymap=root:Packages:geometry:xymap
	Variable px,py, i,j
	for (j=0;j<N;j+=1)
		py = y0 + j*dy
		for (i=0;i<N;i+=1)
			px = x0 + i*dx
			xyDistortion[i][j] = peakcorrection2(xymap,px,py)		// returns cmplx(dx,dy)
		endfor
	endfor
	totalDistortion = cabs(xyDistortion[p][q])

	if (strlen(WinList("GraphDistortionMap", "","WIN:1"))>1)
		DoWindow/F GraphDistortionMap
		SetDrawLayer /K UserFront	
	else
		Display/K=1/W=(168,55,798,615)
		DoWindow/C GraphDistortionMap
		AppendImage totalDistortion
		ModifyImage totalDistortion ctab= {0,7.1,Terrain,1}
		ModifyGraph margin(right)=72,width={Aspect,1},mirror=2
		ModifyGraph lblMargin(left)=7,lblMargin(bottom)=5,axOffset(left)=-2.57143,lblLatPos(left)=-44,lblLatPos(bottom)=47
		SetAxis/A/R left
		ColorScale/N=text0/F=0/S=3/A=RC/X=-13.63/Y=0.40 image=totalDistortion
	endif
	xyDistortion *= 20							// scale up to make lines longer
	for (j=0;j<N;j+=1)
		py = y0 + j*dy
		for (i=0;i<N;i+=1)
			px = x0 + i*dx
			SetDrawLayer UserFront
			SetDrawEnv xcoord=bottom,ycoord=left,arrow=1, arrowlen=7.0,arrowfat=0.3	// draw the arrow
			DrawLine px,py,px+real(xyDistortion[i][j]),py+imag(xyDistortion[i][j])
		endfor
	endfor
	KillWaves/Z xyDistortion
End

//=======================================================================================
//=======================================================================================

//	Make/O root:Packages:micro:X1correctionWave={0.31,0.32,0.2,0.09,-0.1,-0.3,-0.44,-0.47,-0.42,-0.37,-0.32,-0.19,-0.12,-0.07,0.01,0.17,0.35,0.52,0.67,0.73,0.69,0.53,0.32,0.13,-0.01,-0.13,-0.18,-0.18,-0.18,-0.16,-0.26,-0.25,-0.29,-0.38,-0.41,-0.38,-0.33,-0.18,0.01,0.23,0.42}
//	Make/O root:Packages:micro:Y1correctionWave={0.37,0.38,0.35,0.26,0.1,-0.06,-0.22,-0.31,-0.36,-0.35,-0.62,-0.51,-0.39,-0.28,-0.19,-0.03,0.11,0.25,0.34,0.38,0.38,0.34,0.28,0.17,0.04,-0.09,-0.22,-0.33,-0.4,-0.41,-0.39,-0.36,-0.28,-0.19,-0.08,0.06,0.2,0.3,0.4,0.54,0.66}
//	Make/O root:Packages:micro:Z1correctionWave={-0.09,-0.17,-0.21,-0.24,-0.25,-0.22,-0.14,-0.01,0.16,0.28,0.34,0.38,0.37,0.33,0.3,0.27,0.23,0.1,-0.05,-0.22,-0.18,-0.24,-0.28,-0.29,-0.28,-0.25,-0.18,-0.09,0,0.1,0.28,0.34,0.33,0.32,0.32,0.32,0.31,0.26,0.19,0.09,-0.02}
//	SetScale/P x 0,1,"µm", root:Packages:micro:X1correctionWave,root:Packages:micro:Y1correctionWave,root:Packages:micro:Z1correctionWave
//	SetScale d 0,0,"µm", root:Packages:micro:X1correctionWave,root:Packages:micro:Y1correctionWave,root:Packages:micro:Z1correctionWave
//
Function X1correctedKeyence(X1)		// takes PM500 X1 and returns the "real" X1
	Variable X1							// input the PM500 X1 reading
	Wave deltaW=root:Packages:micro:X1correctionWave
	Variable zWidth = DimDelta(deltaW,0)*(numpnts(deltaW)-1)	// width of the correction wave (µm)
	Variable dX1
	dX1 = mod(X1,zWidth)
	dX1 = dX1 < 0 ? dX1+zWidth : dX1
	return deltaW(dX1)+X1				// return corrected X1, this should be the true X1 value
End
//
Function Y1correctedKeyence(Y1)		// takes PM500 Y1 and returns the "real" Y1
	Variable Y1							// input the PM500 Y1 reading
	Variable Yoff=NumVarOrDefault("Y1offKeyence",0)	// a correction added Apr 23, 2009
	Wave deltaW=root:Packages:micro:Y1correctionWave
	Variable zWidth = DimDelta(deltaW,0)*(numpnts(deltaW)-1)	// width of the correction wave (µm)
	Variable dY1
	dY1 = mod(Y1+Yoff,zWidth)
	dY1 = dY1 < 0 ? dY1+zWidth : dY1
	return deltaW(dY1)+Y1				// return corrected Y1, this should be the true Y1 value
End
//
Function Z1correctedKeyence(Z1)		// takes PM500 Z1 and returns the "real" Z1
	Variable Z1							// input the PM500 Z1 reading
	Variable Zoff=NumVarOrDefault("Z1offKeyence",0)	// a correction added Apr 23, 2009
	Wave deltaW=root:Packages:micro:Z1correctionWave
	Variable zWidth = DimDelta(deltaW,0)*(numpnts(deltaW)-1)	// width of the correction wave (µm)
	Variable dZ1
	dZ1 = mod(Z1+Zoff,zWidth)
	dZ1 = dZ1 < 0 ? dZ1+zWidth : dZ1
	return deltaW(dZ1)+Z1				// return corrected Z1, this should be the true Z1 value
End

//=======================================================================================
//=======================================================================================

Static Function/S MakeUnique3Vector(preset)		// make a new 3-vector, and optionally set it to preset[]
	Wave preset							// optional wave, if none call with $""
	String wName = UniqueName("vec3_",1,0)
	Make/N=3/D $wName
	Wave vec=$wName
	if (WaveExists(preset))
		vec = preset
	endif
	return GetWavesDataFolder(vec,2)		// return full path for assingment to a wave
End
//
// Used to make temporary internal waves, a typical use would be:    Wave xxx=$MakeUnique3Vector($"")
// You will probably want to KillWaves/Z xxx, at the end of the routine.
Static Function/S MakeUnique3x3Mat(preset)		// make a new 3x3 matrix, and optionally set it to preset[]
	Wave preset							// optional wave, if none call with $""
	String wName = UniqueName("mat3x3_",1,0)
	Make/N=(3,3)/D $wName
	Wave mat=$wName
	if (WaveExists(preset))
		mat = preset
	endif
	return GetWavesDataFolder(mat,2)		// return full path for assingment to a wave
End

// moved to Utility_JZT,ipf,  Oct 3, 2013
//Function angleVec2Vec(a,b)			// return the angle between two vectors (degree)
//	Wave a,b
//	Variable dot = MatrixDot(a,b) / (norm(a)*norm(b))
//	dot = limit(dot,-1,1)			// ensure that the acos will exist
//	return acos(dot)*180/PI
//End
//
// moved to Utility_JZT,ipf,  Oct 3, 2013
//Function rotationAngleOfMat(rot)				// returns the total rotation angle of a matrix 'rot'
//	Wave rot									// the rotation matrix
//	Variable trace = MatrixTrace(rot)			// trace = 1 + 2*cos(theta)
//	Variable cosine = (trace-1)/2				// cosine of the rotation angle
//	cosine = (cosine>1) ? (2-cosine) : cosine
//	return acos(cosine)*180/PI				// rotation angle in degrees
//End


//=======================================================================================
// ================================== Start of Grain Math Items ==================================

// compute angle and axis of a rotation matrix
// Aug 2, 2007, this was giving the wrong sign for the rotation, so I reversed the "curl" in defn of axis.  JZT
//		changed "axis[0] = rr[1][2] - rr[2][1]"   -->   "axis[0] = rr[2][1] - rr[1][2]"
Function axisOfMatrix(rot,axis)
	// returns total rotation angle, and sets axis to the axis of the total rotation
	Wave rot										// rotation matrix
	Wave axis										// axis of the rotation (angle is returned)

	Variable sumd = notRotationMatrix(rot)	// accept positive sumd that are less than 1e-4
	if (sumd<0 || sumd>1e-4)
		DoAlert 0, "'"+NameOfWave(rot)+"' is not a rotation matrix in axisOfMatrix()"
		axis = NaN
		return NaN
	elseif (0<sumd)								// close enough to a roation matix, but tidy it up first
		Make/N=(3,3)/O/D axisMat__
		axisMat__ = rot
		if (SquareUpMatrix(axisMat__))
			DoAlert 0, "cannot square up '"+NameOfWave(rot)+"' in axisOfMatrix()"
			axis = NaN
			return NaN
		endif
		Wave rr = axisMat__
	else
		Wave rr = rot
	endif

	Variable cosine = (MatrixTrace(rr)-1)/2	// trace = 1 + 2*cos(theta)
	cosine = limit(cosine,-1,1)
	if (cosine<= -1)								// special for 180¡ rotation,
		axis[0] = sqrt((rr[0][0]+1)/2)
		axis[1] = sqrt((rr[1][1]+1)/2)
		axis[2] = sqrt((rr[2][2]+1)/2)		// always assume z positive
		axis[0] = (rr[0][2]+rr[2][0])<0 ? -axis[0] : axis[0]
		axis[1] = (rr[1][2]+rr[2][1])<0 ? -axis[1] : axis[1]
	else											// rotaion < 180¡, usual formula works
		axis[0] = rr[2][1] - rr[1][2]
		axis[1] = rr[0][2] - rr[2][0]
		axis[2] = rr[1][0] - rr[0][1]
		axis /= 2
	endif
	normalize(axis)
	KillWaves /Z axisMat__
	return acos(cosine)*180/PI					// rotation angle in degrees
End
//	else											// rotaion < 180¡, usual formula works
//		axis[0] = rr[1][2] - rr[2][1]
//		axis[1] = rr[2][0] - rr[0][2]
//		axis[2] = rr[0][1] - rr[1][0]
//	endif


Function AngleBetweenRotationVectors(R1,R2)	// returns rotation between orientations defined by Rodriques vectors (degrees)
	Wave R1,R2								// two Rodriques vectors

	Make/N=(3,3)/O/D mat1_JZT, mat2_JZT
	Wave mat1=mat1_JZT, mat2=mat2_JZT
	rotationMatAboutAxis(R1,NaN,mat1)
	rotationMatAboutAxis(R2,NaN,mat2)
	MatrixOp/O temp_angle_11_JZT = Trace(mat2 x mat1^t)
	Variable cosine = (temp_angle_11_JZT[0][0]-1) / 2
	cosine = min(max(-1,cosine),1)			// ensure cosine in range [-1,1]
	return acos(cosine)*180/PI
//	KillWaves/Z temp_angle_11_JZT, mat1_JZT, mat2_JZT
End


// set mat to be a rotation matrix about axis with angle
Function rotationMatAboutAxis(axis,angle,mat)
	Wave axis				// axis about which to rotate (or possibly Rodriques vector)
	Variable angle			// angle to rotate (degrees), assumes axis is true Rodriques vector if angle invalid
	Wave mat				// desired rotation matrix

	Variable len = norm(axis)
	angle = numtype(angle) ? 2*atan(len) : angle*PI/180	// the rotation angle (rad)

	if (angle==0)			// zero angle rotation is just the identity matrix
		mat = (p==q)
		return 0
	endif

	Variable nx=axis[0]/len, ny=axis[1]/len, nz=axis[2]/len
	Variable cosa=cos(angle), sina=sin(angle)
	Variable c1 = 1-cosa
	if (DimSize(mat,0)!=3 && DimSize(mat,1)!=3)
		Abort "in rotationMatAboutAxis(), mat must be (3,3)"
	endif
//	Redimension/N=(3,3) mat
	// from		http://mathworld.wolfram.com/RodriguesRotationFormula.html (I double checked this too.)
	mat[0][0] = cosa+nx*nx*c1;			mat[0][1] =nx*ny*c1-nz*sina;			mat[0][2] = nx*nz*c1+ny*sina;
	mat[1][0] = nx*ny*c1+nz*sina;		mat[1][1] = cosa+ny*ny*c1;			mat[1][2] =ny*nz*c1-nx*sina;
	mat[2][0] = nx*nz*c1-ny*sina;		mat[2][1] = ny*nz*c1+nx*sina;		mat[2][2] = cosa+nz*nz*c1;
	return 0
End


Function notRotationMatrix(mat)			// false if mat is a rotation matrix, i.e. mat^T = mat^-1
	Wave mat
	String wName = UniqueName("rot",1,0)
	MatrixOp/O $wName = Abs(( mat x mat^t ) - Identity(3))
	Wave diff = $wName
	Variable sumd = sum(diff)
	// printWave(diff)
	// print "sum(diff) = ",sumd
	KillWaves/Z $wName

//	Variable returnVal, thresh=2e-10
	Variable returnVal, thresh=2e-5
	if (sumd>thresh)
		returnVal = sumd						// not a rotation matrix
	elseif (MatrixDet(mat)<0)
		returnVal = -1						// an improper rotation matrix
	else
		returnVal = 0							// yes it is a rotation matrix
	endif
	if (ItemsInList(GetRTStackInfo(0))<2)
		printf "'%s'  is %s matrix\r",NameOfWave(mat),SelectString(returnVal ,"an improper rotation","a rotation","NOT a rotation")
		if (sumd>=thresh && sumd<1e-4)
			printf "	but, it is almost a rotation matrix, sum(|diff|) = %g\r",sumd
		endif
	endif
	return returnVal
End


// fix up a matrix so that it is exactly a rotation matrix, not just close to one
Function SquareUpMatrix(mat)
	Wave mat

	Wave vec0=$MakeUnique3Vector($"")
	Wave vec1=$MakeUnique3Vector($"")
	Wave vec2=$MakeUnique3Vector($"")
	Variable err
	vec0 = mat[p][0]
	vec1 = mat[p][1]
	vec2 = mat[p][2]
	Variable norm0=norm(vec0), norm1=norm(vec1), norm2=norm(vec2)

	// Start with the longest vector, and assume it is correct
	if (norm0>=norm1 && norm0>=norm2)	// X is longest
		Cross vec0,vec1						// Z = X x Y
		Wave W_Cross=W_Cross
		vec2 = W_Cross
		Cross vec2,vec0						// Y = Z x X
		vec1 = W_Cross
	elseif (norm1>=norm0 && norm1>=norm2)// Y is longest
		Cross vec1,vec2						// X = Y x Z
		Wave W_Cross=W_Cross
		vec0 = W_Cross
		Cross vec0,vec1						// Z = X x Y
		vec2 = W_Cross
	else											// Z is longest
		Cross vec2,vec0						// Y = Z x X
		Wave W_Cross=W_Cross
		vec1 = W_Cross
		Cross vec1,vec2						// X = Y x Z
		vec0 = W_Cross
	endif

	err = normalize(vec0)
	err += normalize(vec1)
	err += normalize(vec2)
	mat[][0] = vec0[p]
	mat[][1] = vec1[p]
	mat[][2] = vec2[p]

	err = ( notRotationMatrix(mat) != 0 )		// 0 is a rotation matrix
	KillWaves/Z vec0,vec1,vec2,W_Cross
	return numtype(err)
End

// ================================== End of Grain Math Items ===================================
//=======================================================================================


Function Init_microGeo()
	if (exists("HFSToPosix")<3 && stringmatch(igorInfo(2),"Macintosh"))
		DoAlert 0,"HFSToPosix does not exist.  To index anything you must first install 'HFSAndPosix.xop' and re-start Igor.  Follow instructions printed in history"
		print "=============================================================================="
		print "\r		To install: 	HFSAndPosix.xop\r\r"
		print "Navigate to your Igor application folder, it is probably   'Macintosh HD:Applications:Igor Pro Folder:'."
		print "Open the folder 'More Extensions' in your 'Igor Pro Folder', and then open 'Utilities' folder."
		print " inside Utilities, there should be a file named HFSAndPosix.xop."
		print "Make an alias of this xop, and move the alias to your Desktop, or some place you can find it later."
		print "Navigate to your 'Igor Extensions' folder, it is probably 'Macintosh HD:Applications:Igor Pro Folder:Igor Extensions:'."
		print "Drag the alis you just created 'HFSAndPosix.xop alias' into your 'Igor Extensions' folder'."
		print "Finally re-start Igor."
		print "=============================================================================="
	endif
	if (NumVarOrDefault("root:Packages:geometry:geoInited",0))
		return 0
	endif
	InitLatticeSymPackage()
	NewDataFolder/O root:Packages						// ensure Packages exists
	NewDataFolder/O root:Packages:geometry			// ensure geometry exists
	Variable/G root:Packages:geometry:geoInited=0		// flag that initialization was called successsfully
	NVAR geoInited=root:Packages:geometry:geoInited
	NewDataFolder/O root:Packages:micro				// ensure micro exists
	Make/N=3/O/D root:Packages:geometry:pixel2q_ki, root:Packages:geometry:pixel2q_kout
	Make/N=3/O/D root:Packages:geometry:q2pixel_qhat,root:Packages:geometry:q2pixel_ki
	Make/N=3/O/D root:Packages:geometry:q2pixel_kout,root:Packages:geometry:q2pixel_vecB
	Make/N=(3,3)/O/D root:Packages:geometry:q2pixel_mat
	Make/N=3/O/D root:Packages:geometry:pixel2depth_ccd
	Make/N=3/O/D root:Packages:geometry:GeoUpdate_ki
	Make/N=(3,3)/O/D root:Packages:geometry:GeoUpdate_rde

	Make/O root:Packages:micro:X1correctionWave={0.31,0.32,0.2,0.09,-0.1,-0.3,-0.44,-0.47,-0.42,-0.37,-0.32,-0.19,-0.12,-0.07,0.01,0.17,0.35,0.52,0.67,0.73,0.69,0.53,0.32,0.13,-0.01,-0.13,-0.18,-0.18,-0.18,-0.16,-0.26,-0.25,-0.29,-0.38,-0.41,-0.38,-0.33,-0.18,0.01,0.23,0.42}
	Make/O root:Packages:micro:Y1correctionWave={0.37,0.38,0.35,0.26,0.1,-0.06,-0.22,-0.31,-0.36,-0.35,-0.62,-0.51,-0.39,-0.28,-0.19,-0.03,0.11,0.25,0.34,0.38,0.38,0.34,0.28,0.17,0.04,-0.09,-0.22,-0.33,-0.4,-0.41,-0.39,-0.36,-0.28,-0.19,-0.08,0.06,0.2,0.3,0.4,0.54,0.66}
	Make/O root:Packages:micro:Z1correctionWave={-0.09,-0.17,-0.21,-0.24,-0.25,-0.22,-0.14,-0.01,0.16,0.28,0.34,0.38,0.37,0.33,0.3,0.27,0.23,0.1,-0.05,-0.22,-0.18,-0.24,-0.28,-0.29,-0.28,-0.25,-0.18,-0.09,0,0.1,0.28,0.34,0.33,0.32,0.32,0.32,0.31,0.26,0.19,0.09,-0.02}
	SetScale/P x 0,1,"µm", root:Packages:micro:X1correctionWave,root:Packages:micro:Y1correctionWave,root:Packages:micro:Z1correctionWave
	SetScale d 0,0,"µm", root:Packages:micro:X1correctionWave,root:Packages:micro:Y1correctionWave,root:Packages:micro:Z1correctionWave

	if (exists("root:Packages:geometry:useDistortion")!=2)	// create global, it does not exist
		Variable/G root:Packages:geometry:useDistortion=USE_DISTORTION_DEFAULT
	endif
	NVAR useDistortion=root:Packages:geometry:useDistortion
	geoInited = 1
	if (useDistortion && exists("root:Packages:geometry:xymap")!=1)
		geoInited = LoadStandardDistortion() ? 0 : geoInited
	endif
	return 0
End