#pragma rtGlobals=1		// Use modern global access method.
#pragma ModuleName=microGeo
#pragma version = 1.63
#include  "LatticeSym", version>=4.13
//#define MICRO_VERSION_N
//#define MICRO_GEOMETRY_EXISTS

Constant USE_DISTORTION_DEFAULT = 0		// default is TO USE distortion
Constant MAX_Ndetectors = 3					// maximum number of detectors to permitted
StrConstant LaueGoMainMenuName = "LaueGo (micro)"
Strconstant EPICS_PREFIX="34ide:geometryN:"	// prefix for the geometry PVs
Constant FIRST_PIXEL=0						// use for zero-based pixels
StrConstant defaultGeoTagName="geoN"			// default start of name for geometry files

//Static StrConstant GeoWebServer = "www.uni.aps.anl.gov/34ide"
//Static StrConstant GeoWebServer = "top.uni.aps.anl.gov/34ide"
//Static StrConstant GeoWebServer = "www.aps.anl.gov/Sectors/33_34/34ide"
Static StrConstant GeoWebServer = "sector34.xray.aps.anl.gov/34ide"


// All my calculations assume that an image has zero based pixels, both in the image, and in specifying the ROI
//
// to convert from a ROI pixel (zero based) to full unbinned chip use:
//
//	px = startx + px*groupx + (groupx-1)/2	// pixel is zero based here & startx is zero based
//	py = starty + py*groupy + (groupy-1)/2	// groupx=1 is un-binned


//=======================================================================================
// ====================================== Start of Menus ======================================

Menu "Help"
	"LaueGo",DisplayHelpTopic/K=1 "LaueGo"
End
Menu LaueGoMainMenuName, dynamic
//	"-"
	"ÐÐ> LueGo Panel", /Q, MakeMicroPanel(-1)
	"LaueGo Help",DisplayHelpTopic/K=1  "LaueGo"
//	"-"
	"Hwire calculator",PanelHwireMake() 
	"-"
	"Print Current Geo to History",PrintCurrentGeometry()	
	"Set Geometry Parameters...",MakeGeometryParametersPanel("")
	"   Set Default Geo to Reference values",SetDefaultGeo2Reference()
	help={"Show panel used to set (or view) all of the geometry parameters"}
	Submenu DistortionCorrectionMenuItem(-1)
		DistortionCorrectionMenuItem(0) ,/Q, setUseDistortion(0)		// 0 is OFF
		DistortionCorrectionMenuItem(1) ,/Q, setUseDistortion(1)		// 1 is ON
	End
End
//
Function/S DistortionCorrectionMenuItem(m)		// note, "!\022(" causes a menu item to be checked and disabled
	Variable m
	Init_microGeo()
	NVAR useDistortion=root:Packages:geometry:useDistortion
	if (m==-1)										// -1 is special for the menu header
		return "Distortion Correction is "+SelectString(useDistortion,"OFF","ON")
	elseif (m %^ useDistortion)						// only one is true
		return SelectString(useDistortion,"Start","Stop")+" Using Distortion Correction"
	else
		return "!\022("+SelectString(useDistortion,"Not ","")+"Using Distortion Correction"
	endif
End

Function MICRO_VERSION_N_Func()
	return 0
End

// ====================================== End of Menus =======================================
//=======================================================================================



//=======================================================================================
// =================================== Start of micro Panel ====================================

Function MakeMicroPanel(tab)							// makes the main microPanel
	Variable tab
	if (WinType("micropanel")==7)
		DoWindow/F microPanel
		if (!(tab>=0 && tab<=4))
			ControlInfo/W=microPanel tabMicroA
			tab = V_Value
			if (tab<0)
				ControlInfo/W=microPanel tabMicroB
				tab = V_Value
			endif
		endif
	else
//		NewPanel /W=(672,60,991,695)/K=1/N=microPanel
//		TabControl tabMicroA,pos={0,0},size={317,18},proc=microGeo#microGeneralTabProc
		NewPanel /W=(672,60,966,696)/K=1/N=microPanel as "LaueGo"
		TabControl tabMicroA,pos={-10,0},size={305,18},proc=microGeo#microGeneralTabProc
		TabControl tabMicroA,help={"The main micro-beam functions"},userdata(tabNum)="3"
		TabControl tabMicroA,tabLabel(0)="Index",tabLabel(1)="E scans"
		TabControl tabMicroA,tabLabel(2)="3D-arrays",tabLabel(3)="Geo"
		TabControl tabMicroA,tabLabel(4)="Xtal"

//		TabControl tabMicroB,pos={0,20},size={317,18},proc=microGeo#microGeneralTabProc
		TabControl tabMicroB,pos={-10,20},size={305,18},proc=microGeo#microGeneralTabProc
		TabControl tabMicroB,help={"The main micro-beam functions"},userdata(tabNum)="-1"
		TabControl tabMicroB,tabLabel(0)="Details",tabLabel(1)="Laue Sim."
		TabControl tabMicroB,tabLabel(2)="Calibration",tabLabel(3)="Help"
		SetWindow microPanel,userdata(bottom)=num2str(18+20)

		MoveWindowToCorner("microPanel","tr")			// re-position to the top right corner
	endif
	tab = (tab>=0 && tab<=4) ? tab : 3						// default to geometry
	TabControl tabMicroA,value=tab							// set control to tab
	STRUCT WMTabControlAction  tca						// open panel to the geometry tab
	tca.ctrlName="tabMicroA"
	tca.win="microPanel"
	tca.eventCode=2
	tca.tab = tab
	microGeo#microGeneralTabProc(tca)
	SetWindow kwTopWin,hook(closeWin)=microGeo#microPanelHook	// used to save values if window is killed
End
//
Static Function microGeneralTabProc(tca) : TabControl		// changes the panel for each tab when a tab is selected
	STRUCT WMTabControlAction &tca
	switch( tca.eventCode )
		case -1:											// control killed
			return 1
		case 2:												// mouse up
			Variable tab = tca.tab
			break
	endswitch

	ControlInfo/W=$(tca.win) $(tca.ctrlName)
	if (stringmatch(S_Value,"Help"))
		DisplayHelpTopic/K=1  "LaueGo"
		return 0
	endif
	if (NumberByKey("tabNum",GetUserData(tca.win,tca.ctrlName,"tabNum"),"=")==tab)
		return 0
	endif
	String otherTabs = SelectString(stringmatch(tca.ctrlName,"tabMicroA"),"tabMicroA","tabMicroB")
	TabControl $otherTabs, win=$tca.win,userdata(tabNum)="-1"	// ensure other tab control is -1
	TabControl $otherTabs,value=-1						// set control to tab
	String win = GetUserData(tca.win,otherTabs,"panelName")
	if (WinType(win)==7 && strlen(win))
		KillWindow $win
	endif

	TabControl $tca.ctrlName, win=$tca.win,userdata(tabNum)=num2istr(tab)
	win = GetUserData(tca.win,tca.ctrlName,"panelName")
	if (WinType(win)==7 && strlen(win))
		KillWindow $win
	endif

	String tabList
	if (exists("initEnergyWireScans")==6)
		tabList = "Index;EWscan;Arrays3d;Geometry;Lattice"
	else
		tabList = "Index;Escan;Arrays3d;Geometry;Lattice"
	endif
	String fillFunctionList = SelectString(stringmatch(tca.ctrlName,"tabMicroA"),"Detail;LaueSim;Calibration",tabList)

	String funcName = "Fill"+StringFromList(tab,fillFunctionList)+"ParametersPanel"
	if (exists(funcName)!=6)
		funcName = "microGeo#Load"+funcName[4,Inf]
	endif
	if (exists(funcName)==6)
		FUNCREF FillGeometryParametersPanel FillParametersPanel = $funcName
		Variable bottom = str2num(GetUserData("microPanel","","bottom"))
		win = tca.win+FillParametersPanel("","microPanel",45,3+bottom)
	else
		win = ""
	endif
	TabControl $tca.ctrlName, win=$tca.win,userdata(panelName)=win
	SetActiveSubwindow microPanel
	return 0
End
//
Static Function microPanelHook(s)
	STRUCT WMWinHookStruct &s
	if (s.eventCode !=2 && s.eventCode !=14)					// only process kill events (main or sub window)
		return 0
	endif
	// either main or sub-window killed

#if (NumberByKey("IGORFILEVERSION",IgorInfo(3)) == 6.32)
	print "BUG location (Igor version "+StringByKey("IGORFILEVERSION",IgorInfo(3))+") in microPanelHook(),  window =",s.winName
	if (WinType("s.winName")==7)							// this if by-passes the bug, bu I also loose the hookFn functionality
		String hookFn = GetUserData(s.winName,"","KillHookFn")		// get the hook function (if it exists)
		if (exists(hookFn)==6)									// hook exists, execute it
			FUNCREF SetGeoPanelHook  hook=$hookFn
			hook(s)
		endif
	endif
#else
	String hookFn = GetUserData(s.winName,"","KillHookFn")	// get the hook function (if it exists)
	if (exists(hookFn)==6)										// hook exists, execute it
		FUNCREF SetGeoPanelHook  hook=$hookFn
		hook(s)
	endif
#endif
	return 0
End
//
// This is no longer necessary
//
//Static Function/T LoadLatticeParametersPanel(strStruct,hostWin,left,top)
//	String strStruct									// optional passed value of xtal structure, this is used if passed
//	String hostWin										// name of home window
//	Variable left, top									// offsets from the left and top
//
//	SetWindow kwTopWin,userdata(LatticePanelName)=hostWin+"#LatticePanel"
//	NewPanel/K=1/W=(left,top,left+221,top+365)/HOST=$hostWin
//	ModifyPanel frameStyle=0, frameInset=0
//	RenameWindow #,LatticePanel
//	Button buttonInitLattice,pos={33,31},size={150,50},title="Init Lattice\rpackage",proc=microGeo#LoadPackageButtonProc
//	Button buttonInitLattice,help={"Load lattice symmetry procedures"}
//	return "#LatticePanel"
//End
//
Static Function/T LoadIndexParametersPanel(strStruct,hostWin,left,top)
	String strStruct									// optional passed value of xtal structure, this is used if passed
	String hostWin										// name of home window
	Variable left, top									// offsets from the left and top

	SetWindow kwTopWin,userdata(IndexPanelName)=hostWin+"#IndexPanel"
	NewPanel/K=1/W=(left,top,left+221,top+365)/HOST=$hostWin
	ModifyPanel frameStyle=0, frameInset=0
	RenameWindow #,IndexPanel
	Button buttonInitIndex,pos={33,31},size={150,50},title="Init Indexing\rpackage",proc=microGeo#LoadPackageButtonProc
	Button buttonInitIndex,help={"Load Package for Fitting Peaks and Indexing patterns"}
	return "#IndexPanel"
End
//


// OLD  OLD  OLD  OLD  OLD  OLD  OLD  OLD  OLD  OLD  OLD  OLD  OLD  OLD  OLD
Static Function/T LoadEWscanParametersPanel(strStruct,hostWin,left,top)
	String strStruct									// optional passed value of xtal structure, this is used if passed
	String hostWin										// name of home window
	Variable left, top									// offsets from the left and top

	SetWindow kwTopWin,userdata(EWscanPanelName)=hostWin+"#EWscanPanel"
	NewPanel/K=1/W=(left,top,left+221,top+365)/HOST=$hostWin
	ModifyPanel frameStyle=0, frameInset=0
	RenameWindow #,EwscanPanel
	Button buttonInitEWscan,pos={33,31},size={150,50},title="Init E-W scan\rpackage",proc=microGeo#LoadPackageButtonProc
	Button buttonInitEWscan,help={"Load Package for processing E-W scans"}
	return "#EWscanPanel"
End
// NEW  NEW  NEW  NEW  NEW  NEW  NEW  NEW  NEW  NEW  NEW  NEW  NEW  NEW  NEW
Static Function/T LoadEscanParametersPanel(strStruct,hostWin,left,top)
	String strStruct									// optional passed value of xtal structure, this is used if passed
	String hostWin										// name of home window
	Variable left, top									// offsets from the left and top

	SetWindow kwTopWin,userdata(EscanPanelName)=hostWin+"#EscanPanel"
	NewPanel/K=1/W=(left,top,left+221,top+365)/HOST=$hostWin
	ModifyPanel frameStyle=0, frameInset=0
	RenameWindow #,EscanPanel
	Button buttonInitEscan,pos={33,31},size={150,50},title="Init Energy scan\rpackage",proc=microGeo#LoadPackageButtonProc
	Button buttonInitEscan,help={"Load Package for processing Energy scans"}
	return "#EscanPanel"
End


//
Static Function/T LoadArrays3dParametersPanel(strStruct,hostWin,left,top)
	String strStruct									// optional passed value of xtal structure, this is used if passed
	String hostWin										// name of home window
	Variable left, top									// offsets from the left and top

	SetWindow kwTopWin,userdata(Arrays3dPanelName)=hostWin+"#Arrays3dPanel"
	NewPanel/K=1/W=(left,top,left+221,top+365)/HOST=$hostWin
	ModifyPanel frameStyle=0, frameInset=0
	RenameWindow #,Arrays3dPanel
	if (strlen(WinList("ArrayOf3dOrients*", "", "WIN:128")))
		TitleBox titleLoadedArray3d,pos={18,37},size={180,40},title="\\JCmacros already loaded,\rsee the menubar"
		TitleBox titleLoadedArray3d,fSize=16,frame=0
	else
		Button buttonInitArrays3d,pos={33,31},size={150,50},title="Init 3d Arrays\rpackage",proc=microGeo#LoadPackageButtonProc
		Button buttonInitArrays3d,help={"Processes 3D array of matricies into viewable of rotataion angles, dislocation tensors, GND, etc."}
	endif
	return "#Arrays3dPanel"
End
//
Static Function/T LoadDetailParametersPanel(strStruct,hostWin,left,top)
	String strStruct									// optional passed value of xtal structure, this is used if passed
	String hostWin										// name of home window
	Variable left, top									// offsets from the left and top

	SetWindow kwTopWin,userdata(DetailPanelName)=hostWin+"#DetailPanel"
	NewPanel/K=1/W=(left,top,left+221,top+365)/HOST=$hostWin
	ModifyPanel frameStyle=0, frameInset=0
	RenameWindow #,DetailPanel
	Button buttonInitDetail,pos={33,31},size={150,50},title="Init Detail\rpackage",proc=microGeo#LoadPackageButtonProc
	Button buttonInitDetail,help={"Load Package for Details"}
	return "#DetailPanel"
End
//
Static Function/T LoadLaueSimParametersPanel(strStruct,hostWin,left,top)
	String strStruct									// optional passed value of xtal structure, this is used if passed
	String hostWin										// name of home window
	Variable left, top									// offsets from the left and top

	SetWindow kwTopWin,userdata(LaueSimPanelName)=hostWin+"#LaueSimPanel"
	NewPanel/K=1/W=(left,top,left+221,top+365)/HOST=$hostWin
	ModifyPanel frameStyle=0, frameInset=0
	RenameWindow #,LaueSimPanel
	Button buttonInitLaueSim,pos={33,31},size={150,50},title="Init LaueSim\rpackage",proc=microGeo#LoadPackageButtonProc
	Button buttonInitLaueSim,help={"Load Package for Laue Simulation"}
	return "#LaueSimPanel"
End
//
Static Function/T LoadCalibrationParametersPanel(strStruct,hostWin,left,top)
	String strStruct									// optional passed value of xtal structure, this is used if passed
	String hostWin										// name of home window
	Variable left, top									// offsets from the left and top

	SetWindow kwTopWin,userdata(CalibrationPanelName)=hostWin+"#CalibrationPanel"
	NewPanel/K=1/W=(left,top,left+221,top+365)/HOST=$hostWin
	ModifyPanel frameStyle=0, frameInset=0
	RenameWindow #,CalibrationPanel
	Button buttonInitCalibration,pos={33,31},size={150,50},title="Init Calibration\rpackage",proc=microGeo#LoadPackageButtonProc
	Button buttonInitCalibration,help={"Load Package for Calibration"}
	return "#CalibrationPanel"
End
//
Static Function LoadPackageButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	if (ba.eventCode!=2)
		return 0
	endif

	TabControl tabMicroA, win=microPanel,userdata(tabNum)="-1"
	TabControl tabMicroA,value=-1
	TabControl tabMicroB, win=microPanel,userdata(tabNum)="-1"
	TabControl tabMicroB,value=-1
	STRUCT WMTabControlAction  tca							// open panel to the geometry tab
	tca.win = "microPanel"
	tca.eventCode = 2

	String ctrlName
	Variable tab

	if (stringmatch(ba.ctrlName,"buttonInitIndex"))
		Execute/P "INSERTINCLUDE  \"IndexingN\", version>=4.41"
		Execute/P "COMPILEPROCEDURES "
		Execute/P/Q "InitIndexingPackage()"
		tab = 0
		ctrlName = "tabMicroA"
	elseif (stringmatch(ba.ctrlName,"buttonInitIndexLots"))
		Execute/P "DELETEINCLUDE  \"IndexingN\""
		Execute/P "INSERTINCLUDE  \"IndexLotsN\", version>=2.27"
		Execute/P "COMPILEPROCEDURES "
		Execute/P/Q "initSymmetryOperations()"
		tab = 0
		ctrlName = "tabMicroA"


	elseif (stringmatch(ba.ctrlName,"buttonInitEWscan"))	// OLD  OLD  OLD  OLD  OLD
		Execute/P "INSERTINCLUDE  \"EnergyWireScansN\", version>=1.31"
		Execute/P "COMPILEPROCEDURES "
		Execute/P/Q "initEnergyWireScans()"
		tab = 1
		ctrlName = "tabMicroA"
	elseif (stringmatch(ba.ctrlName,"buttonInitEscan"))		// NEW  NEW  NEW  NEW  NEW
		Execute/P "INSERTINCLUDE  \"EnergyScansN\", version>=2.01"
		Execute/P "COMPILEPROCEDURES "
		Execute/P/Q "initEnergyScans()"
		tab = 1
		ctrlName = "tabMicroA"


	elseif (stringmatch(ba.ctrlName,"buttonInitArrays3d"))
		Execute/P "INSERTINCLUDE  \"ArrayOf3dOrientsN\", version>=2.58"
		tab = 2
		ctrlName = "tabMicroA"
	elseif (stringmatch(ba.ctrlName,"buttonInitDetail"))
		Execute/P "INSERTINCLUDE  \"DepthResolvedQueryN\", version>=1.52"
		Execute/P "COMPILEPROCEDURES "
		Execute/P/Q "initDepthResolve()"
		tab = 0
		ctrlName = "tabMicroB"
	elseif (stringmatch(ba.ctrlName,"buttonInitLaueSim"))
		Execute/P "INSERTINCLUDE  \"LaueSimulation\", version>=1.09"
		tab = 1
		ctrlName = "tabMicroB"
	elseif (stringmatch(ba.ctrlName,"buttonInitCalibration"))
		Execute/P "INSERTINCLUDE  \"DetectorCalibration\", version>=0.79"
		tab = 2
		ctrlName = "tabMicroB"
	else
		return 1
	endif

	String cmd
	sprintf cmd,"microGeo#changeTab(\"%s\",%d)",ctrlName,tab
	Execute/P/Q "COMPILEPROCEDURES "
	Execute/P/Q cmd
	return 0
End
//
Static Function changeTab(ctrlName,tab)
	String ctrlName
	Variable tab

	STRUCT WMTabControlAction  tca							// open panel to the geometry tab
	tca.win = "microPanel"
	tca.eventCode = 2
	tca.ctrlName = ctrlName
	tca.tab = tab
	TabControl $ctrlName,value=tab
	microGeo#microGeneralTabProc(tca)
	DoUpdate
End
//
Static Function MoveWindowToCorner(win,corner)
	String win
	String corner													//  one of "TL", "TR", "BL", "BR"
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
		String list = StringByKey("SCREEN1",IgorInfo(0))	// for a Mac, use the whole screen (except menubar)
		Variable i=strsearch(list,"RECT=",0)
		sscanf list[i+5,Inf], "%g,%g,%g,%g", leftScreen,topScreen,rightScreen,bottomScreen
		topScreen += 44											// leave room for menubar on Mac
	endif

	GetWindow $win,wsize										// size of window to move
	Variable width=V_right-V_left, height=V_bottom-V_top	// get size of the window
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


// ==================================== End of micro Panel ====================================
//=======================================================================================



//=======================================================================================
// ====================================== Start of Tests=======================================

Function testBasics(pxIn,pyIn)
	Variable pxIn,pyIn
	Variable px=pxIn, py=pyIn							// save input values

	STRUCT microGeometry g
	GeoReferenceOrientation(g)							// set geometry to reference values
	g.d[0].R[0] += 0.01									// add some imprecision for this example
	g.d[0].R[1] -= 0.02
	g.d[0].R[2] += 0.03
	g.d[1].P[2] = -40e3
	g.d[2].P[2] = -45e3
	GeometryUpdateCalc(g)								// tidy up and do all pre-calculations

	String/G root:Packages:geometry:geoStructStr=""	// use default location
	SVAR geoStructStr = root:Packages:geometry:geoStructStr
	StructPut/S/B=2 g, geoStructStr
	if (pxIn==1024 && pyIn==1024)
		printGeometry(g)
		print "**************************************************************************\r\r"
	endif

	// run the test example
	Make/N=3/O/D xyz
	Variable i
	for (i=0;i<g.Ndetectors;i+=1)
		if (i==1)
			pxIn /= 2;	pyIn/=2						// detectors 1 & 2  are smaller, so reduce pixel position by x2
		endif
		printf "for detector %d\r",i
		pixel2XYZ(g.d[i],pxIn,pyIn,xyz)			// convert pixel position to the beam line coordinate system
		printf "	pixel = (%g, %g) -->  XYZ=(%g, %g, %g)(µm)\r",pxIn,pyIn,xyz[0],xyz[1],xyz[2]

		xyz *= 1.5									// extend point so that it does not lie on detector, I will find the intersection
		XYZ2pixel(g.d[i],xyz,px,py)					// find pixel position where vector xyz will intercept detector
		printf "		Inverting:   --> pixel=(%g, %g),  Æpixel=(%.2g, %.2g)\r",px,py,px-pxIn,py-pyIn
	endfor
	KillWaves/Z xyz
End

// ======================================= End of Tests=======================================
//=======================================================================================



//=======================================================================================
// ===================================== Start of Geometry =====================================

Structure microGeometry							// structure definition
	int16	Ndetectors									// number of detectors in use, must be <= MAX_Ndetectors
	STRUCT detectorGeometry d[MAX_Ndetectors]	// geometry parameters for each detector
	STRUCT wireGeometry wire
	STRUCT sampleGeometry s
EndStructure
//
// Reference position for detector has its face lying in the Z=0 plane exactly centered on the origin, and the X & Y directions on the detector are along the X & Y axes.

// Position of detector: first translate by P, and then rotate detector around R.  Since rho is rotation matrix calculated from R:
//		{X,Y,Z} = rho x [P+{x',y',z'}] ,   where XYZ are beam line coords, and {x'y'z'} are coordinates in detector reference frame.
// Size of detector is measured to the outer edge of the outer most pixels.  So conversion from position to pixel for the x direction is:
//		x' = ( pixel - (Nx-1)/2 )*pitch,   where sizeX = Nx*pitch.  This puts the coordinate of pixel (i,j) at the center of the pixel.
//
Structure detectorGeometry			// structure definition for a detector
	int16 used						// TRUE=detector used, FALSE=detector un-used
	int32 Nx, Ny					// # of un-binned pixels in full detector
	double sizeX,sizeY				// outside size of detector (sizeX = Nx*pitchX), measured to outer edge of outer pixels (micron)
	double R[3]						// rotation vector (length is angle in radians)
	double P[3]						// translation vector (micron)

	uchar timeMeasured[100]		// when this geometry was calculated
	uchar geoNote[100]				// note
	uchar detectorID[100]			// unique detector ID
	uchar distortionMapFile[100]	// name of file with distortion map

	double rho00, rho01, rho02	// rotation matrix internally calculated from R[3]
	double rho10, rho11, rho12
	double rho20, rho21, rho22
EndStructure
//
Structure wireGeometry			// structure definition
	double origin[3]				// raw PM500 coordinates that would put wire center at Origin, (the Si position), (micron)
	double dia						// diameter of wire (micron)
	int16 knife						// true if wire on a knife edge, false for free-standing wire
	double F						// F of the wire for the wire scan in raw PM500 units (not very important) (micron)
	double axis[3]					// unit vector in direction of wire afis in positioner frame, e.g. (0,0,1) is along positioner x-axis
	double axisR[3]					// vector axis rotated by R, now direction of wire in beam line frame (calculated internally)
	double R[3]						// rotation vector for the wire positioner (length is angle in radians)
	double Rmag					// magnitude of | R[3] |, in degrees (computed internally)
	double R00, R01, R02			// rotation matrix from R[3] internally calculated
	double R10, R11, R12
	double R20, R21, R22
EndStructure
//
Structure sampleGeometry
	double O[3]						// raw PM500 coordinates where sample is at origin, (the Si position), (micron)
	double R[3]						// rotation vector for the sample positioner (length is angle in radians)
	double Rmag					// magnitude of | Rs[3] |, in degrees (computed internally)
	double R00, R01, R02			// rotation matrix from R[3] internally calculated
	double R10, R11, R12
	double R20, R21, R22
EndStructure



Function DetectorGeometryLocate()					// This only used by FunctionPath("DetectorGeometryLocate")
	return 0										//   to find the path to this file
End


Function PrintCurrentGeometry()					// prints the current default geometry to the history
	STRUCT microGeometry g
	FillGeometryStructDefault(g)					//fill the geometry structure with current values
	printGeometry(g)
End
//
Function printGeometry(g)							// print the details for passed geometry to the history window
	STRUCT microGeometry &g
	printf "current geomery parameters  (using %d detectors)\r",g.Ndetectors
	if (!SampleBad(g.s))
		printf "Sample\r"
		printf "	Origin = {%g, %g, %g}					// raw PM500 coordinates where sample is at origin (micron)\r",g.s.O[0],g.s.O[1],g.s.O[2]
		printf "	R = {%g, %g, %g}	// (= %g¡) sample positioner rotation vector (radian)\r",g.s.R[0],g.s.R[1],g.s.R[2],g.s.Rmag
		if (NumVarOrDefault("root:Packages:geometry:printVerbose",0))
			printf "			{%+.6f, %+.6f, %+.6f}	// rotation matrix from sample R\r",g.s.R00, g.s.R01, g.s.R02
			printf "	Rs =	{%+.6f, %+.6f, %+.6f}\r",g.s.R10, g.s.R11, g.s.R12
			printf "			{%+.6f, %+.6f, %+.6f}\r",g.s.R20, g.s.R21, g.s.R22
		endif
	endif

	Variable i
	for (i=0;i<MAX_Ndetectors;i+=1)				// info about all of the detectors
		if (g.d[i].used)
			printf "Detector %d\r",i
			printDetector(g.d[i])
		endif
	endfor

	if (WireBad(g.wire))
		printf "Wire UN-Defined *****\r"			// info about the wire
	else
		printf "Wire:\r"								// info about the wire
		printf "	Origin = {%.2f, %.2f, %.2f}	// raw PM500 coordinates to put wire at Origin (Si position) (µm)\r",g.wire.origin[0],g.wire.origin[1],g.wire.origin[2]
		printf "	diameter=%.2f				// diameter of wire (µm)\r",g.wire.dia
		print "\t"+SelectString(g.wire.knife,"free standing wire","wire mounted on a knife edge")
		printf "	wire axis direction = {%.6f, %.6f, %.6f}	// direction of wire axis in PM500 wire coordinates (µm)\r",g.wire.axis[0],g.wire.axis[1],g.wire.axis[2]
		if (g.wire.Rmag > 0)
			printf "	R = {%g, %g, %g}, a rotation of %g¡	// wire positioner rotation vector\r",g.wire.R[0],g.wire.R[1],g.wire.R[2],g.wire.Rmag
			if (NumVarOrDefault("root:Packages:geometry:printVerbose",0))
				printf "			{%+.6f, %+.6f, %+.6f}	// rotation matrix from wire Rw\r",g.wire.R00, g.wire.R01, g.wire.R02
				printf "	Rw =	{%+.6f, %+.6f, %+.6f}\r",g.wire.R10, g.wire.R11, g.wire.R12
				printf "			{%+.6f, %+.6f, %+.6f}\r",g.wire.R20, g.wire.R21, g.wire.R22
				printf "	wire axis direction = {%.6f, %.6f, %.6f}	// direction of wire axis in Beam LIne coordinates (µm)\r",g.wire.axisR[0],g.wire.axisR[1],g.wire.axisR[2]
			endif
		endif
		if (!numtype(g.wire.F))
			printf "	F=%.2f						// F of wire for a wire scan raw PM500 units (µm)\r",g.wire.F
		endif
	endif
End
//
Function printDetector(d)							// print the details for passed detector geometry to the history window
	STRUCT detectorGeometry &d
	if (!(d.used))
		return 1
	endif

	printf "	Nx=%d, Ny=%d			// number of un-binned pixels in detector\r",d.Nx,d.Ny
	printf "	sizeX=%g, sizeY=%g		// size of detector (mm)\r",(d.sizeX/1000), (d.sizeY/1000)
	printf "	R = {%.7g, %.7g, %.7g}, a rotation of %.7g¡	// rotation vector\r",d.R[0],d.R[1],d.R[2],sqrt(d.R[0]*d.R[0] + d.R[1]*d.R[1] + d.R[2]*d.R[2])*180/PI
	printf "	P = {%g, %g, %g}					// translation vector (mm)\r",(d.P[0])/1000,(d.P[1])/1000,(d.P[2])/1000

	printf "	geometry measured on  '%s'\r",d.timeMeasured
	if (strlen(d.geoNote))
		printf "	detector note = '%s'\r",d.geoNote
	endif
	if (strlen(d.distortionMapFile))
		printf "	detector distortion file = '%s'\r",d.distortionMapFile
	endif
	printf "	detector ID = '%s'\r"d.detectorID
	if (NumVarOrDefault("root:Packages:geometry:printVerbose",0))
		printf "			{%+.6f, %+.6f, %+.6f}	// rotation matrix from R\r",d.rho00, d.rho01, d.rho02
		printf "	rho =	{%+.6f, %+.6f, %+.6f}\r",d.rho10, d.rho11, d.rho12
		printf "			{%+.6f, %+.6f, %+.6f}\r",d.rho20, d.rho21, d.rho22
	endif
	Variable chiMin, chiMax, tthMin, tthMax
	Variable/C chiZ, tthZ
	chiZ = chiRange(d)
	chiMin = min(chiMin,real(chiZ))
	chiMax = max(chiMax,imag(chiZ))
	tthZ =  tthRange(d)
	tthMin = min(tthMin,real(tthZ))
	tthMax = max(tthMax,imag(tthZ))
	printf "\tangle ranges:  chi = [%g, %g¡],  2th = [%g, %g¡]\r", real(chiZ),imag(chiZ),real(tthZ),imag(tthZ)
	return 0
End
//
Static Function/C chiRange(d)					// chi is rotation about the z-axis
	STRUCT detectorGeometry &d
	Variable chi,chiMin=+inf, chiMax=-inf
	Variable px,py
	Make/N=3/D/FREE xyz

	px = 0;			py = 0
	pixel2XYZ(d,px,py,xyz)					// convert pixel position to the beam line coordinate system
	chi = atan2(xyz[0],xyz[1]) * 180/PI
	chiMin = min(chiMin,chi)
	chiMax = max(chiMax,chi)

	px = 0;			py = d.Ny-1
	pixel2XYZ(d,px,py,xyz)
	chi = atan2(xyz[0],xyz[1]) * 180/PI
	chiMin = min(chiMin,chi)
	chiMax = max(chiMax,chi)

	px = d.Nx-1;	py = 0
	pixel2XYZ(d,px,py,xyz)
	chi = atan2(xyz[0],xyz[1]) * 180/PI
	chiMin = min(chiMin,chi)
	chiMax = max(chiMax,chi)

	px = d.Nx-1;	py = d.Ny-1
	pixel2XYZ(d,px,py,xyz)
	chi = atan2(xyz[0],xyz[1]) * 180/PI
	chiMin = min(chiMin,chi)
	chiMax = max(chiMax,chi)

	return cmplx(chiMin,chiMax)
End
//
Static Function/C tthRange(d)					// tth is the usual 2-theta
	STRUCT detectorGeometry &d
	Variable tth,tthMin=+inf, tthMax=-inf
	Variable px,py
	Make/N=3/D/FREE xyz

	px = 0;			py = 0
	pixel2XYZ(d,px,py,xyz)					// convert pixel position to the beam line coordinate system
	tth = acos(xyz[2]/norm(xyz)) * 180/PI
	tthMin = min(tthMin,tth)
	tthMax = max(tthMax,tth)

	px = 0;			py = d.Ny-1
	pixel2XYZ(d,px,py,xyz)
	tth = acos(xyz[2]/norm(xyz)) * 180/PI
	tthMax = max(tthMax,tth)

	px = d.Nx-1;	py = 0
	pixel2XYZ(d,px,py,xyz)
	tth = acos(xyz[2]/norm(xyz)) * 180/PI
	tthMin = min(tthMin,tth)
	tthMax = max(tthMax,tth)

	px = d.Nx-1;	py = d.Ny-1
	pixel2XYZ(d,px,py,xyz)
	tth = acos(xyz[2]/norm(xyz)) * 180/PI
	tthMin = min(tthMin,tth)
	tthMax = max(tthMax,tth)

	return cmplx(tthMin,tthMax)
End

Function CopymicroGeometry(f,i)					// copy a microGeometry structure
	STRUCT microGeometry &f, &i					// f is the destination, i is source
	CopySampleGeometry(f.s,i.s)					// copy Sample geometry
	f.Ndetectors = i.Ndetectors
	Variable m
	for (m=0;m<MAX_Ndetectors;m+=1)			// copy the detectors, even copy un-used detectors
		CopyDetectorGeometry(f.d[m],i.d[m])
	endfor
	CopyWireGeometry(f.wire,i.wire)				// copy the wire
End
//
Function CopyDetectorGeometry(f,i)					// copy a detector structure
	STRUCT detectorGeometry &f, &i					// f is the destination, i is source
	f.used = i.used
	f.Nx = i.Nx;			f.Ny = i.Ny
	f.sizeX = i.sizeX;		f.sizeY = i.sizeY
	f.R[0]=i.R[0];		f.R[1]=i.R[1];			f.R[2]=i.R[2];
	f.P[0]=i.P[0];		f.P[1]=i.P[1];			f.P[2]=i.P[2];
	f.timeMeasured = i.timeMeasured
	f.geoNote = i.geoNote
	f.detectorID = i.detectorID
	f.distortionMapFile = i.distortionMapFile
	f.rho00=i.rho00;		f.rho01=i.rho01;		f.rho02=i.rho02
	f.rho10=i.rho10;		f.rho11=i.rho11;		f.rho12=i.rho12
	f.rho20=i.rho20;		f.rho21=i.rho21;		f.rho22=i.rho22
End
//
Function CopyWireGeometry(f,i)					// copy a wire geometry structure, set f = i
	STRUCT wireGeometry &f, &i					// f is the destination, i is source
	f.origin[0]=i.origin[0];	f.origin[1]=i.origin[1];	f.origin[2]=i.origin[2];
	f.F = i.F
	f.dia = i.dia
	f.knife = i.knife
	f.axis[0]=i.axis[0];		f.axis[1]=i.axis[1];		f.axis[2]=i.axis[2];
	f.axisR[0]=i.axisR[0];	f.axisR[1]=i.axisR[1];	f.axisR[2]=i.axisR[2];
	f.R[0] = i.R[0];			f.R[1] = i.R[1];			f.R[2] = i.R[2]
	f.Rmag = i.Rmag
	f.R00=i.R00;				f.R01=i.R01;				f.R02=i.R02
	f.R10=i.R10;				f.R11=i.R11;				f.R12=i.R12
	f.R20=i.R20;				f.R21=i.R21;				f.R22=i.R22
End
Function CopySampleGeometry(f,i)					// copy a Sample geometry structure, set f = i
	STRUCT sampleGeometry &f, &i					// f is the destination, i is source
	f.O[0] = i.O[0];		f.O[1] = i.O[1];		f.O[2] = i.O[2]
	f.R[0] = i.R[0];		f.R[1] = i.R[1];		f.R[2] = i.R[2]
	f.Rmag = i.Rmag
	f.R00=i.R00;			f.R01=i.R01;			f.R02=i.R02
	f.R10=i.R10;			f.R11=i.R11;			f.R12=i.R12
	f.R20=i.R20;			f.R21=i.R21;			f.R22=i.R22
End

Function MicroGeometryBad(g)						// check for a valid or Invalid structure
	STRUCT microGeometry &g						// f is the destination, i is source
	Variable bad = SampleBad(g.s)
	Variable N = limit(round(g.Ndetectors),1,3)	// Ndetectors must be 1, 2, or 3
	bad += (g.Ndetectors != N) || numtype(g.Ndetectors)
	Variable m, i
	for (m=0,i=0;m<MAX_Ndetectors;m+=1)		// check the detectors
		if (g.d[m].used)
			bad += DetectorBad(g.d[m])
			i += 1
		endif
	endfor
	bad += (i != N)
	bad += WireBad(g.wire)							// check the wire
	return (bad>0)
End
//
Static Function DetectorBad(d)
	STRUCT detectorGeometry &d
	if (!(d.used))
		return 1
	endif
	Variable bad = (numtype(d.Nx + d.Nx + d.sizeX + d.sizeY + d.R[0] + d.R[1] + d.R[2] + d.P[0] + d.P[1] + d.P[2])>0)
	bad += (d.Nx<1 || d.Nx>5000)													// detector cannot have more than 5000 pixels along one edge
	bad += (d.Ny<1 || d.Ny>5000)
	bad += (d.sizeX<1e3 || d.sizeX>1e6)												// detector cannot be larger than 1m
	bad += (d.sizeY<1e3 || d.sizeY>1e6)
	bad += (abs(d.R[0])>2*PI || abs(d.R[1])>2*PI || abs(d.R[2])>2*PI)		// rotation cannot be more than 2¹
	bad += (abs(d.P[0])>2e6 || abs(d.P[0])>2e6 || abs(d.P[0])>2e6)			// P cannot be more than 2m in any direction
	return (bad>0)
End
//
Static Function WireBad(w)
	STRUCT wireGeometry &w
	Variable sxyz = abs(w.origin[0])+abs(w.origin[1])+abs(w.origin[2])
	Variable bad = numtype(sxyz) || abs(sxyz)>50e3
	bad += (w.dia<5 || w.dia>5100)
	bad += !(w.knife==0 || w.knife==1)			// knife must be 0 or 1
	bad += (abs(w.axis[0])+abs(w.axis[1])+abs(w.axis[2])) > 2
//	bad += (numtype(w.R[0] + w.R[1] + w.R[2])>0)
//	bad += ( abs(w.R[0])>6.3 || abs(w.R[1])>6.3 || abs(w.R[2])>6.3 )
	return (bad>0)
End
//
Static Function SampleBad(s)
	STRUCT sampleGeometry &s					// sample strucure
	Variable sxyz = abs(s.O[0])+abs(s.O[1])+abs(s.O[2])
	Variable bad = numtype(sxyz) || abs(sxyz)>50e3
	bad += (numtype(s.R[0] + s.R[1] + s.R[2]) > 0)
	return (bad>0)
End




// convert qvector to pixel position
Function/C q2pixel(d,qvec,[depth])					// returns pixel position as a complex number cmplx(px,py)
	STRUCT detectorGeometry &d
	Wave qvec										// qvec need not be normalized
	Variable depth									// sample depth measured along the beam
	depth = ParamIsDefault(depth) ? 0 : depth		// default is 0, the origin

	Variable px,py									// final pixel position, full chip unbinned 0 based pixels
	Wave qhat=root:Packages:geometry:q2pixel_qhat, ki=root:Packages:geometry:q2pixel_ki
	Wave kout=root:Packages:geometry:q2pixel_kout
	qhat = qvec
	normalize(qhat)

	ki = {0,0,1}									// ki = geo.ki[p],  incident beam direction
	//	normalize(ki)
	Variable qLen = -2*MatrixDot(qhat,ki)			// length of qhat, note (q^ dot -ki) always positive
	if (qLen<0)										// this occurs for theta<0, (we do not want to reflect from the back side)
		return cmplx(NaN,NaN)
	endif
	kout = qhat*qLen + ki							// kf - ki = q

	if (ParamIsDefault(depth))
		XYZ2pixel(d,kout,px,py)
	else
		XYZ2pixel(d,kout,px,py,depth=depth)
	endif

	return cmplx(px,py)
End



// convert px,py positions on detector into Q vector, assumes ki={0,0,1}
Function pixel2q(d,px,py,qhat,[depth])				// returns theta (rad)
	STRUCT detectorGeometry &d
	Variable px,py									// pixel position, 0 based, first pixel is (0,0), NOT (1,1)
	Wave qhat										// q-vector in beam line coords (optional)
	Variable depth									// sample depth measured along the beam

	Wave ki=root:Packages:geometry:pixel2q_ki, kout=root:Packages:geometry:pixel2q_kout
	ki = {0,0,1}									//	ki = geo.ki[p],  incident beam direction

	pixel2XYZ(d,px,py,kout)						// kout is in direction of pixel in beam line coords
	if (!ParamIsDefault(depth))					// a depth was passed, offset the sample position by depth*ki[]
		kout -= depth*ki							// koutDepth = d*ki + koutZero
	endif
	normalize(kout)

	Variable theta = acos(MatrixDot(kout,ki))/2	// ki.kf = cos(2theta), (radians)
	if (WaveExists(qhat))
		qhat = kout-ki								// qhat bisects kout and -ki
		normalize(qhat)
	endif
	return theta
End


Function pixel2XYZ(d,px,py,xyz)					// convert pixel position to the beam line coordinate system
	STRUCT detectorGeometry, &d
	Variable px,py									// pixel position on detector (full chip & zero based)
	Wave xyz											// 3-vector to receive the result, position in beam line coords (micron)

	peakCorrect(d,px,py)							// convert pixel on detector to an undistorted distance in pixels (takes zero based pixels)

	Variable xp,yp, zp								// x' and y' (requiring z'=0), detector starts centered on origin and perpendicular to z-axis
	xp = (px - 0.5*(d.Nx-1)) * d.sizeX/d.Nx		// (x' y' z'), position on detector
	yp = (py - 0.5*(d.Ny-1)) * d.sizeY/d.Ny

	xp += d.P[0]										// translate by P
	yp += d.P[1]
	zp = d.P[2]

	xyz[0] = d.rho00*xp + d.rho01*yp + d.rho02*zp	// xyz = rho x [ (x' y' z') + P ]
	xyz[1] = d.rho10*xp + d.rho11*yp + d.rho12*zp	// rho is pre-calculated from vector d.R
	xyz[2] = d.rho20*xp + d.rho21*yp + d.rho22*zp
End


Function/C XYZ2pixel(d,xyz,px,py,[depth])			// find pixel position (px,py) where vector xyz will intercept detector
	STRUCT detectorGeometry, &d
	Wave xyz										// 3-vector, a point in beam line coords giving the direction of ray (micron)
	Variable &px,&py								// pixel position on detector where ray xyz intercepts detector (full chip & zero based)
	Variable depth									// sample depth measured along the beam

	Variable xp,yp,zp
	Variable x=xyz[0], y=xyz[1], z=xyz[2]			// remember (xyz) = rho x [ (x' y' z') + P ]
	xp = d.rho00*x + d.rho10*y + d.rho20*z -d.P[0]	// so   xyz' = rho x xyz - P
	yp = d.rho01*x + d.rho11*y + d.rho21*z -d.P[1]	// this is now the coordinate of the pixel transformed into detector (i.e. prime) space
	zp = d.rho02*x + d.rho12*y + d.rho22*z -d.P[2]

	//	Variable xp0, yp0, zp0						// coordinate of origin transformed into prime space (the detector reference coordinates)
	//	xp0 = -d.P[0]								// coordinate of origin in prime space,  0 = rho x [ {x'y'z'} + P ]  -->   {x''y'z'} = -P
	//	yp0 = -d.P[1]
	//	zp0 = -d.P[2]	

	// for the line that goes from  {x',y',z'} to {x0',y'0',z0'}, where does it intersect the z'=0 plane?  NOTE: line goes from detector to sample
	// the condition is:  z' + t*(z0'-z') = 0,  from vector equation of a line r = t*Ær + r',  where Ær = (ro'-r')
	// I have NOT trapped (z'-z0')==0, this should only occur when detector surface contains the origin.
	// Note this line runs backwards from what is normally used. Going from r to r0, this makes t closer to 0 than 1, and so should improve accuracy

	Variable dxp, dyp, dzp							// Ær', subtract origin from (xp,yp,zp), this is direction of -kf
	dxp = -d.P[0] - xp
	dyp = -d.P[1] - yp
	dzp = -d.P[2] -  zp

	if (!ParamIsDefault(depth))					// a depth was passed, so modify (xp,yp,zp) to shift the origin by depth (shifts r')
		xp += d.rho20*depth						// remember:   xyz' = ( rho x xyz ) - P
		yp += d.rho21*depth						// shifted (xp,yp,zp) by depth*ki, all in prime space
		zp += d.rho22*depth						// This assumes that ki = {0,0,1}
	endif

	// find t,  where the line intercepts the detector
	Variable t = -zp/dzp							// for t==0, the point (x'y'z') already lies on the plane
	if (t>1)											// ray pointing backwards through other side of origin, so going away from detector
		px = NaN
		py = NaN
	else
		x = xp + t*dxp
		y = yp + t*dyp
//		z = zp + t*dzp								// z should be zero, that was the whole point! So this line not needed.
		px = x/d.sizeX*d.Nx + 0.5*(d.Nx-1)		// convert (xy0) in prime space to pixels
		py = y/d.sizeY*d.Ny + 0.5*(d.Ny-1)
	endif

	peakUncorrect(d,px,py)	// go from true peak position to the real peak position (put distortion back in), just invert peakCorrect()
	return cmplx(px,py)
End




//	Distortion corection, take a pixel position from detector and returns an ideal pixel position
// This assumes a zero based coordinate on input, origin is (0,0), not (1,1)
Static Function/C peakCorrect(d,px,py)				// returns cmplx(dx,dy), the change
	STRUCT detectorGeometry, &d
	Variable &px,&py								// un-binned full frame zero based pixel location on input, changed to distorttion corrected value at end

	Variable useDistortion = NumVarOrDefault("root:Packages:geometry:useDistortion",USE_DISTORTION_DEFAULT)
	useDistortion = useDistortion && Exists("root:Packages:geometry:xymap")==1
	if (useDistortion)
		Wave xymap = root:Packages:geometry:xymap	// distortion map
		return peakcorrection2(xymap,px,py)			// returns cmplx(dx,dy)
	else
		return cmplx(0,0)										// no distortion
	endif
End
//Function/C peakcorrection2(xymap,px,py)		// returns cmplx(dx,dy)
//
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


Static Function/C peakUncorrect(d,px,py)			// go from true peak to the measured peak position (put distortion back in), just invert peakCorrect()
	STRUCT detectorGeometry, &d
	Variable &px,&py								// un-binned full frame 0 based pixel with distortion correction, changed to un-distorted value at end

	Wave xymap=$""								// distortion map

	Variable useDistortion = NumVarOrDefault("root:Packages:geometry:useDistortion",USE_DISTORTION_DEFAULT) && WaveExists(xymap)
	if (!useDistortion)								// do not distort
		return cmplx(px,py)
	endif

	Variable maxIter=6, tol=1e-4					// at most 6 iterations, and a tolerance of 1e-4

	Variable pxd=px, pyd=py						// save the starting point (the required distorted position)
	Variable/C dz
	Variable i=0,err = Inf
	for (i=0; i<maxIter && err>tol; i+=1)			// end when it changes by less than tol or after maxIter iterations
		dz = peakCorrect(d,px,py)					// returns changed px,py
		err = abs(px-pxd)+abs(py-pyd)
		px = pxd - real(dz)
		py = pyd - imag(dz)
	endfor
	return cmplx(px,py)
End



Function GeometryUpdateCalc(g)						// update all internally calculated things in the structure
	STRUCT microGeometry &g

	Init_microGeo()

	Variable Rx, Ry, Rz								// used to make the rotation matrix rho from vector R
	Variable theta, c, s, c1
	if (g.Ndetectors>MAX_Ndetectors)
		String str
		sprintf str, "ERROR, g.Ndetectors is %d, which is bigger than max value of %d. Reducing it",g.Ndetectors,MAX_Ndetectors
		DoAlert 0, str
		printf "\r%s\r",str
		g.Ndetectors = MAX_Ndetectors
	endif

	SampleUpdateCalc(g.s)							// update all internally calculated things in the sample structure

	Variable i, N=0
	for (i=0;i<MAX_Ndetectors;i+=1)
		if (g.d[i].used)
			DetectorUpdateCalc(g.d[i])				// update all internally calculated things in the detector structures
			N += 1
		endif
	endfor
	g.Ndetectors = N

	if (Exists("root:Packages:geometry:xymap")==1)
		Wave xymap=root:Packages:geometry:xymap
		resetCenterDistortionMap(xymap,(g.d[0].Nx)/2,(g.d[0].Ny)/2)	// update distortion map so that correction at (xc,yc) goes to zero
	endif

	WireUpdateCalc(g.wire)							// update all internally calculated things in the wire structure
End
//
Static Function resetCenterDistortionMap(xymap,xc,yc)		// set distortion correction of xc & yc into note of xymap, used so correction to (xc,yc) is zero
	Wave xymap											// distortion map
	Variable xc,yc											// values that should get zero distortion correction, should be geo.xcent,geo.ycent
	if (!WaveExists(xymap))								// wave does not exist, do nothing
		return 1
	endif

	String wnote = note(xymap)
	wnote = ReplaceNumberByKey("dxCenter",wnote,0,"=")	// set correction to zero
	wnote = ReplaceNumberByKey("dyCenter",wnote,0,"=")
	Note/K xymap, wnote									// re-set value in wave note
	if (numtype(xc+yc))
		return 1											// invalid (xc,yc) passed so leave corrections at zero
	endif

	Variable px=xc-1, py=yc-1
	Variable/C dp = peakcorrection2(xymap,px,py)					// returns cmplx(dx,dy), px,py are changed
	wnote = ReplaceNumberByKey("dxCenter",wnote,real(dp),"=")	// add this to corrected to get proper value
	wnote = ReplaceNumberByKey("dyCenter",wnote,imag(dp),"=")
	Note/K xymap, wnote									// store values in wave note
End
//
Function DetectorUpdateCalc(d)						// update all internally calculated things in the detector structure
	STRUCT detectorGeometry &d
	if (!(d.used))
		return 1
	endif

	Variable Rx, Ry, Rz								// used to make the rotation matrix rho from vector R
	Variable theta, c, s, c1
	Variable i
	Rx=d.R[0]; Ry=d.R[1]; Rz=d.R[2]				// make the rotation matrix rho from vector R
	theta = sqrt(Rx*Rx+Ry*Ry+Rz*Rz)
	if (theta==0)										// no rotation, set to identity matrix
		d.rho00 = 1;		d.rho01 = 0;		d.rho02 = 0
		d.rho10 = 0;		d.rho11 = 1;		d.rho12 = 0
		d.rho20 = 0;		d.rho21 = 0;		d.rho22 = 1
		return 0
	endif

	c=cos(theta)
	s=sin(theta)
	c1 = 1-c
	Rx /= theta;	Ry /= theta;	Rz /= theta		// make |{Rx,Ry,Rz}| = 1

	d.rho00 = c + Rx*Rx*c1;		d.rho01 = Rx*Ry*c1 - Rz*s;	d.rho02 = Ry*s + Rx*Rz*c1		// this is the Rodrigues formula from:
	d.rho10 = Rz*s + Rx*Ry*c1;	d.rho11 = c + Ry*Ry*c1;		d.rho12 = -Rx*s + Ry*Rz*c1	// http://mathworld.wolfram.com/RodriguesRotationFormula.html
	d.rho20 = -Ry*s + Rx*Rz*c1;	d.rho21 = Rx*s + Ry*Rz*c1;	d.rho22 = c + Rz*Rz*c1
	return 0
End
//
Static Function SampleUpdateCalc(sa)				// update all internally calculated things in the structure
	STRUCT sampleGeometry &sa

	Variable Rx, Ry, Rz								// used to make the rotation matrix rho from vector R
	Variable theta, c, s, c1
	Rx=sa.R[0];	Ry=sa.R[1];	Rz=sa.R[2]
	theta = sqrt(Rx*Rx+Ry*Ry+Rz*Rz)			// angle in radians

	if (numtype(theta) || theta==0)					// rotation of sample stage
		sa.R[0] = 0;	sa.R[1] = 0;	sa.R[2] = 0;
		sa.Rmag = 0
		sa.R00 = 1;		sa.R01 = 0;		sa.R02 = 0	// matrix is identity
		sa.R10 = 0;		sa.R11 = 1;		sa.R12 = 0
		sa.R20 = 0;		sa.R21 = 0;		sa.R22 = 1
	else
		sa.Rmag = theta * 180/PI
		c=cos(theta)
		s=sin(theta)
		c1 = 1-c
		Rx /= theta;	Ry /= theta;	Rz /= theta	// make |{Rx,Ry,Rz}| = 1
		sa.R00 = c + Rx*Rx*c1;			sa.R01 = Rx*Ry*c1 - Rz*s;	sa.R02 = Ry*s + Rx*Rz*c1		// this is the Rodrigues formula from:
		sa.R10 = Rz*s + Rx*Ry*c1;		sa.R11 = c + Ry*Ry*c1;		sa.R12 = -Rx*s + Ry*Rz*c1	// http://mathworld.wolfram.com/RodriguesRotationFormula.html
		sa.R20 = -Ry*s + Rx*Rz*c1;	sa.R21 = Rx*s + Ry*Rz*c1;	sa.R22 = c + Rz*Rz*c1
	endif
End
//
Static Function WireUpdateCalc(w)
	STRUCT wireGeometry &w

	// normalize wire.axis
	Variable len = sqrt(w.axis[0]*w.axis[0] + w.axis[1]*w.axis[1] + w.axis[2]*w.axis[2])
	if (numtype(len)==0 && len>0)
		w.axis[0] /= len ;		w.axis[1] /= len;		w.axis[2] /= len
	else
		w.axis[0] = 1 ; 		w.axis[1] = 0 ;		w.axis[2] = 0
	endif

	// compute wire coordinates that would put wire center at the Si position (micron)
//	Variable ir2=1/sqrt(2), dZ
//	Variable H0=w.H0, Hyc=w.Hyc, F=w.F
//
//	w.xyz0[0] = w.X									// wire coordinates when wire is on the incident beam, computed (micron)
//	w.xyz0[1] = (H0-F)*ir2
//	w.xyz0[2] = (H0+F)*ir2
//
//	w.xyzSi[2] = (Hyc+F)*ir2						// Z of the Si position in wire units & beam line coordinates (micron)
//	dZ = (Hyc-H0)*ir2 								// dist (along z) from where wire intersects beam to Si

////		***********************************************************************************
//// NEEDS WORK, THIS IS WRONG!!!
//	w.xyzSi[1] = (H0-F)*ir2						// Y of the Si in wire units & beam line coordinates (micron)
//	w.xyzSi[0] = w.X									// X of the Si in wire units & beam line coordinates (micron)
//
////	w.xyzSi[1] = (H0-F)*ir2 - dZ*ki[2]/ki[1]	// Y of the Si in wire units & beam line coordinates (micron)
////	w.xyzSi[0] = w.X - dZ*ki[0]/ki[1]				// X of the Si in wire units & beam line coordinates (micron)
////		***********************************************************************************

	// fix up the rotation of wire frame
	Variable Rx, Ry, Rz								// used to make the rotation matrix rho from vector R
	Variable theta, c, s, c1
	Rx=w.R[0];	Ry=w.R[1];	Rz=w.R[2]
	theta = sqrt(Rx*Rx+Ry*Ry+Rz*Rz)			// angle in radians

	if (numtype(w.Rmag) || theta==0)				// rotation of wire stage
		w.R[0] = 0;		w.R[1] = 0;		w.R[2] = 0;
		w.Rmag = 0
		w.R00 = 1;		w.R01 = 0;		w.R02 = 0	// matrix is identity
		w.R10 = 0;		w.R11 = 1;		w.R12 = 0
		w.R20 = 0;		w.R21 = 0;		w.R22 = 1
		w.axisR[0] = w.axis[0]						// both axes the same, no rotation
		w.axisR[1] = w.axis[1]
		w.axisR[2] = w.axis[2]
	else
		w.Rmag = theta * 180/PI
		c=cos(theta)
		s=sin(theta)
		c1 = 1-c
		Rx /= theta;	Ry /= theta;	Rz /= theta	// make |{Rx,Ry,Rz}| = 1
		w.R00 = c + Rx*Rx*c1;			w.R01 = Rx*Ry*c1 - Rz*s;	w.R02 = Ry*s + Rx*Rz*c1	// this is the Rodrigues formula from:
		w.R10 = Rz*s + Rx*Ry*c1;		w.R11 = c + Ry*Ry*c1;		w.R12 = -Rx*s + Ry*Rz*c1	// http://mathworld.wolfram.com/RodriguesRotationFormula.html
		w.R20 = -Ry*s + Rx*Rz*c1;	w.R21 = Rx*s + Ry*Rz*c1;	w.R22 = c + Rz*Rz*c1

		Variable xx=w.axis[0], yy=w.axis[1], zz=w.axis[2]
		w.axisR[0] = w.R00*xx + w.R01*yy + w.R02*zz	// axisR = w.Rij x axis,   rotate by R (a small rotation)
		w.axisR[1] = w.R10*xx + w.R11*yy + w.R12*zz
		w.axisR[2] = w.R20*xx + w.R21*yy + w.R22*zz
	endif
End




//
//// calculate rotation matrix from detector tilts angles
//// note: in IDL matrices are defined as transpose of the usual matrices
//// 
//// usese geo to set ki[3], rded[3][3], and rde[3][3]
//Function GeometryUpdateCalc(geo)			// formerly callled detectortilt()
//	STRUCT microGeometry &geo
//
//	Init_microGeo()
//	// first convert all angles from degrees --> radians
//	Variable xbr = -geo.xbet*PI/180, xcr=-geo.xgam*PI/180
//	Variable xard=geo.xalfd*PI/180, xbrd=geo.xbetd*PI/180
//
//	// rotation matrix for incident beam, use it to set ki[]
//	//
//	// rot(xbr) =	1		0			0
//	//				0	  cos(xbr)	sin(xbr)
//	//				0	-sin(xbr)	cos(xbr)
//	// rot(xbr) rotates a vector about x negative. So +xbetd move the beam upward as it should.
//	//
//	// rot(xcr) =	  cos(xcr)	sin(xcr)	0
//	//				-sin(xcr)	cos(xcr)	0
//	//				  	0			0		1
//	// rot(xcr) rotates a vector about z (which points upward) in a negative direction.  So +xgam is backward.
//	//
//	// rde = rot(xbr) x root(xcr)
//	//
//	Wave rde = root:Packages:geometry:GeoUpdate_rde
//	rde[0][0] = cos(xcr)							// rotation matrix for incident beam
//	rde[1][0] = -sin(xcr)*cos(xbr)
//	rde[2][0] = sin(xbr)*sin(xcr)
//	rde[0][1] = sin(xcr)
//	rde[1][1] = cos(xbr)*cos(xcr)
//	rde[2][1] = -cos(xcr)*sin(xbr)
//	rde[0][2] = 0
//	rde[1][2] = sin(xbr)
//	rde[2][2] = cos(xbr)
//	Wave ki = root:Packages:geometry:GeoUpdate_ki
////	ki = -rde[1][p]									// ki = (0,-1,0) x rde
//	ki = {0,-1,0}
//	MatrixMultiply ki/T, rde
//	Wave M_product=M_product
//	ki = M_product
//	KillWaves M_product
//	normalize(ki)
//	geo.ki[0] = ki[0]								// in Wenge coordinates (X=out door, Y=upstream, Z=up)
//	geo.ki[1] = ki[1]								// note, cannot automatically iterate a vector assignment when it is in a structure
//	geo.ki[2] = ki[2]
//
//	geo.rde00 =   rde[0][0]							// rotation matrix that goes from Wenge system to IDEAL Beam-Line
//	geo.rde10 =   rde[2][0]							// this is like rde[][], but includes the swap of Y=-z and Z=y
//	geo.rde20 = -rde[1][0]							// NOTE, this is a rotation matrix, so its Transpose is its inverse
//	geo.rde01 =   rde[0][1]
//	geo.rde11 =   rde[2][1]
//	geo.rde21 = -rde[1][1]
//	geo.rde02 =   rde[0][2]
//	geo.rde12 =   rde[2][2]
//	geo.rde22 = -rde[1][2]
//	//				// an example of how to use geo.rdexx
//	//			qBL[0] = qW[0]*geo.rde00 + qW[1]*geo.rde01 + qW[2]*geo.rde02	// transform Wenge to Ideal Beam-line
//	//			qBL[1] = qW[0]*geo.rde10 + qW[1]*geo.rde11 + qW[2]*geo.rde12	//   with incident beam along {0,0,1}
//	//			qBL[2] = qW[0]*geo.rde20 + qW[1]*geo.rde21 + qW[2]*geo.rde22
//
//	// rotation matrix for detector
//	//
//	// rot(a) =  cos(a)		0	  sin(a)			// rotation about y axis
//	//			  	0			1		1
//	//			-sin(a)		0	  cos(a)
//	//
//	// rot(b) =	1		0			0				// rotation about x axis
//	//				0	  cos(xbr)	sin(xbr)
//	//				0	-sin(xbr)	cos(xbr)
//	//
//	// rded = rot(xard) x root(xbrd)
//	//
//	geo.rded00 = cos(xard)
//	geo.rded10 = 0
//	geo.rded20 = -sin(xard)
//	geo.rded01 = -sin(xard)*sin(xbrd)
//	geo.rded11 = cos(xbrd)
//	geo.rded21 = -sin(xbrd)*cos(xard)
//	geo.rded02 = sin(xard)*cos(xbrd)
//	geo.rded12 = sin(xbrd)
//	geo.rded22 = cos(xard)*cos(xbrd)
//
//	Variable len = sqrt(geo.wire.axis[0]*geo.wire.axis[0] + geo.wire.axis[1]*geo.wire.axis[1] + geo.wire.axis[2]*geo.wire.axis[2])
//	if (numtype(len)==0)
//		geo.wire.axis[0] /= len  ;  geo.wire.axis[1] /= len  ;  geo.wire.axis[2] /= len
//	else
//		geo.wire.axis[0] = 1  ;  geo.wire.axis[1] = 0  ;  geo.wire.axis[2] = 0
//	endif
//
//	// compute wire coordinates that would put wire center at the Si position (µm)
//	Variable ir2=1/sqrt(2), dZ
//	Variable H0=geo.wire.H0, Hyc=geo.wire.Hyc, F=geo.wire.F
//
//	geo.wire.xyz0[0] = geo.wire.X						// wire coordinates when wire is on the incident beam, computed (µm)
//	geo.wire.xyz0[1] = (H0-F)*ir2
//	geo.wire.xyz0[2] = (H0+F)*ir2
//
//	geo.wire.xyzSi[2] = (Hyc+F)*ir2					// Z of the Si position in wire units & beam line coordinates (µm)
//	dZ = (Hyc-H0)*ir2 								// dist (along z) from where wire intersects beam to Si
//	geo.wire.xyzSi[1] = (H0-F)*ir2 - dZ*ki[2]/ki[1]	// Y of the Si in wire units & beam line coordinates (µm)
//	geo.wire.xyzSi[0] = geo.wire.X - dZ*ki[0]/ki[1]	// X of the Si in wire units & beam line coordinates (µm)
//
//	Wave xymap=root:Packages:geometry:xymap
/////	resetCenterDistortionMap(xymap,geo.xcent,geo.ycent)	// update distortion map so that correction at (xc,yc) goes to zero
//End


Function UpdateDefaultGeometryStruct(g,[local])			// Update the default location with values in g
	STRUCT microGeometry &g								// returns 0 if something set, 0 is nothing done
	Variable local												// if true will also update geoStructStr in the current folder if it already exists (It will NOT create local geoStructStr)
	local = ParamIsDefault(local) ? 0 : local

	String/G root:Packages:geometry:geoStructStr=""		// use default location
	SVAR geoStructStr = root:Packages:geometry:geoStructStr
	GeometryUpdateCalc(g)
	StructPut/S/B=2 g, geoStructStr
	SavePackagePreferences/FLSH=1 "microGeo","microGeoNPrefs",0,g

	if (exists(":geoStructStr")==2 && local)				// if geoStructStr also exists in current folder, update it too
		SVAR geoStructStr = :geoStructStr
		StructPut/S/B=2 g, geoStructStr
	endif
	return 0
End
//Function UpdateDefaultGeometryStruct(g)					// Update the default location with values in g
//	STRUCT microGeometry &g								// returns 0 if something set, 0 is nothing done
//
//	String/G root:Packages:geometry:geoStructStr=""		// use default location
//	SVAR geoStructStr = root:Packages:geometry:geoStructStr
//	GeometryUpdateCalc(g)
//	StructPut/S/B=2 g, geoStructStr
//	return 0
//End
//
Function FillGeometryStructDefault(g)						//fill the geometry structure with test values
	STRUCT microGeometry &g								// returns 0 if something set, 0 is nothing done

	String strStruct=StrVarOrDefault(":geoStructStr","")	// set to values in current directory
	if (strlen(strStruct)<1)
		strStruct=StrVarOrDefault("root:Packages:geometry:geoStructStr","")	// try the default values
	endif
	if (strlen(strStruct)>1)
		StructGet/S/B=2 g, strStruct						// found structure information, load into g
	else
		LoadPackagePreferences/MIS=1 "microGeo","microGeoNPrefs",0,g
		if (V_flag)
			return 1											// did nothing, nothing found
		endif
		StructPut/S/B=2 g, strStruct							// keep a local copy
		String/G root:Packages:geometry:geoStructStr = strStruct
	endif
	GeometryUpdateCalc(g)
	return 0
End


Function detectorNumFromID(ID)							// returns detector number {0,1,2} or -1 for an error, given the detector ID string
	String ID													// detector ID input, (something like "PE0820, 763-1807")
	if (strlen(ID)<1)
		return -1
	endif

	ID = ReplaceString(",", ID, "")
	ID = ReplaceString("_", ID, "-")
	ID = ReplaceString("  ", ID, " ")

	STRUCT microGeometry g
	FillGeometryStructDefault(g)
	Variable i
	for (i=0;i<MAX_Ndetectors;i+=1)						// search each of the g.d[i].detectorID and return at the first match
		if (strsearch(g.d[i].detectorID,ID,0)==0)			// does g.d[i].detectorID start with the passed ID?
			return i
		endif
	endfor
	return -1													// nothing matched, return error
End



Function WriteGeoToFile(fileName,path,g,fileNote,type)
	String fileName						// full path name to the file
	String path							// name of an Igor path to use
	STRUCT microGeometry &g			// the structure to fill from the file
	String fileNote						// a note about this file
	Variable type						// 0=key file,  1=xml file

	if ( !(type==0 || type==1) )
		type = 1						// currently default to old style "$key value" files
		Prompt type, "type of geometry file",popup,"key file (old);xml file (new)"
		DoPrompt "geometry file",type
		if (V_flag)
			return 1
		endif
		type -= 1						// 0=key file,  1=xml file
	endif

	if( (strlen(fileName)+strlen(fileNote))==0 )
		Prompt fileNote, "add a descriptive note to this file?"
		DoPrompt "descriptive note",fileNote
		if (V_flag)
			return 1
		endif
	endif

	String out="", fileType=""
	if (type==0)
		out = Geo2KeyValueStr(g,fileNote)
		fileType = "TEXT"
	elseif (type==1)
		out = Geo2xmlStr(g,fileNote)
		fileType = ".xml"
	endif
	if (strlen(out)<1)
		return 1
	endif

	Variable f
	if (strlen(fileName)<1)				// no file name passed
		fileName = SelectString(strlen(fileName),NewGeoFileName(g),fileName)
		Open/D/C="R*ch"/M="new geometry parameters file"/T=fileType f as fileName
		fileName = S_fileName
	endif
	Open/C="R*ch"/M="new geometry parameters file"/P=$path/T=fileType/Z f as fileName
	fileName = S_fileName
	if (V_flag)
		DoAlert 0, "nothing written to file"
		return 1
	endif
	FBinWrite f, out
	Close f
	String str="wrote geometry to file '"+fileName+"'"
	print str
	DoAlert 0,str
	return 0
End
//
Static Function/T NewGeoFileName(g)
	STRUCT microGeometry &g
	Variable month,day,year,hour,minute,second,TZ, i, epoch=-1		// not using time zone
	String smonth
	for (i=0;i<g.Ndetectors;i+=1)
		if (g.d[i].used)
			// g.d[i].timeMeasured looks like:   "Thu, Jul 28, 2011, 16:17:44 (-5)"
			sscanf g.d[i].timeMeasured, "%3s, %3s %d, %d, %02d:%02d:%02d (%g)",smonth,smonth,day,year,hour,minute,second,TZ
			if (V_flag==8)
				month = WhichListItem(smonth, "Jan;Feb;Mar;Apr;May;Jun;Jul;Aug;Sep;Oct;Nov;Dec")+1
				epoch = max(epoch, date2secs(year,month,day) + 3600*hour+60*minute+second)
			endif
		endif
	endfor
	epoch = epoch<date2secs(2005,1,1) ? DateTime : epoch				// no valid times before 2005, use current time
	return "geoN_"+Secs2Date(epoch,-2,"-")+"_"+ReplaceString(":",Secs2Time(epoch,3),"-")
End


Function ReadGeoFromfile(fileName,path,g)
	String fileName							// full path name to the file
	String path								// name of an Igor path to use
	STRUCT microGeometry &g				// the structure to fill from the file

	Variable f
	#if (NumberByKey("IGORVERS", IgorInfo(0))<6.1)
		Open/M="geoN file"/P=$path/R/T="????"/Z=2 f as fileName
	#else
		String fileFilters = "Text Files (*.xml,*.txt):.txt,.xml;All Files:.*;"
		Open/F=fileFilters/M="geoN file"/P=$path/R/Z=2 f as fileName
	#endif
	if (strlen(S_fileName)<1 || !f)
		return 1
	endif
	fileName = S_fileName
	FStatus f
	Variable len = min(V_logEOF,500)
	String buf=""
	buf = PadString(buf,len,0x20)
	FBinRead f, buf
	Close f
	Variable err
	if (strsearch(buf, "<?xml",0)>=0)
		err = ReadGeoFromXMLfile(fileName,path,g)
	else
		err = ReadGeoFromKeyFile(fileName,path,g)
	endif
End



// string of geometry structure suitable for writing to a geometry keyword file
Static Function/T Geo2KeyValueStr(g,fileNote)
	STRUCT microGeometry &g			// the structure to fill from the file
	String fileNote						// a note about this file

	GeometryUpdateCalc(g)				// calculate other values
	String out, str
	out = "$filetype		geometryFileN\n"
	sprintf str, "$dateWritten	%s\n", date();	out += str
	sprintf str, "$timeWritten	%s (%g)\n", Secs2Time(DateTime,3,1),date2secs(-1,-1,-1)/3600;	out += str
	sprintf str, "$EPOCH			%.0f					// seconds from midnight January 1, 1904\n",DateTime;	out += str
	if (strlen(fileNote))
		sprintf str,"$fileNote		%s\n",ReplaceString("\r",fileNote,"\\r"); out += str
	endif

	if (!SampleBad(g.s))
		out += "\n// Sample\n"
		sprintf str,"$SampleOrigin	{%.2f,%.2f,%.2f}			// sample origin in raw PM500 units (micron)\n",g.s.O[0],g.s.O[1],g.s.O[2];	out += str
		sprintf str,"$SampleRot		{%.8f,%.8f,%.8f}	// sample positioner rotation vector (length is angle in radians)\n",g.s.R[0],g.s.R[1],g.s.R[2]; out += str
	endif

	out += "\n// Detectors\n"
	sprintf str,"$Ndetectors		%d							// number of detectors in use, must be <= MAX_Ndetectors\n",g.Ndetectors;	out += str
	Variable i
	String pre
	for (i=0;i<MAX_Ndetectors;i+=1)
		if (!(g.d[i].used))
			continue
		endif
		sprintf pre,"d%d_",i	;	out += "\n"
		sprintf str,"$%sNx			%d						// number of un-binned pixels in full detector\n",pre,g.d[i].Nx;	out += str
		sprintf str,"$%sNy			%d\n",pre,g.d[i].Ny;	out += str
		sprintf str,"$%ssizeX		%.3f						// size of CCD (mm)\n",pre,(g.d[i].sizeX)/1000;	out += str
		sprintf str,"$%ssizeY		%.3f\n",pre,(g.d[i].sizeY/1000); 		out += str
		sprintf str,"$%sR			{%.8f,%.8f,%.8f}	// rotation vector (length is angle in radians)\n",pre,g.d[i].R[0],g.d[i].R[1],g.d[i].R[2]; 		out += str
		sprintf str,"$%sP			{%.3f,%.3f,%.3f}		// translation vector (mm)\n",pre,(g.d[i].P[0]/1000),(g.d[i].P[1]/1000),(g.d[i].P[2]/1000);  out += str
		if (strlen(g.d[i].timeMeasured))
			sprintf str,"$%stimeMeasured	%s	// when this geometry was calculated\n",pre,g.d[i].timeMeasured; out += str
		endif
		if (strlen(g.d[i].geoNote))
			sprintf str,"$%sgeoNote	%s\n",pre,g.d[i].geoNote; out += str
		endif
		if (strlen(g.d[i].detectorID))
			sprintf str,"$%sdetectorID	%s			// unique detector ID\n",pre,g.d[i].detectorID; out += str
		endif
		if (strlen(g.d[i].distortionMapFile))
			sprintf str,"$%sdistortionMapFile	%s			// name of file with distortion map\n",pre,g.d[i].distortionMapFile; out += str
		endif
	endfor

	if (!WireBad(g.wire))
		out += "\n// Wire\n"
		sprintf str,"$wireDia		%.2f						// diameter of wire (micron)\n",g.wire.dia;	 out += str
		sprintf str,"$wireKnife		%g							// true if wire on a knife edge, false for free-standing wire\n",g.wire.knife; out += str
		sprintf str,"$wireOrigin		{%.2f,%.2f,%.2f}			// wire origin in raw PM500 frame (micron)\n",g.wire.origin[0],g.wire.origin[1],g.wire.origin[2]; out += str
		if (g.wire.Rmag>0)
			sprintf str,"$wireRot		{%.8f,%.8f,%.8f}	// wire positioner rotation vector (length is angle in radians)\n",g.wire.R[0],g.wire.R[1],g.wire.R[2]; out += str
		endif
		sprintf str,"$wireAxis		{%.6f,%.6f,%.6f}	// unit vector along wire axis, usually close to (1,0,0)\n",g.wire.axis[0],g.wire.axis[1],g.wire.axis[2]; out += str
		if (!numtype(g.wire.F))
			sprintf str,"$wireF			%.2f						// F of wire for a constant F wire scan (raw PM500 units)\n",g.wire.F; out += str
		endif
	endif
	return out
End
//
//// write the geometry structure to a keyword file
//Static Function WriteGeoToKeyFile(fileName,path,g,fileNote)
//	String fileName						// full path name to the file
//	String path							// name of an Igor path to use
//	STRUCT microGeometry &g			// the structure to fill from the file
//	String fileNote						// a note about this file
//
//	GeometryUpdateCalc(g)				// calculate other values
//	Variable f								// file id
//	if (strlen(path)<1 || strlen(fileName)<1)
//		Open/D/C="R*ch"/M="new geometry parameters file"/T="TEXT" f as fileName
//		fileName = S_fileName
//	endif
//	Open/C="R*ch"/M="new geometry parameters file"/P=$path/T="TEXT"/Z f as fileName
//	fileName = S_fileName
////	Open/C="R*ch"/M="new geometry parameters file"/P=$path/T="TEXT"/Z=2 f as fileName
//	if (V_flag)
//		DoAlert 0, "nothing written to file"
//		return 1
//	endif
//	if( (strlen(path)+strlen(fileNote))==0 )
//		Prompt fileNote, "add a descriptive note to this file?"
//		DoPrompt "descriptive note",fileNote
//		if (V_flag)
//			fileNote=""
//		endif
//	endif
//	fprintf f,"$filetype		geometryFileN\n"
//	fprintf f, "$dateWritten	%s\n", date()
//	fprintf f, "$timeWritten	%s (%g)\n", Secs2Time(DateTime,3,1),date2secs(-1,-1,-1)/3600
//	fprintf f, "$EPOCH			%.0f					// seconds from midnight January 1, 1904\n",DateTime
//	if (strlen(fileNote))
//		fprintf f,"$fileNote		%s\n",ReplaceString("\r",fileNote,"\\r")
//	endif
//
//	if (!SampleBad(g.s))
//		fprintf f,"\n// Sample\n"
//		fprintf f,"$SampleOrigin	{%.2f,%.2f,%.2f}			// sample origin in raw PM500 units (micron)\n",g.s.O[0],g.s.O[1],g.s.O[2]
//		fprintf f,"$SampleRot		{%.8f,%.8f,%.8f}	// sample positioner rotation vector (length is angle in radians)\n",g.s.R[0],g.s.R[1],g.s.R[2]
//	endif
//
//	fprintf f,"\n// Detectors\n"
//	fprintf f,"$Ndetectors		%d							// number of detectors in use, must be <= MAX_Ndetectors\n",g.Ndetectors
//	Variable i
//	String pre
//	for (i=0;i<MAX_Ndetectors;i+=1)
//		if (!(g.d[i].used))
//			continue
//		endif
//		sprintf pre,"d%d_",i
//		fprintf f,"\n"
//		fprintf f,"$%sNx			%d						// number of un-binned pixels in full detector\n",pre,g.d[i].Nx
//		fprintf f,"$%sNy			%d\n",pre,g.d[i].Ny
//		fprintf f,"$%ssizeX		%.3f						// size of CCD (mm)\n",pre,(g.d[i].sizeX)/1000
//		fprintf f,"$%ssizeY		%.3f\n",pre,(g.d[i].sizeY/1000)
//		fprintf f,"$%sR			{%.8f,%.8f,%.8f}	// rotation vector (length is angle in radians)\n",pre,g.d[i].R[0],g.d[i].R[1],g.d[i].R[2]
//		fprintf f,"$%sP			{%.3f,%.3f,%.3f}		// translation vector (mm)\n",pre,(g.d[i].P[0]/1000),(g.d[i].P[1]/1000),(g.d[i].P[2]/1000)
//		if (strlen(g.d[i].timeMeasured))
//			fprintf f,"$%stimeMeasured	%s	// when this geometry was calculated\n",pre,g.d[i].timeMeasured
//		endif
//		if (strlen(g.d[i].geoNote))
//			fprintf f,"$%sgeoNote	%s\n",pre,g.d[i].geoNote
//		endif
//		if (strlen(g.d[i].detectorID))
//			fprintf f,"$%sdetectorID	%s			// unique detector ID\n",pre,g.d[i].detectorID
//		endif
//		if (strlen(g.d[i].distortionMapFile))
//			fprintf f,"$%sdistortionMapFile	%s			// name of file with distortion map\n",pre,g.d[i].distortionMapFile
//		endif
//	endfor
//
//	if (!WireBad(g.wire))
//		fprintf f,"\n// Wire\n"
//		fprintf f,"$wireDia		%.2f						// diameter of wire (micron)\n",g.wire.dia
//		fprintf f,"$wireKnife		%g							// true if wire on a knife edge, false for free-standing wire\n",g.wire.knife
//		fprintf f,"$wireOrigin		{%.2f,%.2f,%.2f}			// wire origin in raw PM500 frame (micron)\n",g.wire.origin[0],g.wire.origin[1],g.wire.origin[2]
//		if (g.wire.Rmag>0)
//			fprintf f,"$wireRot		{%.8f,%.8f,%.8f}	// wire positioner rotation vector (length is angle in radians)\n",g.wire.R[0],g.wire.R[1],g.wire.R[2]
//		endif
//		fprintf f,"$wireAxis		{%.6f,%.6f,%.6f}	// unit vector along wire axis, usually close to (1,0,0)\n",g.wire.axis[0],g.wire.axis[1],g.wire.axis[2]
//		if (!numtype(g.wire.F))
//			fprintf f,"$wireF			%.2f						// F of wire for a constant F wire scan (raw PM500 units)\n",g.wire.F
//		endif
//	endif
//	Close f
//	String str
//	sprintf str, "finished wrting geometry to file '%s'",fileName
//	print str
//	DoAlert 0,str
//	return 0
//End


// read in the geometry structure from a keyword file
Static Function ReadGeoFromKeyFile(fileName,path,g)
	String fileName							// full path name to the file
	String path								// name of an Igor path to use
	STRUCT microGeometry &g				// the structure to fill from the file

	String list = keyStrFromFile(fileName,"geometryFileN",path)// read in all of the tagged values into a keyword list
	Variable err = GeoFromKeyValueList(list,g)
	if (err==0)
		printf "Loaded geometry information from '%s'\r",StringByKey("keyStrFileName",list,"=")
		printf "   geometry file was set on  %s, %s\r",StringByKey("dateWritten",list,"="),StringByKey("timeWritten",list,"=")
		String str = StringByKey("fileNote",list,"=")
		if (strlen(str))
			printf ",   file note='%s'\r",str
		endif
	endif
	return err
End
//
Static Function GeoFromKeyValueList(list,g)
	String list								// "key=value" list
	STRUCT microGeometry &g				// the structure to fill from the file

	if (strlen(list)<1)
		return 1
	endif

	Variable N = NumberByKey("Ndetectors",list,"=")
	N = N>0 ? round(N) : 0
	if (N<0  || N>3)
		return 1
	endif
	g.Ndetectors = N

	Variable i, value, xx,yy,zz
	sscanf StringByKey("SampleOrigin",list,"="),"{%g,%g,%g}",xx,yy,zz
	if (numtype(xx+yy+zz)==0 && V_flag==3)
		g.s.O[0] = xx;		g.s.O[1] = yy;		g.s.O[2] = zz
	else
		g.s.O[0] = NaN;	g.s.O[1] = NaN;	g.s.O[2] = NaN
	endif
	sscanf StringByKey("SampleRot",list,"="),"{%g,%g,%g}",xx,yy,zz
	if (numtype(xx+yy+zz)==0 && V_flag==3)
		g.s.R[0] = xx;		g.s.R[1] = yy;		g.s.R[2] = zz
	else
		g.s.R[0] = 0;		g.s.R[1] = 0;		g.s.R[2] = 0
	endif

	String pre, str
	for (i=0;i<MAX_Ndetectors;i+=1)
		sprintf pre,"d%d_",i
		g.d[i].used = !numtype(NumberByKey(pre+"sizeX",list,"="))	// if no di_sizeX, then invalid
		if (!(g.d[i].used))
			continue
		endif
		value = NumberByKey(pre+"Nx",list,"=");			g.d[i].Nx = numtype(value) ? g.d[i].Nx : value
		value = NumberByKey(pre+"Ny",list,"=");			g.d[i].Ny = numtype(value) ? g.d[i].Ny : value
		value = NumberByKey(pre+"sizeX",list,"=");		g.d[i].sizeX = numtype(value) ? g.d[i].sizeX : value*1e3 // file uses mm, I need µm
		value = NumberByKey(pre+"sizeY",list,"=");		g.d[i].sizeY = numtype(value) ? g.d[i].sizeY : value*1e3
		sscanf StringByKey(pre+"R",list,"="),"{%g,%g,%g}",xx,yy,zz
		if (numtype(xx+yy+zz)==0 && V_flag==3)
			g.d[i].R[0] = xx
			g.d[i].R[1] = yy
			g.d[i].R[2] = zz
		endif
		sscanf StringByKey(pre+"P",list,"="),"{%g,%g,%g}",xx,yy,zz
		if (numtype(xx+yy+zz)==0 && V_flag==3)
			g.d[i].P[0] = xx*1e3				 	// file uses mm, I need µm
			g.d[i].P[1] = yy*1e3
			g.d[i].P[2] = zz*1e3
		endif
		str = StringByKey(pre+"timeMeasured",list,"=");		g.d[i].timeMeasured = SelectString(strlen(str),g.d[i].timeMeasured,str)
		str = StringByKey(pre+"geoNote",list,"=");				g.d[i].geoNote = SelectString(strlen(str),g.d[i].geoNote,str)
		str = StringByKey(pre+"detectorID",list,"=");			g.d[i].detectorID = SelectString(strlen(str),g.d[i].detectorID,str)
		str = StringByKey(pre+"distortionMapFile",list,"=");	g.d[i].distortionMapFile = SelectString(strlen(str),g.d[i].distortionMapFile,str)
	endfor

	value = NumberByKey("wireDia",list,"=");	g.wire.dia = numtype(value) ? g.wire.dia : value
	value = NumberByKey("wireKnife",list,"=");	g.wire.knife = numtype(value) ? g.wire.knife : value
	value = NumberByKey("wireF",list,"=");		g.wire.F = numtype(value) ? g.wire.F : value
	sscanf StringByKey("wireOrigin",list,"="),"{%g,%g,%g}",xx,yy,zz
	if (numtype(xx+yy+zz)==0 && V_flag==3)
		g.wire.origin[0]=xx;	g.wire.origin[1]=yy;	g.wire.origin[2]=zz
	else
		g.wire.origin[0]=0;	g.wire.origin[1]=0;	g.wire.origin[2]=0
	endif
	sscanf StringByKey("wireAxis",list,"="),"{%g,%g,%g}",xx,yy,zz
	if (numtype(xx+yy+zz)==0 && V_flag==3)
		g.wire.axis[0]=xx;	g.wire.axis[1]=yy;	g.wire.axis[2]=zz
	endif
	sscanf StringByKey("wireRot",list,"="),"{%g,%g,%g}",xx,yy,zz
	if (numtype(xx+yy+zz)==0 && V_flag==3)
		g.wire.R[0] = xx;		g.wire.R[1] = yy;		g.wire.R[2] = zz
	else
		g.wire.R[0] = 0;		g.wire.R[1] = 0;		g.wire.R[2] = 0
	endif

	GeometryUpdateCalc(g)							// calculate other values
	return 0
End
//
//Function twrite()
//	STRUCT microGeometry g
//	GeoReferenceOrientation(g)						//fill the geometry structure with reference values
////	FillGeometryStructDefault(g)				//fill the geometry structure with current values
//	Variable err
//	err = WriteGeoToKeyFile("testMe","home",g,"this is the first test file N")
//	if (err)
//		print "WriteGeoToKeyFile() returned error = ",err
//	endif
//	err = ReadGeoFromKeyFile("testMe","home",g)
//	if (err)
//		print "ReadGeoFromKeyFile() returned error = ",err
//	endif
//	printGeometry(g)
//End
//
Static Function/S keyStrFromFile(fname,ftype,path)	// read in all of the tagged values from a file into a keyword list
	String fname											// full path name to file with tagged geometry values
	String ftype											// the required file identifier, included as a tag (ftype is optional)
	String path											// name of Igor path

	Variable refNum
	Open/M="file containing tagged values"/P=$path/R/Z=2 refNum as fname
	if (strlen(S_fileName)<1 || !refNum)
		return ""
	endif

	Variable i
	String line
	if (strlen(ftype))									// an ftype is present, so check the file
		Variable OK=0
		FReadLine refNum, line
		line[strlen(line)-1,Inf] = " "					// change the trailing term char to a space
		if (strsearch(line,"$filetype",0)==0)			// starts with "$filetype"
			OK = char2num(line[9])<=32				// $filetype ends with whitespace
			i = strsearch(line,"//",0)					// search for comment identifier
			i = (i<0) ? strlen(line) : i					// points to start of comment id (or end of string)
			line = line[0,i-1]							// trim off comment
			line = ReplaceString("\t",line,";")			// change all tabs, commas, and spaces to semi-colons
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
		if (char2num(line)==dollar)					// if it starts with a $, probable tag so process it
			i = strsearch(line,"//",0)					// strip off comments
			if (i>=0)
				line = line[0,i-1]
			endif
			for (i=0;char2num(line[i+1])>32;i+=1)	// find end of tag, it ends with a space or lower
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
	while (strlen(line))									// go until EOF
	Close refNum
	return list
End


// string of geometry structure suitable for writing to an geometry xml file
Static Function/T Geo2xmlStr(g,fileNote)
	STRUCT microGeometry &g			// the structure to fill from the file
	String fileNote						// a note about this file

	GeometryUpdateCalc(g)				// calculate other values
	String out,str
	out = "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n\n"
	out += "<geoN xmlns=\"http://sector34.xray.aps.anl.gov/34ide/geoN\">\n"
	sprintf str, "	<dateWritten>%s</dateWritten>\n", date()	; out += str
	sprintf str, "	<timeWritten>%s (%g)</timeWritten>\n", Secs2Time(DateTime,3,1),date2secs(-1,-1,-1)/3600	; out += str
	sprintf str, "	<EPOCH start=\"midnight Jan 1, 1904\" unit=\"sec\">%.0f</EPOCH>\n", DateTime ; out += str
	if (strlen(fileNote))
		sprintf str, "	<fileNote>%s</fileNote>\n",ReplaceString("\r",fileNote,"\\r")
		out += str
	endif

	if (!SampleBad(g.s))
		out += "\n	<Sample>\n"
		sprintf str, "		<Origin unit=\"micron\">%.8g %.8g %.8g</Origin>	<!-- sample origin in raw PM500 units (micron) -->\n",g.s.O[0],g.s.O[1],g.s.O[2] ; out += str
		sprintf str, "		<R unit=\"radian\">%.9g %.9g %.9g</R>\n",g.s.R[0],g.s.R[1],g.s.R[2]; 	out += str
		out += "	</Sample>\n"
	endif

	sprintf str, "\n	<Detectors Ndetectors=\"%g\">			<!-- Ndetectors, number of detectors in use, must be <= MAX_Ndetectors -->\n",g.Ndetectors
	out += str
	Variable i, comment=1
	String pre
	for (i=0;i<MAX_Ndetectors;i+=1)
		if (!(g.d[i].used))
			continue
		endif
		sprintf str, "		<Detector N=\"%d\">\n",i
		out += str
		sprintf str, "			<Npixels>%d %d</Npixels>%s\n",g.d[i].Nx,g.d[i].Ny,SelectString(comment,"","		<!-- Nx,Ny is number of un-binned pixels in full detector -->")
		out += str
		sprintf str, "			<size unit=\"mm\">%.7g %.7g</size>%s\n",(g.d[i].sizeX)/1000,(g.d[i].sizeY)/1000,SelectString(comment,"","	<!-- sizeX,sizeY otuside size of full detector -->")
		out += str
		sprintf str, "			<R unit=\"radian\">%.9g %.9g %.9g</R>%s\n",g.d[i].R[0],g.d[i].R[1],g.d[i].R[2],SelectString(comment,"","		<!-- Rotation and translation vector -->")
		out += str
		sprintf str, "			<P unit=\"mm\">%.3f %.3f %.3f</P>\n",(g.d[i].P[0]/1000),(g.d[i].P[1]/1000),(g.d[i].P[2]/1000)
		out += str
		if (strlen(g.d[i].timeMeasured))
			sprintf str, "			<timeMeasured>%s</timeMeasured>%s\n",g.d[i].timeMeasured,SelectString(comment,"","	<!-- when this geometry was calculated -->")
			out += str
		endif
		if (strlen(g.d[i].geoNote))
			sprintf str, "			<note>%s</note>\n",g.d[i].geoNote
			out += str
		endif
		if (strlen(g.d[i].detectorID))
			sprintf str, "			<ID>%s</ID>%s\n",g.d[i].detectorID,SelectString(comment,"","				<!-- unique detector ID -->")
			out += str
		endif
		if (strlen(g.d[i].distortionMapFile))
			sprintf str, "			<distortionMap>%s</distortionMap>%s\n",g.d[i].distortionMapFile,SelectString(comment,"","				<!-- file with distortion map -->")
			out += str
		endif
		comment = 0
		out += "		</Detector>\n"
	endfor
	out += "	</Detectors>\n"

	if (!WireBad(g.wire))
		out += "\n	<Wire>\n"
		sprintf str,"		<dia unit=\"micron\">%.5g</dia>\n",g.wire.dia
		out += str
		sprintf str, "		<Knife>%g</Knife>				<!-- true if wire on a knife edge, false for free-standing wire -->\n",g.wire.knife
		out += str
		sprintf str, "		<Origin unit=\"micron\">%.8g %.8g %.8g</Origin>		<!-- wire origin in raw PM500 frame (micron) -->\n",g.wire.origin[0],g.wire.origin[1],g.wire.origin[2]
		out += str
		if (g.wire.Rmag>0)
			sprintf str,"		<R unit=\"radian\">%.9g %.9g %.9g</R>\n",g.wire.R[0],g.wire.R[1],g.wire.R[2]
			out += str
		endif
		sprintf str, "		<Axis>%.8g %.8g %.8g</Axis>		<!-- unit vector along wire axis, usually close to (1,0,0) -->\n",g.wire.axis[0],g.wire.axis[1],g.wire.axis[2]
		out += str
		if (!numtype(g.wire.F))
			sprintf str, "		<F unit=\"micron\">%.8g</F>		<!-- F of wire for a constant F wire scan (raw PM500 units) -->\n",g.wire.F
			out += str
		endif
		out += "	</Wire>\n"
	endif
	out += "</geoN>\n"
	return out
End
//// write the geometry structure to a keyword file
//Static Function WriteGeoToXMLfile(fileName,path,g,fileNote)
//	String fileName						// full path name to the file
//	String path							// name of an Igor path to use
//	STRUCT microGeometry &g			// the structure to fill from the file
//	String fileNote						// a note about this file
//
//	GeometryUpdateCalc(g)				// calculate other values
//	Variable f								// file id
//	if (strlen(path)<1 || strlen(fileName)<1)
//		Open/D/C="R*ch"/M="new geometry parameters file"/T=".xml" f as fileName
//		fileName = S_fileName
//	endif
//	Open/C="R*ch"/M="new geometry parameters file"/P=$path/T=".xml"/Z f as fileName
//	fileName = S_fileName
//	if (V_flag)
//		DoAlert 0, "nothing written to file"
//		return 1
//	endif
//	if( (strlen(path)+strlen(fileNote))==0 )
//		Prompt fileNote, "add a descriptive note to this file?"
//		DoPrompt "descriptive note",fileNote
//		if (V_flag)
//			fileNote=""
//		endif
//	endif
//
//	fprintf f, "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n\n"
//	fprintf f, "<geoN xmlns=\"http://sector34.xray.aps.anl.gov/34ide/geoN\">\n"
//	fprintf f, "	<dateWritten>%s</dateWritten>\n", date()
//	fprintf f, "	<timeWritten>%s (%g)</timeWritten>\n", Secs2Time(DateTime,3,1),date2secs(-1,-1,-1)/3600
//	fprintf f, "	<EPOCH start=\"midnight Jan 1, 1904\" unit=\"sec\">%.0f</EPOCH>\n", DateTime
//	if (strlen(fileNote))
//		fprintf f, "	<fileNote>%s</fileNote>\n",ReplaceString("\r",fileNote,"\\r")
//	endif
//
//	if (!SampleBad(g.s))
//		fprintf f, "\n	<Sample>\n"
//		fprintf f, "		<Origin unit=\"micron\">%.8g %.8g %.8g</Origin>	<!-- sample origin in raw PM500 units (micron) -->\n",g.s.O[0],g.s.O[1],g.s.O[2]
//		fprintf f, "		<R unit=\"radian\">%.9g %.9g %.9g</R>\n",g.s.R[0],g.s.R[1],g.s.R[2]
//		fprintf f, "	</Sample>\n"
//	endif
//
//	fprintf f, "\n	<Detectors Ndetectors=\"%g\">			<!-- Ndetectors, number of detectors in use, must be <= MAX_Ndetectors -->\n",g.Ndetectors
//	Variable i, comment=1
//	String pre
//	for (i=0;i<MAX_Ndetectors;i+=1)
//		if (!(g.d[i].used))
//			continue
//		endif
//		fprintf f, "		<Detector N=\"%d\">\n",i
//		fprintf f, "			<Npixels>%d %d</Npixels>%s\n",g.d[i].Nx,g.d[i].Ny,SelectString(comment,"","		<!-- Nx,Ny is number of un-binned pixels in full detector -->")
//		fprintf f, "			<size unit=\"mm\">%.7g %.7g</size>%s\n",(g.d[i].sizeX)/1000,(g.d[i].sizeY)/1000,SelectString(comment,"","	<!-- sizeX,sizeY otuside size of full detector -->")
//		fprintf f, "			<R unit=\"radian\">%.9g %.9g %.9g</R>%s\n",g.d[i].R[0],g.d[i].R[1],g.d[i].R[2],SelectString(comment,"","		<!-- Rotation and translation vector -->")
//		fprintf f, "			<P unit=\"mm\">%.3f %.3f %.3f</P>\n",(g.d[i].P[0]/1000),(g.d[i].P[1]/1000),(g.d[i].P[2]/1000)
//		if (strlen(g.d[i].timeMeasured))
//			fprintf f, "			<timeMeasured>%s</timeMeasured>%s\n",g.d[i].timeMeasured,SelectString(comment,"","	<!-- when this geometry was calculated -->")
//		endif
//		if (strlen(g.d[i].geoNote))
//			fprintf f, "			<note>%s</note>\n",g.d[i].geoNote
//		endif
//		if (strlen(g.d[i].detectorID))
//			fprintf f, "			<ID>%s</ID>%s\n",g.d[i].detectorID,SelectString(comment,"","				<!-- unique detector ID -->")
//		endif
//		if (strlen(g.d[i].distortionMapFile))
//			fprintf f, "			<distortionMap>%s</distortionMap>%s\n",g.d[i].distortionMapFile,SelectString(comment,"","				<!-- file with distortion map -->")
//		endif
//		comment = 0
//		fprintf f, "		</Detector>\n"
//	endfor
//	fprintf f, "	</Detectors>\n"
//
//	if (!WireBad(g.wire))
//		fprintf f, "\n	<Wire>\n"
//		fprintf f,"		<dia unit=\"micron\">%.5g</dia>\n",g.wire.dia
//		fprintf f, "		<Knife>%g</Knife>				<!-- true if wire on a knife edge, false for free-standing wire -->\n",g.wire.knife
//		fprintf f, "		<Origin unit=\"micron\">%.8g %.8g %.8g</Origin>		<!-- wire origin in raw PM500 frame (micron) -->\n",g.wire.origin[0],g.wire.origin[1],g.wire.origin[2]
//		if (g.wire.Rmag>0)
//			fprintf f,"		<R unit=\"radian\">%.9g %.9g %.9g</R>\n",g.wire.R[0],g.wire.R[1],g.wire.R[2]
//		endif
//		fprintf f, "		<Axis>%.8g %.8g %.8g</Axis>		<!-- unit vector along wire axis, usually close to (1,0,0) -->\n",g.wire.axis[0],g.wire.axis[1],g.wire.axis[2]
//		if (!numtype(g.wire.F))
//			fprintf f, "		<F unit=\"micron\">%.8g</F>		<!-- F of wire for a constant F wire scan (raw PM500 units) -->\n",g.wire.F
//		endif
//		fprintf f, "	</Wire>\n"
//	endif
//	fprintf f, "</geoN>\n"
//	Close f
//	String str
//	sprintf str, "wrote geometry to file '%s'",fileName
//	print str
//	DoAlert 0,str
//	return 0
//End

Static Function ReadGeoFromXMLfile(fileName,path,g)
	String fileName							// full path name to the file
	String path								// name of an Igor path to use
	STRUCT microGeometry &g				// the structure to fill from the file

	Variable f
	Open/M="geoN xml file"/P=$path/R/T=".xml"/Z=2 f as fileName
	if (strlen(S_fileName)<1 || !f)
		return 1
	endif
	fileName = S_fileName
	FStatus f
	String buf=""
	buf = PadString(buf,V_logEOF,0x20)
	FBinRead f, buf
	Close f

	buf = xmlContents(buf)
	Variable err = GeoFromXML(buf,g)
	return err
End
//
Static Function/S xmlContents(buf)
	String buf

	Variable i1, i0=strsearch(buf,"<?xml",0)	// find start of header tag
	if (i0<0)
		return ""
	endif
	i0 = strsearch(buf,"?>",0)					// find end of header tag
	if (i0<0)
		return ""
	endif
	buf = buf[i0+2,Inf]

	return buf
End
//
Static Function GeoFromXML(buf,g)
	String buf								// contents of the <geoN> xml
	STRUCT microGeometry &g				// the structure to fill from the file

	buf = XMLremoveComments(buf)			// remove all of the comments

	String str, nodes = XMLNodeList(buf)
	if (WhichListItem("geoN",nodes)<0)
		print "cannot find geoN"
		return 1
	endif
	buf = XMLtagContents("geoN",buf)
	if (strlen(buf)<1)
		print "cannot find <geoN> in xml"
		return 1
	endif

	Variable Ndetectors=NumberByKey("Ndetectors",XMLattibutes2KeyList("Detectors",buf),"=")
	Ndetectors = Ndetectors>0 ? round(Ndetectors) : 0
	if (Ndetectors<=0  || Ndetectors>MAX_Ndetectors)
		sprintf str,"This file contains %g, but only %g are allowed",Ndetectors,MAX_Ndetectors
		print str
		DoAlert 0,str
		return 1
	endif
	g.Ndetectors = Ndetectors
	Variable i
	for(i=0;i<MAX_Ndetectors;i+=1)
		g.d[i].used = 0
	endfor

	String dateWritten = XMLtagContents("dateWritten",buf)
	String timeWritten = XMLtagContents("timeWritten",buf)
	String fileNote = XMLtagContents("fileNote",buf)

	String list
	String sample = XMLtagContents("Sample",buf)
	if (strlen(sample))
		list = XMLtagContents2List("R",sample)
		g.s.R[0] = str2num(StringFromList(0,list))
		g.s.R[1] = str2num(StringFromList(1,list))
		g.s.R[2] = str2num(StringFromList(2,list))
		list = XMLtagContents2List("Origin",sample)
		g.s.O[0] = str2num(StringFromList(0,list))
		g.s.O[1] = str2num(StringFromList(1,list))
		g.s.O[2] = str2num(StringFromList(2,list))
	endif

	String wire = XMLtagContents("Wire",buf)
	if (strlen(wire))
		list = XMLtagContents2List("Origin",wire)
		g.wire.origin[0] = str2num(StringFromList(0,list))
		g.wire.origin[1] = str2num(StringFromList(1,list))
		g.wire.origin[2] = str2num(StringFromList(2,list))

		list = XMLtagContents2List("R",wire)
		g.wire.R[0] = str2num(StringFromList(0,list))
		g.wire.R[1] = str2num(StringFromList(1,list))
		g.wire.R[2] = str2num(StringFromList(2,list))

		list = XMLtagContents2List("Axis",wire)
		g.wire.axis[0] = str2num(StringFromList(0,list))
		g.wire.axis[1] = str2num(StringFromList(1,list))
		g.wire.axis[2] = str2num(StringFromList(2,list))

		g.wire.dia = str2num(XMLtagContents("dia",wire))
		g.wire.knife = str2num(XMLtagContents("Knife",wire))
		g.wire.F = str2num(XMLtagContents("F",wire))
	endif

	Variable N
	String detector, detectors=XMLtagContents("Detectors",buf)
	for (i=0;i<Ndetectors;i+=1)
		N = NumberByKey("N",XMLattibutes2KeyList("Detector",detectors,occurance=i),"=")
		if (N>=MAX_Ndetectors)
			sprintf str,"This file contains detector number %g, but max value is only %g",N,MAX_Ndetectors-1
			print str
			DoAlert 0,str
			continue
		endif
		detector = XMLtagContents("Detector",buf,occurance=i)

		g.d[N].used = 1
		g.d[N].timeMeasured = XMLtagContents("timeMeasured",detector)
		g.d[N].geoNote = XMLtagContents("note",detector)
		g.d[N].detectorID = XMLtagContents("ID",detector)
		g.d[N].distortionMapFile = XMLtagContents("distortionMap",detector)

		list =  XMLtagContents2List("Npixels",detector)
		g.d[N].Nx = str2num(StringFromList(0,list))
		g.d[N].Ny = str2num(StringFromList(1,list))

		list =  XMLtagContents2List("size",detector)
		g.d[N].sizeX = str2num(StringFromList(0,list))*1e3	// file uses mm, I need µm
		g.d[N].sizeY = str2num(StringFromList(1,list))*1e3

		list =  XMLtagContents2List("R",detector)
		g.d[N].R[0] = str2num(StringFromList(0,list))
		g.d[N].R[1] = str2num(StringFromList(1,list))
		g.d[N].R[2] = str2num(StringFromList(2,list))

		list =  XMLtagContents2List("P",detector)
		g.d[N].P[0] = str2num(StringFromList(0,list))*1e3	// file uses mm, I need µm
		g.d[N].P[1] = str2num(StringFromList(1,list))*1e3
		g.d[N].P[2] = str2num(StringFromList(2,list))*1e3
	endfor

	GeometryUpdateCalc(g)							// calculate other values
	printf "   geometry file was set on  %s, %s",dateWritten,timeWritten
	if (strlen(fileNote))
		printf ",   file note='%s'",fileNote
	endif
	print " "
	return 0
End
//
//Function testReadWriteXML()
//	STRUCT microGeometry g				// the structure to fill from the file
//	microGeo#ReadGeoFromXMLfile("geo.xml","home",g)
////	printGeometry(g)
//	Variable err = microGeo#WriteGeoToXMLfile("test.xml","home",g,"this is just a test write")
//	print "err = ",err
//	microGeo#ReadGeoFromXMLfile("test.xml","home",g)
////	printGeometry(g)
//End



Function SetDefaultGeo2Reference()							// set default geometry to the reference values
	STRUCT microGeometry g
	GeoReferenceOrientation(g)								// set geometry to reference values

	String/G root:Packages:geometry:geoStructStr=""		// use default location
	SVAR geoStructStr = root:Packages:geometry:geoStructStr
	StructPut/S/B=2 g, geoStructStr
	if (exists("root:Packages:geometry:PanelValues:Ndetectors")==2 && exists("root:Packages:geometry:PanelValues:knife")==2)
		SetGeometryPanelGlobals(g)							// set value on geoPanel if it is up
		DoUpdate
		if (WinType("microPanel#geoPanel")==7)
			GeoPanelDetectorDisable("microPanel#geoPanel")
			GeoPanelDirtyUpdate("microPanel#geoPanel",1)
		endif
	endif
End

Function GeoReferenceOrientation(g[,simple])				// sets g to the reference orientation (sort of an ideal set of values)
	STRUCT microGeometry &g
	Variable simple											// true means all positioners are square to beam,  false means beam tilted by mirros 6mrad
	g.Ndetectors = 3											// defining 3 detectors
	simple = ParamIsDefault(simple) ? 0 : simple

	// define Detector 0, located 500mm directly above sample (Orange)
	g.d[0].used = 1
	g.d[0].Nx = 2048 ;			g.d[0].Ny = 2048			// number of un-binned pixels in whole detector
	g.d[0].sizeX = 409.6e3;		g.d[0].sizeY = 409.6e3	// outside size of detector (micron)
	Variable Rval = -2/3*PI/sqrt(3)
	g.d[0].R[0]=Rval;			g.d[0].R[1]=Rval;		g.d[0].R[2]=Rval			// angle of detector, theta = -120¡ about (111)
	g.d[0].P[0]=25e3;			g.d[0].P[1]=0;			g.d[0].P[2]=510e3		// offset to detector (micron)
	g.d[0].timeMeasured = "Dec 4, 2008, 3:33pm"
	g.d[0].geoNote = "reference orientation"
	g.d[0].detectorID = "PE1621 723-3335"
	g.d[0].distortionMapFile = ""

	// define Detector 1, located ~400mm from sample, out along +X direction and up from horizontal by 45¡ (Yellow)
	g.d[1].used = 1
	g.d[1].Nx = 1024 ;			g.d[1].Ny = 1024			// number of un-binned pixels in whole detector
	g.d[1].sizeX = 204.8e3;		g.d[1].sizeY = 204.8e3	// outside size of detector (micron)
	g.d[1].R[0] = -1.77549569095761											// angle of detector
	g.d[1].R[1] = -0.735434395129632
	g.d[1].R[2] = -1.74815646188668
//	g.d[1].P[0]=-187e3;		g.d[1].P[1]=0;			g.d[1].P[2]=400e3	// offset to detector (micron)
	g.d[1].P[0]=-144e3;		g.d[1].P[1]=-1e3;		g.d[1].P[2]=410e3	// offset to detector (micron)
	g.d[1].timeMeasured = "Dec 5, 2008, 11:00am"
	g.d[1].geoNote = "reference orientation"
	g.d[1].detectorID = "PE0820 763-1807"
	g.d[1].distortionMapFile = ""

	// define Detector 2, located ~400mm from sample, out along -X direction and up from horizontal by 45¡ (Purple)
	g.d[2].used = 1
	g.d[2].Nx = 1024 ;			g.d[2].Ny = 1024			// number of un-binned pixels in whole detector
	g.d[2].sizeX = 204.8e3;		g.d[2].sizeY = 204.8e3	// outside size of detector (micron)
 	g.d[2].R[0] = -0.619944391960603		 			// angle of detector
 	g.d[2].R[1] = -1.49667815898843
 	g.d[2].R[2] = -0.610398437087624
	g.d[2].P[0]=-187e3;		g.d[2].P[1]=0;			g.d[2].P[2]=400e3	// offset to detector (micron)
	g.d[2].timeMeasured = "Dec 5, 2008, 11:00am"
	g.d[2].geoNote = "reference orientation"
	g.d[2].detectorID = "PE0820 763-1850"
	g.d[2].distortionMapFile = ""

	// define Wire
	g.wire.origin[0] = 0;			g.wire.origin[1] = 0;		g.wire.origin[2] = 0
	g.wire.dia = 52.
	g.wire.knife = 1											// true if wire on  a knife edge, false for free-standing wire
	g.wire.axis[0] = 1;			g.wire.axis[1] = 0;		g.wire.axis[2] = 0
	g.wire.R[0] = -0.006;		g.wire.R[1] = 0.006;		g.wire.R[2] = -1.8e-5	// 6mrad in X and Y (from the mirrors)
	g.wire.F = 3200											// this should be ~200µm above sample surface

	// define Sample
	g.s.O[0] = 0;					g.s.O[1] = 0;				g.s.O[2] = 0				// sample origin (micron)
	g.s.R[0] = -0.006;			g.s.R[1] = 0.006;			g.s.R[2] = -1.8e-5		// 6mrad in X and Y (from the mirrors)

	if (simple)
		g.wire.R[0] = 0;			g.wire.R[1] = 0;			g.wire.R[2] = 0			// no rotation
		g.s.R[0] = 0;				g.s.R[1] = 0;				g.s.R[2] = 0
	endif

	GeometryUpdateCalc(g)									// tidy up and do all pre-calculations
End
//
//	// detector rotation
//	Variable c135=-1/sqrt(2), s135=-1/sqrt(2)
//	Variable c90=0, s90=1, c45=1/sqrt(2), s45=1/sqrt(2)
//	if (dNum==0)										// detector 0
//		Make/N=(3,3)/O/D x135						// rotate around x-axis by -135¡
//		x135[0][0]=1;		x135[0][1]=0;			x135[0][2]=0
//		x135[1][0]=0;		x135[1][1]=c135;		x135[1][2]=-s135
//		x135[2][0]=0;		x135[2][1]=s135;		x135[2][2]=c135
//		MatrixOp/O dRot = x135						// detector rotation matrix
//		KillWaves/Z x135
//	elseif (dNum==1)									// detector 1
//		Make/N=(3,3)/D/O y90, z45					// z45 x y90 is a sort of random rotation for the lattice
//		y90[0][0]=c90;	y90[0][1]=0;			y90[0][2]=s90
//		y90[1][0]=0;		y90[1][1]=1;			y90[1][2]=0
//		y90[2][0]=-s90;	y90[2][1]=0;			y90[2][2]=c90
//		z45[0][0]=c45;	z45[0][1]=-s45;		z45[0][2]=0
//		z45[1][0]=s45;	z45[1][1]=c45;		z45[1][2]=0
//		z45[2][0]=0;		z45[2][1]=0;			z45[2][2]=1
//		MatrixOp/O dRot = z45 x y90					// detector rotation matrix
//		KillWaves/Z y90,z45
//	elseif (dNum==2)									// detector 2
//		Make/N=(3,3)/D/O y90, z45					// z45 x y90 is a sort of random rotation for the lattice
//		s90 *= -1										// reverse rotation direction for both y90 and z45
//		s45 *= -1
//		y90[0][0]=c90;	y90[0][1]=0;			y90[0][2]=s90
//		y90[1][0]=0;		y90[1][1]=1;			y90[1][2]=0
//		y90[2][0]=-s90;	y90[2][1]=0;			y90[2][2]=c90
//		z45[0][0]=c45;	z45[0][1]=-s45;		z45[0][2]=0
//		z45[1][0]=s45;	z45[1][1]=c45;		z45[1][2]=0
//		z45[2][0]=0;		z45[2][1]=0;			z45[2][2]=1
//		MatrixOp/O dRot = z45 x y90					// detector rotation matrix
//		KillWaves/Z y90,z45
//	endif

// ====================================== End of Geometry =====================================
//=======================================================================================



//=======================================================================================
// ======================================= Start of XHF ======================================

// Beam line corodinate system (X=out door, Y=up, Z=downstream), theta is the angle between Z and H
// definition of of angle: H wire scan is horiz at angle=0,   at angle=0,  F=-Y, H=+Z

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
		H = round(H*1e3)/1e3				// trim to 3 places
		F = round(F*1e3)/1e3
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
		H = round(H*1e3)/1e3				// trim to 3 places
		F = round(F*1e3)/1e3
		list = ReplaceNumberByKey("H2",list,H,"=")
		list = ReplaceNumberByKey("F2",list,F,"=")
	endif

	return list
End

// ======================================== End of XHF ======================================
//=======================================================================================



//=======================================================================================
// =========================== Start of Sample & Wire Positioner Correction ===========================

// Sample Position
Function PM500X1toBeamLineX1(s,XYZ1)		// change from Sample PM500 numbers to Beam Line coords
	STRUCT sampleGeometry &s
	Wave XYZ1									// a 3-vector, PM500 coords on input, change to BeamLine on output

	Variable X1, Y1, Z1
	X1 = X1correctedKeyence(XYZ1[0])			// first correct of errors in the encoder
	Y1 = Y1correctedKeyence(XYZ1[1])
	Z1 = Z1correctedKeyence(XYZ1[2])

	X1 -= X1correctedKeyence(s.O[0])			// translate to origin, subtract off sample origin
	Y1 -= Y1correctedKeyence(s.O[1])	
	Z1 -= Z1correctedKeyence(s.O[2])	

	XYZ1[0] = s.R00*X1 + s.R01*Y1 + s.R02*Z1	// XYZ1 = s.Rij x {X1,Y1,Z1},   rotate by R (a small rotation)
	XYZ1[1] = s.R10*X1 + s.R11*Y1 + s.R12*Z1
	XYZ1[2] = s.R20*X1 + s.R21*Y1 + s.R22*Z1
End


// Wire Position
Function PM500X2toBeamLineX2(w,XYZ2)		// change Wire position from raw PM500 numbers to Beam Line coords
	STRUCT wireGeometry &w
	Wave XYZ2									// a 3-vector, PM500 coords on input, change to BeamLine on output

	Variable X2, Y2, Z2
	X2 = X2correctedKeyence(XYZ2[0])			// first correct of errors in the encoder
	Y2 = Y2correctedKeyence(XYZ2[1])
	Z2 = Z2correctedKeyence(XYZ2[2])

	X2 -= w.origin[0]							// translate to origin (distance of wire from origin)
	Y2 -= w.origin[1]
	Z2 -= w.origin[2]

	XYZ2[0] = w.R00*X2 + w.R01*Y2 + w.R02*Z2	// XYZ2 = w.Rij x {X2,Y2,Z2},   rotate by R (a small rotation)
	XYZ2[1] = w.R10*X2 + w.R11*Y2 + w.R12*Z2
	XYZ2[2] = w.R20*X2 + w.R21*Y2 + w.R22*Z2
End

// ==================================== Keyence Correction =====================================

//	Make/O root:Packages:micro:X1correctionWave={0.31,0.32,0.2,0.09,-0.1,-0.3,-0.44,-0.47,-0.42,-0.37,-0.32,-0.19,-0.12,-0.07,0.01,0.17,0.35,0.52,0.67,0.73,0.69,0.53,0.32,0.13,-0.01,-0.13,-0.18,-0.18,-0.18,-0.16,-0.26,-0.25,-0.29,-0.38,-0.41,-0.38,-0.33,-0.18,0.01,0.23,0.42}
//	Make/O root:Packages:micro:Y1correctionWave={0.37,0.38,0.35,0.26,0.1,-0.06,-0.22,-0.31,-0.36,-0.35,-0.62,-0.51,-0.39,-0.28,-0.19,-0.03,0.11,0.25,0.34,0.38,0.38,0.34,0.28,0.17,0.04,-0.09,-0.22,-0.33,-0.4,-0.41,-0.39,-0.36,-0.28,-0.19,-0.08,0.06,0.2,0.3,0.4,0.54,0.66}
//	Make/O root:Packages:micro:Z1correctionWave={-0.09,-0.17,-0.21,-0.24,-0.25,-0.22,-0.14,-0.01,0.16,0.28,0.34,0.38,0.37,0.33,0.3,0.27,0.23,0.1,-0.05,-0.22,-0.18,-0.24,-0.28,-0.29,-0.28,-0.25,-0.18,-0.09,0,0.1,0.28,0.34,0.33,0.32,0.32,0.32,0.31,0.26,0.19,0.09,-0.02}
//	SetScale/P x 0,1,"µm", root:Packages:micro:X1correctionWave,root:Packages:micro:Y1correctionWave,root:Packages:micro:Z1correctionWave
//	SetScale d 0,0,"µm", root:Packages:micro:X1correctionWave,root:Packages:micro:Y1correctionWave,root:Packages:micro:Z1correctionWave
//
Static Function X1correctedKeyence(X1)		// takes PM500 X1 and returns the "real" X1
	Variable X1									// input the PM500 X1 reading
	Wave deltaW=root:Packages:micro:X1correctionWave
	Variable zWidth = DimDelta(deltaW,0)*(numpnts(deltaW)-1)	// width of the correction wave (µm)
	Variable dX1
	dX1 = mod(X1,zWidth)
	dX1 = dX1 < 0 ? dX1+zWidth : dX1
	return deltaW(dX1)+X1						// return corrected X1, this should be the true X1 value
End
//
Static Function Y1correctedKeyence(Y1)		// takes PM500 Y1 and returns the "real" Y1
	Variable Y1									// input the PM500 Y1 reading
	Wave deltaW=root:Packages:micro:Y1correctionWave
	Variable zWidth = DimDelta(deltaW,0)*(numpnts(deltaW)-1)	// width of the correction wave (µm)
	Variable dY1
	dY1 = mod(Y1,zWidth)
	dY1 = dY1 < 0 ? dY1+zWidth : dY1
	return deltaW(dY1)+Y1						// return corrected Y1, this should be the true Y1 value
End
//
Static Function Z1correctedKeyence(Z1)		// takes PM500 Z1 and returns the "real" Z1
	Variable Z1									// input the PM500 Z1 reading
	Wave deltaW=root:Packages:micro:Z1correctionWave
	Variable zWidth = DimDelta(deltaW,0)*(numpnts(deltaW)-1)	// width of the correction wave (µm)
	Variable dZ1
	dZ1 = mod(Z1,zWidth)
	dZ1 = dZ1 < 0 ? dZ1+zWidth : dZ1
	return deltaW(dZ1)+Z1						// return corrected Z1, this should be the true Z1 value
End


Static Function X2correctedKeyence(X2)		// takes PM500 X2 and returns the "real" X2	****	DUMMY	DUMMY	DUMMY	****
	Variable X2									// input the PM500 X2 reading
	return X2										// return corrected X2, this should be the true X2 value
End
Static Function Y2correctedKeyence(Y2)		// takes PM500 Y2 and returns the "real" Y2	****	DUMMY	DUMMY	DUMMY	****
	Variable Y2									// input the PM500 Y2 reading
	return Y2										// return corrected Y2, this should be the true Y2 value
End
Static Function Z2correctedKeyence(Z2)		// takes PM500 Z2 and returns the "real" Z2	****	DUMMY	DUMMY	DUMMY	****
	Variable Z2									// input the PM500 Z2 reading
	return Z2										// return corrected Z2, this should be the true Z2 value
End

// ============================ End of Sample & Wire Positioner Correction ===========================
//=======================================================================================



//=======================================================================================
// ===================================== Start of Distortion ====================================

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


////	xymap,Nx,Ny,cornerX0,cornerY0,cornerX1,cornerY1
//Function/S loadCCDCorrectionTable(mName)
//	String mName						// name for xymap file
//	String CCDxyMapfile = "CCD_distorMay03_corr.dat"
//	mName = SelectString(strlen(mName)>0,"xymap",CleanupName(mName,0))
//
//	if (ItemsInList(GetRTStackInfo(0))<2)
//		print "Loading the correction table ..."
//	endif
//	Variable refNum
//	Open/M="distortion table"/P=home/R/T="TEXT"/Z=2 refNum as CCDxyMapfile
//	if (V_flag)
//		return ""
//	endif
//	String line
//	FReadLine refNum, line
//
//	Variable Nx,Ny,cornerX0,cornerY0,cornerX1,cornerY1
//	sscanf line, "%d %d %d %d %d %d",Nx,Ny,cornerX0,cornerY0,cornerX1,cornerY1
//	if (V_flag!=6)
//		Close refNum
//		return ""
//	endif
//
//	String noteStr=""
//	noteStr = ReplaceNumberByKey("cornerX0",noteStr,cornerX0,"=")
//	noteStr = ReplaceNumberByKey("cornerY0",noteStr,cornerY0,"=")
//	noteStr = ReplaceNumberByKey("cornerX1",noteStr,cornerX1,"=")
//	noteStr = ReplaceNumberByKey("cornerY1",noteStr,cornerY1,"=")
//	noteStr = ReplaceStringByKey("CCDname",noteStr,"White2084x2084","=")
//
//	Make/N=(Nx,Ny,4)/O $mName
//	Wave xymap = $mName
//	Note/K xymap
//	Note xymap, noteStr
//
//	Variable i,j,x1,y1,x2,y2
//	for (j=0;j<Ny;j+=1)
//		for (i=0;i<Ny;i+=1)
//			FReadLine refNum, line
//			sscanf line, "%g %g %g %g",x1,y1,x2,y2
//			if (V_flag!=4)
//				Close refNum
//				return ""
//			endif
//			xymap[i][j][0] = x1
//			xymap[i][j][1] = y1
//			xymap[i][j][2] = x2
//			xymap[i][j][3] = y2
//		endfor
//	endfor
//	Close refNum
//	return GetWavesDataFolder(xymap,2)
//End
////
//Function LoadStandardDistortion()
//	String str= FunctionPath("LoadStandardDistortion")
//	if (strlen(str)<2)
//		return 1
//	endif
//	String fileName = ParseFilePath(1, str, ":", 1, 0)+"xymap.itx"
//	LoadWave/O/T fileName
//	if (V_flag!=1)
//		DoAlert 0, "Unable to load xymap.itx"
//		return 1
//	elseif (!stringmatch(StringFromList(0,S_waveNames),"xymap"))
//		DoAlert 0,"Loaded something other than xyamp, see history"
//		return 1
//	endif
//
//	// move to root:Packages:geometry if not already there
//	if (stringmatch(GetDataFolder(1),"root:Packages:geometry:"))
//		return 0
//	endif
//	Duplicate/O xymap root:Packages:geometry:xymap
//	KillWaves xymap
//	printf "¥¥loaded standard distortion from '%s' into file root:Packages:geometry:xymap\r",fileName
//	return 0
//End


//Function MakeDistortionMap()
//	Variable N=20
//	Make/N=(N,N)/O totalDistortion
//	Make/N=(N,N)/O/C xyDistortion
//	SetScale/I x 0,2083,"pixels", totalDistortion
//	SetScale/I y 0,2083,"pixels", totalDistortion
//	SetScale d 0,0,"pixels", totalDistortion
//	Variable x0=DimOffset(totalDistortion,0), dx=DimDelta(totalDistortion,0)
//	Variable y0=DimOffset(totalDistortion,1), dy=DimDelta(totalDistortion,1)
//
//	Wave xymap=root:Packages:geometry:xymap
//	Variable px,py, i,j
//	for (j=0;j<N;j+=1)
//		py = y0 + j*dy
//		for (i=0;i<N;i+=1)
//			px = x0 + i*dx
//			xyDistortion[i][j] = peakcorrection2(xymap,px,py)		// returns cmplx(dx,dy)
//		endfor
//	endfor
//	totalDistortion = cabs(xyDistortion[p][q])
//
//	if (strlen(WinList("GraphDistortionMap", "","WIN:1"))>1)
//		DoWindow/F GraphDistortionMap
//		SetDrawLayer /K UserFront	
//	else
//		Display/K=1/W=(168,55,798,615)
//		DoWindow/C GraphDistortionMap
//		AppendImage totalDistortion
//		ModifyImage totalDistortion ctab= {0,7.1,Terrain,1}
//		ModifyGraph margin(right)=72,width={Aspect,1},mirror=2
//		ModifyGraph lblMargin(left)=7,lblMargin(bottom)=5,axOffset(left)=-2.57143,lblLatPos(left)=-44,lblLatPos(bottom)=47
//		SetAxis/A/R left
//		ColorScale/N=text0/F=0/S=3/A=RC/X=-13.63/Y=0.40 image=totalDistortion
//	endif
//	xyDistortion *= 20							// scale up to make lines longer
//	for (j=0;j<N;j+=1)
//		py = y0 + j*dy
//		for (i=0;i<N;i+=1)
//			px = x0 + i*dx
//			SetDrawLayer UserFront
//			SetDrawEnv xcoord=bottom,ycoord=left,arrow=1, arrowlen=7.0,arrowfat=0.3	// draw the arrow
//			DrawLine px,py,px+real(xyDistortion[i][j]),py+imag(xyDistortion[i][j])
//		endfor
//	endfor
//	KillWaves/Z xyDistortion
//End


//Static Function resetCenterDistortionMap(xymap,xc,yc)	// set distortion correction of xc & yc into note of xymap, used so correction to (xc,yc) is zero
//	Wave xymap												// distortion map
//	Variable xc,yc												// values that should get zero distortion correction, should be geo.xcent,geo.ycent
//	if (!WaveExists(xymap))								// wave does not exist, do nothing
//		return 1
//	endif
//
//	String wnote = note(xymap)
//	wnote = ReplaceNumberByKey("dxCenter",wnote,0,"=")			// set correction to zero
//	wnote = ReplaceNumberByKey("dyCenter",wnote,0,"=")
//	Note/K xymap, wnote													// re-set value in wave note
//	if (numtype(xc+yc))
//		return 1															// invalid (xc,yc) passed so leave corrections at zero
//	endif
//
//	Variable px=xc-1, py=yc-1
//	Variable/C dp = peakcorrection2(xymap,px,py)					// returns cmplx(dx,dy), px,py are changed
//	wnote = ReplaceNumberByKey("dxCenter",wnote,real(dp),"=")	// add this to corrected to get proper value
//	wnote = ReplaceNumberByKey("dyCenter",wnote,imag(dp),"=")
//	Note/K xymap, wnote													// store values in wave note
//End

// ===================================== End of Distortion =====================================
//=======================================================================================



//=======================================================================================
//==================================== Start of Wire & Depth ===================================

// Given the three points, P=pixel, S=source, W=wire:
// The next three routines do calculations for the unknown point when the other two are known


// Given a ray starting from depth (on the beam), and ending at the point xyzPixel (probably at a pixel), finds the wire H to put wire
Function DepthPixel2WireH(g,depth,xyzPixel,edge)	// calc H of wire that is tangent to line from depth on beam to pixel on detector (in beam line coords)
	STRUCT microGeometry &g			// g.wire.F provides the F
	Variable depth							// distance along incident beam from origin (the Si positino) (µm)
	Wave xyzPixel						// location of pixel in beam line coordinates
	Variable edge							// 1=leading edge (usual),  0=trailing edge

	Wave Rvec = root:Packages:geometry:depth_Rvec	// rotation vector that puts axis of wire along {1,0,0}
	Wave rho = root:Packages:geometry:depth_rhoX		// rotate to frame where wire axis is along x-axis
	Wave pixel = root:Packages:geometry:depth_pixel	// position of pixel in rotated coords
	Wave S = root:Packages:geometry:depth_S			// position of source in rotated coords
	Wave sigma = root:Packages:geometry:depth_sigma	// normalized vector pointing from Pixel to Source (rotated coords)
	Wave delta = root:Packages:geometry:depth_delta	// scanned direction of wire motion in rotated coords
	Wave a = root:Packages:geometry:depth_a			// direction of wire axis in rotated coords
	Wave nhat = root:Packages:geometry:depth_nhat	// normal to plane
	Wave wo = root:Packages:geometry:depth_wo		// point on the plane of the wire center (beam line coords)
	Wave C = root:Packages:geometry:depth_C			// constant term in linear equation (beam line coords)
	Wave rhoW = root:Packages:geometry:depth_rhoW	// rotation of wire PM500
	Wave mat = root:Packages:geometry:depth_mat		// matrix sent to linear solver

	Variable R=(g.wire.dia)/2							// wire radius
	Variable dF = g.wire.F - YZ2F(g.wire.origin[1],g.wire.origin[2])
	Variable kix=0, kiy=0, kiz=1						// direction of incident beam in beam line coords (do not need kix)
	Variable cosTheta = NumVarOrDefault("root:Packages:geometry:cosThetaWire",cos(PI/4))
	Variable sinTheta = NumVarOrDefault("root:Packages:geometry:sinThetaWire",sin(PI/4))

	Rvec = {0,g.wire.axisR[2], -g.wire.axisR[1]}		// cross product,   axisR x {1,0,0}
	Variable theta = asin(normalize(Rvec))				// length of Rvec = sin(theta)
	Rvec *= theta
	rotationMatFromAxis(Rvec,NaN,rho)				// rotation matrix makes wire axis lie along {1,0,0}

	rhoW[0][0]=g.wire.R00;	rhoW[0][1]=g.wire.R01;	rhoW[0][2]=g.wire.R02
	rhoW[1][0]=g.wire.R10;	rhoW[1][1]=g.wire.R11;	rhoW[1][2]=g.wire.R12
	rhoW[2][0]=g.wire.R20;	rhoW[2][1]=g.wire.R21;	rhoW[2][2]=g.wire.R22

	MatrixOp/O pixel = rho x xyzPixel					// rotate pxiel to new coordinate system

//	S = { kix/kiz, kiy/kiz, 1 }							// source point for the ray (in beam line coords)
//	S *= depth
	S = {kix,kiy,kiz}										// source point for the ray (in beam line coords)
	normalize(S)
	S *= depth
	MatrixOp/O S = rho x S

	sigma = S - pixel										// direction of ray from Pixel to Source in rotated frame
	normalize(sigma)
	nhat = {0, sigma[2], -sigma[1]}					// perpendicular to sigma, pointing in +Z direction (rotated frame)
	normalize(nhat)

	delta[0] = {0, sinTheta, cosTheta}					// direction of wire motion (PM500 coords)
	MatrixOp/O delta = rho x rhoW x delta				// delta = wire.Rij x {0, sinTheta, cosTheta},   rotate delta to beam line coords, then to rotated frame

	a = g.wire.axis[p]
	MatrixOp/O a = rho x rhoW x a						// a = rho x wire.Rij x wire.axis,   should be {1,0,0}

	wo = {0, -dF*cosTheta, dF*sintheta}				// off set to plane with center of wire at dF from origin
	MatrixOp/O wo = rho x rhoW x wo

	edge = edge ? -1 : +1									// use - for leading edge, + for trailing edge
	C = S + edge*R*nhat - wo

	mat[][0] = sigma[p]									// Solve the equation: (S-P)*t + delta*u + axisR*v = C
	mat[][1] = delta[p]
	mat[][2] = a[p]
	MatrixLinearSolve/M=1/Z/O mat C
	KillWaves/Z W_IPIV
	//printf "{t,u,v} = {%g, %g, %g}\r",C[0],C[1],C[2]
	if (V_flag)
		return NaN
	endif
	return C[1]											// the H in PM500 coords
End


// Returns depth (starting point) of ray with one end point xyzPixel (probably a pixel location) that is tangent
// to leading edge of the wire and intersects the incident beam.  The returned depth is relative to the Si position (µm)
// depth is the measured along the incident beam from the origin (the Si position), not just the Z value
Function PixelxyzWire2depth(g,xyzPixel,Xw,edge)	// returns depth (µm)
	STRUCT microGeometry &g			// g.wire.F provides the F
	Wave xyzPixel						// location of pixel in beam line coordinates, origin at Si (µm)
	Wave Xw								// center of wire in raw PM500 coordinates (µm)
	Variable edge							// 1=leading edge (usual),  0=trailing edge

	Wave Rvec = root:Packages:geometry:depth_Rvec	// rotation vector that puts axis of wire along {1,0,0}
	Wave rho = root:Packages:geometry:depth_rhoX		// rotate to frame where wire axis is along x-axis
	Wave pixel = root:Packages:geometry:depth_pixel	// position of pixel (beam line coords)
	Wave wc = root:Packages:geometry:depth_wc		// position of wire center
	Wave ki = root:Packages:geometry:depth_ki			// incident beam direction (in rotated by rho)
	Wave S = root:Packages:geometry:depth_S			// position of source rotated (beam line coords)

	Rvec = {0,g.wire.axisR[2], -g.wire.axisR[1]}		// cross product,   axisR x {,0,0}
	Variable theta = asin(normalize(Rvec))				// length of Rvec = sin(theta)
	Rvec *= theta
	rotationMatFromAxis(Rvec,NaN,rho)				// rotation matrix makes wire axis lie along {1,0,0}

	wc = Xw
	PM500X2toBeamLineX2(g.wire,wc)					// transform Xw into Beam Line coords (subtract origin & rotate)
	MatrixOp/O wc = rho x wc							// rotate wire center to new coordinate
	MatrixOp/O pixel = rho x xyzPixel					// rotate pxiel to new coordinate system
	ki = {0,0,1}											// ki must be normalized
	MatrixOp/O ki = rho x ki								// ki in rotated frame

	Variable R=(g.wire.dia)/2							// wire radius
	Variable vy,vz, v										// v is vector from pixel to the wire center
	vy = wc[1] - pixel[1]
	vz = wc[2] - pixel[2]
	v = sqrt(vy*vy + vz*vz)								// length of vector (0,vy,vz) = wc[]-pixle[]

	Variable phi0 = atan2(vz,vy)						// angle from yhat to V (V is vector from pixel to wire center)
	Variable dphi = asin(R/v)							// angle between lines to wire center and tangent point of wire (from detector pixel)
	Variable tanphi = tan((edge ? phi0-dphi : phi0+dphi))	// phi is angle from yhat to V (measured at the pixel)

	// line from pixel to tangent point is:   z = y*tan(phio±dphi) + b
	Variable b = pixel[2] - pixel[1]*tanphi

	// line of incident beam is:   y = kiy/kiz * z		Thiis line goes through origin, so intercept is 0
	// find intersection of this line and line from pixel to tangent point
	S[2] = b / (1-tanphi*ki[1]/ki[2])					// intersection of two lines at this z value
	S[1] = ki[1]/ki[2] * S[2]							// corresponding y of point on incident beam
	S[0] = ki[0]/ki[2] * S[2]							// corresponding z of point on incident beam

	return MatrixDot(ki,S)								// depth measured along incident beam (remember that ki is normalized)
//	MatrixOp/O S = rho^t x S								// rotate back to beam line frame
//	return S[2]
End

// set mat to be a rotation matrix about axis with angle
Static Function rotationMatFromAxis(axis,angle,mat)
	Wave axis				// axis about which to rotate (or possibly Rodriques vector)
	Variable angle			// angle to rotate (degrees), assumes axis is true Rotation vector if angle invalid
	Wave mat				// desired rotation matrix

	Variable len = norm(axis)
	angle = numtype(angle) ? len : angle*PI/180	// the rotation angle (rad)

	if (angle==0)									// zero angle rotation is just the identity matrix
		mat = (p==q)
		return 0
	endif

	Variable nx=axis[0]/len, ny=axis[1]/len, nz=axis[2]/len
	Variable cosa=cos(angle), sina=sin(angle)
	Variable c1 = 1-cosa
	if (DimSize(mat,0)!=3 && DimSize(mat,1)!=3)
		Abort "in rotationMatFromAxis(), mat must be (3,3)"
//		DoAlert 0, "in rotationMatFromAxis(), mat must be (3,3)"
//		return 1
	endif
	// from		http://mathworld.wolfram.com/RodriguesRotationFormula.html (I double checked this too.)
	mat[0][0] = cosa+nx*nx*c1;			mat[0][1] =nx*ny*c1-nz*sina;			mat[0][2] = nx*nz*c1+ny*sina;
	mat[1][0] = nx*ny*c1+nz*sina;		mat[1][1] = cosa+ny*ny*c1;			mat[1][2] =ny*nz*c1-nx*sina;
	mat[2][0] = nx*nz*c1-ny*sina;		mat[2][1] = ny*nz*c1+nx*sina;		mat[2][2] = cosa+nz*nz*c1;
	return 0
End


//Menu "Test_Wire"
//	"Many Loop Test", testManyDepthWirePixel(NaN,NaN,simple=0)
//	"One Loop Test", test_Depth_Loop(NaN)
//	"Graph of error vs px,py,depth",GraphCorrelationErrors("")
//	SubMenu "2D error surface"
//		"py vs px",graph_px_py_err()
//		"depth vs px",graph_px_depth_err()
//		"depth vs py",graph_py_depth_err()
//	End
//End
////
//Function testManyDepthWirePixel(N,range[,simple,edge])
//	Variable N
//	Variable range
//	Variable simple										// true means all positioners are square to beam,  false means beam tilted by mirros 6mrad
//	Variable edge											// leading edge=1, trailing=0
//	simple = ParamIsDefault(simple) ? 0 : simple
//	edge = ParamIsDefault(edge) ? 1 : edge
//	range = (range>=0) ? range : 1000					// default to 1000 µm
//	String simpleList = "both mirrors tilt by 6mr;all square;use current default"
//	String edgeList = "trailing edge;leading edge"
//	if (!(N>0))
//		simple += 1
//		edge += 1
//		N = (N>0) ? N : 5000							// default to 5000
//		range = (range>0) ? range : 50					// default to ±50µm
//
//		Prompt N,"# of loops"
//		Prompt range "range of depths (µm)"
//		Prompt simple,"geometry",popup, simpleList
//		Prompt edge,"edge of wire",popup, edgeList
//		DoPrompt "test loops",N,range,simple,edge
//		if (V_flag)
//			return 1
//		endif
//		simple -= 1
//		edge -= 1
//	endif
//	printf "testManyDepthWirePixel(%g,%g, simple=%g, edge=%g)		// %s,   %s\r",N,range,simple,edge,StringFromList(simple,simpleList),StringFromList(edge,edgeList)
//
//	Variable i
//	SetRandomSeed 1										// use this if calling testOptimizeAll() repeatedly to test
//	Make/N=(N,4)/O manyErrs=NaN
//	SetDimLabel 1,0,depth,manyErrs
//	SetDimLabel 1,1,px,manyErrs
//	SetDimLabel 1,2,py,manyErrs
//	SetDimLabel 1,3,err,manyErrs
//	Variable depth,px,py,err,errMax=0, imax=-1
//	for (i=0;i<N;i+=1)
//		depth = enoise(range)
//		px = enoise(1022)+1023
//		py = enoise(1022)+1023
//		err = test_Depth_Loop(depth,px=px,py=py,simple=simple,edge=edge)
//		if (err>errMax)
//			errMax = err
//			imax = i
//		endif
//		manyErrs[i][0] = depth
//		manyErrs[i][1] = px
//		manyErrs[i][2] = py
//		manyErrs[i][3] = err
//	endfor
//	WaveStats/Q manyErrs
//
//	ImageStats/G={0,N-1,3,3} manyErrs
//	printf "worst point at i=%d,  pixel = (%g, %g),  depth = %g µm,  err = %g µm\r",imax,manyErrs[imax][1],manyErrs[imax][2],manyErrs[imax][0],manyErrs[imax][3]
//	printf "the error range is [%.2g, %.2g]µm, <error>=%gµm,    using a depth range of ±%gµm\r",V_min,V_max,V_avg,range
//End
////
//Function test_Depth_Loop(depth0[,simple,edge,px,py])
//	Variable depth0
//	Variable simple										// true means all positioners are square to beam,  false means beam tilted by mirros 6mrad
//	Variable edge											// leading edge=1, trailing=0
//	Variable px, py										// detector center at (1023.5, 1023.5)
//	simple = ParamIsDefault(simple) ? 0 : simple
//	edge = ParamIsDefault(edge) ? 1 : edge
//	px = ParamIsDefault(px) ? 1023.5 : px
//	py = ParamIsDefault(py) ? 1023.5 : py
//	depth0 = numtype(depth0) ? 0 : depth0
//	Variable printing = strlen(GetRTStackInfo(2))==0
//	if (printing)
//		simple += 1
//		edge += 1
//		String simpleList = "both mirrors tilt by 6mr;all square;use current default"
//		String edgeList = "trailing edge;leading edge"
//		Prompt depth0,"depth µm"
//		Prompt simple,"geometry",popup, simpleList
//		Prompt px,"X pixel"
//		Prompt py,"Y pixel"
//		Prompt edge,"edge of wire",popup, edgeList
//		DoPrompt "test",depth0,simple,edge,px,py
//		if (V_flag)
//			return NaN
//		endif
//		simple -= 1
//		edge -= 1
//		printf "starting with pixel = {%g, %g},  '%s',   '%s'\r",px,py,StringFromList(simple,simpleList),StringFromList(edge,edgeList)
//	endif
//
//	STRUCT microGeometry g								// returns 0 if something set, 0 is nothing done
//	if (simple==2)
//		FillGeometryStructDefault(g)
//	else
//		GeoReferenceOrientation(g,simple=simple)
//	endif
//
//	Make/N=3/O/D xyzPixel, Xw
//	pixel2XYZ(g.d[0],px,py,xyzPixel)
//	xyzPixel = abs(xyzPixel)<1e-9 ? 0 : xyzPixel
//
//	Variable H = DepthPixel2WireH(g,depth0,xyzPixel,edge)
//
//	Variable F = g.wire.F
//	Xw[0] = {0, HF2Y(H,F), HF2Z(H,F)}
//	Variable depth = PixelxyzWire2depth(g,xyzPixel,Xw,edge)	// returns depth (µm)
//
//	if (printing)
//		printf "pixel = {%g, %g, %g}\r",xyzPixel[0],xyzPixel[1],xyzPixel[2]
//		printf "H(depth=%.10gµm) = %.10g µm\r",depth0,H
//		printf "depth(H=%.10g, F=%.10g) = %.10g µm\r",H,F,depth
//		printf "Ædepth = %.2g µm (= %.2W1Pm)\r",depth-depth0,(depth-depth0)/1e6
//	endif
//	KillWaves/Z W_Cross, xyzPixel, Xw
//	return (depth-depth0)
//End
////
//Function GraphCorrelationErrors(column)
//	String column
//	column = LowerStr(column)
//
//	String colList = "depth;px;py"
//	if (strlen(column)<1)
//		Prompt column,"X-axis",popup, colList
//		Doprompt "Column",column
//		if (V_flag)
//			return 1
//		endif
//	endif
//	Variable ix = WhichListItem(column,colList)
//	if (ix<0)
//		return 1
//	endif
//	Display/K=1/W=(656,44,1205,387) manyErrs[*][3] vs manyErrs[*][ix]
//	ModifyGraph gfMult=110, mode=3, marker=19, msize=3, tick=2, mirror=1, minor=1
//	ModifyGraph zero(left)=2, lowTrip=0.001, standoff=0
//	if (ix==0)
//		ModifyGraph rgb=(2,39321,1)
//	elseif (ix==1)
//		ModifyGraph rgb=(65535,0,0)
//	elseif (ix==2)
//		ModifyGraph rgb=(0,0,65535)
//	endif
//	Label left "error";DelayUpdate
//	Label bottom column+SelectString(stringmatch(column,"depth"),""," (µm)")
//	ShowInfo
//End
////
//Window graph_px_py_err() : Graph
//	PauseUpdate; Silent 1		// building window...
//	Display/K=1/W=(859,44,1438,575) manyErrs[*][2] vs manyErrs[*][1]
//	ModifyGraph gfMult=110,width={Aspect,1}, mode=3, marker=19, msize=3
//	ModifyGraph zColor(manyErrs)={manyErrs[*][3],-2.5e-10,2.5e-10,RedWhiteBlue256}
//	ModifyGraph tick=2, mirror=1, minor=1, lowTrip=0.001, standoff=0
//	Label left "py"
//	Label bottom "px"
//	SetAxis left 0,2048
//	SetAxis bottom 0,2048
//	Cursor/P A manyErrs 69
//	ShowInfo
//	SetDrawLayer UserFront
//	SetDrawEnv xcoord= prel,ycoord= left,dash= 1
//	DrawLine 0,1023.5,1,1023.5
//	SetDrawEnv xcoord= bottom,ycoord= prel,dash= 1
//	DrawLine 1023.5,0,1023.5,1
//EndMacro
////
//Window graph_px_depth_err() : Graph
//	PauseUpdate; Silent 1		// building window...
//	Display/K=1/W=(277,44,856,575) manyErrs[*][0] vs manyErrs[*][1]
//	ModifyGraph gfMult=110,width={Aspect,1}, mode=3, marker=19, msize=3
//	ModifyGraph zColor(manyErrs)={manyErrs[*][3],-2.5e-10,2.5e-10,RedWhiteBlue256}
//	ModifyGraph tick=2, zero(left)=2, mirror=1, minor=1, lowTrip=0.001, standoff=0
//	Label left "depth (µm)"
//	Label bottom "px"
//	SetAxis/A/N=1/E=2 left
//	SetAxis bottom 0,2048
//	Cursor/P A manyErrs 69
//	ShowInfo
//	SetDrawLayer UserFront
//	SetDrawEnv xcoord= bottom,ycoord= prel,dash= 1
//	DrawLine 1023.5,0,1023.5,1
//EndMacro
////
//Window graph_py_depth_err() : Graph
//	PauseUpdate; Silent 1		// building window...
//	Display/K=1/W=(763,44,1342,575) manyErrs[*][0] vs manyErrs[*][2]
//	ModifyGraph gfMult=110,width={Aspect,1}, mode=3, marker=19, msize=3
//	ModifyGraph zColor(manyErrs)={manyErrs[*][3],-2.5e-10,2.5e-10,RedWhiteBlue256}
//	ModifyGraph tick=2, zero(left)=2, mirror=1, minor=1, lowTrip=0.001, standoff=0
//	Label left "depth (µm)"
//	Label bottom "py"
//	SetAxis/A/N=1/E=2 left
//	SetAxis bottom 0,2048
//	Cursor/P A manyErrs 69
//	ShowInfo
//	SetDrawLayer UserFront
//	SetDrawEnv xcoord= bottom,ycoord= prel,dash= 1
//	DrawLine 1023.5,0,1023.5,1
//EndMacro


//===================================== OLD   OLD   OLD   OLD ===================================

// for a ray coming from depth, and tangent to wire, where will it hit detector, returns detector(y,z) position as complex value
Function/C depthWire2pixelXYZ(g,depth,Xw,xyzMin,xyzMax,edge)
	STRUCT microGeometry &g
	Variable depth										// depth along incident beam, relative to the Si position (µm)
	Wave Xw											// wire (x,y,z) beam line coordinates in wire units (µm)
	Wave xyzMin,xyzMax								// two points on the detector (assumed to be ends) relative to Si
	Variable edge										// 1=leading edge (usual),  0=trailing edge
	Variable testing = (WhichListItem("TestdepthWire2pixelXYZ", GetRTStackInfo(0))>=0)

DoAlert 0,"called depthWire2pixelXYZ(),  this is probably wrong, has not been checked"
	Variable R = (g.wire.dia)/2						// radius of wire (µm)
	Variable kiy=0, kiz=1							// direction of incident beam in beam line coords (do not need kix)
	Variable sz = depth, sy = kiy/kiz * depth		// source position, relative to the Si position (beam line coords)

	Variable vy,vz,v									// (vy,vz) vector from source (at depth) to center of wire
	vy = (Xw[1] - g.wire.origin[1]) - sy
	vz = (Xw[2] - g.wire.origin[2]) - sz
	v = sqrt(vy*vy + vz*vz)							// length of v

	Variable phi0 = atan2(vy,vz)					// angle from positive z to tangent ray measured from z toward y
	Variable dphi =asin(R/v)						// angle between v and ray at tangent point (always positive)
	Variable phi = phi0 + (edge ? -dphi : dphi)
	Variable cotPhi = 1/tan(phi)

	// line from P to S is z = y*cot(phi) + b
	Variable b = sz - sy*cotPhi

	// line along surface of detector, line between xyzMin[] & xyzMax[],   y = m*z + beta1
	Variable m = (xyzMax[1]-xyzMin[1])/(xyzMax[2]-xyzMin[2])
	Variable beta1 = xyzMin[1] - m*xyzMin[2]

	Variable detY, detZ								// y, z parts of the pixel location relative to Si (beam line coords)
	detY = (m*b + beta1) / (1-m*cotPhi)			// y value of intersection
	detZ = detY*cotPhi + b							// from line P to S,   z = y*cot(phi) + b

	if (testing)
		Variable cosTheta = NumVarOrDefault("root:Packages:geometry:cosThetaWire",cos(PI/4))
		Variable sinTheta = NumVarOrDefault("root:Packages:geometry:sinThetaWire",sin(PI/4))
		printf "Hwire=%g, Fwire=%g,   Xw=(%g, %g, %g)\r",Xw[1]*sinTheta + Xw[2]*cosTheta, -Xw[1]*cosTheta + Xw[2]*sinTheta,Xw[0], Xw[1], Xw[2]
		printf "s = (%g, %g),   v = (%g, %g),   R=%gµm,  |wire-sample| = %gµm\r",sy,sz,vy,vz,R,v
		printf "line from source to tangent is:   z = %g*y + %g\r",cotPhi,b
		printf "line between xyzMin & xyzMax is:   y = %g*z + %g\r",m,beta1
		// check that detY,detZ satisfy both equations
		printf " is z-y*cot(phi) = %g == %g=b?, Æ=%g\r",detZ-detY*cotPhi,b,detZ-detY*cotPhi-b
		printf " is y-m*z = %g == %g=beta?,  Æ=%g\r",detY-m*detZ,beta1,detY-m*detZ-beta1
	endif
	return cmplx(detY,detZ)
End
//
//
// returns depth (starting point) of ray ending at pixel (px,py) that is tangent to leading edge of the wire
// the returned depth is relative to the Si position (µm)
Function pixel2depth(g,dm,px,py,Xw,edge)
	STRUCT microGeometry &g
	Variable dm										// detector number {0,1,2}
	Variable px,py									// uncorrected unbinned pixel position (pixels) (zero based)
	Wave Xw											// wire (x,y,z) beam line coordinates in wire units (µm)
	Variable edge										// 1=leading edge (usual),  0=trailing edge
	Variable depth=NaN								// computed depth relative to the Si position (µm)

//	Wave xyzDetector=root:Packages:geometry:pixel2depth_xyzDetector
	pixel2XYZ(g.d[dm],px,py,xyzDetector)			// returns xyz of pixel on detector (µm) relative to the Si position
//	depth = pixelXYZ2depth(g,xyzDetector[1],xyzDetector[2],Xw,edge)
	return depth										// relative to the Si position (µm)
End


Function PanelHwireMake() : Panel
	if (ItemsInList(WinList("PanelWireH",";","WIN:64"))==1)
		DoWindow/F PanelWireH
		return 0
	endif

	STRUCT microGeometry g
	FillGeometryStructDefault(g)
	Variable dnum, py = g.d[dnum].Ny/2
	for (dnum=0;dnum<g.Ndetectors && !(g.d[dnum].used);dnum+=1)
	endfor
	NewPanel /W=(128,166,438,253)/K=1/N=PanelWireH as "H wire calculator"
	PopupMenu detectorNum,pos={10,4},size={71,20},proc=microGeo#PanelH_PopMenuProc
	PopupMenu detectorNum,mode=1,popvalue="Orange",value= #"\"Orange;Yellow;Purple\""
	SetVariable setvar_px,pos={104,3},size={80,19},proc=microGeo#PanelH_SetVarProc,title="px"
	SetVariable setvar_px,fSize=12,limits={0,2047,10},value= _NUM:100
	SetVariable setvar_py,pos={190,3},size={80,19},proc=microGeo#PanelH_SetVarProc,title="py"
	SetVariable setvar_py,fSize=12,limits={0,2047,10},value= _NUM:py
	SetVariable setvar_depth,pos={10,62},size={140,19},proc=microGeo#PanelH_SetVarProc,title="depth (µm) "
	SetVariable setvar_depth,fSize=12,limits={-1000,1000,5},format="%.1f",value= _NUM:0
	PopupMenu popupEdge,pos={10,32},size={107,20},proc=microGeo#PanelH_PopMenuProc,fSize=10
	PopupMenu popupEdge,mode=1,popvalue="Leading Edge",value= #"\"Leading Edge;Trailng Edge\""
	SetVariable setvar_Fwire,pos={162,32},size={130,21},proc=microGeo#PanelH_SetVarProc,title="F\\Bwire\\M (µm) "
	SetVariable setvar_Fwire,fSize=12,limits={-1000,1000,5},value= _NUM:g.wire.F
	SetVariable setvar_Hwire,pos={162,62},size={140,21},proc=microGeo#PanelH_SetVarProc,title="H\\Bwire\\M (µm) "
	SetVariable setvar_Hwire,fSize=12,limits={-1000,1000,5},format="%.2f"
	SetWindow PanelWireH userdata(changed)="2;1;"

	calcHofWire("PanelWireH",0)
	return 0
End
//
Static Function PanelH_PopMenuProc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa
	if (pa.eventCode==2)
		calcHofWire(pa.win,0)			// none of the popups count for whoChanged
	endif
	return 0
End
//
Static Function PanelH_SetVarProc(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable whoChanged=NumberByKey(sva.ctrlName,"setvar_px:1;setvar_depth:2;setvar_Hwire:4")
			whoChanged = numtype(whoChanged) ? 0 : whoChanged
//print "ctrlName =",sva.ctrlName,"    whoChanged = ",whoChanged
			calcHofWire(sva.win,whoChanged)
	endswitch
	return 0
End
//
Static Function calcHofWire(win,whoChanged)
	String win
	Variable whoChanged				// who changed, 0=none, 1=px, 2=depth, 4=Hwire

	Variable i0=NaN, i1=NaN
	if (strlen(win) && whoChanged>0)
		String str=GetUserData(win,"","changed")
		i0=str2num(StringFromList(0,str))
		i1=str2num(StringFromList(1,str))
		if (!(whoChanged==i0))		// did not change the most recent one
			i1 = i0						// move most recent to oldest
			i0 = whoChanged			// set current to most recetn
		endif
		str = num2istr(i0)+";"+num2istr(i1)+";"
		SetWindow $win userdata(changed)=str
	endif
	if (numtype(i0+i1))
		i0 = 1
		i1 = 2
	endif
	Variable dnum						// detector number, {0, 1, 2}
	Variable px,py						// requested pixel on detector
	Variable depth						// requested depth (µm)
	Variable edge						// 1=leading edge (usual),  0=trailing edge
	Variable Fwire						// F of the wire, usually -200
	Variable Hwire						// H of the wire (µm)
	ControlInfo/W=$win setvar_px		;	px = V_Value
	ControlInfo/W=$win setvar_py		;	py = V_Value
	ControlInfo/W=$win setvar_Fwire	;	Fwire = V_Value
	ControlInfo/W=$win setvar_depth	;	depth = V_Value
	ControlInfo/W=$win setvar_Hwire	;	Hwire = V_Value
	ControlInfo/W=$win detectorNum	;	dnum = V_Value-1
	ControlInfo/W=$win popupEdge		;	edge = V_Value==1 ? 1 : 0
	STRUCT microGeometry g
	FillGeometryStructDefault(g)
	edge = edge==0 ? 0 : 1				// default for bad input is leading
	g.wire.F = numtype(Fwire) ? g.wire.F : Fwire
	Make/N=3/D/FREE xyz=NaN

	SetDrawLayer/K UserBack			// re-draw colored box
	Make/N=(3,3)/U/W/FREE rgb
	rgb[0][0]= {65535,65535,65535}
	rgb[0][1]= {43688,65535,30000}
	rgb[0][2]= {32768,0,65535}
	SetDrawEnv linefgc=(rgb[dnum][0],rgb[dnum][1],rgb[dnum][2]),fillfgc=(rgb[dnum][0],rgb[dnum][1],rgb[dnum][2]) // Orange
	DrawRect 1,3,90,23

	SetVariable setvar_px,format=""
	if ((i0+i1)==3)					// calc Hwire from depth & px
		pixel2XYZ(g.d[dnum],px,py,xyz)				// convert pixel position to the beam line coordinate system
		Hwire = DepthPixel2WireH(g,depth,xyz,edge)	// calc H of wire that is tangent to line from depth on beam to pixel on detector (in beam line coords)
		SetVariable setvar_Hwire,value= _NUM:Hwire
	elseif ((i0+i1)==5)				// calc depth from px & Hwire
		Make/N=3/D/FREE Xw=0						// position of wire, assume that the Xwire=0
		Xw[1] = HF2Y(Hwire,g.wire.F)					// Y =  H*sin(angle) - F*cos(angle)
		Xw[2] = HF2Z(Hwire,g.wire.F)					// Z =  H*cos(angle) + F*sin(angle)
		pixel2XYZ(g.d[dnum],px,py,xyz)				// convert pixel position to the beam line coordinate system
		depth = PixelxyzWire2depth(g,xyz,Xw,edge)	// returns depth (µm)
		SetVariable setvar_depth,value= _NUM:depth
	elseif ((i0+i1)==6)				// calc px from depth & Hwire
		Make/N=3/D/FREE Xw=0						// position of wire, assume that the Xwire=0
		Xw[1] = HF2Y(Hwire,g.wire.F)					// Y =  H*sin(angle) - F*cos(angle)
		Xw[2] = HF2Z(Hwire,g.wire.F)					// Z =  H*cos(angle) + F*sin(angle)
		px = WireDepth2Pixel(g,dnum,py,Xw,depth,edge)	// find px, given the wire and depth
		SetVariable setvar_px,format="%.1f",value= _NUM:px
	endif
	return 0
End

// this routine works, but it is not really ready for prime time, and it is not that fast either
Static Function WireDepth2Pixel(g,dnum,py,Xw,depth,edge)	// find px, given the wire, depth, & py
	STRUCT microGeometry &g			// g.wire.F provides the F
	Variable dnum						// detector number, {0, 1, 2}
	Variable py							// y pixel position on detector
	Wave Xw							// center of wire in raw PM500 coordinates (µm)
	Variable depth						// distance along incident beam from origin (the Si positino) (µm)
	Variable edge						// 1=leading edge (usual),  0=trailing edge

	Make/N=3/D/FREE xyz=NaN
	Variable px = g.d[dnum].Nx / 2		// start with px = d.Nx/2, center of image
	Variable dpx=g.d[dnum].Nx/10		// px step, start with 1/10 of full image
	Variable dtest,dlast=Inf				// temporary depth
	Variable dstep=30					// depth step, start with 30µm

	Variable counter=0
	do
		pixel2XYZ(g.d[dnum],px,py,xyz)		// convert pixel position to the beam line coordinate system
		dtest = PixelxyzWire2depth(g,xyz,Xw,edge)// returns depth for this loop
		if (abs(dtest-depth)<0.001)
			break
		elseif (abs(depth-dtest)<abs(depth-dlast))	// getting better, continue doing this
			px += dpx
		else
			dpx /= -2							// reverse step size and divide it by 2
			px += dpx
		endif
		dlast = dtest
		counter += 1
	while(abs(dtest-depth)>0.001 && counter<1000)
	px = (abs(dtest-depth)>0.001) ? NaN : px
	return px
End

//===================================== End of Wire & Depth ===================================
//=======================================================================================



//=======================================================================================
//==================================== Start of Utility Items ===================================

Static Function/S MakeUnique3Vector(preset)		// make a new 3-vector, and optionally set it to preset[]
	Wave preset										// optional wave, if none call with $""
	String wName = UniqueName("vec3_",1,0)
	Make/N=3/D $wName
	Wave vec=$wName
	if (WaveExists(preset))
		vec = preset
	endif
	return GetWavesDataFolder(vec,2)				// return full path for assingment to a wave
End
//
// Used to make temporary internal waves, a typical use would be:    Wave xxx=$MakeUnique3Vector($"")
// You will probably want to KillWaves/Z xxx, at the end of the routine.
Static Function/S MakeUnique3x3Mat(preset)		// make a new 3x3 matrix, and optionally set it to preset[]
	Wave preset										// optional wave, if none call with $""
	String wName = UniqueName("mat3x3_",1,0)
	Make/N=(3,3)/D $wName
	Wave mat=$wName
	if (WaveExists(preset))
		mat = preset
	endif
	return GetWavesDataFolder(mat,2)				// return full path for assingment to a wave
End

Static Function/S keyStrFromBuffer(buf)			// read in all of the tagged values from a file into a keyword list
	String buf

	buf = ReplaceString("\r",buf,"\n")
	buf = ReplaceString("\t",buf," ")
	buf += "\n"

	Variable i, i0=0,i1, imax=strlen(buf)
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
			for (i=0;char2num(line[i+1])>32;i+=1)	// find end of tag, it ends with a space or lower
			endfor
			tagName = line[1,i]
			if (strlen(tagName)<1)
				continue
			endif

			// check if tag is already in list, only take the first instance of a tag, not the last
			if (keyInList(tagName,list,"=",""))		// check if this key is already in list
				continue
			endif

			for (i=i+1;char2num(line[i])<=32;i+=1)	// find first non-white space
			endfor
			value = line[i,Inf]						// value associated with tagName

			for (i=strlen(value)-1;char2num(value[i])<=32 && i>0;i-=1)	// strip off trailing whitespace
			endfor
			value = ReplaceString(";",value[0,i],":")
			list = ReplaceStringByKey(tagName,list,value,"=")
		endif
		i0 = i1 + 1
	while(i0<imax)
	return list
End

// ===================================== End of Utility Items ===================================
//=======================================================================================



//=======================================================================================
// ==================================== Start of Math Items ====================================

Function angleVec2Vec(a,b)						// return the angle between two vectors (degree)
	Wave a,b
	Variable dot = MatrixDot(a,b) / (norm(a)*norm(b))
	dot = limit(dot,-1,1)							// ensure that the acos will exist
	return acos(dot)*180/PI
End


Function rotationAngleOfMat(rot)				// returns the total rotation angle of a matrix 'rot'
	Wave rot										// the rotation matrix
	Variable trace = MatrixTrace(rot)			// trace = 1 + 2*cos(theta)
	Variable cosine = (trace-1)/2				// cosine of the rotation angle
	cosine = (cosine>1) ? (2-cosine) : cosine
	return acos(cosine)*180/PI					// rotation angle in degrees
End


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
	Wave R1,R2									// two Rodriques vectors

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
ThreadSafe Function rotationMatAboutAxis(axis,angle,mat)
	Wave axis				// axis about which to rotate (or possibly Rodriques vector)
	Variable angle			// angle to rotate (degrees), assumes axis is true Rodriques vector if angle invalid
	Wave mat				// desired rotation matrix

	Variable len = norm(axis)
	angle = numtype(angle) ? 2*atan(len) : angle*PI/180	// the rotation angle (rad)

	if (angle==0)									// zero angle rotation is just the identity matrix
		mat = (p==q)
		return 0
	endif

	Variable nx=axis[0]/len, ny=axis[1]/len, nz=axis[2]/len
	Variable cosa=cos(angle), sina=sin(angle)
	Variable c1 = 1-cosa
	if (DimSize(mat,0)!=3 && DimSize(mat,1)!=3)
		return NaN
	endif
	// from		http://mathworld.wolfram.com/RodriguesRotationFormula.html (I double checked this too.)
	mat[0][0] = cosa+nx*nx*c1;			mat[0][1] =nx*ny*c1-nz*sina;			mat[0][2] = nx*nz*c1+ny*sina;
	mat[1][0] = nx*ny*c1+nz*sina;		mat[1][1] = cosa+ny*ny*c1;			mat[1][2] =ny*nz*c1-nx*sina;
	mat[2][0] = nx*nz*c1-ny*sina;		mat[2][1] = ny*nz*c1+nx*sina;		mat[2][2] = cosa+nz*nz*c1;
	return 0
End


ThreadSafe Function/WAVE matString2mat33(str)			// returns a FREE wave with (3x3) matrix
	String str
	Variable a0,a1,a2, b0,b1,b2, c0,c1,c2
	sscanf str, "{{%g,%g,%g}{%g,%g,%g}{%g,%g,%g}}", a0,a1,a2, b0,b1,b2, c0,c1,c2
	if (V_flag==9)
		Make/N=(3,3)/D/FREE mat
		mat[0][0]= {a0,a1,a2}					// the reciprocal lattice
		mat[0][1]= {b0,b1,b2}
		mat[0][2]= {c0,c1,c2}
		return mat
	else
		return $""
	endif
End


ThreadSafe Function/WAVE getRLfrom3DWave(scatter,point)
	Wave scatter
	Variable point					// point number, only used when scatter is a list of triplets (not for 3d wave)

	if (!WaveExists(scatter))
		Wave RL=$""
	elseif (WaveDims(scatter)==2 && DimSize(scatter,1)==3 && numtype(point)==0)	// a list of triplets
		String path=GetWavesDataFolder(scatter,1)
		Wave RX = $(path+"RX"), RY = $(path+"RY"), RZ = $(path+"RZ")				// Use Rx,Ry,Rz,   NOT gm
		if (WaveExists(RX) && WaveExists(RY) && WaveExists(RZ))
			Make/N=3/D/FREE 	Rvec = {RX[point], RY[point], RZ[point]}
			Make/N=(3,3)/D/FREE rot
			rotationMatAboutAxis(Rvec,NaN,rot)			// rotation of reference recip to this point
			Wave RL = matString2mat33(StringByKey("recipRef",note(RX),"="))	// starts as reference recip
			MatrixOP/FREE/O RL = rot x RL					// RL is now the reciprocal lattice for this point
		endif
	elseif (WaveDims(scatter)==3)						// a 3D volume (probably k-space)
		String recip_lattice0=StringByKey("recip_lattice0",note(scatter),"=")
		Wave RL = matString2mat33(recip_lattice0)
	else
		Wave RL=$""
	endif
	return RL
End


Function notRotationMatrix(mat)				// false if mat is a rotation matrix, i.e. mat^T = mat^-1
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
		returnVal = sumd							// not a rotation matrix
	elseif (MatrixDet(mat)<0)
		returnVal = -1							// an improper rotation matrix
	else
		returnVal = 0								// yes it is a rotation matrix
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
		Cross vec0,vec1							// Z = X x Y
		Wave W_Cross=W_Cross
		vec2 = W_Cross
		Cross vec2,vec0							// Y = Z x X
		vec1 = W_Cross
	elseif (norm1>=norm0 && norm1>=norm2)// Y is longest
		Cross vec1,vec2							// X = Y x Z
		Wave W_Cross=W_Cross
		vec0 = W_Cross
		Cross vec0,vec1							// Z = X x Y
		vec2 = W_Cross
	else											// Z is longest
		Cross vec2,vec0							// Y = Z x X
		Wave W_Cross=W_Cross
		vec1 = W_Cross
		Cross vec1,vec2							// X = Y x Z
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


Function EulerMatrix(alpha,bet,gam)			// Make the rotation matrix M_Euler from Euler angles (degree)
	Variable alpha,bet,gam						// Euler angles, passed as degrees
	alpha *= PI/180								// need radians internally
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


Function rotationAboutX(theta,rot)				// make a rotation matrix about X axis
	Variable theta									// rotation angle (radian)
	Wave rot										// a 3x3 matrix
	Variable c=cos(theta), s=sin(theta)
	rot[0][0] = 1;		rot[0][1] = 0;		rot[0][2] = 0
	rot[1][0] = 0;		rot[1][1] = c;		rot[1][2] = -s
	rot[2][0] = 0;		rot[2][1] = s;		rot[2][2] = c
End
//
Function rotationAboutY(theta,rot)				// make a rotation matrix about Y axis
	Variable theta									// rotation angle (radian)
	Wave rot										// a 3x3 matrix
	Variable c=cos(theta), s=sin(theta)
	rot[0][0] = c;			rot[0][1] = 0;		rot[0][2] = s
	rot[1][0] = 0;		rot[1][1] = 1;		rot[1][2] = 0
	rot[2][0] = -s;		rot[2][1] = 0;		rot[2][2] = c
End
//
Function rotationAboutZ(theta,rot)				// make a rotation matrix about Z axis
	Variable theta									// rotation angle (radian)
	Wave rot										// a 3x3 matrix
	Variable c=cos(theta), s=sin(theta)
	rot[0][0] = c;			rot[0][1] = -s;		rot[0][2] = 0
	rot[1][0] = s;			rot[1][1] = c;			rot[1][2] = 0
	rot[2][0] = 0;		rot[2][1] = 0;		rot[2][2] = 1
End

// ==================================== End of Math Items =====================================
//=======================================================================================



//=======================================================================================
// ================================ Start of Geometry Sub-Panel =================================

Function MakeGeometryParametersPanel(strStruct)
	String strStruct													// optional passed value of geo structure, this is used if passed
	if (strlen(WinList("GeometrySet",";","WIN:64")))			// window alreay exits, bring it to front
		DoWindow/F GeometrySet
		return 0
	endif
	NewPanel /K=1 /W=(675,60,895,653)
	DoWindow/C GeometrySet
	SetWindow kwTopWin,hook(closeWin)=microPanelHook	// used to save values if window is killed
	FillGeometryParametersPanel(strStruct,"GeometrySet",0,0)
End
//
Function/T FillGeometryParametersPanel(strStruct,hostWin,left,top)
	String strStruct													// optional passed value of geo structure, this is used if passed
	String hostWin													// name of home window
	Variable left, top													// offsets from the left and top

	Init_microGeo()
	NewDataFolder/O root:Packages:geometry:PanelValues			// ensure that the needed data folders exist

	STRUCT microGeometry g
	if (strlen(strStruct)>0)
		StructGet/S/B=2 g, strStruct								// found structure information, load into geo
		Variable/G root:Packages:geometry:PanelValues:dirty=1	// for passed structure
	else
		FillGeometryStructDefault(g)								//fill the geometry structure with current values
		Variable/G root:Packages:geometry:PanelValues:dirty=0
	endif

	// create globals for the panel
	Variable/G root:Packages:geometry:PanelValues:SampleOrigX = g.s.O[0], root:Packages:geometry:PanelValues:SampleOrigY = g.s.O[1], root:Packages:geometry:PanelValues:SampleOrigZ = g.s.O[2]
	Variable/G root:Packages:geometry:PanelValues:SampleRotX = g.s.R[0], root:Packages:geometry:PanelValues:SampleRotY = g.s.R[1], root:Packages:geometry:PanelValues:SampleRotZ = g.s.R[2]

	Variable/G root:Packages:geometry:PanelValues:Ndetectors = g.Ndetectors
	NVAR Ndetectors = root:Packages:geometry:PanelValues:Ndetectors	

	Variable/G root:Packages:geometry:PanelValues:iDetector = 0		// detector displayed {0,1,2}
	NVAR iDetector = root:Packages:geometry:PanelValues:iDetector	

	Variable/G root:Packages:geometry:PanelValues:used0 = g.d[0].used
	Variable/G root:Packages:geometry:PanelValues:Nx0 = g.d[0].Nx,				root:Packages:geometry:PanelValues:Ny0 = g.d[0].Ny
	Variable/G root:Packages:geometry:PanelValues:sizeX0 = g.d[0].sizeX/1e3,	root:Packages:geometry:PanelValues:sizeY0 = g.d[0].sizeY/1e3	// panel uses mm, g.d[m] uses µm
	Variable/G root:Packages:geometry:PanelValues:Rx0 = g.d[0].R[0], root:Packages:geometry:PanelValues:Ry0 = g.d[0].R[1], root:Packages:geometry:PanelValues:Rz0 = g.d[0].R[2]
	Variable/G root:Packages:geometry:PanelValues:Px0 = g.d[0].P[0]/1e3
	Variable/G root:Packages:geometry:PanelValues:Py0 = g.d[0].P[1]/1e3
	Variable/G root:Packages:geometry:PanelValues:Pz0 = g.d[0].P[2]/1e3

	Variable/G root:Packages:geometry:PanelValues:used1 = g.d[1].used
	Variable/G root:Packages:geometry:PanelValues:Nx1 = g.d[1].Nx,				root:Packages:geometry:PanelValues:Ny1 = g.d[1].Ny
	Variable/G root:Packages:geometry:PanelValues:sizeX1 = g.d[1].sizeX/1e3,	root:Packages:geometry:PanelValues:sizeY1 = g.d[1].sizeY/1e3	// panel uses mm, g.d[m] uses µm
	Variable/G root:Packages:geometry:PanelValues:Rx1 = g.d[1].R[0], 	root:Packages:geometry:PanelValues:Ry1 = g.d[1].R[1], root:Packages:geometry:PanelValues:Rz1 = g.d[1].R[2]
	Variable/G root:Packages:geometry:PanelValues:Px1 = g.d[1].P[0]/1e3
	Variable/G root:Packages:geometry:PanelValues:Py1 = g.d[1].P[1]/1e3
	Variable/G root:Packages:geometry:PanelValues:Pz1 = g.d[1].P[2]/1e3

	Variable/G root:Packages:geometry:PanelValues:used2 = g.d[2].used
	Variable/G root:Packages:geometry:PanelValues:Nx2 = g.d[2].Nx,				root:Packages:geometry:PanelValues:Ny2 = g.d[2].Ny
	Variable/G root:Packages:geometry:PanelValues:sizeX2 = g.d[2].sizeX/1e3,	root:Packages:geometry:PanelValues:sizeY2 = g.d[2].sizeY/1e3	// panel uses mm, g.d[m] uses µm
	Variable/G root:Packages:geometry:PanelValues:Rx2 = g.d[2].R[0], root:Packages:geometry:PanelValues:Ry2 = g.d[2].R[1], root:Packages:geometry:PanelValues:Rz2 = g.d[2].R[2]
	Variable/G root:Packages:geometry:PanelValues:Px2 = g.d[2].P[0]/1e3
	Variable/G root:Packages:geometry:PanelValues:Py2 = g.d[2].P[1]/1e3
	Variable/G root:Packages:geometry:PanelValues:Pz2 = g.d[2].P[2]/1e3

	Variable/G root:Packages:geometry:PanelValues:dia = g.wire.dia
	Variable/G root:Packages:geometry:PanelValues:knife = g.wire.knife
	Variable/G root:Packages:geometry:PanelValues:wireF = g.wire.F
	Variable/G root:Packages:geometry:PanelValues:wireX0 = g.wire.origin[0], root:Packages:geometry:PanelValues:wireY0 = g.wire.origin[1], root:Packages:geometry:PanelValues:wireZ0 = g.wire.origin[2]
	Variable/G root:Packages:geometry:PanelValues:axisX = g.wire.axis[0], root:Packages:geometry:PanelValues:axisY = g.wire.axis[1], root:Packages:geometry:PanelValues:axisZ = g.wire.axis[2]
	Variable/G root:Packages:geometry:PanelValues:wireRotX = g.wire.R[0], root:Packages:geometry:PanelValues:wireRotY = g.wire.R[1], root:Packages:geometry:PanelValues:wireRotZ = g.wire.R[2]

	String/G root:Packages:geometry:PanelValues:timeMeasured0 = g.d[0].timeMeasured
	String/G root:Packages:geometry:PanelValues:geoNote0 = g.d[0].geoNote
	String/G root:Packages:geometry:PanelValues:detectorID0 = g.d[0].detectorID
	String/G root:Packages:geometry:PanelValues:distortionMapFile0 = g.d[0].distortionMapFile
	String/G root:Packages:geometry:PanelValues:timeMeasured1 = g.d[1].timeMeasured
	String/G root:Packages:geometry:PanelValues:geoNote1 = g.d[1].geoNote
	String/G root:Packages:geometry:PanelValues:detectorID1 = g.d[1].detectorID
	String/G root:Packages:geometry:PanelValues:distortionMapFile1 = g.d[1].distortionMapFile
	String/G root:Packages:geometry:PanelValues:timeMeasured2 = g.d[2].timeMeasured
	String/G root:Packages:geometry:PanelValues:geoNote2 = g.d[2].geoNote
	String/G root:Packages:geometry:PanelValues:detectorID2 = g.d[2].detectorID
	String/G root:Packages:geometry:PanelValues:distortionMapFile2 = g.d[2].distortionMapFile

	SetWindow kwTopWin,userdata(GeoPanelName)=hostWin+"#geoPanel"
	NewPanel/K=1/W=(left,top,left+221,top+603)/HOST=$hostWin
	ModifyPanel frameStyle=0, frameInset=0
	RenameWindow #,geoPanel
	SetDrawLayer UserBack
	SetDrawEnv fname= "Lucida Grande",fsize= 16
	DrawText 20,20,"Geometry"

	ValDisplay Ndetectors,pos={120,8},size={70,13},title="detectors",value= #"root:Packages:geometry:PanelValues:Ndetectors"
	ValDisplay Ndetectors,labelBack=(56797,56797,56797),format="%g",frame=0, limits={0,0,0},barmisc={0,1000}
	ValDisplay Ndetectors,help={"Number of detectors in use [1-3]"}

	// Detectors
	Variable v=55
	SetDrawEnv linethick= 2
	DrawLine 0,v-28,220,v-28
	PopupMenu detectorPopup,pos={10,v-22},size={203,20},proc=GeoDetectorPopMenuProc,title="Detector "
	PopupMenu detectorPopup,help={"Select detector to define or display"}
	PopupMenu detectorPopup,font="Lucida Grande",fSize=14,mode=1,value=#popStrForDetector(g)

	CheckBox used,pos={10,v},size={100,18},font="Lucida Grande",fSize=14
	CheckBox used,proc=GeoPanelUsedCheckBoxProc,variable= root:Packages:geometry:PanelValues:used0
	SetVariable Nx,pos={19,21+v},size={94,15},proc=GeoPanelVarChangedProc,title="Nx"
	SetVariable Nx,help={"x size of detectors (pixels)"},format="%d"
	SetVariable Nx,limits={1,5000,0},value= root:Packages:geometry:PanelValues:Nx0,bodyWidth= 60
	SetVariable Ny,pos={112,21+v},size={94,15},proc=GeoPanelVarChangedProc,title="Ny"
	SetVariable Ny,help={"y size of detector (pixels)"},format="%d"
	SetVariable Ny,limits={1,5000,0},value= root:Packages:geometry:PanelValues:Ny0,bodyWidth= 60
	SetVariable sizeX,pos={0,39+v},size={112,15},proc=GeoPanelVarChangedProc,title="X (mm)"
	SetVariable sizeX,help={"outside x size of detector (mm)"}
	SetVariable sizeX,limits={0,5100,0},value= root:Packages:geometry:PanelValues:sizeX0,bodyWidth= 60
	SetVariable sizeY,pos={94,39+v},size={112,15},proc=GeoPanelVarChangedProc,title="Y (mm)"
	SetVariable sizeY,help={"outside y size of detector (mm)"}
	SetVariable sizeY,limits={0,500,0},value= root:Packages:geometry:PanelValues:sizeY0,bodyWidth= 60
	SetVariable Rx,pos={8,58+v},size={75,15},proc=GeoPanelVarChangedProc,title="R"
	SetVariable Rx,help={"rotation axis for detector 0"},format="%.6f"
	SetVariable Rx,limits={-3.2,3.2,0},value= root:Packages:geometry:PanelValues:Rx0,bodyWidth= 50
	SetVariable Ry,pos={85,58+v},size={55,15},proc=GeoPanelVarChangedProc,title=" "
	SetVariable Ry,help={"direction of wire axis, beam line coords"},format="%.6f"
	SetVariable Ry,limits={-3.2,3.2,0},value= root:Packages:geometry:PanelValues:Ry0,bodyWidth= 50
	SetVariable Rz,pos={149,58+v},size={50,15},proc=GeoPanelVarChangedProc,title=" "
	SetVariable Rz,help={"direction of wire axis, beam line coords"},format="%.6f"
	SetVariable Rz,limits={-3.2,3.2,0},value= root:Packages:geometry:PanelValues:Rz0,bodyWidth= 50
	SetVariable Px,pos={8,76+v},size={75,15},proc=GeoPanelVarChangedProc,title="P"
	SetVariable Px,help={"rotation axis for detector 0"},format="%.3f"
	SetVariable Px,limits={-1000,1000,0},value= root:Packages:geometry:PanelValues:Px0,bodyWidth= 50
	SetVariable Py,pos={85,76+v},size={55,15},proc=GeoPanelVarChangedProc,title=" "
	SetVariable Py,help={"direction of wire axis, beam line coords"},format="%.3f"
	SetVariable Py,limits={-1000,1000,0},value= root:Packages:geometry:PanelValues:Py0,bodyWidth= 50
	SetVariable Pz,pos={149,76+v},size={50,15},proc=GeoPanelVarChangedProc,title=" "
	SetVariable Pz,help={"direction of wire axis, beam line coords"},format="%.3f"
	SetVariable Pz,limits={-1000,1000,0},value= root:Packages:geometry:PanelValues:Pz0,bodyWidth= 50
	SetVariable timeMeasured,pos={10,96+v},size={202,15},title="time"
	SetVariable timeMeasured,value= root:Packages:geometry:PanelValues:timeMeasured0
	SetVariable timeMeasured,help={"date/time when detector calibrated"}
	SetVariable geoNote,pos={10,116+v},size={202,15},title="note"
	SetVariable geoNote,value= root:Packages:geometry:PanelValues:geoNote0
	SetVariable geoNote,help={"optional note about this detector"}
	SetVariable distortionMapFile,pos={10,136+v},size={202,15},title="distortion file"
	SetVariable distortionMapFile,value= root:Packages:geometry:PanelValues:distortionMapFile0
	SetVariable distortionMapFile,help={"name of file with the distortion map for this detector"}
	SetVariable detectorID,pos={10,156+v},size={202,15},title="detector ID",proc=SetDetectorColorProc
	SetVariable detectorID,value= root:Packages:geometry:PanelValues:detectorID0, labelBack=0
	SetVariable detectorID,help={"ID for this detector"}

	// Wire
	v = 253
	SetDrawEnv linethick= 2
	DrawLine 0,v-20,220,v-20
	SetDrawEnv fname= "Lucida Grande",fsize= 16
	DrawText 10,v,"Scanning wire"
	SetVariable dia,pos={10,v+2},size={81,15},proc=GeoPanelVarChangedProc,title="dia (µm)",bodyWidth= 40
	SetVariable dia,limits={1,5100,0},value= root:Packages:geometry:PanelValues:dia
	SetVariable dia,help={"diameter of wire (micron)"}
	PopupMenu knifePopup,pos={110,v},size={107,20}, proc=GeoKnifePopMenuProc
	PopupMenu knifePopup,mode=(1+g.wire.knife),value= #"\"free-standing;knife edge;\""
	PopupMenu knifePopup,help={"select if scanning wire is a free standing wire or a true knife edge"}
	SetVariable wireAxisX,pos={9,22+v},size={68,15},proc=GeoPanelVarChangedProc,title="axis"
	SetVariable wireAxisX,help={"direction of wire axis, beam line coords"},format="%.6f"
	SetVariable wireAxisX,limits={-1000,1000,0},value= root:Packages:geometry:PanelValues:axisX,bodyWidth= 45
	SetVariable wireAxisY,pos={84,22+v},size={45,15},proc=GeoPanelVarChangedProc,title=" "
	SetVariable wireAxisY,help={"direction of wire axis, beam line coords"},format="%.6f"
	SetVariable wireAxisY,limits={-1000,1000,0},value= root:Packages:geometry:PanelValues:axisY,bodyWidth= 45
	SetVariable wireAxisZ,pos={138,22+v},size={45,15},proc=GeoPanelVarChangedProc,title=" "
	SetVariable wireAxisZ,help={"direction of wire axis, beam line coords"},format="%.6f"
	SetVariable wireAxisZ,limits={-1000,1000,0},value= root:Packages:geometry:PanelValues:axisZ,bodyWidth= 45
	SetDrawEnv fsize= 9, textxjust= 1
	DrawText 102,51+v,"wire positioner origin (µm)"
	SetVariable wXO,pos={17,51+v},size={60,15},proc=GeoPanelVarChangedProc,title=" ",bodyWidth= 60
	SetVariable wXO,help={"X of wire when centered on origin"}
	SetVariable wXO,limits={-50000,50000,0},value= root:Packages:geometry:PanelValues:wireX0
	SetVariable wYO,pos={76,51+v},size={60,15},proc=GeoPanelVarChangedProc,title=" ",bodyWidth= 60
	SetVariable wYO,help={"Y of wire when centered on origin"}
	SetVariable wYO,limits={-50000,50000,0},value= root:Packages:geometry:PanelValues:wireY0
	SetVariable wZO,pos={135,51+v},size={60,15},proc=GeoPanelVarChangedProc,title=" ",bodyWidth= 60
	SetVariable wZO,help={"Z of wire when centered on origin"}
	SetVariable wZO,limits={-50000,50000,0},value= root:Packages:geometry:PanelValues:wireZ0
	SetDrawEnv fsize= 9, textxjust= 1
	DrawText 102,79+v,"wire positioner rotation R (rad)"
	SetVariable wRotX,pos={17,79+v},size={60,15},proc=GeoPanelVarChangedProc,title=" ",bodyWidth= 60
	SetVariable wRotX,help={"rotation vector for wire positioner, X-component (rad)"}
	SetVariable wRotX,limits={-50000,50000,0},value= root:Packages:geometry:PanelValues:wireRotX
	SetVariable wRotY,pos={76,79+v},size={60,15},proc=GeoPanelVarChangedProc,title=" ",bodyWidth= 60
	SetVariable wRotY,help={"rotation vector for wire positioner, Y-component (rad)"}
	SetVariable wRotY,limits={-50000,50000,0},value= root:Packages:geometry:PanelValues:wireRotY
	SetVariable wRotZ,pos={135,79+v},size={60,15},proc=GeoPanelVarChangedProc,title=" ",bodyWidth= 60
	SetVariable wRotZ,help={"rotation vector for wire positioner, Z-component (rad)"}
	SetVariable wRotZ,limits={-50000,50000,0},value= root:Packages:geometry:PanelValues:wireRotZ

	// Sample
	v = 374
	SetDrawEnv linethick= 2
	DrawLine 0,v-20,220,v-20
	SetDrawEnv fname= "Lucida Grande",fsize= 16
	DrawText 10,v,"Sample"
	SetDrawEnv fsize= 9, textxjust= 1
	DrawText 102,12+v,"sample positioner origin (µm)"
	SetVariable sXO,pos={17,12+v},size={60,15},proc=GeoPanelVarChangedProc,title=" ",bodyWidth= 60
	SetVariable sXO,help={"X of sample when centered on origin"}
	SetVariable sXO,limits={-50000,50000,0},value= root:Packages:geometry:PanelValues:SampleOrigX
	SetVariable sYO,pos={76,12+v},size={60,15},proc=GeoPanelVarChangedProc,title=" ",bodyWidth= 60
	SetVariable sYO,help={"Y of sample when centered on origin"}
	SetVariable sYO,limits={-50000,50000,0},value= root:Packages:geometry:PanelValues:SampleOrigY
	SetVariable sZO,pos={135,12+v},size={60,15},proc=GeoPanelVarChangedProc,title=" ",bodyWidth= 60
	SetVariable sZO,help={"Z of sample when centered on origin"}
	SetVariable sZO,limits={-50000,50000,0},value= root:Packages:geometry:PanelValues:SampleOrigZ

	SetDrawEnv fsize= 9, textxjust= 1
	DrawText 102,40+v,"sample positioner rotation R (rad)"
	SetVariable sRotX,pos={17,40+v},size={60,15},proc=GeoPanelVarChangedProc,title=" ",bodyWidth= 60
	SetVariable sRotX,help={"rotation vector for sample positioner, X-component (rad)"}
	SetVariable sRotX,limits={-50000,50000,0},value= root:Packages:geometry:PanelValues:SampleRotX
	SetVariable sRotY,pos={76,40+v},size={60,15},proc=GeoPanelVarChangedProc,title=" ",bodyWidth= 60
	SetVariable sRotY,help={"rotation vector for sample positioner, Y-component (rad)"}
	SetVariable sRotY,limits={-50000,50000,0},value= root:Packages:geometry:PanelValues:SampleRotY
	SetVariable sRotZ,pos={135,40+v},size={60,15},proc=GeoPanelVarChangedProc,title=" ",bodyWidth= 60
	SetVariable sRotZ,help={"rotation vector for sample positioner, Z-component (rad)"}
	SetVariable sRotZ,limits={-50000,50000,0},value= root:Packages:geometry:PanelValues:SampleRotZ

	// Buttons
	v = 435
	SetDrawEnv linethick= 2
	DrawLine 0,v,220,v
	Button buttonSaveDefault,pos={35,440},size={150,20},proc=microGeo#GeometryPanelButtonProc,title="Save as Default"
	Button buttonSaveDefault,help={"saves these values as default in the current data folder"}
	Button buttonSaveLocal,pos={35,462},size={150,20},proc=microGeo#GeometryPanelButtonProc,title="Save in Current Fldr"
	Button buttonSaveLocal,help={"saves these values as default in the current data folder"}
	Button buttonEPICS,pos={145,484},size={40,20},proc=microGeo#GeometryPanelButtonProc,title="EPICS"
	Button buttonEPICS,help={"Fill the values in this panel from EPICS"}
	Button buttonFillFromWeb,pos={35,484},size={105,20},proc=microGeo#GeometryPanelButtonProc,title="Load from web"
	Button buttonFillFromWeb,help={"Fill the values in this panel from the web"}
	Button buttonFillFromFile,pos={35,506},size={150,20},proc=microGeo#GeometryPanelButtonProc,title="Load from a file"
	Button buttonFillFromFile,help={"Fill the values in this panel from a file"}
	if (exists("EPICS_put_PV_num")==6 && exists("EPICS_put_PV_str")==6)
		PopupMenu SaveGeoPopup,pos={35,530},size={150,20},bodyWidth=150,proc=microGeo#GeoSaveFilePopMenuProc
		PopupMenu SaveGeoPopup,mode=1,popvalue="Save Geometry to...",value= #"\"Write Geo to a File;Put Geo to EPICS;\""
		PopupMenu SaveGeoPopup,help={"Write values in this panel to a file or to EPICS"}
	else
		Button buttonWriteToFile,pos={35,528},size={150,20},proc=microGeo#GeometryPanelButtonProc,title="Write to a file"
		Button buttonWriteToFile,help={"Write values in this panel to a file"}
	endif
	Button buttonPrintGeo,pos={35,550},size={150,20},proc=microGeo#GeometryPanelButtonProc,title="Print Geo to History"
	Button buttonPrintGeo,help={"print geometry values to the history"}
	Button buttonCancel,pos={35,572},size={150,20},proc=microGeo#GeometryPanelButtonProc,title="Cancel"
	Button buttonCancel,help={"Close this panel, do not save any changes"}
	SetWindow kwTopWin,userdata(KillHookFn)="SetGeoPanelHook"
	ConnectControlsToGlobals(hostWin+"#geoPanel",0)
	GeoPanelDetectorDisable("#")
	GeoPanelDirtyUpdate(hostWin+"#geoPanel",NaN)
	return "#geoPanel"
End
//
// this is called by the pop up menu to select a particular detector
Function GeoDetectorPopMenuProc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa
	if (pa.eventCode != 2)
		return 0
	endif

	NVAR iDetector = root:Packages:geometry:PanelValues:iDetector		// detector displayed {0,1,2}
	iDetector = pa.popNum - 1
	ConnectControlsToGlobals(pa.win,iDetector)
	GeoPanelDetectorDisable(pa.win)
	return 0
End
//
Static Function/T popStrForDetector(g)
	STRUCT microGeometry &g

	String id0 = SelectString(g.d[0].used,"Detector 0",g.d[0].detectorID)
	String id1 = SelectString(g.d[1].used,"Detector 1",g.d[1].detectorID)
	String id2 = SelectString(g.d[2].used,"Detector 2",g.d[2].detectorID)
	id0 = SelectString(strlen(id0),"Detector 0",id0)
	id1 = SelectString(strlen(id1),"Detector 1",id1)
	id2 = SelectString(strlen(id2),"Detector 2",id2)

	String popStr
	sprintf popStr "\"%s;%s;%s;\"",id0,id1,id2
	return popStr
End
//
Static Function ConnectControlsToGlobals(win,iDetector)
	String win
	Variable iDetector				//  0, 1, or 2 only
	if (iDetector!=0 && iDetector!=1 && iDetector!=2)
		DoAlert 0,"IDetector must be 0, 1, or 2.   iDetector = "+num2str(iDetector)
		return 1
	endif

	String path = "root:Packages:geometry:PanelValues:", is = num2istr(iDetector)
	Variable used = NumVarOrDefault(path+"used"+is,0)
	String id = StrVarOrDefault("root:Packages:geometry:PanelValues:detectorID"+num2istr(iDetector),"")
	id = SelectString(strlen(id) && used, "Detector "+num2istr(iDetector),id)
	CheckBox used,win=$win,title=id, variable= $(path+"used"+is)
	SetVariable Nx,win=$win, value= $(path+"Nx"+is)
	SetVariable Ny,win=$win,value= $(path+"Ny"+is)
	SetVariable sizeX,win=$win,value= $(path+"sizeX"+is)
	SetVariable sizeY,win=$win,value= $(path+"sizeY"+is)
	SetVariable Rx,win=$win,value= $(path+"Rx"+is)
	SetVariable Ry,win=$win,value= $(path+"Ry"+is)
	SetVariable Rz,win=$win,value= $(path+"Rz"+is)
	SetVariable Px,win=$win,value= $(path+"Px"+is)
	SetVariable Py,win=$win,value= $(path+"Py"+is)
	SetVariable Pz,win=$win,value= $(path+"Pz"+is)
	SetVariable timeMeasured,win=$win,value= $(path+"timeMeasured"+is)
	SetVariable geoNote,win=$win,value= $(path+"geoNote"+is)
	SetVariable distortionMapFile,win=$win,value= $(path+"distortionMapFile"+is)
	SetVariable detectorID,win=$win,value= $(path+"detectorID"+is)
	return 0
End
//
// this is called by all of the check boxes, and it just sets the dirty flag for the panel
Function GeoPanelUsedCheckBoxProc(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba
	if (cba.eventCode != 2)
		return 0
	endif
	GeoPanelDetectorDisable(cba.win)
	GeoPanelDirtyUpdate(cba.win,1)					// set dirty = 1
	return 0
End
//
// this is called by the wire knife edge pop up menu to set the dirty flag when it is changed
Function GeoKnifePopMenuProc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa
	if (pa.eventCode != 2)
		return 0
	endif
	NVAR knife = root:Packages:geometry:PanelValues:knife		// 0=free-standing,  1=knife-edge
	knife = pa.popNum - 1
	GeoPanelDirtyUpdate(pa.win,1)					// set dirty = 1
	return 0
End
//
// this is called by all of the varables, and it just sets the dirty flag for the panel
Function GeoPanelVarChangedProc(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	if ( sva.eventCode == 3 || sva.eventCode<1)			// return unless number changed
		return 0
	endif
	GeoPanelDirtyUpdate(sva.win,1)					// set dirty = 1
	return 0
End
//
//Function SetDetectorColorProc(sva) : SetVariableControl
//	STRUCT WMSetVariableAction &sva
//	if (sva.eventCode!=2)			// only process "Enter key"
//		return 0
//	endif
//
//	SetDrawLayer/W=microPanel#geoPanel/K ProgBack
//	strswitch(ReplaceString(",",sva.sval,""))
//		case "PE1621 723-3335":				// Orange
//			SetVariable $(sva.ctrlName) win=$(sva.win), labelBack=(65535,43688,32768)
//			SetDrawEnv/W=microPanel#geoPanel linethick=8, linefgc= (65535,43688,32768)
//			DrawLine/W=microPanel#geoPanel 1,30,1,230
//			break
//		case "PE0820 763-1807":				// Yellow
//			SetVariable $(sva.ctrlName) win=$(sva.win), labelBack=(65535,65535,0)
//			SetDrawEnv/W=microPanel#geoPanel linethick=8, linefgc= (65535,65535,0)
//			DrawLine/W=microPanel#geoPanel 1,30,1,230
//			break
//		case "PE0820 763-1850":				// Purple
//			SetVariable $(sva.ctrlName) win=$(sva.win), labelBack=(65535,30000,65535)
//			SetDrawEnv/W=microPanel#geoPanel linethick=8, linefgc= (65535,30000,65535)
//			DrawLine/W=microPanel#geoPanel 1,30,1,230
//			break
//		default:
//			SetVariable $(sva.ctrlName) win=$(sva.win), labelBack=0
//	endswitch
//	return 0
//End
//
Function SetDetectorColorProc(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	if (sva.eventCode!=2)			// only process "Enter key"
		return 0
	endif

	String color =  detectorID2color(ReplaceString(",",sva.sval,""))
	Variable r=40000,g=40000,b=40000				// default color is gray
	if (stringmatch(color,"Orange"))
		r = 65535;		g = 43688;		b = 32768	
	elseif (stringmatch(color,"Yellow"))
		r = 65535;		g = 65535;		b = 0
	elseif (stringmatch(color,"Purple"))
		r = 65535;		g = 30000;		b = 65535
	else
		color = ""
	endif

	SetDrawLayer/W=microPanel#geoPanel/K ProgBack
	SetVariable $(sva.ctrlName) win=$(sva.win), labelBack=(r,g,b)
	if (strlen(color))
		SetDrawEnv/W=microPanel#geoPanel linethick=8, linefgc= (r,g,b)
		DrawLine/W=microPanel#geoPanel 1,30,1,230
	endif
	return 0
End
//
Function/T detectorID2color(detectorID)
	String detectorID
	Variable i
	for (i=0;i<MAX_Ndetectors;i+=1)
		strswitch(detectorID)
			case "PE1621 723-3335":				// Orange
				return "Orange"
				break
			case "PE0820 763-1807":				// Yellow
				return "Yellow"
				break
			case "PE0820 763-1850":				// Purple
				return "Purple"
		endswitch
	endfor
	return ""
End
//
Static Function GeoPanelDetectorDisable(win)			// enable/disable detector fields based on check box
	String win

	Variable iDetector =NumVarOrDefault("root:Packages:geometry:PanelValues:iDetector",0)
	Variable disable = NumVarOrDefault("root:Packages:geometry:PanelValues:used"+num2istr(iDetector),0) ? 0 : 2

	SetVariable Nx,win=$win, disable=disable;		SetVariable Ny,win=$win, disable=disable
	SetVariable sizeX,win=$win, disable=disable;	SetVariable sizeY,win=$win, disable=disable
	SetVariable Rx,win=$win, disable=disable;		SetVariable Ry,win=$win, disable=disable;		SetVariable Rz,win=$win, disable=disable
	SetVariable Px,win=$win, disable=disable;		SetVariable Py,win=$win, disable=disable;		SetVariable Pz,win=$win, disable=disable
	SetVariable timeMeasured,win=$win, disable=disable
	SetVariable geoNote,win=$win, disable=disable
	SetVariable distortionMapFile,win=$win, disable=disable
	SetVariable detectorID,win=$win, disable=disable

	String id = StrVarOrDefault("root:Packages:geometry:PanelValues:detectorID"+num2istr(iDetector),"")
	id = SelectString(strlen(id) && !disable, "Detector "+num2istr(iDetector),id)
	CheckBox used,win=$win,title=id

	STRUCT microGeometry g
	SetGeometryStructToPanelGlobals(g)				// set structure to the values in the geo panel globals
	PopupMenu detectorPopup,win=$win,value=#popStrForDetector(g)

	// set detector color
	ControlInfo/W=$win detectorID
	STRUCT WMSetVariableAction sva
	sva.sval = StrVarOrDefault(S_DataFolder+S_Value,"")
	sva.ctrlName = "detectorID"
	sva.win = win
	sva.eventCode = 2
	SetDetectorColorProc(sva)

	return 0
End
//
Function SetGeoPanelHook(s)
	STRUCT WMWinHookStruct &s
	if (s.eventCode !=2 && s.eventCode !=14)			// only process kill events (both windows & sub-windows)
		return 0
	endif
	Variable dirty = NumVarOrDefault("root:Packages:geometry:PanelValues:dirty",1)
	if (!dirty)
		return 0											// nothing was changed, do not ask any questions
	endif

	// values were changed in the panel and not saved, ask the user what to do.
	Variable action=1
	Prompt action, "choose an action, do not press cancel",Popup,"Save in default data folder;Save in current data folder;Close window & Discard these values"
	do
		DoPrompt/Help="3D-Xray Diffraction[geometry]" "changes were made, but not saved",action
	while(V_flag)
	// action = V_flag ? 3 : action						// if cancel was pushed, do not save, treat it like a "just quit"
	STRUCT WMButtonAction B_Struct
	B_Struct.ctrlName = StringFromList(action-1,"buttonSaveDefault;buttonSaveLocal;buttonCancel")
	B_Struct.win = s.winName
	GeometryPanelButtonProc(B_Struct)
//	GeometryPanelButtonProc(StringFromList(action-1,"buttonSaveDefault;buttonSaveLocal;buttonCancel"))
End
//
Static Function GeoPanelDirtyUpdate(win,dirty)		// change dirty global, and update buttons to reflect dirty status
	String win
	Variable dirty

	if (!(dirty>=0))
		dirty = NumVarOrDefault("root:Packages:geometry:PanelValues:dirty",0)
	endif
	Variable/G root:Packages:geometry:PanelValues:dirty = dirty	// make sure the global is set
	if (dirty)
		Button buttonSaveDefault disable=0,fColor=(10000,50000,20000),title="Save Default  (Needed)",win=$win
		Button buttonSaveLocal disable=0,fColor=(65535,65535,40000),win=$win
	else
		Button buttonSaveDefault disable=2,fColor=(0,0,0),title="Save Default",win=$win
		Button buttonSaveLocal disable=2,fColor=(0,0,0),win=$win
	endif
	Button buttonEPICS disable=(exists("get_mult_PV")==6 ? 0 : 2),win=$win
	NVAR Ndetectors = root:Packages:geometry:PanelValues:Ndetectors	
	NVAR used0=root:Packages:geometry:PanelValues:used0
	NVAR used1=root:Packages:geometry:PanelValues:used1
	NVAR used2=root:Packages:geometry:PanelValues:used2
	Ndetectors = !(!used0) + !(!used1) + !(!used2)
End
//
// this is called by the save geometry to File/EPICS pop up menu to save the geometry.
Static Function GeoSaveFilePopMenuProc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa
	if (pa.eventCode != 2)
		return 0
	endif
	PopupMenu $pa.ctrlName,win=$pa.win,mode=1,popvalue="Save Geometry to...",value= #"\"Write Geo to a File;Put Geo to EPICS;\""
	STRUCT microGeometry g
	SetGeometryStructToPanelGlobals(g)					// set structure to the values in the geo panel globals
	if (stringmatch(pa.popStr,"Put Geo to EPICS"))
		putGeo2EPICS(g)									// put values of g into EPICS
	elseif (stringmatch(pa.popStr,"Write Geo to a File"))
		WriteGeoToFile("","",g,"",NaN)
	else
		return 1											// This is an ERROR
	endif
	return 0
End
//
Static Function GeometryPanelButtonProc(B_Struct) : ButtonControl
	STRUCT WMButtonAction &B_Struct
	if (B_Struct.eventCode != 2)							// only accept mouse-up
		return 0
	endif
	String ctrlName= B_Struct.ctrlName

	STRUCT microGeometry g								// used to connect to the globals in the geoPanel

	if (stringmatch(ctrlName,"buttonCancel"))			// close window, do not save numbers
		GeoPanelDirtyUpdate(B_Struct.win,0)			// set dirty = 0
		String win = StringFromLIst(0,WinList("*",";","WIN:64"))
		SetWindow $win,hook(closeWin)= $""
		DoWindow/K $win
		return 0
	elseif (stringmatch(ctrlName,"buttonSaveDefault"))
		String/G root:Packages:geometry:geoStructStr=""	// use default location
		SVAR geoStructStr = root:Packages:geometry:geoStructStr
	elseif (stringmatch(ctrlName,"buttonSaveLocal"))
		String/G geoStructStr=""							// save it in local folder
		SVAR geoStructStr = geoStructStr
	elseif (stringmatch(ctrlName,"buttonEPICS"))
		if (GeoFromEPICS(g))
			return 1
		endif
		SetGeometryPanelGlobals(g)						// set the globals in the geo panel from values in g
		GeoPanelDirtyUpdate(B_Struct.win,1)			// set dirty = 1
		GeoPanelDetectorDisable(B_Struct.win)
		return 0
	elseif (stringmatch(ctrlName,"buttonFillFromWeb"))
		Variable useLocalImage=0, epoch=NaN
		String  list=WaveListClass("rawImage*","*","DIMS:2")
		if (ItemsInList(list)>0)						// some loaded images to choose
			Prompt useLocalImage,"image for date/time",popup,"Loaded Image;File with Image"
			DoPrompt "Image Source",useLocalImage
			if (V_flag)
				return 0
			endif
			useLocalImage = useLocalImage==1
		endif
		if (useLocalImage)								// get epoch from a loaded image
			String imageName=StringFromLIst(0,list)
			if (ItemsInList(list)>1)
				Prompt imageName,"Image for Date",popup,list
				DoPrompt "Image",imageName
				if (V_flag)
					return 0
				endif
			endif
			Wave image = $imageName
			if (!WaveExists(image))
				return 1
			endif
			FUNCREF fileTime2EpochProto func=$("fileTime2Epoch"+SelectString(exists("fileTime2Epoch")==6,"Proto",""))
			epoch = func(StringByKey("file_time",note(image),"="))
//			#if (exists("fileTime2Epoch")==6)
//				epoch = fileTime2Epoch(StringByKey("file_time",note(image),"="))
//			#endif
		else
			Variable fnum
			String GrandImageFileFilter="HDF Files (*.h5, *.hdf5, *.hdf):.h5,.hdf5,.hdf,;SPE Files (*.SPE):.SPE;TIFF Files (*.tif,*.tiff):.tif,.tiff;All Files:.*;"
			Open/D/M="Select raw, NOT reconstructed image"/P=home/R/F=GrandImageFileFilter fnum
			if (strlen(S_fileName)>0)
				GetFileFolderInfo/Q/Z S_fileName
				if (V_Flag)
					return 1
				endif
				epoch = V_modificationDate
			endif
		endif
		if (numtype(epoch))		// one last chance to find the time
			DoAlert 1,"Use Current time"
			if (V_flag==1)
				epoch = DateTime
			endif
		endif
		if (numtype(epoch))
			return 1
		elseif (GeoFromWeb(epoch,g))
			String errStr=""
			if (epoch>0)
				errStr="Invalid Date/Time from image, '"+Secs2Date(epoch,2)+"  "+Secs2Time(epoch,1)+"', must be in [1995-Now]."
			else
				errStr = "Failed to get valid epoch from loaded image, epoch = "+num2str(epoch)+"."
			endif
			DoAlert 0, errStr+"\rFailed to get geometry from web."
			return 1
		endif
		SetGeometryPanelGlobals(g)						// set the globals in the geo panel from values in g
		GeoPanelDirtyUpdate(B_Struct.win,1)			// set dirty = 1
		GeoPanelDetectorDisable(B_Struct.win)
		return 0
	elseif (stringmatch(ctrlName,"buttonFillFromFile"))
		// if (ReadGeoFromKeyFile("","home",g))
		if (ReadGeoFromfile("","home",g))				// fill g with values from file
			return 1
		endif
		SetGeometryPanelGlobals(g)						// set the globals in the geo panel from values in g
		GeoPanelDirtyUpdate(B_Struct.win,1)			// set dirty = 1
		GeoPanelDetectorDisable(B_Struct.win)
		return 0
	elseif (stringmatch(ctrlName,"buttonWriteToFile"))
		SetGeometryStructToPanelGlobals(g)			// set structure to the values in the geo panel globals
		WriteGeoToFile("","",g,"",NaN)
		return 0
	elseif (stringmatch(ctrlName,"buttonPrintGeo"))	// print shown geometry parameters to the history
		SetGeometryStructToPanelGlobals(g)			// set structure to the values in the geo panel globals
		GeometryUpdateCalc(g)							// calculate other values
		printGeometry(g)									// print them all out
		return 0
	else
		DoAlert 0, "GeometryPanelButtonProc() called from unsupported control = '"+ctrlName+"'"
		return 1
	endif

	// store the global values from panel into the 'geoStructStr''
	SetGeometryStructToPanelGlobals(g)				// set structure to the values in the geo panel globals
	GeometryUpdateCalc(g)
	StructPut/S/B=2 g, geoStructStr
	GeoPanelDirtyUpdate(B_Struct.win,0)				// set dirty = 0
	SavePackagePreferences/FLSH=1 "microGeo","microGeoNPrefs",0,g
End
//
Function fileTime2EpochProto(fileTime,[UTC])
	String fileTime
	Variable UTC
	return NaN
End
//
Static Function SetGeometryPanelGlobals(g)			// set the globals in the geo panel from values in g
	STRUCT microGeometry &g							// connect to the globals from this structure
	NVAR Ndetectors = root:Packages:geometry:PanelValues:Ndetectors	
	NVAR used0=root:Packages:geometry:PanelValues:used0
	NVAR Nx0=root:Packages:geometry:PanelValues:Nx0,			Ny0=root:Packages:geometry:PanelValues:Ny0
	NVAR sizeX0=root:Packages:geometry:PanelValues:sizeX0,		sizeY0=root:Packages:geometry:PanelValues:sizeY0
	NVAR Rx0=root:Packages:geometry:PanelValues:Rx0,			Ry0=root:Packages:geometry:PanelValues:Ry0,		Rz0=root:Packages:geometry:PanelValues:Rz0
	NVAR Px0=root:Packages:geometry:PanelValues:Px0,			Py0=root:Packages:geometry:PanelValues:Py0,		Pz0=root:Packages:geometry:PanelValues:Pz0
	NVAR used1=root:Packages:geometry:PanelValues:used1
	NVAR Nx1=root:Packages:geometry:PanelValues:Nx1,			Ny1=root:Packages:geometry:PanelValues:Ny1
	NVAR sizeX1=root:Packages:geometry:PanelValues:sizeX1,		sizeY1=root:Packages:geometry:PanelValues:sizeY1
	NVAR Rx1=root:Packages:geometry:PanelValues:Rx1,			Ry1=root:Packages:geometry:PanelValues:Ry1,		Rz1=root:Packages:geometry:PanelValues:Rz1
	NVAR Px1=root:Packages:geometry:PanelValues:Px1,			Py1=root:Packages:geometry:PanelValues:Py1, 		Pz1=root:Packages:geometry:PanelValues:Pz1
	NVAR used2=root:Packages:geometry:PanelValues:used2
	NVAR Nx2=root:Packages:geometry:PanelValues:Nx2,			Ny2=root:Packages:geometry:PanelValues:Ny2
	NVAR sizeX2=root:Packages:geometry:PanelValues:sizeX2,		sizeY2=root:Packages:geometry:PanelValues:sizeY2
	NVAR Rx2=root:Packages:geometry:PanelValues:Rx2,			Ry2=root:Packages:geometry:PanelValues:Ry2,		Rz2=root:Packages:geometry:PanelValues:Rz2
	NVAR Px2=root:Packages:geometry:PanelValues:Px2,			Py2=root:Packages:geometry:PanelValues:Py2,		Pz2=root:Packages:geometry:PanelValues:Pz2
	NVAR dia	= root:Packages:geometry:PanelValues:dia
	NVAR knife = root:Packages:geometry:PanelValues:knife
	NVAR wireF= root:Packages:geometry:PanelValues:wireF
	NVAR wireX0=root:Packages:geometry:PanelValues:wireX0
	NVAR wireY0=root:Packages:geometry:PanelValues:wireY0
	NVAR wireZ0=root:Packages:geometry:PanelValues:wireZ0
	NVAR axisX=root:Packages:geometry:PanelValues:axisX,		axisY=root:Packages:geometry:PanelValues:axisY,	axisZ=root:Packages:geometry:PanelValues:axisZ
	NVAR wireRotX = root:Packages:geometry:PanelValues:wireRotX
	NVAR wireRotY = root:Packages:geometry:PanelValues:wireRotY
	NVAR wireRotZ = root:Packages:geometry:PanelValues:wireRotZ

	Ndetectors = !(!used0) + !(!used1) + !(!used2)			// set the panel globals from the struct
	used0	= g.d[0].used
	Nx0	= g.d[0].Nx;			Ny0	= g.d[0].Ny
	sizeX0	= g.d[0].sizeX/1e3;	sizeY0	= g.d[0].sizeY/1e3	// panel uses mm, g.d[m] uese µm
	Rx0	= g.d[0].R[0];			Ry0	= g.d[0].R[1];			Rz0	= g.d[0].R[2]
	Px0	= g.d[0].P[0]/1e3;	Py0	= g.d[0].P[1]/1e3;	Pz0	= g.d[0].P[2]/1e3
	used1	= g.d[1].used
	Nx1	= g.d[1].Nx;			Ny1	= g.d[1].Ny
	sizeX1	= g.d[1].sizeX/1e3;	sizeY1	= g.d[1].sizeY/1e3
	Rx1	= g.d[1].R[0];			Ry1	= g.d[1].R[1];			Rz1	= g.d[1].R[2]
	Px1	= g.d[1].P[0]/1e3;	Py1	= g.d[1].P[1]/1e3;	Pz1	= g.d[1].P[2]/1e3
	used2	= g.d[2].used
	Nx2	= g.d[2].Nx;			Ny2	= g.d[2].Ny
	sizeX2	= g.d[2].sizeX/1e3;	sizeY2	= g.d[2].sizeY/1e3
	Rx2	= g.d[2].R[0];			Ry2	= g.d[2].R[1];			Rz2	= g.d[2].R[2]
	Px2	= g.d[2].P[0]/1e3;	Py2	= g.d[2].P[1]/1e3;	Pz2	= g.d[2].P[2]/1e3
	dia = g.wire.dia
	wireX0=g.wire.origin[0];	wireY0=g.wire.origin[1];	wireZ0=g.wire.origin[2]
	axisX = g.wire.axis[0];		axisY = g.wire.axis[1];		axisZ = g.wire.axis[2]
	wireRotX = g.wire.R[0];		wireRotY = g.wire.R[1];		wireRotZ = g.wire.R[2]
	knife = g.wire.knife
	wireF = g.wire.F

	Ndetectors = !(!used0) + !(!used1) + !(!used2)

	SVAR timeMeasured0 = root:Packages:geometry:PanelValues:timeMeasured0
	SVAR geoNote0 = root:Packages:geometry:PanelValues:geoNote0
	SVAR detectorID0 = root:Packages:geometry:PanelValues:detectorID0
	SVAR distortionMapFile0 = root:Packages:geometry:PanelValues:distortionMapFile0

	SVAR timeMeasured1 = root:Packages:geometry:PanelValues:timeMeasured1
	SVAR geoNote1 = root:Packages:geometry:PanelValues:geoNote1
	SVAR detectorID1 = root:Packages:geometry:PanelValues:detectorID1
	SVAR distortionMapFile1 = root:Packages:geometry:PanelValues:distortionMapFile1

	SVAR timeMeasured2 = root:Packages:geometry:PanelValues:timeMeasured2
	SVAR geoNote2 = root:Packages:geometry:PanelValues:geoNote2
	SVAR detectorID2 = root:Packages:geometry:PanelValues:detectorID2
	SVAR distortionMapFile2 = root:Packages:geometry:PanelValues:distortionMapFile2
	timeMeasured0		= g.d[0].timeMeasured
	geoNote0				= g.d[0].geoNote
	detectorID0			= g.d[0].detectorID
	distortionMapFile0	= g.d[0].distortionMapFile
	timeMeasured1		= g.d[1].timeMeasured
	geoNote1				= g.d[1].geoNote
	detectorID1			= g.d[1].detectorID
	distortionMapFile1	= g.d[1].distortionMapFile
	timeMeasured2		= g.d[2].timeMeasured
	geoNote2				= g.d[2].geoNote
	detectorID2			= g.d[2].detectorID
	distortionMapFile2	= g.d[2].distortionMapFile

	NVAR SampleRotX = root:Packages:geometry:PanelValues:SampleRotX
	NVAR SampleRotY = root:Packages:geometry:PanelValues:SampleRotY
	NVAR SampleRotZ = root:Packages:geometry:PanelValues:SampleRotZ
	SampleRotX = g.s.R[0];	SampleRotY = g.s.R[1];		SampleRotZ = g.s.R[2]

	NVAR SampleOrigX=root:Packages:geometry:PanelValues:SampleOrigX
	NVAR SampleOrigY=root:Packages:geometry:PanelValues:SampleOrigY
	NVAR SampleOrigZ=root:Packages:geometry:PanelValues:SampleOrigZ
	SampleOrigX = g.s.O[0];	SampleOrigY = g.s.O[1];		SampleOrigZ = g.s.O[2]

	return 0
End
//
Static Function SetGeometryStructToPanelGlobals(g)		// set structure to the values in the geo panel globals
	STRUCT microGeometry &g								// connect to the globals

	NVAR Ndetectors = root:Packages:geometry:PanelValues:Ndetectors	
	NVAR used0=root:Packages:geometry:PanelValues:used0
	NVAR Nx0=root:Packages:geometry:PanelValues:Nx0,			Ny0=root:Packages:geometry:PanelValues:Ny0
	NVAR sizeX0=root:Packages:geometry:PanelValues:sizeX0,	sizeY0=root:Packages:geometry:PanelValues:sizeY0
	NVAR Rx0=root:Packages:geometry:PanelValues:Rx0,			Ry0=root:Packages:geometry:PanelValues:Ry0,		Rz0=root:Packages:geometry:PanelValues:Rz0
	NVAR Px0=root:Packages:geometry:PanelValues:Px0,			Py0=root:Packages:geometry:PanelValues:Py0,		Pz0=root:Packages:geometry:PanelValues:Pz0
	NVAR used1=root:Packages:geometry:PanelValues:used1
	NVAR Nx1=root:Packages:geometry:PanelValues:Nx1,			Ny1=root:Packages:geometry:PanelValues:Ny1
	NVAR sizeX1=root:Packages:geometry:PanelValues:sizeX1,	sizeY1=root:Packages:geometry:PanelValues:sizeY1
	NVAR Rx1=root:Packages:geometry:PanelValues:Rx1,			Ry1=root:Packages:geometry:PanelValues:Ry1,		Rz1=root:Packages:geometry:PanelValues:Rz1
	NVAR Px1=root:Packages:geometry:PanelValues:Px1,			Py1=root:Packages:geometry:PanelValues:Py1, 		Pz1=root:Packages:geometry:PanelValues:Pz1
	NVAR used2=root:Packages:geometry:PanelValues:used2
	NVAR Nx2=root:Packages:geometry:PanelValues:Nx2,			Ny2=root:Packages:geometry:PanelValues:Ny2
	NVAR sizeX2=root:Packages:geometry:PanelValues:sizeX2,	sizeY2=root:Packages:geometry:PanelValues:sizeY2
	NVAR Rx2=root:Packages:geometry:PanelValues:Rx2,			Ry2=root:Packages:geometry:PanelValues:Ry2,		Rz2=root:Packages:geometry:PanelValues:Rz2
	NVAR Px2=root:Packages:geometry:PanelValues:Px2,			Py2=root:Packages:geometry:PanelValues:Py2,		Pz2=root:Packages:geometry:PanelValues:Pz2
	NVAR dia	= root:Packages:geometry:PanelValues:dia
	NVAR knife = root:Packages:geometry:PanelValues:knife
	NVAR wireF = root:Packages:geometry:PanelValues:wireF
	NVAR wireX0=root:Packages:geometry:PanelValues:wireX0
	NVAR wireY0=root:Packages:geometry:PanelValues:wireY0
	NVAR wireZ0=root:Packages:geometry:PanelValues:wireZ0
	NVAR axisX=root:Packages:geometry:PanelValues:axisX,		axisY=root:Packages:geometry:PanelValues:axisY,	axisZ=root:Packages:geometry:PanelValues:axisZ
	Ndetectors = !(!used0) + !(!used1) + !(!used2)

	g.Ndetectors = Ndetectors										// set all the geo values from global values
	g.d[0].used = used0
	g.d[0].Nx = Nx0;			g.d[0].Ny = Ny0
	g.d[0].sizeX = sizeX0*1e3;	g.d[0].sizeY = sizeY0*1e3	// panel uses mm, g.d[m] uese µm
	g.d[0].R[0] = Rx0;			g.d[0].R[1] = Ry0;			g.d[0].R[2] = Rz0
	g.d[0].P[0] = Px0*1e3;	g.d[0].P[1] = Py0*1e3;	g.d[0].P[2] = Pz0*1e3
	g.d[1].used = used1
	g.d[1].Nx = Nx1;			g.d[1].Ny = Ny1
	g.d[1].sizeX = sizeX1*1e3;	g.d[1].sizeY = sizeY1*1e3
	g.d[1].R[0] = Rx1;			g.d[1].R[1] = Ry1;			g.d[1].R[2] = Rz1
	g.d[1].P[0] = Px1*1e3;	g.d[1].P[1] = Py1*1e3;	g.d[1].P[2] = Pz1*1e3
	g.d[2].used = used2
	g.d[2].Nx = Nx2;			g.d[2].Ny = Ny2
	g.d[2].sizeX = sizeX2*1e3;	g.d[2].sizeY = sizeY2*1e3
	g.d[2].R[0] = Rx2;			g.d[2].R[1] = Ry2;			g.d[2].R[2] = Rz2
	g.d[2].P[0] = Px2*1e3;	g.d[2].P[1] = Py2*1e3;	g.d[2].P[2] = Pz2*1e3
	g.wire.dia = dia
	g.wire.origin[0] = wireX0;	g.wire.origin[1] = wireY0;	g.wire.origin[2] = wireZ0
	g.wire.axis[0] = axisX;		g.wire.axis[1] = axisY;		g.wire.axis[2] = axisZ
	g.wire.knife = knife
	g.wire.F = wireF


	SVAR timeMeasured0 = root:Packages:geometry:PanelValues:timeMeasured0
	SVAR geoNote0 = root:Packages:geometry:PanelValues:geoNote0
	SVAR detectorID0 = root:Packages:geometry:PanelValues:detectorID0
	SVAR distortionMapFile0 = root:Packages:geometry:PanelValues:distortionMapFile0
	SVAR timeMeasured1 = root:Packages:geometry:PanelValues:timeMeasured1
	SVAR geoNote1 = root:Packages:geometry:PanelValues:geoNote1
	SVAR detectorID1 = root:Packages:geometry:PanelValues:detectorID1
	SVAR distortionMapFile1 = root:Packages:geometry:PanelValues:distortionMapFile1
	SVAR timeMeasured2 = root:Packages:geometry:PanelValues:timeMeasured2
	SVAR geoNote2 = root:Packages:geometry:PanelValues:geoNote2
	SVAR detectorID2 = root:Packages:geometry:PanelValues:detectorID2
	SVAR distortionMapFile2 = root:Packages:geometry:PanelValues:distortionMapFile2
	g.d[0].timeMeasured	= timeMeasured0
	g.d[0].geoNote			= geoNote0
	g.d[0].detectorID		= detectorID0
	g.d[0].distortionMapFile= distortionMapFile0
	g.d[1].timeMeasured	= timeMeasured1
	g.d[1].geoNote			= geoNote1
	g.d[1].detectorID		= detectorID1
	g.d[1].distortionMapFile= distortionMapFile1
	g.d[2].timeMeasured	= timeMeasured2
	g.d[2].geoNote			= geoNote2
	g.d[2].detectorID		= detectorID2
	g.d[2].distortionMapFile= distortionMapFile2

	NVAR wireRotX = root:Packages:geometry:PanelValues:wireRotX
	NVAR wireRotY = root:Packages:geometry:PanelValues:wireRotY
	NVAR wireRotZ = root:Packages:geometry:PanelValues:wireRotZ
	g.wire.R[0] = wireRotX;	g.wire.R[1] = wireRotY;	g.wire.R[2] = wireRotZ;

	NVAR SampleRotX = root:Packages:geometry:PanelValues:SampleRotX
	NVAR SampleRotY = root:Packages:geometry:PanelValues:SampleRotY
	NVAR SampleRotZ = root:Packages:geometry:PanelValues:SampleRotZ
	g.s.R[0] = SampleRotX;	g.s.R[1] = SampleRotY;	g.s.R[2] = SampleRotZ;

	NVAR SampleOrigX = root:Packages:geometry:PanelValues:SampleOrigX
	NVAR SampleOrigY = root:Packages:geometry:PanelValues:SampleOrigY
	NVAR SampleOrigZ = root:Packages:geometry:PanelValues:SampleOrigZ
	g.s.O[0] = SampleOrigX;	g.s.O[1] = SampleOrigY;	g.s.O[2] = SampleOrigZ

	return 0
End

//=======================================================================================

// With an Igor epoch, retrieve the geo data form the web server
Function GeoFromWeb(epoch,gIn)
	Variable epoch						// Igor epoch, also value of V_modificationDate from GetFileFolderInfo
	STRUCT microGeometry &gIn
	if (!(epoch>date2secs(1995,1,1)) || !(epoch<=DateTime))	// invalid epoch
		return 1
	endif

	String ddate, ttime				// ddate="2007-03-16",ttime="14:32:55"
	Variable im,id,iy
	sscanf Secs2Date(epoch,-1),"%d/%d/%d",id,im,iy
	sprintf ddate "%04d-%02d-%02d",iy,im,id
	ttime = Secs2Time(epoch,3)

	String url
	sprintf url "http://%s/index.php?tag=geoN&date=%sT%s",GeoWebServer,ddate,ttime
	//	printf "url = %s\r\r",url
	String result = FetchURL(url)
	if (strsearch(result,"ERROR",0)==0)
		return 1
	endif

	Variable N=strlen(result), dquote=char2num("\"")			// ASCII value of a double-quote
	if (char2num(result[0])==dquote)							// remove leading double-quote from Apple Script
		result = result[1,N-1]
		N -= 1
	endif
	if (char2num(result[N-1])==dquote)						// remove trailing double-quote from Apple Script
		result = result[0,N-2]
		N -= 1
	endif
	String list = keyStrFromBuffer(result)
	String type = StringByKey("filetype",list,"=")
	String xml = xmlContents(result)
	STRUCT microGeometry g

	if (strsearch(type,"geometryFileN",0,2)==0)				// this is a "$key value" file
		if (GeoFromKeyValueList(list,g))
			g.Ndetectors = -1									// ensure that geometry gets flagged as bad
		endif
	elseif (strlen(xml))										// try as xml
		if (GeoFromXML(xml,g))
			g.Ndetectors = -1									// ensure that geometry gets flagged as bad
		endif
	else															// the input is unknown
		DoAlert 0, "ERROR, GeoFromWeb(), got filetype='"+type+"',  it must be 'geometryFileN'"
		return 1
	endif

	if (MicroGeometryBad(g))									// check for a valid or Invalid structure
		print "Not loading geometry from web, values are invalid"
		DoAlert 0, "Not loading geometry from web, values are invalid"
	else
		GeometryUpdateCalc(g)									// calculate other values
		CopymicroGeometry(gIn,g)								// new geometry is valid, so copy it in
		if (ItemsInList(GetRTStackInfo(0))<3)					// print everything if run from command line
			printf "Loaded geometry information from web using date: %s,  %s\r",Secs2Date(epoch,2),Secs2Time(epoch,1)
//			if (strlen(list))
//				printf "     this geometry created on   %s  %s\r",StringByKey("dateWritten",list,"="),StringByKey("timeWritten",list,"=")
//				printf "     geometry file was set on  %s, %s\r",StringByKey("dateWritten",list,"="),StringByKey("timeWritten",list,"=")
//				String	str = StringByKey("fileNote",list,"=")
//				if (strlen(str))
//					printf "     file note='%s'\r",str
//				endif
//			else
//				printf "     detector0 was set on  %s\r",g.d[0].timeMeasured
//				str = g.d[0].geoNote
//				if (strlen(str))
//					printf "     with detector note='%s'\r",str
//				endif
//			endif
		endif
	endif
	return 0
End

Function GeoFromEPICS(gIn)	//fill the geometry structure from EPICS (uses caget)
	STRUCT microGeometry &gIn

	String geoList= "Ndetectors.RVAL;"
	geoList += "used0;Nx0;Ny0;sizeX0;sizeY0;Rx0;Ry0;Rz0;Px0;Py0;Pz0;timeMeasured0;geoNote0;detectorID0;distortionMapFile0;"
	geoList += "used1;Nx1;Ny1;sizeX1;sizeY1;Rx1;Ry1;Rz1;Px1;Py1;Pz1;timeMeasured1;geoNote1;detectorID1;distortionMapFile1;"
	geoList += "used2;Nx2;Ny2;sizeX2;sizeY2;Rx2;Ry2;Rz2;Px2;Py2;Pz2;timeMeasured2;geoNote2;detectorID2;distortionMapFile2;"
	geoList += "wireDia;wireOriginX;wireOriginY;wireOriginZ;wireKnife;wireAxisX;wireAxisY;wireAxisZ;"
	geoList += "wireRotX;wireRotY;wireRotZ;wireF;"
	geoList += "SampleRotX;SampleRotY;SampleRotZ;SampleOriginX;SampleOriginY;SampleOriginZ;"

	String item,pvList=""
	Variable i
	// make list of pv names for each item in geoList (items in geo that have EPICS pv's)
	for (i=0,item=StringFromList(0,geoList); strlen(item); i+=1,item=StringFromList(i,geoList))
		pvList += EPICS_PREFIX+item+";"
	endfor

	// use FUNCREF incase this computer does not have the function get_mult_PV()
	String funcSpec = SelectString(exists("get_mult_PV")==6 , "protoEPICS_geo", "get_mult_PV")
	FUNCREF protoEPICS_geo  func=$funcSpec
	String epicsOut = func(pvList)								// get values from EPICS
	if (strlen(epicsOut)<2)
		return 1
	endif
	// turn result of get_mult_PV() into a keyword=value list
	Variable n=strlen(EPICS_PREFIX)
	pvList = ""
	for (i=0,item=StringFromList(0,epicsOut); strlen(item); i+=1,item=StringFromList(i,epicsOut))
		if (strsearch(item,EPICS_PREFIX,0)==0)				// remove the prefix
			item = item[n,Inf]
		endif
		item = ReplaceString(".RVAL=",item,"=")
//		Variable len = strlen(item)
//		if (strsearch(item[len-5,len-1],".RVAL",0)==0)		// remove trailing ".RVAL"
//			item = item[0,len-6]
//		endif
		pvList += item+";"
	endfor

	Variable xx,yy,zz
	String si
	STRUCT microGeometry g
	g.s.R[0] = NumberByKey("SampleRotX",pvList,"=")
	g.s.R[1] = NumberByKey("SampleRotY",pvList,"=")
	g.s.R[2] = NumberByKey("SampleRotZ",pvList,"=")
	if (numtype(g.s.R[0] + g.s.R[1] + g.s.R[2]))
		g.s.R[0] = 0;	g.s.R[1] = 0;	g.s.R[2] = 0
	endif
	g.s.O[0] = NumberByKey("SampleOriginX",pvList,"=")
	g.s.O[1] = NumberByKey("SampleOriginY",pvList,"=")
	g.s.O[2] = NumberByKey("SampleOriginZ",pvList,"=")
//	if (numtype(g.s.O[0] + g.s.O[1] + g.s.O[2]))
//		g.s.O[0] = 0;	g.s.O[1] = 0;	g.s.O[2] = 0
//	endif

	g.Ndetectors = round(NumberByKey("Ndetectors",pvList,"="))
	for (i=0;i<MAX_Ndetectors;i+=1)
		si = num2istr(i)
		g.d[i].used = !numtype(NumberByKey("sizeX"+si,pvList,"="))
		if (!(g.d[i].used))											// skip un-used detectors
			continue
		endif
		g.d[i].Nx = NumberByKey("Nx"+si,pvList,"=")
		g.d[i].Ny = NumberByKey("Ny"+si,pvList,"=")
		g.d[i].sizeX = NumberByKey("sizeX"+si,pvList,"=")*1e3	// EPICS uses mm, structure uses µm
		g.d[i].sizeY = NumberByKey("sizeY"+si,pvList,"=")*1e3
		g.d[i].R[0] = NumberByKey("Rx"+si,pvList,"=")
		g.d[i].R[1] = NumberByKey("Ry"+si,pvList,"=")
		g.d[i].R[2] = NumberByKey("Rz"+si,pvList,"=")
		g.d[i].P[0] = NumberByKey("Px"+si,pvList,"=")*1e3	// EPICS uses mm, structure uses µm
		g.d[i].P[1] = NumberByKey("Py"+si,pvList,"=")*1e3
		g.d[i].P[2] = NumberByKey("Pz"+si,pvList,"=")*1e3
		g.d[i].timeMeasured = StringByKey("timeMeasured"+si,pvList,"=")
		g.d[i].geoNote = StringByKey("geoNote"+si,pvList,"=")
		g.d[i].detectorID = StringByKey("detectorID"+si,pvList,"=")
		g.d[i].distortionMapFile = StringByKey("distortionMapFile"+si,pvList,"=")
	endfor

	g.wire.dia			= NumberByKey("wireDia",pvList,"=")
	g.wire.origin[0]	= NumberByKey("wireOriginX",pvList,"=")
	g.wire.origin[1]	= NumberByKey("wireOriginY",pvList,"=")
	g.wire.origin[2]	= NumberByKey("wireOriginZ",pvList,"=")
	g.wire.knife		= NumberByKey("wireKnife",pvList,"=")
	g.wire.F			= NumberByKey("wireF",pvList,"=")
	g.wire.axis[0]		= NumberByKey("wireAxisX",pvList,"=")
	g.wire.axis[1]		= NumberByKey("wireAxisY",pvList,"=")
	g.wire.axis[2]		= NumberByKey("wireAxisZ",pvList,"=")
	g.wire.R[0]		= NumberByKey("wireRotX",pvList,"=")
	g.wire.R[1]		= NumberByKey("wireRotY",pvList,"=")
	g.wire.R[2]		= NumberByKey("wireRotZ",pvList,"=")
	if (numtype(g.wire.R[0] + g.wire.R[1] + g.wire.R[2]))	// for bad (or no) results, use zero rotation
		g.wire.R[0] = 0;	g.wire.R[1] = 0;	g.wire.R[2] = 0
	endif

	if (MicroGeometryBad(g))
		print "Not loading geometry from EPICS, values are invalid"
		DoAlert 0, "Not loading geometry from EPICS, values are invalid"
		return 1
	else
		GeometryUpdateCalc(g)									// calculate other values
		CopymicroGeometry(gIn,g)								// new geometry is valid, so copy it in
		if (ItemsInList(GetRTStackInfo(0))<3)					// print everything if run from command line
			printf "Loaded geometry information from EPICS on  %s,  %s\r",date(),time()
		endif
	endif
	return 0
End
//
Function/S protoEPICS_geo(str)									// proto function for get_mult_PV(), also contains default values for testing
	String str
	DoAlert 0, "EPICS not available, goemetry unchanged"
	return ""
End


Static Function putGeo2EPICS(gIn)	//fill the geometry structure from EPICS (uses caget)
	STRUCT microGeometry &gIn

	if (exists("EPICS_put_PV_num")!=6 || exists("EPICS_put_PV_str")!=6)
		DoAlert 0, "EPICS not available, nothing changed"
		return 1
	endif
	STRUCT microGeometry g										// a working copy of gIn
	CopymicroGeometry(g,gIn)										// copy gIn to g (the working copy)
	GeometryUpdateCalc(g)											// calculate other values
	if (MicroGeometryBad(g))
		printGeometry(g)
		print "The given microgeometry is BAD, not sending it to EPICS"
		DoAlert 0, "The given microgeometry is BAD, not sending it to EPICS"
		return 1
	endif
	if (g.Ndetectors < 1)
		DoAlert 0,"No detectors defined"
		return 1
	endif
	Variable i
	for (i=0;i<g.Ndetectors && !(g.d[i].used);i+=1)
	endfor
	if (i>=g.Ndetectors)
		DoAlert 0,"No detectors defined"
		return 1
	endif

	if (geoLocal_putNum(EPICS_PREFIX+"Ndetectors.RVAL",g.Ndetectors))
		return 1
	endif

	if (geoLocal_putNum(EPICS_PREFIX+"used0",g.d[0].used))
		return 1
	endif
	if (geoLocal_putNum(EPICS_PREFIX+"used1",g.d[1].used))
		return 1
	endif
	if (geoLocal_putNum(EPICS_PREFIX+"used2",g.d[2].used))
		return 1
	endif

	if (g.d[0].used)						// Detector 0 parameters
		putDetector2EPICS(g.d[0],0)	// put detector 0 values
	endif
	if (g.d[1].used)						// Detector 1 parameters
		putDetector2EPICS(g.d[1],1)	// put detector 1 values
	endif
	if (g.d[2].used)						// Detector 2 parameters
		putDetector2EPICS(g.d[2],2)	// put detector 2 values
	endif

	// Wire parameters
	if (geoLocal_putNum(EPICS_PREFIX+"wireDia",g.wire.dia))
		return 1
	endif
	if (geoLocal_putNum(EPICS_PREFIX+"wireOriginX",g.wire.origin[0]))
		return 1
	endif
	if (geoLocal_putNum(EPICS_PREFIX+"wireOriginY",g.wire.origin[1]))
		return 1
	endif
	if (geoLocal_putNum(EPICS_PREFIX+"wireOriginZ",g.wire.origin[2]))
		return 1
	endif
	if (geoLocal_putNum(EPICS_PREFIX+"wireAxisX",g.wire.axis[0]))
		return 1
	endif
	if (geoLocal_putNum(EPICS_PREFIX+"wireAxisY",g.wire.axis[1]))
		return 1
	endif
	if (geoLocal_putNum(EPICS_PREFIX+"wireAxisZ",g.wire.axis[2]))
		return 1
	endif
	if (geoLocal_putNum(EPICS_PREFIX+"wireKnife",g.wire.knife))
		return 1
	endif
	if (geoLocal_putNum(EPICS_PREFIX+"wireRotX",g.wire.R[0]))
		return 1
	endif
	if (geoLocal_putNum(EPICS_PREFIX+"wireRotY",g.wire.R[1]))
		return 1
	endif
	if (geoLocal_putNum(EPICS_PREFIX+"wireRotZ",g.wire.R[2]))
		return 1
	endif
	if (geoLocal_putNum(EPICS_PREFIX+"wireF",g.wire.F))
		return 1
	endif

	// Sample parameters
	if (geoLocal_putNum(EPICS_PREFIX+"SampleRotX",g.s.R[0]))
		return 1
	endif
	if (geoLocal_putNum(EPICS_PREFIX+"SampleRotY",g.s.R[1]))
		return 1
	endif
	if (geoLocal_putNum(EPICS_PREFIX+"SampleRotZ",g.s.R[2]))
		return 1
	endif
	if (geoLocal_putNum(EPICS_PREFIX+"SampleOriginX",g.s.O[0]))
		return 1
	endif
	if (geoLocal_putNum(EPICS_PREFIX+"SampleOriginY",g.s.O[1]))
		return 1
	endif
	if (geoLocal_putNum(EPICS_PREFIX+"SampleOriginZ",g.s.O[2]))
		return 1
	endif

	return 0
End
//
Static Function putDetector2EPICS(d,dnum)	// put detector values into EPICS, does not check d.used
	STRUCT detectorGeometry &d
	Variable dnum		// 0, 1, or 2

	String sd=num2istr(dnum)
	if (geoLocal_putNum(EPICS_PREFIX+"Nx"+sd,d.Nx))
		return 1
	endif
	if (geoLocal_putNum(EPICS_PREFIX+"Ny"+sd,d.Ny))
		return 1
	endif
	if (geoLocal_putNum(EPICS_PREFIX+"sizeX"+sd,d.sizeX * 1e-3))
		return 1
	endif
	if (geoLocal_putNum(EPICS_PREFIX+"sizeY"+sd,d.sizeY * 1e-3))
		return 1
	endif
	if (geoLocal_putNum(EPICS_PREFIX+"Rx"+sd,d.R[0]))
		return 1
	endif
	if (geoLocal_putNum(EPICS_PREFIX+"Ry"+sd,d.R[1]))
		return 1
	endif
	if (geoLocal_putNum(EPICS_PREFIX+"Rz"+sd,d.R[2]))
		return 1
	endif
	if (geoLocal_putNum(EPICS_PREFIX+"Px"+sd,d.P[0] * 1e-3))
		return 1
	endif
	if (geoLocal_putNum(EPICS_PREFIX+"Py"+sd,d.P[1] * 1e-3))
		return 1
	endif
	if (geoLocal_putNum(EPICS_PREFIX+"Pz"+sd,d.P[2] * 1e-3))
		return 1
	endif
	if (geoLocal_putStr(EPICS_PREFIX+"timeMeasured"+sd,d.timeMeasured))
		return 1
	endif
	if (geoLocal_putStr(EPICS_PREFIX+"geoNote"+sd,d.geoNote))
		return 1
	endif
	if (geoLocal_putStr(EPICS_PREFIX+"detectorID"+sd,d.detectorID))
		return 1
	endif
	if (geoLocal_putStr(EPICS_PREFIX+"distortionMapFile"+sd,d.distortionMapFile))
		return 1
	endif
End
//
Static Function  geoLocal_putNum(pv,value)	// utility function for EPICS_put_PV_num()
	String pv								// full PV name
	Variable value							// new number to set
	if (exists("EPICS_put_PV_num")!=6)
		DoAlert 0, "EPICS_put_PV_num not available, nothing changed"
		return 1
	endif
	FUNCREF protoEPICS_putNum  func=$"EPICS_put_PV_num"	//	Function/T EPICS_put_PV_num(pv,value)
	String errStr = func(pv,value)
	if (strlen(errStr))
		DoAlert 0, errStr
		print errStr
		return 1
	endif
	return 0
End
//
Function geoLocal_putStr(pv,value)	// EPICS_put_PV_num function for EPICS_put_PV_str()
	String pv								// full PV name
	String value							// new string value to set
	if (exists("EPICS_put_PV_str")!=6)
		DoAlert 0, "EPICS_put_PV_str not available, nothing changed"
		return 1
	endif
	FUNCREF protoEPICS_putStr  func=$"EPICS_put_PV_str"	//	Function/T EPICS_put_PV_num(pv,value)
	String errStr = func(pv,value)
	if (strlen(errStr))
		DoAlert 0, errStr
		print errStr
		return 1
	endif
	return 0
End
//
Function/S protoEPICS_putNum(pv,value)						// proto function for EPICS_put_PV_num()
	String pv							// full PV name
	Variable value						// new number to set
	DoAlert 0, "EPICS not available, nothing changed"
	return "ERROR -- EPICS not available, nothing changed"
End
//
Function/S protoEPICS_putStr(pv,str)							// proto function for EPICS_put_PV_str()
	String pv							// full PV name
	String str							// new string value to set
	DoAlert 0, "EPICS not available, nothing changed"
	return "ERROR -- EPICS not available, nothing changed"
End

// ================================= End of Geometry Sub-Panel ==================================
//=======================================================================================



//=======================================================================================
// ======================================= Start of Init =======================================

Function Init_microGeo()
	if (NumVarOrDefault("root:Packages:geometry:geoInited",0))
		return 0
	endif

	InitLatticeSymPackage()

if (stringmatch(systemUserName(),"tischler"))
	DoAlert 0,"Note to JZT:  need more consistency in wire.origin. it should be after keyence correction by before rotation"
endif
	NewDataFolder/O root:Packages						// ensure Packages exists
	NewDataFolder/O root:Packages:geometry			// ensure geometry exists
	NewDataFolder/O root:Packages:micro				// ensure micro exists
	Variable/G root:Packages:geometry:geoInited=0		// flag that initialization was called successsfully
	NVAR geoInited=root:Packages:geometry:geoInited

	Make/N=3/O/D root:Packages:geometry:pixel2q_ki, root:Packages:geometry:pixel2q_kout
	Make/N=3/O/D root:Packages:geometry:q2pixel_qhat,root:Packages:geometry:q2pixel_ki
	Make/N=3/O/D root:Packages:geometry:q2pixel_kout

//	Make/N=3/O/D root:Packages:geometry:pixel2depth_xyzDetector
//	Make/N=3/O/D root:Packages:geometry:GeoUpdate_ki

	Make/N=3/O/D root:Packages:geometry:depth_S, root:Packages:geometry:depth_delta
	Make/N=3/O/D root:Packages:geometry:depth_a, root:Packages:geometry:depth_nhat
	Make/N=3/O/D root:Packages:geometry:depth_wo, root:Packages:geometry:depth_C
	Make/N=3/O/D root:Packages:geometry:depth_sigma
	Make/N=(3,3)/O/D root:Packages:geometry:depth_rhoW, root:Packages:geometry:depth_mat
	Make/N=3/O/D root:Packages:geometry:depth_Rvec, root:Packages:geometry:depth_ki
	Make/N=3/O/D root:Packages:geometry:depth_pixel, root:Packages:geometry:depth_wc
	Make/N=(3,3)/O/D root:Packages:geometry:depth_rhoX

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
//Function LoadStandardDistortion()
//	String str= FunctionPath("LoadStandardDistortion")
//	if (strlen(str)<2)
//		return 1
//	endif
//	String fileName = ParseFilePath(1, str, ":", 1, 0)+"xymap.itx"
//
//	DoAlert 0, "Distortion not yet implemented"
//	printf "¥¥ Distortion not yet implemented\r"
//	return 1
//End

Function MakeDistortionMap()
	if (Exists("root:Packages:geometry:xymap")!=1)
		return 1
	endif
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
	return 0
End
//
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

// ======================================== End of Init =======================================
//=======================================================================================