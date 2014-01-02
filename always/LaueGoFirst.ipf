#pragma rtGlobals= 2
#pragma version = 2.03
#pragma ModuleName = LaueGoFirst
#include "Utility_JZT", version>=3.16
#pragma hide = 1
Static Constant JZTalwaysFirst_Version_Min=2.4	// minimum required vesion of "always first.ipf"

StrConstant MICRO_GEOMETRY_VERSION_PATH = "root:Packages:MICRO_GEOMETRY_VERSION"
//		bit flags for MICRO_GEOMETRY_VERSION
//			0 -> OLD
//			1 -> HDF5
//			2 -> TIFF
//			4 -> SPE


Static Function IgorStartOrNewHook(IgorApplicationNameStr)
	String IgorApplicationNameStr
	InitLaueGoFirst()
	return 0
End



//Strconstant browserLeadJZT="\r--> moved to ", browserTailJZT="\r\r"
//ModifyBrowser echoCommands=0, command3="print browserLeadJZT,GetDataFolder(1),browserTailJZT"
Menu "Analysis",dynamic
	Submenu "Packages"
		SubMenu "<B<ULaueGo [� beam]"

			// none
			microMenuShowPreliminary(LaueGoFirst#arrow("B")+"Select Image Types"), LaueGoFirst#SetMicroPreferedImageTypes()
			help = {"Select image types for LaueGo"}

			// NEW
			microMenuShowN(LaueGoFirst#arrow("")+"LaueGo Panel"),MakeMicroPanel(-1)
			help = {"new starting point for the generic micro beam analysis (LaueGo"}
			microMenuShowN("Indexing & Peak Searching"),Execute/P "INSERTINCLUDE  \"IndexingN\", version>=4.00";Execute/P "COMPILEPROCEDURES ";Execute/P "InitIndexingPackage()"
			//	help = {"Find and Fit Peaks in an image, and then Index the peaks"}
			microMenuShowN("Energy-Wire Scans"),Execute/P "INSERTINCLUDE  \"EnergyWireScansN\", version>=1.000";Execute/P "COMPILEPROCEDURES "
			//	help = {"Processes Energy-Wire scans into Depth-Q plots"}
			microMenuShowN("Index Lots of Images"),Execute/P "INSERTINCLUDE  \"IndexLotsN\", version>=2.27";Execute/P "COMPILEPROCEDURES "
			//	help = {"index a large group of images, like from a wire-scan"}
			microMenuShowN("3d-Grains with Gizmo"),Execute/P "INSERTINCLUDE  \"GizmoGrains\", version>=1.3";Execute/P "COMPILEPROCEDURES ";Execute/P "initGrainGizmo()"
			//	help = {"Processes 3D array of matricies into a viewable volume"}
			microMenuShowN("3d-Array of Orientation Mats"),Execute/P "INSERTINCLUDE  \"ArrayOf3dOrientsN\", version>=2.56";Execute/P "COMPILEPROCEDURES "
			//	help = {"Processes 3D array of matricies into viewable of rotataion angles, dislocation tensors, GND, etc."}
			microMenuShowN("Pole Figure"),Execute/P "INSERTINCLUDE  \"PoleFigure\"";Execute/P "COMPILEPROCEDURES "
			//	help = {"Loads Pole figure Macros"}
			microMenuShowN("Details of Depth Resolved Wire Scans"),Execute/P "INSERTINCLUDE  \"DepthResolvedQueryN\"";Execute/P "COMPILEPROCEDURES "
			//	help = {"used to analyze a depth resolved wire scan, mostly for debugging"}

			// OLD
			microMenuShowOLD("��> LaueGo Panel"),MakeMicroPanel(-1)
			help = {"old starting point for the generic micro beam analysis"}
			microMenuShowOLD("Indexing & Peak Searching"),Execute/P "INSERTINCLUDE  \"Indexing\", version>=2.59";Execute/P "COMPILEPROCEDURES ";Execute/P "InitIndexingPackage()"
			//	help = {"Find and Fit Peaks in an image, and then Index the peaks"}
			microMenuShowOLD("Energy-Wire Scans"),Execute/P "INSERTINCLUDE  \"EnergyWireScans\", version>=0.9944";Execute/P "COMPILEPROCEDURES "
			//	help = {"Processes Energy-Wire scans into Depth-Q plots"}
			microMenuShowOLD("Index Lots of Images"),Execute/P "INSERTINCLUDE  \"IndexLots\", version>=2.17";Execute/P "COMPILEPROCEDURES "
			//	help = {"index a large group of images, like from a wire-scan"}
			microMenuShowOLD("3d-Grains with Gizmo"),Execute/P "INSERTINCLUDE  \"GizmoGrains\", version>=1.3";Execute/P "COMPILEPROCEDURES ";Execute/P "initGrainGizmo()"
			//	help = {"Processes 3D array of matricies into a viewable volume"}
			microMenuShowOLD("3d-Array of Orientation Mats"),Execute/P "INSERTINCLUDE  \"ArrayOf3dOrients\", version>=2.41";Execute/P "COMPILEPROCEDURES "
			//	help = {"Processes 3D array of matricies into viewable of rotataion angles, dislocation tensors, GND, etc."}
			microMenuShowOLD("Pole Figure"),Execute/P "INSERTINCLUDE  \"PoleFigure\"";Execute/P "COMPILEPROCEDURES "
			//	help = {"Loads Pole figure Macros"}
			microMenuShowOLD("Details of Depth Resolved Wire Scans"),Execute/P "INSERTINCLUDE  \"DepthResolvedQuery\"";Execute/P "COMPILEPROCEDURES "
			//	help = {"used to analyze a depth resolved wire scan, mostly for debugging"}
			microMenuShowN("Laue Simulation"),Execute/P "INSERTINCLUDE  \"LaueSimulation\", version>=1.06";Execute/P "COMPILEPROCEDURES "

			Submenu "Calibration Procedures"
				microMenuShowN("Detector Calibration"),Execute/P "INSERTINCLUDE  \"DetectorCalibration\"";Execute/P "COMPILEPROCEDURES "
				help = {"Used for the Initial Calibration of the Perkin-Elmer detectors"}
				"micro-mono Calibration",Execute/P "INSERTINCLUDE  \"monoCalibrate\", version>=1.6";Execute/P "COMPILEPROCEDURES ";Execute/P "monoCalibrateInitPackage()"
				help = {"Used to calibrate energy of mono at 34ID-D"}
			End

		End
	End
End


Menu "Analysis"
	Submenu "Packages"
		// "� beam, Sector 34" SubMenu gets insterted below
		"Stereographic Projections",Execute/P "INSERTINCLUDE  \"StereographicProjection\", version>=2.81";Execute/P "COMPILEPROCEDURES ";Execute/P "InitStereoGraphicPackage()"
		SubMenu "EPICS"
			"EPICS("
			"PV I\O", Execute/P "INSERTINCLUDE \"epics\"" ; 	Execute/P "COMPILEPROCEDURES ";Execute/P "epicsInitPackage()"
			help = {"Load procedures for talking to PVs via EPICS (only useful at APS)"}
			"Load Scan Record", Execute/P "INSERTINCLUDE \"LoadEPICSscans\"" ; 	Execute/P "COMPILEPROCEDURES "
			help = {"Load procedures for Loading the dump of a scan record using TkGrab, does not need EPICS support"}
		End
	End
End
//
Function/T microMenuShowN(str)
	String str
	return SelectString(NumVarOrDefault(MICRO_GEOMETRY_VERSION_PATH,256)&5,"",str)		// for both spe AND hdf5
End
Function/T microMenuShowOLD(str)
	String str
	return SelectString(NumVarOrDefault(MICRO_GEOMETRY_VERSION_PATH,NaN)==0,"",str)
End
Function/T microMenuShowPreliminary(str)
	String str
	return SelectString(exists(MICRO_GEOMETRY_VERSION_PATH)!=2,"",str)
End


Menu "Data"
	Submenu "Packages"
		"-"
		"EPICS Load Scan Record", Execute/P "INSERTINCLUDE \"LoadEPICSscans\"" ; 	Execute/P "COMPILEPROCEDURES "
		help = {"Load procedures for Loading the dump of a scan record using TkGrab, does not need EPICS support"}
		"3d-Grains with Gizmo",Execute/P "INSERTINCLUDE  \"GizmoGrains\", version>=1.3";Execute/P "COMPILEPROCEDURES ";Execute/P "initGrainGizmo()"
		help = {"Load procedures looking at 3-d grain structures"}
	End
End
Menu "Load Waves"
	Submenu "Packages"
	"-"
	"EPICS Load Scan Record", Execute/P "INSERTINCLUDE \"LoadEPICSscans\"" ; 	Execute/P "COMPILEPROCEDURES "
	help = {"Load procedures for Loading the dump of a scan record using TkGrab, does not need EPICS support"}
	End
End


Menu "Graph"
	"Set Aspect Ratio to Get Square Pixels \ Range",SetAspectToSquarePixels("")
End

Menu "Edit"
	LaueGoFirst#MenuItemIfWindowTypes("Copy Window Info",1+2+64+4096), /Q,  LaueGoFirst#GetWindowInfo2Scrap()
	SubMenu LaueGoFirst#MenuItemIfScrapValidWindowInfo("  Paste Window","type")
		LaueGoFirst#MenuItemIfScrapValidWindowInfo("Paste Window Size","left"), /Q, LaueGoFirst#PutSizeWindow("")
		LaueGoFirst#MenuItemIfScrapValidWindowInfo("Paste Graph Axis Ranges","axis_left_min"), /Q, LaueGoFirst#PutScrapGraphAxis("")
		"-"
		LaueGoFirst#MenuItemIfScrapValidWindowInfo("Paste Gizmo Quaternion","quaternion"), /Q, LaueGoFirst#PutGizmoQuaternion("")
	End
End


Static Structure microImageTypePrefs
	int16	old
	int16	hdf5
	int16	tiff
	int16	spe
EndStructure

Static Function SetMicroPreferedImageTypes()
	STRUCT microImageTypePrefs prefs
	prefs.old = 0
	prefs.hdf5 = 0
	prefs.tiff = 0
	prefs.spe = 0
	LoadPackagePreferences/MIS=1 "microGeo","microGeoNPrefs",1,prefs
	if (V_bytesRead<8)
		prefs.old = 0
		prefs.hdf5 = 0
		prefs.tiff = 0
		prefs.spe = 0
	endif
	Variable old,hdf5,tiff,spe, repeat
	Prompt old,"Old Style",popup,"---;Old Style, for OLD data"
	Prompt hdf5,"HDF5 images",popup,"---;HDF5 Images"
	Prompt tiff,"TIFF images",popup,"---;TIFF Images in 4500S"
	Prompt spe,"WinView images",popup,"---;WinView (spe) Images"
	//DoAlert 0,"Select either 'Old Style' or some combination of the others"
	do
		old = prefs.old + 1
		hdf5 = prefs.hdf5 + 1
		tiff = prefs.tiff + 1
		spe = prefs.spe + 1
		DoPrompt "Desired Image Types",old,hdf5,tiff,spe
		if (V_flag)
			return -1
		endif
		prefs.old = old-1
		prefs.hdf5 = hdf5-1
		prefs.tiff = tiff-1
		prefs.spe = spe-1
		repeat = 0
		if ((prefs.hdf5 + prefs.tiff + prefs.spe) && prefs.old)
			DoAlert 0,"You cannot select 'Old Style' and anything else"
			repeat = 1
		elseif (!(prefs.hdf5 + prefs.tiff + prefs.spe) && !prefs.old)
			DoAlert 0, "You did not select any image type"
			repeat = 1
		endif
	while (repeat)
	SavePackagePreferences/FLSH=1 "microGeo","microGeoNPrefs",1,prefs
	setMICRO_GEOMETRY_VERSION_PATH()
End

Static Function setMICRO_GEOMETRY_VERSION_PATH([preset])
	Variable preset						// Optional, a bit flag 0=OLD, 1=HDF5, 2=Tiff, 3=SPE, if old, then nothing else allowed
	if (exists(MICRO_GEOMETRY_VERSION_PATH)==2)
		return 0
	endif

	STRUCT microImageTypePrefs prefs
	if (ParamIsDefault(preset))				// no preset given, use stored default values
		LoadPackagePreferences/MIS=1 "microGeo","microGeoNPrefs",1,prefs
		if (V_bytesRead<8)
			DoAlert 0,"Unable to read prefered images types, defaulting to only HDF5"
			prefs.old = 0
			prefs.hdf5 = 1
			prefs.tiff = 0
			prefs.spe = 0
		endif
	else
		prefs.old = (preset == 0)				// old cannot coesist with any of the others
		prefs.hdf5 = (preset & 1)
		prefs.tiff = (preset & 2)
		prefs.spe = (preset & 4)
//		SavePackagePreferences/FLSH=1 "microGeo","microGeoNPrefs",1,prefs	// save the defaults
	endif

	NewDataFolder/O root:Packages
	Variable/G $MICRO_GEOMETRY_VERSION_PATH
	NVAR type=$MICRO_GEOMETRY_VERSION_PATH
	if (prefs.old)
		type = 0
	else
		type += prefs.hdf5 ? 1 : 0
		type += prefs.tiff ? 2 : 0
		type += prefs.spe ? 4 : 0
	endif
	Variable/G $MICRO_GEOMETRY_VERSION_PATH=type
	String insert

	if (type==0)
		insert = "INSERTINCLUDE  \"microGeometry\", version>=2.61"
	elseif (type&2)
		insert = "INSERTINCLUDE  \"microGeometryN\", version>=1.61"
	elseif (type&5)
		insert = "INSERTINCLUDE  \"microGeometryN\", version>=1.61"
	else
		return -1
	endif
	Execute/P insert
	Execute/P "COMPILEPROCEDURES "
	Execute/P/Q "Init_microGeo()"
	Execute/P/Q "MakeMicroPanel(-1)"
	return 0
End


Static Function/T arrow(fmt)
	String fmt
	if (strsearch(IgorInfo(2),"Macintosh",0))
		return " -->"			// for Windows, note leading space in string is necessary
	endif
	String out=""
	Variable i
	for (i=0;i<strlen(fmt);i+=1)
		out += "<"+fmt[i]
	endfor
	out += "��> "
	return out
End






//  ====================================================================================  //
//  ========================== Start of Window Copy/Paste Info =========================  //

Static Function/T MenuItemIfWindowTypes(item,flags)
	String item
	Variable flags
	flags = numtype(flags) ? 0 : flags
	if (ItemsInLIst(WinList("*",";","WIN:"+num2istr(flags)))<1)
		return "("+item					// top graph does not contain an image, so disable menu item
	endif
	return item
End
//
Static Function/T MenuItemIfScrapValidWindowInfo(item,requiredKey)
	String item
	String requiredKey
	String scrap=GetScrapText()
	String type=StringByKey("type",scrap)
	if (!stringmatch(type,"*WindowInfo"))
		return "("+item					// Clipboard does not contain right kind of text
	endif

	if (strlen(requiredKey))
		if (strlen(StringByKey(requiredKey,scrap))<1)
			return "("+item					// the required key is missing
		endif
	endif

	Variable flag = NumberByKey("flag",scrap)
	String gName = StringFromList(0,WinList("*",";","WIN:"+num2istr(1+2+64+4096)))
	Variable topType = WinType(gName)
	topType = topType==7 ? 64 : topType		// Panel should be 64
	topType = topType==13 ? 4096 : topType	// XOP window should be 4096

	if (flag == topType)
		return item
	endif
	return "("+item					// top window does not match kind of information in clipboard
End


Static Function/T GetWindowInfo2Scrap()
	Variable flags = 1+2+64+ (exists("NewGizmo")==4 ? 4096 : 0)
	String wName=StringFromLIst(0,WinList("*",";","WIN:"+num2istr(flags))), type=""
	Variable graph=0, gizmo=0, table=0, panel=0, flag=0
	if (WhichListItem(wName, WinList("*",";","WIN:1"))>=0)
		graph = 1
		flag = 1
		type = "Graph"
	elseif (WhichListItem(wName, WinList("*",";","WIN:2"))>=0)
		table = 1
		flag = 2
		type = "Table"
	elseif (WhichListItem(wName, WinList("*",";","WIN:64"))>=0)
		panel = 1
		flag = 64
		type = "Panel"
	elseif (WhichListItem(wName, WinList("*",";","WIN:4096"))>=0 && flags>=4096)	// a Gizmo/XOP
		type = "Gizmo"
		flag = 4096
		Execute "GetGizmo gizmoNameList"
		String gList=StrVarOrDefault("S_GizmoNames","")
		KillStrings/Z S_GizmoNames
		if (WhichListItem(wName, gList)<0)
			wName = ""
		else
			gizmo = 1
		endif
	endif
	if (strlen(wName)<1)
		DoAlert 0,"There are no Graphs or Gizmos to get info from"
		return ""
	endif

	String keys=ReplaceStringByKey("type","",type+"WindowInfo")
	keys = ReplaceStringByKey("wName",keys,wName)
	keys = ReplaceNumberByKey("flag",keys,flag)
	if (gizmo)
		Execute "GetGizmo winPixels"			// get window position & size
		Variable left=NumVarOrDefault("V_left",NaN), right=NumVarOrDefault("V_right",NaN)
		Variable top=NumVarOrDefault("V_top",NaN), bottom=NumVarOrDefault("V_bottom",NaN)
		KillVariables/Z V_left, V_right, V_top, V_bottom
		if (numtype(top+bottom+left+right))
			DoAlert 0, "Unable to get Size of Gizmo Window"
			return ""
		endif
		keys = ReplaceNumberByKey("top",keys,top)
		keys = ReplaceNumberByKey("left",keys,left)
		keys = ReplaceNumberByKey("right",keys,right)
		keys = ReplaceNumberByKey("bottom",keys,bottom)

		Execute "GetGizmo curRotation"
		Variable euler0=NumVarOrDefault("GizmoEulerA",NaN), euler1=NumVarOrDefault("GizmoEulerB",NaN), euler2=NumVarOrDefault("GizmoEulerC",NaN)
		KillVariables/Z GizmoEulerA,GizmoEulerB,GizmoEulerC
		if (numtype(euler0+euler1+euler2))
			DoAlert 0, "Unable to get Orientation from Gizmo"
			return ""
		endif
		Wave quaternion=Euer2Quaternion(euler0*PI/180, euler1*PI/180,euler2*PI/180)
		String str
		sprintf str, "{%.6f,%.6f,%.6f,%.6f}",quaternion[0],quaternion[1],quaternion[2],quaternion[3]
		keys = ReplaceStringByKey("quaternion",keys,str)

	elseif(graph || table || panel)
		GetWindow/Z kwTopWin, wsize
		if (numtype(V_top+V_bottom+V_left+V_right) || V_flag)
			DoAlert 0, "Unable to get Size of Gizmo Window"
			return ""
		endif
		keys = ReplaceNumberByKey("top",keys,V_top)
		keys = ReplaceNumberByKey("left",keys,V_left)
		keys = ReplaceNumberByKey("right",keys,V_right)
		keys = ReplaceNumberByKey("bottom",keys,V_bottom)
		if(graph)				// for a graph, also get the axis ranges for top, bottom, left, right
			String axis, aList=AxisList("")
			Variable i
			for (i=0;i<ItemsInLIst(aList);i+=1)
				axis = StringFromList(i,aList)
				GetAxis/Q $axis
				if (numtype(V_min+V_max)==0)
					keys = ReplaceNumberByKey("axis_"+axis+"_min",keys,V_min)
					keys = ReplaceNumberByKey("axis_"+axis+"_max",keys,V_max)
				endif
			endfor
		endif
	else
		keys = ""
	endif

	if (strlen(keys)>0)
		PutScrapText keys
	endif
	return keys
End
//
Static Function/WAVE Euer2Quaternion(euler0, euler1, euler2)		// make Quaternion from Euler angles
	Variable euler0, euler1, euler2						// Euler angles (radian)
	// Basically we create 3 Quaternions, one for each euler angle
	// and multiply those together.  The calculation below does the same, just shorter
	Variable cos0=cos(euler0/2), sin0=sin(euler0/2)
	Variable cos1=cos(euler1/2), sin1=sin(euler1/2)
	Variable cos2=cos(euler2/2), sin2=sin(euler2/2)
	Make/N=4/FREE quaternion
	quaternion[0] = sin0*cos1*cos2 - cos0*sin1*sin2
	quaternion[1] = cos0*sin1*cos2 + sin0*cos1*sin2
	quaternion[2] = cos0*cos1*sin2 - sin0*sin1*cos2
	quaternion[3] = cos0*cos1*cos2 + sin0*sin1*sin2
	quaternion[0,2] *= -1
	normalize(quaternion)
	return quaternion
End


Static Function PutSizeWindow(scrap)
	String scrap
	if (strlen(scrap)==0)
		scrap=GetScrapText()
	endif
	Variable graph=0, table=0, panel=0, gizmo=0, flag=0
	String type=""
	strswitch(StringByKey("type",scrap))
		case "GraphWindowInfo":
			type = "Graph"
			graph = 1
			flag = 1
			break
		case "TableWindowInfo":
			type = "Table"
			table = 1
			flag = 2
			break
		case "PanelWindowInfo":
			type = "Panel"
			panel = 1
			flag = 64
			break
		case "GizmoWindowInfo":
			type = "Gizmo"
			gizmo = 1
			flag = 4096
			break
		default:
			DoAlert 0, "Did Not find valid window info in Clipboard"
			return 1
	endswitch
 
	String topWin=StringFromList(0,WinList("*",";","")) 
	String gName = StringFromList(0,WinList("*",";","WIN:"+num2istr(flag)))
	if (!stringmatch(gName,topWin))
		DoAlert 0, "Clipboard info is from "+type+", which does not match top window type"
		return 1
	endif 
  
 	Variable left=NumberByKey("left",scrap), right=NumberByKey("right",scrap), top=NumberByKey("top",scrap), bottom=NumberByKey("bottom",scrap)
	if (numtype(top+bottom+left+right))
		DoAlert 0,"Unable to retrieve desired Window Size, It was not yet saved."
		return 1
	endif
	Variable width=right-left, height=bottom-top		// get desired width and height

	top = NaN
	left = NaN
	if (gizmo)
		Execute "GetGizmo gizmoName"
		gName=StrVarOrDefault("S_GizmoName","")
		KillStrings/Z S_GizmoName
		if (strlen(gName)<1)
			DoAlert 0,"No Gizmo Visible to Re-Size"
			return 1
		endif
		Execute "GetGizmo winPixels"			// get current window position
		left = NumVarOrDefault("V_left",NaN)
		top = NumVarOrDefault("V_top",NaN)
		KillVariables/Z V_left, V_right, V_top, V_bottom
	else
		if (strlen(gName)<1)
			DoAlert 0,"No "+type+" Visible to Re-Size"
			return 1
		endif
		GetWindow/Z kwTopWin, wsize
		if (V_flag==0)
			top = V_top
			left = V_left
		endif
	endif
	if (numtype(left+top))
		DoAlert 0,"Unable to retrieve position of top "+type+"!"
		return 1
	endif
	MoveWindow/W=$gName left, top, left+width, top+height	// change window size leaving top/left corner where it is
	return 0
End


Static Function PutScrapGraphAxis(scrap)
	String scrap
	if (strlen(scrap)==0)
		scrap=GetScrapText()
	endif

	if (!stringmatch(StringByKey("type",scrap),"GraphWindowInfo"))
		DoAlert 0, "Did Not find valid window info in Clipboard"
		return 1
	endif

	String topWin=StringFromList(0,WinList("*",";","")) 
	String gName = StringFromList(0,WinList("*",";","WIN:1"))
	if (!stringmatch(gName,topWin))
		DoAlert 0, "Clipboard info is from 'Graph', which does not match top window type"
		return 1
	endif 
	String axis, aList=AxisList("")
	Variable N=ItemsInList(aList)
	if (N<1)
		DoAlert 0,"Unable to set any Graph Axes!"
		return 1
	endif
  
	Variable i, lo,hi
	for (i=0;i<ItemsInLIst(aList);i+=1)			// for each axis on this graph
		axis = StringFromList(i,aList)
		lo = NumberByKey("axis_"+axis+"_min",scrap)
		hi = NumberByKey("axis_"+axis+"_max",scrap)
		GetAxis/Q $axis
		if (numtype(V_min+V_max+lo+hi)==0)	// true when scrap has info about this axis
			SetAxis $axis lo,hi
		endif
	endfor
	return 0
End


Static Function PutGizmoQuaternion(scrap)		// Set Quaternion
	String scrap
	if (strlen(scrap)==0)
		scrap=GetScrapText()
	endif
	Variable gizmo=stringmatch(StringByKey("type",scrap),"GizmoWindowInfo")
	String quaternion=StringByKey("quaternion",scrap)
	if (!gizmo || strlen(quaternion)<1)
		DoAlert 0,"Clipboard does not contain a Gizmo Quaternion Information"
		return 1
	endif
	Execute "GetGizmo gizmoName"
	String gName=StrVarOrDefault("S_GizmoName","")
	KillStrings/Z S_GizmoName
	if (strlen(gName)<1)
		DoAlert 0,"No Gizmo Visible to Rotate"
		return 1
	endif
	Execute "ModifyGizmo/N="+gName+" SETQUATERNION="+quaternion
	return 0
End


Function SquareUpGizmo(gName)
	String gName
	if (strlen(gName)<1)
		Execute "GetGizmo gizmoName"
		gName=StrVarOrDefault("S_GizmoName","")
		KillStrings/Z S_GizmoName
	endif
	if (strlen(gName)<1)
		return 1
	endif

	Execute "GetGizmo winPixels"			// get window position & size
	Variable left=NumVarOrDefault("V_left",NaN), right=NumVarOrDefault("V_right",NaN)
	Variable top=NumVarOrDefault("V_top",NaN), bottom=NumVarOrDefault("V_bottom",NaN)
	KillVariables/Z V_left, V_right, V_top, V_bottom
	Variable height=bottom-top
	// Variable width=right-left, height=bottom-top
	// printf "[top=%g,  bottom=%g,  �=%g],   [left=%g,  right=%g,  �=%g]\r"top,bottom,height,  left,right,width
	top = max(44,top)
	MoveWindow/W=$gName left, top, left+height, top+height
End

//  ========================== End of Window Copy/Paste Info ===========================  //
//  ====================================================================================  //



//  ====================================================================================  //
//  ============================== Start of Square Pixels ==============================  //

Function SetAspectToSquarePixels(gName)
	String gName										// name of the graph, use "" for the top graph
	Variable printIt = strlen(GetRTStackInfo(2))<=0
	if (strlen(gName)<1)
		gName = StringFromList(0,WinList("*",";","WIN:1"))
	endif
	if (WinType(gName)!=1)
		if (printIt)
			DoAlert 0, "ERROR, in SetAspectToSquarePixels(), '"+gName+"' is not an graph"
		endif
		return NaN										// if no image on graph, do not try to set aspect ratio
	endif

	GetAxis/W=$gName/Q bottom
	if (V_flag)											// if no bottom, try top
		GetAxis/W=$gName/Q top
	endif
	if (V_flag)
		if (printIt)
			DoAlert 0, "ERROR, SetAspectToSquarePixels(), unable to get size of vertical axis"
		endif
		return NaN
	endif
	Variable width = abs(V_max-V_min)

	GetAxis/W=$gName/Q left
	if (V_flag)											// if no left, try right
		GetAxis/W=$gName/Q right
	endif
	if (V_flag)
		if (printIt)
			DoAlert 0, "ERROR, SetAspectToSquarePixels(), unable to get size of horizontal axis"
		endif
		return NaN
	endif
	Variable height = abs(V_max-V_min)

	Variable aspect = height / width
	//	printf "size ,  width = %g,  height = %g,   aspect = height/width = %g\r",width,height,aspect
	if (numtype(aspect) || aspect<=0)
		return NaN
	elseif (aspect<1)
		ModifyGraph/W=$gName height={Aspect,aspect}, width=0
	elseif (aspect>=1)
		ModifyGraph/W=$gName width={Aspect,1/aspect}, height=0
	endif

	return aspect
End

//  =============================== End of Square Pixels ===============================  //
//  ====================================================================================  //



//  ====================================================================================  //
//  ============================== Start of Wave Printing ==============================  //

ThreadSafe Function printWave(w,[name,brief])		// print a wave (vector or matrix) to history
	Wave w
	String name										// optional user supplied name to use
	Variable brief										// print in briefer form
	if (ParamIsDefault(name))
		name = NameOfWave(w)
	endif
	brief = ParamIsDefault(brief) ? 0 : !(!brief)
	if (!WaveExists(w))
		print "in 'printWave', wave does not exist"
		// DoAlert 0, "in 'printWave', wave does not exist"
		return 1
	endif

	if (DimSize(w, 1)<=1)		// for vectors
		printvec(w,name=name)
	elseif (DimSize(w, 2)==0)	// for 2-d matrix
		if (DimSize(w,0)<=1 || DimSize(w,1)<=1)
			printvec(w,name=name)
		else
			return printmat(w,name=name,brief=brief)
		endif
	else
		print "cannot yet handle dimensions 3 or 4"
	endif
	return 0
End
//
ThreadSafe Static Function printvec(w,[name])		// print a vector to screen
	Wave w
	String name										// optional user supplied name to use
	if (ParamIsDefault(name))
		name = NameOfWave(w)
	endif

	if (strlen(name))
		printf "%s = %s\r", name,vec2str(w)
	else
		printf "%s\r", vec2str(w)
	endif
End
//
ThreadSafe Static Function printmat(m,[name,brief,rowMax])
	Wave m
	String name										// optional user supplied name to use
	Variable brief										// print in briefer form
	Variable rowMax									// maximum number of rows to print
	if (ParamIsDefault(name))
		name = NameOfWave(m)
	endif
	rowMax = ParamIsDefault(rowMax) ? 50 : rowMax
	rowMax = (rowMax>0) ? rowMax : 50
	if (DimSize(m, 1)==0 || DimSize(m,2)!=0)	// for 2-d matrix only
		print "Can only print 2-d matricies with printmat"
		// DoAlert 0, "Can only print 2-d matricies with printmat"
		return 1
	endif

	if (brief && strlen(name))
		printf "%s:\r",name
	endif
	Variable Nrow=DimSize(m,0), row
	Nrow = min(Nrow,rowMax)
	for (row=0;row<Nrow;row+=1)
		if (WaveType(m) %& 0x01)					// true for complex numbers
			print printmatOneListComplex(m,row,name=name,brief=brief)
		else
			print printmatOneListReal(m,row,name=name,brief=brief)
		endif
	endfor
	if (DimSize(m,0)>Nrow)
		print "      ."
		print "      ."
		printf "      .\t\t printed only %d of the %d rows\r",Nrow,DimSize(m,1)
	endif
//	if (DimSize(mw,0)>Nrow || DimSize(mw,1)>Nrow)
//		printf "Only printed part of the (%d x %d) matrix\r",DimSize(mw,0),DimSize(mw,1)
//	endif
	return 0
End
//
ThreadSafe Static Function/T printmatOneListReal(m,row,[name,brief])// print one line for real (not complex) matricies
	Wave m
	Variable row										// row number (starts with 0)
	String name										// optional user supplied name to use
	Variable brief										// print in briefer form
	if (ParamIsDefault(name))
		name = NameOfWave(m)
	endif

	String line="", str
	Variable j, Ncol=DimSize(m,1)
	for (j=0;j<Ncol;j+=1)
		if (strlen(line)>100)
			line += "  ..."
			break
		elseif (brief)
			sprintf str, "%g    ",m[row][j]
		else
			sprintf str, "%s[%d][%d] = %g;    ",name,row,j,m[row][j]
		endif
		line += str
	endfor
	line = line[0,strlen(line)-4-1]					// strip off trailing 4 spaces
	return line
End
//
ThreadSafe Static Function/T printmatOneListComplex(m,row,[name,brief])// print one line for complex (not real) matricies
	Wave/C m
	Variable row										// row number (starts with 0)
	String name										// optional user supplied name to use
	Variable brief										// print in briefer form
	if (ParamIsDefault(name))
		name = NameOfWave(m)
	endif
	String line="", str
	Variable j, Ncol=DimSize(m,1)
	for (j=0;j<Ncol;j+=1)
		if (strlen(line)>100)
			line += "  ..."
			break
		elseif (brief)
			sprintf str, "(%g,%g)    ",real(m[row][j]),imag(m[row][j])
		else
			sprintf str, "%s[%d][%d] = (%g,%g);    ",name,row,j,real(m[row][j]),imag(m[row][j])
		endif
		line += str
	endfor
	line = line[0,strlen(line)-4-1]					// strip off trailing 4 spaces
	return line
End
//
ThreadSafe Function/T vec2str(w,[places,maxPrint,bare,sep])		// convert vector to s string suitable for printing, does not include name
	Wave w										// 1d wave to print
	Variable places								// number of places, for default, use negative or NaN
	Variable maxPrint							// maximum number of elements to print, defaults to 20
	Variable bare								// if bare is TRUE, then suppress the "{}" in the output
	String sep									// optional separator, default is ",  "   a comma and 2 spaces

	maxPrint = ParamIsDefault(maxPrint) ? 20 : maxPrint
	maxPrint = maxPrint>0 ? maxPrint : 20
	places = ParamIsDefault(places) ? -1 : places
	bare = ParamIsDefault(bare) ? 0 : !(!bare)
	sep = SelectString(ParamIsDefault(sep),sep,",  ")

	if (!WaveExists(w))
		return SelectString(bare,"{}","")
	endif

	Wave/T tw=$GetWavesDataFolder(w,2)
	Wave/C cw=$GetWavesDataFolder(w,2)
	Variable waveIsComplex = WaveType(w) %& 0x01
	Variable numeric = (WaveType(w)!=0)

	String fmt
	if (waveIsComplex)
		places = places>=0 ? min(20,places) : 5	// default to 5 for unacceptable values
		sprintf fmt,"(%%.%dg, %%.%dg)",places,places
	elseif (numeric)
		places = places>=0 ? min(20,places) : 5	// default to 5 for unacceptable values
		sprintf fmt,"%%.%dg",places
	elseif (places>0)								// must be string, and a maximum length given
		sprintf fmt, "\"%d%%s\"",places
	else												// string with no preferred length
		fmt = "\"%%s\""
	endif

	Variable i=0, n
	n = numpnts(w)
	maxPrint = min(n,maxPrint)
	String str, out=SelectString(bare,"{","")

	do
		if (waveIsComplex)						// a complex wave
			sprintf str,fmt, real(cw[i]),imag(cw[i])
		elseif (numeric && (!waveIsComplex))	// a simple number wave
			sprintf str,fmt, w[i]
		elseif (!numeric)						// a text wave
			sprintf str,"\"%s\"", tw[i]
		endif
		out += str
		if (i<(n-1))
			sprintf str,sep
			out += str
		endif
		i += 1
	while (i<maxPrint)
	if (n>maxPrint)
		sprintf str,"...}\ronly printed %d of %d values\r",maxPrint,n
		out += str
	else
		out += SelectString(bare,"}","")
	endif
	return out
End

//  =============================== End of Wave Printing ===============================  //
//  ====================================================================================  //



//  ====================================================================================  //
//  =========================== Start of Generic Graph Style ===========================  //

Proc Generic_Graph_Style() : GraphStyle
	Silent 1
//	LaueGoFirst#DoGenericGraphStyle("")
	GenericGraphStyle("")
EndMacro
//
Function GenericGraphStyle(gName)
	String gName					// graph name, usually "", for top graph
	if (ItemsInList(WinList("*",";","WIN:1"))<1)
		return 1
	elseif (exists("GenericGraphStyleLocal")==5)
		Execute "GenericGraphStyleLocal(\""+gName+"\")"		// in case users GenericGraphStyleLocal("") is a Macro
	else
		FUNCREF GenericGraphStyleTemplate styleFunc=$"GenericGraphStyleLocal"
		styleFunc(gName)											// defaults to GenericGraphStyleTemplate(gName)
	endif
	return 0
End
//
//	This can be overridden by making a Function named  GenericGraphStyleLocal(gName)
//
Function GenericGraphStyleTemplate(gName)				// default generic graph style funciton
	String gName					// graph name, usually ""
		//	if (NumberByKey("IGORVERS",IgorInfo(0))>=5.01)
		//	//	ModifyGraph gfRelSize=4
		//	//		ModifyGraph/W=$gName gfMult=130
		//	//		ModifyGraph/W=$gName gfMult=110		// just remove all of the gfMult, use the Igor defaults
		//	else
		//		ModifyGraph/W=$gName gfSize=18
		//	endif
	ModifyGraph/Z/W=$gName tick=2, minor=1, standoff=1
	ModifyGraph/W=$gName lowTrip=0.001
	if (strlen(AxisInfo(gName, "top"))<1)
		ModifyGraph/W=$gName/Z mirror(bottom)=1
	endif  
  	if (strlen(AxisInfo(gName, "right"))<1)
		ModifyGraph/W=$gName/Z mirror(left)=1
	endif

	// find the top wave on graph, and try to label the axes
	String wList = ImageNameList(gName,";")		// first assum an image
	if (strlen(wList))
		Wave w=ImageNameToWaveRef(gName,StringFromList(0,wList) )
	else														// next try a line plot
		wList = TraceNameList(gName,";",1)
		Wave w=TraceNameToWaveRef(gName,StringFromList(0,wList))
	endif
	if (WaveExists(w))
		String wnote=note(w)
		String left = StringByKey("GraphAxisLabelVert", wnote ,"=" )
		String bot = StringByKey("GraphAxisLabelHoriz", wnote ,"=" )
		String infoStr=ImageInfo(gName,NameOfWave(w),0)
		if (strlen(infoStr)<1)
			infoStr = TraceInfo(gName,NameOfWave(w),0)
		endif
		String leftName = StringByKey("YAXIS",infoStr)
		String botName = StringByKey("XAXIS",infoStr)
		Variable setHoriz = strlen(AxisLabelFromGraph(gName,w,botName))==0	// do not re-set if alreay labeled
		Variable setVert = strlen(AxisLabelFromGraph(gName,w,leftName))==0
		if (strlen(left) && setVert && strlen(leftName))
			Label $leftName left
		endif
		if (strlen(bot) && setHoriz && strlen(botName))
			Label $botName bot
		endif
	endif
End

//  ============================ End of Generic Graph Style ============================  //
//  ====================================================================================  //



//  ====================================================================================  //
//  =============================== Start of LaueGo Init ===============================  //

Static Function InitLaueGoFirst()
	if (JZTalwaysFirst_Version<JZTalwaysFirst_Version_Min)
		String str
		sprintf str,"\"always first.ipf\" is only version %g, but we require version>=%g\rTime to update",JZTalwaysFirst_Version,JZTalwaysFirst_Version_Min
		DoAlert 0,str
	endif
	if(exists("NewGizmo")==4)			// Only if the Gizmo XOP is available.
		Execute/Q/Z "GizmoMenu AppendItem={JZTcpSize,\"Square Up Gizmo\", \"SquareUpGizmo(\\\"\\\")\"}"
	endif
End

//  ================================ End of LaueGo Init ================================  //
//  ====================================================================================  //
