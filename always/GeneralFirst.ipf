#pragma rtGlobals= 2
#pragma version = 3.15
#pragma ModuleName = JZTgeneral
#pragma hide = 1
#include "Utility_JZT", version>=3.51
//	DefaultFont "Consolas"		// This is in "JonFirst.ipf", that is enough

Static Function IgorStartOrNewHook(IgorApplicationNameStr)
	String IgorApplicationNameStr
	if(exists("NewGizmo")==4)			// If the Gizmo XOP is available, alwalys put in this menu item.
		Execute/Q/Z "GizmoMenu AppendItem={JZTcpSize,\"Square Up Gizmo\", \"SquareUpGizmo(\\\"\\\")\"}"
	endif
	return 0
End


Menu "Analysis"
	Submenu "Packages"
		"-"
		"Physical Constants",Execute/P "INSERTINCLUDE  \"PhysicalConstants\"";Execute/P "COMPILEPROCEDURES "
		help = {"Load Package of Physical Constants."}
		"FWHM",Execute/P "INSERTINCLUDE  \"FWHM\"";Execute/P "COMPILEPROCEDURES "
		help = {"Load commands for getting FWHM."}
		"Image Display Range", Execute/P "INSERTINCLUDE \"ImageDisplayScaling\"" ; 	Execute/P "COMPILEPROCEDURES "
		help = {"Load procedures that set range of color table based on a Marquee region"}
		Submenu "Gizmo"
			"Gizmo Zoom & Translate", Execute/P "INSERTINCLUDE \"GizmoZoomTranslate\", version>=2.00" ; 	Execute/P "COMPILEPROCEDURES ";Execute/P "GZoomTrans#InitGizmoZoomTranslate()"
			help = {"Provide support of Zooming and translating Gizmos"}
			"Gizmo Clip Planes", Execute/P "INSERTINCLUDE \"GizmoClip\", version>=2.00" ; 	Execute/P "COMPILEPROCEDURES ";Execute/P "GClipPlanes#InitGizmoClipPlanes()"
			help = {"Provide support of adding clip planes to a Gizmo"}
			"Gizmo Markers", Execute/P "INSERTINCLUDE \"GizmoMarkers\", version>=2.03" ; 	Execute/P "COMPILEPROCEDURES ";Execute/P "GMarkers#InitGizmoMarkers()"
			help = {"Provide support of putting 'Cursors' on a Gizmo"}
			"Gizmo Movies", Execute/P "INSERTINCLUDE \"GizmoMovies\", version>=2.02" ; 	Execute/P "COMPILEPROCEDURES ";Execute/P "GizmoMovies#InitGizmoMovies()"
			help = {"Provide support of making movies of a Gizmo"}
			"  Load All of the above Gizmo Utiliites",JZTgeneral#LoadAllGizmoUtilities()
		End

		SubMenu "Scattering..."
			"X-ray Data",Execute/P "INSERTINCLUDE  \"Xray\", version>=2.21";Execute/P "COMPILEPROCEDURES ";Execute/P "ElementDataInitPackage()"
			help = {"Load procedures for providing X-ray data"}
			"  Elements Data",Execute/P "INSERTINCLUDE  \"Elements\", version>=1.69";Execute/P "COMPILEPROCEDURES ";Execute/P "ElementDataInitPackage()"
			help = {"Load procedures for providing data on the elements"}
			"  Cromer-Liberman",Execute/P "INSERTINCLUDE  \"CromerLiberman\", version>=1.7";Execute/P "COMPILEPROCEDURES "
			help = {"Load procedures for providing Cromer-Liberman"}
			"Ion Chamber",Execute/P "INSERTINCLUDE  \"IonChamber\", version>=3.2";Execute/P "COMPILEPROCEDURES ";Execute/P "ionChamberInitPackage()"
			help = {"Load procedures for evaluating ion chamber output"}
			"Lattices",Execute/P "INSERTINCLUDE  \"LatticeSym\", version>=3.77";Execute/P "COMPILEPROCEDURES ";Execute/P "InitLatticeSymPackage()"
			help = {"Load lattice symmetry procedures"}
			"-"
			"Neutron Scattering Data",Execute/P "INSERTINCLUDE  \"Neutron\", version>=1.1";Execute/P "COMPILEPROCEDURES ";Execute/P "NeutronDataInitPackage()"
			help = {"Load procedures for providing Neutron Scattering Data"}
		End
		"-"
	End
End


Menu "File"
	"-"
	SubMenu "Add Local User Package -->"
		JZTgeneral_CheckForUserPackages(""), /Q,JZTgeneral#JZTgeneral_StartUpLocalPackage()
	End
End


Menu "Graph"
	"Set Aspect Ratio to Get Square Pixels \ Range",SetAspectToSquarePixels("")
End


Menu "Edit"
	JZTgeneral#MenuItemIfWindowTypes("Copy Window Info",1+2+64+4096), /Q,  JZTgeneral#GetWindowInfo2Scrap()
	SubMenu JZTgeneral#MenuItemIfScrapValidWindowInfo("  Paste Window","type")
		JZTgeneral#MenuItemIfScrapValidWindowInfo("Paste Window Size","left"), /Q, JZTgeneral#PutSizeWindow("")
		JZTgeneral#MenuItemIfScrapValidWindowInfo("Paste Graph Axis Ranges","axis_left_min"), /Q, JZTgeneral#PutScrapGraphAxis("")
		JZTgeneral#MenuItemIfScrapValidWindowInfo("  Paste Both Window Size & Axis Ranges","left,axis_left_min"), /Q, JZTgeneral#PutSizeWindow("");JZTgeneral#PutScrapGraphAxis("")
		"-"
		JZTgeneral#MenuItemIfScrapValidWindowInfo("Paste Gizmo Quaternion","quaternion"), /Q, JZTgeneral#PutGizmoQuaternion("")
	End
End


Menu "Misc"
	"StopAllTimers",StopAllTimers()
End


Menu "Analysis"
	Submenu "Packages"
		"(APS specific"
		Submenu "spec"
			"spec with images",Execute/P "INSERTINCLUDE  \"specImages\", version>=0.44";Execute/P "COMPILEPROCEDURES ";Execute/P "init_specImage(\"\")"
			help = {"Support for spec files plus support for spec with image files (like a Pilatus)."}
			"spec files",Execute/P "INSERTINCLUDE  \"spec\", version>=2.27";Execute/P "COMPILEPROCEDURES ";Execute/P "specInitPackage()"
			help = {"Load everything for reading spec files."}
		End
		"MDA files",Execute/P "INSERTINCLUDE  \"mdaFiles\"";Execute/P "COMPILEPROCEDURES "
		help = {"Load mda files (from APS)."}
		SubMenu "EPICS"
			"EPICS("
			"PV I\O", Execute/P "INSERTINCLUDE \"epics\"" ; 	Execute/P "COMPILEPROCEDURES ";Execute/P "epicsInitPackage()"
			help = {"Load procedures for talking to PVs via EPICS (only useful at APS)"}
			"Load Scan Record", Execute/P "INSERTINCLUDE \"LoadEPICSscans\"" ; 	Execute/P "COMPILEPROCEDURES "
			help = {"Load procedures for Loading the dump of a scan record using TkGrab, does not need EPICS support"}
		End
		"BURT Files",Execute/P "INSERTINCLUDE  \"BurtFiles\"";Execute/P "COMPILEPROCEDURES "
		help = {"Load Burt files (from APS)."}
		"WinView Reader", Execute/P "INSERTINCLUDE \"WinView\", version>=2.01"
		help = {"Load procedures reading and looking at WinView images"}
	End
End

Menu "Data"
	SubMenu "Packages"
		"MDA files from APS",Execute/P "INSERTINCLUDE  \"mdaFiles\"";Execute/P "COMPILEPROCEDURES "
		help = {"Load mda files (from APS)."}
		"EPICS Load Scan Record", Execute/P "INSERTINCLUDE \"LoadEPICSscans\"" ; 	Execute/P "COMPILEPROCEDURES "
		help = {"Load procedures for Loading the dump of a scan record using TkGrab, does not need EPICS support"}
		"BURT files from APS",Execute/P "INSERTINCLUDE  \"BurtFiles\"";Execute/P "COMPILEPROCEDURES "
		help = {"Load Burt files (from APS)."}
	End
End

Menu "Load Waves"
	SubMenu "Packages"
		"MDA files from APS",Execute/P "INSERTINCLUDE  \"mdaFiles\"";Execute/P "COMPILEPROCEDURES "
		help = {"Load mda files (from APS)."}
		"BURT files from APS",Execute/P "INSERTINCLUDE  \"BurtFiles\"";Execute/P "COMPILEPROCEDURES "
		help = {"Load Burt files (from APS)."}
		"EPICS Load Scan Record", Execute/P "INSERTINCLUDE \"LoadEPICSscans\"" ; 	Execute/P "COMPILEPROCEDURES "
		help = {"Load procedures for Loading the dump of a scan record using TkGrab, does not need EPICS support"}
	End
End


Static Function LoadAllGizmoUtilities()
	Execute/P "INSERTINCLUDE \"GizmoZoomTranslate\", version>=2.00"
	Execute/P "COMPILEPROCEDURES "
	Execute/P "GZoomTrans#InitGizmoZoomTranslate()"

	Execute/P "INSERTINCLUDE \"GizmoClip\", version>=2.00"
	Execute/P "COMPILEPROCEDURES "
	Execute/P "GClipPlanes#InitGizmoClipPlanes()"

	Execute/P "INSERTINCLUDE \"GizmoMarkers\", version>=2.00"
	Execute/P "COMPILEPROCEDURES "
	Execute/P "GMarkers#InitGizmoMarkers()"

	Execute/P "INSERTINCLUDE \"GizmoMovies\", version>=1.00"
	Execute/P "COMPILEPROCEDURES "
	Execute/P "GizmoMovies#InitGizmoMovies()"
End


Menu "Help"
	"-"
	"LaueGo Version Info", /Q, print LaueGoVersion(1)
	"Utility_JZT", /Q, DisplayHelpTopic/K=1/Z "JZT Utility functions in \"Utility_JZT.ipf\""
End

Function/T LaueGoVersion(nice)		// returns info from the VersionStatus.xml file
	Variable nice

	String VersionStatusPath = ParseFilePath(1,FunctionPath("JZTgeneral#LaueGoVersion"),":",1,1)+"VersionStatus.xml"
	GetFileFolderInfo/Q/Z=1 VersionStatusPath
	if (V_flag)
		return ""
	endif

	String buf=PadString("",1000,0x20)
	Variable f
	Open/R/Z=1 f as VersionStatusPath	// read top of VersionStatus.xml
	if (V_flag)
		return ""
	endif
	FBinRead f, buf
	Close f

	// extract date & time info from buf
	Variable i=strsearch(buf,"<written ",0)
	if (i<0)
		return ""
	endif
	buf = buf[i+9,Inf]
	i = strsearch(buf,"></written>",0)

	String dateStr="", timeStr="", isoStr="", out=""
	i = strsearch(buf,"date=\"",0)
	dateStr = buf[i+6,Inf]
	i = strsearch(dateStr,"\"",0)
	dateStr = dateStr[0,i-1]

	i = strsearch(buf,"time=\"",0)
	timeStr = buf[i+6,Inf]
	i = strsearch(timeStr,"\"",0)
	timeStr = timeStr[0,i-1]

	i = strsearch(buf,"isoTime=\"",0)
	isoStr = buf[i+9,Inf]
	i = strsearch(isoStr,"\"",0)
	isoStr = isoStr[0,i-1]

	// extract fileCount from buf
	i = strsearch(buf,"<fileCount>",0)
	Variable fileCount=str2num(buf[i+11,i+22])

	if (nice)		// for printting to History
		sprintf out, "LaueGo -- Version Created:  %s,  %s,   %g files,  info from VersionStatus:\r    %s\r",dateStr,timeStr,fileCount,VersionStatusPath
	else				// for use by other routines
		out = ReplaceStringByKey("date",out,dateStr,"=")
		out = ReplaceStringByKey("time",out,timeStr,"=")
		out = ReplaceStringByKey("isoTime",out,isoStr,"=")
		out = ReplaceNumberByKey("fileCount",out,fileCount,"=")
		out = ReplaceStringByKey("VersionStatus",out,VersionStatusPath,"=")
	endif
	return out
End

//  ====================================================================================  //
//  ========================== Start of User Procedures Menu  ==========================  //

//	This provides for an easy way to load any package in the folder "Documents:WaveMetrics:Igor Pro 6 User Files:User Procedures/LocalPackages"
//	just put the ipf file in there and it will be available to this function for loading.
//	if the first 50 lines of the file contain a line such as   "#requiredPackages "HDF5images;microGeometryN;", then it will only load that file if the 
//	ipf files in the list are already loaded.
//	if the package contains the line "#excludeFromPackageMenu", then it will not show up in the list of packages to load (useful for subroutines)
//	Note if the ipf file contains both #requiredPackages and #excludeFromPackageMenu, then the behavior depends upon order, but that situation is stupid.
//
Function/T JZTgeneral_CheckForUserPackages(dirPath)
	String dirPath				// full path to a directory

	String rootPath = SpecialDirPath("Documents", 0, 0, 0 )+"WaveMetrics:Igor Pro 6 User Files:User Procedures:LocalPackages:"
	dirPath = SelectString(strlen(dirPath),rootPath,dirPath)
	String pre = dirPath[strlen(rootPath),Inf]
	pre += SelectString(strlen(pre),"",":")
	String pathName = UniqueName("lpath",12,0)
	NewPath/O/Q/Z $pathName ,  dirPath
	if (V_flag)
		return ""
	endif
	String fname,menuList="",required
	Variable i=0, dirLast
	do
		fname = IndexedFile($pathName,i,".ipf")
		if (strlen(fname)==0)
			break
		endif
		if (allProcsArePresent(ReplaceString(".ipf",fname,"")))
			i += 1													// skip if fname is already included
			continue
		endif
		required = getListOfRequiredProcs(pathName,fname)
		if (allProcsArePresent(required))							// check if required procs are already included
			menuList = AddListItem(pre+fname,menuList,";",Inf)
		endif
		i += 1
	while (1)

	//	check subdirectories, but exclude those labeled as 'always', 'old', or 'subroutine*'
	String dirName="", dirList="", addMenus
	for (i=0,dirName=IndexedDir($pathName,0,1); strlen(dirName); i+=1, dirName=IndexedDir($pathName,i,1))
		dirLast = strlen(dirName)-1
		if (strsearch(dirName,":.",0,2)>0)										// contains a folder starting with "."
			continue
		elseif (strsearch(dirName,"subroutine",0,2)>=0)					// contains "subroutine"
			continue
		elseif (strsearch(dirName,":old",dirLast,3)==dirLast-3)		// ends in ":old"
			continue
		elseif (strsearch(dirName," old",dirLast,3)==dirLast-3)		// ends in " old"
			continue
		elseif (strsearch(dirName,"(old",0,2)>=0)							// contains "(old"
			continue
		elseif (strsearch(dirName,":always",dirLast,3)==dirLast-6)	// ends in ":always"
			continue
		elseif (strsearch(dirName,"(no copy)",dirLast,3)>=0)			// contains "(no copy)"
			continue
		endif
		addMenus = JZTgeneral_CheckForUserPackages(dirName)	// check inside this sub-directory recursively
		menuList += SelectString(strlen(addMenus),"","-;") + addMenus
	endfor
	KillPath/Z $pathName
	return menuList
End


Static Function/T getListOfRequiredProcs(path,ipf)
	String path, ipf
	String line, list=""
	Variable f, i,i0,i1
	Open/R/P=$path/T=".ipf"/Z f as ipf
	for (i=0;i<50;i+=1)							// check only first 50 lines
		FReadLine/N=500 f, line
		if (strlen(line)==0)						// end of file
			Close f
			return list
		endif
		if (strsearch(line,"#requiredPackages",0)==0)
			i0 = strsearch(line,"\"",0)
			i1 = strsearch(line,"\"",i0+1)
			if (i0==-1 || i1==-1 || i1<=i0)
				break
			endif
			list = line[i0+1,i1-1]
			break
		endif
		if (strsearch(line,"#excludeFromPackageMenu",0)==0)
			list = Hash("NON-EXISTANT",1)		// this "package" will not exist
			break
		endif
	endfor
	Close f
	list = ReplaceString(",",list,";")				// change commas to semi-colons
	return list
End
//
Static Function JZTgeneral_StartUpLocalPackage()
	GetLastUserMenuInfo							// sets S_value, V_value, etc.
	String pkg=S_value
	if (V_flag || !stringmatch(pkg,"*.ipf" ))
		return 1
	endif
	String initFuncs=getInitFunctionsName(pkg), cmd
	pkg = ReplaceString(".ipf",pkg,"")
	Variable i=strsearch(pkg,":",inf,1)+1				// start of name part, strip off any leading path part
	if (i>=strlen(pkg))
		DoAlert 0, "ERROR:\rCould not include package '"+pkg+"'"
		return 0
	endif
	pkg = pkg[i,inf]									// trim off any path part
	print "\r  adding local package  ",pkg
	sprintf cmd, "Execute/P \"INSERTINCLUDE  \\\"%s\\\"\";Execute/P \"COMPILEPROCEDURES \"", pkg
	if (strlen(initFuncs))
		cmd += ";Execute/P \""+initFuncs+" \""
	endif
	//	print "cmd = ",cmd
	Execute cmd
	return 0
End
//
Static Function allProcsArePresent(list)
	String list
	String ipf
	Variable i,N=ItemsInList(list)
	for (i=0;i<N;i+=1)
		ipf = StringFromList(i,list)+".ipf"
		if (ItemsInList(WinList(ipf,";","WIN:128"))<1)
			return 0
		endif
	endfor
	return 1
End
//
Static Function/T getInitFunctionsName(ipf)
	String ipf										// name of name of ipf file
	String pathName = UniqueName("lpath",12,0)
	NewPath/O/Q/Z $pathName ,  SpecialDirPath("Documents", 0, 0, 0 )+"WaveMetrics:Igor Pro 6 User Files:User Procedures:LocalPackages:"
	GetFileFolderInfo/P=$pathName/Q/Z=1 ipf
	if (!V_isFile && strsearch(ipf,".ipf",0,2)<0)
		ipf += ".ipf"
	endif

	String line, initFunc=""
	Variable f, i,i0,i1
	Open/R/P=$pathName/T=".ipf"/Z f as ":"+ipf
	if (V_flag)
		DoAlert 0, "ERROR:\rCould not re-open file '"+ipf+"'"
		KillPath/Z $pathName
		return ""
	endif
	for (i=0;i<50;i+=1)
		FReadLine/N=500 f, line
		if (strlen(line)==0)						// end of file
			Close f
			return initFunc
		endif
		if (strsearch(line,"#initFunctionName",0,2)==0)
			i0 = strsearch(line,"\"",0)
			i1 = strsearch(line,"\"",i0+1)
			if (i0==-1 || i1==-1 || i1<=i0)
				break
			endif
			initFunc = line[i0+1,i1-1]
			break
		endif
	endfor
	Close f
	KillPath/Z $pathName
	return initFunc
End
//
//  =========================== End of User Procedures Menu  ===========================  //
//  ====================================================================================  //



//  ====================================================================================  //
//  ========================== Start of Window Copy/Paste Info =========================  //
//
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
Static Function/T MenuItemIfScrapValidWindowInfo(item,requiredKeys)
	String item
	String requiredKeys					// comma separated list of required keys (all must be there)
	String scrap=GetScrapText()
	String type=StringByKey("type",scrap)
	if (!stringmatch(type,"*WindowInfo"))
		return "("+item					// Clipboard does not contain right kind of text
	endif

	Variable i, Nkey=ItemsInList(requiredKeys,",")
	if (Nkey>0)
		for (i=0;i<ItemsInList(requiredKeys,",");i+=1)
			String key = StringFromList(i,requiredKeys,",")
			if (strlen(StringByKey(key,scrap))<1)
				return "("+item					// the required keys are missing
			endif
		endfor
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
//
//  ========================== End of Window Copy/Paste Info ===========================  //
//  ====================================================================================  //



//  ====================================================================================  //
//  =========================== Start of Generic Graph Style ===========================  //
//
Proc Generic_Graph_Style() : GraphStyle
	Silent 1
	GenericGraphStyle("")
EndMacro
//
Function GenericGraphStyle(gName)
	String gName														// graph name, usually "", for top graph
	if (ItemsInList(WinList("*",";","WIN:1"))<1)
		return 1
	elseif (exists("GenericGraphStyleLocal")==5)
		Execute "GenericGraphStyleLocal(\""+gName+"\")"	// in case users GenericGraphStyleLocal("") is a Macro
	else
		FUNCREF GenericGraphStyleTemplate styleFunc=$"GenericGraphStyleLocal"
		styleFunc(gName)												// defaults to GenericGraphStyleTemplate(gName)
	endif
	return 0
End
//
//	This can be overridden by making a Function named  GenericGraphStyleLocal(gName)
//
Function GenericGraphStyleTemplate(gName)	// default generic graph style function
	String gName										// graph name, usually ""
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
//
//  ============================ End of Generic Graph Style ============================  //
//  ====================================================================================  //
