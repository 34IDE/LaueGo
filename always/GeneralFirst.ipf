#pragma rtGlobals= 2
#pragma version = 3.24
#pragma ModuleName = JZTgeneral
#pragma hide = 1
#include "Utility_JZT", version>=3.65
//	DefaultFont "Consolas"		// This is in "JonFirst.ipf", that is enough

#if (IgorVersion()<7)
Constant GIZMO_WIN_TYPE = 13			// numbers for Igor 6 and under
Constant GIZMO_WIN_BIT = 4096
Static Function IgorStartOrNewHook(IgorApplicationNameStr)
	String IgorApplicationNameStr
	if(exists("NewGizmo")==4)			// If the Gizmo XOP is available, alwalys put in this menu item.
		Execute/Q/Z "GizmoMenu AppendItem={JZTcpSize,\"Square Up Gizmo\", \"SquareUpGizmo(\\\"\\\")\"}"
	endif
	return 0
End

#else
Constant GIZMO_WIN_TYPE = 17			// numbers for Igor 7
Constant GIZMO_WIN_BIT = 65536	
Menu "Gizmo"
	"Square Up Gizmo", SquareUpGizmo("")
End
#endif


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
			"Lattices",Execute/P "INSERTINCLUDE  \"LatticeSym\", version>=3.77";Execute/P "COMPILEPROCEDURES ";Execute/P "InitLatticeSymPackage(showPanel=1)"
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
	JZTgeneral#MenuItemIfWindowTypes("Copy Window Info",1+2+64+GIZMO_WIN_BIT), /Q,  JZTgeneral#GetWindowInfo2Scrap()
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
	SubMenu "LaueGo"
		"LaueGo Version Info", /Q, CheckLaueGoVersion(1)
		"   LaueGo Version Info Deep Check...", /Q, DeepCheckActualFiles("")
		"Open LaueGo Web Page", /Q, BrowseHelpFile("http://sector33.xray.aps.anl.gov/~tischler")
		"Utility_JZT", /Q, DisplayHelpTopic/K=1/Z "JZT Utility functions in \"Utility_JZT.ipf\""
	End
End


//  ====================================================================================  //
//  ========================== Start of Help & Version Info  ===========================  //

// This is used to put up a help file (either http: or file:) using the user's web browser
Function BrowseHelpFile(urlStr)
	// search order: User/Docs/WaveMetrics/LocalPackages/doc, User/Docs/WaveMetrics/LaueGo/doc, Applications/Igor Pro/LaueGo/doc
	String urlStr

	if (!StringMatch(urlStr,"http:*") && !StringMatch(urlStr,"file:*")) // not a complete url
		// assume that urlStr is just a file name in doc's folder
		String filePath = SpecialDirPath("Igor Pro User Files",0,0,0)+"User Procedures:LocalPackages:doc:"+urlStr
		GetFileFolderInfo/Q/Z=1 filePath
		if (!V_isFile)								// not in LocalPackages/doc, try Users/LaueGo/doc
			filePath = SpecialDirPath("Igor Pro User Files",0,0,0)+"User Procedures:LaueGo:doc:"+urlStr
			GetFileFolderInfo/Q/Z=1 filePath
			if (!V_isFile)							// not in Users/LaueGo, try Apps/Igor Pro/LaueGo/doc
				filePath = SpecialDirPath("Igor Application",1,0,0)+"User Procedures:LaueGo:doc:"+urlStr
				GetFileFolderInfo/Q/Z=1 filePath
				if (!V_isFile)						// give up, cannot find file
					return 1
				endif
			endif
		endif
		filePath = ParseFilePath(5,filePath,"/",0,0)	// convert from Mac to POSIX
		urlStr = "file://"+filePath							// complete the URL
	endif
	BrowseURL/Z urlStr
	if (V_flag)
		printf "ERROR -- BrowseHelpFile, unable to open   \"%s\"\r",urlStr
	endif
End


Function CheckLaueGoVersion(alert)		// Check if this version is the most recent
	Variable alert			// also put up the alert dialog

	String thisVers = VersionStatusFromDisk("")		// the entire VersionStatus.xml
	thisVers = VSbuf2infoStr(thisVers)
	String VersionStatusPath = ParseFilePath(1,FunctionPath("JZTgeneral#VersionStatusFromDisk"),":",1,1)+"VersionStatus.xml"
	thisVers = ReplaceStringByKey("VersionStatus",thisVers,VersionStatusPath,"=")

	if (strlen(thisVers)<1)
		String str
		sprintf str, "ERROR -- Could not find VersionStatus.xml"
		printf "\r  "+str+"\r"
		DoAlert 0, str
	endif
	String dateStr, timeStr, out, thisHash, latestVers, latestHash
	Variable fileCount, thisEpoch, latestEpoch
	dateStr = StringByKey("date",thisVers,"=")
	timeStr = StringByKey("time",thisVers,"=")
	fileCount = NumberByKey("fileCount",thisVers,"=")
	thisEpoch = NumberByKey("epoch",thisVers,"=")
	thisHash = StringByKey("gitHash",thisVers,"=")

	printf "\r  This version of LaueGo was created:  %s,  %s,   %g files,  commit %s\r",dateStr,timeStr,fileCount,thisHash[0,6]
	sprintf str, "This version of LaueGo with %g files was created:\r  %s,  %s\r",fileCount,dateStr,timeStr
	out = str

	latestVers = VersionStatusFromWeb()			// the entire VersionStatus.xml
	latestVers = VSbuf2infoStr(latestVers)
	latestVers = ReplaceStringByKey("VersionStatus",latestVers,"Web","=")
	latestEpoch = NumberByKey("epoch",latestVers,"=")
	latestHash = StringByKey("gitHash",latestVers,"=")

	if (strlen(latestVers)<1)
		print "\r  Unable to find most recent version info from Web"
		out += "\r  Unable to find most recent version info from Web,\r  check your network connection."
	elseif (StringMatch(thisHash,latestHash))
			str = "  an exact match to the most recent version."
			print str
			out += "\r"+str
	elseif (abs(latestEpoch-thisEpoch) <= 1)
			str = "The hashes do not match, but the most recent version has the same date & time."
			print str
			out += "\r"+str
	else
		dateStr = StringByKey("date",latestVers,"=")
		timeStr = StringByKey("time",latestVers,"=")
		fileCount = NumberByKey("fileCount",latestVers,"=")
		printf "The most recent version was created:  %s,  %s,   %g files\r",dateStr,timeStr,fileCount
		sprintf str, "\rThe most recent version with %g files was created:\r  %s,  %s\r",fileCount,dateStr,timeStr
		out += str
	endif
	print " "
	if (alert)
		DoAlert 0,out
	endif
End


Function/T DeepCheckActualFiles(source,[printIt])
	// check each of  the actual (local) files against contents of VersionStatus.xml
	// this runs slowly since it must read and compute the hash for every file in VersionStatus.xml
	String source		// must be "Web" or "Local" (source of VersionStatus.xml)
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? 0 : !(!printIt)

	source = SelectString(StringMatch(source,"Web*"),source,"Web")
	source = SelectString(StringMatch(source,"Local*"),source,"Local")
	if (WhichListItem(source,"Web;Local;")<0)
		Prompt source "Source of VersionStatus.xml", popup, "Web;Local;"
		DoPrompt "VersionStatus", source
		if (V_flag)
			return ""
		endif
		printIt = 1
	endif

	String buf=""
	if (StringMatch(source,"Web*"))
		buf = VersionStatusFromWeb()
	elseif (StringMatch(source,"Local*"))
		buf = VersionStatusFromDisk("")
	else
		return ""
	endif
	String mismatch = CheckVersionStatusAgainstDisk(buf,progress=1)

	if (printIt)
		printf "DeepCheckActualFiles(\"%s\"", source
		if (!ParamIsDefault(printIt))
			printf ", printIt=%g",printIt
		endif
		printf (")\r")
		if (strlen(mismatch))
			printf "mismatched files are:  %s\r",mismatch
		else
			print "All files match."
		endif
	endif
	return mismatch
End
//
Static Function/T CheckVersionStatusAgainstDisk(buf,[progress])
	String buf
	Variable progress
	progress = ParamIsDefault(progress) || numtype(progress) ? 1 : progress

	Variable fileCount = str2num(XMLtagContents("fileCount",buf))
	progress = numtype(fileCount) || fileCount<1 ? 0 : progress	// no progress bar without valid fileCount
	String fldr = ParseFilePath(1,FunctionPath("JZTgeneral#VersionStatusFromDisk"),":",1,1)
	String fileName="x", list, mismatch=""
	Variable i
	if (progress)
		String progressWin = ProgressPanelStart("",stop=1,showTime=1,status="Checking "+num2str(fileCount)+" files.")
	endif
	for (i=0; strlen(fileName);i+=1)
		if (progress && mod(i,5)==0)
			if (ProgressPanelUpdate(progressWin,i/fileCount*100))
				break
			endif
		endif
		fileName = XMLtagContents("file",buf,occurance=i)
		if (StringMatch(fileName,"*.app"))		// cannot check .apps (they are folders)
			continue
		endif

		list = XMLattibutes2KeyList("file",buf,occurance=i)
		if ( strlen(fileName) && CheckOneFile(list,fldr + ReplaceString("/",fileName,":")) )
			mismatch += fileName+";"
		endif
	endfor
	if (progress)
		printf "total execution time = %s\r",Secs2Time(SecondsInProgressPanel(progressWin),5,0)
		DoWindow/K $progressWin
	endif
	return mismatch
End
//
Static Function CheckOneFile(list,fileName)		// returns 1 on mismatch
	String list			// list of info from VersionStatus
	String fileName	// full file path to the file

	Variable f			// read in full contents of fileName
	Open/R/Z=1 f as fileName
	if (V_flag)
		return 1
	endif
	FStatus f
	String buf = PadString("",V_logEOF,0x20)
	FBinRead f, buf		// read in contents of fileName
	Close f
	if (StringMatch(StringByKey("sha256",list,"="),Hash(buf,1)))
		return 0				// sha256 hashes match, return 0
	endif
	return 1					// hashes do not match, return 1
End



//Static Function/T LaueGoLatestVersionInfo()// get latest verion info from web site
//	String VersionStatusURL = "http://sector33.xray.aps.anl.gov/~tischler/igor/VersionStatus.xml"
//	String vs = FetchURL(VersionStatusURL)
//	if (GetRTError(0))
//		print "***",GetRTErrMessage(), "  Check your network connection."
//		Variable err = GetRTError(1)		// clears the error
//	endif
//	if (!(strlen(vs)>0))
//		return ""
//	elseif (strsearch(vs,"404",0)>0 && strsearch(vs,"Not Found",0,2)>0)
//		return ""
//	endif
//	vs = VSbuf2infoStr(vs)
//	vs = ReplaceStringByKey("VersionStatus",vs,"Web","=")
//	return vs
//End
//
//Static Function/T ThisLaueGoVersion()			// returns info from the Users VersionStatus.xml file
//	String VersionStatusPath = ParseFilePath(1,FunctionPath("JZTgeneral#ThisLaueGoVersion"),":",1,1)+"VersionStatus.xml"
//	Variable f
//	Open/R/Z=1 f as VersionStatusPath			// read top of VersionStatus.xml
//	if (V_flag)
//		return ""
//	endif
//	FStatus f
//	String buf=PadString("",V_logEOF,0x20)	// buf set to full size of file
//	FBinRead f, buf
//	Close f
//	String out = VSbuf2infoStr(buf)
//	out = ReplaceStringByKey("VersionStatus",out,VersionStatusPath,"=")
//	return out
//End
//
Static Function/T VSbuf2infoStr(buf)	// convert VersionStatus.xml to a key=value string
	String buf									// input xml

	buf = XMLtagContents("VersionStatus",buf)
	if (strlen(buf)<1)
		return ""
	endif
	buf = XMLremoveComments(buf)

	String written = XMLattibutes2KeyList("written",buf)
	String isoStr = StringByKey("isoTime",written,"=")
	Variable epoch = ISOtime2IgorEpoch(isoStr)
	String sourceFolder = XMLtagContents("sourceFolder",buf)
	String gitHash = XMLtagContents("gitHash",buf)

	Variable fileCount = str2num(XMLtagContents("fileCount",buf))

	String out=""
	out = ReplaceStringByKey("date",out,StringByKey("date",written,"="),"=")
	out = ReplaceStringByKey("time",out,StringByKey("time",written,"="),"=")
	out = ReplaceStringByKey("isoTime",out,isoStr,"=")
	out = ReplaceNumberByKey("fileCount",out,fileCount,"=")
	if (numtype(epoch)==0 && epoch>0)
		out = ReplaceNumberByKey("epoch",out,epoch,"=")
	endif
	if (strlen(gitHash))
		out = ReplaceStringByKey("gitHash",out,gitHash,"=")
	endif
	return out
End
//Static Function/T VS2infoStr(buf)		// convert top part of VersionStatus.xml to a key=value string
//	String buf									// input xml
//
//	// extract date & time info from buf
//	Variable i=strsearch(buf,"<written ",0,2)
//	if (i<0)
//		return ""
//	endif
//	buf = buf[i+9,Inf]
//	String dateStr = getDelimitedString(buf[strsearch(buf,"date=\"",0,2),Inf])
//	String timeStr = getDelimitedString(buf[strsearch(buf,"time=\"",0,2),Inf])
//	String isoStr = getDelimitedString(buf[strsearch(buf,"isoTime=\"",0,2),Inf])
//	Variable epoch = ISOtime2IgorEpoch(isoStr)
//
//	// extract fileCount from buf
//	Variable fileCount=str2num(buf[strsearch(buf,"<fileCount>",0,2)+11,Inf])
//
//	i = strsearch(buf,"<gitHash>",0,2)	// extract gitHash from top of xml
//	String gitHash = buf[i+9,i+200]
//	i = strsearch(gitHash,"<",0,2)
//	gitHash = gitHash[0,i-1]
//
//	String out=""
//	out = ReplaceStringByKey("date",out,dateStr,"=")
//	out = ReplaceStringByKey("time",out,timeStr,"=")
//	out = ReplaceStringByKey("isoTime",out,isoStr,"=")
//	out = ReplaceNumberByKey("fileCount",out,fileCount,"=")
//	if (numtype(epoch)==0 && epoch>0)
//		out = ReplaceNumberByKey("epoch",out,epoch,"=")
//	endif
//	if (strlen(gitHash))
//		out = ReplaceStringByKey("gitHash",out,gitHash,"=")
//	endif
//	return out
//End
//
//Static Function/T getDelimitedString(buf,[delim])		// returns first occurance of a string delimited by delim in buf
//	//	so, getDelimitedString("date=\"Monday, March 23, 2015\" "), returns "Monday, March 23, 2015"
//	//	or, getDelimitedString("date=_Monday, March 23, 2015_ ",delim="_"), returns "Monday, March 23, 2015"
//	// or, getDelimitedString("date=_Monday, March 23, 2015_ "), returns ""
//	String buf			// input string
//	String delim		// delimiter, defaults to double-quote
//	delim = SelectString(ParamIsDefault(delim),delim,"\"")
//	delim = SelectString(strlen(delim),"\"",delim)
//
//	Variable i0,i1
//	i0 = strsearch(buf,delim,0,2)
//	i1 = strsearch(buf,delim,i0+1,2)
//	if (i0<0 || i1<=i0 || strlen(buf)<=2)
//		return ""
//	endif
//	return buf[i0+1,i1-1]
//End


Static Function/T VersionStatusFromWeb()
	// return the full body of VersionStatus.xml from the web site
	String VersionStatusURL = "http://sector33.xray.aps.anl.gov/~tischler/igor/VersionStatus.xml"
	String vs = FetchURL(VersionStatusURL)
	if (GetRTError(0))
		print "***",GetRTErrMessage(), "  Check your network connection."
		Variable err = GetRTError(1)		// clears the error
	endif
	if (!(strlen(vs)>0))
		return ""
	elseif (strsearch(vs,"404",0)>0 && strsearch(vs,"Not Found",0,2)>0)
		return ""
	endif
	return vs
End
//
Static Function/T VersionStatusFromDisk(VersionStatusPath)
	// returns info from the User Procedures VersionStatus.xml file
	String VersionStatusPath			// full path to VersionStatus.xml (includes the file name)
	if (strlen(VersionStatusPath)<1)
		VersionStatusPath = ParseFilePath(1,FunctionPath("JZTgeneral#VersionStatusFromDisk"),":",1,1)+"VersionStatus.xml"
	endif

	Variable f
	Open/R/Z=1 f as VersionStatusPath			// read all of VersionStatus.xml
	if (V_flag)
		return ""
	endif
	FStatus f
	String buf=PadString("",V_logEOF,0x20)	// buf set to full size of file
	FBinRead f, buf
	Close f
	return buf
End


//  =========================== End of Help & Version Info  ============================  //
//  ====================================================================================  //



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
	String gName = StringFromList(0,WinList("*",";","WIN:"+num2istr(1+2+64+GIZMO_WIN_BIT)))
	Variable topType = WinType(gName)
	topType = topType==7 ? 64 : topType		// Panel should be 64
	topType = topType==GIZMO_WIN_TYPE ? GIZMO_WIN_BIT : topType	// Gizmo window should be GIZMO_WIN_BIT

	if (flag == topType)
		return item
	endif
	return "("+item					// top window does not match kind of information in clipboard
End


Static Function/T GetWindowInfo2Scrap()
	Variable flags = 1+2+64+ (exists("NewGizmo")==4 ? GIZMO_WIN_BIT : 0)
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
	elseif (WhichListItem(wName, WinList("*",";","WIN:"+num2istr(GIZMO_WIN_BIT)))>=0 && flags>=GIZMO_WIN_BIT)	// a Gizmo/XOP
		type = "Gizmo"
		flag = GIZMO_WIN_BIT
#if (IgorVersion()<7)
		Execute "GetGizmo gizmoNameList"
		String gList=StrVarOrDefault("S_GizmoNames","")
		KillStrings/Z S_GizmoNames
#else
		GetGizmo gizmoNameList
		String gList=S_GizmoNames
#endif
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
#if (IgorVersion()<7)
		Execute "GetGizmo winPixels"			// get window position & size
#else
		GetGizmo winPixels						// get window position & size
#endif
		Variable left=NumVarOrDefault("V_left",NaN), right=NumVarOrDefault("V_right",NaN)
		Variable top=NumVarOrDefault("V_top",NaN), bottom=NumVarOrDefault("V_bottom",NaN)
		KillVariables/Z V_left, V_right, V_top, V_bottom
		if (numtype(top + bottom + left + right))
			DoAlert 0, "Unable to get Size of Gizmo Window"
			return ""
		endif
		keys = ReplaceNumberByKey("top",keys,top)
		keys = ReplaceNumberByKey("left",keys,left)
		keys = ReplaceNumberByKey("right",keys,right)
		keys = ReplaceNumberByKey("bottom",keys,bottom)

#if (IgorVersion()<7)
		Execute "GetGizmo curRotation"
		Variable euler0=NumVarOrDefault("GizmoEulerA",NaN), euler1=NumVarOrDefault("GizmoEulerB",NaN), euler2=NumVarOrDefault("GizmoEulerC",NaN)
		KillVariables/Z GizmoEulerA,GizmoEulerB,GizmoEulerC
#else
		GetGizmo curRotation
		Variable euler0=GizmoEulerA, euler1=GizmoEulerB, euler2=GizmoEulerC
#endif
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
			flag = GIZMO_WIN_BIT
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
#if (IgorVersion()<7)
		Execute "GetGizmo gizmoName"
		gName=StrVarOrDefault("S_GizmoName","")
		KillStrings/Z S_GizmoName
#else
		GetGizmo gizmoName
		gName = S_GizmoName
#endif
		if (strlen(gName)<1)
			DoAlert 0,"No Gizmo Visible to Re-Size"
			return 1
		endif
#if (IgorVersion()<7)
		Execute "GetGizmo winPixels"			// get current window position
		left = NumVarOrDefault("V_left",NaN)
		top = NumVarOrDefault("V_top",NaN)
		KillVariables/Z V_left, V_right, V_top, V_bottom
#else
		Variable V_left,V_right,V_top,V_bottom
		GetGizmo winPixels						// get current window position
		left = V_left
		top = V_top
#endif
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
#if (IgorVersion()<7)
	Execute "GetGizmo gizmoName"
	String gName=StrVarOrDefault("S_GizmoName","")
	KillStrings/Z S_GizmoName
#else
	GetGizmo gizmoName
	String gName=S_GizmoName
#endif
	if (strlen(gName)<1)
		DoAlert 0,"No Gizmo Visible to Rotate"
		return 1
	endif
#if (IgorVersion()<7)
	Execute "ModifyGizmo/N="+gName+" SETQUATERNION="+quaternion
#else
	Wave qw = str2vec(quaternion)
	ModifyGizmo/N=$gName SETQUATERNION={qw[0], qw[1], qw[2], qw[3]}
#endif
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
