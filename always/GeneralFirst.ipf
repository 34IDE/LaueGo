#pragma rtGlobals= 2
#pragma version = 3.50
#pragma ModuleName = JZTgeneral
#pragma hide = 1
#include "Utility_JZT", version>=4.64
//	DefaultFont "Consolas"		// This is in "JonFirst.ipf", that is enough

strConstant TEST_TEST="test test"

#if NumVarOrDefault("root:Packages:PackagesJZT:SHOW_PACKAGES_MASK_JZT",-1) & 2		// want to include Gizmo support
#if (IgorVersion()<7)
Static Function IgorStartOrNewHook(IgorApplicationNameStr)
	String IgorApplicationNameStr
	JZTgeneral#PeriodicCheckLaueGoVersion()
	if(exists("NewGizmo")==4)			// If the Gizmo XOP is available, alwalys put in this menu item.
		Execute/Q/Z "GizmoMenu AppendItem={JZT_Line0,\"-\", \"\"}"
		Execute/Q/Z "GizmoMenu AppendItem={JZTcpSize,\"Square Up Gizmo\", \"SquareUpGizmo(\\\"\\\")\"}"
		Execute/Q/Z "GizmoMenu AppendItem={JZTcpAsptec1,\"Set Aspect 1, True\", \"ModifyGizmo aspectRatio = 1\"}"
		Execute/Q/Z "GizmoMenu AppendItem={JZTcpAsptec0,\"Set Aspect 0, Std\", \"ModifyGizmo aspectRatio = 0\"}"
	endif
	return 0
End
#else
Menu "Gizmo"
	SubMenu "Aspect"
		"Square Up Gizmo", SquareUpGizmo("")
		"Set Aspect 1, True", ModifyGizmo aspectRatio = 1
		"Set Aspect 0, Std", ModifyGizmo aspectRatio = 0
	End
	SubMenu "Zoom"
		"Set Zoom 1", SetGizmoZoom(1)
		"Set Zoom 0.8", SetGizmoZoom(0.8)
		"Set Zoom Other...", SetGizmoZoom(NaN)
	End
End
Static Function IgorStartOrNewHook(IgorApplicationNameStr)
	String IgorApplicationNameStr
	JZTgeneral#PeriodicCheckLaueGoVersion()
	return 0
End
#endif
#endif


#if NumVarOrDefault("root:Packages:PackagesJZT:SHOW_PACKAGES_MASK_JZT",-1) & 1
Menu "Analysis"
	Submenu "Packages"
		"-"
		"Physical Constants",Execute/P "INSERTINCLUDE  \"PhysicalConstants\"";Execute/P "COMPILEPROCEDURES "
		help = {"Load Package of Physical Constants."}
		"FWHM",Execute/P "INSERTINCLUDE  \"FWHM\"";Execute/P "COMPILEPROCEDURES "
		help = {"Load commands for getting FWHM."}
		"Image Display Range", Execute/P "INSERTINCLUDE \"ImageDisplayScaling\"" ; 	Execute/P "COMPILEPROCEDURES "
		help = {"Load procedures that set range of color table based on a Marquee region"}
		"Elements Data",Execute/P "INSERTINCLUDE  \"Elements\", version>=1.69";Execute/P "COMPILEPROCEDURES ";Execute/P "ElementDataInitPackage()"
		help = {"Load procedures for providing data on the elements"}
	End
End
#endif

#if NumVarOrDefault("root:Packages:PackagesJZT:SHOW_PACKAGES_MASK_JZT",-1) & 2
#if (IgorVersion()<7)			// GizmoZoomTranslate only works for Igor 7 & 8
Menu "Analysis"
	Submenu "Packages"
		Submenu "Gizmo"
			"Gizmo Zoom & Translate", Execute/P "INSERTINCLUDE \"GizmoZoomTranslate\", version>=2.00" ; 	Execute/P "COMPILEPROCEDURES ";Execute/P "GZoomTrans#InitGizmoZoomTranslate()"
			help = {"Provide support of Zooming and translating Gizmos"}
		End
	End
End
#endif
Menu "Analysis"
	Submenu "Packages"
		Submenu "Gizmo"
			"Gizmo Clip Planes", Execute/P "INSERTINCLUDE \"GizmoClip\", version>=2.00" ; 	Execute/P "COMPILEPROCEDURES ";Execute/P "GClipPlanes#InitGizmoClipPlanes()"
			help = {"Provide support of adding clip planes to a Gizmo"}
			"Gizmo Markers", Execute/P "INSERTINCLUDE \"GizmoMarkers\", version>=2.03" ; 	Execute/P "COMPILEPROCEDURES ";Execute/P "GMarkers#InitGizmoMarkers()"
			help = {"Provide support of putting 'Cursors' on a Gizmo"}
			"Gizmo Movies", Execute/P "INSERTINCLUDE \"GizmoMovies\", version>=2.02" ; 	Execute/P "COMPILEPROCEDURES ";Execute/P "GizmoMovies#InitGizmoMovies()"
			help = {"Provide support of making movies of a Gizmo"}
			"  Load All of the above Gizmo Utiliites",JZTgeneral#LoadAllGizmoUtilities()
		End
	End
End
#endif


#if NumVarOrDefault("root:Packages:PackagesJZT:SHOW_PACKAGES_MASK_JZT",-1) & 32
Menu "File"
	"-"
	SubMenu "Add Local User Package -->"
		JZTgeneral#GetUserPackagesMenuStr(), /Q, JZTgeneral#StartUpLocalPackage()
	End
End
#endif


#if NumVarOrDefault("root:Packages:PackagesJZT:SHOW_PACKAGES_MASK_JZT",-1) & 1
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
#endif


#if NumVarOrDefault("root:Packages:PackagesJZT:SHOW_PACKAGES_MASK_JZT",-1) & 16
Menu "Data"
	SubMenu "Packages"
		"MDA files from APS",Execute/P "INSERTINCLUDE  \"mdaFiles\"";Execute/P "COMPILEPROCEDURES "
		help = {"Load mda files (from APS)."}
		"VTI files from APS",Execute/P "INSERTINCLUDE  \"VTIfiles\"";Execute/P "COMPILEPROCEDURES "
		help = {"Load vti files (from APS)."}
		"EPICS Load Scan Record", Execute/P "INSERTINCLUDE \"LoadEPICSscans\"" ; 	Execute/P "COMPILEPROCEDURES "
		help = {"Load procedures for Loading the dump of a scan record using TkGrab, does not need EPICS support"}
		"BURT files from APS",Execute/P "INSERTINCLUDE  \"BurtFiles\"";Execute/P "COMPILEPROCEDURES "
		help = {"Load Burt files (from APS)."}
		"xye files from APS",Execute/P "INSERTINCLUDE  \"xyeFiles\"";Execute/P "COMPILEPROCEDURES "
		help = {"Load xye files (from APS)."}
	End
End


Menu "Load Waves"
	SubMenu "Packages"
		"MDA files from APS",Execute/P "INSERTINCLUDE  \"mdaFiles\"";Execute/P "COMPILEPROCEDURES "
		help = {"Load mda files (from APS)."}
		"VTI files from APS",Execute/P "INSERTINCLUDE  \"VTIfiles\"";Execute/P "COMPILEPROCEDURES "
		help = {"Load vti files (from APS)."}
		"BURT files from APS",Execute/P "INSERTINCLUDE  \"BurtFiles\"";Execute/P "COMPILEPROCEDURES "
		help = {"Load Burt files (from APS)."}
		"EPICS Load Scan Record", Execute/P "INSERTINCLUDE \"LoadEPICSscans\"" ; 	Execute/P "COMPILEPROCEDURES "
		help = {"Load procedures for Loading the dump of a scan record using TkGrab, does not need EPICS support"}
		"xye files from APS",Execute/P "INSERTINCLUDE  \"xyeFiles\"";Execute/P "COMPILEPROCEDURES "
		help = {"Load xye files (from APS)."}
	End
End
#endif


#if NumVarOrDefault("root:Packages:PackagesJZT:SHOW_PACKAGES_MASK_JZT",-1) & 2		// want to include Gizmo support
Static Function LoadAllGizmoUtilities()
#if (IgorVersion()<7)
	Execute/P "INSERTINCLUDE \"GizmoZoomTranslate\", version>=2.00"
	Execute/P "COMPILEPROCEDURES "
	Execute/P "GZoomTrans#InitGizmoZoomTranslate()"
#endif

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
#endif


#if (IgorVersion()<7)
Menu "Help"
	"-"
	SubMenu "LaueGo"
		"LaueGo Version Info", /Q, CheckLaueGoVersion(1)
		"   LaueGo Version Info Deep Check...", /Q, JZTgeneral#DeepCheckActualFiles("")
		"   Schedule Regular Update Checks...", JZTgeneral#SetPeriodicCheckLaueGoVersion(NaN,printIt=1)
		"Open LaueGo Web Page", /Q, BrowseHelpFile("http://sector33.xray.aps.anl.gov/~tischler")
		"Utility_JZT", /Q, DisplayHelpTopic/K=1/Z "JZT Utility functions in \"Utility_JZT.ipf\""
		SubMenu "Roman"
			"*CHARACTER*(Times New Roman,36)", /Q, GetLastUserMenuInfo; printf "\r  Roman('%s')  -->  d/%d = o/%o = x/%x\r", S_value,V_value,V_value,V_value
		End
		SubMenu "Symbol"
			"*CHARACTER*(Symbol,36)", /Q, GetLastUserMenuInfo; printf "\r  Symbol('%s')  -->  d/%d = o/%o = x/%x\r", S_value,V_value,V_value,V_value
		End
	End
End
#else
Menu "Help"
	"-"
	SubMenu "LaueGo"
		"LaueGo Version Info", /Q, CheckLaueGoVersion(1)
		"   LaueGo Version Info Deep Check...", /Q, JZTgeneral#DeepCheckActualFiles("")
		"   Schedule Regular Update Checks...", JZTgeneral#SetPeriodicCheckLaueGoVersion(NaN,printIt=1)
		"Open LaueGo Web Page", /Q, BrowseHelpFile("http://sector33.xray.aps.anl.gov/~tischler")
		"Utility_JZT", /Q, DisplayHelpTopic/K=1/Z "JZT Utility functions in \"Utility_JZT.ipf\""
		SubMenu "Characters"
			"*CHARACTER*(,36)", /Q, GetLastUserMenuInfo; printf "\r  Character('%s')  -->  d/%d = o/%o = x/%x\r", S_value,V_value,V_value,V_value
		End
	End
End
#endif



//  ====================================================================================  //
//  ================================== Start of Help ===================================  //
//
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
//
//  =================================== End of Help ====================================  //
//  ====================================================================================  //



//  ====================================================================================  //
//  ============================== Start of Version Info  ==============================  //
//
Static Structure LaueGoVersion2Struct
	int16		valid				// all are valid
	int16		days				// how freqently to auto-chceck for updates (use <0 to suppress checking)
	int16		match				// disk matches web (True or False)
	double	lastChecked		// epoch of last check
	STRUCT LaueGoVersionStruct disk
	STRUCT LaueGoVersionStruct web
EndStructure
//
Static Structure LaueGoVersionStruct
	int16		valid				// values are valid
	double	epoch				// the Igor epoch (seconds)
	int32		count				// number of files
	char		hash[40]			// python hash from VersionStatus.xml
EndStructure
//
Function CheckLaueGoVersion(alert)		// Check if this version is the most recent
	Variable alert								// also put up the alert dialog

	STRUCT LaueGoVersion2Struct ss
	FillLaueGoVersionStructUpdate(ss)	// get basic version info from local file and web
	if (!(ss.valid))
		String str
		sprintf str, "ERROR -- Could not find VersionStatus.xml"
		printf "\r  "+str+"\r"
		DoAlert 0, str
	endif
	String diskHash = ss.disk.hash
	String webHash = ss.web.hash

	String dateStr = Secs2Date(ss.disk.epoch,1)
	String timeStr = Secs2Time(ss.disk.epoch,1)
	printf "\r  This version of LaueGo was created:  %s,  %s,   %g files,  commit %s\r",dateStr,timeStr,ss.disk.count,diskHash[0,6]
	sprintf str, "This version of LaueGo with %g files was created:\r  %s,  %s\r",ss.disk.count,dateStr,timeStr
	String out = str

	if (strlen(webHash)<1)
		print "\r  Unable to find most recent distribution info from Web"
		out += "\r  Unable to find most recent distribution info from Web,\r  check your network connection."
	elseif (StringMatch(diskHash,webHash))
			str = "  an exact match to the most recent distribution."
			print str
			out += "\r"+str
	elseif (abs(ss.web.epoch - ss.disk.epoch) <= 1)
			str = "The hashes do not match, but the most recent distribution has the same date & time."
			print str
			out += "\r"+str
	else
		dateStr = Secs2Date(ss.web.epoch,1)
		timeStr = Secs2Time(ss.web.epoch,1)
		printf "The most recent distribution was created:  %s,  %s,   %g files\r",dateStr,timeStr,ss.web.count
		sprintf str, "\rThe most recent distribution with %g files was created:\r  %s,  %s\r",ss.web.count,dateStr,timeStr
		out += str
		Variable dt = ss.web.epoch - ss.disk.epoch
		String deltaString = ElapsedTime2Str(abs(dt))
		if (dt>0)
			printf "%s behind the current distribution.\r",deltaString
			sprintf str, "  %s behind the current distribution.",deltaString
			out += str
		else
			printf "%s **AHEAD** of the current distribution.\r",deltaString
			sprintf str, "  %s **AHEAD** of the current distribution.",deltaString
			out += str
		endif
	endif
	print " "
	if (alert)
		DoAlert 0,out
	endif
End


// set the frequency for checking for LaueGo updates.
Static Function SetPeriodicCheckLaueGoVersion(days,[printIt])
	Variable days
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt

	STRUCT LaueGoVersion2Struct ss
	if (numtype(days) || days <=0)
		LoadPackagePreferences/MIS=1 "LaueGo","VersionInfo",0,ss
		days = ss.days
		Prompt days, "Days between Auto update check (0 or -1 means never)"
		DoPrompt "Days",days
		if (V_flag)
			return 1
		endif
	endif
	days = numtype(days) || days <=0 ? -1 : days
	ss.days = days
	SavePackagePreferences/FLSH=1 "LaueGo","VersionInfo",0,ss
End


// If last check was more than 30 days ago, check for updates, notify user if new update available.
// call this when Igor restarts to see if updates are available
Static Function PeriodicCheckLaueGoVersion()
	STRUCT LaueGoVersion2Struct ss
	LoadPackagePreferences/MIS=1 "LaueGo","VersionInfo",0,ss
	Variable days = ss.days
	days = numtype(days) || days <=0 ? Inf : days
	if ((datetime - ss.lastChecked) > days*24*3600)		// last checked more than "days" ago
		FillLaueGoVersionStructUpdate(ss)
		if (ss.valid && !(ss.match))								// does not match, notify user
			CheckLaueGoVersion(1)
		endif
	endif
End
//
Static Function FillLaueGoVersionStructUpdate(ss)		// get basic version info from local file and web
	STRUCT LaueGoVersion2Struct &ss
	LoadPackagePreferences/MIS=1 "LaueGo","VersionInfo",0,ss
	Variable days = ss.days
	FillLaueGoVersionStruct(ss.disk,"disk")					// get basic version info from local LaueGo folder
	FillLaueGoVersionStruct(ss.web,"web")						// get basic version info from web site
	ss.lastChecked = datetime
	ss.match = StringMatch(ss.disk.hash,ss.web.hash)
	ss.valid = ss.disk.valid && ss.web.valid
	if (ss.valid)														// if values are valid, then update
		SavePackagePreferences/FLSH=1 "LaueGo","VersionInfo",0,ss		// update the saved structure
	endif
	return 0
End
//
Static Function FillLaueGoVersionStruct(s1,source)		// get basic version info from local file and web
	STRUCT LaueGoVersionStruct &s1
	String source														// should be "web" or "disk"

	String Version=""
	if (StringMatch(source,"disk"))
		Version = JZTgeneral#VersionStatusFromDisk("")	// the entire VersionStatus.xml from local LaueGo folder
	elseif (StringMatch(source,"web"))
		Version = JZTgeneral#VersionStatusFromWeb()		// the entire VersionStatus.xml from web site
	else
		Version = ""
	endif
	Version = JZTgeneral#VSbuf2infoStr(Version)

	Variable count = NumberByKey("fileCount",Version,"=")
	Variable epoch = NumberByKey("epoch",Version,"=")
	String hashStr = StringByKey("gitHash",Version,"=")
	if (numtype(count+epoch)==0 && count>0 && epoch>0 && strlen(hashStr)==40)
		s1.count = count
		s1.epoch = epoch
		s1.hash = hashStr
		s1.valid = 1
	else
		s1.valid = 0
	endif
	return 0
End

Static Function printLaueGoVersionStruct()
	STRUCT LaueGoVersion2Struct ss
	LoadPackagePreferences/MIS=1 "LaueGo","VersionInfo",0,ss
	Variable epoch = ss.lastChecked
	printf "at last check:  valid=%g, days=%g, match=%g,  epoch=%d,   %s   %s\r", ss.valid, ss.days, ss.match, ss.lastChecked,Secs2Time(epoch,1),Secs2Date(epoch,1)
	printf "Info at last check  (ocured:  %s ago):\r",ElapsedTime2Str(datetime - ss.lastChecked)
	epoch = ss.disk.epoch
	printf "  disk: valid=%g, count=%g, epoch=%d,   %s   %s   %s\r", ss.disk.valid, ss.disk.count, ss.disk.epoch, Secs2Time(epoch,1),Secs2Date(epoch,1), ss.disk.hash
	epoch = ss.web.epoch
	printf "  web:  valid=%g, count=%g, epoch=%d,   %s   %s   %s\r", ss.web.valid, ss.web.count, ss.web.epoch, Secs2Time(epoch,1),Secs2Date(epoch,1), ss.web.hash
End


Static Function/T DeepCheckActualFiles(source,[printIt])
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
			print "mismatched files are:"
			PrintLongStrings(mismatch,sep=";")
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

	String VS = XMLtagContents("VersionStatus",buf)
	Variable fileCount = str2num(XMLtagContents("fileCount",VS))
	progress = numtype(fileCount) || fileCount<1 ? 0 : progress	// no progress bar without valid fileCount
	String fldr = ParseFilePath(1,FunctionPath("JZTgeneral#VersionStatusFromDisk"),":",1,1)
	String fileName="x", list, mismatch=""
	Variable i, start1=0, start2=0
	if (progress)
		String progressWin = ProgressPanelStart("",stop=1,showTime=1,status="Checking "+num2str(fileCount)+" files.")
	endif
	for (i=0; strlen(fileName);i+=1)
		if (progress && mod(i,20)==0)
			if (ProgressPanelUpdate(progressWin,i/fileCount*100))
				break
			endif
		endif

		fileName = XMLtagContents("file",VS, start=start1)
		list = XMLattibutes2KeyList("file",VS, start=start2)
		start1 = max(start1,start2)
		start2 = max(start1,start2)

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

	GetFileFolderInfo/Q/Z=1 fileName
	if (V_isFolder)
		return 0			// do not check folders, just assume all OK (needed for *.app folders)
	endif

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
//
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
//
//  ====================================================================================  //
//  =============================== End of Version Info  ===============================  //



#if NumVarOrDefault("root:Packages:PackagesJZT:SHOW_PACKAGES_MASK_JZT",-1) & 32		// want to include support for User LocalPackages
//  ====================================================================================  //
//  =========================== Start of User Packages Menu  ===========================  //
//
//	This provides a self-configuring Menu for loading packages. It does NOT automatically load any *.ipf files.
//	It just provides a Menu to load the *.ipf files.
//
//	It looks for all Igor *.ipf files in the folder:
//			":Documents:WaveMetrics:Igor Pro N User Files:User Procedures/LocalPackages"
//	If you don't want to use a separate folder called "LocalPackages", then just execute the command:
//		String root:Packages:PackagesFolder_JZT = "My Packages"
//	where you change "My Packages" to whatever folder name that you want to use.
//	If you set:    String root:Packages:PackagesFolder_JZT="",
//	then every *.ipf file in your "User Procedures" folder will show up on the Menu.
//
//	Put any Igor packages in the LocalPackages folder (or a sub-folder) that you want available for loading.
//	They will automatically be available for including by going to the Igor MenuBar   "File:Add Local User Package-->"
//	It ignores all files in the following sub-folders:
//		if the folder name starts with a "."			// system folders
//		if the folder name contains the word "subroutine"	// hides subroutines
//		if the folder name ends in "always"			// always things
//		if the folder name is "old"						// old things�
//		if the folder name ends in " old"
//		if the folder name contains "(old"
//	
//	If in the first 50 lines of the ipf file there is a line such as:
//	#requiredPackages "HDF5images;microGeometryN;"
//	then this ipf file will only appear in the Menu if all of the ipf files in the list are already loaded.
//	In this example, the ipf file will only appear in the "File:Add Local User Package-->" Menu if the ipf 
//	files HDF5images and microGeometryN have both been PREVIOUSLY loaded.
//	This provides a way to suppress menu items that are not appropriate.
//	
//	If there is a line:
//	#excludeFromPackageMenu
//	then this file will not show up in the Menu. This is useful for subroutine type files.
//	
//	Also, if the ipf file contains a line like:
//	#initFunctionName "Init_abcPackage()"
//	
//	Then Init_abcPackage() will be run after the ipf file is loaded.  Actually, everything in the double quotes is 
//	passed to an Execute command. There is nothing special about the name.
//	
//	So if you include other packages in the package that you are loading, don't forget to include the packages 
//	init_func in you init_func.
//	
//	e.g. if you load abc.ipf (and abc.ipf has #include "xyz"), then Init_abcPackage() should call Init_xyzPackage(). 
//	Otherwise xyz will not be inited.  Note, the init function is not required, it is just available.


Static Function/T GetUserPackagesMenuStr()
	SVAR/Z str = root:Packages:PackagesJZT:LocalUserPackagesMenu_JZT
	if (!SVAR_Exists(str))					// does not exist
		UserPackagesMenuStr("")			// re-make LocalUserPackagesMenu_JZT
	endif
	SVAR str = root:Packages:PackagesJZT:LocalUserPackagesMenu_JZT
	if (cmpstr(str,"=refresh menu=;")==0)
		UserPackagesMenuStr("")			// user request to re-make menu string
	endif
	return str
End

Static Function/T UserPackagesMenuStr(dirPath)
	String dirPath				// OPTIONAL, full path to a directory, usually "" when called
	String str = "=refresh menu=;-;" + CheckForUserPackages(dirPath)
	NewDataFolder/O root:Packages
	NewDataFolder/O root:Packages:PackagesJZT
	String/G root:Packages:PackagesJZT:LocalUserPackagesMenu_JZT = str
	return str
End


Static Function/T CheckForUserPackages(dirPath)
	String dirPath				// full path to a directory
	String PackagesFolder = StrVarOrDefault("root:Packages:PackagesFolder_JZT","LocalPackages")
	String rootPath = SpecialDirPath("Igor Pro User Files", 0, 0, 0 )+"User Procedures:"+PackagesFolder
	dirPath = SelectString(strlen(dirPath),rootPath,dirPath)
	String pre = dirPath[strlen(rootPath),Inf]
	pre += SelectString(strlen(pre),"",":")
	String pathName = UniqueName("lpath",12,0)
	NewPath/O/Q/Z $pathName, dirPath
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
			i += 1															// skip if fname is already included
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
		elseif ( (strsearch(dirName,".app",dirLast,3)+3)==dirLast )	// skip ends with ".app"
			continue
		endif
		addMenus = CheckForUserPackages(dirName)	// check inside this sub-directory recursively
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
	for (i=0;i<50;i+=1)								// check only first 50 lines
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
	list = ReplaceString(",",list,";")			// change commas to semi-colons
	return list
End
//
Static Function StartUpLocalPackage()
	GetLastUserMenuInfo								// sets S_value, V_value, etc.
	String pkg=S_value

	if (!V_flag && cmpstr(pkg,"=refresh menu=")==0)
		UserPackagesMenuStr("")					// re-set the root:Packages:PackagesJZT:LocalUserPackagesMenu_JZT
		return 0
	endif

	if (V_flag || !stringmatch(pkg,"*.ipf" ))
		return 1
	endif
	String initFuncs=getInitFunctionsName(pkg), cmd
	pkg = ReplaceString(".ipf",pkg,"")
	Variable i=strsearch(pkg,":",inf,1)+1		// start of name part, strip off any leading path part
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
	Execute/P "JZTgeneral#UserPackagesMenuStr(\"\")"		// update menu string after loading
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
	String PackagesFolder = StrVarOrDefault("root:Packages:PackagesFolder_JZT","LocalPackages")
	NewPath/O/Q/Z $pathName, SpecialDirPath("Igor Pro User Files", 0, 0, 0 )+"User Procedures:"+PackagesFolder
	GetFileFolderInfo/P=$pathName/Q/Z=1 ipf
	if (!V_isFile && strsearch(ipf,".ipf",0,2)<0)
		ipf += ".ipf"
	endif

	String line, initFunc=""
	Variable f, i,i0,i1
	Open/R/P=$pathName/T=".ipf"/Z f as ipf
	if (V_flag)
		DoAlert 0, "ERROR:\rCould not re-open file '"+ipf+"'"
		KillPath/Z $pathName
		return ""
	endif
	for (i=0;i<50;i+=1)
		FReadLine/N=500 f, line
		if (strlen(line)==0)					// end of file
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

//  ============================ End of User Packages Menu  ============================  //
//  ====================================================================================  //
#endif



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
		Variable left=NumVarOrDefault("V_left",NaN), right=NumVarOrDefault("V_right",NaN)
		Variable top=NumVarOrDefault("V_top",NaN), bottom=NumVarOrDefault("V_bottom",NaN)
		KillVariables/Z V_left, V_right, V_top, V_bottom
#else
		GetGizmo winPixels							// get window position & size
		Variable left=V_left, right=V_right, top=V_top, bottom=V_bottom
#endif
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
