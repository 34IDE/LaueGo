#pragma rtGlobals= 2
// Constant JZTalwaysFirst_Version=2.7
#pragma version = 2.79
#pragma ModuleName=JZTalwaysFirst
#pragma hide = 1

#if NumVarOrDefault("root:Packages:SHOW_PACKAGES_MASK_JZT",-1) & 1
#include "GeneralFirst", version>=3.40
#endif

#if NumVarOrDefault("root:Packages:SHOW_PACKAGES_MASK_JZT",-1) & (4+8+16)	// LauGo or Scattering or APS Specific
#include "LaueGoFirst", version>=2.10
#endif

#if NumVarOrDefault("root:Packages:SHOW_PACKAGES_MASK_JZT",-1) & 32
//#include "LocalPackagesFirst", version>=2.07
#include "LocalPackagesFirst", optional
#endif



Menu "Misc"
	"Change Startup PrefsJZT...", JZTalwaysFirst#SetStartUpPrefsJZT()
End


Static Function AfterFileOpenHook(refNum,file,pathName,type,creator,kind)
	Variable refNum, kind
	String file,pathName,type,creator

	CheckStartupPrefs()

	if ((kind==1) || (kind==2))		// an experiment (packed or unpacked)
		PathInfo $pathName				// expand the path name, "/Users/name/data/Copper" is better than "home'
		pathName = SelectString(V_flag,pathName+":",S_path)
		String strLG = VersionStatusHash(0,"LaueGo")
		strLG = SelectString(strlen(strLG), "", "  ["+strLG+"]")
		String strLocal = VersionStatusHash(0,"LocalPackages")
		strLocal = SelectString(strlen(strLocal), "", "  ["+strLocal+"]")
		printf "\r%s  %s  restarting this file on '%s' from '%s%s'%s%s\r\r",date(),time(),getHostName(1),pathName,file,strLG,strLocal
		ExperimentModified 0			// mark this experiment as still unmodified
	endif
End
Static Function IgorStartOrNewHook(IgorApplicationNameStr)
	String IgorApplicationNameStr
	String strLG = VersionStatusHash(0,"LaueGo")
	strLG = SelectString(strlen(strLG), "", "  ["+strLG+"]")
	String strLocal = VersionStatusHash(0,"LocalPackages")
	strLocal = SelectString(strlen(strLocal), "", "  ["+strLocal+"]")
	printf "%s  %s  starting new Untitled Igor Experiment on '%s'%s%s\r\r",date(),time(),getHostName(1),strLG,strLocal

	CheckStartupPrefs()
	ExperimentModified 0				// mark this experiment as still unmodified
	return 0
	//	GetWindow kwCmdHist wsize
	//	V_right = (V_right>900) ? V_right-100 : V_right
	//	MoveWindow /C V_left+20,V_top-40,V_right,V_bottom-2
End

Static Function/S getHostName(local)	// return unix hostname
	Variable local							// 0-> network name (bob.world.com), 1-> local name (Bob's MacBook)
	if (!stringmatch(StringByKey("OS",IgorInfo(3)),"Macintosh OS X"))
		print "Only know how to get user name from a Mac"
		return ""							// cannot get answer
	endif
	if (local)								// get the name in double quotes
		ExecuteScriptText "do shell script \"scutil --get ComputerName\""	// gets local computer name
	else									// get the network name (hostname) in double quotes
		ExecuteScriptText "do shell script \"hostname\""
	endif
	Variable i=strlen(S_value)-2
	if (i<=0)								// no name,
		DoAlert 0, "Unable to get user name, message is '"+S_value+"'"
	endif
	return S_value[1,i]					// first & last char are double-quotes, remove
End

Static Function/S VersionStatusHash(full,fldr)
	Variable full							// if true returns full hash, otherwise only first 6 characters of buf
	String fldr								// usually either "LaueGo" or "LocalPackages"

	Variable f
	String fname = SpecialDirPath("Igor Pro User Files",0,0,0) + "User Procedures:"+fldr+":VersionStatus.xml"
	Open/R/Z f as fname
	if (V_flag)								// check in 6, Open has trouble following links
		fname = ReplaceString("Igor Pro 7 User",fname,"Igor Pro 6 User")
		Open/R/Z f as fname
	endif
	if (V_flag)								// check in 6, Open has trouble following links
		fname = SpecialDirPath("Igor Application",0,0,0) + "User Procedures:"+fldr+":VersionStatus.xml"
		Open/R/Z f as fname
	endif
	if (V_flag)								// check in 6, Open has trouble following links
		fname = ReplaceString("Igor Pro 7 Folder",fname,"Igor Pro 6.3 Folder")
		Open/R/Z f as fname
	endif
	if (V_flag)
		return ""							// could not open a file
	endif
	FStatus f
	if (V_logEOF<2000)					// way too small
		Close f
		return ""
	endif

	String buf=PadString("",2000,0x20)
	FBinRead f, buf
	Close f
	Variable i0=strsearch(buf,"<gitHash>",0)
	if (i0<0)
		return ""
	endif
	i0 += strlen("<gitHash>")
	Variable i1=strsearch(buf,"</gitHash>",i0)
	if (i1<i0)
		return ""
	endif
	return SelectString(full, buf[i0,i0+5], buf[i0,i1-1])
End

//  ====================================================================================  //
//  ============================== Startup Configuration ===============================  //

Static Structure StartUpPrefsJZT		// the packages available after startup
	uchar General			// General: math and string functions
	uchar Gizmo				// Gizmo: Zoom&Translate, Clip planes, Markers, and Movies
	uchar LaueGo			// LaueGo: Lattices, X-ray, and all of the LaueGo capabilities
	uchar Scattering		// X-ray Data and other scattering related functions, includes Lattice
	uchar APSspecific		// APS Specific: MDA files, BURT files, VTI files, EPICS functions, and spec files
	uchar LocalPackages	// menus showing Local Packages
EndStructure


Static Function CheckStartupPrefs()
	if (exists("root:Packages:SHOW_PACKAGES_MASK_JZT"))
		return 0									// nothing to do
	endif
	STRUCT StartUpPrefsJZT prefs
	LoadPackagePreferences/MIS=1 "JZT","startup",0,prefs
	if (V_bytesRead<5)						// looks like valid values
		SetStartUpPrefsJZT()				// run the dialog to change the startup preferences
	else
		SetMaskFromPrefsStruct(prefs)	// read valid prefs from "~/Library:Preferences/WaveMetrics/JZT/startup"
	endif
End


Static Function ShowStartUpPrefsJZT()	// load the startup preferences and print them
	STRUCT StartUpPrefsJZT prefs
	LoadPackagePreferences/MIS=1 "JZT","startup",0,prefs
	PrintStartUpPrefsJZT(prefs)
End


Static Function SetStartUpPrefsJZT()	// run the dialog to change the startup preferences
	STRUCT StartUpPrefsJZT prefs
	LoadPackagePreferences/MIS=1 "JZT","startup",0,prefs
	if (V_bytesRead<5)
		prefs.General = 1						// Default statup settings (for a new installation)
		prefs.Gizmo = 1
		prefs.LaueGo = 1
		prefs.Scattering = 1
		prefs.APSspecific = 1
		prefs.LocalPackages = 1
	endif
	prefs.General = (prefs.Gizmo || prefs.LaueGo || prefs.Scattering || prefs.APSspecific || prefs.LocalPackages) ? 1 : prefs.General
	prefs.Scattering = prefs.LaueGo ? 1 : prefs.Scattering

	Variable General, Gizmo, LaueGo, Scattering, APSspecific, LocalPackages
	General = !(! prefs.General) + 1	// set values for pop-up menus
	Gizmo = !(! prefs.Gizmo) + 1
	LaueGo = !(! prefs.LaueGo) + 1
	Scattering = !(! prefs.Scattering) + 1
	APSspecific = !(! prefs.APSspecific) + 1
	LocalPackages = !(! prefs.LocalPackages) + 1
	Prompt General,"General purpose Functions",popup,"--;General purpose functions"
	Prompt Gizmo,"Gimzo Functions",popup,"--;Gizmo functions"
	Prompt LaueGo,"Laue Go Functions",popup,"--;LauGo system"
	Prompt Scattering,"X-ray Data Functions & Lattice",popup,"--;Scattering & X-ray data"
	Prompt APSspecific,"APS Specific functions",popup,"--;APS specific functions"
	Prompt LocalPackages,"Local Packages menu",popup,"--;Show Local Packages Menu"
	DoPrompt "Capabilities to have available",General,Gizmo,LaueGo,Scattering,APSspecific,LocalPackages
	if (V_flag)
		return -1
	endif
	prefs.General = General == 2			// set values from pop-up menus
	prefs.Gizmo = Gizmo == 2
	prefs.LaueGo = LaueGo == 2
	prefs.Scattering = Scattering == 2
	prefs.APSspecific = APSspecific == 2
	prefs.LocalPackages = LocalPackages == 2
	prefs.General = (prefs.Gizmo || prefs.LaueGo || prefs.Scattering || prefs.APSspecific || prefs.LocalPackages) ? 1 : prefs.General
	prefs.Scattering = prefs.LaueGo ? 1 : prefs.Scattering

	SavePackagePreferences/FLSH=1 "JZT","startup",0,prefs
	PrintStartUpPrefsJZT(prefs)
	Variable mask = SetMaskFromPrefsStruct(prefs)
	return mask
End

Static Function SetMaskFromPrefsStruct(prefs)
	STRUCT StartUpPrefsJZT &prefs
	Variable mask = prefs.General		// make sure that local copy of mask is correctly set.
	mask += 2*(prefs.Gizmo)
	mask += 4*(prefs.LaueGo)
	mask += 8*(prefs.Scattering)
	mask += 16*(prefs.APSspecific)
	mask += 32*(prefs.LocalPackages)
	NewDataFolder/O root:Packages
	Variable/G root:Packages:SHOW_PACKAGES_MASK_JZT = mask

	Execute "SetIgorOption poundDefine=DOESNTMATTER"		// mark all procedures as needing compile
	Execute "SetIgorOption poundUnDefine=DOESNTMATTER"	// don't leave this defined. note, DOESNTMATTER is an arbitrary string
	Execute/P/Q/Z "COMPILEPROCEDURES "				// re-compile (all)
	return mask
End

Static Function PrintStartUpPrefsJZT(prefs)
	STRUCT StartUpPrefsJZT &prefs

	String labels = "General:General, includes math and string functions.;"
	labels += "Gizmo:Gizmo, includes Zoom&Translate, Clip planes, Markers, and Movies.;"
	labels += "Scattering:X-ray Data & Lattice, includes X-ray data capabilities & Lattice;"
	labels += "LaueGo:LaueGo: includes Lattices, X-ray, and all of the LaueGo capabilities.;"
	labels += "APS:APS Specific, includes MDA file, BURT file, VTI file, EPICS, and spec file capabilities.;"
	labels += "LocalPackages:Show Menus for Local Packages.;"
	Variable someOn  = (prefs.General || prefs.Gizmo || prefs.LaueGo || prefs.Scattering || prefs.APSspecific || prefs.LocalPackages)
	Variable someOff = !(prefs.General && prefs.Gizmo && prefs.LaueGo && prefs.Scattering && prefs.APSspecific && prefs.LocalPackages)

	if (someOn)
		print " "
		print "		Capabilities Available (but not Loaded) at start are:" 
		if (prefs.General)
			print StringByKey("General",labels)
		endif
		if (prefs.Gizmo)
			print StringByKey("Gizmo",labels)
		endif
		if (prefs.LaueGo)
			print StringByKey("LaueGo",labels)
		endif
		if (prefs.Scattering)
			print StringByKey("Scattering",labels)
		endif
		if (prefs.APSspecific)
			print StringByKey("APS",labels)
		endif
		if (prefs.LocalPackages)
			print StringByKey("APS",labels)
		endif
	endif

	if (someOff)
		print " "
		print "		Capabilities NOT Available:" 
		if (! prefs.General)
			print "NO ",StringByKey("General",labels)
		endif
		if (! prefs.Gizmo)
			print "NO ",StringByKey("Gizmo",labels)
		endif
		if (! prefs.LaueGo)
			print "NO ",StringByKey("LaueGo",labels)
		endif
		if (! prefs.Scattering)
			print "NO ",StringByKey("Scattering",labels)
		endif
		if (! prefs.APSspecific)
			print "NO ",StringByKey("APS",labels)
		endif
		if (! prefs.LocalPackages)
			print "NO ",StringByKey("LocalPackages",labels)
		endif
	endif
End


// Set those items that are passed in the funciton
Static Function ChangeStartupPrefsJZT([General, Gizmo, LaueGo, Scattering, APS, LocalPackages, printIt])
	Variable General, Gizmo, LaueGo, Scattering, APS, LocalPackages
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : 0

	STRUCT StartUpPrefsJZT prefs
	LoadPackagePreferences/MIS=1 "JZT","startup",0,prefs

	// make the requested changes
	if (!ParamIsDefault(General) && numtype(General)==0)
		prefs.General = !(!General)
	endif
	if (!ParamIsDefault(Gizmo) && numtype(Gizmo)==0)
		prefs.Gizmo = !(!Gizmo)
	endif
	if (!ParamIsDefault(LaueGo) && numtype(LaueGo)==0)
		prefs.LaueGo = !(!LaueGo)
	endif
	if (!ParamIsDefault(Scattering) && numtype(Scattering)==0)
		prefs.Scattering = !(!Scattering)
	endif
	if (!ParamIsDefault(APS) && numtype(APS)==0)
		prefs.APSspecific = !(!APS)
	endif
	if (!ParamIsDefault(LocalPackages) && numtype(LocalPackages)==0)
		prefs.LocalPackages = !(!LocalPackages)
	endif

	prefs.General = (prefs.Gizmo || prefs.LaueGo || prefs.Scattering || prefs.APSspecific || prefs.LocalPackages) ? 1 : prefs.General
	prefs.Scattering = prefs.LaueGo ? 1 : prefs.Scattering

	SavePackagePreferences/FLSH=1 "JZT","startup",0,prefs
	if (printIt)
		PrintStartUpPrefsJZT(prefs)
	endif
	Variable mask = SetMaskFromPrefsStruct(prefs)
	return mask
End
