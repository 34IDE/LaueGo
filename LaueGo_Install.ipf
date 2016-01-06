#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName=LaueGoInstall
#pragma version = 0.11

Constant LaueGo_Install_Test = 0
Static strConstant gitHubArchiveURL = "http://github.com/34IDE/LaueGo/archive/master.zip"
Static strConstant JZTArchiveURL = "http://sector33.xray.aps.anl.gov/~tischler/igor/LaueGo.zip"
Static Constant minimumArchiveLen = 16777216	// 16MB (16*1024*1024), it must be larger than this.
#if (IgorVersion()<7)
	strConstant xopList = "HDF5.xop;HDF5 Help.ihf;MultiPeakFit.xop;MultiPeakFit Help.ihf;"	// list of xop's to install, this can be overidden
	strConstant IgorExtensionsFolder = "Igor Extensions"
#elif (StringMatch(StringByKey("IGORKIND",IgorInfo(0)),"pro64*"))
	strConstant xopList = "HDF5-64.xop;HDF5 Help.ihf;"	// list of xop's to install, this can be overidden
	strConstant IgorExtensionsFolder = "Igor Extensions (64-bit)"
#elif (StringMatch(StringByKey("IGORKIND",IgorInfo(0)),"pro*"))
	strConstant xopList = "HDF5.xop;HDF5 Help.ihf;"	// list of xop's to install, this can be overidden
	strConstant IgorExtensionsFolder = "Igor Extensions"
#else
	Abort "Cannot Identify Kind of Igor"
#endif


Menu "__LaueGo Install__"
	"Install or Reinstall LaueGo",LaueGo_Install()
End


Function LaueGo_Install()
	DoWindow/H/F			// bring history window to front so you can watch it.
	if (LaueGo_Install_Test)
		DoAlert 0, "This is operating in Testing mode"
	endif
	if (Append2Log("",1))
		print "Cannot open the Log File"
		return 1
	endif

	Variable isEmpty = ExperimentEmpty()
	String str = "***"+SelectString(isEmpty,"This Expermiment appears to NOT be empty","This is considered to be an EMPTY Experiment")
	Append2Log(str,0)
	String UserPath = SpecialDirPath("Igor Pro User Files",0,0,0)+"User Procedures:"

	if (!FileFolderExists(UserPath,folder=1))
		DoAlert 0,"The \"User Path\", does not exist, quit Igor and try again"
		return 1
	endif
	String LaueGoFullPath = UserPath+"LaueGo"		// official location of the LaueGo folder

	Variable existing = FileFolderExists(LaueGoFullPath,folder=1)
	String oldFolderDestination=""
	if (existing)
		sprintf oldFolderDestination, "LaueGo_%sT%s",Secs2Date(DateTime,-2),Secs2Time(DateTime,3)
		oldFolderDestination = SpecialDirPath("Desktop",0,0,0)+CleanupName(oldFolderDestination,0)
		if (FileFolderExists(oldFolderDestination,folder=1))
			sprintf str, "Old folder named \"%s\"  already exists!",oldFolderDestination
			DoAlert 0, str
			Append2Log(str,0)
			return 1
		endif
	endif

	// identify the new LaueGo folder
	Variable downLoad=1
	Prompt download, "source of New LaueGo Folder:", popup, "Download from APS web site;Download from GitHub;Locate new LaueGo folder on this computer"
	DoPrompt "NEW LaueGo", download
	String newFullPath=""							// full path to the expanded archive to install, i.e. the new LaueGo folder
	if (download==1)									// download LaueGo.zip, from APS web site, expand it, and put it on Desktop
		newFullPath = FetchLaueGoArchive(JZTArchiveURL)
	elseif (download==2)							// download LaueGo.zip, from GitHub site, expand it, and put it on Desktop
		newFullPath = FetchLaueGoArchive(gitHubArchiveURL)
	else													// ask user to find the new (unzipped) LaueGo folder, probably on the Desktop
		NewPath/M="Find NEW LaueGo folder"/O/Q/Z NewLaueGo
		if (V_flag)
			return 1
		endif
		PathInfo NewLaueGo
		newFullPath=S_path[0,strlen(S_path)-2]
	endif
	if (StringMatch(LaueGoFullPath,newFullPath))
		str = "The New LaueGo folder is the same as the existing LaueGo folder!"
		Append2Log(str,0)
		DoAlert 0,str
		return 1
	endif

	String name = ParseFilePath(0,newFullPath,":",1,0)
	if (!StringMatch(name,"LaueGo"))
		sprintf str, "The new folder is named \"%s\", it must be named \"LaueGo\".  Rename the folderor try again",name
		DoAlert 0, str
		Append2Log(str,0)
		return 1
	endif

	Append2Log("------------------------------",0)
	sprintf str, "The new LaueGo folder is \"%s\"",newFullPath
	Append2Log(str,0)
	if (existing)
		sprintf str, "Old folder will be moved to  \"%s\"", oldFolderDestination
		Append2Log(str,0)
	else
		Append2Log("There is no existing LaueGo folder to move away",0)
	endif
	Append2Log("------------------------------",0)

	DoAlert 1,"Proceed?\r(See History for what folders will be moved)"
	if (V_flag!=1)
		Append2Log("stopping",0)
		return 1
	endif

	if (existing)
		//	move LaueGo folder to somewere else
		if (!LaueGo_Install_Test)
			MoveFolder/I=0/M="old LaueGo folder"/S="Destination for old LaueGo folder"/Z=1 LaueGoFullPath as oldFolderDestination
			if (V_flag)
				str = "The old LaueGo folder was not moved"
				DoAlert 0, str
				Append2Log(str,0)
				return 1
			endif
		endif
		sprintf str, "Old folder has been moved to:  \"%s\"", ParseFilePath(0,oldFolderDestination,":",1,1)+":"+ParseFilePath(0,oldFolderDestination,":",1,0)
		Append2Log(str,0)
 	endif

	// move new LaueGo to "User Procedures"
	//	put new LaueGo folder into LaueGoFullPath
	if (!LaueGo_Install_Test)
		MoveFolder/I=0/M="new LaueGo folder"/S="Official destination for new LaueGo folder"/Z=1  newFullPath as LaueGoFullPath
		if (V_flag)
			str = "The new LaueGo folder was not moved"
			DoAlert 0, str
			Append2Log(str,0)
			return 1
		endif
	endif
	Append2Log("Succesfully moved the new LaueGo folder to 'User Procedures'",0)

	if (AddLocalPackagesFolder(UserPath))		// Create LocalPackages folder (if needed)
		return 1
	endif

	// install the alias to "always first.ipf"
	String alwaysFirstSource=LaueGoFullPath+":always:always first.ipf"
	String alwaysFirstDestination = SpecialDirPath("Igor Pro User Files",0,0,0)+"Igor Procedures:always first.ipf"
	if (AliasAlreadyThere(alwaysFirstSource,alwaysFirstDestination))
		Append2Log("The alias to \":always:always first.ipf\" is already in \"Igor Procedures\", no action taken",0)
	else
		Append2Log("  about to create the alias",0)
		sprintf str, "  alwaysFirstSource   \"%s\"",alwaysFirstSource
		Append2Log(str,0)
		sprintf str,"  always first destination   \"%s\"",alwaysFirstDestination
		Append2Log(str,0)
		CreateAliasShortcut/I=0/O/Z=1 alwaysFirstSource as alwaysFirstDestination
		if (V_flag)
			str = "CreateAliasShortcut failed to make alias to \"always first.ipf\""
			DoAlert 0,str
			Append2Log(str,0)
			return 1
		endif
	endif

	Append2Log("Create the Igor Help file Aliases",0)
	if (InstallAliasesToHelpFiles(LaueGoFullPath))
		Append2Log("ERROR -- failed to properly install Help File Aliases,  Continuing...",0)
	endif

	sprintf str, "ready to try and install the following xop's:  \"%s\"",xopList
	Append2Log(str,0)
	Variable count = InstallXOPsAsNeeded(xopList,0)
	if (count)
		sprintf str, "Installed %g xop's",count
		Append2Log(str,0)
	elseif (count==0)
		Append2Log("No xop's were installed, they were probably already there.",0)
	else
		Append2Log("ERROR installing xop's",0)
	endif
	Append2Log("finished",0)
	Append2Log("********************************************************************************",0)

	if (isEmpty)
		DoAlert 1, "Done with Installing LaueGo, you need to re-start Igor\rQuit (without saving) Now?"
		if (V_flag==1)
			Execute/P "Quit/N"
		endif
	else
		DoAlert 1, "Done with Installing LaueGo, you need to re-start Igor\rSave and then Quit Now?"
		if (V_flag==1)
			Execute/P "Quit/Y"
		endif
	endif
End



Static Function AddLocalPackagesFolder(UserPath)	// Create LocalPackages folder if needed, returns 0=success, 1=error
	String UserPath

	String FullPath = UserPath+"LocalPackages"
	if (FileFolderExists(FullPath,folder=1))			// LocalPackages already exists, return no error
		return 0
	endif

	String AboutTextBody = "About LocalPackages Folder.txt\r\r"	// text to insert in the file
	AboutTextBody += "Put any Igor packages in this folder (or a sub-folder)that you want available for loading.  They will automatically be available for including by going to the Igor MenuBar   \"File:Add Local User Package-->\"\r"
	AboutTextBody += "It ignores all files in the following sub-folders:\r"
	AboutTextBody += "	if the folder name starts with a \".\"			// system folders\r"
	AboutTextBody += "	if the folder name contains the word \"subroutine\"	// hides subroutines\r"
	AboutTextBody += "	if the folder name ends in \"always\"			// always things\r"
	AboutTextBody += "	if the folder name is \"old\"				// old thingsÉ\r"
	AboutTextBody += "	if the folder name ends in \" old\"\r"
	AboutTextBody += "	if the folder name contains \"(old\"\r\r"
	AboutTextBody += "If in the first 50 lines of the ipf file there is a line such as:\r\r"
	AboutTextBody += "#requiredPackages \"HDF5images;microGeometryN;\"\r\r"
	AboutTextBody += "then this ipf file will only appear in the menu if all of the ipf files in the list are already loaded.  In this example, the ipf file will only appear in the \"File:Add Local User Package-->\" menu if the ipf files HDF5images and microGeometryN have both been previously loaded.\r\r"
	AboutTextBody += "Also, if there is a line like:\r\r"
	AboutTextBody += "#initFunctionName \"Init_xyzPackage()\"\r\r"
	AboutTextBody += "Then the function Init_xyzPackage() will be run after the ipf file is loaded.  Actually, everything in the double quotes is passed to an Execute command.\r"

	String LPpath=UniqueName("localPackages",12,0), str
	NewPath /C/O/Q/Z $LPpath, FullPath		// create the LocalPackages folder
	KillPath/Z $LPpath
	if (V_flag)
		str = "Could not create the 'LocalPackages' folder named \""+FullPath+"\", the new LaueGo folder was not moved"
		DoAlert 0, str
		Append2Log(str,0)
		return 1
	endif
	Append2Log("Created LocalPackages folder \""+FullPath+"\"",0)

	Variable f															// add the "About LocalPackages Folder.txt" file
	Open/Z f as FullPath+":About LocalPackages Folder.txt"
	if (V_flag)
		Append2Log("Minor Error -- Failed to write \"About LocalPackages Folder.txt\" into the LocalPackages folder",0)
	else
		FBinWrite f, AboutTextBody
		Close f
		Append2Log("wrote \"About LocalPackages Folder.txt\" into the LocalPackages folder",0)
	endif
	return 0
End



Static Function/T FetchLaueGoArchive(urlStr)			// download and expand new LaueGo zip file to the Desktop
	String urlStr

	String str
	String destinationZip = SpecialDirPath("Desktop",0,0,0)+"LaueGo.zip"
	String DesktopFldr = SpecialDirPath("Desktop",0,0,0)
	String destinationFldr = DesktopFldr+"LaueGo"
	if (FileFolderExists(destinationZip,file=1))
		str = "the file 'LaueGo.zip' already exists on the desktop\rCannot overwrite existing zip file.\rRemove LaueGo.zip from Desktop and try again."
		Append2Log(str,0)
		DoAlert 0, str
		return ""
	endif
	if (FileFolderExists(destinationFldr))
		str = "the folder 'LaueGo' already exists on the Desktop\rCannot overwrite existing folder.\rTry changing name of existing folder and try again."
		Append2Log(str,0)
		DoAlert 0, str
		return ""
	endif

//	if (strlen(	urlStr)<1)
//		Prompt urlStr, "URL to use for Downloading LaueGo", popup, "APS;gitHub"
//		DoPrompt "Download Site",urlStr
//		if (V_flag)
//			return ""
//		endif
//		urlStr = SelectString(StringMatch(urlStr,"gitHub"),urlStr,gitHubArchiveURL)
//		urlStr = SelectString(StringMatch(urlStr,"APS"),urlStr,JZTArchiveURL)
//	endif
	if (strlen(urlStr)<1)
		str = "Could not get the URL to the LaueGo.zip file.\rTry again."
		Append2Log(str,0)
		DoAlert 0, str
		return ""
	endif

	sprintf str, "Downloading LaueGo from \"%s\"\r",urlStr
	Append2Log(str,0)
	if (StringMatch(urlStr,gitHubArchiveURL))
		str = "gitHub, must use a web browser to download the zip archive.\r"
		str += "On the Web Page, Push the 'Download ZIP' button (on the right)\rThen Follow instructions in Igor History"
		Append2Log(str,0)
		DoAlert 0, str
		BrowseURL/Z "http://github.com/34IDE/LaueGo/"
		if (V_flag)
			str = "Could NOT get gitHub zip archive using Web Browser"
			Append2Log(str,0)
			return ""
		endif
		str = "\r\rYou must now:\r1) find the downloaded file called \"LaueGo-master.zip\"\r2) un-zip it (double click on it)\r3) rename the folder \"LaueGo-master\" to \"LaueGo\"\r"
		str += "You will then have to re-start the install and this time do NOT choose to download, but just select this new \"LaueGo\" folder"
		Append2Log(str,0)
		return ""
	endif

	String response = FetchURL(urlStr)
	if (strlen(response)< minimumArchiveLen)
		sprintf str, "\rDownloaded zip archive is too short, length = %d\rProbably a download error.\rTry again.",strlen(response)
		Append2Log(str,0)
		DoAlert 0, str
		return ""
	endif
	Variable f
	Open/M="Save the LaueGo.zip file"/T=".zip"/Z=1f as destinationZip
	FBinWrite f, response
	Close f

	Variable err = UnZipFiles(destinationZip, DesktopFldr, deleteZip=1, overWrite=0, printIt=1)
	if (err)
		return ""
	endif

//	String cmd=""
//	sprintf str, "set zipPath to quoted form of POSIX path of (\"%s\" as alias)\n",destinationZip
//	cmd += str
//	sprintf str, "set outPath to quoted form of POSIX path of (\"%s\" as alias)\n",DesktopFldr
//	cmd += str
//	cmd += "tell application \"System Events\"\n"
//	cmd += "do shell script \"unzip \" & zipPath & \" -d \" & outPath\n"
//	cmd += "end tell\n"
//	sprintf str, "set LGfolder to \"%s\" as alias\n",destinationFldr
//	cmd += str
//	sprintf str, "	tell application \"Finder\" to set label index of LGfolder to 6"	// 6=Green
//	cmd += str
//	//	print ReplaceString("\n",cmd,"\r")
//	ExecuteScriptText/Z cmd
//	if (V_flag)
//		sprintf str, "failure in unzipping \"%s\"\rcmd = \r'%s'\r\r",destinationZip,ReplaceString("\n",cmd,"\r")
//		Append2Log(str,0)
//		Append2Log(S_value+"\r\r",0)				// there is no valid file path
//	endif
//	DeleteFile /Z destinationZip					// done unzipping, delete the zip file
//	if (V_flag)
//		sprintf str, "Could NOT delete the zip file  \"%s\"\r\tProceeding..."
//		Append2Log(str,0)
//	endif

	if (FileFolderExists(destinationFldr,folder=1))
		sprintf str, "Done Downloading, the new folder   \"%s\"",destinationFldr
		Append2Log(str,0)
	else
		sprintf str, "ERROR --Done Downloading, but NO NEW FOLDER!!   \"%s\"",destinationFldr
		Append2Log(str,0)
	endif
	return destinationFldr
End



// ============================================================================================= //
// =================================== Start of XOP Aliases ==================================== //

//		InstallXOPsAsNeeded("HDF5.xop;HDF5 Help.ihf;",1)
//		InstallXOPsAsNeeded("HDF5.xop;HDF5 Help.ihf;HFSAndPosix.xop;HFSAndPosix Help.ihf;",1)
//
Function InstallXOPsAsNeeded(list,restart)
	String list										// list of .xop's and .ihf's to put in Users "Igor Extensions"
	Variable restart								// 1=ask to restart, 0=don't ask to restart
	restart = numtype(restart) ? 1 : !(!restart)

	PathInfo Igor
	String IgorRoot=S_path
	if (!FileFolderExists(IgorRoot,folder=1))
		Abort "Cannot find the Igor path, quit Igor and try again"
	endif
	if (!FileFolderExists(SpecialDirPath("Igor Pro User Files",0,0,0)+IgorExtensionsFolder,folder=1))
		Abort "Cannot find path to the Users 'Igor Pro User Files' in Users Documents folder, quit Igor and try again"
	endif

	String name, item, str, listFull=""
	Variable i,N=ItemsInList(list)
	for (i=0;i<N;i+=1)
		name = StringFromList(i,list)
		item = FindFileInDirPath(IgorRoot,name)
		if (strlen(item))
			listFull += item+":"+name+";"
		else
			sprintf str, "ERROR -- InstallXOPsAsNeeded(), Could not find path to '%s'",name
			Append2Log(str,0)
			return NaN
		endif
	endfor

	String errStr, fullName
	Variable count=0, ic
	for (i=0;i<ItemsInList(listFull);i+=1)
		fullName = StringFromList(i,listFull)
		try
			ic = CopyAliases2IgorExtensionsLocal(fullName)
			count += numtype(ic) ? 0 : ic
		catch
			sprintf errStr, "ERROR -- in CopyAliases2IgorExtensionsLocal(), could not make alias of \"%s\"",fullName
			if (strlen(GetRTErrMessage()))
				errStr += "\r"+GetRTErrMessage()
			endif
			Variable err = GetRTError(1)
			DoAlert 0, errStr
			Append2Log(errStr,0)
			break
		endtry
	endfor

	if (count>0 && restart)
		DoAlert 1, "Installed "+num2str(count)+" xops & ihf files, you need to re-start Igor\rSave and Quit Now?"
		if (V_flag==1)
			Execute/P "Quit"
		endif
	endif

	count = strlen(errStr) ? -count : count
	return count
End
//
Static Function CopyAliases2IgorExtensionsLocal(source)
	String source

	String extension = ParseFilePath(4,source,":",0,0)
	Variable file=0, folder=0, madeOne=0

	if (StringMatch(IgorInfo(2),"Macintosh") && Stringmatch(extension,"xop"))
		file = 0					// remember, an xop is a folder on the Mac
		folder = 1
	elseif (Stringmatch(extension,"xop"))
		file = 1					// an xop is just a file on Windows
		folder = 0
	elseif (Stringmatch(extension,"ihf"))
		file = 1
		folder = 0
	else
		file = 1					// otherwise, search for anything
		folder = 1
	endif
	if (!FileFolderExists(source,file=file,folder=folder))
		Abort "The source '"+source+"' does not exist"
	endif

	String destination0 = SpecialDirPath("Igor Pro User Files",0,0,0)+IgorExtensionsFolder
	if (!FileFolderExists(destination0,folder=1))
		Abort "The destination folder '"+destination0+"' does not exist"
	endif
	String destination = destination0 + ":"+ParseFilePath(0,source,":",1,0)

	if (!AliasAlreadyThere(source,destination))
		String name=ParseFilePath(0,source,":",1,0), str
		sprintf str, "Creating alias of  \"%s\"  in  \"%s\"",name,destination
		Append2Log(str,0)
		CreateAliasShortcut/I=0/O/Z=1 source as destination
		if (V_flag)
			Abort "CreateAliasShortcut failed to make Alias"
		endif
		madeOne = 1
	endif
	return madeOne
End

// ==================================== End of XOP Aliases ===================================== //
// ============================================================================================= //



// ============================================================================================= //
// ================================ Start of Help File Aliases ================================= //

// find *.ihf files in LaueGo:doc, and put aliases to them in "...:Igor Pro 6 User Files:Igor Help Files:"
Function InstallAliasesToHelpFiles(LaueGoFullPath)
	String LaueGoFullPath
	// probably "MacComputer:Users:userName:Documents:WaveMetrics:Igor Pro 6 User Files:User Procedures:LaueGo"

	Variable n=strlen(LaueGoFullPath)
	if (n<1)
		Append2Log("ERROR -- InstallAliasesToHelpFiles, LaueGoFullPath is empty",0)
		return 1
	endif

	// make path to source help files
	String sourcePath = LaueGoFullPath + SelectString(char2num(LaueGoFullPath[n-1])==58, ":doc:", "doc:")
	String sPath = UniqueName("docPathLaueGo",12,0)
	NewPath/Q/O/Z $sPath, sourcePath
	if (V_flag)
		Append2Log("ERROR -- Unable to make Igor Help Aliases, could not make "+sPath+"  from:",0)
		Append2Log("      "+sourcePath,0)
		return 1
	endif

	// make path to destination help files
	String destPath = SpecialDirPath("Igor Pro User Files",0,0,0)+"Igor Help Files:"
	String dPath = UniqueName("UserHelpPath",12,0)
	NewPath/Q/O/Z $dPath, destPath
	if (V_flag)
		Append2Log("ERROR -- Unable to make Igor Help Aliases, could not make "+dPath+"  from:",0)
		Append2Log("      "+destPath,0)
		return 1
	endif

	String list="", name, str
	sprintf str, "Create Help file Aliases in  \"%s\"",destPath
	Append2Log(str,0)
	Variable i=0
	do													// make list of aliases that are not already there
		name = IndexedFile($sPath,i,".ihf")
		if (AliasAlreadyThere(sourcePath+name,destPath+name))
			sprintf str, "  alias is already present for    \"%s\"",name
			Append2Log(str,0)
		elseif (strlen(name))
			list += name+";"						// alias not there, add it to list
		endif
		i += 1
	while (strlen(name))


	for (i=0;i<ItemsInList(list);i+=1)		// make the aliases that are not already there
		name = StringFromList(i,list)
		sprintf str, "  Create alias for  \"%s\"   from    \"%s\"",name,sourcePath
		Append2Log(str,0)
		CreateAliasShortcut/I=0/O/Z=1 sourcePath+name as destPath+name
		if (V_flag)
			str = "CreateAliasShortcut failed to make alias to \""+name+"\""
			DoAlert 0,str
			Append2Log(str,0)
			return 1
		endif
	endfor

	return 0
End

// ================================= End of Help File Aliases ================================== //
// ============================================================================================= //




// ============================================================================================= //
// ================================= Start of Helpful Utilities ================================ //

Static Function Append2Log(line,header)
	String line
	Variable header		// flag to print header

	print line
	String name = StrVarOrDefault("LogFileName",SpecialDirPath("Desktop",0,0,0)+"LaueGoInstallLog.txt")
	Variable f
	if (!FileFolderExists(name,file=1))
		Open/D=1/M="Log file for install"/T=".txt" f as name
		name = S_fileName
	endif
	Open/A/M="Log file for install"/T=".txt" f as name

	if (strlen(S_fileName)<1)
		return 1
	endif
	String/G LogFileName = S_fileName

	FStatus f
	if (V_logEOF<2 || header)			// print a header
		String str, head=SelectString(V_logEOF>2,"","\r\r\r\r")
		head += "********************************************************************************\r"
		head += "********************************************************************************\r"
		sprintf str,"            Starting LaueGo Install   %s, %s\r\r",Secs2Date(DateTime,1), Secs2Time(DateTime,1)
		head += str
		sprintf str, "  Logging everything in history window to:\r\t\t\"%s\"\r",LogFileName
		head += str
		fprintf f,head
		print head
	endif

	fprintf f, line+"\r"
	Close f
	return 0
End


Static Function ExperimentEmpty()	// returns TRUE=1 if this is a new (empty) experiment
	Variable empty = 1
	empty = empty && StringMatch(IgorInfo(1),"Untitled")
	empty = empty && ItemsInList(FunctionList("*",";","WIN:Procedure"))<1
	empty = empty && ItemsInList(MacroList("*",";","WIN:Procedure"))<1

	String dirList = DataFolderDir(-1,$"root:")
	dirList = ReplaceString("\r",dirList,"")
	String Variables = StringByKey("VARIABLES",dirList)
	Variables = RemoveFromList("Gizmo_Error",Variables,",")
	String Strings = StringByKey("STRINGS",dirList)
	Strings = RemoveFromList("LogFileName",Strings,",")
	String Waves= StringByKey("WAVES",dirList)
	String Folders= StringByKey("FOLDERS",dirList)
	empty = empty && strlen(StringByKey("FOLDERS",dirList))<1
	empty = empty && strlen(StringByKey("WAVES",dirList))<1
	empty = empty && strlen(Variables)<1		// no pre-existing variables
	empty = empty && strlen(Strings)<1			// no pre-existing strings
	empty = empty && strlen(Waves)<1			// no pre-existing waves
	empty = empty && strlen(Folders)<1			// no pre-existing data folders

	String optStr="WIN:"+SelectString(IgorVersion()<7,"69840","4304")		// 69840 = 65536 + 4096 + 128 + 64 + 16
	String wins=WinList("*",";",optStr), str=FunctionPath("LaueGo_Install")
	str = ParseFilePath(0,str,":",1,0)
	wins = RemoveFromList(str,wins)
	wins = RemoveFromList("Procedure",wins)
	wins = RemoveFromList("Utility_JZT.ipf",wins)
	empty = empty && ItemsInList(wins)<1		// no pre-existing windows

	String MainProcedure = ProcedureText("",0,"Procedure")	// contents of main Procedure window
	Variable i=-1, returns=-1
	do
		i = strsearch(MainProcedure,"\r",i+1)
		returns += 1
	while(i>=0)
	empty = empty && strsearch(MainProcedure,"#include",0)<0
	empty = empty && returns<2

	return empty
End


Static Function AliasAlreadyThere(source,destination)
	// returns 1 if source is an alias and it is in destination
	String source				// original file, not an alias
	String destination		// desired new alias (not just a folder)
	GetFileFolderInfo/Q/Z=1 destination
	if (V_isAliasShortcut)
		Variable i
		i = char2num(source[strlen(source)-1])
		source = SelectString(i==58, source, source[0,strlen(source)-2])	// 58=":"
		i = char2num(S_aliasPath[strlen(S_aliasPath)-1])
		S_aliasPath = SelectString(i==58, S_aliasPath, S_aliasPath[0,strlen(S_aliasPath)-2])
		return StringMatch(source,S_aliasPath)
	endif
	return 0
End


Static Function/T FindFileInDirPath(start,name)		// returns full path to name searching from start
	String start			// starting point of search, a full path where to start searching
	String name				// name of file (or folder) to find

	String FullPath=""
	String pathName = UniqueName("fpath",12,0)
	String extension = ParseFilePath(4,name,":",0,0)
	extension = SelectString(strlen(extension),"",".") + extension
	NewPath/Q/Z $pathName,start
	if (V_flag)
		KillPath/Z $pathName
		return ""
	endif

	String list=IndexedFile($pathName,-1,extension), fldr
	Variable i, N=ItemsInList(list)
	for (i=0;i<N;i+=1)						// check all files at start
		if (StringMatch(IndexedFile($pathName,i,extension), name))
			KillPath/Z $pathName			// found it, all done
			return start
		endif
	endfor

	list = IndexedDir($pathName,-1,1)
	KillPath/Z $pathName
	N = ItemsInList(list)
	for (i=0;i<N;i+=1)						// check in all folders at start
		fldr = StringFromList(i,list)
		if (StringMatch(ParseFilePath(0,fldr,":",1,0), name))
			KillPath/Z $pathName			// found a match, done
			FullPath = start
			break
		endif
		if (StringMatch(ParseFilePath(4,fldr,":",0,0),"xop"))
			continue								// don't look inside of xop's (on the Mac an xop is a folder)
		endif
		FullPath = FindFileInDirPath(fldr,name)	// search inside of this folder recursively
		if (strlen(FullPath))				// found it in this folder, go back up the chain
			break
		endif
	endfor
	return FullPath
End


Static Function FileFolderExists(name,[path,file,folder])	// returns 1=exists, 0=does not exist
	String name					// partial or full file name or folder name
	String path					// optional path name, e.g. "home"
	Variable file,folder	// flags, if both set or both unset, it checks for either
	path = SelectString(ParamIsDefault(path),path,"")
	file = ParamIsDefault(file) ? 0 : file
	file = numtype(file) ? 0 : !(!file)
	folder = ParamIsDefault(folder) ? 0 : folder
	folder = numtype(folder) ? 0 : !(!folder)

	if (!file && !folder)	// check for either
		file = 1
		folder = 1
	endif

	if (strlen(path))
		PathInfo $path
		if (V_flag==0)
			return 0
		endif
		name = S_path+name	// add the path to name
	endif

	GetFileFolderInfo/Q/Z=1 name
	Variable found=0
	found = found || (file ? V_isFile : 0)
	found = found || (folder ? V_isFolder : 0)
	return found
End

// ================================== End of Helpful Utilities ================================= //
// ============================================================================================= //




// ============================================================================================= //
// ==================================== Start of Un-Zipping ==================================== //

// Short routine to select the right function to un-zip on either Macintosh or Windows
// NOTE, the Windows only unzips to the Desktop!
Function UnZipFiles(zipFile,DestFolder,[deleteZip,overWrite,printIt])
	String zipFile				// name of zip file to expand
	String DestFolder			// folder to put results (defaults to same folder as zip file"
	Variable deleteZip		// if True, delete the zip file when done (default is NO delete)
	Variable overWrite		// if True, over write existing files when un-zipping (default is NOT overwite)
	Variable printIt
	deleteZip = ParamIsDefault(deleteZip) || numtype(deleteZip) ? 0 : !(!deleteZip)
	overWrite = ParamIsDefault(overWrite) || numtype(overWrite) ? 0 : !(!overWrite)
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt

	Variable err
	if (StringMatch(IgorInfo(2),"Macintosh"))
		err = UnZipOnMac(zipFile,DestFolder, deleteZip=deleteZip, overWrite=overWrite, printIt=printIt)
	elseif (StringMatch(IgorInfo(2),"Windows"))
		err = UnzipFileOnDesktopWindows(zipFile, deleteZip)
	else
		Append2Log("ERROR -- UnZipOnMac() only works on Macintosh or Windows",0)
		err = 1
	endif
	return err
End
//
Static Function UnZipOnMac(zipFile,DestFolder,[deleteZip,overWrite,printIt])
	String zipFile				// name of zip file to expand
	String DestFolder			// folder to put results (defaults to same folder as zip file"
	Variable deleteZip		// if True, delete the zip file when done (default is NO delete)
	Variable overWrite		// if True, over write existing files when un-zipping (default is NOT overwite)
	Variable printIt
	deleteZip = ParamIsDefault(deleteZip) || numtype(deleteZip) ? 0 : deleteZip
	overWrite = ParamIsDefault(overWrite) || numtype(overWrite) ? 0 : overWrite
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt
	if (!StringMatch(IgorInfo(2),"Macintosh"))
		print "ERROR -- UnZipOnMac() only works on Macintosh"
		return 1
	endif

	// check for valid input zip file
	GetFileFolderInfo/P=home/Q/Z=1 zipFile
	if (V_Flag || !V_isFile)
		if (printIt)
			Append2Log("ERROR -- file not found, nothing done",0)
		endif
		return 1
	endif
	printIt = StringMatch(S_Path,zipFile) ? printIt : 1

	String str=""
	zipFile = S_Path
	if (!StringMatch(ParseFilePath(4,zipFile,":",0,0),"zip"))
		if (printIt)
			sprintf str, "ERROR -- \"%s\" is not a zip file\r",zipFile
			Append2Log(str,0)
		endif
		return 1
	endif

	// check for valid destination folder
	if (strlen(DestFolder)<1)
		DestFolder = ParseFilePath(1,zipFile,":",1,0)
	endif
	GetFileFolderInfo/P=home/Q/Z=1 DestFolder
	if (V_Flag || !V_isFolder)
		if (printIt)
			Append2Log("ERROR -- destination folder not found, nothing done",0)
		endif
		return 1
	endif
	DestFolder = S_Path
	printIt = StringMatch(S_Path,DestFolder) ? printIt : 1

	// get POSIX versions of paths for the shell script
	String zipFilePOSIX = ParseFilePath(5,zipFile,"/",0,0)
	String DestFolderPOSIX = ParseFilePath(5,DestFolder,"/",0,0)

	// create the shell script and execute it
	String cmd, switches=SelectString(overWrite,""," -o")
	sprintf cmd, "do shell script \"unzip %s \\\"%s\\\" -d \\\"%s\\\"\"", switches, zipFilePOSIX,DestFolderPOSIX
	ExecuteScriptText/Z cmd						//returns something only on error
	if (V_flag)
		sprintf str, "\r  ERROR -unzipping,  V_flag =",V_flag
		Append2Log(str,0)
		sprintf str, "cmd = ",ReplaceString("\n",cmd,"\r")
		Append2Log(str,0)
		sprintf str, "\r  S_value =",ReplaceString("\n",S_value,"\r")
		Append2Log(str,0)
		return V_flag									// all done, to not consider deleting the zip file
	elseif (printIt)
		sprintf str, "unzipping \"%s\"  -->  \"%s\"\r", zipFilePOSIX, DestFolderPOSIX
		Append2Log(str,0)
	endif

	// optionally delete the zip file if requested
	if (deleteZip)
		DeleteFile/M="Delete the zip file"/Z zipFile
		if (V_flag==0 && printIt)
			sprintf str, "Deleted:  \"%s\"\r", zipFile
			Append2Log(str,0)
		endif
	endif
	return V_flag
End
//
Static Function UnzipFileOnDesktopWindows(ZipFileName, deleteSource)
	string ZipFileName			// name of zip file on Desktop
	variable deleteSource		// also delete the source zip file if this is TRUE
	if (!StringMatch(IgorInfo(2),"Windows"))
		Append2Log("ERROR -- UnzipFileOnDesktopWindows() only works on Windows",0)
		return 1
	endif

	//create Path to Desktop
	NewPath/O Desktop, SpecialDirPath("Desktop", 0, 0, 0 )
	//check that the file exists
	GetFileFolderInfo /P=Desktop /Q/Z=1 ZipFileName
	if(V_Flag!=0)
		Abort "Zip file not found on Desktop"
	endif	

	//create the command file on the desktop, this is Zipjs.bat per 
	//from <a href="https://github.com/npocmaka/batch.scripts/blob/master/hybrids/jscript/zipjs.bat" title="https://github.com/npocmaka/batch.scripts/blob/master/hybrids/jscript/zipjs.bat" rel="nofollow">https://github.com/npocmaka/batch.scripts/blob/master/hybrids/jscript/zi...</a>	
	DoWindow zipjsbat
	if(V_Flag==0)
		CreateZipjsbat()	
	endif
	SaveNotebook/O/P=Desktop zipjsbat as "zipjs.bat"
	DoWindow/K zipjsbat
	//created the zipjs.bat command file which will unzip the file for us, note must kill the internal Notebook
	//or Igor will held the file and Windows will throw errors
	//now create cmd in line with
	//             zipjs.bat unzip -source C:\myDir\myZip.zip -destination C:\MyDir -keep no -force no
	// the destination folder is created by the script. 
	// -keep yes will keep the content of the zip file, -force yes will overwrite the tempfolder for the data if exists
	// be careful, -force yes will wipe out the destination, if exists, so make sure the data are directed to non-existing folder.
	string strToDesktop = SpecialDirPath("Desktop", 0, 1, 0 )
	string cmd = strToDesktop+"zipjs.bat unzip -source "
	cmd +=strToDesktop+ZipFileName+" -destination "
	cmd +=strToDesktop+"ZipFileTempFldr  -keep yes -force yes"
	ExecuteScriptText cmd
	//delete the batch file to clean up...
	DeleteFile /P=Desktop /Z  "zipjs.bat"
	if(deleteSource)
		DeleteFile /P=Desktop /Z  ZipFileName		
	endif
	return 0
End	
//
Static Function CreateZipjsbat()
	//from https://github.com/npocmaka/batch.scripts/blob/master/hybrids/jscript/zipjs.bat
	//how to use see
	//http://stackoverflow.com/questions/28043589/how-can-i-compress-zip-and-uncompress-unzip-files-and-folders-with-bat
	// this is short summary of the description there. Can unzpi, zip and do much more... 
	//
	//// unzip content of a zip to given folder.content of the zip will be not preserved (-keep no).Destination will be not overwritten (-force no)
	//call zipjs.bat unzip -source C:\myDir\myZip.zip -destination C:\MyDir -keep no -force no
	//
	//// lists content of a zip file and full paths will be printed (-flat yes)
	//call zipjs.bat list -source C:\myZip.zip\inZipDir -flat yes
	//
	//// lists content of a zip file and the content will be list as a tree (-flat no)
	//call zipjs.bat list -source C:\myZip.zip -flat no
	//
	//// prints uncompressed size in bytes
	//zipjs.bat getSize -source C:\myZip.zip
	//
	//// zips content of folder without the folder itself
	//call zipjs.bat zipDirItems -source C:\myDir\ -destination C:\MyZip.zip -keep yes -force no
	//
	//// zips file or a folder (with the folder itslelf)
	//call zipjs.bat zipItem -source C:\myDir\myFile.txt -destination C:\MyZip.zip -keep yes -force no
	//
	//// unzips only part of the zip with given path inside
	//call zipjs.bat unZipItem -source C:\myDir\myZip.zip\InzipDir\InzipFile -destination C:\OtherDir -keep no -force yes
	//call zipjs.bat unZipItem -source C:\myDir\myZip.zip\InzipDir -destination C:\OtherDir 
	//
	//// adds content to a zip file
	//call zipjs.bat addToZip -source C:\some_file -destination C:\myDir\myZip.zip\InzipDir -keep no
	//call zipjs.bat addToZip -source  C:\some_file -destination C:\myDir\myZip.zip


	String nb = "zipjsbat"
	NewNotebook/N=$nb/F=0/V=0/K=1/W=(321,81.5,820.5,376.25)
	Notebook $nb defaultTab=20, statusWidth=252
	Notebook $nb font="Arial", fSize=10, fStyle=0, textRGB=(0,0,0)
	Notebook $nb text="@if (@X)==(@Y) @end /* JScript comment\r"
	Notebook $nb text="\t@echo off\r"
	Notebook $nb text="\t\r"
	Notebook $nb text="\trem :: the first argument is the script name as it will be used for proper help message\r"
	Notebook $nb text="\tcscript //E:JScript //nologo \"%~f0\" \"%~nx0\" %*\r"
	Notebook $nb text="\r"
	Notebook $nb text="\texit /b %errorlevel%\r"
	Notebook $nb text="\t\r"
	Notebook $nb text="@if (@X)==(@Y) @end JScript comment */\r"
	Notebook $nb text="\r"
	Notebook $nb text="\r"
	Notebook $nb text="/*\r"
	Notebook $nb text="Compression/uncompression command-line tool that uses Shell.Application and WSH/Jscript -\r"
	Notebook $nb text="http://msdn.microsoft.com/en-us/library/windows/desktop/bb774085(v=vs.85).aspx\r"
	Notebook $nb text="\r"
	Notebook $nb text="Some resources That I've used:\r"
	Notebook $nb text="http://www.robvanderwoude.com/vbstech_files_zip.php\r"
	Notebook $nb text="https://code.google.com/p/jsxt/source/browse/trunk/js/win32/ZipFile.js?r=161\r"
	Notebook $nb text="\r"
	Notebook $nb text="\r"
	Notebook $nb text="\r"
	Notebook $nb text="\r"
	Notebook $nb text="UPDATE *17-03-15*\r"
	Notebook $nb text="\r"
	Notebook $nb text="Devnullius Plussed noticed a bug in ZipDirItems  and ZipItem functions (now fixed)\r"
	Notebook $nb text="And also following issues (at the moment not handled by the script):\r"
	Notebook $nb text="- if there's not enough space on the system drive (usually C:\\) the script could produce various errors "
	Notebook $nb text=", most often the script halts.\r"
	Notebook $nb text="- Folders and files that contain unicode symbols cannot be handled by Shell.Application object.\r"
	Notebook $nb text="\r"
	Notebook $nb text="UPDATE *24-03-15*\r"
	Notebook $nb text="\r"
	Notebook $nb text="Error messages are caught in waitforcount method and if shuch pops-up the script is stopped.\r"
	Notebook $nb text="As I don't know hoe to check the content of the pop-up the exact reason for the failure is not given\r"
	Notebook $nb text="but only the possible reasons.\r"
	Notebook $nb text="\r"
	Notebook $nb text="------\r"
	Notebook $nb text="It's possible to be ported for C#,Powershell and JScript.net so I'm planning to do it at some time.\r"
	Notebook $nb text="\r"
	Notebook $nb text="For sure there's a lot of room for improvements and optimization and I'm absolutely sure there are some "
	Notebook $nb text="bugs\r"
	Notebook $nb text="as the script is big enough to not have.\r"
	Notebook $nb text="\r"
	Notebook $nb text="\r"
	Notebook $nb text="\r"
	Notebook $nb text="!!!\r"
	Notebook $nb text="For suggestions contact me at - npocmaka@gmail.com\r"
	Notebook $nb text="!!!\r"
	Notebook $nb text="\r"
	Notebook $nb text="*/\r"
	Notebook $nb text="\r"
	Notebook $nb text="\r"
	Notebook $nb text="//////////////////////////////////////\r"
	Notebook $nb text="//   CONSTANTS\r"
	Notebook $nb text="\r"
	Notebook $nb text="// TODO - Shell.Application and Scripting.FileSystemObject objects could be set as global variables to a"
	Notebook $nb text="void theit creation\r"
	Notebook $nb text="// in every method.\r"
	Notebook $nb text="\r"
	Notebook $nb text="//empty zip character sequense\r"
	Notebook $nb text="var ZIP_DATA= \"PK\" + String.fromCharCode(5) + String.fromCharCode(6) + \"\\0\\0\\0\\0\\0\\0\\0\\0\\0\\0\\0\\0\\0\\0\\0\\0"
	Notebook $nb text="\\0\\0\";\r"
	Notebook $nb text="\r"
	Notebook $nb text="var SLEEP_INTERVAL=200;\r"
	Notebook $nb text="\r"
	Notebook $nb text="//copy option(s) used by Shell.Application.CopyHere/MoveHere\r"
	Notebook $nb text="var NO_PROGRESS_BAR=4;\r"
	Notebook $nb text="\r"
	Notebook $nb text="\r"
	Notebook $nb text="//oprions used for zip/unzip\r"
	Notebook $nb text="var force=true;\r"
	Notebook $nb text="var move=false;\r"
	Notebook $nb text="\r"
	Notebook $nb text="//option used for listing content of archive\r"
	Notebook $nb text="var flat=false;\r"
	Notebook $nb text="\r"
	Notebook $nb text="var source=\"\";\r"
	Notebook $nb text="var destination=\"\";\r"
	Notebook $nb text="\r"
	Notebook $nb text="var ARGS = WScript.Arguments;\r"
	Notebook $nb text="var scriptName=ARGS.Item(0);\r"
	Notebook $nb text="\r"
	Notebook $nb text="//\r"
	Notebook $nb text="//////////////////////////////////////\r"
	Notebook $nb text="\r"
	Notebook $nb text="//////////////////////////////////////\r"
	Notebook $nb text="//   ADODB.Stream extensions\r"
	Notebook $nb text="\r"
	Notebook $nb text="if ( ! this.ADODB ) {\r"
	Notebook $nb text="\tvar ADODB = {};\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="if ( ! ADODB.Stream ) {\r"
	Notebook $nb text="\tADODB.Stream = {};\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="// writes a binary data to a file\r"
	Notebook $nb text="if ( ! ADODB.Stream.writeFile ) {\r"
	Notebook $nb text="\tADODB.Stream.writeFile = function(filename, bindata)\r"
	Notebook $nb text="\t{\r"
	Notebook $nb text="        var stream = new ActiveXObject(\"ADODB.Stream\");\r"
	Notebook $nb text="        stream.Type = 2;\r"
	Notebook $nb text="        stream.Mode = 3;\r"
	Notebook $nb text="        stream.Charset =\"ASCII\";\r"
	Notebook $nb text="        stream.Open();\r"
	Notebook $nb text="        stream.Position = 0;\r"
	Notebook $nb text="        stream.WriteText(bindata);\r"
	Notebook $nb text="        stream.SaveToFile(filename, 2);\r"
	Notebook $nb text="        stream.Close();\r"
	Notebook $nb text="\t\treturn true;\r"
	Notebook $nb text="\t};\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="//\r"
	Notebook $nb text="//////////////////////////////////////\r"
	Notebook $nb text="\r"
	Notebook $nb text="//////////////////////////////////////\r"
	Notebook $nb text="//   common\r"
	Notebook $nb text="\r"
	Notebook $nb text="if ( ! this.Common ) {\r"
	Notebook $nb text="\tvar Common = {};\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="if ( ! Common.WaitForCount ) {\r"
	Notebook $nb text="\tCommon.WaitForCount = function(folderObject,targetCount,countFunction){\r"
	Notebook $nb text="\t\tvar shell = new ActiveXObject(\"Wscript.Shell\");\r"
	Notebook $nb text="\t\twhile (countFunction(folderObject) < targetCount ){\r"
	Notebook $nb text="\t\t\tWScript.Sleep(SLEEP_INTERVAL);\r"
	Notebook $nb text="\t\t\t//checks if a pop-up with error message appears while zipping\r"
	Notebook $nb text="\t\t\t//at the moment I have no idea how to read the pop-up content\r"
	Notebook $nb text="\t\t\t// to give the exact reason for failing\r"
	Notebook $nb text="\t\t\tif (shell.AppActivate(\"Compressed (zipped) Folders Error\")) {\r"
	Notebook $nb text="\t\t\t\tWScript.Echo(\"Error While zipping\");\r"
	Notebook $nb text="\t\t\t\tWScript.Echo(\"\");\r"
	Notebook $nb text="\t\t\t\tWScript.Echo(\"Possible reasons:\");\r"
	Notebook $nb text="\t\t\t\tWScript.Echo(\" -source contains filename(s) with unicode characters\");\r"
	Notebook $nb text="\t\t\t\tWScript.Echo(\" -produces zip exceeds 8gb size (or 2,5 gb for XP and 2003)\");\r"
	Notebook $nb text="\t\t\t\tWScript.Echo(\" -not enough space on system drive (usually C:\\\\)\");\r"
	Notebook $nb text="\t\t\t\tWScript.Quit(432);\r"
	Notebook $nb text="\t\t\t}\r"
	Notebook $nb text="\t\t\t\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="if ( ! Common.getParent ) {\r"
	Notebook $nb text="\tCommon.getParent = function(path){\r"
	Notebook $nb text="\t\tvar splitted=path.split(\"\\\\\");\r"
	Notebook $nb text="\t\tvar result=\"\";\r"
	Notebook $nb text="\t\tfor (var s=0;s<splitted.length-1;s++){\r"
	Notebook $nb text="\t\t\tif (s==0) {\r"
	Notebook $nb text="\t\t\t\tresult=splitted[s];\r"
	Notebook $nb text="\t\t\t} else {\r"
	Notebook $nb text="\t\t\t\tresult=result+\"\\\\\"+splitted[s];\r"
	Notebook $nb text="\t\t\t}\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t\treturn result;\r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="if ( ! Common.getName ) {\r"
	Notebook $nb text="\tCommon.getName = function(path){\r"
	Notebook $nb text="\t\tvar splitted=path.split(\"\\\\\");\r"
	Notebook $nb text="\t\treturn splitted[splitted.length-1];\r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="//file system object has a problem to create a folder with slashes at the end\r"
	Notebook $nb text="if ( ! Common.stripTrailingSlash ) {\r"
	Notebook $nb text="\tCommon.stripTrailingSlash = function(path){\r"
	Notebook $nb text="\t\twhile (path.substr(path.length - 1,path.length) == '\\\\') {\r"
	Notebook $nb text="\t\t\tpath=path.substr(0, path.length - 1);\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t\treturn path;\r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="}\r"
	Notebook $nb text="//\r"
	Notebook $nb text="//////////////////////////////////////\r"
	Notebook $nb text="\r"
	Notebook $nb text="//////////////////////////////////////\r"
	Notebook $nb text="//   Scripting.FileSystemObject extensions\r"
	Notebook $nb text="\r"
	Notebook $nb text="if (! this.Scripting) {\r"
	Notebook $nb text="\tvar Scripting={};\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="if (! Scripting.FileSystemObject) {\r"
	Notebook $nb text="\tScripting.FileSystemObject={};\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="if ( ! Scripting.FileSystemObject.DeleteItem ) {\r"
	Notebook $nb text="\tScripting.FileSystemObject.DeleteItem = function (item) \r"
	Notebook $nb text="\t{\r"
	Notebook $nb text="\t\tvar FSOObj= new ActiveXObject(\"Scripting.FileSystemObject\");\r"
	Notebook $nb text="\t\tif (FSOObj.FileExists(item)){\r"
	Notebook $nb text="\t\t\tFSOObj.DeleteFile(item);\r"
	Notebook $nb text="\t\t\treturn true;\r"
	Notebook $nb text="\t\t} else if (FSOObj.FolderExists(item) ) {\r"
	Notebook $nb text="\t\t\tFSOObj.DeleteFolder(Common.stripTrailingSlash(item));\r"
	Notebook $nb text="\t\t\treturn true;\r"
	Notebook $nb text="\t\t} else {\r"
	Notebook $nb text="\t\t\treturn false;\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="\r"
	Notebook $nb text="if ( ! Scripting.FileSystemObject.ExistsFile ) {\r"
	Notebook $nb text="\tScripting.FileSystemObject.ExistsFile = function (path)\r"
	Notebook $nb text="\t{\r"
	Notebook $nb text="\t\tvar FSOObj= new ActiveXObject(\"Scripting.FileSystemObject\");\r"
	Notebook $nb text="\t\treturn FSOObj.FileExists(path);\r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="}\r"
	Notebook $nb text="if ( !Scripting.FileSystemObject.ExistsFolder ) {\r"
	Notebook $nb text="\tScripting.FileSystemObject.ExistsFolder = function (path){\r"
	Notebook $nb text="\tvar FSOObj= new ActiveXObject(\"Scripting.FileSystemObject\");\r"
	Notebook $nb text="\t\treturn FSOObj.FolderExists(path);\r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="if ( ! Scripting.FileSystemObject.isFolder ) {\r"
	Notebook $nb text="\tScripting.FileSystemObject.isFolder = function (path){\r"
	Notebook $nb text="\tvar FSOObj= new ActiveXObject(\"Scripting.FileSystemObject\");\r"
	Notebook $nb text="\t\treturn FSOObj.FolderExists(path);\r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="\r"
	Notebook $nb text="\r"
	Notebook $nb text="if ( ! Scripting.FileSystemObject.isEmptyFolder ) {\r"
	Notebook $nb text="\tScripting.FileSystemObject.isEmptyFolder = function (path){\r"
	Notebook $nb text="\tvar FSOObj= new ActiveXObject(\"Scripting.FileSystemObject\");\r"
	Notebook $nb text="\t\tif(FSOObj.FileExists(path)){\r"
	Notebook $nb text="\t\t\treturn false;\r"
	Notebook $nb text="\t\t}else if (FSOObj.FolderExists(path)){\t\r"
	Notebook $nb text="\t\t\tvar folderObj=FSOObj.GetFolder(path);\r"
	Notebook $nb text="\t\t\tif ((folderObj.Files.Count+folderObj.SubFolders.Count)==0){\r"
	Notebook $nb text="\t\t\t\treturn true;\r"
	Notebook $nb text="\t\t\t}\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t\treturn false;\t\r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="if ( ! Scripting.FileSystemObject.CreateFolder) {\r"
	Notebook $nb text="\tScripting.FileSystemObject.CreateFolder = function (path){\r"
	Notebook $nb text="\t\tvar FSOObj= new ActiveXObject(\"Scripting.FileSystemObject\");\r"
	Notebook $nb text="\t\tFSOObj.CreateFolder(path);\r"
	Notebook $nb text="\t\treturn FSOObj.FolderExists(path);\r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="if ( ! Scripting.FileSystemObject.ExistsItem) {\r"
	Notebook $nb text="\tScripting.FileSystemObject.ExistsItem = function (path){\r"
	Notebook $nb text="\t\tvar FSOObj= new ActiveXObject(\"Scripting.FileSystemObject\");\r"
	Notebook $nb text="\t\treturn FSOObj.FolderExists(path)||FSOObj.FileExists(path);\r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="if ( ! Scripting.FileSystemObject.getFullPath) {\r"
	Notebook $nb text="\tScripting.FileSystemObject.getFullPath = function (path){\r"
	Notebook $nb text="\t\tvar FSOObj= new ActiveXObject(\"Scripting.FileSystemObject\");\r"
	Notebook $nb text="        return FSOObj.GetAbsolutePathName(path);\r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="}\r"
	Notebook $nb text="//\r"
	Notebook $nb text="//////////////////////////////////////\r"
	Notebook $nb text="\r"
	Notebook $nb text="//////////////////////////////////////\r"
	Notebook $nb text="//   Shell.Application extensions\r"
	Notebook $nb text="if ( ! this.Shell ) {\r"
	Notebook $nb text="\tvar Shell = {};\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="if (! Shell.Application ) {\r"
	Notebook $nb text="\tShell.Application={};\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="if ( ! Shell.Application.ExistsFolder ) {\r"
	Notebook $nb text="\tShell.Application.ExistsFolder = function(path){\r"
	Notebook $nb text="\t\tvar ShellObj=new ActiveXObject(\"Shell.Application\");\r"
	Notebook $nb text="\t\tvar targetObject = new Object;\r"
	Notebook $nb text="\t\tvar targetObject=ShellObj.NameSpace(path);\r"
	Notebook $nb text="\t\tif (typeof targetObject === 'undefined' || targetObject == null ){\r"
	Notebook $nb text="\t\t\treturn false;\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t\treturn true;\r"
	Notebook $nb text="\t\r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="if ( ! Shell.Application.ExistsSubItem ) {\r"
	Notebook $nb text="\tShell.Application.ExistsSubItem = function(path){\r"
	Notebook $nb text="\t\tvar ShellObj=new ActiveXObject(\"Shell.Application\");\r"
	Notebook $nb text="\t\tvar targetObject = new Object;\r"
	Notebook $nb text="\t\tvar targetObject=ShellObj.NameSpace(Common.getParent(path));\r"
	Notebook $nb text="\t\tif (typeof targetObject === 'undefined' || targetObject == null ){\r"
	Notebook $nb text="\t\t\treturn false;\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t\t\r"
	Notebook $nb text="\t\tvar subItem=targetObject.ParseName(Common.getName(path));\r"
	Notebook $nb text="\t\tif(subItem === 'undefined' || subItem == null ){\r"
	Notebook $nb text="\t\t\treturn false;\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t\treturn true;\r"
	Notebook $nb text="\t\r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="\r"
	Notebook $nb text="if ( ! Shell.Application.ItemCounterL1 ) {\r"
	Notebook $nb text="\tShell.Application.ItemCounterL1 = function(path){\r"
	Notebook $nb text="\t\tvar ShellObj=new ActiveXObject(\"Shell.Application\");\r"
	Notebook $nb text="\t\tvar targetObject = new Object;\r"
	Notebook $nb text="\t\tvar targetObject=ShellObj.NameSpace(path);\r"
	Notebook $nb text="\t\tif (targetObject != null ){\r"
	Notebook $nb text="\t\t\treturn targetObject.Items().Count;\t\r"
	Notebook $nb text="\t\t} else {\r"
	Notebook $nb text="\t\t\treturn 0;\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="}\r"
	Notebook $nb text="// shell application item.size returns the size of uncompressed state of the file.\r"
	Notebook $nb text="if ( ! Shell.Application.getSize ) {\r"
	Notebook $nb text="\tShell.Application.getSize = function(path){\r"
	Notebook $nb text="\t\tvar ShellObj=new ActiveXObject(\"Shell.Application\");\r"
	Notebook $nb text="\t\tvar targetObject = new Object;\r"
	Notebook $nb text="\t\tvar targetObject=ShellObj.NameSpace(path);\r"
	Notebook $nb text="\t\tif (! Shell.Application.ExistsFolder (path)){\r"
	Notebook $nb text="\t\t\tWScript.Echo(path + \"does not exists or the file is incorrect type.Be sure you are using full path to"
	Notebook $nb text=" the file\");\r"
	Notebook $nb text="\t\t\treturn 0;\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t\tif (typeof size === 'undefined'){\r"
	Notebook $nb text="\t\t\tvar size=0;\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t\tif (targetObject != null ){\r"
	Notebook $nb text="\t\t\t\r"
	Notebook $nb text="\t\t\tfor (var i=0; i<targetObject.Items().Count;i++){\r"
	Notebook $nb text="\t\t\t\tif(!targetObject.Items().Item(i).IsFolder){\r"
	Notebook $nb text="\t\t\t\t\tsize=size+targetObject.Items().Item(i).Size;\r"
	Notebook $nb text="\t\t\t\t} else if (targetObject.Items().Item(i).Count!=0){\r"
	Notebook $nb text="\t\t\t\t\tsize=size+Shell.Application.getSize(targetObject.Items().Item(i).Path);\r"
	Notebook $nb text="\t\t\t\t}\r"
	Notebook $nb text="\t\t\t}\r"
	Notebook $nb text="\t\t} else {\r"
	Notebook $nb text="\t\t\treturn 0;\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t\treturn size;\r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="}\r"
	Notebook $nb text="if ( ! Shell.Application.TakeAction ) {\r"
	Notebook $nb text="\tShell.Application.TakeAction = function(destination,item, move ,option){\r"
	Notebook $nb text="\t\tif(typeof destination != 'undefined' && move){\r"
	Notebook $nb text="\t\t\tdestination.MoveHere(item,option);\r"
	Notebook $nb text="\t\t} else if(typeof destination != 'undefined') {\r"
	Notebook $nb text="\t\t\tdestination.CopyHere(item,option);\r"
	Notebook $nb text="\t\t} \r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="//ProcessItem and  ProcessSubItems can be used both for zipping and unzipping\r"
	Notebook $nb text="// When an item is zipped another process is ran and the control is released\r"
	Notebook $nb text="// but when the script stops also the copying to the zipped file stops.\r"
	Notebook $nb text="// Though the zipping is transactional so a zipped files will be visible only after the zipping is done\r"
	Notebook $nb text="// and we can rely on items count when zip operation is performed. \r"
	Notebook $nb text="// Also is impossible to compress an empty folders.\r"
	Notebook $nb text="// So when it comes to zipping two additional checks are added - for empty folders and for count of item"
	Notebook $nb text="s at the \r"
	Notebook $nb text="// destination.\r"
	Notebook $nb text="\r"
	Notebook $nb text="if ( ! Shell.Application.ProcessItem ) {\r"
	Notebook $nb text="\tShell.Application.ProcessItem = function(toProcess, destination  , move ,isZipping,option){\r"
	Notebook $nb text="\t\tvar ShellObj=new ActiveXObject(\"Shell.Application\");\r"
	Notebook $nb text="\t\tdestinationObj=ShellObj.NameSpace(destination);\r"
	Notebook $nb text="\t\t\t\r"
	Notebook $nb text="\t\tif (destinationObj!= null ){\r"
	Notebook $nb text="\t\t\tif (isZipping && Scripting.FileSystemObject.isEmptyFolder(toProcess)) {\r"
	Notebook $nb text="\t\t\t\tWScript.Echo(toProcess +\" is an empty folder and will be not processed\");\r"
	Notebook $nb text="\t\t\t\treturn;\r"
	Notebook $nb text="\t\t\t}\r"
	Notebook $nb text="\t\t\tShell.Application.TakeAction(destinationObj,toProcess, move ,option);\r"
	Notebook $nb text="\t\t\tvar destinationCount=Shell.Application.ItemCounterL1(destination);\r"
	Notebook $nb text="\t\t\tvar final_destination=destination + \"\\\\\" + Common.getName(toProcess);\r"
	Notebook $nb text="\t\t\t\r"
	Notebook $nb text="\t\t\tif (isZipping && !Shell.Application.ExistsSubItem(final_destination)) {\r"
	Notebook $nb text="\t\t\t\tCommon.WaitForCount(destination\r"
	Notebook $nb text="\t\t\t\t\t,destinationCount+1,Shell.Application.ItemCounterL1);\r"
	Notebook $nb text="\t\t\t} else if (isZipping && Shell.Application.ExistsSubItem(final_destination)){\r"
	Notebook $nb text="\t\t\t\tWScript.Echo(final_destination + \" already exists and task cannot be completed\");\r"
	Notebook $nb text="\t\t\t\treturn;\r"
	Notebook $nb text="\t\t\t}\r"
	Notebook $nb text="\t\t}\t\r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="\r"
	Notebook $nb text="\r"
	Notebook $nb text="if ( ! Shell.Application.ProcessSubItems ) {\r"
	Notebook $nb text="\tShell.Application.ProcessSubItems = function(toProcess, destination  , move ,isZipping ,option){\r"
	Notebook $nb text="\t\tvar ShellObj=new ActiveXObject(\"Shell.Application\");\r"
	Notebook $nb text="\t\tvar destinationObj=ShellObj.NameSpace(destination);\r"
	Notebook $nb text="\t\tvar toItemsToProcess=new Object;\r"
	Notebook $nb text="\t\ttoItemsToProcess=ShellObj.NameSpace(toProcess).Items();\r"
	Notebook $nb text="\t\t\t\r"
	Notebook $nb text="\t\tif (destinationObj!= null ){\r"
	Notebook $nb text="\t\t\t\t\t\t\r"
	Notebook $nb text="\t\t\tfor (var i=0;i<toItemsToProcess.Count;i++) {\r"
	Notebook $nb text="\t\t\t\t\r"
	Notebook $nb text="\t\t\t\tif (isZipping && Scripting.FileSystemObject.isEmptyFolder(toItemsToProcess.Item(i).Path)){\r"
	Notebook $nb text="\r"
	Notebook $nb text="\t\t\t\t\tWScript.Echo(\"\");\r"
	Notebook $nb text="\t\t\t\t\tWScript.Echo(toItemsToProcess.Item(i).Path + \" is empty and will be not processed\");\r"
	Notebook $nb text="\t\t\t\t\tWScript.Echo(\"\");\r"
	Notebook $nb text="\r"
	Notebook $nb text="\t\t\t\t} else {\r"
	Notebook $nb text="\t\t\t\t\tShell.Application.TakeAction(destinationObj,toItemsToProcess.Item(i),move,option);\r"
	Notebook $nb text="\t\t\t\t\tvar destinationCount=Shell.Application.ItemCounterL1(destination);\r"
	Notebook $nb text="\t\t\t\t\tif (isZipping) {\r"
	Notebook $nb text="\t\t\t\t\t\tCommon.WaitForCount(destination,destinationCount+1,Shell.Application.ItemCounterL1);\r"
	Notebook $nb text="\t\t\t\t\t}\r"
	Notebook $nb text="\t\t\t\t}\r"
	Notebook $nb text="\t\t\t}\t\r"
	Notebook $nb text="\t\t}\t\r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="if ( ! Shell.Application.ListItems ) {\r"
	Notebook $nb text="\tShell.Application.ListItems = function(parrentObject){\r"
	Notebook $nb text="\t\tvar ShellObj=new ActiveXObject(\"Shell.Application\");\r"
	Notebook $nb text="\t\tvar targetObject = new Object;\r"
	Notebook $nb text="\t\tvar targetObject=ShellObj.NameSpace(parrentObject);\r"
	Notebook $nb text="\r"
	Notebook $nb text="\t\tif (! Shell.Application.ExistsFolder (parrentObject)){\r"
	Notebook $nb text="\t\t\tWScript.Echo(parrentObject + \"does not exists or the file is incorrect type.Be sure the full path the"
	Notebook $nb text=" path is used\");\r"
	Notebook $nb text="\t\t\treturn;\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t\tif (typeof initialSCount == 'undefined') {\r"
	Notebook $nb text="\t\t\tinitialSCount=(parrentObject.split(\"\\\\\").length-1);\r"
	Notebook $nb text="\t\t\tWScript.Echo(parrentObject);\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t\t\r"
	Notebook $nb text="\t\tvar spaces=function(path){\r"
	Notebook $nb text="\t\t\tvar SCount=(path.split(\"\\\\\").length-1)-initialSCount;\r"
	Notebook $nb text="\t\t\tvar s=\"\";\r"
	Notebook $nb text="\t\t\tfor (var i=0;i<=SCount;i++) {\r"
	Notebook $nb text="\t\t\t\ts=\" \"+s;\r"
	Notebook $nb text="\t\t\t}\r"
	Notebook $nb text="\t\t\treturn s;\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t\t\r"
	Notebook $nb text="\t\tvar printP = function (item,end){\r"
	Notebook $nb text="\t\t\tif (flat) {\r"
	Notebook $nb text="\t\t\t\tWScript.Echo(targetObject.Items().Item(i).Path+end);\r"
	Notebook $nb text="\t\t\t}else{\r"
	Notebook $nb text="\t\t\t\tWScript.Echo( spaces(targetObject.Items().Item(i).Path)+targetObject.Items().Item(i).Name+end);\r"
	Notebook $nb text="\t\t\t}\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\r"
	Notebook $nb text="\t\tif (targetObject != null ){\r"
	Notebook $nb text="\t\t\tvar folderPath=\"\";\r"
	Notebook $nb text="\r"
	Notebook $nb text="\t\t\t\tfor (var i=0; i<targetObject.Items().Count;i++) {\r"
	Notebook $nb text="\t\t\t\t\tif(targetObject.Items().Item(i).IsFolder && targetObject.Items().Item(i).Count==0 ){\r"
	Notebook $nb text="\t\t\t\t\t\tprintP(targetObject.Items().Item(i),\"\\\\\");\r"
	Notebook $nb text="\t\t\t\t\t} else if (targetObject.Items().Item(i).IsFolder){\r"
	Notebook $nb text="\t\t\t\t\t\tfolderPath=parrentObject+\"\\\\\"+targetObject.Items().Item(i).Name;\r"
	Notebook $nb text="\t\t\t\t\t\tprintP(targetObject.Items().Item(i),\"\\\\\")\r"
	Notebook $nb text="\t\t\t\t\t\tShell.Application.ListItems(folderPath);\t\t\t\t\t\t\r"
	Notebook $nb text="\r"
	Notebook $nb text="\t\t\t\t\t} else {\r"
	Notebook $nb text="\t\t\t\t\t\tprintP(targetObject.Items().Item(i),\"\")\r"
	Notebook $nb text="\t\t\t\t\t\t\r"
	Notebook $nb text="\t\t\t\t\t}\r"
	Notebook $nb text="\t\t\t\t}\r"
	Notebook $nb text="\r"
	Notebook $nb text="\t\t\t}\r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="//\r"
	Notebook $nb text="//////////////////////////////////////\r"
	Notebook $nb text="\r"
	Notebook $nb text="\r"
	Notebook $nb text="//////////////////////////////////////\r"
	Notebook $nb text="//     ZIP Utils\r"
	Notebook $nb text="\r"
	Notebook $nb text="if ( ! this.ZIPUtils ) {\r"
	Notebook $nb text="\tvar ZIPUtils = {};\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="\r"
	Notebook $nb text="if ( ! this.ZIPUtils.ZipItem) {\t\r"
	Notebook $nb text="\tZIPUtils.ZipItem = function(source, destination ) {\r"
	Notebook $nb text="\t\tif (!Scripting.FileSystemObject.ExistsFolder(source)) {\r"
	Notebook $nb text="\t\t\tWScript.Echo(\"\");\r"
	Notebook $nb text="\t\t\tWScript.Echo(\"file \" + source +\" does not exist\");\r"
	Notebook $nb text="\t\t\tWScript.Quit(2);\t\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t\t\r"
	Notebook $nb text="\t\tif (Scripting.FileSystemObject.ExistsFile(destination) && force ) {\r"
	Notebook $nb text="\t\t\tScripting.FileSystemObject.DeleteItem(destination);\r"
	Notebook $nb text="\t\t\tADODB.Stream.writeFile(destination,ZIP_DATA);\r"
	Notebook $nb text="\t\t} else if (!Scripting.FileSystemObject.ExistsFile(destination)) {\r"
	Notebook $nb text="\t\t\tADODB.Stream.writeFile(destination,ZIP_DATA);\r"
	Notebook $nb text="\t\t} else {\r"
	Notebook $nb text="\t\t\tWScript.Echo(\"Destination \"+destination+\" already exists.Operation will be aborted\");\r"
	Notebook $nb text="\t\t\tWScript.Quit(15);\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t\tsource=Scripting.FileSystemObject.getFullPath(source);\r"
	Notebook $nb text="\t\tdestination=Scripting.FileSystemObject.getFullPath(destination);\r"
	Notebook $nb text="\t\t\r"
	Notebook $nb text="\t\tShell.Application.ProcessItem(source,destination,move,true ,NO_PROGRESS_BAR);\r"
	Notebook $nb text="\r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="\r"
	Notebook $nb text="if ( ! this.ZIPUtils.ZipDirItems) {\t\r"
	Notebook $nb text="\tZIPUtils.ZipDirItems = function(source, destination ) {\r"
	Notebook $nb text="\t\tif (!Scripting.FileSystemObject.ExistsFolder(source)) {\r"
	Notebook $nb text="\t\t\tWScript.Echo();\r"
	Notebook $nb text="\t\t\tWScript.Echo(\"file \" + source +\" does not exist\");\r"
	Notebook $nb text="\t\t\tWScript.Quit(2);\t\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t\tif (Scripting.FileSystemObject.ExistsFile(destination) && force ) {\r"
	Notebook $nb text="\t\t\tScripting.FileSystemObject.DeleteItem(destination);\r"
	Notebook $nb text="\t\t\tADODB.Stream.writeFile(destination,ZIP_DATA);\r"
	Notebook $nb text="\t\t} else if (!Scripting.FileSystemObject.ExistsFile(destination)) {\r"
	Notebook $nb text="\t\t\tADODB.Stream.writeFile(destination,ZIP_DATA);\r"
	Notebook $nb text="\t\t} else {\r"
	Notebook $nb text="\t\t\tWScript.Echo(\"Destination \"+destination+\" already exists.Operation will be aborted\");\r"
	Notebook $nb text="\t\t\tWScript.Quit(15);\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t\t\r"
	Notebook $nb text="\t\tsource=Scripting.FileSystemObject.getFullPath(source);\r"
	Notebook $nb text="\t\tdestination=Scripting.FileSystemObject.getFullPath(destination);\r"
	Notebook $nb text="\t\t\r"
	Notebook $nb text="\t\tShell.Application.ProcessSubItems(source, destination, move ,true,NO_PROGRESS_BAR);\r"
	Notebook $nb text="\t\t\r"
	Notebook $nb text="\t\tif (move){\r"
	Notebook $nb text="\t\t\tScripting.FileSystemObject.DeleteItem(source);\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="if ( ! this.ZIPUtils.Unzip) {\t\r"
	Notebook $nb text="\tZIPUtils.Unzip = function(source, destination ) {\r"
	Notebook $nb text="\t\tif(!Shell.Application.ExistsFolder(source) ){\r"
	Notebook $nb text="\t\t\tWScript.Echo(\"Either the target does not exist or is not a correct type\");\r"
	Notebook $nb text="\t\t\treturn;\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t\t\r"
	Notebook $nb text="\t\tif (Scripting.FileSystemObject.ExistsItem(destination) && force ) {\r"
	Notebook $nb text="\t\t\tScripting.FileSystemObject.DeleteItem(destination);\r"
	Notebook $nb text="\t\t} else if (Scripting.FileSystemObject.ExistsItem(destination)){\r"
	Notebook $nb text="\t\t\tWScript.Echo(\"Destination \" + destination + \" already exists\");\r"
	Notebook $nb text="\t\t\treturn;\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t\tScripting.FileSystemObject.CreateFolder(destination);\r"
	Notebook $nb text="\t\tsource=Scripting.FileSystemObject.getFullPath(source);\r"
	Notebook $nb text="\t\tdestination=Scripting.FileSystemObject.getFullPath(destination);\r"
	Notebook $nb text="\t\tShell.Application.ProcessSubItems(source, destination, move ,false,NO_PROGRESS_BAR);\t\r"
	Notebook $nb text="\t\t\r"
	Notebook $nb text="\t\tif (move){\r"
	Notebook $nb text="\t\t\tScripting.FileSystemObject.DeleteItem(source);\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="    }\t\t\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="if ( ! this.ZIPUtils.AddToZip) {\r"
	Notebook $nb text="\tZIPUtils.AddToZip = function(source, destination ) {\r"
	Notebook $nb text="\t\tif(!Shell.Application.ExistsFolder(destination)) {\r"
	Notebook $nb text="\t\t\tWScript.Echo(destination +\" is not valid path to/within zip.Be sure you are not using relative paths\""
	Notebook $nb text=");\r"
	Notebook $nb text="\t\t\tWscript.Exit(\"101\");\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t\tif(!Scripting.FileSystemObject.ExistsItem(source)){\r"
	Notebook $nb text="\t\t\tWScript.Echo(source +\" does not exist\");\r"
	Notebook $nb text="\t\t\tWscript.Exit(\"102\");\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t\tsource=Scripting.FileSystemObject.getFullPath(source);\r"
	Notebook $nb text="\t\tShell.Application.ProcessItem(source,destination,move,true ,NO_PROGRESS_BAR); \r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="if ( ! this.ZIPUtils.UnzipItem) {\t\r"
	Notebook $nb text="\tZIPUtils.UnzipItem = function(source, destination ) {\r"
	Notebook $nb text="\r"
	Notebook $nb text="\t\tif(!Shell.Application.ExistsSubItem(source)){\r"
	Notebook $nb text="\t\t\tWScript.Echo(source + \":Either the target does not exist or is not a correct type\");\r"
	Notebook $nb text="\t\t\treturn;\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t\t\r"
	Notebook $nb text="\t\tif (Scripting.FileSystemObject.ExistsItem(destination) && force ) {\r"
	Notebook $nb text="\t\t\tScripting.FileSystemObject.DeleteItem(destination);\r"
	Notebook $nb text="\t\t} else if (Scripting.FileSystemObject.ExistsItem(destination)){\r"
	Notebook $nb text="\t\t\tWScript.Echo(destination+\" - Destination already exists\");\r"
	Notebook $nb text="\t\t\treturn;\r"
	Notebook $nb text="\t\t} \r"
	Notebook $nb text="\t\t\r"
	Notebook $nb text="\t\tScripting.FileSystemObject.CreateFolder(destination);\r"
	Notebook $nb text="\t\tdestination=Scripting.FileSystemObject.getFullPath(destination);\r"
	Notebook $nb text="\t\tShell.Application.ProcessItem(source, destination, move ,false,NO_PROGRESS_BAR);\r"
	Notebook $nb text="\t\t                            \r"
	Notebook $nb text="    }\t\t\r"
	Notebook $nb text="}\r"
	Notebook $nb text="if ( ! this.ZIPUtils.getSize) {\t\r"
	Notebook $nb text="\tZIPUtils.getSize = function(path) {\r"
	Notebook $nb text="\t\t// first getting a full path to the file is attempted\r"
	Notebook $nb text="\t\t// as it's required by shell.application\r"
	Notebook $nb text="\t\t// otherwise is assumed that a file within a zip \r"
	Notebook $nb text="\t\t// is aimed\r"
	Notebook $nb text="\t\t\r"
	Notebook $nb text="\t\t//TODO - find full path even if the path points to internal for the \r"
	Notebook $nb text="\t\t// zip directory\r"
	Notebook $nb text="\t\t\r"
	Notebook $nb text="\t\tif (Scripting.FileSystemObject.ExistsFile(path)){\r"
	Notebook $nb text="\t\t\tpath=Scripting.FileSystemObject.getFullPath(path);\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t\tWScript.Echo(Shell.Application.getSize(path));\r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="if ( ! this.ZIPUtils.list) {\t\r"
	Notebook $nb text="\tZIPUtils.list = function(path) {\r"
	Notebook $nb text="\t\t// first getting a full path to the file is attempted\r"
	Notebook $nb text="\t\t// as it's required by shell.application\r"
	Notebook $nb text="\t\t// otherwise is assumed that a file within a zip \r"
	Notebook $nb text="\t\t// is aimed\r"
	Notebook $nb text="\t\t\r"
	Notebook $nb text="\t\t//TODO - find full path even if the path points to internal for the \r"
	Notebook $nb text="\t\t// zip directory\r"
	Notebook $nb text="\t\t\r"
	Notebook $nb text="\t\t// TODO - optional printing of each file uncompressed size\r"
	Notebook $nb text="\t\t\r"
	Notebook $nb text="\t\tif (Scripting.FileSystemObject.ExistsFile(path)){\r"
	Notebook $nb text="\t\t\tpath=Scripting.FileSystemObject.getFullPath(path);\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t\tShell.Application.ListItems(path);\r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="//\r"
	Notebook $nb text="//////////////////////////////////////\r"
	Notebook $nb text="\r"
	Notebook $nb text="/////////////////////////////////////\r"
	Notebook $nb text="//   parsing'n'running\r"
	Notebook $nb text="function printHelp(){\r"
	Notebook $nb text="\r"
	Notebook $nb text="\tWScript.Echo( scriptName + \" list -source zipFile [-flat yes|no]\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\tlist the content of a zip file\");\r"
	Notebook $nb text="\tWScript.Echo( \"\tzipFile - absolute path to the zip file\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\t\tcould be also a directory or a directory inside a zip file or\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\t\tor a .cab file or an .iso file\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t-flat - indicates if the structure of the zip will be printed as tree\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\t\tor with absolute paths (-flat yes).Default is yes.\");\r"
	Notebook $nb text="\tWScript.Echo( \"Example:\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\" + scriptName + \" list -source C:\\\\myZip.zip -flat no\" );\r"
	Notebook $nb text="\tWScript.Echo( \"\t\" + scriptName + \" list -source C:\\\\myZip.zip\\\\inZipDir -flat yes\" );\r"
	Notebook $nb text="\t\r"
	Notebook $nb text="\tWScript.Echo( scriptName + \" getSize -source zipFile\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\tprints uncompressed size of the zipped file in bytes\");\r"
	Notebook $nb text="\tWScript.Echo( \"\tzipFile - absolute path to the zip file\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\t\tcould be also a directory or a directory inside a zip file or\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\t\tor a .cab file or an .iso file\");\r"
	Notebook $nb text="\tWScript.Echo( \"Example:\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\" + scriptName + \" getSize -source C:\\\\myZip.zip\" );\r"
	Notebook $nb text="\t\r"
	Notebook $nb text="\tWScript.Echo( scriptName + \" zipDirItems -source source_dir -destination destination.zip [-force yes|no"
	Notebook $nb text="] [-keep yes|no]\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\tzips the content of given folder without the folder itself \");\r"
	Notebook $nb text="\tWScript.Echo( \"\tsource_dir - path to directory which content will be compressed\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\tEmpty folders in the source directory will be ignored\");\r"
	Notebook $nb text="\tWScript.Echo( \"\tdestination.zip - path/name  of the zip file that will be created\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t-force - indicates if the destination will be overwritten if already exists.\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\tdefault is yes\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t-keep - indicates if the source content will be moved or just copied/kept.\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\tdefault is yes\");\r"
	Notebook $nb text="\tWScript.Echo( \"Example:\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\" + scriptName + \" zipDirItems -source C:\\\\myDir\\\\ -destination C:\\\\MyZip.zip -keep yes"
	Notebook $nb text=" -force no\" );\r"
	Notebook $nb text="\t\r"
	Notebook $nb text="\tWScript.Echo( scriptName + \" zipItem -source item -destination destination.zip [-force yes|no] [-keep y"
	Notebook $nb text="es|no]\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\tzips file or directory to a destination.zip file \");\r"
	Notebook $nb text="\tWScript.Echo( \"\titem - path to file or directory which content will be compressed\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\tIf points to an empty folder it will be ignored\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\tIf points to a folder it also will be included in the zip file alike zipdiritems comma"
	Notebook $nb text="nd\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\tEventually zipping a folder in this way will be faster as it does not process every el"
	Notebook $nb text="ement one by one\");\r"
	Notebook $nb text="\tWScript.Echo( \"\tdestination.zip - path/name  of the zip file that will be created\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t-force - indicates if the destination will be overwritten if already exists.\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\tdefault is yes\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t-keep - indicates if the source content will be moved or just copied/kept.\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\tdefault is yes\");\r"
	Notebook $nb text="\tWScript.Echo( \"Example:\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\" + scriptName + \" zipItem -source C:\\\\myDir\\\\myFile.txt -destination C:\\\\MyZip.zip -ke"
	Notebook $nb text="ep yes -force no\" );\r"
	Notebook $nb text="\t\r"
	Notebook $nb text="\tWScript.Echo( scriptName + \" unzip -source source.zip -destination destination_dir [-force yes|no] [-ke"
	Notebook $nb text="ep yes|no]\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\tunzips the content of a zip file to a given directory \");\r"
	Notebook $nb text="\tWScript.Echo( \"\tsource - path to the zip file that will be expanded\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\tEventually .iso , .cab or even an ordinary directory can be used as a source\");\r"
	Notebook $nb text="\tWScript.Echo( \"\tdestination_dir - path to directory where unzipped items will be stored\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t-force - indicates if the destination will be overwritten if already exists.\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\tdefault is yes\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t-keep - indicates if the source content will be moved or just copied/kept.\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\tdefault is yes\");\r"
	Notebook $nb text="\tWScript.Echo( \"Example:\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\" + scriptName + \" unzip -source C:\\\\myDir\\\\myZip.zip -destination C:\\\\MyDir -keep no -"
	Notebook $nb text="force no\" );\r"
	Notebook $nb text="\t\r"
	Notebook $nb text="\tWScript.Echo( scriptName + \" unZipItem -source source.zip -destination destination_dir [-force yes|no] "
	Notebook $nb text="[-keep yes|no]\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\tunzips  a single within a given zip file to a destination directory\");\r"
	Notebook $nb text="\tWScript.Echo( \"\tsource - path to the file/folcer within a zip  that will be expanded\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\tEventually .iso , .cab or even an ordinary directory can be used as a source\");\r"
	Notebook $nb text="\tWScript.Echo( \"\tdestination_dir - path to directory where unzipped item will be stored\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t-force - indicates if the destination directory will be overwritten if already exists.\""
	Notebook $nb text=");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\tdefault is yes\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t-keep - indicates if the source content will be moved or just copied/kept.\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\tdefault is yes\");\r"
	Notebook $nb text="\tWScript.Echo( \"Example:\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\" + scriptName + \" unZipItem -source C:\\\\myDir\\\\myZip.zip\\\\InzipDir\\\\InzipFile -destina"
	Notebook $nb text="tion C:\\\\OtherDir -keep no -force yes\" );\r"
	Notebook $nb text="\tWScript.Echo( \"\t\" + scriptName + \" unZipItem -source C:\\\\myDir\\\\myZip.zip\\\\InzipDir -destination C:\\\\Ot"
	Notebook $nb text="herDir \" );\r"
	Notebook $nb text="\t\r"
	Notebook $nb text="\tWScript.Echo( scriptName + \" addToZip -source sourceItem -destination destination.zip  [-keep yes|no]\")"
	Notebook $nb text=";\r"
	Notebook $nb text="\tWScript.Echo( \"\t\tadds file or folder to already exist zip file\");\r"
	Notebook $nb text="\tWScript.Echo( \"\tsource - path to the item that will be processed\");\r"
	Notebook $nb text="\tWScript.Echo( \"\tdestination_zip - path to the zip where the item will be added\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t-keep - indicates if the source content will be moved or just copied/kept.\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\tdefault is yes\");\r"
	Notebook $nb text="\tWScript.Echo( \"Example:\");\r"
	Notebook $nb text="\tWScript.Echo( \"\t\" + scriptName + \" addToZip -source C:\\\\some_file -destination C:\\\\myDir\\\\myZip.zip\\\\In"
	Notebook $nb text="zipDir -keep no \" );\r"
	Notebook $nb text="\tWScript.Echo( \"\t\" + scriptName + \" addToZip -source  C:\\\\some_file -destination C:\\\\myDir\\\\myZip.zip \" "
	Notebook $nb text=");\r"
	Notebook $nb text="\t\r"
	Notebook $nb text="\tWScript.Echo( \"\tby Vasil \\\"npocmaka\\\" Arnaudov - npocmaka@gmail.com\" );\r"
	Notebook $nb text="\tWScript.Echo( \"\tver 0.1 \" );\r"
	Notebook $nb text="\tWScript.Echo( \"\tlatest version could be found here https://github.com/npocmaka/batch.scripts/blob/maste"
	Notebook $nb text="r/hybrids/jscript/zipjs.bat\" );\r"
	Notebook $nb text="\t\r"
	Notebook $nb text="\t\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="function parseArguments(){\r"
	Notebook $nb text="\tif (WScript.Arguments.Length==1 || WScript.Arguments.Length==2 || ARGS.Item(1).toLowerCase() == \"-help\""
	Notebook $nb text=" ||  ARGS.Item(1).toLowerCase() == \"-h\" ) {\r"
	Notebook $nb text="\t\tprintHelp();\r"
	Notebook $nb text="\t\tWScript.Quit(0);\r"
	Notebook $nb text="   }\r"
	Notebook $nb text="   \r"
	Notebook $nb text="   //all arguments are key-value pairs plus one for script name and action taken - need to be even numbe"
	Notebook $nb text="r\r"
	Notebook $nb text="\tif (WScript.Arguments.Length % 2 == 1 ) {\r"
	Notebook $nb text="\t\tWScript.Echo(\"Illegal arguments \");\r"
	Notebook $nb text="\t\tprintHelp();\r"
	Notebook $nb text="\t\tWScript.Quit(1);\r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="\t\r"
	Notebook $nb text="\t//ARGS\r"
	Notebook $nb text="\tfor(var arg = 2 ; arg<ARGS.Length-1;arg=arg+2) {\r"
	Notebook $nb text="\t\tif (ARGS.Item(arg) == \"-source\") {\r"
	Notebook $nb text="\t\t\tsource = ARGS.Item(arg +1);\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t\tif (ARGS.Item(arg) == \"-destination\") {\r"
	Notebook $nb text="\t\t\tdestination = ARGS.Item(arg +1);\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t\tif (ARGS.Item(arg).toLowerCase() == \"-keep\" && ARGS.Item(arg +1).toLowerCase() == \"no\") {\r"
	Notebook $nb text="\t\t\tmove=true;\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t\tif (ARGS.Item(arg).toLowerCase() == \"-force\" && ARGS.Item(arg +1).toLowerCase() == \"no\") {\r"
	Notebook $nb text="\t\t\tforce=false;\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t\tif (ARGS.Item(arg).toLowerCase() == \"-flat\" && ARGS.Item(arg +1).toLowerCase() == \"yes\") {\r"
	Notebook $nb text="\t\t\tflat=true;\r"
	Notebook $nb text="\t\t}\r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="\t\r"
	Notebook $nb text="\tif (source == \"\"){\r"
	Notebook $nb text="\t\tWScript.Echo(\"Source not given\");\r"
	Notebook $nb text="\t\tprintHelp();\r"
	Notebook $nb text="\t\tWScript.Quit(59);\r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="var checkDestination=function(){\r"
	Notebook $nb text="\tif (destination == \"\"){\r"
	Notebook $nb text="\t\tWScript.Echo(\"Destination not given\");\r"
	Notebook $nb text="\t\tprintHelp();\r"
	Notebook $nb text="\t\tWScript.Quit(65);\r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="}\r"
	Notebook $nb text="\r"
	Notebook $nb text="var main=function(){\r"
	Notebook $nb text="\tparseArguments();\r"
	Notebook $nb text="\tswitch (ARGS.Item(1).toLowerCase()) {\r"
	Notebook $nb text="\tcase \"list\":\r"
	Notebook $nb text="\t\tZIPUtils.list(source);\r"
	Notebook $nb text="\t\tbreak;\r"
	Notebook $nb text="\tcase \"getsize\":\r"
	Notebook $nb text="\t\tZIPUtils.getSize(source);\r"
	Notebook $nb text="\t\tbreak;\r"
	Notebook $nb text="\tcase \"zipdiritems\":\r"
	Notebook $nb text="\t\tcheckDestination();\r"
	Notebook $nb text="\t\tZIPUtils.ZipDirItems(source,destination);\r"
	Notebook $nb text="\t\tbreak;\r"
	Notebook $nb text="\tcase \"zipitem\":\r"
	Notebook $nb text="\t\tcheckDestination();\r"
	Notebook $nb text="\t\tZIPUtils.ZipDirItems(source,destination);\r"
	Notebook $nb text="\t\tbreak;\r"
	Notebook $nb text="\tcase \"unzip\":\r"
	Notebook $nb text="\t\tcheckDestination();\r"
	Notebook $nb text="\t\tZIPUtils.Unzip(source,destination);\r"
	Notebook $nb text="\t\tbreak;\r"
	Notebook $nb text="\tcase \"unzipitem\":\r"
	Notebook $nb text="\t\tcheckDestination();\r"
	Notebook $nb text="\t\tZIPUtils.UnzipItem(source,destination);\r"
	Notebook $nb text="\t\tbreak;\r"
	Notebook $nb text="\tcase \"addtozip\":\r"
	Notebook $nb text="\t\tcheckDestination();\r"
	Notebook $nb text="\t\tZIPUtils.AddToZip(source,destination);\r"
	Notebook $nb text="\t\tbreak;\r"
	Notebook $nb text="\tdefault:\r"
	Notebook $nb text="\t\tWScript.Echo(\"No valid switch has been passed\");\r"
	Notebook $nb text="\t\tprintHelp();\r"
	Notebook $nb text="\t\t\r"
	Notebook $nb text="\t}\r"
	Notebook $nb text="\t\r"
	Notebook $nb text="}\r"
	Notebook $nb text="main();\r"
	Notebook $nb text="//\r"
	Notebook $nb text="//////////////////////////////////////"
end

// ===================================== End of Un-Zipping ===================================== //
// ============================================================================================= //

