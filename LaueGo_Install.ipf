#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName=LaueGoInstall
#pragma version = 0.03
// #pragma hide = 1
Constant LaueGo_Install_Test = 0
Static strConstant gitHubArchiveURL = "http://github.com/34IDE/LaueGo/archive/master.zip"
Static strConstant JZTArchiveURL = "http://sector33.xray.aps.anl.gov/~tischler/igor/LaueGo.zip"
Static Constant minimumArchiveLen = 16777216	// 16MB (16*1024*1024)

strConstant xopList = "HDF5.xop;HDF5 Help.ihf;"		// a list of xop's to install, this can be overidden
//Static strConstant = xopList = "HDF5.xop;HDF5 Help.ihf;HFSAndPosix.xop;HFSAndPosix Help.ihf;"


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

	Variable empty = ExperimentEmpty()
	String str = "***"+SelectString(empty,"This Expermiment appears to NOT be empty","This is considered to be an EMPTY Experiment")
	Append2Log(str,0)
	String UserPath = SpecialDirPath("Igor Pro User Files",0,0,0)+"User Procedures:"

	if (!LaueGoInstall#FileFolderExists(UserPath,folder=1))
		DoAlert 0,"The \"User Path\", does not exist, quit Igor and try again"
		return 1
	endif
	String LaueGoFullPath = UserPath+"LaueGo"		// official location of the LaueGo folder

	Variable existing = LaueGoInstall#FileFolderExists(LaueGoFullPath,folder=1)
	String oldFolderDestination=""
	if (existing)
		sprintf oldFolderDestination, "LaueGo_%sT%s",Secs2Date(DateTime,-2),Secs2Time(DateTime,3)
		oldFolderDestination = SpecialDirPath("Desktop",0,0,0)+CleanupName(oldFolderDestination,0)
		if (LaueGoInstall#FileFolderExists(oldFolderDestination,folder=1))
			sprintf str, "Old folder named \"%s\"  already exists!",oldFolderDestination
			DoAlert 0, str
			Append2Log(str,0)
			return 1
		endif
	endif

	// identify the new LaueGo folder
	DoAlert/T="find NEW LaueGo" 2, "Download newest version from web?"
	String newFullPath="", name
	if (V_flag==1)										// download LaueGo.zip, expand it, and put it on Desktop
		newFullPath = FetchLaueGoArchive("")
	else													// ask user to identify and existing (presumably new) LaueGo folder to use
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

	name = ParseFilePath(0,newFullPath,":",1,0)
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

	// install the alias to "always first.ipf"
	String alwaysFirstSource=LaueGoFullPath+":always:always first.ipf"
	String alwaysFirstDestination = SpecialDirPath("Igor Pro User Files",0,0,0)+"Igor Procedures:always first.ipf"
	if (LaueGoInstall#AliasAlreadyThere(alwaysFirstSource,alwaysFirstDestination))
		Append2Log("The alias to \":always:always first.ipf\" is already in \"Igor Procedures\", no action taken",0)
	else
		Append2Log("  about to create the alias",0)
		sprintf str, "  alwaysFirstSource =",alwaysFirstSource
		Append2Log(str,0)
		sprintf str,"  always first destination   ",alwaysFirstDestination
		Append2Log(str,0)
		CreateAliasShortcut/I=0/O/Z=1 alwaysFirstSource as alwaysFirstDestination
		if (V_flag)
			str = "CreateAliasShortcut failed to make alias to \"always first.ipf\""
			DoAlert 0,str
			Append2Log(str,0)
			return 1
		endif
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

	if (empty)
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
//
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

Static Function ExperimentEmpty()	// returns TRUE=1 if this is a new experiment
	Variable empty = 1
	empty = empty && StringMatch(IgorInfo(1),"Untitled")
	empty = empty && ItemsInList(FunctionList("*",";","WIN:Procedure"))<1
	empty = empty && ItemsInList(MacroList("*",";","WIN:Procedure"))<1
	empty = empty && ItemsInList(WaveList("*",";",""))<1

	String dirList = DataFolderDir(-1 )
	String Variables = StringByKey("VARIABLES",dirList)
	Variables = RemoveFromList("Gizmo_Error",Variables,",")
	String Strings = StringByKey("STRINGS",dirList)
	Strings = RemoveFromList("LogFileName",Strings,",")
	empty = empty && strlen(StringByKey("FOLDERS",dirList))<1
	empty = empty && strlen(StringByKey("WAVES",dirList))<1
	empty = empty && strlen(Variables)<1
	empty = empty && strlen(Strings)<1

	String wins=WinList("*",";","WIN:4311"), str=FunctionPath("LaueGo_Install")
	str = ParseFilePath(0,str,":",1,0)
	wins = RemoveFromList(str,wins)
	wins = RemoveFromList("Procedure",wins)
	empty = empty && ItemsInList(wins)<1

	String MainProcedure = ProcedureText("",0,"Procedure")	// contents of main Procedure window
	empty = empty && strsearch(MainProcedure,"#include",0)<0

	return empty
End



Function/T FetchLaueGoArchive(urlStr)			// download and expand new LaueGo zip file to the Desktop
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

	if (strlen(	urlStr)<1)
		Prompt urlStr, "URL to use for Downloading LaueGo", popup, "APS;gitHub"
		DoPrompt "Download Site",urlStr
		if (V_flag)
			return ""
		endif
		urlStr = SelectString(StringMatch(urlStr,"gitHub"),urlStr,gitHubArchiveURL)
		urlStr = SelectString(StringMatch(urlStr,"APS"),urlStr,JZTArchiveURL)
	endif
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

	String cmd=""
	sprintf str, "set zipPath to quoted form of POSIX path of (\"%s\" as alias)\n",destinationZip
	cmd += str
	sprintf str, "set outPath to quoted form of POSIX path of (\"%s\" as alias)\n",DesktopFldr
	cmd += str
	cmd += "tell application \"System Events\"\n"
	cmd += "do shell script \"unzip \" & zipPath & \" -d \" & outPath\n"
	cmd += "end tell\n"
	sprintf str, "set LGfolder to \"%s\" as alias\n",destinationFldr
	cmd += str
	sprintf str, "	tell application \"Finder\" to set label index of LGfolder to 6"	// 6=Green
	cmd += str
	//	print ReplaceString("\n",cmd,"\r")
	ExecuteScriptText/Z cmd
	if (V_flag)
		sprintf str, "failure in unzipping \"%s\"\rcmd = \r'%s'\r\r",destinationZip,ReplaceString("\n",cmd,"\r")
		Append2Log(str,0)
		Append2Log(S_value+"\r\r",0)				// there is no valid file path
	endif
	DeleteFile /Z destinationZip					// done unzipping, delete the zip file
	if (V_flag)
		sprintf str, "Could NOT delete the zip file  \"%s\"\r\tProceeding..."
		Append2Log(str,0)
	endif

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
// =================================== Start of xop aliases ==================================== //

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
	if (!FileFolderExists(SpecialDirPath("Igor Pro User Files",0,0,0)+"Igor Extensions",folder=1))
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

	String destination0 = SpecialDirPath("Igor Pro User Files",0,0,0)+"Igor Extensions"
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
//
Static Function AliasAlreadyThere(source,destination)		// returns 1 if the alias already exists
	String source
	String destination
	GetFileFolderInfo/Q/Z=1 destination
	if (V_isAliasShortcut)
		Variable i
		i = char2num(source[strlen(source)-1])
		source = SelectString(i==58, source, source[0,strlen(source)-2])
		i = char2num(S_aliasPath[strlen(S_aliasPath)-1])
		S_aliasPath = SelectString(i==58, S_aliasPath, S_aliasPath[0,strlen(S_aliasPath)-2])
		return StringMatch(source,S_aliasPath)
	endif
	return 0
End
//
Static Function/T FindFileInDirPath(start,name)		// returns full path searching from start
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
//
Static Function FileFolderExists(name,[path,file,folder])
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

// ==================================== End of xop aliases ===================================== //
// ============================================================================================= //

