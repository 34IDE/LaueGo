#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma IgorVersion = 6.3
#pragma version = 0.05
#pragma ModuleName=CurrentDetector
#include "ImageDisplayScaling"
#include "HDF5images"


// REQUIRED_SUBNET and KNOWN_ADs can be overridden in your Main Procedure window
StrConstant REQUIRED_SUBNET = "xray.aps.anl.gov"
StrConstant KNOWN_ADs = "34idePE2:image1=Orange;34idePE1:image1=Yellow;34idePE3:image1=Purple;"


Menu "Macros"
	"Show Current Detector", LoadCurrentDetector("")
End


Menu "Load Waves"
	"  Current Detector...", LoadCurrentDetector("")
End


// load an image from AreaDetector
Function/WAVE LoadCurrentDetector(detID,[printIt])
	String detID				// may be one of the colors in KNOWN_ADs, e.g. "Orange" or the id such as "34idePE2:image1"
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : !(!printIt)

	String str
	if (strsearch(sytemHostname(), REQUIRED_SUBNET,0)<0)
		if (printIt)
			sprintf str, "LoadCurrentDetector() only works on the \"%s\" beamline network", REQUIRED_SUBNET
			DoAlert 1, str
			print str
		endif
		return $""
	endif

	detID = RemoveEnding(detID,":")		// if the user supplied a trailing ":", remove it
	str=ad_from_color(detID)				// if detID is a color, then this convert to AreaDetector image PV
	detID = SelectString(strlen(str), detID, str)
	if (strlen(detID)<1)
		Prompt detID, "detector", popup, ReplaceString("=",KNOWN_ADs,"  ")
		DoPrompt "Detector", detID
		if (V_flag)
			return $""
		endif
		detID = StringFromList(0,detID," ")
	endif
	String color = StringByKey(detID,KNOWN_ADs,"=")
	color = SelectString(strlen(color),StringFromList(0,detID,":"),color)
	color = CleanupName(color,0)

	String pythonScript = FunctionPath("LoadCurrentDetector")
	pythonScript = ParseFilePath(1,pythonScript,":",1,0)+"getADimageH5.py"
	GetFileFolderInfo/Q/Z=1 pythonScript
	if (V_Flag || !V_isFile)
		if (printIt)
			printf "LoadCurrentDetector() unable to find the python script\r\t\t\"%s\"",pythonScript
		endif
		return $""
	endif

	String fileName=SpecialDirPath("Temporary",0,1,1)+color+"Current.h5", cmd, result=""
	Variable isMac = stringmatch(igorInfo(2),"Macintosh")
	Variable isWin = stringmatch(igorInfo(2),"Windows")
	if (isMac)
		String exportCmd = exportVariablesCommands("EPICS_BASE;EPICS_HOST_ARCH;PYEPICS_LIBCA;VERSIONER_PYTHON_PREFER_32_BIT;EPICS_CA_MAX_ARRAY_BYTES")
		pythonScript = ParseFilePath(5,pythonScript,"/",0,0)	// convert HFS to POSIX
		sprintf cmd "do shell script \"   %s ;  \\\"%s\\\" \\\"%s\\\" \\\"%s\\\"\"",exportCmd,pythonScript,detID,fileName
		ExecuteScriptText/Z cmd
		result = S_value
	else
		print "Windows not yet implemented"
		return $""
	endif
	Variable unAvailable = StringMatch(result,"cannot connect to *:PortName_RBV")
	Variable err = unAvailable || strsearch(result,"Traceback (most recent call last)",0)>=0
	err = err || strsearch(result,"ERROR --",0,2)==0
	err = err || strsearch(result,"command not found",0,2)>=0
	fileName = Posix2HFS(fileName)		// done with posix file paths, switch to Mac style
	GetFileFolderInfo/Q/Z=1 fileName	// check that hdf5 file exists
	err = err || (V_Flag || !V_isFile)
	printIt = printIt || err
	if (printIt)
		printf "LoadCurrentDetector(\"%s\")",detID
	endif
	if (err)
		printf "\r"
		str = "Failed to load latest image from "+detID
		str += SelectString(unAvailable, "", "\r  '"+detID+"' cannot be reached, is it on?")
		DoAlert 0, str
		if (unAvailable)
			printf "detID = \"%s\" cannot be reached, check if it is on?\r",detId
		else
			print "\r\r  "+cmd+"\r\r  "+result
		endif
		return $""								// there is no image
	endif

	String fldrName = "root:Packages:micro:CurrentImage"
	String fldrSav = GetDataFolder(1)	// save current data folder
	if (DataFolderExists("root:Packages:micro:"))
		NewDataFolder/O/S $fldrName
	endif
	Wave image = $LoadGenericImageFile(fileName)
	SetDataFolder fldrSav					// and restore to original data folder
	if (!WaveExists(image))
		printf "\r\t\tUnable to load image from disk"
		return $""
	elseif (printIt)
		printf "\t\t// Loaded image to %s\r",GetWavesDataFolder(image,2)
	endif

	String win = StringFromList(0,FindGraphsWithWave(image))
	if (strlen(win))
		DoWindow/F $win							// bring graph with image to front
	else
		FUNCREF NewImageGraph fnew= $"Indexing#NewImageGraphLocal"
		fnew(image)
	endif
	return image
End
//
Static Function/T ad_from_color(color)
	String color

	String ad, item
	Variable i
	do
		item = StringFromList(i,KNOWN_ADs)
		if (StringMatch(StringFromList(1,item,"="),color))
			return StringFromList(0,item,"=")
		endif
		i += 1
	while(strlen(item)>0)
	return ""
End
//
Static Function/T exportVariablesCommands(varList)	// make an export command with quoted values
	String varList						// list of current environment variables to get and format
	String env=getEnvironment()	// get key=value list of all environment variables
	String str, var, value, cmd=""
	Variable i, N=ItemsInList(varList)
	for (i=0;i<N;i+=1)
		var = StringFromList(i,varList)
		value = StringByKey(var,getEnvironment(),"=")
		if (strlen(var)>0 && strlen(value)>0)
			sprintf str "%s=\\\"%s\\\" ",var,value
			cmd += str
		endif
	endfor
	cmd = SelectString(strlen(cmd)>1, "", "export ")+cmd
	return cmd
End
//
//	PYEPICS_LIBCA=/usr/local/epics/base/lib/darwin-x86/libca.dylib
//	EPICS_BASE=/usr/local/epics/base
//	EPICS_HOST_ARCH=darwin-x86
//	VERSIONER_PYTHON_PREFER_32_BIT=NO
//	EPICS_CA_MAX_ARRAY_BYTES=33554432
//	PATH=/opt/local/bin:/opt/local/sbin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin:/usr/local/git/bin:/usr/texbin:/usr/local/epics/base/bin/darwin-x86:/usr/local/epics/extensions/bin/darwin-x86:/usr/local/texlive/2010/bin/universal-darwin:/Users/tischler/bin:/Users/tischler/bin/xop2.3
