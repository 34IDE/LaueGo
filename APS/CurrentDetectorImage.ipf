#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma IgorVersion = 6.3
#pragma version = 0.01
#pragma ModuleName=CurrentDetector
#include "ImageDisplayScaling"
#include "HDF5images"

Menu "Macros"
	"Show Current Detector", LoadCurrentDetector("")
End


Menu "Load Waves"
	"  Current Detector...", LoadCurrentDetector("")
End


// load an image from AreaDetector
Function/WAVE LoadCurrentDetector(detID,[printIt])
	String detID				// may be one of "Orange;Yellow;Purple" or the id such as "34idePE2:image1"
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : !(!printIt)

	if (strsearch(sytemHostname(), "xray.aps.anl.gov",0)<0)
		if (printIt)
			print "LoadCurrentDetector() only works at the APS on the beamline network"
		endif
		return $""
	endif

	String colors = "Orange=34idePE2:image1;Yellow=34idePE1:image1;Purple=34idePE3:image1;"
	String crates = "34idePE2:image1=Orange;34idePE1:image1=Yellow;34idePE3:image1=Purple;"

	String str = StringByKey(detID,colors,"=")
	if (strlen(str))
		detID = str
	endif

	if (strlen(detID)<1)
		Prompt detID, "detector", popup, "34idePE2:image1  Orange;34idePE1:image1  Yellow;34idePE3:image1  Purple;"
		DoPrompt "Detector", detID
		if (V_flag)
			return $""
		endif
		detID = StringFromList(0,detID," ")
	endif
	String color = StringByKey(detID,crates,"=")

	if (strlen(color)<1)
		return $""
	endif

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
	fileName = Posix2HFS(fileName)		// done with posix file paths, switch to Mac

	Variable err = strsearch(result,"Traceback (most recent call last)",0)>=0
	err = err || strsearch(result,"ERROR --",0,2)==0
	err = err || strsearch(result,"command not found",0,2)>=0
	if (err)
		DoAlert 0, "Failed to load latest image from "+detID
		print "\r\r  "+cmd+"\r\r  "+result
		return $""								// there is no image
	endif

	String fldrName = "root:Packages:micro:CurrentImage"
	String fldrSav = GetDataFolder(1)	// save current data folder
	if (DataFolderExists("root:Packages:micro:"))
		NewDataFolder/O/S $fldrName
	endif
	Wave image = $LoadGenericImageFile(fileName)
	if (!WaveExists(image))
		print "Unable to load image from disk"
		return $""
	endif

	SetDataFolder fldrSav					// and restore to original data folder
	if (printIt)
		printf "LoadCurrentDetector(\"%s\"",detID
		if (!ParamIsDefault(printIt))
			printf ", printIt=%g", printIt
		endif
		printf ")\t\t// Loaded iamge to %s\r",GetWavesDataFolder(image,2)
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
