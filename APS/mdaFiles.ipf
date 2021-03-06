#pragma rtGlobals=1		// Use modern global access method.
#pragma IgorVersion = 6.11
#pragma version = 1.15
#pragma ModuleName=mdaAPS

StrConstant mdaFilters = "Data Files (*.mda,*.MDA):.mda,.mda;All Files:.*;"

Menu "Load Waves"
	SubMenu "MDA"
		"Load MDA file...",LoadMDA("")
		"Re-Display MDA wave",DisplayMDAresult($"")
	End
End
//
Menu "New"
	"MDA Plot...",DisplayMDAresult($"")
End

#if stringmatch(IgorInfo(2),"Windows")
Function/T LoadMDA(inFile0)
	String inFile0
	DoAlert 0,"MDA files are loaded using a python script 'mdaAscii2Igor.py', which does not work on Windows.\rNothing loaded."
	return ""
End
//
Static Function/T LoadMDAfile(inFile0)
	String inFile0
	DoAlert 0,"MDA files are loaded using a python script 'mdaAscii2Igor.py', which does not work on Windows.\rNothing loaded."
	return ""
End
#else
Function/T LoadMDA(inFile0)
	String inFile0

	if (strlen(inFile0)<1)
		Variable f
		Open/D=1/F=mdaFilters/M="Pick an mda file"/R f
		inFile0 = S_fileName
		printf "LoadMDA(\"%s\")\r",inFile0
	endif
	if (strlen(inFile0)<1)
		return ""
	endif

	String fldrSav0= GetDataFolder(1)
	String fldrName = CleanupName(ParseFilePath(4,inFile0,":",0,0)+"_"+ParseFilePath(3,inFile0,":",0,0),0)
	NewDataFolder/O/S $fldrName
	String waveNames=LoadMDAfile(inFile0)
	Wave mda = $StringFromList(0,waveNames)
	if (!WaveExists(mda))
		SetDataFolder fldrSav0
		return ""
	endif
	Variable N=ItemsInList(waveNames)
	if (strlen(GetRTStackInfo(2))<1)			// plot if called from comand line
		printf "Loaded %d waves from mda file %s into folder '%s'\r",N,inFile0,fldrName
		DisplayMDAresult($"")
	endif
	SetDataFolder fldrSav0
	return GetDataFolder(2)
End
//
Static Function/T LoadMDAfile(inFile0)
	String inFile0

	if (strlen(inFile0)<1)
		Variable f
		Open/D=1/F=mdaFilters/M="Pick an mda file"/R f
		inFile0 = S_fileName
	endif
	String PosixFileIn=""
	if (isFile("home",inFile0))
		PathInfo home
		PosixFileIn = ParseFilePath(5,S_path+inFile0,"/",0,0)
	elseif (isFile("",inFile0))
		PosixFileIn = ParseFilePath(5,inFile0,"/",0,0)
	endif
	if (strlen(PosixFileIn)<1)
		return ""
	endif

	String str, outFile
	sprintf str,"%.3f",DateTime
	outFile = SpecialDirPath("Temporary",0,1,1)+Hash(str,1)
	//	outFile = "/Users/tischler/Desktop/mda2Igor/test.txt"
	String pythonScript = ParseFilePath(1, FunctionPath("LoadMDAfile"), ":", 1, 0)+"mdaAscii2Igor.py"
	pythonScript = ParseFilePath(5,pythonScript,"/",0,0)
	if (strlen(pythonScript)<1)
		DoAlert 0,"Cannot find 'mdaAscii2Igor.py'"
		return ""
	endif

	String cmd
	sprintf cmd "do shell script \"\\\"%s\\\" \\\"%s\\\" \\\"%s\\\"\"",pythonScript,PosixFileIn,outFile
#if (IgorVersion()>7)
	ExecuteScriptText/UNQ cmd
#else
	ExecuteScriptText cmd
#endif
	Variable err=0
	err = strsearch(S_value,"Traceback (most recent call last)",0)>=0
	err = err || strsearch(S_value,"ERROR --",0)==0
	if (err)
		DoAlert 0, "failure in LoadMDAfile()"
		print "\r\r"
		print cmd
		print "\r\r"
		print S_value
		return ""									// there is no index file
	endif
#if (IgorVersion()<7)
	S_value = TrimBoth(S_value,chars="\"")
#endif
	String mdaFile = Posix2HFS(TrimBoth(S_value))
	if (strlen(mdaFile)<1)
		return ""
	endif
	// printf "mdaFile =   %s\r",mdaFile

	LoadWave/T/M/Q mdaFile
	// V_flag			Number of waves loaded.
	// S_waveNames	Semicolon-separated list of the names of loaded waves.
	DeleteFile/Z mdaFile						// done with with this file

	String/G S_mdaWavesLoaded = S_waveNames
	String class
	Variable i,N=ItemsInList(S_mdaWavesLoaded)
	for (i=0;i<N;i+=1)
		Wave mda = $StringFromList(i,S_mdaWavesLoaded)
		if (WaveExists(mda))
			class = StringByKey("waveClass",note(mda),"=")
			if (strlen(class)<1)
				class = SelectString(WaveDims(mda)==2,"mda","mda,rawimage")
				if (WaveDims(mda)==1 && strsearch(NameOfWave(mda),"P",0)==0)
					class += "Positioner"
				endif
				Note/K mda, ReplaceStringByKey("waveClass",note(mda),class,"=")
			endif
		endif
	endfor

	if (strlen(GetRTStackInfo(2))<1)			// plot if called from comand line
		printf "Loaded %d waves from mda file %s\r",N,inFile0
		DisplayMDAresult($"")
	endif
	return S_waveNames
End

Static Function isFile(path,name)		// returns TRUE if path & name is a file
	String path			// optional path name, e.g. "home"
	String name			// partial or full file name

	if (strlen(path))
		PathInfo $path
		if (V_flag==0)
			return 0
		endif
		name = S_path+name	// add the path to name
	endif
	GetFileFolderInfo/Q/Z=1 name
	return V_isFile && (V_Flag==0)
End
#endif



Function/T DisplayMDAresult(mda)
	Wave mda

	String fldrSav0= GetDataFolder(1)
	NewDataFolder/O root:Packages
	NewDataFolder/O root:Packages:mda
	String defaultPVmda=StrVarOrDefault("root:Packages:mda:defaultPVmda",""), PVname
	String mdaFldr=StrVarOrDefault("root:Packages:mda:defaultFOLDERmda","")

	String mdaList = WaveListClass("mda","*","")
	if (!WaveExists(mda))
		if (ItemsInList(mdaList)<1)			// need to pick a sub-folder
			String fList0=DataFolderDir(1), fList="", fldr
			fList0 = StringByKey("FOLDERS",fList0)
			Variable i
			for (i=0;i<ItemsInList(fList0,",");i+=1)
				fldr = StringFromList(i,fList0,",")
				SetDataFolder $fldr
				if (ItemsInList(WaveListClass("mda","*",""))>0)
					fList += fldr+";"
				endif
				SetDataFolder fldrSav0
			endfor
			if (ItemsInList(fList)<1)
				return ""
			elseif (ItemsInList(fList)==1)
				mdaFldr = StringFromList(0,fList)
			else
				Prompt mdaFldr,"Folder with mda waves",popup,fList
				DoPrompt "mda wave",mdaFldr
				if (V_flag)
					return ""
				endif
			endif
			SetDataFolder $mdaFldr
			mdaList = WaveListClass("mda","*","")
			SetDataFolder fldrSav0
		else
			mdaFldr = ""
		endif

		String wName, item, dlist=""
		Variable k
		SetDataFolder $mdaFldr
		for (k=0;k<ItemsInList(mdaList);k+=1)
			Wave ww=$StringFromList(k,mdaList)
			if (WaveExists(ww))
				PVname = StringByKey("pv",note(ww),"=")
				PVname = SelectString(strlen(PVname),"unknown PV",PVname)
				item = StringFromList(k,mdaList) + " --> " + PVname
				wName = SelectString(stringmatch(PVname,defaultPVmda),wName,item)
				dlist += item+";"
			endif
		endfor
		SetDataFolder fldrSav0
		Prompt wName,"mda wave to display",popup,dlist
		DoPrompt "mda wave",wName
		if (V_flag)
			return ""
		endif
		k = strsearch(wName, " -->",0)
		wName = wName[0,k-1]
		wName = SelectString(strlen(mdaFldr),wName,":"+mdaFldr+":"+wName)
		Wave mda = $wName
	endif
	if (!WaveExists(mda))
		return ""
	endif

	String/G root:Packages:mda:defaultFOLDERmda=GetWavesDataFolder(mda,0)
	PVname = StringByKey("pv",note(mda),"=")
	PVname = SelectString(strlen(PVname),"unknown PV",PVname)
	String/G root:Packages:mda:defaultPVmda=PVname

	if (WaveDims(mda)==1)
		DisplayMDAresult1D(mda)
	elseif (WaveDims(mda)==2)
		DisplayMDAresult2D(mda)
	elseif (WaveDims(mda)==3)
		DisplayMDAresult3D(mda)
	else
		DoAlert 0,"Only know how to display 1D, 2D, and 3D mda waves"
	endif
	return GetWavesDataFolder(mda,2)
End
//
Static Function/T DisplayMDAresult2D(image)
	Wave image

	if (!WaveExists(image))
		String wName, str
		String list = WaveListClass("mda*","*",""), dlist=""
		Variable k
		for (k=0;k<ItemsInList(list);k+=1)
			Wave ww=$StringFromList(k,list)
			str = ""
			if (WaveExists(ww))
				str = StringByKey("pv",note(ww),"=")
			endif
			if (strlen(str))
				dlist += StringFromList(k,list) + " --> " + str+";"
			else
				dlist = StringFromList(k,list) + " --> unknown PV;"
			endif
		endfor
		Prompt wName,"mda image to display",popup,dlist
		DoPrompt "mda image",wName
		if (V_flag)
			return ""
		endif
		k = strsearch(wName, " -->",0)
		wName = wName[0,k-1]
		Wave image = $wName
	endif
	if (!WaveExists(image))
		return ""
	endif

	String win = StringFromList(0,WindowsWithWave(image,1))
	if (strlen(win))
		DoWindow/F $win			// bring window with image to the top
		return GetWavesDataFolder(image,2)
	endif

	wName = NameOfWave(image)
	String wNote=note(image)
	String pv=StringByKey("pv",wNote,"=")
	String xLabel=StringByKey("xLabel",wNote,"=")
	String yLabel=StringByKey("yLabel",wNote,"=")
	String fileTime = ISOtime2niceStr(StringByKey("file_time",wNote,"="))
	xLabel += SelectString(strlen(WaveUnits(image,0))>0,"","  (\U)")
	yLabel += SelectString(strlen(WaveUnits(image,1))>0,"","  (\U)")

	Variable i = strsearch(wName,"_",0)
	String title=""
	if (i>=0)
		title = "\\Zr075file:\\M "+wName[i+1,Inf]
		title += SelectString(strlen(pv),"","\r\\Zr075PV:\\M "+pv)
		title += SelectString(strlen(fileTime),"","\r\\Zr075"+fileTime+"\\M")
	endif

	Display /W=(35,44,563,530)
	AppendImage image
	ModifyImage $wName ctab= {*,*,Rainbow256,0}
	ModifyGraph mirror=2, mirror=1,minor=1, lowTrip=0.001
	Label bottom xLabel
	Label left yLabel
	TextBox/C/N=text0/F=0/B=1/A=LT/X=3/Y=3 title
	return GetWavesDataFolder(image,2)
End
//
Static Function/T DisplayMDAresult1D(line)
	Wave line

	if (!WaveExists(line))
		String wName, str
		String list = WaveListClass("mda","*",""), dlist=""
		Variable k
		for (k=0;k<ItemsInList(list);k+=1)
			Wave ww=$StringFromList(k,list)
			str = ""
			if (WaveExists(ww))
				str = StringByKey("pv",note(ww),"=")
			endif
			if (strlen(str))
				dlist += StringFromList(k,list) + " --> " + str+";"
			else
				dlist = StringFromList(k,list) + " --> unknown PV;"
			endif
		endfor
		Prompt wName,"mda wave to display",popup,dlist
		DoPrompt "mda wave",wName
		if (V_flag)
			return ""
		endif
		k = strsearch(wName, " -->",0)
		wName = wName[0,k-1]
		Wave line = $wName
	endif
	if (!WaveExists(line))
		return ""
	endif

	String win = StringFromList(0,WindowsWithWave(line,1))
	if (strlen(win))
		DoWindow/F $win			// bring window with line to the top
		return GetWavesDataFolder(line,2)
	endif

	wName = NameOfWave(line)
	String wNote=note(line)
	String pv=StringByKey("pv",wNote,"=")
	String xLabel=StringByKey("xLabel",wNote,"="), yLabel=""
	String fileTime = ISOtime2niceStr(StringByKey("file_time",wNote,"="))
	xLabel += SelectString(strlen(WaveUnits(line,0))>0,"","  (\U)")
	if (strlen(pv))
		yLabel = pv + SelectString(strlen(WaveUnits(line,1)),"","  (\U)")
	endif
	Variable i = strsearch(wName,"_",0)
	String title=""
	if (i>=0)
		title = "\\Zr075file:\\M "+wName[i+1,Inf]
		title += SelectString(strlen(pv),"","\r\\Zr075PV:\\M "+pv)
		title += SelectString(strlen(fileTime),"","\r\\Zr075"+fileTime+"\\M")
	endif

	Display /W=(35,44,563,530) line
	ModifyGraph mirror=2, mirror=1,minor=1, lowTrip=0.001
	Label bottom xLabel
	Label left yLabel
	TextBox/C/N=text0/F=0/B=1/A=LT/X=3/Y=3 title
	return GetWavesDataFolder(line,2)
End


Static Function/T DisplayMDAresult3D(planes)
	Wave planes
	if (!WaveExists(planes))
		return ""
	endif
	Variable maxLayer = DimSize(planes,2) - 1
	String wName = NameOfWave(planes)
	if (maxLayer<1)
		printf "'%s' is only a 2D (or less) wave, don't use this routine\r", wName
	endif

	String win=StringFromList(0,WindowsWithWave(planes,1))
	if (strlen(win))
		String str
		sprintf str, "Graph of '%s' already exists.\rPut up another?", wName
		DoAlert 1, str
		if (V_flag == 2)			// NO clicked
			DoWindow/F $win
			return GetWavesDataFolder(planes,2)
		endif
	endif

	Variable layer=0
	String wnote=note(planes)
	String pv=StringByKey("pv",wNote,"=")
	String fileTime = ISOtime2niceStr(StringByKey("file_time",wNote,"="))
	String xLabel=StringByKey("xAxisName",wNote,"=")
	String yLabel=StringByKey("yAxisName",wNote,"=")
	String zLabel = StringByKey("zAxisName",wNote,"=")
	zLabel = SelectString(CmpStr(zLabel, "dummy Z"),"","  "+zLabel) + zPositionStr(planes,layer)
	xLabel += SelectString(strlen(xLabel)>0,"","  (\U)")
	yLabel += SelectString(strlen(yLabel)>0,"","  (\U)")

	String title = "\Zr120layer "+num2istr(layer) + zLabel
	title += "\r\\Zr075file:\\M "+wName
	title += SelectString(strlen(pv),"","\r\\Zr075PV:\\M "+pv)
	title += SelectString(strlen(fileTime),"","\r\\Zr075"+fileTime+"\\M")

	Display/W=(133,45,704,587)/K=1
	AppendImage planes
	ModifyImage $wName ctab= {*,*,Rainbow256,0}
	ModifyImage $wName plane=layer			// always start with layer=0
	ModifyGraph margin(top)=30,width={Aspect,1}, tick=2, mirror=1, minor=1, lowTrip=0.001
	Label left yLabel
	Label bottom xLabel

	TextBox/C/N=title/F=0/B=(65535,65535,65535,25000)/X=3/Y=3/A=LT title

	SetVariable setvar_Layer,pos={25,5},size={68,19},bodyWidth=45,proc=mdaAPS#SetVarProcLayer,title="Layer"
	SetVariable setvar_Layer,fSize=13,format="%d",limits={0,maxLayer,1},value= _NUM:layer
	return GetWavesDataFolder(planes,2)
End
//
Static Function SetVarProcLayer(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	if (sva.eventCode == 1 || sva.eventCode == 2)		// mouse up or Enter key
		String win=sva.win
		Variable layer = round(sva.dval)
		String wName = StringFromList(0,ImageNameList(win,";"))
		Wave planes = ImageNameToWaveRef(win,wName)
		ModifyImage/W=$win $wName plane=layer

		String wnote=note(planes)
		String pv=StringByKey("pv",wnote,"=")
		String zLabel = StringByKey("zAxisName",wNote,"=")
		zLabel = SelectString(CmpStr(zLabel, "dummy Z"),"","  "+zLabel)

		String fileTime = ISOtime2niceStr(StringByKey("file_time",wnote,"="))
		String title = "\Zr120layer "+num2istr(layer) + zLabel + zPositionStr(planes,layer)
		title += "\r\\Zr075file:\\M "+wName
		title += SelectString(strlen(pv),"","\r\\Zr075PV:\\M "+pv)
		title += SelectString(strlen(fileTime),"","\r\\Zr075"+fileTime+"\\M")
		TextBox/C/N=title title
	endif
	return 0
End
//
Static Function/T zPositionStr(wav,layer)
	Wave wav
	Variable layer

	if (WaveDims(wav)<3)
		return ""
	endif
	Variable dz=DimDelta(wav,2), z0=DimOffset(wav,2)
	String unit = WaveUnits(wav,2), str
	unit = SelectString(strlen(unit), "", " ") + unit
	if (dz==1 && z0==0 && strlen(unit)<1)	// don't show position if un-scaled
		return ""
	endif
	sprintf str, " (= %g%s)", (dz*Layer+z0), unit
	return str
End

// moved to Utility_JZT.ipf
//
//Function/T ISOtime2niceStr(iso)
//	String iso
//	Variable year=NaN,month=NaN,day=NaN, hr=NaN,mn=NaN,se=NaN
//	sscanf iso,"%4d-%2d-%2dT%2d:%2d:%2d", year,month,day,hr,mn,se
//	Variable N = V_flag
//
//	if (N<3 || numtype(year+month+day))
//		return ""
//	endif
//	Variable epoch = date2secs(year, month, day )
//	String out = Secs2Date(epoch,2)
//	if (N>=5)
//		epoch += hr*3600
//		epoch += mn*60
//		se = (N>=6) ? se : 0
//		epoch += se
//		Variable fmt = (N>=6) ? 1 : 0
//		out += SelectString(numtype(epoch),"  "+Secs2Time(epoch,fmt),"")
//	endif
//	return out
//End


