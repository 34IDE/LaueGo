#pragma rtGlobals=1		// Use modern global access method.
#pragma IgorVersion = 6.11
#pragma version = 1.08
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
	if (strlen(inFile0)<1)
		return ""
	endif
	PathInfo home
	String path = SelectString(V_flag,"","home")

	String inFile = HFSToPosix(path,inFile0,1,1)
	if (strlen(inFile)<1)
		return ""
	endif

	String str, outFile
	sprintf str,"%.3f",DateTime
	outFile = SpecialDirPath("Temporary",0,1,1)+Hash(str,1)
	//	outFile = "/Users/tischler/Desktop/mda2Igor/test.txt"
	String pythonScript = ParseFilePath(1, FunctionPath("LoadMDAfile"), ":", 1, 0)+"mdaAscii2Igor.py"
	pythonScript = HFSToPosix("",pythonScript,1,1)
	if (strlen(pythonScript)<1)
		DoAlert 0,"Cannot find 'mdaAscii2Igor.py'"
		return ""
	endif

	String cmd
	sprintf cmd "do shell script \"\\\"%s\\\" \\\"%s\\\" \\\"%s\\\"\"",pythonScript,inFile,outFile
	ExecuteScriptText cmd
	Variable err=0
	err = strsearch(S_value,"Traceback (most recent call last)",0)>=0
	err = err || strsearch(S_value,"ERROR --",0)==0
	if (err)
		DoAlert 0, "failure in LoadMDAfile()"
		print "\r\r"
		print cmd
		print "\r\r"
		print S_value
		return ""								// there is no index file
	endif
	Variable i1,i2,i3
	i1 = strsearch(S_value,"\"",0)+1			// remove leading quote
	i2 = strsearch(S_value,"\n",0)-1
	i2 = i2<0 ? Inf : i2
	i3 = strsearch(S_value,"\"",i1)-1			// remove trailing quote
	i3 = i3<0 ? Inf : i3
	i2 = min(i2,i3)
	if (i1<0 || i2<0)
		return ""
	endif
	String mdaFile = PosixToHFS("",S_value[i1,i2],1,1)
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
	else
		DoAlert 0,"Only know how to display 1D and 2D mda waves"
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

	String win = StringFromList(0,FindGraphsWithWave(image))
	if (strlen(win))
		DoWindow/F $win			// bring window with image to the top
		return GetWavesDataFolder(image,2)
	endif

	wName = NameOfWave(image)
	String wNote=note(image)
	String pv=StringByKey("pv",wNote,"=")
	String xLabel=StringByKey("xLabel",wNote,"=")
	String yLabel=StringByKey("yLabel",wNote,"=")
	String unit
	unit = WaveUnits(image,0)
	xLabel += SelectString(strlen(xLabel)>0,"","  (\U)")
	unit = WaveUnits(image,1)
	yLabel += SelectString(strlen(yLabel)>0,"","  (\U)")

	Variable i = strsearch(wName,"_",0)
	String title=""
	if (i>=0)
		title = "mda file:\r3"+wName[i+1,Inf]
		title += SelectString(strlen(pv),"","\rPV = "+pv)
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

	String win = StringFromList(0,FindGraphsWithWave(line))
	if (strlen(win))
		DoWindow/F $win			// bring window with line to the top
		return GetWavesDataFolder(line,2)
	endif

	wName = NameOfWave(line)
	String wNote=note(line)
	String pv=StringByKey("pv",wNote,"=")
	String xLabel=StringByKey("xLabel",wNote,"="), yLabel=""
	String unit=WaveUnits(line,0)
	xLabel += SelectString(strlen(xLabel)>0,"","  (\U)")
	if (strlen(pv))
		yLabel = pv + SelectString(strlen(WaveUnits(line,1)),"","  (\U)")
	endif
	Variable i = strsearch(wName,"_",0)
	String title=""
	if (i>=0)
		title = "mda file:\r3"+wName[i+1,Inf]
		title += SelectString(strlen(pv),"","\rPV = "+pv)
	endif

	Display /W=(35,44,563,530) line
	ModifyGraph mirror=2, mirror=1,minor=1, lowTrip=0.001
	Label bottom xLabel
	Label left yLabel
	TextBox/C/N=text0/F=0/B=1/A=LT/X=3/Y=3 title
	return GetWavesDataFolder(line,2)
End
