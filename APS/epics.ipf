#pragma rtGlobals=1		// Use modern global access method.
#pragma IgorVersion = 6.0
#pragma version=1.7

Menu "Data"
	Submenu "EPICS"
	"spec staus for Inelastic",InelasticSpecStatusPanel()
	"Show a PV",ShowPV("")
	End
End


//	String/G root:Packages:EPICS:EPICS_BASE_PATH = getEpicsBasePath()		// this gets executed whenever the experiment restarts


//	##############################################################
//	####################### Start of Inelastic Status Panel ######################
Function InelasticSpecStatusPanel() : Panel
	if (strlen(WinList("SpecStatusPanel","","WIN:64")))
		DoWindow/F SpecStatusPanel
	endif
	String str = StringByKey("SCREEN1", IgorInfo(0))
	Variable i = strsearch(str,"RECT=",0)+5					//	RECT=left,top,right,bottom ;
	Variable right = str2num(StringFromList(2,str[i,Inf],","))-10
	NewPanel/K=1/W=(right-345-10,50,right-10,195)		//	/W=(left,top,right,bottom)
	DoWindow/C SpecStatusPanel

	TitleBox currentScan,pos={19,4},size={70,17},title="scan #NaN"
	TitleBox currentScan,font="Lucida Grande",fSize=12,frame=2,anchor= LT
	Button buttonScanning,pos={205,1},size={133,22},disable=2,title="off-line"
	Button buttonScanning,labelBack=(40000,40000,40000),fSize=14,fStyle=1
	TitleBox statusStr,pos={57,54},size={218,23},title="status"
	TitleBox statusStr,font="Lucida Grande",fSize=12,frame=2,anchor= MC
	TitleBox timeStr,pos={100,80},size={127,12},title="time stamp"
	TitleBox timeStr,font="Lucida Grande",fSize=10,frame=0,anchor= MC
	TitleBox userDirStr,pos={139,97},size={84,12},title="user dir"
	TitleBox userDirStr,font="Courier",fSize=10,frame=0
	TitleBox userFileStr,pos={5,97},size={66,12},title="file name"
	TitleBox userFileStr,font="Courier",fSize=10,frame=0
	TitleBox TitleStr,pos={104,31},size={115,21},title="title"
	TitleBox TitleStr,font="Lucida Grande",fSize=18,frame=0,fStyle=0,anchor= MT
	TitleBox Qvalue,pos={227,116},size={86,23},title="Q = NaN �\\S-1\\M"
	TitleBox Qvalue,font="Lucida Grande",fSize=12,frame=2,anchor= LB
	TitleBox hklValue,pos={8,116},size={193,23},title="hkl = (�x.xxx, �x.xxx, �x.xxx)"
	TitleBox hklValue,font="Lucida Grande",fSize=12,frame=2,anchor= LB

	SetWindow kwTopWin,hook(stopEpicsUpdates)=InelasticSpecStatusHook
	SetBackground InelasticPanelBackgroundTask()
	CtrlBackground start,period=15*60.15,dialogsOK=0,noBurst=1
End

Function InelasticSpecStatusHook(H_Struct)
	STRUCT WMWinHookStruct &H_Struct
	if (H_Struct.eventCode == 2)				// window was killed
		KillBackground
		return 1
	else
		return 0
	endif
End

Function InelasticPanelBackgroundTask()		// this is the background task that updates the panel
	if (!strlen(WinList("SpecStatusPanel","","WIN:64")))	// check if panel is up
		return 2
	endif
	String epicsOut = specInelasticEpicsStatus()	// first get the epics values
	if (strlen(epicsOut)==0)
		return 2								// error in background task
	endif
	Variable scanning,scanNum,H,K,L,Energy,wavelength,tth,Q_Angstrom
	String specStatus,specTime,specUserDir,specFileName,specTitle,specUser,specHeading
	scanning = NumberByKey("scanning",epicsOut,"=")
	scanning = numtype(scanning) ? -1 : scanning
	scanNum = NumberByKey("scanNum",epicsOut,"=")
	H = NumberByKey("H",epicsOut,"=")
	K = NumberByKey("K",epicsOut,"=")
	L = NumberByKey("L",epicsOut,"=")
	specStatus = StringByKey("specStatus",epicsOut,"=")
	specTime = StringByKey("specTime",epicsOut,"=")
	specUserDir = StringByKey("specUserDir",epicsOut,"=")
	specFileName = StringByKey("specFileName",epicsOut,"=")
	specTitle = StringByKey("specTitle",epicsOut,"=")
	specUser = StringByKey("specUser",epicsOut,"=")
	specHeading = StringByKey("specHeading",epicsOut,"=")
	Energy = NumberByKey("Energy",epicsOut,"=")
	wavelength = NumberByKey("wavelength",epicsOut,"=")
	tth = NumberByKey("tth",epicsOut,"=")
	Q_Angstrom = NumberByKey("Q_Angstrom",epicsOut,"=")

	Variable scanNumOld,scanningOld			// previous values used for the beep
	Variable i,needBeep=0
	String str
	ControlInfo/W=SpecStatusPanel currentScan
	i = strsearch(S_recreation,",title=\"scan #",0)
	scanNumOld = str2num(S_recreation[i+14,Inf])
	needBeep = (scanNumOld != scanNum) & (numtype(scanNumOld)==0 || numtype(scanNum)==0)
	ControlInfo/W=SpecStatusPanel buttonScanning
	i = strsearch(S_recreation,",title=\"",0)
	S_recreation = S_recreation[i+8,Inf]
	i = strsearch(S_recreation,"\"",0)
	S_recreation = S_recreation[0,i-1]
	strswitch(S_recreation)
		case "off-line":
			scanningOld = -1
			break
		case "scanning":
			scanningOld = 1
			break
		case "stopped":
			scanningOld = 0
			break
	endswitch
	needBeep += (scanningOld!= scanning && scanning!=1)
	if (needBeep)
		beep
	endif

	if (scanning<0)
		Button buttonScanning,win=SpecStatusPanel,title="off-line", labelBack=(40000,40000,40000)
	elseif (scanning==1)
		Button buttonScanning,win=SpecStatusPanel,title="scanning", labelBack=(65534,20000,20000)
	else
		Button buttonScanning,win=SpecStatusPanel,title="stopped", labelBack=(40000,40000,65534)
	endif
	TitleBox statusStr,win=SpecStatusPanel,title=specStatus
	TitleBox timeStr,win=SpecStatusPanel,title=specTime
	TitleBox userDirStr,win=SpecStatusPanel,title=specUserDir
	TitleBox userFileStr,win=SpecStatusPanel,title=specFileName
	TitleBox TitleStr,win=SpecStatusPanel,title=specTitle
	sprintf str,"Q = %.3f �\\S-1\\M",Q_Angstrom
	TitleBox Qvalue,win=SpecStatusPanel,title=str
	sprintf str,"scan #%d",scanNum
	TitleBox currentScan,win=SpecStatusPanel,title=str
	sprintf str,"hkl = (%.3f, %.3f, %.3f)",H,K,L
	TitleBox hklValue,win=SpecStatusPanel,title=str
	return 0					// no error
End
Function/S specInelasticEpicsStatus()
	String pvList="jfk:string11.DESC;jfk:string11;jfk:string12.DESC;jfk:string12;jfk:string13.DESC;jfk:string13;jfk:string14.DESC;jfk:string14;jfk:string15.DESC;jfk:string15;jfk:string16.DESC;jfk:string16;jfk:string17.DESC;jfk:string17;jfk:long11.DESC;jfk:long11;jfk:bit11;jfk:float1.DESC;jfk:float2.DESC;jfk:float3.DESC;jfk:float1;jfk:float2;jfk:float3;"
	pvList += "iad:BraggERdbkAO;iad:BraggLambdaRdbkAO;nwp:m4k:c0:m5;"
	String epicsOut = get_mult_PV(pvList)

	if (!stringmatch(StringByKey("jfk:float1.DESC",epicsOut,"="),"spec H"))
		DoAlert 0,"PV mismatch for 'spec H'"
		return ""
	elseif (!stringmatch(StringByKey("jfk:float2.DESC",epicsOut,"="),"spec K"))
		DoAlert 0,"PV mismatch for 'spec K'"
		return ""
	elseif (!stringmatch(StringByKey("jfk:float3.DESC",epicsOut,"="),"spec L"))
		DoAlert 0,"PV mismatch for 'spec L'"
		return ""
	elseif (!stringmatch(StringByKey("jfk:string11.DESC",epicsOut,"="),"spec status"))
		DoAlert 0,"PV mismatch for 'spec status'"
		return ""
	elseif (!stringmatch(StringByKey("jfk:string12.DESC",epicsOut,"="),"spec time stamp"))
		DoAlert 0,"PV mismatch for 'spec time stamp'"
		return ""
	elseif (!stringmatch(StringByKey("jfk:string13.DESC",epicsOut,"="),"spec userdir"))
		DoAlert 0,"PV mismatch for 'spec userdir'"
		return ""
	elseif (!stringmatch(StringByKey("jfk:string14.DESC",epicsOut,"="),"spec file"))
		DoAlert 0,"PV mismatch for 'spec file'"
		return ""
	elseif (!stringmatch(StringByKey("jfk:string15.DESC",epicsOut,"="),"spec title"))
		DoAlert 0,"PV mismatch for 'spec title'"
		return ""
	elseif (!stringmatch(StringByKey("jfk:string16.DESC",epicsOut,"="),"spec user"))
		DoAlert 0,"PV mismatch for 'spec user'"
		return ""
	elseif (!stringmatch(StringByKey("jfk:string17.DESC",epicsOut,"="),"spec heading"))
		DoAlert 0,"PV mismatch for 'spec heading'"
		return ""
	elseif (!stringmatch(StringByKey("jfk:long11.DESC",epicsOut,"="),"spec scan number"))
		DoAlert 0,"PV mismatch for 'spec scan number'"
		return ""
	endif

	Variable scanning = NumberByKey("jfk:bit11",epicsOut,"=")
	Variable scanNum = NumberByKey("jfk:long11",epicsOut,"=")
	Variable H=NumberByKey("jfk:float1",epicsOut,"=")
	Variable K=NumberByKey("jfk:float2",epicsOut,"=")
	Variable L=NumberByKey("jfk:float3",epicsOut,"=")
	String specStatus = StringByKey("jfk:string11",epicsOut,"=")
	String specTime = StringByKey("jfk:string12",epicsOut,"=")
	String specUserDir = StringByKey("jfk:string13",epicsOut,"=")
	String specFileName = StringByKey("jfk:string14",epicsOut,"=")
	String specTitle = StringByKey("jfk:string15",epicsOut,"=")
	String specUser = StringByKey("jfk:string16",epicsOut,"=")
	String specHeading = StringByKey("jfk:string17",epicsOut,"=")
	Variable Energy = NumberByKey("iad:BraggERdbkAO",epicsOut,"=")
	Variable wavelength = NumberByKey("iad:BraggLambdaRdbkAO",epicsOut,"=")
	Variable tth = NumberByKey("nwp:m4k:c0:m5",epicsOut,"=")
	Variable Q_Angstrom = 4*PI*sin(tth/2*PI/180)/wavelength

	String result = ""
	result = ReplaceNumberByKey("scanNum",result,scanNum,"=")
	result = ReplaceNumberByKey("scanning",result,scanning,"=")
	result = ReplaceStringByKey("specStatus",result,specStatus,"=")
	result = ReplaceStringByKey("specTime",result,specTime,"=")
	result = ReplaceStringByKey("specUserDir",result,specUserDir,"=")
	result = ReplaceStringByKey("specFileName",result,specFileName,"=")
	result = ReplaceStringByKey("specTitle",result,specTitle,"=")
	result = ReplaceStringByKey("specUser",result,specUser,"=")
	result = ReplaceStringByKey("specHeading",result,specHeading,"=")
	result = ReplaceNumberByKey("Energy",result,Energy,"=")
	result = ReplaceNumberByKey("wavelength",result,wavelength,"=")
	result = ReplaceNumberByKey("tth",result,tth,"=")
	result = ReplaceNumberByKey("Q_Angstrom",result,Q_Angstrom,"=")
	result = ReplaceNumberByKey("H",result,H,"=")
	result = ReplaceNumberByKey("K",result,K,"=")
	result = ReplaceNumberByKey("L",result,L,"=")

	if (ItemsInList(GetRTStackInfo(0))<2)
		printf "%s\r",specTitle
		printf "status is '%s'     at     %s\r",specStatus,specTime
		if (scanning)
			printf "\t\tspec is currently doing scan %d '%s'\r",scanNum,specHeading
		else
			printf "\t\tspec is not scanning,  scan number is %d\r",scanNum
		endif
		printf "\t\tuserdir = %s,   and data file = %s,   by user '%s'\r",specUserDir,specFileName,specUser
		printf "\t\t(hkl) = (%.3f  %.3f  %.3f)  at  Q=%.3f (1/�)\r",H,K,L,Q_Angstrom
	endif
	return result
End
//	######################## End of Inelastic Status Panel ######################
//	##############################################################


//	##############################################################
//	######################### Generic EPICS Functions #######################
Function/T ShowPV(pv)
	String pv
	if (strlen(pv)<1)
		Prompt pv,"name of the PV"
		DoPrompt "PV?",pv
		if (V_flag)
			return ""
		endif
	endif
	String val = get_PV_str(pv)
	if (strlen(val)<90)
		String str
		sprintf str, "'%s'  -->  '%s'\r",pv,val
		print str
		DoAlert 0, str
	else
		printf "'%s'  -->\r%s\r",pv,val
	endif
	return val
End

Function/T get_mult_PV(pvList)	// return a keyed list of values, one for each PV, e.g. "key1=value2;key2=value2;"
	String pvList							// semi-colon separated list of PVs
	if (strlen(pvList)<1)
		return ""
	endif
	if (!atAPS())
		DoAlert 0, "cannot make EPICS calls unless you are at the APS"
		return ""
	endif

//	String cmd = "do shell script \""+EPICS_BASE+"caget -t -n"
//	SVAR EPICS_BASE_PATH=root:Packages:EPICS:EPICS_BASE_PATH
	String cmd = "do shell script \""+getEpicsBasePath()+"caget -t -n -g 15"
	Variable j,i=0
	String pv=StringFromList(i,pvList)
	do
		cmd += " "+pv
		i += 1
		pv = StringFromList(i,pvList)
	while(strlen(pv))
	cmd += "\""							// add trailing double quote
	ExecuteScriptText cmd
	S_value = S_value[1,Inf]				// remove leading double quote
	i = strlen(S_value)
	S_value = S_value[0,i-2]				// remove trailing double quote

	String value,result = ""
	i = 0
	pv=StringFromList(i,pvList)
	do
		value = StringFromList(i,S_value,num2char(13))
		j = strlen(value)
		if (char2num(value[j-1])==32)	// remove the trailing space
			value = value[0,j-2]
		endif
		if (stringmatch(value," "+pv+" --- Invalid channel name"))
			value = "Invalid channel name"
		endif
		result = ReplaceStringByKey(pv,result,value,"=")
		i += 1
		pv = StringFromList(i,pvList)
	while(strlen(pv))
	return result
End

Function get_PV_num(pv)
	String pv
	if (!atAPS())
		DoAlert 0, "cannot make EPICS calls unless you are at the APS"
		return NaN
	endif
	String cmd
//	sprintf cmd "do shell script \"/Users/tischler/dev/EPICS/extensions/bin/darwin-ppc/caget -t -n %s\"", pv
//	SVAR EPICS_BASE_PATH=root:Packages:EPICS:EPICS_BASE_PATH
	sprintf cmd "do shell script \"%s/caget -t -n -g 15 %s\"", getEpicsBasePath(),pv
	ExecuteScriptText cmd
	return str2num(S_value[1,100])
End

Function/T get_PV_str(pv)
	String pv
	if (!atAPS())
		DoAlert 0, "cannot make EPICS calls unless you are at the APS"
		return ""
	endif
	String cmd
//	sprintf cmd "do shell script \"/Users/tischler/dev/EPICS/extensions/bin/darwin-ppc/caget -t %s\"", pv
//	SVAR EPICS_BASE_PATH=root:Packages:EPICS:EPICS_BASE_PATH
	sprintf cmd "do shell script \"%s/caget -t %s\"", getEpicsBasePath(),pv
	ExecuteScriptText cmd
	Variable i1,i2
	i1 = strsearch(S_value, "\"", 0)+1
	i2 = strsearch(S_value, "\"", Inf,1)-1
	return S_value[i1,i2]
End

Function/T get_PV_wave(pv,n)		// returns the name of the wave created
	String pv
	Variable n						// number of values in wave to get, �0 means get all
	if (!atAPS())
		DoAlert 0, "cannot make EPICS calls unless you are at the APS"
		return ""
	endif
	n = (numtype(n) || n<1) ? 0 : round(n)
//	SVAR EPICS_BASE_PATH=root:Packages:EPICS:EPICS_BASE_PATH
	String cmd
	if (n>0)
//		sprintf cmd "do shell script \"/Users/tischler/dev/EPICS/extensions/bin/darwin-ppc/caget -t -#%d %s\"", n,pv
		sprintf cmd "do shell script \"%s/caget -t -#%d -g 15 %s\"", getEpicsBasePath(),n,pv
	else
//		sprintf cmd "do shell script \"/Users/tischler/dev/EPICS/extensions/bin/darwin-ppc/caget -t %s\"", pv
		sprintf cmd "do shell script \"%s/caget -t -g 15 %s\"", getEpicsBasePath(),pv
	endif
	ExecuteScriptText cmd
	Variable i1,i2,i
	i1 = strsearch(S_value, "\"", 0)+1
	i2 = strsearch(S_value, "\"", Inf,1)-1
	for (;char2num(S_value[i1])<=32 && i1<i2;i1+=1)	// trim leading white space
	endfor
	String list = S_value[i1,i2]
	if (ItemsInList(list," ")<2)
		DoAlert 0, pv+" is not a wave form"
		return ""
	endif
	n = str2num(StringFromList(0, list," "))				// number of points in the wave

	String wName = pv										// name of wave to create
	i = strsearch(pv,".VAL",0)
	if ((strlen(pv)-i)==4)								// ends with ".VAL"
		wName = pv[0,i-1]									// remove the trailing ".VAL"
	endif
	wName = CleanupName(wName,0)						// make wave name based upon PV
	Make/N=(n)/O $wName
	Wave ww=$wName
	for (i=0;i<n;i+=1)
		ww[i] = str2num(StringFromList(i+1, list," "))
	endfor
	return wName
End


Function/T EPICS_put_PV_num(pv,value)
	String pv							// full PV name
	Variable value						// new number to set
	if (!atAPS())
		DoAlert 0, "cannot make EPICS calls unless you are at the APS"
		return "EPICS not Available"
	endif
	String cmd
//	sprintf cmd "do shell script \"%scaput -t -n %s %.15g\"", EPICS_BASE,pv,value
//	SVAR EPICS_BASE_PATH=root:Packages:EPICS:EPICS_BASE_PATH
	sprintf cmd "do shell script \"%scaput -t -n %s %.15g\"", getEpicsBasePath(),pv,value
	ExecuteScriptText/Z cmd
	if (V_flag==0)
		return ""
	endif
	return S_value[0,strlen(S_value)-2]				// remove leading and trailing double quote
End
//
Function/T EPICS_put_PV_str(pv,str)
	String pv							// full PV name
	String str							// new string value to set
	if (!atAPS())
		DoAlert 0, "cannot make EPICS calls unless you are at the APS"
		return "EPICS not Available"
	endif
	str = "\\\""+str+"\\\""
	String cmd
//	sprintf cmd "do shell script \"%scaput -t -s %s %s\"", EPICS_BASE,pv,str
//	SVAR EPICS_BASE_PATH=root:Packages:EPICS:EPICS_BASE_PATH
	sprintf cmd "do shell script \"%scaput -t -s %s %s\"", getEpicsBasePath(),pv,str
	ExecuteScriptText/Z cmd
	if (V_flag==0)
		return ""
	endif
	return S_value[0,strlen(S_value)-2]				// remove leading and trailing double quote
End


Function atAPS()											// returns true if you are currently at the APS
	if (!stringmatch(StringByKey("OS",IgorInfo(3)),"Macintosh OS X"))
		DoAlert 0, "Only know how to get hostname from a Mac"
		return 0								// cannot get answer
	endif
	ExecuteScriptText "do shell script \"hostname\""						//returns something like: 	"tischler.uni.aps.anl.gov"
	String ipAddress = ReplaceString("\"",S_value,"")
	Variable aps = (strsearch(ipAddress, ".aps.anl.gov",0)>0)
	return aps
End

//
//
//	�print get_PV_num("yum:status")
//	  1
//
//
//	�print get_PV_num("Mt:S:NumBucketsFilledAI")
//	  24
//	�print get_PV_wave("Mt:S:FilledBucketsWF",24)
//	  Mt_S_FilledBucketsWF
//	�print Mt_S_FilledBucketsWF
//	  Mt_S_FilledBucketsWF[0]= {0,54,108,162,216,270,324,378,432,486,540,594,648,702,756,810,864,918,972,1026,1080,1134,1188,1242}
//
//
//Function get_PV(pv)
//	String pv
//	String cmd
//	sprintf cmd "do shell script \"/Users/tischler/dev/EPICS/extensions/bin/darwin-ppc/caget -t %s\"", pv
//	ExecuteScriptText cmd
//	return str2num(S_value[1,100])
//End
//
//
//
//Function  xxx()
//	String cmd = "do shell script \"/Users/tischler/dev/EPICS/extensions/bin/darwin-ppc/caget yum:status\""
//	//String cmd = "do shell script \"echo my-command | /bin/tcsh\""
//	ExecuteScriptText cmd
//	print S_value
//End


Function/T getEpicsBasePath()
	String env = getEnvironment()
	String base = StringByKey("EPICS_BASE", env,"=")
	String path = StringByKey("PATH", env,"=")
	String str
	Variable i
	for (i=0, str=StringFromList(i,path,":"); strlen(str); i+=1, str=StringFromList(i,path,":"))
		if (strsearch(str, base,0)>=0)
			return str+"/"
		endif
	endfor
	return ""
End


Function epicsInitPackage()
//	NewDataFolder/O root:Packages
//	NewDataFolder/O root:Packages:EPICS
//	String/G root:Packages:EPICS:EPICS_BASE_PATH = getEpicsBasePath()
	return 0
End