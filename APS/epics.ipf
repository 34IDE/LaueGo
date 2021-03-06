#pragma rtGlobals=1		// Use modern global access method.
#pragma IgorVersion = 6.0
#pragma version=1.84

Menu "Data"
	Submenu "EPICS"
	"Show a PV",ShowPV("")
	End
End




//	#########################################################################
//	######################## Generic EPICS Functions ########################
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
	String pvList						// semi-colon separated list of PVs
	if (strlen(pvList)<1)
		return ""
	endif
	if (!atAPS())
		DoAlert 0, "cannot make EPICS calls unless you are at the APS"
		return ""
	endif

	String cmd = "do shell script \""+getEpicsBasePath()+"caget -t -n -g 15"
	Variable j,i=0
	String pv=StringFromList(i,pvList)
	do
		cmd += " "+pv
		i += 1
		pv = StringFromList(i,pvList)
	while(strlen(pv))
	cmd += "\""									// add trailing double quote
#if (IgorVersion()<7)
	ExecuteScriptText cmd
	S_value = TrimBoth(S_value,chars="\"")
#else
	ExecuteScriptText/UNQ cmd
#endif
	String out = S_value

	String value,result = ""
	i = 0
	pv=StringFromList(i,pvList)
	do
		value = StringFromList(i,out,num2char(13))
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
	sprintf cmd "do shell script \"%s/caget -t -n -g 15 %s\"", getEpicsBasePath(),pv
#if (IgorVersion()<7)
	ExecuteScriptText cmd
	S_value = S_value[1,50]
#else
	ExecuteScriptText/UNQ cmd
#endif
	return str2num(S_value)
End

Function/T get_PV_str(pv)
	String pv
	if (!atAPS())
		DoAlert 0, "cannot make EPICS calls unless you are at the APS"
		return ""
	endif
	String cmd
	sprintf cmd "do shell script \"%s/caget -t %s\"", getEpicsBasePath(),pv
#if (IgorVersion()<7)
	ExecuteScriptText cmd
	S_value = TrimBoth(S_value,chars="\"")
#else
	ExecuteScriptText/UNQ cmd
#endif
	return S_value
End

Function/T get_PV_wave(pv,n)		// returns the name of the wave created
	String pv
	Variable n						// number of values in wave to get, �0 means get all
	if (!atAPS())
		DoAlert 0, "cannot make EPICS calls unless you are at the APS"
		return ""
	endif
	n = (numtype(n) || n<1) ? 0 : round(n)
	String cmd
	if (n>0)
		sprintf cmd "do shell script \"%s/caget -t -#%d -g 15 %s\"", getEpicsBasePath(),n,pv
	else
		sprintf cmd "do shell script \"%s/caget -t -g 15 %s\"", getEpicsBasePath(),pv
	endif
#if (IgorVersion()<7)
	ExecuteScriptText cmd
	S_value = TrimBoth(S_value,chars="\"")
#else
	ExecuteScriptText/UNQ cmd
#endif
	String list = TrimBoth(S_value)
	if (ItemsInList(list," ")<2)
		DoAlert 0, pv+" is not a wave form"
		return ""
	endif
	n = str2num(StringFromList(0, list," "))				// number of points in the wave

	String wName = pv										// name of wave to create
	Variable i = strsearch(pv,".VAL",0)
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


Function/T EPICS_put_PV_num(pv,value,[fmt])
	String pv							// full PV name
	Variable value						// new number to set
	String fmt							// OPTIONAL format specifier for value, default is "%.15g"
	fmt = SelectString(ParamIsDefault(fmt) || strlen(fmt)<2, fmt, "%.15g")
	if (!atAPS())
		DoAlert 0, "cannot make EPICS calls unless you are at the APS"
		return "EPICS not Available"
	endif
	String cmd
//	sprintf cmd "do shell script \"%scaput -t -n %s %.15g\"", getEpicsBasePath(),pv,value
	sprintf cmd "do shell script \"%scaput -t -n %s "+fmt+"\"", getEpicsBasePath(),pv,value
#if (IgorVersion()<7)
	ExecuteScriptText/Z cmd
	S_value = TrimBoth(S_value,chars="\"")
#else
	ExecuteScriptText/Z/UNQ cmd
#endif
	if (V_flag==0)
		return ""
	endif
	return S_value
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
	sprintf cmd "do shell script \"%scaput -t -s %s %s\"", getEpicsBasePath(),pv,str
#if (IgorVersion()<7)
	ExecuteScriptText/Z cmd
	S_value = TrimBoth(S_value,chars="\"")
#else
	ExecuteScriptText/UNQ/Z cmd
#endif
	if (V_flag==0)
		return ""
	endif
	return S_value
End


Function atAPS()											// returns true if you are currently at the APS
	Variable aps = (strsearch(sytemHostname(), ".aps.anl.gov",0)>0)
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
	String base = getEnvironment("EPICS_BASE")
	String path = getEnvironment("PATH")
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