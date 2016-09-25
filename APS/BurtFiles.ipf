#pragma rtGlobals=1		// Use modern global access method.
#pragma IgorVersion = 6.20
#pragma version = 0.04
#pragma ModuleName=BurtFiles

Menu "Load Waves"
	SubMenu "Burt Files"
		"Load PV from Burt Files...",LoadPVfromBurt("","")
		"Re-Display PV from Burt",DisplayBurtResult($"")
	End
End




Function/T DisplayBurtResult(burt)
	Wave/T burt

	String burtName

	String burtList = WaveListClass("BurtPV","*","TEXT:1,DIMS:2")
	if (!WaveExists(wBurt))
		if (ItemsInList(burtList)<1)
			return ""
		elseif (ItemsInList(burtList)==1)
			Wave/T burt = $StringFromList(0,burtList)
		else
			Prompt burtName,"Burt wave to display",popup,burtList
			DoPrompt "Burt wave",burtName
			if (V_flag)
				return ""
			endif
			Wave/T burt = $burtName
		endif
	endif
	if (!WaveExists(burt))
		return ""
	endif
	Variable N=DimSize(burt,0)
	if (N<1)
		return ""
	endif
	Make/N=(N)/FREE/T burt1			// check if values are number (cannot plot strings)
	burt1 = burt[p][1]
	if (!isTextWaveNumbers(burt1))
		return ""
	endif

	String pvName = CleanupName(NameOfWave(burt)+"_PV",0)
	String epochName = CleanupName(NameOfWave(burt)+"_Epoch",0)
	Make/N=(N)/O/D $epochName/WAVE=epoch
	epoch = dateStr2epoch(burt[p][0])
	SetScale d 0,0,"dat", epoch
	Make/N=(N)/O/D $pvName/WAVE=pv
	pv = Text2Number(burt[p][1])
	String wnote=note(burt)
	Note/K epoch, AddClassToWaveNote(wnote,"BurtXvalues")
	Note/K pv, AddClassToWaveNote(wnote,"BurtYvalues")

	String gName=StringFromLIst(0,WindowsWithWave(pv,1))
	if (strlen(gName)>0)
		DoWindow/F $gName
		return gName
	endif

	pvName=StringByKey("PVname",wnote,"=")
	Display pv vs epoch
	ModifyGraph tick=2, mirror=1, minor=1, lowTrip=0.001
	ModifyGraph dateInfo(bottom)={0,0,0}
	Label bottom "time"
	if (strlen(pvName))
		Label left PVname
	endif
	return StringFromList(0,WinList("*", ";", "WIN:1"))
End




#if stringmatch(IgorInfo(2),"Windows")
Function/WAVE LoadPVfromBurt(PVname,start,[period,final])
	String PVname
	String start						// first date (YYYY-MM-DD)
	Variable period					// sampling period (hours), cannot be less than 0.25hr
	String final						// final date (YYYY-MM-DD), defaults to start, if final is "now" or "today" then today's date is used
	DoAlert 0,"Burt files are loaded using a python script 'burt_tool.py', which does not work on Windows.\rNothing loaded."
	return $""
End
//
Static Function/WAVE LoadPVfrom34IDBurtFiles(PVname,start,[period,final])
	String PVname
	String start						// first date (YYYY-MM-DD)
	Variable period					// sampling period (hours), cannot be less than 0.25hr
	String final						// final date (YYYY-MM-DD), defaults to start, if final is "now" or "today" then today's date is used
	DoAlert 0,"Burt files are loaded using a python script 'burt_tool.py', which does not work on Windows.\rNothing loaded."
	return $""
End


#else
Function/WAVE LoadPVfromBurt(PVname,start,[period,final])
	String PVname
	String start						// first date (YYYY-MM-DD)
	Variable period					// sampling period (hours), cannot be less than 0.25hr
	String final						// final date (YYYY-MM-DD), defaults to start, if final is "now" or "today" then today's date is used

	if (strlen(PVname)<1 || strlen(start)<1)
		Prompt PVname,"Name of epics PV"
		Prompt start,"start date in  year-month-day format, e.g. '2012-11-03'"
		Prompt final,"final date, this can also be blank or \"today\" or \"now\""
		Prompt period,"sampling period (hours), cannot be less the 0.25 (15 minutes)"

		period = ParamIsDefault(period) ? 1 : period
		if (strlen(start)<1)
			start = Secs2Date(DateTime,-2)
		endif
		DoPrompt "PV from Burt" PVname,start,final,period
		if(V_flag)
			return $""
		endif
	endif

	if (strlen(PVname)<1 || strlen(start)<1)
		return $""
	endif
	Variable gotFinal = strlen(final)>0
	Variable gotPeriod = gotFinal || !(period==1)
	if (gotFinal)
		Wave/T burt=LoadPVfrom34IDBurtFiles(PVname,start,period=period,final=final)
	elseif(gotPeriod)
		Wave/T burt=LoadPVfrom34IDBurtFiles(PVname,start,period=period)
	else
		Wave/T burt=LoadPVfrom34IDBurtFiles(PVname,start)
	endif
	if (WaveExists(burt))
		DisplayBurtResult(burt)
	endif
	return burt
End
//
Static Function/WAVE LoadPVfrom34IDBurtFiles(PVname,start,[period,final])
	String PVname
	String start						// first date (YYYY-MM-DD)
	Variable period					// sampling period (hours), cannot be less than 0.25hr
	String final						// final date (YYYY-MM-DD), defaults to start, if final is "now" or "today" then today's date is used
	period = ParamIsDefault(period) ? NaN : period
	period = limit(period,0.25,1000)
	final = SelectString(ParamIsDefault(final),final,"")

	if (strlen(PVname)<1)
		return $""
	endif

	Variable year,month,day, hour,minute,second
	sscanf start,"%4d-%02d-%02d",year,month,day
	if (V_flag!=3)
		return $""
	endif
	if (!(year>1990 && year<2050) && (month>=1 && month<=12) && (day>=1 && day<=31))
		printf "start date = '%s' is INVALID, must be of form \"2012-09-20\"\r",start
		return $""
	endif
	sprintf start,"%4d-%02d-%02d",year,month,day

	if (stringmatch(final,"now") || stringmatch(final,"today"))
		final = Secs2Date(DateTime,-2)
	elseif (strlen(final))		// final was passed, check it
		sscanf final,"%4d-%02d-%02d",year,month,day
		if (V_flag!=3)
			return $""
		endif
		if (!(year>1990 && year<2050) && (month>=1 && month<=12) && (day>=1 && day<=31))
			printf "final date = '%s' is INVALID, must be of form \"2012-09-20\"\r",final
			return $""
		endif
		sprintf final,"%4d-%02d-%02d",year,month,day
		period = numtype(period) ? 1 : period
	endif
	if (!(ISOtime2IgorEpoch(final)>ISOtime2IgorEpoch(start)))
		DoAlert 0, "final date is before start date"
		return $""
	endif

	String pythonScript = ParseFilePath(1, FunctionPath("BurtFiles#LoadPVfrom34IDBurtFiles"), ":", 1, 0)+"burt_tool.py"
	pythonScript = ParseFilePath(5,pythonScript,"/",0,0)
	if (strlen(pythonScript)<1)
		DoAlert 0,"Cannot find 'burt_tool.py'"
		return $""
	endif
	pythonScript = "\\\""+pythonScript+"\\\""		// the path to burt_tool.py may have spaces, so it needs to be quoted
	pythonScript += " "+PVname+" "+start			// ./burt_tool.py 34ide:userTran1.J 2012-10-20 1 2012-11-03 > x.txt

	if (numtype(period)==0)
		pythonScript += " "+num2str(period)
	endif
	if (strlen(final))
		pythonScript += " "+final
	endif
	ExecuteScriptText  "do shell script \""+pythonScript+"\""

	String buf = S_value
	buf = ReplaceString("\n", buf, "\r")
	buf = ReplaceString("\r\r", buf, "\r")
	Variable i=strlen(buf)
	buf = buf[1,i-2]				// remove bounding double-quotes
	if (strlen(buf)<1 || stringmatch(buf,"No such file or directory"))
		print "error"
		return $""
	endif

	String line, str
	Variable Nmax=1000, N=0
	Make/N=(Nmax,2)/T/O $CleanupName(PVname,0)/WAVE=pv
	SetDimLabel 1,0,DateTimeStr,pv
	SetDimLabel 1,1,pvValue,pv
	Variable i0,i1
	for (i0=0,i1=0,N=0; i1>=0;)
		i1 = strsearch(buf,"\r",i0)
		line = buf[i0,i1-1]
		if (N>=Nmax)
			Nmax += 1000
			Redimension/N=(Nmax,-1) pv
		endif
		sscanf StringFromList(1,line,"\t"),"%3s %3s %d %02d:%02d:%02d %4d",str,str,day,hour,minute,second,year
		if (V_flag!=7)
			continue
		endif
		month = WhichListItem(str,"Jan;Feb;Mar;Apr;May;Jun;Jul;Aug;Sep;Oct;Nov;Dec")+1
		sprintf str,"%4d-%02d-%02dT%02d:%02d:%02d",year,month,day,hour,minute,second
		pv[N][0] = str
		pv[N][1] = StringFromList(3,line,"\t")
		N += 1
		i0 = i1 + 1
	endfor
	Redimension/N=(N,-1) pv

	String wnote="waveClass=BurtPV"
	wnote = ReplaceStringByKey("PVname",wnote,PVname,"=")
	wnote = ReplaceStringByKey("startDate",wnote,start,"=")

	if (numtype(period)==0 && period != 0)
		wnote = ReplaceNumberByKey("period",wnote,period,"=")
	endif
	if (strlen(final)>0)
		wnote = ReplaceStringByKey("finalDate",wnote,final,"=")
	endif
	Note/K pv, wnote
	return pv
End
#endif



Static Function dateStr2epoch(str)				// convert date-time string of form "2012-11-03T00:12:00" to the Igor epoch
	String str
	Variable year,month,day,hour,minute,second, epoch
	sscanf str,"%04d-%02d-%02dT%02d:%02d:%02d", year,month,day,hour,minute,second
	epoch = (V_flag==6) ? date2secs(year,month,day) + hour*3600 + minute*60 + second : NaN
	return epoch
End
//
Static Function isTextWaveNumbers(wStr)
	Wave/T wStr
	Make/N=(numpnts(wStr))/FREE wFloat
	wFloat = Text2Number(wStr)
	return numtype(sum(wFloat))==0
End
//	43,45,46,48-57
//	43 is "
//	45 is -
//	46 is period
//	48-57 is 0-9
Static Function Text2Number(strIN,[int])		// convert str to a number only if str is ONLY a number
	String strIN
	Variable int									// only integers allowed
	int = ParamIsDefault(int) ? 0 : int
	int = numtype(int) ? 0 : !(!int)

	String str = TrimFrontBackWhiteSpace(strIN)
	Variable DP = int ? NaN : 46				// I use Inf as an invalid character

	if (!int)
		Variable i = strsearch(str,".",0)		// only one decimal point allowed
		if (i>=0)
			if (strsearch(str,".",i+1)>=0)		// found two decimal points
				return NaN
			endif
		endif
	endif
	i=char2num(str[0])
	if (! (i==43 || i==45 || i==DP || (i>=48 && i<=57)))	// first char must be in {+-.0123456789}
		return NaN
	elseif (i<48)
		str = str[1,Inf]							// trim off first char if not a digit
	endif

	// check for exponenet
	str = LowerStr(str)
	i = strsearch(str,"e",0)
	if (strsearch(str,"e",i+1)>=0)				// found two exponentials
		return NaN
	endif

	if (i>=0)
		if (numtype(Text2Number(str[i+1,Inf],int=1)))
			return NaN
		endif
		str = str[0,i-1]
	endif

	// remaining chars can only be {.0123456789}
	Make/N=(strlen(str))/FREE num
	num = char2num(str[p])
	num = num==DP ? 0 : num
	num = (num[p]>=48 && num[p]<=57) ? 0 : num[p]
	if (sum(num))
		return NaN
	endif
	return str2num(strIN)
End
