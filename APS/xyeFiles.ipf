#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName=loadXYE
#pragma IgorVersion = 6.3
#pragma version = 1.00
#include "Utility_JZT", version>=4.21


Static StrConstant xyeFilters = "xye Files (*.xye):.xye,;All Files:.*;"


Menu "Load Waves"
	"Load xye File...", Load_xye_file("","")
	MenuItemIfWaveClassExists("  Graph loaded xye data","xye*","DIMS:1"), Graph_xye($"")
End


Function/S Load_xye_file(fileName,path,[printIt])
	String fileName		// if none given, a dialog will come up.
	String path				// Optional, you can use ""
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRtStackInfo(2))==0 : printIt
	path = SelectString(strlen(path)<1, path, "")			// ensure an empty string, not a NULL string
//	path = SelectString(strlen(path)<1 && strlen(PathList("home","","")), path, "home")	// use home if no path given
 
 	Variable f
	Open/F=xyeFilters/P=$path/M="xye file"/R/Z=2 f as fileName
	String fullFileName = S_fileName				// save the full file name
	if (V_flag)
		return ""
	endif
	FStatus f
	String buf = PadString("",V_logEOF,0x20)
	FBinRead f, buf
	Close f
	buf = ReplaceString("\r",buf,"\n")
	do
		buf = ReplaceString("\n\n",buf,"\n")
	while (strsearch(buf,"\n\n",0)>=0)

	String wnote="waveClass=xye;"
	wnote = addString2WaveNote(buf,wnote,"SPECfile")
	wnote = addString2WaveNote(buf,wnote,"Project")
	wnote = addString2WaveNote(buf,wnote,"Scan_numbers")
	wnote = addString2WaveNote(buf,wnote,"Instrument_configuration")
	wnote = addString2WaveNote(buf,wnote,"Detector_configuration")
	wnote = addString2WaveNote(buf,wnote,"Bad_pixel_file")
	wnote = addNumber2WaveNote(buf,wnote,"Wavelength")
	wnote = addNumber2WaveNote(buf,wnote,"Temperature")
	Variable x0 = get_xye_number(buf,"x_min")
	Variable dx = get_xye_number(buf,"x_step")

	LoadWave/D/G/K=1/W/A/O/Q fullFileName
	if (ItemsInList(S_waveNames)<3)
		return ""
	endif
	String xAxisName=StringFromList(0,S_waveNames)
	Wave wx = $xAxisName
	Wave wy = $StringFromList(1,S_waveNames)
	Wave we = $StringFromList(2,S_waveNames)
	Variable xlo=WaveMin(wx), xhi=WaveMax(wx)
	KillWaves/Z wx

	if (printIt)
		printf "Load_xye_file(\"%s\")\r",fullFileName
		if (WaveExists(we))
			printf "	loaded waves:  %s, %s\r",NameOfWave(wy), NameofWave(we)
		else
			printf "	loaded wave:  %s\r",NameOfWave(wy)
		endif
		printf "	x-axis is \"%s\"  range=[%g, %g]\r",xAxisName, xlo,xhi
	endif

	wnote = ReplaceStringByKey("xAxisName",wnote,xAxisName,"=")
	wnote = ReplaceStringByKey("yAxisName",wnote,StringFromList(1,S_waveNames),"=")
	String allUnits="Q:1/"+ARING+";d:"+ARING+";"

	String xunits = SelectString(Stringmatch(xAxisName,"*theta*"),"",DEGREESIGN)
	xunits += SelectString(Stringmatch(xAxisName,"*chi*"),"",DEGREESIGN)
	xunits += SelectString(Stringmatch(xAxisName,"*phi*"),"",DEGREESIGN)
	if (strlen(xunits)<1)
		xunits = StringByKey(xAxisName,allUnits)
	endif

	SetScale/P x x0,dx,StringByKey(xAxisName,allUnits), wy,we
	Note/K wy, ReplaceStringByKey("errWave",wnote,NameOfWave(we),"=")
	Note/K we, ReplaceStringByKey("waveClass",wnote,"xyeErr","=")

	if (strlen(GetRTStackInfo(2))==0)
		Graph_xye(wy)		// Graph the an xye wave
	endif

	return StringFromList(1,S_waveNames) + ";" + StringFromList(2,S_waveNames)
End
//
Static Function/S addNumber2WaveNote(buf,wnote,key)
	String buf
	String wnote
	String key

	Variable val = get_xye_number(buf,key)
	if (numtype(val)<2)
		wnote = ReplaceNumberByKey(key,wnote,val,"=")
	endif
	return wnote
End
//
Static Function/S addString2WaveNote(buf,wnote,key)
	String buf
	String wnote
	String key

	String str = get_xye_string(buf,key)
	if (strlen(str))
		wnote = ReplaceStringByKey(key,wnote,str,"=")
	endif
	return wnote
End
//
Static Function get_xye_number(buf,key)
	String buf
	String key
	return str2num(get_xye_string(buf,key))
End
//
Static Function/S get_xye_string(buf,key)
	String buf
	String key

	Variable i0, i1
	String find = "# "+key+"="
	i0 = strsearch(buf,find,0)
	i1 = strsearch(buf,"\n",i0)
	if (i0<0 || i1<0)
		return ""
	endif
	return TrimBoth(buf[i0+strlen(find),i1])
End


Function/S Graph_xye(ww,[forceNew])		// Graph the an xye wave
	Wave ww
	Variable forceNew								// if true forces creation of new graph, not just bring existing graph to front
	forceNew = ParamIsDefault(forceNew) || numtype(forceNew) ? 0 : forceNew

	String name
	if (!WaveExists(ww))

		String list=WaveListClass("xye","*","DIMS:1")
		if (ItemsInList(list)==1)
			Wave ww=$StringFromList(0,list)
		else
			Prompt name,"Wave signal", popup, list
			DoPrompt "xye wave", name
			if (V_flag)
				return ""
			endif
			Wave ww=$name
		endif
	endif
	if (!WaveExists(ww))
		return ""
	endif

	String win=StringFromList(0,WindowsWithWave(ww,1))
	if (strlen(win) && !forceNew)			// if ww is already on a graph bring to front and quit
		DoWindow/F $win
		return win
	endif

	String wnote=note(ww)
	String xAxisName = StringByKey("xAxisName",wnote,"=")
	String yAxisName = StringByKey("yAxisName",wnote,"=")
	String xunits=WaveUnits(ww,0), yunits=WaveUnits(ww,1)

	Display /W=(35,45,430,253) ww
	ModifyGraph tick=2, mirror=1, minor=1, lowTrip=0.001
	Label bottom NameUnits2Label(xAxisName,xunits)
	Label left NameUnits2Label(yAxisName,yunits)

	Wave err = $StringByKey("errWave",wnote,"=")
	if (WaveExists(err))
		ErrorBars $NameOfWave(ww) Y,wave=(err,err)
	endif
	DoUpdate

	String str, title						// build the title for graph
	title = StringByKey("SPECfile",wnote,"=")
	title = RemoveTrailingString(title,".spec",1)
	str = StringByKey("Scan_numbers",wnote,"=")
	title += SelectString(strlen(title) && strlen(str), "", ": ") + str
	str = StringByKey("Project",wnote,"=")
	title += SelectString(strlen(title) && strlen(str), "", "\r\\Zr075") + str
	TextBox/C/N=title/F=0/A=LT/X=5.42/Y=5.52 title

	return StringFromList(0,WinList("*",";","WIN:1"))
End
//
Static Function/S NameUnits2Label(name,units)
	String name
	String units

	String out=name
	if (strlen(name) && strlen(units))
		out = name + SelectString(StringMatch(units,"¡"), "  (\\U)", "\\U")
	elseif (strlen(units))
		out = "(\\U)"
	endif
	return out
End
