#pragma rtGlobals=1		// Use modern global access method.
Strconstant EOF = "__EOF__"
#pragma version = 0.02


Menu "Data"
	"Load EPICS Scan",LoadEPICSscan("","",0)
	"Info about EPICS Scan",EPICSscanInfo($"")
End
Menu "Load Waves"
	"Load EPICS Scan",LoadEPICSscan("","",0)
End


// The following is an example of an EPICS_GraphStyle() function.  Use this as a template for making your own.
//
//	Function EPICS_GraphStyle() : GraphStyle
//		if (NumberByKey("IGORVERS",IgorInfo(0))>=5.01)
//			ModifyGraph gfMult=130
//		else
//			ModifyGraph gfSize=18
//		endif
//		ModifyGraph/Z tick=2, minor=1, standoff=0, lowTrip=0.001, mirror=1
//		String yName = StringFromList(0,TraceNameList("",";",1))
//		Wave wy = TraceNameToWaveRef("",yName)
//		Wave wx = XWaveRefFromTrace("", yName )
//		yName = StringByKey("description",note(wy),"=")				// name to use for axis label
//		yName = SelectString(strlen(yName),NameOfWave(wy),yName)	
//		String units = WaveUnits(wy, 1)
//		Label left yName+SelectString(strlen(units),"  \\U","  (\\U)")
//		String xName
//		if (WaveExists(wx))
//			xName = StringByKey("description",note(wx),"=")			// name to use for axis label
//			xName = SelectString(strlen(xName),NameOfWave(wx),xName)	
//			units = WaveUnits(wx, 1)
//		else
//			xName="points"
//			units = ""
//		endif
//		Label bottom xName+SelectString(strlen(units),"  \\U","  (\\U)")
//		String noteStr = note(wy)
//		TextBox/N=text0/F=0/S=3/B=1/A=MT StringByKey("file",noteStr,"=")+"\r"+StringByKey("date",noteStr,"=")
//		return 0
//	End


// Note, this routine is intended solely as a user interface for looking around. If you just want a particular value from a 
// particular wave, use the Igor built in function as         result = StringByKey(keyStr,note(wav),"=")
// or optionally if the result will be a number,           result = NumberByKey(keyStr,note(wav),"=")
Function/S EPICSscanInfo(wav)
	Wave wav
	if (!WaveExists(wav))
		String wList=ListOfWavesForPopup()+WaveList("*",";","DIMS:1")	// wave list including subdirectories
		String wName=StringFromList(0,wList)
		if (ItemsInList(wList)>1)
			Prompt wName, "scan wave to use for values",popup,wList
			DoPrompt "select scan",wName
			if (V_flag)
				return ""
			endif
		endif
		Wave wav=$wName
	endif

	String noteStr = note(wav)
	String pvList=""			// a list of all of the PVs (first part with no extensions)
	String pvListValues = ""		// keyword list pv1a=EXT1:avalue1,EXT1b:value1b;pv2=EXT2a:value2a,EXT2b:value2b;...

	// make list of all PV's,  a PV will be is any key with a colon in it, only save up to (but no including the '.')
	String pv,item,ext,value
	Variable i,j, N = ItemsInList(noteStr)
	for (i=0;i<N;i+=1)
		item = StringFromList(i,noteStr)
		pv = StringFromList(0,item,"=")
		value = StringFromList(1,item,"=")
		if (strsearch(pv,":",0)<0)										// skip items without a colon
			continue
		endif
		j = strsearch(pv,".",0)
		j = j>0 ? j-1 : Inf
		ext = SelectString(numtype(j),pv[j+2,Inf],"VAL")
		pv = pv[0,j]

		if (WhichListItem(pv,pvList)<0)
			pvList += pv+";"
		endif
		item = StringByKey(pv,pvListValues,"=")
		item = ReplaceStringByKey(ext,item,value,"�",",")				// add the ext:value pair to this item
		pvListValues = ReplaceStringByKey(pv,pvListValues,item,"=")	// update value of pv in pvListValues
	endfor

	// alphabetize the pvList
	Variable Npv = ItemsInList(pvList)
	String tempPvList = UniqueName("tempPvList",1,0)
	Make/T/N=(Npv) $tempPvList
	Wave/T tempList = $tempPvList
	tempList = StringFromList(p,pvList)
	Sort/A tempList, tempList
	pvList = ""
	for (i=0;i<Npv;i+=1)
		pvList += tempList[i]+";"
	endfor
	KillWaves/Z tempList

	// now include all of the keys in the list that were not considered a PV
	String moreKeys=""
	for (i=0;i<N;i+=1)
		item = StringFromList(i,noteStr)
		pv = StringFromList(0,item,"=")
		j = strsearch(pv,".",0)
		j = j>0 ? j-1 : Inf
		pv = pv[0,j]
		if (WhichListItem(pv,pvList)>=0)			// already in PVlist
			continue
		endif
		moreKeys += pv+";"
	endfor

	// select a pv
	pv = ""
	Prompt pv, "select a pv", popup,moreKeys+pvList
	DoPrompt "select a pv",pv
	if (V_flag)
		return ""
	endif
	item = ReplaceString("�", StringByKey(pv, pvListValues,"="), "=")
	item = ReplaceString(",", item, ";")
	if (!strlen(item))
		item = StringByKey(pv, noteStr,"=")
	endif
	item = RemoveEnding(item,";")
	if (ItemsInList(GetRTStackInfo(0))<2)
		printf "for the PV '%s', value%s   %s\r",pv,SelectString(strsearch(item,";",0)<0,"s are"," is"),item
	endif
	return item
End
Static Function/S ListOfWavesForPopup()
	String fullList = ""
	String fullName
	Variable m,i
	String fldrSav = GetDataFolder(1)
	String wList,wName,fldr,fldrList=StringByKey("FOLDERS", DataFolderDir(1))
	do
		fldr = StringFromList(m,fldrList,",")
		if (!strlen(fldr))
			break
		endif
		SetDataFolder $fldr
		wList = StringByKey("WAVES", DataFolderDir(2))
		SetDataFolder $fldrSav
		for (i=0,wName="x";strlen(wList)&&strlen(wName);i+=1)
			wName = StringFromList(i,wList,",")
			fullName = ":"+fldr+":"+wName
			Wave wav = $(":"+fldr+":"+wName)
			if (strsearch(note(wav),"scanRecord=",0)==0)		// found one, add it to the list
				fullList += GetWavesDataFolder(wav,4)+";"
				break
			endif
		endfor
		m += 1
	while(1)
	return fullList
End




Function/S LoadEPICSscan(path,fileName,reuseOld)
	String path					// igor file path to the file
	String fileName				// name of the file (full or partial relative to path)
	Variable reuseOld			// if true, then just use the old existing data, do not ask about re-reading

	if (strlen(path)<1)
		PathInfo scan
		if (V_flag)
			path = "scan"
		else
			path = "home"
		endif
	endif
	Variable refNum=0
	Open/Z=1/R/P=$path refNum as fileName								// try to open the file as passed
	if (V_flag)
		Open/Z=2/R/T=".dat"/M="scan file from TkGrab"/P=$path refNum	// select a fie and open it
		if (V_flag==-1)
			return ""
		endif
	endif
	if (strlen(S_fileName)<1)
		Abort "could not open the file '"+fileName+"'"
	endif
	S_fileName = ParseFilePath(5,S_fileName,":",0,0)						// convert to mac style separators
	fileName = ParseFilePath(3,S_fileName,":",0,0)						// get just the filename part, none of the path
	NewPath /O/Q/Z scan, ParseFilePath(1,S_fileName,":",1,0)
	String line,str

	FReadLine refNum, line
	line = RemoveEnding(line,"\r")
//	Variable i = strsearch(line,"# EPICS scan record from ",0)
	Variable i = strsearch(line,"# EPICS scan record:",0)
	if (i<0)
		DoAlert 0, "Invalid scan file"
		Close refNum
		return ""
	endif

	Variable N=strlen(line)
	for (i=strsearch(line,":",0)+1;char2num(line[i])<=32 && i<N;i+=1)
	endfor
	String scanRecord = line[i,Inf]
//	String scanRecord = line[25,Inf]
	String noteStr = "", scanPVs = "", item
	noteStr = ReplaceStringByKey("scanRecord",noteStr,scanRecord,"=")
	scanPVs = "user;date;file;scanPV;drive_positioner;number_of_points_requested;number_of_points_completed;"
	scanPVs += "positioner_settling_delay;detector_settling_delay;status_message"
	scanPVs += "drive_1_positioner;drive_1_units;drive_2_positioner;drive_2_units;"
	scanPVs += "drive_3_positioner;drive_3_units;drive_4_positioner;drive_4_units"
	i = 0
	do
		item = StringFromList(i,scanPVs)
		if (!strlen(item))
			break
		endif
		str = RetrieveFromFile(refNum,item,0)
		if (strlen(str) && !stringmatch(str,EOF))
			noteStr = ReplaceStringByKey(item,noteStr,str,"=")
		endif
		i += 1
	while(1)

	Variable completed,requested
	requested = NumberByKey("number_of_points_requested",noteStr,"=")
	completed = NumberByKey("number_of_points_completed",noteStr,"=")
	if (completed!=requested)
		sprintf str, "in file '', only %d of the requested %d points were measured",fileName,completed,requested
		DoAlert 0,str
		print str
	endif

//	noteStr = ReplaceStringByKey("user",noteStr,RetrieveFromFile(refNum,"user",0),"=")
//	noteStr = ReplaceStringByKey("date",noteStr,RetrieveFromFile(refNum,"date",0),"=")
//	noteStr = ReplaceStringByKey("file",noteStr,RetrieveFromFile(refNum,"file",0),"=")
//	noteStr = ReplaceStringByKey("scanPV",noteStr,RetrieveFromFile(refNum,"scanPV",0),"=")
//	noteStr = ReplaceStringByKey("drive_positioner",noteStr,RetrieveFromFile(refNum,"drive_positioner",0),"=")
//	noteStr = ReplaceStringByKey("number_of_points_requested",noteStr,RetrieveFromFile(refNum,"number_of_points_requested",0),"=")
//	noteStr = ReplaceStringByKey("number_of_points_completed",noteStr,RetrieveFromFile(refNum,"number_of_points_completed",0),"=")
//	noteStr = ReplaceStringByKey("positioner_settling_delay",noteStr,RetrieveFromFile(refNum,"positioner_settling_delay",0),"=")
//	noteStr = ReplaceStringByKey("detector_settling_delay",noteStr,RetrieveFromFile(refNum,"detector_settling_delay",0),"=")
//	noteStr = ReplaceStringByKey("user",noteStr,RetrieveFromFile(refNum,"user",0),"=")
//	noteStr = ReplaceStringByKey("user",noteStr,RetrieveFromFile(refNum,"user",0),"=")
//	noteStr = ReplaceStringByKey("user",noteStr,RetrieveFromFile(refNum,"user",0),"=")
//	noteStr = ReplaceStringByKey("user",noteStr,RetrieveFromFile(refNum,"user",0),"=")
//	str = RetrieveFromFile(refNum,"status_message",0)
//	if (strlen(str))
//		noteStr = ReplaceStringByKey("status_message",noteStr,str,"=")
//	endif
	String tagStr
	str = RetrieveFromFile(refNum,"##",0)
	for (;strlen(str)>0 && !stringmatch(str,EOF);)
		i = strsearch(str," ",0)-1
		if (i>=0)
			tagStr = str[0,i]
			str = TrimLeadingWhiteSpace(str[i+2,Inf])
			// change {str} to just str
			if (stringmatch(str,"{*}"))
				str = str[1,strlen(str)-2]
			endif
			noteStr = ReplaceStringByKey(tagStr,noteStr,str,"=")
		endif
		str = RetrieveFromFile(refNum,"##",-1)
	endfor

	// find column headings:
	String heading
	FSetPos refNum, 0
	do
		heading = line							// save the last line, which will turn out to be the headings
		FReadLine refNum, line
	while(char2num(line[0])==35)			// continue while line starts with a '#'
	heading = RemoveEnding(heading[2,Inf],"\r")

	Variable Ncol = ItemsInList(heading,"\t")	// number of columns
	Variable Nlen = NumberByKey("number_of_points_completed",noteStr,"=")	// length of waves to read
	if (Ncol<1)
		DoAlert 0, "the header is there, but no scan data found"
		Close refNum
		return ""
	endif
	String wNames = ""							// will hold a list of the new wave names
	String pv, pvList = ""						// original pv for each element in wNames
	for (i=0;i<Ncol;i+=1)
		str = StringFromList(i,heading,"\t")
		if (strlen(str)<1)
			DoAlert 0,"an empty column heading"
			Close refNum
			return ""
		endif
		pv = RemoveEnding(str,".VAL")
		pvList = AddListItem(pv,pvList,";",Inf)
		wNames = AddListItem( CleanupName(pv,0),wNames,";",Inf)
	endfor

	// time to make the new waves, so first create the sub-folder
	String fldr, fldrSav = GetDataFolder(1)
	fldr = CleanUpName(fileName,0)
	if (DataFolderExists(fldr))
		if (reuseOld)
			str = fldrSav+fldr+":;"
			SetDataFolder fldr
			i = strsearch(scanRecord,":",0)
			String crate = scanRecord
			if (i>0)
				crate = scanRecord[0,i-1]
			endif
			str += WaveList(crate+"_*",";","DIMS:1")
			SetDataFolder $fldrSav
			return str
		else
			DoAlert 1, "Data folder '"+fldr+"' already exists, re-read the data?"
			if (V_flag==2)
				Close refNum
				return ""
			endif
		endif
	endif
	NewDataFolder/O/S $fldr
	Variable j
	String egu,desc,wNote,pv0
	for (i=0;i<Ncol;i+=1)
		str = StringFromList(i,wNames)
		Make/N=(Nlen)/O $str
		Wave wav=$str
		Note/K wav
		pv = StringFromList(i,pvList)
		if (!strlen(pv))
			Note wav, noteStr
			continue
		endif

		wNote = ReplaceStringByKey("PV",noteStr,pv,"=")
		j = strsearch(pv,".",0)
		pv0 = SelectString(j>=1,pv,pv[0,j-1])	// part of pv preceeding possibl '.'
		if (stringmatch(pv,"*:vsc:c*.S*"))			// a counter channel, do something special
			j = str2num(pv[ strsearch(pv,".S",0)+2])
			desc = StringByKey(pv0+".NM"+num2istr(j),wNote,"=")
			egu = "counts"
		else
			desc = StringByKey(pv0+".DESC",wNote,"=")
			egu = StringByKey(pv0+".EGU",wNote,"=")
		endif
		if (strlen(desc))
			wNote = ReplaceStringByKey("description",wNote,desc,"=")
		endif
		if (strlen(egu))
			wNote = ReplaceStringByKey("units",wNote,egu,"=")
			SetScale d 0,0,egu, wav
		endif
		Note wav, wNote
	endfor

	Variable jline,ind
	for (jline=0;jline<Nlen;jline+=1)			// for each line
		ind = 0
		for (i=0;i<Ncol;i+=1)
			Wave wav=$StringFromList(i,wNames)
			wav[jline] = str2num(line[ind,Inf])
			ind= strsearch(line,"\t",ind)+1
		endfor
		FReadLine refNum, line
	endfor
	Close refNum
	SetDataFolder $fldrSav

	str = GetWavesDataFolder(wav,1)+";"+wNames
	if (ItemsInList(GetRTStackInfo(0))<2)
		printf "Loaded scan from scan record '%s',  %d waves %d long in data folder '%s'\r",scanRecord,Ncol,Nlen,GetWavesDataFolder(wav,1)
		// printf "to graph it use:			Display %s vs %s\r",GetWavesDataFolder(yw,2),GetWavesDataFolder(xw,2)
		PutUpGraph(str)
	endif
	return str
End

// returns value after tag.  "" means no value after tag,  but returns EOF if tag not found
Static Function/S RetrieveFromFile(refNum,tagStr,start)
	Variable refNum
	String tagStr
	Variable start		// 0=start at top,  -1=start from where you are

	if (start==0)
		FSetPos refNum, 0
	elseif (start!=-1)
		Abort "in RetrieveFromFile(), start must be only 0 or -1, not "+num2str(start)
	endif

	String line
//	tagStr = SelectString(stringmatch(tagStr,"##"),"# "+tagStr+": ","## ")
	tagStr = SelectString(stringmatch(tagStr,"##"),"# "+tagStr+":","## ")
	do
		FReadLine refNum, line
		line = RemoveEnding(line,"\r")
		if (strlen(line)<1 || char2num(line)!=35)		// EOF or does not start with '#', so all done searching
			return EOF
		endif
		if (strsearch(line,tagStr,0)==0)				// found it
			break
		endif
	while (1)
	line = line[strlen(tagStr),Inf]
	if (char2num(line)==0)
		line = ""
	endif
	Variable i, N=strlen(line)
	for (i=0;char2num(line[i])<=32 && i<N;i+=1)
	endfor
	return line[i,Inf]
End


Static Function/S getPVinfoAboutWave(wav,EXT)
	Wave wav					// wave to get wave information about
	String EXT					// the .XXX part, whole PV is PV.EXT, note EXT does not include the '.'

	String noteStr = note(wav)
	String pvName = StringByKey("PV", noteStr,"=")
	Variable i = strsearch(pvName,".",0)
	pvName = SelectString(i>=1,pvName,pvName[0,i-1])+"."+EXT
	return StringByKey(pvName, noteStr,"=")
End


Static Function/S TrimLeadingWhiteSpace(str)	// remove any leading white space from str
	String str
	Variable i
	i = -1
	do
		i += 1
	while (char2num(str[i])<=32)
	return str[i,strlen(str)-1]
End


Static Function PutUpGraph(list)
	String list

	String fldr = StringFromList(0,list)
	String wList = list[strlen(fldr)+1,Inf]
	String xName="", yName=""
	Variable N=ItemsInList(wList)
	Variable ix=0,iy=0
	if (N<1)
		DoAlert 0, "No waves to plot"
		return 1
	elseif (N==1)
		yName = StringFromList(0,wList)
	else
		xName = StringFromList(0,wList)
		yName = StringFromList(1,wList)
		ix = 0
		iy = 1
	endif
	if (stringmatch(yName,"*_T") && N>2)
		yName = StringFromList(2,wList)
		iy = 2
	endif

	String desc
	if (N>1)
		String item,popList=""
		Variable i
		for (i=0;i<N;i+=1)
			item = StringFromList(i,wList)
			Wave wav=$(fldr+item)
			desc = StringByKey("description",note(wav),"=")
			item += SelectString(strlen(desc),"","  "+desc)
			popList += item+";"
		endfor
		xName = StringFromList(ix,popList)
		yName = StringFromList(iy,popList)
		Prompt xName, "name of X wave", popup,"_none_;"+popList
		Prompt yName, "name of Y wave", popup,popList
		DoPrompt "pick waves to plot",xName,yName
		if (V_flag)
			return 1
		endif
		xName = SelectString(stringmatch(xName,"_none_"),xName,"")
		xName = StringFromList(0,xName," ")
	endif

	Variable useX = strlen(xName)						// flags that an x-wave is present
	Wave wy = $(fldr+yName)
	Wave wx = $SelectString(useX,UniqueName("xxx",1,0),(fldr+xName))
	if (!WaveExists(wy) || (useX && !WaveExists(wx)))
		DoAlert 0, "'"+yName+"' or '"+xName+"' does not exist"
		return 1
	endif

	if (useX)
		Display wy vs wx
	else
		Display wy
	endif
	if (exists("EPICS_GraphStyle")==5 || exists("EPICS_GraphStyle")==6)
		Execute "EPICS_GraphStyle()"
		return 0
	endif

	if (NumberByKey("IGORVERS",IgorInfo(0))>=5.01)
		ModifyGraph gfMult=130
	else
		ModifyGraph gfSize=18
	endif
	ModifyGraph/Z tick=2, minor=1, standoff=0, lowTrip=0.001, mirror=1

	desc = StringByKey("description",note(wy),"=")				// name to use for axis label
	yName = SelectString(strlen(desc),NameOfWave(wy),desc)	
	String units = WaveUnits(wy, 1)
	Label left yName+SelectString(strlen(units),"  \\U","  (\\U)")
	if (useX)
		desc = StringByKey("description",note(wx),"=")			// name to use for axis label
		xName = SelectString(strlen(desc),NameOfWave(wx),desc)	
		units = WaveUnits(wx, 1)
	else
		xName="points"
		units = ""
	endif
	Label bottom xName+SelectString(strlen(units),"  \\U","  (\\U)")
	String noteStr = note(wy)
	TextBox/N=text0/F=0/S=3/B=1/A=MT StringByKey("file",noteStr,"=")+"\r"+StringByKey("date",noteStr,"=")
	return 0
End

//	# EPICS scan record from hnl:scan1
//	# user: usaxs (Ultra-Small-Angle X-ray Scattering user,33ID-D)
//	# date: Tue Oct 12 10:31:56 2004
//	# file: ./hnl_scan1_20041012_103152.dat
//	# scanPV: hnl:scan1
//	# drive_positioner: hnl:md:energy.VAL
//	# readback_positioner: hnl:md:energy.RBV
//	# status_message: Scanning ...
//	# number_of_points_requested: 402
//	# number_of_points_completed: 134
//	# positioner_settling_delay:  0.5
//	# detector_settling_delay:    1
//	## ID34:Energy.EGU  keV
//	## ID34:Energy.VAL  9.089093208312988e+00
//	## ID34:Gap.EGU  mm
//	## ID34:Gap.VAL  2.032250000000000e+01
//	## ID34:TaperGap.EGU  mm
//	## ID34:TaperGap.VAL  2.784231475196464e-03
//		.
//		.
//		.
//	## yum:wire:center.DESC  spot center
//	## yum:wire:center.EGU  micron
//	## yum:wire:center.VAL  -4.051681526136727e+03
//	## yum:wire:dia.DESC  wire dia.
//	## yum:wire:dia.EGU  micron
//	## yum:wire:dia.VAL  1.089695452651313e+01
//	# hnl:md:energy.VAL	hnl:md:energy.RBV	koa:vsc:c0_calc2.VAL	koa:vsc:c0_calc3.VAL	koa:vsc:c0_calc4.VAL	koa:vsc:c0_calc5.VAL	koa:vsc:c0_calc6.VAL
//	8.183	46284.9	150627	192440	-355943	25.9559
//	8.188	46114.9	151177	193643	-420103	25.9579
//	8.193	45928.9	152022	195388	-377396	25.9579
//	8.198	46410.9	155143	200076	-387315	25.9559
//	8.203	46253.9	156043	201959	-504820	25.9559
//	8.208	46195.9	157329	204393	-278077	25.9559
//	8.213	46085.9	158442	206537	-455424	25.9559
//	8.218	46138.9	160136	209494	-157456	25.9559
//	8.223	45965.9	161002	211405	-263875	25.9579
//		.
//		.
//		.
