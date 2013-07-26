#pragma rtGlobals=1		// Use modern global access method.
#pragma ModuleName=ExamineImageSet
#pragma version = 1.00


Menu "micro", dynamic
	"Scanned Ranges from Images",GetScannedRangesFromImages("","","","")
	help={"find the rnage and order ot what was scanned"}
	"  Find moving axes from image files in a folder",WhatsMovingInFolder("")
	help={"identify which axes were scanned by looking at all of the image files"}
End


Function/T ReadGenericImageHeader(name)
	String name
#if (Exists("WinViewReadHeader")==6)
	return WinViewReadHeader(name)
#elif (Exists("HDFReadHeader")==6)
	return HDFReadHeader(fileName)
#else
	return ""
#endif
End


Function/S WhatsMovingInFolder(pathName)							// returns a list of the moving axes
	String pathName

	if (strlen(pathName)<1 && ItemsInList(GetRTStackInfo(0))<2)				// ask user to find path
		Prompt pathName, "path to use",popup,PathList("!materials", ";", "" )+"-- make new path --"
		DoPrompt "path",pathName
		if (V_flag)
			return ""
		endif
		if (stringmatch(pathName,"-- make new path --"))
			pathName = "imagePath"
			NewPath/M="path to image files"/Z imagePath
			if (V_flag)
				return ""
			endif
		endif
	endif
//	pathName = SelectString(strlen(pathName),"imagePath",pathName)	// for empty pathName, use "imagePath"
//	PathInfo $pathName
//	if (!V_flag && stringmatch(pathName,"imagePath"))
//		NewPath/M="path to image files"/Z $pathName
//	endif
	PathInfo $pathName
	if (!V_flag)
		DoAlert 0, "path  '"+pathName+"'  does not exists"
		return ""
	endif
	String list = directory(pathName)							// get list of all files in this path
	Variable m, Nm=ItemsInLIst(list)

	String name,tags = "X1;Y1;Z1;H1;F1;X2;Y2;Z2;H2;F2;keV;depthSi;CCDy"
	Variable i,Ni=ItemsInLIst(tags)
	for (i=0;i<Ni;i+=1)
		name = StringFromList(i,tags)+"moving"
		Make/N=(Nm)/O/D $name
		Wave w=$name
		w = NaN
	endfor

	// all possible axes to check
	String header
	for (m=0;m<Nm;m+=1)
		header = ReadGenericImageHeader(S_path+StringFromList(m,list))
		for (i=0;i<Ni;i+=1)
			name = StringFromList(i,tags)
			Wave w=$(name+"moving")
			w[m] = NumberByKey(name,header,"=")
		endfor
	endfor

	// remove unused axes
	for (i=Ni-1;i>=0;i-=1)
		name = StringFromList(i,tags)
		Wave w=$(name+"moving")
		WaveStats/M=1/Q w
		if ((V_max-V_min)<1e-5 || numtype(V_min))
			tags = RemoveFromList(name,tags)		// axis not used, remove it
			KillWaves/Z w
		elseif (strlen(name)==2 && !isdigit(name[0]) && isdigit(name[1]) && (V_max-V_min)<0.15)
			tags = RemoveFromList(name,tags)		// axis not used, remove it
			KillWaves/Z w
		endif
	endfor

	String removeList=""							// remove rotated coordinates where appropriate
	if (WhichListItem("H1",tags)>=0 &&WhichListItem("F1",tags)<0)		// pure H1, delete Y1, Z1, and F1
		removeList += "Y1;Z1;F1;"
	elseif (WhichListItem("H1",tags)<0 && WhichListItem("F1",tags)>=0)	// pure F1, delete Y1, Z1, and H1
		removeList += "Y1;Z1;H1;"
	elseif (WhichListItem("Y1",tags)>=0 && WhichListItem("Z1",tags)<0)	// pure Y1, delete F1, H1, and Z1
		removeList += "H1;F1;Z1;"
	elseif (WhichListItem("Z1",tags)<0 && WhichListItem("Z1",tags)>=0)	// pure Z1, delete H1, F1, and Y1
		removeList += "H1;F1;Y1;"
	endif
	if (WhichListItem("H2",tags)>=0 &&WhichListItem("F2",tags)<0)		// pure H2, delete Y2, Z2, and F2
		removeList += "Y2;Z2;F2;"
	elseif (WhichListItem("H2",tags)<0 && WhichListItem("F2",tags)>=0)	// pure F2, delete Y2, Z2, and H2
		removeList += "Y2;Z2;H2;"
	elseif (WhichListItem("Y2",tags)>=0 && WhichListItem("Z2",tags)<0)	// pure Y2, delete F2, H2, and Z2
		removeList += "H2;F2;Z2;"
	elseif (WhichListItem("Z2",tags)<0 && WhichListItem("Z2",tags)>=0)	// pure Z2, delete H2, F2, and Y2
		removeList += "H2;F2;Y2;"
	endif
	for (i=0;i<ItemsInList(removeList);i+=1)
		name = StringFromList(i,removeList)
		tags = RemoveFromList(name,tags)			// axis not used, remove it
		KillWaves/Z $(name+"moving")
	endfor
	Ni = ItemsInLIst(tags)

	for (i=0;i<Ni;i+=1)
		name = StringFromList(i,tags)
		Wave w=$(name+"moving")
		WaveStats/M=1/Q w
		if (ItemsInList(GetRTStackInfo(0))<2)
			printf "for %s,  range = [%g,  %g],  Æ = %g\r",name,V_min,V_max,V_max-V_min
		endif
	endfor
	return tags
End



//Function test()
//	print GetScannedRangesFromImages("imagePath","BTO-F-Eshort_","1-1078","")
//End

Function/T GetScannedRangesFromImages(pathName,namePart,range1,range2)
	String pathName	// either name of path to images, or the full expliiciit path, i.e. "Macintosh HD:Users:tischler:data:cal:recon:"
	String namePart	// the first part of file name, something like  "EW5_"
	String range1		// range of first number after file root (designates energy)
	String range2		// range of second numbers after range1 (designates depth)

	Variable printIt=0
	String str

	pathName = SelectString(strlen(pathName),"imagePath",pathName)
	PathInfo $pathName
	if (!V_flag || strlen(namePart)<1)					// path does not exist or no namePart, ask user
		String pathPart
		str = requestFileRoot(pathName,2)
		pathPart = StringFromList(0,str)
		namePart = StringFromList(1,str)
		if (strlen(pathPart)<1 || strlen(namePart)<1)
			return ""										// invalid inputs
		endif
		if (!stringmatch(pathPart,S_path))				// path was changed
			if (stringmatch(pathName,"imagePath"))
				NewPath/O/M="path to reconstructed image files" imagePath pathPart	// for iamgePath, automatically reassign
			else
				NewPath/M="path to reconstructed image files" $pathName pathPart	// for other names, ask
			endif
		endif
		printIt = 1
	endif
	PathInfo $pathName
	if (strlen(S_path)<1 || strlen(namePart)<1)
		return ""											// invalid inputs
	endif
	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
		printf "using data from files starting with '%s'\r",S_path+namePart
	endif

	if (strlen(range1)<=0 || strlen(range2)<=0)		// if either range is empty, get the full range
		str = get_ranges(pathName,namePart)
		if (strlen(range1)<=0)
			range1=StringFromList(0,str)
			if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
				printf "using outer range of '%s'\r",range1[0,370]	// cannot print lines longer than 400 chars
			endif
			printIt = 1
		endif
		if (strlen(range2)<=0)
			range2=StringFromList(1,str)
			if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
				printf "using depth range of '%s'\r",range2[0,370]	// cannot print lines longer than 400 chars
			endif
			printIt = strlen(range2) ? 1 : printIt
		endif
	endif

	String fileRoot										// complete path up to and including the underscore
	PathInfo $pathName
	if (V_flag)											// if path exists, reset pathName to the explicit path, otherwise assume it is the explicit path
		fileRoot = S_path
	elseif (strsearch(pathName,":",Inf,1) == strlen(pathName)-1)
		fileRoot = pathName+":"							// ensure terminating colon
	endif
	fileRoot += namePart

	String scanList = FindAllScannedRanges(fileRoot,range1,range2)
	if (printIt)
		printf "found %d scanned variables,   scanList = '%s'\r",ItemsInList(scanList), scanList
	endif
	return scanList
End
//
//
//
Static Function/T FindAllScannedRanges(fileRoot,range1,range2)	// for the given set of images, find all of the scanned variables and their ranges, returns them ordered outer->inner
	String fileRoot	// COMPLETE root to file, i.e. the path+name_
	String range1		// range of first number after file root (designates energy)
	String range2		// range of second numbers after range1 (designates depth), this range can be empty for just a single range (non-reconstructed data)

	String name, wNote, varList="X1;Y1;Z1;H1;F1;X2;Y2;Z2;H2;F2;keV;depthSi;"
	varList="X1;H1;F1;keV;depthSi;"
	Variable Nk=ItemsInList(varList), NN=1000
	Make/N=(NN,Nk)/O scannedVariables
	Variable i, j, k, N, once
	for (j=str2num(range1),N=0; numtype(j)==0; j=NextInRange(range1,j))	// loop over range1, energies
		for (i=str2num(range2),once=1; numtype(i)==0 || once; i=NextInRange(range2,i),once=0)	// loop over range2, depths, always do at least once
			name = fileNameFromRangeNumbers(fileRoot,j,i)
			wnote=ReadGenericImageHeader(name)		// wave note to add to file read in
			if (strlen(wnote)<1)
				continue
			endif
			if (N>=NN)
				NN += 500
				Redimension/N=(NN,Nk) scannedVariables
			endif
			for (k=0;k<Nk;k+=1)
				scannedVariables[N][k] = NumberByKey(StringFromList(k,varList), wnote,"=")
			endfor
			N += 1
		endfor
	endfor
	Redimension/N=(N,Nk) scannedVariables

	Make/N=(N)/O tempVariable
	String var, str, scanList=""
	for (i=0;i<Nk;i+=1)									// check each variable
		tempVariable = scannedVariables[p][i]
		var = StringFromList(i,varList)
		if (constantWave(tempVariable))				// remove un-scanned variables
			varList = RemoveFromList(var,varList)
			scannedVariables[][i,Nk-2] = scannedVariables[p][q+1]
			Nk -= 1
			i -= 1
		else
			str = ReplaceString(";",periodOfScannedVariable(tempVariable),",")	// add to result
			scanList += var+"="+str+";"
		endif
	endfor
	Redimension/N=(N,Nk) scannedVariables

	// sort from outermost to innermost
	Make/N=(Nk)/T/O tempStringSort
	Make/N=(Nk)/O tempVariableSort
	tempStringSort = StringFromList(p,scanList)
	tempVariableSort = str2num(StringFromList(3,tempStringSort[p],","))
	Sort/R tempVariableSort, tempStringSort

	// now rebuild scanList based on the sort
	scanList = ""
	for (i=0;i<Nk;i+=1)									// check each variable
		scanList += StringFromList(0,tempStringSort[i],",")+","
		scanList += StringFromList(1,tempStringSort[i],",")+","
		scanList += StringFromList(2,tempStringSort[i],",")+";"
	endfor
	KillWaves/Z tempVariable, tempStringSort,tempVariableSort, scannedVariables
	return scanList
End
//
Static Function constantWave(w)
	Wave w
	WaveStats/Q w
	return (V_max==V_min || V_npnts==0)
End
//
Static Function/T periodOfScannedVariable(w)
	Wave w

	Variable start = w[0]==0 ? 0 : w[0]
	Variable N = numpnts(w)-1
	if (N<1)
		return ""
	endif
	Make/N=(N)/O periodOfScannedVariable_dw_
	Wave dw=periodOfScannedVariable_dw_

	Variable direction										// +1=forward, -1=backward
	dw = w[p+1]-w[p]
	WaveStats/Q/M=1 dw
	Variable maxVal = abs(V_max) > abs(V_min) ? V_max : V_min
	if (maxVal==0)
		return ""
	endif
	direction = maxVal>0 ? -1 : 1
	direction *= (V_max*V_min >=0) ? -1 : 1			// all in one direction, no re-trace

	// get the step size
	dw *= direction										// make the step I care about is the positive one
	dw = dw<=0 ? NaN : dw
	WaveStats/M=1/Q dw
	V_avg /= 10
	dw = dw<=V_avg ? NaN : dw							// this is to remove tiny jitter
	WaveStats/M=1/Q dw
	Variable step = V_avg*direction

	// count the number of steps	in one forward scan
	dw = (w[p+1]-w[p])
	Variable is=0, i, xstart=w[0], maxRange=0
	Variable NN=0										// NN is number of images in one scan of this variable, it will be used to find outer most & innermost loops
	for (i=0;i<N;i+=1)
		if ((direction>0 && dw[i]<-step) || (direction<0 && dw[i]>-step)) // on re-trace
			NN = max(NN,is)
			is = 0
			maxRange = max(maxRange,abs(w[i]-xstart))
			//printf "retrace %d:   [%g,   %g],  Æ=%g\r",i, xstart,w[i],w[i]-xstart
			xstart = w[i+1]
		else												// took one regular step
			is += 1
		endif
	endfor
	maxRange = (maxRange==0) ? abs(w[N-1]-xstart) : maxRange
	NN = max(NN,is)
	Variable Npnts = round(abs(maxRange/step))+1
	if (step==0 || Npnts<1)
		return ""
	endif

	//	printf "maxRange = %g,   start = %g,   direction = %g,    step = %g, Npnts=%g\r",maxRange,start,direction,step,Npnts
	KillWaves/Z periodOfScannedVariable_dw_
	String str
	sprintf str,"%d;%g;%g;%g", Npnts,start,step,NN
	return str
End



Static Function/T requestFileRoot(pathName,suffixes)
	String pathName
	Variable suffixes					// maximum number of suffixes to remove, each suffix looks like "_12", leave the underscore
										// suffixes are ALWAYS preceeded with an "_"
	PathInfo $pathName
	pathName = SelectString(V_flag,"",pathName)
	String message="pick any reconstructed image file in range,  using datafolder = "+GetDataFolder(0)
	Variable refNum
	Open/T=IMAGE_FILE_EXT/D/M=message/P=$pathName/R refNum
	if (strlen(S_fileName)<1)
		return ""
	endif
	String fullPath = ParseFilePath(1,S_fileName,":",1,0)
	String name = ParseFilePath(3,S_fileName,":",0,0)
	String fileRoot=name
	if (suffixes<1)						// nothing to do
		return fullPath+";"+fileRoot
	endif

	Variable i=1
	for (;suffixes>=1 && i>0;suffixes-=1)		// strip off all but last suffix
		i = strsearch(fileRoot,"_",Inf,-1)
		if (i>0)
			fileRoot = fileRoot[0,i-1]
		endif
	endfor
	return fullPath+";"+fileRoot+"_"
End


// get the range of file index numbers by looking at all of the files in a directory
// this only works for the last index on a file, not for depth sorted images
Static Function/T get_FilesIndexRange(pathName,namePart)
	String pathName					// probably 'imagePath'
	String namePart					// name part of file = "EW5_"

	String list = directory(pathName)

	Variable i,m,N = ItemsInList(list)
	Make/N=(N)/O/T fileList_get_FilesIndexRange
	Wave/T fileList = fileList_get_FilesIndexRange
	String name
	for (i=0,m=0;i<N;i+=1)
		name = StringFromList(i,list)
		if (strsearch(name,namePart,0)==0)
			fileList[m] = name			// file starts with namePart, save it in list
			m += 1						// increment pointer into fileList[]
		endif
	endfor
	N = m
	Redimension/N=(N) fileList			// set to correct length

	Make/N=(N)/U/I/O nw_get_FilesIndexRange
	Wave nw=nw_get_FilesIndexRange
	Variable k = strlen(namePart)
	for (m=0;m<N;m+=1)
		name = fileList[m]
		name = name[k,Inf]
		nw[m] = str2num(name)
		i = strsearch(name,"_",0)
	endfor
	Sort nw,nw

	list = ""								// make a list of all of the first number (with no repeats)
	k = NaN
	for (m=0;m<N;m+=1)
		i = nw[m]
		if (i!=k && i<4.2e+9)
			list += num2istr(i)+";"
			k = i
		endif
	endfor
	String range = compressRange(list,";")

	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
		print "range	",ItemsInRange(range),"		",range
	endif
	KillWaves/Z fileList_get_FilesIndexRange, nw_get_FilesIndexRange
	return range
End


Static Function/T get_ranges(pathName,namePart)
	String pathName					// probably 'imagePath'
	String namePart					// name part of file = "EW5_"

	String list = directory(pathName)

	Variable i,m,N = ItemsInList(list)
	Make/N=(N)/O/T fileList_get_ranges
	Wave/T fileList = fileList_get_ranges
	String name
	for (i=0,m=0;i<N;i+=1)
		name = StringFromList(i,list)
		if (strsearch(name,namePart,0)==0)
			fileList[m] = name			// file starts with namePart, save it in list
			m += 1						// increment pointer into fileList[]
		endif
	endfor
	N = m
	Redimension/N=(N) fileList			// set to correct length

	Make/N=(N)/U/I/O n1_get_ranges,n2_get_ranges
	Wave n1=n1_get_ranges, n2=n2_get_ranges
	Variable k = strlen(namePart)
	for (m=0;m<N;m+=1)
		name = fileList[m]
		name = name[k,Inf]
		n1[m] = str2num(name)
		i = strsearch(name,"_",0)
		n2[m] = i>0 ? str2num(name[i+1,Inf]) : 2.1e9
	endfor
	Sort n1,n1
	Sort n2,n2

	list = ""							// make a list of all of the first number (with no repeats)
	k = NaN
	for (m=0;m<N;m+=1)
		i = n1[m]
		if (i!=k && i<4.2e+9)
			list += num2istr(i)+";"
			k = i
		endif
	endfor
	String range1 = compressRange(list,";")

	list = ""							// make a list of all of the second number (with no repeats)
	k = NaN
	for (m=0;m<N;m+=1)
		i = n2[m]
		if (i!=k && i<2e+9)
			list += num2istr(i)+";"
			k = i
		endif
	endfor
	String range2 = compressRange(list,";")
	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
		print "range1	",ItemsInRange(range1),"		",range1
		print "range2	",ItemsInRange(range2),"		",range2
	endif
	KillWaves/Z fileList_get_ranges, n1_get_ranges, n2_get_ranges
	return range1+";"+range2
End


Static Function/T fileNameFromRangeNumbers(fileRoot,jkeV,iDepth)	// return full file name using jkeV and iDepth
	String fileRoot
	Variable jkeV,iDepth			// first and second indicies

	String name
	if (!(jkeV>=0 && jkeV<2e9))
		sprintf name,"%s%d%s",fileRoot,iDepth,IMAGE_FILE_EXT
	elseif (!(iDepth>=0 && iDepth<2e9))
		sprintf name,"%s%d%s",fileRoot,jkeV,IMAGE_FILE_EXT
	else
		sprintf name,"%s%d_%d%s",fileRoot,jkeV,iDepth,IMAGE_FILE_EXT
	endif
	return name
End


//	returns the directory of the path in a semi-colon separated list (works for both Mac & Windows)
Static Function/T directory(pathName)
	String pathName							// path name that I want the directory of
	PathInfo $pathName	

	String list=""
	if (stringmatch(IgorInfo(2),"Macintosh"))	// this section for Mac
		PathInfo $pathName
		String script = "get POSIX path of file \"" + S_path + "\""
		ExecuteScriptText/Z script
		if (V_Flag)
			return ""
		endif
		String POSIXpath = S_value[1,StrLen(S_value)-2] // trim quotes
		POSIXpath = ReplaceString(" ",POSIXpath,"\\\\ ")	// change spaces to escaped spaces
		String cmd
		sprintf cmd, "do shell script \"ls %s\"",POSIXpath
//		sprintf cmd, "do shell script \"ls \\\"%s\\\"\"",POSIXpath
		ExecuteScriptText/Z cmd
		if (V_Flag==0)
			// print POSIXpath, "        ",cmd
			list = S_value
			Variable i=strlen(list)
			list = list[1,i-2]						// remove leading and trailing double quote
			list = ReplaceString("\r", list, ";" )		// change to a semicolon separated list
		else
			print " the 'ls' command failed, you may want to add anther mount point to the directory() function to speed things up"
			list = IndexedFile($pathName,-1,IMAGE_FILE_EXT)
		endif
	else															// this sectio for Windows
		list = IndexedFile($pathName,-1,IMAGE_FILE_EXT)	// list = IndexedFile($pathName,-1,"????")

	endif
	return list
End
//
//Static Function/T directory(pathName)
//	String pathName							// path name that I want the directory of
//	PathInfo $pathName	
//
//	String pathStr,list=""
//	if (stringmatch(IgorInfo(2),"Macintosh"))	// this section for Mac
//		PathInfo $pathName
//		String cmd
//		pathStr = ReplaceString(":", S_path, "/" )
//		pathStr = pathStr[strsearch(pathStr,"/",0),Inf]
//
//		if (strsearch(S_path, "RAID:",0)==0)
//			pathStr = "/Volumes/red.uni.aps.anl.gov"+pathStr
//		elseif (strsearch(S_path,"J:",0)==0)
//			pathStr = "/Volumes/xmd"+pathStr
//		elseif (strsearch(S_path,"AUG06:",0)==0)
//			pathStr = "/Volumes/xmd"+pathStr
//		elseif (strsearch(S_path,"micro2TB:",0)==0)
//			pathStr = "/Volumes/micro2TB"+pathStr
//		endif
//
//		sprintf cmd, "do shell script \"ls \\\"%s\\\"\"",pathStr
//		ExecuteScriptText/Z cmd
//		if (V_Flag==0)
//			// print pathStr, "        ",cmd
//			list = S_value
//			Variable i=strlen(list)
//			list = list[1,i-2]						// remove leading and trailing double quote
//			list = ReplaceString("\r", list, ";" )		// change to a semicolon separated list
//		else
//			print " the 'ls' command failed, you may want to add anther mount point to the directory() function to speed things up"
//			list = IndexedFile($pathName,-1,".SPE")
//		endif
//	else											// this sectio for Windows
////		pathStr = ParseFilePath(5, S_path, "\\", 0, 0)
////		if (strlen(pathStr)<1)
////			return ""
////		endif
////		PathInfo Igor
////		String line,fname=S_path+"User Procedures:JZTbatchFile.bat"
////		Variable fid
////		Open/Z=1 fid as fname
////		if (V_flag)
////			return ""
////		endif
////		fprintf fid,"dir/b \"%s\" > \"C:\\Program Files\\WaveMetrics\\Igor Pro Folder\\User Procedures\\JZTbatchOut.txt\"\r\n",pathStr
////		Close fid
////		ExecuteScriptText/W=15 "\"C:\\Program Files\\WaveMetrics\\Igor Pro Folder\\User Procedures\\JZTbatchFile.bat\""
////		PathInfo Igor
////		fname=S_path+"User Procedures:JZTbatchOut.txt"
////		Open/R/Z=1 fid as fname
////		if (V_flag)
////			return ""
////		endif
////		list = ""
////		do
////			FReadLine fid, line
////			line = ReplaceString("\r",line,"")
////			list += ReplaceString("\r",line,"")+SelectString(strlen(line),"",";")
////		while (strlen(line))
////		Close fid
////		DeleteFile/Z "C:\\Program Files\\WaveMetrics\\Igor Pro Folder\\User Procedures\\JZTbatchFile.bat"
////		DeleteFile/Z "C:\\Program Files\\WaveMetrics\\Igor Pro Folder\\User Procedures\\JZTbatchOut.txt"
//
//		list = IndexedFile($pathName,-1,".SPE")		// list = IndexedFile($pathName,-1,"????")
//
//	endif
//	return list
//End