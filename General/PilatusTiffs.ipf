#pragma rtGlobals=1		// Use modern global access method.
#pragma version = 2.10
#pragma IgorVersion = 6.22					// actually needs 6.12 or greater, but using "= 6.12:" does not work
#pragma ModuleName=pilatus
// #initFunctionName "Init_PilatusImagesPackage()"
//  with version 2.00, It starts to use Christian's newer method, the old "*.pilatuslog" files are gone.
//  with version 2.02, uses only extras string to pass quiet & multiple
//	with version 2.04, changed how I deal with sitution where LoadPilatusTiffImage() but image already exists

strConstant PilatusfileFilters = "TIFF Files (*.tif, *.tiff):.tif,.tiff;All Files:.*;"
// #define UseOLDpilatuslogFiles		// Put this line in your main "Procedure Window" to enable the old way
// #definef EnablePilatusUnWrapping	// to enable the unwrappng functions



#ifndef UseOLDpilatuslogFiles
Menu "Load Waves"
	"Load Pilatuf spec Image...",LoadPilatusTiffImageSpec(NaN,NaN)
	"  Load Pilatus Tiff Image...",LoadPilatusTiffImage("")
End
#else
Menu "Load Waves"
	"Load Pilatus Tiff Image...",LoadPilatusTiffImage("")
End
#endif

//	//	example
//Function PilatusFixBadPixels(image)
//	Wave image
//	Variable bad = (image[255][27] + image[257][27] + image[256][26] + image[256][28])/4
//	image[256][27] = bad
//End

#ifndef UseOLDpilatuslogFiles
//  ============================================================================  //
//  ===================== Start of NEW style with pilatus info in spec file =====================  //

Function/T LoadPilatusTiffImage(fullName[,extras])
	String fullName
	String extras												// example		"quiet:1;multiple:0;"
	extras = SelectString(ParamIsDefault(extras),extras,"")
	Variable quiet = NumberByKey("quiet",extras)				// for quiet=1, if any problems just return
	Variable multiple = NumberByKey("multiple",extras)		// load multiple images, NOT YET WORKING
	quiet = numtype(quiet) ? 0 : !(!quiet)
	multiple = numtype(multiple) ? 0 : !(!multiple)

	GetFileFolderInfo/Q/Z=1 fullName							// check if file exists
	Variable printIt=0
	if (!V_isFile)
		if (quiet)
			return ""
		else
			Variable f
			Open /D/M="Pilatus tiff"/R/F=PilatusfileFilters f	// use /D to get full path name
			fullName = S_filename
			printIt = 1
			quiet = 0
		endif
	endif
	if (strlen(fullName)<1)									// no file name, quit
		return ""
	endif

	String fName = ParseFilePath(3,fullName,":",0,0), name	// only the name part, no path, no extension
	Variable i=strsearch(fName,"_S",Inf,3)
	if (i<0)
		return ""
	endif

	String wName
	name = CleanupName(fName[0,i-1],0)
	Variable scanNum,pindex
	sscanf fname[i+2,Inf], "%d_%d",scanNum,pindex
	if (V_flag==2)
		sprintf wName,"%s_S%d_%d",name,scanNum,pindex
		if (exists(wName)==1 && multiple)
			wName += "_"+num2istr(pindex)
		endif
		wName = CleanupName(wName,0)
	else
		scanNum = NaN
		pindex = NaN
		wName = CleanupName(fName,0)
	endif
	if (exists(wName)==1 || strlen(wName)<1)			// Do NOT allow overwriting
		if (!quiet)
			printf "'%s'    already exists,   nothing done\r",wName
		endif
		return "preExists:"+wName							// passes name of existing file, let higher routine decide what to do
		return ""
	endif

	String wnote=PilatusTiffReadHeader(fullName)			// get info from inside the TIFF file, does not load an image
	if (scanNum>0)
		wnote = ReplaceNumberByKey("scanNum",wnote,scanNum,"=")
	endif
	if (pindex>=0)
		wnote = ReplaceNumberByKey("scanPoint",wnote,pindex,"=")
	endif

	ImageLoad/Z/O/Q/T=tiff fullName						// this overwrites existing image
	if (!V_flag || V_numImages!=1)
		return ""
	endif
	Wave image = $StringFromList(0,S_waveNames)
	Rename image $wName
	if (printIt)
		printf "LoadPilatusTiffImage(\"%s\"",fullName
		if (multiple)
			printf "multiple=%g",multiple
		endif
		printf ") // wave = '%s'\r",wName
	endif

	String funcName="PilatusFixBadPixels"+SelectString(exists("PilatusFixBadPixels")==6,"Proto","")
	FUNCREF PilatusFixBadPixelsProto fixPixels=$funcName
	fixPixels(image)

	SetScale/P x 0,1,"pixel", image
	SetScale/P y 0,1,"pixel", image
	SetScale d 0,0,"<counts>", image
	Note/K image, wnote
	return GetWavesDataFolder(image,2)
End
//
Function PilatusFixBadPixelsProto(image)
	Wave image
End


Function/T PilatusTiffReadROI(fullName,i0,i1,j0,j1,[extras])
	String fullName
	Variable i0,i1,j0,j1			// pixel range of ROI (if i1 or j1<0 then use whole image)
	String extras					// example		"quiet:1;multiple:0;"
	extras = SelectString(ParamIsDefault(extras),extras,"")

	if (ParamIsDefault(extras))
		Wave image = $LoadPilatusTiffImage(fullName)
	else
		Wave image = $LoadPilatusTiffImage(fullName,extras=extras)
	endif
	String wnote = note(image)

	Variable ydim, xdim									// x,y size of whole array
	xdim = DimSize(image,0)
	ydim = DimSize(image,1)
	i1 = (i1<1) ? xdim-1 : i1							// -1 flags use whole range
	j1 = (j1<1) ? ydim-1 : j1

	i0 = max(round(i0),0)
	i1 = max(round(i1),0)
	j0 = max(round(j0),0)
	j1 = max(round(j1),0)

	Variable nx = i1-i0+1
	Variable ny = j1-j0+1
	if (nx<1 || ny<1)									// nothing to read
		KillWaves/Z image
		return ""
	endif

	Duplicate/R=[i0,i1][j0,j1]/FREE image, ROIwave	// save the ROI I want
	Redimension/N=(nx,ny,-1) image
	image = ROIwave									// re-set image to ROI values
	WAVEClear ROIwave

	wnote = ReplaceNumberByKey("startx",wnote,i0,"=")
	wnote = ReplaceNumberByKey("endx",wnote,i1,"=")
	wnote = ReplaceNumberByKey("starty",wnote,j0,"=")
	wnote = ReplaceNumberByKey("endy",wnote,j1,"=")	
	Note/K image,wnote
	return GetWavesDataFolder(image,2)
End


//Function/T PilatusTiffReadHeader(tiffFileName,[extras])			// get info from inside the TIFF file, does not load an image
//	String tiffFileName
//	String extras
//	extras = SelectString(ParamIsDefault(extras),extras,"")		// not used here (yet)
//
Function/T PilatusTiffReadHeader(tiffFileName)			// get info from inside the TIFF file, does not load an image
	String tiffFileName

	DFREF dfrSave = GetDataFolderDFR()
	DFREF dfr = NewFreeDataFolder()
	SetDataFolder dfr
	ImageLoad/Z/O/Q/T=tiff/RTIO tiffFileName
	if (!V_flag)
		SetDataFolder dfrSave
		KillDataFolder/Z dfr
		return ""
	endif

	Variable i,val,N=ItemsInList(S_info,"\n"), j
	String item, list="waveClass=rawImage;detectorID=Pilatus 100K;"
	list= ReplaceStringByKey("file_name",list,tiffFileName,"=")
	for (i=0;i<N;i+=1)
		item = StringFromList(i, S_info,"\n")
		item = ReplaceString("# ",item,"#")

		strswitch(StringFromList(0,item," "))
			case "#Pixel_size":		// Pixel_size 172e-6 m x 172e-6 m
				j = strsearch(item, " ",strlen("#Pixel_size "),1)
				if (j<0)
					break
				endif
				item = TrimFrontBackWhiteSpace(item)
				item = ReplaceString(" x ", item[j+1,Inf],";")
				Variable dPixelX = str2num(StringFromList(0,item))
				val = convert2mm(StringFromList(1,StringFromList(0,item)," "))
				dPixelX *= numtype(val) ? 1 : val
				Variable dPixelY = str2num(StringFromList(1,item))
				val = convert2mm(StringFromList(1,StringFromList(1,item)," "))
				dPixelY *= numtype(val) ? 1 : val
				if (numtype(dPixelX+dPixelY)==0)
					list = ReplaceNumberByKey("dPixelX",list,dPixelX,"=")
					list = ReplaceNumberByKey("dPixelY",list,dPixelY,"=")
				endif
				break
			case "#Exposure_time":
				val = str2num(StringFromList(1,item," "))
				if (numtype(val)==0)
					list = ReplaceNumberByKey("Exposure_time",list,val,"=")
				endif
				break
			case "#Exposure_period":
				val = str2num(StringFromList(1,item," "))
				if (numtype(val)==0)
					list = ReplaceNumberByKey("Exposure_period",list,val,"=")
				endif
				break
			case "#Threshold_setting":
				val = str2num(StringFromList(1,item," "))
				if (numtype(val)==0)
					list = ReplaceNumberByKey("PilatusThreshold_eV",list,val,"=")
				endif
				break
		endswitch
	endfor

	Wave/T T_Tags=:Tag0:T_Tags
	Variable Nid = DimSize(T_Tags,0)
	Make/N=(Nid)/FREE IDs,index=NaN
	IDs = round(str2num(T_Tags[p][0]))
	IDs = numtype(IDs) ? Inf : IDs					// make all bad numbers inf
	MakeIndex IDs, index
	IndexSort index, IDs
	Make/FREE indexList={256,257,272,306}		// list of tiff tags to process
	Variable ii, xdim=NaN, ydim=NaN
	for (i=0;i<numpnts(indexList);i+=1)
		ii = index[BinarySearch(IDs,indexList[i])]
		switch(indexList[i])
			case 256:								// IMAGEWIDTH  --> xdim
				val = str2num(T_Tags[ii][4])
				if (numtype(val)==0 && val>0)
					xdim = val
					list = ReplaceNumberByKey("xdim",list,val,"=")
				endif
				break
			case 257:								// IMAGELENGTH  --> ydim
				val = str2num(T_Tags[ii][4])
				if (numtype(val)==0 && val>0)
					ydim = val
					list = ReplaceNumberByKey("ydim",list,val,"=")
				endif
				break
			case 272:								// MODEL--> detectorID
				item = TrimFrontBackWhiteSpace(T_Tags[ii][4])
				if (strlen(item)>0)
					list = ReplaceStringByKey("detectorID",list,item,"=")
				endif
				break
			case 306:								// DATETIME
				item = TrimFrontBackWhiteSpace(T_Tags[ii][4])
				item = ReplaceString(":",item,"-",0,2)
				item = ReplaceString(" ",item,"T",0,1)
				if (strlen(item)>0)
					list = ReplaceStringByKey("tiffDate",list,item,"=")
				endif
				break
		endswitch
	endfor
	WaveClear T_Tags
	SetDataFolder dfrSave
	KillDataFolder/Z dfr
	if (numtype(xdim+ydim)==0)				// These are default value, they can be reset later if we get good values
		Variable startx=0,endx=xdim-1,groupx=1
		Variable starty=0,endy=ydim-1,groupy=1
		list= ReplaceNumberByKey("xDimDet",list,xdim,"=")
		list= ReplaceNumberByKey("yDimDet",list,ydim,"=")
		list= ReplaceNumberByKey("startx",list,round(startx),"=")
		list= ReplaceNumberByKey("endx",list,round(endx),"=")
		list= ReplaceNumberByKey("groupx",list,round(groupx),"=")
		list= ReplaceNumberByKey("xdim",list,round(xdim),"=")
		list= ReplaceNumberByKey("starty",list,round(starty),"=")
		list= ReplaceNumberByKey("endy",list,round(endy),"=")
		list= ReplaceNumberByKey("groupy",list,round(groupy),"=")
		list= ReplaceNumberByKey("ydim",list,round(ydim),"=")
	endif
	return list
End
//
Static Function convert2mm(str)
	String str
	strswitch(str)
		case "m":
			return 1e3
		case "cm":
			return 10
		case "mm":
			return 1
		case "in":
		case "inch":
			return 25.4
		case "pica":
			return 25.4/6
		case "point":
			return 25.4/144
	endswitch
	return 1
End



// DEPRECATED			DEPRECATED			DEPRECATED			DEPRECATED			DEPRECATED
// in future, this is not the right way to do this
Function/S fetchPilatusInfoExtraSpec(fullName,[specFile,scanNum,pindex])	// only get information from the log file
	String fullName			// full path+name to tiff file
	String specFile			// path to the specFile
	Variable scanNum		// scanNumber in spec file
	Variable pindex			// point in scanNum

//	String list = ReadPilatusTiffHeader(fullName)		// get info from the tiff file itself (but do not load an image)
	String list = ""
	specFile = SelectString(ParamIsDefault(specFile),specFile,StrVarOrDefault("root:Packages:spec:specDefaultFile",""))
	scanNum = ParamIsDefault(scanNum) ? NaN : round(scanNum)
	scanNum = scanNum>0 ? scanNum : NaN
	pindex = ParamIsDefault(pindex) ? NaN : round(pindex)
	pindex = pindex>=0 ? pindex : NaN
	if (strlen(specFile)<1 || numtype(scanNum+pindex))
		return list
	endif

//	data1:EFRC:Nov 2012 33ID:raw:Cu111_Irradiated_1.spc
//	Cu111_Irradiated_1_S002_00000.tif

	String line="", buf=""
	Variable f
	Open/R/Z=1 f specFile
	if (V_flag)
		return list
	endif
	String fname = ParseFilePath(3,fullName,":",0,0)	// only the name part, no path, no extension
//	line=FindPFstart(f,fname)
	do
		buf += line
		FReadLine f, line
	while(strsearch(line,"#PSTATE",0)<0 && strlen(line))
	Close f
	buf = ReplaceString("\r",buf,"\n")
//	list=buf2Positions(buf)

	Variable i0=strsearch(buf, "\n#PD",0), i1=strsearch(buf, "\n",i0+1)
	if (i0>=0 && i1>i0)
		list = ReplaceStringByKey("file_time",list,buf[i0+5,i1-1],"=")
	endif

	i0=strsearch(buf, "\n#PT",0)
	i1=strsearch(buf, "\n",i0+1)
	if (i0>=0 && i1>i0)
		Variable exposure = str2num(buf[i0+5,i1-1])
		if (numtype(exposure)==0)
			list = ReplaceNumberByKey("exposure",list,exposure,"=")
		endif
	endif

	i0=strsearch(buf, "\n#V0",0)
	i1=strsearch(buf, "\n",i0+1)
	if (i0>=0 && i1>i0)
		Variable ringCurrent = str2num(buf[i0+5,i1-1])
		if (numtype(ringCurrent)==0)
			list = ReplaceNumberByKey("ringCurrent",list,ringCurrent,"=")
		endif
	endif

	i0=strsearch(buf, "\n#V1",0)
	i1=strsearch(buf, "\n",i0+1)
	if (i0>=0 && i1>i0)
		Variable keV = str2num(buf[i0+5,i1-1])
		if (numtype(keV)==0)
			list = ReplaceNumberByKey("keV",list,keV,"=")
		endif
	endif

	i0=strsearch(buf, "\n#PS Current Counts",0)
	i1=strsearch(buf, "\n",i0+1)
	if (i0>=0 && i1>i0)
		i0 = i1+1
		i1=strsearch(buf, "#",i0+1)-2
		String counters = buf[i0,i1], scalerList
//		scalerList = specTable2keyVals(counters)
		Variable j, value
		String name, item
		for (j=0;j<ItemsInList(scalerList);j+=1)
			item = StringFromList(j,scalerList)
			i0 = strsearch(item, "=",0)
			name = item[0,i0-1]
			value = str2num(item[i0+1,Inf])
			if (strlen(name)>0 && numtype(value)==0)
				list = ReplaceNumberByKey("scaler_"+name,list,value,"=")
			endif
		endfor
	endif

	if (strlen(list))
		i0=strsearch(buf, "\n#C",0)					// do not include this if nothing else found
		i1=strsearch(buf, "\n",i0+1)
		if (i0>=0 && i1>i0)
			list = ReplaceStringByKey("title",list,buf[i0+4,i1-1],"=")
		endif
	endif

	return list
End

//  ===================== End of NEW style with pilatus info in spec file ======================  //
//  ============================================================================  //



#else

//  ========================== Start of OLD style pilatuslog files ==========================  //
//  ============================================================================  //
// This section is for handling the old style of separate "*.pilatuslog" files
//

Function/S LoadPilatusTiffImage(fullName,[multiple,quiet])
	String fullName
	Variable multiple
	Variable quiet												// for quiet=1, if any problems just return
	multiple = ParamIsDefault(multiple) ? 1 : multiple
	multiple = numtype(multiple) ? 0 : !(!multiple)
	quiet = ParamIsDefault(quiet) ? 0 : quiet
	quiet = numtype(quiet) ? 0 : !(!quiet)

	Variable printIt=0
	GetFileFolderInfo/Q/Z=1 fullName							// check if file exists
	if (!V_isFile)
		if (quiet)
			return ""
		else
			Variable f
			Open /D/M="Pilatus tiff"/R/F=PilatusfileFilters f	// use /D to get full path name
			fullName = S_filename
			printIt = 1
		endif
	endif
	if (strlen(fullName)<1)									// no file name, quit
		return ""
	endif

	String fName = ParseFilePath(3,fullName,":",0,0), name	// only the name part, no path, no extension
	Variable i=strsearch(fName,".",0)
	name = CleanupName(fName[0,i-1],0)
	Variable scanNum,pntNum,eee,pindex
	sscanf fName[i,Inf],".%03d.%03d.%05d_%05d",scanNum,pntNum,eee,pindex
	if (V_flag!=4)
		if (!quiet)
			print "Cannot interpret Pilatus tif file name ='"+fName+"'"
		endif
		return ""
	endif
	String wName
	sprintf wName,"%s_%d_%d",name,scanNum,pntNum
	if (exists(wName)==1 || multiple)
		wName += "_"+num2istr(pindex)
	endif
	wName = CleanupName(wName,0)
	if (exists(wName)==1)
		return ""
	endif

	ImageLoad/Z/O/Q/T=tiff fullName
	if (!V_flag || V_numImages!=1)
		return ""
	endif
	Wave image = $StringFromList(0,S_waveNames)
	Rename image $wName
	if (printIt)
		printf "LoadPilatusTiffImage(\"%s\"",fullName
		if (!multiple)
			printf "multiple=%g",multiple
		endif
		printf ") // wave = '%s'\r",wName
	endif

	String funcName="PilatusFixBadPixels"+SelectString(exists("PilatusFixBadPixels")==6,"Proto","")
	FUNCREF PilatusFixBadPixelsProto fixPixels=$funcName
	fixPixels(image)

	SetScale/P x 0,1,"pixel", image
	SetScale/P y 0,1,"pixel", image
	SetScale d 0,0,"<counts>", image
	String wnote="waveClass=rawImage;"+fetchPilatusInfo(fullName)
	wnote = ReplaceNumberByKey("scanNum",wnote,scanNum,"=")
	wnote= ReplaceStringByKey("file_name",wnote,fullName,"=")
	Variable xdim=DimSize(image,0),ydim=DimSize(image,1)
	Variable startx=0,endx=xdim-1,groupx=1
	Variable starty=0,endy=ydim-1,groupy=1
	wnote= ReplaceNumberByKey("xDimDet",wnote,xdim,"=")
	wnote= ReplaceNumberByKey("yDimDet",wnote,ydim,"=")
	wnote= ReplaceNumberByKey("startx",wnote,round(startx),"=")
	wnote= ReplaceNumberByKey("endx",wnote,round(endx),"=")
	wnote= ReplaceNumberByKey("groupx",wnote,round(groupx),"=")
	wnote= ReplaceNumberByKey("xdim",wnote,round(xdim),"=")
	wnote= ReplaceNumberByKey("starty",wnote,round(starty),"=")
	wnote= ReplaceNumberByKey("endy",wnote,round(endy),"=")
	wnote= ReplaceNumberByKey("groupy",wnote,round(groupy),"=")
	wnote= ReplaceNumberByKey("ydim",wnote,round(ydim),"=")
	wnote= ReplaceStringByKey("file_name",wnote,fullName,"=")

	Note/K image, wnote
	return GetWavesDataFolder(image,2)
End


Function/T PilatusTiffReadROI(fileName,i0,i1,j0,j1)
	String fileName					// fully qualified name of file to open (will not prompt)
	Variable i0,i1,j0,j1			// pixel range of ROI (if i1 or j1<0 then use whole image)

	String wName = LoadPilatusTiffImage(fileName)
	Wave image=$wName
	String wnote = note(image)

	Variable ydim, xdim									// x,y size of whole array
	xdim = DimSize(image,0)
	ydim = DimSize(image,1)
	i1 = (i1<1) ? xdim-1 : i1							// -1 flags use whole range
	j1 = (j1<1) ? ydim-1 : j1

	i0 = max(round(i0),0)
	i1 = max(round(i1),0)
	j0 = max(round(j0),0)
	j1 = max(round(j1),0)

	Variable nx = i1-i0+1
	Variable ny = j1-j0+1
	if (nx<1 || ny<1)									// nothing to read
		KillWaves/Z image
		return ""
	endif

	Duplicate/R=[i0,i1][j0,j1]/FREE image, ROIwave	// save the ROI I want
	Redimension/N=(nx,ny,-1) image
	image = ROIwave									// re-set image to ROI values
	WAVEClear ROIwave

	wnote = ReplaceNumberByKey("startx",wnote,i0,"=")
	wnote = ReplaceNumberByKey("endx",wnote,i1,"=")
	wnote = ReplaceNumberByKey("starty",wnote,j0,"=")
	wnote = ReplaceNumberByKey("endy",wnote,j1,"=")	
	Note/K image,wnote
	return GetWavesDataFolder(image,2)
End


Function/S fetchPilatusInfo(fullName,[logFile])		// only get information from the log file
	String fullName			// full path+name to file
	String logFile			// path to the logfile

	if (ParamIsDefault(logFile))
		logFile = getPilatusLogFileName(fullName)
	endif
	if (strlen(logFile)<1)
		return ""
	endif

	String line, buf=""
	Variable f
	Open/R/Z=1 f logFile
	if (V_flag)
		return ""
	endif
	String fname = ParseFilePath(3,fullName,":",0,0)	// only the name part, no path, no extension
	line=FindPFstart(f,fname)
	do
		buf += line
		FReadLine f, line
	while(strsearch(line,"#PSTATE",0)<0 && strlen(line))
	Close f
	buf = ReplaceString("\r",buf,"\n")
	String list=buf2Positions(buf)

	Variable i0=strsearch(buf, "\n#PD",0), i1=strsearch(buf, "\n",i0+1)
	if (i0>=0 && i1>i0)
		list = ReplaceStringByKey("file_time",list,buf[i0+5,i1-1],"=")
	endif

	i0=strsearch(buf, "\n#PT",0)
	i1=strsearch(buf, "\n",i0+1)
	if (i0>=0 && i1>i0)
		Variable exposure = str2num(buf[i0+5,i1-1])
		if (numtype(exposure)==0)
			list = ReplaceNumberByKey("exposure",list,exposure,"=")
		endif
	endif

	i0=strsearch(buf, "\n#C",0)
	i1=strsearch(buf, "\n",i0+1)
	if (i0>=0 && i1>i0)
		list = ReplaceStringByKey("title",list,buf[i0+4,i1-1],"=")
	endif

	i0=strsearch(buf, "\n#V0",0)
	i1=strsearch(buf, "\n",i0+1)
	if (i0>=0 && i1>i0)
		Variable ringCurrent = str2num(buf[i0+5,i1-1])
		if (numtype(ringCurrent)==0)
			list = ReplaceNumberByKey("ringCurrent",list,ringCurrent,"=")
		endif
	endif

	i0=strsearch(buf, "\n#V1",0)
	i1=strsearch(buf, "\n",i0+1)
	if (i0>=0 && i1>i0)
		Variable keV = str2num(buf[i0+5,i1-1])
		if (numtype(keV)==0)
			list = ReplaceNumberByKey("keV",list,keV,"=")
		endif
	endif

	i0=strsearch(buf, "\n#PS Current Counts",0)
	i1=strsearch(buf, "\n",i0+1)
	if (i0>=0 && i1>i0)
		i0 = i1+1
		i1=strsearch(buf, "#",i0+1)-2
		String counters = buf[i0,i1], scalerList
		scalerList = specTable2keyVals(counters)
		Variable j, value
		String name, item
		for (j=0;j<ItemsInList(scalerList);j+=1)
			item = StringFromList(j,scalerList)
			i0 = strsearch(item, "=",0)
			name = item[0,i0-1]
			value = str2num(item[i0+1,Inf])
			if (strlen(name)>0 && numtype(value)==0)
				list = ReplaceNumberByKey("scaler_"+name,list,value,"=")
			endif
		endfor
	endif

///////////////
	String PilatusID="Pilatus 100K"
	if (strlen(PilatusID))
		list= ReplaceStringByKey("detectorID",list,PilatusID,"=")
	endif
///////////////

	return list
End
//
Static Function/T buf2Positions(buf)
	String buf

	Variable i0=strsearch(buf, "#PP Current Positions",0), i1=strsearch(buf, "\n#",i0+1),j
	if (i0<0 || i1<0)
		return ""
	endif
	i0= strsearch(buf, "\n",i0+1)+1		// skip a new-line
	buf = buf[i0,i1]
	return specTable2keyVals(buf)
End
//
Static Function/T specTable2keyVals(buf)
	String buf

	String names,values,list=""
	Variable i0=0, i1=1, j
	do
		names = StringFromList(i0,buf,"\n")
		values = StringFromList(i1,buf,"\n")
		do
			names = ReplaceString("  ",names," ")
		while(strsearch(names,"  ",0)>=0)
		do
			values = ReplaceString("  ",values," ")
		while(strsearch(values,"  ",0)>=0)
		names = TrimFrontBackWhiteSpace(names)
		values = TrimFrontBackWhiteSpace(values)
		names = ReplaceString(" ",names,";")
		values = ReplaceString(" ",values,";")

		for (j=0;j<ItemsInList(names);j+=1)
			list = ReplaceStringByKey(StringFromList(j,names),list,StringFromList(j,values),"=")
		endfor
		i0 += 2
		i1 += 2
	while(strlen(values)>1)
	return list
End
//
Static Function/T FindPFstart(fileVar,fname)		// search for start of image info
	Variable fileVar		// file ref number
	String fname

	if (strlen(fname)<0)
		return ""
	endif

	String search="#PF PILATUS file name: "+fname

	Variable bufSize=5e5
	Variable Nread=strlen(search)+bufSize
	String buffer=PadString("",Nread,0)		// set up large buffer for reading data file

	FStatus fileVar
	Variable fpos0=V_filePos
	Variable more=1,i,j=-bufSize

	do
		j += bufSize
		FSetPos fileVar, (j+fpos0)
		FStatus fileVar
		if ((V_logEOF-V_filePos)<Nread)		// not enough data to fill buffer
			buffer = ""
			buffer = PadString(buffer,V_logEOF-V_filePos, 0)
			more = 0
		endif
		FBinRead fileVar, buffer				// read next buffer full of data
		i = strsearch(buffer, search, 0)		// check
	while ((i<0) && more)						// continue if not found and more data in file

	if (i<0)										// not found
		return ""
	endif

	FSetPos fileVar, (i+j+fpos0)
	String line=" "
	FReadLine fileVar, line
	return line
End
//
Static Function/S getPilatusLogFileName(fullName)
	String fullName								// full path name to image file
	String logFile=""							// path to the logfile
	if (strlen(fullName)<1)					// none passed, use the default location, or ask user
		logFile=StrVarOrDefault("root:Packages:pilatus:logFile","")			// path to the logfile
		if (!isValidPilatusLogFile(logFile))
			String PilatusLogFilters = "Pilatus Log Files (*.pilatuslog):.pilatuslog;All Files:.*;"
			Variable f
			Open/R/F=PilatusLogFilters/M="Pilatus Log File"/P=home/Z=2 f
			logFile = S_fileName
		endif
	else											// get log file from image file name & location
		String fName = ParseFilePath(3,fullName,":",0,0)	// only the name part, no path, no extension
		String path = ParseFilePath(1,fullName,":",1,0)	// path to the folder with image
		Variable i=strsearch(fName,".",0)
		logFile = path+fName[0,i-1]+".pilatuslog"
	endif

	if (isValidPilatusLogFile(logFile))			// check that logFile is valid, and update default
		NewDataFolder/O root:Packages
		NewDataFolder/O root:Packages:pilatus
		String/G root:Packages:pilatus:logFile=logFile
	else
		logFile = ""
	endif
	return logFile
End
//
Static Function isValidPilatusLogFile(logFile)
	String logFile

	GetFileFolderInfo/Q/Z=1 logFile
	if (!V_isFile)
		return 0
	endif
	String buf=PadString("",1000,0x20 )
	Variable f
	Open/R/Z=1 f as logFile
	if (V_flag)
		return 0
	endif
	FBinRead f,buf
	Close(f)
	return (strsearch(buf,"\n#PF ",0)>=0)
End

//Function test()
//	String fullName="data1:EFRC:July 2012 33ID:tiffs:Cu111_Ni:Cu111_Ni.065.013.46913_01619.tif"
//	fullName="data1:EFRC:July 2012 33ID:tiffs:Cu111_He:Cu111_He.017.000.81017_00536.tif"
//
//	KillWaves/Z Cu111_He_17_0_536
//	Wave image = $LoadPilatusTiffImage(fullName)
//
//	print "keV = ",NumberByKey("keV",note(image),"=")
//	print "ringCurrent = ",NumberByKey("ringCurrent",note(image),"=")
//	print "tth = ",NumberByKey("tth",note(image),"=")
//
//End

//		#PF PILATUS file name: Cu111_He.017.000.81017_00536.tif
//		#PD Fri Jul 27 16:14:55 2012
//		#PT 5  (seconds)
//		#G0 0 0 0 0 0 1 0 0 0 0 0 0 50 0 0.1 0 68 68 50 -1 1 1 3.13542 3.13542 0 463.6 838.8 0
//		#G1 1.54 1.54 1.54 90 90 90 4.079990459 4.079990459 4.079990459 90 90 90 1 0 0 0 1 0 60 30 0 0 0 0 60 30 0 -90 0 0 1.54 1.54
//		#G3 4.079990459 2.19415563e-16 2.498273628e-16 -3.041179982e-17 -4.079990459 2.498273628e-16 0 0 -4.079990459
//		#G4 -0.008316799033 -0.008238379486 0 1.23984244 0 0 0 -90 0 0 0 0 0 0 0 0 -180 -180 -180 -180 -180 -180 -180 -180 -180 0
//		#PQ -0.0083168 -0.00823838 0
//		#PP Current Positions
//		      tth        th       chi       phi       kth       kap      kphi     dummy
//		  -0.5400 -179.9986    0.0000  135.0000  180.0014    0.0000  135.0000    0.0000
//		   slitwt    slitwl    slitwb    slitwr     dcmth    slitmt    slitml    slitmb
//		   0.8000    1.2000    0.8000    0.5000   11.4027    6.5002    1.9999   -3.4999
//		   slitmr        mu        nu       ath      atth      nty1      nty2      nty3
//		   2.0002   17.2885   34.8600 -135.9866    0.0050   13.6363   11.7993   10.1953
//		    sampx     sampy     sampz    post1V      camz     wirex       osx      camy
//		   0.2598    0.0000    0.0033   16.2498   18.2999   11.6051   50.1890   -3.0002
//		     camx       s2x       s2y     s0hap     s0vce     s0vap     s0hce      s00t
//		  -3.0007    0.0006    0.0006    0.0999    0.0000    0.3001   -0.0000    1.5000
//		     s00l      s00b      s00r     s1hce     s1hap     s1vce     s1vap       s1t
//		   2.0000    1.5000    2.0000    0.0000   20.0000    0.0000   20.0000   10.0000
//		      s1b       s1l       s1r     s2hap     s2hce     s2vap     s2vce       s2l
//		  10.0000   10.0000   10.0000    0.0000    0.0000    0.0000    0.0000    0.0000
//		      s2r       s2t       s2b     s3hap     s3hce     s3vap     s3vce       s3l
//		   0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000
//		      s3r       s3t       s3b      m1y1      m1y2      m2y1      m2y2    kohchi
//		   0.0000    0.0000    0.0000   17.8500   21.0500   28.8680   31.9300   -4.2167
//		    ksamx     ksamy     ksamz     atto1     atto2     atto3        hx        hy
//		   0.0400    0.0000    2.5000    0.0000    0.0000    0.0000    0.0000    0.0000
//		       hz        hu        hv        hw       zpx       zpy       zpz     tubez
//		   0.0000    0.0000    0.0000    0.0000   -4.7700   -9.7200    2.6000    0.0000
//		   s01hap    s01hce    s01vap    s01vce      s01l      s01r      s01t      s01b
//		   0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000
//		      ty1       ty2       ty3       tx1       tx2       hbu       hbd       hpu
//		  48.0466   42.4663   45.6355   58.1804   13.9709    0.0000    0.0000    1.3484
//		      hpd    htrans       vbu       vbd       vpu       vpd    vtrans     kbvcu
//		   2.4716   14.5000    0.0000    0.0000    0.6884    1.8716    4.9999    0.0000
//		    kbvel      kbvz     kbvth     kbhcu     kbhel      kbhz     kbhth      hexx
//		   0.0000    1.2800    5.1001    0.0000    0.0000    1.9100    3.6002    0.8000
//		     hexy      hexz      hexu      hexv      hexw 
//		  -2.7000   -1.6200   -0.8670   -0.1000    0.0000 
//		#PX
//		#C Cu(111) He implant
//		#V0 102.238 18 1 4 1 1 1  982.319
//		#V1 10.0002 10.1487 -0.000307083
//		#V2 22.0182 10.1487 1 -0.000569699 -0.000307083
//		#V3 0 0 0 0 1 10.0002 1.23981 11.4027 11.4025 0.2 1 1 -0.0164527 -10 10 3.57265
//		#V4 0 -0.03487 0
//		#V5 4.00001 19.45 3
//		#V6 3.82743 30.399 3
//		#V7 1 -7 1.0e+10 0 7 0 0 -7e-10
//		#V8 1 -0.69 -0.692011 1 -0.015 -0.015575
//		#V9 0.4852 0 0.0348579 1 1 0 0 0 0 0 0
//		#V10 10 10 10 10 20 20 0 0
//		#V11 0 0 0 0 0 0 0 0
//		#V12 0 0 0 0 0 0 0 0
//		#V13 0 0 0 0 0 0 0 0
//		#V14 0.002442 10000000 7 5 -1.27e-10 1.01042 2.442e-10
//		#V15 0.0276188 1000000 7 0 5.1e-11 824.854 2.76188e-08
//		#V16 0.041514 100000 1 5 1.8e-10 1.01042 4.1514e-07
//		#V17 0.041514 -0.021978 0.01221 0.002442 0.510379 0.041514 0.00732601 -0.026862
//		#V18 -0.105006 -0.0708181 -0.114774 -0.197802 -0.01221 -0.017094 -0.01221 -0.01221
//		#V19 0 0 10 0
//		#V20 1000000 0 0 1.01042 1.01042 824.854 0 0
//		#V21 0 0 0 0 0 0 951.042 1.01042
//		#V22 96 97 100
//		#PS Current Counts
//		      sec        I0       det      harm       IC1   scanbar   corrint    ccdtot
//		        5    250831         0         0      4123 1.65908e+06         0    190213
//		   piltot   ccdsca1   ccdsca2   ccdsca3    ccdmax   filters     trans   corrdet
//		7.8206e+07 1.81906e+07    643317    138975    732836       300 0.0348579         0
//		   Energy 
//		       10 
//		#PEPOCH 181079
//		#PE
//		#PROI MinX SizeX MinY SizeY
//		     200       80       51       70
//		     234       10       96       10
//		     258       10      109       10
//		     240       12       30       12
//		#PROIE
//		
//		#PSTAT Total Net Min Max Mean Sigma CentX YCent SigmaX SigmaY SigmaXY
//		      13   732836 18190566 18190566     3248    23737       39       31       15       15        0
//		    1621    97185   643317   643317     6433    14180        4        2        1        2        0
//		    1295     1497   138975   138975     1389       39        4        4        2        2        0
//		    1195     1948   190213   190213     1320      107        5        5        3        3        0
//		      13   732836 78206016 78206016      823     5804      239       93      100       50        0
//		#PSTATE
//  ========================== End of OLD style pilatuslog files ===========================  //
//  ============================================================================  //
#endif



#ifdef EnablePilatusUnWrapping

//  ============================================================================  //
//  ============================== Start of Unwrap Image ==============================  //

Function UnWrapImage(image)
	Wave image

	Variable unWrapped=0
	Variable unwrap_one=1
	do
		unwrap_one = UnWrapImageStep(image)
		unWrapped = unWrapped || unwrap_one
	while (unwrap_one)
	return unWrapped
End
//
Static Function UnWrapImageStep(image,[printIt])
	Wave image
	Variable printIt
	printIt = ParamIsDefault(printIt) ? 0 : (!(!printIt))

	Variable Nx=DimSize(image,0), Ny=DimSize(image,1), i,j
	Duplicate/FREE image wraps
	wraps = -1
	wraps[0,2][] = 0				// set edges to 0 (border 3 wide
	wraps[][0,2] = 0
	wraps[Nx-3,Nx-1][] = 0
	wraps[][Ny-3,Ny-1] = 0

	Variable i2,i1,i0, wrap
	for (j=3;j<Ny-3;j+=1)				// first from the left
		for (i=3;i<Nx-3;i+=1)
			wrap = wrapForOnePixel(image,i,j)
			wraps[i][j] = wrap
			if (wrap)
				break
			endif
		endfor
	endfor
	for (j=3;j<Ny-3;j+=1)				// again from the rigtht
		if (wraps[Nx-4][j]>=0)
			continue						// already did this row
		endif
		for (i=Nx-4;i>=3;i-=1)
			wrap = wrapForOnePixel(image,i,j)
			wraps[i][j] = wrap
			if (wrap)
				break
			endif
		endfor
	endfor

	WaveStats/Q/M=1 wraps
	Variable unwrap_flag = (V_max==1)
	if (unwrap_flag)						// Crawl around the peak, and fill in
		if (printIt)
			DoUpdate
		endif
		i = V_maxRowLoc
		j = V_maxColLoc
		if (printIt)
			print "istart =",i,"  jstart = ",j,"   starting point for going around","   wraps[i][j]=",wraps[i][j]
		endif

		Make/N=(8,2)/FREE boxij
		boxij[0][0]= {-1,-1,0,1,1,1,0,-1}
		boxij[0][1]= {0,1,1,1,0,-1,-1,-1}

		Duplicate/FREE wraps, done
		done = 0
		Variable mlast,mindex,m=0
		do
			wraps[i][j] = wrapForOnePixel(image,i,j,printIt=printIt)
			done[i][j] = 1
			if (printIt)
				Debugger
				Duplicate/O wraps, wrapsView
				Cursor/P/I/W=Graph6 A wrapsView i,j
				DoUpdate
			endif
			for (mindex=0;mindex<8;mindex+=1)
				m = mod(mindex+mlast,8)
				if (wraps[i+boxij[m][0]][j+boxij[m][1]] != 0 && !done[i+boxij[m][0]][j+boxij[m][1]])
					i += boxij[m][0]
					j += boxij[m][1]
					mlast = m-1
					break
				endif
			endfor
		while (atEdge(wraps,i,j))
	endif

	wraps = wraps<0 ? 1 : wraps
	Variable twenty = 2^20						// Pilatus is 20 bits deep
	image += wraps*twenty
//	Duplicate/O wraps, wrapsView
	return unwrap_flag
End
//
Static Function atEdge(wraps,i,j)
	Wave wraps
	Variable i,j
	Make/N=(3,3)/FREE m33
	m33 = wraps[i+p-1][j+q-1]
	Redimension/N=9 m33
	Variable aZero = (m33[0]*m33[1]*m33[2]*m33[3]*m33[4]*m33[5]*m33[6]*m33[7]*m33[8]==0)
	Variable aNeg = WaveMin(m33) < 0
	return aZero && aNeg
End
//
Static Function wrapForOnePixel(image,px,py,[printIt])
	Wave image
	Variable px,py
	Variable printIt
	printIt = ParamIsDefault(printIt) ? 0 : printIt

	Variable ij00=image[px][py]
	Variable i1m=image[px-1][py], i2m=image[px-2][py], i3m=image[px-3][py]
	Variable i1=image[px+1][py], i2=image[px+2][py], i3=image[px+3][py]
	Variable j1m=image[px][py-1], j2m=image[px][py-2], j3m=image[px][py-3]
	Variable j1=image[px][py+1], j2=image[px][py+2], j3=image[px][py+3]

	Variable del=50e3

	if (i2m+del>i3m && i1m+del>i2m && ij00+del>i1m  &&  j2m+del>j3m && j1m+del>j2m && ij00+del>j1m)	// at top left, just increasing
		if (printIt)
			printf "[%g, %g], at top left, just increasing\r",px,py
		endif
		return 0
	elseif (i2m+del>i3m && i1m+del>i2m && ij00+del>i1m  &&  j2+del>j3 && j1+del>j2 && ij00+del>j1)		// at lower left, just increasing
		if (printIt)
			printf "[%g, %g], at lower left, just increasing\r",px,py
		endif
		return 0
	elseif (i2+del>i3 && i1+del>i2 && ij00+del>i1  &&  j2m+del>j3m && j1m+del>j2m && ij00+del>j1m)		// at top right, just increasing
		if (printIt)
			printf "[%g, %g], at top right, just increasing\r",px,py
		endif
		return 0
	elseif (i2+del>i3 && i1+del>i2 && ij00+del>i1  &&  j2+del>j3 && j1+del>j2 && ij00+del>j1)				// at lower right, just increasing
		if (printIt)
			printf "[%g, %g], at lower right, just increasing\r",px,py
		endif
		return 0
	endif

	if (printIt)
		printf "[%g, %g], wrap =1\r",px,py
	endif
	return 1

End

//  =============================== End of Unwrap Imag ==============================  //
//  ============================================================================  //

#endif
