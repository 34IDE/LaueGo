#pragma rtGlobals=1		// Use modern global access method.
#pragma version = 1.17
#pragma ModuleName=TiffProc
#include "ImageDisplayScaling", version>=1.98

//
// Routines for reading in and looking at TIFF files. Eliot Specht, ORNL, based on WinView.ipf code by Jon Tischler, Oak Ridge National Lab
// TischlerJZ@ornl.gov
//
// It should be made more general, to handle other image files

//Menu "BackReflection",  dynamic
//	MenuItemIfWaveClassExists("Display image plot with buttons...","tif*","DIMS:2"), Graph_imageMake($"",0)
//	help={"Display an image and put buttons on the plot for listing spots"}
//	MarqueeTiffMenuItem("Center a Gaussian and Add to pkList"), AddGaussianToList()
//	MarqueeTiffMenuItem("Report Center of a Gaussian"), GaussianCenter(1)
//	help={"fit a gaussian to an ROI an image, and print results.  The region is selected with a marquee in an image plot"}
//	MarqueeTiffMenuItem("Report Center of Mass of a ROI"), GetCenterOfMass(NaN,NaN,NaN,NaN)
//	help={"Compute the center of mass of a region on an image.  The region is selected with a marquee in an image plot"}
//	"-"
//	pkLIstTiffMenuItem("Fit List of Gaussians..."), FitListOfGaussians("","pkList",-1)
//	GcoefTiffMenuItem("Reset pkList from Gaussians"), pkList = Reset_pkList_from_Gaussians(Gauss_coef)
//	"-"
//	"Load TIFF File...", LoadWinViewFile("")
//	help={"Load a TIFF Image from the file, and then display the image"}
////	MenuItemIfWaveClassExists("Get TIFF Info","speImage*","DIMS:2"), WinViewInfo($"","")
////	help={"Get the other information stored with the winview image"}
////	"Show WinView Procedures", DisplayProcedure/B=WinViewProcedures "CenterOfMass"
//	     MenuItemIfWaveClassExists("   Find Z range of image","*","DIMS:2"), getZrange($"",NaN)
//	help={"for an image, find the z range for the specified % range"}
//End

// Nov  15, 2010,  with version 1.13, added the extras string to the file loaders
// Jun  20, 2014,  with version 1.16, added the extras string to TiffReadHeader()

Static Constant TAG_IMAGEWIDTH=256, TAG_IMAGELENGTH=257, TAG_BITSPERSAMPLE=258
Static Constant TAG_MODEL=272, TAG_XRESOLUTION=282, TAG_YRESOLUTION=283, TAG_DATETIME=306
Static Constant TAG_RESOLUTIONUNIT=296, TAG_SOFTWARE=305



Menu "Load Waves"
	"   Load Tiff File...", LoadTiffFile("")
	help={"Load a Tiff Image from the file, and then display the image"}
End


Function/S MarqueeTiffMenuItem(item)
	String item
	if (strlen(WinList("*","","WIN:1"))<1)
		return "("+item
	endif
	Variable V_flag
	GetMarquee/Z								// These menu items requre a Marquee present
	if (V_flag==0)
		return "("+item
	endif
	Wave image = ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";")))
	if (!WaveExists(image))
		return "("+item
	endif
	String class = StringByKey("waveClass", note(image),"=")
	if (!stringmatch(class,"spe*"))				// is not an editable mas, disable (only works for spe images
		item = "("+item
		// item = "("
	endif
	return item
End
//
//Function/S GcoefTiffMenuItem(item)
//	String item
//	if (WaveExists(Gauss_coef))
//		return item
//	endif
//	return "("+item
//End
//
//Function/S pkLIstTiffMenuItem(item)
//	String item
//	if (strlen(StrVarOrDefault("pkList","")))
//		return item
//	endif
//	return "("+item
//End


Function/T TiffSaveImage(fileName,image)
	String fileName					// fully qualified name of file to open (will not prompt)
	Wave image
	if (!WaveExists(image))
		print "ERROR, TiffSaveImage(), the image wave does not exist"
		return ""
	elseif (WaveDims(image)!=2)
		printf "ERROR,TiffSaveImage(), only writes 2-d waves, '%s' is %d-d\r",WaveDims(image),NameOfWave(image)
		return ""
	endif
	Variable TiffType = WaveType(image)

	if (!(TiffType && 0x18))
		printf "ERROR, WinViewSaveImage(), cannot write TIFF files for WaveType(%s) = 0x%x\r",NameOfWave(image),TiffType
		return ""
	endif

#if (exists("Init_microGeo")==6)				// test for presence of micro geometry
	if (exists("Init_microGeo")==6)				// test for presence of micro geometry
		STRUCT microGeometry g
		FillGeometryStructDefault(g)				//fill the geometry structure with current values
		Variable id=detectorNumFromID("Photonic Sciences")
		id = id<0 ? 0 : id								// when no id, just use the first one, zero
		Variable dx = (g.d[id].sizeX)/(g.d[id].Nx)	// pixel size (micron)
		Variable dy = (g.d[id].sizeY)/(g.d[id].Ny)
		Wave/T tagWave = MakeTagsWave(dx,dy)
	endif
#endif

	if (WaveExists(tagWave))
		if (TiffType & 0x10)						// 16 bit image
			ImageSave /D=16/IGOR/T="tiff"/U/WT=tagWave image as fileName
		elseif (TiffType & 0x08)					// 8 bit image
			ImageSave /D=40/IGOR/T="tiff"/U/WT=tagWave image as fileName
		endif
	else
		if (TiffType & 0x10)						// 16 bit image
			ImageSave /D=16/IGOR/T="tiff"/U image as fileName
		elseif (TiffType & 0x08)					// 8 bit image
			ImageSave /D=40/IGOR/T="tiff"/U image as fileName
		else
			return ""
		endif
	endif
	return fileName
End
//
Static Function/WAVE MakeTagsWave(dx,dy)
	Variable dx,dy			// pixel size (micron)

	Make/N=(5,5)/FREE/T T_Tags=""

	String str = num2str(1e4/dx)					// x resolution (pixels/cm)
	T_Tags[0][0] = num2istr(TAG_XRESOLUTION)
	T_Tags[0][2] = "5"
	T_Tags[0][3] = "1"
	T_Tags[0][4] = str

	str = num2str(1e4/dy)							// y resolution (pixels/cm)
	T_Tags[1][0] = num2istr(TAG_YRESOLUTION)
	T_Tags[1][2] = "5"
	T_Tags[1][3] = "1"
	T_Tags[1][4] = str

	T_Tags[2][0] = num2istr(TAG_RESOLUTIONUNIT)	// resolution units, "pixels/cm"
	T_Tags[2][2] = "3"
	T_Tags[2][3] = "1"
	T_Tags[2][4] = "3"								// 3 means "cm"

	str = IgorInfo(1)+", "+	"Igor "+StringByKey("IGORFILEVERSION",IgorInfo(3))+", "+StringByKey("OS",IgorInfo(3))
	T_Tags[3][0] = num2istr(TAG_SOFTWARE)	// name of software, here it is Igor
	T_Tags[3][2] = "2"
	T_Tags[3][3] = num2istr(strlen(str))
	T_Tags[3][4] = str

	Variable seconds=datetime
	str = Secs2Date(DateTime,-2,":")+" " + Secs2Time(seconds,3)
	T_Tags[4][0] = num2istr(TAG_DATETIME)		// date & time file is written
	T_Tags[4][2] = "2"
	T_Tags[4][3] = num2istr(strlen(str))
	T_Tags[4][4] = str

	return T_Tags
End


Function/S LoadTiffFile(fName,[extras])
	String fName											// fully qualified name of file to open
	String extras											// not doing anything here yet
	extras = SelectString(ParamIsDefault(extras),extras,"")

	Variable refNum
//	if (strlen((OnlyTiffFileName(fName)))<1)		// call dialog if no file name passed
	if (strlen((ParseFilePath(0,fName,":",1,0)))<1)	// call dialog if no file name passed
		Open /D/M="image file"/R/T="????" refNum	// use /D to get full path name
		fName = S_filename
	endif
	if (strlen(fName)<1)									// no file name, quit
		return ""
	endif

	String wName =TiffReadROI(fName,0,-1,0,-1)	// load file into wName
	if (strlen(wName)<1)
		return ""
	endif
	Wave image = $wName
	
	if (ItemsInList(GetRTStackInfo(0))<=1)
		String wnote = note(image)
		Variable xdim=DimSize(image,0)
		Variable ydim=DimSize(image,1)
//		String bkgFile = StringByKey("bkgFile", wnote,"=")
		printf "for file '"+fName+"'"
//		if (strlen(bkgFile)>0)
//			printf ",             background file = '%s'",  bkgFile
//		endif
		printf "\r"
		printf "total length = %d x %d  = %d points\r", xdim,ydim,xdim*ydim
//		print "number type is  '"+TiffFileTypeString(NumberByKey("numType", wnote,"="))+"'"
		print "Created a 2-d wave    '"+wName+"'"
		DoAlert 1, "Display this image"
		if (V_Flag==1)
			NewImageGraph(image)
//			Graph_imageMake(image,NaN)
		endif
	endif
	
//	print "Tif/TifLoadFile: Note(image) = ", note(image)
	return GetWavesDataFolder(image,2)
End


Function/T TiffReadROI(fileName,i0,i1,j0,j1,[extras])
	String fileName					// fully qualified name of file to open (will not prompt)
	Variable i0,i1,j0,j1			// pixel range of ROI (if i1 or j1<0 then use whole image)
	String extras											// not doing anything here yet
	extras = SelectString(ParamIsDefault(extras),extras,"")

	String wName = TIffLoadROI(fileName,i0,i1,j0,j1)	// name of ROI read in
	if (strlen(wName)<1)
		return ""												// nothing read in
	endif
	Wave image = $wName
	SetScale/P x 0,1,"pixel", image
	SetScale/P y 0,1,"pixel", image

	String wnote = note(image)
	if (i0!=0 || j0!=0)										// sub-region of image does not start at (0,0)
		wnote = ReplaceNumberByKey("startxRead",wnote,i0,"=")
		wnote = ReplaceNumberByKey("startyRead",wnote,j0,"=")
	endif
	wnote = ReplaceStringByKey("waveClass",wnote,"rawImage","=")
	wnote = ReplaceStringByKey("DetectorID",wnote,"Photonic Sciences","=")
	Note /K image, wnote
	return GetWavesDataFolder(image,2)
End


Function/T TiffLoadROI(fileName,i0,i1,j0,j1,[extras])
	String fileName				// fully qualified name of file to open (will not prompt)
	Variable i0,i1,j0,j1			// pixel range of ROI
	String extras											// not doing anything here yet
	extras = SelectString(ParamIsDefault(extras),extras,"")

	String wnote=""
	if (strlen(extras))
		wnote = TiffReadHeader(fileName, extras=extras)
	else
		wnote = TiffReadHeader(fileName)
	endif
	if (strlen(wnote)<1)
		return ""
	endif

	ImageLoad/T=tiff/Q fileName
	if (V_flag != 1)
		return ""											// could not open file
	endif
	String inWaveName=StringFromList(0,S_waveNames)

	String name=ParseFilePath(3,fileName,":",0,0)	// rename image based on the file name
	name = CleanupName(name,0)						// wave name based on the file name
	if (exists(name))										// if wave already exists, create unique name
		name = name+"_"
		name = UniqueName(name, 1, 1)
	endif
	String wName=ParseFilePath(1,inWaveName,":",1,0)+name
	Rename $inWaveName $wName
	Wave wav = $wName

	Variable ydim, xdim									// x,y size of whole array
	xdim = DimSize(wav,0)
	ydim = DimSize(wav,1)
	i1 = (i1<1) ? xdim-1 : i1							// -1 flags use whole range
	j1 = (j1<1) ? ydim-1 : j1

	i0 = max(round(i0),0)
	i1 = max(round(i1),0)
	j0 = max(round(j0),0)
	j1 = max(round(j1),0)
	Variable nx = i1-i0+1
	Variable ny = j1-j0+1
	if (nx<1 || ny<1)										// nothing to read
		return ""
	endif

	Duplicate/R=[i0,i1][j0,j1]/FREE wav, ROIwave
	Duplicate/O ROIwave, wav

	String flatFile="", bkgFile=""
	wnote = ReplaceNumberByKey("startx",wnote,i0,"=")
	wnote = ReplaceNumberByKey("endx",wnote,i1,"=")
	wnote = ReplaceNumberByKey("starty",wnote,j0,"=")
	wnote = ReplaceNumberByKey("endy",wnote,j1,"=")	
	Note/K wav,wnote
	return GetWavesDataFolder(wav,2)
End


Function/T TiffReadHeader(fName,[extras])
	String fName					// fully qualified name of file to open (will not prompt)
	String extras					// optional switches (only supports EscanOnly in this routine)
	extras = SelectString(ParamIsDefault(extras),extras,"")
	if (strlen(fName)<1)
		return ""
	endif

	String fldrSav= GetDataFolder(1)
	String dataFolderSpec=UniqueName("TAGfolder",11,0)
	NewDataFolder/S $dataFolderSpec

	ImageLoad/Q/T=tiff/RTIO/Z fName
	if (V_Flag!=1)
		return ""						// could not open file
	endif
	Wave tagsIN=$(":Tag0:T_Tags")
	if (WaveExists(tagsIN))
		Duplicate/FREE tagsIN, tags
	endif
	SetDataFolder fldrSav
	KillDataFolder/Z $dataFolderSpec
	if (!WaveExists(tags))
		return ""
	endif

	String wnote=""
	wnote = ReplaceStringByKey("waveClass",wnote,"rawImage","=")
	wnote = ReplaceStringByKey("imageFilePath",wnote,ParseFilePath(1,fName,":",1,0),"=")
	wnote = ReplaceStringByKey("imageFileName",wnote,ParseFilePath(0,fName,":",1,0),"=")
	wnote= ReplaceStringByKey("detectorID",wnote,"Photonic Sciences","=")

	// detector pixels & ROI
	Variable xdim,ydim
	xdim = str2num(strFromTagWave(tags,TAG_IMAGEWIDTH))
	ydim = str2num(strFromTagWave(tags,TAG_IMAGELENGTH))
	if (numtype(xdim+ydim)==0)
		Variable startx=0,endx,groupx=1
		Variable starty=0,endy,groupy=1
		endx = xdim
		endy = ydim
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
	endif

	String model=strFromTagWave(tags,TAG_MODEL)
	if (strlen(model))
		wnote= ReplaceStringByKey("detectorModel",wnote,model,"=")
	endif

	String file_time=strFromTagWave(tags,TAG_DATETIME)
	if (strlen(file_time))
		//		reformat "2011:04:21 20:44:56" --> "2010-04-07T17:41:26-05:00"
		file_time[0,8] = ReplaceString(":",file_time[0,8],"-")
		file_time[10,10] = "T"
		wnote= ReplaceStringByKey("file_time",wnote,file_time,"=")
	endif

	//	Variable exposure = str2num(strFromTagWave(tags,TAG_EXPOSURE))
	//	if (exposure>0)
	//		wnote= ReplaceNumberByKey("exposure",wnote,exposure,"=")
	//	endif
	return wnote
End
//
Static Function/T strFromTagWave(tags,tagNum)
	Wave/T tags
	Variable tagNum				// desired tag number
	if (!WaveExists(tags))
		return ""
	elseif (WaveDims(tags)!=2 || DimSize(tags,1)!=5 || WaveType(tags,1)!=2)
		return ""
	endif

	String value=""				// tag value
	Variable i
	for (i=0;i<DimSize(tags,0);i+=1)
		if (tagNum==round(str2num(tags[i][0])))
			value = TrimTrailingWhiteSpace(tags[i][4])
			break					// only use first occurance
		endif
	endfor
	return value
End




//Function/T OnlyTiffFileName(full)
//	String full
//
//	String name=full
//	Variable ii
//
//	ii = -1
//	do
//		ii = strsearch(name, ":", 0)			// remove the path part
//		if (ii>=0)
//			name= name[ii+1,inf]
//		endif
//	while(ii>0)
//
//	ii = -1
//	do
//		ii = strsearch(name, "\\", 0)			// remove the path part
//		if (ii>=0)
//			name= name[ii+1,inf]
//		endif
//	while(ii>=0)
//
//	return name
//End
