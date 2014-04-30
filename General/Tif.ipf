#pragma rtGlobals=1		// Use modern global access method.
#pragma version = 1.02
#pragma ModuleName=TiffProc
#include "ImageDisplayScaling", version>=1.64
#include "Reorient"

//
// Routines for reading in and looking at16-bit TIF files. Eliot Specht, ORNL, based on WinView.ipf code by Jon Tischler, Oak Ridge National Lab
// TischlerJZ@ornl.gov
//
// It should be made more general, to handle other image files

// Nov  15, 2012, with version 1.01 added extras to the loaders, although is doesn't do anything here yet.

Menu "BackReflection",  dynamic
	MenuItemIfWaveClassExists("Display image plot with buttons...","tif*","DIMS:2"), Graph_imageMake($"",0)
	help={"Display an image and put buttons on the plot for listing spots"}
	MarqueeTiffMenuItem("Center a Gaussian and Add to pkList"), AddGaussianToList()
	MarqueeTiffMenuItem("Report Center of a Gaussian"), GaussianCenter(1)
	help={"fit a gaussian to an ROI an image, and print results.  The region is selected with a marquee in an image plot"}
	MarqueeTiffMenuItem("Report Center of Mass of a ROI"), GetCenterOfMass(NaN,NaN,NaN,NaN)
	help={"Compute the center of mass of a region on an image.  The region is selected with a marquee in an image plot"}
	"-"
	pkLIstTiffMenuItem("Fit List of Gaussians..."), FitListOfGaussians("","pkList",-1)
	GcoefTiffMenuItem("Reset pkList from Gaussians"), pkList = Reset_pkList_from_Gaussians(Gauss_coef)
	"-"
	"Load TIFF File...", LoadWinViewFile("")
	help={"Load a TIFF Image from the file, and then display the image"}
//	MenuItemIfWaveClassExists("Get TIFF Info","speImage*","DIMS:2"), WinViewInfo($"","")
//	help={"Get the other information stored with the winview image"}
//	"Show WinView Procedures", DisplayProcedure/B=WinViewProcedures "CenterOfMass"
	     MenuItemIfWaveClassExists("   Find Z range of image","*","DIMS:2"), getZrange($"",NaN)
	help={"for an image, find the z range for the specified % range"}
End

Menu "Load Waves"
	"Load Tiff File...", LoadTiffFile("")
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
Function/S GcoefTiffMenuItem(item)
	String item
	if (WaveExists(Gauss_coef))
		return item
	endif
	return "("+item
End
//
Function/S pkLIstTiffMenuItem(item)
	String item
	if (strlen(StrVarOrDefault("pkList","")))
		return item
	endif
	return "("+item
End

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
	if (numtype(TiffType))
		printf "ERROR, WinViewSaveImage(), cannot write TIFF files for WaveType(%s) = 0x%x\r",NameOfWave(image),WaveType(image)
		return ""
	endif

// wnote is not used for TIFF files
//	String wnote=note(image)
//	wnote = ReplaceNumberByKey("numType",wnote,igorType2WinView(WaveType(image)),"=")
//	wnote = ReplaceNumberByKey("xdim",wnote,DimSize(image,0),"=")
//	wnote = ReplaceNumberByKey("ydim",wnote,DimSize(image,1),"=")

	ImageSave /D=16/IGOR/T="tiff"/U image as fileName

//	Variable fid								// file id
//	Open/T="????"/A/Z=1 fid as fileName	// this acutally opens file
//	fileName = S_fileName
//	if (V_flag)
//		return ""
//	endif
//	Close fid
	return fileName
End

Function/S LoadTiffFile(fName,[extras])
	String fName											// fully qualified name of file to open
	String extras											// does nothing here yet
	extras = SelectString(ParamIsDefault(extras),extras,"")

	Variable refNum
	if (strlen((OnlyTiffFileName(fName)))<1)				// call dialog if no file name passed
		Open /D/M="image file"/R/T="????" refNum		// use /D to get full path name
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
			NewImageGraph(image,NaN)
//			Graph_imageMake(image,NaN)
		endif
	endif
	
//	print "Tif/TifLoadFile: Note(image) = ", note(image)
	
	return GetWavesDataFolder(image,2)
End

Function/T TiffReadROI(fileName,i0,i1,j0,j1,[extras])
	String fileName					// fully qualified name of file to open (will not prompt)
	Variable i0,i1,j0,j1				// pixel range of ROI (if i1 or j1<0 then use whole image)
	String extras						// does nothing here yet
	extras = SelectString(ParamIsDefault(extras),extras,"")

	Variable xdim
	Variable ydim
	Variable itype
	String wnote

	String inWaveName = TIffLoadROI(fileName,itype,xdim,i0,i1,j0,j1)	// name of ROI read in
	if (strlen(inWaveName)<1)
		return ""											// nothing read in
	endif

	String wName = OnlyTiffFileName(fileName)				// rename image based on the file name
	String separator = SelectString(stringmatch(IgorInfo(2),"Macintosh"),"\\",":")
//	String wName = ParseFilePath(3, fileName,separator,0,0)
	wName = CleanupName(wName,0)						// wave name based on the file name

	if (exists(wName))										// if wave already exists, create unique name
		wName = wName+"_"
		wName = UniqueName(wName, 1, 1)
	endif
	Rename $inWaveName $wName
	Wave image = $wName
	
	wnote = note(image)
	
	SetScale/P x 0,1,"pixel", image
	SetScale/P y 0,1,"pixel", image
	
//	String wnote=TiffReadHeader(fileName,image)			// wave note to add to file read in
	if (!strlen(wnote))
		return ""
	endif
	
	if (i0!=0 || j0!=0)										// sub-region of image does not start at (0,0)
		wnote = ReplaceNumberByKey("startxRead",wnote,i0,"=")
		wnote = ReplaceNumberByKey("startyRead",wnote,j0,"=")
	endif
	
	wnote = ReplaceStringByKey("waveClass",wnote,"speImage","=")
	wnote = ReplaceStringByKey("DetectorID",wnote,"Photonic Sciences","=")
	note /K image, wnote
	
//	print "Tif/TiffReadROI: note(image) = ", note(image)

	return GetWavesDataFolder(image,2)
End

Function/T TiffLoadROI(fileName,itype,xdim,i0,i1,j0,j1,[extras])
	String fileName					// fully qualified name of file to open (will not prompt)
	Variable itype					// WinView file type
									//  0	"float (4 byte)"
									//  1	"long integer (4 byte)"
									//  2	"integer (2 byte)"
									//  3	"unsigned integer (2 byte)"
									//  4	"string/char (1 byte)"
									//  5	"double (8 byte)"
									//  6	"signed int8 (1 byte)"
									//  7	"unsigned int8 (1 byte)"
	Variable xdim					// x size of whole array
	Variable i0,i1,j0,j1			// pixel range of ROI
									// for whole image use 0,xdim,0,ydim
	String extras					// does nothing here yet
	extras = SelectString(ParamIsDefault(extras),extras,"")

	variable ydim
	ImageLoad /T=tiff fileName
	
	if (V_flag != 1)
		return ""							// could not open file
	endif
	//SVAR S_waveNames=S_waveNames
	String wName = S_waveNames[0,strsearch(S_waveNames, ";", 0)-1]
	Wave wav = $wName

	String info = WaveInfo(wav,0)
	if (NumberByKey("NUMTYPE",info) != 16+64)	// unsigned 16-bit integer
		print "TIF file was not 16-bit"
		KillWaves /Z wav
		return ""
	endif

	String wnote=""											// wave note to add to file read in
	wnote = ReplaceStringByKey("imageFileName", wnote, ParseFilePath(0,fileName,":",1,0),"=")
	wnote = ReplaceStringByKey("imageFilePath", wnote, ParseFilePath(1,fileName,":",1,0),"=")
	String flatFile="", bkgFile=""

//	wnote = ReplaceStringByKey("num_Type",wnote,WinViewFileTypeString(header.datatype),"=")
//	wnote = ReplaceNumberByKey("numType",wnote,header.datatype,"=")

	xdim = DimSize(wav,0)
	ydim = DimSize(wav,1)
	
	wnote = ReplaceNumberByKey("xdim",wnote,xdim,"=")			// x size of image in file
	wnote = ReplaceNumberByKey("ydim",wnote,ydim,"=")			// y size of image in file
	
	
	i1 = (i1<1) ? xdim-1 : i1								// -1 flags use whole range
	j1 = (j1<1) ? ydim-1 : j1

	i0 = max(round(i0),0)
	i1 = max(round(i1),0)
	j0 = max(round(j0),0)
	j1 = max(round(j1),0)
	Variable nx = i1-i0+1
	Variable ny = j1-j0+1
	if (nx<1 || ny<1)							// nothing to read
		return ""
	endif
	Variable fType = 0x50 // Igor number type for 16-bit
	

	Duplicate /R=[i0,i1][j0,j1] wav, ROIwave
	Duplicate/O ROIwave, wav
	KillWaves ROIwave
	
	wnote = ReplaceNumberByKey("startx",wnote,i0,"=")
	wnote = ReplaceNumberByKey("endx",wnote,i1,"=")
	wnote = ReplaceNumberByKey("groupx",wnote,1,"=")
	wnote = ReplaceNumberByKey("starty",wnote,j0,"=")
	wnote = ReplaceNumberByKey("endy",wnote,j1,"=")	
	wnote = ReplaceNumberByKey("groupy",wnote,1,"=")	
	
	note /K wav, wnote
	
//	print "Tif/TifLoadROI: note(wav) = ", note(wav)
	
	return GetWavesDataFolder(wav,2)
End

Function/T OnlyTiffFileName(full)
	String full

	String name=full
	Variable ii

	ii = -1
	do
		ii = strsearch(name, ":", 0)			// remove the path part
		if (ii>=0)
			name= name[ii+1,inf]
		endif
	while(ii>0)

	ii = -1
	do
		ii = strsearch(name, "\\", 0)			// remove the path part
		if (ii>=0)
			name= name[ii+1,inf]
		endif
	while(ii>=0)

	return name
End
