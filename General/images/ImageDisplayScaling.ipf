#pragma rtGlobals=1		// Use modern global access method.
#pragma version = 2.12
#pragma ModuleName=ImageDisplayScaling
//
// Routines for rescaling the color table for images, by Jon Tischler, Oak Ridge National Lab
// TischlerJZ@ornl.gov
//
// to use these macros in your Igor experiment, put this file in the "User Procedures" folder
// in the "Igor Pro Folder", and uncomment the following line and put it in the top of the procedure
//folder
//
//	at verion 2.07, added the ImageROIstruct Structure and associated functions
//	at verion 2.12, bad pixel support, and also ROIofImage()


// =============================================================================================
// ==================================== Start of Menu Items ====================================

Menu "Graph"
//	"Set Aspect Ratio to Get Square Pixels \ Range",SetAspectToSquarePixels("")
	// MenuItemIfTopGraphImage("Set Aspect Ratio to Get Square Pixels"),SetAspectToSquarePixels("")
	MenuItemIfTopGraphImage("Sum Image Horizontally & Vertically"),sumImageToBothEdges($"",printIt=1)
	MenuItemIfWaveClassExists("Get Image Info...","speImage*;rawImage*","DIMS:2"), GenericWaveNoteInfo($"","",class="speImage*;rawImage*",options="DIMS:2",type="image")
End
Menu "Data"
	MenuItemIfWaveClassExists("Get Image Wave Info...","speImage*;rawImage*","DIMS:2"), GenericWaveNoteInfo($"","",class="speImage*;rawImage*",options="DIMS:2",type="image")
	"Load Bad Pixels...", LoadDefaultBadPixelImage($"","")
End

Menu "New"
	MenuItemIfWaveClassExists("Display image plot with buttons...","spe*;rawImage*","DIMS:2"), NewImageGraph($"")
End

Menu "Load Waves"
	"Load Jon's Generic Image File...",LoadGenericImageFile("")
End

Menu "Save Waves"
	"Save Jon's Generic Image File...",GenericSaveImage("",$"")
End


Menu "GraphMarquee", dynamic
	"-"
	MarqueeImageDisplayMenuItem("Add Gaussian to pkList"),/Q, AddGaussianToPeakList()
	help={"fit a gaussian to the Marquee, and add it to the peak list"}
	MarqueeImageDisplayMenuItem("Add COM to pkList"),/Q, Add_COM_ToPeakList()
	help={"find peak from Center Of Mass in the Marquee, and add it to the peak list"}
	MarqueeImageDisplayMenuItem("Remove Peak From pkList"),/Q, RemovePeakFromPeakList()		// RemovePeakFromList()
	help={"remove the peak in the marquee from the peak list"}
	MarqueeCsrImageDisplayMenuItem("Add Peak at Cursor A"),/Q, AddPeakAtCursorA()
	help={"add position of cursor A (the round one) to the peak list"}
	MarqueeImageDisplayMenuItem("Report Center of Gaussian"),/Q, Gaussian2dPeakCenter()
	help={"do a gaussian fit and report the result (does not add to peak list)"}
	MarqueeImageDisplayMenuItem("Report Center of Mass"),/Q, ShowCenterOfMass2d(1)
	help={"report the center of mass of the marquee (does not add to peak list)"}
	"-"
	MarqueeImageMenuItem("5%-95% scaling"),/Q, fivePercentScaling()
	help={"reset the z scaling of the image based on pixel range of this marquee, excluding the top and bottom 5%"}
	MarqueeImageMenuItem("full range scaling"),/Q, fullRangeScaling("")
	help={"reset the z scaling of the image based on pixel range of this marquee"}
	MarqueeImageMenuItem("Statistics"),/Q,statsOfROI()
	help={"write to the history some statistics about the pixels within this marquee"}
	"Make Marquee Match Plot Aspect",/Q,MatchMarqueeAspect()
	help={"resize this marquee so that its aspect ration matches that of the plot, used before an expand or shrink"}
	MarqueeImageMenuItem("Sum Image Horizontally & Vertically"), /Q,sumImageToBothEdges(ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";"))),printIt=1)
	help={"reset the z scaling of the image based on pixel range of this marquee, excluding the top and bottom 5%"}
End
//
Function/S MarqueeImageMenuItem(item)
	String item
	if (strlen(ImageNameList("",";"))<1)		// there are no images anywhere, so do not show
		return ""
		//	return "("+item
	endif
	return item
End
//
//Function/S MenuItemIfWaveClassExists(item,classes,options)	// moved to Utility_JZT.ipf
//	String item
//	String classes
//	String options
//	String list = WaveListClass(classes,"*",options)
////	String list = WaveListClass(classes,"*","DIMS:2")
//	return SelectString(strlen(list),"(","")+item
//End
//Menu "GraphMarquee"
//	"-"
//	"5%-95% scaling",/Q,fivePercentScaling()	
//	help={"reset the z scaling of the image based on pixel range of this marquee, excluding the top and bottom 5%"}
//	"full range scaling",/Q,fullRangeScaling()
//	help={"reset the z scaling of the image based on pixel range of this marquee"}
//	"Statistics",/Q,statsOfROI()
//	help={"write to the history some statistics about the pixels within this marquee"}
//	"Make Marquee Match Plot Aspect",/Q,MatchMarqueeAspect()
//	help={"resize this marquee so that its aspect ration matches that of the plot, used before an expand or shrink"}
//End
//
//Function/S MenuItemIfTopGraphImage(item)
//	String item
//	if (strlen(ImageNameList("",""))<1)
//		return "("+item					// top graph does not contain an image, so disable menu item
//	endif
//	return item
//End
//
//// for menus, enable menu item when a particular wave class is present on the graphName
//Function/S MenuItemsWaveClassOnGraph(item,classes,graphName)
//	String item					// item name to return, or "(item" to disable
//	String classes					// list of classes to check for (semi-colon separated)
//	String graphName				// name of graph, use "" for top graph
//
//	Variable i
//	String list = imageNameList("",";")	// first check the images
//	for (i=0;i<ItemsInList(list);i+=1)
//		Wave ww = ImageNameToWaveRef(graphName, StringFromList(i,list))
//		if (WaveInClass(ww,classes))
//			return item
//		endif
//	endfor
//
//	list = TraceNameList("",";",3)		// second check ordinary traces and contours
//	for (i=0;i<ItemsInList(list);i+=1)
//		Wave ww = TraceNameToWaveRef(graphName, StringFromList(i,list))
//		if (WaveInClass(ww,classes))
//			return item
//		endif
//	endfor
//	return "("+item						// top graph does not contain a wave of desired class
//End

Function/S MarqueeImageDisplayMenuItem(item)
	String item
	if (strlen(WinList("*","","WIN:1"))<1)
		return "("+item
	endif
	Variable V_flag
	GetMarquee/Z								// These menu items requre a Marquee present
	if (V_flag==0)
		return "("+item
	endif

	if (strlen(ImageNameList("",";"))<1)		// there are no images anywhere, so do not show
		return ""
	endif

	Wave image = ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";")))
	if (WaveInClass(image,"spe*;rawImage*;"))
		return item
	endif
	return "("+item							// it is an image, but the wrong kind of image
End

Function/S MarqueeCsrImageDisplayMenuItem(item)
	String item

	String wName = StringFromList(0,ImageNameList("",";"))
	if (!strlen(wName))
		return ""								// there are no images anywhere, so do not show
		//	return "("+item
	endif
	if (!WaveExists(ImageNameToWaveRef("",wName)))
		return ""								// this is not an image, do not show
		//	return "("+item
	endif
	if (strlen(CsrWave(A))<=0)				// These menu items requre that cursor A be on the image
		return "("+item						// We are on image, but no cursor
	endif

	Wave image = ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";")))
	if (WaveInClass(image,"spe*;rawImage*;"))	// wrong kind of image
		return item
	endif
		return ""
	return  "("+item
End

// ===================================== End of Menu Items =====================================
// =============================================================================================


// =============================================================================================
// ============================== Start of Getting Bad Pixel Info ==============================

Function/WAVE GetBadPixelsImage(wnote)
	// returns the bad pixels image, if it exists, wnote is from a loaded image
	String wnote

	String id = StringByKey("detectorID",wnote,"=")
	Variable Nx = round(NumberByKey("xDimDet",wnote,"="))
	Variable Ny = round(NumberByKey("yDimDet",wnote,"="))
	if (numtype(Nx+Nx) || Nx<1 || Ny<1 || strlen(id)<1)
		return $""
	endif

	String name =  "root:Packages:imageDisplay:"+ReplaceString("__",CleanupName(id+"_BadPixels",0),"_")
	Wave badPixels = $name
	if (Nx==DimSize(badPixels,0) && Ny==DimSize(badPixels,1) && WaveInClass(badPixels,"badPixelsImage*"))
		return badPixels
	endif
	return $""
End


Function/WAVE LoadDefaultBadPixelImage(image,fileName,[printIt])
	Wave image
	String fileName	// name of file with the bad pixels, may be ""
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt

	String wnote=""
	if (WaveExists(image))
		wnote=note(image)		// a generic image wave note with detectorID and size of image (xDimDet,yDimDet)
	else
		String infoSource, menuStr=""
		FUNCREF geoDetectorInfoWaveProto detectorInfoFunc=$"geoDetectorInfoWave"
		Wave/T dInfo = detectorInfoFunc()
		if (WaveExists(dInfo))
			menuStr = ReplaceString(":", StringByKey("menuStr",note(dInfo),"="),";")
			menuStr += SelectString(strlen(menuStr),""," ;")
		endif
		menuStr += reverseList(WaveListClass("speImage*;rawImage*","*","DIMS:2"))
		if (!strlen(menuStr))
			print "No geometrys or images found, nothing done."
			DoAlert 0, "No geometrys or images found, nothing done."
			return $""
		endif
		Prompt infoSource, "source of detector info", popup, menuStr
		DoPrompt "image",infoSource
		if (V_flag)
			return $""
		endif

		if (numpnts(dInfo))
			Variable i = WhichListItem(infoSource,menuStr)
			if (i>=0)
				wnote = dInfo[i]
			endif
		endif
		if (!strlen(wnote))
			Wave image = $infoSource
			if (WaveExists(image))
				wnote = note(image)
			endif
		endif
		printIt = 1
	endif
	if (printIt)
		printf "LoadDefaultBadPixelImage(%s,\"%s\")\r",NameOfWave(image),fileName
	endif

	String id = TrimBoth(StringByKey("detectorID",wnote,"="))
	Variable Nx = round(NumberByKey("xDimDet",wnote,"="))
	Variable Ny = round(NumberByKey("yDimDet",wnote,"="))
	if (numtype(Nx+Nx) || Nx<1 || Ny<1 || strlen(id)<1)
		return $""
	endif

	Wave pxyBad = ImageDisplayScaling#readBadPixelFile(fileName,id)
	Wave badImage = ImageDisplayScaling#badPixelList2image(Nx,Ny,pxyBad)
	if (!WaveExists(badImage))
		print "Error loading bad pixel file, nothing done."
		DoAlert 0, "Error loading bad pixel file, nothing done."
		return $""
	endif

	String name = ReplaceString("__",CleanupName(id+"_BadPixels",0),"_")

	if (printIt)
		fileName = StringByKey("badPixelFile",note(pxyBad),"=")
		printf "Loaded Bad Pixels from \"%s\", into image \"%s\"[%d][%d]\r", fileName,name, Nx,Ny
	endif

	name = "root:Packages:imageDisplay:"+name
	Duplicate/O badImage, $name
	return $name
End

Function/WAVE geoDetectorInfoWaveProto()
	return $""
End


Static Function/WAVE badPixelList2image(Nx,Ny,pxy)
	// take wave with bad pixels (probably from readBadPixelFile()) and return image with bad pixels
	// good pixels are 1, bad pixels are 0
	Variable Nx,Ny				// number of pixels in image, used to make badImage[][]
	Wave pxy						// list of bad pixels pxy[][2]

	if (!WaveExists(pxy) || numtype(Nx+Ny) || Nx<1 || Ny<1)
		return $""
	endif

	Variable N=DimSize(pxy,0)	// number of bad pixels in pxy[][2]
	Make/N=(Nx,Ny)/B/U/FREE badImage=1
	Variable  px,py, i
	for (i=0;i<N;i+=1)
		px = pxy[i][0]
		py = pxy[i][1]
		if (px<0 || px>=Nx || py<0 || py>=Ny)
			String str
			sprintf str, "Bad pixel = (%d, %d) is outside of detector size [%d, %d]",px,py,Nx,Ny
			print str
			DoAlert 0, str
			return $""
		endif
		badImage[px][py] = 0
	endfor

	String wnote=note(pxy)
	wnote = ReplaceStringByKey("waveClass",wnote,"badPixelsImage","=")
	Note/K badImage, wnote
	return badImage
End


Static Function/WAVE readBadPixelFile(fileName, id)
	// read a bad pixel file, and return free wave with the bad pixels
	String fileName
	String id

	PathInfo home
	String BadPixelFilters = "Bad Pixel Files (*.txt):.txt;All Files:.*;"
	Variable f
	if (V_flag)
		Open/F=BadPixelFilters/M="Bad Pixel File"/P=home/R/Z=2 f as fileName
	else
		Open/F=BadPixelFilters/M="Bad Pixel File"/R/Z=2 f as fileName
	endif
	fileName = S_fileName
	if (f==0 || V_flag)
		return $""
	endif

	FStatus f
	String buf = PadString("",V_logEOF,0x20)
	FBinRead f, buf
	Close f

	buf = ReplaceString("\t",buf," ")
	buf = ReplaceString("\r",buf,"\n")
	for ( ; strsearch(buf,"\n\n",0)>=0; )
		buf = ReplaceString("\n\n",buf,"\n")
	endfor
	buf = TrimBoth(buf)
	if (!StringMatch(StringFromList(0,buf,"\n"),"Bad Pixel File*"))
		return $""						// wrong file type
	endif

	String str = TrimBoth(StringFromList(1,buf,"\n"))
	if (strsearch(str,"$ID",0,2)!=0)
		return $""
	endif
	str = TrimBoth(str[3,Inf])
	id = TrimBoth(id)
	if (!StringMatch(str,id))
		return $""						// detector ID does not match
	endif

	buf = ReplaceString(",",buf," ")	// change to ";" separated list of pixel x y
	buf = ReplaceString("(",buf," ")
	buf = ReplaceString(")",buf," ")
	buf = ReplaceString("\n",buf,";")
	String pair
	Variable i, N=ItemsInlist(buf), px,py
	Variable Npix, Nalloc=100
	Make/N=(Nalloc,2)/I/FREE pxy
	for (i=2,Npix=0;i<N;i+=1)
		pair = StringFromList(i,buf)
		sscanf pair, "%d %d", px,py
		if (V_flag==2)
			if (Npix >= Nalloc)
				Nalloc += 10000
				Redimension/N=(Nalloc,-1) pxy
			endif
			pxy[Npix][0] = px
			pxy[Npix][1] = py
			Npix += 1
		endif
	endfor

	String wnote = "waveClass=badPixelsXY;"
	wnote = ReplaceStringByKey("detectorID",wnote,id,"=")
	wnote = ReplaceStringByKey("badPixelFile",wnote,S_fileName,"=")
	Note/K pxy, wnote

	if (Npix<1)
		WAVEClear pxy
	else
		Redimension/N=(Npix,-1) pxy
	endif
	return pxy
End

	//	Function test_LoadBadPixels()
	//		Make/N=(3,3)/FREE im
	//		Note im, "detectorID=PE1621 723-3335;xDimDet=2048;yDimDet=2048;"
	//		LoadDefaultBadPixelImage(im,"",printIt=1)
	//	End

	//	Function/WAVE test_GetBadPixelsImage()
	//		String wnote = "detectorID=PE1621 723-3335;xDimDet=2048;yDimDet=2048;"
	//		Wave badIm = GetBadPixelsImage(wnote)
	//		if (WaveExists(badIm))
	//			print GetWavesDataFolder(badIm,2)
	//		else
	//			print "could not find bad image"
	//		endif
	//	End

// =============================== End of Getting Bad Pixel Info ===============================
// =============================================================================================


// =============================================================================================
// ============================== Start of Generic Image Reading ===============================

Function imageLoadersAvailable()		// returns true if a real image loader is available
	Variable haveHDF5 = exists("LoadHDF5imageFile")==6	// flags indicating which type of images I know how to read
	Variable haveWinView = exists("LoadWinViewFile")==6
	Variable haveTiff = exists("LoadTiffFile")==6
	Variable havePilatusTiff = exists("LoadPilatusTiffImage")==6
	return (haveHDF5 || haveWinView || haveTiff || havePilatusTiff)
End



// This function understands HDF5, SPE, TIFF, and PilatusTIFF files, returns full path name to the image
Function/S LoadGenericImageFile(fileName,[extras])		// returns full path name to loaded image wave
	String fileName
	String extras
	extras = SelectString(ParamIsDefault(extras), extras,"")

	Variable haveHDF5 = exists("LoadHDF5imageFile")==6	// flags indicating which type of images I know how to read
	Variable haveWinView = exists("LoadWinViewFile")==6
	Variable haveTiff = exists("LoadTiffFile")==6
	Variable havePilatusTiff = exists("LoadPilatusTiffImage")==6

	if (strlen(ParseFilePath(3,fileName,":",0,0))<1)		// call dialog if no file name passed
		initImageDisplayScaling()
		setImageFileFilters()
		SVAR fileFilters=root:Packages:imageDisplay:ImageFileFilters
		Variable f
		Open /D/M="image file"/R/F=fileFilters f			// use /D to get full path name
		fileName = S_filename
	endif
	if (strlen(fileName)<1)								// no file name, quit
		return ""
	endif

	String extension = ParseFilePath(4,fileName,":",0,0), str=""
	if (haveHDF5 && (stringmatch(extension,"h5") || stringmatch(extension,"HDF")))
		FUNCREF LoadFileProtoShort func = $"LoadHDF5imageFile"
	elseif (haveWinView && stringmatch(extension,"SPE"))	// Load an SPE image
		FUNCREF LoadFileProtoShort func = $"LoadWinViewFile"
	elseif(haveTiff && (stringmatch(extension,"TIF") || stringmatch(extension,"TIFF")))
		FUNCREF LoadFileProtoShort func = $"LoadTiffFile"	// Load a TIFF image
	elseif(havePilatusTiff && (stringmatch(extension,"TIF") || stringmatch(extension,"TIFF")))
		FUNCREF LoadFileProtoShort func = $"LoadPilatusTiffImage"	// Load a Pilatus TIFF image
	else
		return ""
	endif

	if (strlen(extras))
		str = func(fileName,extras=extras)				// Load image with extras
	else
		str = func(fileName)								// Load an simple image
	endif		

	Wave image = $str
	if (WaveExists(image))									// make sure note of image contains waveClass of "rawImage"
		String wnote = AddClassToWaveNote(note(image),"rawImage")
		Note/K image, wnote
	endif

	String/G root:Packages:imageDisplay:imageExtension="."+extension
	return str													// return full path name to the image
End
//
Function/T LoadFileProtoShort(fName,[extras])
	String fName
	String extras
	return ""
End



Function/S ReadGenericROI(fileName,i0,i1,j0,j1,[extras])	// returns full path name to loaded image wave
	String fileName											// fully qualified name of file to open (will not prompt)
	Variable i0,i1,j0,j1									// pixel range of ROI (if i1 or j1<0 then use whole image)
	String extras
	extras = SelectString(ParamIsDefault(extras), extras,"")

	Variable f
	if (strlen((ParseFilePath(3,fileName,":",0,0)))<1)	// call dialog if no file name passed
		Open /D/M="image file"/R/T="????" f				// use /D to get full path name
		fileName = S_filename
	endif
	if (strlen(fileName)<1)								// no file name, quit
		return ""
	endif

	String extension = ParseFilePath(4,fileName,":",0,0), str=""
	if (exists("HDF5ReadROI")==6 && (stringmatch(extension,"h5") || stringmatch(extension,"HDF")))
		FUNCREF ReadGenericROIshortProto func = $"HDF5ReadROI"
	elseif (exists("WinViewReadROI")==6 && stringmatch(extension,"SPE"))
		FUNCREF ReadGenericROIshortProto func = $"HDF5ReadROI"
	elseif (exists("TiffReadROI")==6 && (stringmatch(extension,"TIF") || stringmatch(extension,"TIFF")))
		FUNCREF ReadGenericROIshortProto func = $"TiffReadROI"
	elseif (exists("PilatusTiffReadROI")==6 && (stringmatch(extension,"TIF") || stringmatch(extension,"TIFF")))
		FUNCREF ReadGenericROIshortProto func = $"PilatusTiffReadROI"
	else
		return ""
	endif

	if (strlen(extras))
		str = func(fileName,i0,i1,j0,j1,extras=extras)	// an image with extras
	else
		str = func(fileName,i0,i1,j0,j1)						// a simple image
	endif

	String/G root:Packages:imageDisplay:imageExtension="."+extension
	return str
End
//
Function/T ReadGenericROIshortProto(fileName,i0,i1,j0,j1,[extras])
	String fileName
	Variable i0,i1,j0,j1
	String extras
	return ""
End



Function/T ReadGenericHeader(fName,[extras])
	String fName					// fully qualified name of file to open (will not prompt)
	String extras					// optional switches
	extras = SelectString(ParamIsDefault(extras),extras,"")
	GetFileFolderInfo/Q /Z=1 fName							// check if file exists
	if (!V_isFile || V_Flag)									// file not there
		return ""
	endif

	String extension = ParseFilePath(4,fName,":",0,0), str="", funcName=""
	if (exists("HDF5ReadROI")==6 && (stringmatch(extension,"h5") || stringmatch(extension,"HDF")))
		funcName = "ReadHDF5header"
	elseif (exists("PilatusTiffReadHeader")==6 && (stringmatch(extension,"TIF") || stringmatch(extension,"TIFF")))
		funcName = "PilatusTiffReadHeader"
	elseif (exists("WinViewReadHeader")==6 && stringmatch(extension,"SPE"))
		funcName = "WinViewReadHeader"
	elseif (exists("TiffReadHeader")==6 && (stringmatch(extension,"TIF") || stringmatch(extension,"TIFF")))
		funcName = "TiffReadHeader"
	else
		return ""
	endif
	if (strlen(funcName))
		String/G root:Packages:imageDisplay:imageExtension="."+extension
		FUNCREF ReadGenericHeader func = $funcName
		if (strlen(extras))
			str = func(fName,extras=extras)
		else
			str = func(fName)
		endif
	endif
	return str
End


Static Function/S setImageDefaultFileExtension(hard)// note, the extension contains a leading "."
	Variable hard									// hard=1 forces a reset even if image file extension is already set, if hard==0, only does something if no extension
	String ext=StrVarOrDefault("root:Packages:imageDisplay:imageExtension","")
	if (strlen(ext)>0 && !hard)
		return ext
	endif

	String imageTypes="", imageType
	imageTypes += SelectString(Exists("LoadHDF5imageFile")==6,"","HDF5;")
	imageTypes += SelectString(Exists("LoadWinViewFile")==6,"","SPE (WinView);")
	imageTypes += SelectString(Exists("LoadTiffFile")==6,"","TIFF (only for 4500S);")
	imageTypes += SelectString(Exists("LoadPilatusTiffImage")==6,"","Pilatus TIFF;")
	Prompt imageType,"Default image Type",popup,imageTypes
	DoPrompt "Default Image Type",imageType
	if (V_flag)
		return ""
	endif

	strswitch(UpperStr(imageType[0,2]))
		case "HDF":
			ext = ".h5"
			break
		case "SPE":
			ext = ".SPE"
			break
		case "TIF":
			ext = ".tif"
			break
		case "Pil":
			ext = ".tif"
			break
	endswitch
	//	SVAR ImageFileFilters=root:Packages:imageDisplay:ImageFileFilters
	//	Variable f
	//	Open/R/D/F=ImageFileFilters/M="Pick an image file" f
	//	ext = ParseFilePath(4,S_fileName,":",0,0)
	//	ext = SelectString(strlen(ext),"",".") + ext		// add the leading "."
	NewDataFolder/O root:Packages
	NewDataFolder/O root:Packages:imageDisplay
	String/G root:Packages:imageDisplay:imageExtension=ext
	return ext
End

Static Function/S setImageFileFilters()
	NewDataFolder/O root:Packages
	NewDataFolder/O root:Packages:imageDisplay
	String/G root:Packages:imageDisplay:ImageFileFilters=""
	SVAR filters=root:Packages:imageDisplay:ImageFileFilters
	filters += SelectString(Exists("LoadHDF5imageFile")==6,"","HDF Files (*.h5, *.hdf5, *.hdf):.h5,.hdf5,.hdf,;")
	filters += SelectString(Exists("LoadWinViewFile")==6,"","SPE Files (*.SPE):.SPE;")
	Variable aTiff = Exists("LoadTiffFile")==6 || Exists("LoadPilatusTiffImage")==6
	filters += SelectString(aTiff,"","TIFF Files (*.tif,*.tiff):.tif,.tiff;")
	filters += "All Files:.*;"
	return filters
End

// ================================ End of Generic Image Reading ===============================
// =============================================================================================



// =============================================================================================
// =============================== Start of Generic Image Saving ===============================

Function/T GenericSaveImage(fName,image,[extras])
	String fName					// fully qualified name of file to open (will not prompt)
	Wave image
	String extras
	extras = SelectString(ParamIsDefault(extras),extras,"")

	Variable printIt=NumberByKey("printIt",extras)
	printIt = numtype(printIt) ? strlen(GetRTStackInfo(2))<1 : !(!printIt)
	if (!WaveExists(image))
		String classes = StringByKey("classes",extras)
		classes = SelectString(ItemsInList(classes),"speImage*;rawImage*;darkImage;imageMask",classes)
		String imageName="", imageList = WaveListClass(classes,"*","DIMS:2")
		if (ItemsInList(imageList)==1)
			imageName = StringFromLIst(0,imageList)
		elseif (ItemsInList(imageList)>1)
			Prompt imageName,"image to save",popup,imageList
			DoPrompt "image",imageName
			if (V_flag)
				return ""
			endif
		endif
		Wave image = $imageName
		printIt = 1
	endif
	if (!WaveExists(image))
		print "ERROR -- GenericSaveImage(), the image wave does not exist"
		return ""
	elseif (WaveDims(image)!=2)
		printf "ERROR -- GenericSaveImage(), only writes 2-d waves, '%s' is %d-d\r",WaveDims(image),NameOfWave(image)
		return ""
	endif

	Variable haveHDF5 = exists("LoadHDF5imageFile")==6	// flags indicating which type of images I know how to read
	Variable haveWinView = exists("LoadWinViewFile")==6
	Variable haveTiff = exists("LoadTiffFile")==6
	Variable havePilatusTiff = exists("LoadPilatusTiffImage")==6

	String fileFilters=""
	fileFilters += SelectString(haveHDF5,"","HDF Files (*.h5):.h5,;")
	fileFilters += SelectString(haveTiff || havePilatusTiff,"","Tiff Files (*.tif):.tif,;")
	fileFilters += SelectString(haveWinView,"","SPE Files (*.spe):.spe,;")

	Variable f
	if (strlen(ParseFilePath(3,fName,":",0,0))<1)				// call dialog if no file name passed
		Open /D/M="Output image file"/F=fileFilters f			// use /D to get full path name
		fName = S_filename
	endif
	fName = TrimTrailingWhiteSpace(fName)					// remove trailing spaces
	if (strlen(fName)<1)										// no file name, quit
		return ""
	endif
	if (strsearch(fName, ".",0)<0)
		fName += StrVarOrDefault("root:Packages:imageDisplay:imageExtension",".h5")
	endif
	if (printIt)
		printf "%sGenericSaveImage(\"%s\",%s\r",SelectString(strlen(GetRTStackInfo(2)),"","¥"),fName,NameOfWave(image)
		if (strlen(extras))
			printf ", extras=\"%s\"",extras
		endif
		printf ")\r"
	endif

	String funcName="", extension=ParseFilePath(4,fName,":",0,0)
	if (haveHDF5 && stringmatch(extension,"h5"))
		funcName = "HDFSaveImage"
	elseif ((haveTiff || havePilatusTiff) && stringmatch(extension,"tif"))
		funcName = "TiffSaveImage"
	elseif (haveWinView && stringmatch(extension,"spe"))
		funcName = "WinViewSaveImage"
	endif

	if (exists(funcName))
		FUNCREF protoSaveImage funcSaveImage=$"HDFSaveImage"
		fName = funcSaveImage(fName,image)
	endif
	return fName
End
//
Function/T protoSaveImage(fName,image)
	String fName							// fully qualified name of file to open (will not prompt)
	Wave image
	return ""
End

// ================================ End of Generic Image Saving ================================
// =============================================================================================



// =============================================================================================
// ================================ Start of Graph with Buttons ================================

Function/S NewImageGraph(image,[withButtons,kill])
	Wave image
	Variable withButtons				// ignored if ButtonBoxesProc() does ot exist
	Variable kill							// optionally set the /K flag in Display command
	withButtons = ParamIsDefault(withButtons) ? 0 : withButtons
	kill = ParamIsDefault(kill) || numtype(kill) ? 0 : limit(round(kill),0,3)
	Variable NOButtons = (exists("ButtonBoxesProc")!=6)

	if (!WaveExists(image) || numtype(withButtons))
		String wList = reverseList(WaveListClass("speImage*;rawImage*","*","DIMS:2"))
		String wName=SelectString(WaveExists(image),StringFromList(0,wList),NameOfWave(image))
		if (strlen(wName)<1)
			return ""
		endif
		withButtons = numtype(withButtons) ? 0 : withButtons
		withButtons = withButtons ? 1 : 2
		Prompt withButtons,"put the buttons on the graph",popup,"Yes;No"
		Prompt wName, "image to show", popup, wList
		if (NOButtons)
			if (ItemsInList(wList)>1)
				DoPrompt "choose image", wName
			endif
		else
			DoPrompt "choose image", wName, withButtons
			withButtons = (withButtons==1) ? 1 : 0
		endif
		if(V_flag)
			return ""
		endif
		Wave image = $wName
	endif
	if (!WaveExists(image))
		DoAlert 0, "image does not exist"
		return ""
	endif

	String win=StringFromList(0,WindowsWithWave(image,1))
	if (strlen(win))
		DoWindow/F $win
		return ""
	else
		Display/W=(345,44,822,440)/K=(kill)
		AppendImage image
		if (withButtons && !NOButtons)
			Button boxes,pos={0,25},size={60,21},proc=ButtonBoxesProc,title="±boxes"
			FUNCREF ImageButtonRGBproto  func = $"ImageButtonRGB"
			Wave rgb = func(image)
			if (WaveExists(rgb))
				Button boxes fColor=(rgb[0],rgb[1],rgb[2])		// allow other routines to color the button
			endif
			SetWindow kwTopWin,userdata(boxes)=  "boxesOn=0;"
		endif
		if (exists("getFittedPeakInfoHook")==6)
			SetWindow kwTopWin hook(peakInfo)=getFittedPeakInfoHook
		endif
		Execute "GraphImageStyle()"
		SetAspectToSquarePixels("")
	endif
	return GetWavesDataFolder(image,2)
End
//
Function/WAVE ImageButtonRGBproto(image)	// used by NewImageGraph() in ImageDisplayScaling.ipf
	Wave image											// default is none
	return $""
End


//	TextBox/C/N=textTitle/F=0/B=1/X=1/Y=1 imageTitleText(PMN_10_Wire1_1_20)
Static Function/T imageTitleText(image)	// test suitable for a title text box
	Wave image
	if (!WaveExists(image))
		return ""
	endif
	String wnote = note(image)
	if (strlen(wnote)<1)
		return ""
	endif
	String title=StringByKey("title",wnote,"="), user=StringByKey("userName",wnote,"="), sample=StringByKey("sampleName",wnote,"=")
	String path=StringByKey("imageFilePath",wnote,"=")+StringByKey("imageFileName",wnote,"=")
	String file_time=StringByKey("file_time",wnote,"=")
	Variable scanNum = NumberByKey("scanNum",wnote,"=")
	Variable depth = NumberByKey("depth",wnote,"=")
	Variable keV = stringmatch(StringByKey("MonoMode",wnote,"="),"monochromatic") ? NumberByKey("keV",wnote,"=") : NaN
	Variable exposure = NumberByKey("exposure",wnote,"="), add=0
	Variable X1 = NumberByKey("X1",wnote,"="), Y1=NumberByKey("Y1",wnote,"="), Z1=NumberByKey("Z1",wnote,"=")
	Variable H1=NumberByKey("H1",wnote,"="), F1=NumberByKey("F1",wnote,"=")

	String str="", titleStr=""
	str += SelectString(strlen(title),"",title)
	if (strlen(sample))
		str += SelectString(strlen(str),"",", ")
		str += sample
	endif
	titleStr = str
	Variable shrink=0

	if (strlen(file_time))
		str = "\\Zr077"+niceify_file_time(file_time)
		shrink = 1
		titleStr += SelectString(strlen(titleStr),"","\r")+str
	endif

	if (strlen(user) || scanNum>0)
		str = SelectString(shrink,"\\Zr077","") + user
		shrink = 1
		if (scanNum>0)
			str += SelectString(strlen(user),"",", ")
			str += "scan# = "+num2str(scanNum)
		endif
		titleStr += SelectString(strlen(titleStr),"","\r")+str
	endif

	if (exists("imageTitleTextMore")==6)
		FUNCREF imageTitleTextMoreProto  func=$"imageTitleTextMore"
		str = func(image)
		if (strlen(str))
			str = SelectString(shrink,"\\Zr077","") + str
			shrink = 1
			titleStr += SelectString(strlen(titleStr),"","\r")+str
		endif
	endif

	add = 0
	str = ""
	if (numtype(depth)==0)
		str += SelectString(add,"",", ")
		add = 1
		str += "depth="+num2str(depth)+"µm"
	endif
	if (numtype(keV)==0)
		str += SelectString(add,"",", ")
		add = 1
		str += "keV="+num2str(keV)+"keV"
	endif
	if (numtype(exposure)==0)
		str += SelectString(add,"",", ")
		add = 1
		str += "exposure="+num2str(exposure)+"s"
	endif
	if (strlen(str))
		str = SelectString(shrink,"\\Zr077","") + str
		shrink = 1
		titleStr += SelectString(strlen(titleStr),"","\r")+str
	endif

	add = 0
	str = ""
	if (numtype(X1)==0)
		str += SelectString(add,"",", ")
		add = 1
		str += "X="+num2str(X1)
	endif
	if (numtype(H1+F1)==0)
		str += SelectString(add,"",", ")
		add = 1
		str += "H="+num2str(H1)
		str += SelectString(add,"",", ")
		add = 1
		str += "F="+num2str(F1)
	else
		if (numtype(Y1)==0)
			str += SelectString(add,"",", ")
			add = 1
			str += "Y="+num2str(Y1)
		endif
		if (numtype(Z1)==0)
			str += SelectString(add,"",", ")
			add = 1
			str += "Z="+num2str(Z1)
		endif
	endif
	if (strlen(str))
		str = SelectString(shrink,"\\Zr077","") + str
		shrink = 1
		titleStr += SelectString(strlen(titleStr),"","\r")+str
	endif
	str = SelectString(strlen(path),"","\\Zr060"+path)
	if (strlen(str))
		titleStr += SelectString(strlen(titleStr),"","\r")+str
	endif

	titleStr = TrimFrontBackWhiteSpace(titleStr)	// remove trailing <CR>
	return titleStr
End
//
Static Function/T niceify_file_time(timeStr)
	String timeStr

	Variable i=strsearch(timeStr,"T",0)
	if (i>0)
		String tz=""					// fix up the time zone
		Variable iz=strsearch(timeStr,"+",i)
		if (iz<0)
			iz = strsearch(timeStr,"-",i)
		endif
		if (iz>0)						// found time zone
			tz = timeStr[iz,Inf]		// time zone
			tz = ReplaceString(":00",tz,"")
			tz = ReplaceString("+0",tz,"+")
			tz = ReplaceString("-0",tz,"-")
			timeStr[iz,Inf] = " ("+tz+")"
		endif
	endif
	timeStr = ReplaceString("T",timeStr,"   ")
	return timeStr
End
//
Function/T imageTitleTextMoreProto(image)
	Wave image
	return ""
End
//
//		Function/T imageTitleTextMore(image)
//			Wave image
//			String str,title="\[9VO\B2\M\]9(", wnote=note(image)
//			Variable R = NumberByKey("VO2_Reistance",wnote,"=")
//			Variable T = NumberByKey("VO2_Temp",wnote,"=")
//			Variable ep = NumberByKey("VO2_Epoch",wnote,"=")
//			Variable add=0
//			if (numtype(R) && numtype(T))
//				return ""
//			endif
//			if (numtype(T)==0)
//				sprintf str,"%.2f¡C",T
//				title += SelectString(add,"",",  ")+str
//				add = 1
//			endif
//			if (numtype(R)==0)
//				sprintf str,"%.3f\F'Symbol'W\F]0",R
//				title += SelectString(add,"",",  ")+str
//				add = 1
//			endif
//			if (numtype(ep)==0)
//				sprintf str,"Æt=%.1f",ep
//				title += SelectString(add,"",",  ")+str
//				add = 1
//			endif
//			return title+")\r"
//		End

// ================================= End of Graph with Buttons =================================
// =============================================================================================



// =============================================================================================
// =================================== Start of roi Structure ==================================

Structure ImageROIstruct	// a ROI on an image
	int16		empty				// 0 is not empty
	int32		xLo				// must be >= 0
	int32		xHi				// must be >= xLo
	int32		yLo				// must be >= 0
	int32		yHi				// must be >= yLo
	int32		Nx					// size of this roi
	int32		Ny					// computed from yHi-yLo+1
EndStructure
//
Function ImageROIstructInit(roi)
	STRUCT ImageROIstruct &roi
	roi.empty = 0
	roi.xLo = 0	;	roi.xHi = 0
	roi.yLo = 0	;	roi.yHi = 0
	roi.Nx = 1	;	roi.Ny = 1
End
//
Function/T ImageROIstruct2str(roi)
	STRUCT ImageROIstruct &roi
	String out

	if (roi.empty)
		out = "ROI is empty"
	elseif (!ImageROIstructValid(roi))
		out = "INVALID roi structure"
	else
		sprintf out, "ix=[%g, %g] (N=%g),  iy=[%g, %g] (N=%g)\r", roi.xLo,roi.xHi,roi.Nx, roi.yLo,roi.yHi,roi.Ny
	endif
	return out
End
//
Function ImageROIstructValid(roi)
	// check for obvious errors and update Nx,Ny, return 1 if empty, 2 if really invalid
	STRUCT ImageROIstruct &roi

	if (roi.empty)					// if empty, the values do not matter
		return 1
	endif

	Variable err=0
	if (roi.xLo <= roi.xHi && roi.xLo >= 0)
		roi.Nx = roi.xHi - roi.xLo + 1
	else
		roi.Nx = 0
		err = 1
	endif

	if (roi.yLo <= roi.yHi && roi.yLo >= 0)
		roi.Ny = roi.yHi - roi.yLo + 1
	else
		roi.Ny = 0
		err = 1
	endif

	return !err
End

// returns TRUE if roi does NOT fit in the image
Function roiNOTinImage(image, roi)
	Wave image				// a 2D array
	STRUCT ImageROIstruct &roi
	if (!ImageROIstructValid(roi))
		return 1				// roi not even valid
	elseif (roi.xHi >= DimSize(image,0))
		return 1
	elseif (roi.yHi >= DimSize(image,1))
		return 1
	endif	
	return 0
End

//Function aaa()
//	STRUCT ImageROIstruct roi
//	initImageROIstruct(roi)
//	print ImageROIstructValid(roi)
//	print ImageROIstruct2str(roi)
//
//	roi.empty = 1
//	print ImageROIstructValid(roi)
//	print ImageROIstruct2str(roi)
//
//	roi.xLo = -3
//	roi.empty = 0
//	print ImageROIstructValid(roi)
//	print ImageROIstruct2str(roi)
//End

// ==================================== End of roi Structure ===================================
// =============================================================================================



// =============================================================================================
// =================================== Start of Generic Util ===================================

Function/WAVE ROIofImage(image, startx,groupx,endx, starty,groupy,endy)
	// This never returns image, but alwasy a FREE copy or subset of image
	Wave image
	Variable startx,groupx,endx
	Variable starty,groupy,endy

	if (!WaveExists(image))
		return $""
	elseif (WaveDims(image)!=2)				// only works on image (i.e. 2D arrays)
		return $""
	endif

	Duplicate/FREE image, ROI
	Variable Nx=DimSize(image,0), Ny=DimSize(image,1)
	if (startx==0 && starty==0 && groupx==1 && groupy==1 && endx==(Nx-1) && endy==(Ny-1))
		return ROI									// ROI is same size as image
	endif

	// need to extract an ROI
	Nx = floor( (endx-startx+1)/groupx )	// size of requested ROI
	Ny = floor( (endy-starty+1)/groupy )
	Redimension/N=(Nx,Ny) ROI
	ROI = image[p*groupx + startx][q*groupy + starty]	// does not do grouping, just pick first pixel
	return ROI
End


Proc Graph_ImageStyle() : GraphStyle
	GraphImageStyle()
EndMacro
//
Function GraphImageStyle()
	DoUpdate
	String wName = StringFromList(0,ImageNameList("",";"))
	Wave image = ImageNameToWaveRef("",wName)
	if (!WaveExists(image))
		Abort "Unable to find image on top graph"
	endif
	ModifyGraph/Z width={Aspect,DimSize(image,0)/DimSize(image,1)}
	ModifyGraph/Z grid=1, tick=2, mirror=1, minor=1, gridStyle=1
	ModifyGraph/Z lowTrip=0.001, standoff=0, axOffset(bottom)=-1
	ModifyGraph grid=1,gridRGB=(45000,45000,65535)
	ControlInfo  button0
	if (stringmatch(IgorInfo(2),"Macintosh" ))
		ModifyGraph axOffset(left)=(V_flag==1) ? 0.9 : -1.1 		// mac
	else
		ModifyGraph axOffset(left)=(V_flag==1) ? -1.7 : -2.1 		// pc
	endif
	Variable/C lohi = getZrange(image,5)
	Variable lo = real(lohi), hi = imag(lohi)
	if (!(WaveType(image) & 0x40) && lo<0)		// if all values positive, use Terrain256
		lo = abs(lo)
		hi = abs(hi)
		hi = max(lo,hi)
		lo = -hi
		ModifyImage $wName ctab= {lo,hi,RedWhiteBlue256,0}
	else
		ModifyImage $wName ctab= {lo,hi,Terrain256,1}
	endif
	Label/Z left SelectString(strlen(WaveUnits(image,0)),"y","y  (\\U)")
	Label/Z bottom SelectString(strlen(WaveUnits(image,0)),"x","x  (\\U)")
	SetAxis/A/R left
	ShowInfo

	String str = imageTitleText(image)
	if (strlen(str))
		TextBox/C/N=textTitle/F=0/B=1/X=1/Y=1 str
	endif
End

Function/C getZrange(image,delta,[plane,lo,hi])
	Wave image
	Variable delta				//=5,  clip amount (in %), not fraction
	Variable plane
	Variable lo,hi				// these will override delta for the lo and hi ends separately
	plane = ParamIsDefault(plane) ? NaN : plane
	plane = plane>=0 ? round(plane) : NaN
	Variable dLo=lo, dHi=hi
	dLo = ParamIsDefault(lo) ? NaN : dLo
	dLo = numtype(dLo) ? NaN : limit(dLo,0,99.9)
	dHi = ParamIsDefault(hi) ? NaN : dHi
	dHi = numtype(dHi) ? NaN : limit(dHi,0,99.9)

	Variable printIt=0
	if (!WaveExists(image) || numtype(delta) || delta<0)
		delta = (numtype(delta) || delta<0) ? 5 : delta
		String sImage=StringFromList(0, WaveList("*",";","DIMS:2,WIN:" ))
		Prompt sImage, "source image", popup, WaveList("*",";","DIMS:2" )
		Prompt delta, "percent to trim from each end"
		DoPrompt "pick image", sImage,delta
		if (V_flag)
			return NaN
		endif
		Wave image = $sImage
		printIt=1
	endif
	if (!WaveExists(image))
		Abort "input wave does not exist"
	endif
	if (numtype(delta) || delta<0 || delta > 99.9)
		Abort "delta ="+num2str(delta)+", which is too small"
	endif
	delta /= 100
	dLo = numtype(dLo) ? delta : dLo/100			// lo is the lower fraction cut-off
	dHi = numtype(dHi) ? delta : dHi/100			// hi is the upper fraction cut-off
	plane = DimSize(image,2)>1 ? plane : NaN	// only use plane for 3D data

	Variable N
	if (DimSize(image,0)*DimSize(image,1) < 2e5)	// not too big
		Duplicate/FREE image, values
		Variable Nx=DimSize(image,0), Ny=DimSize(image,1)
		N = Nx*Ny
		Redimension/N=(Nx,Ny) values
		if (plane>0)
			values = image[p][q][plane]
		endif
		Redimension/N=(N) values
		Sort values, values
		WaveStats/M=1/Q values
		N = V_npnts
		Redimension/N=(N) values		// trim off any NaN's

		lo = values[round(dLo*(N-1))]
		hi = values[round((1-dHi)*(N-1))]
		if (printIt)
			printf "the [%g%%, %g%%] range = [%g, %g]\r",dLo*100,100*(1-dHi),lo,hi
		endif
		return cmplx(lo,hi)


	else							// fancy way for bigger images
		if (plane>0)
			ImageStats/M=1/P=(plane) image	// check if no range in image (avoids Histogram error)
		else
			ImageStats/M=1 image				// check if no range in image (avoids Histogram error)
		endif
		if (V_max==V_min)
			return cmplx(V_min,V_max)
		endif

		if (plane>0)
			N = DimSize(image,0)*DimSize(image,1)
			Duplicate/FREE image, image0
			Redimension/N=(DimSize(image,0)*DimSize(image,1)) image0
			image0 = image[p][q][plane]
		else
			N = numpnts(image)/2
			Wave image0=image
		endif
		Make/N=(N)/O hh_temp_
		SetScale/P x 0,1,"", hh_temp_
		Histogram/B=1 image0,hh_temp_
		Wave hist=hh_temp_
		Integrate hist

		Variable hmin = hist[0]
		Variable dh = hist[N-1]-hmin
		lo = BinarySearchInterp(hist, dLo*dh+hmin)
		hi = BinarySearchInterp(hist, (1-dHi)*dh+hmin)
		lo = DimOffset(hist,0) + lo *DimDelta(hist,0)
		hi = DimOffset(hist,0) + hi *DimDelta(hist,0)
		if (printIt)
			printf "the [%g%%, %g%%] range = [%g, %g]\r",dLo*100,100*(1-dHi),lo,hi
		endif
		KillWaves/Z hist
		return cmplx(lo,hi)
	endif
End
//Function/C getZrange(image,delta,[plane,lo,hi])
//	Wave image
//	Variable delta				//=5,  clip amount (in %), not fraction
//	Variable plane
//	Variable lo,hi				// these will override delta for the lo and hi ends separately
//	plane = ParamIsDefault(plane) ? NaN : plane
//	plane = plane>=0 ? round(plane) : NaN
//	Variable dLo=lo, dHi=hi
//	dLo = ParamIsDefault(lo) ? NaN : dLo
//	dLo = numtype(dLo) ? NaN : limit(dLo,0,99.9)
//	dHi = ParamIsDefault(hi) ? NaN : dHi
//	dHi = numtype(dHi) ? NaN : limit(dHi,0,99.9)
//
//	Variable printIt=0
//	if (!WaveExists(image) || numtype(delta) || delta<0)
//		delta = (numtype(delta) || delta<0) ? 5 : delta
//		String sImage=StringFromList(0, WaveList("*",";","DIMS:2,WIN:" ))
//		Prompt sImage, "source image", popup, WaveList("*",";","DIMS:2" )
//		Prompt delta, "percent to trim from each end"
//		DoPrompt "pick image", sImage,delta
//		if (V_flag)
//			return NaN
//		endif
//		Wave image = $sImage
//		printIt=1
//	endif
//	if (!WaveExists(image))
//		Abort "input wave does not exist"
//	endif
//	if (numtype(delta) || delta<0 || delta > 99.9)
//		Abort "delta ="+num2str(delta)+", which is too small"
//	endif
//	delta /= 100
//	dLo = numtype(dLo) ? delta : dLo/100
//	dHi = numtype(dHi) ? delta : dHi/100
//	plane = DimSize(image,2)>1 ? plane : NaN	// only use plane for 3D data
//
//	if (plane>0)
//		ImageStats/M=1/P=(plane) image	// check if no range in image (avoids Histogram error)
//	else
//		ImageStats/M=1 image				// check if no range in image (avoids Histogram error)
//	endif
//	if (V_max==V_min)
//		return cmplx(V_min,V_max)
//	endif
//
//	Variable N
//	if (plane>0)
//		N = DimSize(image,0)*DimSize(image,1)
//		Duplicate/FREE image, image0
//		Redimension/N=(DimSize(image,0)*DimSize(image,1)) image0
//		image0 = image[p][q][plane]
//	else
//		N = numpnts(image)/2
//		Wave image0=image
//	endif
//	Make/N=(N)/O hh_temp_
//	SetScale/P x 0,1,"", hh_temp_
//	Histogram/B=1 image0,hh_temp_
//	Wave hist=hh_temp_
//	Integrate hist
//
//	Variable hmin = hist[0]
//	Variable dh = hist[N-1]-hmin
//	lo = BinarySearchInterp(hist, dLo*dh+hmin)
//	hi = BinarySearchInterp(hist, (1-dHi)*dh+hmin)
//	lo = DimOffset(hist,0) + lo *DimDelta(hist,0)
//	hi = DimOffset(hist,0) + hi *DimDelta(hist,0)
//	if (printIt)
//		printf "the [%g%%, %g%%] range = [%g, %g]\r",dLo*100,100*(1-dHi),lo,hi
//	endif
//	KillWaves/Z hist
//	return cmplx(lo,hi)
//End


Structure box
	Variable xlo,xhi, ylo,yhi		// edges of a box
EndStructure

Function fivePercentScaling()
	// change range of current color table to span [5%-95%] of the values in the marquee
	String imageName = StringFromList(0,ImageNameList("",";"))
	Wave image = ImageNameToWaveRef("",imageName)

	STRUCT box b
	if (pointRangeOfMarquee("",imageName,b))
		return 1
	endif

	String wName = UniqueName("roi",1,0)
	Make/N=(DimSize(image,0),DimSize(image,1))/B/U $wName
	Wave roi = $wName
	roi = 1
	roi[b.xlo,b.xhi][b.ylo,b.yhi] = 0

	ImageHistogram/R=roi/I image
	Wave W_ImageHist=W_ImageHist
	Integrate W_ImageHist
	Variable nn = W_ImageHist[Inf]
	W_ImageHist /= nn
	Variable lo, hi
	lo = BinarySearchInterp(W_ImageHist,0.05)
	hi = BinarySearchInterp(W_ImageHist,0.95)
	lo = (numtype(lo)==2) ?0 : lo
	hi = (numtype(lo)==2) ? numpnts(hist)-1 : hi
	lo = pnt2x(W_ImageHist,lo)
	hi = pnt2x(W_ImageHist,hi)

	ModifyOnly_ctab_range("",imageName,lo,hi)
	KillWaves/Z W_ImageHist, roi
End

// change range of current color table to exactly span the range of values in the marquee
Function fullRangeScaling(gName)
	String gName				// graph name, use "" for the top graph
	String imageName = StringFromList(0,ImageNameList(gName,";"))
	Wave image = ImageNameToWaveRef(gName,imageName)

	Variable lo,hi
	String color = get_ctab_range(gName,imageName,lo,hi)
	Variable sym = WhichListItem(color,"RedWhiteBlue;RedWhiteBlue256;BlueRedGreen;RedWhiteGreen")>=0 && -lo==hi	// a symmtric scaling, preserve symmetry

	STRUCT box b
	if (pointRangeOfMarquee("",imageName,b)==0)
		ImageStats/M=1/G={b.xlo,b.xhi,b.ylo,b.yhi} image
	else
		ImageStats/M=1 image
	endif
	lo = V_min  ;  hi = V_max
	if (sym)
		hi = max(abs(V_min),abs(V_max))
		lo = -hi
	endif
	ModifyOnly_ctab_range(gName,imageName,lo,hi)
End

Function/T statsOfROI()		// print out the statics of a Marquee
	String imageName = StringFromList(0,ImageNameList("",";"))
	Wave image = ImageNameToWaveRef("",imageName)
	STRUCT box b
	if (pointRangeOfMarquee("",imageName,b))
		return ""
	endif
	if (pointRangeOfMarquee("",imageName,b)==0)
		ImageStats/G={b.xlo,b.xhi,b.ylo,b.yhi} image
	else
		ImageStats image
	endif
	Variable Vsum=V_avg*V_npnts							// ImageStats does not give the sum
	Variable saturated=0, type = NumberByKey("NUMTYPE",WaveInfo(image,0))
	saturated = (type==80) ? 65535 : saturated		// unsigned 16 bit int
	saturated = (type==16) ? 32767 : saturated		// signed 16 bit int
	saturated = (type==74) ? 255 : saturated		// unsigned 8 bit int
	saturated = (type==8) ? 127 : saturated			// signed 8 bit int
	if (saturated)
		String wName = UniqueName("mat",1,0)
		Variable i0=b.xlo,i1=b.xhi,j0=b.ylo,j1=b.yhi
		Make/N=(i1-i0+1,j1-j0+1)/B $wName
		Wave roi=$wName
		roi = image[p+i0][q+j0]>=saturated
		saturated = sum(roi)
		KillWaves/Z roi
	endif

	if (ItemsInList(GetRTStackInfo(0))<2)
		Variable N=(b.xhi-b.xlo+1)*(b.yhi-b.ylo+1)
		printf "\rfor  '%s'[%d,%d][%d,%d]\r",imageName,b.xlo,b.xhi,b.ylo,b.yhi
		if (V_npnts>0)
			printf "	min at [%d,%d] = %g,   max at [%d,%d] = %g\r",V_minRowLoc,V_minColLoc,V_min,V_maxRowLoc,V_maxColLoc,V_max
			printf "	avg = %g, with an std dev = %g,   and an avg deviation = %g\r",V_avg,V_sdev,V_adev
		endif
		if (N==V_npnts)
			printf "	All %d points are valid, with a sum = %g,  and an rms = %g\r",V_npnts,Vsum,V_rms
		elseif (V_npnts>0)
			printf "	There are %d valid points (out of %d), with a sum = %g,  and an rms = %g\r",V_npnts,N,Vsum,V_rms
		else
			printf "	There are %d valid points (out of %d)\r",V_npnts,N
		endif
		if (saturated)
			printf "	There are %d saturated pixels\r",saturated
		endif
	endif
	if (V_npnts<1)
		return ""											// none of the following statistics mean anything
	endif
	String str=""
	str = ReplaceNumberByKey("V_adev",str,V_adev,"=")
	str = ReplaceNumberByKey("V_avg",str,V_avg,"=")
	str = ReplaceNumberByKey("V_sum",str,Vsum,"=")
	str = ReplaceNumberByKey("V_kurt",str,V_kurt,"=")
	str = ReplaceNumberByKey("V_min",str,V_min,"=")
	str = ReplaceNumberByKey("V_minColLoc",str,V_minColLoc,"=")
	str = ReplaceNumberByKey("V_minRowLoc",str,V_minRowLoc,"=")
	str = ReplaceNumberByKey("V_max",str,V_max,"=")
	str = ReplaceNumberByKey("V_maxColLoc",str,V_maxColLoc,"=")
	str = ReplaceNumberByKey("V_maxRowLoc",str,V_maxRowLoc,"=")
	str = ReplaceNumberByKey("V_npnts",str,V_npnts,"=")
	str = ReplaceNumberByKey("V_rms",str,V_rms,"=")
	str = ReplaceNumberByKey("V_sdev",str,V_sdev,"=")
	str = ReplaceNumberByKey("V_skew",str,V_skew,"=")
	str = ReplaceNumberByKey("N_saturated",str,saturated,"=")
	return str
End


Static Function ModifyOnly_ctab_range(gName,imageName,lo,hi)
	// modify the range of the current color table, does not change the color table, or reverse
	String gName				// optional graph name
	String imageName			// name of particular image on graph
	Variable lo,hi				// range of color table to set (if NaN, then use auto-scale)

	String infoStr = ImageInfo(gName,imageName,0)
	String ctab = StringByKey("ctab",StringByKey("RECREATION",infoStr),"=")
	String colorTable = StringFromList(2,ctab,",")
	Variable colorRev = str2num(StringFromList(3,ctab,","))

	if (numtype(lo) && !numtype(hi))
		ModifyImage/W=$gName $imageName ctab= {*,hi,$colorTable,colorRev}
	elseif (!numtype(lo) && numtype(hi))
		ModifyImage/W=$gName $imageName ctab= {lo,*,$colorTable,colorRev}
	elseif (numtype(lo) && numtype(hi))
		ModifyImage/W=$gName $imageName ctab= {*,*,$colorTable,colorRev}
	else
		ModifyImage/W=$gName $imageName ctab= {lo,hi,$colorTable,colorRev}
	endif
End

Function/T get_ctab_range(gName,imageName,lo,hi)
	// modify the range of the current color table, does not change the color table, or reverse
	String gName				// optional graph name
	String imageName			// name of particular image on graph
	Variable &lo,&hi			// range of color table to set
	String infoStr = ImageInfo(gName,imageName,0)
	String ctab = StringByKey("ctab",StringByKey("RECREATION",infoStr),"=")
	ctab = ctab[strsearch(ctab, "{",0)+1,Inf]
	lo = str2num(StringFromList(0,ctab,","))
	hi = str2num(StringFromList(1,ctab,","))
	return StringFromList(2,ctab,",")
End


Function MatchMarqueeAspect()
	// change the size of a marquee to match the aspect ratio of the plot (gives square pixels if you zoom)
	GetWindow kwTopWin psize
	Variable aspectPlot, aspectMarquee, add
	aspectPlot = (V_bottom-V_top) / (V_right-V_left)			// aspect > 1 --> tall skinny
	GetMarquee/Z
	aspectMarquee = (V_bottom-V_top) / (V_right-V_left)
	if (aspectMarquee > aspectPlot)								// make marquee wider, increase width symmetrically
		add = ((V_bottom-V_top)/aspectPlot+V_left-V_right)/2
		V_left -= add
		V_right += add
	else															// make marquee taller, increase height symmetrically
		add = (aspectPlot*(V_right-V_left)+V_top-V_bottom)/2
		V_top -= add
		V_bottom += add
	endif
	SetMarquee V_left,V_top,V_right,V_bottom
End


Function pointRangeOfMarquee(gName,imageName,b)
	String gName
	String imageName
	STRUCT box &b

	Wave image = ImageNameToWaveRef(gName,imageName)
	String infoStr = ImageInfo("",imageName,0)
	GetMarquee/Z $StringByKey("YAXIS",infoStr), $StringByKey("XAXIS",infoStr)
	if (!V_flag)						// box not set, bad marquee
		b.xlo=NaN ; b.xhi=NaN ; b.ylo=NaN ; b.yhi=NaN
		return 1
	endif

	Variable Nx=DimSize(image,0), Ny=DimSize(image,1)
	V_left = round(V_left-DimOffset(image,0))/DimDelta(image,0)
	V_right = round((V_right-DimOffset(image,0))/DimDelta(image,0))
	V_top = round((V_top-DimOffset(image,1))/DimDelta(image,1))
	V_bottom = round((V_bottom-DimOffset(image,1))/DimDelta(image,1))
	V_left = limit(V_left,0,Nx-1)
	V_right = limit(V_right,0,Nx-1)
	V_top = limit(V_top,0,Ny-1)
	V_bottom = limit(V_bottom,0,Ny-1)

	b.xlo=min(V_left,V_right)		;	b.xhi=max(V_left,V_right)
	b.ylo=min(V_top,V_bottom)	;	b.yhi=max(V_top,V_bottom)
	return 0
End


Function/T sumImageToBothEdges(image,[printIt,planeNum])
	Wave image									// 2D wave to sum left and down
	Variable printIt
	Variable planeNum
	printIt = ParamIsDefault(printIt) ? 0 : printIt
	planeNum = ParamIsDefault(planeNum) ? -1 : planeNum

	if (!WaveExists(image))
		String iName = StringFromList(0,ImageNameList("",";"))
		if (strlen(iName))						// first try to get image from top graph
			DoAlert 2, "creates vertical and horizontal sums using '"+iName+"'?"
			if (V_flag==1)
				Wave image = ImageNameToWaveRef("",iName)
			elseif (V_flag==3)
				return ""
			else
				iName = ""
			endif
		endif
		printIt = 1
	endif
	if (!WaveExists(image))					// next try the current data folder
		if (strlen(WaveList("*",";","DIMS:2"))<1)
			return ""
		endif
		Prompt iName,"image for vertical and horizontal sums",popup,WaveList("*",";","DIMS:2")
		DoPrompt "image",iName
		if (V_flag)
			return ""
		endif
		Wave image = $iName
		printIt = 1
	endif
	if (!WaveExists(image))					// give up
		return ""
	endif
	if (printIt)
		printf "¥sumImageToBothEdges(%s)\r",GetWavesDataFolder(image,2)
	endif

	if (DimSize(image,2)>1 && !(planeNum>=0))
		planeNum = imagePlaneDisplayed("",image)
		planeNum = planeNum>=0 ? planeNum : -1
	endif

	if (planeNum>=0)
		ImageTransform/P=(planeNum) sumAllCols image
		ImageTransform/P=(planeNum) sumAllRows image
	else
		ImageTransform sumAllCols image
		ImageTransform sumAllRows image
	endif
	SetScale/P x DimOffset(image,0),DimDelta(image,0),WaveUnits(image,0), W_sumRows
	SetScale/P x DimOffset(image,1),DimDelta(image,1),WaveUnits(image,1), W_sumCols
	String str=AxisLabelFromGraph("",image,"HORIZONTAL")
	if (strlen(str))
		str = ReplaceString(";",str,"_")										// a semi-colon will mess up the list
		Note/K W_sumRows, ReplaceStringByKey("GraphAxisLabelHoriz",note(W_sumRows),str,"=")
	endif
	str=AxisLabelFromGraph("",image,"VERTICAL")
	if (strlen(str))
		str = ReplaceString(";",str,"_")										// a semi-colon will mess up the list
		String wnote = note(W_sumCols)
		wnote = RemoveByKey("GraphAxisLabelVert",wnote,"=")				// vert label on image
		wnote = ReplaceStringByKey("GraphAxisLabelHoriz",wnote,str,"=")	// is horizontal label on W_sumCols
		Note/K W_sumCols, wnote
	endif
	String fldr = GetWavesDataFolder(image,1)
	String wname = fldr+NameOfWave(image)+"_SumAllY"
	Duplicate/O W_sumRows $wname
	String out=wname+";"
	wname = fldr+NameOfWave(image)+"_SumAllX"
	Duplicate/O W_sumCols $wname
	out += wname+";"
	KillWaves/Z W_sumRows,W_sumCols
	if(ItemsInList(GetRTStackInfo(0))<2 || printIt)
		printf "   summed the image  '%s'",NameOfWave(image)
		if (planeNum>=0)
			printf "[%g]", planeNum
		endif
		printf "  left and down to get waves '%s' ,  and  '%s'\r",StringFromList(0,out),StringFromList(1,out)
	endif
	return out
End


Function imagePlaneDisplayed(win,image)		// returns plane displayed when showing one layer of a 3D array as an image
	String win									// window, use "" for top graph window
	Wave image									// 2D wave to sum left and down

	Variable found=0
	if (WaveExists(image))
		if (strlen(win)<1)
			win = StringFromList(0,WindowsWithWave(image,1))
		endif
		if (strlen(win)<1)						// image was specified, but it is not displayed
			return NaN
		endif
	endif

	if (!WaveExists(image))					// look for first image on win, probably image was not specified
		String iList = ImageNameList(win,";")
		Wave image = ImageNameToWaveRef(win, StringFromList(0,iList) )
	endif
	if (!WaveExists(image))					// could not find image
		return NaN
	endif

	Variable planeNum=0
	if (DimSize(image,2)>1)					// extra planes present
		String str = WinRecreation("", 0), find = "ModifyImage "+NameOfWave(image)+" plane="
		Variable i=strsearch(str,find,0,2)
		if (i>0)
			i += strlen(find)
			planeNum = str2num(str[i,Inf])
		endif
		planeNum = planeNum>=0 ? planeNum : 0
	endif
	return planeNum
End


// This was moved to Utility_JZT
//
//Function DrawMarker(x0,y0,dx,dy,style,[color,thick,dash,win,layer])
//	Variable x0,y0,dx,dy
//	String style
//	String color							// one of "red", "blue", "green", "yellow", "magenta", "cyan", or a triplet like "1,1,0" or "65535,0,0"
//	Variable thick						// line thickness, defaults to 0.50 (thin)
//	Variable dash
//	String win							// optional name of window
//	String layer						// Drawing layer to use.  Default is "UserFront"
//	if (numtype(x0++y0+dx+dy) || dx<=0 || dy<=0)
//		return 1
//	endif
//	if (ParamIsDefault(color))
//		color = SelectString(stringmatch(style,"BoxWithTicks"),"black","30583,30583,30583")	// special default color for BoxWithTicks
//	endif
////	color = SelectString(ParamIsDefault(color),color,"black")	// default color to black
////	color = SelectString(ParamIsDefault(color) && stringmatch(style,"BoxWithTicks"),color,"30583,30583,30583")	// special default color for BoxWithTicks
//	thick = ParamIsDefault(thick) ? 0.5 : thick
//	dash = ParamIsDefault(dash) ? 0 : dash
//	if (ParamIsDefault(win))
//		win = ""
//	endif
//	if (ParamIsDefault(layer))
//		layer = "UserFront"
//	endif
//	String layerList = "ProgBack;UserBack;ProgAxes;UserAxes;ProgFront;UserFront"
//	layer = SelectString(strlen(layer),"UserFront",layer)
//	if (WhichListItem(layer,layerList)<0)
//		printf "ERROR -- layer = '%s' is invlaid, must be one of '%s'\r",layer,layerList
//		return 1
//	endif
//
//	String rgb = StringByKey(color,"red:1,0,0;blue:0,0,1;green:0,1,0;yellow:1,1,0;magenta:1,0,1;cyan:0,1,1;black:0,0,0;white:1,1,1")	// name -> numbers
//	rgb = SelectString(strlen(rgb),color,rgb)
//	Variable r=str2num(StringFromList(0,rgb,",")), g=str2num(StringFromList(1,rgb,",")), b=str2num(StringFromList(2,rgb,","))
//	r = !(r>=0) ? 0 : r											// default invalid to black
//	g = !(g>=0) ? 0 : g
//	b = !(b>=0) ? 0 : b
//	if (r+g+b<=3)												// if numbers are all small change from [0,1] -> [0,65535]
//		r *= 65535
//		g *= 65535
//		b *= 65535
//	endif
//
//	style = SelectString(strlen(style),"cross",style)			// defaults to Cross
//	String list = ImageInfo(win,"",0)
//	if (strlen(list)<1)											// perhaps no image, only a graph
//		list = TraceInfo(win,StringFromList(0,TraceNameList(win,";",1+4)),0)
//	endif
//	String xaxis=StringByKey("XAXIS",list), yaxis=StringByKey("YAXIS",list)
//	String aList = AxisList(win)
//	aList = RemoveFromList(yaxis,aList)
//	if (strlen(xaxis)<1)										// I am getting desperate, last thing to try
//		Variable i
//		for (i=0;i<ItemsInList(aList);i+=1)
//			xaxis = StringFromList(i,aList)
//			if (strsearch(xaxis,"bottom",0,2)>=0)
//				break
//			elseif (strsearch(xaxis,"top",0,2)>=0)
//				break
//			endif
//		endfor
//	endif
//	if (strlen(yaxis)<1)										// I am getting desperate, last thing to try
//		aList = RemoveFromList(xaxis,aList)
//		for (i=0;i<ItemsInList(aList);i+=1)
//			yaxis = StringFromList(i,aList)
//			if (strsearch(yaxis,"left",0,2)>=0)
//				break
//			elseif (strsearch(yaxis,"right",0,2)>=0)
//				break
//			endif
//		endfor
//	endif
//
//	dy = abs(dy/2)		// change FW to HW
//	dx = abs(dx/2)
//	SetDrawLayer/W=$win $layer
//
//	Variable err = 0
//	if (stringmatch(style,"cross gap"))
//		SetDrawEnv/W=$win  xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick, dash=dash
//		DrawLine/W=$win  x0,(y0-dy),x0,(y0-0.1*dy)
//		SetDrawEnv/W=$win  xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick, dash=dash
//		DrawLine/W=$win  x0,(y0+0.1*dy),x0,(y0+dy)
//		SetDrawEnv/W=$win xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick, dash=dash
//		DrawLine/W=$win (x0-dx),y0,(x0-0.1*dx),y0
//		SetDrawEnv/W=$win xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick, dash=dash
//		DrawLine/W=$win (x0+0.1*dx),y0,(x0+dx),y0
//	elseif (stringmatch(style,"cross"))
//		SetDrawEnv/W=$win  xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick, dash=dash
//		DrawLine/W=$win  (x0),(y0-dy),x0,(y0+dy)
//		SetDrawEnv/W=$win xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick, dash=dash
//		DrawLine/W=$win (x0-dx),y0,(x0+dx),y0
//	elseif (stringmatch(style,"X"))
//		SetDrawEnv/W=$win xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick, dash=dash
//		DrawLine/W=$win (x0-dx),(y0-dy),(x0+dx),(y0+dy)
//		SetDrawEnv/W=$win xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick, dash=dash
//		DrawLine/W=$win (x0+dx),(y0-dy),(x0-dx),(y0+dy)
//	elseif (stringmatch(style,"BoxWithTicks"))
//		Variable hw							// half width of box in pixels, usually ~10
//		SetDrawEnv/W=$win xcoord= $xaxis,ycoord= $yaxis,fillpat= 0, linethick= thick, dash=dash
//		DrawRect x0-dx,y0-dy,x0+dx,y0+dy
//		SetDrawEnv/W=$win xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick, dash=dash
//		DrawLine/W=$win x0,y0-dy,x0,y0-dy/2
//		SetDrawEnv/W=$win xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick, dash=dash
//		DrawLine/W=$win x0,y0+dy/2,x0,y0+dy
//		SetDrawEnv/W=$win xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick, dash=dash
//		DrawLine/W=$win x0-dx,y0,x0-dx/2,y0
//		SetDrawEnv/W=$win xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick, dash=dash
//		DrawLine/W=$win x0+dx/2,y0,x0+dx,y0
//	else
//		err = 1
//	endif
//	SetDrawLayer/W=$win UserFront
//	return err
//End

// ==================================== End of Generic Util ===================================
// =============================================================================================



// =============================================================================================
// ================================= Start of Peak Evaluation ==================================

Function AddGaussianToPeakList()
	GetMarquee/Z
	if (V_flag==0)
		DoAlert 0, "no Marquee selected, do nothing"
		return 1
	endif

	Wave image = ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";")))
	Wave peakList = $StringByKey("FullPeakList",GetUserData("","","FitPeaks"),"=")
#if exists("FindPeakListForImage")
	if (!WaveExists(peakList))						// not in user data, search for one
		Wave peakList = FindPeakListForImage(image)
	endif
#endif
#if exists("MakeEmptyPeakListForImage")
	if (!WaveExists(peakList))						// There is no peakList, make a new one
		Wave peakList = MakeEmptyPeakListForImage(image,ask=1,keyVals="FittedPeakShape=Gaussian")
	endif
#endif
	if (!WaveExists(peakList) || !WaveExists(image))
		return 1
	endif
	if (!WaveInClass(peakList,"FittedPeakList*;"))
		return 1
	endif
	if (Gaussian2dPeakCenter())
		return 1
	endif
	Wave W_coef=W_coef, W_sigma=W_sigma
	Variable fw0= ( 2*sqrt(2*ln(2)) )			// this factor corrects for 2-d sigma to FWHM, apply to K3, K5, and the corresponding errors
	Variable px=W_coef[2], py=W_coef[4], amp=W_coef[1]
	Variable fwx=W_coef[3]*fw0, fwy=W_coef[5]*fw0, xcorr=W_coef[6]

	// check for a duplicate within 1.5 pixel
	Variable N=DimSize(peakList,0), ipeak
	for (ipeak=0;ipeak<N;ipeak+=1)
		if (abs(peakList[ipeak][0]-px)<1.5 && abs(peakList[ipeak][1]-py)<1.5)
			break
		endif
	endfor
	String str
	if (ipeak>=N)
		sprintf str,"Add peak at (%.1f, %.1f)  with FWHM=(%.1f, %.1f) & x-corr=%.2f, Amp=%g",px,py,fwx,fwy,xcorr,amp
	else
		sprintf str,"REPLACE peak %d with:\rpeak at (%.1f, %.1f)  with FWHM=(%.1f, %.1f) & x-corr=%.2f, Amp=%g",ipeak,px,py,fwx,fwy,xcorr,amp
	endif
	Cursor/I/H=1 J $NameOfWave(image)  px,py
	DoAlert 1,str
	Cursor/K J
	if (V_flag!=1)
		return 1
	endif

	if (ipeak>=N)
		Redimension/N=(ipeak+1,-1) peakList
	endif
	peakList[ipeak][0] = px;				peakList[ipeak][1] = py
	peakList[ipeak][2] = W_sigma[2];		peakList[ipeak][3] = W_sigma[4]
	peakList[ipeak][4] = fwx;				peakList[ipeak][5] = fwy
	peakList[ipeak][6] = W_sigma[3]*fw0;	peakList[ipeak][7] = W_sigma[5]*fw0
	peakList[ipeak][8] = xcorr;				peakList[ipeak][9] = W_sigma[6]
	peakList[ipeak][10] = amp* fwx*fwy

	if (NumberByKey("boxesOn"GetUserData("","","boxes"),"="))
		DrawMarker(px,py,20,20,"BoxWithTicks")
	endif
	KillWaves/Z W_coef,W_sigma
	return 0
End
//Function AddGaussianToPeakList()
//	GetMarquee/Z
//	if (V_flag==0)
//		DoAlert 0, "no Marquee selected, do nothing"
//		return 1
//	endif
//
//	Wave image = ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";")))
//	Wave peakList = $StringByKey("FullPeakList",GetUserData("","","FitPeaks"),"=")
//	if (!WaveExists(peakList))						// not in user data, search for one
//		Wave peakList = FindPeakListForImage(image)
//	endif
//	if (!WaveExists(peakList) || !WaveExists(image))
//		return 1
//	endif
//	if (!WaveInClass(peakList,"FittedPeakList*;"))
//		return 1
//	endif
//	if (Gaussian2dPeakCenter())
//		return 1
//	endif
//	Wave W_coef=W_coef, W_sigma=W_sigma
//
//	Variable fw= ( 2*sqrt(2*ln(2)) )			// this factor corrects for 2-d sigma to FWHM, apply to K3, K5, and the corresponding errors
//	String str
//	sprintf str,"Add peak at (%.1f, %.1f)  with FWHM=(%.1f, %.1f) & x-corr=%.2f, Amp=%g",W_coef[2],W_coef[4],W_coef[3]*fw,W_coef[5]*fw,W_coef[6],W_coef[1]
//	Cursor/I/H=1 J $NameOfWave(image)  W_coef[2],W_coef[4]
//	DoAlert 1,str
//	Cursor/K J
//	if (V_flag!=1)
//		return 1
//	endif
//
//	Variable N=DimSize(peakList,0), M=DimSize(peakList,1)
//	Redimension/N=(N+1,M) peakList
//	peakList[N][0] = W_coef[2];		peakList[N][1] = W_coef[4]
//	peakList[N][2] = W_sigma[2];		peakList[N][3] = W_sigma[4]
//	peakList[N][4] = W_coef[3] * fw;	peakList[N][5] = W_coef[5] * fw
//	peakList[N][6] = W_sigma[3]*fw;	peakList[N][7] = W_sigma[5]*fw
//	peakList[N][8] = W_coef[6];		peakList[N][9] = W_sigma[6]
//	peakList[N][10] = W_coef[1]*(W_coef[3]*fw)*(W_coef[5]*fw)
//
//	if (NumberByKey("boxesOn"GetUserData("","","boxes"),"="))
//		DrawMarker(W_coef[2],W_coef[4],20,20,"BoxWithTicks")
//	endif
//	KillWaves/Z W_coef,W_sigma
//	return 0
//End
//
Function AddPeakAtCursorA()
	if (strlen(CsrWave(A))<=0)				// This requres that cursor A be on the image
		DoAlert 0, "Cursor A not on the Image, do nothing"
		return 1
	endif

	Wave image = ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";")))
	Wave peakList = $StringByKey("FullPeakList",GetUserData("","","FitPeaks"),"=")

#if exists("FindPeakListForImage")
	if (!WaveExists(peakList))						// not in user data, search for one
		Wave peakList = FindPeakListForImage(image)
	endif
#endif
	if (!WaveExists(peakList) || !WaveExists(image))
		return 1
	endif
	if (!WaveInClass(peakList,"FittedPeakList*;"))
		return 1
	endif

	Variable x0=hcsr(A), y0=vcsr(A), z0=zcsr(A)
	if (numtype(x0+y0))
		return 1
	endif

	String str
	sprintf str,"Add peak at (%d, %d),   Amp=%g",x0,y0,z0
	DoAlert 1,str
	if (V_flag!=1)
		return 1
	endif

	Variable N=DimSize(peakList,0), M=DimSize(peakList,1)
	Redimension/N=(N+1,M) peakList
	peakList[N][] = NaN									// most of the columns cannot be set without fitting a peak
	peakList[N][0] = x0
	peakList[N][1] = y0
	peakList[N][10] = z0
	if (NumberByKey("boxesOn"GetUserData("","","boxes"),"="))
		DrawMarker(x0,y0,20,20,"BoxWithTicks")
	endif
	return 0
End
//
Function RemovePeakFromPeakList()				// returns row number removed, -1 if nothing done
	Variable V_left, V_right, V_bottom, V_top
	String xaxis=StringByKey("XAXIS",ImageInfo("","",0)), yaxis=StringByKey("YAXIS",ImageInfo("","",0))
	GetMarquee/Z $yaxis, $xaxis
	if (V_flag==0)
		DoAlert 0, "no Marquee selected, do nothing"
		return -1
	endif

	Wave image = ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";")))
	Wave peakList = $StringByKey("FullPeakList",GetUserData("","","FitPeaks"),"=")
#if exists("FindPeakListForImage")
	if (!WaveExists(peakList))						// not in user data, search for one
		Wave peakList = FindPeakListForImage(image)
	endif
#endif
	if (!WaveExists(peakList) || !WaveExists(image))
		return -1
	endif

	Variable tol, x0,y0
	tol = round(min(V_right-V_left,V_bottom-V_top)/2)	// ± tolerance in pixels
	x0 = (V_left+V_right)/2
	y0 = (V_bottom+V_top)/2

	// find closest peak
	Variable i,iBest=-1,distBest=Inf, dist
	for (i=0;i<DimSize(peakList,0);i+=1)
		dist = sqrt((peakList[i][0]-x0)^2 + (peakList[i][1]-y0)^2)
		if (dist < distBest && dist<tol)
			iBest = i
			distBest = dist
		endif
	endfor
	if (i<0)
		return -1
	endif

	String str
	sprintf str, "Remove the point at (%d, %d)?",peakList[iBest][0],peakList[iBest][1]
	Cursor/P/I/H=1 J $NameOfWave(image)  peakList[iBest][0],peakList[iBest][1]
	DoAlert 1, str
	Cursor/K J
	if (V_flag!=1)
		return -1
	endif
	DeletePoints/M=0 iBest, 1, peakList
	return iBest
End


Function Gaussian2dPeakCenter([printIt])	// returns 1 if error,  0 is OK
	// the result is passed to calling function by W_coef or {K0, K1, ... K6}
	// this fits a 2-d Gaussian to a region on an image.  The region is selected with a marquee in an image plot, or specified
	Variable printIt							// enables print out of result
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? (strlen(GetRTStackInfo(2))==0) : !(!printIt)

	if (WinType("")!=1)
		Abort "Top window must be a graph, with region selected"
	endif
	String xaxis=StringByKey("XAXIS",ImageInfo("","",0)), yaxis=StringByKey("YAXIS",ImageInfo("","",0))
	GetMarquee/Z $yaxis, $xaxis					// this sets left, right, . . .
	Variable left=V_left, right=V_right, bottom=V_bottom, top=V_top
	String win=S_marqueeWin						// name of top graph window
	if (V_flag==0)
		DoAlert 0, "No region selected on the graph"
		return 1
	endif
	if (right==0 && left==0 || top==0 && bottom==0)
		DoAlert 0, "Unable get range for horizonal or vertical axis info for the selected graph"
		return 1
	endif

	// Find the wave that is the image (assume that it is the only 2-d wave on the graph)
	String imageName=StringFromList(0,ImageNameList(win,";"))	// name of wave with image
	if (strlen(imageName)<1)
		DoAlert 0, "Unable to find image on top graph"
		return 1
	endif
	Wave imageDisplayed=ImageNameToWaveRef(win,imageName)
	Variable Nx=DimSize(imageDisplayed,0), Ny=DimSize(imageDisplayed,1)

	if (WaveDims(imageDisplayed)>2)
		Duplicate/FREE imageDisplayed, image
		Redimension/N=(Nx,Ny) image
		image = imageDisplayed[p][q][0]
	else
		Wave image = imageDisplayed
	endif
	String info=ImageInfo(win,imageName,0)
	Wave xPlot=$(StringByKey("XWAVEDF",info)+StringByKey("XWAVE",info))
	Wave yPlot=$(StringByKey("YWAVEDF",info)+StringByKey("YWAVE",info))
	if (WaveExists(xPlot))
		Duplicate/FREE xPlot, xw
		Redimension/N=(Nx) xw
		xw = (xPlot[p]+xPlot[p+1]) / 2
	else
		Wave xw=$""
	endif
	if (WaveExists(yPlot))
		Duplicate/FREE yPlot, yw
		Redimension/N=(Ny) yw
		yw = (yPlot[p]+yPlot[p+1]) / 2
	else
		Wave yw=$""
	endif

	// do the Gaussian fit
	Variable i = FitOneGaussianPeak(image,left,right,bottom,top, xw=xw, yw=yw)	// returns 1 on error,  0 is OK
	if (i)
		return i
	endif

	Wave W_coef=W_coef, W_sigma=W_sigma
	Variable x0=W_coef[2], y0=W_coef[4], cor=W_coef[6]
	Variable fw=( 2*sqrt(2*ln(2)) )			// this factor corrects for 2-d sigma to FWHM, apply to K3, K5, and the corresponding errors
	Variable dx=W_coef[3]*fw, dy=W_coef[5]*fw
	if (printIt)
		printf "Gaussian peak in '%s' at (%s, %s) with FWHM of (%s, %s) and correlation of %s  Amp=%.3g\r",ImageName, ValErrStr(x0,W_sigma[2]), ValErrStr(y0,W_sigma[4]), ValErrStr(dx,W_sigma[3]*fw), ValErrStr(dy,W_sigma[5]*fw), ValErrStr(cor,W_sigma[6]), W_coef[1]
		DoAlert 1, "Draw a cross at result"
		if (V_Flag==1)
			DrawMarker(x0,y0,abs(right-left),abs(top-bottom),"Cross",win=win)
		endif
	endif
	return 0
End
//Function Gaussian2dPeakCenter(printON)	// returns 1 if error,  0 is OK
//// the result is passed to calling function by W_coef or {K0, K1, ... K6}
//// this fits a 2-d gaussina to a region on an image.  The region is selected with a marquee in an image plot, or specified
//	Variable printON							// enables print out of result
//
//	Variable V_left, V_right, V_bottom, V_top
//	if (WinType("")!=1)
//		Abort "Top window must be a graph, with region selected"
//	endif
//	String xaxis=StringByKey("XAXIS",ImageInfo("","",0)), yaxis=StringByKey("YAXIS",ImageInfo("","",0))
//	GetMarquee/Z $yaxis, $xaxis					// this sets V_left, V_right, . . .
//	String win = S_marqueeWin						// name of top graph window
//	if (V_flag==0)
//		DoAlert 0, "No region selected on the graph"
//		return 1
//	endif
//	if (V_right==0 && V_left==0 || V_top==0 && V_bottom==0)
//		DoAlert 0, "Unable get range for horizonal or vertical axis info for the selected graph"
//		return 1
//	endif
//
//	// Find the wave that is the image (assume that it is the only 2-d wave on the graph
//	String imageName = StringFromList(0,ImageNameList(win,";"))	// name of wave with image
//	if (strlen(ImageName)<1)
//		DoAlert 0, "Unable to find image on top graph"
//		return 1
//	endif
//	Wave image = ImageNameToWaveRef(win,imageName)
//
//	// set up for the gaussian fit
//	Variable i = FitOneGaussianPeak(image,V_left,V_right,V_bottom,V_top)	// returns 1 if error,  0 is OK
//	if (i)
//		return i
//	endif
//
//	Wave W_coef=W_coef
//	Wave W_sigma=W_sigma
//	Variable x0,y0,dx,dy,cor
//	Variable fw= ( 2*sqrt(2*ln(2)) )			// this factor corrects for 2-d sigma to FWHM, apply to K3, K5, and the corresponding errors
//	x0 = W_coef[2]
//	dx = W_coef[3] * fw
//	y0 = W_coef[4]
//	dy = W_coef[5] * fw
//	cor = W_coef[6]
//	if (printON)
//		Variable pixels = (DimDelta(image,0)==1 && DimDelta(image,1)==1 && DimOffset(image,0)==0 && DimOffset(image,1)==0)
//		Variable places1, places2, places3
//		if (pixels)
//			places1 = 2
//			places2 = 2
//			places3 = 2
//		else
//			Variable placesX = limit(ceil(log(limit(1/abs(W_sigma[2]),0.1,1e12))),0,12)
//			Variable placesY = limit(ceil(log(limit(1/abs(W_sigma[4]),0.1,1e12))),0,12)
//			places1 = max(placesX,placesY)
//			placesX = limit(ceil(log(limit(1/abs(W_sigma[3]*fw),0.1,1e12))),0,12)
//			placesY = limit(ceil(log(limit(1/abs(W_sigma[5]*fw),0.1,1e12))),0,12)
//			places2 = max(placesX,placesY)
//			places3 = limit(ceil(log(limit(1/abs(W_sigma[6]),0.1,1e12))),0,12)
//		endif
//		String fmt
//		sprintf fmt, "Gaussian peak in '%%s' at (%%.%df±%%.%df, %%.%df±%%.%df) with FWHM of (%%.%df±%%.%df, %%.%df±%%.%df) and correlation of %%.%df±%%.%df  Amp=%%.3g\r",places1,places1,places1,places1,places2,places2,places2,places2,places3,places3
//		printf fmt,ImageName,x0,W_sigma[2],y0,W_sigma[4],dx,W_sigma[3]*fw,dy,W_sigma[5]*fw,cor,W_sigma[6],W_coef[1]
//		DoAlert 1, "Draw a cross at result"
//		if (V_Flag==1)
//			DrawMarker(x0,y0,abs(V_right-V_left),abs(V_top-V_bottom),"Cross",win=win)
//		endif
//	endif
//	return 0
//End


Function FitOneGaussianPeak(image,left,right,bottom,top,[xw,yw,printIt])	// returns 1 if error,  0 is OK
// the result is passed to calling function by W_coef or {K0, K1, ... K6}
// this fits a 2-d Gaussian to a region on an image.  The region is selected with a marquee in an image plot, or specified
	Wave image									// wave with image to fit
	Variable left, right, bottom, top	// range to use for fit
	Wave xw,yw									// x & y waves giving the x & y axes for the image
	Variable printIt							// enables print out of result
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? (strlen(GetRTStackInfo(2))==0) : !(!printIt)
	if (!WaveExists(image))
		DoAlert 0, "Unable to find image to fit"
		return -1
	endif

	Variable Nx=DimSize(image,0), Ny=DimSize(image,1)
	if (WaveExists(xw))
		if (numpnts(xw)!=Nx)
			Wave xw = $""						// for invalid xw, set to null wave
		endif
	endif
	if (WaveExists(yw))
		if (numpnts(yw)!=Ny)
			Wave yw = $""						// for invalid yw, set to null wave
		endif
	endif

	Variable pixels=1							// using image in pixels
	pixels = (DimDelta(image,0)==1 && DimDelta(image,1)==1 && DimOffset(image,0)==0 && DimOffset(image,1)==0)
	if (pixels)
		Variable swap
		if (top<bottom)						// insist that bottom < top
			swap = top
			top = bottom
			bottom = swap
		endif
		if (right<left)						// insist that left < right
			swap = left
			left = right
			right = swap
		endif
		left = limit(floor(left),0,DimSize(image,0)-1)	// these are pixels, so they need to be integers
		right = limit(ceil(right),0,DimSize(image,0)-1)
		bottom = limit(floor(bottom),0,DimSize(image,1)-1)
		top = limit(ceil(top),0,DimSize(image,1)-1)
	endif

	// set up for the Gaussian fit
	Make/N=7/O/D W_coef, W_sigma=W_sigma
	Variable V_FitOptions=4, V_FitError=0,V_FitQuitReason=0
	//	Execute "SetIgorOption UseVeclib=0"
	Make/T/FREE T_Constraints={"K3 > 0","K5 > 0","K6 > 0","K6 < 1"}
	if (pixels)
		CurveFit/Q/N Gauss2D image[left,right][bottom,top]/C=T_Constraints				// pixels
	elseif (WaveExists(xw) && WaveExists(yw))
		CurveFit/Q/N Gauss2D image(left,right)(bottom,top)/C=T_Constraints/X=xw/Y=yw	// scaled with x & y waves
	elseif (WaveExists(xw))
		CurveFit/Q/N Gauss2D image(left,right)(bottom,top)/C=T_Constraints/X=xw			// scaled with only x wave
	elseif (WaveExists(yw))
		CurveFit/Q/N Gauss2D image(left,right)(bottom,top)/C=T_Constraints/Y=yw			// scaled with only y wave
	else
		CurveFit/Q/N Gauss2D image(left,right)(bottom,top)/C=T_Constraints				// scaled with image scaling
	endif
	if (printIt && V_FitError)
		print FitErrorString(V_FitError,V_FitQuitReason)
	endif
	if (!V_FitError)
		Wave W_coef = W_coef
		Variable zero = -1e-15, one=1+1e-12					// almost zero, almost one
		V_FitError += (W_coef[3]<zero || W_coef[5]<zero || W_coef[6]<zero || W_coef[6]>one) ? 32 : 0			//  constraints failed, set bit 5
	endif
	return V_FitError
End
//Function FitOneGaussianPeak(Image,V_left,V_right,V_bottom,V_top,[Q])	// returns 1 if error,  0 is OK
//// the result is passed to calling function by W_coef or {K0, K1, ... K6}
//// this fits a 2-d gaussina to a region on an image.  The region is selected with a marquee in an image plot, or specified
//	Wave image									// wave with image to fit
//	Variable V_left, V_right, V_bottom, V_top	// range to use for fit
//	Variable Q									// quiet flag
//
//	Q = ParamIsDefault(Q) ? 0 : Q
//	if (!WaveExists(image))
//		DoAlert 0, "Unable to find image to fit"
//		return -1
//	endif
//
//	Variable pixels=1					// using image in pixels
//	pixels = (DimDelta(image,0)==1 && DimDelta(image,1)==1 && DimOffset(image,0)==0 && DimOffset(image,1)==0)
//	if (pixels)
//		Variable swap
//		if (V_top<V_bottom)			// insist that bottom < top
//			swap = V_top
//			V_top = V_bottom
//			V_bottom = swap
//		endif
//		if (V_right<V_left)				// insist that left < right
//			swap = V_left
//			V_left = V_right
//			V_right = swap
//		endif
//		V_left = limit(floor(V_left),0,DimSize(image,0)-1)	// these are pixels, so they need to be integers
//		V_right = limit(ceil(V_right),0,DimSize(image,0)-1)
//		V_bottom = limit(floor(V_bottom),0,DimSize(image,1)-1)
//		V_top = limit(ceil(V_top),0,DimSize(image,1)-1)
//	endif
//
//	// set up for the gaussian fit
//	Make/N=7/O/D W_coef, W_sigma=W_sigma
//	Variable V_FitOptions=4, V_FitError=0,V_FitQuitReason=0
//	V_FitOptions=4
//	//	Execute "SetIgorOption UseVeclib=0"
//	Make/O/T T_Constraints={"K3 > 0","K5 > 0","K6 > 0","K6 < 1"}
//	if (pixels)
//		CurveFit/Q/N Gauss2D image[V_left,V_right][V_bottom,V_top]/C=T_Constraints 
//	else
//		CurveFit/Q/N Gauss2D image(V_left,V_right)(V_bottom,V_top)/C=T_Constraints 
//	endif
//	KillWaves/Z T_Constraints
////	CurveFit/Q/N Gauss2D image[V_left,V_right][V_bottom,V_top]		// constraints on K6 are not needed for Igor 5
//// however, if constraints are not used, then K6 does not remain in [0,1]
////	if (V_FitError && WhichListItem("reFit_GaussianPkList",GetRTStackInfo(0))<0)
////	if (V_FitError && WhichListItem("reFit_GaussianPkList",GetRTStackInfo(0))<0 && WhichListItem("NewFitPeaks",GetRTStackInfo(0))<0)
//	if (V_FitError && WhichListItem("reFit_GaussianPkList",GetRTStackInfo(0))<0 && WhichListItem("NewFitPeaks",GetRTStackInfo(0))<0 && WhichListItem("FitPeaksStepWise",GetRTStackInfo(0))<0)
//		String errString = ""
//		if (V_FitError)
//			errString += "V_FitError = '"
//			errString += SelectString(V_FitError-2,"Singular Matrix","Out of memory","Function return NaN or INF")
//			errString += "'     "
//		endif
//		if (V_FitQuitReason)
//			errString += "V_FitQuitReason = '"
//			errString += SelectString(V_FitQuitReason-2,"iteration limit was reached","user stopped the fit","limit of passes without decreasing chi-square")
//			errString += "'"
//		endif
//		if (!Q)
//			printf "%s\r",errString
//		endif
//	endif
//	if (!V_FitError)
//		Wave W_coef = W_coef
//		Variable zero = -1e-15, one=1+1e-12					// almost zero, almost one
//		V_FitError += (W_coef[3]<zero || W_coef[5]<zero || W_coef[6]<zero || W_coef[6]>one) ? 32 : 0			//  constraints failed, set bit 5
//	endif
//	return V_FitError
//End


Function Add_COM_ToPeakList()
	GetMarquee/Z
	if (V_flag==0)
		DoAlert 0, "no Marquee selected, do nothing"
		return 1
	endif

	Wave image = ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";")))
	Wave peakList = $StringByKey("FullPeakList",GetUserData("","","FitPeaks"),"=")
#if exists("FindPeakListForImage")
	if (!WaveExists(peakList))						// not in user data, search for one
		Wave peakList = FindPeakListForImage(image)
	endif
#endif
#if exists("MakeEmptyPeakListForImage")
	if (!WaveExists(peakList))						// There is no peakList, make a new one
		Wave peakList = MakeEmptyPeakListForImage(image,ask=1,keyVals="FittedPeakShape=COM")
	endif
#endif
	if (!WaveExists(peakList) || !WaveExists(image))
		return 1
	endif
	if (!WaveInClass(peakList,"FittedPeakList*;"))
		return 1
	endif


	STRUCT PeakParams2dStruct pk
	GetCenterOfMass2d(pk,NaN,NaN,NaN,NaN,0)
	Variable px=pk.x0, py=pk.y0, net = numtype(pk.net) ? 1 : pk.net
	if (numtype(px+py))
		return 1
	endif

	// check for a duplicate within 1.5 pixel
	Variable N=DimSize(peakList,0), ipeak
	for (ipeak=0;ipeak<N;ipeak+=1)
		if (abs(peakList[ipeak][0]-px)<1.5 && abs(peakList[ipeak][1]-py)<1.5)
			break
		endif
	endfor
	String str
	if (ipeak>=N)
		sprintf str,"Add peak at (%.1f, %.1f)  Integral=%g",px,py,net
	else
		sprintf str,"REPLACE peak %d with:\rpeak at (%.1f, %.1f)  Integral=%g",ipeak,px,py,net
	endif
	Cursor/I/H=1 J $NameOfWave(image)  px,py
	DoAlert 1,str
	Cursor/K J
	if (V_flag!=1)
		return 1
	endif

	if (ipeak>=N)
		Redimension/N=(ipeak+1,-1) peakList
	endif
	peakList[ipeak][] = NaN;
	peakList[ipeak][0] = px;				peakList[ipeak][1] = py
	peakList[ipeak][10] = net

	if (NumberByKey("boxesOn"GetUserData("","","boxes"),"="))
		DrawMarker(px,py,20,20,"BoxWithTicks")
	endif
	return 0
End

// this gets the center of mass of a region on an image.  The region is selected with a marquee in an image plot
Function ShowCenterOfMass2d(printON)
	Variable printON					// enables print out of result
	STRUCT PeakParams2dStruct pk
	GetCenterOfMass2d(pk,NaN,NaN,NaN,NaN,printON)
End
//
// this gets the center of mass of a region on an image.  The region is selected with a marquee in an image plot
Function GetCenterOfMass2d(pk,V_left,V_right,V_bottom,V_top,printON)
	STRUCT PeakParams2dStruct &pk
	Variable V_left, V_right, V_bottom, V_top
	Variable printON					// enables print out of result

	initPeakParamsStruct(pk)
	pk.func = "COM"
	if (WinType("")!=1)
		Abort "Top window must be a graph, with region selected"
	endif
	if (!(V_left<V_right && V_bottom<V_top))
		String xaxis=StringByKey("XAXIS",ImageInfo("","",0)), yaxis=StringByKey("YAXIS",ImageInfo("","",0))
		GetMarquee/Z $yaxis, $xaxis
		if (V_flag==0)
			DoAlert 0, "No region selected on the graph"
			return 1
		endif
	endif

	// Find the wave that is the image (assume that it is the only 2-d wave on the graph
	String ImageName = StringFromList(0,WaveList("*",";","DIMS:2,WIN:"))	// name of wave with image
	if (strlen(ImageName)<1)
		DoAlert 0, "Unable to find image on top graph"
		return 1
	endif
	Wave image = $ImageName
	pk.units = WaveUnits(image,0)

	// find integer values for V_left,V_right,V_bottom,V_top
	Variable x0=DimOffset(image,0), dx=DimDelta(image,0), y0=DimOffset(image,1), dy=DimDelta(image,1)
	Variable i0,i1,j0,j1
	i0 = (V_left-x0)/dx
	i1 = (V_right-x0)/dx
	j0 = (V_bottom-y0)/dy
	j1 = (V_top-y0)/dy
	Variable swap
	if (j1<j0)						// order vertical indicies
		swap = j1
		j1 = j0
		j0 = swap
	endif
	if (i1<i0)						// order horizontal indicies
		swap = i1
		i1 = i0
		i0 = swap
	endif
	i0 = limit(floor(i0),0,DimSize(image,0)-1)	// these are pixels, so they need to be integers
	i1= limit(ceil(i1),0,DimSize(image,0)-1)
	j0 = limit(floor(j0),0,DimSize(image,1)-1)
	j1 = limit(ceil(j1),0,DimSize(image,1)-1)
	V_left = x0 + i0*dx				// reset to actual region used
	V_right = x0 + i1*dx
	V_bottom = y0 + j0*dy	
	V_top = y0 + j1*dy
	// printf "i=[%g, %g],  j=[%g, %g],    X=[%g, %g],  Y=[%g, %g]\r",i0,i1,j0,j1,V_left, V_right,V_bottom, V_top

	Variable nx=i1-i0+1, ny=j1-j0+1
	Make/N=(nx,ny)/FREE/D roi,xy
	SetScale/I x V_left,V_right,"", roi, xy	// values along X-direction
	SetScale/I y V_bottom,V_top,"", roi, xy	// values along Y-direction

	// find the background,  use average of the edge pixels
	roi = image[p+i0][q+j0]*(p==0 || p==(nx-1) || q==0 || q==(ny-1))	// 1 on edges, 0 inside
	Variable total = 2*(nx+ny-2)				// number of points in along edges
	Variable bkg = sum(roi)/total				// average background / pixel
	// print "bkg = ",bkg,"   total = ",total,"   ",nx,ny
	pk.bkg = bkg

	roi = image[p+i0][q+j0]
	ImageStats/M=1 roi
	total = V_avg*V_npnts
	pk.gross = total
	pk.net = (V_avg-bkg)*V_npnts

	xy = x
	MatrixOP/FREE/O val = sum(roi*xy)
	pk.x0 = val[0]/total

	xy = y
	MatrixOP/FREE/O val = sum(roi*xy)
	pk.y0 = val[0]/total

//	MatrixOP/FREE/O oneSum = sumRows(roi)
//	SetScale/I x V_left,V_right,"",  oneSum		// values along X-direction
//	oneSum *= x									// calculate the first moment
//	pk.x0 = sum(oneSum)/total
//
//	MatrixOP/FREE/O oneSum = sumCols(roi)^t
//	SetScale/I x V_bottom,V_top,"", oneSum	// values along Y-direction
//	oneSum *= x									// calculate the first moment
//	pk.y0 = sum(oneSum)/total

	if (printON)
		printf "Center of Mass of  '%s'[%g,%g][%g,%g] = (%g,%g)\r",ImageName,V_left, V_right, V_bottom, V_top, pk.x0,pk.y0
		if (numtype(pk.x0+pk.y0)==0)
			DoAlert 1, "Draw a cross at result"
			if (V_Flag==1)
				DrawMarker(pk.x0,pk.y0,V_right-V_left,V_top-V_bottom,"Cross")
			endif
		endif
	endif
	return 0
End



Structure PeakParams2dStruct
	double x0,y0, dx0,dy0						// center of peak, & errors
	double fwx, fwy, dfwx, dfwy					// fwhm of peak, & errors
	double xcorr, dxcorr							// x-corr, & error
	double amp, damp								// amplitude of peak, & error
	double bkg, dbkg								// background, & error
	double	gross, net, dgross, dnet				// gross and net integrals, with errors
	char func[100]								// string with functional from of fit, Gaussian, Lorentzian, COM, simple
	char units[100]								// units
EndStructure
//
Function initPeakParamsStruct(pk)
	STRUCT PeakParams2dStruct &pk
	pk.x0 = NaN;		pk.dx0 = NaN			// center of peak, & errors
	pk.y0 = NaN;		pk.dy0 = NaN
	pk.fwx = NaN;		pk.dfwx = NaN			// fwhm of peak, & errors
	pk.fwy = NaN;		pk.dfwy = NaN
	pk.xcorr = NaN;	pk.dxcorr = NaN		// x-corr, & error
	pk.amp = NaN;		pk.damp = NaN		// amplitude of peak, & error
	pk.bkg = NaN;		pk.dbkg = NaN			// background, & error
	pk.gross = NaN;	pk.dgross = NaN		// gross and net integrals, with errors
	pk.net = NaN;		pk.dnet = NaN
	pk.func="";		pk.units=""
End
//
Function/T printPeakParamsStruct(pk)
	STRUCT PeakParams2dStruct &pk

	if (numtype(pk.x0 + pk.y0))
		return "Invalid Peak"
	endif
	String str, out=pk.func
	out += SelectString(strlen(out),""," ")+"peak: "
	out += "("+valueAndError2Str(pk.x0,pk.dx0)+", "+valueAndError2Str(pk.y0,pk.dy0)+")"
	out += SelectString(strlen(pk.units),""," "+pk.units)
	if (numtype(pk.fwx + pk.fwy)==0)
		out += ",   FWHM = ("+valueAndError2Str(pk.fwx,pk.dfwx)+", "+valueAndError2Str(pk.fwy,pk.dfwy)+")"
		str  = valueAndError2Str(pk.xcorr,pk.dxcorr)
		out += SelectString(strlen(str),"",", Xcorr = "+str)
	endif
	out += SelectString(numtype(pk.amp),",  Amp = "+valueAndError2Str(pk.amp,pk.damp),"")
	out += SelectString(numtype(pk.bkg),",  Bkg = "+valueAndError2Str(pk.bkg,pk.dbkg),"")
	out += SelectString(numtype(pk.net),",  Net = "+valueAndError2Str(pk.net,pk.dnet),"")
	out += SelectString(numtype(pk.gross),",  Gross = "+valueAndError2Str(pk.gross,pk.dgross),"")
	return out
End
//
Static Function/T valueAndError2Str(val,err)
	Variable val,err
	String str
	if (numtype(val + err)==0)
		sprintf str "%g±%.2g",val,err
	elseif (numtype(val)==0)
		sprintf str "%g",val
	else
		str = ""
	endif
	return str
End
//
//Function test_PeakParams2d()
//	STRUCT PeakParams2dStruct pk
//
//	initPeakParamsStruct(pk)
//	print printPeakParamsStruct(pk)
//
//	pk.func = "Gaussian";		pk.units = "mm"
//	pk.x0=30;					pk.dx0=.023
//	pk.y0=50;					pk.dy0=0.5
//	pk.fwx=2.1234;			pk.dfwx=.015
//	pk.fwy=2.5234;			pk.dfwy=0.0234567
//	pk.xcorr=-0.7776;		pk.dxcorr=0.001
//	pk.amp = 1234567;		pk.damp=1234567/10
//	pk.bkg = 2345;			pk.dbkg=23459987/10
//	pk.gross= 1e9;			pk.net= 2e8
//	print printPeakParamsStruct(pk)
//
//	pk.func = "";				pk.units = ""
//	print printPeakParamsStruct(pk)
//
//	pk.func = "Gaussian";		pk.units = "mm"
//	pk.amp = NaN
//	print printPeakParamsStruct(pk)
//
//	pk.amp = 1234567;		pk.damp=1234567/10
//	print printPeakParamsStruct(pk)
//End

// ================================== End of Peak Evaluation ===================================
// =============================================================================================



// =============================================================================================
// ======================================== DEPRECATED =========================================

// This id DEPRECATED    DEPRECATED    DEPRECATED    DEPRECATED    DEPRECATED    DEPRECATED    
// Use  GenericWaveNoteInfo() directly
Function/T ShowImageInfo(image,key)
	Wave image
	String key
	return GenericWaveNoteInfo(image,key,class="speImage*;rawImage*",options="DIMS:2",type="image")
End

// ===================================== End of DEPRECATED =====================================
// =============================================================================================



// =============================================================================================
// ====================================== Initialization =======================================

Function initImageDisplayScaling()
	NewDataFolder/O root:Packages
	NewDataFolder/O root:Packages:imageDisplay
	String/G root:Packages:imageDisplay:imageExtension
	SVAR imageExtension=root:Packages:imageDisplay:imageExtension
	imageExtension = SelectString(strlen(imageExtension),".h5",imageExtension)
	setImageFileFilters()
End

// =================================== End of Initialization ===================================
// =============================================================================================
