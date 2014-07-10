#pragma rtGlobals=1		// Use modern global access method.
#pragma version = 1.57
#pragma IgorVersion = 6.0
#pragma ModuleName=depthResolve
//#include "microGeometry", version>=2.48
#include "ImageDisplayScaling", version>=1.81
#if (NumVarOrDefault("root:Packages:MICRO_GEOMETRY_VERSION",0)&2)
#include "tiff"
#else
#if (Exists("HDF5OpenFile")==4)
#include "HDF5images", version>=0.212
#include "FlyImageExtract"
#endif
#if (NumVarOrDefault("root:Packages:MICRO_GEOMETRY_VERSION",0)&4)
#include "WinView", version>=1.98
#else
#endif


// units for the metadata that usually accompanies an image, use this in case HDF5_MetaUnits and WinView_MetaUnits do not exist
Static Strconstant localKeyUnits = "X1:µm;Y1:µm;Z1:µm;H1:µm;F1:µm;X2:µm;Y2:µm;Z2:µm;H2:µm;F2:µm;depth:µm;keV:keV;exposure:s;ringCurrent:mA;"


Menu LaueGoMainMenuName
	SubMenu "DepthResolved Query"
		"Plot Average Intensity vs Depth", IntensityVsDepth("","","","",0,NaN,0,NaN)
		"Make Movie of one type of Scan ",MovieOfOneScan()
		help={"movie of 1 wire scan, some sequential images, or first white image of each wire scan"}
		"Load recontruction 'summary' file", LoadIntegralFile("")
		"  Re-Display intensity from summary",DisplayReconstructionIntegral($"",$"")
		"Sum a ROI in many Images",SumInManyROIs("","","",$"")
		"Make Waves With Metadata from Image Range",ImageMetaData2Waves("","","","")	
		"-"
		MenuItemIfWaveClassExists("Re-Scale Movie Pixel Intensities","PixelIntensities","Dims:2"), RescaleMoviePixelIntensities($"")
		MenuItemIfWaveClassExists("Graph Movie Pixel Intensities","PixelIntensities*","Dims:2"), GraphMoviePixelIntensities($"")
	End
End

Menu "GraphMarquee",dynamic
	"-"
	MarqueeDepthImageMenuItem("Show Average Intensity vs Depth"), MarqueeGetIntensityVsDepth()
	help={"sum up intensity withing marquee for each image of a depth reconstruction"}
End




//Function/T ReadGenericHeader(fName)
//	String fName					// fully qualified name of file to open (will not prompt)
//	String str
//#if (Exists("LoadHDF5imageFile")==6)
//	str = ReadHDF5header(fName)
////	DoAlert 0,"LoadHDF5imageFile Header is not yet done"
////	return ""
//#elif (Exists("LoadWinViewFile")==6)
//	str = WinViewReadHeader(fName)
//#endif
//	return str
//End



//Static Function/S LoadGenericImageFile(fileName,[dark])
//	String fileName
//	Wave dark													// optional darkImage to subtract
//
//	String str=""
//#if (Exists("LoadHDF5imageFile")==6)
//	if (ParamIsDefault(dark))
//		str = LoadHDF5imageFile(fileName)
//	elseif(!WaveExists(dark))
//		str = LoadHDF5imageFile(fileName)
//	else
//		str = LoadHDF5imageFile(fileName,dark=darkImageOrange)
//	endif
//#elif (Exists("LoadWinViewFile")==6)
//	str = LoadWinViewFile("")
//#endif
//	return str
//End


Function/S MarqueeDepthImageMenuItem(item)
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
	if (!WaveExists(image))					// has to be an image
		return "("+item
	endif
	String class = StringByKey("waveClass", note(image),"=")
	if (!stringmatch(class,"spe*"))				// only works for spe images
		item = "("+item
		// item = "("
	endif

	if (!stringmatch(class,"spe*"))				// only works for spe images
		item = "("+item
	endif

	String namePart = ParseFilePath(3,StringByKey("imageFileName",note(image),"="), ":", 0, 1)
	Variable OK = 0
	OK = OK || numtype(NumberByKey("depth",note(image),"="))==0
	OK = OK || numtype(NumberByKey("depthSi",note(image),"="))==0	// also check the old name
	OK = OK || stringmatch(namePart,"*depth*")
	OK = OK || isdigit(namePart[strlen(namePart)-1])
	if (!OK)
		item = "("+item						// only want depth sorted images, or probably a wire scan
	endif
//	if (numtype(NumberByKey("depthSi",note(image),"=")) && !stringmatch(StringByKey("imageFileName",note(image),"="),"*depth*"))
//		item = "("+item						// only want depth sorted images
//	endif
	return item
End
//Function/S MarqueeDepthImageMenuItem(item)
//	String item
//	if (strlen(WinList("*","","WIN:1"))<1)
//		return "("+item
//	endif
//	Variable V_flag
//	GetMarquee/Z								// These menu items requre a Marquee present
//	if (V_flag==0)
//		return "("+item
//	endif
//	Wave image = ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";")))
//	if (!WaveExists(image))					// has to be an image
//		return "("+item
//	endif
//	String class = StringByKey("waveClass", note(image),"=")
//	if (!stringmatch(class,"spe*"))				// only works for spe images
//		item = "("+item
//		// item = "("
//	endif
//
//	if (!stringmatch(class,"spe*"))				// only works for spe images
//		item = "("+item
//	endif
//	if (numtype(NumberByKey("depthSi",note(image),"=")) && !stringmatch(StringByKey("imageFileName",note(image),"="),"*depth*"))
//		item = "("+item						// only want depth sorted images
//	endif
//	return item
//End
//
Function MarqueeGetIntensityVsDepth()
	if (WinType("")!=1)
		Abort "Top window must be a graph, with region selected"
	endif
	Variable V_left, V_right, V_bottom, V_top
	GetMarquee/Z left, bottom
	if (V_flag==0)
		DoAlert 0, "No Marquee selected on the graph, do nothing"
		return 1
	endif

	//	Find the wave that is the image (assume that it is the 1st 2-d wave on the graph
	String ImageName =StringFromList(0,ImageNameList("",";"))	// name of wave with image
	if (strlen(ImageName)<1)
		DoAlert 0, "Unable to find image on top graph"
		return 1
	endif
	Wave image = ImageNameToWaveRef("",ImageName)
	String wnote = note(image)

	V_left = round((V_left-DimOffset(image,0))/DimDelta(image,0))		// convert from scaled to index
	V_right = round((V_right-DimOffset(image,0))/DimDelta(image,0))
	V_top = round((V_top-DimOffset(image,1))/DimDelta(image,1))
	V_bottom = round((V_bottom-DimOffset(image,1))/DimDelta(image,1))
	Variable ilo=min(V_left,V_right), ihi=max(V_left,V_right)				//
	Variable jlo=min(V_top,V_bottom), jhi=max(V_top,V_bottom)

	String filePrefix=StringByKey("imageFileName",wnote,"=")			// first part of file name, e.g. "WH_"
	String pathStr=StringByKey("imageFilePath",wnote,"=")
	String pathName="", list = PathList("*",";","")
	Variable i,N=ItemsInList(list)
	for (i=0;i<N;i+=1)
		PathInfo $StringFromList(i,list)
		if (stringmatch(S_path,pathStr))
			pathName = StringFromList(i,list)
			break
		endif
	endfor
	if (strlen(pathName)<1)
		pathName = "reconPath"
		NewPath/M=""/O $pathName, pathStr
	endif
	filePrefix = ParseFilePath(3,filePrefix,":",1,0)						// remove file extension
	for (i=strlen(filePrefix)-1; i>0 && isdigit(filePrefix[i,i]) ; i-=1)	// remove trailing digits
	endfor
	filePrefix = filePrefix[0,i]
	//	printf " x = [%g, %g],  y = [%g, %g]\r",ilo,ihi,jlo,jhi
	//	printf "'%s'       '%s'\r",pathStr,filePrefix
	IntensityVsDepth(pathName,filePrefix,"","",ilo,ihi,jlo,jhi)
End


// =========================================================================
// =========================================================================
//	Start of Detail set panel

Function/T FillDetailParametersPanel(strStruct,hostWin,left,top)
	String strStruct									// optional passed value of xtal structure, this is used if passed
	String hostWin										// name of home window
	Variable left, top									// offsets from the left and top

	SetWindow kwTopWin,userdata(DetailPanelName)=hostWin+"#DetailPanel"
	NewPanel/K=1/W=(left,top,left+221,top+445+30)/HOST=$hostWin
	ModifyPanel frameStyle=0, frameInset=0
	RenameWindow #,DetailPanel

	Button buttonMakeMovie,pos={29,5},size={160,20},proc=depthResolve#DetailButtonProc,title="Make a Movie"
	Button buttonMakeMovie,help={"Make a Movie of one type of scan, lots of options"}

	Button buttonAvgIntvsDepth,pos={29,40},size={160,20},proc=depthResolve#DetailButtonProc,title="Plot Intensity vs Depth"
	Button buttonAvgIntvsDepth,help={"Plot Average Intensity vs Depth for a Reconstruction"}

	Button buttonLoadSummary,pos={29,75},size={160,20},proc=depthResolve#DetailButtonProc,title="Load Recon Summary"
	Button buttonLoadSummary,help={"Load a reconstruction summary file and plot the intensity"}

	Button buttonRePlotSummary,pos={49,100},size={140,20},proc=depthResolve#DetailButtonProc,title="Re-Plot Summary"
	Button buttonRePlotSummary,help={"Re-Plot intensity vs depth loaded with the button above"}

	Button buttonSumImages,pos={29,135},size={160,20},proc=depthResolve#DetailButtonProc,title="Sum Many Images"
	Button buttonSumImages,help={"Sum a set of images with the option of using an ROI"}

	Button buttonMetaDataImages,pos={29,170},size={190,20},proc=depthResolve#DetailButtonProc,title="Metadata from Image Range"
	Button buttonMetaDataImages,help={"Extract header information from a range of images"}

	//rxadd
	Button buttonFlyExtract,pos={29,205},size={190,20},proc=depthResolve#DetailButtonProc,title="Batch Extract Flyscan Images"
	Button buttonFlyExtract,help={"Extract single images from flyscans and export them as independent HDF files"}
	//end rxadd

	EnableDisableDetailControls(hostWin+"#DetailPanel")
	return "#DetailPanel"
End
//
Static Function EnableDisableDetailControls(win)				// here to enable/disable
	String win												// window (or sub window) to use
	Variable d

	Button buttonMakeMovie,win=$win,disable=0
	Button buttonAvgIntvsDepth,win=$win,disable=0
	Button buttonLoadSummary,win=$win,disable=0
	Button buttonSumImages,win=$win,disable=0
	Button buttonMetaDataImages,win=$win,disable=0
	Button buttonFlyExtract,win=$win,disable=0 		//rxadd
	d = strlen(WaveListClass("intensityVsDepth,intensityRecon","*","DIMS:1"))<1 ? 2 : 0
	Button buttonRePlotSummary,win=$win,disable= d
End
//
Static Function DetailButtonProc(B_Struct) : ButtonControl
	STRUCT WMButtonAction &B_Struct
	if (B_Struct.eventCode != 2)
		return 0
	endif
	String ctrlName=B_Struct.ctrlName

	STRUCT microGeometry g
	FillGeometryStructDefault(g)

	if (stringmatch(ctrlName,"buttonMakeMovie"))
		MovieOfOneScan()
	elseif (stringmatch(ctrlName,"buttonAvgIntvsDepth"))
		IntensityVsDepth("","","","",0,NaN,0,NaN)
	elseif (stringmatch(ctrlName,"buttonLoadSummary"))
		LoadIntegralFile("")
	elseif (stringmatch(ctrlName,"buttonRePlotSummary"))
		DisplayReconstructionIntegral($"",$"")
	elseif (stringmatch(ctrlName,"buttonSumImages"))
		SumInManyROIs("","","",$"")
	elseif (stringmatch(ctrlName,"buttonMetaDataImages"))
		ImageMetaData2Waves("","","","")
	//rxadd
	elseif (stringmatch(ctrlName,"buttonFlyExtract"))
		flyImageExtract()
	//end rxadd
	endif
	EnableDisableDetailControls(GetUserData("microPanel","","DetailPanelName"))
End

//	End of Detail set panel
// =========================================================================
// =========================================================================


//  =================================== Start of Movies ===================================
//  ==================================================================================

Function MovieOfOneScan()		// make the movie frame and fill it
	String gName=FindMovieGraph()
	if (strlen(gName))
		DoWindow/F $gName
		FillMovieOfOneScan("","","",NaN,NaN,"")
	else
		gName = MakeMovieWindow()
		if (strlen(gName)<1)
			return 1
		endif
		DoUpdate
		DoAlert 0,"Adjust the color scale and other setting on this image, then start the movie by pushing the 'Continue' button"
		Button ContinueButton,pos={1,2},size={80,30},proc=MakeMovieContinueButtonProc,title="Continue"
		Button ContinueButton,fColor=(16386,65535,16385)
	endif
	return 0
End
//
Function MakeMovieContinueButtonProc(ctrlName) : ButtonControl
	String ctrlName
	KillControl ContinueButton
	DoUpdate
	FillMovieOfOneScan("","","",NaN,NaN,"")
End

Function/T MakeMovieWindow()
	// pick an image to use for the movie
	String wList = reverseList(WaveListClass("speImage*","*","DIMS:2"))
	wList = wList+"read a new image;"
	String wName = "imageOnMovie"
	if (ItemsInList(wList)>1)
		Prompt wName, "sample image to use for the movie", popup, wList
		DoPrompt "choose an image", wName
		if(V_flag)
			return ""
		endif
	else
		wName = StringFromList(0,wList)
	endif
	if (stringmatch(wName,"read a new image"))		// read in new image to use for imageOnMovie
		wName = LoadGenericImageFile("")
		if (strlen(wName)<1)
			return ""
		endif
		Wave image = $wName
		Duplicate/O image imageOnMovie
		KillWaves/Z image
	elseif(!stringmatch(wName,"imageOnMovie"))		// use existing image to create imageOnMovie
		Wave image = $wName
		if (WaveExists(image))
			Duplicate/O image imageOnMovie
		endif
	else													// use existsing imageOnMovie
		Wave imageOnMovie=imageOnMovie
	endif
	if (!WaveExists(imageOnMovie))
		DoAlert 0, "cannot find image to use for setting up movie"
		return ""
	endif

	// finally have the image to use for the movie, make the graph
	Display /W=(448,130,1180,788)
	AppendImage :imageOnMovie

	WaveStats/Q imageOnMovie
	if (V_avg>100)
		ModifyImage imageOnMovie ctab= {*,*,Terrain,1}
	else
		ModifyImage imageOnMovie ctab= {0,100,Terrain,1}
	endif
	ModifyGraph grid=1, tick=2, mirror=1, minor=1, lowTrip=0.001, standoff=0
	ModifyGraph axOffset(left)=-1.1,axOffset(bottom)=-1
	ModifyGraph gridStyle=1, gridRGB=(45000,45000,65535)
	Label left "y  (\\U)"
	Label bottom "x  (\\U)"
	SetAxis/A/R left
//	Cursor/P/I A imageOnMovie 88,551
//	ShowInfo
//	TextBox/C/N=textParameters ""
//	TextBox/C/N=textFileName/F=0/B=1/A=LT ""
	TextBox/C/N=textParameters/B=1/A=RB ""
	TextBox/C/N=textFileName/F=0/B=1/A=LB ""
	SetWindow kwTopWin,userdata(windowType)=  "OneMovieWindow"
	DoUpdate
	SetAspectToSquarePixels("")	
//	DoAlert 0,"Adjust color scale on this image and adjust other items before starting movie"
	return WinName(0,1,1)
End
//
Function FillMovieOfOneScan(pathName,filePrefix,range,surface,absorpLength,type,[makeSum,movieName,flatten,skipZeros,moreFunc,NoMovie])
	String pathName					// name of path to use
	String filePrefix					// first part of file name, e.g. "WH_4174_"
	String range						// range of files to use
	Variable surface					// depth of surface relative to depthSi(µm), used to calculate boost
	Variable absorpLength			// absorption length (µm)
	String type							// type of movie, images of a wire scan, sequenial images, first image of each wire scan
	Variable makeSum
	String movieName					// name to use for movie (optional)
	Variable flatten					// flag, TRUE make flat movie for Windows
	Variable skipZeros				// flat, TRUE means skip empty images, default is TRUE
	String moreFunc					// name of FillMovieMore function
	Variable NoMovie					// suppress actual creation of Movie
	makeSum = ParamIsDefault(makeSum) ? 0 : (numtype(makeSum)==0 && makeSum != 0)
	movieName = SelectString(ParamIsDefault(movieName),movieName,"")
	flatten = ParamIsDefault(flatten) ? 0 : (numtype(flatten)==0 && flatten != 0)
	skipZeros = ParamIsDefault(skipZeros) ? 1 : (numtype(skipZeros)==0 && skipZeros != 0)
	moreFunc = SelectString(ParamIsDefault(moreFunc),moreFunc,"")
	NoMovie = ParamIsDefault(NoMovie) ? 0 : (numtype(NoMovie)==0 && NoMovie)
	String typeList="single wire scan;sequential images;first image of each wirescan;first white image of each flyScan;all of one flyScan"
	Variable useLog=0

	String gName = FindMovieGraph()						// name of suitable graph with imageOnMovie
	DoWindow/F $gName											// and bring it to the front

	Variable itype = WhichListItem(type, typeList)
	if (itype<0)												// no type passed, or invalid
		Prompt type, "type of movie to make", popup, typeList
		Prompt makeSum,"Also make a sum image",popup,"No Sum;Make Sum Too"
		Prompt flatten,"Flatten movie for windows",popup,"Regular;Flatten for Windows"
		Prompt skipZeros,"Flatten movie for windows",popup,"Show All Images;Skip Empty Images"
		Prompt moreFunc,"Name of Function to run at each step",popup,"-none-;"+FunctionList("FillMovieMore*",";","NPARAMS:2,KIND:2")
		Prompt NoMovie,"Do NOT make the movie",popup,"Make the Movie;NO Movie"
		flatten += 1
		skipZeros += 1
		NoMovie += 1
		DoPrompt "movie type",type,makeSum,flatten,skipZeros,moreFunc,NoMovie
		if (V_flag)
			return 1
		endif
		makeSum = makeSum==2
		flatten = flatten==2
		skipZeros = skipZeros==2
		NoMovie = NoMovie == 2
		moreFunc = SelectString(StringMatch(moreFunc,"-none-"),moreFunc,"")
	endif
	itype = WhichListItem(type, typeList)

	String imageExtension=StrVarOrDefault("root:Packages:imageDisplay:imageExtension",".h5")
	GetFileFolderInfo/P=$pathName/Q/Z filePrefix			// this may be a single flyscan file, itype==4
	if (!V_isFile)
		GetFileFolderInfo/P=$pathName/Q/Z filePrefix+num2istr(str2num(range))+imageExtension
	endif
	if (!V_isFile)
		GetFileFolderInfo/P=$pathName/Q filePrefix+num2istr(str2num(range))+imageExtension
	endif
	if (!V_isFile)
		return 1
	endif
	pathName = SelectString(strlen(pathName),"reconPath",pathName)	// use reconPath as the default path name, if none passed
	String fname=S_Path
	String path = ParseFilePath(1,S_Path,":",1,0)
	NewPath/O/Q $pathName, ParseFilePath(1,S_Path,":",1,0)
	filePrefix = ParseFilePath(3,S_Path,":",1,0)			// remove fle extension
	Variable i
	for (i=strlen(filePrefix)-1; i>0 && isdigit(filePrefix[i,i]) ; i-=1)// remove trailing digits
	endfor
	filePrefix = filePrefix[0,i]
	String fileRoot = path+filePrefix

	String header=ReadGenericHeader(fname)
	Variable epoch=DateTime, Nslices=NumberByKey("Nslices",header,"=")
	Nslices = Nslices>=2 ? Nslices : 0
	if (strlen(range)<1)											// find the range here
		if (itype==4)													// one whole flyscan, itype==4
			if (Nslices<2)
				print "asked for a movie of a fly scan, but only one image!"
				return 1
			endif
			range = "0-"+num2istr(Nslices-1)
			Prompt range,"range in fly scan"
			DoPrompt "Range of Fly Scan",range
			if (V_flag)
				return 1
			endif
			if (lastInRange(range)>=Nslices || str2num(range)<0)
				DoAlert 0, "You CANNOT ask for more images than there are in fly scan"
				return 1
			endif
		elseif (itype <=1 && itype>=0 || itype==3)			// for itype == 0 or 1 or 3
			range = get_FilesIndexRange(pathName,filePrefix)
		elseif (itype ==2 )											// first image of each wirescan, itype==2
			range = get_FirstWhiteFilesIndexRange(pathName,filePrefix)
		else
			DoAlert 0, "in FillMovieOfOneScan() itype = "+num2str(itype)
			print "ERROR -- in FillMovieOfOneScan() itype = "+num2str(itype)
			return 1
		endif
	endif
	Variable totalTime = DateTime-epoch						// time to figure out which scan numbers to use
	Variable Nimages = ItemsInRange(range)					// number of images to put in movie
	if (Nimages<1)
		return 1
	endif
	if ((numtype(surface)==2 || !(absorpLength>0)) && itype==0)
		absorpLength = !(absorpLength>0) ? 15 : absorpLength
		Prompt surface, "depth of surface relative to depthSi(µm), (use NaN to ignore)"
		Prompt absorpLength, "absorption length in material (µm), (use NaN to ignore)"
		DoPrompt "absorption correction",surface,absorpLength
		if (V_flag)
			return 1
		endif
	endif
	if (exists(moreFunc)==6)										// hook for a user supplied function to do some more
		FUNCREF MoreInFillMovieProto func=$moreFunc
		if (exists(moreFunc+"_Init")==6)						// an optional init file associated with moreFunc
			FUNCREF MoreInFillMovieInitProto funcInit=$(moreFunc+"_Init")
		endif
	else
		FUNCREF MoreInFillMovieProto func=MoreInFillMovieProto
		FUNCREF MoreInFillMovieInitProto funcInit=MoreInFillMovieInitProto
	endif

	String ff = SelectString(itype==4,filePrefix,ParseFilePath(0,fname,":",1,0))
	printf "FillMovieOfOneScan(\"%s\",\"%s\",\"%s\",%g,%g,\"%s\"",pathName,ff,range[0,100],surface,absorpLength,type
	if (!ParamIsDefault(makeSum) || makeSum)
		printf ", makeSum=%g",makeSum
	endif
	if (!ParamIsDefault(movieName) && strlen(movieName))
		printf ", movieName=\"%s\"",movieName
	endif
	if (!ParamIsDefault(flatten) || flatten)
		printf ", flatten=%g",flatten
	endif
	if (!ParamIsDefault(skipZeros) || !skipZeros)
		printf ", skipZeros=%g",skipZeros
	endif
	if (strlen(moreFunc))
		printf ", moreFunc=\"%s\"",moreFunc
	endif
	if (NoMovie)
		printf ", NoMovie=1"
	endif
	printf ")\r"
	printf "running in data folder  '%s',    %g images in range\r",GetDataFolder(1),Nimages
	String inputs=""
	inputs = ReplaceStringByKey("pathName",inputs,pathName,"=")
	inputs = ReplaceStringByKey("filePrefix",inputs,filePrefix,"=")
	inputs = ReplaceStringByKey("range",inputs,range,"=")
	inputs = ReplaceStringByKey("fileRoot",inputs,fileRoot,"=")
	inputs = ReplaceStringByKey("imageExtension",inputs,imageExtension,"=")
	inputs = ReplaceNumberByKey("surface",inputs,surface,"=")
	inputs = ReplaceNumberByKey("absorpLength",inputs,absorpLength,"=")
	inputs = ReplaceNumberByKey("itype",inputs,itype,"=")
	inputs = ReplaceStringByKey("movieName",inputs,movieName,"=")
	inputs = ReplaceNumberByKey("flatten",inputs,flatten,"=")
	inputs = ReplaceNumberByKey("skipZeros",inputs,skipZeros,"=")
	
	epoch = DateTime
	Variable boost														// boost factor from absorption
	Variable depthSi,integral
	Variable X1,H1,keV, X10,H10,keV0
	Wave imageOnMovie=imageOnMovie
	String wnote = note(imageOnMovie)
	Variable i0=0, i1=-1, j0=0, j1=-1
	// check if we are only using a sub-image for the movie
	if (DimSize(imageOnMovie,0)!=NumberByKey("xdim", wnote,"=") || DimSize(imageOnMovie,1)!=NumberByKey("ydim", wnote,"="))
		i0 = NumberByKey("startxRead", wnote,"=")
		j0 = NumberByKey("startyRead", wnote,"=")
		i0 = i0>0 ? i0 : 0
		j0 = j0>0 ? j0 : 0
		i1 = i0 + DimSize(imageOnMovie,0) - 1
		j1 = j0 + DimSize(imageOnMovie,1) - 1
		printf "using a sub-region of the images [%d, %d] [%d, %d]\r",i0,i1,j0,j1
	endif
	inputs = ReplaceNumberByKey("i0",inputs,i0,"=")
	inputs = ReplaceNumberByKey("j0",inputs,j0,"=")
	inputs = ReplaceNumberByKey("i1",inputs,i1,"=")
	inputs = ReplaceNumberByKey("j1",inputs,j1,"=")

	Make/N=(Nimages)/O MovieIntegral=NaN						// holds the integral for each frame of the movie
	Variable ii

	String str,textStr, imageName
	//	sprintf fname, "%s%d%s",fileRoot,str2num(range),imageExtension
	// Variable Nslices = NumberByKey("Nslices",ReadGenericHeader(fname),"=")
	//	Nslices = Nslices>=2 ? Nslices : 0
	inputs = ReplaceNumberByKey("Nslices",inputs,Nslices,"=")
	totalTime += DateTime-epoch
	String moreFuncInitStr = funcInit(inputs)
	if (StringMatch(moreFuncInitStr,"**Quit**"))
		return 1
	endif

	epoch = DateTime
	if (!NoMovie)
		String cmd="NewMovie/P=home"
		cmd += SelectString(flatten,"","/L")
		cmd += SelectString(strlen(movieName),""," as "+movieName)
		Execute cmd														// start the movie
	endif

	for (i=str2num(range),ii=0; !numtype(i); i=NextInRange(range,i),ii+=1)
		if (itype==4)
			imageName = ReadGenericROI(fname,i0,i1,j0,j1,extras="slice:"+num2istr(i))	// load one image from the file
		elseif (itype==3)
			sprintf fname, "%s%d%s",fileRoot,i,imageExtension
			imageName = ReadGenericROI(fname,i0,i1,j0,j1,extras="slice:0")	// only load the first image in file
		elseif (Nslices==0)
			sprintf fname, "%s%d%s",fileRoot,i,imageExtension
			imageName = ReadGenericROI(fname,i0,i1,j0,j1)	// returns full path name to loaded image wave
		else																// get all images from one multi image file
			sprintf fname, "%s%d%s",fileRoot,str2num(range),imageExtension
			imageName = ReadGenericROI(fname,i0,i1,j0,j1,extras=ReplaceNumberByKey("slice","",i))	// returns full path name to loaded image wave
		endif
		if (exists(imageName)!=1)
			continue
		endif
		Wave image = $imageName
		if (ii==0 && !(WaveType(image)&0x40))					// if the first frame is signed, ensure imageOnMovie is signed
			Duplicate/O image imageOnMovie
		endif
		if (ii==0 && makeSum)
			Duplicate/O imageOnMovie imageOnMovieSum
			Variable wtype = WaveType(imageOnMovieSum)
			if (wtype & 0x18)											// 8 or 16 bit integer, change to 32 bit
				wtype = (wtype & 0x40) | 0x20					// preserve the unsigned bit
				Redimension/Y=(wtype) imageOnMovieSum
			endif
			imageOnMovieSum = 0
		endif
		wnote = note(image)
	 	depthSi = NumberByKey("depth", wnote,"=")			// values read from image file
	 	depthSi = numtype(depthSi) ? NumberByKey("depthSi", wnote,"=") : depthSi	// try the old name
		X1 = NumberByKey("X1", wnote,"=")
		H1 = NumberByKey("H1", wnote,"=")
		keV = NumberByKey("keV", wnote,"=")
		if (ii==0)														// save first values for comparison
			X10=X1 ;  H10=H1 ;  keV0=KeV
		endif
		integral = sum(image)
		if (integral==0 && skipZeros)							// image is empty, skip it
			printf "skipping '%s',  slice = %g,   it is empty\r",imageName,i
			KillWaves/Z image
			continue
		endif
		MovieIntegral[ii] = integral
		boost = exp(max(depthSi-surface,0)/absorpLength)
		boost = numtype(boost) ? 1 : boost
		boost = itype==0 ? boost : 1
		imageOnMovie = image
		if (makeSum)
			imageOnMovieSum += imageOnMovie
		endif
		wnote = ReplaceNumberByKey("boost",wnote,boost,"=")
		wnote = ReplaceNumberByKey("movieIndex",wnote,ii,"=")		// movie frame index
		wnote = ReplaceNumberByKey("movieFileIndex",wnote,i,"=")	// index from range
		wnote = ReplaceNumberByKey("movieLength",wnote,Nimages,"=")
		Note/K imageOnMovie, wnote

		KillWaves/Z image												// done with image

		textStr = ""
		if (numtype(depthSi)==0)
			sprintf textStr, "%+.1f µm"+SelectString(useLog,"","\rLog")+"\r",depthSi
		endif
		if (numtype(X1+H1+keV)==0 && (itype>0 || ii==0))
			sprintf str, "X1 = %.2f\rH1 = %.2f\rkeV = %.4f\r",X1,H1,keV
			textStr += str
			str = ""
		endif
		if (integral==0)
			sprintf str, "\\Zr070(%d), º=0",i
		elseif (boost>1)
			sprintf str, "\\Zr070(%d), º=%.2fE6\rboost=%.1f",i,integral/1e6,boost
		else
			sprintf str, "\\Zr070(%d), º=%.2fE6",i,integral/1e6
		endif
		textStr += str
		TextBox/C/N=textParameters textStr
	 	str = StringByKey("imageFileName", wnote,"=")
		textStr = SelectString(strlen(str),"",str)
	 	str = StringByKey("imageFilePath", wnote,"=")
		if (strlen(str))
			textStr += SelectString(strlen(textStr),"","\r") + "\\Zr070" + str
		endif
		TextBox/C/N=textFileName textStr						// "filesname\r\\Zr070file path"
		func(imageOnMovie,moreFuncInitStr)
		imageOnMovie *= boost
		if (useLog)
			imageOnMovie = 5000*log(imageOnMovie)				// need 5000 because image is integer
		endif
		if (!NoMovie)
			DoUpdate
			AddMovieFrame
		endif
	endfor
	if (!NoMovie)
		CloseMovie
	endif
	totalTime += DateTime-epoch
	print "done, total execution time is  ",Secs2Time(totalTime,5,0)
	return 0
End

Function MoreInFillMovieProto(image,initStr)
	Wave image
	String initStr					// a generic mechanism to pass specialized info from MoreInFillMovieInitProto() to here
End
Function/T MoreInFillMovieInitProto(list)
	String list
	return ""
End
//Function MoreInFillMovieProto(image,[index])
//	Wave image
//	Variable index
//End
//		// This is an example:
//Function FillMovieMoreTest(image,initStr)
//	Wave image
//	String initStr			// this is not used
//	FitPeaksWithSeedFill(image,1.2,70,50,20,maxNu=15)
//	ReDrawBoxes("")
//	Wave FullPeakList=FullPeakList
//	if (DimSize(FullPeakList,0)<4)
//		ClearhklOffGraph("")
//		return 1
//	endif	
//
//	if (!IndexAndDisplay(FullPeakList,18,39,0.5, 0,1,1,30))
//		DisplayResultOfIndexing(FullPeakIndexed,0)
//	else
//		ClearhklOffGraph("")
//	endif
//End



Function FillMovieMorePixelsIntens(image,initStr)	// fills pixel intensities along a movie
	Wave image
	String initStr				// string with extra info from FillMovieMorePixelsIntens_Init

	Variable hw=NumVarOrDefault(":MoviePixelsIntensHW",0)
	hw = numtype(hw) || hw<0 ? 0 : round(hw)
	String wnote=note(image)
	Variable ii=NumberByKey("movieIndex",wnote,"="), iLen=NumberByKey("movieLength",wnote,"="), depth, Npixels, saveWires=0,X2,Y2,Z2
	Wave PeakList=$StrVarOrDefault(":MoviePixelsIntensPeakListName","")
	Wave MoviePixelIntensities=MoviePixelIntensities
	if (ii>0)											// NOT first frame in movie
		if (!WaveExists(PeakList) || !WaveExists(MoviePixelIntensities))
			return 1										// cannot find waves, not doing anything
		elseif (DimSize(PeakList,0)<1 || DimSize(PeakList,1)<2)
			return 1										// PeakList is invalid
		endif
	endif

	if (ii==0)											// first time through, ask for PeakList & make MoviePixelIntensities
		String PeakListName=NameOfWave(PeakList)
		String wlist=WaveList("*",";","DIMS:2,MAXCOLS:15,MINCOLS:2")
		Prompt PeakListName,"Wave with list of peaks to track",popup, wlist
		Prompt hw,"HW of area around each listed peak to save"
		if (ItemsInList(wlist)<1 && WaveExists(PeakList))
			DoPrompt "Record Peak Intensities", hw
		else
			DoPrompt "Record Peak Intensities", PeakListName,hw
		endif
		if (V_flag)
			DoAlert 0, "Will not be recording pixel intensitites in this movie"
			return 1
		endif
		hw = numtype(hw) || hw<0 ? 0 : round(hw)
		Wave PeakList = $PeakListName
		if (!WaveExists(PeakList))
			DoAlert 0, "Will not be recording pixel intensitites in this movie"
			return 1
		endif
		Variable/G :MoviePixelsIntensHW = hw		// must be saved in a global
		String/G :MoviePixelsIntensPeakListName = NameOfWave(PeakList)
		Npixels = DimSize(PeakList,0)
		Wave MoviePixelIntensities=MoviePixelIntensities

		if (!WaveExists(MoviePixelIntensities))
			Make/N=(iLen,Npixels)/O MoviePixelIntensities=NaN
		else
			Redimension/N=(iLen,Npixels) MoviePixelIntensities
		endif
		SetScale/P x 0,1,"", MoviePixelIntensities

		String pixelStr="",str						// need to save pixelStr in wave note for plotting
		Variable i, px,py
		for (i=0;i<Npixels;i+=1)
			px = PeakList[i][0]
			py = PeakList[i][1]
			sprintf str, "%g:%g,",round(px),round(py)
			pixelStr += str
		endfor

		depth = NumberByKey("depth",wnote,"=")
		if (numtype(depth)==0)
			wnote = ReplaceNumberByKey("firstDepth",wnote,depth,"=")
		else
			X2 = NumberByKey("X2",wnote,"=")	// so no depth but wire position available
			Y2 = NumberByKey("Y2",wnote,"=")
			Z2 = NumberByKey("Z2",wnote,"=")
			saveWires = (numtype(X2+Y2+Z2)==0)
		endif
		wnote = ReplaceStringByKey("PixelPositions",wnote,pixelStr,"=")
		wnote = ReplaceNumberByKey("MoviePixelsIntensHW",wnote,hw,"=")
		wnote = ReplaceNumberByKey("saveWires",wnote,saveWires,"=")
		wnote = ReplaceStringByKey("waveClass",wnote,"PixelIntensities","=")
		Note/K MoviePixelIntensities, wnote
		printf "Creating wave \"%s\" containing intensity around pixels (±%g) listed in \"%s\"\r",NameOfWave(MoviePixelIntensities),hw,NameOfWave(PeakList)
	else
		String pnote=note(MoviePixelIntensities)
		Npixels = DimSize(PeakList,0)
		Variable x0 = NumberByKey("firstDepth",pnote,"=")
		depth = NumberByKey("depth",wnote,"=")
		if (numtype(x0+depth))
			SetScale/P x 0,1,"", MoviePixelIntensities				// no depth scaling
		else
			SetScale/I x x0,depth,"µm", MoviePixelIntensities	// can set depth scaling
		endif
		if (NumberByKey("saveWires",pnote,"="))
			str = StringByKey("X2",pnote,"=")+","+StringByKey("X2",wnote,"=")
			pnote = ReplaceStringByKey("X2",pnote,str,"=")
			str = StringByKey("Y2",pnote,"=")+","+StringByKey("Y2",wnote,"=")
			pnote = ReplaceStringByKey("Y2",pnote,str,"=")
			str = StringByKey("Z2",pnote,"=")+","+StringByKey("Z2",wnote,"=")
			pnote = ReplaceStringByKey("Z2",pnote,str,"=")
			Note/K MoviePixelIntensities, pnote
		endif
	endif
	if (!WaveExists(PeakList) || !WaveExists(MoviePixelIntensities))
		return 1										// cannot find waves, not doing anything
	elseif (DimSize(PeakList,0)<1 || DimSize(PeakList,1)<2)
		return 1										// PeakList is invalid
	endif

	if (hw>0)
		Variable m,Nx=DimSize(image,0), Ny=DimSize(image,1), xlo,xhi,ylo,yhi
		for (m=0;m<Npixels;m+=1)
			px = round(PeakList[m][0])
			py = round(PeakList[m][1])
			xlo = limit(px-hw,0,Nx-1)
			xhi = limit(px+hw,0,Nx-1)
			ylo = limit(py-hw,0,Ny-1)
			yhi = limit(py+hw,0,Ny-1)
			ImageStats/M=1/G={xlo,xhi, ylo,yhi} image
			MoviePixelIntensities[ii][m] = V_avg
		endfor
	else
		MoviePixelIntensities[ii][] = image[PeakList[q][0]][PeakList[q][1]]
	endif
	return 0
End
Function/T FillMovieMorePixelsIntens_Init()
	return ""
End


// ****************************  Start of Fill Movie More Stuff ****************************
Function RescaleMoviePixelIntensities(MoviePixelIntensities,[printIt])
	Wave MoviePixelIntensities
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : !(!printIt)

	if (!WaveExists(MoviePixelIntensities))
		String wlist = WaveListClass("PixelIntensities","*","Dims:2")
		if (ItemsInList(wlist)==1)
			Wave MoviePixelIntensities = $StringFromList(0,wlist)
		else
			String name
			Prompt name,"Pick MoviePixelIntensities wave",popup,wlist
			DoPrompt "Movie Pixel Intensities", name
			if (V_flag)
				return 1
			endif
			Wave MoviePixelIntensities = $name
		endif
		printIt = 1
	endif
	name = NameOfWave(MoviePixelIntensities)
	if (!WaveExists(MoviePixelIntensities))
		return 1
	endif
	if (printIt)
		printf "RescaleMoviePixelIntensities(%s",name
		if (!ParamIsDefault(printIt))
			printf ", printIt=%g",printIt
		endif
		printf ")\r"
	endif

	name = GetWavesDataFolder(MoviePixelIntensities,1)+CleanupName(name[0,31-6],0)+"Scaled"
	Duplicate/O MoviePixelIntensities, $name/WAVE=scaled
	Note/K scaled, ReplaceStringByKey("waveClass",note(MoviePixelIntensities),"PixelIntensitiesScaled","=")
	scaled = NaN
	Variable Npixels=DimSize(MoviePixelIntensities,1), iLen=DimSize(MoviePixelIntensities,0)
	Make/FREE/D/N=(iLen) ww
	Variable m,maxVal
	for (m=0;m<Npixels;m+=1)
		ww = MoviePixelIntensities[p][m]
		maxVal = WaveMax(ww)
		ww /= maxVal
		scaled[][m] = ww[p]
	endfor

	if (printIt)
		printf "Created re-scaled wave \"%s\"\r",NameOfWave(scaled)
	endif
End

Function GraphMoviePixelIntensities(MoviePixelIntensities)
	Wave MoviePixelIntensities
	if (!WaveExists(MoviePixelIntensities))
		String wlist = WaveListClass("PixelIntensities*","*","Dims:2")
		if (ItemsInList(wlist)==1)
			Wave MoviePixelIntensities = $StringFromList(0,wlist)
		else
			String name
			Prompt name,"Pick MoviePixelIntensities wave",popup,wlist
			DoPrompt "Movie Pixel Intensities", name
			if (V_flag)
				return 1
			endif
			Wave MoviePixelIntensities = $name
		endif
	endif
	if (!WaveExists(MoviePixelIntensities))
		return 1
	endif
	name = NameOfWave(MoviePixelIntensities)

	String wnote = note(MoviePixelIntensities)
	Display/W=(16,66,752,761)
	Variable Npixels=DimSize(MoviePixelIntensities,1), m, px,py
	String Ltext="",str
	String pixelStr = StringByKey("PixelPositions",wnote,"=")
	for (m=0;m<Npixels;m+=1)
		AppendToGraph MoviePixelIntensities[*][m]

		str = StringFromList(m,pixelStr,",")
		px = str2num(StringFromList(0,str,":"))
		py = str2num(StringFromList(1,str,":"))
		if (m==0)
			sprintf str,"\\s(%s) [%g, %g]",name,px,py
		else
			sprintf str, "\r\\s(%s#%d) [%g, %g]",name,m,px,py
		endif
		Ltext += str
	endfor
	ModifyGraph tick=2, mirror=1, minor=1, lowTrip=0.001, 	zero(left)=2
	Label left "Intensity"
	Label bottom "Depth  (\\U)"

	Make/N=(7,3)/U/W/FREE colorCycle			// red, green blue, magenta, cyan, yellow, black
	colorCycle[0][0]= {65535,26205,1,65535,0,65535,0}
	colorCycle[0][1]= {0,52428,16019,0,65535,65535,0}
	colorCycle[0][2]= {0,1,65535,52428,65535,0,0}

	Make/W/FREE lSizes = {2,2,2,2,2,3,2}
	Variable i7, istyle
	for (m=0;m<Npixels;m+=1)
		i7 = mod(m,7)
		istyle = mod(floor(m/7),2) ? 11 : 0		// cycle between 11 & 0 every 7 traces
		ModifyGraph rgb[m]=(colorCycle[i7][0],colorCycle[i7][1],colorCycle[i7][2]), lsize[m]=lSizes[i7], lstyle[m]=istyle
	endfor
	Legend/C/N=text0/J/S=3/B=1/X=3.84/Y=2.76 Ltext
	return 0
End
// *****************************  End of Fill Movie More Stuff *****************************



Function FillMovieMoreCrosses(image,initStr)	// finds pixel intensities along a movie
	Wave image
	String initStr												// not used in this function

	Variable width=-1

	SetDrawLayer /K UserFront								// always clear first
	#if Exists("FixUpHolesInPeak")
		FixUpHolesInPeak(image)
	#endif
	width = width<0 ? round(min(DimSize(image,0),DimSize(image,1)) * 0.15) : width
	if (!FindAndFitGaussianCenter(image,width))	// include a box ±width/2 around the peak
		Wave W_coef=W_coef
		DrawMarker(W_coef[2],W_coef[4],width,width,"cross")
	endif

	return 0
End
//
// the result is passed to calling function by W_coef or {K0, K1, ... K6}
// this fits a 2-d gaussina to a region on an image.  The region is a square box of height=width
Function FindAndFitGaussianCenter(imageIn,width)	// returns 1 if error,  0 is OK
	Wave imageIn
	Variable width											// width of box to use for fitting
	if (!WaveExists(imageIn) || !(width>2))
		return 1
	endif

	Variable imax = DimSize(imageIn,0)-1, jmax = DimSize(imageIn,1)-1
	Variable/C zc=FindThePeakOnImage(imageIn)
	Variable xc=real(zc), yc=imag(zc)
	if (numtype(xc+yc))									// nothing there
		return 1
	endif

	Variable maxValue = imageIn(xc)(yc)*1.01		// maximum value to consider as valid when fitting the peak
	ImageStats/M=1/Q imageIn
	if (V_max>maxValue)									// there are bad hot pixels, remove them before fitting
		ImageThreshold/T=(maxValue) imageIn
		Wave M_ImageThresh=M_ImageThresh
		Duplicate/O imageIn imageFitTemp_JZT
		Redimension/D imageFitTemp_JZT
		Wave image = imageFitTemp_JZT
		image = M_ImageThresh[p][q] ? NaN : image[p][q]	// delete all hot pixels
		KillWaves/Z M_ImageThresh
	else
		Wave image = imageIn
	endif

	Variable V_left, V_right, V_bottom, V_top, w2=width/2
	V_left = limit(xc-w2,0,imax)
	V_right = limit(xc+w2,0,imax)
	V_top = limit(yc+w2,0,jmax)
	V_bottom = limit(yc-w2,0,jmax)

	// set up for the gaussian fit
	Variable i = FitOneGaussianPeak(image,V_left,V_right,V_bottom,V_top,Q=1)	// returns 1 if error,  0 is OK
	if (i)
		KillWaves/Z imageFitTemp_JZT
		return i
	endif

	Wave W_coef=W_coef
	Wave W_sigma=W_sigma
	Variable x0,y0,dx,dy,cor
	Variable fw= ( 2*sqrt(2*ln(2)) )				// this factor corrects for 2-d sigma to FWHM, apply to K3, K5, and the corresponding errors
	x0 = W_coef[2]
	dx = W_coef[3] * fw
	y0 = W_coef[4]
	dy = W_coef[5] * fw
	cor = W_coef[6]

	if (ItemsInList(GetRTStackInfo(0))<2)			// enables print out of result
		Variable pixels = (DimDelta(image,0)==1 && DimDelta(image,1)==1 && DimOffset(image,0)==0 && DimOffset(image,1)==0)
		Variable places1, places2, places3
		if (pixels)
			places1 = 2
			places2 = 2
			places3 = 2
		else
			Variable placesX = limit(ceil(log(limit(1/abs(W_sigma[2]),0.1,1e12))),0,12)
			Variable placesY = limit(ceil(log(limit(1/abs(W_sigma[4]),0.1,1e12))),0,12)
			places1 = max(placesX,placesY)
			placesX = limit(ceil(log(limit(1/abs(W_sigma[3]*fw),0.1,1e12))),0,12)
			placesY = limit(ceil(log(limit(1/abs(W_sigma[5]*fw),0.1,1e12))),0,12)
			places2 = max(placesX,placesY)
			places3 = limit(ceil(log(limit(1/abs(W_sigma[6]),0.1,1e12))),0,12)
		endif
		String fmt
		sprintf fmt, "Gaussian peak in '%%s' at (%%.%df±%%.%df, %%.%df±%%.%df) with FWHM of (%%.%df±%%.%df, %%.%df±%%.%df) and correlation of %%.%df±%%.%df  Amp=%%.3g\r",places1,places1,places1,places1,places2,places2,places2,places2,places3,places3
		printf fmt,NameOfWave(imageIn),x0,W_sigma[2],y0,W_sigma[4],dx,W_sigma[3]*fw,dy,W_sigma[5]*fw,cor,W_sigma[6],W_coef[1]
	endif
	KillWaves/Z imageFitTemp_JZT
	return 0
End
//
Static Function/C FindThePeakOnImage(image)					// locate the starting point for finding the peak, this should be pretty noise resistant
	Wave image
	if (!WaveExists(image))
		return cmplx(NaN,NaN)
	endif
	Duplicate/O image findPeakImage
	ImageFilter/N=5 median findPeakImage
	ImageStats/M=1 findPeakImage
	if (V_max<2)
		V_maxRowLoc = NaN
		V_maxColLoc = NaN
	endif
	KillWaves/Z findPeakImage
	return cmplx(V_maxRowLoc,V_maxColLoc)
End




//
// get the range of file index numbers by looking at all of the files in a directory


///
// get range of file numbers that are the first image of each wire scan
// the range of file index numbers by looking at all of the files in a directory
Static Function/T get_FirstWhiteFilesIndexRange(pathName,namePart)
	String pathName					// probably 'imagePath'
	String namePart					// name part of file = "WH_4174__"

	String list = directory(pathName)
	Variable i,m,N = ItemsInList(list)
	if (N<1)
		return ""
	endif
	Make/N=(N)/O/T fileList_get_FilesIndexRange
	Wave/T fileList = fileList_get_FilesIndexRange
	String name
	for (i=0,m=0;i<N;i+=1)
		name = StringFromList(i,list)
		if (strsearch(name,namePart,0)==0)
			fileList[m] = name						// file starts with namePart, save it in list
			m += 1									// increment pointer into fileList[]
		endif
	endfor
	N = m
	Redimension/N=(N) fileList						// set to correct length
	Make/N=(N)/I/O nw_get_FilesIndexRange
	Wave nw=nw_get_FilesIndexRange
	Variable k = strlen(namePart)
	for (m=0;m<N;m+=1)							// get wave with the file numbers
		name = fileList[m]
		name = name[k,100]
		nw[m] = str2num(name)
		i = strsearch(name,"_",0)
	endfor
	Sort nw,nw										// Sort the file numbers
	for (m=N-1;m>0;m-=1)						// remove duplicate numbers
		if (nw[m]==nw[m-1] && nw[m]<2.1e9)	// a duplicate, so delete it
			DeletePoints m,1,nw
		endif
	endfor
	N = numpnts(nw)

	// find which ones are first images of each wire scan, this is a lot of reading of headers, mark by changing new[m] to -1
	PathInfo $pathName
	String wnote, fileRoot=S_path+namePart
	Variable X1,H1,keV,H2
	Variable X1o,H1o,keVo,H2o
	String imageExtension=StrVarOrDefault("root:Packages:imageDisplay:imageExtension",".h5")
	wnote=ReadGenericHeader(fileRoot+num2istr(nw[0])+imageExtension)		// wave note to add to file read in
	if (!strlen(wnote))
		DoAlert 0,"file '"+name+"' does not have a complete header"
		KillWaves/Z fileList_get_FilesIndexRange, nw_get_FilesIndexRange
		return ""
	endif
	X1o=NumberByKey("X1", wnote,"=")
	H1o=NumberByKey("H1", wnote,"=")
	H2o=NumberByKey("H2", wnote,"=")
	keVo=NumberByKey("keV", wnote,"=")
	list = num2istr(nw[0])+";"					// first one is always first one of a scan
	String progressWin = ProgressPanelStart("",stop=1,showTime=1,status="finding first image of each wire scan")
	for (m=1;m<N;m+=1)
		if (mod(m,100)==0)
			if (ProgressPanelUpdate(progressWin,m/N*100))
				break											//   and break out of loop
			endif
		endif
		name= fileRoot+num2istr(nw[m])+imageExtension
		wnote=ReadGenericHeader(name)			// wave note to add to file read in
		if (!strlen(wnote))
			DoAlert 0,"file '"+name+"' does not have a complete header"
			KillWaves/Z fileList_get_FilesIndexRange, nw_get_FilesIndexRange
			return ""
		endif
		X1=NumberByKey("X1", wnote,"=")
		H1=NumberByKey("H1", wnote,"=")
		H2=NumberByKey("H2", wnote,"=")
		keV=NumberByKey("keV", wnote,"=")
		if ((abs(X1-X1o)+abs(H1-H1o))>0.12  || abs(keV-keVo) > 0.0002 || H2<H2o)	// start of new wire scan
			X1o = X1
			H1o = H1
			H2o = H2
			keVo = keV
			list += num2istr(nw[m])+";"
		endif
	endfor
	printf "time to find files = %s\r",Secs2Time(SecondsInProgressPanel(progressWin),5,0)
	DoWindow/K $progressWin
	String range = compressRange(list,";")
	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EscanButtonProc"))
		print "range	",ItemsInRange(range),"		",range
	endif
	KillWaves/Z fileList_get_FilesIndexRange, nw_get_FilesIndexRange
	return range
End
//
// get the range of file index numbers by looking at all of the files in a directory
Static Function/T get_FilesIndexRange(pathName,namePart)
	String pathName					// probably 'imagePath'
	String namePart					// name part of file = "WH_4174__"

	String list = directory(pathName)
	Variable i,m,N = ItemsInList(list)
	if (N<1)
		return ""
	endif
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

	list = ""							// make a list of all of the first number (with no repeats)
	k = NaN
	for (m=0;m<N;m+=1)
		i = nw[m]
		if (i!=k && i<4.2e+9)
			list += num2istr(i)+";"
			k = i
		endif
	endfor
	String range = compressRange(list,";")
	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EscanButtonProc"))
		print "range	",ItemsInRange(range),"		",range
	endif
	KillWaves/Z fileList_get_FilesIndexRange, nw_get_FilesIndexRange
	return range
End
//
//	returns the directory of the path in a semi-colon separated list (works for both Mac & Windows)
Static Function/T directory(pathName)
	String pathName							// path name that I want the directory of
	PathInfo $pathName	

	String list=""
	String imageExtension=StrVarOrDefault("root:Packages:imageDisplay:imageExtension",".h5")
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
			list = IndexedFile($pathName,-1,imageExtension)
		endif
	else											// this sectio for Windows
		list = IndexedFile($pathName,-1,imageExtension)		// list = IndexedFile($pathName,-1,"????")
	endif
	return list
End

Static Function/T FindMovieGraph()
	// Find name of suitable graph with imageOnMovie, and bring it to top
	Variable printIt = (ItemsInList(GetRTStackInfo(0))<2)
	Wave imageOnMovie=imageOnMovie
	String gList=FindGraphsWithWave(imageOnMovie), gName=""
	if (strlen(gList)<1)
		if (printIt)
			DoAlert 0, "No window ready for making the movie, first run\r'MakeMovieWindow()'"
		endif
		return ""
	endif
	Variable i
	for (i=0;i<ItemsInList(gList);i+=1)
		if (stringmatch(GetUserData(StringFromList(i,gList),"","windowType"),"OneMovieWindow"))
			gName = StringFromList(i,gList)
			break
		endif
	endfor
	if (strlen(gName)<1)
		if (printIt)
			DoAlert 0, "No suitable window ready for making the movie, first run\r'MakeMovieWindow()'"
		endif
		return ""
	endif
	return gName
End

//  ==================================== End of Movies ==================================== 
//  =======================================================================================



//  =======================================================================================
//	 ============================ Start of Intensity with Depth ============================

Function/T IntensityVsDepth(pathName,filePrefix,positions,depths,i0,i1,j0,j1)
	String pathName					// name of path to use
	String filePrefix					// first part of file name, e.g. "WH_"
	String positions, depths				// string ranges
	Variable i0,i1,j0,j1				// region of interest of image to sum

	String list = GetFileRootAndDoubleRange(pathName,filePrefix,positions,depths)
	if (strlen(list)<1)
		return ""
	endif
	pathName = StringByKey("pathName",list,"=")
	filePrefix = StringByKey("filePrefix",list,"=")
	positions = StringByKey("positions",list,"=")
	depths = StringByKey("depths",list,"=")
	String fileRoot = StringByKey("fileRoot",list,"=")
	Variable printIt = !(!NumberByKey("printIt",list,"="))
	if ((i0>i1 && (numtype(i0) || i1>=0)) || (j0>j1 && (numtype(j0) || j1>=0)))
		i0 = numtype(i0) ? 0 : max(round(i0),0)
		j0 = numtype(j0) ? 0 : max(round(j0),0)
		i1 = numtype(i1) ? -1 : i1
		j1 = numtype(j1) ? -1 : j1
		i1 = i1<i0 ? -1 : max(i0,round(i1))
		j1 = j1<j0 ? -1 : max(j0,round(j1)) 
		Prompt i0, "first index of ROI in x-direction"
		Prompt j0, "first index of ROI in y-direction"
		Prompt i1, "last index of ROI in x-direction (use -1 to indicate end)"
		Prompt j1, "last index of ROI in y-direction (use -1 to indicate end)"
		DoPrompt "ROI",i0,i1,j0,j1
		if (V_flag)
			return ""
		endif
		printIt = 1
	endif
	i1 = numtype(i1) ? -1 : i1
	j1 = numtype(j1) ? -1 : j1
	if (printIt)
		printf "IntensityVsDepth(\"%s\",\"%s\",\"%s\",\"%s\", %d,%d, %d,%d)\r",pathName,filePrefix,positions,depths,i0,i1,j0,j1
	endif
	if ((i0>i1 && (numtype(i0) || i1>=0)) || (j0>j1 && (numtype(j0) || j1>=0)))
		DoAlert 0, "Invalid ROI"
		return ""
	endif

	Variable epoch = DateTime, sec,secLast=0
	DoWindow/K DisplayLoopStatus
	Display /W=(371,50,678,183)/K=1
	DoWindow/C DisplayLoopStatus
	TextBox/N=file/F=0/A=LT/X=14.01/Y=16.54 "\\Z18starting"
	DoUpdate

	Variable Ndeep = ItemsInRange(depths)
	String wName = "intensVsDepth_full"
	Variable useROI=0
	if (!(i0==0 && j0==0 && i1<0 && j1<0))
		sprintf wName,"intensVsDepth_%d_%d_%d_%d",i0,j0,i1-i0+1,j1-j0+1
		useROI = 1
	endif
	Make/N=(Ndeep,1)/O $wName
	Wave intensVsDepth = $wName
//	Make/N=(Ndeep,1)/O intensVsDepth
	Variable ideep0=str2num(depths), ideep1=lastInRange(depths)
	SetScale/I x,ideep0,ideep1,"",intensVsDepth
	SetScale/P y,0,1,"",intensVsDepth
	Note/K intensVsDepth, ReplaceStringByKey("waveClass",note(intensVsDepth),"intensityVsDepth","=")

	Variable depthSi,keV,X1,Y1,Z1,H1
	String wnote, str
	String fname,imageName
	Variable pos,ipos, j,depth, N=0

	String imageExtension=StrVarOrDefault("root:Packages:imageDisplay:imageExtension",".h5")
	pos = strlen(positions) ? str2num(positions) :  Inf
	for (ipos=0; numtype(pos)<2; pos=NextInRange(positions,pos))
		if (numtype(pos))
			sprintf fname, "%s%d%s",fileRoot,str2num(depths),imageExtension
		else
			sprintf fname, "%s%d_%d%s",fileRoot,pos,str2num(depths),imageExtension
		endif
		GetFileFolderInfo/P=$pathName/Q/Z fname
		if (V_flag)								// file not found
			continue
		endif
		if (DimSize(intensVsDepth,1)<=ipos)								// wave is not long enough
			Redimension/N=(Ndeep,ipos+1) intensVsDepth					// increase size
		endif
		intensVsDepth[][ipos] = NaN

		for (depth=str2num(depths),j=0; !numtype(depth); depth=NextInRange(depths,depth))
			N += 1
			if (numtype(pos))
				sprintf fname, "%s%d%s",fileRoot,depth,imageExtension
			else
				sprintf fname, "%s%d_%d%s",fileRoot,pos,depth,imageExtension
			endif
//			imageName = WinViewReadROI(fName,i0,i1,j0,j1)			// load image file
			imageName = ReadGenericROI(fName,i0,i1,j0,j1)				// load image file
			if (exists(imageName)!=1)
				continue
			endif
			Wave image = $imageName
			wnote = note(image)
			intensVsDepth[j][ipos] = sum(image)/numpnts(image)
			KillWaves image
			sec = DateTime-epoch
			if (sec>(secLast+1))
				TextBox/C/W=DisplayLoopStatus/N=file "\\Z18"+imageName+"\relapsed time "+Secs2Time(sec,5,0)+"\r"+num2str(sec/N)+" sec/image"
				DoUpdate
				secLast = sec
			endif
		 	depthSi = NumberByKey("depth", wnote,"=")
		 	depthSi = numtype(depthSi) ? NumberByKey("depthSi", wnote,"=") : depthSi			// try the old name
	 		keV = NumberByKey("keV", wnote,"=")
			X1 = NumberByKey("X1", wnote,"=")
			Y1 = NumberByKey("Y1", wnote,"=")
			Z1 = NumberByKey("Z1", wnote,"=")
			H1 = YZ2H(Y1,Z1)
			j += 1
		endfor
		ipos += 1
	endfor
	if (useROI)
		wnote = ReplaceNumberByKey("DepthROI_i0",wnote,i0,"=")
		wnote = ReplaceNumberByKey("DepthROI_i1",wnote,i1,"=")
		wnote = ReplaceNumberByKey("DepthROI_j0",wnote,j0,"=")
		wnote = ReplaceNumberByKey("DepthROI_j1",wnote,j1,"=")
	endif
	wnote = ReplaceStringByKey("waveClass",wnote,"intensityVsDepth","=")
	Note/K intensVsDepth, wnote
	if (printIt || ItemsInList(GetRTStackInfo(0))<2 || stringmatch(GetRTStackInfo(0),"MarqueeGetIntensityVsDepth;*"))
		printf "filled '%s'  in   %s\r",GetWavesDataFolder(intensVsDepth,2),Secs2Time(DateTime-epoch,5,1)
	endif
	DoWindow/K DisplayLoopStatus
	if(ItemsInList(GetRTStackInfo(0))<2)		// only display result if called from command line or menu
		String win = StringFromList(0,FindGraphsWithWave(intensVsDepth))
		if (strlen(win)<1)
			if (DimSize(intensVsDepth,1)<2)
				Display /W=(5,44,497,346) intensVsDepth
				ModifyGraph gfMult=130,tick=2,lowTrip=0.001,standoff=0
			else
				Display /W=(418,68,961,466)
				AppendImage intensVsDepth
				ModifyImage intensVsDepth ctab= {*,*,Terrain,1}
			endif
			ModifyGraph mirror=1,minor=1
			ShowInfo
		else
			DoWindow/F $win
		endif
	endif
	return GetWavesDataFolder(intensVsDepth,2)
End
//
Static Function/T GetFileRootAndDoubleRange(pathName,filePrefix,positions,depths)
	String pathName					// name of path to use
	String filePrefix					// first part of file name, e.g. "WH_"
	String positions, depths				// string ranges

	String fileRoot=""
	Variable printIt=0
	if (strlen(pathName)<1 || strlen(positions)<1 || strlen(depths)<1)
		pathName = SelectString(strlen(pathName),"reconPath",pathName)
		Prompt pathName,"name of path that will point to the reconstructed images"
		Prompt positions,"range of positions (1st index) (non-existant are skipped), ['' -> no positions]"
		Prompt depths, "range of depths, the second index"
		DoPrompt "index ranges",pathName,positions,depths
		if (V_flag)
			return ""
		endif
		printIt = 1
	endif
	pathName = SelectString(stringmatch(pathName,CleanupName(pathName,0)),"",pathName)
	if (strlen(depths)<1)										// if only one range, make it depths
		depths = positions										// this way I don't depend upon user to get it right
		positions = ""
	endif
	if (!strlen(pathName) || ItemsInRange(depths)<1)
		DoAlert 0, "nothing to do, positions='"+positions+"',   depths='"+depths+"'  path='"+pathName+"'"
		return ""
	endif

	PathInfo $pathName
	String path = S_path
	if (strlen(filePrefix)<1 || !V_flag || strlen(S_path)<1)	// un-assigned path, or no file prefix
		Variable refNum
		String ImageFileFilters=StrVarOrDefault("root:Packages:imageDisplay:ImageFileFilters","Image Files:.*;")
		Open/F=ImageFileFilters/D/M="pick any one of the images"/R refNum
		if (strlen(S_fileName)<1)
			DoAlert 0,"no file specified"
			return ""
		endif
		filePrefix = ParseFilePath(3,S_fileName, ":", 1, 0)
		path = ParseFilePath(1,S_fileName,":",1,0)
		NewPath/O/Z $pathName,path
		if (strlen(positions))
			filePrefix = filePrefix[0,strsearch(filePrefix,"_",Inf,1)-1]
		endif
		filePrefix = filePrefix[0,strsearch(filePrefix,"_",Inf,1)]
		if (!strlen(filePrefix))
			DoAlert 0,"unable to identify file prefix"
			return ""
		endif
		printIt = 1
	endif
	String list = ReplaceStringByKey("fileRoot","",path+filePrefix,"=")
	list = ReplaceStringByKey("positions",list,positions,"=")
	list = ReplaceStringByKey("depths",list,depths,"=")
	list = ReplaceStringByKey("pathName",list,pathName,"=")
	list = ReplaceStringByKey("filePrefix",list,filePrefix,"=")
	list = ReplaceNumberByKey("printIt",list,printIt,"=")
	return list
End

//  =============================== End of Intensity with Depth ================================
//  ===================================================================================



//	===================================================================================
//	=============================== Start of Sum Images Using ROI ==============================

Function SumInManyROIs(pathName,namePart,range,mask,[dark])
	String pathName			// probably "imagePath"
	String namePart			//	 = "EW5_"
	String range
	Wave mask
	Wave dark												// an optional background wave

	pathName = SelectString(strlen(pathName),"imagePath",pathName)
	PathInfo $pathName
	String str
	if (!V_flag || strlen(namePart)<1)			// path does not exist or no namePart, ask user
		String pathPart
		str = requestFileRoot(pathName,1)
		pathPart = StringFromList(0,str)
		namePart = StringFromList(1,str)
		if (strlen(pathPart)<1 || strlen(namePart)<1)
			return 1							// invalid inputs
		endif
		if (!stringmatch(pathPart,S_path))		// path was changed
			if (stringmatch(pathName,"imagePath"))
				NewPath/O/M="path to reconstructed image files" imagePath pathPart	// for iamgePath, automatically reassign
			else
				NewPath/M="path to reconstructed image files" $pathName pathPart	// for other names, ask
			endif
		endif
	endif
	PathInfo $pathName
	if (strlen(S_path)<1 || strlen(namePart)<1)
		return 1								// invalid inputs
	endif
	String name,fileRoot=S_path+namePart

	Prompt range,"range of file numbers"
	DoPrompt "range", range
	if (V_flag)
		return 1
	endif
	String maskList = WaveListClass("imageMask","*","")
	if (!WaveExists(mask) && ItemsInList(maskList)>0)
		String maskName
		Prompt maskName,"Mask",popup,"_none_;"+maskList
		DoPrompt "Mask",maskName
		if (V_flag)
			return 1
		endif
		if (cmpstr(maskName,"_none_"))
			Wave mask = $maskName
		else
			Wave mask = $""
		endif
	endif
	String darkList = WaveListClass("imageDark;rawImageDark","*","")
	if (!WaveExists(dark) && ItemsInList(darkList)>0)
		String darkName
		Prompt darkName,"Dark Image",popup,"_none_;"+darkList
		DoPrompt "Background",darkName
		if (V_flag)
			return 1
		endif
		if (cmpstr(darkName,"_none_"))
			Wave dark = $darkName
		else
			Wave dark = $""
		endif
	endif
	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"DetailButtonProc"))
		if (WaveExists(mask))
			maskName = NameOfWave(mask)
		else
			maskName = "$\"\""
		endif
		if (WaveExists(dark))
			darkName = NameOfWave(dark)
		else
			darkName = "$\"\""
		endif
		printf "SumInManyROIs(\"%s\",\"%s\",\"%s\",%s",pathName,namePart,range,maskName
		if (WaveExists(dark))
			printf ",dark=%s",darkName
		endif
		printf ")\r"
		printf "using data from files starting with '%s'\r",namePart
	endif

	Variable Nimage = max(ItemsInRange(range),1),  i=str2num(range)
	String imageExtension=StrVarOrDefault("root:Packages:imageDisplay:imageExtension",".h5")
	sprintf name,"%s%d%s",fileRoot,i,imageExtension
	Wave image = $(LoadGenericImageFile(name))
	if (!WaveExists(image))
		printf "could not load first image named '%s'\r",name
		Abort "could not load first image"
	endif
	Variable Nx=DimSize(image,0),  Ny=DimSize(image,1)

	Make/N=(Nimage)/O/D sumIntensDepth=0
	Make/N=(Nimage)/O/D Resistance=NaN, VO2_Epoch=NaN, VO2_T=NaN
	if (WaveExists(mask))
		Make/N=(Nimage)/O/D sumIntensDepthROI=0
	endif
	SetScale/P x 0,1,"", sumIntensDepth

	name = namePart+"_Sum"+num2istr(ItemsInRange(range))
	name = CleanupName(ReplaceString("__",name,"_"),0)
	Duplicate/O image, $name
	Wave imageSum = $name
	Redimension/D imageSum
	imageSum = 0

	Variable j
	String wnote
	String progressWin = ProgressPanelStart("",stop=1,showTime=1)	// display a progress bar
	for (i=0,j=str2num(range); i<Nimage && !numtype(j); i+=1,j=NextInRange(range,j))
		if (ProgressPanelUpdate(progressWin,i/Nimage*100))	// update progress bar
			break												//   and break out of loop
		endif
		sprintf name,"%s%d%s",fileRoot,j,imageExtension
		Wave image = $(LoadGenericImageFile(name))
		wnote = note(image)
		Resistance[i] = NumberByKey("VO2_Resistance",wnote,"=")
		VO2_Epoch[i] = NumberByKey("VO2_Epoch",wnote,"=")
		VO2_T[i] = NumberByKey("VO2_Temp",wnote,"=")
		if (WaveExists(dark))
			Redimension/D image
			MatrixOP/O image = image - dark
		endif
		sumIntensDepth[i] = sum(image)
		if (WaveExists(mask))
			sumIntensDepthROI[i] = SumImageWithMask(image,mask)
		endif
		imageSum += image
		KillWaves/Z image
	endfor
	printf "total execution time = %s\r",Secs2Time(SecondsInProgressPanel(progressWin),5,0)
	DoWindow/K $progressWin
	WaveStats/Q/M=1 VO2_Epoch
	if (V_numNans==Nimage)
		KillWaves/Z Resistance, VO2_Epoch, VO2_T
	endif

	if (WaveExists(mask))
		MakeROIsumPlot()
	endif
	// beep
	return 0
End
//
Function MakeROIsumPlot()
	Wave sumIntensDepthROI=sumIntensDepthROI
	if (!WaveExists(sumIntensDepthROI))
		DoAlert 0, "You must run SumInManyROIs() with a mask first"
		return 1
	endif
	String wName = FindGraphsWithWave(sumIntensDepthROI)
	if (strlen(wName))
		wName = StringFromList(0,wName)
		DoWindow/F $wName
	else
		Display /W=(35,44,607,396) sumIntensDepthROI
		ModifyGraph tick=2,mirror=1,minor=1,lowTrip=0.001,standoff=0
	endif
	DoUpdate
	return 0
End
//
Static Function SumImageWithMask(image,mask)
	Wave image
	Wave mask					// This is a mask, pixels that I want are 1, unwanted are 0
	MatrixOP/FREE roi = greater(1,mask)	// create temp ROI, inverse of mask
	ImageStats/M=1/R=roi  image
	return V_npnts*V_avg
End


// This allows the user to select a file in pathName (or anywhere else), and removes the suffexes to get the proper file root
Static Function/T requestFileRoot(pathName,suffixes)
	String pathName
	Variable suffixes					// maximum number of suffixes to remove, each suffix looks like "_12", leave the underscore, use NaN for auto-determine
										// suffixes are ALWAYS preceeded with an "_"
	Variable killMe=0
	PathInfo $pathName
	if (V_flag==0)							// if pathname does not exist, use Desktop
		pathName = UniqueName("path",12,0)
		NewPath /Q/Z $pathName, SpecialDirPath("Desktop",0,0,0)
		killMe = 1
	endif
	String message="pick any reconstructed image file in range,  using datafolder = "+GetDataFolder(0)
	String ImageFileFilters=StrVarOrDefault("root:Packages:imageDisplay:ImageFileFilters","Image Files:.*;")
	Variable refNum
	Open/F=ImageFileFilters/D/M=message/P=$pathName/R refNum
	if (killMe)
		KillPath/Z $pathName
	endif
	if (strlen(S_fileName)<1)
		return ""
	endif

	String list = findFileRoot(S_fileName,suffixes)

	return StringFromList(1,list)+";"+StringFromList(2,list)
End


//	********** This next routine is really good, use it for all new stuff (and old stuff too) **********
//
// This takes a full path name, and parses the first N suffixes, it returns the list:
//	path  ;  root  ;  s1  ;  s2  ;  extenstion  ;  numberOfsuffixes
//	numbweOfSuffixes  ;  path  ;  root  ;  extension  ;  s1  ;  s2 ...
//	If you pass suffixex==NaN, then is auto determines the number of suffixes,  example:
//	findFileRoot("hd:folder:abc_123_666.h5",nan) -->   "2;hd:folder:;abc_;h5;123;666;"
Static Function/T findFileRoot(fullFilePath,suffixes)
	String fullFilePath
	Variable suffixes					// maximum number of suffixes to remove, each suffix looks like "_12", leave the underscore, use NaN for auto-determine
										// suffixes are ALWAYS preceeded with an "_"
	String path = ParseFilePath(1,fullFilePath,":",1,0)
	String fileRoot = ParseFilePath(3,fullFilePath,":",0,0)
	String extension = ParseFilePath(4,fullFilePath,":",0,0)

	// search to find number of suffixes
	suffixes = suffixes>=0 ? suffixes : Inf
	String sList=""
	Variable i0,i1, ns=0
	i1 = strlen(fileRoot)-1			// position of end
	do
		i0 = strsearch(fileRoot, "_",i1-1,1)
		if (i0<0)
			break
		endif
		Make/N=(i1-i0)/FREE/O flags=1	// check if last part is an integer
		flags = !isdigit(fileRoot[p+I0+1])
		if (sum(flags))
			break
		endif
		sList = fileRoot[i0+1,i1]+";"+sList
		ns += 1
		fileRoot = fileRoot[0,I0]
		i1 = i0-1
	while (ns<suffixes)
	return num2istr(ns)+";" + path+";" +  fileRoot+";" + extension+";" + sList
End



//	================================ End of Sum Images Using ROI ==============================
//	===================================================================================



// ==============================================================================================
// =========  The following group is for reading the ASCII intensity file "*_summary.txt" written after each reconstruction =========

Function/T LoadIntegralFile(fName)
	String fName					// full path name to the file

	Variable refNum
	PathInfo reconPath
	if (V_flag==0)
		Open/M="integral file from reconstruction"/R/T=".txt"/Z=2 refNum as fName
	else
		Open/M="integral file from reconstruction"/P=reconPath/R/T=".txt"/Z=2 refNum as fName
	endif
	if (V_flag)						// could not open file
		return ""
	endif
	close refNum
	fName = S_fileName
	PathInfo reconPath
	String path = ParseFilePath(1, S_fileName, ":", 1, 0)
	if (!stringmatch(S_path,path ))
		NewPath/O reconPath path								// force a re-set reconPath based on this
	endif
	String callingFn=GetRTStackInfo(2)
	Variable printIt = stringmatch(callingFn,"DetailButtonProc") || strlen(callingFn)==0

	String list = microGeo#keyStrFromFile(fName,"depthSortedInfo","reconPath")// read in all of the tagged values into a keyword list

	String item = StringByKey("array0",list,"=")
	if (strlen(item)<1)
		DoAlert 0,"$array0 not found in "+fName
		return ""
	endif

	Variable Nwaves,wavesLen
	Nwaves = str2num(StringFromList(0,item,","))
	wavesLen = str2num(StringFromList(1,item,","))

	Variable i,j
	String wNames="",wUnitsIn="",wUnits="", str, unit
	for (i=2;i<Nwaves+2;i+=1)
		unit = ""
		str = StringFromList(i,item,",")
		if (strsearch(str,")",strlen(str),1) == strlen(str)-1)		// ends with ")"
			j = strsearch(str,"(",strlen(str),1)
			unit = str[j+1,strlen(str)-2]
			str = str[0,j-1]
		endif

		str = CleanupName(str,0)
		if (CheckName(str,1))
			str = UniqueName(str,1,1)
		endif
		wNames += str+";"
		wUnitsIn += unit+";"
	endfor

	String columnInfoStr=""
	for (i=0;i<Nwaves;i+=1)
		str = StringFromList(i,wNames)
		if (stringmatch(str,"index"))
			columnInfoStr += "N='_skip_';  "
			continue
		endif
		columnInfoStr += "F=0,T=2,N="+str+";  "
		wUnits += StringFromList(i,wUnitsIn)+";"
	endfor

	Variable skipLines = countLinesToTag(fname,"array0")				// number of lines in file before $array0
	LoadWave/A/J/B=columnInfoStr/L={0,skipLines,0,0,0}/P=reconPath/Q fName
	if (V_flag<1)
		if (printIt)
			print "nothing loaded"
		endif
		return ""
	endif
	String listOfWaves=S_waveNames
	if (printIt)
		printf "loaded the %d waves  '%s'  from file '%s' in the path '%s'\r",V_flag,listOfWaves,S_fileName,S_path
	endif

	Wave wXaxis = $StringFromList(0,listOfWaves)
	if (!WaveExists(wXaxis))
		return ""
	endif

	for (i=0;i<ItemsInList(listOfWaves);i+=1)
		Wave ww = $StringFromList(i,listOfWaves)
		SetScale d 0,0,StringFromList(i,wUnits), ww
		Note/K ww, list
	endfor

	if (uniformlySpaced(wXaxis))
		if (printIt)
			printf "using the wave '%s' to set the x-scaling for all of the other waves (and then deleting it)\r",NameOfWave(wXaxis)
		endif
		unit = WaveUnits(wXaxis,-1)
		unit = SelectString(stringmatch(unit,"micron"),unit,"µm")
		unit = SelectString(stringmatch(unit,"degree"),unit,"¡")
		unit = SelectString(stringmatch(unit,"deg"),unit,"¡")
		Variable lo=wXaxis[0], hi=wXaxis[numpnts(wXaxis)-1]
		for (i=0;i<ItemsInList(listOfWaves);i+=1)
			Wave ww = $StringFromList(i,listOfWaves)
			SetScale/I x lo,hi,unit, ww
		endfor
		listOfWaves = RemoveFromList(NameOfWave(wXaxis),listOfWaves)
		KillWaves/Z wXaxis
		Wave wy = $StringFromList(0, listOfWaves)
		Wave wx = $""
	else
		Wave wy = $StringFromList(1, listOfWaves)
		Wave wx = $StringFromList(0, listOfWaves)
	endif

	if (stringmatch(WaveUnits(wy,0),"µm") && strsearch(NameOfWave(wy),"Intensity",0,2)==0)
		Note/K wy,ReplaceStringByKey("waveClass",note(wy),"intensityVsDepth","=")
	elseif ((stringmatch(WaveUnits(wx,-1),"µm") || stringmatch(WaveUnits(wx,-1),"micron")) && strsearch(NameOfWave(wy),"Intensity",0,2)==0)
		String wnote = ReplaceStringByKey("waveClass",note(wy),"intensityRecon","=")
		wnote = ReplaceStringByKey("Xwave",wnote,GetWavesDataFolder(wx,2),"=")
		Note/K wy, wnote
	endif
	if (printIt && strlen(FindGraphsWithWave(wy))<1)
		DoAlert 1, "Display the integral for this reconstruction?"
		if (V_flag==1)
			DisplayReconstructionIntegral(wy,wx)
		endif
	endif
	return listOfWaves
End
//
Static Function countLinesToTag(fname,tagName)	// returns the number of lines in file before $tag
	String fname										// full name of file
	String tagName									// tag to search for (does not include $)

	Variable refNum
	Open/M="file containing tagged values"/R/Z=1 refNum as fname
	FStatus refNum
	if (strlen(S_fileName)<1 || !refNum || !V_flag)
		return NaN
	endif
	String buf = PadString("",V_logEOF,0)
	FBinRead refNum, buf
	Close refNum

	Variable i, N=-1
	do
		N = strsearch(buf,"$"+tagName, N+1)
		if (char2num(buf[N+strlen(tagName)+1])<= 32)
			break
		endif
	while (N>0)
	if (N<0)											// $tag not found
		return NaN
	endif
	buf = ReplaceString("\r\n",buf[0,N-1],"\n")
	buf = ReplaceString("\n\r",buf,"\n")
	buf = ReplaceString("\r",buf,"\n")

	// count up the newlines
	for (i=0,N=0; i>=0; i = strsearch(buf,"\n",i+1),N+=1)
	endfor
	return N
End
//
Static Function uniformlySpaced(ww)		// returns true if values of ww are uniformly spaced
	Wave ww

	Variable thresh, diff = ww[1]-ww[0]
	switch(WaveType(ww) & 0x3F)
		case 0x00:					// not available for text
		case 0x01:					// not available for complex
			return 0
		case 0x02:
			thresh = diff*1e-5		// single precision float
			break
		case 0x04:
			thresh = diff*1e-12	// double precision float
			break
		case 0x08:					// all integer types
		case 0x10:
		case 0x20:
			thresh = 0.5			// for all integer types
			break
		default:						// should never reach this
			return 0
	endswitch
	Make/N=(numpnts(ww)-2)/B/FREE wtest
	wtest = abs((ww[p+1]-ww[p])-diff) > thresh
	Variable uniform = sum(wtest)<1
	return uniform
End
//
Function flattenIntegral(integWave)
	Wave integWave
	//	SetScale/P x -50,1,"µm", integWave
	Variable iend = numpnts(integWave)-6
	Variable i1=integWave[2], i2=integWave[iend]
	Variable slope = (i2-i1)/(iend-2)
	Variable b = i1-2*slope
	integWave -= p*slope + b
End

Function DisplayReconstructionIntegral(wy,wx)
	Wave wy,wx
	if (!WaveExists(wy))
		String list = WaveListClass("IntensityVsDepth,intensityRecon","*","DIMS:1")
		String wName = StringFromList(0,list)
		if (ItemsInList(list)<1)
			DoAlert 0, "No intensity vs depth profiles available, try:\r'IntensityVsDepth()' or 'LoadIntegralFile()'"
			return 1
		elseif (ItemsInList(list)>1)
			Prompt wName,"Intensity vs Depth",popup,list
			DoPrompt "Intensity vs Depth",wName
			if (V_flag)
				return 1
			endif
		endif
		Wave wy = $wName
	endif

	if (!WaveExists(wy))
		return 1
	endif
	String win = StringFromList(0,FindGraphsWithWave(wy))		// wave already plotted, bring to front
 	if (strlen(win))
		DoWindow/F $win
		return 0
 	endif

	String wnote=note(wy)
	if (!WaveExists(wx))
		Wave wx = $StringByKey("Xwave",wnote,"=")
	endif

	if (WaveExists(wx))
		Display /W=(3,365,472,622) wy vs wx
	else
		Display /W=(3,365,472,622) wy
	endif
	ModifyGraph gfMult=130, lSize=2, rgb=(0,5871,57708)
	ModifyGraph tick=2, zero(bottom)=2, mirror=1, minor=1, lowTrip=0.001, standoff=0
	ModifyGraph axOffset(left)=-0.857143,axOffset(bottom)=-0.538462
	Label left "Integral  (\\U)"
	Label bottom "Depth  (\\U)"
	WaveStats/Q/M=1 wy
	if (V_min<0)
		SetAxis/A/N=1 left
		ModifyGraph zero(left)=2
	else
		SetAxis/A/E=1 left
	endif

	String str, text="\\Z09"
	Variable i,  j0,j1,percent
	str = StringByKey("ws_infile",wnote,"=")
	i = strlen(ParseFilePath(1, str, "/", 1, 2))
	str = str[i,Inf]
	if (strlen(str))
		text += "in:\t"+str
	endif

	str = StringByKey("ws_outfile",wnote,"=")
	i = strlen(ParseFilePath(1, str, "/", 1, 2))
	str = str[i,Inf]
	if (strlen(str))
		text += "\rout:\t"+str
	endif

	str = StringByKey("ws_geofile",wnote,"=")
	i = strlen(ParseFilePath(1, str, "/", 1, 1))
	str = str[i,Inf]
	if (strlen(str))
		text += "\rgeo:\t"+str
	endif

	j0 = NumberByKey("ws_firstInputIndex",wnote,"=")
	j1 = NumberByKey("ws_lastInputIndex",wnote,"=")
	percent = NumberByKey("ws_percentOfPixels",wnote,"=")
	sprintf str,"\rusing images [%d, %d],  %g%% of pixels",j0,j1,percent
	text += str
	TextBox/C/N=label0/F=0/X=3/Y=6/B=1 text
	Cursor/P A $NameOfWave(wy) 0
	ShowInfo
	return 0
End


// ==============================================================================================
// ==============================================================================================
//  make vectors of most numbers in the header of a group of images
//
Function AllParams2Waves(range)
	String range
	Variable N = ItemsInRange(range)
	String fileRoot = "Macintosh HD:Users:tischler:dev:reconstructXcode_Mar07:FakeData:raw:W_"

	String fname, wnote
	String imageExtension=StrVarOrDefault("root:Packages:imageDisplay:imageExtension",".h5")
	sprintf fname,"%s%d%s",fileRoot,str2num(range),imageExtension		// open first image
//	wnote=WinViewReadHeader(fname)
	wnote=ReadGenericHeader(fname)

	String keyList="", wList="", item, key
	Variable i, value
	String skipList = "controllerType;imageFileName;imageFilePath;num_Type;numType;xdim;ydim;xDimDet;yDimDet;dateExposed;"
	skipList += "ADCrate;ADCtype;ADCresolution;geo_rotate;geo_reverse;geo_flip;startx;endx;groupx;starty;endy;groupy;x'"

	for (item=StringFromList(i,wnote),i=0; strlen(item); i+=1, item=StringFromList(i,wnote))
		value = str2num(StringFromList(1,item,"="))
		key = StringFromList(0,item,"=")
		if (WhichListItem(key,skipList)<0 && numtype(value)==0)
			keyList += key+";"
			wList += CleanupName(key,0)+";"
		endif
	endfor
	// printf "making waves    %s\r",wList
	for (i=0;i<ItemsInList(wList);i+=1)
		item = StringFromList(i,wList)
		Make/N=(N)/O $item = NaN
	endfor

	Variable m,j
	for(j=0,m=str2num(range); numtype(m)==0; m=NextInRange(range,m),j+=1)
		sprintf fname,"%s%d%s",fileRoot,m,imageExtension
//		wnote=WinViewReadHeader(fname)
		wnote=ReadGenericHeader(fname)
		if (!strlen(wnote))
			continue
		endif
		for (i=0;i<ItemsInList(wList);i+=1)
			item = StringFromList(i,wList)
			key = StringFromList(i,keyList)
			Wave ww = $item
			ww[j] = NumberByKey(key,wnote,"=")
		endfor
	endfor

	// printf "remove constant waves    %s\r",wList
	String killList=""
	for (i=0;i<ItemsInList(wList);i+=1)
		item = StringFromList(i,wList)
		Wave ww = $item
		WaveStats/Q ww
		if (V_max==V_min)
			KIllWaves/Z ww
			killList += item+";"
		endif
	endfor
	wLIst = RemoveFromList(killList, wList)

	if (WhichListItem("H2",wList)>=0)			// H2 is in list, so make dH2
		Make/N=(N-1)/O dH2=NaN
		Wave H2=H2
		dH2 = H2[p+1]-H2[p]
	endif

	if (ItemsInList(wList)>0)
		print "the waves that were changing are:   ",wList
	else
		print "NONE of the waves changed, they were all constant"
	endif
End



Function/T ImageMetaData2Waves(path,imageRoot,range,keyList,[mask,dark])	// returns list of wave names with the meta data from many images
	String path												// optional path
	String imageRoot										// root of image file
	String range											// a range of number to go with imageRoot
	String keyList											// space, comma, or semicolon separated key names
	Wave mask
	Wave dark												// an optional background wave

	String ImageFileFilters=StrVarOrDefault("root:Packages:imageDisplay:ImageFileFilters","Image Files:.*;")
	String imageExtension=StrVarOrDefault("root:Packages:imageDisplay:imageExtension",".h5")
	String fName=imageRoot+num2str(str2num(range))+imageExtension
	Variable f=0, printIt=0
	GetFileFolderInfo/P=$path/Q/Z fName
	if (V_Flag || !V_isFile)									// file not found
		Open/D/R/P=$path/F=ImageFileFilters/M="Select ANY image file in range" f
		if (strlen(S_fileName)<1)
			return ""
		endif
		String str = findFileRoot(S_fileName,1)
		imageRoot = StringFromList(1,str)+StringFromList(2,str)
		String extension = StringFromList(3,str)
		printIt = 1
	endif
	keyList = TrimFrontBackWhiteSpace(keyList)			// trim leading white space
	if (ItemsInRange(range)<1 || ItemsInList(keyList)<1)	// no range or no keys
		Prompt range,"range  e.g. 1-200"
		Prompt keyList,"list of keys (use spaces commas or semi-colons"
		DoPrompt "Range & Keys",range,keyList
		if (V_flag)
			return ""
		endif
		printIt = 1
	endif
	keyList = TrimFrontBackWhiteSpace(keyList)			// trim leading white space
	keyList = ReplaceString(" ",keyList,";")				// change sapces to semicolons
	keyList = ReplaceString(",",keyList,";")				// change commas to semicolons
	if (ItemsInList(keyList)<1)								// no keys, quit
		return ""
	endif
	fName=imageRoot+num2str(str2num(range))+imageExtension
	GetFileFolderInfo/P=$path/Q/Z fName
	if (V_Flag || !V_isFile)									// file not found, this is a test for range too
		return ""
	endif
	String maskList = WaveListClass("imageMask","*","")
	if (!WaveExists(mask) && ItemsInList(maskList)>0)
		String maskName
		Prompt maskName,"Mask",popup,"_none_;"+maskList
		DoPrompt "Mask",maskName
		if (V_flag)
			return ""
		endif
		if (cmpstr(maskName,"_none_"))
			Wave mask = $maskName
		else
			Wave mask = $""
		endif
	endif
	String darkList = WaveListClass("imageDark;rawImageDark","*","")
	if (WaveExists(mask) && !WaveExists(dark) && ItemsInList(darkList)>0)
		String darkName
		Prompt darkName,"Background Image",popup,"_none_;"+darkList
		DoPrompt "Background",darkName
		if (V_flag)
			return ""
		endif
		if (cmpstr(darkName,"_none_"))
			Wave dark = $darkName
		else
			Wave dark = $""
		endif
	endif
	if (printit)
		printf "ImageMetaData2Waves(\"%s\",\"%s\",\"%s\",\"%s\"",path,imageRoot,range,keyList
		if (WaveExists(mask))
			printf ",mask=%s",NameOfWave(mask)
		endif
		if (WaveExists(dark))
			printf ",dark=%s",NameOfWave(dark)
		endif
		printf ")\r"
	endif

	String key, wName, wList="",wListFull="", wnote, unit
	Variable N=ItemsInRange(range)
	Variable m,Nkeys = ItemsInList(keyList)
	Make/N=(Nkeys)/Wave/FREE keyWaves				// this holds the list of wave names
	for (m=0;m<Nkeys;m+=1)								// create and assign new waves
		key = StringFromList(m,keyList)
		wName = CleanupName(key+"_"+num2str(N),0)
		Make/N=(N)/D/O $wName
		Wave ww = $wName
#if Exists("LoadHDF5imageFile")==6
		unit = StringByKey(key,HDF5_MetaUnits)
#elif Exists("LoadWinViewFile")==6
		unit = StringByKey(key,WinView_MetaUnits)
#else
		unit = StringByKey(key,localKeyUnits)
#endif
		if (strlen(unit)>1)
			SetScale d 0,0,unit, ww
		endif
		keyWaves[m] = ww
		wListFull += GetWavesDataFolder(ww,2)+";"
		wList += NameOfWave(ww)+";"
		wnote = ""											// set wave note for each of the keys wave
		wnote = ReplaceStringByKey("keyName", wnote, key,"=")
		wnote = ReplaceStringByKey("range", wnote, range,"=")
		wnote = ReplaceStringByKey("imageRoot", wnote, imageRoot,"=")
		if (strlen(path))
			wnote = ReplaceStringByKey("path", wnote, path,"=")
		endif
		Note/K ww, wnote
	endfor
	if (WaveExists(mask))									// set up wave to store integral of mask
		wName = CleanupName("ROI_"+NameoFWave(mask),0)
		Make/N=(N)/D/O $wName=NaN						// store sum of roi in here
		Wave roi = $wName
		wListFull += GetWavesDataFolder(roi,2)+";"
		wList += NameOfWave(roi)+";"
		Duplicate/FREE mask, maskInvert					// do this because ImageStats only evaluates pixels=0
		maskInvert = !mask
	else
		Wave maskInvert=$""
	endif

	String progressWin = ProgressPanelStart("",stop=1,showTime=1)
	Variable i,j, Nmod=max(N/200,2)
	for (i=0,j=-Inf; i<N; i+=1)							// for each of the images
		if (mod(i,Nmod)==0)
			if (ProgressPanelUpdate(progressWin,i/N*100))
				break
			endif
		endif
		j = NextInRange(range,j)
		fName=imageRoot+num2str(j)+imageExtension
		if (WaveExists(maskInvert))						// a mask, get header info AND image
			Wave image = $LoadGenericImageFile(fName)
			if (WaveExists(image))
				wnote = note(image)
				if (WaveExists(dark))
					Redimension/I image
					MatrixOP/O image = image - dark
				endif
				ImageStats/M=1/R=maskInvert image
				roi[i] = V_avg*V_npnts
				KillWaves/Z image
			else
				wnote = ""
			endif
		else
			wnote = ReadGenericHeader(fName)				// get all header info from image
//			wnote = ReadHDF5header(fName)				// get all header info from image
		endif
		for (m=0;m<Nkeys;m+=1)							// for each of the keys
			Wave ww = keyWaves[m]						// wave to hold key value
			ww[i] = NumberByKey(StringFromList(m,keyList),wnote,"=")
		endfor
	endfor
	printf "total execution time = %s\r",Secs2Time(SecondsInProgressPanel(progressWin),5,0)
	DoWindow/K $progressWin
	if (printIt ||1)
		printf "Created & Filled waves:   %s\r",wList
	endif
	return wListFull										// return list of waves with keys
End



// ==============================================================================================
// =============================================  Init =============================================

Function initDepthResolve()
	initImageDisplayScaling()
End


