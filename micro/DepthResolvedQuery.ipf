#pragma rtGlobals=1		// Use modern global access method.
#pragma version = 1.20
#pragma IgorVersion = 6.0
#pragma ModuleName=depthResolve
#include "microGeometry", version>=2.48
#include "WinView", version>=1.90


Menu "micro"
	SubMenu "DepthResolved Query"
		"Plot Average Intensity vs Depth", IntensityVsDepth("","","","",0,NaN,0,NaN)
		"Make Movie of one type of Scan ",MovieOfOneScan()
		help={"movie of 1 wire scan, some sequential images, or first white image of each wire scan"}
		"Load recontruction 'summary' file", LoadIntegralFile("")
		"  Re-Display intensity from summary",DisplayReconstructionIntegral($"")
	End
End

Menu "GraphMarquee",dynamic
	"-"
	MarqueeDepthImageMenuItem("Show Average Intensity vs Depth"), MarqueeGetIntensityVsDepth()
	help={"sum up intensity withing marquee for each image of a depth reconstruction"}
End

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
	OK = OK || numtype(NumberByKey("depthSi",note(image),"="))==0
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
	filePrefix = ParseFilePath(3,filePrefix,":",1,0)						// remove extension (.SPE)
	for (i=strlen(filePrefix)-1; i>0 && isdigit(filePrefix[i,i]) ; i-=1)	// remove trailing digits
	endfor
	filePrefix = filePrefix[0,i]
	//	printf " x = [%g, %g],  y = [%g, %g]\r",ilo,ihi,jlo,jhi
	//	printf "'%s'       '%s'\r",pathStr,filePrefix
	IntensityVsDepth(pathName,filePrefix,"","",ilo,ihi,jlo,jhi)
End



Function MovieOfOneScan()		// make the movie frame and fill it
	String gName=FindMovieGraph()
	if (strlen(gName))
		DoWindow/F $gName
		FillMovieOfOneScan("","","",NaN,NaN,"",0)
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
	FillMovieOfOneScan("","","",NaN,NaN,"",0)
End

Function/T MakeMovieWindow()
	// pick a WinView image to use for the movie
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
		wName = LoadWinViewFile("")
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
	TextBox/C/N=textParameters/A=RB ""
	TextBox/C/N=textFileName/F=0/B=1/A=LB ""
	SetWindow kwTopWin,userdata(windowType)=  "OneMovieWindow"
	DoUpdate
	SetAspectToSquarePixels("")	
//	DoAlert 0,"Adjust color scale on this image and adjust other items before starting movie"
	return WinName(0,1,1)
End
//
Function FillMovieOfOneScan(pathName,filePrefix,range,surface,absorpLength,type,crosses,[makeSum])
	String pathName					// name of path to use
	String filePrefix					// first part of file name, e.g. "WH_4174_"
	String range						// range of files to use
	Variable surface					// depth of surface relative to depthSi(µm), used to calculate boost
	Variable absorpLength			// absorption length (µm)
	String type						// type of movie, images of a wire scan, sequenial images, first white image of each wire scan
	Variable crosses					// flag, TRUE put cross on peak, FALSE no cross
	Variable makeSum
	makeSum = ParamIsDefault(makeSum) ? 0 : (numtype(makeSum)==0 && makeSum != 0)

	String typeList="single wire scan;sequential images;first white image of each wirescan"
	Variable useLog=0

	String gName = FindMovieGraph()								// name of suitable graph with imageOnMovie
	DoWindow/F $gName												// and bring it to the front

	Variable itype = WhichListItem(type, typeList)
	if (itype<0 || numtype(crosses))								// no type passed, or invalid
		Prompt type, "type of movie to make", popup, typeList
		Prompt crosses,"Put fitted crosses on the bigest spot",popup,"NO crosses;show fitted crosses"
		crosses = crosses==0 || crosses==1 ? crosses+1 : 0
		Prompt makeSum,"Also make a sum image",popup,"No Sum;Make Sum Too"
		DoPrompt "movie type",type,crosses,makeSum
		if (V_flag)
			return 1
		endif
		makeSum -= 1
		crosses -= 1
	endif
	itype = WhichListItem(type, typeList)

	GetFileFolderInfo/P=$pathName/Q/Z filePrefix+num2istr(str2num(range))+".SPE"
	if (!V_isFile)
		GetFileFolderInfo/P=$pathName/Q filePrefix+num2istr(str2num(range))+".SPE"
	endif
	if (!V_isFile)
		return 1
	endif
	pathName = SelectString(strlen(pathName),"reconPath",pathName)	// use reconPath as the default path name, if none passed
	String path = ParseFilePath(1,S_Path,":",1,0)
	NewPath/O/Q $pathName, ParseFilePath(1,S_Path,":",1,0)
	filePrefix = ParseFilePath(3,S_Path,":",1,0)							// remove extension (.SPE)
	Variable i
	for (i=strlen(filePrefix)-1; i>0 && isdigit(filePrefix[i,i]) ; i-=1)// remove trailing digits
	endfor
	filePrefix = filePrefix[0,i]
	String fileRoot = path+filePrefix

	if (ItemsInRange(range)<1)												// find the range here
		if (itype <=1 )														// for itype == 0 or 1
			range = get_FilesIndexRange(pathName,filePrefix)
		elseif (itype ==2 )													// first white image of each wirescan, itype==2
			range = get_FirstWhiteFilesIndexRange(pathName,filePrefix)
		else
			Abort "in FillMovieOfOneScan() itype = "+num2str(itype)
		endif
	endif
	if (ItemsInRange(range)<1)
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
	if (exists("FillMovieMore")==6)						// hook for a user supplied function to do some more
		DoAlert 2,"FillMovieMore exists, use it?"
		if (V_flag==1)									// "Yes" was clicked
			FUNCREF MoreInFillMovieProto func=$"FillMovieMore"
		elseif(V_flag==2)									// "No" was clicked
			FUNCREF MoreInFillMovieProto func=MoreInFillMovieProto
		else													// "Cancel was clicked
			return 1
		endif
	endif
	printf "FillMovieOfOneScan(\"%s\",\"%s\",\"%s\",%g,%g,\"%s\",%d)\r",pathName,filePrefix,range[0,100],surface,absorpLength,type,crosses

	NewMovie/P=home
	printf "running in data folder  '%s'\r",GetDataFolder(1)
	Variable boost													// boost factor from absorption
	Variable depthSi,integral, epoch=DateTime
	Variable X1,H1,keV
	Variable X10,H10,keV0
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

	Make/N=(ItemsInRange(range))/O MovieIntegral=NaN	// holds the integral for each frame of the movie
	Variable ii

	String str,textStr, fname,imageName
	Variable width=-1
	for (i=str2num(range),ii=0; !numtype(i); i=NextInRange(range,i),ii+=1)
		sprintf fname, "%s%d.SPE",fileRoot,i
		imageName = WinViewReadROI(fname,i0,i1,j0,j1)	// load file into wName
		if (exists(imageName)!=1)
			continue
		endif
		Wave image = $imageName
		if (ii==0 && !(WaveType(image)&0x40))				// if the first frame is signed, ensure imageOnMovie is signed
			Duplicate/O image imageOnMovie
			if (makeSum)
				Duplicate/O imageOnMovie imageOnMovieSum
				Redimension/D imageOnMovieSum
				imageOnMovieSum = 0
			endif
		endif
		if (ii==0)
			Duplicate/O imageOnMovie imageOnMovieSum
			imageOnMovieSum = 0
			Redimension/D imageOnMovieSum
		endif
		wnote = note(image)
	 	depthSi = NumberByKey("depthSi", wnote,"=")			// values read from image file
		X1 = NumberByKey("X1", wnote,"=")
		H1 = NumberByKey("H1", wnote,"=")
		keV = NumberByKey("keV", wnote,"=")
		if (ii==0)													// save first values for comparison
			X10=X1 ;  H10=H1 ;  keV0=KeV
		endif
		integral = sum(image)
		if (integral==0)											// image is empty
			KillWaves image
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
		Note/K imageOnMovie, ReplaceNumberByKey("boost",note(imageOnMovie),boost,"=")
//		imageOnMovie = image*boost
//		if (useLog)
//			imageOnMovie = 5000*log(imageOnMovie)			// need 5000 because image is integer
//		endif
		if (crosses)
			SetDrawLayer /K UserFront							// always clear first
			#if Exists("FixUpHolesInPeak")
				FixUpHolesInPeak(image)
			#endif
			width = width<0 ? round(min(DimSize(image,0),DimSize(image,1)) * 0.15) : width
			if (!FindAndFitGaussianCenter(image,width))		// include a box ±width/2 around the peak
				Wave W_coef=W_coef
				DrawCross(W_coef[2],W_coef[4],width,width)
			endif
		endif
		KillWaves/Z image

		textStr = ""
		if (numtype(depthSi)==0)
			sprintf textStr, "%+.1f µm"+SelectString(useLog,"","\rLog")+"\r",depthSi
		endif
//		if (itype==2 && numtype(X1+H1+keV)==0)
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
//		TextBox/C/N=textFileName/F=0/B=1/A=LT textStr	// "filesname\r\\Zr070file path"
		TextBox/C/N=textFileName/B=1 textStr				// "filesname\r\\Zr070file path"
//		func(boost)
		func(imageOnMovie)
		imageOnMovie *= boost
		if (useLog)
			imageOnMovie = 5000*log(imageOnMovie)			// need 5000 because image is integer
		endif
		DoUpdate
		AddMovieFrame
	endfor
	CloseMovie
	print "done, total execution time is  ",Secs2Time(DateTime-epoch,5,0)
	return 0
End
//
//Function FillMovieOfOneScan(pathName,filePrefix,range,surface,absorpLength,type,crosses)
//	String pathName					// name of path to use
//	String filePrefix					// first part of file name, e.g. "WH_4174_"
//	String range						// range of files to use
//	Variable surface					// depth of surface relative to depthSi(µm), used to calculate boost
//	Variable absorpLength			// absorption length (µm)
//	String type						// type of movie, images of a wire scan, sequenial images, first white image of each wire scan
//	Variable crosses					// flag, TRUE put cross on peak, FALSE no cross
//
//	String typeList="single wire scan;sequential images;first white image of each wirescan"
//	Variable useLog=0
//
//	String gName = FindMovieGraph()								// name of suitable graph with imageOnMovie
//	DoWindow/F $gName												// and bring it to the front
//
//	Variable itype = WhichListItem(type, typeList)
//	if (itype<0 || numtype(crosses))								// no type passed, or invalid
//		Prompt type, "type of movie to make", popup, typeList
//		Prompt crosses,"Put fitted crosses on the bigest spot",popup,"NO crosses;show fitted crosses"
//		crosses = crosses==0 || crosses==1 ? crosses+1 : 0
//		DoPrompt "movie type",type,crosses
//		if (V_flag)
//			return 1
//		endif
//		crosses -= 1
//	endif
//	itype = WhichListItem(type, typeList)
//
//	GetFileFolderInfo/P=$pathName/Q/Z filePrefix+num2istr(str2num(range))+".SPE"
//	if (!V_isFile)
//		GetFileFolderInfo/P=$pathName/Q filePrefix+num2istr(str2num(range))+".SPE"
//	endif
//	if (!V_isFile)
//		return 1
//	endif
//	pathName = SelectString(strlen(pathName),"reconPath",pathName)	// use reconPath as the default path name, if none passed
//	String path = ParseFilePath(1,S_Path,":",1,0)
//	NewPath/O/Q $pathName, ParseFilePath(1,S_Path,":",1,0)
//	filePrefix = ParseFilePath(3,S_Path,":",1,0)							// remove extension (.SPE)
//	Variable i
//	for (i=strlen(filePrefix)-1; i>0 && isdigit(filePrefix[i,i]) ; i-=1)// remove trailing digits
//	endfor
//	filePrefix = filePrefix[0,i]
//	String fileRoot = path+filePrefix
//
//	if (ItemsInRange(range)<1)												// find the range here
//		if (itype <=1 )														// for itype == 0 or 1
//			range = get_FilesIndexRange(pathName,filePrefix)
//		elseif (itype ==2 )													// first white image of each wirescan, itype==2
//			range = get_FirstWhiteFilesIndexRange(pathName,filePrefix)
//		else
//			Abort "in FillMovieOfOneScan() itype = "+num2str(itype)
//		endif
//	endif
//	if (ItemsInRange(range)<1)
//		return 1
//	endif
//	if ((numtype(surface)==2 || !(absorpLength>0)) && itype==0)
//		absorpLength = !(absorpLength>0) ? 15 : absorpLength
//		Prompt surface, "depth of surface relative to depthSi(µm), (use NaN to ignore)"
//		Prompt absorpLength, "absorption length in material (µm), (use NaN to ignore)"
//		DoPrompt "absorption correction",surface,absorpLength
//		if (V_flag)
//			return 1
//		endif
//	endif
//	if (exists("FillMovieMore")==6)						// hook for a user supplied function to do some more
//		DoAlert 2,"FillMovieMore exists, use it?"
//		if (V_flag==1)									// "Yes" was clicked
//			FUNCREF MoreInFillMovieProto func=$"FillMovieMore"
//		elseif(V_flag==2)									// "No" was clicked
//			FUNCREF MoreInFillMovieProto func=MoreInFillMovieProto
//		else													// "Cancel was clicked
//			return 1
//		endif
//	endif
//	printf "FillMovieOfOneScan(\"%s\",\"%s\",\"%s\",%g,%g,\"%s\",%d)\r",pathName,filePrefix,range[0,100],surface,absorpLength,type,crosses
//
//	NewMovie/P=home
//	printf "running in data folder  '%s'\r",GetDataFolder(1)
//	Variable boost													// boost factor from absorption
//	Variable depthSi,integral, epoch=DateTime
//	Variable X1,H1,keV
//	Variable X10,H10,keV0
//	Wave imageOnMovie=imageOnMovie
//	String wnote = note(imageOnMovie)
//	Variable i0=0, i1=-1, j0=0, j1=-1
//	// check if we are only using a sub-image for the movie
//	if (DimSize(imageOnMovie,0)!=NumberByKey("xdim", wnote,"=") || DimSize(imageOnMovie,1)!=NumberByKey("ydim", wnote,"="))
//		i0 = NumberByKey("startxRead", wnote,"=")
//		j0 = NumberByKey("startyRead", wnote,"=")
//		i0 = i0>0 ? i0 : 0
//		j0 = j0>0 ? j0 : 0
//		i1 = i0 + DimSize(imageOnMovie,0) - 1
//		j1 = j0 + DimSize(imageOnMovie,1) - 1
//		printf "using a sub-region of the images [%d, %d] [%d, %d]\r",i0,i1,j0,j1
//	endif
//
//	Make/N=(ItemsInRange(range))/O MovieIntegral=NaN	// holds the integral for each frame of the movie
//	Variable ii
//
//	String str,textStr, fname,imageName
//	Variable width=-1
//	for (i=str2num(range),ii=0; !numtype(i); i=NextInRange(range,i),ii+=1)
//		sprintf fname, "%s%d.SPE",fileRoot,i
//		imageName = WinViewReadROI(fname,i0,i1,j0,j1)	// load file into wName
//		if (exists(imageName)!=1)
//			continue
//		endif
//		Wave image = $imageName
//		if (ii==0 && !(WaveType(image)&0x40))				// if the first frame is signed, ensure imageOnMovie is signed
//			Duplicate/O image imageOnMovie
//		endif
//		wnote = note(image)
//	 	depthSi = NumberByKey("depthSi", wnote,"=")			// values read from image file
//		X1 = NumberByKey("X1", wnote,"=")
//		H1 = NumberByKey("H1", wnote,"=")
//		keV = NumberByKey("keV", wnote,"=")
//		if (ii==0)													// save first values for comparison
//			X10=X1 ;  H10=H1 ;  keV0=KeV
//		endif
//		integral = sum(image)
//		if (integral==0)											// image is empty
//			KillWaves image
//			continue
//		endif
//		MovieIntegral[ii] = integral
//		boost = exp(max(depthSi-surface,0)/absorpLength)
//		boost = numtype(boost) ? 1 : boost
//		boost = itype==0 ? boost : 1
//		imageOnMovie = image
//		Note/K imageOnMovie, ReplaceNumberByKey("boost",note(imageOnMovie),boost,"=")
////		imageOnMovie = image*boost
////		if (useLog)
////			imageOnMovie = 5000*log(imageOnMovie)			// need 5000 because image is integer
////		endif
//		if (crosses)
//			SetDrawLayer /K UserFront							// always clear first
//			#if Exists("FixUpHolesInPeak")
//				FixUpHolesInPeak(image)
//			#endif
//			width = width<0 ? round(min(DimSize(image,0),DimSize(image,1)) * 0.15) : width
//			if (!FindAndFitGaussianCenter(image,width))		// include a box ±width/2 around the peak
//				Wave W_coef=W_coef
//				DrawCross(W_coef[2],W_coef[4],width,width)
//			endif
//		endif
//		KillWaves/Z image
//
//		textStr = ""
//		if (numtype(depthSi)==0)
//			sprintf textStr, "%+.1f µm"+SelectString(useLog,"","\rLog")+"\r",depthSi
//		endif
////		if (itype==2 && numtype(X1+H1+keV)==0)
//		if (numtype(X1+H1+keV)==0 && (itype>0 || ii==0))
//			sprintf str, "X1 = %.2f\rH1 = %.2f\rkeV = %.4f\r",X1,H1,keV
//			textStr += str
//			str = ""
//		endif
//		if (integral==0)
//			sprintf str, "\\Zr070(%d), º=0",i
//		elseif (boost>1)
//			sprintf str, "\\Zr070(%d), º=%.2fE6\rboost=%.1f",i,integral/1e6,boost
//		else
//			sprintf str, "\\Zr070(%d), º=%.2fE6",i,integral/1e6
//		endif
//		textStr += str
//		TextBox/C/N=textParameters textStr
//	 	str = StringByKey("imageFileName", wnote,"=")
//		textStr = SelectString(strlen(str),"",str)
//	 	str = StringByKey("imageFilePath", wnote,"=")
//		if (strlen(str))
//			textStr += SelectString(strlen(textStr),"","\r") + "\\Zr070" + str
//		endif
////		TextBox/C/N=textFileName/F=0/B=1/A=LT textStr	// "filesname\r\\Zr070file path"
//		TextBox/C/N=textFileName/B=1 textStr				// "filesname\r\\Zr070file path"
////		func(boost)
//		func(imageOnMovie)
//		imageOnMovie *= boost
//		if (useLog)
//			imageOnMovie = 5000*log(imageOnMovie)			// need 5000 because image is integer
//		endif
//		DoUpdate
//		AddMovieFrame
//	endfor
//	CloseMovie
//	print "done, total execution time is  ",Secs2Time(DateTime-epoch,5,0)
//	return 0
//End
Function MoreInFillMovieProto(image)
	Wave image
End
//		// This is an example:
//Function FillMovieMore(image)
//	Wave image
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
	wnote=WinViewReadHeader(fileRoot+num2istr(nw[0])+".SPE")			// wave note to add to file read in
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
	for (m=1;m<N;m+=1)
		name= fileRoot+num2istr(nw[m])+".SPE"
		wnote=WinViewReadHeader(name)			// wave note to add to file read in
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
	String range = compressRange(list,";")
	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
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
	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
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
			list = IndexedFile($pathName,-1,".SPE")
		endif
	else											// this sectio for Windows
		list = IndexedFile($pathName,-1,".SPE")		// list = IndexedFile($pathName,-1,"????")

	endif
	return list
End

Static Function/T FindMovieGraph()
	// Find name of suitable graph with imageOnMovie, and bring it to top
	Variable printIt = (ItemsInList(GetRTStackInfo(0))<2)
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







Function/T IntensityVsDepth(pathName,filePrefix,positions,depths,i0,i1,j0,j1)
	String pathName					// name of path to use
	String filePrefix					// first part of file name, e.g. "WH_"
	String positions, depths				// string ranges
	Variable i0,i1,j0,j1				// region of interest of image to sum

	String list = GetFileRootAndDoubleRange(pathName,filePrefix,positions,depths)
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

	pos = strlen(positions) ? str2num(positions) :  Inf
	for (ipos=0; numtype(pos)<2; pos=NextInRange(positions,pos))
		if (numtype(pos))
			sprintf fname, "%s%d.SPE",fileRoot,str2num(depths)
		else
			sprintf fname, "%s%d_%d.SPE",fileRoot,pos,str2num(depths)
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
				sprintf fname, "%s%d.SPE",fileRoot,depth
			else
				sprintf fname, "%s%d_%d.SPE",fileRoot,pos,depth
			endif
			imageName = WinViewReadROI(fName,i0,i1,j0,j1)				// load image file
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
	 		depthSi = NumberByKey("depthSi", wnote,"=")
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
		Prompt pathName,"name of path that points to the reconstructed images"
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
		depths = positions
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
		Open/D/M="pick any one of the images"/R/T=".spe" refNum
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

	Variable printIt = ItemsInList(GetRTStackInfo(0))<2

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
	endif

	Wave ww = $StringFromList(0, listOfWaves)
	if (stringmatch(WaveUnits(ww,0),"µm") && strsearch(NameOfWave(ww),"Intensity",0,2)==0)
		Note/K ww,ReplaceStringByKey("waveClass",note(ww),"intensityVsDepth","=")
	endif
	if (printIt && strlen(FindGraphsWithWave(ww))<1)
		DoAlert 1, "Display the integral for this reconstruction?"
		if (V_flag==1)
			DisplayReconstructionIntegral(ww)
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

	String wName = UniqueName("uniformTest",1,0)
	Make/N=(numpnts(ww)-1)/B $wName
	Wave wtest = $wName
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
	wtest = abs((ww[p+1]-ww[p])-diff) > thresh
	Variable uniform = sum(wtest)<1
	KillWaves/Z wtest
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

Function DisplayReconstructionIntegral(ww)
	Wave ww
	if (!WaveExists(ww))
		String list = WaveListClass("IntensityVsDepth","*","DIMS:1")
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
		Wave ww = $wName
	endif

	if (!WaveExists(ww))
		return 1
	endif
	String win = StringFromList(0,FindGraphsWithWave(ww))		// wave already plotted, bring to front
 	if (strlen(win))
		DoWindow/F $win
		return 0
 	endif

	WaveStats/Q/M=1 ww

	Display /W=(3,365,472,622) ww
	ModifyGraph gfMult=130, lSize=2, rgb=(0,5871,57708)
	ModifyGraph tick=2, zero(bottom)=2, mirror=1, minor=1, lowTrip=0.001, standoff=0
	ModifyGraph axOffset(left)=-0.857143,axOffset(bottom)=-0.538462
	Label left "Integral  (\\U)"
	Label bottom "Depth  (\\U)"
	if (V_min<0)
		SetAxis/A/N=1 left
		ModifyGraph zero(left)=2
	else
		SetAxis/A/E=1 left
	endif

	String str, wnote = note(ww)
	String text="\\Z09"
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
	percent = NumberByKey("ws_percentOfPixells",wnote,"=")
	sprintf str,"\rusing images [%d, %d],  %g%% of pixels",j0,j1,percent
	text += str
	TextBox/C/N=label0/F=0/X=3/Y=6 text
	Cursor/P A $NameOfWave(ww) 0
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
	sprintf fname,"%s%d.SPE",fileRoot,str2num(range)		// open first image
	wnote=WinViewReadHeader(fname)

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
		sprintf fname,"%s%d.SPE",fileRoot,m
		wnote=WinViewReadHeader(fname)
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