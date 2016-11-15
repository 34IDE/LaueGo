#pragma rtGlobals=1		// Use modern global access method.
#pragma version = 2.06
#pragma ModuleName=WinViewProc
#include "ImageDisplayScaling", version>=1.98
//
// Routines for reading in and looking at Princeton Instruments CCD, by Jon Tischler, Oak Ridge National Lab
// TischlerJZ@ornl.gov
//
// to use these macros in your Igor experiment, put this file in the "User Procedures" folder
// in the "Igor Pro Folder", and uncomment the following line and put it in the top of the procedure
//folder
//
//#include "WinViewProcedures"
//
//
// version 1.3 (changed Mar 1, 2006)
//		Changed WinViewReadROI().  Separated the reading of the header into a new function named WinViewReadHeader() that
//		returns a key list of all of the header information.  Functionally WinViewReadROI() is unchanged, but the new routine
//		WinViewReadHeader() can be called separately if you only want the header information without reading in any image data.
//
// version 1.4 (changed Mar 27, 2006)
//		Changed FitOneGaussianPeak() a better check for imageName, before using it.
//		Changed default color table from PlanetEarth to Terrain
//		Changed default graph to be without buttons (is was with buttons)
//
// version 1.5 (changed April 7, 2006)
//		Changed WinViewInfo() to use wave input, not string input
//
// version 1.6 (changed April 28, 2006)
//		Changed WinViewReadHeader() to use a struct to read in header information, rather than reading in a block of uint16
//		Added the ability to write an spe file
//		Changed FitOneGaussianPeak() and GetCenterOfMass() so that the ROI is constrained to be withing the image
//		Removed the constraint on the cross-correlation coefficient for gauss2D, not needed for Igor 5
//
//	version 1.7 (changed August 8, 2006)
//		Added the 'waveClass' keyword to the wave note of the loaded image
//		Removed the OnlyFileName() function, replaced it with ParseFilePath() and CleanupName(), so now windows & Mac compatible
//		Addded additional constraints on Gauss2D fit, K3 & K5 > 0
//
//		Also changed WinViewInfo() to use WaveListClass("speImage*","*","DIMS:2")
//
//	version 1.8 (changed October 20, 2006)
//		changed FitOneGaussianPeak()
//		It now understands images that are scaled or not scaled
//
//		Also changed WinViewInfo() to use WaveListClass("speImage*","*","DIMS:2")
//
//	version 1.81 (changed November 1, 2006)
//		changed WinView menu & marquee menu, made them dynamic so only appropriate items show
//		depending upon the top image's waveClass
//
//	version 1.82 (changed December 21, 2006)
//		changed WinViewReadHeader() and WinViewLoadROI() so that it returns early if the file is too short
//		changed getZrange() to avoid error when the min and max value are the same
//
//	version 1.83(changed February 25, 2007)
//		added imageFilePath to the file path when reading in images in WinViewReadHeader()
//
//	version 1.84(changed February 26, 2007)
//		fixed WinViewLoadROI() so that it works right if j0==j1
//
//	version 1.85(changed June 6, 2007)
//		changed FitOneGaussianPeak() so that it also reports an error if fit parameters look crazy
//
//	version 1.86(changed June 18, 2007)
//		changed GraphImageStyle() to deal properly with signed images
//
//	version 1.87(changed July 18, 2007)
//		changed Gaussian fitting to handle case of top axis instead of bottom axis, changed GaussianCenter(),
//		and also had to change DrawCross() and DrawBoxWithTicksAtXY()
//
//	version 1.88(changed Sept. 20, 2007)
//		changed FitOneGaussianPeak() so that it avoids errors when coefficients are really close to their limits of 1 or 0, but not quite (off by 1e-12)
//		FitOneGaussianPeak(), added optional Q argument to suppress printout
//
//	version 1.89(changed Sept. 25, 2007)
//		changed WinViewReadROI(), it now includes key in the wave note for when the region of the image read does not start at (0,0), callled startxRead & startyRead
//
//	version 1.90(changed Oct. 1, 2007)
//		replaced  DrawCross() and DrawBoxWithTicksAtXY(), with DrawMarker().  (I kept stubs to the old routines around for bakcward compatibility).
//
//	version 1.90(changed Oct. 23, 2007)
//		chnaged DrawMarker() to accept an optional color input
//		added ReDrawBoxes(), and changed ButtonProc_Boxes() to use ReDrawBoxes()
//
//	version 1.92(changed Nov. 17, 2007)
//		changed GraphImageStyle(),  added the gridRGB=(45000,45000,65535)
//		changed ReDrawBoxes(), made boxes half as big
//		changed DrawMarker() to use skinnier lines
//
//	version 1.93(changed Jan. 16, 2008)
//		changed DrawMarker,  added optional argument for line thickness, and optional window name
//
//	version 1.94(changed Nov. 6, 2008)
//		changed DrawMarker,  added optional argument for line thickness, and optional window name
//
//	version 1.95(changed Nov. 18, 2008)
//		added IMAGE_FILE_EXT
//
//	version 1.96(changed March 19, 2010)
//		removed DrawMarker(), just use the version in ImageDisplayScaling.ipf
//
//	version 1.97(changed Dec 15, 2011)
//		added a detectorID for Roper
//
//	version 1.98(changed Feb 19, 2012)
//		added WinView_MetaUnits
//
//	version 2.00(changed Nov 15, 2012)
//		changed the load functions to have the optional extras string
//
//	version 2.01(changed Nov 27, 2012)
//		changed FitOneGaussianPeak(), to take a wave ref, not a string.  Also removed FitOneGaussianPeak(), it is in ImageDisplayScaling.ipf
//		also remove getZrange(), from here, use the version in ImageDisplayScaling.ipf
//
//	version 2.02(changed Apr 30, 2014)
//		changed line terminations, CR -> LF
//
//	version 2.03(changed Jun 20, 2014)
//		added optional extras argument to WinViewReadHeader()
//
//	version 2.04(changed Apr 2, 2015)
//		use the GraphImageStyle() in ImageDisplayScaling.ipf, remove GraphImageStyle() from this file
//
//	version 2.05(changed Jun 24, 2015)
//		changed FitOneGaussianPeak(), it now uses the optional printIt variable.

StrConstant IMAGE_FILE_EXT = ".SPE"

// units for the metadata that usually accompanies an image
Strconstant WinView_MetaUnits ="X1:�m;Y1:�m;Z1:�m;H1:�m;F1:�m;X2:�m;Y2:�m;Z2:�m;H2:�m;F2:�m;exposure:s;"

Menu "WinView",  dynamic
	MenuItemIfWaveClassExists("Display image plot with buttons...","spe*","DIMS:2"), Graph_imageMake($"",0)
	help={"Display an image and put buttons on the plot for listing spots"}
	MarqueeWinViewMenuItem("Center a Gaussian and Add to pkList"), AddGaussianToList()
	MarqueeWinViewMenuItem("Report Center of a Gaussian"), GaussianCenter(1)
	help={"fit a gaussian to an ROI an image, and print results.  The region is selected with a marquee in an image plot"}
	MarqueeWinViewMenuItem("Report Center of Mass of a ROI"), GetCenterOfMass(NaN,NaN,NaN,NaN)
	help={"Compute the center of mass of a region on an image.  The region is selected with a marquee in an image plot"}
	"-"
	pkLIstWinViewMenuItem("Fit List of Gaussians..."), FitListOfGaussians("","pkList",-1)
	GcoefWinViewMenuItem("Reset pkList from Gaussians"), pkList = Reset_pkList_from_Gaussians(Gauss_coef)
	"-"
	"Load WinView File...", LoadWinViewFile("")
	help={"Load a WinView Image from the file, and then display the image"}
	MenuItemIfWaveClassExists("Get WinView Info","speImage*","DIMS:2"), WinViewInfo($"","")
	help={"Get the other information stored with the winview image"}
//	"Show WinView Procedures", DisplayProcedure/B=WinViewProcedures "CenterOfMass"
	     MenuItemIfWaveClassExists("   Find Z range of image","*","DIMS:2"), getZrange($"",NaN)
	help={"for an image, find the z range for the specified % range"}
End
//Menu "WinView"
//	"Display image plot with buttons...", Graph_imageMake($"",0)
//	help={"Display an image and put buttons on the plot for listing spots"}
//	"Center a Gaussian and Add to pkList", AddGaussianToList()
//	"Report Center of a Gaussian ", GaussianCenter(1)
//	help={"fit a gaussian to an ROI an image, and print results.  The region is selected with a marquee in an image plot"}
//	"Report Center of Mass of a ROI", GetCenterOfMass(NaN,NaN,NaN,NaN)
//	help={"Compute the center of mass of a region on an image.  The region is selected with a marquee in an image plot"}
//	"-"
//	"Fit List of Gaussians...", FitListOfGaussians("","pkList",-1)
//	"Reset pkList from Gaussians", pkList = Reset_pkList_from_Gaussians(Gauss_coef)
//	"-"
//	"Load WinView File...", LoadWinViewFile("")
//	help={"Load a WinView Image from the file, and then display the image"}
//	"Get WinView Info", WinViewInfo($"","")
//	help={"Get the other information stored with the winview image"}
////	"Show WinView Procedures", DisplayProcedure/B=WinViewProcedures "CenterOfMass"
//	"Choose Z range", getZrange($"",NaN)
//End
Menu "Load Waves"
	"Load WinView File...", LoadWinViewFile("")
	help={"Load a WinView Image from the file, and then display the image"}
End
Menu "Analysis"
	"Report Center of a Gaussian ", GaussianCenter(1)
	help={"fit a gaussian to an ROI an image, and print results.  The region is selected with a marquee in an image plot"}
	"Report Center of Mass of a ROI", GetCenterOfMass(NaN,NaN,NaN,NaN)
	help={"Compute the center of mass of a region on an image.  The region is selected with a marquee in an image plot"}
End
//
Function/S MarqueeWinViewMenuItem(item)
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
Function/S GcoefWinViewMenuItem(item)
	String item
	if (WaveExists(Gauss_coef))
		return item
	endif
	return "("+item
End
//
Function/S pkLIstWinViewMenuItem(item)
	String item
	if (strlen(StrVarOrDefault("pkList","")))
		return item
	endif
	return "("+item
End

Menu "GraphMarquee",dynamic
	"-"
	MarqueeWinViewMenuItem("Add Gaussian to pkList"),/Q, AddGaussianToList()
	help={"fit a gaussian to the Marquee, and add it to the peak list"}
	MarqueeWinViewMenuItem("Remove Peak From pkList"),/Q, RemovePeakFromList()
	help={"remove the peak in the marquee from the peak list"}
	MarqueeCursorWinViewMenuItem("Add Peak at Cursor A"),/Q, AddPeakAtCursor()
	help={"add position of cursor A (the round one) to the peak list"}
	MarqueeWinViewMenuItem("Report Center of Gaussian"),/Q, ReportCenterOfGaussian()
	help={"do a gaussian fit and report the result (does not add to peak list)"}
	MarqueeWinViewMenuItem("Report Center of Mass"),/Q, ReportCenterOfMass()
	help={"report the center of mass of the marquee (does not add to peak list)"}
End
//
Function/S MarqueeCursorWinViewMenuItem(item)
	String item
	if (strlen(CsrWave(A))<=0)				// These menu items requre that cursor A be on the image
		return "("+item
	endif
	String wName = StringFromList(0,ImageNameList("",";"))
	if (!strlen(wName))
		return "("+item
	endif
	if (!WaveExists(ImageNameToWaveRef("",wName)))
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


Static Strconstant ksControllerTypes="new120 (TYPE II);old120 (TYPE I);ST130;ST121;ST138;DC131 (Pentamax);ST133 (MicroMAX/SpectroMAX);ST135 (GPIB);VICCD;ST116 (GPIB);OMA3 (GPIB);OMA4"



// use the GraphImageStyle from ImageDisplayScaling.ipf
//
//Proc GraphImageStyle() : GraphStyle
//	PauseUpdate; Silent 1		// modifying window...
//	String wName = StringFromList(0,WaveList("*",";","DIMS:2,WIN:"))
//	if (strlen(wName)<1)
//		Abort "Unable to find image on top graph"
//	endif
//	ModifyGraph/Z width={Aspect,DimSize($wName,0)/DimSize($wName,1)}
//	ModifyGraph/Z grid=1, tick=2, mirror=1, minor=1, gridStyle=1
//	ModifyGraph/Z lowTrip=0.001, standoff=0, axOffset(bottom)=-1
//	ModifyGraph grid=1,gridRGB=(45000,45000,65535)
//	ControlInfo  button0
//	if (stringmatch(IgorInfo(2),"Macintosh" ))
//		ModifyGraph axOffset(left)=(V_flag==1) ? 0.9 : -1.1 		// mac
//	else
//		ModifyGraph axOffset(left)=(V_flag==1) ? -1.7 : -2.1 		// pc
//	endif
//	Variable/C lohi = getZrange($wName,5)
//	Variable lo = real(lohi), hi = imag(lohi)
//	if (!(WaveType($wName) & 0x40))
//		lo = abs(lo)
//		hi = abs(hi)
//		hi = max(lo,hi)
//		lo = -hi
//		ModifyImage $wName ctab= {lo,hi,RedWhiteBlue,0}
//	else
//		ModifyImage $wName ctab= {real(lohi),imag(lohi),Terrain,1}
//	endif
//
////	ModifyImage $wName ctab= {real(lohi),imag(lohi),Terrain,1}
////	ModifyImage $wName ctab= {real(lohi),imag(lohi),PlanetEarth,1}
////
////	ImageStats /M=1  $wName
////	ModifyImage $wName ctab= {*,V_max/10,PlanetEarth,1}
//	Label/Z left "y  (\\U)"
//	Label/Z bottom "x  (\\U)"
//	SetAxis/A/R left
//	ShowInfo
//EndMacro


Function AddGaussianToList()
	GetMarquee/Z
	if (V_flag==0)
		DoAlert 0, "no Marquee selected, do nothing"
		return 0
	endif
	GaussianCenter(0)
	AddPoint2pkList(K2,K4)

	String fldr = GetWavesDataFolder(ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";"))),1 )
	NVAR boxesOn = $(fldr+"boxesOn")
	if (!NVAR_Exists(boxesOn))
		Variable/G $(fldr+"boxesOn")
	endif
	if (boxesOn)
		DrawMarker(K2,K4,20,20,"BoxWithTicks")
	endif
End
Function/T RemovePeakFromList()
	Variable tol=5										// � tolerance in pixels
	Variable V_left, V_right, V_bottom, V_top
	String xaxis=StringByKey("XAXIS",ImageInfo("","",0)), yaxis=StringByKey("YAXIS",ImageInfo("","",0))
	GetMarquee/Z $yaxis, $xaxis
//	GetMarquee/Z left, bottom
	if (V_flag==0)
		DoAlert 0, "no Marquee selected, do nothing"
		return ""
	endif
	String fldr = GetWavesDataFolder(ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";"))),1 )
	SVAR pkList = $(fldr+"pkList")
	if (!SVAR_Exists(pkList))							// list empty, nothing to remove
		return ""
	endif

	Variable i,ix,iy,x0,y0
	x0 = (V_left+V_right)/2
	y0 = (V_bottom+V_top)/2
	i = ItemsInList(pkList)
	pkList = RemovePairFromList(pkList,x0,y0,tol)
	if (i == ItemsInList(pkList))
		DoAlert 0, "no peaks removed"
	endif
	return pkList
End
Function AddPeakAtCursor()
	if (strlen(CsrWave(A))<=0)
		DoAlert 0, "Cursor A not on the Image, do nothing"
		return 0
	endif
	AddPoint2pkList(pcsr(A),qcsr(A))
//	ButtonProc_Include("") 
	String fldr = GetWavesDataFolder(ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";"))),1 )
	NVAR boxesOn = $(fldr+"boxesOn")
	if (!NVAR_Exists(boxesOn))
		Variable/G $(fldr+"boxesOn")
	endif
	if (boxesOn)
		DrawMarker(pcsr(A),qcsr(A),20,20,"BoxWithTicks")
	endif
End
Function ReportCenterOfGaussian()
	GaussianCenter(1)
End
Function ReportCenterOfMass()
	GetCenterOfMass(NaN,NaN,NaN,NaN)
End




Function GaussianCenter(printON)	// returns 1 if error,  0 is OK
// the result is passed to calling function by W_coef or {K0, K1, ... K6}
// this fits a 2-d gaussina to a region on an image.  The region is selected with a marquee in an image plot, or specified
	Variable printON							// enables print out of result

	Variable V_left, V_right, V_bottom, V_top
	if (WinType("")!=1)
		Abort "Top window must be a graph, with region selected"
	endif
	String xaxis=StringByKey("XAXIS",ImageInfo("","",0)), yaxis=StringByKey("YAXIS",ImageInfo("","",0))
	GetMarquee/Z $yaxis, $xaxis
	if (V_flag==0)
		DoAlert 0, "No region selected on the graph"
		return 1
	endif
	if (V_right==0 && V_left==0 || V_top==0 && V_bottom==0)
		DoAlert 0, "Unable get range for horizonal or vertical axis info for the selected graph"
		return 1
	endif
//	GetMarquee/Z left, bottom
//	if (V_flag==0)
//		DoAlert 0, "No region selected on the graph"
//		return 1
//	endif
//	if (V_right==0 && V_left==0)
//		GetMarquee/Z left, top
//	endif
//	if (V_right==0 && V_left==0 || V_top==0 && V_bottom==0)
//		DoAlert 0, "Unable to horizonal axis info on the selected graph"
//		return 1
//	endif

//	Find the wave that is the image (assume that it is the only 2-d wave on the graph
	String ImageName = StringFromList(0,WaveList("*",";","DIMS:2,WIN:"))	// name of wave with image
	if (strlen(ImageName)<1)
		DoAlert 0, "Unable to find image on top graph"
		return 1
	endif
	Wave image = $ImageName

	// set up for the gaussian fit
	Variable i = FitOneGaussianPeak(Image,V_left,V_right,V_bottom,V_top,printIt=1)	// returns 1 if error,  0 is OK
	if (i)
		return i
	endif

	Wave W_coef=W_coef
	Wave W_sigma=W_sigma
	Variable x0,y0,dx,dy,cor
	Variable fw= ( 2*sqrt(2*ln(2)) )			// this factor corrects for 2-d sigma to FWHM, apply to K3, K5, and the corresponding errors
	x0 = W_coef[2]
	dx = W_coef[3] * fw
	y0 = W_coef[4]
	dy = W_coef[5] * fw
	cor = W_coef[6]
	if (printON)
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
		sprintf fmt, "Gaussian peak in '%%s' at (%%.%df�%%.%df, %%.%df�%%.%df) with FWHM of (%%.%df�%%.%df, %%.%df�%%.%df) and correlation of %%.%df�%%.%df  Amp=%%.3g\r",places1,places1,places1,places1,places2,places2,places2,places2,places3,places3
		printf fmt,ImageName,x0,W_sigma[2],y0,W_sigma[4],dx,W_sigma[3]*fw,dy,W_sigma[5]*fw,cor,W_sigma[6],W_coef[1]
//		printf "Gaussian peak in '%s' at (%.2f�%.2f, %.2f�%.2f) with FWHM of (%.2f�%.2f, %.2f�%.2f) and correlation of %.2f�%.2f\r",ImageName,x0,W_sigma[2],y0,W_sigma[4],dx,W_sigma[3]*fw,dy,W_sigma[5]*fw,cor,W_sigma[6]
		DoAlert 1, "Draw a cross at result"
		if (V_Flag==1)
			DrawMarker(x0,y0,abs(V_right-V_left),abs(V_top-V_bottom),"Cross")
		endif
	endif
	return 0
End


Function FitListOfGaussians(ImageName,pkList,halfWidth)		// build Gauss_coef and Gauss_sigma from pkList by fitting gaussians
	String ImageName						// name of image to use
	String pkList							// integer list of rough centers
	Variable halfWidth						// half width of box for fit range

	if (halfWidth<1 || strlen(pkList)<7 || DimSize($ImageName,1)<1)	// need to prompt
		halfWidth = (halfWidth<1) ? 10 : halfWidth
		Prompt halfWidth, "half width of box for fit region"
		Prompt pkList,"list of strings with lists of point", popup,StringList("*",";")
		Prompt ImageName, "name of image to use", popup, WaveList("*",";","DIMS:2")
		DoPrompt "Fit a list of Gaussian Peaks (results go into)", ImageName, pkList, halfWidth
		if (V_flag)
			return 1
		endif
		if (DimSize($ImageName,1)<1 || exists(pkList)!=2)
			return 1
		endif
		SVAR pk = $pkList
		pkList = pk
	endif

	Wave image = $ImageName
	if (!WaveExists(image))
		return 1
	elseif (halfWidth<1 || DimSize(image,1)<1 || strlen(pkList)<5)
		return 1
	endif

	Variable n=ItemsInList(pkList)
	Make/N=(7,n)/O Gauss_coef,Gauss_sigma
	Gauss_coef = NaN
	Gauss_sigma = NaN
	String item
	Variable i,x0,y0,xmin,xmax,ymin,ymax
	Variable nx=DimSize(image,0)
	Variable ny=DimSize(image,1)
	for (i=0;i<n;i+=1)
		item = StringFromList(i,pkList)
		x0 = str2num(StringFromList(0,item,","))
		y0 = str2num(StringFromList(1,item,","))
		xmin=max(0,x0-halfWidth)
		xmax=min(nx-1,x0+halfWidth)
		ymin=max(0,y0-halfWidth)
		ymax=min(ny-1,y0+halfWidth)
		if (FitOneGaussianPeak(image,xmin,xmax,ymin,ymax,printIt=1))
			x0=NaN
			y0=NaN
		else
			Wave W_coef=W_coef
			Wave W_sigma=W_sigma
			x0 =W_coef[2]
			y0 =W_coef[4]
			Gauss_coef[][i] = W_coef[p]
			Gauss_sigma[][i] = W_sigma[p]
		endif
		sprintf item, "%g,%g", x0,y0
	endfor
	return 0
End


//	Function FitOneGaussianPeak(ImageName,V_left,V_right,V_bottom,V_top,[Q])	// returns 1 if error,  0 is OK
//	// the result is passed to calling function by W_coef or {K0, K1, ... K6}
//	// this fits a 2-d gaussina to a region on an image.  The region is selected with a marquee in an image plot, or specified
//		String ImageName 							// name of wave with image
//		Variable V_left, V_right, V_bottom, V_top	// range to use for fit
//		Variable Q										// quiet flag
//	
//		Q = ParamIsDefault(Q) ? 0 : Q
//		Wave image = $ImageName
//		if (!WaveExists(image))
//			DoAlert 0, "Unable to find image named "+ImageName
//			return -1
//		endif
//	
//		Variable pixels=1					// using image in pixels
//		pixels = (DimDelta(image,0)==1 && DimDelta(image,1)==1 && DimOffset(image,0)==0 && DimOffset(image,1)==0)
//		if (pixels)
//			Variable swap
//			if (V_top<V_bottom)			// insist that bottom < top
//				swap = V_top
//				V_top = V_bottom
//				V_bottom = swap
//			endif
//			if (V_right<V_left)				// insist that left < right
//				swap = V_left
//				V_left = V_right
//				V_right = swap
//			endif
//			V_left = limit(floor(V_left),0,DimSize(image,0)-1)	// these are pixels, so they need to be integers
//			V_right = limit(ceil(V_right),0,DimSize(image,0)-1)
//			V_bottom = limit(floor(V_bottom),0,DimSize(image,1)-1)
//			V_top = limit(ceil(V_top),0,DimSize(image,1)-1)
//		endif
//	
//		// set up for the gaussian fit
//		Make/N=7/O/D W_coef, W_sigma=W_sigma
//		Variable V_FitOptions=4, V_FitError=0,V_FitQuitReason=0
//		V_FitOptions=4
//		//	Execute "SetIgorOption UseVeclib=0"
//		Make/O/T T_Constraints={"K3 > 0","K5 > 0","K6 > 0","K6 < 1"}
//		if (pixels)
//			CurveFit/Q/N Gauss2D image[V_left,V_right][V_bottom,V_top]/C=T_Constraints 
//		else
//			CurveFit/Q/N Gauss2D image(V_left,V_right)(V_bottom,V_top)/C=T_Constraints 
//		endif
//		KillWaves/Z T_Constraints
//	//	CurveFit/Q/N Gauss2D image[V_left,V_right][V_bottom,V_top]		// constraints on K6 are not needed for Igor 5
//	// however, if constraints are not used, then K6 does not remain in [0,1]
//	//	if (V_FitError && WhichListItem("reFit_GaussianPkList",GetRTStackInfo(0))<0)
//	//	if (V_FitError && WhichListItem("reFit_GaussianPkList",GetRTStackInfo(0))<0 && WhichListItem("NewFitPeaks",GetRTStackInfo(0))<0)
//		if (V_FitError && WhichListItem("reFit_GaussianPkList",GetRTStackInfo(0))<0 && WhichListItem("NewFitPeaks",GetRTStackInfo(0))<0 && WhichListItem("FitPeaksStepWise",GetRTStackInfo(0))<0)
//			String errString = ""
//			if (V_FitError)
//				errString += "V_FitError = '"
//				errString += SelectString(V_FitError-2,"Singular Matrix","Out of memory","Function return NaN or INF")
//				errString += "'     "
//			endif
//			if (V_FitQuitReason)
//				errString += "V_FitQuitReason = '"
//				errString += SelectString(V_FitQuitReason-2,"iteration limit was reached","user stopped the fit","limit of passes without decreasing chi-square")
//				errString += "'"
//			endif
//			if (!Q)
//				printf "%s\r",errString
//			endif
//		endif
//		if (!V_FitError)
//			Wave W_coef = W_coef
//			Variable zero = -1e-15, one=1+1e-12					// almost zero, almost one
//			V_FitError += (W_coef[3]<zero || W_coef[5]<zero || W_coef[6]<zero || W_coef[6]>one) ? 32 : 0			//  constraints failed, set bit 5
//		endif
//		return V_FitError
//	End


Function AddPoint2pkList(ix,iy)
	Variable ix,iy			// pixel coordinates of point to add

	String fldr = GetWavesDataFolder(ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";"))),1 )
	SVAR pkList = $(fldr+"pkList")
	if (!SVAR_Exists(pkList))
		String/G $(fldr+"pkList")
	endif
	if (numtype(ix) || numtype(iy))
		return 0
	endif
	pkList = AddListItem(num2istr(ix)+","+num2istr(iy),pkList,";",0)	// add to front of list
	pkList = RemoveDuplicatePairsFromList(pkList)			// no duplicates
	pkList = SortList(pklist,";",2)
End


Function/T Reset_pkList_from_Gaussians(Gcoef)
	Wave Gcoef
	if (DimSize(Gcoef,0)!=7)
		return ""
	endif
	Variable i,n=DimSize(Gcoef,1)
	String list=""
	for (i=0;i<n;i+=1)
		list = AddListItem(num2istr(Gcoef[2][i])+","+num2istr(Gcoef[4][i]), list,";",i)
	endfor
	return list
end



Function/C GetCenterOfMass(V_left,V_right,V_bottom,V_top)
// this gets the center of mass of a region on an image.  The region is selected with a marquee in an image plot
	Variable V_left, V_right, V_bottom, V_top

	Variable printON=0					// enables print out of result
	Variable pixels=1					// using image in pixels
	if (WinType("")!=1)
		Abort "Top window must be a graph, with region selected"
	endif
	if (!(V_left<V_right && V_bottom<V_top))
		String xaxis=StringByKey("XAXIS",ImageInfo("","",0)), yaxis=StringByKey("YAXIS",ImageInfo("","",0))
		GetMarquee/Z $yaxis, $xaxis
//		GetMarquee/Z left, bottom
		if (V_flag==0)
			DoAlert 0, "No region selected on the graph"
			return cmplx(0,0)
		endif
		printON = 1
	endif

//	Find the wave that is the image (assume that it is the only 2-d wave on the graph
	String ImageName = StringFromList(0,WaveList("*",";","DIMS:2,WIN:"))	// name of wave with image
	if (strlen(ImageName)<1)
		DoAlert 0, "Unable to find image on top graph"
		return 1
	endif
	Wave image = $ImageName

	if (pixels)
		Variable swap
		if (V_top<V_bottom)			// insist that bottom < top
			swap = V_top
			V_top = V_bottom
			V_bottom = swap
		endif
		if (V_right<V_left)				// insist that left < right
			swap = V_left
			V_left = V_right
			V_right = swap
		endif
		V_left = limit(floor(V_left),0,DimSize(image,0)-1)	// these are pixels, so they need to be integers
		V_right = limit(ceil(V_right),0,DimSize(image,0)-1)
		V_bottom = limit(floor(V_bottom),0,DimSize(image,1)-1)
		V_top = limit(ceil(V_top),0,DimSize(image,1)-1)
	endif

	Variable nx = abs(V_right - V_left) + 1
	Variable ny = abs(V_top - V_bottom) + 1

	Make/O/N=(nx) temp_sum_region_x_
	SetScale/P x V_left,1,"", temp_sum_region_x_
	Make/O/N=(ny) temp_sum_region_y_
	SetScale/P x V_bottom,1,"", temp_sum_region_y_
	Variable total, xMoment, yMoment
	Variable xx, ix, yy, iy

	// find the background,  for a constant background, just sum up the edge of the frame
	Variable bkg=0, useBkg=1
	if (useBkg)
		yy = pnt2x(temp_sum_region_y_,0)
		temp_sum_region_x_ = image(x)(yy)  ;  temp_sum_region_x_[0] = 0
		bkg = sum(temp_sum_region_x_,-inf,inf)

		yy = pnt2x(temp_sum_region_y_,ny-1)
		temp_sum_region_x_ = image(x)(yy)  ;  temp_sum_region_x_[0] = 0
		bkg += sum(temp_sum_region_x_,-inf,inf)

		xx = pnt2x(temp_sum_region_x_,0)
		temp_sum_region_y_ = image(xx)(y)  ;  temp_sum_region_y_[0] = 0
		bkg += sum(temp_sum_region_y_,-inf,inf)

		xx = pnt2x(temp_sum_region_x_,nx-1)
		temp_sum_region_y_ = image(xx)(y)  ;  temp_sum_region_y_[0] = 0
		bkg += sum(temp_sum_region_y_,-inf,inf)

		bkg /= 2*(nx+ny-2)						// divide by number of points used, bkg per point
	endif

	temp_sum_region_x_ = 0						// collapse onto x-axis
	iy = 0
	do
		yy = pnt2x(temp_sum_region_y_,iy)
		temp_sum_region_x_ += image(x)(yy) - bkg
		iy += 1
	while (iy<ny)

	total = sum(temp_sum_region_x_,-inf,inf)		// calculate the first moment
	temp_sum_region_x_ *= x
	xMoment = sum(temp_sum_region_x_,-inf,inf)/total

	temp_sum_region_y_ = 0						// collapse onto y-axis
	ix = 0
	do
		xx = pnt2x(temp_sum_region_x_,ix)
		temp_sum_region_y_ += image(xx)(x) - bkg
		ix += 1
	while (ix<nx)
	total = sum(temp_sum_region_y_,-inf,inf)		// calculate the first moment
	temp_sum_region_y_ *= x
	yMoment = sum(temp_sum_region_y_,-inf,inf)/total

	if (printON)
//		print "bkg = ",bkg
		printf "Center of Mass of  '%s'(%g,%g)(%g,%g) = (%g,%g)\r",ImageName,V_left, V_right, V_bottom, V_top, xMoment,yMoment
		DoAlert 1, "Draw a cross at result"
		if (V_Flag==1)
			DrawMarker(xMoment,yMoment,V_right-V_left,V_top-V_bottom,"Cross")
		endif
	endif

	KillWaves/Z temp_sum_region_x_, temp_sum_region_y_
	return cmplx(xMoment,yMoment)
End




//Static Function/C getZrange(image,delta[plane,lo,hi])
//	Wave image
//	Variable delta				//=5,  clip amount (in %), not fraction
//	Variable plane				// This is not implemented here
//	Variable lo,hi				// This is not implemented here
//
//	Variable printIt=0
//	if (!WaveExists(image) || numtype(delta) || delta<=0)
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
//	if (numtype(delta) || delta<0 || delta > 99)
//		Abort "delta ="+num2str(delta)+", which is too small"
//	endif
//	delta /= 100
//
//	ImageStats/M=1 image					// check if no range in image (avoids Histogram error)
//	if (V_max==V_min)
//		return cmplx(V_min,V_max)
//	endif
//
//	Variable N = numpnts(image)/2
//	Make/N=(N)/O hh_temp_
//	SetScale/P x 0,1,"", hh_temp_
//	Histogram/B=1 image,hh_temp_
//	Wave hist=hh_temp_
//	Integrate hist
//
//	Variable hmin = hist[0]
//	Variable dh = hist[N-1]-hmin
//	lo = BinarySearchInterp(hist, delta*dh+hmin)
//	hi = BinarySearchInterp(hist, (1-delta)*dh+hmin)
//	lo = DimOffset(hist,0) + lo *DimDelta(hist,0)
//	hi = DimOffset(hist,0) + hi *DimDelta(hist,0)
//	if (printIt)
//		printf "the [%g%%, %g%%] range = [%g, %g]\r",delta*100,100*(1-delta),lo,hi
//	endif
//	KillWaves/Z hist
//	return cmplx(lo,hi)
//End


// =======================================================================


//	Function/S test_WinViewSaveImage(redim)
//		Variable redim										// if true, do a redimension
//		PathInfo/S home
//		String fName = S_path+"xtest.spe"					// fully qualified name of file to open
//		printf "writing the test file: '%s'\r",fName
//	
//		Wave image1 = $(LoadWinViewFile("Macintosh HD:Users:tischler:data:Sector 34:Si Apr06:BA1a1746_2788:WhWi_1746.SPE"))
//		String str,listIn = note(image1)
//		print ""
//		print "initial list:"
//		print listIn
//		print ""
//	
//		if (redim)
//			Redimension/I image1
//		endif
//	
//		str = WinViewSaveImage(fName,image1)
//		if (strlen(str)<1)
//			print "WinViewSaveImage() failed, halting"
//			return ""
//		endif
//	
//		Wave image2 = $(LoadWinViewFile(fName))		// then re-read the spe file to check it
//		String listOut = note(image2)
//	
//		if (stringmatch(listIn,listOut))
//			print "The header strings match, before and after, SUCCESS"
//		else
//			beep
//			print "The header strings do not match, FAILURE"
//			print "final list:"
//			print listOut
//			print ""
//		endif
//	
//		Variable i,N=numpnts(image1)
//		for (i=0;i<N;i+=1)
//			if (image1[i]!=image2[i])
//				break
//			endif
//		endfor
//		if (i<N)
//			beep
//			print "images differ at i =",i
//		else
//			print "images are identical"
//		endif
//		return fName										// returns the name of the file writen
//	End
Function/T WinViewSaveImage(fileName,image)
	String fileName					// fully qualified name of file to open (will not prompt)
	Wave image
	if (!WaveExists(image))
		print "ERROR, WinViewSaveImage(), the image wave does not exist"
		return ""
	elseif (WaveDims(image)!=2)
		printf "ERROR, WinViewSaveImage(), only writes 2-d waves, '%s' is %d-d\r",WaveDims(image),NameOfWave(image)
		return ""
	endif
	Variable WinViewType = igorType2WinView(WaveType(image))
	if (numtype(WinViewType))
		printf "ERROR, WinViewSaveImage(), cannot write spe files for WaveType(%s) = 0x%x\r",NameOfWave(image),WaveType(image)
		return ""
	endif

	String wnote=note(image)
	wnote = ReplaceNumberByKey("numType",wnote,igorType2WinView(WaveType(image)),"=")
	wnote = ReplaceNumberByKey("xdim",wnote,DimSize(image,0),"=")
	wnote = ReplaceNumberByKey("ydim",wnote,DimSize(image,1),"=")

	STRUCT speHeader header				// structure containing all of header
	WinViewList2Header(wnote,header)		// convert list of values (from wnote) into a header
	fileName = WinViewWriteHeader(fileName,header)	// write header to file

	Variable fid								// file id
	Open/T="????"/A/Z=1 fid as fileName	// this acutally opens file
	fileName = S_fileName
	if (V_flag)
		return ""
	endif
	FSetPos fid, 4100
	FBinWrite/B=3 fid,image				// writes image
	Close fid
	return fileName
End

Function/S LoadWinViewFile(fName,[extras])
	String fName											// fully qualified name of file to open
	String extras											// not used here for anything yet
	extras = SelectString(ParamIsDefault(extras),extras,"")

	Variable refNum
	if (strlen((OnlyWinFileName(fName)))<1)				// call dialog if no file name passed
		Open /D/M=".spe file"/R/T="????" refNum		// use /D to get full path name
		fName = S_filename
	endif
	if (strlen(fName)<1)									// no file name, quit
		return ""
	endif

	String wName = WinViewReadROI(fName,0,-1,0,-1)	// load file into wName
	if (strlen(wName)<1)
		return ""
	endif
	Wave image = $wName
	if (ItemsInList(GetRTStackInfo(0))<=1)
		String wnote = note(image)
		Variable xdim=NumberByKey("xdim", wnote,"=")
		Variable ydim=NumberByKey("ydim", wnote,"=")
		String bkgFile = StringByKey("bkgFile", wnote,"=")
		printf "for file '"+fName+"'"
		if (strlen(bkgFile)>0)
			printf ",             background file = '%s'",  bkgFile
		endif
		printf "\r"
		printf "total length = %d x %d  = %d points\r", xdim,ydim,xdim*ydim
		print "number type is  '"+WinViewFileTypeString(NumberByKey("numType", wnote,"="))+"'"
		print "Created a 2-d wave    '"+wName+"'"
		DoAlert 1, "Display this image"
		if (V_Flag==1)
			Graph_imageMake(image,NaN)
		endif
	endif
	return GetWavesDataFolder(image,2)
End
//Function/S LoadWinViewFile(S_filename)
//	String S_filename											// fully qualified name of file to open
//
//	Variable refNum
//	if (strlen((OnlyWinFileName(S_filename)))<1)			// call dialog if no file name passed
//		Open /D/M=".spe file"/R/T="????" refNum			// use /D to get full path name
//	endif
//	if (strlen(S_filename)<1)									// no file name, quit
//		return ""
//	endif
//
//	String wName = WinViewReadROI(S_filename,0,-1,0,-1)	// load file into wName
//	Wave image = $wName
//	String wnote = note(image)
//	String bkgFile = StringByKey("bkgFile", wnote,"=")
//	Variable xdim=NumberByKey("xdim", wnote,"=")
//	Variable ydim=NumberByKey("ydim", wnote,"=")
//	if (ItemsInList(GetRTStackInfo(0))<=1)
//		printf "for file '"+S_filename+"'"
//		if (strlen(bkgFile)>0)
//			printf ",             background file = '%s'",  bkgFile
//		endif
//		printf "\r"
//		printf "total length = %d x %d  = %d points\r", xdim,ydim,xdim*ydim
//		print "number type is  '"+WinViewFileTypeString(NumberByKey("numType", wnote,"="))+"'"
//		print "Created a 2-d wave    '"+wName+"'"
//	endif
//	DoAlert 1, "Display this image"
//	if (V_Flag!=1)
//		return GetWavesDataFolder(image,2)
//	endif
//	return Graph_imageMake(image,NaN)
//End


//Function test()
//	Variable width
//	Variable n0, n1
//	String baseName=""											// fully qualified name of file to open
//
//	String S_filename=""										// fully qualified name of file to open
//	Variable refNum
//	if (strlen((OnlyWinFileName(S_filename)))<1)			// call dialog if no file name passed
//		Open /D/M=".spe file"/R/T="????" refNum			// use /D to get full path name
//	endif
//	if (strlen(S_filename)<1)									// no file name, quit
//		return 1
//	endif
//
//	String wname
////	wname = WinViewReadROI(S_filename,547,640,497,600)
//	wname = WinViewReadROI(S_filename,0,-1,0,-1)
//	print "wname=",wname
//End
//
//		or just use:
//		print WinViewReadROI("Macintosh HD:Users:tischler:GeNoise:Ge_inten_1.SPE",0,-1,0,-1)
//
Function/T WinViewReadROI(fileName,i0,i1,j0,j1,[extras])
	String fileName					// fully qualified name of file to open (will not prompt)
	Variable i0,i1,j0,j1			// pixel range of ROI (if i1 or j1<0 then use whole image)
	String extras											// not used here for anything yet
	extras = SelectString(ParamIsDefault(extras),extras,"")

	String wnote=""
	if (strlen(extras))							// wave note to add to file read in
		wnote = WinViewReadHeader(fileName,extras=extras)
	else	
		wnote = WinViewReadHeader(fileName)
	endif
	if (!strlen(wnote))
		return ""
	endif
	Variable xdim=NumberByKey("xdim", wnote,"=")
	Variable ydim=NumberByKey("ydim", wnote,"=")
	Variable itype=NumberByKey("numType", wnote,"=")
	i1 = (i1<1) ? xdim-1 : i1								// -1 flags use whole range
	j1 = (j1<1) ? ydim-1 : j1

	String inWaveName
	if (ParamIsDefault(extras))
		inWaveName = WinViewLoadROI(fileName,itype,xdim,i0,i1,j0,j1)	// name of ROI read in
	else
		inWaveName = WinViewLoadROI(fileName,itype,xdim,i0,i1,j0,j1,extras=extras)	// name of ROI read in
	endif
	if (strlen(inWaveName)<1)
		return ""											// nothing read in
	endif

//	String wName = OnlyFileName(fileName)				// rename image based on the file name
	String separator = SelectString(stringmatch(IgorInfo(2),"Macintosh"),"\\",":")
	String wName = ParseFilePath(3, fileName,separator,0,0)
	wName = CleanupName(wName,0)						// wave name based on the file name

	if (exists(wName))										// if wave already exists, create unique name
		wName = wName+"_"
		wName = UniqueName(wName, 1, 1)
	endif
	Rename $inWaveName $wName
	Wave image = $wName
	SetScale/P x 0,1,"pixel", image
	SetScale/P y 0,1,"pixel", image
	if (i0!=0 || j0!=0)										// sub-region of image does not start at (0,0)
		wnote = ReplaceNumberByKey("startxRead",wnote,i0,"=")
		wnote = ReplaceNumberByKey("startyRead",wnote,j0,"=")
	endif
	wnote = ReplaceStringByKey("waveClass",wnote,"speImage","=")
	Note/K image, wnote
	return GetWavesDataFolder(image,2)
End



//Function/S test_WinViewReadHeader(fName)
//	String fName											// fully qualified name of file to open
//	Variable refNum
//	if (strlen((OnlyWinFileName(fName)))<1)				// call dialog if no file name passed
//		Open /D/M=".spe file"/R/T=".SPE" refNum			// use /D to get full path name
//		fName = S_filename
//	endif
//	if (strlen(fName)<1)									// no file name, quit
//		return ""
//	endif
//	String wnote =  WinViewReadHeader(fName)
//	if (ItemsInList(GetRTStackInfo(0))<2)
//		print wnote
//	endif
//	return wnote
//End
//
Function/T WinViewReadHeader(fileName,[extras])
	String fileName					// fully qualified name of file to open (will not prompt)
	String extras					// optional switches (only supports EscanOnly in this routine)
	extras = SelectString(ParamIsDefault(extras),extras,"")

	Variable fid												// file id (file is assumed already opened)
	Open /Z/M=".spe file"/R/T="????" fid as fileName		// this acutally opens file
	if (V_flag)
		return ""											// could not open file
	endif
	fileName = S_fileName
	FStatus fid
	if (V_logEOF<4100)
		Close fid											// file too short to be interpreted
		return ""
	endif

	NewDataFolder/O root:Packages
	Make/N=2104/U/B/O root:Packages:WinViewHeaderExtra	// additional 4100-1996=2104 bytes, to reach start of image
	Wave buffer=root:Packages:WinViewHeaderExtra

	STRUCT speHeader header								// structure containing all of header
	FBinRead/B=3 fid,header								// read header in one gulp
	FBinRead/B=0 fid, buffer
	Close fid
	// printf "in header, version = '%s'\r",header.version
	String/G root:Packages:WinViewHeaderStruct
	SVAR WinViewHeaderStruct = root:Packages:WinViewHeaderStruct
	StructPut/S/B=2 header, WinViewHeaderStruct

	String wnote=""											// wave note to add to file read in
	wnote = ReplaceStringByKey("imageFileName", wnote, ParseFilePath(0,fileName,":",1,0),"=")
	wnote = ReplaceStringByKey("imageFilePath", wnote, ParseFilePath(1,fileName,":",1,0),"=")
	String flatFile="", bkgFile=""
	flatFile = header.flatFile1
	if (strlen(flatFile)>=60)
		flatFile += header.flatFile2
	endif

	bkgFile = header.bkgFile1
	if (strlen(bkgFile)>=60)
		bkgFile += header.bkgFile2
	endif

	String str = UniqueName("userStrings",1,0)			// change the header.userStrj into strings
	Make/T/N=5/O $str
	Wave/T userString = $str
	userString[0] = header.userStr0
	userString[1] = header.userStr1
	userString[2] = header.userStr2
	userString[3] = header.userStr3
	userString[4] = header.userStr4
	if (header.bkgApplied)
		wnote = ReplaceNumberByKey("bkgApplied",wnote,1,"=")
		bkgFile = OnlyWinFileName(bkgFile)
		if (strlen(bkgFile)>0)
			wnote = ReplaceStringByKey("bkgFile",wnote,bkgFile,"=")
		endif
	endif

	wnote = ReplaceStringByKey("num_Type",wnote,WinViewFileTypeString(header.datatype),"=")
	wnote = ReplaceNumberByKey("numType",wnote,header.datatype,"=")
	wnote = ReplaceNumberByKey("xdim",wnote,header.xdim,"=")			// x size of image in file
	wnote = ReplaceNumberByKey("ydim",wnote,header.ydim,"=")			// y size of image in file
	wnote = ReplaceNumberByKey("xDimDet",wnote,header.xDimDet,"=")	// detector x dimension of chip
	wnote = ReplaceNumberByKey("yDimDet",wnote,header.yDimDet,"=")	// detector y dimension of chip
	if (strlen(header.edate))
		wnote = ReplaceStringByKey("dateExposed",wnote,header.edate,"=")	// date that exposure was taken
	endif
	if (header.ehour!=0 || header.eminute!=0)
		str = num2istr(header.ehour)+":"+num2istr(header.eminute)
		wnote = ReplaceStringByKey("timeExposed",wnote,str,"=")			// time that exposure was taken
	endif

	if (header.flatFieldApplied)												// a flat field was applied
		wnote = ReplaceNumberByKey("flatFieldApplied",wnote,1,"=")
		flatFile = OnlyWinFileName(flatFile)
		if (strlen(flatFile)>0)
			wnote = ReplaceStringByKey("flatFile",wnote,bkgFile,"=")
		endif
	endif
	wnote = ReplaceNumberByKey("exposure",wnote,header.exposure,"=")	// exposure duration
	wnote = ReplaceNumberByKey("ADCrate",wnote,header.ADCrate,"=")
	wnote = ReplaceNumberByKey("ADCtype",wnote,header.ADCtype,"=")
	wnote = ReplaceNumberByKey("ADCresolution",wnote,header.ADCresolution,"=")
	wnote = ReplaceNumberByKey("geo_rotate",wnote,!!(header.geometric&1),"=")	// geometric 1==rotate, 2==reverse, 4==flip
	wnote = ReplaceNumberByKey("geo_reverse",wnote,!!(header.geometric&2),"=")
	wnote = ReplaceNumberByKey("geo_flip",wnote,!!(header.geometric&4),"=")

	Variable i,j
	i = header.controllerType	-1											// controller type
	wnote = AddListItem("controllerType="+StringFromList(i,ksControllerTypes),wnote)
	if (i==6)
		wnote = AddListItem("detectorID="+"Roper ST135 (GPIB)",wnote)
	endif

	header.NumROI = max(header.NumROI,1)								// zero means one
	if (header.NumROI>1)
		Abort "this file contains multiple ROI's"
		wnote = ReplaceNumberByKey("NumROI",wnote,header.NumROI,"=")
	endif
	if (header.NumROI<=1)
		wnote = ReplaceNumberByKey("startx",wnote,header.ROIinfo[0].startx,"=")
		wnote = ReplaceNumberByKey("endx",wnote,header.ROIinfo[0].endx,"=")
		wnote = ReplaceNumberByKey("groupx",wnote,header.ROIinfo[0].groupx,"=")
		wnote = ReplaceNumberByKey("starty",wnote,header.ROIinfo[0].starty,"=")
		wnote = ReplaceNumberByKey("endy",wnote,header.ROIinfo[0].endy,"=")
		wnote = ReplaceNumberByKey("groupy",wnote,header.ROIinfo[0].groupy,"=")
	else
		for (i=0; i<header.NumROI; i+=1)
			str="_"+num2istr(i+1)
			wnote = ReplaceNumberByKey("startx"+str,wnote,header.ROIinfo[i].startx,"=")
			wnote = ReplaceNumberByKey("endx"+str,wnote,header.ROIinfo[i].endx,"=")
			wnote = ReplaceNumberByKey("groupx"+str,wnote,header.ROIinfo[i].groupx,"=")
			wnote = ReplaceNumberByKey("starty"+str,wnote,header.ROIinfo[i].starty,"=")
			wnote = ReplaceNumberByKey("endy"+str,wnote,header.ROIinfo[i].endy,"=")
			wnote = ReplaceNumberByKey("groupy"+str,wnote,header.ROIinfo[i].groupy,"=")
		endfor
	endif

	String item
	for (j=0;j<5;j+=1)					// go through each userString[j] and add to note
		str = userString[j]
		i = strsearch(str,"VAL ",0)			// strip off the PV name, yum:userStringCalc1.SVAL, or 34ide:userStringCalc1.SVAL
		if (!strsearch(userString[0],"yum:",0) && i>0)// special for microdiffraction on 33ID-E
			str = str[i+4,Inf]
		elseif (!strsearch(userString[0],"34ide:",0) && i>0)// special for microdiffraction on 33ID-E
			str = str[i+4,Inf]
		endif
		for (item=StringFromList(0,str),i=0; strlen(item); i+=1, item=StringFromList(i,str))
				wnote = AddListItem(item,wnote)	// add each semicolon separated part of userString[j] to note
		endfor
	endfor
	KillWaves/Z userString

	// this section is special for Sample & Wire positioner's,  if Y1 & Z1 present, then add H1 & F1 (ditto for Y2 & Z2)
	Variable Y1,Z1,H1,F1, Y2,Z2,H2,F2
	Y1 = NumberByKey("Y1",wnote,"=")
	Z1 = NumberByKey("Z1",wnote,"=")
	H1 = NumberByKey("H1",wnote,"=")
	F1 = NumberByKey("F1",wnote,"=")
	Y2 = NumberByKey("Y2",wnote,"=")
	Z2 = NumberByKey("Z2",wnote,"=")
	H2 = NumberByKey("H2",wnote,"=")
	F2 = NumberByKey("F21",wnote,"=")
	Variable cosTheta = NumVarOrDefault("root:Packages:geometry:cosThetaWire",cos(PI/4))
	Variable sinTheta = NumVarOrDefault("root:Packages:geometry:sinThetaWire",sin(PI/4))
	if (numtype(Y1+Z1)==0 && numtype(H1)==2 && numtype(F1)==2)	// Y1&Z1 exist, H1&F1 do not
		H1 =  Y1*sinTheta + Z1*cosTheta
		F1 = -Y1*cosTheta + Z1*sinTheta
		wnote = ReplaceNumberByKey("H1",wnote,H1,"=")
		wnote = ReplaceNumberByKey("F1",wnote,F1,"=")
	endif
	if (numtype(Y2+Z2)==0 && numtype(H2)==2 && numtype(F2)==2)	// Y2&Z2 exist, H2&F2 do not
		H2 =  Y2*sinTheta + Z2*cosTheta
		F2 = -Y2*cosTheta + Z2*sinTheta
		wnote = ReplaceNumberByKey("H2",wnote,H2,"=")
		wnote = ReplaceNumberByKey("F2",wnote,F2,"=")
	endif
//
//	printf "for file '"+fileName+"'"
//	if (header.bkgApplied && (strlen(bkgFile)>0))
//		printf ",             background file = '%s'\r",  bkgFile
//	else
//		printf "\r"
//	endif
//	printf "xdim = %d      ydim = %d      ", header.xdim,header.ydim
//	i = header.xdim * header.ydim
//	printf "total length = %d x %d  = %d points\r", header.xdim,header.ydim,i
//	print "number type is  '"+WinViewFileTypeString(header.datatype)+"'"
	return wnote
End
//
//
//Function/S test_WinViewWriteHeader()
//	PathInfo/S home
//	String fName = S_path+"xtest.spe"					// fully qualified name of file to open
//	printf "writing the test file: '%s'\r",fName
//	String listIn = test_WinViewReadHeader("Macintosh HD:Users:tischler:data:Sector 34:Si Apr06:BA1a1746_2788:WhWi_1746.SPE")
//	print ""
//	print "initial list:"
//	print listIn
//	print ""
//
//	STRUCT speHeader header							// structure containing all of header
//	WinViewList2Header(listIn,header)					// convert list of values (from wnote) into a header
//	WinViewWriteHeader(fName,header)
//
//	// then re-read the header to check
//	String listOut = test_WinViewReadHeader("Macintosh HD:Users:tischler:data:Sector 34:Si Apr06:BA1a1746_2788:WhWi_1746.SPE")
//
//	if (stringmatch(listIn,listOut))
//		print "The header strings match, before and after, SUCCESS"
//	else
//		print "The header strings do not match, FAILURE"
//		print "final list:"
//		print listOut
//		print ""
//	endif
//	return listOut										// returns the key word pair list
//End
//

Function/T WinViewWriteHeader(fileName,header)
	String fileName					// fully qualified name of file to open (will not prompt)
	STRUCT speHeader &header							// structure containing all of header

	Variable fid											// file id (file is assumed already opened)
	Open/M="new .spe file"/T="????" fid as fileName	// this acutally opens file
	fileName = S_fileName
	if (strlen(fileName)<1)							// no file opened
		return ""
	endif
	FBinWrite/B=3 fid,header							// writes structure part of header

	if (exists("root:PackagesWinViewHeaderExtra"))
		Wave buffer=root:PackagesWinViewHeaderExtra
		FBinWrite/B=0 fid,buffer						// writes 1100 more bytes of zero, to fill to 4100
	else
		String wName = UniqueName("buffer",1,0)		// header is 1996 long, image starts at byte 4100
		Make/N=2104/U/B $wName					// so need an additional 2104
		Wave buffer=$wName
		buffer = 0
		FBinWrite/B=0 fid,buffer						// writes 1100 more bytes of zero, to fill to 4100
		KillWaves/Z buffer
	endif
	Close fid
	return fileName
End
//Function/T WinViewWriteHeader(fileName,header)
//	String fileName					// fully qualified name of file to open (will not prompt)
//	STRUCT speHeader &header							// structure containing all of header
//
//	String wName = UniqueName("buffer",1,0)			// header is 1872 long, image starts at byte 4100
//	Make/N=2228/U/B $wName						// so need an additional 4100-1872=2228 bytes
//	Wave buffer=$wName
//	buffer = 0
//
//	Variable fid											// file id (file is assumed already opened)
//	Open/M="new .spe file"/T="????" fid as fileName	// this acutally opens file
//	fileName = S_fileName
//	if (strlen(fileName)<1)							// no file opened
//		return ""
//	endif
//	FBinWrite/B=3 fid,header							// writes structure part of header
//	FBinWrite/B=3 fid,buffer							// writes 2228 more bytes of zero, to fill to 4100
//	Close fid
//
//	KillWaves/Z buffer
//	return fileName
//End
//
//
Structure speHeader		// version 1.6 spe header
	uint16	dioden					//     0
	int16	avgexp					//     2
	int16	msecExpose				//     4
	uint16	xDimDet				//     6, Detector x dimension of chip
	int16	mode					//     8
	float	exposure				//   10, exposure in seconds
	int16	asyavg					//   14
	int16	asyseq					//   16
	uint16	yDimDet				//   18, Detector y dimension of chip
	char	edate[10]				//   20, DDMonYYYY, null terminated
	uint16	ehour					//   30, experiment time hours
	uint16	eminute					//   32, experiment time minutes
	int16	noscan					//   34
	int16	fastacc					//   36
	int16	esecond					//   38, experiment time seconds
	int16	DetType					//   40
	uint16	xdim					//   42, actual # of pixels on x axis
	int16	stdiode					//   44
	float	nanox					//   46
	float	calibdio[10]			//   50
	char	fastfile[16]				//   90
	int16	asynen					// 106
	uint16	datatype				// 108, image data type, 0=float, 1=long, 2=int, 3=uint
	float	calibnan[10]			// 110
	uint16	bkgApplied				// 150, true if a background file was substracted
	int16	astdiode					// 152
	uint16	minblk					// 154
	uint16	numinblk				// 156
	double	calibpol	[4]				// 158
	uint16	ADCrate					// 190, ADC rate
	uint16	ADCtype				// 192, ADC type
	uint16	ADCresolution			// 194, ADC resolution
	uint16	ADCbitAdjust			// 196
	uint16	gain					// 198
	char	userStr0[80]			// 200, experiment remarks
	char	userStr1[80]			// 280
	char	userStr2[80]			// 380
	char	userStr3[80]			// 440
	char	userStr4[80]			// 520
	uint16	geometric				// 600, 1==rotate, 2==reverse, 4==flip
	char	xlabel[16]				// 602
	uint16	cleans					// 618
	uint16	NumSkpPerCln			// 620
	char	califile[16]				// 622
	char	bkgdfile[16]			// 638
	int16	srccmp					// 654
	uint16	ydim					// 656, actual # of pixels on y axis
	int16	scramble				// 658
	int32	lexpos					// 660, long exposure in millisec, used if exposure=-1
	int32	lnoscan					// 664
	int32	lavgexp					// 668
	char	stripfil[16]			// 672
	char	version[16]			// 688, version & date: 01.000 02/01/90"
	int16	controllerType			// 704		1 = new 120 (Type II)
//												2 = old120 (Type I)
//												3 = ST130
//												4 = ST121
//												5 = ST138
//												6 = DC131 (PentaMAX)
//												7 = ST133 (MicroMAX/SpectroMAX)
//												8 = ST135 (GPIB)
//												9 = VICCD
//												10 = ST116 (GPIB)
//												11 = OMA3 (GPIB)
//												12 = OMA4
	int16	flatFieldApplied			// 706, set to 1 if flat field was applied
	double dummy10[100]			//		uses up 800 bytes
	uint16	dummy11				//		uses up 2 bytes
	int16	NumROI					// 1510
	struct speROIinfoblk ROIinfo[10]// 1512
	char	flatFile1[60]			// 1652, 120 long, break into two parts
	char	flatFile2[60]
	char	bkgFile1[60]			// 1752, 120 long, break into two parts
	char	bkgFile2[60]
	char	blemishFile1[60]
	char	blemishFile2[60]
	float	softwave_version
EndStructure
//
Structure speROIinfoblk				// the ROIinfoblk used with version 1.6 spe header
	uint16	startx					// left x start value
	uint16	endx					// right x value.
	uint16	groupx					// amount x is binned/grouped in hw
	uint16	starty					// top y start value
	uint16	endy					// bottom y value
	uint16	groupy					// amount y is binned/grouped in hw
EndStructure
//
//
Function WinViewList2Header(list,header)				// take list of values (from wnote) put into header
	String list											// list with the values used to fill header
	STRUCT speHeader &header							// structure containing all of header

	String WinViewHeaderStruct=StrVarOrDefault("root:Packages:WinViewHeaderStruct","")
	StructGet/S/B=2 header, WinViewHeaderStruct

	header.version = "2.5.12 Sept2002"
	header.exposure = NumberByKey("exposure",list,"=")
	header.xDimDet = NumberByKey("xDimDet",list,"=")
	header.yDimDet = NumberByKey("yDimDet",list,"=")
	header.ehour = str2num(StringByKey("timeExposed",list,"="))
	header.eminute = str2num(StringFromList(1,StringByKey("timeExposed",list,"="),":"))
	header.ehour = header.ehour==65535 ? 0 : header.ehour
	header.eminute = header.eminute==65535 ? 0 : header.eminute
	header.edate = StringByKey("dateExposed",list,"=")
	header.xdim = NumberByKey("xdim",list,"=")
	header.ydim = NumberByKey("ydim",list,"=")
	header.geometric   = NumberByKey("geo_rotate",list,"=")
	header.geometric += 2*NumberByKey("geo_reverse",list,"=")
	header.geometric += 4*NumberByKey("geo_flip",list,"=")
	header.datatype = NumberByKey("numType",list,"=")
	header.ADCrate = NumberByKey("ADCrate",list,"=")
	header.ADCtype = NumberByKey("ADCtype",list,"=")
	header.ADCresolution = NumberByKey("ADCresolution",list,"=")
	header.controllerType = WhichListItem(StringByKey("controllerType",list,"="),ksControllerTypes)+1

	Variable c
	String bkgFile=StringByKey("bkgFile",list,"="), flatFile=StringByKey("flatFile",list,"=")
	header.bkgApplied = NumberByKey("bkgApplied",list,"=")==1
	if (header.bkgApplied)
		header.bkgFile1 = bkgFile[0,59]
		c = char2num(bkgFile[60,60])
		if (strlen(bkgFile)>59)
			header.bkgFile1[59] = c
			header.bkgFile2 = bkgFile[61,119]
		endif
	endif
	header.flatFieldApplied = NumberByKey("flatFieldApplied",list,"=")==1
	if (header.flatFieldApplied)
		header.flatFile1 = flatFile[0,59]
		c = char2num(flatFile[60,60])
		if (strlen(flatFile)>59)
			header.flatFile1[59] = c
			header.flatFile2 = flatFile[61,119]
		endif
	endif

	header.NumROI = 1								// only implemented for one image
	header.ROIinfo[0].startx = NumberByKey("startx",list,"=")
	header.ROIinfo[0].endx    = NumberByKey("endx",list,"=")
	header.ROIinfo[0].groupx= NumberByKey("groupx",list,"=")
	header.ROIinfo[0].starty = NumberByKey("starty",list,"=")
	header.ROIinfo[0].endy    = NumberByKey("endy",list,"=")
	header.ROIinfo[0].groupy= NumberByKey("groupy",list,"=")

	String str
	str = ReplaceStringByKey("X1",""  , StringByKey("X1",list,"="),"=")
	str = ReplaceStringByKey("Y1",str, StringByKey("Y1",list,"="),"=")
	str = ReplaceStringByKey("Z1",str, StringByKey("Z1",list,"="),"=")
	header.userStr0 = str

	if (numtype(NumberByKey("X2",list,"=")) && !numtype(NumberByKey("depthSi",list,"=")))
		str = ReplaceStringByKey("yc","", StringByKey("yc",list,"="),"=")
		str = ReplaceStringByKey("depthSi",str, StringByKey("depthSi",list,"="),"=")
		str = ReplaceStringByKey("depth0",str  , StringByKey("depth0",list,"="),"=")
	else
		str = ReplaceStringByKey("X2",""  , StringByKey("X2",list,"="),"=")
		str = ReplaceStringByKey("Y2",str, StringByKey("Y2",list,"="),"=")
		str = ReplaceStringByKey("Z2",str, StringByKey("Z2",list,"="),"=")
	endif
	header.userStr1 = str

	str = ReplaceStringByKey("cnt1",""  , StringByKey("cnt1",list,"="),"=")
	str = ReplaceStringByKey("cnt2",str, StringByKey("cnt2",list,"="),"=")
	str = ReplaceStringByKey("cnt3",str, StringByKey("cnt3",list,"="),"=")
	header.userStr2 = str

	str = ReplaceStringByKey("taper",""  , StringByKey("taper",list,"="),"=")
	str = ReplaceStringByKey("keV",str, StringByKey("keV",list,"="),"=")
	str = ReplaceStringByKey("x'",str, StringByKey("x'",list,"="),"=")
	header.userStr3 = str

	str = ReplaceStringByKey("CCDy",""  , StringByKey("CCDy",list,"="),"=")
	str = ReplaceStringByKey("mA",str, StringByKey("mA",list,"="),"=")
	str = ReplaceStringByKey("gap",str, StringByKey("gap",list,"="),"=")
	header.userStr4 = str
	return 0
End


//	Function/T WinViewReadHeader(fileName)
//		String fileName					// fully qualified name of file to open (will not prompt)
//	
//		Variable fid						// file id (file is assumed already opened)
//		Open /Z/M=".spe file"/R/T="????" fid as fileName		// this acutally opens file
//		if (V_flag)
//			return ""					// could not open file
//		endif
//		fileName = S_fileName
//		Make /O/W/U/N=1000 buffer_temp__
//		buffer_temp__ = -1
//		FBinRead /F=2/U/B=3 fid,buffer_temp__
//		Close fid
//		Variable xdim,ydim,itype
//		Variable i,j,offset,exposure
//		String wnote=""					// wave note to add to file read in
//	
//		Variable bkgApplied = buffer_temp__(150/2)
//		String bkgFile=""				// name of bkg file if it exists
//		if (bkgApplied)
//			bkgFile = wave2string(buffer_temp__,1752,120)		// extract backgroud file name it starts 1752 bytes into the file
//			bkgFile = OnlyWinFileName(bkgFile)
//		endif
//		if (bkgApplied %& (strlen(bkgFile)>0))
//			wnote = AddListItem("bkgFile="+bkgFile,wnote)
//		endif
//	
//		itype = buffer_temp__(108/2)
//		wnote = AddListItem("num_Type="+WinViewFileTypeString(itype),wnote)
//		wnote = AddListItem("numType="+num2istr(itype),wnote)
//	
//		xdim = buffer_temp__(21)				// x size of image in file
//		wnote = AddListItem("xdim="+num2istr(xdim),wnote)
//		ydim = buffer_temp__(328)			// y size of image in file
//		wnote = AddListItem("ydim="+num2istr(ydim),wnote)
//	
//		i = buffer_temp__(6/2)				// detector x dimension of chip
//		wnote = AddListItem("xDimDet="+num2istr(i),wnote)
//		i = buffer_temp__(18/2)				// detector y dimension of chip
//		wnote = AddListItem("yDimDet="+num2istr(i),wnote)
//	
//		String fileDate = wave2string(buffer_temp__,20,10)	// date that exposure was taken
//		if (strlen(fileDate)>0)
//			wnote = AddListItem("dateExposed="+fileDate,wnote)
//		endif
//		i = buffer_temp__(30/2)				// hour of the exposure
//		j = buffer_temp__(32/2)				// minute of the exposure
//		if (i!=0 || j!=0)
//			wnote = AddListItem("timeExposed="+num2istr(i)+":"+num2istr(j),wnote)
//		endif
//	
//		if (buffer_temp__(706/2))			// a flat field was applied
//			wnote = AddListItem("flatFieldApplied=1",wnote)
//			String flatFile = wave2string(buffer_temp__,1652,120)		// extract flat file name it starts 1652 bytes into the file
//			flatFile = OnlyWinFileName(flatFile)
//			if (strlen(flatFile)>0)
//				wnote = AddListItem("flatFile="+flatFile,wnote)
//			endif
//	
//		endif
//	
//		i = buffer_temp__(190/2)				// ADC rate
//		wnote = AddListItem("ADCrate="+num2istr(i),wnote)
//		i = buffer_temp__(192/2)				// ADC type
//		wnote = AddListItem("ADCtype="+num2istr(i),wnote)
//		i = buffer_temp__(194/2)				// ADC resolution
//		wnote = AddListItem("ADCresolution="+num2istr(i),wnote)
//	
//		i = buffer_temp__(600/2)				// geometric 1==rotate, 2==reverse, 4==flip
//		if (i&1)
//			wnote = AddListItem("geo_rotate=1",wnote)
//		endif
//		if (i&2)
//			wnote = AddListItem("geo_reverse=1",wnote)
//		endif
//		if (i&4)
//			wnote = AddListItem("geo_flip=1",wnote)
//		endif
//	
//		//	for exposure use  float (4bytes starting at byte10)
//		Make /O/R/N=2 bufferFloat_temp__
//		Open /Z/M=".spe file"/R/T="????" fid as fileName		// this acutally opens file
//		if (V_flag)
//			return ""					// could not open file
//		endif
//		FSetPos fid, 10
//		FBinRead /F=4/B=3 fid,bufferFloat_temp__
//		Close fid
//		exposure = bufferFloat_temp__[0]
//		KillWaves/Z bufferFloat_temp__
//		wnote = AddListItem("exposure="+num2str(exposure),wnote)
//	
//		// String controllerTypes="new120 (TYPE II);old120 (TYPE I);ST130;ST121;ST138;DC131 (Pentamax);ST133 (MicroMAX/SpectroMAX);ST135 (GPIB);VICCD;ST116 (GPIB);OMA3 (GPIB);OMA4"
//		i = buffer_temp__(704/2)	-1			// controller type
//		wnote = AddListItem("controllerType="+StringFromList(i,ksControllerTypes),wnote)
//	
//		i = buffer_temp__(1510/2)			// number of ROI's
//		Variable nROI = max(i,1)				// zero means one
//		if (nROI>1)
//			Abort "this file contains multiple ROI's"
//			wnote = AddListItem("NumROI="+num2istr(nROI),wnote)
//		endif
//		String str, ROInames = "startx;endx;groupx;starty;endy;groupy"
//		for (j=1; j<=nROI; j+=1)
//			offset = (1512/2) + (j-1)*12/2
//			for (i=0;i<6;i+=1)
//				if (nROI>1)
//					sprintf str,"%s_%d=%d",StringFromList(i,ROInames),j,buffer_temp__(offset+i)
//				else
//					sprintf str,"%s=%d",StringFromList(i,ROInames),buffer_temp__(offset+i)
//				endif
//				wnote = AddListItem(str,wnote)
//			endfor
//		endfor
//	
//		// read in each of the 5 lists they start at 200, every 80 (200,280, 360, 440, 520)
//		String item,list
//		Variable n
//		for (j=200;j<521;j+=80)							// byte index into the file
//			list = wave2string(buffer_temp__,j,80)		// extract a null terminated string from buf starting at index, and going no more than maxLen characters
//			n = strsearch(list," ",0)
//			list = list[n+1,inf]
//			n=Itemsinlist(list)
//			for (i=0;i<n;i+=1)
//				item = StringFromList(i,list)
//				if (strlen(item)>1)
//					wnote = AddListItem(item,wnote)
//				endif
//			endfor
//		endfor
//	
//		KillWaves/Z buffer_temp__
//	
//	//	printf "for file '"+fileName+"'"
//	//	if (bkgApplied %& (strlen(bkgFile)>0))
//	//		printf ",             background file = '%s'\r",  bkgFile
//	//	else
//	//		printf "\r"
//	//	endif
//	//	printf "xdim = %d      ydim = %d      ", xdim,ydim
//	//	printf "total length = %d x %d  = %d points\r", xdim,ydim,xdim*ydim
//	//	print "number type is  '"+WinViewFileTypeString(itype)+"'"
//	
//		return wnote
//	End
//
//Function/S test_WinViewReadHeader(fName)
//	String fName											// fully qualified name of file to open
//	Variable refNum
//	if (strlen((OnlyWinFileName(fName)))<1)				// call dialog if no file name passed
//		Open /D/M=".spe file"/R/T=".SPE" refNum			// use /D to get full path name
//		fName = S_filename
//	endif
//	if (strlen(fName)<1)									// no file name, quit
//		return ""
//	endif
//	String wnote =  WinViewReadHeader(fName)
//	print wnote
//	return wnote
//End



Function/T WinViewLoadROI(fileName,itype,xdim,i0,i1,j0,j1,[extras])
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
	String extras					// not used here for anything yet
	extras = SelectString(ParamIsDefault(extras),extras,"")
	if (!(itype>=0) || !(xdim>0))			// invalid values, so read them from file
		Variable fid							// file id
		Open /Z/M=".spe file"/R/T="????" fid as fileName		// this acutally opens file
		if (V_flag)
			return ""							// could not open file
		endif
		FStatus fid
		if (V_logEOF<4100)					// file is too short, do nothing
			Close fid
			return ""
		endif
		FSetPos fid, 42
		FBinRead/B=3/f=2/U fid,xdim		// xdim is at byte 42
		FSetPos fid, 108
		FBinRead/B=3/f=2/U fid,itype		// itype is at byte 108
		Close fid
	endif

	Variable bytes=1							// length of image nuber type in bytes
	switch(itype)
		case 5:
			bytes *=2
		case 0:
		case 1:
			bytes *=2
		case 2:
		case 3:
			bytes *=2
	endswitch

	i0 = max(round(i0),0)
	i1 = max(round(i1),0)
	j0 = max(round(j0),0)
	j1 = max(round(j1),0)
	Variable nx = i1-i0+1
	Variable ny = j1-j0+1
	if (nx<1 || ny<1)							// nothing to read
		return ""
	endif
	Variable skip = 4100+bytes*j0*xdim	// offset (bytes) to start of roi, 4100 is to start of image
	skip += ny==1 ? bytes*i0 : 0
//	if (ny==1)
//		skip += bytes*i0
//	endif
	Variable fType = iType2fType(itype)	// Igor number type
	String command
	if (ny>1)
		sprintf command, "GBLoadWave/Q/B/T={%d,%d}/S=%d/W=1/U=%ld \"%s\"",fType,fType,skip,xdim*ny,fileName
	else
		sprintf command, "GBLoadWave/Q/B/T={%d,%d}/S=%d/W=1/U=%ld \"%s\"",fType,fType,skip,nx*ny,fileName
	endif
	command += " ; Variable/G WinView_Nwaves=V_flag"
	Execute command
	NVAR WinView_Nwaves=WinView_Nwaves
	if (WinView_Nwaves<1)
		return ""							// nothing loaded
	endif
	SVAR S_waveNames=S_waveNames
	String wName = S_waveNames[0,strsearch(S_waveNames, ";", 0)-1]
	Wave wav = $wName

	if (nx<xdim && ny>1)					// compress up in the x direction
		Variable i
		Variable next=0						// points to starting point in reduced array to store
		Variable ystart=0					// points to start of this y-value
		for (i=0;i<ny;i+=1)				// loop over each y-value of the array
			wav[next,next+nx-1] = wav[p-next+ystart+i0]	// move x range to next next point in array
			ystart += xdim;					// points to start next y-value
			next += nx;						// points to the next place in reduced array
		endfor
	endif

	Redimension /N=(nx,ny) wav
	return GetWavesDataFolder(wav,2)
End



Function WinViewInfo(image,key)
	Wave image
	String key

	String imageName
	Variable i
	if (!WaveExists(image))
		String item=""
		Prompt imageName,"image to use", popup, WaveListClass("speImage*","*","DIMS:2")
		DoPrompt "pick an image", imageName
		if(V_flag)
			return NaN
		endif
		Wave image=$imageName
	endif
	imageName = NameOfWave(image)
	if (!WaveExists(image))
		return NaN
	endif
	String list = note(image)
//	 possibly extend list to include H and F (from Y and Z)
	Variable Y,Z,H,F
	Y = NumberByKey("Y1",list,"=")
	Z = NumberByKey("Z1",list,"=")
	H = NumberByKey("H1",list,"=")			// get H & F to test for their existance
	F = NumberByKey("F1",list,"=")
	if (!numtype(Y+Z) && numtype(H+F))		// if Y1 & Z1 exist, and H1 & F1 do not
		H = (Z+Y)/sqrt(2)
		F = (Z-Y)/sqrt(2)
		list = ReplaceNumberByKey("H1",list,H,"=")
		list = ReplaceNumberByKey("F1",list,F,"=")
	endif

	Y = NumberByKey("Y2",list,"=")
	Z = NumberByKey("Z2",list,"=")
	H = NumberByKey("H2",list,"=")
	F = NumberByKey("F2",list,"=")
	if (!numtype(Y+Z) && numtype(H+F))		// if Y2 & Z2 exist, and H2 & F2 do not
		H = (Z+Y)/sqrt(2)
		F = (Z-Y)/sqrt(2)
		list = ReplaceNumberByKey("H2",list,H,"=")
		list = ReplaceNumberByKey("F2",list,F,"=")
	endif

	if (strlen(key)>1)
		item = StringByKey(key,list,"=")
		if (strlen(item)>1 && numtype(str2num(item)))
			printf "for '%s',  %s = '%s'\r",imageName,key,item
		endif
	else
		Prompt item,"item",popup,list
		DoPrompt "info about "+imageName, item
		if(V_flag)
			return NaN
		endif
		key = StringFromList(0,item,"=")
		item = StringFromList(1,item,"=")
		printf "for '%s',  %s = %s\r",imageName,key,item
	endif
	return str2num(item)
End



//Function/T OnlyFileName(full)
//	String full
//
//	String name=full
//	Variable ii=-1
//
//	do
//		ii = strsearch(name, ":", 0)			// remove the path part
//		if (ii>0)
//			name= name[ii+1,inf]
//		endif
//	while(ii>0)
//
//	ii = strsearch(name, ".SPE", 0)			// trim off trailing '.SPE'
//	if (ii<0)
//		ii = strsearch(name, ".spe", 0)		// trim off trailing '.spe'
//	endif
//	if (ii>1)
//		name = name[0,ii-1]
//	endif
//
////	do
////		ii = strsearch(name, ".", 0)			// change all '.' to '_'
////		if (ii>=0)
////			name[ii,ii]="_"
////		endif
////	while(ii>=0)
////
////	if (char2num(name[0,0])<char2num("A"))	// must start with a letter
////		name = "X"+name
////	endif
//	name = CleanupName(name,0)
//	return name
//End



Function/T OnlyWinFileName(full)
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



// no longer needed, used by the old version of WinViewReadHeader(fileName), with structs, this is not needed.
//
//	Function/T wave2string(buf,offset,maxLen)		// extract a null terminated string from buf starting at index, and going no more than maxLen characters
//		Wave buf
//		Variable offset					// offset in bytes to the start of string (bytes)
//		Variable maxLen				// maximum length of string (bytes)
//	
//		Variable itype = NumberByKey("NUMTYPE", WaveInfo(buf,0))
//		Variable bytes=1				// number of bytes per index into buf
//		itype = itype & 62
//		switch(itype)
//			case 4:
//				bytes *=2
//			case 2:
//			case 32:
//				bytes *=2
//			case 16:
//				bytes *=2
//		endswitch
//	
//		offset /= bytes
//		String str=""									// result
//		Variable shift,int, i=0
//		do
//			int = buf(offset + floor(i/bytes))
//			if (bytes==1)
//				int = int
//			endif
//			if (bytes==2)
//				if (mod(i,2)==0)
//					int = int & 255						// use low byte
//				else
//					int = (int & (255*256))/256		// use next byte
//				endif
//			endif
//			if (bytes==4)
//				if (mod(i,4)==0)
//					int = int & 255						// use low byte
//				elseif (mod(i,4)==1)
//					shift = 256
//					int = (int & (255*shift))/shift		// use next byte
//				elseif (mod(i,4)==2)
//					shift = 256^2
//					int = (int & (255*shift))/shift		// use next byte
//				else
//					shift = 256^3
//					int = (int & (255*shift))/shift		// use next byte
//				endif
//			endif
//			if (int<1)
//				break
//			endif
//			str[i,i] = num2char(int)
//			i += 1
//		while(i<maxLen)
//		return str
//	End



Function/T WinViewFileTypeString(itype)
	Variable itype

	String stype=""
	if (itype==0)
		stype = "float (4 byte)"
	endif
	if (itype==1)
		stype = "long integer (4 byte)"
	endif
	if (itype==2)
		stype = "integer (2 byte)"
	endif
	if (itype==3)
		stype = "unsigned integer (2 byte)"
	endif
	if (itype==4)
		stype = "string/char (1 byte)"
	endif
	if (itype==5)
		stype = "double (8 byte)"
	endif
	if (itype==6)
		stype = "signed int8 (1 byte)"
	endif
	if (itype==7)
		stype = "unsigned int8 (1 byte)"
	endif
	if (strlen(stype)<1)
		stype="unknown number type"
	endif
	return stype
End


Function iType2fType(itype)		// converts winview number type to Igor number type
	Variable itype

	//fType   in Igor
	//	2		single-precision floating point
	//	4		double-precision floating point
	//	32		32 bit signed integer
	//	16		16 bit signed integer
	//	8		8 bit signed integer
	//	32+64	32 bit signed integer
	//	16+64	16 bit signed integer
	//	8+64		8 bit signed integer
	//
	//
	//		for .spe itype==
	//  0	"float (4 byte)"
	//  1	"long integer (4 byte)"
	//  2	"integer (2 byte)"
	//  3	"unsigned integer (2 byte)"
	//  4	"string/char (1 byte)"
	//  5	"double (8 byte)"
	//  6	"signed int8 (1 byte)"
	//  7	"unsigned int8 (1 byte)"

	Variable ftype=-1
	if (itype==0)
		ftype = 2
	endif
	if (itype==1)
		ftype = 32
	endif
	if (itype==2)
		ftype = 16
	endif
	if (itype==3)
		ftype = 16+64
	endif
	if (itype==4)
		ftype = 8+64
	endif
	if (itype==5)
		ftype = 4
	endif
	if (itype==6)
		ftype = 8
	endif
	if (itype==7)
		ftype = 8+64
	endif
	if (ftype<0)
		DoAlert 0, "unknown number type"
	endif
	return ftype
End

Function igorType2WinView(itype)	// convert Igor WaveType() to WinView number types
	Variable itype					// an igor number type from WaveType()
	switch(itype)
		case 0x02:					// 4 byte float
			return 0
		case 0x04:					// double (8 byte)
			return 5
		case 0x08:					// signed int8 (1 byte)
			return 6
		case 0x10:					// integer (2 byte)
			return 2
		case 0x20:					// long integer (4 byte)
			return 1
		case 0x48:					// unsigned int8 (1 byte)
			return 7
		case 0x50:					// unsigned integer (2 byte)
			return 3
	endswitch
	return NaN						// invalid wave type (complex, unsigned floats, unsigned long)
End



// =======================================================================


Function/S Graph_imageMake(image,withButtons)
	Wave image
	Variable withButtons
	if (!WaveExists(image) || numtype(withButtons))
		String wList = reverseList(WaveListClass("speImage*","*","DIMS:2"))
//		String wName=SelectString(WaveExists(image),StringFromList(0,WaveList("*",";","DIMS:2")),GetWavesDataFolder(image,2))
		String wName=SelectString(WaveExists(image),StringFromList(0,wList),NameOfWave(image))
		withButtons = numtype(withButtons) ? 0 : withButtons
		withButtons = withButtons ? 1 : 2
		Prompt withButtons,"put the buttons on the graph",popup,"Yes;No"
//		Prompt wName, "image to show", popup, WaveList("*",";","DIMS:2")
//		Prompt wName, "image to show", popup, WaveListClass("speImage;speImageNoBkg","*","DIMS:2")
		Prompt wName, "image to show", popup, wList
		DoPrompt "choose image", wName, withButtons
		if(V_flag)
			return ""
		endif
		withButtons = (withButtons==1) ? 1 : 0
		Wave image = $wName
	endif
	if (!WaveExists(image))
		DoAlert 0, "image wave does not exist"
		return ""
	endif

	Display /W=(345,44,822,440)
	AppendImage image
	if (withButtons)
		String/G pkList
		String command
		Button button0,pos={0,26},size={65,20},proc=ButtonProc_Include,title="include"
		PopupMenu popup0,pos={0,4},size={72,20},proc=PopMenuProc,mode=12
		sprintf command, "	PopupMenu popup0,value=#\"%s\"", GetDataFolder(1)+"pkList"
		Execute command
		Button boxes,pos={0,48},size={65,21},proc=ButtonProc_Boxes,title="�boxes"
	endif
	Execute "GraphImageStyle()"
	if (exists("getIndexedPeakInfoHook")==6)
		SetWindow kwTopWin hook(peakInfo)=getIndexedPeakInfoHook
	endif
	return GetWavesDataFolder(image,2)
//	return StringFromList(0,WinList("*",";","WIN:1"))
End
//Function/T Graph_imageMake(wName)
//	String wName
//	if (strlen(wName)<1)
//		wName = StringFromList(0,WaveList("*",";","DIMS:2"))
//		Prompt wName, "image to show", popup, WaveList("*",";","DIMS:2")
//		DoPrompt "choose image", wName
//		if(V_flag)
//			return ""
//		endif
//	endif
//	PauseUpdate
//
//	Display /W=(345,44,822,440)
//	AppendImage $wName
//	Execute "GraphImageStyle()"
//
//	String/G pkList
//	String command
//	Button button0,pos={0,26},size={65,20},proc=ButtonProc_Include,title="include"
//	PopupMenu popup0,pos={0,4},size={72,20},proc=PopMenuProc,mode=12
//	sprintf command, "	PopupMenu popup0,value=#\"%s\"", GetDataFolder(1)+"pkList"
//	Execute command
//	Button boxes,pos={0,48},size={65,21},proc=ButtonProc_Boxes,title="�boxes"
//	return StringFromList(0,WinList("*",";","WIN:1"))
//End


Function ButtonProc_Include(ctrlName) : ButtonControl
	String ctrlName
	Variable i=NaN,j=NaN
	GetMarquee/Z
	if (V_flag)								// use marquee
		AddGaussianToList()
		i = K2
		j= K4
	elseif (strlen(CsrWave(A))>1)
		i = pcsr(A)							// use cursor
		j= qcsr(A)
		AddPoint2pkList(i,j)
	else
		DoAlert 0, "no Marquee selected and cursor A not on graph, do nothing"
	endif
End
Function PopMenuProc(ctrlName,popNum,popStr) : PopupMenuControl
	String ctrlName
	Variable popNum
	String popStr

	SVAR pkList=pkList
	NVAR pkIndex=pkIndex
	Variable i,n=ItemsInList(pkList)
	for (i=0;i<n;i+=1)
		if (stringmatch(StringFromList(i,pkList),popStr))
			break
		endif
	endfor
	if (i<n)
		pkIndex = i
	endif

	String wName = StringFromList(0,ImageNameList("",";"))
	Variable xx=str2num(StringFromList(0,popStr,","))
	Variable yy=str2num(StringFromList(1,popStr,","))
	Cursor/P/S=2/A=1/L=1/I A $wName xx,yy
End
//Function OLDButtonProc_Boxes(ctrlName) : ButtonControl
//	String ctrlName
//
//	String fldr = GetWavesDataFolder(ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";"))),1 )
//	NVAR boxesOn = $(fldr+"boxesOn")
//	SVAR pkList = $(fldr+"pkList")
//	if (!NVAR_Exists(boxesOn))
//		Variable/G $(fldr+"boxesOn")
//	endif
//	if (!SVAR_Exists(pkList))
//		String/G $(fldr+"pkList")
//	endif
//
//	boxesOn = !boxesOn
//	SetDrawLayer /K UserFront				// always clear first
//	if (!boxesOn)
//		return 0
//	endif
//
//	Variable x0,y0,xlo,xhi,ylo,yhi
//	Variable i
//	String item
//	Variable wid=5*2
//	item=StringFromList(0,pkList)
//	SetDrawLayer UserFront
//	for (i=0; strlen(item)>0;  i+=1,item=StringFromList(i,pkList))
//		x0=str2num(StringFromList(0,item,","))
//		y0=str2num(StringFromList(1,item,","))
//		if (numtype(x0) || numtype(y0))
//			continue
//		endif
//		SetDrawEnv xcoord= bottom,ycoord= left,fillpat= 0
//		DrawRect x0-wid,y0-wid,x0+wid,y0+wid
//		SetDrawEnv xcoord= bottom,ycoord= left,linefgc=(30583,30583,30583)
//		DrawLine x0,y0-wid,x0,y0-wid/2
//		SetDrawEnv xcoord= bottom,ycoord= left,linefgc=(30583,30583,30583)
//		DrawLine x0,y0+wid/2,x0,y0+wid
//		SetDrawEnv xcoord= bottom,ycoord= left,linefgc=(30583,30583,30583)
//		DrawLine x0-wid,y0,x0-wid/2,y0
//		SetDrawEnv xcoord= bottom,ycoord= left,linefgc=(30583,30583,30583)
//		DrawLine x0+wid/2,y0,x0+wid,y0
//	endfor
//	return 1
//End
Function ButtonProc_Boxes(ctrlName) : ButtonControl
	String ctrlName

	String fldr = GetWavesDataFolder(ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";"))),1 )
	NVAR boxesOn = $(fldr+"boxesOn")
	SVAR pkList = $(fldr+"pkList")
	if (!NVAR_Exists(boxesOn))
		Variable/G $(fldr+"boxesOn")
	endif
	if (!SVAR_Exists(pkList))
		String/G $(fldr+"pkList")
	endif

	boxesOn = !boxesOn
	SetDrawLayer /K UserFront				// always clear first
	if (!boxesOn)
		String list = GetUserData("","","Indexing")
		SetWindow kwTopWin,userdata(FitPeaks)=RemoveByKey("FullPeakList",list,"=")
		return 0
	endif
	ReDrawBoxes(pkList)
	return 1
End
//Function ButtonProc_Boxes(ctrlName) : ButtonControl
//	String ctrlName
//
//	String fldr = GetWavesDataFolder(ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";"))),1 )
//	NVAR boxesOn = $(fldr+"boxesOn")
//	SVAR pkList = $(fldr+"pkList")
//	if (!NVAR_Exists(boxesOn))
//		Variable/G $(fldr+"boxesOn")
//	endif
//	if (!SVAR_Exists(pkList))
//		String/G $(fldr+"pkList")
//	endif
//
//	boxesOn = !boxesOn
//	SetDrawLayer /K UserFront				// always clear first
//	if (!boxesOn)
//		String list = GetUserData("","","Indexing")
//		SetWindow kwTopWin,userdata(FitPeaks)=RemoveByKey("FullPeakList",list,"=")
//		return 0
//	endif
//
//	Variable x0,y0,xlo,xhi,ylo,yhi
//	Variable i
//	String item
//	Variable wid=5*4
//	item=StringFromList(0,pkList)
//	SetDrawLayer UserFront
//	for (i=0; strlen(item)>0;  i+=1,item=StringFromList(i,pkList))
//		x0=str2num(StringFromList(0,item,","))
//		y0=str2num(StringFromList(1,item,","))
//		DrawMarker(x0,y0,wid,wid,"BoxWithTicks")
//	endfor
//	return 1
//End

Function ReDrawBoxes(pks)
	String pks

	if (strlen(pks)<1)
		String fldr = GetWavesDataFolder(ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";"))),1 )
		SVAR pkList = $(fldr+"pkList")
		if (!SVAR_Exists(pkList))
			return 1
		endif
		pks = pkList
	endif

	Variable x0,y0,xlo,xhi,ylo,yhi
	Variable i
	String item
	Variable wid=5*2				// 5*4
	item=StringFromList(0,pks)
	SetDrawLayer /K UserFront				// always clear first
	SetDrawLayer UserFront
	for (i=0; strlen(item)>0;  i+=1,item=StringFromList(i,pks))
		x0=str2num(StringFromList(0,item,","))
		y0=str2num(StringFromList(1,item,","))
		DrawMarker(x0,y0,wid,wid,"BoxWithTicks")
	endfor
	return 0
End



Static Function DrawMarker(x0,y0,dx,dy,style,[color,thick,win])
	Variable x0,y0,dx,dy
	String style
	String color							// one of "red", "blue", "green", "yellow", "magenta", "cyan", or a triplet like "1,1,0" or "65535,0,0"
	Variable thick						// line thickness, defaults to 0.50 (thin)
	String win							// optional name of window
	if (numtype(x0++y0+dx+dy) || dx<=0 || dy<=0)
		return 1
	endif

	if (ParamIsDefault(color))
		color = SelectString(stringmatch(style,"BoxWithTicks"),"black","30583,30583,30583")	// special default color for BoxWithTicks
	endif
//	color = SelectString(ParamIsDefault(color),color,"black")	// default color to black
//	color = SelectString(ParamIsDefault(color) && stringmatch(style,"BoxWithTicks"),color,"30583,30583,30583")	// special default color for BoxWithTicks
	thick = ParamIsDefault(thick) ? 0.5 : thick
	if (ParamIsDefault(win))
		win = ""
	endif
	String rgb = StringByKey(color,"red:1,0,0;blue:0,0,1;green:0,1,0;yellow:1,1,0;magenta:1,0,1;cyan:0,1,1;black:0,0,0;white:1,1,1")	// name -> numbers
	rgb = SelectString(strlen(rgb),color,rgb)
	Variable r=str2num(StringFromList(0,rgb,",")), g=str2num(StringFromList(1,rgb,",")), b=str2num(StringFromList(2,rgb,","))
	r = !(r>=0) ? 0 : r											// default invalid to black
	g = !(g>=0) ? 0 : g
	b = !(b>=0) ? 0 : b
	if (r+g+b<=3)												// if numbers are all small change from [0,1] -> [0,65535]
		r *= 65535
		g *= 65535
		b *= 65535
	endif

	style = SelectString(strlen(style),"cross",style)			// defaults to Cross
	String list = ImageInfo(win,"",0)
	String xaxis=StringByKey("XAXIS",list), yaxis=StringByKey("YAXIS",list)
	dy = abs(dy/2)		// change FW to HW
	dx = abs(dx/2)
	SetDrawLayer/W=$win UserFront

	if (stringmatch(style,"cross"))
		SetDrawEnv/W=$win  xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick
		DrawLine/W=$win  (x0),(y0-dy),x0,(y0+dy)
		SetDrawEnv/W=$win xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick
		DrawLine/W=$win (x0-dx),y0,(x0+dx),y0
	elseif (stringmatch(style,"X"))
		SetDrawEnv/W=$win xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick
		DrawLine/W=$win (x0-dx),(y0-dy),(x0+dx),(y0+dy)
		SetDrawEnv/W=$win xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick
		DrawLine/W=$win (x0+dx),(y0-dy),(x0-dx),(y0+dy)
	elseif (stringmatch(style,"BoxWithTicks"))
		Variable hw							// half width of box in pixels, usually ~10
		SetDrawEnv/W=$win xcoord= $xaxis,ycoord= $yaxis,fillpat= 0, linethick= thick
		DrawRect x0-dx,y0-dy,x0+dx,y0+dy
		SetDrawEnv/W=$win xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick
		DrawLine/W=$win x0,y0-dy,x0,y0-dy/2
		SetDrawEnv/W=$win xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick
		DrawLine/W=$win x0,y0+dy/2,x0,y0+dy
		SetDrawEnv/W=$win xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick
		DrawLine/W=$win x0-dx,y0,x0-dx/2,y0
		SetDrawEnv/W=$win xcoord= $xaxis,ycoord= $yaxis,linefgc=(r,g,b), linethick= thick
		DrawLine/W=$win x0+dx/2,y0,x0+dx,y0
	else
		return 1
	endif
	return 0
End
//
Function DrawBoxWithTicksAtXY(x0,y0,hw)		// Deprecated, do not use, call DrawMarker() directly
	Variable x0,y0
	Variable hw							// half width of box in pixels, usually ~10
	return DrawMarker(x0,y0,2*hw,2*hw,"BoxWithTicks")
End
//
//Function DrawBoxWithTicksAtXY(x0,y0,hw)
//	Variable x0,y0
//	Variable hw							// half width of box in pixels, usually ~10
//	if (numtype(x0) || numtype(y0) || numtype(hw) || hw<=0)
//		return 1
//	endif
//	String list = ImageInfo("","",0)
//	String xaxis=StringByKey("XAXIS",list), yaxis=StringByKey("YAXIS",list)
//	SetDrawLayer UserFront
//	SetDrawEnv xcoord= $xaxis,ycoord= $yaxis,fillpat= 0
//	DrawRect x0-hw,y0-hw,x0+hw,y0+hw
//	SetDrawEnv xcoord= $xaxis,ycoord= $yaxis,linefgc=(30583,30583,30583)
//	DrawLine x0,y0-hw,x0,y0-hw/2
//	SetDrawEnv xcoord= $xaxis,ycoord= $yaxis,linefgc=(30583,30583,30583)
//	DrawLine x0,y0+hw/2,x0,y0+hw
//	SetDrawEnv xcoord= $xaxis,ycoord= $yaxis,linefgc=(30583,30583,30583)
//	DrawLine x0-hw,y0,x0-hw/2,y0
//	SetDrawEnv xcoord= $xaxis,ycoord= $yaxis,linefgc=(30583,30583,30583)
//	DrawLine x0+hw/2,y0,x0+hw,y0
//	return 0
//End
//
Function DrawCross(x0,y0,dx,dy)
	Variable x0,y0,dx,dy
	return DrawMarker(x0,y0,dx,dy,"Cross")
End
//Function DrawCross(xcm,ycm,dx,dy)
//	Variable xcm,ycm,dx,dy
//	dy = abs(dy/2)		// change FW to HW
//	dx = abs(dx/2)
//	String list = ImageInfo("","",0)
//	String xaxis=StringByKey("XAXIS",list), yaxis=StringByKey("YAXIS",list)
//	SetDrawLayer UserFront
//	SetDrawEnv xcoord= $xaxis,ycoord= $yaxis
//	DrawLine (xcm),(ycm-dy),xcm,(ycm+dy)
//	SetDrawEnv xcoord= $xaxis,ycoord= $yaxis
//	DrawLine (xcm-dx),ycm,(xcm+dx),ycm
End



Function/T RemoveDuplicatePairsFromList(list)
	String list

	Variable tolerance=5
	Variable i,j
	Variable n=ItemsInList(list)
	Variable x0,y0,xx,yy
	String item

	for(j=0;j<n-1;j+=1)						// check all pairs
		item = StringFromList(j,list)
		x0 = str2num(StringFromList(0,item,","))
		y0 = str2num(StringFromList(1,item,","))
		for(i=n-1;i>j;i-=1)						// check all following pairs, starting from the end
			item = StringFromList(i,list)
			xx = str2num(StringFromList(0,item,","))
			yy = str2num(StringFromList(1,item,","))
			if ((abs(xx-x0)<tolerance) && (abs(yy-y0)<tolerance))
				list = RemoveListItem(i,list)
				n -= 1
			endif
		endfor
	endfor
	return list
End

Function/T RemovePairFromList(list,x0,y0,tol)
	String list
	Variable x0,y0
	Variable tol										// � tolerance in pixels, usually~5

	Variable ix,iy,flag=1
	String item
	Variable i=ItemsInList(list)
	do
		i -= 1
		item = StringFromList(i,list)
		ix = str2num(StringFromList(0,item,","))
		iy = str2num(StringFromList(1,item,","))
		if (abs(ix-x0)<=tol && abs(iy-y0)<=tol)
			flag=0
			list = RemoveFromList(item,list)
		endif
	while (i>0)
	return list
End