#pragma rtGlobals=1		// Use modern global access method.
#pragma ModuleName=EnergyScans
#pragma version = 2.51

// version 2.00 brings all of the Q-distributions in to one single routine whether depth or positioner
// version 2.10 cleans out a lot of the old stuff left over from pre 2.00
// version 2.14 speed up Fill_Q_Positions(), by only reading part of image headers
// version 2.30 change Fill1_3DQspace() to allow processing only an roi
// version 2.43 Fill_Q_Positions() and fix Fill1_3DQspace() now automatically use the badPixels wave
// version 2.44 Cleaned out all the old code
// version 2.45 Fill_Q_Positions(), allow for variable ROI in the histograming (only for |Q|)
// version 2.46 Fill_Q_Positions(), allow for variable ROI in the histograming (added 3D Q too)
// version 2.48 added FindStepSizeInVec(vec,threshold), use in FindScalingFromVec(), Fill_Q_Positions(), and Fill1_3DQspace()
// version 2.49 moved FindScalingFromVec() and FindStepSizeInVec() to Utillity_JZT
// version 2.50 Fill_Q_Positions(), added optional dQ, and a DEFAULT_dQfactor to allow user tweaking of dQ
// version 2.51 changed all special Mac character, should print better with Windows & Igor7

#include "ImageDisplayScaling", version>=2.11
#if (Exists("HDF5OpenFile")==4)
#include "HDF5images", version>=0.35
#endif
#if (NumVarOrDefault("root:Packages:MICRO_GEOMETRY_VERSION",0)&4)
#include "WinView", version>=2.03
#endif
#include "Masking", version>1.03
#include "GizmoZoomTranslate", version>=2.03
#include "GizmoClip", version>=2.03
#include "GizmoMarkers", version>=2.21
#include "QspaceVolumesView",  version>=1.21
#include "microGeometryN", version>=1.95

Static Constant hc = 1.239841857					// hc (keV-nm)
Static Constant secPerPixelFixed = 18.8e-6		// the fixed time is takes to process one pixel (sec) after any distortion
Static Constant secPerPixelDistort = 15.24e-3	// time is takes to process distortion for one pixel (sec)  (measured with a 2GHz clock)
Static Constant DEFAULT_I0_GAIN = 1e9
Static Constant DEFAULT_I0_SCALING = 1e5
Static Constant DEFAULT_dQfactor = 1.1
#if (IgorVersion()<7)
	Static strConstant Gmu = "\265"		// Mac option-m, Greek mu
#if StringMatch(IgorInfo(2),"Windows")
	Static strConstant PLUSMINUS = "\241"	// MS Alt 241, plus-minus sign
#else
	Static strConstant PLUSMINUS = "\261"	// Mac option-+, plus-minus sign
#endif
#else
	Static strConstant Gmu = "\xCE\xBC"	// UTF8, Greek mu
	Static strConstant PLUSMINUS = "\xC2\xB1"	// UTF8, plus-minus sign
#endif


//	**	To get the ful chip unbinned pixel position from ROI data use the following:
//	**	fullchip x  =  startx + x*groupx + (groupx-1)/2
//	**	fullchip y  =  starty + y*groupy + (groupy-1)/2
//		py = j*groupy + (groupy-1)/2 + starty
//		px = i*groupx + (groupx-1)/2 + startx
//
//
//				at angle=0,  F=-Y, H=+Z		taken right out of focus5.5.st
//			F = -Y*cos(angle) + Z*sin(angle)
//			H =  Y*sin(angle) + Z*cos(angle)
//			Y =  H*sin(angle) - F*cos(angle)
//			Z =  H*cos(angle) + F*sin(angle)
//
//
//	new for May 3, 2006
//	sz = sample z position
//	s = Æz/Æy of ray from sample to pixel on CCD
//
//	z = (  zs - s*Fo/cos(angle) - R*sqrt(1-s*s)  ) / (1-s*tan(angle) )
//	H = [z/cos(angle) - Fo*tan(angle)]] = [z - Fo*sin(angle)] / cos(angle)
//
//	for a Gaussian,   fwhm = 2*K3*sqrt(ln(2))
//	for a Lorentzian,   fwhm = 2*sqrt(K3)


Menu LaueGoMainMenuName, dynamic
	Submenu "Energy Scans"
		"New Energy scan...",NewEnergyScan("","",NaN)
		help={"Create a new data folder for an energy(+ optional wire) scan, and prompt user for title and d0"}
		"run:  Fill Q-Distribution",Fill_Q_Positions(NaN,"imagePath","","","",$"",printIt=1)
		help={"read the E-scanned images which were optionaly also scanned over positioner(s) or depth sorted, and make a N-d plot of Q vs position"}
			MenuItemIfWaveClassExists("  Re-Plot Q-distribution [re-fit params too]","Qhistogram;QdistAtPositions","MINCOLS:1,MAXCHUNKS:1"),MakeOtherQhistWaves($"")
			MenuItemIfWaveClassExists("  ReCalc an RGB for a Q-distribution","Q_Positions",""),MakeRGBforQdistribution($"",$"")
			MenuItemIfWaveExists("  Graph 1 Q histogram","Qhist"),MakeGraph_Qhist($"")
			help={"re-plot the single slice from the Q-depth surface"}
		MenuItemIfWaveClassExists("Q-Distribution of One Loaded image","rawImage*","DIMS:2"), Fill_Q_1image(NaN,$"")
		help={"just show the Q-distribution from a single image that is alerady loaded"}
		MenuItemIfGraphNamesExist("Layout [Q x Position] + Q-histogram","*_Q"),MakeLayoutQ_Positions_hist("")
		"-"
			"Graph Random Points Q-dist", EnergyScans#MakeGraph_Q_RandomPositions($"",NaN,NaN)
		"-"
		"Make Pseudo White Images",MakePseudoWhiteImages("","")
		      help={"read all the depth sorted images and make the pseudo white beam plot"}
		     MenuItemIfWaveExists("    Make movie of Energy Sums vs depth","imageEsum"),MakeEsumMovie()
		      help={"take the pseudo white data already read in and make the depth movie"}
		     MenuItemIfWaveExists("    Re-Plot Energy Sums vs depth","imageEsum"),MakeEsumPlot()
		      help={"take the pseudo white data already read in and re-plot it"}

		"-"
		"run:  Fill 3D Q-space at One Position",Fill1_3DQspace(0,"imagePath","","",printIt=1)
		"  Make 3D Q-space Gizmo", MakeGizmoQspace3d($"")

		"-"
		MenuItemIfWaveClassExists("Make Mask","rawImage*,speImage*;ImageSummed*","DIMS:2"),MakeMask($"")
		      help={"make a mask from an image (often the pseudo white image)"}
		MarqueeMaskMenuItem("    Clear the Mask"),ClearMask()
		      help={"clear the current mask (usually accessed by clicking in marquee)"}

		"-"
		"Re-scale Q wave to strain",reScaleQwaveToStrain($"",NaN)
		help={"make a nice layout for printing showing the Q-position surface plot and the plot of a cut through it"}
		"-"
		"Select Lineshape Function...", EnergyScans#Escan_SelectLineShapeFunction("")
		"d[hkl] --> d0,     use xtal values",EnergyScans#update_d0()
		help={"reset the d0 used to compute strains using d(nm) from 'Xtal' panel"}
	End
//	help={"the commands associated with the Energy Scans package"}

	Submenu "Set Igor Data Folder"
		DataFolderListForMenus(),/Q, SetDataFolderFromMenu1()		// set to an-other data folder
	End
	help={"used to move between top level folders, also writes a note to history"}
	"Find moving axes from image files in a folder",WhatsMovingInFolder("")
	help={"identify which axes were scanned by looking at all of the image files"}
End

// ===============================================================================================================
// ===============================================================================================================
// =============================================== Start of Menus  ===============================================

Function/S DataFolderListForMenus()				// note, "!\022(" causes a menu item to be checked and disabled
	String fldrSav= GetDataFolder(1)
	SetDataFolder root:
	String list = DataFolderDir(1)
	SetDataFolder fldrSav

	list = StringByKey("FOLDERS",list)				// get list of top level folders
	list = ReplaceString(",",list,";")				// change separator from comma to semi-colnn
	list = RemoveFromList("Packages",list)		// remove "Packages"

	String fldr =  GetDataFolder(0)					// name of current folder
	if (stringmatch(fldrSav,"root:"+fldr+":"))		// we are at one of the top level data folders, so add check & disable
		Variable i = WhichListItem(fldr, list)
		list = RemoveFromList(fldr,list)
		list = AddListItem("!\022("+fldr,list,";",i)
	endif
	return list
End
Function SetDataFolderFromMenu1()
	GetLastUserMenuInfo		// sets S_value, V_value, etc.
	if (strlen(S_value)<1)
		return 1
	endif
	String from, to, fldr = "root:"+S_value+":"
	from=GetDataFolder(1)
	SetDataFolder $fldr
	to=GetDataFolder(1)
	printf "\r\r%s%s changed data folder from  '%s'  -->  '%s'\r",BULLET,BULLET,from,to
End

//Function/S MenuItemIfWaveExists(item,wName)
//	String item
//	String wName
//	return SelectString(WaveExists($wName) ,"(","")+item
//End
//
Function/S MenuItemIfGraphNamesExist(item,nameList)
	String item
	String nameList
	String name
	Variable i, N=ItemsInList(nameList)
	for (i=0;i<N;i+=1)
		name = StringFromlist(i,nameList)
		if (!strlen(WinList(name,";","WIN:1")))
			return "("+item
		endif
	endfor
	return item
End

// ================================================ End of Menus  ================================================
// ===============================================================================================================
// ===============================================================================================================




// ===============================================================================================================
// ===============================================================================================================
// ===============================================================================================================
// =========================================== Start of Q-histograms  ============================================

Function NewEnergyScan(fldr,title,d0)		// set up a folder for a new Energy Scan
	String fldr										// new folder name, it will be at the top level
	String title									// a title to be used on plots
	Variable d0										// d-spacing of unstrained material (nm)

	fldr = StringFromList(ItemsInList(fldr,":")-1, fldr,":")	// remove leading part of path
	fldr = CleanupName(fldr,0)					// igorize the name
	String fldrSav= GetDataFolder(1)
	SetDataFolder root:
	fldr = SelectString(CheckName(fldr,11),fldr,"")	// avoid existing folders
	SetDataFolder $fldrSav
	d0 = (d0>0 && numtype(d0)==0) ? d0 : NaN
	if (strlen(fldr)<1 || numtype(d0) || strlen(title)<1)	// need user input
		Prompt fldr, "new folder name (no colons)"
		Prompt d0,"reference d-spacing (nm)"
		Prompt title, "title to use for plots"
		Variable m,N=CountObjects("root:", 4 )
		String used,list="\\M1:(:Unavailable Names;"
		for (m=0;m<N;m+=1)
			list += "\\M1:(:"+GetIndexedObjName("root:",4,m)+";"
		endfor
		Prompt used, "Unavailable Names",popup,list
		DoPrompt "new scan",fldr,used,d0,title
		if (V_flag)
			return 1
		endif
		d0 = (d0>0 && numtype(d0)==0) ? d0 : NaN
		printf "\r%s|\r%s|\r",BULLET,BULLET
		printf "NewEnergyScan(\"%s\",\"%s\",%.9g)\r",fldr,title,d0
	endif
	fldr = StringFromList(ItemsInList(fldr,":")-1, fldr,":")	// remove leading part of path
	fldr = CleanupName(fldr,0)					// igorize the name
	SetDataFolder root:
	fldr = SelectString(CheckName(fldr,11),fldr,"")// avoid existing folders
	SetDataFolder $fldrSav
	if (strlen(fldr)<1)
		DoAlert 0,"the folder 'root:"+fldr+"' is invalid"
		return 1
	endif
	fldr = "root:"+fldr							// full name of new folder to create
	printf "d0 = %g,   new folder = '%s',    title='%s'\r",d0,fldr,title
	NewDataFolder/S $fldr
	printf "%s%s changed data folder from  '%s'  -->  '%s'\r",BULLET,BULLET,fldrSav,GetDataFolder(1)
	Variable/G $(fldr+":d0")=d0
	if (strlen(title))
		String/G $(fldr+":title")=title
	endif
	return 0
End



// ===============================================================================================================
// ========================= Start of Q distributions any Position or Depth, NEW Section =========================

//	Process many images all at the same depth, but in an array of x-y (possible X-H) positions. One of the positions may be depth,
//	which would occur when there was a wire scan in the process.
//	and it is assumed that the sample is thin so that the yc in geo is correct for all images.
Function Fill_Q_Positions(d0,pathName,nameFmt,range1,range2,mask,[depth,maskNorm,dark,I0normalize,dQ,printIt])	// does not assume depth in image
	Variable d0				// d-spacing of the strain=0 material (nm)
	String pathName		// either name of path to images, or the full explicit path, i.e. "Macintosh HD:Users:tischler:data:cal:recon:"
	String nameFmt			// file name format string (not path info), something like  "EW5_%d.h5", or "EW1_%d_%d.h5"
	String range1			// range of first number in file name (likely energy)
	String range2			// (OPTIONAL) range of second numbers file name (likely depth)
	Wave mask				// optional mask to limit the pixels that get processed (use pixel when mask true)
	Variable depth			// depth measured from origin, needed if images are not the result of a reconstruction (micron), default=0
	Variable maskNorm		// use pixels outside of mask to normalize whole image (suppresses pedestal)
	Wave dark				// an optional background wave
	Variable I0normalize// a Flag, if True then normalize data (default is True)
	Variable dQ				// if passed, this will be the bin size used for the histogram, if bad then a calculated value will be used
	Variable printIt
	depth = ParamIsDefault(depth) ? NaN : depth
	maskNorm = ParamIsDefault(maskNorm) ? 0 : maskNorm
	maskNorm = numtype(maskNorm) ? 0 : !(!maskNorm)
	I0normalize = ParamIsDefault(I0normalize) ? 1 : I0normalize
	I0normalize = numtype(I0normalize) ? NaN : !(!I0normalize)	// a NaN forces a prompt to user
	dQ = numtype(dQ) || dQ<=0 ? NaN : dQ				// if dQ is bad, then a calculated value will be used
	printIt = ParamIsDefault(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt

	if (!((d0>0)) && exists("d0")!=2)					// there is no d0, check with the user
		DoAlert 2, "There is no d0, so you cannot get strain, just Q.  Do you want to set d0?"
		if (V_flag==1)											// chose "Yes"
			Variable dnm = NumVarOrDefault("d0",NaN)
			dnm = numtype(dnm) ? NumVarOrDefault("root:Packages:Lattices:PanelValues:dspace_nm",NaN) : dnm
			Prompt dnm, "referencing d-spacing (nm)"
			DoPrompt "d0",dnm
			if (V_flag)
				return 1
			endif
			Variable/G $(GetDataFolder(1)+"d0")=dnm
			d0 = dnm
		elseif (V_flag==3)									// chose "Cancel"
			return 1												// quitting
		endif
	endif
	if (!((d0>0)))												// invalid d0, check for a local value
		d0 = NumVarOrDefault("d0",NaN)
		d0 = numtype(d0) ? NumVarOrDefault("root:Packages:Lattices:PanelValues:dspace_nm",NaN) : d0
	endif

	String str
	PathInfo $pathName
	if (!V_flag || strlen(nameFmt)<1)					// path does not exist or no nameFmt, ask user
		String pathPart
		str = requestPathFileRootFmt(pathName,2)	// look for at most 2 ranges
		pathPart = StringFromList(0,str)
		nameFmt = StringFromList(1,str)
		if (strlen(pathPart)<1 || strlen(nameFmt)<1)
			return 1												// invalid inputs
		endif
		if (!stringmatch(pathPart,S_path))				// path was changed
			if (stringmatch(pathName,"imagePath") || strsearch(pathName,":",0))
				NewPath/O/M="path to reconstructed image files" imagePath pathPart	// for imagePath, automatically reassign
			else
				NewPath/M="path to reconstructed image files" $pathName pathPart		// for other names, ask
			endif
		endif
		printIt = 1
	endif
	PathInfo $pathName
	if (strlen(S_path)<1 || strlen(nameFmt)<1)
		return 1													// invalid inputs
	endif
	if (printIt)
		printf "using data from files starting with '%s'\r",S_path+nameFmt
	endif

	Variable Nranges = calcNranges(nameFmt)			// number of ranges needed
	String list = WaveListClass("imageMask","*","DIMS:2,BYTE:1")	// setup to use optional mask
	String maskName = ""
	if (WaveExists(mask))
		maskName = GetWavesDataFolder(mask,2)
	elseif (strlen(list))
		maskName = ""
		Prompt maskName, "mask to use with image",popup,"-none-;"+list
		Prompt maskNorm, "subtract dark using pixels outside of mask",popup,"Nothing;Suppress"
		maskNorm += 1
		DoPrompt "mask wave",maskName,maskNorm
		if (V_flag)
			return 1
		endif
		maskNorm = maskNorm==2
		maskName = SelectString(stringmatch(maskName,"-none-"),maskName,"")
		Wave mask = $maskName								// do not check if wave exists, that a valid option
		printIt = 1
	endif
	maskName = SelectString(WaveExists(mask),"$\"\"",maskName)
	String darkList = WaveListClass("imageDark;rawImageDark","*","")
	if (!WaveExists(dark) && ItemsInList(darkList)>0)
		String darkName
		Prompt darkName,"Background Image",popup,"-none-;"+darkList
		DoPrompt "Background",darkName
		if (V_flag)
			return 1
		endif
		if (cmpstr(darkName,"-none-"))
			Wave dark = $darkName
		else
			Wave dark = $""
		endif
	endif

	String progressWin = ProgressPanelStart("",stop=1,showTime=1,status="Determining the range(s) to scan")	// display a progress bar
	Variable ask = (ItemsInRange(range1)<1 || (Nranges>1 && ItemsInRange(range2)<1))	// range is empty, need to ask
	ask = ask || numtype(I0normalize)
	if (ask) 	// if range1 is empty, get the full range from the directory
		if (ItemsInRange(range1)<1 || (Nranges>1 && ItemsInRange(range2)<1))					// range is empty, need to ask
			str = getNranges(pathName,nameFmt, progressWin=progressWin, printIt=1)
			range1 = StringFromList(0,str)
			if (Nranges>1)
				range2 = StringFromList(1,str)
			endif
		endif
		String range1Prompt = SelectString(strlen(range1)>253,range1,"(range 1 is too long to show)")
		String range2Prompt = SelectString(strlen(range2)>253,range2,"(range 2 is too long to show)")
		Prompt range1Prompt,"range1 of image file numbers to use"
		Prompt range2Prompt,"range2 of image file numbers to use"
		I0normalize = numtype(I0normalize) ? 1 : I0normalize
		I0normalize += 1
		Prompt I0normalize,"Normalize the Q-hists",popup,"Do NOT Normalize to Io or exposure;Normalize to Io & exposure"
		if (Nranges>1)
			DoPrompt "range(s)",range1Prompt,range2Prompt, I0normalize
		else
			DoPrompt "range(s)",range1Prompt, I0normalize
		endif
		if (V_flag)
			return ERROR_Fill_Q_Positions("", progressWin)
		endif
		I0normalize -= 1
		printIt = 1
		range1 = SelectString(StringMatch(range1Prompt,"(range 1 is too long to show)"),range1Prompt,range1)
		range2 = SelectString(StringMatch(range2Prompt,"(range 2 is too long to show)"),range2Prompt,range2)
	endif

	if (printIt)
		sprintf str,"Fill_Q_Positions(%g,\"%s\",\"%s\",\"%s\",\"%s\",%s",d0,pathName,nameFmt,range1,range2,maskName
		if (numtype(depth)==0)
			str += ", depth="+num2str(depth)
		endif
		if (maskNorm)
			str += ", maskNorm=1"
		endif
		if (WaveExists(dark))
			darkName = NameOfWave(dark)
			str += ", dark="+darkName
		endif
		if (!ParamIsDefault(I0normalize) || !I0normalize)
			str += ", I0normalize="+num2str(I0normalize)
		endif
		if (numtype(dQ)==0)
			str += ", dQ="+num2str(dQ)
		endif
		str += ")"
		print str[0,390]
	endif
	if (WaveExists(mask))
		if (sum(mask)==0)
			return ERROR_Fill_Q_Positions("You picked a mask that is all zero, stopping", progressWin)
		endif
	endif

	Variable N1=ItemsInRange(range1)					// number of images to be processed
	Variable N2=ItemsInRange(range2)					// number of images to be processed
	Variable N1N2=(N1<1 ? 1:N1) * (N2<1 ? 1:N2)	// total number of images
	if (N1<1 || (Nranges>1 && N2<1))
		return ERROR_Fill_Q_Positions("range is empty",progressWin)
	endif
	if (N1>0)
		printf "range1  %g  %s\r",N1,range1[0,300]
	endif
	if (N2>0)
		printf "range2  %g  %s\r",N2,range2[0,300]
	endif

	String fileFullFmt										// complete path up to and including the underscore
	PathInfo $pathName
	if (V_flag)													// if path exists, reset pathName to the explicit path, otherwise assume it is the explicit path
		fileFullFmt = S_path
	elseif (strsearch(pathName,":",Inf,1) == strlen(pathName)-1)
		fileFullFmt = pathName+":"						// ensure terminating colon
	endif
	fileFullFmt += nameFmt

	Variable useDistortion = NumVarOrDefault("root:Packages:geometry:useDistortion",USE_DISTORTION_DEFAULT)
 	if (printIt)
		printf "processing with distotion "+SelectString(useDistortion,"OFF","on")
		print "   "+SelectString(WaveExists(mask),"and no mask","and using  '"+NameOfWave(mask)+"'  for the mask")
		if (maskNorm)
			print "   Renormalizing pixels using pixels outside of the mask"
		endif
		print SelectString(I0normalize,"   NOT normalizing","   Normalizing")+" to the ion chamber and exposure time"
	endif

	// open the first file in the range, to get info about the images
	String name = fullNameFromFmt(fileFullFmt,str2num(range1),str2num(range2),NaN)
	String wnote = ReadGenericHeader(name)				// wave note of the image file
	if (strlen(wnote)<1)
		sprintf str, "could not load very first image header from '%s'\r",name
		return ERROR_Fill_Q_Positions(str,progressWin)
	endif
	STRUCT imageROIstruct roiAll							// an ROI that will contain all of the images being processed
	imageROIstructInit(roiAll, wnote=wnote)			// initialie roiAll with roi of first image
	
	Variable keV = NumberByKey("keV", wnote,"=")	// energy for this image
	if (numtype(keV))
		return ERROR_Fill_Q_Positions("invalid keV in image '"+name+"'",progressWin)
	endif

	Wave badPixelsAll = GetBadPixelsImage(wnote)	// get full image with bad pixels if badPixels exists
	if (WaveExists(badPixelsAll))
		if ( round(NumberByKey("xDimDet",wnote,"=")) != DimSize(badPixelsAll,0) || round(NumberByKey("yDimDet",wnote,"=")) != DimSize(badPixelsAll,1) )
			return ERROR_Fill_Q_Positions("badPixelsAll and image have different dimensions than the whole detector",progressWin)
		endif
	endif

	STRUCT microGeometry geo
	FillGeometryStructDefault(geo)
	Variable dNum = detectorNumFromID(geo, StringByKey("detectorID", wnote,"="))
	if (!(dNum>=0 && dNum<MAX_Ndetectors))
		return ERROR_Fill_Q_Positions("could not get detector number from from wave note using detector ID",progressWin)
	endif

	Variable Q0=2*PI/d0										// Q of unstrained material (1/nm)
	range2 = SelectString(ItemsInRange(range2),"-1",range2)

	// read header from each of the images, and store it to figure out what was done.
	// hold the positions and energies
	Make/N=(N1N2)/FREE/D X_FillQvsPositions=NaN, H_FillQvsPositions=NaN, depth_FillQvsPositions=NaN, keV_FillQvsPositions=NaN
	Make/N=(N1N2)/FREE/U/I m1_Fill_QHistAt1Depth=0, m2_Fill_QHistAt1Depth=0
	Make/N=(N1N2)/FREE/D Y_FillQvsPositions=NaN, Z_FillQvsPositions=NaN
	Variable seconds
	if (printIt)
		printf "about to examine %d image headers...   ",N1N2
	endif
	// for all the N1*N2 files (go over range), store the energies, X position, H position, and indicies
	Variable microSec0 = stopMSTimer(-2)									// used for timing, number of micro-sec
	Variable microSec = microSec0, Npixels=0
	ProgressPanelUpdate(progressWin,0,status="examining "+num2istr(N1N2)+" image headers",resetClock=1)
	Variable i,m1,m2, varyROI=0, iv
	for (m1=str2num(range1), i=0; numtype(m1)==0; m1=NextInRange(range1,m1))	// loop over range1
		for (m2=str2num(range2); numtype(m2)==0; m2=NextInRange(range2,m2))		// loop over range2
			if (mod(i,50)==0)
				if (ProgressPanelUpdate(progressWin,i/N1N2*100))		// update progress bar
					return ERROR_Fill_Q_Positions("User abort", progressWin)
				endif
			endif
			name = fullNameFromFmt(fileFullFmt,m1,m2,NaN)
			wnote = ReadGenericHeader(name,extras="EscanOnly:1")	// short wave note to add to file read in

			iv = ExtendROIfromWaveNote(roiAll,wnote)					// find ROI that holds ALL the images
			if (iv<0)
				continue																// really bad ROI in wave note, skip this one
			endif
			varyROI = iv ? 1 : varyROI										// set varyROI to true if roiAll is changing

			Npixels += NumberByKey("xdim", wnote,"=") * NumberByKey("ydim", wnote,"=")
			X_FillQvsPositions[i] = NumberByKey("X1",wnote,"=")		// list of X positioner values
			H_FillQvsPositions[i] = NumberByKey("H1",wnote,"=")		// list of H positioner values
			depth_FillQvsPositions[i] = NumberByKey("depth",wnote,"=")	// list of depths
			keV_FillQvsPositions[i] = NumberByKey("keV",wnote,"=")	// list of energies
			m1_Fill_QHistAt1Depth[i] = m1									// file id number for each image
			m2_Fill_QHistAt1Depth[i] = m2									// file id number for each image
			Y_FillQvsPositions[i] = NumberByKey("Y1",wnote,"=")		// list of Y positioner values
			Z_FillQvsPositions[i] = NumberByKey("Z1",wnote,"=")		// list of Z positioner values
			i += 1
		endfor
	endfor
	seconds = (stopMSTimer(-2)-microSec)/1e6
	if (printIt)
		printf "took %s\r",Secs2Time(seconds,5,3)
	endif
	if (varyROI && WaveExists(mask))
		return ERROR_Fill_Q_Positions("You cannot use a mask when you are also varying the ROI size", progressWin)
	elseif (WaveExists(dark) && varyROI)
		return ERROR_Fill_Q_Positions("You cannot use a dark image when you are also varying the ROI size", progressWin)
	endif
	ProgressPanelUpdate(progressWin,0,status="done with headers",resetClock=1)

	if (WaveExists(badPixelsAll))
		Wave badPixels = ExtractROIofImage(badPixelsAll, roiAll)	// badPixels is now a free wave, with size of roiAll
	else
		Wave badPixels = $""
	endif

	microSec = stopMSTimer(-2)												// for timeing bulk of processing
	Variable dx,Nx, dy,Ny, ddep,Nz, off
	FindScalingFromVec(X_FillQvsPositions,0.1,off,dx,Nx)			// get step size and number of points for the X scan
	FindScalingFromVec(H_FillQvsPositions,0.1,off,dy,Ny)			//  "    "    "   "     "   "    "     "   "  H scan
	FindScalingFromVec(depth_FillQvsPositions,0.1,off,ddep,Nz)	//  "    "    "   "     "   "    "     "   "  depths
	Variable depth0 = WaveMin(depth_FillQvsPositions)

	// let the user correct the scan ranges as necessary
	dx = abs(dx)
	dy = abs(dy)
	ddep = abs(ddep)
	Prompt dx,"step size in sample X motion (micron)"
	Prompt Nx,"# of points (NOT intervals) in Xsample"
	Prompt dy,"step size in sample H motion (micron)"
	Prompt Ny,"# of points (NOT intervals) in Hsample"
	Prompt ddep,"step size in depth (micron)"
	Prompt Nz,"# of points (NOT intervals) in depth"
	Prompt depth,"depth to assume (micron)"
	if (Nz<=1 && numtype(depth0))
		DoPrompt "scan sizes",dx,Nx,dy,Ny,depth
		depth = numtype(depth) ? 0 : depth
		Nz = 0
		depth0 = depth
	else
		DoPrompt "scan sizes",dx,Nx,dy,Ny,ddep,Nz
		ddep = (ddep==0) ? 1 : ddep
	endif
	if (V_flag)
		return 1
	endif
	dx = (dx==0) ? 1 : abs(dx)
	dy = (dy==0) ? 1 : abs(dy)
	ddep = (ddep==0) ? 1 : abs(ddep)
	Nx = numtype(Nx) || Nx<=1 ? 0 : Nx
	Ny = numtype(Ny) || Ny<=1 ? 0 : Ny
	Nz = numtype(Nz) || Nz<=1 ? 0 : Nz

	WaveStats/M=1/Q X_FillQvsPositions
	Variable X0 = V_min
	WaveStats/M=1/Q H_FillQvsPositions
	Variable H0 = V_min
	WaveStats/M=1/Q depth_FillQvsPositions
	depth0 = V_min
	depth0 = numtype(depth0) ? depth : depth0
	depth0 = numtype(depth0) ? 0 : depth0
	WaveStats/M=1/Q keV_FillQvsPositions
	if (V_numNans)
		DoAlert 0, "There were "+num2istr(V_numNans)+" bad images found, they will be skipped"
	endif
	Variable ikeVlo=V_minloc, ikeVhi=V_maxloc				// save location of min and max energies
	Variable keVmax=keV_FillQvsPositions[ikeVhi], keVmin=keV_FillQvsPositions[ikeVlo]
	Variable dkeV = FindStepSizeInVec(keV_FillQvsPositions, 0.0002)	// returns energy step, threshold=0.2 eV

	// if varyROI is True, then maskLocal will be size of roiAll
	if (WaveExists(badPixels) && WaveExists(mask))
		Duplicate/FREE mask, maskLocal							// maskLocal is combination of mask and badPixels
		maskLocal = badPixels || mask
	elseif (WaveExists(badPixels))
		Wave maskLocal = badPixels
	elseif (WaveExists(mask))
		Wave maskLocal = mask
	else
		Wave maskLocal = $""
	endif

	// note that Q range only depends upon image size and energy range, not on X, H or depth
	Variable Qmin, Qmax, NQ, dQfactor=NaN						// get range of Q
	name = fullNameFromFmt(fileFullFmt,m1_Fill_QHistAt1Depth[ikeVlo],m2_Fill_QHistAt1Depth[ikeVlo],NaN)	// load one image to get its size & Q range
	String wnoteFull = ReadGenericHeader(name)				// a full typical wave note

	Variable/C thetaZ = thetaRange(geo.d[dNum],roiAll,maskIN=maskLocal,depth=depth0)		// returns the range of theta spanned by roi
	Qmin = 4*PI * sin( real(thetaZ) ) * keVmin/hc			// min Q (1/nm)
	Qmax = 4*PI * sin( imag(thetaZ) ) * keVmax/hc			// max Q (1/nm)
	if (numtype(dQ) || dQ<=0)
		// determine dQ (1/nm), the Q resolution to use.  Base it on the distance between two adjacent pixels
		Variable px,py													// the current pixel to analyze (unbinned full chip pixels)
		px = (roiAll.xLo + roiAll.xHi)/2						// approximate full chip pixel position in center or the image (UN-binned)
		px = (roiAll.yLo + roiAll.yHi)/2
		// dQ is max of Q change in x+1, y+1, or energy step
		dQ = 4*PI*abs(sin(pixel2q(geo.d[dNum],px,py,$""))-sin(pixel2q(geo.d[dNum],px+(roiAll.binx),py,$"")))*keVmax/hc
		dQ = max(dQ, 4*PI*abs(sin(pixel2q(geo.d[dNum],px,py,$""))-sin(pixel2q(geo.d[dNum],px,py+(roiAll.biny),$"")))*keVmax/hc)
		dQ = max(dQ, 4*PI*sin(pixel2q(geo.d[dNum],px,py,$""))/hc * dkeV)	// also consider the energy step size
		dQfactor = NumVarOrDefault("root:Packages:micro:Escan:dQfactor", NaN)
		dQfactor = numtype(dQfactor) || dQfactor<=0 ? NaN : dQfactor
		dQ *= numtype(dQfactor) ? DEFAULT_dQfactor : dQfactor	// allows user to modify dQ used for histogram bins
		if (numtype(dQfactor)==0)
			wnoteFull = ReplaceNumberByKey("dQfactor",wnoteFull,dQfactor,"=")
		endif
	endif
	NQ = round((Qmax-Qmin)/dQ) + 1								// number of Qs

	if (printIt)
		if (Nx>1)
			printf "X range = [%.2f, %.2f] (%sm), dX=%.2f,  %d points\r",X0, X0+(Nx-1)*dx,Gmu, dx,Nx
		else
			printf "only one X value %g %sm\r",X0,Gmu
		endif
		if (Ny>1)
			printf "H range = [%.2f, %.2f] (%sm), dH=%.2f,  %d points\r",H0, H0+(Ny-1)*dy,Gmu, dy,Ny
		else
			printf "only one H value %g %sm\r",H0,Gmu
		endif
		if (Nz>1)
			printf "Depth range = [%.2f, %.2f] (%sm), dDepth=%.2f,  %d points\r",depth0, depth0+(Nz-1)*ddep,Gmu, ddep,Nz
		else
			printf "only one depth value\r"
		endif
		printf "E range = [%g, %g] (keV)\r", keVmin, keVmax
		if (numtype(dQfactor))
			printf "Q range = [%g, %g] (1/nm),  %sQ=%.2g(1/nm)       ",Qmin, Qmax,GDELTA,dQ
		else
			printf "Q range = [%g, %g] (1/nm),  %sQ=%.2g(1/nm), dQfactor=%g       ",Qmin, Qmax,GDELTA,dQ,dQfactor
		endif
		printf "theta range = [%g, %g]%s\r",real(thetaZ)*180/PI,imag(thetaZ)*180/PI,DEGREESIGN
		if (d0>0)
			printf "   d0 = %g (nm),   Q0 = %g (1/nm)\r",d0,Q0
		endif
	endif

	// determine scanned & un-scanned coordinates
	Make/N=4/FREE Nm={Nx,Ny,Nz,-1}, Idims={0,1,2,3}
	Sort/R Nm, Nm, Idims
	Nm = Nm<=1 ? 0 : Nm												// make unscanned dimensions 0, not 1

	// info on scanned coordinates
	Make/N=3/T/FREE LablesWAV={"Xsample", "Hsample", "Depth"}, TagsWAV={"X1", "H1", "Depth"}, UnitsWAV={Gmu+"m", Gmu+"m", Gmu+"m"}
	Make/N=3/D/FREE offsetsWAV={X0, H0, depth0}, deltasWAV={dx, dy, ddep}
	Make/N=4/T/FREE DimTags="",LabelNames="", LabelUnits=""
	Make/N=4/D/FREE del, dim0
	Variable nDim														// number of dimensions in Q_Positions [1-4]
	for (nDim=0; nDim<3 && Nm[nDim]>0; nDim+=1)
		DimTags[nDim] = TagsWAV[Idims[nDim]]
		LabelNames[nDim] = LablesWAV[Idims[nDim]]
		LabelUnits[nDim] = UnitsWAV[Idims[nDim]]
		dim0[nDim] = offsetsWAV[Idims[nDim]]
		del[nDim] = deltasWAV[Idims[nDim]]
	endfor
	Nm[nDim] = NQ														// last one is always Q
	DimTags[nDim] = "Q"
	LabelNames[nDim] = "Q"
	LabelUnits[nDim] = "1/nm"
	dim0[nDim] = Qmin
	del[nDim] = abs((Qmax-Qmin)/NQ)
	nDim += 1 

	Make/N=(Nm[0],Nm[1],Nm[2],Nm[3])/O/D Q_Positions		// array to save Q scans as function of position (X, H, depth)
	Make/N=(Nm[0],Nm[1],Nm[2],Nm[3])/FREE/D Q_PositionsNorm=0	// holds number of pixels contributing to each element in Q_Positions[], use to normalize
	Q_Positions = 0
	Make/N=(NQ)/FREE/D Qhist										// array to hold Q's from one image
	Make/N=(NQ)/FREE/D QhistNorm									// number of pixels used for each Qhist, use to normalize
	SetScale/P x dim0[0],del[0],LabelUnits[0], Q_Positions,Q_PositionsNorm
	SetScale/P y dim0[1],del[1],LabelUnits[1], Q_Positions,Q_PositionsNorm
	SetScale/P z dim0[2],del[2],LabelUnits[2], Q_Positions,Q_PositionsNorm
	SetScale/P t dim0[3],del[3],LabelUnits[3], Q_Positions,Q_PositionsNorm
	SetScale/I x Qmin,Qmax,"1/nm", Qhist

	// get list of positions
	Duplicate/FREE Nm, Ntemp
	Ntemp[nDim] = 0									// not interested in the Q's
	Ntemp = Nm[p]>=1 ? Nm[p] : 1
	Variable Ncoords = Ntemp[0]*Ntemp[1]*Ntemp[2]*Ntemp[3]	// this is actually more than will be neede, trim it later
	WaveClear Ntemp
	Make/N=(Ncoords,9)/D/FREE Q_CoordinatesBig=NaN// array to coordinates of each Qscan {XX,YY,ZZ,HH,FF,depth,i,j,k}

	seconds = (stopMSTimer(-2)-microSec)/1e6
	if (printIt && seconds>0.2)
		printf "setting up arrays took %s\r",Secs2Time(seconds,5,3)
	endif

	// convert coordinates from positioner to sample for use by Q_Coordinates
	Make/N=(i)/D/FREE xSample,ySample,zSample, depthRaw, HSample, FSample
	depthRaw = depth_FillQvsPositions[p]
	depthRaw = numtype(depthRaw) ? 0 : depthRaw
	xSample = -X_FillQvsPositions[p]
	ySample = -Y_FillQvsPositions[p]
	zSample = -Z_FillQvsPositions[p] + depthRaw[p]
	HSample = YZ2H(ySample[p],zSample[p])
	FSample = YZ2F(ySample[p],zSample[p])
	WaveClear X_FillQvsPositions, Y_FillQvsPositions, Z_FillQvsPositions, H_FillQvsPositions, depthRaw

	// done with the setup part, now actually compute something
	if (useDistortion)												// if using the distortion, precompute for all images  here
		Abort "Fill_Q_Positions(), filling the distortion map needs serious fixing"
		Wave DistortionMap = GetDistortionMap(roiAll)		// if using the distortion, precompute for all images here
	endif																	// distortion map now ready

	ProgressPanelUpdate(progressWin,0,status="making sin(theta) array",resetClock=1)
	// This sin(theta) array is the size of roiAll
//print "roiAll =",imageROIstruct2str(roiAll)
	Wave sinThetaAll = MakeSinThetaArray(roiAll,geo.d[dNum],wnote,depth=depth)	// make an array the same size as roiAll, but filled with sin(theta) for this energy
//print "sinThetaAll",NumberByKey("startx",note(sinThetaAll),"="), NumberByKey("endx",note(sinThetaAll),"="), NumberByKey("starty",note(sinThetaAll),"="), NumberByKey("endy",note(sinThetaAll),"="),"  ",DimSize(sinThetaAll,0), DimSize(sinThetaAll,1)

	STRUCT imageROIstruct ROIsinTheta
	ROIsinTheta.empty = 1											// roi of current sinTheta, starts empty to forces a calculation first time
	if (printIt)
		printf "starting the actual Q histogramming of %d images...   ",N1N2
	endif
	microSec = stopMSTimer(-2)									// timing bulk of processing
	Variable sec3=0,timer3
	ProgressPanelUpdate(progressWin,0,status="processing "+num2istr(N1N2)+" images",resetClock=1)	// update progress bar
	// for all the N1N2 files (go over range), compute Qhist for each image
	Variable ipnt, j, k, Nimage
	for (m1=str2num(range1),ipnt=0,Ncoords=0; numtype(m1)==0; m1=NextInRange(range1,m1))	// loop over range, all the images
		for (m2=str2num(range2); numtype(m2)==0; m2=NextInRange(range2,m2),ipnt+=1)					// loop over range2
			if (mod(ipnt,100)==0)
				if (ProgressPanelUpdate(progressWin,ipnt/N1N2*100))// update progress bar
					break													//   and break out of loop
				endif
			endif
			name = fullNameFromFmt(fileFullFmt,m1,m2,NaN)	// image name to process
			Wave image = $(LoadGenericImageFile(name, extras="EscanOnly:1"))	// load next image to histogram
			if (!WaveExists(image))
				printf "\r  could not load image named '%s'\r",name
				continue
			endif

			if (!imageMatchesROI(image,ROIsinTheta))		// this image has a new roi, so re-set sinTheta & maskLocalSub
//print "\r image,",NumberByKey("startx",note(image),"="), NumberByKey("endx",note(image),"="), NumberByKey("starty",note(image),"="), NumberByKey("endy",note(image),"=")
				imageROIstructInit(ROIsinTheta, wnote=note(image))	// re-set roi_Sin(theta) to match current image
//print "ROI struct = ",imageROIstruct2str(ROIsinTheta)
				Wave sinTheta = ExtractROIofImage(sinThetaAll, ROIsinTheta)
//print "sinTheta",NumberByKey("startx",note(sinTheta),"="), NumberByKey("endx",note(sinTheta),"="), NumberByKey("starty",note(sinTheta),"="), NumberByKey("endy",note(sinTheta),"=")
				Nimage = numpnts(sinTheta)
				Redimension/N=(Nimage) sinTheta
				Make/N=(Nimage)/I/FREE indexWaveQ = p
				Sort sinTheta, sinTheta,indexWaveQ				// sort so indexWaveQ[0] is index to lowest sin(theta), indexWaveQ[inf] is greatest
				Wave maskLocalSub = ExtractROIofImage(maskLocal, ROIsinTheta)
			endif

			// accumulate the Q histogram for only ONE image into Qhist
			Qhist = 0													// needed because FillQhist1image() accumulates into Qhist
			QhistNorm = 0

if (!EqualWaves(sinTheta,indexWaveQ,512) || numpnts(image)!=numpnts(sinTheta))
Debugger
endif


			wnote = FillQhist1image(image,sinTheta,indexWaveQ,Qhist,QhistNorm,maskLocalSub,dark=dark,maskNorm=maskNorm,I0normalize=I0normalize)
			KillWaves/Z image											// done with the image
			timer3=startMSTimer
			if (NumberByKey("V_min", wnote,"=")==0  && NumberByKey("V_max", wnote,"=")==0)
				sec3 += stopMSTimer(timer3)/1e6
				continue													// no intensity here, so continue
			endif
			i = round((NumberByKey(DimTags[0], wnote,"=")-dim0[0])/del[0]) ;	i = numtype(i) ? 0 : i
			j = round((NumberByKey(DimTags[1], wnote,"=")-dim0[1])/del[1]) ;	j = numtype(j) ? 0 : j
			k = round((NumberByKey(DimTags[2], wnote,"=")-dim0[2])/del[2]) ;	k = numtype(k) ? 0 : k

			if (nDim==1)
				Q_Positions[] += Qhist[p]							// add Intensity from this image
				Q_PositionsNorm[] += QhistNorm[p]				// accumulate number of pixels contibuting (not pixel intensity)
			elseif (nDim==2)
				Q_Positions[i][] += Qhist[q]						// add Intensity from this image
				Q_PositionsNorm[i][] += QhistNorm[q]			// accumulate number of pixels contibuting (not pixel intensity)
			elseif (nDim==3)
				Q_Positions[i][j][] += Qhist[r]					// add Intensity from this image
				Q_PositionsNorm[i][j][] += QhistNorm[r]		// accumulate number of pixels contibuting (not pixel intensity)
			elseif (nDim==4)
				Q_Positions[i][j][k][] += Qhist[s]				// add Intensity from this image
				Q_PositionsNorm[i][j][k][] += QhistNorm[s]	// accumulate number of pixels contibuting (not pixel intensity)
			endif
			sec3 += stopMSTimer(timer3)/1e6

			Q_CoordinatesBig[Ncoords][0] = xSample[ipnt-1]				// 0  1  2  3  4   5
			Q_CoordinatesBig[Ncoords][1] = ySample[ipnt-1]				//	XX,YY,ZZ,HH,FF,depth, i,j,k
			Q_CoordinatesBig[Ncoords][2] = zSample[ipnt-1]				// note that ipnt has already been incremented
			Q_CoordinatesBig[Ncoords][3] = HSample[ipnt-1]
			Q_CoordinatesBig[Ncoords][4] = FSample[ipnt-1]
			Q_CoordinatesBig[Ncoords][5] = depth_FillQvsPositions[ipnt-1]
			Q_CoordinatesBig[Ncoords][6] = nDim>1 ? i : -1	// use -1 for unused dimensions
			Q_CoordinatesBig[Ncoords][7] = nDim>2 ? j : -1
			Q_CoordinatesBig[Ncoords][8] = nDim>3 ? k : -1
			Ncoords += 1
		endfor
	endfor
	Redimension/N=(Ncoords,-1) Q_CoordinatesBig

	// remove all duplicates from Q_CoordinatesBig
	Make/N=(Ncoords,9)/O/D Q_Coordinates=NaN// array to coordinates of each Qscan without duplicates {XX,YY,ZZ,HH,FF,depth,i,j,k}
	Variable y0,z0
	for (j=0,ipnt=0; j<Ncoords; j+=1)
		x0 = Q_CoordinatesBig[j][0]
		y0 = Q_CoordinatesBig[j][1]
		z0 = Q_CoordinatesBig[j][2]
		if (numtype(x0+y0+z0))
			continue
		endif

		Q_Coordinates[ipnt][] = Q_CoordinatesBig[j][q]
		ipnt += 1
		for (i=j+1;i<Ncoords;i+=1)				 // remove any other points with these coordinates
			if (abs(x0-Q_CoordinatesBig[i][0])+abs(y0-Q_CoordinatesBig[i][1])+abs(z0-Q_CoordinatesBig[i][2]) < 1e-9)	// the same
				Q_CoordinatesBig[i][0] = NaN
				Q_CoordinatesBig[i][1] = NaN
				Q_CoordinatesBig[i][2] = NaN
			endif
		endfor
	endfor
	WaveClear Q_CoordinatesBig
	Ncoords = ipnt
	Redimension/N=(Ncoords,-1) Q_Coordinates

	SetDimLabel 1,0, X, Q_Coordinates			// first 6 are coordinates
	SetDimLabel 1,1, Y, Q_Coordinates
	SetDimLabel 1,2, Z, Q_Coordinates
	SetDimLabel 1,3, H, Q_Coordinates
	SetDimLabel 1,4, F, Q_Coordinates
	SetDimLabel 1,5, depth, Q_Coordinates
	SetDimLabel 1,6, i, Q_Coordinates			// last 3 are indicies into Q_Positions[i][j][k]
	SetDimLabel 1,7, j, Q_Coordinates
	SetDimLabel 1,8, k, Q_Coordinates
	SetScale d 0,0,LabelUnits[0], Q_Coordinates

	if (useDistortion)
		Note/K DistortionMap, "use=0"
	endif
	Q_Positions /= Q_PositionsNorm								// do the normalization
	Q_Positions = numtype(Q_Positions) ? NaN : Q_Positions
	seconds = (stopMSTimer(-2)-microSec)/1e6
	if (printIt)
		printf "took %s	(%.3g %ss/pixel)\r",Secs2Time(seconds,5,1),1e6*seconds/(N1N2*Npixels),Gmu
		printf "		the accumulation/assignment part took %s\r",Secs2Time(sec3,5,2)
	endif

	FindValue/TEXT="Depth"/TXOP=4/Z DimTags
	wnote = wnoteFull													// start with a full typical wave note
	if (V_value>=0)													// if Depth is an axis (meaning it varies), remove it from note
		wnote = RemoveByKey("Depth",wnote,"=")
	endif
	wnote = RemoveByKey("X2",wnote,"=")						// do not need any wire position information
	wnote = RemoveByKey("Y2",wnote,"=")
	wnote = RemoveByKey("Z2",wnote,"=")
	wnote = RemoveByKey("H2",wnote,"=")
	wnote = RemoveByKey("F2",wnote,"=")
	wnote = RemoveByKey("wirebaseX",wnote,"=")
	wnote = RemoveByKey("wirebaseY",wnote,"=")
	wnote = RemoveByKey("wirebaseZ",wnote,"=")
	wnote = RemoveByKey("AerotechH",wnote,"=")
	wnote = RemoveByKey("V_min",wnote,"=")					// the V_min & V_max are left over from the Qhist fitting
	wnote = RemoveByKey("V_max",wnote,"=")
	microSec = stopMSTimer(-2)									// timing final stuff

	if (WaveExists(mask))
		wnote = ReplaceStringByKey("maskWave",wnote,GetWavesDataFolder(mask,2),"=")
	endif
	if (numtype(Q0)==0 && Q0>0)
		wnote = ReplaceNumberByKey("Q0",wnote,Q0,"=")
	endif
	wnote = ReplaceStringByKey("DimTags",wnote,DimTags[0]+","+DimTags[1]+","+DimTags[2]+","+DimTags[3],"=")
	wnote = ReplaceNumberByKey("I0normalize",wnote,I0normalize,"=")
	String title = StrVarOrDefault("title","")
	if (strlen(title))
		title = ReplaceString("=",title,"_")					// in the wave note the string cannot have "=" or ";"
		title = ReplaceString(";",title,"_")
		wnote = ReplaceStringByKey("title",wnote,title,"=")
	endif
	Note/K Q_Positions, wnote
	CompressEmptyDimensions(Q_Positions,1)					// remove empty dimensions from Q_Positions wave (1's & 0's are empty)
	if (WaveDims(Q_Positions)==1)
		Note/K Q_Positions, ReplaceStringByKey("waveClass",note(Q_Positions),"Qhistogram","=")
	endif

	// done with processing, clean up
	DoWindow/K $progressWin										// done with status window
	seconds = (stopMSTimer(-2)-microSec)/1e6
	if (printIt)
		printf "  final stuff took %s\r",Secs2Time(seconds,5,2)
	endif
	seconds = (stopMSTimer(-2)-microSec0)/1e6
	if (printIt && seconds>2.0)
		printf "entire process took %s\r",Secs2Time(seconds,5,3)
	endif

	if (Ncoords>1)
		String cnote=note(Q_Positions)
		cnote = ReplaceStringByKey("waveClass",cnote,"CoordinatesQdist","=")
		cnote = ReplaceStringByKey("sourceWave",cnote,GetWavesDataFolder(Q_Positions,2),"=")
		Note/K Q_Coordinates, cnote
	endif

	if (printIt)
		print " "
	endif
	if (WaveDims(Q_Positions)==1)
		MakeGraph_Qhist(Q_Positions)								// just a single Q-distribution, plot it
	else
		MakeOtherQhistWaves(Q_Positions,printIt=printIt)	// make waves for plotting and put up plot(s)
	endif
	if (Ncoords>1)
		Note/K Q_Coordinates, ReplaceStringByKey("waveClass",note(Q_Positions),"CoordinatesQdist","=")
	endif
	beep
	return 0
End
//
Static Function ERROR_Fill_Q_Positions(message,progressWin)
	String message
	String progressWin
	if (strlen(message))											// if there is a message, show it
		DoAlert 0, message
		print message
	endif
	if (strlen(progressWin))										// if there is a progressWin, kill it
		DoWindow/K $progressWin
	endif
	return 1
End



Function MakeOtherQhistWaves(Q_Positions,[printIt,Qlo,Qhi])	// make waves for plotting and put up plot(s)
	Wave Q_Positions
	Variable printIt
	Variable Qlo, Qhi
	printIt = ParamIsDefault(printIt) ? NaN : printIt
	printIt = numtype(printIt) ? 0 : !(!printIt)
	Qlo = ParamIsDefault(Qlo) ? NaN : Qlo
	Qhi = ParamIsDefault(Qhi) ? NaN : Qhi

	if (!WaveExists(Q_Positions))
		String qName, wlist=WaveListClass("QdistAtPositions","*","MINCOLS:1")+WaveListClass("Qhistogram","*","")
		if (ItemsInList(wlist)==1)
			Wave Q_Positions=$StringFromList(0,wlist)
		elseif (ItemsInList(wlist)>1)
			Prompt qName,"Q_Positions wave",popup,wlist
			DoPrompt "Q_Positions",qName
			if (V_flag)
				return 1
			endif
			Wave Q_Positions=$qName
		endif
		printIt = 1
	endif
	if (!WaveExists(Q_Positions))
		return 1
	endif
	if (printIt)
		printf "MakeOtherQhistWaves(%s",NameOfWave(Q_Positions)
		if (!ParamIsDefault(printIt))
			printf ", printIt=%g",printIt
		endif
		printf ")\r"
	endif

	if (WaveDims(Q_Positions)==1)
		Note/K Q_Positions, ReplaceStringByKey("waveClass",note(Q_Positions),"Qhistogram","=")
		Wave Qhist = Q_Positions
	else
		if (ParamIsDefault(Qlo) && ParamIsDefault(Qhi))
			Wave RGB=reFitAllQdistributions(Q_Positions)		// re-fit all of the Q-distributions, and make waves for plotting
		else
			Wave RGB=reFitAllQdistributions(Q_Positions,Qlo=Qlo,Qhi=Qhi)
		endif

		// put up or refresh the plot of Qhist
		WaveStats/M=1/Q Q_Positions
		Variable imax, jmax, kmax
		if (abs(V_max)>abs(V_min))
			imax = (V_maxRowLoc-DimOffset(Q_Positions,0))/DimDelta(Q_Positions,0)
			jmax = (V_maxColLoc-DimOffset(Q_Positions,1))/DimDelta(Q_Positions,1)
			kmax = (V_maxLayerLoc-DimOffset(Q_Positions,2))/DimDelta(Q_Positions,2)
		else
			imax = (V_minRowLoc-DimOffset(Q_Positions,0))/DimDelta(Q_Positions,0)
			jmax = (V_minColLoc-DimOffset(Q_Positions,1))/DimDelta(Q_Positions,1)
			kmax = (V_minLayerLoc-DimOffset(Q_Positions,2))/DimDelta(Q_Positions,2)
		endif
		Wave Qhist = QhistFromQpositions(Q_Positions,imax,jmax,kmax)	// re-make Qhist, and fit it
	endif
	if (printIt)
		String list = fitOneQhist(Qhist,quiet=0,printIt=0)
		printf "at peak, Q = %.3f(1/nm),   fwhm = %.2g(1/nm)",NumberByKey("Qcenter", list,"="),NumberByKey("Qfwhm", list,"=")
		Variable strainPeak = NumberByKey("strain", list,"=")
		if (numtype(strainPeak)==0)
			printf ",   %sQ = %.3f,   strain = %.1e",GDELTA,NumberByKey("ÆQ", list,"="),strainPeak
		endif
		printf "\r"
		String win = MakeGraph_Qhist(Qhist)
		Variable i = MoveWinToTopRight(win,-1,-1)

		win = MakeGraph_Q_NPositioners(Q_Positions,RGB,imax,jmax,kmax)
		MoveWinToTopRight(win,-1,i-2)
	endif
	return 0
End
//
Static Function/T MakeGraph_Q_NPositioners(Q_positions,image,i0,j0,k0)	// display a Q_Positions[][][][] array, THREE positioners + Q
	Wave Q_positions									// usually Q_Positions
	Wave image											// usually Q_Positions_RGB, RGB wave associated with Q_Positions
	Variable i0,j0	,k0								// default position of cursor
	if (!WaveExists(Q_positions))
		DoAlert 0,"no 'Q_Positions' wave"
		return ""
	endif
	Variable dim = WaveDims(Q_positions)-1	// number of positioners
	if (dim>1 && !WaveExists(image))			// not just one simple Energy-Wirescan, so need an image
		DoAlert 0,"no 'Q_Positions' image wave"
		return ""
	endif
	Variable Nx=DimSize(Q_positions,0), Ny=DimSize(Q_positions,1), Nz=DimSize(Q_positions,2)
	i0 = (i0>=0 && i0<Nx) ? i0 : round(Nx/2)
	j0 = (j0>=0 && j0<Ny) ? j0 : round(Ny/2)
	k0 = (k0>=0 && k0<Nz) ? k0 : round(Nz/2)
	String graphName=CleanupName("Graph_"+GetWavesDataFolder(Q_Positions,0)+"_Q",0)
	graphName = ReplaceString("__",graphName,"_")
	String wnote=note(Q_positions)
	Variable Q0=getCurrentQ0(wnote)
	String labelX=StringByKey("labelX",wnote,"="), labelY=StringByKey("labelY",wnote,"="), labelZ=StringByKey("labelZ",wnote,"=")
	labelX += SelectString(strlen(labelX),"",FixUnitsInLabel(WaveUnits(Q_positions,0)))
	labelY += SelectString(strlen(labelY),"",FixUnitsInLabel(WaveUnits(Q_positions,1)))
	labelZ += SelectString(strlen(labelZ),"",FixUnitsInLabel(WaveUnits(Q_positions,2)))
	String title=StringByKey("title",wnote,"=")	// first try to get title from wave note
	if (strlen(title)<1)
		title = StrVarOrDefault("title","")	// nothing in wavenote, try the title global string
	endif
	title = SelectString(strlen(title),"",title+"\r\\Zr067")
	title += GetWavesDataFolder(Q_positions,1)
	if (strlen(StringByKey("dateExposed",wnote,"=")))
		title += "\r"+StringByKey("dateExposed",wnote,"=")
	endif

	if (strlen(WinList(graphName,"","WIN:"+num2istr(1+GIZMO_WIN_BIT))))	// graph or Gizmo
		DoWindow/F $graphName						// graph is already up, so just bring it to front
	elseif (exists(graphName)==5)	
		Execute graphName+"()"						// a recreation macro exists, run it

	elseif (dim==1)									// a single Energy-Wirescan
		Display/W=(292,67,728,479)/K=1			// nothing exists, create the graph
		DoWindow/C $graphName
		AppendImage Q_Positions
		ModifyImage Q_Positions ctab= {*,*,Terrain,1}
		ModifyGraph tick=2, mirror=1, minor=1, lowTrip=0.001, standoff=0, axOffset(left)=-1.5
		Label bottom labelX
		Label left labelY
		TextBox/C/N=textTitle/F=0 title
		Cursor/P/I/L=1/H=2 A Q_Positions i0,j0
		ShowInfo
//		Button SumQhistOverPositionsButton,pos={100,0},size={90,16},fSize=9,proc=EnergyScans#SumQhistOverPositionsButtonProc,title="Sum over Positions"
		Button SumQhistOverPositionsButton,pos={0,0},size={90,16},fSize=9,proc=EnergyScans#SumQhistOverPositionsButtonProc,title="Sum over Positions"
		if (Q0>0)
			SetDrawLayer UserFront					// draw line at Q0
			SetDrawEnv xcoord=prel,ycoord=left,dash=1
			DrawLine 0,Q0,1,Q0
		endif
		SetWindow kwTopWin hook(QfromCursor)=EnergyScans#QfromCursorHook	// this forces a re-calculation of Qhist when cursor A is moved

	elseif (dim==2)									// Energy-Wirescans at multiple positions
		Display/W=(292,67,728,479)/K=1			// nothing exists, create the graph
		DoWindow/C $graphName
		AppendImage image								// likely an RGB image
		if (numtype(Q0)==0 && DimSize(image,2)!=3)
			WaveStats/M=1/Q Q_positions
			Variable zrange = max(abs(V_min),abs(V_max))
			ModifyImage $NameOfWave(image) ctab= {-zrange,zrange,RedWhiteBlue,0}
		elseif (DimSize(image,2)!=3)
			ModifyImage $NameOfWave(image) ctab= {*,*,Terrain,1}
		endif
		ModifyGraph tick=0, mirror=1, minor=1, lowTrip=0.001, standoff=0, axOffset(left)=-1.5
		Label bottom labelX
		Label left labelY
		TextBox/C/N=textTitle/F=0 title
		Cursor/P/I/L=1/H=1 A $NameOfWave(image) i0,j0
		// ShowInfo
//		ImageStats/M=1 Q_positions
//		Cursor/P/I A $NameOfWave(image) V_maxRowLoc,V_maxColLoc
		WaveStats/M=1/Q Q_positions
		Cursor/I A $NameOfWave(image) V_maxRowLoc,V_maxColLoc
		SetAspectToSquarePixels(graphName)
//		Button SumQhistOverPositionsButton,pos={100,0},size={90,16},fSize=9,proc=EnergyScans#SumQhistOverPositionsButtonProc,title="Sum over Positions"
		Button SumQhistOverPositionsButton,pos={0,0},size={90,16},fSize=9,proc=EnergyScans#SumQhistOverPositionsButtonProc,title="Sum over Positions"
		SetWindow kwTopWin hook(QfromCursor)=EnergyScans#QfromCursorHook	// this forces a re-calculation of Qhist when cursor A is moved

	elseif (dim>2)
		print "Do not yet know how to display 3-Positioner Q-distributions"
		Abort "Do not yet know how to display 3-Positioner Q-distributions"
	endif

	return graphName
End
//
Static Function SumQhistOverPositionsButtonProc(B_Struct) : ButtonControl
	STRUCT WMButtonAction &B_Struct

	if (B_Struct.eventCode != 2)			// only process button-up
		return 0
	endif
	if (stringmatch(B_Struct.ctrlName,"SumQhistOverPositionsButton"))
		String wName = StringFromList(0,ImageNameList(B_Struct.win,";"))
		Wave Q_pos = ImageNameToWaveRef(B_Struct.win, wName )
	endif
	if (!WaveExists(Q_pos))
		return 1
	endif
	Q_SumAllPositions(Q_pos)
	return 0
End

Function/WAVE Q_SumAllPositions(Qwav)			// sum over all the positions in a Q vs Position (including Depth) wave
	Wave Qwav

	String list = WaveListClass("Qhistogram,QdistAtPositions","*","MINCOLS:2")
	if (WaveExists(Qwav))
		String wname = NameOfWave(Qwav)
	elseif (ItemsInList(list)==1)
		Wave Qwav = $StringFromList(0,list)
	elseif (ItemsInList(list)>1)					// need to ask
		Prompt wname,"wave of Q vs Position",popup,list
		DoPrompt "Q vs Positioner",wname
		if (V_flag)
			return $""
		endif
		Wave Qwav = $wname
	endif
	if (!WaveExists(Qwav))
		return $""
	endif

	String wnote=note(Qwav)
	String labelX="Q"
	Variable nDims=WaveDims(Qwav)
	Variable Nx=DimSize(Qwav,0), Ny=DimSize(Qwav,1), Nz=DimSize(Qwav,2), Nw=DimSize(Qwav,3), i, iQ
	if (nDims<=1)
		return Qwav
	elseif (nDims==2)					// remember, last dimension is Q
		iQ = 1
		Duplicate/FREE Qwav, sumTemp1a
		sumTemp1a = numtype(Qwav) ? 0 : Qwav
		MatrixOp/FREE sumTemp1 = sumCols(sumTemp1a)
		Make/N=(Ny)/D/FREE sumTemp=sumTemp1
		WaveClear sumTemp1, sumTemp1a
		labelX = StringByKey("labelY",wnote,"=")
		wnote = RemoveByKey("labelY",wnote,"=")
		wnote = ReplaceStringByKey("labelX",wnote,labelX,"=")
	elseif (nDims==3)
		iQ = 2
		Make/N=(Nx,Ny)/FREE/D sumTemp2 = 0
		Make/N=(Nz)/D/FREE sumTemp=0
		for (i=0;i<Nz;i+=1)
			sumTemp2 = Qwav[p][q][i]
			sumTemp2 = numtype(sumTemp2) ? 0 : sumTemp2
			sumTemp[i] = sum(sumTemp2)
		endfor
		WaveClear sumTemp2
		labelX = StringByKey("labelZ",wnote,"=")
		wnote = RemoveByKey("labelZ",wnote,"=")
		wnote = RemoveByKey("labelY",wnote,"=")
		wnote = ReplaceStringByKey("labelX",wnote,labelX,"=")
	elseif (nDims==4)
		iQ = 3
		Make/N=(Nx,Ny,Nz)/FREE/D sumTemp3 = 0
		Make/N=(Nw)/D/FREE sumTemp=0
		for (i=0;i<Nw;i+=1)
			sumTemp3 = Qwav[p][q][r][i]
			sumTemp3 = numtype(sumTemp3) ? 0 : sumTemp3
			sumTemp[i] = sum(sumTemp3)
		endfor
		WaveClear sumTemp3
		labelX = StringByKey("labelW",wnote,"=")
		wnote = RemoveByKey("labelW",wnote,"=")
		wnote = RemoveByKey("labelZ",wnote,"=")
		wnote = RemoveByKey("labelY",wnote,"=")
		wnote = ReplaceStringByKey("labelX",wnote,labelX,"=")
	endif

	wname = NameOfWave(Qwav)+"_SumPositions"
	Duplicate/O sumTemp, $wname
	WaveClear sumTemp
	Wave ww=$wname
	SetScale/P x, DimOffset(Qwav,iQ),DimDelta(Qwav,iQ),WaveUnits(Qwav,iQ),ww
	String class = StringFromList(0,StringByKey("waveClass",wnote,"="),",")+"_SumPos"
	wnote = ReplaceStringByKey("waveClass",wnote,class,"=")
	Note/K ww,wnote

	String graphName=StringFromList(0,WindowsWithWave(ww,1))
	if (strlen(graphName))									// bring graph with ww to the front
		DoWindow/F $graphName
		return ww
	endif

	// graph the wave
	graphName=CleanupName("Graph_"+GetWavesDataFolder(Qwav,1)+"_Q_SumPos",0)
	graphName = ReplaceString("__",graphName,"_")
//	Display/W=(732,69,1262,429)/K=1 ww
	Display/W=(609,353,1118,616)/K=1 ww
	DoWindow/C $graphName
	ModifyGraph gfMult=130, tick=2, mirror=1, minor=1, lowTrip=0.001, standoff=0
	ModifyGraph rgb=(1,16384,41928), marker=19,msize=3,lsize=1, mode=4
	String labelValue = StringByKey("labelValue",wnote,"=")
	labelValue = SelectString(strlen(labelValue),"Intensity",labelValue)	// default label for Y-axis is "Intensity"
	Label left labelValue+FixUnitsInLabel(WaveUnits(ww,-1))
//	String labelX = StringByKey("labelX",wnote,"=")
	labelX = SelectString(strlen(labelX),"Q",labelX)			// default label for X-axis is "Q"
	Label bottom labelX+FixUnitsInLabel(WaveUnits(ww,0))
	SetAxis/N=1 left 0,*
	ShowInfo
	String str = StringByKey("title",wnote,"=")// first try to get title from wave note
	if (strlen(str)<1)
		str = StrVarOrDefault("title","")	// nothing in wavenote, try the title global string
	endif
	str = SelectString(strlen(str),"",str+"\r\\Zr075")
	str += "\\JR"
	str += GetWavesDataFolder(ww,1)
	if (strlen(StringByKey("dateExposed",wnote,"=")))
		str += "\r"+StringByKey("dateExposed",wnote,"=")
	endif
	TextBox/A=RT/C/N=textTitle/F=0 str

	Variable Q0=getCurrentQ0(wnote)
	if (Q0>0)
		SetDrawLayer UserFront
		SetDrawEnv xcoord= bottom,ycoord= prel,dash= 1
		DrawLine Q0,0,Q0,1
	endif
	return ww
End

Function MakeLayoutQ_Positions_hist(prefix)
	String prefix
	String gName="",gName2=""
	if (strlen(prefix)<1)
		String wlist=WinList("Graph*_Q",";","WIN:1")
		if (ItemsInList(wlist)<1)
			return 1
		elseif (ItemsInList(wlist)==1)
			gName = StringFromList(0,wlist)		
		else
			Prompt gName, "name of graph for layout",popup,WinList("Graph*_Q",";","WIN:1")
			DoPrompt "Pick Graph",gName
			if (V_flag)
				return 1
			endif
		endif
		Variable i
		i = strsearch(gName, "Graph_",0)
		if (i!=0)
			return 1
		endif
		gName = gName[6,Inf]
		i = strsearch(gName,"_Q",Inf,3)
		prefix = gName[0,i-1]
	endif
	gName = "Graph_"+prefix+"_Q"
	gName2 = "Graph_"+prefix+"_Qhist"
	if (WinType(gName)!=1)
		return 1
	endif
	NewLayout/C=1/W=(29,67,530,580)/P=portrait
	AppendLayoutObject/R=(53,23,553,503) graph $gName
	if (WinType(gName2)==1)
		AppendLayoutObject/R=(65,505,555,739) graph  $gName2
	endif
	TextBox/N=stamp0/F=0/A=RB/X=0.17/Y=0.14 "\\Z06\\{\"%s %s\",date(), time()}"
	TextBox/N=stamp1/F=0/A=LB/X=0.17/Y=0.14 "\\Z06\\{\"%s\",CornerStamp1_()}:"+gName
End


Static Function/T MakeGraph_Q_RandomPositions(Q_Coordinates,ix,iy,[rgb,printIt])	// display a Q_Coordinates[][9] array, THREE positioners + Q
	Wave Q_Coordinates								// usually Q_Coordinates
	Variable ix,iy										// columns for x and y, in range [0,5]
	Wave rgb												// the three-column rgb wave
	Variable printIt
	if (ParamIsDefault(rgb))
		Wave rgb = $""
	endif
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : !(!printIt)

	if (!WaveExists(Q_Coordinates))
		String Q_C_list=WaveListClass("CoordinatesQdist","*","DIMS:2,MINCOLS:9,MAXCOLS:9")
		if (ItemsInLIst(Q_C_list)<1)
			return ""
		elseif (ItemsInLIst(Q_C_list)==1)
			Wave Q_Coordinates = $StringFromList(0,Q_C_list)
		else
			String name
			Prompt name,"Q_Coordinates", popup, Q_C_list
			DoPrompt "Q_Coordinates", name
			if (V_flag)
				return ""
			endif
			Wave Q_Coordinates = $name
			printIt = 1
		endif
	endif
	if (!WaveExists(Q_Coordinates))
		DoAlert 0,"no 'Q_Coordinates' wave"
		return ""
	endif
	Variable Ncoord=DimSize(Q_Coordinates,0)
	// 0  1  2  3  4   5
	//	XX,YY,ZZ,HH,FF,depth, i,j,k

	String xName, yName, allAxes="X;Y;Z;H;F;depth"
	if ( !(ix==limit(round(ix),0,5) && iy==limit(round(iy),0,5) && ix != iy) )
		Make/N=(Ncoord)/FREE/D vec
		String axes=""
		vec = Q_Coordinates[p][0]
		WaveStats/M=1/Q vec
		axes += SelectString(V_max-V_min > 0.05, "", "X;")
		vec = Q_Coordinates[p][1]
		WaveStats/M=1/Q vec
		axes += SelectString(V_max-V_min > 0.05, "", "Y;")
		vec = Q_Coordinates[p][2]
		WaveStats/M=1/Q vec
		axes += SelectString(V_max-V_min > 0.05, "", "Z;")
		vec = Q_Coordinates[p][3]
		WaveStats/M=1/Q vec
		axes += SelectString(V_max-V_min > 0.05, "", "H;")
		vec = Q_Coordinates[p][4]
		WaveStats/M=1/Q vec
		axes += SelectString(V_max-V_min > 0.05, "", "F;")
		vec = Q_Coordinates[p][5]
		WaveStats/M=1/Q vec
		axes += SelectString(V_max-V_min > 0.05, "", "depth;")

		xName = StringFromList(0,axes)
		xName = StringFromList(1,axes)
		Prompt xName, "X-axis", popup, axes
		Prompt yName, "Y-axis", popup, axes
		DoPrompt "Axes", xName, yName
		if (V_flag)
			return ""
		endif
		ix = WhichListItem(xName,allAxes)
		iy = WhichListItem(yName,allAxes)
		printIt = 1
	endif
	if ( !(ix==limit(round(ix),0,5) && iy==limit(round(iy),0,5) && ix != iy) )
		return ""
	endif
	xName = StringFromList(ix,allAxes)
	yName = StringFromList(iy,allAxes)
	String graphName=CleanupName("Graph_"+GetWavesDataFolder(Q_Coordinates,0)+"_"+xName+yName,0)
	graphName = ReplaceString("__",graphName,"_")

	if (printIt)
		printf "MakeGraph_Q_RandomPositions(%s, %g,%g",NameOfWave(Q_Coordinates),ix,iy
		if (!ParamIsDefault(printIt))
			printf ", printIt=%g", printIt
		endif
		printf ")\t\t// creating %s\r",graphName
	endif

	String wnote=note(Q_Coordinates)
	Variable Q0=EnergyScans#getCurrentQ0(wnote)

	String labelX = GetDimLabel(Q_Coordinates,1,ix)
	String labelY = GetDimLabel(Q_Coordinates,1,iy)
	labelX += SelectString(strlen(labelX),"",FixUnitsInLabel(WaveUnits(Q_Coordinates,-1)))
	labelY += SelectString(strlen(labelY),"",FixUnitsInLabel(WaveUnits(Q_Coordinates,-1)))

	String title=StringByKey("title",wnote,"=")	// first try to get title from wave note
	if (strlen(title)<1)
		title = StrVarOrDefault("title","")	// nothing in wavenote, try the title global string
	endif
	title = SelectString(strlen(title),"",title+"\r\\Zr067")
	title += GetWavesDataFolder(Q_Coordinates,1)
	if (strlen(StringByKey("dateExposed",wnote,"=")))
		title += "\r"+StringByKey("dateExposed",wnote,"=")
	endif

	if (strlen(WinList(graphName,"","WIN:1")))	// graph
		DoWindow/F $graphName						// graph is already up, so just bring it to front
	elseif (exists(graphName)==5)	
		Execute graphName+"()"						// a recreation macro exists, run it
	else													// Energy-Wirescans at multiple positions, find rgb if not passed
		if (!WaveExists(rgb))						// find rgb if not passed
			String list=WaveListClass("Q_Coordinates_RGB","*","DIMS:2,MINCOLS:3,MAXCOLS:3")
			if (ItemsInLIst(list)<1)
				Wave rgb = $""
			elseif (ItemsInLIst(list)==1)
				Wave rgb = $StringFromList(0,list)
			else
				Prompt name,"Q_Coordinates_RGB wave", popup, list
				DoPrompt "rgb", name
				if (V_flag)
					return ""
				endif
				Wave rgb = $name
			endif
		endif
		if (!WaveExists(rgb))
			print "    Unable to find three-column Q_Coordinates_RGB wave in ",GetWavesDataFolder(Q_Coordinates,1)
			print "    You will have to make it to see anything useful, use the 'Recalc RGB for Q-Distribution' button"
			print "    Use the 'Recalc RGB for Q-Distribution' button"
			print " "
		endif
		Display/W=(292,67,728,479)/K=1			// nothing exists, create the graph
		DoWindow/C $graphName
		AppendToGraph Q_Coordinates[*][iy] vs Q_Coordinates[*][ix]
		ModifyGraph mode=3, marker=18, tick=2, mirror=1, minor=1, lowTrip=0.001, standoff=1
		if (WaveExists(rgb))
			name = NameOfWave(Q_Coordinates)
			ModifyGraph zColor($name)={rgb,*,*,directRGB}
		endif
		Label bottom labelX
		Label left labelY
		TextBox/C/N=textTitle/A=RT/X=2/Y=2/F=0 title
		SetWindow kwTopWin, hook(QfromCursor)=EnergyScans#QfromCursorHook	// this forces a re-calculation of Qhist when cursor A is moved
		SetAxis/A left
		if (ix==2)
			SetAxis/A/R bottom
		else
			SetAxis/A bottom
		endif
		Cursor/P/L=1/H=1 A Q_Coordinates round(Ncoord/2)
		ShowInfo
		DoUpdate/W=$graphName
		SetAspectToSquarePixels(graphName)
	endif
	return graphName
End


// This is normally only called by MakeOtherQhistWaves()
Static Function/WAVE reFitAllQdistributions(Q_Positions,[Qlo,Qhi])// re-fit all of the Q-distributions, and make waves for plotting
	Wave Q_Positions
	Variable Qlo, Qhi				// scales color of _dQ wave, or if no _dQ, scales _Qc
	Qlo = ParamIsDefault(Qlo) ? NaN : Qlo
	Qlo = numtype(Qlo) ? NaN : Qlo
	Qhi = ParamIsDefault(Qhi) ? NaN : Qhi
	Qhi = numtype(Qhi) ? NaN : Qhi
	if (!WaveExists(Q_Positions))
		return $""
	endif

	Variable nDim = WaveDims(Q_Positions)						// number of dimensions
	String wnote = note(Q_Positions)
	Variable Q0=getCurrentQ0(wnote)
	Q0 = Q0>0 ? Q0 : 0
	Variable I0normalize = NumberByKey("I0normalize",wnote,"=")
	I0normalize = numtype(I0normalize) ? 0 : !(!I0normalize)

	Make/N=(4)/D/FREE Nm
	Nm = DimSize(Q_Positions,p)
	Make/N=4/T/FREE DimTags
	String str=StringByKey("DimTags",wnote,"=")
	DimTags = StringFromList(p,str,",")
	// all information accumulated, make auxiallary waves suitable for display

	Nm[nDim-1] = 0					// dimension index for the Q, the last one, zero out the Q dimension
	// create the waves for displaying
	String name=NameOfWave(Q_Positions)
	Make/N=(Nm[0],Nm[1],Nm[2])/O $(name+"_Qwidth")/WAVE=Qwidth =NaN	// holds the image to plot, width of Q peak (1/nm)
	Make/N=(Nm[0],Nm[1],Nm[2])/O $(name+"_Intens")/WAVE=Intens =NaN	// holds the image to plot, max intensity
	Make/N=(Nm[0],Nm[1],Nm[2])/O $(name+"_Qcenter")/WAVE=Qc =NaN		// holds the image to plot, center of Q peak (1/nm)
	Make/N=(Nm[0],Nm[1],Nm[2])/O $(name+"_QpkArea")/WAVE=QpkArea =NaN	// holds the image to plot, center of Q peak (1/nm)
	Make/N=(Nm[0],Nm[1],Nm[2])/O $(name+"_QcenterErr")/WAVE=QcErr =NaN// errors
	Make/N=(Nm[0],Nm[1],Nm[2])/O $(name+"_QwidthErr")/WAVE=QwidthErr =NaN
	Make/N=(Nm[0],Nm[1],Nm[2])/O $(name+"_Qbkg")/WAVE=Qbkg =NaN
	Make/N=(Nm[0],Nm[1],Nm[2])/O $(name+"_QbkgErr")/WAVE=QbkgErr =NaN
	if (numtype(Q0)==0)
		Make/N=(Nm[0],Nm[1],Nm[2])/O $(name+"_dQ")/WAVE=dQ =NaN			// holds the image to plot, delta peak position (1/nm)
		Make/N=(Nm[0],Nm[1],Nm[2])/O $(name+"strain")/WAVE=strain =NaN		// strain is available
		Make/N=(Nm[0],Nm[1],Nm[2])/O $(name+"strainErr")/WAVE=strainErr =NaN
	endif
	Nm[nDim-1] = 0
	Variable Nx=Nm[0], Ny=Nm[1], Nz=Nm[2]
	String list
	Variable quiet=0												// first time through, not quiet
	Variable i,j,k, center
	for (k=0;k<(Nz<1 ? 1 : Nz);k+=1)						// execute z-loop at least once
		for (j=0;j<(Ny<1 ? 1 : Ny);j+=1)					// execute y-loop at least once
			for (i=0;i<(Nx<1 ? 1 : Nx);i+=1)				// execute x-loop at least once
				Wave Qhist = QhistFromQpositions(Q_Positions,i,j,k)		// from a Q_positions type array, re-make Qhist, and fit it
				list = fitOneQhist(Qhist,quiet=quiet,printIt=0)
				quiet = 1
				center = NumberByKey("Qcenter",list,"=")
				Qc[i][j][k] = center											// Q of peak
				Qwidth[i][j][k] = NumberByKey("Qfwhm",list,"=")		// Q width
				Intens[i][j][k] = WaveMax(Qhist)							// max intensity
				QpkArea[i][j][k] = NumberByKey("QpkArea",list,"=")	// area peak
				Qbkg[i][j][k] = NumberByKey("QpkBkg",list,"=")		// bkg of the Lorentzian
				QwidthErr[i][j][k] = NumberByKey("fwhmErr",list,"=")// the errors
				QcErr[i][j][k] = NumberByKey("QcErr",list,"=")
				QbkgErr[i][j][k] = NumberByKey("QpkBkgErr",list,"=")
			endfor
		endfor
	endfor
	CopyScales/I Q_Positions, Intens, Qc,QcErr, Qwidth,QwidthErr, QpkArea, Qbkg, QbkgErr
	if (numtype(Q0)==0)
		dQ = Qc-Q0													// delta Q
		strain = -dQ/Q0
		strainErr = QcErr/Q0
		CopyScales/I Q_Positions, dQ, strain, strainErr
		WaveStats/M=1/Q Qc
		if (Qhi<=Qlo || numtype(Qlo+Qhi))					// figure out ÆQrange if valid number not passed
			Qhi = max(abs(V_min),abs(V_max))				// used to scale the color of _dQ waves, a symmetric Q range
			Qlo = -Qhi
		endif
		SetScale d Qlo,Qhi,"1/nm", dQ
	else
		WaveStats/M=1/Q Qc
		Qlo = V_min													// RGB scaling uses Qc, set range of color scale
		Qhi = V_max
		SetScale d Qlo,Qhi,"1/nm", Qc
	endif

	WaveStats/M=1/Q Intens
	Variable ImaxPnt=V_maxloc, maxIntens=V_max
	Variable XofMax=NaN, HofMax=NaN, DepthofMax=NaN
	Make/N=3/D/FREE maxLocs={V_maxRowLoc,V_maxColLoc,V_maxLayerLoc}
	Make/N=3/T/FREE prefix={"X","Y","Z","W"}
	for (i=0;i<(nDim-1);i+=1)
		if (StringMatch(DimTags[i],"X1"))
			XofMax = maxLocs[i]
			wnote = ReplaceStringByKey("label"+prefix[i],wnote,"X","=")
		elseif (StringMatch(DimTags[i],"H1"))
			HofMax = maxLocs[i]
			wnote = ReplaceStringByKey("label"+prefix[i],wnote,"H","=")
		elseif (StringMatch(DimTags[i],"Depth"))
			DepthofMax = maxLocs[i]
			wnote = ReplaceStringByKey("label"+prefix[i],wnote,"Depth","=")
		endif
	endfor
	wnote = ReplaceStringByKey("label"+prefix[nDim-1],wnote,"Q","=")

	wnote = ReplaceStringByKey("waveClass",wnote,"Q_Positions","=")
	wnote = ReplaceStringByKey("sourceWave",wnote,GetWavesDataFolder(Q_Positions,2),"=")
	wnote = ReplaceNumberByKey("ImaxPnt",wnote,ImaxPnt,"=")
	wnote = ReplaceNumberByKey("maxIntens",wnote,maxIntens,"=")
	wnote = ReplaceNumberByKey("XofMax",wnote,XofMax,"=")		// X of the max point (micron)
	wnote = ReplaceNumberByKey("HofMax",wnote,HofMax,"=")		// H of the max point (micron)
	if (numtype(DepthofMax)==0)
		wnote = ReplaceNumberByKey("DepthofMax",wnote,DepthofMax,"=") // depth of the max point (micron)
	endif
	str = SelectString(I0normalize,"Intensity","Intensity / I\\B0\\M")
	wnote = ReplaceStringByKey("labelValue",wnote,str,"=")
	Note/K Q_Positions, ReplaceStringByKey("waveClass",wnote,"QdistAtPositions","=")
	Note/K Intens, wnote
	Note/K Qc, wnote
	Note/K Qwidth, wnote
	Note/K QpkArea, wnote
	Note/K Qbkg, wnote
	Note/K QcErr, wnote
	Note/K QwidthErr, wnote
	Note/K QbkgErr, wnote
	if (numtype(Q0)==0)
		Note/K dQ, wnote
		Note/K strain, wnote
		Note/K strainErr, wnote
		Wave RGB=MakeRGBforQdistribution(dQ,Intens,Qlo=Qlo,Qhi=Qhi,printIt=1)	// recalc the RGB for fancy plots
	else
		Wave RGB=MakeRGBforQdistribution(Qc,Intens,Qlo=Qlo,Qhi=Qhi,printIt=1)
	endif
	return RGB
End
//
Static Function getCurrentQ0(wnote)		// get current Q0 (1/nm), first try from global Variable d0, then wave note
	String wnote

	Variable d0=NumVarOrDefault("d0",NaN), Q0
	d0 = d0<=0 || numtype(d0) ? NumVarOrDefault("root:Packages:Lattices:PanelValues:dspace_nm",NaN) : d0
	Q0 = 2*PI/d0
	if (Q0<=0 || numtype(Q0))
		Q0 = NumberByKey("Q0",wnote,"=")
		if (Q0<=0 || numtype(Q0))
			Q0 = 2*PI/NumberByKey("d0",wnote,"=")
		endif
	endif
	Q0 = numtype(Q0) || Q0<=0 ? NaN : Q0
	return Q0
End



Function/WAVE MakeRGBforQdistribution(dQ,IntensIN,[Qlo,Qhi,power,printIt])				// this is for user input
	Wave dQ, IntensIN
	Variable Qlo,Qhi				// used to scale the color of _Qc wave, uses a symmetric Q range
	Variable power					// power for scaling intensity, e.g. using IntensIN^power
	Variable printIt
	Qlo = ParamIsDefault(Qlo) ? NaN : Qlo
	Qlo = numtype(Qlo) ? NaN : Qlo
	Qhi = ParamIsDefault(Qhi) ? NaN : Qhi
	Qhi = numtype(Qhi) ? NaN : Qhi
	power = ParamIsDefault(power) ? 1 : power
	power = numtype(power) ? 1 : power
	printIt = ParamIsDefault(printIt) ? NaN : printIt
	printIt = numtype(printIt) ? 0 : !(!printIt)

	String wlist=WaveListClass("Q_positions","*","")
	if (!WaveExists(dQ) || !WaveExists(IntensIN))
		String sdQ=NameOfWave(dQ), sIntens=NameOfWave(IntensIN)
		sdQ = SelectString(strlen(sdQ),"Q_Positions_dQ",sdQ)
		sIntens = SelectString(strlen(sIntens),"Q_Positions_Intens",sIntens)
		Prompt sdQ,"dQ wave (for color)",popup,wlist
		Prompt sIntens,"Intensity wave (for blackness)",popup,wlist+";-none-"
		Prompt Qlo, "Lower scaling range of Q (use NaN for default)"
		Prompt Qhi, "High scaling range of Q (use NaN for default)"
		Prompt power, "power for scaling intensity (Intens^power)"
		DoPrompt "Waves for RGB", sdQ,sIntens,Qlo,Qhi,power
		if (V_flag)
			return $""
		endif
		Wave dQ=$sdQ, IntensIN=$sIntens
		printIt = 1
	endif
	if (!WaveExists(dQ))										//	if (!WaveExists(dQ) || !WaveExists(IntensIN))
		return $""
	endif
	Variable nDim=WaveDims(dQ)								// number of dimensions
	if (nDim<1 || nDim>=4)										// want to make the RGB wave
		return $""
	elseif (Qhi<=Qlo || numtype(Qlo+Qhi))					// figure out Qrange if valid numbers not passed
		WaveStats/M=1/Q dQ
		if (StringMatch(NameOfWave(dQ),"*_dQ"))			// for dQ, prefer symmetric range
			Qhi = max(abs(V_min),abs(V_max))				// used to scale the color of _dQ waves
			Qlo = -Qhi												// need a symmetric Q range
		else
			Qhi = V_max												// an asymmetric Q range
			Qlo = V_min
		endif
	endif
	power = numtype(power) ? 1 : power
	if (printIt)
		if (WaveExists(IntensIN))
			sIntens = NameOfWave(IntensIN)
		else
			sIntens = "$\"\""
		endif
		printf "MakeRGBforQdistribution(%s, %s, Qlo=%g, Qhi=%g",NameOfWave(dQ),sIntens,Qlo, Qhi
		if (!ParamIsDefault(printIt))
			printf ", printIt=%g",printIt
		endif
		if (power!=1)
			printf ", power=%g",power
		endif
		printf ")\r"
	endif

	Make/N=(4)/D/FREE Nm
	Nm = DimSize(dQ,p)
	Nm[nDim] = 3													// dimension index for the RGB, the last one, zero out the Q dimension
	String sourceName = StringByKey("sourceWave",note(dQ),"=")
	String name = NameOfWave($sourceName)
	Make/N=(Nm[0],Nm[1],Nm[2],Nm[3])/O/W/U $(name+"_RGB")/WAVE=RGB =0	// holds the image to plot, color
	Nm[nDim] = 0
	Variable cMax=65535, useIntensity=WaveExists(IntensIN)
	if (useIntensity)
		Duplicate/FREE IntensIN, Intens
		Intens = abs(IntensIN)^power							// scaled intensity to use
		WaveStats/M=1/Q Intens
		Variable maxIntens=V_max, minIntens=V_min
	endif
	// set the RGB wave based on dQ and Intensity
	//	minIntens = 0
	Variable i,j,k, Nx=Nm[0], Ny=Nm[1], Nz=Nm[2]
	Variable red,green,blue, maxc, iScaling=1
	for (k=0;k<(Nz<1 ? 1 : Nz);k+=1)						// execute z-loop at least once
		for (j=0;j<(Ny<1 ? 1 : Ny);j+=1)					// execute y-loop at least once
			for (i=0;i<(Nx<1 ? 1 : Nx);i+=1)				// execute x-loop at least once
				blue = (dQ[i][j][k] - Qlo)/(Qhi-Qlo)		// set the hue, was blue = (dQ[i][j][k] + Qrange)/(2*Qrange)
				red = 1-blue
				green= min(red,blue)
				maxC = max(max(red,green),blue)				// max RGB value
				if (useIntensity)
					iScaling = (Intens[i][j][k]-minIntens)/maxIntens// intensity scaling of this point
				endif
				red = cMax * red/maxC * iScaling			// set to saturated values scaled down by iScaling
				green = cMax * green/maxC * iScaling
				blue = cMax * blue/maxC * iScaling
				if (nDim==3)
					RGB[i][j][k][0] = limit(red,0,cMax)	// red
					RGB[i][j][k][1] = limit(green,0,cMax)	// green
					RGB[i][j][k][2] = limit(blue,0,cMax)	// blue
				elseif (nDim==2)
					RGB[i][j][0] = limit(red,0,cMax)
					RGB[i][j][1] = limit(green,0,cMax)
					RGB[i][j][2] = limit(blue,0,cMax)
				elseif (nDim==1)
					RGB[i][0] = limit(red,0,cMax)
					RGB[i][1] = limit(green,0,cMax)
					RGB[i][2] = limit(blue,0,cMax)
				endif
			endfor
		endfor
	endfor

	String wnote=note(dQ)
	String class=StringByKey("waveClass",wnote,"=")
	class = ReplaceString("Q_Positions",class,"Q_Positions_RGB")
	wnote = ReplaceStringByKey("waveClass",wnote,class,"=")
	wnote = AddClassToWaveNote(wnote,"RGB")
	wnote = ReplaceStringByKey("HueFromWave",wnote,NameOfWave(dQ),"=")
	if (useIntensity)
		wnote = ReplaceStringByKey("IntensFromWave",wnote,NameOfWave(IntensIN),"=")
		wnote = ReplaceNumberByKey("IntensityScalPower",wnote,power,"=")
	endif
	wnote = ReplaceNumberByKey("minColor",wnote,Qlo,"=")	// max and min values, used for the color scaling
	wnote = ReplaceNumberByKey("maxColor",wnote,Qhi,"=")
	Note/K RGB, wnote
	CopyScales/I dQ, RGB
	if (printIt && WaveExists(RGB))
		printf "Created new RGB wave  '%s'\r",NameOfWave(RGB)
	endif

	Wave Q_Coordinates = $StringFromList(0,WaveListClass("CoordinatesQdist","*","DIMS:2,MINCOLS:9,MAXCOLS:9"))
	if (WaveExists(Q_Coordinates))
		Variable Ncoords=DimSize(Q_Coordinates,0)
		name = NameOfWave(Q_Coordinates)
		Make/N=(Ncoords,3)/O/W/U $(name+"_RGB")/WAVE=Q_Coordinates_RGB =0	// holds the image to plot, color
		print "  also creating ",NameOfWave(Q_Coordinates_RGB)
		if (WaveDims(RGB)==4)
			Q_Coordinates_RGB = RGB[Q_Coordinates[p][6]][Q_Coordinates[p][7]][Q_Coordinates[p][8]][q]
		elseif (WaveDims(RGB)==3)
			Q_Coordinates_RGB = RGB[Q_Coordinates[p][6]][Q_Coordinates[p][7]][q]
		elseif (WaveDims(RGB)==2)
			Q_Coordinates_RGB = RGB[Q_Coordinates[p][6]][q]
		else
			return RGB
		endif
		Note/K Q_Coordinates_RGB, ReplaceStringByKey("waveClass",note(RGB),"Q_Coordinates_RGB,RGB","=")
	endif

	return RGB
End



Static Function/S FillQhist1image(image,sinTheta,indexWaveQ,Qhist,QhistNorm,mask,[maskNorm,dark,I0normalize])// fill Qhist from one file, F is for fast, sin(theta) is precomputed
	Wave image							// *******   NOTE, this wave will be altered!!!!  *******
	Wave sinTheta											// array of sin(theta)'s that is same size as the image to be loaded	
	Wave indexWaveQ										// index sorted wave of sinTheta
	Wave Qhist
	Wave QhistNorm											// array that contains number of pixels contributing to each bin in Qhist
	Wave mask												// optional mask, only use pixels that are true, ONLY 0 or 1 (do not use 255)
	Variable maskNorm										// use pixels outside of mask to normalize whole image (suppresses pedestal)
	Wave dark												// an optional background wave
	Variable I0normalize								// a Flag, if True then normalize data (default is True)
	maskNorm = ParamIsDefault(maskNorm) ? 0 : maskNorm
	maskNorm = numtype(maskNorm) ? 0 : !(!maskNorm)
	I0normalize = ParamIsDefault(I0normalize) ? 0 : I0normalize	// default to old way (un-normalized)
	I0normalize = numtype(I0normalize) ? 0 : !(!I0normalize)

	Variable timer = startMSTimer					// start timer for initial processing
	if (!WaveExists(image))
		print "input image does not exist"
		timer = stopMSTimer(timer)					// stop timer
		return ""
	endif

	if (WaveExists(dark))
		if (WaveType(image) & 0x58)					// either 8bit (0x08) or 16bit (0x10) or un-signed (0x40)
			Redimension/I image
		endif
		MatrixOP/FREE image = image - dark
	endif
	String wnote = note(image)
	Variable keV = NumberByKey("keV", wnote,"=")
	if (numtype(keV))
		DoAlert 0,"invalid keV in image wave note"
		timer = stopMSTimer(timer)					// stop timer
		return ""
	endif
	ImageStats/M=1/Q image								// check for empty images, skip them
	if (V_min==0 && V_max==0)
		timer = stopMSTimer(timer)					// stop timer
		wnote = ReplaceNumberByKey("V_min",wnote,V_min,"=")	// flag that this is empty
		wnote = ReplaceNumberByKey("V_max",wnote,V_max,"=")
		return wnote										// image is empty, do not waste time adding zeros
	endif

	Variable inorm = 	I0normalize ? I0normalizationFromNote(wnote) : 1	// normalization by intenstiy & exposure time

	NVAR secExtra=secExtra, secExtra1=secExtra1
	if (NVAR_Exists(secExtra))
		Variable timeExtra = startMSTimer			// start timer for intermediate processing
	endif
	Variable i, j, Ni=DimSize(image,0), Nj=DimSize(image,1)
	Variable k												// index into Qhist for each pixel
	Variable NQ=numpnts(Qhist)
	Variable Qmin = DimOffset(Qhist,0)
	Variable dQ = DimDelta(Qhist,0)
	Variable factor = 4*PI*keV/hc
	MatrixOp/FREE kQimage = (factor*sinTheta-Qmin)/dQ	// this method is 4.5 x faster
	MatrixOp/FREE kQimage = round(kQimage)		// this is now the index into Qhist for each point
	if (WaveExists(mask))
		Variable avgDarkPixel = 0
		if (maskNorm)
			ImageStats/M=1/R=mask image				// note, mask is 1 where we want, so it is correct for the ROI!
			avgDarkPixel = V_avg
		endif
		MatrixOp/FREE image = (image-avgDarkPixel)*mask	// mask values needs to be 1 or 0
	endif
	Variable N=numpnts(image)
	Redimension/N=(N) image,kQimage					// unwrap for faster processing
	SetScale/P x 0,1,"", image, kQimage

	if (NVAR_Exists(secExtra1))
		Variable timeExtra1 = startMSTimer			// start timer for intermediate processing
	endif
	IndexSort indexWaveQ, image						// sinTheta was already sorted, so kQimage is already sorted
	//	MatrixOp/FREE image = waveMap2(image,indexWaveQ,N)	// this is only slightly faster

	if (NVAR_Exists(secExtra1))
		secExtra1 += stopMSTimer(timeExtra1)/1e6
	endif
	// for each pixel in image, accumulate intensity into Qhist
	Variable i0=0,i1										// range of image to sum
	i0 = BinarySearch(kQimage,-0.1)
	i0 = (i0==-1) ? 0 : i0
	i0 = (i0==-2) ? NQ-1 : i0
	i0 += (kQimage[i0]<0)
	k = kQimage[i0]
	for (;k<NQ;)											// for each bin in Qhist
		i1 = BinarySearch(kQimage,k)					// this returns the last index with a value of k
		if (i1<i0)
			break
		endif
		Qhist[k] = sum(image,i0,i1)
		QhistNorm[k] += i1-i0+1						// number of pixels used for this Q bin
		i0 = I1 + 1											// reset start point of region
		k = i0<N ? kQimage[i0] : NQ					// do not index kQimage out of its range
	endfor
	if (NVAR_Exists(secExtra))
		secExtra += stopMSTimer(timeExtra)/1e6
	endif
	if (inorm!=1)
		Qhist *= inorm										// apply the normalization
	endif
	Variable seconds = stopMSTimer(timer)/1e6	// stop timer here
	if (ItemsInList(GetRTStackInfo(0))<3 && seconds>2.0)
		printf "	processing image took %s\r",Secs2Time(seconds,5,2)
	endif
	return wnote
End
//Static Function/S FillQhist1imageFile(fullName,sinTheta,indexWaveQ,Qhist,QhistNorm,mask,[maskNorm,dark,I0normalize])// fill Qhist from one file, F is for fast, sin(theta) is precomputed
//	String fullName
//	Wave sinTheta										// array of sin(theta)'s that is same size as the image to be loaded	
//	Wave indexWaveQ									// index sorted wave of sinTheta
//	Wave Qhist
//	Wave QhistNorm										// array that contains number of pixels contributing to each bin in Qhist
//	Wave mask											// optional mask, only use pixels that are true, ONLY 0 or 1 (do not use 255)
//	Variable maskNorm									// use pixels outside of mask to normalize whole image (suppresses pedestal)
//	Wave dark											// an optional background wave
//	Variable I0normalize							// a Flag, if True then normalize data (default is True)
//	maskNorm = ParamIsDefault(maskNorm) ? 0 : maskNorm
//	maskNorm = numtype(maskNorm) ? 0 : !(!maskNorm)
//	I0normalize = ParamIsDefault(I0normalize) ? 0 : I0normalize	// default to old way (un-normalized)
//	I0normalize = numtype(I0normalize) ? 0 : !(!I0normalize)
//
//	Variable timer = startMSTimer				// start timer for initial processing
//	Wave image = $(LoadGenericImageFile(fullName, extras="EscanOnly:1"))	// load image
//	if (!WaveExists(image))
//		printf "could not load image named '%s'\r",fullName
//		timer = stopMSTimer(timer)				// stop timer
//		return ""
//	endif
//	if (WaveExists(dark))
//		if (WaveType(image) & 0x58)				// either 8bit (0x08) or 16bit (0x10) or un-signed (0x40)
//			Redimension/I image
//		endif
//		MatrixOP/O image = image - dark
//	endif
//	String wnote = note(image)
//	Variable keV = NumberByKey("keV", wnote,"=")
//	if (numtype(keV))
//		DoAlert 0,"invalid keV in image '"+fullName+"'"
//		timer = stopMSTimer(timer)				// stop timer
//		return ""
//	endif
//	ImageStats/M=1/Q image							// check for empty images, skip them
//	if (V_min==0 && V_max==0)
//		timer = stopMSTimer(timer)				// stop timer
//		wnote = ReplaceNumberByKey("V_min",wnote,V_min,"=")	// flag that this is empty
//		wnote = ReplaceNumberByKey("V_max",wnote,V_max,"=")
//		KillWaves/Z image
//		return wnote									// image is empty, do not waste time adding zeros
//	endif
//
//	Variable inorm = 	I0normalize ? I0normalizationFromNote(wnote) : 1	// normalization by intenstiy & exposure time
//
//	NVAR secExtra=secExtra, secExtra1=secExtra1
//	if (NVAR_Exists(secExtra))
//		Variable timeExtra = startMSTimer		// start timer for intermediate processing
//	endif
//	Variable i, j, Ni=DimSize(image,0), Nj=DimSize(image,1)
//	Variable k											// index into Qhist for each pixel
//	Variable NQ=numpnts(Qhist)
//	Variable Qmin = DimOffset(Qhist,0)
//	Variable dQ = DimDelta(Qhist,0)
//	Variable factor = 4*PI*keV/hc
//	MatrixOp/O kQimage = (factor*sinTheta-Qmin)/dQ	// this method is 4.5 x faster
//	MatrixOp/O kQimage = round(kQimage)		// this is now the index into Qhist for each point
//	if (WaveExists(mask))
//		Variable avgDarkPixel = 0
//		if (maskNorm)
//			ImageStats/M=1/R=mask image			// note, mask is 1 where we want, so it is correct for the ROI!
//			avgDarkPixel = V_avg
//		endif
//		MatrixOp/O image = (image-avgDarkPixel)*mask	// mask values needs to be 1 or 0
//	endif
//	Variable N=numpnts(image)
//	Redimension/N=(N) image,kQimage				// unwrap for faster processing
//	SetScale/P x 0,1,"", image, kQimage
//
//	if (NVAR_Exists(secExtra1))
//		Variable timeExtra1 = startMSTimer		// start timer for intermediate processing
//	endif
//	IndexSort indexWaveQ, image					// sinTheta was already sorted, so kQimage is already sorted
//	//	MatrixOp/O image = waveMap2(image,indexWaveQ,N)	// this is only slightly faster
//
//	if (NVAR_Exists(secExtra1))
//		secExtra1 += stopMSTimer(timeExtra1)/1e6
//	endif
//	// for each pixel in image, accumulate intensity into Qhist
//	Variable i0=0,i1									// range of image to sum
//
//	i0 = BinarySearch(kQimage,-0.1)
//	i0 = (i0==-1) ? 0 : i0
//	i0 = (i0==-2) ? NQ-1 : i0
//	i0 += (kQimage[i0]<0)
//	k = kQimage[i0]
//
//	for (;k<NQ;)										// for each bin in Qhist
//		i1 = BinarySearch(kQimage,k)				// this returns the last index with a value of k
//		if (i1<i0)
//			break
//		endif
//		Qhist[k] = sum(image,i0,i1)
//		QhistNorm[k] += i1-i0+1					// number of pixels used for this Q bin
//		i0 = I1 + 1										// reset start point of region
//		k = kQimage[i0]
//	endfor
//	if (NVAR_Exists(secExtra))
//		secExtra += stopMSTimer(timeExtra)/1e6
//	endif
//	if (inorm!=1)
//		Qhist *= inorm									// apply the normalization
//	endif
//	Variable seconds = stopMSTimer(timer)/1e6	// stop timer here
//	if (ItemsInList(GetRTStackInfo(0))<3 && seconds>2.0)
//		printf "	processing image '%s' took %s\r",NameOfWave(image),Secs2Time(seconds,5,2)
//	endif
//	KIllWaves/Z image, kQimage
//	return wnote
//End

Static Function I0normalizationFromNote(wnote)	// normalization by intenstiy & exposure time
	String wnote

	Variable I0cnts,I0gain,I0gainBase,ScalerCountTime,exposure,I0scaling
	I0gainBase = NumVarOrDefault("root:Packages:micro:DEFAULT_I0_GAIN",DEFAULT_I0_GAIN)
	I0cnts = NumberByKey("I0",wnote,"=")
	I0gain = NumberByKey("I0gain",wnote,"=")
	I0gain = numtype(I0gain) ? I0gainBase : I0gain
	ScalerCountTime = NumberByKey("ScalerCountTime",wnote,"=")
	exposure = NumberByKey("exposure",wnote,"=")
	exposure = (numtype(exposure)==0 && exposure>0) ? exposure : 1
	I0scaling = NumVarOrDefault("root:Packages:micro:DEFAULT_I0_SCALING",DEFAULT_I0_SCALING)
	I0scaling = (numtype(I0scaling)==0 && I0scaling>0) ? I0scaling : 1
	Variable inorm = I0scaling/exposure
	if (numtype(I0cnts+I0gain+ScalerCountTime)==0 && I0cnts>0 && I0gain>0 && ScalerCountTime>0)
		inorm *= (ScalerCountTime/I0cnts) * (I0gainBase/I0gain)
	endif
	inorm = numtype(inorm)==0 && inorm>0 ? inorm : 1	// always return a valid inorm
	return inorm
End


Static Function/T requestPathFileRootFmt(pathName,suffixes)
	// returns a full path and fmt string as a list, e.g. "HD:Users:name:data:;EW1_%d_%d.h5"
	// returns the list "path;fileFormat"
	String pathName
	Variable suffixes				// maximum number of suffixes to remove, each suffix looks like "_12", leave the underscore, use NaN for auto-determine
										// suffixes are ALWAYS preceeded with an "_"
	PathInfo $pathName
	pathName = SelectString(V_flag,"",pathName)
	String message="pick any reconstructed image file in range,  using datafolder = "+GetDataFolder(0)
	Variable refNum
	initEnergyScans()
	SVAR ImageFileFilters=root:Packages:imageDisplay:ImageFileFilters
	Open/F=ImageFileFilters/D/M=message/P=$pathName/R refNum
	if (strlen(S_fileName)<1)
		return ""
	endif
	String fullPath = ParseFilePath(1,S_fileName,":",1,0)
	String name = ParseFilePath(3,S_fileName,":",0,0)
	String extension = ParseFilePath(4,S_fileName,":",0,0)
	extension = SelectString(strlen(extension),"",".")+extension		// extension now starts with a dot, or is empty
	String nameFmt = ""

	suffixes = suffixes>0 ? suffixes : Inf
	Variable i
	if (suffixes<1)						// nothing to do
		nameFmt = name
		return fullPath+";"+nameFmt+extension
	else										// desired number is given, or it is really big
		for (i=1;suffixes>=1 && i>0;suffixes-=1)	// strip off all but last suffix
			i = findTrailingIndex(name)
			if (i>=0)
				nameFmt = "_%d" + nameFmt
				name = name[0,i]
			endif
		endfor
		nameFmt = name + nameFmt
		return fullPath+";"+nameFmt+extension
	endif
End
//
Static Function findTrailingIndex(name)
	String name
	Variable N=strlen(name), i=-1
	if (GrepString(name,"_[0-9]+$"))		// true if name ends in "_123" where 123 is any simple unsigned integer
		i = strsearch(name,"_",N-1,1)-1
	endif
	return i
End


Static Function/T getNranges(pathName,nameFmt,[progressWin,printIt])
	String pathName		// probably 'imagePath'
	String nameFmt			// file name format string (not path info), something like  "EW5_%d.h5", or "EW1_%d_%d.h5"
	String progressWin	// if present, then update progress window.
	Variable printIt
	progressWin = SelectString(ParamIsDefault(progressWin),progressWin,"")
	printIt = ParamIsDefault(printIt) ? NaN : printIt
	printIt = numtype(printIt) ? strlen(GetRTStackInfo(2))<1 : !(!printIt)

	Variable Nranges = calcNranges(nameFmt)	// number of ranges needed
	if (Nranges>3)
		print "ERROR -- getNranges(), only know how to deal with up to 3 ranges"
		return ""
	elseif (Nranges<1)
		return ""
	endif
	progressWin = SelectString(ItemsInList(WinList(progressWin,";","WIN:64")),"",progressWin)

	String dlist=directory(pathName), name, subStr
	Variable p0,p1	, iN							// pointers into dlist, iN is number of Items in subStr
	Variable dlen = strlen(dlist)
	Variable im,i,m,N=ItemsInList(dlist), i0,i1,i2
	Make/N=(N,Nranges)/U/I/FREE ns
	for (p0=0,m=0,im=0; p0<dlen; p0=p1+1)	// need to break up dlist into chunks because of StringFromList()
		if (mod(im,50)==0 && strlen(progressWin))	// every 50 chunks
			if (ProgressPanelUpdate(progressWin,im/N*100))
				return ""							//   and break out of loop, and return nothing
			endif
		endif
		p1 = p0 + 10000							// work in chunks of 10000
		p1 = strsearch(dlist,";",p1)
		p1 = p1<0 ? dlen-1 : p1
		subStr = dlist[p0,p1]
		iN = ItemsInList(subStr)
		for (i=0; i<iN; i+=1,im+=1)			// fill this chunk
			name = StringFromList(i,subStr)	// StringFromList() is REALLY slow for big lists
			if (Nranges==3)
				sscanf name,nameFmt, i0,i1,i2
			elseif (Nranges==2)
				sscanf name,nameFmt, i0,i1
			else
				sscanf name,nameFmt, i0
			endif
			if (V_flag==Nranges)
				switch(Nranges)
					case 3:
						ns[m][2] = i2
					case 2:
						ns[m][1] = i1
					case 1:
						ns[m][0] = i0
				endswitch
				m += 1								// increment pointer into fileList[]
			endif
		endfor
	endfor
	if (strlen(progressWin))					// set 1to 100%
		ProgressPanelUpdate(progressWin,100)
	endif
	N = m
	Redimension/N=(N,-1) ns				// trim to correct length

	Make/N=(N)/U/I/FREE ni
	String rangeOut="", list
	Variable k,ii
	for (i=0;i<Nranges;i+=1)				// for each of the ranges
		list = ""
		ni = ns[p][i]
		Sort ni, ni								// sort each of the columns independently
		for (m=0,k=NaN; m<N; m+=1)		// make a list of all of the ith number in this column (with no repeats)
			ii = ni[m]
			if (ii!=k && ii<4.2e+9)
				list += num2istr(ii)+";"	// list of unique numbers in this column
				k = ii
			endif
		endfor
		rangeOut += compressRange(list,";")+";"
	endfor

	if (printIt)
		String range
		for (i=0;i<Nranges;i+=1)
			range = StringFromList(i,rangeOut)
			printf "range%d	 %d 		 %s\r",i+1,ItemsInRange(range),range
		endfor
	endif
	return rangeOut
End
//	Function test()
//		String nameFmt = "FS-slice3-EW008_%d_%d.h5"
//		String pathName = "imagePath"
//	
//		String progressWin=""
//		progressWin = ProgressPanelStart("",stop=1,showTime=1)
//		Variable msStart = stopMSTimer(-2)
//		print EnergyScans#getNranges(pathName,nameFmt, progressWin=progressWin, printIt=1)
//		print (stopMSTimer(-2)-msStart)*1e-6," seconds"
//		DoWindow/K $progressWin
//	End


Static Function/T fullNameFromFmt(fmt,i1,i2,i3)
	String fmt
	Variable i1,i2,i3

	String name=""
	Variable Nranges = calcNranges(fmt)	// number of ranges needed
	if (Nranges==1)
		sprintf name,fmt,i1
	elseif (Nranges==2)
		sprintf name,fmt,i1,i2
	elseif (Nranges==3)
		sprintf name,fmt,i1,i2,i3
	endif
	return name
End


Static Function calcNranges(nameFmt)
	String nameFmt
	Variable i, Nranges=0						// number of ranges needed
	for (i=1; i>0; i+=1)
		i = strsearch(nameFmt,"_%d",i,2)
		Nranges += (i >= 0) ? 1 : 0
	endfor
	return Nranges
End


Static Function/WAVE GetDistortionMap(roi)	// if using the distortion, precompute for all images here
	STRUCT imageROIstruct &roi
	if (imageROIstructBad(roi))
		return $""
	endif

	Variable startx = roi.xLo, starty = roi.yLo
	Variable groupx = roi.binx, groupy = roi.biny
	Variable Ni = roi.Nx,  Nj = roi.Ny
	Abort "GetDistortionMap(), filling the distortion map needs serious fixing"
	// check if distiortion map already exists
	Variable distortionOK=0
	if(exists("root:Packages:geometry:tempCachedDistortionMap")==1)	// wave exists, so check its size and location
		Wave DistortionMap = root:Packages:geometry:tempCachedDistortionMap
		distortionOK = (DimOffset(DistortionMap,0)==startx) && 1
		distortionOK = (DimOffset(DistortionMap,1)==starty) && 1
		distortionOK = (DimDelta(DistortionMap,0)==groupx) && 1
		distortionOK = (DimDelta(DistortionMap,1)==groupy) && 1
		distortionOK = (DimSize(DistortionMap,0)==Ni) && 1
		distortionOK = (DimSize(DistortionMap,1)==Nj) && 1
	endif
	if (!distortionOK)									// existing map (if it exists) is not OK, so make a new one
		Variable microSec = stopMSTimer(-2)		// for timing creation of local distortion map
		KIllWaves/Z root:Packages:geometry:tempCachedDistortionMap	// ensure that this is gone!
		Make/N=(Ni,Nj,2)/O root:Packages:geometry:tempCachedDistortionMapTemp
		Wave distortionMap = root:Packages:geometry:tempCachedDistortionMapTemp
		SetScale/P x startx,groupx,"", distortionMap	// scale distortionMap[][] to full chip pixels
		SetScale/P y starty,groupy,"", distortionMap
		Variable/C dxy
		Wave xymap = root:Packages:geometry:xymap
		Variable i,j, px,py
		for (j=0;j<Nj;j+=1)
			for (i=0;i<Ni;i+=1)
				py = starty + j*groupy
				px = startx + i*groupx
//					dxy = peakcorrection2(xymap,px,py)		// returns cmplx(dx,dy)
				distortionMap[i][j][0] = real(dxy)		// accumulate all the distortions that I will need
				distortionMap[i][j][1] = imag(dxy)
			endfor
		endfor
		Rename root:Packages:geometry:tempCachedDistortionMapTemp, tempCachedDistortionMap
		Variable seconds = (stopMSTimer(-2)-microSec)/1e6
		if (ItemsInList(GetRTStackInfo(0))<3 && seconds>0.2)
			printf "creating local copy of distortion map took %s\r",Secs2Time(seconds,5,3)
		endif
	endif
	Note/K DistortionMap, "use=1"												// distortion map now ready
	return DistortionMap
End


Function/WAVE QhistFromQpositions(source,i,j,k)		// from a Q_posiitions type array, extract the Qhist (a 1d array)
	Wave source
	Variable i,j,k												// index into source to id where to get the Qhist
	Variable dimLast = WaveDims(source)-1				// dimension of Q dimension (always the last one)
	if (dimLast<0)
		return $""
	endif
	String Qname = GetWavesDataFolder(source,1)+"Qhist"
	Make/N=(DimSize(source,dimLast))/O/D $Qname		// re-make array to hold Q's from one image
	Wave Qhist = $Qname
	SetScale/P x DimOffset(source,dimLast),DimDelta(source,dimLast),WaveUnits(source,dimLast), Qhist
	SetScale d 0,0,WaveUnits(source,-1), Qhist

	k = dimLast<3 || numtype(k)>0 || k<0 ? NaN : k
	j = dimLast<2 || numtype(j)>0 || j<0 ? NaN : j
	i = dimLast<1 || numtype(i)>0 || i<0 ? NaN : i
	if (numtype(k)==0)
		Qhist = source[i][j][k][p]
	elseif (numtype(j)==0)
		Qhist = source[i][j][p]
	elseif (numtype(i)==0)
		Qhist = source[i][p]
	else
		Qhist = source[p]
	endif

	// re-set wave note for Qhist
	String wnote=note(source), labelQ=""
	wnote = ReplaceStringByKey("waveClass",wnote,"Qhistogram","=")
	wnote = ReplaceStringByKey("sourceWave",wnote,GetWavesDataFolder(source,2),"=")
	wnote = RemoveByKey("sourceWave_row",wnote,"=")		// clean out old source_row, _col, _layer, _chunk
	wnote = RemoveByKey("sourceWave_col",wnote,"=")
	wnote = RemoveByKey("sourceWave_layer",wnote,"=")
	wnote = RemoveByKey("sourceWave_chunk",wnote,"=")
	String labelX = StringByKey("labelX",wnote,"=")		// save labelX, Y, Z, W
	String labelY = StringByKey("labelY",wnote,"=")
	String labelZ = StringByKey("labelZ",wnote,"=")
	String labelW = StringByKey("labelW",wnote,"=")
	wnote = RemoveByKey("labelX",wnote,"=")					// clean out old labelX, Y, Z, W
	wnote = RemoveByKey("labelY",wnote,"=")
	wnote = RemoveByKey("labelZ",wnote,"=")
	wnote = RemoveByKey("labelW",wnote,"=")
	switch(dimLast)
		case 3:
			wnote = ReplaceNumberByKey("sourceWave_layer",wnote,k,"=")
			if (strlen(labelW))										// move old labelW --> source_labelW
				labelQ = SelectString(strlen(labelQ),labelW,labelQ)	// label of Q
				wnote = ReplaceStringByKey("source_labelW",wnote,labelW,"=")
			endif
		case 2:
			wnote = ReplaceNumberByKey("sourceWave_col",wnote,j,"=")
			if (strlen(labelZ))										// move old labelZ --> source_labelZ
				labelQ = SelectString(strlen(labelQ),labelZ,labelQ)
				wnote = ReplaceStringByKey("source_labelZ",wnote,labelZ,"=")
			endif
		case 1:
			wnote = ReplaceNumberByKey("sourceWave_row",wnote,i,"=")
			if (strlen(labelY))
				labelQ = SelectString(strlen(labelQ),labelY,labelQ)
				wnote = ReplaceStringByKey("source_labelY",wnote,labelY,"=")
			endif
		case 1:
			if (strlen(labelX))
				labelQ = SelectString(strlen(labelQ),labelX,labelQ)
				wnote = ReplaceStringByKey("source_labelX",wnote,labelX,"=")
			endif
			break
		default:
			return $""
	endswitch
	labelQ = SelectString(strlen(labelQ),"Q",labelQ)		// default value is "Q"
	wnote = ReplaceStringByKey("labelX",wnote,labelQ,"=")
	String labelValue = StringByKey("labelValue",wnote,"=")
	labelValue = SelectString(strlen(labelValue),"Intensity",labelValue)// default value is "Intensity"
	wnote = ReplaceStringByKey("labelValue",wnote,labelValue,"=")

	Note/K Qhist, wnote
	// fitOneQhist(Qhist)
	return Qhist
End


Function/T fitOneQhist(Qhist, [quiet,printIt,lo,hi,d0])
	Wave Qhist						// a 1-d Q distribution
	Variable quiet					// only fit, NO fit_Qhist wave and NO ouput, NO graph updating (goes much faster too)
	Variable printIt				// controls print out of peak information to history
	Variable lo,hi					// range of Qhist to use for fitting (you can pass -Inf,Inf for whole range
	Variable d0						// unstrained d-spacing (nm)
	quiet = ParamIsDefault(quiet) ? NaN : quiet
	quiet = numtype(0) ? 0 : !(!quiet)			// default is NOT quiet
	printIt = ParamIsDefault(printIt) ? NaN : printIt
	printIt = numtype(0) ? 0 : !(!printIt)	// default is NO print out
	lo = ParamIsDefault(lo) ? NaN : lo
	hi = ParamIsDefault(hi) ? NaN : hi
	if (!WaveExists(Qhist))
		return ""
	endif
	d0 = ParamIsDefault(d0) ? NaN : d0
	d0 = numtype(d0) || d0<=0 ? NumVarOrDefault("d0",NaN) : d0
	WaveStats/M=1/Q Qhist
	if (V_npnts<6)
		return ""
	endif
	String fldrSav= GetDataFolder(1)
	SetDataFolder  $GetWavesDataFolder(Qhist,1)	// data folder with the data

	lo = max(lo,0)										// lo and hi must fit on wave
	hi = min(DimSize(Qhist,0)-1,hi)
	if (numtype(lo+hi)==2)							// need to find fitting range, NaNs are not acceptable
		if(getPercentOfPeak(Qhist,0.5,lo,hi,5))	// find the 50% points
			lo = NaN										// getPercentOfPeak failed, set lo,hi to bad values
			hi = NaN
		endif
	endif

//		Make/T/O T_Constraints={"K0 > 0","K3 > 0"}	// ensure width>0, and base line positive
//		WaveStats/Q Qhist
//		T_Constraints[0] = "K0 > "+num2str(min(0,V_min))
//		CurveFit/Q lor Qhist/D /C=T_Constraints 
//		KillWaves/Z T_Constraints
	Variable V_FitOptions=4, V_FitError=0		//suppress window & ignore errors, I will trap them
	Variable fwhm=NaN, Qc=NaN, QpkArea=NaN
	Variable fwhmErr=NaN, QcErr=NaN
	Variable bkg=NaN, bkgErr=NaN
	Variable VoigtShape=NaN, VoigtShapeErr=NaN
	Variable sl2 = sqrt(ln(2))
	WaveStats/M=1/Q/R=[lo,hi]/Z Qhist			// check for at least 4 valid points in range
	if (numtype(lo+hi)==0 && V_npnts>4)		// for a valid range, fit the peak
		Variable fitLineShape=NumVarOrDefault("root:Packages:micro:Escan:fitLineShape",0)	// 0=Lorentz, 1=Gauss, 2=Voigt
		if (fitLineShape==0)							// Lorentzian
			CurveFit/Q/N=(quiet) lor Qhist[lo,hi]/AD=(!quiet)
		elseif (fitLineShape==1)						// Gaussian
			CurveFit/Q/N=(quiet) gauss Qhist[lo,hi]/AD=(!quiet)
		elseif (fitLineShape==2)						// Voigt
			CurveFit/Q/N gauss Qhist[lo,hi]			// this is just done to get the starting point
			if (!V_FitError)
				Wave W_coef=W_coef
				Redimension/N=5 W_coef					// reset to fit one Voigt function
				W_coef[2] = 1/K3							// width parameter
				W_coef[3] = K2								// center of peak
				W_coef[4] = 0.02							// start the fit near a Gaussian (easy and works OK)
				FuncFit/Q/N=(quiet) VoigtFit W_coef Qhist[lo,hi]/AD=(!quiet)
			endif
		endif
		if (!V_FitError)
			if (fitLineShape==0)
				fwhm = 2*sqrt(K3)							// for Lorentzian
				Qc = K2
				QpkArea = K1*2*PI/fwhm
				Wave W_sigma=W_sigma
				fwhmErr = 2*abs(W_sigma[3])/2
				QcErr = abs(W_sigma[2])
				bkg = K0 ; 	bkgErr = abs(W_sigma[0])
			elseif (fitLineShape==1)
				fwhm = 2*K3*sl2							// for Gaussian
				Qc = K2
				QpkArea = K1 * abs(K3)*sqrt(PI)
				Wave W_sigma=W_sigma
				fwhmErr = 2*W_sigma[3]*sl2
				QcErr = abs(W_sigma[2])
				bkg = K0 ; 	bkgErr = abs(W_sigma[0])
			elseif (fitLineShape==2)
				Variable c2=W_coef[2], c4=W_coef[4]
				Variable hwg=sl2/c2
				Variable hwl=c4/c2 
				fwhm = 2* (hwl/2 + sqrt( hwl^2/4 + hwg^2) )	// for Voigt
				Qc = W_coef[3]
				QpkArea = W_coef[1]*sqrt(PI)/c2
				VoigtShape = c4
				Wave W_sigma=W_sigma
				Variable denomRoot = sqrt( c4^2/(4*c2^2) + ln(2)/c2^2 )
				Variable e4 = 1/abs(2*c2) + abs(c4)/ ( 4*c2^2 * denomRoot)		// d(hw)/d(c4)
				Variable e2 = -1/c2^3 * ( c4^2/2 + 2*ln(2) )
				e2 /= 2*denomRoot
				e2 -= c4/(2*c2^2)
				fwhmErr = 2 * (abs(e4*W_sigma[4]) + abs(e2*W_sigma[2]) )
				QcErr = abs(W_sigma[3])
				bkg = K0 ; 	bkgErr = abs(W_sigma[0])
				VoigtShapeErr = W_sigma[4]
			endif
		endif
	endif

	String list = ReplaceNumberByKey("Qcenter","",Qc,"=")
	list = ReplaceNumberByKey("Qfwhm",list,fwhm,"=")
	list = ReplaceNumberByKey("QpkArea",list,QpkArea,"=")
	list = ReplaceNumberByKey("QcErr",list,QcErr,"=")
	list = ReplaceNumberByKey("fwhmErr",list,fwhmErr,"=")
	String wnote = note(Qhist)
	wnote = ReplaceNumberByKey("Qcenter",wnote,Qc,"=")
	wnote = ReplaceNumberByKey("Qfwhm",wnote,fwhm,"=")
	if (numtype(QpkArea)==0)
		wnote = ReplaceNumberByKey("QpkArea",wnote,QpkArea,"=")
	endif
	wnote = ReplaceNumberByKey("QcErr",wnote,QcErr,"=")
	wnote = ReplaceNumberByKey("fwhmErr",wnote,fwhmErr,"=")
	list = ReplaceNumberByKey("QpkBkg",list,bkg,"=")
	list = ReplaceNumberByKey("QpkBkgErr",list,bkgErr,"=")
	list = ReplaceNumberByKey("chisq",list,V_chisq,"=")

	String caller=StringFromList(0,GetRTStackInfo(0))
//	Variable printIt = !quiet && (stringmatch(caller,"EscanButtonProc") || stringmatch(caller,"refitQhistogramButtonProc"))
//	Variable printIt = !quiet && (stringmatch(caller,"EscanButtonProc"))
	if (printIt)
		printf "at peak, Q = %s(1/nm),   fwhm = %s(1/nm),   area = %.3g",ValErrStr(Qc,QcErr),ValErrStr(fwhm,fwhmErr),QpkArea
		if (fitLineShape==0)
			printf ",  Lorentzian shape"
		elseif (fitLineShape==1)
			printf ",  Gaussian shape"
		elseif (fitLineShape==2)
			printf ",  Voigt shape = %s",ValErrStr(VoigtShape,VoigtShapeErr)
		endif
		printf ",      chisq=%.3g\r",V_chisq
	endif
	if (numtype(d0)==0 && d0>0)
		Variable Q0 = 2*PI/d0							// Q of unstrained material (1/nm)
		Variable strain = -(Qc-Q0)/Qc					// d is inverse of strain, so a negative
		Variable strainErr = QcErr/Q0
		list = ReplaceNumberByKey("d0",list,d0,"=")
		list = ReplaceNumberByKey("Q0",list,Q0,"=")
		list = ReplaceNumberByKey("strain",list,strain,"=")
		list = ReplaceNumberByKey("strainErr",list,strainErr,"=")
		if (printIt)
			printf "   d0 = %g (nm),   Q0 = %g (1/nm),   strain = %s\r",d0,Q0,ValErrStr(strain,strainErr)
		endif
		list = ReplaceNumberByKey("ÆQ",list,Qc-Q0,"=")
		wnote = ReplaceNumberByKey("ÆQ",wnote,Qc-Q0,"=")
		wnote = ReplaceNumberByKey("strain",wnote,strain,"=")
		wnote = ReplaceNumberByKey("strainErr",wnote,strainErr,"=")
	endif
	Note/K Qhist, wnote

	// if Qhist is on a graph, update textbox and line, so find a plot (if one exists, containing Qhist)
	String win = StringFromList(0,WindowsWithWave(Qhist,1))	// find top most graph containing Qhist
	if (!quiet && strlen(win))
		if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EscanButtonProc"))
			DoWindow/F $win
		endif
		if (WhichListItem("textQ", AnnotationList(win),";")>=0)		// textQ exists, update it
			TextBox/W=$win/C/N=textQ textQ_fromWnote(wnote)
		endif
		if (Q0>0)
			SetDrawLayer/W=$win/K UserFront		// draw dotted line at Q0
			SetDrawEnv/W=$win xcoord= bottom,ycoord= prel,dash= 1
			DrawLine/W=$win Q0,0,Q0,1
		endif
	endif
	SetDataFolder fldrSav
	return list
End


Static Function Escan_SelectLineShapeFunction(shape)
	String shape					// can only be "gauss" or "lorentzian", otherwise user is prompted

	Variable fitLineShapeOld=NumVarOrDefault("root:Packages:micro:Escan:fitLineShape",0)
	Variable fitLineShape=fitLineShapeOld
	String fList="Lorentzian;Gaussian;Voigt"

	if (StringMatch(shape,"gauss*"))
		fitLineShape = 1
	elseif (StringMatch(shape,"lor*"))
		fitLineShape = 0
	elseif (StringMatch(shape,"voigt*"))
		fitLineShape = 0
	else
		fitLineShape += 1
		Prompt fitLineShape,"Line Shape Fitting Function",popup,fList
		DoPrompt "Line Shape",fitLineShape
		if (V_flag)
			return NaN
		endif
		fitLineShape -= 1
	endif

	NewDataFolder/O root:Packages
	NewDataFolder/O root:Packages:micro
	NewDataFolder/O root:Packages:micro:Escan
	Variable/G root:Packages:micro:Escan:fitLineShape = fitLineShape

	printf "\r\r\r"
	if (fitLineShape==fitLineShapeOld)
		printf "Line Shape Function UNCHANGED, still %s\r",StringFromList(fitLineShape,fList)
	else
		printf "Line Shape Function CHANGED from '%s' --> '%s'\r",StringFromList(fitLineShapeOld,fList),StringFromList(fitLineShape,fList)
	endif
	printf "\r\r\r"
	return fitLineShape 
End


Static Function QfromCursorHook(s)		// updates the Q distribution plot by the cursor on a Q_positions image
	STRUCT WMWinHookStruct &s
	if (s.eventCode!=7 || !stringmatch(s.cursorName,"A"))	// only trap cursor A moved events
		return 0									// a zero lets other hooks process this event
	endif

	Variable i,j,k=NaN
	String win=s.winName
	Wave waveDisplayed = ImageNameToWaveRef(win,StringFromList(0,ImageNameList(win,";")))
	if (WaveExists(waveDisplayed))
		i = s.pointNumber
		j = s.yPointNumber
	else
		Wave waveDisplayed = TraceNameToWaveRef(win, StringFromList(0,TraceNameList(win,";",1)))
		if (!WaveExists(waveDisplayed))
			return 0
		endif
		i = waveDisplayed[s.pointNumber][6]
		j = waveDisplayed[s.pointNumber][7]
	endif
	String wnote = note(waveDisplayed)
	Wave source = $StringByKey("sourceWave", wnote,"=")
	if (!WaveExists(source))
		return 0
	endif
	Wave Qhist = QhistFromQpositions(source,i,j,NaN)
	fitOneQhist(Qhist,quiet=0,printIt=0)
	return 1										// a one stops other hooks from processing too
End
//
//	SetWindow kwTopWin hook(QfromCursor)=QfromCursorHook	// this forces a re-calculation of Qhist when cursor A is moved
//	SetWindow kwTopWin userdata(arrows)="value=root:test:EsumDepth;proc=SetEsumDepthProc;stepSize=1;hiLim=71;loLim=-1;"

// ========================== End of Q distributions any Position or Depth, NEW Section ==========================
// ===============================================================================================================



// ===============================================================================================================
// ========================================= Start of Single Qhist Plot  =========================================

Function/S MakeGraph_Qhist(Qhist,[printIt])	// display a Qhist[] wave
	Wave Qhist
	Variable printIt
	printIt = ParamIsDefault(printIt) ? NaN : printIt
	printIt = numtype(printIt) ? 0 : !(!printIt)

	if (!WaveExists(Qhist))
		String QhistName = StringFromLIst(0,WaveListClass("Qhistogram","*","DIMS:1"))
		if (strlen(QhistName))
			Wave Qhist = $QhistName
		endif
		printIt = 1
	endif
	if (!WaveExists(Qhist))
		DoAlert 0, "no 'Qhist' wave exists in this data folder"
		return ""
	endif
	if (printIt)
		printf "MakeGraph_Qhist(%s",NameOfWave(Qhist)
		if (!ParamIsDefault(printIt))
			printf ", printIt=%g",printIt
		endif
		printf ")\r"
	endif
	String graphName=CleanupName("Graph_"+GetDataFolder(0)+"_Qhist",0)
	if (strlen(WinList(graphName,"","WIN:1")))
		DoWindow/F $graphName									// graph is already up, so just bring it to front
		return graphName
	elseif (exists(graphName)==5)	
		Execute graphName+"()"									// a recreation macro exists, run it
		return graphName
	endif
	String wnote = note(Qhist)

	String fitName = GetWavesDataFolder(Qhist,1)+"fit_"+NameOfWave(Qhist)
	Display/K=1/W=(732,69,1262,429) Qhist
	DoWindow/C $graphName
	if (exists(fitName)==1)
		AppendToGraph $fitName
		ModifyGraph lSize($NameOfWave($fitName))=2
	endif
	ModifyGraph tick=2, mirror=1, minor=1, lowTrip=0.001, standoff=0, zero(left)=3
	ModifyGraph mode[0]=4, marker[0]=19, lStyle[0]=1, rgb[0]=(16385,16388,65535), msize[0]=2, hbFill[0]=4
	ModifyGraph axOffset(left)=-3.85714,axOffset(bottom)=-0.461538
	String labelValue = StringByKey("labelValue",wnote,"=")
	labelValue = SelectString(strlen(labelValue),"Intensity",labelValue)	// default label for Y-axis is "Intensity"
	Label left labelValue+FixUnitsInLabel(WaveUnits(Qhist,-1))
	String labelX = StringByKey("labelX",wnote,"=")
	labelX = SelectString(strlen(labelX),"Q",labelX)	// default label for X-axis is "Q"
	Label bottom labelX+FixUnitsInLabel(WaveUnits(Qhist,0))
	SetAxis/N=1 left 0,*
	ShowInfo
	String str = StringByKey("title",wnote,"=")		// first try to get title from wave note
	if (strlen(str)<1)
		str = StrVarOrDefault("title","")					// nothing in wavenote, try the title global string
	endif
	str = SelectString(strlen(str),"",str+"\r\\Zr075")
	str += "\\JR"
	String singleImage=StringByKey("singleImage",wnote,"=")
	if (strlen(singleImage))
		str += singleImage
	else
		str += GetWavesDataFolder(Qhist,1)
	endif
	if (strlen(StringByKey("dateExposed",note(Qhist),"=")))
		str += "\r"+StringByKey("dateExposed",note(Qhist),"=")
	endif
	String iso = StringByKey("file_time",note(Qhist),"=")
	if (strlen(iso))
		if (StringMatch(iso[10],"T"))
			iso = ISOtime2niceStr(iso)
		endif
		str += "\r"+iso
	endif
	TextBox/A=RT/C/N=textTitle/F=0 str

	Variable Qc = NumberByKey("Qcenter",wnote,"="), fwhm = NumberByKey("Qfwhm",wnote,"=")
	Variable dQ = NumberByKey("ÆQ",wnote,"="), strain = NumberByKey("strain",wnote,"=")
	Variable Q0 = NumberByKey("Q0",wnote,"=")
	TextBox/C/N=textQ/F=0/A=LT/B=1 SelectString(numtype(Qc+fwhm),textQ_fromWnote(wnote),"")
	if (Q0>0)
		SetDrawLayer UserFront
		SetDrawEnv xcoord= bottom,ycoord= prel,dash= 1
		DrawLine Q0,0,Q0,1
	endif

	Variable single = strlen(singleImage) || single ? 1 : NaN
	if (numtype(single))
		Wave sourceWave=$StringByKey("sourceWave",wnote,"=")
		if (WaveExists(sourceWave))
			Variable i,Nmax=0
			for (i=0;i<WaveDims(sourceWave);i+=1)
				Nmax = max(Nmax,DimSize(sourceWave,i))
			endfor
			single = numpnts(sourceWave) == Nmax
		else
			single = 1
		endif
	endif
	if (single)
		Button refitButtonQhist,pos={41,1},size={50,20},proc=EnergyScans#refitQhistogramButtonProc,title="re fit"
	else
		Button SaveQhistTraceButton,pos={1,20},size={40,30},proc=SaveQhistTraceButtonProc,title="Keep\rTrace",fSize=10
		Button SaveQhistTraceButton,help={"duplicate the generic trace with a new name and leave it on plot"}
	endif
	return graphName
End
//
Static Function refitQhistogramButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	if (ba.eventCode!=2 || !stringmatch(ba.ctrlName,"refitButtonQhist"))	// only execute on button up
		return 0								// a zero lets other hooks process this event
	endif

	String tList = TraceNameList("",";", 1 )
	Variable i
	for (i=0;i<ItemsInList(tList);i+=1)
		Wave Qhist=TraceNameToWaveRef(ba.win,StringFromList(i,tList))
		if (stringmatch(StringByKey("waveClass",note(Qhist),"="),"Qhistogram"))
			break
		endif
	endfor
	if (WaveExists(Qhist))
		fitOneQhist(Qhist,quiet=0,printIt=1)
	endif
	return 0
End
//
Function SaveQhistTraceButtonProc(ctrlName) : ButtonControl
	String ctrlName

	if (stringmatch(ctrlName,"SaveQhistTraceButton"))
		Wave w = TraceNameToWaveRef("","Qhist")
	endif
	if (!WaveExists(w))
		return 1
	endif

	String wnote = note(w)
	Variable row,col,layer,chunk=NaN
	row = NumberByKey("sourceWave_row",wnote,"=")
	col = NumberByKey("sourceWave_col",wnote,"=")
	layer = NumberByKey("sourceWave_layer",wnote,"=")
	Wave source = $StringByKey("sourceWave",wnote,"=")
	String trace
	if (WaveExists(source) && numtype(row)==0)			// using new style
		trace = ArrayTraceOnGraph("",source,row,col,layer,chunk)	// finds if an array trace (e.g. A[4][*][2]) is on graph, and returns trace name
		if (strlen(trace)==0)								// not on graph, put it there
			String sName = GetWavesDataFolder(source,2)
			if (numtype(layer)==0)
				AppendToGraph $sName[row][col][layer][*]
			elseif (numtype(col)==0)
				AppendToGraph $sName[row][col][*]
			elseif (numtype(row)==0)
				AppendToGraph $sName[row][*]
			endif
			trace = ArrayTraceOnGraph("",source,row,col,layer,chunk)	// get trace name again so I can modify it color, etc.
			if (strlen(trace))								// this should always be true, we just put it there
				ModifyGraph mode($trace)=4,marker($trace)=16,msize($trace)=2,lstyle($trace)=1,rgb($trace)=(9,37540,0)
			endif
		endif

	else													// the OLD way, should be able to get rid of this
		row = NumberByKey("Qhist_row",wnote,"=")
		String qName = GetWavesDataFolder(w,1)+"Q_Depth"
		if (exists(qName)!=1)
			return 1
		endif

		trace = QDepthCutOnGraph("",$qName,row)
		if (strlen(trace)<1)				// it is is not on the graph, put it there
			AppendToGraph $qName[row][*]
		endif
		trace =QDepthCutOnGraph("",$qName,row)
		if (strlen(trace))				// it is is not on the graph, put it there
			ModifyGraph mode($trace)=4,marker($trace)=16,msize($trace)=2,lstyle($trace)=1,rgb($trace)=(9,37540,0)
		endif
	endif

	return 0
End
//
Static Function/S QDepthCutOnGraph(graphName,qw,row)
	String graphName	// name of the graph
	Wave qw			// the Q_Depth 2-d wave
	Variable row		// row in Q_Depth to check for

	String traceList = TraceNameList(graphName,";",1)

	Variable i, N=ItemsInList(traceList)
	String trace
	String info=" "

	String path = GetWavesDataFolder(qw,2)
	Variable instance,rr
	for (i=0;i<N;i+=1)
		trace = StringFromList(i,traceList)
		Wave wt = TraceNameToWaveRef(graphName,trace)
		if (!stringmatch(GetWavesDataFolder(wt,2),path))
			continue
		endif

		sscanf trace, "Q_Depth#%d",instance
		instance = (V_flag<1) ? 0 : instance
		info = StringByKey("YRANGE",TraceInfo(graphName,"Q_Depth",instance))
		sscanf info, "[%d][*]",rr

		if (V_flag==1 && rr==row)
			return trace
		endif
	endfor
	return ""
End
//
Static Function/S ArrayTraceOnGraph(graphName,source,row,col,layer,chunk)	// finds if an array trace (e.g. A[4][*][2]) is on graph, and returns trace name
	String graphName	// name of the graph
	Wave source		// the multi-dim wave
	Variable row,col,layer,chunk	// used to id where trace comes from, at least one must be NaN

	Variable gotVariable=0
	String info0=""
	if (numtype(row)==0)
		info0 += "["+num2istr(row)+"]"
	else
		info0 += "[*]"
		gotVariable = 1
	endif
	if (numtype(col)==0)
		info0 += "["+num2istr(col)+"]"
	elseif (!gotVariable)
		info0 += "[*]"
		gotVariable = 1
	endif
	if (numtype(layer)==0)
		info0 += "["+num2istr(layer)+"]"
	elseif (!gotVariable)
		info0 += "[*]"
		gotVariable = 1
	endif
	if (numtype(chunk)==0)
		info0 += "["+num2istr(chunk)+"]"
	elseif (!gotVariable)
		info0 += "[*]"
		gotVariable = 1
	endif
	if (!gotVariable)
		info0 += "[*]"
		gotVariable = 1
	endif

	String traceList = TraceNameList(graphName,";",1)
	Variable i, N=ItemsInList(traceList)
	String trace, info,fmt
	String name = NameOfWave(source)
	String path = GetWavesDataFolder(source,2)
	Variable instance
	for (i=0;i<N;i+=1)
		trace = StringFromList(i,traceList)
		Wave wt = TraceNameToWaveRef(graphName,trace)
		if (!stringmatch(GetWavesDataFolder(wt,2),path))
			continue
		endif
		sscanf trace, name+"#%d",instance
		instance = (V_flag<1) ? 0 : instance
		info = StringByKey("YRANGE",TraceInfo(graphName,name,instance))
		if (cmpstr(info,info0)==0)
			return trace
		endif
	endfor
	return ""
End
//
Static Function/S textQ_fromWnote(wnote)
	String wnote

	String str, textQ=""
	Variable Qc = NumberByKey("Qcenter",wnote,"=")
	Variable fwhm = NumberByKey("Qfwhm",wnote,"=")
	Variable strain = NumberByKey("strain",wnote,"=")
	Variable QcErr = NumberByKey("QcErr",wnote,"=")
	Variable fwhmErr = NumberByKey("fwhmErr",wnote,"=")
	Variable strainErr = NumberByKey("strainErr",wnote,"=")
	if (numtype(QcErr)==0)
		sprintf str, "Q\\Bc\\M = %g%s%.2g (nm\\S-1\\M)",Qc,PLUSMINUS,QcErr
	else
		sprintf str, "Q\\Bc\\M = %g (nm\\S-1\\M)",Qc
	endif
	textQ = str
	if (numtype(fwhmErr)==0)
		sprintf str, "\rFWHM = %.2g%s%.1g",fwhm,PLUSMINUS,fwhmErr
	else
		sprintf str, "\rFWHM = %.2g",fwhm
	endif
	textQ += str
	if (numtype(strain)==0 && numtype(strainErr)==0)
//		sprintf str,"\rstrain = %.1e %s%.1e",strain,PLUSMINUS,strainErr
		sprintf str,"\r\\F'Symbol'e\\F]0 = %.1e %s%.1e",strain,PLUSMINUS,strainErr
		textQ += str
	elseif (numtype(strain)==0)
//		sprintf str,"\rstrain = %.1e",strain
		sprintf str,"\r\\F'Symbol'e\\F]0 = %.1e",strain
		textQ += str
	endif
	Wave wav = $StringByKey("sourceWave",wnote,"=")
	if (WaveExists(wav))
		Variable m = NumberByKey("sourceWave_row",wnote,"=")
		if (numtype(m)==0)
			textQ += "\r"+StringByKey("source_labelX",wnote,"=")+" = "+num2str(DimOffset(wav,0) + m*DimDelta(wav,0))+" "+WaveUnits(wav,0)
			m = NumberByKey("sourceWave_col",wnote,"=")
			if (numtype(m)==0)
				textQ += "\r"+StringByKey("source_labelY",wnote,"=")+" = "+num2str(DimOffset(wav,1) + m*DimDelta(wav,1))+" "+WaveUnits(wav,1)
				m = NumberByKey("sourceWave_layer",wnote,"=")
				if (numtype(m)==0)
					textQ += "\r"+StringByKey("source_labelZ",wnote,"=")+" = "+num2str(DimOffset(wav,2) + m*DimDelta(wav,2))+" "+WaveUnits(wav,2)
					m = NumberByKey("sourceWave_chunk",wnote,"=")
					if (numtype(m)==0)
						textQ += "\r"+StringByKey("source_labelW",wnote,"=")+" = "+num2str(DimOffset(wav,3) + m*DimDelta(wav,3))+" "+WaveUnits(wav,3)
					endif
				endif
			endif
		endif
	endif
	return textQ
End

// ========================================== End of Single Qhist Plot  ==========================================
// ===============================================================================================================



// ===============================================================================================================
// ================================= Start of Q distribution from 1 loaded image =================================

//	Process a single loaded image.  This pays no attention to any positioners, just process one loaded image
//	a depth of 0 is assumed if image wave note does not say otherwise.
Function Fill_Q_1image(d0,image,[depth,mask,maskNorm,dark,I0normalize,printIt,ask])	// does not assume depth in image
	Variable d0			// d-spacing of the strain=0 material (nm)
	Wave image			// image to process
	Variable depth
	Wave mask			// optional mask to limit the pixels that get processed (use pixel when mask true)
	Variable maskNorm	// use pixels outside of mask to normalize whole image (suppresses pedestal)
	Wave dark			// an optional background wave
	Variable I0normalize // a Flag, if True then normalize data (default is True)
	Variable printIt
	Variable ask		// force a prompt
	depth = ParamIsDefault(depth) ? 0 : depth
	depth = numtype(depth) ? 0 : depth
	maskNorm = ParamIsDefault(maskNorm) ? 0 : maskNorm
	maskNorm = numtype(maskNorm) ? 0 : !(!maskNorm)
	I0normalize = ParamIsDefault(I0normalize) ? 1 : I0normalize
	I0normalize = numtype(I0normalize) ? 1 : !(!I0normalize)
	printIt = ParamIsDefault(printIt) ? NaN : !(!printIt)
	printIt = numtype(printIt) ? strlen(GetRTStackInfo(2))<1 : printIt
	ask = ParamIsDefault(ask) ? NaN : ask
	ask = numtype(ask) ? 0 : !(!ask)

	if (!((d0>0)) && exists("d0")!=2)				// there is no d0, check with the user
		DoAlert 2, "There is no d0, so you cannot get strain, just Q.  Do you want to set d0?"
		if (V_flag==1)										// chose "Yes"
			Variable dnm = NumVarOrDefault("d0",NaN)
			Prompt dnm, "referencing d-spacing (nm)"
			DoPrompt "d0",dnm
			if (V_flag)
				return 1
			endif
			Variable/G $(GetDataFolder(1)+"d0")=dnm
			d0 = dnm
		elseif (V_flag==3)								// chose "Cancel"
			return 1											// quitting
		endif
	endif
	if (!((d0>0)))											// invalid d0, check for a local value
		d0 = NumVarOrDefault("d0",NaN)
		printIt = 1
	endif

	String str, maskName = "", imageName=""
	if (!WaveExists(image))							// no image passed, ask for one
		String iList = WaveListClass("rawImage*","*","DIMS:2")
		String maskList = "-none-;"+WaveListClass("imageMask","*","DIMS:2,BYTE:1")	// setup to use optional mask
		maskList = SelectString(WaveExists(mask),maskList,"")
		if (ItemsInList(iList)==0)
			return 1
		elseif (ItemsInList(iList)>1 || ask)
			Prompt imageName,"image to process",popup,iList
			Prompt depth,"depth to assume (micron)"
			Prompt maskName, "mask to use with image",popup,maskList
			Prompt maskNorm, "subtract dark using pixels outside of mask",popup,"Nothing;Suppress"
			Prompt I0normalize,"Normalize the Q-hists",popup,"Do NOT Normalize to Io or exposure;Normalize to Io & exposure"
			I0normalize += 1
			maskNorm += 1
			DoPrompt "image",imageName,depth,I0normalize,maskName,maskNorm
			I0normalize -= 1
			maskNorm = maskNorm==2
			if (V_flag)
				return 1
			endif
			Wave image = $imageName
			maskName = SelectString(stringmatch(maskName,"-none-"),maskName,"")
			if (exists(maskName)==1)
				Wave mask = $maskName					// do not check if wave exists, that a valid option
			endif
			printIt = 1
		else
			imageName = StringFromList(0,iList)
			Wave image = $imageName
			printIt = 1
		endif
	endif
	if (!WaveExists(image) || numtype(depth))	// no image passed or bad depth, failure
		return 1
	endif
	imageName = NameOfWave(image)
	I0normalize = numtype(I0normalize) ? 1 : I0normalize

	String darkList = WaveListClass("imageDark;rawImageDark","*","")
	if (!WaveExists(dark) && ItemsInList(darkList)>0)
		String darkName=""
		Prompt darkName,"Background Image",popup,"-none-;"+darkList
		DoPrompt "Background",darkName
		if (V_flag)
			return 1
		endif
		if (cmpstr(darkName,"-none-"))
			Wave dark = $darkName
		else
			Wave dark = $""
		endif
	endif

	if (printIt)
		if (WaveExists(mask))
			maskName = GetWavesDataFolder(mask,2)
		endif
		maskName = SelectString(WaveExists(mask),"$\"\"",maskName)
		printf "Fill_Q_1image(%g, %s",d0,imageName
		if (!ParamIsDefault(depth))
			printf ", depth=%g",depth
		endif
		if (WaveExists(mask))
			printf ", mask=%s",maskName
		endif
		if (maskNorm>0)
			printf ", maskNorm=%g",maskNorm
		endif
		if (WaveExists(dark))
			printf ", dark=%s",NameOfWave(dark)
		endif
		if (!ParamIsDefault(I0normalize))
			printf ", I0normalize=%g",I0normalize
		endif
		if (!ParamIsDefault(printIt))
			printf ", printIt=%g",printIt
		endif
		printf ")\r"
	endif
	if (WaveExists(mask))
		if (sum(mask)==0)
			DoAlert 0, "You picked a mask that is all zero, stopping"
			return 1
		endif
	endif

	STRUCT microGeometry geo
	FillGeometryStructDefault(geo)
	Variable useDistortion = NumVarOrDefault("root:Packages:geometry:useDistortion",USE_DISTORTION_DEFAULT)
 	if (printIt)
		printf "processing with distotion "+SelectString(useDistortion,"OFF","on")
		print "   "+SelectString(WaveExists(mask),"and no mask","and using  '"+NameOfWave(mask)+"'  for the mask")
		if (maskNorm)
			print "   Renormalizing pixels using pixels outside of the mask"
		endif
		print SelectString(I0normalize,"   NOT normalizing","   Normalizing")+" to the ion chamber and exposure time"
	endif

	// open the first file in the range, to get info about the images
	Variable Ni,Nj, Npixels=numpnts(image)		// number of pixels in one image
	Ni = DimSize(image,0)
	Nj = DimSize(image,1)
	String wnote = note(image)
	Variable keV, startx, endx, starty, endy, groupx, groupy
	keV = NumberByKey("keV", wnote,"=")					// energy for this image
	Variable X1=NumberByKey("X1", wnote,"="), H1=NumberByKey("H1", wnote,"="), F1=NumberByKey("F1", wnote,"=")
	startx = NumberByKey("startx", wnote,"=");	endx = NumberByKey("endx", wnote,"=");	groupx = NumberByKey("groupx", wnote,"=")
	starty = NumberByKey("starty", wnote,"=");	endy = NumberByKey("endy", wnote,"=");	groupy = NumberByKey("groupy", wnote,"=")
	if (numtype(startx+endx+starty+endy+groupx+groupy))
		DoAlert 0,"could not get ROI from wave note of image '"+imageName+"'"
		return 1
	elseif (numtype(keV))
		DoAlert 0,"invalid keV in image '"+imageName+"'"
		return 1
	endif
	if (numtype(d0)==0 && d0>0)
		wnote = ReplaceNumberByKey("d0",wnote,d0,"=")
	endif

	// note that Q range only depends upon image size and energy
	Variable Qmin, Qmax, dQ, NQ										// get range of Q
	Variable dNum = detectorNumFromID(geo, StringByKey("detectorID", note(image),"="))
	if (!(dNum>=0 && dNum<MAX_Ndetectors))
		DoAlert 0,"could not get detector number from from wave note of image '"+imageName+"' using detector ID"
		return 1
	endif

	STRUCT imageROIstruct roi
	imageROIstructInit(roi, wnote=note(image))					// returns 0=OK, non-zero if problem
	Variable/C thetaZ = thetaRange(geo.d[dNum],roi,maskIN=mask,depth=depth)		// theta range on this image
	Qmin = 4*PI*sin(real(thetaZ))*keV/hc							// min Q (1/nm)
	Qmax = 4*PI*sin(imag(thetaZ))*keV/hc							// max Q (1/nm)

	// determine dQ (1/nm), the Q resolution to use.  Base it on the distance between two adjacent pixels
	Variable px,py															// the current pixel to analyze (unbinned full chip pixels)
	py = round((starty+endy)/2)										// approximate full chip pixel position in center or the image
	px = round((startx+endx)/2)
	dQ = 4*PI*abs(sin(pixel2q(geo.d[dNum],px,py,$""))-sin(pixel2q(geo.d[dNum],px+groupx,py+groupy,$"")))*keV/hc
	//		dQ /= 1.5
	NQ = round((Qmax-Qmin)/dQ) + 1									// number of Qs

	if (printIt)
		printf "Q range = [%g, %g] (1/nm),  %sQ=%.2g(1/nm)       ",Qmin, Qmax,GDELTA,dQ
		printf "theta range = [%g, %g]%s\r",real(thetaZ)*180/PI,imag(thetaZ)*180/PI,DEGREESIGN
		if (d0>0)
			printf "   d0 = %g (nm),   Q0 = %g (1/nm)\r",d0,2*PI/d0
		endif
	endif

	// done with the setup part, now actually compute something
	Variable i, j
	if (useDistortion)													// if using the distortion, precompute for all images  here
		Abort "Fill_Q_Positions(), filling the distortion map needs serious fixing"
		// check if distiortion map already exists
		Variable distortionOK=0
		if(exists("root:Packages:geometry:tempCachedDistortionMap")==1)	// wave exists, so check its size and location
			Wave DistortionMap = root:Packages:geometry:tempCachedDistortionMap
			distortionOK = (DimOffset(DistortionMap,0)==startx) && 1
			distortionOK = (DimOffset(DistortionMap,1)==starty) && 1
			distortionOK = (DimDelta(DistortionMap,0)==groupx) && 1
			distortionOK = (DimDelta(DistortionMap,1)==groupy) && 1
			distortionOK = (DimSize(DistortionMap,0)==Ni) && 1
			distortionOK = (DimSize(DistortionMap,1)==Nj) && 1
		endif
		if (!distortionOK)												// existing map (if it exists) is not OK, so make a new one
			KIllWaves/Z root:Packages:geometry:tempCachedDistortionMap	// ensure that this is gone!
			Make/N=(Ni,Nj,2)/O root:Packages:geometry:tempCachedDistortionMapTemp
			Wave distortionMap = root:Packages:geometry:tempCachedDistortionMapTemp
			SetScale/P x startx,groupx,"", distortionMap		// scale distortionMap[][] to full chip pixels
			SetScale/P y starty,groupy,"", distortionMap
			Variable/C dxy
			Wave xymap = root:Packages:geometry:xymap
			for (j=0;j<Nj;j+=1)
				for (i=0;i<Ni;i+=1)
					py = starty + j*groupy
					px = startx + i*groupx
//					dxy = peakcorrection2(xymap,px,py)				// returns cmplx(dx,dy)
					distortionMap[i][j][0] = real(dxy)				// accumulate all the distortions that I will need
					distortionMap[i][j][1] = imag(dxy)
				endfor
			endfor
			Rename root:Packages:geometry:tempCachedDistortionMapTemp, tempCachedDistortionMap
				if (printIt)
				printf "creating local copy of distortion map\r"
			endif
		endif
		Note/K DistortionMap, "use=1"
	endif																		// distortion map now ready

	Variable timer0=startMSTimer, seconds							// start timer for entire prooess
	Variable timer=startMSTimer										// start timer for initial processing
	Make/N=(Ni,Nj)/D/FREE sinTheta
	Note/K sinTheta, wnote
	CopyScales/I image, sinTheta
	sinTheta = NaN
	// for each pixel in image, compute sin(theta) and save it in sinTheta[][]
	Make/N=(Ni)/I/FREE Xpixels
	Make/N=(Nj)/I/FREE Ypixels
	Ypixels = p*groupy + (groupy-1)/2 + starty
	Xpixels = p*groupx + (groupx-1)/2 + startx
	sinTheta = pixel2q(geo.d[dNum],Xpixels[p],Ypixels[q],$"",depth=depth)	// first set to theta
	sinTheta = sin(sinTheta)											// convert theta -> sin(theta)
	WAVEClear Xpixels,Ypixels
	seconds = stopMSTimer(timer)/1e6								// stop timer here
	if (printIt)
		printf "filling array of sin(theta) from  '%s'  at depth = %g,  took %s\r",NameOfWave(image),depth,Secs2Time(seconds,5,1)
	endif

	Make/N=(Ni*Nj)/I/FREE indexWaveQ
	indexWaveQ = p
	Sort sinTheta, sinTheta,indexWaveQ								// sort so indexWaveQ[0] is index to lowest sin(theta), indexWaveQ[inf] is greatest

	print "starting bulk of processing"
	timer = startMSTimer												// start a timer on bulk of processing
	// accumulate the Q histogram for one image into Qhist
	Make/N=(NQ)/O/D Qhist=0										// array to hold Q's from one image
	Make/N=(NQ)/FREE/D QhistNorm=0								// number of pixels used for each Qhist, use to normalize
	SetScale/I x Qmin,Qmax,"1/nm", Qhist

	Duplicate/FREE image, imageTemp								// needed because FillQhist1image() changes the input image
	wnote = FillQhist1image(imageTemp,sinTheta,indexWaveQ,Qhist,QhistNorm,mask,dark=dark,maskNorm=maskNorm,i0normalize=i0normalize)	// fill Qhist from one file
	wnote = ReplaceStringByKey("waveClass",wnote,"Qhistogram","=")
	wnote = ReplaceStringByKey("singleImage",wnote,GetWavesDataFolder(image,2),"=")
	Note/K Qhist, wnote

	if (useDistortion)
		Note/K DistortionMap, "use=0"
	endif
	Qhist /= QhistNorm									// do the normalization
	Qhist = numtype(Qhist) ? NaN : Qhist
	seconds = stopMSTimer(timer)/1e6
	if (printIt)
		printf "  processing image took %s	(%.3g %ss/pixel)\r",Secs2Time(seconds,5,1),1e6*seconds/Npixels,Gmu
	endif

	String QfitInfo = fitOneQhist(Qhist,d0=d0,quiet=1,printit=1)
	wnote = MergeKeywordLists(wnote,QfitInfo,0,"=",";")

	String title = StrVarOrDefault("title","")
	if (strlen(title))
		title = ReplaceString("=",title,"_")							// in the wave note the string cannot have "=" or ";"
		title = ReplaceString(";",title,"_")
		wnote = ReplaceStringByKey("title",wnote,title,"=")
	endif
	if (WaveExists(mask))
		wnote = ReplaceStringByKey("maskWave",wnote,GetWavesDataFolder(mask,2),"=")
	endif
	str = SelectString(i0normalize,"Intensity","Intensity / I\\B0\\M")
	wnote = ReplaceStringByKey("labelValue",wnote,str,"=")
	Note/K Qhist, wnote

	// put up or refresh the plot of Qhist
	if (printIt)
		printf "at peak, Q = %.3f(1/nm),   fwhm = %.2g(1/nm)",NumberByKey("Qcenter", QfitInfo,"="),NumberByKey("Qfwhm", QfitInfo,"=")
		Variable strain = NumberByKey("strain", QfitInfo,"=")
		if (numtype(strain)==0)
			printf ",   %sQ = %.3f,   strain = %.1e",GDELTA,NumberByKey("ÆQ", QfitInfo,"="),strain
		endif
		printf "\r"
		String win = MakeGraph_Qhist(Qhist)
		i = MoveWinToTopRight(win,-1,-1)
	endif

	// done with processing, clean up
	seconds = stopMSTimer(timer0)/1e6
	if (printIt && seconds>0.2)
		printf "entire process took %s\r",Secs2Time(seconds,5,3)
	endif
	beep
	return 0
End
//
//	Static Function/S FillQhist1image(imageIn,sinTheta,indexWaveQ,Qhist,QhistNorm,mask,[maskNorm,dark,I0normalize])// fill Qhist from one file, F is for fast, sin(theta) is precomputed
//		Wave imageIn
//		Wave sinTheta											// array of sin(theta)'s that is same size as the image to be loaded	
//		Wave indexWaveQ										// index sorted wave of sinTheta
//		Wave Qhist
//		Wave QhistNorm											// array that contains number of pixels contributing to each bin in Qhist
//		Wave mask												// optional mask, only use pixels that are true, ONLY 0 or 1 (do not use 255)
//		Variable maskNorm										// use pixels outside of mask to normalize whole image (suppresses pedestal)
//		Wave dark												// an optional background wave
//		Variable I0normalize								// a Flag, if True then normalize data (default is True)
//		maskNorm = ParamIsDefault(maskNorm) ? 0 : maskNorm
//		maskNorm = numtype(maskNorm) ? 0 : !(!maskNorm)
//		I0normalize = ParamIsDefault(I0normalize) ? 0 : I0normalize	// default to old way (un-normalized)
//		I0normalize = numtype(I0normalize) ? 0 : !(!I0normalize)
//		if (!WaveExists(imageIn))
//			print "input image does not exist"
//			return ""
//		endif
//	
//		Duplicate/FREE imageIn, image
//		if (WaveExists(dark))
//			if (WaveType(image) & 0x58)					// either 8bit (0x08) or 16bit (0x10) or un-signed (0x40)
//				Redimension/I image
//			endif
//			MatrixOP/O image = image - dark
//		endif
//		String wnote = note(image)
//		Variable keV = NumberByKey("keV", wnote,"=")
//		keV = NumberByKey("keV", wnote,"=")
//		if (numtype(keV))
//			DoAlert 0,"invalid keV in image '"+NameOfWave(image)+"'"
//			return ""
//		endif
//		ImageStats/M=1/Q image							// check for empty images, skip them
//		if (V_min==0 && V_max==0)
//			wnote = ReplaceNumberByKey("V_min",wnote,V_min,"=")	// flag that this is empty
//			wnote = ReplaceNumberByKey("V_max",wnote,V_max,"=")
//			return wnote									// image is empty, do not waste time adding zeros
//		endif
//	
//		Variable inorm = 	I0normalize ? I0normalizationFromNote(wnote) : 1	// normalization by intenstiy & exposure time
//	
//		Variable timer = startMSTimer							// start timer for initial processing
//		Variable i, j, Ni=DimSize(image,0), Nj=DimSize(image,1)
//		Variable k														// index into Qhist for each pixel
//		Variable NQ=numpnts(Qhist)
//		Variable Qmin = DimOffset(Qhist,0)
//		Variable dQ = DimDelta(Qhist,0)
//		Variable factor = 4*PI*keV/hc
//		MatrixOp/FREE kQimage = (factor*sinTheta-Qmin)/dQ	// this method is 4.5 x faster
//		MatrixOp/FREE/O kQimage = round(kQimage)			// this is now the index into Qhist for each point
//		if (WaveExists(mask))
//			Variable avgDarkPixel = 0
//			if (maskNorm)
//				ImageStats/M=1/R=mask image						// note, mask is 1 where we want, so it is correct for the ROI!
//				avgDarkPixel = V_avg
//			endif
//			MatrixOp/O/FREE image = (image-avgDarkPixel)*mask	// mask values needs to be 1 or 0
//		endif
//		Variable N=numpnts(image)
//		Redimension/N=(N) image,kQimage							// unwrap for faster processing
//		SetScale/P x 0,1,"", image, kQimage
//	
//		IndexSort indexWaveQ, image								// sinTheta was already sorted, so kQimage is already sorted
//		//	MatrixOp/O/FREE image = waveMap2(image,indexWaveQ,N)	// this is slightly faster
//	
//		// for each pixel in image, accumulate intensity into Qhist
//		Variable i0=0,i1												// range of image to sum
//		i0 = BinarySearch(kQimage,-0.1)
//		i0 = (i0==-1) ? 0 : i0
//		i0 = (i0==-2) ? NQ-1 : i0
//		i0 += (kQimage[i0]<0)
//		k = kQimage[i0]
//		for (;k<NQ;)													// for each bin in Qhist
//			i1 = BinarySearch(kQimage,k)							// this returns the last index with a value of k
//			if (i1<i0)
//				break
//			endif
//			Qhist[k] = sum(image,i0,i1)
//			QhistNorm[k] += i1-i0+1								// number of pixels used for this Q bin
//			i0 = I1 + 1													// reset start point of region
//			k = i0<N ? kQimage[i0] : NQ							// do not index kQimage out of its range
//		endfor
//		if (inorm!=1)
//			Qhist *= inorm												// apply the normalization
//		endif
//		Variable seconds = stopMSTimer(timer)/1e6			// stop timer here
//		if (ItemsInList(GetRTStackInfo(0))<3 && seconds>2.0)
//			printf "	processing image took %s\r",Secs2Time(seconds,5,2)
//		endif
//		return wnote
//	End

// ================================== End of Q distribution from 1 loaded image ==================================
// ===============================================================================================================



// ===============================================================================================================
// ====================================== Start of Q-histograming Utility  =======================================

Static Function/WAVE MakeSinThetaArray(roi,d,wnote,[depth])	// make a free array the same size as roi, but filled with sin(theta)'s
	STRUCT imageROIstruct &roi
	STRUCT detectorGeometry &d					// the input detector geometry
	String wnote
	Variable depth										// usually determined from image file
	depth = ParamIsDefault(depth) ? NaN : depth

	Variable timer = startMSTimer						// start timer for initial processing
	if (imageROIstructBad(roi))
		print "ERROR -- MakeSinThetaArray(), invalid ROI"
		DoAlert 0, "ERROR -- MakeSinThetaArray(), invalid ROI"
		return $""
	endif
	if (numtype(depth))
		depth = NumberByKey("depth", wnote,"=")		// depth for this image
	endif
	depth = numtype(depth) ? 0 : depth					// default depth to zero

	Variable startx = roi.xLo, groupx = roi.binx
	Variable starty = roi.yLo, groupy = roi.biny
	Variable Ni = roi.Nx, Nj = roi.Ny					// number of BINNED pixels
	Make/N=(Ni,Nj)/D/FREE sinThetaCached = NaN
	SetScale/P x, 0, 1,"", sinThetaCached
	SetScale/P y, 0, 1,"", sinThetaCached

	wnote = ReplaceNumberByKey("startx",wnote,startx,"=")
	wnote = ReplaceNumberByKey("endx",wnote,roi.xHi,"=")
	wnote = ReplaceNumberByKey("groupx",wnote,groupx,"=")
	wnote = ReplaceNumberByKey("starty",wnote,starty,"=")
	wnote = ReplaceNumberByKey("endy",wnote,roi.yHi,"=")
	wnote = ReplaceNumberByKey("groupy",wnote,groupy,"=")
	Note/K sinThetaCached, wnote

	// for each pixel in image, compute sin(theta) and save it in sinThetaCached[][]
	Make/N=(Ni)/O/I/FREE Xpixels
	Make/N=(Nj)/O/I/FREE Ypixels
	Ypixels = p*groupy + (groupy-1)/2 + starty
	Xpixels = p*groupx + (groupx-1)/2 + startx
	sinThetaCached = pixel2q(d,Xpixels[p],Ypixels[q],$"",depth=depth)	// first set to theta
	sinThetaCached = sin(sinThetaCached)				// convert theta -> sin(theta)
	WAVEClear Xpixels,Ypixels

	Variable seconds = stopMSTimer(timer)/1e6		// stop timer here
	if (ItemsInList(GetRTStackInfo(0))<3 && seconds>5)
		printf "filling array of sin(theta) from  at depth = %g,  took %s\r",depth,Secs2Time(seconds,5,1)
	endif
	return sinThetaCached
End
//Static Function/S MakeSinThetaArray(fullName,geo,[depth])	// make an array the same size as an image, but filled with sin(theta)'s
//	String fullName										// complete path name to image
//	STRUCT microGeometry &geo						// the input geo
//	Variable depth										// usually determined from image file
//	depth = ParamIsDefault(depth) ? NaN : depth
//
//	Variable timer = startMSTimer						// start timer for initial processing
//	Wave image = $(LoadGenericImageFile(fullName))	// load image
//	if (!WaveExists(image))
//		printf "could not load image named '%s'\r",fullName
//		DoAlert 0,"could not load image"
//		return ""
//	endif
//	String wnote = note(image)
//	Variable startx, starty, groupx, groupy, dNum
//	startx = NumberByKey("startx", wnote,"=")
//	groupx = NumberByKey("groupx", wnote,"=")
//	starty = NumberByKey("starty", wnote,"=")
//	groupy = NumberByKey("groupy", wnote,"=")
//	if (numtype(depth))
//		depth = NumberByKey("depth", wnote,"=")		// depth for this image
//	endif
//	depth = numtype(depth) ? 0 : depth					// default depth to zero
//	dNum = detectorNumFromID(geo, StringByKey("detectorID", wnote,"="))
//	if (!(dNum>=0 && dNum<=2))
//		DoAlert 0,"could not get detector number from from wave note of image '"+NameOfWave(image)+"' using detector ID"
//		return ""
//	elseif (numtype(startx+starty+groupx+groupy))
//		DoAlert 0,"could not get ROI from wave note of image '"+fullName+"'"
//		return ""
//	endif
//
//	Variable Ni=DimSize(image,0), Nj=DimSize(image,1)
//	Make/N=(Ni,Nj)/O/D sinThetaCached
//	Note/K sinThetaCached, wnote
//	CopyScales/I image, sinThetaCached
//	sinThetaCached = NaN
//
//	// for each pixel in image, compute sin(theta) and save it in sinThetaCached[][]
//	Make/N=(Ni)/O/I sinThetaCachedPX
//	Make/N=(Nj)/O/I sinThetaCachedPY
//	sinThetaCachedPY = p*groupy + (groupy-1)/2 + starty
//	sinThetaCachedPX = p*groupx + (groupx-1)/2 + startx
//	sinThetaCached = pixel2q(geo.d[dNum],sinThetaCachedPX[p],sinThetaCachedPY[q],$"",depth=depth)	// first set to theta
//	sinThetaCached = sin(sinThetaCached)				// convert theta -> sin(theta)
//
//	KillWaves/Z sinThetaCachedPX,sinThetaCachedPY
//
//	Variable seconds = stopMSTimer(timer)/1e6		// stop timer here
//	if (ItemsInList(GetRTStackInfo(0))<3 && seconds>5)
//		printf "filling array of sin(theta) from  '%s'  at depth = %g,  took %s\r",NameOfWave(image),depth,Secs2Time(seconds,5,1)
//	endif
//	KIllWaves/Z image
//	return GetWavesDataFolder(sinThetaCached,2)
//End


// returns true if ROI in list is the same as the other ROI
Static Function sameROI(list,startx, endx, starty, endy, groupx, groupy)
	String list											// list with current numbers
	Variable startx, endx, starty, endy, groupx, groupy	// desired numbers

	Variable ierr = abs(startx-NumberByKey("startx", list,"="))
	ierr += abs(endx-NumberByKey("endx", list,"="))
	ierr += abs(groupx-NumberByKey("groupx", list,"="))
	ierr += abs(starty-NumberByKey("starty", list,"="))
	ierr += abs(endy-NumberByKey("endy", list,"="))
	ierr += abs(groupy-NumberByKey("groupy", list,"="))
	return (ierr < 1e-9)								// true if the ROI's are equal
End
//
Static Function/C thetaRange(d,roi,[maskIN,depth])		// returns the range of theta spanned by image
	STRUCT detectorGeometry &d
	STRUCT imageROIstruct &roi									// an ROI of interest
	Wave maskIN															// mask must match roi size exactly
	Variable depth														// optional passed depth for when depth is not zero & not in wave note
	depth = ParamIsDefault(depth) ? NaN : depth			// if depth not passed, set to invalid value

	if (imageROIstructBad(roi))
		DoAlert 0, "ERROR -- thetaRange() ROI is invalid"
		return cmplx(NaN,NaN)
	endif
	Variable groupx = roi.binx, groupy = roi.biny
	Variable Nx = roi.Nx, Ny = roi.Ny

	if (WaveExists(maskIN) && !ParamIsDefault(maskIN))
		Wave mask = maskIN					// check all un masked pixels on the detector
		if (Nx!=DimSize(mask,0) || Ny!=DimSize(mask,1))
			DoAlert 0,"mask size does not match roi size"
			return cmplx(NaN,NaN)
		endif
	else											// if no mask passed, make one that is all ones
		Make/N=(Nx,Ny)/U/B/FREE maskOnes=1
		Wave mask = maskOnes
	endif

	Variable theta,thetaMin=Inf, thetaMax=-Inf	// Bragg angle (rad) and range of theta, test all points in the roi
	Variable ix,iy, px,py					// loop over all pixels
	for (iy=0;iy<Ny;iy+=1)
		py = roi.yLo + groupy*iy
		for (ix=0;ix<Nx;ix+=1)
			if (!mask[ix][iy])				// skip points masked off
				continue
			endif
			px = roi.xLo + groupx*ix
			theta = pixel2q(d,px,py,$"",depth=depth)
			thetaMin = min(theta,thetaMin)
			thetaMax = max(theta,thetaMax)
		endfor
	endfor
	return cmplx(thetaMin,thetaMax)
End
//	Static Function/C thetaRange(d,image,[maskIN,depth])		// returns the range of theta spanned by image
//		STRUCT detectorGeometry &d
//		Wave image
//		Wave maskIN
//		Variable depth											// optional passed depth for when depth is not zero & not in wave note
//		depth = ParamIsDefault(depth) ? NaN : depth			// if depth not passed, set to invalid value
//	
//		String wnote = note(image)
//		Variable startx, groupx, starty, groupy
//		startx = NumberByKey("startx", wnote,"=")
//		starty = NumberByKey("starty", wnote,"=")
//		groupx = NumberByKey("groupx", wnote,"=")
//		groupy = NumberByKey("groupy", wnote,"=")
//		depth = numtype(depth) ? NumberByKey("depth", wnote,"=") : depth	// if depth not passed, try wavenote
//		Variable Nx=DimSize(image,0), Ny=DimSize(image,1)
//		if (numtype(startx+groupx+starty+groupy))
//			DoAlert 0, "could not get ROI from wave note of image '"+NameOfWave(image)+"'"
//			return cmplx(NaN,NaN)
//		endif
//		if (WaveExists(maskIN) && !ParamIsDefault(maskIN))
//			Wave mask = maskIN					// check all un masked pixels on the detector
//			if (Nx!=DimSize(mask,0) || Ny!=DimSize(mask,1))
//				DoAlert 0,"mask size does not match image size"
//				return cmplx(NaN,NaN)
//			endif
//		else											// if no mask passed, make one that is all ones
//			Make/N=(Nx,Ny)/U/B/FREE maskOnes=1
//			Wave mask = maskOnes
//		endif
//	
//		Variable theta,thetaMin=Inf, thetaMax=-Inf	// Bragg angle (rad) and range of theta, test all 4 corners of the image
//		Variable ix,iy, px,py						// loop over all pixels
//		for (iy=0;iy<Ny;iy+=1)
//			py = starty + groupy*iy
//			for (ix=0;ix<Nx;ix+=1)
//				if (!mask[ix][iy])					// skip points masked off
//					continue
//				endif
//				px = startx + groupx*ix
//				theta = pixel2q(d,px,py,$"",depth=depth)
//				thetaMin = min(theta,thetaMin)
//				thetaMax = max(theta,thetaMax)
//			endfor
//		endfor
//		return cmplx(thetaMin,thetaMax)
//	End


#ifdef OLD_WAY_EnergyScans
// find the step size in a coordinate by examining the values
OverRide Function FindScalingFromVec(vec,threshold,first,stepSize,dimN)
	Wave vec
	Variable threshold					// a change greater than this is intentional (less is jitter)
	Variable &first						// use 	SetScale/P x first,stepSize,"",waveName
	Variable &stepSize					// returned |step size|
	Variable &dimN							// returned number of points

	threshold = abs(threshold)						// no negative thresholds
	first = NaN												// init to bad values
	stepSize = NaN
	dimN = 0

	// get the step size rounded to show no significant digits past threshold/10
	stepSize = FindStepSizeInVec(vec,threshold,signed=0)	// returns the |step size| from values in vec
	if (numtype(stepSize) || !WaveExists(vec))	// failed
		return 1
	elseif (stepSize==0)								// only one point, no actual scan
		first = vec[0]
		dimN = 1
		return 0
	endif

	// find the starting point, the median of all the points within threshold of the lowest
	// note that first is always the smallest, regardless of the direction of scan, stepSize is always positive
	Duplicate/FREE vec vecSort
	SetScale/P x 0,1,"", vecSort
	Sort vecSort, vecSort
	Variable i = floor(BinarySearchInterp(vecSort, vecSort[0]+threshold))
	first = vecSort(i/2)							// this gives median of the points within threshold of the lowest

	// find number of points in "one scan"
	Variable N = numpnts(vec)
	Variable tm=0,tn=0, m=0						// tm=# of scan found
	for (i=1;i<N;i+=1)
		if (abs(vec[i]-vec[i-1])<threshold)	// this step is a repeat, do not count it
			continue
		elseif (abs(abs(vec[i]-vec[i-1])-stepSize)<(abs(stepSize)+threshold))	// this an OK step
			m += 1
		else
			tm += m+1									// sum of number of points in "a scan"
			tn += 1										// number of scans found
			m = 0
		endif
	endfor
	if (m>0)
		tm += m+1										// sum of number of points in "a scan"
		tn += 1											// number of scans found
	endif
	dimN = round(tm/ tn)
	return 0
End
#endif


Static Function getPercentOfPeak(pkWave,fraction,lo,hi,minWid)	// return index in to pkWave of the ±50% of peak points, or use cursors A & B
	Wave pkWave					// wave with peak
	Variable fraction			// fraction of peak to go on each side
	Variable &lo, &hi
	Variable minWid			// minimum width, min value of (hi-lo)
	lo = NaN  ;  hi = NaN
	if (!WaveExists(pkWave) || !(0<fraction && fraction<1))
		return 1
	endif
	Variable N=numpnts(pkWave)-1

	String win = StringFromList(0,WindowsWithWave(pkWave,1))	// find top most graph containing pkWave
	if (strlen(win))
		lo = NumberByKey("POINT",CsrInfo(A,win))
		hi = NumberByKey("POINT",CsrInfo(B,win))
	else
		lo = NaN
		hi = NaN
	endif
	if (numtype(lo+hi)==0)
		if (lo>hi)
			Variable swap=lo
			lo = hi
			hi = swap
		endif
		return 0
	endif

	WaveStats/M=1/Q pkWave

	Variable imax,i
	imax = (V_maxloc-DimOffset(pkWave,0))/DimDelta(pkWave,0)// index of the max
	imax = round(imax)
	if (!(1<imax && imax<N-1))										// too close to the edges
		return 1
	endif

	Variable half = (V_max+max(V_min,0))*fraction					// half way point
	if (!(half<V_max))
		return 1
	endif

	for (lo=imax;lo>0;lo-=1)
		if (pkWave[lo]<half)
			break
		endif
	endfor
	lo = max(0,lo)
	for (hi=imax;hi<N;hi+=1)
		if (pkWave[hi]<half)
			break
		endif
	endfor

	// check that width is at least minWid
	if (N<minWid)
		return 1
	elseif((hi-lo)<minWid)
		for (; (hi-lo)<minWid; )
			lo = max(0,lo-1)
			hi = min(N,hi+1)
		endfor
	endif
	return 0
End


Function/S FixUnitsInLabel(units)
	String units
	if (strlen(units)<1)
		return ""
	endif
	return "  "+SelectString(strsearch(units,"1/",0)==0,"(\\U)","(\\E"+units[2,Inf]+"\\S-1\\M)")
End


Static Function MoveWinToTopRight(win,tt,rr)			// returns left edge of window
	String win
	Variable tt,rr													// top right corner of window (use NaN for screen edge)
	Variable left,top,right,bottom
	if (cmpstr(IgorInfo(2),"Windows") == 0)
		GetWindow kwFrameInner, wsize  // get main frame size in Windows
		left = V_left; right = V_right; top = V_top; bottom = V_bottom
	else   // get screen size on Mac
		if (strlen(win))
			win = StringFromList(0,WinList("*",";", "WIN:4183"))
		endif
		String str = StringByKey("SCREEN1", IgorInfo(0))	// RECT=left,top,right,bottom
		Variable i = strsearch(str,"RECT=",0)+5
		sscanf str[i,inf], "%d,%d,%d,%d",left,top,right,bottom	// get screen size
		if (V_flag!=4)
			return NaN
		endif
	endif
	right = !(rr>0) ? right : rr								// desired right side
	top = max(!(tt>52) ? top : tt,45)						// desired top
	GetWindow $win, wsize									// get current window size
	left = max(0,right-V_right+V_left)						// try to preserve width
	bottom = top+V_bottom-V_top
	MoveWindow/W=$win left,top,right,bottom
	return left
End


//	returns the directory of the path in a semi-colon separated list (works for both Mac & Windows)
Static Function/T directory(pathName)
	String pathName							// path name that I want the directory of
	PathInfo $pathName	

	SVAR imageExtension=root:Packages:imageDisplay:imageExtension
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
		ExecuteScriptText/Z cmd
		if (V_Flag==0)
			// print POSIXpath, "        ",cmd
			list = S_value
			Variable i=strlen(list)
			list = list[1,i-2]						// remove leading and trailing double quote
			list = ReplaceString("\r", list, ";" )		// change to a semicolon separated list
		else
			print " the 'ls' command failed, you may want to add anther mount point to the directory() function to speed things up"
			list = IndexedFile($pathName,-1,"."+imageExtension)
		endif
	else											// this section for Windows
		list = IndexedFile($pathName,-1,imageExtension)	// list = IndexedFile($pathName,-1,"????")
	endif
	return list
End



Function/T SortNumericalList(list,sep,increasing)		// return a list of number sorted numerically
	String list								// a list of numbers
	String sep								// separator used for this list
	Variable increasing						// true=increaing order,   false=decreasing order

	Variable i, N = ItemsInList(list,sep)
	Make/N=(N)/D/FREE nw	// need a wave to do sorting
	for (i=0;i<N;i+=1)
		nw[i] = str2num(StringFromList(i,list,sep))
	endfor
	if (increasing)
		Sort nw,nw
	else
		Sort/R nw,nw
	endif

	list = ""
	for (i=0;i<N;i+=1)
		list += num2str(nw[i])+sep
	endfor
	return list
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
		header = ReadGenericHeader(S_path+StringFromList(m,list))
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
			printf "for %s,  range = [%g,  %g],  %s = %g\r",name,V_min,V_max,GDELTA,V_max-V_min
		endif
	endfor
	return tags
End


Static Function CompressEmptyDimensions(Q_Positions,printIt)	// remove empty dimensions from Q_Positions wave (1's & 0's are empty)
	Wave Q_Positions
	Variable printIt
	printIt = !(!printIt)

	if (!WaveExists(Q_Positions))
		printIt = -1
		String list=WaveListClass("QdistAtPositions","*",""), QposName
		if (ItemsInList(list)==1)
			QposName = StringFromList(0,list)
			Wave Q_Positions = $QposName
		elseif (ItemsInList(list)>1)
			Prompt QposName,"Q-distributions array", popup,list
			DoPrompt "Q-distributions",QposName
			if (V_flag)
				return 1
			endif
			Wave Q_Positions = $QposName
		endif
		printf "CompressEmptyDimensions(%s, %g)\r",NameOfWave(Q_Positions),printIt
	endif
	if (!WaveExists(Q_Positions))
		return 1
	endif

	Variable Ndim=WaveDims(Q_Positions)
	Make/N=(Ndim)/FREE idims=DimSize(Q_Positions,p)
	String wnote=note(Q_Positions)

	Make/N=(Ndim)/T/FREE labels="", labelNames={"labelX","labelY","labelZ","labelW"}
	labels = StringByKey(labelNames[p],wnote,"=")		// save labelX, Y, Z, W
	labels = SelectString(p<Ndim,"",labels[p])			// trim off un-needed labels

	Make/N=4/T/FREE/I Nnew=0
	Variable i, m
	for (i=0,m=0; i<Ndim; i+=1)
		if (idims[i]>1)
			labels[m] = labels[i]
			Nnew[m] = idims[i]
			m += 1
		endif
	endfor
	labels = SelectString(p<m,"",labels[p])				// trim off un-needed labels
	Nnew = p<m ? Nnew[p] : 0
	for (i=0;i<Ndim;i+=1)										// re-set the labels
		if (strlen(labels[i]))
		wnote = ReplaceStringByKey(labelNames[i],wnote,labels[i],"=")
		else
			wnote = RemoveByKey(labelNames[i],wnote,"=")
		endif
	endfor

	Duplicate/FREE Q_Positions,Qtemp						// change the dimensions
	Redimension/N=(Nnew[0],Nnew[1],Nnew[2],Nnew[3]) Q_Positions
	Q_Positions = Qtemp
	Note/K Q_Positions,wnote

	if (printIt==-1 || (printIt && !EqualWaves(idims,Nnew,1)))
		printf "changed dimensions of '%s' from %s to %s\r",NameOfWave(Q_Positions),vec2str(idims),vec2str(Nnew)
	endif
	return 0
End


Function reScaleQwaveToStrain(w,d0)
	Wave w						// wave to rescale, assume x-axis is in 1/nm
	Variable d0					// reference d-spacing (nm)
	d0 = !(d0>0) ? NumVarOrDefault("d0",NaN) : d0
	String units
	if (!WaveExists(w) || numtype(d0))
		String wName
		String wlist=WaveList("*",";","")
		Variable i
		for (i=ItemsInLIst(wlist)-1;i>=0;i-=1)
			wName = StringFromList(i,wList)
			units = WaveUnits($wName,0)
			if (!stringmatch(units,"1/nm") && !stringmatch(units,"1/"+ARING))
				wList = RemoveFromList(wName,wList)
			endif
		endfor
		wName = ""
		Prompt wName,"wave to re-scale",popup,wList
		Prompt d0, "d-spacing (nm)"
		DoPrompt "wave to re-scale",wName,d0
		if (V_flag)
			return 1
		endif
		Wave w=$wName
	endif
	if (!WaveExists(w) || !(d0>0))
		return 1
	endif

	units = WaveUnits(w,0)
	if (!stringmatch(units,"1/nm") && !stringmatch(units,"1/"+ARING))
		DoAlert 0,"wave units not 1/nm or 1/"+ARING+", they are '"+units+"'"
		return 1
	endif
	Variable N=numpnts(w), Qstart=DimOffset(w,0), dQ=DimDelta(w,0), Qc=2*PI/d0
	SetScale/P x -(Qstart-Qc)/Qc,-dQ/Qc,"", w
End

// ======================================= End of Q-histograming Utility  ========================================
// ===============================================================================================================

// ============================================ End of Q-histograms  =============================================
// ===============================================================================================================
// ===============================================================================================================
// ===============================================================================================================




// ===============================================================================================================
// ===============================================================================================================
// ===============================================================================================================
// =========================================== Start of Pseudo White  ============================================

Function MakePseudoWhiteImages(pathName,namePart,[range1,range2])
	String pathName			// probably "imagePath"
	String namePart			//	 = "EW5_"
	String range1,range2		// range1=energies, range2=depths, use range1 if only one range

	pathName = SelectString(strlen(pathName),"imagePath",pathName)
	PathInfo $pathName
	String str
	if (!V_flag || strlen(namePart)<1)			// path does not exist or no namePart, ask user
		String pathPart
		str = requestFileRoot(pathName,NaN)
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
		printf "MakePseudoWhiteImages(\"%s\",\"%s\")\r",pathName,namePart
	endif
	PathInfo $pathName
	if (strlen(S_path)<1 || strlen(namePart)<1)
		return 1								// invalid inputs
	endif
	String name,fileRoot=S_path+namePart

	Variable needRangePrompt = ParamIsDefault(range1) || ParamIsDefault(range2)
	if (needRangePrompt)
		str = get_ranges(pathName,namePart)
		range1=StringFromList(0,str)
		range2=StringFromList(1,str)
	endif
	if (ParamIsDefault(range1) || ParamIsDefault(range2) || strlen(range1)==0)
		String range1Prompt = SelectString(strlen(range1)>253,range1,"(range 1 is too long)")
		String range2Prompt = SelectString(strlen(range2)>253,range2,"(range 2 is too long)")
		Prompt range1Prompt,"range of energies, file numbers  (primary)"
		Prompt range2Prompt,"range of depths (or other motion), file numbers  (secondary)"
		DoPrompt "ranges", range1Prompt,range2Prompt
		if (V_flag)
			return 1
		endif
		range1 = SelectString(StringMatch(range1Prompt,"(range 1 is too long)"),range1Prompt,range1)
		range2 = SelectString(StringMatch(range2Prompt,"(range 2 is too long)"),range2Prompt,range2)
	endif
	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EscanButtonProc"))
		printf "MakePseudoWhiteImages(\"%s\",\"%s\",range1=\"%s\",range2=\"%s\")\r",pathName,namePart,range1,range2
		printf "using data from files starting with '%s'\r",namePart
		printf "setting energy range to '%s'  (range1)\r",range1
		printf "setting depth range to '%s'  (range2)\r",range2
	endif

	Variable NkeV = max(ItemsInRange(range1),1)
	Variable Ndepth = max(ItemsInRange(range2),1)

	SVAR imageExtension=root:Packages:imageDisplay:imageExtension
	Variable iDepth=str2num(range2), jkeV=str2num(range1)
	if (!(jkeV>=0))
		sprintf name,"%s%d%s",fileRoot,iDepth,imageExtension
	elseif (!(iDepth>=0))
		sprintf name,"%s%d%s",fileRoot,jkeV,imageExtension
	else
		sprintf name,"%s%d_%d%s",fileRoot,jkeV,iDepth,imageExtension
	endif
	Wave image = $(LoadGenericImageFile(name))
	if (!WaveExists(image))
		printf "could not load first image named '%s'\r",name
		Abort "could not load first image"
	endif
	String wnote = note(image)
	Variable Nx=DimSize(image,0),  Ny=DimSize(image,1)
	Variable tooBig = (Nx*Ny*Ndepth > 1e8)
	tooBig = Ndepth>1 ? tooBig : 1					// only one range, no need to save intermediate sums
	if (tooBig)
		if (Ndepth>1)
			print "There is not enough memory available to make the individual sum images, so no wave 'imageEsum'"
		endif
		Make/N=(Nx,Ny)/O/D imageSumAll
		imageSumAll = 0
	else
		Make/N=(Nx,Ny,Ndepth)/O/I imageEsum
		if (WaveType(image)<=4)				// use double
			Redimension/D imageEsum
		endif
		imageEsum = 0
	endif
	if (Ndepth>1)
		Make/N=(Ndepth)/O/D sumIntensDepth
		SetScale/P x 0,1,"", sumIntensDepth
		sumIntensDepth = 0
	endif

	String depthName = SelectString(numtype(NumberByKey("depth", wnote,"=")),"depth","depthSi")
	Variable depthHi,depthLo = NumberByKey(depthName, wnote,"=")	// depth
	KillWaves/Z imagePlane
	Rename image imagePlane					// need this for the movie
	if (WaveType(imagePlane ) & 0x18)
		Redimension/S imagePlane
	endif

	String progressWin = ProgressPanelStart("",stop=1,showTime=1)	// display a progress bar
	Variable j, i
	for (j=0,jkeV=str2num(range1); j<NkeV && !numtype(jkeV); j+=1,jkeV=NextInRange(range1,jkeV))
		if (ProgressPanelUpdate(progressWin,j/NkeV*100))	// update progress bar
			break												//   and break out of loop
		endif
		for (i=0,iDepth=str2num(range2);i<Ndepth && !numtype(iDepth) || i==0; i+=1,iDepth=NextInRange(range2,iDepth))
			if (!(jkeV>=0))
				sprintf name,"%s%d%s",fileRoot,iDepth,imageExtension
			elseif (!(iDepth>=0))
				sprintf name,"%s%d%s",fileRoot,jkeV,imageExtension
			else
				sprintf name,"%s%d_%d%s",fileRoot,jkeV,iDepth,imageExtension
			endif
			Wave image = $(LoadGenericImageFile(name))
			if (tooBig)
				imageSumAll += image
			else
				imageEsum[][][i] += image[p][q]
			endif
			if (WaveExists(sumIntensDepth))
				sumIntensDepth[i] += sum(image)
			endif

			depthHi = NumberByKey(depthName, note(image),"=")	// depth
			KillWaves/Z image
		endfor
	endfor
	printf "total execution time = %s\r",Secs2Time(SecondsInProgressPanel(progressWin),5,0)
	DoWindow/K $progressWin
	if (!tooBig)
		ImageTransform sumPlanes imageEsum		// sum over all energies
		Duplicate/O M_SumPlanes imageSumAll
		KillWaves/Z M_SumPlanes
		SetScale/P x 0,1,"pixel", imageEsum
		SetScale/P y 0,1,"pixel", imageEsum
		if (numtype(depthHi+depthLo)==0)
			SetScale/I z depthLo,depthHi,Gmu+"m", imageEsum
		endif
	endif
	SetScale/P x 0,1,"pixel", imageSumAll
	SetScale/P y 0,1,"pixel", imageSumAll
	if (numtype(depthHi+depthLo)==0 && WaveExists(sumIntensDepth))
		SetScale/I x depthLo,depthHi,Gmu+"m", sumIntensDepth
	endif
	if (strlen(range1))
		wnote = ReplaceStringByKey("range1",wnote,range1,"=")
	endif
	if (strlen(range2))
		wnote = ReplaceStringByKey("range2",wnote,range2,"=")
	endif
	Note/K imageSumAll, ReplaceStringByKey("waveClass",wnote,"rawImageSum","=")
	if (Ndepth>1)
		MakeEsumPlot()
	else
		KillWaves/Z imagePlane
		NewImageGraph(imageSumAll)
	endif
	beep
	return 0
End
//
Static Function/S FindGraphWithImage(image)	// find the graph window which contains the specified image
	Wave image
	if (!WaveExists(image))
		return ""
	endif
	String name, name0=GetWavesDataFolder(image,2)
	String ilist, win,wlist = WinList("*",";","WIN:1")
	Variable i,Ni, m,Nm=ItemsInList(wlist)
	for (m=0;m<Nm;m+=1)
		win = StringFromList(m,wlist)
		ilist = ImageNameList(win,";")		// list of images in this window
		Ni = ItemsInList(ilist)
		for (i=0;i<ItemsInList(ilist);i+=1)
			name = GetWavesDataFolder(ImageNameToWaveRef(win,StringFromList(i,ilist)),2)
			if (stringmatch(name,name0))
				return win
			endif
		endfor
	endfor
	return ""
End
//
Static Function/T requestFileRoot(pathName,suffixes)
	String pathName
	Variable suffixes					// maximum number of suffixes to remove, each suffix looks like "_12", leave the underscore, use NaN for auto-determine
										// suffixes are ALWAYS preceeded with an "_"
	PathInfo $pathName
	pathName = SelectString(V_flag,"",pathName)
	String message="pick any reconstructed image file in range,  using datafolder = "+GetDataFolder(0)
	Variable refNum
	initEnergyScans()
	SVAR ImageFileFilters=root:Packages:imageDisplay:ImageFileFilters
	Open/F=ImageFileFilters/D/M=message/P=$pathName/R refNum
	if (strlen(S_fileName)<1)
		return ""
	endif
	String fullPath = ParseFilePath(1,S_fileName,":",1,0)
	String name = ParseFilePath(3,S_fileName,":",0,0)
	String fileRoot=name
	if (suffixes<1)						// nothing to do
		return fullPath+";"+fileRoot
	elseif (suffixes>=1)
		Variable i
		for (i=1;suffixes>=1 && i>0;suffixes-=1)	// strip off all but last suffix
			i = strsearch(fileRoot,"_",Inf,-1)
			if (i>0)
				fileRoot = fileRoot[0,i-1]
			endif
		endfor
		return fullPath+";"+fileRoot+"_"

	else									// auto determine the suffixes
		Variable i0,i1
		i1 = strlen(fileRoot)-1			// position of end
		i0 = strsearch(fileRoot, "_",i1,1)
		if (i0<0)
			return fullPath+";"+fileRoot
		endif

		Make/N=(i1-i0)/FREE flags=1	// check if last part is an integer
		flags = !isdigit(fileRoot[p+I0+1])
		if (sum(flags))
			return fullPath+";"+fileRoot
		endif
		fileRoot = fileRoot[0,I0]

		i1 = strlen(fileRoot)-1			// position of end
		i0 = strsearch(fileRoot, "_",i1-1,1)
		if (i0<0)
			return fullPath+";"+fileRoot
		endif

		Make/N=(i1-i0-1)/FREE flags=1	// check if last part is an integer
		flags = !isdigit(fileRoot[p+I0+1])
		if (sum(flags))
			return fullPath+";"+fileRoot
		endif
		fileRoot = fileRoot[0,I0]
		return fullPath+";"+fileRoot
	endif
End

Function MakeEsumMovie()
	Wave imageEsum=imageEsum, imageplane=imageplane, imageSumAll=imageSumAll
	if (!WaveExists(imageEsum) || !WaveExists(imageplane) || !WaveExists(imageSumAll))
		DoAlert 0, "You must run MakePseudoWhiteImages(\"\",\"\") first"
		return 1
	endif

	MakeEsumPlot()
	DoUpdate
	NewMovie/F=10
	SetEsumDepthProc("EsumDepthDisp",-1,"","")
	WaveStats/M=1/Q imageSumAll
	ImageDisplayScaling#ModifyOnly_ctab_range("","imagePlane",V_min,V_max*0.7)
	AddMovieFrame

	WaveStats/M=1/Q imageEsum
	ImageDisplayScaling#ModifyOnly_ctab_range("","imagePlane",0,V_max/3)


	Variable i, N=DimSize(imageEsum,2)
	for (i=0;i<N;i+=1)
		SetEsumDepthProc("EsumDepthDisp",i,"","")
		AddMovieFrame
	endfor
	SetEsumDepthProc("EsumDepthDisp",-1,"","")
	WaveStats/M=1/Q imageSumAll
	ImageDisplayScaling#ModifyOnly_ctab_range("","imagePlane",V_min,V_max*0.7)
	AddMovieFrame
	CloseMovie
	SetEsumDepthProc("EsumDepthDisp",NumVarOrDefault("EsumDepth",-1),"","")
	return 0
End

Function MakeEsumPlot()
	Wave imageEsum=imageEsum, imageplane=imageplane, imageSumAll=imageSumAll
	if (!WaveExists(imageEsum) || !WaveExists(imageplane) || !WaveExists(imageSumAll))
		DoAlert 0, "You must run MakePseudoWhiteImages(\"\",\"\") first"
		return 1
	endif
	String wName = FindGraphWithImage(imagePlane)	
	if (strlen(wName))
		DoWindow/F $wName
	else
		Variable/G EsumDepth=-1
		Display /W=(40,44,736,505)
		AppendImage/T imagePlane
		ModifyImage imagePlane ctab= {*,*,Terrain,1}

		WaveStats/M=1/Q imageSumAll
		ImageDisplayScaling#ModifyOnly_ctab_range("","imagePlane",V_min,V_max*0.7)

		ModifyGraph margin(left)=14,margin(bottom)=14,margin(top)=14,margin(right)=14
		ModifyGraph mirror=2,nticks(left)=3,minor=1,fSize=9,standoff=0
		ModifyGraph tkLblRot(left)=90,btLen=3,tlOffset=-2
		SetAxis/A/R left
		TextBox/N=textEsumDepth/F=0/S=3/A=LT/X=3.14/Y=3.46 "\\F'symbol'\\Zr200S\\Zr050\\F]0(Energys)"
		SetVariable EsumDepthDisp,pos={2,2},size={121,18},proc=SetEsumDepthProc,title="Depth Slice"
		SetVariable EsumDepthDisp,fSize=12
		SetVariable EsumDepthDisp,limits={-1,DimSize(imageEsum,2),1},value=EsumDepth,bodyWidth= 50
		SetWindow kwTopWin hook(arrows)=arrowHook	// this arrow hook allows you to use arrow keys
		String str
		sprintf str, "value=%s;proc=SetEsumDepthProc;stepSize=1;loLim=-1;hiLim=%g",GetDataFolder(1)+"EsumDepth",DimSize(imageEsum,2)
		SetWindow kwTopWin userdata(arrows)=str
	endif
	DoUpdate
	SetAspectToSquarePixels("")
	imageplane = 0
	SetEsumDepthProc("EsumDepthDisp",-1,"","")
	return 0
End
//
Function SetEsumDepthProc(ctrlName,i,varStr,varName) : SetVariableControl
	String ctrlName
	Variable i
	String varStr
	String varName

	Variable depth

	Wave imageplane = ImageNameToWaveRef("","imageplane")
	if (!WaveExists(imageplane))
		return 1
	endif
	String path = GetWavesDataFolder(imageplane,1)
	Wave imageEsum=$(path+"imageEsum"), imageSumAll=$(path+"imageSumAll")
	if (!WaveExists(imageEsum) || !WaveExists(imageSumAll))
		return 1
	endif

	String list
	Variable lo,hi
	if (i<0)									// moving to i==-1, so save current image z limits in window note
		list = GetUserData("","","ctab")
		if (NumberByKey("sum",list,"=")!=1)		// coming off of a regular frame, so store current values
			get_ctab_range("",NameOfWave(imageplane),lo,hi)
			sprintf list,"lo=%g;hi=%g;sum=1",lo,hi
			SetWindow kwTopWin userdata(ctab)=list
		endif
		imageplane = imageSumAll
		WaveStats/M=1/Q imageSumAll
		ImageDisplayScaling#ModifyOnly_ctab_range("","imagePlane",V_min,V_max*0.7)
		TextBox/C/N=textEsumDepth/F=0/S=3/A=LT/X=3.14/Y=3.46 "\\F'symbol'\\Zr200S\\Zr050\\F]0(Energys)   \\F'symbol'\\Zr200S\\Zr050\\F]0(Depths)"
	else
		list = GetUserData("","","ctab")
		if (NumberByKey("sum",list,"="))		// coming off of a sum frame, so reset ctab limits to stored values
			lo = NumberByKey("lo",list,"=")
			hi = NumberByKey("hi",list,"=")
			list = ReplaceNumberByKey("sum",list,0,"=")
			SetWindow kwTopWin userdata(ctab)=list
			ImageDisplayScaling#ModifyOnly_ctab_range("","imagePlane",lo,hi)
		endif
		imageplane = imageEsum[p][q][i]
		depth = DimOffset(imageEsum,2) + i*DimDelta(imageEsum,2)
		TextBox/C/N=textEsumDepth/F=0/S=3/A=LT/X=3.14/Y=3.46 "\\F'symbol'\\Zr200S\\Zr050\\F]0(Energys)  Depth = "+num2str(depth)+" "+Gmu+"m"
	endif
	String title = StrVarOrDefault(path+"title","")
	if (strlen(title))
		AppendText/N=textEsumDepth "\\Zr130"+title
	endif
	AppendText/N=textEsumDepth "\\Zr060"+path
End



//Function aa()
//	print get_ranges("imagePath","EW5_")
//End

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
	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EscanButtonProc"))
		print "range1	",ItemsInRange(range1),"		",range1
		print "range2	",ItemsInRange(range2),"		",range2
	endif
	KillWaves/Z fileList_get_ranges, n1_get_ranges, n2_get_ranges
	return range1+";"+range2
End

Function arrowHook(s)				// for using arrow keys in addition to the little up/down arrows on a SetVariable control
	STRUCT WMWinHookStruct &s
	if (s.eventMod != 8)						// change cursor when command key is held down
		s.cursorCode = 0						// no command key, so set to ordinary cursor
		s.doSetCursor = 0
		return 0
	endif

	s.doSetCursor= 1							// command key is down, change cursor
	s.cursorCode = 22							// vert bar with left-right arrows, or use cursorCode=5
	Variable step=1
	switch(s.keycode)
		case 29:								// right arrow
		case 30:								// up arrow
			step = 1
			break
		case 28:								// left arrow
		case 31:								// down arrow
			step = -1
			break
		default:
			return 0
	endswitch

	String win=s.winName						// get data from window so I know what to do
	String udata = GetUserData(win,"","arrows")
	String funcSpec=StringByKey("proc", udata,"=")
	NVAR value = $StringByKey("value", udata,"=")
	Variable lo=NumberByKey("loLim", udata,"="), hi=NumberByKey("hiLim", udata,"=")
	Variable stepSize = NumberByKey("stepSize", udata,"=")
	if (!NVAR_Exists(value) || exists(funcSpec)!=6 || numtype(stepSize))
		return 0
	endif
	Variable newvalue =  step*stepSize+value
	if (numtype(lo)<2 && numtype(hi)<2)
		newvalue = limit(newvalue,lo,hi)
	endif
	if (numtype(newvalue))
		return 0
	endif
	value = newvalue
	FUNCREF SetEsumDepthProc  arrowProc= $funcSpec
	arrowProc(StringByKey("ctrlName",udata,"="),value,StringByKey("varStr",udata,"="),StringByKey("varName",udata,"="))
	return 1
End

// ============================================ End of Pseudo White  =============================================
// ===============================================================================================================
// ===============================================================================================================
// ===============================================================================================================



// ===============================================================================================================
// ===============================================================================================================
// ===============================================================================================================
// ======================================== Start of 3D Q-space, No Depth ========================================

//	Process many images all at one position, same depth and same x-y (actually X-H) positions.
//	There is no wire scan or any positioners looked at. It is assumed that the images are all from same volume element.
//	No scanning in x, and all at one depth (probably from a thin sample).
// does not require that depth be in image
Function/WAVE Fill1_3DQspace(recipSource,pathName,nameFmt,range,[depth,mask,dark,Qbox, NQx,NQy,NQz, I0normalize,doConvex,FilterFunc,roi,autoGo,printIt])
	Variable recipSource	// 0=beam-line,  1=rceip from an indexation
	String pathName			// either name of path to images, or the full explicit path, i.e. "Macintosh HD:Users:tischler:data:cal:recon:"
	String nameFmt				// the first part of file name, something like  "EW5_%d.h5"
	String range				// range of file indicies to use
	Variable depth				// OPTIONAL, in case image does not contain depth, default is 0
	Wave mask					// OPTIONAL, optional mask to limit the pixels that get processed (use pixel when mask true)
	Wave dark					// OPTIONAL, an optional background wave
	STRUCT boundingVolume &Qbox	// OPTIONAL, do not specify both Qbox AND {NQx,NQy,NQz}
	Variable NQx,NQy,NQz	// OPTIONAL, dmensions of final Qspace3D array
	Variable I0normalize 	// a Flag, if True then normalize data (default is True)
	Variable doConvex			// if True, make the convex hull
	String FilterFunc			// OPTIONAL, name of a filter for the image
	STRUCT imageROIstruct &roi
	Variable autoGo			// OPTIONAL, if true, then do not ask user about Nx,Ny,Nz, just use the auto values
	Variable printIt
	depth = ParamIsDefault(depth) ? 0 : depth
	depth = numtype(depth) ? NaN : depth
	NQx = ParamIsDefault(NQx) || NQx<=0 ? NaN : NQx		// default to NaN, the program figures it out
	NQy = ParamIsDefault(NQy) || NQy<=0 ? NaN : NQy
	NQz = ParamIsDefault(NQz) || NQz<=0 ? NaN : NQz
	I0normalize = ParamIsDefault(I0normalize) || numtype(I0normalize) ? 1 : I0normalize
	doConvex = ParamIsDefault(doConvex) || numtype(doConvex) ? 0 : !(!doConvex)
	FilterFunc = SelectString(ParamIsDefault(FilterFunc),FilterFunc,"")
	autoGo = ParamIsDefault(autoGo) || numtype(autoGo) ? 0 : !(!autoGo)
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : !(!printIt)
	String extras=""
	if (!ParamIsDefault(roi))						// an roi was given, check to see if it is empty or invalid
		if (!imageROIstructBad(roi) && !(roi.empty))
			sprintf extras,"roi=%d,%d,%d,%d;",roi.xLo, roi.xHi, roi.yLo, roi.yHi
		endif
	endif
	if (!(recipSource==0 || recipSource==1))
		recipSource = NumVarOrDefault("root:Packages:micro:Escan:recipSource",NaN)
		Prompt recipSource, "source of reciprocal lattice", popup, "Beamline Coords;Indexing Result"
		recipSource += 1
		DoPrompt "Reciprocal Lattice",recipSource
		if (V_flag)
			return $""
		endif
		recipSource -= 1
		printIt = 1
	endif
	if (!(recipSource==0 || recipSource==1))
		return $""
	endif
	String recipLattice=""
	if (recipSource==1)
		Make/N=(3,3)/D/FREE recip=(p==q)
		String indexList=WaveListClass("IndexedPeakList","*","")
		if (ItemsInList(indexList)==1)
			Wave FullPeakIndexed=$StringFromList(0,indexList)
		elseif (ItemsInList(indexList)>1)
			String indexName
			Prompt indexName,"Source for reciprocal lattice",popup,indexList
			DoPrompt "Reciprocal Lattice",indexName
			if (V_flag)
				return $""
			endif
			printIt = 1
			Wave FullPeakIndexed=$indexName
		endif
		if (WaveExists(FullPeakIndexed))			// if FullPeakIndexed exists, get the indexed reciprocal lattice from it
			String indexNote=note(FullPeakIndexed), recipList="", item
			Variable i
			for (i=0;;i+=1)								// find all reciprocal latticies in wave note
				item = StringByKey("recip_lattice"+num2istr(i),indexNote,"=")
				if (strlen(item)<1)
					break
				endif
				recipList += "recip_lattice"+num2istr(i)+";"
			endfor
			if (ItemsInList(recipList)==1)			// get reciprocal lattice from the wave note
				item = StringFromList(0,recipList)
			elseif (ItemsInList(recipList)>1)
				Prompt item,"reciprocal lattice",popup,recipList
				DoPrompt "Reciprocal Lattice",item
				if (V_flag)
					return $""
				endif
				printIt = 1
			else
				item = ""
			endif
			if (strlen(item)>0)							// save reciprocal lattice for latter inclusion into another wave note
				recipLattice = item +"="+ StringByKey(item,indexNote,"=")
			endif
		endif
	elseif (recipSource==0)
		recipLattice = ""
		if (printIt)
			print "No reciprocal lattice given, just using beam line coordinates"
		endif
//		DoAlert 0,"Not yet implemented use of measured reciprocal lattice"
//		return $""
	else
		DoAlert 0,"recipSource must be only 0 or 1"
		return $""
	endif
	NewDataFolder/O root:Packages
	NewDataFolder/O root:Packages:micro
	NewDataFolder/O root:Packages:micro:Escan
	Variable/G root:Packages:micro:Escan:recipSource=recipSource

	String str
	PathInfo $pathName
	if (!V_flag || strlen(nameFmt)<1)				// path does not exist or no nameFmt, ask user
		String pathPart
		str = requestPathFileFmt(pathName)
		pathPart = StringFromList(0,str)
		nameFmt = StringFromList(1,str)
		if (strlen(pathPart)<1 || strlen(nameFmt)<1)
			return $""										// invalid inputs
		endif
		if (!stringmatch(pathPart,S_path))			// path was changed
			if (stringmatch(pathName,"imagePath"))
				NewPath/O/M="path to reconstructed image files" imagePath pathPart	// for imagePath, automatically reassign
			else
				NewPath/M="path to reconstructed image files" $pathName pathPart		// for other names, ask
			endif
		endif
		printIt = 1
	endif
	PathInfo $pathName
	if (strlen(S_path)<1 || strlen(nameFmt)<1)
		return $""											// invalid inputs
	endif
	if (printIt)
		printf "using data from files starting with '%s'\r",S_path+nameFmt
	endif

	String list = WaveListClass("imageMask","*","DIMS:2,BYTE:1")	// setup to use optional mask
	String maskName = ""
	if (WaveExists(mask))
		maskName = GetWavesDataFolder(mask,2)
	elseif (strlen(list) && ParamIsDefault(mask))
		maskName = ""
		Prompt maskName, "mask to use with image",popup,"-none-;"+list
		DoPrompt "mask wave",maskName
		if (V_flag)
			return $""
		endif
		maskName = SelectString(stringmatch(maskName,"-none-"),maskName,"")
		Wave mask = $maskName							// do not check if wave exists, that a valid option
		printIt = 1
	endif
	maskName = SelectString(WaveExists(mask),"$\"\"",maskName)
	String darkList = WaveListClass("imageDark;rawImageDark","*","")
	if (!WaveExists(dark) && ItemsInList(darkList)>0 && ParamIsDefault(dark))
		String darkName
		Prompt darkName,"Background Image",popup,"-none-;"+darkList
		DoPrompt "Background",darkName
		if (V_flag)
			return $""
		endif
		if (cmpstr(darkName,"-none-"))
			Wave dark = $darkName
		else
			Wave dark = $""
		endif
		printIt = 1
	endif

	String filterFuncList=FunctionList("*_ImageFilter",";","KIND:2,NPARAMS:1,VALTYPE:1")
	filterFuncList = RemoveFromList("Proto_ImageFilter",filterFuncList)
	if (ParamIsDefault(FilterFunc) && strlen(FilterFunc)==0 && ItemsInList(filterFuncList))
		Prompt FilterFunc,"Filter Function for each image",popup, "-none-;"+filterFuncList
		DoPrompt "Filter for Image",FilterFunc
		if (V_flag)
			return $""
		endif
	endif
	FilterFunc = SelectString(StringMatch(FilterFunc,"-none-"),FilterFunc,"")	// remove "-none-"
	FilterFunc = SelectString(StringMatch(FilterFunc,"Proto_ImageFilter"),FilterFunc,"")	// remove proto func too.
	FUNCREF Proto_ImageFilter Local_ImageFilter = $FilterFunc		// if FilterFunc="", then the proto will be run.

	String progressWin = ProgressPanelStart("",stop=1,showTime=1,status="Determining the range to scan")	// display a progress bar
	if (ItemsInRange(range)<1) 						// if range is empty, get the full range from the directory
		range = get_FilesIndexRange(pathName,nameFmt)
		Prompt range,"range of image file numbers to use"
		DoPrompt "range",range
		if (V_flag)
			return $""
		endif
		if (ItemsInRange(range)<1)
			if (printIt)
				printf "range of images indicies is '%s'\r",range[0,370]	// cannot print lines longer than 400 chars
			endif
		endif
		printIt = 1
	endif

	String fileRootFmt									// complete path up to and including the underscore
	PathInfo $pathName
	if (V_flag)												// if path exists, reset pathName to the explicit path, otherwise assume it is the explicit path
		fileRootFmt = S_path
	elseif (strsearch(pathName,":",Inf,1) == strlen(pathName)-1)
		fileRootFmt = pathName+":"					// ensure terminating colon
	endif
	fileRootFmt += nameFmt
	if (printIt)
		sprintf str,"Fill1_3DQspace(%g,\"%s\",\"%s\",\"%s\"",recipSource,pathName,nameFmt,range
		if (WaveExists(mask))
			maskName = NameOfWave(mask)
			str += ",mask="+maskName
		endif
		if (WaveExists(dark))
			darkName = NameOfWave(dark)
			str += ",dark="+darkName
		endif
		str += SelectString(strlen(FilterFunc), "", ",FilterFunc=\""+FilterFunc+"\"")
		str += SelectString(ParamIsDefault(NQx), ",NQx="+num2str(NQx), "")
		str += SelectString(ParamIsDefault(NQy), ",NQy="+num2str(NQy), "")
		str += SelectString(ParamIsDefault(NQz), ",NQz="+num2str(NQz), "")
		str += ")"
		print str[0,390]
	endif
	if (WaveExists(mask))
		if (sum(mask)==0)
			ERROR_Fill_Q_Positions("You picked a mask that is all zero, stopping",progressWin)
			return $""
		endif
	endif
	if (ItemsInRange(range)<1) 						// if range is empty, get the full range from the directory
		ERROR_Fill_Q_Positions("range is empty",progressWin)
		return $""
	endif

	Variable timer0 = stopMSTimer(-2)/1e6			// starting time, used to time prooess (sec)
	Variable useDistortion = NumVarOrDefault("root:Packages:geometry:useDistortion",USE_DISTORTION_DEFAULT)
 	if (printIt)
		printf "processing with distotion %s", SelectString(useDistortion,"OFF","on")
		print "     "+SelectString(WaveExists(mask),"and no mask","and using  '"+NameOfWave(mask)+"'  for the mask")
	endif

	// open the first file in the range, to get info about the images
	ProgressPanelUpdate(progressWin,0,status="setting up")
	String name
	sprintf name, fileRootFmt, str2num(range)
	String wnote = ReadGenericHeader(name)					// wave note of the image file
	if (strlen(wnote)<1)
		sprintf str, "could not load very first image header from '%s'\r",name
		ERROR_Fill_Q_Positions(str,progressWin)
		return $""
	endif
	STRUCT imageROIstruct roiAll									// an ROI that will contain all of the images being processed
	imageROIstructInit(roiAll, wnote=wnote)					// initialie roiAll with roi of first image
	
	String MonoMode = StringByKey("MonoMode",wnote,"=")// energy for this image
	Variable mono = StringMatch(MonoMode,"mono*")			// if not mono, then must be undulator scan
	FUNCREF gap2keV_A gap2keV_Func = $SelectString(exists("gap2keV")==6,"gap2keV_A","gap2keV")
	Variable keV, gap
	keV = NumberByKey("keV", wnote,"=")						// energy for this image
	gap = NumberByKey("undulatorGap", wnote,"=")			// undulator gap for this image
	keV = mono ? keV : gap2keV_Func(gap)						// convert undulator gap to energy for undulator scan
	if (numtype(keV))
		ERROR_Fill_Q_Positions("invalid keV in image '"+name+"'",progressWin)
		return $""
	endif

	// read in the badPixels (if they exist) and combine with mask (1=OK, 0=bad)
	Wave badPixelsAll = GetBadPixelsImage(wnote)			// get full image with bad pixels if badPixels exists
	if (WaveExists(badPixelsAll))
		if ( round(NumberByKey("xDimDet",wnote,"=")) != DimSize(badPixelsAll,0) || round(NumberByKey("yDimDet",wnote,"=")) != DimSize(badPixelsAll,1) )
			ERROR_Fill_Q_Positions("badPixelsAll and image have different dimensions than the whole detector",progressWin)
			return $""
		endif
	endif

	STRUCT microGeometry geo
	FillGeometryStructDefault(geo)
	Variable dNum = detectorNumFromID(geo, StringByKey("detectorID", wnote,"="))
	if (!(dNum>=0 && dNum<MAX_Ndetectors))
		ERROR_Fill_Q_Positions("could not get detector number from from wave note using detector ID",progressWin)
		return $""
	endif
	String wnoteFull = wnote

	// read header from each of the images, and store it to figure out what was done.
	Variable N = ItemsInRange(range)							// number of images to be processed
	Make/N=(N)/FREE/D keV_FillQvsPositions=NaN				// hold the Energies
	Make/N=(N)/FREE/U/I m_Fill_QHistAt1Depth=0
	Make/N=(N)/FREE/D Ic0_FillQvsPositions=NaN				// hold the I0 for normalizations
	Make/N=(N)/FREE/S currents=NaN, normalization=NaN	// save for possible normalizations
	// for all the N files (go over range), store the energies, X position, H position, and indicies
	ProgressPanelUpdate(progressWin,0,status="examining "+num2istr(N)+" image headers",resetClock=1)
	Variable m, j, varyROI=0, iv, Npixels=0
	for (m=str2num(range), i=0; numtype(m)==0; m=NextInRange(range,m), i+=1)		// loop over range
		if (mod(m,20)==0)
			if (ProgressPanelUpdate(progressWin,i/N*100))	// update progress bar
				ERROR_Fill_Q_Positions("User abort",progressWin)
				return $""
			endif
		endif
		sprintf name, fileRootFmt, m
		wnote=ReadGenericHeader(name, extras=extras)		// wave note to add to file read in

		iv = ExtendROIfromWaveNote(roiAll,wnote)			// find ROI that holds ALL the images
		if (iv<0)
			continue														// really bad ROI in wave note, skip this one
		endif
		varyROI = iv ? 1 : varyROI								// set varyROI to true if roiAll is changing

		Npixels += NumberByKey("xdim", wnote,"=") * NumberByKey("ydim", wnote,"=")
		if (mono)
			keV = NumberByKey("keV", wnote,"=")				// list of energies
		else
			gap = NumberByKey("undulatorGap", wnote,"=")	// undulator gap for this image
			keV = gap2keV_Func(gap)
		endif
		keV_FillQvsPositions[i] = keV							// list of energies
		m_Fill_QHistAt1Depth[i] = m								// file id number for each image
		Ic0_FillQvsPositions[i] = NumberByKey("I0", wnote,"=")	// Ic0 for normalization
		currents[i] = NumberByKey("ringCurrent",wnote,"=")
	endfor
	if (varyROI && WaveExists(mask))
		ERROR_Fill_Q_Positions("You cannot use a mask when you are also varying the ROI size", progressWin)
		return $""
	elseif (WaveExists(dark) && varyROI)
		ERROR_Fill_Q_Positions("You cannot use a dark image when you are also varying the ROI size", progressWin)
		return $""
	endif

	if (WaveExists(badPixelsAll))
		Wave badPixels = ExtractROIofImage(badPixelsAll, roiAll)	// badPixels is now a free wave, with size of roiAll
	else
		Wave badPixels = $""
	endif
	// if varyROI is True, then maskLocal will be size of roiAll
	if (WaveExists(badPixels) && WaveExists(mask))
		Duplicate/FREE mask, maskLocal
		maskLocal = badPixels || mask							// maskLocal is combination of mask and badPixels
	elseif (WaveExists(badPixels))
		Wave maskLocal = badPixels								// badPixels is already a local free wave
	elseif (WaveExists(mask))
		Wave maskLocal = mask
	else
		Wave maskLocal = $""
	endif

	Make/N=(N,3)/O/D Qspace3DIndexEnergy
	Qspace3DIndexEnergy[][0] = m_Fill_QHistAt1Depth[p]	// permanently store the index & energy for evey image
	Qspace3DIndexEnergy[][1] = keV_FillQvsPositions[p]
	Qspace3DIndexEnergy[][2] = Ic0_FillQvsPositions[p]	// also store I0 for reference
	Variable useNormalization=0
	WaveStats/Q/M=1 currents
	if (V_numNans==0 && V_numINFs==0 && V_min>1)			// valid currents, so create normalization[]
		V_avg = roundSignificant(V_avg,3)
		normalization = V_avg / currents						// scale to value of average current
		useNormalization = 1
		if (printIt)
			printf "Normalizing by the beam current, scaled to %g mA\r",V_avg
		endif
	endif
	WaveClear currents
	Variable seconds = stopMSTimer(-2)/1e6 - timer0
	if (printIt && seconds>0.2)
		printf "first pass to examine headers took %s\r",Secs2Time(seconds,5,3)
	endif
	ProgressPanelUpdate(progressWin,0,status="done with headers",resetClock=1)
	if (ParamIsDefault(depth) || numtype(depth))
		depth = NumberByKey("depth",wnoteFull,"=")			// no depth given, try from file
	endif

	Variable dkeV = FindStepSizeInVec(keV_FillQvsPositions, 2e-4)	// abs(step size) in the energy scan, threshold=0.2 eV
	Variable ask = (printIt || strlen(GetRTStackInfo(2))==0)	// ask if called from command line
	ask = ask && !autoGo											// don't ask when autoGo=True
	ask = ask || numtype(dkeV)>0									// always ask for invalid dkeV
	if (ask)
		// correct the scan ranges as necessary
		Prompt depth,"depth to use (micron)"
		Prompt dkeV,"step size in E (eV)"
		Prompt I0normalize,"Normalize the Q-hists",popup,"Do NOT Normalize to Io or exposure;Normalize to Io & exposure"
		Prompt doConvex,"Also compute Convex Hull (can take a while)",popup,"No Convex Hull;Compute Convex Hull"
		doConvex += 1
		dkeV *= 1000
		I0normalize += 1
		DoPrompt "scan sizes",dkeV,depth,doConvex,I0normalize
		if (V_flag)
			DoWindow/K $progressWin								// done with status window
			return $""
		endif
		I0normalize -= 1
		dkeV = abs(dkeV/1000)										// sign is not used
		doConvex -= 1
	endif
	depth = numtype(depth) ? 0 : depth							// 0 is better than NaN
	I0normalize = numtype(I0normalize) ? 1 : I0normalize	// default is to normalize
	if (printIt && I0normalize)
		print "Normalizing to Io & exposure"
	endif
	WaveStats/M=1/Q keV_FillQvsPositions
	Variable ikeVlo=V_minloc, ikeVhi=V_maxloc				// save location of min and max energies
	if (V_numNans)
		DoAlert 0, "There were "+num2istr(V_numNans)+" bad images found, they will be skipped"
	endif

	// note that Q range only depends upon image size and energy range, not on X or H position
	// determine dQ (1/nm), the Q resolution to use.  Base it on the distance between two adjacent pixels
	Variable Elo=keV_FillQvsPositions[ikeVlo], Ehi=keV_FillQvsPositions[ikeVhi]
	Variable px,py														// the current pixel to analyze (unbinned full chip pixels)
	px = (roiAll.xLo + roiAll.xHi)/2							// approximate full chip pixel position in center or the image (UN-binned)
	px = (roiAll.yLo + roiAll.yHi)/2
	// choose dQ as largest of pixelX+1, pixelY+1, or dEnergy
	Variable dQ = 4*PI*abs(sin(pixel2q(geo.d[dNum],px,py,$""))-sin(pixel2q(geo.d[dNum],px+(roiAll.binx),py,$"")))*Elo/hc
	dQ = max(dQ, 4*PI*abs(sin(pixel2q(geo.d[dNum],px,py,$""))-sin(pixel2q(geo.d[dNum],px,py+(roiAll.biny),$"")))*Elo/hc )
	dQ = max(dQ, 4*PI*abs(sin(pixel2q(geo.d[dNum],px,py,$""))) * dkeV/hc )
	if (printIt)
		printf "E range = [%g, %g] (keV),  %sE=%.2g eV", keV_FillQvsPositions[ikeVlo], keV_FillQvsPositions[ikeVhi],GDELTA,dkeV
		printf "    %s Scan\r",SelectString(mono,"Undulator","Monochromator")
	endif

	seconds = stopMSTimer(-2)/1e6 - timer0
	if (printIt && seconds>0.2)
		printf "setting up arrays took %s\r",Secs2Time(seconds,5,3)
	endif
	// done with the setup part, now actually compute something

	if (useDistortion)												// if using the distortion, precompute for all images  here
		Abort "Fill1_3DQspace(), filling the distortion map needs serious fixing"
		Wave DistortionMap = GetDistortionMap(roiAll)		// if using the distortion, precompute for all images here
	endif																	// distortion map now ready

	ProgressPanelUpdate(progressWin,0,status="pre-computing Qvecs array, may take 90 sec.",resetClock=1)
	// make an array the same size as an roiAll, but filled with Q's at 1keV for this depth, |Q| = 4*¹*sin(theta) * E/hc
	// print "roiAll =",imageROIstruct2str(roiAll)
	Wave Qvecs1keVall = MakeQvecsArray(roiAll,geo.d[dNum],wnoteFull,recip,Elo,Ehi,mask=maskLocal,depth=depth,printIt=printIt)
	// print "Qvecs1keVall",NumberByKey("startx",note(Qvecs1keVall),"="), NumberByKey("endx",note(Qvecs1keVall),"="), NumberByKey("starty",note(Qvecs1keVall),"="), NumberByKey("endy",note(Qvecs1keVall),"="),"  ",DimSize(Qvecs1keVall,0), DimSize(Qvecs1keVall,1)

	// make 3D volume to hold  the 3D q-histogram
	// find range of Qx,Qy,Qz
	wnote = note(Qvecs1keVall)
	Variable QxLo = NumberByKey("QxLo",wnote,"="), QxHi = NumberByKey("QxHi",wnote,"=")
	Variable QyLo = NumberByKey("QyLo",wnote,"="), QyHi = NumberByKey("QyHi",wnote,"=")
	Variable QzLo = NumberByKey("QzLo",wnote,"="), QzHi = NumberByKey("QzHi",wnote,"=")
	if (printIt)
		printf "Q data range = [%g, %g],  [%g, %g],  [%g, %g]\r",QxLo,QxHi, QyLo,QyHi, QzLo,QzHi
	endif

	if (!ParamIsDefault(Qbox))									// Qbox was passed, use it. Qbox overrides any passed NQx, NQy, or NQz
		if (QxLo > Qbox.xhi || QxHi < Qbox.xlo || QyLo > Qbox.yhi || QyHi < Qbox.ylo || QzLo > Qbox.zhi || QzHi < Qbox.zlo)
			print "Requested Q-volume (Qbox) does not contain any data"
			printf "Q Histogram range = [%g, %g],  [%g, %g],  [%g, %g]\r",QxLo,QxHi, QyLo,QyHi, QzLo,QzHi
			print "Requested Qbox = ",boundingVolumeStruct2str(Qbox)
			DoWindow/K $progressWin								// done with status window
			return $""
		endif
		QxLo = Qbox.xlo ;		QxHi = Qbox.xhi
		QyLo = Qbox.ylo ;		QyHi = Qbox.yhi
		QzLo = Qbox.zlo ;		QzHi = Qbox.zhi
		NQx = Qbox.Nx
		NQy = Qbox.Ny
		NQz = Qbox.Nz

	elseif (numtype(NQx+NQy+NQz))								// don't have valid NQx,NQy,NQz
		dQ *= 4
		NQx = floor((QxHi-QxLo)/dQ)								// calculate size of 3D histogram
		NQy = floor((QyHi-QyLo)/dQ)
		NQz = floor((QzHi-QzLo)/dQ)
		Variable NQ = NQx*NQy*NQz
		Variable Ntot=NQx*NQy*NQz
		if (Ntot>1e7)													// Q volume has too many voxels, coarsen so only 1e7 voxels
			Variable factor = 215/(Ntot^0.333333)				// 215 is about 3rd root of 1e7
			NQx = round(factor*NQx)
			NQy = round(factor*NQy)
			NQz = round(factor*NQz)
			Ntot = NQx*NQy*NQz
		endif

		if (!autoGo)
			Prompt NQx,"no. of Qx points in 3D Q histogram"
			Prompt NQy,"no. of Qy points in 3D Q histogram"
			Prompt NQz,"no. of Qz points in 3D Q histogram"
			DoPrompt "3D Q histogram size",NQx, NQy, NQz
			if (V_flag)
				DoWindow/K $progressWin
				return $""
			endif
			printIt = 1
		endif
		NQx = round(NQx)
		NQy = round(NQy)
		NQz = round(NQz)
		if (printIt)
			printf "NQ = {%g, %g, %g},   (NQx*NQy*NQz) = %g\r",NQx,NQy,NQz, NQ
		endif
	endif

	if (!(NQx>0 && NQy>0 && NQz>0))								// test for valid NQ's
		ERROR_Fill_Q_Positions("Volume has a zero dimension",progressWin)
		return $""
	elseif (printIt)
		printf "Q Histogram range = [%g, %g],  [%g, %g],  [%g, %g]\r",QxLo,QxHi, QyLo,QyHi, QzLo,QzHi
		printf "Make array Qspace3D[%g][%g][%g] to hold Q-space\r",NQx,NQy,NQz
	endif

	Make/N=3/D/FREE Qc={(QxLo+QxHi)/2, (QyLo+QyHi)/2, (QzLo+QzHi)/2}
	Make/N=(NQx, NQy, NQz)/O Qspace3D=0						// array to hold Q-space histogram
	Make/N=(NQx, NQy, NQz)/FREE/I/U Qspace3DNorm=0		// used for normalizing Qspace3D
	String Qunits = SelectString(recipSource,"1/nm","rlu")
	SetScale/I x QxLo,QxHi,Qunits, Qspace3D, Qspace3DNorm
	SetScale/I y QyLo,QyHi,Qunits, Qspace3D, Qspace3DNorm
	SetScale/I z QzLo,QzHi,Qunits, Qspace3D, Qspace3DNorm

	STRUCT imageROIstruct ROIQvecs1keV
	ROIQvecs1keV.empty = 1											// roi of current ROIQvecs1keV, starts empty to forces a calculation first time

	// Starting bulk of processing"
	// for all the N files (go over range), accumulate Qhist from all images
	ProgressPanelUpdate(progressWin,0,status="processing "+num2istr(N)+" images",resetClock=1)	// update progress bar
	Wave i0=$""															// "previous" image in undulator scans
	Variable lowThresh=15, maxPixel=NaN						// only used by undulator gap scan
	Variable ipnt, timer3=stopMSTimer(-2)/1e6, skipUpdate=1	// was 10
	for (m=str2num(range),ipnt=0; numtype(m)==0; m=NextInRange(range,m),ipnt+=1)	// loop over range, all the images
		if (mod(ipnt,skipUpdate)==0)
			if (ProgressPanelUpdate(progressWin,ipnt/N*100))		// update progress bar every 100 images
				break														//   and break out of loop
			endif
		endif
		sprintf name, fileRootFmt, m
		Wave image = $(LoadGenericImageFile(name,extras=extras))
		if (!WaveExists(image))
			printf "could not load image named '%s'\r",name
			continue
		endif
		Local_ImageFilter(image)

		if (!imageMatchesROI(image,ROIQvecs1keV))		// this image has a new roi, so re-set Qvecs1keV & maskLocalSub
			//print "\r image,",NumberByKey("startx",note(image),"="), NumberByKey("endx",note(image),"="), NumberByKey("starty",note(image),"="), NumberByKey("endy",note(image),"=")
			imageROIstructInit(ROIQvecs1keV, wnote=note(image))	// re-set roi_Qvecs1keV to match current image
			//print "ROI struct = ",imageROIstruct2str(ROIQvecs1keV)
			Wave Qvecs1keV = ExtractROIofImage(Qvecs1keVall, ROIQvecs1keV)
			//print "Qvecs1keV",NumberByKey("startx",note(Qvecs1keV),"="), NumberByKey("endx",note(Qvecs1keV),"="), NumberByKey("starty",note(Qvecs1keV),"="), NumberByKey("endy",note(Qvecs1keV),"=")
			Wave maskLocalSub = ExtractROIofImage(maskLocal, ROIQvecs1keV)
		endif


if ( Dimsize(image,0)!=DimSize(Qvecs1keV,0) || Dimsize(image,1)!=DimSize(Qvecs1keV,1) )
Debugger
endif


		// monochromator scan (not undulator gap)
		if (mono)
			wnote = Fill3DQhist1image(image,Qvecs1keV,Qspace3D,Qspace3DNorm,mask=maskLocalSub,dark=dark)// fill Qhist from one file, Qhats are precomputed
			KillWaves/Z image

		// undulator scan (monochromator not used)
		elseif (!WaveExists(i0))									// FIRST image cannot be processed since there is no "previous" image yet
			Wave i0 = image											// save for use as "previous" image in undulator scans
			maxPixel = maxValueOfType(image)					// set max value that a pixel can be
			if (useNormalization)									// using normalized images (by beam current)
				Redimension/S i0										// change i0 to float
				i0 *= normalization[0]
			endif
		else																// SUBSEQUENT, undulator gap scan, use difference between this image and "previous" image
			if (useNormalization)									// valid normalization, so normalize
				Redimension/S image
				image *= normalization[ipnt]
			endif
			// dimage is (image-i0), but zero out any pixels greater than maxPixel or pixels in image less than lowThresh
			// then clip dimage to the range [0,maxPixel]
			MatrixOP/FREE dimage = clip((image-i0) * (greater(maxPixel,i0) && greater(maxPixel,image) && greater(image,lowThresh)),0,maxPixel)
			keV = gap2keV_Func(NumberByKey("undulatorGap",note(image),"="))
			Note/K dimage, ReplaceNumberByKey("keV",note(image),keV,"=")
			wnote = Fill3DQhist1image(dimage,Qvecs1keV,Qspace3D,Qspace3DNorm,mask=maskLocalSub,dark=dark,I0normalize=I0normalize,printIt=printIt)// fill Qhist from one file, Qhats are precomputed
			KillWaves/Z i0												// delete old "previous" image
			Wave i0 = image											// save current image as i0, the new "previous" image
		endif
	endfor
	KillWaves/Z i0, image											// delete left overs
	WaveStats/M=1/Q Qspace3D										// check for empty voume
	Variable allZero = (V_min==0 && V_max==0)
	Variable sec3 = stopMSTimer(-2)/1e6 - timer3
	Qspace3D = Qspace3D / Qspace3DNorm							// do the normalization
	Qspace3D = Qspace3DNorm ? Qspace3D : 0					// fix divide by zeros
	TrimZerosOff3D(Qspace3D)										// trim off zeros on outside of histogram

	Variable NusedPixels = sum(Qspace3DNorm)				// number of pixels that contributed to this histogram
	Qspace3DNorm = Qspace3DNorm ? 1 : 0
	Variable NusedVoxels = sum(Qspace3DNorm)				// number of voxels in Qspace that were filled (some don't get used)
	WaveClear Qspace3DNorm

	if (useDistortion)
		Note/K DistortionMap, "use=0"							// turn off distortion map
	endif

	seconds = stopMSTimer(-2)/1e6 - timer0
	if (printIt)
		printf "Processed %d pixels into %d voxels\r",NusedPixels,NusedVoxels
		printf "\r  processing all %d images took %s	(%.3g %ss/pixel)\r",N,Secs2Time(seconds,5,1),1e6*seconds/(N*Npixels),Gmu
		printf "		the accumulation/assignment part took %s\r",Secs2Time(sec3,5,2)
	endif

	wnote = wnoteFull
	wnote = ReplaceStringByKey("waveClass",wnote,"GizmoXYZ,Qspace3D,Histogram3D","=")
	wnote = ReplaceNumberByKey("depth",wnote,depth,"=")
	wnote = ReplaceStringByKey("Qcenter",wnote,vec2str(Qc,bare=1,sep=","),"=")
	wnote = ReplaceStringByKey("range",wnote,range,"=")
	wnote = ReplaceStringByKey("nameFmt",wnote,nameFmt,"=")
	wnote = ReplaceNumberByKey("NusedPixels",wnote,NusedPixels,"=")
	wnote = ReplaceNumberByKey("NusedVoxels",wnote,NusedVoxels,"=")
	wnote = ReplaceNumberByKey("I0normalize",wnote,I0normalize,"=")
	if (WaveExists(mask))
		wnote = ReplaceStringByKey("maskWave",wnote,GetWavesDataFolder(mask,2),"=")
	endif
	if (strlen(FilterFunc))
		wnote = ReplaceStringByKey("FilterFunc",wnote,FilterFunc,"=")
	endif
	if (strlen(recipLattice)>0)
		wnote += recipLattice+";"
	endif
	wnote = ReplaceNumberByKey("executionSeconds",wnote,seconds,"=")
	Note/K Qspace3D, wnote
	if (seconds>(12*60))											// execution took more than 12 minutes, save experiment
		print "Processing took more than 12 minutees, Saving this Igor Experiemnt"
		SaveExperiment
	endif
	DoWindow/K $progressWin										// done with status window
	if (printIt)
		beep
	endif

	// now display Q-space volume as a Gizmo
	if (allZero)
		printf "'%s' is ALL zeros, nothing to display\r",NameOfWave(Qspace3D)
	elseif (printIt)
		MakeGizmocubeCorners(Qspace3D)
		if (doConvex)
			QspaceVolumesView#MakeConvexHullFrom3D(Qspace3D)
			beep
		endif

		String wName=StringFromlist(0,WindowsWithWave(Qspace3d,4))
		if (strlen(wName)==0)
			MakeGizmoQspace3d(Qspace3D)
		endif
	endif
	return Qspace3d
End
//
Function Proto_ImageFilter(image)				// proto for a filter run on each image before histograming it
	Wave image
	return 0
End
//
Function Median_ImageFilter(image)				// this one is good for removing random hot pixels
	Wave image
	MatrixFilter/N=3 median image				// does a 3x3 median filter over the whole image
	return 0
End
//
Static Function/S Fill3DQhist1image(image,Qvecs1keV,Qhist,QhistNorm,[mask,dark,I0normalize,printIt])// fill Qhist from one file, F is for fast, sin(theta) is precomputed
	Wave image
	Wave Qvecs1keV										// pre-computed array of sin(theta)'s that is same size as the image to be loaded
	Wave Qhist											// wave of Q vectors (at 1keV), accumulate into here
	Wave QhistNorm										// array that contains number of pixels contributing to each bin in Qhist
	Wave mask											// optional mask to limit the pixels that get processed (use pixel when mask true)
	Wave dark											// optional background wave
	Variable I0normalize							// optional Flag, if True then normalize data (default is True)
	Variable printIt
	printIt = ParamIsDefault(printIt) ? NaN : printIt
	printIt = numtype(printIt) ? ItemsInList(GetRTStackInfo(0))<3 : !(!printIt)
	I0normalize = ParamIsDefault(I0normalize) || numtype(I0normalize) ? 1 : I0normalize	// to Normalize

	Variable timer = stopMSTimer(-2)/1e6			// start timer for initial processing
	if (!WaveExists(image))
		return ""
	endif
	Variable Ni=DimSize(image,0),Nj=DimSize(image,1)

	if (WaveExists(mask))								// mask is array that was passed (or not passed)
		Wave maskIn=mask									// there is ALWAYS a maskIn, it is either equal to mask or "1"
	else
		Make/N=(Ni,Nj)/FREE/B/U maskLocal=1
		Wave maskIn=maskLocal							// no mask, so make maskLocal and set to 1 (use all pixels)
	endif
	if (WaveExists(dark))
		if (WaveType(image) & 0x58)					// either 8bit (0x08) or 16bit (0x10) or un-signed (0x40)
			Redimension/I image
		endif
		MatrixOP/O image = image - dark
	endif
	String wnote = note(image)
	Variable keV = NumberByKey("keV", wnote,"=")
	keV = NumberByKey("keV", wnote,"=")
	if (numtype(keV))
		DoAlert 0,"invalid keV in image '"+NameOfWave(image)+"'"
		return ""
	endif
	ImageStats/M=1/Q image								// check for empty images, skip them
	if (V_min==0 && V_max==0)
		wnote = ReplaceNumberByKey("V_min",wnote,V_min,"=")	// flag that this is empty
		wnote = ReplaceNumberByKey("V_max",wnote,V_max,"=")
		return wnote										// image is empty, do not waste time adding zeros
	endif
	Variable inorm = 	I0normalize ? I0normalizationFromNote(wnote) : 1	// normalization by intenstiy & exposure time

	MatrixOP/FREE Qvecs = Qvecs1keV * keV			// actual Q-vector for each pixel

	// for each point in image, sort into Q vector
	Variable Qx0=DimOffset(Qhist,0), Qy0=DimOffset(Qhist,1), Qz0=DimOffset(Qhist,2)
	Variable dQx=DimDelta(Qhist,0), dQy=DimDelta(Qhist,1), dQz=DimDelta(Qhist,2)
	Variable ip,jp, i,j,k
	for (jp=0;jp<Nj;jp+=1)								// for each pixel in the detector
		for (ip=0;ip<Ni;ip+=1)
			if (maskIn[ip][jp])
				i = round((Qvecs[ip][jp][0] - Qx0) / dQx)
				j = round((Qvecs[ip][jp][1] - Qy0) / dQy)
				k = round((Qvecs[ip][jp][2] - Qz0) / dQz)
				Qhist[i][j][k] += image[ip][jp] * inorm
				QhistNorm[i][j][k] += 1
			endif
		endfor
	endfor

	Variable seconds = stopMSTimer(-2)/1e6 - timer	// stop timer here
	if (printIt && seconds>8.0)
		printf "	processing image '%s' took %s\r",NameOfWave(image),Secs2Time(seconds,5,2)
	endif
	return wnote
End


// Make an array the same size as an image, but filled with Q^ to each pixel
// for a 2K x 2K image, this takes 1 minute.
Static Function/WAVE MakeQvecsArray(roi,d,wnote,recip,Elo,Ehi,[depth,mask,printIt])
	STRUCT imageROIstruct &roi						// region of detector of interest
	STRUCT detectorGeometry &d					// the input detector geometry
	String wnote
	Wave recip											// reciprocal lattice, used to find range of Q's
	Variable Elo, Ehi									// Energy Range (keV)
	Variable depth										// usually determined from image file
	Wave mask											// optional mask to limit the pixels that get processed (use pixel when mask true)
	Variable printIt
	depth = ParamIsDefault(depth) ? NaN : depth
	printIt = ParamIsDefault(printIt) ? NaN : printIt
	printIt = numtype(printIt) ? ItemsInList(GetRTStackInfo(0))<3 : !(!printIt)

	if (imageROIstructBad(roi) || !(d.used))
		print "ERROR -- MakeQvecsArray(), invalid ROI"
		DoAlert 0, "ERROR -- MakeQvecsArray(), invalid ROI"
		return $""
	endif
	Variable seconds = stopMSTimer(-2)/1e6			// start timer for initial processing
	Variable startx = roi.xLo, groupx = roi.binx
	Variable starty = roi.yLo, groupy = roi.biny
	Variable Ni = roi.Nx, Nj = roi.Ny					// number of BINNED pixels
	Variable NiNj = Ni*Nj

	Make/N=(Ni,Nj,3)/D/FREE Qvecs1keV=NaN				// Qvector at 1keV

	// for each pixel in image, compute sin(theta) and save it in Qvecs1keV[][][3]
	Variable useMat = 0										// only use recip when it is not an identity matrix
	if (WaveExists(recip))
		MatrixOP/FREE recipInv = Inv(recip)
		MatrixOP/FREE diff = sum(magsqr(recipInv - recip^t))
		useMat = diff[0]>1e-9								// only use recip when it is not an identity matrix
	endif
	Make/N=3/D/FREE qhat, qij, hkl

	Make/N=(NiNj,2)/I/FREE pxpy							// all the pixels in the image, these will be un-binned & zero based
	pxpy[][0] = mod(p,Ni)*groupx + (groupx-1)/2 + startx	// set the X-pixel values [0,Ni-1]
	pxpy[][1] = floor(p/Ni)*groupy + (groupy-1)/2 + starty	// set the Y-pixel values [0,Nj-1]
	Wave qvecs = pixel2qVEC(d,pxpy,depth=depth)	// qvecs in beam line coords, at this point qvecs are normalized
	WaveClear pxpy

	Make/N=3/D/FREE ki={0,0,1}
	Variable factor = 4*PI/hc
	MatrixOP/FREE qLens = -(qvecs x ki) * factor	// length of each q-vector, qvecs^ x ki == -sin(theta)
	qvecs *= qLens[p]											// re-scale qhats to real q-vectors
	Qvecs1keV = qvecs[p+q*Ni][r]							// Qvectors at 1keV, just multiply by E to get Q
	WaveClear qLens

	// find range of Qx,Qy,Qz
	if (useMat)													// if recip is identity, this saves 1 minute for large images
		MatrixOP/FREE/O hkl = (recipInv x (qvecs^t))^t	// this is really hkl/E
	else
		Duplicate/FREE qvecs, hkl
	endif
	WaveClear qvecs

	if (WaveExists(mask))
		Make/N=(NiNj)/B/U/FREE maskLocal				// un-wrap mask into a 1-D wave
		maskLocal = mask[mod(p,Ni)][floor(p/Ni)] ? 1 : 0		// set maskLocal to passed values, 1=OK, 0=BAD
		hkl = maskLocal[p] ? hkl[p][q] : NaN			// remove values that are masked off (set to NaN)
		WaveClear maskLocal
	endif

	MatrixOp/FREE QxLo0 = minVal(ReplaceNaNs(col(hkl,0),Inf))	// minVal() & maxVal don't support NaN's, but ±Inf is OK
	MatrixOp/FREE QyLo0 = minVal(ReplaceNaNs(col(hkl,1),Inf))
	MatrixOp/FREE QzLo0 = minVal(ReplaceNaNs(col(hkl,2),Inf))
	MatrixOp/FREE QxHi0 = maxVal(ReplaceNaNs(col(hkl,0),-Inf))
	MatrixOp/FREE QyHi0 = maxVal(ReplaceNaNs(col(hkl,1),-Inf))
	MatrixOp/FREE QzHi0 = maxVal(ReplaceNaNs(col(hkl,2),-Inf))
	Variable QxLo=QxLo0[0], QyLo=QyLo0[0], QzLo=QzLo0[0], QxHi=QxHi0[0], QyHi=QyHi0[0], QzHi=QzHi0[0]

	QxLo *= QxLo>0 ? Elo : Ehi	;	QxHi *= QxHi>0 ? Ehi : Elo
	QyLo *= QyLo>0 ? Elo : Ehi	;	QyHi *= QyHi>0 ? Ehi : Elo
	QzLo *= QzLo>0 ? Elo : Ehi	;	QzHi *= QzHi>0 ? Ehi : Elo
	wnote = ReplaceNumberByKey("QxLo",wnote,QxLo,"=")
	wnote = ReplaceNumberByKey("QyLo",wnote,QyLo,"=")
	wnote = ReplaceNumberByKey("QzLo",wnote,QzLo,"=")
	wnote = ReplaceNumberByKey("QxHi",wnote,QxHi,"=")
	wnote = ReplaceNumberByKey("QyHi",wnote,QyHi,"=")
	wnote = ReplaceNumberByKey("QzHi",wnote,QzHi,"=")
	wnote = ReplaceNumberByKey("startx",wnote,startx,"=")
	wnote = ReplaceNumberByKey("endx",wnote,roi.xHi,"=")
	wnote = ReplaceNumberByKey("groupx",wnote,groupx,"=")
	wnote = ReplaceNumberByKey("starty",wnote,starty,"=")
	wnote = ReplaceNumberByKey("endy",wnote,roi.yHi,"=")
	wnote = ReplaceNumberByKey("groupy",wnote,groupy,"=")
	Note/K Qvecs1keV, wnote

	seconds = stopMSTimer(-2)/1e6 - seconds
	if (printIt && seconds>5)
		printf "filling array of Q's at depth = %g,  took %s\r", depth, Secs2Time(seconds,5,1)
	endif
	return Qvecs1keV
End


// Conversions for Undulator A
Function gap2keV_A(mm)		// convert undulator gap to energy, This is the default, if you have one named gap2keV(mm) it will be used instead
	Variable mm					// undulator gap (mm)
	Make/D/FREE coef={-16.3368622327141,6.6359060831135,-1.06146366183031,0.0935193681319296,-0.00456055338213102,0.000126380924476218,-1.8912343101741e-06,1.19898530997386e-08}
	return poly(coef,mm)	// return energy (keV)
End
//
Function keV2gap_A(keV)	// convert undulator energy to gap, This is the default, if you have one named keV2gap(keV) it will be used instead
	Variable keV				// energy (keV)
	Make/D/FREE coef={-15.1766965984314,22.1471572343715,-8.44915267136723,1.92769824888557,-0.263034439599908,0.0213128570948232,-0.000946415059657785,1.78104843201671e-05}
	return poly(coef,keV)	// return gap (keV)
End


Static Function/T requestPathFileFmt(pathName)
	// returns a full path and fmt string as a list, e.g. "HD:Users:name:data:;EW1_%d.h5"
	// the user figures it out with some help from the Open
	// returns the list "path;fileFormat"
	String pathName
	PathInfo $pathName
	pathName = SelectString(V_flag,"",pathName)
	String message="pick any reconstructed image file in range,  using datafolder = "+GetDataFolder(0)
	Variable refNum
	initEnergyScans()
	SVAR ImageFileFilters=root:Packages:imageDisplay:ImageFileFilters
	Open/F=ImageFileFilters/D/M=message/P=$pathName/R refNum
	if (strlen(S_fileName)<1)
		return ""
	endif
	String fullPath = ParseFilePath(1,S_fileName,":",1,0)
	String fileRoot = ParseFilePath(0,S_fileName,":",1,0)

	Prompt fileRoot,"put in '%d' where desired"
	DoPrompt "Set filename format", fileRoot
	if (V_flag)
		return ""
	endif
	return fullPath+";"+fileRoot
End


Function TrimZerosOff3D(volIN)
	Wave volIN				// a 3D wave

	if (!WaveExists(volIN))
		return 1
	elseif (!WaveDims(volIN)==3)
		return 1
	endif
	Variable Nx=DimSize(volIN,0), Ny=DimSize(volIN,1), Nz=DimSize(volIN,2)
	Variable X0=DimOffset(volIN,0), Y0=DimOffset(volIN,1), Z0=DimOffset(volIN,2)
	Variable dX=DimDelta(volIN,0), dY=DimDelta(volIN,1), dZ=DimDelta(volIN,2)
	String xUnits=WaveUnits(volIN,0), yUnits=WaveUnits(volIN,1), zUnits=WaveUnits(volIN,2)

	Make/N=(Nx,Ny,Nz)/FREE/B/U vol
	if (WaveType(volIN) & 7)			// a floating point type
		vol = numtype(volIN) ? 0 : volIN!=0
	else
		vol = volIN==0 ? 0 : 1
	endif
	if (sum(vol)==0)
		return 0						// all zero, trim nothing
	endif

	Make/N=(Ny,Nz)/FREE xsurf
	Make/N=(Nx,Nz)/FREE ysurf
	Make/N=(Nx,Ny)/FREE zsurf

	Variable ixlo,ixhi,iylo,iyhi,izlo,izhi	// new integer ranges
	for (ixlo=0;ixlo<Nx;ixlo+=1)		// find range of x
		xsurf = vol[ixlo][p][q]
		if (sum(xsurf))
			break
		endif
	endfor
	for (ixhi=Nx-1;ixhi>=ixlo;ixhi-=1)
		xsurf = vol[ixlo][p][q]
		if (sum(xsurf))
			break
		endif
	endfor

	for (iylo=0;iylo<Ny;iylo+=1)		// find range of y
		ysurf = vol[p][iylo][q]
		if (sum(ysurf))
			break
		endif
	endfor
	for (iyhi=Ny-1;iyhi>=iylo;iyhi-=1)
		ysurf = vol[p][iylo][q]
		if (sum(ysurf))
			break
		endif
	endfor

	for (izlo=0;izlo<Nz;izlo+=1)		// find range of z
		zsurf = vol[p][q][izlo]
		if (sum(zsurf))
			break
		endif
	endfor
	for (izhi=Nz-1;izhi>=izlo;izhi-=1)
		zsurf = vol[p][q][izhi]
		if (sum(zsurf))
			break
		endif
	endfor

	if (ixlo<0 && ixhi>=Nx  &&  iylo<0 && iyhi>=Ny  &&  izlo<0 && izhi>=Nz)
		return 1
	elseif (ixlo==0 && ixhi==(Nx-1)  &&  iylo==0 && iyhi==(Ny-1)  &&  izlo==0 && izhi==(Nz-1))
		return 0						// nothing to trim off, done
	endif
	Nx = ixhi-ixlo+1
	Ny = iyhi-iylo+1
	Nz = izhi-izlo+1
	X0 += ixlo*dX
	Y0 += iylo*dY
	Z0 += izlo*dZ
	if (Nx<1 || Ny<1 || Nz<1)
		return 1
	endif
	Duplicate/FREE volIN, volTemp
	Redimension/N=(Nx,Ny,Nz) volIN
	SetScale/P x, X0, dX, xUnits, volIN
	SetScale/P y, Y0, dY, yUnits, volIN
	SetScale/P z, Z0, dZ, zUnits, volIN
	volIN = volTemp[p+ixlo][q+iylo][r+izlo]
	return 0
End


// get the range of file index numbers by looking at all of the files in a directory
// this only works for the last index on a file, not for depth sorted images
Static Function/T get_FilesIndexRange(pathName,nameFmt,[printIt])
	String pathName					// probably 'imagePath'
	String nameFmt						// name part of file = "EW5_%d.h5"
	Variable printIt
	printIt = ParamIsDefault(printIt) ? NaN : printIt
	printIt = numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : !(!printIt)

	String list = directory(pathName)

	Variable i1 = strsearch(nameFmt,"%d",0)
	Variable jlen=strlen(nameFmt[i1+2,Inf])
	if (i1<0)
		return ""
	endif

	String nameFmtMatch = ReplaceString("%d",nameFmt,"*")
	String name, rlist=""
	Variable m,i, N=ItemsInList(list)
	for (i=0;i<N;i+=1)
		name = StringFromList(i,list)
		if (StringMatch(name,nameFmtMatch))
			m = strlen(name)-jlen
			m = str2num(name[i1,m])
			if (numtype(m)==0 && m>=0 && m<2e9)
				rlist += num2istr(m)+";"
			endif
		endif
	endfor
	rlist = SortList(rlist,";",2)
	String range = compressRange(rlist,";")
	if (printIt)
		printf "range = [%s],    %d points\r",range,ItemsInRange(range)
	endif
	return range
End



Function FillMovieMoreGapScan(image,initStr)	// finds pixel intensities along a movie
	Wave image
	String initStr												// not used, output from FillMovieMoreGapScan_Init()

	Variable lowThresh=15, maxPixel=maxValueOfType(image)	// regular undulator gap scan, use difference between this image and previous image
	String wnote = note(image)
	Variable index = NumberByKey("movieIndex",wnote,"=")		// movie frame index
	if (index<1)												// the first time through frame, save it to start differences
		Duplicate/O image, FillMovieMoreGapScan_i0
	endif

	Wave i0 = FillMovieMoreGapScan_i0					// the previously saved image
	MatrixOP/FREE dimage = clip((image-i0) * (greater(maxPixel,i0) && greater(maxPixel,image) && greater(image,lowThresh)),0,maxPixel)

	i0 = image													// reset i0 to image for the next time through
	image = dimage												// set image to the difference, what we want in movie frame
	return 0
End
//
Function/T FillMovieMoreGapScan_Init(list)
	String list
	Wave i0 = FillMovieMoreGapScan_i0
	if (WaveExists(i0))
		KillWaves/Z i0											// this should never be here, delete in case previous movie stopped early
	endif
	return ""
End
//
Function/T FillMovieMoreGapScan_CleanUp(list)
	String list
	FillMovieMoreGapScan_Init(list)						// no longer need FillMovieMoreGapScan_i0, kill it.
End

// ========================================= End of 3D Q-space, No Depth =========================================
// ===============================================================================================================
// ===============================================================================================================
// ===============================================================================================================




// ===============================================================================================================
// ===============================================================================================================
// ===============================================================================================================
// =========================================== Start of E scans Panel ============================================

Function/T FillEscanParametersPanel(strStruct,hostWin,left,top)
	String strStruct									// optional passed value of xtal structure, this is used if passed
	String hostWin										// name of home window
	Variable left, top								// offsets from the left and top

	NewDataFolder/O root:Packages:micro				// ensure that the needed data folders exist
	NewDataFolder/O root:Packages:micro:Escan
	Variable small=22, big=30, start

	Variable/G root:Packages:micro:Escan:d0			// d-spacing of unstrained material (nm)
	String/G root:Packages:micro:Escan:fldr			// new folder name, it will be at the top level
	String/G root:Packages:micro:Escan:title		// a title to be used on plots

	SetWindow kwTopWin,userdata(EscanPanelName)=hostWin+"#EscanPanel"
//	NewPanel/K=1/W=(left,top,left+221,top+430)/HOST=$hostWin
	NewPanel/K=1/W=(left,top,left+221,top+482)/HOST=$hostWin
	ModifyPanel frameStyle=0, frameInset=0
	RenameWindow #,EscanPanel

	Button buttonNewEscan,pos={13,5},size={195,20},proc=EnergyScans#EscanButtonProc,title="Set Up New Energy scan..."
	Button buttonNewEscan,help={"Create a new data folder for an energy scan, and prompt user for title and d0"}
	Button buttonNewEscan,userdata(active)="0"

	start = 30
	SetVariable setvarTitle,pos={13,start},size={190,18},title="title",fSize=12,value= root:Packages:micro:Escan:title
	SetVariable setvarTitle,help={"title for the plots"},font="Lucida Grande",noedit=1,frame=0
	SetVariable setvard0,pos={45,start+20},size={128,18},title="d0 (nm)",noedit=1,frame=0
	SetVariable setvard0,fSize=12,format="%.7f",limits={0,inf,0},value= root:Packages:micro:Escan:d0
	SetVariable setvard0,help={"referene d-spacing (nm)"},font="Lucida Grande"
	SetVariable setvarfldr,pos={40,start+40},size={162,18},title="folder name",font="Lucida Grande",fSize=12,noedit=1,frame=0
	SetVariable setvarfldr,help={"name of the igor data folder to make"},value= root:Packages:micro:Escan:fldr

	start = 100
	Button buttonFillQ_Distn,pos={13,start},size={195,20},proc=EnergyScans#EscanButtonProc,title="Fill Q-Distribution"
	Button buttonFillQ_Distn,help={"read images containing an Energy scan and make Q-distribution(s)"}
	Button buttonRePlotQDistn,pos={53,start+1*small},size={155,20},proc=EnergyScans#EscanButtonProc,title="Plot Q-Distribution (re-fit too)",fSize=9
	Button buttonRePlotQDistn,help={"re-plot the Q-distribution at all positions or depths"}
	Button buttonRGBforQDistn,pos={53,start+2*small},size={155,20},proc=EnergyScans#EscanButtonProc,title="Recalc RGB for Q-Distribution",fSize=9
	Button buttonRGBforQDistn,help={"Recalculate RGB for the Q-distribution at all positions or depths"}
	Button buttonFitSingleQdistn,pos={53,start+3*small},size={155,20},proc=EnergyScans#EscanButtonProc,title="Display & fit Q at 1 depth",fSize=9
	Button buttonFitSingleQdistn,help={"you should not need this, this is usually done by changing the line on the Q-depth plot"}
	Button buttonlayoutQDistn,pos={53,start+4*small},size={155,20},proc=EnergyScans#EscanButtonProc,title="Make Layout of Q-Distribution",fSize=9
	Button buttonlayoutQDistn,help={"make a nice layout for printing showing the Q-distribution surface plot and the plot of a cut through it"}

	start = start+4*small + big
	Button buttonFillQ1image,pos={13,start},size={195,20},proc=EnergyScans#EscanButtonProc,title="Fill Q dist'n of 1 image"
	Button buttonFillQ1image,help={"make the Q distribution from a single loaded image"}
	Button buttonTransferd0,pos={13,start+big},size={195,20},disable=2,proc=EnergyScans#EscanButtonProc,title="d(hkl) --> d0     (xtal values)",fSize=10
	Button buttonTransferd0,help={"reset the d0 used to compute strains using d(nm) from 'Xtal' panel"}

	start = start+big+big
	Button buttonPseudoWhite,pos={13,start},size={195,20},proc=EnergyScans#EscanButtonProc,title="Make Pseudo White Image"
	Button buttonPseudoWhite,help={"read all the depth sorted images and make the pseudo white beam plot"}
	Button buttonRePlotEsumDepth,pos={53,start+small},size={155,20},disable=2,proc=EnergyScans#EscanButtonProc,title="Re-Plot Energy Sums vs Depth"
	Button buttonRePlotEsumDepth,help={"take the pseudo white data already read in and re-plot it"},fSize=9
	Button buttonMovieEsumsDepth,pos={53,start+2*small},size={155,20},disable=2,proc=EnergyScans#EscanButtonProc,title="Movie of Energy Sums vs Depth"
	Button buttonMovieEsumsDepth,help={"take the pseudo white data already read in and make the depth movie"},fSize=9

	start = start+2*small + big
	Button buttonMakeMask,pos={13,start},size={195,20},proc=EnergyScans#EscanButtonProc,title="Make a Mask"
	Button buttonMakeMask,help={"make a mask from an image (often the pseudo white image)"}
	Button buttonClearMask,pos={53,start+small},size={155,20},disable=2,proc=EnergyScans#EscanButtonProc,title="Clear the Mask"
	Button buttonClearMask,help={"clear the current mask (usually accessed by clicking in marquee)"},fSize=9

	start = start+small + big
	Button buttonFill3DQspace,pos={13,start},size={195,20},disable=0,proc=EnergyScans#EscanButtonProc,title="Fill a 3D Q-space",fSize=12
	Button buttonFill3DQspace,help={"Read In the Energy scan data & create a 3D Q-space"}
	Button buttonGizmoOf3DQspace,pos={53,start+small},size={155,20},disable=2,proc=EnergyScans#EscanButtonProc,title="make 3D Q-space Gizmo",fSize=10
	Button buttonGizmoOf3DQspace,help={"show a Gizmo from a previously loaded 3D Q-space"}

	start = start+small + big
	Button buttonLoadBadPixels,pos={16,start},size={150,20},proc=EnergyScans#EscanButtonProc,title="Load Bad Pixel File"
	Button buttonLoadBadPixels,help={"Load in the bad pixels for a detector"}
	Button buttonLoadBadPixels,fSize=10

	EnableDisableEscanControls(hostWin+"#EscanPanel")
	return "#EscanPanel"
End
//
Static Function EscanButtonProc(ctrlName) : ButtonControl
	String ctrlName

	String win = GetUserData("","","EscanPanelName")
	if (stringmatch(ctrlName,"buttonNewEscan"))
		Variable active = str2num(GetUserData(win,"buttonNewEscan","active"))
		NVAR d0=root:Packages:micro:Escan:d0
		if (!active)										// let user enter values and make new data folder
			SetVariable setvarTitle,win=$win,noedit=0,frame=1,disable=0
			SetVariable setvard0,win=$win,noedit=0,frame=1,disable=0
			SetVariable setvarfldr,win=$win,noedit=0,frame=1,disable=0
			Button buttonNewEscan,win=$win,title="Enter Values, then click here",fColor=(65535,65534,49151),userdata(active)="1"
			d0 = !(d0>0) ? NumVarOrDefault("root:Packages:Lattices:PanelValues:dspace_nm",d0) : d0
		else												// get values and create the new folder
			SVAR fldr=root:Packages:micro:Escan:fldr
			SVAR title=root:Packages:micro:Escan:title
			if (DataFolderExists("root:"+fldr))
				DoAlert 0, "the data folder 'root:"+fldr+"' alreay exists, try another name"
				return 1
			elseif (CheckName(fldr,11))
				DoAlert 0, "the data folder 'root:"+fldr+"' is not a legal folder name, try another name"
				return 1
			endif
			SetVariable setvarTitle,win=$win,noedit=1,frame=0,disable=2
			SetVariable setvard0,win=$win,noedit=1,frame=0,disable=2
			SetVariable setvarfldr,win=$win,noedit=1,frame=0,disable=2
			Button buttonNewEscan,win=$win,title="Set Up New Energy scan...",fColor=(0,0,0),userdata(active)="0"
			if (NewEnergyScan(fldr,title,d0))		// set up a folder for a new Energy Scan
				DoAlert 0,"New folder not made"
				fldr = ""
			endif
		endif
	elseif (stringmatch(ctrlName,"buttonFillQ_Distn"))
		printf "%sFill_Q_Positions(NaN,\"imagePath\",\"\",\"\",\"\",$\"\")\r",BULLET
		Fill_Q_Positions(NaN,"imagePath","","","",$"",printIt=1)
	elseif (stringmatch(ctrlName,"buttonFitSingleQdistn"))
//		printf "%sMakeGraph_Qhist($\"\")\r",BULLET
		MakeGraph_Qhist($"")
	elseif (stringmatch(ctrlName,"buttonRePlotQDistn") && strlen(WaveListClass("QdistAtPositions","*","MINCOLS:1")))
//	elseif (stringmatch(ctrlName,"buttonRePlotQDistn") && strlen(WaveListClass("QdistAtPositions,Qhistogram","*","MINCOLS:1")))
//	elseif (stringmatch(ctrlName,"buttonRePlotQDistn") && strlen(WaveListClass("QdistAtPositions,Qhistogram","*","")))
		MakeOtherQhistWaves($"")
	elseif (stringmatch(ctrlName,"buttonRGBforQDistn") && strlen(WaveListClass("Q_Positions,Qhistogram","*","")))
		MakeRGBforQdistribution($"",$"")
	elseif (stringmatch(ctrlName,"buttonlayoutQDistn") && strlen(WinList("Graph*_Q",";","WIN:1")))
		printf "%sMakeLayoutQ_Positions_hist(\"\")\r",BULLET
		MakeLayoutQ_Positions_hist("")
	elseif (stringmatch(ctrlName,"buttonPseudoWhite"))
		printf "%sMakePseudoWhiteImages(\"\",\"\")\r",BULLET
		MakePseudoWhiteImages("","")
	elseif (stringmatch(ctrlName,"buttonRePlotEsumDepth") && WaveExists(imageEsum))
		printf "%sMakeEsumPlot()\r",BULLET
		MakeEsumPlot()
	elseif (stringmatch(ctrlName,"buttonMovieEsumsDepth") && WaveExists(imageEsum))
		printf "%sMakeEsumMovie()\r",BULLET
		MakeEsumMovie()
	elseif (stringmatch(ctrlName,"buttonMakeMask") && strlen(WaveListClass("rawImage*;ImageSummed*","*","DIMS:2")))
		printf "%sMakeMask($\"\")\r",BULLET
		MakeMask($"")
	elseif (stringmatch(ctrlName,"buttonClearMask") && MaskToClearExists())
		ClearMask()
	elseif (stringmatch(ctrlName,"buttonFillQ1image"))
		printf "%sFill_Q_1image(NaN,$\"\")\r",BULLET
		Fill_Q_1image(NaN,$"")
	elseif (stringmatch(ctrlName,"buttonTransferd0"))
		update_d0()
	elseif (stringmatch(ctrlName,"buttonFill3DQspace"))
		Fill1_3DQspace(0,"imagePath","","",printIt=1)
	elseif (stringmatch(ctrlName,"buttonGizmoOf3DQspace"))
		MakeGizmoQspace3d($"")
	elseif (stringmatch(ctrlName,"buttonLoadBadPixels"))
		LoadDefaultBadPixelImage($"","")
	endif

	EnableDisableEscanControls(GetUserData("microPanel","","EscanPanelName"))
End
//
Static Function EnableDisableEscanControls(win)			// here to enable/disable
	String win												// window (or sub window) to use
	Variable d

	Button buttonNewEscan,win=$win,disable=0			// always OK to make a new data folder
	Button buttonFillQ_Distn,win=$win,disable=0			// always OK to load a Q-Depth surface
	Button buttonPseudoWhite,win=$win,disable=0			// always OK to make a Pseudo White Image

	d = strlen(WaveListClass("Qhistogram","*","DIMS:1")) ? 0 : 2
	Button buttonFitSingleQdistn,win=$win,disable= d

	d = WaveExists(imageEsum) ? 0 : 2						// only if imageEsum exists
	Button buttonRePlotEsumDepth,win=$win,disable= d
	Button buttonMovieEsumsDepth,win=$win,disable= d

	d = strlen(WaveListClass("rawImage*;ImageSummed*","*","DIMS:2")) ? 0 : 2
	Button buttonMakeMask,win=$win,disable=d			// OK to make mask if an imge is around
	d = MaskToClearExists() ? 0 : 2
	Button buttonClearMask,win=$win,disable= d

	d = strlen(WaveListClass("QdistAtPositions","*","MINCOLS:1")) ? 0 : 2
	Button buttonRePlotQDistn,win=$win,disable= d

	d = strlen(WaveListClass("Q_Positions","*","")) ? 0 : 2
	Button buttonRGBforQDistn,win=$win,disable= d

	d = strlen(WinList("Graph*_Q",";","WIN:1")) ? 0 : 2
	Button buttonlayoutQDistn,win=$win,disable= d

	d = strlen(WaveListClass("rawImage*","*","DIMS:2")) ? 0 : 2
	Button buttonFillQ1image,win=$win,disable= d

	d = exists("d0")==2 ? 0 : 2								// only if a local d0 exists
	Button buttonTransferd0,win=$win,disable= d

	d = ItemsInList(WaveListClass("GizmoXYZ;Qspace3D","*","DIMS:3"))>0 ? 0 : 2
	Button buttonGizmoOf3DQspace,win=$win,disable= d

	if (!str2num(GetUserData(win,"buttonNewEscan","active")))	// new Escan NOT active, set to current folder values
		NVAR d0=root:Packages:micro:Escan:d0
		SVAR title=root:Packages:micro:Escan:title, fldr=root:Packages:micro:Escan:fldr
		String fldrLocal = GetDataFolder(0)
		if (stringmatch(GetDataFolder(1),"root:"+fldrLocal+":"))
			NVAR d0Local = :d0
			SVAR titleLocal = :title
			if (SVAR_Exists(titleLocal) || NVAR_Exists(d0Local))
				d0 = d0Local
				title = titleLocal
			endif
		endif
	endif
End

//
Static Function MaskToClearExists()
	String imageName = StringFromList(0,ImageNameList("",";"))
	if (strlen(imageName)<1)
		return 0
	endif
	Wave image = ImageNameToWaveRef("",imageName)
	if (!WaveExists(image))
		return 0
	endif
	String class = StringByKey("waveClass", note(image),"=")
	if (!stringmatch(class,"image&Mask"))	// is not an editable mask, disable
		return 0
	endif
	return 1
End
//
//
Static Function update_d0()	// take current calculation for d(hkl) from Xtal panel and update local value of d0
	Variable dhkl = NumVarOrDefault("root:Packages:Lattices:PanelValues:dspace_nm",NaN)
	if (numtype(dhkl)==0 && dhkl>0)
		NVAR d0default=root:Packages:micro:Escan:d0
		if (NVAR_Exists(d0default))					// d-spacing of unstrained material (nm)
			d0default = dhkl
		endif
		NVAR d0Local = :d0
		if (NVAR_Exists(d0Local))						// d-spacing in current folder
			Variable h=NumVarOrDefault("root:Packages:Lattices:PanelValues:h",NaN)
			Variable k=NumVarOrDefault("root:Packages:Lattices:PanelValues:k",NaN)
			Variable l=NumVarOrDefault("root:Packages:Lattices:PanelValues:l",NaN)
			printf "in folder  '%s', using hkl=(%s),  set d0 to  %8g nm\r",GetDataFolder(0),hkl2str(h,k,l),dhkl
			d0Local = dhkl
		endif
	endif
End

// ============================================ End of E scans Panel =============================================
// ===============================================================================================================
// ===============================================================================================================
// ===============================================================================================================




// ===============================================================================================================
// ===============================================================================================================
// ================================================ Start of Init ================================================

Function initEnergyScans()
	if (exists("root:Packages:micro:DEFAULT_I0_GAIN")!=2)
		Variable/G root:Packages:micro:DEFAULT_I0_GAIN=DEFAULT_I0_GAIN
	endif
	if (exists("root:Packages:micro:DEFAULT_I0_SCALING")!=2)
		Variable/G root:Packages:micro:DEFAULT_I0_SCALING=DEFAULT_I0_SCALING
	endif
	initImageDisplayScaling()
End

// ================================================= End of Init =================================================
// ===============================================================================================================
// ===============================================================================================================
