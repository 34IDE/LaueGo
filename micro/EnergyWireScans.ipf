#pragma rtGlobals=1		// Use modern global access method.
#pragma ModuleName=EnergyWireScans
#pragma version = 0.9948

#include "WinView", version>=1.82
#include "microGeometry", version>=2.56
#include "ImageDisplayScaling", version>=1.44
#include "ExamineImageSet", version>=1.00

Static Constant hc = 1.239841857				// hc (keV-nm)
Static Constant secPerPixelFixed = 18.8e-6		// the fixed time is takes to process one pixel (sec) after any distortion
Static Constant secPerPixelDistort = 15.24e-3	// time is takes to process distortion for one pixel (sec)  (measured with a 2GHz clock)


//	**	To get the ful chip unbinned pixel position from ROI data use the following:
//	**	fullchip x  =  groupx*x + (groupx-1)/2 + (startx - 1)
//	**	fullchip y  =  groupy*y + (groupy-1)/2 + (starty - 1)
//		py =  j*groupy + (groupy-1)/2 + (starty-1)
//		px = i*groupx + (groupx-1)/2 + (startx-1)
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
//	s = ∆z/∆y of ray from sample to pixel on CCD
//
//	z = (  zs - s*Fo/cos(angle) - R*sqrt(1-s*s)  ) / (1-s*tan(angle) )
//	H = [z/cos(angle) - Fo*tan(angle)]] = [z - Fo*sin(angle)] / cos(angle)
//
//	for a Gaussian,   fwhm = 2*K3*sqrt(ln(2))
//	for a Lorentzian,   fwhm = 2*sqrt(K3)


//Menu "GraphMarquee"
//	"-"
//	"Add Marquee to Mask",/Q, AddMarqueeToMask()
//	help={"add the area of this marquee to the mask, (this will include the area of the marquee in the analysis"}
//	"Subtract Marquee From Mask",/Q,SubtractMarqueeFromMask()
//	help={"remove the area of this marquee from the mask, (all pixels within this area will not be used in the analysis"}
//	"Clear the Mask",ClearMask()
//	help={"reset the whole mask to empty, no points will be analyzed.  Use this to clear a mask so you can start addding points again"}
//End
Menu "GraphMarquee",dynamic
	"-"
	MarqueeMaskMenuItem("Add Marquee to Mask"),/Q, AddMarqueeToMask()
	help={"add the area of this marquee to the mask, (this will include the area of the marquee in the analysis"}
	MarqueeMaskMenuItem("Subtract Marquee From Mask"),/Q,SubtractMarqueeFromMask()
	help={"remove the area of this marquee from the mask, (all pixels within this area will not be used in the analysis"}
	MarqueeMaskMenuItem("Clear the Mask"),ClearMask()
	help={"reset the whole mask to empty, no points will be analyzed.  Use this to clear a mask so you can start addding points again"}
End
//
Function/S MarqueeMaskMenuItem(item)
	String item

	String imageName = StringFromList(0,ImageNameList("",";"))
	if (strlen(imageName)<1)
		return "("+item
	endif
	Wave image = ImageNameToWaveRef("",imageName)
	if (!WaveExists(image))
		return "("+item
	endif
	String class = StringByKey("waveClass", note(image),"=")
	if (!stringmatch(class,"image&Mask"))	// is not an editable mas, disable
		item = "("+item
		// item = "("
	endif
	return item
End


Menu "micro", dynamic
	Submenu "Energy-Wire Scans"
		"New Energy-Wire scan...",NewEnergyWireScan("","",NaN)
		help={"Create a new data folder for an energy-wire scan, and prompt user for title and d0"}
		"run:  Fill E vs Depth",Fill_EvsDepth("imagePath","")
		help={"read the depth sorted images, and make a 2-d plot of energy vs depth, this command is not used much"}
		"run:  Fill Q vs Depth",Fill_QvsDepth(NaN,"imagePath","","","",$"")
		help={"read the depth sorted images, and make a 2-d plot of Q vs depth"}
		     MenuItemIfWaveExists("    compute & fit Q at 1 depth","Q_Depth"),fitQat1depth(NaN,":")
		help={"you should not need this, this is usually done by changing the line on the Q-depth plot"}
		"Make Pseudo White Images",MakePseudoWhiteImages("","")
		      help={"read all the depth sorted images and make the pseudo white beam plot"}
		     MenuItemIfWaveExists("    Make movie of Energy Sums vs depth","imageEsum"),MakeEsumMovie()
		      help={"take the pseudo white data already read in and make the depth movie"}
		     MenuItemIfWaveExists("    Re-Plot Energy Sums vs depth","imageEsum"),MakeEsumPlot()
		      help={"take the pseudo white data already read in and re-plot it"}
		"run:  Fill Q at Various Positions",Fill_Q_Positions(NaN,"imagePath","","",$"")
	      help={"read all energy scans (with no wire) and at many positions, and make a 2-d plot of Q at various positions"}
		"-"
		MenuItemIfWaveClassExists("Make Mask","speImage*;ImageSummed*","DIMS:2"),MakeMask($"")
		      help={"make a mask from an image (often the pseudo white image)"}
		//"    Add Marquee to Mask",/Q, AddMarqueeToMask()
		//      help={"add area of marquee to mask (usually accessed by clicking in marquee)"}
		//"    Subtract Marquee From Mask",/Q,SubtractMarqueeFromMask()
		//      help={"subtract area of marquee from mask (usually accessed by clicking in marquee)"}
		MarqueeMaskMenuItem("    Clear the Mask"),ClearMask()
		      help={"clear the current mask (usually accessed by clicking in marquee)"}
		"-"
		MenuItemIfWaveExists("Graph Intensity on [keV x Depth]","keV_Depth"),MakeGraph_keV_Depth()
		help={"re-plot the energy-depth info read in using Fill_EvsDepth()"}
		MenuItemIfWaveExists("Graph Intensity on [Q x Depth]","Q_Depth"),MakeGraph_Q_Depth()
		help={"re-plot the Q-depth info read in using Fill_QvsDepth()"}
//		MenuItemIfWaveClassExists("re-plot x-y positions [to locate Q distns]","Q_Positions*","DIMS:2"),MakeGraph_Q_Positions(Q_Positions_RGB)
		MenuItemIfWaveExists("re-plot x-y positions [to locate Q distns]","Q_Positions_RGB"),MakeGraph_Q_Positions(Q_Positions_RGB)
		     MenuItemIfWaveExists("    Graph 1 Q historam","Qhist"),MakeGraph_Qhist($"")
		     help={"re-plot the single slice from the Q-depth surface"}
		MenuItemIfGraphNamesExist("Layout [Q x Depth] & Q histogram","*_Q_Depth;*_Qhist"),MakeLayoutQ_depth_hist("")
		"Re-scale Q wave to strain",reScaleQwaveToStrain($"",NaN)
		help={"make a nice layout for printing showing the Q-depth surface plot and the plot of a cut through it"}
	End
//	help={"the commands associated with the Energy-Wire Scans package"}

	Submenu "Set Igor Data Folder"
		DataFolderListForMenus(),/Q, SetDataFolderFromMenu1()		// set to an-other data folder
	End
	help={"used to move between top level folders, also writes a note to history"}
End
//
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
	printf "\r\r•• changed data folder from  '%s'  -->  '%s'\r",from,to
End

Function/S MenuItemIfWaveExists(item,wName)
	String item
	String wName
	return SelectString(WaveExists($wName) ,"(","")+item
End
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



Function/S fitQat1depth(row,dataFolder)
	Variable row					// row index in Q_Depth[row][], the first index is depth
	String dataFolder				// data folder with the data
	if (!(row>=0))
		row = 8
		Prompt row, "row index for constant depth cut"
		DoPrompt "row index",row
		if (V_flag)
			return ""
		endif
	endif
	dataFolder = SelectString(strlen(dataFolder)<1,dataFolder,":")	// default to current folder
	String fldrSav= GetDataFolder(1)
	SetDataFolder $dataFolder

	Wave Q_Depth=Q_Depth, Qhist=Qhist
	Qhist = Q_Depth[row][p]							// select the row
	Variable depth = row*DimDelta(Q_Depth,0) + DimOffset(Q_Depth,0)
	Variable V_FitOptions=4, V_FitError=0				//suppress window & ignore errors
	Variable lo,hi, fwhm=NaN, Qc=NaN, Qc_precision=6,fwhm_precision=2
	if(!getPercentOfPeak(Qhist,0.5,lo,hi))				// find the 50% points
		CurveFit/Q lor Qhist[lo,hi]/D 
		if (!V_FitError)
			fwhm = 2*sqrt(K3)						// for Lorentzian
			Qc = K2
			Wave W_coef=W_coef, W_sigma=W_sigma
			Qc_precision = ceil(-log(W_sigma[2]/W_coef[2]))+1
			fwhm_precision = ceil(-log(sqrt(W_sigma[3]/W_coef[3])))+1
			fwhm_precision = numtype(fwhm_precision) ? 2 : limit(fwhm_precision,1,15)
		endif
	endif
	String list = ReplaceNumberByKey("Qcenter","",Qc,"=")
	list = ReplaceNumberByKey("Qfwhm",list,fwhm,"=")

	String wnote = note(Qhist)
	wnote = ReplaceNumberByKey("Qhist_row",wnote,row,"=")
	wnote = ReplaceNumberByKey("Qhist_depth",wnote,depth,"=")
	wnote = ReplaceNumberByKey("Qcenter",wnote,Qc,"=")
	wnote = ReplaceNumberByKey("Qfwhm",wnote,fwhm,"=")
	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
		printf "at peak, Q = %.5g(1/nm),   fwhm = %.2g(1/nm),   fwhm/Q = %.1e,      chisq=%.3g\r",Qc,fwhm, fwhm/Qc,V_chisq
	endif
	Variable d0 = NumVarOrDefault("d0",NaN)			// unstrained d-spacing (nm)
	if (numtype(d0)==0)
		Variable Q0 = 2*PI/d0							// Q of unstrained material (1/nm)
		Variable strain = -(Qc-Q0)/Qc					// d is inverse of strain, so a negative
		list = ReplaceNumberByKey("d0",list,d0,"=")
		list = ReplaceNumberByKey("Q0",list,Q0,"=")
		list = ReplaceNumberByKey("strain",list,strain,"=")
		if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
			printf "   d0 = %g (nm),   Q0 = %g (1/nm),   strain = %.1e\r",d0,Q0,strain
		endif
		wnote = ReplaceNumberByKey("∆Q",wnote,Qc-Q0,"=")
		wnote = ReplaceNumberByKey("strain",wnote,strain,"=")
	endif
	Note/K Qhist, wnote

	// if Qhist is on a graph, update textbox and line, so find a plot (if one exists, containing Qhist)
	String qwName = GetWavesDataFolder(Qhist,2)
	String win, wList = WinList("*",";","WIN:1")
	Variable j, i, N=ItemsInList(wList)
	Variable foundIt=0
	for (i=0;i<N && !foundIt;i+=1)						// search for graph with Qist on it
		win = StringFromList(i,wList)					// name of current window to test
		j = 0
		do
			Wave ww=WaveRefIndexed(win,j,1)
			if (!WaveExists(ww))
				break
			endif
			if (stringmatch(GetWavesDataFolder(ww,2),qwName))
				foundIt = 1
				break
			endif
			j += 1
		while(1)
	endfor
	if (foundIt)											// found the graph with Qhist
		if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
			DoWindow/F $win
		endif
		if (WhichListItem("textQ", AnnotationList(win),";")>=0)	// textQ exists, update it
			String str
			String fmt
			sprintf fmt "Q\\Bc\\M = %%.%dg (nm\\S-1\\M)\rFWHM = %%.%dg",Qc_precision, fwhm_precision
			sprintf str, fmt,Qc,fwhm
			TextBox/W=$win/C/N=textQ str
			if (numtype(strain)==0)
				sprintf str,"strain = %.1e",strain
				AppendText/W=$win/N=textQ str
			endif
			sprintf str,"depth = %g µm",depth
			AppendText/W=$win/N=textQ str
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
//
Static Function getPercentOfPeak(pkWave,fraction,lo,hi)		// return index in to pkWave of the ±50% of peak points, or use cursors A & B
	Wave pkWave					// wave with peak
	Variable fraction				// fraction of peak to go on each side
	Variable &lo, &hi
	lo = NaN  ;  hi = NaN
	if (!WaveExists(pkWave) || !(0<fraction && fraction<1))
		return 1
	endif
	Variable N=numpnts(pkWave)-1

	String win = StringFromList(0,FindGraphsWithWave(pkWave))	// find top most graph containing pkWave
	lo = NumberByKey("POINT",CsrInfo(A,win))
	hi = NumberByKey("POINT",CsrInfo(B,win))
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
	if (!(1<imax && imax<N-1))									// too close to the edges
		return 1
	endif

	Variable half = (V_max+max(V_min,0))*fraction				// half way point
	if (!(half<V_max))
		return 1
	endif

	for (lo=imax;lo>0;lo-=1)
		if (pkWave[lo]<half)
			break
		endif
	endfor
	for (hi=imax;hi<N;hi+=1)
		if (pkWave[hi]<half)
			break
		endif
	endfor
	return 0
End
//Static Function getPercentOfPeak(pkWave,fraction,lo,hi)		// return index in to pkWave of the ±50% of peak points
//	Wave pkWave					// wave with peak
//	Variable fraction				// fraction of peak to go on each side
//	Variable &lo, &hi
//	lo = NaN  ;  hi = NaN
//	if (!WaveExists(pkWave) || !(0<fraction && fraction<1))
//		return 1
//	endif
//	Variable N=numpnts(pkWave)-1
//	WaveStats/M=1/Q pkWave
//
//	Variable imax,i
//	imax = (V_maxloc-DimOffset(pkWave,0))/DimDelta(pkWave,0)// index of the max
//	if (!(1<imax && imax<N-1))									// too close to the edges
//		return 1
//	endif
//
//	Variable half = (V_max+max(V_min,0))*fraction				// half way point
//	if (!(half<V_max))
//		return 1
//	endif
//
//	for (lo=imax;lo>0;lo-=1)
//		if (pkWave[lo]<half)
//			break
//		endif
//	endfor
//	for (hi=imax;hi<N;hi+=1)
//		if (pkWave[hi]<half)
//			break
//		endif
//	endfor
//	return 0
//End
//Function fitQat1depth(row)
//	Variable row					// row index in Q_Depth[row][], the first index is depth
//	if (!(row>=0))
//		row = 8
//		Prompt row, "row index for constant depth cut"
//		DoPrompt "row index",row
//		if (V_flag)
//			return 1
//		endif
//	endif
//	Wave Q_Depth=Q_Depth, Qhist=Qhist
//	Qhist = Q_Depth[row][p]							// select the row
//	CurveFit/Q lor Qhist/D 
//	Variable fwhm = 2*sqrt(K3)						// for Lorentzian
//	printf "at peak, Q = %.5g(1/nm),   fwhm = %.2g(1/nm),   fwhm/Q = %.1e,      chisq=%.3g\r",K2,fwhm, fwhm/K2,V_chisq
//End






//	fileRoot="Macintosh HD:Users:tischler:data:Sector 34:July 4, 2006:EW5:recon:EW5_1_"
//	FillQhist("Macintosh HD:Users:tischler:data:Sector 34:July 4, 2006:EW5:recon:EW5_1_","0-2")
//	Fill_QvsDepth("Macintosh HD:Users:tischler:data:Sector 34:July 4, 2006:EW5:recon:EW5_","1,22,43,64,85,106,127,148,169,190,211,232,253,274,295,316,337,358,379","0-60")
//	Fill_QvsDepth("Macintosh HD:Users:tischler:data:Sector 34:July 4, 2006:EW5:recon:EW5_","211,232","15-22")





// process an E-W scan that has been already depth sorted
Function Fill_QvsDepth(d0,pathName,namePart,range1,range2,mask,[depthSi])
	Variable d0		// d-spacing of the strain=0 material (nm)
	String pathName	// either name of path to images, or the full expliiciit path, i.e. "Macintosh HD:Users:tischler:data:cal:recon:"
	String namePart	// the first part of file name, something like  "EW5_"
	String range1		// range of first number after file root (designates energy)
	String range2		// range of second numbers after range1 (designates depth)
	Wave mask		// optional mask to limit the pixels that get processed (use pixel when mask true)
	Variable depthSi	// depth measured from Si, needed if images are not the result of a reconstruction (micron)
	depthSi = ParamIsDefault(depthSi) ? NaN : depthSi

	Variable printIt=0
	String str
	PathInfo $pathName
	if (!V_flag || strlen(namePart)<1)					// path does not exist or no namePart, ask user
		String pathPart
		str = ExamineImageSet#requestFileRoot(pathName,2)
		pathPart = StringFromList(0,str)
		namePart = StringFromList(1,str)
		if (strlen(pathPart)<1 || strlen(namePart)<1)
			return 1										// invalid inputs
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
		return 1											// invalid inputs
	endif
	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
		printf "using data from files starting with '%s'\r",S_path+namePart
	endif

	String list = WaveListClass("imageMask","*","DIMS:2,BYTE:1")	// setup to use optional mask
	String maskName = ""
	if (WaveExists(mask))
		maskName = GetWavesDataFolder(mask,2)
	elseif (strlen(list))
		maskName = ""
		Prompt maskName, "mask to use with image",popup,"_none_;"+list
		DoPrompt "mask wave",maskName
		if (V_flag)
			return 1
		endif
		maskName = SelectString(stringmatch(maskName,"_none_"),maskName,"")
		Wave mask = $maskName							// do not check if wave exists, that a valid option
		printIt = 1
	endif
	maskName = SelectString(WaveExists(mask),"$\"\"",maskName)

	Display/W=(49,105,687,231)/K=1				// put up temp window to view staus of processing
	String gName = S_name
	ModifyGraph gfSize=18
	TextBox/N=text0/F=0/S=3/A=LC/X=2.51/Y=4.76 "Determining the range to scan"  ;  DoUpdate

	if (strlen(range1)<=0 || strlen(range2)<=0)		// if either range is empty, get the full range
		str = ExamineImageSet#get_ranges(pathName,namePart)
		if (strlen(range1)<=0)
			range1=StringFromList(0,str)
			if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
				printf "setting energy range to '%s'\r",range1[0,370]	// cannot print lines longer than 400 chars
			endif
		endif
		if (strlen(range2)<=0)
			range2=StringFromList(1,str)
			if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
				printf "setting depth range to '%s'\r",range2[0,370]	// cannot print lines longer than 400 chars
			endif
		endif
		printIt = 1
		if (WaveExists(mask))							// prompt for user input on ranges is there is a mask
			Prompt range1,"range of wire scans (probably energy)"
			Prompt range2,"range of depths, (index numbers not micron)"
			DoPrompt "rnages",range1,range2
			if (V_flag)
				return 1
			endif
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
	if (!((d0>0)))										// invalid d0, check for a local value
		d0 = NumVarOrDefault("d0",NaN)
		printIt = 1
	endif
	if (printIt)
		if (numtype(depthSi))
			sprintf str,"Fill_QvsDepth(%g,\"%s\",\"%s\",\"%s\",\"%s\",%s)",d0,pathName,namePart,range1,range2,maskName
		else
			sprintf str,"Fill_QvsDepth(%g,\"%s\",\"%s\",\"%s\",\"%s\",%s,depthSi=%g)",d0,pathName,namePart,range1,range2,maskName,depthSi
		endif
		print str[0,390]
	endif
	if (WaveExists(mask))
		if (sum(mask)==0)
			DoAlert 0, "You picked a mask that is all zero, stopping"
			DoWindow/K $gName							// done with status window
			return 1
		endif
	endif

	Variable timer0 = startMSTimer						// start timer for entire prooess
	Variable timer = startMSTimer						// start timer for initial processing
	STRUCT microGeometry geo
	FillGeometryStructDefault(geo)
	Variable useDistortion = NumVarOrDefault("root:Packages:geometry:useDistortion",USE_DISTORTION_DEFAULT)
 	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
		printf "processing with distotion "+SelectString(useDistortion,"OFF","on")
		if (WaveExists(mask))
			printf "   and using  '%s'  for the mask\r",NameOfWave(mask)
		else
			printf "   and no mask\r"
		endif
	endif

	TextBox/C/N=text0/W=$gName  "setting up"  ;  DoUpdate
	String name
	name = ExamineImageSet#fileNameFromRangeNumbers(fileRoot,str2num(range1),str2num(range2))
	Wave image = $(LoadWinViewFile(name))			// load first image in both ranges
	if (!WaveExists(image))
		printf "could not load very first image named '%s'\r",name
		DoAlert 0,"could not load very first image"
		DoWindow/K $gName								// done with status window
		timer=stopMSTimer(timer)  ;  timer0=stopMSTimer(timer0)
		return 1
	endif
	Variable Ni,Nj, Npixels = numpnts(image)			// number of pixels in one image
	Ni = DimSize(image,0)
	Nj = DimSize(image,1)
	String wnote = note(image)
	KillWaves/Z image
	Variable ddLocal=NumberByKey("CCDy", wnote,"=") + geo.ddOffset
	geo.dd = !(ddLocal>0) ? geo.dd : ddLocal					// use dd from the CCDy if available
	Variable keV, yc, startx, endx, starty, endy, groupx, groupy
	keV = NumberByKey("keV", wnote,"=")				// energy for this image
	yc = NumberByKey("yc", wnote,"=")				// yc for this image
	if (numtype(yc) && numtype(depthSi)!=0)		// ask for depth
		Prompt depthSi,"depth of sample (measured from Si in µm)"
		DoPrompt "depth", depthSi
		if (V_flag)
			DoWindow/K $gName								// done with status window
			timer=stopMSTimer(timer)  ;  timer0=stopMSTimer(timer0)
			return 1
		endif
		printf "		using depthSi=%g\r",depthSi
	endif
	yc = numtype(yc) ? geo.ycent +depthSi*(geo.NyCCD/geo.dpsy/1000) : yc	// for when there is no energy scan
	yc = numtype(yc)&&(strlen(range2)==0) ? geo.ycent : yc	// if no depthSi is given by user, then just use value in geo
	geo.ycent = yc											// geo is only local, and never gets saved.  Use this yc if all else fails
	startx = NumberByKey("startx", wnote,"=");	endx = NumberByKey("endx", wnote,"=");	groupx = NumberByKey("groupx", wnote,"=")
	starty = NumberByKey("starty", wnote,"=");	endy = NumberByKey("endy", wnote,"=");	groupy = NumberByKey("groupy", wnote,"=")
	if (numtype(startx+endx+starty+endy+groupx+groupy))
		DoAlert 0,"could not get ROI from wave note of image '"+name+"'"
		DoWindow/K $gName								// done with status window
		timer=stopMSTimer(timer)  ;  timer0=stopMSTimer(timer0)
		return 1
	elseif (numtype(yc+keV))
		DoAlert 0,"invalid yc or keV in image '"+name+"'"
		DoWindow/K $gName								// done with status window
		timer=stopMSTimer(timer)  ;  timer0=stopMSTimer(timer0)
		return 1
	endif
	String strStruct										// string to hold geo for an easy copying
	StructPut/S/B=0 geo, strStruct						// this is done to copy geo into geoLocal
	STRUCT microGeometry geoLocal						// the local geo with yc of this image
	StructGet/S/B=0 geoLocal, strStruct				// now we can change yc in geoLocal without making trouble above

	Variable NkeV=ItemsInRange(range1), Ndepth=max(1,ItemsInRange(range2))
	Variable N = Ndepth*NkeV							// number of images to be processed
	Variable Q0 = 2*PI/d0								// Q of unstrained material (1/nm)

	Make/N=(N)/O/D depth_FillQvsDepth, keV_FillQvsDepth
	depth_FillQvsDepth = NaN
	keV_FillQvsDepth = NaN
	Make/N=(N)/O/U/I m1_FillQvsDepth,m2_FillQvsDepth
	m1_FillQvsDepth = 0
	m2_FillQvsDepth = 0
	Variable seconds
	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
		seconds = Npixels*N*secPerPixelFixed
		if (useDistortion)
			seconds += Npixels*secPerPixelDistort
		endif
		seconds *= numtype(cpuFrequency()) ? 1 : 2e9/cpuFrequency()
		printf "about to examine %d images, which should take %s,  and end at %s\r",N,Secs2Time(seconds,5,1),Secs2Time(DateTime+seconds,1)
	endif
	// for all the N files (go over both range1 & range2), store the energies, depths, and indicies
	TextBox/C/N=text0/W=$gName  "examining "+num2istr(N)+" image headers\r(to determine what to do)"  ;  DoUpdate
	Variable m1,m2, i, j
	for (m2=str2num(range2), i=0; numtype(m2)==0 || i==0; m2=NextInRange(range2,m2))// loop over range2, depths
		for (m1=str2num(range1); numtype(m1)==0; m1=NextInRange(range1,m1), i+=1)	// loop over range1, energies
			name = ExamineImageSet#fileNameFromRangeNumbers(fileRoot,m1,m2)
			wnote=WinViewReadHeader(name)			// wave note to add to file read in
			if (!sameROI(wnote,startx, endx, starty, endy, groupx, groupy))	// skip bad ROI's
				continue
			endif
			depth_FillQvsDepth[i] = NumberByKey("depthSi", wnote,"=")
			keV_FillQvsDepth[i] = NumberByKey("keV", wnote,"=")
			m1_FillQvsDepth[i] = m1
			m2_FillQvsDepth[i] = m2
		endfor
	endfor
	seconds = stopMSTimer(timer)/1e6
	if (ItemsInList(GetRTStackInfo(0))<2 && seconds>0.2)
		printf "first pass to examine headers took %s\r",Secs2Time(seconds,5,3)
	endif

	timer = startMSTimer								// start a timer on bulk of processing
	WaveStats/M=1/Q keV_FillQvsDepth
	Variable ikeVlo=V_minloc, ikeVhi=V_maxloc		// save location of min and max energies
	if (V_numNans)
		DoAlert 0, "There were "+num2istr(V_numNans)+" bad images found, they will be skipped"
	endif
	WaveStats/M=1/Q depth_FillQvsDepth
	Variable depthMin=V_min, depthMax=V_max		// depth range
	depthMin = numtype(depthMin) ? 0 : depthMin
	depthMax = numtype(depthMax) ? 0 : depthMax
	Variable dDepth = (depthMax-depthMin)/(max(Ndepth,2)-1)

	Variable theta, Qmin, Qmax, dQ, NQ					// get range of Q
	name = ExamineImageSet#fileNameFromRangeNumbers(fileRoot,m1_FillQvsDepth[ikeVlo],m2_FillQvsDepth[ikeVlo])
	Wave image = $(LoadWinViewFile(name))			// load image
	if (!WaveExists(image))
		printf "could not load image named '%s'\r",name
		DoAlert 0,"could not load image for finding theta range"
		DoWindow/K $gName								// done with status window
		timer=stopMSTimer(timer)  ;  timer0=stopMSTimer(timer0)
		return 1
	endif
	geoLocal.ycent = NumberByKey("yc", note(image),"=")	// change yc in local copy of geo
	geoLocal.ycent = numtype(geoLocal.ycent) ? geo.ycent : geoLocal.ycent
	theta = real(thetaRange(geoLocal,image))			// min theta on this image
	KillWaves/Z image
	keV = keV_FillQvsDepth[ikeVlo]
	Qmin = 4*PI*sin(theta)*keV/hc					// min Q (1/nm)
	name = ExamineImageSet#fileNameFromRangeNumbers(fileRoot,m1_FillQvsDepth[ikeVhi],m2_FillQvsDepth[ikeVhi])
	Wave image = $(LoadWinViewFile(name))			// load image
	geoLocal.ycent = NumberByKey("yc", note(image),"=")	// change yc in local copy of geo
	geoLocal.ycent = numtype(geoLocal.ycent) ? geo.ycent : geoLocal.ycent
	theta = imag(thetaRange(geoLocal,image))			// max theta on this image
	KillWaves/Z image
	keV = keV_FillQvsDepth[ikeVhi]
	Qmax = 4*PI*sin(theta)*keV/hc					// max Q (1/nm)

	// determine dQ (1/nm), the Q resolution to use.  Base it on the distance between two adjacent pixels
	Make/N=3/O/D qhat
	Variable px,py										// the current pixel to analyze (unbinned full chip pixels)
	py = round((starty+endy)/2)						// approximate full chip pixel position in center or the image
	px = round((startx+endx)/2)
	dQ = 4*PI*abs(sin(pixel2q(geo,px,py,qhat))-sin(pixel2q(geo,px,py+groupy,qhat)))*keV/hc
//		dQ /= 1.5
	NQ = round((Qmax-Qmin)/dQ) + 1					// number of Qs

	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
		Variable deV = 1000*(keV_FillQvsDepth[ikeVhi]-keV_FillQvsDepth[ikeVlo])/(NkeV-1)
		printf "depth range = [%g, %g] (µm),  ∆z=%.2g µm\r",depthMin,depthMax,dDepth
		printf "E range = [%g, %g] (keV),  ∆E=%.2g eV\r", keV_FillQvsDepth[ikeVlo], keV_FillQvsDepth[ikeVhi],deV
		printf "Q range = [%g, %g] (1/nm),  ∆Q=%.2g(1/nm)\r",Qmin, Qmax,dQ
		if (d0>0)
			printf "   d0 = %g (nm),   Q0 = %g (1/nm)\r",d0,Q0
		endif
	endif

	Make/N=(Ndepth,NQ)/O/D Q_Depth					// make the array to save the Q vs Depth array
	Make/N=(Ndepth,NQ)/O/D Q_DepthNorm			// holds number of pixels contributing to each element in Q_Depth[], use to normalize
	Q_Depth = 0
	Q_DepthNorm = 0
	Make/N=(NQ)/O/D Qhist								// array to hold Q's from one image
	Make/N=(NQ)/O/D QhistNorm						// number of pixels used for each Qhist, use to normalize
	SetScale/I x depthMin,depthMax,"µm", Q_Depth, Q_DepthNorm
	SetScale/I y Qmin,Qmax,"1/nm", Q_Depth, Q_DepthNorm
	SetScale/I x Qmin,Qmax,"1/nm", Qhist
	seconds = stopMSTimer(timer)/1e6
	if (ItemsInList(GetRTStackInfo(0))<2 && seconds>0.2)
		printf "setting up arrays took %s\r",Secs2Time(seconds,5,3)
	endif
	// done with the setup part, now actually compute something

	if (useDistortion)										// if using the distortion, precompute for all images  here
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
			timer = startMSTimer						// start a timer on creation of local distortion map
			KIllWaves/Z root:Packages:geometry:tempCachedDistortionMap	// ensure that this is gone!
			Make/N=(Ni,Nj,2)/O root:Packages:geometry:tempCachedDistortionMapTemp
			Wave distortionMap = root:Packages:geometry:tempCachedDistortionMapTemp
			SetScale/P x startx,groupx,"", distortionMap	// scale distortionMap[][] to full chip pixels
			SetScale/P y starty,groupy,"", distortionMap
			Variable/C dxy
			Wave xymap = root:Packages:geometry:xymap
			TextBox/C/N=text0/W=$gName  "making local copy of the distortion map"  ;  DoUpdate
			for (j=0;j<Nj;j+=1)
				for (i=0;i<Ni;i+=1)
					py = starty + j*groupy
					px = startx + i*groupx
					dxy = peakcorrection2(xymap,px,py)	// returns cmplx(dx,dy)
					distortionMap[i][j][0] = real(dxy)		// accumulate all the distortions that I will need
					distortionMap[i][j][1] = imag(dxy)
				endfor
			endfor
			Rename root:Packages:geometry:tempCachedDistortionMapTemp, tempCachedDistortionMap
			seconds = stopMSTimer(timer)/1e6
			if (ItemsInList(GetRTStackInfo(0))<2 && seconds>0.2)
				printf "creating local copy of distortion map took %s\r",Secs2Time(seconds,5,3)
			endif
		endif
		Note/K DistortionMap, "use=1"
	endif

	print "starting bulk of processing"
	timer = startMSTimer								// start a timer on bulk of processing
	Variable sec1=0, timer1, sec2=0,timer2, sec3=0,timer3
	Variable/G secExtra=0, secExtra1=0
	Variable depth
	// for all the N files (go over both range1 & range2), compute Qhist for each image
	for (m2=str2num(range2),j=0; numtype(m2)==0 || j==0; m2=NextInRange(range2,m2))		// loop over range2 (depth)
		name = ExamineImageSet#fileNameFromRangeNumbers(fileRoot,str2num(range1),m2)
		timer1=startMSTimer
		Wave sinTheta = $(MakeSinThetaArray(name,geo))// make an array the same size as an image, but filled with sin(theta) for this energy
		sec1 += stopMSTimer(timer1)/1e6
		Redimension/N=(Ni*Nj) sinTheta
		Make/N=(Ni*Nj)/I/O indexWaveQ
		indexWaveQ = p
		Sort sinTheta, sinTheta,indexWaveQ				// sort so indexWaveQ[0] is index to lowest sin(theta), indexWaveQ[inf] is greatest
		for (m1=str2num(range1); numtype(m1)==0; m1=NextInRange(range1,m1),j+=1)	// loop over range1 (energy)
			name = ExamineImageSet#fileNameFromRangeNumbers(fileRoot,m1,m2)
			if (mod(j,10)==0)
				TextBox/C/N=text0/W=$gName  "processing image "+num2istr(j)+" (of "+num2istr(N)+")"  ;  DoUpdate
			endif
			Qhist = 0										// needed because FillQhist1image() accumulates into Qhist
			QhistNorm = 0
		timer2=startMSTimer
			wnote = FillQhist1image(name,sinTheta,indexWaveQ,Qhist,QhistNorm,mask)	// fill Qhist from one file
		sec2 += stopMSTimer(timer2)/1e6
		timer3=startMSTimer
			if (NumberByKey("V_min", wnote,"=")==0  && NumberByKey("V_max", wnote,"=")==0)
				sec3 += stopMSTimer(timer3)/1e6
				continue									// no intensity here, so continue
			endif
			depth = NumberByKey("depthSi", wnote,"=")
			depth = numtype(depth) ? 0 : depth
			i = limit(round((depth-depthMin)/dDepth),0,Ndepth-1)
			i = numtype(i) ? 0 : i
			Q_Depth[i][] += Qhist[q]					// add Intensity from this image
			Q_DepthNorm[i][] += QhistNorm[q]		// accumulate number of pixels contibuting (not pixel intensity)
		sec3 += stopMSTimer(timer3)/1e6
			if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
				DoUpdate
			endif
		endfor
	endfor
	Note/K DistortionMap, "use=0"
	Q_Depth /= Q_DepthNorm							// do the normalization
	Q_Depth = numtype(Q_Depth[p][q]) ? NaN : Q_Depth[p][q]
	seconds = stopMSTimer(timer)/1e6
	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
		printf "\r  processing all %d images took %s	(%.3g µs/pixel)\r",N,Secs2Time(seconds,5,1),1e6*seconds/(N*Npixels)
		printf "		the sinTheta part took %s\r",Secs2Time(sec1,5,2)
		printf "		the FillQhist1image part took %s\r",Secs2Time(sec2,5,2)
		printf "		the accumulation/assignment part took %s\r",Secs2Time(sec3,5,2)
		printf "		the secExtra part took %s\r",Secs2Time(secExtra,5,2)
		printf "		the secExtra1 part took %s\r",Secs2Time(secExtra1,5,2)
	endif

	timer = startMSTimer								// start a timer on final stuff
	WaveStats/Q Q_depth
	wnote = ReplaceNumberByKey("ImaxPnt",wnote,V_maxloc,"=")
	ImageStats/Q Q_depth
	depth = V_maxRowLoc*DimDelta(Q_Depth,0) + DimOffset(Q_Depth,0)
	wnote = ReplaceNumberByKey("DepthOfMax",wnote,depth,"=")	// depth with the max point (µm)
	if (numtype(Q0)==0 && Q0>0)
		wnote = ReplaceNumberByKey("Q0",wnote,Q0,"=")
	endif
	String title = StrVarOrDefault("title","")
	if (strlen(title))
		title = ReplaceString("=",title,"_")				// in the wave note the string cannot have "=" or ";"
		title = ReplaceString(";",title,"_")
		wnote = ReplaceStringByKey("title",wnote,title,"=")
	endif

	if (WaveExists(mask))
		wnote = ReplaceStringByKey("maskWave",wnote,GetWavesDataFolder(mask,2),"=")
	endif
	wnote = ReplaceStringByKey("waveClass",wnote,"Qhistogram","=")
	// add geometry to the wave note
	wnote = ReplaceNumberByKey("geo_NxCCD",wnote,geo.NxCCD,"=")
	wnote = ReplaceNumberByKey("geo_NyCCD",wnote,geo.NyCCD,"=")
	wnote = ReplaceNumberByKey("geo_dpsx",wnote,geo.dpsx,"=")
	wnote = ReplaceNumberByKey("geo_dpsy",wnote,geo.dpsy,"=")
	wnote = ReplaceNumberByKey("geo_xcent",wnote,geo.xcent,"=")
	wnote = ReplaceNumberByKey("geo_ycent",wnote,geo.ycent,"=")
	wnote = ReplaceNumberByKey("geo_ddOffset",wnote,geo.ddOffset,"=")
	wnote = ReplaceNumberByKey("geo_dd",wnote,geo.dd,"=")
	wnote = ReplaceNumberByKey("geo_xbet",wnote,geo.xbet,"=")
	wnote = ReplaceNumberByKey("geo_xgam",wnote,geo.xgam,"=")
	wnote = ReplaceNumberByKey("geo_xalfd",wnote,geo.xalfd,"=")
	wnote = ReplaceNumberByKey("geo_xbetd",wnote,geo.xbetd,"=")
	wnote = ReplaceNumberByKey("geo_wire_F",wnote,geo.wire.F,"=")
	wnote = ReplaceNumberByKey("geo_wire_H0",wnote,geo.wire.H0,"=")
	wnote = ReplaceNumberByKey("geo_wire_Hyc",wnote,geo.wire.Hyc,"=")
	wnote = ReplaceNumberByKey("geo_wire_X",wnote,geo.wire.X,"=")
	wnote = ReplaceNumberByKey("geo_wire_dia",wnote,geo.wire.dia,"=")
	Note/K Q_depth, wnote
	ImageStats/M=1/Q Q_Depth
	Qhist = Q_Depth[V_maxRowLoc][p]					// add Intensity from this image
	Note/K Qhist, wnote
	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
		Make/T/O T_Constraints={"K0 > 0","K3 > 0"}	// ensure width>0, and base line positive
		WaveStats/Q Qhist
		T_Constraints[0] = "K0 > "+num2str(min(0,V_min))
		CurveFit/Q lor Qhist/D /C=T_Constraints 
		KillWaves/Z T_Constraints
		Variable fwhm = 2*sqrt(K3)					// for Lorentzian
		printf "at peak, Q = %.3f(1/nm),   fwhm = %.2g(1/nm)",K2,fwhm
		wnote = ReplaceNumberByKey("Qcenter",wnote,K2,"=")
		wnote = ReplaceNumberByKey("Qfwhm",wnote,fwhm,"=")
		if (Q0>0)
			Variable strain = -(K2-Q0)/Q0				// d is inverse of strain, so a negative
			printf ",   ∆Q = %.3f,   strain = %.1e",K2-Q0,strain
			wnote = ReplaceNumberByKey("∆Q",wnote,K2-Q0,"=")
			wnote = ReplaceNumberByKey("strain",wnote,strain,"=")
		endif
		Note/K Qhist, wnote
		printf "\r"
		MakeGraph_Q_Depth()
		MakeGraph_Qhist(Qhist)
	endif

	DoWindow/K $gName									// done with status window
	KillWaves/Z depth_FillQvsDepth, keV_FillQvsDepth, m1_FillQvsDepth, m2_FillQvsDepth, qhat, sinTheta, indexWaveQ, Q_DepthNorm, QhistNorm

	seconds = stopMSTimer(timer)/1e6
	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
		printf "  final stuff took %s\r",Secs2Time(seconds,5,2)
	endif
	seconds = stopMSTimer(timer0)/1e6
	if (ItemsInList(GetRTStackInfo(0))<2 && seconds>0.2)
		printf "entire process took %s\r",Secs2Time(seconds,5,3)
	endif
	beep
	return 0
End
//
Function MakeGraph_Q_Depth()				// display a Q_Depth[][] array
	Wave Q_Depth=Q_Depth
	if (!WaveExists(Q_Depth))
		DoAlert 0,"no 'Q_Depth' in this data folder"
		return 1
	endif
	String graphName=CleanupName("Graph_"+GetDataFolder(0)+"_Q_Depth",0)
	String wnote = note(Q_Depth)
	if (strlen(WinList(graphName,"","WIN:1")))
		DoWindow/F $graphName			// graph is already up, so just bring it to front
	elseif (exists(graphName)==5)	
		Execute graphName+"()"			// a recreation macro exists, run it
	else
		Display /W=(357,62,744,474)	// nothing exists, create the graph
		DoWindow/C $graphName
		AppendImage Q_Depth
		ModifyImage Q_Depth ctab= {*,*,Terrain,1}
		ModifyGraph gfMult=130, tick=2, mirror=1, minor=1, lowTrip=0.001, standoff=0
		ModifyGraph axOffset(left)=-1.5
//		String units = WaveUnits(Q_Depth,1)
//		if (strsearch(units,"1/",0)==0)
//			Label left "Q  (\\E"+units[2,Inf]+"\\S-1\\M)"
//		else
//			Label left "Q  (\\U)"
//		endif
//		Label bottom "Depth  (\\U)"
		Label left "Q"+FixUnitsInLabel(WaveUnits(Q_Depth,1))
		Label bottom "Depth"+FixUnitsInLabel(WaveUnits(Q_Depth,0))
		String str = StringByKey("title",wnote,"=")// first try to get title from wave note
		if (strlen(str)<1)
			str = StrVarOrDefault("title","")	// nothing in wavenote, try the title global string
		endif
		str = SelectString(strlen(str),"",str+"\r\\Zr067")
		str += GetWavesDataFolder(Q_Depth,1)+"\r"+StringByKey("dateExposed",wnote,"=")
		TextBox/C/N=text0/F=0 str
		Button SumQhistOverDepthButton,pos={100,0},size={70,16},fSize=9,proc=SumQhistOverDepthButtonProc,title="Sum Depths"
	endif

	Variable iTag = NumberByKey("ImaxPnt",wnote,"=")	// point in Q_Depth[][] of the max
	Variable Q0 = NumberByKey("Q0",wnote,"=")	// un-strained Q (1/nm)
	if (numtype(Q0)==0 && Q0>0)			// always redo the line showing Q0 if valid value is passed
		SetWindow $graphName userdata(Q0)="Q0="+num2str(Q0)+";"
	endif

	Variable x0=DimOffset(Q_Depth,0),x1,dx=DimDelta(Q_Depth,0), Nx=DimSize(Q_Depth,0)
	x1 = x0 + dx*(Nx-1)

	String fullName = GetWavesDataFolder(Q_Depth,1)+"Depth_onQplot"
	Variable/G $fullName = NumberByKey("DepthOfMax",wnote,"=")	// depth where the max is (µm)
	NVAR Depth_onQplot = $fullName
//	Variable/G :Depth_onQplot=NumberByKey("DepthOfMax",wnote,"=")	// depth where the max is (µm)
//	NVAR Depth_onQplot = :Depth_onQplot
	Depth_onQplot = numtype(Depth_onQplot) ? (x0 + round(Nx/2)*dx) : Depth_onQplot
	SetVariable QcutMarker,pos={1,2},size={95,15},proc=SetVarProcMakeQcut,title="Depth"
	SetVariable QcutMarker,limits={x0,x1,dx},value= $fullName
	SetVariable QcutMarker fSize=12, font="Lucida Grande"
	SetWindow kwTopWin hook(arrows)=arrowHook	// this arrow hook allows you to use arrow keys
	sprintf str, "value=%s;proc=SetVarProcMakeQcut;ctrlName=QcutMarker;stepSize=%g;loLim=%g;hiLim=%g",fullName,dx,x0,x1
	SetWindow kwTopWin userdata(arrows)=str
	if (numtype(iTag)==0 && iTag>0)			// for a valid number, redo the peak marker
		if (WhichListItem("peak0",AnnotationList(graphName))>=0 && stringmatch(StringByKey("TYPE",AnnotationInfo(graphName,"peak0")),"Tag"))
			Tag/C/N=peak0 Q_Depth, iTag		// simply move existing tag to peak
		else
			Tag/N=peak0/F=0/A=MB/X=-10/Y=20 Q_Depth, iTag, "\\OZ"
		endif
	endif
	SetVarProcMakeQcut("QcutMarker",Depth_onQplot,"","Depth_onQplot")
End
//
Function SetVarProcMakeQcut(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName
	if (!stringmatch(ctrlName,"QcutMarker"))
		return 1
	endif
	String wName = StringFromList(0,WinList("*",";","WIN:1"))
	String imageName = StringFromList(0,ImageNameList(wName,";"))
	Wave image = ImageNameToWaveRef(wName,imageName)
	String path = GetWavesDataFolder(image,1)
	Variable Q0 = NumberByKey("Q0",GetUserData(wName,"","Q0"),"=")

	SetDrawLayer/W=$wName/K UserFront
	SetDrawLayer/W=$wName UserFront
	SetDrawEnv/W=$wName xcoord= bottom,ycoord= prel,dash= 8
	DrawLine/W=$wName varNum,0,varNum,1
	if (numtype(Q0)==0)
		SetDrawEnv/W=$wName xcoord= prel,ycoord= left,dash= 1
		DrawLine/W=$wName 0,Q0,1,Q0
	endif
	Variable row = (varNum - DimOffset(image,0))/DimDelta(image,0)
	fitQat1depth(row,path)
	KillWaves/Z $(path+"W_coef")
	KillWaves/Z $(path+"W_sigma")
	KillWaves/Z $(path+"W_ParamConfidenceInterval")
End
//
Function SumQhistOverDepthButtonProc(B_Struct) : ButtonControl
	STRUCT WMButtonAction &B_Struct

	if (B_Struct.eventCode != 2)			// only process button-up
		return 0
	endif
	if (stringmatch(B_Struct.ctrlName,"SumQhistOverDepthButton"))
		String imageName = StringFromList(0,ImageNameList(B_Struct.win,";"))
		Wave image = ImageNameToWaveRef(B_Struct.win, imageName )
	endif
	if (!WaveExists(image))
		return 1
	endif
	Q_SumAllDepths(image)
	return 0
End
//
Function/T Q_SumAllDepths(image)			// sum all the depths in a Q vs Depth image
	Wave image

	String wname
	String list = WaveListClass("Qhistogram","*","DIMS:2")
	if (ItemsInList(list)==1)
		Wave image = $StringFromList(0,list)
	endif
	if (!WaveInClass(image,"Qhistogram"))
		Prompt wname,"2D wave of Q vs Depth",popup,WaveListClass("Qhistogram","*","DIMS:2")
		DoPrompt "Q vs Depth",wname
		if (V_flag)
			return ""
		endif
		Wave image = $wname
	endif
	wname=NameOfWave(image)
	if (!WaveInClass(image,"Qhistogram" ) || WaveDims(image)!=2)
		DoAlert 0,"cannot sum depths in "+wname+", it is not a 2D Q-histogram"
		return ""
	endif
	Duplicate/O image, Q_SumAllDepths_image
	Q_SumAllDepths_image = numtype(image) ? 0 : image

	ImageTransform sumAllCols Q_SumAllDepths_image
	KillWaves/Z Q_SumAllDepths_image
	Wave W_sumCols=W_sumCols
	SetScale/P x, DimOffset(image,1),DimDelta(image,1),WaveUnits(image,1),W_sumCols

	wname += "_SumDepths"
	Duplicate/O W_sumCols, $wname
	KillWaves/Z W_sumCols
	String wnote=ReplaceStringByKey("waveClass",note(image),"Qhistogram_SumDepth","=")
	Wave ww=$wname
	Note/K ww,wnote

	String graphName=StringFromList(0,FindGraphsWithWave(ww))
	if (strlen(graphName))									// bring graph with ww to the front
		DoWindow/F $graphName
		return GetWavesDataFolder(ww,2)
	endif

	// graph the wave
	graphName=CleanupName("Graph_"+GetDataFolder(0)+"_Q_SumDepths",0)
	Display/K=1/W=(732,69,1262,429) ww
	DoWindow/C $graphName
	ModifyGraph gfMult=130, tick=2, mirror=1, minor=1, lowTrip=0.001, standoff=0
	ModifyGraph rgb=(1,16384,41928), marker=19,msize=3,lsize=1, mode=4
	String labelValue = StringByKey("labelValue",wnote,"=")
	labelValue = SelectString(strlen(labelValue),"Intensity",labelValue)	// default label for Y-axis is "Intensity"
	Label left labelValue+FixUnitsInLabel(WaveUnits(ww,-1))
	String labelX = StringByKey("labelX",wnote,"=")
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
	str += GetWavesDataFolder(ww,1)+"\r"+StringByKey("dateExposed",wnote,"=")
	TextBox/A=RT/C/N=text0/F=0 str

	Variable Q0 = NumberByKey("Q0",wnote,"=")
	if (Q0>0)
		SetDrawLayer UserFront
		SetDrawEnv xcoord= bottom,ycoord= prel,dash= 1
		DrawLine Q0,0,Q0,1
	endif
	return GetWavesDataFolder(ww,2)
End
//
Function/S MakeGraph_Qhist(Qhist)	// display a Qhist[] wave
	Wave Qhist
	if (!WaveExists(Qhist))
		String QhistName = StringFromLIst(0,WaveListClass("Qhistogram","*","DIMS:1"))
		if (strlen(QhistName))
			Wave Qhist = $QhistName
		endif
	endif
	if (!WaveExists(Qhist))
		DoAlert 0, "no 'Qhist' wave exists in this data folder"
		return ""
	endif
	String graphName=CleanupName("Graph_"+GetDataFolder(0)+"_Qhist",0)
	if (strlen(WinList(graphName,"","WIN:1")))
		DoWindow/F $graphName			// graph is already up, so just bring it to front
		return graphName
	elseif (exists(graphName)==5)	
		Execute graphName+"()"			// a recreation macro exists, run it
		return graphName
	endif
	String fitName = GetWavesDataFolder(Qhist,1)+"fit_"+NameOfWave(Qhist)
	Display/K=1/W=(732,69,1262,429) Qhist
	DoWindow/C $graphName
	if (exists(fitName)==1)
		AppendToGraph $fitName
		ModifyGraph lSize($NameOfWave($fitName))=2
	endif
	String wnote = note(Qhist)
	ModifyGraph gfMult=130, tick=2, mirror=1, minor=1, lowTrip=0.001, standoff=0
	ModifyGraph mode(Qhist)=4, marker(Qhist)=19, lStyle(Qhist)=1, rgb(Qhist)=(16385,16388,65535), msize(Qhist)=4, hbFill(Qhist)=4
	ModifyGraph axOffset(left)=-3.85714,axOffset(bottom)=-0.461538
	String labelValue = StringByKey("labelValue",wnote,"=")
	labelValue = SelectString(strlen(labelValue),"Intensity",labelValue)	// default label for Y-axis is "Intensity"
	Label left labelValue+FixUnitsInLabel(WaveUnits(Qhist,-1))
	String labelX = StringByKey("labelX",wnote,"=")
	labelX = SelectString(strlen(labelX),"Q",labelX)			// default label for X-axis is "Q"
	Label bottom labelX+FixUnitsInLabel(WaveUnits(Qhist,0))
//	SetAxis/A/E=1/N=1 left
	SetAxis/N=1 left 0,*
	ShowInfo
	String str = StringByKey("title",wnote,"=")// first try to get title from wave note
	if (strlen(str)<1)
		str = StrVarOrDefault("title","")	// nothing in wavenote, try the title global string
	endif
	str = SelectString(strlen(str),"",str+"\r\\Zr075")
	str += "\\JR"
	str += GetWavesDataFolder(Qhist,1)+"\r"+StringByKey("dateExposed",note(Qhist),"=")
	TextBox/A=RT/C/N=text0/F=0 str

	Variable Qc = NumberByKey("Qcenter",wnote,"="), fwhm = NumberByKey("Qfwhm",wnote,"=")
	Variable dQ = NumberByKey("∆Q",wnote,"="), strain = NumberByKey("strain",wnote,"=")
	Variable Q0 = NumberByKey("Q0",wnote,"=")
	TextBox/C/N=textQ/F=0/A=LT/B=1 SelectString(numtype(Qc+fwhm),textQ_fromWnote(wnote),"")
//	if (numtype(Qc+fwhm)==0)				// found valid Q values
//		TextBox/C/N=textQ/F=0/A=LT/B=1 textQ_fromWnote(wnote)
//	endif
	if (Q0>0)
		SetDrawLayer UserFront
		SetDrawEnv xcoord= bottom,ycoord= prel,dash= 1
		DrawLine Q0,0,Q0,1
	endif
	Button SaveQhistTraceButton,pos={1,20},size={40,30},proc=SaveQhistTraceButtonProc,title="Keep\rTrace",fSize=10
	Button SaveQhistTraceButton,help={"duplicate the generic trace with a new name and leave it on plot"}
	return graphName
End
//
//Function MakeGraph_Qhist()				// display a Qhist[] wave
//	if (!WaveExists(Qhist))
//		DoAlert 0, "no 'Qhist' wave exists in this data folder"
//		return 1
//	endif
//	String graphName=CleanupName("Graph_"+GetDataFolder(0)+"_Qhist",0)
//	if (strlen(WinList(graphName,"","WIN:1")))
//		DoWindow/F $graphName			// graph is already up, so just bring it to front
//		return 0
//	elseif (exists(graphName)==5)	
//		Execute graphName+"()"			// a recreation macro exists, run it
//		return 0
//	endif
//	Variable fitExists = exists("fit_Qhist")==1
//	Display /W=(732,69,1262,429) Qhist
//	DoWindow/C $graphName
//	if (fitExists)
//		AppendToGraph fit_Qhist
//		ModifyGraph lSize(fit_Qhist)=2
//	endif
//	ModifyGraph gfMult=130, tick=2, mirror=1, minor=1, lowTrip=0.001, standoff=0
//	ModifyGraph mode(Qhist)=4, marker(Qhist)=19, lStyle(Qhist)=1, rgb(Qhist)=(16385,16388,65535), msize(Qhist)=4, hbFill(Qhist)=4
//	ModifyGraph axOffset(left)=-3.85714,axOffset(bottom)=-0.461538
//	Label left "Intensity"
//	String units = WaveUnits(Qhist,0)
//	if (strsearch(units,"1/",0)==0)
//		Label bottom "Q  (\\E"+units[2,Inf]+"\\S-1\\M)"
//	else
//		Label bottom "Q  (\\U)"
//	endif
//	SetAxis/A/E=1/N=1 left
//	ShowInfo
//	String wnote = note(Qhist)
//	String str = StringByKey("title",wnote,"=")// first try to get title from wave note
//	if (strlen(str)<1)
//		str = StrVarOrDefault("title","")	// nothing in wavenote, try the title global string
//	endif
//	str = SelectString(strlen(str),"",str+"\r\\Zr075")
//	str += GetWavesDataFolder(Qhist,1)+"\r"+StringByKey("dateExposed",note(Qhist),"=")
//	TextBox/A=RT/C/N=text0/F=0 str
//
//	Variable Qc = NumberByKey("Qcenter",wnote,"="), fwhm = NumberByKey("Qfwhm",wnote,"=")
//	Variable dQ = NumberByKey("∆Q",wnote,"="), strain = NumberByKey("strain",wnote,"=")
//	Variable Q0 = NumberByKey("Q0",wnote,"=")
//	if (numtype(Qc+fwhm)==0)				// found valid Q values
//		sprintf str, "Q\\Bc\\M = %g (nm\\S-1\\M)\rFWHM = %.2g",Qc,fwhm
//		TextBox/C/N=textQ/F=0/A=LT str
//		if (numtype(dQ+strain)==0)
//			sprintf str,"strain = %.1e",strain
//			AppendText/N=textQ str
//		endif
//	endif
//	if (Q0>0)
//		SetDrawLayer UserFront
//		SetDrawEnv xcoord= bottom,ycoord= prel,dash= 1
//		DrawLine Q0,0,Q0,1
//	endif
////	Button SaveQhistTraceButton,pos={1,2},size={75,20},proc=SaveQhistTraceButtonProc,title="Keep Trace"
////	Button SaveQhistTraceButton,help={"duplicate the generic trace with a new name and leave it on plot"}
//	Button SaveQhistTraceButton,pos={1,20},size={40,30},proc=SaveQhistTraceButtonProc,title="Keep\rTrace",fSize=10
//	Button SaveQhistTraceButton,help={"duplicate the generic trace with a new name and leave it on plot"}
//	return 0
//End
//
Function SaveQhistTraceButtonProc(ctrlName) : ButtonControl
	String ctrlName

	if (stringmatch(ctrlName,"SaveQhistTraceButton"))
		Wave w = TraceNameToWaveRef("","Qhist")
	endif
	if (!WaveExists(w))
		return 1
	endif

//	String win = WinList("*",";","WIN:1")

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
//Function SaveQhistTraceButtonProc(ctrlName) : ButtonControl
//	String ctrlName
//
//	if (stringmatch(ctrlName,"SaveQhistTraceButton"))
//		Wave w = TraceNameToWaveRef("","Qhist")
//	endif
//	if (!WaveExists(w))
//		return 1
//	endif
//	Variable row = NumberByKey("Qhist_row",note(w),"=")
//	String qName = GetWavesDataFolder(w,1)+"Q_Depth"
//	if (exists(qName)!=1)
//		return 1
//	endif
//
//	String trace = QDepthCutOnGraph("",$qName,row)
//	if (strlen(trace)<1)				// it is is not on the graph, put it there
//		AppendToGraph $qName[row][*]
//	endif
//	trace = QDepthCutOnGraph("",$qName,row)
//	if (strlen(trace))				// it is is not on the graph, put it there
//		ModifyGraph mode($trace)=4,marker($trace)=16,msize($trace)=2,lstyle($trace)=1,rgb($trace)=(9,37540,0)
//	endif
//	return 0
//End
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

	String str, textQ
	Variable Qc = NumberByKey("Qcenter",wnote,"=")
	Variable fwhm = NumberByKey("Qfwhm",wnote,"=")
	Variable strain = NumberByKey("strain",wnote,"=")

	sprintf textQ, "Q\\Bc\\M = %g (nm\\S-1\\M)\rFWHM = %.2g",Qc,fwhm
	if (numtype(strain)==0)
		sprintf str,"\rstrain = %.1e",strain
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


Static Function/S MakeSinThetaArray(fullName,geo)	// make an array the same size as an image, but filled with sin(theta)'s
	String fullName										// complete path name to image
	STRUCT microGeometry &geo							// the input geo with the yc of Si calibration standard

	Variable timer = startMSTimer						// start timer for initial processing
	Wave image = $(LoadWinViewFile(fullName))		// load image
	if (!WaveExists(image))
		printf "could not load image named '%s'\r",fullName
		DoAlert 0,"could not load image"
		return ""
	endif
	String wnote = note(image)
	Variable startx, starty, groupx, groupy, yc, ddLocal
	startx = NumberByKey("startx", wnote,"=")
	groupx = NumberByKey("groupx", wnote,"=")
	starty = NumberByKey("starty", wnote,"=")
	groupy = NumberByKey("groupy", wnote,"=")
	yc = NumberByKey("yc", wnote,"=")				// yc for this image
	ddLocal=NumberByKey("CCDy", note(image),"=") + geo.ddOffset
	if (numtype(startx+starty+groupx+groupy))
		DoAlert 0,"could not get ROI from wave note of image '"+fullName+"'"
		return ""
	endif
	String strStruct										// string to hold geo for an easy copying
	StructPut/S/B=0 geo, strStruct						// this is done to copy geo into geoLocal
	STRUCT microGeometry geoLocal						// the local geo with yc of this image
	StructGet/S/B=0 geoLocal, strStruct
	geoLocal.ycent = numtype(yc) ? geo.ycent : yc		// change yc in local copy of geo
	geoLocal.dd = !(ddLocal>0) ? geo.dd : ddLocal

	Variable Ni=DimSize(image,0), Nj=DimSize(image,1)
	Make/N=(Ni,Nj)/O/D sinThetaCached
	Note/K sinThetaCached, wnote
	CopyScales/I image, sinThetaCached
	sinThetaCached = NaN

	// for each pixel in image, compute sin(theta) and save it in sinThetaCached[][]
	Make/N=(Ni)/O/I sinThetaCachedPX
	Make/N=(Nj)/O/I sinThetaCachedPY
	sinThetaCachedPY = p*groupy + (groupy-1)/2 + (starty-1)
	sinThetaCachedPX = p*groupx + (groupx-1)/2 + (startx-1)
	sinThetaCached = pixel2sinEWspecial(geoLocal,sinThetaCachedPX[p],sinThetaCachedPY[q])
	KillWaves/Z sinThetaCachedPX,sinThetaCachedPY

	Variable seconds = stopMSTimer(timer)/1e6		// stop timer here
	if (ItemsInList(GetRTStackInfo(0))<3 && seconds>5)
		Variable depthSi = NumberByKey("depthSi", wnote,"=")
		printf "filling array of sin(theta) from  '%s'  at depth = %g,  took %s\r",NameOfWave(image),depthSi,Secs2Time(seconds,5,1)
	endif
	KIllWaves/Z image
	return GetWavesDataFolder(sinThetaCached,2)
End
//
// convert xp,yp positions on screen into absolute coordinates x,y,z (Wenge coords)
Static Function pixel2sinEWspecial(geo,px,py)
	STRUCT microGeometry &geo
	Variable px,py							// pixel position, 0 based, first pixel is (0,0), NOT (1,1)

	// convert px,py positions on detector into x,y,z (Wenge coords) with origin at Si position
	Wave distortionMap=root:Packages:geometry:tempCachedDistortionMap
	if (NumVarOrDefault("root:Packages:geometry:useDistortion",USE_DISTORTION_DEFAULT))// do not distort
		Variable i,j, Ni = DimSize(distortionMap,0), Nj = DimSize(distortionMap,1)
		i = (px+1 - DimOffset(distortionMap,0)) / DimDelta(distortionMap,0)
		j = (py+1 - DimOffset(distortionMap,1)) / DimDelta(distortionMap,1)
		if (i<0 || j<0 || i>=Ni || j>=Nj)
			Abort "peakcorrection2EWspecial(), 'root:Packages:geometry:tempCachedDistortionMap' has a scaling that does not fit input pixels"
		endif
		px += distortionMap[i][j][0]		// add correction and convert origin from (1,1) to (0,0)
		py += distortionMap[i][j][1]
	endif
	Variable xxst,yyst
	xxst = (1+px-geo.xcent)*geo.dpsx/geo.NxCCD	// the "1+" is to make 1 based pixels
	xxst = CCD_X_REVERSE ? -xxst : xxst
	yyst = -(1+py-geo.ycent)*geo.dpsy/geo.NyCCD
	Variable xcc,ycc,ddc,zcc
	xcc = geo.rded00*xxst + geo.rded10*yyst
	ycc = geo.rded01*xxst + geo.rded11*yyst
	zcc = geo.rded02*xxst + geo.rded12*yyst + geo.dd	// dd is up
	Variable dot = (xcc*geo.ki[0] + ycc*geo.ki[1] + zcc*geo.ki[2] ) / sqrt(xcc*xcc+ycc*ycc+zcc*zcc)
	return sqrt( (1-dot)/2 )				// sin(Bragg angle), uses half-angle formula
End
//Static Function/S MakeSinThetaArray(fullName,geo)		// make an array the same size as an image, but filled with sin(theta)'s
//	String fullName										// complete path name to image
//	STRUCT microGeometry &geo						// the input geo with the yc of Si calibration standard
//
//	Variable timer = startMSTimer						// start timer for initial processing
//	Wave image = $(LoadWinViewFile(fullName))		// load image
//	if (!WaveExists(image))
//		printf "could not load image named '%s'\r",fullName
//		DoAlert 0,"could not load image"
//		return ""
//	endif
//	String wnote = note(image)
//	Variable startx, starty, groupx, groupy, yc
//	startx = NumberByKey("startx", wnote,"=")
//	groupx = NumberByKey("groupx", wnote,"=")
//	starty = NumberByKey("starty", wnote,"=")
//	groupy = NumberByKey("groupy", wnote,"=")
//	yc = NumberByKey("yc", wnote,"=")			// yc for this image
//	if (numtype(startx+starty+groupx+groupy))
//		DoAlert 0,"could not get ROI from wave note of image '"+fullName+"'"
//		return ""
//	elseif (numtype(yc))
//		DoAlert 0,"invalid yc in image '"+fullName+"'"
//		return ""
//	endif
//	String strStruct									// string to hold geo for an easy copying
//	StructPut/S/B=0 geo, strStruct						// this is done to copy geo into geoLocal
//	STRUCT microGeometry geoLocal						// the local geo with yc of this image
//	StructGet/S/B=0 geoLocal, strStruct
//	geoLocal.ycent = yc									// change yc in local copy of geo
//
//	Make/N=(DimSize(image,0),DimSize(image,1))/O/D sinThetaCached
//	Note/K sinThetaCached, wnote
//	CopyScales/I image, sinThetaCached
//	sinThetaCached = NaN
//	// for each pixel in image, compute sin(theta) and save it in sinThetaCached[][]
//	Variable i, j, Ni=DimSize(image,0), Nj=DimSize(image,1)
//	Variable px,py										// the current pixel to analyze (unbinned full chip pixels)
////	Make/N=3/O/D qhat
//	for (j=0;j<Nj;j+=1)
//		py =  j*groupy + (groupy-1)/2 + (starty-1)
//		for (i=0;i<Ni;i+=1)
//			px = i*groupx + (groupx-1)/2 + (startx-1)
//			sinThetaCached[i][j] = sin(pixel2q(geoLocal,px,py,$""))
//		endfor
//	endfor
//	Variable seconds = stopMSTimer(timer)/1e6		// stop timer here
//	if (ItemsInList(GetRTStackInfo(0))<3 && seconds>5)
//		Variable depthSi = NumberByKey("depthSi", wnote,"=")
//		printf "filling array of sin(theta) from  '%s'  at depth = %g,  took %s\r",NameOfWave(image),depthSi,Secs2Time(seconds,5,1)
//	endif
//	KIllWaves/Z image, qhat
//	return GetWavesDataFolder(sinThetaCached,2)
//End



//=======================================================================================
//=======================================================================================
//=======================================================================================


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
Static Function/C thetaRange(geo,image)		// returns the range of theta spanned by image
	STRUCT microGeometry &geo
	Wave image

	String wnote = note(image)
	Variable startx, endx, starty, endy
	startx = NumberByKey("startx", wnote,"=")
	endx = NumberByKey("endx", wnote,"=")
	starty = NumberByKey("starty", wnote,"=")
	endy = NumberByKey("endy", wnote,"=")
	if (numtype(startx+endx+starty+endy))
		DoAlert 0, "could not get ROI from wave note of image '"+NameOfWave(image)+"'"
		return cmplx(NaN,NaN)
	endif

//	Make/N=3/O/D qhat
	Variable theta								// Bragg angle (radian)
	Variable thetaMin, thetaMax					// range of theta, test all 4 corners of the image

//	theta = pixel2q(geo,startx-1,starty-1,qhat)
	theta = pixel2q(geo,startx-1,starty-1,$"")
	thetaMin = theta
	thetaMax = theta

//	theta = pixel2q(geo,endx-1,starty-1,qhat)
	theta = pixel2q(geo,endx-1,starty-1,$"")
	thetaMin = min(theta,thetaMin)
	thetaMax = max(theta,thetaMax)

//	theta = pixel2q(geo,startx-1,endy-1,qhat)
	theta = pixel2q(geo,startx-1,endy-1,$"")
	thetaMin = min(theta,thetaMin)
	thetaMax = max(theta,thetaMax)

//	theta = pixel2q(geo,endx-1,endy-1,qhat)
	theta = pixel2q(geo,endx-1,endy-1,$"")
	thetaMin = min(theta,thetaMin)
	thetaMax = max(theta,thetaMax)

	//	print thetaMin*180/PI,thetaMax*180/PI
//	KillWaves/Z qhat
	return cmplx(thetaMin,thetaMax)
End


Static Function/S FillQhist1image(fullName,sinTheta,indexWaveQ,Qhist,QhistNorm,mask)// fill Qhist from one file, F is for fast, sin(theta) is precomputed
	String fullName
	Wave sinTheta										// array of sin(theta)'s that is same size as the image to be loaded	
	Wave indexWaveQ									// index sorted wave of sinTheta
	Wave Qhist
	Wave QhistNorm									// array that contains number of pixels contributing to each bin in Qhist
	Wave mask											// optional mask, only use pixels that are true, ONLY 0 or 1 (do not use 255)

	Variable timer = startMSTimer						// start timer for initial processing
	Wave image = $(LoadWinViewFile(fullName))		// load image
	if (!WaveExists(image))
		printf "could not load image named '%s'\r",fullName
//		DoAlert 0,"could not load image"
		timer = stopMSTimer(timer)					// stop timer
		return ""
	endif
	String wnote = note(image)
	Variable keV = NumberByKey("keV", wnote,"=")
	keV = NumberByKey("keV", wnote,"=")
	if (numtype(keV))
		DoAlert 0,"invalid keV in image '"+fullName+"'"
		timer = stopMSTimer(timer)					// stop timer
		return ""
	endif
	ImageStats/M=1/Q image							// check for empty images, skip them
	if (V_min==0 && V_max==0)
		timer = stopMSTimer(timer)					// stop timer
		wnote = ReplaceNumberByKey("V_min",wnote,V_min,"=")	// flag that this is empty
		wnote = ReplaceNumberByKey("V_max",wnote,V_max,"=")
		KillWaves/Z image
		return wnote									// image is empty, do not waste time adding zeros
	endif

	NVAR secExtra=secExtra, secExtra1=secExtra1
	if (NVAR_Exists(secExtra))
		Variable timeExtra = startMSTimer				// start timer for intermediate processing
	endif
	Variable i, j, Ni=DimSize(image,0), Nj=DimSize(image,1)
	Variable k											// index into Qhist for each pixel
	Variable NQ=numpnts(Qhist)
	Variable Qmin = DimOffset(Qhist,0)
	Variable dQ = DimDelta(Qhist,0)
	Variable factor = 4*PI*keV/hc
	MatrixOp/O kQimage = (factor*sinTheta-Qmin)/dQ	// this method is 4.5 x faster
	MatrixOp/O kQimage = round(kQimage)				// this is now the index into Qhist for each point
	if (WaveExists(mask))
		MatrixOp/O image = image*mask				// mask values needs to be 1 or 0
	endif
	Variable N=numpnts(image)
	Redimension/N=(N) image,kQimage					// unwrap for faster processing
	SetScale/P x 0,1,"", image, kQimage

	if (NVAR_Exists(secExtra1))
		Variable timeExtra1 = startMSTimer			// start timer for intermediate processing
	endif
	IndexSort indexWaveQ, image						// sinTheta was already sorted, so kQimage is already sorted
//	MatrixOp/O image = waveMap2(image,indexWaveQ,N)// this is slightly faster

	if (NVAR_Exists(secExtra1))
		secExtra1 += stopMSTimer(timeExtra1)/1e6
	endif
	// for each pixel in image, accumulate intensity into Qhist
	Variable i0=0,i1									// range of image to sum
	for (k=kQimage[0];k<NQ;)							// for each bin in Qhist
		i1 = BinarySearch(kQimage,k)					// this returns the last index with a value of k
		if (i1<i0)
			break
		endif
		Qhist[k] = sum(image,i0,i1)
//		QhistNorm[k] += i1-i0-1
		QhistNorm[k] += i1-i0+1						// number of pixels used for this Q bin
		i0 = I1 + 1										// reset start point of region
		k = kQimage[i0]
	endfor
	if (NVAR_Exists(secExtra))
		secExtra += stopMSTimer(timeExtra)/1e6
	endif

	Variable seconds = stopMSTimer(timer)/1e6		// stop timer here
	if (ItemsInList(GetRTStackInfo(0))<3 && seconds>0.5)
		printf "	processing image '%s' took %s\r",NameOfWave(image),Secs2Time(seconds,5,2)
	endif
	KIllWaves/Z image, kQimage
	return wnote
End
//Static Function/S xxFillQhist1image(fullName,geo,sinTheta,Qhist,QhistNorm,mask)// fill Qhist from one file, F is for fast, sin(theta) is precomputed
//	String fullName
//	STRUCT microGeometry &geo						// the input geo with the yc of Si calibration standard
//	Wave sinTheta										// array of sin(theta)'s that is same size as the image to be loaded	
//	Wave Qhist
//	Wave QhistNorm									// array that contains number of pixels contributing to each bin in Qhist
//	Wave mask											// optional mask, only use pixels that are true
//
//	Variable timer = startMSTimer						// start timer for initial processing
//	Variable sinThetaExists = WaveExists(sinTheta)		// use sinTheta[][] if it exists
//	Wave image = $(LoadWinViewFile(fullName))		// load image
//	if (!WaveExists(image))
//		printf "could not load image named '%s'\r",fullName
//		DoAlert 0,"could not load image"
//		timer = stopMSTimer(timer)					// stop timer
//		return ""
//	endif
//	String wnote = note(image)
//	Variable keV = NumberByKey("keV", wnote,"=")
//	Variable startx, starty, groupx, groupy, yc
//	startx = NumberByKey("startx", wnote,"=")
//	groupx = NumberByKey("groupx", wnote,"=")
//	starty = NumberByKey("starty", wnote,"=")
//	groupy = NumberByKey("groupy", wnote,"=")
//	yc = NumberByKey("yc", wnote,"=")				// yc for this image
//	keV = NumberByKey("keV", wnote,"=")
//	if (numtype(startx+starty+groupx+groupy))
//		DoAlert 0,"could not get ROI from wave note of image '"+fullName+"'"
//		timer = stopMSTimer(timer)					// stop timer
//		return ""
//	elseif (numtype(yc+keV))
//		DoAlert 0,"invalid yc or keV in image '"+fullName+"'"
//		timer = stopMSTimer(timer)					// stop timer
//		return ""
//	endif
//	ImageStats/M=1/Q image							// check for empty images, skip them
//	if (V_min==0 && V_max==0)
//		timer = stopMSTimer(timer)					// stop timer
//		wnote = ReplaceNumberByKey("V_min",wnote,V_min,"=")	// flag that this is empty
//		wnote = ReplaceNumberByKey("V_max",wnote,V_max,"=")
//		KillWaves/Z image
//		return wnote									// image is empty, do not waste time adding zeros
//	endif
//
//	Variable NQ=numpnts(Qhist)
//	Variable Qmin = DimOffset(Qhist,0)
//	Variable dQ = DimDelta(Qhist,0)
//	String strStruct									// string to hold geo for an easy copying
//	StructPut/S/B=0 geo, strStruct						// this is done to copy geo into geoLocal
//	STRUCT microGeometry geoLocal						// the local geo with yc of this image
//	StructGet/S/B=0 geoLocal, strStruct
//	geoLocal.ycent = yc									// change yc in local copy of geo
//
//	// for each pixel in image, accumulate intensity into Qhist
//	Variable haveMask = WaveExists(mask), useIt
//	Variable i, j, Ni=DimSize(image,0), Nj=DimSize(image,1)
//	Variable Qval										// Q value for one pixel
//	Variable k											// index into Qhist for each pixel
//	Variable px,py										// the current pixel to analyze (unbinned full chip pixels)
//	Make/N=3/O/D qhat
//	for (j=0;j<Nj;j+=1)
//		py =  j*groupy + (groupy-1)/2 + (starty-1)
//		for (i=0;i<Ni;i+=1)
////			useIt = image[i][j]>0						// use this pixel, check mask and note zero
////			useIt = (useIt && haveMask) ? mask[i][j] : useIt
//			useIt = haveMask ? mask[i][j] : 1			// do not use maked off pixels, but process even if pixel value is 0
//			if (useIt)									// do not spend time accumulating zeros
//				px = i*groupx + (groupx-1)/2 + (startx-1)
//				if (sinThetaExists)
//					Qval = 4*PI*sinTheta[i][j]*keV/hc	// sinTheta[][] exists, use it
//				else
//					Qval = 4*PI*sin(pixel2q(geoLocal,px,py,qhat))*keV/hc	// does not exist, so compute theta for each point
//				endif
////				k = limit(round((Qval-Qmin)/dQ),0,NQ-1)
//				k = round((Qval-Qmin)/dQ)
//				if (k<0 || k>=NQ)						// if k is out of range do not use it
//					continue
//				endif
//				Qhist[k] += image[i][j]
//				QhistNorm[k] += 1
//			endif
//		endfor
//	endfor
//	Variable seconds = stopMSTimer(timer)/1e6		// stop timer here
//	if (ItemsInList(GetRTStackInfo(0))<3 && seconds>0.5)
//		printf "	processing image '%s' took %s\r",NameOfWave(image),Secs2Time(seconds,5,2)
//	endif
//	KIllWaves/Z image, qhat
//	return wnote
//End


//=======================================================================================
//=======================================================================================
//=======================================================================================

// make the wave keV_Depth[][] which contains the average intensity/pixel for each image as a fn of energy and depth.
// it also returns the depth range (as two indices) that contain some intensity (returns NaN for no intensity at all)
Function/C Fill_EvsDepth(pathName,namePart)
	String pathName			// probably 'imagePath'
	String namePart			//	 = "EW5_"

	String str
	PathInfo $pathName
	if (!V_flag || strlen(namePart)<1)			// path does not exist or no namePart, ask user
		String pathPart
		str = ExamineImageSet#requestFileRoot(pathName,2)
		pathPart = StringFromList(0,str)
		namePart = StringFromList(1,str)
		if (strlen(pathPart)<1 || strlen(namePart)<1)
			return cmplx(NaN,NaN)			// invalid inputs
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
		return cmplx(NaN,NaN)				// invalid inputs
	endif

	str = ExamineImageSet#get_ranges(pathName,namePart)
	String range1=StringFromList(0,str), range2=StringFromList(1,str)
	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
		printf "using data from files starting with '%s'\r",namePart
		printf "setting energy range to '%s'\r",range1
		printf "setting depth range to '%s'\r",range2
	endif
	Variable NkeV = ItemsInRange(range1)
	Variable Ndepth = ItemsInRange(range2)

	Make/N=(Ndepth,NkeV)/O keV_Depth		// array to recieve the energy vs depth intensity
	keV_Depth = 0
	String fileRoot
	Variable i,m
	for (m=str2num(range1), i=0; numtype(m)==0; m=NextInRange(range1,m), i+=1)	// loop over all energies which is range1
		fileRoot = S_path+namePart+num2istr(m)+"_"
		ImageIntensityVsDepth(fileRoot,range2)
		DoUpdate
		Wave avgIntensity=avgIntensity
		keV_Depth[][i] += avgIntensity[p]
	endfor

	String wnote, name
	i = str2num(range2)
	sprintf name, "%s%s%d_%d%s",S_path,namePart,str2num(range1),i,IMAGE_FILE_EXT
	wnote=WinViewReadHeader(name)			// wave note to add to file read in
	Variable Elo = NumberByKey("keV", wnote,"=")	// energy for fiirst image
	sprintf name, "%s%s%d_%d%s",S_path,namePart,lastInRange(range1),i,IMAGE_FILE_EXT
	wnote=WinViewReadHeader(name)			// wave note to add to file read in
	Variable Ehi = NumberByKey("keV", wnote,"=")// energy for last image

	Wave depthIntensity=depthIntensity
	SetScale/I x depthIntensity[0],depthIntensity[Inf],"µm", keV_Depth
	SetScale/I y Elo,Ehi,"keV", keV_Depth

	// find range of depth in keV_Depth[][] that contains all the intensity
	ImageTransform sumAllRows keV_Depth		// this produces the wave W_sumRows[]
	Wave W_sumRows=W_sumRows
	// find the range in W_sumRows where there is intensity,  outside the range [i1,i2], W_sumRows[] is zero
	Variable i1,i2
	for (i1=0;i1<Ndepth && W_sumRows[i1]==0;i1+=1)		// find the first non-zero value at the beginning of W_sumRows[]
	endfor
	for (i2=Ndepth-1;i2>=i1 && W_sumRows[i2]==0;i2-=1)	// find the last non-zero values at the end of W_sumRows[]
	endfor
	Variable x1=i1*DimDelta(keV_Depth,0)+DimOffset(keV_Depth,0), x2=i2*DimDelta(keV_Depth,0)+DimOffset(keV_Depth,0)
	if (i1>i2)									// special for all zeros
		i1 = -1  ;  i2 = -1
		x1 = NaN  ;  x2 = NaN
	endif
	wnote = ReplaceNumberByKey("iDepthRangeLo",wnote,x1,"=")
	wnote = ReplaceNumberByKey("iDepthRangeHi",wnote,x2,"=")

	String title = StrVarOrDefault("title","")
	if (strlen(title))
		title = ReplaceString("=",title,"_")		// in the wave note the string cannot have "=" or ";"
		title = ReplaceString(";",title,"_")
		wnote = ReplaceStringByKey("title",wnote,title,"=")
	endif
	wnote = ReplaceStringByKey("waveClass",wnote,"EvsDepth","=")
	Note/K keV_Depth, wnote

	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
		printf "'%s', contains non-zero intensity in the depth range of [%d, %d] == (%g, %g)(µm)\r",GetWavesDataFolder(keV_Depth,2),i1,i2,x1,x2
		MakeGraph_keV_Depth()
	endif
	KillWaves/Z W_sumRows
	return cmplx(i1,i2)
End
//
Function MakeGraph_keV_Depth()				// display a keV_Depth[][] array
	Wave keV_Depth=keV_Depth
	if (!WaveExists(keV_Depth))
		DoAlert 0,"no energy depth data found"
		return 1
	endif
	String graphName=CleanupName("Graph_"+GetDataFolder(0)+"_keV_Depth",0)
	String wnote = note(keV_Depth)
	if (strlen(WinList(graphName,"","WIN:1")))
		DoWindow/F $graphName				// graph is already up, so just bring it to front
	elseif (exists(graphName)==5)	
		Execute graphName+"()"				// a recreation macro exists, run it
	else
		Display /W=(8,55,461,422)			// nothing exists, create the graph
		DoWindow/C $graphName
		AppendImage keV_Depth
		ModifyImage keV_Depth ctab= {*,*,Terrain,1}
		ModifyGraph gfMult=130, tick=2, mirror=1, minor=1, lowTrip=0.001, standoff=0
		ModifyGraph axOffset(left)=-1.71429,axOffset(bottom)=-0.384615
//		Label left "Energy  (\\U)"
//		Label bottom "Depth  (\\U)"
		Label left "Energy"+FixUnitsInLabel(WaveUnits(keV_Depth,1))
		Label bottom "Depth"+FixUnitsInLabel(WaveUnits(keV_Depth,0))
		String str = StringByKey("title",wnote,"=")// first try to get title from wave note
		if (strlen(str)<1)
			str = StrVarOrDefault("title","")	// nothing in wavenote, try the title global string
		endif
		str = SelectString(strlen(str),"",str+"\r\\Zr067")
		str += GetWavesDataFolder(keV_Depth,1)+"\r"+StringByKey("dateExposed",wnote,"=")
		TextBox/C/N=text0/F=0 str
	endif

	Variable x1 = NumberByKey("iDepthRangeLo",wnote,"=")	// range in depth that contains non-zero values
	Variable x2 = NumberByKey("iDepthRangeHi",wnote,"=")
	if (numtype(x1+x2)==0)					// always redo the lines showing the range, if x1 & x2 are valid
		SetDrawLayer/W=$graphName/K UserFront
		SetDrawLayer/W=$graphName UserFront
		SetDrawEnv/W=$graphName xcoord= bottom,ycoord= prel,dash= 1
		DrawLine/W=$graphName x1,0,x1,1
		SetDrawEnv/W=$graphName xcoord= bottom,ycoord= prel,dash= 1
		DrawLine/W=$graphName x2,0,x2,1
	endif
End

Function MakeLayoutQ_depth_hist(prefix)
	String prefix
	String gName="",gName2=""
	if (strlen(prefix)<1)
		Prompt gName, "name of graph for layout",popup,WinList("Graph*_Q_Depth",";","WIN:1")
		DoPrompt "Pick Graph",gName
		if (V_flag)
			return 1
		endif
		Variable i
		i = strsearch(gName, "Graph_",0)
		if (i!=0)
			return 1
		endif
		gName = gName[6,Inf]
		i = strsearch(gName,"_Q_Depth",Inf,3)
		prefix = gName[0,i-1]
	endif
	gName = "Graph_"+prefix+"_Q_Depth"
	gName2 = "Graph_"+prefix+"_Qhist"
	if (WinType(gName)!=1 || WinType(gName2)!=1)
		return 1
	endif
	NewLayout/C=1/W=(29,67,530,580)/P=portrait
	AppendLayoutObject/R=(65,505,555,739) graph  $gName2
	AppendLayoutObject/R=(53,23,553,503) graph $gName
	TextBox/N=stamp0/F=0/A=RB/X=0.17/Y=0.14 "\\Z06\\{\"%s %s\",date(), time()}"
	TextBox/N=stamp1/F=0/A=LB/X=0.17/Y=0.14 "\\Z06\\{\"%s\",CornerStamp1_()}:"+gName
End



// for a set of images over a range of depth indicies, make waves showing actual depth, peak intensity, and average intensity
Function ImageIntensityVsDepth(fileRoot,range)
	String fileRoot
	String range

	Variable N=ItemsInRange(range)			// number of points in range
	if (N<1)
		DoAlert 0, "There is nothing in this range,  range = '"+range+"'"
		return 0
	endif

	Variable depth1=NaN,depth2=NaN
	String name
	sprintf name, "%s%d%s",fileRoot,str2num(range),IMAGE_FILE_EXT	// open first file in range
	Wave image = $(LoadWinViewFile(name))
	if (!WaveExists(image))
		printf "could not load first image named '%s'\r",name
		Abort "could not load first image"
	endif
	depth1 = NumberByKey("depthSi", note(image),"=")
	KillWaves/Z image
	sprintf name, "%s%d%s",fileRoot,lastInRange(range),IMAGE_FILE_EXT	// open first file in range
	Wave image = $(LoadWinViewFile(name))
	if (!WaveExists(image))
		printf "could not load last image named '%s'\r",name
		Abort "could not load first image"
	endif
	depth2 = NumberByKey("depthSi", note(image),"=")
	KillWaves/Z image
	if (numtype(depth1+depth2))
		DoAlert 1, "Could not load depth information, this is probably not depth sorted images, contiuing?"
		if (V_flag!=1)
			return 0
		endif
	endif

	Make/N=(N)/O/D maxIntensity,avgIntensity,depthIntensity	// wave to hold the result
	maxIntensity = NaN
	avgIntensity = NaN
	depthIntensity = NaN

	Variable m, i
	for (m=str2num(range), i=0; numtype(m)==0; m=NextInRange(range,m), i+=1)	// loop over all values in range
		sprintf name, "%s%d%s",fileRoot,m,IMAGE_FILE_EXT
		Wave image = $(LoadWinViewFile(name))
		if (!WaveExists(image))
			continue							// just skip over missing images
			printf "could not load image named '%s',  continuing ...\r",name
		endif
		ImageStats/M=1/Q image
		maxIntensity[i] = V_max
		avgIntensity[i] = V_avg
		depthIntensity[i] = NumberByKey("depthSi", note(image),"=")
		KIllWaves/Z image
	endfor
	return N
End


Function NewEnergyWireScan(fldr,title,d0)		// set up a folder for a new Energy-Wire Scan
	String fldr									// new folder name, it will be at the top level
	String title									// a title to be used on plots
	Variable d0									// d-spacing of unstrained material (nm)

	fldr = StringFromList(ItemsInList(fldr,":")-1, fldr,":")	// remove leading part of path
	fldr = CleanupName(fldr,0)					// igorize the name
	String fldrSav= GetDataFolder(1)
	SetDataFolder root:
	fldr = SelectString(CheckName(fldr,11),fldr,"")	// avoid existing folders
	SetDataFolder $fldrSav
	d0 = (d0>0) ? d0 : NaN
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
		printf "\r•|\r•|\r"
		printf "NewEnergyWireScan(\"%s\",\"%s\",%.9g)\r",fldr,title,d0
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
	printf "•• changed data folder from  '%s'  -->  '%s'\r",fldrSav,GetDataFolder(1)
	if (d0>0)									// for a valid d0, create the global
		Variable/G $(fldr+":d0")=d0
	endif
	if (strlen(title))
		String/G $(fldr+":title")=title
	endif
	return 0
End


//=======================================================================================
//=======================================================================================
//===================  This section is used to make an image mask by hand

Function MakeMask(image)
	Wave image
	if (!WaveExists(image) || WaveDims(image)!=2)
		String imageName = ""
		String list = WaveListClass("speImage*;ImageSummed*","*","DIMS:2")
		if (ItemsInList(list)==1)
			imageName = StringFromList(0,list)
		else
			Prompt imageName,"image to use",popup,list
			DoPrompt/Help="3D-Xray Diffraction[Fit Peaks]" "image to mask",imageName
			if (V_flag)
				return 1
			endif
		endif
		Wave image = $imageName
		printf "MakeMask(%s)\r",imageName
	endif
	if (!WaveExists(image) || WaveDims(image)!=2)
		return 1
	endif

	Variable Nx=DimSize(image,0), Ny=DimSize(image,1)
	String wnote=ReplaceStringByKey("sourceImage","",GetWavesDataFolder(image,2),"=")
	String imName = GetWavesDataFolder(image,1)+CleanupName(NameOfWave(image)+"_IM",0)
	Duplicate/O image $imName
	Wave both = $imName
	wnote = ReplaceStringByKey("waveClass",wnote,"image&Mask","=")
	Note/K both, wnote

	String maskName = GetWavesDataFolder(image,1)+CleanupName(NameOfWave(image)+"Mask",0)
	Variable overWrite
	do
		overWrite=1
		if (exists(maskName)>0)
			DoAlert 1, " The wave '"+maskName+"' alreay exists, Edit it?"
			overWrite = (V_flag==1)
			if (!overWrite)
				Prompt maskName, "name of wave with mask"
				DoPrompt "mask wave",maskName
				if (V_flag)
					KillWaves/Z $imName
					return 1
				endif
			endif
		endif
	while(!overWrite)
	Make/N=(Nx,Ny)/B/U/O $maskName
	Wave mask = $maskName
	wnote = ReplaceStringByKey("waveClass",wnote,"imageMask","=")
	Note/K mask, wnote

	//	WaveStats/Q image
	//	Variable delta = (V_max-V_min)/30			// set both to show existing mask (if one had been set)
	Variable delta = deltaForMaskOffset(image)
	both = mask[p][q] ? image[p][q]+delta : image[p][q]

	Display/W=(345,44,822,440)/K=1
	AppendImage both
	Execute "GraphImageStyle()"
	ModifyGraph axOffset(left)=1
	Button NewFreeHandButton,pos={3,2},size={55,30},proc=MaskFreeHandButtonProc,title="New\rFreehand"
	Button NewFreeHandButton,help={"make a new freehand region to add (or subtract) from mask"}
	Button NewFreeHandButton,fSize=10
	Button AddFreeHandButton,pos={3,34},size={55,30},disable=1,proc=MaskFreeHandButtonProc,title="Add\rFreehand"
	Button AddFreeHandButton,help={"add the freehand region to mask"},fSize=10
	Button SubtractFreeHandButton,pos={3,66},size={55,30},disable=1,proc=MaskFreeHandButtonProc,title="Subtract\rFreehand"
	Button SubtractFreeHandButton,help={"subtract the freehand region to mask"},fSize=10
	Button KillFreeHandButton,pos={3,98},size={55,30},disable=1,proc=MaskFreeHandButtonProc,title="Kill\rFreehand"
	Button KillFreeHandButton,help={"kill the freehand region, do nothing to mask"},fSize=10
	DoWindow/T $WinName(0,1),"Close Me When Done"
	SetWindow kwTopWin hook(setMask)=setMaskImageHook
	SetWindow kwTopWin userdata(images)="image="+GetWavesDataFolder(image,2)+";mask="+maskName+";both="+imName+";"
	SetWindow kwTopWin userdata(deleteOnClose)=imName+";"
	return 0
End
//
Function setMaskImageHook(s)		// things to do when the mask edit window is closed
	STRUCT WMWinHookStruct &s
	if (s.eventCode!=2)
		return 0
	endif
	String win= s.winName
	String item,killList = GetUserData(s.winName,"","deleteOnClose")
	Variable i, N=ItemsInList(killList)
	for (i=0;i<N;i+=1)
		item = StringFromList(i,killList)
		Wave image = $item
		RemoveImage/Z/W=$win $NameOfWave(image)
		KillWaves/Z $item
	endfor

	Wave image = $StringByKey("image",GetUserData("","","images"),"=")
	Wave mask = $StringByKey("mask",GetUserData("","","images"),"=")
	if (!WaveExists(image) || !WaveExists(mask))
		return 1
	endif
	String wnote = note(image)
	wnote = ReplaceStringByKey("mask",wnote,GetWavesDataFolder(mask,2),"=")
	Note/K image,wnote
	return 1
End
//
Function AddMarqueeToMask()
	Wave image = $StringByKey("image",GetUserData("","","images"),"=")
	Wave mask = $StringByKey("mask",GetUserData("","","images"),"=")
	Wave both = $StringByKey("both",GetUserData("","","images"),"=")
	if (!WaveExists(image) || !WaveExists(mask) || !WaveExists(both))
		DoAlert 0, "Nothing done, this is probably the wrong kind of image, see 'Make Mask'"
		return 1
	endif
	GetMarquee/K left,bottom
	if (!V_flag)
		return 1
	endif
	Variable x0 = DimOffset(both,0), dx = DimDelta(both,0)
	Variable y0 = DimOffset(both,1), dy = DimDelta(both,1)
	Variable ilo = round(min((V_left-x0)/dx, (V_right)/dx))
	Variable ihi = round(max((V_left-x0)/dx, (V_right)/dx))
	Variable jlo = round(min((V_top-x0)/dx, (V_bottom)/dx))
	Variable jhi = round(max((V_top-x0)/dx, (V_bottom)/dx))
	mask[ilo,ihi][jlo,jhi] = 1
	Variable delta = deltaForMaskOffset(image)
	both = mask[p][q] ? image[p][q]+delta : image[p][q]
End
//
Function SubtractMarqueeFromMask()
	Wave image = $StringByKey("image",GetUserData("","","images"),"=")
	Wave mask = $StringByKey("mask",GetUserData("","","images"),"=")
	Wave both = $StringByKey("both",GetUserData("","","images"),"=")
	if (!WaveExists(image) || !WaveExists(mask) || !WaveExists(both))
		DoAlert 0, "Nothing done, this is probably the wrong kind of image, see 'Make Mask'"
		return 1
	endif
	GetMarquee/K left,bottom
	if (!V_flag)
		return 1
	endif
	Variable x0 = DimOffset(both,0), dx = DimDelta(both,0)
	Variable y0 = DimOffset(both,1), dy = DimDelta(both,1)
	Variable ilo = round(min((V_left-x0)/dx, (V_right)/dx))
	Variable ihi = round(max((V_left-x0)/dx, (V_right)/dx))
	Variable jlo = round(min((V_top-x0)/dx, (V_bottom)/dx))
	Variable jhi = round(max((V_top-x0)/dx, (V_bottom)/dx))
	mask[ilo,ihi][jlo,jhi] = 0
	Variable delta = deltaForMaskOffset(image)
	both = mask[p][q] ? image[p][q]+delta : image[p][q]
End
//
Function ClearMask()
	String win = StringFromList(0,WinList("*",";","WIN:1"))
	Wave image = $StringByKey("image",GetUserData(win,"","images"),"=")
	Wave mask = $StringByKey("mask",GetUserData(win,"","images"),"=")
	Wave both = $StringByKey("both",GetUserData(win,"","images"),"=")
	if (!WaveExists(image) || !WaveExists(mask) || !WaveExists(both))
		DoAlert 0, "Nothing done, this is probably the wrong kind of image, see 'Make Mask'"
		return 1
	endif
	GetMarquee/Z/K
	mask = 0
	both = image
End

//		Make/N=(Nx,Ny)/O/U/B imageSeedTemp
//
//		GraphWaveDraw/W=Graph0/O yWavName, xWavName
//		GraphWaveEdit/W=Graph0 yWavName
//		GraphNormal/W=Graph0

Function MaskFreeHandButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	if (ba.eventCode!=2)									// only process on mouse-up
		return 0
	endif
	String win=ba.win

	String xName,yName
	if (stringmatch(ba.ctrlName,"NewFreeHandButton"))
		Button NewFreeHandButton win=$win, disable=2	// disable new freehand button
		Button AddFreeHandButton win=$win, disable=0		// enable all other buttons
		Button SubtractFreeHandButton win=$win, disable=0
		Button KillFreeHandButton win=$win, disable=0
		Make/N=1/O yWavFreeHandMask=NaN, xWavFreeHandMask=NaN
		xName=GetWavesDataFolder(xWavFreeHandMask,2)
		yName=GetWavesDataFolder(yWavFreeHandMask,2)

		String list = GetUserData(win,"","deleteOnClose")
		list += SelectString(WhichListItem(xName,list)<0,"",xName+";")
		list += SelectString(WhichListItem(yName,list)<0,"",yName+";")
		SetWindow kwTopWin,userdata(deleteOnClose)=list

		list = GetUserData(win,"","images")
		list += SelectString(WhichListItem(xName,list)<0,"","xFreeHand="+xName+";")
		list += SelectString(WhichListItem(yName,list)<0,"","yFreeHand="+yName+";")
		SetWindow kwTopWin,userdata(images)=list

		GraphWaveDraw/W=$win/O yWavFreeHandMask, xWavFreeHandMask
		DoAlert 0,"Draw and edit a freehand curve on the graph to add/subtract from mask"
		return 0
	endif

	GraphNormal/W=$win
	Button NewFreeHandButton win=$win, disable=0		// enable new freehand button
	Button AddFreeHandButton win=$win, disable=1			// hide all other buttons
	Button SubtractFreeHandButton win=$win, disable=1
	Button KillFreeHandButton win=$win, disable=1

	Wave xw=$StringByKey("xFreeHand",GetUserData(win,"","images"),"=")
	Wave yw=$StringByKey("yFreeHand",GetUserData(win,"","images"),"=")
	yName = NameOfWave(yw)
	if (WhichListItem(yName,TraceNameList(win,";",1))>=0)
		RemoveFromGraph/W=$win $yName	
	endif

	if (stringmatch(ba.ctrlName,"KillFreeHandButton"))
		return 0
	endif

	Wave image = $StringByKey("image",GetUserData(win,"","images"),"=")
	Wave mask = $StringByKey("mask",GetUserData(win,"","images"),"=")
	Wave both = $StringByKey("both",GetUserData(win,"","images"),"=")
	if (!WaveExists(image) || !WaveExists(mask) || !WaveExists(both))
		DoAlert 0, "Nothing done to mask, this is probably the wrong kind of image, see 'Make Mask'"
		return 0
	endif
	Wave tempMask = $PixelInsideCurve(image,yw,xw)				// get temp mask, need to add/subtract this
	if (!WaveExists(tempMask))
		DoAlert 0,"Mask un-changed"
		return 0
	endif
	if (stringmatch(ba.ctrlName,"AddFreeHandButton"))
		mask = mask || tempMask							// combine the two masks
	elseif (stringmatch(ba.ctrlName,"SubtractFreeHandButton"))
		mask = mask && !tempMask
	endif
	KillWaves/Z tempMask
	Variable delta = deltaForMaskOffset(image)
	both = mask ? image+delta : image
	return 0
End
//
Static Function/T PixelInsideCurve(image,yw,xw)
	Wave image
	Wave yw,xw		// curve defining the ROI to add or subtract

	Variable Nx=DimSize(image,0), Ny=DimSize(image,1)
	if( !(Nx>2 && Ny>2) )
		return ""
	endif

	Variable i,N=DimSize(yw,0)
	if ((abs(xw[0]-xw[N-1])+abs(yw[0]-yw[N-1]))>4)	// add extra end point to close region
		Redimension/N=(N+1) yw,xw
		xw[N] = xw[0]
		yw[N] = yw[0]
		N += 1
	else												// ensure that region exactly closes
		xw[N-1] = xw[0]
		yw[N-1] = yw[0]
	endif

	Make/N=(Nx,Ny)/O/U/B PixelInsideMaskJZT=0
	Wave mask=PixelInsideMaskJZT

	Make/N=(N+1)/O PixelInsideYcross
	Wave yCross=PixelInsideYcross					// list of y-crossing values from line segments
	Variable nCross=0								// number of crossings saved in yCross

	WaveStats/M=1/Q xw
	Variable pxLo=floor(V_min), pxHi=ceil(V_max)
	pxLo = limit(pxLo,0,Nx-1)
	pxHi = limit(pxHi,0,Nx-1)
	Variable px,py									// pixel in image
	Variable x0,y0, x1,y1							// ends of a line segment
	Variable slope, yy, dx							// slope of line segement
	Variable outside, i0,i1

	for (px=pxLo;px<=pxHi;px+=1)
		yCross = NaN								// find evry line segment that crosses this px value
		nCross = 0
		x1=xw[0]+enoise(0.01);	 y1=yw[0]
		for (i=1;i<N;i+=1)							// check all line segments to see if they intersect this px
			x0 = x1;		y0 = y1	
			x1 = xw[i]+enoise(0.01);	y1=yw[i]
			if ((x0-px)*(x1-px)<0)				// segment crosses px, solve for y, and save it
				dx = x1-x0
				if (dx==0)							// skip vertical lines
					continue
				endif
				slope = (y1-y0)/dx				// slope of line segment
				yy = slope*px+(y0-slope*x0)		// segement intersects px,  y0-slope*x0==b
				yCross[nCross] = round(yy)
				nCross += 1
			endif
		endfor

		Sort yCross, yCross							// the NaN's are at the end
		outside = 0									// start from inside
		for (i=0,x1=-1;i<nCross;i+=1)
			i0 = i1
			i1 = yCross[i]
			if (outside)								// was outside, so this region is inside
				mask[px][i0,i1] = 1
			endif
			outside = !outside
		endfor
	endfor

	for (i=0,x1=-1;i<N;i+=1)						// set all verticies and set all vertical line segments
		x0 = x1;		y0 = y1	
		x1 = xw[i];		y1=yw[i]
		dx = x1-x0
		if (abs(dx)<0.5)							// set all of one vertical line segement
			px = round(x0)
			i0 = round(y0)
			i1 = round(y1)
			mask[px][i0,i1] = 1
		else	
			mask[round(x1)][round(y1)] = 1		// set every vertex
		endif
	endfor

	KillWaves/Z PixelInsideYcross
	return GetWavesDataFolder(PixelInsideMaskJZT,2)
End
//
Static Function deltaForMaskOffset(image)		// set both to show existing mask (if one had been set)
	Wave image
	if (!WaveExists(image))
		return NaN
	endif
	WaveStats/Q image
	Variable delta = V_avg>0 ? V_sdev/3 : -V_sdev/3	//	was:     delta = (V_max-V_min)/30
	return delta
End

//=======================================================================================
//=======================================================================================


Function MakePseudoWhiteImages(pathName,namePart)
	String pathName			// probably "imagePath"
	String namePart			//	 = "EW5_"

	pathName = SelectString(strlen(pathName),"imagePath",pathName)
	PathInfo $pathName
	String str
	if (!V_flag || strlen(namePart)<1)			// path does not exist or no namePart, ask user
		String pathPart
		str = ExamineImageSet#requestFileRoot(pathName,2)
		pathPart = StringFromList(0,str)
		namePart = StringFromList(1,str)
		if (strlen(pathPart)<1 || strlen(namePart)<1)
			return 1							// invalid inputs
		endif
		if (!stringmatch(pathPart,S_path))		// path was changed
			if (stringmatch(pathName,"imagePath"))
				NewPath/O/M="path to reconstructed image files" imagePath pathPart	// for iamgePath, automatically reassign
			else
				NewPath/M="path to reconstructed spe files" $pathName pathPart	// for other names, ask
			endif
		endif
		printf "MakePseudoWhiteImages(\"%s\",\"%s\")\r",pathName,namePart
	endif
	PathInfo $pathName
	if (strlen(S_path)<1 || strlen(namePart)<1)
		return 1								// invalid inputs
	endif
	String name,fileRoot=S_path+namePart

	Variable timer=startMSTimer, timer1, seconds
	str = ExamineImageSet#get_ranges(pathName,namePart)
	String range1=StringFromList(0,str), range2=StringFromList(1,str)
	Prompt range1,"range of energies, file numbers"
	Prompt range2,"range of depths (or other motion), file numbers"
	DoPrompt "ranges",range1, range2
	if (V_flag)
		return 1
	endif
	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
		printf "using data from files starting with '%s'\r",namePart
		printf "setting energy range to '%s'\r",range1
		printf "setting depth range to '%s'\r",range2
	endif

	Variable NkeV = max(ItemsInRange(range1),1)
	Variable Ndepth = max(ItemsInRange(range2),1)

	Variable iDepth=str2num(range2), jkeV=str2num(range1)
	if (!(jkeV>=0))
		sprintf name,"%s%d%s",fileRoot,iDepth,IMAGE_FILE_EXT
	elseif (!(iDepth>=0))
		sprintf name,"%s%d%s",fileRoot,jkeV,IMAGE_FILE_EXT
	else
		sprintf name,"%s%d_%d%s",fileRoot,jkeV,iDepth,IMAGE_FILE_EXT
	endif
	Wave image = $(LoadWinViewFile(name))
	if (!WaveExists(image))
		printf "could not load first image named '%s'\r",name
		Abort "could not load first image"
	endif
	Variable Nx=DimSize(image,0),  Ny=DimSize(image,1)
	Variable tooBig = (Nx*Ny*Ndepth > 1e8)
	if (tooBig)
		print "There is not enough memory available to make the individual sum images, so no wave 'imageEsum'"
		Make/N=(Nx,Ny)/O/D imageSumAll
		imageSumAll = 0
	else
		Make/N=(Nx,Ny,Ndepth)/O/U/I imageEsum
		if (WaveType(image)<=4)				// use double
			Redimension/D imageEsum
		endif
		imageEsum = 0
	endif
	Make/N=(Ndepth)/O/D sumIntensDepth
	SetScale/P x 0,1,"", sumIntensDepth
	sumIntensDepth = 0
	Variable depthHi,depthLo = NumberByKey("depthSi", note(image),"=")	// depth
	KillWaves/Z imagePlane
	Rename image imagePlane					// need this for the movie
	if (WaveType(imagePlane ) & 0x18)
		Redimension/S imagePlane
	endif

	Variable j, i
	for (j=0,jkeV=str2num(range1); j<NkeV; j+=1,jkeV=NextInRange(range1,jkeV))
		for (i=0,iDepth=str2num(range2);i<Ndepth;i+=1,iDepth=NextInRange(range2,iDepth))
			if (i==0 && j==0)
				timer1 = startMSTimer			// time one image to get estimated execution time
			endif
			if (!(jkeV>=0))
				sprintf name,"%s%d%s",fileRoot,iDepth,IMAGE_FILE_EXT
			elseif (!(iDepth>=0))
				sprintf name,"%s%d%s",fileRoot,jkeV,IMAGE_FILE_EXT
			else
				sprintf name,"%s%d_%d%s",fileRoot,jkeV,iDepth,IMAGE_FILE_EXT
			endif
			Wave image = $(LoadWinViewFile(name))
			if (tooBig)
				imageSumAll += image
			else
				imageEsum[][][i] += image[p][q]
			endif
			sumIntensDepth[i] += sum(image)
			depthHi = NumberByKey("depthSi", note(image),"=")	// depth
			KillWaves/Z image
			if (i==0 && j==0)
				seconds = stopMSTimer(timer1)*1e-6
				seconds *= (NkeV*Ndepth-1)	// total estimated time
				printf "should run for an additional %s, which will be at %s\r",Secs2Time(seconds,5,0),Secs2Time(60*round((DateTime+seconds)/60),1,0)
			endif
		endfor
	endfor
	if (!tooBig)
		ImageTransform sumPlanes imageEsum		// sum over all energies
		Duplicate/O M_SumPlanes imageSumAll
		KillWaves/Z M_SumPlanes
		SetScale/P x 0,1,"pixel", imageEsum
		SetScale/P y 0,1,"pixel", imageEsum
		if (numtype(depthHi+depthLo)==0)
			SetScale/I z depthLo,depthHi,"µm", imageEsum
		endif
	endif
	SetScale/P x 0,1,"pixel", imageSumAll
	SetScale/P y 0,1,"pixel", imageSumAll
	if (numtype(depthHi+depthLo)==0)
		SetScale/I x depthLo,depthHi,"µm", sumIntensDepth
	endif
	seconds = stopMSTimer(timer)*1e-6
	printf "total execution time was %s\r",Secs2Time(seconds,5,1)

	MakeEsumPlot()
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
	ImageDisplayScaling# ModifyOnly_ctab_range("","imagePlane",V_min,V_max*0.7)
	AddMovieFrame

	WaveStats/M=1/Q imageEsum
	ImageDisplayScaling# ModifyOnly_ctab_range("","imagePlane",0,V_max/3)


	Variable i, N=DimSize(imageEsum,2)
	for (i=0;i<N;i+=1)
		SetEsumDepthProc("EsumDepthDisp",i,"","")
		AddMovieFrame
	endfor
	SetEsumDepthProc("EsumDepthDisp",-1,"","")
	WaveStats/M=1/Q imageSumAll
	ImageDisplayScaling# ModifyOnly_ctab_range("","imagePlane",V_min,V_max*0.7)
	AddMovieFrame
	CloseMovie
	SetEsumDepthProc("EsumDepthDisp",NumVarOrDefault("EsumDepth",-1),"","")

//	DoAlert 1,"play the movie?"
//	if (V_flag==1)
//		PlayMovie [/I/M/P=pathName /W=(left, top, right, bottom) /Z] [as fileNameStr ]
	//endif
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
		ImageDisplayScaling# ModifyOnly_ctab_range("","imagePlane",V_min,V_max*0.7)

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
		get_ctab_range("",NameOfWave(imageplane),lo,hi)
		sprintf list,"lo=%g;hi=%g;sum=1",lo,hi
		SetWindow kwTopWin userdata(ctab)=list
		imageplane = imageSumAll
		WaveStats/M=1/Q imageSumAll
		ImageDisplayScaling# ModifyOnly_ctab_range("","imagePlane",V_min,V_max*0.7)
		TextBox/C/N=textEsumDepth/F=0/S=3/A=LT/X=3.14/Y=3.46 "\\F'symbol'\\Zr200S\\Zr050\\F]0(Energys)   \\F'symbol'\\Zr200S\\Zr050\\F]0(Depths)"
	else
		list = GetUserData("","","ctab")
		if (NumberByKey("sum",list,"="))		// coming off of a sum frame, so reset ctab limits to stored values
			lo = NumberByKey("lo",list,"=")
			hi = NumberByKey("hi",list,"=")
			list = ReplaceNumberByKey("sum",list,0,"=")
			SetWindow kwTopWin userdata(ctab)=list
			ImageDisplayScaling# ModifyOnly_ctab_range("","imagePlane",lo,hi)
		endif
		imageplane = imageEsum[p][q][i]
		depth = DimOffset(imageEsum,2) + i*DimDelta(imageEsum,2)
		TextBox/C/N=textEsumDepth/F=0/S=3/A=LT/X=3.14/Y=3.46 "\\F'symbol'\\Zr200S\\Zr050\\F]0(Energys)  Depth = "+num2str(depth)+" µm"
	endif
	String title = StrVarOrDefault(path+"title","")
	if (strlen(title))
		AppendText/N=textEsumDepth "\\Zr130"+title
	endif
	AppendText/N=textEsumDepth "\\Zr060"+path
End


//=======================================================================================
//=======================================================================================


Function arrowHook(s)				// for using arrow keys in addition to the little up/down arrows on a SetVariable control
	STRUCT WMWinHookStruct &s
	if (s.eventCode!=11 || s.eventMod)
		return 0
	endif

	Variable step=1
	if (s.keycode==29 || s.keycode==30)		// 29 is right,  30 is up	make bigger
		step = 1
	elseif (s.keycode==28 || s.keycode==31)	// 29 is leftt,  31 is down	make smaller
		step = -1
	else
		return 0
	endif
	String win=s.winName
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
//
//	SetWindow kwTopWin hook(arrows)=arrowHook
//	SetWindow kwTopWin userdata(arrows)="value=root:test:EsumDepth;proc=SetEsumDepthProc;stepSize=1;hiLim=71;loLim=-1;"



Function/T SortNumericalList(list,sep,increasing)		// return a list of number sorted numerically
	String list								// a list of numbers
	String sep								// separator used for this list
	Variable increasing						// true=increaing order,   false=decreasing order

	Variable i, N = ItemsInList(list,sep)

	Make/N=(N)/D N_SortNumericalList	// need a wave to do sorting
	Wave nw = N_SortNumericalList

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
	KillWaves/Z N_SortNumericalList
	return list
End





//Function FillGeometryStructDefault(geo)				//fill the geometry structure with test values
//	STRUCT microGeometry &geo
//
//	String strStruct=StrVarOrDefault(":geoStructStr","")// set to values in current directory
//	if (strlen(strStruct)<1)
//		strStruct=StrVarOrDefault("root:Packages:geometry:geoStructStr","")	// try the default values
//	endif
//	if (strlen(strStruct)>1)
//		StructGet/S/B=2 geo, strStruct					// found structure information, load into geo
//	else
//		LoadPackagePreferences/MIS=1 "microGeo","microGeoPrefs",0,geo
//		if (V_flag)
//			strStruct=""								// flags an error
//		endif
//	endif
//	if (strlen(strStruct)<1)							// nothing found, just use some reasonable looking default values
//		geo.xalfd = 0.15605
//		geo.xbetd = 0.19635
//		geo.dd = 40.0143
//		geo.ddOffset = 0.0143
//		geo.NxCCD = 2084
//		geo.NyCCD = 2084
//		geo.dpsx = 50.01
//		geo.dpsy = 50.168
//		geo.xcent = 1089.7
//		geo.ycent = 942.901							// yc of Si calibration standard
//		geo.xbet = 0.30943
//		geo.xgam = 0.75192
//		geo.wire.dia = 52
//		geo.wire.H0 = -5305.9
//		geo.wire.Hyc = -4889.03
//		geo.wire.F = 3830.
//	endif
//	GeometryUpdateCalc(geo)
//End



//=======================================================================================
//=======================================================================================
//=======================================================================================
//=============================== Start of Q distributions, No Depth ================================

//	Process many images all at the same depth, but in an array of x-y (actually X-H) positions.  There is no wire scan associated with this,
//	and it is assumed that the sample is thin so that the yc in geo is correct for all images.
Function Fill_Q_Positions(d0,pathName,namePart,range,mask)	// does not assume depth
	Variable d0			// d-spacing of the strain=0 material (nm)
	String pathName	// either name of path to images, or the full expliiciit path, i.e. "Macintosh HD:Users:tischler:data:cal:recon:"
	String namePart	// the first part of file name, something like  "EW5_"
	String range		// range of file indicies to use
	Wave mask			// optional mask to limit the pixels that get processed (use pixel when mask true)

	Variable printIt=0
	if (!((d0>0)) && exists("d0")!=2)						// there is no d0, check with the user
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
		elseif (V_flag==3)									// chose "Cancel"
			return 1											// quitting
		endif
	endif
	if (!((d0>0)))											// invalid d0, check for a local value
		d0 = NumVarOrDefault("d0",NaN)
		printIt = 1
	endif

	String str
	PathInfo $pathName
	if (!V_flag || strlen(namePart)<1)						// path does not exist or no namePart, ask user
		String pathPart
		str = ExamineImageSet#requestFileRoot(pathName,1)
		pathPart = StringFromList(0,str)
		namePart = StringFromList(1,str)
		if (strlen(pathPart)<1 || strlen(namePart)<1)
			return 1											// invalid inputs
		endif
		if (!stringmatch(pathPart,S_path))					// path was changed
			if (stringmatch(pathName,"imagePath"))
				NewPath/O/M="path to reconstructed spe files" imagePath pathPart	// for imagePath, automatically reassign
			else
				NewPath/M="path to reconstructed spe files" $pathName pathPart	// for other names, ask
			endif
		endif
		printIt = 1
	endif
	PathInfo $pathName
	if (strlen(S_path)<1 || strlen(namePart)<1)
		return 1												// invalid inputs
	endif
	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
		printf "using data from files starting with '%s'\r",S_path+namePart
	endif

	String list = WaveListClass("imageMask","*","DIMS:2,BYTE:1")	// setup to use optional mask
	String maskName = ""
	if (WaveExists(mask))
		maskName = GetWavesDataFolder(mask,2)
	elseif (strlen(list))
		maskName = ""
		Prompt maskName, "mask to use with image",popup,"_none_;"+list
		DoPrompt "mask wave",maskName
		if (V_flag)
			return 1
		endif
		maskName = SelectString(stringmatch(maskName,"_none_"),maskName,"")
		Wave mask = $maskName								// do not check if wave exists, that a valid option
		printIt = 1
	endif
	maskName = SelectString(WaveExists(mask),"$\"\"",maskName)

	Display/W=(49,105,687,231)/K=1					// put up temp window to view staus of processing
	String gName = S_name
	ModifyGraph gfSize=18
	TextBox/N=text0/F=0/S=3/A=LC/X=2.51/Y=4.76 "Determining the index range to scan"  ;  DoUpdate

	if (ItemsInRange(range)<1) 								// if range is empty, get the full range from the directory
		range = ExamineImageSet#get_FilesIndexRange(pathName,namePart)
		Prompt range,"range of spe file numbers to use"
		DoPrompt "range",range
		if (V_flag)
			return 1
		endif
		if (ItemsInRange(range)<1)
			if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
				printf "range of images indicies is '%s'\r",range[0,370]	// cannot print lines longer than 400 chars
			endif
		endif
		printIt = 1
	endif

	String fileRoot											// complete path up to and including the underscore
	PathInfo $pathName
	if (V_flag)												// if path exists, reset pathName to the explicit path, otherwise assume it is the explicit path
		fileRoot = S_path
	elseif (strsearch(pathName,":",Inf,1) == strlen(pathName)-1)
		fileRoot = pathName+":"								// ensure terminating colon
	endif
	fileRoot += namePart
	if (printIt)
		sprintf str,"Fill_Q_Positions(%g,\"%s\",\"%s\",\"%s\",%s)",d0,pathName,namePart,range,maskName
		print str[0,390]
	endif
	if (WaveExists(mask))
		if (sum(mask)==0)
			DoAlert 0, "You picked a mask that is all zero, stopping"
			DoWindow/K $gName								// done with status window
			return 1
		endif
	endif
	if (ItemsInRange(range)<1) 								// if range is empty, get the full range from the directory
		DoAlert 0, "range is empty"
		DoWindow/K $gName									// done with status window
		return 1
	endif

	Variable timer0 = startMSTimer							// start timer for entire prooess
	Variable timer = startMSTimer							// start timer for initial processing
	STRUCT microGeometry geo
	FillGeometryStructDefault(geo)
	Variable useDistortion = NumVarOrDefault("root:Packages:geometry:useDistortion",USE_DISTORTION_DEFAULT)
 	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
		printf "processing with distotion "+SelectString(useDistortion,"OFF","on")
		print "   "+SelectString(WaveExists(mask),"and no mask","and using  '"+NameOfWave(mask)+"'  for the mask")
	endif

	// open the first file in the range, to get info about the images
	TextBox/C/N=text0/W=$gName  "setting up"  ;  DoUpdate
	String name
	sprintf name, "%s%d%s",fileRoot,str2num(range),IMAGE_FILE_EXT
	Wave image = $(LoadWinViewFile(name))				// load first image in range
	if (!WaveExists(image))
		printf "could not load very first image named '%s'\r",name
		DoAlert 0,"could not load very first image"
		DoWindow/K $gName									// done with status window
		timer=stopMSTimer(timer)  ;  timer0=stopMSTimer(timer0)
		return 1
	endif
	Variable Ni,Nj, Npixels = numpnts(image)				// number of pixels in one image
	Ni = DimSize(image,0)
	Nj = DimSize(image,1)
	String wnote = note(image)
	KillWaves/Z image
	Variable ddLocal=NumberByKey("CCDy", wnote,"=") + geo.ddOffset
	geo.dd = !(ddLocal>0) ? geo.dd : ddLocal					// use dd from the CCDy if available
	Variable keV, startx, endx, starty, endy, groupx, groupy
	keV = NumberByKey("keV", wnote,"=")					// energy for this image
	startx = NumberByKey("startx", wnote,"=");	endx = NumberByKey("endx", wnote,"=");	groupx = NumberByKey("groupx", wnote,"=")
	starty = NumberByKey("starty", wnote,"=");	endy = NumberByKey("endy", wnote,"=");	groupy = NumberByKey("groupy", wnote,"=")
	if (numtype(startx+endx+starty+endy+groupx+groupy))
		DoAlert 0,"could not get ROI from wave note of image '"+name+"'"
		DoWindow/K $gName									// done with status window
		timer=stopMSTimer(timer)  ;  timer0=stopMSTimer(timer0)
		return 1
	elseif (numtype(keV))
		DoAlert 0,"invalid keV in image '"+name+"'"
		DoWindow/K $gName									// done with status window
		timer=stopMSTimer(timer)  ;  timer0=stopMSTimer(timer0)
		return 1
	endif

	Variable Q0 = 2*PI/d0									// Q of unstrained material (1/nm)
	Variable N = ItemsInRange(range)						// number of images to be processed

	// read header from each of the images, and store it to figure out what was done.
	Make/N=(N)/O/D X_FillQvsPositions,H_FillQvsPositions, keV_FillQvsPositions	// hold the sample positions and energies
	X_FillQvsPositions = NaN
	H_FillQvsPositions = NaN
	keV_FillQvsPositions = NaN
	Make/N=(N)/O/U/I m_Fill_QHistAt1Depth
	m_Fill_QHistAt1Depth = 0
	Variable seconds
	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
		seconds = Npixels*N*secPerPixelFixed
		if (useDistortion)
			seconds += Npixels*secPerPixelDistort
		endif
		seconds *= numtype(cpuFrequency()) ? 1 : 2e9/cpuFrequency()
		printf "about to examine %d images, which should take %s,  and end at %s\r",N,Secs2Time(seconds,5,1),Secs2Time(DateTime+seconds,1)
	endif
	// for all the N files (go over range), store the energies, X position, H position, and indicies
	TextBox/C/N=text0/W=$gName  "examining "+num2istr(N)+" image headers\r(to determine what to do)"  ;  DoUpdate
	Variable m, i, j
	for (m=str2num(range), i=0; numtype(m)==0; m=NextInRange(range,m), i+=1)		// loop over range
		sprintf name, "%s%d%s",fileRoot,m,IMAGE_FILE_EXT
		wnote=WinViewReadHeader(name)					// wave note to add to file read in
		if (!sameROI(wnote,startx, endx, starty, endy, groupx, groupy))	// skip bad ROI's
			continue
		endif
		X_FillQvsPositions[i] = NumberByKey("X1", wnote,"=")		// list of X positions
		H_FillQvsPositions[i] = NumberByKey("H1", wnote,"=")		// list of H positions
		keV_FillQvsPositions[i] = NumberByKey("keV", wnote,"=")	// list of energies
		m_Fill_QHistAt1Depth[i] = m									// file id number for each image
	endfor
	seconds = stopMSTimer(timer)/1e6
	if (ItemsInList(GetRTStackInfo(0))<2 && seconds>0.2)
		printf "first pass to examine headers took %s\r",Secs2Time(seconds,5,3)
	endif

	timer = startMSTimer											// start a timer on bulk of processing
	Variable dkeV,NkeV, dx,Nx, dy,Ny, off
	FindScalingFromVec(keV_FillQvsPositions,2e-4,off,dkeV,NkeV)// get step size and number of points for the energy scan
	FindScalingFromVec(X_FillQvsPositions,0.1,off,dx,Nx)		// get step size and number of points for the X scan
	FindScalingFromVec(H_FillQvsPositions,0.1,off,dy,Ny)		// get step size and number of points for the H scan
	// correct the scan ranges as necessary
	dkeV = abs(dkeV)
	dx = abs(dx)
	dy = abs(dy)
	Prompt dkeV,"energy step size (keV)"
	Prompt NkeV,"# of points (NOT intervals) in 1 energy scan"
	Prompt dx,"step size in sample X motion (micron)"
	Prompt Nx,"# of points (NOT intervals) in 1 sample X scan"
	Prompt dy,"step size in sample H motion (micron)"
	Prompt Ny,"# of points (NOT intervals) in 1 sample H scan"
	DoPrompt "scan sizes",dkeV,NkeV,dx,Nx,dy,Ny
	if (V_flag)
		return 1
	endif

	WaveStats/M=1/Q X_FillQvsPositions
	Variable X0 = V_min
	WaveStats/M=1/Q H_FillQvsPositions
	Variable H0 = V_min
	WaveStats/M=1/Q keV_FillQvsPositions
	Variable ikeVlo=V_minloc, ikeVhi=V_maxloc			// save location of min and max energies
	if (V_numNans)
		DoAlert 0, "There were "+num2istr(V_numNans)+" bad images found, they will be skipped"
	endif

	// note that Q range only depends upon image size and energy range, not on X or H position
	Variable thetaLo,thetaHi, Qmin, Qmax, dQ, NQ			// get range of Q
	sprintf name, "%s%d%s",fileRoot,m_Fill_QHistAt1Depth[ikeVlo],IMAGE_FILE_EXT	// load one image to get its size & Q range
	Wave image = $(LoadWinViewFile(name))				// load image
	thetaLo = real(thetaRange(geo,image))					// min theta on this image
	KillWaves/Z image
	keV = keV_FillQvsPositions[ikeVlo]
	Qmin = 4*PI*sin(thetaLo)*keV/hc						// min Q (1/nm)
	sprintf name, "%s%d%s",fileRoot,m_Fill_QHistAt1Depth[ikeVhi],IMAGE_FILE_EXT
	Wave image = $(LoadWinViewFile(name))				// load image
	thetaHi = imag(thetaRange(geo,image))					// max theta on this image
	KillWaves/Z image
	keV = keV_FillQvsPositions[ikeVhi]
	Qmax = 4*PI*sin(thetaHi)*keV/hc						// max Q (1/nm)

	// determine dQ (1/nm), the Q resolution to use.  Base it on the distance between two adjacent pixels
	Variable px,py											// the current pixel to analyze (unbinned full chip pixels)
	py = round((starty+endy)/2)							// approximate full chip pixel position in center or the image
	px = round((startx+endx)/2)
	dQ = 4*PI*abs(sin(pixel2q(geo,px,py,$""))-sin(pixel2q(geo,px,py+groupy,$"")))*keV/hc
//		dQ /= 1.5
	NQ = round((Qmax-Qmin)/dQ) + 1						// number of Qs

	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
		Variable deV = 1000*(keV_FillQvsPositions[ikeVhi]-keV_FillQvsPositions[ikeVlo])/(NkeV-1)
		if (Nx>1)
			printf "X range = [0, %.2f] (µm), dX=%.2f,  %d points\r",(Nx-1)*dx, dx,Nx
		else
			printf "only one X value\r"
		endif
		if (Ny>1)
			printf "H range = [0, %.2f] (µm), dH=%.2f,  %d points\r",(Ny-1)*dy, dy,Ny
		else
			printf "only one H value\r"
		endif
		printf "E range = [%g, %g] (keV),  ∆E=%.2g eV\r", keV_FillQvsPositions[ikeVlo], keV_FillQvsPositions[ikeVhi],deV
		printf "Q range = [%g, %g] (1/nm),  ∆Q=%.2g(1/nm)       ",Qmin, Qmax,dQ
		printf "theta range = [%g, %g]°\r",thetaLo*180/PI,thetaHi*180/PI
		if (d0>0)
			printf "   d0 = %g (nm),   Q0 = %g (1/nm)\r",d0,Q0
		endif
	endif

	Make/N=(Nx,Ny,NQ)/O/D Q_Positions					// array to save Q scans as function of position (X and H)
	Make/N=(Nx,Ny,NQ)/O/D Q_PositionsNorm			// holds number of pixels contributing to each element in Q_Positions[], use to normalize
	Q_Positions = 0
	Q_PositionsNorm = 0
	Make/N=(NQ)/O/D Qhist									// array to hold Q's from one image
	Make/N=(NQ)/O/D QhistNorm							// number of pixels used for each Qhist, use to normalize
	SetScale/P x X0,dx,"µm", Q_Positions,Q_PositionsNorm
	SetScale/P y H0,dy,"µm", Q_Positions,Q_PositionsNorm
	SetScale/I z Qmin,Qmax,"1/nm", Q_Positions, Q_PositionsNorm
	SetScale/I x Qmin,Qmax,"1/nm", Qhist

	seconds = stopMSTimer(timer)/1e6
	if (ItemsInList(GetRTStackInfo(0))<2 && seconds>0.2)
		printf "setting up arrays took %s\r",Secs2Time(seconds,5,3)
	endif
	// done with the setup part, now actually compute something

	if (useDistortion)											// if using the distortion, precompute for all images  here
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
		if (!distortionOK)										// existing map (if it exists) is not OK, so make a new one
			timer = startMSTimer							// start a timer on creation of local distortion map
			KIllWaves/Z root:Packages:geometry:tempCachedDistortionMap	// ensure that this is gone!
			Make/N=(Ni,Nj,2)/O root:Packages:geometry:tempCachedDistortionMapTemp
			Wave distortionMap = root:Packages:geometry:tempCachedDistortionMapTemp
			SetScale/P x startx,groupx,"", distortionMap	// scale distortionMap[][] to full chip pixels
			SetScale/P y starty,groupy,"", distortionMap
			Variable/C dxy
			Wave xymap = root:Packages:geometry:xymap
			TextBox/C/N=text0/W=$gName  "making local copy of the distortion map"  ;  DoUpdate
			for (j=0;j<Nj;j+=1)
				for (i=0;i<Ni;i+=1)
					py = starty + j*groupy
					px = startx + i*groupx
					dxy = peakcorrection2(xymap,px,py)	// returns cmplx(dx,dy)
					distortionMap[i][j][0] = real(dxy)		// accumulate all the distortions that I will need
					distortionMap[i][j][1] = imag(dxy)
				endfor
			endfor
			Rename root:Packages:geometry:tempCachedDistortionMapTemp, tempCachedDistortionMap
			seconds = stopMSTimer(timer)/1e6
			if (ItemsInList(GetRTStackInfo(0))<2 && seconds>0.2)
				printf "creating local copy of distortion map took %s\r",Secs2Time(seconds,5,3)
			endif
		endif
		Note/K DistortionMap, "use=1"
	endif														// distortion map now ready

	Wave sinTheta = $(MakeSinThetaArray(name,geo))	// make an array the same size as an image, but filled with sin(theta) for this energy
	Redimension/N=(Ni*Nj) sinTheta
	Make/N=(Ni*Nj)/I/O indexWaveQ
	indexWaveQ = p
	Sort sinTheta, sinTheta,indexWaveQ						// sort so indexWaveQ[0] is index to lowest sin(theta), indexWaveQ[inf] is greatest

	print "starting bulk of processing"
	timer = startMSTimer									// start a timer on bulk of processing
	Variable sec3=0,timer3
	// for all the N files (go over range), compute Qhist for each image
	for (m=str2num(range); numtype(m)==0; m=NextInRange(range,m))	// loop over range, all the images
		if (mod(m,10)==0)
			TextBox/C/N=text0/W=$gName  "processing image "+num2istr(m)+" (of "+num2istr(N)+")"  ;  DoUpdate
		endif
		sprintf name, "%s%d%s",fileRoot,m,IMAGE_FILE_EXT

		// accumulate the Q histogram for one image into Qhist
		Qhist = 0												// needed because FillQhist1image() accumulates into Qhist
		QhistNorm = 0
		wnote = FillQhist1image(name,sinTheta,indexWaveQ,Qhist,QhistNorm,mask)	// fill Qhist from one file
		timer3=startMSTimer
		if (NumberByKey("V_min", wnote,"=")==0  && NumberByKey("V_max", wnote,"=")==0)
			sec3 += stopMSTimer(timer3)/1e6
			continue											// no intensity here, so continue
		endif
		i = round((NumberByKey("X1", wnote,"=")-X0)/abs(dx))
		j = round((NumberByKey("H1", wnote,"=")-H0)/abs(dy))
		i = numtype(i) ? 0 : i
		j = numtype(j) ? 0 : j
		Q_Positions[i][j][] += Qhist[r]						// add Intensity from this image
		Q_PositionsNorm[i][j][] += QhistNorm[r]			// accumulate number of pixels contibuting (not pixel intensity)
		sec3 += stopMSTimer(timer3)/1e6
		if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
			DoUpdate
		endif
	endfor
	Note/K DistortionMap, "use=0"
	Q_Positions /= Q_PositionsNorm						// do the normalization
	Q_Positions = numtype(Q_Positions[p][q]) ? NaN : Q_Positions[p][q]
	seconds = stopMSTimer(timer)/1e6
	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
		printf "\r  processing all %d images took %s	(%.3g µs/pixel)\r",N,Secs2Time(seconds,5,1),1e6*seconds/(N*Npixels)
		printf "		the accumulation/assignment part took %s\r",Secs2Time(sec3,5,2)
	endif

	timer = startMSTimer									// start a timer on final stuff

	// create the 2d wave to display
	Make/N=(Nx,Ny)/O Q_Positions_dQ	=NaN				// holds the image to plot, delta peak position (1/nm)
	Make/N=(Nx,Ny)/O Q_Positions_Intens	=NaN			// holds the image to plot, max intensity
	Make/N=(Nx,Ny,3)/O/W/U Q_Positions_RGB=0		// holds the RGB image to plot
	Q0 = !(Q0>0) ? 0 : Q0
	for (j=0;j<Ny;j+=1)
		for (i=0;i<Nx;i+=1)
			if (0)												// get position of peak from the max or fit?
				Qhist = Q_Positions[i][j][p]
				WaveStats/Q/M=1 Qhist
				Q_Positions_dQ[i][j] = V_maxloc-Q0		// set delta Q from position of max (no fitting)
			else
				Wave Qhist = $QhistFromQposiitions(Q_Positions,i,j,NaN)		// from a Q_posiitions type array, re-make Qhist, and fit it
				list = fitOneQhist(Qhist,quiet=1)
				Q_Positions_dQ[i][j] = NumberByKey("Qcenter",list,"=")-Q0	// delta Q
				WaveStats/Q/M=1 Qhist
				Q_Positions_Intens[i][j] = V_max			// max intensity
			endif
		endfor
	endfor
	ImageStats/M=1 Q_Positions_dQ
	Variable Qrange = max(abs(V_min),abs(V_max))		// used to scale the color of _Qpk waves, a symmetric Q range
	Variable imax, jmax
	if (abs(V_max)>abs(V_min))
		imax=V_maxRowLoc; jmax=V_maxColLoc
	else
		imax=V_minRowLoc; jmax=V_minColLoc
	endif
	WaveStats/M=1/Q Q_Positions_dQ
	Variable ImaxPnt = V_maxloc
	Variable X1 = V_maxRowLoc*DimDelta(Q_Positions,0) + DimOffset(Q_Positions,0)
	Variable H1 = V_maxColLoc*DimDelta(Q_Positions,1) + DimOffset(Q_Positions,1)
	CopyScales/I Q_Positions, Q_Positions_dQ,Q_Positions_Intens,Q_Positions_RGB
	SetScale d -Qrange,Qrange,"1/nm", Q_Positions_dQ,Q_Positions_RGB
	WaveStats/M=1/Q Q_Positions_Intens
	Variable maxIntens = V_max

	// set the RGB wave based on dQ and Intensity
	Variable cMax = 65535
	Q_Positions_RGB = 0
	Variable red,green,blue, maxc, intens
	for (j=0;j<Ny;j+=1)
		for (i=0;i<Nx;i+=1)
			blue = (Q_Positions_dQ[i][j] + Qrange)/(2*Qrange)	// set the hue
			red = 1-blue
			green= min(red,blue)
			maxC = max(max(red,green),blue)						// max RGB value
			intens = Q_Positions_Intens[i][j]						// intensity of this point
			red = cMax * red/maxC * intens/maxIntens				// set to saturated values scaled down by intensity
			green = cMax * green/maxC * intens/maxIntens
			blue = cMax * blue/maxC * intens/maxIntens
			Q_Positions_RGB[i][j][0] = limit(red,0,cMax)		// red
			Q_Positions_RGB[i][j][1] = limit(green,0,cMax)		// green
			Q_Positions_RGB[i][j][2] = limit(blue,0,cMax)		// blue
		endfor
	endfor

//	String wnote=""
	wnote = ReplaceStringByKey("sourceWave",wnote,GetWavesDataFolder(Q_Positions,2),"=")
	wnote = ReplaceStringByKey("waveClass",wnote,"Q_Positions","=")
	Note/K Q_Positions_Intens, wnote
	wnote = ReplaceNumberByKey("minColor",wnote,-Qrange,"=")		// max and min values, used for the color scaling
	wnote = ReplaceNumberByKey("maxColor",wnote,Qrange,"=")
	Note/K Q_Positions_dQ, wnote
	Note/K Q_Positions_RGB, wnote
	Wave Qhist = $QhistFromQposiitions(Q_Positions,imax,jmax,NaN)	// from a Q_posiitions type array, re-make Qhist, and fit it
	fitOneQhist(Qhist)

	if (numtype(Q0)==0 && Q0>0)
		wnote = ReplaceNumberByKey("Q0",wnote,Q0,"=")
	endif
	String title = StrVarOrDefault("title","")
	if (strlen(title))
		title = ReplaceString("=",title,"_")									// in the wave note the string cannot have "=" or ";"
		title = ReplaceString(";",title,"_")
		wnote = ReplaceStringByKey("title",wnote,title,"=")
	endif
	if (WaveExists(mask))
		wnote = ReplaceStringByKey("maskWave",wnote,GetWavesDataFolder(mask,2),"=")
	endif
	wnote = ReplaceNumberByKey("ImaxPnt",wnote,ImaxPnt,"=")
	wnote = ReplaceNumberByKey("XofMax",wnote,X1,"=")				// X with the max point (µm)
	wnote = ReplaceNumberByKey("HofMax",wnote,H1,"=")				// H with the max point (µm)
	wnote = ReplaceStringByKey("labelX",wnote,"X","=")
	wnote = ReplaceStringByKey("labelY",wnote,"Y","=")
	wnote = ReplaceStringByKey("labelZ",wnote,"Q","=")
	wnote = ReplaceStringByKey("labelValue",wnote,"Intensity","=")
	wnote = ReplaceNumberByKey("geo_NxCCD",wnote,geo.NxCCD,"=")	// add geometry to the wave note
	wnote = ReplaceNumberByKey("geo_NyCCD",wnote,geo.NyCCD,"=")
	wnote = ReplaceNumberByKey("geo_dpsx",wnote,geo.dpsx,"=")
	wnote = ReplaceNumberByKey("geo_dpsy",wnote,geo.dpsy,"=")
	wnote = ReplaceNumberByKey("geo_xcent",wnote,geo.xcent,"=")
	wnote = ReplaceNumberByKey("geo_ycent",wnote,geo.ycent,"=")
	wnote = ReplaceNumberByKey("geo_ddOffset",wnote,geo.ddOffset,"=")
	wnote = ReplaceNumberByKey("geo_dd",wnote,geo.dd,"=")
	wnote = ReplaceNumberByKey("geo_xbet",wnote,geo.xbet,"=")
	wnote = ReplaceNumberByKey("geo_xgam",wnote,geo.xgam,"=")
	wnote = ReplaceNumberByKey("geo_xalfd",wnote,geo.xalfd,"=")
	wnote = ReplaceNumberByKey("geo_xbetd",wnote,geo.xbetd,"=")
	wnote = ReplaceNumberByKey("geo_wire_F",wnote,geo.wire.F,"=")
	wnote = ReplaceNumberByKey("geo_wire_H0",wnote,geo.wire.H0,"=")
	wnote = ReplaceNumberByKey("geo_wire_Hyc",wnote,geo.wire.Hyc,"=")
	wnote = ReplaceNumberByKey("geo_wire_X",wnote,geo.wire.X,"=")
	wnote = ReplaceNumberByKey("geo_wire_dia",wnote,geo.wire.dia,"=")
	Note/K Q_Positions, ReplaceStringByKey("waveClass",wnote,"QdistAtPositions","=")

	wnote = ReplaceStringByKey("sourceWave",wnote,GetWavesDataFolder(Q_Positions,2),"=")
	Note/K Q_Positions_Intens, wnote
	wnote = ReplaceNumberByKey("minColor",wnote,-Qrange,"=")		// max and min values, used for the color scaling
	wnote = ReplaceNumberByKey("maxColor",wnote,Qrange,"=")
	Note/K Q_Positions_dQ, wnote
	Note/K Q_Positions_RGB, wnote

	// put up or refresh the plot of Qhist
	Wave Qhist = $QhistFromQposiitions(Q_Positions,imax,jmax,NaN)	// re-make Qhist, and fit it
	list = fitOneQhist(Qhist)
	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
		// list = note(Qhist)
		printf "at peak, Q = %.3f(1/nm),   fwhm = %.2g(1/nm)",NumberByKey("Qcenter", list,"="),NumberByKey("Qfwhm", list,"=")
		Variable strain = NumberByKey("strain", list,"=")
		if (numtype(strain)==0)
			printf ",   ∆Q = %.3f,   strain = %.1e",NumberByKey("∆Q", list,"="),strain
		endif
		printf "\r"
////		MakeGraph_Q_Positions(Q_Positions_dQ)
//		MakeGraph_Q_Positions(Q_Positions_RGB)
//		MakeGraph_Qhist(Qhist)
		String win = MakeGraph_Qhist(Qhist)
		i = MoveWinToTopRight(win,-1,-1)
		win = MakeGraph_Q_Positions(Q_Positions_RGB)
		MoveWinToTopRight("Graph_root_Q_Positions_RGB",-1,i-2)
	endif

	// done with processing, clean up
	DoWindow/K $gName														// done with status window
	KillWaves/Z X_FillQvsPositions,H_FillQvsPositions, keV_FillQvsPositions, m_Fill_QHistAt1Depth, sinTheta, indexWaveQ, Q_PositionsNorm, QhistNorm
	seconds = stopMSTimer(timer)/1e6
	if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
		printf "  final stuff took %s\r",Secs2Time(seconds,5,2)
	endif
	seconds = stopMSTimer(timer0)/1e6
	if (ItemsInList(GetRTStackInfo(0))<2 && seconds>0.2)
		printf "entire process took %s\r",Secs2Time(seconds,5,3)
	endif
	beep
	return 0
End
//
Static Function MoveWinToTopRight(win,tt,rr)			// returns left edge of window
	String win
	Variable tt,rr											// top right corner of window (use NaN for screen edge)
	if (strlen(win))
		win = StringFromList(0,WinList("*",";", "WIN:4183"))
	endif
	String str = StringByKey("SCREEN1", IgorInfo(0))		// RECT=left,top,right,bottom
	Variable i = strsearch(str,"RECT=",0)+5
	Variable left,top,right,bottom
	sscanf str[i,inf], "%d,%d,%d,%d",left,top,right,bottom	// get screen size
	if (V_flag!=4)
		return NaN
	endif
	right = !(rr>0) ? right : rr								// desired right side
	top = max(!(tt>52) ? top : tt,45)						// desired top
	GetWindow $win, wsize									// get current window size
	left = max(0,right-V_right+V_left)					// try to preserve width
	bottom = top+V_bottom-V_top
	MoveWindow/W=$win left,top,right,bottom
	return left
End



Function/S QhistFromQposiitions(source,i,j,k)				// from a Q_posiitions type array, make the Qhist (a 1d array)
	Wave source
	Variable i,j,k												// index into source to id where to get the Qhist
	Variable dimLast = WaveDims(source)-1				// dimension of Q dimension (always the last one)
	String Qname = GetWavesDataFolder(source,1)+"Qhist"
	Make/N=(DimSize(source,dimLast))/O/D $Qname		// re-make array to hold Q's from one image
	Wave Qhist = $Qname
	SetScale/P x DimOffset(source,dimLast),DimDelta(source,dimLast),WaveUnits(source,dimLast), Qhist
	SetScale d 0,0,WaveUnits(source,-1), Qhist

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
	String wnote = note(source)
	wnote = ReplaceStringByKey("waveClass",wnote,"Qhistogram","=")
	wnote = ReplaceStringByKey("sourceWave",wnote,GetWavesDataFolder(source,2),"=")
	wnote = ReplaceNumberByKey("sourceWave_row",wnote,i,"=")

	String labelX=""
	if (dimLast==0)
		wnote = RemoveByKey("sourceWave_row",wnote,"=")
		wnote = RemoveByKey("sourceWave_col",wnote,"=")
		wnote = RemoveByKey("sourceWave_layer",wnote,"=")
		wnote = RemoveByKey("sourceWave_chunk",wnote,"=")
		labelX = StringByKey("labelX",wnote,"=")
	elseif (dimLast==1)
		wnote = RemoveByKey("sourceWave_col",wnote,"=")
		wnote = RemoveByKey("sourceWave_layer",wnote,"=")
		wnote = RemoveByKey("sourceWave_chunk",wnote,"=")
		labelX = StringByKey("labelY",wnote,"=")
	elseif (dimLast==2)
		wnote = ReplaceNumberByKey("sourceWave_col",wnote,j,"=")
		wnote = RemoveByKey("sourceWave_layer",wnote,"=")
		wnote = RemoveByKey("sourceWave_chunk",wnote,"=")
		labelX = StringByKey("labelZ",wnote,"=")
	elseif (dimLast==3)
		wnote = ReplaceNumberByKey("sourceWave_col",wnote,j,"=")
		wnote = ReplaceNumberByKey("sourceWave_layer",wnote,NumberByKey("plane", wnote,"="),"=")
		wnote = RemoveByKey("sourceWave_chunk",wnote,"=")
		labelX = StringByKey("labelW",wnote,"=")
	endif
	if (strsearch(";"+wnote,";labelX=",0)>=0)
		wnote = ReplaceStringByKey("source_labelX",wnote,StringByKey("labelX",wnote,"="),"=")
		wnote = RemoveByKey("labelX",wnote,"=")
	endif
	if (strsearch(";"+wnote,";labelY=",0)>=0)
		wnote = ReplaceStringByKey("source_labelY",wnote,StringByKey("labelY",wnote,"="),"=")
		wnote = RemoveByKey("labelY",wnote,"=")
	endif
	if (strsearch(";"+wnote,";labelZ=",0)>=0)
		wnote = ReplaceStringByKey("source_labelZ",wnote,StringByKey("labelZ",wnote,"="),"=")
		wnote = RemoveByKey("labelZ",wnote,"=")
	endif
	if (strsearch(";"+wnote,";labelW=",0)>=0)
		wnote = ReplaceStringByKey("source_labelW",wnote,StringByKey("labelW",wnote,"="),"=")
		wnote = RemoveByKey("labelW",wnote,"=")
	endif
	labelX = SelectString(strlen(labelX),"Q",labelX)		// default value is "Q"
	wnote = ReplaceStringByKey("labelX",wnote,labelX,"=")
	String labelValue = StringByKey("labelValue",wnote,"=")
	labelValue = SelectString(strlen(labelValue),"Intensity",labelValue)// default value is "Intensity"
	wnote = ReplaceStringByKey("labelValue",wnote,labelValue,"=")

	Note/K Qhist, wnote
	// fitOneQhist(Qhist)
	return GetWavesDataFolder(Qhist,2)
End


Function/S fitOneQhist(Qhist, [quiet])
	Wave Qhist						// a 1-d Q distribution
	Variable quiet						// only fit, NO fit_Qhist wave and NO ouput, NO graph updating (goes much faster too)
	if (!WaveExists(Qhist))
		return ""
	endif
	quiet = ParamIsDefault(quiet) ? 0 : quiet
	String fldrSav= GetDataFolder(1)
	SetDataFolder  $GetWavesDataFolder(Qhist,1)		// data folder with the data

//		Make/T/O T_Constraints={"K0 > 0","K3 > 0"}	// ensure width>0, and base line positive
//		WaveStats/Q Qhist
//		T_Constraints[0] = "K0 > "+num2str(min(0,V_min))
//		CurveFit/Q lor Qhist/D /C=T_Constraints 
//		KillWaves/Z T_Constraints
	Variable V_FitOptions=4, V_FitError=0			//suppress window & ignore errors
	Variable lo,hi, fwhm=NaN, Qc=NaN
	if(!getPercentOfPeak(Qhist,0.5,lo,hi))				// find the 50% points
		if (quiet)
			CurveFit/Q/N lor Qhist[lo,hi]
		else
			CurveFit/Q lor Qhist[lo,hi]/D 
		endif
		if (!V_FitError)
			fwhm = 2*sqrt(K3)							// for Lorentzian
			Qc = K2
		endif
	endif
	String list = ReplaceNumberByKey("Qcenter","",Qc,"=")
	list = ReplaceNumberByKey("Qfwhm",list,fwhm,"=")

	String wnote = note(Qhist)
	wnote = ReplaceNumberByKey("Qcenter",wnote,Qc,"=")
	wnote = ReplaceNumberByKey("Qfwhm",wnote,fwhm,"=")
	if (!quiet && (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc")))
		printf "at peak, Q = %.5g(1/nm),   fwhm = %.2g(1/nm),   fwhm/Q = %.1e,      chisq=%.3g\r",Qc,fwhm, fwhm/Qc,V_chisq
	endif
	Variable d0 = NumVarOrDefault("d0",NaN)			// unstrained d-spacing (nm)
	if (numtype(d0)==0)
		Variable Q0 = 2*PI/d0							// Q of unstrained material (1/nm)
		Variable strain = -(Qc-Q0)/Qc					// d is inverse of strain, so a negative
		list = ReplaceNumberByKey("d0",list,d0,"=")
		list = ReplaceNumberByKey("Q0",list,Q0,"=")
		list = ReplaceNumberByKey("strain",list,strain,"=")
		if (!quiet && (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc")))
			printf "   d0 = %g (nm),   Q0 = %g (1/nm),   strain = %.1e\r",d0,Q0,strain
		endif
		list = ReplaceNumberByKey("∆Q",list,Qc-Q0,"=")
		wnote = ReplaceNumberByKey("∆Q",wnote,Qc-Q0,"=")
		wnote = ReplaceNumberByKey("strain",wnote,strain,"=")
	endif
	Note/K Qhist, wnote

	// if Qhist is on a graph, update textbox and line, so find a plot (if one exists, containing Qhist)
	String win = StringFromList(0,FindGraphsWithWave(Qhist))	// find top most graph containing Qhist
	if (!quiet && strlen(win))
		if (ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"EWscanButtonProc"))
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



Function/S MakeGraph_Q_Positions(image)			// display a Q_Positions[][][] array
	Wave image										// usually Q_Positions_dQ
	if (!WaveExists(image))
		DoAlert 0,"no 'Q_Positions' in this data folder"
		return ""
	endif
	String graphName=CleanupName("Graph_"+GetDataFolder(0)+"_"+NameOfWave(image),0)
	String wnote = note(image)
	if (strlen(WinList(graphName,"","WIN:1")))
		DoWindow/F $graphName					// graph is already up, so just bring it to front
	elseif (exists(graphName)==5)	
		Execute graphName+"()"					// a recreation macro exists, run it
	else
		Display/K=1/W=(292,67,728,479)			// nothing exists, create the graph
		DoWindow/C $graphName
		AppendImage image
		if (numtype(NumberByKey("Q0",wnote,"="))==0)
			WaveStats/M=1/Q image
			Variable zrange = max(abs(V_min),abs(V_max))
			ModifyImage $NameOfWave(image) ctab= {-zrange,zrange,RedWhiteBlue,0}
		else
			ModifyImage $NameOfWave(image) ctab= {*,*,Terrain,1}
		endif
		ModifyGraph gfMult=130, tick=0, mirror=1, minor=1, lowTrip=0.001, standoff=0
		ModifyGraph axOffset(left)=-1.5
		String units = WaveUnits(image,1)

		String str = StringByKey("labelX",wnote,"=")
		if (strlen(str))
			Label bottom str+FixUnitsInLabel(WaveUnits(image,0))
		endif
		str = StringByKey("labelY",wnote,"=")
		if (strlen(str))
			Label left str+FixUnitsInLabel(WaveUnits(image,1))
		endif
		str = StringByKey("title",wnote,"=")		// first try to get title from wave note
		if (strlen(str)<1)
			str = StrVarOrDefault("title","")		// nothing in wavenote, try the title global string
		endif
		str = SelectString(strlen(str),"",str+"\r\\Zr067")
		str += GetWavesDataFolder(image,1)
		if (strlen(StringByKey("dateExposed",wnote,"=")))
			str += "\r"+StringByKey("dateExposed",wnote,"=")
		endif
		TextBox/C/N=text0/F=0 str
			Cursor/P/I/L=1/H=1 A $NameOfWave(image) 0,0
		// ShowInfo
	endif
	ImageStats/M=1 image
	Cursor/P/I A $NameOfWave(image) V_maxRowLoc,V_maxColLoc

	Variable Q0 = NumberByKey("Q0",wnote,"=")	// un-strained Q (1/nm)
	if (numtype(Q0)==0 && Q0>0)					// always redo the line showing Q0 if valid value is passed
		SetWindow $graphName userdata(Q0)="Q0="+num2str(Q0)+";"
	endif
	SetAspectToSquarePixels(graphName)

	SetWindow kwTopWin hook(QfromCursor)=QfromCursorHook	// this forces a re-calculation of Qhist when cursor A is moved
//	sprintf str, "value=%s;proc=SetVarProcViewQdistn;ctrlName=QcutMarker;stepSize=%g;loLim=%g;hiLim=%g",fullName,dx,x0,x1
//	SetWindow kwTopWin userdata(QfromCursor)=str
	MakeGraph_Qhist(Qhist)							// also put up the associated plot of Q distributions
	return graphName
End


Function QfromCursorHook(s)					// updates the Q distribution plot by the cursor on a Q_positions image
	STRUCT WMWinHookStruct &s
	if (s.eventCode!=7 || !stringmatch(s.cursorName,"A"))	// only trap cursor A moved events
		return 0								// a zero lets other hooks process this event
	endif

	String win=s.winName
	Wave image = ImageNameToWaveRef(win,StringFromList(0,ImageNameList(win,";")))
	if (!WaveExists(image))
		return 0
	endif
	String wnote = note(image)
	Wave source = $StringByKey("sourceWave", wnote,"=")
	if (!WaveExists(source))
		return 0
	endif
	Wave Qhist = $QhistFromQposiitions(source,s.pointNumber,s.yPointNumber,NaN)
	fitOneQhist(Qhist)
	return 1									// a one stops other hooks from processing too
//	Variable step=1
//	if (s.keycode==29 || s.keycode==30)		// 29 is right,  30 is up	make bigger
//		step = 1
//	elseif (s.keycode==28 || s.keycode==31)	// 29 is leftt,  31 is down	make smaller
//		step = -1
//	else
//		return 0
//	endif
//	String win=s.winName
//	String udata = GetUserData(win,"","QfromCursor")
//	String funcSpec=StringByKey("proc", udata,"=")
//	NVAR value = $StringByKey("value", udata,"=")
///	Variable lo=NumberByKey("loLim", udata,"="), hi=NumberByKey("hiLim", udata,"=")
//	Variable stepSize = NumberByKey("stepSize", udata,"=")
//	if (!NVAR_Exists(value) || exists(funcSpec)!=6 || numtype(stepSize))
//		return 0
//	endif
//	Variable newvalue =  step*stepSize+value
//	if (numtype(lo)<2 && numtype(hi)<2)
//		newvalue = limit(newvalue,lo,hi)
//	endif
//	if (numtype(newvalue))
//		return 0
//	endif
//	value = newvalue
//	FUNCREF SetEsumDepthProc  arrowProc= $funcSpec
//	arrowProc(StringByKey("ctrlName",udata,"="),value,StringByKey("varStr",udata,"="),StringByKey("varName",udata,"="))
//	return 1									// a one stops other hooks from processing too
End
//
//	SetWindow kwTopWin hook(QfromCursor)=QfromCursorHook	// this forces a re-calculation of Qhist when cursor A is moved
//	SetWindow kwTopWin userdata(arrows)="value=root:test:EsumDepth;proc=SetEsumDepthProc;stepSize=1;hiLim=71;loLim=-1;"



// find the step size in a coordinate by examining the values
Static Function FindScalingFromVec(vec,threshold,first,stepSize,dimN)
	Wave vec
	Variable threshold					// a change greater than this is intentional (less is jitter)
	Variable &first						// use 	SetScale/P x first,stepSize,"",waveName
	Variable &stepSize
	Variable &dimN
	threshold = threshold>0 ? threshold : 0.1	// no negative thresholds

	Variable last,i, Ndim=1
	String sortName = UniqueName("vecSort",1,0)

	// first, the step size
	Duplicate/O vec $sortName
	Wave vecSort = $sortName
	Variable N = numpnts(vecSort)
	SetScale/P x 0,1,"", vecSort

	vecSort = vecSort[p+1]-vecSort[p]					// make vecSort the differences
	vecSort = abs(vecSort[p])<threshold ? NaN : vecSort[p]
	Sort vecSort, vecSort

	WaveStats/Q vecSort
	N = V_npnts
	Redimension/N=(N) vecSort						// remove all NaN's (small steps) 
	if (N<1)
		KillWaves/Z vecSort
		first = vec[0]									// no steps, just use first point
		stepSize = 0
		dimN = 1
		return 0
	endif
	stepSize = vecSort[N/2]							// take the median value (avoids problems with average)

	// second, the starting point
	Duplicate/O vec $sortName
	Wave vecSort = $sortName
	SetScale/P x 0,1,"", vecSort
	Sort vecSort, vecSort
	i = floor(BinarySearchInterp(vecSort, vecSort[0]+threshold))
	first = faverage(vecSort,0,i)

	// third, find number of points in "one scan"
	N = numpnts(vec)
	Variable tm=0,tn=0, m=0
	for (i=1;i<N;i+=1)
		if (abs(vec[i]-vec[i-1])<threshold)				// this step is a repeat, do not count it
			continue
		elseif (abs(vec[i]-vec[i-1]-stepSize)<(abs(stepSize)+threshold))	// this step is the right size
			m += 1
		else
			tm += m+1										// sum of number of points in "a scan"
			tn += 1											// number of scans found
			m = 0
		endif
	endfor
	if (m>0)
		tm += m+1											// sum of number of points in "a scan"
		tn += 1												// number of scans found
	endif
	dimN = round(tm/ tn)
	KillWaves/Z vecSort
	return 0
End
//Static Function FindScalingFromVec(vec,threshold,first,stepSize,dimN)
//	Wave vec
//	Variable threshold					// a change greater than this is intentional (less is jitter)
//	Variable &first						// use 	SetScale/P x first,stepSize,"",waveName
//	Variable &stepSize
//	Variable &dimN
//	threshold = threshold>0 ? threshold : 0.1
//
//	Variable last,i, Ndim=1
//	String sortName = UniqueName("vecSort",1,0)
//
//	// first, the step size
//	Duplicate/O vec $sortName
//	Wave vecSort = $sortName
//	Variable N = numpnts(vecSort)
//	SetScale/P x 0,1,"", vecSort
//	Sort vecSort, vecSort
//	vecSort = vecSort[p+1]-vecSort[p]					// make vecSort the differences
//	vecSort = vecSort[p]<threshold ? NaN : vecSort		// set all small steps to NaN
//	Sort vecSort, vecSort								// sort by increasing step size
//	WaveStats/Q vecSort
//	N = V_npnts
//	Redimension/N=(N) vecSort						// remove all NaN's (small steps) 
//	if (N<1)
//		KillWaves/Z vecSort
//		return 1
//	endif
//	i = BinarySearch(vecSort,vecSort[0]+threshold)
//	i = (i==-2) ? N-1 : i
//	stepSize = faverage(vecSort,0,i)					// the most appropriate step size
//
//	// second, the starting point
//	Duplicate/O vec $sortName
//	Wave vecSort = $sortName
//	SetScale/P x 0,1,"", vecSort
//	Sort vecSort, vecSort
//	i = floor(BinarySearchInterp(vecSort, vecSort[0]+threshold))
//	first = faverage(vecSort,0,i)
//
//	// third, find number of points in "one scan"
//	N = numpnts(vec)
//	Variable tm=0,tn=0, m=0
//	for (i=1;i<N;i+=1)
//		if (abs(vec[i]-vec[i-1])<threshold)				// this step is a repeat, do not count it
//			continue
//		elseif (abs(vec[i]-vec[i-1]-stepSize)<abs(stepSize+threshold))	// this step is the right size
//			m += 1
//		else
//			tm += m+1										// sum of number of points in "a scan"
//			tn += 1											// number of scans found
//			m = 0
//		endif
//	endfor
//	if (m>0)
//		tm += m+1											// sum of number of points in "a scan"
//		tn += 1												// number of scans found
//	endif
//	dimN = round(tm/ tn)
//	KillWaves/Z vecSort
//	return 0
//End
//// find the step size in a coordinate by examining the values
//Static Function FindScalingFromVec(vec,threshold,first,stepSize,dimN)
//	Wave vec
//	Variable threshold					// a change greater than this is intentional (less is jitter)
//	Variable &first						// use 	SetScale/P x first,stepSize,"",waveName
//	Variable &stepSize
//	Variable &dimN
//	threshold = threshold>0 ? threshold : 0.1
//
//	Variable last,i, Ndim=1
//	String sortName = UniqueName("vecSort",1,0)
//	Duplicate/O vec $sortName
//	Wave vecSort = $sortName
//	Variable N = numpnts(vecSort)
//	SetScale/P x 0,1,"", vecSort
//	Sort vecSort, vecSort
//
//	// first, get the best guess at the stepSize
//	Variable delta, stepSizeSum, Nsum, currentStep
//	stepSize=Inf
//	for (i=0;i<(N-1);i+=1)
//		delta = vecSort[i+1]-vecSort[i]
//		if (delta<(stepSize-threshold) && delta>threshold)// found a significantly smaler step (but not too small)
//			stepSizeSum = delta
//			Nsum = 1
//			stepSize = delta
//		elseif (abs(delta-stepSize)<threshold)				// found more of this step size, accumulate into sum
//			stepSizeSum += delta
//			Nsum += 1
//			stepSize = stepSizeSum/Nsum
//		endif
//	endfor
//	if (numtype(stepSize)==1)								// no steps, all are the same
//		WaveStats/M=1/Q vec
//		first = V_avg
//		stepSize = 0
//		dimN=1
//		KillWaves/Z vecSort
//		return 1
//	endif
//
//	// second, the starting point
//	i = floor(BinarySearchInterp(vecSort, vecSort[0]+threshold))
//	first = faverage(vecSort,0,i)
//
//	// third, find number of points in "a scan"
//	Variable tm=0,tn=0, m=0
//	for (i=1;i<N;i+=1)
//		if (abs(vec[i]-vec[i-1])<threshold)				// this step is a repeat, do not count it
//			continue
//		elseif (abs(vec[i]-vec[i-1]-stepSize)<abs(stepSize+threshold))	// this step is the right size
//			m += 1
//		else
//			tm += m+1										// sum of number of points in "a scan"
//			tn += 1											// number of scans found
//			m = 0
//		endif
//	endfor
//	if (m>0)
//		tm += m+1											// sum of number of points in "a scan"
//		tn += 1												// number of scans found
//	endif
//	dimN = round(tm/ tn)
//
//	KillWaves/Z vecSort
//	return 0
//End


//Static Function FixAspectRatio(graphName)		// set aspect ratio of graph to match the axes
//	String graphName
//	if (strlen(graphName)<1)
//		graphName = StringFromList(0,WinList("*",";","WIN:1"))
//	endif
//	DoUpdate
//	if (WinType(graphName)!=1)
//		DoAlert 0, "FixAspectRatio(), unable to get size of bottom axis"
//		return 1
//	endif
//	GetAxis /W=$graphName/Q bottom
//	if (V_flag)
//		DoAlert 0, "FixAspectRatio(), unable to get size of bottom axis"
//		return 1
//	endif
//	Variable xRange=abs(V_max-V_min)
//	GetAxis /W=$graphName/Q left
//	if (V_flag)
//		DoAlert 0, "FixAspectRatio(), unable to get size of left axis"
//		return 1
//	endif
//	Variable yRange=abs(V_max-V_min)
//	Variable aspect = yRange/xRange
//	//	printf "xRange=%g, yRange=%g,  aspect=%g\r",xRange,yRange,aspect
//	if (aspect>=1 && !numtype(aspect))
//		ModifyGraph width={Aspect,1/aspect},height=0		// tall & thin
//	elseif (!numtype(aspect))
//		ModifyGraph height={Aspect,aspect},width=0		// short & fat 		
//	else
//		DoAlert 0, "Unable to determie aspect ratio of graph"
//	endif
//End


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
			if (!stringmatch(units,"1/nm") && !stringmatch(units,"1/Å"))
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
	if (!stringmatch(units,"1/nm") && !stringmatch(units,"1/Å"))
		DoAlert 0,"wave units not 1/nm or 1/Å, they are '"+units+"'"
		return 1
	endif
	Variable N=numpnts(w), Qstart=DimOffset(w,0), dQ=DimDelta(w,0), Qc=2*PI/d0
	SetScale/P x -(Qstart-Qc)/Qc,-dQ/Qc,"", w
End


Function/S FixUnitsInLabel(units)
	String units
	if (strlen(units)<1)
		return ""
	endif
	return "  "+SelectString(strsearch(units,"1/",0)==0,"(\\U)","(\\E"+units[2,Inf]+"\\S-1\\M)")
End

//=============================== End of Q distributions, No Depth =================================
//=======================================================================================
//=======================================================================================



//Function/T FillEWscanParametersPanel(strStruct,hostWin,left,top)
//	String strStruct									// optional passed value of xtal structure, this is used if passed
//	String hostWin										// name of home window
//	Variable left, top									// offsets from the left and top
//
//	NewDataFolder/O root:Packages:micro				// ensure that the needed data folders exist
//	NewDataFolder/O root:Packages:micro:EWscan
//
//
//	SetWindow kwTopWin,userdata(EWscanPanelName)=hostWin+"#EWscanPanel"
//	NewPanel/K=1/W=(left,top,left+221,top+365)/HOST=$hostWin
//	ModifyPanel frameStyle=0, frameInset=0
//	RenameWindow #,EwscanPanel
//
//	Button buttonInitEWscan,pos={33,31},size={150,50},title="ready to go"
//
//	return "#EWscanPanel"
//End
Function/T FillEWscanParametersPanel(strStruct,hostWin,left,top)
	String strStruct									// optional passed value of xtal structure, this is used if passed
	String hostWin										// name of home window
	Variable left, top									// offsets from the left and top

	NewDataFolder/O root:Packages:micro				// ensure that the needed data folders exist
	NewDataFolder/O root:Packages:micro:EWscan

	Variable/G root:Packages:micro:EWscan:d0				// d-spacing of unstrained material (nm)
	String/G root:Packages:micro:EWscan:fldr				// new folder name, it will be at the top level
	String/G root:Packages:micro:EWscan:title				// a title to be used on plots

	SetWindow kwTopWin,userdata(EWscanPanelName)=hostWin+"#EWscanPanel"
	NewPanel/K=1/W=(left,top,left+221,top+465)/HOST=$hostWin
	ModifyPanel frameStyle=0, frameInset=0
	RenameWindow #,EwscanPanel

	Button buttonNewEWscan,pos={13,5},size={195,20},proc=EWscanButtonProc,title="Set Up New E-W scan..."
	Button buttonNewEWscan,help={"Create a new data folder for an energy-wire scan, and prompt user for title and d0"}
	Button buttonNewEWscan,userdata(active)="0"

	SetVariable setvarTitle,pos={13,30},size={190,18},title="title",fSize=12,value= root:Packages:micro:EWscan:title
	SetVariable setvarTitle,help={"title for the plots"},font="Lucida Grande",noedit=1,frame=0
	SetVariable setvard0,pos={45,50},size={128,18},title="d0 (nm)",noedit=1,frame=0
	SetVariable setvard0,fSize=12,format="%.7f",limits={0,inf,0},value= root:Packages:micro:EWscan:d0
	SetVariable setvard0,help={"referene d-spacing (nm)"},font="Lucida Grande"
	SetVariable setvarfldr,pos={40,70},size={162,18},title="folder name",font="Lucida Grande",fSize=12,noedit=1,frame=0
	SetVariable setvarfldr,help={"name of the igor data folder to make"},value= root:Packages:micro:EWscan:fldr
//	SetVariable setvarTitle,pos={13,30},size={190,18},title="title",fSize=12,value= root:Packages:micro:EWscan:title
//	SetVariable setvarTitle,help={"title for the plots"},font="Lucida Grande",disable=2
//	SetVariable setvard0,pos={45,50},size={128,18},title="d0 (nm)",disable=2
//	SetVariable setvard0,fSize=12,format="%.7f",limits={0,inf,0},value= root:Packages:micro:EWscan:d0
//	SetVariable setvard0,help={"referene d-spacing (nm)"},font="Lucida Grande"
//	SetVariable setvarfldr,pos={40,70},size={162,18},title="folder name",font="Lucida Grande",fSize=12,disable=2
//	SetVariable setvarfldr,help={"name of the igor data folder to make"},value= root:Packages:micro:EWscan:fldr

	Button buttonFillQvsDepth,pos={13,100},size={195,20},proc=EWscanButtonProc,title="Make Q-Depth Surface"
	Button buttonFillQvsDepth,help={"read the depth sorted images, and make a 2-d plot of Q vs depth"}
	Button buttonFitQAt1Depth,pos={53,123},size={155,20},proc=EWscanButtonProc,title="Display & fit Q at 1 depth",fSize=9
	Button buttonFitQAt1Depth,help={"you should not need this, this is usually done by changing the line on the Q-depth plot"}
	Button buttonRePlotQDepth,pos={53,146},size={155,20},proc=EWscanButtonProc,title="Re-Plot Q-Depth Surface",fSize=9
	Button buttonRePlotQDepth,help={"re-plot the Q-depth info read in using Fill_QvsDepth()"}
	Button buttonlayoutQDepth,pos={53,169},size={155,20},proc=EWscanButtonProc,title="Make Layout of Q-Depth",fSize=9
	Button buttonlayoutQDepth,help={"make a nice layout for printing showing the Q-depth surface plot and the plot of a cut through it"}

	Button buttonFillEvsDepth,pos={13,200},size={195,20},proc=EWscanButtonProc,title="Make E-Depth Surface"
	Button buttonFillEvsDepth,help={"read the depth sorted images, and make a 2-d plot of energy vs depth, this command is not used much"}
	Button buttonRePlotEDepthSurf,pos={53,223},size={155,20},disable=2,proc=EWscanButtonProc,title="Re-Plot E-Depth Surface"
	Button buttonRePlotEDepthSurf,help={"re-plot the energy-depth info read in using Fill_EvsDepth()"},fSize=9

	Button buttonPseudoWhite,pos={13,250},size={195,20},proc=EWscanButtonProc,title="Make Pseudo White Image"
	Button buttonPseudoWhite,help={"read all the depth sorted images and make the pseudo white beam plot"}
	Button buttonRePlotEsumDepth,pos={53,273},size={155,20},disable=2,proc=EWscanButtonProc,title="Re-Plot Energy Sums vs Depth"
	Button buttonRePlotEsumDepth,help={"take the pseudo white data already read in and re-plot it"},fSize=9
	Button buttonMovieEsumsDepth,pos={53,296},size={155,20},disable=2,proc=EWscanButtonProc,title="Movie of Energy Sums vs Depth"
	Button buttonMovieEsumsDepth,help={"take the pseudo white data already read in and make the depth movie"},fSize=9

	Button buttonMakeMask,pos={13,325},size={195,20},proc=EWscanButtonProc,title="Make a Mask"
	Button buttonMakeMask,help={"make a mask from an image (often the pseudo white image)"}
	Button buttonClearMask,pos={53,348},size={155,20},disable=2,proc=EWscanButtonProc,title="Clear the Mask"
	Button buttonClearMask,help={"clear the current mask (usually accessed by clicking in marquee)"},fSize=9

	Button buttonFillQNoDepth,pos={13,379},size={195,20},proc=EWscanButtonProc,title="Make Q dist'n at one Depth"
	Button buttonFillQNoDepth,help={"read raw images at many positions (but all from one depth and mult. energies) and get Q disributions"}

	Button buttonPlotQpositionsRGBNoDepth,pos={13,410},size={195,20},disable=2,proc=EWscanButtonProc,title="Display positions of Q dist'n"
	Button buttonPlotQpositionsRGBNoDepth,help={"make RGB plot of postions used to get Q distributions (no depths)"}

//		     MenuItemIfWaveExists("    Graph 1 Q historam","Qhist"),MakeGraph_Qhist()
//		     help={"re-plot the single slice from the Q-depth surface"}
	EnableDisableEWscanControls(hostWin+"#EWscanPanel")
	return "#EWscanPanel"
End


Function EWscanButtonProc(ctrlName) : ButtonControl
	String ctrlName

	String win = GetUserData("","","EWscanPanelName")
	if (stringmatch(ctrlName,"buttonNewEWscan"))
		Variable active = str2num(GetUserData(win,"buttonNewEWscan","active"))
		NVAR d0=root:Packages:micro:EWscan:d0
		if (!active)										// let user enter values and make new data folder
			SetVariable setvarTitle,win=$win,noedit=0,frame=1,disable=0
			SetVariable setvard0,win=$win,noedit=0,frame=1,disable=0
			SetVariable setvarfldr,win=$win,noedit=0,frame=1,disable=0
			Button buttonNewEWscan,win=$win,title="Enter Values, then click here",fColor=(65535,65534,49151),userdata(active)="1"
			d0 = !(d0>0) ? NumVarOrDefault("root:Packages:Lattices:PanelValues:dspace_nm",d0) : d0
		else												// get values and create the new folder
			SVAR fldr=root:Packages:micro:EWscan:fldr
			SVAR title=root:Packages:micro:EWscan:title
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
			Button buttonNewEWscan,win=$win,title="Set Up New E-W scan...",fColor=(0,0,0),userdata(active)="0"
			if (NewEnergyWireScan(fldr,title,d0))		// set up a folder for a new Energy-Wire Scan
				DoAlert 0,"New folder not made"
				fldr = ""
			endif
		endif
	elseif (stringmatch(ctrlName,"buttonFillQvsDepth"))
		printf "•Fill Q vs Depth array\r"
		Fill_QvsDepth(NaN,"imagePath","","","",$"")
	elseif (stringmatch(ctrlName,"buttonFitQAt1Depth") && WaveExists(Q_Depth))
		printf "•MakeGraph_Qhist($\"\")\r"
		MakeGraph_Qhist($"")
//		fitQat1depth(NaN,":")
	elseif (stringmatch(ctrlName,"buttonRePlotQDepth") && WaveExists(Q_Depth))
		printf "•MakeGraph_Q_Depth()\r"
		MakeGraph_Q_Depth()
	elseif (stringmatch(ctrlName,"buttonlayoutQDepth") && strlen(WinList("*_Q_Depth",";","WIN:1")) && strlen(WinList("*_Qhist",";","WIN:1")))
		MakeLayoutQ_depth_hist("")
	elseif (stringmatch(ctrlName,"buttonFillEvsDepth"))
		printf "•Fill_EvsDepth(\"imagePath\",\"\")\r"
		Fill_EvsDepth("imagePath","")
	elseif (stringmatch(ctrlName,"buttonRePlotEDepthSurf") && WaveExists(keV_Depth))
		printf "•MakeGraph_keV_Depth()\r"
		MakeGraph_keV_Depth()
	elseif (stringmatch(ctrlName,"buttonPseudoWhite"))
		printf "•MakePseudoWhiteImages(\"\",\"\")\r"
		MakePseudoWhiteImages("","")
	elseif (stringmatch(ctrlName,"buttonRePlotEsumDepth") && WaveExists(imageEsum))
		printf "•MakeEsumPlot()\r"
		MakeEsumPlot()
	elseif (stringmatch(ctrlName,"buttonMovieEsumsDepth") && WaveExists(imageEsum))
		printf "•MakeEsumMovie()\r"
		MakeEsumMovie()
	elseif (stringmatch(ctrlName,"buttonMakeMask") && strlen(WaveListClass("speImage*;ImageSummed*","*","DIMS:2")))
		printf "•MakeMask($\"\")\r"
		MakeMask($"")
	elseif (stringmatch(ctrlName,"buttonClearMask") && MaskToClearExists())
		ClearMask()
	elseif (stringmatch(ctrlName,"buttonFillQNoDepth"))
		printf "•Fill_Q_Positions(NaN,\"imagePath\",\"\",\"\",$\"\")\r"
		Fill_Q_Positions(NaN,"imagePath","","",$"")
	elseif (stringmatch(ctrlName,"buttonPlotQpositionsRGBNoDepth"))
		printf "•MakeGraph_Q_Positions(Q_Positions_RGB)\r"
		MakeGraph_Q_Positions(Q_Positions_RGB)
	endif
	EnableDisableEWscanControls(GetUserData("microPanel","","EWscanPanelName"))
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
Static Function EnableDisableEWscanControls(win)				// here to enable/disable
	String win												// window (or sub window) to use
	Variable d

	Button buttonNewEWscan,win=$win,disable=0			// always OK to make a new data folder
	Button buttonFillQvsDepth,win=$win,disable=0			// always OK to load a Q-Depth surface
	Button buttonFillEvsDepth,win=$win,disable=0			// always OK to load a Energy-Depth surface
	Button buttonPseudoWhite,win=$win,disable=0			// always OK to make a Pseudo White Image

	d = WaveExists(Q_Depth)||WaveExists(Qhist) ? 0 : 2	// only if Q_Depth exists
	Button buttonFitQAt1Depth,win=$win,disable= d

	d = WaveExists(imageEsum) ? 0 : 2						// only if imageEsum exists
	Button buttonRePlotEsumDepth,win=$win,disable= d
	Button buttonMovieEsumsDepth,win=$win,disable= d

	d = strlen(WaveListClass("speImage*;ImageSummed*","*","DIMS:2")) ? 0 : 2
	Button buttonMakeMask,win=$win,disable=d			// OK to make mask if an imge is around
	d = MaskToClearExists() ? 0 : 2
	Button buttonClearMask,win=$win,disable= d

	d = WaveExists(Q_Depth) ? 0 : 2						// only if Q_Depth exists
	Button buttonRePlotQDepth,win=$win,disable= d

	d = WaveExists(keV_Depth) ? 0 : 2						// only if keV_Depth exists
	Button buttonRePlotEDepthSurf,win=$win,disable= d

	d = WaveExists(Q_Positions_RGB) ? 0 : 2				// only if Q_Positions_RGB exists
	Button buttonPlotQpositionsRGBNoDepth,win=$win,disable= d

	d = (strlen(WinList("*_Q_Depth",";","WIN:1")) && strlen(WinList("*_Qhist",";","WIN:1"))) ? 0 : 2
	Button buttonlayoutQDepth,win=$win,disable= d

	if (!str2num(GetUserData(win,"buttonNewEWscan","active")))	// new EWscan NOT active, set to current folder values
		NVAR d0=root:Packages:micro:EWscan:d0
		SVAR title=root:Packages:micro:EWscan:title, fldr=root:Packages:micro:EWscan:fldr
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