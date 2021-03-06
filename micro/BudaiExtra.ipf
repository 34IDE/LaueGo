#pragma rtGlobals=1		// Use modern global access method.
#pragma version = 0.02
#include "WinView", version>=1.81
#include  "EnergyWireScans", version>=0.8
#include <WaveSelectorWidget>


Menu "Macros"
	"Sum First Image of Each wire scan from an X-W scan", sumReconstructedImages("","",NaN,NaN,NaN)
	"Make X-Depth surface of intensity in a mask",MakeXvsDepthImageWithMask("","",NaN,NaN,NaN,NaN,$"")
	"-"
	"Sum all raw spe images in a range",sumAllImagesInRange("","","")
	"Make X-H surface of intensity in a mask",MakeXvsHimageWithMask("","","",NaN,$"")
	"-"
	"RGB composite image from separate images",MakeRGBcompositeImage()
End




// make a 2d wave where x-y are the X1 and H1, and the value is the sum of the image using mask
Function MakeXvsHimageWithMask(pathName,fileRoot,range,NX,mask)
	String pathName				// name of path that points to raw (un-reconstructed) images
	String fileRoot
	String range					// range of images to sum
	Variable NX						// number of X positions in 2d scan
	Wave mask
	String path
	if (strlen(pathName)<1 || strlen(fileRoot)<1 || numtype(str2num(range)) || !(NX>0) || !WaveExists(mask))	// bad inputs, ask
		pathName = SelectString(strlen(pathName),"rawPath",pathName)
		Prompt pathName, "name of path to raw images, probably 'rawPath'"
		Prompt fileRoot,"the file root, something like 'EYMO_WB2_'"
		Prompt range, "range of indicies to sum, in form \"1-7300\""
		Prompt NX, "number of X positions in X-H scan"
		String maskName								// name of mask wave
		String list = WaveListClass("imageMask","*","DIMS:2,BYTE:1")	// setup to use optional mask
		if (strlen(list)<1)
			Abort "no masks found in this data folder"
		endif
		Prompt maskName, "mask to use with image",popup,list
		DoPrompt "inputs", pathName,fileRoot,range,NX,maskName
		if (V_flag)
			return 1
		endif
		printf "MakeXvsHimageWithMask(\"%s\",\"%s\",\"%s\",%g,%s)\r",pathName,fileRoot,range,NX,maskName
		PathInfo $pathName
		if (!V_flag)
			NewPath/M="select folder containing raw images"/Q/Z $pathName
			if (V_flag)
				Abort "no path selected, stopping"
			endif
			PathInfo $pathName
			printf "set  '%s'  to  '%s'\r",pathName, S_path
		endif
		Wave mask = $maskName
	endif
	if (strlen(pathName)<1 || strlen(fileRoot)<1 || numtype(str2num(range)) || !(NX>0) || !WaveExists(mask))
		Abort "invalid inputs, doing nothing"
	endif
	PathInfo $pathName
	if (!V_flag)
		Abort "the path '"+pathName+"' does not exist, stopping"
	endif
	path = S_path
	String fullName, wName

	Variable N = ItemsInRange(range)
	Variable NH = round(N/NX)
	Make/N=(NX,NH)/O sumMaskedImage			// iH=[0,NH-1],  iX=[0,NX-1]
	String wnote = "waveClass=speImageSum"
	wnote = ReplaceStringByKey("mask",wnote,GetWavesDataFolder(mask,2),"=")
	Note/K sumMaskedImage,wnote
	sumMaskedImage = NaN

	Variable Xlo,Hlo,Xhi,Hhi
	Variable i,iX,iH, i0=str2num(range), mA
	Variable timer = startMSTimer
	for (i=str2num(range); !numtype(i); i=NextInRange(range,i))
		iX = mod(i-i0,NX)
		iH = round((i-i0-iX)/NX)
		sprintf fullName, "%s%s%d.SPE",path,fileRoot,i
		wName = WinViewReadROI(fullName,0,-1,0,-1)
		if (strlen(wName)<1)
			Abort "could not read image   "+fullName
		endif
		Wave image = $wName
		wnote = note(image)
		if (i==str2num(range))							// first point
			Xlo = NumberByKey("X1",wnote,"=")
			Hlo = (NumberByKey("Y1",wnote,"=")+NumberByKey("Z1",wnote,"="))/sqrt(2)
		elseif (i==lastInRange(range))						// last point
			Xhi = NumberByKey("X1",wnote,"=")
			Hhi = (NumberByKey("Y1",wnote,"=")+NumberByKey("Z1",wnote,"="))/sqrt(2)
		endif
		MatrixOp/O image = image * mask
		ImageStats/M=1 image
		mA = NumberByKey("mA",wnote,"=")				// normalize by incident beam current (relative to 100 mA)
		sumMaskedImage[iX][iH] = V_avg*V_npnts*(numtype(mA) ? 1 : 100/mA)
//		sumMaskedImage[iX][iH] = V_avg*V_npnts
		KillWaves/Z image
	endfor
	SetScale/I x Xlo,Xhi,"�m", sumMaskedImage
	SetScale/I y Hlo,Hhi,"�m", sumMaskedImage
	printf "execution took %s\r",Secs2Time(stopMSTimer(timer )/1e6,5,1)
	print "Made an X-H surface of intensity in a mask, and put into wave named  \"sumMaskedImage\""
	String win = StringFromList(0,FindGraphsWithWave(sumMaskedImage))
	if (strlen(win))
		DoWindow/F $win
	else
		Display /W=(577,218,1265,753)
		AppendImage sumMaskedImage
		ModifyImage sumMaskedImage ctab= {*,*,Rainbow,1}
		ModifyGraph mirror=2, minor=1
		Label bottom "X pm500  (\\U)"
		Label left "H pm500  (\\U)"
		SetAxis/A/R
	endif
End



// make a 2d wave where x-y are the X1 and depth, and the value is the sum of the image using mask
Function MakeXvsDepthImageWithMask(pathName,fileRoot,i0,i1,Nwire,Ndepth,mask)
	String pathName				// name of path that points to reconstructed images
	String fileRoot
	Variable i0,i1					// first and last file numbers for the original X-W scan scan, i0=1, i1=7381
	Variable Nwire					// number of points in one wire scan 121
	Variable Ndepth				// number of depths
	Wave mask
	String path
	if (strlen(pathName)<1 || strlen(fileRoot)<1 || !(i0>0) || !(i1>i0) || !(Nwire>0) || !(Ndepth>0) || !WaveExists(mask))	// bad inputs, ask
		pathName = SelectString(strlen(pathName),"reconPath",pathName)
		i0 = !(i0>0) ? 1 : i0
		i1 = !(i1>i0) ? 7381 : i1
		Prompt pathName, "name of path to reconstructed images, probably 'reconPath'"
		Prompt fileRoot,"the file root, something like 'EYMO_WB2_'"
		Prompt i0, "first index number in RAW X-Wire scan"
		Prompt i1, "last index number in RAW X-Wire scan"
		Prompt Nwire,"number of points in each wire scan"
		Prompt Ndepth,"number of depths in the reconstruction"

		String maskName								// name of mask wave
		String list = WaveListClass("imageMask","*","DIMS:2,BYTE:1")	// setup to use optional mask
		if (strlen(list)<1)
			Abort "no masks found in this data folder"
		endif
		Prompt maskName, "mask to use with image",popup,list
		DoPrompt "inputs", pathName,fileRoot,i0,i1,Nwire,Ndepth,maskName
		if (V_flag)
			return 1
		endif
		printf "MakeXvsDepthImageWithMask(\"%s\",\"%s\",%g,%g,%g,%g,%s)\r",pathName,fileRoot,i0,i1,Nwire,Ndepth,maskName
		PathInfo $pathName
		if (!V_flag)
			NewPath/M="select folder containing reconstructed images"/Q/Z $pathName
			if (V_flag)
				Abort "no path selected, stopping"
			endif
			PathInfo $pathName
			printf "set  '%s'  to  '%s'\r",pathName, S_path
		endif
		Wave mask = $maskName
	endif
	if (strlen(pathName)<1 || strlen(fileRoot)<1 || !(i0>0) || !(i1>i0) || !(Nwire>0) || !(Ndepth>0) || !WaveExists(mask))
		Abort "invalid inputs, doing nothing"
	endif
	PathInfo $pathName
	if (!V_flag)
		Abort "the path '"+pathName+"' does not exist, stopping"
	endif
	path = S_path
	String fullName, wName

	Variable N = i1-i0+1
	Variable NX = round(N/Nwire)
	Make/N=(Ndepth,NX)/O sumMaskedImage		// depth=[0,Ndepth-1],  x=[0,NX-1]
	String wnote = "waveClass=speImageSum"
	wnote = ReplaceStringByKey("mask",wnote,GetWavesDataFolder(mask,2),"=")
	Note/K sumMaskedImage,wnote

	Variable idepth,ix,i
	for (i=0,ix=i0;ix<i1;ix+=Nwire,i+=1)
		printf "starting ix=%d (of %d)s\r",ix,i1
		for (idepth=0;idepth<Ndepth;idepth+=1)
			sprintf fullName, "%s%s%d_%d.SPE",path,fileRoot,ix,idepth
			wName = WinViewReadROI(fullName,0,-1,0,-1)
			if (strlen(wName)<1)
				print fullName
				Abort "could not read image   "+fullName
			endif
			Wave image = $wName
			MatrixOp/O image = image * mask
			ImageStats/M=1 image
			sumMaskedImage[idepth][i] = V_avg*V_npnts
			KillWaves/Z image
		endfor
	endfor
	print "Made an X-Depth surface of intensity in a mask, and put into wave named  \"sumMaskedImage\""
End



// sum the raw images over range
Function sumAllImagesInRange(pathName,fileRoot,range)
	String pathName				// name of path that points to raw images
	String fileRoot
	String range					// range of images to sum
	String path
	if (strlen(pathName)<1 || strlen(fileRoot)<1 || strlen(range)<1)	// bad inputs, ask
		pathName = SelectString(strlen(pathName),"rawPath",pathName)
		Prompt pathName, "name of path to raw (unreconstructed) images, probably 'rawPath'"
		Prompt fileRoot,"the file root, something like 'EYMO_'"
		Prompt range, "range of indicies to sum, in form \"1-7300\""
		DoPrompt "inputs",pathName,fileRoot,range
		if (V_flag)
			return 1
		endif
		printf "sumAllImagesInRange(\"%s\",\"%s\",\"%s\")\r",pathName,fileRoot,range
		PathInfo $pathName
		if (!V_flag)
			NewPath/M="select folder containing raw images"/Q/Z $pathName
			if (V_flag)
				Abort "no path selected, stopping"
			endif
			PathInfo $pathName
			printf "set  '%s'  to  '%s'\r",pathName, S_path
		endif
	endif
	if (strlen(pathName)<1 || strlen(fileRoot)<1 || strlen(range)<1)
		Abort "invalid inputs, doing nothing"
	endif
	PathInfo $pathName
	if (!V_flag)
		Abort "the path '"+pathName+"' does not exist, stopping"
	endif
	path = S_path
	String fullName, wName

	KillWaves/Z sumImage
	if (exists("sumImage")!=1)
		sprintf fullName, "%s%s%d.SPE",path,fileRoot,str2num(range)
		wName = WinViewReadROI(fullName,0,-1,0,-1)	// load file into wName
		if (strlen(wName)<1)
			Abort "could not read first image    "+fullName
		endif
		Rename $wName sumImage
	endif
	Wave sumImage=sumImage
	Redimension/D sumImage
	String wnote = note(sumImage), class
	class = StringByKey("waveClass", note(sumImage),"=")
	wnote = ReplaceStringByKey("waveClass", wnote, class+"Sum","=")
	Note/K sumImage, wnote

	sumImage = 0
	Variable i, mA
	for (i=str2num(range); !numtype(i); i=NextInRange(range,i))
		sprintf fullName, "%s%s%d.SPE",path,fileRoot,i
		wName = WinViewReadROI(fullName,0,-1,0,-1)	// load file into wName
		if (strlen(wName)<1)
			Abort "could not read image   "+fullName
		endif
		Wave image = $wName
		mA = NumberByKey("mA",note(image),"=")		// normalize by incident beam current
		mA = numtype(mA) ? 1 : mA/100					// normalize to 100 mA
		image /= mA
		sumImage += image
		KillWaves/Z image
	endfor
	printf "Summed all raw spe images within range [%s] and put results into an image named \"sumImage\"\r",range
End




// sum first Image of each wire scan from an X-W scan, this is equivalent to summing all of the reconstructed images from an X-W scan
Function sumReconstructedImages(pathName,fileRoot,i0,i1,Nwire)
	String pathName				// name of path to raw data
	String fileRoot
	Variable i0,i1					// first and last file numbers for the original X-W scan scan, i0=1, i1=7381
	Variable Nwire					// number of points in one wire scan 121
	String path
	if (strlen(pathName)<1 || strlen(fileRoot)<1 || !(i0>0) || !(i1>i0) || !(Nwire>0))	// bad inputs, ask
		pathName = SelectString(strlen(pathName),"rawPath",pathName)
		i0 = !(i0>0) ? 1 : i0
		i1 = !(i1>i0) ? 7381 : i1
		Prompt pathName, "name of path to raw (NOT reconstructed) images, probably 'rawPath'"
		Prompt fileRoot,"the file root, something like 'EYMO_WB2_'"
		Prompt i0, "first index number in RAW X-Wire scan"
		Prompt i1, "last index number in RAW X-Wire scan"
		Prompt Nwire,"number of points in each wire scan"
		DoPrompt "inputs",pathName,fileRoot,i0,i1,Nwire
		if (V_flag)
			return 1
		endif
		printf "sumReconstructedImages(\"%s\",\"%s\",%g,%g,%g)\r",pathName,fileRoot,i0,i1,Nwire
		PathInfo $pathName
		if (!V_flag)
			NewPath/M="select folder containing RAW images"/Q/Z $pathName
			if (V_flag)
				Abort "no path selected, stopping"
			endif
			PathInfo $pathName
			printf "set  '%s'  to  '%s'\r",pathName, S_path
		endif
	endif
	if (strlen(pathName)<1 || strlen(fileRoot)<1 || !(i0>0) || !(i1>i0) || !(Nwire>0))
		Abort "invalid inputs, doing nothing"
	endif
	PathInfo $pathName
	if (!V_flag)
		Abort "the path '"+pathName+"' does not exist, stopping"
	endif
	path = S_path						// path to raw data
	String fullName, wName

	KillWaves/Z sumImage
	if (exists("sumImage")!=1)
		sprintf fullName, "%s%s%d.SPE",path,fileRoot,i0
		wName = WinViewReadROI(fullName,0,-1,0,-1)	// load file into wName
		if (strlen(wName)<1)
			Abort "could not read first image    "+fullName
		endif
		Rename $wName sumImage
	endif
	Wave sumImage=sumImage
	Redimension/D sumImage
	String wnote = note(sumImage), class
	class = StringByKey("waveClass", note(sumImage),"=")
	wnote = ReplaceStringByKey("waveClass", wnote, class+"sum","=")
	Note/K sumImage, wnote

	sumImage = 0
	Variable i
	for (i=i0;i<=i1;i+= Nwire)
		sprintf fullName, "%s%s%d.SPE",path,fileRoot,i
		wName = WinViewReadROI(fullName,0,-1,0,-1)	// load file into wName
		if (strlen(wName)<1)
			Abort "could not read image   "+fullName
		endif
		Wave image = $wName
		sumImage += image
		KillWaves/Z image
	endfor
	print "Summed Reconstructedimages from an X-Wire scan, and put results into an image named \"sumImage\""
End




//  ============================================================================  //
//  ============================ Start of RGB composite image ===========================  //

Function MakeRGBcompositeImage()
	NewDataFolder/O root:Packages
	NewDataFolder/O root:Packages:RGBcomposite
	String/G root:Packages:RGBcomposite:rgbName
	DoWindow/K RGBcompositeImagePanel
	NewPanel/K=1/W=(150,50,446,620)/N=RGBcompositeImagePanel
	Listbox listR win=RGBcompositeImagePanel, size={287,150}
	Listbox listG win=RGBcompositeImagePanel, size={287,150}
	Listbox listB win=RGBcompositeImagePanel, size={287,150}
	MakeListIntoWaveSelector("RGBcompositeImagePanel", "listR",selectionMode=WMWS_SelectionSingle,content=WMWS_Waves,listoptions="DIMS:2",nameFilterProc="FilterProc_speImageSum")
	MakeListIntoWaveSelector("RGBcompositeImagePanel", "listG",selectionMode=WMWS_SelectionSingle,content=WMWS_Waves,listoptions="DIMS:2",nameFilterProc="FilterProc_speImageSum")
	MakeListIntoWaveSelector("RGBcompositeImagePanel", "listB",selectionMode=WMWS_SelectionSingle,content=WMWS_Waves,listoptions="DIMS:2",nameFilterProc="FilterProc_speImageSum")
	ListBox listR,pos={1,20},size={287,150}
	ListBox listG,pos={1,190},size={287,150}
	ListBox listB,pos={1,360},size={287,150}
	SetDrawLayer UserBack
	SetDrawEnv fstyle= 1
	DrawText 2,19,"Red"
	SetDrawEnv fstyle= 1
	DrawText 2,190,"Green"
	SetDrawEnv fstyle=1
	DrawText 2,360,"Blue"
	CheckBox checkRedInvert,pos={50,5},size={65,14},title="invert Red",value= 0
	CheckBox checkGreenInvert,pos={50,176},size={75,14},title="invert Green",value= 0
	CheckBox checkBlueInvert,pos={50,346},size={68,14},title="invert Blue",value= 0
	SetVariable setRGBname,pos={13,515},size={218,18},title="RGB wave name",fSize=12
	SetVariable setRGBname,value=root:Packages:RGBcomposite:rgbName,bodyWidth=120
	Button RGBbuttonDoIt,pos={74,539},size={50,20},proc=MakeRGBButtonProc,title="DoIt"
	Button RGBbuttonCancel,pos={148,539},size={50,20},proc=MakeRGBButtonProc,title="Cancel"
	PauseForUser RGBcompositeImagePanel
End
//
Function FilterProc_speImageSum(aName, contents)
	String aName		// object name with full data folder path
	Variable contents	// content code as described for the content parameter
	if (contents!=WMWS_Waves)
		return 0		// reject if not a wave
	endif
	return stringmatch(StringByKey("waveClass",note($aName),"="),"speImageSum*")	//	return 0 to reject or 1 to accept
End
//
Function MakeRGBButtonProc(ctrlName) : ButtonControl
	String ctrlName
	if (stringmatch(ctrlName,"RGBbuttonCancel"))
		DoWindow/K RGBcompositeImagePanel
	elseif (stringmatch(ctrlName,"RGBbuttonDoIt"))
  		Wave Rimage = $StringFromList(0,WS_SelectedObjectsList("RGBcompositeImagePanel","listR"))
 		Wave Gimage = $StringFromList(0,WS_SelectedObjectsList("RGBcompositeImagePanel","listG"))
		Wave Bimage = $StringFromList(0,WS_SelectedObjectsList("RGBcompositeImagePanel","listB"))
		ControlInfo/W=RGBcompositeImagePanel checkRedInvert
		Variable redInvert = V_Value
		ControlInfo/W=RGBcompositeImagePanel checkGreenInvert
		Variable greenInvert = V_Value
		ControlInfo/W=RGBcompositeImagePanel checkBlueInvert
		Variable blueInvert = V_Value
		SVAR rgbName = root:Packages:RGBcomposite:rgbName
		String outName=CalcRGBcompositeImage(Rimage,Gimage,Bimage,redInvert,greenInvert,blueInvert,rgbName)
		printf "made RGB wave named  '%s'\r",outName
		DoWindow/K RGBcompositeImagePanel
		if (ItemsInList(GetRTStackInfo(0))<3 && strlen(outName) && !strlen(FindGraphsWithWave($outName)))
			Display /W=(58,77,570,565)
			AppendImage $outName
			ModifyGraph mirror=2
			SetAxis/A/R left
		endif
	endif
End
//
Function/S CalcRGBcompositeImage(Rimage,Gimage,Bimage,Ri,Gi,Bi,rgbName)
	Wave Rimage,Gimage,Bimage
	Variable Ri,Gi,Bi					// invert effect of red, green, or blue layers
	String rgbName
	if (!strlen(rgbName))
		DoAlert 0, "no name given for RGB wave"
		return ""
	endif
	rgbName = CleanupName(rgbName,0)
	if (WaveExists(rgbName))
		DoAlert 1, "the wave '"+rgbName+"' already exists, replace it?"
		if (V_flag!=1)
			return ""
		endif
	endif
	String imageName=""
	if (WaveExists(Rimage))
		imageName = GetWavesDataFolder(Rimage,2)
	elseif(WaveExists(Gimage))
		imageName = GetWavesDataFolder(Gimage,2)
	elseif(WaveExists(Bimage))
		imageName = GetWavesDataFolder(Bimage,2)
	else
		DoAlert 0,"no waves specified for Red, Green, or Blue planes"
		return ""
	endif
	Wave image = $imageName
	String wnote = note(image)
	Variable ix=DimSize(image,0), iy=DimSize(image,1)
	Make/N=(ix,iy,3)/O/W/U $rgbName
	Wave rgb = $rgbName
	CopyScales/I image,rgb
	rgb = 0

	if (WaveExists(Rimage))
		ImageStats/M=1  Rimage
		rgb[][][0] = round(65535*Rimage[p][q]/V_max)
	endif
	if (WaveExists(Gimage))
		ImageStats/M=1  Gimage
		rgb[][][1] = round(65535*Gimage[p][q]/V_max)
	endif
	if (WaveExists(Bimage))
		ImageStats/M=1  Bimage
		rgb[][][2] = round(65535*Bimage[p][q]/V_max)
	endif
	if (Ri)
		rgb[][][0] = 65535 - rgb[p][q][0]
	endif
	if (Gi)
		rgb[][][1] = 65535 - rgb[p][q][1]
	endif
	if (Bi)
		rgb[][][2] = 65535 - rgb[p][q][2]
	endif

	String class = StringByKey("waveClass",wnote,"=")+"RGB"
	wnote = ReplaceStringByKey("waveClass", wnote,class,"=")
	Note/K rgb, wnote
	return GetWavesDataFolder(rgb,2)
End

//  ============================= End of RGB composite image ===========================  //
//  ============================================================================  //
