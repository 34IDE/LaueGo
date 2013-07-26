#pragma rtGlobals=1		// Use modern global access method.
#pragma ModuleName=Masking
#pragma version = 1.01


Menu "GraphMarquee",dynamic
	"-"
	MarqueeMaskMenuItem("Invert this Mask"),/Q, InvertMask()
	help={"Invert the entire mask"}
	MarqueeMaskMenuItem("Add Marquee to Mask"),/Q, AddMarqueeToMask()
	help={"add the area of this marquee to the mask, (this will include the area of the marquee in the analysis)"}
	MarqueeMaskMenuItem("Subtract Marquee From Mask"),/Q,SubtractMarqueeFromMask()
	help={"remove the area of this marquee from the mask, (all pixels within this area will not be used in the analysis)"}
		MarqueeMaskMenuItem("Add Cursor A to mask",csr="A0"),/Q,AddCursorToMask("A")
		help={"add point at cursor A to the mask,  used to remove  a bad pixel from analysis"}
		MarqueeMaskMenuItem("Remove Cursor A from mask",csr="A1"),/Q,RemoveCursorFromMask("A")
		help={"remove point at cursor A from the mask,  undo the add pixel at cursor command"}
		MarqueeMaskMenuItem("Add Cursor B to mask",csr="B0"),/Q,AddCursorToMask("B")
		help={"add point at cursor B to the mask,  used to remove  a bad pixel from analysis"}
		MarqueeMaskMenuItem("Remove Cursor B from mask",csr="B1"),/Q,RemoveCursorFromMask("B")
		help={"remove point at cursor B from the mask,  undo the add pixel at cursor command"}
	MarqueeMaskMenuItem("Clear the Mask"),ClearMask()
	help={"reset the whole mask to empty, no points will be analyzed.  Use this to clear a mask so you can start addding points again"}
End


OverRide Function/S MarqueeMaskMenuItem(item,[csr])
	String item
	String csr

	Variable disable = 0											// 0=show,  1=show disabled,  2=do not show

	String class = ""
	String imageName = StringFromList(0,ImageNameList("",";"))
	if (strlen(imageName)<1)
		disable = 2
	endif
	Wave image = ImageNameToWaveRef("",imageName)
	if (!WaveExists(image))
		disable = 2
	else
		class = StringByKey("waveClass", note(image),"=")
		disable = stringmatch(class,"image&Mask") ? 0 : 1		// not an editable mask
	endif

	csr = SelectString(ParamIsDefault(csr),csr,"")
	Variable ic = char2num(UpperStr(csr[0]))
	if (strlen(csr) && disable)									// not good mask, and a cursor, always no show
		disable = 2
	elseif ((ic>=char2num("A") && ic<=char2num("J")))		// have a valid cursor letter
		Variable flag=!(!str2num(csr[1]))
		csr = num2char(ic)
		String list = CsrInfo($num2char(ic))
		Wave mask = $StringByKey("mask",GetUserData("","","images"),"=")
		Variable i=NumberByKey("POINT",list), j=NumberByKey("YPOINT",list)
		if (strlen(list)<1)
			disable = 2											// do not show, cursor not on image
		elseif(!WaveExists(mask))
			disable = 2											// do not show, cannot find mask
		elseif (!(i>=0 && j>=0))
			disable = 2											// do not show, cannot find point on image	
		else
			disable = flag==mask[i][j] ? 0 : 2
		endif
	endif

	if (disable==1)
		item = "("+item
	elseif (disable==2)
		item = ""
	endif
	return item
End
//OverRide Function/S MarqueeMaskMenuItem(item,[csr])
//	String item
//	String csr
//
//	Variable disable = 0
//	csr = SelectString(ParamIsDefault(csr),csr,"")
//	Variable ic = char2num(UpperStr(csr[0]))
//	if ((ic>=char2num("A") && ic<=char2num("J")))	// have a valid cursor letter
//		disable = strlen(CsrInfo($num2char(ic)))<1
//	endif
//
//	String imageName = StringFromList(0,ImageNameList("",";"))
//	if (strlen(imageName)<1)
//		return "("+item
//	endif
//	Wave image = ImageNameToWaveRef("",imageName)
//	if (!WaveExists(image))
//		return "("+item
//	endif
//	String class = StringByKey("waveClass", note(image),"=")
//
//	disable = disable || !stringmatch(class,"image&Mask")	// not an editable mask
//	if (disable)
//		item = "("+item
//		// item = "("
//	endif
//	return item
//End


//=======================================================================================
//=======================================================================================
//===================  This section is used to make an image mask by hand

Function MakeMask(image)
	Wave image
	if (!WaveExists(image) || WaveDims(image)!=2)
		String imageName = ""
		String list = WaveListClass("rawImage*;ImageSummed*","*","DIMS:2")
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
Function InvertMask()
	Wave image = $StringByKey("image",GetUserData("","","images"),"=")
	Wave mask = $StringByKey("mask",GetUserData("","","images"),"=")
	Wave both = $StringByKey("both",GetUserData("","","images"),"=")
	if (!WaveExists(image) || !WaveExists(mask) || !WaveExists(both))
		DoAlert 0, "Nothing done, this is probably the wrong kind of image, see 'Make Mask'"
		return 1
	endif
	mask = !mask
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
//
Function AddCursorToMask(cursorName)
	String cursorName
	Wave image = $StringByKey("image",GetUserData("","","images"),"=")
	Wave mask = $StringByKey("mask",GetUserData("","","images"),"=")
	Wave both = $StringByKey("both",GetUserData("","","images"),"=")

	if (!WaveExists(image) || !WaveExists(mask) || !WaveExists(both))
		DoAlert 0, "Nothing done, this is probably the wrong kind of image, see 'Make Mask'"
		return 1
	endif

	Variable ic = char2num(UpperStr(cursorName[0]))
	if (! (ic>=char2num("A") && ic<=char2num("J")) )
		DoAlert 0, "AddCursorToMask(), Nothing done, Cursor must be in {A, B, ... J}"
		return 1
	endif
	cursorName = num2char(ic)
	Variable i=NumberByKey("POINT",CsrInfo($cursorName))
	Variable j=NumberByKey("YPOINT",CsrInfo($cursorName))
	if (!(i>=0 && j>=0))
		DoAlert 0, "AddCursorToMask(), Could not determine pixel location of Cursor "+cursorName
		return 1
	endif

	if (mask[i][j])
		String str
		sprintf str, "Pixel [%d, %d] at cursor %s is already masked.  No change",i,j,cursorName
		DoAlert 0, str
	else
		mask[i][j] = 1
	endif
	Variable delta = masking#deltaForMaskOffset(image)
	both = mask[p][q] ? image[p][q]+delta : image[p][q]
	printf "Added [%d, %d] to mask\r",i,j
End
//
Function RemoveCursorFromMask(cursorName)
	String cursorName
	Wave image = $StringByKey("image",GetUserData("","","images"),"=")
	Wave mask = $StringByKey("mask",GetUserData("","","images"),"=")
	Wave both = $StringByKey("both",GetUserData("","","images"),"=")

	if (!WaveExists(image) || !WaveExists(mask) || !WaveExists(both))
		DoAlert 0, "Nothing done, this is probably the wrong kind of image, see 'Make Mask'"
		return 1
	endif

	Variable ic = char2num(UpperStr(cursorName[0]))
	if (! (ic>=char2num("A") && ic<=char2num("J")) )
		DoAlert 0, "RemoveCursorFromMask(), Nothing done, Cursor must be in {A, B, ... J}"
		return 1
	endif
	cursorName = num2char(ic)
	Variable i=NumberByKey("POINT",CsrInfo($cursorName))
	Variable j=NumberByKey("YPOINT",CsrInfo($cursorName))
	if (!(i>=0 && j>=0))
		DoAlert 0, "RemoveCursorFromMask(), Could not determine pixel location of Cursor "+cursorName
		return 1
	endif

	if (!mask[i][j])
		String str
		sprintf str, "Pixel [%d, %d] at cursor %s was NOT masked.  No change",i,j,cursorName
		DoAlert 0, str
	else
		mask[i][j] = 0
	endif
	Variable delta = masking#deltaForMaskOffset(image)
	both = mask[p][q] ? image[p][q]+delta : image[p][q]
	printf "Remove [%d, %d] from mask\r",i,j
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
