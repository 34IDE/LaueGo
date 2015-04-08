#pragma rtGlobals=3		// Use modern global access method.
#pragma version = 2.08
#pragma IgorVersion = 6.3
#pragma ModuleName=GMarkers
#include "GizmoUtility", version>=0.16


Static strConstant Q_UNITS_LIST = "1/nm;nm\\S-1\\M"
Static Constant NUMBER_OF_GIZMO_MARKERS = 8


Static Function AfterFileOpenHook(refNum,file,pathName,type,creator,kind)
	Variable refNum, kind
	String file,pathName,type,creator
	if ((kind==1) || (kind==2))		// an experiment (packed or unpacked)
		GMarkers#InitGizmoMarkers()
	endif
	return 0
End
Static Function IgorStartOrNewHook(IgorApplicationNameStr)
	String IgorApplicationNameStr
	GMarkers#InitGizmoMarkers()
	return 0
End



//  ======================================================================================  //
//  ================================== Start of Markers ==================================  //

Static Structure GizmoMarkerInfoStruct
	String win
	Variable MarkerNum			// current marker number slected on panel
	STRUCT GizmoMarkerInfoStruct1 M[NUMBER_OF_GIZMO_MARKERS]
EndStructure
//
Static Structure GizmoMarkerInfoStruct1
	int16 used
	String MarkerName					// something like 1 red
	String wName						// name of wave that the marker is on
	String xUnit,yUnit,zUnit		// units for xyz axes
	double x0,y0,z0					// position in xyz
	double point						// point number that marker is on
	double intensity					// value of wName at point
	double ix,iy,iz					// indecies into wName to point number
	double h,k,l						// hkl corresponding to xyz (assuming recip is valid)
	double recip[9]					// reciprocal lattice (if it exists)
EndStructure




// provides external access to the marker information, provides key:value, for all the information
Function/S GizmoMarkerInfo2keyVals()
	STRUCT GizmoMarkerInfoStruct info
	info.win = "GizmoScatterMarkerPanel"
	if (GizmoMarkerInfo(info)<1)			// no markers dissplayed
		return ""
	endif

	String out=""
	out = ReplaceStringByKey("win",out,info.win)
	out = ReplaceNumberByKey("ActiveMarkerNum",out,info.MarkerNum)

	String str, sn
	Variable i, h,k,l, intensity, ix,iy,iz, x0,y0,z0
	for (i=0;i<NUMBER_OF_GIZMO_MARKERS;i+=1)
		if (info.M[i].used)
			sn = num2istr(i)
			x0=info.M[i].x0 ; 	y0=info.M[i].y0 ;	 z0=info.M[i].z0
			ix=info.M[i].ix;		iy=info.M[i].iy;		iz=info.M[i].iz

			out = ReplaceStringByKey("name"+sn,out,info.M[i].MarkerName)
			out = ReplaceStringByKey("wave"+sn,out,info.M[i].wName)

			if (numtype(x0+y0+z0)==0)
				sprintf str,"%g,%g,%g",x0,y0,z0
				out = ReplaceStringByKey("xyz"+sn,out,str)
	
				if (numtype(ix+iy+iz)==0)
					sprintf str,"%g,%g,%g",ix,iy,iz
					out = ReplaceStringByKey("ixyz"+sn,out,str)
				elseif (numtype(info.M[i].point)==0)
					out = ReplaceNumberByKey("point"+sn,out,info.M[i].point)
				endif
	
				intensity=info.M[i].intensity
				if (numtype(intensity)==0)
					out = ReplaceNumberByKey("intensity"+sn,out,intensity)
				endif
				h=info.M[i].h;	k=info.M[i].k;	l=info.M[i].l
				if (numtype(h+k+l)==0)
					sprintf str,"%g,%g,%g",h,k,l
					out = ReplaceStringByKey("hkl"+sn,out,str)
				endif
			endif
		endif
	endfor
	return out
End




Function MakeGizmoScatterMarkerPanel() : Panel	// Create the Cut Plane Panel.  Doesn't create duplicate panels.
	if (ItemsInList(WinList("GizmoScatterMarkerPanel",";","WIN:64")))
		DoWindow/F GizmoScatterMarkerPanel				// If Marker Panel exists, just bring it to front and quit
		return 0
	endif

	Variable left=1094,top=245,width=150
	NewPanel/K=1/W=(left,top,left+width,top+340)/N=GizmoScatterMarkerPanel  as "Marker"
	SetDrawLayer UserBack
	SetDrawEnv fstyle= 1
	DrawText 7,19,"Scatter Marker (Gizmo)"
	SetDrawEnv fsize= 10
	DrawText 10,40,"absolute position"
	SetDrawEnv fsize= 10
	DrawText 10,200,"along wave"
	SetDrawEnv linethick= 3
	DrawLine 0,132,width,132
	SetDrawEnv linethick= 3
	DrawLine 0,181,width,181
	SetDrawEnv linethick= 3
	DrawLine 0,280,width,280
	SetDrawEnv linethick= 3
	DrawLine 0,280,width,280

	SetVariable Xabsolute,pos={20,40},size={101,19},bodyWidth=90,title="X",fSize=12
	SetVariable Xabsolute,proc=GMarkers#GizmoScatterMarkerSetVarProc,value= _NUM:NaN
	SetVariable Yabsolute,pos={20,63},size={101,19},bodyWidth=90,title="Y",fSize=12
	SetVariable Yabsolute,proc=GMarkers#GizmoScatterMarkerSetVarProc,value= _NUM:NaN
	SetVariable Zabsolute,pos={20,86},size={101,19},bodyWidth=90,title="Z",fSize=12
	SetVariable Zabsolute,proc=GMarkers#GizmoScatterMarkerSetVarProc,value= _NUM:NaN
	SetVariable xyzStep,pos={19,108},size={103,19},bodyWidth=70,proc=GMarkers#GizmoScatterMarkerSetVarProc,title="Æxyz"
	SetVariable xyzStep,fSize=12,limits={0,inf,1},value= _NUM:1

	CheckBox OnPointCheckBox,pos={6,141},size={61,32},title="stay on\rpoints"
//	CheckBox OnPointCheckBox,fSize=12,value=0, proc=GMarkers#OnPointCheckProc
	CheckBox OnPointCheckBox,fSize=12,value=0
	Button MarkerInfoButton,pos={76,143},size={25,30},proc=GMarkers#GizmoMarkerInfoButtonProc
	Button MarkerInfoButton,fSize=24,fStyle=1,font="Bookman Old Style",title="i" // font "Hoefler Text" doesn't exist on Windows
//	Button MarkerInfoButton,pos={80,143},size={20,30},proc=GMarkers#GizmoMarkerInfoButtonProc
//	Button MarkerInfoButton,fSize=24,title="?"
	Button FitPeakButton,pos={110,138},size={35,40},proc=GMarkers#GizmoScatterMarkerButtonProc,title="Fit\rpeak"

	SetVariable alongWave,pos={10,204},size={124,19},bodyWidth=90,proc=GMarkers#GizmoScatterMarkerSetVarProc,title="point"
	SetVariable alongWave,fSize=12,limits={0,inf,1},value= _NUM:0

	TitleBox colorTitleBox,pos={0,234},size={width,20},title="                                        "
	TitleBox colorTitleBox,labelBack=(0,0,0),fSize=6,frame=0
	TitleBox colorTitleBox,anchor= LC,fixedSize=1
	PopupMenu MarkerNumberPopup,pos={33,235},size={60,20},proc=GMarkers#GizmoChangeMarkerPopMenuProc
	PopupMenu MarkerNumberPopup,fSize=18
	PopupMenu MarkerNumberPopup,mode=1,popvalue="0 black",value= #"\"0 black;1 red;2 green;3 blue;4 magenta;5 cyan;6 yellow;7 white\""
	SetVariable MarkerSizeSetVar,pos={7,259},size={74,16},bodyWidth=50,proc=GMarkers#GizmoMarkerSizeSetVarProc,title="Size"
	SetVariable MarkerSizeSetVar,fSize=10,limits={0,10,0.1},value= _NUM:0.5
	Button ShowHideMarkerButton,pos={97,257},size={50,20},proc=GMarkers#GizmoMarkerShowHideButtonProc,title="Show"
	Button ShowHideMarkerButton,fSize=12

	PopupMenu waveSelectPopup,pos={3,287},size={100,20}
	PopupMenu waveSelectPopup,mode=1,popvalue="",value=""
	Button removeMarkerButton,pos={17,313},size={110,21},proc=GMarkers#GizmoScatterMarkerButtonProc,title="add marker"

	ValDisplay IntensDisp,pos={12,340},size={109,17},title=" I",fSize=12			// Only show hkl & Intensity when a reciprocal lattice is present
	ValDisplay IntensDisp,limits={0,0,0},barmisc={0,1000},value= _NUM:0
	ValDisplay Hvaldisp,pos={12,360},size={109,17},title="H",fSize=12
	ValDisplay Hvaldisp,limits={0,0,0},barmisc={0,1000},value= _NUM:0
	ValDisplay Kvaldisp,pos={12,380},size={109,17},title="K",fSize=12
	ValDisplay Kvaldisp,limits={0,0,0},barmisc={0,1000},value= _NUM:0
	ValDisplay Lvaldisp,pos={12,400},size={109,17},title="L",fSize=12
	ValDisplay Lvaldisp,limits={0,0,0},barmisc={0,1000},value= _NUM:0

	SetWindow kwTopWin,hook(update)=GMarkers#GizmoScatterMarkerUpdateHook
	return 0
End


Function ShowHideGizmoMarker(MarkerNum,show,[info])			// show/hide a marker
	Variable MarkerNum					// marker number to change [0,NUMBER_OF_GIZMO_MARKERS-1]
	Variable show							// true=show, false=hide
	STRUCT GizmoMarkerInfoStruct &info

	MarkerNum = round(MarkerNum)
	if (!(MarkerNum==limit(MarkerNum,0,NUMBER_OF_GIZMO_MARKERS-1)))
		return 1
	endif

	if (ParamIsDefault(info))
		STRUCT GizmoMarkerInfoStruct infoLocal
		infoLocal.win = "GizmoScatterMarkerPanel"
		GizmoMarkerInfo(infoLocal)
		Wave scatter = $(infoLocal.M[MarkerNum].wName)
	else
		Wave scatter = $(info.M[MarkerNum].wName)
	endif
	if (!WaveExists(scatter))
		return 1
	endif

	String path=GetWavesDataFolder(scatter,1)
	Wave sizeW=$(path+"gizmoScatterMarkerArraySize"), sizeW0=$(path+"gizmoScatterMarkerArraySize0")
	sizeW[MarkerNum][] = show ? sizeW0[MarkerNum] : NaN

	return 0
End


Function moveGizmoMarkerToPoint(point,scatter,MarkerNum)	// used by external routines to move the marker
	Variable point
	Wave scatter
	Variable MarkerNum

	if (!GizmoScatterMarkerDisplayed())
		return 1
	elseif (!WaveExists(scatter))
		return 1
	elseif (WaveDims(scatter)!=2 || DimSize(scatter,1)!=3)
		return 1
	elseif (!(point>=0 && point<DimSize(scatter,0)))
		return 1
	endif

	Wave marker=$GizmoGetScatterMarkerWave()
	if (!WaveExists(marker))
		return 1
	endif

	Variable mx=scatter[point][0], my=scatter[point][1], mz=scatter[point][2]
	if (numtype(mx+my+mz))
		// print "NaN"
	endif
	marker[MarkerNum][0] = mx
	marker[MarkerNum][1] = my
	marker[MarkerNum][2] = mz
	GizmoMarkerPanelUpdate()
	return 0
End


Function moveGizmoMarkerToXYZ(mx,my,mz,scatter,MarkerNum)	// used by external routines to move the marker
	Variable mx,my,mz
	Wave scatter							// This can be defaulted
	Variable MarkerNum					// This defaults to the current marker

	if (!GizmoScatterMarkerDisplayed())
		return 1
	endif
	if (numtype(mx+my+mz))
		return 1
	endif

	if (!(MarkerNum>=0 && MarkerNum<=7))
		ControlInfo/W=GizmoScatterMarkerPanel MarkerNumberPopup
		MarkerNum = V_Value-1
	endif
	if (!(MarkerNum>=0 && MarkerNum<=7))
		return 1
	endif

	if (!WaveExists(scatter))
		ControlInfo/W=GizmoScatterMarkerPanel waveSelectPopup
		Wave scatter=$S_Value
	endif
	if (!WaveExists(scatter))
		return 1
	elseif (!(WaveDims(scatter)==2 && DimSize(scatter,1)==3) && WaveDims(scatter)!=3)
		return 1
	endif

	Wave marker=$GizmoGetScatterMarkerWave()
	if (!WaveExists(marker))
		return 1
	endif
	marker[MarkerNum][0] = mx
	marker[MarkerNum][1] = my
	marker[MarkerNum][2] = mz
	GizmoMarkerPanelUpdate()
	return 0
End


Static Function GizmoScatterMarkerSetVarProc(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	if (sva.eventCode!=1 && sva.eventCode!=2)	// mouse up and Enter key
		return 0
	endif

	ControlInfo/W=$(sva.win) waveSelectPopup
	Wave scatter=$S_Value
	Variable dx=1, dy=1, dz=1, xLo=-inf,yLo=-inf,zLo=-inf, xHi=Inf,yHi=Inf,zHi=Inf, Nx,Ny,Nz
	if (WaveExists(scatter))
		if (WaveDims(scatter)==3)			// scatter is a 3D array, so use DimDelta for step sizes
			dx = DimDelta(scatter,0)
			dy = DimDelta(scatter,1)
			dz = DimDelta(scatter,2)
			Nx=DimSize(scatter,0)
			Ny=DimSize(scatter,1)
			Nz=DimSize(scatter,2)
			xLo = DimOffset(scatter,0)
			yLo = DimOffset(scatter,1)
			zLo = DimOffset(scatter,2)
			xHi = xLo + dx*(Nx-1)
			yHi = yLo + dy*(Ny-1)
			zHi = zLo + dz*(Nz-1)
		endif
	endif

	if (stringmatch(sva.ctrlName,"xyzStep"))			// xyzStep changed, use its value to set step sizes
		ControlInfo/W=$(sva.win) xyzStep
		Variable step = ((numtype(V_Value)==0 && V_Value>0) ? V_Value : NaN
		if (numtype(step)==0)									// a valid step was given, use it
			SetVariable Xabsolute, limits={xLo,xHi,step}, win=$(sva.win)
			SetVariable Yabsolute, limits={yLo,yHi,step}, win=$(sva.win)
			SetVariable Zabsolute, limits={zLo,zHi,step}, win=$(sva.win)
		elseif (numtype(step) && WaveDims(scatter)==3)// step is bad, but know 3D scaling
			SetVariable Xabsolute, limits={xLo,xHi,dx}, win=$(sva.win)
			SetVariable Yabsolute, limits={yLo,yHi,dy}, win=$(sva.win)
			SetVariable Zabsolute, limits={zLo,zHi,dz}, win=$(sva.win)
		else
			step = 1													// give up, use a deault step of 1
			SetVariable Xabsolute, limits={xLo,xHi,step}, win=$(sva.win)
			SetVariable Yabsolute, limits={yLo,yHi,step}, win=$(sva.win)
			SetVariable Zabsolute, limits={zLo,zHi,step}, win=$(sva.win)
		endif
	endif

	ControlInfo/W=$(sva.win) MarkerNumberPopup
	Variable MarkerNum = V_Value-1
	Wave marker=$GizmoGetScatterMarkerWave()
	Variable point=NaN, mx=NaN, my=NaN, mz=NaN
	if (!WaveExists(marker) || numtype(sva.dval))
		return 1
	elseif (stringmatch(sva.ctrlName,"Xabsolute"))
		marker[MarkerNum][0] = sva.dval
	elseif (stringmatch(sva.ctrlName,"Yabsolute"))
		marker[MarkerNum][1] = sva.dval
	elseif (stringmatch(sva.ctrlName,"Zabsolute"))
		marker[MarkerNum][2] = sva.dval
	elseif (stringmatch(sva.ctrlName,"alongWave"))
		if (!WaveExists(marker) || !WaveExists(scatter))
			return 1
		endif
		Variable i = sva.dval
		if (!(i>=0))
			return 1
		endif
		point = i
	endif

	ControlInfo/W=$(sva.win) OnPointCheckBox
	if (V_Value && numtype(point))
		point = GizmoScatterMarkerGetNearPoint(marker[MarkerNum][0],marker[MarkerNum][1],marker[MarkerNum][2])
	endif
	if (numtype(point)==0 || V_Value)								// given a point, or told to stay on points
		Variable N
		if (WaveDims(scatter)==2 && DimSize(scatter,1)==3)	// triplets
			N = DimSize(scatter,0)
			i = limit(point,0,N-1)
			mx = scatter[i][0]
			my = scatter[i][1]
			mz = scatter[i][2]
		else
			N = numpnts(scatter)
			i = limit(point,0,N-1)
			Variable ix = mod(i,Nx)
			Variable iy = mod(floor(i/Nx),Ny)
			Variable iz = floor(i/(Nx*Ny))
			mx = limit(ix*dx+xLo,xLo,xHi)
			my = limit(iy*dy+yLo,yLo,yHi)
			mz = limit(iz*dz+zLo,zLo,zHi)
		endif
		if (!(i>=0 && i<N))
			return 1
		endif

		marker[MarkerNum][0] = mx
		marker[MarkerNum][1] = my
		marker[MarkerNum][2] = mz
		if (i!=sva.dval)
			SetVariable alongWave,value= _NUM:i, win=$(sva.win)
		endif
	endif
	GizmoMarkerPanelUpdate()

	return 0
End
//
Static Function GizmoScatterMarkerButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	if (ba.eventCode != 2)	// mouse up event
		return 0
	elseif (itemsInList(WinList("*",";","WIN:4096"))<1)
		return 0
	endif

	ControlInfo/W=$(ba.win) MarkerNumberPopup
	Variable MarkerNum = V_Value-1

	Variable displayed=GizmoMarkerPanelUpdate()
	Wave marker=$GizmoGetScatterMarkerWave()

	Variable anAdd=0
	if (stringmatch(ba.ctrlName,"removeMarkerButton"))
		ControlInfo/W=$(ba.win) removeMarkerButton
		String title = StringByKey("title", S_recreation,"=",",")
		title = ReplaceString("\"",title,"")
		if (stringmatch(title,"add marker*"))			// determine desired action
			anAdd = 1
		endif
	endif
	if (!WaveExists(marker) && !anAdd)
		return 0
	endif

	ControlInfo/W=$(ba.win) waveSelectPopup
	Wave scatter=$S_Value
	if (stringmatch(ba.ctrlName,"FitPeakButton"))	// fit a 3D Gaussian peak at the marker position
		if (WaveDims(scatter)==3)							// this only works for 3D waves (NO triplit xyz)
			Make/N=3/D/FREE fitXYZ
			fitXYZ = marker[MarkerNum][p]
			FitPeakAt3Dmarker(scatter,fitXYZ,NaN,printIt=1)
		endif

	elseif (stringmatch(ba.ctrlName,"removeMarkerButton"))
		ControlInfo/W=$(ba.win) removeMarkerButton
		title = StringByKey("title", S_recreation,"=",",")
		title = ReplaceString("\"",title,"")
		if (stringmatch(title,"add marker*"))			// determine desired action
			GizmoAddScatterMarker(scatter=scatter)
		elseif (stringmatch(title,"remove marker*"))
			GizmoRemoveScatterMarkers()
		endif
	endif

	GizmoMarkerPanelUpdate()
	return 0
End
//
Static Function GizmoChangeMarkerPopMenuProc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa
	if (pa.eventCode !=2)
		return 0
	endif

	String win=pa.win
	ControlInfo/W=$win waveSelectPopup
	Wave scatter=$S_Value
	if (!WaveExists(scatter))
		return 0
	endif
	String path=GetWavesDataFolder(scatter,1)
	GizmoMarkerPanelUpdate()
	return 0
End
//
//Static Function OnPointCheckProc(cba) : CheckBoxControl
//	STRUCT WMCheckboxAction &cba
//	if (cba.eventCode !=2)
//		return 0
//	endif
//	return 0
//End
//
Static Function GizmoMarkerSizeSetVarProc(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable size = sva.dval
			break
		default:
			return 0
	endswitch

	String win=sva.win
	ControlInfo/W=$win waveSelectPopup
	Wave scatter=$S_Value
	if (!WaveExists(scatter))
		return 0
	endif
	String path=GetWavesDataFolder(scatter,1)
	Wave sizeW=$(path+"gizmoScatterMarkerArraySize"), sizeW0=$(path+"gizmoScatterMarkerArraySize0")
	if (WaveExists(sizeW) && WaveExists(sizeW0))
		ControlInfo/W=$win MarkerNumberPopup
		Variable MarkerNum = V_Value-1
		sizeW0[MarkerNum] = size
//		sizeW[MarkerNum][] = numtype(sizeW[MarkerNum][q]) ? NaN : size
		sizeW[MarkerNum][] = size
	endif
	return 0
End
//
Static Function GizmoMarkerShowHideButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	if (ba.eventCode != 2)
		return 0
	endif

	String win=ba.win
	ControlInfo/W=$win waveSelectPopup
	Wave scatter=$S_Value
	if (!WaveExists(scatter))
		return 0
	endif

	String path=GetWavesDataFolder(scatter,1)
	Wave sizeW=$(path+"gizmoScatterMarkerArraySize"), sizeW0=$(path+"gizmoScatterMarkerArraySize0")
	if (WaveExists(sizeW) && WaveExists(sizeW0))
		ControlInfo/W=$win $(ba.ctrlName)
		Variable i=strsearch(S_recreation, ",title=",0)
		Variable show = strsearch(S_recreation[i+8,Inf], "Show",0)==0
		ControlInfo/W=$win MarkerNumberPopup
		Variable MarkerNum = V_Value-1
		sizeW[MarkerNum][] = show ? sizeW0[MarkerNum] : NaN
	endif

	Wave gizmoScatter=$(path+"gizmoScatterMarkerArray")
	if (WaveExists(gizmoScatter))				// start new markers at center if its position is NaN
		if (numtype(gizmoScatter[MarkerNum][0]+gizmoScatter[MarkerNum][1]+gizmoScatter[MarkerNum][2]))
			Wave center=centerOf3Ddata(scatter)
			if (WaveExists(center))				// put up new marker at center of volume
				gizmoScatter[MarkerNum][] = center[q]
			endif
		endif
	endif

	GizmoMarkerPanelUpdate()
	return 0
End
//
Static Function GizmoMarkerInfoButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	if (ba.eventCode != 2)
		return 0
	endif

	STRUCT GizmoMarkerInfoStruct info
	info.win = ba.win
	if (GizmoMarkerInfo(info)<1)	
		print "No Markers Displayed"
		return 0
	endif
	Variable MarkerNum=info.MarkerNum
	if (ba.eventMod & (2|8) || !(info.M[MarkerNum].used))	// print all info when (Shift or Comamnd Key), or for an unsued marker
		print "¥PrintGizmoMarkerInfoAll()"
		PrintGizmoMarkerInfoStruct(info)
		return 0
	endif
	if (!(info.M[MarkerNum].used))
		return 0
	endif
	String wName=info.M[MarkerNum].wName
	String MarkerName=info.M[MarkerNum].MarkerName

	if (strlen(GetRTStackInfo(2))<1)
		print "¥GizmoMarkerInfoPrint()"
	endif
	printf "for marker%d,  '%s' on wave '%s'\r",MarkerNum,MarkerName,wName

	Variable x0=info.M[MarkerNum].x0, y0=info.M[MarkerNum].y0, z0=info.M[MarkerNum].z0
	Variable point=info.M[MarkerNum].point
	Variable ix=info.M[MarkerNum].ix, iy=info.M[MarkerNum].iy, iz=info.M[MarkerNum].iz
	String xUnit=info.M[MarkerNum].xUnit, yUnit=info.M[MarkerNum].yUnit, zUnit=info.M[MarkerNum].zUnit
	Variable intensity=info.M[MarkerNum].intensity
	if (numtype(x0+y0+z0)==0)
		if (stringmatch(xUnit,yUnit) && stringmatch(xUnit,zUnit))
			printf "   at {%g, %g, %g}%s",x0,y0,z0,xUnit
		else
			printf "   at {%g (%s), %g (%s), %g (%s)}",x0,xUnit,y0,yUnit,z0,zUnit
		endif
		if (numtype(ix+iy+iz)==0)
			printf "   nearest actual point is [%g, %g, %g],   point#=%d",ix,iy,iz,point
		else
			printf "   nearest actual point is point#=%d",point
		endif
		if (numtype(intensity)==0)
			printf ",   intensity=%g",intensity
		endif
		printf "\r"
	endif

	Variable h=info.M[MarkerNum].h, k=info.M[MarkerNum].k, l=info.M[MarkerNum].l
	if (numtype(h+k+l)==0)
		printf "   (hkl) = (%g, %g, %g)\r",h,k,l
	endif

	String list=""
	Variable i
	for (i=0;i<NUMBER_OF_GIZMO_MARKERS;i+=1)
		if (i!=info.MarkerNum && info.M[i].used)
			list += info.M[i].MarkerName+";"
		endif
	endfor
	if (strlen(list)>0)					// I can show differences from info.MarkerNum too!
		String name2=StringFromList(0,list)
		if (itemsInList(list)>1)
			Prompt name2,"Other Marker to compare to",popup,list
			DoPrompt "Other Marker Compare", name2
			if (V_flag)
				return 1
			endif
		endif
		Variable num2 = str2num(name2)
		Make/N=3/D/FREE distance
		distance[0] = info.M[MarkerNum].x0 - info.M[num2].x0
		distance[1] = info.M[MarkerNum].y0 - info.M[num2].y0
		distance[2] = info.M[MarkerNum].z0 - info.M[num2].z0
		printf "from '%s'  -->  '%s':\r",info.M[num2].MarkerName,info.M[MarkerNum].MarkerName
		if (stringmatch(xUnit,yUnit) && stringmatch(xUnit,zUnit))
			printf "   distance = {%g, %g, %g}%s,  |distance|=%g\t\t",distance[0],distance[1],distance[2],xUnit,norm(distance)
		else
			printf "   distance = {%g (%s), %g (%s), %g (%s)}%s,  |distance|=%g\t\t",distance[0],xUnit,distance[1],yUnit,distance[2],zUnit,norm(distance)
		endif
		if (numtype(ix+iy+iz + info.M[num2].ix + info.M[num2].iy + info.M[num2].iz)==0)
			printf "   Æ(ijk) = [%g, %g, %g],   Æpoint = %g\r",(ix - info.M[num2].ix), (ix - info.M[num2].iy), (iz - info.M[num2].iz), (point - info.M[num2].point)
		else
			printf "   Æpoint = %g\r",(point - info.M[num2].point)
		endif
		if (numtype(h+k+l + info.M[num2].h + info.M[num2].k + info.M[num2].l)==0)
			Variable dh = h - info.M[num2].h, dk = k - info.M[num2].k, dl = l - info.M[num2].l
			printf "   Æ(hkl) = (%.3g, %.3g, %.3g),  |Æ(hkl)| = %.3g (rlu)\r",dh, dk, dl,sqrt(dh^2+dk^2+dl^2)
		endif

		if (exists("InitLatticeSymPackage")==6)
			Make/N=9/D/FREE rl0, rl1						// consider printing the rotation between two reicprocal lattices
			rl1 = info.M[MarkerNum].recip[p]
			rl0 = info.M[num2].recip[p]
			Redimension/N=(3,3) rl0, rl1
			list =  Compare2ReciprocalLatticesList(rl0,rl1)
			Variable angle=NumberByKey("axisRotationAngle",list,"=")
			Wave scatter=$(info.M[MarkerNum].wName)
			if (angle>0)
				printf "   rotation is about the axis {XYZ} = %s  by  %.3g¡\r",StringByKey("axisDirectionXYZ",list,"="),angle
				printf "   in both crystals, the axis is an hkl={%s},    which is %.2f¡ away from the hkl={%s}\r",StringByKey("axisHKL",list,"="),NumberByKey("axisHKLdeviation",list,"="),StringByKey("axisHKL0",list,"=")
			elseif (WaveDims(scatter)==1 && angle==0)		// this only makes sense for 1D, 3D waves will always give zero
				printf "   the reciprocal lattices are identical\r"
			endif
		endif
	endif
	return 0
End
//
Function GizmoMarkerInfoPrint()
	STRUCT WMButtonAction ba
	ba.eventCode = 2
	ba.win = "GizmoScatterMarkerPanel"
	GizmoMarkerInfoButtonProc(ba)
End
//
Static Function GizmoScatterMarkerUpdateHook(s)		// update all values on Marker Panel. Needed since Marker Panel can serve many Gizmos
	STRUCT WMWinHookStruct &s

	Variable hookResult = 0
	switch(s.eventCode)
		case 0:				// Activate
		case 6:				// resize
//		case 12:			// moved
			GizmoMarkerPanelUpdate()
			break
	endswitch
	return hookResult		// 0 if nothing done, else 1
End
//
Static Function GizmoMarkerPanelUpdate()
	if (strlen(WinList("*","","WIN:4096"))<1)	// if no gizmos, disable everything
		PopupMenu waveSelectPopup,disable=2, win=GizmoScatterMarkerPanel
		SetVariable xyzStep,disable=2, win=GizmoScatterMarkerPanel
		Button FitPeakButton,disable=1, win=GizmoScatterMarkerPanel
		SetVariable alongWave,disable=2, win=GizmoScatterMarkerPanel
		Button removeMarkerButton,disable=2, win=GizmoScatterMarkerPanel
		CheckBox OnPointCheckBox,disable=2, win=GizmoScatterMarkerPanel
		SetVariable MarkerSizeSetVar,disable=2, win=GizmoScatterMarkerPanel
		Button ShowHideMarkerButton ,disable=2, win=GizmoScatterMarkerPanel
		SetVariable Xabsolute,disable=2, win=GizmoScatterMarkerPanel
		SetVariable Yabsolute,disable=2, win=GizmoScatterMarkerPanel
		SetVariable Zabsolute,disable=2, win=GizmoScatterMarkerPanel
		SetVariable alongWave,disable=2, win=GizmoScatterMarkerPanel
		Button MarkerInfoButton,disable=2, win=GizmoScatterMarkerPanel
		updatePanelHKL($"",NaN,NaN,NaN)
		SetVariable MarkerSizeSetVar,disable=2, win=GizmoScatterMarkerPanel
		PopupMenu MarkerNumberPopup,disable=2, win=GizmoScatterMarkerPanel
		Button FitPeakButton,disable=2				// only show this button when marking a 3D wave
		return 0
	endif

	Variable displayed=GizmoScatterMarkerDisplayed()
	String list = GizmoListIsoSurfaceWaves()+GizmoListScatterWaves(), mstr="", str
	Variable i,N=ItemsInLIst(list)
	for (i=0;i<N;i+=1)
		str = StringFromLIst(i,list)
		str = StringFromLIst(0,str,"=")
		if (strlen(str))
			mstr += str+";"
		endif
	endfor
	mstr = "\""+mstr+"\""
	PopupMenu waveSelectPopup,disable=0, value=#mstr

	STRUCT GizmoMarkerInfoStruct info
	info.win = "GizmoScatterMarkerPanel"
	Variable Ndisplay = GizmoMarkerInfo(info)	// get info and number displayed
	Variable MarkerNum = info.MarkerNum
	Variable dval = displayed ? 0 : 2				// 0=show, 2=disable
	dval = info.M[MarkerNum].used ? dval : 2
	dval = Ndisplay ? dval : 2						// always disabled when none displayed

	Button MarkerInfoButton,disable=dval
	SetVariable xyzStep,disable=dval, win=GizmoScatterMarkerPanel
	Button FitPeakButton,disable=(dval? 1 : 0)
	SetVariable alongWave,disable=dval, win=GizmoScatterMarkerPanel
	if (displayed)
		Button removeMarkerButton,disable=0, title="remove marker"
	else
		Button removeMarkerButton,disable=0, title="add marker"
	endif
	CheckBox OnPointCheckBox,disable=dval
	SetVariable MarkerSizeSetVar,disable=dval, win=GizmoScatterMarkerPanel
	Button ShowHideMarkerButton ,disable=(displayed ? 0 : 2)
	Button ShowHideMarkerButton, title=SelectString(info.M[MarkerNum].used,"Show","Hide")
	if (Ndisplay<1)
		return 0												// nothing more to do
	endif

	Wave marker=$GizmoGetScatterMarkerWave()
	Variable mx=NaN,my=NaN,mz=NaN
	if (WaveExists(marker))
		mx = marker[MarkerNum][0]
		my = marker[MarkerNum][1]
		mz = marker[MarkerNum][2]
	endif
	SetVariable Xabsolute,disable=dval,value= _NUM:mx, win=GizmoScatterMarkerPanel
	SetVariable Yabsolute,disable=dval,value= _NUM:my, win=GizmoScatterMarkerPanel
	SetVariable Zabsolute,disable=dval,value= _NUM:mz, win=GizmoScatterMarkerPanel
	Variable point = GizmoScatterMarkerGetNearPoint(mx,my,mz)
	SetVariable alongWave,value= _NUM:point, win=GizmoScatterMarkerPanel

	ControlInfo/W=GizmoScatterMarkerPanel waveSelectPopup
	Wave scatter=$S_Value
	updatePanelHKL(scatter,mx,my,mz)
	if (dval==0)
		updatePanelXYZlimits(scatter)
	endif

	String path=GetWavesDataFolder(scatter,1)
	Wave sizeW0=$(path+"gizmoScatterMarkerArraySize0")
	Variable size=sizeW0[MarkerNum]
	SetVariable MarkerSizeSetVar,value= _NUM:size, win=GizmoScatterMarkerPanel

	PopupMenu MarkerNumberPopup,mode=(MarkerNum+1),disable=0
	Wave RGBA=$(path+"gizmoScatterMarkerArrayRGBA")
	if (WaveExists(RGBA))
		Variable r,g,b
		r = RGBA[MarkerNum][0]*65535
		g = RGBA[MarkerNum][1]*65535
		b = RGBA[MarkerNum][2]*65535
		TitleBox colorTitleBox,labelBack=(r,g,b),win=GizmoScatterMarkerPanel
	endif

	if (WaveDims(scatter)!=3)
		Button FitPeakButton,disable=1		// only show this button when marking a 3D wave
	endif
	return displayed
End


Static Function updatePanelXYZlimits(scatter)
	Wave scatter

	if (!WaveExists(scatter))
		return 1
	endif
	Variable dx=1, dy=1, dz=1, xLo=-inf,yLo=-inf,zLo=-inf, xHi=Inf,yHi=Inf,zHi=Inf, Nx,Ny,Nz
	if (WaveDims(scatter)==3)
		dx = DimDelta(scatter,0)
		dy = DimDelta(scatter,1)
		dz = DimDelta(scatter,2)
		Nx=DimSize(scatter,0)
		Ny=DimSize(scatter,1)
		Nz=DimSize(scatter,2)
		xLo = DimOffset(scatter,0)
		yLo = DimOffset(scatter,1)
		zLo = DimOffset(scatter,2)
		xHi = xLo + dx*(Nx-1)
		yHi = yLo + dy*(Ny-1)
		zHi = zLo + dz*(Nz-1)
	endif

	ControlInfo/W=GizmoScatterMarkerPanel xyzStep
	Variable step = ((numtype(V_Value)==0 && V_Value>0) ? V_Value : NaN
	if (numtype(step)==0)					// a valid step was given, use it
		SetVariable Xabsolute, limits={xLo,xHi,step}, win=GizmoScatterMarkerPanel
		SetVariable Yabsolute, limits={yLo,yHi,step}, win=GizmoScatterMarkerPanel
		SetVariable Zabsolute, limits={zLo,zHi,step}, win=GizmoScatterMarkerPanel
	elseif (numtype(step) && WaveDims(scatter)==3)	// step is bad, but know 3D scaling
		SetVariable Xabsolute, limits={xLo,xHi,dx}, win=GizmoScatterMarkerPanel
		SetVariable Yabsolute, limits={yLo,yHi,dy}, win=GizmoScatterMarkerPanel
		SetVariable Zabsolute, limits={zLo,zHi,dz}, win=GizmoScatterMarkerPanel
	else
		step = 1									// give up, use a deault step of 1
		SetVariable Xabsolute, limits={xLo,xHi,step}, win=GizmoScatterMarkerPanel
		SetVariable Yabsolute, limits={yLo,yHi,step}, win=GizmoScatterMarkerPanel
		SetVariable Zabsolute, limits={zLo,zHi,step}, win=GizmoScatterMarkerPanel
	endif
End


Static Function updatePanelHKL(scatter,qx,qy,qz)
	Wave scatter
	Variable qx,qy,qz

	Variable bottom=340, calc=1
	if (!WaveExists(scatter) || numtype(qx+qy+qz))
		calc = 0
	elseif (WhichListItem(WaveUnits(scatter,0),Q_UNITS_LIST)<0)
		calc = 0
	endif
	Variable h=NaN,k=NaN,l=NaN
	if (calc)
		FUNCREF getRLfrom3DWaveProto getRL=$"getRLfrom3DWave"
		Wave RL = getRl(scatter,NaN)
		if (WaveExists(RL))
			Make/N=3/D/FREE Qc={qx,qy,qz}
			MatrixOP/FREE hkl = Inv(RL) x Qc
			h=hkl[0]
			k=hkl[1]
			l=hkl[2]
		endif
	endif

	GetWindow/Z GizmoScatterMarkerPanel wsize
	if (V_flag)
		return 1
	endif

	if (calc && numtype(h+k+l)==0)						// show hkl and intensity
		bottom += 80
		if (WaveDims(scatter)==3)
			Variable intens = scatter(qx)(qy)(qz)
			ValDisplay IntensDisp, disable=0, value= _NUM:intens
		else
			ValDisplay IntensDisp disable=1
		endif
		ValDisplay Hvaldisp,value= _NUM:h
		ValDisplay Kvaldisp,value= _NUM:k
		ValDisplay Lvaldisp,value= _NUM:l
		MoveWindow/W=GizmoScatterMarkerPanel V_left, V_top, V_right, V_top+bottom
	else
		MoveWindow/W=GizmoScatterMarkerPanel V_left, V_top, V_right, V_top+bottom
	endif
	return 0
End


Static Function/T GizmoAddScatterMarker([scatter,rgba,alpha])		// adds marker to gizmo if needed, will also create the marker wave if needed
	Wave scatter				// scatter wave that associated with the marker
	String rgba					// red, green, blue, or "" is black, or you can give your own rgba as "1,1,0,0.5"
	Variable alpha
	rgba = SelectString(ParamIsDefault(rgba),rgba,"")
	alpha = ParamIsDefault(alpha) ? NaN : alpha

	Execute "ModifyGizmo stopRotation"
	String wname=GizmoGetScatterMarkerWave()
	if (strlen(wname)<1)								// scatterMarker wave is not in object list, add it
		if (!WaveExists(scatter))						// need this to know where the scatterMarker waves are located
			return ""
		endif
		String path=GetWavesDataFolder(scatter,1)// path to the scatter marker waves
		Wave gizmoScatterMarkerArray=$(path+"gizmoScatterMarkerArray"), gizmoScatterMarkerArrayRGBA=$(path+"gizmoScatterMarkerArrayRGBA")
		Wave gizmoScatterMarkerArraySize=$(path+"gizmoScatterMarkerArraySize"), gizmoScatterMarkerArraySize0=$(path+"gizmoScatterMarkerArraySize0")
		if (!WaveExists(gizmoScatterMarkerArray))// old for compatibility
			Wave gizmoScatterMarkerArray=$(path+"gizmoScatterMarker")
		endif
		if (!WaveExists(gizmoScatterMarkerArray))
			Make/N=(8,3)/O $(path+"gizmoScatterMarkerArray")/WAVE=gizmoScatterMarkerArray =NaN
			Make/N=(8,3)/O $(path+"gizmoScatterMarkerArraySize")/WAVE=gizmoScatterMarkerArraySize =NaN
			Make/N=(8)/O $(path+"gizmoScatterMarkerArraySize0")/WAVE=gizmoScatterMarkerArraySize0 =0.5
			Make/N=(8,4)/O $(path+"gizmoScatterMarkerArrayRGBA")/WAVE=gizmoScatterMarkerArrayRGBA =NaN
			gizmoScatterMarkerArraySize[][0] = 0.5
			Note/K gizmoScatterMarkerArray,"waveClass=gizmoScatterMarkerXYZ"
			Note/K gizmoScatterMarkerArraySize,"waveClass=gizmoScatterMarkerSize"
			Note/K gizmoScatterMarkerArraySize0,"waveClass=gizmoScatterMarkerSize0"
			Note/K gizmoScatterMarkerArrayRGBA,"waveClass=gizmoScatterMarkerRGBA"
			gizmoScatterMarkerArrayRGBA[0][0]= {0,1,0,0,1,0,1,1}	// preset to "red;green;blue;magenta;cyan;yellow;white"
			gizmoScatterMarkerArrayRGBA[0][1]= {0,0,0.85,0,0,1,1,1}
			gizmoScatterMarkerArrayRGBA[0][2]= {0,0,0,1,1,1,0,1}
			gizmoScatterMarkerArrayRGBA[0][3]= {1,1,1,1,1,1,1,1}
		endif

		String CrossObject=""
		Execute "GetGizmo/Z objectItemExists=CrossPathGroup0"
		if (NumVarOrDefault("V_Flag",0))
			CrossObject = "CrossPathGroup0"		// use existing cross object
		else
			CrossObject = AddGizmoMarkerGroup("CrossPathGroup", rgba=rgba,alpha=alpha)
		endif

		Execute "ModifyGizmo startRecMacro"
		Execute "AppendToGizmo Scatter="+GetWavesDataFolder(gizmoScatterMarkerArray,2)+",name=scatterMarkerArray"
		Execute "ModifyGizmo ModifyObject=scatterMarkerArray property={ markerType,0}"
		Execute "ModifyGizmo ModifyObject=scatterMarkerArray property={ rotationType,0}"
		Execute "ModifyGizmo ModifyObject=scatterMarkerArray property={ Shape,7}"
		if (WaveExists(gizmoScatterMarkerArrayRGBA))
			Execute "ModifyGizmo ModifyObject=scatterMarkerArray property={ scatterColorType,1}"
			Execute "ModifyGizmo ModifyObject=scatterMarkerArray property={ colorWave,"+GetWavesDataFolder(gizmoScatterMarkerArrayRGBA,2)+"}"
		else
			Execute "ModifyGizmo ModifyObject=scatterMarkerArray property={ scatterColorType,0}"
			Execute "ModifyGizmo ModifyObject=scatterMarkerArray property={ color,0,0,0,.5}"
		endif
		if (WaveExists(gizmoScatterMarkerArraySize))
			Execute "ModifyGizmo ModifyObject=scatterMarkerArray property={ sizeType,1}"
			Execute "ModifyGizmo ModifyObject=scatterMarkerArray property={ sizeWave,"+GetWavesDataFolder(gizmoScatterMarkerArraySize,2)+"}"
		else
			Execute "ModifyGizmo ModifyObject=scatterMarkerArray property={ sizeType,0}"
			Execute "ModifyGizmo ModifyObject=scatterMarkerArray property={ size,0.5}"
		endif
		Execute "ModifyGizmo ModifyObject=scatterMarkerArray property={ objectName,"+CrossObject+"}"
		Execute "ModifyGizmo endRecMacro"
		wname = GetWavesDataFolder(gizmoScatterMarkerArray,2)
	endif

	Execute "GetGizmo displayItemExists=scatterMarkerArray"
	if (!NumVarOrDefault("V_Flag",1))		// object not in display list, add it before other objects
		Execute "GetGizmo displayList"		// list of all displayed objects
		Wave/T TW_DisplayList=TW_DisplayList
		Variable i,N=DimSize(TW_DisplayList,0)
		for (i=N-1;i>=0;i-=1)
			if (strsearch(TW_DisplayList[i],"ModifyGizmo setDisplayList="+num2istr(i)+", object",0)<0)
				do
					i += 1							// put marker after the "light" objects
				while(i<N && strsearch(TW_DisplayList[i],", object=light",0,2)>0)
				break
			endif
		endfor
		if (i>=0)									// position it before other objects
			Execute "ModifyGizmo insertDisplayList="+num2istr(i)+", object=scatterMarkerArray"
		else											// just put it at end
			Execute "ModifyGizmo setDisplayList=-1, object=scatterMarkerArray"
		endif
	endif

	return wname
End
//
Static Function GizmoRemoveScatterMarkers()		// removes marker from gizmo if needed, it does not remove it from the object list
	Execute "ModifyGizmo stopRotation"
	Execute "GetGizmo displayItemExists=scatterMarkerArray"
	NVAR V_flag=V_flag
	if (V_flag)										// object not in display list, remove it
		Execute "RemoveFromGizmo displayItem=scatterMarkerArray"
	endif
End
//
Static Function GizmoScatterMarkerDisplayed()		// returns true if scatter marker is displayed
	Execute "GetGizmo displayItemExists=scatterMarkerArray"
	NVAR V_Flag=V_Flag
	return V_Flag
End
//
Static Function/T GizmoGetScatterMarkerWave()	// gets name of wave with position of scatter marker
	String wname=""

	Execute "ModifyGizmo stopRotation"
	Execute "GetGizmo objectItemExists=scatterMarkerArray"
	NVAR V_Flag=V_Flag
	if (!V_Flag)									// not in object list, so stop here
		return ""
	endif

	Execute "GetGizmo objectList"			// find name of wave with marker position
	String str
	Wave/T TW_gizmoObjectList
	Variable i,i0,i1
	for (i=0;i<numpnts(TW_gizmoObjectList);i+=1)
		str = TW_gizmoObjectList[i]
		if (strsearch(str,"AppendToGizmo Scatter=",0)>=0 && strsearch(str,"name=scatterMarkerArray",0)>0)
			i0 = strsearch(str,"=",0)+1
			i1 = strsearch(str,",name=scatterMarkerArray",0)-1
			wname = str[i0,i1]
			break
		endif
	endfor
	KillWaves/Z TW_gizmoObjectList
	KIllStrings/Z S_gizmoObjectList

	return wname
End
//
Static Function GizmoScatterMarkerGetNearPoint(mx,my,mz)
	Variable mx, my, mz
	if (numtype(mx+my+mz))
		return NaN
	endif
	ControlInfo/W=GizmoScatterMarkerPanel waveSelectPopup
	if (V_flag!=3)
		return NaN
	endif
	Wave scatter=$S_Value
	if (!WaveExists(scatter))
		return NaN
	endif

	if (WaveDims(scatter)==2 && DimSize(scatter,1)==3)	// a triplet scatter wave
		Variable N=DimSize(scatter,0)
		Make/N=(N)/FREE delta
		delta = (scatter[p][0]-mx)^2 + (scatter[p][1]-my)^2 + (scatter[p][2]-mz)^2
		WaveStats/M=1/Q delta
		return V_minloc

	elseif (WaveDims(scatter)==3)								// a scaled 3d wave
		Variable ix,iy,iz
		Variable Nx=DimSize(scatter,0), Ny=DimSize(scatter,1), Nz=DimSize(scatter,2)
		ix = round( (mx-DimOffset(scatter,0)) / DimDelta(scatter,0) )
		iy = round( (my-DimOffset(scatter,1)) / DimDelta(scatter,1) )
		iz = round( (mz-DimOffset(scatter,2)) / DimDelta(scatter,2) )
		ix = limit(ix,0,DimSize(scatter,0)-1)
		iy = limit(iy,0,DimSize(scatter,1)-1)
		iz = limit(iz,0,DimSize(scatter,2)-1)
		return (iz*Nx*Ny) + iy*Nx + ix
	endif
	return NaN
End



Static Function/WAVE centerOf3Ddata(ww)	// finds center of data, works for triplets and for a 3D array
	Wave ww
	if (!WaveExists(ww))
		return $""
	endif

	if (WaveDims(ww)==2 && DimSize(ww,1)==3 && DimSize(ww,0)>0)	// triplets
		MatrixOP/FREE center = sumCols(ww)/numRows(ww)
		Redimension/N=3 center
	elseif (WaveDims(ww)==3 && numpnts(ww)>0)							// 3D array
		Make/N=3/D/FREE ccc
		ccc = DimOffset(ww,p) + DimDelta(ww,p)*DimSize(ww,p)/2
		Wave center = ccc
	else
		WAVE center = $""
	endif
	return center
End



Function/T FitPeakAt3Dmarker(space3D,Qc,QxHW,[QyHW,QzHW,printIt])
	Wave space3D
	Wave Qc								// center of sub-volume to fit
	Variable QxHW,QyHW,QzHW		// half widths dQz, dQy, dQz for the sub volume
	Variable printIt

	printIt = ParamIsDefault(printIt) ? 0 : !(!printIt)
	QyHW = ParamIsDefault(QyHW) ? NaN : QyHW	// QyHW & QzHW default to QwhX
	QzHW = ParamIsDefault(QzHW) ? NaN : QzHW
	if (WaveExists(space3D))
		if (WaveDims(space3D)!=3)
			return ""					// quit if really bad input passed
		endif
	endif
	if (WaveExists(Qc))
		if (numpnts(Qc)!=3)
			return ""					// quit if really bad input passed
		endif
	endif
	String spaceList = WaveListClass("GizmoXYZ;Qspace3D*","*","DIMS:3")
	if (!WaveExists(space3D))					// no space3D passed
		if (ItemsInList(spaceList)==1)		// only 1 choice, use it
			Wave space3D = $StringFromList(0,spaceList)
		elseif (ItemsInList(spaceList)<1)	// no acceptable choices, quit
			return ""
		endif
	endif
	String QcList=WaveList("*",";","DIMS:2,MAXROWS:1,MINROWS:1,MAXCOLS:3,MINCOLS:3" )+WaveList("*",";","DIMS:1,MAXROWS:3,MINROWS:3" )
	if (!WaveExists(Qc))						// no Qc passed
		if (ItemsInList(QcList)==1)			// only 1 choice, use it
			Wave Qc = $StringFromList(0,QcList)
		elseif (ItemsInList(QcList)<1)		// no acceptable choices, quit
			return ""
		endif
	endif
	if (numpnts(Qc)==3 && WaveType(Qc,2)==2)	// a free wave was passed, deal with it
		QcList = "_free_;"+QcList
	endif

	if (!WaveExists(space3D) || !WaveExists(Qc) || !(QxHW>0))
		printIt = 1
		String QcName=NameOfWave(Qc), spaceName = NameOfWave(space3D)
		QcName = SelectString(strlen(QcName),"gizmoScatterMarker",QcName)
		QxHW = QxHW>0 ? QxHW : 1
		Prompt spaceName,"3D space",popup,spaceList
		Prompt QcName,"Center of 3D space",popup,QcList
		Prompt QxHW,"X-HW in Volume to Use"
		Prompt QyHW,"Y-HW in Volume to Use (NaN defaults to X-HW)"
		Prompt QzHW,"Z-HW in Volume to Use (NaN defaults to X-HW)"
		DoPrompt "Fit 3D Peak",spaceName,QcName,QxHW,QyHW,QzHW
		if (V_flag)
			return ""
		endif
		printf "FitPeakAt3Dmarker(%s, %s, %g",spaceName,QcName,QxHW
		if (QyHW>0)
			printf ", QyHW=%g",QyHW
		endif
		if (QzHW>0)
			printf ", QzHW=%g",QzHW
		endif
		printf ")\r"
		Wave space3D=$spaceName
		if (!stringmatch(QcName,"_free_"))
			Wave Qc=$QcName						// this is for a free wave which does not have a name
		endif
	endif
	QyHW = QyHW>0 ? QyHW : QxHW
	QzHW = QzHW>0 ? QzHW : QxHW
	if (!WaveExists(space3D) || !WaveExists(Qc) || !(QxHW>0) || !(QyHW>0) || !(QzHW>0))
		return ""
	endif

	Make/N=3/D/FREE Qhw={QxHW,QyHW,QzHW}	// half widths dQz, dQy, dQz for the sub volume
	if (!WaveExists(space3D) || !WaveExists(Qc))
		return ""
	elseif (numtype(sum(Qhw)+sum(Qc)))
		return ""
	endif

	Make/N=3/D/FREE Qlo, Qhi, iLo,iHi,Np
	Qlo = Qc[p]-Qhw[p]
	Qhi = Qc[p]+Qhw[p]
	Np = DimSize(space3D,p)
	iLo = (Qlo[p]-DimOffSet(space3D,p))/DimDelta(space3D,p)
	iHi = (QHi[p]-DimOffSet(space3D,p))/DimDelta(space3D,p)
	iLo = limit(floor(iLo),0,Np)
	iHi = limit(ceil(iHi),0,Np)
	Np = iHi[p]-iLo[p]+1
	if (0 && printIt)
		printWave(Qlo,Name="Qlo")
		printWave(Qhi,Name="Qhi")
		printWave(iLo,Name="iLo")
		printWave(iHi,Name="iHi")
		printWave(Np,Name="Np")
	endif
	if (WaveMin(Np)<=0)
		return ""
	endif
	Make/N=(Np[0],Np[1],Np[2])/FREE/D sub3D, stdDev
	sub3D = space3D[iLo[0]+p][iLo[1]+q][iLo[2]+r]
	SetScale/P x, iLo[0]*DimDelta(space3D,0)+DimOffset(space3D,0), DimDelta(space3D,0),"",sub3D
	SetScale/P y, iLo[1]*DimDelta(space3D,1)+DimOffset(space3D,1), DimDelta(space3D,1),"",sub3D
	SetScale/P z, iLo[2]*DimDelta(space3D,2)+DimOffset(space3D,2), DimDelta(space3D,2),"",sub3D
	stdDev = sub3D					// this assumes that one count in detector == 1 photon

	// set the starting point for the fit
	Variable maxVal=WaveMax(sub3D),minVal=WaveMin(sub3D)
	Make/N=8/D/O W_coef = {NaN,NaN,Qc[0],Qhw[0]/4,Qc[1],Qhw[1]/4,Qc[2],Qhw[2]/4}
	if (maxVal>-minVal)			// a positive peak
		W_coef[1] = maxVal
		W_coef[0] = minVal==0 ? 1 : minVal
	else								// a negative peak
		W_coef[1] = minVal
		W_coef[0] = maxVal==0 ? 1 : maxVal
	endif

//	Variable V_FitOptions=2, V_FitError=0, V_FitQuitReason=0		// V_FitOptions=2 means robust fitting
	Variable V_FitError=0, V_FitQuitReason=0
	FuncFitMD/Q GMarkers#Gaussian3DFitFunc, W_coef, sub3D/W=stdDev/I=1
	Variable chisq = V_chisq
	Wave W_sigma=W_sigma
	Make/N=3/T/FREE units=WaveUnits(space3D,p)
	Make/N=3/D/FREE Qo=W_coef[2*p + 2]

	Variable err = V_FitError
	err = err || !( abs(Qc[0]-Qo[0]) < Qhw[0] )
	err = err || !( abs(Qc[1]-Qo[1]) < Qhw[0] )
	err = err || !( abs(Qc[2]-Qo[2]) < Qhw[0] )
	if (err)
		print "ERROR -- Fit failed, peak not inside of fitting volume"
		print FitErrorString(V_FitError,V_FitQuitReason)
		return ""
	endif

	Wave hkl=$""
#if exists("diffractometer#sample2crystal")
	STRUCT sampleStructure sa	
	String str = ParseFilePath(1,GetWavesDataFolder(space3D,1),":",1,0)+"sampleStructStr"
	String strStruct=StrVarOrDefault(str,"")				// fill the sample structure with values in spec scan directory
	StructGet/S/B=2 sa, strStruct								// found structure information, load into s
	Wave hkl = diffractometer#sample2crystal(sa,Qo)		// rotate qvec from sample-location frame into crystal based frame, the hkl
#endif
	if (!WaveExists(hkl))
		FUNCREF getRLfrom3DWaveProto getRL=$"getRLfrom3DWave"
		Wave RL = getRl(space3D,NaN)
		if (WaveExists(RL))
			MatrixOP/FREE hkl = Inv(RL) x Qo
		endif
	endif

	if (printIt)
		units = SelectString(stringmatch(units[p],"nm\\S-1*"),units[p],"1/nm")		// two different ways of writing 1/nm
		units = SelectString(strlen(units[p]),""," ("+units[p]+")")						// add () when something present
		Make/N=3/T/FREE Qstr=SelectString(strsearch(units[p],"1/nm",0)>=0,"","Q")	// is it Q?
		printf "Gaussian Fit has offset = %s,  amp = %s,  chisq = %g\r",ValErrStr(W_coef[0],W_sigma[0],sp=1),ValErrStr(W_coef[1],W_sigma[1],sp=1),chisq
		printf "  <%sx> = %s%s,   FWHMx = %s%s\r",Qstr[0],ValErrStr(Qo[0],W_sigma[2],sp=1),units[0],ValErrStr(W_coef[3],W_sigma[3],sp=1),units[0]
		printf "  <%sy> = %s%s,   FWHMy = %s%s\r",Qstr[1],ValErrStr(Qo[1],W_sigma[4],sp=1),units[1],ValErrStr(W_coef[5],W_sigma[5],sp=1),units[1]
		printf "  <%sz> = %s%s,   FWHMz = %s%s\r",Qstr[2],ValErrStr(Qo[2],W_sigma[6],sp=1),units[1],ValErrStr(W_coef[7],W_sigma[7],sp=1),units[2]
		if (stringmatch(units[0],units[1]) && stringmatch(units[0],units[2]))
			Variable Qmag=norm(Qo), dQmag=sqrt(W_sigma[2]^2 + W_sigma[4]^2 + W_sigma[6]^2)
			printf "  <|%s|> = %s%s\r",Qstr[0],ValErrStr(Qmag,dQmag,sp=1),units[0]
		endif
		if (WaveDims(space3D)==3)
			Np = DimSize(space3D,p)
			Make/N=3/D/FREE ijk
			ijk = (Qo[p]-DimOffset(space3D,p)) / DimDelta(space3D,p)
			ijk = limit(round(ijk[p]),0,Np[p]-1)
			Variable i=ijk[0], j=ijk[1], k=ijk[2]
			Variable N = k*Np[0]*Np[1] + j*Np[0] + i
			printf "  Closest point to peak is the %s[%g, %g, %g] or [%d] = %g\r",NameOfWave(space3D),i,j,k,N,space3D[i][j][k]
		endif
		if (WaveExists(hkl))
			printf "hkl Peak Center = %s\r",vec2str(hkl)
		endif
	endif

	String keyVals=""
	keyVals = ReplaceNumberByKey("offset",keyVals,W_coef[0],"=")
	keyVals = ReplaceNumberByKey("offsetErr",keyVals,W_sigma[0],"=")
	keyVals = ReplaceNumberByKey("amp",keyVals,W_coef[1],"=")
	keyVals = ReplaceNumberByKey("ampErr",keyVals,W_sigma[1],"=")
	keyVals = ReplaceNumberByKey("Xc",keyVals,W_coef[2],"=")
	keyVals = ReplaceNumberByKey("XcErr",keyVals,W_sigma[2],"=")
	keyVals = ReplaceNumberByKey("Yc",keyVals,W_coef[4],"=")
	keyVals = ReplaceNumberByKey("YcErr",keyVals,W_sigma[4],"=")
	keyVals = ReplaceNumberByKey("Zc",keyVals,W_coef[6],"=")
	keyVals = ReplaceNumberByKey("ZcErr",keyVals,W_sigma[6],"=")
	keyVals = ReplaceNumberByKey("FWX",keyVals,W_coef[3],"=")
	keyVals = ReplaceNumberByKey("FWXErr",keyVals,W_sigma[3],"=")
	keyVals = ReplaceNumberByKey("FWY",keyVals,W_coef[5],"=")
	keyVals = ReplaceNumberByKey("FWYErr",keyVals,W_sigma[5],"=")
	keyVals = ReplaceNumberByKey("FWZ",keyVals,W_coef[7],"=")
	keyVals = ReplaceNumberByKey("FWZErr",keyVals,W_sigma[7],"=")
	keyVals = ReplaceNumberByKey("chisq",keyVals,V_chisq,"=")
	if (strlen(units[0]))
		keyVals = ReplaceStringByKey("Xunit",keyVals,units[0],"=")
	endif
	if (strlen(units[1]))
		keyVals = ReplaceStringByKey("Yunit",keyVals,units[1],"=")
	endif
	if (strlen(units[2]))
		keyVals = ReplaceStringByKey("Zunit",keyVals,units[2],"=")
	endif
	if (WaveExists(hkl))
		keyVals = ReplaceStringByKey("hklPeakCenter",keyVals,vec2str(hkl,sep=","),"=")
	endif
	return keyVals
End
//
Static Function Gaussian3DFitFunc(w,xx,yy,zz) : FitFunc
	Wave w
	Variable xx
	Variable yy
	Variable zz

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ Variable ss = ((xx-x0)/xFWHM)^2 + ((yy-y0)/yFWHM)^2 + ((zz-z0)/zFWHM)^2
	//CurveFitDialog/ ss *= 4*ln(2)
	//CurveFitDialog/ ss = max(-500,ss)
	//CurveFitDialog/ f(xx,yy,zz) = offset + amp*exp(-ss)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 3
	//CurveFitDialog/ xx
	//CurveFitDialog/ yy
	//CurveFitDialog/ zz
	//CurveFitDialog/ Coefficients 8
	//CurveFitDialog/ w[0] = offset
	//CurveFitDialog/ w[1] = amp
	//CurveFitDialog/ w[2] = x0
	//CurveFitDialog/ w[3] = xFWHM
	//CurveFitDialog/ w[4] = y0
	//CurveFitDialog/ w[5] = yFWHM
	//CurveFitDialog/ w[6] = z0
	//CurveFitDialog/ w[7] = zFWHM

	Variable ss = ((xx-w[2])/w[3])^2 + ((yy-w[4])/w[5])^2 + ((zz-w[6])/w[7])^2
	ss *= 4*ln(2)
	ss = max(-500,ss)
	return w[0] + w[1]*exp(-ss)
End


Static Function GizmoMarkerInfo(info)		// returns number of markers displayed
	STRUCT GizmoMarkerInfoStruct &info

	String win=info.win
	ControlInfo/W=$win MarkerNumberPopup
	info.MarkerNum = V_Value-1

	Variable Ndisplay=0, i				// i must lie in range [0,7]
	for (i=0;i<NUMBER_OF_GIZMO_MARKERS;i+=1)
		Ndisplay += (GizmoMarkerInfo1(info.M[i],win,i)==0)
	endfor
	return Ndisplay
End
//
Static Function GizmoMarkerInfo1(M,win,MarkerNum)
	STRUCT GizmoMarkerInfoStruct1 &M
	String win
	Variable MarkerNum					// must lie in range [0, NUMBER_OF_GIZMO_MARKERS-1]
	MarkerNum = round(MarkerNum)
	Init_GizmoMarkerInfoStruct1(M)	// set all values for unsued

	if (strlen(win)<0)
		return 1
	elseif (!(MarkerNum>=0 && MarkerNum<=7))
		return 1
	endif

	String S_Value = ""
	ControlInfo/W=$win waveSelectPopup
	Wave scatter=$S_Value
	if (!WaveExists(scatter))
		return 1
	endif
	String path=GetWavesDataFolder(scatter,1)
	Wave wSize=$(path+"gizmoScatterMarkerArraySize")
	if (!WaveExists(wSize))
		return 1
	endif
	M.wName = GetWavesDataFolder(scatter,2)
	ControlInfo/W=$win MarkerNumberPopup
	Variable i = strsearch(S_recreation, "value= #\"",0)
	String list = S_recreation[i+11,Inf]
	i = strsearch(list, "\"",0)
	list = list[0,i-2]
	M.MarkerName = StringFromList(MarkerNum,list)
	Variable point = NaN

	Wave gizmoScatter=$(path+"gizmoScatterMarkerArray")
	if (WaveExists(gizmoScatter))
		Variable x0=gizmoScatter[MarkerNum][0], y0=gizmoScatter[MarkerNum][1], z0=gizmoScatter[MarkerNum][2]
		M.used = numtype(wSize[MarkerNum][0]+wSize[MarkerNum][1]+wSize[MarkerNum][2]+x0+y0+z0)==0
	else
		M.used = 0
		return 1
	endif

	if (M.used)
		String xUnit="", yUnit="", zUnit=""
		if (WaveDims(scatter)==3)
			xUnit=WaveUnits(scatter,0)
			yUnit=WaveUnits(scatter,1)
			zUnit=WaveUnits(scatter,2)
			if (strsearch(xUnit, "\S",0)>0)
				xUnit = ReplaceString("\S", xUnit,"^(")
				xUnit = ReplaceString("\M", xUnit,")")
			endif
			if (strsearch(yUnit, "\S",0)>0)
				yUnit = ReplaceString("\S", yUnit,"^(")
				yUnit = ReplaceString("\M", yUnit,")")
			endif
			if (strsearch(zUnit, "\S",0)>0)
				zUnit = ReplaceString("\S", zUnit,"^(")
				zUnit = ReplaceString("\M", zUnit,")")
			endif
		else
			xUnit = WaveUnits(scatter,-1)
			yUnit = xUnit
			zUnit = xUnit
		endif
		M.xUnit = xUnit;		M.yUnit = yUnit;		M.zUnit = zUnit

		M.x0 = x0;		M.y0 = y0;		M.z0 = z0
		point = GMarkers#GizmoScatterMarkerGetNearPoint(x0,y0,z0)
		M.point = point
		if (numtype(point))
			M.used = 0
			M.intensity = NaN
			M.ix = NaN	;	M.iy = NaN	;	M.iz = NaN
		else
			M.intensity = WaveDims(scatter)==3 ? scatter[point] : NaN
			if (WaveDims(scatter)==3)
				M.ix=round((x0-DimOffset(scatter,0))/DimDelta(scatter,0))
				M.iy=round((y0-DimOffset(scatter,1))/DimDelta(scatter,1))
				M.iz=round((z0-DimOffset(scatter,2))/DimDelta(scatter,2))
			endif
		endif
	endif

	FUNCREF getRLfrom3DWaveProto getRL=$"getRLfrom3DWave"
	Wave RL = getRl(scatter,point)
	if (WaveDims(scatter)==3 && WaveExists(RL) && WhichListItem(WaveUnits(scatter,0),Q_UNITS_LIST)>=0)
		Make/N=3/D/FREE Qc={x0,y0,z0}		// This is probably Q-space, so try for an hkl
		MatrixOP/FREE hkl = Inv(RL) x Qc
		M.h = hkl[0];		M.k = hkl[1];		M.l = hkl[2]
	endif
	if (WaveExists(RL))							// save the reciprocal lattice
		Redimension/N=9 RL
		for (i=0;i<9;i+=1)
			M.recip[i] = RL[i]
		endfor
	endif
	return 0
End

ThreadSafe Function/WAVE getRLfrom3DWaveProto(scatter,point)
	Wave scatter
	Variable point
	return $""
End


Static Function Init_GizmoMarkerInfoStruct1(M)		// zero out the structure, set everything to "Not Used"
	STRUCT GizmoMarkerInfoStruct1 &M
	M.used = 0
	M.MarkerName = ""
	M.wName = ""
	M.xUnit = "";	M.yUnit = "";	M.zUnit = ""
	M.x0 = NaN;	M.y0 = NaN;	M.z0 = NaN
	M.point = NaN
	M.intensity = NaN
	M.ix = NaN;		M.iy = NaN;		M.iz = NaN
	M.h = NaN;		M.k = NaN;		M.l = NaN
	M.recip[0] = NaN;		M.recip[1] = NaN;		M.recip[2] = NaN;
	M.recip[3] = NaN;		M.recip[4] = NaN;		M.recip[5] = NaN;
	M.recip[6] = NaN;		M.recip[7] = NaN;		M.recip[8] = NaN;
End


Function PrintGizmoMarkerInfoAll()
	STRUCT GizmoMarkerInfoStruct info
	info.win = "GizmoScatterMarkerPanel"
	if (GizmoMarkerInfo(info)<1)			// no markers dissplayed
		print "No Markers Displayed"
		return 0
	endif
	PrintGizmoMarkerInfoStruct(info)
End
//
Static Function PrintGizmoMarkerInfoStruct(info)
	STRUCT GizmoMarkerInfoStruct &info
	printf "Currently slected Marker %g  on %s\r",info.MarkerNum,info.win
	Variable x0,y0,z0, point, ix,iy,iz, intensity, h,k,l
	String xUnit,yUnit,zUnit
	Variable i
	for (i=0;i<NUMBER_OF_GIZMO_MARKERS;i+=1)
		if (info.M[i].used)
			print " "
			printf "for marker%d,  '%s' on wave '%s'%s\r",i,info.M[i].MarkerName,info.M[i].wName,SelectString(i==info.MarkerNum,"","\t\t**** Selected Marker ****")
			x0=info.M[i].x0 ; 	y0=info.M[i].y0 ;	 z0=info.M[i].z0
			point=info.M[i].point
			ix=info.M[i].ix;		iy=info.M[i].iy;		iz=info.M[i].iz
			xUnit=info.M[i].xUnit;	yUnit=info.M[i].yUnit;	zUnit=info.M[i].zUnit
			intensity=info.M[i].intensity
			h=info.M[i].h;	k=info.M[i].k;	l=info.M[i].l
			if (numtype(x0+y0+z0)==0)
				if (stringmatch(xUnit,yUnit) && stringmatch(xUnit,zUnit))
					printf "   at {%g, %g, %g}%s\r",x0,y0,z0,xUnit
				else
					printf "   at {%g (%s), %g (%s), %g (%s)}\r",x0,xUnit,y0,yUnit,z0,zUnit
				endif
				if (numtype(ix+iy+iz)==0)
					printf "   nearest actual point is [%g, %g, %g],   point#=%d,   intensity=%g\r",ix,iy,iz,point,intensity
				else
					printf "   nearest actual point is point#=%d,   intensity=%g\r",point,intensity
				endif
			endif
			if (numtype(h+k+l)==0)
				printf "   (hkl) = (%g, %g, %g)\r",h,k,l
			endif
		endif
	endfor
	return 0
End

#if strlen(WinList("LatticeSym.ipf","","WIN:128"))
//#if exists("InitLatticeSymPackage")==6
Static Function/T Compare2ReciprocalLatticesList(rl0,rl1)	// provide key=value list comparing two reciprocal lattices
	Wave rl0,rl1											// two 3x3 reciprocal lattices

	if (numtype(sum(rl0)+sum(rl1)))					// invalid reciprocal lattices
		return 	""
	elseif (EqualWaves(rl0,rl1,1))					// identical reciprocal lattices
		return 	"axisRotationAngle=0;axisHKLdeviation=0"
	endif

	MatrixOP/FREE rot01 = rl1 x Inv(rl0)			// rotation matrix from 0 to 1
	Make/N=3/D/FREE axis
	Variable angle = axisOfMatrix(rot01,axis)	// total rotation angle between point0 and point1 (degrees)

	MatrixOP/FREE hkl = Normalize(Inv(rl0) x axis)
	Variable h=round(24*hkl[0]), k=round(24*hkl[1]), l=round(24*hkl[2])
	lowestOrderHKL(h,k,l)
	Make/FREE hkl0 = {h,k,l}
	MatrixOP/FREE axis0 = Normalize(rl0 x hkl0)// axis of hkl0
	Variable dot = limit(MatrixDot(axis,axis0),-1,1)// ensure that the acos will exist
	Variable axisAnlge = acos(dot)*180/PI

	String list=""
	list = ReplaceStringByKey("axisDirectionXYZ",list,vec2str(axis,bare=1),"=")
	list = ReplaceNumberByKey("axisRotationAngle",list,angle,"=")
	list = ReplaceStringByKey("axisHKL",list,vec2str(hkl,bare=1),"=")
	list = ReplaceStringByKey("axisHKL0",list,hkl2str(h,k,l),"=")
	list = ReplaceNumberByKey("axisHKLdeviation",list,axisAnlge,"=")
	return list
End
#else
Static Function/T Compare2ReciprocalLatticesList(rl0,rl1)	
	Wave rl0,rl1
	return ""
End
#endif

//  ======================================================================================  //
//  =================================== End of Markers ===================================  //



//  ======================================================================================  //
//  =============================== Start of Movie Support ===============================  //

#if strlen(WinList("GizmoMovies.ipf","","WIN:128"))
//#if exists("AddGizmoMovieFrame_Proto")==6
Function AddGizmoMovieFrame_Markers(gm,data)
	STRUCT gizmoMovieStruct &gm		// needed to actually add frames to the movie
	String data								// string with data
	// data keys are:
	//		num			marker number [0,7]
	//		color			you can also use black, red, ... to identify marker
	//		action		"on" or "off"  on shows marker, off hides marker
	//							also, you can just use "show", "hide", "on", or "off" in place of action
	//		xyz			move marker to position given by "x y z"
	//		point			move marker to point number
	//		Nframes		add Nframes to the movie, if Nframes<1 or not present, then NO frames are added
	//
	// you can only work on one marker at a time
	//	structure of data		"num=1,action=on"  or  "num=1,action=off"  or  "num=1,xyz=2.2 -3.7 22"
	//			or  "num=1,point=135,Nframes=1"  or  "#=1,show"  or  "red,show"

	if (!GMarkers#GizmoScatterMarkerDisplayed())	// no markers, cannot do anything
		return 1
	endif
	String str

	Variable MarkerNum=NumberByKey("num",data,"=",",")				// marker number [0,7]
	if (numtype(MarkerNum))
		String colors="black;red;green;blue;magenta;cyan;yellow;white"
		Make/N=8/I/FREE iwave=WhichListItem(StringFromList(p,colors),data,",")
		WaveStats/Q/M=1 iwave
		MarkerNum = V_maxloc
	endif
	MarkerNum = MarkerNum==limit(MarkerNum,0,7) ? round(MarkerNum) : 0	// default to 0

	Variable show=-1
	String action=StringByKey("action",data,"=",",")
	if (strlen(action))
		if (WhichListItem(action,"on;show;true;1")>=0)
			show = 1
		elseif (WhichListItem(action,"off;hide;false;0")>=0)
			show = 0
		endif	
	else
		show = WhichListItem("hide",data,",")>=0 ? 0 : show
		show = WhichListItem("off",data,",")>=0 ? 0 : show
		show = WhichListItem("show",data,",")>=0 ? 1 : show
		show = WhichListItem("on",data,",")>=0 ? 1 : show
	endif

	Variable point=NumberByKey("point",data,"=",",")
	if (numtype(point))						// point is bad, check xyz
		str=StringByKey("xyz",data,"=",",")
		str = ReplaceString("\t",str," ")
		str = ReplaceString("  ",str," ")
		Make/N=3/D/FREE xyzw
		xyzw = str2num(StringFromList(p,str," "))
		if (numtype(sum(xyzw)))			// xyzw is invalid
			WaveClear xyzw
		endif
	endif
	Variable Nframes=NumberByKey("Nframes",data,"=",",")	// add some frames
	Nframes = Nframes>0 ? round(Nframes) : 0

	if (point>=0 || WaveExists(xyzw))	// need wave scatter
		STRUCT GizmoMarkerInfoStruct infoLocal
		infoLocal.win = "GizmoScatterMarkerPanel"
		GizmoMarkerInfo(infoLocal)
		Wave scatter = $(infoLocal.M[MarkerNum].wName)
		if (!WaveExists(scatter))			// cannot find scatter wave
			point = NaN
			WaveClear xyzw
		endif
	endif

	// move the marker
	if (point>=0)								// move marker to point number
		moveGizmoMarkerToPoint(point,scatter,MarkerNum)
	elseif (WaveExists(xyzw))				// or maybe move marker to xyz
		moveGizmoMarkerToXYZ(xyzw[0],xyzw[1],xyzw[2],scatter,MarkerNum)
	endif

	// show/hide the marker
	if (show>=0)								// show or hide the marker
		ShowHideGizmoMarker(MarkerNum,show)
	endif

	// Sleep/B/S 3								// give the gizmo a chance to update
	// add some frames to the movie
	if (Nframes>0)
		GizmoMovies#AddGizmoFrames2Movie(gm,Nframes)
	endif
	return 0
End
#endif

//  ================================ End of Movie Support ================================  //
//  ======================================================================================  //




//  ======================================================================================  //
//  =================================== Start of Init ====================================  //

Static Function InitGizmoMarkers()
	GizmoUtil#InitGizmoUtilityGeneral()
	Execute/Q/Z "GizmoMenu AppendItem={JZTm0,\"Marker Panel\", \"MakeGizmoScatterMarkerPanel()\"}"
	Execute/Q/Z "GizmoMenu AppendItem={JZTm1,\"   Marker Info\", \"PrintGizmoMarkerInfoAll()\"}"
End

//  ==================================== End of Init =====================================  //
//  ======================================================================================  //



