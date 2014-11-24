#pragma rtGlobals=2		// Use modern global access method.
#pragma ModuleName=GizmoUtil
#pragma IgorVersion = 6.20
#pragma version = 2.01
#include "ColorNames"

Static Constant GIZMO_MARKER_END_SIZE = 0.07		// puts boxes on ends of 3D marker (you can OverRide this in the Main procedure)
Static Constant GIZMO_MARKER_END_TYPE = 1			// 0=box on ends of lines, 1=sphere on ends of lines (you can OverRide this in the Main procedure)
Static Constant MAX_N_OBJECTS=500					// maximum number of objects that can be created using these routines (such as scatter, groups, titles, ...)
Static Constant GIZMO_SCALE_BAR_LEFT_EDGE = -1.9	// left edge of scale bar on a Gizmo

//	MakeGizmocubeCorners(xyz)	Make a scatter wave that contains xyz and is CUBICAL
//	AddGizmoTitleGroup()			Add the title to a gizmo. Returns name of item to include in the display list
//	setGizmoAxisLabels()			set the three axis labels
//
//	isGizmoObjectDisplayed()		returns -1 if object no found on displayList, if found return place in displayList
//	GetGizmoObjects()				returns list of objects of type "type", from top gizmo or named gizmo
//	GizmoListScatterWaves()		return list of all scatter plots, except 'scatterMarker0, scatterMarker1, ...
//	GizmoListIsoSurfaceWaves()	return list of all iso surface waves
//	GizmoListSurfaceWaves([gizmo])	return list of all Surface plots on gizmo

//Questionable:

//	AppendParametricXYZ2Gizmo()	Adds a parametric XYZ wave to the gizmo with its associated RGBA wave, returns name of object for parametricXYZ
//	GizmosWithWave()				Only used by AppendParametricXYZ2Gizmo(), Maybe this should not be here:




//  ======================================================================================  //
//  ============================ Start of Make Gizmo Cubical =============================  //

//	ActivateCornerCubeOnGizmo(), Show corner cubes on Gizmo, returns true if error displaying corner cubes
//	DeActivateCornerCubeOnGizmo(), de-activeate corner cubes on gizmo, does not delete the scatter object, just don't display it
//	FindMakeCubeCornerWaves(),  finds (or makes) the corner cube wave for a gizmo
//	isCornerCubeDisplayedOnGizmo(), returns True/False whether corener cubes are displayed
//	AddGizmoCornerCubesObject(), adds a scatter for corner cubes, returns name of scatter object (or existing scatter corner cubes object)
//	getCornerCubeObjectNameOnGizmo(), returns name of corner cube object on the gizmo, returns "" if no corner cube displayed

Function ActivateCornerCubeOnGizmo([GizmoName,forceCalc])			// Show corner cubes on Gizmo, returns true if error displaying corner cubes
	String GizmoName			// optional name of gizmo, defaults to top gizmo
	Variable forceCalc		// force recalculation of corners
	GizmoName = SelectString(ParamIsDefault(GizmoName),GizmoName,"")
	forceCalc = ParamIsDefault(forceCalc) || numtype(forceCalc) ? 0 : !(!forceCalc)

	String objectName=getCornerCubeObjectNameOnGizmo(GizmoName=GizmoName)	// if "", then add the corner cubes
	if (strlen(objectName)<1)
		Wave wCorners = FindMakeCubeCornerWaves(GizmoName=GizmoName, forceCalc=forceCalc)
		objectName = AddGizmoCornerCubesObject(wCorners,GizmoName=GizmoName)

	elseif (forceCalc)		// get name of corners wave and recalculate cube corners
		Wave wCorners = GizmoObjectWave(objectName,"Scatter")
		Wave ww = $(StringByKey("sourceWavePath",note(wCorners),"=")+StringByKey("sourceWave",note(wCorners),"="))
		Wave wCorners = MakeGizmocubeCorners(ww)
	endif
	if (strlen(objectName)<1)
		return 1
	endif

	if (!isCornerCubeDisplayedOnGizmo(GizmoName=GizmoName))
		String Nswitch=SelectString(strlen(GizmoName),"","/N="+GizmoName)
		Execute "ModifyGizmo"+Nswitch+" setDisplayList=-1, object="+objectName
		Execute "ModifyGizmo"+Nswitch+" compile"
	endif
	return 0
End


Function DeActivateCornerCubeOnGizmo([GizmoName])			// returns true if the gizmo corner cubes are un-displayed
	String GizmoName			// optional name of gizmo, defaults to top gizmo
	GizmoName = SelectString(ParamIsDefault(GizmoName),GizmoName,"")

	if (isCornerCubeDisplayedOnGizmo(GizmoName=GizmoName))
		String objectName=getCornerCubeObjectNameOnGizmo(GizmoName=GizmoName)	// if "", then add the corner cubes
		String Nswitch=SelectString(strlen(GizmoName),"","/N="+GizmoName)
		Execute "RemoveFromGizmo"+Nswitch+" displayItem="+objectName
		Execute "ModifyGizmo"+Nswitch+" compile"
	endif
	return 0
End


Static Function/WAVE FindMakeCubeCornerWaves([GizmoName,forceCalc])		// finds (or makes) the corner cube wave for a gizmo
	String GizmoName			// optional name of gizmo, defaults to top gizmo
	Variable forceCalc		// force recalculation of corners
	GizmoName = SelectString(ParamIsDefault(GizmoName),GizmoName,"")
	forceCalc = ParamIsDefault(forceCalc) || numtype(forceCalc) ? 0 : !(!forceCalc)

	String cornerList=WaveListClass("GizmoCorners","*","DIMS:2,MINROWS:8,MAXROWS:8,MINCOLS:3,MAXCOLS:3")
	String gizmoScatterList=GizmoListScatterWaves(gizmo=GizmoName)
	String gizmoIsoSurfaceList=GizmoListIsoSurfaceWaves(gizmo=GizmoName)		// get list of all iso surface waves
	String gizmoSurfaceList=GizmoListSurfaceWaves(gizmo=GizmoName)		// get list of all surface waves
	String gizmoListAll = gizmoScatterList+gizmoSurfaceList

	String str, list3Dobjects=""
	Variable i
	for (i=0;i<ItemsInList(gizmoListAll);i+=1)
		str = StringFromList(i, gizmoListAll)
		str = StringFromList(0,str,"=")
		Wave ww=$str
		if (WaveExists(ww))
			if (!WaveInClass(ww,"GizmoCorners;gizmoScatterMarkerXYZ"))		// make list of scatter waves that are NOT cube corners
				list3Dobjects += str + ";"
			endif
		endif
	endfor
	if (itemsInList(list3Dobjects)<1)					// no scatter waves, nothing to do
		return $""
	endif

	String list=""
	for (i=0;i<ItemsInList(list3Dobjects);i+=1)
		Wave ww = $StringFromList(i,list3Dobjects)
		list += WavesWithMatchingKeyVals(cornerList,"sourceWave="+NameOfWave(ww)+";")
	endfor

	i = ItemsInList(list)
	if (i<=1)
		str = StringFromList(0,list)						// only one choice, take it
	else
		Prompt str,"Cube Corners for Gizmo",popup,list
		DoPrompt "Cube Corners",str						// multiple choices, choose one
		if (V_flag)
			return $""
		endif
	endif
	Wave corners = $str

	if (!WaveExists(corners))								// no corner waves, make one
		if (ItemsInList(list3Dobjects)==1)
			Wave scat = $StringFromList(0,list3Dobjects)	// only one choice, take it
		else
			Prompt str,"3D object waves in Gizmo",popup,list3Dobjects
			DoPrompt "Scatter XYZ",str					// multiple choices, choose one
			if (V_flag)
				return $""
			endif
			Wave scat = $str
		endif
		Wave corners = MakeGizmocubeCorners(scat)

	elseif (forceCalc)
		Wave scat = $(StringByKey("sourceWavePath",note(corners),"=")+StringByKey("sourceWave",note(corners),"="))
		Wave corners = MakeGizmocubeCorners(scat)
	endif

	return corners
End


Static Function isCornerCubeDisplayedOnGizmo([GizmoName])		// returns true if the gizmo corner cubes are displayed
	String GizmoName			// optional name of gizmo, defaults to top gizmo
	GizmoName = SelectString(ParamIsDefault(GizmoName),GizmoName,"")
	String objectName=getCornerCubeObjectNameOnGizmo(GizmoName=GizmoName)	// probably "CubeCornersScatter"
	Variable i = isGizmoObjectDisplayed(objectName,gizmo=GizmoName)		// returns -1 if object not found on displayList, if found return place in displayList
	return (i>=0)
End


// Add a Corner Cubes Scatter object to a Gizmo. Returns name of objectName if successful.
Static Function/T AddGizmoCornerCubesObject(wCorners,[GizmoName])
	Wave wCorners				// wave with corner cube positions
	String GizmoName			// optional name of gizmo, defaults to top gizmo
	GizmoName = SelectString(ParamIsDefault(GizmoName),GizmoName,"")

	String Nswitch=""
	if (!WaveExists(wCorners))								// check that wCorners is valid
		return ""
	elseif(! (WaveDims(wCorners)==2 && DimSize(wCorners,0)==8 && DimSize(wCorners,1)==3) )
		return ""
	elseif (!ParamIsDefault(GizmoName) && strlen(GizmoName))
		if (WinType(GizmoName)!=13)
			return ""
		endif
		Nswitch = "/N="+GizmoName
	endif

	String objectName = getCornerCubeObjectNameOnGizmo(GizmoName=GizmoName)
	objectName = SelectString(strlen(objectName),"CubeCornersScatter",objectName)

	String scatterList=GetGizmoObjects("scatter")		// list of current scatter objects
	if (WhichListItem(objectName,scatterList)>=0)		// group already exists
		return objectName
	endif

	if (WhichListItem(objectName,scatterList)<0 && strlen(objectName))	// group does not exist, create it
		Execute "ModifyGizmo startRecMacro"
		Execute "AppendToGizmo"+Nswitch+" Scatter="+GetWavesDataFolder(wCorners,2)+",name="+objectName
		Execute "ModifyGizmo"+Nswitch+" ModifyObject="+objectName+" property={ scatterColorType,0}"
		Execute "ModifyGizmo"+Nswitch+" ModifyObject="+objectName+" property={ markerType,0}"
		Execute "ModifyGizmo"+Nswitch+" ModifyObject="+objectName+" property={ sizeType,0}"
		Execute "ModifyGizmo"+Nswitch+" ModifyObject="+objectName+" property={ rotationType,0}"
		Execute "ModifyGizmo"+Nswitch+" ModifyObject="+objectName+" property={ Shape,1}"
		Execute "ModifyGizmo"+Nswitch+" ModifyObject="+objectName+" property={ size,1}"
		Execute "ModifyGizmo"+Nswitch+" ModifyObject="+objectName+" property={ color,0,0,0,1}"
		Execute "ModifyGizmo"+Nswitch+" userString={CubeCorners,\""+objectName+"\"}"	// save name of cube corner object
		Execute "ModifyGizmo endRecMacro"
	endif
	return objectName
End


Static Function/T getCornerCubeObjectNameOnGizmo([GizmoName])		// returns name of corner cube object on gizmo, it may not yet be displayed
	String GizmoName			// optional name of gizmo, defaults to top gizmo

	GizmoName = SelectString(ParamIsDefault(GizmoName),GizmoName,"")
	String Nswitch=""
	if (!ParamIsDefault(GizmoName) && strlen(GizmoName))
		if (WinType(GizmoName)!=13)
			return ""
		endif
		Nswitch = "/N="+GizmoName
	endif

	Execute "GetGizmo/Z"+Nswitch+" userString=CubeCorners"
	String objectName=StrVarOrDefault("S_GizmoUserString","")
	KillStrings/Z S_GizmoUserString

	if (strlen(objectName)<1)			// in case of legacy gizmos
		String scatterList = GetGizmoObjects("scatter",gizmo=GizmoName)
		if (WhichListItem("CubeCornersScatter",scatterList)>=0)
			objectName = "CubeCornersScatter"
		elseif (WhichListItem("CubeCornersScatter",scatterList)>=0)
			objectName = "CubeCorners"
		endif
	endif

	objectName = SelectString(WhichListItem(objectName, GetGizmoObjects("scatter"))<0,objectName,"")	// check if object exists
	return objectName
End


// Make a scatter wave that contains xyz and is CUBICAL
Function/WAVE MakeGizmocubeCorners(xyz)
	Wave xyz							// xyz[N][3], or xyz[nx][ny][nz]
	if (!WaveExists(xyz))
		return $""
	endif

	Variable Xlo=NaN, Xhi=NaN, Ylo=NaN, Yhi=NaN, Zlo=NaN, Zhi=NaN
	String units
	if (WaveDims(xyz)==3)				// a scaled 3D array
		units = WaveUnits(xyz,0)
		Xlo = DimOffset(xyz,0)
		Xhi = DimDelta(xyz,0)*(DimSize(xyz,0)-1) + Xlo
		Ylo = DimOffset(xyz,1)
		Yhi = DimDelta(xyz,1)*(DimSize(xyz,1)-1) + Ylo
		Zlo = DimOffset(xyz,2)
		Zhi = DimDelta(xyz,2)*(DimSize(xyz,2)-1) + Zlo
		Variable swap
		if (Xlo>Xhi)
			swap = Xlo
			Xlo = Xhi
			Xhi = swap
		endif
		if (Ylo>Yhi)
			swap = Ylo
			Ylo = Yhi
			Yhi = swap
		endif
		if (Zlo>Zhi)
			swap = Zlo
			Zlo = Zhi
			Zhi = swap
		endif
	else
		Variable N=DimSize(xyz,0)		// a list of (xyz) triplets, xyz[N][3]
		units = WaveUnits(xyz,-1)
		Make/N=(N)/FREE/D cubeCorners_All
		cubeCorners_All = xyz[p][0]
		WaveStats/Q cubeCorners_All
		Xlo=V_min; Xhi=V_max
		cubeCorners_All = xyz[p][1]
		WaveStats/Q cubeCorners_All
		Ylo=V_min; Yhi=V_max
		cubeCorners_All = xyz[p][2]
		WaveStats/Q cubeCorners_All
		Zlo=V_min; Zhi=V_max
	endif

	Variable dX=Xhi-Xlo, dY=Yhi-Ylo, dZ=Zhi-Zlo
	Variable d = max(max(dX,dY),dZ)/2, mid

	Variable printIt = strlen(GetRTStackInfo(2))<=0
	if (printIt)
		printf "X = [%g, %g], 	Æ=%g\r",Xlo,Xhi,dX
		printf "Y = [%g, %g], 	Æ=%g\r",Ylo,Yhi,dY
		printf "Z = [%g, %g], 	Æ=%g\r",Zlo,Zhi,dZ
	endif
	mid = (Xlo+Xhi)/2;		Xlo = mid-d;	Xhi = mid+d
	mid = (Ylo+Yhi)/2;		Ylo = mid-d;	Yhi = mid+d
	mid = (Zlo+zhi)/2;		Zlo = mid-d;	Zhi = mid+d
	if (printIt)
		print " "
		printf "X = [%g, %g], 	Æ=%g\r",Xlo,Xhi,Xhi-Xlo
		printf "Y = [%g, %g], 	Æ=%g\r",Ylo,Yhi,Yhi-Ylo
		printf "Z = [%g, %g], 	Æ=%g\r",Zlo,Zhi,Zhi-Zlo
	endif

	String name=CleanupName(NameOfWave(xyz)+"Corners",0)
	Make/N=(8,3)/O $name=NaN
	Wave corners = $name
	SetScale d 0,0,units, corners
	corners[0][0] = Xlo;	corners[0][1] = Ylo; 	corners[0][2] = Zlo
	corners[1][0] = Xhi;	corners[1][1] = Ylo; 	corners[1][2] = Zlo
	corners[2][0] = Xlo;	corners[2][1] = Yhi; 	corners[2][2] = Zlo
	corners[3][0] = Xhi;	corners[3][1] = Yhi; 	corners[3][2] = Zlo
	corners[4][0] = Xlo;	corners[4][1] = Ylo; 	corners[4][2] = Zhi
	corners[5][0] = Xhi;	corners[5][1] = Ylo; 	corners[5][2] = Zhi
	corners[6][0] = Xlo;	corners[6][1] = Yhi; 	corners[6][2] = Zhi
	corners[7][0] = Xhi;	corners[7][1] = Yhi; 	corners[7][2] = Zhi
	String wnote="waveClass=GizmoCorners;"
	wnote = ReplaceStringByKey("sourceWave",wnote,NameOfWave(xyz),"=")
	wnote = ReplaceStringByKey("sourceWavePath",wnote,GetWavesDataFolder(xyz,1),"=")
	Note/K corners, wnote
	return corners
End

//  ============================= End of Make Gizmo Cubical ==============================  //
//  ======================================================================================  //



//  ======================================================================================  //
//  =========================== Start of Add to orModify Gizmo ===========================  //

// Add the title to a gizmo. Returns name of item to include in the display list
Function/T AddGizmoTitleGroup(groupName,title1,[title2,title3,title4,pos])
	String groupName			// probably "groupTitle", If this is empty, then a unique name will be assigned
	String title1
	String title2,title3,title4
	String pos
	title2 = SelectString(ParamIsDefault(title2),title2,"")
	title3 = SelectString(ParamIsDefault(title3),title3,"")
	title4 = SelectString(ParamIsDefault(title4),title4,"")
	pos = SelectString(ParamIsDefault(pos),pos,"TL")
	pos = SelectString(strlen(pos),"TL",pos)	// defaults to top left
	String font="Geneva"

	if (WhichListItem(pos,"TL;BL;")<0)
		return ""
	elseif (strlen(title1)<1)
		return ""
	elseif (CheckName(groupName,5))
		if (strlen(groupName)<1)						// create a unique group name
			NewDataFolder/O root:Packages				// ensure Packages exists
			NewDataFolder/O root:Packages:JZT_GizmoUtility	// ensure geometry exists
			if (exists("root:Packages:JZT_GizmoUtility:stringNumber")!=2)
				Variable/G root:Packages:JZT_GizmoUtility:stringNumber=-1
			endif
			NVAR stringNumber = root:Packages:JZT_GizmoUtility:stringNumber
			stringNumber = numtype(stringNumber) ? -1 : limit(round(stringNumber),-1,Inf)
			stringNumber += 1
			groupName = "gizmoStringGroup"+num2istr(stringNumber)
		endif
		groupName = CleanupName(groupName,0)
	endif

	Variable i
	for (i=0;i<3 && strlen(title2)==0;i+=1)	// don't leave empty strings between full ones
		title2 = title3							// e.g. "aa","","cc","dd"  --> "aa,"cc","dd",""
		title3 = title4
		title4 = ""
	endfor
	for (i=0;i<2 && strlen(title3)==0;i+=1)
		title3 = title4
		title4 = ""
	endfor

	// ************************* Group Object Start *******************
	Execute "AppendToGizmo group,name="+groupName
	Execute "ModifyGizmo currentGroupObject=\""+groupName+"\""
	Execute "AppendToGizmo string=\""+title1+"\",strFont=\""+font+"\",name=Title1"
	Execute "ModifyGizmo modifyObject=Title1 property={Clipped,0}"
	if (strlen(title2))
		Execute "AppendToGizmo string=\""+title2+"\",strFont=\""+font+"\",name=Title2"
		Execute "ModifyGizmo modifyObject=Title2 property={Clipped,0}"
	endif
	if (strlen(title3))
		Execute "AppendToGizmo string=\""+title3+"\",strFont=\""+font+"\",name=Title3"
		Execute "ModifyGizmo modifyObject=Title3 property={Clipped,0}"
	endif
	if (strlen(title4))
		Execute "AppendToGizmo string=\""+title4+"\",strFont=\""+font+"\",name=Title4"
		Execute "ModifyGizmo modifyObject=Title4 property={Clipped,0}"
	endif

	if (stringmatch(pos,"TL"))
		Execute "ModifyGizmo setDisplayList=0, opName=translateTitle, operation=translate, data={-1.9,1.9,0}"
	elseif (stringmatch(pos,"BL"))
		Execute "ModifyGizmo setDisplayList=0, opName=translateTitle, operation=translate, data={-1.9,-1.95,0}"
	endif
	Execute "ModifyGizmo setDisplayList=1, opName=rotateTitle, operation=rotate, data={180,1,0,0}"
	Execute "ModifyGizmo setDisplayList=2, opName=scaleTitle, operation=scale, data={0.1,0.1,0.1}"
	Execute "ModifyGizmo setDisplayList=3, object=Title1"
	if (strlen(title2))
		Execute "ModifyGizmo setDisplayList=4, opName=translateTitle2, operation=translate, data={0,1,0}"
		Execute "ModifyGizmo setDisplayList=5, opName=scaleTitle2, operation=scale, data={0.8,0.8,0.8}"
		Execute "ModifyGizmo setDisplayList=6, object=Title2"
		if (strlen(title3))
			Execute "ModifyGizmo setDisplayList=7, opName=translateTitle3, operation=translate, data={0,1,0}"
			Execute "ModifyGizmo setDisplayList=8, opName=scaleTitle3, operation=scale, data={0.7,0.7,0.7}"
			Execute "ModifyGizmo setDisplayList=9, object=Title3"
			if (strlen(title4))
				Execute "ModifyGizmo setDisplayList=10, opName=translateTitle4, operation=translate, data={0,1.5,0}"
				Execute "ModifyGizmo setDisplayList=11, object=Title4"
			endif
		endif
	endif
	Execute "ModifyGizmo currentGroupObject=\"::\""
	// ************************* Group Object End *******************
	return groupName
End


Function/T AddScaleBarGroup(groupName,maxLength,units,[scaleFactor,font])
	String groupName			// probably "ScaleBarGroup0", If this is empty, then a unique name will be assigned
	Variable maxLength		// maximum length in Gizmo, scale bar will be less than this
	String units				// units for maxLength and used in scale bar label
	Variable scaleFactor
	String font
	scaleFactor = ParamIsDefault(scaleFactor) ? 1 : scaleFactor
	scaleFactor = numtype(scaleFactor) || scaleFactor<=0 ? 1 : scaleFactor
	font = SelectString(ParamIsDefault(font),font,"Geneva")
	font = SelectString(strlen(font),"Geneva",font)
	if (maxLength<=0 || numtype(maxLength))
		return ""
	endif

	// for scale bar use multipliers of 1, 2, or 5 ONLY
	Variable BarLength = 10^floor(log(maxLength))
	if (5*BarLength < maxLength)
		BarLength = 5*BarLength
	elseif (2*BarLength < maxLength)
		BarLength = 2*BarLength
	endif
	String unitStr = num2str(BarLength)+" "+units
	Variable dxGizmoLine = BarLength/maxLength * scaleFactor * 2

	if (CheckName(groupName,5))						// invalid groupName passed, create one
		if (strlen(groupName)<1)						// create a unique group name
			NewDataFolder/O root:Packages			// ensure Packages exists
			NewDataFolder/O root:Packages:JZT_GizmoUtility	// ensure geometry exists
			if (exists("root:Packages:JZT_GizmoUtility:ScaleBarNumber")!=2)
				Variable/G root:Packages:JZT_GizmoUtility:ScaleBarNumber=-1
			endif
			NVAR ScaleBarNumber = root:Packages:JZT_GizmoUtility:ScaleBarNumber
			ScaleBarNumber = numtype(ScaleBarNumber) ? -1 : limit(round(ScaleBarNumber),-1,Inf)
			ScaleBarNumber += 1
			groupName = "ScaleBarGroup"+num2istr(ScaleBarNumber)
		endif
		groupName = CleanupName(groupName,0)
	endif
	if (CheckName(groupName,5))						// invalid groupName passed, give up
		return ""
	endif

	String cmd
	// ************************* Group Object Start *******************
	Execute "AppendToGizmo group,name="+groupName
	Execute "ModifyGizmo currentGroupObject=\""+groupName+"\""
	sprintf cmd "AppendToGizmo string=\"%s\",strFont=\"%s\",name=stringScaleBar",unitStr,font
	Execute cmd
	Execute "ModifyGizmo modifyObject=stringScaleBar property={Clipped,0}"
	sprintf cmd "AppendToGizmo line={%g,-1.9,0,%g,-1.9,0}, name=line0",GIZMO_SCALE_BAR_LEFT_EDGE,dxGizmoLine+GIZMO_SCALE_BAR_LEFT_EDGE
	Execute cmd
	Execute "AppendToGizmo attribute lineWidth=5, name=lineWidth0"
	Execute "ModifyGizmo setDisplayList=0, opName=pushMatrix0, operation=pushMatrix"
	Execute "ModifyGizmo setDisplayList=1, attribute=lineWidth0"
	Execute "ModifyGizmo setDisplayList=2, object=line0"
	Execute "ModifyGizmo setDisplayList=3, opName=translateScaleBarText, operation=translate, data={"+num2str(GIZMO_SCALE_BAR_LEFT_EDGE)+",-1.85,0}"
	Execute "ModifyGizmo setDisplayList=4, opName=scaleScaleBarText, operation=scale, data={0.1,0.1,0.1}"
	Execute "ModifyGizmo setDisplayList=5, opName=rotateScaleBarText, operation=rotate, data={180,1,0,0}"
	Execute "ModifyGizmo setDisplayList=6, object=stringScaleBar"
	Execute "ModifyGizmo setDisplayList=7, opName=popMatrix0, operation=popMatrix"
	Execute "ModifyGizmo currentGroupObject=\"::\""
	// ************************* Group Object End *******************
	return groupName
End
//
Static Function GzimoReSetScaleBarHookProc(s)
	STRUCT WMGizmoHookStruct &s
	if (!StringMatch(s.eventName,"scale"))
		return 0
	endif

	Variable scaleFactor=1
	String units = "nm"
	String win=s.winName
	Execute "GetGizmo/N="+win+"/Z objectNameList"
	String listObjectNames=StrVarOrDefault("S_ObjectNames",""), scaleBarName="", str, font="Geneva"
	KillStrings/Z S_ObjectNames
	Variable i
	for (i=0;i<ItemsInList(listObjectNames);i+=1)
		if (StringMatch(StringFromList(i,listObjectNames),"ScaleBarGroup*"))
			scaleBarName = StringFromList(i,listObjectNames)
			break
		endif
	endfor
	if (strlen(scaleBarName)<1)						// no scale bar, do nothing
		return 0
	endif

	Execute "GetGizmo/N="+win+"/Z displayItemExists=scaleBase"
	if (NumVarOrDefault("V_flag",0))
		Execute "GetGizmo/N="+win+"/Z displayList"
		Wave/T TW_DisplayList=TW_DisplayList
		for (i=0;i<numpnts(TW_DisplayList);i+=1)
			if (strsearch(TW_DisplayList[i],"opName=scaleBase, operation=scale,",0,2)>0)
				str = TW_DisplayList[i]
  					i = strsearch(str,"data={",0,2)
 					scaleFactor = str2num(str[i+6,Inf])
 					scaleFactor = numtype(scaleFactor) ? 1 : scaleFactor
				break
			endif

		endfor
		KillWaves/Z TW_DisplayList
		KillStrings/Z S_DisplayList
		KillVariables/Z V_Flag
	endif

	Execute "GetGizmo/N="+win+"/Z dataLimits"
	Variable dx = NumVarOrDefault("GizmoXmax",1)-NumVarOrDefault("GizmoXmin",1)
	Variable dy = NumVarOrDefault("GizmoYmax",1)-NumVarOrDefault("GizmoYmin",1)
	Variable dz = NumVarOrDefault("GizmoZmax",1)-NumVarOrDefault("GizmoZmin",1)
	KillVariables/Z GizmoXmin,GizmoXmax,GizmoYmin,GizmoYmax,GizmoZmin,GizmoZmax
	Variable maxLength = max(max(dx,dy),dz)
	Variable BarLength = 10^floor(log(maxLength))
	if (5*BarLength < maxLength)					// for scale bar use multipliers of 1, 2, or 5 ONLY
		BarLength = 5*BarLength
	elseif (2*BarLength < maxLength)
		BarLength = 2*BarLength
	endif
	Variable dxGizmoLine = BarLength/maxLength * scaleFactor * 2
	//	printf "maxLength = %g,  BarLength = %g,   dxGizmoLine = %g,  dxGizmoLine-1.75 = %g\r",maxLength,BarLength,dxGizmoLine, dxGizmoLine-1.75

	Execute "GetGizmo/N="+win+"/Z objectList"
	Wave/T TW_gizmoObjectList=TW_gizmoObjectList
	Variable i0,i1
	for (i0=0,i=0; i<numpnts(TW_gizmoObjectList); i+=1)
		if (strsearch(TW_gizmoObjectList[i], "AppendToGizmo group,name="+scaleBarName,0,2)>0)
			i0 = i
			break
		endif
	endfor
	for (i1=0,i=i0; i<numpnts(TW_gizmoObjectList); i+=1)
		if (strsearch(TW_gizmoObjectList[i], "ModifyGizmo currentGroupObject=\"::\"",0,2)>0)
			i1 = i
			break
		endif
	endfor

	for (i=i0;i<i1; i+=1)
		if (strsearch(TW_gizmoObjectList[i], "name=stringScaleBar",0,2)>0)
			str = TW_gizmoObjectList[i]
			i = strsearch(str, "AppendToGizmo string=\"",0,2)
			if (i<0)
				break
			endif
			str = str[i+22,Inf]
			i = strsearch(str,"\",",0)
			if (i<0)
				break
			endif
			str = str[0,i-1]
			i = strsearch(str," ",Inf,1)
			if (i<0)
				break
			endif
			str = str[i+1,Inf]
			units = str
			break
		endif
	endfor
	//	printf "units = '%s'\r",units
	KillWaves/Z TW_gizmoObjectList
	KillStrings/Z S_gizmoObjectList
	String unitStr = num2str(BarLength)+" "+units

	Execute "ModifyGizmo/N="+win+"/Z startRecMacro"
	Execute "ModifyGizmo/N="+win+"/Z currentGroupObject=\""+scaleBarName+"\""
	Execute "ModifyGizmo/N="+win+"/Z modifyObject=stringScaleBar, property={string,\""+unitStr+"\"}"
	sprintf str "ModifyGizmo/N=%s/Z modifyObject=line0, property={vertex,%g,-1.9,0,%g,-1.9,0}", win,GIZMO_SCALE_BAR_LEFT_EDGE,dxGizmoLine+GIZMO_SCALE_BAR_LEFT_EDGE
	Execute str
	Execute "ModifyGizmo/N="+win+"/Z currentGroupObject=\"::\""
	Execute "ModifyGizmo/N="+win+"/Z endRecMacro"
	return 0	 
End



// Add a marker to a gizmo. Returns name of item to include in the display list.
// Note this just creates a marker group, the marker group needs to be part of a scatter plot to display it
Function/T AddGizmoMarkerGroup(groupName,[rgba,alpha])
	String groupName			// probably "groupTitle", If this is empty, then a unique name will be assigned
	String rgba					// red, green, blue, or "" is black, or you can give your own rgba as "1,1,0,0.5"
	Variable alpha
	rgba = SelectString(ParamIsDefault(rgba),rgba,"")
	alpha = ParamIsDefault(alpha) ? 1 : alpha
	alpha = numtype(alpha) ? 1 : limit(alpha,0,1)


	if (CheckName(groupName,5))						// invalid groupName passed, create one
		if (strlen(groupName)<1)						// create a unique group name
			NewDataFolder/O root:Packages				// ensure Packages exists
			NewDataFolder/O root:Packages:JZT_GizmoUtility	// ensure geometry exists
			if (exists("root:Packages:JZT_GizmoUtility:MarkerNumber")!=2)
				Variable/G root:Packages:JZT_GizmoUtility:MarkerNumber=-1
			endif
			NVAR MarkerNumber = root:Packages:JZT_GizmoUtility:MarkerNumber
			MarkerNumber = numtype(MarkerNumber) ? -1 : limit(round(MarkerNumber),-1,Inf)
			MarkerNumber += 1
			groupName = "CrossPathGroup"+num2istr(MarkerNumber)
		endif
		groupName = CleanupName(groupName,0)
	endif
	if (CheckName(groupName,5))						// invalid groupName passed, create one
		return ""
	endif
	Variable r,g,b,a
	sscanf rgba,"%g,%g,%g,%g", r,g,b,a
	if (V_flag!=4)
		Wave rgb=color2RGB(rgba,1)					// a wide range of color names is available
		if (rgb[0]>0 && rgb[1]>0 && rgb[2]>0 && alpha>=0)
			sprintf rgba,"%g,%g,%g,%g", rgb[0],rgb[1],rgb[2],alpha
		else
			rgba = ""
		endif
	endif

	// ************************* Group Object Start *******************
	Execute "AppendToGizmo group,name="+groupName
	Execute "ModifyGizmo currentGroupObject=\""+groupName+"\""
	String str
	Variable v = GIZMO_MARKER_END_SIZE>0 ? 1-2*GIZMO_MARKER_END_SIZE : 1
	sprintf str, "AppendToGizmo line={%g,0,0,%g,0,0}, name=lineX",-v,v	;	Execute str
	sprintf str, "AppendToGizmo line={0,%g,0,0,%g,0}, name=lineY",-v,v	;	Execute str
	sprintf str, "AppendToGizmo line={0,0,%g,0,0,%g}, name=lineZ",-v,v	;	Execute str
	if (GIZMO_MARKER_END_SIZE>0)
		Variable endType=limit(GIZMO_MARKER_END_TYPE,0,1)
		if (endType==0)							// boxes
			sprintf str, "AppendToGizmo box={%g,%g,%g},name=endMarkers",GIZMO_MARKER_END_SIZE,GIZMO_MARKER_END_SIZE,GIZMO_MARKER_END_SIZE
			Execute str
			Execute "ModifyGizmo modifyObject=endMarkers property={calcNormals,1}"
		elseif (endType==1)						// spheres
			sprintf str, "AppendToGizmo sphere={%g,25,25},name=endMarkers",GIZMO_MARKER_END_SIZE
			Execute str
		endif
	endif
	Execute "AppendToGizmo attribute lineWidth=2, name=lineWidthCross"
	if (strlen(rgba))
		Execute "AppendToGizmo attribute color={"+rgba+"},name=colorLine"
	endif
	Execute "ModifyGizmo setDisplayList=0, opName=pushAttribute0, operation=pushAttribute, data=5"
	Execute "ModifyGizmo setDisplayList=1, attribute=lineWidthCross"
	if (strlen(rgba))
		Execute "ModifyGizmo setDisplayList=2, attribute=colorLine"
	endif
	Execute "ModifyGizmo setDisplayList=-1, object=lineX"
	Execute "ModifyGizmo setDisplayList=-1, object=lineY"
	Execute "ModifyGizmo setDisplayList=-1, object=lineZ"
	if (GIZMO_MARKER_END_SIZE>0)
		v = 1 - GIZMO_MARKER_END_SIZE
		sprintf str, "ModifyGizmo setDisplayList=-1, opName=translateXminus, operation=translate, data={%g,0,0}",-v	;	Execute str
		Execute "ModifyGizmo setDisplayList=-1, object=endMarkers"
		sprintf str, "ModifyGizmo setDisplayList=-1, opName=translateXplus, operation=translate, data={%g,0,0}",2*v	;	Execute str
		Execute "ModifyGizmo setDisplayList=-1, object=endMarkers"
		sprintf str, "ModifyGizmo setDisplayList=-1, opName=translateYminus, operation=translate, data={%g,%g,0}",-v,-v ; Execute str
		Execute "ModifyGizmo setDisplayList=-1, object=endMarkers"
		sprintf str, "ModifyGizmo setDisplayList=-1, opName=translateYplus, operation=translate, data={0,%g,0}",2*v	;	Execute str
		Execute "ModifyGizmo setDisplayList=-1, object=endMarkers"
		sprintf str, "ModifyGizmo setDisplayList=-1, opName=translateZminus, operation=translate, data={0,%g,%g}",-v,-v ; Execute str
		Execute "ModifyGizmo setDisplayList=-1, object=endMarkers"
		sprintf str, "ModifyGizmo setDisplayList=-1, opName=translateZplus, operation=translate, data={0,0,%g}",2*v	;	Execute str
		Execute "ModifyGizmo setDisplayList=-1, object=endMarkers"
	endif
	Execute "ModifyGizmo setDisplayList=-1, opName=popAttribute0, operation=popAttribute"
	Execute "ModifyGizmo currentGroupObject=\"::\""
	// ************************* Group Object End *******************

	return groupName
End


Function/T AddGizmoBeamLineAxesGroup(groupName)
	String groupName			// probably "BeamLineAxesGroup0", If this is empty, then a unique name will be assigned
	String font="Geneva"

	if (CheckName(groupName,5))						// invalid groupName passed, create one
		if (strlen(groupName)<1)						// create a unique group name
			NewDataFolder/O root:Packages				// ensure Packages exists
			NewDataFolder/O root:Packages:JZT_GizmoUtility	// ensure geometry exists
			if (exists("root:Packages:JZT_GizmoUtility:BeamLineAxesNumber")!=2)
				Variable/G root:Packages:JZT_GizmoUtility:BeamLineAxesNumber=-1
			endif
			NVAR BeamLineAxesNumber = root:Packages:JZT_GizmoUtility:BeamLineAxesNumber
			BeamLineAxesNumber = numtype(BeamLineAxesNumber) ? -1 : limit(round(BeamLineAxesNumber),-1,Inf)
			BeamLineAxesNumber += 1
			groupName = "BeamLineAxesGroup"+num2istr(BeamLineAxesNumber)
		endif
		groupName = CleanupName(groupName,0)
	endif
	if (CheckName(groupName,5))						// invalid groupName passed, give up
		return ""
	endif

	// ************************* Group Object Start *******************
	Execute "AppendToGizmo group,name="+groupName
	Execute "ModifyGizmo currentGroupObject=\""+groupName+"\""
	Execute "AppendToGizmo freeAxesCue={0,0,0,1},name=freeAxesCue_BeamLine"
	Execute "AppendToGizmo string=\"X(BL)\",strFont=\""+font+"\",name=stringX"
	Execute "ModifyGizmo modifyObject=stringX property={Clipped,0}"
	Execute "AppendToGizmo string=\"Y(BL)\",strFont=\""+font+"\",name=stringY"
	Execute "ModifyGizmo modifyObject=stringY property={Clipped,0}"
	Execute "AppendToGizmo string=\"Z(BL)\",strFont=\""+font+"\",name=stringZ"
	Execute "ModifyGizmo modifyObject=stringZ property={Clipped,0}"
	Execute "AppendToGizmo cylinder={0,0.2,1,25,25},name=cylinderArrow"
	Execute "AppendToGizmo attribute lineWidth=2, name=lineWidthCross"
	Execute "AppendToGizmo attribute color={1,0,0,1},name=colorXred"
	Execute "AppendToGizmo attribute color={0,1,0,1},name=colorYgreen"
	Execute "AppendToGizmo attribute color={0,0,1,1},name=colorZblue"
	Execute "ModifyGizmo setDisplayList=0, opName=pushAttribute0, operation=pushAttribute, data=5"
	Execute "ModifyGizmo setDisplayList=1, attribute=lineWidthCross"
	Execute "ModifyGizmo setDisplayList=2, opName=pushMatrix0, operation=pushMatrix"
	Execute "ModifyGizmo setDisplayList=3, opName=rotate45, operation=rotate, data={-45,1,0,0}"
	Execute "ModifyGizmo setDisplayList=4, object=freeAxesCue_BeamLine"
	Execute "ModifyGizmo setDisplayList=5, opName=pushMatrix_X, operation=pushMatrix"
	Execute "ModifyGizmo setDisplayList=6, opName=translateX, operation=translate, data={1,0,0}"
	Execute "ModifyGizmo setDisplayList=7, opName=scaleX, operation=scale, data={0.07,0.07,0.07}"
	Execute "ModifyGizmo setDisplayList=8, attribute=colorXred"
	Execute "ModifyGizmo setDisplayList=9, opName=rotateX, operation=rotate, data={-90,0,1,0}"
	Execute "ModifyGizmo setDisplayList=10, object=stringX"
	Execute "ModifyGizmo setDisplayList=11, object=cylinderArrow"
	Execute "ModifyGizmo setDisplayList=12, opName=popMatrix_X, operation=popMatrix"
	Execute "ModifyGizmo setDisplayList=13, opName=pushMatrix_Y, operation=pushMatrix"
	Execute "ModifyGizmo setDisplayList=14, opName=translateY, operation=translate, data={0,1,0}"
	Execute "ModifyGizmo setDisplayList=15, attribute=colorYgreen"
	Execute "ModifyGizmo setDisplayList=16, opName=scaleY, operation=scale, data={0.07,0.07,0.07}"
	Execute "ModifyGizmo setDisplayList=17, opName=rotateY, operation=rotate, data={90,1,0,0}"
	Execute "ModifyGizmo setDisplayList=18, object=stringY"
	Execute "ModifyGizmo setDisplayList=19, object=cylinderArrow"
	Execute "ModifyGizmo setDisplayList=20, opName=popMatrix_Y, operation=popMatrix"
	Execute "ModifyGizmo setDisplayList=21, opName=pushMatrix_Z, operation=pushMatrix"
	Execute "ModifyGizmo setDisplayList=22, opName=translateZ, operation=translate, data={0,0,1}"
	Execute "ModifyGizmo setDisplayList=23, attribute=colorZblue"
	Execute "ModifyGizmo setDisplayList=24, opName=scaleZ, operation=scale, data={0.07,0.07,0.07}"
	Execute "ModifyGizmo setDisplayList=25, opName=rotateZ, operation=rotate, data={180,1,0,0}"
	Execute "ModifyGizmo setDisplayList=26, object=stringZ"
	Execute "ModifyGizmo setDisplayList=27, object=cylinderArrow"
	Execute "ModifyGizmo setDisplayList=28, opName=popMatrix_Z, operation=popMatrix"
	Execute "ModifyGizmo setDisplayList=29, opName=popMatrix0, operation=popMatrix"
	Execute "ModifyGizmo setDisplayList=30, opName=popAttribute0, operation=popAttribute"
	Execute "ModifyGizmo currentGroupObject=\"::\""
	// ************************* Group Object End *******************

	return groupName
End


Function setGizmoAxisLabels(xlabel,ylabel,zlabel)
	String xlabel,ylabel,zlabel
	if (itemsInList(WinList("*",";","WIN:4096"))==0)
		return 1
	endif

	Execute "GetGizmo objectItemExists=axes0"
	NVAR V_Flag=V_Flag
	if (!V_flag)
		Execute "AppendToGizmo Axes=boxAxes,name=axes0"
	endif
	String cmd

	if (stringmatch(xlabel,"empty"))
		Execute "ModifyGizmo ModifyObject=axes0,property={0,axisLabel,0}"
	elseif (strlen(xlabel))
		Execute "ModifyGizmo ModifyObject=axes0,property={0,ticks,3}"
		Execute "ModifyGizmo ModifyObject=axes0,property={0,axisLabel,1}"
		sprintf cmd, "ModifyGizmo ModifyObject=axes0,property={0,axisLabelText,\"%s\"}",xlabel
		Execute cmd
		Execute "ModifyGizmo ModifyObject=axes0,property={0,axisLabelDistance,0.05}"
		Execute "ModifyGizmo ModifyObject=axes0,property={0,axisLabelScale,0.5}"
	endif

	if (stringmatch(ylabel,"empty"))
		Execute "ModifyGizmo ModifyObject=axes0,property={1,axisLabel,0}"
	elseif (strlen(ylabel))
		Execute "ModifyGizmo ModifyObject=axes0,property={1,ticks,3}"
		Execute "ModifyGizmo ModifyObject=axes0,property={1,axisLabel,1}"
		sprintf cmd, "ModifyGizmo ModifyObject=axes0,property={1,axisLabelText,\"%s\"}",ylabel
		Execute cmd
		Execute "ModifyGizmo ModifyObject=axes0,property={1,axisLabelDistance,0.05}"
		Execute "ModifyGizmo ModifyObject=axes0,property={1,axisLabelScale,0.5}"
	endif

	if (stringmatch(zlabel,"empty"))
		Execute "ModifyGizmo ModifyObject=axes0,property={2,axisLabel,0}"
	elseif (strlen(zlabel))
		Execute "ModifyGizmo ModifyObject=axes0,property={2,ticks,3}"
		Execute "ModifyGizmo ModifyObject=axes0,property={2,axisLabel,1}"
		sprintf cmd, "ModifyGizmo ModifyObject=axes0,property={2,axisLabelText,\"%s\"}",zlabel
		Execute cmd
		Execute "ModifyGizmo ModifyObject=axes0,property={2,axisLabelDistance,0.2}"
		Execute "ModifyGizmo ModifyObject=axes0,property={2,axisLabelScale,0.5}"
	endif
End


//	for clipPlane,  equation is:		A*x + B*y + C*z = D
//	ModifyGizmo setDisplayList=12, opName=ClipPlane0, operation=ClipPlane, data={i,A,B,C,D}

// Add a set of clip planes to a gizmo. Returns name of groupName if successful.
Function/T AddGizmoClipPlaneGroup(groupName)
	String groupName			// probably "groupClipPlane", If this is empty, then a unique name will be assigned

	String groupList=GetGizmoObjects("group")					// list of current group objects
	Variable groupExists=(WhichListItem(groupName,groupList)>=0)	// group already exists

	groupName = SelectString(strlen(groupName),"",CleanupName(groupName,0))
	if (CheckName(groupName,5))
		if (strlen(groupName)<1)						// create a unique group name
			Variable i
			for (i=0;i<MAX_N_OBJECTS;i+=1)			// at most MAX_N_OBJECTS groups
				groupName = "gizmoClipPlaneGroup"+num2istr(i)
				if (WhichListItem(groupName,groupList)<0)
					break
				else
					groupName = ""
				endif
			endfor
		endif
	endif
	if (WhichListItem(groupName,groupList)<0 && strlen(groupName))	// group does not exist, create it
		// ************************* Group Object Start *******************
		Execute "AppendToGizmo group,name="+groupName
		Execute "ModifyGizmo currentGroupObject=\""+groupName+"\""
		Execute "ModifyGizmo setDisplayList=0, opName=enable0, operation=enable, data=12288"	// enable all of the 6 clip planes
		Execute "ModifyGizmo setDisplayList=1, opName=enable1, operation=enable, data=12289"
		Execute "ModifyGizmo setDisplayList=2, opName=enable2, operation=enable, data=12290"
		Execute "ModifyGizmo setDisplayList=3, opName=enable3, operation=enable, data=12291"
		Execute "ModifyGizmo setDisplayList=4, opName=enable4, operation=enable, data=12292"
		Execute "ModifyGizmo setDisplayList=5, opName=enable5, operation=enable, data=12293"
		Execute "ModifyGizmo currentGroupObject=\"::\""
		// ************************* Group Object End *******************
		Execute "ModifyGizmo userString={"+groupName+",\"\"}"	// there are no clip planes yet defined, so make empty
	endif
	return groupName
End

// Modify a clip plane in a gizmo. Returns 1 on error, uses clipVal=NaN to remove a clip plane
Function ModifyGizmoClipPlaneGroup(groupName,action,clipVal)
	String groupName			// probably "groupClipPlane", If this is empty, then a unique name will be assigned
	String action				// "add;delete;remove;+X;-X;+Y;-Y;+Z;-Z;update"
	Variable clipVal				// value of clip plane, only used when adding planes (if valid when removeX, then can be used to select which Xplane gets removed)

	String groupList=GetGizmoObjects("group")					// list of current group objects
	if (WhichListItem(groupName,groupList)<0)				// group must already exists 
		return 1
	endif
	Execute "GetGizmo/Z userString="+groupName
	String keyVals=StrVarOrDefault("S_GizmoUserString","")
	KillStrings/Z S_GizmoUserString

	Variable err=0
	strswitch(action)
		case "delete":
		case "remove":
			Execute "RemoveFromGizmo/Z object="+groupName
			Execute "ModifyGizmo userString={"+groupName+",\"\"}"	// there are no clip planes, so make empty
			break

		case "update":
			String item
			Variable val, i
			for (i=0;i<ItemsInlist(keyVals);i+=1)
				item = StringFromList(i,keyVals)
				ModifyGizmoClipPlaneGroup(groupName,StringFromList(0,item,":"),str2num(StringFromList(1,item,":")))
			endfor
			break

		case "+X":													// change, add, or remove +X clip plane
		case "-X":
		case "+Y":
		case "-Y":
		case "+Z":
		case "-Z":
			action = UpperStr(action)
			String clipPlane="ClipPlane"+action[1]+SelectString(char2num(action[0])==43,"n","p")
			if (numtype(clipVal))
				keyVals = RemoveByKey(action,keyVals)
				Execute "ModifyGizmo currentGroupObject=\""+groupName+"\""
				Execute "RemoveFromGizmo/Z displayItem="+clipPlane
				Execute "ModifyGizmo currentGroupObject=\"::\""
			else
				keyVals = ReplaceNumberByKey(action,keyVals,clipVal)
				String box=GetGizmoBoxDisplayed("")
				Make/N=6/D/FREE xyzLoHi = str2num(StringFromList(p,box))
				Variable clipBox
				if (strsearch(action,"X",0,2)==1)							// convert X clipVal to range [-1,1]
					clipBox = (2*clipVal - xyzLoHi[0] - xyzLoHi[1]) / (xyzLoHi[1]-xyzLoHi[0])
				elseif (strsearch(action,"Y",0,2)==1)						// convert Y clipVal to range [-1,1]
					clipBox = (2*clipVal - xyzLoHi[2] - xyzLoHi[3]) / (xyzLoHi[3]-xyzLoHi[2])
				elseif (strsearch(action,"Z",0,2)==1)						// convert Z clipVal to range [-1,1]
					clipBox = (2*clipVal - xyzLoHi[4] - xyzLoHi[5]) / (xyzLoHi[5]-xyzLoHi[4])
				endif
				String dataStr=""
				if (strsearch(action,"-X",0,2)>=0)
					dataStr = "{0,1,0,0,"
					clipBox *= -1
				elseif (strsearch(action,"+X",0,2)>=0)
					dataStr = "{1,-1,0,0,"
				elseif (strsearch(action,"-Y",0,2)>=0)
					dataStr = "{2,0,1,0,"
					clipBox *= -1
				elseif (strsearch(action,"+Y",0,2)>=0)
					dataStr = "{3,0,-1,0,"
				elseif (strsearch(action,"-Z",0,2)>=0)
					dataStr = "{4,0,0,1,"
					clipBox *= -1
				elseif (strsearch(action,"+Z",0,2)>=0)
					dataStr = "{5,0,0,-1,"
				endif
				dataStr += num2str(clipBox)+"}"

				Execute "ModifyGizmo currentGroupObject=\""+groupName+"\""
				Execute "GetGizmo displayNameList"
				String names=StrVarOrDefault("S_DisplayNames","")
				KillStrings/Z S_DisplayNames
				if (WhichListItem(clipPlane,names)<0)				// not present, add it
					Execute "ModifyGizmo setDisplayList=-1, opName="+clipPlane+", operation=ClipPlane, data="+dataStr
				else													// just modify existing
					Execute "ModifyGizmo/Z opName="+clipPlane+", operation=ClipPlane, data="+dataStr
				endif
				Execute "ModifyGizmo currentGroupObject=\"::\""
			endif
			Execute "ModifyGizmo userString={"+groupName+",\""+keyVals+"\"}"
			Execute "ModifyGizmo compile"
			break

		default:
			err = 1													// unknown action
	endswitch
	return err
End

//  ============================ End of Add to orModify Gizmo ============================  //
//  ======================================================================================  //




//  ======================================================================================  //
//  =================== Start of Examining things Displayed in a Gizmo ===================  //

Function isGizmoObjectDisplayed(object,[gizmo])		// returns -1 if object no found on displayList, if found return place in displayList
	String object						// name of a gizmo object
	String gizmo
	String Nswitch = ""
	if (!ParamIsDefault(gizmo) && strlen(gizmo))
		if (WinType(gizmo)!=13)
			return -1
		endif
		Nswitch = "/N="+gizmo
	endif

	Execute "GetGizmo/Z"+Nswitch+" displayList"
	String displayList=StrVarOrDefault("S_DisplayList","")
	KillWaves/Z TW_DisplayList
	KillStrings/Z S_DisplayList

	String item, find="*, object="+object
	Variable i,i0, num=-1
	for (i=0;i<ItemsInList(displayList);i+=1)
		item = StringFromList(i,displayList)
		if (stringmatch(item,find))
			i0 = strsearch(item," setDisplayList=",0,2)
			if (i0>0)
				i0 += strlen(" setDisplayList=")
				num = str2num(item[i0,Inf])
			endif
			break
		endif
	endfor
	num = num>=0 ? num : -1
	return num
End



Function/T GetGizmoObjects(type,[gizmo,exclude])		// returns list of objects of type "type", from top gizmo or named gizmo
	String type						// type of object to look for, e.g. "scatter"
	String gizmo
	String exclude					// an optional list of objects to exclude

	String Nswitch = ""
	if (!ParamIsDefault(gizmo) && strlen(gizmo))
		if (WinType(gizmo)!=13)
			return ""
		endif
		Nswitch = "/N="+gizmo
	endif
	if (ParamIsDefault(exclude))
		exclude = ""
	endif

	Execute "GetGizmo/Z"+Nswitch+" objectList"
	String objectList=StrVarOrDefault("S_gizmoObjectList","")
	KillWaves/Z TW_gizmoObjectList
	KillStrings/Z S_gizmoObjectList
	String item, objects="", objName
	String typeSep=type+SelectString(stringmatch(type,"group"),"=",",")
	Variable i
	for (i=0;i<ItemsInList(objectList,";");i+=1)
		item = StringFromList(i,objectList,";")
		if (strsearch(item,"AppendToGizmo "+typeSep,0,2)==1)
			objName = TrimFrontBackWhiteSpace(StringByKey("name",item,"=",","))
			if (WhichListItem(objName,exclude)<0)
				objects += objName+";"
			endif
		endif
	endfor
	return objects
End


Function/T GizmoListScatterWaves([gizmo])		// get list of all scatter plots, except 'scatterMarker0, scatterMarker1, ...
	String gizmo									// optional name of gizmo
	String Nswitch = ""
	gizmo = SelectString(ParamIsDefault(gizmo),gizmo,"")
	if (strlen(gizmo))
		if (WinType(gizmo)!=13)
			return ""
		endif
		Nswitch = "/N="+gizmo
	endif

	Execute "ModifyGizmo"+Nswitch+" stopRotation"
	Execute "GetGizmo"+Nswitch+" displayList"		// list of all displayed objects
	Execute "GetGizmo"+Nswitch+" objectList"		// find name of wave with marker position
	Wave/T TW_DisplayList=TW_DisplayList, TW_gizmoObjectList=TW_gizmoObjectList

	String keyList=""
	String str, wname, list, scatterName
	Variable i,i0,i1
	for (i=0;i<numpnts(TW_gizmoObjectList);i+=1)
		str = TW_gizmoObjectList[i]

		i0 = strsearch(str,"AppendToGizmo ",0)
		if (i0<0 || i0>1)							// only process AppendToGizmo commands
			continue
		endif
		list = str[i0+strlen("AppendToGizmo "),Inf]

		scatterName = StringByKey("name",list,"=",",")
		if (GrepString(scatterName,"scatterMarker[0-9]+"))// skip scatterMarker0, scatterMarker1, ...
			continue
		elseif (stringmatch(scatterName,"scatterMarkerJZT"))	// skip scatterMarkerJZT too.
			continue
		endif

		wname = StringByKey("Scatter",list,"=",",")
		if (strlen(wname)<1)
			continue
		endif
		keyList = ReplaceStringByKey(wname,keyList,scatterName,"=")
	endfor
	KillWaves/Z TW_DisplayList, TW_gizmoObjectList
	KIllStrings/Z S_DisplayList, S_gizmoObjectList
	return keyList
End


Function/T GizmoListIsoSurfaceWaves([gizmo])		// get list of all iso surface waves
	String gizmo									// optional name of gizmo
	String Nswitch = ""
	gizmo = SelectString(ParamIsDefault(gizmo),gizmo,"")
	if (strlen(gizmo))
		if (WinType(gizmo)!=13)
			return ""
		endif
		Nswitch = "/N="+gizmo
	endif

	Execute "ModifyGizmo"+Nswitch+" stopRotation"
	Execute "GetGizmo"+Nswitch+" displayList"		// list of all displayed objects
	Execute "GetGizmo"+Nswitch+" objectList"		// find name of wave with marker position
	Wave/T TW_DisplayList=TW_DisplayList, TW_gizmoObjectList=TW_gizmoObjectList

	String keyList=""
	String str, wname, list, scatterName
	Variable i,i0,i1
	for (i=0;i<numpnts(TW_gizmoObjectList);i+=1)
		str = TW_gizmoObjectList[i]

		i0 = strsearch(str,"AppendToGizmo ",0)
		if (i0<0 || i0>1)							// only process AppendToGizmo commands
			continue
		endif
		list = str[i0+strlen("AppendToGizmo "),Inf]

		scatterName = StringByKey("name",list,"=",",")
		wname = StringByKey("isoSurface",list,"=",",")
		if (strlen(wname)<1)
			continue
		endif
		keyList = ReplaceStringByKey(wname,keyList,scatterName,"=")
	endfor
	KillWaves/Z TW_DisplayList, TW_gizmoObjectList
	KIllStrings/Z S_DisplayList, S_gizmoObjectList
	return keyList
End


Function/T GizmoListSurfaceWaves([gizmo])		// get list of all Surface plots on gizmo
	String gizmo									// optional name of gizmo
	String Nswitch = ""
	gizmo = SelectString(ParamIsDefault(gizmo),gizmo,"")
	if (strlen(gizmo))
		if (WinType(gizmo)!=13)
			return ""
		endif
		Nswitch = "/N="+gizmo
	endif

	Execute "ModifyGizmo"+Nswitch+" stopRotation"
	Execute "GetGizmo"+Nswitch+" displayList"		// list of all displayed objects
	Execute "GetGizmo"+Nswitch+" objectList"		// find name of wave with marker position
	Wave/T TW_DisplayList=TW_DisplayList, TW_gizmoObjectList=TW_gizmoObjectList

	String keyList=""
	String str, wname, list, surfaceName
	Variable i,i0,i1
	for (i=0;i<numpnts(TW_gizmoObjectList);i+=1)
		str = TW_gizmoObjectList[i]

		i0 = strsearch(str,"AppendToGizmo ",0)
		if (i0<0 || i0>1)							// only process AppendToGizmo commands
			continue
		endif
		list = str[i0+strlen("AppendToGizmo "),Inf]

		surfaceName = StringByKey("name",list,"=",",")

		wname = StringByKey("Surface",list,"=",",")
		if (strlen(wname)<1)
			continue
		endif
		keyList = ReplaceStringByKey(wname,keyList,surfaceName,"=")
	endfor
	KillWaves/Z TW_DisplayList, TW_gizmoObjectList
	KIllStrings/Z S_DisplayList, S_gizmoObjectList
	return keyList
End


Function/T GetGizmoBoxDisplayed(gizName)		// returns list with XYZ range of gizmo in USER units
	String gizName				// use empty string for top gizmo
	if (itemsInList(WinList("*",";","WIN:4096"))==0)
		return ""
	endif

	String Nstr=SelectString(strlen(gizName),"","/N="+gizName)

	Execute "GetGizmo"+Nstr+" userBoxLimits"
	NVAR GizmoBoxXmin=GizmoBoxXmin, GizmoBoxXmax=GizmoBoxXmax
	NVAR GizmoBoxYmin=GizmoBoxYmin, GizmoBoxYmax=GizmoBoxYmax
	NVAR GizmoBoxZmin=GizmoBoxZmin, GizmoBoxZmax=GizmoBoxZmax
	//	printf "GizmoBox = [%g, %g] [%g, %g] [%g, %g]\r",GizmoBoxXmin, GizmoBoxXmax,GizmoBoxYmin, GizmoBoxYmax,GizmoBoxZmin, GizmoBoxZmax
	Variable Xlo=GizmoBoxXmin,Xhi=GizmoBoxXmax,Ylo=GizmoBoxYmin,Yhi=GizmoBoxYmax,Zlo=GizmoBoxZmin,Zhi=GizmoBoxZmax
	KillVariables/Z GizmoBoxXmin, GizmoBoxXmax, GizmoBoxYmin, GizmoBoxYmax, GizmoBoxZmin, GizmoBoxZmax
	if (numtype(Xlo+Xhi+Ylo+Yhi+Zlo+Zhi))
		Execute "GetGizmo"+Nstr+" dataLimits"
		NVAR GizmoXmin=GizmoXmin, GizmoXmax=GizmoXmax
		NVAR GizmoYmin=GizmoYmin, GizmoYmax=GizmoYmax
		NVAR GizmoZmin=GizmoZmin, GizmoZmax=GizmoZmax
		Xlo=GizmoXmin;Xhi=GizmoXmax; Ylo=GizmoYmin;Yhi=GizmoYmax; Zlo=GizmoZmin;Zhi=GizmoZmax
		//	printf "GizmoRange = [%g, %g] [%g, %g] [%g, %g]\r",GizmoXmin, GizmoXmax,GizmoYmin, GizmoYmax,GizmoZmin, GizmoZmax
		KillVariables/Z GizmoXmin, GizmoXmax, GizmoYmin, GizmoYmax, GizmoZmin, GizmoZmax
	endif
	String out=""
	if (numtype(Xlo+Xhi+Ylo+Yhi+Zlo+Zhi)==0)
		sprintf out,"%g;%g;%g;%g;%g;%g",Xlo,Xhi,Ylo,Yhi,Zlo,Zhi
	endif
	return out
End


Function/WAVE GizmoObjectWave(objectName,objType,[gizmo])
	// get the wave associated with objectName
	String objectName					// name of object, e.g. "scatter0", "CubeCornersScatter", ...
	String objType						// object type, "Scatter", "Surface", "voxelgram", ...
	String gizmo						// optional name of gizmo
	gizmo = SelectString(ParamIsDefault(gizmo),gizmo,"")
	if (strlen(gizmo))
		if (WinType(gizmo)!=13)
			return $""
		endif
	endif

	// get all ojects
	Execute "GetGizmo"+SelectString(strlen(gizmo),"","/N="+gizmo)+" objectList"
	Wave/T TW_gizmoObjectList=TW_gizmoObjectList

	Wave ww = $""
	Variable i,i0,i1,N=DimSize(TW_gizmoObjectList,0)
	String list, name
	for (i=0;i<N;i+=1)
		list = TrimFrontBackWhiteSpace(TW_gizmoObjectList[i])
		if (strsearch(list,"AppendToGizmo ",0)!=0)
			continue
		endif
		list = ReplaceString("AppendToGizmo ",list,"")
		list = TrimFrontBackWhiteSpace(list)
		name = StringByKey("name",list,"=",",")
		if (StringMatch(name,objectName))
			name = TrimFrontBackWhiteSpace(StringByKey(objType,list,"=",","))
			Wave ww = $StringByKey(objType,list,"=",",")
			break
		endif
	endfor
	KillWaves/Z TW_gizmoObjectList
	return ww
End

//  ==================== End of Examining things Displayed in a Gizmo ====================  //
//  ======================================================================================  //



//  ======================================================================================  //
//  ====================== Start of QUESTIONABLE items for this file =====================  //

// Will add a parametric XYZ wave to the gizmo with its associated RGBA wave, this is a pretty generic routine
// returns name of object for parametricXYZ
Static Function/T AppendParametricXYZ2Gizmo(parametricXYZ,class)
	Wave parametricXYZ
	String class							// wave class

	if(exists("NewGizmo")!=4)			// Do nothing if the Gizmo XOP is not available.
		DoAlert 0, "Gizmo XOP must be installed"
		return ""
	endif

	if (!WaveExists(parametricXYZ))
		class = SelectString(strlen(class),"*",class)
		String wStr="",wList = WaveListClass(class,"*","DIMS:3,MAXLAYERS:3,MINLAYERS:3")
		if (ItemsInList(wList)<1)
			DoAlert 0,"Parametric XYZ wave not found"
			return ""
		elseif (ItemsInList(wList)==1)
			wStr = StringFromLIst(0,wList)
		else
			Prompt wStr, "Parametric XYZ wave",popup, wList
			DoPrompt "Parametric XYZ ("+class+")",wStr
			if (V_flag)
				return ""
			endif
		endif
		Wave parametricXYZ = $wStr
		printf "AppendParametricXYZ2Gizmo(%s)\r",NameOfWave(parametricXYZ)
	endif
	if (!WaveExists(parametricXYZ))
		DoAlert 0,"Parametric XYZ wave not found"
		return ""
	endif
	if (WaveDims(parametricXYZ)!=3 || DimSize(parametricXYZ,2)!=3)
		DoAlert 0,"Parametric XYZ wave "+NameOfWave(parametricXYZ)+" is not of dimension [N][M][3]"
		return ""
	endif

	Variable i,j
	String win=GizmosWithWave(parametricXYZ)
	String xyzName = GetWavesDataFolder(parametricXYZ,2)
	if (strlen(win)>0)
		DoWindow/F $win				// for wave already in gizmo, bring it to front
		Execute "GetGizmo/N="+win+"/Z objectList"
		Wave/T TW_gizmoObjectList=TW_gizmoObjectList
		for (i=0,j=-1; i<numpnts(TW_gizmoObjectList) && j<0; i+=1)
			j = strsearch(TW_gizmoObjectList[i],xyzName,0,2)
		endfor
		return StringByKey("name",TW_gizmoObjectList[i-1],"=",",")
	endif

	String wnote=note(parametricXYZ)
	Wave sourceWave = $StringByKey("source", wnote,"=")
	if (!WaveExists(sourceWave))
		return ""
	endif
	win = GizmosWithWave(sourceWave)
	if (strlen(win)<1)
		return ""
	endif
	DoWindow/F $win					// bring gizmo to front, and then append the surface
	DoUpdate

	// Find unsued name for surface, starting with "ParametricSurf0"
	Execute "GetGizmo/N="+win+"/Z objectNameList"
	SVAR S_ObjectNames=S_ObjectNames
	String rgbaName=StringByKey("RGBA",wnote,"="), surfaceName
	for (i=0,j=0; j>=0; i+=1)
		surfaceName = "ParametricSurf"+num2istr(i)
		j = WhichListItem(surfaceName, S_ObjectNames)
	endfor
	Execute "AppendToGizmo Surface="+xyzName+",name="+surfaceName
	Execute "ModifyGizmo ModifyObject="+surfaceName+" property={ srcMode,4}"
	Execute "ModifyGizmo modifyObject="+surfaceName+" property={calcNormals,1}"
	if (exists(rgbaName)==1)
		Execute "ModifyGizmo ModifyObject="+surfaceName+" property={ surfaceColorType,3}"
		Execute "ModifyGizmo ModifyObject="+surfaceName+" property={ surfaceColorWave,"+rgbaName+"}"
	endif
	Execute "ModifyGizmo setDisplayList=-1, object="+surfaceName
	return surfaceName
End



// Maybe this should not be here:
Static Function/T GizmosWithWave(scatter)
	Wave scatter

	if (!WaveExists(scatter))
		return ""
	endif
	String list = WinList("*", ";", "WIN:4096" )	// list of possible gizmos
	Variable N=ItemsInlist(list)
	String str, wName=GetWavesDataFolder(scatter,2)
	Variable i
	for (i=N-1;i>=0;i-=1)
		str = WinRecreation(StringFromList(i,list), 0)
		if (strsearch(str,wName,0)<0)
			list = RemoveListItem(i,list)
		endif
	endfor
	N=ItemsInlist(list)

	if (N<1)
		return ""
	elseif (N==1)
		return StringFromList(0,list)
	endif

	String wGizmo
	Prompt wGizmo,"Which Gizmo?",popup,list
	DoPrompt "Gizmo?",wGizmo
	if (V_flag)
		return ""
	endif
	return wGizmo
End

//  ======================= End of QUESTIONABLE items for this file ======================  //
//  ======================================================================================  //




Static Function InitGizmoUtilityGeneral()
	Execute/Q/Z "GizmoMenu AppendItem={JZT0,\"-\", \"\"}"
	Execute/Q/Z "GizmoMenu AppendItem={JZT1,\"Put Cube Corners on Gizmo\", \"ActivateCornerCubeOnGizmo(forceCalc=1)\"}"
	Execute/Q/Z "GizmoMenu AppendItem={JZT2,\"De-Activate Cube Corners on Gizmo\", \"DeActivateCornerCubeOnGizmo()\"}"
	Execute/Q/Z "GizmoMenu AppendItem={JZT3,\"-\", \"\"}"
	if (strlen(WinList("microGeometryN.ipf", ";","WIN:128")))
		Execute/Q/Z "GizmoMenu AppendItem={JZTr0,\"Gizmo X-H plane\", \"ModifyGizmo stopRotation ; ModifyGizmo SETQUATERNION={0,0,0,1}\"}"	// X-H
		Execute/Q/Z "GizmoMenu AppendItem={JZTr1,\"Gizmo X-F plane\", \"ModifyGizmo stopRotation ; ModifyGizmo SETQUATERNION={0.707106781186547,0,0,0.707106781186547}\"}"	// X-F
		Execute/Q/Z "GizmoMenu AppendItem={JZTr2,\"Gizmo H-F plane\", \"ModifyGizmo stopRotation ; ModifyGizmo SETQUATERNION={0.5,0.5,0.5,0.5}\"}"	// H-F
		Execute/Q/Z "GizmoMenu AppendItem={JZTr3,\"Gizmo Beamline\", \"ModifyGizmo stopRotation ; ModifyGizmo SETQUATERNION={-0.270598,0.653282,-0.270598,0.653281}\"}"	// ModifyGizmo euler={0,-90,45}
	else
		Execute/Q/Z "GizmoMenu AppendItem={JZTr0,\"Gizmo X-Y plane [along beam]\", \"ModifyGizmo stopRotation ; ModifyGizmo SETQUATERNION={0.0,0.0,0.0,1.0}\"}"
		Execute/Q/Z "GizmoMenu AppendItem={JZTr1,\"Gizmo Y-Z plane [side view]\", \"ModifyGizmo stopRotation ; ModifyGizmo SETQUATERNION={0.5,0.5,0.5,0.5}\"}"
		Execute/Q/Z "GizmoMenu AppendItem={JZTr2,\"Gizmo X-Z plane\", \"ModifyGizmo stopRotation ; ModifyGizmo SETQUATERNION={0.707107,0.0,0.0,0.707107}\"}"
		Execute/Q/Z "GizmoMenu AppendItem={JZTr3,\"Gizmo X-Z plane [top view]\", \"ModifyGizmo stopRotation ; ModifyGizmo SETQUATERNION={-0.707107,0.0,0.0,0.707107}\"}"
		Execute/Q/Z "GizmoMenu AppendItem={JZTr4,\"Gizmo 111 vertical [0-11 to right]\", \"ModifyGizmo stopRotation ; ModifyGizmo SETQUATERNION={0.060003,-0.540625,-0.455768,0.704556}\"}"
		//	ModifyGizmo euler = {45, acos(1/sqrt(3))*180/PI, 90}
	endif
	Execute/Q/Z "GizmoMenu AppendItem={JZT4,\"-\", \"\"}"
End
