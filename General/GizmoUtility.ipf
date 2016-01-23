#pragma rtGlobals=2		// Use modern global access method.
#pragma ModuleName=GizmoUtil
#pragma IgorVersion = 6.21
#pragma version = 2.13
#include "ColorNames"

Static Constant GIZMO_MARKER_END_SIZE = 0.07		// puts boxes on ends of 3D marker (you can OverRide this in the Main procedure)
Static Constant GIZMO_MARKER_END_TYPE = 1			// 0=box on ends of lines, 1=sphere on ends of lines (you can OverRide this in the Main procedure)
Static Constant MAX_N_OBJECTS=500					// maximum number of objects that can be created using these routines (such as scatter, groups, titles, ...)
#if (IgorVersion()<7)
Static Constant GIZMO_SCALE_BAR_LEFT_EDGE = -1.9	// left edge of scale bar on a Gizmo
#else
Static Constant GIZMO_SCALE_BAR_LEFT_EDGE = -0.95	// left edge of scale bar on a Gizmo
#endif


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
#if (IgorVersion()<7)
		String Nswitch=SelectString(strlen(GizmoName),"","/N="+GizmoName)
		Execute "ModifyGizmo"+Nswitch+" setDisplayList=-1, object="+objectName
		Execute "ModifyGizmo"+Nswitch+" compile"
#else
		ModifyGizmo/N=$GizmoName setDisplayList=-1, object=$objectName
		ModifyGizmo/N=$GizmoName compile
#endif
	endif
	return 0
End


Function DeActivateCornerCubeOnGizmo([GizmoName])			// returns true if the gizmo corner cubes are un-displayed
	String GizmoName			// optional name of gizmo, defaults to top gizmo
	GizmoName = SelectString(ParamIsDefault(GizmoName),GizmoName,"")

	if (isCornerCubeDisplayedOnGizmo(GizmoName=GizmoName))
		String objectName=getCornerCubeObjectNameOnGizmo(GizmoName=GizmoName)	// if "", then add the corner cubes
#if (IgorVersion()<7)
		String Nswitch=SelectString(strlen(GizmoName),"","/N="+GizmoName)
		Execute "RemoveFromGizmo"+Nswitch+" displayItem="+objectName
		Execute "ModifyGizmo"+Nswitch+" compile"
#else
		RemoveFromGizmo/N=$GizmoName displayItem=$objectName
		ModifyGizmo/N=$GizmoName compile
#endif
	endif
	return 0
End


Static Function/WAVE FindMakeCubeCornerWaves([GizmoName,forceCalc])		// finds (or makes) the corner cube wave for a gizmo
	String GizmoName			// optional name of gizmo, defaults to top gizmo
	Variable forceCalc		// force recalculation of corners
	GizmoName = SelectString(ParamIsDefault(GizmoName),GizmoName,"")
	forceCalc = ParamIsDefault(forceCalc) || numtype(forceCalc) ? 0 : !(!forceCalc)

	String cornerList=WaveListClass("GizmoCorners","*","DIMS:2,MINROWS:2,MAXROWS:2,MINCOLS:3,MAXCOLS:3")
	cornerList += WaveListClass("GizmoCorners","*","DIMS:2,MINROWS:8,MAXROWS:8,MINCOLS:3,MAXCOLS:3")
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
	elseif(! (WaveDims(wCorners)==2 && (DimSize(wCorners,0)==8 || DimSize(wCorners,0)==2) && DimSize(wCorners,1)==3) )
		return ""
	elseif (!ParamIsDefault(GizmoName) && strlen(GizmoName))
		if (WinType(GizmoName)!=GIZMO_WIN_TYPE)
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
#if (IgorVersion()<7)
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
#else
		ModifyGizmo startRecMacro
		AppendToGizmo/N=$GizmoName Scatter=$GetWavesDataFolder(wCorners,2), name=$objectName
		ModifyGizmo/N=$GizmoName ModifyObject=$objectName objectType=scatter property={ scatterColorType,0}
		ModifyGizmo/N=$GizmoName ModifyObject=$objectName objectType=scatter property={ markerType,0}
		ModifyGizmo/N=$GizmoName ModifyObject=$objectName objectType=scatter property={ sizeType,0}
		ModifyGizmo/N=$GizmoName ModifyObject=$objectName objectType=scatter property={ rotationType,0}
		ModifyGizmo/N=$GizmoName ModifyObject=$objectName objectType=scatter property={ Shape,1}
		ModifyGizmo/N=$GizmoName ModifyObject=$objectName objectType=scatter property={ size,1}
		ModifyGizmo/N=$GizmoName ModifyObject=$objectName objectType=scatter property={ color,0,0,0,1}
//		ModifyGizmo/N=$GizmoName userString={CubeCorners, objectName}	// save name of cube corner object
		SetWindow $GizmoName userdata(CubeCorners)="AtomViewCubeCorners"
		ModifyGizmo endRecMacro
#endif
	endif
	return objectName
End


Static Function/T getCornerCubeObjectNameOnGizmo([GizmoName])		// returns name of corner cube object on gizmo, it may not yet be displayed
	String GizmoName			// optional name of gizmo, defaults to top gizmo

	GizmoName = SelectString(ParamIsDefault(GizmoName),GizmoName,"")
	String Nswitch=""
	if (!ParamIsDefault(GizmoName) && strlen(GizmoName))
		if (WinType(GizmoName)!=GIZMO_WIN_TYPE)
			return ""
		endif
		Nswitch = "/N="+GizmoName
	endif

#if (IgorVersion()<7)
	Execute "GetGizmo/Z"+Nswitch+" userString=CubeCorners"
	String objectName=StrVarOrDefault("S_GizmoUserString","")
	KillStrings/Z S_GizmoUserString
#else
//	GetGizmo/Z/N=$GizmoName userString=CubeCorners
	String objectName = GetUserData(GizmoName,"","CubeCorners")
#endif

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
	Make/N=(2,3)/O $name/WAVE=corners =NaN
	SetScale d 0,0,units, corners
	corners[0][0] = Xlo;	corners[0][1] = Ylo; 	corners[0][2] = Zlo
	corners[1][0] = Xhi;	corners[1][1] = Yhi; 	corners[1][2] = Zhi
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
	pos = SelectString(ParamIsDefault(pos),pos,"LT")
	pos = SelectString(strlen(pos),"LT",pos)	// defaults to top left
	String font="Geneva"
	Variable fontSize=18

	if (WhichListItem(pos,"LT;LB;")<0)
		return ""
	elseif (strlen(title1)<1)
		return ""
	endif

	// remove empty strings, slide down
	Make/N=4/T/FREE titleWave={title1,title2,title3,title4}
	Variable i,j
	for (j=0;j<3;j+=1)									// remove any empty strings
		for (i=0;i<(4-j);i+=1)
			if (strlen(titleWave[0]))
				break
			endif
			titleWave[j,2] = titleWave[p+1]		// slide down 1
			titleWave[3] = ""
		endfor
	endfor
	Make/N=4/W/U/FREE useStr = ( strlen(titleWave[p]) > 0 )
	if (sum(useStr)<1)									// nothing to do
		return ""
	endif

#if (IgorVersion()>=7)
	groupName = "textTitle"
	String title = "\Z18"+titleWave[0]
	title += SelectString(useStr[1],"", "\r\Zr080"+titleWave[1])
	title += SelectString(useStr[2],"", "\r\Zr070"+titleWave[2])
	title += SelectString(useStr[3],"", "\r"+titleWave[3])
	TextBox/C/N=$groupName/F=0/B=1/A=$pos/X=2/Y=2 title

#else
	if (CheckName(groupName,5))
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

	if (stringmatch(pos,"LT"))
		Execute "ModifyGizmo setDisplayList=0, opName=translateTitle, operation=translate, data={-1.9,1.9,0}"
	elseif (stringmatch(pos,"LB"))
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
#endif

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
	font = SelectString(ParamIsDefault(font),font,"")
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
	Variable dxGizmoLine = BarLength/maxLength * scaleFactor

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

	// ************************* Group Object Start *******************
#if (IgorVersion()<7)
	font = SelectString(strlen(font),"Geneva",font)
	Execute "AppendToGizmo group,name="+groupName
	Execute "ModifyGizmo currentGroupObject=\""+groupName+"\""
	String cmd
	sprintf cmd "AppendToGizmo string=\"%s\",strFont=\"%s\",name=stringScaleBar",unitStr,font
	Execute cmd
	Execute "ModifyGizmo modifyObject=stringScaleBar property={Clipped,0}"
	dxGizmoLine *= 2
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
#else
	AppendToGizmo group,name=$groupName
	ModifyGizmo currentGroupObject=groupName
	AppendToGizmo line={0,0, 0,dxGizmoLine, 0,0}, name=lineScaleBar
	AppendToGizmo attribute lineWidth=5, name=widthScaleBar
	ModifyGizmo setDisplayList=0, opName=translate0, operation=translate, data={GIZMO_SCALE_BAR_LEFT_EDGE,GIZMO_SCALE_BAR_LEFT_EDGE,0}
	ModifyGizmo setDisplayList=1, attribute=widthScaleBar
	ModifyGizmo setDisplayList=2, object=lineScaleBar
	ModifyGizmo currentGroupObject="::"

	String str=ReplaceStringByKey("object","",groupName,"=")
	str = ReplaceNumberByKey("scaleFactor",str,scaleFactor,"=")
	str = ReplaceNumberByKey("BarLength",str,BarLength,"=")
	str = ReplaceStringByKey("units",str,units,"=")
	if (strlen(font))
		str = ReplaceStringByKey("font",str,font,"=")
	endif
	SetWindow kwTopWin,userdata(ScaleBar)=str

	Variable fontSize=18
	unitStr = SelectString(strlen(units),""," ")+units
	if (strlen(font))
		sprintf str "\\F'%s'\\Z%d%g%s", font,fontSize,BarLength,unitStr
	else
		sprintf str "\\Z%d%g%s", fontSize,BarLength,unitStr
	endif
	TextBox/C/N=textScaleBar/F=0/B=1/A=LB/X=3/Y=3 str
#endif
	// ************************* Group Object End *******************
	return groupName
End
//Function/T AddScaleBarGroup(groupName,maxLength,units,[scaleFactor,font])
//	String groupName			// probably "ScaleBarGroup0", If this is empty, then a unique name will be assigned
//	Variable maxLength		// maximum length in Gizmo, scale bar will be less than this
//	String units				// units for maxLength and used in scale bar label
//	Variable scaleFactor
//	String font
//	scaleFactor = ParamIsDefault(scaleFactor) ? 1 : scaleFactor
//	scaleFactor = numtype(scaleFactor) || scaleFactor<=0 ? 1 : scaleFactor
//	font = SelectString(ParamIsDefault(font),font,"Geneva")
//	font = SelectString(strlen(font),"Geneva",font)
//	if (maxLength<=0 || numtype(maxLength))
//		return ""
//	endif
//
//	// for scale bar use multipliers of 1, 2, or 5 ONLY
//	Variable BarLength = 10^floor(log(maxLength))
//	if (5*BarLength < maxLength)
//		BarLength = 5*BarLength
//	elseif (2*BarLength < maxLength)
//		BarLength = 2*BarLength
//	endif
//	String unitStr = num2str(BarLength)+" "+units
//	Variable dxGizmoLine = BarLength/maxLength * scaleFactor * 2
//
//	if (CheckName(groupName,5))						// invalid groupName passed, create one
//		if (strlen(groupName)<1)						// create a unique group name
//			NewDataFolder/O root:Packages			// ensure Packages exists
//			NewDataFolder/O root:Packages:JZT_GizmoUtility	// ensure geometry exists
//			if (exists("root:Packages:JZT_GizmoUtility:ScaleBarNumber")!=2)
//				Variable/G root:Packages:JZT_GizmoUtility:ScaleBarNumber=-1
//			endif
//			NVAR ScaleBarNumber = root:Packages:JZT_GizmoUtility:ScaleBarNumber
//			ScaleBarNumber = numtype(ScaleBarNumber) ? -1 : limit(round(ScaleBarNumber),-1,Inf)
//			ScaleBarNumber += 1
//			groupName = "ScaleBarGroup"+num2istr(ScaleBarNumber)
//		endif
//		groupName = CleanupName(groupName,0)
//	endif
//	if (CheckName(groupName,5))						// invalid groupName passed, give up
//		return ""
//	endif
//
//	// ************************* Group Object Start *******************
//#if (IgorVersion()<7)
//	Execute "AppendToGizmo group,name="+groupName
//	Execute "ModifyGizmo currentGroupObject=\""+groupName+"\""
//	String cmd
//	sprintf cmd "AppendToGizmo string=\"%s\",strFont=\"%s\",name=stringScaleBar",unitStr,font
//	Execute cmd
//	Execute "ModifyGizmo modifyObject=stringScaleBar property={Clipped,0}"
//	sprintf cmd "AppendToGizmo line={%g,-1.9,0,%g,-1.9,0}, name=line0",GIZMO_SCALE_BAR_LEFT_EDGE,dxGizmoLine+GIZMO_SCALE_BAR_LEFT_EDGE
//	Execute cmd
//	Execute "AppendToGizmo attribute lineWidth=5, name=lineWidth0"
//	Execute "ModifyGizmo setDisplayList=0, opName=pushMatrix0, operation=pushMatrix"
//	Execute "ModifyGizmo setDisplayList=1, attribute=lineWidth0"
//	Execute "ModifyGizmo setDisplayList=2, object=line0"
//	Execute "ModifyGizmo setDisplayList=3, opName=translateScaleBarText, operation=translate, data={"+num2str(GIZMO_SCALE_BAR_LEFT_EDGE)+",-1.85,0}"
//	Execute "ModifyGizmo setDisplayList=4, opName=scaleScaleBarText, operation=scale, data={0.1,0.1,0.1}"
//	Execute "ModifyGizmo setDisplayList=5, opName=rotateScaleBarText, operation=rotate, data={180,1,0,0}"
//	Execute "ModifyGizmo setDisplayList=6, object=stringScaleBar"
//	Execute "ModifyGizmo setDisplayList=7, opName=popMatrix0, operation=popMatrix"
//	Execute "ModifyGizmo currentGroupObject=\"::\""
//#else
//	AppendToGizmo group,name=$groupName
//	ModifyGizmo currentGroupObject=groupName
//	AppendToGizmo string=unitStr,strFont=font,name=stringScaleBar
//// xxxxxxxxxxxxxxxxxxxxxxxxxxx
////	ModifyGizmo modifyObject=stringScaleBar objectType=string, property={Clipped,0}
//	AppendToGizmo line={GIZMO_SCALE_BAR_LEFT_EDGE,-1.9,0,dxGizmoLine+GIZMO_SCALE_BAR_LEFT_EDGE,-1.9,0}, name=line0
//	AppendToGizmo attribute lineWidth=5, name=lineWidth0
//	ModifyGizmo setDisplayList=0, opName=pushMatrix0, operation=pushMatrix
//	ModifyGizmo setDisplayList=1, attribute=lineWidth0
//	ModifyGizmo setDisplayList=2, object=line0
//	ModifyGizmo setDisplayList=3, opName=translateScaleBarText, operation=translate, data={GIZMO_SCALE_BAR_LEFT_EDGE,-1.85,0}
////	ModifyGizmo setDisplayList=4, opName=scaleScaleBarText, operation=scale, data={0.1,0.1,0.1}
//	ModifyGizmo setDisplayList=4, opName=rotateScaleBarText, operation=rotate, data={180,1,0,0}
//	ModifyGizmo setDisplayList=5, object=stringScaleBar
//	ModifyGizmo setDisplayList=6, opName=popMatrix0, operation=popMatrix
//	ModifyGizmo currentGroupObject="::"
//#endif
//	// ************************* Group Object End *******************
//	return groupName
//End
//


#if (IgorVersion()<7)
Static Function GzimoReSetScaleBarHookProc(s)
	STRUCT WMGizmoHookStruct &s
	if (!StringMatch(s.eventName,"scale"))
		return 0
	endif

	Variable scaleFactor=1
	String units = "nm"
	String win=s.winName
	String scaleBarName="", str, font="Geneva"
	Execute "GetGizmo/N="+win+"/Z objectNameList"
	String listObjectNames=StrVarOrDefault("S_ObjectNames","")
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
	Variable Vflag = (NumVarOrDefault("V_flag",0))
	KillVariables/Z V_Flag
	if (Vflag)
		Execute "GetGizmo/N="+win+"/Z displayList"
		KillStrings/Z S_DisplayList
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
	KillStrings/Z S_gizmoObjectList
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
#else
Static Function GzimoReSetScaleBarHookProc(s)
	STRUCT WMGizmoHookStruct &s
	if (!StringMatch(s.eventName,"scale"))
		return 0
	endif
	ResetScaleBarLength(win=s.winName)
End

Function ResetScaleBarLength([units,scaleFactor,font,fontSize,win])
	String units
	Variable scaleFactor
	String font
	Variable fontSize
	String win
	win = SelectString(ParamIsDefault(win),win,WinName(0,GIZMO_WIN_BIT))
	if (strlen(win)<1)
		return NaN
	endif
	String wnote=GetUserData(win,"","ScaleBar")
	String groupName = StringByKey("object",wnote,"=")
	if (strlen(groupName)<1)
		return NaN
	endif
	units = SelectString(ParamIsDefault(units), units, StringByKey("units",wnote,"="))
	scaleFactor = ParamIsDefault(scaleFactor) || numtype(scaleFactor) ? NumberByKey("scaleFactor",wnote,"=") : scaleFactor
	scaleFactor = numtype(scaleFactor) ? 1 : scaleFactor
	font = SelectString(ParamIsDefault(font), font, StringByKey("font",wnote,"="))
	fontSize = ParamIsDefault(fontSize) || numtype(fontSize) ? NumberByKey("fontSize",wnote,"=") : fontSize
	fontSize = numtype(fontSize) ? 18 : limit(round(fontSize),6,64)

	GetGizmo/N=$win/Z dataLimits
	Variable dX=GizmoXmax-GizmoXmin, dY=GizmoYmax-GizmoYmin, dZ=GizmoZmax-GizmoZmin
	if (dX<=0 || dY<=0 || dZ<=0 || numtype(dX+dY+dZ))
		return NaN
	endif
	Variable maxLen = max(dZ, max(dY, dX))
	maxLen = maxLen <=0 || numtype(maxLen) ? NaN : maxLen
//	print "   Gizmo XYZ range",GizmoXmax,GizmoXmin, GizmoYmax,GizmoYmin, GizmoZmax,GizmoZmin
//	print "   dX,dY,dZ,",dX,dY,dZ, "max=",maxLen

	// for scale bar use multipliers of 1, 2, or 5 ONLY
	Variable BarLength = 10^floor(log(maxLen))
	if (5*BarLength < maxLen)
		BarLength = 5*BarLength
	elseif (2*BarLength < maxLen)
		BarLength = 2*BarLength
	endif
	Variable dxGizmoLine = BarLength/maxLen * scaleFactor
//	printf "   BarLength = %g '%s',   dxGizmoLine = %g,   maxLen = %g\r",BarLength,units,dxGizmoLine,maxLen
	String object = win+":"+groupName
//printf "   object =\"%s\"\r",object

	ModifyGizmo/N=$win currentGroupObject=object
	Modifygizmo/N=$win modifyobject=lineScaleBar,objectType=line,property={vertex,0,0,0,dxGizmoLine,0,0}
	ModifyGizmo/N=$win currentGroupObject="::"

	units = SelectString(strlen(units),""," ")+units
	String str
	if (strlen(font))
		sprintf str "\\F'%s'\\Z%d%g%s", font,fontSize,BarLength,units
	else
		sprintf str "\\Z%d%g%s", fontSize,BarLength,units
	endif
	TextBox/C/N=textScaleBar/F=0/B=1/A=LB/X=3/Y=3 str
	return BarLength
End

#endif
//Static Function GzimoReSetScaleBarHookProc(s)
//	STRUCT WMGizmoHookStruct &s
//	if (!StringMatch(s.eventName,"scale"))
//		return 0
//	endif
//
//	Variable scaleFactor=1
//	String units = "nm"
//	String win=s.winName
//	String scaleBarName="", str, font="Geneva"
//#if (IgorVersion()<7)
//	Execute "GetGizmo/N="+win+"/Z objectNameList"
//	String listObjectNames=StrVarOrDefault("S_ObjectNames","")
//	KillStrings/Z S_ObjectNames
//#else
//	GetGizmo/N=$win/Z objectNameList
//	String listObjectNames=S_ObjectNames
//#endif
//	Variable i
//	for (i=0;i<ItemsInList(listObjectNames);i+=1)
//		if (StringMatch(StringFromList(i,listObjectNames),"ScaleBarGroup*"))
//			scaleBarName = StringFromList(i,listObjectNames)
//			break
//		endif
//	endfor
//	if (strlen(scaleBarName)<1)						// no scale bar, do nothing
//		return 0
//	endif
//
//#if (IgorVersion()<7)
//	Execute "GetGizmo/N="+win+"/Z displayItemExists=scaleBase"
//	Variable Vflag = (NumVarOrDefault("V_flag",0))
//	KillVariables/Z V_Flag
//#else
//	GetGizmo/N=$win/Z displayItemExists=scaleBase
//	Variable Vflag = V_flag
//#endif
//	if (Vflag)
//#if (IgorVersion()<7)
//		Execute "GetGizmo/N="+win+"/Z displayList"
//		KillStrings/Z S_DisplayList
//#else
//		GetGizmo/N=$win/Z displayList
//#endif
//		Wave/T TW_DisplayList=TW_DisplayList
//
//		for (i=0;i<numpnts(TW_DisplayList);i+=1)
//			if (strsearch(TW_DisplayList[i],"opName=scaleBase, operation=scale,",0,2)>0)
//				str = TW_DisplayList[i]
//  					i = strsearch(str,"data={",0,2)
// 					scaleFactor = str2num(str[i+6,Inf])
// 					scaleFactor = numtype(scaleFactor) ? 1 : scaleFactor
//				break
//			endif
//
//		endfor
//		KillWaves/Z TW_DisplayList
//	endif
//
//#if (IgorVersion()<7)
//	Execute "GetGizmo/N="+win+"/Z dataLimits"
//	Variable dx = NumVarOrDefault("GizmoXmax",1)-NumVarOrDefault("GizmoXmin",1)
//	Variable dy = NumVarOrDefault("GizmoYmax",1)-NumVarOrDefault("GizmoYmin",1)
//	Variable dz = NumVarOrDefault("GizmoZmax",1)-NumVarOrDefault("GizmoZmin",1)
//	KillVariables/Z GizmoXmin,GizmoXmax,GizmoYmin,GizmoYmax,GizmoZmin,GizmoZmax
//#else
//	GetGizmo/N=$win/Z dataLimits
//	Variable dx=GizmoXmax-GizmoXmin, dy=GizmoYmax-GizmoYmin, dz=GizmoZmax-GizmoZmin
//#endif
//
//	Variable maxLength = max(max(dx,dy),dz)
//	Variable BarLength = 10^floor(log(maxLength))
//	if (5*BarLength < maxLength)					// for scale bar use multipliers of 1, 2, or 5 ONLY
//		BarLength = 5*BarLength
//	elseif (2*BarLength < maxLength)
//		BarLength = 2*BarLength
//	endif
//	Variable dxGizmoLine = BarLength/maxLength * scaleFactor * 2
//	//	printf "maxLength = %g,  BarLength = %g,   dxGizmoLine = %g,  dxGizmoLine-1.75 = %g\r",maxLength,BarLength,dxGizmoLine, dxGizmoLine-1.75
//
//#if (IgorVersion()<7)
//	Execute "GetGizmo/N="+win+"/Z objectList"
//	KillStrings/Z S_gizmoObjectList
//#else
//	GetGizmo/N=$win/Z objectList
//#endif
//	Wave/T TW_gizmoObjectList=TW_gizmoObjectList
//
//	Variable i0,i1
//	for (i0=0,i=0; i<numpnts(TW_gizmoObjectList); i+=1)
//		if (strsearch(TW_gizmoObjectList[i], "AppendToGizmo group,name="+scaleBarName,0,2)>0)
//			i0 = i
//			break
//		endif
//	endfor
//	for (i1=0,i=i0; i<numpnts(TW_gizmoObjectList); i+=1)
//		if (strsearch(TW_gizmoObjectList[i], "ModifyGizmo currentGroupObject=\"::\"",0,2)>0)
//			i1 = i
//			break
//		endif
//	endfor
//
//	for (i=i0;i<i1; i+=1)
//		if (strsearch(TW_gizmoObjectList[i], "name=stringScaleBar",0,2)>0)
//			str = TW_gizmoObjectList[i]
//			i = strsearch(str, "AppendToGizmo string=\"",0,2)
//			if (i<0)
//				break
//			endif
//			str = str[i+22,Inf]
//			i = strsearch(str,"\",",0)
//			if (i<0)
//				break
//			endif
//			str = str[0,i-1]
//			i = strsearch(str," ",Inf,1)
//			if (i<0)
//				break
//			endif
//			str = str[i+1,Inf]
//			units = str
//			break
//		endif
//	endfor
//	//	printf "units = '%s'\r",units
//	KillWaves/Z TW_gizmoObjectList
//	String unitStr = num2str(BarLength)+" "+units
//
//#if (IgorVersion()<7)
//	Execute "ModifyGizmo/N="+win+"/Z startRecMacro"
//	Execute "ModifyGizmo/N="+win+"/Z currentGroupObject=\""+scaleBarName+"\""
//	Execute "ModifyGizmo/N="+win+"/Z modifyObject=stringScaleBar, property={string,\""+unitStr+"\"}"
//	sprintf str "ModifyGizmo/N=%s/Z modifyObject=line0, property={vertex,%g,-1.9,0,%g,-1.9,0}", win,GIZMO_SCALE_BAR_LEFT_EDGE,dxGizmoLine+GIZMO_SCALE_BAR_LEFT_EDGE
//	Execute str
//	Execute "ModifyGizmo/N="+win+"/Z currentGroupObject=\"::\""
//	Execute "ModifyGizmo/N="+win+"/Z endRecMacro"
//#else
//	ModifyGizmo/N=$win/Z startRecMacro
//	ModifyGizmo/N=$win/Z currentGroupObject=scaleBarName
//	ModifyGizmo/N=$win/Z modifyObject=stringScaleBar, objectType=string, property={string,unitStr}
//	ModifyGizmo/N=$win/Z modifyObject=line0, objectType=line property={vertex,GIZMO_SCALE_BAR_LEFT_EDGE,-1.9,0,dxGizmoLine+GIZMO_SCALE_BAR_LEFT_EDGE,-1.9,0}
//	ModifyGizmo/N=$win/Z currentGroupObject="::"
//	ModifyGizmo/N=$win/Z endRecMacro
//#endif
//	return 0	 
//End



// Add a marker to a gizmo. Returns name of item to include in the display list.
// Note this just creates a marker group, the marker group needs to be part of a scatter plot to display it
Function/T AddGizmoMarkerGroup(groupName,[rgba,alpha,scale])
	String groupName			// probably "groupTitle", If this is empty, then a unique name will be assigned
	String rgba					// red, green, blue, or "" is black, or you can give your own rgba as "1,1,0,0.5"
	Variable alpha
	Variable scale				// optional scale factor, default is 1
	rgba = SelectString(ParamIsDefault(rgba),rgba,"")
	alpha = ParamIsDefault(alpha) ? 1 : alpha
	alpha = numtype(alpha) ? 1 : limit(alpha,0,1)
	scale = ParamIsDefault(scale) || numtype(scale) ? 1 : scale

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
#if (IgorVersion()<7)
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
	if (!(scale==1))
		sprintf str, "ModifyGizmo setDisplayList=-1, opName=scaleAll, operation=scale, data={%g,%g,%g}",scale,scale,scale	;	Execute str
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
#else
	AppendToGizmo group,name=groupName
	ModifyGizmo currentGroupObject=groupName
	String str
	Variable v = GIZMO_MARKER_END_SIZE>0 ? 1-2*GIZMO_MARKER_END_SIZE : 1
	AppendToGizmo line={-v,0,0,v,0,0}, name=lineX
	AppendToGizmo line={0,-v,0,0,v,0}, name=lineY
	AppendToGizmo line={0,0,-v,0,0,v}, name=lineZ
	if (GIZMO_MARKER_END_SIZE>0)
		Variable endType=limit(GIZMO_MARKER_END_TYPE,0,1)
		if (endType==0)							// boxes
			AppendToGizmo box={GIZMO_MARKER_END_SIZE,GIZMO_MARKER_END_SIZE,GIZMO_MARKER_END_SIZE},name=endMarkers
			ModifyGizmo modifyObject=endMarkers objectType=box property={calcNormals,1}
		elseif (endType==1)						// spheres
			AppendToGizmo sphere={GIZMO_MARKER_END_SIZE,25,25},name=endMarkers
		endif
	endif
	AppendToGizmo attribute lineWidth=2, name=lineWidthCross
	if (strlen(rgba))
		Wave wrgba = str2vec(str)
		AppendToGizmo attribute color={wrgba[0],wrgba[1],wrgba[2],wrgba[3]},name=colorLine
	endif
	ModifyGizmo setDisplayList=0, opName=pushAttribute0, operation=pushAttribute, data=5
	ModifyGizmo setDisplayList=1, attribute=lineWidthCross
	if (strlen(rgba))
		ModifyGizmo setDisplayList=2, attribute=colorLine
	endif
	if (!(scale==1))
		ModifyGizmo setDisplayList=-1, opName=scaleAll, operation=scale, data={scale,scale,scale}
	endif
	ModifyGizmo setDisplayList=-1, object=lineX
	ModifyGizmo setDisplayList=-1, object=lineY
	ModifyGizmo setDisplayList=-1, object=lineZ
	if (GIZMO_MARKER_END_SIZE>0)
		v = 1 - GIZMO_MARKER_END_SIZE
		ModifyGizmo setDisplayList=-1, opName=translateXminus, operation=translate, data={-v,0,0}
		ModifyGizmo setDisplayList=-1, object=endMarkers
		ModifyGizmo setDisplayList=-1, opName=translateXplus, operation=translate, data={2*v,0,0}
		ModifyGizmo setDisplayList=-1, object=endMarkers
		ModifyGizmo setDisplayList=-1, opName=translateYminus, operation=translate, data={-v,-v,0}
		ModifyGizmo setDisplayList=-1, object=endMarkers
		ModifyGizmo setDisplayList=-1, opName=translateYplus, operation=translate, data={0,2*v,0}
		ModifyGizmo setDisplayList=-1, object=endMarkers
		ModifyGizmo setDisplayList=-1, opName=translateZminus, operation=translate, data={0,-v,-v}
		ModifyGizmo setDisplayList=-1, object=endMarkers
		ModifyGizmo setDisplayList=-1, opName=translateZplus, operation=translate, data={0,0,2*v}
		ModifyGizmo setDisplayList=-1, object=endMarkers
	endif
	ModifyGizmo setDisplayList=-1, opName=popAttribute0, operation=popAttribute
	ModifyGizmo currentGroupObject="::"
#endif
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
#if (IgorVersion()<7)
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
#else
	AppendToGizmo group,name=groupName
	ModifyGizmo currentGroupObject=groupName
	AppendToGizmo freeAxesCue={0,0,0,1},name=freeAxesCue_BeamLine
	AppendToGizmo string="X(BL)",strFont=font,name=stringX
// xxxxxxxxxxxxxxxxxxxxxxxxxxx
//	ModifyGizmo modifyObject=stringX property={Clipped,0}
	AppendToGizmo string="Y(BL)",strFont=font,name=stringY
// xxxxxxxxxxxxxxxxxxxxxxxxxxx
//	ModifyGizmo modifyObject=stringY property={Clipped,0}
	AppendToGizmo string="Z(BL)",strFont=font,name=stringZ
// xxxxxxxxxxxxxxxxxxxxxxxxxxx
//	ModifyGizmo modifyObject=stringZ property={Clipped,0}
	AppendToGizmo cylinder={0,0.2,1,25,25},name=cylinderArrow
	AppendToGizmo attribute lineWidth=2, name=lineWidthCross
	AppendToGizmo attribute color={1,0,0,1},name=colorXred
	AppendToGizmo attribute color={0,1,0,1},name=colorYgreen
	AppendToGizmo attribute color={0,0,1,1},name=colorZblue
	ModifyGizmo setDisplayList=0, opName=pushAttribute0, operation=pushAttribute, data=5
	ModifyGizmo setDisplayList=1, attribute=lineWidthCross
	ModifyGizmo setDisplayList=2, opName=pushMatrix0, operation=pushMatrix
	ModifyGizmo setDisplayList=3, opName=rotate45, operation=rotate, data={-45,1,0,0}
	ModifyGizmo setDisplayList=4, object=freeAxesCue_BeamLine
	ModifyGizmo setDisplayList=5, opName=pushMatrix_X, operation=pushMatrix
	ModifyGizmo setDisplayList=6, opName=translateX, operation=translate, data={1,0,0}
	ModifyGizmo setDisplayList=7, opName=scaleX, operation=scale, data={0.07,0.07,0.07}
	ModifyGizmo setDisplayList=8, attribute=colorXred
	ModifyGizmo setDisplayList=9, opName=rotateX, operation=rotate, data={-90,0,1,0}
	ModifyGizmo setDisplayList=10, object=stringX
	ModifyGizmo setDisplayList=11, object=cylinderArrow
	ModifyGizmo setDisplayList=12, opName=popMatrix_X, operation=popMatrix
	ModifyGizmo setDisplayList=13, opName=pushMatrix_Y, operation=pushMatrix
	ModifyGizmo setDisplayList=14, opName=translateY, operation=translate, data={0,1,0}
	ModifyGizmo setDisplayList=15, attribute=colorYgreen
	ModifyGizmo setDisplayList=16, opName=scaleY, operation=scale, data={0.07,0.07,0.07}
	ModifyGizmo setDisplayList=17, opName=rotateY, operation=rotate, data={90,1,0,0}
	ModifyGizmo setDisplayList=18, object=stringY
	ModifyGizmo setDisplayList=19, object=cylinderArrow
	ModifyGizmo setDisplayList=20, opName=popMatrix_Y, operation=popMatrix
	ModifyGizmo setDisplayList=21, opName=pushMatrix_Z, operation=pushMatrix
	ModifyGizmo setDisplayList=22, opName=translateZ, operation=translate, data={0,0,1}
	ModifyGizmo setDisplayList=23, attribute=colorZblue
	ModifyGizmo setDisplayList=24, opName=scaleZ, operation=scale, data={0.07,0.07,0.07}
	ModifyGizmo setDisplayList=25, opName=rotateZ, operation=rotate, data={180,1,0,0}
	ModifyGizmo setDisplayList=26, object=stringZ
	ModifyGizmo setDisplayList=27, object=cylinderArrow
	ModifyGizmo setDisplayList=28, opName=popMatrix_Z, operation=popMatrix
	ModifyGizmo setDisplayList=29, opName=popMatrix0, operation=popMatrix
	ModifyGizmo setDisplayList=30, opName=popAttribute0, operation=popAttribute
	ModifyGizmo currentGroupObject="::"
#endif
	// ************************* Group Object End *******************

	return groupName
End


Function setGizmoAxisLabels(xlabel,ylabel,zlabel)
	String xlabel,ylabel,zlabel
	if (itemsInList(WinList("*",";","WIN:"+num2istr(GIZMO_WIN_BIT)))==0)
		return 1
	endif

#if (IgorVersion()<7)
	Execute "GetGizmo objectItemExists=axes0"
	NVAR V_Flag=V_Flag
	if (!V_flag)
		Execute "AppendToGizmo Axes=boxAxes,name=axes0"
	endif
	String cmd
#else
	GetGizmo objectItemExists=axes0
	if (!V_flag)
		AppendToGizmo Axes=boxAxes,name=axes0
	endif
#endif

#if (IgorVersion()<7)
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
#else
	if (stringmatch(xlabel,"empty"))
		ModifyGizmo ModifyObject=axes0, objectType=axes property={0,axisLabel,0}
	elseif (strlen(xlabel))
		ModifyGizmo ModifyObject=axes0, objectType=axes property={0,ticks,3}
		ModifyGizmo ModifyObject=axes0, objectType=axes property={0,axisLabel,1}
		ModifyGizmo ModifyObject=axes0, objectType=axes property={0,axisLabelText,xlabel}
		ModifyGizmo ModifyObject=axes0, objectType=axes property={0,axisLabelDistance,0.05}
		ModifyGizmo ModifyObject=axes0, objectType=axes property={0,axisLabelScale,0.5}
	endif

	if (stringmatch(ylabel,"empty"))
		ModifyGizmo ModifyObject=axes0, objectType=axes property={1,axisLabel,0}
	elseif (strlen(ylabel))
		ModifyGizmo ModifyObject=axes0, objectType=axes property={1,ticks,3}
		ModifyGizmo ModifyObject=axes0, objectType=axes property={1,axisLabel,1}
		ModifyGizmo ModifyObject=axes0, objectType=axes property={1,axisLabelText,ylabel}
		ModifyGizmo ModifyObject=axes0, objectType=axes property={1,axisLabelDistance,0.05}
		ModifyGizmo ModifyObject=axes0, objectType=axes property={1,axisLabelScale,0.5}
	endif

	if (stringmatch(zlabel,"empty"))
		ModifyGizmo ModifyObject=axes0, objectType=axes property={2,axisLabel,0}
	elseif (strlen(zlabel))
		ModifyGizmo ModifyObject=axes0, objectType=axes property={2,ticks,3}
		ModifyGizmo ModifyObject=axes0, objectType=axes property={2,axisLabel,1}
		ModifyGizmo ModifyObject=axes0, objectType=axes property={2,axisLabelText,zlabel}
		ModifyGizmo ModifyObject=axes0, objectType=axes property={2,axisLabelDistance,0.2}
		ModifyGizmo ModifyObject=axes0, objectType=axes property={2,axisLabelScale,0.5}
	endif
#endif
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
#if (IgorVersion()<7)
		Execute "AppendToGizmo group,name="+groupName
		Execute "ModifyGizmo currentGroupObject=\""+groupName+"\""
		Execute "ModifyGizmo setDisplayList=0, opName=enable0, operation=enable, data=12288"	// enable all of the 6 clip planes
		Execute "ModifyGizmo setDisplayList=1, opName=enable1, operation=enable, data=12289"
		Execute "ModifyGizmo setDisplayList=2, opName=enable2, operation=enable, data=12290"
		Execute "ModifyGizmo setDisplayList=3, opName=enable3, operation=enable, data=12291"
		Execute "ModifyGizmo setDisplayList=4, opName=enable4, operation=enable, data=12292"
		Execute "ModifyGizmo setDisplayList=5, opName=enable5, operation=enable, data=12293"
		Execute "ModifyGizmo currentGroupObject=\"::\""
#else
		AppendToGizmo group,name=groupName
		ModifyGizmo currentGroupObject=groupName
		ModifyGizmo setDisplayList=0, opName=enable0, operation=enable, data={12288,}	// enable all of the 6 clip planes
		ModifyGizmo setDisplayList=1, opName=enable1, operation=enable, data={12289,}
		ModifyGizmo setDisplayList=2, opName=enable2, operation=enable, data={12290,}
		ModifyGizmo setDisplayList=3, opName=enable3, operation=enable, data={12291,}
		ModifyGizmo setDisplayList=4, opName=enable4, operation=enable, data={12292,}
		ModifyGizmo setDisplayList=5, opName=enable5, operation=enable, data={12293,}
		ModifyGizmo currentGroupObject="::"
#endif
		// ************************* Group Object End *******************

#if (IgorVersion()<7)
		Execute "ModifyGizmo userString={"+groupName+",\"\"}"	// there are no clip planes yet defined, so make empty
#else
		ModifyGizmo userString={$groupName,""}						// there are no clip planes yet defined, so make empty
#endif
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
#if (IgorVersion()<7)
	Execute "GetGizmo/Z userString="+groupName
	String keyVals=StrVarOrDefault("S_GizmoUserString","")
	KillStrings/Z S_GizmoUserString
#else
	GetGizmo/Z userString=groupName
	String keyVals=S_GizmoUserString
#endif

	Variable err=0
	strswitch(action)
		case "delete":
		case "remove":
#if (IgorVersion()<7)
			Execute "RemoveFromGizmo/Z object="+groupName
			Execute "ModifyGizmo userString={"+groupName+",\"\"}"	// there are no clip planes, so make empty
#else
			RemoveFromGizmo/Z object=$groupName
			ModifyGizmo userString={$groupName,""}						// there are no clip planes, so make empty
#endif
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
#if (IgorVersion()<7)
				Execute "ModifyGizmo currentGroupObject=\""+groupName+"\""
				Execute "RemoveFromGizmo/Z displayItem="+clipPlane
				Execute "ModifyGizmo currentGroupObject=\"::\""
#else
				ModifyGizmo currentGroupObject=groupName
				RemoveFromGizmo/Z displayItem=$clipPlane
				ModifyGizmo currentGroupObject="::"
#endif
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
#if (IgorVersion()<7)
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
#else
				Variable plane, A,B,C,D
				if (strsearch(action,"-X",0,2)>=0)
					plane = 0
					A=1; B=0; C=0		//		dataStr = "{0,1,0,0,"
					clipBox *= -1
				elseif (strsearch(action,"+X",0,2)>=0)
					plane = 1
					A=-1; B=0; C=0		//		dataStr = "{1,-1,0,0,"
				elseif (strsearch(action,"-Y",0,2)>=0)
					plane = 2
					A=0; B=1; C=0		//		dataStr = "{2,0,1,0,"
					clipBox *= -1
				elseif (strsearch(action,"+Y",0,2)>=0)
					plane = 3
					A=0; B=-1; C=0		//		dataStr = "{3,0,-1,0,"
				elseif (strsearch(action,"-Z",0,2)>=0)
					plane = 4
					A=0; B=0; C=1		//		dataStr = "{4,0,0,1,"
					clipBox *= -1
				elseif (strsearch(action,"+Z",0,2)>=0)
					plane = 5
					A=0; B=0; C=-1		//		dataStr = "{5,0,0,-1,"
				endif
//				dataStr += num2str(clipBox)+"}"

#endif

#if (IgorVersion()<7)
				Execute "ModifyGizmo currentGroupObject=\""+groupName+"\""
				Execute "GetGizmo displayNameList"
				String names=StrVarOrDefault("S_DisplayNames","")
				KillStrings/Z S_DisplayNames
#else
				ModifyGizmo currentGroupObject=groupName
				GetGizmo displayNameList
				String names=S_DisplayNames
#endif
				if (WhichListItem(clipPlane,names)<0)				// not present, add it
#if (IgorVersion()<7)
					Execute "ModifyGizmo setDisplayList=-1, opName="+clipPlane+", operation=ClipPlane, data="+dataStr
#else
					ModifyGizmo setDisplayList=-1, opName=$clipPlane, operation=ClipPlane, data={plane,A,B,C,D}
#endif
				else													// just modify existing
#if (IgorVersion()<7)
					Execute "ModifyGizmo/Z opName="+clipPlane+", operation=ClipPlane, data="+dataStr
#else
					ModifyGizmo/Z opName=$clipPlane, operation=ClipPlane, data={plane,A,B,C,D}
#endif
				endif
#if (IgorVersion()<7)
				Execute "ModifyGizmo currentGroupObject=\"::\""
#else
				ModifyGizmo currentGroupObject="::"
#endif
			endif
#if (IgorVersion()<7)
			Execute "ModifyGizmo userString={"+groupName+",\""+keyVals+"\"}"
			Execute "ModifyGizmo compile"
#else
			ModifyGizmo userString={$groupName,keyVals}
			ModifyGizmo compile
#endif
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
	gizmo = SelectString(ParamIsDefault(gizmo),gizmo,"")

	if (strlen(gizmo))
		if (WinType(gizmo)!=GIZMO_WIN_TYPE)
			return -1
		endif
	endif

#if (IgorVersion()<7)
	String Nswitch = SelectString(strlen(gizmo),"","/N="+gizmo)
	Execute "GetGizmo/Z"+Nswitch+" displayList"
	String displayList=StrVarOrDefault("S_DisplayList","")
	KillWaves/Z TW_DisplayList
	KillStrings/Z S_DisplayList
#else
	GetGizmo/Z/N=$gizmo displayList
	String displayList=S_DisplayList
	KillWaves/Z TW_DisplayList
#endif

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
	gizmo = SelectString(ParamIsDefault(gizmo),gizmo,"")

	if (strlen(gizmo))
		if (WinType(gizmo)!=GIZMO_WIN_TYPE)
			return ""
		endif
	endif
	exclude = SelectString(ParamIsDefault(exclude),exclude,"")

#if (IgorVersion()<7)
	String Nswitch = SelectString(strlen(gizmo),"","/N="+gizmo)
	Execute "GetGizmo/Z"+Nswitch+" objectList"
	KillStrings/Z S_gizmoObjectList
#else
	GetGizmo/Z/N=$gizmo objectList
#endif
	Wave/T TW_gizmoObjectList=TW_gizmoObjectList

	String item, objects="", objName
	String typeSep=type+SelectString(stringmatch(type,"group"),"=",",")
	Variable i, N=DimSize(TW_gizmoObjectList,0)
	for (i=0;i<N;i+=1)
		item = TW_gizmoObjectList[i]
		if (strsearch(item,"AppendToGizmo "+typeSep,0,2)==1)
			objName = TrimFrontBackWhiteSpace(StringByKey("name",item,"=",","))
			if (WhichListItem(objName,exclude)<0)
				objects += objName+";"
			endif
		endif
	endfor
	KillWaves/Z TW_gizmoObjectList
	return objects
End


Function/T GizmoListScatterWaves([gizmo])		// get list of all scatter plots, except 'scatterMarker0, scatterMarker1, ...
	String gizmo									// optional name of gizmo
	gizmo = SelectString(ParamIsDefault(gizmo),gizmo,"")
	if (strlen(gizmo))
		if (WinType(gizmo)!=GIZMO_WIN_TYPE)
			return ""
		endif
	endif

#if (IgorVersion()<7)
	String Nswitch = SelectString(strlen(gizmo),"","/N="+gizmo)
	Execute "ModifyGizmo"+Nswitch+" stopRotation"
	Execute "GetGizmo"+Nswitch+" displayList"		// list of all displayed objects
	Execute "GetGizmo"+Nswitch+" objectList"		// find name of wave with marker position
	KIllStrings/Z S_DisplayList, S_gizmoObjectList
#else
	ModifyGizmo/N=$gizmo stopRotation
	GetGizmo/N=$gizmo displayList			// list of all displayed objects
	GetGizmo/N=$gizmo objectList				// find name of wave with marker position
#endif
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
	return keyList
End


Function/T GizmoListIsoSurfaceWaves([gizmo])		// get list of all iso surface waves
	String gizmo									// optional name of gizmo
	gizmo = SelectString(ParamIsDefault(gizmo),gizmo,"")
	if (strlen(gizmo))
		if (WinType(gizmo)!=GIZMO_WIN_TYPE)
			return ""
		endif
	endif

#if (IgorVersion()<7)
	String Nswitch = SelectString(strlen(gizmo),"","/N="+gizmo)
	Execute "ModifyGizmo"+Nswitch+" stopRotation"
	Execute "GetGizmo"+Nswitch+" displayList"		// list of all displayed objects
	Execute "GetGizmo"+Nswitch+" objectList"		// find name of wave with marker position
	KIllStrings/Z S_DisplayList, S_gizmoObjectList
#else
	ModifyGizmo/N=$gizmo stopRotation
	GetGizmo/N=$gizmo displayList			// list of all displayed objects
	GetGizmo/N=$gizmo objectList				// find name of wave with marker position
#endif
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
	return keyList
End


Function/T GizmoListSurfaceWaves([gizmo])		// get list of all Surface plots on gizmo
	String gizmo									// optional name of gizmo
	gizmo = SelectString(ParamIsDefault(gizmo),gizmo,"")
	if (strlen(gizmo))
		if (WinType(gizmo)!=GIZMO_WIN_TYPE)
			return ""
		endif
	endif

#if (IgorVersion()<7)
	String Nswitch = SelectString(strlen(gizmo),"","/N="+gizmo)
	Execute "ModifyGizmo"+Nswitch+" stopRotation"
	Execute "GetGizmo"+Nswitch+" objectList"		// find name of wave with marker position
//	Execute "GetGizmo"+Nswitch+" displayList"		// list of all displayed objects
//	KIllStrings/Z S_DisplayList
//	KillWaves/Z TW_DisplayList
	KIllStrings/Z S_gizmoObjectList
#else
	ModifyGizmo/N=$gizmo stopRotation
	GetGizmo/N=$gizmo objectList				// find name of wave with marker position
#endif
	Wave/T TW_gizmoObjectList=TW_gizmoObjectList

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
	KillWaves/Z TW_gizmoObjectList
	return keyList
End


Function/T GetGizmoBoxDisplayed(gizmo)		// returns list with XYZ range of gizmo in USER units
	String gizmo				// use empty string for top gizmo
	if (itemsInList(WinList("*",";","WIN:"+num2istr(GIZMO_WIN_BIT)))==0)
		return ""
	endif

#if (IgorVersion()<7)
	//	String Nstr=SelectString(strlen(gizmo),"","/N="+gizmo)
	String Nswitch = SelectString(strlen(gizmo),"","/N="+gizmo)
	Execute "GetGizmo"+Nswitch+" userBoxLimits"
	NVAR GizmoBoxXmin=GizmoBoxXmin, GizmoBoxXmax=GizmoBoxXmax
	NVAR GizmoBoxYmin=GizmoBoxYmin, GizmoBoxYmax=GizmoBoxYmax
	NVAR GizmoBoxZmin=GizmoBoxZmin, GizmoBoxZmax=GizmoBoxZmax
	Variable Xlo=GizmoBoxXmin,Xhi=GizmoBoxXmax,Ylo=GizmoBoxYmin,Yhi=GizmoBoxYmax,Zlo=GizmoBoxZmin,Zhi=GizmoBoxZmax
	KillVariables/Z GizmoBoxXmin, GizmoBoxXmax, GizmoBoxYmin, GizmoBoxYmax, GizmoBoxZmin, GizmoBoxZmax
#else
	GetGizmo/N=$gizmo userBoxLimits
	Variable Xlo=GizmoBoxXmin,Xhi=GizmoBoxXmax,Ylo=GizmoBoxYmin,Yhi=GizmoBoxYmax,Zlo=GizmoBoxZmin,Zhi=GizmoBoxZmax
#endif
	//	printf "GizmoBox = [%g, %g] [%g, %g] [%g, %g]\r",Xlo,Xhi, Ylo,Yhi, Zlo,Zhi

	if (numtype(Xlo+Xhi+Ylo+Yhi+Zlo+Zhi))
#if (IgorVersion()<7)
		Execute "GetGizmo"+Nswitch+" dataLimits"
		NVAR GizmoXmin=GizmoXmin, GizmoXmax=GizmoXmax
		NVAR GizmoYmin=GizmoYmin, GizmoYmax=GizmoYmax
		NVAR GizmoZmin=GizmoZmin, GizmoZmax=GizmoZmax
		Xlo=GizmoXmin;Xhi=GizmoXmax; Ylo=GizmoYmin;Yhi=GizmoYmax; Zlo=GizmoZmin;Zhi=GizmoZmax
		KillVariables/Z GizmoXmin, GizmoXmax, GizmoYmin, GizmoYmax, GizmoZmin, GizmoZmax
#else
		GetGizmo/N=$gizmo dataLimits
		Xlo=GizmoXmin;Xhi=GizmoXmax; Ylo=GizmoYmin;Yhi=GizmoYmax; Zlo=GizmoZmin;Zhi=GizmoZmax
#endif
		//	printf "GizmoRange = [%g, %g] [%g, %g] [%g, %g]\r",Xlo,Xhi, Ylo,Yhi, Zlo,Zhi
	endif
	String out=""
	if (numtype(Xlo+Xhi+Ylo+Yhi+Zlo+Zhi)==0)
		sprintf out,"%g;%g;%g;%g;%g;%g",Xlo,Xhi,Ylo,Yhi,Zlo,Zhi
	endif
	return out
End



Function/T GetItemFromGizmoObject(win,object,name)
	String win
	String object
	String name

#if (IgorVersion()<7)
	String Nswitch = SelectString(strlen(win),"","/N="+win)
	Execute "GetGizmo"+Nswitch+"/Z objectList"
#else
	GetGizmo/N=$win/Z objectList
#endif
	Wave/T TW_gizmoObjectList=TW_gizmoObjectList
	KillStrings/Z S_gizmoObjectList
	if (!WaveExists(TW_gizmoObjectList))
		return ""
	endif

	String line, str, findAppend, findModify, out=""
	sprintf findAppend, "*AppendToGizmo %s=*,name=%s", name,object
	sprintf findModify, "*ModifyGizmo modifyObject=%s,*%s*", object,name

	Variable i0,i1,i, N=DimSize(TW_gizmoObjectList,0)
	for (i=0;i<N;i+=1)
		line = TW_gizmoObjectList[i]
		if (StringMatch(line,findAppend))
			sprintf str, "AppendToGizmo %s=",name
			i0 = strsearch(line,str,0,2)+strlen(str)
			if (i0<0)
				return ""
			endif
			i1 = strsearch(line,",",i0,2)
			out = line[i0,i1-1]
			break
		elseif (StringMatch(line,findModify))
			i0 = strsearch(line,",property={",0,2)
			if (i0<1)
				return ""
			endif
			i0 += strlen(",property={")
			i0 = strsearch(line,name,i0,2) + strlen(name) + 1
			i1 = strsearch(line,"}",i0,2)
			out = TrimFrontBackWhiteSpace(line[i0,i1-1])
			break
		endif

	endfor
	KillWaves/Z TW_gizmoObjectList
	return out
End




Function/WAVE GizmoObjectWave(objectName,objType,[gizmo])
	// get the wave associated with objectName
	String objectName					// name of object, e.g. "scatter0", "CubeCornersScatter", ...
	String objType						// object type, "Scatter", "Surface", "voxelgram", ...
	String gizmo						// optional name of gizmo
	gizmo = SelectString(ParamIsDefault(gizmo),gizmo,"")
	if (strlen(gizmo))
		if (WinType(gizmo)!=GIZMO_WIN_TYPE)
			return $""
		endif
	endif

	// get all ojects
#if (IgorVersion()<7)
	String Nswitch = SelectString(strlen(gizmo),"","/N="+gizmo)
	Execute "GetGizmo"+Nswitch+" objectList"
#else
	GetGizmo/N=$gizmo objectList
#endif
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

#if (IgorVersion()<7)
		String Nswitch = SelectString(strlen(win),"","/N="+win)
		Execute "GetGizmo"+Nswitch+"/Z objectList"
#else
		GetGizmo/N=$win/Z objectList
#endif
		Wave/T TW_gizmoObjectList=TW_gizmoObjectList

		for (i=0,j=-1; i<numpnts(TW_gizmoObjectList) && j<0; i+=1)
			j = strsearch(TW_gizmoObjectList[i],xyzName,0,2)
		endfor
		KillWaves/Z TW_gizmoObjectList
		return StringByKey("name",TW_gizmoObjectList[i-1],"=",",")
	endif

	String wnote=note(parametricXYZ)
	Wave sourceWave = $StringByKey("source", wnote,"=")
	if (!WaveExists(sourceWave))
		KillWaves/Z TW_gizmoObjectList
		return ""
	endif
	win = GizmosWithWave(sourceWave)
	if (strlen(win)<1)
		KillWaves/Z TW_gizmoObjectList
		return ""
	endif
	DoWindow/F $win					// bring gizmo to front, and then append the surface
	DoUpdate

	// Find unsued name for surface, starting with "ParametricSurf0"
#if (IgorVersion()<7)
	Nswitch = SelectString(strlen(win),"","/N="+win)
	Execute "GetGizmo"+Nswitch+"/Z objectNameList"
	String ObjectNames = StrVarOrDefault("S_ObjectNames","")
	KillStrings/Z S_ObjectNames
	KillWaves/Z TW_gizmoObjectList
#else
	GetGizmo/N=$win/Z objectNameList
	String ObjectNames = S_ObjectNames
#endif
	String rgbaName=StringByKey("RGBA",wnote,"="), surfaceName
	for (i=0,j=0; j>=0; i+=1)
		surfaceName = "ParametricSurf"+num2istr(i)
		j = WhichListItem(surfaceName, ObjectNames)
	endfor

#if (IgorVersion()<7)
	Execute "AppendToGizmo"+Nswitch+" Surface="+xyzName+",name="+surfaceName
	Execute "ModifyGizmo"+Nswitch+" ModifyObject="+surfaceName+" property={ srcMode,4}"
	Execute "ModifyGizmo"+Nswitch+" modifyObject="+surfaceName+" property={calcNormals,1}"
	if (exists(rgbaName)==1)
		Execute "ModifyGizmo"+Nswitch+" ModifyObject="+surfaceName+" property={ surfaceColorType,3}"
		Execute "ModifyGizmo"+Nswitch+" ModifyObject="+surfaceName+" property={ surfaceColorWave,"+rgbaName+"}"
	endif
	Execute "ModifyGizmo"+Nswitch+" setDisplayList=-1, object="+surfaceName
#else
	AppendToGizmo/N=$win Surface=$xyzName,name=$surfaceName
	ModifyGizmo/N=$win ModifyObject=$surfaceName objectType=Surface, property={ srcMode,4}
	ModifyGizmo/N=$win modifyObject=$surfaceName objectType=Surface, property={calcNormals,1}
	if (exists(rgbaName)==1)
		ModifyGizmo/N=$win ModifyObject=$surfaceName objectType=Surface, property={ surfaceColorType,3}
		ModifyGizmo/N=$win ModifyObject=$surfaceName objectType=Surface, property={ surfaceColorWave,$rgbaName}
	endif
	ModifyGizmo/N=$win setDisplayList=-1, object=$surfaceName
#endif
	return surfaceName
End



// Maybe this should not be here:
Static Function/T GizmosWithWave(scatter)
	Wave scatter

	if (!WaveExists(scatter))
		return ""
	endif
	String list = WinList("*", ";", "WIN:"+num2istr(GIZMO_WIN_BIT) )	// list of possible gizmos
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



#if (IgorVersion()>=7)
Menu "Gizmo"
	"-"
	"Put Cube Corners on Gizmo", ActivateCornerCubeOnGizmo(forceCalc=1)
	"De-Activate Cube Corners on Gizmo", DeActivateCornerCubeOnGizmo()
	"-"
End

#if (strlen(WinList("microGeometryN.ipf", ";","WIN:128")))
Menu "Gizmo"
	"Gizmo X-H plane", ModifyGizmo stopRotation ; ModifyGizmo SETQUATERNION={0,0,0,1}												// X-H
	"Gizmo X-F plane", ModifyGizmo stopRotation ; ModifyGizmo SETQUATERNION={0.707106781186547,0,0,0.707106781186547}	// X-F
	"Gizmo H-F plane", ModifyGizmo stopRotation ; ModifyGizmo SETQUATERNION={0.5,0.5,0.5,0.5}									// H-F
	"Gizmo Beamline", ModifyGizmo stopRotation ; ModifyGizmo SETQUATERNION={-0.270598,0.653282,-0.270598,0.653281}		// ModifyGizmo euler={0,-90,45}
	"-"
End
#else
Menu "Gizmo"
	"Gizmo X-Y plane [along beam, Z-axis]", ModifyGizmo stopRotation ; ModifyGizmo SETQUATERNION={0.0,0.0,0.0,1.0}
	"Gizmo Y-Z plane [side view, X-axis]", ModifyGizmo stopRotation ; ModifyGizmo SETQUATERNION={0.5,0.5,0.5,0.5}
	"Gizmo X-Z plane [-Y axis]", ModifyGizmo stopRotation ; ModifyGizmo SETQUATERNION={0.707107,0.0,0.0,0.707107}
	"Gizmo X-Z plane [top view]", ModifyGizmo stopRotation ; ModifyGizmo SETQUATERNION={-0.707107,0.0,0.0,0.707107}
	"Gizmo 111 vertical [0-11 to right]", ModifyGizmo stopRotation ; ModifyGizmo SETQUATERNION={0.060003,-0.540625,-0.455768,0.704556}
	//	ModifyGizmo euler = {45, acos(1/sqrt(3))*180/PI, 90}
	"-"
End
#endif
#endif



#if (IgorVersion()<7)
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
		Execute/Q/Z "GizmoMenu AppendItem={JZTr0,\"Gizmo X-Y plane [along beam, Z-axis]\", \"ModifyGizmo stopRotation ; ModifyGizmo SETQUATERNION={0.0,0.0,0.0,1.0}\"}"
		Execute/Q/Z "GizmoMenu AppendItem={JZTr1,\"Gizmo Y-Z plane [side view, X-axis]\", \"ModifyGizmo stopRotation ; ModifyGizmo SETQUATERNION={0.5,0.5,0.5,0.5}\"}"
		Execute/Q/Z "GizmoMenu AppendItem={JZTr2,\"Gizmo X-Z plane [-Y axis]\", \"ModifyGizmo stopRotation ; ModifyGizmo SETQUATERNION={0.707107,0.0,0.0,0.707107}\"}"
		Execute/Q/Z "GizmoMenu AppendItem={JZTr3,\"Gizmo X-Z plane [top view]\", \"ModifyGizmo stopRotation ; ModifyGizmo SETQUATERNION={-0.707107,0.0,0.0,0.707107}\"}"
		Execute/Q/Z "GizmoMenu AppendItem={JZTr4,\"Gizmo 111 vertical [0-11 to right]\", \"ModifyGizmo stopRotation ; ModifyGizmo SETQUATERNION={0.060003,-0.540625,-0.455768,0.704556}\"}"
		//	ModifyGizmo euler = {45, acos(1/sqrt(3))*180/PI, 90}
	endif
	Execute/Q/Z "GizmoMenu AppendItem={JZT4,\"-\", \"\"}"
End
#else
Static Function InitGizmoUtilityGeneral()
End
#endif
