#pragma rtGlobals=3		// Use modern global access method.
#pragma version = 2.03
#pragma IgorVersion = 6.2
#pragma ModuleName=GClipPlanes
#include "GizmoUtility", version>=0.16

// with version 1.30, changed how clipping is done

Static Function AfterFileOpenHook(refNum,file,pathName,type,creator,kind)
	Variable refNum, kind
	String file,pathName,type,creator
	if ((kind==1) || (kind==2))		// an experiment (packed or unpacked)
		GClipPlanes#InitGizmoClipPlanes()
	endif
	return 0
End
Static Function IgorStartOrNewHook(IgorApplicationNameStr)
	String IgorApplicationNameStr
	GClipPlanes#InitGizmoClipPlanes()
	return 0
End


#if (IgorVersion()>=7)
Menu "Gizmo"
	"Clip Plane", MakeGizmoClipPlanePanel()
End
#endif


//  ======================================================================================  //
//  ================================ Start of Clip Plane =================================  //

Static StrConstant GizmoClipPlanePanelName = "GizmoClipPlanePanel"
Static StrConstant GizmoClipPlanePanelGroupName = "gizmoClipPlaneGroupPanel"

Function MakeGizmoClipPlanePanel() : Panel				// Creates the Clip Plane Panel, or brings it to the front
	if (ItemsInList(WinList(GizmoClipPlanePanelName,";","WIN:64")))
		DoWindow/F $GizmoClipPlanePanelName			// It already exists, bring to front (don't create a duplicate)
		return 0
	endif

	Variable left=1094,top=245,width=150
	NewPanel/K=1/W=(left,top,left+width,top+130)/N=$GizmoClipPlanePanelName  as "Clip Planes"
	SetDrawLayer UserBack
	SetDrawEnv fstyle= 1
	DrawText 7,19,"Clip Planes (Gizmo)"
	SetDrawEnv linethick= 3
	DrawLine 0,98,150,98

	SetVariable clipValue,pos={10,24},size={114,19},bodyWidth=80,title="value",fSize=12
	SetVariable clipValue,proc=GClipPlanes#GizmoClipPlaneSetVarProc,value= _NUM:NaN
	PopupMenu popupXYZ title="clip plane",pos={10,48},size={102,20},proc=GClipPlanes#GizmoClipPlanePopMenuProc
	PopupMenu popupXYZ ,mode=1,fSize=12, value=#"\"+X;-X;+Y;-Y;+Z;-Z\""
	SetVariable xyzStep,pos={19,73},size={103,19},bodyWidth=70,proc=GClipPlanes#GizmoClipPlaneSetVarProc,title="Æxyz"
	SetVariable xyzStep,fSize=12,value= _NUM:1
	Button removeClipPlaneButton,pos={3,104},size={130,21},proc=GClipPlanes#GizmoClipPlaneButtonProc,title="remove clip planes"
	SetWindow kwTopWin,hook(update)=GClipPlanes#GizmoClipPlanelUpdateHook
	return 0
End
//
Static Function GizmoClipPlanelUpdateHook(s)		// update values on Clip Panel. Needed since panel can server multiple Gizmos
	STRUCT WMWinHookStruct &s
	switch(s.eventCode)
		case 0:				// Activate
//		case 6:				// resize
//		case 12:			// moved
			GizmoClipPlanePanelUpdate(s.winName)	// up dates all values on Clip Plane Panel, to reflect status of top Gizmo
			return 1
			break
		endswitch
	return 0		// 0 if nothing done, else 1
End
//
Static Function GizmoClipPlaneSetVarProc(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	if (sva.eventCode!=1 && sva.eventCode!=2)		// mouse up andq Enter key
		return 0
	endif

	if (stringmatch(sva.ctrlName,"xyzStep"))
		SetLimitsClipPlanePanel(sva.win)
	elseif (stringmatch(sva.ctrlName,"clipValue"))
		ControlInfo/W=$(sva.win) popupXYZ
		ModifyGizmoClipPlaneGroup(GizmoClipPlanePanelGroupName,S_Value,sva.dval)
		GizmoClipPlanePanelUpdate(sva.win)
	endif
	return 0
End
//
Static Function GizmoClipPlanePopMenuProc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa
	if (pa.eventCode != 2)		// break out if NOT mouse-up
		return 0
	endif
#if (IgorVersion()<7)
	Execute "GetGizmo/Z userString="+GizmoClipPlanePanelGroupName
	String keyVals=StrVarOrDefault("S_GizmoUserString","")
	KillStrings/Z S_GizmoUserString
#else
	GetGizmo/Z userString=$GizmoClipPlanePanelGroupName
	String keyVals=S_GizmoUserString
#endif
	Variable val=NumberByKey(pa.popStr,keyVals)
	SetVariable clipValue,value= _NUM:val
	SetLimitsClipPlanePanel(pa.win)
	return 0
End
//
Static Function GizmoClipPlaneButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	if (ba.eventCode != 2)									// mouse up event
		return 0
	elseif (itemsInList(WinList("*",";","WIN:"+num2istr(GIZMO_WIN_BIT)))<1)	// no gizmos available, nothing to do
		return 0
	elseif (!stringmatch(ba.ctrlName,"removeClipPlaneButton"))
		return 1											// unknown buttong
	endif

	String gizName = StringFromList(0,WinList("*",";","WIN:"+num2istr(GIZMO_WIN_BIT)))		// name of top gizmo
	String groupList=GetGizmoObjects("group",gizmo=gizName)			// list of current group objects
	Variable clipON=(WhichListItem(GizmoClipPlanePanelGroupName,groupList)>=0)		// group exists for clip planes
	if (clipON)
		ModifyGizmoClipPlaneGroup(GizmoClipPlanePanelGroupName,"delete",NaN)
	else
		AddGizmoClipPlaneGroup(GizmoClipPlanePanelGroupName)
	endif
	GizmoClipPlanePanelUpdate(ba.win)
	return 0
End
//
Static Function GizmoClipPlanePanelUpdate(win)	// update the Panel to reflect the current Gizmo
	String win

	String gizName = StringFromList(0,WinList("*",";","WIN:"+num2istr(GIZMO_WIN_BIT)))		// name of top gizmo
	String Nstr=SelectString(strlen(gizName),"","/N="+gizName)
	String groupList=""	
	if (strlen(gizName)>0)
		groupList = GetGizmoObjects("group",gizmo=gizName)				// list of current group objects
	endif
	Variable clipON=(WhichListItem(GizmoClipPlanePanelGroupName,groupList)>=0)		// group exists for clip planes
	if (!clipON)										// nothing here, disable controls and re-set removeClipPlaneButton
		Button removeClipPlaneButton, win=$win, title="add clip planes"
		Variable dval = 2
		SetVariable xyzStep,win=$win, disable=dval
		SetVariable ClipValue, win=$win, disable=dval,value= _NUM:NaN
		PopupMenu popupXYZ disable=dval
		return 0
	endif

	// otherwise clip planes group is present
#if (IgorVersion()<7)
	Execute "GetGizmo/Z"+Nstr+" userString="+GizmoClipPlanePanelGroupName
	String keyVals=StrVarOrDefault("S_GizmoUserString","")
	KillStrings/Z S_GizmoUserString
#else
	GetGizmo/Z/N=$gizName userString=$GizmoClipPlanePanelGroupName
	String keyVals=S_GizmoUserString
#endif
	Button removeClipPlaneButton, win=$win, title="remove clip planes"

	if (strlen(keyVals)<1)
#if (IgorVersion()<7)
		Execute "RemoveFromGizmo/Z displayItem="+GizmoClipPlanePanelGroupName
#else
		RemoveFromGizmo/Z displayItem=$GizmoClipPlanePanelGroupName
#endif
	else
		// check if GizmoClipPlanePanelGroupName is displayed, if not then display it
		Variable index = isGizmoObjectDisplayed(GizmoClipPlanePanelGroupName)
		if (index<0)
			index = IndexForClipPlanesInDisplayList()
#if (IgorVersion()<7)
			if (index>=0)
				Execute "ModifyGizmo insertDisplayList="+num2istr(index)+", object="+GizmoClipPlanePanelGroupName
			else
				Execute "ModifyGizmo setDisplayList=-1, object="+GizmoClipPlanePanelGroupName
			endif
#else
			if (index>=0)
				ModifyGizmo insertDisplayList=index, object=$GizmoClipPlanePanelGroupName
			else
				ModifyGizmo setDisplayList=-1, object=$GizmoClipPlanePanelGroupName
			endif
#endif
		endif
	endif

	dval = 0
	PopupMenu popupXYZ disable=dval
	ControlInfo/W=$win popupXYZ
	Variable val=NumberByKey(S_Value, keyVals)
	SetVariable ClipValue, win=$win, disable=dval,value= _NUM:val
	SetVariable xyzStep,win=$win, disable=dval
	return 1
End
//
Static Function SetLimitsClipPlanePanel(win)			// reset value and limits of clip value based on the clip plane
	String win					// panel window name

	ControlInfo/W=$win popupXYZ
	Variable i = 2*WhichListItem(S_Value[1], "X;Y;Z")
	String box=GetGizmoBoxDisplayed("")
	Variable lo=str2num(StringFromList(i,box)), hi=str2num(StringFromList(i+1,box))

	ControlInfo/W=$win xyzStep
	Variable delta = V_Value

	ControlInfo/W=$win clipValue
	Variable value = limit(V_Value,lo,hi)
	if (numtype(hi+lo+delta)==0)
		SetVariable clipValue limits={lo,hi,delta}
		if (!(value == V_Value))
			SetVariable clipValue, win=$win,value= _NUM:value
		endif
	endif
End
//
// use with:
//	Variable index = IndexForClipPlanesInDisplayList()
//	ModifyGizmo insertDisplayList=index, object=gizmoClipPlaneGroupPanel
//
// finds position in display list for the clip planes,  before any "real" objects, but after everything else
Static Function IndexForClipPlanesInDisplayList()	// put it right before the first "real" object, ignore groups
	String realObjects=""							//	list of "real" objects, isoSurface, scatter, sphere, ...
#if (IgorVersion()<7)
	Execute "GetGizmo objectList"				// find name of wave with marker position
	KillStrings/Z S_gizmoObjectList
#else
	GetGizmo objectList								// find name of wave with marker position
#endif
	Wave/T TW_gizmoObjectList=TW_gizmoObjectList

	String realTypes=LowerStr("isoSurface;Scatter;Path;surface;ribbon;voxelgram;line;triangle;quad;box;sphere;cylinder;disk;")
	String item, type, name
	Variable i,m
	for (m=0;m<numpnts(TW_gizmoObjectList);m+=1)
		item = TW_gizmoObjectList[m]
		i = strsearch(item,"AppendToGizmo ",0)
		if (i==1)
			item = item[14,Inf]
			i = strsearch(item,"=",0)
			type = item[1,i-1]
			i = strsearch(item,"name=",0)
			name = item[i+5,Inf]
			if (WhichListItem(LowerStr(type), realTypes)>=0)
				realObjects += name+";"
			endif
		endif
	endfor

#if (IgorVersion()<7)
	Execute "GetGizmo displayList"				// list of all displayed objects
	KillStrings/Z S_DisplayList
#else
	GetGizmo displayList							// list of all displayed objects
#endif
	Wave/T TW_DisplayList=TW_DisplayList
	Variable index=-1
	for (m=0;m<numpnts(TW_DisplayList) && index<0;m+=1)		// find first realObject in TW_DisplayList
		item = TW_DisplayList[m]
		i = strsearch(item,"object=",0,2)
		if (i>0)
			item = item[i+7,Inf]
			if (WhichListItem(item,realObjects)>=0)
				item = TW_DisplayList[m]
				i = strsearch(item," setDisplayList=",0,2)
				if (i>0)
					index = str2num(item[i+16,Inf])
				endif
			endif
		endif
	endfor
	KillWaves/Z TW_gizmoObjectList, TW_DisplayList
	return index
End


Function getClipPlaneValue(plane)	// Utility to get the current clip value for a plane.
	String plane					// plane id, one of {"-X", "+X",  "-Y", "+Y",  "-Z", "+Z"}
#if (IgorVersion()<7)
	Execute "GetGizmo/Z userString="+GizmoClipPlanePanelGroupName
	String keyVals=StrVarOrDefault("S_GizmoUserString","")
	KillStrings/Z S_GizmoUserString
#else
	GetGizmo/Z userString=$GizmoClipPlanePanelGroupName
	String keyVals=S_GizmoUserString
#endif
	return NumberByKey(plane,keyVals)
End

//  ======================================================================================  //
//  ================================= End of Clip Plane ==================================  //


Static Function InitGizmoClipPlanes()
	GizmoUtil#InitGizmoUtilityGeneral()
#if (IgorVersion()<7)
	Execute/Q/Z "GizmoMenu AppendItem={JZTc0,\"Clip Plane\", \"MakeGizmoClipPlanePanel()\"}"
#endif
End


