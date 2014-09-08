#pragma rtGlobals=3		// Use modern global access method.
#pragma version = 2.00
#pragma IgorVersion = 6.2
#pragma ModuleName=GZoomTrans
#include "GizmoUtility", version>=0.16


Static Function AfterFileOpenHook(refNum,file,pathName,type,creator,kind)
	Variable refNum, kind
	String file,pathName,type,creator
	if ((kind==1) || (kind==2))		// an experiment (packed or unpacked)
		GZoomTrans#InitGizmoZoomTranslate()
	endif
	return 0
End
Static Function IgorStartOrNewHook(IgorApplicationNameStr)
	String IgorApplicationNameStr
	GZoomTrans#InitGizmoZoomTranslate()
	return 0
End



//  ======================================================================================  //
//  ============================== Start of Zoom/Translate ===============================  //

Function MakeGizmoZoomTransPanel() : Panel
	if (ItemsInList(WinList("GizmoZoomTransPanel", ";", "WIN:64")))
 		DoWindow/F GizmoZoomTransPanel
		return 0
	endif
	NewDataFolder/O root:Packages
	NewDataFolder/O root:Packages:GizmoZoomTranslate

	NewPanel/K=1/W=(1094,50,1250,220)/N=GizmoZoomTransPanel  as "for Gizmo"
	SetDrawLayer UserBack
	SetDrawEnv fstyle= 1
	DrawText 7,19,"Zoom & Translate Gizmo"
	SetVariable zoomVal,pos={18,25},size={101,19},bodyWidth=65,proc=GZoomTrans#GizmoZoomTranslateSetVarProc,title="zoom"
	SetVariable zoomVal,fSize=12,limits={1e-12,inf,0.1},value= _NUM:1
	SetVariable xTranslate,pos={21,50},size={90,19},bodyWidth=70,proc=GZoomTrans#GizmoZoomTranslateSetVarProc,title="∆X"
	SetVariable xTranslate,fSize=12,value= _NUM:0
	SetVariable yTranslate,pos={21,73},size={90,19},bodyWidth=70,proc=GZoomTrans#GizmoZoomTranslateSetVarProc,title="∆Y"
	SetVariable yTranslate,fSize=12,value= _NUM:0
	SetVariable zTranslate,pos={21,95},size={90,19},bodyWidth=70,proc=GZoomTrans#GizmoZoomTranslateSetVarProc,title="∆Z"
	SetVariable zTranslate,fSize=12,value= _NUM:0
	Button autoButton,pos={28,120},size={90,20},proc=GZoomTrans#GizmoZoomTransButtonProc,title="full scale"
	Button cubifyButton,pos={28,143},size={90,20},proc=GZoomTrans#GizmoZoomTransButtonProc,title="cubify"
	SetWindow GizmoZoomTransPanel hook(update)=GZoomTrans#GizmoZoomTranslateUpdateHook
	return 0
End

Function showBox()
	if (itemsInList(WinList("*",";","WIN:4096"))==0)
		return 1
	endif

	String fldrSav= GetDataFolder(1)
	SetDataFolder root:Packages:GizmoZoomTranslate
	Execute "GetGizmo userBoxLimits"
	NVAR GizmoBoxXmin=GizmoBoxXmin, GizmoBoxXmax=GizmoBoxXmax
	NVAR GizmoBoxYmin=GizmoBoxYmin, GizmoBoxYmax=GizmoBoxYmax
	NVAR GizmoBoxZmin=GizmoBoxZmin, GizmoBoxZmax=GizmoBoxZmax
	printf "GizmoBox = [%g, %g] [%g, %g] [%g, %g]\r",GizmoBoxXmin, GizmoBoxXmax,GizmoBoxYmin, GizmoBoxYmax,GizmoBoxZmin, GizmoBoxZmax

	Variable GizmoXmin,GizmoXmax, GizmoYmin,GizmoYmax, GizmoZmin,GizmoZmax
	getGizmoDataLimits(GizmoXmin,GizmoXmax, GizmoYmin,GizmoYmax, GizmoZmin,GizmoZmax)
	printf "Gizmo Data = [%g, %g] [%g, %g] [%g, %g]\r",GizmoXmin, GizmoXmax,GizmoYmin, GizmoYmax,GizmoZmin, GizmoZmax

	Variable zoom,xc,yc,zc, xhw,yhw,zhw
	getGizmoCurrentZoomCent(zoom,xc,yc,zc, xhw,yhw,zhw)
	printf "center offset = {%g, %g, %g},   half-widths = {%g, %g, %g}\r",xc,yc,zc,xhw,yhw,zhw
	SetDataFolder fldrSav
End


Function translateGizmo(dx,dy,dz)
	Variable dx,dy,dz
	if (itemsInList(WinList("*",";","WIN:4096"))==0)
		return 1
	endif

	if (numtype(dx+dy+dz))
		Execute/Z/Q "ModifyGizmo autoScale"
		return 0
	endif

	Variable Xmin,Xmax, Ymin,Ymax, Zmin,Zmax
	getGizmoDataLimits(Xmin,Xmax, Ymin,Ymax, Zmin,Zmax)
	Variable XminB, XmaxB, YminB, YmaxB, ZminB, ZmaxB
	getGizmoCurrentBox(XminB, XmaxB, YminB, YmaxB, ZminB, ZmaxB)

	dx -= (XminB+XmaxB)/2 - (Xmin+Xmax)/2
	dy -= (YminB+YmaxB)/2 - (Ymin+Ymax)/2
	dz -= (ZminB+ZmaxB)/2 - (Zmin+Zmax)/2
	Execute "ModifyGizmo scalingOption=0"
	String cmd
	sprintf cmd "ModifyGizmo setOuterBox={%g,%g,%g,%g,%g,%g}",XminB+dx,XmaxB+dx, YminB+dy,YmaxB+dy, ZminB+dz,ZmaxB+dz
	Execute cmd
	Execute "ModifyGizmo scalingMode=8"
End


Function zoomGizmo(zoom)
	Variable zoom
	if (itemsInList(WinList("*",";","WIN:4096"))==0)
		return 1
	endif

	if (zoom<=0 || numtype(zoom))
		Execute/Z/Q "ModifyGizmo autoScale"
		return 0
	endif

	Variable Xmin,Xmax, Ymin,Ymax, Zmin,Zmax
	getGizmoDataLimits(Xmin,Xmax, Ymin,Ymax, Zmin,Zmax)

	Variable XminB, XmaxB, YminB, YmaxB, ZminB, ZmaxB
	getGizmoCurrentBox(XminB, XmaxB, YminB, YmaxB, ZminB, ZmaxB)

	Variable xc=(XmaxB+XminB)/2, yc=(YmaxB+YminB)/2, zc=(ZmaxB+ZminB)/2
	Variable xWidth=Xmax-Xmin, yWidth=Ymax-Ymin, zWidth=Zmax-Zmin
	Variable dx = zoom*xWidth/2, dy = zoom*yWidth/2, dz = zoom*zWidth/2

	Execute "ModifyGizmo scalingOption=0"
	String cmd
	sprintf cmd "ModifyGizmo setOuterBox={%g,%g,%g,%g,%g,%g}",xc-dx, xc+dx, yc-dy, yc+dy, zc-dz, zc+dz
	Execute cmd
	Execute "ModifyGizmo scalingMode=8"
	//	ModifyGizmo scalingMode=8
	//	ModifyGizmo ScalingMode=2
End


Static Function GizmoZoomTranslateUpdateHook(s)
	STRUCT WMWinHookStruct &s

	Variable hookResult = 0
	switch(s.eventCode)
		case 0:				// Activate
		case 6:				// resize
		case 12:			// moved
			Variable zoom,xc,yc,zc, v
			hookResult = !getGizmoCurrentZoomCent(zoom,xc,yc,zc,v,v,v)
			SetVariable zoomVal,value= _NUM:zoom
			SetVariable xTranslate,value= _NUM:xc
			SetVariable yTranslate,value= _NUM:yc
			SetVariable zTranslate,value= _NUM:zc
			break
	endswitch
	return hookResult		// 0 if nothing done, else 1
End
//
Static Function GizmoZoomTranslateSetVarProc(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
			Variable dval = sva.dval
			String ctrlName = sva.ctrlName
			break
		default:
			return 0
	endswitch

	Variable zoom,xc,yc,zc, xhw,yhw,zhw
	getGizmoCurrentZoomCent(zoom,xc,yc,zc,xhw,yhw,zhw)
	strswitch(ctrlName)
		case "zoomVal":
			zoomGizmo(dval)
			break
		case "xTranslate":
			translateGizmo(dval,yc,zc)
			break
		case "yTranslate":
			translateGizmo(xc,dval,zc)
			break
		case "zTranslate":
			translateGizmo(xc,yc,dval)
			break
	endswitch

	return 0
End
//
Static Function GizmoZoomTransButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	if (ba.eventCode != 2)	// mouse up event
		return 0
	elseif (itemsInList(WinList("*",";","WIN:4096"))<1)
		return 0
	endif

	if (stringmatch(ba.ctrlName,"autoButton"))
		Execute/Z/Q "ModifyGizmo scalingOption=0"
		Execute/Z/Q "ModifyGizmo autoScale"
		Variable zoom,xc,yc,zc, xhw,yhw,zhw
		getGizmoCurrentZoomCent(zoom,xc,yc,zc, xhw,yhw,zhw)
		SetVariable zoomVal,value= _NUM:zoom
		SetVariable xTranslate,value= _NUM:xc
		SetVariable yTranslate,value= _NUM:yc
		SetVariable zTranslate,value= _NUM:zc
	elseif (stringmatch(ba.ctrlName,"cubifyButton"))
		Variable Xlo,Xhi,Ylo,Yhi,Zlo,Zhi
		getGizmoCurrentBox(Xlo,Xhi,Ylo,Yhi,Zlo,Zhi)
		Variable hw = max(max(Xhi-Xlo,Yhi-Ylo),Zhi-Zlo)/2
		Variable dx,dy,dz
		dx = hw - (Xhi-Xlo)/2
		dy = hw - (Yhi-Ylo)/2
		dz = hw - (Zhi-Zlo)/2
		String cmd
		sprintf cmd "ModifyGizmo setOuterBox={%g,%g,%g,%g,%g,%g}",Xlo-dx,Xhi+dx, Ylo-dy,Yhi+dy, Zlo-dz,Zhi+dz
		Execute cmd
		Execute "ModifyGizmo scalingMode=8"
	endif

	return 0
End
//
Static Function getGizmoCurrentZoomCent(zoom,xc,yc,zc,xhw,yhw,zhw)
	Variable &zoom, &xc,&yc,&zc, &xhw,&yhw,&zhw
	if (itemsInList(WinList("*",";","WIN:4096"))==0)
		zoom = NaN
		xc = NaN
		yc = NaN
		zc = NaN
		return 1
	endif

	Variable Xmin,Xmax, Ymin,Ymax, Zmin,Zmax
	getGizmoDataLimits(Xmin,Xmax, Ymin,Ymax, Zmin,Zmax)
	Variable XminB, XmaxB, YminB, YmaxB, ZminB, ZmaxB
	getGizmoCurrentBox(XminB, XmaxB, YminB, YmaxB, ZminB, ZmaxB)
	zoom = (XmaxB-XminB) / (Xmax-Xmin)
	xc = (XmaxB+XminB)/2 - (Xmax+Xmin)/2
	yc = (YmaxB+YminB)/2 - (Ymax+Ymin)/2
	zc = (ZmaxB+ZminB)/2 - (Zmax+Zmin)/2
	xhw = (XmaxB-XminB)/2
	yhw = (YmaxB-YminB)/2
	zhw = (ZmaxB-ZminB)/2
	return 0
End
//
Static Function getGizmoCurrentBox(Xlo,Xhi,Ylo,Yhi,Zlo,Zhi)
	Variable &Xlo,&Xhi,&Ylo,&Yhi,&Zlo,&Zhi

	String fldrSav= GetDataFolder(1)
	SetDataFolder root:Packages:GizmoZoomTranslate
	Execute "GetGizmo userBoxLimits"
	SetDataFolder fldrSav
	Xlo = NumVarOrDefault("root:Packages:GizmoZoomTranslate:GizmoBoxXmin",NaN)
	Xhi = NumVarOrDefault("root:Packages:GizmoZoomTranslate:GizmoBoxXmax",NaN)
	Ylo = NumVarOrDefault("root:Packages:GizmoZoomTranslate:GizmoBoxYmin",NaN)
	Yhi = NumVarOrDefault("root:Packages:GizmoZoomTranslate:GizmoBoxYmax",NaN)
	Zlo = NumVarOrDefault("root:Packages:GizmoZoomTranslate:GizmoBoxZmin",NaN)
	Zhi = NumVarOrDefault("root:Packages:GizmoZoomTranslate:GizmoBoxZmax",NaN)
	if (numtype(Xlo+Xhi+Ylo+Yhi+Zlo+Zhi))
		getGizmoDataLimits(Xlo,Xhi,Ylo,Yhi,Zlo,Zhi)	// use full data range if there  are no box values
	endif
End
//
Static Function getGizmoDataLimits(Xlo,Xhi,Ylo,Yhi,Zlo,Zhi)
	Variable &Xlo,&Xhi,&Ylo,&Yhi,&Zlo,&Zhi

	String fldrSav= GetDataFolder(1)
	SetDataFolder root:Packages:GizmoZoomTranslate
	Execute "GetGizmo dataLimits"
	SetDataFolder fldrSav
	Xlo = NumVarOrDefault("root:Packages:GizmoZoomTranslate:GizmoXmin",NaN)
	Xhi = NumVarOrDefault("root:Packages:GizmoZoomTranslate:GizmoXmax",NaN)
	Ylo = NumVarOrDefault("root:Packages:GizmoZoomTranslate:GizmoYmin",NaN)
	Yhi = NumVarOrDefault("root:Packages:GizmoZoomTranslate:GizmoYmax",NaN)
	Zlo = NumVarOrDefault("root:Packages:GizmoZoomTranslate:GizmoZmin",NaN)
	Zhi = NumVarOrDefault("root:Packages:GizmoZoomTranslate:GizmoZmax",NaN)
End

//  ======================================================================================  //
//  =============================== End of Zoom/Translate ================================  //



Static Function InitGizmoZoomTranslate()
	GizmoUtil#InitGizmoUtilityGeneral()
	Execute/Q/Z "GizmoMenu AppendItem={JZTz0,\"Zoom & Translate\", \"MakeGizmoZoomTransPanel()\"}"
	NewDataFolder/O root:Packages
	NewDataFolder/O root:Packages:GizmoZoomTranslate
End



