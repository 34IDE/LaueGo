#pragma rtGlobals=1		// Use modern global access method.
#pragma version = 1.42
// #include "vector-math"	// was needed for normalize()
Constant debugN = 0


Menu "Grains"
	"(Info"
	"Statistics of a Boundary",BoundaryStats($"",NaN,NaN,NaN)
	"Check for All Sigma Boundarys",checkAllSigmaBoundarys($"",NaN,NaN,NaN)
	"List Neighbors to a Grain",makeBorderList($"",NaN)
	"Angle Between two Grains",angleBetweenGrains($"",NaN,NaN)
	"   Angle Between two Voxels",AskAngleBetweenPointsSymReduced(NaN,NaN)
	"Find Number of Grains Above a Size", NumOfGrainsBiggerThan($"",NaN)
	"-"
	"Table of grain centers & orientations",Table_GrainOrientations($"")
	"Table of a multigrain 3d wave",Table_AllGrains($"")
	"Graph of the grain volumes", Graph_grain_Vol($"")
	"Table of XYZ indicies",Table_XYZindex($"")
	"Table of all grain statistics",Table_GrainStatistics($"")
	"Histogram of total angles",Graph_HistTotalAngle($"")
	"-"
	"(New Data"
	"Prepare\Read In data",PrepareData()
	"Calculate Grains from Neighbor Rotations",findGrainsByComparingNeighbors(NaN,NaN,0)
End
Menu "Load Waves"
	"-"
	"(Gizmo Grain Data"
	"Prepare\Read In data",PrepareData()
	"Calculate Grains from Neighbor Rotations",findGrainsByComparingNeighbors(NaN,NaN,0)
End


//Menu "Grains"
//	"Zoomer Panel for a Gizmo",PanelGizmoView()
//	"-"
//End

// when putting up a sphere, subtract (1,1,1) from it so that it matches the iso-surfaces
//	e.g. if I want a sphere at (2,3,4), put it at (1,2,3),  see MoveGizmoCursorProc() for example
//	angle ~ sqrt(axisX^2+axisY^2+axisZ^2)*180/PI/2 * 1.04386


Window PanelGizmoView() : Panel
	PauseUpdate; Silent 1		// building window...
	String gizWin = StringFromList(0,WinList("*",";","WIN:4096"))		// gizmo window for positioning
	if (strlen(WinList("Panel_GizmoViewer","","WIN:64"))>1)
		DoWindow/F Panel_GizmoViewer
		AutoPositionWindow /M=0/R=$gizWin Panel_GizmoViewer
		return
	endif
	initGrainGizmo()
	NewPanel /K=1 /W=(690,50,859,620)					// W=(left, top, right, bottom )
	DoWindow /C Panel_GizmoViewer
	AutoPositionWindow /M=0/R=$gizWin Panel_GizmoViewer	// position panel by top most gizmo

	SetDrawLayer UserBack
	SetDrawEnv xcoord= rel,ycoord= abs,linethick= 3
	DrawLine 0,292,1,292
	SetDrawEnv xcoord= rel,ycoord= abs,linethick= 3
	DrawLine 0,443,1,443
	SetDrawEnv xcoord= rel,ycoord= abs,linethick= 3
	DrawLine 0,536,1,536
	String micron = StrVarOrDefault("root:Packages:Grains:micron","µm")

	Button buttonHome,pos={22,99},size={50,20},proc=ButtonHomeProc,title="Home"
	Button buttonStopRotation,pos={85,99},size={50,20},proc=ButtonStopProc,title="Stop"
	Button buttonStopRotation help={"stp rotating"}
	Button buttonZoomIn,pos={10,3},size={65,20},proc=ButtonZoomProc,title="Zoom In"
	Button buttonZoomOut,pos={85,3},size={75,20},proc=ButtonZoomProc,title="Zoom Out"
	Button buttonXplus,pos={43,27},size={30,20},proc=ButtonTranslateProc,title="X+"
	Button buttonXminus,pos={83,27},size={30,20},proc=ButtonTranslateProc,title="X-"
	Button buttonYplus,pos={43,51},size={30,20},proc=ButtonTranslateProc,title="Y+"
	Button buttonYminus,pos={83,51},size={30,20},proc=ButtonTranslateProc,title="Y-"
	Button buttonZplus,pos={43,75},size={30,20},proc=ButtonTranslateProc,title="Z+"
	Button buttonZminus,pos={83,75},size={30,20},proc=ButtonTranslateProc,title="Z-"
	PopupMenu popupGrainDim,pos={34,123},size={100,20},proc=PopGrainProc,title="dim grain"
	PopupMenu popupGrainDim,mode=1,popvalue="none",value= #"\"none;\"+getGrainNameList()"
	PopupMenu popupGrainBrite,pos={22,147},size={111,20},proc=PopGrainProc,title="briten grain"
	PopupMenu popupGrainBrite,mode=1,popvalue="none",value= #"\"none;\"+getGrainNameList()"
	PopupMenu popupGrainRestore,pos={16,171},size={118,20},proc=PopGrainProc,title="restore grain"
	PopupMenu popupGrainRestore,mode=1,popvalue="none",value= #"\"none;\"+getGrainNameList()"
	Button buttonAllDim,pos={40,195},size={85,20},proc=ButtonGrainColorProc,title="Dim All"
	Button buttonAllSolid,pos={40,219},size={85,20},proc=ButtonGrainColorProc,title="All Solid"
	Button buttonAllRestore,pos={40,243},size={85,20},proc=ButtonGrainColorProc,title="Restore All"
	PopupMenu GizmoPopup,pos={48,267},size={64,20},proc=PopGizmoProc,fSize=12
	PopupMenu GizmoPopup,mode=1,popvalue="Gizmo",value= #"\"Gizmo;Show Axes;Hide Axes;Re-Compile;Show Info Window\""

	PopupMenu popupBndryGrain1,pos={19,299},size={87,20},proc=PopGrainBoundryProc,title="grain1"
	PopupMenu popupBndryGrain1,help={"select first grain of boundary, usually central grain"}
	PopupMenu popupBndryGrain1,mode=1,popvalue="none",value= #"\"none;\"+getGrainNameList()"
	PopupMenu popupBndryGrain2,pos={19,322},size={87,20},proc=PopGrainBoundryProc,title="grain2"
	PopupMenu popupBndryGrain2,help={"select secpmd grain of boundary"}
	PopupMenu popupBndryGrain2,mode=1,popvalue="none",value= #"\"none;\"+getGrainNameList()"
	SetVariable setRadius,pos={10,346},size={140,18},title="radius ("+micron+")",fSize=12,format="%.2f"
	SetVariable setRadius,limits={0,inf,.5},value= root:Packages:Grains:radiusBndry
	SetVariable setRadius,help={"distance from cursor to accept when finding grain boundary points ("+micron+")"}
	Button buttonCalcBndry,pos={11,367},size={136,20},proc=ButtonGrainBndryProc,title="Update Boundary"
	Button buttonCalcBndry,help={"calculate the grain boundary parameters"}
	Button buttonCenterBndry,pos={11,392},size={136,20},proc=ButtonGrainBndryProc,title="Center on Boundary"
	Button buttonCenterBndry,help={"re-center plot on the grain boundary"}
	PopupMenu popupScatterBndry,pos={10,417},size={143,20},proc=PopGrainBoundryPointsProc,title="boundary points"
	PopupMenu popupScatterBndry,help={"turns scatter boundary points on/off"}
	PopupMenu popupScatterBndry,fSize=12,mode=1,popvalue="Off",value= #"\"Off;On\""

	SetVariable Xcursor,pos={1,452},size={70,18},proc=MoveGizmoCursorProc,title="X"
	SetVariable Xcursor,fSize=12
	SetVariable Xcursor,limits={-10,inf,1},value= root:Packages:Grains:Xcursor
	SetVariable Ycursor,pos={1,472},size={70,18},proc=MoveGizmoCursorProc,title="Y"
	SetVariable Ycursor,fSize=12
	SetVariable Ycursor,limits={-10,inf,1},value= root:Packages:Grains:Ycursor
	SetVariable Zcursor,pos={1,492},size={70,18},proc=MoveGizmoCursorProc,title="Z"
	SetVariable Zcursor,fSize=12
	SetVariable Zcursor,limits={-10,inf,1},value= root:Packages:Grains:Zcursor
	Button CursorOn,pos={75,459},size={90,20},proc=ButtonGizmoCursorProc,title="Cursor ON"
	Button CursorOff,pos={75,482},size={90,20},proc=ButtonGizmoCursorProc,title="Cursor OFF"
	Button GetGrainNum,pos={11,513},size={90,20},proc=ButtonGizmoCursorProc,title="Get Grain #"
	SetVariable grainNumber,pos={108,513},size={50,18},proc=SetGrainNumProc,title=" "
	SetVariable grainNumber,help={"grain number works with 3d cursor"},fSize=12
	SetVariable grainNumber,fStyle=1
	SetVariable grainNumber,limits={0,inf,1},value= root:Packages:Grains:cursorGrainNum

	SetVariable explode,pos={19,544},size={122,18},proc=ExplodeProc,title="Explode"
	SetVariable explode,help={"set exploded view level [0,100]"},fSize=12,fStyle=1
	SetVariable explode,limits={0,100,1},value= root:Packages:Grains:explode
EndMacro

// These two fucntions are needed to control the cursor marker on the gizmo plot
Function MoveGizmoCursorProc(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName
	NVAR xx = root:Packages:Grains:Xcursor
	NVAR yy = root:Packages:Grains:Ycursor
	NVAR zz = root:Packages:Grains:Zcursor
	String cmd
	sprintf cmd, "ModifyGizmo opName=translateCursor, operation=translate, data={%g,%g,%g}",xx,yy,zz
	Execute cmd
	sprintf cmd, "ModifyGizmo opName=translateCursorBack, operation=translate, data={%g,%g,%g}",-xx,-yy,-zz
	Execute cmd

	// this next section makes sure that the range of the control is set to a bit larger than the size of the grains
	if (cmpstr(ctrlName,"Xcursor") || cmpstr(ctrlName,"Ycursor") || cmpstr(ctrlName,"Zcursor"))
		String gWaveName=StringFromList(0,getGrainWaveList())		// name of first grain wave in gizmo
		gWaveName = StringByKey("grainsWave", note($gWaveName),"=")
		if (exists(gWaveName)!=1)
			DoAlert 0, "unable to find 3d grains wave in ButtonGizmoCursorProc()"
			return 1
		endif
		Wave grains=$gWaveName
		Variable xlo,xhi, ylo,yhi, zlo,zhi
		xlo = DimOffset(grains,0)
		ylo = DimOffset(grains,1)
		zlo = DimOffset(grains,2)
		xhi = xlo + DimDelta(grains,0)*DimSize(grains,0)
		yhi = ylo + DimDelta(grains,1)*DimSize(grains,1)
		zhi = zlo + DimDelta(grains,2)*DimSize(grains,2)
		SetVariable Xcursor,limits={xlo-2,xhi+2,1}
		SetVariable Ycursor,limits={ylo-2,yhi+2,1}
		SetVariable Zcursor,limits={zlo-2,zhi+2,1}
	endif
End
//Function MoveGizmoCursorProc(ctrlName,varNum,varStr,varName) : SetVariableControl
//	String ctrlName
//	Variable varNum
//	String varStr
//	String varName
//	NVAR xx = root:Packages:Grains:Xcursor
//	NVAR yy = root:Packages:Grains:Ycursor
//	NVAR zz = root:Packages:Grains:Zcursor
//	String cmd
//	sprintf cmd, "ModifyGizmo opName=translateCursor, operation=translate, data={%g,%g,%g}",xx,yy,zz
//	Execute cmd
//	sprintf cmd, "ModifyGizmo opName=translateCursorBack, operation=translate, data={%g,%g,%g}",-xx,-yy,-zz
//	Execute cmd
//End


Function ButtonGizmoCursorProc(ctrlName) : ButtonControl
	String ctrlName
	Execute "GetGizmo gizmoName"		// check if a gizmo is up
	SVAR S_GizmoName=S_GizmoName
	if (strlen(S_GizmoName)<1)
		DoAlert 0, "No Gizmo Windows Open"
		return 1
	endif

	if (stringmatch(ctrlName,"CursorOn"))
		Execute "ModifyGizmo modifyObject=sphereCursor, property={radius,0.7}"
		return 0
	elseif (stringmatch(ctrlName,"CursorOff"))
		Execute "ModifyGizmo modifyObject=sphereCursor, property={radius,0.0}"
		return 0
	elseif (stringmatch(ctrlName,"GetGrainNum"))	// take the 3d cursor position,  and update the grain number in the panel
		String gWaveName=StringFromList(0,getGrainWaveList())		// name of first grain wave in gizmo
		gWaveName = StringByKey("grainsWave", note($gWaveName),"=")
		if (exists(gWaveName)!=1)
			DoAlert 0, "unable to find 3d grains wave in ButtonGizmoCursorProc()"
			return 1
		endif
		Wave grains=$gWaveName
		NVAR xx = root:Packages:Grains:Xcursor
		NVAR yy = root:Packages:Grains:Ycursor
		NVAR zz = root:Packages:Grains:Zcursor
		if (!NVAR_Exists(xx) || !NVAR_Exists(yy) || !NVAR_Exists(zz))
			DoAlert 0, "unable to retrieve cursor positions ButtonGizmoCursorProc()"
			return 1
		endif
		Variable grainNum = xyzPositionToGrainNum(grains,xx,yy,zz)
		if (exists("root:Packages:Grains:cursorGrainNum")==2)
			NVAR cursorGrainNum = root:Packages:Grains:cursorGrainNum
			cursorGrainNum = grainNum
		endif
		return 0
	endif
	DoAlert 0, "ERROR, ButtonGizmoCursorProc() called with ctrlName= '"+ctrlName+"'"
	return 1
End
//Function ButtonGizmoCursorProc(ctrlName) : ButtonControl
//	String ctrlName
//	if (stringmatch(ctrlName,"CursorOn"))
//		Execute "ModifyGizmo modifyObject=sphereCursor, property={radius,0.7}"
//	elseif (stringmatch(ctrlName,"CursorOff"))
//		Execute "ModifyGizmo modifyObject=sphereCursor, property={radius,0.0}"
//	else
//		DoAlert 0, "ERROR, ButtonGizmoCursorProc() called with ctrlName= '"+ctrlName+"'"
//	endif
//End

// take the 3d position (µm),  and return the grain number
Function xyzPositionToGrainNum(grains,xx,yy,zz)
	Wave grains						// 3d wave with all of the grain numbers (e.g. gNeighbor)
	Variable xx,yy,zz
	if (!WaveExists(grains) || numtype(xx+yy+zz))
		String grainName=NameOfWave(grains)
		Prompt grainName,"3d wave with all grains",popup,ListMultiGrainMats()
		xx = numtype(xx) ? 0 : xx
		yy = numtype(yy) ? 0 : yy
		zz = numtype(zz) ? 0 : zz
		String micron = StrVarOrDefault("root:Packages:Grains:micron","micron")
		Prompt xx,"X position in grain ("+micron+")"
		Prompt yy,"Y position in grain ("+micron+")"
		Prompt zz,"Z position in grain ("+micron+")"
		DoPrompt "choose grains",grainName,xx,yy,zz
		if (V_flag)
			return -1
		endif
		Wave grains = $grainName
	endif
	if (!WaveExists(grains))
		return -1
	endif
	Variable i,j,k,gNum=NaN
	i = (xx-DimOffset(grains,0))/DimDelta(grains,0)
	j = (yy-DimOffset(grains,1))/DimDelta(grains,1)
	k = (zz-DimOffset(grains,2))/DimDelta(grains,2)
	if (numtype(i+j+k)==0 && i>=0 && j>=0 && k>=0 && i<DimSize(grains,0) && j<DimSize(grains,1) && k<DimSize(grains,2))
		gNum = grains[i][j][k]
	endif
	if (ItemsInList(GetRTStackInfo(0))<2)
		String name=NameOfWave(grains)
		if (numtype(gNum))
			printf "   3D cursor at %s(%g)(%g)(%g)  ==  %s[%g][%g][%g]  is unidentified\r",name,xx,yy,zz,name,i,j,k
		else
			printf "   3D cursor at %s(%g)(%g)(%g)  ==  %s[%g][%g][%g]  is grain #%g\r",name,xx,yy,zz,name,i,j,k,gNum
		endif
	endif
	return gNum
End

Function SetGrainNumProc(ctrlName,grainNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable grainNum
	String varStr
	String varName
	if (!stringmatch(ctrlName,"grainNumber"))
		return 1
	endif

	// check that cursor variables exist
	if (exists("root:Packages:Grains:Xcursor")!=2 || exists("root:Packages:Grains:Ycursor")!=2 || exists("root:Packages:Grains:Zcursor")!=2)
		DoAlert 0, "XYZ cursor variables do not exist, do a 'initGrainGizmo()'"
		print "   do this ->		initGrainGizmo()"
		return 1
	endif
	String gWaveName=StringFromList(0,getGrainWaveList())		// name of first grain wave in gizmo
	gWaveName = StringByKey("grainsWave", note($gWaveName),"=")
	if (exists(gWaveName)!=1)
		DoAlert 0, "unable to find 3d grains wave in SetGrainNumProc()"
		return 1
	endif
	Wave grains=$gWaveName										// finally this is the 3d wave of all grains
	Wave com=$(GetWavesDataFolder(grains,2)+"_COM")			// wave with single index to com of each grain
	Variable index = com[grainNum]								// single index to com of grainNum
	if (numtype(index))
		DoAlert 0, "got back a point index of "+""+""
	endif
	String dataPath = GetWavesDataFolder(grains,1)+"raw:"
	Wave XX=$(dataPath+"XX")										// x,y,z positions read in from file
	Wave YY=$(dataPath+"YY")
	Wave ZZ=$(dataPath+"ZZ")
	if (!WaveExists(XX) || !WaveExists(YY) || !WaveExists(ZZ))
		DoAlert 0, "unable to find XX, YY, ZZ waves in SetGrainNumProc()"
		return 1
	endif
	NVAR Xcursor = root:Packages:Grains:Xcursor
	NVAR Ycursor = root:Packages:Grains:Ycursor
	NVAR Zcursor = root:Packages:Grains:Zcursor
	Xcursor = XX[index]											// set cursor xyz to com of grainNum
	Ycursor = YY[index]
	Zcursor = ZZ[index]
	MoveGizmoCursorProc("SetGrainNumProc",0,"","")				// move the cursor to this position
	ButtonGizmoCursorProc("CursorOn")	// make sure cursor is on
	return 0
End


Function PopGrainBoundryPointsProc(ctrlName,popNum,popStr) : PopupMenuControl
	String ctrlName
	Variable popNum
	String popStr

	if (!stringmatch(ctrlName,"popupScatterBndry"))
		return 1
	endif
	String cmd
	sprintf cmd, "ModifyGizmo ModifyObject=scatterBndry property={ size,%d}",(stringmatch(popStr,"On")
	Execute cmd
End

Function PopGizmoProc(ctrlName,popNum,popStr) : PopupMenuControl
	String ctrlName
	Variable popNum
	String popStr
	if (stringmatch(popStr,"Show Axes"))
//		Execute "ModifyGizmo showAxisCue=1"
		Execute "ModifyGizmo insertDisplayList=2, object=freeAxesCue0"
	elseif (stringmatch(popStr,"Hide Axes"))
		Execute "ModifyGizmo showAxisCue=0"
		Execute "RemoveFromGizmo displayItemNumber=2"
	elseif (stringmatch(popStr,"Re-Compile"))
		Execute "ModifyGizmo compile"
	elseif (stringmatch(popStr,"Show Info Window"))
		Execute "ModifyGizmo showInfo"
	else
		return 1
	endif
	PopupMenu GizmoPopup mode=1
End



Function ExplodeProc(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName

	if (!stringmatch(ctrlName,"explode") || varNum<0 || numtype(varNum))
		return 1
	endif
	Variable explode = varNum
	if (explode!=0)								// when exploded, hide dots and cursor
		Execute "ModifyGizmo modifyObject=sphereCursor, property={radius,0.0}"
		if (strlen(getGrainPropertyFromObject("scatterBndry","size")))
			Execute "ModifyGizmo ModifyObject=scatterBndry property={ size,0}"
		endif
	endif

	Variable xc0=NaN,yc0=NaN,zc0=NaN
	String gName, noteStr, gizmoNote = GetGizmoNote()
	gName = StringByKey("grainBndry1",gizmoNote,"=",":")
	Wave grain =$getGrainWaveFromName(gName)

	if (WaveExists(grain))
		noteStr = note(grain)
		xc0 = NumberByKey("grain_xc",noteStr,"=")		// center of explosion in 'grain'
		yc0 = NumberByKey("grain_yc",noteStr,"=")
		zc0 = NumberByKey("grain_zc",noteStr,"=")
	endif
	if (numtype(xc0+yc0+zc0))			// get center from gizmo note
		xc0 = NumberByKey("xc",gizmoNote,"=",":")			// center of explosion is mean center
		yc0 = NumberByKey("yc",gizmoNote,"=",":")
		zc0 = NumberByKey("zc",gizmoNote,"=",":")
	endif
	if (numtype(xc0+yc0+zc0))
		DoAlert 0, "failed to find center of mass in ExplodeProc()"
		return 1
	endif

	String grainList = getGrainNameList()				// list of grain names from the Gizmo iso-surface plot
	Variable N=ItemsInLIst(grainList)
	if (N<1)
		DoAlert 0, "no grains found in ExplodeProc()"
		return 1
	endif
//	printf "in ExplodeProc, the control '%s' sent a value of %g\r",ctrlName,varNum
//	print " "
//	print "gName = ",gName, "    explode=",explode
//	print "xc0,yc0,zc0 =",xc0,yc0,zc0,"µm"
//	print "grainList =",grainList

	Variable i,xc,yc,zc,x0,y0,z0,factor=.1			// move 0.1µm for every 1µm of Æxc
	for (i=0;i<N;i+=1)
		Wave grain =$getGrainWaveFromName(StringFromList(i,grainList))
		if (!WaveExists(grain))
			continue
		endif
		noteStr = note(grain)
		xc = NumberByKey("grain_xc",noteStr,"=")	// center of mass from wave note
		yc = NumberByKey("grain_yc",noteStr,"=")
		zc = NumberByKey("grain_zc",noteStr,"=")
		x0 = NumberByKey("x0",noteStr,"=")			// original offset from wave note
		y0 = NumberByKey("y0",noteStr,"=")
		z0 = NumberByKey("z0",noteStr,"=")

//		printf "for %s, xc=%g, yc=%g, zc=%g µm,   starting point is (%g, %g, %g)µm",NameOfWave(grain),xc,yc,zc,x0,y0,z0
//		printf "     Æ to COM is (%g, %g, %g)µm",xc-xc0,yc-yc0,zc-zc0
		x0 += explode*factor*(xc-xc0)				// new offset in exploding about (xc0,yc0,zc0)
		y0 += explode*factor*(yc-yc0)
		z0 += explode*factor*(zc-zc0)
//		printf "explode will put offset to (%g, %g,%g)\r",x0,y0,z0

		SetScale/P x x0,DimDelta(grain,0),WaveUnits(grain,0), grain	// explosion is done by changing offsets
		SetScale/P y y0,DimDelta(grain,1),WaveUnits(grain,1), grain
		SetScale/P z z0,DimDelta(grain,2),WaveUnits(grain,2), grain
	endfor
End


Function PopGrainBoundryProc(ctrlName,popNum,popStr) : PopupMenuControl
	String ctrlName
	Variable popNum
	String popStr
	String tagStr, noteStr = GetGizmoNote()
	if (stringmatch(ctrlName,"popupBndryGrain1"))
		tagStr = "grainBndry1"
	elseif (stringmatch(ctrlName,"popupBndryGrain2"))
		tagStr = "grainBndry2"
	else
		return 1
	endif
	noteStr = ReplaceStringByKey(tagStr,noteStr,popStr,"=",":")
	SetGizmoNote(noteStr)
End
Function ButtonGrainBndryProc(ctrlName) : ButtonControl
	String ctrlName

	Execute "GetGizmo gizmoName"		// check if a gizmo is up
	SVAR S_GizmoName=S_GizmoName
	if (strlen(S_GizmoName)<1)
		DoAlert 0, "No Gizmo Windows Open"
		return 1
	endif
	Variable gb1, gb2
	String gName1,gName2, str = GetGizmoNote()
	gName1 = StringByKey("grainBndry1",str,"=",":")
	gName2 = StringByKey("grainBndry2",str,"=",":")
	Wave grain=$getGrainPropertyFromObject(gName1,"isoSurface")
	gb1 = WaveExists(grain) ? NumberByKey("grainNum",note(grain),"=") : NaN
	Wave grain=$getGrainPropertyFromObject(gName2,"isoSurface")
	gb2 = WaveExists(grain) ? NumberByKey("grainNum",note(grain),"=") : NaN
	if (numtype(gb1+gb2))
		DoAlert 0, "unable to determine grains for boundary, set yet?"
		return 1
	endif

	Variable i, xc,yc,zc
	String cmd
	if (!cmpstr(ctrlName,"buttonCalcBndry"))
		Wave grainWave = $getGrainPropertyFromObject(gName1,"isoSurface")
		Wave grains = $StringByKey("grainsWave",note(grainWave),"=")	// get the wave name
		if (!WaveExists(grains))
			DoAlert 0, "unable to find original grains wave"
			return 1
		endif
		Variable radius = NumVarOrDefault("root:Packages:Grains:radiusBndry", Inf)
		String stats = BoundaryStats(gNeighbor,gb1,gb2,radius)			// xc,yc,zc is the center of the boundary
		xc=NumberByKey("xc",stats,"=")
		yc=NumberByKey("yc",stats,"=")
		zc=NumberByKey("zc",stats,"=")

		str = GetGizmoNote()												// write the center to the gizmo note
		str = ReplaceNumberByKey("xGBc",str,xc,"=",":")
		str = ReplaceNumberByKey("yGBc",str,yc,"=",":")
		str = ReplaceNumberByKey("zGBc",str,zc,"=",":")
		SetGizmoNote(str)

		Execute "ModifyGizmo startRecMacro"
		// position the line and sphere at the center
		sprintf cmd, "ModifyGizmo opName=NormalTranslate, operation=translate, data={%g,%g,%g}",xc,yc,zc
		Execute cmd
		sprintf cmd, "ModifyGizmo opName=NormalTranslateBack, operation=translate, data={%g,%g,%g}",-xc,-yc,-zc
		Execute cmd

		Variable xN,yN,zN
		str = StringByKey("bndryNormal",stats,"=")						// change the line to be in the right direction
		xN = 10*str2num(StringFromList(0,str,","))
		yN = 10*str2num(StringFromList(1,str,","))
		zN = 10*str2num(StringFromList(2,str,","))
		Wave boundaryNormal=$(GetWavesDataFolder(grains,1)+"GrainBoundaryPath")
		boundaryNormal[0][0]=-xN  ;  boundaryNormal[1][0]=xN
		boundaryNormal[0][1]=-yN  ;  boundaryNormal[1][1]=yN
		boundaryNormal[0][2]=-zN  ;  boundaryNormal[1][2]=zN
		boundaryNormal += 1
		Execute cmd
		Execute "ModifyGizmo endRecMacro"
		Execute "ModifyGizmo compile"
	elseif (!cmpstr(ctrlName,"buttonCenterBndry"))						// request to center at boundary
		str = GetGizmoNote()
		xc=NumberByKey("xGBc",str,"=",":")								// read center position from gizmo note
		yc=NumberByKey("yGBc",str,"=",":")
		zc=NumberByKey("zGBc",str,"=",":")
		if (numtype(xc+yc+zc))
			DoAlert 0, "grain boundary center not yet defined"
		endif
		sprintf cmd, "ModifyGizmo opName=translateSampleCenter, operation=translate, data={%g,%g,%g}",-xc,-yc,-zc
		Execute cmd
		Execute "ModifyGizmo compile"
	endif
	return 0
End


Function ButtonTranslateProc(ctrlName) : ButtonControl
	String ctrlName

	Execute "GetGizmo displayList"
	SVAR S_DisplayList=S_DisplayList

	Variable i, x1,x2,tx,ty,tz,  step=0.1
	i = strsearch(S_DisplayList,"opName=ortho0, operation=ortho, data={",0)
	String cmd, str = S_DisplayList[i,Inf]
	sscanf str,"opName=ortho0, operation=ortho, data={%g,%g,",x1,x2
	if (V_flag<2)
		DoAlert 0, "unable to get the otho width in ButtonTranslateProc"
		print str
		return 1
	endif
	step *= abs(x2-x1)
	i = strsearch(S_DisplayList,"opName=translateSampleCenter, operation=translate, data={",0)
	str = S_DisplayList[i,Inf]
	sscanf str," opName=translateSampleCenter, operation=translate, data={%g,%g,%g}",tx,ty,tz
	if (V_flag<3)
		DoAlert 0, "unable to get the otho width in ButtonTranslateProc"
		print str
		return 1
	endif
	if (!cmpstr(ctrlName,"buttonXminus"))
		tx -= step
	elseif (!cmpstr(ctrlName,"buttonXplus"))
		tx += step
	elseif (!cmpstr(ctrlName,"buttonYminus"))
		ty -= step
	elseif (!cmpstr(ctrlName,"buttonYplus"))
		ty += step
	elseif (!cmpstr(ctrlName,"buttonZminus"))
		tz -= step
	elseif (!cmpstr(ctrlName,"buttonZplus"))
		tz += step
	else
		KillWaves/Z TW_DisplayList
		KillStrings/Z S_DisplayList
		return 0
	endif
	sprintf cmd, "ModifyGizmo opName=translateSampleCenter, operation=translate, data={%g,%g,%g}",tx,ty,tz
	Execute cmd
	return 0
End
Function ButtonZoomProc(ctrlName) : ButtonControl
	String ctrlName

	Execute "GetGizmo displayList"
	SVAR S_DisplayList=S_DisplayList

	Variable factor=1.2
	Variable i = strsearch(S_DisplayList,"opName=ortho0, operation=ortho, data={",0)
	Variable x1,x2,y1,y2,z1,z2
	String str = S_DisplayList[i,Inf]
	sscanf str,"opName=ortho0, operation=ortho, data={%g,%g,%g,%g,%g,%g}*",x1,x2,y1,y2,z1,z2
	if (V_flag<6)
		DoAlert 0, "unable to get the zoom factor, only found "+num2istr(V_flag)+" of 6 values"
		print S_DisplayList[0,100]
		return 1
	endif
	if (cmpstr(ctrlName,"buttonZoomIn"))
		x1 *= factor
		x2 *= factor
		y1 *= factor
		y2 *= factor
	elseif (cmpstr(ctrlName,"buttonZoomOut"))
		x1 /= factor
		x2 /= factor
		y1 /= factor
		y2 /= factor
	else
		KillWaves/Z TW_DisplayList
		KillStrings/Z S_DisplayList
		return 1
	endif
	String cmd
	sprintf cmd, "ModifyGizmo opName=ortho0, operation=ortho, data={%g,%g,%g,%g,%g,%g}",x1,x2,y1,y2,z1,z2
	Execute cmd
//	sprintf cmd, "ModifyGizmo modifyObject=freeAxesCue0, property={0,0,0,%g}",(x2-x1)/6
//	Execute cmd	
	KillWaves/Z TW_DisplayList
	KillStrings/Z S_DisplayList
	return 0
End
Function ButtonHomeProc(ctrlName) : ButtonControl
	String ctrlName
	Variable or, ci,cj,ck
	String cmd, noteStr = GetGizmoNote()
	or = NumberByKey("ortho0",noteStr,"=",":")
	sprintf cmd, "ModifyGizmo opName=ortho0, operation=ortho, data={%g,%g,%g,%g,%g,%g}",-or,or,-or,or,-or,or
	Execute cmd
//	sprintf cmd, "ModifyGizmo modifyObject=freeAxesCue0, property={0,0,0,%g}",or/3
//	Execute cmd	

	ci = NumberByKey("ci",noteStr,"=",":")
	cj = NumberByKey("cj",noteStr,"=",":")
	ck = NumberByKey("ck",noteStr,"=",":")
	sprintf cmd "ModifyGizmo opName=translateSampleCenter, operation=translate, data={%g,%g,%g}",-ci,-cj,-ck
	Execute cmd
End
Function ButtonStopProc(ctrlName) : ButtonControl
	String ctrlName
	Execute "ModifyGizmo stopRotation"
End

Function PopGrainProc(ctrlName,popNum,popStr) : PopupMenuControl
	String ctrlName
	Variable popNum
	String popStr

	Execute "GetGizmo gizmoName"		// check if a gizmo is up
	SVAR S_GizmoName=S_GizmoName
	if (strlen(S_GizmoName)<1)
		DoAlert 0, "No Gizmo Windows Open"
		return 0
	endif

	// first get the rgb for this grain
	Execute "GetGizmo objectList"
	SVAR S_gizmoObjectList=S_gizmoObjectList
	String item,str
	str = "ModifyGizmo ModifyObject="+popStr+" property={ backColor"
	Variable i=0
	do
		item = StringFromList(i,S_gizmoObjectList)
		if (strlen(item)<1)
			return 1
		elseif (strsearch(item,str,0)>=0)
			break
		endif
		i += 1
	while(1)

	i = strsearch(item,"property={ backColor,",0)+21
	item = item[i,Inf]
	Variable r,g,b,alpha
	sscanf item,"%g,%g,%g,%g",r,g,b,alpha

	String cmd
	if (!cmpstr("popupGrainDim",ctrlName))		// dim the grain 'popStr'
		sprintf cmd, "ModifyGizmo ModifyObject=%s property={ backColor,%g,%g,%g,%g}",popStr,r,g,b,0.1
		Execute cmd
	elseif (!cmpstr("popupGrainBrite",ctrlName))	// brighten the grain 'popStr'
		sprintf cmd, "ModifyGizmo ModifyObject=%s property={ backColor,%g,%g,%g,%g}",popStr,r,g,b,0.9
		Execute cmd
	elseif (!cmpstr("popupGrainRestore",ctrlName))	// restore the grain 'popStr'
		// find the correct alpha
		String wName = getGrainWaveFromName(popStr)
		Wave wav = $wName
		alpha = NumberByKey("alpha",note(wav),"=")
		alpha = numtype(alpha) ? 0.6 : alpha
		sprintf cmd, "ModifyGizmo ModifyObject=%s property={ backColor,%g,%g,%g,%g}",popStr,r,g,b,alpha
		Execute cmd
	endif
	Killwaves/Z TW_gizmoObjectList
	KillStrings/Z S_gizmoObjectList
End
Function ButtonGrainColorProc(ctrlName) : ButtonControl
	String ctrlName

	String grainList = getGrainNameList()				// list of grain names from the Gizmo iso-surface plot
	if (strlen(grainList)<1)
		DoAlert 0, "no grains found in ButtonGrainColorProc"
		return 1
	endif
	Variable alpha
	if (!cmpstr(ctrlName,"buttonAllDim"))
		alpha = 0.05
	elseif (!cmpstr(ctrlName,"buttonAllSolid"))
		alpha = 1
	elseif (!cmpstr(ctrlName,"buttonAllRestore"))
		alpha = -1										// flag to use the original grains alpha
	else
		return 1
	endif
	// now apply alpha to all of the grains
	String grain, cmd, wName
	Variable r,g,b,a
	Variable i,N=ItemsInList(grainList)
	Execute "ModifyGizmo startRecMacro"
	for (i=0;i<N;i+=1)
		grain =StringFromList(i,grainList)
		// first get the color for the grain
		sscanf getGrainPropertyFromObject(grain,"backColor"), "%g,%g,%g,",r,g,b
		if (V_flag!=3)
			continue
		endif
		a = alpha
		if (alpha<0)										// get original alpha
			wName = getGrainWaveFromName(grain)		// returns the wave name for object objName
			Wave wav = $wName
			a = NumberByKey("alpha",note(wav),"=")
			a = numtype(a) ? 0.6 : a
		endif
		sprintf cmd, "ModifyGizmo ModifyObject=%s property={backColor,%g,%g,%g,%g}",grain,r,g,b,a
		Execute cmd
	endfor
	Execute "ModifyGizmo endRecMacro"
	return 0
End


Function/S Get_GrainNo_Dialog(grainNo)
	Wave grainNo
	if (!WaveExists(grainNo))
		String PromptList = WaveList("*",";","DIMS:3")
		if (cmpstr(GetDataFolder(1),"root:raw:"))
			String fldrSav=GetDataFolder(1)
			SetDataFolder root:raw:
			String tempWlist = WaveList("*",";","DIMS:3")
			SetDataFolder $fldrSav

			Variable i = 0
			for (i=0;i<ItemsInList(tempWlist);i+=1)
				PromptList = AddListItem("root:raw:"+StringFromList(i,tempWlist),PromptList,";",Inf)
			endfor
		endif
		String grainNoWave = SelectString(WhichListItem("grainNo",PromptList)>=0,"root:raw:grainNo","grainNo")
		Prompt grainNoWave, "3D wave with all grain numbers", popup,PromptList
		DoPrompt "pick grainNo",grainNoWave
		if (V_flag)
			grainNoWave = ""
		endif
		Wave grainNo = $grainNoWave
	endif
	return GetWavesDataFolder(grainNo,2)
End

Function/S getGrainNameList()			// returns the list of iso-surface names from a Gizmo iso-surface plot
	Execute "GetGizmo gizmoName"		// check if a gizmo is up
	SVAR S_GizmoName=S_GizmoName
	if (strlen(S_GizmoName)<1)
		DoAlert 0, "No Gizmo Windows Open"
		return ""
	endif
	Execute "GetGizmo objectList"
	SVAR S_gizmoObjectList=S_gizmoObjectList
	String grainList="", item
	Variable i=0
	do
		item = StringFromList(i,S_gizmoObjectList)
		if (strlen(item)<1)
			break
		elseif (strsearch(item,"AppendToGizmo isoSurface=",0)>=0)
			grainList += StringByKey("name",item,"=",",")+";"
		endif
		i += 1
	while(1)
	Killwaves/Z TW_gizmoObjectList
	KillStrings/Z S_gizmoObjectList
	return grainList
End
Function/S getGrainNumbersList()			// returns the list of waves making the iso-surface from a Gizmo iso-surface plot
	String grainList = ""
	String wName, wList = getGrainWaveList()	
	Variable i,j
	for(i=0,wName="x";strlen(wName);i+=1)
		wName = StringFromList(i,wList)
		if (exists(wName)!=1)
			continue
		endif
		j=NumberByKey("grainNum",note($wName),"=")
		if (!numtype(j))
			grainList = AddListItem(num2istr(j),grainList,";",Inf)
		endif
	endfor
	return grainList
End
Function/S getGrainWaveList()			// returns the list of waves making the iso-surface from a Gizmo iso-surface plot
	Execute "GetGizmo gizmoName"		// check if a gizmo is up
	SVAR S_GizmoName=S_GizmoName
	if (strlen(S_GizmoName)<1)
		DoAlert 0, "No Gizmo Windows Open"
		return ""
	endif
	Execute "GetGizmo objectList"
	SVAR S_gizmoObjectList=S_gizmoObjectList
	String grainList="", item
	Variable i=0
	do
		item = StringFromList(i,S_gizmoObjectList)
		if (strsearch(item,"AppendToGizmo isoSurface=",0)>=0)
			item = RemoveFromList("\tAppendToGizmo",item," ")
			item = StringByKey("isoSurface",item,"=",",")		// the file name
			if (strlen(item)>=1)
				grainList += item+";"
			endif
		endif
		i += 1
	while(strlen(item)>=1)
	Killwaves/Z TW_gizmoObjectList
	KillStrings/Z S_gizmoObjectList
	return grainList
End
Function/S getGrainWaveFromName(objName)		// returns the wave name for object objName
	String objName									// name of object to get wave name for
	Execute "GetGizmo gizmoName"		// check if a gizmo is up
	SVAR S_GizmoName=S_GizmoName
	if (strlen(S_GizmoName)<1)
		DoAlert 0, "No Gizmo Windows Open"
		return ""
	endif
	Execute "GetGizmo objectList"
	SVAR S_gizmoObjectList=S_gizmoObjectList
	String wName="", item
	Variable i=0
	do
		item = StringFromList(i,S_gizmoObjectList)+","
		if (strlen(item)<2)
			break
		elseif (strsearch(item,"AppendToGizmo isoSurface=",0)<0)
			i += 1
			continue
		elseif (strsearch(item,"name="+objName+",",0)<0)
			i += 1
			continue
		endif
		// found it
		i = strsearch(item,"isoSurface=",0)
		wName = StringByKey("isoSurface",  item[i,Inf],"=",",")
		break
	while(1)
	Killwaves/Z TW_gizmoObjectList
	KillStrings/Z S_gizmoObjectList
	return wName
End

Function/S getGrainPropertyFromObject(objName,property) // returns the value of property for object objName
	String objName							// name of object to get property of
	String property							// property to get
	Execute "GetGizmo gizmoName"			// check if a gizmo is up
	SVAR S_GizmoName=S_GizmoName
	if (strlen(S_GizmoName)<1)
		return ""
	endif
	Execute "GetGizmo objectList"
	SVAR S_gizmoObjectList=S_gizmoObjectList
	String search, value = ""
	sprintf search, "ModifyGizmo ModifyObject=%s property={ %s,",objName,property
	Variable i = strsearch(S_gizmoObjectList,search,0)
	if (i>=0)								// found it
		i += strlen(search)
		value = S_gizmoObjectList[i,Inf]
		i = strsearch(value,"}",0)
		value = value[0,i-1]
	endif
	if (strlen(value)==0)
		String item
		i=0
		do
			item = StringFromList(i,S_gizmoObjectList)
			if (strsearch(item,"\tAppendToGizmo ",0)>=0)
				item = RemoveFromList("\tAppendToGizmo",item," ")
				if (stringmatch(StringByKey("name",item,"=",","),objName)) // check for matching name
					value =StringByKey(property,item,"=",",")	
					break
				endif
			endif
			i += 1
		while(strlen(item)>=1)
	endif
	Killwaves/Z TW_gizmoObjectList
	KillStrings/Z S_gizmoObjectList
	return value
End



//	=========================================================================
//	=========================================================================
//	=========================================================================
//								Start of Grain Boundary Analysis
//	=========================================================================

// This routine goes through each grain and checks each of its neighbors, if a pair matches the criterion
// for a sigma N boundary then it prints it out.  To match a sigma N, the axis has to be aligned to within
// angleTol, and the total rotation also has to be within angleTol.  Note, this routine skips grains that are
// smaller than minVol.
Function checkAllSigmaBoundarys(grains,minVol,minBndryPoints,angleTol)
	Wave grains						// 3d wave with all of the grain numbers (e.g. gNeighbor)
	Variable minVol						// ignore grains smaller than this
	Variable minBndryPoints			// mininum number of points in common needed for inclusion of a neighbor
	Variable angleTol					// angle tolerance (degrees)
	if (!WaveExists(grains) || WaveDims(grains)!=3 || !(minVol>=0)  || !(angleTol>0) || !(minBndryPoints>=0))
		minVol = minVol>0 ? minVol : 20
		angleTol = angleTol>0 ? angleTol : 15
		minBndryPoints = !(minBndryPoints>=0) ? 1 : minBndryPoints
		String grainName=NameOfWave(grains)
		Prompt grainName,"3d wave with all grains",popup,ListMultiGrainMats()
		Prompt minVol, "min volume for a grain"
		Prompt angleTol, "angular tolerance for rotation and axis alignment (degreee), for Brandon Criterion use 15"
		Prompt minBndryPoints, "only include neighbors with at least this many bordering voxels"
		DoPrompt "grain array",grainName,minVol,minBndryPoints,angleTol
		if (V_flag)
			return 1
		endif
		Wave grains = $grainName
		if (ItemsInList(GetRTStackInfo(0))<2)
			printf "		checkAllSigmaBoundarys(%s,%g,%g,%g)\r",GetWavesDataFolder(grains,2),minVol,minBndryPoints,angleTol
		endif
	endif
	if (!WaveExists(grains) || WaveDims(grains)!=3 || !(angleTol>=0)  || !(angleTol>0) || !(minBndryPoints>=0))
		return 1
	endif
	String gName = GetWavesDataFolder(grains,2)

	Wave gVol=$(gName+"_Vol")	// number of voxels in the grain
	Variable N = numpnts(gVol)
	Variable totalAngle				// total rotation angle between grains (degrees)
	Variable j,i,xx,yy,zz,Sigma
	Variable angleErr,axisErr		// error in angle and axis alignment for a particular SigmaN
	String range,list
	String wName
	wName = UniqueName("sigmaAxis",1,0) ; Make/N=3/D $wName ; Wave sigmaAxis=$wName
	String degree = StrVarOrDefault("root:Packages:Grains:degree","¡")
	printf "Sigma\t\t i\tj\t\t  angle(%s)\t\t\t\t\t(hkl)\t\t\t\t\t\t Æaxis%s\t\t\t Æangle%s\r",degree,degree,degree
	for (i=0;i<N;i+=1)
		if (gVol[i]<minVol)
			break
		endif
		range = ReplaceString(";",makeBorderList(grains,i),",")
		for(j = NextInRange(range,i);!numtype(j)&&gVol[j]>=minVol;j = NextInRange(range,j))
			if (j==i)
				continue
			endif
			list = BoundaryStats(gNeighbor,j,i,Inf)
			if (NumberByKey("Nboundary1",list,"=")<minBndryPoints && NumberByKey("Nboundary2",list,"=")<minBndryPoints)
				continue
			endif
			sscanf StringByKey("hklAxis",list,"="), "%g,%g,%g",xx,yy,zz
			sigmaAxis = {xx,yy,zz}
			totalAngle = NumberByKey("totalAngle",list,"=")

			Sigma = -1
			do						// this loop checks for all SigmaN that match totalAngle and sigmaAxis
				sigmaType(totalAngle,sigmaAxis,angleTol,angleErr,axisErr,Sigma)
				if (Sigma>0)
					integerDirectionExact(sigmaAxis,12)
					sigmaAxis = abs(sigmaAxis)
					Sort sigmaAxis, sigmaAxis
					printf "%g   \t\t%d\t%d\t\t%8.4f\t\t(%8.5f, %8.5f, %8.5f)\t\t%7.3f\t\t\t%+7.3f\r"Sigma,i,j,totalAngle,sigmaAxis[0],sigmaAxis[1],sigmaAxis[2],axisErr, angleErr
				endif
			while(Sigma>0)

		endfor
	endfor
	KillWaves/Z sigmaAxis
End
//Function checkAllSigmaBoundarys(grains,minVol,angleTol)
//	Wave grains						// 3d wave with all of the grain numbers (e.g. gNeighbor)
//	Variable minVol						// ignore grains smaller than this
//	Variable angleTol					// angle tolerance (degrees)
//	if (!WaveExists(grains) || WaveDims(grains)!=3 || !(angleTol>=0)  || !(angleTol>0))
//		minVol = minVol>0 ? minVol : 20
//		angleTol = angleTol>0 ? angleTol : 1
//		String grainName=NameOfWave(grains)
//		Prompt grainName,"3d wave with all grains",popup,ListMultiGrainMats()
//		Prompt minVol, "min volume for a grain"
//		Prompt angleTol, "angular tolerance for rotation and axis alignment (degreee)"
//		DoPrompt "grain array",grainName,minVol,angleTol
//		if (V_flag)
//			return 1
//		endif
//		Wave grains = $grainName
//		if (ItemsInList(GetRTStackInfo(0))<2)
//			printf "		checkAllSigmaBoundarys(%s,%g,%g)\r",GetWavesDataFolder(grains,2),minVol,angleTol
//		endif
//	endif
//	if (!WaveExists(grains) || WaveDims(grains)!=3 || !(angleTol>=0)  || !(angleTol>0))
//		return 1
//	endif
//	String gName = GetWavesDataFolder(grains,2)
//
//	Wave gVol=$(gName+"_Vol")	// number of voxels in the grain
//	Variable N = numpnts(gVol)
//	Variable totalAngle				// total rotation angle between grains (degrees)
//	Variable j,i,xx,yy,zz,sig
//	String range,list
//	String wName
//	wName = UniqueName("sigmaAxis",1,0) ; Make/N=3/D $wName ; Wave sigmaAxis=$wName
//	String degree = StrVarOrDefault("root:Packages:Grains:degree","¡")
//	printf "sigma\t i\tj\t\t  angle(%s)\t\t\t(hkl)\r",degree
//	for (i=0;i<N;i+=1)
//		if (gVol[i]<minVol)
//			break
//		endif
//		range = ReplaceString(";",makeBorderList(grains,i),",")
//		for(j = NextInRange(range,i);!numtype(j)&&gVol[j]>=minVol;j = NextInRange(range,j))
//			if (j==i)
//				continue
//			endif
//			list = BoundaryStats(gNeighbor,j,i,Inf)
//			sscanf StringByKey("hklAxis",list,"="), "%g,%g,%g",xx,yy,zz
//			sigmaAxis = {xx,yy,zz}
//			totalAngle = NumberByKey("totalAngle",list,"=")
//			sig = sigmaType(totalAngle,sigmaAxis,angleTol)
//
//			integerDirectionExact(sigmaAxis,12)
//			sigmaAxis = abs(sigmaAxis)
//			Sort sigmaAxis, sigmaAxis
//			if (sig>0)
//				printf "%d   \t%d\t%d\t\t%8.4f\t\t(%8.5f, %8.5f, %8.5f)\r"sig,i,j,totalAngle,sigmaAxis[0],sigmaAxis[1],sigmaAxis[2]
//			endif
//		endfor
//	endfor
//	KillWaves/Z sigmaAxis
//End



// This routine is only called by checkAllSigmaBoundarys().  It looks at the given axis and angle and 
// returns the sigma N type.  If it is not a match to any type, it returns 0.
Static Function sigmaType(angle,axisIn,angleTol,angleErr,axisErr,Sigma)
	Variable angle						// measured total rotation angle (degrees)
	Wave axisIn							// measured axis of rotation given
	Variable angleTol					// angle tolerance (degrees), actual tolerance is angleTol/sqrt(SigmaN)
										// for Brandon Criterion, angleTol=15¡
	Variable &angleErr, &axisErr		// actual error in angle and axis for this Sigma
	Variable &Sigma					// on entry, it is the last Sigma, on exit the result
										// first time use a value of Sigma²0, for subsequent calls, Sigma=lastSigma
	Variable lastSigma=Sigma+1e-4
	Sigma = 0
	angle = abs(angle)
	if (angle<angleTol && lastSigma<1)	// the Sigma 1 is sort of simple, no need for axis
		angleErr = angle
		axisErr = NaN
		Sigma = 1
		return 1
	endif
	String wName
	wName = UniqueName("hkl",1,0) ; Make/N=3/D $wName ; Wave hkl=$wName
	wName = UniqueName("axis",1,0) ; Make/N=3/D $wName ; Wave axis=$wName
	axis = abs(axisIn)					// make a "normal form" of axis, all positive and ordered
	Sort axis, axis
	normalize(axis)					// and normalized too

	Wave sigmaCSL = root:Packages:Grains:sigmaCSL
	if (!WaveExists(sigmaCSL))
		initGrainGizmo()				// this will make sigmaCSL
	endif

	Variable Brandon					// Brandon Criterion for angular tolerance, usually 15/sqrt(sigmaN)  (degree)
	// for use of Brandon's Criterion of a CSL, see:  http://www.amc.anl.gov/ANLSoftwareLibrary/02-MMSLib/DIFF/FINDCSL/FindCSL2ang2.for
	// or Hanada, et al, Acta Metall, 34(1986) No.1, pp 13-21.
	Variable i, delta
	for (i=0;i<DimSize(sigmaCSL,0);i+=1)
		if (lastSigma > sigmaCSL[i][0])	// skip SigmaN that were already checked
			continue
		endif
		hkl = sigmaCSL[i][2+p]
		normalize(hkl)
		Brandon = angleTol/sqrt(sigmaCSL[i][0])
		angleErr = (angle - sigmaCSL[i][1])
		axisErr = acos(max(-1,min(MatrixDot(axis,hkl),1))) * 180/PI
		delta = 2*atan(  sqrt(tan(angleErr/2*PI/180)^2 + tan(axisErr/2*PI/180)^2)  )* 180/PI
		if (delta < Brandon)
			Sigma = sigmaCSL[i][0]
			Sigma = round(Sigma*100)/100
			break
		endif
	endfor
	KillWaves/Z hkl,axis
	return Sigma
End
//Static Function sigmaType(angle,axisIn,angleTol)
//	Variable angle						// total rotation angle (degrees)
//	Wave axisIn							// axis of rotation given
//	Variable angleTol					// angle tolerance (degrees)
//
//	angle = abs(angle)
//	if (angle<angleTol)					// the Sigma 1 is sort of simple
//		return 1
//	endif
//	Variable cosMin = cos(angleTol*PI/180)	// only accept is cos > cosMin
//	String wName
//	wName = UniqueName("hkl",1,0) ; Make/N=3/D $wName ; Wave hkl=$wName
//	wName = UniqueName("axis",1,0) ; Make/N=3/D $wName ; Wave axis=$wName
//	axis = abs(axisIn)					// make a "normal form" of axis, all positive and ordered
//	Sort axis, axis
//	normalize(axis)					// and normalized too
//
//	Wave sigmaCSL = root:Packages:Grains:sigmaCSL
//	if (!WaveExists(sigmaCSL))
//		initGrainGizmo()				// this will make sigmaCSL
//	endif
//
//	Variable i, sigma = 0
//	for (i=0;i<DimSize(sigmaCSL,0);i+=1)
//		hkl = sigmaCSL[i][2+p]
//		normalize(hkl)
//		if (MatrixDot(axis,hkl)>cosMin && abs(angle-sigmaCSL[i][1])<angleTol)
//			sigma = sigmaCSL[i][0]
//			break
//		endif
//	endfor
//	KillWaves/Z hkl,axis
//	return sigma
//End


Function/S BoundaryStats(grains,i1,i2,radius)// calc all information about a boundary bwteen two neighboring grains
	Wave grains						// 3d wave with all of the grain numbers (e.g. gNeighbor)
	Variable i1,i2						// grain numbers
	Variable radius						// only consider points within radius of cursor

	String str
	String micron = StrVarOrDefault("root:Packages:Grains:micron","µm")
	if (!WaveExists(grains) || numtype(i1+i2) || i1<0 || i2<0 || numtype(radius)==2 || radius<=0)
		String grainName=NameOfWave(grains)
		Prompt grainName,"3d wave with all grains",popup,ListMultiGrainMats()
		String gizmoNote=GetGizmoNote()
		if (numtype(i1))
			str = StringByKey("grainBndry1",gizmoNote,"=",":")
			i1 = str2num(str[5,Inf])
			i1 = numtype(i1) ? 0 : max(i1,0)
		endif
		if (numtype(i2))
			str = StringByKey("grainBndry2",gizmoNote,"=",":")
			i2 = str2num(str[5,Inf])
			i2 = numtype(i2) ? 0 : max(i2,0)
		endif
		Prompt i1, "first grain number"
		Prompt i2, "second grain number"
		radius = (numtype(radius)==2 || radius<0) ? NumVarOrDefault("root:Packages:Grains:radiusBndry", Inf) : radius
		Prompt radius, "distance from cursor to accept when finding grain boundary points ("+micron+")"
		DoPrompt "choose grains",grainName,i1,i2,radius
		if (V_flag)
			return ""
		endif
		Wave grains = $grainName
		printf "   	BoundaryStats(%s,%g,%g,%g)\r",NameOfWave(grains),i1,i2,radius
	endif
	if (!WaveExists(grains))
		return ""
	endif
	WaveStats/Q grains
	Variable maxGrain = V_max
	if (numtype(i1+i2) || i1<0 || i2<0 || i1>maxGrain || i2>maxGrain || i1==i2 || numtype(radius)==2 || radius<=0)
		if (ItemsInList(GetRTStackInfo(0))<2)
			DoAlert 0,"invalid inputs"
		endif
		return ""
	endif
	// check that grains are neighbors
	String list = makeBorderList(grains,i2)
	if (WhichListItem(num2istr(i1),list)<0)
		if (ItemsInList(GetRTStackInfo(0))<2)
			sprintf str, "grains %d and %d are not neighbors\r",i1,i2
			printf str
			DoAlert 0, str
		endif
		return ""
	endif
	// all checks are done, now do the work

	String fName = GetWavesDataFolder(grains,1)
	Wave indicies = $(fname+"raw:indicies")
	Wave XYZindex = $(fname+"raw:XYZindex")
	Wave OrientMat = $(fname+"raw:OrientMat")
	String gName = GetWavesDataFolder(grains,2)
	Wave COMsingle = $(gName+"_COM")

	Variable Nx,Ny,Nz
	Nx=DimSize(XYZindex,0)
	Ny=DimSize(XYZindex,1)
	Nz=DimSize(XYZindex,2)
	String out=""
	out=ReplaceNumberByKey("i1",out,i1,"=")
	out=ReplaceNumberByKey("i2",out,i2,"=")

	Wave g1 = $getGrainCoords(grains,i1)		// list of grain coords for each grain
	Wave g2 = $getGrainCoords(grains,i2)
	Variable N1=DimSize(g1,0), N2=DimSize(g2,0)
	Variable single1 = grainCOMsingle(g1)		// single index to a measured point nearest to COM
	Variable single2 = grainCOMsingle(g2)

	Variable totalAngle = angleBetweenPointsSymReduced(single1,single2)
	out=ReplaceNumberByKey("totalAngle",out,totalAngle,"=")

	Make/N=(N2)/T/O g21neighgors			// contains a list of neighbors to each grain in g2 (neighbors are in g2)
	g21neighgors = ""
	Variable i,j, ix0,iy0,iz0,  ix,iy,iz
	Variable Nbndry1=0, Nbndry2=0, Nadd

	str=UniqueName("offsets",1,0)
	// offsets[][] contains the 6 voxels to the 1st nearest neighbors (faces in common)
	Make/N=(6,3) $str
	Wave offsets=$str
	offsets[0][0]= {-1,1,  0,0,  0,0}
	offsets[0][1]= {  0,0,-1,1,  0,0}
	offsets[0][2]= {  0,0,  0,0,-1,1}
	Variable Noffset = DimSize(offsets,0)

	for (i=0;i<N2;i+=1)
		Nadd = 0
		for (j=0;j<Noffset;j+=1)				// loop through the offsets
			ix = g2[i][0] + offsets[j][0]
			iy = g2[i][1] + offsets[j][1]
			iz = g2[i][2] + offsets[j][2]
			if (numtype(ix+iy+iz) || ix>=Nx || iy>=Ny || iz>=Nz || ix<0 || iy<0 || iz<0)
				continue						// point is NaN or out of range
			elseif (grains[ix][iy][iz]==i1)
				sprintf str "%d,%d,%d",ix,iy,iz
				g21neighgors[i] = AddListItem(str, g21neighgors[i])
				Nadd = 1
			endif
		endfor
		Nbndry2 += Nadd
	endfor

	Variable jlast=-1							// remove all non-surface points from g2 and compress
	for (i=0;i<Nbndry2;i+=1)
		for (j=jlast+1;j<N2;j+=1)
			if (strlen(g21neighgors[j])>0)
				g2[i][0,2] = g2[j][q]
				jlast = j
				break
			endif
		endfor
	endfor
	Redimension/N=(Nbndry2,3) g2

	Make/N=(N1)/O g1_flag_temp__
	g1_flag_temp__ = 0
	Variable k
	for (i=0;i<N2;i+=1)
		if (strlen(g21neighgors[i])>0)
			str = g21neighgors[i]
			for (j=0;j<ItemsInList(str);j+=1)	// check each triplet in g21neighgors[i]
				ix = str2num(StringFromList(0,str,","))
				iy = str2num(StringFromList(1,str,","))
				iz = str2num(StringFromList(2,str,","))
				for (k=0;k<N1;k+=1)
					if (g1[k][0]==ix && g1[k][1]==iy && g1[k][2]==iz)
						g1_flag_temp__[k] = 1
						break
					endif
				endfor
			endfor
		endif
	endfor
	Nbndry1 = sum(g1_flag_temp__)

	jlast=-1									// remove all non-surface points from g1 and compress
	for (i=0;i<Nbndry1;i+=1)
		for (j=jlast+1;j<N1;j+=1)
			if (g1_flag_temp__[j])
				g1[i][0,2] = g1[j][q]
				jlast = j
				break
			endif
		endfor
	endfor
	Redimension/N=(Nbndry1,3) g1
	KillWaves g1_flag_temp__
	out=ReplaceNumberByKey("Nboundary1",out,Nbndry1,"=")
	out=ReplaceNumberByKey("Nboundary2",out,Nbndry2,"=")
	//
	// g1 and g2 now contain only points in the interface (i.e. points that border the other grain by a face)

	// convert g1 and g2 from indicies to microns
	Variable x0 = DimOffset(grains,0), 	y0 = DimOffset(grains,1), z0 = DimOffset(grains,2)
	Variable dx = DimDelta(grains,0), dy = DimDelta(grains,1), dz = DimDelta(grains,2)
	g1[][0] = g1[p][0]*dx + x0
	g1[][1] = g1[p][1]*dy + y0
	g1[][2] = g1[p][2]*dz + z0
	g2[][0] = g2[p][0]*dx + x0
	g2[][1] = g2[p][1]*dy + y0
	g2[][2] = g2[p][2]*dz + z0

	// find the center of mass of the grain boundary (xc,yc,zc)
	Variable xc=0,yc=0,zc=0

	// here is the rule, if the cursor is on a boundary voxel, then use the cursor, otherwise use com.
	str = getGrainPropertyFromObject("sphereCursor","sphere")
	if (str2num(str[1,Inf])>1e-5)				// cursor is active (since radius is >0)
		xc = NumVarOrDefault("root:Packages:Grains:Xcursor",0)
		yc = NumVarOrDefault("root:Packages:Grains:Ycursor",0)
		zc = NumVarOrDefault("root:Packages:Grains:Zcursor",0)
	else
		if (Nbndry1>=3)
			Make/N=(Nbndry1)/O temp_sum_up_
			temp_sum_up_ = g1[p][0]  ;  xc = faverage(temp_sum_up_)
			temp_sum_up_ = g1[p][1]  ;  yc = faverage(temp_sum_up_)
			temp_sum_up_ = g1[p][2]  ;  zc = faverage(temp_sum_up_)
		endif
		if (Nbndry2>=3)
			Make/N=(Nbndry2)/O temp_sum_up_
			temp_sum_up_ = g2[p][0]  ;  xc = (xc+faverage(temp_sum_up_))/2
			temp_sum_up_ = g2[p][1]  ;  yc = (yc+faverage(temp_sum_up_))/2
			temp_sum_up_ = g2[p][2]  ;  zc = (zc+faverage(temp_sum_up_))/2
		endif
	endif
	out=ReplaceNumberByKey("xc",out,xc,"=")
	out=ReplaceNumberByKey("yc",out,yc,"=")
	out=ReplaceNumberByKey("zc",out,zc,"=")

	// refilter g1 and g2 removing all points more than  distance radius from (xc,yc,zc)
	Variable rad2 = radius*radius
	for (i=0;i<Nbndry1;i+=1)
		if (((xc-g1[i][0])^2 + (yc-g1[i][1])^2 + (zc-g1[i][2])^2)>rad2)
			g1[i][] = NaN									// if a point is farther than radius from (xc,yc,zc) set to NaN
		endif
	endfor
	DeleteNaNs(g1)											// remove the NaNs just added, and compress
	for (i=0;i<Nbndry2;i+=1)
		if (((xc-g2[i][0])^2 + (yc-g2[i][1])^2 + (zc-g2[i][2])^2)>rad2)
			g2[i][] = NaN									// if a point is farther than radius from (xc,yc,zc) set to NaN
		endif
	endfor
	DeleteNaNs(g2)											// remove the NaNs just added, and compress

	Make/N=3/D/O normal1,normal2, normal12			// vector in the average surface normal for each grain (beam line coords)
	NormalOfPlane(g1,normal1)							// normal1 is in beam line coords
	normalize(normal1)
	NormalOfPlane(g2,normal2)							// normal2 is in beam line coords
	normalize(normal2)
	Duplicate/O g1 bndryScatterPoints

	Variable normalAngleErr = acos(max(min(MatrixDot(normal1,normal2),1),-1))
	if (!numtype(normalAngleErr))
		out=ReplaceNumberByKey("normalAngleErr",out,normalAngleErr,"=")
	endif
	normal12 = (normal1+normal2)/2						// average normal to boundary
	String degree = StrVarOrDefault("root:Packages:Grains:degree","¡")
	if (normalAngleErr*180/PI>5 && ItemsInList(GetRTStackInfo(0))<2)
//	if (normalAngleErr*180/PI>5 && !stringmatch(StringFromList(0, GetRTStackInfo(0)),"GizmoAroundOneGrain"))
		DoAlert 0, "WARNING in BoundaryStats(), angle between normals > 5"+degree+", consider reducing radius"
		normal12 = normal1
	endif
	if (!numtype(normal12[0]))
		sprintf str, "%g,%g,%g",normal12[0],normal12[1],normal12[2]
		out=ReplaceStringByKey("bndryNormal",out,str,"=")
	endif
	Variable maxhkl=12
	Make/N=(3,3)/O mat1,mat2
	i = COMsingle[i1]
	mat1 = OrientMat[i][p][q]
	i = COMsingle[i2]
	mat2 = OrientMat[i][p][q]
	Make/N=3/O/D axisOfMat							// get the axis of the rotation
	Make/N=(3,3)/O/D BASi							// BASi, total rotation matrix from A to B, cubic reduced
	symReducedRotation(mat1,mat2,BASi)				// BASi = (B Inv(A) S^i) = (mat2 Inv(mat1) S^i),   A=mat1, B=mat2

	axisOfMatrix(BASi,axisOfMat,squareUp=1)				// make the axis of the rotation matrix BASi (grain2 coordinates)
	Make/N=3/O/D hklAxis								// hkl of total rotation matrix axis, in terms of grains1
	hklAxis = -axisOfMat								// change from grain2 to grain1 (for the axis they are just opposite)
	integerDirectionExact(hklAxis,maxhkl)				// hkl of rotation axis in terms of grain1
	MatrixOp/O tempAxisOfMat = Inv(mat2) x axisOfMat// transform from grain2(hkl) coords to sample (xyz)
	axisOfMat = tempAxisOfMat
	TransformSample2Beamline(axisOfMat)				// transform axisOfMat from sample (xyz) -> beamline (XYZ)
	sprintf str, "%g,%g,%g",hklAxis[0],hklAxis[1],hklAxis[2]
	out=ReplaceStringByKey("hklAxis",out,str,"=")
	sprintf str, "%g,%g,%g",axisOfMat[0],axisOfMat[1],axisOfMat[2]
	out=ReplaceStringByKey("axisOfMat",out,str,"=")
//	axisOfMatrix(BASi,hklAxis,squareUp=1)				// make the axis of the rotation matrix BASi (sample coords, xyz)
///	hklAxis *= -1										// change from grain2 to grain1
//	integerDirectionApprox(hklAxis,maxhkl)			// hkl of rotation axis in terms of grain1
//	axisOfMatrix(BASi,axisOfMat,squareUp=1)			// make the axis of the rotation matrix BASi (sample coords, xyz)
//	TransformSample2Beamline(axisOfMat)				// transform axisOfMat from sample (xyz) -> beamline (XYZ)
	Variable twist=naN,tilt=NaN, angleBetweenAxes=NaN// twist and tilt angles for the boundary (degree)
	if (normalize(normal12)>0.1)						// if normal12 is valid, then proceed to calculate the hkl1 & hkl2
		// next find the hkl for each of the two grains that lies along normal12
		Make/N=3/O/D normalXHF
		// XYZ -> xyz,  	convert normal from beamline coords (XYZ) to the sample (xyz) == (X -H -F)
		normalXHF = normal12
		TransformBeamline2Sample(normalXHF)
		MatrixOp/O hkl1 = mat1 x normalXHF			//  A*(xyz) = (hkl)
		integerDirectionApprox(hkl1,maxhkl)
		MatrixOp/O hkl2 = mat2 x normalXHF
		integerDirectionApprox(hkl2,maxhkl)
		Variable cosine = abs(MatrixDot(axisOfMat,normal12))
		angleBetweenAxes = acos(min(cosine,1))*180/PI
		TiltTwistOfBoundary(mat1,mat2,normalXHF,twist,tilt)	// find the twist and tilt of the boundary (degree))
		sprintf str, "%g,%g,%g",hkl1[0],hkl1[1],hkl1[2]
		out=ReplaceStringByKey("hkl1",out,str,"=")
		sprintf str, "%g,%g,%g",hkl2[0],hkl2[1],hkl2[2]
		out=ReplaceStringByKey("hkl2",out,str,"=")
	endif

	if (ItemsInList(GetRTStackInfo(0))<2)
		printf "   for the boundary between grains %g and %g in '%s':\r",i1,i2,NameOfWave(grains)
		printf "   angle between grains %d and %d is %g%s,  rotated about hkl=(%g, %g, %g) of grain1,",i1,i2,totalAngle,degree,hklAxis[0],hklAxis[1],hklAxis[2]
		printf "   --> (beamLine coords)=(%g, %g, %g)\r", axisOfMat[0],axisOfMat[1],axisOfMat[2]
		if (!numtype(angleBetweenAxes))
			printf "   angle between axis of rotation and interface normal is %g%s (after cubic symmetry operations),",angleBetweenAxes,degree
			printf "   tilt=%g%s,  twist=%g%s",tilt,degree,twist,degree
			Variable tiltTwistPure = 5							// threshold to decide whether to flag boundary as twist or tilt
			if (abs(angleBetweenAxes-90)<tiltTwistPure)
				printf ",   close to a pure tilt boundary"
			elseif (angleBetweenAxes<tiltTwistPure || abs(angleBetweenAxes-180)<tiltTwistPure)
				printf ",   close to a pure twist boundary"
			endif
			printf "\r"
		endif
		printf "   grain 1 has %d points in the surface, and grain 2 has %d in the surface\r",Nbndry1,Nbndry2
		if (!numtype(normal1[0]) || !numtype(normal2[0]))
			printf "   surface normal (average) between grains 1 & 2 is  {%g, %g, %g}   angle between individual normals is %.3g deg\r",normal12[0],normal12[1],normal12[2],normalAngleErr*180/PI
			printf "         surface normals are: grain1-> {%g, %g, %g}   and   grain2-> {%g, %g, %g}\r",normal1[0],normal1[1],normal1[2],normal2[0],normal2[1],normal2[2]
			printf "   normals use voxels out to %g %s away from point (%g, %g, %g)%s\r",radius,micron,xc,yc,zc,micron
			printf "   approx surface normals are grain1->(%g %g %g)   and   grain2->(%g %g %g) (with an hkl limit of %g)\r",hkl1[0],hkl1[1],hkl1[2],hkl2[0],hkl2[1],hkl2[2],maxhkl
			hkl1 = abs(hkl1) ; Sort hkl1, hkl1
			hkl2 = abs(hkl2) ; Sort hkl2, hkl2
			printf "   							or grain1->(%g %g %g)   and   grain2->(%g %g %g)\r",hkl1[0],hkl1[1],hkl1[2],hkl2[0],hkl2[1],hkl2[2]
		endif
	endif
	KillWaves/Z hkl1,hkl2,mat1,mat2, axisOfMat,tempAxisOfMat, rotAB, normalXHF,BASi,hklAxis
	KillWaves/Z W_coef,W_sigma,W_ParamConfidenceInterval,temp_sum_up_
	KillWaves/Z g1,g2,g21neighgors,normal1,normal12,normal2,offsets
	return out
End
Static Function DeleteNaNs(g)
	Wave g
	Variable N=DimSize(g,0), M=DimSize(g,1)
	Variable iold, inew
	for (iold=0, inew=0;iold<N;iold+=1)				// check over every point in array
		if (numtype(g[iold][0])<2)						// a valid point
			g[inew][] = g[iold][q]
			inew += 1
		endif
	endfor

	N = inew
	Redimension/N=(N,M) g
	return N
End
//
//Function/S BoundaryStats(grains,i1,i2)// calc all information about a boundary bwteen two neighboring grains
//	Wave grains						// 3d wave with all of the grain numbers (e.g. gNeighbor)
//	Variable i1,i2						// grain numbers
//
//	String str
//	if (!WaveExists(grains) || numtype(i1+i2) || i1<0 || i2<0)
//		String grainName=NameOfWave(grains)
//		Prompt grainName,"3d wave with all grains",popup,ListMultiGrainMats()
//		String gizmoNote=GetGizmoNote()
//		if (numtype(i1))
//			str = StringByKey("grainBndry1",gizmoNote,"=",":")
//			i1 = str2num(str[5,Inf])
//			i1 = numtype(i1) ? 0 : max(i1,0)
//		endif
//		if (numtype(i2))
//			str = StringByKey("grainBndry2",gizmoNote,"=",":")
//			i2 = str2num(str[5,Inf])
//			i2 = numtype(i2) ? 0 : max(i2,0)
//		endif
//		Prompt i1, "first grain number"
//		Prompt i2, "second grain number"
//		DoPrompt "choose grains",grainName,i1,i2
//		if (V_flag)
//			return ""
//		endif
//		Wave grains = $grainName
//		printf "   	BoundaryStats(%s,%g,%g)\r",NameOfWave(grains),i1,i2
//	endif
//	if (!WaveExists(grains))
//		return ""
//	endif
//	WaveStats/Q grains
//	Variable maxGrain = V_max
//	if (numtype(i1+i2) || i1<0 || i2<0 || i1>maxGrain || i2>maxGrain || i1==i2)
//		if (ItemsInList(GetRTStackInfo(0))<2)
//			DoAlert 0,"invalid inputs"
//		endif
//		return ""
//	endif
//	// check that grains are neighbors
//	String list = makeBorderList(grains,i2)
//	if (WhichListItem(num2istr(i1),list)<0)
//		if (ItemsInList(GetRTStackInfo(0))<2)
//			sprintf str, "grains %d and %d are not neighbors\r",i1,i2
//			printf str
//			DoAlert 0, str
//		endif
//		return ""
//	endif
//	// all checks are done, now do the work
//
//	String fName = GetWavesDataFolder(grains,1)
//	Wave indicies = $(fname+"raw:indicies")
//	Wave XYZindex = $(fname+"raw:XYZindex")
//	Wave OrientMat = $(fname+"raw:OrientMat")
//	String gName = GetWavesDataFolder(grains,2)
//	Wave COMsingle = $(gName+"_COM")
//
//	Variable Nx,Ny,Nz
//	Nx=DimSize(XYZindex,0)
//	Ny=DimSize(XYZindex,1)
//	Nz=DimSize(XYZindex,2)
//	String out=""
//	out=ReplaceNumberByKey("i1",out,i1,"=")
//	out=ReplaceNumberByKey("i2",out,i2,"=")
//
//	Wave g1 = $getGrainCoords(grains,i1)		// list of grain coords for each grain
//	Wave g2 = $getGrainCoords(grains,i2)
//	Variable N1=DimSize(g1,0), N2=DimSize(g2,0)
//	Variable single1 = grainCOMsingle(g1)		// single index to a measured point nearest to COM
//	Variable single2 = grainCOMsingle(g2)
//
//	Variable totalAngle = angleBetweenPointsSymReduced(single1,single2)
//	out=ReplaceNumberByKey("totalAngle",out,totalAngle,"=")
//
//	Make/N=(N2)/T/O g21neighgors			// contains a list of neighbors to each grain in g2 (neighbors are in g2)
//	g21neighgors = ""
//	Variable i,j, ix0,iy0,iz0,  ix,iy,iz
//	Variable Nbndry1=0, Nbndry2=0, Nadd
//
//	str=UniqueName("offsets",1,0)
//	// offsets[][] contains the 6 voxels to the 1st nearest neighbors (faces in common)
//	Make/N=(6,3) $str
//	Wave offsets=$str
//	offsets[0][0]= {-1,1,  0,0,  0,0}
//	offsets[0][1]= {  0,0,-1,1,  0,0}
//	offsets[0][2]= {  0,0,  0,0,-1,1}
//	Variable Noffset = DimSize(offsets,0)
//
//	for (i=0;i<N2;i+=1)
//		Nadd = 0
//		for (j=0;j<Noffset;j+=1)				// loop through the offsets
//			ix = g2[i][0] + offsets[j][0]
//			iy = g2[i][1] + offsets[j][1]
//			iz = g2[i][2] + offsets[j][2]
//			if (numtype(ix+iy+iz) || ix>=Nx || iy>=Ny || iz>=Nz || ix<0 || iy<0 || iz<0)
//				continue						// point is NaN or out of range
//			elseif (grains[ix][iy][iz]==i1)
//				sprintf str "%d,%d,%d",ix,iy,iz
//				g21neighgors[i] = AddListItem(str, g21neighgors[i])
//				Nadd = 1
//			endif
//		endfor
//		Nbndry2 += Nadd
////if (Nadd)
////print i,"    ",g2[i][0],g2[i][1],g2[i][2],"    ",g21neighgors[i]
////endif
//	endfor
//
//	Variable jlast=-1							// remove all non-surface points from g2 and compress
//	for (i=0;i<Nbndry2;i+=1)
//		for (j=jlast+1;j<N2;j+=1)
//			if (strlen(g21neighgors[j])>0)
//				g2[i][0,2] = g2[j][q]
//				jlast = j
//				break
//			endif
//		endfor
//	endfor
//	Redimension/N=(Nbndry2,3) g2
//
//	Make/N=(N1)/O g1_flag_temp__
//	g1_flag_temp__ = 0
//	Variable k
//	for (i=0;i<N2;i+=1)
//		if (strlen(g21neighgors[i])>0)
//			str = g21neighgors[i]
//			for (j=0;j<ItemsInList(str);j+=1)	// check each triplet in g21neighgors[i]
//				ix = str2num(StringFromList(0,str,","))
//				iy = str2num(StringFromList(1,str,","))
//				iz = str2num(StringFromList(2,str,","))
//				for (k=0;k<N1;k+=1)
//					if (g1[k][0]==ix && g1[k][1]==iy && g1[k][2]==iz)
//						g1_flag_temp__[k] = 1
//						break
//					endif
//				endfor
//			endfor
//		endif
//	endfor
//	Nbndry1 = sum(g1_flag_temp__)
//
//	jlast=-1									// remove all non-surface points from g1 and compress
//	for (i=0;i<Nbndry1;i+=1)
//		for (j=jlast+1;j<N1;j+=1)
//			if (g1_flag_temp__[j])
//				g1[i][0,2] = g1[j][q]
//				jlast = j
//				break
//			endif
//		endfor
//	endfor
//	Redimension/N=(Nbndry1,3) g1
//	KillWaves g1_flag_temp__
//	out=ReplaceNumberByKey("Nboundary1",out,Nbndry1,"=")
//	out=ReplaceNumberByKey("Nboundary2",out,Nbndry2,"=")
//	//
//	// g1 and g2 now contain only points in the interface (i.e. points that border the other grain by a face)
//
//	Make/N=3/D/O normal1,normal2, normal12			// vector in the average surface normal for each grain
//	normal1 = 0
//	normal2 = 0
//	Variable x0,y0,z0, dx,dy,dz,xc,yc,zc
//	x0 = DimOffset(grains,0)
//	y0 = DimOffset(grains,1)
//	z0 = DimOffset(grains,2)
//	dx = DimDelta(grains,0)
//	dy = DimDelta(grains,1)
//	dz = DimDelta(grains,2)
//	if (Nbndry1>=3)
//		g1[][0] = g1[p][0]*dx + x0						// first do grain1
//		g1[][1] = g1[p][1]*dy + y0
//		g1[][2] = g1[p][2]*dz + z0
//		Make/N=(Nbndry1)/O temp_sum_up_
//		temp_sum_up_ = g1[p][0]  ;  xc = faverage(temp_sum_up_)
//		temp_sum_up_ = g1[p][1]  ;  yc = faverage(temp_sum_up_)
//		temp_sum_up_ = g1[p][2]  ;  zc = faverage(temp_sum_up_)
//		NormalOfPlane(g1,normal1)
////print NameOfWave(g1),xc,yc,zc,g1,normal1
//		normalize(normal1)
//		Duplicate/O g1 bndryScatterPoints
//	endif
//	if (Nbndry2>=3)
//		g2[][0] = g2[p][0]*dx + x0						// second do grain2
//		g2[][1] = g2[p][1]*dy + y0
//		g2[][2] = g2[p][2]*dz + z0
//		Make/N=(Nbndry2)/O temp_sum_up_
//		temp_sum_up_ = g2[p][0]  ;  xc = (xc+faverage(temp_sum_up_))/2
//		temp_sum_up_ = g2[p][1]  ;  yc = (yc+faverage(temp_sum_up_))/2
//		temp_sum_up_ = g2[p][2]  ;  zc = (zc+faverage(temp_sum_up_))/2
//		NormalOfPlane(g2,normal2)
//		normalize(normal2)
////print NameOfWave(g2),xc,yc,zc,g2,normal1
//		out=ReplaceNumberByKey("xc",out,xc,"=")
//		out=ReplaceNumberByKey("yc",out,yc,"=")
//		out=ReplaceNumberByKey("zc",out,zc,"=")
//	endif
//	Variable normalAngleErr = acos(max(min(MatrixDot(normal1,normal2),1),-1))
//	out=ReplaceNumberByKey("normalAngleErr",out,normalAngleErr,"=")
//	normal12 = (normal1+normal2)/2						// average normal to boundary
//	if (normalAngleErr*180/PI>5 && !stringmatch(StringFromList(0, GetRTStackInfo(0)),"GizmoAroundOneGrain"))
//		DoAlert 0, "WARNING in BoundaryStats(), angle between normals > 5¡"
//		normal12 = normal1
//	endif
//	sprintf str, "%g,%g,%g",normal12[0],normal12[1],normal12[2]
//	out=ReplaceStringByKey("bndryNormal",out,str,"=")
//
//	Variable maxhkl=12
//	if (normalize(normal12)>0.1)							// if normal12 is valid, then proceed to calculate the hkl1 & hkl2
//		// next find the hkl for each of the two grains that lies along normal12
//		Make/N=3/O/D hkl1,hkl2, normalXHF
//		Make/N=(3,3)/O mat1,mat2
//		// XYZ -> xyz,  	convert normal from beamline coords (XYZ) to the sample (xyz) == (X -H -F)
//		normalXHF[0] = normal12[0]
//		normalXHF[1] = -(normal12[1]+normal12[2])/sqrt(2)
//		normalXHF[2] = (normal12[1]-normal12[2])/sqrt(2)
//		i = COMsingle[i1]
//		//	print "total rotaion angle for mat1 is ",anglefromMatrix(i)
//		mat1 = OrientMat[i][p][q]
//		i = COMsingle[i2]
//		//	print "total rotaion angle for mat2 is ",anglefromMatrix(i)
//		mat2 = OrientMat[i][p][q]
//		MatrixOp/O hkl1 = mat1 x normalXHF			//  A*(xyz) = (hkl)
//		integerDirectionApprox(hkl1,maxhkl)
//		MatrixOp/O hkl2 = mat2 x normalXHF
//		integerDirectionApprox(hkl2,maxhkl)
//		Make/N=3/O/D axisOfMat						// get the axis of the rotation
//		MatrixOp/O rot = Inv(mat1) x mat2				// rot is the total rotation matrix
//		Variable angle = axisOfMatrix(rot,axisOfMat,squareUp=1)	// compute axis and angle
//		// xyz -> XYZ,  	convert from sample (xyz) to beamline coords (XYZ), and again (xyz) == (X -H -F)
//		Variable ysample = axisOfMat[1], zsample = axisOfMat[2]
//		axisOfMat[1] = (-ysample+zsample)/sqrt(2)
//		axisOfMat[2] = (-ysample-zsample)/sqrt(2)
//		Variable angleBetweenAxes = acos(MatrixDot(axisOfMat,normal12))*180/PI
//		angleBetweenAxes = min(angleBetweenAxes, 180-angleBetweenAxes)
//		sprintf str, "%g,%g,%g",hkl1[0],hkl1[1],hkl1[2]
//		out=ReplaceStringByKey("hkl1",out,str,"=")
//		sprintf str, "%g,%g,%g",hkl2[0],hkl2[1],hkl2[2]
//		out=ReplaceStringByKey("hkl2",out,str,"=")
//	endif
//
//	if (ItemsInList(GetRTStackInfo(0))<2)
//		printf "   for the boundary between grains %g and %g in '%s':\r",i1,i2,NameOfWave(grains)
//		printf "   angle between grains %d and %d is %g¡,  rotated about the (%g, %g, %g)\r",i1,i2,totalAngle,axisOfMat[0],axisOfMat[1],axisOfMat[2]
//		printf "   angle between axis of rotation and interface normal is %g¡\r",angleBetweenAxes
//		printf "   grain 1 has %d points in the surface, and grain 2 has %d in the surface\r",Nbndry1,Nbndry2
//		printf "   surface normal (average) between grains 1 & 2 is  {%g, %g, %g}   angle between individual normals is %.3g deg\r",normal12[0],normal12[1],normal12[2],normalAngleErr*180/PI
//		printf "         surface normals are: grain1-> {%g, %g, %g}   and   grain2-> {%g, %g, %g}\r",normal1[0],normal1[1],normal1[2],normal2[0],normal2[1],normal2[2]
//		printf "   approx surface normals are grain1->(%g %g %g)   and   grain2->(%g %g %g) (with an hkl limit of %g)\r",hkl1[0],hkl1[1],hkl1[2],hkl2[0],hkl2[1],hkl2[2],maxhkl
//		hkl1 = abs(hkl1) ; Sort hkl1, hkl1
//		hkl2 = abs(hkl2) ; Sort hkl2, hkl2
//		printf "   							or grain1->(%g %g %g)   and   grain2->(%g %g %g)\r",hkl1[0],hkl1[1],hkl1[2],hkl2[0],hkl2[1],hkl2[2]
//	endif
//	KillWaves/Z hkl1,hkl2,mat1,mat2, axisOfMat, rot, normalXHF
//	KillWaves/Z W_coef,W_sigma,W_ParamConfidenceInterval,temp_sum_up_
//	KillWaves/Z g1,g2,g21neighgors,normal1,normal12,normal2,offsets
//	return out
//End


Function TiltTwistOfBoundary(Ain,B,normal,twist,tilt)	// return the twist and tilt angles for this boundary
	Wave Ain,B								// the two orientation matrices (in sample representation, NOT beamline)
	Wave normal							// normal to interface (sample representation, NOT beamline!)
	Variable &twist,&tilt					// angles returned (degree)

	String wName
	wName = UniqueName("hkl",1,0) ; Make/N=3/D $wName ; Wave hA=$wName
	wName = UniqueName("hkl",1,0) ; Make/N=3/D $wName ; Wave hB=$wName
	wName = UniqueName("SymMatrix",1,0)	;  Make/N=(3,3)/D/O $wName  ;  Wave symOp = $wName
	wName = UniqueName("Amat",1,0)	;  Make/N=(3,3)/D/O $wName  ;  Wave A = $wName
	A = Ain
	BestSymOp(A,B,symOp)
	MatrixOp/O A = Inv(symOp) x A			// the new A, symmetry reduced

	MatrixOp/O hA = Normalize(A x normal)
	MatrixOp/O hB = Normalize(B x normal)

	if (MatrixDot(hA,hB)<0)
		hB *= -1
	endif
	Cross hA, hB
	Wave W_Cross=W_Cross
	Variable asine = normalize(W_Cross)
	tilt = asin(asine)*180/PI

	wName = UniqueName("rot",1,0)
	Make/N=3/D $wName
	Wave tiltMat=$wName
	rotationMatAboutAxis(W_Cross,-tilt,tiltMat)

	MatrixOp/O hA = tiltMat x hA		// this brings hA == hB
	wName = UniqueName("matrix1D",1,0)  ;  Make/N=(1)/D/O $wName  ;  Wave cosine11 = $wName
	MatrixOp/O cosine11 = (Trace(tiltMat x A x Inv(B))-1)/2
	Variable cosine = min(max(-1,cosine11[0]),1)			// ensure cosine in range [-1,1]
	twist = acos(cosine)*180/PI
	KillWaves/Z hA,hB, W_Cross, cosine11, tiltMat, symOp, A
End



Function integerDirectionExact(hkl,maxInt)		// changes vector hkl to be an approxiamtely integer vector, like integerDirectionApprox(), but exact
	Wave hkl
	Variable maxInt
	String wName = UniqueName("hkl",1,0)
	Make/N=(numpnts(hkl))/O/D $wName
	Wave hkli = $wName
	WaveStats/Q/M=1 hkl
	Variable biggest = max(abs(V_max),abs(V_min))
	hkli = round(hkl[p]*maxInt/biggest)
	Variable factor = gcf(GetWavesDataFolder(hkli,2))
	hkl /= factor
	KillWaves/Z hkli
End
//
Function integerDirectionApprox(hkl,maxInt)		// changes vector hkl to be an integer vector, approximately in same direction
	Wave hkl
	Variable maxInt
	WaveStats/Q/M=1 hkl
	Variable biggest = max(abs(V_max),abs(V_min))
	hkl = round(hkl[p]*maxInt/biggest)
	Variable factor = gcf(GetWavesDataFolder(hkl,2))
	hkl /= factor
End
//
//Function integerDirectionExact(hkl,maxInt)		// changes vector hkl to be an approxiamtely integer vector, like integerDirectionApprox(), but exact
//	Wave hkl
//	Variable maxInt
//	Make/N=3/O/D test_hkl_
//	normalize(hkl)
//	test_hkl_ = hkl
//	hkl = round(maxInt*hkl[p])
//	Variable factor = gcf(GetWavesDataFolder(hkl,2))
//	hkl /= factor
//	factor = norm(hkl)/norm(test_hkl_)
//	hkl = test_hkl_*factor
//	KillWaves/Z test_hkl_
//End
//	//
//Function integerDirectionApprox(hkl,maxInt)		// changes vector hkl to be an integer vector, approximately in same direction
//	Wave hkl
//	Variable maxInt
//	normalize(hkl)
//	hkl = round(maxInt*hkl[p])
//	Variable factor = gcf(GetWavesDataFolder(hkl,2))
//	hkl /= factor
//End
//
// this revision of gcf uses the the "PrimeFactors" command available with Igor 5
Static Function gcf(factors)				// find greatest common factor of a wave whose values are all integers
	String factors					// either a list of factors or the name of a wave containing the factors

	String tempName=""
	Variable N,i
	if (exists(factors)==1 && strsearch(factors,";",0)<0)		// passed a wave name
		Wave wav=$factors
	else
		N = ItemsInList(factors)								// a list was passed
		if (N<1)
			return NaN
		endif
		tempName = UniqueName("tempPrimes",1,0)
		Make/N=(N) $tempName
		Wave wav=$tempName
		for (i=0;i<N;i+=1)
			wav[i] = str2num(StringFromList(i,factors))
		endfor
	endif
	N=numpnts(wav)
	if (N<1)
		return NaN
	endif

	String primeListsName = UniqueName("ListOfPrimes",1,0)
	Make/N=(N)/T $primeListsName
	Wave/T primeLists=$primeListsName
	primeLists = ""
	String allPrimes="", item
	Variable j,gcf
	for (i=0;i<N;i+=1)
		PrimeFactors/Q  wav[i]
		Wave W_PrimeFactors=W_PrimeFactors
		for (j=0;j<numpnts(W_PrimeFactors);j+=1)
			primeLists[i]=AddListItem(num2istr(W_PrimeFactors[j]), primeLists[i])
		endfor
		allPrimes+= primeLists[i]
		if (wav[i]==0)									// ensure "0" is there for a 0
			primeLists[i]=AddListItem("0", primeLists[i])
		endif
	endfor

	Variable m,found,jN = ItemsInList(allPrimes)
	j = 0
	gcf=1
	for (j=0;j<jN;j+=1)								// loop over all of the primes
		item = StringFromList(j,allPrimes)
		found = 1
		for (i=0;i<N;i+=1)
			m = WhichListItem(item,primeLists[i])
			if (m>=0)
				primeLists[i] = RemoveListItem(m,primeLists[i])
			elseif (cmpstr(primeLists[i],"0;"))
				found=0
			endif
		endfor
		if (found)
			gcf *= str2num(item)
		endif
	endfor
	KillWaves/Z primeLists,W_PrimeFactors, $tempName
	return gcf
End



Function NormalOfPlane(xyzPos,normal)
	Wave xyzPos
	Wave normal

	if (DimSize(xyzPos,0)<3)								// need at least 3 points to define a plane
		normal = NaN
		return 1
	endif
	Variable X1,Y1,Z1,XX,YY,ZZ,XY,XZ,YZ
	Variable N = DimSize(xyzPos,0)
	Make/N=(N)/O/D temp_for_plane_
	Wave sumTemp = temp_for_plane_
	sumTemp = xyzPos[p][0]  ;  X1 = sum(sumTemp)
	sumTemp = xyzPos[p][1]  ;  Y1 = sum(sumTemp)
	sumTemp = xyzPos[p][2]  ;  Z1 = sum(sumTemp)
	sumTemp = xyzPos[p][0]*xyzPos[p][0]  ;  XX = sum(sumTemp)
	sumTemp = xyzPos[p][1]*xyzPos[p][1]  ;  YY = sum(sumTemp)
	sumTemp = xyzPos[p][2]*xyzPos[p][2]  ;  ZZ = sum(sumTemp)
	sumTemp = xyzPos[p][0]*xyzPos[p][1]  ;  XY = sum(sumTemp)
	sumTemp = xyzPos[p][0]*xyzPos[p][2]  ;  XZ = sum(sumTemp)
	sumTemp = xyzPos[p][1]*xyzPos[p][2]  ;  YZ = sum(sumTemp)

	Make/N=3/D/O con1={X1,Y1,Z1}
	Make/N=(3,3)/D/O mat2
	mat2[0][0]=XX  ;  mat2[0][1]=XY  ;  mat2[0][2]=XZ	// note, mat2 is symmetric
	mat2[1][0]=XY  ;  mat2[1][1]=YY  ;  mat2[1][2]=YZ
	mat2[2][0]=XZ  ;  mat2[2][1]=YZ  ;  mat2[2][2]=ZZ

	MatrixLLS  mat2  con1
	Wave M_B = M_B
	Redimension/N=3/D normal
	normal = M_B

//	print xyzPos
//	print "x1,y1,z1=",x1,y1,z1
//	printf "X1=%g,   Y1=%g,   Z1=%g\r",x1,y1,z1
//	printf "XX=%g,   YY=%g,   ZZ=%g,      XY=%g,   XZ=%g,   YZ=%g\r",XX,YY,ZZ,XY,XZ,YZ
//	print mat2
//	MatrixMultiply mat2, M_B
//	Wave M_product=M_product
//	print con1,M_product
//	KillWaves/Z M_product
	KillWaves/Z M_B,M_A,temp_for_plane_,mat2,con1
	return 0
End



Function/S makeBorderList(grains,grainNum)	// it includes grainNum in the list too
	Variable grainNum							// find the grains that border grainsNum
	Wave grains								// 3d wave with all of the grain numbers (e.g. gNeighbor)

	if (!WaveExists(grains) || numtype(grainNum) || grainNum<0)
		String grainName=NameOfWave(grains)
		Prompt grainName,"3d wave with all grains",popup,ListMultiGrainMats()
		grainNum = numtype(grainNum)==2 ? 0 : max(grainNum,0)
		Prompt grainNum, "grain number"
		DoPrompt "choose grain",grainName,grainNum
		if (V_flag)
			return ""
		endif
		Wave grains = $grainName
	endif
	if (!WaveExists(grains))
		return ""
	endif
	WaveStats/Q grains
	Variable maxGrain = V_max
	if (numtype(grainNum) || grainNum<0 || grainNum>maxGrain)
		return ""
	endif

	Make/N=(maxGrain+1)/O tempGrainFoundList_,tempGrainFoundIndex_
	tempGrainFoundList_ = 0
	tempGrainFoundIndex_ = p
	String wName = getGrainCoords(grains,grainNum)		// create a wave containing the [ix][iy][iz] coords of each point in grain
	Wave xyz = $wName
	Variable Nx,Ny,Nz
	Nx=DimSize(grains,0)
	Ny=DimSize(grains,1)
	Nz=DimSize(grains,2)

	// check the neighbor of each xyz[i][3] for new neighbors
	Variable Nxyz = DimSize(xyz,0)
	Variable i,j,ix,iy,iz
	for (i=0;i<Nxyz;i+=1)
		ix = xyz[i][0]
		iy = xyz[i][1]
		iz = xyz[i][2]

		j = grains[max(ix-1,0)][iy][iz]
		j = numtype(j) ? maxGrain : j
		tempGrainFoundList_[j] += 1

		j = grains[min(ix+1,Nx-1)][iy][iz]
		j = numtype(j) ? maxGrain : j
		tempGrainFoundList_[j] += 1

		j = grains[ix][max(iy-1,0)][iz]
		j = numtype(j) ? maxGrain : j
		tempGrainFoundList_[j] += 1

		j = grains[ix][min(iy+1,Ny-1)][iz]
		j = numtype(j) ? maxGrain : j
		tempGrainFoundList_[j] += 1

		j = grains[ix][iy][max(iz-1,0)]
		j = numtype(j) ? maxGrain : j
		tempGrainFoundList_[j] += 1

		j = grains[ix][iy][min(iz+1,Nz-1)]
		j = numtype(j) ? maxGrain : j
		tempGrainFoundList_[j] += 1
	endfor
	tempGrainFoundList_[maxGrain]= -1

	Sort/R tempGrainFoundList_, tempGrainFoundList_, tempGrainFoundIndex_
	Variable N = BinarySearch(tempGrainFoundList_,0.5)+1
	ReDimension/N=(N) tempGrainFoundIndex_
	Sort tempGrainFoundIndex_,tempGrainFoundIndex_
	String str=""
	for (i=0;i<N;i+=1)
		str = AddListItem(num2istr(tempGrainFoundIndex_[i]),str,";",Inf)
	endfor
	if (ItemsInList(GetRTStackInfo(0))<2)
		printf "grain %d has %d neighbors, they are  {%s}\r",grainNum,N-1,str
	endif
	KillWaves/Z tempGrainFoundList_, xyz, tempGrainFoundIndex_
	return str
End

//	=========================================================================
//								End of Grain Boundary Analysis
//	=========================================================================
//	=========================================================================
//	=========================================================================

//  ===================================================================
//  ===================================================================
//  ===================================================================
//		graph and table routines
//  ===================================================================

Function Table_GrainOrientations(grains) : Table	// put up a table of a multigrain 3d wave
	Wave grains								// 3d wave with all of the grain numbers (e.g. gNeighbor)
	if (!WaveExists(grains))
		String grainName=NameOfWave(grains)
		Prompt grainName,"3d wave with all grains",popup,ListMultiGrainMats()
		DoPrompt "choose grains",grainName
		if (V_flag)
			return 1
		endif
		Wave grains = $grainName
	endif
	if (!WaveExists(grains))
		return 1
	endif
	String tableName = "GrainOrientations_"+NameOfWave(grains)
	if (strlen(WinList(tableName,"","WIN:2")))
		DoWindow/F $tableName
		return 0
	endif

	String name=GetWavesDataFolder(grains,2)
	String  rawfldr=GetWavesDataFolder(grains,1)+"raw:"
	Wave vol = $(name+"_Vol")
	String com = name+"_COM"	// REMOVE
	Wave xc = $(name+"_xc")					// single index to COM voxel (micron)
	Wave yc = $(name+"_yc")
	Wave zc = $(name+"_zc")
	Wave Rodrigues = $(name+"_Rodrigues")	// Rodrigues vector

	Edit/W=(3,63,536,316)/K=1 vol,xc,yc,zc,Rodrigues
	DoWindow/C $tableName
	ModifyTable width=70
	ModifyTable width(Point)=40
	ModifyTable width(Rodrigues)=64,format(Rodrigues)=3,digits(Rodrigues)=5
End


Function Table_AllGrains(grains) : Table	// put up a table of a multigrain 3d wave
	Wave grains						// 3d wave with all of the grain numbers (e.g. gNeighbor)
	if (!WaveExists(grains))
		String grainName=NameOfWave(grains)
		Prompt grainName,"3d wave with all grains",popup,ListMultiGrainMats()
		DoPrompt "choose grains",grainName
		if (V_flag)
			return 1
		endif
		Wave grains = $grainName
	endif
	if (!WaveExists(grains))
		return 1
	endif
	Edit/W=(100,47,1003,793) grains
	ModifyTable elements(gNeighbor)=(-3,6,-2,-1)
End

Function Graph_grain_Vol(grainVol) : Table	// put up a table of the grain volumes
	Wave grainVol						// 3d wave with all of the grain numbers (e.g. gNeighbor)
	if (!WaveExists(grainVol))
		String grainVolName=NameOfWave(grainVol)
		Prompt grainVolName,"3d wave with all grainVol waves",popup,WaveList("*_Vol",";","DIMS:1")
		DoPrompt "choose grainVol",grainVolName
		if (V_flag)
			return 1
		endif
		Wave grainVol = $grainVolName
	endif
	if (!WaveExists(grainVol))
		DoAlert 1, "make a wave temporarily?"
		if (V_flag)
			Make/N=2 gNeighbor_Vol
			Wave grainVol=gNeighbor_Vol
		else
			return 1
		endif
	endif
	Display /W=(542,44,1150,440) grainVol
	ModifyGraph gfMult=130,mode=4,marker=19,lStyle=1,rgb=(0,0,65000),msize=2
	ModifyGraph log(left)=1,tick=2,mirror=1,minor=1,lowTrip=0.001,standoff=0
	SetAxis/N=1/E=1/A left
	Cursor/P A $StringFromList(0,TraceNameList("",";",1)) numpnts(grainVol)/4
	ShowInfo
End

Function Graph_HistTotalAngle(HistTotalAngle) : Graph
	Wave HistTotalAngle					// 1d wave with a histogram of the angles
	if (!WaveExists(HistTotalAngle))
		String histName=NameOfWave(HistTotalAngle)
		Prompt histName,"3d wave with all HistTotalAngle",popup,WaveList("*HistTotalAngle",";","DIMS:1")
		DoPrompt "choose HistTotalAngle",histName
		if (V_flag)
			return 1
		endif
		Wave HistTotalAngle = $histName
	endif
	if (!WaveExists(HistTotalAngle))
		Make/N=2 HistTotalAngle
	endif
	Display /W=(5,44,540,441) HistTotalAngle
	ModifyGraph gfMult=130,tick=2,mirror=1,minor=1,lowTrip=0.001,standoff=0
	SetAxis/E=1/A left
	ShowInfo
	ValDisplay sumOfHistogram,pos={75,1},size={80,14},title="sum",fSize=10
	ValDisplay sumOfHistogram,format="%d",frame=2,limits={0,0,0},barmisc={0,1000}
	String str
	sprintf str,"ValDisplay sumOfHistogram,value= %d",sum(HistTotalAngle)
	Execute str
End

Function Table_XYZindex(XYZIndex) : Table	// put up a table of XYZ indicies
	Wave XYZIndex							// 2d wave with all of the grain numbers (e.g. gNeighbor)
	if (!WaveExists(XYZIndex))
		Wave XYZIndex = $"XYZIndex"
	endif
	if (!WaveExists(XYZIndex))
		Wave XYZIndex = $":raw:XYZIndex"
	endif
	if (!WaveExists(XYZIndex))
		return 1
	endif
	Edit/W=(12,46,888,636) XYZindex
	ModifyTable width(Point)=34,elements(XYZindex)=(-3,0,-2,1)
End

Function Table_GrainStatistics(grains) : Table	// put up a table of a multigrain 3d wave
	Wave grains								// 3d wave with all of the grain numbers (e.g. gNeighbor)
	if (!WaveExists(grains))
		String grainName=NameOfWave(grains)
		Prompt grainName,"3d wave with all grains",popup,ListMultiGrainMats()
		DoPrompt "choose grains",grainName
		if (V_flag)
			return 1
		endif
		Wave grains = $grainName
	endif
	if (!WaveExists(grains))
		return 1
	endif

	String name=GetWavesDataFolder(grains,2)
	String vol = name+"_Vol"
	String com = name+"_COM"
	String xc = name+"_xc"
	String yc = name+"_yc"
	String zc = name+"_zc"
	String ixlo = name+"_ixlo"
	String ixhi = name+"_ixhi"
	String iylo = name+"_iylo"
	String iyhi = name+"_iyhi"
	String izlo = name+"_izlo"
	String izhi = name+"_izhi"

	Edit/W=(27,163,704,849) $vol,$com,$xc,$yc,$zc,$ixlo,$ixhi,$iylo,$iyhi,$izlo,$izhi
	ModifyTable width(Point)=40,width[1]=68,width[2]=64,width[3]=54
	ModifyTable width[4]=54,width[5]=54,width[6]=54,width[7]=54
	ModifyTable width[8]=54,width[9]=54,width[10]=54,width[11]=54
End

//  ===================================================================
//		End of graph and table routines
//  ===================================================================
//  ===================================================================
//  ===================================================================



//Window Table_grainNo() : Table
//	PauseUpdate; Silent 1
//	if (strlen(WinList("TableGrainNo","","WIN:2")))
//		DoWindow /F TableGrainNo
//		return
//	endif
//	String fldrSav= GetDataFolder(1)
//	SetDataFolder root:raw:
//	Edit/W=(5,44,877,682)/K=1 grainNo
//	ModifyTable width(Point)=36,elements(grainNo)=(-3,9,-2)
//	DoWindow /C TableGrainNo
//	SetDataFolder fldrSav
//EndMacro


/// assign a grain number to each voxel in the volume, returns the full path name of the 3d array of grain id's.  
//// number of id'd grains is in the note of grainNo
//Function/S AssignGrainNums(orients3d,thresh)
//	Wave orients3d
//	Variable thresh				// angle threhold for a grain boundary
//	if (WaveDims(orients3d)!=4)
//		Abort  " dim of orient3d is not 4"
//	endif
//
//	if (numtype(thresh) || thresh<=0)
//		thresh = (numtype(thresh)||(thresh<=0)) ? 10 : thresh
//		Prompt thresh, "minimum angle between two grains (¡)"
//		DoPrompt "min angle between grains",thresh
//		if (V_flag)
//			return ""
//		endif
//	endif
//	if (numtype(thresh) || thresh<=0)
//		DoAlert 0, "thresh = "+num2str(thresh)+" not allowed"
//		return ""
//	endif
//	DoAlert 1, "this will take ~1/2 hour, proceed"
//	if (V_flag!=1)
//		return ""
//	endif
//	printf "		AssignGrainNums(%s,%g)\r",GetWavesDataFolder(orients3d,2),thresh
//
//	Variable nx, ny, nz
//	nx = DimSize(orients3d,0)
//	ny = DimSize(orients3d,1)
//	nz = DimSize(orients3d,2)
//	String grainNoStr = GetWavesDataFolder(orients3d,1)+"grainNo"
//	Make/N=(nx,ny,nz)/O $grainNoStr,checkNeighbor_temp_
//	Wave grainNo = $grainNoStr		// this contains the grain number for each voxel
//	Wave checkNeighbor=checkNeighbor_temp_		// 3d array of flags.  true flags a voxel whose grain id was set, but its neighbors have not been checked
//	grainNo  = -1						// set all voxels to un assigned values (good values are all ³0)
//	checkNeighbor = 0					// set all check neighbor flags to unset
//	Variable i,j,k						// indicies for one grain
//	Variable gindex = 0					// grain index
//	Variable tickStart = ticks
//	do
//		if (GetNextGrainToIdentify(grainNo,i,j,k))	// return the ijk of a grain that has not been id'd, if true, none found
//			break
//		endif
//		if (numtype(orients3d[i][j][k]))			// invalid orientation at ijk, skip it
//			grainNo[i][j][k] = NaN
//			checkNeighbor[i][j][k] = 0
//			continue
//		endif
//		grainNo[i][j][k] = gindex
//		checkNeighbor[i][j][k] = 1					// neighbors for this grain have not been checked, obviously
//
//printf "starting to fill out grain %d at (ijk)=(%d %d %d)     after %s\r",gindex,i,j,k,num2sexigesmal((ticks-tickStart)/60.15,0)
//		// this loop always starts with an ijk of a ligitimate grain (with grain number of gindex)
//		do											// expanding out from ijk, find the rest of this grain
//			if (AssignGrainNeighbors(orients3d,grainNo,checkNeighbor,i,j,k,thresh))	// check 6 neighbors of ijk & set flags
//				Abort "We should never get a bad return from 'AssignGrainNeighbors()'"
//				return ""
//			endif
//			if (GetNextUnCheckedNeighbor(checkNeighbor,i,j,k)	)// return the ijk of a grain that has the checkNeighbor set (=1)
//				break
//			endif
//		while(1)
//		gindex += 1									// increment grain index
//	while(1)
//	KillWaves/Z checkNeighbor_temp_
//	printf "identified %d grains, execution took %s\r",gindex,num2sexigesmal((ticks-tickStart)/60.15,0)
//
//	Note/K $grainNoStr
//	String noteStr = note(orients3d)
//	noteStr = ReplaceNumberByKey("numGrains",noteStr,gindex,"=")
//	Note $grainNoStr noteStr
//
//	return GetWavesDataFolder(grainNo,2)
//End


//Function/S GetNbiggestGrains(grainNo,minSize) // returns list of values in order of grain size, down to grains of size minSize
//	Wave grainNo
//	Variable minSize
//	if (numtype(minSize) || minSize<1 || !WaveExists(grainNo))
//		minSize = (minSize>=1) ? minSize : 20
//		String grainNumName=""
//		Prompt grainNumName, "array of gain id's",popup,WaveList("*grainNo",";","DIMS:3")
//		Prompt minSize, "min no. of voxels in a grain"
//		if (WaveExists(grainNo))
//			DoPrompt "min grain size to count",minSize
//		else
//			DoPrompt "grain to count",grainNumName,minSize
//			Wave grainNo = $grainNumName
//		endif
//		if (V_flag)
//			return ""
//		endif
//		printf "         GetNbestGrains(%s,%d)\r", GetWavesDataFolder(grainNo,2),minSize
//	endif
//
//	String noteStr = note(grainNo)
//	Variable nGrains = NumberByKey("numGrains",noteStr,"=")
//
//	String Ssizes=UniqueName("VoxSizes",1,0)
//	String Sindex=UniqueName("index",1,0)
//	Make/O/N=(nGrains) $Ssizes,$Sindex
//	Wave sizes = $Ssizes
//	Wave index = $Sindex
//	Histogram/B={0,1,nGrains} grainNo,sizes
//	index = p
//	Sort/R sizes, sizes,index
//
//	String grainList=""
//	Variable i = 0
//	do
//		if (sizes[i]<minSize)
//			break
//		endif
//		grainList += num2istr(index[i])+","+num2istr(sizes[i])+";"
//		i += 1
//	while(1)
//	KillWaves/Z sizes,index
//	return grainList
//End



//Static Function GetNextGrainToIdentify(grainNo,i,j,k)		// return the ijk of a grain that has not yet been id'd
//	Wave grainNo				// 3d array  the number that identifies the grain (-1 flags not yet id'd grains)
//	Variable &i,&j,&k
//	Variable nx,ny,nz
//	nx = DimSize(grainNo,0)
//	ny = DimSize(grainNo,1)
//	nz = DimSize(grainNo,2)
//	if (nz<1)
//		return -1				// error, grainNo has the wrong dimensions
//	endif
//
//	for (k=0;k<nz;k+=1)
//		for (j=0;j<ny;j+=1)
//			for (i=0;i<nx;i+=1)
//				if (grainNo[i][j][k]<0)
//					return 0		// ijk is the location of a grain whose id is still negative
//				endif
//			endfor
//		endfor
//	endfor
//	return 1						// all grains are identified with a grain number
//End


//Static Function GetNextUnCheckedNeighbor(checkNeighbor,i,j,k)		// return the ijk of a grain that has the checkNeighbor set (=1)
//	Wave checkNeighbor		// 3d array of flags.  true flags a voxel whose grain id was set, but its neighbors have not been checked
//	Variable &i,&j,&k
//	Variable nx,ny,nz
//	nx = DimSize(checkNeighbor,0)
//	ny = DimSize(checkNeighbor,1)
//	nz = DimSize(checkNeighbor,2)
//	if (nz<1)
//		return -1				// error, checkNeighbor has the wrong dimensions
//	endif
//
//	for (k=0;k<nz;k+=1)
//		for (j=0;j<ny;j+=1)
//			for (i=0;i<nx;i+=1)
//				if (checkNeighbor[i][j][k])
//					return 0		// ijk is the location of a grain whose checkNeighbor flag is set.
//				endif
//			endfor
//		endfor
//	endfor
//	return 1						// no grains have a checkNeighbor flag set
//End








//  ===================================================================
//  ===================================================================
//  ===================================================================
//		read data file and process
//  ===================================================================

Function PrepareData()
	String fldrSav= GetDataFolder(1)

	initGrainGizmo()
	String fullFilePath=""
	Variable i,readin=1
	if (exists(":raw:XX")==1)
		SetDataFolder :raw
		DoAlert 2, "re-read in the data"
		if (V_flag==3)		// "cancel" clicked
			SetDataFolder fldrSav
			return 1
		endif
		readin = (V_flag==1)
		if (readin)
			printf "re-reading and overwriting existing data\r"
		endif
	else
		String dataFolderStr="", folderList=StringByKey("FOLDERS", DataFolderDir(1))
		folderList = RemoveFromList("Packages",folderList,",")
		folderList += SelectString(strlen(folderList),"",",")+"new folder"
		folderList = ReplaceString(",",folderList,";")
		Prompt dataFolderStr,"data folder with data",popup,folderList
		DoPrompt "data folder",dataFolderStr
		if (V_flag)
			SetDataFolder fldrSav
			return 1
		endif
		if (stringmatch(dataFolderStr,"new folder"))
			Variable refNum
			Open/D/M="data file"/P=home/R refNum		// S_fileName will contain the full file name
			fullFilePath = S_fileName
			if (strlen(fullFilePath)<1)
				DoAlert 0, "No file selected"
				SetDataFolder fldrSav
				return 1
			endif
			dataFolderStr = ParseFilePath(3, fullFilePath, ":", 0, 0)
			i = strsearch(dataFolderStr, "_Orients", Inf,1)
			if (i>0)
				dataFolderStr = dataFolderStr[0,i-1]
			endif
			dataFolderStr = CleanupName(dataFolderStr,0)
			if (DataFolderExists(":"+dataFolderStr))
				DoAlert 0, "Data folder "+dataFolderStr+" already exists"
				SetDataFolder fldrSav
				return 1
			endif
			NewDataFolder/S $(":"+dataFolderStr)
			NewDataFolder/O/S :raw
			readin = 1
		else
			SetDataFolder (":"+dataFolderStr)
			NewDataFolder/O/S :raw
			DoAlert 2, "re-read in the data"
			if (V_flag==3)		// "cancel" clicked
				SetDataFolder fldrSav
				return 1
			endif
			readin = (V_flag==1)
			if (readin)
				printf "re-reading and overwriting existing data\r"
			endif
		endif
	endif
	if (readin)
		String columnInfoStr
		columnInfoStr    =	"C=1,N=XX;C=1,N=YY;C=1,N=ZZ;"
		columnInfoStr +=	"C=1,N=Hx;C=1,N=Hy;C=1,N=Hz;"
		columnInfoStr +=	"C=1,N=Kx;C=1,N=Ky;C=1,N=Kz;"
		columnInfoStr +=	"C=1,N=Lx;C=1,N=Ly;C=1,N=Lz;"
//		columnInfoStr +=	"C=1,N=Hx;C=1,N=Kx;C=1,N=Lx;"	// this done because IDL says A[col,row], backwards, WRONG!
//		columnInfoStr +=	"C=1,N=Hy;C=1,N=Ky;C=1,N=Ly;"
//		columnInfoStr +=	"C=1,N=Hz;C=1,N=Kz;C=1,N=Lz;"
		columnInfoStr +=	"C=1,N=axisX;C=1,N=axisY;C=1,N=axisZ;C=1,N=angle;"
		LoadWave/A/B=columnInfoStr/G/O/P=home fullFilePath
		String wName,noteStr = "dataFile="+S_path+S_fileName
		for (i=0;i<V_flag;i+=1)
			wName = StringFromList(i,S_waveNames)
			Note/K $wName
			Note $wName, noteStr
		endfor
	endif

	Wave XX=XX, YY=YY, ZZ=ZZ
	Wave Hx=Hx,Hy=Hy,Hz=Hz, Kx=Kx,Ky=Ky,Kz=Kz, Lx=Lx,Ly=Ly,Lz=Lz
	Wave axisX=axisX,axisY=axisY,axisZ=axisZ,angle=angle

	Variable N=numpnts(XX)
	for (i=0;i<N;i+=1)
		if (Hx[i]==0 && Hy[i]==0 && Hz[i]==0 && Kx[i]==0 && Ky[i]==0 && Kz[i]==0 && Lx[i]==0 && Ly[i]==0 && Lz[i]==0)
			XX[i]=NaN ; YY[i]=NaN; ZZ[i]=NaN
			Hx[i]=NaN ; Hy[i]=NaN ; Hz[i]=NaN
			Kx[i]=NaN ; Ky[i]=NaN ; Kz[i]=NaN
			Lx[i]=NaN ; Ly[i]=NaN ; Lz[i]=NaN
			axisX[i]=NaN ; axisY[i]=NaN ; axisZ[i]=NaN ; angle[i] = NaN
		endif
	endfor
	Variable dX, dY, dZ
	dX = findDelta(XX)
	dY = findDelta(YY)
	dZ = findDelta(ZZ)

	Variable Xmin, Xmax, Ymin,Ymax, Zmin, Zmax
	WaveStats/Q XX
	Xmin = V_min  ;  Xmax = V_max
	printf "range of X is [%g, %g], delta=%g\r",Xmin,Xmax,dX
	WaveStats/Q YY
	Ymin = V_min  ;  Ymax = V_max
	printf "range of Y is [%g, %g], delta=%g\r",Ymin,Ymax,dY
	WaveStats/Q ZZ
	Zmin = V_min  ;  Zmax = V_max
	printf "range of Z is [%g, %g], delta=%g\r",Zmin,Zmax,dZ

	Variable Nx,Ny,Nz
	Nx = round((Xmax-Xmin)/dX + 1)
	Ny = round((Ymax-Ymin)/dY + 1)
	Nz = round((Zmax-Zmin)/dZ + 1)
	printf "(%d, %d, %d)\r",Nx,Ny,Nz
	Make/N=(Nx,Ny,Nz)/O XYZindex
	Make/N=(N,3)/O indicies					// contains the [ix][iy][iz] into XYZindex for each point
	String micron = StrVarOrDefault("root:Packages:Grains:micron","micron")
	SetScale/I x Xmin,Xmax,micron, XYZindex
	SetScale/I y Ymin,Ymax,micron, XYZindex
	SetScale/I z Zmin,Zmax,micron, XYZindex
	if (readin)
		Note/K XYZindex
		Note XYZindex,noteStr
		Note/K indicies
		Note indicies,noteStr
	endif
	indicies = NaN

	Variable ix,iy,iz
	for (i=0;i<N;i+=1)
		ix = round( (XX[i]-DimOffset(XYZindex,0))/DimDelta(XYZindex,0) )
		iy = round( (YY[i]-DimOffset(XYZindex,1))/DimDelta(XYZindex,1) )
		iz = round( (ZZ[i]-DimOffset(XYZindex,2))/DimDelta(XYZindex,2) )
		if (numtype(ix+iy+iz)==0)
			XYZindex[ix][iy][iz] = i
			indicies[i][0] = ix
			indicies[i][1] = iy
			indicies[i][2] = iz
		endif
	endfor

//	Add a section to combine the (Hx,Hy,Hz)(Kx,Ky,Kz)(Lx,Ly,Lz) into one wave OrientMat[N][3][3]
//	Note that (Hx,Hy,Hz) is normalized vector
	wName = UniqueName("vec",1,0)  ;  Make/N=3 $wName  ;  Wave vec=$wName
//	Make/N=3/O vec
	Make/N=(N,3,3)/O OrientMat
	if (readin)
		Note/K OrientMat
		Note OrientMat,noteStr
	endif
	for (i=0;i<N;i+=1)
		vec[0]=Hx[i] ; vec[1]=Hy[i] ; vec[2]=Hz[i]	// doing this only makes sense for cubic
		normalize(vec)
		OrientMat[i][0][] = vec[r]
		vec[0]=Kx[i] ; vec[1]=Ky[i] ; vec[2]=Kz[i]
		normalize(vec)
		OrientMat[i][1][] = vec[r]
		vec[0]=Lx[i] ; vec[1]=Ly[i] ; vec[2]=Lz[i]
		normalize(vec)
		OrientMat[i][2][] = vec[r]
	endfor

	Make/N=(N,3)/O axisXYZ							// make official Rodriques vectors for each voxel
	for (i=0;i<N;i+=1)
		vec[0] = axisX[i]  ;  vec[1] = axisY[i]  ;  vec[2] = axisZ[i]
		normalize(vec)
		vec *= tan(angle[i]/2*PI/180)
		axisXYZ[i][] = vec[q]
	endfor
	//	axisXYZ[][] = numtype(axisXYZ[p][q]) ? 0 : axisXYZ[p][q]
	make3dColorTable(axisXYZ,1,"axisXYZ_Colors")	// make a gizmo color wave to go with each voxel

	KillWaves/Z vec
	SetDataFolder "::"
//	SetDataFolder fldrSav
End
Static Function findDelta(wav)
	Wave wav

	Variable N=numpnts(wav)
	String wName=UniqueName("tempDeltaWave",1,0 )
	Make/N=(N) $wName
	Wave dxw=$wName
	dxw = wav
	Sort dxw,dxw
	dxw = dxw[p+1]-dxw[p]				// a list of delta's

	WaveStats/Q dxw
	Variable thresh = V_max*1e-5			// set a minimum delta to consider
	dxw = (dxw[p]<thresh)?NaN : dxw[p]	// mark for trimimng all delta too close to zero
	Sort dxw,dxw
	WaveStats/Q dxw
	Redimension/N=(V_npnts) dxw			// with all NaN's at end, trim them off

	Variable delta= V_avg
	if ((V_max-V_min) > thresh)
		delta = findMode(dxw,thresh)
		delta = (numtype(delta) || (delta<thresh)) ? V_avg : delta
		String str
		sprintf str,"for wave '%s' give the delta, min=%g, max=%g, avg=%g",NameOfWave(wav),V_min,V_max,V_avg
		Prompt delta,str
		DoPrompt "input delta",delta
	endif
	KillWaves/Z dxw
	return delta
End
Static Function findMode(wav,thresh)
	Wave wav
	Variable thresh
	thresh = max(0,thresh)
	String wName=UniqueName("tempModeWave",1,0 )
	Duplicate wav $wName
	Wave cnt=$wName
	wName=UniqueName("tempModeWave",1,0 )
	Duplicate wav $wName
	Wave val=$wName
	val = NaN
	cnt = 0
	Variable i, N=numpnts(wav), j, m=0
	for(i=0;i<N;i+=1)
		for (j=0;j<m;j+=1)		// is wav[i] in val[] ?
			if (abs(val[j]-wav[i])<thresh)
				break
			endif
		endfor
		val[j] = wav[i]
		cnt[j] += 1
		m = max(j+1,m)
	endfor
	WaveStats/Q cnt
	Variable mode=val[V_maxloc]
	mode = roundSignificant(mode,round(log(mode/thresh)))	// round val to N significant figures
	KillWaves/Z cnt,val
	return mode
End



Function MakeAllGrainStatistics(grains)
	Wave grains						// 3d wave with all of the grain numbers (e.g. gNeighbor)
	if (!WaveExists(grains) || WaveDims(grains)!=3)
		String grainName
		Prompt grainName,"3d wave with all grains",popup,ListMultiGrainMats()
		DoPrompt "grain array",grainName
		if (V_flag)
			return 1
		endif
		Wave grains = $grainName
	endif
	if (!WaveExists(grains) || WaveDims(grains)!=3)
		return 1
	endif
	String gName = GetWavesDataFolder(grains,2)

	Wave indicies=:raw:indicies
	Wave XYZindex=:raw:XYZindex
	Wave axisX=:raw:axisX, axisY=:raw:axisY, axisZ=:raw:axisZ
	Wave angle=:raw:angle
	Wave gVol=$(gName+"_Vol")	// number of voxels in the grain
	Variable N = numpnts(gVol)
	Make/N=(N)/O $(gName+"_COM")	;	Wave grain_COM=$(gName+"_COM")	// single index to the center of mass
	Make/N=(N)/O $(gName+"_xc")	;	Wave grain_xc=$(gName+"_xc")		// scaled postiion of the COM
	Make/N=(N)/O $(gName+"_yc")	;	Wave grain_yc=$(gName+"_yc")
	Make/N=(N)/O $(gName+"_zc")	;	Wave grain_zc=$(gName+"_zc")
	Make/N=(N)/O $(gName+"_ixlo")	;	Wave grain_ixlo=$(gName+"_ixlo")	// total extent of grain in [ix][iy][iz] indicies
	Make/N=(N)/O $(gName+"_ixhi")	;	Wave grain_ixhi=$(gName+"_ixhi")
	Make/N=(N)/O $(gName+"_iylo")	;	Wave grain_iylo=$(gName+"_iylo")
	Make/N=(N)/O $(gName+"_iyhi")	;	Wave grain_iyhi=$(gName+"_iyhi")
	Make/N=(N)/O $(gName+"_izlo")	;	Wave grain_izlo=$(gName+"_izlo")
	Make/N=(N)/O $(gName+"_izhi")	;	Wave grain_izhi=$(gName+"_izhi")
	Make/N=(N,3)/O $(gName+"_Rodrigues")	;	Wave Rodrigues=$(gName+"_Rodrigues")

	grain_ixlo = 1e10
	grain_ixhi = -1e10
	grain_iylo = 1e10
	grain_iyhi = -1e10
	grain_izlo = 1e10
	grain_izhi = -1e10
	String wName = UniqueName("vec",1,0)
	Make/N=3 $wName
	Wave vec=$wName
	Variable i,j, ix,iy,iz
	for (i=0;i<N;i+=1)
		wName = getGrainCoords(grains,i)		// create a wave containing the [ix][iy][iz] coords of each point in grain
		if (exists(wName)!=1)
			continue
		endif
		Wave xyz = $wName
		j = grainCOMsingle(xyz)
		grain_COM[i] = j
		ix = indicies[j][0]
		iy = indicies[j][1]
		iz = indicies[j][2]
		grain_xc[i] = DimDelta(XYZindex,0)*ix + DimOffset(XYZindex,0)
		grain_yc[i] = DimDelta(XYZindex,1)*iy + DimOffset(XYZindex,1)
		grain_zc[i] = DimDelta(XYZindex,2)*iz + DimOffset(XYZindex,2)
		vec[0] = axisX[j]  ;  vec[1] = axisY[j]  ;  vec[2] = axisZ[j]
		normalize(vec)
		vec *= tan(angle[j]/2*PI/180)
		Rodrigues[i][] = vec[q]					// these are now the Rodrigues vectors
		// now find the extent of each grain in [][][] index space
		for (j=0;j<gVol[i];j+=1)
			grain_ixlo[i] = min(xyz[j][0],grain_ixlo[i])
			grain_ixhi[i] = max(xyz[j][0],grain_ixhi[i])
			grain_iylo[i] = min(xyz[j][1],grain_iylo[i])
			grain_iyhi[i] = max(xyz[j][1],grain_iyhi[i])
			grain_izlo[i] = min(xyz[j][2],grain_izlo[i])
			grain_izhi[i] = max(xyz[j][2],grain_izhi[i])
		endfor
		KillWaves/Z xyz
	endfor
	KillWaves/Z vec
	return 0
End



Proc from_data_file()
	Abort "cannot run this"

	Hx Hy Hz	X	 H
	Kx Ky Kz  *	Y = 	 K
	Lx Ly Lz	Z	 L

	H = (Hx,Hy,Hz)
	K = (Kx,Ky,Kz)
	L = (Lx,Ly,Lz)

	from Al_Anneal350C_Orients.dat
		XX			YY			ZZ			Hx			Hy			Hz			Kx			Ky			Kz			Lx			Ly			Lz		axisX		axisY		axisZ		angle
	0.0000		9.0000		9.0000		0.0000		0.0000		0.0000		0.0000		0.0000		0.0000		0.0000		0.0000		0.0000		0.0000		0.0000		0.0000		0.0000
	0.0000		9.0000		10.0000	0.0000		0.0000		0.0000		0.0000		0.0000		0.0000		0.0000		0.0000		0.0000		0.0000		0.0000		0.0000		0.0000
	0.0000		9.0000		11.0000	0.0000		0.0000		0.0000		0.0000		0.0000		0.0000		0.0000		0.0000		0.0000		0.0000		0.0000		0.0000		0.0000
	0.0000		9.0000		12.0000	0.8805		0.4518		0.1437		-0.4316	0.8893		-0.1516	-0.1963	0.0714		0.9779		-0.2230	0.3400		0.8834		29.0933
	0.0000		9.0000		13.0000	0.8806		0.4515		0.1438		-0.4312	0.8894		-0.1518	-0.1964	0.0717		0.9779		-0.2235	0.3402		0.8827		29.0798
	0.0000		9.0000		14.0000	0.8814		0.4498		0.1443		-0.4294	0.8902		-0.1523	-0.1970	0.0723		0.9778		-0.2245	0.3412		0.8793		28.9965
	0.0000		9.0000		15.0000	0.8904		0.4289		0.1522		-0.4494	0.8816		0.1445		-0.0722	-0.1971	0.9777		0.3417		0.2244		0.8783		28.9715
	0.0000		9.0000		16.0000	0.8904		0.4290		0.1523		-0.4495	0.8815		0.1447		-0.0722	-0.1973	0.9777		0.3420		0.2245		0.8785		28.9823
	0.0000		9.0000		17.0000	0.8902		0.4293		0.1523		-0.4498	0.8813		0.1447		-0.0721	-0.1973	0.9777		0.3421		0.2244		0.8791		28.9997
EndMacro




//	=========================================================================
//	=========================================================================
//	=========================================================================

//	here is the way tog get rotation axis and angle.
//	From rotation matrices A, B, you get C=A*B^(-1)
//	The angle = acos((C11+C22+C33-1)/2)
			// cos(angle) = { Trace(A*Binverse)-1 } / 2
//	axis direction (x,y,z): x=C23-C32, y=C31-C13, z=C12-C21
//	Make sure to normalize A and B first.
//		Wenge


// make a wave that has the color for each of the grains.  The color comes from the orientation of the COM 
// voxel for each grain.  The color is the scaled RGB using the Rodrigues vectors
Function makeGrainsColorTable(grains)
	Wave grains						// 3d wave with all of the grain numbers (e.g. gNeighbor)
	if (!WaveExists(grains) || WaveDims(grains)!=3)
		String grainName
		Prompt grainName,"3d wave with all grains",popup,ListMultiGrainMats()
		DoPrompt "grain array for making color table",grainName
		if (V_flag)
			return 1
		endif
		Wave grains = $grainName
	endif
	if (!WaveExists(grains) || WaveDims(grains)!=3)
		DoAlert 0, "in makeGrainsColorTable(), the wave 'grains[][][]' does not exist or has wrong dimension"
		return 1
	endif
	String gName = GetWavesDataFolder(grains,2)
	Wave Rodrigues=$(gName+"_Rodrigues")	// Rodrigues vector

	if (0)			// this section scales the colors to range over all grain angles
		//	make3dColorTable(Rodrigues,65535,gName+"_Colors")
		make3dColorTable(Rodrigues,1,gName+"_Colors")
	else				// this section sets the grain colors to exactly match the individual voxel colors
		Variable N=DimSize(Rodrigues,0)
		Wave com=$(gName+"_COM")					// wave with single index to com of each grain's COM
		String wName = GetWavesDataFolder(grains,1)+"raw:axisXYZ_Colors"
		Wave axisXYZ_Colors=$wName	// wave with single index to color of each voxel
		wName = gName+"_Colors"
		Make/O/N=(N,4) $wName
		Wave colors = $wName
		colors[][] = axisXYZ_Colors[com[p]][q]
		colors[][3] = numtype(colors[p][1]+colors[p][2]+colors[p][3]) ? 0 : colors[p][3]
	endif
End
Function make3dColorTable(vecList,maxRGB,outName)
	Wave vecList						// list of 3-vectors (i.e. points in Rodrigues space)
	Variable maxRGB					// max value for an RGB, usually either 1 or 65535
	String outName						// name of output wave, usually something like "vecList_Colors"
	if (!WaveExists(vecList) || WaveDims(vecList)!=2 || DimSize(vecList,1)!=3 || !(maxRGB>0) || strlen(outName)<1)
		DoAlert 0, "in make3dColorTable(), the wave 'vecList[][3]' does not exist or has wrong dimension, or other invalid input"
		return 1
	endif

	Variable N = DimSize(vecList,0)
	Variable xlo=1e10,xhi=-1e10,ylo=1e10,yhi=-1e10,zlo=1e10,zhi=-1e10
	Variable i
	for (i=0;i<N;i+=1)
		if (numtype(vecList[i][1]+vecList[i][2]+vecList[i][3]))
			continue
		endif
		xlo = min(xlo,vecList[i][0])
		xhi = max(xhi,vecList[i][0])
		ylo = min(ylo,vecList[i][1])
		yhi = max(yhi,vecList[i][1])
		zlo = min(zlo,vecList[i][2])
		zhi = max(zhi,vecList[i][2])
	endfor

	Variable dx=xhi-xlo, dy=yhi-ylo, dz=zhi-zlo	// find range of axes, needed to scale the color table
	//	printf "x0=%g, dx=%g,    y0=%g, dy=%g,    z0=%g, dz=%g\r",xlo,dx,ylo,dy,zlo,dz
	if (abs(maxRGB-65535)<100)
		Make/O/W/U/N=(N,4) $outName
	else
		Make/O/N=(N,4) $outName
	endif
	Wave colors = $outName
	colors = 0										// default is black (maxRGB is white)
	colors[][3] = 1									// default alpha is 1
	for (i=0;i<N;i+=1)
		colors[i][2] = maxRGB*(vecList[i][0]-xlo)/dx
		colors[i][1] = maxRGB*(vecList[i][1]-ylo)/dy
		colors[i][0] = maxRGB*(vecList[i][2]-zlo)/dz
		colors[i][3] = numtype(colors[i][1]+colors[i][2]+colors[i][3]) ? 0 : colors[i][3]
	endfor
End
///
//Function makeGrainsColorTable(grains)
//	Wave grains						// 3d wave with all of the grain numbers (e.g. gNeighbor)
//	if (!WaveExists(grains) || WaveDims(grains)!=3)
//		String grainName
//		Prompt grainName,"3d wave with all grains",popup,ListMultiGrainMats()
//		DoPrompt "grain array for making color table",grainName
//		if (V_flag)
//			return 1
//		endif
//		Wave grains = $grainName
//	endif
//	if (!WaveExists(grains) || WaveDims(grains)!=3)
//		DoAlert 0, "in makeGrainsColorTable(), the wave 'grains[][][]' does not exist or has wrong dimension"
//		return 1
//	endif
//
//	String gName = GetWavesDataFolder(grains,2)
//	Wave grain_COM=$(gName+"_COM")	// single index to the center of mass
//	Variable Ngrains = numpnts(grain_COM)
//	String fldr = GetWavesDataFolder(grains,1)
//	Wave Rodrigues=$(gName+"_Rodrigues")	// Rodrigues vector
//	Variable xlo=1e10,xhi=-1e10,ylo=1e10,yhi=-1e10,zlo=1e10,zhi=-1e10
//	Variable i
//	for (i=0;i<Ngrains;i+=1)
//		xlo = min(xlo,Rodrigues[i][0])
//		xhi = max(xhi,Rodrigues[i][0])
//		ylo = min(ylo,Rodrigues[i][1])
//		yhi = max(yhi,Rodrigues[i][1])
//		zlo = min(zlo,Rodrigues[i][2])
//		zhi = max(zhi,Rodrigues[i][2])
//	endfor
//
//	// find range of axes, needed to scale the color table
//	Variable dx=xhi-xlo, dy=yhi-ylo, dz=zhi-zlo
//	//	printf "x0=%g, dx=%g,    y0=%g, dy=%g,    z0=%g, dz=%g\r",xlo,dx,ylo,dy,zlo,dz
//	Make/O/W/U/N=(Ngrains,3) $(gName+"_Colors")
//	Wave grainColors = $(gName+"_Colors")
//	Variable maxColor=0, minColor=3*65535, colorSum
//	grainColors = 0						// default is black (65535 is whire)
//	for (i=0;i<Ngrains;i+=1)
//		grainColors[i][2] = 65535*(Rodrigues[i][0]-xlo)/dx	// color is scaled by Rodrigues space
//		grainColors[i][1] = 65535*(Rodrigues[i][1]-ylo)/dy
//		grainColors[i][0] = 65535*(Rodrigues[i][2]-zlo)/dz
//		colorSum = grainColors[i+1][0]+grainColors[i+1][1]+grainColors[i+1][2]
//		maxColor = max(maxColor,colorSum)
//		minColor = min(minColor,colorSum)
//	endfor
//	maxColor = round(maxColor/3)
//	minColor = round(minColor/3)
//	printf "color sum in range [%g,%g]\r",maxColor,minColor
//End




Function NumOfGrainsBiggerThan(grains,minSize)
	Wave grains				// 3d wave with all of the grain numbers (e.g. gNeighbor)
	Variable minSize			// mininum size to include
	if (!WaveExists(grains) || WaveDims(grains)!=3 || numtype(minSize) || minSize<1)
		minSize = (minSize>=1) ? minSize : 20
		Prompt minSize, "min no. of voxels in a grain"
		String grainNoWave=NameOfWave(grains)
		Prompt grainNoWave,"3d wave with all grains",popup,ListMultiGrainMats()
		DoPrompt "grains to make" grainNoWave,minSize
		if (V_flag)
			return 1
		endif
		Wave grains = $grainNoWave
		if (ItemsInList(GetRTStackInfo(0))<2)
			printf "     NumOfGrainsBiggerThan(%s,%g)\r",NameOfWave(grains),minSize
		endif
	endif
	if (numtype(minSize) || minSize<1)
		return 0
	endif
	if (!WaveExists(grains) || WaveDims(grains)!=3)
		DoAlert 0, "'grains' wave does not exist"
		return 0
	endif

	String gName = GetWavesDataFolder(grains,2)
	Wave gVol=$(gName+"_Vol")

	Variable N=BinarySearch(gvol,minSize)+1
	if (ItemsInList(GetRTStackInfo(0))<2)
		printf "from array '%s' there are %d grains of volume %d or bigger\r",NameOfWave(grains),N,minSize
	endif
	return N
End


// the contents of the supplied wave is index (i,j,k) not microns
Function/S getGrainCoords(grains,grainNum)	// create a wave containing the [ix][iy][iz] coords of each point in grain
	Wave grains								// 3d wave with all of the grain numbers (e.g. gNeighbor)
	Variable grainNum							// grain number

	if (!WaveExists(grains) || numtype(grainNum) || grainNum<0)
		String grainName=NameOfWave(grains)
		Prompt grainName,"3d wave with all grains",popup,ListMultiGrainMats()
		grainNum = numtype(grainNum)==2 ? 0 : max(grainNum,0)
		Prompt grainNum, "grain number"
		DoPrompt "choose grain",grainName,grainNum
		if (V_flag)
			return ""
		endif
		Wave grains = $grainName
	endif
	if (!WaveExists(grains))
		return ""
	endif
	WaveStats/Q grains
	Variable maxGrain = V_max
	if (numtype(grainNum) || grainNum<0 || grainNum>maxGrain)
		return ""
	endif

	Variable Nx,Ny,Nz
	Nx=DimSize(grains,0)
	Ny=DimSize(grains,1)
	Nz=DimSize(grains,2)

	String str="xyz"+num2istr(grainNum)
	Make/N=(100,3)/O $str
	Wave xyz=$str

	Variable ix,iy,iz, N=0
	for (ix=0;ix<Nx;ix+=1)
		for (iy=0;iy<Ny;iy+=1)
			for (iz=0;iz<Nz;iz+=1)
				if (grains[ix][iy][iz]==grainNum)
					if (N>=DimSize(xyz,0))
						Redimension/N=(DimSize(xyz,0)+100,3) xyz		// extend xyz[][3]
					endif
					xyz[N][0] = ix
					xyz[N][1] = iy
					xyz[N][2] = iz
					N += 1
				endif
			endfor
		endfor
	endfor
	Redimension/N=(N,3) xyz		// extend xyz[][3]
	if (ItemsInList(GetRTStackInfo(0))<2)
		printf "   for grain %d, made a list of the voxel indicies in '%s'\r",grainNum,GetWavesDataFolder(xyz,2)
	endif
	return GetWavesDataFolder(xyz,2)
End


Function grainCOMsingle(xyz)	// returns the single index to the COM of the grain
	Wave xyz					// xyz[N][3], a list of ix,iy,iz for each point in the grain
	if (!WaveExists(xyz))
		String xyzName
		Prompt xyzName, "wave with gain points 'xyz*'",popup,WaveList("xyz*",";", "DIMS:2,MAXCOLS:3,MINCOLS:3")
		DoPrompt "grain list",xyzName
		if (V_flag)
			return NaN
		endif
		Wave xyz=$xyzName
	endif
	if (!WaveExists(xyz))
		return NaN
	endif

	Wave XYZindex=:raw:XYZindex		// 3d array with single index for each point
	Variable N=DimSize(xyz,0)

	Variable ix,iy,iz			// results, must actually land on a point in xyz
	Make/N=(N)/O temp_find_avg
	temp_find_avg = xyz[p][0]
	ix = sum(temp_find_avg)/N
	temp_find_avg = xyz[p][1]
	iy = sum(temp_find_avg)/N
	temp_find_avg = xyz[p][2]
	iz = sum(temp_find_avg)/N
	KillWaves/Z temp_find_avg
	// now ix,iy,iz are the average values, next make sure that we have a real voxel

	Variable i,d2,min2=1e100,imin
	for (i=0;i<N;i+=1)
		d2 = (xyz[i][0]-ix)^2 + (xyz[i][1]-iy)^2 + (xyz[i][2]-iz)^2
		if (d2 < min2)
			imin = i
			min2 = d2
		endif
	endfor
	ix = xyz[imin][0]
	iy = xyz[imin][1]
	iz = xyz[imin][2]
	Variable index = XYZindex[ix][iy][iz]
	if (ItemsInList(GetRTStackInfo(0))<2)
		printf "   Center of Mass for grain list %s is point %d\r",NameOfWave(xyz),index
	endif
	return index
End


Function angleBetweenGrains(grains,i1,i2)
	Wave grains						// 3d wave with all of the grain numbers (e.g. gNeighbor)
	Variable i1,i2						// grain numbers

	if (!WaveExists(grains) || numtype(i1+i2) || i1<0 || i2<0)
		String grainName=NameOfWave(grains)
		Prompt grainName,"3d wave with all grains",popup,ListMultiGrainMats()
		i1 = numtype(i1)==2 ? 0 : max(i1,0)
		i2 = numtype(i2)==2 ? 1 : max(i2,0)
		Prompt i1, "first grain number"
		Prompt i2, "second grain number"
		DoPrompt "choose grains",grainName,i1,i2
		if (V_flag)
			return NaN
		endif
		Wave grains = $grainName
	endif
	if (!WaveExists(grains))
		return NaN
	endif
	WaveStats/Q grains
	Variable maxGrain = V_max
	if (numtype(i1+i2) || i1<0 || i2<0 || i1>maxGrain || i2>maxGrain)
		if (ItemsInList(GetRTStackInfo(0))<2)
			DoAlert 0,"invalid inputs"
		endif
		return NaN
	endif

	Wave g1 = $getGrainCoords(grains,i1)		// list of grain coords for each grain
	Wave g2 = $getGrainCoords(grains,i2)
	Variable j1 = grainCOMsingle(g1)
	Variable j2 = grainCOMsingle(g2)
	Variable angle = angleBetweenPointsSymReduced(j1,j2)
	if (ItemsInList(GetRTStackInfo(0))<2)
		String degree = StrVarOrDefault("root:Packages:Grains:degree","¡")
		printf "   angle between grains %d and %d is %g%s\r",i1,i2,angle,degree
	endif
	KillWaves/Z g1,g2
	return angle
End


Function anglefromMatrix(index)			// returns total rotaion angle (degree)
	// function to determine total rotation angle for one rotation matrix
	Variable index						// index to point
	Wave OrientMat=:raw:OrientMat
	Variable trace
	trace = OrientMat[index][0][0]+OrientMat[index][1][1]+OrientMat[index][2][2]
	Variable cosine = (trace-1) / 2
	cosine = min(max(-1,cosine),1)	// ensure valid range for cosine
	return acos(cosine)*180/PI
End


Function/S ListMultiGrainMats()	// creates a list of the matricies that contain all of the grains
	String list=WaveList("g*",";","DIMS:3")
	// next remove all items ending with "_number" from list
	String item,out=""
	Variable j,i=0
	for (i=0;i<ItemsInList(list);i+=1)
		item = StringFromList(i,list)
		j = strsearch(item,"_",Inf,1)
		if ((j<0 || numtype(str2num(item[j+1,Inf]))) && strlen(note($item)))
			out = AddListItem(item,out)
		endif
	endfor
	return out
End


Function locateValueInWave(wav,val)	// returns the point number i, where wav[i]==val
	Wave wav
	Variable val

	Variable i,N=numpnts(wav)
	for (i=0;i<N;i+=1)
		if (wav[i]==val)
			return i
		endif
	endfor
	return NaN
End

Function AskAngleBetweenPointsSymReduced(index1,index2)	// returns angle between 1 and 2 in degrees
	Variable index1,index2			// index to points 1 and 2
	if (numtype(index1+index2) || index1<0 || index2<0)
		index1 = (numtype(index1) || index1<0) ? 0 : index1
		index2 = (numtype(index2) || index2<0) ? index1+1 : index2
		Prompt index1, "index to first voxel"
		Prompt index2, "index to second voxel"
		DoPrompt "voxel numbers",index1,index2
		if (V_flag)
			return NaN
		endif
	endif
	angleBetweenPointsSymReduced(index1,index2)
End
// function to determine the angle between two rotation matricies
// this will return the smallest angle after a check of all the symmetry operations
Function angleBetweenPointsSymReduced(index1,index2)	// returns angle between 1 and 2 in degrees
	Variable index1,index2			// index to points 1 and 2
	String wName
	wName = UniqueName("OrientMatrix",1,0)  ;  Make/N=(3,3)/D/O $wName  ;  Wave A = $wName
	wName = UniqueName("OrientMatrix",1,0)  ;  Make/N=(3,3)/D/O $wName  ;  Wave B = $wName
	wName = UniqueName("SymMatrix",1,0)  ;  Make/N=(3,3)/D/O $wName  ;  Wave symOp = $wName
	wName = UniqueName("matrix1D",1,0)  ;  Make/N=(1)/D/O $wName  ;  Wave cosine11 = $wName
	Wave OrientMat=:raw:OrientMat
	A = OrientMat[index1][p][q]
	B = OrientMat[index2][p][q]

	Wave symOps=$StrVarOrDefault("root:Packages:Grains:SymmetryOpsPath","")	// wave with symmetry ops
	Variable N = NumberByKey("Nproper",note(symOps),"=")	// number of symmetyry operations to check
	N = numtype(N) ? DimSize(symOps,0) : N
	Variable i, cosMax = -2
	Variable iBest
	for (i=0;i<N;i+=1)
		symOp = symOps[i][p][q]				// symmetry operation
		//	rotation matrix from A to B is just (B * A^-1 * S) with a symmetry op
		MatrixOp/O cosine11 = (Trace(B x Inv(A) x symOp)-1)/2
		iBest = cosMax<cosine11[0] ? i : iBest
		cosMax = max(cosMax,cosine11[0])
	endfor
	cosMax = min(max(-1,cosMax),1)			// ensure cosine in range [-1,1]
	Variable angle = acos(cosMax)*180/PI

	String stackStr = GetRTStackInfo(0)
	if (ItemsInList(stackStr)<2 || strsearch(stackStr,"AskAngleBetweenPointsSymReduced;", 0)==0)
		if (numtype(angle))
			print "one of the matrices is invalid"
		else
			printf "taking the transpose of both matrices before doing anthing, maximize [%s x S x Inv(%s)]\r",NameOfWave(B),NameOfWave(A)
			print3x3(A) ;  print3x3(B)
			symOp = symOps[iBest][p][q] ; print3x3(symOp)
			String degree = StrVarOrDefault("root:Packages:Grains:degree","¡")
			printf "total angle between %d and %d minimized with symmetry op %d is %g%s\r",index1,index2,iBest,angle,degree
		endif
	endif
	KillWaves/Z A,B,symOp, cosine11
	return angle
End
//
Static Function BestSymOp(A,B,Si)				// return best cubic symmetry op for going from A to B
	Wave A,B									// two orientation matricies
	Wave Si										// best cubic symmetry op (returned)
	String wName
	wName = UniqueName("SymMatrix",1,0)	;  Make/N=(3,3)/D/O $wName  ;  Wave symOp = $wName
	wName = UniqueName("matrix1D",1,0)		;  Make/N=(1)/D/O $wName	;  Wave trace11 = $wName

	Wave symOps=$StrVarOrDefault("root:Packages:Grains:SymmetryOpsPath","")	// wave with symmetry ops
	Variable N = NumberByKey("Nproper",note(symOps),"=")	// number of symmetyry operations to check
	N = numtype(N) ? DimSize(symOps,0) : N
	Variable i, traceMax = -4
	for (i=0;i<N;i+=1)
		symOp = symOps[i][p][q]					// symmetry operation
		MatrixOp/O trace11 = Trace(B x Inv(A) x symOp)	// Tr[ B(A^-1)(S^i) ]
		if (trace11[0] > traceMax)					// save value for optimum rotation matrix
			traceMax = trace11[0]
			Si = symOp
		endif
	endfor
	KillWaves/Z symOp, trace11
End
//
// function to determine the angle between two rotation matricies
//	rotation matrix from a to b is just (Rb * Ra^-1), and don't forget that Ra^-1 = Rtranspose
Function angleBetweenMatricies(mat1,mat2)	// returns angle between 1 and 2 in degrees
	Wave mat1,mat2
	MatrixOp/O temp_angle_11_ = Trace(mat2 x mat1^t)
	Variable cosine = (temp_angle_11_[0][0]-1) / 2
	KillWaves/Z temp_angle_11_
	cosine = min(max(-1,cosine),1)			// ensure cosine in range [-1,1]
	return acos(cosine)*180/PI
End



// returns rot, rotation matrix from A to B, which is symmetry reduced to smallest total angle
Static Function symReducedRotation(A,B,rot)
	Wave A,B									// two orientation matricies
	Wave rot									// resulting symmetry reduced rotation matrix
	String wName
	wName = UniqueName("SymMatrix",1,0)	;  Make/N=(3,3)/D/O $wName  ;  Wave symOp = $wName
	wName = UniqueName("rotMat",1,0)		;  Make/N=(3,3)/D/O $wName  ;  Wave BASi_ = $wName
	Wave symOps=$StrVarOrDefault("root:Packages:Grains:SymmetryOpsPath","")	// wave with symmetry ops
	Variable N = NumberByKey("Nproper",note(symOps),"=")	// number of symmetyry operations to check
	N = numtype(N) ? DimSize(symOps,0) : N
	Variable i, trace, traceMax = -4
	for (i=0;i<N;i+=1)
		symOp = symOps[i][p][q]				// symmetry operation
		MatrixOp/O BASi_ = B x Inv(A) x symOp	// B(A^-1)(S^i)
		trace = MatrixTrace(BASi_)
		if (trace > traceMax)					// save value for optimum rotation matrix
			traceMax = trace
			rot = BASi_
		endif
	endfor
	Variable cosine = (traceMax-1)/2
	cosine = max(min(cosine,1),-1)
	Variable angle = acos(cosine)*180/PI
	KillWaves/Z symOp, BASi_
	return angle
End


Function TransformSample2Beamline(vec)	// change vec[3] from sample (xyz) -> beamline (XYZ) coordinates
	Wave vec
	Variable ys = vec[1], zs=vec[2]			// save sample system y and z
	vec[1] = (-ys+zs)/sqrt(2)	
	vec[2] = (-ys-zs)/sqrt(2)
	return 0
End
Function TransformBeamline2Sample(vec)	// change vec[3] from beamline (XYZ) -> sample (xyz) coordinates
	Wave vec
	Variable Ybl = vec[1], Zbl=vec[2]		// save beamline system Y and Z
	vec[1] = (-Ybl-Zbl)/sqrt(2)
	vec[2] = (  Ybl-Zbl)/sqrt(2)
	return 0
End




Function point_in_list(ix,iy,iz,list3,N)			// returns true if (ix,iy,iz) is a valid point in list3
	Variable ix,iy,iz
	Wave list3									// list3[N][3], list of triplets to check
	Variable N									// number of points to check in list3
	Variable i
	for (i=0;i<N;i+=1)
		if (list3[i][0]==ix)
			if (list3[i][1]==iy)
				if (list3[i][2]==iz)			// ix, iy, & iz  all match, point is in list3
					return 1
				endif
			endif
		endif
	endfor
	return 0
End




//	=========================================================================
//	=========================================================================
//	=========================================================================




Function rot3d2xyz(rot3d)				// turn 3-d matrix in to separate waves for a scatter plot
	Wave rot3d
	Variable N=numpnts(rot3d)			// total number of scatter points
	Variable nj=DimSize(rot3d,0)
	Variable nk=DimSize(rot3d,1)
	Variable nm=DimSize(rot3d,2)
	Make/N=(N,3)/O xyz				// total number of scatter points
	Variable j,k,m,i=0
	for (m=0;m<nm;m+=1)
		for (k=0;k<nk;k+=1)
			for (j=0;j<nj;j+=1)
				xyz[i][0] = j
				xyz[i][1] = k
				xyz[i][2] = m
				i += 1
			endfor
		endfor
	endfor
End


Function GrainCOM(grainNoW,ci,cj,ck)	// returns the grain Center Of Mass in ci,cj,ck
	Wave grainNoW						// 3d wave of numbered grains
	Variable &ci,&cj,&ck				// returns with the COM (Center of Mass)
	if (!WaveExists(grainNoW))
		Wave grainNoW = $Get_GrainNo_Dialog(grainNoW)
	endif
	ci=0 ; cj=0 ; ck=0
	Variable i,j,k, n=0
	for (k=0;k<DimSize(grainNoW,2);k+=1)
		for (j=0;j<DimSize(grainNoW,1);j+=1)
			for (i=0;i<DimSize(grainNoW,0);i+=1)
				if(grainNoW[i][j][k]>-0.5)
					ci += i
					cj += j
					ck += k
					n += 1
				endif
			endfor
		endfor
	endfor
	ci /= n
	cj /= n
	ck /= n
End
//Function test()
//	Variable i,j,k
//	GrainCOM($"",i,j,k)
//	print i,j,k
//End


//Function EditRefMatrix() : Table
//	if (WinType("ref_matrix_win")==2)
//		DoWindow /F ref_matrix_win
//		return 1
//	endif
//	Edit/W=(127,45,437,170)/K=1 root:ref_matrix
//	Execute "ModifyTable width(Point)=42,format(root:ref_matrix)=3"
//	DoWindow/C ref_matrix_win
//End


Function/S grainRGBlist(i)
	Variable i
	i = mod(abs(i),7)
	switch(i)
		case 0:
			return "1,0,0"
		case 1:
			return "1,1,0"
		case 2:
			return "0,1,0"
		case 3:
			return "0,1,1"
		case 4:
			return "0,0,1"
		case 5:
			return "1,0,1"
		case 6:
			return "0.1,0.1,0.1"
	endswitch
	return ""
End

Function MakeWiderRainbow()
	Variable N=101
	Make/N=(N,3)/O M_colors
	Variable i, frac
	Make/N=3/O color100_, color110_, color010_, color011_, color001_, color101_
	color100_ = {1,0,0}		// red
	color110_ = {1,1,0}		// yellow
	color010_ = {0,1,0}		// green
	color011_ = {0,1,1}		// cyan
	color001_ = {0,0,1}		// blue
	color101_ = {1,0,1}		// magenta
	Variable step = (N-1)/5
	for (i=0;i<N;i+=1)
		if (i<=step)
			frac = i/step
			M_colors[i][] = (1-frac)*color100_[q] + frac*color110_[q]
		elseif (i<=2*step)
			frac = (i-step)/step
			M_colors[i][] = (1-frac)*color110_[q] + frac*color010_[q]
		elseif (i<=3*step)
			frac = (i-2*step)/step
			M_colors[i][] = (1-frac)*color010_[q] + frac*color011_[q]
		elseif (i<=4*step)
			frac = (i-3*step)/step
			M_colors[i][] = (1-frac)*color011_[q] + frac*color001_[q]
		else
			frac = (i-4*step)/step
			M_colors[i][] = (1-frac)*color001_[q] + frac*color101_[q]
		endif
	endfor
	M_colors *= 65535
	KillWaves/Z color100_, color110_, color010_, color011_, color001_, color101_
End

Function BadWave(a)		// returns true if any of the wave elements are NaN or an Inf
	Wave a
	Variable N=numpnts(a)
	Make/N=(N)/O/B BadWaves_temp_
	BadWaves_temp_ = numtype(a)
	Variable bad = sum(BadWaves_temp_,0,N-1)>0
	KillWaves/Z BadWaves_temp_
	return bad
End





//	here is the way to get rotation axis and angle.
//	From rotation matrices A, B, you get C=A*B^(-1)
//	The angle = acos((C11+C22+C33-1)/2)
			// cos(angle) = { Trace(A*Binverse)-1 } / 2
//	axis direction (x,y,z): x=C23-C32, y=C31-C13, z=C12-C21
//	Make sure to normalize A and B first.
Function AngleBetweenMats24andSwaps(orient,ref,Sswap)
	Wave orient									// matrix of all of the orientation matricies
	String Sswap								// name of matrix to get the swap vector, if "" then don't make Sswap
	Wave ref								// reference orientation matrix
	// From the orientation matricies, make swap vectors that will swap a*, b*, c* so that they are all
	// as close together as possible (less than 90¡)

	if (DimSize(orient,0)!=3 || DimSize(orient,1)!=3 || WaveDims(orient)!=2)
		Abort "MakeSwapVectorsReference, orient matrix is illegal size"
	endif
	if (!WaveExists(ref))
		Abort "no refrence orientation matrix given"
	endif

	Make/N=3/O oneaxis__

	Duplicate/O orient orient_Temp_Norm_					// a normalized version of orient (lattice vectors all length 1)
	oneaxis__ = orient[p][0]
	normalize(oneaxis__)
	orient_Temp_Norm_[][0] = oneaxis__[p]
	oneaxis__ = orient[p][1]
	normalize(oneaxis__)
	orient_Temp_Norm_[][1] = oneaxis__[p]
	oneaxis__ = orient[p][2]
	normalize(oneaxis__)
	orient_Temp_Norm_[][2] = oneaxis__[p]

	Duplicate/O ref ref_Temp_Norm_					// a normalized version of ref (lattice vectors all length 1)
	oneaxis__ = ref[p][0]
	normalize(oneaxis__)
	ref_Temp_Norm_[][0] = oneaxis__[p]
	oneaxis__ = ref[p][1]
	normalize(oneaxis__)
	ref_Temp_Norm_[][1] = oneaxis__[p]
	oneaxis__ = ref[p][2]
	normalize(oneaxis__)
	ref_Temp_Norm_[][2] = oneaxis__[p]
	MatrixTranspose ref_Temp_Norm_

	Variable i,j,k,eps
	Variable i0=1,j0=2,k0=3
	Make/N=(3,3)/O rmat__
	Variable cosine
	Variable maxCos=-Inf
	for(k=-3;k<4;k+=1)								// check all 48 of the 90¡ type rotations
		for(j=-3;j<4;j+=1)
			for(i=-3;i<4;i+=1)
				eps = epsilon_ijk(i,j,k)				// one possible swapping
				if (eps==1)							// only 24 right handed swappings are checked
					rmat__ = orient_Temp_Norm_
					SwapXYZ2xyz(rmat__,i,j,k)		// rmat__ is now swapped version of normalized orient
					MatrixMultiply rmat__,ref_Temp_Norm_		// M_product now contains product
					cosine = (MatrixTrace(M_product)-1)/2		// cosine of the rotation angle
					cosine = (cosine>1) ? (2-cosine) : cosine
					if (cosine>maxCos)
						maxCos = cosine
						i0=i  ;  j0=j  ;  k0=k
					endif
				endif
			endfor
		endfor
	endfor

	if (strlen(Sswap)>0)
		Make/N=3/O $Sswap
		Wave swaps = $Sswap
		swaps[0,2]={i0,j0,k0}								// save best swap
	endif
	Killwaves/Z rmat__,M_product
	KillWaves/Z oneaxis__, ref_Temp_Norm_,orient_Temp_Norm_
	return acos(maxCos)*180/PI							// rotation angle in degrees
End




Function SwapXYZ2xyz(mat,x2,y2,z2)		// permute axes so X->x, Y->y, Z->z
	Wave mat				// matrix to permute (this is NOT a rotation)
	Variable x2,y2,z2		// destinations of x,y,z (use ±1, ±2, and ±3 (cannot use 0, since 0==±0)

	Variable eps = epsilon_ijk(x2,y2,z2)
	if (!eps)
		DoAlert 0, "SwapXYZ2xyz x2,y2,z2, must be a permutation of ±1,2, or 3"
	endif

	Make/N=(3,3)/O axses__
	axses__ = mat
	mat[][0] = axses__[p][abs(x2)-1] * sign(x2)
	mat[][1] = axses__[p][abs(y2)-1] * sign(y2)
	mat[][2] = axses__[p][abs(z2)-1] * sign(z2)

	Killwaves/Z axses__
	return eps
End



Function MakeRotationMatFromAB(a,b,rot)	// makes rotation mat between a and b
	Wave a,b							// input matricies
	Wave rot							// rotation matrix

	// makes rmat, which is the  rotation matricies between succesive grains
	// rmat[i] is the rotation needed to go from orient[i-1] to orient[i]

	if (numtype(a[0][0]) || numtype(b[0][0]))
		rot = NaN
		return 1
	endif

	Make/N=(3,3)/O temp__, invA_temp__
	rot = a
	temp__ = (p==q)
	MatrixLLS /O/Z  rot temp__		// temp__ is now a^-1
	invA_temp__ = temp__				// invA_temp__ is now a^-1
	MatrixMultiply b,invA_temp__		// rotation matrix is now in M_product
	Wave M_product=M_product

	Variable det = MatrixDet(M_product)
	if (det>1)
		M_product /= det
	endif
	rot = M_product
	Killwaves/Z temp__,M_product, invA_temp__
	return 0
EndMacro





// ============================================================================
// ==============  this section are matrix math routines from the old PoleFigure.ipf   ===============


Function rotationCosineOfMat(rot)				// returns cos(total rotation angle of matrix) for matrix 'rot'
	Wave rot
	Variable trace = MatrixTrace(rot)			// trace = 1 + 2*cos(theta)
	return (trace-1)/2						// cosine of the rotation angle
End


Function rotationMatAboutAxis(axis,angle,mat)
	Wave axis				// axis about which to rotate
	Variable angle			// angle to rotate (passed as degrees, but convert to radians)
	Wave mat				// desired rotation matrix
	angle *= PI/180

	Make/N=3/O/D xhat_rotate__, yhat_rotate__, zhat_rotate__
	Wave xhat=xhat_rotate__
	Wave yhat=yhat_rotate__
	Wave zhat=zhat_rotate__

	zhat = axis
	if (normalize(zhat) <=0)
		KillWaves/Z xhat_rotate__, yhat_rotate__, zhat_rotate__
		return 1			// error, cannot rotate about a zero length axis
	endif

	Variable i
	i = ( abs(zhat[0])<= abs(zhat[1]) ) ? 0 : 1
	i = ( abs(zhat[2])< abs(zhat[i]) ) ? 2 : i
	xhat = zhat
	xhat[i] = 2							// choose x-z plane
	normalize(xhat)
	Variable dotxz = MatrixDot(xhat,zhat)
	xhat -= dotxz*zhat					// xhat is now perpendicular to zhat
	normalize(xhat)					// xhat is now normalized vector perp to zhat
	crossV4(zhat,xhat,yhat)
	normalize(yhat)					// yhat is now normalized vector perp to zhat and xhat
//	xhat, yhat, zhat is now an orthonormal system with zhat || to axis

//	now rotate around zhat in the xhat -> yhat direction, an angle of angle
	Variable cosa = cos(angle)
	Variable sina = sin(angle)
	Redimension/N=(3,3) mat
	for(i=0;i<3;i+=1)
		mat[][i] = (xhat[i]*cosa+yhat[i]*sina)*xhat[p] + (-xhat[i]*sina + yhat[i]*cosa)*yhat[p] + zhat[i]*zhat[p]
	endfor
	KillWaves/Z xhat_rotate__, yhat_rotate__, zhat_rotate__
End
// cross supplanted by the Cross operation in Igor 5
Static Function crossV4(a,b,c)	// vector cross product
	Wave a,b,c
	c[2] = a[0]*b[1] - a[1]*b[0]
	c[0] = a[1]*b[2] - a[2]*b[1]
	c[1] = a[2]*b[0] - a[0]*b[2]
End


Function rotation(mat,axis,angle)		// rotate mat about x, y, or z by angle
	Wave mat				// matrix to rotate
	Variable axis			// axis to rotate about 0=x, 1=y, 2=z
	Variable angle			// angle to rotate (radians)

	Variable cosa = cos(angle)
	Variable sina = sin(angle)
	Make/N=(3,3)/O mat_rot_
	mat_rot_ = (p==q)

	switch(axis)
		case 0:			// rotate about x-axis
			mat_rot_[1][1] = cosa
			mat_rot_[2][2] = cosa
			mat_rot_[1][2] = -sina
			mat_rot_[2][1] =   sina
			break
		case 1:			// rotate about y-axis
			mat_rot_[0][0] = cosa
			mat_rot_[2][2] = cosa
			mat_rot_[0][2] =   sina
			mat_rot_[2][0] = -sina
			break
		case 2:			// rotate about z-axis
			mat_rot_[0][0] = cosa
			mat_rot_[1][1] = cosa
			mat_rot_[0][1] = -sina
			mat_rot_[1][0] =   sina
			break
		default:
	endswitch
	Killwaves/Z mat_rot_
End



//	here is the way to get rotation axis and angle.
//	From rotation matrices A, B, you get C=A*B^(-1)
//	The angle = acos((C11+C22+C33-1)/2)
			// cos(angle) = { Trace(A*Binverse)-1 } / 2
//	axis direction (x,y,z): x=C23-C32, y=C31-C13, z=C12-C21
//	Make sure to normalize A and B first.
//		Wenge

Function AngleBetweenMats(am,bm)
	Wave am,bm										// the two matricies

	Make/N=(3,3)/O am__, bm__
	Make/N=3/O oneaxis__

	oneaxis__ = am[p][0]
	normalize(oneaxis__)
	am__[][0] = oneaxis__[p]
	oneaxis__ = am[p][1]
	normalize(oneaxis__)
	am__[][1] = oneaxis__[p]
	oneaxis__ = am[p][2]
	normalize(oneaxis__)
	am__[][2] = oneaxis__[p]

	oneaxis__ = bm[p][0]
	normalize(oneaxis__)
	bm__[][0] = oneaxis__[p]
	oneaxis__ = bm[p][1]
	normalize(oneaxis__)
	bm__[][1] = oneaxis__[p]
	oneaxis__ = bm[p][2]
	normalize(oneaxis__)
	bm__[][2] = oneaxis__[p]

	MatrixTranspose bm__
	MatrixMultiply am__,bm__		// M_product now contain product
	KillWaves/Z am__, bm__, oneaxis__
	return rotationAngleOfMat(M_product)
//	return acos((MatrixTrace(M_product)-1)/2)*180/PI
End



Static Function normalize(a)	// normalize a and return the initial magnitude
	Wave a
	if (WaveDims(a)==1)
		Variable norm_a = norm(a)
		if (norm_a==0)
			return 0
		endif
		a /= norm_a
		return norm_a
	elseif(WaveDims(a)==2 && DimSize(a,0)==DimSize(a,1))
		Variable det = MatrixDet(a)
		a /= det
		return det
	endif
End



Function epsilon_ijk(i,j,k)
	Variable i,j,k					// returns +1 for even permutaion, -1 for odd permutaions

	Variable ps = sign(i)*sign(j)*sign(k)
	i = round(abs(i))
	j = round(abs(j))
	k = round(abs(k))

	if (i*j*k ==0)					// no zeros allowed
		return 0
	endif
	if (i>3 || j>3 || k>3)			// must be in range [1,3]
		return 0
	endif
	if (i==j || j==k || k==i)		// do duplicates
		return 0
	endif
	Variable perm = (i==1) + 2*(j==2)  + 4*(k==3)
	perm = ((perm==0) || (perm==7)) ? 1 : -1
	return ps*perm
End
//Function test48()
//	Variable i,j,k,eps,n=1
//	Variable pos=0,neg=0
//	for(k=-3;k<4;k+=1)
//		for(j=-3;j<4;j+=1)
//			for(i=-3;i<4;i+=1)
//				eps = epsilon_ijk(i,j,k)
//				if (eps)
//					printf "%3d    %+d, %+d, %+d     %+d\r",n,i,j,k,eps
//					pos = (eps>0) ? pos+1 : pos
//					neg = (eps>0) ? neg+1 : neg
//					n += 1
//				endif
//			endfor
//		endfor
//	endfor
//	print "pos = ",pos,"     neg=",neg
//End



//	=========================================================================
//	=========================================================================
//	=========================================================================



Function/S GetGizmoNote()
	Execute "GetGizmo gizmoName"		// check if a gizmo is up
	SVAR S_GizmoName=S_GizmoName
	if (strlen(S_GizmoName)<1)
		return ""
	endif
	Execute "GetGizmo objectList"
	Wave/T TW_gizmoObjectList=TW_gizmoObjectList
	String str
	Variable i, i2=0, i1=2
	for (i=0;i<numpnts(TW_gizmoObjectList);i+=1)
		str = TW_gizmoObjectList[i]
		if (stringmatch(str," AppendToGizmo string=*,name=GizmoNoteString*"))
			i1 = strsearch(str,"\"",0)+1
			i2 = strsearch(str,"\"",i1+1)-1
			break
		endif
	endfor
	KillWaves/Z TW_gizmoObjectList
	KillStrings/Z S_gizmoObjectList
	return SelectString(i2<i1,str[i1,i2],"")
End
Function SetGizmoNote(noteStr)
	String noteStr
	// first see if noteStr exists
	if (strsearch(noteStr,";",0)>=0)
		Abort "this note trick for Gizmos does not work with semi-colons"
	endif

	Execute "GetGizmo gizmoName"		// check if a gizmo is up
	SVAR S_GizmoName=S_GizmoName
	if (strlen(S_GizmoName)<1)
		return 1
	endif
	Execute "GetGizmo objectList"
	SVAR S_gizmoObjectList=S_gizmoObjectList
	Variable i
	String cmd
	if (!(stringmatch(S_gizmoObjectList,"*AppendToGizmo string=*,name=GizmoNoteString*")))
		sprintf cmd, "AppendToGizmo string=\"%s\",strFont=\"\",name=GizmoNoteString",noteStr
	else
		sprintf cmd, "ModifyGizmo modifyObject=GizmoNoteString, property={string,\"%s\"}",noteStr
	endif
	Execute cmd
	KillWaves/Z TW_gizmoObjectList
	KillStrings/Z S_gizmoObjectList
End



Function initGrainGizmo()
	Variable doit = 0
	doit = (exists("root:Packages:Grains:explode")!=2) ? 1 : doit
	doit = (exists("root:Packages:Grains:Xcursor")!=2) ? 1 : doit
	doit = (exists("root:Packages:Grains:Ycursor")!=2) ? 1 : doit
	doit = (exists("root:Packages:Grains:Zcursor")!=2) ? 1 : doit
	doit = (exists("root:Packages:Grains:cursorGrainNum")!=2) ? 1 : doit
	doit = (exists("root:Packages:Grains:radiusBndry")!=2) ? 1 : doit
	doit = (exists("root:Packages:Grains:SymmetryOpsPath")!=2) ? 1 : doit
	doit = (exists("root:Packages:Grains:CubicSymmetryOps")!=1) ? 1 : doit
	doit = (exists("root:Packages:Grains:sigmaCSL")!=1) ? 1 : doit
	// always set micorn and degree to make it easier to move an experiment back and forth between Mac and Windows
	String/G root:Packages:Grains:micron=SelectString(stringmatch(igorInfo(2),"Macintosh"),"micron","µm")
	String/G root:Packages:Grains:degree=SelectString(stringmatch(igorInfo(2),"Macintosh"),"degree","¡")

	if (!doit)
		return 0
	endif
	NewDataFolder/O root:Packages
	NewDataFolder/O root:Packages:Grains
	Variable explode = NumVarOrDefault("root:Packages:Grains:explode",0)
	Variable Xcursor = NumVarOrDefault("root:Packages:Grains:Xcursor",0)
	Variable Ycursor = NumVarOrDefault("root:Packages:Grains:Ycursor",0)
	Variable Zcursor = NumVarOrDefault("root:Packages:Grains:Zcursor",0)
	Variable cursorGrainNum = NumVarOrDefault("root:Packages:Grains:cursorGrainNum",NaN)
	Variable radiusBndry = NumVarOrDefault("root:Packages:Grains:radiusBndry",Inf)
	String SymmetryOpsPath = StrVarOrDefault("root:Packages:Grains:SymmetryOpsPath","root:Packages:Grains:CubicSymmetryOps")

	Variable/G root:Packages:Grains:explode=explode
	Variable/G root:Packages:Grains:Xcursor=Xcursor
	Variable/G root:Packages:Grains:Ycursor=Ycursor
	Variable/G root:Packages:Grains:Zcursor=Zcursor
	Variable/G root:Packages:Grains:cursorGrainNum=cursorGrainNum
	Variable/G root:Packages:Grains:radiusBndry=radiusBndry
	String/G root:Packages:Grains:SymmetryOpsPath=SymmetryOpsPath
	MakeCubicSymmetryOps("root:Packages:Grains:")
	Make_sigmaCSL("root:Packages:Grains:")
End

Function Make_sigmaCSL(fldr)			// make a wave with data on Sigma boundaries for Cubic CSL
	// columns for sigmaCSL are "Sigma, angle¡, h,k,l
	String fldr							// folder used to put CubicSymmetryOps
	String wName = fldr+"sigmaCSL"
	Make/N=(1,1)/O $wName
	Wave sigmaCSL=$wName
	// This table from H. Mykura, H. (1980).  Grain-Boundary Structure and Kinetics, pp.  445-456.  Metals Park: American Society for Metals.  
	//	 Faxed to me by Budai Apr. 4, 2005
	sigmaCSL[0][0]= {3,5,7,9,11,13.1,13.2,15,17.1,17.2,19.1,19.2,21.1,21.2,23,25.1,25.2,27.1,27.2,29.1,29.2,31.1,31.2,33.1,33.2,33.3,35.1,35.2,37.1,37.2,37.3,39.1,39.2,41.1,41.2,41.3,43.1,43.2,43.3,45.1,45.2}
	sigmaCSL[41][0]= {45.3,47.1,47.2,49.1,49.2,49.3}
	sigmaCSL[0][1]= {60,36.87,38.21,38.94,50.48,22.62,27.8,48.19,28.07,61.93,26.53,46.83,21.79,44.4,40.45,16.25,51.68,31.58,35.42,43.6,46.39,17.9,52.19,20.05,33.55,58.98,34.04,43.23,18.92,43.13,50.57,32.21}
	sigmaCSL[32][1]= {50.13,12.68,40.88,55.88,15.18,27.91,60.77,28.62,36.87,53.13,37.07,43.66,43.58,43.58,49.22}
	sigmaCSL[0][2]= {1,0,1,0,0,0,1,0,0,1,0,1,1,1,1,0,1,0,0,0,1,1,1,0,1,0,1,1,0,0,1,1,1,0,0,0,1,0,2,1,1,1,1,0,1,1,2}
	sigmaCSL[0][3]= {1,0,1,1,1,0,1,1,0,2,1,1,1,1,1,0,3,1,1,0,2,1,1,1,1,1,1,3,0,1,1,1,2,0,1,1,1,1,3,1,2,2,3,2,1,1,2}
	sigmaCSL[0][4]= {1,1,1,1,1,1,1,2,1,2,1,1,1,2,3,1,3,1,2,1,2,1,2,1,3,1,2,3,1,3,1,1,3,1,2,1,1,2,3,3,2,2,3,3,1,5,3}

	wName = fldr+"hkl_temp__"				// do this section to make sure that hkl are ordered properly (123) is OK,  not (3,2,1)
	Make/N=3/O $wName
	Wave hkl=$wName
	Variable i
	for (i=0;i<DimSize(sigmaCSL,0);i+=1)
		hkl = sigmaCSL[i][p+2]
		hkl = abs(hkl)
		Sort hkl,hkl
		sigmaCSL[i][2,4] = hkl[q-2]
	endfor
	KillWaves/Z hkl
End

Function/S MakeCubicSymmetryOps(fldr)		// make a wave with the 48 cubic symmetry operations
//	the first 24 are proper rotations (determinant==1) and the second 24 involve mirrors (determinant==-1)
	String fldr							// folder used to put CubicSymmetryOps

	String wName = fldr+"CubicSymmetryOps"
	Make/N=(48,3,3)/D/O $wName
	Wave ops = $wName
	wName = UniqueName("mat",1,0)
	Make/N=(3,3)/D/O $wName
	Wave mat=$wName
	wName = UniqueName("axis",1,0)
	Make/N=3/D/O $wName
	Wave axis=$wName

	ops = 0
	// first the identity (1 operation)
	ops[0][0][0] = 1						// first is identity
	ops[0][1][1] = 1
	ops[0][2][2] = 1

	// the 3 4-fold axes (9 operations)
	axis = {1,0,0}
	rotationMatAboutAxis(axis,90,mat)		// 90¡ about x-axis
	ops[1][][] = mat[q][r]
	rotationMatAboutAxis(axis,180,mat)	// 180¡ about x-axis
	ops[2][][] = mat[q][r]
	rotationMatAboutAxis(axis,-90,mat)	// -90¡ about x-axis
	ops[3][][] = mat[q][r]

	axis = {0,1,0}
	rotationMatAboutAxis(axis,90,mat)		// 90¡ about y-axis
	ops[4][][] = mat[q][r]
	rotationMatAboutAxis(axis,180,mat)	// 180¡ about y-axis
	ops[5][][] = mat[q][r]
	rotationMatAboutAxis(axis,-90,mat)	// -90¡ about y-axis
	ops[6][][] = mat[q][r]

	axis = {0,0,1}
	rotationMatAboutAxis(axis,90,mat)		// 90¡ about y-axis
	ops[7][][] = mat[q][r]
	rotationMatAboutAxis(axis,180,mat)	// 180¡ about y-axis
	ops[8][][] = mat[q][r]
	rotationMatAboutAxis(axis,-90,mat)	// -90¡ about y-axis
	ops[9][][] = mat[q][r]

	// the 4 3-fold axes (8 operations)
	axis = {1,1,1}
	rotationMatAboutAxis(axis,120,mat)	// 120¡ about (1 1 1)
	ops[10][][] = mat[q][r]
	rotationMatAboutAxis(axis,-120,mat)	// -120¡ about (1 1 1)
	ops[11][][] = mat[q][r]

	axis = {-1,1,1}
	rotationMatAboutAxis(axis,120,mat)	// 120¡ about (-1 1 1)
	ops[12][][] = mat[q][r]
	rotationMatAboutAxis(axis,-120,mat)	// -120¡ about (-1 1 1)
	ops[13][][] = mat[q][r]

	axis = {1,-1,1}
	rotationMatAboutAxis(axis,120,mat)	// 120¡ about (1 -1 1)
	ops[14][][] = mat[q][r]
	rotationMatAboutAxis(axis,-120,mat)	// -120¡ about (1 -1 1)
	ops[15][][] = mat[q][r]

	axis = {-1,-1,1}
	rotationMatAboutAxis(axis,120,mat)	// 120¡ about (-1 -1 1)
	ops[16][][] = mat[q][r]
	rotationMatAboutAxis(axis,-120,mat)	// -120¡ about (-1 -1 1)
	ops[17][][] = mat[q][r]

	// the 6 2-fold axes (8 operations)
	axis = {1,0,1}
	rotationMatAboutAxis(axis,180,mat)	// 180¡ about (1 0 1)
	ops[18][][] = mat[q][r]

	axis = {0,1,1}
	rotationMatAboutAxis(axis,180,mat)	// 180¡ about (0 1 1)
	ops[19][][] = mat[q][r]

	axis = {-1,0,1}
	rotationMatAboutAxis(axis,180,mat)	// 180¡ about (-1 0 1)
	ops[20][][] = mat[q][r]

	axis = {0,-1,1}
	rotationMatAboutAxis(axis,180,mat)	// 180¡ about (0 -1 1)
	ops[21][][] = mat[q][r]

	axis = {1,1,0}
	rotationMatAboutAxis(axis,180,mat)	// 180¡ about (1 1 0)
	ops[22][][] = mat[q][r]

	axis = {1,-1,0}
	rotationMatAboutAxis(axis,180,mat)	// 180¡ about (1 -1 0)
	ops[23][][] = mat[q][r]

	Variable i
	for (i=0;i<24;i+=1)					// add the mirror operations
		ops[i+24][][] = -ops[i][q][r]
	endfor

	ops = abs(ops[p][q][r])<1e-13 ? 0 : ops[p][q][r]
	ops = abs(1-ops[p][q][r])<1e-13 ? 1 : ops[p][q][r]
	ops = abs(1+ops[p][q][r])<1e-13 ? -1 : ops[p][q][r]
	KillWaves/Z mat,axis
	Note/K ops
	Note ops, "Nproper=24;"				// only the first 24 operations are proper
	return GetWavesDataFolder(ops,2)
End
//	CubicSymmetryOps[0][0][0]= {1,1,1,1,0,-1,0,0,-1,0,0,0,0,0,0,0,0,0,0,-1,0,-1,0,0,-1,-1,-1,-1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0}
//	CubicSymmetryOps[0][1][0]= {0,0,0,0,0,0,0,-1,0,1,0,1,-1,0,-1,0,0,1,0,0,0,0,1,-1,0,0,0,0,0,0,0,1,0,-1,0,-1,1,0,1,0,0,-1,0,0,0,0,-1,1}
//	CubicSymmetryOps[0][2][0]= {0,0,0,0,1,0,-1,0,0,0,1,0,0,-1,0,1,-1,0,1,0,-1,0,0,0,0,0,0,0,-1,0,1,0,0,0,-1,0,0,1,0,-1,1,0,-1,0,1,0,0,0}
//	CubicSymmetryOps[0][0][1]= {0,0,0,0,0,0,0,1,0,-1,1,0,0,-1,0,-1,1,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,-1,0,1,-1,0,0,1,0,1,-1,0,0,0,0,0,-1,1}
//	CubicSymmetryOps[0][1][1]= {1,0,-1,0,1,1,1,0,-1,0,0,0,0,0,0,0,0,0,-1,0,-1,0,0,0,-1,0,1,0,-1,-1,-1,0,1,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0}
//	CubicSymmetryOps[0][2][1]= {0,-1,0,1,0,0,0,0,0,0,0,1,1,0,-1,0,0,-1,0,1,0,-1,0,0,0,1,0,-1,0,0,0,0,0,0,0,-1,-1,0,1,0,0,1,0,-1,0,1,0,0}
//	CubicSymmetryOps[0][0][2]= {0,0,0,0,-1,0,1,0,0,0,0,1,-1,0,1,0,0,-1,1,0,-1,0,0,0,0,0,0,0,1,0,-1,0,0,0,0,-1,1,0,-1,0,0,1,-1,0,1,0,0,0}
//	CubicSymmetryOps[0][1][2]= {0,1,0,-1,0,0,0,0,0,0,1,0,0,1,0,-1,-1,0,0,1,0,-1,0,0,0,-1,0,1,0,0,0,0,0,0,-1,0,0,-1,0,1,1,0,0,-1,0,1,0,0}
//	CubicSymmetryOps[0][2][2]= {1,0,-1,0,0,-1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,-1,0,1,0,0,1,0,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,1,1}
//
//	The follow is a test routine for the cubic symmetry ops
//
// check the following for each symmetry op
//	1 the determinant is ±1
//	2 the sum of the absolute value is 3 (i.e. 3 non-zero elements of ±1)
//	3 each operation is unique
Function checkCubicSymOps(printAll)
	Variable printAll
	if (numtype(printAll) || printAll<0 || printAll>1)
		printAll = (printAll>0) ? 2 : 1
		Prompt printAll, "print all the symmetry matrices?", popup, "only errors;print matrices"
		DoPrompt "print all",printAll
		if (V_flag)
			return 1
		endif
		printAll -= 1
	endif
	Wave cubicOps = root:Packages:Grains:CubicSymmetryOps

	Make/N=(3,3)/O/D op,opj,tempSum33
	Variable i,j,N=DimSize(cubicOps,0)
	Variable det,sumAbs,err=0
	for (i=0;i<N;i+=1)
		op = cubicOps[i][p][q]
		if (printAll)
			printSymmetryMatrix(i)
			if (i==24)
				print "***********************"
			endif
			print ""
		endif
		det = MatrixDet(op)
		tempSum33 = abs(op)
		sumAbs = sum(tempSum33)
		if ((abs(det-1)>1e-12 && abs(det+1)>1e-12) || abs(sumAbs-3)>1e-12)
			printf "checking op %d, Det(op)=%g,  Sum|op[i][j]|=%g\r",i,det,sumAbs
			print3x3(op)
			err += 1
		endif
		for (j=i+1;j<N;j+=1)				// check for two equal matrices
			MatrixOp/O temp33 = abs( op - opj )
			if (sum(temp33)<.1)			// matricies match
				printf "checking op %d, Det(op)=%g,  Sum|op[i][j]|=%g\r",i,det,sumAbs
				printf "op %d and %d match\r",i,j
				print3x3(op)
				print3x3(opj)
				err += 1
			endif
		endfor
	endfor
	if (err>0)
		printf "fournd %d errors\r",err
	else
		printf "no errors found\r"
	endif
	KIllWaves/Z op,opj,tempSum33,temp33
End
Static Function printSymmetryMatrix(i)
	Variable i				// index, of which cubic matrix to print i in [0,47]
	Wave symOps=$StrVarOrDefault("root:Packages:Grains:SymmetryOpsPath","")	// wave with symmetry ops

	Make/N=(3,3)/O/D temp_printSymmetryMatrix
	Wave op = temp_printSymmetryMatrix
	op = symOps[i][p][q]
	Make/N=3/O/D temp_axis_
	Wave axis = temp_axis_
	Variable angle = axisOfMatrix(op,axis,squareUp=1)
	integerDirectionApprox(axis,12)

	Variable det = MatrixDet(op)
	Variable a0=op[0][0], a1=op[1][1], a2=op[2][2]	// elements on trace

	String str
	if (det==1 && angle==0)
		str = "Identity matrix"
	else
		String degree = StrVarOrDefault("root:Packages:Grains:degree","¡")
		sprintf str, "a %g%s rotation of about the   (%g, %g, %g)",angle,axis[0],axis[1],axis[2]
	endif
	printf "symOp[%02d] =\t%+d, %+d, %+d   %s\r",i,symOps[i][0][0],symOps[i][0][1],symOps[i][0][2],str
	printf "\t\t\t\t%+d, %+d, %+d\r",symOps[i][1][0],symOps[i][1][1],symOps[i][1][2]
	printf "\t\t\t\t%+d, %+d, %+d\r",symOps[i][2][0],symOps[i][2][1],symOps[i][2][2]

	KillWaves/Z temp_printSymmetryMatrix,temp_axis_
	return angle
End
Static Function print3x3(mat3x3)
	Wave mat3x3
	printf "%s = \r",NameOfwave(mat3x3)
	printf "\t%+.4f, %+.4f, %+.4f\r",mat3x3[0][0],mat3x3[0][1],mat3x3[0][2]
	printf "\t%+.4f, %+.4f, %+.4f\r",mat3x3[1][0],mat3x3[1][1],mat3x3[1][2]
	printf "\t%+.4f, %+.4f, %+.4f\r",mat3x3[2][0],mat3x3[2][1],mat3x3[2][2]
End
