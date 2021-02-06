#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma IgorVersion = 8.0
#include "GizmoMarkers", version>=2.04
//#if (IgorVersion()<7)
//#include "GizmoZoomTranslate", version>=2.00
//#endif


#if defined(ZONE_TESTING) || defined(QS_TESTING) || defined(ZONE_QS_TESTING) || defined(DOTS_TESTING)
Menu "Zones"
	MenuItemIfWaveClassExists("Make Gizmo Waves","Zones*","MINCOLS:4"),GizmoMakeWavesForZones($"",$"",$"")
	MenuItemIfWaveExists("Gizmo of Zones",":GizmoWaves:GizmoGhats"), GizmoOfZones()
	MenuItemIfWaveExists("Gizmo of all Pair Rotations",":GizmoWaves:PairRotationsView"), GizmoOfMakePairRotations()
	"  View Gizmo from Detector", ModifyGizmo stopRotation, SETQUATERNION={-0.653030,-0.653030,-0.271204,0.271204}
	"  View Gizmo from Side", ModifyGizmo stopRotation, SETQUATERNION={0,-1/sqrt(2),0,1/sqrt(2)}
	"  View Gizmo at Beamline", ModifyGizmo stopRotation, SETQUATERNION={-0,1/sqrt(2),0,1/sqrt(2)}
	"-"
End
#endif


//  ====================================================================================  //
//  ================================== Start of Gizmo ==================================  //

Function GizmoMakeWavesForZones(GhatsMeasured,ZoneAxes,ZoneZone)
	// for a Gizmo, make the waves: GizmoGhats, GizmoZoneLines, GizmoZoneCircles, GizmoZoneZonePoints
	Wave GhatsMeasured, ZoneAxes			// only uses first 3 columns
	Wave ZoneZone

	if (!WaveExists(GhatsMeasured) || !WaveExists(ZoneAxes))
		String GhatList = WaveListClass("QsMeasured","*","DIMS:2,MINCOLS:3,MINROWS:1")
		String ZAxesList = WaveListClass("Zones","*","DIMS:2,MINCOLS:3,MINROWS:1")
		String ZZList = WaveListClass("ZoneZone","*","DIMS:2,MINCOLS:3,MINROWS:1")
		if (ItemsInList(GhatList)<1 || ItemsInList(ZAxesList)<1)
			return 1
		elseif (ItemsInList(GhatList)==1 && ItemsInList(ZAxesList)==1)
			Wave GhatsMeasured = $StringFromList(0,GhatList)
			Wave ZoneAxes = $StringFromList(0,ZAxesList)
		else
			String Gname, Zname, ZZname
			Prompt Gname,"Ghats Measured",popup, GhatList
			Prompt Zname,"Zone Axes",popup, ZAxesList
			Prompt ZZname,"ZoneZone",popup, ZZList
			DoPrompt "Waves for Gizmo", Gname, Zname, ZZname
			if (V_flag)
				return 1
			endif
			Wave GhatsMeasured = $Gname
			Wave ZoneAxes = $Zname
			Wave ZoneZone = $ZZname
		endif
	endif

	Variable Ntest=DimSize(GhatsMeasured,0), Nz=DimSize(ZoneAxes,0), Nzz=DimSize(ZoneZone,0)
	if (!(Nz>0 && Ntest>0))
		return 1
	endif
	Variable i,Ncir=50						// number of intervals in a circle

	Make/FREE RGBA0={0.2,0.2,0.2, 0.2}	// default color (RGBA)
	Make/N=(7,4)/FREE RGBAfirst
	RGBAfirst[0][0]= {1,0,0,1,0,1,0}	// special colors for the first Ncolors zones
	RGBAfirst[0][1]= {0,1,0,0,1,1,0}	// red, green, blue, magenta, cyan, yellow, dark
	RGBAfirst[0][2]= {0,0,1,1,1,0,0}
	RGBAfirst[0][3]= {1,1,1,1,1,1,0.5}
	Variable Ncolors=DimSize(RGBAfirst,0)	// only color the first Ncolors

	String fldrSav0= GetDataFolder(1)
	NewDataFolder/O/S GizmoWaves			// folder for Gizmo Waves
	Make/N=(Ntest,3)/O GizmoGhats=NaN
	Make/N=(Nz*3-1,3)/O GizmoZoneLines=NaN
	Make/N=(Nz*3-1,4)/O GizmoZoneLinesRGBA=NaN
	Make/N=(Nz*(Ncir+2)-1,3)/O GizmoZoneCircles=NaN
	Make/N=(Nz*(Ncir+2)-1,4)/O GizmoZoneCirclesRGBA=NaN
	if (Nzz>0)
		Make/N=(Nzz,3)/O GizmoZoneZonePoints=NaN
	endif
	SetDataFolder fldrSav0

	GizmoGhats = GhatsMeasured[p][q]	// first 3 columns of GhatsMeasured is G^

	for (i=0;i<Nz;i+=1)						// lines of the Zone axes
		GizmoZoneLines[3*i][] = -ZoneAxes[i][q]
		GizmoZoneLines[3*i+1][] = ZoneAxes[i][q]
		if (i<Ncolors)
			GizmoZoneLinesRGBA[3*i][] = RGBAfirst[i][q]
			GizmoZoneLinesRGBA[3*i+1][] = RGBAfirst[i][q]
		else
			GizmoZoneLinesRGBA[3*i][] = RGBA0[q]
			GizmoZoneLinesRGBA[3*i+1][] = RGBA0[q]
		endif
	endfor

	Make/N=3/D/FREE axis
	Variable mz, rad
	for (mz=0;mz<Nz;mz+=1)					// for each zone, make a great circle, with axis GizmoZoneLines
		axis[] = ZoneAxes[mz][p]
		Wave zero = perpVector(axis)
		Cross axis, zero
		Wave W_Cross=W_Cross				// vector perpendicular to both axis and zero
		normalize(W_Cross)
		for (i=0;i<=Ncir;i+=1)
			rad = i*2*PI/Ncir
			GizmoZoneCircles[52*mz+i][] = cos(rad)*zero[q] + sin(rad)*W_Cross[q]
			if (mz<Ncolors)
				GizmoZoneCirclesRGBA[52*mz+i][] = RGBAfirst[mz][q]
			else
				GizmoZoneCirclesRGBA[52*mz+i][] = RGBA0[q]
			endif
		endfor
	endfor
	KillWaves/Z W_Cross

	if (Nzz>0)									// the Zone of Zones
		GizmoZoneZonePoints[][] = ZoneZone[p][q]
	endif
	return 0
End
//
// called by GizmoMakeWavesForZones()
Static Function/WAVE perpVector(axisIN)
	Wave axisIN

	Make/N=3/D/FREE perp, axis=axisIN
	normalize(axis)
	perp = abs(axis)
	WaveStats/Q/M=1 perp
	Variable imax=V_maxloc

	if (abs(axis[imax])>0.99999)		// axis is real close to a (001) type
		imax = mod(imax+1,DimSize(axis,0))
	endif
	perp = (p==imax) - axis[imax]*axis
	normalize(perp)
	return perp
End


Function GizmoOfZones()
	if (ItemsInList(WinList("GizmoZones",";","WIN:"+num2istr(GIZMO_WIN_BIT))))
		DoWindow/F GizmoZones
		return 1
	endif

	String fldr = GetDataFolder(1)
	Wave GhatScatter = $(fldr + "GizmoWaves:GizmoGhats")
	Wave GizmoZoneLines = $(fldr + "GizmoWaves:GizmoZoneLines")
	Wave GizmoZoneLinesRGBA = $(fldr + "GizmoWaves:GizmoZoneLinesRGBA")
	Wave GizmoZoneCircles = $(fldr + "GizmoWaves:GizmoZoneCircles")
	Wave GizmoZoneCirclesRGBA = $(fldr + "GizmoWaves:GizmoZoneCirclesRGBA")

	Wave GizmoZoneZonePoints = $(fldr + "GizmoWaves:GizmoZoneZonePoints")
	Wave testZoneHatsGizmo = $(fldr + "GizmoWaves:testZoneHatsGizmo")

	if (!WaveExists(GhatScatter) || !WaveExists(GizmoZoneLines) || !WaveExists(GizmoZoneCircles))
		return 1
	endif


	NewGizmo/N=GizmoZones/T="GizmoZones"/W=(1140,49,1593,502)
	ModifyGizmo startRecMacro=700
	ModifyGizmo scalingOption=63
	ModifyGizmo keepPlotSquare=1
	AppendToGizmo Scatter=GhatScatter,name=GhatScatter
	ModifyGizmo ModifyObject=GhatScatter,objectType=scatter,property={ scatterColorType,0}
	ModifyGizmo ModifyObject=GhatScatter,objectType=scatter,property={ markerType,0}
	ModifyGizmo ModifyObject=GhatScatter,objectType=scatter,property={ sizeType,0}
	ModifyGizmo ModifyObject=GhatScatter,objectType=scatter,property={ rotationType,0}
	ModifyGizmo ModifyObject=GhatScatter,objectType=scatter,property={ Shape,2}
	ModifyGizmo ModifyObject=GhatScatter,objectType=scatter,property={ size,0.15}
	ModifyGizmo ModifyObject=GhatScatter,objectType=scatter,property={ color,0,0,0.8,1}

	AppendToGizmo Axes=boxAxes,name=axes0
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={-1,axisScalingMode,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,axisColor,0.5,0.5,0.5,0.4}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={1,axisColor,0.5,0.5,0.5,0.4}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={2,axisColor,0.5,0.5,0.5,0.4}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={3,axisColor,0.5,0.5,0.5,0.4}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={4,axisColor,0.5,0.5,0.5,0.4}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={5,axisColor,0.5,0.5,0.5,0.4}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={6,axisColor,0.5,0.5,0.5,0.4}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={7,axisColor,0.5,0.5,0.5,0.4}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={8,axisColor,0.5,0.5,0.5,0.4}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={9,axisColor,0.5,0.5,0.5,0.4}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={10,axisColor,0.5,0.5,0.5,0.4}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={11,axisColor,0.5,0.5,0.5,0.4}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,ticks,3}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={1,ticks,3}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={2,ticks,3}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,labelColor,0.5,0.5,0.5,0.4}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={1,labelColor,0.5,0.5,0.5,0.4}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={2,labelColor,0.5,0.5,0.5,0.4}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,axisLabel,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={1,axisLabel,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={2,axisLabel,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,axisLabelText,"X"}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={1,axisLabelText,"Y"}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={2,axisLabelText,"Z"}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,axisLabelDistance,0}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={1,axisLabelDistance,0}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={2,axisLabelDistance,0}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,axisLabelScale,1.5}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={1,axisLabelScale,1.5}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={2,axisLabelScale,1.5}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,axisLabelRGBA,0.5,0.5,0.5,0.5}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={1,axisLabelRGBA,0.5,0.5,0.5,0.5}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={2,axisLabelRGBA,0.5,0.5,0.5,0.5}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,labelBillboarding,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={1,labelBillboarding,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={2,labelBillboarding,1}
	ModifyGizmo modifyObject=axes0,objectType=Axes,property={-1,Clipped,0}

	AppendToGizmo Path=GizmoZoneLines,name=ZoneLines
	ModifyGizmo ModifyObject=ZoneLines,objectType=path,property={pathColorType,2}
	ModifyGizmo ModifyObject=ZoneLines,objectType=path,property={lineWidthType,1}
	ModifyGizmo ModifyObject=ZoneLines,objectType=path,property={lineWidth,2}
	if (WaveExists(GizmoZoneLinesRGBA))
		ModifyGizmo ModifyObject=ZoneLines,objectType=path,property={ pathColorWave,GizmoZoneLinesRGBA}
	endif

	AppendToGizmo Path=GizmoZoneCircles,name=ZoneCirclesPath
	ModifyGizmo ModifyObject=ZoneCirclesPath,objectType=path,property={pathColorType,2}
	ModifyGizmo ModifyObject=ZoneCirclesPath,objectType=path,property={lineWidthType,1}
	ModifyGizmo ModifyObject=ZoneCirclesPath,objectType=path,property={pathColorWave,GizmoZoneCirclesRGBA}

	if (WaveExists(GizmoZoneZonePoints))
		AppendToGizmo Scatter=GizmoZoneZonePoints,name=ZoneZoneSscatter
		ModifyGizmo ModifyObject=ZoneZoneSscatter,objectType=scatter,property={scatterColorType,0}
		ModifyGizmo ModifyObject=ZoneZoneSscatter,objectType=scatter,property={markerType,0}
		ModifyGizmo ModifyObject=ZoneZoneSscatter,objectType=scatter,property={sizeType,0}
		ModifyGizmo ModifyObject=ZoneZoneSscatter,objectType=scatter,property={rotationType,0}
		ModifyGizmo ModifyObject=ZoneZoneSscatter,objectType=scatter,property={Shape,2}
		ModifyGizmo ModifyObject=ZoneZoneSscatter,objectType=scatter,property={size,0.5}
		ModifyGizmo ModifyObject=ZoneZoneSscatter,objectType=scatter,property={color,0.25,0.75,1,0.3}
	endif

	if (WaveExists(testZoneHatsGizmo))
		AppendToGizmo Scatter=testZoneHatsGizmo,name=TestZoneHatsScatter
		ModifyGizmo ModifyObject=TestZoneHatsScatter,objectType=scatter,property={scatterColorType,0}
		ModifyGizmo ModifyObject=TestZoneHatsScatter,objectType=scatter,property={markerType,0}
		ModifyGizmo ModifyObject=TestZoneHatsScatter,objectType=scatter,property={sizeType,0}
		ModifyGizmo ModifyObject=TestZoneHatsScatter,objectType=scatter,property={rotationType,0}
		ModifyGizmo ModifyObject=TestZoneHatsScatter,objectType=scatter,property={Shape,2}
		ModifyGizmo ModifyObject=TestZoneHatsScatter,objectType=scatter,property={size,0.07}
		ModifyGizmo ModifyObject=TestZoneHatsScatter,objectType=scatter,property={color,0,0,0,1}
	endif

	// line with arrow showing incident beam
	AppendToGizmo line={0,0,-0.9,0,0,0}, name=IncidentBeamLine
	ModifyGizmo ModifyObject=IncidentBeamLine,objectType=line,property={ arrowMode,2}
	ModifyGizmo ModifyObject=IncidentBeamLine,objectType=line,property={ endArrowHeight,0.06}
	ModifyGizmo ModifyObject=IncidentBeamLine,objectType=line,property={ endArrowBase,0.02}
	ModifyGizmo ModifyObject=IncidentBeamLine,objectType=line,property={ cylinderStartRadius,0.005}
	ModifyGizmo ModifyObject=IncidentBeamLine,objectType=line,property={ cylinderEndRadius,0.005}

	AppendToGizmo attribute blendFunc={770,771},name=blendFunc0
	ModifyGizmo setDisplayList=0, opName=enable0, operation=enable, data=3042
	ModifyGizmo setDisplayList=1, opName=scale0, operation=scale, data={2,2,2}
	ModifyGizmo setDisplayList=2, opName=ortho0, operation=ortho, data={-2.5,2.5,-2.5,2.5,-3,3}
	ModifyGizmo setDisplayList=3, attribute=blendFunc0
	ModifyGizmo setDisplayList=4, object=IncidentBeamLine
	ModifyGizmo setDisplayList=5, object=GhatScatter
	ModifyGizmo setDisplayList=6, object=axes0
	ModifyGizmo setDisplayList=7, object=ZoneLines
	ModifyGizmo setDisplayList=8, object=ZoneCirclesPath

	if (GizmoObjExists("ZoneZoneSscatter", name="GizmoZones"))
		ModifyGizmo setDisplayList=-1, object=ZoneZoneSscatter
	endif
	if (GizmoObjExists("TestZoneHatsScatter", name="GizmoZones"))
		ModifyGizmo setDisplayList=-1, object=TestZoneHatsScatter
	endif
	ModifyGizmo autoscaling=1
	ModifyGizmo currentGroupObject=""
	ModifyGizmo endRecMacro
	ModifyGizmo SETQUATERNION={0,1/sqrt(2),0,1/sqrt(2)}
End


Function GizmoOfMakePairRotations()
	if (ItemsInList(WinList("GizmoPairRotations",";","WIN:"+num2istr(GIZMO_WIN_BIT))))
		DoWindow/F GizmoPairRotations
		return 0
	endif

	Wave PairRotationsView=:GizmoWaves:PairRotationsView
	Wave PairRotationsViewRGBA=:GizmoWaves:PairRotationsViewRGBA
	Wave PairRotationsViewSize=:GizmoWaves:PairRotationsViewSize
	if (!WaveExists(PairRotationsViewSize))
		return 1
	endif
	Wave corners = MakeGizmocubeCorners(PairRotationsView)
	corners = { {-PI,PI},{-PI,PI},{-PI,PI} }				// want to see total rotation space

	NewGizmo/N=GizmoPairRotations/T="GizmoPairRotations"/W=(1140,528,1593,981)
	ModifyGizmo startRecMacro=700
	ModifyGizmo scalingOption=63
	ModifyGizmo keepPlotSquare=1

	AppendToGizmo Scatter=PairRotationsView,name=PairRotationScatter
	ModifyGizmo ModifyObject=PairRotationScatter,objectType=scatter,property={ markerType,0}
	ModifyGizmo ModifyObject=PairRotationScatter,objectType=scatter,property={ rotationType,0}
	ModifyGizmo ModifyObject=PairRotationScatter,objectType=scatter,property={ Shape,2}
	if (WaveExists(PairRotationsViewSize))
		ModifyGizmo ModifyObject=PairRotationScatter,objectType=scatter,property={sizeType,1}
		ModifyGizmo ModifyObject=PairRotationScatter,objectType=scatter,property={sizeWave,PairRotationsViewSize}
	else
		ModifyGizmo ModifyObject=PairRotationScatter,objectType=scatter,property={sizeType,0}
		ModifyGizmo ModifyObject=PairRotationScatter,objectType=scatter,property={size,0.05}
	endif

	if (WaveExists(PairRotationsViewRGBA))
		ModifyGizmo ModifyObject=PairRotationScatter,objectType=scatter,property={ scatterColorType,1}
		ModifyGizmo ModifyObject=PairRotationScatter,objectType=scatter,property={ colorWave,PairRotationsViewRGBA}
	else
		ModifyGizmo ModifyObject=PairRotationScatter,objectType=scatter,property={ scatterColorType,0}
		ModifyGizmo ModifyObject=PairRotationScatter,objectType=scatter,property={ color, 1,0,0,1}
	endif

	AppendToGizmo Axes=boxAxes,name=axes0
	setGizmoAxisLabels("X","Y","Z")

	String title1="green: rotAxis from FindPeakInRots()", title2="red: symmetry reduced rotAxis", title3="black: MakeSimulatedTestSpots()"
	String Pink3dCross = AddGizmoMarkerGroup("",rgba="1,0.6,1,1",scale=0.15)

	AppendToGizmo Scatter=corners,name=CubeCornersScatter
	ModifyGizmo ModifyObject=CubeCornersScatter, objectType=scatter, property={markerType,0}
	ModifyGizmo ModifyObject=CubeCornersScatter, objectType=scatter, property={sizeType,0}
	ModifyGizmo ModifyObject=CubeCornersScatter, objectType=scatter, property={Shape,1}
	ModifyGizmo ModifyObject=CubeCornersScatter, objectType=scatter, property={size,1}

	AppendToGizmo attribute blendFunc={770,771},name=blendFunc0
	ModifyGizmo setDisplayList=0, opName=enable0, operation=enable, data=3042
	ModifyGizmo setDisplayList=1, opName=scale0, operation=scale, data={1,1,1}
	ModifyGizmo setDisplayList=2, opName=ortho0, operation=ortho, data={-2,2,-2,2,-3,3}
	ModifyGizmo setDisplayList=3, attribute=blendFunc0
	ModifyGizmo setDisplayList=4, object=$Pink3dCross
	ModifyGizmo setDisplayList=5, object=axes0
	ModifyGizmo setDisplayList=6, object=PairRotationScatter
	ModifyGizmo setDisplayList=7, object=CubeCornersScatter

	ModifyGizmo autoscaling=1
	ModifyGizmo currentGroupObject=""
	ModifyGizmo userString={cubecorners,"CubeCornersScatter"}
	ModifyGizmo endRecMacro
	ModifyGizmo SETQUATERNION={-0.347238,0.475219,-0.653693,0.475728}

	TextBox/C/N=textTitle/F=0/B=1/A=LT/X=2.05/Y=1.90 "\\Z18" + title1
	AppendText "\\Zr080red: symmetry reduced rotAxis" + title2
	AppendText "\\Zr080" + title3
	return 0
End


Function GizmoOfAllRotations()
	if (ItemsInList(WinList("GizmoAllRotations",";","WIN:"+num2istr(GIZMO_WIN_BIT))))
		DoWindow/F GizmoAllRotations
		return 0
	endif
	if(exists("NewGizmo")!=4)	// Do nothing if the Gizmo XOP is not available.
		DoAlert 0, "Gizmo XOP must be installed"
		return 1
	endif

	Wave AllRotations = AllRotationsView
	Wave AllRotationsSize = AllMatchesViewSize
	if (!WaveExists(AllRotations) || !WaveExists(AllRotationsSize))
		return 1
	endif
	Wave corners = MakeGizmocubeCorners(AllRotations)
	corners = { {-PI,PI},{-PI,PI},{-PI,PI} }				// want to see total rotation space

	Execute "NewGizmo/N=GizmoAllRotations/T=\"GizmoAllRotations\" /W=(1217,536,1670,989)"
	Execute "ModifyGizmo startRecMacro"
	Execute "AppendToGizmo Scatter="+GetWavesDataFolder(AllRotations,2)+",name=scatter0"
	Execute "ModifyGizmo ModifyObject=scatter0 property={ scatterColorType,0}"
	Execute "ModifyGizmo ModifyObject=scatter0 property={ markerType,0}"
	Execute "ModifyGizmo ModifyObject=scatter0 property={ sizeType,1}"
	Execute "ModifyGizmo ModifyObject=scatter0 property={ rotationType,0}"
	Execute "ModifyGizmo ModifyObject=scatter0 property={ Shape,1}"
	Execute "ModifyGizmo ModifyObject=scatter0 property={ size,2}"
	Execute "ModifyGizmo ModifyObject=scatter0 property={ color,0,0,0,1}"
	Execute "ModifyGizmo ModifyObject=scatter0 property={ sizeWave,"+GetWavesDataFolder(AllRotationsSize,2)+"}"
	Execute "AppendToGizmo Axes=boxAxes,name=axes0"

	setGizmoAxisLabels("Rx","Ry","Rz")

	Execute "AppendToGizmo Scatter="+GetWavesDataFolder(corners,2)+",name=CubeCornersScatter"
	Execute "ModifyGizmo ModifyObject=CubeCornersScatter property={ markerType,0}"
	Execute "ModifyGizmo ModifyObject=CubeCornersScatter property={ sizeType,0}"
	Execute "ModifyGizmo ModifyObject=CubeCornersScatter property={ Shape,1}"
	Execute "ModifyGizmo ModifyObject=CubeCornersScatter property={ size,1}"

	String titleGroupName = AddGizmoTitleGroup("","rotation vectors ±¹")
	
	Execute "ModifyGizmo  setDisplayList=0, object="+titleGroupName
	Execute "ModifyGizmo  setDisplayList=1, opName=MainTransform, operation=mainTransform"
	Execute "ModifyGizmo  setDisplayList=2, object=scatter0"
	Execute "ModifyGizmo  setDisplayList=3, object=axes0"
	Execute "ModifyGizmo  setDisplayList=4, object=CubeCornersScatter"

	Execute "ModifyGizmo  SETQUATERNION={0.842023,-0.117635,-0.006120,0.526424}"
	Execute "ModifyGizmo  autoscaling=1"
	Execute "ModifyGizmo  currentGroupObject=\"\""
	Execute "ModifyGizmo  compile"
	Execute "ModifyGizmo  userString={cubecorners,\"CubeCornersScatter\"}"

	Execute "ModifyGizmo  showAxisCue=1"
	Execute "ModifyGizmo  endRecMacro"
End

//  =================================== End of Gizmo ===================================  //
//  ====================================================================================  //

