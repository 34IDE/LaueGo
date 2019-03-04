#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include "GizmoMarkers", version>=2.04
#if (IgorVersion()<7)
#include "GizmoZoomTranslate", version>=2.00
#endif


#if defined(ZONE_TESTING) || defined(QS_TESTING) || defined(ZONE_QS_TESTING) || defined(DOTS_TESTING)
Menu "Zones"
	MenuItemIfWaveClassExists("Make Gizmo Waves","Zones*","MINCOLS:4"),GizmoMakeWavesForZones($"",$"",$"")
	MenuItemIfWaveExists("Gizmo of Zones",":GizmoWaves:GizmoGhats"), GizmoOfZones()
	MenuItemIfWaveExists("Gizmo of all Pair Rotations",":GizmoWaves:PairRotationsView"), GizmoOfMakePairRotations()
	"  View Gizmo from Detector", ModifyGizmo SETQUATERNION={-0.653030,-0.653030,-0.271204,0.271204}
	"  View Gizmo from Side", ModifyGizmo SETQUATERNION={0,-1/sqrt(2),0,1/sqrt(2)}
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
	if(exists("NewGizmo")!=4)	// Do nothing if the Gizmo XOP is not available.
		DoAlert 0, "Gizmo XOP must be installed"
		return 1
	endif

	String fldr = GetDataFolder(1)

	Execute "NewGizmo/N=GizmoZones/T=\"GizmoZones\" /W=(1140,49,1593,502)"
	Execute "ModifyGizmo startRecMacro"
	Execute "AppendToGizmo Scatter="+fldr+"GizmoWaves:GizmoGhats"+",name=GhatScatter"
	Execute "ModifyGizmo ModifyObject=GhatScatter property={ scatterColorType,0}"
	Execute "ModifyGizmo ModifyObject=GhatScatter property={ markerType,0}"
	Execute "ModifyGizmo ModifyObject=GhatScatter property={ sizeType,0}"
	Execute "ModifyGizmo ModifyObject=GhatScatter property={ rotationType,0}"
	Execute "ModifyGizmo ModifyObject=GhatScatter property={ Shape,1}"
	Execute "ModifyGizmo ModifyObject=GhatScatter property={ size,5}"
	Execute "ModifyGizmo ModifyObject=GhatScatter property={ color,0,6.10361e-05,0.8,1}"
	Execute "AppendToGizmo Axes=boxAxes,name=axes0"
	Execute "ModifyGizmo ModifyObject=axes0,property={0,axisRange,-1,-1,-1,1,-1,-1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={1,axisRange,-1,-1,-1,-1,1,-1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={2,axisRange,-1,-1,-1,-1,-1,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={3,axisRange,-1,1,-1,-1,1,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={4,axisRange,1,1,-1,1,1,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={5,axisRange,1,-1,-1,1,-1,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={6,axisRange,-1,-1,1,-1,1,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={7,axisRange,1,-1,1,1,1,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={8,axisRange,1,-1,-1,1,1,-1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={9,axisRange,-1,1,-1,1,1,-1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={10,axisRange,-1,1,1,1,1,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={11,axisRange,-1,-1,1,1,-1,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={-1,axisScalingMode,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={-1,axisColor,0.5,0.5,0.5,0.4}"
	Execute "ModifyGizmo ModifyObject=axes0,property={0,ticks,3}"
	Execute "ModifyGizmo ModifyObject=axes0,property={1,ticks,3}"
	Execute "ModifyGizmo ModifyObject=axes0,property={2,ticks,3}"
	Execute "ModifyGizmo ModifyObject=axes0,property={0,labelColor,0.5,0.5,0.5,0.4}"
	Execute "ModifyGizmo ModifyObject=axes0,property={1,labelColor,0.5,0.5,0.5,0.4}"
	Execute "ModifyGizmo ModifyObject=axes0,property={2,labelColor,0.5,0.5,0.5,0.4}"
	Execute "ModifyGizmo ModifyObject=axes0,property={0,axisLabel,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={1,axisLabel,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={2,axisLabel,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={0,axisLabelText,\"X\"}"
	Execute "ModifyGizmo ModifyObject=axes0,property={1,axisLabelText,\"Y\"}"
	Execute "ModifyGizmo ModifyObject=axes0,property={2,axisLabelText,\"Z\"}"
	Execute "ModifyGizmo ModifyObject=axes0,property={0,axisLabelCenter,0}"
	Execute "ModifyGizmo ModifyObject=axes0,property={1,axisLabelCenter,0}"
	Execute "ModifyGizmo ModifyObject=axes0,property={2,axisLabelCenter,0}"
	Execute "ModifyGizmo ModifyObject=axes0,property={0,axisLabelDistance,0}"
	Execute "ModifyGizmo ModifyObject=axes0,property={1,axisLabelDistance,0}"
	Execute "ModifyGizmo ModifyObject=axes0,property={2,axisLabelDistance,0}"
	Execute "ModifyGizmo ModifyObject=axes0,property={0,axisLabelScale,0.5}"
	Execute "ModifyGizmo ModifyObject=axes0,property={1,axisLabelScale,0.5}"
	Execute "ModifyGizmo ModifyObject=axes0,property={2,axisLabelScale,0.5}"
	Execute "ModifyGizmo ModifyObject=axes0,property={0,axisLabelRGBA,0.5,0.5,0.5,0.5}"
	Execute "ModifyGizmo ModifyObject=axes0,property={1,axisLabelRGBA,0.5,0.5,0.5,0.5}"
	Execute "ModifyGizmo ModifyObject=axes0,property={2,axisLabelRGBA,0.5,0.5,0.5,0.5}"
	Execute "ModifyGizmo modifyObject=axes0 property={Clipped,0}"

	Execute "AppendToGizmo Path="+fldr+"GizmoWaves:GizmoZoneLines"+",name=ZoneLines"
	Execute "ModifyGizmo ModifyObject=ZoneLines property={ pathColorType,2}"
	Execute "ModifyGizmo ModifyObject=ZoneLines property={ lineWidthType,1}"
	Execute "ModifyGizmo ModifyObject=ZoneLines property={ lineWidth,2}"
	Execute "ModifyGizmo ModifyObject=ZoneLines property={ pathColorWave,"+fldr+"GizmoWaves:GizmoZoneLinesRGBA"+"}"
	Execute "AppendToGizmo Path="+fldr+"GizmoWaves:GizmoZoneCircles"+",name=ZoneCirclesPath"
	Execute "ModifyGizmo ModifyObject=ZoneCirclesPath property={ pathColorType,2}"
	Execute "ModifyGizmo ModifyObject=ZoneCirclesPath property={ lineWidthType,1}"
	Execute "ModifyGizmo ModifyObject=ZoneCirclesPath property={ lineWidth,1}"
	Execute "ModifyGizmo ModifyObject=ZoneCirclesPath property={ pathColorWave,"+fldr+"GizmoWaves:GizmoZoneCirclesRGBA"+"}"

	if (exists(GetDataFolder(1)+"GizmoWaves:GizmoZoneZonePoints"))
		Execute "AppendToGizmo Scatter="+fldr+"GizmoWaves:GizmoZoneZonePoints"+",name=ZoneZoneSscatter"
		Execute "ModifyGizmo ModifyObject=ZoneZoneSscatter property={ scatterColorType,0}"
		Execute "ModifyGizmo ModifyObject=ZoneZoneSscatter property={ markerType,0}"
		Execute "ModifyGizmo ModifyObject=ZoneZoneSscatter property={ sizeType,0}"
		Execute "ModifyGizmo ModifyObject=ZoneZoneSscatter property={ rotationType,0}"
		Execute "ModifyGizmo ModifyObject=ZoneZoneSscatter property={ Shape,2}"
		Execute "ModifyGizmo ModifyObject=ZoneZoneSscatter property={ size,0.5}"
		Execute "ModifyGizmo ModifyObject=ZoneZoneSscatter property={ color,0.250019,0.748074,1,0.3}"
	endif


	Wave testZoneHatsGizmo=$(fldr+"GizmoWaves:testZoneHatsGizmo")
	Variable showTestZoneHats = WaveExists(testZoneHatsGizmo)
	if (showTestZoneHats)
		Execute "AppendToGizmo Scatter="+fldr+"GizmoWaves:testZoneHatsGizmo"+",name=TestZoneHats"
		Execute "ModifyGizmo ModifyObject=TestZoneHats property={ scatterColorType,0}"
		Execute "ModifyGizmo ModifyObject=TestZoneHats property={ markerType,0}"
		Execute "ModifyGizmo ModifyObject=TestZoneHats property={ sizeType,0}"
		Execute "ModifyGizmo ModifyObject=TestZoneHats property={ rotationType,0}"
		Execute "ModifyGizmo ModifyObject=TestZoneHats property={ Shape,1}"
		Execute "ModifyGizmo ModifyObject=TestZoneHats property={ size,6}"
		Execute "ModifyGizmo ModifyObject=TestZoneHats property={ color,0,0,0,1}"
	endif

	Execute "AppendToGizmo attribute blendFunc={770,771},name=blendFunc0"
	Execute "ModifyGizmo setDisplayList=0, opName=enable0, operation=enable, data=3042"
	Execute "ModifyGizmo setDisplayList=1, opName=scale0, operation=scale, data={1.5,1.5,1.5}"
	Execute "ModifyGizmo setDisplayList=2, opName=ortho0, operation=ortho, data={-2.5,2.5,-2.5,2.5,-3,3}"
	Execute "ModifyGizmo setDisplayList=3, attribute=blendFunc0"
	Execute "ModifyGizmo setDisplayList=4, object=GhatScatter"
	Execute "ModifyGizmo setDisplayList=5, object=axes0"
	Execute "ModifyGizmo setDisplayList=6, object=ZoneLines"
	Execute "ModifyGizmo setDisplayList=7, object=ZoneCirclesPath"
	if (showTestZoneHats)
		Execute "ModifyGizmo setDisplayList=-1, object=TestZoneHats"
	endif
	if (exists(GetDataFolder(1)+"GizmoWaves:GizmoZoneZonePoints"))
		Execute "ModifyGizmo setDisplayList=-1, object=ZoneZoneSscatter"
	endif

	Execute "ModifyGizmo SETQUATERNION={0.350994,-0.719064,0.202320,-0.564657}"
	Execute "ModifyGizmo autoscaling=1"
	Execute "ModifyGizmo currentGroupObject=\"\""
	Execute "ModifyGizmo compile"

	Execute "ModifyGizmo endRecMacro"
End


Function GizmoOfMakePairRotations()
	if (ItemsInList(WinList("GizmoPairRotations",";","WIN:"+num2istr(GIZMO_WIN_BIT))))
		DoWindow/F GizmoPairRotations
		return 0
	endif
	if(exists("NewGizmo")!=4)	// Do nothing if the Gizmo XOP is not available.
		DoAlert 0, "Gizmo XOP must be installed"
		return 1
	endif

	Wave PairRotationsView=:GizmoWaves:PairRotationsView
	Wave PairRotationsViewRGBA=:GizmoWaves:PairRotationsViewRGBA
	Wave PairRotationsViewSize=:GizmoWaves:PairRotationsViewSize
	if (!WaveExists(PairRotationsViewSize))
		return 1
	endif
	Wave corners = MakeGizmocubeCorners(PairRotationsView)
	corners = { {-PI,PI},{-PI,PI},{-PI,PI} }				// want to see total rotation space

	Execute "NewGizmo/N=GizmoPairRotations/T=\"GizmoPairRotations\" /W=(1140,528,1593,981)"
	Execute "ModifyGizmo startRecMacro"

	Execute "AppendToGizmo Scatter="+GetWavesDataFolder(PairRotationsView,2)+",name=PairRotationScatter"
	Execute "ModifyGizmo ModifyObject=PairRotationScatter property={ markerType,0}"
	Execute "ModifyGizmo ModifyObject=PairRotationScatter property={ rotationType,0}"
	Execute "ModifyGizmo ModifyObject=PairRotationScatter property={ Shape,1}"
	if (WaveExists(PairRotationsViewSize))
		Execute "ModifyGizmo ModifyObject=PairRotationScatter property={ sizeType,1}"
		Execute "ModifyGizmo ModifyObject=PairRotationScatter property={ sizeWave,"+GetWavesDataFolder(PairRotationsViewSize,2)+"}"
	else
		Execute "ModifyGizmo ModifyObject=PairRotationScatter property={ sizeType,0}"
		Execute "ModifyGizmo ModifyObject=PairRotationScatter property={ size,2}"
	endif
	if (WaveExists(PairRotationsViewRGBA))
		Execute "ModifyGizmo ModifyObject=PairRotationScatter property={ scatterColorType,1}"
		Execute "ModifyGizmo ModifyObject=PairRotationScatter property={ colorWave,"+GetWavesDataFolder(PairRotationsViewRGBA,2)+"}"
	else
		Execute "ModifyGizmo ModifyObject=PairRotationScatter property={ scatterColorType,1}"
		Execute "ModifyGizmo ModifyObject=PairRotationScatter property={ color,0,0,0,1}"
	endif

	Execute "AppendToGizmo Axes=boxAxes,name=axes0"
	setGizmoAxisLabels("X","Y","Z")
	String titleGroupName = AddGizmoTitleGroup("","green: rotAxis from FindPeakInRots()",title2="red: symmetry reduced rotAxis",title3="black: MakeSimulatedTestSpots()")
	String Pink3dCross = AddGizmoMarkerGroup("",rgba="1,0.6,1,1",scale=0.15)

	Execute "AppendToGizmo Scatter="+GetWavesDataFolder(corners,2)+",name=CubeCornersScatter"
	Execute "ModifyGizmo ModifyObject=CubeCornersScatter property={ markerType,0}"
	Execute "ModifyGizmo ModifyObject=CubeCornersScatter property={ sizeType,0}"
	Execute "ModifyGizmo ModifyObject=CubeCornersScatter property={ Shape,1}"
	Execute "ModifyGizmo ModifyObject=CubeCornersScatter property={ size,1}"

	Execute "AppendToGizmo attribute blendFunc={770,771},name=blendFunc0"
	Execute "ModifyGizmo setDisplayList=0, object="+titleGroupName
	Execute "ModifyGizmo setDisplayList=1, opName=MainTransform, operation=mainTransform"
	Execute "ModifyGizmo setDisplayList=2, opName=enable0, operation=enable, data=3042"
	Execute "ModifyGizmo setDisplayList=3, opName=scale0, operation=scale, data={1.25,1.25,1.25}"
	Execute "ModifyGizmo setDisplayList=4, opName=ortho0, operation=ortho, data={-2,2,-2,2,-3,3}"
	Execute "ModifyGizmo setDisplayList=5, attribute=blendFunc0"
	Execute "ModifyGizmo setDisplayList=6, object="+Pink3dCross
	Execute "ModifyGizmo setDisplayList=7, object=axes0"
	Execute "ModifyGizmo setDisplayList=8, object=PairRotationScatter"
	Execute "ModifyGizmo setDisplayList=9, object=CubeCornersScatter"
	Execute "ModifyGizmo SETQUATERNION={-0.917130,-0.171883,-0.272546,-0.234738}"
	Execute "ModifyGizmo autoscaling=1"
	Execute "ModifyGizmo currentGroupObject=\"\""
	Execute "ModifyGizmo compile"

	Execute "ModifyGizmo userString={cubecorners,\"CubeCornersScatter\"}"
	Execute "ModifyGizmo endRecMacro"
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

