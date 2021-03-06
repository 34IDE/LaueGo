#pragma rtGlobals=1		// Use modern global access method.
#pragma version = 1.31

//	#include "GizmoGrains"
#include "NeighborByWalking", version>=1.0
#include "HistogramAxes", version>=1.0
#include "GrainMath", version>=1.4
#include "CutUpGrains", version>=1.0

Menu "Grains"
	"Gizmo of All Grains",MakeGizmoMultiGrain($"","",0,NaN)
	"Gizmo Around One Grain",GizmoAroundOneGrain($"",Nan,NaN,NaN)
	"Zoomer Panel for a Gizmo",PanelGizmoView()
	"Gizmo of All Rodriques Vectors", GizmoAllRodriquesVectors($"")
 	"-"
//	"(Info"
//	"List Neighbors to a Grain",makeBorderList($"",NaN)
//	"Statistics of a Boundary",BoundaryStats($"",NaN,NaN)
//	"Angle Between two Grains",angleBetweenGrains($"",NaN,NaN)
//	"Find Number of Grains Above a Size", NumOfGrainsBiggerThan($"",NaN)
//	"-"
//	"(New Data"
//	"Prepare\Read In data",PrepareData()
//	"Calculate Grains from Neighbor Rotations",findGrainsByComparingNeighbors(NaN,NaN,0)
End




//	=========================================================================
//	=========================================================================
//	=========================================================================

//	ModifyGizmo outputResFactor=n		// n=2, default

Function GizmoAroundOneGrain(grains,grainNum,minSize,minBndryPoints)		// it includes grainNum in the list too
	Wave grains						// 3d wave with all of the grain numbers (e.g. gNeighbor)
	Variable grainNum					// find the grains that border grainsNum
	Variable minSize					// mininum size to include
	Variable minBndryPoints			// mininum number of points in common needed for inclusion of a neighbor

	if (!WaveExists(grains) || WaveDims(grains)!=3 || numtype(grainNum+minSize) || grainNum<0 || !(minBndryPoints>=0))
		String grainName=NameOfWave(grains)
		Prompt grainName,"3d wave with all grains",popup,ListMultiGrainMats()
		grainNum = numtype(grainNum)==2 ? 0 : max(grainNum,0)
		Prompt grainNum, "grain number"
		minSize = numtype(minSize) ? 200 : minSize
		Prompt minSize, "min no. of voxels in an acceptable grain"
		minBndryPoints = !(minBndryPoints>=0) ? 3 : minBndryPoints
		Prompt minBndryPoints, "only include neighbors with at least this many bordering voxels"
		DoPrompt "choose grain",grainName,grainNum,minSize,minBndryPoints
		if (V_flag)
			return 1
		endif
		Wave grains = $grainName
		if (ItemsInList(GetRTStackInfo(0))<2)
			printf "   GizmoAroundOneGrain(%s,%g,%g,%g)\r",NameOfWave(grains),grainNum,minSize,minBndryPoints
		endif
	endif
	if (!WaveExists(grains))
		return 1
	endif
	WaveStats/Q grains
	Variable maxGrain = V_max
	if (numtype(grainNum) || grainNum<0 || grainNum>maxGrain)
		return 2
	endif
	if (!(minBndryPoints>=0))
		return 3
	endif

	String list = makeBorderList(grains,grainNum)
	Variable N=ItemsInList(list)
	Make/N=(N)/O very_temp_wave_
	Wave abc = very_temp_wave_
	abc = str2num(StringFromList(p,list))
	Sort abc,abc

	String gName = GetWavesDataFolder(grains,2)
	Wave vol=$(gName+"_Vol")
	Variable i,size
	for (i=0;i<N;i+=1)
		size = vol[abc[i]]
		if (size<3)
			abc[i] = NaN
		endif
	endfor
	list = ""
	for (i=0;i<N;i+=1)
		if (numtype(abc[i]))
			continue
		endif
		list = AddListItem(num2istr(abc[i]),list,",",Inf)
	endfor
	i = strlen(list)
	list[i-1,i]=""
	N = ItemsInList(list,",")

	// remove grains from list that do not have at least minBndryPoints points in common
	String str
	Variable j
	i = 1
	do
		j = str2num(StringFromList(i,list,","))
		str = BoundaryStats(grains,grainNum,j,Inf)
		if (NumberByKey("Nboundary1",str,"=")<minBndryPoints || NumberByKey("Nboundary2",str,"=")<minBndryPoints)
			list = RemoveFromList(num2istr(j),list,",")
			N -= 1
		else
			i += 1
		endif
	while (i<N)

//	if (N>1)
//		list = RemoveFromList(num2istr(grainNum), list,",")
//		list = num2istr(grainNum)+","+list
//	endif
	KillWaves/Z very_temp_wave_
	if (ItemsInList(GetRTStackInfo(0))<2)
		printf "   from '%s'. using %d grains: {%s}\r",NameOfWave(grains),itemsinrange(list),list
	endif
	if (N<1)
		return 1
	endif
	DoAlert 1, "Put Up a Gizmo for this List?"
	if (V_flag==1)
		i = MakeGizmoMultiGrain(grains,list,0.1,grainNum)
		if (i==0)
			PopGrainProc("popupGrainBrite",0,"grain"+num2istr(grainNum))
		endif
		// DoAlert 0, "Need to move the origin to center of mass of "+num2istr(grainNum)
	endif
	return 0
End



//Window GizmoMultiGrain() : GizmoPlot
//	PauseUpdate; Silent 1	// Building window...
//	MakeGizmoMultiGrain($"","",1,NaN)
//EndMacro
Function MakeGizmoMultiGrain(grainNos,range,showScaleCube,firstGrain) : GizmoPlot
	Wave grainNos					// 3d wave with all of the grain numbers (e.g. gNeighbor)
	String range					// list of grain numbers, possibly non-contiguous range, ie "1-20,25,30-40"
	Variable showScaleCube			// =1
	Variable firstGrain				// if a special first grain is desired, this one will be first in the list, and 'brite'
	Variable scaleCubeSize=5

	if (!WaveExists(grainNos) || WaveDims(grainNos)!=3 || numtype(NextInRange(range,-inf)) || numtype(showScaleCube) || abs(showScaleCube-.1)<1e-5)
		//	note, do not prompt for firstGrain
		range = SelectString(strlen(range),"0-Inf",range)
		Prompt range,"range of grain numbers"
		Prompt showScaleCube, "Show a Scale Cube", popup, "Show Scale Cube;No Scale Cube"
		showScaleCube = round(showScaleCube) ? 1 : 2
		String grainNoWave=NameOfWave(grainNos)
		Prompt grainNoWave,"3d wave with all grains",popup,ListMultiGrainMats()
		DoPrompt "grains to make",range,grainNoWave,showScaleCube
		if (V_flag)
			return 1
		endif
		Wave grainNos = $grainNoWave
		showScaleCube = (showScaleCube==1)
		if (ItemsInList(GetRTStackInfo(0))<2)
			printf "\tMakeGizmoMultiGrain(%s,\"%s\",%d,%g)\r",NameOfWave(grainNos),range,showScaleCube,firstGrain
		endif
	endif
	if (!WaveExists(grainNos))
		DoAlert 0, "'grainNos' wave does not exist"
		return 1
	elseif(exists("NewGizmo")!=4)						// Do nothing if the Gizmo XOP is not available.
		DoAlert 0, "Gizmo XOP must be installed"
		return 1
	endif
	Variable i
	if (lastInRange(range)==Inf)
		Wave vol=$(GetWavesDataFolder(grainNos,2)+"_Vol")
		i = strsearch(range,"-Inf",Inf,3)			// do a special check to see if this ends in '-Inf', note this is dash inf, not minus
		if (i+4 == strlen(range))
			range[i+1,Inf] = num2istr(numpnts(vol)-1)
		else
			DoAlert 0, "The range '"+range+"' makes no sense, max range is \"0-Inf\""
			return 1
		endif
	endif

	range = expandRange(range,",")
	i = WhichListItem(num2istr(firstGrain),range,",")
	firstGrain = numtype(i) ? NaN : firstGrain			// ensure invalid if grain not in list
	if (i>=0)											// firstGrain is in list
		range = RemoveFromList(num2istr(firstGrain), range,",")	// put firstGrain first in the list
		range = AddListItem(num2istr(firstGrain), range,",")
	endif

	String gList = makeRange3dGrain(grainNos,range) // list of individual grain arrays, one array for each grain
	if (ItemsInList(gList)<1)
		DoAlert 0, "No grains created by makeRange3dGrain()"
		return 1
	endif
	String cmd, str
	Variable nx,ny,nz, nz2
	nx = DimSize(grainNos,0)
	ny = DimSize(grainNos,1)
	nz = DimSize(grainNos,2)
	nz2 = nz/2
	String noteStr=""

	str = StringByKey("SCREEN1",IgorInfo(0))
	str = str[strsearch(str,"RECT=",0),Inf]
	Variable size = 0.85*str2num(StringFromList(3,str[strsearch(str,"RECT=",0),Inf],","))
	sprintf cmd, "NewGizmo/N=MultiGrainGizmo/W=(%d,%d,%d,%d)/K=1",90,50,90+size,50+size
	Execute cmd										//	Execute "NewGizmo/N=MultiGrainGizmo/W=(90,50,690,650)/K=1"
	Execute "ModifyGizmo startRecMacro"
	Execute "ModifyGizmo outputResFactor=2"		// n=2, default
	Execute "ModifyGizmo scalingMode=1"
	SetGizmoNote(" ")

	sprintf cmd, "AppendToGizmo freeAxesCue={0,0,0,%g},name=freeAxesCue0",nz2/3
	Execute cmd	
	Execute "AppendToGizmo attribute lineWidth=2, name=lineWidthAxesCue"
	Execute "AppendToGizmo light=Directional,name=light0"
	Execute "ModifyGizmo light=light0 property={ position,0.683157,-0.000063,-0.730271,0}"
	Execute "ModifyGizmo light=light0 property={ direction,0.683157,-0.000063,-0.730271}"
	Execute "ModifyGizmo light=light0 property={ specular,0.8,0.8,0.8,1.0}"
	Execute "ModifyGizmo light=light0 property={ diffuse,0.5,0.55,0.55,1.0}"
	Execute "AppendToGizmo light=Directional,name=light1"
	Execute "ModifyGizmo light=light1 property={ position,-0.791350,-0.503586,-0.346650,0}"
	Execute "ModifyGizmo light=light1 property={ direction,-0.791350,-0.503586,-0.346650}"
	Execute "ModifyGizmo light=light1 property={ specular,0.8,0.8,0.8,1.0}"
	Execute "ModifyGizmo light=light1 property={ diffuse,0.55,0.55,0.55,1.0}"
	Execute "AppendToGizmo attribute blendFunc={770,771},name=blendFunc0"
//	Execute "AppendToGizmo sphere={0.0,25,25},name=sphereCursor"			// initially size is 0, so invisible (25,25, gives a sphere)
	Execute "AppendToGizmo sphere={0.0,4,2},name=sphereCursor"				// initially size is 0, so invisible (4,2 gives an octahedron)
	Execute "AppendToGizmo attribute diffuse={1,0,0,1,1032},name=diffuseCursor"
	Execute "AppendToGizmo attribute specular={0.5,0.5,0.5,1,1029},name=specularGrains"

	// ensure it exists and append 'GrainBoundaryPath' which is the surface normal of a boundary if needed
	str = GetWavesDataFolder(grainNos,1)+"GrainBoundaryPath"
	if (exists(str)!=1)
		Make/N=(2,3) $str													// just a dummy if if does not exist, it will be updated later
	endif
	Wave GrainBoundaryPath=$str
	GrainBoundaryPath = 1														// why do I need a 1. why not 0 ?
	sprintf cmd, "AppendToGizmo Path=%s,name=GrainBoundaryPath",GetWavesDataFolder(GrainBoundaryPath,2)
	Execute cmd
	Execute "ModifyGizmo ModifyObject=GrainBoundaryPath property={ pathColorType,1}"
	Execute "ModifyGizmo ModifyObject=GrainBoundaryPath property={ lineWidthType,1}"
	Execute "ModifyGizmo ModifyObject=GrainBoundaryPath property={ lineWidth,4}"
	Execute "ModifyGizmo ModifyObject=GrainBoundaryPath property={ pathColor,0,0.6,0,1}"
	Execute "ModifyGizmo modifyObject=GrainBoundaryPath property={calcNormals,1}"
	Execute "AppendToGizmo sphere={0.5,25,25},name=sphereNormal"
	Execute "AppendToGizmo attribute diffuse={0.1,0.9,0.1,1,1032},name=diffuseNormal"

	// ensure it exists and append 'bndryScatterPoints' which will be the points showing a boundary if needed
	str = GetWavesDataFolder(grainNos,1)+"bndryScatterPoints"
	if (exists(str)!=1)
		Make/N=(2,3) $str													// if it does not exist, just create a dummy version
	endif
	Wave bndryScatterPoints=$str
	sprintf cmd, "AppendToGizmo Scatter=%s,name=scatterBndry",GetWavesDataFolder(bndryScatterPoints,2)
	Execute cmd
	Execute "ModifyGizmo ModifyObject=scatterBndry property={ scatterColorType,0}"
	Execute "ModifyGizmo ModifyObject=scatterBndry property={ markerType,0}"
	Execute "ModifyGizmo ModifyObject=scatterBndry property={ sizeType,0}"
	Execute "ModifyGizmo ModifyObject=scatterBndry property={ rotationType,0}"
	Execute "ModifyGizmo ModifyObject=scatterBndry property={ Shape,2}"
	Execute "ModifyGizmo ModifyObject=scatterBndry property={ size,0}"		// marker size=0, initiall invisible
	Execute "ModifyGizmo ModifyObject=scatterBndry property={ color,0,0,0,1}"

	if (showScaleCube)
		scaleCubeSize = round(min(nx,min(ny,nz))/2)
		sprintf cmd, "AppendToGizmo box={%g,%g,%g},name=scaleBox",scaleCubeSize,scaleCubeSize,scaleCubeSize
		Execute cmd
		Execute "ModifyGizmo modifyObject=scaleBox property={calcNormals,1}"
		Execute "AppendToGizmo attribute color={0.5,0.5,0.5,0.2},name=scaleBoxcolor"
		Execute "AppendToGizmo attribute diffuse={0.5,0.5,0.5,0.5,1028},name=scaleBoxdiffuse"
		Execute "AppendToGizmo attribute specular={0.2,0.2,0.2,1,1028},name=scaleBoxspecular"
		Execute "AppendToGizmo attribute ambient={1,1,1,1,1032},name=scaleBoxambient"
		sprintf cmd, "AppendToGizmo string=\"%g%s\",strFont=\"Times\",name=scaleString",scaleCubeSize,SelectString(stringmatch(igorInfo(2),"Macintosh"),"micron","�m")
		Execute cmd
		Execute "AppendToGizmo attribute color={0,0,0,1},name=scaleStringColor"
	endif

	String description,grainName,displayList=""
	Variable id, isoValue, N=ItemsInList(gList)
	for (i=0;i<N;i+=1)
		grainName = StringFromList(i,gList)
		str = note($grainName)
		id = NumberByKey("grainNum", str,"=")
		isoValue = NumberByKey("isoValue",str,"=")
		description = ""
		description = ReplaceStringByKey("name",description,"grain"+num2istr(id),"=")
		description = ReplaceStringByKey("WaveName",description,grainName,"=")
		description = ReplaceStringByKey("RGBA",description,StringByKey("RGBA",str,"="),"=")
		description = ReplaceNumberByKey("isoValue",description,isoValue,"=")
		AddIsoSurfaceGrainToObjectList(description)
		displayList += StringByKey("name",description,"=")+";"
	endfor

//	sprintf cmd "ModifyGizmo setDisplayList=0, opName=ortho0, operation=ortho, data={%g,%g,%g,%g,%g,%g}",-nz2,nz2,-nz2,nz2,-nz2,nz2
	sprintf cmd "ModifyGizmo setDisplayList=0, opName=ortho0, operation=ortho, data={%g,%g,%g,%g,%g,%g}",-nz2,nz2,-nz2,nz2,-nz,nz
	Execute cmd
	Execute "ModifyGizmo setDisplayList=1, attribute=lineWidthAxesCue"
	Execute "ModifyGizmo setDisplayList=2, object=freeAxesCue0"
	Variable ci,cj,ck
	GrainCOM(grainNos,ci,cj,ck)		// moves the grain Center Of Mass in ci,cj,ck
	sprintf cmd "ModifyGizmo setDisplayList=3, opName=translateSampleCenter, operation=translate, data={%g,%g,%g}",-ci,-cj,-ck
	Execute cmd
	Execute "ModifyGizmo setDisplayList=4, object=light0"
	Execute "ModifyGizmo setDisplayList=5, object=light1"
	Execute "ModifyGizmo setDisplayList=6, attribute=diffuseCursor"
	Execute "ModifyGizmo setDisplayList=7, opName=translateCursor, operation=translate, data={0,0,0}"
	Execute "ModifyGizmo setDisplayList=8, object=sphereCursor"
	Execute "ModifyGizmo setDisplayList=9, opName=translateCursorBack, operation=translate, data={0,0,0}"

	Variable xc=0, yc=0, zc=0
	if (!numtype(firstGrain))
//		Wave grain = $StringFromList(i,gList)	// changed April 1, 2005, JZT
		Wave grain = $StringFromList(0,gList)	// first grain in list is the 'central' grain from firstGrain
		str = note(grain)
		xc=NumberByKey("grain_xc",str,"=")
		yc=NumberByKey("grain_yc",str,"=")
		zc=NumberByKey("grain_zc",str,"=")
		xc = numtype(xc) ? 0 : xc
		yc = numtype(yc) ? 0 : yc
		zc = numtype(zc) ? 0 : zc
	endif
	sprintf cmd, "ModifyGizmo setDisplayList=10, opName=NormalTranslate, operation=translate, data={%g,%g,%g}",xc,yc,zc
	Execute cmd
	Execute "ModifyGizmo setDisplayList=11, attribute=diffuseNormal"
	Execute "ModifyGizmo setDisplayList=12, object=sphereNormal"
	Execute "ModifyGizmo setDisplayList=13, object=GrainBoundaryPath"
	sprintf cmd, "ModifyGizmo setDisplayList=14, opName=NormalTranslateBack, operation=translate, data={%g,%g,%g}",-xc,-yc,-zc
	Execute cmd
	Execute "ModifyGizmo setDisplayList=15, attribute=specularGrains"
	Execute "ModifyGizmo setDisplayList=16, attribute=blendFunc0"
	Execute "ModifyGizmo setDisplayList=17, opName=translate111, operation=translate, data={1,1,1}"
	Execute "ModifyGizmo setDisplayList=18, object=scatterBndry"

//	Execute "ModifyGizmo setDisplayList=17, opName=enableLineSmooth, operation=enable, data=2848"
//	Execute "ModifyGizmo setDisplayList=18, opName=enable0, operation=enable, data=3042"
	noteStr = ReplaceNumberByKey("ortho0",noteStr,nz2,"=",":")
	noteStr = ReplaceNumberByKey("ci",noteStr,ci,"=",":")
	noteStr = ReplaceNumberByKey("cj",noteStr,cj,"=",":")
	noteStr = ReplaceNumberByKey("ck",noteStr,ck,"=",":")
	noteStr = ReplaceNumberByKey("xc",noteStr,DimDelta(grainNos,0)*ci + DimOffset(grainNos,0),"=",":")
	noteStr = ReplaceNumberByKey("yc",noteStr,DimDelta(grainNos,1)*cj + DimOffset(grainNos,1),"=",":")
	noteStr = ReplaceNumberByKey("zc",noteStr,DimDelta(grainNos,2)*ck + DimOffset(grainNos,2),"=",":")
	// for firstGrain=NaN, this loops just runs indicies down from (N-1) to 0
	// for firstGrain=number, then it does index 0 first, then backwards from (N-1) to 1
	for (i=0;i<N;i+=1)
		// this strange thing, 
		grainName= StringFromList(mod(N-1-i+!numtype(firstGrain),N),displayList)
		if (strlen(grainName)<1)
			break
		endif
		sprintf cmd, "ModifyGizmo setDisplayList=-1, object=%s", grainName
		Execute cmd
	endfor
//		for (i=0;i<N;i+=1)
//			grainName= StringFromList(N-1-i,displayList)
//			if (strlen(grainName)<1)
//				break
//			endif
//			sprintf cmd, "ModifyGizmo setDisplayList=-1, object=%s", grainName
//			Execute cmd
//		endfor

	if (showScaleCube)
		sprintf cmd, "ModifyGizmo setDisplayList=-1, opName=scaleBoxtranslate, operation=translate, data={%g,%g,%g}",2*nx,-2*ny,nz/4
		Execute cmd
		Execute "ModifyGizmo setDisplayList=-1, attribute=scaleBoxspecular"
		Execute "ModifyGizmo setDisplayList=-1, attribute=scaleBoxambient"
		Execute "ModifyGizmo setDisplayList=-1, attribute=scaleBoxdiffuse"
		Execute "ModifyGizmo setDisplayList=-1, attribute=scaleBoxcolor"
		Execute "ModifyGizmo setDisplayList=-1, object=scaleBox"

		Execute "ModifyGizmo setDisplayList=-1, attribute=scaleStringColor"
		sprintf cmd, "ModifyGizmo setDisplayList=-1, opName=translateScaleBoxText, operation=translate, data={0,-%g,0}",scaleCubeSize
		Execute cmd
		Execute "ModifyGizmo setDisplayList=-1, opName=scaleStringRotate, operation=rotate, data={-90,0,1,0}"
		sprintf cmd,"ModifyGizmo setDisplayList=-1, opName=scaleStringscale, operation=scale, data={%g,%g,%g}",scaleCubeSize,scaleCubeSize,scaleCubeSize
		Execute cmd
		Execute "ModifyGizmo setDisplayList=-1, object=scaleString"
	endif
//	Execute "ModifyGizmo SETQUATERNION={-0.813679,0.072385,-0.562105,-0.129330}"
//	Execute "ModifyGizmo SETQUATERNION={0,0.707107,0,0.707107}"
	Execute "ModifyGizmo SETQUATERNION={-0.07,0.6,0.06,0.78}"
	Execute "ModifyGizmo autoscaling=1"
	Execute "ModifyGizmo compile"
//	Execute "ModifyGizmo showInfo"
//	Execute "ModifyGizmo infoWindow={504,264,999,529}"
	Execute "ModifyGizmo bringToFront"
//	Execute "ModifyGizmo showAxisCue=1"			// the axis is too small to be seen, so don't bother putting it up
	Execute "ModifyGizmo endRecMacro"
	SetGizmoNote(noteStr)
	Execute "PanelGizmoView()"
//	if (strlen(WinList("Panel_Gizmo_Viewer","","WIN:64"))>1)
//		DoWindow/F Panel_Gizmo_Viewer
//	else
//		Execute "PanelGizmoView()"
//	endif
	return 0
End


Function AddIsoSurfaceGrainToObjectList(description)
	String description

	String name = StringByKey("name",description,"=")
	String wname =  StringByKey("WaveName",description,"=")
	String RGBA = StringByKey("RGBA",description,"=")
	Wave grain = $wName
	Variable isoValue
	String cmd
//	isoValue = (NumberByKey("lo",note(grain),"=")+NumberByKey("hi",note(grain),"="))/4 - .5
	isoValue = NumberByKey("isoValue", note(grain),"=")
	if (numtype(isoValue))
//		isoValue = (NumberByKey("grainNum",note(grain),"=")+1)/4 - 1
		isoValue = NumberByKey("grainNum",note(grain),"=")/2 - 0.5
//		isoValue = (NumberByKey("grainNum",note(grain),"=")+1)*0.9 - 1
	endif

	// write these values to the wave note
	String noteStr=note(grain)
	noteStr = ReplaceNumberByKey("red",noteStr,str2num(StringFromList(0,RGBA,",")),"=")
	noteStr = ReplaceNumberByKey("green",noteStr,str2num(StringFromList(1,RGBA,",")),"=")
	noteStr = ReplaceNumberByKey("blue",noteStr,str2num(StringFromList(2,RGBA,",")),"=")
	noteStr = ReplaceNumberByKey("alpha",noteStr,str2num(StringFromList(3,RGBA,",")),"=")
	noteStr = ReplaceNumberByKey("isoValue",noteStr,isoValue,"=")
	Note/K grain
	Note grain, noteStr

	Execute "AppendToGizmo isoSurface="+wName+",name="+name
	Execute "ModifyGizmo ModifyObject="+name+" property={surfaceColorType,1}"
	Execute "ModifyGizmo ModifyObject="+name+" property={lineColorType,0}"
	Execute "ModifyGizmo ModifyObject="+name+" property={lineWidthType,0}"
	Execute "ModifyGizmo ModifyObject="+name+" property={fillMode,2}"
	Execute "ModifyGizmo ModifyObject="+name+" property={lineWidth,1}"
	sprintf cmd, "ModifyGizmo ModifyObject=%s property={isoValue,%g}",name,isoValue
	Execute cmd
//	Execute "ModifyGizmo ModifyObject="+name+" property={frontColor,"+RGBA+"}"
	Execute "ModifyGizmo ModifyObject="+name+" property={frontColor,1,1,1,0}"		// a volume surrounded by small values
	Execute "ModifyGizmo ModifyObject="+name+" property={backColor,"+RGBA+"}"
	Execute "ModifyGizmo modifyObject="+name+" property={calcNormals,1}"
	Execute "ModifyGizmo modifyObject="+name+" property={useCube,0}"			// 1 will be faster, but not right
End




Function GizmoAllRodriquesVectors(grains)
	Wave grains								// 3d wave with all of the grain numbers (e.g. gNeighbor)
	// Do nothing if the Gizmo XOP is not available.
	if(exists("NewGizmo")!=4)
		DoAlert 0, "Gizmo XOP must be installed"
		return 1
	endif
	if (exists(":raw:axisXYZ")!=1)
		DoAlert 0, "cannot find ':raw:axisXYZ' or ':raw:axisXYZ_Colors'"
		return 1
	endif
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

	String str, cmd
	str = StringByKey("SCREEN1",IgorInfo(0))
	str = str[strsearch(str,"RECT=",0),Inf]
	Variable size = 0.85*str2num(StringFromList(3,str[strsearch(str,"RECT=",0),Inf],","))
	size = round(size*0.6)
	sprintf cmd, "NewGizmo/N=RodriquesGizmo/W=(%d,%d,%d,%d)/K=1",90,90,90+size,90+size
	Execute cmd
	Execute "ModifyGizmo startRecMacro"
	Execute "ModifyGizmo outputResFactor=2"		// n=2, default

	Execute "AppendToGizmo Scatter=:raw:axisXYZ,name=scatterAll"
	Execute "ModifyGizmo ModifyObject=scatterAll property={ scatterColorType,1}"
	Execute "ModifyGizmo ModifyObject=scatterAll property={ markerType,0}"
	Execute "ModifyGizmo ModifyObject=scatterAll property={ sizeType,0}"
	Execute "ModifyGizmo ModifyObject=scatterAll property={ rotationType,0}"
	Execute "ModifyGizmo ModifyObject=scatterAll property={ Shape,2}"
	Execute "ModifyGizmo ModifyObject=scatterAll property={ size,0.3}"
	Execute "ModifyGizmo ModifyObject=scatterAll property={ colorWave,:raw:axisXYZ_Colors}"
	sprintf cmd, "AppendToGizmo Scatter=%s_Rodrigues,name=scatterGrains",GetWavesDataFolder(grains,2)
	Execute cmd
	Execute "ModifyGizmo ModifyObject=scatterGrains property={ scatterColorType,1}"
	Execute "ModifyGizmo ModifyObject=scatterGrains property={ markerType,0}"
	Execute "ModifyGizmo ModifyObject=scatterGrains property={ sizeType,0}"
	Execute "ModifyGizmo ModifyObject=scatterGrains property={ rotationType,0}"
	Execute "ModifyGizmo ModifyObject=scatterGrains property={ Shape,14}"
	Execute "ModifyGizmo ModifyObject=scatterGrains property={ size,1}"
	sprintf cmd,"ModifyGizmo ModifyObject=scatterGrains property={ colorWave,%s_Colors}",GetWavesDataFolder(grains,2)
	Execute cmd
	Execute "AppendToGizmo light=Directional,name=light0"
	Execute "ModifyGizmo light=light0 property={ position,0.0,0.0,-1.0,0.0}"
	Execute "ModifyGizmo light=light0 property={ direction,0.0,0.0,-1.0}"
	Execute "ModifyGizmo light=light0 property={ ambient,0.866667,0.866667,0.866667,1.0}"
	Execute "ModifyGizmo light=light0 property={ specular,0.266667,0.266667,0.266667,1.0}"
	Execute "ModifyGizmo light=light0 property={ diffuse,0.266667,0.266667,0.266667,1.0}"
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
	Execute "ModifyGizmo ModifyObject=axes0,property={0,axisType,4194305}"
	Execute "ModifyGizmo ModifyObject=axes0,property={1,axisType,4194306}"
	Execute "ModifyGizmo ModifyObject=axes0,property={2,axisType,4194308}"
	Execute "ModifyGizmo ModifyObject=axes0,property={3,axisType,4194312}"
	Execute "ModifyGizmo ModifyObject=axes0,property={4,axisType,4194320}"
	Execute "ModifyGizmo ModifyObject=axes0,property={5,axisType,4194336}"
	Execute "ModifyGizmo ModifyObject=axes0,property={6,axisType,4194368}"
	Execute "ModifyGizmo ModifyObject=axes0,property={7,axisType,4194432}"
	Execute "ModifyGizmo ModifyObject=axes0,property={8,axisType,4194560}"
	Execute "ModifyGizmo ModifyObject=axes0,property={9,axisType,4194816}"
	Execute "ModifyGizmo ModifyObject=axes0,property={10,axisType,4195328}"
	Execute "ModifyGizmo ModifyObject=axes0,property={11,axisType,4196352}"
	Execute "ModifyGizmo ModifyObject=axes0,property={-1,axisScalingMode,1}"
	Execute "ModifyGizmo ModifyObject=axes0,property={-1,axisColor,0,0,0,1}"
	Execute "ModifyGizmo setDisplayList=0, object=light0"
	Execute "ModifyGizmo setDisplayList=1, object=scatterAll"
	Execute "ModifyGizmo setDisplayList=2, object=scatterGrains"
	Execute "ModifyGizmo setDisplayList=3, object=axes0"
	//	Execute "ModifyGizmo SETQUATERNION={-0.082902,-0.728104,-0.095385,0.673720}"
	//	Execute "ModifyGizmo SETQUATERNION={-0.042222,0.662365,0.033769,0.747230}"
	Execute "ModifyGizmo SETQUATERNION={0.5,0.4,0.5,0.6}"
	Execute "ModifyGizmo autoscaling=1"
	Execute "ModifyGizmo currentGroupObject=\"\""
	Execute "ModifyGizmo compile"

	Execute "ModifyGizmo showAxisCue=1"
	Execute "ModifyGizmo endRecMacro"
End


//	=========================================================================
//	=========================================================================
//	=========================================================================




