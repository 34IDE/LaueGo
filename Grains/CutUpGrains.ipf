#pragma rtGlobals=1		// Use modern global access method.
#pragma version = 1.01
#include "GrainMath", version>=1.0



Menu "Grains"
	"-"
	"(2d Cuts"
	"Graph a cut array", GraphCut($"")
	"Graph All cut arrays",Graph_All_Cuts()
	"make all cuts",MakeAllCuts($"")
	"  make one cut",oneCut($""<NaN)
End



//  ===================================================================
//  ===================================================================
//  ===================================================================
//		cut up and display grains[][][]
//  ===================================================================

Function MakeAllCuts(grains)
	Wave grains			// 3d wave with all of the grain numbers (e.g. gNeighbor)
	if (!WaveExists(grains))
		String grainName
		Prompt grainName,"3d wave with all grains",popup,ListMultiGrainMats()
		DoPrompt "grain array to take cuts from",grainName
		if (V_flag)
			return 1
		endif
		Wave grains = $grainName
	endif
	if (!WaveExists(grains) || WaveDims(grains)!=3)
		DoAlert 0, "in MakeAllCuts(), the wave 'grains[][][]' does not exist or has wrong dimension"
		return 1
	endif

	Variable Ny = DimSize(grains,1)
	Variable i
	String str
	for (i=0;i<Ny;i+=1)
		str = oneCut(grains,i)
		if (strlen(str)<1)
			DoAlert 0, "failed to make cut"+num2istr(i)
		endif
	endfor
	makeGrainsColorTable(grains)
	return 0
End
Function/S oneCut(grains,column)
	Wave grains					// 3d wave with all of the grain numbers (e.g. gNeighbor)
	Variable column					// column index to select
	String str
	if (!WaveExists(grains) || WaveDims(grains)!=3 || numtype(column) || column<0)	// prompt for column
		String grainName=NameOfWave(grains)
		Prompt grainName,"3d wave with all grains",popup,ListMultiGrainMats()
		column = (numtype(column) || column<0) ? 0 : column
//		sprintf str, "colum (y) in grains[%d][%d][%d] to cut must be in range [0,%d]",Nx,Ny,Nz,Ny-1
		Prompt column,"colum (y) in to take cut from"
		DoPrompt "grain array to take cut from & column",grainName,column
		if (V_flag)
			return ""
		endif
		Wave grains = $grainName
	endif
	if (!WaveExists(grains) || WaveDims(grains)!=3)
		DoAlert 0, "in oneCut(), the wave 'grains[][][]' does not exist or has wrong dimension"
		return ""
	endif
	Variable Nx = DimSize(grains,0), Ny = DimSize(grains,1), Nz = DimSize(grains,2)
	Variable calledDirectly = (ItemsInList(GetRTStackInfo(0))<2)
	if (numtype(column) || column<0 || column>=Ny)
		if (calledDirectly)
			sprintf str, "column %d of '%s'  is invalid, it must be in [0,%d]",column,NameOfWave(grains),Ny-1
			DoAlert 0, str
		endif
		return ""
	endif

	String gName = GetWavesDataFolder(grains,2)
	sprintf str,"%s_cut%d",gName,column
	Make/N=(Nz,Nx)/O $str
	Wave cut = $str
	cut = grains[q][column][p]
	String micron = StrVarOrDefault("root:Packages:Grains:micron","µm")
	SetScale/P x DimOffset(grains,2),DimDelta(grains,2),micron, cut
	SetScale/P y DimOffset(grains,0),DimDelta(grains,0),micron, cut

	String noteStr = note(grains)
	noteStr = ReplaceNumberByKey("iycut",noteStr,column,"=")
	noteStr = ReplaceNumberByKey("ycut",noteStr,column*DimDelta(grains,1)+DimOffset(grains,1),"=")
	Note/K cut
	Note cut,noteStr
	if (calledDirectly)
		printf "made cut:   %s[%d][%d] from column %d of %s[%d][%d][%d]\r",NameOfWave(cut),Nz,Nx,column,NameofWave(grains),Nx,Ny,Nz
	endif
	return GetWavesDataFolder(cut,2)
End
Function GraphCut(cut) : Graph
	Wave cut
	Variable calledFromTop = (ItemsInList(GetRTStackInfo(0))<2 )

	if ((!WaveExists(cut) || WaveDims(cut)!=2) && calledFromTop)
		String cutName
		Prompt cutName,"2d cut array",popup,WaveList("*_cut*",";","DIMS:2")
		DoPrompt "cut array", cutName
		if (V_flag)
			return 1
		endif
		Wave cut = $cutName
	endif
	if ((!WaveExists(cut) || WaveDims(cut)!=2))
		DoAlert 0, "in GraphCut(), the wave 'cut[][]' does not exist or has wrong dimension"
		return 1
	endif

	String grName="Graph_"+NameOfWave(cut)
	if (strlen(WinList(grName,"","WIN:1"))>2)
		DoWindow/F $grName
		return 0
	endif
	Display /W=(61,44,991,166)
	DoWindow/C $grName
	AppendImage cut
	ModifyImage $NameOfWave(cut) ctab= {0,300,Fiddle,0}
	ModifyGraph mirror=2, minor=1
	String str
	String micron = StrVarOrDefault("root:Packages:Grains:micron","µm")
	sprintf str, "%s\rcol %g\rY=%g%s",NameOfWave(cut),NumberByKey("iycut",note(cut),"="),NumberByKey("iycut",note(cut),"="),micron
	TextBox/N=text0/F=0/S=3/A=LT/X=0.69/Y=1.27 str
	return 0
End



Function Graph_All_Cuts(grains) : Graph
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
	String graphName = "Graph_All_Cuts_"+NameOfWave(grains)

	if (strlen(WinList(graphName,"","WIN:1")))
		DoWindow/F $graphName
		return 1
	endif
	Display /W=(113,146,792,781)
	DoWindow/C $graphName
	Variable grainMax = numpnts(:Anneal350C:gNeighbor_Vol)-1
	String cut, axis, cmd
	Variable i=0
	for (i=0;i<10;i+=1)
		cut = "cut"+num2istr(i)
		axis = SelectString(i,"left",cut)
		Wave cutWave = $(gName+"_"+cut)
		AppendImage/L=$axis cutWave
		ModifyImage $NameOfWave(cutWave) cindex= $(gName+"_Colors")
		ModifyGraph freePos($axis)=0
		ModifyGraph axisEnab($axis)={(9-i)/10+.005,(9-i+1)/10 - 0.005}
		Label $axis "\\e"
	endfor
	ModifyGraph gfSize=12
	ModifyGraph mirror=2,mirror=2
	ModifyGraph mirror(bottom)=1
	ModifyGraph lblPos(left)=44, lblLatPos(left)=453
	Label left "X  (\\U)"
	Label bottom "Z  (\\U)"
	ShowInfo
	TextBox/N=text0/F=0/S=3/A=LT/X=0.81/Y=90.13 "\\Z09cut9\rcol 9\rY=9µm"
	SetDrawLayer UserFront
	i = 1
	do
		SetDrawEnv linethick= 2
		DrawLine 0,i/10,1,i/10
		i += 1
	while(i<10)
End


