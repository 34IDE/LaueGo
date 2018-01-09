#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName=vti
#pragma IgorVersion = 6.3
#pragma version = 0.10
#include "ColorNames"

Menu "VTI"
	"Load VTI file...",LoadVTIfile("")
	"Display VTI Gizmo",MakeGizmoIsoVolume($"")
End

Menu "Load Waves"
	"Load VTI file...",LoadVTIfile("")
End

Static StrConstant VTIfilters = "VTI Files (*.vti):.vti,;All Files:.*;"



//  ======================================================================================  //
//  ============================== Start of VTI file loading =============================  //

//	http://www.cacr.caltech.edu/~slombey/asci/vtk/vtk_formats.simple.html

Function/WAVE LoadVTIfile(fileNameFull,[printIt])
	String fileNameFull
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt

	Variable f=0
	Open/R/P=home/F=VTIfilters/M="VTI files"/Z=2 f as fileNameFull
	if (strlen(S_fileName)<1 || !f)
		return $""
	endif
	fileNameFull = S_fileName
	Variable tick0=stopMSTimer(-2)
	FStatus f
	String vti = PadString("",V_logEOF,0x20)
	FBinRead f, vti
	Close f
	String wName = CleanupName(ParseFilePath(3,fileNameFull,":",0,0),0)

	if (StringMatch(wName,"*_ascii") && strlen(wName)>6)
		wName = wName[0, strlen(wName)-7]
	endif

	if (Exists(wName))
		wName = UniqueName(wName,1,1)
	endif
	if (printIt)
		printf "read from \"%s\", %d bytes long\r",fileNameFull,strlen(vti)
	endif

	Variable istart=strsearch(vti,"<DataArray ",0), line0=0, i=0, numLines
	if (istart<0)
		return $""
	endif
	do 					// find line0, the line number where data starts
		i = strsearch(vti,"\n",i+1)
		line0 += 1
	while(i<istart)

	istart = i			// find the number of data lines to read, start where we just ended
	i = strsearch(vti,"</DataArray>",istart)	// find the end
	String str = vti[istart,i]
	numLines = strlen(str)
	str = ReplaceString("\n",str,"")				// remove all new line characters
	numLines -= strlen(str)+1							// number of data lines

	str = XMLattibutes2KeyList("ImageData",vti)
	Wave xyz0 = str2vec(StringByKey("Origin",str,"="))		// starting point
	Wave dxyz = str2vec(StringByKey("Spacing",str,"="))		// increment

	str = XMLattibutes2KeyList("Piece",vti)
	Wave Ns = str2vec(StringByKey("Extent",str,"="))			// array size

	str = XMLattibutes2KeyList("DataArray",vti)
	String format = StringByKey("format",str,"=")				// format must be "ascii"
	if (!StringMatch(format,"ascii"))
		sprintf str, "This routine only understands \"ascii\" format, not \"%s\"", format
		if (printIt)
			DoAlert 0, str
			print str
		endif
		return $""
	endif

	LoadWave/A/G/H/K=1/M/L={line0-1,line0,numLines,0,0}/P=home/Q fileNameFull
	if (V_flag!=1)
		if (printIt)
			printf "ERROR -- in LoadWave of \"%s\",  V_flag = %g\r",S_waveNames, V_flag
		endif
		return $""
	endif
	Wave wav = $StringFromList(0,S_waveNames)

	Variable Nx, Ny=-1, Nz=-1, N						// -1 causes dimension to be unchanged in the Redimension
	if (WaveExists(Ns))
		Nx = Ns[1]-Ns[0] + 1
		Ny = DimSize(Ns,0)>2 ? Ns[3]-Ns[2] + 1 : Ny
		Nz = DimSize(Ns,0)>4 ? Ns[5]-Ns[4] + 1 : Nz
		N = abs(Nx * Ny * Nz)							// abs needed in case Nz=-1
	else
		Nx = numpnts(wav)
		N = Nx
	endif
	if (printIt)
		printf "Loaded wave:  %s[%d][%d][%d],  (%d total points), range=[%g, %g],  took %.2f sec\r",wName,Nx,Ny,Nz,N, WaveMin(wav),WaveMax(wav), (stopMSTimer(-2)-tick0)*1e-6
	endif

	Rename wav, $wName
	MatrixOP/O wav = wav^t								// transpose in place
	Redimension/N=(N) wav								// make 1D
	if (Ny>0)
		Redimension/N=(Nx,Ny,Nz) wav					// make correct 3D shape
	endif

	if (WaveExists(xyz0) && WaveExists(dxyz))
		SetScale/P x xyz0[0],dxyz[0],"", wav
		if (numpnts(xyz0)>=2 && Ny>0)
			SetScale/P y xyz0[1],dxyz[1],"", wav
			if (numpnts(xyz0)>=3 && Nz>0)
				SetScale/P z xyz0[2],dxyz[2],"", wav
			endif
		endif
	endif

	String wnote="waveClass=vtiPointData;"
	wnote = ReplaceStringByKey("fullFile",wnote,fileNameFull,"=")
	Note/K wav, wnote
	return wav
End

	//		# Alternative python code
	//		import vtk
	//		from vtk.util.numpy_support import vtk_to_numpy
	//
	//		import sys
	//		import numpy as np
	//
	//		if sys.argv < 3:
	//			print "Usage python vti2txt inputVtiName outputTxtName"
	//			exit(0)
	//		filename = sys.argv[1]
	//		outfilename = sys.argv[2]
	//
	//		reader = vtk.vtkXMLImageDataReader()
	//		reader.SetFileName(filename)
	//		reader.Update()
	//
	//		output = reader.GetOutput()
	//
	//		dataOrigin = output.GetOrigin()
	//		dataSpacing = output.GetSpacing()
	//		dataExtent = output.GetExtent()
	//
	//		vtkData = output.GetPointData().GetScalars()
	//
	//		npData = vtk_to_numpy(vtkData)
	//
	//		header = "origin: " + str(dataOrigin) + "\nspacing: " + str(dataSpacing) + "\nextent: " + str(dataExtent)
	//		np.savetxt(outfilename, npData, header=header)

//  =============================== End of VTI file loading ==============================  //
//  ======================================================================================  //




//  ======================================================================================  //
//  =============================== Start of Make the Gizmo ==============================  //

//	Make/FREE isoValues = {0.2, 0.1, 0.01, 0.002, 0.001}
//	isoValues *= WaveMax(volWave)
//	Make/FREE/T isoColors = {"0.5,0,0,1","0.75,0.25,0,.8","0.75,0.25,0,0.5","0.75,0.25,0,0.3","0.75,0.25,0,0.05"}

Function MakeGizmoIsoVolume(volWave,[isoMin,isoMax,Niso,ColorTable,revColors,isoValues,isoColors])
	Wave volWave
	// The following are all OPTIONAL
	Variable isoMin					// used to set smallest iso value
	Variable isoMax					// used to set biggest iso value
	Variable Niso						// number of requested iso levels
	String ColorTable					// color table to use for setting iso levels
	Variable revColors				// optionaly reverse the color table
	Wave isoValues						// user supplied values for iso surfaces
	Wave isoColors						// user supplied rgba for iso surfaces
	isoMax = ParamIsDefault(isoMax) ? NaN : isoMax
	isoMin = ParamIsDefault(isoMin) ? NaN : isoMin
	Niso = ParamIsDefault(Niso) ? NaN : Niso
	if (ParamIsDefault(ColorTable))
		ColorTable = ""
	endif
	revColors = ParamIsDefault(revColors) ? NaN : revColors
	if(exists("NewGizmo")!=4)			// Do nothing if the Gizmo XOP is not available.
		DoAlert 0, "Gizmo XOP must be installed"
		return 1
	endif

	Variable printIt=0
	if (!WaveExists(volWave))
		String wList=WaveListClass("vtiPointData*","*","DIMS:3"), wName=""
		if (ItemsInList(wList)<2)
			wName = StringFromList(0,wList)
		else
			wName = StrVarOrDefault("root:Packages:QspaceVolumes:Gizmo:volWaveName","")
			Prompt wName,"3D Qspace",popup,wList
			DoPrompt "Pick 3D Qspace",wName
			if (V_flag)
				return 1
			endif
		endif
		Wave volWave = $wName
		printIt = 1
	endif
	if (!WaveExists(volWave))
		DoAlert 0, "Cannot find 3D Qspace Wave"
		return 1
	endif
	String gizName=StringFromlist(0,WindowsWithWave(volWave,4))	// volWave already on a Gizmo, just bring it to front
	if (strlen(gizName)>0)
		DoWindow/F $gizName
		return 0
	endif

	gizName = "Gizmo"+CleanupName(NameOfWave(volWave),0)
	Variable QxLo=NaN,QxHi=NaN, QyLo=NaN,QyHi=NaN, QzLo=NaN,QzHi=NaN
	String str
	Variable i, isoMaxInternal=NaN, isoMinInternal=NaN
	if (WaveExists(isoValues) && WaveExists(isoColors))	// use the user supplied iso values and colors
		if (numpnts(isoValues) != DimSize(isoColors,0))
			DoAlert 0,"The user supplied isoValues[ ] and isoColors[ ] waves have different sizes."
			return 1
		endif
//		for (i=0;i<numpnts(isoValues);i+=1)
//			sprintf str,"%g:%s;", isoValues[i],isoColors[i]
//			isoStr += str
//		endfor
	else
		if (numtype(isoMax))
			// isoMaxInternal = FindIsoLevelByFraction(volWave,0.0005)
			isoMaxInternal = FindIsoLevelByFraction(volWave,numpnts(volWave)^-0.6)
			isoMaxInternal = roundSignificant(isoMaxInternal,2)
			isoMax = numtype(isoMax) ? NumVarOrDefault("root:Packages:QspaceVolumes:Gizmo:isoMax",NaN) : isoMax
			isoMax = numtype(isoMax) ? isoMaxInternal : isoMax
		endif
		if (numtype(isoMin))
			isoMinInternal = FindIsoLevelByFraction(volWave,0.25)
			isoMinInternal = roundSignificant(isoMinInternal,2)
			isoMin = numtype(isoMin) ? NumVarOrDefault("root:Packages:QspaceVolumes:Gizmo:isoMin",NaN) : isoMin
			isoMin = numtype(isoMin) ? isoMinInternal : isoMin
		endif
		if (WhichListItem(ColorTable,CTabList())<0)
			ColorTable = StrVarOrDefault("root:Packages:QspaceVolumes:Gizmo:ColorTable","Rainbow256")
		endif
		Niso = (Niso>0 &&Niso<50) ? Niso : NumVarOrDefault("root:Packages:QspaceVolumes:Gizmo:Niso",5)
		revColors = numtype(revColors) ? NumVarOrDefault("root:Packages:QspaceVolumes:Gizmo:revColors",0) : !(!revColors)

		Prompt ColorTable,"Color Table",popup,CTabList()
		Prompt Niso,"Number of Iso Levels",popup,"1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;"
		Prompt isoMax,"Max value for an Iso Level"
		Prompt isoMin,"Min value for an Iso Level"
		Prompt revColors,"Reverse Colors",popup,"Direct Colors;ReverseColors"
		revColors += 1
		DoPrompt "Iso Levels & Colors",ColorTable,revColors,Niso,isoMax,isoMin
		if (V_flag)
			return 1
		endif
		revColors -= 1
		isoMax = numtype(isoMax) ? isoMaxInternal : isoMax
		isoMin = numtype(isoMin) ? isoMinInternal : isoMin
		if (!(Niso>0 &&Niso<50)  ||  WhichListItem(ColorTable,CTabList())<0 || numtype(isoMax) || numtype(isoMin))
			return 1
		endif
		Wave isoValues = EqualVolumeIsoLevels(volWave,Niso,isoMin,isoMax)
		Wave isoColors = IsoColorsWave(Niso,ColorTable,revColors)		// returns a wave with rgba for each of N levels
		printIt = 1
	endif

	if (printIt)
		printf "MakeGizmoIsoVolume(%s",NameOfWave(volWave)
		if (WaveExists(isoValues) && WaveExists(isoColors))	// use the user supplied iso values and colors
			printf ", isoValues=%s,  isoColors=%s",NameOfWave(isoValues),NameOfWave(isoColors)
		else
			printf ", isoMax=%g, isoMin=%g, Niso=%g, ColorTable=\"%s\", revColors=%g",isoMax,isoMin,Niso,ColorTable,revColors
		endif
		printf ")\r"
	endif

	if (numtype(sum(isoValues)+sum(isoColors)))
		return 1
	endif
//	String cornerList = WaveListClass("GizmoCorners","*","DIMS:2,MINROWS:2,MAXROWS:2,MINCOLS:3,MAXCOLS:3"), cName=""
//	cornerList += WaveListClass("GizmoCorners","*","DIMS:2,MINROWS:8,MAXROWS:8,MINCOLS:3,MAXCOLS:3")
//	cornerList = WavesWithMatchingKeyVals(cornerList,"sourceWave="+NameOfWave(volWave))
//	if (ItemsInList(cornerList)==1)
//		cName = StringFromList(0,cornerList)
//	elseif (ItemsInList(cornerList)>1)
//		cName = StrVarOrDefault("root:Packages:QspaceVolumes:Gizmo:cornersName","")
//		Prompt cName,"3D Qspace CORNERS",popup,cornerList
//		DoPrompt "Pick 3D Qspace CORNERS",cName
//		if (V_flag)
//			return 1
//		endif
//	endif
//
//	Wave corners = $cName

	NewDataFolder/O root:Packages
	NewDataFolder/O root:Packages:QspaceVolumes
	NewDataFolder/O root:Packages:QspaceVolumes:Gizmo
	String/G root:Packages:QspaceVolumes:Gizmo:volWaveName = NameOfWave(volWave)
//	if (WaveExists(corners))
//		String/G root:Packages:QspaceVolumes:Gizmo:cornersName = NameOfWave(corners)
//	endif
	Variable/G root:Packages:QspaceVolumes:Gizmo:isoMax = (isoMax==isoMaxInternal) ? NaN : isoMax
	Variable/G root:Packages:QspaceVolumes:Gizmo:isoMin = (isoMin==isoMinInternal) ? NaN : isoMin
	Variable/G root:Packages:QspaceVolumes:Gizmo:Niso = Niso
	String/G root:Packages:QspaceVolumes:Gizmo:ColorTable = ColorTable
	Variable/G root:Packages:QspaceVolumes:Gizmo:revColors = revColors

#if (IgorVersion()<7)
	Execute "NewGizmo/N="+gizName+"/T=\""+gizName+"\" /W=(154,44,922,812)"
	Execute "ModifyGizmo startRecMacro"
	Execute "AppendToGizmo attribute blendFunc={770,771},name=blendingFunction"
#else
	Variable Billboarding = 1
	NewGizmo/N=$gizName/T=gizName/W=(154,44,922,812)
	ModifyGizmo startRecMacro
	AppendToGizmo attribute blendFunc={770,771},name=blendingFunction
#endif
	String object,displayObjectList=""

//		if (WaveExists(corners))
//			displayObjectList += "scatterCubeCorners;"
//			Variable icor=DimSize(corners,0)-1			// this will be 2 or 7 (7 is old way)
//			ImageStats/M=1/G={0,icor, 0,0} corners;		QxLo = V_min	;	QxHi = V_max	// need these later for setting a*, b*, c*
//			ImageStats/M=1/G={0,icor, 1,1} corners;		QyLo = V_min	;	QyHi = V_max
//			ImageStats/M=1/G={0,icor, 2,2} corners;		QzLo = V_min	;	QzHi = V_max
//	#if (IgorVersion()<7)
//			Execute "AppendToGizmo Scatter="+GetWavesDataFolder(corners,2)+",name=scatterCubeCorners"
//			Execute "ModifyGizmo ModifyObject=scatterCubeCorners property={ scatterColorType,0}"
//			Execute "ModifyGizmo ModifyObject=scatterCubeCorners property={ markerType,0}"
//			Execute "ModifyGizmo ModifyObject=scatterCubeCorners property={ sizeType,0}"
//			Execute "ModifyGizmo ModifyObject=scatterCubeCorners property={ rotationType,0}"
//			Execute "ModifyGizmo ModifyObject=scatterCubeCorners property={ Shape,1}"
//			Execute "ModifyGizmo ModifyObject=scatterCubeCorners property={ size,1}"
//			Execute "ModifyGizmo ModifyObject=scatterCubeCorners property={ color,0,0,0,0}"
//	#else
//			AppendToGizmo Scatter=$GetWavesDataFolder(corners,2),name=scatterCubeCorners
//			ModifyGizmo ModifyObject=scatterCubeCorners objectType=scatter property={ scatterColorType,0}
//			ModifyGizmo ModifyObject=scatterCubeCorners objectType=scatter property={ markerType,0}
//			ModifyGizmo ModifyObject=scatterCubeCorners objectType=scatter property={ sizeType,0}
//			ModifyGizmo ModifyObject=scatterCubeCorners objectType=scatter property={ rotationType,0}
//			ModifyGizmo ModifyObject=scatterCubeCorners objectType=scatter property={ Shape,1}
//			ModifyGizmo ModifyObject=scatterCubeCorners objectType=scatter property={ size,1}
//			ModifyGizmo ModifyObject=scatterCubeCorners objectType=scatter property={ color,0,0,0,0}
//	#endif
//		else
//			QxLo = DimOffset(volWave,0)	;	QxHi = QxLo + (DimSize(volWave,0)-1) * DimDelta(volWave,0)
//			QyLo = DimOffset(volWave,1)	;	QyHi = QyLo + (DimSize(volWave,1)-1) * DimDelta(volWave,1)
//			QzLo = DimOffset(volWave,2)	;	QzHi = QzLo + (DimSize(volWave,2)-1) * DimDelta(volWave,2)
//		endif
	QxLo = DimOffset(volWave,0)	;	QxHi = QxLo + (DimSize(volWave,0)-1) * DimDelta(volWave,0)
	QyLo = DimOffset(volWave,1)	;	QyHi = QyLo + (DimSize(volWave,1)-1) * DimDelta(volWave,1)
	QzLo = DimOffset(volWave,2)	;	QzHi = QzLo + (DimSize(volWave,2)-1) * DimDelta(volWave,2)

	Variable val, r,g,b,a
	printf " Iso Values\t\t\t\tColors  {r,g,b,a}\r"
	for (i=0;i<DimSize(isoValues,0);i+=1)
		val = isoValues[i]
		r = isoColors[i][0]	;	g = isoColors[i][1]	;	b = isoColors[i][2]	;	a = isoColors[i][3]
		sprintf str,"%g,%g,%g,%g", r,g,b,a
		displayObjectList += AddIso2Gizmo(volWave,"isoSurface"+num2istr(i),val,str,"")+";"
		if (printIt)
			printf "%g   \t\t'%s'     {%s}\r",val,RGBA2name(r,g,b,a,1),str
		endif
	endfor

	Variable Qunits = StringMatch(WaveUnits(volWave,0),"*nm")
#if (IgorVersion()<7)
	Execute "AppendToGizmo Axes=boxAxes,name=axesBeamLineQ"
	Execute "ModifyGizmo ModifyObject=axesBeamLineQ,property={-1,axisScalingMode,1}"
	Execute "ModifyGizmo ModifyObject=axesBeamLineQ,property={-1,axisColor,0,0,0,1}"
	Execute "ModifyGizmo ModifyObject=axesBeamLineQ,property={0,ticks,3}"
	Execute "ModifyGizmo ModifyObject=axesBeamLineQ,property={1,ticks,3}"
	Execute "ModifyGizmo ModifyObject=axesBeamLineQ,property={2,ticks,3}"
	Execute "ModifyGizmo ModifyObject=axesBeamLineQ,property={0,fontScaleFactor,0.8}"
	Execute "ModifyGizmo ModifyObject=axesBeamLineQ,property={1,fontScaleFactor,0.8}"
	Execute "ModifyGizmo ModifyObject=axesBeamLineQ,property={2,fontScaleFactor,0.8}"
	Execute "ModifyGizmo ModifyObject=axesBeamLineQ,property={0,axisLabel,1}"
	Execute "ModifyGizmo ModifyObject=axesBeamLineQ,property={1,axisLabel,1}"
	Execute "ModifyGizmo ModifyObject=axesBeamLineQ,property={2,axisLabel,1}"
	if (Qunits)
		Execute "ModifyGizmo ModifyObject=axesBeamLineQ,property={0,axisLabelText,\"Qx  (1/nm)\"}"
		Execute "ModifyGizmo ModifyObject=axesBeamLineQ,property={1,axisLabelText,\"Qy  (1/nm)\"}"
		Execute "ModifyGizmo ModifyObject=axesBeamLineQ,property={2,axisLabelText,\"Qz  (1/nm)\"}"
	endif
	Execute "ModifyGizmo ModifyObject=axesBeamLineQ,property={0,axisLabelCenter,-0.1}"
	Execute "ModifyGizmo ModifyObject=axesBeamLineQ,property={1,axisLabelCenter,-0.1}"
	Execute "ModifyGizmo ModifyObject=axesBeamLineQ,property={2,axisLabelCenter,-0.1}"
	Execute "ModifyGizmo ModifyObject=axesBeamLineQ,property={0,axisLabelDistance,0.05}"
	Execute "ModifyGizmo ModifyObject=axesBeamLineQ,property={1,axisLabelDistance,0.05}"
	Execute "ModifyGizmo ModifyObject=axesBeamLineQ,property={2,axisLabelDistance,0.3}"
	Execute "ModifyGizmo ModifyObject=axesBeamLineQ,property={0,axisLabelScale,0.4}"
	Execute "ModifyGizmo ModifyObject=axesBeamLineQ,property={1,axisLabelScale,0.4}"
	Execute "ModifyGizmo ModifyObject=axesBeamLineQ,property={2,axisLabelScale,0.4}"
	Execute "ModifyGizmo ModifyObject=axesBeamLineQ,property={0,axisLabelRGBA,0,0,0,1}"
	Execute "ModifyGizmo ModifyObject=axesBeamLineQ,property={1,axisLabelRGBA,0,0,0,1}"
	Execute "ModifyGizmo ModifyObject=axesBeamLineQ,property={2,axisLabelRGBA,0,0,0,1}"
	Execute "ModifyGizmo modifyObject=axesBeamLineQ property={Clipped,0}"
#else
	AppendToGizmo Axes=boxAxes,name=axesBeamLineQ
	ModifyGizmo ModifyObject=axesBeamLineQ, objectType=axes property={-1,axisScalingMode,1}
	ModifyGizmo ModifyObject=axesBeamLineQ, objectType=axes property={-1,axisColor,0,0,0,1}
	ModifyGizmo ModifyObject=axesBeamLineQ, objectType=axes property={0,ticks,3}
	ModifyGizmo ModifyObject=axesBeamLineQ, objectType=axes property={1,ticks,3}
	ModifyGizmo ModifyObject=axesBeamLineQ, objectType=axes property={2,ticks,3}
	ModifyGizmo ModifyObject=axesBeamLineQ, objectType=axes property={0,fontScaleFactor,1}
	ModifyGizmo ModifyObject=axesBeamLineQ, objectType=axes property={1,fontScaleFactor,1}
	ModifyGizmo ModifyObject=axesBeamLineQ, objectType=axes property={2,fontScaleFactor,1}
	ModifyGizmo ModifyObject=axesBeamLineQ, objectType=axes property={0,axisLabel,1}
	ModifyGizmo ModifyObject=axesBeamLineQ, objectType=axes property={1,axisLabel,1}
	ModifyGizmo ModifyObject=axesBeamLineQ, objectType=axes property={2,axisLabel,1}
	if (Qunits)
		ModifyGizmo ModifyObject=axesBeamLineQ, objectType=axes property={0,axisLabelText,"Qx  (1/nm)"}
		ModifyGizmo ModifyObject=axesBeamLineQ, objectType=axes property={1,axisLabelText,"Qy  (1/nm)"}
		ModifyGizmo ModifyObject=axesBeamLineQ, objectType=axes property={2,axisLabelText,"Qz  (1/nm)"}
	endif
	ModifyGizmo ModifyObject=axesBeamLineQ, objectType=axes property={0,axisLabelCenter,0}
	ModifyGizmo ModifyObject=axesBeamLineQ, objectType=axes property={1,axisLabelCenter,0}
	ModifyGizmo ModifyObject=axesBeamLineQ, objectType=axes property={2,axisLabelCenter,0}
	ModifyGizmo ModifyObject=axesBeamLineQ, objectType=axes property={0,axisLabelDistance,0}
	ModifyGizmo ModifyObject=axesBeamLineQ, objectType=axes property={1,axisLabelDistance,0}
	ModifyGizmo ModifyObject=axesBeamLineQ, objectType=axes property={2,axisLabelDistance,0}
	ModifyGizmo ModifyObject=axesBeamLineQ, objectType=axes property={0,axisLabelScale,1}
	ModifyGizmo ModifyObject=axesBeamLineQ, objectType=axes property={1,axisLabelScale,1}
	ModifyGizmo ModifyObject=axesBeamLineQ, objectType=axes property={2,axisLabelScale,1}
	ModifyGizmo ModifyObject=axesBeamLineQ, objectType=axes property={0,axisLabelRGBA,0,0,0,1}
	ModifyGizmo ModifyObject=axesBeamLineQ, objectType=axes property={1,axisLabelRGBA,0,0,0,1}
	ModifyGizmo ModifyObject=axesBeamLineQ, objectType=axes property={2,axisLabelRGBA,0,0,0,1}
	ModifyGizmo modifyObject=axesBeamLineQ, objectType=axes property={-1,Clipped,0}
	ModifyGizmo ModifyObject=axesBeamLineQ,objectType=Axes,property={0,labelBillboarding,Billboarding}
	ModifyGizmo ModifyObject=axesBeamLineQ,objectType=Axes,property={1,labelBillboarding,Billboarding}
	ModifyGizmo ModifyObject=axesBeamLineQ,objectType=Axes,property={2,labelBillboarding,Billboarding}
#endif
	displayObjectList += "axesBeamLineQ;"

#if (IgorVersion()<7)
	Execute "ModifyGizmo setDisplayList=-1, opName=ortho0, operation=ortho, data={-2,2,-2,2,-3,3}"
	Execute "ModifyGizmo setDisplayList=-1, opName=scale0, operation=scale, data={1.25,1.25,1.25}"
	Execute "ModifyGizmo setDisplayList=-1, attribute=blendingFunction"
	Execute "ModifyGizmo setDisplayList=-1, opName=enableBlend, operation=enable, data=3042"
#else
	ModifyGizmo setDisplayList=-1, opName=ortho0, operation=ortho, data={-2,2,-2,2,-3,3}
	ModifyGizmo setDisplayList=-1, opName=scale0, operation=scale, data={1.25,1.25,1.25}
	ModifyGizmo setDisplayList=-1, attribute=blendingFunction
	ModifyGizmo setDisplayList=-1, opName=enableBlend, operation=enable, data=3042
#endif

	for (i=0;i<ItemsInList(displayObjectList);i+=1)
		object = StringFromList(i,displayObjectList)
		if (strlen(object))
#if (IgorVersion()<7)
			Execute "ModifyGizmo setDisplayList=-1, object="+object
#else
			ModifyGizmo setDisplayList=-1, object=$object
#endif
		endif
	endfor

#if (IgorVersion()<7)
//	Execute "ModifyGizmo SETQUATERNION={0.628289,0.154694,-0.151222,0.747298}"
	Execute "ModifyGizmo autoscaling=1"
	Execute "ModifyGizmo aspectRatio=1"
	Execute "ModifyGizmo currentGroupObject=\"\""
	Execute "ModifyGizmo compile"
	Execute "ModifyGizmo bringToFront"
	Execute "ModifyGizmo endRecMacro"
#else
//	ModifyGizmo SETQUATERNION={0.628289,0.154694,-0.151222,0.747298}
	ModifyGizmo autoscaling=1
	ModifyGizmo aspectRatio=1
	ModifyGizmo currentGroupObject=""
	ModifyGizmo compile
	ModifyGizmo bringToFront
	ModifyGizmo endRecMacro
#endif
	return 0
End
//
Static Function/T AddIso2Gizmo(isoWave,isoName,isoValue,front,back)
	Wave isoWave				// 3d wave
	String isoName				// probably astar, bstar, or cstar, or maybe h,k,l
	Variable isoValue
	String front,back			// if back is empty, make back the complement of front

	if (!WaveExists(isoWave) || strlen(isoName)<1 || numtype(isoValue))
		return ""
	endif
	Variable r,g,b,a
	sscanf front,"%g,%g,%g,%g",r,g,b,a
	if (V_flag!=4)
		return ""
	endif

	Variable Mr, Mg, Mb, Ma
	sscanf back,"%g,%g,%g,%g",Mr,Mg,Mb,Ma
	if (V_flag!=4)
		Mr=1-r; Mg=1-g; Mb=1-b; Ma=a		// use complement color for back side
		sprintf back,"%g,%g,%g,%g",Mr,Mg,Mb,Ma
	endif

#if (IgorVersion()<7)
	Execute "AppendToGizmo isoSurface="+GetWavesDataFolder(isoWave,2)+",name="+isoName
	Execute "ModifyGizmo ModifyObject="+isoName+" property={ surfaceColorType,1}"
	Execute "ModifyGizmo ModifyObject="+isoName+" property={ lineColorType,0}"
	Execute "ModifyGizmo ModifyObject="+isoName+" property={ lineWidthType,0}"
	Execute "ModifyGizmo ModifyObject="+isoName+" property={ fillMode,2}"
	Execute "ModifyGizmo ModifyObject="+isoName+" property={ lineWidth,1}"
	Execute "ModifyGizmo ModifyObject="+isoName+" property={ isoValue,"+num2str(isoValue)+"}"
	Execute "ModifyGizmo ModifyObject="+isoName+" property={ frontColor,"+back+"}"
	Execute "ModifyGizmo ModifyObject="+isoName+" property={ backColor,"+front+"}"
	Execute "ModifyGizmo modifyObject="+isoName+" property={calcNormals,1}"
	//	Execute "ModifyGizmo setDisplayList=-1, object=isoSurface0"
#else
	AppendToGizmo isoSurface=$GetWavesDataFolder(isoWave,2),name=$isoName
	ModifyGizmo ModifyObject=$isoName objectType=isosurface property={ surfaceColorType,1}
	ModifyGizmo ModifyObject=$isoName objectType=isosurface property={ lineColorType,0}
	ModifyGizmo ModifyObject=$isoName objectType=isosurface property={ lineWidthType,0}
	ModifyGizmo ModifyObject=$isoName objectType=isosurface property={ fillMode,2}
	ModifyGizmo ModifyObject=$isoName objectType=isosurface property={ lineWidth,1}
	ModifyGizmo ModifyObject=$isoName objectType=isosurface property={ isoValue,isoValue}
	ModifyGizmo ModifyObject=$isoName objectType=isosurface property={ frontColor,Mr,Mg,Mb,Ma}
	ModifyGizmo ModifyObject=$isoName objectType=isosurface property={ backColor,r,g,b,a}
	ModifyGizmo ModifyObject=$isoName objectType=isosurface property={calcNormals,1}
	//	ModifyGizmo setDisplayList=-1, object=isoSurface0
#endif
	return isoName
End

//  ================================ End of Make the Gizmo ===============================  //
//  ======================================================================================  //



//  ======================================================================================  //
//  ============================= Start of Gizmo Iso Scaling  ============================  //

Static Function/WAVE IsoColorsWave(N,colorName,revColors)			// returns a wave with rgba for each of N levels
	Variable N							// number of iso values
	String colorName
	Variable revColors				// reverse the color table
	if (!(N>0 &&N<50)  ||  WhichListItem(colorName,CTabList())<0 || numtype(revColors))
		return $""
	endif

	ColorTab2Wave $colorName
	Wave M_colors=M_colors
	if (!WaveExists(M_colors))
		return $""
	endif
	Duplicate/FREE M_colors, colors
	Redimension/D colors
	colors = M_colors/65535
	KillWaves/Z M_colors
	if (revColors)
		SetScale/I x N,1,"", colors
	else
		SetScale/I x 1,N,"", colors
	endif

	Make/N=(N,4)/D/FREE isoColors=NaN
	Variable a, i, base
	base = N<2 ? 1 : exp(-2.6/(N-1))		// -2.6 gives a min value of 0.07,  -2.6 = ln(0.07)
	for (i=1;i<=N;i+=1)
//		a = 1 - 0.2375*(i-1)
//		a = 0.6^(i-1)
		a = base^(i-1)
		isoColors[i-1][0,2] = colors(i)[q]
		isoColors[i-1][3] = a
	endfor
	return isoColors
End
//
Static Function/WAVE EqualVolumeIsoLevels(volWave,N,iMin,iMax,[include0])
	Wave volWave						// wave with all the values
	Variable N							// number of iso values
	Variable iMin						// intensity minimum in 3D data, one step below the minimum iso-value
	Variable iMax						// intensity maximum value, if N==1, then this is not used
	Variable include0					// when true, include zeros

	if (!(N>0 &&N<50) || numtype(iMin+iMax))
		return $""
	endif

	if (iMin>iMax)						// ensure that iMin <= iMax
		Variable val = iMin
		iMin = iMax
		iMax = val
	endif

	Make/N=(N)/D/FREE isoValues=NaN
	if (N==1)
		isoValues = {iMin}
		return isoValues
	elseif (N==2)
		isoValues = {iMax, iMin}
		return isoValues
	endif
	isoValues[0] = iMax
	isoValues[N-1] = iMin

	Variable fLo, fHi, df
	fLo = FindFractionAboveLevel(volWave,iMin,include0=include0)	// return level s.t. fraction of points in volWave are above level
	fHi = FindFractionAboveLevel(volWave,iMax,include0=include0)
	df = (fHi-fLo) / (N-1)

	Variable i
	for (i=1;i<(N-1);i+=1)
		isoValues[i] = FindIsoLevelByFraction(volWave,(fHi-i*df),include0=include0)
	endfor
	return isoValues
End
//
Static Function FindFractionAboveLevel(volWave,level,[include0])	// return level s.t. fraction of points in volWave are above level
	Wave volWave
	Variable level						// given level
	Variable include0					// when true, include zeros
	if (!WaveExists(volWave) || numtype(level))
		return NaN
	endif

	WaveStats/Q/M=1 volWave
	Variable Lmax=V_max, Lmin=V_min, N=V_npnts, Ntot=numpnts(volWave)
	if (N<1)								// invalid N
		return NaN
	elseif (level<=Lmin)			// fraction is zero
		return 0
	elseif (level>=Lmax)			// fraction is 1
		return 1
	endif
 	include0 = ParamIsDefault(include0) || numtype(include0) ?  (Lmin*Lmax)<0 : include0

	Duplicate/FREE volWave,tt
	Redimension/N=(Ntot) tt
	if (!include0)
		MatrixOP/FREE tt = Replace(tt,0,NaN)	// excludes zeros
		WaveStats/Q/M=1 tt
		N = V_npnts
		Lmax = V_max
		Lmin = V_min
		if (N<1)							// invalid number of points
			return NaN
		endif
	endif

	MatrixOP/FREE/O sumAbove = sum(greater(tt,level))
	return sumAbove[0]/N			// the fraction of points above level
End
//
Static Function FindIsoLevelByFraction(volWave,fraction,[include0])	// return level s.t. fraction of points in volWave are above level
	Wave volWave
	Variable fraction					// fraction of points above returned level
	Variable include0					// when true, include zeros
	if (!WaveExists(volWave))
		return NaN
	endif

	WaveStats/Q/M=1 volWave
	Variable Lmax=V_max, Lmin=V_min, N=V_npnts, Ntot=numpnts(volWave)
	if (numtype(fraction))			// invalid fraction
		return NaN
	elseif (fraction<=0)			// fraction is zero, use min value
		return Lmin
	elseif (fraction>=1)			// fraction too big, use max value
		return Lmax
	endif
 	include0 = ParamIsDefault(include0) || numtype(include0) ?  (Lmin*Lmax)<0 : include0

	Duplicate/FREE volWave,tt
	Redimension/N=(Ntot) tt
	if (!include0)
		MatrixOP/FREE tt = Replace(tt,0,NaN)	// excludes zeros
		WaveStats/Q/M=1 tt
		N = V_npnts
		Lmax = V_max
		Lmin = V_min
	endif

	Make/N=1/FREE sumAbove=-1
	Variable sumAboveLast, thresh, targetAbove=fraction*N
	do
		sumAboveLast = sumAbove[0]	// this loop avoids the need to sum the wave
		thresh = (Lmax+Lmin)/2
		MatrixOP/FREE/O sumAbove = sum(greater(tt,thresh))
		if (sumAbove[0]>targetAbove)
			Lmin = thresh
		else
			Lmax = thresh
		endif
	while(abs(sumAboveLast-sumAbove[0]))	// stop when no change
	return thresh
End

//  ============================== End of Gizmo Iso Scaling  =============================  //
//  ======================================================================================  //
