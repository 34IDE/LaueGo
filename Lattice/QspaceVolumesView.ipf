#pragma rtGlobals=1		// Use modern global access method.
#pragma ModuleName=QspaceVolumesView
#pragma version = 1.19
#include "ImageDisplayScaling", version>= 2.06
#include "ColorNames"
#include "GizmoUtility" version>= 2.15

Menu "Qspace"
	MenuItemIfWaveClassExists("Gizmo of Qspace ...","Qspace3D*","DIMS:3"), MakeGizmoQspace3D($"")
	MenuItemIfWaveClassExists("ReScale by Qn ...","Qspace3D","DIMS:3"), RescaleQspace3DbyQn($"",NaN)
	MenuItemIfWaveClassExists("Make a Radial Line for Gizmo ...","GizmoCorners","DIMS:2,MAXCOLS:3,MINCOLS:3"), MakeRadialLine($"")
	SubMenu "Utility"
		"  Make a Test 3D Qspace",QspaceVolumesView#MakeTestQspaceVolume()
	End
End


// This routine is for PLOTTING 3D Q-space arrays.   (It is NOT for creating the 3D arrays from measurements !!)


//  ======================================================================================  //
//  ================================ Start of Re-Scaling =================================  //

Function/WAVE RescaleQspace3DbyQn(Qspace3D,n)
	// do a Q^n scaling on this wave.  The Q comes from the x-y scaling
	Wave Qspace3D
	Variable n				// the power, usually 4

	if (!WaveExists(Qspace3D) || !(n>0))
		Variable ask = 1
		String wList=WaveListClass("Qspace3D","*","DIMS:3"), wName=""
		wList = ExcludeWavesInClass(wList,"Qpower")		// cannot re-scale a rescaled wave
		if (WaveExists(Qspace3D))
			wName = NameOfWave(Qspace3D)
			ask = n>0 ? 0 : 1
		else
			wName = StringFromList(0,wList)
			ask = n>0 ?  ItemsInList(wList)>1 : 1
		endif
		Wave Qspace3D = $wName
		if (ask && ItemsInList(wList)<1)
			return $""
		endif
		if (ask)
			n = n>0 ? n : 4
			Prompt wName,"3D Qspace",popup,wList
			Prompt n,"Power of Q"
			DoPrompt "Pick 3D Qspace",wName,n
			if (V_flag)
				return $""
			endif
		endif
		printf "RescaleQspace3DbyQn(%s,%g)\r",wName,n
		Wave Qspace3D = $wName
	endif
	if (!WaveExists(Qspace3D) || !(n>0))
		DoAlert 0, "Cannot find 3D Qspace Wave, or invalid power"
		return $""
	endif

	Variable qx,qy,qz
	sscanf StringByKey("qpeak", note(Qspace3D),"="), "%g, %g, %g", qx,qy,qz
	if (V_flag!=3)
		return $""
	endif

	Variable dx=DimDelta(Qspace3D,0), dy=DimDelta(Qspace3D,1), dz=DimDelta(Qspace3D,2)
	Variable x0=DimOffset(Qspace3D,0), y0=DimOffset(Qspace3D,1), z0=DimOffset(Qspace3D,2)
	Duplicate/FREE Qspace3D, Qn
	Qn = Qspace3D[p][q][r] * ( (p*dx+x0-qx)^2 + (q*dy+y0-qy)^2 + (r*dz+z0-qz)^2) ^2
	wName = NameOfWave(Qspace3D)+"Qn"
	Duplicate/O Qn, $wName
	WaveClear Qn
	Wave Qn = $wName

	String wnote=note(Qspace3D)
	wnote = AddClassToWaveNote(wnote,"Qpower")
	wnote = ReplaceNumberByKey("Qpower",wnote,n,"=")
	Note/K Qn, wnote
	printf "rescaled '%s' by q^%g, using a Bragg peak at {%g, %g, %g}(1/nm), result in '%s'\r",NameOfWave(Qspace3D),n,qx,qy,qz,wName
	return Qn
End


Function/WAVE BraggPeak(Qspace3D)
	Wave Qspace3D
	WaveStats/M=1/Q Qspace3D
	Make/N=3/D/FREE qvec
	qvec = {V_maxRowLoc, V_maxColLoc, V_maxLayerLoc}
	if (strlen(GetRTStackInfo(2))==0)
		printf "At maximum, Q = %s\r",vec2str(qvec)
	endif
	return qvec
End


#ifdef NOT_YET
Function combineQspaceVolumes(Qspace1,Qspace2)
	Wave Qspace1,Qspace2

	Variable Nx1=DimSize(Qspace1,0), Ny1=DimSize(Qspace1,1), Nz1=DimSize(Qspace1,2)
	Variable Nx2=DimSize(Qspace2,0), Ny2=DimSize(Qspace2,1), Nz2=DimSize(Qspace2,2)
	Variable dx1=DimDelta(Qspace1,0), dy1=DimDelta(Qspace1,1), dz1=DimDelta(Qspace1,2)
	Variable dx2=DimDelta(Qspace2,0), dy2=DimDelta(Qspace2,1), dz2=DimDelta(Qspace2,2)
	Variable dQ1V = dx1*dy1*dz1, dQ2V = dx2*dy2*dz2

	STRUCT boundingVolume v1
	STRUCT boundingVolume v2
	STRUCT boundingVolume v
	v1.xlo=DimOffset(Qspace1,0) ;		v1.xhi=DimOffset(Qspace1,0) + (Nx1-1)*dx1
	v1.ylo=DimOffset(Qspace1,1) ;		v1.yhi=DimOffset(Qspace1,1) + (Ny1-1)*dy1
	v1.zlo=DimOffset(Qspace1,2) ;		v1.zhi=DimOffset(Qspace1,2) + (Nz1-1)*dz1
	updateBoundingVolumeStruct(v1)
	v2.xlo=DimOffset(Qspace2,0) ;		v2.xhi=DimOffset(Qspace2,0) + (Nx2-1)*dx2
	v2.ylo=DimOffset(Qspace2,1) ;		v2.yhi=DimOffset(Qspace2,1) + (Ny2-1)*dy2
	v2.zlo=DimOffset(Qspace2,2) ;		v2.zhi=DimOffset(Qspace2,2) + (Nz2-1)*dz2
	updateBoundingVolumeStruct(v2)
	combineBoundingVolumeStructs(v1,v2,v)

	print "Measured Q-Volume1 is: ",boundingVolumeStruct2str(v1)
	print "Measured Q-Volume2 is: ",boundingVolumeStruct2str(v2)
	print "Combined Q-Volume is: ",boundingVolumeStruct2str(v)

	String wnote=note(Qspace1), wnote2=note(Qspace2)
	String scanNum=StringByKey("scanNum",wnote,"=")+","+StringByKey("scanNum",wnote2,"=")
	wnote = ReplaceStringByKey("scanNum",wnote,scanNum,"=")

	Variable qx1,qy1,qz1, qx2,qy2,qz2
	sscanf StringByKey("qpeak", wnote,"="), "%g, %g, %g", qx1,qy1,qz1
	if (V_flag!=3)
		qx1=NaN; qy1=NaN; qz1=NaN
	endif
	sscanf StringByKey("qpeak",wnote2,"="), "%g, %g, %g", qx2,qy2,qz2
	if (V_flag!=3)
		qx2=NaN; qy2=NaN; qz2=NaN
	endif
	Variable qx,qy,qz
	if (numtype(qx1+qx2)==0)			// both are good
		Variable qmax1, qmax2
		qmax1 = Qspace1(qx1)(qy1)(qz1)
		qmax2 = Qspace2(qx2)(qy2)(qz2)
		if (qmax1>=qmax2)
			qx = qx1
			qy = qy1
			qz = qz1
		else
			qx = qx2
			qy = qy2
			qz = qz2
		endif
	elseif(numtype(qx1)==0)				// only qx1 is good
		qx=qx1
		qy=qy1
		qz=qz1
	elseif(numtype(qx2)==0)				// only qx2 is good
		qx=qx2
		qy=qy2
		qz=qz2
	else										// both are bad
		qx=NaN
		qy=NaN
		qz=NaN
	endif
	String str
	sprintf str,"%g, %g, %g",qx,qy,qz
	wnote = ReplaceStringByKey("qpeak",wnote,str,"=")
	Variable dx=max(dx1,dx2), dy=max(dy1,dy2), dz=max(dz1,dz2)
	Variable dQV = dx*dy*dz
	Variable Nx=ceil(v.xW / dx)+1, Ny=ceil(v.yW / dy)+1, Nz=ceil(v.zW / dz)+1
	Make/N=(Nx,Ny,Nz)/D/O QspaceCombined=0
	SetScale/P x v.xlo,dx,"nm\S-1\M", QspaceCombined
	SetScale/P y v.ylo,dy,"nm\S-1\M", QspaceCombined
	SetScale/P z v.zlo,dz,"nm\S-1\M", QspaceCombined
	printf "dQ1vol = %g,   dQ2vol = %g,   dQvol = %g\r",dQ1V,dQ2V,dQV

	Make/N=(Nx1)/D/FREE ix
	Make/N=(Ny1)/D/FREE iy
	Make/N=(Nz1)/D/FREE iz
	ix = DimOffset(Qspace1,0) + p*dx1							// at this point ix is qx value of points in Qspace1
	iy = DimOffset(Qspace1,1) + p*dy1
	iz = DimOffset(Qspace1,2) + p*dz1
	ix = round((ix[p]-DimOffset(QspaceCombined,0))/dx)		// now ix is index into QspcaceCombined
	iy = round((iy[p]-DimOffset(QspaceCombined,1))/dy)
	iz = round((iz[p]-DimOffset(QspaceCombined,2))/dz)
	Variable i,j,k, Qval
	for (k=0;k<Nz1;k+=1)
		for (j=0;j<Ny1;j+=1)
			for (i=0;i<Nx1;i+=1)
				Qval = Qspace1[i][j][k] * (dQ1V/dQV)
				Qval = numtype(Qval) ? 0 : Qval
				QspaceCombined[ix[i]][iy[j]][iz[k]] += Qval
			endfor
		endfor
	endfor
print "max = ",WaveMax(QspaceCombined),"   max of Qspace1 = ",WaveMax(Qspace1),"   max of Qspace2 = ",WaveMax(Qspace2)
print " "

	Make/O/N=(Nx2)/D/FREE ix
	Make/O/N=(Ny2)/D/FREE iy
	Make/O/N=(Nz2)/D/FREE iz
	ix = DimOffset(Qspace2,0) + p*dx2		// at this point ix is qx value of points in Qspace2
	iy = DimOffset(Qspace2,1) + p*dy2
	iz = DimOffset(Qspace2,2) + p*dz2
	ix = round((ix[p]-DimOffset(QspaceCombined,0))/dx)		// now ix is index into QspcaceCombined
	iy = round((iy[p]-DimOffset(QspaceCombined,1))/dy)
	iz = round((iz[p]-DimOffset(QspaceCombined,2))/dz)
	for (k=0;k<Nz2;k+=1)
		for (j=0;j<Ny2;j+=1)
			for (i=0;i<Nx2;i+=1)
				Qval = Qspace2[i][j][k] *(dQ2V/dQV)
				Qval = numtype(Qval) ? 0 : Qval
				QspaceCombined[ix[i]][iy[j]][iz[k]] += Qval
			endfor
		endfor
	endfor
WaveStats/M=1/Q QspaceCombined
printf "at max:   QspaceCombined(%g)(%g)(%g) = %g\r",V_maxRowLoc,V_maxColLoc,V_maxLayerLoc,QspaceCombined(V_maxRowLoc)(V_maxColLoc)(V_maxLayerLoc)
print " "

	MatrixOP/FREE eq0 = equal(QspaceCombined,0)		// set to 1 if point is zero
	Variable zeros=sum(eq0)						// number of zero points in QspaceCombined
	wnote = ReplaceNumberByKey("zeros", wnote,zeros,"=")

	Variable Nqtot = Nx*Ny*Nz
	WaveStats/Q/M=1 QspaceCombined
	printf "V_npnts = %g, V_numNaNs = %g, V_numINFs = %g\r",V_npnts, V_numNaNs, V_numINFs
	printf "V_avg = %g,  V_Sum = %g,  V_sdev = %g\r",V_avg, V_Sum, V_sdev
	printf "there are %g zeros, and %g non-zero points  (total = %g)\r",zeros,Nqtot-zeros,Nqtot
	printf "min[%g, %g, %g] = %g\r",V_minRowLoc,V_minColLoc,V_minLayerLoc,V_min
	printf "max[%g, %g, %g] = %g\r",V_maxRowLoc,V_maxColLoc,V_maxLayerLoc,V_max

	Note/K QspaceCombined, wnote
	MakeSampledPoints(QspaceCombined,name="SampledPointsCombinedXYZ")
	MakeGizmocubeCorners(QspaceCombined)
	QspaceCombined = QspaceCombined==0 ? NaN : QspaceCombined
End
#endif

//  ================================= End of Re-Scaling ==================================  //
//  ======================================================================================  //



//  ======================================================================================  //
//  =============================== Start of Make the Gizmo ==============================  //

Function MakeGizmoQspace3D(Qspace3D,[isoMax,isoMin,Niso,ColorTable,revColors,isoValues,isoColors,volumeScaling])
	Wave Qspace3D
	Variable isoMax						// used to set biggest iso value
	Variable isoMin						// used to set smallest iso value
	Variable Niso						// number of requested iso levels
	String ColorTable					// color table to use for setting iso levels
	Variable revColors					// optionaly reverse the color table
	Wave isoValues						// user supplied values for iso surfaces
	Wave/T isoColors					// user supplied rgba for iso surfaces
	Variable volumeScaling				// 1=volume scaling, 0= linear scaling, default is volume
	isoMax = ParamIsDefault(isoMax) ? NaN : isoMax
	isoMin = ParamIsDefault(isoMin) ? NaN : isoMin
	Niso = ParamIsDefault(Niso) ? NaN : Niso
	if (ParamIsDefault(ColorTable))
		ColorTable = ""
	endif
	revColors = ParamIsDefault(revColors) ? NaN : revColors
	volumeScaling = ParamIsDefault(volumeScaling) ? 1 : volumeScaling
	volumeScaling = numtype(volumeScaling) ? 1 : !(!volumeScaling)
	if(exists("NewGizmo")!=4)			// Do nothing if the Gizmo XOP is not available.
		DoAlert 0, "Gizmo XOP must be installed"
		return 1
	endif

	Variable printIt=0
	if (!WaveExists(Qspace3D))
		String wList=WaveListClass("Qspace3D*","*","DIMS:3"), wName=""
		if (ItemsInList(wList)<2)
			wName = StringFromList(0,wList)
		else
			wName = StrVarOrDefault("root:Packages:QspaceVolumes:Gizmo:Qspace3Dname","")
			Prompt wName,"3D Qspace",popup,wList
			DoPrompt "Pick 3D Qspace",wName
			if (V_flag)
				return 1
			endif
		endif
		Wave Qspace3D = $wName
		printIt = 1
	endif
	if (!WaveExists(Qspace3D))
		DoAlert 0, "Cannot find 3D Qspace Wave"
		return 1
	endif
	String gizName=StringFromlist(0,FindGizmosWithWave(Qspace3D))	// Qspace3D already on a Gizmo, just bring it to front
	if (strlen(gizName)>0)
		DoWindow/F $gizName
		return 0
	endif

	String wnote=note(Qspace3D), str
	String userName = StringByKey("userName",wnote,"=")
	String sampleName = StringByKey("sampleName",wnote,"=")
	Variable Qpower=NumberByKey("Qpower",wnote,"=")
	String scanRange=StringByKey("scanRange",wnote,"=")
	String command=StringByKey("command",wnote,"=")
	String scanNumStr=StringByKey("scanNum",wnote,"=")
	scanNumStr = CleanupName(scanNumStr,0)
	if (strsearch(scanNumStr,"X",0)==0)		// remove leading "X" from CleanUpName
		scanNumStr = scanNumStr[1,Inf]
	endif
	Variable scanNum = str2num(scanNumStr)
	scanNum = scanNum>=0 ? scanNum : NaN
	scanNumStr = SelectString(scanNum>=0,"",num2istr(scanNum))
	scanRange = SelectString(strlen(scanRange),scanNumStr,scanRange)
	Variable X1=NumberByKey("X1",wnote,"="), Y1=NumberByKey("Y1",wnote,"="), Z1=NumberByKey("Z1",wnote,"=")
	String beamLine = StringByKey("beamLine",wnote,"=")
	String file_time = StringByKey("file_time",wnote,"=")
	String file_name = StringByKey("file_name",wnote,"=")
	Make/N=5/T/FREE wTitle=""
	wTitle[0] = StringByKey("title",wnote,"=")
	wTitle[0] += SelectString(strlen(wTitle[0]),"",", ")+sampleName
//	if (strlen(wTitle[0])<1 && strlen(command))
//		wTitle[0] = command
//	endif
	if (Qpower!=0 && numtype(Qpower)==0)
		wTitle[0] += ", Q^"+num2istr(Qpower)+" scaled"
	endif
	wTitle[1] = SelectString(strlen(scanRange),"","#"+scanRange)
	if (strlen(userName))
		wTitle[1] += SelectString(strlen(wTitle[1]),"",", ")+userName
	endif
	if (strlen(beamLine))
		wTitle[1] += SelectString(strlen(wTitle[1]),"",", ")+beamLine
	endif
	if (strlen(file_time))
		wTitle[1] += SelectString(strlen(wTitle[1]),"",", ")+file_time
	endif
	String sep = SelectString(strsearch(file_name,"\\",0)>0, ":", "\\")
	wTitle[2] = ParseFilePath(0,file_name,sep,1,0)
	if (numtype(X1+Y1+Z1)==0)
		sprintf str, "XYZ = {%g, %g, %g}",X1,Y1,Z1
		wTitle[3] = str
	endif
	if (strlen(command) && ItemsInRange(scanRange)<2)
		wTitle[4] = command
	endif
//	gizName = "Gizmo"+CleanupName(NameOfWave(Qspace3D),0)+"_"+scanNumStr
	gizName = "Gizmo"+CleanupName(NameOfWave(Qspace3D),0)

	String recip_lattice0=StringByKey("recip_lattice0",wnote,"=")
#if Exists("str2recip")
	Wave RL = str2recip(recip_lattice0)
#else
	Wave RL=$""
#endif
	Variable QxLo=NaN,QxHi=NaN, QyLo=NaN,QyHi=NaN, QzLo=NaN,QzHi=NaN, showRecip
	if (WaveExists(RL))
		Make/N=3/D/FREE astar,bstar,cstar
		astar = RL[p][0]
		bstar = RL[p][1]
		cstar = RL[p][2]
		showRecip = numtype(sum(RL))==0
	else
		Wave astar=$"", bstar=$"", cstar=$""
		showRecip = 0
	endif

	String convexList = WaveListClass("ConvexHullTriangles","*","DIMS:2,MINROWS:1,MINCOLS:3,MAXCOLS:3"), sName=""
	convexList = WavesWithMatchingKeyVals(convexList,"sourceWave="+NameOfWave(Qspace3D))
	convexList = "- none -;"+convexList
	sName = StrVarOrDefault("root:Packages:QspaceVolumes:Gizmo:ConvexHullName","")
	Prompt sName,"Show the 3D Sampled Volume",popup,convexList

//	Make/FREE isoValues = {0.2, 0.1, 0.01, 0.002, 0.001}
//	isoValues *= WaveMax(Qspace3D)
//	Make/FREE/T isoColors = {"0.5,0,0,1","0.75,0.25,0,.8","0.75,0.25,0,0.5","0.75,0.25,0,0.3","0.75,0.25,0,0.05"}
	String isoStr = ""						// make the isoStr
	Variable i, isoMaxInternal=NaN, isoMinInternal=NaN
	if (WaveExists(isoValues) && WaveExists(isoColors))	// use the user supplied iso values and colors
		if (!EqualWaves(isoValues,isoColors,512))
			DoAlert 0,"The user supplied isoValues[ ] and isoColors[ ] waves have different sizes."
			return 1
		endif
		for (i=0;i<numpnts(isoValues);i+=1)
			sprintf str,"%g:%s;", isoValues[i],isoColors[i]
			isoStr += str
		endfor
		DoPrompt "Pick 3D Qspace Sampled Points",sName
		if (V_flag)
			return 1
		endif
	else
		if (numtype(isoMax))
			// isoMaxInternal = FindIsoLevelByFraction(Qspace3D,0.0005)
			isoMaxInternal = FindIsoLevelByFraction(Qspace3D,numpnts(Qspace3D)^-0.6)
			isoMaxInternal = roundSignificant(isoMaxInternal,2)
			isoMax = numtype(isoMax) ? NumVarOrDefault("root:Packages:QspaceVolumes:Gizmo:isoMax",NaN) : isoMax
			isoMax = numtype(isoMax) ? isoMaxInternal : isoMax
		endif
		if (numtype(isoMin))
			isoMinInternal = FindIsoLevelByFraction(Qspace3D,0.25)
			isoMinInternal = roundSignificant(isoMinInternal,2)
			isoMin = numtype(isoMin) ? NumVarOrDefault("root:Packages:QspaceVolumes:Gizmo:isoMin",NaN) : isoMin
			isoMin = numtype(isoMin) ? isoMinInternal : isoMin
		endif
		if (WhichListItem(ColorTable,CTabList())<0)
			ColorTable = StrVarOrDefault("root:Packages:QspaceVolumes:Gizmo:ColorTable","Rainbow256")
		endif
		Niso = (Niso>0 &&Niso<50) ? Niso : NumVarOrDefault("root:Packages:QspaceVolumes:Gizmo:Niso",5)
		revColors = numtype(revColors) ? NumVarOrDefault("root:Packages:QspaceVolumes:Gizmo:revColors",0) : !(!revColors)

		volumeScaling = volumeScaling ? 1 : 2
		Prompt volumeScaling,"linear or volume scaling",popup,"Volume Scaling;Linear Scaling"
		Prompt ColorTable,"Color Table",popup,CTabList()
		Prompt Niso,"Number of Iso Levels",popup,"1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;"
		Prompt isoMax,"Max value for an Iso Level"
		Prompt isoMin,"Min value for an Iso Level"
		revColors += 1
		Prompt revColors,"Reverse Colors",popup,"Direct Colors;ReverseColors"
		DoPrompt "Iso Levels & Colors",ColorTable,revColors,Niso,volumeScaling,isoMax,isoMin,sName
		if (V_flag)
			return 1
		endif
		revColors -= 1
		volumeScaling = (volumeScaling==1)
		isoMax = numtype(isoMax) ? isoMaxInternal : isoMax
		isoMin = numtype(isoMin) ? isoMinInternal : isoMin
		if (!(Niso>0 &&Niso<50)  ||  WhichListItem(ColorTable,CTabList())<0 || numtype(isoMax) || numtype(isoMin))
			return 1
		endif
		isoStr = AutoMakeIsoColors(ColorTable,Niso,isoMax,isoMin,revColors,volumeScaling)
		printIt = 1
	endif

	if (printIt)
		printf "MakeGizmoQspace3D(%s",NameOfWave(Qspace3D)
		if (WaveExists(isoValues) && WaveExists(isoColors))	// use the user supplied iso values and colors
			printf ", isoValues=%s,  soColors=%s",NameOfWave(isoValues),NameOfWave(isoColors)
		else
			printf ", isoMax=%g, isoMin=%g, Niso=%g, ColorTable=\"%s\", revColors=%g",isoMax,isoMin,Niso,ColorTable,revColors
		endif
		if (!ParamIsDefault(volumeScaling) || !volumeScaling)
			printf ", volumeScaling=%g",volumeScaling
		endif
		printf ")\r"
	endif
	if (strlen(isoStr)==0)
		return 1
	endif

	String cornerList = WaveListClass("GizmoCorners","*","DIMS:2,MINROWS:2,MAXROWS:2,MINCOLS:3,MAXCOLS:3"), cName=""
	cornerList += WaveListClass("GizmoCorners","*","DIMS:2,MINROWS:8,MAXROWS:8,MINCOLS:3,MAXCOLS:3")
	cornerList = WavesWithMatchingKeyVals(cornerList,"sourceWave="+NameOfWave(Qspace3D))
	if (ItemsInList(cornerList)==1)
		cName = StringFromList(0,cornerList)
	elseif (ItemsInList(cornerList)>1)
		cName = StrVarOrDefault("root:Packages:QspaceVolumes:Gizmo:cornersName","")
		Prompt cName,"3D Qspace CORNERS",popup,cornerList
		DoPrompt "Pick 3D Qspace CORNERS",cName
		if (V_flag)
			return 1
		endif
	endif

	Wave corners = $cName
	Wave SampledVolumeXYZ = $SelectString(stringmatch(sName,"- none -"),sName,"")
	Wave gizmoScatterMarker=gizmoScatterMarker

	NewDataFolder/O root:Packages
	NewDataFolder/O root:Packages:QspaceVolumes
	NewDataFolder/O root:Packages:QspaceVolumes:Gizmo
	String/G root:Packages:QspaceVolumes:Gizmo:Qspace3Dname = NameOfWave(Qspace3D)
	if (WaveExists(corners))
		String/G root:Packages:QspaceVolumes:Gizmo:cornersName = NameOfWave(corners)
	endif
	String/G root:Packages:QspaceVolumes:Gizmo:ConvexHullName = sName
	Variable/G root:Packages:QspaceVolumes:Gizmo:isoMax = (isoMax==isoMaxInternal) ? NaN : isoMax
	Variable/G root:Packages:QspaceVolumes:Gizmo:isoMin = (isoMin==isoMinInternal) ? NaN : isoMin
	Variable/G root:Packages:QspaceVolumes:Gizmo:Niso = Niso
	String/G root:Packages:QspaceVolumes:Gizmo:ColorTable = ColorTable
	Variable/G root:Packages:QspaceVolumes:Gizmo:revColors = revColors

	String groupTitle = ""
#if (IgorVersion()<7)
	Execute "NewGizmo/N="+gizName+"/T=\""+gizName+"\" /W=(154,44,922,812)"
	Execute "ModifyGizmo startRecMacro"
	groupTitle = AddGizmoTitle(wTitle,"groupTitleQsp")
	Execute "AppendToGizmo attribute blendFunc={770,771},name=blendingFunction"
#else
	Variable Billboarding = 1
	NewGizmo/N=$gizName/T=gizName/W=(154,44,922,812)
	ModifyGizmo startRecMacro
	AddGizmoTitleGroup("groupTitleQsp",wTitle[0],title2=wTitle[1],title3=wTitle[2],title4=wTitle[3],pos="LT")
	AppendToGizmo attribute blendFunc={770,771},name=blendingFunction
#endif
	String object,displayObjectList=""

	if (WaveExists(gizmoScatterMarker))
		String CrossPathGroup=AddGizmoMarkerGroup("")		// Creates the marker for use by gizmoScatterMarker
#if (IgorVersion()<7)
		Execute "AppendToGizmo Scatter="+GetWavesDataFolder(gizmoScatterMarker,2)+",name=scatterMarker0"
		Execute "ModifyGizmo ModifyObject=scatterMarker0 property={ scatterColorType,0}"
		Execute "ModifyGizmo ModifyObject=scatterMarker0 property={ markerType,0}"
		Execute "ModifyGizmo ModifyObject=scatterMarker0 property={ sizeType,0}"
		Execute "ModifyGizmo ModifyObject=scatterMarker0 property={ rotationType,0}"
		Execute "ModifyGizmo ModifyObject=scatterMarker0 property={ Shape,7}"
		Execute "ModifyGizmo ModifyObject=scatterMarker0 property={ size,0.5}"
		Execute "ModifyGizmo ModifyObject=scatterMarker0 property={ color,0,0,0,0.5}"
		Execute "ModifyGizmo ModifyObject=scatterMarker0 property={ objectName,"+CrossPathGroup+"}"
#else
		AppendToGizmo Scatter=$GetWavesDataFolder(gizmoScatterMarker,2),name=scatterMarker0
		ModifyGizmo ModifyObject=scatterMarker0 objectType=scatter property={ scatterColorType,0}
		ModifyGizmo ModifyObject=scatterMarker0 objectType=scatter property={ markerType,0}
		ModifyGizmo ModifyObject=scatterMarker0 objectType=scatter property={ sizeType,0}
		ModifyGizmo ModifyObject=scatterMarker0 objectType=scatter property={ rotationType,0}
		ModifyGizmo ModifyObject=scatterMarker0 objectType=scatter property={ Shape,7}
		ModifyGizmo ModifyObject=scatterMarker0 objectType=scatter property={ size,0.5}
		ModifyGizmo ModifyObject=scatterMarker0 objectType=scatter property={ color,0,0,0,0.5}
		ModifyGizmo ModifyObject=scatterMarker0 objectType=scatter property={ objectName,$CrossPathGroup}
#endif
		displayObjectList += "scatterMarker0;"
	endif

	if (WaveExists(corners))
		displayObjectList += "scatterCubeCorners;"
		Variable icor=DimSize(corners,0)-1			// this will be 2 or 7 (7 is old way)
		ImageStats/M=1/G={0,icor, 0,0} corners;		QxLo = V_min	;	QxHi = V_max	// need these later for setting a*, b*, c*
		ImageStats/M=1/G={0,icor, 1,1} corners;		QyLo = V_min	;	QyHi = V_max
		ImageStats/M=1/G={0,icor, 2,2} corners;		QzLo = V_min	;	QzHi = V_max
#if (IgorVersion()<7)
		Execute "AppendToGizmo Scatter="+GetWavesDataFolder(corners,2)+",name=scatterCubeCorners"
		Execute "ModifyGizmo ModifyObject=scatterCubeCorners property={ scatterColorType,0}"
		Execute "ModifyGizmo ModifyObject=scatterCubeCorners property={ markerType,0}"
		Execute "ModifyGizmo ModifyObject=scatterCubeCorners property={ sizeType,0}"
		Execute "ModifyGizmo ModifyObject=scatterCubeCorners property={ rotationType,0}"
		Execute "ModifyGizmo ModifyObject=scatterCubeCorners property={ Shape,1}"
		Execute "ModifyGizmo ModifyObject=scatterCubeCorners property={ size,1}"
		Execute "ModifyGizmo ModifyObject=scatterCubeCorners property={ color,0,0,0,0}"
#else
		AppendToGizmo Scatter=$GetWavesDataFolder(corners,2),name=scatterCubeCorners
		ModifyGizmo ModifyObject=scatterCubeCorners objectType=scatter property={ scatterColorType,0}
		ModifyGizmo ModifyObject=scatterCubeCorners objectType=scatter property={ markerType,0}
		ModifyGizmo ModifyObject=scatterCubeCorners objectType=scatter property={ sizeType,0}
		ModifyGizmo ModifyObject=scatterCubeCorners objectType=scatter property={ rotationType,0}
		ModifyGizmo ModifyObject=scatterCubeCorners objectType=scatter property={ Shape,1}
		ModifyGizmo ModifyObject=scatterCubeCorners objectType=scatter property={ size,1}
		ModifyGizmo ModifyObject=scatterCubeCorners objectType=scatter property={ color,0,0,0,0}
#endif
	else
		QxLo = DimOffset(Qspace3D,0)	;	QxHi = QxLo + (DimSize(Qspace3D,0)-1) * DimDelta(Qspace3D,0)
		QyLo = DimOffset(Qspace3D,1)	;	QyHi = QyLo + (DimSize(Qspace3D,1)-1) * DimDelta(Qspace3D,1)
		QzLo = DimOffset(Qspace3D,2)	;	QzHi = QzLo + (DimSize(Qspace3D,2)-1) * DimDelta(Qspace3D,2)
	endif

	if (showRecip)
		Make/N=3/D/FREE Qc = {(QxLo+QxHi)/2, (QyLo+QyHi)/2, (QzLo+QzHi)/2}		// center of volume
		printf "Qx=[%g, %g],  Qy=[%g, %g],  Qz=[%g, %g], Qc={%s}\r",QxLo,QxHi, QyLo,QyHi, QzLo,QzHi,vec2str(Qc,bare=1)
		MatrixOP/FREE hklCenter = Inv(RL) x Qc
		Variable dLo,dHi

		scaleRange3D(astar,Qc,QxLo,QxHi,QyLo,QyHi,QzLo,QzHi,dLo,dHi)
		printf "H range = [%g, %g],  Æ=%g,  centered at %s\r",hklCenter[0]+dLo,hklCenter[0]+dHi,dHi-dLo, vec2str(hklCenter)
		displayObjectList += AddGizmoRecipAxis(astar,Qc,dLo,dHi,hklCenter[0],"astarRecipAxis","red")+";"

		scaleRange3D(bstar,Qc,QxLo,QxHi,QyLo,QyHi,QzLo,QzHi,dLo,dHi)
		printf "K range = [%g, %g],  Æ=%g\r",hklCenter[1]+dLo,hklCenter[1]+dHi,dHi-dLo
		displayObjectList += AddGizmoRecipAxis(bstar,Qc,dLo,dHi,hklCenter[1],"bstarRecipAxis","green")+";"

		scaleRange3D(cstar,Qc,QxLo,QxHi,QyLo,QyHi,QzLo,QzHi,dLo,dHi)
		printf "L range = [%g, %g],  Æ=%g\r",hklCenter[2]+dLo,hklCenter[2]+dHi,dHi-dLo
		displayObjectList += AddGizmoRecipAxis(cstar,Qc,dLo,dHi,hklCenter[2],"cstarRecipAxis","blue")+";"
	endif

	Variable val, r,g,b,a
	printf " Iso Values\t\t\t\tColors  {r,g,b,a}\r"
	for (i=0;i<ItemsInList(isoStr);i+=1)
		str = StringFromList(i,isoStr)
		val = str2num(StringFromList(0,str,":"))			// isovalue
		str = StringFromList(1,str,":")						// iso color, "rgba"
		displayObjectList += AddIso2Gizmo(Qspace3D,"isoSurface"+num2istr(i),val,str,"")+";"
		if (printIt)
			sscanf str,"%g,%g,%g,%g",r,g,b,a
			printf "%g   \t\t'%s'     {%s}\r",val,RGBA2name(r,g,b,a,1),str
		endif
	endfor

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
	Execute "ModifyGizmo ModifyObject=axesBeamLineQ,property={0,axisLabelText,\"Qx  (1/nm)\"}"
	Execute "ModifyGizmo ModifyObject=axesBeamLineQ,property={1,axisLabelText,\"Qy  (1/nm)\"}"
	Execute "ModifyGizmo ModifyObject=axesBeamLineQ,property={2,axisLabelText,\"Qz  (1/nm)\"}"
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
	ModifyGizmo ModifyObject=axesBeamLineQ, objectType=axes property={0,axisLabelText,"Qx  (1/nm)"}
	ModifyGizmo ModifyObject=axesBeamLineQ, objectType=axes property={1,axisLabelText,"Qy  (1/nm)"}
	ModifyGizmo ModifyObject=axesBeamLineQ, objectType=axes property={2,axisLabelText,"Qz  (1/nm)"}
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

	if (WaveExists(SampledVolumeXYZ))
#if (IgorVersion()<7)
		Execute "AppendToGizmo Surface="+GetWavesDataFolder(SampledVolumeXYZ,2)+",name=sampledVolumeSurface"
		Execute "ModifyGizmo ModifyObject=sampledVolumeSurface property={ surfaceColorType,1}"
		Execute "ModifyGizmo ModifyObject=sampledVolumeSurface property={ srcMode,1}"
		Execute "ModifyGizmo ModifyObject=sampledVolumeSurface property={ frontColor,0.749996,0.811093,1,0.1}"
		Execute "ModifyGizmo ModifyObject=sampledVolumeSurface property={ backColor,0.749996,0.811093,1,0.1}"
#else
		AppendToGizmo Surface=$GetWavesDataFolder(SampledVolumeXYZ,2),name=sampledVolumeSurface
		ModifyGizmo ModifyObject=sampledVolumeSurface objectType=surface property={ surfaceColorType,1}
		ModifyGizmo ModifyObject=sampledVolumeSurface objectType=surface property={ srcMode,1}
		ModifyGizmo ModifyObject=sampledVolumeSurface objectType=surface property={ frontColor,0.749996,0.811093,1,0.1}
		ModifyGizmo ModifyObject=sampledVolumeSurface objectType=surface property={ backColor,0.749996,0.811093,1,0.1}
#endif
		displayObjectList += "sampledVolumeSurface;"
	endif

	if (strlen(groupTitle))
#if (IgorVersion()<7)
		Execute "ModifyGizmo setDisplayList=-1, object="+groupTitle
		Execute "ModifyGizmo setDisplayList=-1, opName=MainTransform, operation=mainTransform"
#endif
	endif
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
	Execute "ModifyGizmo SETQUATERNION={0.628289,0.154694,-0.151222,0.747298}"
	Execute "ModifyGizmo autoscaling=1"
	Execute "ModifyGizmo currentGroupObject=\"\""
	Execute "ModifyGizmo compile"
//	Execute "ModifyGizmo showInfo"
//	Execute "ModifyGizmo infoWindow={987,235,1488,603}"
	Execute "ModifyGizmo bringToFront"
	Execute "ModifyGizmo endRecMacro"
#else
	ModifyGizmo SETQUATERNION={0.628289,0.154694,-0.151222,0.747298}
	ModifyGizmo autoscaling=1
	ModifyGizmo currentGroupObject=""
	ModifyGizmo compile
//	ModifyGizmo showInfo
//	ModifyGizmo infoWindow={987,235,1488,603}
	ModifyGizmo bringToFront
	ModifyGizmo endRecMacro
#endif
	return 0
End
//
// //   This was moved to LatticeSym.ipf
//Static Function/WAVE str2recip(str)		// returns a FREE wave with reciprocal lattice
//	String str
//	Variable as0,as1,as2,bs0,bs1,bs2,cs0,cs1,cs2
//	sscanf str, "{{%g,%g,%g}{%g,%g,%g}{%g,%g,%g}}",as0,as1,as2,bs0,bs1,bs2,cs0,cs1,cs2
//	if (V_flag==9)
//		Make/N=(3,3)/D/FREE RL
//		RL[0][0]= {as0,as1,as2}						// the reciprocal lattice
//		RL[0][1]= {bs0,bs1,bs2}
//		RL[0][2]= {cs0,cs1,cs2}
//		return RL
//	else
//		return $""
//	endif
//End
//
Static Function/T AutoMakeIsoColors(colorName,N,iMax,iMin,revColors,volumeScaling)	// uses ~volume scaling (not log, not linear)
	String colorName
	Variable N							// number of iso values
	Variable iMax						// maximum value
	Variable iMin						// minimum in 3D data, one step below the minimum iso-value, ONLY used with linear scaling (volumeScaling=0)
	Variable revColors					// reverse the color table
	Variable volumeScaling				// 1=volume scaling, 0= linear scaling, default is volume

	if (!(N>0 &&N<50)  ||  WhichListItem(colorName,CTabList())<0 || numtype(revColors))
		return ""
	endif

	Make/N=(N)/O/FREE isoValues
	Make/N=(N)/T/O/FREE isoColors
	volumeScaling = (iMax<=0 || iMin<=0) ? 0 : volumeScaling	// volume scaling does not work right with negatives or zeros

	if (volumeScaling)					// volume scaling
		// isoValues = iMax * 1 / ( (p+1)^4 )
		Variable m = ln(iMax/iMin)/ln(N)		// power
		isoValues = iMax * 1 / ( (p+1)^m )
	else
		Variable delta = (iMax-iMin)/N
		isoValues = iMax - p*delta
	endif

	ColorTab2Wave $colorName
	Wave M_colors=M_colors
	if (!WaveExists(M_colors))
		return ""
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
	String str
	Variable r,g,b,a, i, base
	base = N<2 ? 1 : exp(-2.6/(N-1))		// -2.6 gives a min value of 0.07,  -2.6 = ln(0.07)
	for (i=1;i<=N;i+=1)
//		a = 1 - 0.2375*(i-1)
//		a = 0.6^(i-1)
		a = base^(i-1)
		sprintf str,"%g,%g,%g,%g", colors(i)[0], colors(i)[1], colors(i)[2],a
		isoColors[i-1] = str
	endfor

	String out=""
	for (i=0;i<N;i+=1)
		out += num2str(isoValues[i])+":"+isoColors[i]+";"
	endfor

	return out
End
//
Static Function FindIsoLevelByFraction(Qspace3D,fraction)	// return level s.t. fraction of points in Qspace3D are above level
	Wave Qspace3D
	Variable fraction				// fraction of points above returned level
	if (!WaveExists(Qspace3D))
		return NaN
	endif

	WaveStats/Q/M=1 Qspace3D
	Variable Lmax=V_max, Lmin=V_min, N=V_npnts, Ntot=numpnts(Qspace3D)
	if (numtype(fraction))			// invalid fractoin
		return NaN
	elseif (fraction<=0)			// fraction is zero, use min value
		return Lmin
	elseif (fraction>=1)			// fraction too big, use max value
		return Lmax
	endif
 
	Duplicate/FREE Qspace3D,tt
	Redimension/N=(Ntot) tt
	Make/N=1/FREE sumAbove=-Inf
	Variable sumAboveLast, thresh, targetAbove=fraction*N
	do
		sumAboveLast = sumAbove[0]
		thresh = (Lmax+Lmin)/2
		MatrixOP/FREE/O sumAbove = sum(greater(tt,thresh))
		if (sumAbove[0]>targetAbove)
			Lmin = thresh
		else
			Lmax = thresh
		endif
	while(abs(sumAboveLast-sumAbove[0]))
	return thresh
End
//
Static Function/T AddIso2Gizmo(isoWave,isoName,isoValue,front,back)
	Wave isoWave						// 3d wave
	String isoName			// probably astar, bstar, or cstar, or maybe h,k,l
	Variable isoValue
	String front,back

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
//
Static Function/T AddGizmoRecipAxis(rlv,Qc,dLo,dHi,hCenter,name,rgba,[alpha])
	Wave rlv					// reciprocal lattice vector
	Wave Qc					// center of cube (in 1/nm)
	Variable dLo,dHi,hCenter	// length of rlv to show (in rlv units, e.g. h, k, or l), range is [hCenter+dLo, hCenter+dHi]
	String name					// probably astarRecipAxis, bstarRecipAxis, or cstarRecipAxis
	String rgba					// red, green, blue, or "" is black, or you can give your own rgba as "1,1,0,0.5"
	Variable alpha
	alpha = ParamIsDefault(alpha) ? 1 : alpha
	alpha = numtype(alpha) ? 1 : limit(alpha,0,1)

	if (!WaveExists(rlv) || !WaveExists(Qc))		// either a*, b*, c*
		return ""
	elseif (numtype(dLo+dHi))
		return ""
	endif

	String colors="red:1,0,0;green:0,0.6,0;blue:0,0,1;cyan:0,1,1;magenta:1,0,1;yellow:1,1,0;black:0,0,0;white:1,1,1;", str
	Variable r,g,b,a
	sscanf rgba,"%g,%g,%g,%g", r,g,b,a
	if (V_flag!=4)
		rgba = StringByKey(rgba,colors)+","+num2str(alpha)
	endif
	if (strlen(rgba)<1)
		rgba = "0,0,0,1"
	endif


	if (0)
		Make/N=3/D/FREE vec=rlv, xhat={1,0,0}
		normalize(vec)
		Cross xhat, vec
		Wave W_Cross=W_Cross
		Make/N=3/D/FREE axis = W_Cross				// axis of rotation
		KillWaves/Z W_Cross
		Variable angle = asin(normalize(axis))*180/PI// rotation angle (degree)
		printf "axis = %s,   angle = %g¡\r",vec2str(axis),angle
	endif


	Make/N=(2,3)/FREE/D qs=NaN
	qs[0][] = Qc[q] + dLo*rlv[q]
	qs[1][] = Qc[q] + dHi*rlv[q]

	// print "rlv =",vec2str(rlv,bare=1)," hkl range = ",hCenter+dLo,hCenter+dHi
#if (IgorVersion()<7)
	sprintf str, "ModifyGizmo ModifyObject=%s,property={0,axisRange,%g,%g,%g,%g,%g,%g}", name,qs[0][0],qs[0][1],qs[0][2], qs[1][0],qs[1][1],qs[1][2]
	// print str
	// print " "
	Execute "AppendToGizmo Axes=CustomAxis,name="+name
	Execute str
	Execute "ModifyGizmo ModifyObject="+name+", property={0,lineWidth,2}"
	Execute "ModifyGizmo ModifyObject="+name+", property={0,axisScalingMode,1}"
	Execute "ModifyGizmo ModifyObject="+name+", property={0,axisColor,"+rgba+"}"
	Execute "ModifyGizmo ModifyObject="+name+", property={0,ticks,3}"
	Execute "ModifyGizmo ModifyObject="+name+", property={0,labelRotationAxis,5,1,1}"
	Execute "ModifyGizmo ModifyObject="+name+", property={0,labelRotationAngle,43.443}"
	Execute "ModifyGizmo ModifyObject="+name+", property={0,labelOffset,-0.15,-0.05,0.05}"
	Execute "ModifyGizmo ModifyObject="+name+", property={0,labelColor,"+rgba+"}"
	Execute "ModifyGizmo ModifyObject="+name+", property={0,fontScaleFactor,0.7}"
	Execute "ModifyGizmo ModifyObject="+name+", property={0,numTicks,7}"
	Execute "ModifyGizmo ModifyObject="+name+", property={0,axisMinValue,"+num2str(hCenter+dLo)+"}"
	Execute "ModifyGizmo ModifyObject="+name+", property={0,axisMaxValue,"+num2str(hCenter+dHi)+"}"
	Execute "ModifyGizmo modifyObject="+name+", property={Clipped,0}"
#else
	Wave rgbaVec = str2vec(rgba)
	AppendToGizmo Axes=CustomAxis,name=$name
	ModifyGizmo ModifyObject=$name, objectType=axes property={0,axisRange,qs[0][0],qs[0][1],qs[0][2], qs[1][0],qs[1][1],qs[1][2]}
	ModifyGizmo ModifyObject=$name, objectType=axes property={0,lineWidth,2}
	ModifyGizmo ModifyObject=$name, objectType=axes property={0,axisScalingMode,1}
	ModifyGizmo ModifyObject=$name, objectType=axes property={0,axisColor, rgbaVec[0],rgbaVec[1],rgbaVec[2],rgbaVec[3]}
	ModifyGizmo ModifyObject=$name, objectType=axes property={0,ticks,3}
	ModifyGizmo ModifyObject=$name, objectType=axes property={0,labelRotationAxis,5,1,1}
	ModifyGizmo ModifyObject=$name, objectType=axes property={0,labelRotationAngle,43.443}
	ModifyGizmo ModifyObject=$name, objectType=axes property={0,labelOffset,-0.15,-0.05,0.05}
	ModifyGizmo ModifyObject=$name, objectType=axes property={0,labelColor, rgbaVec[0],rgbaVec[1],rgbaVec[2],rgbaVec[3]}
	ModifyGizmo ModifyObject=$name, objectType=axes property={0,fontScaleFactor,0.7}
	ModifyGizmo ModifyObject=$name, objectType=axes property={0,numTicks,7}
	ModifyGizmo ModifyObject=$name, objectType=axes property={0,axisMinValue, hCenter+dLo}
	ModifyGizmo ModifyObject=$name, objectType=axes property={0,axisMaxValue, hCenter+dHi}
	ModifyGizmo ModifyObject=$name, objectType=axes property={-1,Clipped,0}
#endif
	return name
End
//
//
// Example of how to use AddGizmoTitle()
//
//	Make/N=1/T/FREE wTitle="Title to Display"
//	String groupTitle = AddGizmoTitle(wTitle,"")
//	if (strlen(groupTitle))
//		Execute "ModifyGizmo setDisplayList=-1, object="+groupTitle
//		Execute "ModifyGizmo setDisplayList=-1, opName=MainTransform, operation=mainTransform"
//	endif
//
//	*********************************************************
//	NOTE, this whole function is really only useful for Igor6.
//	*********************************************************
//
Static Function/T AddGizmoTitle(wTitle,groupName)
	Wave/T wTitle
	String groupName							// defaults to "groupTitle"
	if (WaveType(wTitle,1)!=2)				// wTitle must be a text wave
		return ""
	endif
	groupName = SelectString(strlen(groupName),"groupTitle",groupName)

	Variable i,m,N=max(numpnts(wTitle),10)	// at most 10 lines
	for (i=0;i<N;i+=1)
		if (strlen(wTitle[i]))
			wTitle[m] = wTitle[i]
			if (m<i)
				wTitle[i] = ""
			endif
			m += 1
		endif
	endfor
	N = m										// number of non-blank lines
	if (N<1)
		return ""								// wTitle is empty
	endif

	String name, cmd
	Variable scale=1,trans=0

	// ************************* Group Object Start *******************
#if (IgorVersion()<7)
	Execute "AppendToGizmo group,name="+groupName
	Execute "ModifyGizmo currentGroupObject=\""+groupName+"\""
#else
	AppendToGizmo group,name=$groupName
	ModifyGizmo currentGroupObject=groupName
#endif

	for (i=0;i<N;i+=1)
		name = "Title"+num2str(i)
#if (IgorVersion()<7)
		Execute "AppendToGizmo string=\""+wTitle[i]+"\",strFont=\""+GenevaEquivFont+"\",name="+name
		Execute "ModifyGizmo modifyObject="+name+" property={Clipped,0}"
#else
		AppendToGizmo string=wTitle[i],strFont=GenevaEquivFont,name=$name
// xxxxxxxxxxxxxxxxxxxxxxxxxxx
//		ModifyGizmo modifyObject=$name objectType=string property={Clipped,0}
#endif
	endfor

#if (IgorVersion()<7)
	Execute "ModifyGizmo setDisplayList=0, opName=translateTitle, operation=translate, data={-1.9,1.9,0}"
	Execute "ModifyGizmo setDisplayList=1, opName=scaleTitle, operation=scale, data={0.1,0.1,0.1}"
	Execute "ModifyGizmo setDisplayList=2, opName=rotateTitle, operation=rotate, data={180,1,0,0}"
	Execute "ModifyGizmo setDisplayList=3, object=Title0"
#else
	ModifyGizmo setDisplayList=0, opName=translateTitle, operation=translate, data={-1.9,1.9,0}
	ModifyGizmo setDisplayList=1, opName=scaleTitle, operation=scale, data={0.1,0.1,0.1}
	ModifyGizmo setDisplayList=2, opName=rotateTitle, operation=rotate, data={180,1,0,0}
	ModifyGizmo setDisplayList=3, object=Title0
#endif

	for (i=1;i<N;i+=1)
		if (i==1)
			scale = 0.8
			trans = 1
		elseif (i==2)
			scale = 0.7
			trans = 1
		else
			scale = 1
			trans = 1.5
		endif
		name = "Title"+num2str(i)

#if (IgorVersion()<7)
		sprintf cmd, "ModifyGizmo setDisplayList=-1, opName=translate%s, operation=translate, data={0,%g,0}",name,trans
		Execute cmd
		if (scale<1)
			sprintf cmd, "ModifyGizmo setDisplayList=-1, opName=scale%s, operation=scale, data={%g,%g,%g}",name,scale,scale,scale
			Execute cmd
		endif
		Execute "ModifyGizmo setDisplayList=-1, object="+name
#else
		ModifyGizmo setDisplayList=-1, opName=$("translate"+name), operation=translate, data={0,trans,0}
		if (scale<1)
			ModifyGizmo setDisplayList=-1, opName=$("scale"+name), operation=scale, data={scale,scale,scale}
		endif
		ModifyGizmo setDisplayList=-1, object=$name
#endif
	endfor
#if (IgorVersion()<7)
	Execute "ModifyGizmo currentGroupObject=\"::\""
#else
	ModifyGizmo currentGroupObject="::"
#endif
	// ************************* Group Object End *******************
	return groupName
End
//
//
// find scale factor for scale*vec that starts at center and first hits a box wall, returns ± scale factors
Static Function scaleRange3D(vec,center,Xlo,Xhi,Ylo,Yhi,Zlo,Zhi,lo,hi)
	Wave vec							// vector direction
	Wave center						// vec starts at center
	Variable Xlo,Xhi,Ylo,Yhi,Zlo,Zhi		// box limits
	Variable &lo, &hi

	// find intersection of line and plane
	Variable vn, s0
	Make/N=3/D/FREE planeLow={Xlo,Ylo,Zlo}, planeHigh={Xhi,Yhi,Zhi}, dplane, n

	lo = -Inf 
	hi = Inf

	n={1,0,0}
	vn = MatrixDot(vec,n)
	if (vn!=0)
		dplane = planeLow - center				// plane at x == Xlo
		s0 = MatrixDot(dplane,n)/vn
		if (s0<0)
			lo = max(lo,s0)
		else
			hi = min(hi,s0)
		endif
		dplane = planeHigh - center				// plane at x == Xhi
		s0 = MatrixDot(dplane,n)/vn
		if (s0<0)
			lo = max(lo,s0)
		else
			hi = min(hi,s0)
		endif
	endif

	n = {0,1,0}
	vn = MatrixDot(vec,n)
	if (vn!=0)
		dplane = planeLow - center				// plane at y == Ylo
		s0 = MatrixDot(dplane,n)/vn
		if (s0<0)
			lo = max(lo,s0)
		else
			hi = min(hi,s0)
		endif
		dplane = planeHigh - center				// plane at y == Yhi
		s0 = MatrixDot(dplane,n)/vn
		if (s0<0)
			lo = max(lo,s0)
		else
			hi = min(hi,s0)
		endif
	endif

	n = {0,0,1}
	vn = MatrixDot(vec,n)
	if (vn!=0)
		dplane = planeLow - center				// plane at z == Zlo
		s0 = MatrixDot(dplane,n)/vn
		if (s0<0)
			lo = max(lo,s0)
		else
			hi = min(hi,s0)
		endif
		dplane = planeHigh - center				// plane at z == Zhi
		s0 = MatrixDot(dplane,n)/vn
		if (s0<0)
			lo = max(lo,s0)
		else
			hi = min(hi,s0)
		endif
	endif

	return numtype(lo+hi)
End



Static Function/WAVE MakeConvexHullFrom3D(Qspace3D,[name,tolerance])
	Wave Qspace3D						// a 3D wave that we want to find the outside of,  assume that smallest value defines outside (probably 0)
	String name							// user can supply the name for convex hull wave
	Variable tolerance					// This is needed because the ConvexHull command has trouble with too many points with the same value, 1e-4 is good for Q-distn
	name = SelectString(ParamIsDefault(name),name,"")
	tolerance = ParamIsDefault(tolerance) ? 1e-4 : tolerance
	tolerance = tolerance>=0 ? tolerance : 1e-4

	if (WaveDims(Qspace3D)!=3)
		String wName, wlist=WaveListClass("Qspace3D*","*","DIMS:3")
		if (ItemsInLIst(wlist)<1)
			return $""
		elseif (ItemsInLIst(wlist)==1)
			Wave Qspace3D = $StringFromList(0,wlist)
		else
			Prompt wName,"Wave of 3D space for Convex Hull",popup,wlist
			DoPrompt "Pick 3D Wave",wName
			if (V_flag)
				return $""
			endif
			Wave Qspace3D = $wName
		endif
	endif
	if (WaveDims(Qspace3D)!=3)
		return $""
	endif

	Variable Nx=DimSize(Qspace3D,0), Ny=DimSize(Qspace3D,1), Nz=DimSize(Qspace3D,2)
	Variable x0=DimOffset(Qspace3D,0), y0=DimOffset(Qspace3D,1), z0=DimOffset(Qspace3D,2)
	Variable dx=DimDelta(Qspace3D,0), dy=DimDelta(Qspace3D,1), dz=DimDelta(Qspace3D,2)
	Make/N=(Nx*Ny*Nz,3)/FREE/D xyzTemp=NaN
	Make/N=3/D/FREE xyzi
	WaveStats/Q/M=1 Qspace3D
	Variable ix,iy,iz,N, threshold=V_min, xx,yy,zz
	for (N=0,iz=0; iz<Nz; iz+=1)
		zz =  iz*dz + z0
		for (iy=0;iy<Ny;iy+=1)
			yy =  iy*dy + y0
			for (ix=0;ix<Nx;ix+=1)
				if (Qspace3D[ix][iy][iz]>threshold)		// only save inside values, this also skips NaN
					xx =  ix*dx + x0
					xyzi = {xx,yy,zz}
					xyzTemp [N][] = xyzi[q]
					N += 1
				endif
			endfor
		endfor
	endfor
	Redimension/N=(N,3) xyzTemp

	xyzTemp += gnoise(tolerance)
	ConvexHull xyzTemp
	WaveClear xyzTemp
	Wave M_Hull=M_Hull

	if (strlen(name)<1)
		name = "ConvexHullTrianlges"+NameOfWave(Qspace3D)
	endif
	name =  CleanupName(name,0)
	Duplicate/O M_Hull, $name
	Wave hull=$name
	KillWaves/Z M_Hull

	String wnote=note(Qspace3D)
	wnote = ReplaceStringByKey("waveClass",wnote,"ConvexHullTriangles","=")
	wnote = ReplaceStringByKey("sourceWave",wnote,NameOfWave(Qspace3D),"=")
	Note/K hull, wnote
	return hull
End



Static Function MakeTestQspaceVolume()				// Make a test volume
	Variable Nx=21, Ny=21,Nz=21
	Variable Nqtot = Nx*Ny*Nz

	STRUCT boundingVolume v
	initBoundingVolumeStruct(v)
	v.xlo = -20 ;		v.xhi = +20
	v.ylo = -20 ;		v.yhi = +20
	v.zlo = 140 ;		v.zhi = 180
	updateBoundingVolumeStruct(v)
	print "Measured Q-Volume is: ",boundingVolumeStruct2str(v)
	Variable dQ = (v.xW)/(Nx-1)
	printf "dQ = %g 1/nm,   N = [%g, %g, %g],   %g points\r",dQ,Nx,Ny,Nz,Nqtot

	Make/N=(Nx,Ny,Nz)/D/O Qspace3D=0
	SetScale/P x v.xlo,dQ,"nm\S-1\M", Qspace3D
	SetScale/P y v.ylo,dQ,"nm\S-1\M", Qspace3D
	SetScale/P z v.zlo,dQ,"nm\S-1\M", Qspace3D

	Variable xc = (v.xlo + v.xhi)/2
	Variable yc = (v.ylo + v.yhi)/2
	Variable zc = (v.zlo + v.zhi)/2
	Variable wid = min(min(v.xW,v.yW),v.zW)/2
	Qspace3D = exp(-((x-xc)^2 + (y-yc)^2 + (z-zc)^2)/wid)

	String wnote = ReplaceStringByKey("waveClass","","Qspace3D","=")
	wnote = ReplaceStringByKey("title",wnote,"Test","=")
	wnote = ReplaceStringByKey("userName",wnote,"JZT","=")
	wnote = ReplaceStringByKey("sampleName",wnote,"Si","=")
	wnote = ReplaceNumberByKey("scanNum",wnote,1,"=")
	wnote = ReplaceNumberByKey("kev",wnote,8.05092114285714,"=")
	Variable as = 2*PI/0.154
	String str
	sprintf str, "{{%g,%g,%g}{%g,%g,%g}{%g,%g,%g}}",as,0,0,0,as,0,0,0,as
	wnote = ReplaceStringByKey("recip_lattice0",wnote,str,"=")
	Wave qvec = BraggPeak(Qspace3D)
	str = vec2str(qvec,bare=1)
	str = ReplaceString(" ",str,"")
	wnote = ReplaceStringByKey("qPeak", wnote,str,"=")
	Note/K Qspace3D, wnote

	WaveStats/Q/M=1 Qspace3D
	printf "V_npnts = %g, V_numNaNs = %g, V_numINFs = %g\r",V_npnts, V_numNaNs, V_numINFs
	printf "V_avg = %g,  V_Sum = %g,  V_sdev = %g\r",V_avg, V_Sum, V_sdev
	printf "min[%g, %g, %g] = %g\r",V_minRowLoc,V_minColLoc,V_minLayerLoc,V_min
	printf "max[%g, %g, %g] = %g\r",V_maxRowLoc,V_maxColLoc,V_maxLayerLoc,V_max

	Variable dx=DimDelta(Qspace3D,0), dy=DimDelta(Qspace3D,1), dz=DimDelta(Qspace3D,2)
	Variable x0=DimOffset(Qspace3D,0), y0=DimOffset(Qspace3D,1), z0=DimOffset(Qspace3D,2)
	Variable xx,yy,zz
	Make/N=(Nqtot,3)/O SampledPointsXYZ=NaN
	Variable i,j,k, m, edge=3
	for (k=edge,m=0;k<Nz-edge;k+=1)
		zz = k*dz + z0
		for (j=edge;j<Ny-edge;j+=1)
			yy = j*dy + y0
			for (i=edge;i<Nx-edge;i+=1)
				xx = i*dx + x0
				SampledPointsXYZ[m][0] = xx
				SampledPointsXYZ[m][1] = yy
				SampledPointsXYZ[m][2] = zz
				m += 1
			endfor
		endfor
	endfor
	Redimension/N=(m,-1) SampledPointsXYZ
	wnote="waveClass=sampledPoints;"
	wnote = ReplaceStringByKey("sourceWave",wnote,NameOfWave(Qspace3D),"=")
	wnote = ReplaceStringByKey("sourceWaveFullPath",wnote,GetWavesDataFolder(Qspace3D,2),"=")
	Note/K SampledPointsXYZ, wnote

	MakeConvexHullFrom3D(Qspace3D,name="ConvexHullTrianlges")
//	Duplicate/FREE SampledPointsXYZ,xyzTemp
//	Redimension/D xyzTemp
//	xyzTemp += gnoise(1e-4)
//	ConvexHull xyzTemp
//	WaveClear xyzTemp
//	Wave M_Hull=M_Hull
//	String name = CleanupName("ConvexHullTrianlges",0)
//	Duplicate M_Hull, $name
//	Wave hull=$name
//	KillWaves/Z M_Hull
//	wnote = ReplaceStringByKey("waveClass",wnote,"ConvexHullTriangles","=")
//	Note/K hull, wnote

	MakeGizmocubeCorners(Qspace3D)
	Qspace3D = Qspace3D==0 ? NaN : Qspace3D

	return 0
End

//  =============================== End of Make the Gizmo ================================  //
//  ======================================================================================  //




//  ======================================================================================  //
//  ============================ Start of Gizmo Enhancements =============================  //

Function/WAVE MakeRadialLine(corners,[point,printIt])
	Wave corners							// corners of the volume, used to set ends of line
	String point							// string with a 3-vector that is a point on the radial line
	Variable printIt
	point = SelectString(ParamIsDefault(point),point,"")
	printIt = ParamIsDefault(printIt) ? NaN : printIt
	printIt = numtype(printIt) ? (strlen(GetRTStackInfo(2)) < 1) : !(!printIt)

	Variable ask=!WaveExists(corners)
	Variable pointOK=1
	if (!ParamIsDefault(point))		// if a default, assume the center of volume
		Wave vcTemp = str2vec(point)
		pointOK = numpnts(vcTemp)==3
		ask = ask || !pointOK
	endif

	String name=""
	if (ask)
		String clist=WaveListClass("GizmoCorners","*","DIMS:2,MAXCOLS:3,MINCOLS:3")
		if (ItemsInList(clist)==1 && pointOK)
			Wave corners=$StringFromList(0,clist)	// only one choice, take it
		elseif (ItemsInList(clist)>=1)				// multiple choices, ask user
			Prompt name,"Wave with Corners of Gizmo",popup,clist
			Prompt point,"point on radial line, use \"\" for default"
			DoPrompt "Corners Wave",name, point
			if (V_flag)
				print "Could not find the corners wave for input, nothing done"
				return $""
			endif
			Wave corners=$name
			Wave vcTemp = str2vec(point)
			pointOK = numpnts(vcTemp)==3
		else
			return $""						// cannot find corners wave, give up
		endif
		printIt = 1							// force printing, since no input wave was provided
	endif
	if (!WaveExists(corners))
		return $""
	endif

	Make/N=3/D/FREE lo, hi				// lowest and highest values in x,y,z
	Make/N=(DimSize(corners,0))/D/FREE v=corners[p][0]
	lo[0] = WaveMin(v) ;	hi[0] = WaveMax(v)
	v = corners[p][1]
	lo[1] = WaveMin(v) ;	hi[1] = WaveMax(v)
	v = corners[p][2]
	lo[2] = WaveMin(v) ;	hi[2] = WaveMax(v)
	if (numpnts(vcTemp)==3)
		Wave vc = vcTemp
	else
		MatrixOP/FREE vc = (lo+hi)/2	// center
	endif
	if (norm(vc)==0)
		return $""							// center is origin, cannot draw a line from origin to origin
	endif

	name=CleanupName(NameOfWave(corners)+"_radialLine",0)
	Make/N=(2,3)/O $name/WAVE=line = NaN
	String wnote=ReplaceStringByKey("waveClass",note(corners),"GizmoRadialLine,GizmoPath","=")
	Note/K line, wnote
	if (printIt)
		printf "Created a radial line named  '%s'  in folder  '%s'\r",NameOfWave(line),GetWavesDataFolder(line,1)
		print "To add this to a Gizmo, add a new Path with this wave as the source"
	endif
	line[1][] = vc[q]						// one end of line always lies at center of Gizmo
	if (lo[0]<=0 && hi[0]<=0 && lo[1]<=0 && hi[1]<=0 && lo[2]<=0 && hi[2]<=0)
		line[0][] = 0						// origin is inside box, so start of line is {0,0,0}
	else
		Make/N=3/D/FREE towardZero, ps
		towardZero = vc[p]>0 ? lo[p] : hi[p]	// limit that projects towards origin
		ps = towardZero/vc
		Variable psMin = WaveMax(ps)	// closest to 1, when ps==1, line has zero length
		line[0][] = psMin*vc[q]
	endif
	return line
End

//Function Make_gLinePath()
//	Wave qvec = BraggPeak()
//
//	Wave corners=QspaceCorners
//	STRUCT boundingVolume v
//	v.xlo = corners[0][0] ;	v.xhi = corners[1][0]
//	v.ylo = corners[0][1] ;	v.yhi = corners[2][1]
//	v.zlo = corners[0][2] ;	v.zhi = corners[4][2]
//	updateBoundingVolumeStruct(v)
//
//	if (!insideBoundingVolumeStruct(v,qvec))
//		print "Bragg point is not inside of volume"
//	endif
//
//	Variable flox = (v.xlo)/qvec[0], fhix = (v.xhi)/qvec[0], swap
//	if (flox>fhix)
//		swap = flox
//		flox = fhix
//		fhix = swap
//	endif
//	Variable floy = (v.ylo)/qvec[1], fhiy = (v.yhi)/qvec[1]
//	if (floy>fhiy)
//		swap = floy
//		floy = fhiy
//		fhiy = swap
//	endif
//	Variable floz = (v.zlo)/qvec[2], fhiz = (v.zhi)/qvec[0]
//	if (floz>fhiz)
//		swap = floz
//		floz = fhiz
//		fhiz = swap
//	endif
//
//	Variable flo = max(max(flox,floy),floz)
//	Variable fhi = min(min(fhix,fhiy),fhiz)
//	Make/N=(2,3)/O gLinePath
//	gLinePath[0][]= flo*qvec[q]
//	gLinePath[1][]= fhi*qvec[q]
//End

//  ============================= End of QGizmo Enhancements =============================  //
//  ======================================================================================  //



//  ======================================================================================  //
//  ========================= Start of Bounding Volume Structure =========================  //

Structure boundingVolume			// a generic bounding volume, probably in k-space
	double	xlo, xhi					// range of x, these 6 form the corners of a box
	double	ylo, yhi					// range of y
	double	zlo, zhi					// range of z
	double	xW, yW, zW				// box size is xW x yW x zW
	int16		Nx, Ny, Nz				// dimensions of the array
	double	dx, dy, dz				// step size along box
	double	vol						// volume of box = xW*yW*zW
EndStructure
//Structure boundingVolume				// a generic bounding volume, probably in k-space
//	Variable xlo,xhi, ylo,yhi,zlo,zhi		// corners of a box
//	Variable dx,dy,dz					// box size is dx x dy x dz
//	Variable vol
//EndStructure
//
//ThreadSafe Function initBoundingVolumeStruct(v)
Function initBoundingVolumeStruct(v)
	STRUCT boundingVolume &v
	v.xlo = Inf  ;	v.ylo = Inf ;	v.zlo = Inf
	v.xhi = -Inf ;	v.yhi = -Inf ;	v.zhi = -Inf
	v.Nx = 0     ;	v.Ny = 0 ;		v.Nz = 0			// default dimension
	updateBoundingVolumeStruct(v)
End
//
//ThreadSafe Function updateBoundingVolumeStruct(v)
Function updateBoundingVolumeStruct(v)
	STRUCT boundingVolume &v
	v.xW = abs(v.xhi - v.xlo)
	v.yW = abs(v.yhi - v.ylo)
	v.zW = abs(v.zhi - v.zlo)
	v.vol = v.xW * v.yW * v.zW

	if (v.Nx > 1)					// prefer to use Nx
		v.dx = v.xW / (v.Nx - 1)
	elseif (v.dx > 0)				// N not good, try to use dx
		v.Nx = ceil(v.xW / v.dx) + 1
	endif

	if (v.Ny > 1)					// prefer to use Ny
		v.dy = v.yW / (v.Ny - 1)
	elseif (v.dy > 0)				// N not good, try to use dy
		v.Ny = ceil(v.yW / v.dy) + 1
	endif

	if (v.Nz > 1)					// prefer to use Nz
		v.dz = v.zW / (v.Nz - 1)
	elseif (v.dz > 0)				// N not good, try to use dz
		v.Nz = ceil(v.zW / v.dz) + 1
	endif

	v.dx = v.Nx>1 ? v.xW / (v.Nx - 1) : 0		// re-update the dx
	v.dy = v.Ny>1 ? v.yW / (v.Ny - 1) : 0
	v.dz = v.Nz>1 ? v.zW / (v.Nz - 1) : 0
End
//
//ThreadSafe Function/S boundingVolumeStruct2str(v)
Function/S boundingVolumeStruct2str(v)
	STRUCT boundingVolume &v
	String str
//	sprintf str,"Vol = X=[%g, %g] ÆX=%g,  Y=[%g, %g] ÆY=%g,  Z=[%g, %g] ÆZ=%g,  Vol=%g\r",v.xlo,v.xhi,v.xW, v.ylo,v.yhi,v.yW, v.zlo,v.zhi,v.zW,v.vol
	sprintf str,"Vol = X=[%g, %g] Nx=%g,  Y=[%g, %g] Ny=%g,  Z=[%g, %g] Nz=%g,  Vol=%g\r",v.xlo,v.xhi,v.Nx, v.ylo,v.yhi,v.Ny, v.zlo,v.zhi,v.Nz,v.vol
	return str
End
//
//ThreadSafe Function copyBoundingVolumeStructs(vin,vout)		// note, vall may be v1 or v2
Function copyBoundingVolumeStructs(vin,vout)		// note, vall may be v1 or v2
	STRUCT boundingVolume &vin, &vout
	vout.xlo = vin.xlo ;	vout.ylo = vin.ylo ;	vout.zlo = vin.zlo
	vout.xhi = vin.xhi ;	vout.yhi = vin.yhi ;	vout.zhi = vin.zhi
	vout.xW = vin.xW   ;	vout.yW = vin.yW   ;	vout.zW = vin.zW
	vout.Nx = vin.Nx   ;	vout.Ny = vin.Ny   ;	vout.Nz = vin.Nz
	vout.dx = vin.dx   ;	vout.dy = vin.dy   ;	vout.dz = vin.dz
	vout.vol = vin.vol
End
//
//ThreadSafe Function extendBoundingVolumeStruct(v,vec)
Function extendBoundingVolumeStruct(v,vec)
	// extend v to include the vector vec, maintain the dx,dy,dz when doing this
	STRUCT boundingVolume &v
	Wave vec

	updateBoundingVolumeStruct(v)									// ensure that dx,dy,dz are calculated
	Variable xlo = min(xlo,vec[0]), xhi = max(xhi,vec[0])	// new ranges
	Variable ylo = min(ylo,vec[1]), yhi = max(yhi,vec[1])
	Variable zlo = min(zlo,vec[2]), zhi = max(zhi,vec[2])
	Variable Nx,Ny,Nz
	Nx = ceil(abs(xhi-xlo) / v.dx) + 1							// new Nx,Ny,Nz
	Ny = ceil(abs(yhi-ylo) / v.dy) + 1
	Nz = ceil(abs(zhi-zlo) / v.dz) + 1
	v.Nx = Nx>0 && numtype(Nx)==0 ? Nx : 0
	v.Ny = Ny>0 && numtype(Ny)==0 ? Ny : 0
	v.Nz = Nz>0 && numtype(Nz)==0 ? Nz : 0

	v.xlo = xlo ;		v.xhi = xhi
	v.ylo = ylo ;		v.yhi = yhi
	v.zlo = zlo ;		v.zhi = zhi
	updateBoundingVolumeStruct(v)
End
//
//ThreadSafe Function insideBoundingVolumeStruct(v,vec)		// returns TRUE if vec is inside v
Function insideBoundingVolumeStruct(v,vec)		// returns TRUE if vec is inside v
	STRUCT boundingVolume &v
	Wave vec
	Variable inside = v.xlo <= vec[0] && vec[0] <= v.xhi
	inside = inside && (v.ylo <= vec[1] && vec[1] <= v.yhi)
	inside = inside && (v.zlo <= vec[2] && vec[2] <= v.zhi)
	return inside
End
//
//ThreadSafe Function combineBoundingVolumeStructs(v1,v2,vall)	// note, vall may be v1 or v2
Function combineBoundingVolumeStructs(v1,v2,vall)	// note, vall may be v1 or v2
	STRUCT boundingVolume &v1, &v2, &vall

	updateBoundingVolumeStruct(v1)
	updateBoundingVolumeStruct(v2)

	Variable xlo=min(v1.xlo, v2.xlo), xhi=max(v1.xhi, v2.xhi)
	Variable ylo=min(v1.ylo, v2.ylo), yhi=max(v1.yhi, v2.yhi)
	Variable zlo=min(v1.zlo, v2.zlo), zhi=max(v1.zhi, v2.zhi)
	Variable dx=max(v1.dx, v2.dx), dy=max(v1.dy, v2.dy), dz=max(v1.dz, v2.dz)
	Variable Nx,Ny,Nz
	Nx = ceil(abs(xhi-xlo) / dx) + 1						// new Nx,Ny,Nz
	Ny = ceil(abs(yhi-ylo) / dy) + 1
	Nz = ceil(abs(zhi-zlo) / dz) + 1
	vall.Nx = Nx>0 && numtype(Nx)==0 ? Nx : 0
	vall.Ny = Ny>0 && numtype(Ny)==0 ? Ny : 0
	vall.Nz = Nz>0 && numtype(Nz)==0 ? Nz : 0
	vall.xlo = xlo ;		vall.xhi = xhi
	vall.ylo = ylo ;		vall.yhi = yhi
	vall.zlo = zlo ;		vall.zhi = zhi
	updateBoundingVolumeStruct(vall)
End

//  ========================== End of Bounding Volume Structure ==========================  //
//  ======================================================================================  //

