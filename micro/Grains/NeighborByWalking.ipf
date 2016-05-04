#pragma rtGlobals=1		// Use modern global access method.
#include "GrainMath"
#pragma version = 1.01

Function findGrainsByComparingNeighbors(maxGrains,da,useEdges)		// find grains by looking at voxel to voxel rotationse
	Variable maxGrains					// maximum number of grains to find (if you put in Inf it wil find all)
	Variable da							// threshold angle between two grains to trigger a grain boundary (deg)
	Variable useEdges					// when finding parts of a grain, used grains with only an edge in common (always use faces)
	if (numtype(maxGrains)==2 || maxGrains<1 || numtype(da) || da<=0 || numtype(useEdges) | useEdges<0)
		maxGrains = (numtype(maxGrains) || maxGrains<1) ? 10 : maxGrains
		da = (numtype(da) || da<=0) ? 1 : da
		useEdges = (numtype(useEdges) || useEdges<0) ? 1 : useEdges+1
		Prompt maxGrains, "maximum number of grains to find"		
		Prompt da, "min angle between two grains for a grain boundary"
		Prompt useEdges, "used both edges and faces when extending a grain",popup,"only faces;edges and faces"
		DoPrompt "grains",maxGrains,da,useEdges
		if (V_flag)
			return 0
		endif
		useEdges -= 1
		printf "    findGrainsByComparingNeighbors(%g,%g,%d)\r",maxGrains,da,useEdges
	endif
	if (numtype(maxGrains)==2 || maxGrains<1 || numtype(da) || da<=0 || numtype(useEdges) | useEdges<0)
		String str
		sprintf str, "invalid inputs, maxGrains=%g,  da=%g, useEdges=$g",maxGrains,da,useEdges
		Abort str
	endif

	Wave indicies=:raw:indicies			// for each point, provides ix,iy,iz into grains[][][]
	Wave XYZindex=:raw:XYZindex		// 3d array with single index for each point
	Variable N=DimSize(indicies,0)
	Make/N=(N)/O totalAngle
	totalAngle = anglefromMatrix(p)	// totalAngle is angle from reference orientation to each grain (ref is no rotation)
	WaveStats/Q totalAngle
	Variable Nmax=V_npnts
	maxGrains = min(Nmax,maxGrains)// cannot have more grains than there are valid voxels
	Variable maxAngle,deltaAngle,minA,maxA,i
	Variable/C tempZ
	String wName = "gNeighbor"
	Duplicate/O XYZindex $wName		// a 3d space to hold the grain numbers called grains
	Wave grains = $wName
	grains = NaN
	String noteStr = note(grains)
	noteStr = ReplaceNumberByKey("dAngle",noteStr,da,"=")
	Note/K grains
	Note grains,noteStr

	Make/N=(maxGrains)/O $(wName+"_Vol")		// volume of grain (voxels)
	Wave gVol=$(wName+"_Vol")
	gVol = NaN
	Variable Ngrains=0
	Variable Nixyz
	Variable j, ix,iy,iz, Nadded,Ntotal=0
	Variable tick0=ticks
	do
		tempZ = MaxHistoramTotalAngle(totalAngle)
		Wave HistTotalAngle=HistTotalAngle 
		if (strlen(GraphWithWaveName("HistTotalAngle")))
			sprintf str, "ValDisplay sumOfHistogram value=#\"%g\",win=%s",sum(HistTotalAngle),GraphWithWaveName("HistTotalAngle")
			Execute str
		endif
		DoUpdate
		maxAngle = real(tempZ)
		deltaAngle = abs(imag(tempZ))
		minA = maxAngle-deltaAngle/2
		maxA = maxAngle+deltaAngle/2
		for (i=0;i<N;i+=1)						// find first point in totalAngle that has value of maxAngle (± deltaAngle/2)
			if (minA<totalAngle[i] && totalAngle[i]<maxA)
				break
			endif
		endfor
		if (i>=N)								// deltaAngle was too tight, loosen tolerance and search again
			minA = maxAngle-deltaAngle
			maxA = maxAngle+deltaAngle
			for (i=0;i<N;i+=1)					// find first point in totalAngle that has value of maxAngle (± deltaAngle/2)
				if (minA<totalAngle[i] && totalAngle[i]<maxA)
					break
				endif
			endfor
		endif
		if (i>=N)								// deltaAngle was too tight, loosen tolerance and search again
			Abort "unable to find an angle that was in the histogram??!!"
		endif

		// identify a list of all points in grains[][][] that are connected to the peak (by a surface, not just edge)
		Make/N=(1,3)/O ixyz					// list of indicies into grains for the particle
		ixyz[0][0] = indicies[i][0]				// put in the starting point of the particle
		ixyz[0][1] = indicies[i][1]
		ixyz[0][2] = indicies[i][2]

		// get new particle, points returned in ixyz, and Nixyz is the number found
		Nixyz = FindParticleFromOrientsMat(ixyz,da,useEdges)	// ixyz[Nixyz] are the points in grains to consider
if (Ngrains<debugN)
//	printf "found a 'grain' in grains, it consits of %d elements in grains[][][]\r",Nixyz
	String degree = StrVarOrDefault("root:Packages:Grains:degree","¡")
	printf "found grain %d in grains with totalAngle=%g%s , it consits of %d elements in grains[][][]\r",Ngrains,totalAngle[i],degree,Nixyz
endif
		// put the Nixyz points in ixyz into grains[][][]
		for (i=0;i<Nixyz;i+=1)
			ix = ixyz[i][0]
			iy = ixyz[i][1]
			iz = ixyz[i][2]
			j = XYZindex[ix][iy][iz]				// the single index that correspoinds to (ix,iy,iz)
			totalAngle[j] = NaN						// delete the corresponding point from totalAngle
			grains[ix][iy][iz] = Ngrains
		endfor
		gVol[Ngrains] = Nixyz
		Ngrains += min(Nixyz,1)
		Ntotal += Nixyz
		if (Ntotal>12000)
			print Ntotal,numpnts(totalAngle)
		endif
	while (Nixyz && Ngrains<maxGrains && Ntotal<Nmax)
	Redimension/N=(Ngrains) gVol
	DoUpdate

	// re-order grain number in order of increasing grain size, ie grain 0 is the largest
	Duplicate/O gVol, temp_grainNum_
	temp_grainNum_ = p
	Sort/R gVol,gVol,temp_grainNum_
	//	next re-label grains[][][] with new grain numbers, use -1 as a flag for the grain to be relabeled
	grains *= -1										// this flags all grains as being wrong (except zero)
	for (i=0;i<Ngrains;i+=1)
		j =  locateValueInWave(temp_grainNum_,i)		// new grain number, (i is old grain number)
		grains = (grains[p][q][r]==-i) ? j : grains[p][q][r]
	endfor

	MakeAllGrainStatistics(grains)
	makeGrainsColorTable(grains)

	if (ItemsInList(GetRTStackInfo(0))<2)
		WaveStats/Q indicies
		Variable validOriginal = V_npnts/3			// number of valid voxels in original data
		WaveStats/Q grains
		Variable valid = V_npnts					// number of valid voxels assigned to a grain
		WaveStats/Q totalAngle
		Variable leftOver = V_npnts					// number of voxels that were indexed but not assigned
		printf "created grains in    %s[%d][%d][%d]\r",NameOfWave(grains),DimSize(XYZindex,0),DimSize(XYZindex,1),DimSize(XYZindex,2)
		printf "identified Ngrains = %d containing %d voxels,  there are %g valid voxels left over (identified %.1f%)\r",Ngrains,valid,leftOver,valid/validOriginal*100
//		printf "execution took %s\r", num2sexigesmal((ticks-tick0)/60.15,1)
		printf "execution took %s\r", Secs2Time((ticks-tick0)/60.15,5,1)
	endif
	//	now grains[][][] contains the identifier for each grain, and XYZindex[][][] contains pointers to their source
	KillWaves/Z totalAngle, ixyz, temp_grainNum_
	return Ngrains
End
//Function doubleCheckIndividuals(da)
//	Variable da							// threshold angle between two grains to trigger a grain boundary (deg)
//	Wave gVol=gVol
//	Wave XYZindex=:raw:XYZindex		// 3d array with single index for each point
//	Wave grains=grains
//	Variable Nx,Ny,Nz, N=numpnts(gVol)
//	Nx=DimSize(XYZindex,0)
//	Ny=DimSize(XYZindex,1)
//	Nz=DimSize(XYZindex,2)
//
//	// offsets[18][3] contains the 18 voxels to the 1st and 2nd nearest neighbors (faces and edges)
//	Make/N=(18,3) offsets
//	offsets[0][0]= {-1,1,  0,0,  0,0,	-1,1,-1,1,  0,0,	-1,  1,-1,  1,  0,  0}
//	offsets[0][1]= {  0,0,-1,1,  0,0,	-1,1,  0,0,-1,1,	  1,-1,  0,  0,-1,  1}
//	offsets[0][2]= {  0,0,  0,0,-1,1,	  0,0,-1,1,-1,1,	  0,  0,  1,-1,  1,-1}
//	Variable Noffset=DimSize(offsets,0)
//
//	Variable i,j, ix,iy,iz
//	Variable ix0,iy0,iz0
//	Variable i1,i2,ang, grainNum
//	for (i=0;i<N;i+=1)							// loop over every singleton
//		if (gVol[i]>1)
//			continue
//		endif
//		if (grainNum2Indicies(i,ix0,iy0,iz0))	// get ix,iy,iz for grain i
//			DoAlert 0, "could not find indicies for grain "+num2istr(i)
//			continue							// error, could not find indicies for this grain number
//		endif
//		i1 = XYZindex[ix0][iy0][iz0]
//		for (j=0;j<Noffset;j+=1)				// loop through the offsets
//			ix = ix0 + offsets[j][0]
//			iy = iy0 + offsets[j][1]
//			iz = iz0 + offsets[j][2]
//			if (numtype(ix+iy+iz) || ix>=Nx || iy>=Ny || iz>=Nz || ix<0 || iy<0 || iz<0)
//				continue						// point is NaN or out of range
//			endif
//			i2 = XYZindex[ix][iy][iz]
//			// compute the angle between ixyzAdd[i] and checkNext[j][3]
//			if (abs(angleBetweenPoints(i1,i2))<da)	// angle is smal enough
//				ang = angleBetweenPoints(i1,i2)
//				grainNum = grains[ix][iy][iz]
//				printf "found that grain%d at [%d][%d][%d] is close to grain%d at [%d][%d][%d],  with angle=%g¡\r",i,ix0,iy0,iz0,grainNum,ix,iy,iz,ang
//				// add this point to ixyz[N][3] and ixyzNEW[iNew][3]
//				// attach i1 to grain
////				break
//			endif
//		endfor									// end of loop of j =[0,Noffset]
//	endfor
//	KillWaves/Z offsets
//End
//Static Function grainNum2Indicies(grainNum,ix,iy,iz)
//	Variable grainNum
//	Variable &ix,&iy,&iz
//
//	Wave grains=grains
//	Variable Nx,Ny,Nz
//	Nx=DimSize(grains,0)
//	Ny=DimSize(grains,1)
//	Nz=DimSize(grains,2)
//	for (ix=0;ix<Nx;ix+=1)
//		for (iy=0;iy<Ny;iy+=1)
//			for (iz=0;iz<Nz;iz+=1)
//				if (grains[ix][iy][iz]==grainNum)
//					return 0
//				endif
//			endfor
//		endfor
//	endfor
//	return 1
//End
Static Function/C MaxHistoramTotalAngle(totalAngle)		// returns the angle for the maximum of the histogram
	// the angle of the max is in the real part, and the bin width for the histogram is the imaginary part
	Wave totalAngle
	WaveStats/Q/M=1 totalAngle			// needed to get number of points and range for histogram
	Make/N=(1.5*sqrt(V_npnts+V_numNans))/D/O HistTotalAngle
	SetScale/I x V_min,V_max,"", HistTotalAngle
	Histogram/B=2 totalAngle,HistTotalAngle	// put histogram into HistTotalAngle[]
	WaveStats/M=1/Q HistTotalAngle		// needed to find the peak of HistTotalAngle[]
	Variable/C val = cmplx(V_maxloc,DimDelta(HistTotalAngle,0))
	KillWaves/Z HistTotalAngle
	return val
End
// NOTE, capitalization is important!!
Static Function/S GraphWithWaveName(waveNameStr)	// Find the topmost graph that shows waveNameStr
	String waveNameStr							// name of wave to search for
	String wList = WinList("Graph*",";","WIN:1")
	String graphNameStr

	Variable i
	for (i=0, graphNameStr=StringFromList(0,wList); strlen(graphNameStr); i+=1, graphNameStr = StringFromList(i,wList))
		if (WhichListItem(waveNameStr,TraceNameList(graphNameStr,";",1))>=0)
			break
		endif
	endfor
	return graphNameStr
End

// search for contiguous grains based on a grain to grain angle change of da
Function FindParticleFromOrientsMat(ixyz,da,useEdges)
	Wave ixyz				// returned list of indicies that are adjacent, first point is starting point
	Variable da				// angle threshold, is angles changes by more  than da, a grain boundary
	Variable useEdges		// when finding parts of a grain, used grains with only an edge in common (always use faces)

	Wave OrientMat = :raw:OrientMat	// input waves that are needed
	Wave XYZindex = :raw:XYZindex
	Wave indicies = :raw:indicies
	if (!WaveExists(OrientMat) || !WaveExists(XYZindex) ||!WaveExists(indicies))
		DoAlert 0, "FindParticleFromOrientsMat() could not find data arrys OrientMat, XYZindex, or indicies"
		return 0
	endif
	if (WaveDims(OrientMat)!=3 || WaveDims(XYZindex)!=3 || WaveDims(indicies)!=2)
		DoAlert 0, "FindParticleFromOrientsMat() OrientMat, XYZindex, or indicies has wrong dimensions"
		return 0
	endif
	Variable N=1			// number of valid voxels adjacent in grain, (also number of valid points in ixyz)
	Variable Nx,Ny,Nz
	Nx=DimSize(XYZindex,0)
	Ny=DimSize(XYZindex,1)
	Nz=DimSize(XYZindex,2)

	String str
	// ixyzNEW is a list of voxels that were added in the last pass, need to check voxels about these
	str=UniqueName("ixyzNEW",1,0)			// new points added to ixyz[][3] that have not been checked
	Make/N=(100,3) $str
	Wave ixyzNEW=$str
	ixyzNEW[0][0] = ixyz[0][0]				// init ixyzNew[][3] to be the starting point
	ixyzNEW[0][1] = ixyz[0][1]
	ixyzNEW[0][2] = ixyz[0][2]

	str=UniqueName("offsets",1,0)
	if (useEdges)
		// offsets[][] contains the 18 voxels to the 1st and 2nd nearest neighbors (faces and edges)
		Make/N=(18,3) $str
		Wave offsets=$str
		offsets[0][0]= {-1,1,  0,0,  0,0,	-1,1,-1,1,  0,0,	-1,  1,-1,  1,  0,  0}
		offsets[0][1]= {  0,0,-1,1,  0,0,	-1,1,  0,0,-1,1,	  1,-1,  0,  0,-1,  1}
		offsets[0][2]= {  0,0,  0,0,-1,1,	  0,0,-1,1,-1,1,	  0,  0,  1,-1,  1,-1}
	else
		// offsets[][] contains the 6 voxels to the 1st nearest neighbors (faces in common)
		Make/N=(6,3) $str
		Wave offsets=$str
		offsets[0][0]= {-1,1,  0,0,  0,0}
		offsets[0][1]= {  0,0,-1,1,  0,0}
		offsets[0][2]= {  0,0,  0,0,-1,1}
	endif

	Variable iNew=1, Ncheck=0, iAdd
	Variable Noffset = DimSize(offsets,0)
	Variable i,j, ix,iy,iz, i1,i2
	do
		Duplicate/O ixyzNEW, ixyzAdd
		iAdd = iNew
		ixyzNEW = NaN
		iNew = 0
		// check the neighborhood of each ixyzAdd[][3] individually
		for (i=0;i<iAdd;i+=1)						// loop over the ixyzAdd[][3]
			i1 = XYZindex[ixyzAdd[i][0]][ixyzAdd[i][1]][ixyzAdd[i][2]]
			for (j=0;j<Noffset;j+=1)				// loop through the offsets
				ix = ixyzAdd[i][0] + offsets[j][0]
				iy = ixyzAdd[i][1] + offsets[j][1]
				iz = ixyzAdd[i][2] + offsets[j][2]
				if (numtype(ix+iy+iz) || ix>=Nx || iy>=Ny || iz>=Nz || ix<0 || iy<0 || iz<0)
					continue						// point is NaN or out of range
				endif
				if  (point_in_list(ix,iy,iz,ixyz,N))
					continue						// point already in ixyz[]
				endif
				i2 = XYZindex[ix][iy][iz]
				// compute the angle between ixyzAdd[i] and ixyzAdd[i]offsets[j]
				if (abs(angleBetweenPointsSymReduced(i1,i2))<da)	// angle is smal enough
					// add this point to ixyz[N][3] and ixyzNEW[iNew][3]
					if (DimSize(ixyz,0)<(N+1))			// need to extend ixyz
						Redimension/N=(DimSize(ixyz,0)+100,3) ixyz
						ixyz[N,Inf][] = NaN
					endif
					if (DimSize(ixyzNEW,0)<(iNew+1))	// need to extend ixyzNEW
						Redimension/N=(DimSize(ixyzNEW,0)+100,3) ixyzNEW
					endif
					ixyz[N][0] = ix
					ixyz[N][1] = iy
					ixyz[N][2] = iz
					N += 1
					ixyzNEW[iNew][0] = ix
					ixyzNEW[iNew][1] = iy
					ixyzNEW[iNew][2] = iz
					iNew += 1
				endif
			endfor									// end of loop of j =[0,Noffset]
		endfor										// end of loop of i =[0, iAdd]
	while(iNew>0)
	Redimension/N=(N,3) ixyz			// trim to exact size

	//	now ixyz[N][3] are the list of indicies into in3d[][][] that form a contiguous particle
	KillWaves/Z offsets, ixyzNEW, ixyzAdd
	return N
End



//  ===================================================================
//  ===================================================================
//  ===================================================================
//		Utility routines
//  ===================================================================


Function/S makeRange3dGrain(grains,range)
	Wave grains					// 3d wave with all of the grain numbers (e.g. gNeighbor)
	String range					// list of grain numbers, possibly non-contiguous range, ie "1-20,25,30-40"
	if (!WaveExists(grains) || WaveDims(grains)!=3 || numtype(NextInRange(range,-inf)))
		String grainName=NameOfWave(grains)
		Prompt grainName,"3d wave with all grains",popup,ListMultiGrainMats()
		Prompt range,"range of grain numbers"
		DoPrompt "range of grains to make",grainName,range
		if (V_flag)
			return ""
		endif
		Wave grains = $grainName
		if (ItemsInList(GetRTStackInfo(0))<2)
			printf "     makeRange3dGrain(%s,\"%s\")\r",NameOfWave(grains),range
		endif
	endif
	if (!WaveExists(grains) || WaveDims(grains)!=3 || numtype(NextInRange(range,-inf)))
		return ""
	endif

	String fullPath, list = "" 						// list of individual grain arrays, one array for each grain
	range = expandRange(range,",")

	Variable i,j, N = ItemsInRange(range)
	for (i=0;i<N;i+=1)
		j = str2num(StringFromList(i,range,","))
		fullPath = makeSingle3dGrain(grains,j)
		if (exists(fullPath)==1)
			list += fullPath+";"
		else
			print "ERROR,  in makeRange3dGrain() for grain no. ",j,"   unable to make the grain array"
		endif
	endfor
	return list
End

Function/S makeSingle3dGrain(grains,grainNum)
	Variable grainNum				// find the grains that border grainsNum
	Wave grains					// 3d wave with all of the grain numbers (e.g. gNeighbor)
//	if (!WaveExists(grains) || numtype(grainNum) || grainNum<0)
//		String grainName=NameOfWave(grains)
//		Prompt grainName,"3d wave with all grains",popup,ListMultiGrainMats()
//		grainNum = numtype(grainNum)==2 ? 0 : max(grainNum,0)
//		Prompt grainNum, "grain number"
//		DoPrompt "choose grain",grainName,grainNum
//		if (V_flag)
//			return ""
//		endif
//		Wave grains = $grainName
//	endif
	if (!WaveExists(grains))
		return ""
	endif
	WaveStats/Q grains
	Variable maxGrain = V_max, totalVol=V_npnts
	if (numtype(grainNum) || grainNum<0 || grainNum>maxGrain)
		return ""
	endif
	String gName = GetWavesDataFolder(grains,2)

	Wave grain_ixlo=$(gName+"_ixlo")	// total extent of grain in [ix][iy][iz] indicies
	Wave grain_ixhi=$(gName+"_ixhi")
	Wave grain_iylo=$(gName+"_iylo")
	Wave grain_iyhi=$(gName+"_iyhi")
	Wave grain_izlo=$(gName+"_izlo")
	Wave grain_izhi=$(gName+"_izhi")
	Variable Nx,Ny,Nz, ix0, iy0, iz0
	ix0 = grain_ixlo[grainNum]
	iy0 = grain_iylo[grainNum]
	iz0 = grain_izlo[grainNum]
	Nx = grain_ixhi[grainNum] - ix0 + 3	// two extra points, one at each end
	Ny = grain_iyhi[grainNum] - iy0 + 3
	Nz = grain_izhi[grainNum] - iz0 + 3
	String noteStr = note(grains)
	noteStr = ReplaceStringByKey("grainsWave",noteStr,GetWavesDataFolder(grains,2),"=")
	Wave grain_xc=$(gName+"_xc")				// scaled postiion of the COM
	Wave grain_yc=$(gName+"_yc")
	Wave grain_zc=$(gName+"_zc")
	Wave grain_COM=$(gName+"_COM")			// single index to the center of mass
	Wave grain_Vol=$(gName+"_Vol")				// size of grain (voxels)
	noteStr = ReplaceNumberByKey("grainNum",noteStr,grainNum,"=")
	noteStr = ReplaceNumberByKey("COM_single",noteStr,grain_COM[grainNum],"=")
	noteStr = ReplaceNumberByKey("volume",noteStr,grain_Vol[grainNum],"=")
	noteStr = ReplaceNumberByKey("grain_xc",noteStr,grain_xc[grainNum],"=")
	noteStr = ReplaceNumberByKey("grain_yc",noteStr,grain_yc[grainNum],"=")
	noteStr = ReplaceNumberByKey("grain_zc",noteStr,grain_zc[grainNum],"=")

	Variable x0,y0,z0,dx,dy,dz
	dx = DimDelta(grains,0)
	dy = DimDelta(grains,1)
	dz = DimDelta(grains,2)
	x0 = ix0*dx + DimOffset(grains,0) - dx		// (x0,y0,z0) is xyz of the low corner of the grain
	y0 = iy0*dy + DimOffset(grains,1)	- dy
	z0 = iz0*dz + DimOffset(grains,2) - dz
	noteStr = ReplaceNumberByKey("x0",noteStr,x0,"=")
	noteStr = ReplaceNumberByKey("y0",noteStr,y0,"=")
	noteStr = ReplaceNumberByKey("z0",noteStr,z0,"=")
//print "for grain 6, corner at ",x0,y0,z0,"(µm)    and delta's are",dx,dy,dz, "(µm)"
//printf "ranges are x=[%g,%g]    x=[%g,%g]    x=[%g,%g] (µm)\r",x0,x0+(Nx-1)*dx,y0,y0+(Ny-1)*dy,z0,z0+(Nz-1)*dz
	String wName = gName+"_"+num2istr(grainNum)
	Make/N=(Nx,Ny,Nz)/O $wName
	Wave grain = $wName
	String micron = StrVarOrDefault("root:Packages:Grains:micron","micron")
	SetScale/P x x0,dx,micron, grain
	SetScale/P y y0,dy,micron, grain
	SetScale/P z z0,dz,micron, grain
	grain = NaN		// the backgound
	grain[1,Nx-2][1,Ny-2][1,Nz-2] = grains[p+ix0-1][q+iy0-1][r+iz0-1]
	grain = (grain[p][q][r]!=grainNum) ? grainNum-1 : grainNum

	String str
	Wave grainColors = $(gName+"_Colors")
	if (!WaveExists(grainColors))
		str = ""
		sprintf str, "'%s_Colors' does not exist, try issuing:  makeGrainsColorTable(%s)",gName,gName
		DoAlert 0, str
		print "     ",str
	Endif
	Wave gVol=gNeighbor_Vol
	if (WaveExists(grainColors) && WaveExists(gVol))
		WaveStats/Q grainColors
		Variable scaleFactor = (V_max<1.1) ? 1 : 65535
		Variable vol = grainNum==0 ? 0 : sum(gVol,0,grainNum-1)	// volume before grain i
		Variable alpha = 1
		alpha = (vol<=0.85*totalVol) ? 0.7 : alpha	// set alpha on the volume, (big grains translucent, small opaque)
		alpha = (vol<=0.75*totalVol) ? 0.5 : alpha
		alpha = (vol<=0.5*totalVol) ? 0.3 : alpha
//		sprintf str, "%g,%g,%g,%g",grainColors(grainNum)[0]/65535,grainColors(grainNum)[1]/65535,grainColors(grainNum)[2]/65535,alpha
		sprintf str, "%g,%g,%g,%g",grainColors[grainNum][0]/scaleFactor,grainColors[grainNum][1]/scaleFactor,grainColors[grainNum][2]/scaleFactor,alpha
		noteStr = ReplaceStringByKey("RGBA",noteStr,str,"=")
	endif

	noteStr = ReplaceNumberByKey("isoValue",noteStr,grainNum-.5,"=")
	Note/K grain
	Note grain, noteStr
	//	print "ix0, iy0, iz0 = ",ix0, iy0, iz0
	return GetWavesDataFolder(grain,2)
End

