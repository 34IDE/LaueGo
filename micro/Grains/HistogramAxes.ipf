#pragma rtGlobals=1		// Use modern global access method.
#pragma version = 1.01


Menu "Macros"
	"-"
	"(Old methods"
	"  Calculate Grains from Axes Histogram",findGrainsByAxesHist(NaN)
End


// Histograming Rodriguiez space is probably a bad idea.  Stick to the total rotation angle.
//	so use the (Hx,Hy,Hz)(Kx,Ky,Kz)(Lx,Ly,Lz) matrix


//  ===================================================================
//  ===================================================================
//  ===================================================================
//		Most of the following functions has to do with histograming Rodriguiez space.
//  ===================================================================

Function  findGrainsByAxesHist(maxGrains)		// find the grains by looking at a histogram of the axes
	Variable maxGrains			// maximum number of grains to find
	if (numtype(maxGrains) || maxGrains<1)
		maxGrains = 10
		Prompt maxGrains, "maximum number of grains to find"		
		DoPrompt "max no. of grains",maxGrains
		if (V_flag)
			return 0
		endif
	endif

	Wave axisX=:raw:axisX,axisY=:raw:axisY,axisZ=:raw:axisZ
	Wave indicies=indicies
	Variable N=numpnts(axisX)
	Make/N=(N,3)/O axes				// make a single matrix to hold the 3 axisN components
	axes[][0] = axisX[p]				// rotaion axis (and magnitude) for each measured voxel
	axes[][1] = axisY[p]
	axes[][2] = axisZ[p]

	Wave XYZindex=:raw:XYZindex
	Duplicate/O XYZindex grains			// a 3d space to hold the grain numbers called grains
	grains = NaN
	Variable Ngrains=0
	Variable delta						// scaled step size in Hist3dAxes[][][] for all 3 dimensions

	STRUCT peak3D peak				// holds peak result from Histogram3dAxisSpace()
	String str
	Variable Nixyz
	STRUCT boxStruct box
	Variable i,j, ix,iy,iz, Nadded,Ntotal=0
	Make/N=(maxGrains)/O Nincrement
	Nincrement = NaN
if (debugN)
WaveStats/Q axes
printf "before starting axes has %g valid triplets\r",V_npnts/3
endif
	Variable tick0=ticks
	do
// Variable timer = startMSTimer
		str = Histogram3dAxisSpace(axes,peak)		// peak position is at (peak.ix,peak.iy,peak.iz)
//printf "Histogram3dAxisSpace() took %g msec\r",stopMSTimer(timer)/1000
		Wave Hist3dAxes = $StringByKey("HIST",str,"=")	// a 3d histogram of axes[N][3]
if (Ngrains<debugN)
printf "\rFor Ngrains = %d\r",Ngrains
		printf "peak of Hist3dAxes(%g)(%g)(%g) = Hist3dAxes[%d][%d][%d]=%d\r",peak.x,peak.y,peak.z,peak.ix,peak.iy,peak.iz,Hist3dAxes[peak.ix][peak.iy][peak.iz]
endif

		// identify a list of all points in Hist3dAxes that are connected to the peak (by a surface, not just edge)
		Make/N=(1,3)/O ixyz						// list of indicies into Hist3dAxes for the particle
		ixyz[0][0] = peak.ix
		ixyz[0][1] = peak.iy
		ixyz[0][2] = peak.iz
//timer = startMSTimer
		Nixyz = FindParticle3d(Hist3dAxes,ixyz)	// ixyz[Nixyz] are the points in Hist3dAxes to consider
//printf "FindParticle3d() took %g msec\r",stopMSTimer(timer)/1000
if (Ngrains<debugN)
		printf "found a 'grain' in Hist3dAxes, it consits of %d elements in Hist3dAxes[][][]\r",Nixyz
endif

//timer = startMSTimer
		delta = DimDelta(Hist3dAxes,0)
		Nadded = 0
		for (i=0;i<Nixyz;i+=1)
			box.xlo = ixyz[i][0]*DimDelta(Hist3dAxes,0) + DimOffset(Hist3dAxes,0) - delta/2
			box.ylo = ixyz[i][1]*DimDelta(Hist3dAxes,1) + DimOffset(Hist3dAxes,1) - delta/2
			box.zlo = ixyz[i][2]*DimDelta(Hist3dAxes,2) + DimOffset(Hist3dAxes,2) - delta/2
			box.xhi = box.xlo + delta
			box.yhi = box.ylo + delta
			box.zhi = box.zlo + delta
if (Ngrains<debugN)
printf "for ixyz[%d][],  box is:   x=[%g, %g],  y=[%g, %g],  z=[%g, %g]\r",i,box.xlo,box.xhi,box.ylo,box.yhi,box.zlo,box.zhi
endif
			for (j=0;j<N;j+=1)			// search each of the measured axes, is it in box?
				ix = indicies[j][0]
				iy = indicies[j][1]
				iz = indicies[j][2]
				if (!vectorInBox(box,axes[j][0],axes[j][1],axes[j][2]) || numtype(ix+iy+iz))
					continue				// invalid, do not assign
				endif
				if (numtype(axes[j][0]+axes[j][1]+axes[j][2]))
					continue
				endif
				grains[ix][iy][iz] = Ngrains
				axes[j][0] = NaN			// mark this axis as  used
				axes[j][1] = NaN
				axes[j][2] = NaN
				Nadded += 1
			endfor
		endfor
//printf "double loop  took %g msec\r",stopMSTimer(timer)/1000
if (Ngrains<debugN)
WaveStats/Q axes
printf "      now axes has %g valid triplets,     Nadded = %d\r",V_npnts/3,Nadded
endif

		Nincrement[Ngrains] = Nadded
		Ngrains += min(Nadded,1)
		Ntotal += Nadded
		DoUpdate
	while (Nadded && Ngrains<maxGrains)

	if (ItemsInList(GetRTStackInfo(0))<2)
		WaveStats/Q axes
		Variable valid = V_npnts/3
		printf "identified Ngrains = %d containing %d voxels,  there are %g valid voxels left over (identified %.1f%)\r",Ngrains,Ntotal,valid,Ntotal/(Ntotal+valid)*100
//		printf "execution took %s\r", num2sexigesmal((ticks-tick0)/60.15,1)
		printf "execution took %s\r", Secs2Time((ticks-tick0)/60.15,5,1)
	endif

	//	now grains[][][] contains the identifier for each grain, and XYZindex[][][] contains pointers to their source
	KillWaves/Z axes,Hist3dAxes
	return Ngrains
End
Function testCompareAngles(len,testRotation)		// this does a very simpleminded test in Rodriguiz Space
	Variable len
	Variable testRotation	// try 5

	Make/N=3/O a,b,ba		// two vectors to play with
	Variable factor = 180/PI	
	Variable d1,d2
	Variable calledDirectly = (ItemsInList(GetRTStackInfo(0))<2)

	a = {len,0,0}
	b = {len*cos(testRotation*PI/180), len*sin(testRotation*PI/180), 0}
	ba = b-a
	d1 = norm(ba)
	String degree = StrVarOrDefault("root:Packages:Grains:degree","¡")
	if (calledDirectly)
		//	print a,b
		printf "for a %g%s rotation of b about z-axis with |a|=|b|=%g, ^a=%g%s, ^b=%g%s      |b-a| = %g\r",testRotation,degree,norm(a),norm(a)*factor,degree,norm(b)*factor,degree,d1
	endif

	a = {len,0,0}
	b = {len,0,0}
	b *= (factor+testRotation)/factor
	ba = b-a
	d2 = norm(ba)
	Variable error = (d2-d1)/(d2+d1)*2
	if (calledDirectly)
		//	print "",a,b
		printf "for a %g%s extension of b with a||b, ^a=%g%s, ^b=%g%s     |b-a| = %g\r",testRotation,degree,norm(a)*factor,degree,norm(b)*factor,degree,d2
		printf "error between the two methods is %.3f%\r",error*100
	endif
	KillWaves/Z a,b,ba
	return error
End
//Function testManyCompareAngles()
////	printf "%.3f%\r",testCompareAngles(len,testRotation)*100
//	printf "%.3f%\r",testCompareAngles(1,.01)*100
//	printf "%.3f%\r",testCompareAngles(1,.1)*100
//	printf "%.3f%\r",testCompareAngles(1,1)*100
//	printf "%.3f%\r",testCompareAngles(1,10)*100
//	printf "%.3f%\r",testCompareAngles(1,100)*100
//	print ""
//	printf "%.3f%\r",testCompareAngles(.01,5)*100
//	printf "%.3f%\r",testCompareAngles(.1,5)*100
//	printf "%.3f%\r",testCompareAngles(1,5)*100
//	printf "%.3f%\r",testCompareAngles(10,5)*100
//	printf "%.3f%\r",testCompareAngles(100,5)*100
//	//
//	// results are:
//	//	  -0.006%
//	//	  0.001%
//	//	  0.001%
//	//	  0.127%
//	//	  13.013%
//	//    
//	//	  0.032%
//	//	  0.032%
//	//	  0.032%
//	//	  0.032%
//	//	  0.032%
//End



Static Function FindParticle3d(in3d,ixyz)
	Wave in3d				// 3d space to search for contiguous particle (usually Hist3dAxes)
							// in3d denotes holes by NaN (or Infs), finite values are all valid
	Wave ixyz				// returned list of indicies that are adjacent, first point is starting point

	if (numtype(in3d[ixyz[0][0]][ixyz[0][1]][ixyz[0][2]]))
		return 0			// starting point not valid
	endif
	Variable N=1			// number of valid voxels adjacent to peak, (also number of valid points in ixyz)

	Variable Nx,Ny,Nz
	Nx=DimSize(in3d,0)
	Ny=DimSize(in3d,1)
	Nz=DimSize(in3d,2)

	String str=UniqueName("checkNext",1,0)
	Make/N=(600,3) $str
	Wave checkNext=$str
	str=UniqueName("ixyzNEW",1,0)			// new points added to ixyz[][3] that have not been checked
	Make/N=(100,3) $str
	Wave ixyzNEW=$str
	checkNext = NaN
	ixyzNEW[0][0] = ixyz[0][0]				// init ixyzNew[][3] to be just the peak
	ixyzNEW[0][1] = ixyz[0][1]
	ixyzNEW[0][2] = ixyz[0][2]
	Variable iNew=1, Ncheck=0
	Variable i,j, ix,iy,iz
	do
		if (DimSize(checkNext,0)<(6*iNew))		// need to extend checkNext
			Redimension/N=(6*iNew,3) checkNext
		endif
		checkNext = NaN
		for (i=0;i<iNew;i+=1)						// loop over the ixyzNEW[][3]
			for (j=0;j<6;j+=1)
				checkNext[j+i][] = ixyzNEW[i][q]
			endfor
			checkNext[i+0][0] += -1
			checkNext[i+1][0] += 1
			checkNext[i+2][1] += -1
			checkNext[i+3][1] += 1
			checkNext[i+4][2] += -1
			checkNext[i+5][2] += 1
		endfor
		Ncheck = 6*iNew
		Valid_checkNext(checkNext,Ncheck,Nx,Ny,Nz,ixyz,N)	// set invalid points in checkNext to -Inf

		// now check the points in checkNext[Ncheck][3] against in3d (are they valid or NaN's ?)
		// append the valid ones to ixyz[][3] and ixyzNEW[][3]
		ixyzNEW = NaN								// reset ixyzNEW[][3] to empty
		iNew = 0
		for (i=0;i<Ncheck;i+=1)
			ix = checkNext[i][0]
			iy = checkNext[i][1]
			iz = checkNext[i][2]
			if (numtype(ix+iy+iz))
				continue
			elseif (numtype(in3d[ix][iy][iz]))
				continue
			endif
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
		endfor
	while(iNew>0)
	Redimension/N=(N,3) ixyz			// trim to exact size

	//	now ixyz[N][3] are the list of indicies into in3d[][][] that form a contiguous particle
	KillWaves/Z checkNext, ixyzNEW
	return N
End


// go through checkNext[][] and remove points out of range, invalid, or already identifiec
Static Function Valid_checkNext(checkNext,Ncheck,Nx,Ny,Nz,ixyz,Nixyz)
	Wave checkNext				// checkNext[Ncheck][3], triplets to check
	Variable Ncheck				// number of points in checkNext to check
	Variable Nx,Ny,Nz			// valid range of ix=[0,Nx-1],  iy=[0,Ny-1],  iz=[0,Nz-1]
	Wave ixyz					// list of found voxels, do not want duplicates
	Variable Nixyz				// length of ixyz[][3]

	Variable i,j,ix,iy,iz

	// check for duplicate in checkNext[Ncheck][3]
	for (i=0;i<Ncheck;i+=1)
		ix = checkNext[i][0]
		iy = checkNext[i][1]
		iz = checkNext[i][2]
		for (j=i+1;j<Ncheck;j+=1)
			if (ix==checkNext[j][0] && iy==checkNext[j][1] && iz==checkNext[j][2])
				checkNext[i][0] = -Inf			// this flags an invalid point to check
				checkNext[i][1] = -Inf
				checkNext[i][2] = -Inf
			endif
		endfor
	endfor

	for (i=0;i<Ncheck;i+=1)
		ix = checkNext[i][0]
		iy = checkNext[i][1]
		iz = checkNext[i][2]
		// check that each triplet has acceptable range
		if (numtype(ix+iy+iz) || ix>=Nx || iy>=Ny || iz>=Nz || ix<0 || iy<0 || iz<0)
			checkNext[i][0] = -Inf			// this flags an invalid point to check
			checkNext[i][1] = -Inf
			checkNext[i][2] = -Inf
		// check that each triplet has already been identified
		elseif  (point_in_list(ix,iy,iz,ixyz,Nixyz))
			checkNext[i][0] = -Inf			// this flags an invalid point to check
			checkNext[i][1] = -Inf
			checkNext[i][2] = -Inf
		endif
	endfor
End





Function runHist3dAxis()
	STRUCT peak3D peak
	String str = Histogram3dAxisSpace(axes,peak)
	print str
	printf "ix=%g;iy=%g;iz=%g;x=%g;y=%g;z=%g\r",peak.ix,peak.iy,peak.iz,peak.x,peak.y,peak.z
End
// take the set of 3vectors (axisX, axisY, axisZ), and histogram them into a 3d space named Hist3dAxes
// each point in Hist3dAxes represents a single orientation.  So the peak in Hist3dAxes is the biggest grain.
// for each vector (axisX[i], axisY[i], axisZ[i]), the length is the angle of rotation, and the direction
// the axis, so that  angle ~ sqrt(axisX^2+axisY^2+axisZ^2)*180/PI/2 * 1.04386
Function/S Histogram3dAxisSpace(axes,peak)
	Wave axes
	STRUCT peak3D &peak

//Variable timer = startMSTimer
	if (DimSize(axes,1)!=3 || DimSize(axes,2))
		DoAlert 0, "dimensiont of axes must be axes[N][3]"
		return ""
	endif
	String fldr = GetWavesDataFolder(axes,1)
	Variable N = DimSize(axes,0)

	String wName = UniqueName("axisTemp",1,0)
	Make/N=(N) $wName
	Wave axis=$wName
	Variable Xmin, Xmax, Ymin,Ymax, Zmin, Zmax
	axis = axes[p][0]
	WaveStats/Q axis
	Xmin = V_min  ;  Xmax = V_max
	axis = axes[p][1]
	WaveStats/Q axis
	Ymin = V_min  ;  Ymax = V_max
	axis = axes[p][2]
	WaveStats/Q axis
	Zmin = V_min  ;  Zmax = V_max
	KillWaves/Z axis
//	printf "range of axisX is [%g, %g],  Æ=%g\r",Xmin,Xmax,Xmax-Xmin
//	printf "range of axisY is [%g, %g],  Æ=%g\r",Ymin,Ymax,Ymax-Ymin
//	printf "range of axisZ is [%g, %g],  Æ=%g\r",Zmin,Zmax,Zmax-Zmin
	Variable amin = min(min(Xmin,Ymin),Zmin)
	Variable amax = max(max(Xmax,Ymax),Zmax)
	Variable Nx,Ny,Nz, NN=20
	NN = 2*round(N^0.333333)
	Variable delta = (amax-amin)/(NN-1)
	Nx = ceil((Xmax-Xmin)/delta) + 1
	Ny = ceil((Ymax-Ymin)/delta) + 1
	Nz = ceil((Zmax-Zmin)/delta) + 1
	if (ItemsInList(GetRTStackInfo(0))<2)
		printf "max range is [%g, %g]    dimensions=(%d,%d,%d)\r",amin,amax,Nx,Ny,Nz
	endif

	String HistName = fldr+"Hist3dAxes"
	Make/N=(Nx,Ny,Nz)/O $HistName
	Wave Hist3dAxes=Hist3dAxes
	Hist3dAxes = 0
	SetScale/P x Xmin,delta, "", Hist3dAxes
	SetScale/P y Ymin,delta, "", Hist3dAxes
	SetScale/P z Zmin,delta, "", Hist3dAxes
	Variable i, ix,iy,iz
//printf "		up to loop in Histogram3dAxisSpace() took %g msec\r",stopMSTimer(timer)/1000
//timer = startMSTimer
	for (i=0;i<N;i+=1)
		if (numtype(axes[i][0]+axes[i][1]+axes[i][2]))
			continue					// skip invalid numbers
		endif
		ix = round( (axes[i][0]-DimOffset(Hist3dAxes,0))/DimDelta(Hist3dAxes,0) )
		iy = round( (axes[i][1]-DimOffset(Hist3dAxes,1))/DimDelta(Hist3dAxes,1) )
		iz = round( (axes[i][2]-DimOffset(Hist3dAxes,2))/DimDelta(Hist3dAxes,2) )
		Hist3dAxes[ix][iy][iz] += 1
	endfor
//printf "		loop in Histogram3dAxisSpace() took %g msec\r",stopMSTimer(timer)/1000
//timer = startMSTimer
	Hist3dAxes = Hist3dAxes[p][q][r]==0 ? NaN : Hist3dAxes[p][q][r]
//printf "		to set zeros to NaN in Histogram3dAxisSpace() took %g msec\r",stopMSTimer(timer)/1000

//timer = startMSTimer
	WaveStats/Q/M=1 Hist3dAxes
	ix = round( (V_maxRowLoc-DimOffset(Hist3dAxes,0))/DimDelta(Hist3dAxes,0) )
	iy = round( (V_maxColLoc-DimOffset(Hist3dAxes,1))/DimDelta(Hist3dAxes,1) )
	iz = round( (V_maxLayerLoc-DimOffset(Hist3dAxes,2))/DimDelta(Hist3dAxes,2) )
	if (ItemsInList(GetRTStackInfo(0))<2)
		printf "Hist3dAxes(%g)(%g)(%g) = Hist3dAxes[%d][%d][%d]=%d,     V_max=%g\r",V_maxRowLoc,V_maxColLoc,V_maxLayerLoc,ix,iy,iz,Hist3dAxes[ix][iy][iz],V_max
	endif
	String str="HIST="+GetWavesDataFolder(Hist3dAxes,2)
	str = ReplaceNumberByKey("ix", str,ix,"=")
	str = ReplaceNumberByKey("iy", str,iy,"=")
	str = ReplaceNumberByKey("iz", str,iz,"=")
	str = ReplaceNumberByKey("x", str,V_maxRowLoc,"=")
	str = ReplaceNumberByKey("y", str,V_maxColLoc,"=")
	str = ReplaceNumberByKey("z", str,V_maxLayerLoc,"=")
	peak.ix = ix
	peak.iy = iy
	peak.iz = iz
	peak.x = V_maxRowLoc
	peak.y = V_maxColLoc
	peak.z = V_maxLayerLoc
//printf "		last part of Histogram3dAxisSpace() took %g msec\r",stopMSTimer(timer)/1000
	return str
End




Function vectorInBox(box,xx,yy,zz)
	STRUCT boxStruct &box
	Variable xx,yy,zz
	if (xx<box.xlo || xx>box.xhi)
		return 0
	endif
	if (yy<box.ylo || yy>box.yhi)
		return 0
	endif
	if (zz<box.zlo || zz>box.zhi)
		return 0
	endif
	return 1
End

Static Structure peak3D
	Variable ix,iy,iz
	Variable x,y,z
EndStructure 

Static Structure boxStruct
	Variable xlo,xhi,ylo,yhi,zlo,zhi
EndStructure

