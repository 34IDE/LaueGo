#pragma rtGlobals=3		// Use modern globala access method and strict wave access.
#pragma ModuleName=IndexingInternal
#pragma version = 0.07
#include "IndexingN", version>=4.70

Static Constant hc = 1.239841857			// keV-nm
Constant threshDivide = 3

//	#define ZONE_TESTING
//	#define QS_TESTING
#if defined(ZONE_TESTING) || defined(QS_TESTING)
Menu "Zones"
	SubMenu "Test"
		"Make Simulated Test Pattern...", MakeSimulatedTestPattern(NaN)
		"Make Test Spots...", IndexingInternal#MakeTestSpots(NaN)
		 MenuItemIfWaveClassExists("  Graph Simulated Pattern...","FittedPeakList*",""), GraphTestPattern($"",$"")
		SubMenu "Browse Zone Info"
			"Browse Zone Axis Definition",BrowseURL "http://reference.iucr.org/dictionary/Zone_axis"
			"Browse Zone Axis Wikipedia", BrowseURL "http://en.wikipedia.org/wiki/Zone_axis"
		End
	End
	"-"
	MenuItemIfWaveClassExists("Make Zones Wave from Ghats...","Ghat*","MINCOLS:3"), IndexingInternal#MakeZonesWave($"",NaN)
	"-"
	MenuItemIfWaveClassExists("Table of ZonesWave","Zones*","MINCOLS:4"),/Q, DisplayTableOfWave($"",classes="Zones*",promptStr="Zones",options="MINCOLS:4",top=52,left=15,colWid=80)
	MenuItemIfWaveClassExists("Table of Ghats","QsMeasured*","MINCOLS:3"),/Q, DisplayTableOfWave($"",classes="QsMeasured*",promptStr="Ghat",options="MINCOLS:3",top=65,left=17,colWid=80)
End
#endif

//  ====================================================================================  //
//  =============================== Start of Qs Indexing ===============================  //

Function/WAVE runIndexingQsIgor(args)
	String args
	Wave FullPeakList = $StringByKey("FullPeakList",args)
	Variable keVmaxCalc = NumberByKey("keVmaxCalc",args)	// 14, maximum energy for peaks to check (keV)
	Variable keVmaxTest = NumberByKey("keVmaxTest",args)	// 30, maximum energy for finding peaks (keV)
	Variable angTol = NumberByKey("angleTolerance",args)	// ~0.5¡, angular matrch in pair angles (degree), as this gets bigger it runs slower
	Variable hp = NumberByKey("hp",args)							// preferred hkl
	Variable kp = NumberByKey("kp",args)	
	Variable lp = NumberByKey("lp",args)	
	Variable cone = NumberByKey("cone",args)					// for possible Qhats, range of allowed tilts from central hkl (degree)
	Variable maxSpots = NumberByKey("maxSpots",args)			// -n max num. of spots from FullPeakList to use, default is 250
	Wave FullPeakList1 = $StringByKey("FullPeakList1",args)
	Wave FullPeakList2 = $StringByKey("FullPeakList2",args)
	Variable printIt = NumberByKey("printIt",args)
	maxSpots = numtype(maxSpots) || maxSpots<3 ? 250 : maxSpots
	printIt = numtype(printIt) ? 0 : printIt
	if (cone<=0 || numtype(cone))
		return $""
	elseif (angTol<=0 || numtype(angTol) || angTol>10)
		return $""
	elseif (numtype(hp+kp+lp))
		return $""
	elseif (keVmaxCalc<=0 || numtype(keVmaxCalc))
		return $""
	elseif (keVmaxTest<=0 || numtype(keVmaxTest))
		return $""
	endif
	Variable sec0=stopMSTimer(-2)*1e-6, sec1=sec0

	Make/N=3/D/FREE hklPrefer={hp,kp,lp}							// hkl near center of data
	Variable NmaxHKL = 200												// max number of Possible hkl to create

	Wave GhatsMeasured=FullPeakList2Ghats(maxSpots,FullPeakList,FullPeakList1=FullPeakList1,FullPeakList2=FullPeakList2)	// convert peaks to a Qlist+intens+dNum
	String noteMeasured=note(GhatsMeasured)
	Variable NG0=DimSize(GhatsMeasured,0)
	Make/N=(NG0,3)/FREE GhatsOnly=GhatsMeasured[p][q]
	MatrixOP/FREE qvec0 = sumCols(GhatsOnly)					// assume that average Ghat is Q-vector to center of detector

#ifdef QS_TESTING
	if (printIt)
		printf "*making GhatsMeasured[%d], elapsed time = %.3f s\r",NG0,stopMSTimer(-2)*1e-6 - sec1 ; sec1 = stopMSTimer(-2)*1e-6
	endif
#endif

	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))
		DoAlert 0, "no lattice structure found, did you forget to set it?"
		return $""
	endif
	Wave direct=directFrom_xtal(xtal)
	Wave recipBase=recipFrom_xtal(xtal)
	Make/N=3/D/FREE ki={0,0,1}			// use this for the incident beam direction
#ifdef QS_TESTING
	Wave PossibleQhats = MakePossibleQhats(recipBase,qvec0,cone*PI/180,hklPrefer,keVmaxCalc,NmaxHKL,kin=ki,printIt=printIt)	// q^ available for checking
	if (printIt)
		printf "*finished setup, %d Possible Qhats,  elapsed time = %.3f s\r",DimSize(PossibleQhats,0),stopMSTimer(-2)*1e-6 - sec1 ; sec1 = stopMSTimer(-2)*1e-6
	endif
#else
	Wave PossibleQhats = MakePossibleQhats(recipBase,qvec0,cone*PI/180,hklPrefer,keVmaxCalc,NmaxHKL,kin=ki)	// q^ available for checking
#endif

	// Possible HKL Pairs & Measured Ghat Pairs are lists of {dot,i,j}
	Wave PossibleHKLPairs = MakePairsList(PossibleQhats)	// all pairs of hkl, {dot, i, j}
	Wave MeasuredGhatPairs = MakePairsList(GhatsMeasured)	//   only uses first 3 columns of argument
#ifdef QS_TESTING
	if (printIt)
		printf "*made both pairs waves, elapsed time = %.3f s\r",stopMSTimer(-2)*1e-6 - sec1 ; sec1 = stopMSTimer(-2)*1e-6
	endif
#endif

	// get all of the rotations between pairs, make a list of all the rotation vectors (radian)
	Variable dotZero = 1-cos(angTol*PI/180)						// used to compare two dot products for equivalence
	Variable Nalloc = 1000
	Make/N=(Nalloc,4)/U/I/FREE indexPairs=0
	Variable iMeasured, nMeasuredPairs=DimSize(MeasuredGhatPairs,0)
	Variable iPossible, nPossiblePairs=DimSize(PossibleHKLPairs,0)
#ifdef QS_TESTING
	if (printIt)
		printf "*there are %d pairs of Measured Ghats,  and %d pairs of Possible hkl's\r",nMeasuredPairs,nPossiblePairs
	endif
#endif
	Variable doti, nPairs
	for (iMeasured=0,nPairs=0; iMeasured<nMeasuredPairs; iMeasured+=1)
		doti = MeasuredGhatPairs[iMeasured][0]
		for (iPossible=0; iPossible<nPossiblePairs; iPossible+=1)
			if (abs(doti-PossibleHKLPairs[iPossible][0])<dotZero)			// do the pairs have the same angular separation?
				if (nPairs>=Nalloc)
					Nalloc += 1000
					Redimension/N=(Nalloc,-1) indexPairs
				endif
				indexPairs[nPairs][0] = MeasuredGhatPairs[iMeasured][1]	// save this set of pairs
				indexPairs[nPairs][1] = MeasuredGhatPairs[iMeasured][2]
				indexPairs[nPairs][2] = PossibleHKLPairs[iPossible][1]
				indexPairs[nPairs][3] = PossibleHKLPairs[iPossible][2]
				nPairs += 1
				if (mod(nPairs,10000)==0)
					print "nPairs = ",nPairs
				endif
			endif
		endfor
	endfor
	Redimension/N=(nPairs,-1) indexPairs
#ifdef QS_TESTING
	if (printIt)
		printf "*made pairs of pairs: indexPairs[%d][4], elapsed time = %.3f s\r",nPairs, stopMSTimer(-2)*1e-6 - sec1 ; sec1 = stopMSTimer(-2)*1e-6
	endif
//	Duplicate/O indexPairs, indexPairsView
#endif
	if (nPairs<1)
		printf "ERROR -- nPairs = %g,  must be ³ 2\r",nPairs
		return $""
	endif
	// each row of indexPairs generates two rotations (+ and -), so there are 32*nPairs of rotations to check
	// what is the best way to check for lumps or groupings in the set of rotation vectors?

	// generate the rotation vectors, store them in PairRotations
	Wave PairRotations = ConstructAllPairRotations(indexPairs,GhatsMeasured,PossibleQhats)
#ifdef QS_TESTING
	if (printIt)
		printf "*made PairRotations[%d][3] the rotation vectors, elapsed time = %.3f s\r",DimSize(PairRotations,0),stopMSTimer(-2)*1e-6 - sec1 ; sec1 = stopMSTimer(-2)*1e-6
	endif
	if (DataFolderExists(":GizmoWaves"))
		Duplicate/O PairRotations, :GizmoWaves:PairRotationsView
	endif
#endif

	Variable threshStart=angTol*PI/180		// largest acceptable error between measured and calculated angles
	Variable threshLimit=0.001*PI/180		// at 0.001¡ stop searching for a better match
#ifdef QS_TESTING
	if (printIt)
		printf "running with a starting threshold of %g¡, and a stopping threshold of %g¡\r",threshStart*180/PI, threshLimit*180/PI
	endif
#endif
	Make/N=3/D/FREE rotAxis, center=0		// start center at {0,0,0}
	Variable hits, hitsLast=Inf, thresh, radius=Inf	// hits is largest number of PairRotations that are the same rotation
	for (thresh=threshStart; thresh>threshLimit; thresh /= threshDivide)
		hits = FindPeakInRots(PairRotations,rotAxis,thresh,center,radius)
		NVAR J_IndexOfMax=J_IndexOfMax
#ifdef QS_TESTING
		if (printIt)
			printf "  hits = %g,   hitsLast = %g,  thresh=%.3g¡,  radius=%.3g¡,   axis=%s  (%d)\r",hits,hitsLast,thresh*180/PI,radius*180/PI,vec2str(rotAxis,zeroThresh=1e-9),J_IndexOfMax
		endif
#endif


		if ((hits+2) > hitsLast)				// no significant reduction
#ifdef QS_TESTING
			if (printIt)
				print "  stop looping because (hits+2) > hitsLast"
			endif
#endif
			break
		endif


		hitsLast = hits
		center = rotAxis							// for next iteration, reset
		radius = 1.1 * thresh					// reduce radius about center to the last threshold
	endfor
	Variable angle = norm(rotAxis)*180/PI
	Make/N=(3,3)/D/FREE rotMat
	rotationMatAboutAxis(rotAxis,angle,rotMat)	// make rotMat the matrix version of rotAxis
	MatrixOP/O/FREE recip = rotMat x recipBase
	//	printWave(recip,name="rotated recip",brief=1,zeroThresh=1e-9)
#ifdef QS_TESTING
	if (printIt)
		printf "*found <rot> = %s  |<rot>|=%.3f¡,  with %d hits, elapsed time = %.3f s\r",vec2str(rotAxis,zeroThresh=1e-9),angle,hits,stopMSTimer(-2)*1e-6-sec1 ; sec1 = stopMSTimer(-2)*1e-6
	endif
#endif

	// find rms error between Ghats and the predicted hkl directions
	STRUCT microGeometry g
	if (FillGeometryStructDefault(g))		//fill the geometry structure with current values
		DoAlert 0, "no geometry structure found, did you forget to set it?"
		return $""
	endif

	Make/N=3/D/FREE vec
	Variable startx=NumberByKey("startx",noteMeasured,"="), groupx=NumberByKey("groupx",noteMeasured,"=")
	Variable starty=NumberByKey("starty",noteMeasured,"="), groupy=NumberByKey("groupy",noteMeasured,"=")
	Variable/C pz
	Variable i, NG1, delta
	Make/N=(NG0,3)/D/FREE hkls=NaN, hklsAll=NaN
	for (i=0,NG1=0; i<NG0; i+=1)				// find the indexed hkl's from this rotation
		vec = GhatsMeasured[i][p]
		Wave hkl = LowestAllowedHKLfromGhat(xtal,recip,vec,keVmaxTest)	// reduce to nearest hkl
		if (!WaveExists(hkl))
			continue
		endif
		MatrixOP/FREE qvec = recip x hkl
		MatrixOP/FREE qhat = Normalize( qvec )
		delta = acos(limit(MatrixDot(qhat,vec),-1,1))*180/PI
		hklsAll[i][] = hkl[q]					// save found hkl for every measured G, used to fill IndexedWave[][]
		if (delta < angTol)
			hkls[NG1][] = hkl[q]				// just save the good ones, this is used for the RefineOrientation()
			NG1 += 1
		endif
	endfor
	Redimension/N=(NG1,-1) hkls

	// do a least-squares optimization on the rotation vector
	Variable err = RefineOrientation(GhatsMeasured,hkls,recipBase,rotAxis)
	angle = norm(rotAxis)*180/PI
	rotationMatAboutAxis(rotAxis,angle,rotMat)		// re-make rotMat
	MatrixOP/O/FREE recip = rotMat x recipBase		//  and re-make recip
#ifdef QS_TESTING
	if (printIt)
		printf "*after non-linear least-squares optimization, <rot> = %s  |<rot>|=%.3f¡,  elapsed time = %.3f s\r",vec2str(rotAxis,zeroThresh=1e-9),angle,stopMSTimer(-2)*1e-6-sec1 ; sec1 = stopMSTimer(-2)*1e-6
	endif
#endif

	// now save all indexed values using the optimized rotation
	Variable Npatterns=1, ipat=0, NG
	String peakListName = ParseFilePath(3,StringByKey("peakListWave",noteMeasured,"="),":",0,0)
	String IndexedName = CleanupName("FullPeakIndexed"+ReplaceString("FullPeakList",peakListName,""),0)
	Make/N=(NG0,12,Npatterns)/O $IndexedName/WAVE=IndexedWave = NaN
	Variable dNum, px,py, keV, intensity, intensityMax=-Inf
	Variable rmsGhat, minError=Inf, maxError=-Inf, goodness0=0
	for (i=0,NG=0; i<NG0; i+=1)							// re-calc with optimized rotation
		vec = GhatsMeasured[i][p]
		dNum = GhatsMeasured[i][4]
		hkl = hklsAll[i][p]
		MatrixOP/FREE qvec = recip x hkl
		MatrixOP/FREE qhat = Normalize( qvec )
		keV = -hc * norm(qvec) / ( 4*PI*MatrixDot(ki,vec) )		// |Q| = 4*¹*sin(theta)/lambda,  MatrixDot(ki,vec) =-sin(theta)
		delta = acos(limit(MatrixDot(qhat,vec),-1,1))*180/PI	// angle error (degree)
		if (delta < angTol)
			maxError = max(delta,maxError)
			minError = min(delta,minError)
			rmsGhat += delta*delta
			intensity = genericIntensity(xtal,qvec,hkl,keV)
			intensityMax = max(intensityMax,intensity)
			goodness0 += intensity
			pz = q2pixel(g.d[dNum],qhat)
			px = limit(real(pz),0,2047)
			py = limit(imag(pz),0,2047)
			px = ( px - (startx-FIRST_PIXEL) - (groupx-1)/2 )/groupx	// change to binned pixels
			py = ( py - (starty-FIRST_PIXEL) - (groupy-1)/2 )/groupy	// pixels are still zero based
			IndexedWave[NG][0,2][ipat] = qhat[q]		// Q^
			IndexedWave[NG][3,5][ipat] = hkl[q-3]		// save the hkl, of this indexed point
			IndexedWave[NG][6][ipat]   = intensity	// intensity
			IndexedWave[NG][7][ipat]   = keV			// energy (keV)
			IndexedWave[NG][8][ipat]   = delta			// angle error (degree)
			IndexedWave[NG][9][ipat]   = px				// pixel X
			IndexedWave[NG][10][ipat]  = py				// pixel Y
			IndexedWave[NG][11][ipat]  = dNum			// detector Number
			NG += 1
		endif
	endfor
	Redimension/N=(NG,-1,-1) IndexedWave
	IndexedWave[][6][] /= intensityMax					// normalize so max intensity is 1
	goodness0 *= NG*NG/intensityMax
	rmsGhat = sqrt(rmsGhat/NG)
	IndexedWave[][8][ipat] = abs(IndexedWave[p][8][ipat])<5e-6 ? 0 : IndexedWave[p][8][ipat]	// precision limit
	Variable executionTime = stopMSTimer(-2)*1e-6 - sec0
#ifdef QS_TESTING
	printf "Ghat rms error (from %d of %d) Ghats = %.3g¡,  range=[%.4f, %.4f¡],  goodness0=%g\r",NG,NG0,rmsGhat,minError,maxError,goodness0
	printf "Total time = %.3f s\r",executionTime
#endif

	SetDimLabel 1,0,Qx,IndexedWave			;	SetDimLabel 1,1,Qy,IndexedWave
	SetDimLabel 1,2,Qz,IndexedWave			;	SetDimLabel 1,3,h,IndexedWave
	SetDimLabel 1,4,k,IndexedWave			;	SetDimLabel 1,5,l,IndexedWave
	SetDimLabel 1,6,Intensity,IndexedWave	;	SetDimLabel 1,7,keV,IndexedWave
	SetDimLabel 1,8,angleErr,IndexedWave	;	SetDimLabel 1,9,pixelX,IndexedWave
	SetDimLabel 1,10,pixelY,IndexedWave	;	SetDimLabel 1,11,detNum,IndexedWave

	String str, wnote=ReplaceStringByKey("waveClass",noteMeasured,"IndexedPeakList","=")
	wnote = ReplaceStringByKey("peakListWave",wnote,GetWavesDataFolder(FullPeakList,2),"=")

	wnote = ReplaceStringByKey("structureDesc",wnote,xtal.desc,"=")
	sprintf str,"{ %g, %g, %g, %g, %g, %g }",xtal.a,xtal.b,xtal.c,xtal.alpha,xtal.beta,xtal.gam
	wnote = ReplaceStringByKey("latticeParameters",wnote,str,"=")
	wnote = ReplaceStringByKey("lengthUnit",wnote,"nm","=")
	wnote = ReplaceNumberByKey("SpaceGroup",wnote,xtal.SpaceGroup,"=")

	wnote = ReplaceNumberByKey("keVmaxCalc",wnote,keVmaxCalc,"=")
	wnote = ReplaceNumberByKey("keVmaxTest",wnote,keVmaxTest,"=")
	wnote = ReplaceStringByKey("hklPrefer",wnote,vec2str(hklPrefer),"=")
	wnote = ReplaceNumberByKey("cone",wnote,cone,"=")
	wnote = ReplaceNumberByKey("angleTolerance",wnote,angTol,"=")
	wnote = ReplaceNumberByKey("Nindexed",wnote,NG,"=")
	wnote = ReplaceNumberByKey("NiData",wnote,NG0,"=")
	wnote = ReplaceNumberByKey("executionTime",wnote,executionTime,"=")
	wnote = ReplaceNumberByKey("NpatternsFound",wnote,Npatterns,"=")
	wnote = ReplaceNumberByKey("rms_error0",wnote,rmsGhat,"=")
	wnote = ReplaceStringByKey("rotation_matrix0",wnote,encodeMatAsStr(rotMat,places=8,vsep=""),"=")
	wnote = ReplaceStringByKey("recip_lattice0",wnote,encodeMatAsStr(recip,places=8,vsep=""),"=")
	wnote = ReplaceNumberByKey("goodness0",wnote,goodness0,"=")
	Wave euler = rotAxis2EulerAngles(rotAxis)	// returns the three Euler angles (radian)
	euler *= 180/PI
	wnote = ReplaceStringByKey("EulerAngles0",wnote,vec2str(euler,places=11),"=")
	Note/K IndexedWave, wnote
	return IndexedWave
End
//
Static Function/WAVE MakePossibleQhats(recip0,qvec0,cone,hkl0,keVmax,Nmax,[kin,printIt])
	// Make a list of q^[][3] that are available for checking (UN-Rotated)
	Wave recip0						// UN-rotated reciprocal lattice
	Wave qvec0						// this q-vector diffracts to center of the detector (only use direction)
	Variable cone					// only consider q-vectors with kf within cone of qvec0+ki (radian)
	Wave hkl0						// preferred hkl, the hkl that we want to diffract near qvec0
	Variable keVmax				// maximum allowed energy (keV)
	Variable Nmax					// max number of test q-vectors
	Wave kin							// optional incident wave direction (usually defaults to 001)
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) || strlen(GetRTStackInfo(2))==0 ? 0 : !(!printIt)

	Variable tol = 0.1								// angle tolerance (degree)
	if (numtype(cone+keVmax+Nmax))
		return $""
	elseif (cone<tol || cone>(PI-tol))			// cone in range [tol, 180¡-tol]
		return $""
	elseif (keVmax<=0)
		return $""
	elseif (Nmax<3 || Nmax>250)					// max number of hkls is in [3,250]
		return $""
	endif

	Make/N=3/D/FREE ki={0,0,1}					// use this for the incident beam direction
	if (!ParamIsDefault(kin))						// one was passed, use it
		ki = kin[p]
		normalize(ki)
	endif
	MatrixOP/FREE qhat0 = Normalize(qvec0)	// direction of q-vector that will scatter toward detector center
	Redimension/N=3 qhat0

	// find rotation that takes (recip0 x hkl0) --> qvec0
	//		rot0 x (recip0 x hkl0) == qhat0
	MatrixOP/FREE q0 = Normalize(recip0 x hkl0)
	Cross q0, qhat0									// find a rotation that puts (recip0 x hkl0) parallel to qhat0
	Wave W_Cross=W_Cross
	Variable sine = normalize(W_Cross)
	Variable cosine = MatrixDot(q0,qhat0)
	Variable angle = atan2(sine,cosine)*180/PI
	Make/N=(3,3)/D/FREE rot0
	rotationMatAboutAxis(W_Cross,angle,rot0)
	KillWaves/Z W_Cross
	MatrixOP/FREE recip = rot0 x recip0
	//		so now:  qhat0 || (rot0 x recip0 x hkl0)

	Variable dot = MatrixDot(ki,qhat0)
	Make/N=3/D/FREE kf0 = (ki - 2*dot*qhat0)//		kf^ = ki^ - 2*(ki^ . q^)*q^
	// I want only kf's that are within cone of kf0
	Variable maxTheta = limit(acos(MatrixDot(kf0,ki))+cone,tol,PI-tol)/2
	Variable maxQLen = 4*PI*sin(maxTheta) * keVmax/hc	// max len of a q-vector = recip0 x hkl)
	MatrixOP/FREE hklMax = floor( rowRepeat(maxQLen,3)/sqrt(sumCols(magSqr(recip)))^t )
	Variable hmax=hklMax[0], kmax=hklMax[1], lmax=hklMax[2], 	maxMax=WaveMax(hklMax)
	if (printIt)
		printf "  find Possible Q-directions using index range: h=[%+d, %+d], k=[%+d, %+d], l=[%+d, %+d], max(theta)=%g¡\r",-hmax,hmax, -kmax,kmax, -lmax,lmax,maxTheta*180/PI
	endif

	Variable minDot = cos(cone)					// minimum allowed value of dot(kf0,kf)
	Variable cosTol=cos(tol), coneMax=-Inf
	Make/N=(Nmax,3)/FREE TestQs=NaN				// contains q^x, q^y, q^z, other has h,k,l
	Make/N=3/D/FREE hkl, kf
	Variable h,k,l, m, h1,k1,l1, N, Qmag
	for (m=1,N=0; m<=maxMax && N<Nmax; m+=1)
		l1 = min(lmax,m)
		for (l=0; l<=l1 && N<Nmax; l=LatticeSym#incrementIndex(l))
			hkl[2] = l
			k1 = min(kmax,l1)
			for (k=0; k<=k1 && N<Nmax; k=LatticeSym#incrementIndex(k))
				hkl[1] = k
				h1 = min(hmax,k1)
				for (h=0; h<=h1 && N<Nmax; h=LatticeSym#incrementIndex(h))
					hkl[0] = h
					MatrixOP/FREE/O qhat = recip x hkl
					Qmag = normalize(qhat)
					dot = MatrixDot(ki,qhat)
					kf = (ki - 2*dot*qhat)				// kf^ = ki^ - 2*(ki^ . q^)*q^
					if (Qmag>maxQLen)
						continue								// length is too long, skip
					elseif (dot > 0)						// (-ki .dot. q^) < 0 is invalid
						continue								// Bragg angle must be < 90¡
					elseif (MatrixDot(kf,kf0)<minDot)
						continue								// kf is too far from kf0
					elseif (!isNewDirection(TestQs,qhat,cosTol))
						continue								// skip already existing directions too
					endif
					TestQs[N][] = qhat[q]				// a new valid & unique q-direction
					N += 1
					coneMax = max(coneMax,acos(MatrixDot(kf,kf0)))
				endfor
			endfor
		endfor
	endfor
	Redimension/N=(N,-1) TestQs
	MatrixOP/FREE/O TestQs = ( Inv(rot0) x TestQs^t )^t	// rotate TestQs back to UN-rotated frame

	SetDimLabel 1,0,Qx,TestQs							// first 3 columns are the normalized qx,qy,qz
	SetDimLabel 1,1,Qy,TestQs
	SetDimLabel 1,2,Qz,TestQs
	String wnote="waveClass=QhatsTest;"
	wnote = ReplaceStringByKey("qvec0",wnote,vec2str(qvec0,places=14,sep=","),"=")
	wnote = ReplaceNumberByKey("cone",wnote,cone*180/PI,"=")
	wnote = ReplaceStringByKey("hkl0",wnote,vec2str(hkl0,places=14,sep=","),"=")
	wnote = ReplaceNumberByKey("keVmax",wnote,keVmax,"=")
	if (Nmax!=250)
		wnote = ReplaceNumberByKey("Nmax",wnote,Nmax,"=")
	endif
	wnote = ReplaceNumberByKey("coneMax",wnote,coneMax*180/PI,"=")
	wnote = ReplaceStringByKey("recip",wnote,encodeMatAsStr(recip0,places=14),"=")
	wnote = ReplaceStringByKey("ki",wnote,vec2str(ki,places=14,sep=","),"=")
	Note/K TestQs, wnote
	return TestQs
End
//Function Test_MakePossibleQhats()
//	STRUCT crystalStructure xtal
//	if (FillCrystalStructDefault(xtal))
//		DoAlert 0, "no lattice structure found, did you forget to set it?"
//		return 1
//	endif
//	Wave recip=recipFrom_xtal(xtal)
//
//	Make/D/FREE qvec0={0,1,-1}
//	Variable cone=40
//	Make/N=3/D/FREE hkl0={0,0,2}
//	Variable keVmax=30, Nmax=250
//
//	Variable sec0=stopMSTimer(-2)*1e-6
//	Wave qhats = MakePossibleQhats(recip,qvec0,cone*PI/180,hkl0,keVmax,Nmax,printIt=1)
//	printf "Total time = %.3f s\r",stopMSTimer(-2)*1e-6 - sec0
//
//	print "coneMax = ",NumberByKey("coneMax",note(qhats),"=")
//	DuplicaTe/O qhats, qhatsView
//	Variable N=DimSize(qhats,0)
//	MatrixOP/O anglesView = acos(NormalizeRows(qhats) x Normalize(qvec0))
//	anglesView *= 180/PI
//End

//  ================================ End of Qs Indexing ================================  //
//  ====================================================================================  //




//  ====================================================================================  //
//  ============================== Start of Zones Indexing =============================  //

// for Zones
//		hu+kv+lw=0
//			[uvw] means:		direct x uvw	= a direct lattice vector
//			(hkl) means:		recip  x hkl	= a reciprocal lattice vector
//			u,v,w, are all integers and are used to describe a direct lattice vector
//			h,k,l, are all integers and are used to describe a reciprocal lattice vector
//
// see:			http://reference.iucr.org/dictionary/Zone_axis
// and also:	http://en.wikipedia.org/wiki/Zone_axis


// Solve with angles as I usualy do, but only use angles of zones, not Ghats.
//	This should let me use a reduced angle tolerance and fewer hkl's.
// Use rotation vectors rather than Euler angles.

Function/WAVE runIndexingZonesIgor(args)
	String args
	Wave FullPeakList = $StringByKey("FullPeakList",args)
	Variable keVmax = NumberByKey("keVmaxTest",args)	// 30, maximum energy for finding peaks (keV)
	Variable angTol = NumberByKey("angleTolerance",args)	// ~0.5¡, angular matrch in pair angles (degree), as this gets bigger it runs slower
	Variable hp = NumberByKey("hp",args)					// preferred hkl
	Variable kp = NumberByKey("kp",args)	
	Variable lp = NumberByKey("lp",args)	
	Variable cone = NumberByKey("cone",args)			// for possible zone axes, range of allowed tilts from central hkl (degree)
	Variable maxSpots = NumberByKey("maxSpots",args)	// -n max num. of spots from FullPeakList to use, default is 250
	Wave FullPeakList1 = $StringByKey("FullPeakList1",args)
	Wave FullPeakList2 = $StringByKey("FullPeakList2",args)
	Variable printIt = NumberByKey("printIt",args)
	maxSpots = numtype(maxSpots) || maxSpots<3 ? 250 : maxSpots
	printIt = numtype(printIt) ? 0 : printIt
	if (cone<=0 || numtype(cone))
		return $""
	elseif (angTol<=0 || numtype(angTol) || angTol>10)
		return $""
	elseif (numtype(hp+kp+lp))
		return $""
	elseif (keVmax<=0 || numtype(keVmax))
		return $""
	endif
	Variable sec0=stopMSTimer(-2)*1e-6, sec1=sec0

	Wave GhatsMeasured=FullPeakList2Ghats(maxSpots,FullPeakList,FullPeakList1=FullPeakList1,FullPeakList2=FullPeakList2)	// convert peaks to a Qlist+intens+dNum
	String noteMeasured=note(GhatsMeasured)
	Wave ZonesWave=MakeZonesWave(GhatsMeasured, 4, tolAngle=angTol)

#ifdef ZONE_TESTING
	if (printIt)
		printf "*making GhatsMeasured and ZonesWave done, elapsed time = %.3f s\r",stopMSTimer(-2)*1e-6 - sec1 ; sec1 = stopMSTimer(-2)*1e-6
	endif
#endif

	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))
		DoAlert 0, "no lattice structure found, did you forget to set it?"
		return $""
	endif
	Wave direct=directFrom_xtal(xtal)
	Wave recipBase=recipFrom_xtal(xtal)
	// perpDir is the q-hat direction that diffract to near detector center (UN-rotated lattice coords)
	Make/N=3/D/FREE hklPrefer={hp,kp,lp}
	MatrixOP/FREE perpDir = recipBase x hklPrefer

#ifdef ZONE_TESTING
	if (printIt)
		printf "*finished setup, elapsed time = %.3f s\r",stopMSTimer(-2)*1e-6 - sec1 ; sec1 = stopMSTimer(-2)*1e-6
	endif
#endif

	Variable maxLatticeLen = ( MatrixDet(direct)^(1/3) ) * 6
	Variable NmaxTestZones = round(8000/DimSize(ZonesWave,0))	// maximum number of test zones to make
	NmaxTestZones = limit(NmaxTestZones,4,1000)					// want N(testZones)*N(measuredZones) ~ 8000
#ifdef ZONE_TESTING
	Wave PossibleZoneAxes=MakeVecHatsTestWave(direct,perpDir,cone*PI/180,maxLatticeLen,NmaxTestZones,printIt=printIt)
	if (printIt)
		printf "*made %d (out of %d) PossibleZoneAxes to test against, elapsed time = %.3f s\r",DimSize(PossibleZoneAxes,0),NmaxTestZones,stopMSTimer(-2)*1e-6 - sec1 ; sec1 = stopMSTimer(-2)*1e-6
	endif
#else
	Wave PossibleZoneAxes=MakeVecHatsTestWave(direct,perpDir,cone*PI/180,maxLatticeLen,NmaxTestZones)
#endif

	// PossibleZonePairs & MeasuredZonePairs are lists of {dot,i,j}
	Wave PossibleZonePairs = MakePairsList(PossibleZoneAxes)	// all pairs of directions, returns {dot,i,j} of each pair
	Wave MeasuredZonePairs = MakePairsList(ZonesWave)			//   only uses first 3 columns of argument
#ifdef ZONE_TESTING
	if (printIt)
		printf "*made both pairs waves, elapsed time = %.3f s\r",stopMSTimer(-2)*1e-6 - sec1 ; sec1 = stopMSTimer(-2)*1e-6
	endif
#endif


	// get all of the rotations between pairs
	// make a list of all the rotation vectors (radian)


	Variable dotTol = 1-cos(angTol*PI/180)
	Variable Nalloc = 1000
	Make/N=(Nalloc,4)/U/I/FREE indexPairs=0

	Variable iMeasured, nMeasuredPairs=DimSize(MeasuredZonePairs,0)
	Variable iPossible, nPossiblePairs=DimSize(PossibleZonePairs,0)
#ifdef ZONE_TESTING
	if (printIt)
		printf "*there are %d pairs of Measured Zones,  and %d pairs of Possible zones\r",nMeasuredPairs,nPossiblePairs
	endif
#endif
	Variable doti, nPairs
	for (iMeasured=0,nPairs=0; iMeasured<nMeasuredPairs; iMeasured+=1)
		doti = MeasuredZonePairs[iMeasured][0]
		for (iPossible=0; iPossible<nPossiblePairs; iPossible+=1)
			if (abs(doti-PossibleZonePairs[iPossible][0])<dotTol)			// do the pairs have the same angular separation?
				if (nPairs>=Nalloc)
					Nalloc += 1000
					Redimension/N=(Nalloc,-1) indexPairs
				endif
				indexPairs[nPairs][0] = MeasuredZonePairs[iMeasured][1]	// save this set of pairs
				indexPairs[nPairs][1] = MeasuredZonePairs[iMeasured][2]
				indexPairs[nPairs][2] = PossibleZonePairs[iPossible][1]
				indexPairs[nPairs][3] = PossibleZonePairs[iPossible][2]
				nPairs += 1
				if (mod(nPairs,10000)==0)
					print "nPairs = ",nPairs
				endif
			endif
		endfor
	endfor
	Redimension/N=(nPairs,-1) indexPairs
#ifdef ZONE_TESTING
	if (printIt)
		printf "*made pairs of pairs: indexPairs[%d][4], elapsed time = %.3f s\r",nPairs, stopMSTimer(-2)*1e-6 - sec1 ; sec1 = stopMSTimer(-2)*1e-6
	endif
#endif
	if (nPairs<1)
		printf "ERROR -- nPairs = %g,  must be ³ 2\r",nPairs
		return $""
	endif

	// each row of indexPairs generate two rotations (+ and -), so there are 1228 rotations to check
	// what is the best way to check for lumps or groupings in the set of rotation vectors?
	//
	// first generate the rotation vectors in PairRotations

	Wave PairRotations = ConstructAllPairRotations(indexPairs,ZonesWave,PossibleZoneAxes)
#ifdef ZONE_TESTING
	if (printIt)
		printf "*made PairRotations[%d][3] the rotation vectors, elapsed time = %.3f s\r",DimSize(PairRotations,0),stopMSTimer(-2)*1e-6 - sec1 ; sec1 = stopMSTimer(-2)*1e-6
	endif
	if (DataFolderExists(":GizmoWaves"))
		Duplicate/O PairRotations, :GizmoWaves:PairRotationsView
	endif
#endif

	Variable threshStart=angTol*PI/180		// largest acceptable error between measured and calculated angles
	Variable threshLimit=0.001*PI/180		// at 0.001¡ stop searching for a better match
#ifdef ZONE_TESTING
	if (printIt)
		printf "running with a starting threshold of %g¡, and a stopping threshold of %g¡\r",threshStart*180/PI, threshLimit*180/PI
	endif
#endif
	Make/N=3/D/FREE rotAxis, center=0		// start center at {0,0,0}
	Variable hits, hitsLast=Inf, thresh, radius=Inf	// hits is largest number of PairRotations that are the same rotation
	for (thresh=threshStart; thresh>threshLimit; thresh /= threshDivide)
		hits = FindPeakInRots(PairRotations,rotAxis,thresh,center,radius)
		NVAR J_IndexOfMax=J_IndexOfMax
#ifdef ZONE_TESTING
		if (printIt)
			printf "  hits = %g,   hitsLast = %g,  thresh=%.3g¡,  radius=%.3g¡,   axis=%s  (%d)\r",hits,hitsLast,thresh*180/PI,radius*180/PI,vec2str(rotAxis,zeroThresh=1e-9),J_IndexOfMax
		endif
#endif


		if ((hits+2) > hitsLast)				// no significant reduction
#ifdef ZONE_TESTING
			if (printIt)
				print "  stop looping because (hits+2) > hitsLast"
			endif
#endif
			break
		endif


		hitsLast = hits
		center = rotAxis							// for next iteration, reset
		radius = 1.1 * thresh					// reduce radius about center to the last threshold
	endfor
	Variable angle = norm(rotAxis)*180/PI
	Make/N=(3,3)/D/FREE rotMat
	rotationMatAboutAxis(rotAxis,angle,rotMat)	// make rotMat the matrix version of rotAxis
	MatrixOP/O/FREE recip = rotMat x recipBase
#ifdef ZONE_TESTING
	if (printIt)
		printf "*found <rot> = %s  |<rot>|=%.3f¡,  with %d hits, elapsed time = %.3f s\r",vec2str(rotAxis,zeroThresh=1e-9),angle,hits,stopMSTimer(-2)*1e-6-sec1 ; sec1 = stopMSTimer(-2)*1e-6
	endif
#endif

	// find rms error in the zones directions
	// rotate PossibleZoneAxes by rotAxis for comparison to measured Zone axes
	MatrixOP/FREE PossibleZoneAxes = ( rotMat x (PossibleZoneAxes^t) )^t
#ifdef ZONE_TESTING
if (DataFolderExists(":GizmoWaves"))
	Duplicate/O PossibleZoneAxes, :GizmoWaves:testZoneHatsGizmo
endif
#endif
	Make/N=3/D/FREE vec
	Variable i, delta, rmsZones, NZ
	for (i=0,NZ=0; i<DimSize(ZonesWave,0); i+=1)
		vec = ZonesWave[i][p]
		MatrixOp/FREE dots = PossibleZoneAxes x vec
		delta = acos(limit(WaveMax(dots),-1,1))*180/PI	// angle to closest zone axis (degree)
		if (delta < angTol)
			rmsZones += delta*delta
			NZ += 1
		endif
	endfor
	rmsZones = sqrt(rmsZones/NZ)
#ifdef ZONE_TESTING
	if (printIt)
		printf "Zones rms error (from %d of %d) zones = %.3g¡\r",NZ,DimSize(ZonesWave,0),rmsZones
	endif
#endif

	// find rms error of the Ghats to predicted hkl directions
	STRUCT microGeometry g
	if (FillGeometryStructDefault(g))	//fill the geometry structure with current values
		DoAlert 0, "no geometry structure found, did you forget to set it?"
		return $""
	endif

	Variable startx=NumberByKey("startx",noteMeasured,"="), groupx=NumberByKey("groupx",noteMeasured,"=")
	Variable starty=NumberByKey("starty",noteMeasured,"="), groupy=NumberByKey("groupy",noteMeasured,"=")

	Variable NG1, NG0=DimSize(GhatsMeasured,0)
	Variable/C pz
	Make/N=3/D/FREE ki={0,0,1}			// use this for the incident beam direction
	Make/N=(NG0,3)/D/FREE hkls=NaN, hklsAll=NaN
	for (i=0,NG1=0; i<NG0; i+=1)
		vec = GhatsMeasured[i][p]
		Wave hkl = LowestAllowedHKLfromGhat(xtal,recip,vec,keVmax)	// reduce to nearest hkl
		if (!WaveExists(hkl))
			continue
		endif
		MatrixOP/FREE qvec = recip x hkl
		MatrixOP/FREE qhat = Normalize( qvec )
		delta = acos(limit(MatrixDot(qhat,vec),-1,1))*180/PI
		hklsAll[i][] = hkl[q]				// save found hkl for every measured G, used to fill IndexedWave[][]
		if (delta < angTol)
			hkls[NG1][] = hkl[q]			// just save the good ones, this is used for the RefineOrientation()
			NG1 += 1
		endif
	endfor
	Redimension/N=(NG1,-1) hkls

	// do a least-squares optimization on the rotation vector
	Variable err = RefineOrientation(GhatsMeasured,hkls,recipBase,rotAxis)
	angle = norm(rotAxis)*180/PI
	rotationMatAboutAxis(rotAxis,angle,rotMat)	// re-make rotMat and recip
	MatrixOP/O/FREE recip = rotMat x recipBase	// re-make recip
#ifdef ZONE_TESTING
	if (printIt)
		printf "*after non-linear least-squares optimization, <rot> = %s  |<rot>|=%.3f¡,  elapsed time = %.3f s\r",vec2str(rotAxis,zeroThresh=1e-9),angle,stopMSTimer(-2)*1e-6-sec1 ; sec1 = stopMSTimer(-2)*1e-6
	endif
#endif

	// now save all indexed values using the optimized rotation
	Variable Npatterns=1
	String peakListName = ParseFilePath(3,StringByKey("peakListWave",noteMeasured,"="),":",0,0)
	String IndexedName = CleanupName("FullPeakIndexed"+ReplaceString("FullPeakList",peakListName,""),0)
	Make/N=(NG0,12,Npatterns)/O $IndexedName/WAVE=IndexedWave = NaN
	Variable dNum, px,py, keV, intensity, intensityMax=-Inf, ipat=0
	Variable NG, rmsGhat, minError=Inf, maxError=-Inf, goodness0=0
	for (i=0,NG=0; i<NG0; i+=1)					// re-calc with optimized rotation
		vec = GhatsMeasured[i][p]
		dNum = GhatsMeasured[i][4]
		hkl = hklsAll[i][p]
		MatrixOP/FREE qvec = recip x hkl
		MatrixOP/FREE qhat = Normalize( qvec )
		keV = -hc * norm(qvec) / ( 4*PI*MatrixDot(ki,vec) )	// |Q| = 4*¹*sin(theta)/lambda
		delta = acos(limit(MatrixDot(qhat,vec),-1,1))*180/PI
		if (delta < angTol)
			maxError = max(delta,maxError)
			minError = min(delta,minError)
			rmsGhat += delta*delta
			intensity = genericIntensity(xtal,qvec,hkl,keV)
			intensityMax = max(intensityMax,intensity)
			goodness0 += intensity
			pz = q2pixel(g.d[dNum],qhat)
			px = limit(real(pz),0,2047)
			py = limit(imag(pz),0,2047)
			px = ( px - (startx-FIRST_PIXEL) - (groupx-1)/2 )/groupx	// change to binned pixels
			py = ( py - (starty-FIRST_PIXEL) - (groupy-1)/2 )/groupy	// pixels are still zero based
			IndexedWave[NG][0,2][ipat] = qhat[q]		// Q^
			IndexedWave[NG][3,5][ipat] = hkl[q-3]		// save the hkl, of this indexed point
			IndexedWave[NG][6][ipat]   = intensity	// intensity
			IndexedWave[NG][7][ipat]   = keV			// energy (keV)
			IndexedWave[NG][8][ipat]   = delta			// angle error (degree)
			IndexedWave[NG][9][ipat]   = px				// pixel X
			IndexedWave[NG][10][ipat]  = py				// pixel Y
			IndexedWave[NG][11][ipat]  = dNum			// detector Number
			NG += 1
		endif
	endfor
	Redimension/N=(NG,-1,-1) IndexedWave
	IndexedWave[][6][] /= intensityMax
	goodness0 *= NG*NG/intensityMax
	rmsGhat = sqrt(rmsGhat/NG)
	Variable executionTime = stopMSTimer(-2)*1e-6 - sec0
#ifdef ZONE_TESTING
	printf "Ghat rms error (from %d of %d) Ghats = %.3g¡,  range=[%.4f, %.4f¡],  goodness0=%g\r",NG,NG0,rmsGhat,minError,maxError,goodness0
	printf "Total time = %.3f s\r",executionTime
#endif

	SetDimLabel 1,0,Qx,IndexedWave			;	SetDimLabel 1,1,Qy,IndexedWave
	SetDimLabel 1,2,Qz,IndexedWave			;	SetDimLabel 1,3,h,IndexedWave
	SetDimLabel 1,4,k,IndexedWave			;	SetDimLabel 1,5,l,IndexedWave
	SetDimLabel 1,6,Intensity,IndexedWave	;	SetDimLabel 1,7,keV,IndexedWave
	SetDimLabel 1,8,angleErr,IndexedWave	;	SetDimLabel 1,9,pixelX,IndexedWave
	SetDimLabel 1,10,pixelY,IndexedWave	;	SetDimLabel 1,11,detNum,IndexedWave

	String str, wnote=ReplaceStringByKey("waveClass",noteMeasured,"IndexedPeakList","=")
	wnote = ReplaceStringByKey("peakListWave",wnote,GetWavesDataFolder(GhatsMeasured,2),"=")

	wnote = ReplaceStringByKey("structureDesc",wnote,xtal.desc,"=")
	sprintf str,"{ %g, %g, %g, %g, %g, %g }",xtal.a,xtal.b,xtal.c,xtal.alpha,xtal.beta,xtal.gam
	wnote = ReplaceStringByKey("latticeParameters",wnote,str,"=")
	wnote = ReplaceStringByKey("lengthUnit",wnote,"nm","=")
	wnote = ReplaceNumberByKey("SpaceGroup",wnote,xtal.SpaceGroup,"=")

	wnote = ReplaceNumberByKey("keVmaxCalc",wnote,keVmax,"=")
	wnote = ReplaceNumberByKey("keVmaxTest",wnote,keVmax,"=")
	wnote = ReplaceStringByKey("hklPrefer",wnote,vec2str(hklPrefer),"=")
	wnote = ReplaceNumberByKey("cone",wnote,cone,"=")
	wnote = ReplaceNumberByKey("angleTolerance",wnote,angTol,"=")
	wnote = ReplaceNumberByKey("Nindexed",wnote,NG,"=")
	wnote = ReplaceNumberByKey("NiData",wnote,NG0,"=")
	wnote = ReplaceNumberByKey("executionTime",wnote,executionTime,"=")
	wnote = ReplaceNumberByKey("NpatternsFound",wnote,Npatterns,"=")
	wnote = ReplaceNumberByKey("rms_error0",wnote,rmsGhat,"=")
	wnote = ReplaceStringByKey("rotation_matrix0",wnote,encodeMatAsStr(rotMat,places=8,vsep=""),"=")
	wnote = ReplaceStringByKey("recip_lattice0",wnote,encodeMatAsStr(recip,places=8,vsep=""),"=")
	wnote = ReplaceNumberByKey("goodness0",wnote,goodness0,"=")
	Wave euler = rotAxis2EulerAngles(rotAxis)	// returns the three Euler angles (radian)
	euler *= 180/PI
	wnote = ReplaceStringByKey("EulerAngles0",wnote,vec2str(euler,places=11),"=")

	// still to add to wnote:
	//		AtomDesctiption1={Si1  0 0 0 1}

	Note/K IndexedWave, wnote
	return IndexedWave
End
//
Static Function/WAVE MakeVecHatsTestWave(lattice,perpDirIN,cone,maxLatticeLen,Nmax,[tol,printIt])
	// Make a list of unit vectors from lattice, contains zone axes to consider
	Wave lattice					// either recip or direct
	Wave perpDirIN					// points toward q-vector that diffracts to detector center (in UN-rotated lattice coords)
	Variable cone					// only consider vectors at least 90-cone from perpDir (radian)
	Variable maxLatticeLen		// max len of a (lattice x ijk) vector
	Variable Nmax					// max number of test directions
	Variable tol					// angle (degree), two directions are the same if they are this close, usually 0.1¡
	Variable printIt
	tol = ParamIsDefault(tol) || numtype(tol) || tol<=0 ? 0.1*PI/180 : tol*PI/180
	printIt = ParamIsDefault(printIt) || numtype(printIt) || strlen(GetRTStackInfo(2))==0 ? 0 : !(!printIt)
	cone = limit(cone,tol,PI/2-tol)				// cone in range [tol, 90-tol]

	Make/N=3/D/FREE perpDir = perpDirIN
	normalize(perpDir)

	Variable imax,jmax,kmax, maxMax, i1,j1,k1
	Make/N=3/D/FREE vec, ijk
	vec = lattice[p][2]
	kmax = floor(maxLatticeLen/norm(vec))
	vec = lattice[p][1]
	jmax = floor(maxLatticeLen/norm(vec))
	vec = lattice[p][0]
	imax = floor(maxLatticeLen/norm(vec))
	maxMax = max(max(imax,jmax),kmax)
	if (printIt)
		printf "  find directions of matrix using indicies in: i=[%d, %d], j=[%d, %d], k=[%d, %d]\r",-imax,imax, -jmax,jmax, -kmax,kmax
	endif

	Variable maxDot = 1-cos(cone)				// maximum allowed value of dot(perpDir,zoneAxis)
	Variable cosTol = cos(tol)					// convert angle to dot, this should be slightly less than 1.0
	Make/N=(Nmax,3)/FREE TestHats=NaN			// contains x,y,z
	Variable i,j,k,m, N
	for (m=1,N=0; m<=maxMax && N<Nmax; m+=1)
		k1 = min(kmax,m)
		for (k=0; k<=k1 && N<Nmax; k=LatticeSym#incrementIndex(k))
			ijk[2] = k
			j1 = min(jmax,k1)
			for (j=0; j<=j1 && N<Nmax; j=LatticeSym#incrementIndex(j))
				ijk[1] = j
				i1 = min(imax,j1)
				for (i=0; i<=i1 && N<Nmax; i=LatticeSym#incrementIndex(i))
					ijk[0] = i
					MatrixOP/FREE/O vec = lattice x ijk
					if (normalize(vec)>maxLatticeLen)
						continue								// length is too long, skip
					elseif (abs(MatrixDot(vec,perpDir))>maxDot)
						continue								// zone axis too close to perpDir
					elseif (!isNewDirection(TestHats,vec,cosTol))
						continue								// skip already existing directions too
					endif
					TestHats[N][] = vec[q]				// save this vec as a new test zone axis
					N += 1
				endfor
			endfor
		endfor
	endfor
	Redimension/N=(N,-1) TestHats
	SetDimLabel 1,0,x,TestHats			// first 3 columns are the zone axis (unit vector)
	SetDimLabel 1,1,y,TestHats
	SetDimLabel 1,2,z,TestHats

	String wnote="waveClass=DirectionsTest;"
	wnote = ReplaceStringByKey("perpDir",wnote,vec2str(perpDir,sep=","),"=")
	wnote = ReplaceNumberByKey("ZoneAngleRange",wnote,cone*180/PI,"=")
	Note/K TestHats, wnote
	return TestHats
End
//
Static Function/WAVE MakeZonesWave(GhatsMeasured,Nmin,[tolAngle,printIt])
	Wave GhatsMeasured	// the measured Ghats (with identified hkls, only used for testing phase)
	Variable Nmin			// minimum number of spot required to indentify a zone
	Variable tolAngle		// angular tolerance used in accepting Q-vec as part of a zone
	Variable printIt
	tolAngle = ParamIsDefault(tolAngle) || numtype(tolAngle) || tolAngle<=0 ? 1 : tolAngle
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : !(!printIt)
	Variable smallestNmin=3

	String GhatList=WaveListClass("QsMeasured*","*","MINCOLS:3")
	if (!WaveExists(GhatsMeasured) && ItemsInList(GhatList)==1)
		Wave GhatsMeasured = $StringFromList(0,GhatList)
	elseif (!WaveExists(GhatsMeasured) && ItemsInList(GhatList)<1)
		return $""
	endif

	if (!WaveExists(GhatsMeasured) || Nmin<smallestNmin || numtype(Nmin))
		if (WaveExists(GhatsMeasured))
			String name=NameOfWave(GhatsMeasured)
		endif
		Nmin = numtype(Nmin) ? 7 : limit(Nmin,smallestNmin,100)
		Prompt name,"Ghats Wave",popup,GhatList
		Prompt Nmin,"min number on a zone [3,100]"
		Prompt tolAngle, "angle tolerance in defining a zone (degree)"
		if (WaveExists(GhatsMeasured))
			DoPrompt "Find Zones", Nmin, tolAngle
		else
			DoPrompt "Find Zones",name, Nmin, tolAngle
			Wave GhatsMeasured = $name
		endif
		if (V_flag)
			return $""
		endif
		printIt = 1
	endif
	if (!WaveExists(GhatsMeasured) || Nmin<smallestNmin || numtype(Nmin))
		return $""
	endif
	Variable N=DimSize(GhatsMeasured,0)
	if (!(N>=Nmin))
		return $""
	endif

	if (printIt)
		printf "MakeZonesWave(%s, %g",NameOfWave(GhatsMeasured),Nmin
		if (!ParamIsDefault(tolAngle) || tolAngle!=1)
			printf ", tolAngle=%g",tolAngle
		endif
		if (!ParamIsDefault(printIt))
			printf ", printIt=%g",printIt
		endif
		printf ")\r"
	endif

	Wave zw = FindZonesWithNspots(GhatsMeasured,Nmin,tolAngle)
	Variable Nz = WaveExists(zw) ? DimSize(zw,0) : 0
	if (Nz<1)
		printf "for Nmin = %d,  No zones found\r",Nmin
		return $""
	endif

	zw = abs(zw)<1e-15 ? 0 : zw

	if (printIt)
		printf "for Nmin = %d,  found %d measured zones\r",Nmin,Nz
	endif
#ifdef ZONE_TESTING
//Duplicate/O zw, ZonesWave
GizmoMakeWavesForZones(GhatsMeasured,zw,$"")
#endif
	return zw
End
//
Static Function/WAVE FindZonesWithNspots(Ghats,Nmin,tolAngle)
	Wave Ghats				// list of measured G^ (first 3 columns contain G^, the rest are ignored)
	Variable Nmin			// minimum number of spot to constiture a zone
	Variable tolAngle		// angular tolerance used in accepting Q-vec as part of a zone

	Variable NG=DimSize(Ghats,0)
	if (NG<Nmin || Nmin<=2)
		return $""
	endif

	Variable Nalloc=100
	Make/N=(Nalloc)/I/FREE nZones=0		// holds number of spots contributing to a zone
	Make/N=(Nalloc,3)/D/FREE zones=0	// direction of each zone axis
																			// the next two waves are only used for testing
	Variable minDot = cos(5.0*PI/180)	// G-vectors must be > 5 degrees apart to define a zone (not too parallel)
	Variable tol = cos(tolAngle*PI/180)	// tolerance on dot product
	Variable sinMax = abs(sin(tolAngle*PI/180))
	Make/N=3/D/FREE g0, g1, axis, axisAvg

	String gIndexList							// holds list of G^'s used to form one zone
	Variable Nz=0								// Nz is number of zones found, m is first spot to start searching
	Variable iz									// number of spots in this zone
	Variable mag, dot
	Variable i,j, m=0							// m is first spot to start searching
	for (m=0;m<(NG-1);m+=1)				// m is first spot that defines this possible zone
		g0 = Ghats[m][p]
		for (j=m+1; j<(NG-Nmin+1); j+=1)// find second vector to define the zone
			g1 = Ghats[j][p]
			if (abs(MatrixDot(g0,g1)) > minDot)
				continue							// vectors too parallel for defining a zone axis, keep searching
			endif

			// have a possible pair, (g0,g1) from (m,j)
			Cross g0,g1
			Wave W_Cross=W_Cross
			axis = W_Cross						// possible new zone axis
			normalize(axis)
			if (!isNewAxis(zones,axis,tol))
				continue							// already have this axis, try again
			endif
			// the (g0,g1) pair defining axis[3], has not been found before, check it out

			// have a new unique zone axis, find the spots belonging to this zone
			axisAvg = 0
			gIndexList = ""
			for (i=0,iz=0; i<NG; i+=1)	// find how many other spots belong to this axis
				g1 = Ghats[i][p]
				Cross g0,g1
				mag = normalize(W_cross)
				dot = MatrixDot(W_cross,axis)
				if (mag<sinMax)				// only look at the cross product when g0 & g1 not too parallel
					iz += 1						// this spot real close to g0, increment iz
					gIndexList += num2istr(i)+";"
					axisAvg += axis
				elseif (abs(dot) > tol)
					iz += 1						// a spot in zone, increment iz
					gIndexList += num2istr(i)+";"
					axisAvg += dot<0 ? -W_cross : W_cross
				endif
			endfor								// for i
			if (iz >= Nmin)					// found another zone axis with at least Nmin spots, save it
				if (Nz>=Nalloc)
					Nalloc += 100				// need to extend zones, nZones
					Redimension/N=(Nalloc,-1) zones
					Redimension/N=(Nalloc) nZones
				endif
				normalize(axisAvg)			// want a unit vector, so don't bother to divide by iz
				zones[Nz][] = axisAvg[q]	// axis of this zone
				nZones[Nz] = iz				// save number of spots contributing to this zone
				Nz += 1
			endif
		endfor									// for j
	endfor										// for m
	KillWaves/Z W_Cross
	if (Nz<1)
		return $""								// failed to find anything
	endif

	Redimension/N=(Nz) nZones				// want to sort zones by decreasing number of nZones
	Duplicate/FREE nZones, indexWave
	MakeIndex/R nZones, indexWave

	// assemble the output wave
	Make/N=(Nz,4)/D/FREE zoneOut		// this is the output wave
	zoneOut[][0,2] = zones[indexWave[p]][q]
	zoneOut[][3] = nZones[indexWave[p]]
	WaveClear zones, nZones
	SetDimLabel 1,0,Zx,zoneOut			// first 3 columns are the zone axis (unit vector)
	SetDimLabel 1,1,Zy,zoneOut
	SetDimLabel 1,2,Zz,zoneOut
	SetDimLabel 1,3,Ns,zoneOut			// column 3 is number of spots used to define each zone axis

	String wnote = note(Ghats)
	wnote = ReplaceStringByKey("waveClass",wnote,"Zones","=")
	wnote = ReplaceNumberByKey("Nmin",wnote,Nmin,"=")
	wnote = ReplaceNumberByKey("tolAngle",wnote,tolAngle,"=")
	wnote = ReplaceStringByKey("GhatSource",wnote,GetWavesDataFolder(Ghats,2),"=")
	Note/K zoneOut, wnote
	return zoneOut
End
//
//	returns True if axis is not parallel (or anti-parallel) to the existing directions in Gtest
Static Function isNewAxis(existingAxes,axis,tol)	// only used by FindZonesWithNspots()
	Wave existingAxes				// array(N,3) list of existing axes
	Wave axis
	Variable tol					// tolerance on dot:  tol = cos(Æangle),  cos(0.1¡) = 0.999998476913288
	if (norm(axis)==0)
		return 0
	endif
	MatrixOP/FREE dots = abs(existingAxes x axis)
	Variable dotMax = WaveMax(dots)
	return (dotMax < tol)
End

//  =============================== End of Zones Indexing ==============================  //
//  ====================================================================================  //




//  ====================================================================================  //
//  ============================= Start of Common Indexing =============================  //

Static Function RefineOrientation(Gmeasured,hkls,recipBase,rotAxis)
	Wave Gmeasured			// list of Gmeasured for ALL measured reflections
	Wave hkls				// provides the hkls
	Wave recipBase			// unrotated recip lattice
	Wave rotAxis			// starting rotation axis, contains answer on return

	Variable i,N=DimSize(hkls,0)
	Make/N=(N,6)/D/FREE wOpt=NaN
	wOpt[][3,5] = hkls[p][q-3]

	Wave rotMat = RotMatAboutAxisOnly(rotAxis)
	MatrixOP/FREE recip = rotMat x recipBase// rotated recip

	Make/N=(DimSize(Gmeasured,0),3)/D/FREE GhatsAll=Gmeasured[p][q]
	Make/N=3/D/FREE qhat, hkl
	for (i=0;i<N;i+=1)								// fill Ghats with measured G^ closest to each indexed Q^
		hkl = wOpt[i][p+3]
		MatrixOP/FREE qhat = Normalize( recip x hkl )
		MatrixOP/FREE dots = GhatsAll x qhat
		WaveStats/Q/M=1 dots
		if (V_maxLoc>=0)
			wOpt[i][0,2] = GhatsAll[V_maxLoc][q]
		endif
	endfor

	Note/K wOpt, encodeMatAsStr(recipBase)
	Variable err = RefineOrientation_Err(wOpt, rotAxis[0],rotAxis[1],rotAxis[2])
	Optimize/Q/R=rotAxis/X=rotAxis/Y=(err) RefineOrientation_Err, wOpt
	err = RefineOrientation_Err(wOpt, rotAxis[0],rotAxis[1],rotAxis[2])
	KillWaves/Z W_OptGradient
	return err
End
//
Function RefineOrientation_Err(wOpt, rx,ry,rz)
	Wave wOpt
	Variable rx,ry,rz

	Variable N=DimSize(wOpt,0)
	Wave recipBase = decodeMatFromStr(note(wOpt))
	Make/N=(N,3)/D/FREE Ghats=wOpt[p][q], hkls=wOpt[p][q+3]

	Make/N=3/D/FREE rotAxis = {rx,ry,rz}
	Wave rotMat = RotMatAboutAxisOnly(rotAxis)
	MatrixOP/FREE qhats = NormalizeCols( rotMat x recipBase x (hkls^t) )^t

	MatrixOP/FREE sumDots = sum( qhats * Ghats )
	return abs(N-sumDots[0])
End


Static Function/WAVE MakePairsList(vecHats)// returns pairDots, whose columns are: {dot, i, j}
	Wave vecHats										// a wave containing unit vectors (only use first 3 columns)

	Variable notParallel = cos(0.05*PI/180)	// = 0.999999619228249, sufficient for parallel
	Variable Nv=DimSize(vecHats,0)
	Variable Np = Nv*(Nv-1)/2						// number of unique parirs of vecHats, init to max possible value
	Make/N=(Np,3)/FREE pairDots
	Make/N=3/D/FREE vecj, veci
	Variable dot, j,i,m=0
	for (j=0,m=0; j<(Nv-1); j+=1)
		vecj = vecHats[j][p]
		for (i=j+1;i<Nv;i+=1)
			veci = vecHats[i][p]
			dot = MatrixDot(vecj,veci)
			if (abs(dot)<notParallel)				// vectors are neither parallel nor anti-parallel
				pairDots[m][0] = dot				// store this pair, save dot and i,j
				pairDots[m][1] = j
				pairDots[m][2] = i
				m += 1			
			endif
		endfor
	endfor
	Np = m
	Redimension/N=(Np,-1) pairDots
	return pairDots
End


// find hkl that is paralllel to gVec, is allowed, and has energy < keVmax
Static Function/WAVE LowestAllowedHKLfromGhat(xtal,recip,gVec,keVmax,[kin])	// in and out can be the same wave
	STRUCT crystalStructure &xtal
	Wave recip						// reciprocal lattice (in actual rotated frame)
	Wave gVec						// q-vector that is known, gVec point in the measured direction
	Variable keVmax				// maximum allowed energy
	Wave kin							// optional incident wave direction (usually defaults to 001)

	Make/N=3/D/FREE ki={0,0,1}			// use this for the incident beam direction
	if (WaveExists(kin))
		ki = kin[p]
	endif
	normalize(ki)
	MatrixOP/FREE ghat = Normalize(gVec)

	Variable Qmax = -4*PI*MatrixDot(ki,ghat)*keVmax/hc	// |Q| = 4*¹*sin(theta)/lambda

	MatrixOP/FREE hklBase = Inv(recip) x ghat				// hkl vector but with maximum value of 1 or -1
	MatrixOP/FREE maxBase = maxVal(abs(hklBase))
	hklBase /= maxBase[0]											// now largest value of |hklBase| == 1

	MatrixOP/FREE qBase = recip x hklBase
	Variable maxInt=floor(Qmax/norm(qBase))		
	Make/N=3/D/FREE hkl=NaN

	Variable smallAdd = min(WaveType(gVec),WaveType(recip)) & 0x04 ? 1e-13 : 1e-6
	Variable i, ibest=0, dot, dotMax=-Inf
	for (i=1;i<=maxInt;i+=1)
		hkl = round(hklBase * i)
		if (!allowedHKL(hkl[0],hkl[1],hkl[2],xtal))
			continue
		endif
		MatrixOP/FREE qhat = Normalize(recip x hkl)
		dot = MatrixDot(qhat,ghat)
		if (dot > (dotMax+smallAdd))
			ibest = i
			dotMax = dot
		endif
	endfor

	if (ibest<1)
		return $""
	endif
	hkl = round(hklBase * ibest)
	hkl = hkl[p]==0 ? 0 : hkl[p]							// avoid "-0"
	return hkl
End


Static Function/WAVE FullPeakList2Ghats(maxSpots,FullPeakList,[FullPeakList1,FullPeakList2])	// convert peaks to a Qlist+intens, peaks file for analysis by Euler
	Variable maxSpots
	Wave FullPeakList
	Wave FullPeakList1
	Wave FullPeakList2

	if (numtype(maxSpots) || maxSpots<3)
		return $""
	endif
	if (ParamIsDefault(FullPeakList1))
		Wave FullPeakList1=$""
	endif
	if (ParamIsDefault(FullPeakList2))
		Wave FullPeakList2=$""
	endif

	if (!WaveExists(FullPeakList))
		DoAlert 0, "input wave for FullPeakList2File() does not exists"
		return $""
	elseif (DimSize(FullPeakList,1)!=11)
		DoAlert 0, "the passed full peak list '"+NameOfWave(FullPeakList)+"' is not the right size"
		return $""
	elseif (WaveExists(FullPeakList1))
		if (DimSize(FullPeakList1,1)!=11)
			DoAlert 0, "Full peak list 1'"+NameOfWave(FullPeakList1)+"' is the wrong size"
			return $""
		endif
	elseif (WaveExists(FullPeakList2))
		if (DimSize(FullPeakList2,1)!=11)
			DoAlert 0, "Full peak list 1'"+NameOfWave(FullPeakList2)+"' is the wrong size"
			return $""
		endif
	endif

	Variable N0=DimSize(FullPeakList,0), N1=0,N2=0
	if (WaveExists(FullPeakList1))
		N1 = DimSize(FullPeakList1,0)
	endif
	if (WaveExists(FullPeakList2))
		N2 = DimSize(FullPeakList2,0)
	endif
	Variable N=max(maxSpots,N0+N1+N2)
	if (N<1)													// nothing to write
		return $""
	endif

	STRUCT microGeometry geo							// note, dd and yc are reset from wave note below if it exists
	if (FillGeometryStructDefault(geo))			//fill the geometry structure with default values
		DoAlert 0, "no geometry structure found, did you forget to set it?"
		return $""
	endif
	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))
		DoAlert 0, "no lattice structure found, did you forget to set it?"
		return $""
	endif

	Make/N=(N,5)/D/FREE GhatMeasured				// temporary waves for the q info, holds Qhat and intensity
	SetDimLabel 1,0,Qx,GhatMeasured	;	SetDimLabel 1,1,Qy,GhatMeasured
	SetDimLabel 1,2,Qz,GhatMeasured	;	SetDimLabel 1,3,intensity,GhatMeasured
	SetDimLabel 1,4,dNum,GhatMeasured

	Make/N=3/D/FREE qhat
	String wnotePeak=note(FullPeakList), wnote
	Variable dNum=limit(detectorNumFromID(StringByKey("detectorID",wnotePeak,"=")),0,MAX_Ndetectors)
	Variable startx,groupx, starty,groupy			// ROI of the original image
	startx = NumberByKey("startx",wnotePeak,"=")
	groupx = NumberByKey("groupx",wnotePeak,"=")
	starty = NumberByKey("starty",wnotePeak,"=")
	groupy = NumberByKey("groupy",wnotePeak,"=")
	startx = numtype(startx) ? FIRST_PIXEL : startx
	groupx = numtype(groupx) ? 1 : groupx
	starty = numtype(starty) ? FIRST_PIXEL : starty
	groupy = numtype(groupy) ? 1 : groupy
	Variable depth = NumberByKey("depth",wnotePeak,"=")
	depth = numtype(depth) ? 0 : depth

	Variable i, px,py
	for (i=0,N=0;i<min(maxSpots,N0);i+=1)
		px = (startx-FIRST_PIXEL) + groupx*FullPeakList[i][0] + (groupx-1)/2		// change to un-binned pixels
		py = (starty-FIRST_PIXEL) + groupy*FullPeakList[i][1] + (groupy-1)/2		// pixels are still zero based
		pixel2q(geo.d[dNum],px,py,qhat,depth=depth)	// was in Wenge-coord system (OLD) or BeamLine(New)
		if (norm(qhat)>0)											// check for a valid Q
			GhatMeasured[N][0,2] = qhat[q]
			GhatMeasured[N][3] = FullPeakList[i][10]	// the intensity
			GhatMeasured[N][4] = dNum							// detector number
			N += 1
		endif
	endfor

	if (WaveExists(FullPeakList1) && N<maxSpots)
		wnote = note(FullPeakList1)
		dNum = limit(detectorNumFromID(StringByKey("detectorID",wnote,"=")),0,MAX_Ndetectors)
		startx = NumberByKey("startx",wnote,"=")
		groupx = NumberByKey("groupx",wnote,"=")
		starty = NumberByKey("starty",wnote,"=")
		groupy = NumberByKey("groupy",wnote,"=")
		startx = numtype(startx) ? FIRST_PIXEL : startx
		groupx = numtype(groupx) ? 1 : groupx
		starty = numtype(starty) ? FIRST_PIXEL : starty
		groupy = numtype(groupy) ? 1 : groupy
		depth = NumberByKey("depth",wnote,"=")
		depth = numtype(depth) ? 0 : depth
		for (i=0; i<N1 && N<maxSpots; i+=1)
			px = (startx-FIRST_PIXEL) + groupx*FullPeakList1[i][0] + (groupx-1)/2		// change to un-binned pixels
			py = (starty-FIRST_PIXEL) + groupy*FullPeakList1[i][1] + (groupy-1)/2		// pixels are still zero based
			pixel2q(geo.d[dNum],px,py,qhat,depth=depth)	// was in Wenge-coord system (OLD) or BeamLine(New)
			if (norm(qhat)>0)											// check for a valid Q
				GhatMeasured[N][0,2] = qhat[q]
				GhatMeasured[N][3] = FullPeakList1[i][10]	// the intensity
				GhatMeasured[N][4] = dNum							// detector number
				N += 1
			endif
		endfor
	endif

	if (WaveExists(FullPeakList2) && N<maxSpots)
		wnote = note(FullPeakList2)
		dNum = limit(detectorNumFromID(StringByKey("detectorID",wnote,"=")),0,MAX_Ndetectors)
		startx = NumberByKey("startx",wnote,"=")
		groupx = NumberByKey("groupx",wnote,"=")
		starty = NumberByKey("starty",wnote,"=")
		groupy = NumberByKey("groupy",wnote,"=")
		startx = numtype(startx) ? FIRST_PIXEL : startx
		groupx = numtype(groupx) ? 1 : groupx
		starty = numtype(starty) ? FIRST_PIXEL : starty
		groupy = numtype(groupy) ? 1 : groupy
		depth = NumberByKey("depth",wnote,"=")
		depth = numtype(depth) ? 0 : depth
		for (i=0; i<N2 && N<maxSpots; i+=1)
			px = (startx-FIRST_PIXEL) + groupx*FullPeakList2[i][0] + (groupx-1)/2		// change to un-binned pixels
			py = (starty-FIRST_PIXEL) + groupy*FullPeakList2[i][1] + (groupy-1)/2		// pixels are still zero based
			pixel2q(geo.d[dNum],px,py,qhat,depth=depth)	// was in Wenge-coord system (OLD) or BeamLine(New)
			if (norm(qhat)>0)											// check for a valid Q
				GhatMeasured[N][0,2] = qhat[q]
				GhatMeasured[N][3] = FullPeakList2[i][10]	// the intensity
				GhatMeasured[N][4] = dNum							// detector number
				N += 1
			endif
		endfor
	endif
	Redimension/N=(N,-1) GhatMeasured							// N is now the actual number of spots in GhatMeasured

	Variable val
	String str, wnoteQs
	if (StringMatch(StringByKey("WaveClass",wnotePeak,"="),"FittedPeakListTest*"))
		wnoteQs = "WaveClass=QsMeasuredTest;"
	else
		wnoteQs = "WaveClass=QsMeasured;"
	endif
	wnoteQs = ReplaceStringByKey("structureDesc",wnoteQs,xtal.desc,"=")
	sprintf str,"{ %g, %g, %g, %g, %g, %g }",xtal.a,xtal.b,xtal.c,xtal.alpha,xtal.beta,xtal.gam
	wnoteQs = ReplaceStringByKey("latticeParameters",wnoteQs,str,"=")
	wnoteQs = ReplaceStringByKey("lengthUnit",wnoteQs,"nm","=")
	wnoteQs = ReplaceNumberByKey("SpaceGroup",wnoteQs,xtal.SpaceGroup,"=")
	for (i=0;i<xtal.N;i+=1)
		sprintf str,"{%s,%g,%g,%g,%g}",xtal.atom[i].name,xtal.atom[i].x,xtal.atom[i].y,xtal.atom[i].z,xtal.atom[i].occ
		wnoteQs = ReplaceStringByKey("AtomDesctiption"+num2istr(i+1),wnoteQs,str,"=")
	endfor

	// copy a bunch of key=value pairs from wnotePeaks to wnote
	String commonKeys="imageFileName;exposure;dateExposed;rawIgorImage;fittedIgorImage;fractionBkg;"
	commonKeys += "minSpotSeparation;minPeakWidth;maxPeakWidth;totalPeakIntensity;totalIntensity;"
	commonKeys += "startx;groupx;endx;starty;groupy;endy;depth;"
//	commonKeys += "detectorID;"
	wnoteQs = MergeKeywordLists(wnoteQs,wnotePeak,0,"=",";",keys=commonKeys)
	wnoteQs = ReplaceStringByKey("peakListWave",wnoteQs,GetWavesDataFolder(FullPeakList,2),"=")

	Note/K GhatMeasured, wnoteQs
	return GhatMeasured
End


Static Function isNewDirection(Ghats,gvec,cosTol)
	// returns True if gvec is not parallel to any of the vectors in Ghat
	Wave Ghats
	Wave gvec
	Variable cosTol						// tolerance on dot:  cosTol = cos(dAnlge*PI/180),  cos(0.1¡)=0.999998476913288
	if (norm(gvec)==0)
		return 0								// skip the (000)
	endif
	MatrixOP/FREE dots = Ghats x gvec
	Variable dotMax = WaveMax(dots)
	return (dotMax < cosTol || numtype(dotMax))	// if Ghats filled with NaN then anything works
End


Static Function FindPeakInRots(PairRotations,rotAxis,thresh,center,radius)	//  returns topSum number of hits at this rot
	Wave PairRotations	// input list of rotation vectors for all possible pairs
	Wave rotAxis			// returned rotation vector (radian), the result
	Variable thresh		// threshold (radian)
	Wave center				// center of volume to check, only check PairRotations[i][3] within radius of center
	Variable radius		// radius around center to check, probably starts with Inf

	radius = numtype(radius) || radius<=0 ? Inf : radius
	Variable N=DimSize(PairRotations,0), i, thresh2=(thresh*thresh)

	//	inVolume[] identifies those rotations within radius of center[3]
	MatrixOP/FREE inVolume = greater(radius,sqrt(sumRows(magSqr(PairRotations-rowRepeat(center,N)))))
	Make/N=(N)/FREE/I summs=0
	Make/N=3/D/FREE rot
	for (i=0;i<N;i+=1)
		if (inVolume[i])							// PairRotations[i][p] is within radius of center
			rot = PairRotations[i][p]			// a single rotation vector from PairRotations[N][3]
			// sumi[0] is the number of rotations in PairRotations[N][3] that are within thresh of rot
			MatrixOP/FREE sumi = sum(greater(thresh2,sumRows(magSqr(PairRotations-rowRepeat(rot,N)))))
			summs[i] = sumi[0]					// number of rotations within thresh of PairRotations[i][3]
		endif
	endfor

	Variable topSum1 = WaveMax(summs)-1	// this ensures that NumTopRotations > 0
	MatrixOP/FREE isFrequentRot = greater(summs,topSum1)	// flags pairs having more than topSum1 matches


if (0)
	Variable NumTopRotations = sum(isFrequentRot)				// number of rotations that occur > topSum1 times
	MatrixOP/FREE rotAvg = sumCols(PairRotations*colRepeat(isFrequentRot,3))
	rotAxis = rotAvg / NumTopRotations		// the result being returned
	//printf "thresh = %g¡,  summs = [%g, %g]   NumTopRotations = %g,   topSum1 = %g\r",thresh*180/PI,WaveMin(summs),WaveMax(summs),NumTopRotations,topSum1
else
	// find position of largest point in summs, not the average of largest points
	WaveStats/M=1/Q summs
	Make/N=3/D/FREE rotMostLikely = PairRotations[V_maxloc][p]

	// sameRot is flags identifying all rotations within thresh of rotMostLikely
	MatrixOP/FREE sameRot = greater(thresh2,sumRows(magSqr(PairRotations-rowRepeat(rotMostLikely,N))))
	Variable NumSameRotations = sum(sameRot)			// number of rotations that occur at rotMostLikely

	MatrixOP/FREE rotAvg = sumCols(PairRotations*colRepeat(sameRot,3))
	rotAxis = rotAvg / NumSameRotations		// the result being returned

#if defined(ZONE_TESTING) || defined(QS_TESTING)
	// set global to hold index of first point at the max
	WaveStats/Q/M=1 sameRot
	Variable/G J_IndexOfMax = V_maxLoc
#endif
endif


#if defined(ZONE_TESTING) || defined(QS_TESTING)
	if (DataFolderExists(":GizmoWaves"))
		//	Duplicate/O summs, summsView
		Duplicate/O PairRotations, :GizmoWaves:PairRotationsViewSize
		Wave PairRotationsViewSize = :GizmoWaves:PairRotationsViewSize
		Redimension/S  PairRotationsViewSize
		PairRotationsViewSize = limit(round(summs[p]/1),2,25)

		Duplicate/O PairRotationsViewSize, :GizmoWaves:PairRotationsViewRGBA
		Wave PairRotationsViewRGBA = :GizmoWaves:PairRotationsViewRGBA
		Redimension/N=(-1,4) PairRotationsViewRGBA
		PairRotationsViewRGBA = 0
		PairRotationsViewRGBA[][3] = inVolume[p] ? 1 : 0.3
		PairRotationsViewRGBA[][1,2] = isFrequentRot[p] ? 1 : 0
	endif
#endif

	return topSum1+1								// largest number of parirs in one rotaiton
End


// each row of indexPairs generate two rotations (+ and -), so there are 1228 rotations to check
// what is the best way to check for lumps or groupings in the set of rotation vectors?
// this returns a list of rotation vectors (length = radian)
Static Function/WAVE ConstructAllPairRotations(indexPairs,MeasuredHats,TestHats)
	Wave indexPairs				// indexes into MeasuredHats and TestHats of pairs to define rotations
	Wave MeasuredHats				// only uses first 3 columns (ignores number and [uvw])
	Wave TestHats					// only uses first 3 columns (actually, only 3 columns are passed)

	Variable nPairs=DimSize(indexPairs,0)
	Make/N=(nPairs,3)/D/FREE g0A, g1A, q0A,q1A
	g0A = MeasuredHats[indexPairs[p][0]][q]
	g1A = MeasuredHats[indexPairs[p][1]][q]
	q0A = TestHats[indexPairs[p][2]][q]		// possible q's from given lattice
	q1A = TestHats[indexPairs[p][3]][q]

	Make/N=(nPairs,3)/D/FREE RotsHalf=NaN
	Make/N=(2*nPairs,3)/D/FREE PairRotations=NaN
	Make/N=(nPairs)/FREE mm						// really, just used to iterate with MultiThread

	//	MultiThread 
	MultiThread mm[] = rotFromPairsIndexed(p,g0A,g1A,q0A,q1A,RotsHalf)
	PairRotations[0,nPairs-1][] = RotsHalf[p][q]						// fill first half

#if defined(ZONE_TESTING) || defined(QS_TESTING)
mm = numtype(mm) ? 1 : 0
if (sum(mm))
	printf "found %d bad values in first half (out of %d)\r",sum(mm),DimSize(mm,0)
endif
#endif

	//	MultiThread 
	MultiThread mm[] = rotFromPairsIndexed(p,g1A,g0A,q0A,q1A,RotsHalf)
	PairRotations[nPairs,2*nPairs-1][] = RotsHalf[p-nPairs][q]	// fill second half

#if defined(ZONE_TESTING) || defined(QS_TESTING)
mm = numtype(mm) ? 1 : 0
if (sum(mm))
	printf "found %d bad values in first half (out of %d)\r",sum(mm),DimSize(mm,0)
endif
#endif

	return PairRotations
End
//
ThreadSafe Static Function rotFromPairsIndexed(m,g0A,g1A,q0A,q1A,rotAxisA)
	Variable m
	Wave g0A,g1A				// measured g's from data
	Wave q0A,q1A				// possible q's from given lattice
	Wave rotAxisA

	Make/N=3/D/FREE g0=g0A[m][p], g1=g1A[m][p]
	Make/N=3/D/FREE q0=q0A[m][p], q1=q1A[m][p]

	normalize(g0)
	normalize(g1)
	normalize(q0)
	normalize(q1)

	Wave rotAxis = rotFromPairs(g0,g1,q0,q1)
	rotAxisA[m][] = rotAxis[q]			// set the rotation
//	if (check_rotFromPairs(g0,g1,q0,q1,rotAxis,6e-4))
//		return NaN
//	endif
	return norm(rotAxis)
End
//
// ThreadSafe
ThreadSafe Static Function/WAVE rotFromPairs(g0,g1,q0,q1) // helper routine for rotFromPairsIndexed()
	Wave g0, g1					// measured g's from data
	Wave q0, q1					// possible q's from given lattice
	// want rotation that takes (q0,q1) --> (g0,g1)
	//		i.e. rot x q0 = g0,   and rot x q1 = g1
	// first make a q2 that forma a non-coplanar triplet with q0 & q1, then the same for g2 from g0 & g1
	// then what we want is:  rot x {q0,q1,q2} = {g0,g1,g2}
	// equivalent to:  rot = {g0,g1,g2} x {q0,q1,q2}^-1

	Make/N=3/D/FREE g2,q2, rotAxis
	Cross g0, g1							// make a g2 that forms three non-coplanar vectors
	Wave W_Cross=W_Cross
	g2 = W_Cross
	normalize(g2)

	Cross q0, q1							// make a q2 that forms three non-coplanar vectors
	q2 = W_Cross
	normalize(q2)
	KillWaves/Z W_Cross

	Make/N=(3,3)/D/FREE gmat, qmat
	gmat[][0] = g0[p]						// gmat = {g0,g1,g2}
	gmat[][1] = g1[p]
	gmat[][2] = g2[p]
	qmat[][0] = q0[p]						// qmat = {q0,q1,q2}
	qmat[][1] = q1[p]
	qmat[][2] = q2[p]
	MatrixOP/FREE rotMat = gmat x Inv(qmat)
//print " "
//printWave(rotMat,name="++ rotMat",brief=1)		// rotMat look right here, problem is in axisOfMatrixLocal()
//print " "

	Variable theta = axisOfMatrixLocal(rotMat,rotAxis,squareUp=1)*PI/180	// PROBLEM is in axisOfMatrixLocal()
	rotAxis *= theta

	return rotAxis
End

//  ============================== End of Common Indexing ==============================  //
//  ====================================================================================  //




//  ====================================================================================  //
//  ================================= Start of Utility =================================  //

Static Function genericIntensity(xtal,qvec,hkl,keV)
	STRUCT crystalStructure &xtal
	Wave qvec							// q-vector (length is actual Q-vector, NOT normalized)
	Wave hkl
	Variable keV

	Variable f1							// approx atomic form factor with Z=1 (this is VERY approximate, the atom position part is exact)
	Variable F2							// Structure factor for the cell squared
	Variable s2 = norm(qvec)^2	// s = 2*PI/d (1/nm))
	Variable/C Fcell

	f1 = 0.18 + 0.7*exp(-s2*0.00027)				// apporx. generic atomic form factor (with Z=1)
	Fcell = Fstruct(xtal,hkl[0],hkl[1],hkl[2],keV=keV)
	F2 = f1*f1 * magSqr(Fcell)

	// from Zachariasen eqn 3.72, Integrated Intensity = [re^2] * [ (1+cos^2(2theta))/2 ]  *  [ F^2 * lambda^4 / (2sin^2(theta)) ]
	// using s = 4*PI*sinTheta/lambda, last term becomes:
	//	F2*(lambda^2)*8*(PI^2)/(s^2)
	Variable lambda = hc/keV
	Variable intens = F2*(lambda*lambda)*8*(PI*PI)/s2	// intens = (Integrated Intensity)/(re^2), assuming polarization term = 1
	return intens
End



ThreadSafe Static Function/WAVE RotMatAboutAxisOnly(axis)		// return rotation matrix
	// see:   http://en.wikipedia.org/wiki/Rotation_matrix
	Wave axis			// axis of rotation, length is rotation angle (radian)

	Make/N=(3,3)/D/FREE rot
	Variable angle = norm(axis)		// start with angle (radian)
	Make/N=3/D/FREE a=axis/angle
	angle *= 180/PI						// convert angle to degree
	rotationMatAboutAxis(a,angle,rot)
	return rot
End
//ThreadSafe Static Function/WAVE RotMatAboutAxisOnly(axis)		// return rotation matrix
//	// see:   http://en.wikipedia.org/wiki/Rotation_matrix
//	Wave axis			// axis of rotation, length is rotation angle (radian)
//
//	Make/N=3/D/FREE a=axis
//	Variable angle = normalize(a)
//	Variable c=cos(angle), s=sin(angle)
//	Variable c1 = 1-c
//
//	Make/N=3/D/FREE rot
//	rot[0][0] = {c+a[0]*a[0]*c1, a[1]*a[0]*c1+a[2]*s, a[2]*a[0]*c1-a[1]*s}
//	rot[0][1] = {a[0]*a[1]*c1-a[2]*s, c+a[1]*a[1]*c1, a[2]*a[1]*c1+a[0]*s}
//	rot[0][2] = {a[0]*a[2]*c1+a[1]*s, a[1]*a[2]*c1-a[0]*s, c+a[2]*a[2]*c1}
//	return rot
//End



//		Test_axisOfMatrixLocal(0)
//		Test_axisOfMatrixLocal(1)
//		Test_axisOfMatrixLocal(2)
//		Test_axisOfMatrixLocal(3)
//		Test_axisOfMatrixLocal(4)
//		Test_axisOfMatrixLocal(5)
//		Test_axisOfMatrixLocal(6)
//		Test_axisOfMatrixLocal(7)
//		Test_axisOfMatrixLocal(8)
//		Test_axisOfMatrixLocal(9)
//Function Test_axisOfMatrixLocal(flag)
//	Variable flag
//
//	Make/N=(3,3)/D/FREE rot
//	Make/N=3/D/FREE correctAxis=NaN
//	Variable c=cos(10*PI/180), s=sin(10*PI/180), ir2=1/sqrt(2), correctAngle=NaN
//
//	if (flag==0)
//		rot = 0								// zero matrix, not a rotation, just fails
//	elseif (flag==1)
//		correctAxis = 0					// identity matrix, 0¡ rotation
//		correctAngle = 0
//		rot = (p==q)
//	elseif (flag==2)
//		correctAxis = {1,1,1}			// 120¡ about (111)
//		correctAngle = 120
//		rot[0][0] = {0,1,0}
//		rot[0][1] = {0,0,1}
//		rot[0][2] = {1,0,0}
//	elseif (flag==3)
//		correctAxis = {1,0,0}			// 180¡ about (001)
//		correctAngle = 180
//		rot[0][0] = {1,0,0}
//		rot[0][1] = {0,-1,0}
//		rot[0][2] = {0,0,-1}
//	elseif (flag==4)
//		correctAxis = {0,0,1}			// 180¡ about (001)
//		correctAngle = 180
//		rot[0][0] = {-1,0,0}
//		rot[0][1] = {0,-1,0}
//		rot[0][2] = {0,0,1}
//	elseif (flag==5)
//		correctAxis = {-1,0,1}			// 180 about (-1,0,1)
//		correctAngle = 180
//		rot[0][0] = {0,0,-1}
//		rot[0][1] = {0,-1,0}
//		rot[0][2] = {-1,0,0}
//	elseif (flag==6)
//		correctAxis = {1,0,0}			// +10¡ rotation about +x-axis
//		correctAngle = 10
//		rot[0][0] = {1,0,0}
//		rot[0][1] = {0,c,s}
//		rot[0][2] = {0,-s,c}
//	elseif (flag==7)
//		correctAxis = {1,0,0}			// +10¡ rotation about -x-axis
//		correctAngle = -10
//		rot[0][0] = {1,0,0}
//		rot[0][1] = {0,c,-s}
//		rot[0][2] = {0,s,c}
//	elseif (flag==8)
//		correctAxis = {0, 1, -1-sqrt(2)}	// 180¡ about this funny axis
//		correctAngle = 180
//		rot[0][0] = {-1,0,0}
//		rot[0][1] = {0,-ir2,-ir2}
//		rot[0][2] = {0,-ir2,ir2}
//	else
//		printf "Unknown flag number %g\r",flag
//		return 1
//	endif
//	if (correctAngle<0)
//		correctAxis = -correctAxis
//		correctAngle = -correctAngle
//	endif
//
//	printWave(rot,name="rot",brief=1)
//	MatrixOP/FREE delta = sqrt(sum(magSqr((rot x (rot^t)) - Identity(3))))
//	Variable determinant = MatrixDet(rot)
//	if (delta[0]>1e-12)
//		printf "ERROR -- rot is NOT a rotation matrix, |rot| = %g,   Æ = %g\r",determinant,delta[0]
//	elseif (determinant<0)
//		printf "WARNING -- |rot| = %g, NOT a proper rotation\r",determinant
//	endif
//
//	Make/N=3/D/FREE axis
//	// ******************************************************************
//	Variable angle = axisOfMatrixLocal(rot,axis)	// THIS IS THE TEST
//	// ******************************************************************
//	Make/N=3/D/FREE axisRad = axis*angle*PI/180
//	printf "axis = %s    angle = %g¡    scaled axis = %s\r",vec2str(axis,sep=", "),angle,vec2str(axisRad,places=12,sep=", ")
//	MatrixOP/FREE axisErr = sqrt(sum(magSqr(Normalize(correctAxis) - axis)))
//	if (abs(180-angle)<57e-12 && axisErr[0]>1e-12)			// 180¡, check the opposite rotation too
//		MatrixOP/FREE axisErrNeg = sqrt(sum(magSqr(Normalize(correctAxis) + axis)))
//		axisErr = min(axisErr,axisErrNeg)
//	endif
//
//	Make/N=(3,3)/D/FREE rotBackCalc
//	rotationMatAboutAxis(axis,angle,rotBackCalc)
//	MatrixOP/FREE rotErr = sqrt(sum(magSqr(rotBackCalc-rot)))
//	if (rotErr[0]<1e-12 && axisErr[0]<1e-12)
//		print "Success, matrix calculated from this axis matches rot above"
//		return 0
//	else
//		print " "
//		printf "ERROR -- Calculated axis = %s,  Correct axis = %s\r", vec2str(axis,sep=", "), vec2str(correctAxis,sep=", ")
//		printWave(rotBackCalc,name="         rotation calculated from axis does NOT match rot from above",brief=1, zeroThresh=1e-12)
//		return 1
//	endif
//End
// ThreadSafe 
ThreadSafe Static Function axisOfMatrixLocal(mat,axis,[squareUp,tol])
	// returns total rotation angle (deg), and sets axis to the axis of the total rotation
	Wave mat									// should be a rotation matrix
	Wave axis								// desired axis of the rotation (normalized)
	Variable squareUp						// optionally square up mat (default is NOT square up)
	Variable tol							// positive threshold for zero
	squareUp = ParamIsDefault(squareUp) ? NaN : squareUp
	squareUp = numtype(squareUp) ? 0 : !(!squareUp)
	tol = ParamIsDefault(tol) || numtype(tol) || tol<=0 ? NaN : tol
	if (numtype(tol))
		tol = WaveType(mat) & 0x04 ? 1e-13 : 1e-6
	endif

	Make/N=(3,3)/FREE/D rot=mat
	if (squareUp)
		if (SquareUpMatrix(rot))
			axis = NaN								// default for error
			return NaN
		endif
	else
		MatrixOp/FREE sumd = sum(Abs((mat x mat^t) - Identity(3)))
		if (sumd[0]<0 || sumd[0]>1e-4)		// not close enough to a rotation mat, an error
			axis = NaN								// default for error
			return NaN
		endif
	endif
	Variable cosine = (MatrixTrace(rot)-1)/2	// trace = 1 + 2*cos(theta)
	Variable angle = acos(limit(cosine,-1,1))	// angle (radian)

	// see:   http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToAngle/index.htm
	axis[0] = rot[2][1] - rot[1][2]			// first try the usual formula
	axis[1] = rot[0][2] - rot[2][0]			// this gives axis in right direction, but length is wrong
	axis[2] = rot[1][0] - rot[0][1]

	if (norm(axis)<tol)							// special for 0¡ or 180¡					
		if (cosine>0)								// special for 0¡ rotation					
			axis = 0									// (111) and zero
			angle = 0
		elseif (rot[0][0] > rot[1][1] && rot[0][0] > rot[2][2])	// special for some 180¡ rotations
			axis[0] = (1 + rot[0][0] - rot[1][1] - rot[2][2])
			axis[1] = (rot[0][1] + rot[1][0])
			axis[2] = (rot[0][2] + rot[2][0])
		elseif (rot[1][1] > rot[2][2])
			axis[0] = (rot[0][1] + rot[1][0])
			axis[1] = (1+ rot[1][1] - rot[0][0] - rot[2][2])
			axis[2] = (rot[1][2] + rot[2][1])
		else
			axis[0] = (rot[0][2] + rot[2][0])
			axis[1] = (rot[1][2] + rot[2][1])
			axis[2] = (1+ rot[2][2] - rot[0][0] - rot[1][1])
		endif
	endif
	normalize(axis)								// returns normalized axis
	return angle*180/PI							// return rotation angle (degree)
End
//Function Test_axisOfMatrix(flag)
//	Variable flag
//
//	Variable c,s
//	Make/N=(3,3)/D/FREE rot
//	if (flag==0)
//		rot[0][0] = {-1,0,0}
//		rot[0][1] = {0,-1,0}
//		rot[0][2] = {0,0,1}
//		printWave(rot,name="rot 180¡ about z-axis",brief=1)
//	elseif (flag==1)
//		rot[0][0] = {0,0,-1}
//		rot[0][1] = {0,-1,0}
//		rot[0][2] = {-1,0,0}
//		printWave(rot,name="Correct rot for this (q0,q1) & (g0,g1)",brief=1)
//	elseif (flag==2)
//		rot = 0
//		printWave(rot,name="rot = 0",brief=1)
//	elseif (flag==3)
//		rot = p==q
//		printWave(rot,name="Identity matrix",brief=1)
//	elseif (flag==4)
//		c=cos(10*PI/180)
//		s=sin(10*PI/180)
//		rot[0][0] = {1,0,0}
//		rot[0][1] = {0,c,s}
//		rot[0][2] = {0,-s,c}
//		printWave(rot,name="rot of 10¡ about x-axis",brief=1)
//	elseif (flag==5)
//		c=cos(-10*PI/180)
//		s=sin(-10*PI/180)
//		rot[0][0] = {1,0,0}
//		rot[0][1] = {0,c,s}
//		rot[0][2] = {0,-s,c}
//		printWave(rot,name="rot of -10¡ about x-axis",brief=1)
//	endif
//
//	print " "
//	Make/N=3/D/FREE axis
//	Variable angle = axisOfMatrix(rot,axis)
//	print "axis = ",vec2str(axis),"   angle = ",angle
//End


// This is just for completeness, DO NOT use Euler angles.
Static Function/WAVE rotAxis2EulerAngles(rotAxis)	// returns a vector[3] containing the three Euler angles (radian)
	Wave rotAxis
	//		Rz(gamma) * zhat = zhat
	//		and A = Rz(gamma)*Ry(beta)*Rz(alpha), and Rz*zhat = zhat
	//		A*zhat = Rz(gamma)*Ry(beta)*zhat
	//		for v=A*zhat, then v[2] gives cos(beta), and v[0],v[[1] gives azimuth
	//		then, Rz(gamma)*Ry(beta)*Rz(alpha) = A,  so
	//		Rz(alpha) = RyInv(beta)*RzInv(gamma) * A
	//		and Rz[0][0] = cos(angle), and Rz[0][1] = sin(angle)

	Variable alpha=NaN,bet=NaN,gam=NaN
	Wave A = RotMatAboutAxisOnly(rotAxis)

	Make/N=3/D/FREE zhat={0,0,1}
	MatrixOP/FREE v = A x zhat				// Rz(gamma) * zhat = zhat
	bet = acos(limit(v[2],-1,1))				// v[2] is v .dot. (001)
	gam = PI - atan2(v[1],v[0])				// azimuthal angle
	Wave RzInv = MatrixRz(-gam)				// Rz(-gamma) is Inverse(Rz(gamma))
	Wave RyInv = MatrixRy(-bet)
	// Rz(alpha) = RyInv(beta)*RzInv(gamma)*A
	MatrixOP/FREE Rzalpha = RyInv x RzInv x A	// Rzalpha = Rz(alpha),  Rzalpha = RyInv x mat
	alpha = atan2(Rzalpha[0][1],Rzalpha[0][0])	// atan2( sin(alpha), cos(alpha) )
	Make/N=3/D/FREE euler={alpha,bet,gam}
	return euler
End
//
Static Function/WAVE MatrixRz(angle)	// rotation matrix about the z axis thru angle (radian)
	Variable angle				// rotation angle (radian)
	Make/N=(3,3)/D/FREE Rz=0
	Variable cosine=cos(angle), sine=sin(angle)
	Rz[0][0] = cosine
	Rz[1][1] = cosine
	Rz[2][2] = 1;
	Rz[0][1] = sine
	Rz[1][0] = -sine
	return Rz
End
//
Static Function/WAVE MatrixRy(angle)	// rotation matrix about the y axis thru angle (radian)
	Variable angle				// rotation angle (radian)
	Make/N=(3,3)/D/FREE Ry=0
	Variable cosine=cos(angle), sine=sin(angle)
	Ry[0][0] = cosine;
	Ry[1][1] = 1;
	Ry[2][2] = cosine;
	Ry[2][0] = sine;
	Ry[0][2] = -sine;
	return Ry
End

//  ================================== End of Utility ==================================  //
//  ====================================================================================  //




//  ====================================================================================  //
//  ================================ Start of Test Data ================================  //

// make the wave GhatsMeasured, GhatsMeasured are the result of finding peaks on a detector.
// This wave is a good set of test data for testing all of this stuff.
Function/WAVE MakeSimulatedTestPattern(Nreq,[cone,lattice,angle,axis,tol,keVmax,printIt])
	Variable Nreq						// number of test spots to make
	Variable cone						// cone angle (degree)
	String lattice						// "C","M1",M2"
	Variable angle						// rotation angle from lattice about axis (degree)
	String axis							// string representation of axis for angle rotation
	Variable tol						// angle (degree), usually 0.1¡, reject hkl's that differ by less than tol
	Variable keVmax					// max energy of a reflection (keV)
	Variable printIt
	cone = ParamIsDefault(cone) || numtype(cone) || cone<=0 ? 20 : cone
	lattice = SelectString(ParamIsDefault(lattice),lattice,"C")
	angle = ParamIsDefault(angle) || numtype(angle) ? 0 : angle
	axis = SelectString(ParamIsDefault(axis),axis,"")
	tol = ParamIsDefault(tol) ? 0.1 : tol
	keVmax = ParamIsDefault(keVmax) || numtype(keVmax) || keVmax<=0 ? 30 : keVmax
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : !(!printIt)

	Wave axis3 = str2vec(axis)
	if (Nreq<=0 || numtype(Nreq))
		Nreq = 98
		Prompt Nreq,"Number Ghat's to make"
		Prompt cone,"cone angle (0 gives all)"
		Prompt lattice,"Test lattice to construct",popup,"Cubic (C);Monoclinic 1 (M1);Monoclinic 2 (M2)"
		Prompt angle, "rotation angle applied to lattice"
		if (strlen(axis)<1)
			axis="1,2,3"
		endif
		Prompt axis, "axis of test rotation applied to lattice"
		DoPrompt "Number of G's",Nreq,cone,lattice,angle,axis
		if (V_flag)
			return $""
		endif

		Variable i=strsearch(lattice,"(",0)
		lattice = lattice[i+1,Inf]
		i = strsearch(lattice,")",0)
		lattice = lattice[0,i-1]

		angle = numtype(angle) ? 0 : angle
		if (angle)
			Wave axis3 = str2vec(axis)
		endif
		printIt = 1
	endif
	if (Nreq<=0 || numtype(Nreq))
		return $""
	endif
	if (printIt)
		printf "MakeSimulatedTestPattern(%g",Nreq
		if (!ParamIsDefault(lattice))
			printf ", lattice=\"%s\"",lattice
		endif
		if (!ParamIsDefault(cone) || cone>0)
			printf ", cone=%g",cone
		endif
		if (!ParamIsDefault(angle) || !(angle==0))
			printf ", angle=%g, axis=\"%s\"",angle,vec2str(axis3,sep=",",bare=1)
		endif
		if (!ParamIsDefault(tol) || tol!=0.1)
			printf ", tol=%g",tol
		endif
		if (!ParamIsDefault(printIt))
			printf ", printIt=%g",printIt
		endif
		printf ")\r"
	endif
	if (numtype(cone+tol) || tol<=0 || abs(cone>180))
		return $""
	elseif (angle)
		if (!WaveExists(axis3))
			return $""
		elseif (!(norm(axis3)>0))
			return $""
		elseif (numpnts(axis3)!=3)
			return $""
		endif
	endif

	Variable cosTol=cos(tol*PI/180)		// convert degree to dot, this should be slightly less than 1.0
	Make/N=(3,3)/D/FREE direct
	String FullPeakListName, desc=""
	if (StringMatch(lattice,"C"))
		direct = (p==q) * 0.54310206		// simple cubic direct lattice with a = 0.3 nm
//		Variable SpaceGroup = 224			// [195,230]
		Variable SpaceGroup = 227			// [195,230], diamond structure
		desc = "Cubic Test"
		FullPeakListName = "FullPeakListTestCubic"
	elseif (StringMatch(lattice,"M1"))
		direct[0][0] = {0.28, 0, 0.076}	// a   a non-cubic direct lattice (Monoclinic)
		direct[0][1] = {0, 0.3, 0}		// b	(a^c !=90, a^b = c^b = 90,  alpha=gamma=90
		direct[0][2] = {0, 0, 0.29}		// c
		SpaceGroup = 14						// [3,15]
		desc = "Monoclinic 1 Test"
		FullPeakListName = "FullPeakListTestM1"
	elseif (StringMatch(lattice,"M2"))
		direct[0][0] = {0.28, 0, 0.075}	// a   a non-cubic direct lattice (Monoclinic)
		direct[0][1] = {0, 0.3, 0}		// b	(a^c !=90, a^b = c^b = 90,  alpha=gamma=90
		direct[0][2] = {0, 0, 0.29}		// c
		SpaceGroup = 14						// [3,15]
		desc = "Monoclinic 2 Test"
		FullPeakListName = "FullPeakListTestM2"
	else
		return $""
	endif

	Variable Vc = MatrixDet(direct)
	LatticeSym#SetSymOpsForSpaceGroup(SpaceGroup)

	MatrixOP/FREE/O recip   = 2*PI * (Inv(direct))^t
	recip = recip==0 ? 0 : recip
	String directStr=encodeMatAsStr(direct), recipStr=encodeMatAsStr(recip)

	if (angle)
		MatrixOP/FREE axisNorm = Normalize(axis3)
		Duplicate/FREE axisNorm, axisTemp
		axisTemp *= angle*PI/180						// length of axisTemp is rot in radians
		Wave rot = RotMatAboutAxisOnly(axisTemp)
		MatrixOP/FREE/O recip = rot x recip
		MatrixOP/FREE/O direct = rot x direct
		direct = direct == 0 ? 0 : direct			// get rid of -0
		recip = recip == 0 ? 0 : recip
		Make/N=3/D/FREE scaledRotAxis = axisNorm*(angle*PI/180)
		if (printIt)
			printf "test rotation is %g¡ about the (%s)  =  %s\r",angle,vec2str(axis3,sep=" ",bare=1),vec2str(scaledRotAxis)
		endif
	endif
	if (!CheckDirectRecip(direct,recip))
		print "Invalid computation of direct & recip"
		Abort "Invalid computation of direct & recip"
		return $""
	endif

	Wave LC = direct2LatticeConstants(direct)
	if (printIt)
		printf "lattice constants: %.9g, %.9g, %.9g,   %g¡, %g¡, %g¡,   SG=%G\r",LC[0],LC[1],LC[2],LC[3],LC[4],LC[5],SpaceGroup
	endif
	STRUCT crystalStructure xtal
	xtal.a = LC[0] ;			xtal.b = LC[1] ;		xtal.c = LC[2]
	xtal.alpha = LC[3] ;	xtal.beta = LC[4] ;	xtal.gam = LC[5]
	xtal.SpaceGroup = SpaceGroup
	xtal.desc = desc
	if (!isValidLatticeConstants(xtal)	)
		print "ERROR -- Lattice Constants are NOT valid"
		return $""
	endif
	xtal.Unconventional00 = NaN
	xtal.sourceFile = ""
	LatticeSym#CleanOutCrystalStructure(xtal)
	LatticeSym#ForceLatticeToStructure(xtal)
	MakeSymmetryOps(xtal)					// make a wave with the symmetry operation
	UpdateCrystalStructureDefaults(xtal)

	Make/N=3/D/FREE qUp=0
	Variable cosCone=-2						// -2 includes everything
	if (cone>0)
		qUp = {0,1,-1}							// direction of q-vector that will diffract straight up
		normalize(qUp)
		cosCone = cos(cone*PI/180)
	endif

	STRUCT microGeometry geo							// note, dd and yc are reset from wave note below if it exists
	if (FillGeometryStructDefault(geo))			//fill the geometry structure with default values
		DoAlert 0, "no geometry structure found, did you forget to set it?"
		return $""
	endif

#if defined(ZONE_TESTING) || defined(QS_TESTING)
	Make/N=(Nreq,3)/D/O GhatsMeasured=0, GhatsMeasuredHKL=0 // this holds the result
#else
	Make/N=(Nreq,3)/D/FREE GhatsMeasured=0, GhatsMeasuredHKL=0
#endif
	Make/N=(Nreq,3)/D/FREE Ghat3=0		// only need this for the isNewDirection() test
	Make/N=3/D/FREE hkl, ki={0,0,1}
	Make/N=(Nreq,11)/D/O $FullPeakListName/WAVE=FullPeakListTest = NaN
	Variable keV, intens, intensMax=0
	Variable/C pz
	Variable h,k,l, N							// N is the actual number of unique spots saved
	Variable Nx=2048, Ny=2048, px,py
	Variable hklMax, hklMaxMax=20
	for (hklMax=0,N=0; hklMax<=hklMaxMax && N<Nreq; hklMax+=1)
		for (l=0; l<=hklMax && N<Nreq; l=LatticeSym#incrementIndex(l))
			hkl[2] = l
			for (k=0; k<=hklMax && N<Nreq; k=LatticeSym#incrementIndex(k))
				hkl[1] = k
				for (h=0; h<=hklMax && N<Nreq; h=LatticeSym#incrementIndex(h))
					hkl[0] = h
					if (!allowedHKL(h,k,l,xtal))
						continue
					endif
					MatrixOP/FREE/O qvec = recip x hkl
					MatrixOP/FREE/O qhat = Normalize(qvec)
					pz = q2pixel(geo.d[0],qhat)
					px = real(pz)
					py = imag(pz)
					if (px<0 || px>=Nx || py<0 || py>=Ny || numtype(px+py))	// misses the detector
						continue
					endif
					MatrixOP/FREE Qmag = sqrt(sum(magSqr(recip x hkl)))
					qvec = qhat * Qmag[0]
					keV = -hc*Qmag[0] / (4*PI*MatrixDot(ki,qhat))	// Q = 4*PI*sin(theta)/lambda
					if (isNewDirection(Ghat3,qhat,cosTol) && MatrixDot(qhat,qUp)>cosCone || keV>keVmax)
						FullPeakListTest[N][0] = px
						FullPeakListTest[N][1] = py
						intens = genericIntensity(xtal,qvec,hkl,keV)
						intensMax = max(intens,intensMax)
						FullPeakListTest[N][10] = intens

						GhatsMeasured[N][] = qhat[q]		// a new direction, save it...
						GhatsMeasuredHKL[N][] = hkl[q]	// hkl associated with each test Ghat
						Ghat3[N][] = qhat[q]			// Ghat3 only needed for isNewDirection()
						N += 1
					endif
				endfor
			endfor
		endfor
	endfor
	Redimension/N=(N,-1) GhatsMeasured,GhatsMeasuredHKL, Ghat3	// trim to actual size
	Redimension/N=(N,-1) FullPeakListTest

	// find the central hkl
	MatrixOP/FREE GhatAvg = Normalize(sumCols(Ghat3))
	Redimension/N=3 GhatAvg
	MatrixOP/FREE/O hkl = Normalize(Inv(recip) x GhatAvg)
	Variable small = smallestNonZeroValue(hkl,tol=1e-8)
	hkl *= 12/small
	removeCommonFactors(hkl)
	FullPeakListTest[][10] /= intensMax

	String wnote="waveClass=FittedPeakListTest;"
	wnote = ReplaceStringByKey("direct",wnote,directStr,"=")
	wnote = ReplaceStringByKey("recip",wnote,recipStr,"=")
	wnote = ReplaceNumberByKey("SpaceGroup",wnote,SpaceGroup,"=")
	wnote = ReplaceNumberByKey("cone",wnote,cone,"=")
	wnote = ReplaceNumberByKey("tol",wnote,tol,"=")
	wnote = ReplaceStringByKey("hklCenter",wnote,vec2str(hkl,bare=1,sep=","),"=")
	if (angle)
		wnote = ReplaceNumberByKey("rotAngle",wnote,angle,"=")
		wnote = ReplaceStringByKey("rotAxis",wnote,vec2str(axis3,bare=1,sep=","),"=")
	endif
	Note/K GhatsMeasured, wnote
	wnote = ReplaceNumberByKey("startx",wnote,0,"=")
	wnote = ReplaceNumberByKey("endx",wnote,2047,"=")
	wnote = ReplaceNumberByKey("groupx",wnote,1,"=")
	wnote = ReplaceNumberByKey("starty",wnote,0,"=")
	wnote = ReplaceNumberByKey("endy",wnote,2047,"=")
	wnote = ReplaceNumberByKey("groupy",wnote,1,"=")
	wnote = ReplaceNumberByKey("xdim",wnote,geo.d[0].Nx,"=")
	wnote = ReplaceNumberByKey("ydim",wnote,geo.d[0].Nx,"=")
	wnote = ReplaceNumberByKey("xDimDet",wnote,geo.d[0].Nx,"=")
	wnote = ReplaceNumberByKey("yDimDet",wnote,geo.d[0].Nx,"=")
	wnote = ReplaceStringByKey("detectorID",wnote,"PE1621 723-3335","=")
	wnote = ReplaceStringByKey("detectorModel",wnote,"XRD0820","=")
	wnote = ReplaceStringByKey("beamLine",wnote,"34ID-E","=")
	wnote = ReplaceNumberByKey("exposure",wnote,1,"=")
	wnote = ReplaceStringByKey("MonoMode",wnote,"white slitted","=")
	Note/K FullPeakListTest,wnote
	SetDimLabel 1,0,x0,FullPeakListTest		;	SetDimLabel 1,1,y0,FullPeakListTest
	SetDimLabel 1,2,x0Err,FullPeakListTest	;	SetDimLabel 1,3,y0Err,FullPeakListTest
	SetDimLabel 1,4,fwx,FullPeakListTest		;	SetDimLabel 1,5,fwy,FullPeakListTest
	SetDimLabel 1,6,fwxErr,FullPeakListTest	;	SetDimLabel 1,7,fwyErr,FullPeakListTest
	SetDimLabel 1,8,correlation,FullPeakListTest ;	SetDimLabel 1,9,correlationErr,FullPeakListTest
	SetDimLabel 1,10,area,FullPeakListTest
	if (printIt)
		printf "maximum hkl went to (%d %d %d)\r",hklMax,hklMax,hklMax
		printf "Made %d simulated peaks, stored in '%s'\r",N,NameOfWave(FullPeakListTest)
	endif
	return FullPeakListTest
End

Function GraphTestPattern(FullPeakList,FullPeakIndexed,[pattern])
	Wave FullPeakList
	Wave FullPeakIndexed
	Variable pattern
	pattern = ParamIsDefault(pattern) || numtype(pattern) ? 0 : pattern

	if (!WaveExists(FullPeakList))
		String wList=WaveListClass("FittedPeakList*","*","MINCOLS:11")
		if (ItemsInList(wList)<1)
			return 1
		elseif (ItemsInList(wList)==1)
			Wave FullPeakList=$StringFromList(0,wList)
		else
			String name
			Prompt name,"FullPeakList",popup,wList
			DoPrompt "FullPeakList", name
			if (V_flag)
				return 1
			endif
			Wave FullPeakList = $name
		endif
	endif
	String win = StringFromList(0,FindGraphsWithWave(FullPeakList))
	if (strlen(win))
		DoWindow/F $win
		return 0
	endif
	if (!WaveExists(FullPeakIndexed))
		Wave FullPeakIndexed = FindIndexListForPeakList(FullPeakList)
	endif

	Display/W=(80,69,794,726)/N=GraphTestPattern/K=1 FullPeakList[*][1] vs FullPeakList[*][0]
	if (WaveExists(FullPeakIndexed))
		AppendToGraph FullPeakIndexed[*][10] vs FullPeakIndexed[*][9]
	endif
	ModifyGraph width={Aspect,1}, tick=2, mirror=1, minor=1, lowTrip=0.001, mode=3, mrkThick=2
	ModifyGraph marker[0]=7, msize[0]=5
	if (WaveExists(FullPeakIndexed))
		ModifyGraph marker[1]=8, rgb[1]=(2,39321,1), msize[1]=8
	endif

	String wnote=note(FullPeakList)
	Variable dNum=limit(detectorNumFromID(StringByKey("detectorID",wnote,"=")),0,MAX_Ndetectors)
	Variable startx=NumberByKey("startx",wnote,"="), starty=NumberByKey("starty",wnote,"=")
	Variable groupx=NumberByKey("groupx",wnote,"="), groupy=NumberByKey("groupy",wnote,"=")
	Variable Nx=NumberByKey("xdim",wnote,"="), Ny=NumberByKey("ydim",wnote,"=")
	SetAxis/R left Nx,starty
	SetAxis bottom startx,Ny

	STRUCT microGeometry geo						// note, dd and yc are reset from wave note below if it exists
	if (!FillGeometryStructDefault(geo))		//fill the geometry structure with default values
		Make/N=3/FREE xyz = {0,1,-1}				// Q of surface normal
		Variable/C pz =  q2pixel(geo.d[dNum],xyz)		// pixel location of surface normal
		Variable px=real(pz), py=imag(pz)	//, Nx=geo.d[dNum].Nx, Ny=geo.d[dNum].Ny
		if (numtype(px+py))								// perhaps the beam is on path of direct beam
			xyz = {0,0, geo.d[dNum].P[2]}
			XYZ2pixel(geo.d[dNum],xyz,px,py)
		endif
		if (px==limit(px,0,Nx-1) && py==limit(py,0,Ny-1))
			Variable NxFW=Nx/10, NyFW=Ny/10
			if (numtype(startx+starty+groupx+groupy)==0 && (groupx>1 || groupy>1))	// binned image
				px = round(( px-startx-(groupx-1)/2 )/groupx)	// pixel is zero based here & startx is zero based
				py = round(( py-starty-(groupy-1)/2 )/groupy)	// groupx=1 is un-binned
				NxFW /= groupx
				NyFW /= groupy
			endif
			DrawMarker(px,py,round(NxFW),round(NyFW),"cross gap",dash=2,layer="UserAxes")
			String str
			sprintf str,"%g,%g",px,py
			SetWindow kwTopWin, userdata(pixelCenter)=str
		endif
	endif

	String list = GetUserData("","","Indexing")
	if (WaveExists(FullPeakIndexed))
		list = ReplaceStringByKey("FullPeakIndexed",list,GetWavesDataFolder(FullPeakIndexed,2),"=")
	endif
	list = ReplaceNumberByKey("patternNum",list,pattern,"=")
	SetWindow kwTopWin,hook(peakInfo)=getFittedPeakInfoHook
	SetWindow kwTopWin,userdata(Indexing)=list
	SetWindow kwTopWin userdata(FitPeaks) = "FullPeakList="+GetWavesDataFolder(FullPeakList,2)
	return 0
End
//
Static Function/WAVE FindIndexListForPeakList(FullPeakList)	// returns wave ref of FullPeakIndex for an image
	Wave FullPeakList
	String wList=WaveListClass("IndexedPeakList*","*","MAXCOLS:12")
	Variable i, N=ItemsInList(wList)
	for (i=0;i<N;i+=1)
		Wave indexWave = $StringFromList(i,wList)
		Wave pkTest = $StringByKey("peakListWave",note(indexWave),"=")
		if (WaveRefsEqual(FullPeakList,pkTest))
			return indexWave
		endif
	endfor
	return $""
End



// make the wave GhatsMeasured, GhatsMeasured are the result of finding peaks on a detector.
// This wave is a good set of test data for testing all of this stuff.
Static Function/WAVE MakeTestSpots(Nreq,[cone,lattice,angle,axis,tol,printIt])
	Variable Nreq						// number of test spots to make
	Variable cone						// cone angle (degree)
	String lattice						// "C","M1",M2"
	Variable angle						// rotation angle from lattice about axis (degree)
	String axis							// string representation of axis for angle rotation
	Variable tol						// angle (degree), usually 0.1¡, reject hkl's that differ by less than tol
	Variable printIt
	cone = ParamIsDefault(cone) || numtype(cone) || cone<=0 ? 20 : cone
	lattice = SelectString(ParamIsDefault(lattice),lattice,"C")
	angle = ParamIsDefault(angle) || numtype(angle) ? 0 : angle
	axis = SelectString(ParamIsDefault(axis),axis,"")
	tol = ParamIsDefault(tol) ? 0.1 : tol
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? 0 : !(!printIt)

	Wave axis3 = str2vec(axis)
	if (Nreq<=0 || numtype(Nreq))
		Nreq = 98
		Prompt Nreq,"Number Ghat's to make"
		Prompt cone,"cone angle (0 gives all)"
		Prompt lattice,"Test lattice to construct",popup,"Cubic (C);Monoclinic 1 (M1);Monoclinic 2 (M2)"
		Prompt angle, "rotation angle applied to lattice"
		if (strlen(axis)<1)
			axis="1,2,3"
		endif
		Prompt axis, "axis of test rotation applied to lattice"
		DoPrompt "Number of G's",Nreq,cone,lattice,angle,axis
		if (V_flag)
			return $""
		endif

		Variable i=strsearch(lattice,"(",0)
		lattice = lattice[i+1,Inf]
		i = strsearch(lattice,")",0)
		lattice = lattice[0,i-1]

		angle = numtype(angle) ? 0 : angle
		if (angle)
			Wave axis3 = str2vec(axis)
		endif
		printIt = 1
	endif
	if (Nreq<=0 || numtype(Nreq))
		return $""
	endif
	if (printIt)
		printf "MakeTestSpots(%g",Nreq
		if (!ParamIsDefault(lattice))
			printf ", lattice=\"%s\"",lattice
		endif
		if (!ParamIsDefault(cone) || cone>0)
			printf ", cone=%g",cone
		endif
		if (!ParamIsDefault(angle) || !(angle==0))
			printf ", angle=%g, axis=\"%s\"",angle,vec2str(axis3,sep=",",bare=1)
		endif
		if (!ParamIsDefault(tol) || tol!=0.1)
			printf ", tol=%g",tol
		endif
		if (!ParamIsDefault(printIt))
			printf ", printIt=%g",printIt
		endif
		printf ")\r"
	endif
	if (numtype(cone+tol) || tol<=0 || abs(cone>180))
		return $""
	elseif (angle)
		if (!WaveExists(axis3))
			return $""
		elseif (!(norm(axis3)>0))
			return $""
		elseif (numpnts(axis3)!=3)
			return $""
		endif
	endif

	Variable cosTol=cos(tol*PI/180)		// convert degree to dot, this should be slightly less than 1.0
	Make/N=(3,3)/D/FREE direct
	String desc=""

	if (StringMatch(lattice,"C"))
		direct = (p==q) * 0.54310206		// simple cubic direct lattice with a = 0.3 nm
		Variable SpaceGroup = 224			// [195,230]
SpaceGroup = 227			// [195,230]

		desc = "Cubic Test"
	elseif (StringMatch(lattice,"M1"))
		direct[0][0] = {0.28, 0, 0.076}	// a   a non-cubic direct lattice (Monoclinic)
		direct[0][1] = {0, 0.3, 0}		// b	(a^c !=90, a^b = c^b = 90,  alpha=gamma=90
		direct[0][2] = {0, 0, 0.29}		// c
		SpaceGroup = 14						// [3,15]
		desc = "Monoclinic 1 Test"
	elseif (StringMatch(lattice,"M2"))
		direct[0][0] = {0.28, 0, 0.075}	// a   a non-cubic direct lattice (Monoclinic)
		direct[0][1] = {0, 0.3, 0}		// b	(a^c !=90, a^b = c^b = 90,  alpha=gamma=90
		direct[0][2] = {0, 0, 0.29}		// c
		SpaceGroup = 14						// [3,15]
		desc = "Monoclinic Test"
	endif

	Variable Vc = MatrixDet(direct)
	LatticeSym#SetSymOpsForSpaceGroup(SpaceGroup)

	MatrixOP/FREE/O recip   = 2*PI * (Inv(direct))^t
	recip = recip==0 ? 0 : recip
	String directStr=encodeMatAsStr(direct), recipStr=encodeMatAsStr(recip)

	if (angle)
		MatrixOP/FREE axisNorm = Normalize(axis3)
		Duplicate/FREE axisNorm, axisTemp
		axisTemp *= angle*PI/180						// length of axisTemp is rot in radians
		Wave rot = RotMatAboutAxisOnly(axisTemp)
		MatrixOP/FREE/O recip = rot x recip
		MatrixOP/FREE/O direct = rot x direct
		direct = direct == 0 ? 0 : direct			// get rid of -0
		recip = recip == 0 ? 0 : recip
		Make/N=3/D/FREE scaledRotAxis = axisNorm*(angle*PI/180)
		if (printIt)
			printf "test rotation is %g¡ about the (%s)  =  %s\r",angle,vec2str(axis3,sep=" ",bare=1),vec2str(scaledRotAxis)
		endif
	endif
	if (!CheckDirectRecip(direct,recip))
		print "Invalid computation of direct & recip"
		Abort "Invalid computation of direct & recip"
		return $""
	endif

	Wave LC = direct2LatticeConstants(direct)
	if (printIt)
		printf "lattice constants: %.9g, %.9g, %.9g,   %g¡, %g¡, %g¡,   SG=%G\r",LC[0],LC[1],LC[2],LC[3],LC[4],LC[5],SpaceGroup
	endif
	STRUCT crystalStructure xtal
	xtal.a = LC[0] ;			xtal.b = LC[1] ;		xtal.c = LC[2]
	xtal.alpha = LC[3] ;	xtal.beta = LC[4] ;	xtal.gam = LC[5]
	xtal.SpaceGroup = SpaceGroup
	xtal.desc = desc
	if (!isValidLatticeConstants(xtal)	)
		print "ERROR -- Lattice Constants are NOT valid"
		return $""
	endif
	xtal.Unconventional00 = NaN
	xtal.sourceFile = ""
	LatticeSym#CleanOutCrystalStructure(xtal)
	LatticeSym#ForceLatticeToStructure(xtal)
	MakeSymmetryOps(xtal)					// make a wave with the symmetry operation
	UpdateCrystalStructureDefaults(xtal)

	Make/N=3/D/FREE qUp=0
	Variable cosCone=-2						// -2 includes everything
	if (cone>0)
		qUp = {0,1,-1}							// direction of q-vector that will diffract straight up
		normalize(qUp)
		cosCone = cos(cone*PI/180)
	endif

	Make/N=(Nreq,3)/D/FREE Ghat3=0		// only need this for the isNewDirection() test
	Make/N=(Nreq,3)/D/O GhatsMeasured=0 // this holds the result
	Make/N=(Nreq,3)/D/O GhatsMeasuredHKL=0 // this holds the result
	Make/N=3/D/FREE hkl
	Variable h,k,l, N	=0						// N is the actual number of unique spots saved
	Variable hklMax, hklMaxMax=20
	for (hklMax=0,N=0; hklMax<=hklMaxMax && N<Nreq; hklMax+=1)
		for (l=0; l<=hklMax && N<Nreq; l=LatticeSym#incrementIndex(l))
			hkl[2] = l
			for (k=0; k<=hklMax && N<Nreq; k=LatticeSym#incrementIndex(k))
				hkl[1] = k
				for (h=0; h<=hklMax && N<Nreq; h=LatticeSym#incrementIndex(h))
					hkl[0] = h
					MatrixOP/FREE/O qvec = Normalize(recip x hkl)
					if (isNewDirection(Ghat3,qvec,cosTol) && MatrixDot(qvec,qUp)>cosCone)
						GhatsMeasured[N][] = qvec[q]		// a new direction, save it...
						GhatsMeasuredHKL[N][] = hkl[q]	// hkl associated with each test Ghat
						Ghat3[N][] = qvec[q]			// Ghat3 only needed for isNewDirection()
						N += 1
					endif
				endfor
			endfor
		endfor
	endfor
	Redimension/N=(N,-1) GhatsMeasured,GhatsMeasuredHKL, Ghat3	// trim to actual size
	if (printIt)
		printf "maximum hkl went to (%d %d %d)\r",hklMax,hklMax,hklMax
		printf "Made %d unique directions, stored in '%s'\r",N,NameOfWave(GhatsMeasured)
	endif
	SetDimLabel 1,0,Gx,GhatsMeasured		// columns are the zone axis direction (unit vector)
	SetDimLabel 1,1,Gy,GhatsMeasured
	SetDimLabel 1,2,Gz,GhatsMeasured
	SetDimLabel 1,0,H,GhatsMeasuredHKL		// columns are hkl
	SetDimLabel 1,1,K,GhatsMeasuredHKL
	SetDimLabel 1,2,L,GhatsMeasuredHKL

	// find the central hkl
	MatrixOP/FREE GhatAvg = Normalize(sumCols(Ghat3))
	Redimension/N=3 GhatAvg
	MatrixOP/FREE/O hkl = Normalize(Inv(recip) x GhatAvg)
	Variable small = smallestNonZeroValue(hkl,tol=1e-8)
	hkl *= 12/small
	removeCommonFactors(hkl)

	String wnote="waveClass=QsMeasuredTest;"
	wnote = ReplaceStringByKey("direct",wnote,directStr,"=")
	wnote = ReplaceStringByKey("recip",wnote,recipStr,"=")
	wnote = ReplaceNumberByKey("SpaceGroup",wnote,SpaceGroup,"=")
	wnote = ReplaceNumberByKey("cone",wnote,cone,"=")
	wnote = ReplaceNumberByKey("tol",wnote,tol,"=")
	wnote = ReplaceStringByKey("hklCenter",wnote,vec2str(hkl,bare=1,sep=","),"=")
	if (angle)
		wnote = ReplaceNumberByKey("rotAngle",wnote,angle,"=")
		wnote = ReplaceStringByKey("rotAxis",wnote,vec2str(axis3,bare=1,sep=","),"=")
	endif
	Note/K GhatsMeasured, wnote
	return GhatsMeasured
End


Static Function removeCommonFactors(vec)
	Wave vec									// vector of integers, remove all common factors

	Variable N=numpnts(vec)
	Make/N=(N)/D/FREE ivec = round(vec)

	MatrixOP/FREE maxAbs = maxVal(Abs(ivec))
	Variable i, maxDiv=maxAbs[0]		// maxDiv is max possible divisor
	Variable tol = 0.1/maxDiv
	for (i=maxDiv;i>=2;i-=1)			// check all divisorts in range [2, maxDiv]
		MatrixOP/FREE/O mods = sum(abs(round(ivec/i)-(ivec/i)))
		if (mods[0]<tol)					// i is a factor of every elemnt in ivec
			ivec /= i
		endif
	endfor
	vec = ivec
End

//  ================================= End of Test Data =================================  //
//  ====================================================================================  //
