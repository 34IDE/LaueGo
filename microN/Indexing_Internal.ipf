#pragma rtGlobals=3		// Use modern globala access method and strict wave access.
#pragma ModuleName=IndexingInternal
#pragma version = 0.24
#include "IndexingN", version>=4.80

#if defined(ZONE_TESTING) || defined(QS_TESTING) || defined(ZONE_QS_TESTING)
#include "Indexing_InternalQZ_Gizmo"
#endif

Static Constant hc = 1.239841857			// keV-nm
Constant threshDivide = 3
Static Constant INDEXING_KEV_MAX = 30		// maximum for keV calc (not test) (kev)
Static Constant INDEXING_KEV_MIN = 2		// minimum energy (kev)


//	#define ZONE_TESTING
//	#define QS_TESTING
//	#define ZONE_QS_TESTING
//	#if defined(ZONE_TESTING) || defined(QS_TESTING) || defined(ZONE_QS_TESTING)
//	Constant DEBUG_LEVEL=3
//	#endif


Menu LaueGoMainMenuName
	MenuItemIfWaveClassExists("  Calc & Show Zones...","FittedPeakList*","DIMS:2,MINCOLS:11"),PutZoneLinesOnGraph($"")
End

Menu "Zones"
	MenuItemIfWaveClassExists("Calc & Show Zones...","FittedPeakList*","DIMS:2,MINCOLS:11"),PutZoneLinesOnGraph($"")
End
#if defined(ZONE_TESTING) || defined(QS_TESTING) || defined(ZONE_QS_TESTING)
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
//  ============================ Start of Zones+Qs Indexing ============================  //

Function/WAVE runIndexingQsZones(args)
	String args

	Wave FullPeakList = $StringByKey("FullPeakList",args)
	Variable keVmaxCalc = NumberByKey("keVmaxCalc",args)	// 14, maximum energy for peaks to check (keV)
	Variable keVmaxTest = NumberByKey("keVmaxTest",args)	// 30, maximum energy for finding peaks (keV)
	Variable angTol = NumberByKey("angleTolerance",args)	// ~0.5¡, angular matrch in pair angles (degree), as this gets bigger it runs slower
	Variable angPrecision = NumberByKey("angPrecision",args)	// ~0.05¡, used to define zones, precision of peaks, detector, & calibration (degree)
		// angPrecision	how well we know measured spot directions, combination of peak searching and calibration (used to find zones)
		// angTol			how well we angles from lattice constants, lattice constants may be poorly known. (used to match angles between Q's)
	Variable hp = NumberByKey("hp",args)							// preferred hkl
	Variable kp = NumberByKey("kp",args)
	Variable lp = NumberByKey("lp",args)
	Variable cone = NumberByKey("cone",args)					// for possible Qhats & zone axes, range of allowed tilts from central hkl (degree)
	Variable maxSpots = NumberByKey("maxSpots",args)			// -n max num. of spots from FullPeakList to use, default is 250
	Wave FullPeakList1 = $StringByKey("FullPeakList1",args)
	Wave FullPeakList2 = $StringByKey("FullPeakList2",args)
	Variable minSpotsForZone = NumberByKey("minSpotsForZone",args)// need at least this many spots to define a zone
	Variable minNQsMatch = NumberByKey("minNQsMatch",args)// number of Q's in a zone that need to match before I will use it to calc a possible angle
	Variable printIt = NumberByKey("printIt",args)
	maxSpots = numtype(maxSpots) || maxSpots<3 ? 250 : maxSpots
	minSpotsForZone = numtype(minSpotsForZone) || minSpotsForZone<3 || minSpotsForZone>30 ? 5 : round(minSpotsForZone)
	minNQsMatch = numtype(minNQsMatch) || minNQsMatch<2 || minNQsMatch>30 ? 3 : round(minNQsMatch)
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
	angPrecision = numtype(angPrecision) || angPrecision<=0 ? 0.05 : angPrecision
	Variable err
	timeIncrement(fresh=1)
	Variable sec0=stopMSTimer(-2)*1e-6
	Make/N=3/D/FREE ki={0,0,1}										// use this for the incident beam direction

	Wave GhatsMeasured=FullPeakList2Ghats(maxSpots,FullPeakList,FullPeakList1=FullPeakList1,FullPeakList2=FullPeakList2)	// convert peaks to a Qlist+intens+dNum
	String noteMeasured=note(GhatsMeasured)
	Variable NG0=DimSize(GhatsMeasured,0)
	Make/N=(NG0,3)/FREE GhatsOnly=GhatsMeasured[p][q]
	MatrixOP/FREE qvec0 = sumCols(GhatsOnly)					// assume that average Ghat is Q-vector to center of detector
	// find the kf for every measurd spot
	MatrixOP/FREE kfs = NormalizeRows( (rowRepeat(ki,NG0) - 2*colRepeat((GhatsOnly x ki),3)*GhatsOnly) )	//		kf^ = ki^ - 2*(ki^ . q^)*q^



	STRUCT allZonesStruct measuredZones							// holds info about the measured zones
	//	In the data find the measured zone axes.
	// measuredZones is the structure containing all the information about the measured zones
	// measuredZones.zn[].axisHat[3] and measuredZones.zn[].g.xyz[3] are in BEAMLINE coordinates
	err = ZoneInfoFromMeasuredSpots(measuredZones,GhatsMeasured, minSpotsForZone, angPrecision)	// find the measured zone axes info
	if (err || measuredZones.Nz < 2)
#ifdef ZONE_QS_TESTING
		printf "Unable to find any measured zones, err=%g\r",err
#endif
		return $""
	endif
	Variable nZmeasured = measuredZones.Nz
#ifdef ZONE_QS_TESTING
	if (printIt)
		printf "*made GhatsMeasured with %d G^'s from the %d zones in measuredZones (require %d Q's to define a zone), %s\r",DimSize(GhatsMeasured,0),nZmeasured,minSpotsForZone,timeIncrement()
		if (DEBUG_LEVEL>=4)
			print ZonesStruct2Str(measuredZones,maxPrint=5)
			if (DEBUG_LEVEL>=5)
				Duplicate GhatsMeasured, GhatsMeasuredView
				GhatsMeasuredView = abs(GhatsMeasuredView)<1e-12 ? 0 : GhatsMeasuredView
			endif
		endif
	endif
	Make/N=(nZmeasured,3)/FREE ZonesWave
	ZonesWave = measuredZones.zn[p].axisHat[q]
	GizmoMakeWavesForZones(GhatsMeasured,ZonesWave,$"")
	Variable iz
	for (iz=0; iz<(measuredZones.Nz); iz+=1)
		MakeZoneLinesOnGraph(measuredZones)
	endfor
#endif

	//	Make a list of Possible zones (this uses hklPrefer, keVmax, cone). and calc the q^s for each zone.
	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))
		DoAlert 0, "no lattice structure found, did you forget to set it?"
		return $""
	endif
	Wave direct=directFrom_xtal(xtal)
	Wave recipBase=recipFrom_xtal(xtal)
	Make/N=3/D/FREE hklPrefer={hp,kp,lp}							// hkl near center of data
#ifdef ZONE_QS_TESTING
	if (printIt)
		printf "*finished setup, have the measured zone info, %s\r",timeIncrement()
		print " "
	endif
#endif
	Variable NmaxTestZones = round(10*STRUCT_ZONES_MAX/nZmeasured)
	NmaxTestZones = limit(NmaxTestZones,4,STRUCT_ZONES_MAX)	// maximum number of possible zones to make
	Variable NminTestZones = 3											// want at least 3
	STRUCT allZonesStruct possibleZones								// holds info about the possible zones
	// possibleZones is the structure containing all the information about the measured zones
	// possibleZones.zn[].axisHat[3] and possibleZones.zn[].g.xyz[3] are in BASE CRYSTAL coordinates

#ifndef ZONE_QS_TESTING
	MakePossibleZoneInfo(possibleZones,xtal,hklPrefer,qvec0,cone*PI/180,keVmaxCalc,NminTestZones,NmaxTestZones,minSpotsForZone)
	Variable nZpossible = possibleZones.Nz
#else
	MakePossibleZoneInfo(possibleZones,xtal,hklPrefer,qvec0,cone*PI/180,keVmaxCalc,NminTestZones,NmaxTestZones,minSpotsForZone,printIt=DEBUG_LEVEL>=4)
	Variable nZpossible = possibleZones.Nz
	if (printIt)
		printf "*made %d (out of %d) Possible Zones to test against,  %s\r",nZpossible,NmaxTestZones,timeIncrement()
		if (DEBUG_LEVEL==4)
			print "\r        zone 0 is (check the first reflection):\r"+OneZoneStruct2Str(possibleZones.zn[0],maxPrint=6)+"\r"
		elseif (DEBUG_LEVEL>=5)
			print ZonesStruct2Str(possibleZones,maxPrint=4)
		endif
	endif
	if (DEBUG_LEVEL>=3)		// make Q^s for display in a Gizmo
		Make/N=((nZpossible)*STRUCT_ZONESPOTS_MAX,3)/O PossibleQhats=NaN
		Make/N=((nZpossible)*STRUCT_ZONESPOTS_MAX,3)/I/O Possiblehkls=NaN
		Variable iq, Nqs, nNs
		for (iq=0,Nqs=0;iq < nZpossible; iq+=1)
			nNs = possibleZones.zn[iq].Ns
			PossibleQhats[Nqs,Nqs+nNs-1][] = possibleZones.zn[iq].g[p-Nqs].xyz[q]
			Possiblehkls[Nqs,Nqs+nNs-1][] = possibleZones.zn[iq].h[p-Nqs].hkl[q]
			Nqs += nNs
		endfor
	endif
	Redimension/N=(Nqs,-1) PossibleQhats, Possiblehkls
#endif


//	pick a possible zone (loop)
//	compare q^(hkl)s with each of the measured zone g^'s to find the rotAxis that works
//			note, this step is like an indexing step, but rotation axis direction is fixed.
//			try to use a correlation method (but sparse) where matching is set by angTol (see sparse1Dcorrelate below)
//		save all that are good enough (based on angTol & NumSpots) along with some stats, {NumSpots, sum(intens), rmsErr, rotAxis}.
//
// given the q^ waves for two zones, what are the rotations that connect them?
	Variable Nalloc=200, NallRots
	Make/N=(Nalloc,3)/D/FREE AllRotations=NaN
	Make/N=(Nalloc)/I/FREE AllMatches=0
	Make/N=(3)/D/FREE rotAxis
	Variable i, ip, im, nrots, angTolRad=angTol*PI/180
	Variable angTolRad2 = angTolRad*angTolRad
	for (im=0,NallRots=0; im<nZmeasured; im+=1)
		Make/N=(measuredZones.zn[im].Ns,3)/FREE qMeas
		qMeas = measuredZones.zn[im].g[p].xyz[q]
		for (ip=0; ip<nZpossible; ip+=1)
			Make/N=(possibleZones.zn[ip].Ns,3)/FREE qPossible
			qPossible = possibleZones.zn[ip].g[p].xyz[q]
			Wave axes = RotsBetweenZones(qPossible,qMeas, minNQsMatch, angTolRad)
			nrots = DimSize(axes,0)
			if (nrots<1)
				continue
			endif
			// append all rotations in axes[][0,3] to AllRotations[][3] & AllMatches[] if not a duplicate rotation
			for (i=0;i<nrots;i+=1)
				rotAxis = axes[i][p]							// check if this rotAxis is already in AllRotations
				MatrixOP/FREE diffs2 = sumRows(magSqr(AllRotations-RowRepeat(rotAxis,Nalloc)))
				if (WaveMin(diffs2) < angTolRad2)			// an existing rotAxis, note: WaveMin(diffs2)==NaN the first time
					WaveStats/Q/M=1 diffs2						// find location of smallest value in diffs2
					AllMatches[V_minloc] += axes[i][3]		//  and inrcrement the matches
				else
					if (NallRots>=Nalloc)
						Nalloc += 200
						Redimension/N=(Nalloc,-1) AllRotations
						Redimension/N=(Nalloc) AllMatches
					endif
					AllRotations[NallRots][] = rotAxis[q]	// save the new rotation
					AllMatches[NallRots] = axes[i][3]		//   and the number of matches
					NallRots += 1
				endif
			endfor										// end of loop i, over rotation axes for this zone pair
		endfor											// end of loop ip, over possible zones
	endfor												// end of loop im, over measured zones
	Redimension/N=(NallRots,-1) AllRotations// trim to correct length
	Redimension/N=(NallRots) AllMatches
	Variable maxMatches = WaveMax(AllMatches)
	if (NallRots<1 || !(maxMatches>1))
		printf "Failed to find any rotation axes that rotate %d spots from the Possible zones to the Measured zones, try increasing keVmaxCalc (=%g keV), or cone (=%g¡).\r", minNQsMatch,keVmaxCalc,cone
		return $""
	endif

#ifdef ZONE_QS_TESTING
MatrixOP/O AllRotationsAngleView = sqrt(sumRows(magsqr(AllRotations)))
AllRotationsAngleView *= 180/PI
printf "AllRotations[%g][3] has a max of %g, and the smallest rotation is %g¡,  %s\r",NallRots,WaveMax(AllMatches),WaveMin(AllRotationsAngleView),timeIncrement()
Duplicate/O AllRotations, AllRotationsView
for (i=0;i<NallRots;i+=1)
	if (AllMatches[i] >= maxMatches)
		rotAxis = AllRotations[i][p]
		printf "#%d,  rot = %s,   |rot|=%.3f¡   %d matches\r",i,vec2str(rotAxis),norm(rotAxis)*180/PI,AllMatches[i]
	endif
endfor
Duplicate/O AllMatches, AllMatchesView, AllMatchesViewSize
Redimension/N=(-1,3) AllMatchesViewSize
AllMatchesViewSize = limit(30*AllMatches[p]^4/(maxMatches^4),1,20)
#endif

	WaveStats/Q/M=1 AllMatches
	rotAxis = AllRotations[V_maxloc][p]
	Variable angle = norm(rotAxis)*180/PI
	Make/N=(3,3)/D/FREE rotMat
	rotationMatAboutAxis(rotAxis,angle,rotMat)	// make rotMat the matrix version of rotAxis
	MatrixOP/O/FREE recip = rotMat x recipBase
	//	printWave(recip,name="rotated recip",brief=1,zeroThresh=1e-9)
#ifdef ZONE_QS_TESTING
	if (printIt)
		printf "*found <rot> = %s  |<rot>|=%.3f¡,  with %d matches, %s\r",vec2str(rotAxis,zeroThresh=1e-9),angle,maxMatches,timeIncrement()
	endif
#endif

	// find rms error between Ghats and the predicted hkl directions
	STRUCT microGeometry geo
	if (FillGeometryStructDefault(geo,alert=1))	//fill the geometry structure with current values
		return $""
	endif

	Make/N=3/D/FREE vec
	Variable startx=NumberByKey("startx",noteMeasured,"="), groupx=NumberByKey("groupx",noteMeasured,"=")
	Variable starty=NumberByKey("starty",noteMeasured,"="), groupy=NumberByKey("groupy",noteMeasured,"=")
	Variable/C pz
	Variable NG1, delta
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
#ifdef ZONE_QS_TESTING
	if (printIt)
		printf "*found %g hkl's that match this rotation, %s\r",NG1,timeIncrement()
	endif
#endif
	if (NG1<3)										// found less than 3 matching hkl's, failure
		return $""
	endif
	Redimension/N=(NG1,-1) hkls

	// do a least-squares optimization on the rotation vector
	err = RefineOrientation(GhatsMeasured,hkls,recipBase,rotAxis)
	angle = norm(rotAxis)*180/PI
	rotationMatAboutAxis(rotAxis,angle,rotMat)		// re-make rotMat
	MatrixOP/O/FREE recip = rotMat x recipBase		//  and re-make recip
#ifdef ZONE_QS_TESTING
	if (printIt)
		printf "*after non-linear least-squares optimization, <rot> = %s  |<rot>|=%.3f¡,  %s\r",vec2str(rotAxis,zeroThresh=1e-9),angle,timeIncrement()
	endif
#endif

	// now save all indexed values using the optimized rotation
	Variable Npatterns=1, ipat=0, NG
	String peakListName = ParseFilePath(3,StringByKey("peakListWave",noteMeasured,"="),":",0,0)
	String IndexedName = CleanupName("FullPeakIndexed"+ReplaceString("FullPeakList",peakListName,""),0)
	Make/N=(NG0,12,Npatterns)/O $IndexedName/WAVE=IndexedWave = NaN
	Variable dNum, px,py, keV, intensity, intensityMax=-Inf
	Variable rmsGhat, minError=Inf, maxError=-Inf
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
			intensityMax = numtype(intensity) ? intensityMax : max(intensityMax,intensity)
			pz = q2pixel(geo.d[dNum],qhat)
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
	Variable goodness0 = goodNessOfPattern(IndexedWave,ipat)
	rmsGhat = sqrt(rmsGhat/NG)
	IndexedWave[][8][ipat] = abs(IndexedWave[p][8][ipat])<5e-6 ? 0 : IndexedWave[p][8][ipat]	// precision limit
	Variable executionTime = stopMSTimer(-2)*1e-6 - sec0



#ifdef ZONE_QS_TESTING
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
	if (strlen(xtal.SpaceGroupID))
		wnote = ReplaceStringByKey("SpaceGroupID",wnote,xtal.SpaceGroupID,"=")
	endif
	if (numtype(xtal.SpaceGroupIDnum))
		wnote = ReplaceNumberByKey("SpaceGroupIDnum",wnote,xtal.SpaceGroupIDnum,"=")
	endif
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
Static Function/WAVE MakeZoneLinesOnGraph(allZones)
	Struct allZonesStruct &allZones
	STRUCT microGeometry g
	FillGeometryStructDefault(g)			//fill the geometry structure with current values

//print ZonesStruct2Str(allZones)

	Make/N=3/D/FREE qvec
	Variable/C pz
	Variable j, i, Nspots, Npnts=0, ipnt=0
	for (j=0;j<allZones.Nz;j+=1)
		Npnts += allZones.zn[j].Ns + 1
	endfor
	Make/N=(Npnts,2)/O ZoneLines=NaN

	for (j=0,ipnt=0; j<allZones.Nz; j+=1)
		Nspots = allZones.zn[j].Ns
		if (Nspots<2)
			continue
		endif
		Make/N=(Nspots,2)/FREE xLine=NaN, yLine=NaN
		for (i=0; i<Nspots; i+=1)
			qvec = allZones.zn[j].g[i].xyz[p]
			pz = q2pixel(g.d[0],qvec)			// convert qvector to pixel position
			xLine[i] = real(pz)
			yLine[i] = imag(pz)
		endfor
		Sort xLine, xLine, yLine
		ZoneLines[ipnt,ipnt+Nspots-1][0] = xLine[p-ipnt]
		ZoneLines[ipnt,ipnt+Nspots-1][1] = yLine[p-ipnt]
		ipnt += Nspots+1
	endfor

	Redimension/N=(Npnts-1,-1) ZoneLines		// remove final NaN
	return ZoneLines
End


Static Function ZoneInfoFromMeasuredSpots(zInfo,Ghats,Nmin,precision)
	// returns zoneInfo that contains zone axis and the list of Ghats contributing to each zone
	// zInfo is the structure containing all the information about the measured zones
	// zInfo.zn[].axisHat[3] and zInfo.zn[].g.xyz[3] are in BEAMLINE (measured) coordinates
	STRUCT allZonesStruct &zInfo	// holds info about the measured zones
	Wave Ghats							// list of measured G^ (first 3 columns contain G^, the rest are ignored)
	Variable Nmin						// minimum number of spots to constitute a valid measured zone
	Variable precision				// angular tolerance used in accepting Q-vec as part of a zone (degree)

	if (!WaveExists(Ghats))
		return 1
	endif
	Variable NG=DimSize(Ghats,0)
	if (NG<Nmin || Nmin<=2)
		return 1
	endif
	if (strsearch(StringByKey("peakListWave",note(Ghats),"="),"FullPeakListTest",0)>0 && WaveInClass(Ghats,"QsMeasuredTest"))
		Wave GhatsMeasuredHKL = $(GetWavesDataFolder(Ghats,1)+"GhatsMeasuredHKL")
	endif
	Make/N=(NG,3)/D/FREE gs = Ghats[p][q]					// the normalized g-vectors
	Make/N=(STRUCT_ZONES_MAX)/I/FREE nQsinZone=0		// holds number of spots contributing to a zone
	Make/N=(STRUCT_ZONES_MAX,3)/D/FREE ZoneAxes=0		// direction of each zone axis
	Make/N=(STRUCT_ZONES_MAX,3)/D/FREE AllAxes=NaN	// holds local info about the measured zones
	Make/N=(STRUCT_ZONES_MAX,STRUCT_ZONESPOTS_MAX)/I/FREE zSpotIndex=-1
	Variable minDot = cos(5.0*PI/180)			// G-vectors must be > 5 degrees apart to define a zone (not too parallel)
	Variable tol = cos(precision*PI/180)		// tolerance on dot product
	Variable perpTol = cos((90-precision)*PI/180)
	Make/N=3/D/FREE g0, g1, axis
	Variable Nz=0								// Nz is number of zones found, m is first spot to start searching
	Variable nQs								// number of spots in this zone
	Variable i,j,m								// m is first spot to start searching
	for (m=0,Nz=0; m<(NG-1) && Nz<STRUCT_ZONES_MAX; m+=1)	// m is first spot that defines this possible zone
		g0 = gs[m][p]
		for (j=m+1; j<(NG-Nmin+1) && Nz<STRUCT_ZONES_MAX; j+=1)// find second vector to define the zone
			g1 = gs[j][p]
			if (abs(MatrixDot(g0,g1)) > minDot)
				continue							// vectors too parallel for defining a zone axis, keep searching
			endif

			// have a possible pair, (g0,g1) from (m,j)
			Cross g0,g1
			Wave W_Cross=W_Cross
			axis = W_Cross						// possible new zone axis
			normalize(axis)
			if (!isNewDirection(ZoneAxes,axis,tol,anti=1))
				continue							// already have this axis, try again
			endif
			// the (g0,g1) pair defining axis[3], has not been found before, check it out

			// have a new unique zone axis, find the spots belonging to this zone
			// find the best axis comparing to all gs that are close to perpendicualr to axis
			MatrixOP/FREE isPerp = greater(perpTol, Abs(gs x axis))	// true when g^ perp to axis
			nQs = sum(isPerp)
			if (nQs < Nmin)					// not enough spost to qualify as a valid zone, skip
				continue
			endif

			Variable ip
			Make/N=(nQs,3)/D/FREE gPlane
			for (i=0,ip=0; i<NG; i+=1)
				if (isPerp[i])
					gPlane[ip][] = gs[i][q]
					ip += 1
				endif
			endfor
			Wave axis = FindVecPerp2Many(gPlane)	// set axis to vec that is most perp to all of gPlane
			if (!WaveExists(axis))
				continue
			endif

			ZoneAxes[Nz][] = axis[q]		// axis of this zone
			nQs = min(nQs,STRUCT_ZONESPOTS_MAX)	// only room for STRUCT_ZONESPOTS_MAX spots
			nQsinZone[Nz] = nQs				// save number of spots contributing to this zone
			for (i=0,ip=0; i<NG; i+=1)
				if (isPerp[i])
					zSpotIndex[Nz][ip] = i
					ip += 1
				endif
			endfor
			Nz += 1

//			axisAvg = 0
//			oneSpotIndexList = -1
//			for (i=0,nQs=0; i<NG; i+=1)	// find how many other spots belong to this axis
//				g1 = Ghats[i][p]
//				Cross g0,g1
//				mag = normalize(W_cross)
//				dot = MatrixDot(W_cross,axis)
//				if (mag<sinMax)				// only look at the cross product when g0 & g1 not too parallel
//					oneSpotIndexList[nQs] = i
//					axisAvg += axis
//					nQs += 1						// this spot real close to g0, increment nQs
//				elseif (abs(dot) > tol)
//					oneSpotIndexList[nQs] = i
//					axisAvg += dot<0 ? -W_cross : W_cross
//					nQs += 1						// a spot in zone, increment nQs
//				endif
//			endfor								// for i
//			if (nQs >= Nmin)					// found another zone with at least Nmin spots, save it
//				normalize(axisAvg)			// want a unit vector, so don't bother dividing by nQs
//				ZoneAxes[Nz][]=axisAvg[q]	// axis of this zone
//				nQsinZone[Nz] = nQs			// save number of spots contributing to this zone
//				zSpotIndex[Nz][0,nQs-1] = oneSpotIndexList[q]
//				Nz += 1
//			endif

		endfor									// for j
	endfor										// for m
	KillWaves/Z W_Cross
	if (Nz<1)
		return 1									// failed to find anything
	endif

	Redimension/N=(Nz) nQsinZone
	Duplicate/FREE nQsinZone, indexWave
	MakeIndex/R nQsinZone, indexWave	// sort zones by decreasing number of nQsinZone

	ZonesStructClear(zInfo)				// put sorted values into the zInfo structure
	zInfo.Nz = Nz
	zInfo.NsMin = WaveMin(nQsinZone)	// smallest number of sponts in a zone
	zInfo.measured = 1
	for (m=0;m<Nz;m+=1)
		j = indexWave[m]
		zInfo.zn[m].Ns = nQsinZone[j]
		zInfo.zn[m].axisHat[0] = ZoneAxes[j][0]
		zInfo.zn[m].axisHat[1] = ZoneAxes[j][1]
		zInfo.zn[m].axisHat[2] = ZoneAxes[j][2]
		for (i=0;i<nQsinZone[j];i+=1)
			zInfo.zn[m].g[i].xyz[0] = Ghats[zSpotIndex[j][i]][0]
			zInfo.zn[m].g[i].xyz[1] = Ghats[zSpotIndex[j][i]][1]
			zInfo.zn[m].g[i].xyz[2] = Ghats[zSpotIndex[j][i]][2]
//			if (WaveExists(GhatsMeasuredHKL))
//				zInfo.zn[m].h[i].hkl[0] = GhatsMeasuredHKL[zSpotIndex[j][0]]
//				zInfo.zn[m].h[i].hkl[1] = GhatsMeasuredHKL[zSpotIndex[j][1]]
//				zInfo.zn[m].h[i].hkl[2] = GhatsMeasuredHKL[zSpotIndex[j][2]]
//			endif
		endfor
	endfor
	zInfo.precision = precision
	zInfo.GhatSource = GetWavesDataFolder(Ghats,2)
	return 0
End
//
Static Function/WAVE FindVecPerp2Many(hats)		// return a vec that is most perp to all of hats
	Wave hats					// (N,3)
	Variable N=DimSize(hats,0)
	if (N<2)
		return $""
	endif

	Make/N=3/D/FREE g0 = hats[0][p]

	MatrixOp/FREE absDots = Abs(hats x g0)
	WaveStats/M=1/Q absDots
	Make/N=3/D/FREE g1 = hats[V_minloc][p]	// g1 is most perpendicular to g0
	Cross g0,g1
	Wave W_Cross=W_Cross
	Make/N=3/D/FREE axis=W_Cross
	normalize(axis)								// approximate answer
	if (N<3)
		return axis									// only two vectors, so this is the answer
	endif

	Make/N=3/D/FREE a1,a0=g0
	Cross axis, a0
	a1 = W_Cross
	KillWaves/Z W_Cross
	normalize(a1)
	// now axis is close to answer and a0 & a1 are perp to axis

	// want to minimize SUM{ (ax*a0 + ay*a1 + axis) . hats }
	// so take derivatives and setting to zero gives:
	//	2*ax*g0i^2 + 2*g0i*ay*g1i + 2*g0i*gai = 0
	//	2*ay*g1i^2 + 2*ax*g0i*g1i + 2*g1i*gai = 0
	//
	//	g0i^2   * ax + g0i*g1i * ay = -g0i*gai
	//	g0i*g1i * ax + g1i^2 * ay   = -g1i*gai
	//
	//	g0i^2,   g0i*g1i   x   ax   = -g0i*gai
	// g0i*g1i, g1i^2         ay     -g1i*gai
	//        mat         x  axay  =  bvec
	MatrixOP/FREE g0i = hats x a0			// hats[i] . a0
	MatrixOP/FREE g1i = hats x a1			// hats[i] . a1
	MatrixOP/FREE gai = hats x axis			// hats[i] . axis

	Make/N=(2,2)/D/FREE mat	// form linear equation mat x axay = bvec
	Make/N=(2)/D/FREE bvec
	MatrixOP/FREE ss = sum(magSqr(g0i))
	mat[0][0] = ss[0]				// SUM{ g0i^2 }
	MatrixOP/FREE ss = sum(g0i*g1i)
	mat[0][1] = ss[0]				// SUM{ g0i*g1i }
	mat[1][0] = ss[0]				// SUM{ g0i*g1i }
	MatrixOP/FREE ss = sum(magSqr(g1i))
	mat[1][1] = ss[0]				// SUM{ g1i^2 }

	MatrixOP/FREE ss = -sum(g0i * gai)
	bvec[0] = ss[0]	
	MatrixOP/FREE ss = -sum(g1i * gai)
	bvec[1] = ss[0]

	// now have mat x axay = bvec, solve for axay
	MatrixOP/FREE axay = Inv(mat) x bvec
	axis = axay[0]*a0[p] + axay[1]*a1[p] + axis[p]
	normalize(axis)				// final answer
	return axis
End


Static Function MakePossibleZoneInfo(zInfo,xtal,hklC,qCin,cone,keVmax,NzMin,NzMax,minSpotsForZone,[maxLatticeLen,precision,printIt])
	// returns zoneInfo that contains zone axis and the list of Ghats & hkls contributing to each zone
	// zInfo is the structure containing all the information about the Possible zones
	// zInfo.zn[].axisHat[3] and zInfo.zn[].g.xyz[3] are in BASE CRYSTAL coordinates
	STRUCT allZonesStruct &zInfo	// holds info about the possible zones
	STRUCT crystalStructure &xtal
	Wave qCin						// direction of q-vector that diffracts to detector center (in Beam Line coords)
	Wave hklC						// hkl that diffracts to detector center
	Variable cone					// limits selection of zone axes and kf, only zones least 90-cone from qC (radian)
	Variable keVmax				// max energy (keV)
	Variable NzMin					// minimum number of zones to accept
	Variable NzMax					// max number of possible zone axes
	Variable minSpotsForZone	// minimum number of spots needed to make a zone acceptable
	Variable maxLatticeLen		// max len of a (lattice x ijk) vector
	Variable precision			// angular tolerance used in accepting Q-vec as part of a zone (radian)
	Variable printIt
	maxLatticeLen = ParamIsDefault(maxLatticeLen) || numtype(maxLatticeLen) || maxLatticeLen<=0 ? NaN : maxLatticeLen
	precision = ParamIsDefault(precision) || numtype(precision) || precision<=0 ? 1e-6 : precision*PI/180
	printIt = ParamIsDefault(printIt) || numtype(printIt) || strlen(GetRTStackInfo(2))==0 ? 0 : !(!printIt)
	cone = limit(cone,0.001,PI-precision)		// cone in range [0.001, ¹-tol] (rad)
	Variable coneZone = limit(cone,0.001,PI/2-precision)	// zone cone in range [0.001, ¹/2-tol] (rad)
	NzMax = limit(NzMax,1,STRUCT_ZONES_MAX)
	NzMin = min(NzMin,NzMax)
	Wave lattice=directFrom_xtal(xtal)
	maxLatticeLen = numtype(maxLatticeLen) ? ( MatrixDet(lattice)^(1/3) ) * 6 : maxLatticeLen
	MatrixOP/FREE qC = Normalize( qCin )		// this needs to be normalized, but don't change passed qCin
	MatrixOP/FREE qC0 = Normalize( 2*PI * (Inv(lattice))^t x hklC )	// similar to qC, but in UN-rotated coordinates

	MatrixOP/FREE ijkMax = floor( rowRepeat(maxLatticeLen,3)/sqrt(sumCols(magSqr(lattice)))^t )
	Variable imax=ijkMax[0], jmax=ijkMax[1], kmax=ijkMax[2], 	maxMax=WaveMax(ijkMax)
	if (printIt)
		printf "  find directions of possible zones using indicies in: i=[%d, %d], j=[%d, %d], k=[%d, %d]\r",-imax,imax, -jmax,jmax, -kmax,kmax
	endif

timeIncrement(fresh=1)
	Variable maxDot = 1-cos(coneZone)			// maximum allowed value of dot(qC,zoneAxis)
	Variable cosTol = cos(precision)			// convert angle to dot, this should be slightly less than 1.0
	Make/N=(NzMax,3)/I/FREE TestAxesijks=0	// contains possible i,j,k for zone axes
	Make/N=(NzMax,3)/D/FREE TestAxesHats=NaN// contains corresponding possible unit vectors, x,y,z

	Wave ijks = hklOrderedList(maxMax)
	Variable Nijks=DimSize(ijks,0)
#ifdef ZONE_QS_TESTING
printf "  -- MakePossibleZoneInfo(), starting with a list of %d ijk's to examine, %s\r",Nijks,timeIncrement()
#endif
	Make/N=3/I/FREE ijk
	Variable m, Nz
	for (m=0,Nz=0; m<Nijks && Nz<NzMax; m+=1)
		ijk = ijks[m][p]
		MatrixOP/FREE/O vec = lattice x ijk
		if (normalize(vec)>maxLatticeLen)			// length is too long, skip
			continue
		elseif (abs(MatrixDot(vec,qC0))>maxDot)	// zone axis too close to qC
			continue
		elseif (!isNewDirection(TestAxesHats,vec,cosTol,anti=1))	// skip already existing directions too
			continue
		endif
		TestAxesHats[Nz][] = vec[q]					// save this vec as a new possible zone axis
		TestAxesijks[Nz][] = ijk[q]					//  and save the ijk that goes with it
		Nz += 1
	endfor
	Redimension/N=(Nz,-1) TestAxesHats, TestAxesijks
#ifdef ZONE_QS_TESTING
printf "  -- MakePossibleZoneInfo(), found %d possible zone axes, %s\r",Nz,timeIncrement()
#endif
	if (Nz<NzMin)											// did not find enough zones
		return 1
	endif
#ifdef ZONE_QS_TESTING
Duplicate/O TestAxesijks, TestAxesijksView
Duplicate/O TestAxesHats, TestAxesHatsView
TestAxesHatsView = abs(TestAxesHatsView)<1e-12 ? 0 : TestAxesHatsView
#endif

	// next fill the structure info for each zone axis
	Variable NsMin, Ns
	zInfo.measured = 0									// 0=calculated (not measured)
	zInfo.precision = precision
	zInfo.GhatSource = ""
	zInfo.desc = ""

	Make/N=(xtal.N)/WAVE/FREE atomWaves=$("root:Packages:Lattices:atom"+num2istr(p))	// cannot reach these waves from within a separate thread
	Make/N=(Nz)/FREE/WAVE AllPossibleZonesW
	MultiThread AllPossibleZonesW[] = PossibleQsInOneZone(p,TestAxesijks,lattice,xtal,qC,hklC,cone,keVmax,atomWaves=atomWaves,printIt=printIt)
#ifdef ZONE_QS_TESTING
print "  -- MakePossibleZoneInfo(), just after calls to PossibleQsInOneZone(...), ",timeIncrement()
#endif

	Variable i, j, NtotalQs=0
	for (i=0,m=0,NsMin=Inf; i<Nz; i+=1)			// for each zone
		ijk = TestAxesijks[i][p]
		vec = TestAxesHats[i][p]
		Wave possibleZonesW = AllPossibleZonesW[i]	// contains result from PossibleQsInOneZone()
		// columns of possibleZonesW are {Qx,Qy,Qz, H,K,L, angle}
		Ns = DimSize(possibleZonesW,0)
		if (Ns<minSpotsForZone)
			continue
		endif

		zInfo.zn[m].Ns = Ns
		NsMin = min(NsMin,Ns)
		zInfo.zn[m].axisHat[0] = vec[0]
		zInfo.zn[m].axisHat[1] = vec[1]
		zInfo.zn[m].axisHat[2] = vec[2]
		zInfo.zn[m].axisijk[0] = ijk[0]				// have ijk, because these are possible (not measured)
		zInfo.zn[m].axisijk[1] = ijk[1]
		zInfo.zn[m].axisijk[2] = ijk[2]
		for (j=0;j<Ns;j+=1)
			zInfo.zn[m].g[j].xyz[0] = possibleZonesW[j][0]	// copy direction and hkl of each spot
			zInfo.zn[m].g[j].xyz[1] = possibleZonesW[j][1]
			zInfo.zn[m].g[j].xyz[2] = possibleZonesW[j][2]
			zInfo.zn[m].h[j].hkl[0] = possibleZonesW[j][3]	// have hkl, because these are possible (not measured)
			zInfo.zn[m].h[j].hkl[1] = possibleZonesW[j][4]
			zInfo.zn[m].h[j].hkl[2] = possibleZonesW[j][5]
			NtotalQs += 1
		endfor
		m += 1
	endfor
	zInfo.Nz = m
	zInfo.NsMin = NsMin
#ifdef ZONE_QS_TESTING
printf "  -- MakePossibleZoneInfo(), at end found %d zone axes, containing %d Q's (require %d spots for a zone),  %s\r",m,NtotalQs,minSpotsForZone,timeIncrement()
#endif
	return 0
End


ThreadSafe Static Function/WAVE PossibleQsInOneZone(pnt,ijks,lattice,xtal,qhat0in,hklC,cone,keVmax,[atomWaves,kin,printIt])
	// returns hkl & G^ for all hkl found on this zone, zone axis is lattice x ijk
	Variable pnt					// selects one ijk from ijks (done this way for MultiThreading)
	Wave ijks						// i,j,k for this zone,   zone axis = lattice x ijk
	Wave lattice					// direct lattice
	STRUCT crystalStructure &xtal
	Wave qhat0in					// direction of q-vector that diffracts to center of detector
	Wave hklC						// hkl that diffracts to near detector center (the preferred hkl)
	Variable cone					// max cone angle from kf (radian)
	Variable keVmax				// max energy (keV)
	Wave/WAVE atomWaves			// contains wave refs to root:Packages:Lattices:atomN
	Wave kin							// optional incident wave direction (usually defaults to 001)
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? 0 : !(!printIt)

	Make/N=3/I/FREE ijk = ijks[pnt][p]
	if (printIt)
		printf "PossibleQsInOneZone(ijk=[%s], direct, q0^=%s, hklC=(%s), %g¡, %gkeV",vec2str(ijk,sep=",",bare=1),vec2str(qhat0in,sep=",",zeroThresh=1e-7),vec2str(hklC,sep=",",bare=1),cone*180/PI,keVmax
		if (!ParamIsDefault(kin))
			printf ", kin=",vec2str(kin,sep=",")
		endif
		printf ")\r"
	endif

	Variable tol = 0.01									// angular tolerance (radian) (= ~0.5¡)
	Make/N=3/D/FREE qhat0=qhat0in[p]				// q-vector that diffracts to center of detector
	normalize(qhat0)

	Make/N=3/D/FREE ki={0,0,1}						// use this for the incident beam direction
	if (!ParamIsDefault(kin))							// one was passed, use it
		ki = kin[p]
		normalize(ki)
	endif
	MatrixOP/FREE axisBase = Normalize( lattice x ijk )	// zone axis
	MatrixOP/FREE recipBase = 2*PI * (Inv(lattice))^t	// reciprocal lattice

	MatrixOP/FREE qC = Normalize(recipBase x hklC)
	Cross qC, qhat0									// find a rotation that puts (recipBase x hkl0) parallel to qhat0
	Wave W_Cross=W_Cross
	Variable sine=normalize(W_Cross), cosine=MatrixDot(qC,qhat0)
	Variable angle = atan2(sine,cosine)*180/PI
	Make/N=(3,3)/D/FREE rot0
	rotationMatAboutAxis(W_Cross,angle,rot0)
	MatrixOP/FREE rot0Inv = Inv(rot0)				// I will want this later
	if (printIt)
		printf "  rotating %g¡ about the %s\r",angle,vec2str(W_Cross,zeroThresh=1e-7)
	endif
	KillWaves/Z W_Cross
	MatrixOP/FREE recip = rot0 x recipBase
	MatrixOP/FREE axis = rot0 x axisBase
	//		so now:  qhat0 || (rot0 x recip0 x hkl0)

	Variable dot = MatrixDot(ki,qhat0)				// points to detector center
	Make/N=3/D/FREE kf0 = (ki - 2*dot*qhat0)	//		kf^ = ki^ - 2*(ki^ . q^)*q^

	Variable maxTheta = limit(acos(MatrixDot(kf0,ki))+cone,tol,PI-tol)/2
	Variable maxQLen = 4*PI*sin(maxTheta) * keVmax/hc	// max len of a q-vector = recip0 x hkl)
	MatrixOP/FREE hklMax = floor( rowRepeat(maxQLen,3)/sqrt(sumCols(magSqr(recip)))^t )

	Variable maxMax=WaveMax(hklMax)
	if (printIt)
		printf "  find Possible Q-directions using index range: h=[%+d, %+d], k=[%+d, %+d], l=[%+d, %+d], max(theta)=%g¡\r",-hklMax[0],hklMax[0], -hklMax[1],hklMax[1], -hklMax[2],hklMax[2],maxTheta*180/PI
	endif

	//		knowing lattice=L and ijk, want hkl s.t.  [ L x ijk ] ¥ [ Inv(L)^t x hkl ] = 0
	//			or  axis ¥ [ recip x hkl ] = 0
	//			for a given k & l, this can be solved for h (thus removing one loop from below)
	//
	//			since:  recip x hkl = recip x [ (h00)+(0kl) ] = recip x (h00) + recip x (0kl) ]
	//			axis ¥ [ recip x h00 ] = - axis ¥ [ recip x 0kl]
	//			h (axis ¥ a*) = - axis ¥ [ recip x 0kl]
	//			h  = -( axis ¥ [ recip x 0kl] ) / (axis ¥ a*)

	Variable minDot=cos(cone)							// minimum allowed value of dot(kf0,kf)
	Variable cosTol=cos(tol)
	Variable coneMax=-Inf, angleMax=-Inf

	// find which of the two {h,k,l} to iterate over, choose the two with smallest
	// of |axis¥a*|, |axis¥b*|, |axis¥c*|
	MatrixOP/FREE absDots = abs(axis^t x recip)^t	// absDots[] contains {|axis¥a*|, |axis¥b*|, |axis¥c*|}, used to make index
	Make/N=3/I/FREE index
	MakeIndex absDots, index
	if (printIt)
		print "  {|axis¥a*|, |axis¥b*|, |axis¥c*|} =",vec2str(absDots,zeroThresh=1e-12)
		printf "  index = %s  solve for hkl[index[2]]=hkl[%d],  loop over hkl[%d] & hkl[%d] (the first two)\r",vec2str(index),index[2],index[0],index[1]
	endif

	Variable Nmax=STRUCT_ZONESPOTS_MAX
	Make/N=3/D/FREE kf, q0=0, star=recip[p][index[2]]		// star is one of a*, b*, or c*
	Make/N=3/I/FREE hkl
	Variable axis_dot_star = MatrixDot(axis,star)
	axis_dot_star = abs(axis_dot_star)<1e-7 ? 0 : axis_dot_star
	//	print "axis_dot_star = ",axis_dot_star
	Make/N=3/I/FREE hklLo, hklHi
	hklLo[index[0]] = -hklMax[index[0]]	;		hklHi[index[0]] = hklMax[index[0]]
	hklLo[index[1]] = -hklMax[index[1]]	;		hklHi[index[1]] = hklMax[index[1]]
	hklLo[index[2]] = 0							;		hklHi[index[2]] = 0
	Wave hkls = hklOrderedList(maxMax, hklLo=hklLo, hklHi=hklHi)

	Make/N=(Nmax,7)/D/FREE out=NaN
	Make/N=(Nmax,3)/D/FREE TestQs=NaN
	Variable Nhkls=DimSize(hkls,0)
	Variable m, Ns, Qmag, keV

	for (m=0,Ns=0; m<Nhkls && Ns<Nmax; m+=1)
		hkl[index[0]] = hkls[m][index[0]]			// hkls[m][]
		hkl[index[1]] = hkls[m][index[1]]
		// solve for hkl2 = hkl[index[2]]
		// do not loop over hkl2, rather solve for hkl2 using:  axis ¥ [ recip x hkl ] = 0
		//	hkl2  = -( axis ¥ [ recip x 0kl] ) / (axis ¥ star)   from above
		hkl[index[2]] = 0
		MatrixOP/FREE/O axis_dot_q0kl = axis . (recip x hkl)
		hkl[index[2]] = round( -axis_dot_q0kl[0] / axis_dot_star )
		MatrixOP/FREE/O qhat = recip x hkl
		Qmag = normalize(qhat)
		dot = MatrixDot(ki,qhat)			// this will be negative for legit reflections
		kf = (ki - 2*dot*qhat)				// kf^ = ki^ - 2*(ki^ . q^)*q^
		keV = -Qmag*hc / (4*PI*dot)

		if (abs(MatrixDot(axis,qhat)) > 1e-5)	// axis ¥ [ recip x hkl ] = 0, not on zone
			continue
		elseif (Qmag>maxQLen)				// length is too long, skip
			continue
		elseif (MatrixDot(kf,kf0)<minDot)
			continue								// kf is too far from kf0
		elseif (cone<PI/2 && dot >= 0)	// for cone angles < 90¡, also check for valid reflection
			continue								// (ki .dot. q^) must be negative!, Bragg angle must be < 90¡
		elseif (cone<PI/3 && keV>keVmax)// for cone angles < 60¡, also check for valid energy
			continue
		elseif (!isNewDirection(TestQs,qhat,cosTol))
			continue								// skip already existing directions too
		elseif (!allowedHKL(hkl[0],hkl[1],hkl[2],xtal,atomWaves=atomWaves))
			continue
		endif
		TestQs[Ns][] = qhat[q]				// a new valid & unique q-direction

		MatrixOP/FREE/O qhatBase = Normalize( recipBase x hkl )
		out[Ns][0,2] = qhatBase[q]		// q^ in crystal base coordinates
		out[Ns][3,5] = hkl[q-3]
		if (Ns==0)
			q0 = qhatBase[p]					// q of first reflection in crystal base coordinates
			out[Ns][6] = 0						// angle from first, measured around the zone axis (radian)
		else
			out[Ns][6] = acos(MatrixDot(q0,qhatBase))		// angle about the axis, measured from first point
		endif
		coneMax = max(coneMax,acos(limit(MatrixDot(kf,kf0),-1,1)))
		angleMax = max(angleMax,acos(limit(MatrixDot(q0,qhatBase),-1,1)))
		Ns += 1
	endfor					// end of loop m
	Redimension/N=(Ns,-1) out
	SetDimLabel 1,0,Qx,out	;	SetDimLabel 1,1,Qy,out	;	SetDimLabel 1,2,Qz,out
	SetDimLabel 1,3,H,out	;	SetDimLabel 1,4,K,out	;	SetDimLabel 1,5,L,out
	SetDimLabel 1,6,angle,out
	if (printIt)
		printf "  for zone axis=%s, found %d spots,   max cone angle = %.2f¡ out of [0,%g¡],  angular range around zone is %.2f¡\r",vec2str(axis,zeroThresh=1e-7,sep=", "),Ns,coneMax*180/PI,cone*180/PI,angleMax*180/PI
	endif
	return out
End
//
//Function test_PossibleQsInOneZone(flag)
//	Variable flag
//
//	STRUCT crystalStructure xtal
//	FillCrystalStructDefault(xtal)
//	Wave direct=directFrom_xtal(xtal)
//
//	Make/N=3/D/FREE qC = {0,1,-1}
//	Make/N=3/D/FREE hklC={0,0,1}						// hkl that diffracts to detector center
//	Variable keVmax = 17
//	Variable cone = 30.*PI/180
//
//	Make/N=(2,3)/D/FREE ijks
//	Make/N=3/I/FREE ijk
//	ijks[0][0]= {1,0}
//	ijks[0][1]= {0,1}
//	ijks[0][2]= {0,0}
//
//	timeIncrement(fresh=1)
//	Variable ip, ip2
//	for (ip=0;ip<10;ip+=1)
//		ip2 = 2^ip
//		if (flag & ip2)
//			Wave out = PossibleQsInOneZone(ip,ijks,direct,xtal,qC,hklC,cone,keVmax,printIt=1)
//			printf "    test took %s\r",timeIncrement()
//			print " "
//		endif
//	endfor
//	printf "End of All Tests,  %s\r",timeIncrement()
//
//	if (WaveExists(out))
//		Duplicate/O out, outView
//		if (DimSize(out,0)>0)
//			outView = abs(outView)<1e-12 ? 0 : outView
//		endif
//		DisplayTableOfWave(outView)
//	endif
//End


Static Function/WAVE RotsBetweenZones(pQsIn,mQsIN, matchMin, tol)
	// returns ALL the rotations that take at least matchMin pQs --> mQs
	Wave pQsIn				// Q^s on a possible zone
	Wave mQsIN				// Q^s on a measured zone
	Variable matchMin		// number of required Q^s to match
	Variable tol			// anglular tolerance for Q^ matching (radian)

	Variable nP=DimSize(pQsIn,0), nM=DimSize(mQsIN,0)
	if (numtype(nP+nM+matchMin+tol) || nP<matchMin || nM<matchMin || matchMin<=0 || tol<=0)
		return $""
	endif

	// find the min acceptable separation between two q's to calculate a rotation
	Make/N=(nM,3)/D/FREE mQs=mQsIN[p][q]
	Make/N=3/D/FREE qvec
	Variable i, minDot=Inf, tooCloseDot=(1-0.01*PI/180)
	for (i=0;i<nM;i+=1)
		qvec = mQs[0][p]
		MatrixOP/FREE dots = mQs x qvec
		dots = abs(dots) > tooCloseDot ? NaN : dots	// 0.999825 = 1-0.01*PI/180
		minDot = min(minDot,WaveMin(dots))
	endfor
	minDot = cos(acos(limit(minDot,-1,1))/10)	// two q's must be searated by at least minDot to calc rot

	Variable Nalloc = 100
	Make/N=(Nalloc,3)/D/FREE AllAxes=NaN			// {rx,ry,rz},   rotation vectors
	Make/N=3/D/FREE qM0, qM1, qP0, qP1, rotAxis
	Variable dotM, angleM, angleP
	Variable m0, m1, p0,p1, Nall
	for (m0=0,Nall=0; m0<(nM-1); m0+=1)			// m0 is first measured direction
		qM0 = mQsIN[m0][p]
		for (m1=1; m1<nM; m1+=1)					// m1 is second measured direction
			qM1 = mQsIN[m1][p]
			dotM = MatrixDot(qM0,qM1)
			if (dotM>minDot)							// qM0 ^ qM1 is too small
				continue
			endif
			angleM = acos(limit(dotM,-1,1))

			for (p0=0; p0<(nP-1); p0+=1)			// p0 is first possible direction
				qP0 = pQsIN[p0][p]
				for (p1=1;p1<nP;p1+=1)				// m1 is second possible direction
					qP1 = pQsIN[p1][p]
					angleP = acos(limit(MatrixDot(qP0,qP1),-1,1))
					if ( abs(angleP-angleM) > tol)
						continue
					endif

					// calc the rotation and accumulate the rotAxis
					Wave rotAxis = rotFromPairs(qM0,qM1,qP0,qP1)
					if (Nall>=Nalloc)
						Nalloc += 300
						Redimension/N=(Nalloc,-1) AllAxes
					endif
					allAxes[Nall][] = rotAxis[q]
					Nall += 1
				endfor
			endfor
		endfor
	endfor
	Redimension/N=(Nall,-1) AllAxes
	if (Nall<=1)
		return $""
	endif
	Make/N=(Nall,4)/D/FREE axes=NaN				// {rx,ry,rz, matches},   rotation vectors & number of matches
	Make/N=3/D/FREE rotSum
	Variable m, match, Naxes, tol2=(tol*tol)
	for (i=0,Naxes=0; i<Nall; i+=1)				// group by matching axis
		rotAxis = AllAxes[i][p]
		AllAxes[i][] = NaN
		MatrixOP/FREE dist2 = sumRows( magSqr(AllAxes - RowRepeat(rotAxis,Nall)) )
		rotSum = rotAxis
		for (m=i+1,match=1; m<Nall; m+=1)
			if (dist2[m]<tol2)
				match += 1
				rotSum += AllAxes[m][p]
				AllAxes[m][] = NaN
			endif
		endfor
		if (match >= matchMin)						// this is a good one, save it
			axes[Naxes][0,2] = rotSum[q]/match		// save the average axis
			axes[Naxes][3] = match						//    and number of matches
			Naxes += 1
		endif
	endfor
	if (Naxes<1)
		return $""
	endif
	Redimension/N=(Naxes,-1) axes
	return axes
End



//       ==========================================================================       //
//       ========================= Start of Zone Structure ========================       //

Static Constant STRUCT_ZONESPOTS_MAX=50	// max number of Ghats in a single zone
//Static Constant STRUCT_ZONES_MAX=100		// max number of zones
Static Constant STRUCT_ZONES_MAX=50		// max number of zones
//
Static Structure allZonesStruct					// structure definition of a measured or calculated zone
	int16 measured								// 1=measured, 0=calculated
	int16 Nz										// actual number of zones defined here
	int16 NsMin									// smallest number of spots in any one zone
	Struct oneZoneStruct zn[STRUCT_ZONES_MAX]	// q vectors in this zone
	double precision							// precision used when filling the structure
	char GhatSource[200]					// name of wave used to create the zone info
	char desc[200]								// an optional description
EndStructure
//
Static Structure oneZoneStruct			// defines ONE zone
	int32 Ns										// number of spots in this zone
	double axisHat[3]							// unit vector pointing in zone axis direction
	int32 axisijk[3]							// ijk for the axis (if it is known)
	Struct double3vecStruct g[STRUCT_ZONESPOTS_MAX]	// g^s for this zone point
	Struct int3vecStruct h[STRUCT_ZONESPOTS_MAX]		// hkls for this zone point
	// the {hkl} are only filled when measured==0, or for test data, for measured spots, hkl is not known.
EndStructure
//
Static Structure double3vecStruct		// structure definition of a single 3-vector (double precision float)
	double xyz[3]
EndStructure
//
Static Structure int3vecStruct			// structure definition of a single 3-vector (4byte int)
	int32 hkl[3]
EndStructure
//
Static Function ZonesStructClear(zones)
	STRUCT allZonesStruct &zones
	zones.measured	 = 0
	zones.Nz	 = 0
	zones.NsMin = 0
	Variable m, i
	for (m=0;m<STRUCT_ZONES_MAX;m+=1)
		zones.zn[m].axisHat[0] = NaN
		zones.zn[m].axisHat[1] = NaN
		zones.zn[m].axisHat[2] = NaN
		zones.zn[m].axisijk[0] = 0
		zones.zn[m].axisijk[1] = 0
		zones.zn[m].axisijk[2] = 0
		zones.zn[m].Ns = 0
		for (i=0;i<STRUCT_ZONESPOTS_MAX;i+=1)
			zones.zn[m].g[i].xyz[0] = NaN
			zones.zn[m].g[i].xyz[1] = NaN
			zones.zn[m].g[i].xyz[2] = NaN
			zones.zn[m].h[i].hkl[0] = 0
			zones.zn[m].h[i].hkl[1] = 0
			zones.zn[m].h[i].hkl[2] = 0
		endfor
	endfor
	zones.precision = NaN
	zones.GhatSource = ""
	zones.desc = ""
	return 0
End
//
Static Function ZonesStructCopy(dest,source)
	STRUCT allZonesStruct &dest
	STRUCT allZonesStruct &source
	ZonesStructClear(dest)			// start fresh
	dest.measured = source.measured
	dest.Nz = source.Nz
	dest.NsMin = source.NsMin
	Variable m, i
	for (m=0;m<STRUCT_ZONES_MAX;m+=1)
		dest.zn[m].axisHat[0] = source.zn[m].axisHat[0]
		dest.zn[m].axisHat[1] = source.zn[m].axisHat[1]
		dest.zn[m].axisHat[2] = source.zn[m].axisHat[2]
		dest.zn[m].axisijk[0] = source.zn[m].axisijk[0]
		dest.zn[m].axisijk[1] = source.zn[m].axisijk[1]
		dest.zn[m].axisijk[2] = source.zn[m].axisijk[2]
		dest.zn[m].Ns = source.zn[m].Ns
		for (i=0;i<STRUCT_ZONESPOTS_MAX;i+=1)
			dest.zn[m].g[i].xyz[0] = source.zn[m].g[i].xyz[0]
			dest.zn[m].g[i].xyz[1] = source.zn[m].g[i].xyz[1]
			dest.zn[m].g[i].xyz[2] = source.zn[m].g[i].xyz[2]
			dest.zn[m].h[i].hkl[0] = source.zn[m].h[i].hkl[0]
			dest.zn[m].h[i].hkl[1] = source.zn[m].h[i].hkl[1]
			dest.zn[m].h[i].hkl[2] = source.zn[m].h[i].hkl[2]
		endfor
	endfor
	dest.precision = source.precision
	dest.GhatSource = source.GhatSource
	dest.desc = source.desc
	return 0
End
//
Static Function/T ZonesStruct2Str(zi,[maxPrint])
	STRUCT allZonesStruct &zi
	Variable maxPrint
	maxPrint = ParamIsDefault(maxPrint) || numtype(maxPrint) || maxPrint<=1 ? 10 : maxPrint

	Variable addTerm=0
	String str, out=""
	if (strlen(zi.desc))
		out += zi.desc+"\r"
	endif
	sprintf str, "%d %s Zones, (smallest zone has %d spots)\r", zi.Nz, SelectString(zi.measured,"Calculated", "Measured"), zi.NsMin
	out += str
	addTerm = 0
	if (numtype(zi.precision)==0 && zi.precision>0)
		sprintf str, "  computed with precision = %g¡   ",zi.precision
		out += str
		addTerm = 1
	endif
	if (strlen(zi.GhatSource))
		out += "  Ghat source = \""+zi.GhatSource+"\""
		addTerm = 1
	endif
	out += SelectString(addTerm,"","\r")
	out += "\r"
	Make/N=3/D/FREE vec
	Variable m, i
	for (m=0;m<zi.Nz;m+=1)
		out += "zone "+num2istr(m)+" "+OneZoneStruct2Str(zi.zn[m])
	endfor
	return out
End
//
Static Function/T OneZoneStruct2Str(zn,[maxPrint])
	STRUCT oneZoneStruct &zn
	Variable maxPrint
	maxPrint = ParamIsDefault(maxPrint) || numtype(maxPrint) || maxPrint<=1 ? 10 : maxPrint

	String str, out=""
	if (zn.Ns < 1)
		sprintf str, "contains 0 spots\r"
		out += str
		return out
	endif

	Make/N=3/D/FREE vec
	vec = zn.axisHat[p]
	sprintf str, "contains %d spots, with an axis = %s",zn.Ns, vec2str(vec,zeroThresh=1e-5,sep=", ")
	out += str
	vec = zn.axisijk[p]
	if (norm(vec)>1e-5)
		sprintf str, " = direct x [%s]",vec2str(vec,bare=1,sep=",")
		out += str
	endif
	out += "\r"

	str = "  G^s = "
	Variable i
	for (i=0;i<min(maxPrint, zn.Ns);i+=1)
		vec = zn.g[i].xyz[p]
		str += SelectString(i,"",", ")
		str += vec2str(vec,zeroThresh=1e-5,maxPrint=maxPrint,sep=",")
	endfor
	out += str + SelectString(zn.Ns > maxPrint, "", "...") + "\r"

	vec = zn.h[0].hkl[p]
	if (norm(vec))
		str = "  hkls = "
		for (i=0;i<min(maxPrint, zn.Ns);i+=1)
			vec = zn.h[i].hkl[p]
			str += SelectString(i,"",", ")
			str += vec2str(vec,zeroThresh=1e-5,maxPrint=maxPrint,sep=",")
		endfor
		out += str + SelectString(zn.Ns > maxPrint, "", "...") + "\r"
	endif
	return out
End
//
//	Function testStruct()
//		STRUCT allZonesStruct zi
//		ZonesStructClear(zi)
//		zi.measured = 1
//		zi.Nz = 2
//		zi.NsMin = 2
//		zi.precision =  0.1
//		zi.GhatSource = "fullWavePath:name"
//		zi.desc = "test print"
//
//		zi.zn[0].Ns	= 3
//		zi.zn[0].axisHat[0] = 0.1	;	zi.zn[0].axisHat[1] = 0.2 ; zi.zn[0].axisHat[2] = 0.3
//		zi.zn[0].g[0].xyz[0] = 0.11	;	zi.zn[0].g[0].xyz[1] = 0.12	;	zi.zn[0].g[0].xyz[2] = 0.13
//		zi.zn[0].g[1].xyz[0] = 0.21	;	zi.zn[0].g[1].xyz[1] = 0.22	;	zi.zn[0].g[1].xyz[2] = 0.23
//		zi.zn[0].g[2].xyz[0] = 0.31	;	zi.zn[0].g[2].xyz[1] = 0.32	;	zi.zn[0].g[2].xyz[2] = 0.33
//
//		zi.zn[1].Ns	= 2
//		zi.zn[1].axisHat[0] = -0.1	;	zi.zn[1].axisHat[1] = -0.2	;	zi.zn[1].axisHat[2] = -0.3
//		zi.zn[1].g[0].xyz[0] = 0.51	;	zi.zn[1].g[0].xyz[1] = 0.52	;	zi.zn[1].g[0].xyz[2] = 0.53
//		zi.zn[1].h[0].hkl[0] = 1		;	zi.zn[1].h[0].hkl[1] = 2		;	zi.zn[1].h[0].hkl[2] = 3
//		zi.zn[1].g[1].xyz[0] = 0.61	;	zi.zn[1].g[1].xyz[1] = 0.62	;	zi.zn[1].g[1].xyz[2] = 0.63
//		zi.zn[1].h[1].hkl[0] = -1		;	zi.zn[1].h[1].hkl[1] = -2		;	zi.zn[1].h[1].hkl[2] = -3
//		print ZonesStruct2Str(zi)
//	End

//       ========================== End of Zone Structure =========================       //
//       ==========================================================================       //

//  ============================= End of Zones+Qs Indexing =============================  //
//  ====================================================================================  //




//  ====================================================================================  //
//  =============================== Start of Qs Indexing ===============================  //

Function/WAVE runIndexingQs(args)
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
	elseif (angTol<=0 || numtype(angTol) || angTol>10)		// angTol must be in (0,10), and 10 is really BIG
		return $""
	elseif (numtype(hp+kp+lp))
		return $""
	elseif (keVmaxCalc<=0 || numtype(keVmaxCalc))
		return $""
	elseif (keVmaxTest<=0 || numtype(keVmaxTest))
		return $""
	endif
	timeIncrement(fresh=1)

	Make/N=3/D/FREE hklPrefer={hp,kp,lp}							// hkl near center of data
	Variable NmaxHKL = 500												// max number of Possible hkl to create
	Make/N=3/D/FREE ki={0,0,1}										// use this for the incident beam direction
	Wave kfsMeasured = FullPeakList2kfHats(maxSpots,FullPeakList,FullPeakList1=FullPeakList1,FullPeakList2=FullPeakList2)	// convert peaks to kf^s
	Wave GhatsMeasured=FullPeakList2Ghats(maxSpots,FullPeakList,FullPeakList1=FullPeakList1,FullPeakList2=FullPeakList2)	// convert peaks to a Qlist+intens+dNum
	String noteMeasured=note(GhatsMeasured)
	Variable NG0=DimSize(GhatsMeasured,0)
	Make/N=(NG0,3)/FREE GhatsOnly=GhatsMeasured[p][q]
	MatrixOP/FREE kf0 = Normalize(sumCols(kfsMeasured))		// assume that average of kfsMeasured is kf-vector to center of detector
	Redimension/N=3 kf0
	Make/N=3/D/FREE qvec0 = (kf0 - ki)								// try the usual definition of q
	if (norm(qvec0)<0.001)												// kf0 & ki are almost parallel
		Wave qvec0 = perp2Vector(ki)									// qvec0 is just perpendicular to ki (in some unspecified direction)
	endif
	Variable measuredSpan = AngularSpanOfVectors(GhatsOnly)

#ifdef QS_TESTING
	if (printIt)
		printf "*made GhatsMeasured[%d][%d] with a angular span of %.2f¡, %s\r",NG0,DimSize(GhatsMeasured,1),measuredSpan,timeIncrement()
		Duplicate/O GhatsOnly, GhatsOnlyView
	endif
#endif

	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))
		DoAlert 0, "no lattice structure found, did you forget to set it?"
		return $""
	endif
	Wave direct=directFrom_xtal(xtal)
	Wave recipBase=recipFrom_xtal(xtal)
	Variable atDetector = hkPreferPointsAtDetector(xtal,kfsMeasured)	// True if hklPrefer points at detector (forward scattering), False for top detector
#ifdef QS_TESTING
	Wave PossibleQhats = MakePossibleQhats(atDetector,recipBase,qvec0,cone*PI/180,hklPrefer,keVmaxCalc,Nmax=NmaxHKL,kin=ki,printIt=printIt)	// q^ available for checking
	Duplicate/O PossibleQhats, PossibleQhatsView
	if (printIt)
		printf "*made PossibleQhats[%d][%d], %s\r",DimSize(PossibleQhats,0),DimSize(PossibleQhats,1),timeIncrement()
	endif
#else
	Wave PossibleQhats = MakePossibleQhats(atDetector,recipBase,qvec0,cone*PI/180,hklPrefer,keVmaxCalc,Nmax=NmaxHKL,kin=ki)	// q^ available for checking
#endif

	// Possible HKL Pairs & Measured Ghat Pairs are lists of {dot,i,j}
	Variable minDot = cos(measuredSpan*PI/180)										// smallest dot products between all of the GhatsMeasured
	Wave MeasuredGhatPairs = MakePairsList(GhatsMeasured)						//   only uses first 3 columns of argument
	Wave PossibleQhatPairs = MakePairsList(PossibleQhats,mindot=minDot)	// all pairs of hkl, {dot, i, j}
	Variable iMeasured, nMeasuredPairs=DimSize(MeasuredGhatPairs,0)
	Variable iPossible, nPossiblePairs=DimSize(PossibleQhatPairs,0)
#ifdef QS_TESTING
	if (printIt)
		printf "*made both pairs waves: PossibleQhatPairs[%d][3], MeasuredGhatPairs[%d][3], %s\r",nPossiblePairs,nMeasuredPairs,timeIncrement()
	endif
#endif

	// get all of the rotations between pairs, make a list of all the rotation vectors (radian)
	Variable angTolRad = angTol*PI/180								// used to compare two directions for equivalence
	Variable Nalloc = 1000
	Make/N=(Nalloc,4)/U/I/FREE indexPairs=0						// holds (i,j of measured pair) &  (i,j of calculated pair)
	Variable angleMeasured, nPairs
	for (iMeasured=0,nPairs=0; iMeasured<nMeasuredPairs; iMeasured+=1)
		angleMeasured = acos(MeasuredGhatPairs[iMeasured][0])
		for (iPossible=0; iPossible<nPossiblePairs; iPossible+=1)
			if (abs(angleMeasured-acos(PossibleQhatPairs[iPossible][0]))<angTolRad)// do the pairs have the same angular separation?
				if (nPairs>=Nalloc)
					Nalloc += 1000
					Redimension/N=(Nalloc,-1) indexPairs
				endif
				indexPairs[nPairs][0] = MeasuredGhatPairs[iMeasured][1]	// save this set of pairs
				indexPairs[nPairs][1] = MeasuredGhatPairs[iMeasured][2]
				indexPairs[nPairs][2] = PossibleQhatPairs[iPossible][1]
				indexPairs[nPairs][3] = PossibleQhatPairs[iPossible][2]
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
		printf "*made pairs of pairs: indexPairs[%d][4], elapsed time = %s\r",nPairs, timeIncrement()
	endif
	//	Duplicate/O indexPairs, indexPairsView
#endif
	if (nPairs<1)
		printf "ERROR -- nPairs = %g,  must be ³ 2\r",nPairs
		return $""
	endif
	// each row of indexPairs generates two rotations (+ and -), so there are 2*nPairs of rotations to check
	//   generate the rotation vectors, store them in PairRotations[2*nPairs][3]
	Wave PairRotations = ConstructAllPairRotations(indexPairs,GhatsMeasured,PossibleQhats)
#ifdef QS_TESTING
	if (printIt)
		printf "*made PairRotations[%d][3] the rotation vectors, %s\r",DimSize(PairRotations,0), timeIncrement()
	endif
	if (DataFolderExists(":GizmoWaves"))
		Duplicate/O PairRotations, :GizmoWaves:PairRotationsView
	endif
#endif

	// check for lumps or groupings in the set of rotation vectors PairRotations[][3]
	Variable threshStart=angTolRad			// largest acceptable error between measured and calculated angles (radian)
	Variable threshStop=0.001*PI/180		// stop searching when thresh reaches 0.001¡
#ifdef QS_TESTING
	if (printIt)
		printf "running with a starting threshold of %g¡, and a stopping threshold of %g¡\r",threshStart*180/PI, threshStop*180/PI
	endif
#endif
	Make/N=3/D/FREE rotAxis, center=0		// start center at {0,0,0}
	Variable hits, hitsLast=Inf, thresh, radius=Inf	// hits is largest number of PairRotations that are the same rotation
	for (thresh=threshStart; thresh>threshStop; thresh /= threshDivide)
		hits = FindPeakInRots(PairRotations,rotAxis,thresh,center,radius)
#ifdef QS_TESTING
		if (printIt)
			Variable indexOfMax = NumberByKey("indexOfMax",note(rotAxis),"=")
			printf "  hits = %g,   hitsLast = %g,  thresh=%.3g¡,  radius=%.3g¡,   axis=%s  (%d), %s\r",hits,hitsLast,thresh*180/PI,radius*180/PI,vec2str(rotAxis,zeroThresh=1e-9),indexOfMax,timeIncrement()
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
		printf "*found <rot> = %s  |<rot>|=%.3f¡,  with %d hits, %s\r",vec2str(rotAxis,zeroThresh=1e-9),angle,hits,timeIncrement()
	endif
#endif

	// find rms error between Ghats and the predicted hkl directions
	STRUCT microGeometry g
	if (FillGeometryStructDefault(g, alert=1))	//fill the geometry structure with current values
		return $""
	endif

	Make/N=3/D/FREE vec, hkl
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
		printf "*after non-linear least-squares optimization, <rot> = %s  |<rot>|=%.3f¡,  %s\r",vec2str(rotAxis,zeroThresh=1e-9),angle,timeIncrement()
	endif
#endif

	// Save all indexed values using the optimized rotation
	Variable Npatterns=1, ipat=0, NG
	String peakListName = ParseFilePath(3,StringByKey("peakListWave",noteMeasured,"="),":",0,0)
	String IndexedName = CleanupName("FullPeakIndexed"+ReplaceString("FullPeakList",peakListName,""),0)
	Make/N=(NG0,12,Npatterns)/O $IndexedName/WAVE=IndexedWave = NaN
	Make/N=3/D/O/FREE hkl
	Variable dNum, px,py, keV, intensity, intensityMax=-Inf
	Variable rmsGhat, minError=Inf, maxError=-Inf
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
			intensityMax = numtype(intensity) ? intensityMax : max(intensityMax,intensity)
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
	Variable goodness0 = goodNessOfPattern(IndexedWave,ipat)
	rmsGhat = sqrt(rmsGhat/NG)
	IndexedWave[][8][ipat] = abs(IndexedWave[p][8][ipat])<5e-6 ? 0 : IndexedWave[p][8][ipat]	// precision limit
#ifdef QS_TESTING
	printf "Ghat rms error (from %d of %d) Ghats = %.3g¡,  range=[%.4f, %.4f¡],  goodness0=%g\r",NG,NG0,rmsGhat,minError,maxError,goodness0
	printf "Total time = %s,   largest Ætime = %.3f sec\r",timeIncrement(),NumVarOrDefault("root:Packages:timeIncrement:maxDelta",NaN)
#endif
	String executionTimeStr = timeIncrement(seconds=1)

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
	if (strlen(xtal.SpaceGroupID))
		wnote = ReplaceStringByKey("SpaceGroupID",wnote,xtal.SpaceGroupID,"=")
	endif
	if (numtype(xtal.SpaceGroupIDnum))
		wnote = ReplaceNumberByKey("SpaceGroupIDnum",wnote,xtal.SpaceGroupIDnum,"=")
	endif

	wnote = ReplaceNumberByKey("keVmaxCalc",wnote,keVmaxCalc,"=")
	wnote = ReplaceNumberByKey("keVmaxTest",wnote,keVmaxTest,"=")
	wnote = ReplaceStringByKey("hklPrefer",wnote,vec2str(hklPrefer),"=")
	wnote = ReplaceNumberByKey("cone",wnote,cone,"=")
	wnote = ReplaceNumberByKey("angleTolerance",wnote,angTol,"=")
	wnote = ReplaceNumberByKey("Nindexed",wnote,NG,"=")
	wnote = ReplaceNumberByKey("NiData",wnote,NG0,"=")
	wnote = ReplaceStringByKey("executionTime",wnote,executionTimeStr,"=")
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
// qvec0 comes from the average of all the measured qvectors, for forward scattering (ki¥qvec0)==0
//
Static Function/WAVE MakePossibleQhats(atDetector,recip0,qvec0,cone,hkl0,keVmax[,keVmin,Nmax,kin,printIt])
	// Make a list of q^[][3] that are available for checking (in UN-Rotated frame, the frame of recip0)
	Variable atDetector			// True if hklPrefer points at detector (forward scattering), false means hklPrefer is qCenter
	Wave recip0						// UN-rotated reciprocal lattice, will be diagonal, for cubic, tetragonal, & orthorhombic
	Wave qvec0						// this q-vector diffracts to center of the detector (only use its direction)
	Variable cone					// only consider q-vectors with kf within cone of qvec0+ki (radian)
	Wave hkl0						// preferred hkl, the hkl specidifed by the user, this defines orientation
	Variable keVmax				// maximum allowed energy (keV)
	Variable keVmin				// minimum    "       "   (keV)
	Variable Nmax					// max number of test q-vectors, defaults to 250
	Wave kin							// optional incident wave direction (defaults to 001)
	Variable printIt
	Nmax = ParamIsDefault(Nmax) || numtype(Nmax) ? 250 : round(Nmax)
	keVmin = ParamIsDefault(keVmin) || keVmin<=0 || numtype(keVmin) ? 0 : keVmin
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : !(!printIt)

	Variable tol = 0.1*PI/180						// angle tolerance (radian)
	cone = min(cone,PI-tol)						// limit cone to a bit less than 180¡
	if (numtype(cone+keVmax+keVmin+Nmax))
		return $""
	elseif (cone<tol || cone>(PI-tol))			// cone in range [tol, 180¡-tol]
		return $""
	elseif (keVmax<INDEXING_KEV_MIN || keVmin>=keVmax || keVmax>INDEXING_KEV_MAX)
		return $""
	elseif (Nmax<3 || Nmax>500)					// max number of hkls is in range [3,500]
		return $""
	endif
	keVmin = keVmin==0 ? keVmin : keVmax/20

	Make/N=3/D/FREE ki={0,0,1}					// use this for the incident beam direction
	if (!ParamIsDefault(kin))						// one was passed, use it
		ki = kin[p]
		normalize(ki)
	endif

	MatrixOP/FREE qhat0 = Normalize(qvec0)	// direction of q-vector that will scatter toward detector center
	Redimension/N=3 qhat0							// changes [1,3] or [3,1] into just [3]
	Variable ki_dot_q = MatrixDot(ki,qhat0)
	Make/N=3/D/FREE kf0 = (ki - 2*ki_dot_q*qhat0)//		kf^ = ki^ - 2*(ki^ . q^)*q^

	// find rotation that takes (recip0 x hkl0) --> qvec0
	//		rot0 x (recip0 x hkl0) == qhat0
	MatrixOP/FREE q0 = Normalize(recip0 x hkl0)
	Variable sine, cosine, angle
	if (atDetector)									// forward scattering
		// compute rotation rot0 so (rot0 x recip0) = recip,  (recip x hklPrefer) points at kf0 (detector center)
		Cross q0, kf0									// find a rotation axis that puts q0=Normalize(recip0 x hkl0) parallel to kf0
		cosine = MatrixDot(q0,kf0)				// cosine of angular distance to rotate q0
	else													// top detector
		// compute rotation rot0 so (rot0 x recip0) = recip,  (recip x hklPrefer) points at q0 (q0 = kf0-ki), i.e. q0 diffracts to detector center
		Cross q0, qhat0								// find a rotation axis that puts q0=Normalize(recip0 x hkl0) parallel to qhat0
		cosine = MatrixDot(q0,qhat0)				// cosine of angular distance to rotate q0
	endif
	Wave W_Cross=W_Cross
	sine = normalize(W_Cross)
	angle = atan2(sine,cosine)*180/PI
	Make/N=(3,3)/D/FREE rot0
	rotationMatAboutAxis(W_Cross,angle,rot0)
	KillWaves/Z W_Cross
	MatrixOP/FREE recip = rot0 x recip0		// recip is now a R.L. that has (recip x hkl0) pointing in the desired direction
	//		so now:  qhat0 || (rot0 x recip0 x hkl0)

	// I only want kf's that are within cone of kf0
	Variable maxTheta = limit((acos(MatrixDot(kf0,ki))+cone)/2,tol,PI/2-tol)
	Variable maxQLen = 4*PI*sin(maxTheta) * keVmax/hc	// max len of a q-vector = (recip0 x hkl)
	Variable minTheta = limit((acos(MatrixDot(kf0,ki))-cone)/2, 0, maxTheta)
	Variable minQLen = 4*PI*sin(minTheta) * keVmin/hc	// min len of a q-vector = (recip0 x hkl)
	MatrixOP/FREE hklMax = floor( rowRepeat(maxQLen,3)/sqrt(sumCols(magSqr(recip)))^t )
	Variable hmax=hklMax[0], kmax=hklMax[1], lmax=hklMax[2]
	if (printIt)
		printf "  find Possible Q-directions using index range: h=[%+d, %+d], k=[%+d, %+d], l=[%+d, %+d]\r",-hmax,hmax, -kmax,kmax, -lmax,lmax
		printf "  hklPrefer %s detector center,  theta: [%.2f¡, %.f¡],  Q: [%.1f, %.1f](1/nm)\r", SelectString(atDetector, "difracts to", "points at"), minTheta*180/PI, maxTheta*180/PI,minQLen,maxQLen
		printf "  kf0 = %s,  hkl0 = (%s)\r",vec2str(kf0),vec2str(hkl0,sep=" ",bare=1)
	endif

	Variable minDot=cos(cone)						// minimum allowed value of dot(kf0,kf)
	Variable cosTol=cos(tol), coneMax=-Inf
	Make/N=(Nmax,3)/D/FREE TestQs=NaN			// contains q^x, q^y, q^z, other has h,k,l
	Make/N=3/D/FREE hkl, kf


//		elseif (cone<PI/3 && keV>keVmax)// for cone angles < 60¡, also check for valid energy

	
	Make/N=3/D/FREE hklLo={-hmax,-kmax,-lmax}, hklHi={hmax,kmax,lmax}	// range of h,k,l
	Wave hkls = hklOrderedList(NaN,hklLo=hklLo,hklHi=hklHi)
	Variable Nhkls=DimSize(hkls,0)
	Variable N, Qmag, m, kf_dot_kf0				// = MatrixDot(kf,kf0)

#if defined(ZONE_TESTING) || defined(QS_TESTING) || defined(ZONE_QS_TESTING)
	if (DEBUG_LEVEL>=3 && printIt)
		printf "\r\t  hkl\t\t  acos(ki¥q^)\tacos(kf¥kfo^)\t\t\t\t\tkf^\r"
	endif
#endif
	for (m=0,N=0; m<Nhkls && N<Nmax; m+=1)
		hkl = hkls[m][p]
		MatrixOP/FREE/O qhat = recip x hkl
		Qmag = normalize(qhat)
		ki_dot_q = MatrixDot(ki,qhat)
		kf = (ki - 2*ki_dot_q*qhat)				// kf^ = ki^ - 2*(ki^ . q^)*q^
		kf_dot_kf0 = MatrixDot(kf,kf0)
		if (Qmag>maxQLen || Qmag<minQLen)
			continue										// length is too long, skip
		elseif (ki_dot_q > 0)						// (-ki .dot. q^) < 0 is invalid, e.g. NEVER need (100) & (-100)
			continue										// Bragg angle must be < 90¡
		elseif (kf_dot_kf0<minDot)
			continue										// kf is too far from kf0
		elseif (!isNewDirection(TestQs,qhat,cosTol))
			continue										// skip existing directions too
		endif
#if defined(ZONE_TESTING) || defined(QS_TESTING) || defined(ZONE_QS_TESTING)
		if (DEBUG_LEVEL>=3 && printIt)
			printf "(%s)   \t%.3f¡\t\t\t%.3f¡\t\t%s\r",vec2str(hkl,bare=1,sep=" "),acos(ki_dot_q)*180/PI, acos(limit(kf_dot_kf0,-1,1))*180/PI,vec2str(kf,fmt="% .5f")
		endif
#endif
		TestQs[N][] = qhat[q]						// a new valid & unique q-direction
		N += 1
		coneMax = max(coneMax,acos(limit(kf_dot_kf0,-1,1)))
	endfor
	Redimension/N=(N,-1) TestQs
	MatrixOP/FREE/O TestQs = ( Inv(rot0) x TestQs^t )^t	// rotate TestQs back to UN-rotated frame
	if (printIt)
		printf "  made %d test Q's with a maximum cone of %.2f¡ and an angular span of %.2f¡\r",N,coneMax*180/PI,AngularSpanOfVectors(TestQs)
	endif

	SetDimLabel 1,0,Qx,TestQs						// first 3 columns are the normalized qx,qy,qz
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
//
//	Function Test_MakePossibleQhats(cone,keVmax)		// try: Test_MakePossibleQhats(40, 12)
//		Variable cone			//	=40
//		Variable keVmax		//	=12
//	
//		STRUCT crystalStructure xtal
//		if (FillCrystalStructDefault(xtal))
//			DoAlert 0, "no lattice structure found, did you forget to set it?"
//			return 1
//		endif
//		Wave recip=recipFrom_xtal(xtal)
//	
//		Make/D/FREE qvec0={0,1,-1}, hkl0={0,0,2}
//	
//		Wave FullPeakList=FullPeakListOrange
//		if (!WaveExists(FullPeakList))
//			DoAlert 0, "FullPeakListOrange could not be found."
//			return 1
//		endif
//		Wave GhatsMeasured=IndexingInternal#FullPeakList2Ghats(250,FullPeakList)	// convert peaks to a Qlist+intens+dNum
//		Variable NG0=DimSize(GhatsMeasured,0)
//		Make/N=(NG0,3)/FREE GhatsOnly=GhatsMeasured[p][q]
//		MatrixOP/FREE qvec0 = sumCols(GhatsOnly)					// assume that average Ghat is Q-vector to center of detector
//		Redimension/N=3 qvec0
//		Variable atDetector = IndexingInternal#hkPreferPointsAtDetector(xtal,GhatsOnly)	// True if hklPrefer points at detector (forward scattering), False for top detector
//	
//		timeIncrement(fresh=1)
//		Wave qhats = IndexingInternal#MakePossibleQhats(atDetector,recip,qvec0,cone*PI/180,hkl0,keVmax,printIt=1)
//		print "Total",timeIncrement()
//		DuplicaTe/O qhats, qhatsView
//	End

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

Function/WAVE runIndexingOnlyZones(args)
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
	Variable threshStop=0.001*PI/180		// at 0.001¡ stop searching for a better match
#ifdef ZONE_TESTING
	if (printIt)
		printf "running with a starting threshold of %g¡, and a stopping threshold of %g¡\r",threshStart*180/PI, threshStop*180/PI
	endif
#endif
	Make/N=3/D/FREE rotAxis, center=0		// start center at {0,0,0}
	Variable hits, hitsLast=Inf, thresh, radius=Inf	// hits is largest number of PairRotations that are the same rotation
	for (thresh=threshStart; thresh>threshStop; thresh /= threshDivide)
		hits = FindPeakInRots(PairRotations,rotAxis,thresh,center,radius)
#ifdef ZONE_TESTING
		if (printIt)
			Variable indexOfMax = NumberByKey("indexOfMax",note(rotAxis),"=")
			printf "  hits = %g,   hitsLast = %g,  thresh=%.3g¡,  radius=%.3g¡,   axis=%s  (%d)\r",hits,hitsLast,thresh*180/PI,radius*180/PI,vec2str(rotAxis,zeroThresh=1e-9),indexOfMax
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
	if (FillGeometryStructDefault(g, alert=1))	//fill the geometry structure with current values
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
	Variable NG, rmsGhat, minError=Inf, maxError=-Inf
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
			intensityMax = numtype(intensity) ? intensityMax : max(intensityMax,intensity)
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
	Variable goodness0 = goodNessOfPattern(IndexedWave,ipat)
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
	if (strlen(xtal.SpaceGroupID))
		wnote = ReplaceStringByKey("SpaceGroupID",wnote,xtal.SpaceGroupID,"=")
	endif
	if (numtype(xtal.SpaceGroupIDnum))
		wnote = ReplaceNumberByKey("SpaceGroupIDnum",wnote,xtal.SpaceGroupIDnum,"=")
	endif

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
		if (printIt)
			printf "for Nmin = %d,  No zones found\r",Nmin
		endif
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


Function/WAVE PutZoneLinesOnGraph(FullPeakList,[Nmin,tolAngle,printIt])
	// Find the zones from fitted peaks, and put lines on plot showing the zones.
	Wave FullPeakList
	Variable Nmin			// minimum number of spot required to indentify a zone
	Variable tolAngle		// angular tolerance used in accepting Q-vec as part of a zone
	Variable printIt
	Nmin = ParamIsDefault(Nmin) || numtype(Nmin) || Nmin<3 ? 7 : Nmin
	tolAngle = ParamIsDefault(tolAngle) || numtype(tolAngle) || tolAngle<=0 ? 1 : tolAngle
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : !(!printIt)

	STRUCT microGeometry g
	if (FillGeometryStructDefault(g, alert=1))	//fill the geometry structure with current values
		return $""
	endif
	String FullPeakName=""
	if (!WaveExists(FullPeakList))		// need to find the FullPeakList
		Wave FullPeakList = $StringByKey("FullPeakList",GetUserData("","","FitPeaks"),"=")
		Prompt Nmin,"min number of spots on a zone [3,100]"
		Prompt tolAngle, "angle tolerance in defining a zone (degree)"
		if (!WaveExists(FullPeakList))
			String list = WaveListClass("FittedPeakList","*","MINCOLS:11")
			if (ItemsInList(list)<1)
				return $""
			elseif (ItemsInList(list)==1)
				Wave FullPeakList = $StringFromList(0,list)
			else
				Nmin = limit(Nmin,3,100)
				Prompt FullPeakName,"Full Peak List",popup, list
				DoPrompt "Peak List",FullPeakName, Nmin, tolAngle
				if (V_flag)
					return $""
				endif
				Wave FullPeakList = $FullPeakName
			endif
		endif
		if (ParamIsDefault(Nmin) || ParamIsDefault(tolAngle) || numtype(Nmin+tolAngle) || Nmin<3 || tolAngle <=0)
			DoPrompt "Zone Finding Parameters",Nmin, tolAngle
			if (V_flag)
				return $""
			endif
		endif
		printIt = 1
	endif
	FullPeakName = NameOfWave(FullPeakList)
	if (printIt)
		printf "PutZoneLinesOnGraph(%s",FullPeakName
		if (!ParamIsDefault(Nmin) || Nmin!=7)
			printf ", Nmin=%g", Nmin
		endif
		if (!ParamIsDefault(tolAngle) || tolAngle!=1)
			printf ", tolAngle=%g", tolAngle
		endif
		printf ")\r"
	endif
	if (!WaveExists(FullPeakList))
		return $""
	endif

	String ending="Other", outName=""
	if (strsearch(FullPeakName,"FullPeakList",0)==0)
		ending = FullPeakName[12,inf]
	endif
	outName = CleanupName("ZoneLines"+ending,0)

	Variable i, N=DimSize(FullPeakList,0)
	Make/N=(N,3)/FREE Ghats=NaN
	Make/N=3/D/FREE qhat
	for (i=0;i<N;i+=1)
		pixel2q(g.d[0],FullPeakList[i][0],FullPeakList[i][1],qhat)
		Ghats[i][] = qhat[q]
	endfor
	Note/K Ghats, "waveClass=QsMeasured;"
	Wave waxes = MakeZonesWave(Ghats, Nmin, tolAngle=tolAngle, printIt=0)
	if (!WaveExists(waxes))
		return $""
	endif
	Wave zl = MakeZoneLinesForGraph(g.d[0],waxes)
#if defined(ZONE_TESTING) || defined(QS_TESTING) || defined(ZONE_QS_TESTING)
	Duplicate/O waxes, $CleanupName("ZoneAxes"+ending,0)
#endif
	Duplicate/O zl, $outName
	Wave zOut = $outName
	if (printIt)
		printf "  created wave '%s' containing lines for %d zones\r",outName,DimSize(waxes,0)
	endif

	Wave image = $StringByKey("fittedIgorImage",note(FullPeakList),"=")
	String win = StringFromList(0,WindowsWithWave(image,1))
	CheckDisplayed/W=$win zOut
	if (strlen(win) && !V_flag)		// image displayed, but zOut NOT on image
		AppendToGraph/W=$win zOut[*][1] vs zOut[*][0]
		ModifyGraph lSize($outName)=2, lStyle($outName)=7, rgb($outName)=(0,50000,0)
	endif
	return zOut
End
//
Static Function/WAVE MakeZoneLinesForGraph(d,axes,[dAngle,depth])
	STRUCT detectorGeometry &d
	Wave axes					// list of zone axes dim(N,3)
	Variable dAngle			// angular step size along lines (degree)
	Variable depth				// depth of sample (rarely if ever used), probably 0
	dAngle = ParamIsDefault(dAngle) || numtype(dAngle) || dAngle<0 ? NaN : dAngle
	depth = ParamIsDefault(depth) ? NaN : depth
	if (!WaveExists(axes))
		return $""
	endif
	Variable N=DimSize(axes,0)

	if (!(dAngle>0))								// find the angular step size
		Make/N=3/D/FREE q1,q2
		pixel2q(d,d.Nx/2,0,q1)
		pixel2q(d,d.Nx/2,d.Ny-1,q2)
		dAngle = acos(MatrixDot(q1,q2))		// anglular size in y-direction (radian)
		pixel2q(d,0,d.Ny/2,q1)
		pixel2q(d,d.Nx-1,d.Ny/2,q2)
		dAngle = min(acos(MatrixDot(q1,q2)),dAngle)	// anglular size in x-direction (radian)
		dAngle /= 50
//		printf "using an angular step size of %g¡\r",dAngle*180/PI
	endif
	Variable NLineMax = ceil(PI/dAngle)+2
	Make/N=(NLineMax)/FREE pxs,pys,angles
	Make/N=(N*(NLineMax+1),2)/FREE ZoneLines=NaN	// (x,y) parirs of each point on the curve
	SetScale d 0,0,"pixel", ZoneLines

	Make/N=3/D/FREE qDetCenter				// q that diffracts to detector center
	pixel2q(d,d.Nx/2,d.Ny/2,qDetCenter,depth=depth)

	Make/N=3/D/FREE q0, qperp, qvec, axis
	Variable/C pz
	Variable plusEnd, minusEnd, Nline
	Variable iaxis, m, i, sine, cosine, angle
	for (iaxis=0,m=0; iaxis<N; iaxis+=1)	// once for each zone axis
		axis = axes[iaxis][p]
		normalize(axis)
		q0 = qDetCenter - MatrixDot(qDetCenter,axis)*axis
		normalize(q0)								// q on zone that is closest to detector center
//		printf "q0[%d] = %s,   angle=%g¡\r",iaxis,vec2str(q0), acos(MatrixDot(q0,qDetCenter))*180/PI

		Cross axis, q0
		Wave W_Cross=W_Cross
		qperp = W_Cross
		normalize(qperp)							// q on zone that is perpendicular to q0
		plusEnd = 1
		minusEnd = 1
		Nline = 0
		angles = Inf
		for (i=0;i<100 && (plusEnd || minusEnd);i+=1)
			angle = i*dAngle
			cosine = cos(angle)
			sine = sin(angle)
			if (plusEnd)
				qvec = cosine*q0 + sine*qperp
				pz = q2pixel(d,qvec,depth=depth)
				if (pixelOnDetector(d,pz))		// true if more in this direction
					pxs[Nline] = real(pz)
					pys[Nline] = imag(pz)
					angles[Nline] = angle
					Nline += 1
				else
					plusEnd = 0							// done with the plus direction
				endif
			endif
			if (minusEnd && i)
				qvec = cosine*q0 - sine*qperp
				pz = q2pixel(d,qvec,depth=depth)
				if (pixelOnDetector(d,pz))		// true if more in this direction
					pxs[Nline] = real(pz)
					pys[Nline] = imag(pz)
					angles[Nline] = -angle
					Nline += 1
				else
					minusEnd = 0						// done with the minus direction
				endif
			endif
		endfor
		Sort angles, pxs, pys
		ZoneLines[m,m+Nline-1][0] = pxs[p-m]
		ZoneLines[m,m+Nline-1][1] = pys[p-m]
		m += Nline + 1
	endfor
	KillWaves/Z W_Cross

	Redimension/N=(m-1,-1) ZoneLines
	return ZoneLines
End
//
Static Function pixelOnDetector(d,pz)
	STRUCT detectorGeometry &d
	Variable/C pz
	Variable px=real(pz), py=imag(pz)

	if (numtype(px+py))
		return 0
	elseif (px<0 || px>=(d.Nx-1))
		return 0
	elseif (py<0 || py>=(d.Ny-1))
		return 0
	endif
	return 1
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


Static Function/WAVE MakePairsList(vecHats,[minDot])// returns pairDots, whose columns are: {dot, i, j}
	Wave vecHats										// a wave containing unit vectors (only use first 3 columns)
	Variable minDot									// do not accept any pairs with a dot less than minDot (-1 gets everything)
	minDot = ParamIsDefault(minDot) || numtype(minDot) ? -2 : minDot

	Variable notParallel = cos(0.05*PI/180)	// = 0.999999619228249, sufficient for parallel
	Variable Nv=DimSize(vecHats,0)
	Variable Np = Nv*(Nv-1)/2						// number of unique parirs of vecHats, init to max possible value
	Make/N=(Np,3)/D/FREE pairDots=NaN
	Make/N=3/D/FREE vecj, veci
	Variable dot, j,i,m=0
	for (j=0,m=0; j<(Nv-1); j+=1)
		vecj = vecHats[j][p]
		for (i=j+1;i<Nv;i+=1)
			veci = vecHats[i][p]
			dot = MatrixDot(vecj,veci)
			if (abs(dot)<notParallel && dot>minDot)	// vectors are neither parallel nor anti-parallel, and not too far apart
				pairDots[m][0] = dot				// store this pair, save dot and i,j
				pairDots[m][1] = j
				pairDots[m][2] = i
				m += 1			
			endif
		endfor
	endfor
	Np = m
	Redimension/N=(Np,-1) pairDots
	pairDots[][0] = limit(pairDots[p][q],-1,1)	// ensure that dot is in [-1,1]
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





// convert px,py positions on detector into kf vectors, assumes ki={0,0,1}
ThreadSafe Function/WAVE pixel2kfVEC(d,pxpy,[depth,DeltaPixelCorrect])	// returns the normalized kf's in beam line coords
	STRUCT detectorGeometry &d
	Wave pxpy								// list of pixels to use (must be in raw un-binned pixels), perhaps an ROI, first dimension is number of pixels, second is x,y
	Variable depth							// sample depth measured along the beam (optional, 0 is default)
	Wave DeltaPixelCorrect				// contains optional pixel corrections dimensions are [N][2]
	depth = ParamIsDefault(depth) || numtype(depth) ? 0 : depth

	if (!ParamIsDefault(DeltaPixelCorrect) && DimSize(DeltaPixelCorrect,1)==2)
		Wave kf = pixel2XYZVEC(d,pxpy, DeltaPixelCorrect=DeltaPixelCorrect)
	else
		Wave kf = pixel2XYZVEC(d,pxpy)	// kf is in direction of pixel in beam line coords
	endif

	Variable N=DimSize(pxpy,0)
	Make/N=3/D/FREE ki={0,0,1}				//	ki = geo.ki[p],  incident beam direction (normalized)
	if (depth==0)
		MatrixOP/FREE/O kf = NormalizeRows(kf)
	else
		Make/N=3/D/FREE depthVec = depth*ki[p]
		MatrixOP/FREE/O kf = NormalizeRows(kf - rowRepeat(depthVec,N))	// kf(depth) = d*ki + kf(depth=0)
	endif

	WaveClear ki, depthVec
	return kf										// normalized kf-vectors in beam line coords
End





Static Function/WAVE FullPeakList2kfHats(maxSpots,FullPeakList,[FullPeakList1,FullPeakList2,DeltaPixelCorrect])	// convert peaks to kf^s
	Variable maxSpots
	Wave FullPeakList
	Wave FullPeakList1
	Wave FullPeakList2
	Wave DeltaPixelCorrect			// contains optional pixel corrections dimensions are [N][2]

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
	elseif (DimSize(FullPeakList,1)<11)
		DoAlert 0, "the passed full peak list '"+NameOfWave(FullPeakList)+"' is not the right size"
		return $""
	elseif (WaveExists(FullPeakList1))
		if (DimSize(FullPeakList1,1)<11)
			DoAlert 0, "Full peak list 1'"+NameOfWave(FullPeakList1)+"' is the wrong size"
			return $""
		endif
	elseif (WaveExists(FullPeakList2))
		if (DimSize(FullPeakList2,1)<11)
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
	maxSpots = min(maxSpots,N0+N1+N2)
	if (maxSpots<1)										// nothing to do
		return $""
	endif
	N0 = min(maxSpots,N0)
	N1 = min(maxSpots,N1)
	N2 = min(maxSpots,N2)

	STRUCT microGeometry geo							// note, dd and yc are reset from wave note below if it exists
	if (FillGeometryStructDefault(geo,alert=1))//fill the geometry structure with default values
		return $""
	endif
	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))
		DoAlert 0, "no lattice structure found, did you forget to set it?"
		return $""
	endif

	Variable N=0											// number of valid kf's set in kfs[][3]
//	Make/N=(maxSpots,4)/D/FREE kfs					// will hold the result from all detectors, normalized kf's & dNum
	Make/N=(maxSpots,3)/D/FREE kfs					// will hold the result from all detectors, normalized kf's
	SetDimLabel 1,0,Qx,kfs	;	SetDimLabel 1,1,Qy,kfs	;	SetDimLabel 1,2,Qz,kfs
//	SetDimLabel 1,3,dNum,kfs

	String wnotePeak=note(FullPeakList), wnote
	Variable dNum=limit(detectorNumFromID(geo, StringByKey("detectorID",wnotePeak,"=")),0,MAX_Ndetectors)
	Variable startx,groupx, starty,groupy			// ROI of the original image
	startx = NumberByKey("startx",wnotePeak,"=")
	groupx = NumberByKey("groupx",wnotePeak,"=")
	starty = NumberByKey("starty",wnotePeak,"=")
	groupy = NumberByKey("groupy",wnotePeak,"=")
	startx = numtype(startx) ? FIRST_PIXEL : startx
	groupx = numtype(groupx) ? 1 : groupx
	starty = numtype(starty) ? FIRST_PIXEL : starty
	groupy = numtype(groupy) ? 1 : groupy
	Variable depth = NumberByKey("depth",wnotePeak,"="), i

	Make/N=(N0,2)/D/FREE pxpy
	pxpy[][0] = (startx-FIRST_PIXEL) + groupx*FullPeakList[p][0] + (groupx-1)/2		// change to un-binned pixels
	pxpy[][1] = (starty-FIRST_PIXEL) + groupy*FullPeakList[p][1] + (groupy-1)/2		// pixels are still zero based
	if (WaveExists(DeltaPixelCorrect))
		Wave kfs0 = pixel2kfVEC(geo.d[dNum],pxpy,depth=depth,DeltaPixelCorrect=DeltaPixelCorrect)	// normalized kf's
	else
		Wave kfs0 = pixel2kfVEC(geo.d[dNum],pxpy,depth=depth)
	endif
	MatrixOP/FREE badRows = greater(0.5,sumRows(MagSqr(magSqr(kfs0))))	// valid if |kf|==1, so must be >0.5
	if (sum(badRows)<1)									// no bad rows, trasfer all values of kfs0
		kfs[0,N0-1][0,2] = kfs0[p][q]
	else														// only transfer the valid rows
		for (i=0;i<N0;i+=1)
			kfs[N+i][0,2] = kfs0[i][q]
		endfor
		N0 -= sum(badRows)								// reduce N0 so that it is only the valid rows
	endif
//	kfs[0,N0-1][3] = dNum
	N += N0
	WaveClear kfs0

	if (WaveExists(FullPeakList1) && N<maxSpots)
		wnote = note(FullPeakList1)
		dNum = limit(detectorNumFromID(geo, StringByKey("detectorID",wnote,"=")),0,MAX_Ndetectors)
		startx = NumberByKey("startx",wnote,"=")
		groupx = NumberByKey("groupx",wnote,"=")
		starty = NumberByKey("starty",wnote,"=")
		groupy = NumberByKey("groupy",wnote,"=")
		startx = numtype(startx) ? FIRST_PIXEL : startx
		groupx = numtype(groupx) ? 1 : groupx
		starty = numtype(starty) ? FIRST_PIXEL : starty
		groupy = numtype(groupy) ? 1 : groupy
		depth = NumberByKey("depth",wnote,"=")

		N1 = min(N1,maxSpots-N)							// ensure that N1 fits in what is left, can add at most maxSpots-N
		Redimension/N=(N1,2) pxpy
		pxpy[][0] = (startx-FIRST_PIXEL) + groupx*FullPeakList1[p][0] + (groupx-1)/2		// change to un-binned pixels
		pxpy[][1] = (starty-FIRST_PIXEL) + groupy*FullPeakList1[p][1] + (groupy-1)/2		// pixels are still zero based
		if (WaveExists(DeltaPixelCorrect))
			Wave kfs1 = pixel2kfVEC(geo.d[dNum],pxpy,depth=depth,DeltaPixelCorrect=DeltaPixelCorrect)	// normalized kf's
		else
			Wave kfs1 = pixel2kfVEC(geo.d[dNum],pxpy,depth=depth)
		endif
		MatrixOP/FREE badRows = greater(0.5,sumRows(MagSqr(magSqr(kfs1))))	// valid if |kf|==1, so must be >0.5

		if (sum(badRows)<1)									// no bad rows, trasfer all values of kfs1
			kfs[N,N+N1-1][0,2] = kfs1[p-N][q]
		else														// only transfer the valid rows
			for (i=0;i<N1;i+=1)
				kfs[N+i][0,2] = kfs1[i][q]
			endfor
			N1 -= sum(badRows)								// reduce N1 so that it is only the valid rows
		endif
//		kfs[N,N+N1-1][3] = dNum
		N += N1
		WaveClear kfs1
	endif

	if (WaveExists(FullPeakList2) && N<maxSpots)
		wnote = note(FullPeakList2)
		dNum = limit(detectorNumFromID(geo, StringByKey("detectorID",wnote,"=")),0,MAX_Ndetectors)
		startx = NumberByKey("startx",wnote,"=")
		groupx = NumberByKey("groupx",wnote,"=")
		starty = NumberByKey("starty",wnote,"=")
		groupy = NumberByKey("groupy",wnote,"=")
		startx = numtype(startx) ? FIRST_PIXEL : startx
		groupx = numtype(groupx) ? 1 : groupx
		starty = numtype(starty) ? FIRST_PIXEL : starty
		groupy = numtype(groupy) ? 1 : groupy
		depth = NumberByKey("depth",wnote,"=")

		N2 = min(N2,maxSpots-N)							// ensure that N2 fits in what is left, can add at most maxSpots-N
		Redimension/N=(N2,2) pxpy
		pxpy[][0] = (startx-FIRST_PIXEL) + groupx*FullPeakList2[p][0] + (groupx-1)/2		// change to un-binned pixels
		pxpy[][1] = (starty-FIRST_PIXEL) + groupy*FullPeakList2[p][1] + (groupy-1)/2		// pixels are still zero based
		if (WaveExists(DeltaPixelCorrect))
			Wave kfs2 = pixel2kfVEC(geo.d[dNum],pxpy,depth=depth,DeltaPixelCorrect=DeltaPixelCorrect)	// normalized kf's
		else
			Wave kfs2 = pixel2kfVEC(geo.d[dNum],pxpy,depth=depth)
		endif
		MatrixOP/FREE badRows = greater(0.5,sumRows(MagSqr(magSqr(kfs2))))	// valid if |kf|==1, so must be >0.5
		if (sum(badRows)<1)									// no bad rows, trasfer all values of kfs2
			kfs[N,N+N2-1][0,2] = kfs2[p-N][q]
		else														// only transfer the valid rows
			for (i=0;i<N2;i+=1)
				kfs[N+i][0,2] = kfs2[i][q]
			endfor
			N2 -= sum(badRows)								// reduce N2 so that it is only the valid rows
		endif
//		kfs[N,N+N2-1][3] = dNum
		N += N2
		WaveClear kfs2
	endif
	Redimension/N=(N,-1) kfs								// N is now the actual number of valid spots in kfs
	WaveClear pxpy, badRows

	Variable val
	String str, wnotekfs
	if (StringMatch(StringByKey("WaveClass",wnotePeak,"="),"FittedPeakListTest*"))
		wnotekfs = "WaveClass=kfsMeasuredTest;"
	else
		wnotekfs = "WaveClass=kfsMeasured;"
	endif
	wnotekfs = ReplaceStringByKey("structureDesc",wnotekfs,xtal.desc,"=")
	sprintf str,"{ %g, %g, %g, %g, %g, %g }",xtal.a,xtal.b,xtal.c,xtal.alpha,xtal.beta,xtal.gam
	wnotekfs = ReplaceStringByKey("latticeParameters",wnotekfs,str,"=")
	wnotekfs = ReplaceStringByKey("lengthUnit",wnotekfs,"nm","=")
	wnotekfs = ReplaceNumberByKey("SpaceGroup",wnotekfs,xtal.SpaceGroup,"=")
	if (strlen(xtal.SpaceGroupID))
		wnote = ReplaceStringByKey("SpaceGroupID",wnote,xtal.SpaceGroupID,"=")
	endif
	if (numtype(xtal.SpaceGroupIDnum))
		wnote = ReplaceNumberByKey("SpaceGroupIDnum",wnote,xtal.SpaceGroupIDnum,"=")
	endif
	for (i=0;i<xtal.N;i+=1)
		sprintf str,"{%s,%g,%g,%g,%g}",xtal.atom[i].name,xtal.atom[i].x,xtal.atom[i].y,xtal.atom[i].z,xtal.atom[i].occ
		wnotekfs = ReplaceStringByKey("AtomDesctiption"+num2istr(i+1),wnotekfs,str,"=")
	endfor

	// copy a bunch of key=value pairs from wnotePeaks to wnote
	String commonKeys="imageFileName;exposure;dateExposed;rawIgorImage;fittedIgorImage;fractionBkg;"
	commonKeys += "minSpotSeparation;minPeakWidth;maxPeakWidth;totalPeakIntensity;totalIntensity;"
	commonKeys += "startx;groupx;endx;starty;groupy;endy;depth;"
//	commonKeys += "detectorID;"
	wnotekfs = MergeKeywordLists(wnotekfs,wnotePeak,0,"=",";",keys=commonKeys)
	wnotekfs = ReplaceStringByKey("peakListWave",wnotekfs,GetWavesDataFolder(FullPeakList,2),"=")

	Note/K kfs, wnotekfs
	return kfs
End
//
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
	elseif (DimSize(FullPeakList,1)<11)
		DoAlert 0, "the passed full peak list '"+NameOfWave(FullPeakList)+"' is not the right size"
		return $""
	elseif (WaveExists(FullPeakList1))
		if (DimSize(FullPeakList1,1)<11)
			DoAlert 0, "Full peak list 1'"+NameOfWave(FullPeakList1)+"' is the wrong size"
			return $""
		endif
	elseif (WaveExists(FullPeakList2))
		if (DimSize(FullPeakList2,1)<11)
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
	if (FillGeometryStructDefault(geo,alert=1))//fill the geometry structure with default values
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

	Variable dNum=limit(detectorNumFromID(geo, StringByKey("detectorID",wnotePeak,"=")),0,MAX_Ndetectors)
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

	Variable i, px,py
	for (i=0,N=0;i<min(maxSpots,N0);i+=1)
		px = (startx-FIRST_PIXEL) + groupx*FullPeakList[i][0] + (groupx-1)/2		// change to un-binned pixels
		py = (starty-FIRST_PIXEL) + groupy*FullPeakList[i][1] + (groupy-1)/2		// pixels are still zero based
		pixel2q(geo.d[dNum],px,py,qhat,depth=depth)	// was in Wenge-coord system (OLD) or BeamLine(New)
		if (norm(qhat)>0)											// check for a valid Q
			GhatMeasured[N][0,2] = qhat[q]
			GhatMeasured[N][3] = FullPeakList[i][10]	// the integral
			GhatMeasured[N][4] = dNum							// detector number
			N += 1
		endif
	endfor

	if (WaveExists(FullPeakList1) && N<maxSpots)
		wnote = note(FullPeakList1)
		dNum = limit(detectorNumFromID(geo, StringByKey("detectorID",wnote,"=")),0,MAX_Ndetectors)
		startx = NumberByKey("startx",wnote,"=")
		groupx = NumberByKey("groupx",wnote,"=")
		starty = NumberByKey("starty",wnote,"=")
		groupy = NumberByKey("groupy",wnote,"=")
		startx = numtype(startx) ? FIRST_PIXEL : startx
		groupx = numtype(groupx) ? 1 : groupx
		starty = numtype(starty) ? FIRST_PIXEL : starty
		groupy = numtype(groupy) ? 1 : groupy
		depth = NumberByKey("depth",wnote,"=")
		for (i=0; i<N1 && N<maxSpots; i+=1)
			px = (startx-FIRST_PIXEL) + groupx*FullPeakList1[i][0] + (groupx-1)/2		// change to un-binned pixels
			py = (starty-FIRST_PIXEL) + groupy*FullPeakList1[i][1] + (groupy-1)/2		// pixels are still zero based
			pixel2q(geo.d[dNum],px,py,qhat,depth=depth)	// was in Wenge-coord system (OLD) or BeamLine(New)
			if (norm(qhat)>0)											// check for a valid Q
				GhatMeasured[N][0,2] = qhat[q]
				GhatMeasured[N][3] = FullPeakList1[i][10]	// the integral
				GhatMeasured[N][4] = dNum							// detector number
				N += 1
			endif
		endfor
	endif

	if (WaveExists(FullPeakList2) && N<maxSpots)
		wnote = note(FullPeakList2)
		dNum = limit(detectorNumFromID(geo, StringByKey("detectorID",wnote,"=")),0,MAX_Ndetectors)
		startx = NumberByKey("startx",wnote,"=")
		groupx = NumberByKey("groupx",wnote,"=")
		starty = NumberByKey("starty",wnote,"=")
		groupy = NumberByKey("groupy",wnote,"=")
		startx = numtype(startx) ? FIRST_PIXEL : startx
		groupx = numtype(groupx) ? 1 : groupx
		starty = numtype(starty) ? FIRST_PIXEL : starty
		groupy = numtype(groupy) ? 1 : groupy
		depth = NumberByKey("depth",wnote,"=")
		for (i=0; i<N2 && N<maxSpots; i+=1)
			px = (startx-FIRST_PIXEL) + groupx*FullPeakList2[i][0] + (groupx-1)/2		// change to un-binned pixels
			py = (starty-FIRST_PIXEL) + groupy*FullPeakList2[i][1] + (groupy-1)/2		// pixels are still zero based
			pixel2q(geo.d[dNum],px,py,qhat,depth=depth)	// was in Wenge-coord system (OLD) or BeamLine(New)
			if (norm(qhat)>0)											// check for a valid Q
				GhatMeasured[N][0,2] = qhat[q]
				GhatMeasured[N][3] = FullPeakList2[i][10]	// the integral
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
	if (strlen(xtal.SpaceGroupID))
		wnote = ReplaceStringByKey("SpaceGroupID",wnote,xtal.SpaceGroupID,"=")
	endif
	if (numtype(xtal.SpaceGroupIDnum))
		wnote = ReplaceNumberByKey("SpaceGroupIDnum",wnote,xtal.SpaceGroupIDnum,"=")
	endif
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
	MultiThread summs[] = summs_Set(p,inVolume,PairRotations,thresh2)

	Variable topSum1 = WaveMax(summs)-1	// this ensures that NumTopRotations > 0
	MatrixOP/FREE isFrequentRot = greater(summs,topSum1)	// flags pairs having more than topSum1 matches

//if (0)
//	Variable NumTopRotations = sum(isFrequentRot)				// number of rotations that occur > topSum1 times
//	MatrixOP/FREE rotAvg = sumCols(PairRotations*colRepeat(isFrequentRot,3))
//	rotAxis = rotAvg / NumTopRotations		// the result being returned
//	//printf "thresh = %g¡,  summs = [%g, %g]   NumTopRotations = %g,   topSum1 = %g\r",thresh*180/PI,WaveMin(summs),WaveMax(summs),NumTopRotations,topSum1
//else
	// find position of largest point in summs, not the average of largest points
	WaveStats/M=1/Q summs
	Make/N=3/D/FREE rotMostLikely = PairRotations[V_maxloc][p]
#if defined(ZONE_TESTING) || defined(QS_TESTING)
	Variable indexOfMax=V_maxLoc					// index of a point at the max
	Note/K rotAxis, ReplaceNumberByKey("indexOfMax",note(rotAxis),indexOfMax,"=")
#endif

	// sameRot is flags identifying all rotations within thresh of rotMostLikely
	MatrixOP/FREE sameRot = greater(thresh2,sumRows(magSqr(PairRotations-rowRepeat(rotMostLikely,N))))
	Variable NumSameRotations = sum(sameRot)// number of rotations that occur at rotMostLikely
	MatrixOP/FREE rotAvg = sumCols(PairRotations*colRepeat(sameRot,3))
	rotAxis = rotAvg / NumSameRotations		// the result being returned
//endif

#if defined(ZONE_TESTING) || defined(QS_TESTING)
	if (DataFolderExists(":GizmoWaves"))
		//	Duplicate/O summs, summsView
		Duplicate/O PairRotations, :GizmoWaves:PairRotationsViewSize
		Wave PairRotationsViewSize = :GizmoWaves:PairRotationsViewSize
		Redimension/S  PairRotationsViewSize
		Variable sizeMin = N>10000 ? 1 : 2
		PairRotationsViewSize = limit(round(summs[p]/1),sizeMin,25)
		Variable sizeMax = limit(WaveMax(PairRotationsViewSize)+6,6,25)
		PairRotationsViewSize[indexOfMax][] = sizeMax	// make closest point bigger

		Duplicate/O PairRotationsViewSize, :GizmoWaves:PairRotationsViewRGBA
		Wave PairRotationsViewRGBA = :GizmoWaves:PairRotationsViewRGBA
		Redimension/N=(-1,4) PairRotationsViewRGBA
		PairRotationsViewRGBA = 0
		PairRotationsViewRGBA[][3] = inVolume[p] ? 1 : 0.3
		PairRotationsViewRGBA[][1,2] = isFrequentRot[p] ? 1 : 0

		Make/FREE red={1,0,0,1}
		PairRotationsViewRGBA[indexOfMax][] = red[q]	// make closest point red
	endif
#endif

	return topSum1+1								// largest number of parirs in one rotaiton
End
//
ThreadSafe Function summs_Set(i,inVolume,PairRotations,thresh2)
	Variable i
	Wave inVolume
	Wave PairRotations
	Variable thresh2
	if (!inVolume[i])				// PairRotations[i][p] is within radius of center
		return 0
	endif
	Variable N=DimSize(PairRotations,0)
	Make/N=3/D/FREE rot=PairRotations[i][p]	// a single rotation vector from PairRotations[N][3]
	// sumi[0] is the number of rotations in PairRotations[N][3] that are within thresh of rot
	MatrixOP/FREE sumi = sum(greater(thresh2,sumRows(magSqr(PairRotations-rowRepeat(rot,N)))))
	return sumi[0]					// number of rotations in PairRotations[all][3] within thresh of PairRotations[i][3]
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
	Make/N=(nPairs)/FREE mm						// really, just use mm to iterate with MultiThread

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
	Wave g0A,g1A			// measured g's from data, assuming that they are all normalized
	Wave q0A,q1A			// possible q's from given lattice, assuming that they are all normalized
	Wave rotAxisA

	Make/N=3/D/FREE g0=g0A[m][p], g1=g1A[m][p]
	Make/N=3/D/FREE q0=q0A[m][p], q1=q1A[m][p]


if (0)
	normalize(g0)
	normalize(g1)
	normalize(q0)
	normalize(q1)
endif

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


Static Function hkPreferPointsAtDetector(xtal,kfHats,[kin])
	// used to determine the meaning of hklPrefer
	// returns 1 if hklPrefer points at detector center (forward scattering)
	// returns 0 if hklPrefer is the q-vector that Diffracts to the detector center (top detector)
	STRUCT crystalStructure &xtal
	Wave kfHats						// contains all of the measured kf^s
	Wave kin							// optional incident wave direction (defaults to 001)

	Make/N=3/D/FREE ki={0,0,1}								// use this for the incident beam direction
	if (!ParamIsDefault(kin))									// one was passed, use it
		ki = kin[p]
		normalize(ki)
	endif

	Variable N=DimSize(kfHats,0)
	MatrixOP/FREE kf0 = Normalize(sumCols(kfHats))	// assume that average kf points to center of detector
	Redimension/N=3 kf0
	if (MatrixDot(ki,kf0)>0.99)								// kf0 & ki are within 8¡, forward scattering
		return 1
	endif

	Cross kf0, ki
	Wave W_Cross=W_Cross
	Duplicate/FREE W_Cross, a				// a[3] is the axis perp to ki, kf, & qhat0
	KillWaves/Z W_Cross
	normalize(a)
	// now make kfPerp the components kf perpendicular to the a[3], and normlized kfPerp
	MatrixOP/FREE kfPerp = NormalizeRows( kfHats - rowRepeat(a,N)*colRepeat((kfHats x a),3) )

	// alpha & beta are angles in the plane of diffraction, this plane contains ki, kf0
	// alpha is the angle between kf0 and the most distant kf (cone angle from detector center detector edge)
	// beta is the angle from ki to kf closest to ki
	MatrixOP/FREE cosAlpha = minVal( kfPerp x kf0 )			// min { kfPerp ¥ kf0 } = cos(alpha)
	MatrixOP/FREE cosBeta = maxVal( kfPerp x ki )			// max { kfPerp ¥ ki } = cos(beta)
	Variable atDetector = (cosAlpha[0] < cosBeta[0])

	//	if beta < alpha return 1		meaning: incident beam is closer to the detector than the size of the detector
	//	if alpha <= beta return 0		meaning: detector is far from the incident beam 
	return atDetector
End
//Static Function hkPreferPointsAtDetector(xtal,Ghats,[kin])
//	// used to determine the meaning of hklPrefer
//	// returns 1 if hklPrefer points at detector center (forward scattering)
//	// returns 0 if hklPrefer is the q-vector that Diffracts to the detector center (top detector)
//	STRUCT crystalStructure &xtal
//	Wave Ghats						// contains all of the measured G^s
//	Wave kin							// optional incident wave direction (defaults to 001)
//
//	Make/N=3/D/FREE ki={0,0,1}								// use this for the incident beam direction
//	if (!ParamIsDefault(kin))									// one was passed, use it
//		ki = kin[p]
//		normalize(ki)
//	endif
//
//	Variable N=DimSize(Ghats,0)
//	MatrixOP/FREE qhat0 = Normalize(sumCols(Ghats))	// assume that average Ghat is Q-vector to center of detector
//	Redimension/N=3 qhat0
//	Variable dot = MatrixDot(ki,qhat0)
//	Make/N=3/D/FREE kf0 = (ki - 2*dot*qhat0)			//		kf^ = ki^ - 2*(ki^ . q^)*q^,   kf0 is normalized
//	if (MatrixDot(ki,kf0)>0.99)								// kf0 & ki are within 8¡, forward scattering
//		return 1
//	endif
//
//	// form a kf for every Ghat
//	MatrixOP/FREE kf = (rowRepeat(ki,N) - 2*colRepeat((Ghats x ki),3)*Ghats)	//		kf^ = ki^ - 2*(ki^ . q^)*q^
//
//	Cross kf0, ki
//	Wave W_Cross=W_Cross
//	Duplicate/FREE W_Cross, a				// a[3] is the axis perp to ki, kf, & qhat0
//	KillWaves/Z W_Cross
//	normalize(a)
//	// now make kfPerp the components kf perpendicular to the a[3], and normlized kfPerp
//	MatrixOP/FREE kfPerp = NormalizeRows( kf - rowRepeat(a,N)*colRepeat((kf x a),3) )
//
//	// alpha & beta are angles in the plane of diffraction, this plane contains ki, kf0 & qhat0
//	// alpha is the angle between kf0 and the most distant kf (cone angle from detector center detector edge)
//	// beta is the angle from ki to kf closest to ki
//	MatrixOP/FREE cosAlpha = minVal( kfPerp x kf0 )			// min { kfPerp ¥ kf0 } = cos(alpha)
//	MatrixOP/FREE cosBeta = maxVal( kfPerp x ki )			// max { kfPerp ¥ ki } = cos(beta)
//	Variable atDetector = (cosAlpha[0] < cosBeta[0])
//
//	//	if beta < alpha return 1		meaning: incident beam is closer to the detector than the size of the detector
//	//	if alpha <= beta return 0		meaning: detector is far from the incident beam 
//	return atDetector
//End
//Function test_hkPreferPointsAtDetector()
//	Wave FullPeakList=FullPeakListOrange
//	STRUCT crystalStructure xtal
//	if (FillCrystalStructDefault(xtal))
//		DoAlert 0, "no lattice structure found, did you forget to set it?"
//		return 1
//	endif
//	Wave GhatsMeasured=IndexingInternal#FullPeakList2Ghats(250,FullPeakList)	// convert peaks to a Qlist+intens+dNum
//	Variable NG0=DimSize(GhatsMeasured,0)
//	Make/N=(NG0,3)/FREE GhatsOnly=GhatsMeasured[p][q]
//	print hkPreferPointsAtDetector(xtal,GhatsOnly)
//End


ThreadSafe Static Function goodNessOfPattern(FullPeakIndexed,ipat)
	Wave FullPeakIndexed
	Variable ipat			// pattern number usually 0

	Variable Nspots=DimSize(FullPeakIndexed,0)
	Make/N=(Nspots)/D/FREE intensity=FullPeakIndexed[p][6]
	WaveStats/M=1/Q intensity
	Variable goodness = V_sum>0 ? V_sum*(Nspots*Nspots) : Nspots
	return goodness
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
	intens = numtype(intens) || intens<0 ? 0 : intens
	return intens
End


ThreadSafe Static Function/WAVE hklOrderedList(hklMax,[hklLo,hklHi,Nmax])
	// If you need an ordered list of hkl's, this is much faster then getting fancy
	// and trying to order the loops.  LatticeSym#incrementIndex(h) is not too bad, but 
	// using an another loop around the hkl loops slows you down a factor of 33 for hklMax=16
	// for hklMax=20 this is faster by x88, for hklMax=30 faster by x420
	// If you only want the smallest Nmax hkl, then do a Redimension/N=(Nmax,-1) hkls
	// Note, if you only want k&l to vary and h fixed, you can say:
	// Wave hkls = hklOrderedList(hklMax,hLo=1,hHi=1)

	Variable hklMax
	Wave hklLo, hklHi								// optional range of h,k,l (if not passed or NaN, use hklMax)

	Variable Nmax									// max number of hkls to make, default is 1e6 (and 1e6 takes ~1minute)
	Nmax = ParamIsDefault(Nmax) || numtype(Nmax) || round(Nmax)<1 ? 1e6 : round(Nmax)
	hklMax = round(hklMax)

	Variable hLo=NaN,hHi=NaN, kLo=NaN,kHi=NaN, lLo=NaN,lHi=NaN		// can be passed in hklLo & hklHi to overRide hklMax
	if (WaveExists(hklLo))
		hLo = hklLo[0]
		kLo = hklLo[1]
		lLo = hklLo[2]
	endif
	if (WaveExists(hklHi))
		hHi = hklHi[0]
		kHi = hklHi[1]
		lHi = hklHi[2]
	endif
	hLo = numtype(hLo) ? -hklMax : round(hLo)	// use ±hklMax if value is invalid
	kLo = numtype(kLo) ? -hklMax : round(kLo)
	lLo = numtype(lLo) ? -hklMax : round(lLo)
	hHi = numtype(hHi) ?  hklMax : round(hHi)
	kHi = numtype(kHi) ?  hklMax : round(kHi)
	lHi = numtype(lHi) ?  hklMax : round(lHi)
	Variable Nh=hHi-hLo+1, Nk=kHi-kLo+1, Nl=lHi-lLo+1
	Variable N=Nh*Nk*Nl
	if (numtype(Nh+Nk+Nl) || Nh<1 || Nk<1 || Nl<1 || N>Nmax)
		return $""
	endif

	Make/N=(N,3)/I/FREE hkls
	Make/N=3/I/FREE hkl
	Variable h,k,l, i
	for (h=hHi,i=0; h>=hLo; h-=1)
		hkl[0] = h
		for (k=kHi; k>=kLo; k-=1)
			hkl[1] = k
			for (l=lHi; l>=lLo; l-=1)
				hkl[2] = l
				hkls[i][] = hkl[q]
				i += 1
			endfor
		endfor
	endfor

	MatrixOP/FREE len2 = sumRows(magSqr(hkls))
	Make/N=(N)/I/FREE index
	MakeIndex len2, index
	WaveClear len2
	Duplicate/FREE hkls, hklsOut
	hklsOut = hkls[index[p]][q]
	return hklsOut
End
//Function check_hklOrderedList(hklMax)
//	Variable hklMax
//
//	Variable h,k,l, N, hklM
//	Variable Nmax = (2*hklMax+1)^3
//	printf "should generate (2*%d+1)^3 = %d^3 = %d hkl's\r",hklMax,2*hklMax+1,Nmax
//
//	Make/N=(Nmax,3)/FREE/I hkls=1e9
//	Make/N=3/I/FREE hkl
//	timeIncrement(fresh=1)
//	for (hklM=0,N=0; hklM<=hklMax; hklM+=1)
//		for (l=0; abs(l)<=hklM; l=LatticeSym#incrementIndex(l))
//			hkl[2] = l
//			for (k=0; k<=hklM; k=LatticeSym#incrementIndex(k))
//				hkl[1] = k
//				for (h=0; h<=hklM; h=LatticeSym#incrementIndex(h))
//					hkl[0] = h
//					MatrixOP/FREE absMax = maxVal(abs(hkl))
//					if (absMax[0]<hklM)
//						continue
//					endif
//					hkls[N][] = hkl[q]
//					N += 1
//				endfor
//			endfor
//		endfor
//	endfor
//	print "N = ",N,"  incrementIndex with outer loop ",timeIncrement()
//	Duplicate/O hkls hklsView
//
//	timeIncrement(fresh=1)
//	for (l=0,N=0; abs(l)<=hklMax; l=LatticeSym#incrementIndex(l))
//		hkl[2] = l
//		for (k=0; k<=hklMax; k=LatticeSym#incrementIndex(k))
//			hkl[1] = k
//			for (h=0; h<=hklMax; h=LatticeSym#incrementIndex(h))
//				hkl[0] = h
//				hkls[N][] = hkl[q]
//				N += 1
//			endfor
//		endfor
//	endfor
//	print "N = ",N,"  only incrementIndex ",timeIncrement()
//
//	timeIncrement(fresh=1)
//	N = 0
//	for (l=-hklMax; l<=hklMax; l+=1)
//		hkl[2] = l
//		for (k=-hklMax; k<=hklMax; k+=1)
//			hkl[1] = k
//			for (h=-hklMax; h<=hklMax; h+=1)
//				hkl[0] = h
//				N += 1
//			endfor
//		endfor
//	endfor
//	print "N = ",N,"  assign hkl,N   ",timeIncrement()
//
//	timeIncrement(fresh=1)
//	N = 0
//	for (l=-hklMax; l<=hklMax; l+=1)
//		for (k=-hklMax; k<=hklMax; k+=1)
//			for (h=-hklMax; h<=hklMax; h+=1)
//				N += 1
//			endfor
//		endfor
//	endfor
//	print "N = ",N,"  just looping   ",timeIncrement()
//
//	timeIncrement(fresh=1)
//	Wave hkls = hklOrderedList(16)
//	print "N = ",N,"  hklOrderedList ",timeIncrement()
//	Duplicate/O hkls hkls2View
//End



ThreadSafe Static Function AngularSpanOfVectors(hats)
	Wave hats
	MatrixOP/FREE minDots = minVal(hats x hats^t)
	return acos(limit(minDots[0],-1,1))*180/PI
End


ThreadSafe Static Function isNewDirection(existingAxes,axis,cosTol,[anti])
	//	returns True if axis does not point toward an existing directions in existingAxes
	// if anti==1, then axis is also not anti-parallel to existingAxes
	Wave existingAxes			// array(N,3) list of existing axes
	Wave axis
	Variable cosTol			// tolerance on dot:  cosTol = cos(Æangle),  cos(0.1¡) = 0.999998476913288
	Variable anti				// is true then parallel and anti-parallel return TRUE (default is FALSE)
	if (norm(axis)==0)
		return 0					// skip the (000)
	endif
	anti = ParamIsDefault(anti) ? 0 : anti

	if (anti)
		MatrixOP/FREE dots = abs(existingAxes x axis)
	else
		MatrixOP/FREE dots = existingAxes x axis
	endif
	Variable dotMax = WaveMax(dots)
	return (dotMax < cosTol || numtype(dotMax))	// if existingAxes filled with NaN then anything works
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


Function/T timeIncrement([fresh,longTime,seconds])
	Variable fresh							// reset all ticks to zero, default is NO restart
	Variable longTime						// lets you set what is a long time (default is 0.2 sec)
	Variable seconds						// returns total execution time in seconds, not the regular string
	fresh = ParamIsDefault(fresh) || numtype(fresh) ? 0 : fresh
	longTime = ParamIsDefault(longTime) || numtype(longTime)==2 ? 0.2 : longTime
	longTime = longTime<=0 ? Inf : longTime		// 0 or negative times eliminate long time flag
	seconds = ParamIsDefault(seconds) || numtype(seconds) ? 0 : seconds

	Variable now=stopMSTimer(-2)
	if (exists("root:Packages:timeIncrement:tick0")!=2 || fresh)
		NewDataFolder/O root:Packages
		NewDataFolder/O root:Packages:timeIncrement
		Variable/G root:Packages:timeIncrement:tick0=now, root:Packages:timeIncrement:tick1=now
		Variable/G root:Packages:timeIncrement:maxDelta=0
		return ""
	endif

	NVAR tick0=root:Packages:timeIncrement:tick0, tick1=root:Packages:timeIncrement:tick1
	NVAR maxDelta=root:Packages:timeIncrement:maxDelta
	String out=""
	if (tick1 > tick0)
		maxDelta = max( maxDelta, (now-tick1)*1e-6 )
		sprintf out, "elapsed time = %.3f s  (Æ = %.3f s)", (now-tick0)*1e-6, (now-tick1)*1e-6
		if ((now-tick1)*1e-6 > longTime)
			out += "\t\t********* This is a long step. *********"
		endif
	else
		sprintf out, "elapsed time = %.3f s", (now-tick0)*1e-6
	endif
	tick1 = now

	if (seconds)
		return num2str( (now-tick0)*1e-6 )
	endif
	return out
End

//  ================================== End of Utility ==================================  //
//  ====================================================================================  //




//  ====================================================================================  //
//  ================================ Start of Test Data ================================  //

// make the wave GhatsMeasured, GhatsMeasured are the result of finding peaks on a detector.
// This wave is a good set of test data for testing all of this stuff.
Function/WAVE MakeSimulatedTestPattern(Nreq,[lattice,angle,axis,tol,keVmax,dNum,printIt])
	Variable Nreq						// number of test spots to make
	String lattice						// "C","M1",M2"
	Variable angle						// rotation angle from lattice about axis (degree)
	String axis							// string representation of axis for angle rotation
	Variable tol						// angle (degree), usually 0.1¡, reject hkl's that differ by less than tol
	Variable keVmax					// max energy of a reflection (keV)
	Variable dNum						// detector number
	Variable printIt
	lattice = SelectString(ParamIsDefault(lattice),lattice,"C")
	angle = ParamIsDefault(angle) || numtype(angle) ? 0 : angle
	axis = SelectString(ParamIsDefault(axis),axis,"")
	tol = ParamIsDefault(tol) ? 0.1 : tol
	keVmax = ParamIsDefault(keVmax) || numtype(keVmax) || keVmax<=0 ? 30 : keVmax
	dNum = ParamIsDefault(dNum) || numtype(dNum) || dNum<=0 ? 0 : dNum
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : !(!printIt)

	STRUCT microGeometry geo							// note, dd and yc are reset from wave note below if it exists
	if (FillGeometryStructDefault(geo,alert=1))//fill the geometry structure with default values
		return $""
	endif

	Wave axis3 = str2vec(axis)
	if (Nreq<=0 || numtype(Nreq))
		Nreq = 98
		Prompt Nreq,"Number Ghat's to make"
		Prompt lattice,"Test lattice to construct",popup,"Cubic (C);Monoclinic 1 (M1);Monoclinic 2 (M2)"
		Prompt angle, "rotation angle applied to lattice"
		if (strlen(axis)<1)
			axis="1,2,3"
		endif
		Prompt axis, "axis of test rotation applied to lattice"
		Prompt dNum, "Detector Number", popup, DetectorMenuList(geo)
		DoPrompt "Number of G's",Nreq,lattice,angle,axis, dNum
		if (V_flag)
			return $""
		endif
		dNum -= 1

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
		if (!ParamIsDefault(angle) || !(angle==0))
			printf ", angle=%g, axis=\"%s\"",angle,vec2str(axis3,sep=",",bare=1)
		endif
		if (!ParamIsDefault(tol) || tol!=0.1)
			printf ", tol=%g",tol
		endif
		if (!ParamIsDefault(keVmax) || keVmax!=30)
			printf ", keVmax=%g",keVmax
		endif
		if (!ParamIsDefault(dNum) || dNum!=0)
			printf ", dNum=%g",dNum
		endif
		if (!ParamIsDefault(printIt))
			printf ", printIt=%g",printIt
		endif
		printf ")\r"
	endif
	if (numtype(tol+keVmax) || tol<=0 || keVmax<=0)
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

	String SpaceGroupID=""
	Variable cosTol=cos(tol*PI/180)		// convert degree to dot, this should be slightly less than 1.0
	Make/N=(3,3)/D/FREE direct
	String FullPeakListName, desc=""
	if (StringMatch(lattice,"C"))
		direct = (p==q) * 0.54310206		// simple cubic direct lattice with a = 0.3 nm
		SpaceGroupID = "227:1"				// conventional diamond structure
		desc = "Cubic Test"
		FullPeakListName = "FullPeakListTestCubic"
	elseif (StringMatch(lattice,"M1"))
		direct[0][0] = {0.28, 0, 0.076}	// a   a non-cubic direct lattice (Monoclinic)
		direct[0][1] = {0, 0.3, 0}		// b	(a^c !=90, a^b = c^b = 90,  alpha=gamma=90
		direct[0][2] = {0, 0, 0.29}		// c
		SpaceGroupID = "14:b1"
		desc = "Monoclinic 1 Test"
		FullPeakListName = "FullPeakListTestM1"
	elseif (StringMatch(lattice,"M2"))
		direct[0][0] = {0.28, 0, 0.075}	// a   a non-cubic direct lattice (Monoclinic)
		direct[0][1] = {0, 0.3, 0}		// b	(a^c !=90, a^b = c^b = 90,  alpha=gamma=90
		direct[0][2] = {0, 0, 0.29}		// c
		SpaceGroupID = "14:b1"
		desc = "Monoclinic 2 Test"
		FullPeakListName = "FullPeakListTestM2"
	else
		return $""
	endif

	Variable SpaceGroup = str2num(SpaceGroupID)
	Variable SpaceGroupIDnum = SpaceGroupID2num(SpaceGroupID)
	Variable Vc = MatrixDet(direct)
	LatticeSym#SetSymOpsForSpaceGroup(SpaceGroupID)

	MatrixOP/FREE/O recip = 2*PI * (Inv(direct))^t
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
		printf "lattice constants: %.9g, %.9g, %.9g,   %g¡, %g¡, %g¡,   SG=%s\r",LC[0],LC[1],LC[2],LC[3],LC[4],LC[5],SpaceGroupID
	endif
	STRUCT crystalStructure xtal
	xtal.a = LC[0] ;			xtal.b = LC[1] ;		xtal.c = LC[2]
	xtal.alpha = LC[3] ;	xtal.beta = LC[4] ;	xtal.gam = LC[5]
	xtal.SpaceGroup = SpaceGroup
	xtal.SpaceGroupID = SpaceGroupID
	xtal.SpaceGroupIDnum = SpaceGroupIDnum
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

	if (dNum<0 || dNum>= geo.Ndetectors)
		DoAlert 0, "dNum="+num2str(dNum)+" is out of the allowed range [0,"+num2str(geo.Ndetectors -1)+"]"
		return $""
	elseif (!(geo.d[dNum].used))
		DoAlert 0, "Detector dNum="+num2str(dNum)+" is not used"
		return $""
	endif


	Variable Nx = geo.d[dNum].Nx, Ny = geo.d[dNum].Ny, px,py
	Make/N=3/D/FREE ki={0,0,1}, qhat0
	XYZ2pixel(geo.d[dNum],ki,px,py)
	Variable atDetector = (px>=0 && px<Nx && py>=0 && py<Ny)	// if True ki hits detector (forward scattering), False for top detector

#if defined(ZONE_TESTING) || defined(QS_TESTING) || defined(ZONE_QS_TESTING)
	Make/N=(Nreq,3)/D/O GhatsMeasured=0, GhatsMeasuredHKL=0 // this holds the result
#else
	Make/N=(Nreq,3)/D/FREE GhatsMeasured=0, GhatsMeasuredHKL=0
#endif
	Make/N=(Nreq,3)/D/FREE kfsFound=0

	Make/N=3/D/FREE hkl
	Make/N=(Nreq,12)/D/O $FullPeakListName/WAVE=FullPeakListTest = NaN
	Variable keV, intens, intensMax=0, hklMax=20
	Variable/C pz

	Wave hkls = hklOrderedList(hklMax)
	Variable Nhkls=DimSize(hkls,0)
	Redimension/N=(Nhkls,-1) hkls
	MatrixOP/FREE qvecs = ( recip x hkls^t )^t
	MatrixOP/FREE qhats = NormalizeRows(qvecs)
	MatrixOP/FREE Qmags = sqrt(sumRows(magSqr(qvecs)))
	MatrixOP/FREE kfs = NormalizeRows( (rowRepeat(ki,Nhkls) - 2*colRepeat((qhats x ki),3)*qhats) )	//		kf^ = ki^ - 2*(ki^ . q^)*q^

	Make/N=3/D/FREE qhat, qvec, hkl, kf
	Variable h,k,l, N							// N is the actual number of unique spots saved
	for (i=0,N=0; i<Nhkls && N<Nreq; i+=1)
		qhat = qhats[i][p]
		pz = q2pixel(geo.d[dNum],qhat)
		px = real(pz)
		py = imag(pz)
		if (px<0 || px>=Nx || py<0 || py>=Ny || numtype(px+py))	// misses the detector
			continue
		endif
		keV = -hc*Qmags[i] / (4*PI*MatrixDot(ki,qhat))	// Q = 4*PI*sin(theta)/lambda
		kf = kfs[i][p]
		hkl = hkls[i][p]

		if (isNewDirection(kfsFound,kf,cosTol) && keV<=keVmax && allowedHKL(hkl[0],hkl[1],hkl[2],xtal))
			qvec = qvecs[i][p]
			FullPeakListTest[N][0] = px
			FullPeakListTest[N][1] = py
			intens = genericIntensity(xtal,qvec,hkl,keV)
			intensMax = max(intens,intensMax)
			FullPeakListTest[N][10] = intens

			GhatsMeasured[N][] = qhat[q]		// a new direction, save it...
			GhatsMeasuredHKL[N][] = hkl[q]	// hkl associated with each test Ghat
			kfsFound[N][] = kf[q]
			N += 1
		endif
	endfor
	Redimension/N=(N,-1) GhatsMeasured,GhatsMeasuredHKL		// trim to actual size
	Redimension/N=(N,-1) FullPeakListTest

	// find the central hkl
	MatrixOP/FREE GhatAvg = Normalize(sumCols(GhatsMeasured))
	Redimension/N=3 GhatAvg
	MatrixOP/FREE/O hkl = Normalize(Inv(recip) x GhatAvg)
	Variable small = smallestNonZeroValue(hkl,tol=1e-8)
	hkl *= 12/small
	removeCommonFactors(hkl)
	FullPeakListTest[][10] /= intensMax
	FullPeakListTest[][11] = FullPeakListTest[p][10]

	String wnote="waveClass=FittedPeakListTest;"
	wnote = ReplaceStringByKey("direct",wnote,directStr,"=")
	wnote = ReplaceStringByKey("recip",wnote,recipStr,"=")
	wnote = ReplaceNumberByKey("SpaceGroup",wnote,SpaceGroup,"=")
	wnote = ReplaceStringByKey("SpaceGroupID",wnote,SpaceGroupID,"=")
	wnote = ReplaceNumberByKey("SpaceGroupIDnum",wnote,SpaceGroupIDnum,"=")
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
	wnote = ReplaceNumberByKey("xdim",wnote,Nx,"=")
	wnote = ReplaceNumberByKey("ydim",wnote,Ny,"=")
	wnote = ReplaceNumberByKey("xDimDet",wnote,Nx,"=")
	wnote = ReplaceNumberByKey("yDimDet",wnote,Ny,"=")
	wnote = ReplaceStringByKey("detectorID",wnote,geo.d[dNum].detectorID,"=")
	wnote = ReplaceStringByKey("detectorModel",wnote,"XRD0820","=")
	wnote = ReplaceStringByKey("beamLine",wnote,"34ID-E","=")
	wnote = ReplaceNumberByKey("exposure",wnote,1,"=")
	wnote = ReplaceStringByKey("MonoMode",wnote,"white slitted","=")
	Note/K FullPeakListTest,wnote
	SetDimLabel 1,0,x0,FullPeakListTest				;	SetDimLabel 1,1,y0,FullPeakListTest
	SetDimLabel 1,2,x0Err,FullPeakListTest			;	SetDimLabel 1,3,y0Err,FullPeakListTest
	SetDimLabel 1,4,fwx,FullPeakListTest				;	SetDimLabel 1,5,fwy,FullPeakListTest
	SetDimLabel 1,6,fwxErr,FullPeakListTest			;	SetDimLabel 1,7,fwyErr,FullPeakListTest
	SetDimLabel 1,8,correlation,FullPeakListTest	;	SetDimLabel 1,9,correlationErr,FullPeakListTest
	SetDimLabel 1,10,area,FullPeakListTest			;	SetDimLabel 1,11,amp,FullPeakListTest
	if (printIt)
		printf "maximum hkl went to (%d %d %d)\r",hklMax,hklMax,hklMax
		printf "Made %d simulated peaks, stored in '%s'\r",N,NameOfWave(FullPeakListTest)
	endif
	return FullPeakListTest
End
//// make the wave GhatsMeasured, GhatsMeasured are the result of finding peaks on a detector.
//// This wave is a good set of test data for testing all of this stuff.
//Function/WAVE MakeSimulatedTestPattern(Nreq,[cone,lattice,angle,axis,tol,keVmax,dNum,printIt])
//	Variable Nreq						// number of test spots to make
//	Variable cone						// cone angle (degree)
//	String lattice						// "C","M1",M2"
//	Variable angle						// rotation angle from lattice about axis (degree)
//	String axis							// string representation of axis for angle rotation
//	Variable tol						// angle (degree), usually 0.1¡, reject hkl's that differ by less than tol
//	Variable keVmax					// max energy of a reflection (keV)
//	Variable dNum						// detector number
//	Variable printIt
//	cone = ParamIsDefault(cone) || numtype(cone) || cone<=0 ? 20 : cone
//	lattice = SelectString(ParamIsDefault(lattice),lattice,"C")
//	angle = ParamIsDefault(angle) || numtype(angle) ? 0 : angle
//	axis = SelectString(ParamIsDefault(axis),axis,"")
//	tol = ParamIsDefault(tol) ? 0.1 : tol
//	keVmax = ParamIsDefault(keVmax) || numtype(keVmax) || keVmax<=0 ? 30 : keVmax
//	dNum = ParamIsDefault(dNum) || numtype(dNum) || dNum<=0 ? 0 : dNum
//	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : !(!printIt)
//
//	Wave axis3 = str2vec(axis)
//	if (Nreq<=0 || numtype(Nreq))
//		Nreq = 98
//		Prompt Nreq,"Number Ghat's to make"
//		Prompt cone,"cone angle (0 gives all)"
//		Prompt lattice,"Test lattice to construct",popup,"Cubic (C);Monoclinic 1 (M1);Monoclinic 2 (M2)"
//		Prompt angle, "rotation angle applied to lattice"
//		if (strlen(axis)<1)
//			axis="1,2,3"
//		endif
//		Prompt axis, "axis of test rotation applied to lattice"
//		Prompt dNum, "Detector Number", popup, "Orange;Yellow;Purple;MARCCD"
//		DoPrompt "Number of G's",Nreq,cone,lattice,angle,axis, dNum
//		if (V_flag)
//			return $""
//		endif
//		dNum -= 1
//
//		Variable i=strsearch(lattice,"(",0)
//		lattice = lattice[i+1,Inf]
//		i = strsearch(lattice,")",0)
//		lattice = lattice[0,i-1]
//
//		angle = numtype(angle) ? 0 : angle
//		if (angle)
//			Wave axis3 = str2vec(axis)
//		endif
//		printIt = 1
//	endif
//	if (Nreq<=0 || numtype(Nreq))
//		return $""
//	endif
//	if (printIt)
//		printf "MakeSimulatedTestPattern(%g",Nreq
//		if (!ParamIsDefault(lattice))
//			printf ", lattice=\"%s\"",lattice
//		endif
//		if (!ParamIsDefault(cone) || cone>0)
//			printf ", cone=%g",cone
//		endif
//		if (!ParamIsDefault(angle) || !(angle==0))
//			printf ", angle=%g, axis=\"%s\"",angle,vec2str(axis3,sep=",",bare=1)
//		endif
//		if (!ParamIsDefault(tol) || tol!=0.1)
//			printf ", tol=%g",tol
//		endif
//		if (!ParamIsDefault(keVmax) || keVmax!=30)
//			printf ", keVmax=%g",keVmax
//		endif
//		if (!ParamIsDefault(dNum) || dNum!=0)
//			printf ", dNum=%g",dNum
//		endif
//		if (!ParamIsDefault(printIt))
//			printf ", printIt=%g",printIt
//		endif
//		printf ")\r"
//	endif
//	if (numtype(cone+tol) || tol<=0 || abs(cone>180))
//		return $""
//	elseif (angle)
//		if (!WaveExists(axis3))
//			return $""
//		elseif (!(norm(axis3)>0))
//			return $""
//		elseif (numpnts(axis3)!=3)
//			return $""
//		endif
//	endif
//
//	Variable cosTol=cos(tol*PI/180)		// convert degree to dot, this should be slightly less than 1.0
//	Make/N=(3,3)/D/FREE direct
//	String FullPeakListName, desc=""
//	if (StringMatch(lattice,"C"))
//		direct = (p==q) * 0.54310206		// simple cubic direct lattice with a = 0.3 nm
////		Variable SpaceGroup = 224			// [195,230]
//		Variable SpaceGroup = 227			// [195,230], diamond structure
//		desc = "Cubic Test"
//		FullPeakListName = "FullPeakListTestCubic"
//	elseif (StringMatch(lattice,"M1"))
//		direct[0][0] = {0.28, 0, 0.076}	// a   a non-cubic direct lattice (Monoclinic)
//		direct[0][1] = {0, 0.3, 0}		// b	(a^c !=90, a^b = c^b = 90,  alpha=gamma=90
//		direct[0][2] = {0, 0, 0.29}		// c
//		SpaceGroup = 14						// [3,15]
//		desc = "Monoclinic 1 Test"
//		FullPeakListName = "FullPeakListTestM1"
//	elseif (StringMatch(lattice,"M2"))
//		direct[0][0] = {0.28, 0, 0.075}	// a   a non-cubic direct lattice (Monoclinic)
//		direct[0][1] = {0, 0.3, 0}		// b	(a^c !=90, a^b = c^b = 90,  alpha=gamma=90
//		direct[0][2] = {0, 0, 0.29}		// c
//		SpaceGroup = 14						// [3,15]
//		desc = "Monoclinic 2 Test"
//		FullPeakListName = "FullPeakListTestM2"
//	else
//		return $""
//	endif
//
//	Variable Vc = MatrixDet(direct)
//	LatticeSym#SetSymOpsForSpaceGroup(SpaceGroup)
//
//	MatrixOP/FREE/O recip = 2*PI * (Inv(direct))^t
//	recip = recip==0 ? 0 : recip
//	String directStr=encodeMatAsStr(direct), recipStr=encodeMatAsStr(recip)
//
//	if (angle)
//		MatrixOP/FREE axisNorm = Normalize(axis3)
//		Duplicate/FREE axisNorm, axisTemp
//		axisTemp *= angle*PI/180						// length of axisTemp is rot in radians
//		Wave rot = RotMatAboutAxisOnly(axisTemp)
//		MatrixOP/FREE/O recip = rot x recip
//		MatrixOP/FREE/O direct = rot x direct
//		direct = direct == 0 ? 0 : direct			// get rid of -0
//		recip = recip == 0 ? 0 : recip
//		Make/N=3/D/FREE scaledRotAxis = axisNorm*(angle*PI/180)
//		if (printIt)
//			printf "test rotation is %g¡ about the (%s)  =  %s\r",angle,vec2str(axis3,sep=" ",bare=1),vec2str(scaledRotAxis)
//		endif
//	endif
//	if (!CheckDirectRecip(direct,recip))
//		print "Invalid computation of direct & recip"
//		Abort "Invalid computation of direct & recip"
//		return $""
//	endif
//
//	Wave LC = direct2LatticeConstants(direct)
//	if (printIt)
//		printf "lattice constants: %.9g, %.9g, %.9g,   %g¡, %g¡, %g¡,   SG=%G\r",LC[0],LC[1],LC[2],LC[3],LC[4],LC[5],SpaceGroup
//	endif
//	STRUCT crystalStructure xtal
//	xtal.a = LC[0] ;			xtal.b = LC[1] ;		xtal.c = LC[2]
//	xtal.alpha = LC[3] ;	xtal.beta = LC[4] ;	xtal.gam = LC[5]
//	xtal.SpaceGroup = SpaceGroup
//	xtal.desc = desc
//	if (!isValidLatticeConstants(xtal)	)
//		print "ERROR -- Lattice Constants are NOT valid"
//		return $""
//	endif
//	xtal.Unconventional00 = NaN
//	xtal.sourceFile = ""
//	LatticeSym#CleanOutCrystalStructure(xtal)
//	LatticeSym#ForceLatticeToStructure(xtal)
//	MakeSymmetryOps(xtal)					// make a wave with the symmetry operation
//	UpdateCrystalStructureDefaults(xtal)
//
//	Make/N=3/D/FREE qUp=0
//	Variable cosCone=-2						// -2 includes everything
//	if (cone>0)
//		qUp = {0,1,-1}							// direction of q-vector that will diffract straight up
//		normalize(qUp)
//		cosCone = cos(cone*PI/180)
//	endif
//
//	STRUCT microGeometry geo							// note, dd and yc are reset from wave note below if it exists
//	if (FillGeometryStructDefault(geo))			//fill the geometry structure with default values
//		DoAlert 0, "no geometry structure found, did you forget to set it?"
//		return $""
//	endif
//	if (dNum<0 || dNum>= geo.Ndetectors)
//		DoAlert 0, "dNum="+num2str(dNum)+" is out of the allowed range [0,"+num2str(geo.Ndetectors -1)+"]"
//		return $""
//	elseif (!(geo.d[dNum].used))
//		DoAlert 0, "Detector dNum="+num2str(dNum)+" is not used"
//		return $""
//	endif
//
//#if defined(ZONE_TESTING) || defined(QS_TESTING)
//	Make/N=(Nreq,3)/D/O GhatsMeasured=0, GhatsMeasuredHKL=0 // this holds the result
//#else
//	Make/N=(Nreq,3)/D/FREE GhatsMeasured=0, GhatsMeasuredHKL=0
//#endif
//	Make/N=(Nreq,3)/D/FREE Ghat3=0		// only need this for the isNewDirection() test
//	Make/N=3/D/FREE hkl, ki={0,0,1}
//	Make/N=(Nreq,11)/D/O $FullPeakListName/WAVE=FullPeakListTest = NaN
//	Variable keV, intens, intensMax=0
//	Variable/C pz
//	Variable Nx = geo.d[dNum].Nx, Ny = geo.d[dNum].Ny
//	Variable hklMax=20, px,py
//
//	Wave hkls = hklOrderedList(hklMax)
//	Variable Nhkls=DimSize(hkls,0)
//	Redimension/N=(Nhkls,-1) hkls
//	MatrixOP/FREE qvecs = ( recip x hkls^t )^t
//	MatrixOP/FREE qhats = NormalizeRows(qvecs)
//	MatrixOP/FREE Qmags = sqrt(sumRows(magSqr(qvecs)))
//
//	Make/N=3/D/FREE qhat, qvec, hkl
//	Variable h,k,l, N							// N is the actual number of unique spots saved
//	for (i=0; i<Nhkls && N<Nreq; i+=1)
//		hkl = hkls[i][p]
//		qhat = qhats[i][p]
//		qvec = qvecs[i][p]
//		pz = q2pixel(geo.d[dNum],qhat)
//		px = real(pz)
//		py = imag(pz)
//		if (px<0 || px>=Nx || py<0 || py>=Ny || numtype(px+py))	// misses the detector
//			continue
//		endif
//		keV = -hc*Qmags[i] / (4*PI*MatrixDot(ki,qhat))	// Q = 4*PI*sin(theta)/lambda
//if (numtype(keV)==0 && allowedHKL(hkl[0],hkl[1],hkl[2],xtal))
//if (!(MatrixDot(qhat,qUp)>cosCone))
//print "* MatrixDot(qhat,qUp)>cosCone  ",MatrixDot(qhat,qUp),cosCone
//else
//print "*****"
//endif
//Debugger
//endif
//
//		if (isNewDirection(Ghat3,qhat,cosTol) && MatrixDot(qhat,qUp)>cosCone && keV<=keVmax && allowedHKL(hkl[0],hkl[1],hkl[2],xtal))
//			FullPeakListTest[N][0] = px
//			FullPeakListTest[N][1] = py
//			intens = genericIntensity(xtal,qvec,hkl,keV)
//			intensMax = max(intens,intensMax)
//			FullPeakListTest[N][10] = intens
//
//			GhatsMeasured[N][] = qhat[q]		// a new direction, save it...
//			GhatsMeasuredHKL[N][] = hkl[q]	// hkl associated with each test Ghat
//			Ghat3[N][] = qhat[q]				// Ghat3 only needed for isNewDirection()
//			N += 1
//		endif
//	endfor
//	Redimension/N=(N,-1) GhatsMeasured,GhatsMeasuredHKL, Ghat3	// trim to actual size
//	Redimension/N=(N,-1) FullPeakListTest
//
//	// find the central hkl
//	MatrixOP/FREE GhatAvg = Normalize(sumCols(Ghat3))
//	Redimension/N=3 GhatAvg
//	MatrixOP/FREE/O hkl = Normalize(Inv(recip) x GhatAvg)
//	Variable small = smallestNonZeroValue(hkl,tol=1e-8)
//	hkl *= 12/small
//	removeCommonFactors(hkl)
//	FullPeakListTest[][10] /= intensMax
//
//	String wnote="waveClass=FittedPeakListTest;"
//	wnote = ReplaceStringByKey("direct",wnote,directStr,"=")
//	wnote = ReplaceStringByKey("recip",wnote,recipStr,"=")
//	wnote = ReplaceNumberByKey("SpaceGroup",wnote,SpaceGroup,"=")
//	wnote = ReplaceNumberByKey("cone",wnote,cone,"=")
//	wnote = ReplaceNumberByKey("tol",wnote,tol,"=")
//	wnote = ReplaceStringByKey("hklCenter",wnote,vec2str(hkl,bare=1,sep=","),"=")
//	if (angle)
//		wnote = ReplaceNumberByKey("rotAngle",wnote,angle,"=")
//		wnote = ReplaceStringByKey("rotAxis",wnote,vec2str(axis3,bare=1,sep=","),"=")
//	endif
//	Note/K GhatsMeasured, wnote
//	wnote = ReplaceNumberByKey("startx",wnote,0,"=")
//	wnote = ReplaceNumberByKey("endx",wnote,2047,"=")
//	wnote = ReplaceNumberByKey("groupx",wnote,1,"=")
//	wnote = ReplaceNumberByKey("starty",wnote,0,"=")
//	wnote = ReplaceNumberByKey("endy",wnote,2047,"=")
//	wnote = ReplaceNumberByKey("groupy",wnote,1,"=")
//	wnote = ReplaceNumberByKey("xdim",wnote,Nx,"=")
//	wnote = ReplaceNumberByKey("ydim",wnote,Ny,"=")
//	wnote = ReplaceNumberByKey("xDimDet",wnote,Nx,"=")
//	wnote = ReplaceNumberByKey("yDimDet",wnote,Ny,"=")
//	wnote = ReplaceStringByKey("detectorID",wnote,"PE1621 723-3335","=")
//	wnote = ReplaceStringByKey("detectorModel",wnote,"XRD0820","=")
//	wnote = ReplaceStringByKey("beamLine",wnote,"34ID-E","=")
//	wnote = ReplaceNumberByKey("exposure",wnote,1,"=")
//	wnote = ReplaceStringByKey("MonoMode",wnote,"white slitted","=")
//	Note/K FullPeakListTest,wnote
//	SetDimLabel 1,0,x0,FullPeakListTest		;	SetDimLabel 1,1,y0,FullPeakListTest
//	SetDimLabel 1,2,x0Err,FullPeakListTest	;	SetDimLabel 1,3,y0Err,FullPeakListTest
//	SetDimLabel 1,4,fwx,FullPeakListTest		;	SetDimLabel 1,5,fwy,FullPeakListTest
//	SetDimLabel 1,6,fwxErr,FullPeakListTest	;	SetDimLabel 1,7,fwyErr,FullPeakListTest
//	SetDimLabel 1,8,correlation,FullPeakListTest ;	SetDimLabel 1,9,correlationErr,FullPeakListTest
//	SetDimLabel 1,10,area,FullPeakListTest
//	if (printIt)
//		printf "maximum hkl went to (%d %d %d)\r",hklMax,hklMax,hklMax
//		printf "Made %d simulated peaks, stored in '%s'\r",N,NameOfWave(FullPeakListTest)
//	endif
//	return FullPeakListTest
//End

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
	String win = StringFromList(0,WindowsWithWave(FullPeakList,1))
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
	STRUCT microGeometry geo
	if (FillGeometryStructDefault(geo,alert=1))//fill the geometry structure with default values
		return 1
	endif
	Variable dNum=limit(detectorNumFromID(geo, StringByKey("detectorID",wnote,"=")),0,MAX_Ndetectors)
	Variable startx=NumberByKey("startx",wnote,"="), starty=NumberByKey("starty",wnote,"=")
	Variable groupx=NumberByKey("groupx",wnote,"="), groupy=NumberByKey("groupy",wnote,"=")
	Variable Nx=NumberByKey("xdim",wnote,"="), Ny=NumberByKey("ydim",wnote,"=")
	SetAxis/R left Nx,starty
	SetAxis bottom startx,Ny

	Variable px,py
	Make/N=3/D/FREE xyz, xyzAbs				// xyz points to detector center
	pixel2XYZ(geo.d[dNum],0.5*(geo.d[dNum].Nx - 1),0.5*(geo.d[dNum].Ny - 1),xyz)
	xyzAbs = abs(xyz)
	WaveStats/M=1/Q xyzAbs						// need to find which value in xyz has the largest magnitude
	xyz = p==V_maxloc ? sign(xyz[p]) : 0	// xyz is now the closest orthogonal direction that points to the detector
	XYZ2pixel(geo.d[dNum],xyz,px,py)		// set pixel where xyz hits detector
	if (px>=0 && px <= geo.d[dNum].Nx && py>=0 && py <= geo.d[dNum].Ny)	// (px,py) IS on detector, draw cross
		Variable NxFW=Nx/10, NyFW=Ny/10							// size of cross
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

	Variable cosTol=cos(tol*PI/180)			// convert degree to dot, this should be slightly less than 1.0
	Make/N=(3,3)/D/FREE direct
	String desc="", SpaceGroupID=""

	if (StringMatch(lattice,"C"))
		direct = (p==q) * 0.54310206			// simple cubic direct lattice with a = 0.3 nm
		SpaceGroupID = "224:1"
SpaceGroupID = "227:1"
		desc = "Cubic Test"
	elseif (StringMatch(lattice,"M1"))
		direct[0][0] = {0.28, 0, 0.076}		// a   a non-cubic direct lattice (Monoclinic)
		direct[0][1] = {0, 0.3, 0}			// b	(a^c !=90, a^b = c^b = 90,  alpha=gamma=90
		direct[0][2] = {0, 0, 0.29}			// c
		SpaceGroupID = "14:b1"
		desc = "Monoclinic 1 Test"
	elseif (StringMatch(lattice,"M2"))
		direct[0][0] = {0.28, 0, 0.075}		// a   a non-cubic direct lattice (Monoclinic)
		direct[0][1] = {0, 0.3, 0}			// b	(a^c !=90, a^b = c^b = 90,  alpha=gamma=90
		direct[0][2] = {0, 0, 0.29}			// c
		SpaceGroupID = "14:b1"
		desc = "Monoclinic Test"
	endif
	Variable SpaceGroup = str2num(SpaceGroupID)
	Variable SpaceGroupIDnum = SpaceGroupID2num(SpaceGroupID)
	Variable Vc = MatrixDet(direct)
	LatticeSym#SetSymOpsForSpaceGroup(SpaceGroupID)

	MatrixOP/FREE/O recip   = 2*PI * (Inv(direct))^t
	recip = recip==0 ? 0 : recip
	String directStr=encodeMatAsStr(direct), recipStr=encodeMatAsStr(recip)

	if (angle)
		MatrixOP/FREE axisNorm = Normalize(axis3)
		Duplicate/FREE axisNorm, axisTemp
		axisTemp *= angle*PI/180				// length of axisTemp is rot in radians
		Wave rot = RotMatAboutAxisOnly(axisTemp)
		MatrixOP/FREE/O recip = rot x recip
		MatrixOP/FREE/O direct = rot x direct
		direct = direct == 0 ? 0 : direct	// get rid of -0
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
		printf "lattice constants: %.9g, %.9g, %.9g,   %g¡, %g¡, %g¡,   SG=%s\r",LC[0],LC[1],LC[2],LC[3],LC[4],LC[5],SpaceGroupID
	endif
	STRUCT crystalStructure xtal
	xtal.a = LC[0] ;			xtal.b = LC[1] ;		xtal.c = LC[2]
	xtal.alpha = LC[3] ;	xtal.beta = LC[4] ;	xtal.gam = LC[5]
	xtal.SpaceGroup = SpaceGroup
	xtal.SpaceGroupID = SpaceGroupID
	xtal.SpaceGroupIDnum = SpaceGroupIDnum
	xtal.desc = desc
	if (!isValidLatticeConstants(xtal)	)
		print "ERROR -- Lattice Constants are NOT valid"
		return $""
	endif
	xtal.Unconventional00 = NaN
	xtal.sourceFile = ""
	LatticeSym#CleanOutCrystalStructure(xtal)
	LatticeSym#ForceLatticeToStructure(xtal)
	MakeSymmetryOps(xtal)						// make a wave with the symmetry operation
	UpdateCrystalStructureDefaults(xtal)

	Make/N=3/D/FREE qUp=0
	Variable cosCone=-2							// -2 includes everything
	if (cone>0)
		qUp = {0,1,-1}								// direction of q-vector that will diffract straight up
		normalize(qUp)
		cosCone = cos(cone*PI/180)
	endif

	Make/N=(Nreq,3)/D/FREE Ghat3=0			// only need this for the isNewDirection() test
	Make/N=(Nreq,3)/D/O GhatsMeasured=0	// this holds the result
	Make/N=(Nreq,3)/D/O GhatsMeasuredHKL=0	// this holds the result
	Make/N=3/D/FREE hkl
	Variable h,k,l, N	=0							// N is the actual number of unique spots saved
	Variable hklMax, hklMaxMax=20
	for (hklMax=0,N=0; hklMax<=hklMaxMax && N<Nreq; hklMax+=1)
		for (l=0; l<=hklMax && N<Nreq; l=LatticeSym#incrementIndex(l))
			hkl[2] = l
			for (k=0; k<=hklMax && N<Nreq; k=LatticeSym#incrementIndex(k))
				hkl[1] = k
				for (h=0; h<=hklMax && N<Nreq; h=LatticeSym#incrementIndex(h))
					hkl[0] = h
					MatrixOP/FREE/O qhkl = Normalize(recip x hkl)
					if (isNewDirection(Ghat3,qhkl,cosTol) && MatrixDot(qhkl,qUp)>cosCone)
						GhatsMeasured[N][] = qhkl[q]		// a new direction, save it...
						GhatsMeasuredHKL[N][] = hkl[q]	// hkl associated with each test Ghat
						Ghat3[N][] = qhkl[q]				// Ghat3 only needed for isNewDirection()
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
	wnote = ReplaceStringByKey("SpaceGroupID",wnote,SpaceGroupID,"=")
	wnote = ReplaceNumberByKey("SpaceGroupIDnum",wnote,SpaceGroupIDnum,"=")
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
