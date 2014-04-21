#pragma rtGlobals=1		// Use modern global access method.
#pragma ModuleName=detectorCalibration
#pragma version = 0.81
#include "microGeometryN", version>=1.63
#include "ImageDisplayScaling", version>=1.97

Static Constant hc = 1.239841857					// keV-nm
// #define OLD_ORIENTATION_METHOD

//#define OptimizeWireAxis				// if defined, optimize wire.axis[], otherwise  wire.R[]
//		for Optimization of Wire Geometry from a wire scan, see section  "=== Start of Wire Calibration Fit ==="

Menu LaueGoMainMenuName
	SubMenu "Optimize - Calibrate"
		MenuItemIfWaveClassExists("Enter Measured Energies for Calibration","IndexedPeakList*",""),EnterMeasuredEnergies("")
		MenuItemIfWaveClassExists("Set Calibration Input Data for Optimize","IndexedPeakList*",""),MakeCalibData(NaN)
		MenuItemIfWaveClassExists("Optimize All 3 Detector Geometrys","DetectorCalibrationList*",""),OptimizeAll($"",$"",$"")
		"-"
		MenuItemIfWaveClassExists("  Graph Calibration Spots on 1 Detector","DetectorCalibrationList*",""),GraphCalibrationDetector($"")
		MenuItemIfWaveClassExists("  Graph Calibration Spots on All Detectors","DetectorCalibrationList*",""),GraphAllCalibrationDetector()
		MenuItemIfWaveClassExists("  Table of Calibration Data","DetectorCalibrationList*",""),TableCalibrationData($"")
		"Write Detector values to EPICS...",WriteDetectorGeo2EPICS(NaN)
		"-"
		"(testing"
		"  Make Fake Calibration Data",MakeFakeCalibration3Detectors(NaN)
		"  test many Optimizations",testManyOptimize(NaN,NaN)
		"  test first of many Optimizations",testOneOptimize(NaN)
		"-"
		"(Optimize Geometry from WireScan Peak TurnOffs"
		MenuItemIfWaveClassExists("  Make Table of Turn Offs","PixelIntensities*","DIMS:2"), MakeWireFittingTable($"")
		MenuItemIfWaveClassExists("  Optimize Wire Geometry","PixelTurnOffs","DIMS:2"), OptimizeWireGeometry($"",NaN)
		MenuItemIfWaveClassExists("  Show Table of Wire Scan Turn Offs","PixelTurnOffs","DIMS:2"), TableOfPixelTurnOffs($"")
		"-"
		"Calculate Wire Origin",calcWireOrigin(NaN,NaN,NaN,NaN,NaN)
		"For Detector in Direct Beam, find P from pixel",findPxPyFromPixel(NaN,NaN,NaN)
	End
End


///	noteStr=ReplaceStringByKey("waveClass",noteStr,"DetectorCalibrationList"+SelectString(dNum,"0",""),"=")


// To calibrate, find the two vectors R and P that best fit the data in CalibrationList[][]
// I do not think tOptimizeDetectorCalibrationhat you can also optimize the sizeX & sizeY, they must be constant since I really only measure angles


// ==============================================================================================
// ==================================== Start of Calibration ====================================

// ==============================================================================================
// ==================================== Start of Wire Origin ====================================

// Takes old calibration values for wire and converts to wire origin
//	example: calcWireOrigin(-6454,-6019,4350,880,1092)
Function/T calcWireOrigin(H0,Hyc,F0,px,py)
	Variable H0									// H where wire is centered on incident beam
	Variable Hyc								// H were wire is directly above the sample origin
	Variable F0									// wire is at constant F
	Variable px,py								// pixel used to measure Hyc

	Variable printIt = strlen(GetRTStackInfo(2))==0
	if (numtype(H0+Hyc+F0+px+py))
		Prompt H0,"H0, PM500 wire H when wire splits incident beam"
		Prompt Hyc,"Hyc, PM500 wire H when wire nearly over sample origin"
		Prompt F0,"F, PM500 wire F for both points"
		Prompt px,"pixel X (used to get Hyc)"
		Prompt py,"pixel Y (used to get Hyc)"
		DoPrompt "Wire Coordinates",H0,Hyc,F0,px,py
		if (V_flag)
			return ""
		endif
		printf "calcWireOrigin(%g, %g, %g, %g, %g)\r",H0,Hyc,F0,px,py
		printIt = 1
	endif
	if (numtype(H0+Hyc+F0))
		return ""
	endif

	Variable X0=0, Y0,Z0
	Variable Xyc=0,Yyc,Zyc
	Y0 = HF2Y(H0,F0)										// position of wire when leading edge cuts incident beam
	Z0 = HF2Z(H0,F0)
	Yyc= HF2Y(Hyc,F0)									// position of wire centered in diffracted beam
	Zyc = HF2Z(Hyc,F0)									// F is constant, so no Fyc, just keep using F0

	STRUCT microGeometry g
	FillGeometryStructDefault(g)
	g.wire.origin[0] = 0									// This is what we are trying to find, set local copy to zero
	g.wire.origin[1] = 0
	g.wire.origin[2] = 0

	Make/N=3/O/D xyzPixel_getWireOrigin, xyz_yc_getWireOrigin, xyz0_getWireOrigin, origin_getWireOrigin
	Wave xyzPixel = xyzPixel_getWireOrigin, xyz_yc=xyz_yc_getWireOrigin, xyz0=xyz0_getWireOrigin,  origin=origin_getWireOrigin
	Make/N=(3,3)/O/D rhoW_getWireOrigin
	Wave rhoW=rhoW_getWireOrigin
	rhoW[0][0]=g.wire.R00;	rhoW[0][1]=g.wire.R01;	rhoW[0][2]=g.wire.R02
	rhoW[1][0]=g.wire.R10;	rhoW[1][1]=g.wire.R11;	rhoW[1][2]=g.wire.R12
	rhoW[2][0]=g.wire.R20;	rhoW[2][1]=g.wire.R21;	rhoW[2][2]=g.wire.R22
	xyz_yc = {Xyc,Yyc,Zyc}

	pixel2XYZ(g.d[0],px,py,xyzPixel)						// convert pixel position to the beam line coordinate system

	Variable depth,i
	origin = 0
	do
		g.wire.origin[0] += origin[0]						// This is what we are trying to find, set local copy to zero
		g.wire.origin[1] += origin[1]
		g.wire.origin[2] += origin[2]

		depth = PixelxyzWire2depth(g,xyzPixel,xyz_yc,1)// returns depth (µm) for leading edge
		depth += PixelxyzWire2depth(g,xyzPixel,xyz_yc,0)// returns depth (µm) for trailing edge
		depth /= 2											// This is Z of the beam from origin

		xyz0 = {X0,Y0,Z0}
		PM500X2toBeamLineX2(g.wire,xyz0)				// beam goes through xyz0 along (001)

		origin[0] = {xyz0[0], xyz0[1], depth}				// in beam line system, just need to rotate back to PM500 system
		MatrixOp/O origin = Inv(rhoW) x origin			// (xyz) = (wire.Rij) x (xyz),   rotate wire position from PM500 coords to beam line coords
		i += 1												// count number of iterations
	while (norm(origin)>1e-9 && i<20)
	if (i>=20)
		DoAlert 0, "ERROR -- Failed to converge after "+num2str(i)+" interations"
		return "nan;nan;nan"
	endif
	origin += g.wire.origin[p]								// Final value

	if (printIt)
		printf "wire at beam	= {%.2f, %.2f, %.2f} (PM500 frame)\r",X0,Y0,Z0
		printf "wire above		= {%.2f, %.2f, %.2f} (PM500 frame)\r",Xyc,Yyc,Zyc
		printf "desired wire Origin	= {%.3f, %.3f, %.3f} (PM500 frame)\r",origin[0],origin[1],origin[2]
	endif
	String str
	sprintf str,"%.3f;%.3f;%.3f",origin[0],origin[1],origin[2]

	// Go back and check the result
	print "\rCheck depth with these new values"
	g.wire.origin[0] = origin[0]							// Set to the computed values
	g.wire.origin[1] = origin[1]
	g.wire.origin[2] = origin[2]

	depth = PixelxyzWire2depth(g,xyzPixel,xyz_yc,1)		// returns depth (µm) for leading edge
	depth += PixelxyzWire2depth(g,xyzPixel,xyz_yc,0)	// returns depth (µm) for trailing edge
	depth /= 2												// This is Z of the beam from origin
	printf "average depth = %g\r",depth
	if (abs(depth)>1e-6)
		DoAlert 0,"ERROR -- Bad calculation, depth not zero, depth="+num2str(depth)
		str = "nan;nan;nan"
	endif
	KillWaves/Z xyz0_getWireOrigin, xyz_yc_getWireOrigin, rhoW_getWireOrigin, origin_getWireOrigin
	Killwaves/Z xyzPixel_getWireOrigin
	return str
End
//// Takes old calibration values for wire and converts to wire origin
////	example: calcWireOrigin(-6454,-5982,-6056,4350,880,1092)
//Function/T calcWireOrigin(H0,Hyc1,Hyc2,F0,px,py)
//	Variable H0									// H where wire is centered on incident beam
//	Variable Hyc1,Hyc2							// H were wire is directly above the sample origin
//	Variable F0									// wire is at constant F
//	Variable px,py								// pixel used to measure Hyc
//
//	Variable printIt = strlen(GetRTStackInfo(2))==0
//	if (numtype(H0+Hyc1+Hyc2+F0+px+py))
//		Prompt H0,"H0, PM500 wire H when wire splits incident beam"
//		Prompt Hyc1,"Hyc Leading Edge, PM500 wire H when wire nearly over sample origin"
//		Prompt Hyc2,"Hyc Trailing Edge, PM500 wire H when wire nearly over sample origin"
//		Prompt F0,"F, PM500 wire F for both points"
//		Prompt px,"pixel X (used to get Hyc)"
//		Prompt py,"pixel Y (used to get Hyc)"
//		DoPrompt "Wire Coordinates",H0,Hyc1,Hyc2,F0,px,py
//		if (V_flag)
//			return ""
//		endif
//		printf "calcWireOrigin(%g, %g, %g, %g, %g, %g)\r",H0,Hyc1,Hyc2,F0,px,py
//		printIt = 1
//	endif
//	if (numtype(H0+Hyc1+Hyc2+F0))
//		return ""
//	endif
//	if (Hyc1>Hyc2)											// in case Hyc was entered backward
//		Variable swap=Hyc1
//		Hyc1 = Hyc2
//		Hyc2 = swap
//	endif
//
//	Variable X0=0, Y0,Z0
//	Variable Xyc=0,Yyc1,Zyc1,Yyc2,Zyc2
//	Y0 = HF2Y(H0,F0)										// position of wire when leading edge cuts incident beam
//	Z0 = HF2Z(H0,F0)
//	Yyc1 = HF2Y(Hyc1,F0)									// position of wire when leading edge cuts diffracted beam
//	Zyc1 = HF2Z(Hyc1,F0)									// F is constant, so no Fyc, just keep using F0
//	Yyc2 = HF2Y(Hyc2,F0)									// position of wire when trailing edge cuts diffracted beam
//	Zyc2 = HF2Z(Hyc2,F0)
//
//	STRUCT microGeometry g
//	FillGeometryStructDefault(g)
//	g.wire.origin[0] = 0									// This is what we are trying to find, set local copy to zero
//	g.wire.origin[1] = 0
//	g.wire.origin[2] = 0
//
//	Make/N=3/O/D xyzPixel_getWireOrigin, xyz_yc_getWireOrigin, xyz0_getWireOrigin, origin_getWireOrigin
//	Wave xyzPixel = xyzPixel_getWireOrigin, xyz_yc=xyz_yc_getWireOrigin, xyz0=xyz0_getWireOrigin,  origin=origin_getWireOrigin
//	Make/N=(3,3)/O/D rhoW_getWireOrigin
//	Wave rhoW=rhoW_getWireOrigin
//	rhoW[0][0]=g.wire.R00;	rhoW[0][1]=g.wire.R01;	rhoW[0][2]=g.wire.R02
//	rhoW[1][0]=g.wire.R10;	rhoW[1][1]=g.wire.R11;	rhoW[1][2]=g.wire.R12
//	rhoW[2][0]=g.wire.R20;	rhoW[2][1]=g.wire.R21;	rhoW[2][2]=g.wire.R22
//
//	pixel2XYZ(g.d[0],px,py,xyzPixel)						// convert pixel position to the beam line coordinate system
//
//	Variable depth, depth1,depth2
//	origin = 0
//	do
//		g.wire.origin[0] += origin[0]						// This is what we are trying to find, set local copy to zero
//		g.wire.origin[1] += origin[1]
//		g.wire.origin[2] += origin[2]
//
//		xyz_yc = {Xyc,Yyc1,Zyc1}
//		depth1 = PixelxyzWire2depth(g,xyzPixel,xyz_yc,1)// returns depth (µm) for leading edge
//		xyz_yc = {Xyc,Yyc2,Zyc2}
//		depth2 = PixelxyzWire2depth(g,xyzPixel,xyz_yc,0)// returns depth (µm) for trailing edge
//		depth = (depth1+depth2)/2							// This is Z of the beam from origin
//
//		xyz0 = {X0,Y0,Z0}
//		PM500X2toBeamLineX2(g.wire,xyz0)				// beam goes through xyz0 along (001)
//
//		origin[0] = {xyz0[0], xyz0[1], depth}				// in beam line system, just need to rotate back to PM500 system
//		MatrixOp/O origin = Inv(rhoW) x origin			// (xyz) = (wire.Rij) x (xyz),   rotate wire position from PM500 coords to beam line coords
//	while (norm(origin)>1e-9)
//
//	origin += g.wire.origin[p]								// Final value
//
//	if (printIt)
//		printf "wire at beam	= {%.2f, %.2f, %.2f} (PM500 frame)\r",X0,Y0,Z0
//		printf "wire above		= {%.2f, %.2f, %.2f} (PM500 frame) Leading Edge\r",Xyc,Yyc1,Zyc1
//		printf "wire above		= {%.2f, %.2f, %.2f} (PM500 frame) Trailing Edge\r",Xyc,Yyc2,Zyc2
//		printf "desired wire Origin	= {%.3f, %.3f, %.3f} (PM500 frame)\r",origin[0],origin[1],origin[2]
//	endif
//	String str
//	sprintf str,"%.3f;%.3f;%.3f",origin[0],origin[1],origin[2]
//
//	// Go back and check the result
//	print "\rCheck depth with these new values"
//	g.wire.origin[0] = origin[0]							// Set to the computed values
//	g.wire.origin[1] = origin[1]
//	g.wire.origin[2] = origin[2]
//
//	xyz_yc = {Xyc,Yyc1,Zyc1}
//	depth1 = PixelxyzWire2depth(g,xyzPixel,xyz_yc,1)	// returns depth (µm) for leading edge
//	xyz_yc = {Xyc,Yyc2,Zyc2}
//	depth2 = PixelxyzWire2depth(g,xyzPixel,xyz_yc,0)	// returns depth (µm) for trailing edge
//	depth = (depth1+depth2)/2								// This is Z of the beam from origin
//	printf "depth from leading edge = %g,  trailing = %g,   average depth = %g\r",depth1,depth2,depth
//	if (abs(depth)>1e-6)
//		DoAlert 0,"Bad calculation, depth not zero, depth="+num2str(depth)
//		str = "nan;nan;nan"
//	endif
//	KillWaves/Z xyz0_getWireOrigin, xyz_yc_getWireOrigin, rhoW_getWireOrigin, origin_getWireOrigin
//	Killwaves/Z xyzPixel_getWireOrigin
//	return str
//End

// ===================================== End of Wire Origin =====================================
// ==============================================================================================



// ==============================================================================================
// ==================================== Start of Optimization ===================================

Function/T MakeCalibData(dNum)					// Make the calibration list needed by OptimizeAll()
	Variable dNum
	dNum = mod(dNum,1) ? NaN : dNum
	dNum = dNum!=limit(dNum,0,2) ? NaN : dNum

	String wlist = reverseList(WaveListClass("FittedPeakList","*","DIMS:2"))
	String peakListName = StringFromList(0,wlist)
	Wave FullPeakList = $peakListName
	Variable getPeakList = (!WaveInClass(FullPeakList,"FittedPeakList") || ItemsInList(wlist)!=1)
	if (numtype(dNum) || getpeakList)
		Prompt dNum,"detector number",popup,"Detector 0, Orange;Detector1, Yellow;Detector 2, Purple"
		Prompt peakListName, "List of Fitted Peaks",popup,wlist
		if (!getpeakList)
			DoPrompt "detector",dNum
			dNum -= 1
		elseif (numtype(dNum)==0)
			DoPrompt "detector",peakListName
		else
			DoPrompt "detector",dNum, peakListName
			dNum -= 1
		endif
		if (V_flag)
			return ""
		endif
		Wave FullPeakList = $peakListName
	endif
	printf "	MakeCalibData(%d),   using FullPeakList='%s'\r",dNum,peakListName
	Variable Npeaks=DimSize(FullPeakList,0), k

	// Find the correct FullPeakIndexed wave, may be a FullPeakIndexedAux for dNum = 1 or 2
	String str, indexList=""
	wlist = reverseList(WaveListClass("IndexedPeakList*","*",""))
	for (k=0;k<ItemsInList(wlist);k+=1)								// find the best choice
		str = StringByKey("peakListWave",note($StringFromList(k,wlist)),"=")
		if (stringmatch(peakListName,ParseFilePath(3,str,":",0,0)))
			indexList = StringFromList(k,wlist)
		endif
	endfor
	Wave FullPeakIndexed=$indexList
	if (!WaveExists(FullPeakIndexed))
		DoAlert 0,"Cannot find the correct FullPeakIndexed wave"
		return ""
	endif
	Variable Nindex=DimSize(FullPeakIndexed,0), m

	// always ask for Ename, even side detectors can have energies, algthough they don't need them
	String Ename=""														// find the list of energies that goes with FullPeakList
	wlist = SelectString(dNum,"","_none_;")+WaveListClass("measuredEnergies","*","MAXCOLS:1,MAXROWS:"+num2istr(Npeaks))
	if (ItemsInList(wlist))
		Prompt Ename,"Measured Energies",popup,wlist
		DoPrompt "Energies",Ename
		if (V_flag)
			return ""
		endif
	endif
	Wave Emeasured=$Ename											// Emeasured should line up with FullPeakList (NOT indexed peaks)
	printf "For detector %d,  using FullPeakIndexed='%s',   with Energies in '%s'\r",dNum,NameOfWave(FullPeakIndexed),NameOfWave(Emeasured)
	if (!WaveExists(Emeasured) && dNum==0)
		DoAlert 0,"Cannot find the correct wave with measured energies"
		return ""
	elseif (WaveExists(Emeasured))
		Duplicate/O Emeasured, $(NameOfWave(Emeasured)+"_dE")	// not really needed, just to look at.
		Wave deltaE = $(NameOfWave(Emeasured)+"_dE")
		SetScale d 0,0,"eV", deltaE
		Note/K deltaE, ReplaceStringByKey("waveClass", note(Emeasured), "measuredEnergies_dE","=")
	endif

	Variable rhox, rhoy, rhoz
	if (dNum==0)															// set sample rotation {rhox,rhoy,rhoz} from the FullPeakIndexed wave
		Wave axisSampleFake = $findSampleAxis(FullPeakIndexed)	// rotation vector for sample
		rhox = axisSampleFake[0]
		rhoy = axisSampleFake[1]
		rhoz = axisSampleFake[2]
		KillWaves/Z axisSampleFake
	else																	// find the sample rotation sample {rhox,rhoy,rhoz} from CalibrationList0
		wlist = reverseList(WaveListClass("DetectorCalibrationList0","*","DIMS:2"))
		String wname=StringFromList(0,wlist)
		if (ItemsInList(wlist)>1)										// ask which one
			Prompt wname,"Calibration List with sample rotation {rhox,rhoy,rhoz}",popup,wlist
			DoPrompt "sample rotation",wname
			if (V_flag)
				wname = ""
			endif
		endif
		if (strlen(wname))												// pre-set rho before asking
			Wave ww=$wname
			wname = note(ww)
			rhox=NumberByKey("rhox",wname,"=")
			rhoy=NumberByKey("rhoy",wname,"=")
			rhoz=NumberByKey("rhoz",wname,"=")
		endif
		Prompt rhox,"rotation axis of sample (X)"
		Prompt rhoy,"rotation axis of sample (Y)"
		Prompt rhoz,"rotation axis of sample (Z)"
		DoPrompt "sample rotation",rhox,rhoy,rhoz
		if (V_flag)
			return ""
		endif
	endif
	printf "	using sample rotation of  {%g, %g, %g}\r",rhox,rhoy,rhoz

	STRUCT microGeometry g
	FillGeometryStructDefault(g)										//fill the geometry structure with current values

	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))									// fill the lattice structure with test values
		DoAlert 0, "no crystal structure found"
		return ""
	endif
	Make/N=(3,3)/O/D recip0_MakeCalib_
	Wave recip0 = recip0_MakeCalib_

	recip0[0][0]=xtal.as0;	recip0[0][1]=xtal.bs0;	recip0[0][2]=xtal.cs0	// reciprocal lattice { a*[], b*[], c*[] }
	recip0[1][0]=xtal.as1;	recip0[1][1]=xtal.bs1;	recip0[1][2]=xtal.cs1
	recip0[2][0]=xtal.as2;	recip0[2][1]=xtal.bs2;	recip0[2][2]=xtal.cs2

	String name = CleanupName("CalibrationList"+ReplaceString("FullPeakList",NameOfWave(FullPeakList),""),0)
	name = name[0,29]+num2istr(dNum)
//	Make/N=(Npeaks,10)/O/D $name=NaN								// list of measured reflections used to calibrate the detector, {px,py,qx,qy,yz,keV,theta}
	Make/N=(Npeaks,12)/O/D $name=NaN								// list of measured reflections used to calibrate the detector, {px,py,qx,qy,yz,keV,theta}
	Wave cList = $name
	SetDimLabel 1,0,px,cList;		SetDimLabel 1,1,py,cList											// measued pixel
	SetDimLabel 1,2,H,cList;			SetDimLabel 1,3,K,cList;			SetDimLabel 1,4,L,cList			// (hkl)
	SetDimLabel 1,5,qx_hkl,cList;	SetDimLabel 1,6,qy_hkl,cList;	SetDimLabel 1,7,qz_hkl,cList	// qvector from (hkl)
	SetDimLabel 1,8,keV,cList														// measured keV
	SetDimLabel 1,9,theta_keV_hkl,cList											// theta calculated from known d(hkl) and measured keV
	SetDimLabel 1,10,calculated_keV,cList											// energy calculated from sample orientation (no measured parameters)
	SetDimLabel 1,11,deltaE_eV,cList												// measured - calculated energy (eV)
	String noteStr = ""
	noteStr=ReplaceStringByKey("waveClass",noteStr,"DetectorCalibrationList"+SelectString(dNum,"0",""),"=")
	noteStr = ReplaceStringByKey("detectorID",noteStr,g.d[dNum].detectorID,"=")// detector ID
	noteStr = ReplaceNumberByKey("Nx",noteStr,g.d[dNum].Nx,"=")				// number of un-binned pixels in whole detector
	noteStr = ReplaceNumberByKey("Ny",noteStr,g.d[dNum].Ny,"=")
	noteStr = ReplaceNumberByKey("sizeX",noteStr,g.d[dNum].sizeX,"=")		// outside size of detector (micron)
	noteStr = ReplaceNumberByKey("sizeY",noteStr,g.d[dNum].sizeY,"=")
	noteStr = ReplaceNumberByKey("rhox",noteStr,rhox,"=")					// rotation vector for sample
	noteStr = ReplaceNumberByKey("rhoy",noteStr,rhoy,"=")
	noteStr = ReplaceNumberByKey("rhoz",noteStr,rhoz,"=")
	Note/K cList, noteStr
	//	CalibrationList[][0,1]	= (px,py) are measured on image
	//	CalibrationList[][2-4]	= hkl are hkl for each measured spot
	//	CalibrationList[][5-7]	= Q-vector, calculated from knowing hkl and the standard reciprocal lattice (NOT rotated from standard orientation)
	//	CalibrationList[][8]	= keV, measured energy of spot
	//	CalibrationList[][9]	= theta (deg), obtained from calculated d(hkl) and measured energy using: lambda = 2d sin(theta)
	//	CalibrationList[][10]	= keV calculated from sample orientation and known lattice (does not use measured spot position)
	//	CalibrationList[][11]	= ÆE (eV),  ( CalibrationList[][8] - CalibrationList[][10] ) *1e3

	Variable tolerance, i									// dpixel tolerance to find indexed spots that match fitted peaks
	for (tolerance=Inf,m=0; m<(Npeaks-1); m+=1)
		for (i=m+1;i<Npeaks;i+=1)
			tolerance = min((FullPeakList[m][0]-FullPeakList[i][0])^2 +(FullPeakList[m][1]-FullPeakList[i][1])^2, tolerance)
		endfor
	endfor
	tolerance = tolerance/4									// in units of pixel^2,  half way between the two nearest peaks

	// Find closest indexed peak to each of the fitted peaks, must also be withing tolerance
	Make/N=3/O/D hkl_MakeCalib_, qvec_MakeCalib_
	Wave hkl=hkl_MakeCalib_, qvec=qvec_MakeCalib_
	Variable px,py, dist, kBest, distBest, lambda
	for (m=0;m<Npeaks;m+=1)								// for each fitted peak, find corresponding indexed peak (if one is close enough)
		px = FullPeakList[m][0]								// fitted peak positions
		py = FullPeakList[m][1]
		for (k=0,distBest=Inf; k<Nindex; k+=1)
			dist = (px-FullPeakIndexed[k][9])^2 + (py-FullPeakIndexed[k][10])^2
			if (dist < distBest)
				distBest = dist
				kBest = k
			endif
		endfor
		if (distBest > tolerance)								// this peak not indexed, it is of no use to the calibration list
			continue
		endif

		hkl = FullPeakIndexed[kBest][p+3]					// hkl for this fitted peak
		cList[m][0] = px										// CalibrationList[][0,1] = (px,py) are measured on image
		cList[m][1] = py
		cList[m][2,4] = hkl[q-2]							// CalibrationList[][2-4] = hkl are hkl for each measured spot

		MatrixOp/O qvec =  recip0 x hkl
		qvec = abs(qvec)<1e-14 ? 0 : qvec
		cList[m][5,7] = qvec[q-5]							// CalibrationList[][5-7] = Q-vector, calculated from knowing hkl and the standard reciprocal lattice (NOT rotated from standard orientation)

		if (WaveExists(Emeasured))
			cList[m][8] = Emeasured[m]					// CalibrationList[][8] = keV, measured energy of spot
			lambda = hc/Emeasured[m]
			cList[m][9] = asin(lambda*norm(qvec)/4/PI)*180/PI // CalibrationList[][9] = theta (deg), obtained from calculated d(hkl) and measured energy using: Q = 4¹ sin(theta)/lambda
			deltaE[m] = (Emeasured[m] - FullPeakIndexed[kBest][7])*1e3	// not really needed, just to look at.
		endif
		cList[m][10] = NaN									// initially there is no calculated energy
		cList[m][11] = NaN
	endfor

	for (m=Npeaks-1;m>=0;m-=1)							// remove all bad lines
		if (numtype(cList[m][0]+cList[m][1]+cList[m][2]+cList[m][3]+cList[m][4]+cList[m][5]+cList[m][6]+cList[m][7]))
			DeletePoints/M=0 m, 1, cList
		endif
	endfor
	KillWaves/Z recip0_MakeCalib_,qvec_MakeCalib_,hkl_MakeCalib_

	str = GetWavesDataFolder(cList,2)
	printf "results are in:  '%s'\r",str
	return str
End
//
Static Function/T findSampleAxis(windex,[perfect])		// find the axis of the sample from FullPeakIndexed waves
	Wave windex
	Variable perfect											// if set to true, then always return the perfect Si rotation
	perfect = ParamIsDefault(perfect) ? 0 : perfect		// default to NOT perfect

	String wName=""
	if (!WaveInClass(windex,"IndexedPeakList") && !perfect)
		String wlist = reverseList(WaveListClass("IndexedPeakList","*",""))+"perfect Si;"
		Prompt wName, "Indexed List to provide sample rotation",popup,wlist
		DoPrompt "Sample Rotation",wName
		if (V_flag)
			return ""
		endif
		Wave windex = $wName
		perfect = stringmatch(wName,"perfect Si")
	endif

	if (perfect)
		Make/N=3/O/D perfectSiAxis = {-2.2256407716336,-0.884037203122317,-0.38914203759791}
		Note/K perfectSiAxis,"waveClass=sampleAxis;"
		return GetWavesDataFolder(perfectSiAxis,2)
	endif

	// get rotation matrix from wave note of windex
	Variable r00,r10,r20, r01,r11,r21, r02,r12,r22
	sscanf StringByKey("rotation_matrix0",note(windex),"="), "{{%g,%g,%g}{%g,%g,%g}{%g,%g,%g}}", r00,r10,r20, r01,r11,r21, r02,r12,r22
	if (V_flag!=9)
		printf "Unable to get sample axis from '"+NameOfWave(windex)+"', cannot get rotation_matrix0 from wave note"
		return ""
	endif

	wName = NameOfWave(windex)+"_SampleAxis"		// name of wave to hold result
	if (strlen(wName)>30)									// make a shorter wave name
		wName = ReplaceString("FullPeakIndexed_",NameOfWave(windex),"")
		wName = ReplaceString("FullPeakIndexed",wName,"")
		wName = wName[0,19]+"_SampleAxis"
	endif

	Make/N=3/O/D $wName
	Wave axis = $wName
	Note/K axis, ReplaceStringByKey("waveClass",note(windex),"sampleAxis","=")

	String rname = UniqueName("rotMat",1,0)				// calculate the axis
	Make/N=(3,3)/O/D $rname=NaN
	Wave rotMat = $rname
	rotMat[0][0]=r00;	rotMat[0][1]=r01;	rotMat[0][2]=r02
	rotMat[1][0]=r10;	rotMat[1][1]=r11;	rotMat[1][2]=r12
	rotMat[2][0]=r20;	rotMat[2][1]=r21;	rotMat[2][2]=r22
	Variable angle = axisOfMatrix(rotMat,axis)				// get axis from rotation matrix
	KillWaves/Z rotMat
	axis *= angle*PI/180									// re-scale axis to length in radians
	return GetWavesDataFolder(axis,2)
End



Function/WAVE EnterMeasuredEnergies(peakID)				// Allows easy entry of Emeasured[]
	String peakID				// method used to id measured peaks can only be "Pixel Position" or "hkl"

	String peakMethods="hkl;Pixel Position"		// list of valid methods
	String peakLists = WaveListClass("FittedPeakList","*","DIMS:2")
	String indexLists = WaveListClass("IndexedPeakList","*","MINCOLS:2")
	if (ItemsInList(peakLists)<1)
		DoAlert 0,"No FullPeakList files found. You need to fit an image first"
	endif
	if (ItemsInLIst(indexLists)<1)					// only possibility when no indexed peaks wave
		peakID = "Pixel Position"
	endif

	Variable i=WhichListItem(peakID,peakMethods)
	if (i<0)
		peakID = SelectString(i<0,peakID,"hkl")
		Prompt peakID "Method for identifying Mesaured Peaks",popup,peakMethods
		DoPrompt "Measured Peaks",peakID
		if (V_flag)
			return $""
		endif
	endif
	Variable Nmeas=-1, ishkl=stringmatch(peakID,"hkl")
	if (WhichListItem(peakID,peakMethods)<0)
		return $""
	endif

	if (ishkl)										// will enter energies using hkl's
		if (ItemsInList(indexLists)<1)
			return $""
		elseif (ItemsInList(indexLists)==1)
			Wave FullPeakIndexed = $StringFromList(0,indexLists)
		else
			String indexName
			Prompt indexName,"Indexing to use",popup,indexLists
			DoPrompt "Indexing",indexName
			if (V_flag)
				return $""
			endif
			Wave FullPeakIndexed = $indexName
		endif
		Wave FullPeakList=$StringByKey("peakListWave",note(FullPeakIndexed),"=")
		Nmeas = DimSize(FullPeakIndexed,0)
	else												// will enter energies using peak positions
		if (ItemsInList(peakLists)<1)
			return $""
		elseif (ItemsInList(peakLists)==1)
			Wave FullPeakList = $StringFromList(0,peakLists)
		else
			String peakListName
			Prompt peakListName,"Peaks list to use",popup,peakLists
			DoPrompt "Peaks",peakListName
			if (V_flag)
				return $""
			endif
			Wave FullPeakList = $peakListName
		endif
		Nmeas = DimSize(FullPeakList,0)
	endif
	if (!WaveExists(FullPeakList))
		DoAlert 0,"Could not find FullPeakList"
		return $""
	endif
	Variable Nfit=DimSize(FullPeakList,0)			// number of fitted peaks
	String detectorID=StringByKey("detectorID",note(FullPeakList),"=")
	String color=detectorID2color(detectorID)		// detector color, {orange, yellow, purple}

	String idName=UniqueName("PeakLabels",1,0), eName=UniqueName("Energies_keV",1,0)
	Make/N=(Nmeas) $eName/WAVE=ew=NaN
	if (ishkl)										// entering energies for hkl
		Make/N=(Nmeas,3) $idName/WAVE=idw
		idw = FullPeakIndexed[p][3+q]
		SetDimLabel 1,0,H,idw
		SetDimLabel 1,1,K,idw
		SetDimLabel 1,2,L,idw
		Edit/K=1/W=(214,245,571,738)/N=EnergyInput idw.ld,ew.y
		ModifyTable format(Point)=1,width(Point)=40,format(ew.y)=3,digits(ew.y)=4
		ModifyTable width(ew.y)=90,width(idw.l)=20,width(idw.d)=50
	else												// entering energies for pixel positions
		Make/N=(Nmeas,2) $idName/WAVE=idw
		idw = FullPeakList[p][q]
		SetDimLabel 1,0,X_pixel,idw
		SetDimLabel 1,1,Y_pixel,idw
		Edit/K=1/W=(214,245,534,738)/N=EnergyInput idw.ld,ew.y
		ModifyTable format(Point)=1,width(Point)=40,format(ew.y)=3,digits(ew.y)=4
		ModifyTable width(ew.y)=90,width(idw.l)=20,width(idw.d)=74
	endif
	DoUpdate
	DoAlert 0,"Kill the table after you have filled it in.\rYou do not need to fill in energies for each line."
	PauseForUser EnergyInput,EnergyInput
	ew = numtype(ew) ? NaN : ew					// ensure that there are no Inf's

	eName=UniqueName("Emeasured"+color,1,0)
	Make/N=(Nfit) $eName/WAVE=Emeasured=NaN
	String wNote="waveClass=measuredEnergies;"
	wNote = ReplaceStringByKey("detectorID",wNote,detectorID,"=")
	Note/K Emeasured, wNote

	Variable px,py, m
	Make/N=(Nfit,2)/FREE pxy=FullPeakList[p][q]
	Make/N=2/FREE pxy0
	if (ishkl)
		Variable Nindex=DimSize(FullPeakIndexed,0)
		Make/N=3/FREE idwm
		Make/N=(Nindex,3)/FREE hkl
		hkl = FullPeakIndexed[p][q+3]
		Make/N=(Nindex,3)/FREE hkl=FullPeakIndexed[p][q+3]	// list of hkl's that are in FullPeakIndexed
		if (DimSize(idw,0)<Nindex)
			i = DimSize(idw,0)
			Redimension/N=(Nindex,-1) idw
			idw[i,Nindex-1][] = NaN					// hkls used to input energies, now same length as hkl[][3]
		endif
		for (m=0;m<DimSize(ew,0);m+=1)			// for each of the measured energies
			if (!(ew[m]>0))
				continue								// skip invalid energies
			endif
			idwm = idw[m][p]							// hkl of one measurement
			MatrixOP/FREE dhkl = sumRows(magSqr(hkl-rowRepeat(idwm,Nindex)))
			WaveStats/M=1/Q dhkl						// find row in FullPeakIndexed with the hkl from idwm
			if (!(V_min<0.01 && V_minloc>=0))
				continue								// skip if I could not find the hkl in hkl[][3]
			endif
			pxy0 = {FullPeakIndexed[V_minloc][9],FullPeakIndexed[V_minloc][10]}	// pixel coordinates for this E measurement
			MatrixOP/FREE dpxy = sumRows(magSqr(pxy-rowRepeat(pxy0,Nfit)))
			WaveStats/M=1/Q dpxy						// find row in pxy which is closest to pxy0, gives row index in FullPeakList
			Emeasured[V_minloc] = ew[m]				// save into Emeasured
		endfor
	else
		for (m=0;m<DimSize(ew,0);m+=1)			// for each measured energy
			if (!(ew[m]>0))
				continue								// skip invalid energies
			endif
			pxy0 = idw[m][p]							// pixel where energy was measured
			MatrixOP/FREE dpxy = sumRows(magSqr(pxy-rowRepeat(pxy0,Nfit)))
			WaveStats/M=1/Q dpxy						// find row in pxy which is closest to pxy0
			Emeasured[V_minloc] = ew[m]				// save into Emeasured
		endfor
	endif
	KillWaves/Z idw,ew
	if (strlen(GetRTStackInfo(2))<1 || stringmatch(GetRTStackInfo(2),"CalibrationButtonProc"))
		printf "Created  '%s'  containing the measured energies\r",NameOfWave(Emeasured)
	endif
	return Emeasured
End



// this uses pre-defined CalibrationList's and takes them as input, it sets up and optimizes each of the three detectors
Function OptimizeAll(calib0,calib1,calib2)
	Wave calib0,calib1,calib2								// calibration lists
	if (WaveExists(calib0))
		if (!WaveInClass(calib0,"DetectorCalibrationList0"))
			DoAlert 0, "The first argument "+NameOfWave(calib0)+" is not of 'DetectorCalibrationList0' class"
			return 1
		endif
	endif
	if (WaveExists(calib1))
		if (!WaveInClass(calib1,"DetectorCalibrationList*"))
			DoAlert 0, "The second argument "+NameOfWave(calib1)+" is not of 'DetectorCalibrationList' class"
			return 1
		endif
	endif
	if (WaveExists(calib2))
		if (!WaveInClass(calib2,"DetectorCalibrationList*"))
			DoAlert 0, "The third argument "+NameOfWave(calib2)+" is not of 'DetectorCalibrationList' class"
			return 1
		endif
	endif

	Variable printing = (!strlen(GetRTStackInfo(2))) || stringmatch(GetRTStackInfo(2),"testOneOptimize")
	String cName0,cName1,cName2
	if (!WaveExists(calib1) || !WaveExists(calib1) || !WaveExists(calib1))
		cName0 = SelectString(WaveExists(calib0),"CalibrationList0",NameOfWave(calib0))
		cName1 = SelectString(WaveExists(calib1),"CalibrationList1",NameOfWave(calib1))
		cName2 = SelectString(WaveExists(calib2),"CalibrationList2",NameOfWave(calib2))
//		String wList0 = WaveListClass("DetectorCalibrationList0","*","DIMS:2;MAXCOLS:10;MINCOLS:10")+"_none_"
//		String wList12 = WaveListClass("DetectorCalibrationList","*","DIMS:2;MAXCOLS:10;MINCOLS:10")+"_none_"
		String wList0 = WaveListClass("DetectorCalibrationList0","*","DIMS:2;MAXCOLS:12;MINCOLS:10")+"_none_"
		String wList12 = WaveListClass("DetectorCalibrationList","*","DIMS:2;MAXCOLS:12;MINCOLS:10")+"_none_"
		Prompt cName0,"top detector",popup,wList0
		Prompt cName1,"side detector 1",popup,wList12
		Prompt cName2,"side detector 2",popup,wList12
		DoPrompt "Calibration Lists",cName0,cName1,cName2
		if (V_flag)
			return 1
		endif
		cName0 = SelectString(exists(cName0)==1,"",cName0)
		cName1 = SelectString(exists(cName1)==1,"",cName1)
		cName2 = SelectString(exists(cName2)==1,"",cName2)
		Wave calib0=$cName0, calib1=$cName1, calib2=$cName2

		printing = 1
		cName0 = SelectString(strlen(cName0),"$\"\"",cName0)
		cName1 = SelectString(strlen(cName1),"$\"\"",cName1)
		cName2 = SelectString(strlen(cName2),"$\"\"",cName2)
 		printf "OptimizeAll(%s,%s,%s)\r",cName0,cName1,cName2
	endif

	String calib0FullName="", calib1FullName="", calib2FullName=""
	if (WaveExists(calib0))
		calib0FullName = GetWavesDataFolder(calib0,2)
	endif
	if (WaveExists(calib1))
		calib1FullName = GetWavesDataFolder(calib1,2)
	endif
	if (WaveExists(calib2))
		calib2FullName = GetWavesDataFolder(calib2,2)
	endif

	if (WaveExists(calib0))
		if (!WaveInClass(calib0,"DetectorCalibrationList0"))
			DoAlert 0, "The first argument "+NameOfWave(calib0)+" is not of 'DetectorCalibrationList0' class"
			return 1
		endif
	endif
	if (WaveExists(calib1))
		if (!WaveInClass(calib1,"DetectorCalibrationList*"))
			DoAlert 0, "The second argument "+NameOfWave(calib1)+" is not of 'DetectorCalibrationList' class"
			return 1
		endif
		if (	stringmatch(calib0FullName,calib1FullName))
			DoAlert 0, "Calibration Lists 0 & 1 are the same wave, '"+NameOfWave(calib0)+"'"
			return 1
		endif
	endif
	if (WaveExists(calib2))
		if (!WaveInClass(calib2,"DetectorCalibrationList*"))
			DoAlert 0, "The third argument "+NameOfWave(calib2)+" is not of 'DetectorCalibrationList' class"
			return 1
		endif
		if (	stringmatch(calib1FullName,calib2FullName))
			DoAlert 0, "Calibration Lists 1 & 2 are the same wave, '"+NameOfWave(calib1)+"'"
			return 1
		endif
	endif
	String cNameList = calib0FullName+";"+calib1FullName+";"+calib2FullName+";"

	STRUCT microGeometry g
	FillGeometryStructDefault(g)							//fill the geometry structure with current values
	//	print "in testOptimizeAll:"
	//	printGeometry(g)

	Variable Rstart, Rend, Rangle
	Variable err, sec, failed=0
	String noteStr, errList=""
	Variable rhox=NaN, rhoy=NaN, rhoz=NaN
	STRUCT detectorGeometry d								// this is the structure that will be optimized
	Variable err3=0											// sum of error from all 3 detectors
	Variable dNum												// detector number {0,1,2}
	for (dNum=0;dNum<MAX_Ndetectors && !failed;dNum+=1)
		if (!(g.d[dNum].used))
			continue											// skip un-used detectors
		endif
		Wave cListN = $StringFromList(dNum,cNameList)	// data to use for optimization
		if (!WaveExists(cListN))
			continue
		endif

		noteStr=note(cListN)
		if (dNum && numtype(rhox+rhoy+rhoz)==0)
			noteStr = ReplaceNumberByKey("rhox",noteStr,rhox,"=")
			noteStr = ReplaceNumberByKey("rhoy",noteStr,rhoy,"=")
			noteStr = ReplaceNumberByKey("rhoz",noteStr,rhoz,"=")
			Note/K cListN, noteStr
		endif
		CopyDetectorGeometry(d,g.d[dNum])				// set new detector structure to start as the old one
		Rstart = sqrt((d.R[0]*d.R[0])+(d.R[1]*d.R[1])+(d.R[2]*d.R[2]))*180/PI
		if (printing)
			printf "Detector %d\r",dNum
			Rangle = sqrt((g.d[dNum].R[0])^2 + (g.d[dNum].R[1])^2 + (g.d[dNum].R[2])^2)*180/PI
			printf "started at  R={%g,%g,%g},   P={%g,%g,%g}mm,   |R| = %g¡\r",d.R[0],d.R[1],d.R[2],(d.P[0])/1000,(d.P[1])/1000,(d.P[2])/1000,Rstart
			Variable angle = sqrt(NumberByKey("rhox",noteStr,"=")^2 + NumberByKey("rhoy",noteStr,"=")^2 + NumberByKey("rhoz",noteStr,"=")^2)*180/PI
			printf "sample started at axis = {%g, %g, %g},   |sample angle| = %g¡\r",NumberByKey("rhox",noteStr,"="),NumberByKey("rhoy",noteStr,"="),NumberByKey("rhoz",noteStr,"="),angle
		endif

		failed = OptimizeDetectorCalibration(d,cListN)
		if (failed && printing)
			printf "-----------------------\rERROR in Optimize\r  %s\r  -----------------------\r",note(cListN)
		endif
		noteStr=note(cListN)
		if ( NumberByKey("V_OptNumIters",noteStr,"=")<50 && dNum==0)
			printf "-----------------------\rPOSSIBLE ERROR (no. of iterations = %d is too few), re-run this\r  -----------------------\r", NumberByKey("V_OptNumIters",noteStr,"=")
		endif
		Rend = sqrt((d.R[0]*d.R[0])+(d.R[1]*d.R[1])+(d.R[2]*d.R[2]))*180/PI
		err = NumberByKey("err",noteStr,"=")
		sec = NumberByKey("exectutionSec",noteStr,"=")
		err3 += err
		errList += num2str(err)+";"
		if (dNum==0)
			rhox = NumberByKey("rhox",noteStr,"=")
			rhoy = NumberByKey("rhoy",noteStr,"=")
			rhoz = NumberByKey("rhoz",noteStr,"=")
		endif
		if (printing)
			if (sec>4)
				printf "\tOptimization toook %s\r",Secs2Time(sec,5,1)
			endif
			printf "ended at  R={%g,%g,%g},   P={%g,%g,%g}mm,   |R| = %g¡\r",d.R[0],d.R[1],d.R[2],(d.P[0])/1000,(d.P[1])/1000,(d.P[2])/1000,Rend
			if (dNum==0)
				printf "  final sample rotation is  {%g,%g,%g},   |sample angle| = %g¡\r",rhox,rhoy,rhoz,sqrt(rhox^2 + rhoy^2 + rhoz^2)*180/PI
			endif
			printf "error started at %g,   reduced to  %g,   after %d iterations\r",NumberByKey("errStart",noteStr,"="),err, NumberByKey("V_OptNumIters",noteStr,"=")
			printf "(final - initial) --> ÆR={%.2g,%.2g,%.2g},   ÆP={%.3f,%.3f,%.3f}\r",(d.R[0]-g.d[dNum].R[0]), (d.R[1]-g.d[dNum].R[1]), (d.R[2]-g.d[dNum].R[2]), (d.P[0]-g.d[dNum].P[0])/1000, (d.P[1]-g.d[dNum].P[1])/1000, (d.P[2]-g.d[dNum].P[2])/1000
		endif
		CopyDetectorGeometry(g.d[dNum],d)				// updagte structure with fitted values for this detector
	endfor
	err3 /= (g.Ndetectors)
	if (printing && numtype(Rend)==0)
		print " "
		printGeometry(g)
		print "errList =",errList,"     total error =",err3
	endif

	if (!failed && printing)
		DoAlert 1,"Update current Geometry with these fitted values, err = "+num2str(err3)
		if (V_flag==1)
			print "Updated the current geometry with these values"
			UpdateDefaultGeometryStruct(g,local=1)
		else
			print "REJECTED these values, no update done"
		endif
	endif
	return (failed ? NaN : err3)
End


Static Function OptimizeDetectorCalibration(d,CalibrationList)	// optimizes values in d to best fit values in CalibrationList
	STRUCT detectorGeometry &d						// detector geometry to optimize
	Wave CalibrationList								// data to use for optimization
	if (!WaveInClass(CalibrationList,"DetectorCalibrationList*"))
		DoAlert 0, "The CalibrationList wave is not of waveClass, 'DetectorCalibrationList*'"
		return -1
	endif

	Variable timer=startMSTimer
	Make/N=3/O/D ki={0,0,1}, kf, qcalc,qmeas
	Make/N=(3,3)/O/D rhoCalib

	String str, noteStr=note(CalibrationList)
	Variable zero = WaveInClass(CalibrationList,"DetectorCalibrationList0")	// a zero detector, fit rho too
	Variable rhox=NumberByKey("rhox",noteStr,"="), rhoy=NumberByKey("rhoy",noteStr,"="), rhoz=NumberByKey("rhoz",noteStr,"=")

	d.R[0] = (d.R[0] == 0) ? 0.02 : d.R[0]				// avoid starting exactly on a zero
	d.R[1] = (d.R[1] == 0) ? 0.02 : d.R[1]				// offset 0 angles by ~1¡
	d.R[2] = (d.R[2] == 0) ? 0.02 : d.R[2]
	d.P[0] = (d.P[0] == 0) ? 500 : d.P[0]				// offset 0 positions by 0.5mm
	d.P[1] = (d.P[1] == 0) ? 500 : d.P[1]
	d.P[2] = (d.P[2] == 0) ? 500 : d.P[2]

//	Note/K CalibrationList,AddListItem("printErrX",noteStr)
	Variable err, errStart
	if (zero)
		errStart=CalibrationErrorRPrho(CalibrationList,d.R[0],d.R[1],d.R[2],d.P[0],d.P[1],d.P[2],rhox,rhoy,rhoz)	// save error at start
	else
		errStart=CalibrationErrorRP(CalibrationList,d.R[0],d.R[1],d.R[2],d.P[0],d.P[1],d.P[2])							// save error at start
	endif
//	Note/K CalibrationList,noteStr

	Variable maxIters=NumVarOrDefault("maxIters",1000)
	if (zero)
		Make/N=9/O/D optimizeStepWave={d.R[0],d.R[1],d.R[2],d.P[0],d.P[1],d.P[2],rhox,rhoy,rhoz}
		optimizeStepWave = abs(optimizeStepWave)
		optimizeStepWave = max(optimizeStepWave[p],0.1)
		Optimize/M={0,0}/X={d.R[0],d.R[1],d.R[2],d.P[0],d.P[1],d.P[2],rhox,rhoy,rhoz}/Y=(errStart)/Q/I=(maxIters)/R=optimizeStepWave CalibrationErrorRPrho, CalibrationList
	else
		Make/N=6/O/D optimizeStepWave={d.R[0],d.R[1],d.R[2],d.P[0],d.P[1],d.P[2]}
		optimizeStepWave = abs(optimizeStepWave)
		optimizeStepWave = max(optimizeStepWave[p],0.1)
		Optimize/M={0,0}/X={d.R[0],d.R[1],d.R[2],d.P[0],d.P[1],d.P[2]}/Y=(errStart)/Q/I=(maxIters)/R=optimizeStepWave CalibrationErrorRP, CalibrationList
	endif
	noteStr = note(CalibrationList)
	Wave W_Extremum=W_Extremum

	if (V_flag)
		noteStr = ReplaceNumberByKey("ERROR",noteStr,V_flag,"=")
	endif
	noteStr = ReplaceNumberByKey("errStart",noteStr,errStart,"=")
	noteStr = ReplaceNumberByKey("err",noteStr,V_min,"=")
	noteStr = ReplaceNumberByKey("V_OptTermCode",noteStr,V_OptTermCode,"=")
	noteStr = ReplaceNumberByKey("V_OptNumIters",noteStr,V_OptNumIters,"=")
	noteStr = ReplaceNumberByKey("V_OptNumFunctionCalls",noteStr,V_OptNumFunctionCalls,"=")
	noteStr = ReplaceNumberByKey("V_min",noteStr,V_min,"=")
	if (zero)												// only update rho if I fit it
		noteStr = ReplaceNumberByKey("rhox",noteStr,W_Extremum[6],"=")
		noteStr = ReplaceNumberByKey("rhoy",noteStr,W_Extremum[7],"=")
		noteStr = ReplaceNumberByKey("rhoz",noteStr,W_Extremum[8],"=")
	endif
	noteStr = ReplaceNumberByKey("exectutionSec",noteStr,stopMSTimer(timer)*1e-6,"=")
	Note/K CalibrationList, noteStr

	Variable dR = 1e-8									// this gives a limit of (~0.5e-6)¡
	d.R[0]= round(W_Extremum[0]/dR)*dR			// round to 8 places
	d.R[1]= round(W_Extremum[1]/dR)*dR
	d.R[2]= round(W_Extremum[2]/dR)*dR
	d.R[0] = (d.R[0]==0) ? 0 : d.R[0]					// change -0 to +0
	d.R[1] = (d.R[1]==0) ? 0 : d.R[1]
	d.R[2] = (d.R[2]==0) ? 0 : d.R[2]

	Variable dP = min(d.sizeX/d.Nx,d.sizeY/d.Ny)/100	// resolution limit in Px,Py,Pz (1% of a pixel)
	dP = 10^round(log(dP))							// make rounding to a fixed number of places (not arbitrary number)
	d.P[0]= round(W_Extremum[3]/dP)*dP
	d.P[1]= round(W_Extremum[4]/dP)*dP
	d.P[2]= round(W_Extremum[5]/dP)*dP
	d.P[0] = (d.P[0]==0) ? 0 : d.P[0]					// change -0 to +0
	d.P[1] = (d.P[1]==0) ? 0 : d.P[1]
	d.P[2] = (d.P[2]==0) ? 0 : d.P[2]

	d.geoNote = "Optimized using "+NameOfWave(CalibrationList)
	sprintf str,"%s, %s (%g)",date(),Secs2Time(DateTime,3),date2secs(-1,-1,-1)/3600
	d.timeMeasured = str

	CalibrationList[][10] = CalculatedQvecEnergy(CalibrationList,p,qcalc)	// fill calculated energy from sample orientation and known lattice constants
	CalibrationList[][11] = ( CalibrationList[p][8] - CalibrationList[p][10] ) *1e3	// ÆE (eV)

	KillWaves/Z W_Extremum, W_OptGradient, optimizeStepWave
	KillWaves/Z ki, kf, qcalc,qmeas, rhoCalib
	return V_flag
End
//
// This routine is for optimizing the detector posiition & rotation and sample orientation simultaneously (only useful for detector 1,2,...   not 0)
Function CalibrationErrorRP(CalibrationList,Rx,Ry,Rz,Px,Py,Pz)	// returns error between measured and calculated spots, only called by OptimizeDetectorCalibration()
	Wave CalibrationList
	Variable Rx,Ry,Rz						// rotation vector for detector
	Variable Px,Py,Pz						// translation vector for detector
	// This routine uses the (px,py) from CalibrationList, and compares it to the angles calculated from (hkl) and lattice, it also
	// compares the measured theta (from px,py) to the theta calculated from the measured energy
	// The routine does not use the {qx,qy,qz}, or the theta columns in CalibrationList
	//
	//	CalibrationList[][0,1]	= (px,py) are measured on image
	//	CalibrationList[][2-4]	= hkl are hkl for each measured spot
	//	CalibrationList[][5-7]	= Q-vector, calculated from knowing hkl and the standard reciprocal lattice (NOT rotated from standard orientation)
	//	CalibrationList[][8]	= keV, measured energy of spot
	//	CalibrationList[][9]	= theta (deg), obtained from calculated d(hkl) and measured energy using: lambda = 2d sin(theta)
	//	CalibrationList[][10]	= keV calculated from sample orientation and known lattice (does not use measured spot position)
	//	CalibrationList[][11]	= ÆE (eV),  ( CalibrationList[][8] - CalibrationList[][10] ) *1e3

	//	printf "CalibrationErrorRP(%s,%g,%g,%g,   %g,%g,%g)\r",NameOfWave(CalibrationList),Rx,Ry,Rz,Px,Py,Pz
	Wave ki=ki, kf=kf									// preallocated waves (speeds things up)
	Wave rhoSample=rhoCalib							// rotation matrix computed from {rhox,rhoy,rhoz}
	Wave qcalc=qcalc, qmeas=qmeas
	String noteStr = note(CalibrationList)

	Variable rhox=NumberByKey("rhox",noteStr,"="), rhoy=NumberByKey("rhoy",noteStr,"="), rhoz=NumberByKey("rhoz",noteStr,"=")
	RotationVec2Matrix(rhox,rhoy,rhoz,NaN,rhoSample)// Compute the sample rotation rhoSample. Apply to qhat's to make them line up with (px,py) pairs

	STRUCT detectorGeometry d							// a local version of detector structure
	d.used = 1
	d.Nx = NumberByKey("Nx",noteStr,"=")			// number of un-binned pixels in whole detector
	d.Ny = NumberByKey("Ny",noteStr,"=")
	d.sizeX = NumberByKey("sizeX",noteStr,"=")		// outside size of detector (micron)
	d.sizeY = NumberByKey("sizeY",noteStr,"=")		// outside size of detector (micron)
	d.R[0] = Rx;	d.R[1] = Ry;	d.R[2] = Rz			// Rotation vector for detector (degree)
	d.P[0] = Px;	d.P[1] = Py;	d.P[2] = Pz			// offset to detector (micron)
	DetectorUpdateCalc(d)								// update all fields in this detector structure (basically rho)
	//print "in CalibrationErrorRP:"
	//	printDetector(d)

	Variable N=DimSize(CalibrationList,0)				// nuber of spots to use in computing error values
	ImageStats/M=1/G={0,N-1,0,1} CalibrationList	// weight is used to more heavily weight the points with valid energy
	Variable Ntot=V_npnts/2								// probable number of usable points
	ImageStats/M=1/G={0,N-1,8,8} CalibrationList
	Variable NeV=V_npnts								// probable number of valid points with energy
	Variable weightLo=Ntot/(2*Ntot-NeV), weightHi=weightLo*Ntot/NeV
	if (Nev<=0 || abs(2*Ntot-NeV)<0.1 || numtype(weightLo+weightHi) || NeV>=Ntot)
		weightLo = 1										// something wrong, no weighting
		weightHi = 1
	endif

	Variable keV, Qlen, dq, i
	Variable errQ=0, nQ=0
	for (i=0;i<N;i+=1)
		if (numtype(CalibrationList[i][0]+CalibrationList[i][1]+CalibrationList[i][5]+CalibrationList[i][6]+CalibrationList[i][7]))
			continue										// skip invalid points
		endif

		// compute ideal Q-vector from just rhoSample and the reciprocal lattice
		qcalc = CalibrationList[i][p+5]					// Q vector for this reflection in stardard orientation
		MatrixOp/O qcalc = rhoSample x qcalc			// rotate Q-vector by rhoSample

		// construct measured Q-vectors from (px,py) & keV
		pixel2XYZ(d,CalibrationList[i][0],CalibrationList[i][1],kf)		// find kf, convert pixel position to the beam line coordinate system
		normalize(kf)
		qmeas = kf - ki									// wrong length, but parallel to correct q
		normalize(qmeas)
		keV = CalibrationList[i][8]
		// compute |Q| using:  sin(theta) = -MatrixDot(ki,qmeas),  Q=4¹ sin(theta)/lambda,  if no energy, then use |Q| of calculated spot
		Qlen = (keV>0) ? -4*PI*MatrixDot(ki,qmeas)*keV/hc : norm(qcalc)
		qmeas *= Qlen

		qcalc -= qmeas
		dq = norm(qcalc)
		dq *= (keV>0) ? weightHi : weightLo				// boost when I have a measured energy
		errQ += dq^2										// error is distance in Q space
		nQ += 1
	endfor
	errQ /= nQ

	if (WhichListItem("printErrX",noteStr)>=0)
		//	printf "{%g,%g,%g,   %g,%g,%g} --> errQ/nQ = %g\r",Rx,Ry,Rz,Px,Py,Pz,errQ
		printf "errQ/nQ = %g,  nQ=%d\r",errQ,nQ
	endif
	return 	errQ										// error between calculated and measured
End
//Function CalibrationErrorRP(CalibrationList,Rx,Ry,Rz,Px,Py,Pz)	// returns error between measured and calculated spots, only called by OptimizeDetectorCalibration()
//	Wave CalibrationList
//	Variable Rx,Ry,Rz						// rotation vector for detector
//	Variable Px,Py,Pz						// translation vector for detector
//	// This routine uses the (px,py) from CalibrationList, and compares it to the angles calculated from (hkl) and lattice, it also
//	// compares the measured theta (from px,py) to the theta calculated from the measured energy
//	// The routine does not use the {qx,qy,qz}, or the theta columns in CalibrationList
//	//
//	//	CalibrationList[][0,1]	= (px,py) are measured on image
//	//	CalibrationList[][2-4]	= hkl are hkl for each measured spot
//	//	CalibrationList[][5-7]	= Q-vector, calculated from knowing hkl and the standard reciprocal lattice (NOT rotated from standard orientation)
//	//	CalibrationList[][8]	= keV, measured energy of spot
//	//	CalibrationList[][9]	= theta (deg), obtained from calculated d(hkl) and measured energy using: lambda = 2d sin(theta)
//	//	CalibrationList[][10]	= keV calculated from sample orientation and known lattice (does not use measured spot position)
//	//	CalibrationList[][11]	= ÆE (eV),  ( CalibrationList[][8] - CalibrationList[][10] ) *1e3
//
//	//	printf "CalibrationErrorRP(%s,%g,%g,%g,   %g,%g,%g)\r",NameOfWave(CalibrationList),Rx,Ry,Rz,Px,Py,Pz
//	Wave ki=ki, kf=kf									// preallocated waves (speeds things up)
//	Wave rhoSample=rhoCalib							// rotation matrix computed from {rhox,rhoy,rhoz}
//	Wave qcalc=qcalc, qmeas=qmeas
//	String noteStr = note(CalibrationList)
//
//	Variable rhox=NumberByKey("rhox",noteStr,"="), rhoy=NumberByKey("rhoy",noteStr,"="), rhoz=NumberByKey("rhoz",noteStr,"=")
//	RotationVec2Matrix(rhox,rhoy,rhoz,NaN,rhoSample)// Compute the sample rotation rhoSample. Apply to qhat's to make them line up with (px,py) pairs
//
//	STRUCT detectorGeometry d							// a local version of detector structure
//	d.used = 1
//	d.Nx = NumberByKey("Nx",noteStr,"=")			// number of un-binned pixels in whole detector
//	d.Ny = NumberByKey("Ny",noteStr,"=")
//	d.sizeX = NumberByKey("sizeX",noteStr,"=")		// outside size of detector (micron)
//	d.sizeY = NumberByKey("sizeY",noteStr,"=")		// outside size of detector (micron)
//	d.R[0] = Rx;	d.R[1] = Ry;	d.R[2] = Rz			// Rotation vector for detector (degree)
//	d.P[0] = Px;	d.P[1] = Py;	d.P[2] = Pz			// offset to detector (micron)
//	DetectorUpdateCalc(d)								// update all fields in this detector structure (basically rho)
//	//print "in CalibrationErrorRP:"
//	//	printDetector(d)
//
//	Variable keV, Qlen, dq
//	Variable i, N=DimSize(CalibrationList,0)			// nuber of spots to use in computing error values
//	Variable errQ=0, nQ=0
//	for (i=0;i<N;i+=1)
//		if (numtype(CalibrationList[i][0]+CalibrationList[i][1]+CalibrationList[i][p+5]+CalibrationList[i][p+6]+CalibrationList[i][p+7]))
//			continue										// skip invalid points
//		endif
//
//		// compute ideal Q-vector from just rhoSample and the reciprocal lattice
//		qcalc = CalibrationList[i][p+5]					// Q vector for this reflection in stardard orientation
//		MatrixOp/O qcalc = rhoSample x qcalc			// rotate Q-vector by rhoSample
//
//		// construct measured Q-vectors from (px,py) & keV
//		pixel2XYZ(d,CalibrationList[i][0],CalibrationList[i][1],kf)		// find kf, convert pixel position to the beam line coordinate system
//		normalize(kf)
//		qmeas = kf - ki									// wrong length, but parallel to correct q
//		normalize(qmeas)
//		keV = CalibrationList[i][8]
//		// compute |Q| using:  sin(theta) = -MatrixDot(ki,qmeas),  Q=4¹ sin(theta)/lambda,  if no energy, then use |Q| of calculated spot
//		Qlen = keV>0 ? -4*PI*MatrixDot(ki,qmeas)*keV/hc : norm(qcalc)
//		qmeas *= Qlen
//
//		qcalc -= qmeas
//		dq = norm(qcalc)
//		errQ += dq^2										// error is distance in Q space
//		nQ += 1
//	endfor
//	errQ /= nQ
//
//	if (WhichListItem("printErrX",noteStr)>=0)
//		//	printf "{%g,%g,%g,   %g,%g,%g} --> errQ/nQ = %g\r",Rx,Ry,Rz,Px,Py,Pz,errQ
//		printf "errQ/nQ = %g,  nQ=%d\r",errQ,nQ
//	endif
//	return 	errQ										// error between calculated and measured
//End
//
//
// This routine is for optimizing the detector posiition & rotation and sample orientation simultaneously (only useful for detector0)
Function CalibrationErrorRPrho(CalibrationList,Rx,Ry,Rz,Px,Py,Pz,rhox,rhoy,rhoz)	// returns error between measured and calculated spots, only called by OptimizeDetectorCalibration()
	Wave CalibrationList
	Variable Rx,Ry,Rz						// rotation vector for detector
	Variable Px,Py,Pz						// translation vector for detector
	Variable rhox,rhoy,rhoz					// rotation vector to apply to qhat's to make them line up with (px,py) pairs
	// This routine uses the (px,py) from CalibrationList, and compares it to the angles calculated from (hkl) and lattice, it also
	// compares the measured theta (from px,py) to the theta calculated from the measured energy
	// The routine does not use the {qx,qy,qz}, or the theta columns in CalibrationList
	//
	//	CalibrationList[][0,1]	= (px,py) are measured on image
	//	CalibrationList[][2-4]	= hkl are hkl for each measured spot
	//	CalibrationList[][5-7]	= Q-vector, calculated from knowing hkl and the standard reciprocal lattice (NOT rotated from standard orientation)
	//	CalibrationList[][8]	= keV, measured energy of spot
	//	CalibrationList[][9]	= theta (deg), obtained from calculated d(hkl) and measured energy using: lambda = 2d sin(theta)
	//	CalibrationList[][10]	= keV calculated from sample orientation and known lattice (does not use measured spot position)
	//	CalibrationList[][11]	= ÆE (eV),  ( CalibrationList[][8] - CalibrationList[][10] ) *1e3

	//	printf "CalibrationErrorRPrho(%s,%g,%g,%g,   %g,%g,%g)\r",NameOfWave(CalibrationList),Rx,Ry,Rz,Px,Py,Pz
	Wave ki=ki, kf=kf									// preallocated waves (speeds things up)
	Wave rhoSample=rhoCalib							// rotation matrix computed from {rhox,rhoy,rhoz}
	Wave qcalc=qcalc, qmeas=qmeas

	RotationVec2Matrix(rhox,rhoy,rhoz,NaN,rhoSample)// compute rhoSample

	STRUCT detectorGeometry d							// a local version of detector structure
	String noteStr = note(CalibrationList)
	d.used = 1
	d.Nx = NumberByKey("Nx",noteStr,"=")			// number of un-binned pixels in whole detector
	d.Ny = NumberByKey("Ny",noteStr,"=")
	d.sizeX = NumberByKey("sizeX",noteStr,"=")		// outside size of detector (micron)
	d.sizeY = NumberByKey("sizeY",noteStr,"=")		// outside size of detector (micron)
	d.R[0] = Rx;	d.R[1] = Ry;	d.R[2] = Rz			// Rotation vector for detector (degree)
	d.P[0] = Px;	d.P[1] = Py;	d.P[2] = Pz			// offset to detector (micron)
	DetectorUpdateCalc(d)								// update all fields in this detector structure (basically rho)
	//print "in CalibrationErrorRPrho:"
	//	printDetector(d)

	Variable N=DimSize(CalibrationList,0)				// nuber of spots to use in computing error values
	ImageStats/M=1/G={0,N-1,0,1} CalibrationList	// weight is used to more heavily weight the points with valid energy
	Variable Ntot=V_npnts/2							// probable number of usable points
	ImageStats/M=1/G={0,N-1,8,8} CalibrationList
	Variable NeV=V_npnts								// probable number of valid points with energy
	Variable weightLo=Ntot/(2*Ntot-NeV), weightHi=weightLo*Ntot/NeV
	if (Nev<=0 || abs(2*Ntot-NeV)<0.1 || numtype(weightLo+weightHi) || NeV>=Ntot)
		weightLo = 1									// something wrong, no weighting
		weightHi = 1
	endif

	Variable keV, Q, dq, i
	Variable errQ=0, nQ=0
//	Variable dqMax=0
	for (i=0;i<N;i+=1)
		if (numtype(CalibrationList[i][0]+CalibrationList[i][1]+CalibrationList[i][5]+CalibrationList[i][6]+CalibrationList[i][7]))
			continue										// skip invalid points
		endif

		// compute ideal Q-vector from just rhoSample and the reciprocal lattice
		qcalc = CalibrationList[i][p+5]					// Q vector for this reflection in stardard orientation
		MatrixOp/O qcalc = rhoSample x qcalc			// rotate Q- vector by rhoSample

		// construct measured Q-vectors from (px,py) & keV
		pixel2XYZ(d,CalibrationList[i][0],CalibrationList[i][1],kf)		// find kf, convert pixel position to the beam line coordinate system
		normalize(kf)
		qmeas = kf - ki									// wrong length, but parallel to correct q
		normalize(qmeas)
		keV = CalibrationList[i][8]
		// compute |Q| using  sin(theta) = -MatrixDot(ki,qmeas),  Q=4¹ sin(theta)/lambda,  if no energy, then use |Q| of calculated spot
		Q = (keV>0) ? -4*PI*MatrixDot(ki,qmeas)*keV/hc : norm(qcalc)
		qmeas *= Q

		qcalc -= qmeas
		dq = norm(qcalc)
		dq *= (keV>0) ? weightHi : weightLo				// boost when I have a measured energy
		errQ += dq^2									// error is distance in Q space
		nQ += 1
	endfor
	errQ /= nQ

	Variable angleOffset=NumberByKey("angleOffset",noteStr,"=")	// rotation about incident beam (usually zero)
	angleOffset = numtype(angleOffset) ? 0 : angleOffset
#ifdef OLD_ORIENTATION_METHOD
	pixel2XYZ(d,(d.Nx-1)/2,(d.Ny-1)/2,kf)			// find detector center position, and try to move it to kf[0]==0
	normalize(kf)
	Variable errOrient =  kf[0]/norm(kf)
#else
	Variable errOrient = errVertAngle(d,angleOffset)	// error to use for orienting system along Y-axs (controls rotation about incident beam)
#endif
	errOrient = 1e4*(errOrient*errOrient)

	if (WhichListItem("printErrX",noteStr)>=0)
		//	print "norm(kf)=",norm(kf),"   ( kf[0]/norm(kf) )^2 = ",( kf[0]/norm(kf) )^2
		printf "errQ/nQ = %g,  nQ=%d,      errOrient = %g\r",errQ,nQ,errOrient
	endif
	return (errQ + errOrient)							// normalization makes E and Q equally important
End
//
Static Function errVertAngle(d,angleOffset)// calculate the error to use for orienting system along Y-axs
	STRUCT detectorGeometry, &d
	Variable angleOffset						// user supplied angle offset¡, when 0 detector face perp to y-axis, rotates about z-axis

	Variable tanSize=min(d.sizeX,d.sizeY)/abs(d.P[2])				// tan(angular size of detector)
	Make/N=3/D/FREE xyz
	pixel2XYZ(d,(d.Nx-1)/2,(d.Ny-1)/2,xyz)						// vector to detector center
	Variable tanCenter = sqrt(xyz[0]^2+xyz[1]^2)/abs(xyz[2])	// tan(angle to detector center)

	Variable angle										// angle between disired vertical and y^ (radian)
	if (tanSize>tanCenter)								// close to forward (or back) scattering
		// find direction along detector pixels (x or y) that is closest to y-axis, and make that the vertical
		Make/N=(3,3)/D/FREE rho
		rho[0][0] = {d.rho00, d.rho10, d.rho20}		// rotation matrix
		rho[0][1] = {d.rho01, d.rho11, d.rho21}
		rho[0][2] = {d.rho02, d.rho12, d.rho22}

		Make/N=(3,4)/D/FREE alongxy					// four vectors along the four directions ±x & ±y
		alongxy[0][0] = {1,0,0}
		alongxy[0][1] = {-1,0,0}
		alongxy[0][2] = {0,1,0}
		alongxy[0][3] = {0,-1,0}

		MatrixOP/FREE alongxy = rho x alongxy			// unit vector pointing along 4 directions of detector in beam line coords
		MatrixOP/FREE dotYs = row(alongxy,1)			// dot of each alongxy with y^
		WaveStats/Q/M=1 dotYs
		angle = atan2(alongxy[0][V_maxloc],alongxy[1][V_maxloc])	// angle between "most vertical direction" and yhat

	else													// far from forward scattering (the usual)
		angle = atan2(d.rho02,d.rho12)					// angle from y-axis to detector normal
		//Make/N=3/D/FREE xyzN						// xyzN = rho x [ (0 0 1) ], rho is pre-calculated from vector d.R
		//xyzN[0] = { d.rho02, d.rho12, d.rho22 }
		// printf "kf^ = %s,  angle from vertical = %g¡\r",vec2str(xyzN),angle*180/PI
	endif

	Variable errOrient = angle - (angleOffset*PI/180)
	// printf "errOrient = %g (only angle in radians)\r",errOrient
	return errOrient									// error to force correct orientation (radian)
End


Static Function CalculatedQvecEnergy(CalibrationList,m,qcalc)	// fills the calculated qcalc for point m in CalibrationList, returns calculated energy
	Wave CalibrationList
	Variable m
	Wave qcalc

	Wave rhoSample=rhoCalib							// rotation matrix computed from {rhox,rhoy,rhoz}
	Wave ki=ki											// preallocated waves (speeds things up)

	String noteStr = note(CalibrationList)
	Variable rhox=NumberByKey("rhox",noteStr,"="), rhoy=NumberByKey("rhoy",noteStr,"="), rhoz=NumberByKey("rhoz",noteStr,"=")
	RotationVec2Matrix(rhox,rhoy,rhoz,NaN,rhoSample)// Compute the sample rotation rhoSample. Apply to qhat's to make them line up with (px,py) pairs

	// compute ideal Q-vector from just rhoSample and the reciprocal lattice
	qcalc = CalibrationList[m][p+5]					// Q vector for this reflection in stardard orientation
	Variable Qlen = norm(qcalc)
	MatrixOp/O qcalc = rhoSample x qcalc				// rotate Q-vector by rhoSample

	Variable sintheta = -MatrixDot(ki,qcalc)/Qlen		// sin(theta) = -ki dot q^ 
	return hc * Qlen / (4*PI*sintheta)					// energy (keV)
End


Function TableCalibrationData(calib)
	Wave calib
//	String wList = WaveListClass("DetectorCalibrationList*","*","DIMS:2;MAXCOLS:10;MINCOLS:10")
	String wList = WaveListClass("DetectorCalibrationList*","*","DIMS:2;MAXCOLS:12;MINCOLS:10")
	if (!WaveExists(calib) && ItemsInList(wList)==1)
		Wave calib = $StringFromList(0,wList)
	elseif (!WaveExists(calib))
		String calibName
		Prompt calibName,"Calibration List wave",popup,wLIst
		DoPrompt "Calibration List Wave",calibName
		if (V_flag)
			return 1
		endif
		Wave calib = $calibName
	endif
	if (!WaveExists(calib))
		DoAlert 0,"cannot find CalibrationList wave"
		return 1
	endif

	String wName = "Table_"+NameOfWave(calib)
	String name = StringFromList(0,WinList(wName,";", "WIN:2"))
	if (strlen(name))
		DoWindow/F $name
		return 0
	endif
	Edit/K=1/N=$wName/W=(255,125,1144,599) calib.ld
	ModifyTable format(Point)=1,width(Point)=28,width(calib.l)=42,width(calib.d)=80
End



Function GraphAllCalibrationDetector()
//	String list = WaveListClass("DetectorCalibrationList*","*","DIMS:2;MAXCOLS:10;MINCOLS:10")
	String list = WaveListClass("DetectorCalibrationList*","*","DIMS:2;MAXCOLS:12;MINCOLS:10")
	Variable i,N=ItemsInList(list)
	for (i=0;i<N;i+=1)
		GraphCalibrationDetector($StringFromList(i,list))
	endfor
End

Function GraphCalibrationDetector(calib)
	Wave calib
//	String calibName,wList = WaveListClass("DetectorCalibrationList*","*","DIMS:2;MAXCOLS:10;MINCOLS:10")
	String calibName,wList = WaveListClass("DetectorCalibrationList*","*","DIMS:2;MAXCOLS:12;MINCOLS:10")
	if (!WaveExists(calib) && ItemsInList(wList)==1)
		Wave calib = $StringFromList(0,wList)
	elseif (!WaveExists(calib))
		Prompt calibName,"Calibration List wave",popup,wLIst
		DoPrompt "Calibration List Wave",calibName
		if (V_flag)
			return 1
		endif
		Wave calib = $calibName
	endif
	if (!WaveExists(calib))
		DoAlert 0,"cannot find CalibrationList wave"
		return 1
	endif

	String noteStr = note(calib)
	String detectorID=StringByKey("detectorID",noteStr,"=")

	calibName = NameOfWave(calib)
	String cName1 = calibName+"#1"
	String wName = "Graph_"+calibName
	String name = StringFromList(0,WinList(wName,";", "WIN:1"))
	if (strlen(name))
		DoWindow/F $name
		return 0
	endif
	Variable Nx=NumberByKey("Nx",noteStr,"="), Ny=NumberByKey("Ny",noteStr,"=")
	Variable sizeX=NumberByKey("sizeX",noteStr,"="), sizeY=NumberByKey("sizeY",noteStr,"=")
	Variable wid = round(max(sizeX,sizeY)/787.7)
	wid = min(wid,1000)
	Variable xoff=190
	xoff += stringmatch(detectorID,"PE0820, 476-1807") ? -190 : 0
	xoff += stringmatch(detectorID,"PE0820, 476-1850") ? 450 : 0

	Display /W=(xoff,40,xoff+wid,40+wid)/K=1/N=$wName calib[*][1] vs calib[*][0]
	ModifyGraph width={Aspect,Nx/Ny}
	AppendToGraph calib[*][1] vs calib[*][0]
	ModifyGraph gfMult=110, tick=2, mirror=1, minor=1, lowTrip=0.001, standoff=0, mode=3
	ModifyGraph rgb($calibName)=(0,0,0), marker($calibName)=8, marker($cname1)=19

	ImageStats/M=1/G={0,DimSize(calib,0)-1,8,8} calib
	Variable NkeV=V_npnts
	ImageStats/M=1/G={0,DimSize(calib,0)-1,11,11} calib
	Variable NdE=V_npnts
	if (NdE>2 && NdE>NkeV/2)
		Variable zrange = max(abs(V_max),abs(V_min))
		ModifyGraph zColor($cname1)={calib[*][11],-zrange,zrange,RedWhiteBlue256}
		ColorScale/C/N=textColorScale "ÆE  (eV)"
	else
		ModifyGraph zColor($cname1)={calib[*][8],*,*,Rainbow256}
		ColorScale/C/N=textColorScale "E  (keV)"
	endif
	ColorScale/C/N=textColorScale/F=0/S=3/M/B=1/A=RC/X=2.94/Y=0.98/E trace=$cname1
	ColorScale/C/N=textColorScale lblMargin=0, minor=1
	//	ColorScale/C/N=textColorScale width=12, fsize=10, lblMargin=0, minor=1

	ModifyGraph axOffset(left)=-2.7,axOffset(bottom)=-0.86
	Label left "Y (pixels)"
	Label bottom "X (pixels)"
	SetAxis left Ny+.5,-.5
	SetAxis bottom -.5,Nx+.5
	Cursor/P/S=2 A $calibName 0
	TextBox/C/N=textID/F=0/S=3/M/B=1/A=LT detectorID
	ShowInfo
End


Static Function RotationVec2Matrix(rx,ry,rz,theta,rho)	// compute the rotation matrix from a rotation vector
	Variable rx,ry,rz									// rotation axis
	Variable theta										// rotation angle (rad) (if NaN, then assume,   | {rx,ry,rz} | = theta (rad)
	Wave rho											// a (3,3) matrix to recieve rotation matrix

	if (DimSize(rho,0)!=3 && DimSize(rho,1)!=3)
		Abort "in rotationVec2Matrix(), mat must be (3,3)"
	endif

	Variable len = sqrt(rx*rx + ry*ry + rz*rz)	// length of {rx,ry,rz}
	theta = numtype(theta) ? len : theta				// if theta invalid, assume len is theta (rad)
	if (theta==0 || len==0)							// zero angle rotation just return the identity matrix
		rho = (p==q)
		return 0
	elseif (numtype(rx+ry+rz))					// totally invalid
		rho = NaN
		return 1
	endif
	Variable c=cos(theta), s=sin(theta)
	Variable c1 = 1-c
	rx /= len;	ry /= len;	rz /= len					// make |{rx,ry,rz}| = 1

	// this is the Rodrigues formula from:   http://mathworld.wolfram.com/RodriguesRotationFormula.html
	rho[0][0] = c + rx*rx*c1;			rho[0][1] = rx*ry*c1 - rz*s;		rho[0][2] = ry*s + rx*rz*c1
	rho[1][0] = rz*s + rx*ry*c1;		rho[1][1] = c + ry*ry*c1;			rho[1][2] = -rx*s + ry*rz*c1
	rho[2][0] = -ry*s + rx*rz*c1;	rho[2][1] = rx*s + ry*rz*c1;		rho[2][2] = c + rz*rz*c1
End

// ==================================== End of Optimization =====================================
// ==============================================================================================


// ==================================== Start of EPICS update ===================================
// ==============================================================================================

Function WriteDetectorGeo2EPICS(dNum)
	Variable dNum

	if (Exists("EPICS_put_PV_num")!=6)
		DoAlert 0,"EPICS not available"
		return 1
	elseif (dNum!=0 && dNum!=1 && dNum!=2)			// dNum must be 0, 1, or 2
		Prompt dNum,"detector number",popup,"0 (Orange);1 (Yellow);2 (Purple);"
		DoPrompt "Detector Number",dNum
		if (V_flag)
			return 1
		endif
		dNum -= 1
	endif
	STRUCT microGeometry g
	if (FillGeometryStructDefault(g))						//fill the geometry structure with current values
		DoAlert 0, "no geometry structure found, did you forget to set it?"
		return 1
	endif
	String si=num2istr(dNum),pv
	printf "Writing values for detector %d to EPICS\r",dNum

	if (putNum(EPICS_PREFIX+"Nx"+si,g.d[dNum].Nx))
		return 1
	endif
	if (putNum(EPICS_PREFIX+"Ny"+si,g.d[dNum].Ny))
		return 1
	endif
	if (putNum(EPICS_PREFIX+"sizeX"+si,(g.d[dNum].sizeX)/1000))
		return 1
	endif
	if (putNum(EPICS_PREFIX+"sizeY"+si,(g.d[dNum].sizeY)/1000))
		return 1
	endif

	if (putNum(EPICS_PREFIX+"Rx"+si,g.d[dNum].R[0]))
		return 1
	endif
	if (putNum(EPICS_PREFIX+"Ry"+si,g.d[dNum].R[1]))
		return 1
	endif
	if (putNum(EPICS_PREFIX+"Rz"+si,g.d[dNum].R[2]))
		return 1
	endif

	if (putNum(EPICS_PREFIX+"Px"+si,(g.d[dNum].P[0])/1000))
		return 1
	endif
	if (putNum(EPICS_PREFIX+"Py"+si,(g.d[dNum].P[1])/1000))
		return 1
	endif
	if (putNum(EPICS_PREFIX+"Pz"+si,(g.d[dNum].P[2])/1000))
		return 1
	endif

	if (putStr(EPICS_PREFIX+"timeMeasured"+si,g.d[dNum].timeMeasured))
		return 1
	endif
	if (putStr(EPICS_PREFIX+"geoNote"+si,g.d[dNum].geoNote))
		return 1
	endif
	if (putStr(EPICS_PREFIX+"detectorID"+si,g.d[dNum].detectorID))
		return 1
	endif
	if (putStr(EPICS_PREFIX+"distortionMapFile"+si,g.d[dNum].distortionMapFile))
		return 1
	endif
	return 0
End
//
Static Function putNum(pv,value)
	String pv
	Variable value

	String err
#if (exists("EPICS_put_PV_num")==6)
	err = EPICS_put_PV_num(pv,value)
#else
	printf "EPICS Unavailable"
	return 1
#endif
	if (strlen(err)<16)
		return 0
	endif
	printf "EPICS Error on PV='%s',    %s\r",pv,err
	err = "For '"+pv+"'  "+err
	DoAlert 0,err
	return 1
End
//
Static Function putStr(pv,str)
	String pv
	String str

	String err
#if (exists("EPICS_put_PV_num")==6)
	err = EPICS_put_PV_str(pv,str)
#else
	printf "EPICS Unavailable"
	return 1
#endif
	if (strlen(err)<16)
		return 0
	endif
	printf "EPICS Error on PV='%s',    %s\r",pv,err
	err = "At '"+pv+"'  "+err
	DoAlert 0,err
	return 1
End

// ==================================== End of EPICS update =====================================
// ==============================================================================================



Function findPxPyFromPixel(detNum,px,py)			// For Detector in the beam, find P from (px,py)
	Variable detNum									// detector number, usually 0
	Variable px,py										// pixel at center
	STRUCT microGeometry g
	FillGeometryStructDefault(g)						//fill the geometry structure with current values
	STRUCT detectorGeometry d
	Variable detNumMax = (g.Ndetectors)-1

	if (!(detNum>=0 && px>=0 && py>=0))			// need to prompt
		String str,popStr
		sprintf str, "Detector Number [0, %d]", detNumMax
		detNum += 1
		Prompt detNum, str, popup, expandRange("0-"+num2str(detNumMax),";")
		Prompt px, "X-pixel where beam hits"
		Prompt py, "Y-pixel where beam hits"
		DoPrompt "Incident Beam Pixel", detNum, px,py
		if (V_flag)
			return 1
		endif
		detNum -= 1
		printf "¥findPxPyFromPixel(%g, %.6g, %.6g)\r",detNum,px,py
	endif

	String alertStr=""
	if (detNum != limit(round(detNum),0,detNumMax))
		sprintf alertStr "Detector number \"%g\" is NOT an integer in the range [0, %d]", detNum,detNumMax
	else
		d = g.d[0]
		if (!(d.used))
			sprintf alertStr "Detector \"%g\" is not used", detNum
		endif
	endif
	if (strlen(alertStr))
		DoAlert 0, alertStr
	endif

	Variable xp,yp										// x' and y' (requiring z'=0), detector starts centered on origin and perpendicular to z-axis
	xp = (px - 0.5*(d.Nx-1)) * d.sizeX/d.Nx	// (x' y' z'), position on detector
	yp = (py - 0.5*(d.Ny-1)) * d.sizeY/d.Ny
	printf "for Detector %g (id=%s) with center at pixel [%.1f, %.1f], use P = {%.6g, %.6g, %.6g}\r",detNum,d.detectorID,px,py,-xp/1000,-yp/1000,(d.P[2])/1000
	return 0
End




// ===================================== End of Calibration =====================================
// ==============================================================================================





// ==============================================================================================
// ================================ Start of Optimization Tests =================================

Function testOneOptimize(noise)
	Variable noise
	noise = (noise>=0) ? noise : 0.25					// default to 0.25
	SetRandomSeed 1										// use this if calling testOptimize() repeatedly to test
	Variable oneErr
	String list
	Variable timer=startMSTimer
	list = MakeFakeCalibration3Detectors(noise)
	Wave c0=$StringFromList(0,list), c1=$StringFromList(1,list), c2=$StringFromList(2,list)
	oneErr =  OptimizeAll(c0,c1,c2)
	printf "the error was  %.2g,  and the total execution time for optimizations was %s\r",oneErr,Secs2Time(stopMSTimer(timer)*1e-6,5,1)
End

Function testManyOptimize(N,noise)
	Variable N
	Variable noise
	noise = (noise>=0) ? noise : 0.25					// default to 0.25
	if (!(N>0))
		N = (N>0) ? N : 3									// default to 3
		Prompt N,"# of loops"
		Prompt noise "noise added to pixels"
		DoPrompt "test loops",N,noise
		if (V_flag)
			return 1
		endif
		printf "testManyOptimize(%g,%g)\r",N,noise
	endif

	Variable i
	SetDefaultGeo2Reference()	
	DoAlert 1, "Shift around the starting point?"
	if (V_flag==1)										// shift around the stasrting point some
		STRUCT microGeometry g	
		FillGeometryStructDefault(g)					//fill the geometry structure with current values
		Variable dP=10e3, dR=0.2
dP=1e3
dR=0.1

		for (i=0;i<MAX_Ndetectors;i+=1)
			if (g.d[i].used)
				g.d[i].R[0] += dR/2;		g.d[i].R[1] += dR;		g.d[i].R[2] += 1.5*dR	// and add imprecision so we do not start at the answer, approximately 10¡ & 10mm
				g.d[i].P[0] += dP;		g.d[i].P[1] -= dP;		g.d[i].P[2] -= 10*dP
			endif
		endfor
		UpdateDefaultGeometryStruct(g)					// save the changes
	endif

	SetRandomSeed 1										// use this if calling testOptimizeAll() repeatedly to test
	Make/N=(N)/O manyErrs=NaN
	String list
	Variable timer=startMSTimer
	for (i=0;i<N;i+=1)
		list = MakeFakeCalibration3Detectors(noise)
		Wave c0=$StringFromList(0,list), c1=$StringFromList(1,list), c2=$StringFromList(2,list)
		manyErrs[i] = OptimizeAll(c0,c1,c2)
	endfor
	WaveStats/Q manyErrs
	Variable fwhm = 2*noise*sqrt(2*ln(2))
	printf "the error range was [%.2g, %.2g], <error>=%g,    using a noise of  sigma=%g (FWHM=%.2g),  and the total execution time for the %d optimizations was %s\r",V_min,V_max,V_avg,noise,fwhm,N,Secs2Time(stopMSTimer(timer)*1e-6,5,1)
End


//Function ShiftStartingPoints()							// shitfs the starting point for the fit by a medium amount
//	STRUCT microGeometry g
//	FillGeometryStructDefault(g)						//fill the geometry structure with current values
//
//	Variable i, dP=1e3, dR=0.2
//	for (i=0;i<g.Ndetectors;i+=1)
//		g.d[i].R[0] += dR/2;	g.d[i].R[1] += dR;		g.d[i].R[2] += 1.5*dR		// and add imprecision so we do not start at the answer, approximately 10¡ & 10mm
//		g.d[i].P[0] += dP;		g.d[i].P[1] -= dP;		g.d[i].P[2] -= 10*dP
//	endfor
//	UpdateDefaultGeometryStruct(g)						// save the changes
//End

Function/T MakeFakeCalibration3Detectors(noise)
	Variable noise											// gaussian noise added to each test pixel position (sigma in pixels)
	Variable printing = (!strlen(GetRTStackInfo(2))) || stringmatch(GetRTStackInfo(2),"testOneOptimize")
	if (!(noise>=0))
		noise = (noise>=0) ? noise : 0.25				// default to 0.25
		Prompt noise "noise added to pixels"
		DoPrompt "pixel noise",noise
		if (V_flag)
			return ""
		endif
		printing = 1
		printf "MakeFakeCalibrationAllDetectors(%g)\r",noise
	endif
	if (!(noise>=0))
		DoAlert 0, "bad inputs"
		return ""
	endif
	if (noise && printing)
		printf "Æpixel has sigma=%g, and ÆeV=±%geV:\r",noise,noise/2
	elseif (printing)
		printf "No noise:\r"
	endif

	Variable dNum											// detector number {0,1,2}
	String cName,cNameLIst=""
	for (dNum=0; dNum<3; dNum+=1)
		cNameList += MakeFakeCalibrationData(noise,dNum,NaN)	// data to use for optimization
		cNameList += ";"
	endfor

	if (printing)
		print "made the following:  ",cNameList
	endif
	return cNameList
End
//
Static Function/T MakeFakeCalibrationData(noise,dNum,N)
	Variable noise											// signa of noise (pixels)
	Variable dNum											// detector number {0,1,2}
	Variable N												// number of spots to make
	dNum = limit(round(dNum),0,2)					// only know about detectors 0,1,2
	N = numtype(N) ? 25 : N
	N = limit(N,1,25)
	Variable printing = (!strlen(GetRTStackInfo(2)))

	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))					//fill the lattice structure with test values
		DoAlert 0, "no crystal structure found"
		return ""
	endif

	Make/N=(3,3)/O/D recip, recip0, rhoSample
	recip0[0][0]=xtal.as0;	recip0[0][1]=xtal.bs0;	recip0[0][2]=xtal.cs0	// reciprocal lattice { a*[], b*[], c*[] }
	recip0[1][0]=xtal.as1;	recip0[1][1]=xtal.bs1;	recip0[1][2]=xtal.cs1
	recip0[2][0]=xtal.as2;	recip0[2][1]=xtal.bs2;	recip0[2][2]=xtal.cs2

	// sample rotation
	Make/N=3/O/D axisSampleFake
	axisSampleFake={-135*PI/180,0,0}
	Variable angle = norm(axisSampleFake)*180/PI	// sample rotation angle (deg)
	RotationVec2Matrix(axisSampleFake[0],axisSampleFake[1],axisSampleFake[2],NaN,rhoSample)
	MatrixOp/O recip = rhoSample x recip0				// rotate recip0
	KillWaves/Z rhoSample

	STRUCT microGeometry g
	GeoReferenceOrientation(g)							// use reference geometry
//	FillGeometryStructDefault(g)						//fill the geometry structure with current values
	if (printing)
		print "for detector",dNum
		printDetector(g.d[dNum])
		printf "sample rotation axis = {%g,%g,%g},  a rotation angle of  %g¡\r",axisSampleFake[0],axisSampleFake[1],axisSampleFake[2],angle
	endif

	// desired pixels for spots on detector (a sort of uniformly spaced array of spots that covers most of the detector)
	String pkList=""
	if (dNum==0)
		pkList="1023.5,1023.5;140,1180;1907,1180;1023.5,166;1023.5,1881;140,86;1907,86;140,1961;1907,1961;631,544;"
		pkList += "602,1501;1513,1492;1542,514;197,598;265,1623;612,1852;1494,1824;1793,1580;1827,730;1480,163;578,191;"
		pkList += "587,1000;1007,1426;1431,1014;1026,580;"
	else
		pkLIst = "511.5,511.5;70,590;940,590;512,83;512,941;70,43;920,63;70,981;950,981;316,272;301,751;757,746;771,257;"
		pkLIst += "99,299;133,812;306,926;747,912;897,790;914,365;740,82;289,96;294,500;504,713;716,507;513,290;"
	endif
	Variable i
	for (i=N;i<25;i+=1)									// remove extra spots from pkList
		pkList = RemoveListItem(N,pkList)
	endfor

	//	CalibrationList[][0,1]	= (px,py) are measured on image
	//	CalibrationList[][2-4]	= hkl are hkl for each measured spot
	//	CalibrationList[][5-7]	= Q-vector, calculated from knowing hkl and the standard reciprocal lattice (NOT rotated from standard orientation)
	//	CalibrationList[][8]	= keV, measured energy of spot
	//	CalibrationList[][9]	= theta (deg), obtained from calculated d(hkl) and measured energy using: lambda = 2d sin(theta)
	//	CalibrationList[][10]	= keV calculated from sample orientation and known lattice (does not use measured spot position)
	//	CalibrationList[][11]	= ÆE (eV),  ( CalibrationList[][8] - CalibrationList[][10] ) *1e3
	String name="CalibrationListTest"+num2istr(dNum)
//	Make/N=(N,10)/O/D $name=NaN					// a list of measured reflections used to calibrate the detector, {px,py,qx,qy,yz,keV,theta}
	Make/N=(N,12)/O/D $name=NaN					// a list of measured reflections used to calibrate the detector, {px,py,qx,qy,yz,keV,theta}
	Wave cList = $name
	SetDimLabel 1,0,px,cList;		SetDimLabel 1,1,py,cList											// measued pixel
	SetDimLabel 1,2,H,cList;		SetDimLabel 1,3,K,cList;		SetDimLabel 1,4,L,cList					// (hkl)
	SetDimLabel 1,5,qx_hkl,cList;	SetDimLabel 1,6,qy_hkl,cList;	SetDimLabel 1,7,qz_hkl,cList	// qvector from (hkl)
	SetDimLabel 1,8,keV,cList														// measured keV
	SetDimLabel 1,9,theta_keV_hkl,cList											// theta calculated from known d(hkl) and measured keV
	SetDimLabel 1,10,calculated_keV,cList											// energy calculated from sample orientation (no measured parameters)
	SetDimLabel 1,11,deltaE_eV,cList												// measured - calculated energy (eV)
	String noteStr = ""
	noteStr=ReplaceStringByKey("waveClass",noteStr,"DetectorCalibrationList"+SelectString(dNum,"0",""),"=")
	noteStr = ReplaceStringByKey("detectorID",noteStr,g.d[dNum].detectorID,"=")// detector ID
	noteStr = ReplaceNumberByKey("Nx",noteStr,g.d[dNum].Nx,"=")				// number of un-binned pixels in whole detector
	noteStr = ReplaceNumberByKey("Ny",noteStr,g.d[dNum].Ny,"=")
	noteStr = ReplaceNumberByKey("sizeX",noteStr,g.d[dNum].sizeX,"=")		// outside size of detector (micron)
	noteStr = ReplaceNumberByKey("sizeY",noteStr,g.d[dNum].sizeY,"=")
	noteStr = ReplaceNumberByKey("rhox",noteStr,axisSampleFake[0],"=")	// rotation vector for sample
	noteStr = ReplaceNumberByKey("rhoy",noteStr,axisSampleFake[1],"=")
	noteStr = ReplaceNumberByKey("rhoz",noteStr,axisSampleFake[2],"=")
	Note/K cList, noteStr
	KillWaves/Z axisSampleFake

	Make/N=3/O/D hkl, ki={0,0,1}, kf, qhat, qvec
	Variable dot, Q, keV
	Variable px,py, h,k,l, hmax, m
	for (i=0,m=0; i<N; i+=1)
		sscanf StringFromLIst(i,pkList),"%g,%g",px,py// first find integral hkl closest to (px,py)
		pixel2XYZ(g.d[dNum],px,py,kf)					// get kf from pixel on detector
		normalize(kf)
		qvec = kf - ki
		MatrixOp/O hkl = Inv(recip) x qvec
		ClosestHKL(hkl)	// in and out can be the same wave
		if (numtype(sum(hkl)))							// invalid point
			continue
		endif
		h = hkl[0]
		k = hkl[1]
		l = hkl[2]
		lowestAllowedHKL(h,k,l)

		// got a valid hkl from the (px,py), use values that are exact for hkl for another point in the cList[m][*]
		hkl = {h,k,l}										// re-compute using exact hkl
		MatrixOp/O qhat =  recip x hkl
		Q = normalize(qhat)
		dot = -MatrixDot(ki,qhat)						// = sin(thetaBragg)
		keV = Q*hc/(4*PI*dot)
		kf = ki +2*dot*qhat
		XYZ2pixel(g.d[dNum],kf,px,py)					// find pixel where kf hits detector 0
		if (keV>99 || px<0 || px>=g.d[dNum].Nx || py<0 || py>=g.d[dNum].Ny)
			continue
		endif

		// add some noise the pretend data to make it more realistic
		px += gnoise(noise)								// add some noise to the pixel position (make things more realistic)
		py += gnoise(noise)
		keV += noise>0 ? enoise(noise/2000) : 0		// when noise present, add ~1eV of noise

		// recalculate qhat to be STANDARD orientation, cList[][5-7] does not need to know actual orientations, only angle between qhat pairs is used
		MatrixOp/O qvec =  recip0 x hkl

		cList[m][0] = px;		cList[m][1] = py											// measured pixel position
		cList[m][2] = h;		cList[m][3] = k;		cList[m][4] = l					// (hkl)
		cList[m][5] = qvec[0];	cList[m][6] = qvec[1];	cList[m][7] = qvec[2]	// calculated from hkl (not measured)
		cList[m][8] =keV																	// keV,  from Q = 4¹ sin(theta)/lambda == Q = 4¹ sin(theta) * (E/hc)
		cList[m][9] = asin(hc*Q/(4*PI*keV))*180/PI								// theta (degree),  calculated from Q = 4¹ sin(theta) * (E/hc)
		cLIst[m][10] = NaN
		cList[m][11] = NaN
		m += 1
	endfor
	N = m
	if (printing)
		print "N=",N
	endif
//	Redimension/N=(N,10) cList
	Redimension/N=(N,12) cList

	if (N>1)
		ImageStats/G={0,N-1,8,8} cList				// pretend that did not measure spot with highest energy
		cList[V_maxRowLoc][8] = NaN
		cList[V_maxRowLoc][9] = NaN
	endif
	KillWaves/Z recip,recip0, qhat, qvec, ki, kf, hkl
	return GetWavesDataFolder(cList,2)
End
//
// find hkl that is paralllel to the input and allowed
Static Function ClosestHKL(hkl)	// in and out can be the same wave
	Wave hkl												// real numbers on intput, set to h,k,l on output

	Variable xin=hkl[0],yin=hkl[1],zin=hkl[2]		// save input values
	hkl = abs(hkl)
	Variable i, maxVal = WaveMax(hkl)
	Variable ibest=0, dotMax=-Inf, dot
	for (i=1;i<=25;i+=1)
		hkl = {xin,yin,zin}
		hkl = round(hkl*i/maxVal)
		dot = (hkl[0]*xin + hkl[1]*yin + hkl[2]*zin)/norm(hkl)
		if (dot > dotMax)
			ibest = i
			dotMax = dot
		endif
	endfor
	hkl = {xin,yin,zin}
	hkl = round(hkl*ibest/maxVal)
	hkl = hkl[p]==0 ? 0 : hkl[p]							// avoid "-0"
	return 0
End



////
//Static Function RemoveCommonFactors(h,k,l)
//	Variable &h,&k,&l
//	h = round(h)
//	k = round(k)
//	l = round(l)
//	Variable i, maxh = max(max(abs(h),abs(k)),abs(l))
//	Variable hh,kk,ll
//	for (i=maxh;i>1;i-=1)
//		if (mod(h/i,1)==0 && mod(k/i,1)==0 && mod(l/i,1)==0)
//			h /= i
//			k /= i
//			l /= i
//			return i
//		endif
//	endfor
//	return 0
//End
////Function test(h,k,l)
////	Variable h,k,l
////	RemoveCommonFactors(h,k,l)
////	print h," ",k," ",l
////End
//

// ================================= End of Optimization Tests ==================================
// ==============================================================================================




// ==============================================================================================
// =============================== Start of Calibration set panel ===============================

Function/T FillCalibrationParametersPanel(strStruct,hostWin,left,top)
	String strStruct									// optional passed value of xtal structure, this is used if passed
	String hostWin										// name of home window
	Variable left, top									// offsets from the left and top

	SetWindow kwTopWin,userdata(CalibrationPanelName)=hostWin+"#CalibrationPanel"
	NewPanel/K=1/W=(left,top,left+221,top+445+30)/HOST=$hostWin
	ModifyPanel frameStyle=0, frameInset=0
	RenameWindow #,CalibrationPanel

	Button buttonEnterMeasuredEnergies,pos={29,5},size={180,20},proc=detectorCalibration#CalibrationButtonProc,title="Enter Measured Energies"
	Button buttonEnterMeasuredEnergies,help={"Set the Calibration Data"}

	Button buttonMakeCalibData,pos={29,5+35},size={160,20},proc=detectorCalibration#CalibrationButtonProc,title="Set Calibration Data"
	Button buttonMakeCalibData,help={"Set the Calibration Data"}

	Button buttonCalibOptimizeAll3,pos={29,40+35},size={160,20},proc=detectorCalibration#CalibrationButtonProc,title="Optimize All"
	Button buttonCalibOptimizeAll3,help={"Run Optimization on some set of detectors"}

	Button buttonCalibGraph1,pos={29,75+35},size={160,20},proc=detectorCalibration#CalibrationButtonProc,title="Graph 1 Detector"
	Button buttonCalibGraph1,help={"Graph Results for one detector"}

	Button buttonCalibGraphAll,pos={29,100+35},size={160,20},proc=detectorCalibration#CalibrationButtonProc,title="Graph All Detectors"
	Button buttonCalibGraphAll,help={"Graph results for all 3 detectors"}

	Button buttonCalibTable,pos={29,135+35},size={160,20},proc=detectorCalibration#CalibrationButtonProc,title="Table of Calib Data"
	Button buttonCalibTable,help={"Show table of Calibration Data"}

	Button buttonCalibWireOrigin,pos={29,170+35},size={160,20},proc=detectorCalibration#CalibrationButtonProc,title="Calc. Wire Origin"
	Button buttonCalibWireOrigin,help={"Calculate Wire Origin from position of a spot"}

	Button buttonCalibWrite2EPICS,pos={29,205+35},size={160,20},proc=detectorCalibration#CalibrationButtonProc,title="Write to EPICS"
	Button buttonCalibWrite2EPICS,help={"Write Optimized results to EPICS"}

	EnableDisableCalibControls(hostWin+"#CalibrationPanel")
	return "#CalibrationPanel"
End
//
Static Function EnableDisableCalibControls(win)				// here to enable/disable
	String win												// window (or sub window) to use
	Variable d

	d = strlen(WaveListClass("FittedPeakList*","*",""))<1 ? 2 : 0
	Button buttonEnterMeasuredEnergies,win=$win,disable=d

	d = strlen(WaveListClass("IndexedPeakList*","*",""))<1 ? 2 : 0
	Button buttonMakeCalibData,win=$win,disable=d

	d = strlen(WaveListClass("DetectorCalibrationList*","*",""))<1 ? 2 : 0
	Button buttonCalibOptimizeAll3,win=$win,disable=d
	Button buttonCalibGraph1,win=$win,disable=d
	Button buttonCalibGraphAll,win=$win,disable= d
	Button buttonCalibTable,win=$win,disable=0

	d = Exists("EPICS_put_PV_num")==6 ? 0 : 2
	Button buttonCalibWrite2EPICS,win=$win,disable=d

	Button buttonCalibWireOrigin,win=$win,disable=0
End
//
Static Function CalibrationButtonProc(B_Struct) : ButtonControl
	STRUCT WMButtonAction &B_Struct
	if (B_Struct.eventCode != 2)
		return 0
	endif
	String ctrlName=B_Struct.ctrlName

	STRUCT microGeometry g
	FillGeometryStructDefault(g)

	if (stringmatch(ctrlName,"buttonEnterMeasuredEnergies"))
		EnterMeasuredEnergies("")
	elseif (stringmatch(ctrlName,"buttonMakeCalibData"))
		MakeCalibData(NaN)
	elseif (stringmatch(ctrlName,"buttonCalibOptimizeAll3"))
		OptimizeAll($"",$"",$"")
	elseif (stringmatch(ctrlName,"buttonCalibGraph1"))
		GraphCalibrationDetector($"")
	elseif (stringmatch(ctrlName,"buttonCalibGraphAll"))
		GraphAllCalibrationDetector()
	elseif (stringmatch(ctrlName,"buttonCalibTable"))
		TableCalibrationData($"")
	elseif (stringmatch(ctrlName,"buttonCalibWrite2EPICS"))
		WriteDetectorGeo2EPICS(NaN)
	elseif (stringmatch(ctrlName,"buttonCalibWireOrigin"))
		calcWireOrigin(NaN,NaN,NaN,NaN,NaN)
	endif
	EnableDisableCalibControls(GetUserData("microPanel","","CalibrationPanelName"))
End

// ================================ End of Calibration set panel ================================
// ==============================================================================================




// ==============================================================================================
// =============================== Start of Wire Calibration Fit ================================


Function/WAVE MakeWireFittingTable(PixelIntensities,[printIt])
	Wave PixelIntensities
	Variable printIt
	printIt = ParamIsDefault(printIt) ? NaN : printIt
	printIt = numtype(printIt) ? 0 : !(!printIt)

	if (!WaveExists(PixelIntensities))
		String wList = reverseList(WaveListClass("PixelIntensities*","*","DIMS:2"))
		if (ItemsInList(wList)<2)
			Wave PixelIntensities = $StringFromList(0,wList)
		else
			String wName
			Prompt wName,"Pixel Intensity Wave",popup,wList
			DoPrompt "Pixel Intensities",wName
			if (V_flag)
				return $""
			endif
			Wave PixelIntensities = $wName
		endif
		printIt = 1
	endif
	if (printIt)
		printf "MakeWireFittingTable(%s",NameOfWave(PixelIntensities)
		if (!ParamIsDefault(printIt))
			printf ", printIt=1"
		endif
		printf ")\r"
	endif
	if (!WaveExists(PixelIntensities))
		return $""
	endif

	String wnote=note(PixelIntensities)
	Variable NpIN=DimSize(PixelIntensities,1), i
	String pxList=StringByKey("PixelPositions",wnote,"="), str
	String X2list=StringByKey("X2",wnote,"="), Y2list=StringByKey("Y2",wnote,"="), Z2list=StringByKey("Z2",wnote,"=")
	Variable detNum=detectorNumFromID(StringByKey("detectorID",wnote,"="))

	Make/N=(NpIN)/FREE/D iTurnOffs=NaN
	Variable N=DimSize(PixelIntensities,0), i0,i1, NpGood=0
	Make/N=(N)/D/FREE vals
	for (i=0;i<NpIN;i+=1)
		vals = PixelIntensities[p][i]
		iTurnOffs[i] = FitWireIntensityProfile(vals)
		if (numtype(iTurnOffs[i]))
			print " *** Could NOT find level for point",i
		else
			NpGood += 1									// increment good ones
		endif
 	endfor
	Make/N=(NpGood,8)/D/O PixelTurnOffs=NaN	// px,py,x,y,z
	wnote = RemoveByKey("X2",wnote,"=")
	wnote = RemoveByKey("Y2",wnote,"=")
	wnote = RemoveByKey("Z2",wnote,"=")
	wnote = ReplaceStringByKey("waveClass",wnote,"PixelTurnOffs","=")
	wnote = ReplaceNumberByKey("detNum",wnote,detNum,"=")
	Note/K PixelTurnOffs, wnote
	SetDimLabel 1,0,px0,PixelTurnOffs		;	SetDimLabel 1,1,py0,PixelTurnOffs
	SetDimLabel 1,2,wireX,PixelTurnOffs	;	SetDimLabel 1,3,wireY,PixelTurnOffs
	SetDimLabel 1,4,wireZ,PixelTurnOffs	;	SetDimLabel 1,5,ip,PixelTurnOffs
	SetDimLabel 1,6,depth0,PixelTurnOffs	;	SetDimLabel 1,7,depth1,PixelTurnOffs

	Variable m, ip
	for (i=0,m=0; i<NpIN; i+=1)
		ip = iTurnOffs[i]
		if (numtype(ip))
			continue
		endif
		str = StringFromList(i,pxList,",")
		PixelTurnOffs[m][0] = str2num(StringFromList(0,str,":"))
		PixelTurnOffs[m][1] = str2num(StringFromList(1,str,":"))
		PixelTurnOffs[m][2] = str2num(StringFromList(ip,X2list,","))
		PixelTurnOffs[m][3] = str2num(StringFromList(ip,Y2list,","))
		PixelTurnOffs[m][4] = str2num(StringFromList(ip,Z2list,","))
		PixelTurnOffs[m][5] = iTurnOffs[i]
		m += 1
	endfor

	Wave depths = CalcDepthsFromPixelTurnOffs(PixelTurnOffs,$"")
	PixelTurnOffs[][6] = depths[p]			// save initial depths
	PixelTurnOffs[][7] = NaN
	//	Duplicate/O depths, depthsView
	WaveStats/Q depths
	printf " at start,  <depth> = %g,   range = [%g, %g],   std. dev = %g\r",V_avg,V_min,V_max,V_sdev
	return PixelTurnOffs
End
//
Static Function FitWireIntensityProfile(wy)
	// used by MakeWireFittingTable() to find where wire clips spot
	Wave wy
	Variable hw=50
	if (!WaveExists(wy))
		return NaN
	endif

	SetScale/P x,0,1,"" wy
	Variable N=DimSize(wy,0)
	Variable i0=2, i1=N-2
	WaveStats/M=1/Q/R=[i0,i1] wy
	i1 = V_minloc
	WaveStats/M=1/Q/R=[i0,i1] wy
	FindLevel/B=3/EDGE=2/P/Q/R=[i0,i1] wy, (V_min+V_max)/2
	V_LevelX = round(V_LevelX)

	i0 = limit(V_LevelX-hw,2,V_LevelX)
	i1 = limit(V_LevelX+hw,i0+5,N-1)
	//	print i0,"  ",i1
	if (i1-i0 < 10)
		return NaN
	endif
	Variable V_fitOptions=4, V_FitError=0
	CurveFit/Q/NTHR=0 Sigmoid wy[i0,i1]
	if (V_FitError)
		return NaN
	endif
	Wave W_coef=W_coef, W_sigma=W_sigma
	Variable mid=W_coef[2], rate=W_coef[3]
	KillWaves/Z W_coef, W_sigma
	//	print "rate = ",rate,"  mid =",mid

	return rate<5 ? mid : NaN		// if rate>5, then step is too broad
	return mid
End



Function/WAVE OptimizeWireGeometry(PixelTurnOffs,Nfit,[printIt])
	Wave PixelTurnOffs
	Variable Nfit
	Variable printIt
	Nfit = numtype(Nfit) ? NaN : round(Nfit)
	Nfit = Nfit==limit(Nfit,1,4) ? Nfit : NaN
	printIt = ParamIsDefault(printIt) ? NaN : printIt
	printIt = numtype(printIt) ? 0 : !(!printIt)

#ifdef OptimizeWireAxis
	String fitList="wire.origin(Y);wire.origin(Y,Z);wire.origin(Y,Z), wire.axis(Z);wire.origin(Y,Z), wire.axis(Y,Z);"
#else
	String fitList="wire.origin(Y);wire.origin(Y,Z);wire.origin(Y,Z), wire.R(Z);wire.origin(Y,Z), wire.R(Y,Z);"
#endif
	if (!WaveExists(PixelTurnOffs) || numtype(Nfit))
		String wList = reverseList(WaveListClass("PixelTurnOffs","*","DIMS:2"))
		if (ItemsInList(wList)<2 && numtype(Nfit)==0)
			Wave PixelTurnOffs = $StringFromList(0,wList)
		else
			Nfit = numtype(Nfit) ? 4 : Nfit
			String wName, fitStr
			Prompt wName,"Pixel Turn Off Wave",popup,wList
			Prompt Nfit, "What to Optimize", popup,fitList
			DoPrompt "Pixel Turn Offs",wName, Nfit
			if (V_flag)
				return $""
			endif
			Wave PixelTurnOffs = $wName
		endif
		printIt = 1
	endif
	if (!WaveExists(PixelTurnOffs))
		return $""
	endif
	if (printIt)
		printf "OptimizeWireGeometry(%s,%g",NameOfWave(PixelTurnOffs),Nfit
		if (!ParamIsDefault(printIt))
			printf ", printIt=1"
		endif
		printf ")\t\t// Optimizing {%s}\r",StringFromList(Nfit-1,fitList)
	endif
	if (!(Nfit==limit(Nfit,1,4)))
		return $""
	endif

	STRUCT microGeometry g
	FillGeometryStructDefault(g)
	Make/N=4/D/FREE xWave, typXWave={1,1,0.001,0.001}
	xWave[0] = g.wire.origin[1]
	xWave[1] = g.wire.origin[2]
#ifdef OptimizeWireAxis
	xWave[2] = g.wire.axis[1]
	xWave[3] = g.wire.axis[2]
#else
	xWave[2] = g.wire.R[1]
	xWave[3] = g.wire.R[2]
#endif
	Redimension/N=(Nfit) xWave,typXWave

	Wave depths = CalcDepthsFromPixelTurnOffs(PixelTurnOffs,xWave)
	WaveStats/Q depths
	Variable errStart=V_sdev
	printf "Before Optimize:  <depth> = %g,   range = [%g, %g],  Æ=%g\t\txWave = %s\r",V_avg,V_min,V_max,errStart,vec2str(xWave)
	String wnote=note(PixelTurnOffs)
	wnote = ReplaceStringByKey("OptimizeXvalues0",wnote,vec2str(xWave),"=")
	wnote = ReplaceNumberByKey("OptimizeXvaluesErr0",wnote,errStart,"=")

	Optimize/Q/A=0/R=typXWave/X=xWave/Y=0.1 detectorCalibration#depthError,PixelTurnOffs
	if (V_flag)
		printf "\r\tERROR -- V_flag = %g,  V_OptTermCode = %g,  V_OptNumIters=%g,  V_OptNumFunctionCalls=%g\r",V_flag,V_OptTermCode,V_OptNumIters,V_OptNumFunctionCalls
		return $""
	endif
	Variable errEnd = V_min

	Wave depths = CalcDepthsFromPixelTurnOffs(PixelTurnOffs,xWave)
	PixelTurnOffs[][7] = depths[p]			// save final depths
	WaveStats/Q depths
	printf "After Optimize:  <depth> = %g,   range = [%g, %g],  Æ=%g\t\txWave = %s\r",V_avg,V_min,V_max,errEnd,vec2str(xWave)
	//	Wave W_OptGradient=W_OptGradient
	//	print "\tGradient =",vec2str(W_OptGradient)
	KillWaves/Z W_OptGradient

	print " "
	printf "\t\twire.origin(Y) changed from  %g  ->  %g\r",g.wire.origin[1], xWave[0]
	if (Nfit>1)
		printf "\t\twire.origin(Z) changed from  %g  ->  %g\r",g.wire.origin[2], xWave[1]
	endif
#ifdef OptimizeWireAxis
	if (Nfit>2)
		printf "\t\twire.axis(Y) changed from  %g  ->  %g\r",g.wire.axis[1], xWave[2]
	endif
	if (Nfit>3)
		printf "\t\twire.axis(Z) changed from  %g  ->  %g\r",g.wire.axis[2], xWave[3]
	endif
#else
	if (Nfit>2)
		printf "\t\twire.R(Y) changed from  %g  ->  %g\r",g.wire.R[1], xWave[2]
	endif
	if (Nfit>3)
		printf "\t\twire.R(Z) changed from  %g  ->  %g\r",g.wire.R[2], xWave[3]
	endif
#endif
	printf "\terror changed from  %g  ->  %g\r",errStart, errEnd

	wnote = ReplaceStringByKey("OptimizeXvalues1",wnote,vec2str(xWave),"=")
	wnote = ReplaceNumberByKey("OptimizeXvaluesErr1",wnote,errEnd,"=")
	Note/K PixelTurnOffs, wnote
	return xWave
End
//
Static Function depthError(PixelTurnOffs,xw)		// function called by Optimize
	Wave PixelTurnOffs
	Wave xw			// things that get adjusted

	Wave depths = CalcDepthsFromPixelTurnOffs(PixelTurnOffs,xw)
	WaveStats/Q depths
	return V_sdev
End
//
Static Function/WAVE CalcDepthsFromPixelTurnOffs(PixelTurnOffs,xw)
	Wave PixelTurnOffs
	Wave xw			// things that get adjusted

	Variable detNum=NumberByKey("detNum",note(PixelTurnOffs),"=")
	STRUCT microGeometry g
	FillGeometryStructDefault(g)

	if (WaveExists(xw))
		Variable Nfit=DimSize(xw,0)
		g.wire.origin[1] = xw[0]		// wire Origin Y
		g.wire.origin[2] = xw[1]		// wire Origin Z

#ifdef OptimizeWireAxis
		if (Nfit>2)
			g.wire.axis[1] = xw[2]		// wire axis[Y]
			if (Nfit>3)
				g.wire.axis[2] = xw[3]	// wire axis[Y]
			endif
		endif
#else
		if (Nfit>2)
			g.wire.R[1] = xw[2]			// wire stage Ry
			if (Nfit>3)
				g.wire.R[2] = xw[3]		// wire stage Rz
			endif
		endif
#endif
		GeometryUpdateCalc(g)
	endif

	Variable px,py,i, Np=DimSize(PixelTurnOffs,0)
	Make/N=3/D/FREE wxyz, dxyz
	Make/N=(Np)/FREE/D depths
	for (i=0;i<Np;i+=1)
		wxyz = PixelTurnOffs[i][2+p]
		px = PixelTurnOffs[i][0]
		py = PixelTurnOffs[i][1]
		pixel2XYZ(g.d[detNum],px,py,dxyz)
		depths[i] = PixelxyzWire2depth(g,dxyz,wxyz,1)
	endfor
	return depths
End



Function/WAVE TableOfPixelTurnOffs(PixelTurnOffs,[printIt])
	Wave PixelTurnOffs
	Variable printIt
	printIt = ParamIsDefault(printIt) ? NaN : printIt
	printIt = numtype(printIt) ? 0 : !(!printIt)

	if (!WaveExists(PixelTurnOffs))
		String wList = reverseList(WaveListClass("PixelTurnOffs","*","DIMS:2"))
		if (ItemsInList(wList)<2)
			Wave PixelTurnOffs = $StringFromList(0,wList)
		else
			String wName
			Prompt wName,"Pixel Turn Off Wave",popup,wList
			DoPrompt "Pixel Turn Offs",wName
			if (V_flag)
				return $""
			endif
			Wave PixelTurnOffs = $wName
		endif
		printIt = 1
	endif
	if (!WaveExists(PixelTurnOffs))
		return $""
	endif
	if (printIt)
		printf "TableOfPixelTurnOffs(%s",NameOfWave(PixelTurnOffs)
		if (!ParamIsDefault(printIt))
			printf ", printIt=1"
		endif
		printf ")\r"
	endif

	String win=StringFromList(0,FindTablesWithWave(PixelTurnOffs))
	if (strlen(win))
		DoWindow/F $win
	else
		Variable left=166, top=60, height, N=DimSize(PixelTurnOffs,0)
		height = limit(16*(N+0.5) + 69, 106, 600)
		Edit/W=(left,top,left+740-8,top+height)/K=1 PixelTurnOffs.ld
		ModifyTable format(Point)=1,width(Point)=30,width(PixelTurnOffs.l)=28
	endif
	return PixelTurnOffs
End


// ================================ End of Wire Calibration Fit =================================
// ==============================================================================================
