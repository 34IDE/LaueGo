#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=1		// Use modern global access method.
#pragma ModuleName=detectorCalibration
#pragma version = 0.99
#include "microGeometryN", version>=1.83
#include "ImageDisplayScaling", version>=2.04
#include "ColorNames", version >=2.05

Static Constant hc = 1.239841984					// keV-nm

//#define OptimizeWireAxis				// if defined, optimize wire.axis[], otherwise  wire.R[]
//		for Optimization of Wire Geometry from a wire scan, see section  "=== Start of Wire Calibration Fit ==="

Menu LaueGoMainMenuName
	SubMenu "Optimize - Calibrate"
		"Print Help", PrintCalibrationListHelp()
		MenuItemIfWaveClassExists("Enter Measured Energies for Calibration","IndexedPeakList*",""),EnterMeasuredEnergies()
		MenuItemIfWaveClassExists("Set Calibration Input Data for Optimize","IndexedPeakList*",""),MakeCalibData(NaN)
		MenuItemIfWaveClassExists("Optimize All 3 Detector Geometrys","DetectorCalibrationList*",""),OptimizeAll($"",$"",$"")
		"Show Beam Line Preference, the azimuth", ShowPreferredBeamLine()
		"Set Beam Line Preference, the azimuth...",SetPreferredBeamLine("")
		"-"
		MenuItemIfWaveClassExists("  Graph Calibration Spots on 1 Detector","DetectorCalibrationList*",""),GraphCalibrationDetector($"")
		MenuItemIfWaveClassExists("  Graph Calibration Spots on All Detectors","DetectorCalibrationList*",""),GraphAllCalibrationDetector()
		MenuItemIfWaveClassExists("  Table of Calibration Data","DetectorCalibrationList*",""), DisplayTableOfWave($"",classes="DetectorCalibrationList*",promptStr="Calibration List Wave",options="DIMS:2;MAXCOLS:7;MINCOLS:7")
		MenuItemIfFunctionExists("Write Detector values to EPICS...","EPICS_put_PV_num"), WriteDetectorGeo2EPICS(NaN)
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

Function PrintCalibrationListHelp()
	print "With waves FullPeakList and FullPeakIndexed. Call MakeCalibData(NaN) to make the CalibrationList wave."
	print "	Use the \"Set Calibration Data\" button to call MakeCalibData(NaN)."
	print "	When it asks for the energy, you can either supply an Emeasured wave, or enter it later before calling OptimizeAll()."
	print " "

	print "in OptimizeAll(calib0,calib1,calib2), the sample rotation is optimized ONLY when calib0 is present."
	print " "
	print "*Required:		CalibrationList[][0,2]	= (h,k,l) the known hkl of the reflection"
	print " Optional:		CalibrationList[][3-4]	= (px,py) are measured peak position on image, if not present, then energy should be"
	print " Optional:		CalibrationList[][5]	= the measured energy of spot (keV), (some are known, some are not)"
	print "   Note:\ta line should have either (px,py), keV, or both:"
	print "		if both are present, then both are used"
	print "		if one of them is present, then that one is used"
	print "		if neither are present, then the line is ignored"
	print " "
	print "		Columns 6 & 7 are ALWAYS filled in by the optimization routine."
	print " Set Afterward:		CalibrationList[][6]	= keV calculated from sample orientation and known lattice (does not use measured spot position)"
	print " Set Afterward:		CalibrationList[][7]	= ∆E (eV),  ( CalibrationList[][5] - CalibrationList[][6] ) *1e3"
End

// ==============================================================================================
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
		Prompt H0,"H0, positioner wire H when wire splits incident beam"
		Prompt Hyc,"Hyc, positioner wire H when wire nearly over sample origin"
		Prompt F0,"F, positioner wire F for both points"
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

	Make/N=3/D/FREE origin
	Make/N=(3,3)/D/FREE rhoW
	Make/N=3/D/FREE xyzPixel, xyz0, xyz_yc
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
		PositionerX2_toBeamLineX2(g.wire,xyz0)		// beam goes through xyz0 along (001)

		origin[0] = {xyz0[0], xyz0[1], depth}				// in beam line system, just need to rotate back to positioner system
		MatrixOp/FREE origin = Inv(rhoW) x origin		// (xyz) = (wire.Rij) x (xyz),   rotate wire position from positioner coords to beam line coords
		i += 1												// count number of iterations
	while (norm(origin)>1e-9 && i<20)
	if (i>=20)
		DoAlert 0, "ERROR -- Failed to converge after "+num2str(i)+" interations"
		return "nan;nan;nan"
	endif
	origin += g.wire.origin[p]								// Final value

	if (printIt)
		printf "wire at beam	= {%.2f, %.2f, %.2f} (positioner frame)\r",X0,Y0,Z0
		printf "wire above		= {%.2f, %.2f, %.2f} (positioner frame)\r",Xyc,Yyc,Zyc
		printf "desired wire Origin	= {%.3f, %.3f, %.3f} (positioner frame)\r",origin[0],origin[1],origin[2]
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
	return str
End

// ===================================== End of Wire Origin =====================================
// ==============================================================================================



// ==============================================================================================
// ==================================== Start of Optimization ===================================

//	The user supplies ONLY columns [0,2], [4-4], [7], the rest are calculated.
//
//	CalibrationList[][0,2]	= [REQUIRED] hkl for each measured spot
//	CalibrationList[][3,4]	= [OPTIONAL] (px,py) are measured on image
//	CalibrationList[][5]	= [OPTIONAL] keV, MEASURED energy of spot (otherwise leave as NaN)
//
//	CalibrationList[][6]	= [CALCULATED] keV calculated from sample orientation and known lattice (does NOT use measured spot position)
//	CalibrationList[][7]	= [CALCULATED] ∆E (eV),  ( CalibrationList[][6] - CalibrationList[][5] ) *1e3

Function/WAVE MakeCalibData(dNum)					// Make the calibration list needed by OptimizeAll()
	Variable dNum
	dNum = mod(dNum,1) ? NaN : dNum
	dNum = dNum!=limit(dNum,0,2) ? NaN : dNum

	String wlist = reverseList(WaveListClass("FittedPeakList","*","DIMS:2"))
	String peakListName = StringFromList(0,wlist)
	Wave FullPeakList = $peakListName
	Variable getPeakList = (!WaveInClass(FullPeakList,"FittedPeakList") || ItemsInList(wlist)!=1)

	STRUCT microGeometry g
	FillGeometryStructDefault(g)										//fill the geometry structure with current values

	if (numtype(dNum) || getpeakList)
		Prompt dNum,"detector number",popup, DetectorMenuList(g)
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
			return $""
		endif
		Wave FullPeakList = $peakListName
	endif
	printf "MakeCalibData(%d)\t\t// using FullPeakList='%s'\r",dNum,peakListName
	Variable Npeaks=DimSize(FullPeakList,0), k

	// Find the correct FullPeakIndexed wave, may be a FullPeakIndexedAux for dNum = 1 or 2
	String str, indexList=""
	wlist = reverseList(WaveListClass("IndexedPeakList*","*",""))
	for (k=0;k<ItemsInList(wlist);k+=1)								// find the best choice
		str = StringByKey("peakListWave",note($StringFromList(k,wlist)),"=")
		str = StringFromList(0,str,",")
		if (stringmatch(peakListName,ParseFilePath(3,str,":",0,0)))
			indexList = StringFromList(k,wlist)
		endif
	endfor
	Wave FullPeakIndexed=$indexList
	if (!WaveExists(FullPeakIndexed))
		DoAlert 0,"Cannot find the correct FullPeakIndexed wave"
		return $""
	endif
	Variable Nindex=DimSize(FullPeakIndexed,0), m

	// always ask for EmeasuredName, even side detectors can have energies, algthough they don't need them
	String EmeasuredName=""										// find the list of energies that goes with FullPeakList
	wlist = WaveListClass("measuredEnergiesEX","*","MINROWS:1,MAXROWS:"+num2istr(Npeaks)) + SelectString(dnum,"...none yet...;","_none_;")
	if (ItemsInList(wlist))
		Prompt EmeasuredName,"Measured Energies",popup,wlist
		DoPrompt "Energies",EmeasuredName
		if (V_flag)
			return $""
		endif
	endif
	Wave Emeasured = $EmeasuredName									// Emeasured should line up with FullPeakList (NOT indexed peaks)
	printf "For detector %d,  using FullPeakIndexed='%s'",dNum,NameOfWave(FullPeakIndexed)
	if (WaveExists(Emeasured))
		printf ",   with Energies from '%s'\r",NameOfWave(Emeasured)
	else
		printf ",   with NO Energies\r"
	endif

	Variable rhox, rhoy, rhoz
	if (dNum==0)															// set sample rotation {rhox,rhoy,rhoz} from the FullPeakIndexed wave
		Wave axis = findSampleAxis(FullPeakIndexed)			// rotation vector for sample
		rhox = axis[0]
		rhoy = axis[1]
		rhoz = axis[2]
	else																		// find the sample rotation sample {rhox,rhoy,rhoz} from CalibrationList0
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
			return $""
		endif
	endif
	printf "	starting with a sample rotation of  {%g, %g, %g}\r",rhox,rhoy,rhoz

	Variable azimuth	// rotation about incident beam (usually 90° for 34-ID-E, 180° for 34-ID-C)
	STRUCT BeamLinePreference BLp										// 1st look in BL preferences
	LoadPackagePreferences "microGeo" , "microBLprefs", 0, BLp
	if (V_bytesRead>1)				// read the BLp
		azimuth = numtype(BLp.azimuth) || abs(BLp.azimuth)>180 ? NaN : BLp.azimuth
	endif
	azimuth = numtype(azimuth) || abs(azimuth)>180 ? ProbableAzimuthalAngleDetector(g.d[0]) : azimuth	// 2nd try to guess from detector0
	azimuth = numtype(azimuth) || abs(azimuth)>180 ? 90 : azimuth	// 3rd give up and just use 90°

	String name = CleanupName("CalibrationList"+ReplaceString("FullPeakList",NameOfWave(FullPeakList),""),0)
	name = name[0,29]+num2istr(dNum)
	Make/N=(Npeaks,8)/O/D $name/Wave=cList = NaN					// list of measured reflections used to calibrate the detector, {px,py,qx,qy,yz,keV,theta}
	SetDimLabel 1,0,H,cList;			SetDimLabel 1,1,K,cList;	SetDimLabel 1,2,L,cList	// (hkl)
	SetDimLabel 1,3,px,cList;			SetDimLabel 1,4,py,cList	// measued pixel
	SetDimLabel 1,5,keV_meas,cList										// measured keV
	SetDimLabel 1,6,calculated_keV,cList								// energy calculated from sample orientation (no measured parameters)
	SetDimLabel 1,7,deltaE_eV,cList										// measured - calculated energy (eV)
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
	noteStr = ReplaceNumberByKey("azimuth",noteStr,azimuth,"=")
	Note/K cList, noteStr

	Variable tolerance, i									// dpixel tolerance to find indexed spots that match fitted peaks
	for (tolerance=Inf,m=0; m<(Npeaks-1); m+=1)
		for (i=m+1;i<Npeaks;i+=1)
			tolerance = min((FullPeakList[m][0]-FullPeakList[i][0])^2 +(FullPeakList[m][1]-FullPeakList[i][1])^2, tolerance)
		endfor
	endfor
	tolerance = tolerance/4								// in units of pixel^2,  half way between the two nearest peaks

	// Find closest indexed peak to each of the fitted peaks, must also be within tolerance
	Make/N=3/D/FREE hkl
	Variable px,py, dist, kBest, distBest, keV
	for (m=0;m<Npeaks;m+=1)								// for each fitted peak, find corresponding indexed peak (if one is close enough)
		px = FullPeakList[m][0]							// fitted peak positions
		py = FullPeakList[m][1]
		for (k=0,distBest=Inf; k<Nindex; k+=1)
			dist = (px-FullPeakIndexed[k][9])^2 + (py-FullPeakIndexed[k][10])^2
			if (dist < distBest)
				distBest = dist
				kBest = k
			endif
		endfor
		if (distBest > tolerance)							// this peak not indexed, it is of no use to the calibration list
			continue
		endif

		hkl = FullPeakIndexed[kBest][p+3]				// hkl for this fitted peak
		cList[m][0,2] = hkl[q]								// CalibrationList[][0-2] = hkl are hkl for each measured spot
		cList[m][3] = px										// CalibrationList[][3,4] = (px,py) are MEASURED from image
		cList[m][4] = py
		clist[m][5] = Emeasured_from_hkl(Emeasured,hkl)	// returns NaN if no measured E (or if Emeasured does not exist)
	endfor
	// columns 6 & 7 will be filled out later at the start of optimization
	cList[][6] = NaN											// initially there is no calculated energy, or error
	cList[][7] = NaN

	Variable haveE, haveHKL, havePixel
	for (m=Npeaks-1;m>=0;m-=1)							// remove all bad lines, must be hkl and (pxiel or energy)
		haveHKL = numtype(cList[m][0]+cList[m][1]+cList[m][2]) == 0
		havePixel = numtype(cList[m][3]+cList[m][4]) == 0
		haveE = numtype(cList[m][5]) == 0
		if (!haveHKL || (!haveE && !havePixel))		// check for valid hkl, and either pixel or energy
			DeletePoints/M=0 m, 1, cList
		endif
	endfor

	if (!WaveExists(Emeasured) && dNum==0)
		sprintf str, "Detector 0\rYou will not be able to run \"Optimize\" until you have entered measured energies into \"%s\"", NameOfWave(cList)
		DoAlert 0, str
		DisplayTableOfWave(cList,classes="DetectorCalibrationList*",promptStr="Calibration List Wave",options="DIMS:2;MAXCOLS:7;MINCOLS:7")
	endif
	printf "Calibration table is:  \"%s\"\r",GetWavesDataFolder(cList,2)
	return cList
End
//
Static Function Emeasured_from_hkl(Emeasured,hkl)	// returns keV from line with hkl = Emeasured[m][1,3]
	Wave Emeasured		// wave with columns {Emeasured, h,k,l, px,py}
	Wave hkl
	if (!WaveExists(Emeasured) || !WaveExists(hkl))
		return NaN
	endif
	Make/N=3/D/FREE dhkl
	Variable i, keV, N=DimSize(Emeasured,0)
	for (i=0,keV=NaN; i<N; i+=1)
		dhkl = hkl[p] - Emeasured[i][p+1]
		if (norm(dhkl) < 0.01)
			keV = Emeasured[i][0]
			break
		endif
	endfor
	return keV
End
//
Static Function/WAVE findSampleAxis(windex,[perfect])		// find the axis of the sample from FullPeakIndexed waves
	Wave windex
	Variable perfect													// if set to true, then always return the perfect Si rotation
	perfect = ParamIsDefault(perfect) ? 0 : perfect		// default to NOT perfect

	String wName=""
	if (!WaveInClass(windex,"IndexedPeakList") && !perfect)
		String wlist = reverseList(WaveListClass("IndexedPeakList","*",""))+"perfect Si;"
		Prompt wName, "Indexed List to provide sample rotation",popup,wlist
		DoPrompt "Sample Rotation",wName
		if (V_flag)
			return $""
		endif
		Wave windex = $wName
		perfect = stringmatch(wName,"perfect Si")
	endif

	if (perfect)
		Make/N=3/D/FREE perfectSiAxis = {-2.2256407716336,-0.884037203122317,-0.38914203759791}
		Note/K perfectSiAxis,"waveClass=sampleAxis;"
		return perfectSiAxis
	endif

	// get rotation matrix from wave note of windex
	Variable r00,r10,r20, r01,r11,r21, r02,r12,r22
	sscanf StringByKey("rotation_matrix0",note(windex),"="), "{{%g,%g,%g}{%g,%g,%g}{%g,%g,%g}}", r00,r10,r20, r01,r11,r21, r02,r12,r22
	if (V_flag!=9)
		printf "Unable to get sample axis from '"+NameOfWave(windex)+"', cannot get rotation_matrix0 from wave note"
		return $""
	endif

	wName = NameOfWave(windex)+"_SampleAxis"		// name of wave to hold result
	if (strlen(wName)>30)									// make a shorter wave name
		wName = ReplaceString("FullPeakIndexed_",NameOfWave(windex),"")
		wName = ReplaceString("FullPeakIndexed",wName,"")
		wName = wName[0,19]+"_SampleAxis"
	endif

	Make/N=3/D/FREE axis
	Note/K axis, ReplaceStringByKey("waveClass",note(windex),"sampleAxis","=")

	Make/N=(3,3)/D/FREE rotMat
	rotMat[0][0]=r00;	rotMat[0][1]=r01;	rotMat[0][2]=r02
	rotMat[1][0]=r10;	rotMat[1][1]=r11;	rotMat[1][2]=r12
	rotMat[2][0]=r20;	rotMat[2][1]=r21;	rotMat[2][2]=r22
	Variable angle = axisOfMatrix(rotMat,axis,squareUp=1)	// get axis from rotation matrix
	axis *= angle*PI/180									// re-scale axis to length in radians
	return axis
End



Function/WAVE EnterMeasuredEnergies()				// Allows easy entry of Emeasured[]
	String indexLists = WaveListClass("IndexedPeakList","*","MINCOLS:2")

	if (ItemsInList(indexLists)<1)
		DoAlert 0,"No IndexedPeakList files found. You need to index an image first"
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

	if (!WaveExists(FullPeakIndexed))
		DoAlert 0,"Could not find FullPeakIndexed"
		return $""
	endif
	Variable Npeaks = DimSize(FullPeakIndexed,0)	// number of peaks
	if (Npeaks<1)
		String str
		sprintf str, "the wave \"%s\" is empty",NameOfWave(FullPeakIndexed)
		DoAlert 0, str
		return $""
	endif

	Variable dNum = FullPeakIndexed[0][11]		// detector number
	STRUCT microGeometry g
	FillGeometryStructDefault(g)						//fill the geometry structure with current values
	String detectorID = g.d[dNum].detectorID
	String color=detectorID2color(detectorID)	// detector color, {orange, yellow, purple}

	String EmeasuredName=UniqueName("Emeasured"+color,1,0)
	Make/N=(Npeaks,6) $EmeasuredName/WAVE=Em = NaN
	String wNote="waveClass=measuredEnergiesEX;"
	wNote = ReplaceStringByKey("detectorID",wNote,detectorID,"=")
	Note/K Em, wNote
	SetDimLabel 1,0,keV_measure,Em
	SetDimLabel 1,1,H,Em	;	SetDimLabel 1,2,K,Em	;	SetDimLabel 1,3,L,Em
	SetDimLabel 1,4,px,Em	;	SetDimLabel 1,5,py,Em
	Em[][1,3] = FullPeakIndexed[p][q+2]			// hkl in columns 3, 4, 5
	Em[][4,5] = round(FullPeakIndexed[p][q+5])	// px,py in columns 9, 10
	DisplayTableOfWave(Em)
	DoWindow/C EnergyInput
	DoUpdate
	DoAlert 0,"Kill the table after you have filled it in.\rYou do not need to fill in energies for every line."
	PauseForUser EnergyInput,EnergyInput
	Em = numtype(Em) ? NaN : Em						// ensure that there are no Inf's

	MatrixOP/FREE NkeV0 = sum(greater(col(Em,0),0))	// number of measured energies
	Variable NkeV = NkeV0[0]
	if (NkeV<1)
		KillWaves/Z Em
		print "There are NO measured energies, so no wave created."
		return $""
	endif

	Variable i
	for (i=Npeaks-1; i>=0; i-=1)
		if (!(Em[i][0]>0))
			DeletePoints/M=0 i, 1, Em					// remove all un-needed lines
		endif
	endfor

	if (strlen(GetRTStackInfo(2))<1 || stringmatch(GetRTStackInfo(2),"CalibrationButtonProc"))
		printf "Created  \"%s\"  containing the %g measured energies\r",NameOfWave(Em),NkeV
	endif
	return Em
End


// this uses pre-defined CalibrationList's and takes them as input, it sets up and optimizes each of the three detectors
Function OptimizeAll(calib0,calib1,calib2,[printIt])
	Wave calib0,calib1,calib2								// calibration lists
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : !(!printIt)
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

	String cName0,cName1,cName2
	if (!WaveExists(calib1) || !WaveExists(calib1) || !WaveExists(calib1))
		cName0 = SelectString(WaveExists(calib0),"CalibrationList0",NameOfWave(calib0))
		cName1 = SelectString(WaveExists(calib1),"CalibrationList1",NameOfWave(calib1))
		cName2 = SelectString(WaveExists(calib2),"CalibrationList2",NameOfWave(calib2))
		String wList0 = WaveListClass("DetectorCalibrationList0","*","DIMS:2;MAXCOLS:7;MINCOLS:7")+"_none_"
		String wList12 = WaveListClass("DetectorCalibrationList","*","DIMS:2;MAXCOLS:7;MINCOLS:7")+"_none_"
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

		printIt = 1
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
		MatrixOP/FREE Nenergies = sum(greater(col(calib0,5),0))
		print "Nenergies =",Nenergies[0]
		if (Nenergies[0]<1)				// must have some measured energies
			DoAlert 0, "The first argument "+NameOfWave(calib0)+" has NO measured energies, it must have at least 1."
			return 1
		endif
	endif
	if (WaveExists(calib1))
		if (!WaveInClass(calib1,"DetectorCalibrationList*"))
			DoAlert 0, "The second argument "+NameOfWave(calib1)+" is not of 'DetectorCalibrationList' class"
			return 1
		endif
		if (stringmatch(calib0FullName,calib1FullName))
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

	Variable Rstart, Rend, Rangle
	Variable err, seconds, failed=0
	String noteStr, errList=""
	Variable rhox=NaN, rhoy=NaN, rhoz=NaN, rhox0=NaN, rhoy0=NaN, rhoz0=NaN
	STRUCT detectorGeometry d								// this is the structure that will be optimized
	Variable err3=0											// sum of error from all 3 detectors
	Variable dNum												// detector number {0,1,2}
	for (dNum=0;dNum<MAX_Ndetectors && !failed;dNum+=1)
		if (!(g.d[dNum].used))
			continue												// skip un-used detectors
		endif
		Wave cList = $StringFromList(dNum,cNameList)	// data to use for optimization
		if (!WaveExists(cList))
			continue
		endif
		FillOutCalibTable(cList)							// fill in columns that can be calculated, this MODIFIES clist

		noteStr=note(cList)
		if (dNum && numtype(rhox+rhoy+rhoz)==0)
			noteStr = ReplaceNumberByKey("rhox",noteStr,rhox,"=")
			noteStr = ReplaceNumberByKey("rhoy",noteStr,rhoy,"=")
			noteStr = ReplaceNumberByKey("rhoz",noteStr,rhoz,"=")
			Note/K cList, noteStr
		endif
		d = g.d[dNum]											// set new detector structure to start as the old one, was:  CopyDetectorGeometry(d,g.d[dNum])
		Rstart = sqrt((d.R[0]*d.R[0])+(d.R[1]*d.R[1])+(d.R[2]*d.R[2]))*180/PI
		rhox = NumberByKey("rhox",noteStr,"=")
		rhoy = NumberByKey("rhoy",noteStr,"=")
		rhoz = NumberByKey("rhoz",noteStr,"=")
		rhox0 = rhox	;	rhoy0 = rhoy	;	rhoz0 = rhoz	// store starting rotation
		if (printIt)
			printf "Detector %d\r",dNum
			Rangle = sqrt((g.d[dNum].R[0])^2 + (g.d[dNum].R[1])^2 + (g.d[dNum].R[2])^2)*180/PI
			printf "started at  R={%g,%g,%g},   P={%g,%g,%g}mm,   |R| = %g°\r",d.R[0],d.R[1],d.R[2],(d.P[0])/1000,(d.P[1])/1000,(d.P[2])/1000,Rstart
			Variable angle = sqrt(rhox^2 + rhoy^2 + rhoz^2)*180/PI
			printf "sample started at axis = {%g, %g, %g},   |sample angle| = %g°\r",rhox,rhoy,rhoz,angle
		endif

		failed = OptimizeDetectorCalibration(d,cList)
		noteStr=note(cList)
		if (failed && printIt)
			print "----------------------- ERROR in Optimize"
			// print " ",note(cList)
			print "---",OptimizeError2str(NumberByKey("ERROR",noteStr,"="),NumberByKey("V_OptTermCode",noteStr,"="),NumberByKey("V_OptNumIters",noteStr,"="))
			printf "-----------------------"
		endif
		if ( NumberByKey("V_OptNumIters",noteStr,"=")<50 && dNum==0)
			printf "-----------------------\rPOSSIBLE ERROR (no. of iterations = %d is too few), re-run this\r  -----------------------\r", NumberByKey("V_OptNumIters",noteStr,"=")
		endif
		Rend = sqrt((d.R[0]*d.R[0])+(d.R[1]*d.R[1])+(d.R[2]*d.R[2]))*180/PI
		err = NumberByKey("err",noteStr,"=")
		seconds = NumberByKey("exectutionSec",noteStr,"=")
		err3 += err
		errList += num2str(err)+";"
		if (dNum==0)
			rhox = NumberByKey("rhox",noteStr,"=")
			rhoy = NumberByKey("rhoy",noteStr,"=")
			rhoz = NumberByKey("rhoz",noteStr,"=")
		endif
		if (printIt)
			if (seconds>4)
				printf "\tOptimization toook %s\r",Secs2Time(seconds,5,1)
			endif
			printf "ended at  R={%g,%g,%g},   P={%g,%g,%g}mm,   |R| = %g°\r",d.R[0],d.R[1],d.R[2],(d.P[0])/1000,(d.P[1])/1000,(d.P[2])/1000,Rend
			if (dNum==0)
				printf "  final sample rotation is  {%g,%g,%g},   |sample angle| = %g°",rhox,rhoy,rhoz,sqrt(rhox^2 + rhoy^2 + rhoz^2)*180/PI
				Make/N=3/D/FREE drho = {rhox-rhox0, rhoy-rhoy0, rhoz-rhoz0}
				if (numtype(sum(drho))==0)
					printf ",    %srho = %s, |%srho|=%.2g°\r",GDELTA,vec2str(drho,fmt="%.2g",sep=", "),GDELTA,norm(drho)*180/PI
				endif
				printf "\r"
			endif
			printf "error started at %g,   reduced to  %g,   after %d iterations\r",NumberByKey("errStart",noteStr,"="),err, NumberByKey("V_OptNumIters",noteStr,"=")
			printf "(final - initial) --> ∆R={%.2g,%.2g,%.2g},   ∆P={%.3f,%.3f,%.3f}\r",(d.R[0]-g.d[dNum].R[0]), (d.R[1]-g.d[dNum].R[1]), (d.R[2]-g.d[dNum].R[2]), (d.P[0]-g.d[dNum].P[0])/1000, (d.P[1]-g.d[dNum].P[1])/1000, (d.P[2]-g.d[dNum].P[2])/1000
		endif
		g.d[dNum] = d											// updagte structure with fitted values for this detector, was:  CopyDetectorGeometry(g.d[dNum],d)
	endfor
	err3 /= (g.Ndetectors)

	Variable m, maxDeltaE=0
	for (m=0;m<ItemsInList(cNameList);m+=1)
		Wave cList = $StringFromList(m,cNameList)	// data to use for optimization
		if (!WaveExists(cList))
			continue
		endif
		ImageStats/G={0,DimSize(cList,0)-1, 7,7} cList
		if (V_max>abs(maxDeltaE) || (-V_min)>abs(maxDeltaE))
			maxDeltaE = abs(V_max) > abs(V_min) ? V_max : V_min
		endif
	endfor

	if (printIt && numtype(Rend)==0)
		print " "
		DoAlert 1, "Print optimized geometry to history?"
		if (V_flag==1)
			printGeometry(g)
		endif
		printf "errList = %s     total error = %g",errList,err3
		if (abs(maxDeltaE) > 0)
			printf ",   max(|%sE|) = %.3f\r",GDELTA,maxDeltaE
		endif
		print "\r "
	endif

	if (!failed && printIt)
		String str
		sprintf str, "Update current Geometry with these fitted values\r  starting err = %.3g  -->  %s", NumberByKey("errStart",noteStr,"="), errList
		DoAlert 1, str
		if (V_flag==1)
			print "Updated the current geometry with these values"
			UpdateDefaultGeometryStruct(g,local=1)
		else
			print "REJECTED these values, no update done"
		endif
	endif
	return (failed ? NaN : err3)
End
//
Static Function FillOutCalibTable(clist)	// fill in columns that can be calculated, this MODIFIES clist
	// modifies calculated_keV=[6], deltaE_eV=[7]
	Wave clist							// table with calibration data
	String wnote=note(clist)

	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))				// fill the lattice structure with current values
		DoAlert 0, "no crystal structure found"
		return 1
	endif
	Wave recip0 = recipFrom_xtal(xtal)					// STANDARD reciprocal lattice, not rotated to match actual sample
	wnote = ReplaceStringByKey("recipSTD",wnote, encodeMatAsStr(recip0), "=")
	Note/K clist, wnote

	Variable rhox=NumberByKey("rhox",wnote,"="), rhoy=NumberByKey("rhoy",wnote,"="), rhoz=NumberByKey("rhoz",wnote,"=")
	Make/N=(3,3)/D/FREE rhoSample
	RotationVec2Matrix(rhox,rhoy,rhoz,NaN,rhoSample)	// Compute rhoSample, the sample rotation mat. Apply to qhat's to make them line up with (px,py) pairs

	Make/N=3/D/FREE hkl, ki={0,0,1}
	Variable keVmeas, keVcalc, theta, eV, m, N=DimSize(clist,0)
	for (m=0;m<N;m+=1)
		hkl = clist[m][p]
		keVmeas = cList[m][5]									// measured energy
		MatrixOp/FREE qvec = recip0 x hkl
		qvec = abs(qvec)<1e-14 ? 0 : qvec
		theta = asin(hc/keVmeas*norm(qvec)/4/PI)*180/PI // CalibrationList[][9] = theta (deg), obtained from calculated d(hkl) and measured energy using: Q = 4π sin(theta)/lambda

		keVcalc = NaN
		eV = NaN
		if (WaveExists(rhoSample))
			MatrixOP/FREE qhat = Normalize(rhoSample x qvec)
			keVcalc = -hc * norm(qvec) / (4*PI*MatrixDot(qhat,ki))
			eV = (keVmeas - keVcalc) * 1000
		endif
		cList[m][6] = keVcalc
		cList[m][7] = eV
	endfor
	return 0
End
//
Static Function OptimizeDetectorCalibration(d,CalibrationList)	// optimizes values in d to best fit values in CalibrationList
	STRUCT detectorGeometry &d						// detector geometry to optimize
	Wave CalibrationList								// data to use for optimization
	if (!WaveInClass(CalibrationList,"DetectorCalibrationList*"))
		DoAlert 0, "The CalibrationList wave is not of waveClass, 'DetectorCalibrationList*'"
		return -1
	endif

	Variable timer=startMSTimer

	String str, noteStr=note(CalibrationList)
	Variable zero = WaveInClass(CalibrationList,"DetectorCalibrationList0")	// a zero detector, fit rho too
	Variable rhox=NumberByKey("rhox",noteStr,"="), rhoy=NumberByKey("rhoy",noteStr,"="), rhoz=NumberByKey("rhoz",noteStr,"=")

	d.R[0] = (d.R[0] == 0) ? 0.02 : d.R[0]			// avoid starting exactly on a zero
	d.R[1] = (d.R[1] == 0) ? 0.02 : d.R[1]			// offset 0 angles by ~1°
	d.R[2] = (d.R[2] == 0) ? 0.02 : d.R[2]
	d.P[0] = (d.P[0] == 0) ? 500 : d.P[0]			// offset 0 positions by 0.5mm
	d.P[1] = (d.P[1] == 0) ? 500 : d.P[1]
	d.P[2] = (d.P[2] == 0) ? 500 : d.P[2]

	// Variable/G printErrX = 1
	Variable err, errStart
	if (zero)
		errStart=CalibrationErrorRPrho(CalibrationList,d.R[0],d.R[1],d.R[2],d.P[0],d.P[1],d.P[2],rhox,rhoy,rhoz)	// save error at start
	else
		errStart=CalibrationErrorRP(CalibrationList,d.R[0],d.R[1],d.R[2],d.P[0],d.P[1],d.P[2])							// save error at start
	endif

	Variable maxIters=NumVarOrDefault("maxIters",1000)
	if (zero)
		Make/N=9/D/FREE optimizeStepWave={d.R[0],d.R[1],d.R[2],d.P[0],d.P[1],d.P[2],rhox,rhoy,rhoz}
		optimizeStepWave = abs(optimizeStepWave)
		optimizeStepWave = max(optimizeStepWave[p],0.1)
		Optimize/M={0,0}/X={d.R[0],d.R[1],d.R[2],d.P[0],d.P[1],d.P[2],rhox,rhoy,rhoz}/Y=(errStart)/Q/I=(maxIters)/R=optimizeStepWave CalibrationErrorRPrho, CalibrationList
	else
		Make/N=6/D/FREE optimizeStepWave={d.R[0],d.R[1],d.R[2],d.P[0],d.P[1],d.P[2]}
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
	if (zero)												// only update rho if itis being fitted
		noteStr = ReplaceNumberByKey("rhox",noteStr,W_Extremum[6],"=")
		noteStr = ReplaceNumberByKey("rhoy",noteStr,W_Extremum[7],"=")
		noteStr = ReplaceNumberByKey("rhoz",noteStr,W_Extremum[8],"=")
	endif
	noteStr = ReplaceNumberByKey("exectutionSec",noteStr,stopMSTimer(timer)*1e-6,"=")
	Note/K CalibrationList, noteStr

	Variable dR = 1e-8									// this gives a limit of (~0.5e-6)°
	d.R[0]= round(W_Extremum[0]/dR)*dR			// round to 8 places
	d.R[1]= round(W_Extremum[1]/dR)*dR
	d.R[2]= round(W_Extremum[2]/dR)*dR
	d.R[0] = (d.R[0]==0) ? 0 : d.R[0]				// change -0 to +0
	d.R[1] = (d.R[1]==0) ? 0 : d.R[1]
	d.R[2] = (d.R[2]==0) ? 0 : d.R[2]

	Variable dP = min(d.sizeX/d.Nx,d.sizeY/d.Ny)/100	// resolution limit in Px,Py,Pz (1% of a pixel)
	dP = 10^round(log(dP))							// make rounding to a fixed number of places (not arbitrary number)
	d.P[0]= round(W_Extremum[3]/dP)*dP
	d.P[1]= round(W_Extremum[4]/dP)*dP
	d.P[2]= round(W_Extremum[5]/dP)*dP
	d.P[0] = (d.P[0]==0) ? 0 : d.P[0]				// change -0 to +0
	d.P[1] = (d.P[1]==0) ? 0 : d.P[1]
	d.P[2] = (d.P[2]==0) ? 0 : d.P[2]

	d.geoNote = "Optimized using "+NameOfWave(CalibrationList)
	sprintf str,"%s, %s (%g)",date(),Secs2Time(DateTime,3),date2secs(-1,-1,-1)/3600
	d.timeMeasured = str

	CalibrationList[][6] = CalculatedQvecEnergy(CalibrationList,p)	// fill calculated energy from sample orientation and known lattice constants
	CalibrationList[][7] = ( CalibrationList[p][5] - CalibrationList[p][6] ) *1e3	// ∆E (eV)

	KillWaves/Z W_Extremum, W_OptGradient
	return V_flag
End
//
// This routine is for optimizing the detector position & rotation and sample orientation simultaneously (only useful for detector 1,2,...   not 0)
Function CalibrationErrorRP(CalibrationList,Rx,Ry,Rz,Px,Py,Pz)	// returns error between measured and calculated spots, only called by OptimizeDetectorCalibration()
	Wave CalibrationList
	Variable Rx,Ry,Rz						// rotation vector for detector, the vectors [R] & [P] are changed to minimize returned rms value
	Variable Px,Py,Pz						// translation vector for detector
	// This routine uses the (px,py) from CalibrationList, and compares it to the angles calculated from (hkl) and lattice,
	// it also compares the measured theta (from px,py) to the theta calculated from the measured energy.
	//
	//	The user supplies ONLY columns [0-2]=hkl, [3,4]=(px,py), [5]=measured_keV, the rest are calculated.

	//	printf "CalibrationErrorRP(%s,%g,%g,%g,   %g,%g,%g)\r",NameOfWave(CalibrationList),Rx,Ry,Rz,Px,Py,Pz
	String noteStr = note(CalibrationList)
	Variable rhox=NumberByKey("rhox",noteStr,"="), rhoy=NumberByKey("rhoy",noteStr,"="), rhoz=NumberByKey("rhoz",noteStr,"=")
	Make/N=(3,3)/D/FREE rhoSample							// rotation matrix computed from {rhox,rhoy,rhoz}
	RotationVec2Matrix(rhox,rhoy,rhoz,NaN,rhoSample)	// Compute rhoSample, the sample rotation mat. Apply to qhat's to make them line up with (px,py) pairs

	Variable rms = CalibrationErrorRPinternal(CalibrationList,Rx,Ry,Rz,Px,Py,Pz,rhoSample)	// returns error between measured and calculated spots, only called by OptimizeDetectorCalibration()

	//	if (NumVarOrDefault("printErrX",0))
	//		printf "rms = %g\r",rms
	//	endif
	return rms														// error between calculated and measured (rad)
End
//
// This routine is for optimizing the detector position & rotation and sample orientation simultaneously (only useful for detector0)
Function CalibrationErrorRPrho(CalibrationList,Rx,Ry,Rz,Px,Py,Pz,rhox,rhoy,rhoz)	// returns error between measured and calculated spots, only called by OptimizeDetectorCalibration()
	Wave CalibrationList
	Variable Rx,Ry,Rz					// rotation vector for detector, these two vectors a changed to minimize returned value
	Variable Px,Py,Pz					// translation vector for detector
	Variable rhox,rhoy,rhoz			// rotation vector to apply to qhat's to make them line up with (px,py) pairs
	// This routine uses the (px,py) from CalibrationList, and compares it to the angles calculated from (hkl) and lattice,
	// it also compares the measured theta (from px,py) to the theta calculated from the measured energy.
	//
	//	The user supplies ONLY columns [0-2]=hkl, [3,4]=(px,py), [5]=measured_keV, the rest are calculated.

	//	printf "CalibrationErrorRPrho(%s,%g,%g,%g,   %g,%g,%g,   %g,%g,%g)\r",NameOfWave(CalibrationList),Rx,Ry,Rz,Px,Py,Pz,,rhox,rhoy,rhoz
	String noteStr = note(CalibrationList)
	Make/N=(3,3)/D/FREE rhoSample							// rotation matrix computed from {rhox,rhoy,rhoz}
	RotationVec2Matrix(rhox,rhoy,rhoz,NaN,rhoSample)	// Compute rhoSample, the sample rotation mat. Apply to qhat's to make them line up with (px,py) pairs

	Variable rms = CalibrationErrorRPinternal(CalibrationList,Rx,Ry,Rz,Px,Py,Pz,rhoSample)	// returns error between measured and calculated spots, only called by OptimizeDetectorCalibration()

	// find error term that sets orientation about the incident beam.
	Variable azimuth=NumberByKey("azimuth",noteStr,"=")	// rotation about incident beam (usually 90° for 34-ID-E, 180° for 34-ID-C)
	azimuth = (numtype(azimuth) || abs(azimuth)>180) ? 90 : azimuth	// should never need this

	// find error term that sets orientation about the incident beam.
	STRUCT detectorGeometry d								// a local version of detector structure
	Assemble_d(d, noteStr, Rx,Ry,Rz, Px,Py,Pz)			// fill d to match this orientation
	Variable errOrient = errAzimuthAngle(d,azimuth)	// error to use for orienting system around Z-axs (controls rotation about incident beam)
	errOrient *= 10												// weight this error heavier

	return sqrt(rms^2 + errOrient^2)
End
//
Static Function CalibrationErrorRPinternal(CalibrationList,Rx,Ry,Rz,Px,Py,Pz,rhoSample)	// returns error between measured and calculated spots, only called by OptimizeDetectorCalibration()
	Wave CalibrationList
	Variable Rx,Ry,Rz						// rotation vector for detector, these two vectors a changed to minimize returned value
	Variable Px,Py,Pz						// translation vector for detector
	Wave rhoSample							// rotation matrix to apply to qhat's to make them line up with (px,py) pairs
	// This routine uses the (px,py) from CalibrationList, and compares it to the angles calculated from (hkl) and lattice, 
	// it also compares the measured theta (from px,py) to the theta calculated from the measured energy
	// The routine does not use the {qx,qy,qz}, or the theta columns in CalibrationList
	//
	//	The user supplies ONLY columns [0-2]=hkl, [3,4]=(px,py), [5]=measured_keV, the rest are calculated.
	//
	//	CalibrationList[][0,2]= [REQUIRED] hkl for each measured spot
	//	CalibrationList[][3,4]= [OPTIONAL] (px,py) are measured on image
	//	CalibrationList[][5]	= [OPTIONAL] keV, MEASURED energy of spot (otherwise leave as NaN)
	//
	//	CalibrationList[][6]	= [CALCULATED] keV calculated from sample orientation and known lattice (does NOT use measured spot position)
	//	CalibrationList[][7]	= [CALCULATED] ∆E (eV),  ( CalibrationList[][6] - CalibrationList[][5] ) *1e3
	//
	// if CalibrationList[][5] is valid, then have energy
	// if CalibrationList[][3,4] are valid, then have pixel position
	// if CalibrationList[][3,4,5] are valid, then have pixel position & energy

	String noteStr = note(CalibrationList)

	STRUCT detectorGeometry d							// a local version of detector structure
	Assemble_d(d, noteStr, Rx,Ry,Rz, Px,Py,Pz)		// fill d to match this orientation

	Wave recip0 = decodeMatFromStr(StringByKey("recipSTD",noteStr, "="))

	Make/N=3/D/FREE ki={0,0,1}, kf, hkl, qmeas
	Make/N=3/D/FREE radialHat							// unit vector in qmeas direction
	Make/N=3/D/FREE dq
	Variable N=DimSize(CalibrationList,0)			// number of measured spots to use in computing error values
	Variable keV, Qlen, i, haveE, haveQhat
	Variable weightSum, dqRad, dqPerp, err2, pixelX, pixelY, dq2

	MatrixOP/FREE Nenergy0 = sum(greater(col(CalibrationList,5),0))	// number of measured energies
	Variable weightE = N/max(Nenergy0[0],1)			// weighting for points with measured energy
	weightE = min(weightE*weightE,40)					//   this is weight for err^2, when no energy use weight=1

	for (i=0,weightSum=0,err2=0; i<N; i+=1)
		pixelX = CalibrationList[i][3]	;	pixelX = pixelX==limit(pixelX,0,d.Nx - 1) ? pixelX : NaN
		pixelY = CalibrationList[i][4]	;	pixelY = pixelY==limit(pixelY,0,d.Ny - 1) ? pixelY : NaN
		keV = CalibrationList[i][5]
		keV = numtype(keV)==0 && keV>0 ? keV : NaN
		haveQhat = numtype(pixelX+pixelY) == 0		// valid pixelX & pixelY are present
		haveE = keV > 0										// valid energy is present

		hkl = CalibrationList[i][p]
		if (numtype(sum(hkl)))
			continue											// Invalid hkl, skip, you must have the hkl
		elseif (!haveQhat && !haveE)
			continue											// neither energy nor direction, nothing was measured
		endif

		// compute CALCULATED Q-vector from just rhoSample and the reciprocal lattice
		MatrixOp/FREE qcalc = rhoSample x recip0 x hkl	// rotate Q-vector by rhoSample, this is the Calculated qvec

		// compute MEASURED Q-vectors from (pixelX,pixelY) & keV
		if (haveQhat)										// pixel is valid
			pixel2XYZ(d,pixelX,pixelY,kf)				// find kf, convert pixel position to the beam line coordinate system
			normalize(kf)
			qmeas = kf - ki									// direction of qMeasured, wrong length, but parallel to correct q
		else
			qmeas = qcalc									// user gave an hkl & keV, but no (pixelX,pixelY), so just use calculated direction
		endif
		normalize(qmeas)
		radialHat = qmeas									// will need this later to calc the error

		// compute |Qmeas| using:  sin(theta) = -MatrixDot(ki,qmeas),  Q=4π sin(theta)/lambda,  if no energy, then use |Q| of calculated instead
		Qlen = haveE ? -4*PI*MatrixDot(ki,qmeas)*keV/hc : norm(qcalc)	// measured |Q|
		if (numtype(Qlen) || Qlen <= 0)
			continue											// inalid |Q|, give up on this spot
		endif
		qmeas *= Qlen										// set length of qmeas using Qlen from either measured keV or calculated energy

		// now have qcalc[][3] & qmeas[][3], calculate the difference appropriately
		dq = qcalc - qmeas
		dqRad = MatrixDot(radialHat,dq)				// length of dq in radial direction, the radial error in dq
		dq -= dqRad*radialHat								// remove radial part of dq, dq --> perpendicular part
		dqPerp = norm(dq)									// length of dq in perpendicular direction, the perpendicular error in dq

		dq2  = haveQhat ? (dqPerp*dqPerp) : 0		// ignore perpendicular error if no measured Q^, only have a measured E
		dq2 += haveE ? weightE*(dqRad*dqRad) : 0	// ignore radial error if no measured energy, only have a measured Q^
		err2 += dq2 / (Qlen*Qlen)						// weight by Qlen, so units of err2 are ~radian^2
		weightSum += haveE ? weightE : 1				// accumulate Sum(weight), weight points with energy more heavily
	endfor
	Variable rms = sqrt(err2 / weightSum)
	return rms
End
//
Static Function errAzimuthAngle(d,azimuth)// calculate the error to use for orienting system around the Z-axs, Y-axis is up
	// this error fixes the rotation of the whole system about the incident beam.
	STRUCT detectorGeometry, &d
	Variable azimuth			// user supplied desired rotaion angle°, sets the rotation about the z-azis, 0=X-axis, 90=Y-axis, ...

	Variable tanSize=min(d.sizeX,d.sizeY)/abs(d.P[2])// tan(angular size of detector)
	Make/N=3/D/FREE xyz
	pixel2XYZ(d,(d.Nx-1)/2,(d.Ny-1)/2,xyz)			// vector to detector center
	Variable tanCenter = sqrt(xyz[0]^2+xyz[1]^2)/abs(xyz[2])	// tan(angle to detector center)

	Variable angle											// desired rotation angle about z-axis, rotation angle from x-axis to center of detector (radian)
	if (tanSize>tanCenter)								// close to forward (or back) scattering
		// find direction along detector pixels (x or y) that is closest to y-axis, and make that the vertical
		Make/N=(3,3)/D/FREE rho
		rho[0][0] = {d.rho00, d.rho10, d.rho20}		// rotation matrix
		rho[0][1] = {d.rho01, d.rho11, d.rho21}
		rho[0][2] = {d.rho02, d.rho12, d.rho22}

		Make/N=(3,4)/D/FREE alongxy						// four vectors along the four directions ±x & ±y
		alongxy[0][0] = {1,0,0}
		alongxy[0][1] = {-1,0,0}
		alongxy[0][2] = {0,1,0}
		alongxy[0][3] = {0,-1,0}

		MatrixOP/FREE alongxy = rho x alongxy		// unit vector pointing along 4 directions of detector in beam line coords
		MatrixOP/FREE dotYs = row(alongxy,1)			// dot of each alongxy with y^
		WaveStats/Q/M=1 dotYs
//		angle = atan2(alongxy[0][V_maxloc],alongxy[1][V_maxloc])	// angle between "most vertical direction" and yhat
		angle = atan2(alongxy[1][V_maxloc],alongxy[0][V_maxloc])	// angle between xhat and "most vertical direction" and yhat

	else														// far from forward scattering (the usual)
		Variable minus = sign(d.P[2])
		angle = atan2(minus * d.rho12, minus * d.rho02)
	endif

	Variable errOrient = angle - (azimuth*PI/180)
	// printf "errOrient = %g (only angle in radians)\r",errOrient
	return errOrient										// error to force correct orientation (radian)
End
//	//
//	Static Function errVertAngle(d,angleOffset)// calculate the error to use for orienting system around the Z-axs, Y-axis is up
//		// this error fixes the rotation of the whole system about the incident beam.
//		STRUCT detectorGeometry, &d
//		Variable angleOffset		// user supplied angle offset°, when 0 detector face perp to y-axis, rotates about z-axis
//	
//		Variable tanSize=min(d.sizeX,d.sizeY)/abs(d.P[2])// tan(angular size of detector)
//		Make/N=3/D/FREE xyz
//		pixel2XYZ(d,(d.Nx-1)/2,(d.Ny-1)/2,xyz)			// vector to detector center
//		Variable tanCenter = sqrt(xyz[0]^2+xyz[1]^2)/abs(xyz[2])	// tan(angle to detector center)
//	
//		Variable angle											// angle between disired vertical and y^ (radian)
//		if (tanSize>tanCenter)								// close to forward (or back) scattering
//			// find direction along detector pixels (x or y) that is closest to y-axis, and make that the vertical
//			Make/N=(3,3)/D/FREE rho
//			rho[0][0] = {d.rho00, d.rho10, d.rho20}		// rotation matrix
//			rho[0][1] = {d.rho01, d.rho11, d.rho21}
//			rho[0][2] = {d.rho02, d.rho12, d.rho22}
//	
//			Make/N=(3,4)/D/FREE alongxy						// four vectors along the four directions ±x & ±y
//			alongxy[0][0] = {1,0,0}
//			alongxy[0][1] = {-1,0,0}
//			alongxy[0][2] = {0,1,0}
//			alongxy[0][3] = {0,-1,0}
//	
//			MatrixOP/FREE alongxy = rho x alongxy		// unit vector pointing along 4 directions of detector in beam line coords
//			MatrixOP/FREE dotYs = row(alongxy,1)			// dot of each alongxy with y^
//			WaveStats/Q/M=1 dotYs
//			angle = atan2(alongxy[0][V_maxloc],alongxy[1][V_maxloc])	// angle between "most vertical direction" and yhat
//	
//		else														// far from forward scattering (the usual)
//			angle = atan2(d.rho02,d.rho12)					// angle from y-axis to detector normal
//		endif
//	
//		Variable errOrient = angle - (angleOffset*PI/180)
//		// printf "errOrient = %g (only angle in radians)\r",errOrient
//		return errOrient										// error to force correct orientation (radian)
//	End
//
Static Function Assemble_d(d, noteStr, Rx,Ry,Rz, Px,Py,Pz)	// fill d to match this orientation
	STRUCT detectorGeometry &d							// a local version of detector structure, this is filled
	String noteStr											// noteStr = note(CalibrationList)
	Variable  Rx,Ry,Rz, Px,Py,Pz
	d.used = 1
	d.Nx = NumberByKey("Nx",noteStr,"=")				// number of un-binned pixels in whole detector
	d.Ny = NumberByKey("Ny",noteStr,"=")
	d.sizeX = NumberByKey("sizeX",noteStr,"=")		// outside size of detector (micron)
	d.sizeY = NumberByKey("sizeY",noteStr,"=")		// outside size of detector (micron)
	d.R[0] = Rx;	d.R[1] = Ry;	d.R[2] = Rz				// Rotation vector for detector (degree)
	d.P[0] = Px;	d.P[1] = Py;	d.P[2] = Pz				// offset to detector (micron)
	DetectorUpdateCalc(d)									// update all fields in this detector structure (basically rho)
End


Static Function CalculatedQvecEnergy(CalibrationList,m)	// fills the calculated qcalc for point m in CalibrationList, returns calculated energy
	Wave CalibrationList
	Variable m

	Make/N=3/D/FREE ki={0,0,1}
	String wnote = note(CalibrationList)

	Variable rhox=NumberByKey("rhox",wnote,"="), rhoy=NumberByKey("rhoy",wnote,"="), rhoz=NumberByKey("rhoz",wnote,"=")
	Make/N=(3,3)/D/FREE rhoSample							// rotation matrix computed from {rhox,rhoy,rhoz}
	RotationVec2Matrix(rhox,rhoy,rhoz,NaN,rhoSample)	// Compute the sample rotation rhoSample. Apply to qhat's to make them line up with (px,py) pairs

	// compute ideal Q-vector from just rhoSample and the reciprocal lattice
	Wave recip0 = decodeMatFromStr(StringByKey("recipSTD",wnote, "="))	// returns a FREE wave defined by str
	Make/N=3/D/FREE hkl = CalibrationList[m][p]
	MatrixOP/FREE qcalc = rhoSample x recip0 x hkl	// rotate Q-vector by rhoSample
	Variable Qlen = norm(qcalc)

	Variable sintheta = -MatrixDot(ki,qcalc)/Qlen		// sin(theta) = -ki dot q^ 
	return hc * Qlen / (4*PI*sintheta)					// energy (keV)
End




Function GraphAllCalibrationDetector()
	String list = WaveListClass("DetectorCalibrationList*","*","DIMS:2;MAXCOLS:7;MINCOLS:7")
	Variable i,N=ItemsInList(list)
	for (i=0;i<N;i+=1)
		GraphCalibrationDetector($StringFromList(i,list))
	endfor
End

Function GraphCalibrationDetector(calib)
	Wave calib
	String calibName,wList = WaveListClass("DetectorCalibrationList*","*","DIMS:2;MAXCOLS:7;MINCOLS:7")
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

	Display /W=(xoff,40,xoff+wid,40+wid)/K=1/N=$wName calib[*][4] vs calib[*][3]
	ModifyGraph width={Aspect,sizeX/sizeY}
	AppendToGraph calib[*][4] vs calib[*][3]
	ModifyGraph gfMult=110, tick=2, mirror=1, minor=1, lowTrip=0.001, standoff=0, mode=3
	ModifyGraph rgb($calibName)=(0,0,0), marker($calibName)=8, marker($cname1)=19

	ImageStats/M=1/G={0,DimSize(calib,0)-1,5,5} calib
	Variable NkeV=V_npnts
	ImageStats/M=1/G={0,DimSize(calib,0)-1,7,7} calib
	Variable NdE=V_npnts
	if (NdE>2 && NdE>NkeV/2)
		Variable zrange = max(abs(V_max),abs(V_min))
		ModifyGraph zColor($cname1)={calib[*][7],-zrange,zrange,RedWhiteBlue256}
		ColorScale/C/N=textColorScale "∆E  (eV)"
	else
		ModifyGraph zColor($cname1)={calib[*][5],*,*,Rainbow256}
		ColorScale/C/N=textColorScale "E  (keV)"
	endif
	ColorScale/C/N=textColorScale/F=0/S=3/M/B=1/A=RC/X=2.94/Y=0.98/E trace=$cname1
	ColorScale/C/N=textColorScale lblMargin=0, minor=1

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
	Variable rx,ry,rz								// rotation axis
	Variable theta									// rotation angle (rad) (if NaN, then assume,   | {rx,ry,rz} | = theta (rad)
	Wave rho											// a (3,3) matrix to recieve rotation matrix

	if (DimSize(rho,0)!=3 && DimSize(rho,1)!=3)
		Abort "in rotationVec2Matrix(), mat must be (3,3)"
	endif

	Variable len = sqrt(rx*rx + ry*ry + rz*rz)	// length of {rx,ry,rz}
	theta = numtype(theta) ? len : theta			// if theta invalid, assume len is theta (rad)
	if (theta==0 || len==0)							// zero angle rotation just return the identity matrix
		rho = (p==q)
		return 0
	elseif (numtype(rx+ry+rz))						// totally invalid
		rho = NaN
		return 1
	endif
	Variable c=cos(theta), s=sin(theta)
	Variable c1 = 1-c
	rx /= len;	ry /= len;	rz /= len				// make |{rx,ry,rz}| = 1

	// this is the Rodrigues formula from:   http://mathworld.wolfram.com/RodriguesRotationFormula.html
	rho[0][0] = c + rx*rx*c1;			rho[0][1] = rx*ry*c1 - rz*s;		rho[0][2] = ry*s + rx*rz*c1
	rho[1][0] = rz*s + rx*ry*c1;		rho[1][1] = c + ry*ry*c1;			rho[1][2] = -rx*s + ry*rz*c1
	rho[2][0] = -ry*s + rx*rz*c1;	rho[2][1] = rx*s + ry*rz*c1;		rho[2][2] = c + rz*rz*c1
End

// ==================================== End of Optimization =====================================
// ==============================================================================================



// ==============================================================================================
// ==================================== Start of EPICS update ===================================

Function WriteDetectorGeo2EPICS(dNum)
	Variable dNum

	STRUCT microGeometry g
	if (FillGeometryStructDefault(g, alert=1))		//fill the geometry structure with current values
		return 1
	endif

	if (Exists("EPICS_put_PV_num")!=6)
		DoAlert 0,"EPICS not available"
		return 1
	elseif (dNum!=0 && dNum!=1 && dNum!=2)			// dNum must be 0, 1, or 2
		Prompt dNum,"detector number",popup, DetectorMenuList(g)
		DoPrompt "Detector Number",dNum
		if (V_flag)
			return 1
		endif
		dNum -= 1
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



Function findPxPyFromPixel(detNum,px,py)	// For Detector in the beam, find P from (px,py)
	Variable detNum									// detector number, usually 0
	Variable px,py										// pixel at center
	STRUCT microGeometry g
	FillGeometryStructDefault(g)					//fill the geometry structure with current values
	STRUCT detectorGeometry d
	Variable detNumMax = (g.Ndetectors)-1

	if (!(detNum>=0 && px>=0 && py>=0))		// need to prompt
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
		printf "•findPxPyFromPixel(%g, %.6g, %.6g)\r",detNum,px,py
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


Static Function detectorNumFromQvec(g,qvec,[depth])	// given Qvec, returns detector number {0,1,2,...} or -1 for an error
	STRUCT microGeometry &g								// contains all the detector info
	Wave qvec													// q-vector, returns id num of detector struck by this qvec
	Variable depth												// sample depth measured along the beam
	depth = ParamIsDefault(depth) ? NaN : depth	// default is NaN, XYZ2pixel() will handle this properly

	Variable i, id=-1, dist=Inf
	for (i=0;i<MAX_Ndetectors;i+=1)						// check each detector, find the closest one that detects the x-ray
		if (! g.d[i].used)
			continue												// skip un-used detectors
		endif
		Wave xyz = q2XYZ(g.d[i],qvec,depth=depth)	// get xyz position where diffraction from qvec hits detector
		if (!WaveExists(xyz))
			continue												// diffracted ray failed to hit detector
		endif
		if (norm(xyz)<dist)
			dist = norm(xyz)									// closer than previous, save it
			id = i
		endif
	endfor
	return id
End
//
Static Function/WAVE q2XYZ(d,qvec,[depth])			// returns xyz, position where diffraction from qvec hits detector
	STRUCT detectorGeometry &d
	Wave qvec													// qvec need not be normalized
	Variable depth												// sample depth measured along the beam
	depth = ParamIsDefault(depth) ? NaN : depth	// default is NaN, XYZ2pixel() will handle this properly

	Make/N=3/D/FREE kout, qhat=qvec, ki={0,0,1}	// ki = geo.ki[p],  incident beam direction
	normalize(qhat)
	//	normalize(ki)

	Variable qLen = -2*MatrixDot(qhat,ki)				// length of qhat, note (q^ dot -ki) always positive
	if (qLen<0)													// this occurs for theta<0, (we do not want to reflect from the back side)
		return $""
	endif
	kout = qhat*qLen + ki									// kf - ki = q

	Variable px,py												// final pixel position, full chip unbinned 0 based pixels
	XYZ2pixel(d,kout,px,py,depth=depth)
	if (px>=0 &&  px < d.Nx && py>=0 && py < d.Ny)
		Make/N=3/D/FREE xyz=NaN							// pixel is on detector
		pixel2XYZ(d,px,py,xyz)
	else
		Wave xyz = $""											// pixel noes not lie on detector
	endif
	return xyz
End
//	Function test_detectorNumFromQvec()
//		STRUCT microGeometry g
//		FillGeometryStructDefault(g)							//fill the geometry structure with current values
//	
//		Make/N=(3,3)/D/FREE recip0 = {{11.5,0,0}, {0,11.5,0}, {0,0,11.5}}
//		Make/N=(3,3)/D/FREE rhoSample
//		Make/N=3/D/FREE hkl={0, 0, 4}, qcalc
//		Variable rhox=3.58149181126276, rhoy=-1.34908451768581, rhoz=-0.611186383166512;
//		detectorCalibration#RotationVec2Matrix(rhox,rhoy,rhoz,NaN,rhoSample)	// Compute rhoSample, the sample rotation mat. Apply to qhat's to make them line up with (px,py) pairs
//		Variable id
//	
//		print "these should come from Orange"
//		MatrixOp/FREE qcalc = rhoSample x recip0 x hkl	// rotate Q- vector by rhoSample
//		id = detectorCalibration#detectorNumFromQvec(g,qcalc)
//		print "	",id,g.d[id].color
//	
//		hkl={3, -1, 13}
//		MatrixOp/FREE qcalc = rhoSample x recip0 x hkl	// rotate Q- vector by rhoSample
//		id = detectorCalibration#detectorNumFromQvec(g,qcalc)
//		print "	",id,g.d[id].color
//	
//		print "this should come from Yellow"
//		hkl={7, -1, 10}
//		MatrixOp/FREE qcalc = rhoSample x recip0 x hkl	// rotate Q- vector by rhoSample
//		id = detectorCalibration#detectorNumFromQvec(g,qcalc)
//		print "	",id,g.d[id].color
//	
//		print "these should come from Purple"
//		hkl={0, 8, 12}
//		MatrixOp/FREE qcalc = rhoSample x recip0 x hkl	// rotate Q- vector by rhoSample
//		id = detectorCalibration#detectorNumFromQvec(g,qcalc)
//		print "	",id,g.d[id].color
//	
//		hkl={0, 6, 14}
//		MatrixOp/FREE qcalc = rhoSample x recip0 x hkl	// rotate Q- vector by rhoSample
//		id = detectorCalibration#detectorNumFromQvec(g,qcalc)
//		print "	",id,g.d[id].color
//	
//		hkl={-2, 6, 12}
//		MatrixOp/FREE qcalc = rhoSample x recip0 x hkl	// rotate Q- vector by rhoSample
//		id = detectorCalibration#detectorNumFromQvec(g,qcalc)
//		print "	",id,g.d[id].color
//	End


Static Function/T OptimizeError2str(flag,OptTermCode,OptNumIters)
	Variable flag
	Variable OptTermCode
	Variable OptNumIters
	if (!flag)
		return ""
	endif

	String flagErrs = "57:User abort;"
	flagErrs += "788:Iteration limit was exceeded.;"
	flagErrs += "789:	Maximum step size was exceeded in five consecutive iterations.;"
	flagErrs += "790:	The number of points in the typical X size wave specified by /R does not match the number of X values specified by the /X flag;"
	flagErrs += "791:	Gradient nearly zero and no iterations taken. This means the starting point is very nearly a critical point. It could be a solution, or it could be so close to a saddle point or a maximum (when searching for a minimum) that the gradient has no useful information. Try a slightly different starting point.;"

	String termErrs = "1:Gradient tolerance was satisfied.;"
	termErrs += "2:Step size tolerance was satisfied.;"
	termErrs += "3:No step was found that was better than the last iteration. This could be because the current step is a solution, or your function may be too nonlinear for Optimize to solve, or your tolerances may be too large (or too small), or finite difference gradients are not sufficiently accurate for this problem.;"
	termErrs += "4:Iteration limit was exceeded.;"
	termErrs += "5:Maximum step size was exceeded in five consecutive iterations. This may mean that the maximum step size is too small, or that the function is unbounded in the search direction (that is, goes to -inf if you are searching for a minimum), or that the function approaches the solution asymptotically (function is bounded but doesn't have a well-defined extreme point).;"
	termErrs += "6:Same as V_flag = 791.;"

	String out = StringByKey(num2istr(flag), flagErrs)
	String str = StringByKey(num2istr(OptTermCode), termErrs)
	if (flag==788 && OptTermCode==4)
		str = ""
	endif

	if (strlen(str))
		if (strlen(out))
			out += "\r   " + str
		else
			out = str
		endif
	endif
	if (strlen(out))
		sprintf str, "   %d iterations", OptNumIters
		out += str
	endif
	return out
End

// ===================================== End of Calibration =====================================
// ==============================================================================================
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
	oneErr =  OptimizeAll(c0,c1,c2, printIt=1)
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
				g.d[i].R[0] += dR/2;		g.d[i].R[1] += dR;		g.d[i].R[2] += 1.5*dR	// and add imprecision so we do not start at the answer, approximately 10° & 10mm
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
		manyErrs[i] = OptimizeAll(c0,c1,c2, printIt=1)
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
//		g.d[i].R[0] += dR/2;	g.d[i].R[1] += dR;		g.d[i].R[2] += 1.5*dR		// and add imprecision so we do not start at the answer, approximately 10° & 10mm
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
		printf "∆pixel has sigma=%g, and ∆eV=±%geV:\r",noise,noise/2
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
	if (FillCrystalStructDefault(xtal))				//fill the lattice structure with test values
		DoAlert 0, "no crystal structure found"
		return ""
	endif
	Wave recip0 = recipFrom_xtal(xtal)					// STANDARD reciprocal lattice, not rotated to match actual sample
	Make/N=(3,3)/D/FREE recip, rhoSample

	// sample rotation
	Make/N=3/D/FREE axisSampleFake
	axisSampleFake={-135*PI/180,0,0}
	Variable angle = norm(axisSampleFake)*180/PI	// sample rotation angle (deg)
	RotationVec2Matrix(axisSampleFake[0],axisSampleFake[1],axisSampleFake[2],NaN,rhoSample)
	MatrixOp/FREE recip = rhoSample x recip0		// rotate recip0

	STRUCT microGeometry g
	GeoReferenceOrientation(g)							// use reference geometry
//	FillGeometryStructDefault(g)							//fill the geometry structure with current values
	if (printing)
		print "for detector",dNum
		printDetector(g.d[dNum])
		printf "sample rotation axis = {%g,%g,%g},  a rotation angle of  %g°\r",axisSampleFake[0],axisSampleFake[1],axisSampleFake[2],angle
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
	for (i=N;i<25;i+=1)										// remove extra spots from pkList
		pkList = RemoveListItem(N,pkList)
	endfor

	//	CalibrationList[][0,2]	= [REQUIRED] hkl for each measured spot
	//	CalibrationList[][3,4]	= [OPTIONAL] (px,py) are measured on image
	//	CalibrationList[][5]	= [OPTIONAL] keV, MEASURED energy of spot (otherwise leave as NaN)
	//
	//	CalibrationList[][6]	= [CALCULATED] keV calculated from sample orientation and known lattice (does NOT use measured spot position)
	//	CalibrationList[][7]	= [CALCULATED] ∆E (eV),  ( CalibrationList[][6] - CalibrationList[][5] ) *1e3

	String name="CalibrationListTest"+num2istr(dNum)
	Make/N=(N,8)/O/D $name/Wave=cList =NaN	// a list of measured reflections used to calibrate the detector, {px,py,qx,qy,yz,keV,theta}
	SetDimLabel 1,0,H,cList;			SetDimLabel 1,1,K,cList;	SetDimLabel 1,2,L,cList	// (hkl)
	SetDimLabel 1,3,px,cList;			SetDimLabel 1,4,py,cList	// measued pixel
	SetDimLabel 1,5,keV_meas,cList										// measured keV
	SetDimLabel 1,6,calculated_keV,cList								// energy calculated from sample orientation (no measured parameters)
	SetDimLabel 1,7,deltaE_eV,cList										// measured - calculated energy (eV)
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
	noteStr = ReplaceStringByKey("recipSTD",noteStr, encodeMatAsStr(recip0), "=")
	Note/K cList, noteStr

	Make/N=3/D/FREE ki={0,0,1}, kf, qvec
	Variable dot, Q, keV
	Variable px,py, h,k,l, hmax, m
	for (i=0,m=0; i<N; i+=1)
		sscanf StringFromLIst(i,pkList),"%g,%g",px,py// first find integral hkl closest to (px,py)
		pixel2XYZ(g.d[dNum],px,py,kf)					// get kf from pixel on detector
		normalize(kf)
		qvec = kf - ki
		MatrixOp/FREE hkl = Inv(recip) x qvec
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
		MatrixOp/FREE qhat = recip x hkl
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
		// MatrixOp/FREE qvec =  recip0 x hkl

		cList[m][0] = h;		cList[m][1] = k;		cList[m][2] = l	// (hkl)
		cList[m][3] = px;		cList[m][4] = py								// measured pixel position
		cList[m][5] = keV															// keV,  from Q = 4π sin(theta)/lambda == Q = 4π sin(theta) * (E/hc)
		cLIst[m][6] = NaN
		cList[m][7] = NaN
		m += 1
	endfor
	N = m
	if (printing)
		print "N=",N
	endif
	Redimension/N=(N,8) cList

	if (N>1)
		ImageStats/G={0,N-1,7,7} cList				// pretend that did not measure spot with highest energy
		cList[V_maxRowLoc][8] = NaN
		cList[V_maxRowLoc][9] = NaN
	endif
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

// ================================= End of Optimization Tests ==================================
// ==============================================================================================



// ==============================================================================================
// =============================== Start of Calibration set panel ===============================

Function/T FillCalibrationParametersPanel(strStruct,hostWin,left,top)
	String strStruct									// optional passed value of xtal structure, this is used if passed
	String hostWin										// name of home window
	Variable left, top								// offsets from the left and top

	SetWindow kwTopWin,userdata(CalibrationPanelName)=hostWin+"#CalibrationPanel"
	NewPanel/K=1/W=(left,top,left+221,top+445+30)/HOST=$hostWin
	ModifyPanel frameStyle=0, frameInset=0
	RenameWindow #,CalibrationPanel

	Button buttonEnterMeasuredEnergies,pos={29,5},size={180,20},proc=detectorCalibration#CalibrationButtonProc,title="Enter Measured Energies"
	Button buttonEnterMeasuredEnergies,help={"Set the Calibration Data"}

	Button buttonMakeCalibData,pos={29,5+35},size={160,20},proc=detectorCalibration#CalibrationButtonProc,title="Set Calibration Data"
	Button buttonMakeCalibData,help={"Set the Calibration Data"}

	Button buttonCalibOptimizeAll3,pos={29,40+35},size={160,20},proc=detectorCalibration#CalibrationButtonProc,title="Optimize..."
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

	Button buttonCalibPrintHelp,pos={29,240+35},size={160,20},proc=detectorCalibration#CalibrationButtonProc,title="Print Help to History"
	Button buttonCalibPrintHelp,help={"Print some Help text to the History"}

	EnableDisableCalibControls(hostWin+"#CalibrationPanel")
	return "#CalibrationPanel"
End
//
Static Function EnableDisableCalibControls(win)	// here to enable/disable
	String win													// window (or sub window) to use
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

	Button buttonCalibPrintHelp,win=$win,disable=0	// these two are always enabled
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
		EnterMeasuredEnergies()
	elseif (stringmatch(ctrlName,"buttonMakeCalibData"))
		MakeCalibData(NaN)
	elseif (stringmatch(ctrlName,"buttonCalibOptimizeAll3"))
		OptimizeAll($"",$"",$"", printIt=1)
	elseif (stringmatch(ctrlName,"buttonCalibGraph1"))
		GraphCalibrationDetector($"")
	elseif (stringmatch(ctrlName,"buttonCalibGraphAll"))
		GraphAllCalibrationDetector()
	elseif (stringmatch(ctrlName,"buttonCalibTable"))
		DisplayTableOfWave($"",classes="DetectorCalibrationList*,measuredEnergies*",promptStr="Calibration List Wave",options="DIMS:2")
	elseif (stringmatch(ctrlName,"buttonCalibWrite2EPICS"))
		WriteDetectorGeo2EPICS(NaN)
	elseif (stringmatch(ctrlName,"buttonCalibPrintHelp"))
		print " "
		PrintCalibrationListHelp()
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

	STRUCT microGeometry geo
	if (FillGeometryStructDefault(geo))			//fill the geometry structure with test values
		return $""
	endif
	Variable detNum=detectorNumFromID(geo,StringByKey("detectorID",wnote,"="))

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
	printf "Before Optimize:  <depth> = %g,   range = [%g, %g],  ∆=%g\t\txWave = %s\r",V_avg,V_min,V_max,errStart,vec2str(xWave)
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
	printf "After Optimize:  <depth> = %g,   range = [%g, %g],  ∆=%g\t\txWave = %s\r",V_avg,V_min,V_max,errEnd,vec2str(xWave)
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

	String win=StringFromList(0,WindowsWithWave(PixelTurnOffs,2))
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



// ==============================================================================================
// ================================ Start of Preferred Beam Line ================================

Function/T ShowPreferredBeamLine()
	// Show current values of Preferred Beam Line struct
	STRUCT BeamLinePreference BLp
	LoadPackagePreferences "microGeo" , "microBLprefs", 0, BLp
	String str=""
	if (V_bytesRead>1)				// read the BLp
		print "Current values of struct:"
		str = desc_BeamLinePreference(BLp)
	else
		str = "NO values were found for 'microBLprefs', use  SetPreferredBeamLine(\"\")  to set some"
	endif
	print str
	DoAlert 0, str
	return str
End
//
Function/T SetPreferredBeamLine(BL0)
	String BL0

	String knownBLs="APS/34-ID-E;APS/34-ID-C", BL=BL0
	if (WhichListItem(BL, knownBLs)<0)
		// try to guess correct Beam Line, based on curret geometry
		STRUCT microGeometry g
		if (!FillGeometryStructDefault(g,alert=0))	//fill geometry structure with current values
			BL = SelectString( CmpStr(g.d[0].detectorID, "PE1621 723-3335")==0, BL, "APS/34-ID-E")
			BL = SelectString( g.d[0].Nx == 514, BL, "APS/34-ID-C")
		endif
		Prompt BL,"Preferred Beam Line", popup, knownBLs
		DoPrompt "Choose Preferred Beam Line", BL
		if (V_flag)
			return ""
		endif
	endif
	if (CmpStr(BL0,BL))
		printf "%sSetPreferredBeamLine(\"%s\")\r",BULLET,BL
	endif

	Variable azimuth=ProbableAzimuthalAngleDetector(g.d[0])	// start at current detector0
	String BLname=BL, id="none", color="black"
	Make/N=3/U/W/FREE rgb=0
	if (CmpStr(BL, "APS/34-ID-E")==0)
		azimuth = 90
		id = "PE1621 723-3335"
		color = SelectString(strlen(g.d[0].color), "Orange", g.d[0].color)
		rgb = {65535,43688,32768}
	elseif (CmpStr(BL, "APS/34-ID-C")==0)
		azimuth = numtype(azimuth) ? 180 : azimuth
	else
		azimuth = numtype(azimuth) ? 90 : azimuth	// know nothing, just set to a top detector
	endif

	Prompt BLname, "Name of Beam Line", popup, knownBLs+";none"
	Prompt azimuth, "Azimuthal angle° around beam for detector 0, 90°→top,  180°→-X"
	Prompt id, "ID of Detector 0"
	Prompt color, "Color for Detector 0"
	DoPrompt "Prefered Beam Line", BLname,id,color,azimuth
	if (V_flag)
		return ""
	endif
	if (CmpStr(color,"Orange"))
		Wave rgb = color2RGB(color,65535)
	endif

	STRUCT BeamLinePreference BLp
	LoadPackagePreferences "microGeo" , "microBLprefs", 0, BLp
	if (V_bytesRead>1)				// read the BLp
		print "OLD values of struct:"
		print desc_BeamLinePreference(BLp)
		print " "
	endif

	// reset BLp to new values and save
	BLp.detectorID = id
	BLp.color = color
	BLp.rgb[0] = rgb[0]
	BLp.rgb[1] = rgb[1]
	BLp.rgb[2] = rgb[2]
	BLp.BL = BLname
	BLp.azimuth = azimuth
	SavePackagePreferences/FLSH=1 "microGeo","microBLprefs",0,BLp
	print desc_BeamLinePreference(BLp)
	return ""
End
//
Static Structure BeamLinePreference		// structure definition for a prefered Beam Line (mainly for azimuth of detector0)
	uchar BL[50]						// beam line name, something like "APS/34-ID-E", or "34-ID-C", etc.
	uchar detectorID[100]			// default detector ID, for this beam line
	uchar color[30]					// name of color, e.g. "orange" "yellow", ...
	uint16 rgb[3]					// default for beamline, rgb of color, based on 65535 = full, color is optional, it will be assigned if it is not specified.
	double azimuth					// detector azimuth angle° (-180,+180], 34-ID-E is +90, 34-ID-C is +180
EndStructure
//
Static Function/T desc_BeamLinePreference(BLp)
	// returns a string describing a BeamLinePreference structure
	STRUCT BeamLinePreference &BLp
	String out="", line
	sprintf line, "Current Beam Line \"%s\"\r    Default Prefered Detector:",BLp.BL
	out = line + "\r"
	sprintf line, "      ID = \"%s\"  ",BLp.detectorID
	out += line + "\r"
	sprintf line, "      %s,  rgb={%d, %d, %d}", BLp.color, BLp.rgb[0],BLp.rgb[1],BLp.rgb[2]
	out += line + "\r"
	sprintf line, "      with a detector azimuth angle of %g°", BLp.azimuth
	out += line
	return out
End


Static Function ProbableAzimuthalAngleDetector(d)
	STRUCT detectorGeometry &d
	// returns the azimuthal angle of the detector (orientation about z-axis)
	// if existing azimuthal angle is within 5° of a multiple of 45°, then it rounds to nearest multiple of 45°.
	// for 34-ID-E, this should give 90°, for 34-ID-C 180°
	Variable angleTol = 5				// if azimuth is within 5° of multiple of 45, lockin to nearest multiple of 45
	Variable minus = sign(d.P[2])
	Variable azimuth=atan2(minus * d.rho12, minus * d.rho02) * 180/PI
	Variable absAngle45=mod(abs(azimuth),45)
	Variable deltaAzimuth = min(absAngle45,45-absAngle45)	// distance from nearst multiple of 45°
	if (deltaAzimuth <= angleTol)
		azimuth = round(azimuth/45) * 45							// round to nearest multiple of 45°
		azimuth = abs(azimuth+180)<0.1 ? 180 : azimuth		// close to -180, prefer +180
	endif
	azimuth = numtype(azimuth) || abs(azimuth)>180 ? NaN : azimuth	// a bad value
	return azimuth
End

// ================================= End of Preferred Beam Line =================================
// ==============================================================================================



