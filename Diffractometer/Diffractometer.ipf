#pragma rtGlobals=1		// Use modern global access method.
#pragma version = 0.50
#pragma ModuleName=diffractometer
#include "LatticeSym", version>=3.76
#initFunctionName "Init_Diffractometer()"

#define HAVE_VECTOR_OPS
Constant MAX_REF_ORIENTATIONS = 20				// maximum number of orientations that may be defined (on first 2 used so far)
StrConstant defaultGeoTagName="geoD"				// default start of name for geometry files
Static Constant MAX_Naxes = 10
Static strConstant DEFAULT_DIFF_NAME = "fourc"	// default diffractometer name
Static strConstant DEFAULT_AXIS_NAMES = "2-theta;theta;chi;phi"
Static Constant DEFAULT_Naxes = 4					// number of axes used by diffractometer (you can use OverRide with this!)
Static Constant hc_keVnm = 1.239841856			// h*c (keV-nm)
Static strConstant DetectorFileFilters =  "Detector Files (*.xml):.xml;All Files:.*;"
Static Constant MAX_DETECTORS = 5					// maximum number of detectors that can be stored, each detector must have a UNIQUE detectorID
Static strConstant DEF_PILATUS_ID = "PILATUS 100K, 1-0005, APS"
Static StrConstant GeoWebServer = "sector34.xray.aps.anl.gov/34ide"


// Routines with names ending with "VEC" are vectorized and will run MUCH faster
// For vectorized operations, remake image into a 1D object, where rows first, then columns (this is how Redimension works)


Menu "Diffractometer"
	SubMenu "Sample\Detector"
		"Print Sample", PrintCurrentSample()
		"Print Detectors",PrintCurrentDetectors()
		"-"
		"Load New Detector Geometry...",LoadDetectorGeoGeneral(NaN)
//		"Load Detector Geometry from a File...",LoadDetectorsFromFile("","")
		"Write Detector Geometry to File...",WriteDetectorToFile("","")
		"Enter Detector Parameters Manually...", diffractometer#SetDetectorParameters(NaN,NaN,NaN,NaN,"","","",0)
		"Re-Set Detector for a  PIlatus 100K...", diffractometer#ResetSetPilatus100Kcalibration("",NaN,NaN,NaN)
		"-"
		"Set Sample Orientation...", diffractometer#SetSampleReferenceReflections("","","","","","","",NaN,NaN)
		"Change Diffractometer Type...",diffractometer#SelectDiffractometer("")
		SubMenu "Set to Reference Values"
			"Set Sample to Reference",SetDefaultSample2Reference()
			"Set Detector to Reference",SetDefaultDetector2Reference(NaN)
		End
	End
	diffractometer#MenuIfFunctionExists("Movie of spec Images ...","specInfo"), MakeMovieSpecScanImages(NaN)
	diffractometer#MenuIfFunctionExists("Movie of Energy Scanned Images ...","GetImageFileParameters"), MakeMovieEnergyImages("","")
End

Static Function/T MenuIfFunctionExists(str,funcName)
	String str
	String funcName
	if (exists(funcName)<=2)
		return "("+str
	endif
	return str
End

// Constant lambda = 0.154		// wavelength (nm), used with default orientations
// primary		hkl = {0,1,0},  reflection:	tth=60, th=30, phi=0, kappa=0, phi=0
// secondary	hkl = {-1,0,0},  reflection:	tth=60, th=30, chi=-90, phi=0

//	MatrixOP/FREE  qout = T x UB x recip x hkl






//ThreadSafe 
Function/WAVE HKLofPixel(s,d,keV,px,py,A,[depth])// returns hkl of pixel in beam CRYSTAL frame, (not sample frame, not BL coords)
	STRUCT sampleStructure &s						// structure defining a sample
	STRUCT detectorGeometry, &d					// structure defining a detector
	Variable keV									// for keV<=0, just return qvec with length 1
	Variable px,py									// pixel on detector
	Wave A											// wave with diffractometer angles
	Variable depth									// depth of point (measured along Z from origin)
	if (!(keV>0))
		return $""
	endif

	Make/N=3/D/FREE qvec
	if (ParamIsDefault(depth))
		QofPixel(d,keV,px,py,A,qvec)				// returns Q of pixel in beam SAMPLE frame, (not crystal frame, not BL coords)
	else
		QofPixel(d,keV,px,py,A,qvec,depth=depth)	// same as above, but with depth passed
	endif
	Wave hkl = sample2crystal(s,qvec)				// rotate qvec from sample-location frame into crystal based frame, the hkl
	WaveClear qvec
	return hkl
End

//ThreadSafe 
Function QofPixel(d,keV,px,py,A,qvecIN,[depth])	// returns Q of pixel in beam SAMPLE frame, (not crystal frame, not BL coords)
	STRUCT detectorGeometry, &d
	Variable keV									// for keV<=0, just return qvec with length 1
	Variable px,py									// pixel on detector
	Wave A											// wave with diffractometer angles
	Wave qvecIN									// 3-vector to recieve the result
	Variable depth									// depth of point (measured along Z from origin)
	if (!(keV>0))
		return NaN
	endif

	Make/N=3/D/FREE qBL,qvec
	Variable theta
	if (ParamIsDefault(depth))
		theta = QofPixelBL(d,keV,px,py,A,qBL)	// returns Q of pixel in beam line frame, This is NOT q in sample
	else
		theta = QofPixelBL(d,keV,px,py,A,qBL,depth=depth)	// returns Q of pixel in beam line frame, This is NOT q in sample
	endif
	beamline2sample(A,qBL,qvec)	// rotate vector from beamline to sample-location frame, (un-does diffractometer)
	if (WaveExists(qvecIN))
		qvecIN = qvec
	endif
	return theta
End
//
// Vectorized version of QofPixel()
Threadsafe Static Function/WAVE QofPixelVEC(d,keV,pxpy,A,[depth])// returns Qs of pixels in SAMPLE frame, (not crystal frame, not BL coords)
	STRUCT detectorGeometry, &d					// one detector geometry for whole image
	Variable keV									// for keV<=0, just return qvec with length 1, one energy for whole image
	Wave pxpy										// array of pixels to use (must be in raw un-binned pixels), perhaps an ROI
	Wave A											// wave with diffractometer angles, one set of angles for whole image
	Variable depth									// depth of point (measured along Z from origin)
	if (!(keV>0))
		return $""
	endif

	if (ParamIsDefault(depth))
		Wave vBLs = QofPixelBLVEC(d,keV,pxpy,A)	// returns Q of pixel in beam line frame, This is NOT q in sample
	else
		Wave vBLs = QofPixelBLVEC(d,keV,pxpy,A,depth=depth)	// returns Q of pixel in beam line frame, This is NOT q in sample
	endif
	Wave Qsamples = beamline2sampleVEC(A,vBLs)	// rotate vector from beamline to sample-location frame, (un-does diffractometer)
	WaveClear vBLs
	return Qsamples
End


ThreadSafe Static Function beamline2sample(A,vBL,vSample)		// rotate vector from beamline to sample-location frame, (un-does diffractometer)
	Wave A
	Wave vBL						// vector in beam line frame
	Wave vSample					// vector in sample frame
	String DiffractometerName = StrVarOrDefault("root:Packages:Diffractometer:DiffractometerName",DEFAULT_DIFF_NAME)
//	String fname = "diffractometer#"+DiffractometerName+"Matrix"
	String fname = DiffractometerName+"Matrix"
	FUNCREF protoMatrix diffractMatrix = $fname
	Wave Tmat = diffractMatrix(A)
	MatrixOP/FREE vec = Inv(Tmat) x vBL			// temporary vector in sample frame
	vSample = vec
	WaveClear Tmat, vec
	return norm(vec)
End
//
// Vectorized version of beamline2sample()
ThreadSafe Static Function/WAVE beamline2sampleVEC(A,vBLs)	// rotate vector from beamline to sample-location frame, (un-does diffractometer)
	Wave A
	Wave vBLs										// vectors in beam line frame

	String DiffractometerName = StrVarOrDefault("root:Packages:Diffractometer:DiffractometerName",DEFAULT_DIFF_NAME)
//	String fname = "diffractometer#"+DiffractometerName+"Matrix"
	String fname = DiffractometerName+"Matrix"
	FUNCREF protoMatrix diffractMatrix = $fname
	Wave Tmat = diffractMatrix(A)

//	MatrixOP/FREE vSamples = Inv(Tmat) x vBLs	// vectors in sample frame
	MatrixOP/FREE vSamples = (Inv(Tmat) x vBLs^t)^t	// vectors in sample frame
	WaveClear Tmat
	return vSamples
End


//ThreadSafe 
Static Function/WAVE sample2crystal(s,vSample)	// rotate vector from sample-location frame into crystal based frame
	STRUCT sampleStructure &s						// structure definition for a sample
	Wave vSample
	Make/N=(3,3)/D/FREE UB
	UB = s.UB[p+3*q]
	MatrixOP/FREE vCrystal = Inv(UB) x vSample
	WaveClear UB
	if (numtype(sum(vCrystal)))
		WaveClear vCrystal							// causes return to pass null wave reference
	endif
	return vCrystal
End


Static Function/WAVE Measured2UB(ref0,ref1)	// compute UB from 2 measured reflections, sets the orientation of crystal w.r.t sample-location
	STRUCT oneOrientation &ref0, &ref1		// measured orientations used to find UB

	String DiffractometerName = StrVarOrDefault("root:Packages:Diffractometer:DiffractometerName",DEFAULT_DIFF_NAME)
//	String fname = "diffractometer#"+DiffractometerName+"Matrix"
	String fname = DiffractometerName+"Matrix"
	FUNCREF protoMatrix diffractMatrix = $fname
//	String dname = "diffractometer#"+DiffractometerName+"DetectorMatrix"
	String dname = DiffractometerName+"DetectorMatrix"
	FUNCREF protoMatrix detectorMatrix = $dname
	Variable Naxes=NumVarOrDefault("root:Packages:Diffractometer:Naxes",DEFAULT_Naxes)

	STRUCT detectorGeometry d
	if (FillDetectorStructDefault(d,""))			//fill the detector structure with current values
		DoAlert 0, "No Detector Structure, please set one"
		return $""
	endif
	Make/N=(Naxes)/D/FREE A					// angles of a measured point

	if (sampleRefBad(ref0) || sampleRefBad(ref1))
		return $""
	endif

	// direction of primary reflection, 2theta = A0[0]
	A = ref0.A[p]
	Make/N=3/D/FREE ki, kf, qhkl
	ki={0,0,2*PI/(ref0.lambda)}
	pixel2xyz(d,ref0.px,ref0.py,A,kf)			// convert pixel position to the beam line coordinate system
//		Wave rot = detectorMatrix(A)			// a point detector
//		MatrixOP/FREE/O kf = rot x ki
	normalize(kf)								// change length of kf to 2¹/lambda
	kf *= 2*PI/(ref0.lambda)
	Wave Tmat = diffractMatrix(A)
	MatrixOP/FREE/O qhklMeasured = Inv(Tmat) x (kf-ki)
	normalize(qhklMeasured)
	hkl2Q(ref0.h,ref0.k,ref0.l, qhkl,normal=1)
	Wave axis = wCross(qhkl,qhklMeasured)
	Variable angle = asin(normalize(axis))
	Wave UB = rotationMatFromAxis(axis,angle)// first rotation in making UB, aligns first measured H0
	WaveClear Tmat

	axis = qhklMeasured						// axis of rotation to use with second reflection

	// direction of second reflection, 2theta = A1[0]
	A = ref1.A[p]
	ki={0,0,2*PI/(ref1.lambda)}
	pixel2xyz(d,ref1.px,ref1.py,A,kf)		// convert pixel position to the beam line coordinate system
//		Wave rot = detectorMatrix(A)			// a point detector
//		MatrixOP/FREE/O kf = rot x ki
	normalize(kf)								// change length of kf to 2¹/lambda
	kf *= 2*PI/(ref1.lambda)
	Wave Tmat = diffractMatrix(A)
	MatrixOP/FREE/O qhklMeasured = Inv(Tmat) x (kf-ki)
	normalize(qhklMeasured)
	hkl2Q(ref1.h,ref1.k,ref1.l, qhkl,normal=1)// note that both qhklMeasured and qhkl are normalized
	MatrixOP/FREE/O qhkl = UB x qhkl			// put in first rotation
	// need to rotate about axis (previous qMeasured) to bring qhkl as close as possible to qhklMeasured

	// change qhkl & qhklMeasured to components that are perpendicular to axis
	MatrixOP/FREE/O qhkl = qhkl - axis x (axis . qhkl)	// subtract off componenet of qhkl that is || to axis
	MatrixOP/FREE/O qhklMeasured = qhklMeasured - axis x (axis . qhklMeasured)

	angle = asin(norm(wCross(qhkl,qhklMeasured)))	// rotate from qhkl to qhklMeasured about axis
	Wave rot2 = rotationMatFromAxis(axis,angle)		// set mat to be a rotation matrix about axis with angle

	MatrixOP/FREE/O UB = rot2 x UB
	UB = abs(UB) < 1e-14 ? 0 : UB
	WaveClear Tmat, axis, rot2, qhklMeasured, qhkl, ki,kf, A
	return UB
End




//  ============================================================================  //
//  ============================== Start of Make a Movie ==============================  //

Function MakeMovieSpecScanImages(scanNum)
	Variable scanNum

	if (!(scanNum>0))
		Prompt scanNum,"Scan Number"
		DoPrompt "Scan Number",scanNum
		if (V_flag)
			return 1
		endif
		scanNum = round(scanNum)
		printf "MakeMovieSpecScanImages(%g)\r",scanNum
	endif
	if (!(scanNum>0))
		return 1
	endif

	String motors = StrVarOrDefault("root:Packages:Diffractometer:DiffractometerAxisNames","")
	Variable Nm=ItemsInList(motors)
	Make/N=(Nm)/D/FREE A

	Variable m
	for (m=0;m<Nm;m+=1)
		A[m] = specInfo(scanNum,StringFromList(m,motors))
	endfor
	Variable keV = specInfo(scanNum,"DCM_energy")
	keV = numtype(keV) ? 10*hc_keVnm/specInfo(scanNum,"wavelength") : keV

	Wave image=$LoadSpecImage(scanNum,0,extras="quiet:1")
	if (!WaveExists(image))
		return 1
	endif
	String wnote = note(image)
	Duplicate/O image, movieFrame
	KillWaves/Z image

	Variable Np=specInfo(scanNum,"Npoints")			// number of points in this scan
	wnote = ReplaceNumberByKey("scanNum",wnote, scanNum ,"=")
	wnote = ReplaceNumberByKey("SCAN_LEN",wnote, Np ,"=")
	Note/K movieFrame, wnote

	print wnote

	String gName=MakeGraphMovieFrame()
	if (strlen(gName)<1)
		return 1
	endif

	Variable maxVal=-Inf, minVal=Inf
	Variable i, theta
	NewMovie
	for (i=0;i<Np;i+=1)
		Wave image=$LoadSpecImage(scanNum,i,extras="quiet:1")
		if (!WaveExists(image))
			DoAlert 0,"Could not load image "+num2istr(i)+ " from scan "+num2istr(scanNum)
			break
		endif
		wnote = note(image)
		movieFrame = image
		Note/K movieFrame, wnote
		maxVal = max(maxVal,WaveMax(image))
		minVal = min(minVal,WaveMin(image))
		KillWaves/Z image
		LabelReLabelGraphMovieFrame(movieFrame)
		DoWindow/F $gName
		DoUpdate
		AddMovieFrame
	endfor
	CloseMovie
	if (i<Np)
		DoAlert 0,"Could not open all of the images"
		return 1		
	endif
	printf "range of images is [%g, %g]\r",		maxVal, minVal
	return 0
End

Static Function/S MakeGraphMovieFrame()
	Wave movieFrame
	if (!WaveExists(movieFrame))
		return ""
	endif

	String gList=FindGraphsWithWave(movieFrame)
	String gName=StringFromList(0,gList)
	if (ItemsInList(gList)>1)
		String item
		Variable i
		for (i=0;i<ItemsInList(gList);i+=1)
			item = StringFromList(i,gList)
			if (strsearch(item,"GraphMovieFrame",0)==0)
				gName = item
				break
			endif
		endfor
	endif
	if (strlen(gName)>0)
		DoWindow/F $gName
		return gName
	endif

	Variable scanNum=NumberByKey("scanNum",note(movieFrame),"=")
	if (numtype(scanNum))
		return ""
	endif

	gName = "GraphMovieFrame"+num2istr(scanNum)
	Variable imageMax = 1e6
	Display /W=(459,182,1389,577)/K=1
	DoWindow /C $gName
	AppendImage movieFrame
	ModifyImage movieFrame ctab= {1,imageMax,Terrain256,1}
	ModifyImage movieFrame log= 1
	ModifyGraph height={Aspect,0.400411}, mirror=1, minor=1
	SetAxis/A/R left

	LabelReLabelGraphMovieFrame(movieFrame)

	STRUCT detectorGeometry d
	if (FillDetectorStructDefault(d,""))				//fill the detector structure with current values
		DoAlert 0, "No Detector Structure, please set one"
		return gName
	endif
	Variable Nx=(d.Nx), Ny=(d.Ny), px=(d.px0), py=(d.py0)
	if (numtype(px+py)==0)						// pixel positionof incident beam for marker
		Variable dMarker=round(max(Nx,Ny)/10)
		dMarker = min(dMarker,min(Nx,Ny))
		DrawMarker(px,py,dMarker,dMarker,"cross gap",dash=2,layer="UserAxes")
	endif

	return gName
End
//
Static Function LabelReLabelGraphMovieFrame(movieFrame)
	Wave movieFrame
	if (!WaveExists(movieFrame))
		return 1
	endif

	String wnote=note(movieFrame)
	Variable scanNum=NumberByKey("scanNum",note(movieFrame),"=")
	String timeWritten = specInfoT(scanNum,"timeWritten")
	String specCommand = specInfoT(scanNum,"specCommand")
	String specComment = StrVarOrDefault("specComment","")

	String file_name = ParseFilePath(0, StringByKey("file_name",wnote,"="), ":", 1,0)
	Variable keV=specInfo(scanNum,"keV")
	if (numtype(keV))
		keV = NumberByKey("keV",wnote,"=")
	endif

	String motors = StrVarOrDefault("root:Packages:Diffractometer:DiffractometerAxisNames","")
	Variable Nm=ItemsInList(motors)
	String mot, str= ""
	Variable i, val
	for (i=0;i<Nm;i+=1)
		mot = StringFromList(i,motors)
		val = specInfo(scanNum,mot)
		if (numtype(val)==0)
			str += "\r"+WordToGreek(mot)+" = "+num2str(val)+"¡"
		endif
	endfor
	if (keV>0)
		str += "\r"+num2str(keV)+" keV"
	endif
	if (strlen(str))
		str = "\\JR"+str
		TextBox/C/N=textPosition/F=0/B=1/X=1.68/Y=1.80 str
	endif

	str = specComment
	str += SelectString(strlen(str) && strlen(file_name),"","\r")
	str += file_name
	str += SelectString(strlen(str) && strlen(specCommand),"","\r")
	str += specCommand
	str += SelectString(strlen(str) && strlen(timeWritten),"","\r")
	str += timeWritten
	if (strlen(str))
		TextBox/C/N=textTitle/F=0/B=1/A=LT/X=1.9/Y=1.2 str
	endif
	return 0
End

//  =============================== End of Make a Movie ==============================  //
//  ============================================================================  //




//  ============================================================================  //
//  ============================ Start of Sample Orientation ============================  //

Structure sampleStructure			// structure definition for a sample
	uchar name[256]				// sample name
	STRUCT crystalStructure xtal	// structure definition for a crystal lattice
	double UB[9]					// a linear version of the UB matrix
									//   for UB[3][3],    UB = s.UB[p+3*q]
	int16 Nrefs						// number of refs that were set (max of MAX_REF_ORIENTATIONS)
	STRUCT oneOrientation refs[MAX_REF_ORIENTATIONS]	// measured orientations used to find UB
EndStructure
//
Structure oneOrientation			// structure definition, one crystal reflection used to make an orientation
	double	h,k,l
	double	A[MAX_Naxes]			// angles where the reflection was found (tth,theta,chi,phi)
	double	lambda					// wavelength used to find this reflection
	double	px,py					// pixel location of peak center (not used if px <= 0, e.g. for a Bicron detector)
EndStructure


// when you have ref0 & ref1, and want to set the Sample structure
Static Function SetSampleStruct(name,xtal,ref0,ref1,s)
	String name										// name of sample (truncated to 256)
	STRUCT crystalStructure &xtal
	STRUCT oneOrientation &ref0, &ref1			// measured orientations used to find UB
	STRUCT sampleStructure &s						// structure definition for a sample

	s.name = name[0,255]
	copy_xtal(s.xtal,xtal)
	s.Nrefs = 2
	copyOneOrientation(s.refs[0],ref0)				// copy a oneOrientation structure
	copyOneOrientation(s.refs[1],ref1)
	UpdateDefaultSampleStruct(s)
End


Function FillSampleStructDefault(s)				//fill the sample structure with current default values
	STRUCT sampleStructure &s						// returns 0 if something set, 0 is nothing done

	String strStruct=StrVarOrDefault(":sampleStructStr","")	// set to values in current directory
	if (strlen(strStruct)<1)
		strStruct=StrVarOrDefault("root:Packages:Diffractometer:sampleStructStr","")	// try the default location
	endif
	if (strlen(strStruct)>1)						// found structure information, load into s
		StructGet/S/B=2 s, strStruct
	else												// not found in this experiment, try the prefs
		LoadPackagePreferences/MIS=1 "Diffractometer","samplePrefs",0,s
		if (V_flag)
			return 1								// did nothing, nothing found
		endif
	endif
	UpdateDefaultSampleStruct(s)					// this ensures that if we loaded it from prefs, that local values are updated
	return 0
End


Function UpdateDefaultSampleStruct(s,[local])		// Update the default location with values in s
	STRUCT sampleStructure &s						// returns 0 if something set, 0 is nothing done
	Variable local									// if true will also update sampleStructStr in the current folder if it already exists (It will NOT create local sampleStructStr)
	local = ParamIsDefault(local) ? 0 : local
	local = numtype(local) ? 0 : !(!local)

	SampleUpdateCalc(s)
	String sStructStr
	StructPut/S/B=2 s, sStructStr

	String/G root:Packages:Diffractometer:sampleStructStr=sStructStr		// always save to default location
	SavePackagePreferences/FLSH=1 "Diffractometer","samplePrefs",0,s	// alwasy update prefs too

	local = local || (exists(":sampleStructStr")==2)// save locally if told to, or if a local already exists in current datafolder
	if (local)
		String/G :sampleStructStr=sStructStr		// make a local copy too
	endif
	return 0
End


Function SetDefaultSample2Reference()				// set default sample to the reference values
	STRUCT sampleStructure s
	if (SampleReferenceOrientation(s))				// set geometry to reference values
		print "ERROR -- in SampleReferenceOrientation(), nothing changed"
		return 1
	endif
	UpdateDefaultSampleStruct(s)
	printSampleStructure(s)
	return 0
End
//
Static Function SampleReferenceOrientation(s)		// sets s to the reference orientation (sort of an ideal set of values)
	STRUCT sampleStructure &s

	String DiffractometerName = StrVarOrDefault("root:Packages:Diffractometer:DiffractometerName",DEFAULT_DIFF_NAME), str=""
	STRUCT detectorGeometry d
	Variable px=0, py=0
	if (!FillDetectorStructDefault(d,""))			//fill the detector structure with test values
		px = d.px0									// center of detector
		py = d.py0
	endif

	String name = DiffractometerName+"ReferenceOrientation"
	FUNCREF protoReferenceOrientation referenceOrientation = $name
	if (referenceOrientation(s,px,py))
		sprintf str,"ERROR -- the Diffractometer '%s' is an unknown type.  Cannot set a Default Reference Orientation", DiffractometerName
		DoAlert 0,str
		print str
		return 1
	endif
	return 0
End


Static Function SampleUpdateCalc(s)
	STRUCT sampleStructure &s						// structure definition for a sample

	Wave UB = Measured2UB(s.refs[0],s.refs[1])
	if (WaveExists(UB))							// if a valid UB matrix, then reset structue, otherwise structure unchanged
		Make/N=9/FREE/D line=UB
		Variable i
		for (i=0;i<9;i+=1)
			s.UB[i] = line[i]
		endfor
	endif
End


Function copyOneOrientation(dest,in)				// copy a oneOrientation structure
	STRUCT oneOrientation &dest, &in				// dest is the destination, in is source
	Variable Naxes=NumVarOrDefault("root:Packages:Diffractometer:Naxes",DEFAULT_Naxes)
	dest.h = in.h
	dest.k = in.k
	dest.l = in.l
	Variable i
	for (i=0;i<Naxes;i+=1)
		dest.A[i] = in.A[i]
	endfor
	dest.lambda = in.lambda
	dest.px = in.px
	dest.py = in.py
End

Function PrintCurrentSample()						// prints the current default sample to the history
	STRUCT sampleStructure s
	FillSampleStructDefault(s)						//fill the geometry structure with current values
	printSampleStructure(s)
End
//
Function printSampleStructure(s)
	STRUCT sampleStructure &s						// structure definition for a sample

	Variable Naxes=NumVarOrDefault("root:Packages:Diffractometer:Naxes",DEFAULT_Naxes)
	if (strlen(s.name))
		printf "sample:  '%s'\r", s.name
	endif

	if (!LatticeSym#LatticeBad(s.xtal))
		String sym = getHMboth(s.xtal.SpaceGroup)
		printf "    for '%s'  lattice is  #%d   %s     %.9gnm, %.9gnm, %.9gnm,   %g¡, %g¡, %g¡\r",s.xtal.desc,s.xtal.SpaceGroup,sym,s.xtal.a,s.xtal.b,s.xtal.c,s.xtal.alpha,s.xtal.beta,s.xtal.gam
	endif

	if (s.Nrefs > 0 && !sampleRefBad(s.refs[0]) && !sampleRefBad(s.refs[1]))
		printf " the %d reference reflections are:\r",s.Nrefs
		Variable i,m, N=min(s.Nrefs, MAX_REF_ORIENTATIONS)
		if (N > s.Nrefs)
			printf "***  ERROR -- s.Nrefs = %g,   which is greater than MAX_REF_ORIENTATIONS (=%g)  ***\r",s.Nrefs, MAX_REF_ORIENTATIONS
		endif
		for (i=0;i< N ; i+=1)
			printf "    hkl[%d]=(%g %g %g),  angles = {%g",i,s.refs[i].h, s.refs[i].k, s.refs[i].l, s.refs[i].A[0]
			for (m=1;m<Naxes;m+=1)
				printf ", %g",s.refs[i].A[m]
			endfor
			printf "},   lambda = %g nm\r", s.refs[i].lambda
		endfor
	else
		print "No reference reflections found, UB probably taken directly from spec"
	endif
	printf "UB =\r"
	printf "  %+5.3f  %+5.3f  %+5.3f\r",s.UB[0],s.UB[3],s.UB[6]
	printf "  %+5.3f  %+5.3f  %+5.3f\r",s.UB[1],s.UB[4],s.UB[7]
	printf "  %+5.3f  %+5.3f  %+5.3f\r",s.UB[2],s.UB[5],s.UB[8]
End



Static Function SetSampleReferenceReflections(name,hklStr0,Astr0,hklStr1,Astr1,pstr0,pstr1,lam0,lam1)
	String name
	String hklStr0, Astr0
	String hklStr1, Astr1
	String pstr0,pstr1
	Variable lam0, lam1

	STRUCT sampleStructure s						// structure definition for a sample
	FillSampleStructDefault(s)

	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))				//fill the lattice structure with current values
		DoAlert 0, "No Lattice, please set one"
		return 1
	endif

	Variable Naxes=NumVarOrDefault("root:Packages:Diffractometer:Naxes",DEFAULT_Naxes)
	String str, fmt="%.13g", fmtSpace="%.13g"
	Variable i
	for (i=1;i<Naxes;i+=1)
		fmt += ", %.13g"
		fmtSpace += " %.13g"
	endfor

	name = SelectString(strlen(name),s.name,name)
	Variable ask=numtype(lam0+lam1)	>0			// flag to prompt the user
	lam0 = lam0>0 ? lam0 : s.refs[0].lambda
	lam1 = lam1>1 ? lam1 : s.refs[1].lambda

	Wave vec = string2wave(hklStr0)
	if (numpnts(vec)!=3)
		ask = 1
		sprintf hklStr0, "%g, %g, %g", s.refs[0].h, s.refs[0].k, s.refs[0].l
	endif
	Wave vec = string2wave(hklStr1)
	if (numpnts(vec)!=3)
		ask = 1
		sprintf hklStr1, "%g, %g, %g", s.refs[1].h, s.refs[1].k, s.refs[1].l
	endif
	Wave vec = string2wave(Astr0)
	if (numpnts(vec)!=Naxes)
		ask = 1
		Make/N=(Naxes)/D/FREE Avec
		Avec = s.refs[0].A[p]
		Astr0 = vec2str(Avec,bare=1,places=13)
	endif
	Wave vec = string2wave(Astr1)
	if (numpnts(vec)!=Naxes)
		ask = 1
		Make/N=(Naxes)/D/FREE Avec
		Avec = s.refs[1].A[p]
		Astr1 = vec2str(Avec,bare=1,places=13)
	endif

	Wave vec = string2wave(pstr0)
	if (numpnts(vec)!=2)
		ask = 1
		sprintf pstr0, "%g, %g", s.refs[0].px, s.refs[0].py
	endif

	Wave vec = string2wave(pstr1)
	if (numpnts(vec)!=2)
		ask = 1
		sprintf pstr1, "%g, %g", s.refs[1].px, s.refs[1].py
	endif

	if (ask)
		Prompt hklStr0,"hkl 0"
		Prompt hklStr1,"hkl 1"
		Prompt Astr0,"Angles 0 (degree)"
		Prompt Astr1,"Angles 1 (degree)"
		Prompt pstr0,"Center Pixel 0"
		Prompt pstr1,"Center Pixel 1"
		Prompt lam0, "Wavelength 0 (nm)"
		Prompt lam1, "Wavelength 1 (nm)"
		Prompt name,"Name of Sample"
		DoPrompt "Sample Orientation",hklStr0,hklStr1,Astr0,Astr1,pstr0,pstr1,lam0,lam1,name
		if (V_flag)
			return 1
		endif
		printf "SetSampleReferenceReflections(\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%g,%g)\r",name,hklStr0,Astr0,hklStr1,Astr1,pstr0,pstr1,lam0,lam1
	endif

	STRUCT sampleStructure sf						// structure definition for a sample
	sf.Nrefs = 2
	sf.name = name
	sf.refs[0].lambda = lam0
	sf.refs[1].lambda = lam1
	if (!LatticeSym#LatticeBad(xtal))
		copy_xtal(sf.xtal,s.xtal)				
	else
		copy_xtal(sf.xtal,xtal)
	endif

	Wave vec = string2wave(hklStr0)
	if (numpnts(vec)!=3)
		return 1
	endif
	sf.refs[0].h = vec[0] ;	sf.refs[0].k = vec[1] ;	sf.refs[0].l = vec[2]

	Wave vec = string2wave(hklStr1)
	if (numpnts(vec)!=3)
		return 1
	endif
	sf.refs[1].h = vec[0] ;	sf.refs[1].k = vec[1] ;	sf.refs[1].l = vec[2]

	Wave vec = string2wave(Astr0)
	if (numpnts(vec)!=Naxes)
		return 1
	endif
	for (i=0;i<Naxes;i+=1)
		sf.refs[0].A[i] = vec[i]
	endfor

	Wave vec = string2wave(Astr1)
	if (numpnts(vec)!=Naxes)
		return 1
	endif
	for (i=0;i<Naxes;i+=1)
		sf.refs[1].A[i] = vec[i]
	endfor

	Wave vec = string2wave(pstr0)
	if (numpnts(vec)!=2)
		return 1
	endif
	sf.refs[0].px = vec[0] ;	sf.refs[0].py = vec[1]

	Wave vec = string2wave(pstr1)
	if (numpnts(vec)!=2)
		return 1
	endif
	sf.refs[1].px = vec[0] ;	sf.refs[1].py = vec[1]

	SampleUpdateCalc(sf)

	print "changed sample from:"
	printSampleStructure(s)
	print " "
	print "to:"
	printSampleStructure(sf)

	UpdateDefaultSampleStruct(s)					// update default locations with this value of s
	return 0
End
//
//Static Function LatticeBad(xtal)
//	STRUCT crystalStructure &xtal
//
//	Variable bad=0
//	bad += numtype(xtal.a + xtal.b + xtal.c)
//	bad += !(xtal.a > 0) || !(xtal.b > 0) || !(xtal.c > 0)
//	bad += numtype(xtal.alpha + xtal.beta + xtal.gam)
//	bad += !(xtal.alpha > 0) || !(xtal.beta > 0) || !(xtal.gam > 0)
//	bad += xtal.alpha >=180 || xtal.beta >= 180 || xtal.gam >= 180
//	bad += numtype(xtal.SpaceGroup) || xtal.SpaceGroup < 1 || xtal.SpaceGroup > 230
//	return (bad>0)
//End
//
Static Function sampleRefBad(ref)
	STRUCT oneOrientation &ref					// measured orientations used to find UB

	Variable N=NumVarOrDefault("root:Packages:Diffractometer:Naxes",NaN)
	if (numtype(ref.h + ref.k + ref.l + ref.lambda + ref.px + ref.py + N))
		return 1
	elseif (ref.px<0 || ref.py<0 || ref.lambda<=0)
		return 1
	elseif (N<=0)
		return 1
	endif	

	Make/N=(N)/FREE/D Asum=ref.A[p]
	if (numtype(sum(Asum)))
		return 1
	endif

	return 0									// all looks OK
End

//Function test(str)
//	String str
//
//	Wave vec = string2wave(str)
//	printWave(vec,name="vec",brief=1)
//End
Static Function/WAVE string2wave(list)
	String list
	list = ReplaceString(",",list," ")		// change all separators to spaces
	list = ReplaceString(";",list," ")
	list = ReplaceString(":",list," ")
	list = TrimTrailingWhiteSpace(list)
	do
		list = ReplaceString("  ",list," ")	// change all multiple spaces to single space
	while (strsearch(list,"  ",0)>=0)
	list = ReplaceString(" ",list,";")		// change space separators to semi-colons

	Variable i, val, N=ItemsInList(list)
	if (N<1)
		return $""
	endif
	Make/N=(N)/D/FREE vec
	for (i=0;i<N;i+=1)
		val = str2num(StringFromList(i,list))
		if (numtype(val)==2)
			return $""
		endif
		vec[i] = val
	endfor
	return vec
End


#if strlen(FunctionList("specInfo","","WIN:spec.ipf"))
// special extra stuff for spec
Function SetLocalSampleStructFromSpec(scanNum)	// sets local sample structure from spec info
	Variable scanNum

	Wave UB = specProc#UBmatrix(scanNum)		// qvec = UB x hkl
	if (numtype(sum(UB)))
		return 1
	endif

	String DiffractometerName = StrVarOrDefault("root:Packages:Diffractometer:DiffractometerName",DEFAULT_DIFF_NAME)
	String name = "diffractometer#"+DiffractometerName+"Spec2BL"
	FUNCREF protoSpec2BL genericSpec2BL = $name
	Wave spec2BL = genericSpec2BL(0)
	MatrixOP/FREE UBnew = spec2BL x UB * 10	// spec uses Angstrom, I use nm
	STRUCT sampleStructure s						// structure definition for a sample
	s.name = ""
	Variable i
	for (i=0;i<9;i+=1)
		s.UB[i] = UBnew[i]
	endfor

	Variable H0=specInfo(scanNum,"H0"), K0=specInfo(scanNum,"K0"), L0=specInfo(scanNum,"L0")
	Variable H1=specInfo(scanNum,"H1"), K1=specInfo(scanNum,"K1"), L1=specInfo(scanNum,"L1")
	if (numtype(H0+K0+L0+H1+K1+L1)==0)
		s.Nrefs = 2
		s.refs[0].h = H0 ;	s.refs[0].k = K0 ;	s.refs[0].l = L0
		s.refs[1].h = H1 ;	s.refs[1].k = K1 ;	s.refs[1].l = L1
	endif

	UpdateDefaultSampleStruct(s)					// update default locations with this value of s
	return 0
End
//
Function/WAVE protoSpec2BL(inverse)
	Variable inverse
//	Make/N=(3,3)/FREE/D mat=(p==q)
	Make/N=(3,3)/FREE/D mat=NaN
	return mat
End
#endif

//  ============================================================================  //
//  ============================ End of Sample Orientation =============================  //




//  ============================================================================  //
//  ========================= Start of Diffractometer Definitions ==========================  //

Structure diffractometerTypeStruct	// structure that identifies a diffractometer type, used with package prefs
	uchar name[100]				// name, e.g. fourc,kappa, ...
	int16 Naxes						// number of axes in this diffractometer
	uchar axesNames[256]			// semicolon separated list of axis names
	double KA						// kappa angle (degree)
EndStructure



// The routines in this section, contain those parts that are unique to each different kind of diffractometer

// CDP:		A = {nu, mu, delta}
//
// Static
ThreadSafe Function/WAVE CDPMatrix(A)
	Wave A
	Variable theta = A[1] * PI/180
	Make/N=(3,3)/D/FREE muRot					// rotates -theta (degree) about y-axis
	Variable c=cos(theta), s=sin(theta)
	muRot[0][0]= {c,0,-s}
	muRot[0][1]= {0,1,0}
	muRot[0][2]= {s,0,c}
	return muRot
End

//	// CDP:		A = {nu, hexu, delta}
////
//Static Constant hexuOffset = 25.1555				// = 50.311 / 2
////
//Static Function/WAVE CDPMatrix(A)
////ThreadSafe Static Function/WAVE CDPMatrix(A)
//	Wave A
//	Variable theta = (A[1] + hexuOffset) * PI/180
//	Make/N=(3,3)/D/FREE hexRot					// rotates -theta (degree) about y-axis
//	Variable c=cos(theta), s=sin(theta)
//	hexRot[0][0]= {c,0,-s}
//	hexRot[0][1]= {0,1,0}
//	hexRot[0][2]= {s,0,c}
//	return hexRot
//End


// Static
ThreadSafe Function/WAVE CDPDetectorMatrix(A)	// returns a rotation matrix describing current detector rotation
	Wave A
	return doubleDetectorMatrix(A[2],A[0])		// delta, nu
End


//ThreadSafe 
Static Function/WAVE CDPSpec2BL(inverse)		// same as for fourc
	Variable inverse								// if true return Inv(spec2BL)
	return fourcSpec2BL(inverse)
End





// fourc:	A = {delta, theta, chi, phi}
// Static
ThreadSafe Function/WAVE fourcMatrix(A)
	Wave A
	Variable theta=A[1], chi=A[2], phi=A[3]						// A[0]=2theta

	Make/N=(3,3)/D/FREE thetaM,chiM,phiM

	Variable ct=cos(theta*PI/180), st=sin(theta*PI/180)			// theta rotation, about x-axis
	thetaM[0][0]= {1, 0, 0}
	thetaM[0][1]= {0, ct, -st}
	thetaM[0][2]= {0, st, ct}

	Variable cc=cos(chi*PI/180), sc=sin(chi*PI/180)				// chi rotation, about z-axis
	chiM[0][0]= {cc, sc, 0}
	chiM[0][1]= {-sc, cc, 0}
	chiM[0][2]= {0, 0,1}

	Variable cp=cos(phi*PI/180), sp=sin(phi*PI/180)			// phi rotation,  about x-axis
	phiM[0][0]= {1, 0, 0}
	phiM[0][1]= {0, cp, -sp}
	phiM[0][2]= {0, sp, cp}

	MatrixOP/FREE sampleRot = thetaM x chiM x phiM
	WaveClear thetaM,chiM,phiM
	return sampleRot
End
//
//ThreadSafe Static Function/WAVE fourcDetectorMatrix(A)	// returns a rotation matrix describing current detector rotation
ThreadSafe Function/WAVE fourcDetectorMatrix(A)			// returns a rotation matrix describing current detector rotation
	Wave A
	return simpleVerticalDetectorMatrix(A[0])	
End


Static Function fourcReferenceOrientation(s,px,py)
	STRUCT sampleStructure &s
	Variable px,py

	s.name = "spec default values"
	s.Nrefs = 2
	s.refs[0].h=0 ;			s.refs[0].k=1 ;	s.refs[0].l=0		// primary reflection
	s.refs[0].A[0] = 60	;	s.refs[0].A[1] = 30	;	s.refs[0].A[2] = 0	;	s.refs[0].A[3] = 0
	s.refs[1].h=-1 ;		s.refs[1].k=0 ;	s.refs[1].l=0	// secondary reflection
	s.refs[1].A[0] = 60	;	s.refs[1].A[1] = 30	;	s.refs[1].A[2] = -90;	s.refs[1].A[3] = 0

	s.xtal.desc = "fourc default"
	s.xtal.SpaceGroup = 221
	s.xtal.a = 0.154
	LatticeSym#ForceLatticeToStructure(s.xtal)
	s.refs[0].lambda = 0.154
	s.refs[1].lambda = 0.154
	s.refs[0].px = px ;		s.refs[0].py = py		// using a point detector, NOT an area detector
	s.refs[1].px = px ;		s.refs[1].py = py		// using a point detector, NOT an area detector
	return 0
End



//ThreadSafe 
Static Function/WAVE fourcSpec2BL(inverse)
	Variable inverse						// if true return Inv(spec2BL)
	Make/N=(3,3)/D/FREE rotMat
	if (inverse)							// sets rotMat = Inv(spec2BL)
		rotMat[0][0]= {0,0,1}				// qSPEC = rotMat x qBL,    qBL = spec2BL x qSPE
		rotMat[0][1]= {1,0,0}
		rotMat[0][2]= {0,1,0}
	else										// set rotMat to spec2BL
		rotMat[0][0] = {0,1,0}				//  qBL = spec2BL x qSPEC,   and rotMat = spec2BL
		rotMat[0][1] = {0,0,1}
		rotMat[0][2] = {1,0,0}
	endif
	return rotMat
End


// kappa:	A = {delta, theta, kappa, phi, mu, nu}
//
//Function testKappa()
//	Variable theta,kappa,phi, tth
//	Variable h=0,k=0,l=1
//
//	Make/N=3/D/FREE A={nan,theta,kappa,phi}
//
//	Make/N=3/FREE/D qvec
//	hkl2Q(h,k,l, qvec)
//	Wave Tmat = kappaMatrix(A)
//
//	MatrixOP/FREE/O qvec = Tmat x qvec
//	qvec = abs(qvec) < 1e-14 ? 0 : qvec
//	printWave(qvec,name="Q rotated")
//End
////
//Static Constant KA = 50
//ThreadSafe Static Function/WAVE kappaMatrix(A)
ThreadSafe Function/WAVE kappaMatrix(A)
	Wave A
	Variable theta=A[1], kappa=A[2], phi=A[3]						// A[0]=2theta

	Variable KA=NumVarOrDefault("root:Packages:Diffractometer:KA",50)// kappa angle
	Make/N=(3,3)/D/FREE thetaM,kappaM,phiM

	Variable ct=cos(theta*PI/180), st=sin(theta*PI/180)			// theta rotation
	thetaM[0][0]= {1, 0, 0}
	thetaM[0][1]= {0, ct, -st}
	thetaM[0][2]= {0, st, ct}

	Variable cKA=cos(KA*PI/180), sKA=sin(KA*PI/180)			// kappa rotation
	Variable ck=cos(kappa*PI/180), sk=sin(kappa*PI/180)
	kappaM[0][0]= {ck+cKA*(1-ck),  cKA*sKA*(1-ck), -sKA*sk}
	kappaM[0][1]= {cKA*sKA*(1-ck), ck+sKA*sKA, cKA*sk}
	kappaM[0][2]= {sKA*sk, -cKA*sk, ck}

	Variable cp=cos(phi*PI/180), sp=sin(phi*PI/180)			// phi rotation
	phiM[0][0]= {1, 0, 0}
	phiM[0][1]= {0, cp, -sp}
	phiM[0][2]= {0, sp, cp}

	MatrixOP/FREE sampleRot = thetaM x kappaM x phiM
	WaveClear thetaM,kappaM,phiM
	return sampleRot
End
//
// Static
ThreadSafe Function/WAVE kappaDetectorMatrix(A)// returns a rotation matrix describing current detector rotation
	Wave A
	return doubleDetectorMatrix(A[0],A[5])			// delta, nu
End
//Static Function/WAVE kappaDetectorMatrix(A)		// returns a rotation matrix describing current detector rotation
//	Wave A
//	return simpleVerticalDetectorMatrix(A[0])			// delta
//End


//ThreadSafe 
Static Function/WAVE kappapec2BL(inverse)			// same as for fourc
	Variable inverse									// if true return Inv(spec2BL)
	return fourcSpec2BL(inverse)
End



// This is a simple detector that moves vertically up in 2-theta, used in many diffractometers
// Static
ThreadSafe Function/WAVE simpleVerticalDetectorMatrix(A0)// returns a rotation matrix describing current detector rotation
	Variable A0										// vertical detector angle (degree)
	Variable tth=A0*PI/180						// 2theta (rad), this is the angle of the arm, not actual scattering angle
	Make/N=(3,3)/D/FREE rot
	Variable c=cos(tth), s=sin(tth)					// positive rotation is up from +z (at 0¡) to +y (at A0=90¡)
	rot[0][0] = {1, 0, 0}
	rot[0][1] = {0, c,-s}
	rot[0][2] = {0, s, c}
	return rot
End


// This is a 2-axis detector that has an upward moving axis mounted on a horizontal table
// Static
ThreadSafe Function/WAVE doubleDetectorMatrix(delta,nu)// returns a rotation matrix describing current detector rotation
	Variable delta							// a horizontal axis that moves detector upward (degree)
	Variable nu								// holds delta rotates about a vertical axis, moves delta sideways (degree)

	Variable c=cos(delta*PI/180), s=sin(delta*PI/180)
	Make/N=(3,3)/D/FREE rotDelta		// positive rotation is up from +z (at 0¡) to +y (at delta=90¡)
	rotDelta[0][0] = {1, 0, 0}				// rotation matrix about x-axis
	rotDelta[0][1] = {0, c,-s}
	rotDelta[0][2] = {0, s, c}

	c = cos(nu*PI/180)					// positive rotation is from +z (at 0¡) to +x (at nu=90¡)
	s = sin(nu*PI/180)
	Make/N=(3,3)/D/FREE rotNu
	rotNu[0][0] = {c, 0,-s}				// rotation matrix about y-axis
	rotNu[0][1] = {0, 1, 0}
	rotNu[0][2] = {s, 0, c}

	MatrixOP/FREE rot = rotNu x rotDelta
	WaveClear rotNu, rotDelta
	return rot
End



ThreadSafe Function/WAVE protoMatrix(A)		// used for both the sample and detector matricies
	Wave A
	return $""
End

Function protoReferenceOrientation(s,px,py)	// fills s with the reference orientation for an ideal sample
	STRUCT sampleStructure &s
	Variable px,py

	s.name = "proto default values"
	s.Nrefs = 0
	s.xtal.desc = "default Si"
	s.xtal.SpaceGroup = 227
	s.xtal.a = 0.54310206
	LatticeSym#ForceLatticeToStructure(s.xtal)
	Variable astar = 2*PI/(s.xtal.a)
	s.UB[0] = astar	;	s.UB[1] = 0		;	s.UB[2] = 0
	s.UB[3] = 0		;	s.UB[4] = astar	;	s.UB[5] = 0
	s.UB[6] = 0		;	s.UB[7] = 0	;		s.UB[8] = astar
	return 1									// flags an error, the values above are kind of pointless
End

//  ============================================================================  //
//  ========================== End of Diffractometer Definitions ==========================  //




//  ============================================================================  //
//  ============================ Start of Detector Orientation ===========================  //

//	Note, for a point detector, just configure the detector to be 1x1.  So in reality, everyhing is treated as an area detector
//
// This section is used with an area detector

//	px = startx + px*groupx + (groupx-1)/2	// pixel is zero based here & startx is zero based
//	py = starty + py*groupy + (groupy-1)/2	// groupx=1 is un-binned


Structure detectorGeometrys		// structure definition for a detector
	int16 N							// number of detectors defined (must be <= MAX_DETECTORS)
	int16 last						// last detector index (probably the one you want)
	uchar diffractometer[100]		// name of diffractometer (e.g. fourc)
	Struct detectorGeometry d[MAX_DETECTORS]
EndStructure

Structure detectorGeometry			// structure definition for a detector
	int16 used						// TRUE=detector used, FALSE=detector un-used
	int32 Nx, Ny					// # of un-binned pixels in full detector
	double sizeX,sizeY				// outside size of detector (sizeX = Nx*pitchX), measured to outer edge of outer pixels (mm)
	double R[3]						// rotation vector (length is angle in radians)
	double P[3]						// translation vector (mm)

	uchar timeMeasured[100]		// when this geometry was calculated
	uchar geoNote[100]				// note
	uchar detectorID[100]			// unique detector ID
	uchar distortionMapFile[100]	// name of file with distortion map

	double rho00, rho01, rho02	// rotation matrix internally calculated from R[3]
	double rho10, rho11, rho12
	double rho20, rho21, rho22

	double px0, py0					// pixel which is hit by direct beam, px0 & py0 are redundant, but convienent to have
EndStructure


Function DetectorGeometryLocate()					// This only used by FunctionPath("DetectorGeometryLocate")
	return 0										//   to find the path to this file
End


// ThreadSafe 
Static Function QofPixelBL(d,keV,px,py,A,qvecIN,[depth])	// returns Q of pixel in beam line frame, This is NOT q in sample
	STRUCT detectorGeometry, &d
	Variable keV									// for keV<=0, just return qvec with length 1
	Variable px,py									// pixel on detector
	Wave A											// wave with diffractometer angles
	Wave qvecIN									// 3-vector to recieve the result
	Variable depth									// depth of point (measured along Z from origin)
	depth = ParamIsDefault(depth) ? 0 : depth
	depth = numtype(depth) ? 0 : depth
	keV = keV>0 ? keV : NaN

	Variable klen = keV>0 ? 2*PI*keV/hc_keVnm : 1// 2*PI/lambda
	Make/N=3/D/FREE kf, ki={0,0,1}				// kf points to position on detector (mm)
	pixel2xyz(d,px,py,A,kf)						// pixel (px,py) on detector
	kf[2] -= depth									// depth correction
	normalize(kf)

	//	MatrixOP/FREE  qout = T x UB x recip x hkl

	Variable theta = acos(MatrixDot(kf,ki))/2	// ki.kf = cos(2theta), (radians)
	if (WaveExists(qvecIN))
		qvecIN = (kf - ki)*klen						// q in beam line frame  (kf - ki)
		if (!(keV>0))
			normalize(qvecIN)						// no energy, so make a unit vector
		endif
	endif
	WaveClear ki,kf
	return theta
End
//
// Vectorized version of QofPixelBL
ThreadSafe Static Function/WAVE QofPixelBLVEC(d,keV,pxpy,A,[depth])	// returns Q of pixel in beam line frame, This is NOT q in sample
	STRUCT detectorGeometry, &d
	Variable keV									// for keV<=0, just return qvecs with length 1
	Wave pxpy										// array of pixels to use (must be in raw un-binned pixels), perhaps an ROI
	Wave A											// wave with diffractometer angles
	Variable depth									// depth of point (measured along Z from origin)
	depth = ParamIsDefault(depth) ? 0 : depth
	depth = numtype(depth) ? 0 : depth
	keV = keV>0 ? keV : NaN

	Variable N=DimSize(pxpy,0)					// number of pixels to process
	Variable klen = keV>0 ? 2*PI*keV/hc_keVnm : 1// 2*PI/lambda
	Make/N=3/D/FREE ki={0,0,klen}				// incident wave vector, length=2PI/lambda

	Wave kf = pixel2xyzVEC(d,pxpy,A)				// position (xyz BL coords) of every pixel
	if (depth==0)
		MatrixOP/FREE/O kf = NormalizeRows(kf)*klen
	else
		Make/N=3/D/FREE depthVec = {0,0,depth}
		MatrixOP/FREE/O kf = NormalizeRows(kf - rowRepeat(depthVec,N))*klen	// shift z component by depth
		WaveClear depthVec
	endif

	if (!(keV>0))
		MatrixOP/FREE qvecs = NormalizeRows(kf - rowRepeat(ki,N))	// no energy, so make a unit vector
	else
		MatrixOP/FREE qvecs = kf - rowRepeat(ki,N)
	endif
	WaveClear ki,kf
	//	MatrixOP/FREE  qout = T x UB x recip x hkl
	return qvecs
End


Static Function pixel2xyz(d,px,py,A,xyz,[DeltaPixelCorrect])	// convert pixel position to the beam line coordinate system
	STRUCT detectorGeometry, &d
	Variable px,py									// pixel position on detector (full chip & zero based)
	Wave A
	Wave xyz										// 3-vector to receive the result, position in beam line coords (mm)
	Wave DeltaPixelCorrect							// contains optional pixel corrections dimensions are [Nx][Ny][2]

	String DiffractometerName = StrVarOrDefault("root:Packages:Diffractometer:DiffractometerName",DEFAULT_DIFF_NAME)
	Variable xp,yp, zp								// x' and y' (requiring z'=0), detector starts centered on origin and perpendicular to z-axis
	xp = (px - 0.5*(d.Nx-1)) * d.sizeX/d.Nx		// (x' y' z'), position on detector
	yp = (py - 0.5*(d.Ny-1)) * d.sizeY/d.Ny

	if (!ParamIsDefault(DeltaPixelCorrect) && DimSize(DeltaPixelCorrect,2)==2)
		Variable ix=limit(round(px),0,(d.Nx-1)),iy=limit(round(py),0,(d.Ny-1))
		px += DeltaPixelCorrect[ix][iy][0]// start as pixel position on detector + correction for distortion
		py += DeltaPixelCorrect[ix][iy][1]
	endif
	//	//	Make/N=3/D/FREE vec= { ( px-(Nx/2) ) * dx, ( py-(Ny/2) ) * dy, dist }	// position of pixel at nu=0
	//		Make/N=3/D/FREE vec= { ( (Nx/2)-px ) * dx, ( py-(Ny/2) ) * dy, dist }	// position of pixel at nu=0

	xp += d.P[0]									// translate by P
	yp += d.P[1]
	zp = d.P[2]

	xyz[0] = d.rho00*xp + d.rho01*yp + d.rho02*zp	// xyz = rho x [ (x' y' z') + P ]
	xyz[1] = d.rho10*xp + d.rho11*yp + d.rho12*zp	// rho is pre-calculated from vector d.R
	xyz[2] = d.rho20*xp + d.rho21*yp + d.rho22*zp

	// xyz is now the position of pixel with all angles at zero, so rotate detector by tth
	String dname = DiffractometerName+"DetectorMatrix"
	FUNCREF protoMatrix detectorMatrix = $dname
	Wave rot = detectorMatrix(A)					// rotation of detector
	MatrixOP/FREE/O vec = rot x xyz				// rotate vec by tth, to the actual detector position
	xyz = vec

	Variable tthTotal, dot=vec[2]/norm(vec)		// calculate total scattering angle
	if (dot >= 1)
		tthTotal = 0
	elseif (dot <= -1)
		tthTotal =180
	else
		tthTotal = acos(dot)*180/PI
	endif
	return tthTotal
End
//
// Vectorized version of pixel2xyz
ThreadSafe Static Function/WAVE pixel2xyzVEC(d,pxpy,A,[DeltaPixelCorrect])		// convert pixel position to the beam line coordinate system
	STRUCT detectorGeometry, &d
	Wave pxpy										// array of pixels to use (must be in raw un-binned pixels), perhaps an ROI, first dimension is number of pixels, second is x,y
	Wave A
	Wave DeltaPixelCorrect							// contains optional pixel corrections dimensions are [N][2]

	String DiffractometerName = StrVarOrDefault("root:Packages:Diffractometer:DiffractometerName",DEFAULT_DIFF_NAME)
	Variable N=DimSize(pxpy,0), Nx=d.Nx, Ny=d.Ny
	Variable ixc=(Nx-1)/2,  iyc=(Ny-1)/2		// center of detector (pixles)

	String dname = DiffractometerName+"DetectorMatrix"
	FUNCREF protoMatrix detectorMatrix = $dname
	Wave rot = detectorMatrix(A)					// rotation of detector
	if (!WaveExists(rot))
		return $""
	endif

	Make/N=3/D/FREE Pw, dxyz, ixyzc={ixc,iyc,0}
	dxyz = {d.sizeX/Nx, d.sizeY/Ny, 0}
	Pw = d.P[p]

	Make/N=(3,3)/D/FREE rho
	rho[0][0] = d.rho00 ;	rho[0][1] = d.rho01 ;	rho[0][2] = d.rho02	// xyz = rho x [ (x' y' z') + P ]
	rho[1][0] = d.rho10 ;	rho[1][1] = d.rho11 ;	rho[1][2] = d.rho12	// rho is pre-calculated from vector d.R
	rho[2][0] = d.rho20 ;	rho[2][1] = d.rho21 ;	rho[2][2] = d.rho22
	// xyz is now the position of pixel with all angles at zero, so rotate detector by tth

	Make/N=(N,3)/FREE xyz						// array of 3-vector to receive the result, position in beam line coords (mm)
	if (!ParamIsDefault(DeltaPixelCorrect) && DimSize(DeltaPixelCorrect,1)==2)
		xyz[][0,1] = pxpy[p][q] + DeltaPixelCorrect[p][q]	// start as pixel position on detector + correction for distortion
	else
		xyz[][0,1] = pxpy[p][q]					// start as pixel position on detector
	endif

	//	xyz = (xyz - ixyzc[q]) * dxyz[q] + Pw[q]		// (x' y' z'), position on detector for all pixels at A[]=0
	//	MatrixOP/FREE/O xyz = (xyz - rowRepeat(ixyzc,N)) * rowRepeat(dxyz,N) +rowRepeat(Pw,N)		// (x' y' z'), position on detector for all pixels at A[]=0
	//	MatrixOP/FREE/O xyz = (rot x rho x xyz^t)^t			// rotate detector by calib, and rotate by tth, to the final detector position
	MatrixOP/FREE/O xyz = (rot x rho x ((xyz - rowRepeat(ixyzc,N)) * rowRepeat(dxyz,N) +rowRepeat(Pw,N))^t)^t		// (x' y' z'), position on detector for all pixels at A[]=0
	WaveClear rho, rot, Pw, dxyz, ixyzc
	return xyz
End



Static Function/WAVE MakePixelListFromImage(image,[pxpyName])
	Wave image
	String pxpyName					// name to use for wave, when pxpyName passed, then create a real (NOT FREE) wave
	pxpyName = SelectString(ParamIsDefault(pxpyName),pxpyName,"")

	String wnote=note(image)
	Variable N=numpnts(image), Nx=DimSize(image,0)
	Variable startx=NumberByKey("startx",wnote,"="), groupx=round(NumberByKey("groupx",wnote,"="))
	Variable starty=NumberByKey("starty",wnote,"="), groupy=round(NumberByKey("groupy",wnote,"="))
	startx = (startx<0 || numtype(startx)) ? 0 : startx
	starty = (starty<0 || numtype(starty)) ? 0 : starty
	groupx = (groupx<1 || numtype(groupx)) ? 1 : groupx
	groupy = (groupy<1 || numtype(groupy)) ? 1 : groupy
	if (strlen(pxpyName)==0)
		Make/N=(N,2)/FREE pxpy					// all of the pixel values desired, FREE wave
	else
		pxpyName = CleanupName(pxpyName,0)
		Make/N=(N,2) $pxpyName/WAVE=pxpy	// all of the pixel values desired, NOT FREE wave
	endif
	pxpy[][0] = mod(p,Nx)*groupx+startx			// units are full-chip unbinned, even if pixel are from a binned ROI
	pxpy[][1] = floor((p+0.4)/Nx)*groupy+starty
	return pxpy
End



Static Function DetectorUpdateCalc(d)				// update all internally calculated things in the detector structure
	STRUCT detectorGeometry &d
	if (!(d.used))
		return 1
	endif

	Variable Rx, Ry, Rz								// used to make the rotation matrix rho from vector R
	Variable theta, c, s, c1
	Variable i
	Rx=d.R[0]; Ry=d.R[1]; Rz=d.R[2]				// make the rotation matrix rho from vector R
	theta = sqrt(Rx*Rx+Ry*Ry+Rz*Rz)
	if (theta==0)									// no rotation, set to identity matrix
		d.rho00 = 1;		d.rho01 = 0;		d.rho02 = 0
		d.rho10 = 0;		d.rho11 = 1;		d.rho12 = 0
		d.rho20 = 0;		d.rho21 = 0;		d.rho22 = 1
		return 0
	endif

	c=cos(theta)
	s=sin(theta)
	c1 = 1-c
	Rx /= theta;	Ry /= theta;	Rz /= theta		// make |{Rx,Ry,Rz}| = 1

	d.rho00 = Rx*Rx*c1 + c;		d.rho01 = Rx*Ry*c1 - Rz*s;	d.rho02 = Rx*Rz*c1 + Ry*s	// this is the Rodrigues formula from:
	d.rho10 = Rx*Ry*c1 + Rz*s;	d.rho11 = Ry*Ry*c1 + c;		d.rho12 = Ry*Rz*c1 - Rx*s	// http://mathworld.wolfram.com/RodriguesRotationFormula.html
	d.rho20 = Rx*Rz*c1 - Ry*s;	d.rho21 = Ry*Rz*c1 + Rx*s ;	d.rho22 = Rz*Rz*c1 + c

	Variable px0,py0
	px0 = 0.5*(d.Nx-1) - ((d.Nx/d.sizeX) * (d.P[0]))			// px0 & py0 are redundant, but convienent to have
	py0 = 0.5*(d.Ny-1) - ((d.Ny/d.sizeY) * (d.P[1]))
	if (px0>=0 && px0<(d.Nx) && py0>=0 && py0<(d.Ny))
		d.px0 = px0
		d.py0 = py0
	else
		d.px0 = NaN
		d.py0 = NaN
	endif
	return 0
End
//
Static Function UpdateDefaultDetectorStruct(ds,[local])			// Update the default location with values in ds
	STRUCT detectorGeometrys &ds								// returns 0 if something set, 0 is nothing done
	Variable local												// if true will also update detStructStr in the current folder if it already exists (It will NOT create local detStructStr)
	local = ParamIsDefault(local) ? 0 : local
	local = numtype(local) ? 0 : !(!local)

	Variable i,N=0
	for (i=0;i<MAX_DETECTORS;i+=1)
		DetectorUpdateCalc(ds.d[i])
		N += !(!(ds.d[i].used))
	endfor
	ds.N = N

	String detStructStr
	StructPut/S/B=2 ds, detStructStr
	String/G root:Packages:Diffractometer:detectorStructStr=detStructStr	// always save to default location
	SavePackagePreferences/FLSH=1 "Diffractometer","detectorPrefs",0,ds	// alwasy update prefs too

	local = local || (exists(":detStructStr")==2)				// save locally if told to, or if a local already exists in current datafolder
	if (local)
		String/G :detectorStructStr=detStructStr				// make a local copy too
	endif
	return 0
End

Function FillDetectorsStruct(ds)								//fill the detector structures with current values
	STRUCT detectorGeometrys &ds								// returns 0 if something set, 0 is nothing done

	String strStruct=StrVarOrDefault(":detectorStructStr","")	// set to values in current directory
	if (strlen(strStruct)<1)
		strStruct=StrVarOrDefault("root:Packages:Diffractometer:detectorStructStr","")	// try the default values
	endif
	if (strlen(strStruct)>1)									// found structure information in this experiment, load into ds
		StructGet/S/B=2 ds, strStruct
	else															// no detector in this experiment, try the prefs
		LoadPackagePreferences/MIS=1 "Diffractometer","detectorPrefs",0,ds
		if (V_flag)
			return 1											// did nothing, nothing found
		endif
	endif
	UpdateDefaultDetectorStruct(ds)							// this ensures that if we loaded it from prefs, that local values are updated
	return 0
End

Function FillDetectorStructDefault(d,id)						//fill the detector structure with current values
	STRUCT detectorGeometry &d								// returns 0 if something set, 0 is nothing done
	String id													// detector id, if "", then use ds.last

	STRUCT detectorGeometrys ds								// returns 0 if something set, 0 is nothing done
	if (FillDetectorsStruct(ds))
		return 1
	endif

	Variable i = ds.last											// default to pass back if no id specified
	if (strlen(id)>0)											// search for detector with matching id
		for (i=0;i<MAX_DETECTORS;i+=1)
			if (ds.d[i].used && StringMatch(ds.d[i].detectorID,id))	// found detector with correct ID
				break
			endif
		endfor
	endif

	if (i>=0 && i<MAX_DETECTORS)
		CopyOneDetectorGeometry(d,ds.d[i])
	else
		return 1
	endif
	return 0
End
//
Static Function FindDetectorIdIndex(ds,id,[used])				// find index into ds.d[i] for id, returns -1 if not found
	STRUCT detectorGeometrys &ds
	String id
	Variable used												// if true, only find used detectors
	used = ParamIsDefault(used) ? 0 : used
	used = numtype(used) ? 0 : !(!used)
	if (strlen(id)<1)
		return -1
	endif
	Variable i
	for (i=0;i<MAX_DETECTORS;i+=1)
		if (ds.d[i].used && StringMatch(ds.d[i].detectorID,id))	// found detector with correct ID and it is used
			return i
		elseif (!used && StringMatch(ds.d[i].detectorID,id))		// correct ID and don't care if it is used
			return 1
		endif
	endfor
	return -1													// failed to find id
End
//
Static Function NextEmptyDetectorIndex(ds)						// find index into ds.d[i] for first empty slot, returns -1 if all are full
	STRUCT detectorGeometrys &ds
	String id

	Variable i
	for (i=0;i<MAX_DETECTORS;i+=1)
		if (!(ds.d[i].used))										// found a detector that is not used
			return i
		endif
	endfor
	return -1													// failed to find id
End



Static Function UpdateDetectorInList(ds,d,force)			// update a detector in a list of detector geometrys
	STRUCT detectorGeometrys &ds						// list of detector geometrys
	STRUCT detectorGeometry &d						// detector geometry to replace
	Variable force										// flag, if true, then overwrite last detector if list is full

	String id = d.detectorID
	Variable i, last=-1, empty=-1
	for (i=0;i<MAX_DETECTORS;i+=1)
		if (StringMatch(ds.d[i].detectorID,id))			// found detector with correct ID
			break
		endif
		if (!(ds.d[i].used) && empty<0)					// save first empty slot
			empty = i
		endif
	endfor

	if (i<MAX_DETECTORS)
		last = i											// overwrite existing detector with same ID
	elseif (empty>=0)
		last = empty									// overwrite first empty slot
	elseif (force)
		last = MAX_DETECTORS-1						// overwrite last detector in list
	endif
	if (last>=0)
		DetectorUpdateCalc(d)
		CopyOneDetectorGeometry(ds.d[last],d)
		ds.last = last
	endif
	return last
End



Function SetDefaultDetector2Reference(point)		// set default detector to the reference values
	Variable point												// if TRUE then a point detector
	if (numtype(point))
		Prompt point,"Default Detector Type",popup,"Area;Point"
		DoPrompt "Detector Type", point
		if (V_flag)
			return 1
		endif
		point -= 1
		printf "SetDefaultDetector2Reference(%g)\r",point
	endif

	STRUCT detectorGeometrys ds							// returns 0 if something set, 0 is nothing done
	if (FillDetectorsStruct(ds))							//fill the detector structures with current values
		return 1
	endif

	STRUCT detectorGeometry d
	DetectorReferenceOrientation(d,point)				// get reference values for a detector
	ds.diffractometer = StrVarOrDefault("root:Packages:Diffractometer:DiffractometerName",DEFAULT_DIFF_NAME)	// name of diffractometer (e.g. fourc)

	UpdateDetectorInList(ds,d,0)							// does a DetectorUpdateCalc(), and stores new values where they will be found
	UpdateDefaultDetectorStruct(ds)						// does a DetectorUpdateCalc(), and stores new values where they will be found
	return 0
End
//
Static Function DetectorReferenceOrientation(d,point)			// sets d to the reference orientation (sort of an ideal set of values)
	STRUCT detectorGeometry &d
	Variable point												// if TRUE then a point detector

	point = numtype(point) ? 0 : !(!point)
	if (point)
		// define a simple point detector 1m from sample, with a size of 1cm square
		d.used = 1
		d.Nx = 1 ;			d.Ny = 1							// number of un-binned pixels in whole detector
		d.sizeX = 10;		d.sizeY = 10						// outside size of detector (mm)
		d.R[0]=0;			d.R[1]=0;			d.R[2]=0		// angle of detector, theta = 0¡
		d.P[0]=0;			d.P[1]=0;			d.P[2]=1000	// offset to detector (mm)
		d.timeMeasured = ""
		d.geoNote = ""
		d.detectorID = "Point Detector"
		d.distortionMapFile = ""
	else
		// define a Pilatus 100K Detector located at 1m from sample with pixels along x and y
		d.used = 1
		d.Nx = 487 ;				d.Ny = 195					// number of un-binned pixels in whole detector
		d.sizeX = 487*0.172;		d.sizeY = 195*0.172		// outside size of detector (mm)
		d.R[0]=0;			d.R[1]=0;			d.R[2]=0		// angle of detector @ theta=0¡
		d.P[0]=0;			d.P[1]=0;			d.P[2]=1000	// offset to detector (mm)
		d.timeMeasured = ""
		d.geoNote = ""
		d.detectorID = DEF_PILATUS_ID
		d.distortionMapFile = ""
	endif
	DetectorUpdateCalc(d)										// tidy up and do all pre-calculations
End



Function PrintCurrentDetectors([more])			// prints the current default detector geometrys to the history
	Variable more
	STRUCT detectorGeometrys ds
	FillDetectorsStruct(ds)							//fill the geometry structure with all known detector values
	printf "using diffractometer '%s'\r",SelectString(strlen(ds.diffractometer),"UNKNOWN",ds.diffractometer)
	Variable i
	for (i=0;i<MAX_DETECTORS;i+=1)
		if (ParamIsDefault(more))
			printOneDetector(ds.d[i])
		else
			printOneDetector(ds.d[i],more=more)
		endif
	endfor
End
//
Static Function printOneDetector(d,[more])			// print the details for passed detector geometry to the history window
	STRUCT detectorGeometry &d
	Variable more
	if (!(d.used))
		return 1
	endif
	more = ParamIsDefault(more) ? 0 : more
	more = numtype(more) ? 0 : !(!more)

	if (DetectorBad(d))
		print "************************************************************"
		print "        BAD   BAD   BAD   BAD   BAD   BAD   BAD   BAD   BAD   BAD   BAD   BAD   BAD"
		print "************************************************************"
	endif

	printf "	Nx=%d, Ny=%d						// number of un-binned pixels in detector\r",d.Nx,d.Ny
	printf "	sizeX=%g, sizeY=%g			// size of detector (mm)\r",(d.sizeX), (d.sizeY)
	printf "	   pixel size is %g x %g mm\r",(d.sizeX)/(d.Nx), (d.sizeY)/(d.Ny)
	printf "	R = {%.8g, %.8g, %.8g}, a rotation of %.7g¡			// rotation vector\r",d.R[0],d.R[1],d.R[2],sqrt(d.R[0]*d.R[0] + d.R[1]*d.R[1] + d.R[2]*d.R[2])*180/PI
	printf "	P = {%g, %g, %g}							// translation vector (mm)\r", d.P[0], d.P[1], d.P[2]
	if (numtype((d.px0)+(d.py0))==0)
		printf "	    incident beam hits pixel {%g, %g}\r", d.px0, d.py0
	endif

	if (strlen(d.timeMeasured))
		printf "	detector set on  '%s'\r",d.timeMeasured
	endif
	if (strlen(d.geoNote))
		printf "	detector note = '%s'\r",d.geoNote
	endif
	if (strlen(d.distortionMapFile))
		printf "	detector distortion file = '%s'\r",d.distortionMapFile
	endif
	printf "	detector ID = '%s'\r"d.detectorID
	Make/N=(3,3)/D/FREE rho
	rho[0][0] = d.rho00	;	rho[0][1] = d.rho01	;		rho[0][2] = d.rho02
	rho[1][0] = d.rho10	;	rho[1][1] = d.rho11	;		rho[1][2] = d.rho12
	rho[2][0] = d.rho20	;	rho[2][1] = d.rho21	;		rho[2][2] = d.rho22
	rho = abs(rho)<1e-11 ? 0 : rho
	if (more)
		printf "| rho | = %g\t\t// rotation matrix from R\r",MatrixDet(rho)
		printWave(rho,name="",brief=1)
	endif
	rho = abs(rho)<1e-6 ? 0 : rho
	Make/N=3/FREE hatx=rho[p][0], haty=rho[p][1], hatD=d.P[p]
	MatrixOP/FREE/O hatD = rho x hatD
	printf "	pixels along +x --> %s,    pixels along +y --> %s\r",vec2str(hatx),vec2str(haty)
	printf "	center of detector = %s,    (at all angles==0)\r",vec2str(hatD)
	return 0
//	Variable chiMin, chiMax, tthMin, tthMax
//	Variable/C chiZ, tthZ
//	chiZ = chiRange(d)
//	chiMin = min(chiMin,real(chiZ))
//	chiMax = max(chiMax,imag(chiZ))
//	tthZ =  tthRange(d)
//	tthMin = min(tthMin,real(tthZ))
//	tthMax = max(tthMax,imag(tthZ))
//	printf "\tangle ranges:  chi = [%g, %g¡],  2th = [%g, %g¡]\r", real(chiZ),imag(chiZ),real(tthZ),imag(tthZ)
//	return 0
End



Function LoadDetectorGeoGeneral(method)		// General routine for loading detector info
	Variable method
	
	String list="Select Local File;from Web;by Date from Local Directory"
	if (!(method>=0 && method<ItemsInList(list)))
		method = 0
		Prompt method "Select a new geomety by...",popup,list
		DoPrompt "Select Geometry", method
		if (V_flag)
			return 1
		endif
		method -= 1
	endif

	Variable err=1
	if (method==0)
		err = LoadDetectorsFromFile("","")
	elseif (method==1)
		err = LoadDetectorFromWeb("","")
	elseif (method==2)
		err = LoadDetectorsByDateLocal("","")
	else
		DoAlert 0,"Unknown, method = "+num2str(method)+"\rNothing done."
	endif
	return err
End


Function LoadDetectorFromWeb(dateStr,timeStr,[tagName,epoch,printIt])
	String dateStr			// something like "4/2/2010"
	String timeStr			// something like "13:22:00"
	String tagName			// name of tag, defaults to "geoD"
	Variable epoch			// desired epoch
	Variable printIt

	tagName = SelectString(ParamIsDefault(tagName),tagName,defaultGeoTagName)	// default to defaultGeoTagName
	tagName = SelectString(strlen(tagName),defaultGeoTagName,tagName)	// default to defaultGeoTagName
	epoch = ParamIsDefault(epoch) ? NaN : epoch
	printIt = ParamIsDefault(printIt) ? (strlen(GetRTStackInfo(2))==0) : printIt
	printIt = numtype(printIt) ? (strlen(GetRTStackInfo(2))==0) : printIt

	dateStr = dateStr2ISO8601Str(dateStr)
	Variable err=1

	Variable year,month,day, hour,minute,second
	sscanf dateStr+"T"+timeStr,"%d/%d/%dT%d:%d:%d",year,month,day,hour,minute,second
	if (!(epoch>0) && V_flag==6)			// no valid epoch passed, try to interpret date & time
		epoch = date2secs(year,month,day)+3600*hour+60*minute+second
	endif

	if (!(epoch>0))							// unable to find desired date & time, so ask
		epoch = getEpochManually(dateStr,timeStr)
	endif
	if (numtype(epoch) || epoch<=0)
		return err
	endif
	dateStr = Secs2Date(epoch,-2)
	timeStr = Secs2Time(epoch,3)

	String url
	sprintf url "http://%s/index.php?tag=%s&date=%sT%s",GeoWebServer,tagName,dateStr,timeStr
	//	printf "url = %s\r\r",url
	String result = FetchURL(url)
	if (strsearch(result,"ERROR",0)==0)
		printf "result = '%s'\r", ReplaceString("\n",result,"")
		return err
	endif

	Variable N=strlen(result), dquote=char2num("\"")			// ASCII value of a double-quote
	if (char2num(result[0])==dquote)							// remove leading double-quote from Apple Script
		result = result[1,N-1]
		N -= 1
	endif
	if (char2num(result[N-1])==dquote)						// remove trailing double-quote from Apple Script
		result = result[0,N-2]
		N -= 1
	endif


	String buf = diffractometer#xmlContents(result)
	STRUCT detectorGeometrys ds
	err = diffractometer#DetectorsFromXML(buf,ds)
	if (err)
		DoAlert 0,"UN-able to load a valid detector\rNothing done"
	else
		printf "using diffractometer '%s'\r",SelectString(strlen(ds.diffractometer),"UNKNOWN",ds.diffractometer)
		Variable i
		for (i=0;i<MAX_DETECTORS;i+=1)
			diffractometer#printOneDetector(ds.d[i])
		endfor
		diffractometer#UpdateDefaultDetectorStruct(ds)
	endif
	return err
End


Function LoadDetectorsByDateLocal(dateStr,timeStr,[tagName,epoch,printIt])
	String dateStr			// something like "4/2/2010"
	String timeStr			// something like "13:22:00"
	String tagName			// name of tag, defaults to "geoD"
	Variable epoch			// desired epoch
	Variable printIt

	tagName = SelectString(strlen(tagName),defaultGeoTagName,tagName)	// default to defaultGeoTagName
	epoch = ParamIsDefault(epoch) ? NaN : epoch
	printIt = ParamIsDefault(printIt) ? (strlen(GetRTStackInfo(2))==0) : printIt
	printIt = numtype(printIt) ? (strlen(GetRTStackInfo(2))==0) : printIt

	dateStr = dateStr2ISO8601Str(dateStr)
	Variable year,month,day, hour,minute,second
	sscanf dateStr+"T"+timeStr,"%d/%d/%dT%d:%d:%d",year,month,day,hour,minute,second
	if (!(epoch>0) && V_flag==6)			// no valid epoch passed, try to interpret date & time
		epoch = date2secs(year,month,day)+3600*hour+60*minute+second
	endif

	if (!(epoch>0))							// unable to find desired date & time, so ask
		epoch = getEpochManually(dateStr,timeStr)
	endif
	if (numtype(epoch) || epoch<=0)
		return 1
	endif

	String fileName=getGeoNameFromDirectory(tagName,epoch)	// full path name to the file
	printf "Loading Geometry from: '%s'\r",fileName
	Variable err = LoadDetectorsFromFile(fileName,"geoXMLpath")
	return err
End
//
Static Function getEpochManually(dateStr,timeStr)	// returns epoch from a user entered date & time
	String dateStr								// something like "4/2/2010", defaults to now
	String timeStr								// something like "13:22:00"

	dateStr = dateStr2ISO8601Str(dateStr)

	Variable year,month,day, hour,minute,second
	sscanf dateStr+"T"+timeStr,"%d/%d/%dT%d:%d:%d",year,month,day,hour,minute,second
	Variable epoch=NaN						// desired epoch
	if (V_flag==6)								// no valid epoch passed, try to interpret date & time
		epoch = date2secs(year,month,day)+3600*hour+60*minute+second
	endif

	if (!(epoch>0))							// unable to determine valid date & time, so ask user
		Variable now=DateTime
		Prompt dateStr,"Date, '4/2/2010', empty means now"
		Prompt timeStr,"Time, '13:22:00', empty means now"
		dateStr = SelectString(strlen(dateStr),Secs2Date(now,-2),dateStr)
		timeStr = SelectString(strlen(timeStr),Secs2Time(now,3),timeStr)
		DoPrompt "Desired Time",dateStr,timeStr
		if (V_flag)
			return NaN
		endif
		if (strlen(dateStr)==0 || stringmatch(dateStr,"now"))
			dateStr = Secs2Date(now,-2)	// default to now
		endif
		if (strlen(timeStr)==0 || stringmatch(timeStr,"now"))
			timeStr = Secs2Time(now,3)	// default to now
		endif
		dateStr = dateStr2ISO8601Str(dateStr)
		sscanf dateStr+"T"+timeStr,"%d-%d-%dT%d:%d:%d",year,month,day,hour,minute,second
		if (V_flag==6)
			epoch = date2secs(year,month,day)+3600*hour+60*minute+second
		endif
	endif
	epoch = epoch<=0 || numtype(epoch) ? NaN : epoch
	return epoch
End
//
// returns full file name to geometry based on epoch, searches all files in geoXMLpath starting with tagName
Function/T getGeoNameFromDirectory(tagName,epoch,[ext,printIt])
	String tagName								// name of tag, defaults to "geoD"
	Variable epoch								// desired epoch
	String ext									// extension, probably xml, but maybe txt
	Variable printIt
	ext = SelectString(ParamIsDefault(ext),ext,"")
	ext = SelectString(strlen(ext),".xml",ext)					// default is '.xml'
	tagName = SelectString(strlen(tagName),defaultGeoTagName,tagName)	// default to defaultGeoTagName
	epoch = epoch>0 ? epoch : DateTime
	printIt = ParamIsDefault(printIt) ? (strlen(GetRTStackInfo(2))==0) : printIt
	printIt = numtype(printIt) ? (strlen(GetRTStackInfo(2))==0) : printIt
	if (strlen(ext)>0 && char2num(ext)!=46)	// ensure ext has leaing ","
		ext = "."+ext
	endif

	PathInfo geoXMLpath						// full path to file, V_flag 0 if path does not exist, 1 if it exists.
	if (V_flag==0)								// path does not exist, set to default near microN
		String path = ParseFilePath(1,FunctionPath("DetectorGeometryLocate"),":",1,1)+"geometrys"
		if (strlen(path))
			NewPath/M="Path to "+tagName+"*"+ext+" files"/O/Q geoXMLpath, path
		endif
		PathInfo geoXMLpath
	endif

	if (V_flag)
		DoAlert 1, "Use default location for Geometry files\r(No means select a new location)"
		V_flag = V_flag==1
	endif
	if (V_flag==0)								// need to ask user to find the path
		NewPath/M="Path to "+tagName+"*"+ext+" files"/O/Q geoXMLpath
	endif
	String dirListAll=IndexedFile(geoXMLpath,-1,ext)
	dirListAll = SortList(dirListAll)	// list of all files in directory
	Variable N=ItemsInList(dirListAll)
	if (N<1)
		if (printIt)
			DoAlert 0,"No "+tagName+"*"+ext+" files found"
		endif
		return ""
	endif

	Make/N=(N)/D/FREE times=inf
	Make/N=(N)/T/FREE files=""
	String fmt = tagName+"_%4d-%02d-%02d_%02d-%02d-%02d"+ext, str
	Variable m,i, year,month,day,hr,minute,second
	for (m=0,i=0; i<N; i+=1)				// check each file for valid file name
		str = StringFromList(i,dirListAll)
		sscanf str, fmt,year,month,day,hr,minute,second
		if (V_flag==6)							// this file name is valid, save it
			files[m] = str
			times[m] = date2secs(year,month,day)+3600*hr+60*minute+second
			m += 1
		endif
	endfor
	N = m
	Redimension/N=(N) times, files
	Sort times,times,files

	i = BinarySearch(times,epoch)
	if (i==-1)									// requested time is before any of the geometrys
		if (printIt)
			print "ERROR -- given time is before all geometrys"
		endif
		return ""
	endif
	i = i==-2 ? N-1 : i						// -2 means after to use last one

	PathInfo geoXMLpath						// get full path to file
	return S_path+files[i]
End


Function LoadDetectorsFromFile(fileName,path)
	String fileName							// full path name to the file
	String path								// name of an Igor path to use

	STRUCT detectorGeometrys ds
	Variable err = ReadDetectorsFromFile(fileName,path,ds)
	if (err)
		DoAlert 0,"UN-able to load a valid detector\rNothing done"
	else
		printf "using diffractometer '%s'\r",SelectString(strlen(ds.diffractometer),"UNKNOWN",ds.diffractometer)
		Variable i
		for (i=0;i<MAX_DETECTORS;i+=1)
			printOneDetector(ds.d[i])
		endfor
		UpdateDefaultDetectorStruct(ds)
	endif
	return err
End
//
Static Function ReadDetectorsFromFile(fileName,path,ds)
	String fileName							// full path name to the file
	String path								// name of an Igor path to use
	STRUCT detectorGeometrys &ds

	Variable f
	Open/F=DetectorFileFilters/M="Detector xml file"/P=$path/R/Z=2 f as fileName
	if (strlen(S_fileName)<1 || !f)
		return 1
	endif
	fileName = S_fileName
	FStatus f
	String buf=""
	buf = PadString(buf,V_logEOF,0x20)
	FBinRead f, buf
	Close f

	buf = xmlContents(buf)
	Variable err = DetectorsFromXML(buf,ds)
	return err
End
//
Static Function DetectorsFromXML(buf,ds)
	String buf								// contents of the <geoN> xml
	STRUCT detectorGeometrys &ds

	buf = XMLremoveComments(buf)		// remove all of the comments

	String str, nodes = XMLNodeList(buf)
	if (WhichListItem("diffractometer",nodes)<0)
		print "cannot find diffractometer"
		return 1
	endif

	String diffractometerName = StringByKey("type",XMLattibutes2KeyList("diffractometer",buf),"=")
	ds.diffractometer = SelectString(strlen(diffractometerName),"",diffractometerName)
	if (char2num(ds.diffractometer)<32)
		ds.diffractometer = ""
	endif
	buf = XMLtagContents("diffractometer",buf)
	if (strlen(buf)<1)
		print "cannot find <diffractometer> in xml"
		return 1
	endif

	Variable Ndetectors=NumberByKey("Ndetectors",XMLattibutes2KeyList("Detectors",buf),"=")
	Ndetectors = Ndetectors>0 ? round(Ndetectors) : 0
	if (Ndetectors<=0  || Ndetectors>MAX_DETECTORS)
		sprintf str,"This file contains %g, but only %g are allowed",Ndetectors,MAX_DETECTORS
		print str
		DoAlert 0,str
		return 1
	endif

	Variable i
	for(i=0;i<MAX_DETECTORS;i+=1)
		ds.d[i].used = 0
	endfor

	String dateWritten = XMLtagContents("dateWritten",buf)
	String timeWritten = XMLtagContents("timeWritten",buf)
//	String fileNote = XMLtagContents("fileNote",buf)

	STRUCT detectorGeometry d
	String list
	Variable N
	String detector, detectors=XMLtagContents("Detectors",buf)
	for (i=0;i<Ndetectors;i+=1)
		N = NumberByKey("N",XMLattibutes2KeyList("Detector",detectors,occurance=i),"=")
		if (N>=MAX_DETECTORS || numtype(N))
			sprintf str,"This file contains detector number %g, but max value is only %g",N,MAX_DETECTORS-1
			print str
			DoAlert 0,str
			continue
		endif
		detector = XMLtagContents("Detector",detectors,occurance=i)

		d.used = 1
		d.timeMeasured = XMLtagContents("timeMeasured",detector)
		d.geoNote = XMLtagContents("note",detector)
		d.detectorID = XMLtagContents("ID",detector)
		d.distortionMapFile = XMLtagContents("distortionMap",detector)

		list =  XMLtagContents2List("Npixels",detector)
		d.Nx = str2num(StringFromList(0,list))
		d.Ny = str2num(StringFromList(1,list))

		list =  XMLtagContents2List("size",detector)
		d.sizeX = str2num(StringFromList(0,list))
		d.sizeY = str2num(StringFromList(1,list))

		list =  XMLtagContents2List("R",detector)
		d.R[0] = str2num(StringFromList(0,list))
		d.R[1] = str2num(StringFromList(1,list))
		d.R[2] = str2num(StringFromList(2,list))

		list =  XMLtagContents2List("P",detector)
		d.P[0] = str2num(StringFromList(0,list))
		d.P[1] = str2num(StringFromList(1,list))
		d.P[2] = str2num(StringFromList(2,list))

		DetectorUpdateCalc(d)					// calculate other values
		CopyOneDetectorGeometry(ds.d[N],d)
	endfor

	printf "   Load detector geometry that was written %s, %s",dateWritten,timeWritten
	print " "
	return DetectorBad(d)
End
//
Static Function/S xmlContents(buf)
	String buf

	Variable i1, i0=strsearch(buf,"<?xml",0)	// find start of header tag
	if (i0<0)
		return ""
	endif
	i0 = strsearch(buf,"?>",0)					// find end of header tag
	if (i0<0)
		return ""
	endif
	buf = buf[i0+2,Inf]
	return buf
End


Function WriteDetectorToFile(fileName,path,[ds])
	String fileName						// full path name to the file
	String path							// name of an Igor path to use
	STRUCT detectorGeometrys &ds		// structure defining a detector

	if (ParamIsDefault(ds))
		STRUCT detectorGeometrys dds
		if (FillDetectorsStruct(dds))		//fill the detector structure with test values
			DoAlert 0,"No current Detector Set\rNothing Done"
			return 1
		endif
	else
		CopyDetectorGeometrys(dds,ds)
	endif
	if (!(dds.N > 0))
		DoAlert 0,"No Valid Detectors\rNothing Done"
		return 1
	endif

	String xml = Detectors2xmlStr(dds)
	if (strlen(xml)<1)
		DoAlert 0,"Detector is Bad\rNothing Done"
		return 1
	endif

	Variable f
	if (strlen(fileName)<1)				// no file name passed
		fileName = SelectString(strlen(fileName),NewDetectorFileName(dds),fileName)
		Open/D/C="R*ch"/M="new detector geometry parameters file"/T=".xml" f as fileName
		fileName = S_fileName
	endif
	Open/C="R*ch"/M="new detector geometry parameters file"/P=$path/T=".xml"/Z f as fileName
	fileName = S_fileName
	if (V_flag)
		DoAlert 0, "nothing written to file"
		return 1
	endif
	FBinWrite f, xml
	Close f
	String str="wrote detector geometry to file '"+fileName+"'"
	print str
	DoAlert 0,str
	return 0
End
//
Static Function/T NewDetectorFileName(ds)
	STRUCT detectorGeometrys &ds
	Variable month,day,year,hour,minute,second,TZ, i, epoch=-1		// not using time zone
	String smonth
	for (i=0;i<MAX_DETECTORS;i+=1)
		if (ds.d[i].used)
			// sscanf d.timeMeasured, "%3s, %3s %d, %d, %02d:%02d:%02d (%g)",smonth,smonth,day,year,hour,minute,second,TZ
			sscanf ds.d[i].timeMeasured, "%3s, %3s %d, %d, %02d:%02d:%02d ",smonth,smonth,day,year,hour,minute,second
			if (V_flag==8)
				month = WhichListItem(smonth, "Jan;Feb;Mar;Apr;May;Jun;Jul;Aug;Sep;Oct;Nov;Dec")+1
				epoch = max(epoch, date2secs(year,month,day) + 3600*hour+60*minute+second)
			endif
			break
		endif
	endfor
	epoch = epoch<date2secs(2005,1,1) ? DateTime : epoch				// no valid times before 2005, use current time
	return "detector_"+Secs2Date(epoch,-2,"-")+"_"+ReplaceString(":",Secs2Time(epoch,3),"-")
End
//
Static Function/T Detectors2xmlStr(ds)
	STRUCT detectorGeometrys, &ds					// structure defining detectors

	Variable now = DateTime
	String str, xml="<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n\n<diffractometer xmlns=\"http://sector34.xray.aps.anl.gov/diffractometer\""
	if (strlen(ds.diffractometer))
		xml += " type=\""+ds.diffractometer+"\""
	endif
	xml += ">\n"
	xml += "\t<dateWritten>"+Secs2Date(now,2)+"</dateWritten>\n"
	sprintf str, "	<timeWritten>%s (%g)</timeWritten>\n", Secs2Time(now,3,1),date2secs(-1,-1,-1)/3600
	xml += str
	sprintf str, "	<EPOCH start=\"midnight Jan 1, 1904\" unit=\"sec\">%.0f</EPOCH>\n", now ; xml += str

	Variable last = ds.last
	if (last>=0 && last<MAX_DETECTORS)
		sprintf str, "\n	<Detectors Ndetectors=\"%d\" last=\"%d\">\t\t\t\t<!-- Ndetectors, number of detectors defined & last one used -->\n",ds.N,last
	else
		sprintf str, "\n	<Detectors Ndetectors=\"%d\">\t\t\t\t\t\t<!-- Ndetectors, number of detectors in defined -->\n",ds.N
	endif
	xml += str
	STRUCT detectorGeometry d
	Variable i, comment=1
	for (i=0;i<MAX_DETECTORS;i+=1,comment=0)
		if (!ds.d[i].used)
			continue
		endif
		CopyOneDetectorGeometry(d,ds.d[i])
		if (DetectorBad(d))
			continue
		endif
		sprintf str, "\t\t<Detector N=\"%d\">\n",i	;	xml += str
		sprintf str, "\t\t\t<Npixels>%d %d</Npixels>%s\n",d.Nx,d.Ny,SelectString(comment,"","\t\t\t\t<!-- Nx,Ny is number of un-binned pixels in full detector -->")
		xml += str
		sprintf str, "\t\t\t<size unit=\"mm\">%.7g %.7g</size>%s\n",d.sizeX,d.sizeY,SelectString(comment,"","\t\t<!-- sizeX,sizeY otuside size of full detector -->")
		xml += str
		sprintf str, "\t\t\t<R unit=\"radian\">%.9g %.9g %.9g</R>%s\n",d.R[0],d.R[1],d.R[2],SelectString(comment,"","\t<!-- Rotation and translation vector -->")
		xml += str
		sprintf str, "\t\t\t<P unit=\"mm\">%.3f %.3f %.3f</P>\n",d.P[0],d.P[1],d.P[2]
		xml += str
		sprintf str, "\t\t\t<pixel0>%.2f %.2f</pixel0>\t\t\t<!-- only provided for user reference -->\n",d.px0,d.py0
		xml += str
		if (strlen(d.timeMeasured))
			sprintf str, "\t\t\t<timeMeasured>%s</timeMeasured>%s\n",d.timeMeasured,SelectString(comment,"","	<!-- when this geometry was calculated -->")
			xml += str
		endif

		if (strlen(d.geoNote))
			sprintf str, "\t\t\t<note>%s</note>\n",d.geoNote
			xml += str
		endif
		if (strlen(d.detectorID))
			sprintf str, "\t\t\t<ID>%s</ID>%s\n",d.detectorID,SelectString(comment,"","\t<!-- unique detector ID -->")
			xml += str
		endif

		if (strlen(d.distortionMapFile))
			sprintf str, "\t\t\t<distortionMap>%s</distortionMap>%s\n",d.distortionMapFile,SelectString(comment,"","				<!-- file with distortion map -->")
			xml += str
		endif
		xml += "\t\t</Detector>\n"
	endfor
	xml +="\t</Detectors>\n"
	xml += "</diffractometer>\n"
	return xml
End



Function CopyDetectorGeometrys(f,i)				// copy a all detector structures
	STRUCT detectorGeometrys &f, &i				// f is the destination, i is source

	f.N = i.N
	f.last = i.last
	f.diffractometer = i.diffractometer
	Variable j
	for (j=0;j<MAX_DETECTORS;j+=1)
		CopyOneDetectorGeometry(f.d[j],i.d[j])
	endfor
End
//
Static Function CopyOneDetectorGeometry(f,i)		// copy a detector structure
	STRUCT detectorGeometry &f, &i				// f is the destination, i is source
	f.used = i.used
	f.Nx = i.Nx;			f.Ny = i.Ny
	f.sizeX = i.sizeX;		f.sizeY = i.sizeY
	f.R[0]=i.R[0];		f.R[1]=i.R[1];			f.R[2]=i.R[2];
	f.P[0]=i.P[0];		f.P[1]=i.P[1];			f.P[2]=i.P[2];
	f.timeMeasured = i.timeMeasured
	f.geoNote = i.geoNote
	f.detectorID = i.detectorID
	f.distortionMapFile = i.distortionMapFile
	f.rho00=i.rho00;		f.rho01=i.rho01;		f.rho02=i.rho02
	f.rho10=i.rho10;		f.rho11=i.rho11;		f.rho12=i.rho12
	f.rho20=i.rho20;		f.rho21=i.rho21;		f.rho22=i.rho22
	f.px0 = i.px0	;		f.py0 = i.py0
End

Static Function DetectorBad(d)
	STRUCT detectorGeometry &d
	if (!(d.used))
		return 1
	endif
	Variable bad = (numtype(d.Nx + d.Nx + d.sizeX + d.sizeY + d.R[0] + d.R[1] + d.R[2] + d.P[0] + d.P[1] + d.P[2])>0)
	bad += (d.Nx<1 || d.Nx>5000)												// detector cannot have more than 5000 pixels along one edge
	bad += (d.Ny<1 || d.Ny>5000)
	bad += (d.sizeX<1 || d.sizeX>2000)											// detector cannot be larger than 2m
	bad += (d.sizeY<1 || d.sizeY>2000)
	bad += (abs(d.R[0])>2*PI || abs(d.R[1])>2*PI || abs(d.R[2])>2*PI)		// rotation cannot be more than 2¹
	bad += (abs(d.P[0])>4000 || abs(d.P[0])>4000 || abs(d.P[0])>4000)		// P cannot be more than 4m in any direction
	bad += InValidPilatus(d)													// can only be bad if ID looks like a Pilatus 100K
	return (!(!bad))
End


Static Function InValidPilatus(d)							// checks if orientation of Pilatus is valid (beam goes in front)
	STRUCT detectorGeometry &d							// returns 0 if something set, 0 is nothing done
	if (strsearch(d.detectorID,"PILATUS 100K",0,2)<0)	// if the ID does not looks like a Pilatus 100K, then return isValid
		return 0
	endif

	Variable dpixel=0.172, Nx=487, Ny=195
	Variable err = sqrt((d.Nx - Nx)^2 + (d.Ny - Ny)^2 + (d.sizeX/dpixel - Nx)^2 + (d.sizeY/dpixel - Ny)^2)
	err = numtype(err) ? Inf : err
	if (err>0.01)											// check size to 0.01 pixels
		return 1
	endif

	Make/N=(3,3)/D/FREE rho
	rho[0][0] = d.rho00 ;	rho[0][1] = d.rho01 ;	rho[0][2] = d.rho02	// xyz = rho x [ (x' y' z') + P ]
	rho[1][0] = d.rho10 ;	rho[1][1] = d.rho11 ;	rho[1][2] = d.rho12	// rho is pre-calculated from vector d.R
	rho[2][0] = d.rho20 ;	rho[2][1] = d.rho21 ;	rho[2][2] = d.rho22

	Variable Naxes=NumVarOrDefault("root:Packages:Diffractometer:Naxes",DEFAULT_Naxes)
	Make/N=(Naxes)/FREE A=0
	Make/N=3/D/FREE xyz0,xhat,yhat

	Variable px0=d.px0, py0=d.py0
	pixel2xyz(d,px0,py0,A,xyz0)
	pixel2xyz(d,px0+10,py0,A,xhat)
	pixel2xyz(d,px0,py0+10,A,yhat)
	xhat = xhat - xyz0
	yhat = yhat - xyz0
	normalize(xhat)
	normalize(yhat)
	Cross xhat, yhat
	Wave W_cross
	Variable zcomp = W_cross[2]
	KillWaves/Z W_cross
	return abs(zcomp-1)<0.1 ? 0 : 1
End


Static Function SetDetectorParameters(Nx,Ny,dx,dy,Rstr,Pstr,dNote,Ptype,[quiet,id,fresh])	// set default detector to the reference values
	Variable Nx,Ny
	Variable dx,dy
	String Rstr, Pstr
	String dNote
	Variable Ptype								// 0="P[0],P[1],P[2]",  1="px,py,dist along Z"
	Variable quiet
	String id									// optional detector id
	Variable fresh								// a new detector (not changing a current one)
	quiet = ParamIsDefault(quiet) ? 0 : !quiet
	quiet = numtype(quiet) ? 0 : !(!quiet)		// default is NOT quiet
	fresh = ParamIsDefault(fresh) ? 0 : fresh
	fresh = numtype(fresh) ? 0 : !(!fresh)		// default is NOT a fresh one

	STRUCT detectorGeometrys ds
	FillDetectorsStruct(ds)
	id = SelectString(ParamIsDefault(id),id,"")
	Variable Nid=FindDetectorIdIndex(ds,id)		// index to this ID, -1 means not found
	if (Nid<0 && fresh)
		Nid = NextEmptyDetectorIndex(ds)
	endif
	Nid = Nid<0 ? ds.last : Nid
	if (fresh)
		Nid = Nid<0 ? MAX_DETECTORS-1 : Nid
	else
		Nid = Nid<0 ? 0 : Nid
	endif
	STRUCT detectorGeometry d
	CopyOneDetectorGeometry(d,ds.d[Nid])
	// set the default values

	Variable ask = ( !(Nx>0) || !(Ny>0) || !(dx>0) || !(dy>0) )
	Nx = Nx>0 ? Nx : d.Nx
	Ny = Ny>0 ? Ny : d.Ny
	dx = dx>0 ? dx : d.sizeX / d.Nx
	dy = dy>0 ? dy : d.sizeY / d.Ny
	dNote = SelectString(strlen(dNote),d.geoNote,dNote)

	Rstr = ReplaceString(",",Rstr," ")
	Rstr = ReplaceString(";",Rstr," ")
	Rstr = ReplaceString(":",Rstr," ")
	Variable xx,yy,zz
	sscanf Rstr,"%g %g %g",xx,yy,zz
	if (V_flag!=3)
		ask = 1
		sprintf Rstr, "%.13g, %.13g, %.13g", d.R[0], d.R[1], d.R[2]
	endif
	Pstr = ReplaceString(",",Pstr," ")
	Pstr = ReplaceString(";",Pstr," ")
	Pstr = ReplaceString(":",Pstr," ")
	sscanf Pstr,"%g %g %g",xx,yy,zz
	if (V_flag!=3)
		ask = 1
		sprintf Pstr, "%.13g, %.13g, %.13g", d.P[0], d.P[1], d.P[2]
	endif
	id = SelectString(strlen(id),d.detectorID,id)

	if (ask)
		Prompt Nx, "number of pixels in X"
		Prompt Ny, "number of pixels in Y"
		Prompt dx, "X size of pixels (mm)"
		Prompt dy, "Y size of pixels (mm)"
		Prompt Rstr, "Rotation vector (rad)"
		Prompt Pstr, "Translation vector (mm)"
		Prompt id,"Unique Detector ID"
		Prompt dNote, "Optional Note for Geometry"
		Prompt Ptype,"Translation vector description", popup,"P[0], P[1], P[2];px, py, P[2]"
		Ptype += 1
		DoPrompt "Detector Size & Location", Nx,Ny,dx,dy,Ptype,Pstr,Rstr,id,dNote
		if (V_flag)
			return 1
		endif
		Ptype -= 1
		printf  "SetDetectorParameters(%g,%g,%g,%g,\"%s\",\"%s\",\"%s\", %g)\r",Nx,Ny,dx,dy,Rstr,Pstr,dNote,Ptype
		quiet = 0
	endif

	STRUCT detectorGeometry df
	CopyOneDetectorGeometry(df,d)
	df.Nx = Nx ;				df.Ny = Ny						// get size and number of pixels
	df.sizeX = dx * Nx ;	df.sizeY = dy * Ny
	df.geoNote = dNote

	Rstr = ReplaceString(",",Rstr," ")						// get the R[3] vector
	Rstr = ReplaceString(";",Rstr," ")
	Rstr = ReplaceString(":",Rstr," ")
	sscanf Rstr,"%g %g %g",xx,yy,zz
	if (V_flag!=3)
		return 1
	endif
	df.R[0] = xx ;	df.R[1] = yy ;	df.R[2] = zz

	Pstr = ReplaceString(",",Pstr," ")						// get the P[3] vector
	Pstr = ReplaceString(";",Pstr," ")
	Pstr = ReplaceString(":",Pstr," ")
	sscanf Pstr,"%g %g %g",xx,yy,zz
	if (V_flag!=3)
		return 1
	endif
	if (Ptype==1)													// user entered center pixel & distance, not P[3]
		xx = (0.5*(Nx - 1) - xx) / (Nx/df.sizeX)
		yy = (0.5*(Ny-1) - yy) / (Ny/df.sizeY)
	endif
	df.P[0] = xx ;	df.P[1] = yy ;	df.P[2] = zz

	df.timeMeasured = date()+", "+time()
	df.detectorID = id
 	df.used = 1
	df.distortionMapFile = ""
	DetectorUpdateCalc(df)										// update all internally calculated things in the detector structure
	if (DetectorBad(df))
		DoAlert 0,"Detector Parameters are BAD\rNothing changed"
		return 1
	endif

	CopyOneDetectorGeometry(ds.d[Nid],df)
	ds.last = Nid
	String DiffractometerName = StrVarOrDefault("root:Packages:Diffractometer:DiffractometerName","")
	if (strlen(DiffractometerName))
		ds.diffractometer = DiffractometerName
	endif
	UpdateDefaultDetectorStruct(ds)							// does a DetectorUpdateCalc(), and stores new values where they will be found
	if (!quiet)
		print "changed detector from:"
		printOneDetector(d)
		print " "
		print "to:"
		printOneDetector(df)
	endif
	return 0
End


Static Function ResetSetPilatus100Kcalibration(xydir,px0,py0,dist)
	String xydir
	Variable px0,py0
	Variable dist

	String id = DEF_PILATUS_ID
	Variable Nx=487, Ny=195, dpixel=0.172			// number of un-binned pixels in detector
	String dirStrings="px -> +X & py -> +Y;px -> -Y & py -> +X;px -> -X & py -> -Y;px -> +Y & py -> -X;"

	STRUCT detectorGeometrys ds
	FillDetectorsStruct(ds)
	Variable Nid=FindDetectorIdIndex(ds,id)				// index to this ID, -1 means not found
	if (Nid<0)
		Nid = NextEmptyDetectorIndex(ds)
	endif
	Nid = Nid<0 ? ds.last : Nid
	Nid = Nid<0 ? MAX_DETECTORS-1 : Nid
	STRUCT detectorGeometry d
	CopyOneDetectorGeometry(d,ds.d[Nid])
	// set the default values

	px0 = px0 > 0 ? px0 : d.px0							// first default to previous
	py0 = py0 > 0 ? py0 : d.py0
	dist = dist > 0 ? dist : abs(d.P[2])

	px0 = px0 > 0 ? px0 : (Nx-1)/2					// reasonable default values, if previous are invalid
	py0 = py0 > 0 ? py0 : (Ny-1)/2
	dist = dist > 0 ? dist : 1000

	xydir = SelectString(WhichListItem(xydir,dirStrings,";",0,0)<0,xydir,"")// invalid xydir, set to ""
	if (strlen(xydir)<1)								// xydir not passed, set to current value
		Make/N=3/D/FREE xhat={d.rho00, d.rho10, d.rho20}
		normalize(xhat)
		if (abs(xhat[0]-1)<1e-3)
			xydir = "px -> +X & py -> +Y"
		elseif (abs(xhat[1]+1)<1e-3)
			xydir = "px -> -Y & py -> +X"
		elseif (abs(xhat[0]+1)<1e-3)
			xydir = "px -> -X & py -> -Y"
		elseif (abs(xhat[1]+1)<1e-3)
			xydir = "px -> -Y & py -> +X"
		else
			xydir = ""
		endif
	endif

	Prompt xydir,"X & Y pixel directions",popup,dirStrings
	Prompt px0, "X-pixel of beam at 2theta=0,  [0,"+num2istr(Nx-1)+"]"
	Prompt py0, "Y-pixel of beam at 2theta=0,  [0,"+num2istr(Ny-1)+"]"
	Prompt dist, "detector distance (mm)"
	DoPrompt "Detector Parameters",xydir,px0,py0,dist
	if (V_flag)
		return 1
	endif
	printf "ResetSetPilatus100Kcalibration(\"%s\", %g,%g, %g)\r",xydir,px0,py0,dist
	if (!(dist>0))
		print "Detector distance = %g, which is invalid"
		return 1
	elseif (numtype(px0+py0) || px0<0 || px0>=Nx || py0<0 || py0>=Ny)
		printf "Position of direct beam on detecter = (%g, %g), does not lie on the detector\r",px0,py0
		return 1
	endif

	Make/N=3/D/FREE 	xhat = NaN, yhat = NaN
	Variable i = WhichListItem(xydir,dirStrings,";",0,0)
	if (i==0)
		xhat = {1,0,0}
		yhat = {0,1,0}
	elseif (i==1)
		xhat = {0,-1,0}
		yhat = {1,0,0}
	elseif (i==2)
		xhat = {-1,0,0}
		yhat = {0,-1,0}
	elseif (i==3)
		xhat = {0,1,0}
		yhat = {-1,0,0}
	else
		return 1
	endif
	// printf "xhat = %s,   yhat = %s\r",vec2str(xhat),vec2str(yhat)
	Make/N=(3,3)/D/FREE rho=(p==q)
	rho[][0] = xhat[p]
	rho[][1] = yhat[p]
	if (MatrixDet(rho)<0)
		rho[2][2] = -1
	endif
	print "| rho | =",MatrixDet(rho)
	printWave(rho,name="",brief=1)

	Make/N=3/D/FREE Rvec
	Variable angle = axisOfMatrix(rho,Rvec)
	printf "rotation axis = %s,   angle = %g¡\r",vec2str(Rvec), angle
	Rvec *= angle*PI/180

	Variable/C pz = PvectorFromPixels(px0,py0)
	Make/N=3/D/FREE Pvec={real(pz),imag(pz),dist}
	MatrixOP/FREE vec = rho x Pvec
	Pvec[2] *= vec[2]<0 ? -1 : 1

	String dNote=d.geoNote, Rstr,Pstr
	Rstr = vec2str(Rvec,bare=1,places=12)
	Pstr = vec2str(Pvec,bare=1,places=12)

	print " "
	print "reset detector to:"
	printf "Using Nx=%g,  Ny=%g,  Æpixel=%g µm\r",Nx,Ny,dpixel
	printf "R[3] = \"%s\"\r",Rstr
	printf "P[3] = \"%s\",  pixel=[%g, %g]\r",Pstr,px0,py0

	Variable reset=1
	String str
	sprintf str,"Reset Detector to:  %s,  dist = %g mm",xydir,dist
	Prompt reset,str,popup,"NO;Reset Detector"
	Prompt dNote,"Note"
	DoPrompt "Reset Detector Parameters",reset,dNote
	if (V_flag)
		return 1
	endif
	print " "
	if (reset==2)
		SetDetectorParameters(Nx,Ny,dpixel,dpixel,Rstr,Pstr,dNote,0,quiet=0,id=id)	// set default detector to the reference values
	else
		print "Detector Calibration UN-changed"
	endif
	return 0
End


Static Function/C PvectorFromPixels(px0,py0)
	Variable px0,py0

	STRUCT detectorGeometry d
	FillDetectorStructDefault(d,"")
	if (!(px0>=0 && px0<(d.Nx) && py0>=0 && py0<(d.Ny)))
		Prompt px0,"X-pixel of incident beam"
		Prompt py0,"Y-pixel of incident beam"
		DoPrompt "Pixel of Incident Beam",px0,py0
		if (V_flag)
			return cmplx(NaN,NaN)
		endif
	endif
	if (!(px0>=0 && px0<(d.Nx) && py0>=0 && py0<(d.Ny)))
		return cmplx(NaN,NaN)
	endif

	Variable P0 = (0.5*(d.Nx-1) - px0) * (d.sizeX/d.Nx)
	Variable P1 = (0.5*(d.Ny-1) - py0) * (d.sizeY/d.Ny)
	if (strlen(GetRTStackInfo(2))<1)
		printf "new P =   \"%.10g,  %.10g,  %.10g\"\r",P0,P1,d.P[2]
	endif
	return cmplx(P0,P1)
End

//  ============================================================================  //
//  ============================ End of Detector Orientation ============================  //




//  ============================================================================  //
//  ============================= Start of Utility Routines =============================  //

//ThreadSafe Static Function/WAVE recip_from_xtal(xtal)
Static Function/WAVE recip_from_xtal(xtal)
	STRUCT crystalStructure &xtal
	Make/N=(3,3)/D/FREE recip
	recip[0][0] = {xtal.as0, xtal.as1, xtal.as2}
	recip[0][1] = {xtal.bs0, xtal.bs1, xtal.bs2}
	recip[0][2] = {xtal.cs0, xtal.cs1, xtal.cs2}
	return recip
End


// set mat to be a rotation matrix about axis with angle
//ThreadSafe Static Function/WAVE rotationMatFromAxis(axis,angle)
Static Function/WAVE rotationMatFromAxis(axis,angle)
	Wave axis				// axis about which to rotate (or possibly Rodriques vector)
	Variable angle			// angle to rotate (degrees), assumes axis is true Rotation vector if angle invalid
	Wave mat				// desired rotation matrix

	if (!WaveExists(axis) || numtype(angle))
		return $""
	endif

	Make/N=(3,3)/FREE mat
	Variable len = norm(axis)
	angle = numtype(angle) ? len : angle*PI/180	// the rotation angle (rad)
	if (angle==0)									// zero angle rotation is just the identity matrix
		mat = (p==q)
		return mat
	endif

	Variable nx=axis[0]/len, ny=axis[1]/len, nz=axis[2]/len
	Variable cosa=cos(angle), sina=sin(angle)
	Variable c1 = 1-cosa
	// from		http://mathworld.wolfram.com/RodriguesRotationFormula.html (I double checked this too.)
	mat[0][0] = nx*nx*c1 + cosa;			mat[0][1] = nx*ny*c1 - nz*sina;		mat[0][2] = nx*nz*c1 + ny*sina
	mat[1][0] = nx*ny*c1 + nz*sina;		mat[1][1] = ny*ny*c1 + cosa;			mat[1][2] = ny*nz*c1 - nx*sina
	mat[2][0] = nx*nz*c1 - ny*sina;		mat[2][1] = ny*nz*c1 + nx*sina;		mat[2][2] = nz*nz*c1 + cosa
	return mat
End



ThreadSafe Static Function/WAVE wCross(a,b)
	Wave a,b

	if (!WaveExists(a) || !WaveExists(b))
		return $""							// both waves must exist
	elseif (DimSize(a,0)!=3 || WaveDims(a)!=1)
		return $""							// a must be 3 long
	elseif (DimSize(b,0)!=3 || WaveDims(b)!=1)
		return $""							// b must be 3 long
	endif
	Make/N=3/D/FREE c
	c[0] = a[1]*b[2] - a[2]*b[1]
	c[1] = a[2]*b[0] - a[0]*b[2]
	c[2] = a[0]*b[1] - a[1]*b[0]
	return c
End


ThreadSafe Static Function angleBetweenVecs(a,b)
	Wave a,b
	Variable dot = MatrixDot(a,b)/(norm(a)*norm(b))
	return acos(dot) * 180/PI
End

//  ============================================================================  //
//  ============================== End of Utility Routines =============================  //



//  ============================================================================  //
//  ================================= Start of Init ==================================  //


Static Function SelectDiffractometer(name)
	String name

	SVAR/Z DiffractometerName = root:Packages:Diffractometer:DiffractometerName
	SVAR/Z DiffractometerAxisNames = root:Packages:Diffractometer:DiffractometerAxisNames
	NVAR/Z Naxes = root:Packages:Diffractometer:Naxes
	NVAR/Z KA = root:Packages:Diffractometer:KA
	if (!NVAR_Exists(Naxes) || !SVAR_Exists(DiffractometerName))
		Init_Diffractometer()
		SVAR/Z DiffractometerName = root:Packages:Diffractometer:DiffractometerName
		SVAR/Z DiffractometerAxisNames = root:Packages:Diffractometer:DiffractometerAxisNames
		NVAR/Z Naxes = root:Packages:Diffractometer:Naxes
		NVAR/Z KA = root:Packages:Diffractometer:KA
	endif
	String listFull = "fourc:33BM,4,2-theta theta chi phi;"				// "name0:where0,Naxes0,axisList0;name1:where1,Naxes1,axisList1; ..."
	listFull += "kappa:33ID,6,delta theta kappa phi mu nu;"
	listFull += "CDP:33ID,3,nu mu delta;"
	if (!keyInList(DiffractometerName,listFull,":",";"))		// current diffractometer is not in fullList, add it
		listFull += DiffractometerName+":"+num2istr(Naxes)+","+ReplaceString(";", DiffractometerAxisNames,",")+";"
	endif

	if (!keyInList(name,listFull,"",""))
		Variable m, i=NaN
		String key,value,item,listPopup=""
		for (i=0;i<ItemsInList(listFull);i+=1)		// build popup string & and set default to current diffractometer name
			item = StringFromList(i,listFull)
			key = StringFromList(0,item,":")
			value = StringByKey(key,item)
			item = StringFromList(0,item,":")+",  used on "+StringFromList(0,value,",")
			listPopup += item+";"
			name = SelectString(stringmatch(key,DiffractometerName),name,item)
		endfor

		item = name
		Prompt item,"Diffractometer Type",popup,listPopup
		DoPrompt "Diffractometer", item
		if (V_flag)
			return 1
		endif
		name = StringFromList(0,item,",")
		printf "SelectDiffractometer(\"%s\")\r",name
	endif
	value = StringByKey(name,listFull)
	Variable N = str2num(StringFromList(1,value,","))
	String axes = ReplaceString(" ",StringFromList(2,value,","),";")

	if (!keyInList(name,listFull,"","") && !(N>0))	// check that name is valid
		return 1
	endif
	DiffractometerName = name						// update global name
	Naxes = N										// update global Naxes
	DiffractometerAxisNames = axes					// update global axis names
	KA = (strsearch(name,"kappa",0,2)<0) ? NaN : 50

	STRUCT diffractometerTypeStruct dt
	dt.name = name									// set structure to new global values
	dt.Naxes = N
	dt.axesNames = axes
	dt.KA = KA
	SavePackagePreferences/FLSH=1 "Diffractometer","diffractometerTypePrefs",0,dt	// update prefs
	return 0
End


Function Init_Diffractometer()
	InitLatticeSymPackage()

	NewDataFolder/O root:Packages
	NewDataFolder/O root:Packages:Diffractometer

	Variable setPackageDefaults=0
	STRUCT diffractometerTypeStruct dt
	LoadPackagePreferences/MIS=0 "Diffractometer","diffractometerTypePrefs",0,dt
	if (V_flag || V_bytesRead<V_structSize || strlen(dt.name)<1)	// could not load diffractometer type, set to default
		setPackageDefaults = 1
		dt.name = DEFAULT_DIFF_NAME
		dt.Naxes = DEFAULT_Naxes
		dt.axesNames = DEFAULT_AXIS_NAMES
		dt.KA = 50
	endif

	if (Exists("root:Packages:Diffractometer:Naxes")!=2)
		Variable/G root:Packages:Diffractometer:Naxes=dt.Naxes
	endif
	if (Exists("root:Packages:Diffractometer:DiffractometerName")!=2)
		String/G root:Packages:Diffractometer:DiffractometerName=dt.name
	endif
	if (Exists("root:Packages:Diffractometer:DiffractometerAxisNames")!=2)
		String/G root:Packages:Diffractometer:DiffractometerAxisNames=dt.axesNames
	endif
	if (Exists("root:Packages:Diffractometer:KA")!=2)
		Variable/G root:Packages:Diffractometer:KA=dt.KA	// kappa angle
	endif
	if (setPackageDefaults)
		dt.name = StrVarOrDefault("root:Packages:Diffractometer:DiffractometerName",DEFAULT_DIFF_NAME)
		dt.Naxes = NumVarOrDefault("root:Packages:Diffractometer:Naxes",DEFAULT_Naxes)
		dt.axesNames = StrVarOrDefault("root:Packages:Diffractometer:DiffractometerAxisNames",DEFAULT_AXIS_NAMES)
		dt.KA = NumVarOrDefault("root:Packages:Diffractometer:KA",50)
		SavePackagePreferences/FLSH=1 "Diffractometer","diffractometerTypePrefs",0,dt	// set package prefs
	endif

	if (strlen(FunctionList("FillDetectorsStruct","","KIND:2")))
		STRUCT detectorGeometrys ds
		if (FillDetectorsStruct(ds))					//fill the detector structure with current values
			SetDefaultDetector2Reference(0)			// there was an error, set to idealized reference values
		endif
	endif

	if (strlen(FunctionList("FillSampleStructDefault","","KIND:2")))
		STRUCT sampleStructure s
		if (FillSampleStructDefault(s))					//fill the sample structure with current values
			SetDefaultSample2Reference()				// there was an error, set to idealized reference values
		endif
	endif
End

//  ============================================================================  //
//  ================================== End of Init ==================================  //



//	#G1
//	1.54			def g_aa        'U[0]'  # a lattice constant (real space)
//	1.54			def g_bb        'U[1]'  # b lattice constant (real space)
//	1.54			def g_cc        'U[2]'  # c lattice constant (real space)
//	90				def g_al        'U[3]'  # alpha lattice angle (real space)
//	90				def g_be        'U[4]'  # beta  lattice angle (real space)
//	90				def g_ga        'U[5]'  # gamma lattice angle (real space)
//	
//	4.079990459		def g_aa_s      'U[6]'  # a lattice constant (reciprocal space)
//	4.079990459		def g_bb_s      'U[7]'  # b lattice constant (reciprocal space)
//	4.079990459		def g_cc_s      'U[8]'  # c lattice constant (reciprocal space)
//	90				def g_al_s      'U[9]'  # alpha lattice angle (reciprocal space)
//	90				def g_be_s      'U[10]' # beta  lattice angle (reciprocal space)
//	90				def g_ga_s      'U[11]' # gamma lattice angle (reciprocal space)
//	
//	1				def g_h0        'U[12]' # H of primary reflection
//	0				def g_k0        'U[13]' # K of primary reflection
//	0				def g_l0        'U[14]' # L of primary reflection
//	
//	0				def g_h1        'U[15]' # H of secondary reflection
//	1				def g_k1        'U[16]' # K of secondary reflection
//	0				def g_l1        'U[17]' # L of secondary reflection
//	
//	60				def g_u00       'U[18]' # tth0, angles of primary reflection
//	30				def g_u01       'U[19]'	th0
//	0				def g_u02       'U[20]'	om0
//	0				def g_u03       'U[21]'	chi0
//	0				def g_u04       'U[22]'	phi0
//	0				def g_u05       'U[23]'
//	
//	60				def g_u10       'U[24]' # tth1, angles of secondary reflection
//	30				def g_u11       'U[25]'	th1
//	0				def g_u12       'U[26]'	om1
//	-90				def g_u13       'U[27]'	chi1
//	0				def g_u14       'U[28]'	phi1
//	0				def g_u15       'U[29]'
//	
//	1.54			def g_lambda0   'U[30]' # lambda when or0 was set
//	1.54			def g_lambda1   'U[31]' # lambda when or1 was set

