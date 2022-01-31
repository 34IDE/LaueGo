#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=1		// Use modern global access method.
#pragma ModuleName=Indexing
#pragma IgorVersion = 6.2
#pragma version = 5.05
#include "LatticeSym", version>=6.28
#include "microGeometryN", version>=1.98
#include "Masking", version>1.04
#include "ImageDisplayScaling", version>=2.16
#if (NumVarOrDefault("root:Packages:MICRO_GEOMETRY_VERSION",0)&2)
#include "tiff"
//#else
#endif
//#if (Exists("HDF5OpenFile")==4)
#if (NumVarOrDefault("root:Packages:MICRO_GEOMETRY_VERSION",0)&1 && Exists("HDF5OpenFile")==4)
#include "HDF5images", version>=0.35
#endif
#if (NumVarOrDefault("root:Packages:MICRO_GEOMETRY_VERSION",0)&4)
#include "WinView", version>=2.06
//#else
#endif
Constant INDEXING_MAX_CALC = 30
Constant INDEXING_MAX_TEST = 51
Constant INDEXING_ASK_BEFORE_DISPLAYING = 0		// This can be overridden in your Main Procedure window

//	with version 2.50, I changed so that yc in the file header always superceedes yc in the geo panel
//	with version 3.00, allow use of both microGeo and microGeoN
//	with version 3.40, change strain refinement to use energy, it also now works in k-spcace, not angle space
//	with version 4.00, add the ability to use both spe and hdf5 with the new geometry
//	with version 4.18, changed default name of FullPeakList to something like FullPeakListOrange
//	with version 4.60, changed how pick_keV_calc() works
//	with version 4.88, changed number of columns in FullPeakList from 11 to 12
//	with version 4.98, added support of showing hkl on an image from a user supplied reciprocal lattice
//	with version 4.99, improved CloseAllowedhkl(), more UTF8 and less symbol font

//	#define USE_ENERGY_STRAIN_REFINE		// paste this line (uncommented) in your Procedure Window to use newer strain refinement

#if (IgorVersion()>7)
	// Igor 7 use UNICODE
	Static strConstant S_DEGREE = "\xC2\xB0"						// UTF8, degree sign
	Static strConstant S_THETA = "\xCE\xB8"							// UTF8, lower case theta
	Static strConstant S_DELTA_THETA = "\xCE\x94\xCE\xB8"		// UTF8, lower case theta
#else
	Static strConstant S_DEGREE = "\\F'Symbol'\260\\F]0"		// degree sign in Symbol font
	Static strConstant S_THETA = "\\F'Symbol'q\\F]0"				// theta in Symbol font
	Static strConstant S_DELTA_THETA = "\\F'Symbol'Dq\\F]0"	// Delta-theta in Symbol font
#endif


Menu LaueGoMainMenuName
	SubMenu "Crystal Lattice"
		"Set Crystal Lattice...", /Q, MakeLatticeParametersPanel("")
		"Show Crystal Lattice",/Q,showLattice()
		help={"Shows the crystal structure and lattice that is currently defined"}
		"d[hkl]",/Q,get_dhkl(NaN,NaN,NaN)
		help={"show d-spacing for the hkl reflection"}
	End
	"Load image File...", LoadGenericImageFile("")
	MenuItemIfWaveClassExists("Save Mask File...","imageMask","DIMS:2"),GenericSaveImage("",$"",extras="classes:imageMask")
	help={"Load an Image from the file, and then display the image.  If you just want to display a loaded file, go to 'HDF or WinView' menu"}
	MenuItemIfWaveClassExists("Remove Background from Image...","speImage;rawImage*","DIMS:2"),RemoveBackgroundImage($"",NaN)
	help={"remove background from an image, do not use on depth sotred images"}
	MenuItemIfWaveClassExists("Fit Peaks on Image...","speImage*;rawImage*","DIMS:2"),FitPeaks($"",NaN,NaN)
	help={"Identify peaks in an Image, does not do any background removal"}
	MenuItemIfWaveClassExists("   TESTING, Fit Peaks step-wise on Image...","speImageNoBkg;speImage;rawImageNoBkg;rawImage","DIMS:2"), FitPeaksStepWise($"",NaN,NaN,NaN,NaN)
	help={"Identify peaks in an Image, does not do any background removal, uses a simple threshold to find peaks"}
	MenuItemIfWaveClassExists("   TESTING, Fit Peaks on Image...","speImageNoBkg;speImage;rawImageNoBkg;rawImage","DIMS:2"),FitPeaksNew($"",NaN,NaN)
//	help={"Identify peaks in an Image, does not do any background removal, uses a simple threshold to find peaks"}
//	     MenuItemIfWaveClassExists("   OLD Remove Bkg & Fit Peaks on Image...","speImage;rawImage","DIMS:2"),OLD_removeBkg_FitPeaks($"",NaN,NaN)
	help={"removes bkg and then finds and fits peaks, do not use anymore"}
//	  pkLIstImageMenuItem("   Reset Fitted Peak list from boxes"),setFittedPeaksFromList($"",NaN,pkList)
//	  help={"reset list of peaks for fitting to the boxes on an image plot"}
	MenuItemIfWaveClassExists("Index and Display...","FittedPeakList*","DIMS:2,MAXCOLS:12,MINCOLS:11"),IndexAndDisplay($"",NaN,NaN,NaN,NaN,NaN,NaN,NaN)
	help={"Index the identified peaks, using the defined lattice, and then show it all on a plot"}
	   MenuItemIfWaveClassExists("  Refine Strain","IndexedPeakList*",""),Indexing#MenuStrainRefine(NaN,"")
//	   MenuItemIfWaveClassExists("  Refine Strain","IndexedPeakList*",""),DeviatoricStrainRefine(NaN,"",printit=1)
//	   MenuItemIfWaveClassExists("  Refine Strain","IndexedPeakList*",""),TotalStrainRefine(NaN,"")
	   help={"compute deviatoric strain for the pattern just idexed"}
	   MenuItemIfWaveClassExists("  Re-Draw hkl tags","IndexedPeakList*",""),DisplayResultOfIndexing($"",NaN)
	   help={"for an image plot, clear old hkl tags (if any) and redraw new ones using results of an indexing"}
	   MenuItemIfWaveClassExists("  Clear hkl tags off Graph","IndexedPeakList*",""),ClearhklOffGraph("")
	   help={"clear the hkl tags off of an image"}
	   MenuItemIfWaveClassExists(SelectString(exists("MakeStereo")==6 , "  Load Stereographic","  Stereographic Projection..."),"IndexedPeakList*",""),StereoOfIndexedPattern($"",NaN)
	   help={"Make a Stereographic projection with the same orientation as the indexed pattern"}
	   MenuItemIfWaveClassExists("  Report Sheet of Indexed Image","IndexedPeakList*",""),MakeIndexingReportSheet()
	   help={"if the top window is an images showing the hkl's, this makes a layout with reflections for printing"}
	   MenuItemIfWaveClassExists("  calc energy of [hkl]","IndexedPeakList*",""),EnergyOfhkl($"",NaN,NaN,NaN,NaN)
	   help={"calculate the energy of an hkl using the current indexing"}
		MenuItemIfWaveClassExists("  predict cone angle for detector","speImage;rawImage*","DIMS:2"),EstimateConeAngle(NaN)
	   help={"if you know the hkl at the center of the detector, this calculates the cone angle to fill the detector"}
	   MenuItemIfWaveClassExists("  Rotation between two indexations...","IndexedPeakList*",""),RotationBetweenTwoIndexations($"",NaN,$"",NaN)
	   help={"calculate the rotation between two different indexation results, useful when indexing finds multiple patterns"}
	//	SubMenu("Tables")
	MenuItemsWaveClassOnGraph("Angle Between Cursors","speImage*;rawImage*",""),/Q,AngleBetweenTwoPointsGeneral($"",$"",nan, nan,nan,nan)
	SubMenu(MenuItemIfWaveClassExists("Tables","FittedPeakList*;IndexedPeakList*",""))
		MenuItemIfWaveClassExists("Fitted Peaks","FittedPeakList*",""),DisplayTableOfWave($"",classes="FittedPeakList*",promptStr="Wave with list of fitted peaks")
		help={"Show table of identified peak positions"}
		MenuItemIfWaveClassExists("Indexed Peaks","IndexedPeakList*",""),DisplayTableOfWave($"",classes="IndexedPeakList*",promptStr="Wave with list of indexed peaks")
		help={"Show table of indexed peaks"}
		MenuItemIfWaveClassExists("Strain Peaks","IndexedPeakList*;StrainPeakList",""),DisplayTableOfWave($"",classes="StrainPeakList*",promptStr="Wave with list of Strain Refined peaks")
		help={"Show table of peaks from the strain refinement"}
		"New Wave for measured Energies...",MakeMeasured_hkl_EnergiesWave(NaN,"")
	End
	"  Choose Custom Indexing Program...",pickIndexingFunction("")
End
// MenuStrainRefine() is only temporary until I decide to delete the function DeviatoricStrainRefine().  We should only be using TotalStrainRefine().
Static Function MenuStrainRefine(num,str)
	Variable num
	String str
#ifdef USE_ENERGY_STRAIN_REFINE
	TotalStrainRefine(num,str)
#else
	DoAlert 0,"You are using the Deviatoric strain refinement, a Total refinement is available, see\r'#define USE_ENERGY_STRAIN_REFINE' in IndexingN.ipf"
	DeviatoricStrainRefine(num,str,printit=1)
#endif
End
Menu "GraphMarquee", dynamic
	"-"
	MenuItemIndexPeak("Indexed Peak Info"),/Q,IndexedPeakInfoFromMarquee()
End
//
Function/S MenuItemIndexPeak(item)
	String item
	GetMarquee/Z
	if (!V_flag)
		return "("+item
	endif
	if (strlen(GetUserData("","","Indexing"))<1)
		return "("+item
	endif
	return item
End



Static Constant hc = 1.239841857			// keV-nm
Static Constant SMALLEST1 = 1.2e-16


Function/WAVE IndexAndDisplay(FullPeakList0,keVmaxCalc,keVmaxTest,angleTolerance,hp,kp,lp,cone,[maxSpots,FullPeakList1,FullPeakList2,depth,printIt])
	Wave FullPeakList0				// contains the result of a peak fitting
	Variable keVmaxCalc			// 17, maximum energy to calculate (keV)
	Variable keVmaxTest			// 30, maximum energy to test (keV)  [-t]
	Variable angleTolerance		// 0.25, angular tolerance (deg)
	Variable hp,kp,lp				// preferred hkl
	Variable cone					// angle from preferred hkl, (0 < cone ≤ 180°)
	Variable maxSpots				// -n max num. of spots from data file to use, default is 250
	Wave FullPeakList1				// contains the result of a peak fitting
	Wave FullPeakList2				// contains the result of a peak fitting
	Variable depth					// optional depth (normally get depth from FullPeakList0 wave note)
	Variable printIt				// forces print out
	maxSpots = ParamIsDefault(maxSpots) ? -1 : maxSpots
	maxSpots = ((maxSpots>2) && numtype(maxSpots)==0) ? maxSpots : -1
	if (ParamIsDefault(FullPeakList1))
		Wave FullPeakList1=$""
	endif
	if (ParamIsDefault(FullPeakList2))
		Wave FullPeakList2=$""
	endif
	depth = ParamIsDefault(depth) || numtype(depth) ? NaN : depth
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? (strlen(GetRTStackInfo(2))==0) : !(!printIt)

	if (NumVarOrDefault("root:Packages:geometry:PanelValues:dirty",0))	// geo waiting
		Variable/G root:Packages:geometry:PanelValues:dirty=0
		DoAlert 0,"You have made changes in the geometry panel, deal with them first"
		MakeGeometryParametersPanel("")
		return $""
	endif
	Variable is4500S=NumVarOrDefault(MICRO_GEOMETRY_VERSION_PATH,0)&2
	Variable badWave = !WaveExists(FullPeakList0) && !WaveExists(FullPeakList1) && !WaveExists(FullPeakList2)
	if (!badWave)
		if (WaveExists(FullPeakList0))
			badWave = badWave && (DimSize(FullPeakList0,0)<1 || DimSize(FullPeakList0,1)<11)
		endif
		if (WaveExists(FullPeakList1))
			badWave = badWave && (DimSize(FullPeakList1,0)<1 || DimSize(FullPeakList1,1)<11)
		endif
		if (WaveExists(FullPeakList2))
			badWave = badWave && (DimSize(FullPeakList2,0)<1 || DimSize(FullPeakList2,1)<11)
		endif
	endif
	Variable badNums= !(keVmaxCalc>1 && keVmaxCalc<INDEXING_MAX_CALC) || !(keVmaxTest>1 && keVmaxTest<INDEXING_MAX_TEST)
	badNums += !(angleTolerance>=0.01 && angleTolerance<10)
	badNums += numtype(hp+kp+lp)
	badNums += !(cone>1 && cone<=180)
	badNums += (abs(hp)+abs(kp)+abs(lp))==0
	if (badWave || badNums)
		keVmaxCalc = (keVmaxCalc>1 && keVmaxCalc<INDEXING_MAX_CALC) ? keVmaxCalc : NaN
		keVmaxCalc = numtype(keVmaxCalc) ? NumVarOrDefault("root:Packages:micro:Index:keVmaxCalc",NaN) : keVmaxCalc
		keVmaxCalc = numtype(keVmaxCalc) ? pick_keV_calc() : keVmaxCalc
		keVmaxTest = (keVmaxTest>1 && keVmaxTest<INDEXING_MAX_TEST) ? keVmaxTest : NaN
		keVmaxTest = numtype(keVmaxTest) ? NumVarOrDefault("root:Packages:micro:Index:keVmaxTest",NaN) : keVmaxTest
		keVmaxTest = numtype(keVmaxTest) ? min(3*keVmaxCalc,30) : keVmaxTest
		angleTolerance = (angleTolerance>=0.01 && angleTolerance<10) ? angleTolerance : NumVarOrDefault("root:Packages:micro:Index:angleTolerance",(is4500S ? 0.5 : 0.1))
		hp = numtype(hp) ? NumVarOrDefault("root:Packages:micro:Index:hp",0) : hp
		kp = numtype(kp) ? NumVarOrDefault("root:Packages:micro:Index:kp",0) : kp
		lp = numtype(lp) ? NumVarOrDefault("root:Packages:micro:Index:lp",2) : lp
		cone = (cone>1 && cone<=180) ? cone : NumVarOrDefault("root:Packages:micro:Index:cone",( is4500S ? 179 : 72 ))
		String hkl
		sprintf hkl,"%d, %d, %d",hp,kp,lp
		Prompt hkl,"preferred hkl of center"
		Prompt cone,"max angle"+DEGREESIGN+" from hkl prefer to spot"
		String pkList=WaveListClass("FittedPeakList*","*","DIMS:2,MAXCOLS:12,MINCOLS:11")
		String peakListStr0="", peakListStr1="", peakListStr2=""
		Variable multi=0
		if (WaveExists(FullPeakList0))
			peakListStr0 = NameOfWave(FullPeakList0)
			multi += 1
		else
			peakListStr0 = StringFromList(0,pkList)
		endif
		if (WaveExists(FullPeakList1))
			peakListStr1 = NameOfWave(FullPeakList1)
			multi += 1
		else
			peakListStr1 = StringFromList(1,pkList)
		endif
		if (WaveExists(FullPeakList2))
			peakListStr2 = NameOfWave(FullPeakList2)
			multi += 1
		else
			peakListStr2 = StringFromList(2,pkList)
		endif
		multi = multi>1 ? 1 : (ItemsInList(pkList)>1)
		pkList = reverseList(pkList)
		Prompt peakListStr0,"wave with fitted peaks",popup,pkList
		Prompt peakListStr1,"wave1 with fitted peaks",popup,"-none-;"+pkList
		Prompt peakListStr2,"wave2 with fitted peaks",popup,"-none-;"+pkList
		Prompt keVmaxCalc,"max energy for searching (keV)"
		Prompt keVmaxTest,"max energy for matching (keV)"
		Prompt angleTolerance "max angle"+DEGREESIGN+" between peak & hkl"
		Prompt depth "depth ("+Gmu+"m) (NaN for default)"
		if (multi)
			DoPrompt/Help="3D-Xray Diffraction[Indexing]" "index",peakListStr0,keVmaxCalc,peakListStr1,hkl,peakListStr2,keVmaxTest,cone,angleTolerance, depth
		else
			DoPrompt/Help="3D-Xray Diffraction[Indexing]" "index",peakListStr0,keVmaxCalc,hkl,keVmaxTest,cone,angleTolerance,depth
		endif
		if (V_flag)
			return $""
		endif
		peakListStr0 = SelectString(strsearch(peakListStr0,"-none-",0,2)>=0,peakListStr0,"")
		peakListStr1 = SelectString(strsearch(peakListStr1,"-none-",0,2)>=0,peakListStr1,"")
		peakListStr2 = SelectString(strsearch(peakListStr2,"-none-",0,2)>=0,peakListStr2,"")
		Wave FullPeakList0=$peakListStr0, FullPeakList1=$peakListStr1, FullPeakList2=$peakListStr2
		badWave = !WaveExists(FullPeakList0) && !WaveExists(FullPeakList1) && !WaveExists(FullPeakList2)
		Variable N0=0,N1=0,N2=0
		if (!badWave)
			if (WaveExists(FullPeakList0))
				badWave = badWave && (DimSize(FullPeakList0,0)<1 || DimSize(FullPeakList0,1)<11)
				N0 = DimSize(FullPeakList0,0)
			endif
			if (WaveExists(FullPeakList1))
				badWave = badWave && (DimSize(FullPeakList1,0)<1 || DimSize(FullPeakList1,1)<11)
				N1 = DimSize(FullPeakList1,0)
			endif
			if (WaveExists(FullPeakList2))
				badWave = badWave && (DimSize(FullPeakList2,0)<1 || DimSize(FullPeakList2,1)<11)
				N2 = DimSize(FullPeakList2,0)
			endif
		endif
		if (!(cone>1 && cone<=180) && WaveExists(FullPeakList0))
			cone = EstimateConeAngle(NaN,image=FullPeakList0)	// user requested auto-guess of cone angle
		endif
		if (ParamIsDefault(maxSpots) && (N0+N1+N2)>250)		// more than 250 spots in FullPeakLIst
			Prompt maxSpots, "max number to try to index, out of "+num2istr(N0+N1+N2)+", -1 uses default of 250"
			DoPrompt "Number of Spots", maxSpots
			if (V_flag)
				return $""
			endif
			maxSpots = ((maxSpots>2) && numtype(maxSpots)==0) ? maxSpots : NumVarOrDefault("root:Packages:micro:Index:maxSpots",-1)
		endif
		badNums= !(keVmaxCalc>1 && keVmaxCalc<INDEXING_MAX_CALC) || !(keVmaxTest>1 && keVmaxTest<INDEXING_MAX_TEST)
		badNums += !(angleTolerance>=0.01 && angleTolerance<10)
		badNums += numtype(hp+kp+lp) + !(cone>1 && cone<=180)
		badNums += (abs(hp)+abs(kp)+abs(lp))==0
		sscanf hkl,"%d, %d, %d", hp,kp,lp
		printf "•IndexAndDisplay(%s,%g,%g,%g, %d,%d,%d,%g",peakListStr0,keVmaxCalc,keVmaxTest,angleTolerance,hp,kp,lp,cone
		if (maxSpots>0)
			printf ",maxSpots=%g",maxSpots
		endif
		if (strlen(peakListStr1))
			printf ", FullPeakList1=%s",peakListStr1
		endif
		if (strlen(peakListStr2))
			printf ", FullPeakList2=%s",peakListStr2
		endif
		if (numtype(depth)==0)
			printf ", depth=%g",depth
		endif
//		if (!ParamIsDefault(printIt))
//			printf ", printIt=%g",printIt
//		endif
		printf ")\r"
		if (badWave || badNums)
			DoAlert 0, "Invalid inputs sent to IndexAndDisplay()"
			return $""
		endif
	endif
	NewDataFolder/O root:Packages
	NewDataFolder/O root:Packages:micro
	NewDataFolder/O root:Packages:micro:Index
	Variable/G root:Packages:micro:Index:angleTolerance=angleTolerance, root:Packages:micro:Index:cone=cone
	Variable/G root:Packages:micro:Index:hp=hp, root:Packages:micro:Index:kp=kp, root:Packages:micro:Index:lp=lp
	Variable/G root:Packages:micro:Index:maxSpots=maxSpots
	Variable/G root:Packages:micro:Index:keVmaxCalc=keVmaxCalc, root:Packages:micro:Index:keVmaxTest=keVmaxTest
	String/G root:Packages:micro:Index:FullPeakList0=NameOfWave(FullPeakList0), root:Packages:micro:Index:FullPeakList1=NameOfWave(FullPeakList1), root:Packages:micro:Index:FullPeakList2=NameOfWave(FullPeakList2)

	depth = numtype(depth) ? NumberByKey("depth",note(FullPeakList0),"=") : depth
	String args=""
	args = ReplaceStringByKey("FullPeakList",args, GetWavesDataFolder(FullPeakList0,2))
	if (WaveExists(FullPeakList1))
		args = ReplaceStringByKey("FullPeakList1",args, GetWavesDataFolder(FullPeakList1,2))
	endif
	if (WaveExists(FullPeakList2))
		args = ReplaceStringByKey("FullPeakList2",args, GetWavesDataFolder(FullPeakList2,2))
	endif
	args = ReplaceNumberByKey("keVmaxCalc",args,keVmaxCalc)
	args = ReplaceNumberByKey("keVmaxTest",args,keVmaxTest)
	args = ReplaceNumberByKey("angleTolerance",args,angleTolerance)
	args = ReplaceNumberByKey("hp",args,hp)
	args = ReplaceNumberByKey("kp",args,kp)
	args = ReplaceNumberByKey("lp",args,lp)
	args = ReplaceNumberByKey("cone",args,cone)
	args = ReplaceNumberByKey("maxSpots",args,maxSpots)
	args = ReplaceNumberByKey("depth",args,depth)
	args = ReplaceNumberByKey("printIt",args,printIt)

	String funcName = StringFromList(0,GetIndexingFuncOrExec())
	FUNCREF runIndexingEulerCommand runIndexingFunc = $funcName
	Wave FullPeakIndexed = runIndexingFunc(args)
	if (!WaveExists(FullPeakIndexed))
		if (printIt)
			print "\tNothing indexed"
		endif
		return $""
	endif
	String wnote = note(FullPeakIndexed)
	if (numtype(depth)==0)									// transfer valid depth to FullPeakIndexed wave note
		Note/K FullPeakIndexed, ReplaceNumberByKey("depth",wnote,depth,"=")
	endif
	Variable NpatternsFound = NumberByKey("NpatternsFound",wnote,"=")
	Variable Nindexed = NumberByKey("Nindexed",wnote,"=")
	Variable NiData = NumberByKey("NiData",wnote,"=")
	Variable executionTime = NumberByKey("executionTime",wnote,"=")
	Variable rms_error = NumberByKey("rms_error0",wnote,"=")

	if (printIt)
		Make/N=(DimSize(FullPeakIndexed,0))/FREE keVs=NaN
		keVs = FullPeakIndexed[p][7][0]
		keVs = numtype(keVs) ? -Inf : keVs
		Variable keVmax = WaveMax(keVs)
		String secStr
		sprintf secStr, "%.3g", executionTime
		String timeStr=SelectString(executionTime>=60,secStr+" sec", Secs2Time(executionTime,5,1)+" ("+secStr+" sec)")
		String str = ReplaceString(";", GetIndexingFuncOrExec(), "")		// removes all ";" is like a concatenate
		str = SelectString(strlen(str),"Euler",str)
		printf "from \"%s\", found %d patterns, indexed %d out of %d spots (Emax=%.3f keV)  with rms=%.3g"+DEGREESIGN+" in %s",str,NpatternsFound,Nindexed,NiData,keVmax,rms_error,timeStr
		printf ",   "+SurfaceNormalString(FullPeakIndexed)+"\r"
		if (NpatternsFound>1)
			Wave RL0=str2recip(StringByKey("recip_lattice0",wnote,"="))
			Variable pat, angle
			for (pat=1;pat<NpatternsFound;pat+=1)
				keVs = FullPeakIndexed[p][7][pat]
				keVs = numtype(keVs) ? -Inf : keVs
				keVmax = WaveMax(keVs)
				rms_error = NumberByKey("rms_error"+num2istr(pat),wnote,"=")
				printf "  also found pattern %d (Emax=%.3f keV)  with rms=%g"+DEGREESIGN+",   %s\r",pat,keVmax,rms_error,SurfaceNormalString(FullPeakIndexed,pattern=pat)
				Wave RL=str2recip(StringByKey("recip_lattice"+num2istr(pat),wnote,"="))
				Make/N=3/D/FREE hkl3
				angle = rotationBetweenRecipLattices(RL0,RL,hkl3,ints=1)
				printf "  pattern%d  is a rotation of %.3f"+DEGREESIGN+" about an (%s)  from pattern0\r",pat,angle,hkl2str(hkl3[0],hkl3[1],hkl3[2])
			endfor
		endif
		if (INDEXING_ASK_BEFORE_DISPLAYING)
			DoAlert 1, "Examine fit on a picture"
			if (V_flag==1)
				DisplayResultOfIndexing(FullPeakIndexed,NaN)
			endif
		else
			DisplayResultOfIndexing(FullPeakIndexed,NaN)
		endif
	endif
	return FullPeakIndexed
End
//
Static Function rotationBetweenRecipLattices(RL0,RL,hkl,[ints])		// returns the hkl, and angle
	Wave RL0, RL				// want rotation from RL0 to RL (both are reciprocal lattices)
	Wave hkl					// returned vector with hkl
	Variable ints				// if true, convert hkl to integers
	ints = ParamIsDefault(ints) ? 0 : !(!ints)

	MatrixOP/FREE/O rot = RL x Inv(RL0)
	Make/N=3/D/FREE axis
	Variable angle = axisOfMatrix(rot,axis,squareUp=1)	// angle (degree)
	MatrixOP/O hkl = Inv(RL0) x axis

	Make/N=3/D/FREE mv, remain
	mv = hkl==0 ? Inf : abs(hkl)
	Variable mm = WaveMin(mv)
	if (mm>0)
		hkl = hkl/mm
	endif
	if (ints)								// convert hkl to nearby integers
		Make/N=3/D/FREE remain
		hkl = round(hkl*40)
		mm = max(-WaveMin(hkl),WaveMax(hkl))
		Variable i
		for (i=mm;i>=2;i-=1)				// remove common factors
			remain = abs(mod(hkl[p],i))	// if i divides evenly into each of hkl, then divide hkl by i
			if (sum(remain)==0)
				hkl /= i
			endif
		endfor
	endif
	return angle
End
//
Static Function/T SurfaceNormalString(FullPeakIndexed,[pattern])
	Wave FullPeakIndexed
	Variable pattern
	if (ParamIsDefault(pattern) || !(pattern>=0))
		pattern = 0
	endif

	String str = StringByKey("recip_lattice"+num2istr(pattern),note(FullPeakIndexed),"=")
	Make/N=(3,3)/FREE/D recip=NaN
	Variable a0,a1,a2, b0,b1,b2, c0,c1,c2
	sscanf str,"{{%g,%g,%g}{%g,%g,%g}{%g,%g,%g}}",a0,a1,a2, b0,b1,b2, c0,c1,c2
	if (V_flag!=9)
		return ""
	endif
	recip[0][0]=a0;	recip[0][1]=b0;	recip[0][2]=c0
	recip[1][0]=a1;	recip[1][1]=b1;	recip[1][2]=c1
	recip[2][0]=a2;	recip[2][1]=b2;	recip[2][2]=c2

	Make/N=3/D/FREE normal=NaN, surf={0,1,-1}
	MatrixOp/FREE/O normal = Inv(recip) x surf
	Make/N=3/D/FREE hkl
	hkl = abs(normal)
	Variable maxVal = WaveMax(hkl)
	hkl = round(normal/maxVal * 24)
	hkl = hkl==0 ? 0 : hkl
	Variable h=hkl[0],k=hkl[1],l=hkl[2]
	lowestOrderHKL(h,k,l)
	hkl = {h,k,l}
	MatrixOp/FREE/O qvec = recip x hkl
	normal = {0,1,-1}
	Variable angle = angleVec2Vec(normal,qvec)
	if (pattern>0)
		sprintf str, "surface normal(pattern=%d) is near (%s), which is %g"+DEGREESIGN+" from assumed normal",pattern,hkl2str(h,k,l),angle
	else
		sprintf str, "surface normal is near (%s), which is %g"+DEGREESIGN+" from assumed normal",hkl2str(h,k,l),angle
	endif
	return str
End
//
Static Function pick_keV_calc()				// given xtal, return a good number for the keV calc to be used in Euler indexing program
	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))
		DoAlert 0, "no lattice structure found, did you forget to set it?"
		Abort "in pick_keV_calc()"
	endif
	//	Variable factor = NumVarOrDefault(MICRO_GEOMETRY_VERSION_PATH,0)&2 ? 3.0 : 5.0
//	Variable factor = 5.0
	Variable factor = 4.0
	Variable keV = factor * (PrimitiveCellFactor(xtal)/(xtal.Vc))^0.33333
	keV = round(keV*10)/10						// round to 0.1 keV
	return keV
End


Function/T MakeIndexedWaveForAuxDetector(dNum,peakList,indexedList)	// create the FullPeakIndexed wave for Aux detectors
	// modified to do all of the indexed patterns
	Variable dNum										// this should probably only be called with dNum = 1 or 2
	Wave peakList										// the FullPeakList wave
	Wave/Z indexedList								// the FullPeakIndexed wave, only used to get the reciprocal lattice

	STRUCT microGeometry g
	if (FillGeometryStructDefault(g, alert=1))	//fill the geometry structure with current values
		return ""
	endif

	String flist = reverseList(WaveListClass("FittedPeakList*","*","DIMS:2"))
	Variable perfect=0
	if (!WaveInClass(indexedList,"IndexedPeakList*") || !WaveInClass(peakList,"FittedPeakList*") || !(dNum==limit(round(dNum),0,2)))
		String indexName="", peakListName=""
		Wave image = ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";" )))
		String win=StringFromList(0,WindowsWithWave(image,1))
		if (!(dNum == limit(round(dNum),0,2)))				// invalid dNum
			dNum = detectorNumFromID(g, StringByKey("detectorID", note(image),"="))
		endif
		dNum = !(dNum==limit(round(dNum),0,2)) ? 1 : dNum
		if (!WaveExists(peakList))
			Wave peakList = FindPeakListForImage(image)
		endif
		if (WaveExists(peakList))
			peakListName = NameOfWave(peakList)
		endif
		if (!WaveExists(indexedList))
			Wave indexedList = $StringByKey("FullPeakIndexed",GetUserData(win,"","Indexing"),"=")
		endif
		if (WaveExists(indexedList))
			indexName = NameOfWave(indexedList)
		endif
		Prompt dNum,"Detector Number (probably 1 or)",popup,DetectorMenuList(g)
		String ilist = reverseList(WaveListClass("IndexedPeakList*","*",""))+"perfect Si;"
		Prompt indexName, "FullPeakIndexed List to provide sample rotation",popup,ilist
		Prompt peakListName, "Fitted Peaks List",popup,flist
		dNum += 1
		DoPrompt "Sample Rotation",dNum,peakListName,indexName
		dNum -= 1
		if (V_flag)
			return ""
		endif
		perfect = stringmatch(indexName,"perfect Si")
		Wave peakList=$peakListName, indexedList=$indexName
		printf "MakeIndexedWaveForAuxDetector(%g,%s,%s)\r",dNum,NameOfWave(peakList),NameOfWave(indexedList)
	endif
	String AuxPeakIndexedName="AuxPeakIndexed"+ReplaceString("FullPeakList",NameOfWave(peakList),"")
	Variable Npatterns = max(1,DimSize(indexedList,2))			// always at least 1
	String indexNote=note(indexedList)
	Make/N=(3,3,Npatterns)/FREE/D recipN=NaN, recipInvN=NaN	// stores a recip for each pattern
	Make/N=(3,3)/FREE/D recip=NaN, recipInv=NaN
	Variable angleTolerance = NumberByKey("angleTolerance",indexNote,"=")
	angleTolerance = 2*( numtype(angleTolerance) ? 0.1 : angleTolerance )
	Variable depth = NumberByKey("depth",indexNote,"=")

	if (perfect)
		Make/N=3/D/FREE axis={-2.2256407716336,-0.884037203122317,-0.38914203759791}	// length is radians
		Variable len = normalize(axis)
		len = tan(len/2)											// length of a true Rodriques vector tan(angle/2)
		axis *= len
		rotationMatAboutAxis(axis,NaN,recip)
		MatrixOp/FREE/O recipInv = Inv(recip)
		recipN[][][0] = recip[p][q]
		recipInvN[][][0] = recipInv[p][q]
	else
		Variable i, r00,r10,r20, r01,r11,r21, r02,r12,r22
		String list
		for (i=0;i<Npatterns;i+=1)
			list = StringByKey("recip_lattice"+num2istr(i),indexNote,"=")
			sscanf list, "{{%g,%g,%g}{%g,%g,%g}{%g,%g,%g}}", r00,r10,r20, r01,r11,r21, r02,r12,r22
			if (V_flag==9)
				recip[0][0]=r00;	recip[0][1]=r01;	recip[0][2]=r02
				recip[1][0]=r10;	recip[1][1]=r11;	recip[1][2]=r12
				recip[2][0]=r20;	recip[2][1]=r21;	recip[2][2]=r22
				MatrixOp/FREE/O recipInv = Inv(recip)
				recipN[][][i] = recip[p][q]
				recipInvN[][][i] = recipInv[p][q]
			endif
		endfor
	endif

	Make/N=3/FREE/D qmeas
	Variable N=DimSize(peakList,0)								// number of fitted peaks to examine
	Make/N=(N,12,Npatterns)/O/D $AuxPeakIndexedName=NaN		// this holds the result
	Wave AuxPeakIndexed=$AuxPeakIndexedName

	Make/N=(Npatterns)/FREE rmsPix=0,rmsAng=0, Nindexed=0
	Variable/C pz
	Variable theta, ipat, errorAng, j, NindexedTotal, found1
	for (NindexedTotal=0,i=0; i<N; i+=1)						// loop over each measured peak
		theta = pixel2q(g.d[dNum],peakList[i][0],peakList[i][1],qmeas,depth=depth)
		found1 = 0
		for (ipat=0;ipat<Npatterns;ipat+=1)
			recipInv = recipInvN[p][q][ipat]
			MatrixOp/FREE/O hkl = recipInv x qmeas				// un-rotate from measured to reference orientation
			normalize(hkl)
			CloseAllowedhkl(hkl)
			recip = recipN[p][q][ipat]
			MatrixOp/FREE/O qcalc = recip x hkl
			normalize(qcalc)
			errorAng = acos(limit(MatrixDot(qcalc,qmeas),-1,1))*180/PI	// angle error (degree)
			if (errorAng<=angleTolerance)							// found a match
				j = Nindexed[ipat]
				AuxPeakIndexed[j][3,5][ipat] = hkl[q-3]
				AuxPeakIndexed[j][7][ipat] = hc / (2*get_dhkl(hkl[0],hkl[1],hkl[2])*sin(theta))	//   hc/E = 2 d sin(theta)
				AuxPeakIndexed[j][6][ipat] = peakList[i][10]	// transfer the integral
				AuxPeakIndexed[j][0,2][ipat] = qmeas[q]			// actual measured Q-vector (normalized)
				AuxPeakIndexed[j][8][ipat] =	errorAng				// angle error (degree)
				rmsAng[ipat] += (errorAng)^2
				pz = q2pixel(g.d[dNum],qcalc)
				AuxPeakIndexed[j][9][ipat] = real(pz)
				AuxPeakIndexed[j][10][ipat] = imag(pz)
				rmsPix[ipat] += (peakList[i][0]-real(pz))^2 + (peakList[i][1]-imag(pz))^2
				Nindexed[ipat] +=1
				found1 = 1
			endif
		endfor
		NindexedTotal += found1
	endfor
	Variable NindexedMax=WaveMax(Nindexed)
	Redimension/N=(NindexedMax,-1,-1) AuxPeakIndexed	// trim off unused end
	AuxPeakIndexed[][11][] = dNum
	rmsPix = sqrt(rmsPix/Nindexed[p])
	rmsAng = sqrt(rmsAng/Nindexed[p])
	if (NindexedTotal<1)
		rmsPix = NaN
		rmsAng = NaN
	endif
	printf "rms error = %spixels,   or   %s%s\r",vec2str(rmsPix,places=3,bare=Npatterns<2),vec2str(rmsAng,places=4,bare=Npatterns<2),DEGREESIGN

	pixel2q(g.d[0],925,1095,qcalc,depth=depth)
	recipInv = recipInvN[p][q][0]
	MatrixOp/FREE/O qcalc = recipInv x qcalc
	CloseAllowedhkl(qcalc)
	printf "on top detector,   hkl = (%g %g %g)",qcalc[0],qcalc[1],qcalc[2]
	recip = recipN[p][q][0]
	MatrixOp/FREE/O qcalc = recip x qcalc
	pz = q2pixel(g.d[0],qcalc)
	printf " --> pixel = (%d, %d)",real(pz),imag(pz)
	if (Npatterns>1)
		printf " for pattern 0"
	endif
	printf "\r"

	SetDimLabel 1,0,Qx,AuxPeakIndexed			;	SetDimLabel 1,1,Qy,AuxPeakIndexed
	SetDimLabel 1,2,Qz,AuxPeakIndexed			;	SetDimLabel 1,3,h,AuxPeakIndexed
	SetDimLabel 1,4,k,AuxPeakIndexed			;	SetDimLabel 1,5,l,AuxPeakIndexed
	SetDimLabel 1,6,Intensity,AuxPeakIndexed;	SetDimLabel 1,7,keV,AuxPeakIndexed
	SetDimLabel 1,8,angleErr,AuxPeakIndexed	;	SetDimLabel 1,9,pixelX,AuxPeakIndexed
	SetDimLabel 1,10,pixelY,AuxPeakIndexed	;	SetDimLabel 1,11,detNum,AuxPeakIndexed
	String wNote="waveClass=IndexedPeakListAux;"
	wNote = ReplaceStringByKey("structureDesc",wnote,StringByKey("structureDesc",indexNote,"="),"=")
	wNote = ReplaceStringByKey("latticeParameters",wnote,StringByKey("latticeParameters",indexNote,"="),"=")
	wNote = ReplaceStringByKey("lengthUnit",wnote,StringByKey("lengthUnit",indexNote,"="),"=")
	wNote = ReplaceStringByKey("SpaceGroup",wnote,StringByKey("SpaceGroup",indexNote,"="),"=")
	String str = StringByKey("SpaceGroupID",indexNote,"=")
	if (strlen(str))
		wNote = ReplaceStringByKey("SpaceGroupID",wnote,str,"=")
	endif
	str = StringByKey("SpaceGroupIDnum",indexNote,"=")
	if (strlen(str))
		wNote = ReplaceStringByKey("SpaceGroupIDnum",wnote,str,"=")
	endif
	str = StringByKey("AtomDesctiption1",indexNote,"=")
	for (i=1;strlen(str);i+=1)
		wNote = ReplaceStringByKey("AtomDesctiption"+num2istr(i),wnote,str,"=")
		str = StringByKey("AtomDesctiption"+num2istr(i+1),indexNote,"=")
	endfor
	wNote = ReplaceStringByKey("peakListWave",wnote, GetWavesDataFolder(peakList,2),"=")
	wNote = ReplaceStringByKey("fittedIgorImage",wnote,StringByKey("fittedIgorImage",note(peakList),"="),"=")
	for (ipat=0;ipat<Npatterns;ipat+=1)
		String d=num2istr(ipat)
		wNote = ReplaceStringByKey("EulerAngles"+d,wnote,StringByKey("EulerAngles"+d,indexNote,"="),"=")
		wNote = ReplaceStringByKey("rotation_matrix"+d,wnote,StringByKey("rotation_matrix"+d,indexNote,"="),"=")
		wNote = ReplaceStringByKey("recip_lattice"+d,wnote,StringByKey("recip_lattice"+d,indexNote,"="),"=")
		wNote = ReplaceStringByKey("structureDesc",wnote,StringByKey("structureDesc",indexNote,"="),"=")
		wNote = ReplaceNumberByKey("rms_error"+d,wnote,rmsAng[ipat],"=")
	endfor
	wNote = ReplaceNumberByKey("NpatternsFound",wnote,Npatterns,"=")
	wNote = ReplaceNumberByKey("Nindexed",wnote,NindexedTotal,"=")
	wNote = ReplaceNumberByKey("NiData",wnote,N,"=")
	Note/K AuxPeakIndexed,wNote
	return GetWavesDataFolder(AuxPeakIndexed,2)		// can not compute (px,py), but still a valid result
End
//
Function ClearhklOffGraph(win)
	String win
	String item,list = AnnotationList("")
	Variable i
	for (i=0,item=StringFromList(0,list); strlen(item); i+=1,item=StringFromList(i,list))
		if (strsearch(item, "hkl",0)==0)
			Tag/W=$win/K/N=$item
		endif
	endfor
	list = GetUserData("","","Indexing")
	SetWindow kwTopWin,userdata(Indexing)=RemoveByKey("FullPeakIndexed",list,"=")
End
//
Function DisplayResultOfIndexing(FullPeakIndexed,pattern)
	Wave FullPeakIndexed				// contains the result of indexing peaks
	Variable pattern						// pattern in FullPeakIndexed to draw

	Variable i, badWave = !WaveExists(FullPeakIndexed)
	if (badWave)							// try getting FullPeakIndexed from the top image
		Wave image = ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";")))
		Wave FullPeakIndexed = FindIndexListForImage(image)
		badWave = !WaveExists(FullPeakIndexed)
	endif
	if (!badWave)
		badWave = (DimSize(FullPeakIndexed,0)<1 || DimSize(FullPeakIndexed,1)<12)
	endif
	if (badWave)
		String wList = reverseList(WaveListClass("IndexedPeakList*","*",""))
		i = ItemsInList(wList)
		if (i==1)							// only one, so choose it
			Wave FullPeakIndexed = $(StringFromList(0,wList))
		elseif(i>0)							// more than 1, so ask
			String wName=""
			Prompt wName,"Indexed Peak List",popup,wList
			DoPrompt "Indexed Peak List",wName
			if (V_flag)
				return 1
			endif
			Wave FullPeakIndexed = $wName
		endif
	endif
	badWave = !WaveExists(FullPeakIndexed)
	if (!badWave)
		badWave = (DimSize(FullPeakIndexed,0)<1 || DimSize(FullPeakIndexed,1)<12)
	endif
	if (badWave)
		return 1
	endif
	String wnote = note(FullPeakIndexed)
	Variable NpatternsFound = max(DimSize(FullPeakIndexed,2),1)
	NpatternsFound = min(NpatternsFound,NumberByKey("NpatternsFound",wnote,"="))
	if (NpatternsFound<1)
		return 1
	elseif (!(NpatternsFound!=1))
		pattern = 0
	elseif (!(pattern>=0))
		pattern = 0
		Prompt pattern, "pattern number to fit [0,"+num2istr(NpatternsFound-1)+",]",popup,expandRange("0-"+num2istr(NpatternsFound-1),";")
		DoPrompt "pattern number",pattern
		if (V_flag)
			return 1
		endif
		pattern = min(NpatternsFound-1,pattern-1)
	endif
	if (!(pattern>=0))
		return 1
	endif
	// get a graph of the image at the front, with boxes on to the front
	Wave image = $""					// init to NULL
 	String imageName					// name of image to look at
	String win							// name of window displaying image
	imageName = StringByKey("fittedIgorImage",wnote,"=")

	Variable err=0
	if (itemsInList(imageName,",")>=1)
		for (i=0;i<itemsInList(imageName,",");i+=1)
			Wave image = $StringFromList(i,imageName,",")
			win = FindGraphWithImage(image)
			err = UpdateIndexingOn1image(win,image,FullPeakIndexed,pattern) || err
		endfor
	endif

	if (err)
		Wave image = $(StringByKey("rawIgorImage",wnote,"="))
		win = FindGraphWithImage(rawimage)
		err = UpdateIndexingOn1image(win,image,FullPeakIndexed,pattern)
	endif
	return err
End
//

Static Function UpdateIndexingOn1image(win,image,FullPeakIndexed,pattern)
	String win								// window name
	Wave image								// image wave in win
	Wave FullPeakIndexed				// contains the result of indexing peaks
	Variable pattern						// pattern in FullPeakIndexed to draw

	if (!WaveExists(image))			// could not find the image
		return 1
	endif
	if (strlen(win))
		DoWindow/F $win					// bring window with image to the top
	else
		NewImageGraphLocal(image)		// no window exists, so make one showing image
  	endif
	STRUCT WMButtonAction B_Struct
	B_Struct.eventCode = 2
	B_Struct.eventMod = 63				// special to force ON or OFF
	B_Struct.win = ""
	B_Struct.ctrlName="ON"				// make sure that boxes are on
	ButtonBoxesProc(B_Struct)
	ClearhklOffGraph("")

	STRUCT microGeometry g
	FillGeometryStructDefault(g)		//fill the geometry structure with current values
	Variable dNum = detectorNumFromID(g, StringByKey("detectorID",note(image),"="))
	dNum = max(0,dNum)					// default to 0

	String wnote = note(FullPeakIndexed)
	String peakListName = StringByKey("peakListWave", wnote,"=")
	if (strlen(peakListName))
		SetWindow kwTopWin userdata(FitPeaks )="FullPeakList="+peakListName
	endif

	Variable wid = (NumberByKey("maxPeakWidth",wnote,"=")+NumberByKey("minSpotSeparation",wnote,"="))/2
	wid = numtype(wid) ? NumberByKey("minSpotSeparation",wnote,"=") : wid
	//	wid = numtype(wid) ? 30 : wid		// width of the cross
	wid = numtype(wid) ? 20 : wid	// width of the cross
	//	Variable tsize = round(75 - N/4)
	//	tsize = limit(tsize,45,75)
	MatrixOP/FREE num_dNum = sum(equal(col(FullPeakIndexed[][][pattern],11),dNum))
	Variable tsize = round(80 - num_dNum[0]/4)
	tsize = limit(tsize,45,80)

	Variable i, px,py, N=DimSize(FullPeakIndexed,0)
	SetDrawLayer UserFront				// draw the X's
	for (i=0;i<N;i+=1)
		if (dNum != FullPeakIndexed[i][11][pattern])
			continue							// wrong detector, skip this spot
		endif
		px = FullPeakIndexed[i][9][pattern]
		py = FullPeakIndexed[i][10][pattern]
		if (!numtype(px+py))
			Variable fwx=wid*DimSize(image,0)/2048, fwy=wid*DimSize(image,1)/2048
			DrawMarker(px,py,fwx,fwy,"X",color="55000,5000,0")
			DrawhklTags(i,FullPeakIndexed[i][3][pattern],FullPeakIndexed[i][4][pattern],FullPeakIndexed[i][5][pattern],px,py,wid,0,tsize)
		endif
	endfor
	TextBox/K/N=indexedPeakInfo
	TextBox/C/N=indexedPeakInfo/F=0 ""
	String list = GetUserData("","","Indexing")
	list = ReplaceStringByKey("FullPeakIndexed",list,GetWavesDataFolder(FullPeakIndexed,2),"=")
	list = ReplaceNumberByKey("patternNum",list,pattern,"=")
	SetWindow kwTopWin,userdata(Indexing)=list
	return 0
End
//
Static Function DrawhklTags(i,h,k,l,px,py,dx,dy,tsize)
	Variable i
	Variable h,k,l				// hkl value
	Variable px,py				// pixel position
	Variable dx,dy				// off set of label from pixel center
	Variable tsize
	tsize = tsize>2 && tsize<150 ? tsize : 75
	String str
	//	sprintf str,"\\F'Comic Sans MS'\\Zr%03d\\[0",tsize
	sprintf str,"\\F'%s'\\Zr%03d",BAR_FONT_ALWAYS,tsize
	str += SelectString (IgorVersion() > 7, "\\[0", "")
	str += hkl2IgorBarStr(h,k,l)							// change to Igor string with negatives --> bars
	Wave image = ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";")))
	px = limit(px,0,DimSize(image,0)-1)
	py = limit(py,0,DimSize(image,1)-1)
	Variable p = round(px + DimSize(image,0)*round(py))
	String name = "hkl"+num2istr(i)
	String wname=NameOfWave(image)
	Tag/C/N=$name/I=1/F=0/G=(55000,5000,0)/B=1/A=MB/L=0/X=1/Y=0.5 $wname,p,str
	return 0
End
//
// use version of hkl2IgorBarStr() that is in Latticesym.ipf
//ThreadSafe Static Function/T hkl2IgorBarStr(h,k,l)	// changes negatives to a bar over the number, only for displaying, not printing
//	Variable h,k,l				// hkl value
//
//	if (numtype(h+k+l))
//		return num2str(h)+","+num2str(k)+","+num2str(l)
//	endif
////	String extra=SelectString(abs(h)>9 || abs(k)>9 || abs(l)>9 || h<0 || k<0 || l<0,""," ")
//	String extra=SelectString(abs(h)>9 || abs(k)>9 || abs(l)>9,""," ")
//	String str=""
//	str += minus2bar(num2istr(h),spaces=floor(log(abs(h)))) + extra
//	str += minus2bar(num2istr(k),spaces=floor(log(abs(k)))) + extra
//	str += minus2bar(num2istr(l),spaces=floor(log(abs(l))))
//	return str
//End
//
//Static Function DrawhklTags(i,h,k,l,px,py,dx,dy,tsize)
//	Variable i
//	Variable h,k,l				// hkl value
//	Variable px,py				// pixel position
//	Variable dx,dy				// off set of label from pixel center
//	Variable tsize
//	tsize = tsize>2 && tsize<150 ? tsize : 75
//	String str
//	sprintf str,"\\F'Comic Sans MS'\\Zr%03d\\[0",tsize		// use for gfMult=100
//	if (stringmatch(IgorInfo(2),"Macintosh"))
//		str += SelectString(h<0,"","\\S \\f01—\\f00\\M\\X0")+" "+num2istr(abs(h))	// mac
//		str += SelectString(k<0,"","\\[1\\S\\f01 —\\f00\\M\\X1")+" "+num2istr(abs(k))
//		str += SelectString(l<0,"","\\[2\\S \\f01—\\f00\\M\\X2")+" "+num2istr(abs(l))
//	else
//		str += SelectString(h<0,"","\\S \\f01Ø\\f00\\M\\X0")+" "+num2istr(abs(h))	// windows
//		str += SelectString(k<0,"","\\[1\\S\\f01 Ø\\f00\\M\\X1")+" "+num2istr(abs(k))
//		str += SelectString(l<0,"","\\[2\\S \\f01Ø\\f00\\M\\X2")+" "+num2istr(abs(l))
//	endif
//	Wave image = ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";")))
//	px = limit(px,0,DimSize(image,0)-1)
//	py = limit(py,0,DimSize(image,1)-1)
//	Variable p = round(px + DimSize(image,0)*round(py))
//	String name = "hkl"+num2istr(i)
//	String wname=NameOfWave(image)
//	Tag/C/N=$name/I=1/F=0/G=(55000,5000,0)/B=1/A=MB/L=0/X=1/Y=0.5 $wname,p,str
//	return 0
//End
//
Function IndexedPeakInfoFromMarquee()	// find the energy of a spot in the marquee
	String imageName = StringFromList(0,ImageNameList("",";"))
	Wave image = ImageNameToWaveRef("",imageName)
	STRUCT box b
	if (pointRangeOfMarquee("",imageName,b)!=0)
		return NaN
	endif

	String list = GetUserData("","","Indexing")
	Variable ip = NumberByKey("patternNum",list,"=")		// pattern number
	Wave FullPeakIndexed=$(StringByKey("FullPeakIndexed",list,"="))
	if (!WaveExists(FullPeakIndexed))
		DoAlert 0,"FullPeakIndexed cannot be found"
		return NaN
	endif
	ip = (ip<0 || numtype(ip)) ? 0 : ip
	ip = limit(ip,0,max(DimSize(FullPeakIndexed,2)-1,0))
	Variable px,py,i,N=DimSize(FullPeakIndexed,0)
	for (i=0;i<N;i+=1)				// first find the peak
		px = FullPeakIndexed[i][9][ip]
		py = FullPeakIndexed[i][10][ip]
		if (b.xlo<px && px<b.xhi && b.ylo<py && py<b.yhi)
			break
		endif
	endfor
	if (i>=N)
		DoAlert 0,"No indexed peak is in the Marquee"
		return NaN
	endif
	Variable h=FullPeakIndexed[i][3][ip], k=FullPeakIndexed[i][4][ip], l=FullPeakIndexed[i][5][ip]
	Variable Intens=FullPeakIndexed[i][6][ip], keV=FullPeakIndexed[i][7][ip],angleErr=FullPeakIndexed[i][8][ip]
	Variable Np = NumberByKey("NpatternsFound",note(FullPeakIndexed),"=")
	if (Np>1)
		printf "Indexed peak %d at pixel(%g, %g),     in pattern %d of %d, hkl=(%d %d %d),   Intensity=%g,  Energy=%.4f keV,  angleErr=%.3f%s\r",i,px,py,ip,Np,h,k,l,Intens,keV,angleErr,DEGREESIGN
	else
		printf "Indexed peak %d at pixel(%g, %g),     hkl=(%d %d %d),   Intensity=%g,  Energy=%.4f keV,  angleErr=%.3f%s\r",i,px,py,h,k,l,Intens,keV,angleErr,DEGREESIGN
	endif
End
// This puts up a window with info about an indexed peak on the current graph, when you shift-click on a spot
//
//	if (exists("getIndexedPeakInfoHook")==6)
//		SetWindow Graph1 hook(peakInfo)=getIndexedPeakInfoHook
//	endif
//
//
//
Function getIndexedPeakInfoHook(s)			// old name, use getFittedPeakInfoHook() for all future use
	STRUCT WMWinHookStruct &s
	return getFittedPeakInfoHook(s)
End
//
Function getFittedPeakInfoHook(s)	// Command=fitted peak,  Shift=Indexed peak,  Command-Shift=strained peak,  &  Command-mouseDown=follow mouse
	STRUCT WMWinHookStruct &s

	// use shift to look at the indexing result, use command to look at fitted peaks
	String win = (s.winName)
	Variable process = (s.eventMod&10 && (s.eventCode==4 || s.eventCode==3)) && (s.keycode==0 || s.keycode==255)		// eventMod==2 is shift key, 8 is command, 1 is mouse, and no regular key held down
	if (s.eventCode==8)
		return 0							// graph modified occurs when the tag is made or changed, ignore this
	elseif (!process)
		Tag/K/N=indexedPeakInfo/W=$win
		s.doSetCursor= 0
		s.cursorCode = 0				// for nothing
		return 0							// shift key not down, ignore
	endif

	Variable useIndex = s.eventMod&2		// eventMod==2 is shift key,  use info from indexed peaks
	Variable useFitted = s.eventMod&8		// eventMod==8 is command key, use info from fitted peaks
	Variable useStrain = useIndex && useFitted
	useFitted = useFitted && !useStrain	// both means use strain
	Variable useMouse = s.eventMod==9		// evendMod==8+1, is command & mouse down, use mouse position
	if (useMouse)
		useIndex = 0
		useFitted = 0
		useStrain = 0
	endif

	s.doSetCursor= 1
	if (useStrain)
		s.cursorCode = 19					// lightning bolt for strain
	elseif (useIndex)
		s.cursorCode = 26					// "?" in box for index
	elseif (useFitted)
		s.cursorCode = 17					// big "?" for fitted
	elseif (useMouse)
		s.cursorCode = 18					// small box with cross for only mouse
	else
		s.cursorCode = 0					// for nothing
		s.doSetCursor= 0
	endif

	Variable i,N, useMissing=0, useUser=0
	if (useFitted)							// check if useFitted also include Missing
		String tname = StringFromLIst(0,TraceNamesInClass("MissingPeakList*",win))
		if (strlen(tname))
			Wave missing = TraceNameToWaveRef(win,tname)
		else
			Wave missing = $""
		endif
		useMissing = WaveExists(missing)

		tname = StringFromLIst(0,TraceNamesInClass("UserRecipPeakList*",win))
		if (strlen(tname))
			Wave userHKL = TraceNameToWaveRef(win,tname)
		else
			Wave userHKL = $""
		endif
		useUser = WaveExists(userHKL)
	endif

	String imageName = StringFromList(0,ImageNameList(win,";")), wnote=""
	String FullPeakName="", FullPeakIndexName="", PeaksForStrainName="", TraceName=""
	Wave image = ImageNameToWaveRef(win,StringFromList(0,ImageNameList(win,";")))
	if (!WaveExists(image))
		FullPeakName = StringFromList(0,TraceNamesInClass("FittedPeakList*",win))
		Wave FullPeakList = TraceNameToWaveRef(win,FullPeakName)
		if (!WaveExists(FullPeakList))
			return 1
		endif
		wnote = note(FullPeakList)
	elseif (WaveDims(image)!=2)
		return 1
	else
		wnote = note(image)
	endif

	STRUCT microGeometry geo
	if (FillGeometryStructDefault(geo, alert=1))	//fill the geometry structure with test values
		return 1
	endif

	Variable startx,groupx, starty,groupy	// ROI of the actual image
	startx = NumberByKey("startx",wnote,"=")
	groupx = NumberByKey("groupx",wnote,"=")
	starty = NumberByKey("starty",wnote,"=")
	groupy = NumberByKey("groupy",wnote,"=")
	startx = numtype(startx) ? FIRST_PIXEL : startx
	groupx = numtype(groupx) ? 1 : groupx
	starty = numtype(starty) ? FIRST_PIXEL : starty
	groupy = numtype(groupy) ? 1 : groupy

	Variable dNum = max(detectorNumFromID(geo, StringByKey("detectorID", wnote,"=")),0)
	Make/N=3/D/FREE Pvec = geo.d[dNum].P[p]	// get theta resolution from detector (resolution to 1/2 pixel in theta)
	Variable perpDist = norm(Pvec)
	Variable dangleX = groupx * geo.d[dNum].sizeX / geo.d[dNum].Nx
	Variable dangleY = groupy * geo.d[dNum].sizeY / geo.d[dNum].Ny
	Variable anglePlaces = -floor( log(min(dangleX,dangleY)/4 / perpDist * 180/PI) )
	anglePlaces = numtype(anglePlaces) ? 4 : limit(anglePlaces,1,8)
	String angFmt = "%."+num2istr(anglePlaces)+"f"+S_DEGREE

	GetWindow $win psize
	Variable vert = limit((s.mouseLoc.v-V_top)/(V_bottom-V_top),0,1)	// fractional position on graph
	Variable horiz = limit((s.mouseLoc.h-V_left)/(V_right-V_left),0,1)

	Variable mx,my										// pixel position at the mouse-click
	GetAxis/W=$win/Q bottom
	if (V_flag)
		return 1
	endif
	mx = (V_max-V_min)*horiz + V_min
	GetAxis/W=$win/Q left
	if (V_flag)
		return 1
	endif
	my = (V_min-V_max)*vert + V_max

	String tagStr="", str, SpaceGroupID="", hklStr=""
	Variable h,k,l, angleErr,keV=NaN, SpaceGroupIDnum
	Variable dist2, m
	Variable px,py
	Variable ip = NumberByKey("patternNum",GetUserData("","","Indexing"),"=")	// pattern number
	ip = (ip<0 || numtype(ip)) ? 0 : ip

	Wave FullPeakIndexed=$(StringByKey("FullPeakIndexed",GetUserData("","","Indexing"),"="))
	if (!WaveExists(FullPeakIndexed))
		FullPeakIndexName = StringFromList(0,TraceNamesInClass("IndexedPeakList*",win))
		Wave FullPeakIndexed = TraceNameToWaveRef(win,FullPeakIndexName)
	else
		FullPeakIndexName = NameOfWave(FullPeakIndexed)
	endif
	if (WaveExists(FullPeakIndexed))
		ip = limit(ip,0,max(DimSize(FullPeakIndexed,2)-1,0))
		Wave PeaksForStrain=$( GetWavesDataFolder(FullPeakIndexed,1)+"PeaksForStrain")
		PeaksForStrainName = NameOfWave(PeaksForStrain)
	endif
	useStrain = useStrain && (modDate(PeaksForStrain) >= modDate(FullPeakIndexed))	// only use strain if not stale
	useStrain = useStrain && (NumberByKey("patternNum", note(PeaksForStrain),"=")==ip)
	useIndex = useIndex && !useStrain
	Variable useDistortion = NumVarOrDefault("root:Packages:geometry:useDistortion",USE_DISTORTION_DEFAULT)
	Variable wIndex=NaN

	if (useFitted || useMouse)				// examining fitted peaks (or mouse position), also includes missing & userHKL
		if (!WaveExists(FullPeakList))
			Wave FullPeakList=$(StringByKey("FullPeakList",GetUserData("","","FitPeaks"),"="))
			FullPeakName = NameOfWave(FullPeakList)
		endif
		if (!WaveExists(FullPeakList))
			Wave FullPeakList=FindPeakListForImage(image)
		endif
		if (!WaveExists(FullPeakList) && useFitted)
			return 0
		endif

		Variable fwx=NaN, fwy=NaN, xc=NaN
		if (useFitted)
			Variable dmiss=Inf, imiss=-1, duser=Inf, iuser=-1
			if (useMissing)
				N = DimSize(missing,0)
				Make/N=2/FREE mxmy={mx,my}								// location from mouse
				Make/N=(N,2)/FREE pxpy = missing[p][3+q] - mxmy[q]
				MatrixOP/FREE dpxpy2 = sumRows(magSqr(pxpy))		// dist^2 from mxmy to all fitted peaks
				WaveStats/M=1/Q dpxpy2
				imiss = V_minLoc
				dmiss = V_min
				wIndex = imiss
				TraceName = NameOfWave(missing)
			endif
			if (useUser)
				N = DimSize(userHKL,0)
				Make/N=2/FREE mxmy={mx,my}								// location from mouse
				Make/N=(N,2)/FREE pxpy = userHKL[p][3+q] - mxmy[q]
				MatrixOP/FREE dpxpy2 = sumRows(magSqr(pxpy))		// dist^2 from mxmy to all fitted peaks
				WaveStats/M=1/Q dpxpy2
				iuser = V_minLoc
				duser = V_min
				wIndex = iuser
				TraceName = NameOfWave(userHKL)
			endif
			dist2=Inf
			m=-1								// find the fitted peak closest to the mouse-click
			N = DimSize(FullPeakList,0)
			Make/N=2/FREE mxmy={mx,my}								// location from mouse
			Make/N=(N,2)/FREE pxpy = FullPeakList[p][q] - mxmy[q]
			MatrixOP/FREE dpxpy2 = sumRows(magSqr(pxpy))		// dist^2 from mxmy to all fitted peaks
			WaveStats/M=1/Q dpxpy2
			m = V_minLoc
			dist2 = V_min
			wIndex = m
			TraceName = FullPeakName
			if (numtype(dist2) && numtype(dmiss))
				return 1						// no fitted or missing peak found
			endif

			if (duser<dmiss && duser<dist2)	// use userHKL
				useFitted = 0			
				useMissing = 0			
				px = userHKL[iuser][3]			// user peaks, already binned
				py = userHKL[iuser][4]
				keV = userHKL[iuser][5]
			elseif (dmiss<dist2)					// use missing peak
				useFitted = 0			
				useUser = 0
				px = missing[imiss][3]			// missing peaks, already binned
				py = missing[imiss][4]
				keV = missing[imiss][5]
			else										// use fitted peaks
				useMissing = 0
				useUser = 0
				px = FullPeakList[m][0]			// binned peaks
				py = FullPeakList[m][1]
				fwx=FullPeakList[m][4]
				fwy=FullPeakList[m][5]
				xc=FullPeakList[m][8]
			endif
		elseif (useMouse)
			px = mx
			py = my
		else
			Abort "not useFitted or useMouse"
		endif

		if (WaveExists(FullPeakList))
			wnote = note(FullPeakList)
		else
			wnote = note(image)		
		endif
		Make/N=3/O/D/FREE qBL
		Variable pxUnb = (startx-FIRST_PIXEL) + groupx*px + (groupx-1)/2	// change to un-binned pixels
		Variable pyUnb = (starty-FIRST_PIXEL) + groupy*py + (groupy-1)/2	// pixels are still zero based
		Variable depth = NumberByKey("depth",wnote,"=")
		if (numtype(keV) || keV<=0)		// try to get the energy
			if (StringMatch(StringByKey("MonoMode",wnote,"="),"monochromatic"))	// mono mode, try to get keV from wnote
				keV = NumberByKey("keV",wnote,"=")
			endif
		endif
		keV = (numtype(keV)==0 && keV>0) ? keV : NaN

		Variable theta = pixel2q(geo.d[dNum],pxUnb,pyUnb,qBL,depth=depth)	// get theta, and q^ in Beam Line system
		Variable/C pz = useDistortion ? microGeo#peakCorrect(geo.d[dNum],pxUnb,pyUnb) : cmplx(0,0)

		Variable Qmag=NaN
		if (useFitted)
			sprintf tagStr,"\\Zr090Fitted peak position (%.2f, %.2f)",px,py
			if (numtype(fwx+fwy)==0)
				sprintf str,"\rFWHM: x=%.2f, y=%.2f,  x-corr=%.3f",fwx,fwy,xc
				tagStr += str
			endif
			sprintf str,"\r"+S_THETA+" = "+angFmt+",   #%d",theta*180/PI,m
			tagStr += str
		elseif (useMissing)
			hklStr = hkl2str(missing[imiss][0],missing[imiss][1],missing[imiss][2], bar=1)
			sprintf tagStr,"\\Zr090Missing peak: (%.2f, %.2f)\r(%s) at %.4f keV,   #%d",px,py,hklStr,keV,imiss
			String desc = StringByKey("xtalDesc",note(missing),"=")
			if (strlen(desc)>0)
				tagStr +="\r"+desc
			endif
			if (ip>0)
				tagStr +="\r pattern "+num2istr(ip)
			endif
			tagStr +="\\F]0"
		elseif (useUser)
			hklStr = hkl2str(userHKL[iuser][0],userHKL[iuser][1],userHKL[iuser][2], bar=1)
			sprintf tagStr,"\\Zr090User Recip peak: (%.2f, %.2f)\r(%s) at %.4f keV,   #%d\r%s\\F]0",px,py,hklStr,keV,iuser,NameOfWave(userHKL)
		elseif (useMouse)
			if (useDistortion)
				sprintf tagStr,"\\Zr090Mouse position (%.2f, %.2f)\r"+S_THETA+" = "+angFmt+",     distort=%.2f px",px,py,theta*180/PI,cabs(pz)
			else
				sprintf tagStr,"\\Zr090Mouse position (%.2f, %.2f)\r"+S_THETA+" = "+angFmt,px,py,theta*180/PI
			endif
			if (numtype(keV))
				sprintf str,"\r\[0\\Zr075q\X0\y+20^\y-20\BBL\\M\\Zr075 = {%.3f, %.3f, %.3f}\M\\Zr090",qBL[0],qBL[1],qBL[2]
			else
				Qmag = keV>0 ? 4*PI*sin(theta)*keV/hc : 1	// 4*PI*sin(theta)/lambda = 4*PI*sin(theta)*E/hc
				sprintf str,"\r\[0\\Zr075\\f01q\f00\BBL\\M\\Zr075 = {%.3f, %.3f, %.3f},  |q| = %.4f\M\\Zr090",qBL[0]*Qmag,qBL[1]*Qmag,qBL[2]*Qmag, Qmag
			endif
			tagStr += str
		endif

		if (WaveExists(FullPeakIndexed) && !useMissing && !useUser)	// an indexation exists, so try to add some hkl info too
			wnote = note(FullPeakIndexed)
			Wave recip = str2recip(StringByKey("recip_lattice"+num2istr(ip),wnote,"="))
			if (WaveExists(recip))
				MatrixOp/O/FREE hkl = Inv(recip) x qBL		// go from q in Ideal beam line to hkl (but wrong length)
				sprintf str, "hkl = (%.3f, %.3f, %.3f)",hkl[0]*24,hkl[1]*24,hkl[2]*24
				tagStr += "\r"+str
				Make/N=3/O/D/FREE hkli=hkl						// try to make a nice integral hkl
				CloseAllowedhkl(hkli)
				MatrixOp/O/FREE xyz = recip x hkli
				keV = hc*norm(xyz)/(4*PI*sin(theta))
				hklStr = hkl2str(hkli[0],hkli[1],hkli[2], bar=1)
				sprintf str, "hkl = (%s),  %.4f keV",hklStr,keV
				tagStr += "\r"+str
				if (useFitted)
					MatrixOp/O/FREE hkl = recip x hkli			// go from integral hkl to Q in BL
					Variable dtheta = angleVec2Vec(hkl,qBL)	// angle between fitted peak and indexed peak
					if (useDistortion)
						sprintf str,"\r"+S_DELTA_THETA+" = "+angFmt+",   distort=%.2f px",dtheta,cabs(pz)
					else
						sprintf str,"\r"+S_DELTA_THETA+" = "+angFmt,dtheta
					endif
					tagStr += str
				endif
			endif
		elseif (WaveExists(image))								// test for recip in image wave note
			wnote = note(image)
			Wave recip = str2recip(StringByKey("recip_lattice"+num2istr(ip),wnote,"="))
			if (!WaveExists(recip))
				Wave recip = str2recip(StringByKey("recipRef",wnote,"="))	
			endif
			if (WaveExists(recip) && WaveExists(qBL))
				MatrixOp/O/FREE hkl = Inv(recip) x qBL		// go from q in Ideal beam line to hkl (but wrong length)
				sprintf str, "hkl = (%.3f, %.3f, %.3f)",hkl[0]*24,hkl[1]*24,hkl[2]*24
				tagStr += "\r"+str
				Make/N=3/O/D/FREE hkli=hkl						// try to make a nice integral hkl
				CloseAllowedhkl(hkli)
				MatrixOp/O/FREE xyz = recip x hkli
				keV = hc*norm(xyz)/(4*PI*sin(theta))
				hklStr = hkl2str(hkli[0],hkli[1],hkli[2], bar=1)
				sprintf str, "hkl = (%s),  %.4f keV",hklStr,keV
				tagStr += "\r"+str
			endif
		endif

	elseif (useIndex)													// examining indexed peaks
		if (!WaveExists(FullPeakIndexed))
			return 1
		endif
		dist2=Inf
		m=-1																// find the indexed peak closest to the mouse-click
		N = DimSize(FullPeakIndexed,0)
		Make/N=2/FREE mxmy={mx,my}								// location from mouse
		Make/N=(N,2)/FREE pxpy = FullPeakIndexed[p][9+q][ip] - mxmy[q]
		MatrixOP/FREE dpxpy2 = sumRows(magSqr(pxpy))		// dist^2 from mxmy to all indexed peaks
		if (DimSize(FullPeakIndexed,1)>11)						// only intersested in pixels on this detector
			dpxpy2 = FullPeakIndexed[p][11][ip]==dNum ? dpxpy2[p] : NaN
		endif
		WaveStats/M=1/Q dpxpy2
		m = V_minLoc
		wIndex = m
		TraceName = FullPeakIndexName
		if (m<0)
			return 1														// no indexed peak found
		endif
		h =   FullPeakIndexed[m][3][ip]
		k =   FullPeakIndexed[m][4][ip]
		l =   FullPeakIndexed[m][5][ip]
		px =  FullPeakIndexed[m][9][ip]
		py =  FullPeakIndexed[m][10][ip]
		keV = FullPeakIndexed[m][7][ip]
		angleErr=FullPeakIndexed[m][8][ip]
		SpaceGroupIDnum = NumberByKey("SpaceGroupIDnum",note(FullPeakIndexed),"=")
		if (numtype(SpaceGroupIDnum))
			Variable SG = NumberByKey("SpaceGroup",note(FullPeakIndexed),"=")
			SpaceGroupID = LatticeSym#FindDefaultIDForSG(SG)
			SpaceGroupIDnum = LatticeSym#FindDefaultIDnumForSG(SG)
		endif
		if (strlen(SpaceGroupID)<1)
			SpaceGroupID = IDnum2SpaceGroupID(SpaceGroupIDnum)
		endif

		hklStr = hkl2str(h,k,l, bar=1)
		if (ip>0)
			sprintf str,"hkl=(%s),   %.4f keV\rpixel(%.2f, %.2f),   #%d, %d\rangleErr="+angFmt,hklStr,keV,px,py,m,ip,angleErr
		else
			sprintf str,"hkl=(%s),   %.4f keV\rpixel(%.2f, %.2f),   #%d\rangleErr="+angFmt,hklStr,keV,px,py,m,angleErr
		endif
		tagStr = "\\Zr090Indexed peak position\r" + str
		tagStr += SelectString(numtype(SpaceGroupIDnum),"\r"+getHMsym2(SpaceGroupIDnum)+"    Space Group "+SpaceGroupID,"")

	elseif (useStrain)												// examining strain refined peaks
		if (!WaveExists(PeaksForStrain))
			return 1
		endif
		dist2=Inf
		m=-1																// find the indexed peak closest to the mouse-click
		N = DimSize(PeaksForStrain,0)
		Make/N=2/FREE mxmy={mx,my}								// location from mouse
		Make/N=(N,2)/FREE pxpy = PeaksForStrain[p][11+q] - mxmy[q]
		MatrixOP/FREE dpxpy2 = sumRows(magSqr(pxpy))		// dist^2 from mxmy to all strain refined peaks
		if (DimSize(PeaksForStrain,1)>14)						// only intersested in pixels on this detector
			dpxpy2 = PeaksForStrain[p][14]==dNum ? dpxpy2[p] : NaN
		endif
		WaveStats/M=1/Q dpxpy2
		m = V_minLoc
		wIndex = m
		TraceName = PeaksForStrainName
		dist2 = V_min
		if (m<0)
			return 1														// no indexed peak found
		endif
		h=PeaksForStrain[m][0]
		k=PeaksForStrain[m][1]
		l=PeaksForStrain[m][2]
		px = PeaksForStrain[m][11]
		py = PeaksForStrain[m][12]
		keV=PeaksForStrain[m][10]
		angleErr=PeaksForStrain[m][13]
#ifdef USE_ENERGY_STRAIN_REFINE
		String angleErrUnit = " nm\S-1\M"
#else
		String angleErrUnit = DEGREESIGN
#endif
		hklStr = hkl2str(h,k,l, bar=1)
		sprintf str,"pixel(%.2f, %.2f)\r%.4f keV\r%s=%.2g%s\rhkl=(%s),   #%d",px,py, keV,GDELTA,angleErr,angleErrUnit,hklStr,m
		tagStr = "\\Zr090Strained peak position\r" + str
	endif
	if (WaveExists(image))
		px = limit(px,0,DimSize(image,0)-1)					// needed in case (px,py) is outside the image
		py = limit(py,0,DimSize(image,1)-1)
		wIndex = round(px) + DimSize(image,0)*round(py)	// convert px,py back to index into image
	endif
	GetAxis/W=$win/Q bottom
	horiz = (px-V_min)/ (V_max-V_min)
	GetAxis/W=$win/Q left
	vert = (py-V_max)/(V_min-V_max)
	Variable x0 = horiz<0.5 ? 5 : -5,  y0 = vert<0.5 ? -5 : 5
	String anchor = SelectString(horiz<0.5,"R","L")+ SelectString(vert<0.5,"B","T")
	str = SelectString(strlen(imageName), TraceName, imageName)
	Tag/C/N=indexedPeakInfo/W=$win/A=$anchor/F=2/L=2/X=(x0)/Y=(y0)/P=1 $str, wIndex, tagStr
	DoUpdate
	return 0																// 1 means do not send back to Igor for more processing, 0 means all OK
End


Function ButtonBoxesProc(B_Struct) : ButtonControl
	STRUCT WMButtonAction &B_Struct
	if (B_Struct.eventCode !=2)						// only process on button up
		return 0
	endif
	Wave peakList = $StringByKey("FullPeakList",GetUserData(B_Struct.win,"","FitPeaks"),"=")
	if (!WaveExists(peakList))						// not in user data, search for one
		Wave image = ImageNameToWaveRef(B_Struct.win,StringFromList(0,ImageNameList(B_Struct.win,";")))
		Wave peakList = FindPeakListForImage(image)
	endif
	if (!WaveExists(peakList))						// give up
		DoAlert 0,"Could not find FullPeakList for the image on this Graph"
		return 0
	endif

	String boxKeys = GetUserData(B_Struct.win,"","boxes")	// save for later
	Variable boxesOn = NumberByKey("boxesOn"boxKeys,"=")
	boxesOn = numtype(boxesOn) ? 1 : !boxesOn
	if (B_Struct.eventMod==63)						// special to force ON or OFF
		boxesOn = stringmatch(B_Struct.ctrlName,"ON") ? 1 : boxesOn
		boxesOn = stringmatch(B_Struct.ctrlName,"OFF") ? 0 : boxesOn
	endif
	boxKeys = ReplaceNumberByKey("boxesOn",boxKeys,boxesOn,"=")
	SetWindow $(B_Struct.win), userdata(boxes)=boxKeys

	SetDrawLayer /K UserFront				// always clear first
	if (!boxesOn)								// turn boxes off, really nothing to do
//		String list = GetUserData(B_Struct.win,"","Indexing")
//		SetWindow $(B_Struct.win),userdata(FitPeaks)=RemoveByKey("FullPeakList",list,"=")
	else										// Re-Draw the Boxes
		Variable x0,y0, i, wid
		GetAxis/W=$(B_Struct.win)/Q bottom
		wid = abs(V_max-V_min)			// wid depends upon the size of the graph & how much of the plot is showing
		GetAxis /W=$(B_Struct.win)/Q left
		wid = round(max(wid,abs(V_max-V_min))/50)	
		wid = numtype(wid) ? -1 : wid
		wid = max(wid,2*5)
		SetDrawLayer UserFront
		Variable Npeaks=DimSize(peakList,0)
		if (Npeaks>60)						// reduce size when there are lots of peaks
			wid *= 0.5
		elseif (Npeaks>40)
			wid *= 0.75
		endif
		for (i=0;i<Npeaks;i+=1)
			x0 = peakList[i][0]
			y0 = peakList[i][1]
			if (numtype(x0+y0)==0)
				DrawMarker(x0,y0,wid,wid,"BoxWithTicks")
			endif
		endfor
	endif
	return 0
End

Structure ROIstructure				// structure definition
	double startx						// starting pixel in raw un-binned pixels
	double endx
	double groupx
	double starty
	double endy
	double groupy
EndStructure


Function RemoveMissingPeaksFromGraph(gName)
	String gName					// name of graph to use
	if (strlen(gName)<1)
		gName=WinName(0,1)				// name of top graph
	endif
	String traceName, traceList = TraceNameList(gName,";",1)
	Variable i, N = ItemsInLIst(traceList)
	for (i=0;i<N;i+=1)
		traceName = StringFromLIst(i,traceList)
		if (WaveInClass(TraceNameToWaveRef(gName,traceName),"MissingPeakList*"))
			RemoveFromGraph/W=$gName $traceName
		endif
	endfor
End


Function AddMissingReflections(FullPeakIndexed,[pattern,detector])	// calculate the missing reflections
	Wave FullPeakIndexed
	Variable pattern									// defaults to 0
	Variable detector									// defaults to 0, probably never use this

	String gName=WinName(0,1)					// name of top graph
	Wave image = ImageNameToWaveRef(gName, StringFromLIst(0,ImageNameList(gName,";")) )
	pattern = ParamIsDefault(pattern) || numtype(pattern) || pattern<0 ? NaN : pattern
	if (!(pattern>=0))
		pattern = NumberByKey("patternNum",GetUserData(gName,"","Indexing"),"=")
	endif
	pattern = pattern>=0 ? pattern : 0
	if (!WaveExists(FullPeakIndexed))
		Wave FullPeakIndexed = $StringByKey("FullPeakIndexed",GetUserData(gName,"","Indexing"),"=")
	endif
	if (!WaveExists(FullPeakIndexed) || !WaveExists(image))
		return 1
	elseif (!WaveInClass(FullPeakIndexed,"IndexedPeakList*"))
		return 1
	endif
	String wnote = note(image)

	detector = ParamIsDefault(detector) ? NaN : detector
	if (!(detector>=0))
		STRUCT microGeometry geo
		if (FillGeometryStructDefault(geo, alert=1))	//fill the geometry structure with test values
			return 1
		endif
		detector = detectorNumFromID(geo, StringByKey("detectorID",wnote,"="))
	endif
	detector = detector>=0 ? detector : 0

	STRUCT ROIstructure roi
	roi.startx = NumberByKey("startx",wnote,"=")
	roi.endx = NumberByKey("endx",wnote,"=")
	roi.groupx = NumberByKey("groupx",wnote,"=")
	roi.starty = NumberByKey("starty",wnote,"=")
	roi.endy = NumberByKey("endy",wnote,"=")
	roi.groupy = NumberByKey("groupy",wnote,"=")
	if (numtype(roi.startx + roi.endx + roi.groupx + roi.starty + roi.endy + roi.groupy))
		return 1
	endif

	String indexNote = note(FullPeakIndexed)
	Variable keVmax = NumberByKey("keVmaxTest",indexNote,"=")
	keVmax = numtype(keVmax) ? 30 : keVmax
	STRUCT crystalStructure xtal
	if (setXtalFromNote(indexNote,xtal))
		if (FillCrystalStructDefault(xtal))
			DoAlert 0, "no lattice structure found, did you forget to set it?"
			Abort "in AddMissingReflections()"
		endif
	endif

	Variable NpatternsFound=NumberByKey("NpatternsFound",indexNote,"=")
	NpatternsFound = numtype(NpatternsFound) ? pattern+1 : NpatternsFound
	pattern = limit(round(pattern),0,NpatternsFound-1)
	String str = StringByKey("recip_lattice"+num2istr(pattern),indexNote,"=")
	Make/N=(3,3)/FREE/D recip=NaN
	Variable a0,a1,a2, b0,b1,b2, c0,c1,c2
	sscanf str,"{{%g,%g,%g}{%g,%g,%g}{%g,%g,%g}}",a0,a1,a2, b0,b1,b2, c0,c1,c2
	if (V_flag!=9)
		return 1
	endif
	xtal.as0 = a0;		xtal.bs0 = b0;		xtal.cs0 = c0				// put actual recip into xtal
	xtal.as1 = a1;		xtal.bs1 = b1;		xtal.cs1 = c1
	xtal.as2 = a2;		xtal.bs2 = b2;		xtal.cs2 = c2

	if (NpatternsFound>1 || pattern>0)
		Wave Missing = $CalcMissingReflections(FullPeakIndexed,roi,xtal,keVmax,pattern=pattern,detector=detector)
	else
		Wave Missing = $CalcMissingReflections(FullPeakIndexed,roi,xtal,keVmax,detector=detector)
	endif
	if (!WaveExists(Missing))
		return 1
	endif
	if (strlen(GetRTStackInfo(2))<1 || WhichListItem("hklTagsPopMenuProc",GetRTStackInfo(0))>=0)
		printf "found %d missing reflections\r",DimSize(Missing,0)
	endif

	if (WhichListItem(NameOfWave(Missing),TraceNameList(gName,";",1))<0)		// Missing not on graph
		AppendToGraph/W=$gName Missing[*][4] vs Missing[*][3]
		String mName = NameOfWave(Missing)
		ModifyGraph mode($mName)=3, marker($mName)=8, rgb($mName)=(0,0,65535)
	endif

	return 0
End
//

//	Function/T test()
//		Variable k, k0=5
//		String str="{"
//		for (k=k0;k<=11; k = (k>k0 ? k0-(k-k0) : (k0-k)+k0+1))
//			str += num2istr(k)+", "
//		endfor
//		return str+"}"
//	End
//
//	result is:		  {5, 6, 4, 7, 3, 8, 2, 9, 1, 10, 0, 11, -1, }
//
Static Function/T CalcMissingReflections(FullPeakIndexed,roi,xtal,keVmax,[pattern,detector])	// calculate the missing reflections
	Wave FullPeakIndexed
	STRUCT ROIstructure &roi
	STRUCT crystalStructure &xtal
	Variable keVmax
	Variable pattern									// defaults to 0
	Variable detector									// defaults to 0, probably never use this
	pattern = ParamIsDefault(pattern) || numtype(pattern) ? 0 : round(pattern)
	detector = ParamIsDefault(detector) ? 0 : round(detector)
	if (!(pattern>=0) || !(detector>=0) || !WaveInClass(FullPeakIndexed,"IndexedPeakList*"))
		return ""
	endif
	Variable is4500s=NumVarOrDefault(MICRO_GEOMETRY_VERSION_PATH,0)&2
	Variable keVmin = is4500s ? 3 : 7
	if (!(keVmin<keVmax))
		return ""
	endif

	Wave recip = recipFrom_xtal(xtal)			// returns a FREE wave with reciprocal lattice

	STRUCT microGeometry g	
	FillGeometryStructDefault(g)					//fill the geometry structure with current values

	Variable xlo, xhi,ylo,yhi
	xlo = (roi.startx-(roi.startx-FIRST_PIXEL)-(roi.groupx-1)/2)/roi.groupx		// get range of this ROI in binned pixels
	ylo = (roi.starty-(roi.starty-FIRST_PIXEL)-(roi.groupy-1)/2)/roi.groupy
	xhi = (roi.endx-(roi.startx-FIRST_PIXEL)-(roi.groupx-1)/2)/roi.groupx
	yhi = (roi.endy-(roi.starty-FIRST_PIXEL)-(roi.groupy-1)/2)/roi.groupy
	Variable depth = NumberByKey("depth",note(FullPeakIndexed),"=")

	// find Qvector to center of roi, and max cone angle
	Variable px=(roi.startx+roi.endx)/2, py=(roi.starty+roi.endy)/2
	Make/N=3/D/FREE q00, qC, q00,q10,q01,q11			// q vectors to center of detector and each of its corners
	pixel2q(g.d[detector],px,py,qC,depth=depth)							// find q vectors to center and all 4 corners too
	pixel2q(g.d[detector],roi.startx,roi.starty,q00,depth=depth)
	pixel2q(g.d[detector],roi.endx,roi.starty,q10,depth=depth)
	pixel2q(g.d[detector],roi.startx,roi.endy,q01,depth=depth)
	pixel2q(g.d[detector],roi.endx,roi.endy,q11,depth=depth)

	Variable hlo,hhi, klo,khi, llo,lhi, Qmag
	Qmag = 4*PI*keVMin/hc * (-qC[2])					// since ki={0,0,1}, MatrixDot(ki,qC) = qC[2]
	MatrixOP/FREE/O hkl = (Inv(recip) x qC) * Qmag
	hlo = hkl[0]	;	klo = hkl[1]	;	llo = hkl[2]
	hhi = hlo		;	khi = klo		;	lhi = llo
	Qmag = 4*PI*keVMax/hc * (-qC[2])
	MatrixOP/FREE/O hkl = (Inv(recip) x qC) * Qmag
	hlo = min(hkl[0],hlo)	;	klo = min(hkl[1],klo)	;	llo = min(hkl[2],llo)
	hhi = max(hkl[0],hhi)	;	khi = max(hkl[1],khi)	;	lhi = max(hkl[2],lhi)

	Qmag = 4*PI*keVMin/hc * (-q00[2])					// since ki={0,0,1}, MatrixDot(ki,q00) = q00[2]
	MatrixOP/FREE/O hkl = (Inv(recip) x q00) * Qmag
	hlo = min(hkl[0],hlo)	;	klo = min(hkl[1],klo)	;	llo = min(hkl[2],llo)
	hhi = max(hkl[0],hhi)	;	khi = max(hkl[1],khi)	;	lhi = max(hkl[2],lhi)
	Qmag = 4*PI*keVMax/hc * (-q00[2])
	MatrixOP/FREE/O hkl = (Inv(recip) x q00) * Qmag
	hlo = min(hkl[0],hlo)	;	klo = min(hkl[1],klo)	;	llo = min(hkl[2],llo)
	hhi = max(hkl[0],hhi)	;	khi = max(hkl[1],khi)	;	lhi = max(hkl[2],lhi)

	Qmag = 4*PI*keVMin/hc * (-q10[2])
	MatrixOP/FREE/O hkl = (Inv(recip) x q10) * Qmag
	hlo = min(hkl[0],hlo)	;	klo = min(hkl[1],klo)	;	llo = min(hkl[2],llo)
	hhi = max(hkl[0],hhi)	;	khi = max(hkl[1],khi)	;	lhi = max(hkl[2],lhi)
	Qmag = 4*PI*keVMax/hc * (-q10[2])
	MatrixOP/FREE/O hkl = (Inv(recip) x q10) * Qmag
	hlo = min(hkl[0],hlo)	;	klo = min(hkl[1],klo)	;	llo = min(hkl[2],llo)
	hhi = max(hkl[0],hhi)	;	khi = max(hkl[1],khi)	;	lhi = max(hkl[2],lhi)

	Qmag = 4*PI*keVMin/hc * (-q01[2])
	MatrixOP/FREE/O hkl = (Inv(recip) x q01) * Qmag
	hlo = min(hkl[0],hlo)	;	klo = min(hkl[1],klo)	;	llo = min(hkl[2],llo)
	hhi = max(hkl[0],hhi)	;	khi = max(hkl[1],khi)	;	lhi = max(hkl[2],lhi)
	Qmag = 4*PI*keVMax/hc * (-q01[2])
	MatrixOP/FREE/O hkl = (Inv(recip) x q01) * Qmag
	hlo = min(hkl[0],hlo)	;	klo = min(hkl[1],klo)	;	llo = min(hkl[2],llo)
	hhi = max(hkl[0],hhi)	;	khi = max(hkl[1],khi)	;	lhi = max(hkl[2],lhi)

	Qmag = 4*PI*keVMin/hc * (-q11[2])
	MatrixOP/FREE/O hkl = (Inv(recip) x q11) * Qmag
	hlo = min(hkl[0],hlo)	;	klo = min(hkl[1],klo)	;	llo = min(hkl[2],llo)
	hhi = max(hkl[0],hhi)	;	khi = max(hkl[1],khi)	;	lhi = max(hkl[2],lhi)
	Qmag = 4*PI*keVMax/hc * (-q11[2])
	MatrixOP/FREE/O hkl = (Inv(recip) x q11) * Qmag
	hlo = min(hkl[0],hlo)	;	klo = min(hkl[1],klo)	;	llo = min(hkl[2],llo)
	hhi = max(hkl[0],hhi)	;	khi = max(hkl[1],khi)	;	lhi = max(hkl[2],lhi)

	hlo=floor(hlo)			;	hhi = ceil(hhi)
	klo=floor(klo)			;	khi = ceil(khi)
	llo=floor(llo)			;	lhi = ceil(lhi)

	Variable m, Nm=100
	String mName="Missing"+num2istr(detector)
	Make/N=(Nm,6)/O $mName=NaN
	Wave Missing=$mName
	SetScale d,0,1,"pixel",Missing
	SetDimLabel 1,0,H,Missing		;	SetDimLabel 1,1,K,Missing		;	SetDimLabel 1,2,L,Missing
	SetDimLabel 1,3,px,Missing		;	SetDimLabel 1,4,py,Missing
	SetDimLabel 1,5,keV,Missing

	String str, noteStr=ReplaceStringByKey("waveClass","","MissingPeakList","=")
	sprintf str,"%d,%d",hlo,hhi
	noteStr=ReplaceStringByKey("Hrange",noteStr,str,"=")
	sprintf str,"%d,%d",klo,khi
	noteStr=ReplaceStringByKey("Krange",noteStr,str,"=")
	sprintf str,"%d,%d",llo,lhi
	noteStr=ReplaceStringByKey("Lrange",noteStr,str,"=")
	noteStr=ReplaceNumberByKey("keVmin",noteStr,keVmin,"=")
	noteStr=ReplaceNumberByKey("keVmax",noteStr,keVmax,"=")
	if (!ParamIsDefault(pattern))
		noteStr=ReplaceNumberByKey("pattern",noteStr,pattern,"=")
	endif
	if (!ParamIsDefault(detector))
		noteStr=ReplaceNumberByKey("detectorNum",noteStr,detector,"=")
	endif
	if (numtype(depth)==0)
		noteStr=ReplaceNumberByKey("depth",noteStr,depth,"=")
	endif
	str = ReplaceString(";",xtal.desc, "_")
	str = ReplaceString("=",xtal.desc, "_")
	noteStr=ReplaceStringByKey("xtalDesc", noteStr, str,"=")
	noteStr=ReplaceStringByKey("recip", noteStr, encodeMatAsStr(recip),"=")
	Note/K Missing, noteStr

	Make/N=3/FREE/D hkli,hklt
	Variable/C pz
	Variable h,k,l,  keV, theta
	Variable i,Ni=DimSize(FullPeakIndexed,0)
	for (l=llo; l<=lhi; l+=1)											// loop symmetrically outward from (h0,k0,l0)
		for (k=klo; k<=khi; k+=1)
			for (h=hlo; h<=hhi; h+=1)
				hklt = {h,k,l}
				for (i=0;i<Ni;i+=1)									// reject if already in FullPeakIndexed
					hkli = FullPeakIndexed[i][p+3][pattern]
					if (MatrixDot(hkli,hklt)/(norm(hklt)*norm(hkli))>0.99999)	// skip if hkl are parallel
						break
					endif
				endfor
				if (i<Ni)												// found hkl in FullPeakIndexed, so skip this and continue
					continue
				endif
				MatrixOP/FREE/O qvec = recip x hklt					// find point on detector
				Qmag = norm(qvec)
				pz = q2pixel(g.d[detector],qvec)	
				px = (real(pz)-(roi.startx-FIRST_PIXEL)-(roi.groupx-1)/2)/roi.groupx	// binned pixels
				py = (imag(pz)-(roi.starty-FIRST_PIXEL)-(roi.groupy-1)/2)/roi.groupy
				if (!(xlo<px && px<xhi && ylo<py && py<yhi))			// must land on detector
					continue
				endif
				theta = pixel2q(g.d[detector],real(pz),imag(pz),qvec,depth=depth)	// get theta to calc energy
				keV = Qmag * hc / (4*PI*sin(theta))					// energy(keV)
				if (!(keVmin<keV && keV<keVmax))
					continue
				endif
				if (!allowedHKL(h,k,l,xtal))							// reject forbidden reflections
					continue
				endif

				// check for duplicates
				for (i=0;i<m;i+=1)
					if ((abs(px-Missing[i][3])+abs(py-Missing[i][4]))<1)
						break
					endif
				endfor
				if (i<m)
					continue											// skip duplicate location
				endif

				if (m>=Nm)											// extend wave if necessary
					Nm += 100
					Redimension/N=(Nm,-1) Missing
				endif
				Missing[m][0,2] = hklt[q]
				Missing[m][3] = px
				Missing[m][4] = py
				Missing[m][5] =keV									// energy(keV)
				m += 1
			endfor
		endfor
	endfor
	Nm = m
	Redimension/N=(Nm,-1) Missing
	return GetWavesDataFolder(Missing,2)
End
//
Static Function setXtalFromNote(wnote,xtal)
	String wnote
	STRUCT crystalStructure &xtal

	String structureDesc = StringByKey("structureDesc",wnote,"=")
	structureDesc = structureDesc[0,99]
	Variable SpaceGroup = NumberByKey("SpaceGroup",wnote,"=")
	if (!(1<=SpaceGroup && SpaceGroup<=230))
		return 1
	endif
	String SpaceGroupID = StringByKey("SpaceGroupID",wnote,"=")
	Variable SpaceGroupIDnum = NumberByKey("SpaceGroupIDnum",wnote,"=")
	if (strlen(SpaceGroupID)<1)
		SpaceGroupID = LatticeSym#FindDefaultIDForSG(SpaceGroup)
	endif
	SpaceGroupIDnum = numtype(SpaceGroupIDnum) ? str2num(SpaceGroupID) : SpaceGroupIDnum

	String str = StringByKey("latticeParameters",wnote,"=")
	Variable a,b,c,alpha,bet,gam
	sscanf str, "{%g, %g, %g, %g, %g, %g}",a,b,c,alpha,bet,gam
	if (V_flag!=6)
		return 1
	endif

	str = StringByKey("lengthUnit",wnote,"=")
	if (stringmatch(str,"Å") || stringmatch(str,"Angstrom)"))
		a /= 10
		b /= 10
		c /= 10
	endif
	xtal.a = a
	xtal.b = b
	xtal.c = c
	xtal.alpha = alpha
	xtal.beta = bet
	xtal.gam = gam
	xtal.SpaceGroup = SpaceGroup
	xtal.SpaceGroupID = SpaceGroupID
	xtal.SpaceGroupIDnum = SpaceGroupIDnum
	xtal.desc = structureDesc
	xtal.Unconventional00 = NaN

	Variable xx,yy,zz,occ,Zatom
	Variable N
	String atom
	for (N=0;N<20;N+=1)
		str = StringByKey("AtomDesctiption"+num2istr(N+1),wnote,"=")
		if (strlen(str)<2)
			break
		endif
		sscanf str, "{%s  %g %g %g %g}",atom,xx,yy,zz,occ
		Zatom = LatticeSym#ZfromLabel(atom)
		if (numtype(xx+yy+zz+Zatom) || V_flag<4)
			break
		endif
		xtal.atom[N].name = atom
		xtal.atom[N].Zatom = LatticeSym#ZfromLabel(atom)
		xtal.atom[N].x = xx
		xtal.atom[N].y = yy
		xtal.atom[N].z = zz
		if (V_flag==5)
			xtal.atom[N].occ = numtype(occ) ? 1 : occ
		endif
	endfor
	xtal.N = N
	LatticeSym#setDirectRecip(xtal)	
	return 0
End





// supprot for adding hkl from a user supplied reicprocal lattice
Function RemoveUserRecipFromGraph(gName)
	String gName					// name of graph to use
	if (strlen(gName)<1)
		gName=WinName(0,1)				// name of top graph
	endif
	String traceName, traceList = TraceNameList(gName,";",1)
	Variable i, N = ItemsInLIst(traceList)
	for (i=0;i<N;i+=1)
		traceName = StringFromLIst(i,traceList)
		if (WaveInClass(TraceNameToWaveRef(gName,traceName),"UserRecipPeakList*"))
			RemoveFromGraph/W=$gName $traceName
		endif
	endfor
End
//
Function AddReflectionsFromUserRecip(recip,[keVmax,detector])	// calculate reflections based on user supplied recip[3][3]
	Wave recip										// a 3x3 matrix of the reciprocal lattice
	Variable keVmax
	Variable detector								// defaults to 0, probably never use this
	keVmax = ParamIsDefault(keVmax) || keVmax<=0 ? NaN : keVmax
	detector = ParamIsDefault(detector) ? NaN : detector

	if (!WaveExists(recip) || WaveDims(recip)!=2 || DimSize(recip,0)!=3 || DimSize(recip,1)!=3)
		String rname, rlist=WaveList("*",";", "DIMS:2,MINROWS:3,MAXROWS:3,MINCOLS:3,MAXCOLS:3")
		if (ItemsInList(rlist)<1)
			return 1
		elseif (ItemsInList(rlist)==1)
			Wave recip = $StringFromList(0,rlist)
		else
			Prompt rname, "Pick reciprocal lattice", popup, WaveList("*",";", "DIMS:2,MINROWS:3,MAXROWS:3,MINCOLS:3,MAXCOLS:3")
			DoPrompt "Reciprocal Lattice", rname
			if (V_flag)
				return 1
			endif
			Wave recip = $rname
		endif
	endif
	if (!WaveExists(recip) || WaveDims(recip)!=2 || DimSize(recip,0)!=3 || DimSize(recip,1)!=3)
		return 1
	endif

	String gName=WinName(0,1)					// name of top graph
	Wave image = ImageNameToWaveRef(gName, StringFromLIst(0,ImageNameList(gName,";")) )
	if (!WaveExists(image))
		return 1
	endif
	String wnote = note(image)

	if (!(detector>=0))
		STRUCT microGeometry geo
		if (FillGeometryStructDefault(geo, alert=1))	//fill the geometry structure with test values
			return 1
		endif
		detector = detectorNumFromID(geo, StringByKey("detectorID",wnote,"="))
	endif
	detector = detector>=0 ? detector : 0

	STRUCT ROIstructure roi
	roi.startx = NumberByKey("startx",wnote,"=")
	roi.endx = NumberByKey("endx",wnote,"=")
	roi.groupx = NumberByKey("groupx",wnote,"=")
	roi.starty = NumberByKey("starty",wnote,"=")
	roi.endy = NumberByKey("endy",wnote,"=")
	roi.groupy = NumberByKey("groupy",wnote,"=")
	if (numtype(roi.startx + roi.endx + roi.groupx + roi.starty + roi.endy + roi.groupy))
		return 1
	endif

	if (numtype(keVmax) || keVmax<=0)						// valid keVmax was not supplied, check FullPeakIndexed[][]
		Wave FullPeakIndexed = $StringByKey("FullPeakIndexed",GetUserData(gName,"","Indexing"),"=")
		if (WaveExists(FullPeakIndexed))
			keVmax = NumberByKey("keVmaxTest",note(FullPeakIndexed),"=")
		endif
	endif
	if (numtype(keVmax) || keVmax<=0)						// still no valid keVmax, ask user
		keVmax = 25
		Prompt keVmax, "Upper cut-off energy (keV)"
		DoPrompt "cut-off energy", keVmax
		if (V_flag)
			return 1
		endif
	endif
	keVmax = numtype(keVmax) ? 30 : keVmax

	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))
		DoAlert 0, "ERROR AddReflectionsFromUserRecip()\r   no lattice structure found,\r   did you forget to set it?"
		return 1
	endif
	xtal.as0 = recip[0][0];		xtal.bs0 = recip[0][1];		xtal.cs0 = recip[0][2]				// put actual recip into xtal
	xtal.as1 = recip[1][0];		xtal.bs1 = recip[1][1];		xtal.cs1 = recip[1][2]
	xtal.as2 = recip[2][0];		xtal.bs2 = recip[2][1];		xtal.cs2 = recip[2][2]

	String mName = "hkl_" + NameOfWave(recip) + SelectString(detector, "", "_"+num2istr(detector))
	mName = CleanupName(mName,0)
	Wave userHKL = CalcReflections(roi,xtal,keVmax,detector=detector, name=mName)
	if (!WaveExists(userHKL))
		return 1
	endif
	if (strlen(GetRTStackInfo(2))<1 || WhichListItem("hklTagsPopMenuProc",GetRTStackInfo(0))>=0)
		printf "found %d reflections from %s,  created wave: \"%s\"\r",DimSize(userHKL,0),NameOfWave(recip), mName
	endif

	if (WhichListItem(NameOfWave(userHKL),TraceNameList(gName,";",1))<0)		// userHKL not on graph
		AppendToGraph/W=$gName userHKL[*][4] vs userHKL[*][3]
		ModifyGraph mode($mName)=3, marker($mName)=7, rgb($mName)=(0,52000,0)
	endif

	return 0
End
//
Static Function/WAVE CalcReflections(roi,xtal,keVmax,[detector,depth,name])	// calculate the missing reflections
	STRUCT ROIstructure &roi
	STRUCT crystalStructure &xtal
	Variable keVmax
	Variable detector									// defaults to 0, probably never use this
	Variable depth
	String name											// name to use for output wave, may be empty
	detector = ParamIsDefault(detector) ? 0 : round(detector)
	depth = ParamIsDefault(depth) || numtype(depth) ? NaN : depth
	name = SelectString(ParamIsDefault(name) , name, "")
	Variable is4500s=NumVarOrDefault(MICRO_GEOMETRY_VERSION_PATH,0)&2
	Variable keVmin = is4500s ? 3 : 7
	if (!(keVmin<keVmax))
		return $""
	endif

	Wave recip = recipFrom_xtal(xtal)				// returns a FREE wave with reciprocal lattice
	STRUCT microGeometry g	
	FillGeometryStructDefault(g)					//fill the geometry structure with current values

	Variable xlo, xhi,ylo,yhi
	xlo = (roi.startx-(roi.startx-FIRST_PIXEL)-(roi.groupx-1)/2)/roi.groupx		// get range of this ROI in binned pixels
	ylo = (roi.starty-(roi.starty-FIRST_PIXEL)-(roi.groupy-1)/2)/roi.groupy
	xhi = (roi.endx-(roi.startx-FIRST_PIXEL)-(roi.groupx-1)/2)/roi.groupx
	yhi = (roi.endy-(roi.starty-FIRST_PIXEL)-(roi.groupy-1)/2)/roi.groupy

	// find Qvector to center of roi, and max cone angle
	Variable px=(roi.startx+roi.endx)/2, py=(roi.starty+roi.endy)/2
	Make/N=3/D/FREE q00, qC, q00,q10,q01,q11			// q vectors to center of detector and each of its corners
	pixel2q(g.d[detector],px,py,qC,depth=depth)							// find q vectors to center and all 4 corners too
	pixel2q(g.d[detector],roi.startx,roi.starty,q00,depth=depth)
	pixel2q(g.d[detector],roi.endx,roi.starty,q10,depth=depth)
	pixel2q(g.d[detector],roi.startx,roi.endy,q01,depth=depth)
	pixel2q(g.d[detector],roi.endx,roi.endy,q11,depth=depth)

	Variable hlo,hhi, klo,khi, llo,lhi, Qmag
	Qmag = 4*PI*keVMin/hc * (-qC[2])					// since ki={0,0,1}, MatrixDot(ki,qC) = qC[2]
	MatrixOP/FREE/O hkl = (Inv(recip) x qC) * Qmag
	hlo = hkl[0]	;	klo = hkl[1]	;	llo = hkl[2]
	hhi = hlo		;	khi = klo		;	lhi = llo
	Qmag = 4*PI*keVMax/hc * (-qC[2])
	MatrixOP/FREE/O hkl = (Inv(recip) x qC) * Qmag
	hlo = min(hkl[0],hlo)	;	klo = min(hkl[1],klo)	;	llo = min(hkl[2],llo)
	hhi = max(hkl[0],hhi)	;	khi = max(hkl[1],khi)	;	lhi = max(hkl[2],lhi)

	Qmag = 4*PI*keVMin/hc * (-q00[2])					// since ki={0,0,1}, MatrixDot(ki,q00) = q00[2]
	MatrixOP/FREE/O hkl = (Inv(recip) x q00) * Qmag
	hlo = min(hkl[0],hlo)	;	klo = min(hkl[1],klo)	;	llo = min(hkl[2],llo)
	hhi = max(hkl[0],hhi)	;	khi = max(hkl[1],khi)	;	lhi = max(hkl[2],lhi)
	Qmag = 4*PI*keVMax/hc * (-q00[2])
	MatrixOP/FREE/O hkl = (Inv(recip) x q00) * Qmag
	hlo = min(hkl[0],hlo)	;	klo = min(hkl[1],klo)	;	llo = min(hkl[2],llo)
	hhi = max(hkl[0],hhi)	;	khi = max(hkl[1],khi)	;	lhi = max(hkl[2],lhi)

	Qmag = 4*PI*keVMin/hc * (-q10[2])
	MatrixOP/FREE/O hkl = (Inv(recip) x q10) * Qmag
	hlo = min(hkl[0],hlo)	;	klo = min(hkl[1],klo)	;	llo = min(hkl[2],llo)
	hhi = max(hkl[0],hhi)	;	khi = max(hkl[1],khi)	;	lhi = max(hkl[2],lhi)
	Qmag = 4*PI*keVMax/hc * (-q10[2])
	MatrixOP/FREE/O hkl = (Inv(recip) x q10) * Qmag
	hlo = min(hkl[0],hlo)	;	klo = min(hkl[1],klo)	;	llo = min(hkl[2],llo)
	hhi = max(hkl[0],hhi)	;	khi = max(hkl[1],khi)	;	lhi = max(hkl[2],lhi)

	Qmag = 4*PI*keVMin/hc * (-q01[2])
	MatrixOP/FREE/O hkl = (Inv(recip) x q01) * Qmag
	hlo = min(hkl[0],hlo)	;	klo = min(hkl[1],klo)	;	llo = min(hkl[2],llo)
	hhi = max(hkl[0],hhi)	;	khi = max(hkl[1],khi)	;	lhi = max(hkl[2],lhi)
	Qmag = 4*PI*keVMax/hc * (-q01[2])
	MatrixOP/FREE/O hkl = (Inv(recip) x q01) * Qmag
	hlo = min(hkl[0],hlo)	;	klo = min(hkl[1],klo)	;	llo = min(hkl[2],llo)
	hhi = max(hkl[0],hhi)	;	khi = max(hkl[1],khi)	;	lhi = max(hkl[2],lhi)

	Qmag = 4*PI*keVMin/hc * (-q11[2])
	MatrixOP/FREE/O hkl = (Inv(recip) x q11) * Qmag
	hlo = min(hkl[0],hlo)	;	klo = min(hkl[1],klo)	;	llo = min(hkl[2],llo)
	hhi = max(hkl[0],hhi)	;	khi = max(hkl[1],khi)	;	lhi = max(hkl[2],lhi)
	Qmag = 4*PI*keVMax/hc * (-q11[2])
	MatrixOP/FREE/O hkl = (Inv(recip) x q11) * Qmag
	hlo = min(hkl[0],hlo)	;	klo = min(hkl[1],klo)	;	llo = min(hkl[2],llo)
	hhi = max(hkl[0],hhi)	;	khi = max(hkl[1],khi)	;	lhi = max(hkl[2],lhi)

	hlo=floor(hlo)			;	hhi = ceil(hhi)
	klo=floor(klo)			;	khi = ceil(khi)
	llo=floor(llo)			;	lhi = ceil(lhi)

	Variable m, Nm=100

	if (strlen(name)<1)
		name = "hkl_"+num2istr(detector)
	endif
	Make/N=(Nm,6)/O $name=NaN
	Wave userHKL=$name
	SetScale d,0,1,"pixel",userHKL
	SetDimLabel 1,0,H,userHKL		;	SetDimLabel 1,1,K,userHKL		;	SetDimLabel 1,2,L,userHKL
	SetDimLabel 1,3,px,userHKL		;	SetDimLabel 1,4,py,userHKL
	SetDimLabel 1,5,keV,userHKL

	String str, noteStr=ReplaceStringByKey("waveClass","","UserRecipPeakList","=")
	sprintf str,"%d,%d",hlo,hhi
	noteStr=ReplaceStringByKey("Hrange",noteStr,str,"=")
	sprintf str,"%d,%d",klo,khi
	noteStr=ReplaceStringByKey("Krange",noteStr,str,"=")
	sprintf str,"%d,%d",llo,lhi
	noteStr=ReplaceStringByKey("Lrange",noteStr,str,"=")
	noteStr=ReplaceNumberByKey("keVmin",noteStr,keVmin,"=")
	noteStr=ReplaceNumberByKey("keVmax",noteStr,keVmax,"=")
	if (!ParamIsDefault(detector))
		noteStr=ReplaceNumberByKey("detectorNum",noteStr,detector,"=")
	endif
	if (numtype(depth)==0)
		noteStr=ReplaceNumberByKey("depth",noteStr,depth,"=")
	endif
	str = ReplaceString(";",xtal.desc, "_")
	str = ReplaceString("=",xtal.desc, "_")
	noteStr=ReplaceStringByKey("xtalDesc", noteStr, str,"=")
	noteStr=ReplaceStringByKey("recip", noteStr, encodeMatAsStr(recip),"=")
	Note/K userHKL, noteStr

	Make/N=3/FREE/D hkli,hklt
	Variable/C pz
	Variable h,k,l,  keV, theta
	Variable i
	for (l=llo; l<=lhi; l+=1)										// loop symmetrically outward from (h0,k0,l0)
		for (k=klo; k<=khi; k+=1)
			for (h=hlo; h<=hhi; h+=1)
				hklt = {h,k,l}
				MatrixOP/FREE/O qvec = recip x hklt				// find point on detector
				Qmag = norm(qvec)
				pz = q2pixel(g.d[detector],qvec)	
				px = (real(pz)-(roi.startx-FIRST_PIXEL)-(roi.groupx-1)/2)/roi.groupx	// binned pixels
				py = (imag(pz)-(roi.starty-FIRST_PIXEL)-(roi.groupy-1)/2)/roi.groupy
				if (!(xlo<px && px<xhi && ylo<py && py<yhi))	// must land on detector
					continue
				endif
				theta = pixel2q(g.d[detector],real(pz),imag(pz),qvec,depth=depth)	// get theta to calc energy
				keV = Qmag * hc / (4*PI*sin(theta))				// energy(keV)
				if (!(keVmin<keV && keV<keVmax))
					continue
				endif
				if (!allowedHKL(h,k,l,xtal))						// reject forbidden reflections
					continue
				endif

				// check for duplicates
				for (i=0;i<m;i+=1)
					if ((abs(px-userHKL[i][3])+abs(py-userHKL[i][4]))<1)
						break
					endif
				endfor
				if (i<m)
					continue												// skip duplicate location
				endif

				if (m>=Nm)												// extend wave if necessary
					Nm += 100
					Redimension/N=(Nm,-1) userHKL
				endif
				userHKL[m][0,2] = hklt[q]
				userHKL[m][3] = px
				userHKL[m][4] = py
				userHKL[m][5] =keV										// energy(keV)
				m += 1
			endfor
		endfor
	endfor
	Nm = m
	Redimension/N=(Nm,-1) userHKL
	return userHKL
End



Function/WAVE FindPeakListForImage(image)	// returns wave ref of FullPeakList for an image
	Wave image
	String wList=WaveListClass("FittedPeakList*","*","MINCOLS:1")
	Variable i, N=ItemsInList(wList)
	for (i=0;i<N;i+=1)
		Wave peakList = $StringFromList(i,wList)
		Wave iTest = $StringByKey("fittedIgorImage",note(peakList),"=")
		if (WaveRefsEqual(image,iTest))
			return peakList
		endif
	endfor
	return $""
End

Function/WAVE FindIndexListForImage(image)	// returns wave ref of FullPeakIndex for an image
	Wave image
	String wList=WaveListClass("IndexedPeakList*","*","MAXCOLS:12"), ilist
	Variable m, i, N=ItemsInList(wList)
	for (i=0;i<N;i+=1)
		Wave indexList = $StringFromList(i,wList)
		ilist = StringByKey("fittedIgorImage",note(indexList),"=")		// this can now be a list (comma separated)
		for (m=0;m<ItemsInLIst(ilist,",");m+=1)
			Wave iTest = $StringFromList(m,ilist,",")
			if (WaveRefsEqual(image,iTest))
				return indexList
			endif
		endfor
	endfor
	return $""
End
//Function/WAVE FindPeakListForImage(image)
//	Wave image
//	if (!WaveExists(image))
//		return ""
//	endif
//	String fullName = GetWavesDataFolder(image,2)
//
//	String wlist = WaveListClass("FittedPeakList*","*","MINCOLS:1")
//	Variable i
//	for (i=0;i<ItemsInList(wlist);i+=1)
//		Wave peakList = $StringFromList(i,wlist)
//		if (stringmatch(fullName,StringByKey("fittedIgorImage", note(peakList),"=")))
//			return GetWavesDataFolder(peakList,2)
//		endif
//	endfor
//	return ""
//End
//
//Function/T FindIndexListForImage(image)
//	Wave image
//	if (!WaveExists(image))
//		return ""
//	endif
//	String fullName = GetWavesDataFolder(image,2)
//
//	String wlist = WaveListClass("IndexedPeakList*","*","MAXCOLS:12")
//	Variable i
//	for (i=0;i<ItemsInList(wlist);i+=1)
//		Wave indexList = $StringFromList(i,wlist)
//		if (stringmatch(fullName,StringByKey("fittedIgorImage", note(indexList),"=")))
//			return GetWavesDataFolder(indexList,2)
//		endif
//	endfor
//	return ""
//End


// make an empty FittedPeakList for the image.  This is particularly useful when selecting "AddGaussianToPeakList"
//	before any peak fitting.  When calling from AddGaussianToPeakLIst, use keyVals="FittedPeakShape=Gaussian".
Function/WAVE MakeEmptyPeakListForImage(image,[ask,keyVals])
	Wave image
	Variable ask							// if True, will ask user before creating new FullPeakList, False just creates.
	String keyVals							// additional "key=value;" pairs to add to wave note
	if (!WaveExists(image))
		return $""
	endif

	Wave peakListWave=FindPeakListForImage(image)
	if (WaveExists(peakListWave))
		DoAlert 0, "FullPeakList already exists, will use it."
		return peakListWave
	endif
	if (ParamIsDefault(keyVals))
		keyVals = ""
	endif
	ask = ParamIsDefault(ask) ? 0 : !(ask==0)
	if (ask)									// ask user whether to create new FullPeakList
		DoAlert 1,"There is no peak list for '"+NameOfWave(image)+"', create one?"
		if (V_flag!=1)
			return $""
		endif
	endif

	String wnote = note(image)
	String peakListName="FullPeakList"+detectorID2color(StringByKey("detectorID",wnote,"="))
	Make/N=(0,12)/O $peakListName/WAVE=FullPeakList
	wnote = ReplaceStringByKey("fittedIgorImage",wnote,GetWavesDataFolder(image,2),"=")
	wnote = ReplaceStringByKey("waveClass",wnote,"FittedPeakList","=")
	String item,key,value
	Variable i
	for (i=0;i<ItemsInList(keyVals);i+=1)	// add any extra keyVals that were passed
		item = StringFromList(0,keyVals)
		key = StringFromList(0,item,"=")
		value = StringFromList(1,item,"=")
		if (strlen(item)>0 && strlen(value)>0)
			wnote = ReplaceStringByKey(key,wnote,value,"=")
		endif
	endfor
	Note/K FullPeakList,wnote
	SetDimLabel 1,0,x0,FullPeakList			;	SetDimLabel 1,1,y0,FullPeakList
	SetDimLabel 1,2,x0Err,FullPeakList		;	SetDimLabel 1,3,y0Err,FullPeakList
	SetDimLabel 1,4,fwx,FullPeakList		;	SetDimLabel 1,5,fwy,FullPeakList
	SetDimLabel 1,6,fwxErr,FullPeakList	;	SetDimLabel 1,7,fwyErr,FullPeakList
	SetDimLabel 1,8,correlation,FullPeakList ;	SetDimLabel 1,9,correlationErr,FullPeakList
	SetDimLabel 1,10,area,FullPeakList		;	SetDimLabel 1,11,amp,FullPeakList
	return FullPeakList
End


Function MakeIndexingReportSheet()	// make a report of the indexation using the top graph
	String gName = StringFromList(0,WinList("*",";","WIN:1"))
	String list = GetUserData(gName,"","Indexing")
	Wave FullPeakIndexed = $(StringByKey("FullPeakIndexed",list,"="))
	Variable ip = NumberByKey("patternNum",list,"=")
	ip = numtype(ip) || ip<0 ? 0 : ip
	if (!WaveExists(FullPeakIndexed))
		DoAlert 0,"cannot find indexed peaks on the top graph"
		return 1
	endif

	String imageName = StringFromList(0,ImageNameList(gName,";"))
	Wave image = ImageNameToWaveRef(gName,imageName)
	String textOut="\\Zr150"+imageName+SelectString(ip,"","  pattern "+num2istr(ip))
	if (WaveExists(image))
		String dateExposed="", imageFileName=""
		dateExposed = StringByKey("dateExposed",note(image),"=")
		imageFileName = StringByKey("imageFileName",note(image),"=")
		textOut += "\\Zr075"
		textOut += SelectString(strlen(imageFileName),"","   from file '"+imageFileName+"'")
		textOut += SelectString(strlen(dateExposed),"","   taken on  "+dateExposed)
	endif

//	Wave FullPeakIndexed=$(GetWavesDataFolder(image,1)+"FullPeakIndexed")
	String str = ParagraphsFromIndexing(FullPeakIndexed,ip)
	String par1 = StringFromLIst(0,str), par2 = StringFromLIst(1,str)
	NewLayout/C=1/K=1/P=Portrait
	AppendLayoutObject/F=0/R=(43,23,582,506) graph $gName
	if (strlen(par1))
		TextBox/N=text0/C/F=0/T={70,110,210}/A=LT/X=6/Y=69.4 par1
	endif
	if (strlen(par2))
		TextBox/N=text1/C/F=0/T={70,110,210}/A=RT/X=6/Y=69.4 par2
	endif
	TextBox/N=text2/F=0/X=4.70/Y=66.80 textOut
	str = StringByKey("structureDesc",note(FullPeakIndexed),"=")
	if (strlen(str))
		TextBox/C/N=textDesc/B=1/F=0/A=RT/X=5/Y=0 "\\Zr150"+str
	endif
	if (countChars(par2,"\r\n") <= 12)			// add reciprocal lattice if there is room
		str = "recip_lattice" + num2istr(ip)
		Wave recip = decodeMatFromStr(StringByKey(str, note(FullPeakIndexed),"="))	
		if (WaveExists(recip))
			Variable i
			String recipStr = "reciprocal lattice vectors:\r     a*           b*           c*\r"
			for (i=0;i<3;i+=1)
				sprintf str,"%+10.6f   %+10.6f   %+10.6f\r",recip[i][0],recip[i][1],recip[i][2]
				recipStr += str
			endfor
			recipStr = TrimBoth(recipStr)
		endif
		TextBox/N=textRecip/F=0/A=RB/X=4/Y=4 recipStr
	endif
	if (strlen(par2))
		SetDrawLayer UserFront
		SetDrawEnv linethick= 0.1
		DrawLine 313,540,313,743
	endif
	Textbox/C/N=stamp0/F=0/A=RB/X=0.1/Y=0.1 "\\Z06\\{\"%s %s\",date(), time()}"
	Textbox/C/N=stamp1/F=0/A=LB/X=0.1/Y=0.1 "\\Z06\\{\"%s\",CornerStamp1_()}"+":"+WinName(0, 1)
End
//
Static Function/S ParagraphsFromIndexing(fpi,ip)
	Wave fpi									// usually FullPeakIndexed
	Variable ip								// pattern number (usually 0)

	Variable maxLines=17					// max number of lines in one paragraph
	String line, par1="",par2=""

	String deg = SelectString(cmpstr(IgorInfo(2),"Windows")==0 && IgorVersion()<7, DEGREESIGN, S_DEGREE)	// use symbol font for old Windows
	String topLine = "\\Z09 (hkl)\t keV\t           pixel\t err"+deg

	Variable i, N = DimSize(fpi,0)
	for (i=0;i<min(N,maxLines);i+=1)	// set par1
		if (numtype(fpi[i][0][ip]))
			continue
		endif
		sprintf line,"\r(%s)	%.3f	(% 6.1f, % 6.1f)	%.3f",hkl2str(fpi[i][3][ip],fpi[i][4][ip],fpi[i][5][ip]),fpi[i][7][ip],fpi[i][9][ip],fpi[i][10][ip],fpi[i][8][ip]
		par1 += line
	endfor
	if (N>i)
		N = min(N,i+maxLines)
		for (;i<N;i+=1)						// set par2
			if (numtype(fpi[i][0][ip]))
				continue
			endif
			sprintf line,"\r(%s)	%.3f	(% 6.1f, % 6.1f)	%.3f",hkl2str(fpi[i][3][ip],fpi[i][4][ip],fpi[i][5][ip]),fpi[i][7][ip],fpi[i][9][ip],fpi[i][10][ip],fpi[i][8][ip]
			par2 += line
		endfor
	endif
	par1 = SelectString(strlen(par1),"",topLine+par1)
	par2 = SelectString(strlen(par2),"",topLine+par2)
	return par1+";"+par2
End


Function/T pickIndexingFunction(path)
	String path						// path to folder on disk to be searched
	if (strlen(path)<1)
		path = getBinaryPath()								// default path to the executables
	endif

	String defaultExe
	Variable isMac = stringmatch(igorInfo(2),"Macintosh")
	if (isMac && stringmatch(igorinfo(4),"PowerPC"))
		defaultExe = "Euler_ppc"
		isMac = 1
	elseif (isMac && stringmatch(igorinfo(4),"Intel"))
		defaultExe = "Euler_i386"
		isMac = 2
	else
		defaultExe = "Euler.exe"						// default for MSWindows
	endif
	String str = GetIndexingFuncOrExec()			// = "indexingIgorFunc;indexingExecutable;"
	String indexFuncOld=StringFromList(0,str), exeOld=StringFromList(1,str)
	String old = exeOld + indexFuncOld				// one is always "", so this works
	old = SelectString(strlen(old), "default", old)
	STRUCT indexingTypeStruct indexPref
	indexPref.executable = exeOld
	indexPref.internal = indexFuncOld

	// get list of possibilities
	NewPath/O/Q/Z IndexingSearchPath, path
	if (V_Flag)
		DoAlert 0,"Unable to set path to "+path
		return ""
	endif
	String exeList=""
	if (isMac)
		exeList = IndexedFile(IndexingSearchPath,-1,"????")
		Variable i, bad
		String fname
		for (i=ItemsInList(exeList)-1;i>=0;i-=1)
			fname = StringFromList(i,exeList)
			GetFileFolderInfo/Q/Z path+fname
			bad = (strlen(S_fileType+S_creator) || V_isInvisible || !V_isFile || V_isAliasShortcut || V_isStationery || V_logEOF<300e3)
			bad = bad || StringMatch(fname,"*.ipf")		// exclude Igor procedure files
			bad = bad || StringMatch(fname,"*.exe")		// on a Mac, exclude .exe files
			if (isMac ==1)											// for ppc, exclude files ending in i386 or intel
				bad = bad || StringMatch(fname,"*i386") || StringMatch(fname,"*intel")
			elseif (isMac==2)										// for intel, exclude files ending with ppc
				bad = bad || StringMatch(fname,"*ppc")	
			endif
			if (bad)
				exeList = RemoveFromList(fname,exeList)
			endif
		endfor
	else
		exeList = IndexedFile(IndexingSearchPath,-1,".exe")	// on MSWindows, only *.exe files
	endif

	String indexFuncList = FunctionList("runIndexing*",";", "KIND:2,VALTYPE:8,NPARAMS:1")
	indexFuncList = RemoveFromList("runIndexingEulerCommand",indexFuncList)
	if (ItemsInList(exeList+indexFuncList)<1)			// nothing to select from
		DoAlert 0,"No executables or indexing functions found.\r  Doing nothing."
		return ""
	endif

	// choose from list of possibilities
	String exe=exeOld
	if (strlen(exe)==0 && strlen(indexFuncOld)>0)
		exe = StringFromList(0,indexFuncOld)
	endif
	String popUpList = "[default binary];"
	if (strlen(exeList) && strlen(indexFuncList))
		popUpList += exeList+" ;"+indexFuncList
	elseif (strlen(exeList))
		popUpList += exeList
	elseif (strlen(indexFuncList))
		popUpList += indexFuncList
	endif
	if (strlen(WinList("Indexing_Internal*",";","WIN:128"))==0 && strlen(indexFuncList)==0)
		popUpList += "[#include internal Igor functions];"	// Indexing_Internal.ipf not usually included, option to include it
	endif
	Prompt exe,"name of indexing program",popup, popupList
	DoPrompt/HELP="Choose a different program for indexing.\r  Euler is fast\r  IndexAll is slow\r  default is usually best." "indexing program",exe
	if (V_flag)
		return ""
	endif
	if (StringMatch(exe,"[#include internal Igor functions]"))	// now include Indexing_Internal.ipf
		Execute/P "INSERTINCLUDE  \"Indexing_Internal\", version>=0.30"
		Execute/P "COMPILEPROCEDURES "
		Execute/P "pickIndexingFunction(\"\")"
		DoAlert 0, "re-running to Select an internal Igor indexing routine."
		return ""															// re-running pickIndexingFunction(""), allows chance for COMPILEPROCEDURES to run
	endif
	exe = SelectString(StringMatch(exe,"[default binary]") || strlen(exe)<2, exe, "")	// change invalids to ""

	String exeNew="", indexFuncNew=""
	if (WhichListItem(exe,indexFuncList)>=0)
		indexFuncNew = exe
	else
		exeNew = exe
	endif

	String final = exeNew + indexFuncNew
	final = SelectString(strlen(final), "default", final)
	printf "\r\r ==========  Changing Indexing Executable from \"%s\"  -->  \"%s\"  ==========\r\r\r", old, final

	if (strlen(exeNew))				// set globals, or delete depending upon what was chosen
		String/G root:Packages:micro:indexingExecutable = exeNew
		indexPref.executable = exeNew
	else
		KillStrings/Z root:Packages:micro:indexingExecutable
		indexPref.executable = ""
	endif
	if (strlen(indexFuncNew))
		String/G root:Packages:micro:indexingIgorFunc = indexFuncNew
		indexPref.internal = indexFuncNew
	else
		KillStrings/Z root:Packages:micro:indexingIgorFunc
		indexPref.internal = ""
	endif
	SavePackagePreferences/FLSH=1 "microGeo","indexing",0,indexPref
	return exe
End
//
Static Function/T GetIndexingFuncOrExec()		// returns "indexingIgorFunc;indexingExecutable;"
	String indexingIgorFunc = StrVarOrDefault("root:Packages:micro:indexingIgorFunc","")
	String indexingExecutable = StrVarOrDefault("root:Packages:micro:indexingExecutable","")

	Variable funcExists = ItemsInList(FunctionList(indexingIgorFunc,";", "KIND:2,VALTYPE:8,NPARAMS:1"))>0
	Variable execExists = strlen(indexingExecutable)>0
	indexingIgorFunc = SelectString(funcExists, "", indexingIgorFunc)
	indexingExecutable = SelectString(execExists, "", indexingExecutable)

	if (strlen(indexingIgorFunc+indexingExecutable)<1)	// nothing defined locally, check preferences
		STRUCT indexingTypeStruct indexPref
		LoadPackagePreferences/MIS=0 "microGeo", "indexing", 0, indexPref
		indexingIgorFunc = indexPref.internal
		indexingExecutable = indexPref.executable
	endif

	funcExists = ItemsInList(FunctionList(indexingIgorFunc,";", "KIND:2,VALTYPE:8,NPARAMS:1"))>0
	indexingIgorFunc = SelectString(funcExists, "", indexingIgorFunc)
	String EulerPath = getBinaryPath()						// path to the Euler executable
	GetFileFolderInfo/Q/Z EulerPath+indexingExecutable
	indexingExecutable = SelectString(V_Flag==0 && V_isFile, "", indexingExecutable)

	return indexingIgorFunc+";"+indexingExecutable+";"
End
//
Static Structure indexingTypeStruct	// structure definition of a measured or calculated zone
	char executable[400]					// name of executable (empty implies the default)
	char internal[400]						// name of Igor function to use for indexing (empty implies an executable)
EndStructure


// using the result from FitPeaks() run Euler, and read in the results from the index file
Function/WAVE runIndexingEulerCommand(args)
	String args
	Wave FullPeakList = $StringByKey("FullPeakList",args)
	Variable keVmaxCalc = NumberByKey("keVmaxCalc",args)				// 17, maximum energy to calculate (keV)
	Variable keVmaxTest = NumberByKey("keVmaxTest",args)				// 26, maximum energy to test (keV)  [-t]
	Variable angleTolerance = NumberByKey("angleTolerance",args)	// 0.25, angular tolerance (deg)
	Variable hp = NumberByKey("hp",args)										// preferred hkl
	Variable kp = NumberByKey("kp",args)
	Variable lp = NumberByKey("lp",args)
	Variable cone = NumberByKey("cone",args)								// angle from preferred hkl, (0 < cone ≤ 180°)
	Variable maxSpots = NumberByKey("maxSpots",args)						// -n max num. of spots from data file to use, default is 250
	Wave FullPeakList1 = $StringByKey("FullPeakList1",args)
	Wave FullPeakList2 = $StringByKey("FullPeakList2",args)
	Variable depth = NumberByKey("depth",args)								// depth (micron)
	Variable printIt = NumberByKey("printIt",args)
	maxSpots = ((maxSpots>2) && numtype(maxSpots)==0) ? maxSpots : -1
	printIt = numtype(printIt) ? (strlen(GetRTStackInfo(2))==0) : printIt

	cone = cone==180 ? 180-0.001 : cone
	Variable badNums= !(keVmaxCalc>1 && keVmaxCalc<INDEXING_MAX_CALC) || !(keVmaxTest>1 && keVmaxTest<INDEXING_MAX_TEST)
	badNums += !(angleTolerance>=0.01 && angleTolerance<10)
	badNums += numtype(hp+kp+lp)
	badNums += !(cone>1 && cone<=180)
	if (badNums)
		DoAlert 0, "Invalid inputs sent to runIndexingEulerCommand()"
		printf "runIndexingEulerCommand(%s,%g,%g,%g)\r",NameOfWave(FullPeakList),keVmaxCalc,keVmaxTest,angleTolerance
		return $""
	elseif (!WaveExists(FullPeakList))
		DoAlert 0, "the input wave does not exist in runIndexingEulerCommand()"
		return $""
	elseif (DimSize(FullPeakList,0)<1 || DimSize(FullPeakList,1)<11)
		DoAlert 0, "Full peak list '"+NameOfWave(FullPeakList)+"' is empty or the wrong size"
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

	Variable isMac = stringmatch(igorInfo(2),"Macintosh")
	Variable isWin = stringmatch(igorInfo(2),"Windows")

	// first write the command file to drive Euler
	String upath=SpecialDirPath("Temporary",0,1,0)// local path (probably unix path) for the command line
	String mpath=SpecialDirPath("Temporary",0,0,0)// mac style path for use only within Igor
	PathInfo home
	if (V_flag && !isWin)						// the path "home" exists, use it (but not on windows, to avoid certain permission issue since 2018)
		mpath = S_path
		if (isMac)
			upath = ParseFilePath(5,S_path,"/",0,0)	// Convert HFS>POSIX only on Mac, not Windows
		endif
	endif
	NewPath/O/Q/Z EulerCalcFolder, mpath				// need a new path name since "home" may not exist
	if (V_flag)
		DoAlert 0, "Unable to create path to do Euler calculation, try saving the experiment first."
		return $""
	endif
	String peakFile="generic_Peaks.txt"				// name of file with input peak positions
	if(FullPeakList2Qfile(FullPeakList,peakFile,"EulerCalcFolder",FullPeakList1=FullPeakList1,FullPeakList2=FullPeakList2,depth=depth))// convert peaks to a Qlist+intens, and write to a file
		return $""												// nothing happened
	endif

	// find the full path name of the Euler executable
	String exe, EulerPath=getBinaryPath()					// path to the Euler executable
	if (isMac && stringmatch(igorinfo(4),"PowerPC"))	// find the default name for this architecture
		exe = "Euler_ppc"
	elseif (isMac && stringmatch(igorinfo(4),"Intel"))
		exe = "Euler_i386"
	elseif (isWin)
		exe = "Euler.exe"
	else
		exe = "Euler"
	endif
	if (strlen(StringFromList(1,GetIndexingFuncOrExec())))
		exe = StringFromList(1,GetIndexingFuncOrExec())
	endif
	GetFileFolderInfo/Q/Z EulerPath+exe
	EulerPath += SelectString(V_Flag==0 && V_isFile,"Euler",exe)	// use just plane Euler if Euler_ppc or Euler_i386 does not exist
	if (stringmatch(igorInfo(2),"Macintosh"))			// on Mac, convert EulerPath from HFS to Posix
		EulerPath = ParseFilePath(5,EulerPath,"/",0,0)
	endif
	if (strlen(EulerPath)<1)
		DoAlert 0, "cannot find the executable '"+exe+"'"
		return $""
	endif

	String cmd, result
	Variable err=0
	if (isMac)
		if (maxSpots > 2)
			sprintf cmd "do shell script \"cd \\\"%s\\\" ; \\\"%s\\\" -k %g -t %g -a %g -h %d %d %d -c %g -n %d -f %s -q\"",upath,EulerPath,keVmaxCalc,keVmaxTest,angleTolerance,hp,kp,lp,cone,maxSpots,peakFile
		else
			sprintf cmd "do shell script \"cd \\\"%s\\\" ; \\\"%s\\\" -k %g -t %g -a %g -h %d %d %d -c %g -f %s -q\"",upath,EulerPath,keVmaxCalc,keVmaxTest,angleTolerance,hp,kp,lp,cone,peakFile
		endif
		ExecuteScriptText/Z cmd
		result = S_value
	elseif (isWin)		// create a .BAT file to run Euler
		variable BatchFileRN
		variable OutputFileRN

		String BatchFileName = "EulerExecute.bat"
		String FullBatchFileName = mpath + BatchFileName
		String WindowsEulerPath = MacToWindows(EulerPath)
		String OutputFile = MacToWindows(mpath + "EulerOutput.txt")
		String FullPeakFileName = MacToWindows(mpath + peakFile)
		
		Open /Z BatchFileRN  as FullBatchFileName
		variable errflag=V_flag
		if (errflag != 0)
			printf "Failed to open %s for writing", FullBatchFileName
			return $""
		endif
		if (maxSpots > 2)
			fprintf BatchFileRN, "\"%s\" -k %g -t %g -a %g -h %d %d %d -c %g -n %d -f \"%s\" > \"%s\"\r\n", WindowsEulerPath,keVmaxCalc,keVmaxTest,angleTolerance,hp,kp,lp,cone,maxSpots,FullPeakFileName,OutputFile
		else
			fprintf BatchFileRN, "\"%s\" -k %g -t %g -a %g -h %d %d %d -c %g -f \"%s\" > \"%s\"\r\n", WindowsEulerPath,keVmaxCalc,keVmaxTest,angleTolerance,hp,kp,lp,cone,FullPeakFileName,OutputFile
		endif
		close BatchFileRN
		sprintf cmd "\"%s\"", MacToWindows(FullBatchFileName)
		ExecuteScriptText cmd
		open /R/Z OutputFileRN as OutputFile
		if (V_flag != 0)
			print "Error opening batch output file ", OutputFile
			result = ""
		else
			FReadLine /T="" OutputFileRN, result
			Close OutputFileRN
		endif
	endif

	if (ItemsInList(GetRTStackInfo(0))<2)		// print everything if run from command line
		print result
	endif

	err = strsearch(result,"writing output to file '",0)<0
	if (err)
		if (!NumVarOrDefault("root:Packages:micro:Index:SKIP_EULER_ERRORS",0))
			if (printIt)
				DoAlert 0, "failure in runIndexingEulerCommand()"
			endif
			print "\r\r"
			print cmd
			print "\r\r"
			print result
		endif
		return $""								// there is no index file
	endif
	Variable i = strsearch(result,"writing output to file '",0)
	String line = result[i,Inf]
	i = strsearch(line,"'",0)+1
	line = line[i,Inf]
	String indexFile = line[0,strsearch(line,"'",0)-1]
	if (!strlen(indexFile))
		DoAlert 0, "in runIndexingEulerCommand(), cannot find the index file"
		print "\r\r"
		print result
		return $""								// there is no index file
	endif
	return readIndexFile(indexFile,"EulerCalcFolder")	// read an index file and create an array to hold the results, return array name
End
//
// 
Static Function/WAVE readIndexFile(indexFile,path)
	String indexFile						// name of index file, the output from Euler
	String path							// name of Igor path to go with indexFile

	Variable f = OpenFileOf_ftype(indexFile,"IndexFile",path)// open an index file
	if (f<1)
		return $""						// could not get any info from indexFile 
	endif
	FStatus f
	String buffer=""
	buffer = PadString(buffer,V_logEOF-V_filePos,0x20)
	FBinRead f, buffer					// read entire file into buffer
	Close f
	buffer = ReplaceString("\r\n",buffer,"\n")				// ensure that line term is a new-line only
	buffer = ReplaceString("\n\r",buffer,"\n")
	buffer = ReplaceString("\r",buffer,"\n")
	if (strsearch(buffer,"$IndexFile\n",0)==0)
		buffer = ReplaceString("$IndexFile\n",buffer,"$fileType	IndexFile\n",0,1)
	endif
	Variable Nbuf = strlen(buffer)

	STRUCT microGeometry geo
	if (FillGeometryStructDefault(geo, alert=1))		//fill the geometry structure with test values
		return $""													// can not compute (px,py), but still a valid result
	endif

	Variable ipattern=0										// number of current patterns
	Variable i1,i0 = strsearch(buffer,"$pattern"+num2istr(ipattern),0)		// put into "key=value" everything up pattern ipattern
	i0= i0<0 ? Inf : i0-1
	String wnote = keyStrFromBuffer(buffer[0,i0])
	String peakListNameList = StringByKey("peakListWave",wnote,"=")
	String peakListName = ParseFilePath(3,StringFromList(0,peakListNameList,","),":",0,0)
	String FullPeakIndexedName = CleanupName("FullPeakIndexed"+ReplaceString("FullPeakList",peakListName,""),0)
	Variable Npatterns = NumberByKey("NpatternsFound",wnote,"=")
	if (!(Npatterns>0))
			return $""
	endif
	Variable startx,groupx, starty,groupy, ddLocal
	startx = NumberByKey("startx",wnote,"=")
	groupx = NumberByKey("groupx",wnote,"=")
	starty = NumberByKey("starty",wnote,"=")
	groupy = NumberByKey("groupy",wnote,"=")
	startx = numtype(startx) ? FIRST_PIXEL : startx
	groupx = numtype(groupx) ? 1 : groupx
	starty = numtype(starty) ? FIRST_PIXEL : starty
	groupy = numtype(groupy) ? 1 : groupy

	Make/N=3/O/D/FREE qBL, qvec
	Variable dNum											// can have multiple detectors at once
	String list1
	Variable/C pz
	Variable px,py
	Variable qx,qy,qz,hh,kk,ll,int,keV,err
	Variable Ni=-1, cols,i,val, nLines
	do
		i0 +=1												// range of header info for this pattern
		i1 = strsearch(buffer,"\n$array"+num2istr(ipattern),i0)
		list1 = keyStrFromBuffer(buffer[i0,i1])
		list1 = RemoveByKey("pattern"+num2istr(ipattern),list1,"=")
		wnote = MergeKeywordLists(wnote,list1,0,"=","")

		// buffer is now the data (and the leading $arrayN line)
		i0 = i1+1
		i1 = strsearch(buffer,"$pattern"+num2istr(ipattern+1),i0)
		i1 = (i1<1 ? Nbuf : i1)-1
		sscanf buffer[i0,i1],"$array%d %d %d",i,cols,nLines
		if (ipattern==0)									// for first pattern, make the output wave
			Ni = nLines
			Make/N=(Ni,cols+2,Npatterns)/O $FullPeakIndexedName// make is long enough to also hold the pixel positions (an extra 3 cols)
			Wave FullPeakIndexed = $FullPeakIndexedName
			FullPeakIndexed = NaN
		elseif (nLines>Ni)
			Ni = nLines
			Redimension/N=(Ni,cols+2,Npatterns) FullPeakIndexed
		endif

		i0 = strsearch(buffer,"\n",i0+1)+1				// start of first data
		for (i=0;i<nLines;i+=1)
			sscanf buffer[i0,i1], "    [ %d]   (%g  %g %g)     ( %d %d  %d)   %g,    %g,    %g",val,qx,qy,qz,hh,kk,ll,int,keV,err
			qBL[0]=qx ; qBL[1]=qy ; qBL[2]=qz
			FullPeakIndexed[i][0][ipattern] = qBL[0]
			FullPeakIndexed[i][1][ipattern] = qBL[1]
			FullPeakIndexed[i][2][ipattern] = qBL[2]
			FullPeakIndexed[i][3][ipattern] = hh
			FullPeakIndexed[i][4][ipattern] = kk
			FullPeakIndexed[i][5][ipattern] = ll
			FullPeakIndexed[i][6][ipattern] = int
			FullPeakIndexed[i][7][ipattern] = keV
			FullPeakIndexed[i][8][ipattern] = err
			i0 = strsearch(buffer,"\n",i0+1)+1			// start of next line
		endfor
		for (i=0;i<nLines;i+=1)							// convert the indexed qhats to pixel positions
			qvec = FullPeakIndexed[i][p][ipattern]
//			pz = q2pixel(geo.d[dNum],qvec)
//			FullPeakIndexed[i][9][ipattern] = (real(pz)-(startx-FIRST_PIXEL)-(groupx-1)/2)/groupx		// change to binned pixels
//			FullPeakIndexed[i][10][ipattern] = (imag(pz)-(starty-FIRST_PIXEL)-(groupy-1)/2)/groupy		// pixels are still zero based
//			pz = q2pixel(geo.d[dNum],qvec)
//			px = real(pz)															// trial pixels
//			py = imag(pz)
			for (dNum=0;dNum<geo.Ndetectors;dNum+=1)
				pz = q2pixel(geo.d[dNum],qvec)
				px = real(pz)														// trial pixels for this detector
				py = imag(pz)
				if (px>=0 && px <geo.d[dNum].Nx && py>=0 && py<geo.d[dNum].Ny)
					FullPeakIndexed[i][9][ipattern] = (px-(startx-FIRST_PIXEL)-(groupx-1)/2)/groupx		// change to binned pixels
					FullPeakIndexed[i][10][ipattern] = (py-(starty-FIRST_PIXEL)-(groupy-1)/2)/groupy	// pixels are still zero based
					FullPeakIndexed[i][11][ipattern] = dNum
					break
				endif
			endfor
		endfor
		ipattern += 1
	while (ipattern<Npatterns && i0<(Nbuf-1))
	for (ipattern=0;ipattern<Npatterns;ipattern+=1)		// remove all (000) for all patterns
		for (i=0;i<Ni;i+=1)
			if (FullPeakIndexed[i][0][ipattern]==0 && FullPeakIndexed[i][1][ipattern]==0 && FullPeakIndexed[i][2][ipattern]==0)
				FullPeakIndexed[i][][ipattern] = NaN
			endif
		endfor
	endfor
	wnote = ReplaceStringByKey("waveClass",wnote,"IndexedPeakList","=")

	// a legacy section
	Variable SpaceGroup = NumberByKey("SpaceGroup",wnote,"=")
	if (numtype(SpaceGroup))					// could not find 'SpaceGroup', see if it is still called 'latticeStructure'
		SpaceGroup = NumberByKey("latticeStructure",wnote,"=")
		if (numtype(SpaceGroup)==0)
			wnote = ReplaceNumberByKey("SpaceGroup",wnote,spaceGroup,"=")
		endif
	endif
	String SpaceGroupID = StringByKey("SpaceGroupID",wnote,"=")
	Variable SpaceGroupIDnum = NumberByKey("SpaceGroupIDnum",wnote,"=")
	if (strlen(SpaceGroupID)<1)
		SpaceGroupID = LatticeSym#FindDefaultIDForSG(SpaceGroup)
	endif
	SpaceGroupIDnum = numtype(SpaceGroupIDnum) ? str2num(SpaceGroupID) : SpaceGroupIDnum

	Note/K FullPeakIndexed, wnote
	SetDimLabel 1,0,Qx,FullPeakIndexed				;	SetDimLabel 1,1,Qy,FullPeakIndexed
	SetDimLabel 1,2,Qz,FullPeakIndexed				;	SetDimLabel 1,3,h,FullPeakIndexed
	SetDimLabel 1,4,k,FullPeakIndexed				;	SetDimLabel 1,5,l,FullPeakIndexed
	SetDimLabel 1,6,Intensity,FullPeakIndexed	;	SetDimLabel 1,7,keV,FullPeakIndexed
	SetDimLabel 1,8,angleErr,FullPeakIndexed	;	SetDimLabel 1,9,pixelX,FullPeakIndexed
	SetDimLabel 1,10,pixelY,FullPeakIndexed		;	SetDimLabel 1,11,detNum,FullPeakIndexed

	KillWaves/Z M_Lower,M_Upper,W_LUPermutation,M_x
	return FullPeakIndexed
End
//
Static Function OpenFileOf_ftype(fname,ftype,path)// open a file only if it is of type ftype
	String fname									// full path name to file with tagged geometry values
	String ftype										// the required file identifier, included as a tag (ftype is optional)
	String path										// name of Igor path

	Variable OK = 0
	Variable f = 0


	if (strlen(ftype)<1)
		OK = 1
	else
		String list = getListOfTypesInFile(fname,path)
		OK = (WhichListItem(ftype,list)>=0)
	endif

	if (OK)
		Open/M="file with a keyword list"/P=$path/R/Z=2 f as fname
		if (strlen(S_fileName)<1 || !f)
			if (f)
				Close f
			endif
			return 0
		endif
	endif
	return f
End
//
Static Function/S keyStrFromBuffer(buffer)// get key list from the buffer
	String buffer 							// long sting with line feeds, each line is of the form:   "$tag   value   // comment\n"

	Variable N=strlen(buffer)
	Variable i0,i1							// indicies into buffer
	String tagName, value					// the found tag and value (if they exist)
	String line								// one line of buffer
	String list=""							// the result
	buffer = ReplaceString("\r\n",buffer,"\n")	// ensure that line term is a new-line only
	buffer = ReplaceString("\n\r",buffer,"\n")
	buffer = ReplaceString("\r",buffer,"\n")
	i0 = 0
	i1 = strsearch(buffer,"\n",i0)
	do
		line = buffer[i0,i1]						// one line from the buffer
		tagValueFromLine(line,tagName,value)	// get tag and value from this line
		if (strlen(tagName) && !keyInList(tagName,list,"=","")) // valid tag, and not yet in list
			list = ReplaceStringByKey(tagName,list,value,"=")
		endif
		i0 = i1+1
		i1 = strsearch(buffer,"\n",i0)
		i1 = (i1<0) ? N-1 : i1
	while(i0<N)
	return list
End
//
Static Function/S tagValueFromLine(lineIn,tagName,value)// find the tag and value from a line
	String lineIn								// line of text to process
	String &tagName							// the found tag (if it exists)
	String &value								// the found value for that tag (if it exists)

	String line=lineIn
	tagName = ""
	value = ""
	Variable i,dollar = char2num("$")
	if (char2num(line)!=dollar)				// if it starts with a $, probable tag so process it
		return ""
	endif

	i = strsearch(line,"//",0)					// strip off comments
	if (i>=0)
		line = line[0,i-1]
	endif

	for (i=0;char2num(line[i+1])>32;i+=1)	// find end of tag, it ends with a space or lower
	endfor
	tagName = line[1,i]
	if (strlen(tagName)<1)						// no tag
		return ""
	endif

	for (i=i+1;char2num(line[i])<=32;i+=1)	// find first non-white space, start of value
	endfor
	value = line[i,Inf]							// value associated with tagName

	for (i=strlen(value)-1;char2num(value[i])<=32 && i>0;i-=1)	// strip off trailing whitespace
	endfor
	value = ReplaceString(";",value[0,i],":")	// cannot have a semicolon or equals in the value
	return ReplaceStringByKey(tagName,"",value,"=")
End
//Function test()
//	Variable f
//	Open/M="test file of tagged lines"/P=home/R/T="TEXT" f as "Macintosh HD:Users:tischler:data:Sector 34:August 2006:calibration:generic_Index.txt"
//	print S_fileName
//	FStatus f
//	print "V_logEOF=",V_logEOF
//	String buffer=""
//	buffer = PadString(buffer,V_logEOF,0x20)
//	FBinRead f, buffer
//	Close f
//	print "strlen(buffer) = ",strlen(buffer)
//	String list = keyStrFromBuffer(buffer)
//	print list
//End
//
//Function test()
//	String tagName,value,line
//	line = "$EulerAngles {-171.29728810, 149.32171846, 229.41456865}	// Euler angles for this pattern (deg)\n"
//	String output = tagValueFromLine(line,tagName,value)
//	printf "tagName = '%s'\r",tagName
//	printf "value = '%s'\r",value
//	printf "output = '%s'\r",output
//End


// remove bkg from an image, identify the peaks, and fit them all
Function/S RemoveBackgroundImage(rawImage,fractionBkg,[printIt])
	Wave rawImage				// the raw image from the image file
	Variable fractionBkg		// fraction of image that is background, use something like 0.99 or 0.995 (even 0.7 works OK)
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt

	if (!WaveExists(rawImage) || WaveDims(rawImage)!=2 || !(fractionBkg==limit(fractionBkg,1e-3,1)))
		String imageName = SelectString(WaveExists(rawImage),"",NameOfWave(rawImage))
		fractionBkg = fractionBkg==limit(fractionBkg,1e-3,1) ? fractionBkg : 0.8
		Prompt imageName,"name of image with peaks to find",popup,reverseList(WaveListClass("speImage;rawImage*","*","DIMS:2"))
		Prompt fractionBkg,"fraction of image that is background, usually in range [0.7, 0.995]"
		DoPrompt/Help="3D-Xray Diffraction[Fit Peaks]" "peak fitting",imageName,fractionBkg
		if (V_flag)
			return ""
		endif
		Wave rawImage = $imageName
		printf "FitPeaks(%s,%g)\r",imageName,fractionBkg
		printIt = 1
	endif
	if (!WaveExists(rawImage) || WaveDims(rawImage)!=2 || !(fractionBkg==limit(fractionBkg,1e-3,1)))
		return ""
	endif
	Variable timer=startMSTimer
//	Variable timer1 = startMSTimer
	Variable thresh = numpnts(rawImage)*fractionBkg

	WaveStats/Q/M=1 rawImage
	Make/N=200/D/O FitPeaks_Hist
	Variable i
	for (i=-1;i<10 || numtype(i);)							// keep expanding the histogram range while i<10,   I want i to be bigger than 10
		SetScale/I x,V_min,V_max,"",FitPeaks_Hist			// lower high end of histogram to zoom in around zero
		Histogram/B=2 rawImage,FitPeaks_Hist
		Integrate/P FitPeaks_Hist/D=FitPeaks_Hist_INT	// integrate the histogram
		i=BinarySearchInterp(FitPeaks_Hist_INT,thresh)	// find point where FitPeaks_Hist_INT[i] == thresh
		V_max = numtype(i) ? V_max/2 : pnt2x(FitPeaks_Hist_INT,i+4)
	endfor
//	print "first part took ",stopMSTimer(timer1)/1e6,"sec"
//	timer1 = startMSTimer

	thresh = DimDelta(FitPeaks_Hist,0)*i + DimOffset(FitPeaks_Hist,0)
	ImageThreshold/Q/T=(thresh)/M=0 rawImage					// create M_ImageThresh
	Wave M_ImageThresh=M_ImageThresh

	ImageStats/M=1 M_ImageThresh
	Variable fracBlack = V_avg*V_npnts/255/V_npnts
	//	printf "number of points in peak is %d out of %d points (%.2f%%)\r",fracBlack*V_npnts,V_npnts,fracBlack*100
	Make/N=(9,9)/O/B/U RemoveBackgroundImage_black
	RemoveBackgroundImage_black = 255
	for(;fracBlack<0.10;)											// loop until fracBlack is atleast 10%
		Imagemorphology/S=RemoveBackgroundImage_black/O BinaryDilation M_ImageThresh
		ImageThreshold/Q/T=(2)/M=0 M_ImageThresh
		ImageStats/M=1 M_ImageThresh
		fracBlack = V_avg*V_npnts/255/V_npnts
		//	printf "number of points in peak is %d out of %d points (%.2f%%)\r",fracBlack*V_npnts,V_npnts,fracBlack*100
	endfor
	ImageTransform/O invert M_ImageThresh						//	M_ImageThresh = !M_ImageThresh

//	print "second part took ",stopMSTimer(timer1)/1e6,"sec"
//	timer1 = startMSTimer

	ImageInterpolate /f={.125,.125} bilinear rawImage
	Redimension/S M_InterpolatedImage

	KillWaves/Z FitPeaks_rawImage_smaller
	Rename M_InterpolatedImage FitPeaks_rawImage_smaller		// make a source image using fewer pixels
//	print "\t\t1st interpolate part took ",stopMSTimer(timer1)/1e6,"sec"

//	timer1 = startMSTimer
 	ImageInterpolate /f={.125,.125} bilinear M_ImageThresh		// make a mask with fewer pixels
	KillWaves/Z FitPeaks_mask_smaller
	Rename M_InterpolatedImage FitPeaks_mask_smaller
	Redimension/B/U FitPeaks_mask_smaller
	ImageRemoveBackground/F/R=FitPeaks_mask_smaller/P=2/W FitPeaks_rawImage_smaller	// creates M_RemovedBackground, uses 3rd order polynomial
	Wave W_BackgroundCoeff=W_BackgroundCoeff
	//	print W_BackgroundCoeff
//	print "\t\t2nd interpolate part took ",stopMSTimer(timer1)/1e6,"sec"

//	timer1 = startMSTimer
	ImageInterpolate /f={8,8} bilinear M_RemovedBackground		// creates M_InterpolatedImage, a bkg with full number of pixels
	Wave M_InterpolatedImage=M_InterpolatedImage

	String outName = NameOfWave(rawImage)+"NoBkg"
	String fldr = GetWavesDataFolder(rawImage,1)
	if (CheckName(fldr+outName,1))
		outName = CleanupName(outName,0)
	endif
	outName = fldr+outName
	Duplicate/O rawImage $outName
	Wave noBkg = $outName

	String wnote=ReplaceStringByKey("rawIgorImage",note(rawImage),GetWavesDataFolder(rawImage,2),"=")
	wnote = ReplaceNumberByKey("fractionBkg",wnote,fractionBkg,"=")
	wnote = ReplaceStringByKey("waveClass",wnote,"rawImageNoBkg","=")
	Note/K noBkg, wnote
	Redimension/I noBkg											// signed int32 version of rawImage
	noBkg = rawImage - M_InterpolatedImage[p][q]					// p,q because M_InterpolatedImage is smaller

	Variable sec0=stopMSTimer(timer)/1e6
	if (printIt)
		printf "removing background  took %g sec\r",sec0
		printf "  result stored in the wave '%s'\r", GetWavesDataFolder(noBkg,2)
	endif

	KillWaves/Z FitPeaks_rawImage_smaller,FitPeaks_mask_smaller, M_ImageThresh, M_RemovedBackground, M_InterpolatedImage, W_BackgroundCoeff
	KillWaves/Z FitPeaks_Hist_INT,FitPeaks_Hist, RemoveBackgroundImage_black
	return GetWavesDataFolder(noBkg,2)
End







Function/T MakeMaskThreshold(image,threshold,[dilate,printIt])
	Wave image						// the image from the image file
	Variable threshold				// mask off all pixels with values above threshold
	Variable dilate
	Variable printIt
	dilate = ParamIsDefault(dilate) ? 0 : dilate
	dilate = dilate>=0 ? dilate : 0
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt

	String imageName
	if (!WaveExists(image) || WaveDims(image)!=2 || !(threshold>0))
		imageName = SelectString(WaveExists(image),"",NameOfWave(image))
		threshold = threshold>0 ? threshold : NaN
		Prompt imageName,"name of image with peaks to find",popup,reverseList(WaveListClass("speImageNoBkg;speImage;rawImageNoBkg;rawImage*","*","DIMS:2,BYTE:0"))
		Prompt threshold,"max pixel, masks off all pixels above this, NaN masks peaks in peakList"
		Prompt dilate,"expand around peaks (pixels)"
		DoPrompt "make mask",imageName,threshold,dilate
		if (V_flag)
			return ""
		endif
		Wave image = $imageName
		printf "•MakeMaskThreshold(%s,%g",imageName,threshold
		if (dilate>0)
			printf ",dilate=%g",dilate
		endif
		printf ")\r"
		printIt = 1
	endif
	if (!WaveExists(image) || WaveDims(image)!=2)
		return ""
	endif
	if (!(threshold>0))
		Wave peakList = FindPeakListForImage(image)
		if (!WaveExists(peakList))
			DoAlert 0,"MakeMaskThreshold(), could not find peakList & threshold is invalid."
			return ""
		endif
		threshold = -1							// this is a flag for making mask
	endif
	dilate = dilate>=0 ? dilate : 0

	String maskName=NameOfWave(image)+"Mask"
	Duplicate/O image,$maskName
	Wave mask=$maskName
	if (threshold>0)							// make mask using a threshold
		ImageThreshold/O/Q/T=(threshold)/M=0 mask
	else											// make mask using FullPeaksList
		// set peaks to 255
		Redimension/B/U mask
		mask = 0
		Variable px,py,i,N=DimSize(peakList,0)
		for (i=0;i<N;i+=1)
			px = round(peakList[i][0])
			py = round(peakList[i][1])
			mask[px][py] = 255
		endfor
	endif

	//	ImageThreshold/O/Q/T=(threshold)/I/M=0 mask
	String wNote=note(image)
	wNote = ReplaceStringByKey("waveClass",wNote,"imageMask","=")
	Note/K mask,wNote
	if (dilate)
		ImageMorphology/E=4/I=(dilate)/O BinaryDilation  mask
		//		ImageMorphology/E=4/I=(dilate)/O BinaryErosion  mask
	endif
	return GetWavesDataFolder(mask,2)
End
//Function/T MakeMaskThreshold(image,threshold,[dilate])
//	Wave image						// the image from the image file
//	Variable threshold				// mask off all pixels with values above threshold
//	Variable dilate
//	dilate = ParamIsDefault(dilate) ? 0 : dilate
//	dilate = dilate>=0 ? dilate : 0
//
//	String imageName
//	Variable printIt = ItemsInList(GetRTStackInfo(0))<2 || stringmatch(StringFromList(0,GetRTStackInfo(0)),"IndexButtonProc")
//	if (!WaveExists(image) || WaveDims(image)!=2 || !(threshold>0))
//		imageName = SelectString(WaveExists(image),"",NameOfWave(image))
//		threshold = threshold>0 ? threshold : NaN
//		Prompt imageName,"name of image with peaks to find",popup,reverseList(WaveListClass("speImageNoBkg;speImage;rawImageNoBkg;rawImage*","*","DIMS:2,BYTE:0"))
//		Prompt threshold,"max pixel value, mask off all pixels above this"
//		Prompt dilate,"expand around peaks (pixels)"
//		DoPrompt "make mask",imageName,threshold,dilate
//		if (V_flag)
//			return ""
//		endif
//		Wave image = $imageName
//		printf "•MakeMaskThreshold(%s,%g",imageName,threshold
//		if (dilate>0)
//			printf ",dilate=%g",dilate
//		endif
//		printf ")\r"
//		printIt = 1
//	endif
//	if (!WaveExists(image) || WaveDims(image)!=2 || !(threshold>0))
//		return ""
//	endif
//	dilate = dilate>=0 ? dilate : 0
//
//	String maskName=NameOfWave(image)+"Mask"
//	Duplicate/O image,$maskName
//	Wave mask=$maskName
//	ImageThreshold/O/Q/T=(threshold)/M=0 mask
//	//	ImageThreshold/O/Q/T=(threshold)/I/M=0 mask
//	String wNote=note(image)
//	wNote = ReplaceStringByKey("waveClass",wNote,"imageMask","=")
//	Note/K mask,wNote
//	if (dilate)
//		ImageMorphology/E=4/I=(dilate)/O BinaryDilation  mask
//		//		ImageMorphology/E=4/I=(dilate)/O BinaryErosion  mask
//	endif
//	return GetWavesDataFolder(mask,2)
//End



//	-b box size (half width)
//	-R maximum R factor
//	-m min size of peak (pixels)
//	-M max number of peaks to examine(default=50)
//	-s minimum separation between two peaks (default=2*boxsize)
//	-t user supplied threshold (optional)
//	-p use -p L for Lorentzian (default), -p G for Gaussian
//	-S smooth the image if present
//	-T threshold ratio, set threshold to (ratio*[std dev] + avg) (optional)
//	-K mask_file_name (use pixels with mask==0)
Function/S FitPeaksWithExternal(image,minPeakWidth,boxSize,maxRfactor,threshAboveAvg,[mask,whoami,peakShape,smoothing,thresholdRatio,maxNum,printIt])
	Wave image						// the image from the file
	Variable minPeakWidth		// minimum width for an acceptable peak (pixels)
	Variable boxSize				// maximum width for an acceptable peak (pixels)
	Variable maxRfactor			// max R-factor
	Variable threshAboveAvg	// threshold above average value,  use 25
	Wave mask						// starting mask for the fit (this mask is unchanged by this routine)
	Variable whoami				// if passed, then just return the name of this routine
	String peakShape				// Lorentzian (default), or Gaussian
	Variable smoothing			// if TRUE, do a smoothing operation before peak search
	Variable thresholdRatio	// set threshold to (ratio*[std dev] + avg) if threshold passed as NaN (optional)
	Variable maxNum				// maximum number of peaks to fit (defaults to unlimited)
	Variable printIt
	if (!ParamIsDefault(whoami))
		return GetRTStackInfo(1)
	endif
	maxNum = ParamIsDefault(maxNum) ? -1 : maxNum
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt

	if (exists("root:Packages:micro:PeakFit:minPeakWidthLast")!=2)
		initPeaksPart()
	endif
	NVAR minPeakWidthLast = root:Packages:micro:PeakFit:minPeakWidthLast
	NVAR boxSizeLast = root:Packages:micro:PeakFit:boxSizeLast
	NVAR maxRfactorLast = root:Packages:micro:PeakFit:maxRfactorLast
	NVAR threshAboveAvgLast = root:Packages:micro:PeakFit:threshAboveAvgLast
	NVAR smoothingLast = root:Packages:micro:PeakFit:smoothingLast
	NVAR thresholdRatioLast = root:Packages:micro:PeakFit:thresholdRatioLast
	NVAR maxNumLast = root:Packages:micro:PeakFit:maxNumLast
	SVAR peakShapeLast = root:Packages:micro:PeakFit:peakShapeLast
	SVAR maskNameLast = root:Packages:micro:PeakFit:maskNameLast
	if (ParamIsDefault(peakShape))
		peakShape = peakShapeLast
	endif
	if (ParamIsDefault(mask))
		Wave mask = $""
	endif
	smoothing = ParamIsDefault(smoothing) ? smoothingLast : smoothing
	smoothing = !(!smoothing)
	thresholdRatio = ParamIsDefault(thresholdRatio) ? thresholdRatioLast : thresholdRatio
	thresholdRatio = thresholdRatio>0 ? thresholdRatio : NaN

	String imageName, line=""
	if (!WaveExists(image) || WaveDims(image)!=2 || !(minPeakWidth>0) || !(boxSize>0) || !(threshAboveAvg>0) || !(maxRfactor>0))
		imageName = SelectString(WaveExists(image),"",NameOfWave(image))
		minPeakWidth = minPeakWidth>.05 ? minPeakWidth : minPeakWidthLast
		boxSize = boxSize>.5 ? boxSize : boxSizeLast
		threshAboveAvg = threshAboveAvg>0 ? threshAboveAvg : threshAboveAvgLast
		maxRfactor = maxRfactor>0 ? maxRfactor : maxRfactorLast
		maxNum = maxNum>0 ? maxNum : maxNumLast
		Prompt imageName,"name of image with peaks to find",popup,reverseList(WaveListClass("speImageNoBkg;speImage;rawImageNoBkg;rawImage","*","DIMS:2"))
		Prompt minPeakWidth,"minimum allowed width of a peak"
		Prompt boxSize,"boxSize for finding peak"
		Prompt maxRfactor,"max R-factor for a peak fit"
		Prompt threshAboveAvg,"min peak height, use NaN for auto detect"
		Prompt maxRfactor,"max R-factor"
		Prompt peakShape, "shape for fitted peak",popup,"Lorentzian;Gaussian"
		Prompt smoothing, "smooth before fitting",popup,"NO smoothing;Smooth"
		Prompt thresholdRatio, "threshold = (ratio*[std dev] + avg)"
		Prompt maxNum,"Max number of peaks to fit, use -1 for unlimited"
		smoothing += 1
		DoPrompt/Help="3D-Xray Diffraction[Fit Peaks]" "peak fitting",imageName,minPeakWidth,boxSize,maxRfactor,threshAboveAvg,thresholdRatio,peakShape,smoothing,maxNum
		if (V_flag)
			return ""
		endif
		smoothing = (smoothing==2)
		Wave image = $imageName
		printIt = 2
	endif
	if (printIt>1)
		sprintf line,"FitPeaksWithExternal(%s, %g, %g, %g, %g, peakShape=\"%s\"",NameOfWave(image),minPeakWidth,boxSize,maxRfactor,threshAboveAvg,peakShape
		if (ParamIsDefault(smoothing) && smoothing!=smoothingLast || smoothing || !ParamIsDefault(smoothing))
			line += ", smoothing="+num2str(smoothing)
		endif
		if (numtype(thresholdRatio)==0 || !ParamIsDefault(thresholdRatio))
			line += ", thresholdRatio="+num2str(thresholdRatio)
		endif
		if (numtype(maxNum)||maxNum<1)
			line += ", maxNum="+num2str(maxNum)
		endif
		if (strlen(GetRTStackInfo(2)))
			line = "•"+line
		endif
	endif
	if (!WaveExists(image) || WaveDims(image)!=2)
		if (printIt && strlen(line))
			printf "%s)\r",line
		endif
		return ""
	endif
	if (!(threshAboveAvg>0))
		threshAboveAvg = NaN
	endif
	if (!(minPeakWidth>0) || !(boxSize>0))
		if (printIt && strlen(line))
			printf "%s)\r",line
		endif
		return ""
	endif
	imageName = GetWavesDataFolder(image,2)
	Variable Nx=DimSize(image,0), Ny=DimSize(image,1)
	if (Nx<2 || Ny<2)
		if (printIt && strlen(line))
			printf "%s)\r",line
		endif
		return ""
	endif
	maxNum = maxNum<1 || numtype(maxNum) ? -1 : maxNum

	String maskName=""
	if (!WaveExists(mask) && printIt)
		String options
		sprintf options "BYTE:1,UNSIGNED:1,MINROWS:%d,MAXROWS:%d,MINCOLS:%d,MAXCOLS:%d,",Nx,Nx,Ny,Ny
		String maskList = WaveList("*",";",options)
		if (itemsInList(maskList)>0)
			maskName=maskNameLast
			Prompt maskName,"mask to use",popup,"_none_;"+maskList
			DoPrompt "mask",maskName
			if (V_flag)
				if (printIt && strlen(line))
					printf "%s)\r",line
				endif
				return ""
			endif
			Wave mask=$maskName
		endif
	endif
	if (WaveExists(mask))										// if the mask wave exists
		if (DimSize(mask,0)!=Nx || DimSize(mask,1)!=Ny)	// check that mask exists is right size
			DoAlert 0,"Requested mask is not same size as image wave"
			if (printIt && strlen(line))
				printf "%s)\r",line
			endif
			return ""
		endif
		line += SelectString(strlen(line), "",", mask="+NameOfWave(mask))
		maskName = NameOfWave(mask)
	else
		maskName = ""
	endif
	if (strlen(line))
		printf "%s)\r",line
	endif
	minPeakWidthLast = minPeakWidth							// save current values for next time
	boxSizeLast = boxSize
	maxRfactorLast = maxRfactor
	threshAboveAvgLast = threshAboveAvg
	smoothingLast = smoothing
	thresholdRatioLast = thresholdRatio
	maxNumLast = maxNum
	peakShapeLast = peakShape
	maskNameLast = maskName
	Variable timer=startMSTimer
	Wave FullPeakList = $runPeakSearchCommand(image,boxSize,maxRfactor,minPeakWidth,threshAboveAvg,maxNum,peakShape=peakShape,smoothing=smoothing,thresholdRatio=thresholdRatio,mask=mask)
	Variable sec=stopMSTimer(timer)/1e6
	if (!WaveExists(FullPeakList))
		return ""
	endif
	Variable Nu = DimSize(FullPeakList,0)

	String win = StringFromList(0,WindowsWithWave(image,1))
	if (strlen(win))
		SetWindow $win userdata(FitPeaks) = "FullPeakList="+GetWavesDataFolder(FullPeakList,2)
	endif
	if (printIt)
		threshAboveAvg = threshAboveAvg>0 ? threshAboveAvg : NumberBykey("threshAboveAvg",note(FullPeakList),"=")
		printf "found and fitted %d peaks (using a threshold of %.1f),  this all took %g sec\r",Nu,threshAboveAvg,sec
		printf "  result stored in the wave '%s'\r", GetWavesDataFolder(FullPeakList,2)
	endif
	return GetWavesDataFolder(FullPeakList,2)
End
//
//#if stringmatch(igorInfo(2),"Macintosh") //commented RX
//
// using the result from FitPeaks() run peaksearch, and read in the results from output file
Static Function/S runPeakSearchCommand(image,boxSize,maxRfactor,minSize,threshold,maxNum,[peakShape,smoothing,thresholdRatio,mask])
	Wave image								// image to use
	Variable boxSize
	Variable maxRfactor
	Variable minSize
	Variable threshold
	Variable maxNum						// maximum number of peaks to fit (defaults to unlimited)
	String peakShape
	Variable smoothing						// if TRUE, do a smoothing operation before peak search
	Variable thresholdRatio					// set threshold to (ratio*[std dev] + avg) if threshold passed as NaN (optional)
	Wave mask								// wave for mask, use pixels with mask==0

	Variable badNums = !(boxSize>0) || !(maxRfactor>0) || !(minSize>0) || numtype(smoothing) || thresholdRatio<0 || numtype(thresholdRatio)==1
	if (badNums)
		DoAlert 0, "Invalid inputs sent to runPeakSearchCommand()"
		printf "runPeakSearchCommand(%s,%g,%g,%g,%g",NameOfWave(image),boxSize,maxRfactor,minSize,threshold
		if (!ParamIsDefault(peakShape))
			printf ",peakShape=%s",peakShape
		endif
		if (!ParamIsDefault(smoothing))
			printf ",smoothing=%g",smoothing
		endif
		if (!ParamIsDefault(thresholdRatio))
			printf ",thresholdRatio=%g",thresholdRatio
		endif
		if (WaveExists(mask))
			printf ",mask=%s",NameOfWave(mask)
		endif
		printf ")\r"
		return ""
	endif

	String upath=SpecialDirPath("Temporary",0,1,0)	// local path (probably unix path) for the command line
	String mpath=SpecialDirPath("Temporary",0,0,0)	// mac style path for use only within Igor
	String wnote = note(image)
	String imageFile=StringByKey("imageFilePath",wnote,"=")+StringByKey("imageFileName",wnote,"=")	// file with image (full path name)
	GetFileFolderInfo/Q/Z=1 imageFile				// check that image file exists
	Variable deleteImage = 0						// flag is true means delete image afterward
	if (V_Flag || !V_isFile)							// could not find original image file, try to save this one from Igor
		imageFile = ""
#if (Exists("Write1HDF5imageFile")==6)
		imageFile = Write1HDF5imageFile(image,mpath+"peakFitFile")	// write image to temporary file
		deleteImage = strlen(imageFile)>1
#endif
	endif
	if (strlen(imageFile)<1)
		DoAlert 0, "image file '"+imageFile+"' does not exist"
		return ""
	endif
	String imageFileMac=imageFile
	Variable isMac = stringmatch(igorInfo(2),"Macintosh")
	Variable isWin = stringmatch(igorInfo(2), "Windows")
	if (isMac) // RX, JZT
		imageFile = ParseFilePath(5,imageFile,"/",0,0)	// on Mac only convert HFS to POSIX
	endif	
	// first write the command file to drive peaksearch program
	PathInfo home
	if (V_flag && !isWin)						// the path "home" exists, use it (but not on windows, to avoid certain permission issue since 2018)
		mpath = S_path
		if (isMac)
			upath = ParseFilePath(5,S_path,"/",0,0)
		endif
	endif
	NewPath/O/Q/Z PeakSearchCalcFolder, mpath	// need a new path name since "home" may not exist
	if (V_flag)
		DoAlert 0, "Unable to create path to run peaksearch, try saving the experiment first."
		return ""
	endif
	String peakFile="generic_XY.txt"			// name of file to receive result

	// find the full path name of the peaksearch executable
	String name,PeakSearchPath=getBinaryPath()	// path to the peaksearch executables

	if (isMac && stringmatch(igorinfo(4),"PowerPC"))	// find the default name for this architecture
		name = "peaksearch_ppc"
	elseif (isMac && stringmatch(igorinfo(4),"Intel"))
		name = "peaksearch_i386"
	elseif(isWin)
		name = "peaksearch_win32.exe" // RX
	else
		name = "peaksearch"
	endif
	name = StrVarOrDefault("root:Packages:micro:PeakFit:PeakSearchExecutableMac",name)	// override default if PeakSearchExecutableMac is set
	GetFileFolderInfo/Q/Z PeakSearchPath+name
	PeakSearchPath += SelectString(V_Flag==0 && V_isFile,"peaksearch",name)	// use just plane peaksearch if peaksearch_ppc or peaksearch_i386 does not exist
	if (isMac)												// RX, JZT,  only on Mac, convert paths from HFS to POSIX
		PeakSearchPath = ParseFilePath(5,PeakSearchPath,"/",0,0)
	endif
	if (strlen(PeakSearchPath)<1)
		DoAlert 0, "cannot find the executable '"+name+"'"
		return ""
	endif

	String maskFileName=""
	SVAR maskHashLast = root:Packages:micro:PeakFit:maskHashLast
	if (WaveExists(mask))								// a mask file was passed, write it to disk, (and do not delete it afterwards)
		if (Exists("Write1HDF5imageFile")!=6)
			DoAlert 0,"External peak fitting works ONLY with hdf5 images"
			return ""
		endif
		String maskHash="NOT"+maskHashLast
		GetFileFolderInfo/P=PeakSearchCalcFolder/Q/Z=1 "generic_mask.h5"		// find out if existing file is the one I want
		if (V_Flag==0 && V_isFile && !V_isAliasShortcut && (modDate(mask)-V_modificationDate)<1)
			maskHash = Hash(S_Path+S_fileType+num2istr(V_modificationDate),1)
		endif
		if (cmpstr(maskHash,maskHashLast,1))			// they differ, write the mask to disk
		#if (Exists("Write1HDF5imageFile")==6)		// only works for HDF5 files
			maskFileName = Write1HDF5imageFile(mask,mpath+"generic_mask.h5")	// name of file to hold mask
		#endif
			if (isMac) 											// RX, JZT		only on Mac, need POSIX paths
				maskFileName = ParseFilePath(5,maskFileName,"/",0,0)	
			endif
			if (strlen(maskFileName)<1)
				DoAlert 0, "Unable to write mask image file to disk"
				return ""		
			endif
			GetFileFolderInfo/P=PeakSearchCalcFolder/Q/Z=1 "generic_mask.h5"
			maskHashLast = Hash(S_Path+S_fileType+num2istr(V_modificationDate),1)	// update last hash
		else
			maskFileName = upath+"generic_mask.h5"	// just use existing file
		endif
	endif
	String cmd,params
	sprintf params "-b %g -R %g -m %g",boxSize,maxRfactor,minSize
	if (threshold>0)
		params += " -t "+num2str(threshold)
	endif
	params += SelectString(smoothing,""," -S")
	params += SelectString(thresholdRatio>0,""," -T "+num2str(thresholdRatio))
	if (maxNum>1)
		params += " -M "+num2str(maxNum)
	endif
	if (!ParamIsDefault(peakShape) && strlen(peakShape)>0)
		params += " -p "+peakShape
	endif
	if (strlen(maskFileName))
		params += " -K \\\""+maskFileName+"\\\""
	endif
	if (isMac)
		sprintf cmd "do shell script \"cd \\\"%s\\\" ; \\\"%s\\\" %s  \\\"%s\\\"  %s\"",upath,PeakSearchPath,params,imageFile,peakFile
		ExecuteScriptText cmd
	elseif (isWin) // create a .bat file //RX
		variable BatchFileRN
		variable OutputFileRN
		String BatchFileName = "PeakSearchExecute.bat"
		String FullBatchFileName = mpath + BatchFileName
		String WinPeakSearchPath = MacToWindows(PeakSearchPath)
		String WinimageFilePath = MacToWindows(imageFile)
		
		//String peakFile = MacToWindows(mpath + peakFile)
		//String FullPeakFileName = MacToWindows(mpath + peakFile)
		Open /Z BatchFileRN  as FullBatchFileName
		variable errflag=V_flag
		if (errflag != 0)
			printf "Failed to open %s for writing", FullBatchFileName
			return ""
		endif
		fprintf BatchFileRN, "%s:\r\n", ParseFilePath(0,mpath,":",0,0)   // drive letter + ":"
		fprintf BatchFileRN, "cd \"%s\"\r\n", MacToWindows(mpath)
		fprintf BatchFileRN, "\"%s\" %s \"%s\" \"%s\"\r\n", WinPeakSearchPath,params,WinimageFilePath,peakFile
		close BatchFileRN
		sprintf cmd "\"%s\"", MacToWindows(FullBatchFileName)
		ExecuteScriptText cmd
		open /R/Z OutputFileRN as peakFile
		if (V_flag != 0)
			print "Error opening peakSearch batch output file ", peakFile
		else
			Close OutputFileRN
		endif
	else   // prompt & exit if other OS
		printf "External peak search only runs on Mac or Windows."
		return ""
	endif
	if (ItemsInList(GetRTStackInfo(0))<2)		// print everything if run from command line
		print S_value
	endif
	if (deleteImage)
		DeleteFile/Z=1 imageFileMac
		if (V_flag)
			print "Unable to delete temporary file '"+imageFileMac+"'"
			DoAlert 0,"Unable to delete temporary file '"+imageFileMac+"', this should never happen"
		endif
	endif
	if (strlen(S_value)>2)
		DoAlert 0, "failure in runPeakSearchCommand()"
		print "\r\r"
		print cmd
		print "\r\r"
		print S_value
		return ""								// there is no peaksearch output file
	endif
	wnote = ReplaceStringByKey("fittedIgorImage",wnote,GetWavesDataFolder(image,2),"=")
	return readPeakXYfile(peakFile,"PeakSearchCalcFolder",wnote)	// read peaksearch results and create an array to hold them, return array name
End
//
// begin block commented RX:
//#elif stringmatch(igorInfo(2),"Windows")			// Windows version of runIndexingEulerCommand()
// //
// //	This is just a stub until I get a final version
// // using the result from FitPeaks() run peaksearch, and read in the results from output file
// Static Function/S runPeakSearchCommand(image,boxSize,maxRfactor,minSize,threshold,maxNum,[peakShape,smoothing,thresholdRatio,mask])
	// Wave image								// image to use
	// Variable boxSize
	// Variable maxRfactor
	// Variable minSize
	// Variable threshold
	// Variable maxNum						// maximum number of peaks to fit (defaults to unlimited)
	// String peakShape
	// Variable smoothing						// if TRUE, do a smoothing operation before peak search
	// Variable thresholdRatio					// set threshold to (ratio*[std dev] + avg) if threshold passed as NaN (optional)
	// Wave mask								// wave for mask, use pixels with mask==0
	// DoAlert 0, "runPeakSearchCommand() not yet implemented for Windows"
	// printf "runPeakSearchCommand() not yet implemented for Windows\r"
	// return ""
// End
// //
//#endif //end block commented RX
//
Static Function/S readPeakXYfile(peakFile,path,wnote)
	String peakFile					// name of index file, the output from Euler
	String path							// name of Igor path to go with peakFile
	String wnote						// info from image file

	Variable f = OpenFileOf_ftype(peakFile,"PixelPeakList",path)// open a peaksearch file
	if (f<1)
		return ""						// could not get any info from indexFile 
	endif

	FStatus f
	String buffer=""
	buffer = PadString(buffer,V_logEOF-V_filePos,0x20)
	FBinRead f, buffer				// read entire file into buffer
	Close f
	buffer = ReplaceString("\r\n",buffer,"\n")				// ensure that line term is a new-line only
	buffer = ReplaceString("\n\r",buffer,"\n")
	buffer = ReplaceString("\r",buffer,"\n")
	Variable Nbuf = strlen(buffer)

	Variable i0=strsearch(buffer,"\n",0)+1, i1=strsearch(buffer,"$peakList",0)-1
	i1 = strsearch(buffer,"\n",i1+3)+1					// extend i1 to include the $peakList line too
	if (i0<2 || i1<10)
		return ""
	endif
	String list = keyStrFromBuffer(buffer[i0,i1])
	Variable Npeaks = NumberByKey("Npeaks",list,"=")		// number of fitted peaks in following table
	Npeaks = Npeaks>0 ? Npeaks : 0
	Variable Ncols = str2num(StringByKey("peakList",list,"="))	// number of columns
	Ncols = Ncols==5 ? 4 : Ncols							// ***** special to fix old error *****

	wnote = ReplaceStringByKey("waveClass",wnote,"FittedPeakList","=")
	wnote = ReplaceNumberByKey("boxsize",wnote,NumberByKey("boxsize",list,"="),"=")				// box size used for peak fitting
	wnote = ReplaceNumberByKey("maxRfactor",wnote,NumberByKey("maxRfactor",list,"="),"=")		// max allowed R-factor
	wnote = ReplaceNumberByKey("minPeakWidth",wnote,NumberByKey("minwidth",list,"="),"=")		// min allowed width of a peak
	wnote = ReplaceNumberByKey("maxPeakWidth",wnote,NumberByKey("maxwidth",list,"="),"=")	// max allowed width of a peak
	wnote = ReplaceNumberByKey("threshAboveAvg",wnote,NumberByKey("threshold",list,"="),"=")	// threshold for blob searching
	wnote = ReplaceNumberByKey("minSpotSeparation",wnote,NumberByKey("minSeparation",list,"="),"=")		// minimum separation between any two peaks
	wnote = ReplaceNumberByKey("totalIntensity",wnote,NumberByKey("totalSum",list,"="),"=")			// sum of all pixels in image
	wnote = ReplaceNumberByKey("totalPeakIntensity",wnote,NumberByKey("sumAboveThreshold",list,"="),"=")	// sum of all pixels in image above threshold
	wnote = ReplaceNumberByKey("NumPixelsAboveThreshold",wnote,NumberByKey("numAboveThreshold",list,"="),"=")	// number of pixels above threshold
	string peakShape = StringByKey("peakShape",list,"=")
	if (!stringmatch(peakShape,"Lorentzian"))
		wnote = ReplaceStringByKey("FittedPeakShape",wnote,StringByKey("peakShape",list,"="),"=")		// shape for peak fit
	endif
	Variable ismooth = NumberByKey("smooth",list,"=")
	if (ismooth)
		wnote = ReplaceNumberByKey("smoothed",wnote,ismooth,"=")	// fit to smoothed image
	endif
	String peakListName="FullPeakList"+detectorID2color(StringByKey("detectorID",wnote,"="))
	Make/N=(Npeaks,12)/O $peakListName/WAVE=FullPeakList = NaN
	SetDimLabel 1,0,x0,FullPeakList			;	SetDimLabel 1,1,y0,FullPeakList
	SetDimLabel 1,2,x0Err,FullPeakList		;	SetDimLabel 1,3,y0Err,FullPeakList
	SetDimLabel 1,4,fwx,FullPeakList		;	SetDimLabel 1,5,fwy,FullPeakList
	SetDimLabel 1,6,fwxErr,FullPeakList	;	SetDimLabel 1,7,fwyErr,FullPeakList
	if (Ncols==8)
		SetDimLabel 1,8,tilt_degree,FullPeakList	;	SetDimLabel 1,9,chisq,FullPeakList
	else
		SetDimLabel 1,8,correlation,FullPeakList	;	SetDimLabel 1,9,correlationErr,FullPeakList
	endif
	SetDimLabel 1,10,area,FullPeakList					;	SetDimLabel 1,11,amp,FullPeakList
	Note/K FullPeakList,wnote
	if (Npeaks<1)
		return GetWavesDataFolder(FullPeakList,2)
	endif

	i0 = i1
	i1 = strsearch(buffer,"\n",i0)
	Variable px,py,maxIntens,integral
	Variable hwhmX,hwhmY,tilt,chisq							// new values added to peakSearch output
	Variable m=0
	do
		if (Ncols<8)
			sscanf buffer[i0,i1],"%g %g %g %g",px,py,maxIntens,integral
			FullPeakList[m][0] = px
			FullPeakList[m][1] = py
			FullPeakList[m][10] = integral					// actual integral of the peak (from image, not fit func)
		else
			sscanf buffer[i0,i1],"%g %g %g %g %g %g %g %g",px,py,maxIntens,integral,hwhmX,hwhmY,tilt,chisq
			FullPeakList[m][0] = px
			FullPeakList[m][1] = py
			FullPeakList[m][10] = integral
			FullPeakList[m][4] = 2*hwhmX
			FullPeakList[m][5] = 2*hwhmY
			FullPeakList[m][8] = tilt
			FullPeakList[m][9] = chisq						// this is not really right, but do it anyhow
			FullPeakList[m][11] = maxIntens					// amplitude of the fit function
		endif
		i0 = i1+1
		i1 = strsearch(buffer,"\n",i0)
		m += 1
	while(m<Npeaks && V_flag==Ncols && i1>i0)
	return GetWavesDataFolder(FullPeakList,2)
End












// First get a background:
//		1) down sample the image by 2 (could be 4)
//		2) do a median filter to get a background level for each pixel
//		3) gaussian filter the background to smooth it
//		4) upsample to get a background the correct size
//	Subtract this background from the original image
//	and median filter (only 5) and gauss smooth to get rid of specks
//	Thresold this image to generate a mask of regions that represent the peak locations (an ROI)
//	Examine each region in the ROI separately, fitting a gaussian
Function/S FitPeaksWithSeedFill(image,minPeakWidth,maxPeakWidth,minSep,threshAboveAvg,[mask,maxNum,whoami,printIt])
	Wave image						// the image from the image file
	Variable minPeakWidth		// minimum width for an acceptable peak (pixels)
	Variable maxPeakWidth		// maximum width for an acceptable peak (pixels)
	Variable minSep				// minimum separation between two peaks (pixels)
	Variable threshAboveAvg	// threshold above average value,  use 25
	Wave mask						// starting mask for the fit (this mask is unchanged by this routine)
	Variable maxNum				// maximum number of peaks to find, normally goes to completion
	Variable whoami				// if passed, then just return the name of this routine
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt
	if (!ParamIsDefault(whoami))
		return GetRTStackInfo(1)
	endif

	if (exists("root:Packages:micro:PeakFit:minPeakWidthLast")!=2)
		indexing#initPeaksPart()
	endif
	NVAR minPeakWidthLast = root:Packages:micro:PeakFit:minPeakWidthLast
	NVAR maxPeakWidthLast = root:Packages:micro:PeakFit:maxPeakWidthLast
	NVAR minSepLast = root:Packages:micro:PeakFit:minSepLast
	NVAR threshAboveAvgLast = root:Packages:micro:PeakFit:threshAboveAvgLast
	NVAR maxNumLast = root:Packages:micro:PeakFit:maxNumLast
	SVAR maskNameLast = root:Packages:micro:PeakFit:maskNameLast
	maxNum = ParamIsDefault(maxNum) ? Inf : maxNum
	if (ParamIsDefault(mask))
		Wave mask = $""
	endif

	String imageName, line=""
	if (!WaveExists(image) || WaveDims(image)!=2 || !(minPeakWidth>0) || !(maxPeakWidth>0) || !(minSep>3) || !(threshAboveAvg>0))
		imageName = SelectString(WaveExists(image),"",NameOfWave(image))
		minPeakWidth = minPeakWidth>.5 ? minPeakWidth : minPeakWidthLast
		maxPeakWidth = maxPeakWidth>.5 ? maxPeakWidth : maxPeakWidthLast
		minSep = minSep>minPeakWidth ? minSep : minSepLast
		threshAboveAvg = threshAboveAvg>0 ? threshAboveAvg : threshAboveAvgLast
		Prompt imageName,"name of image with peaks to find",popup,reverseList(WaveListClass("speImageNoBkg;speImage;rawImageNoBkg;rawImage*","*","DIMS:2"))
		Prompt minPeakWidth,"minimum allowed width of a peak"
		Prompt maxPeakWidth,"maximum allowed width of a peak"
		Prompt minSep,"minimum distance between two peaks, reject peaks closer than this"
		Prompt threshAboveAvg,"min peak height, use NaN for auto detect"
		DoPrompt/Help="3D-Xray Diffraction[Fit Peaks]" "peak fitting",imageName,minPeakWidth,maxPeakWidth,minSep,threshAboveAvg
		if (V_flag)
			return ""
		endif
		Wave image = $imageName
		// printf "FitPeaksWithSeedFill(%s,%g,%g,%g,%g)\r",imageName,minPeakWidth,maxPeakWidth,minSep,threshAboveAvg
		sprintf line, "FitPeaksWithSeedFill(%s,%g,%g,%g,%g",imageName,minPeakWidth,maxPeakWidth,minSep,threshAboveAvg
		line = SelectString(strlen(GetRTStackInfo(2)),"","•")+line
		printIt = 1
	endif
	if (!WaveExists(image) || WaveDims(image)!=2)
		if (strlen(line))
			printf "%s)\r",line
		endif
		return ""
	endif
	if (!(threshAboveAvg>0))
		if (strlen(StringByKey("depthSi",note(image),"=")))
			// threshAboveAvg = 10
			threshAboveAvg = noiseInImage(image)/2
		else
			threshAboveAvg = StatsMedian(image)/4 
			// !(threshAboveAvg>0) = !(threshAboveAvg>0) ? StatsMedian(image)/4 : threshAboveAvg
		endif
		if (printIt)
			printf "   using a threshAboveAvg = %.3g\r",threshAboveAvg
		endif
	endif
	if (!(minPeakWidth>0) || !(maxPeakWidth>0) || !(minSep>3) || !(threshAboveAvg>0))
		if (strlen(line))
			printf "%s)\r",line
		endif
		return ""
	endif
	imageName = GetWavesDataFolder(image,2)
	Variable Nx=DimSize(image,0), Ny=DimSize(image,1)
	if (Nx<2 || Ny<2)
		if (strlen(line))
			printf "%s)\r",line
		endif
		return ""
	endif

	String maskName=""
	if (WaveExists(mask) && strlen(line))
		line += ", mask="+NameOfWave(mask)
	endif
	if (!WaveExists(mask) && printIt)
		String options
		sprintf options "BYTE:1,UNSIGNED:1,MINROWS:%d,MAXROWS:%d,MINCOLS:%d,MAXCOLS:%d,",Nx,Nx,Ny,Ny
		String maskList = WaveList("*",";",options)
		if (itemsInList(maskList)>0)
			maskName=maskNameLast
			Prompt maskName,"mask to use",popup,"_none_;"+maskList
			DoPrompt "mask",maskName
			if (V_flag)
				if (strlen(line))
					printf "%s)\r",line
				endif
				return ""
			endif
			Wave mask=$maskName
			if (WaveExists(mask))
				line += ", mask="+maskName
			endif
		endif
	endif
	if (WaveExists(mask))										// if the mask wave exists
		if (DimSize(mask,0)!=Nx || DimSize(mask,1)!=Ny)	// check that mask exists is right size
			DoAlert 0,"Requested mask '"+NameOfWave(mask)+"'is not same size as image wave"
			if (strlen(line))
				printf "%s)\r",line
			endif
			return ""
		endif
		maskName = NameOfWave(mask)
	else
		maskName = ""
	endif
	minPeakWidthLast = minPeakWidth
	maxPeakWidthLast = maxPeakWidth
	minSepLast = minSep
	threshAboveAvgLast = threshAboveAvg
	maxNumLast = maxNum
	maskNameLast = maskName
	if (strlen(line))
		printf "%s)\r",line
	endif

	Variable timer=startMSTimer
	Variable Nfilt=20
	Variable downSample = NumVarOrDefault("downSample", 2 )					// could be 4 (don't be greedy),  1 also works but is slow

	Duplicate/FREE image, imageFit,imageS						// image to use when finding next hot pixel
	Redimension/D imageS, imageFit

	ImageInterpolate/PXSZ={(downSample),(downSample)} Pixelate imageS			// down sample the image to get background
	Wave M_PixelatedImage=M_PixelatedImage

//	Variable Ndown = Nfilt/downSample
	Variable Ndown = max(7,round(Nfilt/downSample))
	MatrixFilter/N=(Ndown)/M=(Ndown*Ndown*0.5) median M_PixelatedImage
	MatrixFilter/N=(2*Nfilt) gauss M_PixelatedImage								// this is now the background

	ImageInterpolate/F={(downSample),(downSample)} bilinear M_PixelatedImage	// up sample the background
	Wave M_InterpolatedImage=M_InterpolatedImage
	Variable rowsToAdd = DimSize(image,0)-DimSize(M_InterpolatedImage,0)		// redimension to exactly the right size
	Variable colsToAdd = DimSize(image,1)-DimSize(M_InterpolatedImage,1)
	ImageTransform/N={(rowsToAdd),(colsToAdd)}/O padimage M_InterpolatedImage

	ImageStats M_InterpolatedImage
	Variable bkg = V_rms
	MatrixOp/O/FREE imageS = image - M_InterpolatedImage		// imageS = image - M_InterpolatedImage
	MatrixOp/O/FREE imageFit = image - M_InterpolatedImage	// imageFit = image - M_InterpolatedImage
	MatrixFilter/N=5 median imageS
	MatrixFilter/N=10 gauss imageS

	Duplicate/FREE imageS,FitPeaks_ImageMask
	ImageThreshold/T=(threshAboveAvg)/I/O FitPeaks_ImageMask
	if (WaveExists(mask))										// if a mask was given, apply it here
		FitPeaks_ImageMask = FitPeaks_ImageMask | mask
	endif

	Duplicate/FREE imageS weights, ones						// weights to use in the peak fitting
	Redimension/B/U ones
	ones = 1

	String str
	Make/O/T JZT_Constraints={"K3 > 0.1","K5 > 0.1","K6 > 0","K6 < 1","","","",""}
	sprintf str,"K3 > %g",minPeakWidth/100
	JZT_Constraints[0] = str
	sprintf str,"K5 > %g",minPeakWidth/100
	JZT_Constraints[1] = str
	sprintf str,"K2 > %g",-maxPeakWidth
	JZT_Constraints[4] = str
	sprintf str,"K2 < %g",DimOffset(image,0) + (DimSize(image,0)-1)*DimDelta(image,0) + maxPeakWidth
	JZT_Constraints[5] = str
	sprintf str,"K4 > %g",-maxPeakWidth
	JZT_Constraints[6] = str
	sprintf str,"K4 < %g",DimOffset(image,1) + (DimSize(image,1)-1)*DimDelta(image,1) + maxPeakWidth
	JZT_Constraints[7] = str
//	sprintf str,"K1 > %g",-100*threshAboveAvg
//	JZT_Constraints[8] = str

	Variable Nu=0, Nlen=50
	String wnote = note(image)
	String peakListName="FullPeakList"+detectorID2color(StringByKey("detectorID",wnote,"="))
	Make/N=(Nlen,12)/O $peakListName/WAVE=FullPeakList = NaN
	wnote = ReplaceNumberByKey("minPeakWidth",wnote,minPeakWidth,"=")
	wnote = ReplaceNumberByKey("maxPeakWidth",wnote,maxPeakWidth,"=")
	wnote = ReplaceNumberByKey("minSpotSeparation",wnote,minSep,"=")
	wnote = ReplaceNumberByKey("threshAboveAvg",wnote,threshAboveAvg,"=")
	wnote = ReplaceStringByKey("fittedIgorImage",wnote,GetWavesDataFolder(image,2),"=")
	wnote = ReplaceStringByKey("waveClass",wnote,"FittedPeakList","=")
	SetDimLabel 1,0,x0,FullPeakList			;	SetDimLabel 1,1,y0,FullPeakList
	SetDimLabel 1,2,x0Err,FullPeakList		;	SetDimLabel 1,3,y0Err,FullPeakList
	SetDimLabel 1,4,fwx,FullPeakList		;	SetDimLabel 1,5,fwy,FullPeakList
	SetDimLabel 1,6,fwxErr,FullPeakList	;	SetDimLabel 1,7,fwyErr,FullPeakList
	SetDimLabel 1,8,correlation,FullPeakList ;	SetDimLabel 1,9,correlationErr,FullPeakList
	SetDimLabel 1,10,area,FullPeakList		;	SetDimLabel 1,11,amp,FullPeakList
	Note/K FullPeakList,wnote

	Variable fw= ( 2*sqrt(2*ln(2)) )							// this factor corrects for 2-d sigma to FWHM, apply to K3, K5, and the corresponding errors
	Variable xx,yy,fwx,fwy,sigX,sigY, amp						// holds the resutls of gaussian fit
	Variable px,py, left,right,top,bot, err, i
	do
		// DoUpdate
		ImageStats/M=1/R=FitPeaks_ImageMask imageS		// find the highest point
		if (V_max<threshAboveAvg || V_npnts<2)
			break
		endif
		px=V_maxRowLoc ; py=V_maxColLoc					// look for a peak at (px,py)
		// fit a gaussian about px,py
		left = limit(floor(px-maxPeakWidth),0,Nx-1)		// box to use for fitting
		right = limit(ceil(px+maxPeakWidth),0,Nx-1)
		top = limit(floor(py-maxPeakWidth),0,Ny-1)
		bot = limit(ceil(py+maxPeakWidth),0,Ny-1)
		ImageSeedFill seedP=px,seedQ=py,min=0,max=0,target=1,srcwave=FitPeaks_ImageMask
		Wave M_SeedFill=M_SeedFill
		MatrixOp/O/FREE weights = equal(M_SeedFill,ones)	// weights = M_SeedFill==1 ? 1 : 0
		err = FitGaussianPeak(imageFit,left,right,top,bot,weight=weights)	// returns 0 == OK,  non-zero if error
		ImageSeedFill/O seedP=px,seedQ=py,min=0,max=0,target=255,srcwave=FitPeaks_ImageMask	// fil in this region, already used
		if (err)
			continue
		endif
		Wave W_coef=W_coef, W_sigma=W_sigma
		xx = W_coef[2]
		yy = W_coef[4]
		sigX = W_sigma[2]
		sigY = W_sigma[4]
		fwx = abs(W_coef[3])*fw
		fwy = abs(W_coef[5])*fw
		amp = abs(W_coef[1])

//		if (ItemsInList(GetRTStackInfo(0))<2)
//			printf "%d    Gaussian peak in at (%.2f±%.2f, %.2f±%.2f) with FWHM of (%.2f±%.2f, %.2f±%.2f),   amp=%g,  image[x,y]= %g\r",Nu,xx,sigX,yy,sigY,fwx,W_sigma[3]*fw,fwy,W_sigma[5]*fw,amp,image[xx][yy]
//		endif
		err = ( sigX>maxPeakWidth || sigY>maxPeakWidth || numtype(sigX+sigY) )						// errors bars are huge
		if (!(min(fwx,fwy)>minPeakWidth) || !(max(fwx,fwy)<maxPeakWidth) || !(max(amp,image[xx][yy])>bkg) || err)	// cannot have spots too wide or too narrow
			continue
		endif
		for (i=0;i<Nu;i+=1)
			if (((xx-FullPeakList[i][0])^2 + (yy-FullPeakList[i][1])^2) < minsep^2)
				i = -1
				break
			endif
		endfor
		if (i<0)
			continue
		endif

		if (Nu>=Nlen)												// FullPeakList is too short, extend it
			Nlen += 50
			Redimension/N=(Nlen,-1) FullPeakList
		endif
		FullPeakList[Nu][0]=xx				;	FullPeakList[Nu][1]=yy
		FullPeakList[Nu][2]=sigX				;	FullPeakList[Nu][3]=sigY
		FullPeakList[Nu][4]=fwx				;	FullPeakList[Nu][5]=fwy
		FullPeakList[Nu][6]=W_sigma[3]*fw	;	FullPeakList[Nu][7]=W_sigma[5]*fw
		FullPeakList[Nu][8]=W_coef[6]		;	FullPeakList[Nu][9]=W_sigma[6]
		FullPeakList[Nu][10]=W_coef[1]*fwx*fwy
		FullPeakList[Nu][11]=W_coef[1]
		Nu += 1
	while(Nu<maxNum)
	Redimension/N=(Nu,-1) FullPeakList						// set to exact length
	ImageStats/M=1/G={0,Nu-1,10,10}/Q FullPeakList
	wnote = ReplaceNumberByKey("totalPeakIntensity",wnote,V_avg*V_npnts,"=")	// total intensity in all fitted peaks
	wnote = ReplaceNumberByKey("totalIntensity",wnote,sum(image),"=")				// total intensity in the image
	Note/K FullPeakList,wnote

	CheckDisplayed image
	if (V_flag)
		SetWindow kwTopWin userdata(FitPeaks) = "FullPeakList="+GetWavesDataFolder(FullPeakList,2)
	endif

	Variable sec=stopMSTimer(timer)/1e6
	if (printIt)
		printf "found and fitted %d peaks (using a threshold of %.1f),  this all took %g sec\r",Nu,threshAboveAvg,sec
		printf "  result stored in the wave '%s'\r", GetWavesDataFolder(FullPeakList,2)
	endif
	KillWaves/Z M_SeedFill, M_PixelatedImage,M_InterpolatedImage
	KillWaves/Z JZT_Constraints
	KillWaves/Z W_coef, W_sigma
	return GetWavesDataFolder(FullPeakList,2)
End

Function noiseInImage(image)
	Wave image
	if (!WaveExists(image))
		return NaN
	endif
	Duplicate/O/I/FREE image, diff_noiseInImage
	Wave diff = diff_noiseInImage
	diff = image[p+1][q] - image[p][q]
	WaveStats/Q diff

	if (ItemsInList(GetRTStackInfo(0))<2)
		print "V_sdev =",V_sdev
	endif
	return V_sdev
End


// identify the peaks, and fit them all (no background removal), a lot like FitPeaksNew(), but proceeds in a more step wise fashion
Function/S FitPeaksStepWise(image,minPeakWidth,maxPeakWidth,minSep,threshAboveAvg,[mask,maxNum,whoami,printIt])
	Wave image						// the image from the image file
	Variable minPeakWidth			// minimum width for an acceptable peak (pixels)
	Variable maxPeakWidth			// maximum width for an acceptable peak (pixels)
	Variable minSep				// minimum separation between two peaks (pixels)
	Variable threshAboveAvg		// threshold above average value
	Wave mask						// starting mask for the fit (this mask is unchanged by this routine)
	Variable maxNum				// maximum number of peaks to find, normally goes to completion
	Variable whoami				// if passed, then just return the name of this routine
	Variable printIt
	if (!ParamIsDefault(whoami))
		return GetRTStackInfo(1)
	endif
	maxNum = ParamIsDefault(maxNum) ? Inf : maxNum
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt

	String imageName
	if (!WaveExists(image) || WaveDims(image)!=2 || !(minPeakWidth>0) || !(maxPeakWidth>0) || !(minSep>3) || !(threshAboveAvg>0))
		imageName = SelectString(WaveExists(image),"",NameOfWave(image))
		minPeakWidth = minPeakWidth>.5 ? minPeakWidth : 1.2
		maxPeakWidth = maxPeakWidth>.5 ? maxPeakWidth : 10
		minSep = minSep>minPeakWidth ? minSep : 40
		threshAboveAvg = threshAboveAvg>0 ? threshAboveAvg : 100
		Prompt imageName,"name of image with peaks to find",popup,reverseList(WaveListClass("speImageNoBkg;speImage;rawImageNoBkg;rawImage*","*","DIMS:2"))
		Prompt minPeakWidth,"minimum allowed width of a peak"
		Prompt maxPeakWidth,"maximum allowed width of a peak"
		Prompt minSep,"minimum distance between two peaks"
		Prompt threshAboveAvg,"min peak height"
		DoPrompt/Help="3D-Xray Diffraction[Fit Peaks]" "peak fitting",imageName,minPeakWidth,maxPeakWidth,minSep,threshAboveAvg
		if (V_flag)
			return ""
		endif
		Wave image = $imageName
		printf "FitPeaksStepWise(%s,%g,%g,%g,%g)\r",imageName,minPeakWidth,maxPeakWidth,minSep,threshAboveAvg
		printIt = 1
	endif
	if (!WaveExists(image) || WaveDims(image)!=2 || !(minPeakWidth>0) || !(maxPeakWidth>0) || !(minSep>3) || !(threshAboveAvg>0))
		return ""
	endif
	imageName = GetWavesDataFolder(image,2)
	Variable Nx=DimSize(image,0), Ny=DimSize(image,1)
	if (Nx<2 || Ny<2)
		return ""
	endif

	Variable timer0=startMSTimer

	Duplicate/O image FitPeaksStepWise_smth					// image to use when finding next hot pixel
	Wave imageS = FitPeaksStepWise_smth
	Redimension/D imageS
	imageS = image[p][q]<=0 ? NaN : image[p][q]
	MatrixFilter/N=3 min  imageS

	ImageStats/M=1/Q imageS									// find the average value
	Variable threshold = threshAboveAvg+max(V_avg,0)		// the actual threshold, only consider pixels above this level
	if (printIt)
		printf "average for the image = %g,  so only consider pixels values above %g\r",V_avg,threshold
	endif

	Make/N=(Nx,Ny)/O/B/U FitPeaks_ImageMask				// make a mask for this image
	FitPeaks_ImageMask = 0									// set to 0, this includes all of the image
	if (!ParamIsDefault(mask))								// preset FitPeaks_ImageMask with 'mask' if it was passed and is valid
		if (WaveExists(mask) && DimSize(mask,0)==Nx && DimSize(mask,1)==Ny)	// if mask exists and is right size, use it as starting point
			FitPeaks_ImageMask = !(!mask[p][q])
		endif
	endif

	Duplicate/O image image_weights_temp_					// weights to use in the peak fitting
	Wave weights = image_weights_temp_
	weights = FitPeaks_ImageMask<1

	Make/N=(3,3)/O/U/I FitPeaks_3x3
	Make/O/T JZT_Constraints={"K3 > 0","K5 > 0","K6 > 0","K6 < 1"}

	Variable Nu=0, Nlen=50
	String wnote = note(image)
	String peakListName="FullPeakList"+detectorID2color(StringByKey("detectorID",wnote,"="))
	Make/N=(Nlen,12)/O $peakListName/WAVE=FullPeakList = NaN
	wnote = ReplaceNumberByKey("minPeakWidth",wnote,minPeakWidth,"=")
	wnote = ReplaceNumberByKey("maxPeakWidth",wnote,maxPeakWidth,"=")
	wnote = ReplaceNumberByKey("minSpotSeparation",wnote,minSep,"=")
	wnote = ReplaceNumberByKey("threshAboveAvg",wnote,threshAboveAvg,"=")
	wnote = ReplaceStringByKey("fittedIgorImage",wnote,GetWavesDataFolder(image,2),"=")
	wnote = ReplaceStringByKey("waveClass",wnote,"FittedPeakList","=")
	SetDimLabel 1,0,x0,FullPeakList			;	SetDimLabel 1,1,y0,FullPeakList
	SetDimLabel 1,2,x0Err,FullPeakList		;	SetDimLabel 1,3,y0Err,FullPeakList
	SetDimLabel 1,4,fwx,FullPeakList		;	SetDimLabel 1,5,fwy,FullPeakList
	SetDimLabel 1,6,fwxErr,FullPeakList	;	SetDimLabel 1,7,fwyErr,FullPeakList
	SetDimLabel 1,8,correlation,FullPeakList ;	SetDimLabel 1,9,correlationErr,FullPeakList
	SetDimLabel 1,10,area,FullPeakList		;	SetDimLabel 1,10,area,FullPeakList

	Variable fw= ( 2*sqrt(2*ln(2)) )							// this factor corrects for 2-d sigma to FWHM, apply to K3, K5, and the corresponding errors
	Variable xx,yy,fwx,fwy,sigX,sigY,amp,sigAmp				// holds the resutls of gaussian fit
	Variable px,py, left,right,top,bot
	do
		// DoUpdate
		ImageStats/M=1/R=FitPeaks_ImageMask imageS		// find the highest point
		if (V_max<threshold || V_npnts<2)
			break
		endif
		px=V_maxRowLoc ; py=V_maxColLoc					// look for a peak at (px,py)
		//	if (px==491 && py==157)
		//		printf "point = [%d, %d]\r",px,py
		//	endif

		// check if this is just a lone pixel
		FitPeaks_3x3 = image[p+px-1][q+py-1]
		FitPeaks_3x3[1][1] = 0
		ImageStats/M=1 FitPeaks_3x3
		if (!(V_max>threshold))									// this is a lone pixel
			image[px][py] = V_avg								// set lone pixels to surrounding average, and continue
			FitPeaks_ImageMask[px][py] = 255
			weights[px][py] = 0
			continue
		endif

		// fit a gaussian about px,py
		left = limit(floor(px-maxPeakWidth/2),0,Nx-1)		// box to use for fitting
		right = limit(ceil(px+maxPeakWidth/2),0,Nx-1)
		top = limit(floor(py-maxPeakWidth/2),0,Ny-1)
		bot = limit(ceil(py+maxPeakWidth/2),0,Ny-1)
		if (FitGaussianPeak(image,left,right,top,bot,weight=weights))			// returns 0 == OK,  non-zero if error
			FitPeaks_ImageMask[left,right][top,bot] = 255	// extend the FitPeaks_ImageMask to reject spots within maxPeakWidth
			weights[left,right][top,bot] = 0
			continue
		endif
		Wave W_coef=W_coef, W_sigma=W_sigma
		xx = W_coef[2]
		yy = W_coef[4]
		fwx = abs(W_coef[3] * fw)
		fwy = abs(W_coef[5] * fw)
		sigX = W_sigma[2]
		sigY = W_sigma[4]
		amp = abs(W_coef[1])
		sigAmp = W_sigma[2]
		//	if (ItemsInList(GetRTStackInfo(0))<2)
		//		printf "%d    Gaussian peak in at (%.2f±%.2f, %.2f±%.2f) with FWHM of (%.2f±%.2f, %.2f±%.2f),   amp=%g\r",Nu,xx,sigX,yy,sigY,fwx,W_sigma[3]*fw,fwy,W_sigma[5]*fw,W_coef[1]*fwx*fwy
		//	endif

		if (min(fwx,fwy)<minPeakWidth || max(fwx,fwy)>maxPeakWidth || sigX>maxPeakWidth || sigY>maxPeakWidth || amp<threshAboveAvg || amp/sigAmp<3)	// cannot have spots too wide or too narrow
//			left = limit(floor(px-min(minSep,fwx)),0,Nx-1)// fitted peak is invalid, delete area, use a smaller region than minSep
//			right = limit(ceil(px+min(minSep,fwx)),0,Nx-1)
//			top = limit(floor(py-min(minSep,fwy)),0,Ny-1)
//			bot = limit(ceil(py+min(minSep,fwy)),0,Ny-1)
			left = limit(floor(px-min(maxPeakWidth/2,fwx)),0,Nx-1)// fitted peak is invalid, delete area
			right = limit(ceil(px+min(maxPeakWidth/2,fwx)),0,Nx-1)
			top = limit(floor(py-min(maxPeakWidth/2,fwy)),0,Ny-1)
			bot = limit(ceil(py+min(maxPeakWidth/2,fwy)),0,Ny-1)
			FitPeaks_ImageMask[left,right][top,bot] = 255	// extend the FitPeaks_ImageMask to reject spots within maxPeakWidth
			weights[left,right][top,bot] = 0
			continue
		endif

		// exclusion area around (px,py) for future peaks when this peak is good
//px = round(xx)
//py = round(yy)
		left = limit(px-minSep,0,Nx-1)						// fitted peak was valid, delete area so it does not get counted twice
		right = limit(px+minSep,0,Nx-1)
		top = limit(py-minSep,0,Ny-1)
		bot = limit(py+minSep,0,Ny-1)
		FitPeaks_ImageMask[left,right][top,bot] = 255		// extend the FitPeaks_ImageMask to reject spots within maxPeakWidth
		weights[left,right][top,bot] = 0

		if (Nu>=Nlen)												// FullPeakList is too short, extend it
			Nlen += 50
			Redimension/N=(Nlen,-1) FullPeakList
		endif
		FullPeakList[Nu][0]=xx					;	FullPeakList[Nu][1]=yy
		FullPeakList[Nu][2]=sigX				;	FullPeakList[Nu][3]=sigY
		FullPeakList[Nu][4]=fwx				;	FullPeakList[Nu][5]=fwy
		FullPeakList[Nu][6]=W_sigma[3]*fw	;	FullPeakList[Nu][7]=W_sigma[5]*fw
		FullPeakList[Nu][8]=W_coef[6]		;	FullPeakList[Nu][9]=W_sigma[6]
		FullPeakList[Nu][10]=W_coef[1]*fwx*fwy
		FullPeakList[Nu][11]=W_coef[1]
		Nu += 1
	while(Nu<maxNum)
	Redimension/N=(Nu,-1) FullPeakList						// set to exact length
	ImageStats/M=1/G={0,Nu-1,10,10}/Q FullPeakList
	wnote = ReplaceNumberByKey("totalPeakIntensity",wnote,V_avg*V_npnts,"=")	// total intensity in all fitted peaks
	wnote = ReplaceNumberByKey("totalIntensity",wnote,sum(image),"=")			// total intensity in the image
	Note/K FullPeakList,wnote

//	String/G pkLIst=""
//	String str
//	Variable i
//	for (i=0;i<Nu;i+=1)											// re-set pkLIst to the fitted peaks
//		sprintf str,"%.0f,%.0f;",FullPeakList[i][0],FullPeakList[i][1]
//		pkList += str
//	endfor

	Variable sec0=stopMSTimer(timer0)/1e6
	if (printIt)
		printf "found and fitted %d peaks,  this all took %g sec\r",Nu,sec0
		printf "  result stored in the wave '%s'\r", GetWavesDataFolder(FullPeakList,2)
	endif
	KillWaves/Z FitPeaks_ImageMask, FitPeaks_3x3, image_weights_temp_
	KillWaves/Z FitPeaksStepWise_smth, JZT_Constraints
	return GetWavesDataFolder(FullPeakList,2)
End
//
// this function does not appear to be used anywhere.
//		Function/S FitPeaksProto(image,minPeakWidth,maxPeakWidth,minSep,threshAboveAvg,[mask,maxNum,whoami,printIt])
//			Wave image						// the image from the image file
//			Variable minPeakWidth		// minimum width for an acceptable peak (pixels)
//			Variable maxPeakWidth		// maximum width for an acceptable peak (pixels)
//			Variable minSep				// minimum separation between two peaks (pixels)
//			Variable threshAboveAvg	// threshold above average value,  use 25
//			Wave mask						// starting mask for the fit (this mask is unchanged by this routine)
//			Variable maxNum				// maximum number of peaks to find, normally goes to completion
//			Variable whoami				// if passed, then just return the name of this routine
//			Variable printIt
//			return ""
//		End
// 
Static Function FitGaussianPeak(image,left,right,bottom,top,[weight])	// returns 1 if error,  0 is OK
	// the result is passed to calling function by W_coef or {K0, K1, ... K6}
	// this fits a 2-d gaussina to a region on an image.  The region is selected with a marquee in an image plot, or specified
	Wave image
	Variable left, right, bottom, top								// range to use for fit (defaults to full image for bad values)
	Wave weight
	Variable noWeight = ParamIsDefault(weight)

	left = left==limit(left,0,DimSize(image,0)-1) ? left : 0	// for bad values, use whole image
	bottom = bottom==limit(bottom,0,DimSize(image,1)-1) ? bottom : 0
	right = right==limit(right,0,DimSize(image,0)-1) ? right : DimSize(image,0)-1
	top = top==limit(top,0,DimSize(image,1)-1) ? top : DimSize(image,1)-1

	// set up for the gaussian fit
	Wave/T JZT_Constraints=JZT_Constraints					// = {"K3 > 0","K5 > 0","K6 > 0","K6 < 1"}
	Variable V_FitOptions=4, V_FitError=0,V_FitQuitReason=0
	V_FitOptions=4

	if (!noWeight)												// limit box to non-zero weights
		Variable l,r,t,b
		maskBoundingBox(weight,l,r,t,b)
		left = max(left,l)
		right = min(right,r)
		top = min(top,b)
		bottom = max(bottom,t)
	endif
	if ((right-left)<3 || (top-bottom)<3)
		return 1													// need bigger region for a 2d fit
	endif

	if (noWeight)
		CurveFit/Q/N Gauss2D image[left,right][bottom,top]/C=JZT_Constraints
	else
		CurveFit/Q/N Gauss2D image[left,right][bottom,top]/C=JZT_Constraints/W=weight
	endif

	if (!V_FitError)												// if no error, but constraints failed, set bit 5 as error anyhow
		Wave W_coef = W_coef
		V_FitError += (W_coef[3]<-1e-7 || W_coef[5]<-1e-7 || W_coef[6]<-1e-7 || W_coef[6]>1.0001) ? 32 : 0
	endif
	return V_FitError
End
//
Static Function maskBoundingBox(mask,left,right,top,bot)
	Variable &left,&right,&top,&bot
	Wave mask
	left = Inf
	right = -Inf
	top = Inf
	bot = -Inf

	Variable i,Nx=DimSize(mask,0), Ny=DimSize(mask,1)
	for (i=0;i<Nx;i+=1)
		ImageStats/M=1/G={i,i,0,Ny-1} mask
		if (V_max)
			left = min(left,V_maxRowLoc)
			right = max(right,V_maxRowLoc)
		endif
	endfor
	for (i=0;i<Ny;i+=1)
		ImageStats/M=1/G={0,Nx-1,i,i} mask
		if (V_max)
			top = min(top,V_maxColLoc)
			bot = max(bot,V_maxColLoc)
		endif
	endfor
End


//Static Function FitGaussianPeak(image,left,right,bottom,top,[weight])	// returns 1 if error,  0 is OK
//	// the result is passed to calling function by W_coef or {K0, K1, ... K6}
//	// this fits a 2-d gaussina to a region on an image.  The region is selected with a marquee in an image plot, or specified
//	Wave image
//	Variable left, right, bottom, top								// range to use for fit
//	Wave weight
//	Variable noWeight = ParamIsDefault(weight)
//
//	// set up for the gaussian fit
//	Wave JZT_Constraints=JZT_Constraints					// = {"K3 > 0","K5 > 0","K6 > 0","K6 < 1"}
//	Variable V_FitOptions=4, V_FitError=0,V_FitQuitReason=0
//	V_FitOptions=4
//	if (noWeight)
//		CurveFit/Q/N Gauss2D image[left,right][bottom,top]/C=JZT_Constraints
//	else
//		CurveFit/Q/N Gauss2D image[left,right][bottom,top]/C=JZT_Constraints/W=weight
//	endif
//
//	if (!V_FitError)												// if no error, but constraints failed, set bit 5 as error anyhow
//		Wave W_coef = W_coef
//		V_FitError += (W_coef[3]<-1e-7 || W_coef[5]<-1e-7 || W_coef[6]<-1e-7 || W_coef[6]>1.0001) ? 32 : 0
//	endif
//	return V_FitError
//End

// identify the peaks, and fit them all (no background removal)
Function/S FitPeaksNew(image,threshAboveAvg,dist,[printIt])
	Wave image						// the image from the image file
	Variable threshAboveAvg		// threshold above average value
	Variable dist					// minimum distance between spots, spots have to be at least this far apart (pixels)
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt

	if (!WaveExists(image) || WaveDims(image)!=2 || !(threshAboveAvg>0) || !(dist>3))
		String imageName = SelectString(WaveExists(image),"",NameOfWave(image))
		threshAboveAvg = threshAboveAvg>0 ? threshAboveAvg : 100
		dist = dist>3 ? dist : 30
		Prompt imageName,"name of image with peaks to find",popup,reverseList(WaveListClass("speImageNoBkg;speImage;rawImageNoBkg;rawImage*","*","DIMS:2"))
		if (NumVarOrDefault(MICRO_GEOMETRY_VERSION_PATH,0)&2)
			Prompt threshAboveAvg,"Peak Threshold"
			Prompt dist, "minimum distance (pixels) between any two peaks (typically 5)"
		else
			Prompt threshAboveAvg,"fraction of image that is background, usually in range [0.7, 0.995], use 0 to skip bkg removal"
			Prompt dist, "minimum distance (pixels) between any two peaks (use ~30 for unbinned Si)"
		endif
		DoPrompt/Help="3D-Xray Diffraction[Fit Peaks]" "peak fitting",imageName,threshAboveAvg,dist
		if (V_flag)
			return ""
		endif
		Wave image = $imageName
		printf "FitPeaks(%s,%g,%g)\r",imageName,threshAboveAvg,dist
		printIt = 1
	endif
	if (!WaveExists(image) || WaveDims(image)!=2 || !(threshAboveAvg>0) || !(dist>3))
		return ""
	endif
	Variable Nx=DimSize(image,0), Ny=DimSize(image,1)

	Variable timer0=startMSTimer
	ImageStats/M=1/Q image									// find the average value
	Variable threshold = threshAboveAvg+V_avg					// the actual threshold, only consider pixels above this level
if (printIt)
print "average for the image = ",V_avg
endif
	Duplicate/O image, FitPeaks_ImageMask
	Wave FitPeaks_ImageMask=FitPeaks_ImageMask
	Redimension/B/U FitPeaks_ImageMask					// make a mask for the image
	FitPeaks_ImageMask = 0									// set to 0, this includes all of the image
	Variable localAvg											// average of local spots (not including center)

	Make/N=(3,3)/O/U/I FitPeaks_3x3
	Variable xx,yy,fwx,fwy,sigX,sigY
	Variable fw= ( 2*sqrt(2*ln(2)) )							// this factor corrects for 2-d sigma to FWHM, apply to K3, K5, and the corresponding errors

	String wnote = note(image)
	String peakListName="FullPeakList"+detectorID2color(StringByKey("detectorID",wnote,"="))
	Make/N=(50,12)/O $peakListName/WAVE=FullPeakList = NaN
	wnote = ReplaceNumberByKey("minSpotSeparation",wnote,dist,"=")
	wnote = ReplaceNumberByKey("threshAboveAvg",wnote,threshAboveAvg,"=")
	wnote = ReplaceStringByKey("fittedIgorImage",wnote,GetWavesDataFolder(image,2),"=")
	wnote = ReplaceStringByKey("waveClass",wnote,"FittedPeakList","=")
	Note/K FullPeakList,wnote
	SetDimLabel 1,0,x0,FullPeakList			;	SetDimLabel 1,1,y0,FullPeakList
	SetDimLabel 1,2,x0Err,FullPeakList		;	SetDimLabel 1,3,y0Err,FullPeakList
	SetDimLabel 1,4,fwx,FullPeakList		;	SetDimLabel 1,5,fwy,FullPeakList
	SetDimLabel 1,6,fwxErr,FullPeakList	;	SetDimLabel 1,7,fwyErr,FullPeakList
	SetDimLabel 1,8,correlation,FullPeakList ;	SetDimLabel 1,9,correlationErr,FullPeakList
	SetDimLabel 1,10,area,FullPeakList		;	SetDimLabel 1,11,amp,FullPeakList

	Variable i
	Variable px,py, left,right,top,bot
	Variable Nu=0, Nlen=50
	do
		ImageStats/M=1/R=FitPeaks_ImageMask image			// find the highest point
		if (V_max<threshold || V_npnts<2)
			break
		endif
		px = V_maxRowLoc
		py = V_maxColLoc

		// check if this is just a lone pixel
		FitPeaks_3x3 = image[p+px][q+py]
		FitPeaks_3x3[1][1] = 0
		ImageStats/M=1 FitPeaks_3x3
		localAvg = V_avg
		if (!(V_max>threshold))								// this is a lone pixel
			image[px][py] = localAvg							// set lone pixels to local average, and continue
			FitPeaks_ImageMask[px][py] = 255
			continue
		endif

		// fit a gaussian about px,py
		left = limit(px-dist/2,0,Nx-1)
		right = limit(px+dist/2,0,Nx-1)
		top = limit(py-dist/2,0,Ny-1)
		bot = limit(py+dist/2,0,Ny-1)
		if (FitOneGaussianPeak(image,left,right,top,bot))		// returns 0 == OK,  non-zero if error
			FitPeaks_ImageMask[left,right][top,bot] = 255		// extend the FitPeaks_ImageMask to include this spot
			continue
		endif
		Wave W_coef=W_coef, W_sigma=W_sigma
		xx = W_coef[2]
		yy = W_coef[4]
		fwx = W_coef[3] * fw
		fwy = W_coef[5] * fw
		sigX = W_sigma[2]
		sigY = W_sigma[4]
		if (fwx>dist || fwy>dist || fwx<1 || fwy<1 || sigX>dist || sigY>dist)	// cannot have spots wider than dist
			FitPeaks_ImageMask[left,right][top,bot] = 255		// extend the FitPeaks_ImageMask to include this spot
			continue
		endif
		// printf "%d  %d  Gaussian peak in at (%.2f±%.2f, %.2f±%.2f) with FWHM of (%.2f±%.2f, %.2f±%.2f),   amp=%g\r",i,Nu,xx,sigX,yy,sigY,fwx,W_sigma[3]*fw,fwy,W_sigma[5]*fw,W_coef[1]*fwx*fwy
		// check here that the found xx,yy is not too close to existing peaks (if if is stamp out with roi again)
		for (i=0;i<Nu;i+=1)									// check that this is not too close to an existing peak
			if (sqrt((FullPeakList[i][0]-xx)^2 + (FullPeakList[i][1]-yy)^2)<dist)
				break
			endif
		endfor
		if (i<Nu)												// there is already a spot too close to here
			FitPeaks_ImageMask[left,right][top,bot] = 255		// extend the FitPeaks_ImageMask to include this spot
			continue
		endif
		if (Nu>=Nlen)											// FullPeakList is too short, extend it
			Nlen += 50
			Redimension/N=(Nlen,-1) FullPeakList
		endif
if (printIt)
printf "found Gaussian start @ (%d, %d)=%d   with left=%.1f,  right=%.1f,  top=%.1f,  bot=%.1f      result=(%.3f, %.3f)\r",px,py,image[px][py],left,right,top,bot,xx,yy
endif
		FullPeakList[Nu][0]=xx					;	FullPeakList[Nu][1]=yy
		FullPeakList[Nu][2]=sigX				;	FullPeakList[Nu][3]=sigY
		FullPeakList[Nu][4]=fwx				;	FullPeakList[Nu][5]=fwy
		FullPeakList[Nu][6]=W_sigma[3]*fw	;	FullPeakList[Nu][7]=W_sigma[5]*fw
		FullPeakList[Nu][8]=W_coef[6]		;	FullPeakList[Nu][9]=W_sigma[6]
		FullPeakList[Nu][10]=W_coef[1]*fwx*fwy
		FullPeakList[Nu][11]=W_coef[1]
		Nu += 1
		FitPeaks_ImageMask[left,right][top,bot] = 255			// extend the FitPeaks_ImageMask to include this spot
	while(1)
	Redimension/N=(Nu,-1) FullPeakList

//	String/G pkLIst=""
//	String str
//	for (i=0;i<Nu;i+=1)
//		sprintf str,"%.0f,%.0f;",FullPeakList[i][0],FullPeakList[i][1]
//		pkList += str
//	endfor

	Variable sec0=stopMSTimer(timer0)/1e6
	if (printIt)
		printf "found and fitted %d peaks,  this all took %g sec\r",i,sec0
		printf "  result stored in the wave '%s'\r", GetWavesDataFolder(FullPeakList,2)
	endif
	KillWaves/Z FitPeaks_ImageMask, FitPeaks_3x3
	return GetWavesDataFolder(FullPeakList,2)
End


// remove bkg from an image, identify the peaks, and fit them all
Function/S FitPeaks(image,fractionBkg,dist,[printIt])
	Wave image					// the image from the image file (with or without bkg removed)
	Variable fractionBkg			// fraction of image that is background, use something like 0.99 or 0.995 (even 0.7 works OK)
	Variable dist					// minimum distance between spots, spots have to be at least this far apart (pixels)
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt

	if (!WaveExists(image) || WaveDims(image)!=2 || !(fractionBkg==limit(fractionBkg,1e-3,1) || !(dist>0)))
		String imageName = SelectString(WaveExists(image),"",NameOfWave(image))
		fractionBkg = fractionBkg==limit(fractionBkg,1e-3,1) ? fractionBkg : 0.9
		dist = dist>0 ? dist : 30
		Prompt imageName,"name of image with peaks to find",popup,reverseList(WaveListClass("speImage*;rawImage*","*","DIMS:2"))
		Prompt fractionBkg,"fraction of image that is background, usually in range [0.7, 0.995]"
		Prompt dist, "minimum distance (pixels) between any two peaks (use ~30 for unbinned Si)"
		DoPrompt/Help="3D-Xray Diffraction[Fit Peaks]" "peak fitting",imageName,fractionBkg,dist
		if (V_flag)
			return ""
		endif
		Wave image = $imageName
		printf "FitPeaks(%s,%g,%g)\r",imageName,fractionBkg,dist
		printIt = 1
	endif
	if (!WaveExists(image) || WaveDims(image)!=2 || !(fractionBkg==limit(fractionBkg,1e-3,1) || !(dist>0)))
		return ""
	endif
	Variable timer=startMSTimer

//	timer2 = startMSTimer
	Make/N=100/O hist_
	SetScale/P x 0,2,"", hist_
	Histogram/B=2 image,hist_
//	print "\t\tHistogram part took ",stopMSTimer(timer2)/1e6,"sec"

	Variable i = BinarySearchInterp(hist_, hist_(0)/10)
	i = numtype(i) ? numpnts(hist_)/2 : i
 	Variable cutOff = pnt2x(hist_,i)
// 	Variable cutOff = pnt2x(hist_,BinarySearchInterp(hist_, hist_(0)/10))
	printf "search for peaks above a threshold of %g\r",cutOff
//	print "cutOff=",cutOff											// threshold in image wave for finding peaks
//	timer2 = startMSTimer
	Duplicate/O image,FitPeaks_mask_bigger

	ImageThreshold/O/Q/I/T=(cutOff)/M=0 FitPeaks_mask_bigger
//	print "\t\tmaking mask part took ",stopMSTimer(timer2)/1e6,"sec"

	DoUpdate
//	timer2 = startMSTimer
	Variable id = dist>20 ? 5 : 1
//	Imagemorphology/O/E=5 BinaryDilation FitPeaks_mask_bigger
	Imagemorphology/O/E=(id) BinaryDilation FitPeaks_mask_bigger
//	print "\t\tImagemorphology part took ",stopMSTimer(timer2)/1e6,"sec"

//	timer2 = startMSTimer
	ImageAnalyzeParticles/A=5/M=3/Q stats FitPeaks_mask_bigger
	KillWaves/Z M_RawMoments,M_Particle
	Sort/R W_ImageObjArea,W_ImageObjArea,W_SpotX,W_SpotY,W_circularity,W_rectangularity,W_ImageObjPerimeter,W_xmin,W_xmax,W_ymin,W_ymax
//	print "\t\tImageAnalyzeParticles part took ",stopMSTimer(timer2)/1e6,"sec"

//	timer2 = startMSTimer
	i = reFit_GaussianPkList(image,dist)
//	print "\t\t refitting the peaks took ",stopMSTimer(timer2)/1e6,"sec"

	Variable sec0=stopMSTimer(timer)/1e6
	if (printIt)
		printf "fitted %d peaks,  this all took %g sec\r",i,sec0
		printf "  result stored in the wave '%s'\r", GetWavesDataFolder(image,2)
	endif
//	Wave FullPeakList=FullPeakList

	KillWaves/Z hist_, W_sigma,W_coef,W_ParamConfidenceInterval, FitPeaks_mask_bigger
	KillWaves/Z W_xmin,W_xmax,W_ymin,W_ymax,W_circularity,W_rectangularity,W_ImageObjPerimeter,W_ImageObjArea,W_SpotX,W_SpotY
	return GetWavesDataFolder(image,2)
End



//		// remove bkg from an image, identify the peaks, and fit them all
//		Function/S OLD_removeBkg_FitPeaks(rawImage,fractionBkg,dist)
//			Wave rawImage					// the raw image from the image file
//			Variable fractionBkg			// fraction of image that is background, use something like 0.99 or 0.995 (even 0.7 works OK)
//			Variable dist					// minimum distance between spots, spots have to be at least this far apart (pixels)
//		
//			Variable printIt = ItemsInList(GetRTStackInfo(0))<2
//			if (!WaveExists(rawImage) || WaveDims(rawImage)!=2 || !(fractionBkg==limit(fractionBkg,1e-3,1) || !(dist>0)))
//				String imageName = SelectString(WaveExists(rawImage),"",NameOfWave(rawImage))
//				fractionBkg = fractionBkg==limit(fractionBkg,1e-3,1) ? fractionBkg : 0.8
//				dist = dist>0 ? dist : 30
//		//		Prompt imageName,"name of image with peaks to find",popup,WaveList("*",";","DIMS:2")
//		//		Prompt imageName,"name of image with peaks to find",popup,WaveList_Tags("imageFileName","rawIgorImage","*","DIMS:2")
//				Prompt imageName,"name of image with peaks to find",popup,reverseList(WaveListClass("speImage;rawImage*","*","DIMS:2"))
//				Prompt fractionBkg,"fraction of image that is background, usually in range [0.7, 0.995]"
//				Prompt dist, "minimum distance (pixels) between any two peaks (use ~30 for unbinned Si)"
//				DoPrompt/Help="3D-Xray Diffraction[Fit Peaks]" "peak fitting",imageName,fractionBkg,dist
//				if (V_flag)
//					return ""
//				endif
//				Wave rawImage = $imageName
//				printf "FitPeaks(%s,%g,%g)\r",imageName,fractionBkg,dist
//				printIt = 1
//			endif
//			if (!WaveExists(rawImage) || WaveDims(rawImage)!=2 || !(fractionBkg==limit(fractionBkg,1e-3,1) || !(dist>0)))
//				return ""
//			endif
//		//	Variable timer=startMSTimer
//			Variable timer0=startMSTimer
//			Variable thresh = numpnts(rawImage)*fractionBkg
//		
//			WaveStats/Q/M=1 rawImage
//			Make/N=200/D/O FitPeaks_Hist
//			Variable i
//			for (i=-1;i<10 || numtype(i);)							// keep expanding the histogram range while i<10,   I want i to be bigger than 10
//				SetScale/I x,V_min,V_max,"",FitPeaks_Hist			// lower high end of histogram to zoom in around zero
//				Histogram/B=2 rawImage,FitPeaks_Hist
//				Integrate/P FitPeaks_Hist/D=FitPeaks_Hist_INT	// integrate the histogram
//				i=BinarySearchInterp(FitPeaks_Hist_INT,thresh)	// find point where FitPeaks_Hist_INT[i] == thresh
//				V_max = numtype(i) ? V_max/2 : pnt2x(FitPeaks_Hist_INT,i+4)
//			endfor
//		//	print "first part took ",stopMSTimer(timer)/1e6,"sec"
//		//	timer = startMSTimer
//		
//			thresh = DimDelta(FitPeaks_Hist,0)*i + DimOffset(FitPeaks_Hist,0)
//			ImageThreshold/Q/T=(thresh)/M=0 rawImage					// create M_ImageThresh
//			Wave M_ImageThresh=M_ImageThresh
//		
//			ImageStats/M=1 M_ImageThresh
//			Variable fracBlack = V_avg*V_npnts/255/V_npnts
//			//	printf "number of points in peak is %d out of %d points (%.2f%%)\r",fracBlack*V_npnts,V_npnts,fracBlack*100
//			Make/N=(9,9)/O/B/U black
//			black = 255
//			for(;fracBlack<0.10;)											// loop until fracBlack is atleast 10%
//				Imagemorphology/S=black/O BinaryDilation M_ImageThresh
//				ImageThreshold/Q/T=(2)/M=0 M_ImageThresh
//				ImageStats/M=1 M_ImageThresh
//				fracBlack = V_avg*V_npnts/255/V_npnts
//				//	printf "number of points in peak is %d out of %d points (%.2f%%)\r",fracBlack*V_npnts,V_npnts,fracBlack*100
//			endfor
//			ImageTransform/O invert M_ImageThresh						//	M_ImageThresh = !M_ImageThresh
//		
//		//	print "second part took ",stopMSTimer(timer)/1e6,"sec"
//		//	timer = startMSTimer
//		
//		//	Variable timer2 = startMSTimer
//			ImageInterpolate /f={.125,.125} bilinear rawImage
//			Redimension/S M_InterpolatedImage
//		
//			KillWaves/Z FitPeaks_rawImage_smaller
//			Rename M_InterpolatedImage FitPeaks_rawImage_smaller		// make a source image using fewer pixels
//		//	print "\t\t1st interpolate part took ",stopMSTimer(timer2)/1e6,"sec"
//		
//		//	timer2 = startMSTimer
//		 	ImageInterpolate /f={.125,.125} bilinear M_ImageThresh		// make a mask with fewer pixels
//			KillWaves/Z FitPeaks_mask_smaller
//			Rename M_InterpolatedImage FitPeaks_mask_smaller
//			Redimension/B/U FitPeaks_mask_smaller
//			ImageRemoveBackground/F/R=FitPeaks_mask_smaller/P=2/W FitPeaks_rawImage_smaller	// creates M_RemovedBackground, uses 3rd order polynomial
//			Wave W_BackgroundCoeff=W_BackgroundCoeff
//			//	print W_BackgroundCoeff
//		//	print "\t\t2nd interpolate part took ",stopMSTimer(timer2)/1e6,"sec"
//		
//		//	timer2 = startMSTimer
//			ImageInterpolate /f={8,8} bilinear M_RemovedBackground		// creates M_InterpolatedImage, a bkg with full number of pixels
//			Wave M_InterpolatedImage=M_InterpolatedImage
//		
//			String outName = NameOfWave(rawImage)+"NoBkg"
//			String fldr = GetWavesDataFolder(rawImage,1)
//			if (CheckName(fldr+outName,1))
//				outName = CleanupName(outName,0)
//			endif
//			outName = fldr+outName
//			Duplicate/O rawImage $outName
//			Wave noBkg = $outName
//		
//			String wnote=ReplaceStringByKey("rawIgorImage",note(rawImage),GetWavesDataFolder(rawImage,2),"=")
//			wnote = ReplaceNumberByKey("fractionBkg",wnote,fractionBkg,"=")
//			wnote = ReplaceNumberByKey("minSpotSeparation",wnote,dist,"=")
//			wnote = ReplaceStringByKey("waveClass",wnote,"rawImageNoBkg","=")
//			Note/K noBkg, wnote
//			Redimension/I noBkg											// signed int32 version of rawImage
//			noBkg = rawImage - M_InterpolatedImage[p][q]					// p,q because M_InterpolatedImage is smaller
//		//	print "\t\t",DimSize(rawImage,0),DimSize(rawImage,1)
//		//	print "\t\t",DimSize(noBkg,0),DimSize(noBkg,1)
//		//	print "\t\t",DimSize(M_InterpolatedImage,0),DimSize(M_InterpolatedImage,1)
//		//	print "\t\t3rd interpolate part took ",stopMSTimer(timer2)/1e6,"sec"
//		
//		
//		//	timer2 = startMSTimer
//			Make/N=100/O hist_
//			SetScale/P x 0,2,"", hist_
//			Histogram/B=2 noBkg,hist_
//		//	print "\t\tHistogram part took ",stopMSTimer(timer2)/1e6,"sec"
//		
//			i = BinarySearchInterp(hist_, hist_(0)/10)
//			i = numtype(i) ? numpnts(hist_)/2 : i
//		 	Variable cutOff = pnt2x(hist_,i)
//		// 	Variable cutOff = pnt2x(hist_,BinarySearchInterp(hist_, hist_(0)/10))
//			printf "after background removal search for peaks above a threshold of %g\r",cutOff
//		//	print "cutOff=",cutOff											// threshold in NoBkg wave for finding peaks
//		//	timer2 = startMSTimer
//			Duplicate/O noBkg,FitPeaks_mask_bigger
//		
//			ImageThreshold/O/Q/I/T=(cutOff)/M=0 FitPeaks_mask_bigger
//		//	print "\t\tmaking mask part took ",stopMSTimer(timer2)/1e6,"sec"
//		
//			DoUpdate
//		//	timer2 = startMSTimer
//			Variable id = dist>20 ? 5 : 1
//		//	Imagemorphology/O/E=5 BinaryDilation FitPeaks_mask_bigger
//			Imagemorphology/O/E=(id) BinaryDilation FitPeaks_mask_bigger
//		//	print "\t\tImagemorphology part took ",stopMSTimer(timer2)/1e6,"sec"
//		
//		//	timer2 = startMSTimer
//			ImageAnalyzeParticles/A=5/M=3/Q stats FitPeaks_mask_bigger
//			KillWaves/Z M_RawMoments,M_Particle
//			Sort/R W_ImageObjArea,W_ImageObjArea,W_SpotX,W_SpotY,W_circularity,W_rectangularity,W_ImageObjPerimeter,W_xmin,W_xmax,W_ymin,W_ymax
//		//	print "\t\tImageAnalyzeParticles part took ",stopMSTimer(timer2)/1e6,"sec"
//		
//		//	timer2 = startMSTimer
//			i = reFit_GaussianPkList(noBkg,dist)
//		
//		//	print "\t\t refitting the peaks took ",stopMSTimer(timer2)/1e6,"sec"
//		
//		//	print "fancy part took ",stopMSTimer(timer)/1e6,"sec"	
//			Variable sec0=stopMSTimer(timer0)/1e6
//			if (printIt)
//				printf "found and fitted %d peaks,  this all took %g sec\r",i,sec0
//		//		print "peak finding and fitting all took ",sec0,"sec"	
//				printf "  result stored in the wave '%s'\r", GetWavesDataFolder(noBkg,2)
//			endif
//		//	Wave FullPeakList=FullPeakList
//		
//			KillWaves/Z hist_, W_sigma,W_coef,W_ParamConfidenceInterval
//			KillWaves/Z FitPeaks_rawImage_smaller,FitPeaks_mask_smaller, M_ImageThresh, M_RemovedBackground, M_InterpolatedImage, W_BackgroundCoeff
//			KillWaves/Z FitPeaks_Hist_INT,FitPeaks_Hist, black, FitPeaks_mask_bigger
//			KillWaves/Z W_xmin,W_xmax,W_ymin,W_ymax,W_circularity,W_rectangularity,W_ImageObjPerimeter,W_ImageObjArea,W_SpotX,W_SpotY
//			return GetWavesDataFolder(noBkg,2)
//		End



// Fit every peak from the output of ImageAnalyzeParticles,  it sets the string pkList based on the fits
// it also sets FullPeakList[][10] which contains the exact info about the fit of each peak
Static Function reFit_GaussianPkList(image,dist)
	Wave image
	Variable dist						// minimum distance between spots (pixels)

	Wave W_spotX=W_spotX, W_spotY=W_spotY
	Wave W_xmin=W_xmin, W_xmax=W_xmax, W_ymin=W_ymin, W_ymax=W_ymax

	Variable Nx=DimSize(image,0), Ny=DimSize(image,1)
	Variable j,i,N=numpnts(W_spotX),err
	Variable fw= ( 2*sqrt(2*ln(2)) )			// this factor corrects for 2-d sigma to FWHM, apply to K3, K5, and the corresponding errors
	Variable fwx,fwy,sigX,sigY
	Variable xx,yy, tooClose, Nu=0
	Variable left,right,top,bot

	String peakListName="FullPeakList"+detectorID2color(StringByKey("detectorID",note(image),"="))
	Make/N=(N,12)/O $peakListName/WAVE=FullPeakList = NaN
	String wnote = ReplaceStringByKey("fittedIgorImage",note(image),GetWavesDataFolder(image,2),"=")
	wnote = ReplaceStringByKey("waveClass",wnote,"FittedPeakList","=")
	Note/K FullPeakList,wnote
	SetDimLabel 1,0,x0,FullPeakList			;	SetDimLabel 1,1,y0,FullPeakList
	SetDimLabel 1,2,x0Err,FullPeakList		;	SetDimLabel 1,3,y0Err,FullPeakList
	SetDimLabel 1,4,fwx,FullPeakList		;	SetDimLabel 1,5,fwy,FullPeakList
	SetDimLabel 1,6,fwxErr,FullPeakList	;	SetDimLabel 1,7,fwyErr,FullPeakList
	SetDimLabel 1,8,correlation,FullPeakList ;	SetDimLabel 1,9,correlationErr,FullPeakList
	SetDimLabel 1,10,area,FullPeakList		;	SetDimLabel 1,11,amp,FullPeakList
	for (i=0;i<N;i+=1)
		xx = (W_xmin[i]+W_xmax[i])/2
		yy = (W_ymin[i]+W_ymax[i])/2

		tooClose = 0
		for (j=0;j<i;j+=1)						// check only spots already fitted
			if ((xx-W_spotX[j])^2 + (yy-W_spotY[j])^2 <= dist^2)
				tooClose = 1
				break
			endif
		endfor
		if (tooClose)
			continue
		endif
		left = floor(W_xmin[i])
		right = ceil(W_xmax[i])
		top = floor(W_ymin[i])
		bot = ceil(W_ymax[i])
		if (right-left<dist/2)					// fit box cannot be too small
			left = xx - dist/2
			right = xx + dist/2
		endif
		if (bot-top<dist/2)
			top = yy - dist/2
			bot = yy + dist/2
		endif
		left = limit(left,0,Nx-1)
		right = limit(right,0,Nx-1)
		top = limit(top,0,Ny-1)
		bot = limit(bot,0,Ny-1)

		err = FitOneGaussianPeak(image,left,right,top,bot)	// returns 0 == OK,  non-zero if error
		if (err)
			continue
		endif
		Wave W_coef=W_coef, W_sigma=W_sigma
		xx = W_coef[2]
		yy = W_coef[4]
		fwx = W_coef[3] * fw
		fwy = W_coef[5] * fw
		sigX = W_sigma[2]
		sigY = W_sigma[4]
		if (fwx>dist || fwy>dist || fwx<1 || fwy<1 || sigX>dist || sigY>dist)	// cannot have spots wider than dist
			continue
		endif
		// printf "%d  %d  Gaussian peak in at (%.2f±%.2f, %.2f±%.2f) with FWHM of (%.2f±%.2f, %.2f±%.2f),   amp=%g\r",i,Nu,xx,sigX,yy,sigY,fwx,W_sigma[3]*fw,fwy,W_sigma[5]*fw,W_coef[1]*fwx*fwy

		tooClose = 0													// check that this spot is not too close to another already included
		for (j=0;j<Nu-1;j+=1)					// check all spots already fitted
			if ((FullPeakList[j][0]-xx)^2 + (FullPeakList[j][1]-yy)^2 <= dist^2)
				tooClose = 1
				break
			endif
		endfor
		if (tooClose)
			continue
		endif
		FullPeakList[Nu][0]=xx					;	FullPeakList[Nu][1]=yy
		FullPeakList[Nu][2]=sigX				;	FullPeakList[Nu][3]=sigY
		FullPeakList[Nu][4]=fwx				;	FullPeakList[Nu][5]=fwy
		FullPeakList[Nu][6]=W_sigma[3]*fw	;	FullPeakList[Nu][7]=W_sigma[5]*fw
		FullPeakList[Nu][8]=W_coef[6]		;	FullPeakList[Nu][9]=W_sigma[6]
		FullPeakList[Nu][10]=W_coef[1]*fwx*fwy
		FullPeakList[Nu][11]=W_coef[1]		// amplitude
//
//		sprintf str,"%.0f,%.0f;",xx,yy
//		sprintf str,"%d,%d;",round(xx),round(yy)
//		pkList += str
		Nu += 1
	endfor
	Redimension/N=(Nu,-1) FullPeakList
	KillWaves/Z W_sigma,W_coef,W_ParamConfidenceInterval
	return Nu
End


//	This was removed when purging the whole pkList thing
//
//// This routine allows the editing of which peaks to fit
//// Fit every peak in the list 'peaks',  it sets the string pkList based on the fits
//// it also sets FullPeakList[][11] which contains the exact info about the fit of each peak
//Function setFittedPeaksFromList(image,dist,peaks)
//	Wave image
//	Variable dist						// minimum distance between spots (pixels)
//	String peaks						// list of peaks, basically this is pkList
//
//	peaks = SelectString(strlen(peaks),StrVarOrDefault(":pkList",""),peaks)
//	if (strlen(peaks)<1)
//		DoAlert 0,"no list of peaks available in setFittedPeaksFromList(), cannot do anything"
//		return 0
//	endif
//	String str
//	if (!WaveExists(image) || !(dist>0))
//		String wList = reverseList(WaveListClass("speImage*;rawImage*","*","DIMS:2")), wName=""
//		Prompt dist,"minimum distance between peaks (pixels)"
//		Prompt wName,"Indexed Peak List",popup,wList
//		if (dist>0 && ItemsInList(wList)==1)	// only need image, and there is only one image
//			Wave image = $(StringFromList(0,wList))
//		else
//			if (dist>0)					// only need image, and need to ask
//				DoPrompt "image to use",wName
//				Wave image = $wName
//			elseif (WaveExists(image))		// only need dist
//				DoPrompt "dist between peakst",dist
//			else								// need both image and dist
//				DoPrompt "image to use",wName,dist
//				Wave image = $wName
//			endif
//			if (V_flag)
//				return 0
//			endif
//		endif
//		sprintf str, "setFittedPeaksFromList(%s,%g,\"%s\")",GetWavesDataFolder(image,2),dist,peaks
//		print str[0,390]
//	endif
//	if (!WaveExists(image) || !(dist>0))
//		return 0
//	endif
//
//	if (exists("pkLIst")!=2)
//		String/G pkLIst=""
//	endif
//	SVAR pkList = pkList
//	pkList = ""
//	String imageName = GetWavesDataFolder(image,2)
//	Variable Nx=DimSize(image,0), Ny=DimSize(image,1)
//	Variable j,i,N=ItemsInList(peaks),err,xx,yy
//
//	Make/N=(N)/O peakxx_,peakyy_,peakInt_		// re-order the points by peak intensity
//	for (i=0;i<N;i+=1)
//		str = StringFromList(i,peaks)
//		xx = str2num(StringFromList(0,str,","))
//		yy = str2num(StringFromList(1,str,","))
//		peakxx_[i] = xx							// fill arrays with poistions
//		peakyy_[i] = yy
//		peakInt_[i] = image[xx][yy]				// and intensities
//	endfor
//	Sort/R peakInt_, peakxx_,peakyy_,peakInt_	// sort by intensity
//	peaks=""
//	for (i=0;i<N;i+=1)								// and reset the peak list
//		sprintf str,"%.0f,%.0f;",peakxx_[i],peakyy_[i]
//		peaks += str
//	endfor
//	KillWaves/Z peakxx_,peakyy_,peakInt_
//
//	Variable fw= ( 2*sqrt(2*ln(2)) )			// this factor corrects for 2-d sigma to FWHM, apply to K3, K5, and the corresponding errors
//	Variable fwx,fwy,sigX,sigY
//	Variable Nu=0, tooClose
//	Variable left,right,top,bot
//	Make/N=(N,11)/O FullPeakList
//	String wnote = ReplaceStringByKey("fittedIgorImage",note(image),GetWavesDataFolder(image,2),"=")
//	wnote = ReplaceStringByKey("waveClass",wnote,"FittedPeakList","=")
//	Note/K FullPeakList,wnote
//	SetDimLabel 1,0,x0,FullPeakList		;	SetDimLabel 1,1,y0,FullPeakList
//	SetDimLabel 1,2,x0Err,FullPeakList	;	SetDimLabel 1,3,y0Err,FullPeakList
//	SetDimLabel 1,4,fwx,FullPeakList		;	SetDimLabel 1,5,fwy,FullPeakList
//	SetDimLabel 1,6,fwxErr,FullPeakList	;	SetDimLabel 1,7,fwyErr,FullPeakList
//	SetDimLabel 1,8,correlation,FullPeakList ;	SetDimLabel 1,9,correlationErr,FullPeakList
//	SetDimLabel 1,10,area,FullPeakList
//	FullPeakList = NaN
//	for (i=0;i<N;i+=1)
//		str = StringFromList(i,peaks)
//		xx = str2num(StringFromList(0,str,","))
//		yy = str2num(StringFromList(1,str,","))
//		left = limit(xx-dist/2,0,Nx-1)
//		right = limit(xx+dist/2,0,Nx-1)
//		top = limit(yy-dist/2,0,Ny-1)
//		bot = limit(yy+dist/2,0,Ny-1)
//		err = FitOneGaussianPeak(imageName,left,right,top,bot)	// returns 0 == OK,  non-zero if error
//		if (err)
//			continue
//		endif
//		Wave W_coef=W_coef, W_sigma=W_sigma
//		xx = W_coef[2]
//		yy = W_coef[4]
//		fwx = W_coef[3] * fw
//		fwy = W_coef[5] * fw
//		sigX = W_sigma[2]
//		sigY = W_sigma[4]
//		if (fwx>dist*2 || fwy>dist*2 || fwx<1 || fwy<1 || sigX>dist || sigY>dist)	// cannot have spots wider than 2*dist
//			continue
//		endif
//		for (j=0,tooClose=0;j<Nu;j+=1)							// check that this spot is not too close to an existing spot
//			if (abs(FullPeakList[j][0]-xx)<dist && abs(FullPeakList[j][1]-yy)<dist)
//				tooClose = 1
//			endif
//		endfor
//		if (tooClose)
//			continue
//		endif
//		FullPeakList[Nu][0]=xx				;	FullPeakList[Nu][1]=yy
//		FullPeakList[Nu][2]=sigX				;	FullPeakList[Nu][3]=sigY
//		FullPeakList[Nu][4]=fwx				;	FullPeakList[Nu][5]=fwy
//		FullPeakList[Nu][6]=W_sigma[3]*fw	;	FullPeakList[Nu][7]=W_sigma[5]*fw
//		FullPeakList[Nu][8]=W_coef[6]		;	FullPeakList[Nu][9]=W_sigma[6]
//		FullPeakList[Nu][10]=W_coef[1]*fwx*fwy
//		sprintf str,"%.0f,%.0f;",xx,yy
//		pkList += str
//		Nu += 1
//	endfor
//	Redimension/N=(Nu,11) FullPeakList
//	KillWaves/Z W_sigma,W_coef,W_ParamConfidenceInterval
//	return Nu
//End



//Static Function resetpkLIstFromFullPeakList(FullPeakList)	// reset the pkList string from data in FullPeakList
//	Wave FullPeakList
//	if (!WaveInClass(FullPeakList,"FittedPeakList*"))
//		return 1
//	endif
//
//	String fldr = ParseFilePath(1,StringByKey("fittedIgorImage", note(FullPeakList),"="),":",-1,0)
//	if (strlen(fldr)<1)
//		return 1
//	endif
//	SVAR pkList = $(fldr+"pkList")
//	if (!SVAR_Exists(pkList))
//		String/G $(fldr+"pkList")
//	endif
//
//	pkList = ""
//	Variable N=DimSize(FullPeakList,0), i
//	String str
//	Variable px,py
//	for (i=0;i<N;i+=1)
//		px = FullPeakList[i][0]
//		py = FullPeakList[i][1]
//		if (numtype(px+py)==0)
//			sprintf str,"%.0f,%.0f;",px,py
//			pkLIst += str
//		endif
//	endfor
//	return 0
//End



Static Function FullPeakList2Qfile(FullPeakList,fname,pathName,[FullPeakList1,FullPeakList2,depth])	// convert peaks to a Qlist+intens, peaks file for analysis by Euler
	Wave FullPeakList
	String fname
	String pathName
	Wave FullPeakList1
	Wave FullPeakList2
	Variable depth				// optional, depth to use (micron)
	if (ParamIsDefault(FullPeakList1))
		Wave FullPeakList1=$""
	endif
	if (ParamIsDefault(FullPeakList2))
		Wave FullPeakList2=$""
	endif
	depth = ParamIsDefault(depth) || numtype(depth) ? NaN : depth

	fname = SelectString(strlen(fname),"generic_Peaks.txt",fname)
	if (!WaveExists(FullPeakList))
		DoAlert 0, "input wave for FullPeakList2File() does not exists"
		return 1
	elseif (DimSize(FullPeakList,1)<11)
		DoAlert 0, "the passed full peak list '"+NameOfWave(FullPeakList)+"' is not the right size"
		return 1
	elseif (WaveExists(FullPeakList1))
		if (DimSize(FullPeakList1,1)<11)
			DoAlert 0, "Full peak list 1'"+NameOfWave(FullPeakList1)+"' is the wrong size"
			return 1
		endif
	elseif (WaveExists(FullPeakList2))
		if (DimSize(FullPeakList2,1)<11)
			DoAlert 0, "Full peak list 1'"+NameOfWave(FullPeakList2)+"' is the wrong size"
			return 1
		endif
	endif

	Variable N0=DimSize(FullPeakList,0), N1=0,N2=0
	if (WaveExists(FullPeakList1))
		N1 = DimSize(FullPeakList1,0)
	endif
	if (WaveExists(FullPeakList2))
		N2 = DimSize(FullPeakList2,0)
	endif
	Variable N=N0+N1+N2
	if (N<1)												// nothing to write
		return 1
	endif

	STRUCT microGeometry geo								// note, dd and yc are reset from wave note below if it exists
	if (FillGeometryStructDefault(geo, alert=1))	//fill the geometry structure with default values
		return 1
	endif
	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))
		DoAlert 0, "no lattice structure found, did you forget to set it?"
		return 1
	endif

	Make/N=(N,4)/O/FREE Qs								// temporary waves for the q info, holds Qhat and intensity
	Make/N=3/D/FREE qhat
	String wnote = note(FullPeakList)
	Variable dNum = max(detectorNumFromID(geo, StringByKey("detectorID", wnote,"=")),0)
	Variable startx,groupx, starty,groupy			// ROI of the original image
	startx = NumberByKey("startx",wnote,"=")
	groupx = NumberByKey("groupx",wnote,"=")
	starty = NumberByKey("starty",wnote,"=")
	groupy = NumberByKey("groupy",wnote,"=")
	startx = numtype(startx) ? FIRST_PIXEL : startx
	groupx = numtype(groupx) ? 1 : groupx
	starty = numtype(starty) ? FIRST_PIXEL : starty
	groupy = numtype(groupy) ? 1 : groupy
	String fittedIgorImages=StringByKey("fittedIgorImage",wnote,"=")
	String PeakListWaves=GetWavesDataFolder(FullPeakList,2)

	Variable i, px,py
	Variable depthFile = numtype(depth) ? NumberByKey("depth",wnote,"=") : depth
	for (i=0,N=0;i<N0;i+=1)
		px = (startx-FIRST_PIXEL) + groupx*FullPeakList[i][0] + (groupx-1)/2		// change to un-binned pixels
		py = (starty-FIRST_PIXEL) + groupy*FullPeakList[i][1] + (groupy-1)/2		// pixels are still zero based
		pixel2q(geo.d[dNum],px,py,qhat,depth=depthFile)		// was in Wenge-coord system (OLD) or BeamLine(New)
		if (norm(qhat)>0)									// check for a valid Q
			Qs[N][0,2] = qhat[q]
			Qs[N][3] = FullPeakList[i][10]			// the integral
			N += 1
		endif
	endfor

	if (WaveExists(FullPeakList1))
		wnote = note(FullPeakList1)
		dNum = max(detectorNumFromID(geo, StringByKey("detectorID", wnote,"=")),0)
		startx = NumberByKey("startx",wnote,"=")
		groupx = NumberByKey("groupx",wnote,"=")
		starty = NumberByKey("starty",wnote,"=")
		groupy = NumberByKey("groupy",wnote,"=")
		startx = numtype(startx) ? FIRST_PIXEL : startx
		groupx = numtype(groupx) ? 1 : groupx
		starty = numtype(starty) ? FIRST_PIXEL : starty
		groupy = numtype(groupy) ? 1 : groupy
		depthFile = numtype(depth) ? NumberByKey("depth",wnote,"=") : depth
		fittedIgorImages += ","+StringByKey("fittedIgorImage",wnote,"=")
		PeakListWaves += ","+GetWavesDataFolder(FullPeakList1,2)
		for (i=0;i<N1;i+=1)
			px = (startx-FIRST_PIXEL) + groupx*FullPeakList1[i][0] + (groupx-1)/2		// change to un-binned pixels
			py = (starty-FIRST_PIXEL) + groupy*FullPeakList1[i][1] + (groupy-1)/2		// pixels are still zero based
			pixel2q(geo.d[dNum],px,py,qhat,depth=depthFile)		// was in Wenge-coord system (OLD) or BeamLine(New)
			if (norm(qhat)>0)								// check for a valid Q
				Qs[N][0,2] = qhat[q]
				Qs[N][3] = FullPeakList1[i][10]		// the integral
				N += 1
			endif
		endfor
	endif

	if (WaveExists(FullPeakList2))
		wnote = note(FullPeakList2)
		dNum = max(detectorNumFromID(geo, StringByKey("detectorID", wnote,"=")),0)
		startx = NumberByKey("startx",wnote,"=")
		groupx = NumberByKey("groupx",wnote,"=")
		starty = NumberByKey("starty",wnote,"=")
		groupy = NumberByKey("groupy",wnote,"=")
		startx = numtype(startx) ? FIRST_PIXEL : startx
		groupx = numtype(groupx) ? 1 : groupx
		starty = numtype(starty) ? FIRST_PIXEL : starty
		groupy = numtype(groupy) ? 1 : groupy
		depthFile = numtype(depth) ? NumberByKey("depth",wnote,"=") : depth
		fittedIgorImages += ","+StringByKey("fittedIgorImage",wnote,"=")
		PeakListWaves += ","+GetWavesDataFolder(FullPeakList2,2)
		for (i=0;i<N2;i+=1)
			px = (startx-FIRST_PIXEL) + groupx*FullPeakList2[i][0] + (groupx-1)/2		// change to un-binned pixels
			py = (starty-FIRST_PIXEL) + groupy*FullPeakList2[i][1] + (groupy-1)/2		// pixels are still zero based
			pixel2q(geo.d[dNum],px,py,qhat,depth=depthFile)		// was in Wenge-coord system (OLD) or BeamLine(New)
			if (norm(qhat)>0)								// check for a valid Q
				Qs[N][0,2] = qhat[q]
				Qs[N][3] = FullPeakList2[i][10]		// the integral
				N += 1
			endif
		endfor
	endif
	// N is now the actual number of Qs to write

	wnote = note(FullPeakList)						// reset wnote to first peak list
	startx = NumberByKey("startx",wnote,"=")	// and reset the values writting to first peak file too
	groupx = NumberByKey("groupx",wnote,"=")
	starty = NumberByKey("starty",wnote,"=")
	groupy = NumberByKey("groupy",wnote,"=")
	startx = numtype(startx) ? FIRST_PIXEL : startx
	groupx = numtype(groupx) ? 1 : groupx
	starty = numtype(starty) ? FIRST_PIXEL : starty
	groupy = numtype(groupy) ? 1 : groupy
	do
		fittedIgorImages = ReplaceString(",,",fittedIgorImages,"")
	while (strsearch(fittedIgorImages,",,",0)>0)
	fittedIgorImages = TrimBoth(fittedIgorImages,chars=",")
	do
		PeakListWaves = ReplaceString(",,",PeakListWaves,"")
	while (strsearch(PeakListWaves,",,",0)>0)
	PeakListWaves = TrimBoth(PeakListWaves,chars=",")

	Variable refNum
	Open/C="R*ch"/P=$pathName/T="TEXT" refNum as fname
	if (!strlen(S_fileName))
		DoAlert 0, "Cannot open file '"+fname+"'"
		return 1
	endif

	String str
	Variable val
	fprintf refNum,"$PeaksFile\n"
	str = xtal.desc
	if (strlen(str))
		fprintf refNum,"$structureDesc		%s\n",xtal.desc
	endif
	fprintf refNum,"$latticeParameters	{ %g, %g, %g, %g, %g, %g }// using nm and degrees\n",xtal.a,xtal.b,xtal.c,xtal.alpha,xtal.beta,xtal.gam
	fprintf refNum,"$lengthUnit			nm					// length unit for lattice constants a,b,c\n"
	fprintf refNum,"$SpaceGroupID			%s					// Structure number from International Tables\n",xtal.SpaceGroupID
	fprintf refNum,"$SpaceGroupIDnum		%d					// Structure number from International Tables\n",xtal.SpaceGroupIDnum
	fprintf refNum,"$SpaceGroup			%d					// Structure number from International Tables\n",xtal.SpaceGroup
	if (stringmatch(StringFromList(1,GetIndexingFuncOrExec()),"EulerOrig"))
		fprintf refNum,"$latticeStructure		%d					// Structure number from International Tables\n",xtal.SpaceGroup
	endif
	str = "\t// {element  x y z occupancy}"
	for (i=0;i<xtal.N;i+=1)
		fprintf refNum,"$AtomDesctiption%d	{%s  %g %g %g %g}%s\n",i+1,xtal.atom[i].name,xtal.atom[i].x,xtal.atom[i].y,xtal.atom[i].z,xtal.atom[i].occ,str
		str = ""
	endfor

	str = StringByKey("imageFileName",wnote,"=")
	if (strlen(str))
		fprintf refNum,"$imageFileName			%s		// name of image file read in by Igor\n",str
	endif
//	fprintf refNum,"$peakListWave	%s	// Igor wave containng list of peak positions\n",GetWavesDataFolder(FullPeakList,2)
	fprintf refNum,"$peakListWave	%s	// Igor wave containng list of peak positions\n",PeakListWaves
	val = NumberByKey("exposure",wnote,"=")
	if (!numtype(val))
		fprintf refNum,"$exposure			%g					// CCD exposure of original image (sec)\n",val
	endif
	str = StringByKey("dateExposed",wnote,"=")
	if (strlen(str))
		fprintf refNum,"$dateExposed		%s			// date of original CCD image\n",str
	endif
	str = StringByKey("rawIgorImage",wnote,"=")
	if (strlen(str))
		fprintf refNum,"$rawIgorImage		%s		// name of Igor image wave\n",str
	endif
//	str = StringByKey("fittedIgorImage",wnote,"=")
	if (strlen(fittedIgorImages))
		fprintf refNum,"$fittedIgorImage	%s	// name of flattened Igor image wave used to fit peaks\n",fittedIgorImages
	endif
	val = NumberByKey("fractionBkg",wnote,"=")
	if (!numtype(val))
		fprintf refNum,"$fractionBkg		%g					// in peak analysis, the given fraction of image that is bkg\n",val
	endif
	val = NumberByKey("minSpotSeparation",wnote,"=")
	if (!numtype(val))
		fprintf refNum,"$minSpotSeparation	%g					// in peak analysis, min separaion between spots (pixels)\n",val
	endif
	val = NumberByKey("minPeakWidth",wnote,"=")
	if (!numtype(val))
		fprintf refNum,"$minPeakWidth	%g					// in peak analysis, min allowed fwhm for a spot (pixels)\n",val
	endif
	val = NumberByKey("maxPeakWidth",wnote,"=")
	if (!numtype(val))
		fprintf refNum,"$maxPeakWidth	%g					// in peak analysis, max allowed fwhm for a spot (pixels)\n",val
	endif
	val = NumberByKey("totalPeakIntensity",wnote,"=")
	if (!numtype(val))
		fprintf refNum,"$totalPeakIntensity	%g					// sum of the areas of all fitted peaks\n",val
	endif
	val = NumberByKey("totalIntensity",wnote,"=")
	if (!numtype(val))
		fprintf refNum,"$totalIntensity	%g						// sum of all the intensity in the image\n",val
	endif
	val = NumberByKey("startx",wnote,"=")					// write the ROI of the image
	if (!numtype(val))
		fprintf refNum,"$startx	%g								// ROI of the raw image file\n",val
	endif
	val = NumberByKey("groupx",wnote,"=")
	if (!numtype(val))
		fprintf refNum,"$groupx	%g\n",val
	endif
	val = NumberByKey("endx",wnote,"=")
	if (!numtype(val))
		fprintf refNum,"$endx	%g\n",val
	endif
	val = NumberByKey("starty",wnote,"=")
	if (!numtype(val))
		fprintf refNum,"$starty	%g\n",val
	endif
	val = NumberByKey("groupy",wnote,"=")
	if (!numtype(val))
		fprintf refNum,"$groupy	%g\n",val
	endif
	val = NumberByKey("endy",wnote,"=")
	if (!numtype(val))
		fprintf refNum,"$endy	%g\n",val
	endif
	val = NumberByKey("CCDy",wnote,"=")
	if (!numtype(val))
		fprintf refNum,"$CCDy	%g\n",val
	endif
	val = numtype(depth) ? NumberByKey("depth",wnote,"=") : depth
	if (!numtype(val))
		fprintf refNum,"$depth	%g\n",val
	endif

	fprintf refNum,"\n// the following table contains xyz compotnents of G^ and the integral of the peak\n"
	fprintf refNum,"$N_Ghat+Intens 	%d		// number of G^ vectors\n",N
	for (i=0;i<N;i+=1)
		fprintf refNum, "% .7f, % .7f, % .7f,     %g\n",Qs[i][0],Qs[i][1],Qs[i][2],Qs[i][3]
	endfor
	Close refNum
	return 0
End


Function/T MakeMeasured_hkl_EnergiesWave(N,type)
	Variable N
	String type		// type of info measured, one of "keV", "Qnm", "dnm"

	String typeList = "keV;Qnm;dnm"
	if (!(N>0) || WhichListItem(type,typeList)<0)
		if (WaveExists(hklEmeasured))
			N = N>0 ? N : DimSize(hklEmeasured,0)
			if (WhichListItem(type,typeList)<0)
				type = GetDimLabel(hklEmeasured,1,3)
				type = ReplaceString("_measured",type,"")
			endif
		endif
		N = (N>0) ? min(N,100) : 10
		Prompt N,"number of measured energies"
		Prompt type,"type of measured length info",popup,typeList
		DoPrompt SelectString(WaveExists(hklEmeasured),"","EXISTING ")+"measured Energies",N,type
		if (V_flag)
			return ""
		endif
	endif
	if (!(N>0) || WhichListItem(type,typeList)<0)
		return ""
	endif

	Variable fresh = !WaveExists($"hklEmeasured")
	Make/N=(N,7)/O hklEmeasured							// holds hkl,E,px,py,detNum for points that have the energies measured
	if (fresh)
		hklEmeasured = NaN
	endif
	hklEmeasured[][3] = hklEmeasured[p][3]>0 ? hklEmeasured[p][3] : NaN
	hklEmeasured[][4] = hklEmeasured[p][4]>0 ? hklEmeasured[p][4] : NaN
	hklEmeasured[][5] = hklEmeasured[p][5]>0 ? hklEmeasured[p][5] : NaN
	hklEmeasured[][6] = numtype(hklEmeasured[p][6]) ? 0 : hklEmeasured[p][6]
	SetDimLabel 1,0,H,hklEmeasured ;			SetDimLabel 1,1,K,hklEmeasured ;		SetDimLabel 1,2,L,hklEmeasured
	SetDimLabel 1,3,$(type+"_measured"),hklEmeasured
	SetDimLabel 1,4,pixelX,hklEmeasured ;		SetDimLabel 1,5,pixelY,hklEmeasured;	SetDimLabel 1,6,detNum,hklEmeasured
	Note/K hklEmeasured,"waveClass=measuredEnergy;type="+type+";"
	CheckDisplayed/A hklEmeasured
	if (V_flag==0)												// not displayed, make the table
		DisplayTableOfWave(hklEmeasured)
		SetWindow kwTopWin,hook(fill)=MeasuredFillHook
	else
		String win,wins=WinList("*",";","WIN:2")
		Variable i
		for (i=0;i<ItemsInList(wins);i+=1)
			win = StringFromList(i,wins)
			CheckDisplayed/W=$win hklEmeasured
			if (V_flag)
				DoWindow/F $win
				break
			endif
		endfor
	endif
	return GetWavesDataFolder(hklEmeasured,2)
End
//
Function MeasuredFillHook(s)
	STRUCT WMWinHookStruct &s

	Variable hookResult = 0
	if (s.eventCode!=2)
		return 0
	endif
	DoAlert 1,"Should I fill in the rest of the table?"
	if (V_flag!=1)
		return 0
	endif
	String win=s.winName
	Wave ww = WaveRefIndexed(win,0,1)
	if (WaveExists(ww))
		printf "Completing columns in '%s'\r",NameOfWave(ww)
		Complete_hklEmeasured(ww)
	endif
	return hookResult		// 0 if nothing done, else 1
End
//
// Fill in extra columns in hklEmeasured
Function Complete_hklEmeasured(hklEmeasured)
	Wave hklEmeasured
	String wList = WaveListClass("IndexedPeakList*","*","")

	// order by moddate()
	Variable Ni=ItemsInList(wList)
	Make/N=(Ni)/WAVE/FREE waves
	Make/N=(Ni)/D/FREE times
	waves = $StringFromList(p,wList)
	times = modDate(waves[p])
	Sort/R times, times,waves

	Variable col_pixelX = colFromDimLabel(hklEmeasured,"pixelX")
	Variable col_pixelY = colFromDimLabel(hklEmeasured,"pixelY")
	Variable col_detNum = colFromDimLabel(hklEmeasured,"detNum")
  	Make/N=3/D/FREE hkl
	Variable detNum,px,py, ipeak
	Variable m,Nm=Dimsize(hklEmeasured,0), N,i
	for (m=0;m<Nm;m+=1)			// loop over each of the measured spots
		hkl = hklEmeasured[m][p]
		for (i=0;i<Ni;i+=1)
			Wave indexWave = waves[i]
			N = DimSize(indexWave,0)
			Make/N=(N,3)/FREE/D/O dhkls

			dhkls = indexWave[p][q+3] - hkl[q]
			MatrixOP/FREE/O dmags = sumRows(dhkls*dhkls)
			WaveStats/M=1/Q  dmags
			if (V_min<1e-9)
				break
			endif
		endfor
		if (V_min>1e-9)
			continue
		endif

		Wave peaks = $StringByKey("peakListWave",note(waves[i]),"=")
		if (!WaveExists(peaks))
			continue
		endif

		STRUCT microGeometry geo
		if (FillGeometryStructDefault(geo, alert=1))	//fill the geometry structure with test values
			return 1
		endif
		detNum = detectorNumFromID(geo, StringByKey("detectorID",note(peaks),"="))
		if (col_detNum>=0 && detNum>=0)
			hklEmeasured[m][col_detNum]=detNum
		endif
		px = indexWave[V_minloc][9]
		py = indexWave[V_minloc][10]
		ipeak = FindPeakInPeaksList(peaks,px,py,1)
		px = peaks[ipeak][0]
		px = peaks[ipeak][1]
		if ((col_pixelX+col_pixelY)>0 && ipeak>=0)
			hklEmeasured[m][col_pixelX]=px
			hklEmeasured[m][col_pixelY]=py
		endif
		printf "hkl = %s,  indexWave=%s[%d],  detNum=%d, ipeak = %d,  pixel=(%g, %g)\r",vec2str(hkl),NameOfWave(waves[i]),V_minloc,detNum,ipeak,px,py
	endfor
End
//
Static Function FindPeakInPeaksList(peaks,px,py,maxDist)	// find closest peak to (px,py) from the peaks file
	Wave peaks
	Variable px,py
	Variable maxDist										// only accept if pixel is within a distance of maxDist

	Variable N=DimSize(peaks,0)
	Make/N=(N)/FREE xy
	xy[][0] = peaks[p][0]-px
	xy[][1] = peaks[p][1]-py
	MatrixOP/FREE dist2 = sumRows(xy*xy)
	WaveStats/M=1/Q dist2
	if ((V_min)>maxDist^2)
		return NaN
	endif
	return V_minloc
End
//
Static Function colFromDimLabel(w,colLabel)				// return column number with given DimLabel
	Wave w
	String colLabel											// column label to search for
	Variable i,N=DimSize(w,1)
	for (i=0;i<N;i+=1)
		if (stringmatch(GetDimLabel(hklEmeasured,1,i),colLabel))
			return i
		endif
	endfor
	return NaN
End
//
//Function/T MakeMeasured_hkl_EnergiesWave(N,type)
//	Variable N
//	String type
//	if (!(N>0))
//		if (WaveExists(hklEmeasured))
//			N = DimSize(hklEmeasured,0)
//		endif
//		N = (N>0) ? min(N,100) : 10
//		Prompt N,"number of measured energies"
//		DoPrompt "measured Energies",N
//		if (V_flag)
//			return ""
//		endif
//	endif
//	if (!(N>0))
//		return ""
//	endif
//	Make/N=(N,7)/O hklEmeasured								// holds hkl,E,px,py,detNum for points that have the energies measured
//	SetDimLabel 1,0,H,hklEmeasured ;			SetDimLabel 1,1,K,hklEmeasured ;		SetDimLabel 1,2,L,hklEmeasured
//	SetDimLabel 1,3,keVmeasured,hklEmeasured
//	SetDimLabel 1,4,pixelX,hklEmeasured ;		SetDimLabel 1,5,pixelY,hklEmeasured;	SetDimLabel 1,6,detNum,hklEmeasured
//	Note/K hklEmeasured,"waveClass=measuredEnergy;"
//
//	CheckDisplayed/A hklEmeasured
//	if (V_flag==0)												// not displayed, make the table
//		Edit/W=(5,44,666,483)/K=1 hklEmeasured.ld
//		ModifyTable format(Point)=1,width(Point)=28,width(hklEmeasured.l)=40
//	else
//		String win,wins=WinList("*",";","WIN:2")
//		Variable i
//		for (i=0;i<ItemsInList(wins);i+=1)
//			win = StringFromList(i,wins)
//			CheckDisplayed/W=$win hklEmeasured
//			if (V_flag)
//				DoWindow/F $win
//				break
//			endif
//		endfor
//	endif
//	return GetWavesDataFolder(hklEmeasured,2)
//End
//




// find best hkl that is most paralllel to the input and allowed
Static Function CloseAllowedhkl(hklIN,[qmax])	// in and out can be the same wave
	Wave hklIN						// real numbers on intput, set to h,k,l on output
	Variable qmax					// maximum allowed value of |q| (1/nm), used to determine maxFactor
	qmax = ParamIsDefault(qmax) || numtype(qmax) || qmax<=0 ? 230 : qmax	// 230(1/nm) puts 32keV at 2theta=90°
	//	Variable lambda = 4*PI/sqrt(2)/qmax, keV=hc/lambda
	//	print "lambda =",lambda, "   ",keV

	STRUCT crystalStructure xtal			// temporary crystal structure
	FillCrystalStructDefault(xtal)			//fill the lattice structure with default values

	WaveStats/Q/M=1 hklIN
	Variable maxVal=max(V_max,-V_min)
	Make/N=3/D/FREE hkl1 = hklIN/maxVal	// make largest value in |hkl1| = 1

	Wave recip = recipFrom_xtal(xtal)		// FREE wave with reciprocal lattice
	MatrixOP/FREE qvec0 = recip x hkl1	// qvec for hkl1
	Variable maxFactor = max(qmax/norm(qvec0), 2.1)	// I can multiply hkl1 by maxFactor without exceeding qmax
	Variable N=floor(maxFactor)
	MatrixOP/FREE hkls = colRepeat(hkl1,N)
	hkls *= (q+1)

	MatrixOP/FREE dots = sumCols( normalizeCols(recip x round(hkls)) * colRepeat(normalize(qvec0),N) )
	Redimension/N=(N) dots					// want (N), not (1,N)
	dots += (N-1-p)*SMALLEST1				// selects the lowest index when getting V_maxloc
	WaveStats/Q/M=1 dots
	Variable ibest = limit(V_maxloc+1,1,N)
	Variable xx=round(ibest*hkl1[0]), yy=round(ibest*hkl1[1]), zz=round(ibest*hkl1[2])
	lowestAllowedHKL(xx,yy,zz)
	hklIN = {xx,yy,zz}
	return 0
End
//
//// find hkl that is paralllel to the input and allowed
//Static Function CloseAllowedhkl(hkl,[maxMultiply])	// in and out can be the same wave
//	Wave hkl							// real numbers on intput, set to h,k,l on output
//	Variable maxMultiply				// used to determine maxFactor, maxFactor = maxMultiply*60
//	maxMultiply = ParamIsDefault(maxMultiply) ? 1 : maxMultiply
//	maxMultiply = maxMultiply>0 && numtype(maxMultiply)==0 ? maxMultiply : 1
//
//	STRUCT crystalStructure xtal			// temporary crystal structure
//	FillCrystalStructDefault(xtal)			//fill the lattice structure with default values
//	Variable maxFactor = round(max(max(xtal.a,xtal.b),xtal.c)*maxMultiply*60)
//	Variable xin=hkl[0],yin=hkl[1],zin=hkl[2]		// save input values
//	hkl = abs(hkl)
//	Variable i, maxVal = WaveMax(hkl)
//	Variable ibest=0, dotMax=-Inf, dot
//	for (i=1;i<=maxFactor;i+=1)			// maxFactor was 25, which was good for Si
//		hkl = {xin,yin,zin}
//		hkl = round(hkl*i/maxVal)
//		dot = (hkl[0]*xin + hkl[1]*yin + hkl[2]*zin)/norm(hkl)
//		if (dot > (dotMax+SMALLEST1))
//			ibest = i
//			dotMax = dot
//		endif
//	endfor
//	hkl = {xin,yin,zin}
//	hkl = round(hkl*ibest/maxVal)
//	xin = hkl[0]
//	yin = hkl[1]
//	zin = hkl[2]
//	lowestAllowedHKL(xin,yin,zin)
//	hkl = {xin,yin,zin}
//	hkl = hkl[p]==0 ? 0 : hkl[p]		// avoid "-0"
//	return 0
//End
//Static Function CloseAllowedhkl(hkl)	// in and out can be the same wave
//	Wave hkl							// real numbers on intput, set to h,k,l on output
//
//	Variable xin=hkl[0],yin=hkl[1],zin=hkl[2]		// save input values
//	hkl = abs(hkl)
//	Variable i, maxVal = WaveMax(hkl)
//	Variable ibest=0, dotMax=-Inf, dot
//	for (i=1;i<=25;i+=1)
//		hkl = {xin,yin,zin}
//		hkl = round(hkl*i/maxVal)
//		dot = (hkl[0]*xin + hkl[1]*yin + hkl[2]*zin)/norm(hkl)
//		if (dot > (dotMax+SMALLEST1))
//			ibest = i
//			dotMax = dot
//		endif
//	endfor
//	hkl = {xin,yin,zin}
//	hkl = round(hkl*ibest/maxVal)
//	xin = hkl[0]
//	yin = hkl[1]
//	zin = hkl[2]
//	lowestAllowedHKL(xin,yin,zin)
//	hkl = {xin,yin,zin}
//	hkl = hkl[p]==0 ? 0 : hkl[p]		// avoid "-0"
//	return 0
//End
//Static Function CloseAllowedhkl(hkl)	// in and out can be the same wave
//	Wave hkl							// real numbers on intput, set to h,k,l on output
//
//	Variable xin=hkl[0],yin=hkl[1],zin=hkl[2]		// save input values
//	hkl = abs(hkl)
//	Variable i, maxVal = WaveMax(hkl)
//	Variable ibest=0, dotMax=-Inf, dot
//	for (i=1;i<=25;i+=1)
//		hkl = {xin,yin,zin}
//		hkl = round(hkl*i/maxVal)
//		dot = (hkl[0]*xin + hkl[1]*yin + hkl[2]*zin)/norm(hkl)
//		if (dot > (dotMax+SMALLEST1))
//			ibest = i
//			dotMax = dot
//		endif
//	endfor
//	hkl = {xin,yin,zin}
//	hkl = round(hkl*ibest/maxVal)
//	hkl = hkl[p]==0 ? 0 : hkl[p]		// avoid "-0"
//	return 0
//End
//Function smallest()
//	Variable hi=1,lo=0, mid, best
//	do
//		best = 1-mid
//		mid = (hi+lo)/2
//		if ((1-mid) != 1)
//			lo = mid
//		else
//			hi =mid
//		endif
//	while(lo!=hi)
//	printf "1 - %g ,  is as close as you can get to one\r",best
//End
//
//	•smallest()
//	  1 - 1.11022e-16 ,  is as close as you can get to one
//
//
//
//Static  Function CloseAllowedhkl(hkl)	// in and out can be the same wave
//	Wave hkl							// real numbers on intput, set to h,k,l on output
//
//	Variable xin=hkl[0],yin=hkl[1],zin=hkl[2]		// save input values
//	hkl = abs(hkl)
//	Variable minNonZeroVal = WaveMax(hkl)/60	// if it is less than this, then it is a zero
//	hkl = (hkl[p]<minNonZeroVal) ? Inf : hkl[p]
//	WaveStats/Q/M=1 hkl							// do this to get V_minloc
//
//	hkl = {xin,yin,zin}
//	Variable minh = abs(hkl[V_minloc])			// location of smallest non-zero value, e.g. for {0,2,7}, minh=2
//	hkl = round(hkl[p]/minh)
//	hkl = hkl[p]==0 ? 0 : hkl[p]					// get rid of "-0"
//
//	xin=hkl[0]; yin=hkl[1]; zin=hkl[2]
//	lowestAllowedHKL(xin,yin,zin)					// convert hkl to allowed reflection
//	hkl = {xin,yin,zin}
//	return 0
//End
//Static  Function CloseAllowedhkl(in,out)	// in and out can be the same wave
//	Wave in,out								// input and output hkl
//	String wName=""
//	if (stringmatch(GetWavesDataFolder(in,2),GetWavesDataFolder(out,2)))
//		wName = UniqueName("CloseAllowedhkl",1,0)
//		Make/N=3/O/D $wName=in
//		Wave inTemp = $wName
//	else
//		Wave inTemp=in
//	endif
//
//	out = abs(inTemp)
//	Variable minNonZeroVal = WaveMax(out)/60	// if it is less than this, then it is a zero
//	out = (out[p]<minNonZeroVal) ? Inf : out[p]
//	WaveStats/Q/M=1 out
//	Variable minh = abs(inTemp[V_minloc])		// location of smallest non-zero value, e.g. for {0,2,7}, minh=2
//	out = round(inTemp[p]/minh)
//	out = out[p]==0 ? 0 : out[p]				// get rid of "-0"
//	Variable h=out[0], k=out[1], l=out[2]
//	lowestAllowedHKL(h,k,l)					// convert hkl to allowed reflection
//	out = {h,k,l}
//	KillWaves/Z $wName
//	return 0
//End
//Static  Function CloseAllowedhkl(in,out)
//	Wave in,out								// input and output hkl
//
//	out = abs(in)
//	Variable minNonZeroVal = WaveMax(out)/60	// if it is less than this, then it is a zero
//	out = (out[p]<minNonZeroVal) ? Inf : out[p]
//	WaveStats/Q/M=1 out
//	Variable minh = abs(in[V_minloc])		// location of smallest non-zero value, e.g. for {0,2,7}, minh=2
//	out = round(in[p]/minh)
//	out = out[p]==0 ? 0 : out[p]				// get rid of "-0"
//	Variable h=out[0], k=out[1], l=out[2]
//	lowestAllowedHKL(h,k,l)					// convert hkl to allowed reflection
//	out = {h,k,l}
//	return 0
//End
//Static Function CloseAllowedhkl(in,out)
//	Wave in,out								// input and output hkl
//	out = abs(round(in*60)/60)
//	out = out[p]>0 ?out[p] : Inf
//	WaveStats/Q/M=1 out
//	Variable minh = abs(in[V_minloc])		// smallest non-zero value, e.g. for {0,2,7}, minh=2
//	out = round(in[p]/minh*12)
//	out = out[p]==0 ? 0 : out[p]				// get rid of "-0"
//	Variable h=out[0], k=out[1], l=out[2]
//	lowestAllowedHKL(h,k,l)					// convert hkl to allowed reflection
//	out = {h,k,l}
//	return 0
//End
//// find hkl that is paralllel to the input and allowed
//Static Function CloseAllowedhkl(SpaceGroup,in,out)
//	Wave in,out								// input and output hkl
//	Variable SpaceGroup					// Space Group number from International tables
//	out = abs(round(in*60)/60)
//	out = out[p]>0 ?out[p] : Inf
//	WaveStats/Q/M=1 out
//	Variable minh = abs(in[V_minloc])		// smallest non-zero value, e.g. for {0,2,7}, minh=2
//	out = round(in[p]/minh*12)
//	out = out[p]==0 ? 0 : out[p]				// get rid of "-0"
//	Variable h=out[0], k=out[1], l=out[2]
//	lowestAllowedHKL(h,k,l,SpaceGroup)	// convert hkl to allowed reflection
//	out = {h,k,l}
//	return 0
//End



Function EnergyOfhkl(FullPeakIndexed,pattern,h,k,l,[printIt])
	Wave FullPeakIndexed
	Variable pattern
	Variable h,k,l
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt

	Variable i
	if (!WaveExists(FullPeakIndexed))
		String wList = reverseList(WaveListClass("IndexedPeakList*","*",""))
		i = ItemsInList(wList)
			if (i==1)						// only one, so choose it
			Wave FullPeakIndexed = $(StringFromList(0,wList))
		elseif(i>0)						// more than 1, so ask
			String wName=""
			Prompt wName,"Indexed Peak List",popup,wList
			DoPrompt "Indexed Peak List",wName
			if (V_flag)
				return NaN
			endif
			Wave FullPeakIndexed = $wName
		endif
	endif
	if (!WaveExists(FullPeakIndexed))
		if (printIt)
			print "cannot find wave of indexed peak information 'FullPeakIndexed'"
		endif
		return NaN
	endif

	Variable Np = max(DimSize(FullPeakIndexed,2),1)
	pattern = (Np==1) ? 0 : pattern



	// need to fix this so it takes the user data from the graph, not the panel
	if (numtype(pattern))
		pattern = NumberByKey("patternNum",Getuserdata("","","Indexing"),"=")
	endif



	if (Np>1 && (pattern<0 || pattern>=Np || numtype(pattern)))
		pattern = pattern>=0 ? pattern : 0
		Prompt pattern, "pattern to choose [0,"+num2istr(Np-1)+",]",popup,expandRange("0-"+num2istr(Np-1),";")
		DoPrompt "pattern number",pattern
		if (V_flag || !(pattern>=0))
			return NaN
		endif
		pattern -= 1
	endif
	pattern = pattern>=0 ? pattern : 0
	if (numtype(h+k+l))
		Prompt h, "H"
		Prompt k, "K"
		Prompt l, "L"
		DoPrompt "(hkl)",h,k,l
		if (V_flag || numtype(h+k+l))
			return NaN
		endif
	endif
	if (printIt)
		printf "EnergyOfhkl(%s, %g, %g,%g,%g)\r", NameOfWave(FullPeakIndexed),pattern,h,k,l
	endif
	if (numtype(h+k+l))
		return NaN
	endif

	Make/N=3/O/D/FREE ki={0,0,1}
	String wnote = note(FullPeakIndexed)
	Wave recip = decodeMatFromStr(StringByKey("recip_lattice"+num2istr(pattern),wnote,"="))
	if (!WaveExists(recip))
		if (printIt)
			print "Unable to read recip_lattice"
		endif
		return NaN
	endif
	Make/N=3/D/FREE hkl = {h,k,l}
	MatrixOP/FREE qhat = recip x hkl
	Variable d = 2*PI/normalize(qhat)			// d-spacing, and normalized qhat
	Variable sineTheta = -MatrixDot(qhat,ki)// sin(theta) = -ki dot qhat
	if (sineTheta<=0)									// reflection is unreachable
		if (printIt)
			printf "the (%g, %g, %g) reflection is unreachable\r",h,k,l
		endif
		return NaN
	endif
	Variable keV = hc/(2*d*sineTheta)			// energy
	if (printIt)
		Variable px=NaN, py=NaN
		STRUCT microGeometry geo
		if (!FillGeometryStructDefault(geo))	//fill the geometry structure with test values
			Variable dNum = max(detectorNumFromID(geo, StringByKey("detectorID", wnote,"=")),0)
			Variable startx,groupx, starty,groupy, ddLocal, ycLocal
			startx = NumberByKey("startx",wnote,"=")
			groupx = NumberByKey("groupx",wnote,"=")
			starty = NumberByKey("starty",wnote,"=")
			groupy = NumberByKey("groupy",wnote,"=")
			startx = numtype(startx) ? FIRST_PIXEL : startx
			groupx = numtype(groupx) ? 1 : groupx
			starty = numtype(starty) ? FIRST_PIXEL : starty
			groupy = numtype(groupy) ? 1 : groupy
			Variable/C pz = q2pixel(geo.d[dNum],qhat)
			px = (real(pz)-(startx-FIRST_PIXEL)-(groupx-1)/2)/groupx		// change to binned pixels
			py = (imag(pz)-(starty-FIRST_PIXEL)-(groupy-1)/2)/groupy		// pixels are still zero based
		endif
		if (printIt)
			printf "for %s pattern #%d,  d[(%s)] = %.9g nm,   E(%s) = %.4f keV   (theta=%.3f%s)",NameOfWave(FullPeakIndexed),pattern,hkl2str(h,k,l),d,hkl2str(h,k,l),keV,asin(sineTheta)*180/PI,DEGREESIGN
			if (numtype(px+py)==0)
				printf "       should be at pixel [%.2f, %.2f]",px,py
			endif
			printf "\r"
		endif
	endif
	return keV
End



Function EstimateConeAngle(dnum, [image, printIt])
	Variable dnum
	Wave image									// used to get dnum if it is not specified
	Variable printIt
	printIt = ParamIsDefault(printIt) ? NaN : printIt
	printIt = numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : !(!printIt)

	STRUCT microGeometry g
	FillGeometryStructDefault(g)			//fill the geometry structure with current values

	if (dnum==limit(round(dnum),0,MAX_Ndetectors))	// test for a valid dnum
		dnum = g.d[dnum].used ? dnum : NaN
	else
		dnum = NaN
	endif

	if (numtype(dnum))						// I need to find dnum, try from the image ID
		if (!WaveExists(image))
			Wave image = ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";")))
		endif
		dNum = WaveExists(image) ? detectorNumFromID(g, StringByKey("detectorID",note(image),"=")) : NaN
	endif

	if (dnum==limit(round(dnum),0,MAX_Ndetectors))	// test for a valid dnum
		dnum = g.d[dnum].used ? dnum : NaN
	else
		dnum = NaN
	endif

	if (numtype(dnum))						// invalid dnum
		if (printIt)
			printf "Could not determine cone angle\r"
		endif
		return NaN
	endif

	Variable cone=EstimateConeAngleFromDet(g.d[dnum])
	if (printIt)
		printf "cone angle to fill detector %g is probably %.3g%s\r",dnum,cone,DEGREESIGN
	endif
	return cone
End


ThreadSafe Static Function EstimateConeAngleFromDet(d)
	// Returns cone angle (deg) that will fill this detector.
	// This routine assumes that you know the hkl at detector center
	STRUCT detectorGeometry &d			// detector parameters

	Make/N=3/D/FREE q0
	pixel2q(d,(d.Nx-1)/2,(d.Ny-1)/2,q0)// q^ that produces ray to detector center
	if (norm(q0)==0)
		q0 = {0,1,0}
	endif

	Make/D/FREE pxpy = {{0,0.5,1, 0,1, 0,0.5,1}, {0,0,0, 0.5,0.5, 1,1,1}}
	pxpy[][0] *= d.Nx-1						// {pxpy} are 8 points around edge of detector
	pxpy[][1] *= d.Ny-1
	Wave Qs = pixel2qVEC(d,pxpy)			// returns the qhat in beam line coords for the 8 points

	MatrixOP/FREE minDot = minVal(Qs x q0)
	Variable cone = acos(limit(minDot[0],-1,1))*180/PI
	return cone
End


Function RotationBetweenTwoIndexations(windex0,pattern0,windex1,pattern1,[printIt])
	Wave windex0						// FullPeakIndexed wave with first pattern
	Variable pattern0					// pattern number in first FullPeakIndexed wave
	Wave windex1						// FullPeakIndexed wave with second pattern
	Variable pattern1					// pattern number in second FullPeakIndexed wave
	Variable printIt						// used to force printing
	printIt = ParamIsDefault(printIt) ? NaN : !(!printIt)
	printIt = numtype(printIt) ? (strlen(GetRTStackInfo(2))<1 || stringmatch(GetRTStackInfo(2),"hklTagsPopMenuProc")) : printIt

	Variable sameWave=0
	String indexList = WaveListClass("IndexedPeakList*","*","MAXCOLS:12,MINCOLs:12")
	if (!WaveExists(windex0) || !WaveExists(windex1) || !(pattern0>=0) || !(pattern1>=0))
		String wName0=NameOfWave(windex0), wName1=NameOfWave(windex1)
		if (!WaveExists(windex0))
			wName0 = StringFromList(0,indexList)
		endif
		if (!WaveExists(windex1))
			wName1 = StringFromList(0,indexList)
		endif
		pattern0 = pattern0>=0 ? pattern0 : 0
		pattern1 = pattern1>=0 ? pattern1 : 0
		sameWave = stringmatch(wName0,wName1)
		pattern1 += sameWave
		Prompt wName0,"Index Wave with first reciprocal lattice",popup,indexList
		Prompt wName1,"Index Wave with second reciprocal lattice",popup,indexList
		Prompt pattern0,"pattern number in first wave with reciprocal lattice"
		Prompt pattern1,"pattern number in second wave with reciprocal lattice"
		DoPrompt "Pick Indexed Results",wName0,pattern0,wName1,pattern1
		if (V_flag)
			return NaN
		endif
		Wave windex0=$wName0, windex1=$wName1
		printf "•RotationBetweenTwoIndexations(%s,%g, %s,%g)",NameOfWave(windex0),pattern0,NameOfWave(windex1),pattern1
		if (!ParamIsDefault(printIt))
			printf ", printIt=%g",printIt
		endif
		printf ")\r"
	endif
	pattern0 = round(pattern0)
	pattern1 = round(pattern1)
	if (!WaveExists(windex0) || !WaveExists(windex1) || !(pattern0>=0) || !(pattern1>=0))
		return NaN
	endif

	Wave recip0 = matString2mat(StringByKey("recip_lattice"+num2istr(pattern0),note(windex0),"="))
	Wave recip1 = matString2mat(StringByKey("recip_lattice"+num2istr(pattern1),note(windex1),"="))
	if (!WaveExists(recip0) || !WaveExists(recip1))
		return NaN
	endif

	MatrixOP/FREE rot01 = recip1 x Inv(recip0)	// rotation matrix from 0 to 1
	Make/N=3/D/FREE axis
	Variable angle = axisOfMatrix(rot01,axis,squareUp=1)	// total rotation angle between point0 and point1 (degrees)
	if (printIt)
		MatrixOP/FREE hkl = Normalize(Inv(recip0) x axis)
		Variable h=round(24*hkl[0]), k=round(24*hkl[1]), l=round(24*hkl[2])
		lowestOrderHKL(h,k,l)
		Make/FREE hkl0 = {h,k,l}

		printf "rotation from point %s pattern%d  to  %s pattern%d is about the axis {XYZ} = %s  by  %.3g%s\r",NameOfWave(windex0),pattern0,NameOfWave(windex1),pattern1,vec2str(axis),angle,DEGREESIGN
		printf "in both crystals, the axis is an hkl=%s,    approximately {%s} (off by %.2f%s)\r",vec2str(hkl),hkl2str(h,k,l),angleVec2Vec(hkl0,hkl),DEGREESIGN
	endif
	return angle									// total rotation angle (degrees)
End




// find angle between two Qvectors from points on different images
Function AngleBetweenTwoPointsGeneral(image0,image1, px0,py0, px1,py1)
	Wave image0, image1
	Variable px0,py0				// first pixel (in binning of image)
	Variable px1,py1				// second pixel (in binning of image)

	Variable printIt = (ItemsInList(GetRTStackInfo(0))<2)

	if (!WaveExists(image0))		// use image on top graph if not passed
		Wave image0 = ImageNameToWaveRef("",StringFromList(0,ImageNameList("",";")))
	endif
	if (!WaveExists(image0))
		if (printIt)
			DoAlert 0, "could not find the image with cursors"
		endif
		return 1
	endif
	String gName0 = StringFromList(0,WindowsWithWave(image0,1)), gName1
	if (strlen(gName0)<1)
		return 1
	endif
	if (WaveExists(image1))
		gName1 = StringFromList(0,WindowsWithWave(image1,1))
	else
		gName1 = gName0
		Wave image1 = image0
	endif
	String infoA=CsrInfo(A,gName0), infoB=CsrInfo(B,gName1)
	if ((strlen(infoA)>0 || numtype(px0+py0)==0) && (strlen(infoB)<=0 && numtype(px1+py1)))	// no B point on image0, ask
		String iList=WaveListClass("rawImage*","*","")
		iList = OnlyWavesThatAreDisplayed(iList)
		String iname0,iname1
		if (WaveExists(image0))
			iname0 = NameOfWave(image0)
		else
			iname0 = StringFromList(0,iList)
		endif
		if (WaveExists(image1))
			iname1 = NameOfWave(image1)
		else
			iname1 = StringFromList(1,iList)
		endif
		Prompt iname0,"image to use for cursor A", popup, iList
		Prompt iname1,"image to use for cursor B", popup, iList
		DoPrompt "Select Image(s)",iname0,iname1
		if (V_flag)
			return 1
		endif
		Wave image0 = $iname0
		Wave image1 = $iname1
	endif
	gName0 = StringFromList(0,WindowsWithWave(image0,1))
	gName1 = StringFromList(0,WindowsWithWave(image1,1))
	if (strlen(gName0)==0 || strlen(gName1)==0)
		return 1
	endif
	if (strlen(gName0))
		infoA = CsrInfo(A,gName0)
	else
		infoA = ""
	endif
	if (strlen(gName1))
		infoB = CsrInfo(B,gName1)
	else
		infoB = ""
	endif

	if (numtype(px0+py0) && strlen(infoA))			// bad pixels on image0, use cursor A
		px0 = hcsr(A,gName0)
		py0 = vcsr(A,gName0)
	endif
	if (numtype(px1+py1) && strlen(infoB))			// bad pixels on image1, use cursor B
		px1 = hcsr(B,gName1)
		py1 = vcsr(B,gName1)
	endif

	if (numtype(px0+py0+px1+py1))		// bad pixels, use cursors A and B
		if (printIt)
			DoAlert 0, "could not get the two pixel positions"
		endif
		return 1
	endif
	STRUCT microGeometry geo
	if (FillGeometryStructDefault(geo, alert=printIt))	//fill the geometry structure with test values
		return 1
	endif

	Make/N=3/FREE/D qhat0,qhat1
	Variable angle											// angle between points (degree)
	Variable px,py, depth, dNum
	Variable startx,groupx, starty,groupy					// ROI of the original image

	String wnote = note(image0)
	dNum = max(detectorNumFromID(geo, StringByKey("detectorID", wnote,"=")),0)
	startx = NumberByKey("startx",wnote,"=")
	groupx = NumberByKey("groupx",wnote,"=")
	starty = NumberByKey("starty",wnote,"=")
	groupy = NumberByKey("groupy",wnote,"=")
	startx = numtype(startx) ? FIRST_PIXEL : startx
	groupx = numtype(groupx) ? 1 : groupx
	starty = numtype(starty) ? FIRST_PIXEL : starty
	groupy = numtype(groupy) ? 1 : groupy
	depth = NumberByKey("depth",wnote,"=")
	px = (startx-FIRST_PIXEL) + groupx*px0 + (groupx-1)/2		// change to un-binned pixels
	py = (starty-FIRST_PIXEL) + groupy*py0 + (groupy-1)/2		// pixels are still zero based
	pixel2q(geo.d[dNum],px,py,qhat0,depth=depth)						// in Beam LIne Coord system

	wnote = note(image1)
	dNum = max(detectorNumFromID(geo, StringByKey("detectorID", wnote,"=")),0)
	startx = NumberByKey("startx",wnote,"=")
	groupx = NumberByKey("groupx",wnote,"=")
	starty = NumberByKey("starty",wnote,"=")
	groupy = NumberByKey("groupy",wnote,"=")
	startx = numtype(startx) ? FIRST_PIXEL : startx
	groupx = numtype(groupx) ? 1 : groupx
	starty = numtype(starty) ? FIRST_PIXEL : starty
	groupy = numtype(groupy) ? 1 : groupy
	depth = NumberByKey("depth",wnote,"=")
	px = (startx-FIRST_PIXEL) + groupx*px1 + (groupx-1)/2		// change to un-binned pixels
	py = (starty-FIRST_PIXEL) + groupy*py1 + (groupy-1)/2		// pixels are still zero based
	pixel2q(geo.d[dNum],px,py,qhat1,depth=depth)						// in Beam Line Coord system

	Variable cosine = limit(MatrixDot(qhat0,qhat1),-1,1)
	angle = acos(cosine)*180/PI
	if (printIt)
		if (StringMatch(GetWavesDataFolder(image0,2),GetWavesDataFolder(image1,2)))
			printf "Between points [%.2f, %.2f],  and [%.2f, %.2f]  the angle between Q vectors is %.3f%s\r",px0,py0,px1,py1,angle,DEGREESIGN
		else
			printf "Between points %s[%.2f, %.2f],  and %s[%.2f, %.2f]  the angle between Q vectors is %.3f%s\r",NameOfWave(image0),px0,py0,NameOfWave(image0),px1,py1,angle,DEGREESIGN
		endif
	endif
	return angle
End
//
//		DEPRECATED		DEPRECATED		DEPRECATED		DEPRECATED		Use AngleBetweenTwoPointsGeneral() instead
//
// find angle between two Qvectors from points on an image
Function AngleBetweenTwoPoints(image,px0,py0,px1,py1)
	Wave image
	Variable px0,py0				// first pixel (in binning of image)
	Variable px1,py1				// second pixel (in binning of image)
	Variable angle = AngleBetweenTwoPointsGeneral(image,image, px0,py0, px1,py1)
	return angle
End




Static Function/WAVE matString2mat(str)
	String str
	Variable a0,a1,a2, b0,b1,b2, c0,c1,c2
	sscanf str, "{{%g,%g,%g}{%g,%g,%g}{%g,%g,%g}}", a0,a1,a2, b0,b1,b2, c0,c1,c2
	if (V_flag==9)
		Make/N=(3,3)/D/FREE mat
		mat[0][0]=a0;		mat[0][1]=b0;		mat[0][2]=c0
		mat[1][0]=a1;		mat[1][1]=b1;		mat[1][2]=c1
		mat[2][0]=a2;		mat[2][1]=b2;		mat[2][2]=c2
		return mat
	endif
	return $""
End




// if the image cannot be found in existing graphs, it returns empty string
Function/S FindGraphWithImage(image)	// find the graph window which contains the specified image
	Wave image
	if (!WaveExists(image))
		return ""
	endif
	String name, name0=GetWavesDataFolder(image,2)
	String ilist, win,wlist = WinList("*",";","WIN:1")
	Variable i,Ni, m,Nm=ItemsInList(wlist)
	for (m=0;m<Nm;m+=1)
		win = StringFromList(m,wlist)
		ilist = ImageNameList(win,";")		// list of images in this window
		Ni = ItemsInList(ilist)
		for (i=0;i<ItemsInList(ilist);i+=1)
			name = GetWavesDataFolder(ImageNameToWaveRef(win,StringFromList(i,ilist)),2)
			if (stringmatch(name,name0))
				return win
			endif
		endfor
	endfor
	return ""
End

Static Function/S FindTableWithWave(w0)	// find the table window which contains the specified wave
	Wave w0
	if (!WaveExists(w0))
		return ""
	endif
	String name0=GetWavesDataFolder(w0,2)		// full path name of desired wave
	String table,tlist = WinList("*",";","WIN:2")
	Variable i, m,Nm=ItemsInList(tlist)
	for (m=0;m<Nm;m+=1)						// loop over all displayed tables
		table = StringFromList(m,tlist)				// table name
		i = 0
		do											// check each wave in this table
			Wave wav = WaveRefIndexed(table,i,1)	// 1 means table data
			if (!WaveExists(wav))					// no more waves in this table
				break
			endif
			if (stringmatch(GetWavesDataFolder(wav,2),name0))
				return table						// found our wave, return table name
			endif
			i += 1
		while(1)
	endfor
	return ""										// no exising displayed table found
End


// this routine moved to Utility_JZT.ipf
//// return a list of waves in current folder having all tags in 'require', and none of the tags in 'avoid'
//// This is similar to WaveList(), but with a finer selection
//Function/T WaveListClass(waveClassList,search,options)
//	String waveClassList				// a list of acceptable wave classes (semi-colon separated)
//	String search						// same as first argument in WaveList()
//	String options						// same as last argument in WaveList()
//
//	String in = WaveList(search,";",options), out=""
//	String name,key, wnote
//	Variable m
//	for (m=0, name=StringFromList(0,in); strlen(name); m+=1,name=StringFromList(m,in))
//		wnote = note($name)
//		if (WhichListItem(StringByKey("waveClass",wnote,"="),waveClassList)<0)
//			continue
//		endif
//		out += name+";"
//	endfor
//	return out
//End
//
// return a list of waves in current folder having all tags in 'require', and none of the tags in 'avoid'
// This is similar to WaveList(), but with a finer selection
Function/T WaveList_Tags(require,avoid,search,options)
	String require						// list of required tags in wave note
	String avoid							// list of tags to avoid in wave note
	String search						// same as first argument in WaveList()
	String options						// same as last argument in WaveList()

	String in = WaveList(search,";",options), out=""
	String name,key, wnote
	Variable m,i,Na=ItemsInList(avoid),Nr=ItemsInList(require), bad
	for (m=0, name=StringFromList(0,in); strlen(name); m+=1,name=StringFromList(m,in))
		wnote = note($name)
		bad = 0

		for (i=0;i<Nr;i+=1)				// check for required keys
			key = StringFromList(i,require)
			if (strlen(StringByKey(key,wnote,"="))==0)
				bad = 1						// a required key not found, set bad and break
				break
			endif
		endfor

		for (i=0;i<Na;i+=1)				// check for keys to avoid
			key = StringFromList(i,avoid)
			if (strlen(StringByKey(key,wnote,"=")))
				bad = 1						// found a key to avoid, so stop and set bad flag
				break
			endif
		endfor
		if (!bad)
			out += name+";"
		endif
	endfor
	return out
End



// =========================================================================
// =========================================================================
//	Start of Index strain refinement

#ifdef USE_ENERGY_STRAIN_REFINE
Function/T TotalStrainRefine(pattern,constrain,[coords,FullPeakIndexed,FullPeakIndexed1,FullPeakIndexed2,hklEmeasured,printIt])
	Variable pattern									// pattern number, usually 0
	String constrain									// constraint on optimization, "111111", a 1 is refine, a 0 is keep constant
	Variable coords									// coordinate system to pass in return {0=crystal, 1=BL, 2=XHF, 3=Sample (outward surface normal)}
	Wave FullPeakIndexed
	Wave FullPeakIndexed1, FullPeakIndexed2	// from other (Yellow and Purple) detectors)
	Wave hklEmeasured									// measured energies, contains h,k,l,keV, any following columns such as px,py,detectorNum are ignored
	Variable printIt

	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt
	coords = ParamIsDefault(coords) || numtype(coords) ? 1 : round(coords)		// Beam Line system is the default
	if (coords<0 || coords>3)									// coords must be 0, 1, 2, or 3
		return ""
	endif
	if (!ParamIsDefault(FullPeakIndexed) && !WaveExists(FullPeakIndexed))
		return ""
	endif
	if (!ParamIsDefault(FullPeakIndexed) && !WaveExists(FullPeakIndexed))
		return ""
	endif
	if (!ParamIsDefault(FullPeakIndexed1) && !WaveExists(FullPeakIndexed1))
		return ""
	endif
	if (!ParamIsDefault(FullPeakIndexed2) && !WaveExists(FullPeakIndexed2))
		return ""
	endif
	if (!ParamIsDefault(hklEmeasured) && !WaveExists(hklEmeasured))
		return ""
	endif
	String indexWaveList = WaveListClass("IndexedPeakList*","*","")
	String hklElist=WaveListClass("measuredEnergy*","*","")
	if (!WaveExists(FullPeakIndexed))
		if (ItemsInList(indexWaveList)<1)
			return ""
		elseif (ItemsInList(indexWaveList)==1)
			Wave FullPeakIndexed = $StringFromList(0,indexWaveList)
		endif
	endif
	if (!WaveExists(FullPeakIndexed) || !WaveExists(hklEmeasured) && strlen(hklElist))
		String peaksName, indexName="",indexName1="",indexName2="", hklEname=""
		//coords += 1				// change from 0 based to 1 based
		//Prompt coords,"strain coordinates",popup,"Crystal;Beam Line;XHF;Sample (outward normal)"
		Prompt indexName,"index list",popup,indexWaveList
		Prompt indexName1,"index list 2nd detector",popup,"_none_;"+indexWaveList
		Prompt indexName2,"index list 3rd detector",popup,"_none_;"+indexWaveList
		Prompt hklEname,"measured energies",popup,hklElist+"_none_;"
		if (ItemsInList(hklElist)>0)
			DoPrompt "pick index list",indexName,indexName1,indexName2,hklEname
		else
			DoPrompt "pick index list",indexName,indexName1,indexName2
		endif
		if (V_flag)
			return ""
		endif
		//coords -= 1												// change back to 0 based from 1 based
		Wave FullPeakIndexed = $indexName
		Wave FullPeakIndexed1 = $indexName1
		Wave FullPeakIndexed2 = $indexName2
		Wave hklEmeasured = $hklEname
	endif
	if (!WaveExists(FullPeakIndexed))
		return ""
	endif
	STRUCT microGeometry geo
	if (FillGeometryStructDefault(geo, alert=1))		//fill the geometry structure with test values
		return ""
	endif
	Variable Npatterns=DimSize(FullPeakIndexed,2)
	if (WaveExists(FullPeakIndexed1))
		Npatterns = min(Npatterns,DimSize(FullPeakIndexed1,2))
	endif
	if (WaveExists(FullPeakIndexed2))
		Npatterns = min(Npatterns,DimSize(FullPeakIndexed2,2))
	endif
	Variable deviatoric=1
	if (WaveExists(hklEmeasured))
		Make/N=(DimSize(hklEmeasured,0))/FREE energes
		WaveStats/M=1/Q energes
		deviatoric = V_npnts<1									// use deviatoric=TRUE when no energies
	endif
	deviatoric = !(!deviatoric)
	Variable aFit,bFit,cFit,alphaFit,betFit,gamFit=NaN
	sscanf constrain, "%1d%1d%1d%1d%1d%1d",aFit,bFit,cFit,alphaFit,betFit,gamFit
	if (V_flag!=6 || !((aFit+bFit+cFit+alphaFit+betFit+gamFit)>0))
		aFit = NaN													// flags bad input
	endif
	if (deviatoric && str2num(constrain[0,2])>=111)
		aFit = NaN													// when deviatoric, cannot fit all three lengths
	endif
	if (!(pattern>=0) || numtype(aFit))					// need to ask user
		pattern = limit(pattern,0,Npatterns-1)
		if (StartStrainPanel(pattern,Npatterns,aFit,bFit,cFit,alphaFit,betFit,gamFit,deviatoric=deviatoric)<0)
			return ""
		endif
		printIt = 1
	endif
	aFit = aFit ? 1 : 0 ;				bFit = bFit ? 1 : 0 ;			cFit = cFit ? 1 : 0		// ensure only 0 or 1
	alphaFit = alphaFit ? 1 : 0 ;	betFit = betFit ? 1 : 0 ;		gamFit = gamFit ? 1 : 0
	Variable Nfit = aFit+bFit+cFit+alphaFit+betFit+gamFit
	if (Nfit<1 || (aFit&&bFit&&cFit&&deviatoric))
		return ""													// this is deviatoric, you cannot fit all 3 lengths
	endif
	sprintf constrain "%d%d%d%d%d%d", aFit,bFit,cFit,alphaFit,betFit,gamFit		// constrain is for {a,b,c,alpha,beta,gamma}

	Make/WAVE/FREE/N=3 indexWaves,Qwaves
	indexWaves[0] = FullPeakIndexed
	indexWaves[1] = FullPeakIndexed1
	indexWaves[2] = FullPeakIndexed2
	Make/N=(MAX_Ndetectors)/FREE Ns=NaN, dNums=NaN
	Ns=DimSize(indexWaves[p],0)
	Ns = numtype(Ns) ? 0 : Ns
	Variable N = sum(Ns)										// total number of indexed peaks to use for fit

	Variable idetector
	for (idetector=0;idetector<3;idetector+=1)
		if (WaveExists(indexWaves[idetector]))
			Wave FullPeakList = $StringByKey("peakListWave",note(indexWaves[idetector]),"=")
			Wave Qs = FullPeakList2Qwave(FullPeakList)	// measured qhats
			if (!WaveExists(Qs))
				return ""
			endif
			dNums[idetector] = max(detectorNumFromID(geo, StringByKey("detectorID", note(FullPeakList),"=")),0)
			Qwaves[idetector] = Qs
		endif
	endfor
	String funcName = "Indexing#latticeMismatch"+num2istr(Nfit+3)	// need 3 extra for the rotation, which is always there
	KillWaves/Z epsilon,epsilonBL,epsilonXHF,epsilonSample
	//
	// done with I/O, now setup the calculation

	// from initial recip lattice, find lattice constants, and initial rotation (stored as rx,ry,rz)
	Variable as0,as1,as2,bs0,bs1,bs2,cs0,cs1,cs2
	Make/N=(3,3)/FREE/D RLmeas
	sscanf StringByKey("recip_lattice"+num2istr(pattern),note(FullPeakIndexed),"="), "{{%g,%g,%g}{%g,%g,%g}{%g,%g,%g}}",as0,as1,as2,bs0,bs1,bs2,cs0,cs1,cs2
	RLmeas[0][0] = as0 ;		RLmeas[0][1] = bs0 ;		RLmeas[0][2] = cs0		// the original measured RL
	RLmeas[1][0] = as1 ;		RLmeas[1][1] = bs1 ;		RLmeas[1][2] = cs1
	RLmeas[2][0] = as2 ;		RLmeas[2][1] = bs2 ;		RLmeas[2][2] = cs2
	Make/N=6/O/D optimize_LatticeConstantsWave
	Wave LC=optimize_LatticeConstantsWave
	Variable Vc = RL2latticeConstants(RLmeas,LC)		// calc original lattice constants from RL, stored in LC, save original Vc too
	forceLattice(LC,NumberByKey("SpaceGroup",note(FullPeakIndexed),"="),Vc=Vc)	// force lattice constants to match Space Group exactly, preserve Vc

	Make/N=(3,3)/D/FREE DL, RL0
	Make/N=(3,3)/O/D RL_latticeMismatch, rho_latticeMismatch
	Wave RL=RL_latticeMismatch, rho=rho_latticeMismatch
	RLfromLatticeConstants(LC,DL,RL0)						// make reference RL from lattice constants (used to compute rho), this RL exactly matches Space Group
	MatrixOp/O rho = RLmeas x Inv(RL0)						// RLmeas = rho x RL0,  the rotation (or almost a perfect rotation matrix)
	Make/N=3/O/D/FREE axis
	Variable angle = axisOfMatrix(rho,axis,squareUp=1)	// rho is almost a  perfect rotation matrix, by remaking it, it will be a perfect rotation matrix
	if (numtype(angle))
		return ""
	endif
	axis *= angle													// length of axis is now rotation angle (rad)
	rotationMatAboutAxis(axis,angle,rho)					// forces rho to be a perfect rotation matrix
	MatrixOp/O/FREE RLmeas = rho x RL0						// ensure that starting point perfectly agrees with Space Group

	Make/N=(3,3)/FREE/D Ameas, Astart						// Vcartesian = A x Vcell,   used to make the strain tensor
	Astart = DL														// the starting direct lattice
	if (printIt)
		if (stringmatch(StringFromList(0,GetRTStackInfo(0)),"IndexButtonProc"))
			printf "•"
		endif
		printf "TotalStrainRefine(%d,\"%d%d%d%d%d%d\"",pattern,aFit,bFit,cFit,alphaFit,betFit,gamFit
		if (!ParamIsDefault(coords))
			printf ",coords=%g",coords
		endif
		if (WaveExists(indexWaves[0]))
			printf ", FullPeakIndexed=%s",NameOfWave(indexWaves[0])
		endif
		if (WaveExists(indexWaves[1]))
			printf ", FullPeakIndexed1=%s",NameOfWave(indexWaves[1])
		endif
		if (WaveExists(indexWaves[2]))
			printf ", FullPeakIndexed2=%s",NameOfWave(indexWaves[2])
		endif
		if (WaveExists(hklEmeasured))
			printf ", hklEmeasured=%s",NameOfWave(hklEmeasured)
		endif
		printf ")\r"
		if (Npatterns>1)
			printf "For pattern  %d  of  [0, %d]\r",pattern,Npatterns-1
		endif
		printf "starting reciprocal lattice,	a* = {%+.6f, %+.6f, %+.6f} (1/nm)\r",as0,as1,as2
		printf "							b* = {%+.6f, %+.6f, %+.6f}\r",bs0,bs1,bs2
		printf "							c* = {%+.6f, %+.6f, %+.6f}\r",cs0,cs1,cs2
		print " "
		printf "a = %.7g nm,   b = %.7g,   c = %.7g,   alpha = %.8g%s,   beta = %.8g%s,   gamma = %.8g%s,   Vc = %.7g (nm^3)\r",LC[0],LC[1],LC[2],LC[3],DEGREESIGN,LC[4],DEGREESIGN,LC[5],DEGREESIGN,Vc
		printf "a/c = %.7g,   b/c = %.7g,   a/b = %.7g\r"LC[0]/LC[2],LC[1]/LC[2],LC[0]/LC[1]
	endif

	String measuredType = ""									// default to no extra measured information (keV, Q, or d)
	if (WaveExists(hklEmeasured))							// get proper column label for measured keV, Q, or d
		measuredType = GetDimLabel(hklEmeasured,1,3)
		measuredType = SelectString(strlen(measuredType),"keV_measured",measuredType)
		Make/N=(DimSize(hklEmeasured,0))/FREE measured
		measured = hklEmeasured[p][3]
		WaveStats/M=1/Q measured
		if (strsearch(measuredType,"keV",0,2)==0)
		elseif (strsearch(measuredType,"d",0,2)==0)
		elseif (strsearch(measuredType,"Q",0,2)==0)
		else
			V_npnts = -1											// removes use of measured energy, Q, or d
		endif
		if (V_npnts<=0)											// actually nothing extra was measured!
			Wave hklEmeasured=$""	
			measuredType = ""
		endif
	endif
	measuredType = SelectString(strlen(measuredType),"E_NOT_meas",measuredType)	// default to no extra measured information (keV, Q, or d)
	Make/N=(N,17)/O/D PeaksForStrain=NaN
	SetDimLabel 1,0,h,PeaksForStrain ;				SetDimLabel 1,1,k,PeaksForStrain ;				SetDimLabel 1,2,l,PeaksForStrain
	SetDimLabel 1,3,measQx,PeaksForStrain ;		SetDimLabel 1,4,measQy,PeaksForStrain ;		SetDimLabel 1,5,measQz,PeaksForStrain
	SetDimLabel 1,6,QxCalc,PeaksForStrain ;	SetDimLabel 1,7,QyCalc,PeaksForStrain ;	SetDimLabel 1,8,QzCalc,PeaksForStrain
	SetDimLabel 1,9,Intens_meas,PeaksForStrain;SetDimLabel 1,10,keV_calc,PeaksForStrain
	SetDimLabel 1,11,pixX_calc,PeaksForStrain;	SetDimLabel 1,12,pixY_calc,PeaksForStrain;	SetDimLabel 1,13,err_deg,PeaksForStrain
	String measuredTypeColumnStr = SelectString(strlen(measuredType),"none",measuredType)
	SetDimLabel 1,14,detNum,PeaksForStrain;	SetDimLabel 1,15,$measuredTypeColumnStr,PeaksForStrain;	SetDimLabel 1,16,deV,PeaksForStrain
	String wnote = ReplaceNumberByKey("patternNum","",pattern,"=")
	wnote = ReplaceStringByKey("constrain",wnote,constrain,"=")
	wnote = ReplaceStringByKey("waveClass",wnote,"StrainPeakList","=")
	if (strlen(measuredType))
		wnote = ReplaceStringByKey("measuredType",wnote,measuredType,"=")
	endif
	Note/K PeaksForStrain,wnote

	Make/N=3/D/FREE qhat,qcalc, qmeas, hkl
	Variable i,j, m, dNum
	Variable d, sine, keV, qLen
	for (idetector=0,m=0; idetector<MAX_Ndetectors; idetector+=1)
		dNum = dNums[idetector]
		Wave FullPeakIndexed = indexWaves[idetector]
		Wave Qs = Qwaves[idetector]
		if (!WaveExists(Qs))
			continue
		endif
		Make/N=(DimSize(Qs,0),3)/FREE/D/O qqs				// used to find which measured Q^ is closest to an indexed spot
		for (i=0;i<Ns[idetector];i+=1)						// for each indexed peak on this detector
			hkl = FullPeakIndexed[i][p+3][pattern]		// set the hkl
			MatrixOp/O/FREE qcalc = RLmeas x hkl			// predicted g^ from each measured hkl
			qhat = qcalc
			normalize(qhat)
			qqs = Qs[p][q]*qhat[q]
			MatrixOP/FREE/O dots = sumRows(qqs)
			WaveStats/M=1/Q dots
			if (V_max>0.99)										// require that Q's are parallel to within cos(0.99)=8.11°
				PeaksForStrain[m][0,2] = hkl[q]				// set the hkl
				qmeas = Qs[V_maxloc][p]						// closest measured peak to this hkl
				d = dSpacingFromLatticeConstants(hkl[0],hkl[1],hkl[2],LC[0],LC[1],LC[2],LC[3],LC[4],LC[5])
				qmeas = Qs[V_maxloc][p]*2*PI/d				// closest measured peak to this hkl
				PeaksForStrain[m][3,5] = qmeas[q-3]		// store measured q (use a calulated energy for length of Q)
				PeaksForStrain[m][6,8] = qcalc[q-6]		// predicted Qvectors (1/nm)
				PeaksForStrain[m][9] = Qs[V_maxloc][3]	// measured intensity
				PeaksForStrain[m][14] = dNum
				m += 1
			endif
		endfor
	endfor
	N = m
	Redimension/N=(N,-1) PeaksForStrain
	if ((Nfit+3+(!deviatoric))>(2*N))						// there are (Nfit+3) free parameters (Nfit lattice constants + 3 rotations)
		return ""													// each measured spot provides 2 values
	endif

	if (WaveExists(hklEmeasured))							// measured energies exist, set them
		Make/N=(N,3)/D/FREE dhkls
		Make/N=3/D/FREE dhkl
		for (j=0;j<DimSize(hklEmeasured,0);j+=1)		// loop over each measured value, find matching hkl
			hkl = hklEmeasured[j][p]
			dhkls = PeaksForStrain[p][q] - hkl[q]
			MatrixOP/FREE/O deltas = sumRows(dhkls*dhkls)
			WaveStats/Q/M=1 deltas
			if (V_min<1e-5)
				m = V_minloc
				qhat = PeaksForStrain[m][p+3]				// measured Qvector
				normalize(qhat)
				if (strsearch(measuredType,"keV",0,2)==0)
					sine = -qhat[2]								// measured sin(Bragg angle)
					keV = hklEmeasured[j][3]					// measured energy
					qLen = 4*PI*sine*keV/hc					// Q = 4π sin(theta)/lambda
				elseif (strsearch(measuredType,"d",0,2)==0)
					qLen = 2*PI/hklEmeasured[j][3]			// passing d-spacing,  Q=2π/d
				elseif (strsearch(measuredType,"Q",0,2)==0)
					qLen = hklEmeasured[j][3]					// directly passing Q
				endif
				qmeas = qhat * qLen
				PeaksForStrain[m][3,5] = qmeas[q-3]			// re-assign measured Q using measured energy for |Q|
				PeaksForStrain[m][15] = hklEmeasured[j][3]	// save measured energy
			endif
		endfor
	endif

	Variable rmsErr = latticeMismatchAll(PeaksForStrain,axis[0],axis[1],axis[2],LC)
	if (printIt)
		printf "at start, rms error = %.5f%s",rmsErr*180/PI,DEGREESIGN
		printf " \t\t****** %s ******\r",SelectString(deviatoric,"Total Strain Refinement, Vc can vary","Deviatoric Strain Refinement, Vc is fixed")
	endif

	// set up for Optimization
	Make/N=(Nfit+3)/FREE/D typXWave, xWave
	typXWave[0,2] = axis[p]
	for (m=3,i=0;i<6 && m<(Nfit+3);i+=1)
		if (str2num(constrain[i]))
			typXWave[m] = LC[i]
			m += 1
		endif
	endfor
	xWave = typXWave

	Variable/G printOnlyFirstStrainNaN=1
	Optimize/Q/A=0/M={0,0}/R=typXWave/X=xWave/I=1000/Y=(rmsErr) $funcName, PeaksForStrain
	if (V_flag==791 && V_OptTermCode==6)					// too close to a critical point, change starting point and try again
		Variable del=1.001
		xWave[0]*=del
		if (numpnts(xWave)>3)
			xWave[3]*=del
			xWave[6]*=del
		endif
		printOnlyFirstStrainNaN=1
		Optimize/Q/A=0/M={0,0}/R=typXWave/X=xWave/I=1000/Y=(rmsErr) $funcName, PeaksForStrain
	elseif (V_flag==789 ||  V_OptTermCode==5)				// optimization got lost, try a different method
		printf "firtst try at optimization using a 'LIne Search' failed with V_flag=%g, V_OptTermCode=%g,   try again using 'More-Hebdon'\r",V_flag,V_OptTermCode
		xWave = typXWave
		printOnlyFirstStrainNaN=1
		Optimize/Q/A=0/M={2,0}/R=typXWave/X=xWave/I=1000/Y=(rmsErr) $funcName, PeaksForStrain
	endif
	KillVariables/Z printOnlyFirstStrainNaN

	if (V_flag)														// Optimization failed
		if (printIt)
			printf "Optimize failed with V_flag=%d,   V_OptTermCode=%d,   V_OptNumIters=%d,   V_OptNumFunctionCalls=%d\r",V_flag,V_OptTermCode,V_OptNumIters,V_OptNumFunctionCalls
		endif
		KillWaves/Z W_OptGradient, W_Cross
		KillWaves/Z rho_latticeMismatch
		return ""
	endif

	// Optimization finished, compute the lattice constants and all of the strain tensors

	axis = xWave[p]												// axis = {rx,ry,rz}
	fillLC(LC,constrain,xWave[3],xWave[4],xWave[5],xWave[6],xWave[7],xWave[8])
	KillWaves/Z W_OptGradient
	if (deviatoric)
		forceLattice(LC,1,Vc=Vc)								// re-adjusts a,b,c to give correct Vc, (SpaceGroup=1 is triclinic, so no readjusting is done)
	endif
	rmsErr = latticeMismatchAll(PeaksForStrain,axis[0],axis[1],axis[2],LC)		// 3 rotations and all 6 lattice constants

	angle= norm(axis)
	rotationMatAboutAxis(axis,angle,rho)					// calculate rho, the rotation matrix from {rx,ry,rz}
	if (deviatoric)
		RLfromLatticeConstants(LC,DL,RL,Vc=Vc)			// calculate the reciprocal lattice from {a,b,c,alpha,bet,gam}
	else
		Vc=RLfromLatticeConstants(LC,DL,RL)				// calculate the reciprocal lattice from {a,b,c,alpha,bet,gam}
	endif
	MatrixOp/O RL = rho x RL									// rotate the reciprocal lattice by rho, this is the strained RL
	MatrixOp/O/FREE DL = rho x DL
	Ameas = DL														// Vcartesian = A x Vcell,   used to make the strain tensor

	MatrixOp/FREE F = Ameas x Inv(Astart)					// transformed un-strained to strained,  Ameas = F x Ao,        F = R x U
	MatrixOp/FREE U = F^t x F									// U is symmetric, so F^t x F = U^t x R^t x R x U = U x R^-1 x R x U = U x U = U^2
	if (sqrtSymmetricMat(U))									// U starts as F^t x F, and is replaced by sqrt(U^2)
		KillWaves/Z rho_latticeMismatch
		return ""
	endif
	MatrixOp/O epsilon = U - Identity(3)					// U = I + epsilon

	Variable trace = MatrixTrace(epsilon)
	if (deviatoric)
		epsilon -= (p==q)*trace/3								// only deviatoric part, subtract trace/3 from diagonal
	endif
	trace = MatrixTrace(epsilon)								// re-set trace to be just the epsilon part, should be zero
	trace = abs(trace)<1e-15 ? 0 : trace					// set really small numbers to zero
	Make/N=(3,3)/D/FREE epsilonAbs
	epsilonAbs = abs(epsilon)
	WaveStats/Q epsilonAbs
	Variable SumAbsEpsilon=V_Sum
	Variable epsilonvM											// von Mises strain
	epsilonvM = (epsilon[0][0]-epsilon[1][1])^2 + (epsilon[1][1]-epsilon[2][2])^2 + (epsilon[2][2]-epsilon[0][0])^2 
	epsilonvM += 6*( epsilon[0][1]^2 + epsilon[1][2]^2 + epsilon[2][0]^2 )
	epsilonvM = sqrt(epsilonvM/2)

	MatrixOp/O epsilonBL = rho x epsilon x Inv(rho)	// epsilon in beam line coordinates

	Make/N=(3,3)/D/FREE Rlocal
	Rlocal = 0
	Rlocal[0][0] = 1
	Rlocal[1,2][1,2] = 1/sqrt(2)
	Rlocal[2][1] *= -1
	MatrixOp/O epsilonXHF = Rlocal x epsilonBL x Inv(Rlocal)	// get epsilon in the XHF coordinate system
	Rlocal[1,2][1,2] = -1/sqrt(2)							// now Rlocal transforms BL -->  "sample system" (uses outward surface normal for sample at 45°)
	Rlocal[2][1] *= -1
	MatrixOp/O epsilonSample = Rlocal x epsilonBL x Inv(Rlocal)// get epsilon in the Sample coordinate system

	Make/N=3/FREE/D fithat
	Variable startx,groupx, starty,groupy, detNum, keVmeasured
	Variable/C pz
	for (i=0;i<N;i+=1)											// for each peak on the detector
		qcalc = PeaksForStrain[i][p+6]						// strained peak directions
		qmeas = PeaksForStrain[i][p+3]						// measured peak directions
		MatrixOP/FREE/O delta = sqrt(sumSqr(qcalc-qmeas))
		PeaksForStrain[i][13] = delta[0]/norm(qcalc)*180/PI
		sine = -qcalc[2]/norm(qcalc)							// sin(Bragg angle)
		d = dSpacingFromLatticeConstants(PeaksForStrain[i][0],PeaksForStrain[i][1],PeaksForStrain[i][2],LC[0],LC[1],LC[2],LC[3],LC[4],LC[5])
		keV = hc/(2*d*sine)
		PeaksForStrain[i][10] = keV							// calculated energy of strained reflection (at strained position, not fitted position)

		if (strsearch(measuredType,"keV",0,2)==0)		// get the measured energy to calculate dE
			keVmeasured = PeaksForStrain[i][15]
		elseif (strsearch(measuredType,"d",0,2)==0)
			d = PeaksForStrain[i][15]							// d (nm)
			keVmeasured = hc / (2*d*sine)					// lambda = 2 d sin(theta)
		elseif (strsearch(measuredType,"Q",0,2)==0)
			qLen = PeaksForStrain[i][15]						// Q (1/nm)
			keVmeasured = qLen*hc/(4*PI*sine)				// Q = 4π sin(theta)/lambda
		else
			keVmeasured = NaN										// unknown measuredType
		endif

		PeaksForStrain[i][16] = (keV-keVmeasured)*1000// energy error
		detNum = PeaksForStrain[i][14]
		pz = q2pixel(geo.d[detNum],qcalc)
		idetector = BinarySearch(dNums,detNum)
		Wave FullPeakIndexed = indexWaves[idetector]
		if (WaveExists(FullPeakIndexed))
			wnote = note(FullPeakIndexed)
			startx = NumberByKey("startx",wnote,"=")
			groupx = NumberByKey("groupx",wnote,"=")
			starty = NumberByKey("starty",wnote,"=")
			groupy = NumberByKey("groupy",wnote,"=")
			startx = numtype(startx) ? FIRST_PIXEL : startx
			groupx = numtype(groupx) ? 1 : groupx
			starty = numtype(starty) ? FIRST_PIXEL : starty
			groupy = numtype(groupy) ? 1 : groupy
			PeaksForStrain[i][11] = (real(pz)-(startx-FIRST_PIXEL)-(groupx-1)/2)/groupx	// change to binned pixels
			PeaksForStrain[i][12] = (imag(pz)-(starty-FIRST_PIXEL)-(groupy-1)/2)/groupy	// pixels are still zero based
		endif
	endfor

	Variable a=LC[0], b=LC[1], c=LC[2], alpha=LC[3], bet=LC[4], gam=LC[5]
	String str, star1 = SelectString(deviatoric,"","*****  "), star2 = SelectString(deviatoric,"","  *****")
	if (printIt)
		print " "
		printf "final reciprocal lattice,	a* = {%+.6f, %+.6f, %+.6f} (1/nm)\r",RL[0][0],RL[1][0],RL[2][0]
		printf "							b* = {%+.6f, %+.6f, %+.6f}\r",RL[0][1],RL[1][1],RL[2][1]
		printf "							c* = {%+.6f, %+.6f, %+.6f}\r",RL[0][2],RL[1][2],RL[2][2]
		printf "final lattice,	a = {%+.6f, %+.6f, %+.6f} (nm)\r",DL[0][0],DL[1][0],DL[2][0]
		printf "				b = {%+.6f, %+.6f, %+.6f}\r",DL[0][1],DL[1][1],DL[2][1]
		printf "				c = {%+.6f, %+.6f, %+.6f}\r",DL[0][2],DL[1][2],DL[2][2]
		print " "
		printf "a = %.7g nm,   b = %.7g,   c = %.7g,   alpha = %.8g%s,   beta = %.8g%s,   gamma = %.8g%s,   Vc = %.7g (nm^3)\r",a,b,c,alpha,DEGREESIGN,bet,DEGREESIGN,gam,DEGREESIGN,Vc
		printf "a/c = %.7g,   b/c = %.7g,   a/b = %.7g\r",a/c,b/c,a/b
		printf "at end, rms error    = %.5f%s\r",rmsErr*180/PI,DEGREESIGN
		print " "
		sprintf str,"  \t%sTr(epsilon) = %.2g%s",star1,trace,star2
		printf "epsilonCrystal =	{%+.6f, %+.6f, %+.6f}  \tvon Mises strain = %.3g%s\r",epsilon[0][0],epsilon[0][1],epsilon[0][2], epsilonvM, SelectString(abs(trace)>1e-12,"",str)
		printf "					{%+.6f, %+.6f, %+.6f}  \tmax{ | epsilon(ij) | } = %.3g\r",epsilon[1][0],epsilon[1][1],epsilon[1][2], V_max
		printf "					{%+.6f, %+.6f, %+.6f}  \tSum{ | epsilon | } = %.3g\r\r",epsilon[2][0],epsilon[2][1],epsilon[2][2], SumAbsEpsilon
		WaveStats/Q epsilonBL
		sprintf str,"\t%sTr(epsilonBL) = %.2g%s",star1,MatrixTrace(epsilonBL),star2
		printf "epsilonBL =		{%+.6f, %+.6f, %+.6f}  \tmax{ | epsilon(ij) | } = %.3g\r",epsilonBL[0][0],epsilonBL[0][1],epsilonBL[0][2],max(abs(V_max),abs(V_min))
		printf "					{%+.6f, %+.6f, %+.6f}%s\r",epsilonBL[1][0],epsilonBL[1][1],epsilonBL[1][2], SelectString(abs(MatrixTrace(epsilonBL))>1e-12,"",str)
		printf "					{%+.6f, %+.6f, %+.6f}\r\r",epsilonBL[2][0],epsilonBL[2][1],epsilonBL[2][2]
		WaveStats/Q epsilonXHF
		sprintf str,"  \t%sTr(epsilonXHF) = %.2g%s",star1,MatrixTrace(epsilonXHF),star2
		printf "epsilonXHF =		{%+.6f, %+.6f, %+.6f}  \tmax{ | epsilon(ij) | } = %.3g\r",epsilonXHF[0][0],epsilonXHF[0][1],epsilonXHF[0][2],max(abs(V_max),abs(V_min))
		printf "					{%+.6f, %+.6f, %+.6f}%s\r",epsilonXHF[1][0],epsilonXHF[1][1],epsilonXHF[1][2], SelectString(abs(MatrixTrace(epsilonXHF))>1e-12,"",str)
		printf "					{%+.6f, %+.6f, %+.6f}\r\r",epsilonXHF[2][0],epsilonXHF[2][1],epsilonXHF[2][2]
		WaveStats/Q epsilonSample
		sprintf str,"  \t%sTr(epsilonSample) = %.2g%s",star1,MatrixTrace(epsilonSample),star2
		printf "epsilonSample =	{%+.6f, %+.6f, %+.6f}  \tmax{ | epsilon(ij) | } = %.3g\r",epsilonSample[0][0],epsilonSample[0][1],epsilonSample[0][2],max(abs(V_max),abs(V_min))
		printf "					{%+.6f, %+.6f, %+.6f}%s\r",epsilonSample[1][0],epsilonSample[1][1],epsilonSample[1][2], SelectString(abs(MatrixTrace(epsilonSample))>1e-12,"",str)
		printf "					{%+.6f, %+.6f, %+.6f}\r\r",epsilonSample[2][0],epsilonSample[2][1],epsilonSample[2][2]
	endif

	wnote = note(indexWaves[0])
	sprintf str, "{{%.8g,%.8g,%.8g}{%.8g,%.8g,%.8g}{%.8g,%.8g,%.8g}}",RL[0][0],RL[1][0],RL[2][0],  RL[0][1],RL[1][1],RL[2][1],  RL[0][2],RL[1][2],RL[2][2]
	wnote = ReplaceStringByKey("recip_lattice_refined",wnote,str,"=")
	sprintf str, "{%.8g %.8g %.8g %.8g %.8g %.8g}",a,b,c,alpha,bet,gam
	wnote = ReplaceStringByKey("lattice_constants_refined",wnote,str,"=")
	wnote = ReplaceNumberByKey("SumAbsEpsilon",wnote,SumAbsEpsilon,"=")
	wnote = ReplaceNumberByKey("vonMisesStrain",wnote,epsilonvM,"=")
	Note/K epsilon, ReplaceStringByKey("waveClass",wnote,"strain_tensor_STD","=")
	Note/K epsilonBL, ReplaceStringByKey("waveClass",wnote,"strain_tensor_BL","=")
	Note/K epsilonXHF, ReplaceStringByKey("waveClass",wnote,"strain_tensor_XHF","=")
	Note/K epsilonSample, ReplaceStringByKey("waveClass",wnote,"strain_tensor_Sample","=")

	wnote = note(PeaksForStrain)
	sprintf str, "{{%.8g,%.8g,%.8g}{%.8g,%.8g,%.8g}{%.8g,%.8g,%.8g}}",RL[0][0],RL[1][0],RL[2][0],  RL[0][1],RL[1][1],RL[2][1],  RL[0][2],RL[1][2],RL[2][2]
	wnote = ReplaceStringByKey("recip_lattice_refined",wnote,str,"=")
	sprintf str, "{%.8g %.8g %.8g %.8g %.8g %.8g}",a,b,c,alpha,bet,gam
	wnote = ReplaceStringByKey("lattice_constants_refined",wnote,str,"=")
	Note/K PeaksForStrain, wnote

	KillWaves/Z rho_latticeMismatch, RL_latticeMismatch
	KillWaves/Z  optimize_LatticeConstantsWave
	KillWaves/Z M_Lower,M_Upper,W_LUPermutation,M_x, W_Cross

	if (coords==0)
		str = GetWavesDataFolder(epsilon,2)
	elseif (coords==1)
		str = GetWavesDataFolder(epsilonBL,2)
	elseif (coords==2)
		str = GetWavesDataFolder(epsilonXHF,2)
	elseif (coords==3)
		str = GetWavesDataFolder(epsilonSample,2)
	endif
	return str
End
//
// Returns rms error in (radian).  This is new, the old way only worked before we had energies.
// this uses the measured energy if available, also points with energy get weighted more heavily
//	The measured Qvecs columns [3,4,5] use the measured energy when it is available, otherwise a computed energy is used.
Static Function latticeMismatchAll(PeaksForStrain,rx,ry,rz,LC)			// 3 rotations and all 6 lattice constants (uses MatrixOP for speed)
	Wave PeaksForStrain
	Variable rx,ry,rz
	Wave LC																	//	lattice parameters:   LC = {a,b,c,alpha,bet,gam}

	Wave RL=RL_latticeMismatch, LC=optimize_LatticeConstantsWave
	if (!WaveExists(RL) || !WaveExists(LC))
		Abort "some of the waves do not exist in latticeMismatchAll()"
	endif
	Variable Vc = NumberByKey("Vc",note(PeaksForStrain),"=")		// This is to be held constant

	Make/N=3/D/FREE axis = {rx,ry,rz}
	Make/N=(3,3)/D/FREE rho
	Variable angle= norm(axis)
	rotationMatAboutAxis(axis,angle,rho)							// calculate rho, the rotation matrix from {rx,ry,rz}

	RLfromLatticeConstants(LC,$"",RL,Vc=Vc)						// calculate the reciprocal lattice from {a,b,c,alpha,bet,gam}
	MatrixOp/O RL = rho x RL											// rotate the reciprocal lattice by rho

	Variable N = DimSize(PeaksForStrain,0)
	Make/N=(N,3)/D/FREE Qvecs, hkls
	Qvecs = PeaksForStrain[p][q+3]									// measured directions (assumed normalized)
	hkls = PeaksForStrain[p][q]										// hkl of reflections (first 3 columns)
	MatrixOp/O/FREE Gvecs = (RL x hkls^t)^t						// calcualted Gvectors from hkls
	PeaksForStrain[][6,8] = Gvecs[p][q-6]							// save for later checking

	MatrixOP/FREE haveE = greater(col(PeaksForStrain,15),0)	// energy >0 (and not NaN)
	MatrixOP/FREE haveQ = greater(sumRows(magSqr(Qvecs)), 0)

	Variable weightE = N/max(sum(haveE),1)						// enhancemet for weighting points with measured energy
	weightE = min(weightE,6)											// max enhancement factor for an energy point
	MatrixOP/FREE weights = haveE*weightE + haveQ				// weighting for each point
	Variable sumWeights = sum(weights)								// used to calculate rms, the effective number of values
	weights /= sumWeights												// this makes sum(weights) = 1

	MatrixOP/FREE Qlen = sqrt(sumRows(Gvecs*Gvecs))			// |Q| for each of the Q's (or G's), used to normalize
	MatrixOP/FREE radHats = NormalizeRows(Gvecs)				// radial unit vector for each Q
	MatrixOP/FREE deltas = (Gvecs-Qvecs)/colRepeat(Qlen,3)	// scale or longer Qvectors will count for more
	MatrixOp/FREE radDots = sumRows(deltas * radHats)		// the radial componenet of deltas

	MatrixOP/FREE perpDots2 = sumRows(magSqr(deltas - colRepeat(radDots,3)*radHats))	// length^2 of perp part of deltas
	MatrixOp/FREE radDots2 = magSqr(radDots)					// length^2 of radial part of deltas
	MatrixOP/FREE err2 = sum( (haveE*radDots2 + haveQ*perpDots2) * weights )	// Sum(i){ w_i(delta_i^2) }
	Variable rms = sqrt(err2[0])
	rms = numtype(rms) ? Inf : rms
	return rms																// error in radian
End
//
#else				// else for #ifdef USE_ENERGY_STRAIN_REFINE
//
Function/T DeviatoricStrainRefine(pattern,constrain,[coords,FullPeakList,FullPeakIndexed,printIt])
	Variable pattern											// pattern number, usually 0
	String constrain											// constraint on optimization, "111111", a 1 is refine, a 0 is keep constant
	Variable coords											// coordinate system to pass in return {0=crystal, 1=BL, 2=XHF, 3=Sample (outward surface normal)}
	Wave FullPeakList,FullPeakIndexed
	Variable printIt											// force full print out of results
	printIt = ParamIsDefault(printIt) ? 0 : !(!printIt)
	printIt = numtype(printIt) ? 0 : printIt
	coords = ParamIsDefault(coords) || numtype(coords) ? 1 : round(coords)		// Beam Line system is the default
	if (coords<0 || coords>3)								// coords must be 0, 1, 2, or 3
		return ""
	endif

	if (!ParamIsDefault(FullPeakList))
		if (!WaveExists(FullPeakList))
			return ""
		endif
	endif
	if (!ParamIsDefault(FullPeakIndexed))
		if (!WaveExists(FullPeakIndexed))
			return ""
		endif
	endif

	if (!WaveExists(FullPeakIndexed))
		String indexWaveList = WaveListClass("IndexedPeakList*","*","")
		if (ItemsInList(indexWaveList) < 1)
			return ""
		elseif (ItemsInList(indexWaveList) == 1)
			Wave FullPeakIndexed = $StringFromList(0,indexWaveList)
		else
			String indexName
			Prompt coords,"strain coordinates",popup,"Crystal;Beam Line;XHF;Sample (outward normal)"
			Prompt indexName,"index list",popup,indexWaveList
			coords += 1				// change from 0 based to 1 based
			DoPrompt "Pick index list",indexName, coords
			if (V_flag)
				return ""
			endif
			coords -= 1				// change from 1 based to 0 based
			Wave FullPeakIndexed = $indexName
		endif
	endif
	if (!WaveExists(FullPeakIndexed))
		return ""
	endif
	String peakLists = StringByKey("peakListWave",note(FullPeakIndexed),"=")
	Variable Npeaks = ItemsInList(peakLists,",")
	Make/N=(Npeaks)/FREE/WAVE FullPeakListWaves = $StringFromList(p,peakLists,",")
	Variable iPeak, mm=0
	for (iPeak=0;iPeak<Npeaks;iPeak+=1)
		mm += WaveExists(FullPeakListWaves[iPeak])
	endfor
	if (mm<1)						// make sure some FullPeakList's exist
		return ""
	endif

	STRUCT microGeometry geo
	if (FillGeometryStructDefault(geo, alert=1))	//fill the geometry structure with test values
		return ""
	endif
	Variable Npatterns=DimSize(FullPeakIndexed,2)

	Variable aFit,bFit,cFit,alphaFit,betFit,gamFit=NaN
	sscanf constrain, "%1d%1d%1d%1d%1d%1d",aFit,bFit,cFit,alphaFit,betFit,gamFit
	if (V_flag!=6 || !((aFit+bFit+cFit+alphaFit+betFit+gamFit)>0) || str2num(constrain[0,2])>=111)
		aFit = NaN												// flags bad input
	endif

	if (!(pattern>=0) || numtype(aFit))						// need to ask user
		pattern = limit(pattern,0,Npatterns-1)
		if (StartStrainPanel(pattern,Npatterns,aFit,bFit,cFit,alphaFit,betFit,gamFit)<0)
			return ""
		endif
		printIt = 1
	endif
	aFit = aFit ? 1 : 0 ;				bFit = bFit ? 1 : 0 ;			cFit = cFit ? 1 : 0		// ensure only 0 or 1
	alphaFit = alphaFit ? 1 : 0 ;		betFit = betFit ? 1 : 0 ;		gamFit = gamFit ? 1 : 0
	Variable Nfit = aFit+bFit+cFit+alphaFit+betFit+gamFit
	if (Nfit<1 || (aFit&&bFit&&cFit))
		return ""													// this is deviatoric, you cannot fit all 3 lengths
	endif
	sprintf constrain "%d%d%d%d%d%d", aFit,bFit,cFit,alphaFit,betFit,gamFit		// constrain is for {a,b,c,alpha,beta,gamma}

	Variable N=DimSize(FullPeakIndexed,0)
	KillWaves/Z epsilon,epsilonBL,epsilonXHF,epsilonSample
	if ((Nfit+3)>=(2*N))										// there are (Nfit+3) free parameters (Nfit lattice constants + 3 rotations)
		return ""													// each measured spot provides 2 values
	endif

	Wave Qs = $""													// all possible measured Qs, check all used detectors
	Variable iQs
	for (iPeak=0;iPeak<Npeaks;iPeak+=1)
		Wave FullPeakList = FullPeakListWaves[iPeak]
		if (!WaveExists(FullPeakList))
			continue
		endif
		Wave Qis = FullPeakList2Qwave(FullPeakList)	// Qxyz, intens (4 columns)
		if (!WaveExists(Qis) || DimSize(Qis,0)<1)
			continue
		endif
		if (!WaveExists(Qs))
			Duplicate/FREE Qis, Qs								// create Qs
		else
			iQs = DimSize(Qs,0)
			Redimension/N=(iQs+DimSize(Qis,0),-1) Qs	// extend Qs
			Qs[iQs,Inf][] = Qis[p-iQs][q]
		endif
	endfor
	if (!WaveExists(Qs) || DimSize(Qs,0)<1)
		return ""
	endif
	String funcName = "Indexing#latticeMismatch"+num2istr(Nfit+3)	// need 3 extra for the rotation, which is always there
	//
	// done with I/O, now setup the calculation

	// from initial recip lattice, find lattice constants, and initial rotation (stored as rx,ry,rz)
	Wave RLmeas = decodeMatFromStr(StringByKey("recip_lattice0",note(FullPeakIndexed),"="))
	Duplicate/FREE RLmeas, RLmeasStart
	Make/N=3/D/FREE hkl
	Make/N=3/O/D axis_latticeMismatch
	Wave axis=axis_latticeMismatch

	Make/N=6/O/D optimize_LatticeConstantsWave
	Wave LC=optimize_LatticeConstantsWave
	Variable Vc = RL2latticeConstants(RLmeas,LC)		// calc original lattice constants from RL, stored in LC, save original Vc too
	forceLattice(LC,NumberByKey("SpaceGroup",note(FullPeakIndexed),"="),Vc=Vc)	// force lattice constants to match Space Group exactly, preserve Vc

	Make/N=(3,3)/D/FREE RL0
	Make/N=(3,3)/O/D DL_latticeMismatch, RL_latticeMismatch, rho_latticeMismatch
	Wave DL=DL_latticeMismatch, RL=RL_latticeMismatch, rho=rho_latticeMismatch
	RLfromLatticeConstants(LC,DL,RL0)						// make reference RL from lattice constants (used to compute rho), this RL exactly matches Space Group
	MatrixOp/O rho = RLmeas x Inv(RL0)						// RLmeas = rho x RL0,  the rotation (or almost a perfect rotation matrix)
	Variable angle = axisOfMatrix(rho,axis,squareUp=1)		// rho is almost a  perfect rotation matrix, by remaking it, it will be a perfect rotation matrix
	if (numtype(angle))
		return ""
	endif
	axis *= angle													// length of axis is now rotation angle (rad)
	rotationMatAboutAxis(axis,angle,rho)					// forces rho to be a perfect rotation matrix
	MatrixOp/O RLmeas = rho x RL0							// ensure that starting point perfectly agrees with Space Group

	Make/N=(3,3)/D/FREE Ameas, Astart						// Vcartesian = A x Vcell,   used to make the strain tensor
	Astart = DL														// the starting direct lattice
	if (printIt)
		if (stringmatch(StringFromList(0,GetRTStackInfo(0)),"IndexButtonProc"))
			printf "•"
		endif
		printf "DeviatoricStrainRefine(%d,\"%d%d%d%d%d%d\"",pattern,aFit,bFit,cFit,alphaFit,betFit,gamFit
		if (!ParamIsDefault(coords))
			printf ",coords=%g",coords
		endif
		if (!ParamIsDefault(printit))
			printf ",printit=%g",printit
		endif
		printf ")\r"
		if (Npatterns>1)
			printf "For pattern  %d  of  [0, %d]\r",pattern,Npatterns-1
		endif
		Variable as0,as1,as2,bs0,bs1,bs2,cs0,cs1,cs2
		as0 = RLmeasStart[0][0]	;	bs0 = RLmeasStart[0][1]	;	cs0 = RLmeasStart[0][2]	// the original measured RL
		as1 = RLmeasStart[1][0]	;	bs1 = RLmeasStart[1][1]	;	cs1 = RLmeasStart[1][2]
		as2 = RLmeasStart[2][0]	;	bs2 = RLmeasStart[2][1]	;	cs2 = RLmeasStart[2][2]
		printf "starting reciprocal lattice,	a* = {%+.6f, %+.6f, %+.6f} (1/nm)\r",as0,as1,as2
		printf "							b* = {%+.6f, %+.6f, %+.6f}\r",bs0,bs1,bs2
		printf "							c* = {%+.6f, %+.6f, %+.6f}\r",cs0,cs1,cs2
		print " "
		printf "a = %.7g nm,   b = %.7g,   c = %.7g,   alpha = %.8g%s,   beta = %.8g%s,   gamma = %.8g%s,   Vc = %.7g (nm^3)\r",LC[0],LC[1],LC[2],LC[3],DEGREESIGN,LC[4],DEGREESIGN,LC[5],DEGREESIGN,Vc
		printf "a/c = %.7g,   b/c = %.7g,   a/b = %.7g\r"LC[0]/LC[2],LC[1]/LC[2],LC[0]/LC[1]
	endif

	Make/N=(N,15)/O/D PeaksForStrain=NaN
	SetDimLabel 1,0,h,PeaksForStrain ;				SetDimLabel 1,1,k,PeaksForStrain ;			SetDimLabel 1,2,l,PeaksForStrain
	SetDimLabel 1,3,fitQx,PeaksForStrain ;		SetDimLabel 1,4,fitQy,PeaksForStrain ;	SetDimLabel 1,5,fitQz,PeaksForStrain
	SetDimLabel 1,6,indexQx,PeaksForStrain ;	SetDimLabel 1,7,indexQy,PeaksForStrain ;SetDimLabel 1,8,indexQz,PeaksForStrain
	SetDimLabel 1,9,fit_Intens,PeaksForStrain;	SetDimLabel 1,10,keV,PeaksForStrain
	SetDimLabel 1,11,pixelX,PeaksForStrain;		SetDimLabel 1,12,pixelY,PeaksForStrain;	SetDimLabel 1,13,angErr,PeaksForStrain
	SetDimLabel 1,14,detNum,PeaksForStrain
	String wnote = ReplaceNumberByKey("patternNum","",pattern,"=")
	wnote = ReplaceStringByKey("constrain",wnote,constrain,"=")
	wnote = ReplaceStringByKey("waveClass",wnote,"StrainPeakList","=")
	Note/K PeaksForStrain,wnote

	Variable Nj=DimSize(Qs,0), jmax								// find the measured peak that goes with each indexed peak (in Q^)
	Make/N=(Nj,3)/D/FREE Q3s = Qs[p][q]						// only the Q^ vector, without the intensity
	MatrixOP/FREE Q3s = NormalizeRows(Q3s)

	Variable i,m
	for (i=0,m=0;i<N;i+=1)
		hkl = FullPeakIndexed[i][p+3][pattern]				// set the hkl
		MatrixOp/O/FREE qhat = RLmeas x hkl					// predicted g^ from each measured hkl
		normalize(qhat)

		MatrixOP/FREE dots = Q3s x qhat
		WaveStats/M=1/Q dots
		if (V_npnts>0)
			jmax = V_maxloc
			PeaksForStrain[m][0,2] = hkl[q]						// store the hkl
			PeaksForStrain[m][3,5] = Q3s[jmax][q-3]			// store closest fitted q^
			PeaksForStrain[m][6,8] = qhat[q-6]					// store predicted g^
			PeaksForStrain[m][9] = Qs[jmax][3]					// 	and its intensity
			PeaksForStrain[m][14] = FullPeakIndexed[i][11][pattern]	// detector number
			m += 1
		else
			PeaksForStrain[m][] = NaN								// bad point
		endif
	endfor
	N = m
	Redimension/N=(N,-1) PeaksForStrain						// trim to number of valid points in PeaksForStrain
	Variable rmsErr = latticeMismatchAll(PeaksForStrain,axis[0],axis[1],axis[2],LC)
	if (printIt)
		printf "at start, rms error = %.5f%s\r",rmsErr*180/PI,DEGREESIGN	// 3 rotations and all 6 lattice constants
	endif
	if ((Nfit+3)>=(2*N))											// there are (Nfit+3) free parameters (Nfit lattice constants + 3 rotations),  requires 1 more 
		return ""														// each measured spot provides 2 values
	endif

	Make/N=(Nfit+3)/O/D optimize_typXWave, optimize_xWave
	optimize_typXWave[0,2] = axis[p]
	for (m=3,i=0;i<6 && m<(Nfit+3);i+=1)
		if (str2num(constrain[i]))
			optimize_typXWave[m] = LC[i]
			m += 1
		endif
	endfor
	optimize_xWave = optimize_typXWave

	Variable/G printOnlyFirstStrainNaN=1
//	Optimize/Q/A=0/M={0,0}/R=optimize_typXWave/X=optimize_xWave/I=1000 $funcName, PeaksForStrain
	Optimize/Q/A=0/M={0,0}/R=optimize_typXWave/X=optimize_xWave/I=1000/Y=(rmsErr) $funcName, PeaksForStrain
	if (V_flag==791 && V_OptTermCode==6)				// too close to a critical point, change starting point and try again
		Variable del=1.001
		optimize_xWave[0]*=del
		if (numpnts(optimize_xWave)>3)
			optimize_xWave[3]*=del
			optimize_xWave[6]*=del
		endif
		printOnlyFirstStrainNaN=1
		Optimize/Q/A=0/M={0,0}/R=optimize_typXWave/X=optimize_xWave/I=1000/Y=(rmsErr) $funcName, PeaksForStrain
	elseif (V_flag==789 ||  V_OptTermCode==5)			// optimization got lost, try a different method
		printf "firtst try at optimization using a 'LIne Search' failed with V_flag=%g, V_OptTermCode=%g,   try again using 'More-Hebdon'\r",V_flag,V_OptTermCode
		optimize_xWave = optimize_typXWave
		printOnlyFirstStrainNaN=1
		Optimize/Q/A=0/M={2,0}/R=optimize_typXWave/X=optimize_xWave/I=1000/Y=(rmsErr) $funcName, PeaksForStrain
	endif
	KillVariables/Z printOnlyFirstStrainNaN

	if (V_flag)													// Optimization failed
		if (printIt)
			printf "Optimize failed with V_flag=%d,   V_OptTermCode=%d,   V_OptNumIters=%d,   V_OptNumFunctionCalls=%d\r",V_flag,V_OptTermCode,V_OptNumIters,V_OptNumFunctionCalls
		endif
		KillWaves/Z totalTrans
		KillWaves/Z W_Cross
		KillWaves/Z W_OptGradient, optimize_xWave, optimize_typXWave
		KillWaves/Z axis_latticeMismatch, rho_latticeMismatch
		return ""
	endif

	axis = optimize_xWave[p]								// axis = {rx,ry,rz}
	fillLC(LC,constrain,optimize_xWave[3],optimize_xWave[4],optimize_xWave[5],optimize_xWave[6],optimize_xWave[7],optimize_xWave[8])
	KillWaves/Z W_OptGradient, optimize_xWave, optimize_typXWave
	forceLattice(LC,1,Vc=Vc)								// re-adjusts a,b,c to give correct Vc, (SpaceGroup=1 is triclinic, so no readjusting is done)
	rmsErr = latticeMismatchAll(PeaksForStrain,axis[0],axis[1],axis[2],LC)		// 3 rotations and all 6 lattice constants

	angle= norm(axis)
	rotationMatAboutAxis(axis,angle,rho)				// calculate rho, the rotation matrix from {rx,ry,rz}
	RLfromLatticeConstants(LC,DL,RL,Vc=Vc)			// calculate the reciprocal lattice from {a,b,c,alpha,bet,gam}
	MatrixOp/O RL = rho x RL								// rotate the reciprocal lattice by rho, this is the strained RL
	MatrixOp/O DL = rho x DL
	Ameas = DL												// Vcartesian = A x Vcell,   used to make the strain tensor

	MatrixOp/FREE F = Ameas x Inv(Astart)			// transformed un-strained to strained,  Ameas = F x Ao,  F = R x U
	MatrixOp/FREE U = F^t x F							// U is symmetric, so F^t x F = U^t x R^t x R x U = U x R^-1 x R x U = U x U = U^2
	if (sqrtSymmetricMat(U))								// U starts as F^t x F, and is replaced by sqrt(U^2)
		KillWaves/Z axis_latticeMismatch, rho_latticeMismatch
		return ""
	endif
	MatrixOp/O epsilon = U - Identity(3)				// U = I + epsilon

	Variable trace = MatrixTrace(epsilon)
	epsilon -= (p==q)*trace/3							// only deviatoric part, subtract trace/3 from diagonal
	trace = MatrixTrace(epsilon)						// re-set trace to be just the epsilon part, should be zero
	trace = abs(trace)<1e-15 ? 0 : trace				// set really small numbers to zero
	Wave epsilonAbs = $(microGeo#MakeUnique3x3Mat($""))
	epsilonAbs = abs(epsilon)
	WaveStats/Q epsilonAbs
	Variable SumAbsEpsilon=V_Sum
	KillWaves/Z epsilonAbs
	Variable epsilonvM										// von Mises strain
	epsilonvM = (epsilon[0][0]-epsilon[1][1])^2 + (epsilon[1][1]-epsilon[2][2])^2 + (epsilon[2][2]-epsilon[0][0])^2 
	epsilonvM += 6*( epsilon[0][1]^2 + epsilon[1][2]^2 + epsilon[2][0]^2 )
	epsilonvM = sqrt(epsilonvM/2)

	MatrixOp/O epsilonBL = rho x epsilon x Inv(rho)		// epsilon in beam line coordinates

	Make/N=(3,3)/D/FREE Rlocal							// make a new 3x3 matrix,for the rotation from BL --> XHF coordinates (45° about X)
	Rlocal = 0
	Rlocal[0][0] = 1
	Rlocal[1,2][1,2] = 1/sqrt(2)
	Rlocal[2][1] *= -1
	MatrixOp/O epsilonXHF = Rlocal x epsilonBL x Inv(Rlocal)	// get epsilon in the XHF coordinate system
	Rlocal[1,2][1,2] = -1/sqrt(2)						// now Rlocal transforms BL -->  "sample system" (uses outward surface normal for sample at 45°)
	Rlocal[2][1] *= -1
	MatrixOp/O epsilonSample = Rlocal x epsilonBL x Inv(Rlocal)// get epsilon in the Sample coordinate system

	Variable startx,groupx, starty,groupy
	startx = NumberByKey("startx",note(FullPeakIndexed),"=")
	groupx = NumberByKey("groupx",note(FullPeakIndexed),"=")
	starty = NumberByKey("starty",note(FullPeakIndexed),"=")
	groupy = NumberByKey("groupy",note(FullPeakIndexed),"=")
	startx = numtype(startx) ? FIRST_PIXEL : startx
	groupx = numtype(groupx) ? 1 : groupx
	starty = numtype(starty) ? FIRST_PIXEL : starty
	groupy = numtype(groupy) ? 1 : groupy
	Variable/C pz
	Variable sine, d, keV
	Variable dNum = max(detectorNumFromID(geo, StringByKey("detectorID", wnote,"=")),0)
	Make/N=3/O/D fithat_PeaksForStrain
	Wave fithat = fithat_PeaksForStrain
	for (i=0;i<N;i+=1)
		qhat = PeaksForStrain[i][p+6]						// strained peak directions
		fithat = PeaksForStrain[i][p+3]						// fitted peak directions
		PeaksForStrain[i][13] = acos(limit(MatrixDot(qhat,fithat),-1,1))*180/PI
		sine = -qhat[3]										// sin(Bragg angle)
		d = dSpacingFromLatticeConstants(PeaksForStrain[i][0],PeaksForStrain[i][1],PeaksForStrain[i][2],LC[0],LC[1],LC[2],LC[3],LC[4],LC[5])
		keV = hc/(2*d*sine)
		PeaksForStrain[i][10] = keV				// energy of strained reflection (at strained position, not fitted position)
		if (DimSize(PeaksForStrain,1) > 14)
			dNum = PeaksForStrain[i][14]
		endif
		pz = q2pixel(geo.d[dNum],qhat)
		PeaksForStrain[i][11][pattern] = (real(pz)-(startx-FIRST_PIXEL)-(groupx-1)/2)/groupx		// change to binned pixels
		PeaksForStrain[i][12][pattern] = (imag(pz)-(starty-FIRST_PIXEL)-(groupy-1)/2)/groupy	// pixels are still zero based
	endfor

	Variable a=LC[0], b=LC[1], c=LC[2], alpha=LC[3], bet=LC[4], gam=LC[5]
	String str
	if (printIt)
		print " "
		printf "final reciprocal lattice,	a* = {%+.6f, %+.6f, %+.6f} (1/nm)\r",RL[0][0],RL[1][0],RL[2][0]
		printf "							b* = {%+.6f, %+.6f, %+.6f}\r",RL[0][1],RL[1][1],RL[2][1]
		printf "							c* = {%+.6f, %+.6f, %+.6f}\r",RL[0][2],RL[1][2],RL[2][2]
		printf "final lattice,	a = {%+.6f, %+.6f, %+.6f} (nm)\r",DL[0][0],DL[1][0],DL[2][0]
		printf "				b = {%+.6f, %+.6f, %+.6f}\r",DL[0][1],DL[1][1],DL[2][1]
		printf "				c = {%+.6f, %+.6f, %+.6f}\r",DL[0][2],DL[1][2],DL[2][2]
		print " "
		printf "a = %.7g nm,   b = %.7g,   c = %.7g,   alpha = %.8g%s,   beta = %.8g%s,   gamma = %.8g%s,   Vc = %.7g (nm^3)\r",a,b,c,alpha,DEGREESIGN,bet,DEGREESIGN,gam,DEGREESIGN,Vc
		printf "a/c = %.7g,   b/c = %.7g,   a/b = %.7g\r",a/c,b/c,a/b
		printf "at end, rms error    = %.5f%s\r",rmsErr*180/PI,DEGREESIGN
		print " "
		sprintf str,"  \t*****  Tr(epsilon) = %.2g  *****",trace
		printf "epsilonCrystal =	{%+.6f, %+.6f, %+.6f}  \tvon Mises strain = %.3g%s\r",epsilon[0][0],epsilon[0][1],epsilon[0][2], epsilonvM, SelectString(abs(trace)>1e-12,"",str)
		printf "					{%+.6f, %+.6f, %+.6f}  \tmax{ | epsilon(ij) | } = %.3g\r",epsilon[1][0],epsilon[1][1],epsilon[1][2], V_max
		printf "					{%+.6f, %+.6f, %+.6f}  \tSum{ | epsilon | } = %.3g\r\r",epsilon[2][0],epsilon[2][1],epsilon[2][2], SumAbsEpsilon
		WaveStats/Q epsilonBL
		sprintf str,"  \t*****  Tr(epsilonBL) = %.2g  *****",MatrixTrace(epsilonBL)
		printf "epsilonBL =		{%+.6f, %+.6f, %+.6f}  \tmax{ | epsilon(ij) | } = %.3g\r",epsilonBL[0][0],epsilonBL[0][1],epsilonBL[0][2],max(abs(V_max),abs(V_min))
		printf "					{%+.6f, %+.6f, %+.6f}%s\r",epsilonBL[1][0],epsilonBL[1][1],epsilonBL[1][2], SelectString(abs(MatrixTrace(epsilonBL))>1e-12,"",str)
		printf "					{%+.6f, %+.6f, %+.6f}\r\r",epsilonBL[2][0],epsilonBL[2][1],epsilonBL[2][2]
		WaveStats/Q epsilonXHF
		sprintf str,"  \t*****  Tr(epsilonXHF) = %.2g  *****",MatrixTrace(epsilonXHF)
		printf "epsilonXHF =		{%+.6f, %+.6f, %+.6f}  \tmax{ | epsilon(ij) | } = %.3g\r",epsilonXHF[0][0],epsilonXHF[0][1],epsilonXHF[0][2],max(abs(V_max),abs(V_min))
		printf "					{%+.6f, %+.6f, %+.6f}%s\r",epsilonXHF[1][0],epsilonXHF[1][1],epsilonXHF[1][2], SelectString(abs(MatrixTrace(epsilonXHF))>1e-12,"",str)
		printf "					{%+.6f, %+.6f, %+.6f}\r\r",epsilonXHF[2][0],epsilonXHF[2][1],epsilonXHF[2][2]
		WaveStats/Q epsilonSample
		sprintf str,"  \t*****  Tr(epsilonSample) = %.2g  *****",MatrixTrace(epsilonSample)
		printf "epsilonSample =	{%+.6f, %+.6f, %+.6f}  \tmax{ | epsilon(ij) | } = %.3g\r",epsilonSample[0][0],epsilonSample[0][1],epsilonSample[0][2],max(abs(V_max),abs(V_min))
		printf "					{%+.6f, %+.6f, %+.6f}%s\r",epsilonSample[1][0],epsilonSample[1][1],epsilonSample[1][2], SelectString(abs(MatrixTrace(epsilonSample))>1e-12,"",str)
		printf "					{%+.6f, %+.6f, %+.6f}\r\r",epsilonSample[2][0],epsilonSample[2][1],epsilonSample[2][2]
	endif

	wnote = note(FullPeakIndexed)
	sprintf str, "{{%.8g,%.8g,%.8g}{%.8g,%.8g,%.8g}{%.8g,%.8g,%.8g}}",RL[0][0],RL[1][0],RL[2][0],  RL[0][1],RL[1][1],RL[2][1],  RL[0][2],RL[1][2],RL[2][2]
	wnote = ReplaceStringByKey("recip_lattice_refined",wnote,str,"=")
	sprintf str, "{%.8g %.8g %.8g %.8g %.8g %.8g}",a,b,c,alpha,bet,gam
	wnote = ReplaceStringByKey("lattice_constants_refined",wnote,str,"=")
	wnote = ReplaceNumberByKey("SumAbsEpsilon",wnote,SumAbsEpsilon,"=")
	wnote = ReplaceNumberByKey("vonMisesStrain",wnote,epsilonvM,"=")
	Note/K epsilon, ReplaceStringByKey("waveClass",wnote,"strain_tensor_STD","=")
	Note/K epsilonBL, ReplaceStringByKey("waveClass",wnote,"strain_tensor_BL","=")
	Note/K epsilonXHF, ReplaceStringByKey("waveClass",wnote,"strain_tensor_XHF","=")
	Note/K epsilonSample, ReplaceStringByKey("waveClass",wnote,"strain_tensor_Sample","=")

	wnote = note(PeaksForStrain)
	sprintf str, "{{%.8g,%.8g,%.8g}{%.8g,%.8g,%.8g}{%.8g,%.8g,%.8g}}",RL[0][0],RL[1][0],RL[2][0],  RL[0][1],RL[1][1],RL[2][1],  RL[0][2],RL[1][2],RL[2][2]
	wnote = ReplaceStringByKey("recip_lattice_refined",wnote,str,"=")
	sprintf str, "{%.8g %.8g %.8g %.8g %.8g %.8g}",a,b,c,alpha,bet,gam
	wnote = ReplaceStringByKey("lattice_constants_refined",wnote,str,"=")
	Note/K PeaksForStrain, wnote

	KillWaves/Z  axis_latticeMismatch, rho_latticeMismatch
	KillWaves/Z fithat_PeaksForStrain, optimize_LatticeConstantsWave
	KillWaves/Z RL_latticeMismatch, DL_latticeMismatch
	KillWaves/Z W_Cross
	KillWaves/Z M_Lower,M_Upper,W_LUPermutation,M_x

	if (coords==0)
		str = GetWavesDataFolder(epsilon,2)
	elseif (coords==1)
		str = GetWavesDataFolder(epsilonBL,2)
	elseif (coords==2)
		str = GetWavesDataFolder(epsilonXHF,2)
	elseif (coords==3)
		str = GetWavesDataFolder(epsilonSample,2)
	endif
	return str
End
//
Static Function latticeMismatchAll(PeaksForStrain,rx,ry,rz,LC)			// 3 rotations and all 6 lattice constants (uses MatrixOP for speed)
	Wave PeaksForStrain
	Variable rx,ry,rz
	Wave LC									//	lattice parameters:   LC = {a,b,c,alpha,bet,gam}

	Wave axis=axis_latticeMismatch, rho=rho_latticeMismatch
	Wave RL=RL_latticeMismatch, LC=optimize_LatticeConstantsWave
	if (!WaveExists(axis) || !WaveExists(rho) || !WaveExists(RL) || !WaveExists(LC))
		Abort "some of the waves do not exist in latticeMismatchAll()"
	endif
	Variable Vc = NumberByKey("Vc",note(PeaksForStrain),"=")		// This is to be held constant

	axis = {rx,ry,rz}
	Variable angle= norm(axis)
	rotationMatAboutAxis(axis,angle,rho)									// calculate rho, the rotation matrix from {rx,ry,rz}

	RLfromLatticeConstants(LC,$"",RL,Vc=Vc)								// calculate the reciprocal lattice from {a,b,c,alpha,bet,gam}
	MatrixOp/O RL = rho x RL													// rotate the reciprocal lattice by rho

	Variable N = DimSize(PeaksForStrain,0)
	Make/N=(N,3)/D/FREE qhats, ghats
	Make/N=(N)/D/FREE weight=1,  ones=1
	qhats = PeaksForStrain[p][q+3]											// measured directions (assumed normalized)
	ghats = PeaksForStrain[p][q]											// hkl of reflection
	MatrixOp/O/FREE ghats = NormalizeCols(RL x ghats^t)^t			// calcualted normalized direction from hkls
	PeaksForStrain[][6,8] = ghats[p][q-6]								// save for later checking
	// weight = sqrt(sqrt(PeaksForStrain[p][9]))
	MatrixOP/FREE errw = sum((ones-sumRows(ghats*qhats))*weight)	// 1-dot is ~ 0.5*∆theta^2,  using cos(theta) ~ 1-(theta^2)/2!
	Variable err =sqrt(errw[0]*2/sum(weight))							// change to rms
	err = numtype(err) ? Inf : err
	if (numtype(err) && NumVarOrDefault("printOnlyFirstStrainNaN",1))
		printf "err = %g,    LC={%g, %g, %g, %g, %g, %g}\r",err,LC[0],LC[1],LC[2],LC[3],LC[4],LC[5]
		Variable/G printOnlyFirstStrainNaN = 0
	endif
	return err	// units are radians
End
//
#endif			// end of #ifdef USE_ENERGY_STRAIN_REFINE
//
Static Function latticeMismatch4(PeaksForStrain,rx,ry,rz,f1)	// only fit rx,ry,rz, and one others
	Wave PeaksForStrain
	Variable rx,ry,rz
	Variable f1
	Wave LC=optimize_LatticeConstantsWave
	fillLC(LC,StringByKey("constrain",note(PeaksForStrain),"="),f1,NaN,NaN,NaN,NaN,NaN)
	return latticeMismatchAll(PeaksForStrain,rx,ry,rz,LC)
End
//
Static Function latticeMismatch5(PeaksForStrain,rx,ry,rz,f1,f2)	// only fit rx,ry,rz, and two others
	Wave PeaksForStrain
	Variable rx,ry,rz
	Variable f1,f2
	Wave LC=optimize_LatticeConstantsWave
	fillLC(LC,StringByKey("constrain",note(PeaksForStrain),"="),f1,f2,NaN,NaN,NaN,NaN)
	return latticeMismatchAll(PeaksForStrain,rx,ry,rz,LC)
End
//
Static Function latticeMismatch6(PeaksForStrain,rx,ry,rz,f1,f2,f3)	// only fit rx,ry,rz, and three others
	Wave PeaksForStrain
	Variable rx,ry,rz
	Variable f1,f2,f3
	Wave LC=optimize_LatticeConstantsWave
	fillLC(LC,StringByKey("constrain",note(PeaksForStrain),"="),f1,f2,f3,NaN,NaN,NaN)
	return latticeMismatchAll(PeaksForStrain,rx,ry,rz,LC)
End
//
Static Function latticeMismatch7(PeaksForStrain,rx,ry,rz,f1,f2,f3,f4)	// only fit rx,ry,rz, and four others
	Wave PeaksForStrain
	Variable rx,ry,rz
	Variable f1,f2,f3,f4
	Wave LC=optimize_LatticeConstantsWave
	fillLC(LC,StringByKey("constrain",note(PeaksForStrain),"="),f1,f2,f3,f4,NaN,NaN)
	return latticeMismatchAll(PeaksForStrain,rx,ry,rz,LC)
End
//
Static Function latticeMismatch8(PeaksForStrain,rx,ry,rz,f1,f2,f3,f4,f5)	// only fit rx,ry,rz, and five others
	Wave PeaksForStrain
	Variable rx,ry,rz
	Variable f1,f2,f3,f4,f5
	Wave LC=optimize_LatticeConstantsWave
	fillLC(LC,StringByKey("constrain",note(PeaksForStrain),"="),f1,f2,f3,f4,f5,NaN)
	return latticeMismatchAll(PeaksForStrain,rx,ry,rz,LC)
End
//
// This should never be called when optimizing Deviatoric strain, can be called when there are measured energies
Static Function latticeMismatch9(PeaksForStrain,rx,ry,rz,f1,f2,f3,f4,f5,f6)	// only fit rx,ry,rz, and five others
	Wave PeaksForStrain
	Variable rx,ry,rz
	Variable f1,f2,f3,f4,f5,,f6
	Wave LC=optimize_LatticeConstantsWave
	fillLC(LC,StringByKey("constrain",note(PeaksForStrain),"="),f1,f2,f3,f4,f5,f6)
	return latticeMismatchAll(PeaksForStrain,rx,ry,rz,LC)
End
//
Static Function fillLC(LC,constrain,f1,f2,f3,f4,f5,f6)
	Wave LC										// lattice constants
	String constrain
	Variable f1,f2,f3,f4,f5,f6

	Variable i
	i = strsearch(constrain,"1",0)
	LC[i] = f1
	i = strsearch(constrain,"1",i+1)
	if (i<0)
		return 0
	endif
	LC[i] = f2
	i = strsearch(constrain,"1",i+1)
	if (i<0)
		return 0
	endif
	LC[i] = f3
	i = strsearch(constrain,"1",i+1)
	if (i<0)
		return 0
	endif
	LC[i] = f4
	i = strsearch(constrain,"1",i+1)
	if (i<0)
		return 0
	endif
	LC[i] = f5
	i = strsearch(constrain,"1",i+1)
	if (i<0)
		return 0
	endif
	LC[i] = f6
	return 0
End


Static Function forceLattice(LC,SG,[Vc])
	Wave LC
	Variable SG
	Variable Vc										// if present, then preserve Vc

	STRUCT crystalStructure xtal				// this sruct is local to this routine
	xtal.SpaceGroup = SG
	xtal.a = LC[0]
	xtal.b = LC[1]
	xtal.c = LC[2]
	xtal.alpha = LC[3]
	xtal.beta = LC[4]
	xtal.gam = LC[5]
	LatticeSym#ForceLatticeToStructure(xtal)		// also updates xtal.Vc
	LC[0] = xtal.a
	LC[1] = xtal.b
	LC[2] = xtal.c
	LC[3] = xtal.alpha
	LC[4] = xtal.beta
	LC[5] = xtal.gam
	if (!ParamIsDefault(Vc) && (Vc>0))			// passed a valid Vc, preserve it by adjusting a,b,c
		Variable factor = (Vc/xtal.Vc)^(1/3)
		LC[0,2] *= factor
	endif
End


// The choice of DL from lattice constants is same as section Vol B, 3.3 of the International Tables
Static Function RLfromLatticeConstants(LC,DL,RL[,Vc])	// set a*, b*, c* from lattice constants, option to preserve Vc
	Wave LC											// wave with 6 lattice constants
	Wave DL, RL										// direct and recip lattices
	Variable Vc											// rescale so volume of unit cell is this (if valid Vc)
	Vc = ParamIsDefault(Vc) ? NaN : Vc

	Variable a=LC[0], b=LC[1], c=LC[2]
	Variable alpha=LC[3], bet=LC[4], gam=LC[5]

	// calculate the reciprocal lattice from {a,b,c,alpha,bet,gam}
	Variable sa = sin((alpha)*PI/180), ca = cos((alpha)*PI/180)
	Variable cb = cos((bet)*PI/180), cg = cos((gam)*PI/180)
	Variable phi = sqrt(1.0 - ca*ca - cb*cb - cg*cg + 2*ca*cb*cg)	// = Vc/(a*b*c), dimensionless
	if (Vc>0)											// re-scale {a,b,c} to preserve volume of unit cell (deviatoric strain only)
		Variable scale = (Vc/(a*b*c*phi))^(1/3)	// a*b*c*phi = Vc
		a *= scale
		b *= scale
		c *= scale
	else
		Vc = a*b*c * phi								// no Vc passed, use computed Vc
	endif

	Variable pv = (2*PI) / (Vc)						// used for scaling reciprocal lattice
	Variable a0,a1,a2,  b0,b1,b2,  c0,c1,c2		//components of the direct lattice vectors
	Variable rhom = (a==b && b==c && alpha==bet && bet==gam && abs(bet-90)>0.1)	// rhombohedral
	if (rhom)											// special for rhombohedral axes
		Variable pp=a/3*sqrt(1+2*ca), qq=a/3*sqrt(1-ca)
		a0=pp+2*qq	; a1=pp-qq				; a2=pp-qq
		b0=pp-qq		; b1=pp+2*qq				; b2=pp-qq
		c0=pp-qq		; c1=pp-qq				; c2=pp+2*qq
	else												// for all others
		a0=a*phi/sa	; a1=a*(cg-ca*cb)/sa	; a2=a*cb
		b0=0			; b1=b*sa					; b2=b*ca
		c0=0			; c1=0						; c2=c
	endif
	if (WaveExists(DL))								// set direct lattice
		DL[0][0] = a0	;	DL[0][1] = b0	;	DL[0][2] = c0
		DL[1][0] = a1	;	DL[1][1] = b1	;	DL[1][2] = c1
		DL[2][0] = a2	;	DL[2][1] = b2	;	DL[2][2] = c2
	endif
	if (WaveExists(RL))								// set reciprocal lattice
		// a* = (b x c)*2π/Vc				b* = (c x a)*2π/Vc					c*=(a x b)*2π/Vc
		RL[0][0]=(b1*c2-b2*c1)*pv	;	RL[0][1]=(c1*a2-c2*a1)*pv ;	RL[0][2]=(a1*b2-a2*b1)*pv
		RL[1][0]=(b2*c0-b0*c2)*pv ;	RL[1][1]=(c2*a0-c0*a2)*pv ;	RL[1][2]=(a2*b0-a0*b2)*pv
		RL[2][0]=(b0*c1-b1*c0)*pv ;	RL[2][1]=(c0*a1-c1*a0)*pv ;	RL[2][2]=(a0*b1-a1*b0)*pv
	endif
	return Vc
End


Static Function RL2latticeConstants(RL,LC)			// calc lattice constants from RL
	Wave RL											// recip lattice
	Wave LC											// wave to hold 6 lattice constants

	Make/N=3/O/D/FREE astar, bstar, cstar
	astar = RL[p][0]
	bstar = RL[p][1]
	cstar = RL[p][2]

	Cross bstar,cstar
	Wave W_Cross=W_Cross
	Variable pv = 2*PI/MatrixDot(astar,W_Cross)	// real lattice for un-strained material
	Make/N=3/O/D/FREE avec,bvec,cvec
	Cross bstar,cstar ;	avec = pv*W_Cross
	Cross cstar,astar ;	bvec = pv*W_Cross
	Cross astar,bstar ;	cvec = pv*W_Cross
	Cross bvec,cvec
	Variable Vc = MatrixDot(avec,W_Cross)			// volume of unstrained real unit cell
	LC[0] = norm(avec)							// the lattice parameters
	LC[1] = norm(bvec)
	LC[2] = norm(cvec)
	LC[3] = angleVec2Vec(bvec,cvec)
	LC[4] = angleVec2Vec(avec,cvec)
	LC[5] = angleVec2Vec(avec,bvec)
	return Vc
End


Static Function StartStrainPanel(ipat,npat,aIn,bIn,cIn,alphaIn,betIn,gamIn,[deviatoric])
	Variable &ipat, npat							// number of patterns, and initial pattern
	Variable &aIn,&bIn,&cIn, &alphaIn,&betIn,&gamIn
	Variable deviatoric							// if TRUE, cannot set a, b, and c
	deviatoric = ParamIsDefault(deviatoric) ? 1 : deviatoric
	deviatoric = numtype(deviatoric) ? 1 : !(!deviatoric)

	NewDataFolder/O root:Packages:micro
	NewDataFolder/O root:Packages:micro:strainRefine
	String path = "root:Packages:micro:strainRefine:"
	if (exists(path+"aFit")!=2)				// ensure globals exist
		Variable/G $(path+"aFit")=1,$(path+"bFit")=1, $(path+"cFit")=0
		Variable/G $(path+"alphaFit")=1,$(path+"betFit")=1, $(path+"gamFit")=1
		Variable/G $(path+"patFit")=0
		String/G $(path+"patternNumberList")="0"
	endif
	NVAR aFit=$(path+"aFit"), bFit=$(path+"bFit"), cFit=$(path+"cFit")
	NVAR alphaFit=$(path+"alphaFit"), betFit=$(path+"betFit"), gamFit=$(path+"gamFit")
	if (deviatoric)
		cFit = (aFit&&bFit&&cFit) ? 0 : cFit	// cannot set all three lengths
	endif
	SVAR patternNumberList=$(path+"patternNumberList")

	ipat = !(ipat>=0) ? 0 : ipat
	npat = !(npat>ipat) ? ipat+1 : npat
	patternNumberList = expandRange("0-"+num2istr(max(0,npat-1)),";")

	NewPanel /W=(865,44,1020,177)/K=1
	PopupMenu popupPattern,pos={12,2},size={111,20},title="pattern #",fSize=12, disable=(npat==1 ? 2 : 0)
	PopupMenu popupPattern,mode=4,bodyWidth= 52,popvalue=num2istr(ipat),value=StrVarOrDefault("root:Packages:micro:strainRefine:patternNumberList","0")
	CheckBox check_a,pos={22,33},size={28,23},title="a",fSize=18,variable=$(path+"aFit"), proc=Indexing#UpdateStartStrainCheck
	CheckBox check_b,pos={22,55},size={29,23},title="b",fSize=18,variable=$(path+"bFit"), proc=Indexing#UpdateStartStrainCheck
	CheckBox check_c,pos={22,77},size={28,23},title="c",fSize=18,variable=$(path+"cFit"), proc=Indexing#UpdateStartStrainCheck
	CheckBox check_alpha,pos={92,35},size={30,18},title="a",font="Symbol",fSize=18
	CheckBox check_alpha,variable=$(path+"alphaFit"), proc=Indexing#UpdateStartStrainCheck
	CheckBox check_bet,pos={92,57},size={28,18},title="b",font="Symbol",fSize=18
	CheckBox check_bet,variable=$(path+"betFit"), proc=Indexing#UpdateStartStrainCheck
	CheckBox check_gam,pos={92,79},size={26,18},title="g",font="Symbol",fSize=18
	CheckBox check_gam,variable=$(path+"gamFit"), proc=Indexing#UpdateStartStrainCheck
	Button DoIt_button,pos={11,107},size={50,20},proc=Indexing#StartStrainButtonProc,title="Go"
	Button Cancel_button,pos={86,107},size={50,20},proc=Indexing#StartStrainButtonProc,title="Cancel"
	SetWindow kwTopWin,userdata(deviatoric)=num2str(deviatoric)
	DoUpdate
	String win = StringFromList(0,WinList("*",";","WIN:64"))
	EnableDisableStartStrainCheck(win)
	PauseForUser $win

	aIn = NumVarOrDefault(path+"aFit",1)
	bIn = NumVarOrDefault(path+"bFit",1)
	cIn = NumVarOrDefault(path+"cFit",1)
	alphaIn = NumVarOrDefault(path+"alphaFit",1)
	betIn = NumVarOrDefault(path+"betFit",1)
	gamIn = NumVarOrDefault(path+"gamFit",1)
	ipat = NumVarOrDefault(path+"patFit",1)
	Variable Nfit = aIn+bIn+cIn+alphaIn+betIn+gamIn
	if (!(Nfit>0 && (Nfit+deviatoric)<=6))
		ipat = -1
	endif
	aIn = aIn ? 1 : 0				// ensure only 0 or 1
	bIn = bIn ? 1 : 0
	cIn = cIn ? 1 : 0
	alphaIn = alphaIn ? 1 : 0
	betIn = betIn ? 1 : 0
	gamIn = gamIn ? 1 : 0
	return ipat					// <0 is error, >=0 is OK
End
//
Static Function StartStrainButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	if (ba.eventCode != 2)		// not mouse up
		return 0
	endif
	Variable pattern
	Variable/G root:Packages:micro:strainRefine:patFit
	NVAR patFit= root:Packages:micro:strainRefine:patFit
	if (stringmatch(ba.ctrlName,"DoIt_button"))
		ControlInfo/W=$(ba.win) popupPattern
		patFit = V_flag ? str2num(S_Value): 0				// if no control present, then only pattern 0 allowed
	else
		patFit = -1				// flags a "Cancel"
	endif
	DoWindow/K $(ba.win)
	KillStrings/Z patternNumberList
	return 0
End
//
Static Function UpdateStartStrainCheck(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba
	if (cba.eventCode != 2)
		return 0
	endif
	EnableDisableStartStrainCheck(cba.win)
	return 0
End
//
Static Function EnableDisableStartStrainCheck(win)
	String win
	DoUpdate
	String ctrlName, list="check_a;check_b;check_c;check_alpha;check_bet;check_gam"
	Variable i, n, nAll
	Variable deviatoric=str2num(GetUserData(win,"","deviatoric"))
	deviatoric = numtype(deviatoric) ? 1 : !(!deviatoric)	// deviatoric TRUE means you cannot fit a, b, and c at once

	for (n=0,i=0; i<3; i+=1)				//	count how many boxes are checked only a,b,c
		ctrlName = StringFromList(i,list)
		ControlInfo/W=$win $ctrlName
		n += V_Value
	endfor
	if (!deviatoric || n<2)
		for (i=0;i<3;i+=1)					// enable all
			ctrlName = StringFromList(i,list)
			CheckBox $ctrlName,disable=0
		endfor
	elseif (n==3)							// un-check c box, for deviatoric, you cannot fit a,b, & c (unless deviatoric=1)
		CheckBox check_c,value=0,disable=2
		n = 5
	elseif(n==2)							// disable the only un-checked box (unless deviatoric=1)
		for (i=0;i<3;i+=1)
			ctrlName = StringFromList(i,list)
			ControlInfo/W=$win $ctrlName
			CheckBox $ctrlName,disable=(V_Value ? 0 : 2)
		endfor
	endif

	for (nAll=n,i=3; i<6; i+=1)			//	count how many boxes are checked, all lattice constants
		ctrlName = StringFromList(i,list)
		ControlInfo/W=$win $ctrlName
		nAll += V_Value
	endfor
	if (deviatoric)
		Button DoIt_button,disable=((nAll<1 || n>2 || nAll>5) ? 2 : 0)
	else
		Button DoIt_button,disable=((nAll<1) ? 2 : 0)
	endif
	return 0
End

Static Function/WAVE FullPeakList2Qwave(FullPeakList, [printIt]) // convert fitted peaks to a Qlist+intens, peaks file for analysis be Euler
	// returns a FREE wave!!
	Wave FullPeakList
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt
	if (!WaveExists(FullPeakList))
		if (printIt)
			DoAlert 0, "input wave for FullPeakList2File() does not exists"
		endif
		return $""
	elseif (DimSize(FullPeakList,1)<11)
		if (printIt)
			DoAlert 0, "the passed full peak list '"+NameOfWave(FullPeakList)+"' is not the right size"
		endif
		return $""
	endif
	Variable N=DimSize(FullPeakList,0)
	if (N<1)												// nothing to write
		return $""
	endif
	STRUCT microGeometry geo
	if (FillGeometryStructDefault(geo, alert=1))	//fill the geometry structure with test values
		return $""
	endif

	Make/N=(N,4)/FREE/D Qs=NaN					// create temp waves for the q + intensity
	Make/N=3/FREE/D qhat
	String wnote = note(FullPeakList)
	Variable dNum = max(detectorNumFromID(geo, StringByKey("detectorID", wnote,"=")),0)
	Variable startx,groupx, starty,groupy		// ROI of the original image
	startx = NumberByKey("startx",wnote,"=")
	groupx = NumberByKey("groupx",wnote,"=")
	starty = NumberByKey("starty",wnote,"=")
	groupy = NumberByKey("groupy",wnote,"=")
	startx = numtype(startx) ? FIRST_PIXEL : startx
	groupx = numtype(groupx) ? 1 : groupx
	starty = numtype(starty) ? FIRST_PIXEL : starty
	groupy = numtype(groupy) ? 1 : groupy
	Variable depth = NumberByKey("depth",wnote,"=")

	Variable i, m, px,py
	for (i=0,m=0;i<N;i+=1)
		px = (startx-FIRST_PIXEL) + groupx*FullPeakList[i][0] + (groupx-1)/2		// change to un-binned pixels
		py = (starty-FIRST_PIXEL) + groupy*FullPeakList[i][1] + (groupy-1)/2		// pixels are still zero based
		pixel2q(geo.d[dNum],px,py,qhat,depth=depth)// in Beam Line Coord system
		if (norm(qhat)>0)								// check for a valid Q
			normalize(qhat)
			Qs[m][0,2] = qhat[q]
			Qs[m][3] = max(FullPeakList[i][10],0) // save intensity too
			m += 1
		endif
	endfor
	N = m
	Redimension/N=(N,4) Qs
	ImageStats/M=1/G={0,N-1,3,3} Qs
	Qs[][3] /= V_max
	return Qs
End
//Static Function/T FullPeakList2Qwave(FullPeakList) // convert fitted peaks to a Qlist+intens, peaks file for analysis be Euler
//	Wave FullPeakList
//	if (!WaveExists(FullPeakList))
//		DoAlert 0, "input wave for FullPeakList2File() does not exists"
//		return ""
//	elseif (DimSize(FullPeakList,1)!=11)
//		DoAlert 0, "the passed full peak list '"+NameOfWave(FullPeakList)+"' is not the right size"
//		return ""
//	endif
//	Variable N=DimSize(FullPeakList,0)
//	if (N<1)											// nothing to write
//		return ""
//	endif
//	STRUCT microGeometry geo
//	if (FillGeometryStructDefault(geo))			//fill the geometry structure with test values
//		DoAlert 0, "no geometry structure found, did you forget to set it?"
//		return ""
//	endif
//
//	Make/N=(N,4)/O/D FullPeakList2Qfile_Q=NaN // create temp waves for the q + intensity
//	Wave Qs=FullPeakList2Qfile_Q					// holds Qhat and intensity
//	Make/N=3/O/D FullPeakList2Qfile_qhat
//	Wave qhat=FullPeakList2Qfile_qhat
//	String wnote = note(FullPeakList)
//	Variable dNum = max(detectorNumFromID(StringByKey(geo, "detectorID", wnote,"=")),0)
//	Variable startx,groupx, starty,groupy			// ROI of the original image
//	startx = NumberByKey("startx",wnote,"=")
//	groupx = NumberByKey("groupx",wnote,"=")
//	starty = NumberByKey("starty",wnote,"=")
//	groupy = NumberByKey("groupy",wnote,"=")
//	startx = numtype(startx) ? FIRST_PIXEL : startx
//	groupx = numtype(groupx) ? 1 : groupx
//	starty = numtype(starty) ? FIRST_PIXEL : starty
//	groupy = numtype(groupy) ? 1 : groupy
//	Variable i, m, px,py
//	for (i=0,m=0;i<N;i+=1)
//		px = (startx-FIRST_PIXEL) + groupx*FullPeakList[i][0] + (groupx-1)/2		// change to un-binned pixels
//		py = (starty-FIRST_PIXEL) + groupy*FullPeakList[i][1] + (groupy-1)/2		// pixels are still zero based
//		pixel2q(geo.d[dNum],px,py,qhat)			// in Beam Line Coord system
//		if (norm(qhat)>0)							// check for a valid Q
//			normalize(qhat)
//			Qs[m][0,2] = qhat[q]
//			Qs[m][3] = max(FullPeakList[i][10],0) // save intensity too
//			m += 1
//		endif
//	endfor
//	KillWaves/Z FullPeakList2Qfile_qhat
//	N = m
//	Redimension/N=(N,4) FullPeakList2Qfile_Q
//	ImageStats/M=1/G={0,N-1,3,3} FullPeakList2Qfile_Q
//	FullPeakList2Qfile_Q[][3] /= V_max
//	return GetWavesDataFolder(FullPeakList2Qfile_Q,2)
//End
//
Function sqrtSymmetricMat(symmetricMat)			// computes sqrt(symmetricMat),  returns 1=error, 0=OK,   symmetricMat must be symmetric
	Wave symmetricMat									// Note, symmetricMat must be a "positive definite" matrix, or some of the eigen values will be negative

	MatrixEigenV /SYM/EVEC/Z symmetricMat
	if (V_flag)
		return 1
	endif
	Wave eigenValues = W_eigenValues, V = M_eigenVectors

//	// start of check
//	MatrixOp/O diagCheck = Inv(V) x symmetricMat x V
//	diagCheck = abs(diagCheck[p][q]) > 1e-15 ? 1 : 0
//	print diagCheck
//	// end of chek

	String wName = UniqueName("sqrtMat",1,0)
	Duplicate/O symmetricMat $wName
	Wave D = $wName								// holds the diagonal wave
	D = (p==q) ? sqrt(eigenValues[p]) : 0			// put sqrt(eigen values) on the diagonal, eigenValues better be >=0
	MatrixOp/O symmetricMat = V x D x Inv(V)	// transform back, putting root into original matrix
	KillWaves/Z W_eigenValues, M_eigenVectors, D
	return 0
End


Function VolumeOfLattice(a0,a1,a2,b0,b1,b2,c0,c1,c2)		// Volume = a .dot. (b .cross. c), the triple product
	Variable a0,a1,a2,b0,b1,b2,c0,c1,c2
	Variable Vc=0
	Vc += (b1*c2 - b2*c1)*a0
	Vc += (b2*c0 - b0*c2)*a1
	Vc += (b0*c1 - b1*c0)*a2
	return Vc
End

//	End of Index strain refinement
// =========================================================================
// =========================================================================



// =========================================================================
// =========================================================================
//	Start of Stereographic projection

Function StereoOfIndexedPattern(FullPeakIndexed,pattern,[centerType,showDetectors,printIt])
	Wave FullPeakIndexed							// provides the reciprocal lattice
	Variable pattern									// in case more than one, default is 0
	String centerType									// choice for center of stereographic projections, "Surface Normal" or "Detector Center"
	Variable showDetectors							// flag, show detector outlines on seterographic projection
	Variable printIt
	centerType = SelectString(ParamIsDefault(centerType), centerType, "")
	showDetectors = ParamIsDefault(showDetectors) ? NaN : showDetectors
	showDetectors = numtype(showDetectors) ? 1 : !(!showDetectors)
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt
	if (exists("MakeStereo")!=6)
		print "Loading Stereographic package"
		Execute/P "INSERTINCLUDE  \"StereographicProjection\", version>=2.86";Execute/P "COMPILEPROCEDURES "
		return 1
	endif
	if (!WaveExists(FullPeakIndexed))
		String wList=WaveListClass("IndexedPeakList*","*","")
		String wName=StringFromList(0,wList)
		if (ItemsInList(wList)!=1)						// only ask if more than one choice
			Prompt wName, "name of Full Peak Indexed list",popup,WaveListClass("IndexedPeakList*","*","")
			DoPrompt "Full Peak Indexed",wName
			if (V_flag)
				return 1
			endif
		endif
		Wave FullPeakIndexed=$wName
		printIt = 1
	endif
	if (!WaveExists(FullPeakIndexed))
		return 1
	endif
	Variable Npatterns=DimSize(FullPeakIndexed,2)
	pattern = (Npatterns==1) ? 0 : pattern
	pattern = limit(pattern,0,Npatterns-1)
	if (numtype(pattern))
		Prompt pattern, "pattern number",popup,expandRange("0-"+num2istr(Npatterns-1),";")
		DoPrompt "pattern number",pattern
		if (V_flag)
			return 1
		endif
		pattern -= 1
		printIt = 1
	endif

	Make/N=3/O/D/FREE vec
	Make/N=(3,3)/O/D/FREE RL
	Variable as0,as1,as2,bs0,bs1,bs2,cs0,cs1,cs2
	String wnote = note(FullPeakIndexed)
	sscanf StringByKey("recip_lattice"+num2istr(pattern),wnote,"="), "{{%g,%g,%g}{%g,%g,%g}{%g,%g,%g}}",as0,as1,as2,bs0,bs1,bs2,cs0,cs1,cs2
	RL[0][0]= {as0,as1,as2}								// the original measured RL
	RL[0][1]= {bs0,bs1,bs2}
	RL[0][2]= {cs0,cs1,cs2}

	// find (hkl) at center of detector and make that the center of the stereographic pattern
	STRUCT microGeometry geo
	if (FillGeometryStructDefault(geo, alert=1))	//fill the geometry structure with current default values
		return 1
	endif
	Variable startx=NumberByKey("startx",wnote,"="), endx=NumberByKey("endx",wnote,"=")
	Variable starty=NumberByKey("starty",wnote,"="), endy=NumberByKey("endy",wnote,"=")
	Variable depth = NumberByKey("depth",wnote,"=")

	Variable sufNormal=NaN
	if (stringmatch(centerType,"Surface*") || numtype(startx+endx+starty+endy))
		sufNormal = 1
		centerType = "Surface Normal"
	elseif (stringmatch(centerType,"Detector*"))
		sufNormal = 0
		centerType = "Detector Center"
	endif
	if (numtype(sufNormal))
		Prompt centerType,"Choose Center of Stereographic Projection",popup,"Surface Normal;Detector Center"
		Prompt showDetectors,"Show Detector Outlines",popup,"No Outlines;Show Detector Outlines"
		showDetectors += 1
		DoPrompt "Center",centerType,showDetectors
		if (V_flag)
			return 1
		endif
		showDetectors -= 1
		sufNormal = stringmatch(centerType"Surface*")
	endif

	String outlines=""										// string with detector outlines
	if (showDetectors)										// make detector outlines, store in a string named outlines
		// format of outlines is detctorOutline0;detctorOutline1;detctorOutline2;
		// detctorOutline0 = "r=65535:g=43688:b=32768:hkl0=-0.2 0.5 -0.8:hkl1=-0.3 0.7 -0.5:hkl2=0.3 0.8 -0.5:hkl3=0.2 0.5 -0.8:;
		// each of the hkli are a normalized hkl vector
		String str
		Make/N=3/FREE rgb
		Variable i, Nx,Ny
		for (i=0;i<MAX_Ndetectors;i+=1)
			if (geo.d[i].used)								// this detector is used, display it
				Nx = geo.d[i].Nx
				Ny = geo.d[i].Ny

				Wave rgb = detectorRGBwave(geo,i)
				str = ReplaceStringByKey("rgb","",vec2str(rgb,bare=1,sep=" "),"=",":")

				if (!(pixel2q(geo.d[i],0,0, vec,depth=depth)))
					continue
				endif
				MatrixOp/O/FREE hkl = Normalize(Inv(RL) x vec)	// go from q in Ideal beam line to hkl (but wrong length)
				str = ReplaceStringByKey("hkl0",str,vec2str(hkl,bare=1,sep=" "),"=",":")

				if (!(pixel2q(geo.d[i],Nx-1,0, vec,depth=depth)))
					continue
				endif
				MatrixOp/O/FREE hkl = Normalize(Inv(RL) x vec)
				str = ReplaceStringByKey("hkl1",str,vec2str(hkl,bare=1,sep=" "),"=",":")

				if (!(pixel2q(geo.d[i],Nx-1,Ny-1, vec,depth=depth)))
					continue
				endif
				MatrixOp/O/FREE hkl = Normalize(Inv(RL) x vec)
				str = ReplaceStringByKey("hkl2",str,vec2str(hkl,bare=1,sep=" "),"=",":")

				if (!(pixel2q(geo.d[i],0,Ny-1, vec,depth=depth)))
					continue
				endif
				MatrixOp/O/FREE hkl = Normalize(Inv(RL) x vec)
				str = ReplaceStringByKey("hkl3",str,vec2str(hkl,bare=1,sep=" "),"=",":")
				outlines += str+";"
			endif
		endfor
	endif

	vec = NaN													// q-vector at center of stereographic projection
	if (!sufNormal)											// use center of detector for central hkl
		pixel2q(geo.d[0],(startx+endx)/2,(starty+endy)/2, vec,depth=depth)
	endif
	if (numtype(norm(vec)))								// default if other way failed, use hkl of surface normal
		Variable cosTheta = NumVarOrDefault("root:Packages:geometry:cosThetaWire",cos(PI/4))
		Variable sinTheta = NumVarOrDefault("root:Packages:geometry:sinThetaWire",sin(PI/4))
		vec = {0, cosTheta, -sinTheta}					// vec = {0,1,-1}
		sufNormal = 1
	endif

	MatrixOp/O/FREE hkl = Normalize(Inv(RL) x vec)// go from q in Ideal beam line to hkl (but wrong length)
	hkl = round(hkl*10000)/1000
	hkl = hkl==0 ? 0 : hkl
	Variable h=hkl[0], k=hkl[1], l=hkl[2]
	lowestOrderHKL(h,k,l)									// remove common factors
	vec = {1,0,0}												// now get azimuthal orientation
	MatrixOp/O/FREE hklAz = Normalize(Inv(RL) x vec)
	hklAz = round(hklAz*10000)/1000
	if (printIt)
		if (stringmatch(StringFromList(0,GetRTStackInfo(0)),"IndexButtonProc"))
			printf "•"
		endif
		printf "StereoOfIndexedPattern(%s,%d)\r",NameOfWave(FullPeakIndexed),pattern
		printf "Making Stereographic projection with pole at (%g, %g, %g),  and  x toward (%g, %g, %g)",h,k,l,round(hklAz[0]*1000)/100, round(hklAz[1]*1000)/100, round(hklAz[2]*1000)/100
		printf ",  centered on "+SelectString(sufNormal,"surface normal","detector")+"\r"
	endif
	FUNCREF MakeStereoProto func = $"MakeStereo"
	func(h,k,l,  NaN,9, 0,hklPerp=hklAz,outlines=outlines)	// MakeStereo(h,k,l,  NaN,9, 0,hklPerp=hklAz)
	return 0
End
Function MakeStereoProto(Hz,Kz,Lz,hklmax,fntSize,phi,[Qmax,hklPerp, WulffStepIn,WulffPhiIn,markers,outlines])
	Variable Hz,Kz,Lz			// hkl of the pole
	Variable hklmax
	Variable fntSize
	Variable phi				// aximuthal angle (degree)
	Variable Qmax				// maximum Q (1/nm)
	Wave hklPerp				// optional hkl giving the perpendicular direction long=0
	Variable WulffStepIn, WulffPhiIn	// for Wulff Net, use WulffStep=0 for no Wulff Net
	Variable markers			// flag, show markers for each hkl
	String outlines			// string with outlines of detectors to display
End

//	End of Stereographic projection
// =========================================================================
// =========================================================================



// =========================================================================
// =========================================================================
//	Start of Index set panel

Function/T FillIndexParametersPanel(strStruct,hostWin,left,top)
	String strStruct									// optional passed value of xtal structure, this is used if passed
	String hostWin										// name of home window
	Variable left, top									// offsets from the left and top

	NewDataFolder/O root:Packages:micro				// ensure that the needed data folders exist
	NewDataFolder/O root:Packages:micro:Index
	NewDataFolder/O root:Packages:Lattices:PanelValues
	String/G root:Packages:Lattices:PanelValues:desc
	Variable/G root:Packages:micro:Index:h_e
	Variable/G root:Packages:micro:Index:k_e
	Variable/G root:Packages:micro:Index:l_e
	NVAR h=root:Packages:micro:Index:h_e
	NVAR k=root:Packages:micro:Index:k_e
	NVAR l=root:Packages:micro:Index:l_e

	SetWindow kwTopWin,userdata(IndexPanelName)=hostWin+"#IndexPanel"
	NewPanel/K=1/W=(left,top,left+221,top+560+30)/HOST=$hostWin
	ModifyPanel frameStyle=0, frameInset=0
	RenameWindow #,IndexPanel

	Button buttonLoadImage,pos={1,5},size={120,20},proc=IndexButtonProc,title="Load Image"
	Button buttonLoadImage,help={"Load an Image from the file, and then display the image."}
	Button buttonViewImage,pos={130,5},size={90,20},proc=IndexButtonProc,title="View Image"
	Button buttonViewImage,help={"Display an existing Image that is already loaded"}

	Button buttonMakeMask,pos={1,35},size={120,20},proc=IndexButtonProc,title="Make Mask"
	Button buttonMakeMask,help={"Make a Mask of Peaks to Exclude from peak searching"}
	Button buttonEditMask,pos={130,35},size={90,20},proc=IndexButtonProc,title="Edit Mask"
	Button buttonEditMask,help={"Edit an existing mask"}

	Button buttonFitPeaks,pos={1,65},size={120,20},proc=IndexButtonProc,title="Fit Peaks"
	Button buttonFitPeaks,help={"Identify peaks in an Image, does not do any background removal"}
	Button buttonBkgRemove,pos={130,65},size={90,20},proc=IndexButtonProc,title="Background"
	Button buttonBkgRemove,help={"remove background from an image, do not use on depth sotred images"}

	Button buttonIndex,pos={1,95},size={120,20},proc=IndexButtonProc,title="Index & Display"
	Button buttonIndex,help={"Index the identified peaks, using the defined lattice, and then show it all on a plot"}
	Button buttonAuxIndex,pos={130,95},size={90,20},proc=IndexButtonProc,title="Index Aux."
	Button buttonAuxIndex,help={"Index the identified peaks, using the defined lattice, and then show it all on a plot"}

	PopupMenu popuphklTags,pos={1,125},size={120,20},proc=hklTagsPopMenuProc,title="hkl Tags..."
	PopupMenu popuphklTags,help={"re-draw or clear the hkl tags from a plot"}
	PopupMenu popuphklTags,fSize=14
	PopupMenu popuphklTags,mode=0,value= #"\"Re-Draw hkl tags;Clear hkl tags off Graph;Add Missing Peaks;Remove Missing Peaks;Rotation Between Orientations...;\""
	if (strlen(WaveList("*",";", "DIMS:2,MINROWS:3,MAXROWS:3,MINCOLS:3,MAXCOLS:3"))>1)
		PopupMenu popuphklTags,value= #"\"Re-Draw hkl tags;Clear hkl tags off Graph;Add Missing Peaks;Remove Missing Peaks;Rotation Between Orientations...;Add Peaks from a user recip;Remove All Peaks from user recips;\""
	endif
	Button buttonReportIndex,pos={130,125},size={90,20},proc=IndexButtonProc,title="Report"
	Button buttonReportIndex,help={"if the top window is an images showing the hkl's, this makes a layout with reflections for printing"}

	Button buttonStrainRefine,pos={29,155},size={160,20},proc=IndexButtonProc,title="Strain Refine"
	Button buttonStrainRefine,help={"after indexing, use this to refine the lattice and find the strain"}
	if (exists("MakeStereo")==6)
		Button buttonStereographic,pos={29,185},size={160,20},proc=IndexButtonProc,title="Stereographic..."
	else
		Button buttonStereographic,pos={29,185},size={160,20},proc=IndexButtonProc,title="Load Stereographic"
	endif
	Button buttonStereographic,help={"Make a Stereographic projection with the same orientation as the indexed pattern"}

	Button buttonCalcEnergyIndex,pos={29,185+30},size={160,20},proc=IndexButtonProc,title="Calc Energy of (hkl)"
	Button buttonCalcEnergyIndex,help={"calculate the energy of an hkl using the current indexing"}
	SetVariable h_IndexEneVar,pos={27,212+30},size={50,15},fsize=12,proc=CalcE_SetVarProc,title="H"
	SetVariable h_IndexEneVar,value= root:Packages:micro:Index:h_e
	SetVariable k_IndexEneVar,pos={91,212+30},size={50,15},fsize=12,proc=CalcE_SetVarProc,title="K"
	SetVariable k_IndexEneVar,value= root:Packages:micro:Index:k_e
	SetVariable L_IndexEneVar,pos={155,212+30},size={50,15},fsize=12,proc=CalcE_SetVarProc,title="L"
	SetVariable L_IndexEneVar,value= root:Packages:micro:Index:l_e

	PopupMenu popupTables,pos={62,240+30},size={66,20},proc=TablesPopMenuProc,title="Display Tables..."
	PopupMenu popupTables,help={"Show table of identified peak positions indexed peaks, or the peaks for strain refinement"}
	PopupMenu popupTables,fSize=14,mode=0,value= #"\"Fitted Peaks;Indexed Peaks;Strain Peaks;Missing Peaks;Measured Energies...\""
	if (strlen(WaveListClass("UserRecipPeakList","*","")))
		PopupMenu popupTables,value= #"\"Fitted Peaks;Indexed Peaks;Strain Peaks;Missing Peaks;Measured Energies...;Peaks from User recip\""
	endif

	Variable offset = 0
	if (exists("Load3dRecipLatticesFileXML")==6)			// do this if function exists (it is in xmlMultiIndex.ipf)
		offset += 305
		SetDrawLayer UserBack
		SetDrawEnv linethick= 2
		DrawLine 0,offset,221,offset

		Button buttonIndexNewPoleXML,pos={1,15+offset},size={140,20},proc=IndexButtonProc,title="New Pole Figure"
		Button buttonIndexNewPoleXML,help={"Create a Pole Figure from current xml file data."}
		Button buttonIndexRePlotPoleXML,pos={150,15+offset},size={70,20},proc=IndexButtonProc,title="Re-Plot"
		Button buttonIndexRePlotPoleXML,help={"Re-Plot the calculated Pole Figure"}

		Button buttonIndexRefineStrainMany,pos={1,40+offset},size={140,20},proc=IndexButtonProc,title="Strain Refine Many"
		Button buttonIndexRefineStrainMany,help={"Strain Refine Many Indexed Images"}
		Button buttonIndexRefineStrain1,pos={150,40+offset},size={70,20},proc=IndexButtonProc,title="Refine 1"
		Button buttonIndexRefineStrain1,help={"Strain Refine 1 Indexed Image"}

		SetDrawLayer UserBack
		SetDrawEnv linethick= 1
		DrawLine 0,65+offset,221,65+offset

		PopupMenu popupIndexLoadXML,pos={1,70+offset},size={140,20},fSize=14,proc=Indexing#IndexXMLPopUpMenuProc,title="Load 3d XML"
		PopupMenu popupIndexLoadXML,mode=0,help={"Load & Prepare for Display an XML output file"}
		FUNCREF ValidRawXMLdataAvailableProto fvalid = $"multiIndex#ValidRawXMLdataAvailable"
		String mstr="\"Load New 3D XML file;"+SelectString(fvalid(""),"","Load/Append 3D XML file;Process Loaded XML;")+"\""
		PopupMenu popupIndexLoadXML,value=#mstr
		Button buttonIndexLoadOneXML,pos={150,70+offset},size={70,20},proc=IndexButtonProc,title="One Step"
		Button buttonIndexLoadOneXML,help={"Load One <step> from an XML file"}

//		Button buttonIndex2D3DXMLrgba,pos={1,95+offset},size={219,20},proc=IndexButtonProc,title="New 2D/3D RGBA"
		Button buttonIndex2D3DXMLrgba,pos={1,95+offset},size={140,20},proc=IndexButtonProc,title="New 2D/3D RGBA"
		Button buttonIndex2D3DXMLrgba,help={"calculate a new RGBA wave for use with Gizmo"}
		Button buttonIndex2Dplot,pos={150,95+offset},size={70,20},proc=IndexButtonProc,title="2D Plot"
		Button buttonIndex2Dplot,help={"Put up a 2D plot of loaded XML data"}

		Button buttonIndexGizmoXML,pos={1,120+offset},size={140,20},proc=IndexButtonProc,title="Gizmo of 3d XML"
		Button buttonIndexGizmoXML,help={"Make a Gizmo from the XML data"}
		Button buttonIndexColorHexagon,pos={150,120+offset},size={70,20},proc=IndexButtonProc,title="Color Hex"
		Button buttonIndexColorHexagon,help={"Display direction colors"}

		Button buttonIndexHistogram,pos={1,145+offset},size={140,20},proc=IndexButtonProc,title="Histogram 3D array"
		Button buttonIndexHistogram,help={"Create a histogram of an interpolated 3D array"}
		Button buttonIndexHistogramPlot,pos={150,145+offset},size={70,20},proc=IndexButtonProc,title="Re-Plot"
		Button buttonIndexHistogramPlot,help={"Re-Plot the histogram of an interpolated 3D array"}

		Button buttonIndexRotation,pos={1,170+offset},size={190,20},proc=IndexButtonProc,title="Rotation between 2 Points"
		Button buttonIndexRotation,help={"Find the rotation angle between two points, probably identified with cursors A & B"}

		offset += 125 + 30 + 25
	endif
	if (exists("IndexLots")==6)
		offset += 25
		SetDrawLayer UserBack
		SetDrawEnv linethick= 2
		DrawLine 0,offset,221,offset
		Button buttonIndexLots,pos={29,10+offset},size={160,20},proc=IndexButtonProc,title="Index Lots"
		Button buttonIndexLots,help={"index lots of images"}
		Button buttonWrite3dFile,pos={39,35+offset},size={160,20},proc=IndexButtonProc,title="Write 3d Orients file"
		Button buttonWrite3dFile,help={"Write 3d Orientations file"}
		Button buttonFindClosestPoint,pos={1,60+offset},size={140,20},proc=IndexButtonProc,title="Find Closest Point"
		Button buttonFindClosestPoint,help={"find point in original indexing list that is closest to the cursor"}
		Button buttonIndexingListTable,pos={150,60+offset},size={70,20},proc=IndexButtonProc,title="Table"
		Button buttonIndexingListTable,help={"Make a table of the indexing list"}
	else
		Button buttonInitIndexLots,pos={29,340+offset},size={150,50},title="Init IndexLots\rpackage",proc=microGeo#LoadPackageButtonProc
		Button buttonInitIndexLots,help={"Load Package for Indexing Lots of Images"}
	endif

	Indexing#EnableDisableIndexControls(hostWin+"#IndexPanel")
	return "#IndexPanel"
End
//Function/T FillIndexParametersPanel(strStruct,hostWin,left,top)
//	String strStruct									// optional passed value of xtal structure, this is used if passed
//	String hostWin										// name of home window
//	Variable left, top									// offsets from the left and top
//
//	NewDataFolder/O root:Packages:micro				// ensure that the needed data folders exist
//	NewDataFolder/O root:Packages:micro:Index
//	String/G root:Packages:Lattices:PanelValues:desc
//	Variable/G root:Packages:micro:Index:h_e
//	Variable/G root:Packages:micro:Index:k_e
//	Variable/G root:Packages:micro:Index:l_e
//	NVAR h=root:Packages:micro:Index:h_e
//	NVAR k=root:Packages:micro:Index:k_e
//	NVAR l=root:Packages:micro:Index:l_e
//
//	SetWindow kwTopWin,userdata(IndexPanelName)=hostWin+"#IndexPanel"
//	NewPanel/K=1/W=(left,top,left+221,top+445+30)/HOST=$hostWin
//	ModifyPanel frameStyle=0, frameInset=0
//	RenameWindow #,IndexPanel
//
//	Button buttonLoadImage,pos={17,5},size={80,20},proc=IndexButtonProc,title="Load Image"
//	Button buttonLoadImage,help={"Load an Image from the file, and then display the image."}
//	Button buttonViewImage,pos={111,5},size={80,20},proc=IndexButtonProc,title="View Image"
//	Button buttonViewImage,help={"Display an existing Image that is already loaded"}
//
//	Button buttonBkgRemove,pos={29,40},size={160,20},proc=IndexButtonProc,title="Remove Background"
//	Button buttonBkgRemove,help={"remove background from an image, do not use on depth sotred images"}
//	Button buttonFitPeaks,pos={29,65},size={160,20},proc=IndexButtonProc,title="Fit Peaks"
//	Button buttonFitPeaks,help={"Identify peaks in an Image, does not do any background removal"}
//
////	Button buttonPeaksFromBoxes,pos={29,90},size={160,20},proc=IndexButtonProc,title="set peaks from boxes"
////	Button buttonPeaksFromBoxes,help={"reset list of peaks for fitting to the boxes on an image plot"}
//	Button buttonIndex,pos={29,100},size={160,20},proc=IndexButtonProc,title="Index & Display"
//	Button buttonIndex,help={"Index the identified peaks, using the defined lattice, and then show it all on a plot"}
//	Button buttonAuxIndex,pos={29,125},size={160,20},proc=IndexButtonProc,title="Index Aux. Detectors"
//	Button buttonAuxIndex,help={"Index the identified peaks, using the defined lattice, and then show it all on a plot"}
////	Button buttonIndex,pos={29,125},size={160,20},proc=IndexButtonProc,title="Index & Display"
////	Button buttonIndex,help={"Index the identified peaks, using the defined lattice, and then show it all on a plot"}
//
//	PopupMenu popuphklTags,pos={69,150},size={76,20},proc=hklTagsPopMenuProc,title="hkl Tags..."
//	PopupMenu popuphklTags,help={"re-draw or clear the hkl tags from a plot"}
//	PopupMenu popuphklTags,fSize=14
//	PopupMenu popuphklTags,mode=0,value= #"\"Re-Draw hkl tags;Clear hkl tags off Graph;Add Missing Peaks;Remove Missing Peaks\""
//	Button buttonStrainRefine,pos={29,175},size={160,20},proc=IndexButtonProc,title="Strain Refine"
//	Button buttonStrainRefine,help={"after indexing, use this to refine the lattice and find the strain"}
//	if (exists("MakeStereo")==6)
//		Button buttonStereographic,pos={29,205},size={160,20},proc=IndexButtonProc,title="Stereographic..."
//	else
//		Button buttonStereographic,pos={29,205},size={160,20},proc=IndexButtonProc,title="Load Stereographic"
//	endif
//	Button buttonStereographic,help={"Make a Stereographic projection with the same orientation as the indexed pattern"}
//	Button buttonReportIndex,pos={29,205+30},size={160,20},proc=IndexButtonProc,title="Make Report Sheet"
//	Button buttonReportIndex,help={"if the top window is an images showing the hkl's, this makes a layout with reflections for printing"}
//
//	Button buttonCalcEnergyIndex,pos={29,235+30},size={160,20},proc=IndexButtonProc,title="Calc Energy of (hkl)"
//	Button buttonCalcEnergyIndex,help={"calculate the energy of an hkl using the current indexing"}
//	SetVariable h_IndexEneVar,pos={27,257+30},size={50,15},proc=CalcE_SetVarProc,title="H"
//	SetVariable h_IndexEneVar,value= root:Packages:micro:Index:h_e
//	SetVariable k_IndexEneVar,pos={91,257+30},size={50,15},proc=CalcE_SetVarProc,title="K"
//	SetVariable k_IndexEneVar,value= root:Packages:micro:Index:k_e
//	SetVariable L_IndexEneVar,pos={155,257+30},size={50,15},proc=CalcE_SetVarProc,title="L"
//	SetVariable L_IndexEneVar,value= root:Packages:micro:Index:l_e
//
//	PopupMenu popupTables,pos={62,290+30},size={66,20},proc=TablesPopMenuProc,title="Display Tables..."
//	PopupMenu popupTables,help={"Show table of identified peak positions indexed peaks, or the peaks for strain refinement"}
//	PopupMenu popupTables,fSize=14,mode=0,value= #"\"Fitted Peaks;Indexed Peaks;Strain Peaks;Missing Peaks;Measured Energies...\""
//
//	if (exists("IndexLots")==6)
//		SetDrawLayer UserBack
//		SetDrawEnv linethick= 2
//		DrawLine 0,330+30,221,330+30
//		Button buttonIndexLots,pos={29,340+30},size={160,20},proc=IndexButtonProc,title="Index Lots"
//		Button buttonIndexLots,help={"index lots of images"}
//		Button buttonWrite3dFile,pos={39,365+30},size={160,20},proc=IndexButtonProc,title="Write 3d Orients file"
//		Button buttonWrite3dFile,help={"Write 3d Orientations file"}
//		Button buttonFindClosestPoint,pos={39,390+30},size={160,20},proc=IndexButtonProc,title="find closest point"
//		Button buttonFindClosestPoint,help={"find point in original indexing list that is closest to the cursor"}
//		Button buttonIndexingListTable,pos={39,415+30},size={160,20},proc=IndexButtonProc,title="Indexing List Table"
//		Button buttonIndexingListTable,help={"Make a table of the indexing list"}
//	else
//		Button buttonInitIndexLots,pos={29,340+30},size={150,50},title="Init IndexLots\rpackage",proc=microGeo#LoadPackageButtonProc
//		Button buttonInitIndexLots,help={"Load Package for Indexing Lots of Images"}
//	endif
//
//	Indexing#EnableDisableIndexControls(hostWin+"#IndexPanel")
//	return "#IndexPanel"
//End
//
Static Function EnableDisableIndexControls(win)				// here to enable/disable
	String win												// window (or sub window) to use
	Variable d

	Button buttonLoadImage,win=$win,disable=0		// always OK to load an image

	d = strlen(WaveListClass("tif*;spe*;rawImage*;HDF*","*","DIMS:2"))<1 ? 2 : 0
	Button buttonViewImage,win=$win,disable= d

	d = strlen(WaveListClass("tif*;speImage;rawImage;HDFImage","*","DIMS:2"))<1 ? 2 : 0
	Button buttonBkgRemove,win=$win,disable= d

	d = strlen(WaveListClass("tif*;speImage*;rawImage*","*","DIMS:2"))<1 ? 2 : 0
	Button buttonFitPeaks,win=$win,disable= d
	Button buttonMakeMask,win=$win,disable= d

	d = strlen(WaveListClass("imageMask*","*","DIMS:2,BYTE:1,UNSIGNED:1"))<1 ? 2 : 0
	Button buttonEditMask,win=$win,disable= d

// want to get rid of the whole pkList thing.
//	d = strlen(StrVarOrDefault("pkList","")) ? 0 : 2
//	Button buttonPeaksFromBoxes,win=$win,disable= d

	STRUCT microGeometry g
	FillGeometryStructDefault(g)
	d = (strlen(WaveListClass("IndexedPeakList*","*","")) && g.Ndetectors>1) ? 0 : 2
	Button buttonAuxIndex,win=$win,disable= d

	d = strlen(WaveListClass("FittedPeakList*","*","MINCOLS:1"))<1 ? 2 : 0
	Button buttonIndex,win=$win,disable= d

	d = (strlen(WaveListClass("FittedPeakList*","*","DIMS:2,MAXCOLS:12,MINCOLS:11")) && strlen(WaveListClass("IndexedPeakList*","*",""))) ? 0 : 2
	Button buttonStrainRefine,win=$win,disable= d

	d = (strlen(WaveListClass("IndexedPeakList*","*",""))) ? 0 : 2
	Button buttonStereographic,win=$win,disable= d

	d = strlen(WaveListClass("IndexedPeakList*","*",""))<1 ? 2 : 0
	PopupMenu popuphklTags,win=$win,disable= d
	Button buttonReportIndex,win=$win,disable= d
	Button buttonCalcEnergyIndex,win=$win,disable= d
	if (d)
		SetVariable h_IndexEneVar,win=$win,noedit=1,frame=0,disable= d
		SetVariable k_IndexEneVar,win=$win,noedit=1,frame=0,disable= d
		SetVariable L_IndexEneVar,win=$win,noedit=1,frame=0,disable= d
	else
		SetVariable h_IndexEneVar,win=$win,noedit=0,frame=2,disable= d
		SetVariable k_IndexEneVar,win=$win,noedit=0,frame=2,disable= d
		SetVariable L_IndexEneVar,win=$win,noedit=0,frame=2,disable= d
	endif

	d = strlen(WaveListClass("FittedPeakList*;IndexedPeakList*","*",""))<1 ? 2 : 0
	PopupMenu popupTables,win=$win,disable=d

	if (exists("Load3dRecipLatticesFileXML")==6)
		FUNCREF ValidRawXMLdataAvailableProto fvalid = $"multiIndex#ValidRawXMLdataAvailable"
		String mstr="\"Load New 3D XML file;"+SelectString(fvalid(""),"","Load/Append 3D XML file;Process Loaded XML;")+"\""
		String title=SelectString(fvalid(""),"Load 3d XML","Process Loaded XML")
		PopupMenu popupIndexLoadXML,win=$win,disable=0,value=#mstr, title=title
		Button buttonIndexColorHexagon,win=$win,disable= 0

		d = strlen(WaveListClass("Random3dArraysGm","*",""))<1 ? 2 : 0
		Button buttonIndexNewPoleXML,win=$win,disable= d	
		Button buttonIndexRefineStrainMany,win=$win,disable= d	
		Button buttonIndexRefineStrain1,win=$win,disable= d	

		d = strlen(WaveListClass("poleXYpoints","*",""))<1 ? 2 : 0
		Button buttonIndexRePlotPoleXML,win=$win,disable= d	

		d = strlen(WaveListClass("Random3dArrays","*",""))<1 ? 2 : 0
		Button buttonIndexLoadOneXML,win=$win,disable= d
		d = strlen(WaveListClass("Random3dArraysXYZ","*",""))<1 ? 2 : 0
		Button buttonIndexGizmoXML,win=$win,disable= d	

//		d = strlen(WaveListClass("Interpolated3dArraysXYZ","*",""))<1 ? 2 : 0
		d = strlen(WaveListClass("Random3dArrays*","*",""))<1 ? 2 : 0
		Button buttonIndex2D3DXMLrgba,win=$win,disable= d
		d = strlen(WaveListClass("Random3dArrays","*",""))<1 ? 2 : 0
		Button buttonIndex2Dplot,win=$win,disable= d
//	MenuItemIfWaveClassExists("2D plot of loaded XML data","Random3dArrays",""), Make2Dplot_xmlData("")

		d = strlen(WaveListClass("Interpolated3dArrays","*",""))<1 ? 2 : 0
		Button buttonIndexHistogram,win=$win,disable= d
		d = strlen(WaveListClass("HistogramFrom3d","*",""))<1 ? 2 : 0
		Button buttonIndexHistogramPlot,win=$win,disable= d

		d = strlen(WaveListClass("Random3dArrays","*",""))<1 ? 2 : 0
		Button buttonIndexRotation,win=$win,disable= d
	endif

	if (exists("IndexLots")==6)
		Button buttonIndexLots,win=$win,disable= 0					// always OK to do this

		d = strlen(WaveListClass("indexationResultList","*","TEXT:1"))<1 ? 2 : 0
		Button buttonWrite3dFile,win=$win,disable= d

		d = strlen(WaveListClass("OrientationSliceWave*","*",""))<1 ? 2 : 0
		Button buttonFindClosestPoint,win=$win,disable= d

		d = strlen(WaveListClass("indexationResultList*","*","TEXT:1"))<1 ? 2 : 0
		Button buttonIndexingListTable,win=$win,disable= d
	endif
End
Function ValidRawXMLdataAvailableProto(fldr)
	String fldr
	return 0
End


Function IndexButtonProc(B_Struct) : ButtonControl
	STRUCT WMButtonAction &B_Struct
	if (B_Struct.eventCode != 2)
		return 0
	endif
	String ctrlName=B_Struct.ctrlName

	STRUCT microGeometry g
	FillGeometryStructDefault(g)

	if (stringmatch(ctrlName,"buttonLoadImage"))
		Wave image = $LoadGenericImageFile("")
		if (WaveExists(image))
			String wnote=note(image), bkgFile=StringByKey("bkgFile",wnote,"=")
			Variable xdim=NumberByKey("xdim",wnote,"="), ydim=NumberByKey("ydim",wnote,"=")
			printf "for file '"+GetWavesDataFolder(image,2)+"'"
			if (strlen(bkgFile)>0)
				printf ",             background file = '%s'",  bkgFile
			endif
			printf "\r"
			printf "total length = %d x %d  = %d points\r", xdim,ydim,xdim*ydim
			print "number type is  '"+IgorFileTypeString(NumberByKey("NUMTYPE",WaveInfo(image,0)))+"'"
			print "Created a 2-d wave    '"+GetWavesDataFolder(image,2)+"'"
			NewImageGraphLocal(image)
		endif
	elseif (stringmatch(ctrlName,"buttonViewImage") && strlen(WaveListClass("spe*;rawImage*;HDF*","*","DIMS:2")))
		NewImageGraphLocal($"")
//		NewImageGraph($"",withButtons=1)
//		Graph_imageMake($"",1)
	elseif (stringmatch(ctrlName,"buttonMakeMask") && strlen(WaveListClass("speImage;rawImage*","*","DIMS:2")))
		MakeMaskThreshold($"",NaN,printIt=1)
	elseif (stringmatch(ctrlName,"buttonEditMask") && strlen(WaveListClass("imageMask*","*","DIMS:2,BYTE:1,UNSIGNED:1")))
		MakeMask($"")
	elseif (stringmatch(ctrlName,"buttonBkgRemove") && strlen(WaveListClass("speImage;rawImage","*","DIMS:2")))
		RemoveBackgroundImage($"",NaN,printIt=1)
	elseif (stringmatch(ctrlName,"buttonFitPeaks") && strlen(WaveListClass("speImage*;rawImage*","*","DIMS:2")))
		// FitPeaks($"",NaN,NaN)
		//		FunctionList("FitPeak*",";" , "WIN:[Indexing]")
		if (exists("root:Packages:micro:PeakFit:fitPeakFuncLast")!=2)
			indexing#initPeaksPart()
		endif
		SVAR fitPeakFuncLast = root:Packages:micro:PeakFit:fitPeakFuncLast
		String fitPeakFunc=fitPeakFuncLast
		Prompt fitPeakFunc,"peak fitting method",popup,"FitPeaksWithExternal;FitPeaksWithSeedFill;FitPeaks;FitPeaksNew;FitPeaksStepWise"
		DoPrompt "Peak Fit",fitPeakFunc
		if (!V_flag)
			if (stringmatch(fitPeakFunc,"FitPeaks"))
				FitPeaks($"",NaN,NaN,printIt=1)
			elseif (stringmatch(fitPeakFunc,"FitPeaksNew"))
				FitPeaksNew($"",NaN,NaN,printIt=1)
			elseif (stringmatch(fitPeakFunc,"FitPeaksStepWise"))
				FitPeaksStepWise($"",NaN,NaN,NaN,NaN,printIt=1)
			elseif (stringmatch(fitPeakFunc,"FitPeaksWithSeedFill"))
				FitPeaksWithSeedFill($"",NaN,NaN,NaN,NaN,printIt=1)
			elseif (stringmatch(fitPeakFunc,"FitPeaksWithExternal"))
				FitPeaksWithExternal($"",NaN,NaN,NaN,NaN,printIt=1)
			endif
			fitPeakFuncLast = fitPeakFunc
		endif
	elseif (stringmatch(ctrlName,"buttonAuxIndex") && strlen(WaveListClass("IndexedPeakList*;","*","")) && g.Ndetectors>1)
		MakeIndexedWaveForAuxDetector(NaN,$"",$"")
//	elseif (stringmatch(ctrlName,"buttonPeaksFromBoxes") && strlen(StrVarOrDefault("pkList","")))
//		SVAR pkList=pkList
//		setFittedPeaksFromList($"",NaN,pkList)
	elseif (stringmatch(ctrlName,"buttonIndex") && strlen(WaveListClass("FittedPeakList*","*","DIMS:2,MAXCOLS:12,MINCOLS:11")))
		IndexAndDisplay($"",NaN,NaN,NaN,NaN,NaN,NaN,NaN,printit=1)
	elseif (stringmatch(ctrlName,"buttonStrainRefine"))
#ifdef USE_ENERGY_STRAIN_REFINE			//#if (Exists("TotalStrainRefine")==6)
		TotalStrainRefine(NaN,"")
#else
		DeviatoricStrainRefine(NaN,"",printit=1)
#endif
	elseif (stringmatch(ctrlName,"buttonStereographic"))
		StereoOfIndexedPattern($"",NaN,printIt=1)
		if (exists("MakeStereo")!=6)
			Button buttonStereographic,win=$(B_Struct.win),title="Stereographic..."
		endif
	elseif (stringmatch(ctrlName,"buttonReportIndex"))
		MakeIndexingReportSheet()
	elseif (stringmatch(ctrlName,"buttonCalcEnergyIndex"))
		NVAR h=root:Packages:micro:Index:h_e
		NVAR k=root:Packages:micro:Index:k_e
		NVAR l=root:Packages:micro:Index:l_e
		Variable/G Energy_keV = EnergyOfhkl($"",NaN,h,k,l,printIt=1)	// the energy is saved in a global so that it is available by the user

#if (Exists("Load3dRecipLatticesFileXML")==6)
	elseif (stringmatch(ctrlName,"buttonIndexNewPoleXML") && strlen(WaveListClass("Random3dArraysGm","*","")))
		NewPoleFigure("")
	elseif (stringmatch(ctrlName,"buttonIndexRePlotPoleXML") && strlen(WaveListClass("poleXYpoints","*","")))
		Make_GraphPoles($"")
	elseif (stringmatch(ctrlName,"buttonIndexRefineStrainMany") && strlen(WaveListClass("Random3dArraysGm","*","")))
		DeviatoricStrainRefineXML_ALL("","")
	elseif (stringmatch(ctrlName,"buttonIndexRefineStrain1") && strlen(WaveListClass("Random3dArraysGm","*","")))
		DeviatoricStrainRefineXML(-1,NaN,"",printit=1)
	elseif (stringmatch(ctrlName,"buttonIndexLoadOneXML") && strlen(WaveListClass("Random3dArrays","*","")))
		LoadPixelsFromXML(-1)
	elseif (stringmatch(ctrlName,"buttonIndex2D3DXMLrgba") && strlen(WaveListClass("Random3dArraysXYZ,Random3dArrays","*","")))
		Make2D_3D_RGBA($"",$"","",NaN,NaN)
	elseif (stringmatch(ctrlName,"buttonIndex2Dplot") && strlen(WaveListClass("Random3dArraysXYZ,Random3dArrays","*","")))
		Make2Dplot_xmlData("")
	elseif (stringmatch(ctrlName,"buttonIndexColorHexagon"))
		MakColorHexagon()
	elseif (stringmatch(ctrlName,"buttonIndexGizmoXML") && strlen(WaveListClass("Random3dArrays","*","")))
		MakeGizmo_xmlData($"")
	elseif (stringmatch(ctrlName,"buttonIndexHistogram") && strlen(WaveListClass("Interpolated3dArrays","*","")))
		multiIndex#HistogramOf3dArray($"")
	elseif (stringmatch(ctrlName,"buttonIndexHistogramPlot") && strlen(WaveListClass("HistogramFrom3d","*","")))
		multiIndex#DisplayHistogramFrom3d($"")
	elseif (stringmatch(ctrlName,"buttonIndexRotation") && strlen(WaveListClass("Random3dArrays","*","")))
		RotationBetweenTwoPoints(NaN,NaN)
#endif

#if (Exists("IndexLots")==6)
	elseif (stringmatch(ctrlName,"buttonIndexLots") )
		IndexLots("","","","",NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN)
	elseif (stringmatch(ctrlName,"buttonWrite3dFile") && strlen(WaveListClass("indexationResultList","*","TEXT:1")))
		Write3dOrientsFile($"","")
	elseif (stringmatch(ctrlName,"buttonFindClosestPoint") && strlen(WaveListClass("OrientationSliceWave*","*","")))
		findClosestPoint(NaN,NaN)
	elseif (stringmatch(ctrlName,"buttonIndexingListTable") && strlen(WaveListClass("indexationResultList*","*","TEXT:1")))
		TableOfIndexingList($"")
#endif
	endif
	Indexing#EnableDisableIndexControls(GetUserData("microPanel","","IndexPanelName"))
End
//
Function hklTagsPopMenuProc(ctrlName,popNum,popStr) : PopupMenuControl
	String ctrlName
	Variable popNum
	String popStr
	if (strlen(WaveListClass("IndexedPeakList*","*",""))<1)
		return 1
	endif
	if (stringmatch(popStr,"Re-Draw hkl tags"))
		DisplayResultOfIndexing($"",NaN)
	elseif (stringmatch(popStr,"Clear hkl tags off Graph"))
		ClearhklOffGraph("")
	elseif (stringmatch(popStr,"Add Missing Peaks"))
		AddMissingReflections($"")
	elseif (stringmatch(popStr,"Remove Missing Peaks"))
		RemoveMissingPeaksFromGraph("")
	elseif (stringmatch(popStr,"Rotation Between Orientations..."))
		RotationBetweenTwoIndexations($"",NaN,$"",NaN)
	elseif (stringmatch(popStr,"Add Peaks from a user recip"))
		AddReflectionsFromUserRecip($"")
	elseif (stringmatch(popStr,"Remove All Peaks from user recips"))
		RemoveUserRecipFromGraph("")
	endif
	Indexing#EnableDisableIndexControls(GetUserData("microPanel","","IndexPanelName"))
End
//
Function CalcE_SetVarProc(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName
	if (strlen(WaveListClass("IndexedPeakList*","*",""))<1)
		return 1
	endif
	NVAR h=root:Packages:micro:Index:h_e
	NVAR k=root:Packages:micro:Index:k_e
	NVAR l=root:Packages:micro:Index:l_e
	// EnergyOfhkl($"",NaN,h,k,l)
	Indexing#EnableDisableIndexControls(GetUserData("microPanel","","IndexPanelName"))
End
//
Function TablesPopMenuProc(ctrlName,popNum,popStr) : PopupMenuControl
	String ctrlName
	Variable popNum
	String popStr
	if (strlen(WaveListClass("FittedPeakList*;IndexedPeakList*","*",""))<1)
		return 1
	endif
	if (stringmatch(popStr,"Fitted Peaks"))
		DisplayTableOfWave($"",classes="FittedPeakList*",promptStr="Wave with list of fitted peaks")
	elseif (stringmatch(popStr,"Indexed Peaks"))
		DisplayTableOfWave($"",classes="IndexedPeakList*",promptStr="Wave with list of indexed peaks")
	elseif (stringmatch(popStr,"Strain Peaks"))
		DisplayTableOfWave($"",classes="StrainPeakList*",promptStr="Wave with list of Strain Refined peaks")
	elseif (stringmatch(popStr,"Missing Peaks"))
		DisplayTableOfWave($"",classes="MissingPeakList*",promptStr="Wave with list of missing peaks")
	elseif (stringmatch(popStr,"Measured Energies*"))
		MakeMeasured_hkl_EnergiesWave(NaN,"")
	elseif (stringmatch(popStr,"Peaks from User recip"))
		DisplayTableOfWave($"",classes="UserRecipPeakList*",promptStr="Wave with list of User hkl peaks")
	endif
End
//
Static Function IndexXMLPopUpMenuProc(ctrlName,popNum,popStr) : PopupMenuControl
	String ctrlName
	Variable popNum
	String popStr
	Variable err=1		// when no error, also re-process the loaded data
#if (Exists("Load3dRecipLatticesFileXML")==6)
	if (stringmatch(popStr,"Load New 3D XML file"))
		err = Load3dRecipLatticesFileXML("")
	elseif (stringmatch(popStr,"Load/Append 3D XML file"))
		err = Load3dRecipLatticesFileXML("", appendXML=1)
	elseif (stringmatch(popStr,"Process Loaded XML"))
		err = 0
	endif
	if (!err)
		ProcessLoadedXMLfile(Inf,NaN)
	endif
#endif
	Indexing#EnableDisableIndexControls(GetUserData("microPanel","","IndexPanelName"))
End

Static Function/S NewImageGraphLocal(image,[withButtons])
	Wave image
	Variable withButtons
	withButtons = ParamIsDefault(withButtons) ? NaN : withButtons
	withButtons = numtype(withButtons) ? 1 : !(!withButtons)

	String win=StringFromList(0,WindowsWithWave(image,1))
	if (strlen(win))
		DoWindow/F $win
		return GetWavesDataFolder(image,2)
	endif

	String result = NewImageGraph(image,withButtons=withButtons)
	Wave image = $result
	if (!WaveExists(image))
		return ""
	elseif (strlen(GetUserData("","","pixelCenter")))
		return result								// no need to add the cross, it is already there
	endif

	// the following puts a cross where the detector is direction above the sample or on the direct beam
	String wnote=note(image)					// get start & group in case image is binned
	Variable startx=NumberByKey("startx",wnote,"="), starty=NumberByKey("starty",wnote,"=")
	Variable endx=NumberByKey("endx",wnote,"="), endy=NumberByKey("endy",wnote,"=")
	Variable groupx=NumberByKey("groupx",wnote,"="), groupy=NumberByKey("groupy",wnote,"=")

	STRUCT microGeometry g
	FillGeometryStructDefault(g)				//fill the geometry structure with current values
	Variable dnum = detectorNumFromID(g, StringByKey("detectorID",note(image),"="))
	if (dnum<0)
		return result
	endif

	String axes=microGeo#AxesHitDetector(g.d[dnum])		// returns pixel where X, Y, Z, H, or F axes hit detector
	String list, winData=""
	Variable i, px,py
	for (i=0;i<ItemsInList(axes);i+=1)						// for each axis that hits the detector, draw a cross there
		list = StringFromList(i,axes)							// probably only contain 1 item
		px = NumberByKey("px",list,"=",",")
		py = NumberByKey("py",list,"=",",")
		if (startx<=px && px<=endx && starty<=py && py<=endy)	// (px,py) IS on this ROI image of detector, draw cross
			Variable NxFW=(g.d[dnum].Nx)/10, NyFW=(g.d[dnum].Ny)/10						// size of cross
			if (numtype(startx+starty+groupx+groupy)==0 && (groupx>1 || groupy>1))	// a binned image, rescale to binned pixels
				px = round(( px-startx-(groupx-1)/2 )/groupx)	// pixel is zero based here & startx is zero based
				py = round(( py-starty-(groupy-1)/2 )/groupy)	// groupx=1 is un-binned
				NxFW /= groupx
				NyFW /= groupy
			endif
			DrawMarker(px,py,round(NxFW),round(NyFW),"cross gap",dash=2,layer="UserAxes")
			winData += list + ";"
		endif
	endfor
	SetWindow kwTopWin, userdata(pixelCenter)=winData
	return result
End

//	End of Index set panel
// =========================================================================
// =========================================================================



// =========================================================================
// =========================================================================
//	Start of compatibility section

//Function/S LoadGenericImageFile(fileName)
//	String fileName
//	String str=""
//
//#if (Exists("LoadWinViewFile")==6)
//	str = LoadWinViewFile("")
//#elif (Exists("LoadHDF5imageFile")==6)
//	str = LoadHDF5imageFile(fileName)
//#endif
//	return str
//End

//// Get rid of this when HDF is available
//Function/S pkLIstImageMenuItem(item)
//	String item
//	return pkLIstWinViewMenuItem(item)
//End

//	End of compatibility section
// =========================================================================
// =========================================================================



Function InitIndexingPackage()					// used to initialize this package
	InitLatticeSymPackage()					// init Lattice Sym package
	initImageDisplayScaling()
	NewDataFolder/O root:Packages:geometry
	initPeaksPart()
End

Static Function initPeaksPart()
	// presets for peaksearching
	NewDataFolder/O root:Packages:micro
	NewDataFolder/O root:Packages:micro:PeakFit

	if (exists("root:Packages:micro:PeakFit:fitPeakFuncLast")!=2)
		String/G root:Packages:micro:PeakFit:fitPeakFuncLast = ""
	endif

	if (exists("root:Packages:micro:PeakFit:minPeakWidthLast")!=2)
		Variable/G root:Packages:micro:PeakFit:minPeakWidthLast = 1.13
	endif
	if (exists("root:Packages:micro:PeakFit:boxSizeLast")!=2)
		Variable/G root:Packages:micro:PeakFit:boxSizeLast = 18
	endif
	if (exists("root:Packages:micro:PeakFit:maxRfactorLast")!=2)
		Variable/G root:Packages:micro:PeakFit:maxRfactorLast = 0.5
	endif
	if (exists("root:Packages:micro:PeakFit:threshAboveAvgLast")!=2)
		Variable/G root:Packages:micro:PeakFit:threshAboveAvgLast = NaN
	endif
	if (exists("root:Packages:micro:PeakFit:peakShapeLast")!=2)
		String/G root:Packages:micro:PeakFit:peakShapeLast = "Lorentzian"
	endif
	if (exists("root:Packages:micro:PeakFit:maskNameLast")!=2)
		String/G root:Packages:micro:PeakFit:maskNameLast = ""
		String/G root:Packages:micro:PeakFit:maskHashLast = ""
	endif
	if (exists("root:Packages:micro:PeakFit:maskHashLast")!=2)
		String/G root:Packages:micro:PeakFit:maskHashLast = ""
	endif

	if (exists("root:Packages:micro:PeakFit:maxPeakWidthLast")!=2)
		Variable/G root:Packages:micro:PeakFit:maxPeakWidthLast = 70
	endif
	if (exists("root:Packages:micro:PeakFit:minSepLast")!=2)
		Variable/G root:Packages:micro:PeakFit:minSepLast = 50
	endif
	if (exists("root:Packages:micro:PeakFit:smoothingLast")!=2)
		Variable/G root:Packages:micro:PeakFit:smoothingLast = 0
	endif
	if (exists("root:Packages:micro:PeakFit:thresholdRatioLast")!=2)
		Variable/G root:Packages:micro:PeakFit:thresholdRatioLast = NaN
	endif
	if (exists("root:Packages:micro:PeakFit:maxNumLast")!=2)
		Variable/G root:Packages:micro:PeakFit:maxNumLast = Inf
	endif
End


Function Itest()
		STRUCT microGeometry geo
		FillGeometryStructDefault(geo)
#if (Exists("MICRO_VERSION_ORIGINAL_Func")==6)
		print "geo.ddOffset = ",geo.ddOffset
#else
		print "Ndetectors =",geo.Ndetectors
#endif

#if (Exists("MICRO_VERSION_N_Func")==6)
		print "Ndetectors =",geo.Ndetectors
#else
		print "geo.ddOffset = ",geo.ddOffset
#endif
End

// change colons in path to backslash
Function/S MacToWindows(MacPath)
	String MacPath
	String WinPath=""
	variable FirstColon = 0
	variable i
	for(i=0; i<strlen(MacPath); i += 1)
		if ( char2num(MacPath[i]) == char2num(":") )
			if (FirstColon == 0)
				// RX 03.17.2012 - there's no drive letter for a network path on Windows, e.g. "\\computer1\folder1\something\"
				if((char2num(MacPath[0]) == char2num("\\")) && (char2num(MacPath[1]) == char2num("\\"))) 
					WinPath = WinPath + "\\"
				else
					WinPath = WinPath + ":\\"		// The first colon is the drive letter, e.g.. c:\
				endif
				FirstColon = 1
			else
				WinPath = WinPath + "\\"		// The subsequent colons separate folders in the tree
			endif
		else
			WinPath = WinPath + MacPath[i]
		endif
	endfor
	return WinPath
End

Static Function/S getBinaryPath()
	// returns path to the binaries.
	// looks for folder inside of folder containing runIndexingEulerCommand()
	// if "/binary" does not exist, just return path to runIndexingEulerCommand()
	String path = ParseFilePath(1,FunctionPath("runIndexingEulerCommand"),":",1,0)// default path to the executables
	GetFileFolderInfo/Q/Z=1 path+"binary:"
	path = SelectString(V_isFolder && !V_Flag, path, S_Path)
	return path
End

