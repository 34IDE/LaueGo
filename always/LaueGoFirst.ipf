#pragma rtGlobals= 2
#pragma version = 2.09
#pragma ModuleName = LaueGoFirst
#pragma hide = 1
//Static Constant JZTalwaysFirst_Version_Min=2.4	// minimum required vesion of "always first.ipf"

#if NumVarOrDefault("root:Packages:SHOW_PACKAGES_MASK_JZT",-1) & 4		// LaueGo
StrConstant MICRO_GEOMETRY_VERSION_PATH = "root:Packages:MICRO_GEOMETRY_VERSION"
//		bit flags for MICRO_GEOMETRY_VERSION
//			0 -> OLD
//			1 -> HDF5
//			2 -> TIFF
//			4 -> SPE
#endif


#if NumVarOrDefault("root:Packages:SHOW_PACKAGES_MASK_JZT",-1) & 8		// Scattering
Menu "Analysis"
	Submenu "Packages"
		"-"
		SubMenu "Scattering..."
			"X-ray Data",Execute/P "INSERTINCLUDE  \"Xray\", version>=2.21";Execute/P "COMPILEPROCEDURES ";Execute/P "ElementDataInitPackage()"
			help = {"Load procedures for providing X-ray data"}
			"  Elements Data",Execute/P "INSERTINCLUDE  \"Elements\", version>=1.69";Execute/P "COMPILEPROCEDURES ";Execute/P "ElementDataInitPackage()"
			help = {"Load procedures for providing data on the elements"}
			"  Cromer-Liberman",Execute/P "INSERTINCLUDE  \"CromerLiberman\", version>=1.7";Execute/P "COMPILEPROCEDURES "
			help = {"Load procedures for providing Cromer-Liberman"}
			"Ion Chamber",Execute/P "INSERTINCLUDE  \"IonChamber\", version>=3.2";Execute/P "COMPILEPROCEDURES ";Execute/P "ionChamberInitPackage()"
			help = {"Load procedures for evaluating ion chamber output"}
			"Lattices",Execute/P "INSERTINCLUDE  \"LatticeSym\", version>=3.77";Execute/P "COMPILEPROCEDURES ";Execute/P "InitLatticeSymPackage(showPanel=1)"
			help = {"Load lattice symmetry procedures"}
			"-"
			"Neutron Scattering Data",Execute/P "INSERTINCLUDE  \"Neutron\", version>=1.1";Execute/P "COMPILEPROCEDURES ";Execute/P "NeutronDataInitPackage()"
			help = {"Load procedures for providing Neutron Scattering Data"}
		End
	End
End
#endif

#if NumVarOrDefault("root:Packages:SHOW_PACKAGES_MASK_JZT",-1) & 16	// APS Specific
Menu "Analysis"
	Submenu "Packages"
		"(APS specific"
		Submenu "spec"
			"spec with images",Execute/P "INSERTINCLUDE  \"specImages\", version>=0.44";Execute/P "COMPILEPROCEDURES ";Execute/P "init_specImage(\"\")"
			help = {"Support for spec files plus support for spec with image files (like a Pilatus)."}
			"spec files",Execute/P "INSERTINCLUDE  \"spec\", version>=2.27";Execute/P "COMPILEPROCEDURES ";Execute/P "specInitPackage()"
			help = {"Load everything for reading spec files."}
		End
		"MDA files",Execute/P "INSERTINCLUDE  \"mdaFiles\"";Execute/P "COMPILEPROCEDURES "
		help = {"Load mda files (from APS)."}
		SubMenu "EPICS"
			"EPICS("
			"PV I\O", Execute/P "INSERTINCLUDE \"epics\"" ; 	Execute/P "COMPILEPROCEDURES ";Execute/P "epicsInitPackage()"
			help = {"Load procedures for talking to PVs via EPICS (only useful at APS)"}
			"Load Scan Record", Execute/P "INSERTINCLUDE \"LoadEPICSscans\"" ; 	Execute/P "COMPILEPROCEDURES "
			help = {"Load procedures for Loading the dump of a scan record using TkGrab, does not need EPICS support"}
		End
		"BURT Files",Execute/P "INSERTINCLUDE  \"BurtFiles\"";Execute/P "COMPILEPROCEDURES "
		help = {"Load Burt files (from APS)."}
		"WinView Reader", Execute/P "INSERTINCLUDE \"WinView\", version>=2.01"
		help = {"Load procedures reading and looking at WinView images"}
	End
End
#endif

#if NumVarOrDefault("root:Packages:SHOW_PACKAGES_MASK_JZT",-1) & 4		// LaueGo
//Strconstant browserLeadJZT="\r--> moved to ", browserTailJZT="\r\r"
//ModifyBrowser echoCommands=0, command3="print browserLeadJZT,GetDataFolder(1),browserTailJZT"
Menu "Analysis",dynamic
	Submenu "Packages"
		"-"
		SubMenu "<B<ULaueGo [µ beam]"

			// none
			microMenuShowPreliminary(LaueGoFirst#arrow("B")+"Select Image Types"), LaueGoFirst#SetMicroPreferedImageTypes()
			help = {"Select image types for LaueGo"}

			// NEW
			microMenuShowN(LaueGoFirst#arrow("")+"LaueGo Panel"),MakeMicroPanel(-1)
			help = {"new starting point for the generic micro beam analysis (LaueGo"}
			microMenuShowN("Indexing & Peak Searching"),Execute/P "INSERTINCLUDE  \"IndexingN\", version>=4.00";Execute/P "COMPILEPROCEDURES ";Execute/P "InitIndexingPackage()"
			//	help = {"Find and Fit Peaks in an image, and then Index the peaks"}
			microMenuShowN("Energy-Wire Scans"),Execute/P "INSERTINCLUDE  \"EnergyWireScansN\", version>=1.000";Execute/P "COMPILEPROCEDURES "
			//	help = {"Processes Energy-Wire scans into Depth-Q plots"}
			microMenuShowN("Index Lots of Images"),Execute/P "INSERTINCLUDE  \"IndexLotsN\", version>=2.27";Execute/P "COMPILEPROCEDURES "
			//	help = {"index a large group of images, like from a wire-scan"}
			microMenuShowN("3d-Grains with Gizmo"),Execute/P "INSERTINCLUDE  \"GizmoGrains\", version>=1.3";Execute/P "COMPILEPROCEDURES ";Execute/P "initGrainGizmo()"
			//	help = {"Processes 3D array of matricies into a viewable volume"}
			microMenuShowN("3d-Array of Orientation Mats"),Execute/P "INSERTINCLUDE  \"ArrayOf3dOrientsN\", version>=2.56";Execute/P "COMPILEPROCEDURES "
			//	help = {"Processes 3D array of matricies into viewable of rotataion angles, dislocation tensors, GND, etc."}
			microMenuShowN("Pole Figure"),Execute/P "INSERTINCLUDE  \"PoleFigure\"";Execute/P "COMPILEPROCEDURES "
			//	help = {"Loads Pole figure Macros"}
			microMenuShowN("Details of Depth Resolved Wire Scans"),Execute/P "INSERTINCLUDE  \"DepthResolvedQueryN\"";Execute/P "COMPILEPROCEDURES "
			//	help = {"used to analyze a depth resolved wire scan, mostly for debugging"}

			// OLD
			microMenuShowOLD("ÐÐ> LaueGo Panel"),MakeMicroPanel(-1)
			help = {"old starting point for the generic micro beam analysis"}
			microMenuShowOLD("Indexing & Peak Searching"),Execute/P "INSERTINCLUDE  \"Indexing\", version>=2.59";Execute/P "COMPILEPROCEDURES ";Execute/P "InitIndexingPackage()"
			//	help = {"Find and Fit Peaks in an image, and then Index the peaks"}
			microMenuShowOLD("Energy-Wire Scans"),Execute/P "INSERTINCLUDE  \"EnergyWireScans\", version>=0.9944";Execute/P "COMPILEPROCEDURES "
			//	help = {"Processes Energy-Wire scans into Depth-Q plots"}
			microMenuShowOLD("Index Lots of Images"),Execute/P "INSERTINCLUDE  \"IndexLots\", version>=2.17";Execute/P "COMPILEPROCEDURES "
			//	help = {"index a large group of images, like from a wire-scan"}
			microMenuShowOLD("3d-Grains with Gizmo"),Execute/P "INSERTINCLUDE  \"GizmoGrains\", version>=1.3";Execute/P "COMPILEPROCEDURES ";Execute/P "initGrainGizmo()"
			//	help = {"Processes 3D array of matricies into a viewable volume"}
			microMenuShowOLD("3d-Array of Orientation Mats"),Execute/P "INSERTINCLUDE  \"ArrayOf3dOrients\", version>=2.41";Execute/P "COMPILEPROCEDURES "
			//	help = {"Processes 3D array of matricies into viewable of rotataion angles, dislocation tensors, GND, etc."}
			microMenuShowOLD("Pole Figure"),Execute/P "INSERTINCLUDE  \"PoleFigure\"";Execute/P "COMPILEPROCEDURES "
			//	help = {"Loads Pole figure Macros"}
			microMenuShowOLD("Details of Depth Resolved Wire Scans"),Execute/P "INSERTINCLUDE  \"DepthResolvedQuery\"";Execute/P "COMPILEPROCEDURES "
			//	help = {"used to analyze a depth resolved wire scan, mostly for debugging"}
			microMenuShowN("Laue Simulation"),Execute/P "INSERTINCLUDE  \"LaueSimulation\", version>=1.06";Execute/P "COMPILEPROCEDURES "

			Submenu "Calibration Procedures"
				microMenuShowN("Detector Calibration"),Execute/P "INSERTINCLUDE  \"DetectorCalibration\"";Execute/P "COMPILEPROCEDURES "
				help = {"Used for the Initial Calibration of the Perkin-Elmer detectors"}
				"micro-mono Calibration",Execute/P "INSERTINCLUDE  \"monoCalibrate\", version>=1.6";Execute/P "COMPILEPROCEDURES ";Execute/P "monoCalibrateInitPackage()"
				help = {"Used to calibrate energy of mono at 34ID-D"}
			End

		End
	End
End
#endif


#if NumVarOrDefault("root:Packages:SHOW_PACKAGES_MASK_JZT",-1) & 8		// Scattering
Menu "Analysis"
	Submenu "Packages"
		// "µ beam, Sector 34" SubMenu gets insterted below
		"Stereographic Projections",Execute/P "INSERTINCLUDE  \"StereographicProjection\", version>=2.81";Execute/P "COMPILEPROCEDURES ";Execute/P "InitStereoGraphicPackage()"
	End
End
#endif


#if NumVarOrDefault("root:Packages:SHOW_PACKAGES_MASK_JZT",-1) & 4		// LaueGo
Function/T microMenuShowN(str)
	String str
	return SelectString(NumVarOrDefault(MICRO_GEOMETRY_VERSION_PATH,256)&5,"",str)		// for both spe AND hdf5
End
Function/T microMenuShowOLD(str)
	String str
	return SelectString(NumVarOrDefault(MICRO_GEOMETRY_VERSION_PATH,NaN)==0,"",str)
End
Function/T microMenuShowPreliminary(str)
	String str
	return SelectString(exists(MICRO_GEOMETRY_VERSION_PATH)!=2,"",str)
End


Static Structure microImageTypePrefs
	int16	old
	int16	hdf5
	int16	tiff
	int16	spe
EndStructure

Static Function SetMicroPreferedImageTypes()
	STRUCT microImageTypePrefs prefs
	prefs.old = 0
	prefs.hdf5 = 0
	prefs.tiff = 0
	prefs.spe = 0
	LoadPackagePreferences/MIS=1 "microGeo","microGeoNPrefs",1,prefs
	if (V_bytesRead<8)
		prefs.old = 0
		prefs.hdf5 = 0
		prefs.tiff = 0
		prefs.spe = 0
	endif
	Variable old,hdf5,tiff,spe, repeat
	Prompt old,"Old Style",popup,"---;Old Style, for OLD data"
	Prompt hdf5,"HDF5 images",popup,"---;HDF5 Images"
	Prompt tiff,"TIFF images",popup,"---;TIFF Images in 4500S"
	Prompt spe,"WinView images",popup,"---;WinView (spe) Images"
	//DoAlert 0,"Select either 'Old Style' or some combination of the others"
	do
		old = prefs.old + 1
		hdf5 = prefs.hdf5 + 1
		tiff = prefs.tiff + 1
		spe = prefs.spe + 1
		DoPrompt "Desired Image Types",old,hdf5,tiff,spe
		if (V_flag)
			return -1
		endif
		prefs.old = old-1
		prefs.hdf5 = hdf5-1
		prefs.tiff = tiff-1
		prefs.spe = spe-1
		repeat = 0
		if ((prefs.hdf5 + prefs.tiff + prefs.spe) && prefs.old)
			DoAlert 0,"You cannot select 'Old Style' and anything else"
			repeat = 1
		elseif (!(prefs.hdf5 + prefs.tiff + prefs.spe) && !prefs.old)
			DoAlert 0, "You did not select any image type"
			repeat = 1
		endif
	while (repeat)
	SavePackagePreferences/FLSH=1 "microGeo","microGeoNPrefs",1,prefs
	setMICRO_GEOMETRY_VERSION_PATH()
End

Static Function setMICRO_GEOMETRY_VERSION_PATH([preset])
	Variable preset						// Optional, a bit flag 0=OLD, 1=HDF5, 2=Tiff, 3=SPE, if old, then nothing else allowed
	if (exists(MICRO_GEOMETRY_VERSION_PATH)==2)
		return 0
	endif

	STRUCT microImageTypePrefs prefs
	if (ParamIsDefault(preset))				// no preset given, use stored default values
		LoadPackagePreferences/MIS=1 "microGeo","microGeoNPrefs",1,prefs
		if (V_bytesRead<8)
			DoAlert 0,"Unable to read prefered images types, defaulting to only HDF5"
			prefs.old = 0
			prefs.hdf5 = 1
			prefs.tiff = 0
			prefs.spe = 0
		endif
	else
		prefs.old = (preset == 0)				// old cannot coesist with any of the others
		prefs.hdf5 = (preset & 1)
		prefs.tiff = (preset & 2)
		prefs.spe = (preset & 4)
//		SavePackagePreferences/FLSH=1 "microGeo","microGeoNPrefs",1,prefs	// save the defaults
	endif

	NewDataFolder/O root:Packages
	Variable/G $MICRO_GEOMETRY_VERSION_PATH
	NVAR type=$MICRO_GEOMETRY_VERSION_PATH
	if (prefs.old)
		type = 0
	else
		type += prefs.hdf5 ? 1 : 0
		type += prefs.tiff ? 2 : 0
		type += prefs.spe ? 4 : 0
	endif
	Variable/G $MICRO_GEOMETRY_VERSION_PATH=type
	String insert

	if (type==0)
		insert = "INSERTINCLUDE  \"microGeometry\", version>=2.61"
	elseif (type&2)
		insert = "INSERTINCLUDE  \"microGeometryN\", version>=1.61"
	elseif (type&5)
		insert = "INSERTINCLUDE  \"microGeometryN\", version>=1.61"
	else
		return -1
	endif
	Execute/P insert
	Execute/P "COMPILEPROCEDURES "
	Execute/P/Q "Init_microGeo()"
	Execute/P/Q "MakeMicroPanel(-1)"
	return 0
End


Static Function/T arrow(fmt)
	String fmt
	if (strsearch(IgorInfo(2),"Macintosh",0))
		return " -->"			// for Windows, note leading space in string is necessary
	endif
	String out=""
	Variable i
	for (i=0;i<strlen(fmt);i+=1)
		out += "<"+fmt[i]
	endfor
	out += "ÐÐ> "
	return out
End
#endif

// The following were moved to Utility_JZT.ipf  on Dec 10, 2014
//
//	moved MenuItemIfWindowTypes(), MenuItemIfScrapValidWindowInfo(), GetWindowInfo2Scrap(), Euer2Quaternion(), 
//		PutSizeWindow(), PutScrapGraphAxis(), PutGizmoQuaternion()
//
// moved to SquareUpGizmo()
//
// moved to SetAspectToSquarePixels()
//
// moved to printWave(), printvec(), printmat(), printmatOneListReal(), printmatOneListComplex()
//
// moved to GenericGraphStyle(), Generic_Graph_Style(), GenericGraphStyleTemplate()


//  ====================================================================================  //
//  =============================== Start of LaueGo Init ===============================  //
//
//Static Function InitLaueGoFirst()
//	if (JZTalwaysFirst_Version<JZTalwaysFirst_Version_Min)
//		String str
//		sprintf str,"\"always first.ipf\" is only version %g, but we require version>=%g\rTime to update",JZTalwaysFirst_Version,JZTalwaysFirst_Version_Min
//		DoAlert 0,str
//	endif
//End
//
//  ================================ End of LaueGo Init ================================  //
//  ====================================================================================  //
