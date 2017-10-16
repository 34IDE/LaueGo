#pragma TextEncoding = "MacRoman"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName=LatticeSym
#pragma version = 6.35
#include "Utility_JZT" version>=4.44
#include "xtl_Locate"										// used to find the path to the materials files (only contains CrystalsAreHere() )

// #define 	OLD_LATTICE_ORIENTATION					// used to get old direct lattice orientation (pre version 5.00)
//	#define DO_HEXAGONAL_EXTRA							// If this is used, then Hex & Trigonal F's are too big

Static strConstant NEW_LINE="\n"						//	was NL="\r"
Static Constant minPossibleBondLength = 0.050		// 0.050 nm = 50 pm, minimum possible distance between atoms (smallest known bond is 74 pm)
Static Constant ELEMENT_Zmax = 118
//Static strConstant BAR_FONT_ALWAYS = "Arial"		//	unicode Overline only works well for Arial and Tahoma fonts, a Qt problem
strConstant BAR_FONT_ALWAYS = "Tahoma"			//	unicode Overline only works well for Arial and Tahoma fonts, a Qt problem
Static strConstant OVERLINE = "\xCC\x85"			// put this AFTER a character to put a bar over it (unicode U+0305), see:  https://en.wikipedia.org/wiki/Overline


//	remember to execute    InitLatticeSymPackage()
//
//	with version 2.0, added a text description to the crystalStructure structure
//
//	with version 3.0, fixed up and changed the symmetry operators, I now think that they are right
//	with version 3.1, only compute symmetry operations as needed
//	with version 3.2, fixed Trigonal, now can use rhombohedral axes
//	with version 3.3, fixed lowestAllowedHKL() made alpha, beta & gamma non-zero 
//	with version 3.41, fixed Fstruct(), removed multiplicity since it is taken care of with by calculating all atom positions
//	with version 3.44, added PrimitiveCellFactor()
//	with version 3.45, fixed PrimitiveCellFactor(), to deal proplerly with other systems
//	with version 3.46, made Fstruct() not forbidden for non-integral h,k,or l
//	with version 3.47, restrict Rhombohedral to only SG={146,148,155,160,161,166,167} symmetry symbols starting with an "R"
//	with version 3.48, UpdateCrystalStructureDefaults() changed to remove hashID when structure changes
//	with version 3.50, setDirectRecip() changed to handle Unconventional correctly & have at least 1 dummy atom, reMakeAtomXYZs() changed to make one atom when none are given
//	with version 3.51, changed allowed reflection limit from 0.1 electron/atom to 0.01
//	with version 3.52, added MakeSymmetryOps() to LatticeSym.ipf
//	with version 3.53, problem setting lattice parameters for Trigonal (using Hex axes) in LatticePanelParamProc()
//	with version 3.54, changed UpdateCrystalStructureDefaults() to first write structstr to Packages
//	with version 3.55, more changing how crystalStructStr is saved.
//	with version 3.56, changed again how crystalStructStr is saved.  Changed FillCrystalStructDefault() & UpdateCrystalStructureDefaults()
//	with version 3.57, added ability to read xml files.  Modified readCrystalStructure(), and added: readFileXML(), RemoveCommentsXML(), TagAttributesXML(), TagContentsXML(), ConvertUnits2meters()
//	with version 3.58, moved all interpretation of xml to Utilit_JZT.ipf, deleted:  RemoveCommentsXML(), TagAttributesXML(), TagContentsXML()
//	with version 3.59, improved ConvertUnits2meters()
//	with version 3.61, Added ability to write current xtal structure to a file.
//	with version 3.62, Added #include "Utility_JZT", previously it was not explicitly stated.
//	with version 3.63, Included ability of Fstruct() to use Cromer-Liberman values for atomic structure factors (if CromerLiberman.ipf is loaded).
//	with version 3.64, added FindMaterialsFile(), so now will automatically search for local copies of materials files.
//	with version 3.65, made the following functions ThreadSafe.:
//		dSpacing(), ZfromLabel(), Z2symbol(), num2StrFull(), SingleWhiteSpace(), MakeSymmetryOps(), isMatInMats(), getHMsym()
//		latticeSystem(), ALLOW_FC(), ALLOW_BC(), ALLOW_CC(), ALLOW_AC(), ALLOW_RHOM_HEX(), ALLOW_HEXAGONAL()
//		isRhombohedral(), PrimitiveCellFactor(), lowestOrderHKL(), hkl2str(), placesOfPrecision()
//	with version 3.67, fixed up minus2bar() and added hkl2IgorBarStr()
//	with version 3.68, fixed up Rhombohedral/Hex issue in LatticePanelParamProc(), and prints Rhombohedral constants in print_crystalStructure()
//	with version 3.69, changed getSymString() to getHMsym(), added getFullHMSym(), getHMboth(), and getHallSymbol().  DO NOT use getSymString any more.
//	with version 3.73, moved str2recip() to here from microGeometryN and other places
//	with version 3.74, changed setDirectRecip(), added a default name when no atoms are defined
//	with version 3.75, added:  Hex2Rhom(aH,cH)  &  Rhom2Hex(aR,alpha)
//	with version 3.76, changed: FillCrystalStructDefault() so to save the lattice in Pacakges whenever it has to load from preferences
//	with version 3.77, aded: findClosestHKL(), finds the hkl that is closest to the given d-spacing
//	with version 3.78, fixed SERIOUS BUG in readCrystalStructureXML()
//	with version 3.79, allowed xml input files to contain <fract_xyz>0 0 0</fract_xyz> instead of separate fract_x, fract_y, fract_z
//		added a list of bond types & Temperature, this changes crystalStructure, bondTypeStructure, print_crystalStructure(), readFileXML()
//	with version 3.80, added full support for Temperature, including proper handling of Temperature units (supports C, K, R, & F, default is C)
//	with version 3.81, changed "bond_type" to "bond_chemical", also changed bonds to use label0 & label1 (not name0, name1)
//	with version 3.82, added menus for AtomView.ipf and PowderPatterns.ipf
//		also moved menus for AtomView.ipf and PowderPatterns.ipf to the Analysis->Lattice sub menu
//	with version 3.83, modified symmtry2SG() and added SymString2SG(), this makes it easy to find Space Group number when you know the symbol
//	with version 3.90, modified the panel so that user can use wild cards to find the correct space group. LOTS of changes
//	with version 3.92, changed the hash to work properly with atom positions, changed LatticeBad(), added atomBAD(), xtalHashID(), copy_atomType(), modified copy_xtal() & UpdateCrystalStructureDefaults()
//	with version 3.93, changed input of atom positions, to use a table rather than a dialog, changed EditAtomPositions()
//	with version 3.96, added support for Debye-Waller factor, added some optional arguments to Fstruct()
//	with version 3.97, changed Fstruct() to use Q properly, and changed getFstruct(), added inputs for Temperature & keV
//	with version 3.99, added proper support for Biso, Uiso, Uij, and WyckoffSymbol
//		also changed ConvertUnits2meters() to handle units like "nm^2" or "Angstrom^-1"
//	with version 4.00, fixed editing of Biso, Uiso, Uij
//		also fixed EditAtomPositions(), atomBad(), atomThermalInfo(), Fstruct(), readFileXML() to correctly handle thermal parameters.
//	with version 4.01, changed atom editing to use a much nicer Panel.  Changed EditAtomPositions(), added SetAtomEditPanelValues(), ChangeAtomEditPanelSize(), 
//		EditAtomPanelHook(), EditAtomTypeListAction(), AtomN_PopMenuProc(), SetAtomThermalPopMenuProc(), 
//		AtomSetZ_ButtonProc(), AtomEditSetVarProc(), AtomEditWyckoffPopMenuProc(), AttomPanel2PanelStruct()
//	with version 4.03, added reading CIF files, changed readCrystalStructure(), showCrystalStructure(), FindMaterialsFile()
//		added readCrystalStructureCIF(), readFileCIF(), CIF_interpret(), CIF_loop_Labels(), CIF_readNumber(), CIF_readString()
//		also removed ForceLatticeToStructure() and UpdateCrystalStructureDefaults() from readFileXML()
//		and changed STRUCTURE_ATOMS_MAX from 20 to 40
//	with version 4.06, added ability for user to set temperature when there are thermal parameters present
//		modified crystalStructure: added haveDebyeT
//		added xtalHasDebye(), MinimalChemFormula(), xtalHasDebye(), copy_bondType()
//		changed EditAtomPositions(), print_crystalStructure(), getFstruct(), Fstruct(), copy_xtal(), LatticeBad(), ForceLatticeToStructure(), 
//		FillCrystalStructDefault(), FillLatticeParametersPanel(), UpdatePanelLatticeConstControls(), LatticePanelParamProc(), LatticePanelButtonProc()
//	with version 4.07, added valence to the atomTypeStructure, and updated other routines to properly support it 
//	with version 4.08, put valence intot the atom edit panel, removed code for editing atoms positions in a table, added valence to the Fstruct().
//	with version 4.09, added GetWyckoffSymStrings(), ForceXYZtoWyckoff(),  FindWyckoffSymbol()
//	with version 4.10, now show AtomView & PowderPattern on Lattice Panel
//	with version 4.13, changed STRUCTURE_ATOMS_MAX from 40 to 50 (50 is the maximum allowed by Igor)
//	with version 4.14, fixed up how thermal info is handled, changes to:
//		atomThermalInfo(), SetAtomEditPanelValues(), AtomPanel2PanelStruct(), print_crystalStructure(), atomThermalInfo(), 
//		crystalStructure2xml(), Fstruct(),  and added CleanOutCrystalStructure()
//	with version 4.15, more changes along with vers 4.14, changed: setLattice(), FillLatticeParametersPanel(), readFileXML(),
//		LatticePanelButtonProc(), readCrystalStructure_xtl(), readCrystalStructureXML(), CIF_interpret(), read_cri_fileOLD()
//	with version 4.16, added OverOccupyList() and OverOccupyList2Str(), and changed other routines to print a warning when occupy>1
//		changed: print_crystalStructure(), readCrystalStructure90, LatticePanelButtonProc()
//	with version 4.17, changed MinimalChemFormula(), it better removes common factors
//	with version 4.18, fixed an ERROR in Fstruct, it now properly checks for valid xtal.SpaceGroup
//		also added recipFrom_xtal(), a little utility that make a free recip matrix from xtal
//	with version 4.19, small change in reMakeAtomXYZs(), avoids error when debug on "NVAR SVAR WAVE Checking" is on
//	with version 4.20, the xtl file is now carried along in the crystalStructure structure
//	with version 4.21, added directFrom_xtal(), which is a lot like recipFrom_xtal()
//	with version 4.23, added str2hkl(hklStr,h,k,l)
//	with version 4.24, moved ConvertUnits2meters(), ConvertTemperatureUnits() , and placesOfPrecision() to Utility_JZT.ipf
//		also added NearestAllowedHKL(xtal,hkl)
//	with version 4.25, added equivalentHKLs(xtal,hkl0,[noNeg]), get list of hkl's symmetry equivalent to hkl0
//	with version 4.26, added bond info to the AtomView menu
//	with version 4.27, when writing files, use "\n" instead of "\r" for line termination. 
//		also changed crystalStructure2xml() and convertOneXTL2XML() to take a new line argument.
//	with version 4.28, in positionsOfOneAtomType() duplicate atoms are now within 50pm (in x,y,&z) not just 1e-3
//		also positionsOfOneAtomType() uses free waves rather than real temp waves
//	with version 4.30, in str2recip(str), now handles both "}{" and "},{" type strings
//		also added some helpful comments about convert recip <--> direct using MatrixOP
//	with version 4.31, added wave note info in directFrom_xtal(xtal) and in recipFrom_xtal(xtal)
//	with version 4.32, added direct2LatticeConstants(direct)
//	with version 4.33, added isValidLatticeConstants(direct), and imporved formatting in print_crystalStructure()
//	with version 4.34, added DescribeSymOps(direct), CheckDirectRecip()
//	with version 4.35, MatrixOP in Fstruct(), change variable Q --> Qmag, and calc of F uses complex math now.
//		also in positionsOfOneAtomType, make condition based on true atom distances, & use MatrixOP too.
//	with version 4.37, changed allowedHKL(), now a ThreadSafe version of allowedHKL()
//		allowedHKL() now cannot use Cromer, but it is ThreadSafe (done for indexing routines).
//		to use allowedHKL() with MultiThread you MUST pass the wave refs to the atomN waves
//		also added isValidSpaceGroup()
//	with version 4.38, added ELEMENT_Symbols and changed Get_f_proto() to be a better.
//	with version 4.39, Fixed ERROR in Fstruct()
//	with version 4.40, Get_f_proto() was ThreadSafe, it must NOT be ThreadSafe (since Get_F is not)
//	with version 4.41, removed ELEMENT_Symbols (it is now in Utility_JZT.ipf)
//	with version 4.42, fixed problem with Rhombohedal lattices & GetLatticeConstants()
//	with version 4.43, added showPanel to InitLatticeSymPackage()
//	with version 4.44, definition of MenuItemIfWindowAbsent() was changed
//	with version 4.45, added SetToDummyXTAL() for use by FillCrystalStructDefault() for startup bug found by Jan
//	with version 4.46, added SetToDummyATOM() and cleaned up v4.45
//	with version 4.48, fixed a bug in positionsOfOneAtomType() & equivalentHKLs() occured when rowRepeat(vec,1)
//
//	with version 5.00, changed definition of direct lattice to a||x and c*|| z (was b||y), added #define OLD_LATTICE_ORIENTATION to get old way
//	with version 5.04, changed readFileXML()
//	with version 5.05, changed setDirectRecip() so that Rhombohedral sytems orient a,b,c about the 111
//	with version 5.06, small menu fix
//	with version 5.07, improved FindMaterialsFile()
//	with version 5.09, FindMaterialsFile() also looks in Documents for "materials" folder
//	with version 5.10, allow all thermal parameters to be set, and Uij can be negative too
//	with version 5.12, can now also get name from _chemical_name_mineral
//	with version 5.13, can now also get name from _chemical_name_mineral working correctly
//	with version 5.16, added print_crystalStructStr()
//	with version 5.17, ARING should be "\201",  NOT "\305"
//	with version 5.18, structureTitle, now adjusts font size to show all of the string
//	with version 5.19, added DO_HEXAGONAL_EXTRA, Fstruct was too big by factor of 2
//	with version 5.20, added muOfXtal() and get_muOfXtal(), MenuItemIfCromerPresent(), and Hex2RhomFractionalFractonal()
//	with version 5.24, cleaned up Hex2Rhom(), Rhom2Hex(), Rhom2HexFractonal(), Hex2RhomFractonal(), and test_Hex_Rhom()
//	with verison 5.25, added RhomLatticeFromHex(directH), and HexLatticeFromRhom(directR)
//	with verison 5.26, moved DEGREESIGN, BULLET, ARING, BCHAR to Utility_JZT.ipf
//	with verison 5.27, remove materialsXML, change path materials to materialsPath
//	with verison 5.29, fixed bug in FindMaterialsFile() involving materialsPath
//
//	with verison 6.00, Stop using just SpaceGroup [1,230], switch to using the complete set of 530 types
//								Added SpaceGroupID and SpaceGroupIDnum to the crystalStructure structure
//								This is important since some of the space groups have many variants (e.g. SG=15 has 18 different ways of using it)
//	with verison 6.06, added keV to the Lattice Panel, also Get_f_proto(), Svector should not have the "/10"
//	with verison 6.08, changed space_group_ID --> space_group_id (to align with CIF usage)
//	with verison 6.09, in muOfXtal() was getting mult and occ wrong
//	with verison 6.11, fixed the CIF file reading, it now uses NextLineInBuf() from Utility_JZT.ipf
//	with verison 6.12, fixed a bug with valence when reading a CIF file, also tries to read multiplicity.
//	with verison 6.14, changed equivXYZM and equivXYZB to end with SpaceGroupID (not just SG)
//	with verison 6.15, changed slightly the formula for putting fractional coords into range [0,1).
//	with verison 6.16, added reading of sym ops from CIF files, can now get space group id from sym ops too.
//	with verison 6.17, fixes to the reading for sym ops from CIF files, changed ParseOneSymEquation(), modfied setSymLineIDnum() to emphasize duplicates
//	with verison 6.19, changed definition of num2fraction(), now uses tolerance
//	with verison 6.20, removed unnecessay copy_xtal(), copy_atomType(), copy_bondType()
//	with verison 6.21, changed lowestOrderHKL() for better speed, uses gcd(), also now returns the integer divisor used
//	with verison 6.24, added FstructMax to be used to find the minimum F for an allowed reflection, no longer just assumes 0.01 electrons
//	with verison 6.25, fixed error in space group numbering (affected Orthorhombic & Tetragonal)
//	with verison 6.26, fixed sign problem in lowestOrderHKL()
//	with verison 6.27, fixed up reading xml files in readCrystalStructureXML()
//	with version 6.28, simplified minus2bar(), and changed hkl2IgorBarStr()
//	with version 6.30, fixed error writing WyckoffSymbol to an xml file
//	with version 6.31, start to add ability to transform between different settings for the same Space Group.
//	with version 6.32, change all Wyckoff routines to use SpaceGroupID rather than SpaceGroup
//	with version 6.34, Wyckoff routines may be done, and can change the xtal setting

//	Rhombohedral Transformation:
//
//	for Rhombohedral (hkl), and Hexagonal (HKL)
//		h = (2H+K+L)/3
//		k = (-H+K+L)/3
//		l = (-H-2K+L)/3
//		and the condition (-H+K+L) = 3n (where n is an integer) for allowed reflections, three Rhombohedral cells per Hexagonal cell
//
//		H = h=k
//		K = k-l
//		L = h+k+l
//
//	for Hexagonal lattice constants aH,cH (alpha=beta=90, gamma=120)
//	a(Rhom) = sqrt(3*aH^2 + cH^2)/3
//	sin(alpha(Rhom)/2) = 1.5/sqrt(3+(cH/aH)^2)
//	These equations are implented in the routines:  Hex2Rhom(aH,cH)  &  Rhom2Hex(aR,alpha)
//
//
//	Although not used here, note that the following also works:
//		NOTE Inv(A^t) = Inv(A)^t, so the order of Inv() and ^t do not matter
//		MatrixOP recipLatice   = 2*PI * (Inv(directLattice))^t
//		MatrixOP directLattice = 2*PI * Inv(recipLatice^t)
//		Vc      = MatrixDet(directLattice)		// Volume of unit cell
//		VcRecip = MatrixDet(recipLatice)		// Volume of reciprocal lattice cell
//
//		kf^ = ki^ - 2*(ki^ . q^)*q^		note: (ki^ . q^) must be NEGATIVE (both Bragg & Laue)
//
//
//							SG				  idNum			#IDs		#SpaceGroups
//	Cubic				[195,230]		[489,530]		 42			36
//	Hexagonal		[168,194]		[462,488]		 27			27
//	Trigonal			[143,167]		[430,461]		 32			25
//	Tetragonal		[75,142]			[349,429]		 81			68
//	Orthorhombic	[16,74]			[108,348]		241			59
//	Monoclinic		[3,15]			[3,107]			105			13
//	Triclinic		[1,2]				[1,2]				  2			 2
//
//											SG										idNum
//	Rhombohedral	[146,148,155,160,161,166,167]	[434,437,445,451,453,459,461]


Menu "Analysis"
	SubMenu "Lattice"
		"<BLattice Set Panel...",MakeLatticeParametersPanel("")
		help={"Define the crystal structure and lattice to be used"}
		"Show Crystal Structure",showCrystalStructure()
		help={"Shows the crystal structure and lattice that is currently defined"}
		"  Edit the Atom Positions...",EditAtomPositionsMenu()
		"  Change Current xtal Setting...", ChangeSettingCurrentXtal("")
		help={"Manually set/change the atom positions"}
		"Load a new Crystal Structure...",LoadCrystal("")
		help={"load a rystal structure from a fie"}
		"Write Current Crystal to XML File...",writeCrystalStructure2xmlFile("","")
		help={"takes the current crystal structure and writes it to an xml file"}
		"d[hkl]",/Q,get_dhkl(NaN,NaN,NaN,T=NaN)
		help={"show d-spacing for the hkl reflection"}
		"Fstructure [hkl]", /Q ,getFstruct(NaN,NaN,NaN)
		help={"Crude Structure Factor using current lattice"}
		LatticeSym#MenuItemIfCromerPresent("Set Energy..."), setCromerEnergy(NaN)
		help={"The energy that is used by the Cromer-Liberman"}
		LatticeSym#MenuItemIfCromerPresent("calculate mu..."), get_muOfXtal(NaN)
		help={"Calculate mu (absorption factor) using current lattice structure"}
		"Find Closest hkl from d-spacing or Q",findClosestHKL(NaN)
		help={"Knowing either the d-spacing or the Q, find closest hkl's"}
		"\\M0Space Group number <--> symmetry",symmtry2SG("")
		help={"find the Space Group number from symmetry string,  e.g. Pmma, or sym from number"}
		"Describe the Symmetry Operations", DescribeSymOps($"")
		"angle between two hkl's",angleBetweenHKLs(NaN,NaN,NaN,  NaN,NaN,NaN)
		"  Convert old xtl files to new xml files",ConverXTLfile2XMLfile("")
		"-"
		MenuItemIfWindowAbsent("Include Powder Patterns Support","PowderPatterns.ipf","WIN:128"), Execute/P "INSERTINCLUDE  \"PowderPatterns\", version>=0.24";Execute/P "COMPILEPROCEDURES ";Execute/P "Init_PowderPatternLattice()"
		help = {"Load Function used to compute Powder Patterns from Loaded Lattice"}
		MenuItemIfWindowAbsent("Include Atom View Support","AtomView.ipf","WIN:128"), Execute/P "INSERTINCLUDE  \"AtomView\", version>=0.43";Execute/P "COMPILEPROCEDURES ";Execute/P "Init_AtomViewLattice()"
		help = {"Load Function used to Display Atoms from the Loaded Lattice in a Gizmo"}
	End
End

Static Constant STRUCTURE_ATOMS_MAX=50	// max number of atom types in a material structure (for Si only need 1)
//Strconstant CommonPredefinedStructures = "FCC:255;BCC:229;Diamond:227;Perovskite:221;Cubic:195;Hexagonal:194;Wurtzite (B4):186;Sapphire:167;"
Strconstant CommonPredefinedStructures = "FCC:255;BCC:229;Diamond:227;Perovskite:221;Cubic:195;Wurtzite (B4):186;Sapphire:167;"
Static Constant HASHID_LEN = 66				// length of string to hold hashID in crystalStructure
Static Constant MAX_FILE_LEN = 400			// max length of a file name to store

Static Constant amu_eV = 931.494061e6		// energy of one amu (eV),  these 4 numbers only used by Debye-Waller calculation
Static Constant hbar = 6.58211928e-16		// reduced Plank constant (eV - sec)
Static Constant kB = 8.6173324e-5			// Boltzman constant (eV / K)
Static Constant c = 299792458				// speed of light (m/sec)

Static Constant hc_keVnm = 1.2398419739		// h*c (keV-nm),  these two used to calculate mu in muOfXtal
Static Constant re_nm = 2.8179403227e-06	// Thompson radius (nm)

// These constants have been moved to Utility_JZT.ipf
//	#if (IgorVersion()<7)
//		strConstant DEGREESIGN = "\241"			// option-shift-8
//		strConstant BULLET = "\245"				// option-8
//		strConstant ARING = "\201"				// Angstrom sign, option-shift-A
//	#if StringMatch(IgorInfo(2),"Windows")
//		strConstant BCHAR = "\257"
//	#else
//		strConstant BCHAR = "\321"
//	#endif
//	#else
//		strConstant DEGREESIGN = "\xC2\xB0"	// UTF8, DEGREE SIGN
//		strConstant BULLET = "\xE2\x80\xA2"
//		strConstant ARING = "\xC3\x85"			// Aring, Angstrom sign
//		strConstant BCHAR = "\xE2\x80\x94"		// EM DASH
//	#endif

Static Function/S MenuItemIfCromerPresent(item)		// Shows menu item if CromerLiberman is present
	String item			// the string that appears in the menu
	Variable present = strlen(WinList("CromerLiberman.ipf","","WIN:129"))>0
	return SelectString(present,"(","")+item
End


// =========================================================================
// =========================================================================
//	Start of Structure definitions

Structure crystalStructure	// structure definition for a crystal lattice
	char desc[100]					// name or decription of this crystal
	double a,b,c					// lattice constants, length (nm)
	double alpha,beta,gam		// angles (degree)
	int16 SpaceGroup				// Space Group number from international tables, allowed range is [1, 230]
	char SpaceGroupID[12]		// id of SpaceGroup, e.g. "15:-b2", not just a number anymore
	int16 SpaceGroupIDnum		// index to the SpaceGroupID, allowed range is [1, 530]
	double  a0,  b0,  c0		// direct lattice from constants { a[], b[], c[] }
	double  a1,  b1,  c1
	double  a2,  b2,  c2
	double  as0,  bs0,  cs0	// reciprocal lattice { a*[], b*[], c*[] }
	double  as1,  bs1,  cs1	// a*,b*,c* already have the 2PI in them
	double  as2,  bs2,  cs2
	double Vc						// volume of cell, triple product of (a[]xb[]).c
	double density					// calculated density (g/cm^3)
	double Temperature			// Temperature (C)
	double	alphaT				// coef of thermal expansion, a = ao*(1+alphaT*(TempC-22.5))
	int16 N							// number of atoms described here
	Struct atomTypeStructure atom[STRUCTURE_ATOMS_MAX]
	int16 Vibrate					// True if DebyeT, Biso, Uiso, or Uij available for some atom
	int16 haveDebyeT				// True if one of the atoms has a Debye Temperature (a Temperature dependent thermal parameter)
	int16 Nbonds					// number of bonds described here
	Struct bondTypeStructure bond[2*STRUCTURE_ATOMS_MAX]
	double Unconventional00,Unconventional01,Unconventional02	// transform matrix for an unconventional unit cel
	double Unconventional10,Unconventional11,Unconventional12
	double Unconventional20,Unconventional21,Unconventional22
	char sourceFile[MAX_FILE_LEN]	// name of source file
	char hashID[HASHID_LEN]	// hash function for this strucutre (needs to hold at least 64 chars), This MUST be the LAST item
EndStructure
//
Structure crystalStructure5	// structure definition for a crystal lattice
	char desc[100]					// name or decription of this crystal
	double a,b,c					// lattice constants, length (nm)
	double alpha,beta,gam		// angles (degree)
	int16 SpaceGroup				// Space Group number from international tables, allowed range is [1, 230]
	double  a0,  b0,  c0		// direct lattice from constants { a[], b[], c[] }
	double  a1,  b1,  c1
	double  a2,  b2,  c2
	double  as0,  bs0,  cs0	// reciprocal lattice { a*[], b*[], c*[] }
	double  as1,  bs1,  cs1	// a*,b*,c* already have the 2PI in them
	double  as2,  bs2,  cs2
	double Vc						// volume of cell, triple product of (a[]xb[]).c
	double density					// calculated density (g/cm^3)
	double Temperature			// Temperature (C)
	double	alphaT				// coef of thermal expansion, a = ao*(1+alphaT*(TempC-22.5))
	int16 N							// number of atoms described here
	Struct atomTypeStructure atom[STRUCTURE_ATOMS_MAX]
	int16 Vibrate					// True if DebyeT, Biso, Uiso, or Uij available for some atom
	int16 haveDebyeT				// True if one of the atoms has a Debye Temperature (a Temperature dependent thermal parameter)
	int16 Nbonds					// number of bonds described here
	Struct bondTypeStructure bond[2*STRUCTURE_ATOMS_MAX]
	double Unconventional00,Unconventional01,Unconventional02	// transform matrix for an unconventional unit cel
	double Unconventional10,Unconventional11,Unconventional12
	double Unconventional20,Unconventional21,Unconventional22
	char sourceFile[MAX_FILE_LEN]	// name of source file
	char hashID[HASHID_LEN]	// hash function for this strucutre (needs to hold at least 64 chars), This MUST be the LAST item
EndStructure
//
Structure atomTypeStructure	// defines one type of atom in a crystal structure
	char name[60]				// label for this atom, usually starts with atomic symbol
	int16 Zatom					// Z of the atom
	double x						// fractional coord along x lattice vector
	double y
	double z
	double occ					// occupancy
	char WyckoffSymbol[2]	// a single letter, e.g. 'a', also called the Wyckoff letter
	int16 valence				// valence of atom, must be an integer
	int16 mult					// multiplicity
	double DebyeT				// Debye Temperature (K),  for DebyeT, B, Uiso, & U_ij, use only one method
	double Biso					// B-factor (nm^2) using:   exp(-M) = exp(-B * sin^2(theta)/lam^2)		B = 8 * PI^2 * <u^2> =  8 * PI^2 * Uiso, 	exp[-B*q^2 / (16 PI^2) ]
	double Uiso					// isotropic U (nm^2)
	double U11					// anisotropic U(11) (nm^2)
	double U22					// anisotropic U(22) (nm^2)
	double U33					// anisotropic U(33) (nm^2)
	double U12					// anisotropic U(12) (nm^2)
	double U13					// anisotropic U(13) (nm^2)
	double U23					// anisotropic U(23) (nm^2)
EndStructure
//
Structure bondTypeStructure	// defines the type of bond between two atom types
	char label0[60]			// label for first atom, usually starts with atomic symbol
	char label1[60]			// label for second atom, usually starts with atomic symbol
	int16 N						// number of bonds in len (often just 1)
	double len[5]				// length of bond (possibly multiple values) (nm)
EndStructure


Function init_crystalStructure(xtal)		// set all values to empty or invalid values
	STRUCT crystalStructure &xtal
	
	xtal.desc = ""
	xtal.a = NaN		;		xtal.b = NaN		;	xtal.c = NaN
	xtal.alpha = NaN	;		xtal.beta = NaN	;	xtal.gam = NaN
	xtal.SpaceGroup = 0
	xtal.SpaceGroupID = ""	;	xtal.SpaceGroupIDnum = 0
	xtal.a0 = NaN		;	xtal.a1 = NaN	;		xtal.a2 = NaN
	xtal.b0 = NaN		;	xtal.b1 = NaN	;		xtal.b2 = NaN
	xtal.c0 = NaN		;	xtal.c1 = NaN	;		xtal.c2 = NaN

	xtal.Vc = NaN
	xtal.density = NaN
	xtal.Temperature = 0	;	xtal.alphaT = 0
	xtal.Vibrate = 0			;	xtal.haveDebyeT = 0

	xtal.N = 0
	Variable i
	for (i=0; i<STRUCTURE_ATOMS_MAX; i+=1)
		init_atomTypeStructure(xtal.atom[i])
	endfor

	xtal.Nbonds = 0
	for (i=0; i<(2*STRUCTURE_ATOMS_MAX); i+=1)
		init_bondTypeStructure(xtal.bond[i])
	endfor

	xtal.Unconventional00 = NaN	;	xtal.Unconventional01 = NaN	;	xtal.Unconventional02 = NaN
	xtal.Unconventional10 = NaN	;	xtal.Unconventional11 = NaN	;	xtal.Unconventional12 = NaN
	xtal.Unconventional20 = NaN	;	xtal.Unconventional21 = NaN	;	xtal.Unconventional22 = NaN
	xtal.sourceFile = ""			;	xtal.hashID = ""
End

Function init_atomTypeStructure(atom)		// set all values to empty or invalid values
	STRUCT atomTypeStructure &atom
	atom.name = ""
	atom.Zatom = 0
	atom.x = NaN	;	atom.y = NaN	;	atom.z = NaN
	atom.occ = 1
	atom.WyckoffSymbol = ""
	atom.valence = 0
	atom.mult = 1
	atom.DebyeT = NaN	;	atom.Biso = NaN	;	atom.Uiso = NaN
	atom.U11 = NaN		;	atom.U22 = NaN		;	atom.U33 = NaN
	atom.U12 = NaN		;	atom.U13 = NaN		;	atom.U23 = NaN
End

Function init_bondTypeStructure(bond)		// set all values to empty or invalid values
	STRUCT bondTypeStructure &bond
	bond.label0 = ""
	bond.label1 = ""
	bond.N = 0
	bond.len[0] = NaN		;	bond.len[1] = NaN		;	bond.len[2] = NaN
	bond.len[3] = NaN		;	bond.len[4] = NaN
End

//	End of Structure definitions
// =========================================================================
// =========================================================================


// =========================================================================
// =========================================================================
//	Start of setting particular lattice constants

Function MakeLatticeParametersPanel(strStruct)
	String strStruct									// optional passed value of geo structure, this is used if passed
	if (strlen(WinList("LatticeSet",";","WIN:64")))	// window alreay exits, bring it to front
		DoWindow/F LatticeSet
		return 0
	endif
	Variable showEnergy = (exists("Get_f")==6)
	Variable bottom = 60+508 + (showEnergy ? 20 : 0)
	NewPanel /K=1 /W=(675,60,675+220,bottom)
	DoWindow/C LatticeSet
	FillLatticeParametersPanel(strStruct,"LatticeSet",0,0)
End


Function showCrystalStructure()						// prints the structure that is currently being used
	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))					//fill the lattice structure with test values
		DoAlert 0, "no crystal structure found"
		return 1
	endif
	String str, sym = getHMboth(xtal.SpaceGroupIDnum)
	sprintf str, "'%s'  %d atoms,  %s   %s\r%.9g, %.9g, %.9gnm,\r%g%s, %g%s, %g%s",xtal.desc,xtal.N,xtal.SpaceGroupID,sym,xtal.a,xtal.b,xtal.c,xtal.alpha,DEGREESIGN,xtal.beta,DEGREESIGN,xtal.gam,DEGREESIGN
	Variable netCharge = NetChargeCell(xtal)
	if (netCharge)
		str += "\r  *** Charge imbalance in cell = "+num2str(netCharge)+" ***"
	endif
	DoAlert 1,str+"\r\tprint all details to history too?"
	if (V_flag==1)
		print_crystalStructure(xtal)					// prints out the value in a crystalStructure structure
	else
		printf "currently using  '%s'  lattice is  %s   %s     %.9gnm, %.9gnm, %.9gnm,   %g%s, %g%s, %g%s",xtal.desc,xtal.SpaceGroupID,sym,xtal.a,xtal.b,xtal.c,xtal.alpha,DEGREESIGN,xtal.beta,DEGREESIGN,xtal.gam,DEGREESIGN
		if (xtal.N > 0)
			printf ",   %g defined atom types\r",xtal.N
		else
			print ",   NO actual atoms defined"
		endif
		if (xtal.alphaT>0)
			printf "coefficient of thermal expansion is %g\r",xtal.alphaT
		endif
	endif
	return 0
End


// ============================= Start of Atom Editing =============================
Function EditAtomPositionsMenu()		// returns number of atoms defined
	if (!DataFolderExists("root:Packages:Lattices"))
		InitLatticeSymPackage()
	endif
	if (!DataFolderExists("root:Packages:Lattices"))
		DoAlert 0, "ERROR -- Package not initialized!!!"
		return NaN
	endif

	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))				//fill the lattice structure with test values
		DoAlert 0, "ERROR -- no crystal structure found"
		return NaN
	endif

	Variable Natoms=EditAtomPositions(xtal)	// takes and MODIFYS xtal, returns number of atoms defined
	if (Natoms>=0)
		if (ForceLatticeToStructure(xtal))	// sets everything including remaking atom0,atom1,...
			return NaN
		endif
		// xtal.hashID = xtalHashID(xtal)	// this gets done in UpdateCrystalStructureDefaults()
		// reMakeAtomXYZs(xtal)				// this gets done in ForceLatticeToStructure()
		UpdateCrystalStructureDefaults(xtal)
	endif
	return Natoms
End
//

Function EditAtomPositions(xtal_IN)		// Create and Handle the Edit Atoms Panel
	STRUCT crystalStructure &xtal_IN
	if (LatticeBad(xtal_IN))
		DoAlert 0, "ERROR -- Invalid Lattice in EditAtomPositions()"
		return NaN
	endif
	STRUCT crystalStructure xtal					// the working copy
	xtal = xtal_IN										// copy xtal_IN to a working copy, was:  copy_xtal(xtal,xtal_IN)

	SetSymOpsForSpaceGroup(xtal.SpaceGroupID)	// ensure that symmetry ops are right
	Variable Natoms = round(xtal.N)				// number of predefined atoms
	Natoms = numtype(Natoms) ? 0 : limit(Natoms,0,STRUCTURE_ATOMS_MAX)
	String wyckMenuStr = "\" ;"+WyckoffMenuStr(xtal.SpaceGroupID)+"\""

	if (Natoms<1)
		DoAlert 1, "No atoms defined, create one?"
		if (V_flag==1)									// create one dummy atom
			Natoms = 1
			xtal.N = 1
			xtal.atom[0].name = ""
			xtal.atom[0].Zatom = 1
			xtal.atom[0].x=0 ;			xtal.atom[0].y=0 ;			xtal.atom[0].z=0
			xtal.atom[0].valence = 0
			xtal.atom[0].occ = 1
			xtal.atom[0].WyckoffSymbol = ""
			xtal.atom[0].DebyeT=NaN ;	xtal.atom[0].Biso=NaN ;	xtal.atom[0].Uiso=NaN
			xtal.atom[0].U11=NaN ;		xtal.atom[0].U22=NaN ;		xtal.atom[0].U33=NaN
			xtal.atom[0].U12=NaN ;		xtal.atom[0].U13=NaN ;		xtal.atom[0].U23=NaN
		else
			return NaN									// don't change anything
		endif
	endif

	String listName=UniqueName("atomNameList",1,0)
	Make/O/T/N=(Natoms) $listName/WAVE=atomNameList
	listName = GetWavesDataFolder(atomNameList,2)
	atomNameList = xtal.atom[p].name

	NewPanel /W=(74,44,403,301-24)/K=1 as "Edit Atoms,  Close When Done..."
	SetDrawLayer UserBack
	SetDrawEnv linethick= 2
	DrawLine 0,164,330,164
	ListBox pickAtom,pos={10,5},size={174,100},proc=LatticeSym#EditAtomTypeListAction,fSize=12
	ListBox pickAtom,listWave=$listName,mode= 2,selRow=0
	PopupMenu NatomPopup,pos={212,1},size={89,20},proc=LatticeSym#AtomN_PopMenuProc,title="N atoms",fSize=12
	PopupMenu NatomPopup,popvalue="0",value= #"\"0;1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;\""
	PopupMenu NatomPopup,mode=Natoms+1
	SetVariable setAtomLabel,pos={202,23},size={107,19},bodyWidth=76,title="label",proc=LatticeSym#AtomEditSetVarProc
	SetVariable setAtomLabel,fSize=12,value= _STR:""
	SetVariable setAtomZ,pos={203,47},size={60,19},title="Z",fSize=12,proc=LatticeSym#AtomEditSetVarProc
	SetVariable setAtomZ,limits={1,92,1},value= _NUM:0
	PopupMenu setValencePopup,pos={216,70},size={90,19},bodyWidth=40,title="Valence"
	PopupMenu setValencePopup,value= #"\"-2;-1;0;1;2;3;4;5;6;\""
	PopupMenu setValencePopup,fSize=12,mode=3,proc=LatticeSym#AtomEdit_PopMenuProc	// mode=3 sets valence=0
	PopupMenu setAtomWyckoff,mode=1,popvalue=" ",value= #wyckMenuStr
	PopupMenu setAtomWyckoff,pos={213,93},size={87,20},title="Wyckoff",fSize=12,proc=LatticeSym#AtomEdit_PopMenuProc
	PopupMenu useThermalPopUp,pos={211,141},size={113,20},fSize=12, proc=LatticeSym#SetAtomThermalPopMenuProc
	PopupMenu useThermalPopUp,mode=5,popvalue="anisotropic Uij",value= #"\"NO Thermal;Debye T;B-isotropic;U-isotropic;anisotropic Uij\""
	Button setAtomZButton,pos={265,46},size={40,20},proc=LatticeSym#AtomSetZ_ButtonProc,title="set Z"
	Button setAtomZButton,fSize=10
	SetVariable setAtomRelX,pos={14,118},size={90,19},bodyWidth=80,title="x/a",proc=LatticeSym#AtomEditSetVarProc
	SetVariable setAtomRelX,fSize=12,limits={0,1,0},value= _NUM:0
	SetVariable setAtomRelY,pos={123,118},size={91,19},bodyWidth=80,title="y/b",proc=LatticeSym#AtomEditSetVarProc
	SetVariable setAtomRelY,fSize=12,limits={0,1,0},value= _NUM:0
	SetVariable setAtomRelZ,pos={234,118},size={91,19},bodyWidth=80,title="z/c",proc=LatticeSym#AtomEditSetVarProc
	SetVariable setAtomRelZ,fSize=12,limits={0,1,0},value= _NUM:0
	SetVariable setAtomOcc,pos={10,141},size={110,19},title="Occupancy",fSize=12,proc=LatticeSym#AtomEditSetVarProc
	SetVariable setAtomOcc,limits={0,1,0},value= _NUM:1

	SetVariable setAtomDebyeT,pos={17,173},size={127,19},bodyWidth=60,title="Debye ("+DEGREESIGN+"K)"
	SetVariable setAtomDebyeT,fSize=12,limits={0.01,inf,0},value= _NUM:NaN , proc=LatticeSym#AtomEditSetVarProc
	SetVariable setAtomBiso,pos={17,170},size={155,24},bodyWidth=60,title="B-isotropic ("+ARING+"\\S2\\M)"
	SetVariable setAtomBiso,fSize=12,limits={0.0001,100,0},value= _NUM:0 , proc=LatticeSym#AtomEditSetVarProc
	SetVariable setAtomUiso,pos={17,170},size={155,24},bodyWidth=60,title="U-isotropic ("+ARING+"\\S2\\M)"
	SetVariable setAtomUiso,fSize=12,limits={0.0001,100,0},value= _NUM:0 , proc=LatticeSym#AtomEditSetVarProc
	SetVariable setAtomU11,pos={10,173},size={73,21},bodyWidth=50,title="U\\B11\\M"
	SetVariable setAtomU11,fSize=12,limits={0.0001,100,0},value= _NUM:0 , proc=LatticeSym#AtomEditSetVarProc
	SetVariable setAtomU22,pos={113,173},size={73,21},bodyWidth=50,title="U\\B22\\M"
	SetVariable setAtomU22,fSize=12,limits={0.0001,100,0},value= _NUM:0 , proc=LatticeSym#AtomEditSetVarProc
	SetVariable setAtomU33,pos={211,173},size={73,21},bodyWidth=50,title="U\\B33\\M"
	SetVariable setAtomU33,fSize=12,limits={0.0001,100,0},value= _NUM:0 , proc=LatticeSym#AtomEditSetVarProc
	SetVariable setAtomU12,pos={11,204},size={73,21},bodyWidth=50,title="U\\B12\\M"
	SetVariable setAtomU12,fSize=12,limits={0.0001,100,0},value= _NUM:0 , proc=LatticeSym#AtomEditSetVarProc
	SetVariable setAtomU13,pos={114,204},size={73,21},bodyWidth=50,title="U\\B13\\M"
	SetVariable setAtomU13,fSize=12,limits={0.0001,100,0},value= _NUM:0 , proc=LatticeSym#AtomEditSetVarProc
	SetVariable setAtomU23,pos={212,204},size={73,21},bodyWidth=50,title="U\\B23\\M"
	SetVariable setAtomU23,fSize=12,limits={0.0001,100,0},value= _NUM:0 , proc=LatticeSym#AtomEditSetVarProc
	TitleBox titleUij1,pos={290,171},size={14,21},title=ARING+"\\S2\\M",fSize=12,frame=0
	TitleBox titleUij2,pos={291,203},size={14,21},title=ARING+"\\S2\\M",fSize=12,frame=0

	String win=StringFromList(0,WinList("*",";","WIN:64"))
	if (Natoms>0)
		SetAtomEditPanelValues(win,xtal.atom[0])
	endif

	String xtalStrStructName=UniqueName("xtalS",4,0)
	String/G $xtalStrStructName							// create temporary global to hold xtal
	SVAR strStruct=$xtalStrStructName
	StructPut/S/B=2 xtal, strStruct						// keep a copy in the Panel
	SetWindow kwTopWin,hook(TableEnter)=LatticeSym#EditAtomPanelHook,userdata(xtalSname)=xtalStrStructName
	SetWindow kwTopWin,userdata(listWave)=listName
	DoUpdate
	PauseForUser $win											// WAIT here for Panel to close then proceed
	KillWaves/Z atomNameList

	// user is done entering atom type information, save it
	StructGet_xtal(strStruct,xtal)						// retrieve xtal from the global string
	KillStrings/Z strStruct								// done with this global string
	if (LatticeBad(xtal,atomsToo=1))
		DoAlert 0, "ERROR -- Invalid Lattice"
		return NaN
	endif

	DoAlert 1,"Entered Values are Valid,\rRepalce Curret Values with These?"
	if (V_flag!=1)
		return NaN
	endif
	xtal.Vibrate = xtalVibrates(xtal)					// True if some Thermal vibration info present in xtal
	xtal.haveDebyeT = xtalHasDebye(xtal)				// True if some one of the atoms has a Debye Temperature
	xtal.hashID = xtalHashID(xtal)
	xtal_IN = xtal												// copy working copy to the passed copy, was:  copy_xtal(xtal_IN,xtal)
	return Natoms
End
//
Static Function SetAtomEditPanelValues(win,atom)			// set panel values using numbers from atom
	String win									// window to set
	STRUCT atomTypeStructure &atom

	SetVariable setAtomLabel,win=$win,value=_STR:atom.name
	SetVariable setAtomZ,win=$win,value= _NUM:atom.Zatom
	SetVariable setAtomRelX,win=$win,value= _NUM:atom.x
	SetVariable setAtomRelY,win=$win,value= _NUM:atom.y
	SetVariable setAtomRelZ,win=$win,value= _NUM:atom.z
	SetVariable setAtomOcc,win=$win,value= _NUM:atom.occ
	Variable mm=char2num(LowerStr(atom.WyckoffSymbol))-96
	mm = limit(mm,0,Inf)+1
	PopupMenu setAtomWyckoff,win=$win,mode=mm
	PopupMenu setValencePopup,win=$win,mode=(3 + max(atom.valence,-2))
	SetVariable setAtomDebyeT,win=$win,value=_NUM:atom.DebyeT
	SetVariable setAtomBiso,win=$win,value=_NUM:(atom.Biso)*100
	SetVariable setAtomUiso,win=$win,value=_NUM:(atom.Uiso)*100
	SetVariable setAtomU11,win=$win,value=_NUM:(atom.U11)*100
	SetVariable setAtomU22,win=$win,value=_NUM:(atom.U22)*100
	SetVariable setAtomU33,win=$win,value=_NUM:(atom.U33)*100
	SetVariable setAtomU12,win=$win,value=_NUM:(atom.U12)*100
	SetVariable setAtomU13,win=$win,value=_NUM:(atom.U13)*100
	SetVariable setAtomU23,win=$win,value=_NUM:(atom.U23)*100

	Variable mode = atomThermalInfo(atom)
	mode = limit(mode,0,4)								// when mode==5, use 4
	ChangeAtomEditPanelSize(win,mode)
End
//
Static Function ChangeAtomEditPanelSize(win,mode)		// changes size of AtomEditPanel to match Thermal Info
	String win									// window to set
	Variable mode								// Thermal Mode, 0=None, 1=Debye, 2=Biso, 3=Uiso, 4=Uaniso[ij]

	PopupMenu useThermalPopUp,win=$win,mode=(mode+1)		// set the Thermal popup to corrct value
	SetVariable setAtomDebyeT,win=$win, disable=1			// first disable all
	SetVariable setAtomBiso,win=$win, disable=1
	SetVariable setAtomUiso,win=$win, disable=1
	SetVariable setAtomU11,win=$win, disable=1
	SetVariable setAtomU22,win=$win, disable=1
	SetVariable setAtomU33,win=$win, disable=1
	SetVariable setAtomU12,win=$win, disable=1
	SetVariable setAtomU13,win=$win, disable=1
	SetVariable setAtomU23,win=$win, disable=1
	TitleBox titleUij1, disable=1
	TitleBox titleUij2, disable=1

	GetWindow $win wsizeRM
	V_bottom = V_top+257
	switch(mode)				// adjust bottom of window and which SetVariables are visible
		case 0:					// No Thermal
			V_bottom -= 95
			break
		case 1:					// Debye Temperature
			SetVariable setAtomDebyeT,win=$win, disable=0
			V_bottom -= 60
			break
		case 2:					// B-iso
			SetVariable setAtomBiso,win=$win, disable=0
			V_bottom -= 60
			break
		case 3:					// U-iso
			SetVariable setAtomUiso,win=$win, disable=0
			V_bottom -= 60
			break
		case 4:					// U-aniso
		default:
			V_bottom -= 24
			SetVariable setAtomU11,win=$win, disable=0
			SetVariable setAtomU22,win=$win, disable=0
			SetVariable setAtomU33,win=$win, disable=0
			SetVariable setAtomU12,win=$win, disable=0
			SetVariable setAtomU13,win=$win, disable=0
			SetVariable setAtomU23,win=$win, disable=0
			TitleBox titleUij1, disable=0
			TitleBox titleUij2, disable=0
	endswitch
	MoveWindow/W=$win V_left,V_top,V_right,V_bottom
End
//
Static Function EditAtomPanelHook(s)				// Just used to Kill the window
	STRUCT WMWinHookStruct &s
	if ((s.keycode == 119 && s.eventMod == 8) || s.eventCode == 17)	// user entered command-W, or Closed Window
		KillWindow $(s.winName)
	endif
	return 0
End
//
Static Function EditAtomTypeListAction(LB_Struct) : ListboxControl	// proc for the listbox
	// When the user clicks or up/down arrows on the list box, change values in the panel
	STRUCT WMListboxAction &LB_Struct
	if (LB_Struct.eventCode != 4)			// only process when event is cell selection
		return 0
	endif

	SVAR strStruct = $GetUserData(LB_Struct.win,"","xtalSname")	// global string containing xtal
	if (!SVAR_Exists(strStruct))
		print "ERROR -- Cannot find strStruct from panel userdata"
		return 0
	endif
	STRUCT crystalStructure xtal
	StructGet_xtal(strStruct,xtal)											// retrieve xtal from the string strStruct
	SetAtomEditPanelValues(LB_Struct.win,xtal.atom[LB_Struct.row])	// put values from xtal into panel values
	return 0
End
//
Static Function AtomN_PopMenuProc(pa) : PopupMenuControl		// Called when Number of atoms is to change
	STRUCT WMPopupAction &pa
	if (pa.eventCode != 2)		// only act on mouse up
		return 0
	endif

	Variable Natoms = pa.popNum - 1
	if (0 <= Natoms && Natoms<=STRUCTURE_ATOMS_MAX)
		SVAR strStruct = $GetUserData(pa.win,"","xtalSname")
		if (!SVAR_Exists(strStruct))
			print "ERROR -- Cannot find strStruct from panel userdata"
			return 0
		endif
		STRUCT crystalStructure xtal
		StructGet_xtal(strStruct,xtal)							// retrieve xtal from the string strStruct
		xtal.N = Natoms
		StructPut/S/B=2 xtal, strStruct							// update a copy with new xtal.N
		Wave/T atomNameList = $GetUserData(pa.win,"","listWave")	// contains names that are in list box
		Redimension/N=(Natoms) atomNameList
		if (WaveExists(atomNameList))							// update name in atomNameList from label
			atomNameList = xtal.atom[p].name					// re-set all the names
			atomNameList = SelectString(strlen(atomNameList[p]),"no label",atomNameList[p])
		endif
	else
		printf "ERROR -- Unable to set Natoms to %g,  must be in range [%d, %d]\r",Natoms,0,STRUCTURE_ATOMS_MAX
	endif
	return 0
End
//
Static Function SetAtomThermalPopMenuProc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa
	if (pa.eventCode != 2)		// only act on mouse up
		return 0
	endif
	ChangeAtomEditPanelSize(pa.win,pa.popNum - 1)
	return 0
End
//
Static Function AtomSetZ_ButtonProc(ba) : ButtonControl			// "setZ" button, try to set Z from label
	STRUCT WMButtonAction &ba
	if (ba.eventCode != 2)		// only act on mouse up
		return 0
	endif

	String str, win = ba.win
	ControlInfo/W=$win setAtomLabel				// label from control
	if (strlen(S_Value)<1)
		return 0
	endif
	str = SelectString(isdigit(S_Value[1]),S_Value[0,1],S_Value[0])
	Variable Zatom = ZfromLabel(str)
	if (Zatom == limit(round(Zatom),1,92))
		SetVariable setAtomZ,win=$win,value= _NUM:Zatom
		AtomPanel2PanelStruct(win,ba.ctrlName)	// collect all panel values and put them into working copy of xtal
	endif
	return 0
End
//
Static Function AtomEditSetVarProc(sva) : SetVariableControl	// a generic proc used by all the SetVariables in the Edit Atom Panel
	STRUCT WMSetVariableAction &sva
	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			AtomPanel2PanelStruct(sva.win,sva.ctrlName)	// collect all panel values and put them into working copy of xtal
			break
	endswitch
	return 0
End
//
Static Function AtomEdit_PopMenuProc(pa) : PopupMenuControl	// process Wyckoff popup menu
	STRUCT WMPopupAction &pa
	if ( pa.eventCode == 2 )
		AtomPanel2PanelStruct(pa.win,pa.ctrlName)	// collect all panel values and put them into working copy of xtal
	endif
	return 0
End
//Static Function AtomEditWyckoffPopMenuProc(pa) : PopupMenuControl	// process Wyckoff popup menu
//	STRUCT WMPopupAction &pa
//	if ( pa.eventCode == 2 )
//		AtomPanel2PanelStruct(pa.win)	// collect all panel values and put them into working copy of xtal
//	endif
//	return 0
//End
//
Static Function AtomPanel2PanelStruct(win,ctrlName)	// gather panel values and put them into copy of xtal stored in panel
	String win
	String ctrlName

	SVAR strStruct = $GetUserData(win,"","xtalSname")
	if (!SVAR_Exists(strStruct))
		print "ERROR -- Cannot find strStruct from panel userdata"
		return 0
	elseif (strlen(strStruct)<1)
		print "ERROR -- Could not get xtal values from Panel"
		return 0
	endif

	STRUCT crystalStructure xtal
	StructGet_xtal(strStruct,xtal)												// retrieve xtal from the string strStruct

	String symbol
	Variable Natoms, m, xx,yy,zz, mult=0
	ControlInfo/W=$win NatomPopup		;		Natoms = V_Value-1		// Number of atoms
	ControlInfo/W=$win pickAtom			;		m = V_Value					// current atom number
	if (m>=Natoms || Natoms>STRUCTURE_ATOMS_MAX)
		print "Too Many Atoms"
		return 0
	endif
	ControlInfo/W=$win setAtomLabel		;		xtal.atom[m].name = S_Value[0,59]	// gather all the values
	ControlInfo/W=$win setAtomZ			;		xtal.atom[m].Zatom = V_Value
	ControlInfo/W=$win setAtomWyckoff	;		symbol = SelectString(StringMatch(S_Value," "),S_Value,"")
	ControlInfo/W=$win setAtomRelX		;		xx = V_Value
	ControlInfo/W=$win setAtomRelY		;		yy = V_Value
	ControlInfo/W=$win setAtomRelZ		;		zz = V_Value
	if (StringMatch(ctrlName,"setAtomWyckoff") && strlen(symbol)>0)	// changed Wyckoff, force xyz to conform
		if (!ForceXYZtoWyckoff(xtal.SpaceGroupID,symbol,xx,yy,zz))
			SetVariable setAtomRelX,win=$win,value= _NUM:xx
			SetVariable setAtomRelY,win=$win,value= _NUM:yy
			SetVariable setAtomRelZ,win=$win,value= _NUM:zz
		endif
	elseif (StringMatch(ctrlName,"setAtomRel*") && strlen(symbol)==0)	// changed xyz, find Wyckoff symbol
		Variable mm
		symbol = FindWyckoffSymbol(xtal.SpaceGroupID,xx,yy,zz,mult)
		if (char2num(symbol)==65)
			mm = 27+1
		else
			mm=char2num(LowerStr(symbol))-96 + 1
		endif
		if (mm>0)
			mm = limit(mm,0,27)
			PopupMenu setAtomWyckoff,win=$win,mode=mm
		endif
	endif
	xtal.atom[m].x = xx
	xtal.atom[m].y = yy
	xtal.atom[m].z = zz
	xtal.atom[m].WyckoffSymbol = symbol[0]
	if (char2num(symbol)>32 && !(mult>0))
		mult = WyckoffMultiplicity(xtal.SpaceGroupID,symbol)
	endif
	if (mult>0)
		xtal.atom[m].mult = mult
	endif
	ControlInfo/W=$win setAtomOcc		;		xtal.atom[m].occ = V_Value
	ControlInfo/W=$win setValencePopup	;		xtal.atom[m].valence = str2num(S_Value)
	ControlInfo/W=$win setAtomDebyeT	;		xtal.atom[m].DebyeT = V_Value
	ControlInfo/W=$win setAtomBiso		;		xtal.atom[m].Biso = V_Value * 0.01	// Panel uses Angstrom^2, xtal uses nm^2
	ControlInfo/W=$win setAtomUiso		;		xtal.atom[m].Uiso = V_Value * 0.01
	ControlInfo/W=$win setAtomU11		;		xtal.atom[m].U11 = V_Value * 0.01
	ControlInfo/W=$win setAtomU22		;		xtal.atom[m].U22 = V_Value * 0.01
	ControlInfo/W=$win setAtomU33		;		xtal.atom[m].U33 = V_Value * 0.01
	ControlInfo/W=$win setAtomU12		;		xtal.atom[m].U12 = V_Value * 0.01
	ControlInfo/W=$win setAtomU13		;		xtal.atom[m].U13 = V_Value * 0.01
	ControlInfo/W=$win setAtomU23		;		xtal.atom[m].U23 = V_Value * 0.01

	Wave/T atomNameList = $GetUserData(win,"","listWave")
	if (WaveExists(atomNameList))				// update name in atomNameList from label
		atomNameList[m] = xtal.atom[m].name
	endif

//	print "transfer values from panel to xtal and atomNameList"
//	PopupMenu useThermalPopUp,mode=5,popvalue="anisotropic Uij",value= #"\"NO Thermal;Debye T;B-isotropic;U-isotropic;anisotropic Uij\""

	StructPut/S/B=2 xtal, strStruct				// reset copy of xtal in Panel, saving in global string
	return 0
End

// ============================= End of Atom Editing ==============================


Function print_crystalStructStr(strStruct)	// prints the contents of a crystalStructStr in readable form
	String strStruct
	STRUCT crystalStructure xtal					// holds values from strStruct
	StructGet_xtal(strStruct,xtal)				// retrieve xtal from the string strStruct
	print_crystalStructure(xtal)
End


Function print_crystalStructure(xtal, [brief])			// prints out the value in a crystalStructure
	STRUCT crystalStructure &xtal				// this sruct is printed in this routine
	Variable brief										// if true, print short version
	brief = ParamIsDefault(brief) || numtype(brief) ? 0 : brief

	Variable i,N=xtal.N
	if (WhichListItem("LatticePanelButtonProc",GetRTStackInfo(0))>=0)
		printf "%sshowCrystalStructure()\r",BULLET
	endif
	if (strlen(xtal.desc))
		printf "'%s'   \t",xtal.desc
	else
		printf "\t\t\t"
	endif
	Variable SG = xtal.SpaceGroup, idNum=xtal.SpaceGroupIDnum
	String id = xtal.SpaceGroupID
	printf "Space Group=%s   %s   %s       Vc = %g (nm^3)",id,StringFromList(latticeSystem(id),LatticeSystemNames), getHMboth(idNum),xtal.Vc
	if (xtal.density>0)
		printf "       density = %g (g/cm^3)",xtal.density
	endif
	if (xtal.alphaT>0)
		printf "     coefficient of thermal expansion is %g\r",xtal.alphaT
	else
		printf "\r"
	endif
	printf "lattice constants  { %.15gnm, %.15gnm, %.15gnm,   %.15g%s, %.15g%s, %.15g%s }",xtal.a,xtal.b,xtal.c,xtal.alpha,DEGREESIGN,xtal.beta,DEGREESIGN,xtal.gam,DEGREESIGN
	if (numtype(xtal.Temperature)==0)
		Variable Temperature = xtal.Temperature
		String unit = "C"
		if (Temperature < -100)					// display in Kelvin
			unit = "K"
			Temperature = ConvertTemperatureUnits(Temperature,"C",unit)		// xtal.Temperature is always C
		endif
		printf ",   Temperature = %g %s%s\r",Temperature,DEGREESIGN,unit
	else
		printf "\r"
	endif
	if (SG<=167 && SG>=143 && (abs(90-xtal.alpha)+abs(90-xtal.beta)+abs(120-xtal.gam))<1e-6)	// also show rhombohedral lattice constants
		Variable aRhom = sqrt(3*(xtal.a)^2 + (xtal.c)^2)/3
		Variable alphaRhom = 2*asin(1.5/sqrt(3+(xtal.c/xtal.a)^2)) * 180/PI
		printf "   Rhombohedral constants, aRhom = %.8g nm,   alpha(Rhom) = %.8g%s\r",aRhom,alphaRhom,DEGREESIGN
	endif
	if (strlen(xtal.sourceFile)>1)
		printf "xtal read from file: \"%s\"\r",xtal.sourceFile
	endif

	if (N<1)
		print "No Atoms Defined"
	else
		String formula = MinimalChemFormula(xtal)				// this executes reMakeAtomXYZs(xtal), so don't do it again
		if (strlen(formula)>0)
			printf "atom type locations:\t chemical formula = \"%s\"\r",formula
		else
			printf "atom type locations:\r"
		endif
		Variable mult, itemp
		for (i=0;i<xtal.N;i+=1)									// loop over the defined atoms
			String vstr=""
			if (xtal.atom[i].valence)
				sprintf vstr, ", %+d",xtal.atom[i].valence
			endif
			if (strlen(xtal.atom[i].WyckoffSymbol))
				printf "     %s (Z=%g%s)\t%s\t{%g,  %g,  %g}",xtal.atom[i].name,xtal.atom[i].Zatom,vstr,xtal.atom[i].WyckoffSymbol,xtal.atom[i].x,xtal.atom[i].y,xtal.atom[i].z
			else
				printf "     %s (Z=%g%s)\t\t{%g,  %g,  %g}",xtal.atom[i].name,xtal.atom[i].Zatom,vstr,xtal.atom[i].x,xtal.atom[i].y,xtal.atom[i].z
			endif
			if (xtal.atom[i].occ != 1)
				printf "    occ = %g",xtal.atom[i].occ
			endif

			mult = DimSize($("root:Packages:Lattices:atom"+num2istr(i)),0)
			mult = !(mult>0) ? xtal.atom[i].mult : mult
			if (mult>1)
				printf "\tmultiplicity = %d",mult
			endif

			itemp = atomThermalInfo(xtal.atom[i])
			if (itemp==1)					// Thermal vibration information, print at most one kind of info
				printf "\t    Debye Temperature = %g K",xtal.atom[i].DebyeT
			elseif (itemp==2)
				printf "    Isotropic B = %g (%s^2)",(xtal.atom[i].Biso * 100),ARING
			elseif (itemp==3)
				printf "    Isotropic U = %g (%s^2)",(xtal.atom[i].Uiso * 100),ARING
			elseif (itemp>=4)
				printf "    Anisotropic U: U11=%+g, U22=%+g, U33=%+g",(xtal.atom[i].U11 *100),(xtal.atom[i].U22 *100),(xtal.atom[i].U33 *100)
				if (itemp==5)
					printf "    U12=%+g, U13=%+g, U23=%+g",(xtal.atom[i].U12 *100),(xtal.atom[i].U13 *100),(xtal.atom[i].U23 *100)
				endif
				printf " ("+ARING+"^2)"
			endif
			printf "\r"
		endfor

		String list = OverOccupyList(xtal)						// check if some sites have occ>1
		if (strlen(list))
			print "*************************************************"
			print "\t\t Found over occupied sites:"
			print OverOccupyList2Str(list)
			print "*************************************************"
		endif

		Variable netCharge = NetChargeCell(xtal)
		if (netCharge != 0)
			print "\t*************************************************"
			printf "\t\tnetCharge in cell = %g\r",netCharge
			print "\t*************************************************"
		else
			print "\t\tCharge Neutral"
		endif
		if (brief)
			return 0
		endif

		if (xtal.Nbonds > 0)
			printf "defined bonds types (nm):\r"
			for (i=0;i<xtal.Nbonds;i+=1)
				Make/N=(xtal.bond[i].N)/FREE/D lens
				lens = xtal.bond[i].len[p]
				printf "     %s  <-->  %s:  %s\r",xtal.bond[i].label0,xtal.bond[i].label1,vec2str(lens,bare=1)
			endfor
		endif
		printf "\rindividual atom positions:\rtype	frac xyz\r"
		Variable j
		for (i=0;i<xtal.N;i+=1)
			Wave wa = $("root:Packages:Lattices:atom"+num2istr(i))
			if (!WaveExists(wa))
				break
			endif
			for (j=0;j<DimSize(wa,0);j+=1)
				printf "%d\t\t{%g,  %g,  %g}\r"i,wa[j][0],wa[j][1],wa[j][2]
			endfor
		endfor
	endif
	if (numtype(xtal.Unconventional00)==0)					// Unconventional exists, transform the lattice by it
		printf "\rusing an Uncoventional Cell, the transform is:\r"
		printf "\t\t%+6.3f\t\t%+6.3f\t\t%+6.3f\r",xtal.Unconventional00,xtal.Unconventional01,xtal.Unconventional02
		printf "\tT =\t%+6.3f\t\t%+6.3f\t\t%+6.3f\r",xtal.Unconventional10,xtal.Unconventional11,xtal.Unconventional12
		printf "\t\t%+6.3f\t\t%+6.3f\t\t%+6.3f\r",xtal.Unconventional20,xtal.Unconventional21,xtal.Unconventional22
	endif
	printf "\r"
	printf "					a				b				c   (nm)\r"
	printf "direct	\t%+10.7f\t%+10.7f\t%+10.7f\r",xtal.a0,xtal.b0,xtal.c0
	printf "lattice\t%+10.7f\t%+10.7f\t%+10.7f\r",xtal.a1,xtal.b1,xtal.c1
	printf " \t\t\t%+10.7f\t%+10.7f\t%+10.7f\r",xtal.a2,xtal.b2,xtal.c2
	printf "\r					a*				b*				c*  (1/nm)\r"
	printf "recip	\t%+10.6f\t%+10.6f\t%+10.6f\r",xtal.as0,xtal.bs0,xtal.cs0
	printf "lattice\t%+10.6f\t%+10.6f\t%+10.6f\r",xtal.as1,xtal.bs1,xtal.cs1
	printf "	\t\t\t%+10.6f\t%+10.6f\t%+10.6f\r",xtal.as2,xtal.bs2,xtal.cs2
	if (strlen(xtal.hashID)<2)
		printf "hash id = '%s'\r",xtal.hashID
	endif
End


Static Function/S xtalHashID(xtalIN)		// Calculates the hashID for xtal
	STRUCT crystalStructure &xtalIN		// This does NOT get modified

	STRUCT crystalStructure xtal
	xtal = xtalIN								// not changing xtalIN, was:  copy_xtal(xtal,xtalIN)

	Make/FREE waveStruct
	StructPut xtal, waveStruct

	String strAll="", hex
	Variable i, m, N=numpnts(waveStruct)-HASHID_LEN
	for (i=0;i<N;i+=1)
		sprintf hex,"%X", waveStruct[i]
		strAll += hex
	endfor

	for (m=0;m<xtal.N;m+=1)				// for each defined atom type
		StructPut xtal.atom[m], waveStruct
		for (i=0;i<numpnts(waveStruct);i+=1)
			sprintf hex,"%X", waveStruct[i]
			strAll += hex
		endfor
	endfor

	return hash(strAll,1)
End


Static Function LatticeBad(xtal,[atomsToo])
	STRUCT crystalStructure &xtal
	Variable atomsToo			// also check atom definitions
	atomsToo = ParamIsDefault(atomsToo) ? 0 : atomsToo
	atomsToo = numtype(atomsToo) ? 0 : !(!atomsToo)

	Variable i, vibrate=0, bad=0, haveDebyeT=0
	bad += numtype(xtal.a + xtal.b + xtal.c)
	bad += !(xtal.a > 0) || !(xtal.b > 0) || !(xtal.c > 0)
	bad += numtype(xtal.alpha + xtal.beta + xtal.gam)
	bad += !(xtal.alpha > 0) || !(xtal.beta > 0) || !(xtal.gam > 0)
	bad += xtal.alpha >=180 || xtal.beta >= 180 || xtal.gam >= 180
	bad += !isValidSpaceGroup(xtal.SpaceGroup)
	for (i=0;i<xtal.N && atomsToo;i+=1)
		bad += atomBAD(xtal.atom[i])
		vibrate = vibrate || atomThermalInfo(xtal.atom[i])
		haveDebyeT = haveDebyeT || (xtal.atom[i].DebyeT)>0
	endfor
	if (atomsToo)
		bad += (xtal.Vibrate && !vibrate) || (!(xtal.Vibrate) && vibrate)		// if xtal.Vibrate, then one of the atoms must have Thermal Info
		bad += (xtal.haveDebyeT && !haveDebyeT) || (!(xtal.haveDebyeT) && haveDebyeT)
	endif
	return (bad>0)
End
//
Static Function atomBAD(atom)
	Struct atomTypeStructure &atom
	Variable bad=0
	bad += strlen(atom.name)<1
	bad += numtype(atom.x+atom.y+atom.z)>0
	bad += numtype(atom.occ) || atom.occ<0 || atom.occ>1
	bad += atom.Zatom != limit(round(atom.Zatom),1,92)
	bad += !(abs(atom.valence)<10)					// a valence of 10 is too big
	bad += atom.DebyeT < 0								// Debye, and isotropic U or B must be positive
	bad += atom.Biso < 0 || numtype(atom.Biso)==1
	bad += atom.Uiso < 0 || numtype(atom.Uiso)==1

	// Note, anisotropic Uij can be negative
	if (numtype(atom.U11)==0 || numtype(atom.U22)==0 || numtype(atom.U33)==0)	// using anisotropic Uii
		// if one is Uii valid, all symmetric Uii must be valid
		bad += numtype(atom.U11 + atom.U22 + atom.U33)!=0
		if (numtype(atom.U12)==0 || numtype(atom.U13)==0 || numtype(atom.U23)==0)	// using anisotropic Uij
			// if one Uij is valid, all asymmetric Uij must be valid
			bad += numtype(atom.U12 + atom.U13 + atom.U23)!=0
		endif
	endif
	return (bad>0)
End
//
ThreadSafe Static Function atomThermalInfo(atom,[T_K])	// True if atom contains thermal information of some kind
	Struct atomTypeStructure &atom
	Variable T_K											// optional temperature
	T_K = ParamIsDefault(T_K) ? 300.0 : T_K		// only consider T_K if it is passed (0 is valid)
	// returns type of thermal info
	//		0 = No thermal information for this atom
	//		1 = Debye T
	//		2 = Biso
	//		3 = Uiso
	//		4 = U11, U22, U22
	//		5 = U11, U22, U22, U12, U13, U23

	Variable thetaM = atom.DebyeT
	Variable U11_OK = numtype(atom.U11 + atom.U22 + atom.U33) == 0
	Variable U12_OK = numtype(atom.U12 + atom.U13 + atom.U23) == 0

	Variable vibrate=0
	if (thetaM>0 && T_K>0)			// have a valid and Debye Temperature and Temperature
		vibrate = 1
	elseif (atom.Biso > 0)			// have a valid Biso
		vibrate = 2
	elseif (atom.Uiso > 0)			// have a valid Uiso
		vibrate = 3
	elseif (U11_OK && U12_OK)		// both U11 and U12 are valid
		vibrate = 5
	elseif (U11_OK)					// have valid U11, but U12 is bad
		vibrate = 4
	endif
	return vibrate
End


// forces lattice constants to match the Space Group number (e.g. for cubic, forces b and c to be a, and all angles 90)
Static Function ForceLatticeToStructure(xtal)
	STRUCT crystalStructure &xtal					// this sruct is set in this routine

	// Make sure that the SpaceGroup, SpaceGroupID, and SpaceGroupIDnum all agree
	Variable SG, idNum
	String id = xtal.SpaceGroupID	// first try to use the id to identify the material			
	if (isValidSpaceGroupID(id))
		SG = str2num(id)
		idNum = SpaceGroupID2num(id)
	else
		SG = xtal.SpaceGroup			// failed with SpaceGroupID, try just SpaceGroup
		id = FindDefaultIDforSG(SG)	
		idNum = SpaceGroupID2num(id)
		xtal.SpaceGroupID = id			// set id & idNum using values from SpaceGroup
		xtal.SpaceGroupIDnum = idNum
	endif
	if (!isValidSpaceGroupIDnum(idNum))	// SpaceGroupIDnum must be in range [1,530]
		DoAlert 0, "invalid Space Group ID "+id+"     or SpaceGroup #"+num2str(SG)
		return 1
	endif

	//	system			  SG			   idNum
	//	Triclinic		[1,2]				[1,2]					a,b,c,alpha,beta,gamma
	//	Monoclinic		[3,15]			[3,107]				a,b,c,gamma
	//	Orthorhombic	[16,74]			[108,348]			a,b,c
	//	Tetragonal		[75,142]			[349,429]			a,c
	//	Trigonal			[143,167]		[430,461]			Hex: a,c   Rhom: a,alpha
	//	Hexagonal		[168,194]		[462,488]			a,c
	//	Cubic				[195,230]		[489,530]			a
	//	
	//	for Rhombohedral		a,alpha
	//	SG =		[146,148,155,160,161,166,167]
	//	idNum =	[434,437,445,451,453,459,461]

	// force the lattice constant to confomr to system
	Variable system = latticeSystem(id)
	if (system == CUBIC)					// Cubic space groups [195,230]
			xtal.b = xtal.a
			xtal.c = xtal.a
			xtal.alpha=90;  xtal.beta=90;  xtal.gam=90

	elseif (system == HEXAGONAL)			// Hexagonal space groups [168,194]
			xtal.b = xtal.a
			xtal.alpha=90;  xtal.beta=90;  xtal.gam=120

	elseif (system == TRIGONAL)			// Trigonal space groups [143,167]
		// generally use hexagonal cell, for rhomohedral may use rhomohedral cell
		Variable useRhom = 0
		if (isRhombohedral(SG))
			useRhom = abs(xtal.alpha - xtal.beta)<0.1 && abs(xtal.alpha - xtal.gam)<0.1
			useRhom = useRhom | ( strsearch(id,":R",0,2)>0 )
		endif
		if (useRhom)
			// looks like a rhombohedral structure using rhombohedral lattice constants
			xtal.b = xtal.a					// set to rhombohedral a=b=c, alpha=beta=gamma
			xtal.c = xtal.a
			xtal.beta = xtal.alpha
			xtal.gam = xtal.alpha
		else										// set to hexagonal lattice constants
			xtal.b = xtal.a
			xtal.alpha=90  ;  xtal.beta=90  ;  xtal.gam=120
		endif

	elseif (system == TETRAGONAL)		// Tetragonal space groups [75,142]
		xtal.b = xtal.a
		xtal.alpha=90;  xtal.beta=90;  xtal.gam=90

	elseif (system == ORTHORHOMBIC)		// Orthorhombic space groups [16,74]
		xtal.alpha=90;  xtal.beta=90;  xtal.gam=90

	elseif (system == MONOCLINIC)		// Monoclinic space groups [3,15]
		xtal.alpha=90;  xtal.gam=90

	elseif (system == TRICLINIC)			// Triclinic space groups 1,2]
		// nothing to fix
	else
		return 1									// system is invalid, NOT valid
	endif

	// check that lattice constants appear valid (non NaN or negatives or zero...)
	String str
	if (xtal.a<=0 || xtal.b<=0 || xtal.c<=0 || numtype(xtal.a+xtal.b+xtal.c))
		sprintf str,"invalid, (a,b,c) = (%g,%g,%g)",xtal.a,xtal.b,xtal.c
		DoAlert 0, str
		return 1
	endif
	Variable alpha=xtal.alpha,bet=xtal.beta,gam=xtal.gam
	if (alpha<=0 || bet<=0 || gam<=0 || alpha>=180 || bet>=180 || gam>=180 || numtype(alpha+bet+gam))
		sprintf str,"invalid, (alpha,beta,gam) = (%g,%g,%g)",alpha,bet,gam
		DoAlert 0, str
		return 1
	endif
	setDirectRecip(xtal)									// update Vc, direct and recip, density, and also calculates atom positions
	CleanOutCrystalStructure(xtal)
	xtal.Vibrate = xtalVibrates(xtal)					// True if some Thermal vibration info present in xtal
	xtal.haveDebyeT = xtalHasDebye(xtal)				// True if some one of the atoms has a Debye Temperature
	return 0
End
//
#ifndef OLD_LATTICE_ORIENTATION		// added at version<5.00, with a_Parallel_X, a||x and c*||z
Static Function setDirectRecip(xtal)					// set direct and recip lattice vectors from a,b,c,..., also calculates Vc & density
	STRUCT crystalStructure &xtal
	// Although not used here, note that the following also works:
	//		MatrixOP recipLatice = 2*PI * (Inv(directLattice))^t
	//		MatrixOP directLattice = 2*PI * Inv(recipLatice^t)
	//		Vc = MatrixDet(directLattice),    VcRecip = MatrixDet(recipLatice)
	//
	//	see:		https://en.wikipedia.org/wiki/Fractional_coordinates
	//  and		International Tables (2006) Vol. B, chapter 3.3 page 360
	//
	Variable a=xtal.a, b=xtal.b, c=xtal.c
	Variable ca = cos((xtal.alpha)*PI/180), cb = cos((xtal.beta)*PI/180)
	Variable cg = cos((xtal.gam)*PI/180), sg = sin((xtal.gam)*PI/180)
	Variable phi = sqrt(1.0 - ca*ca - cb*cb - cg*cg + 2*ca*cb*cg)	// = Vc/(a*b*c)
	Variable Vc = a*b*c * phi							// volume of unit cell
	Variable pv = (2*PI) / Vc							// used to scale reciprocal lattice

	Variable a0,a1,a2,  b0,b1,b2,  c0,c1,c2		//components of the direct lattice vectors
	if (!isRhombohedralXtal(xtal))
		// conventional choice International Tables (2006) Vol. B, chapter 3.3 page 360
		a0=a				; a1=0						; a2=0
		b0=b*cg			; b1=b*sg					; b2=0
		c0=c*cb			; c1=c*(ca-cb*cg)/sg	; c2=c*phi/sg	// z || c* (or axb)
	else
		// Rhombohedral cell choice International Tables (2006) Vol. B, chapter 3.3 page 360 (top of second column)
		// a,b,c are symmetric about the 111 direction, or (a+b+c) == {x,x,x}
		Variable pp=sqrt(1+2*ca), qq=sqrt(1-ca)
		Variable pmq=(a/3)*(pp-qq), p2q=(a/3)*(pp+2*qq)
		a0=p2q			; a1=pmq					; a2=pmq
		b0=pmq			; b1=p2q					; b2=pmq
		c0=pmq			; c1=pmq					; c2=p2q
	endif

	xtal.Vc = Vc
	xtal.a0=a0		; xtal.a1=a1	; xtal.a2=a2	// assign to xtal values
	xtal.b0=b0		; xtal.b1=b1	; xtal.b2=b2
	xtal.c0=c0		; xtal.c1=c1	; xtal.c2=c2
	xtal.as0=(b1*c2-b2*c1)*pv	; xtal.as1=(b2*c0-b0*c2)*pv	; xtal.as2=(b0*c1-b1*c0)*pv	// (b x c)*2PI/Vc
	xtal.bs0=(c1*a2-c2*a1)*pv	; xtal.bs1=(c2*a0-c0*a2)*pv	; xtal.bs2=(c0*a1-c1*a0)*pv	// (c x a)*2PI/Vc
	xtal.cs0=(a1*b2-a2*b1)*pv	; xtal.cs1=(a2*b0-a0*b2)*pv	; xtal.cs2=(a0*b1-a1*b0)*pv	// (a x b)*2PI/Vc

	Variable allZero = abs(xtal.Unconventional00)+abs(xtal.Unconventional01)+abs(xtal.Unconventional02)
	allZero += abs(xtal.Unconventional10)+abs(xtal.Unconventional11)+abs(xtal.Unconventional12)
	allZero += abs(xtal.Unconventional20)+abs(xtal.Unconventional21)+abs(xtal.Unconventional22)
	xtal.Unconventional00 = allZero==0 ? NaN : xtal.Unconventional00
	if (numtype(xtal.Unconventional00)==0 && xtal.Unconventional00>-100 && xtal.Unconventional00<100)
		Make/N=(3,3)/O/D root:Packages:Lattices:Unconventional
		Wave Unconventional=root:Packages:Lattices:Unconventional
		Unconventional[0][0]=xtal.Unconventional00;	Unconventional[0][1]=xtal.Unconventional01;	Unconventional[0][2]=xtal.Unconventional02
		Unconventional[1][0]=xtal.Unconventional10;	Unconventional[1][1]=xtal.Unconventional11;	Unconventional[1][2]=xtal.Unconventional12
		Unconventional[2][0]=xtal.Unconventional20;	Unconventional[2][1]=xtal.Unconventional21;	Unconventional[2][2]=xtal.Unconventional22
	else
		KillWaves/Z root:Packages:Lattices:Unconventional
	endif

	xtal.density = densityOfCrystalStructure(xtal)
	if (xtal.N==0)								// no atom defined, make one dummy atom
		xtal.N = 1
		xtal.atom[0].x = 0
		xtal.atom[0].y = 0
		xtal.atom[0].z = 0
		xtal.atom[0].Zatom = 1				// Z of the atom
		xtal.atom[0].name = "H1"
		xtal.atom[0].occ = 1
		xtal.atom[0].valence = 0
	endif
	return 0
End
#else													// OLD_LATTICE_ORIENTATION was defined, this gives c||z which is non-standard
Static Function setDirectRecip(xtal)					// set direct and recip lattice vectors from a,b,c,..., also calculates Vc & density
	STRUCT crystalStructure &xtal
	// Although not used here, note that the following also works:
	//		MatrixOP recipLatice = 2*PI * (Inv(directLattice))^t
	//		MatrixOP directLattice = 2*PI * Inv(recipLatice^t)
	//		Vc = MatrixDet(directLattice),    VcRecip = MatrixDet(recipLatice)
	Variable a=xtal.a, b=xtal.b, c=xtal.c
	Variable sa = sin((xtal.alpha)*PI/180), ca = cos((xtal.alpha)*PI/180)
	Variable cb = cos((xtal.beta)*PI/180), cg = cos((xtal.gam)*PI/180)
	Variable phi = sqrt(1.0 - ca*ca - cb*cb - cg*cg + 2*ca*cb*cg)	// = Vc/(a*b*c)
	xtal.Vc = a*b*c * phi								// volume of unit cell
	Variable pv = (2*PI) / (xtal.Vc)				// used for scaling
	Variable a0,a1,a2,  b0,b1,b2,  c0,c1,c2		//components of the direct lattice vectors
	a0=a*phi/sa		; a1=a*(cg-ca*cb)/sa	; a2=a*cb
	b0=0				; b1=b*sa		; b2=b*ca
	c0=0				; c1=0			; c2=c
	xtal.a0=a0		; xtal.a1=a1	; xtal.a2=a2
	xtal.b0=b0		; xtal.b1=b1	; xtal.b2=b2
	xtal.c0=c0		; xtal.c1=c1	; xtal.c2=c2
	xtal.as0=(b1*c2-b2*c1)*pv	; xtal.as1=(b2*c0-b0*c2)*pv	; xtal.as2=(b0*c1-b1*c0)*pv	// (b x c)*2PI/Vc
	xtal.bs0=(c1*a2-c2*a1)*pv	; xtal.bs1=(c2*a0-c0*a2)*pv	; xtal.bs2=(c0*a1-c1*a0)*pv	// (c x a)*2PI/Vc
	xtal.cs0=(a1*b2-a2*b1)*pv	; xtal.cs1=(a2*b0-a0*b2)*pv	; xtal.cs2=(a0*b1-a1*b0)*pv	// (a x b)*2PI/Vc

	Variable allZero = abs(xtal.Unconventional00)+abs(xtal.Unconventional01)+abs(xtal.Unconventional02)
	allZero += abs(xtal.Unconventional10)+abs(xtal.Unconventional11)+abs(xtal.Unconventional12)
	allZero += abs(xtal.Unconventional20)+abs(xtal.Unconventional21)+abs(xtal.Unconventional22)
	xtal.Unconventional00 = allZero==0 ? NaN : xtal.Unconventional00
	if (numtype(xtal.Unconventional00)==0 && xtal.Unconventional00>-100 && xtal.Unconventional00<100)
		Make/N=(3,3)/O/D root:Packages:Lattices:Unconventional
		Wave Unconventional=root:Packages:Lattices:Unconventional
		Unconventional[0][0]=xtal.Unconventional00;	Unconventional[0][1]=xtal.Unconventional01;	Unconventional[0][2]=xtal.Unconventional02
		Unconventional[1][0]=xtal.Unconventional10;	Unconventional[1][1]=xtal.Unconventional11;	Unconventional[1][2]=xtal.Unconventional12
		Unconventional[2][0]=xtal.Unconventional20;	Unconventional[2][1]=xtal.Unconventional21;	Unconventional[2][2]=xtal.Unconventional22
	else
		KillWaves/Z root:Packages:Lattices:Unconventional
	endif

	xtal.density = densityOfCrystalStructure(xtal)
	if (xtal.N==0)								// no atom defined, make one dummy atom
		xtal.N = 1
		xtal.atom[0].x = 0
		xtal.atom[0].y = 0
		xtal.atom[0].z = 0
		xtal.atom[0].Zatom = 1				// Z of the atom
		xtal.atom[0].name = "H1"
		xtal.atom[0].occ = 1
		xtal.atom[0].valence = 0
	endif
	return 0
End
#endif

Function/WAVE direct2LatticeConstants(direct)	// calculate lattice constants angles in degree
	// take three direct lattice vectors and return lattice constants as a free wave[6]
	Wave direct								// 3x3 matrix with direct lattice vectors, {a,b,c}
	Make/N=3/D/FREE a=direct[p][0], b=direct[p][1], c=direct[p][2]
	Make/N=6/D/FREE LatticeConstants=NaN
	LatticeConstants[0] = norm(a)	// a
	LatticeConstants[1] = norm(b)	// b
	LatticeConstants[2] = norm(c)	// c
	LatticeConstants[3] = acos(MatrixDot(b,c)/norm(b)/norm(c)) * 180/PI	// alpha
	LatticeConstants[4] = acos(MatrixDot(a,c)/norm(a)/norm(c)) * 180/PI	// beta
	LatticeConstants[5] = acos(MatrixDot(a,b)/norm(a)/norm(b)) * 180/PI	// gamma
	return LatticeConstants
End

Static Function xtalVibrates(xtal)						// True if some Thermal vibration info present in xtal (for any atom)
	STRUCT crystalStructure &xtal
	Variable i,vibrate
	for (i=0,vibrate=0; i < xtal.N; i+=1)
		vibrate = vibrate || atomThermalInfo(xtal.atom[i])
	endfor
	return vibrate
End

Static Function xtalHasDebye(xtal)
	STRUCT crystalStructure &xtal						// this sruct is filled  by this routine

	Variable i, N=xtal.N
	for (i=0;i<N;i+=1)
		if (xtal.atom[i].DebyeT > 0)
			return 1
		endif
	endfor
	return 0
End

Static Function CleanOutCrystalStructure(xtal)	// clean out all unused values and set to defaults
	STRUCT crystalStructure &xtal					// the structure to fix up
	Variable m, i, itemp, N=xtal.N, Nbonds=xtal.Nbonds

	if (N>=STRUCTURE_ATOMS_MAX || N<0)		// invalid number of atoms
		xtal.N = 0
		return 1
	elseif (Nbonds>=(2*STRUCTURE_ATOMS_MAX) || Nbonds<0)
		xtal.Nbonds = 0
		return 1
	elseif (strlen(xtal.desc)>99)
		String str = xtal.desc
		xtal.desc = str[0,99]
		return 1
	endif

	Variable Zatom
	for (m=0;m<N;m+=1)							// loop over the defined atoms
		Zatom = xtal.atom[m].Zatom
		xtal.atom[m].Zatom = (Zatom<1 || Zatom>103) ? 1 : round(Zatom)
		str = xtal.atom[m].name
		if (strlen(str)>59)
			xtal.atom[m].name = str[0,59]
		endif
		str = xtal.atom[m].WyckoffSymbol
		if (strlen(str)>1)
			xtal.atom[m].WyckoffSymbol = str[0]
		endif
		itemp = atomThermalInfo(xtal.atom[m])
		if (itemp!=1)
			xtal.atom[m].DebyeT = NaN
		endif
		if (itemp!=2)
			xtal.atom[m].Biso = NaN
		endif
		if (itemp!=3)
			xtal.atom[m].Uiso = NaN
		endif
		if (itemp<4)
			xtal.atom[m].U11 = NaN ; xtal.atom[m].U22 = NaN ; xtal.atom[m].U33 = NaN
		endif
		if (itemp!=5)
			xtal.atom[m].U12 = NaN ; xtal.atom[m].U13 = NaN ; xtal.atom[m].U23 = NaN
		endif
	endfor
	for (m=N;m<STRUCTURE_ATOMS_MAX;m+=1)	// loop over the unused atoms
		SetToDummyATOM(xtal.atom[m])
	endfor

	Variable Nlen
	for (m=0;m<Nbonds;m+=1)						// loop over the defined bonds
		Nlen = xtal.bond[m].N						// number of bond lengths
		Nlen = (Nlen>0 && Nlen<5) ? Nlen : 0	// Nlen must be in [0,4], otherwise set to 0
		if (strlen(xtal.bond[m].label0)>59 || strlen(xtal.bond[m].label1)>59)
			xtal.bond[m].label0 = ""				// label is too long, bad bond
			xtal.bond[m].label1 = ""
			Nlen = 0
		endif
		for (i=Nlen;i<5;i+=1)
			xtal.bond[m].len[i] = NaN				// unused bond lengths
		endfor
	endfor
	for (m=Nbonds;m<2*STRUCTURE_ATOMS_MAX;m+=1)	// loop over the unused bonds
		xtal.bond[m].label0 = ""
		xtal.bond[m].label1 = ""
		xtal.bond[m].N = 0
		for (i=0;i<5;i+=1)
			xtal.bond[m].len[i] = NaN				// unused bond lengths
		endfor
	endfor
	return 0
End


Function densityOfCrystalStructure(xtal)		// returns the density (g/cm^3)
	STRUCT crystalStructure &xtal					// this sruct is filled  by this routine
	Variable NA=6.02214199e23							// Avagadro's number
	String name

	reMakeAtomXYZs(xtal)
	Variable m, amuAll									// atomic mass of all atoms in cell
	for (m=0,amuAll=0; m<(xtal.N); m+=1)			// for each atom type
		name="root:Packages:Lattices:atom"+num2istr(m)
		Wave ww = $name
		amuAll += Element_amu(xtal.atom[m].Zatom)*DimSize(ww,0) * xtal.atom[m].occ
	endfor
	return (amuAll/NA)/(xtal.Vc * 1e-21)			// grams / cm^3
End

Static Function NetChargeCell(xtal)				// find the net charge in a cell (from valences), should be zero
	STRUCT crystalStructure &xtal

	Variable i, mult, occ, valence, netCharge=0
	for (i=0;i<xtal.N;i+=1)							// loop over the defined atoms
		mult = DimSize($("root:Packages:Lattices:atom"+num2istr(i)),0)
		mult = mult>0 ? mult : xtal.atom[i].mult
		valence = xtal.atom[i].valence
		occ = xtal.atom[i].occ > 0 ? xtal.atom[i].occ : 1
		netCharge += (mult>=1 && abs(valence)>0) ? mult*occ*valence : 0
	endfor
	return netCharge
End

Function get_dhkl(h,k,l,[T])
	Variable h,k,l
	Variable T										// Temperature, only used if 
	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))				//fill the lattice structure with test values
		DoAlert 0, "no lattice structure found"
		return NaN
	endif
	Variable alphaTbad = !(abs(xtal.alphaT)<0.1)	// if  |alphaT| > 0.1 then it must be wrong
	if (numtype(h+k+l))
		h = numtype(h) ? 0 : h
		k = numtype(k) ? 0 : k
		l = numtype(l) ? 2 : l
		Prompt h,"H"
		Prompt k,"K"
		Prompt l,"L"
		Prompt T,"Temperature, standard is 22.5, (C)"
		Variable askForT = (!ParamIsDefault(T) && !alphaTbad)
		if (askForT)
			DoPrompt "(hkl)",h,k,l,T
		else
			DoPrompt "(hkl)",h,k,l
		endif
		if (V_flag)
			return NaN
		endif
	endif
	if (numtype(h+k+l))
		return NaN
	endif
	Variable d, usingT = !alphaTbad && numtype(T)==0
	if (usingT)											// both alphaT and T are valid
		d = dSpacing(xtal,h,k,l,T=T)
	else
		d = dSpacing(xtal,h,k,l)
	endif
	if (ItemsInList(GetRTStackInfo(0))<2)
		Variable places = placesOfPrecision(xtal.a)
		if (usingT)										// both alphaT and T are valid
			printf "d['%s' (T=%g),  (%s)] = %.9g nm\r",xtal.desc,T,hkl2str(h,k,l),roundSignificant(d,places)
		else
			printf "d['%s',  (%s)] = %.9g nm\r",xtal.desc,hkl2str(h,k,l),roundSignificant(d,places)
		endif
	endif
	return d
End
//
ThreadSafe Function dSpacing(xtal,h,k,l,[T])		// returns d-spacing for the hkl (nm)
	STRUCT crystalStructure &xtal		// this sruct is set in this routine
	Variable h,k,l
	Variable T
	Variable xx,yy,zz
	xx = h*xtal.as0 + k*xtal.bs0 + l*xtal.cs0
	yy = h*xtal.as1 + k*xtal.bs1 + l*xtal.cs1
	zz = h*xtal.as2 + k*xtal.bs2 + l*xtal.cs2
	Variable d = 2*PI/sqrt(xx*xx + yy*yy + zz*zz)
	if (abs(xtal.alphaT)<0.1 && !ParamIsDefault(T))	// do if T passed, and valid alphaT
		T = limit(T,-273.15,Inf)				// limit T to > absolute zero
		d = d*(1+xtal.alphaT*(T-22.5))		// apply temperature correction
	endif
	return d
End

// returns d-spacing (nm) given (hkl) and lattice constants, just a local utility
Function dSpacingFromLatticeConstants(h,k,l,a,b,c,alpha,bet,gam)
	Variable h,k,l
	Variable a,b,c,alpha,bet,gam							// lattice constants, lengths in nm angles in degrees

	STRUCT crystalStructure xtal						// this sruct is only used locally
	xtal.a = a			;	xtal.b = b		;	xtal.c = c	// put values into structure
	xtal.alpha = alpha	;	xtal.beta = bet	;	xtal.gam = gam
	xtal.SpaceGroup = 1
	ForceLatticeToStructure(xtal)

	Variable k0,k1,k2, d					// 	k =  {k0,k1,k2} = h*as + k*bs + l*cs
	k0 = h*xtal.as0 + k*xtal.bs0 + l*xtal.cs0
	k1 = h*xtal.as1 + k*xtal.bs1 + l*xtal.cs1
	k2 = h*xtal.as2 + k*xtal.bs2 + l*xtal.cs2
	d = 2*PI/sqrt(k0*k0 + k1*k1 + k2*k2)	// d = 2*pi/|k|
	return d
End

Function hkl2Q(h,k,l, qvec,[normal])				// compute qvector for (h,k,l), returns |Q|
	Variable h,k,l
	Wave qvec
	Variable normal									// TRUE -> normalize Q, but always return full length
	normal = ParamIsDefault(normal) ? 0 : normal
	if (numtype(h+k+l) || !WaveExists(qvec))
		if (WaveExists(qvec))
			qvec = NaN
		endif
		return NaN
	endif
	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))				//fill the lattice structure with current values
		DoAlert 0, "No Lattice, please set one"
		return 1
	endif
	qvec[0] = h*xtal.as0+k*xtal.bs0+l*xtal.cs0	// qvec =  recip x hkl
	qvec[1] = h*xtal.as1+k*xtal.bs1+l*xtal.cs1
	qvec[2] = h*xtal.as2+k*xtal.bs2+l*xtal.cs2

	Variable Qlen = norm(qvec)
	if (normal)
		qvec /= Qlen
	endif
	return Qlen
End

Function angleBetweenHKLs(h1,k1,l1,  h2,k2,l2, [printIt])	// find angle between (h1,k1,l1) and (h2,k2,l2)
	Variable h1,k1,l1,h2,k2,l2
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt
	if (numtype(h1+k1+l1+h2+k2+l2) && printIt )
		h1 = numtype(h1) ? 0 : h1
		k1 = numtype(k1) ? 0 : k1
		l1 = numtype(l1) ? 2 : l1
		h2 = numtype(h2) ? 1 : h2
		k2 = numtype(k2) ? 1 : k2
		l2= numtype(l2) ? 1 : l2
		String hklStr1, hklStr2
		sprintf hklStr1,"%g %g %g",h1,k1,l1
		sprintf hklStr2,"%g %g %g",h2,k2,l2
		Prompt hklStr1,"(hkl) of the first Q"
		Prompt hklStr2,"(hkl) of the second Q"
		DoPrompt "hkl max & font size",hklStr1,hklStr2
		if (V_flag)
			return 1
		endif
		sscanf hklStr1, "%g %g %g",h1,k1,l1
		sscanf hklStr2, "%g %g %g",h2,k2,l2
		printIt = 1
	endif

	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))					//fill the lattice structure with current values
		DoAlert 0, "No Lattice, please set one"
		return 1
	endif
	Make/N=3/O/D hkl_angleBetweenHKLs1, hkl_angleBetweenHKLs2
	Wave q1=hkl_angleBetweenHKLs1, q2=hkl_angleBetweenHKLs2
	q1[0] = h1*xtal.as0+k1*xtal.bs0+l1*xtal.cs0		// q1 =  recip x hkl1
	q1[1] = h1*xtal.as1+k1*xtal.bs1+l1*xtal.cs1
	q1[2] = h1*xtal.as2+k1*xtal.bs2+l1*xtal.cs2
	q2[0] = h2*xtal.as0+k2*xtal.bs0+l2*xtal.cs0		// q2 =  recip x hkl2
	q2[1] = h2*xtal.as1+k2*xtal.bs1+l2*xtal.cs1
	q2[2] = h2*xtal.as2+k2*xtal.bs2+l2*xtal.cs2
	normalize(q1)
	normalize(q2)
	Variable angle = acos(limit(MatrixDot(q1,q2),-1,1))
	angle *= 180/PI
	KillWaves hkl_angleBetweenHKLs1, hkl_angleBetweenHKLs2
	if (printIt)
		printf "angle between (%s) and (%s) is %g%s\r",hkl2str(h1,k1,l1),hkl2str(h2,k2,l2),angle,DEGREESIGN
	endif
	return angle
End


Function/T findClosestHKL(dIN,[tolerance,usingQ])
	Variable dIN									// desired d(nm),   or perhaps Q(1/nm),  see usingQ
	Variable usingQ									// if true, then dIN is actually Q (1/nm)
	Variable tolerance								// optional tolerance used for more output
	dIN = numtype(dIN) ? NaN : dIN
	usingQ = ParamIsDefault(usingQ) ? 0 : usingQ
	usingQ = numtype(usingQ) ? 0 : usingQ
	tolerance = ParamIsDefault(tolerance) ? 0.0002 : tolerance
	tolerance = (numtype(tolerance) || tolerance<=0) ? 0.0002 : tolerance

	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))				//fill the lattice structure with test values
		DoAlert 0, "ERROR findClosestHKL()\rno lattice structure found"
		return ""
	endif

	Variable printIt = strlen(GetRTStackInfo(2))<1 || !ParamIsDefault(tolerance)
	if (!(dIN>0))									// ask if not enough info or too much
		Prompt usingQ,"d-spacing or Q",popup,"d (nm);Q (1/nm)"
		Prompt dIN,"known d-spacing (nm)  or Q (1/nm)"
		Prompt tolerance,"tolerance for checking d-spacings"
		usingQ += 1
		DoPrompt "d-spacing",dIN, usingQ,tolerance
		if (V_flag)
			return ""
		endif
		usingQ -= 1
		dIN = usingQ ? 2*PI/dIN : dIN				// if user supplied Q, then convert to d
		printIt = 1
	elseif (usingQ)
		dIN = 2*PI/dIN								// if user supplied Q, then convert to d
	endif
	if (!(dIN>0) || !(tolerance>0))
		DoAlert 0, "ERROR findClosestHKL()\rinvalid d-spacing"
		return ""									// invalid d-spacing
	endif

	Variable hmax = ceil(dSpacing(xtal,1,0,0)/dIN) + 1
	Variable kmax = ceil(dSpacing(xtal,0,1,0)/dIN) + 1
	Variable lmax = ceil(dSpacing(xtal,0,0,1)/dIN) + 1

	Variable Nalloc=100, Nsaved=0
	Make/N=(Nalloc,6)/FREE/D saved=NaN
	Variable/C Fhkl
	Variable h,k,l, dhkl, derr, err0=Inf
	for (l=0;l<=lmax;l=incrementIndex(l))
		for (k=0;k<=kmax;k=incrementIndex(k))
			for (h=0;h<=hmax;h=incrementIndex(h))
				Fhkl = Fstruct(xtal,h,k,l)
				if (magsqr(Fhkl)<1 && numtype(real(Fhkl))==0)	// only consider allowed hkl (if invalid, use it)
					continue
				endif
				dhkl = dSpacing(xtal,h,k,l)

				derr = abs(dhkl-dIN)
				if (printIt && derr<tolerance)
					printf "d(%s) = %g nm,  Q = %g (1/nm),  F = %g + i %g,  |F| = %g\r",hkl2str(h,k,l),dhkl,2*PI/dhkl,real(Fhkl),imag(Fhkl),cabs(Fhkl)
				endif
				if (derr<err0)							// found a new best one
					err0 = derr
					saved = NaN
					Nsaved = 0
				endif
				if (derr<=err0)						// a duplicate best one
					if (Nsaved>=Nalloc)
						Nalloc += 100
						Redimension/N=(Nalloc,-1) saved
					endif
					saved[Nsaved][0] = h	;	saved[Nsaved][1] = k	;	saved[Nsaved][2] = l
					saved[Nsaved][3] = dhkl
					saved[Nsaved][4] = real(Fhkl)	;	saved[Nsaved][5] = imag(Fhkl)
					Nsaved += 1
				endif

			endfor
		endfor
	endfor
	if (Nsaved<1)
		DoAlert 0, "ERROR findClosestHKL()\rNothing found?"
		return ""
	endif

	Redimension/N=(Nsaved,-1) saved
	Make/N=3/D/FREE hklj,hkli
	String str,list=""
	Variable i,j
	for (j=0;j<(Nsaved-1);j+=1)
		hklj = saved[j][p]
		sprintf str, "%g,%g,%g,%g,%g,%g;",saved[j][0],saved[j][1],saved[j][2],saved[j][3],saved[j][4],saved[j][5]
		list += str
		for (i=j+1;i<Nsaved;i+=1)
			hkli = -saved[i][p]
			MatrixOp/FREE/O dij = sum(abs(hklj-hkli))
			if (dij[0] < 0.1)
				DeletePoints/M=0 i, 1, saved
				i -= 1
				Nsaved -= 1
			endif
		endfor
	endfor

	if (printIt)
		Variable d0=saved[0][3], Fmag=sqrt(saved[0][4]^2+saved[0][5]^2)
		print " "
		printf "Closest reflection to  %g nm  is (1 of %d):\r",dIN,Nsaved
		printf "d(%s) = %g nm,  Q = %g (1/nm),  F = %g + i %g,  |F| = %g\r",hkl2str(saved[0][0],saved[0][1],saved[0][2]),d0,2*PI/d0,saved[0][4],saved[0][5],Fmag
	endif
	return list
End
//
ThreadSafe Static Function incrementIndex(i)	// this is used for looping symmetrically over + and - indicies (see findClosestHKL for example)
	Variable i
	if (i==0)			// for i==0, increment
		return 1
	elseif (i>0)
		return -i		// for i>0, make negative
	else
		return 1-i		// for i<0, increment and make positive
	endif
End



Function FillCrystalStructDefault(xtal)			// fill the crystal structure with 'current' values
	STRUCT crystalStructure &xtal					// returns 0 if something set, 1 is nothing done

	String strStruct=StrVarOrDefault(":crystalStructStr","")// set to values in current directory
	if (strlen(strStruct)<1)
		strStruct=StrVarOrDefault("root:Packages:Lattices:crystalStructStr","")	// try the values in Packages
	endif

	if (strlen(strStruct)>1)							// found structure information, load into xtal
		StructGet_xtal(strStruct,xtal)				// retrieve xtal from the string strStruct
	else														// last chance, load from generic defaults
		LoadPackagePreferences/MIS=1 "LatticeSym","LatticeSymPrefs",1,xtal
		if (V_flag)
			return 1											// did nothing, nothing found, give up
		endif
		if (!isValidLatticeConstants(xtal))
			SetToDummyXTAL(xtal)						// set to valid dummy values
		endif
		StructPut/S/B=2 xtal, strStruct				// keep a local copy after loading from PackagePreferences
		String/G root:Packages:Lattices:crystalStructStr = strStruct
	endif
	setDirectRecip(xtal)								// ensure that these are set
	CleanOutCrystalStructure(xtal)
	xtal.Vibrate = xtalVibrates(xtal)				// True if some Thermal vibration info present in xtal
	xtal.haveDebyeT = xtalHasDebye(xtal)			// True if some one of the atoms has a Debye Temperature
	return 0
End
//
Static Function SetToDummyXTAL(xtal)				// fill the crystal structure with valid dummy values
	STRUCT crystalStructure &xtal
	xtal.desc = "Dummy Structure"
	xtal.SpaceGroupID = "195"							// simple cubic
	xtal.SpaceGroupIDnum = SpaceGroupID2num(xtal.SpaceGroupID)
	xtal.SpaceGroup = 195
	xtal.a = 0.5		;	xtal.b = 0.5		;	xtal.c = 0.5
	xtal.alpha = 90	;	xtal.beta = 90	;	xtal.gam = 90
	xtal.Vc = 0.125
	xtal.Vibrate = 0	;	xtal.haveDebyeT = 0
	xtal.Nbonds = 0
	xtal.N = 1
	xtal.sourceFile=""
	xtal.hashID = ""

	xtal.N = 1
	SetToDummyATOM(xtal.atom[0])
	xtal.atom[0].name = "H1"
End
//
Static Function SetToDummyATOM(atom)		// fill the atom structure with valid dummy values
	STRUCT atomTypeStructure &atom
	atom.name = ""
	atom.Zatom = 1
	atom.x = 0	;	atom.y = 0	;	atom.z = 0
	atom.occ = 1
	atom.valence = 0
	atom.WyckoffSymbol = ""
	atom.mult = 1
	atom.DebyeT = NaN	;	atom.Biso = NaN	;	atom.Uiso = NaN
	atom.U11 = NaN	;	atom.U22 = NaN	;	atom.U33 = NaN	
	atom.U12 = NaN	;	atom.U13 = NaN	;	atom.U23 = NaN
End


Function UpdateCrystalStructureDefaults(xtal)		// save xtal as a string in local, Packages, and PackagePreferences
	STRUCT crystalStructure &xtal						// returns 0 if something set, 1 is nothing done

	xtal.hashID = xtalHashID(xtal)						// re-set hash function to identify associated waves
//	String strStruct
//	xtal.hashID = ""											// re-set hash function to identify associated waves
//	StructPut/S xtal, strStruct
//	xtal.hashID = hash(strStruct,1)

	Wave atom0 = root:Packages:Lattices:atom0
	if (WaveExists(atom0))
//		String wnote=StringByKey("ID",note(atom0),"=")
		String wnote=note(atom0)
		if (!stringmatch(StringByKey("ID",wnote,"="),xtal.hashID))
			Note/K atom0, ReplaceStringByKey("ID",wnote,"","=")
		endif
	endif

	// save as defaults for future use
	String strStruct
	StructPut/S/B=2 xtal, strStruct						// string to write to globals
	String/G root:Packages:Lattices:crystalStructStr=strStruct	// always save in the Packages
	if (!stringmatch(GetDataFolder(1),"root:"))	// don't save in root dataFolder
		String/G :crystalStructStr=strStruct			// usually save in local data folder too
	endif
	SavePackagePreferences/FLSH=1 "LatticeSym","LatticeSymPrefs",1,xtal
	return 0
//	if (exists(":crystalStructStr")==2)				// use the global in this data folder
//		SVAR crystalStructStr=:crystalStructStr
//		crystalStructStr = strStruct
//	elseif (exists("root:crystalStructStr")==2)	// global exists in root
//		SVAR crystalStructStr=root:crystalStructStr// update global in root
//		crystalStructStr = strStruct
//	else	
//		String/G root:Packages:Lattices:crystalStructStr=strStruct	// neither exist, make global in Packages
//	endif
End


Static Function/T OverOccupyList(xtalIN,[occMax,printIt])	// checks if sites in an xtal have occupancy > 1
	// returns string with those sites with occupancy > occMax, if all OK, then string is empty.
	// if occMax=1 (default) then only over occupied sites are listed, if occMax=-Inf, then all are listed
	// the output is a list of the form "x0,y0,z0,occ0;x1,y1,z1,oc1;..."
	STRUCT crystalStructure &xtalIN
	Variable occMax			// maximum occupancy for printing, only print if (occ > occMax)
	Variable printIt
	occMax = ParamIsDefault(occMax) ? 1 : occMax
	occMax = numtype(occMax)==2 ? 1 : occMax
	printIt = ParamIsDefault(printIt) ? 0 : printIt
	printIt = numtype(printIt) ? 0 : !(!printIt)

	STRUCT crystalStructure xtal
	xtal = xtalIN					// was:  copy_xtal(xtal,xtalIN)
	String str, out=""

	Variable Natoms=xtal.N		// nuber of atom types
	Variable iatom, occ, i
	Variable x0,y0,z0, dist
	for (iatom=0;iatom<Natoms;iatom+=1)
		x0 = xtal.atom[iatom].x
		y0 = xtal.atom[iatom].y
		z0 = xtal.atom[iatom].z
		if (numtype(x0+y0+z0))
			continue
		endif
		for (i=iatom,occ=0; i<Natoms; i+=1)	// now loop over iatom and all of the following atoms
			dist = sqrt((x0-xtal.atom[i].x)^2 + (y0-xtal.atom[i].y)^2 + (z0-xtal.atom[i].z)^2)
			if (dist < 0.01)							// within 0.01 of a cell, assume same position
				occ += xtal.atom[i].occ			// accumulate occupancy
				xtal.atom[i].x = NaN				// this atom has been used, remove from further consideration
			endif
		endfor
		if (occ>occMax)
			sprintf str, "%g,%g,%g,%g;",x0,y0,z0,occ
			out += str
		endif
	endfor

	if (printIt && strlen(out))
		print OverOccupyList2Str(out,occMax=occMax)
	endif
	return out
End
//
Static Function/T OverOccupyList2Str(list,[occMax])
	String list
	Variable occMax			// maximum occupancy for printing, only print if (occ > occMax), use occMax=-Inf for everything
	occMax = ParamIsDefault(occMax) ? 1 : occMax	// default is 1, only print if occ>1
	occMax = numtype(occMax)==2 ? 1 : occMax

	String out="", item, str
	Variable i, N=ItemsInList(list), x,y,z,occ
	for (i=0;i<N;i+=1)
		item = StringFromList(i,list)
		x = str2num(StringFromlist(0,item,","))
		y = str2num(StringFromlist(1,item,","))
		z = str2num(StringFromlist(2,item,","))
		occ = str2num(StringFromlist(3,item,","))
		if (numtype(x+y+z+occ)==0 && occ>occMax)
			sprintf str, "site {%g, %g, %g}, occupancy = %g",x,y,z,occ
			out += SelectString(strlen(out),"","\r") + str
		endif
	endfor
	return out
End
//	Function Test_OverOccupyList(occMax)
//		Variable occMax
//		STRUCT crystalStructure xtal
//		FillCrystalStructDefault(xtal)
//		String list = LatticeSym#OverOccupyList(xtal,occMax=occMax)
//		if (strlen(list))
//			print LatticeSym#OverOccupyList2Str(list,occMax=occMax)
//		else
//			print "Occupancy is OK"
//		endif
//	End



//	DEPRECATED	DEPRECATED	DEPRECATED	DEPRECATED
// The following TWO routines are DEPRECATED, USE THE SET LATTICE PANEL.  It is ONLY called by setLattice()
//
// THIS ROUTINE IS DEPRECATED, use the set lattice panel instead
//		available from MakeLatticeParametersPanel(""),  or from the "Set Crystal Structure..." menu item
//
// prompts the user for lattice info, and set local structure.  Store the information in a structure string in current datafolder
Function setLattice()
	String strName=""
	Variable SpaceGroup									//Space Group number, from International Tables
	Variable a,b,c,alpha,bet,gam
	String desc

	STRUCT crystalStructure xtal
	xtal.a = 0
	FillCrystalStructDefault(xtal)						//fill the lattice structure with default values
	a = xtal.a			;	b = xtal.b		;	c = xtal.c
	alpha = xtal.alpha	;	bet = xtal.beta	;	gam = xtal.gam
	SpaceGroup = xtal.SpaceGroup
	desc = xtal.desc

	if (GetLatticeConstants(xtal,SpaceGroup,strName,a,b,c,alpha,bet,gam,desc))
		return 1
	endif
	xtal.Unconventional00 = NaN
	xtal.sourceFile = ""
	CleanOutCrystalStructure(xtal)
	UpdateCrystalStructureDefaults(xtal)
	return 0
End


// gets lattice constants for known structures via dialogs,  it uses Space Group numbers from the International Tables
// known structures are:
//	FCC			225
//	BCC			229
//	Diamond		227
//	Perovskite		221
//	Simple Cubic	195
//	Hexagonal		194
//	Sapphire		167
//	Triclinic		1
//
//	DEPRECATED	DEPRECATED	DEPRECATED	DEPRECATED
// THIS ROUTINE IS DEPRECATED, use the set lattice panel instead.  This is ONLY called by setLattice()
Static Function GetLatticeConstants(xtal,SpaceGroup,structureName,a,b,c,alpha,bet,gam,desc)
	STRUCT crystalStructure &xtal				// this sruct is set in this routine
	Variable SpaceGroup								// same as for Internationl tables
	String structureName							// see structures below for list of valid names
	Variable a,b,c,alpha,bet,gam					// lattice constants, lengths in nm angles in degrees
	String desc											// description of lattice (i.e. the name)
	a = numtype(a)||a<=0 ? 0.405 : a			// default is aluminum
	b = numtype(b)||b<=0 ? a : b
	c = numtype(c)||c<=0 ? a : c
	alpha = numtype(alpha)||alpha<=0 ? 90 : alpha	// for invalid angles, use 90 (degree)
	bet = numtype(bet)||bet<=0 ? 90 : bet
	gam = numtype(gam)||gam<=0 ? 90 : gam
	//	Cubic				[195,230]	//	a
	//	Hexagonal		[168,194]	//	a,c
	//	Trigonal			[143,167]	//	a,alpha
	//	Tetragonal		[75,142]		//	a,c
	//	Orthorhombic	[16,74]		//	a,b,c
	//	Monoclinic		[3,15]		//	a,b,c,gamma
	//	Triclinic		[1,2]			//	a,b,c,alpha,beta,gamma
//	String item="",structures="FCC:225;BCC:229;Diamond:227;Perovskite:221;Simple Cubic:195;Hexagonal:194;Wurtzite (B4):186;Sapphire:167;Triclinic:1;Space Group #..."
	String item="",structures="FCC:225;BCC:229;Diamond:227;Perovskite:221;Simple Cubic:195;Hexagonal:194;Wurtzite (B4):186;Sapphire:167;Triclinic:1;Space Group #..."
	if (strlen(xtal.desc)>0 && xtal.SpaceGroup==SpaceGroup)
		structures = xtal.desc+":"+num2istr(SpaceGroup)+";"+structures
	endif
	Variable i,N=ItemsInList(structures)
	if (isValidSpaceGroup(SpaceGroup))				// valid SpaceGroup number
		structureName = ""
		for(i=0;i<N;i+=1)
			item = StringFromList(i,structures)
			if (SpaceGroup==str2num(StringFromList(1,item,":")))
				structureName = StringFromList(0,item,":")
			endif
		endfor
		if (strlen(structureName)<1)
			DoAlert 0, "unknown structure number "+num2str(SpaceGroup)
		endif
	endif
	if (strlen(structureName)>1)						// a structureName passed
		SpaceGroup = -1
		for(i=0;i<N;i+=1)
			item = StringFromList(i,structures)
			if (stringmatch(structureName,StringFromList(0,item,":")))
				SpaceGroup = str2num(StringFromList(1,item,":"))
			endif
		endfor
		if (SpaceGroup<1)
			DoAlert 0, "unknown structure name '"+structureName+"'"
		endif
	endif
	item = structureName+":"+num2istr(SpaceGroup)
	Prompt item, "lattice structure", popup, structures
	DoPrompt/HELP="" "lattice structure",item
	if (V_flag)
		return 1
	endif
	if (stringmatch(item,"Space Group*"))
		Prompt SpaceGroup,"Space Group number from International Tables [1,230]"
		DoPrompt "Space Group",SpaceGroup
		if (V_flag)
			return 1
		endif
	else
		SpaceGroup=str2num(StringFromList(1,item,":"))
		structureName = StringFromList(0,item,":")
	endif
	if (!isValidSpaceGroup(SpaceGroup))
		return 1
	endif

	Prompt a, "a (nm)"
	Prompt b, "b (nm)"
	Prompt c, "c (nm)"
	Prompt alpha, "alpha ("+DEGREESIGN+")"
	Prompt bet, "beta ("+DEGREESIGN+")"
	Prompt gam, "gamma ("+DEGREESIGN+")"

	if (SpaceGroup>=195)			// Cubic
		DoPrompt/HELP="" "Cubic lattice constants",a
		b=a;  c=a
		alpha=90;  bet=90;  gam=90
	elseif(SpaceGroup>=168)		// Hexagonal
		c = (a==c) ? 3*a : c
		DoPrompt/HELP="" "Hexagonal lattice constants",a,c
		b = a
		alpha=90;  bet=90;  gam=120
	elseif(SpaceGroup>=143)		// Trigonal
		Variable useHex=1
		if (isRhombohedral(SpaceGroup))
			Prompt useHex, "Use Hex axes", popup, "Hexagonal Axes;Rhombohedral Axes"
			DoPrompt "What kind of Axes?", useHex
			if (V_flag)
				return 1
			endif
			useHex = useHex == 1
		endif

		if (useHex)
			c = (a==c) ? 3*a : c
			DoPrompt/HELP="" "Trigonal lattice constants (hexagoanl cell)",a,c
			b = a
			alpha=90;  bet=90;  gam=120
		else
			DoPrompt/HELP="" "Rhombohedral lattice constants (rhombohedral cell)",a,alpha
			b=a  ;  c=a
			bet=alpha;  gam=alpha
		endif
	elseif(SpaceGroup>=75)			// Tetragonal
		DoPrompt/HELP="" "Tetragonal lattice constants",a,c
		b = a
		alpha=90;  bet=90;  gam=90
	elseif(SpaceGroup>=16)			// Orthorhombic
		DoPrompt/HELP="" "Orthorhombic lattice constants",a,b,c
		alpha=90;  bet=90;  gam=90
	elseif(SpaceGroup>=3)			// Monoclinic
		DoPrompt/HELP="" "Monoclinic lattice constants",a,b,c,gam
//		alpha=90;  bet=90
		alpha=90;  gam=90
	else								// Triclinic
		DoPrompt/HELP="" "Triclinic lattice constants",a,b,c,alpha,bet,gam
	endif

	xtal.a = a			;	xtal.b = b		;	xtal.c = c	// put values into structure
	xtal.alpha = alpha	;	xtal.beta = bet	;	xtal.gam = gam
	xtal.SpaceGroup=SpaceGroup
	xtal.SpaceGroupID = FindDefaultIDforSG(SpaceGroup)
	xtal.SpaceGroupIDnum = SpaceGroupID2num(xtal.SpaceGroupID)
	xtal.Unconventional00=NaN
	xtal.sourceFile = ""
	ForceLatticeToStructure(xtal)
	return 0
End

//	End of setting particular lattice constants
// =========================================================================
// =========================================================================


// =========================================================================
// =========================================================================
//	Start of lattice set panel

Function/T FillLatticeParametersPanel(strStruct,hostWin,left,top)
	String strStruct									// optional passed value of xtal structure, this is used if passed
	String hostWin										// name of home window
	Variable left, top									// offsets from the left and top

	NewDataFolder/O root:Packages:Lattices:PanelValues	// ensure that the needed data folders exist
	Variable new = !(NumVarOrDefault("root:Packages:Lattices:PanelValues:a",-1)>0)
	String/G root:Packages:Lattices:PanelValues:desc		// create, but don't fill in the values for the panel
	String/G root:Packages:Lattices:PanelValues:SpaceGroupID
	Variable/G root:Packages:Lattices:PanelValues:a
	Variable/G root:Packages:Lattices:PanelValues:b
	Variable/G root:Packages:Lattices:PanelValues:c
	Variable/G root:Packages:Lattices:PanelValues:alpha
	Variable/G root:Packages:Lattices:PanelValues:bet
	Variable/G root:Packages:Lattices:PanelValues:gam
	Variable/G root:Packages:Lattices:PanelValues:alphaT
	Variable/G root:Packages:Lattices:PanelValues:dirty		// =1 (xtal.sourceFile bad too),  =2 (values changed, but xtal.sourceFile is OK)
	Variable/G root:Packages:Lattices:PanelValues:T_C
	SVAR SpaceGroupID=root:Packages:Lattices:PanelValues:SpaceGroupID
	SVAR desc=root:Packages:Lattices:PanelValues:desc
	NVAR a=root:Packages:Lattices:PanelValues:a
	NVAR b=root:Packages:Lattices:PanelValues:b
	NVAR c=root:Packages:Lattices:PanelValues:c
	NVAR alpha=root:Packages:Lattices:PanelValues:alpha
	NVAR bet=root:Packages:Lattices:PanelValues:bet
	NVAR gam=root:Packages:Lattices:PanelValues:gam
	NVAR alphaT=root:Packages:Lattices:PanelValues:alphaT
	NVAR T_C=root:Packages:Lattices:PanelValues:T_C
	NVAR dirty=root:Packages:Lattices:PanelValues:dirty

	Variable/G root:Packages:Lattices:PanelValues:h
	Variable/G root:Packages:Lattices:PanelValues:k
	Variable/G root:Packages:Lattices:PanelValues:l
	Variable/G root:Packages:Lattices:PanelValues:dspace_nm
	Variable/G root:Packages:Lattices:PanelValues:Fr
	Variable/G root:Packages:Lattices:PanelValues:Fi
	String/G root:Packages:Lattices:PanelValues:crystalStructStr
	SVAR crystalStructStr = root:Packages:Lattices:PanelValues:crystalStructStr

	STRUCT crystalStructure xtal
	if (strlen(strStruct))								// start using the passed values
		StructGet_xtal(strStruct,xtal)				// found passed structure information, load into xtal
		CleanOutCrystalStructure(xtal)
		dirty = 2											// xtal.sourceFile should be OK
		crystalStructStr = strStruct
	elseif(new)												// no old values present, use usual defaults
		FillCrystalStructDefault(xtal)
		a=xtal.a  ;  b=xtal.b  ;  c=xtal.c
		alpha=xtal.alpha  ;  bet=xtal.beta  ;  gam=xtal.gam
		SpaceGroupID = xtal.SpaceGroupID
		alphaT=xtal.alphaT
		desc=xtal.desc
		T_C = xtal.Temperature
		T_C = numtype(T_C) || T_C<-273.14 ? 22.5 : T_C
		T_C = (xtal.haveDebyeT) ? T_C : NaN		// if no Debye Temperatuers, no temperature needed
		dirty = 0
		StructPut/S xtal, crystalStructStr
	endif

	SetWindow kwTopWin,userdata(LatticePanelName)=hostWin+"#LatticePanel"
	Variable showEnergy = (exists("Get_f")==6)
	Variable bottom = top+60+508 + (showEnergy ? 20 : 0)
	NewPanel/K=1/W=(left,top,left+221,bottom)/HOST=$hostWin
	ModifyPanel frameStyle=0, frameInset=0
	RenameWindow #,LatticePanel
	SetVariable setDesc,pos={4,13},size={211,18},title="name"
	SetVariable setDesc,help={"description of this lattice"},fSize=12,proc=LatticePanelParamProc
	SetVariable setDesc,limits={-inf,inf,0},value= root:Packages:Lattices:PanelValues:desc
	Button buttonFindSpaceGroup,pos={35,36},size={150,20},title="New Space Group...", proc=LatticePanelButtonProc
	TitleBox structureTitle,pos={24,63},size={82,15},title="\\JC"
	TitleBox structureTitle,fSize=14,frame=0,anchor= LC
	SetVariable set_a_nm,pos={45,93},size={128,18},proc=LatticePanelParamProc,title="a (nm)"
	SetVariable set_a_nm,fSize=12,format="%.7f"
	SetVariable set_a_nm,limits={0,inf,0},value= root:Packages:Lattices:PanelValues:a
	SetVariable set_b_nm,pos={45,113},size={128,18},proc=LatticePanelParamProc,title="b (nm)"
	SetVariable set_b_nm,fSize=12,format="%.7f"
	SetVariable set_b_nm,limits={0,inf,0},value= root:Packages:Lattices:PanelValues:b
	SetVariable set_c_nm,pos={45,133},size={128,18},proc=LatticePanelParamProc,title="c (nm)"
	SetVariable set_c_nm,fSize=12,format="%.7f"
	SetVariable set_c_nm,limits={0,inf,0},value= root:Packages:Lattices:PanelValues:c
	SetVariable set_a_nm,help={"lattice constant (nm), length of 'a' vector"}
	SetVariable set_b_nm,help={"lattice constant (nm), length of 'b' vector"}
	SetVariable set_c_nm,help={"lattice constant (nm), length of 'c' vector"}
	SetVariable set_alpha,help={"lattice constant (degree), angle between b & c vectors"}
	SetVariable set_beta,help={"lattice constant (degree), angle between a & c vectors"}
	SetVariable set_gamma,help={"lattice constant (degree), angle between a & b vectors"}

	SetVariable set_alpha,pos={45,163},size={100,17},proc=LatticePanelParamProc,format="%10.5f"
	SetVariable set_beta,pos={45,183},size={100,17},proc=LatticePanelParamProc,format="%10.5f"
	SetVariable set_gamma,pos={45,203},size={100,17},proc=LatticePanelParamProc,format="%10.5f"
	SetVariable set_alpha,limits={0,180,0},value= root:Packages:Lattices:PanelValues:alpha
	SetVariable set_beta,limits={0,180,0},value= root:Packages:Lattices:PanelValues:bet
	SetVariable set_gamma,limits={0,180,0},value= root:Packages:Lattices:PanelValues:gam
#if (IgorVersion()<7)
	SetVariable set_alpha,font="Symbol", fSize=14, title="a\260"	// "\260" = option-5, degree sign in Symbol font
	SetVariable set_beta,font="Symbol", fSize=14, title="b\260"
	SetVariable set_gamma,font="Symbol", fSize=14, title="g\260"
#else
	SetVariable set_alpha, fSize=12, title="\xCE\xB1"+DEGREESIGN
	SetVariable set_beta, fSize=12, title="\xCE\xB2"+DEGREESIGN
	SetVariable set_gamma, fSize=12, title="\xC9\xA3"+DEGREESIGN	// there are two gamma's, this (latin small letter gamma) looks better than "\xCE\xB3"
#endif

	Button buttonLatticeSave,pos={35,233},size={150,20},proc=LatticePanelButtonProc,title="Save"
	Button buttonLatticeSave,help={"Save these values as the current values"}
	Button buttonLatticeRevert,pos={35,258},size={150,20},proc=LatticePanelButtonProc,title="Revert"
	Button buttonLatticeRevert,help={"revert to the current default"}

	Button buttonLatticeFromFile,pos={35,283},size={150,20},proc=LatticePanelButtonProc,title="Load from a file"
	Button buttonLatticeFromFile,help={"Fill the values in this panel from a file"}
	Button buttonPrintLattice,pos={35,308},size={150,20},proc=LatticePanelButtonProc,title="Print Lattice to History"
	Button buttonPrintLattice,help={"print lattice values to the history"}
	Button buttonWriteLattice,pos={35,333},size={150,20},proc=LatticePanelButtonProc,title="Write Lattice to File"
	Button buttonWriteLattice,help={"write current lattice values to an xml file"}

	PopupMenu popupLatticeEdit,pos={35,358},size={150,20},proc=LatticeSym#LatticeEditPopMenuProc,title="Utility & Atom Edit"
	PopupMenu popupLatticeEdit,help={"Provides supprort for editing atoms and other crystalographic changes"}
	String mstr
	sprintf mstr, "\"Edit Atom Positions%s;Change Current Setting%s;Find Closest hkl%s;Space Group number <--> symmetry%s;Describe the Symmetry Operations%s;angle between two hkl's%s\"", HORIZ_ELLIPSIS,HORIZ_ELLIPSIS,HORIZ_ELLIPSIS,HORIZ_ELLIPSIS,HORIZ_ELLIPSIS,HORIZ_ELLIPSIS
	PopupMenu popupLatticeEdit,fSize=14,mode=0,value= #mstr

	Button buttonAtomView,pos={35,383},size={150,20},proc=LatticePanelButtonProc,title="Add Atom View..."
	Button buttonAtomView,help={"Provides supprort for viewing the lattice as a 3D Gizmo"}
	PopupMenu popupAtomView,pos={35,383},size={150,20},proc=AtomView#AtomViewPopMenuProc,title="Atom View...",disable= 1
	PopupMenu popupAtomView,help={"Provides supprort for viewing the lattice as a 3D Gizmo"}
	PopupMenu popupAtomView,fSize=14,mode=0,value= #"\"Make Cells of Atoms...;  Bond Info;  Gizmo of Atoms;    Atom Type at Cursor\""
	Button buttonPowderPattern,pos={35,408},size={150,20},proc=LatticePanelButtonProc,title="Add Powder Patterns..."
	Button buttonPowderPattern,help={"Provides supprort creating a powder pattern simulation for this lattice"}
	PopupMenu popupPowderPattern,pos={35,408},size={150,20},proc=powder#PowderPatternPopMenuProc,title="Powder Pattern...",disable= 1
	PopupMenu popupPowderPattern,help={"Provides supprort creating a powder pattern simulation for this lattice"}
	PopupMenu popupPowderPattern,fSize=14,mode=0,value= #"\"Calculate Powder Lines ...;  Make Powder Pattern from Lines ...;  Graph of Powder Pattern;  Graph of Powder Lines;  Table of Powder Lines\""

	SetVariable h_LatticeVar,pos={12,439},size={60,15},proc=LatticePanelParamProc,title="H",font="Lucida Grande",fSize=12
	SetVariable h_LatticeVar,value= root:Packages:Lattices:PanelValues:h
	SetVariable k_LatticeVar,pos={84,439},size={60,15},proc=LatticePanelParamProc,title="K",font="Lucida Grande",fSize=12
	SetVariable k_LatticeVar,value= root:Packages:Lattices:PanelValues:k
	SetVariable L_LatticeVar,pos={156,439},size={60,15},proc=LatticePanelParamProc,title="L",font="Lucida Grande",fSize=12
	SetVariable L_LatticeVar,value= root:Packages:Lattices:PanelValues:l
	ValDisplay d_LatticeDisp,pos={12,460},size={123,16},title="d(nm)",font="Lucida Grande",fSize=12,format="%.7f"
	ValDisplay d_LatticeDisp,limits={0,0,0},value=#"root:Packages:Lattices:PanelValues:dspace_nm"
	ValDisplay d_LatticeDisp,help={"d-spacing (nm) calculated using the lattice"},frame=0
	SetVariable T_LatticeVar,pos={143,459},size={70,18},proc=LatticePanelParamProc,title="T("+DEGREESIGN+"C)",font="Lucida Grande",fSize=12
	SetVariable T_LatticeVar,limits={-273.15,inf,0},value=root:Packages:Lattices:PanelValues:T_C
	SetVariable T_LatticeVar,help={"Temperature (C) used when Thermal factors are given"}
	ValDisplay Fr_3atticeDisp,pos={20,483},size={80,17},title="F =",value= #"root:Packages:Lattices:PanelValues:Fr"
	ValDisplay Fr_3atticeDisp,font="Lucida Grande",fSize=12,format="%.3f",limits={0,0,0},barmisc={0,1000}
	ValDisplay Fr_3atticeDisp,help={"real part of F(hkl) calulated from the lattice"},frame=0
	ValDisplay Fi_3atticeDisp,pos={104,483},size={80,17},title="+i",value= #"root:Packages:Lattices:PanelValues:Fi"
	ValDisplay Fi_3atticeDisp,font="Lucida Grande",fSize=12,format="%.3f",limits={0,0,0},barmisc={0,1000}
	ValDisplay Fi_3atticeDisp,help={"imag part of F(hkl) calulated from the lattice"},frame=0
	if (showEnergy)
		SetVariable setvarEnergy,pos={37,504},size={140,20},title="Energy (keV)"
		SetVariable setvarEnergy,font="Lucida Grande",fSize=12,format="%.4f"
		SetVariable setvarEnergy,limits={0,inf,0},value= root:Packages:Lattices:keV
		SetVariable setvarEnergy,proc=LatticePanelParamProc
	endif

	String subWin = GetUserData("","","LatticePanelName")
	UpdatePanelLatticeConstControls(subWin,SpaceGroupID)

	STRUCT WMSetVariableAction sva
	sva.win = subWin
	sva.eventCode = 2
	LatticePanelParamProc(sva)
	return "#LatticePanel"
End
//
Static Function UpdatePanelLatticeConstControls(subWin,SpaceGroupID)
	// update a,b,c, alpha,beta,gamma in LatticeSet Panel, changes who is enabled, not the values
	String subWin
	String SpaceGroupID

	if (!isValidSpaceGroupID(SpaceGroupID))
		return 1
	endif

	NVAR a=root:Packages:Lattices:PanelValues:a
	NVAR b=root:Packages:Lattices:PanelValues:b
	NVAR c=root:Packages:Lattices:PanelValues:c
	NVAR alpha=root:Packages:Lattices:PanelValues:alpha
	NVAR bet=root:Packages:Lattices:PanelValues:bet
	NVAR gam=root:Packages:Lattices:PanelValues:gam
	NVAR T_C=root:Packages:Lattices:PanelValues:T_C

	Variable SG = str2num(SpaceGroupID)		// first part of id is space group number
	String titleStr="\\JC"+SpaceGroupID+" "
	if (SG>=195)															// Cubic, a
		SetVariable set_a_nm,noedit=0,frame=1,win=$subWin	// enable
		SetVariable set_b_nm,noedit=1,frame=0,win=$subWin	// disable
		SetVariable set_c_nm,noedit=1,frame=0,win=$subWin
		SetVariable set_alpha,noedit=1,frame=0,win=$subWin
		SetVariable set_beta,noedit=1,frame=0,win=$subWin
		SetVariable set_gamma,noedit=1,frame=0,win=$subWin
		b=a  ;  c=a  ;  alpha=90 ; bet=90 ; gam=90
		titleStr += "Cubic"
	elseif (SG>=168)														// Hexagonal, a, c
		SetVariable set_a_nm,noedit=0,frame=1,win=$subWin	// enable
		SetVariable set_c_nm,noedit=0,frame=1,win=$subWin
		SetVariable set_b_nm,noedit=1,frame=0,win=$subWin	// disable
		SetVariable set_alpha,noedit=1,frame=0,win=$subWin
		SetVariable set_beta,noedit=1,frame=0,win=$subWin
		SetVariable set_gamma,noedit=1,frame=0,win=$subWin
		b=a ; alpha=90 ; bet=90 ; gam=120
		titleStr += "Hexagonal"
	elseif (isRhombohedral(SG) && !((abs(90-alpha)+abs(90-bet)+abs(120-gam))<1e-6))	// Rhombohedral, with rhombohedral cell
		SetVariable set_a_nm,noedit=0,frame=1,win=$subWin	// enable
		SetVariable set_alpha,noedit=0,frame=1,win=$subWin
		SetVariable set_b_nm,noedit=1,frame=0,win=$subWin	// disable
		SetVariable set_c_nm,noedit=1,frame=0,win=$subWin
		SetVariable set_beta,noedit=1,frame=0,win=$subWin
		SetVariable set_gamma,noedit=1,frame=0,win=$subWin
		b=a  ;  c=a
		bet=alpha  ;  gam=alpha
		titleStr += "Rhombohedral"
	elseif (SG>=143)														// Trigonal, with hexagonal cell
		SetVariable set_a_nm,noedit=0,frame=1,win=$subWin	// enable
		SetVariable set_c_nm,noedit=0,frame=1,win=$subWin
		SetVariable set_b_nm,noedit=1,frame=0,win=$subWin	// disable
		SetVariable set_alpha,noedit=1,frame=0,win=$subWin
		SetVariable set_beta,noedit=1,frame=0,win=$subWin
		SetVariable set_gamma,noedit=1,frame=0,win=$subWin
		b=a
		alpha=90  ;  bet=90  ;  gam=120
		titleStr += "Trigonal"
	elseif (SG>=75)														// Tetragonal, a, c
		SetVariable set_a_nm,noedit=0,frame=1,win=$subWin	// enable
		SetVariable set_c_nm,noedit=0,frame=1,win=$subWin
		SetVariable set_b_nm,noedit=1,frame=0,win=$subWin	// disable
		SetVariable set_alpha,noedit=1,frame=0,win=$subWin
		SetVariable set_beta,noedit=1,frame=0,win=$subWin
		SetVariable set_gamma,noedit=1,frame=0,win=$subWin
		b=a ; alpha=90 ; bet=90 ; gam=90
		titleStr += "Tetragonal"
	elseif (SG>=16)														// Orthorhombic, a, b, c
		SetVariable set_a_nm,noedit=0,frame=1,win=$subWin	// enable
		SetVariable set_b_nm,noedit=0,frame=1,win=$subWin
		SetVariable set_c_nm,noedit=0,frame=1,win=$subWin
		SetVariable set_alpha,noedit=1,frame=0,win=$subWin	// disable
		SetVariable set_beta,noedit=1,frame=0,win=$subWin
		SetVariable set_gamma,noedit=1,frame=0,win=$subWin
		alpha=90 ; bet=90 ; gam=90
		titleStr += "Orthorhombic"
	elseif (SG>=3)															// Monoclinic, a, b, c, gamma
		SetVariable set_a_nm,noedit=0,frame=1,win=$subWin	// enable
		SetVariable set_b_nm,noedit=0,frame=1,win=$subWin
		SetVariable set_c_nm,noedit=0,frame=1,win=$subWin
		SetVariable set_beta,noedit=0,frame=1,win=$subWin
		SetVariable set_alpha,noedit=1,frame=0,win=$subWin	// disable
		SetVariable set_gamma,noedit=1,frame=0,win=$subWin
		alpha=90 ; gam=90										//	alpha=90 ; bet=90
		titleStr += "Monoclinic"
	else															// Triclinic, a, b, c, alpha, beta, gamma
		SetVariable set_a_nm,noedit=0,frame=1,win=$subWin	// enable all
		SetVariable set_b_nm,noedit=0,frame=1,win=$subWin
		SetVariable set_c_nm,noedit=0,frame=1,win=$subWin
		SetVariable set_alpha,noedit=0,frame=1,win=$subWin
		SetVariable set_beta,noedit=0,frame=1,win=$subWin
		SetVariable set_gamma,noedit=0,frame=1,win=$subWin
		titleStr += "Triclinic"
	endif
//	titleStr += "   \\F'Courier'"+getHMboth(SpaceGroupID2num(SpaceGroupID))
	titleStr += "   \\F'"+BAR_FONT_ALWAYS+"'"+getHMboth(SpaceGroupID2num(SpaceGroupID))

	titleStr = minus2bar(titleStr, single=1)								// change all minuses to a bar over following character
	Variable/C sizeLeft = titleStrLength(titleStr)
	TitleBox structureTitle,pos={imag(sizeLeft),63},title=titleStr, fSize=real(sizeLeft), win=$subWin
	SetVariable T_LatticeVar disable=(numtype(T_C)>0),win=$subWin

	if (exists("Init_AtomViewLattice")==6)
		Button buttonAtomView,win=$subWin,disable=1			// hide the button
		PopupMenu popupAtomView,win=$subWin,disable=0	// show the popup
	else
		Button buttonAtomView,win=$subWin,disable=0			// show the button
		PopupMenu popupAtomView,win=$subWin,disable=1	// hide the popup
	endif
	if (exists("Init_PowderPatternLattice")==6)
		Button buttonPowderPattern,win=$subWin,disable=1	// hide the button
		PopupMenu popupPowderPattern,win=$subWin,disable=0	// show the popup
	else
		Button buttonPowderPattern,win=$subWin,disable=0	// show the button
		PopupMenu popupPowderPattern,win=$subWin,disable=1	// hide the popup
	endif
	return 0
End
//
Static Function/C titleStrLength(str)
	String str
	Make/FREE/I sizes={14,12,11,10,9}
	Variable i=0, isize=sizes[0]
	for (i=0; i<numpnts(sizes); i+=1)
		isize=sizes[i]
		MeasureStyledText/F=GenevaEquivFont/SIZE=(isize) str
		if (i==0 && V_width < 194)
			break
		elseif (V_width < 214)
			break
		endif
	endfor
	Variable left = isize>=14 ? 24 : 4
	return cmplx(isize,left)
End
//
Function LatticePanelParamProc(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	if (sva.eventCode != 2 && sva.eventCode!=1)		// NOT Enter key or mouse up
		return 0
	endif

	NVAR dirty=root:Packages:Lattices:PanelValues:dirty
	strswitch(sva.ctrlName)
		case "set_a_nm":
		case "set_b_nm":
		case "set_c_nm":
		case "set_alpha":
		case "set_beta":
		case "set_gamma":
		case "setDesc":
		case "T_C":
			dirty = 1								// a big change, xtal.sourceFile no longer valid
			break
	endswitch

	if (!dirty)
		NVAR h=root:Packages:Lattices:PanelValues:h, k=root:Packages:Lattices:PanelValues:k, l=root:Packages:Lattices:PanelValues:l
		NVAR dspace_nm=root:Packages:Lattices:PanelValues:dspace_nm
		NVAR Fr=root:Packages:Lattices:PanelValues:Fr, Fi=root:Packages:Lattices:PanelValues:Fi
		NVAR T_C=root:Packages:Lattices:PanelValues:T_C
		STRUCT crystalStructure xtal			// returns 0 if something set, 1 is nothing done
		FillCrystalStructDefault(xtal)		//fill the crystal structure with 'current' values
		dspace_nm = dSpacing(xtal,h,k,l)
//		Variable/C Fc = Fstruct(xtal,h,k,l,T_K=(xtal.Temperature)+273.15)
		Variable/C Fc = Fstruct(xtal,h,k,l,T_K=T_C+273.15)
		Fr = real(Fc)
		Fi = imag(Fc)
	endif

	// set buttons to enable/disable, based upon 'dirty'
	String subWin=sva.win
	if (dirty)
		Button buttonLatticeSave disable=0,fColor=(57346,65535,49151),title="Save  (Needed)",win=$subWin
		Button buttonLatticeRevert disable=0,fColor=(65535,65534,49151),win=$subWin
	else
		Button buttonLatticeSave disable=2,fColor=(0,0,0),title="Saved",win=$subWin
		Button buttonLatticeRevert disable=2,fColor=(0,0,0),win=$subWin
	endif
	return 0
End
//
Static Function/T SelectNewSG(find)
	String find
	if (strlen(find)<1)
		find = StrVarOrDefault("root:Packages:Lattices:PanelValues:SpaceGroupSearch","")
		Prompt find, "Space Group Search, use * for wild card"
		DoPrompt "Search String", find
		if (V_flag)
			return ""
		endif
	endif
	String/G root:Packages:Lattices:PanelValues:SpaceGroupSearch=find

	String str,list0="", symList=""
	Variable i
	for (i=0;i<ItemsInList(CommonPredefinedStructures);i+=1)
		str = StringFromList(i,CommonPredefinedStructures)
		list0 += StringFromList(0,str,":")+";"
	endfor
	list0 = ListMatch(list0, find)

	for (i=0;i<ItemsinList(list0);i+=1)
		str = StringFromList(i,list0)
		sprintf str, "%d  %s;",NumberByKey(str,CommonPredefinedStructures),str
		symList += str
	endfor	

	String list = symmtry2SG(find,types=-1,printIt=0), sym
	list = RemoveDuplicatesFromList(list)
	String system, systemNames="Triclinic\t;Monoclinic\t;Orthorhombic;Tetragonal\t;Trigonal\t;Hexagonal\t;Cubic\t\t"
	String id, allIDs=MakeAllIDs()
	Variable Nlist=ItemsInList(list), idNum
	for (i=0; i<Nlist; i+=1)
		idNum = str2num(StringFromList(i,list))
		if (isValidSpaceGroupIDnum(idNum))
			id = StringFromList(idNum-1,allIDs)
			system = StringFromList(latticeSystem(id),systemNames)
			sprintf str, "%s  %s  %s  [%s];", id,system,getHMsym2(idNum),getHallSymbol(idNum)
			symList += str
		endif
	endfor
	Variable N=ItemsInList(symList)
	if (strsearch(find,"*",0)<0 && N<1)
		symList = SelectNewSG(find+"*")
		N = ItemsInList(symList)
	endif
	if (N<1)
		symmtry2SG("")							// this will print stuff to history, hopefully helpful
		return ""
	elseif (N==1)
		sym = StringFromList(0,symList)
		return StringFromList(0,sym," ")
	endif
	Prompt sym,"Space Group",popup,addDefaults2symList(symList)
	DoPrompt "Space Group",sym
	if (V_flag)
		return ""
	endif
	return StringFromList(0,sym," ")
End
//
Static Function/S addDefaults2symList(in)	// add "Default" to the default Space Groups
	String in
	String line, id, out=""
	Variable SG, i,N=ItemsInList(in)
	for (i=0;i<N;i+=1)
		line = StringFromList(i,in)
		id = StringFromList(0,line," ")
		SG = str2num(id)
		if (strsearch(id,":",0)<0)
			out += line+";"					// a simple SG with no options
		elseif (StringMatch(id,FindDefaultIDforSG(SG)))
			out += line+"  <--Default;"		// a simple SG with no options
		else
			out += line+";"					// a simple SG with no options
		endif
	endfor
	return out
End


Function LatticePanelButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	if (ba.eventCode != 2)			// not mouse up
		return 0
	endif

	String ctrlName=ba.ctrlName
	SVAR desc=root:Packages:Lattices:PanelValues:desc
	SVAR SpaceGroupID=root:Packages:Lattices:PanelValues:SpaceGroupID
	NVAR a=root:Packages:Lattices:PanelValues:a
	NVAR b=root:Packages:Lattices:PanelValues:b
	NVAR c=root:Packages:Lattices:PanelValues:c
	NVAR alpha=root:Packages:Lattices:PanelValues:alpha
	NVAR bet=root:Packages:Lattices:PanelValues:bet
	NVAR gam=root:Packages:Lattices:PanelValues:gam
	NVAR alphaT=root:Packages:Lattices:PanelValues:alphaT
	NVAR T_C=root:Packages:Lattices:PanelValues:T_C
	NVAR dirty = root:Packages:Lattices:PanelValues:dirty
	SVAR crystalStructStr = root:Packages:Lattices:PanelValues:crystalStructStr
	STRUCT crystalStructure xtal
	StructGet_xtal(crystalStructStr,xtal)					// pre-load with current information
	CleanOutCrystalStructure(xtal)
	STRUCT WMSetVariableAction sva
	sva.win = ba.win
	sva.eventCode = 2

	if (stringmatch(ctrlName,"buttonLatticeRevert"))	// close window, do not save numbers
		FillCrystalStructDefault(xtal)					//fill the lattice structure with default values
		a=xtal.a  ;  b=xtal.b  ;  c=xtal.c
		alpha=xtal.alpha  ;  bet=xtal.beta  ;  gam=xtal.gam
		SpaceGroupID = xtal.SpaceGroupID
		alphaT = xtal.alphaT
		desc = xtal.desc
		T_C = xtal.Temperature
		dirty = 0
		UpdatePanelLatticeConstControls(ba.win,SpaceGroupID)
		LatticePanelParamProc(sva)
		StructPut/S xtal, crystalStructStr
	elseif (stringmatch(ctrlName,"buttonLatticeSave") && dirty)
		xtal.a = a  ;  xtal.b = b  ;  xtal.c = c
		xtal.alpha = alpha  ;  xtal.beta = bet  ;  xtal.gam = gam
		xtal.SpaceGroupID = SpaceGroupID
		xtal.alphaT = alphaT
		xtal.Temperature = T_C
		xtal.desc = desc[0,99]								// desc is limited to 100 chars
		if (dirty==1)											// dirty==1, means xtal.sourceFile is bad
			xtal.sourceFile = ""
		endif
		ForceLatticeToStructure(xtal)
		UpdateCrystalStructureDefaults(xtal)
		dirty = 0
		UpdatePanelLatticeConstControls(ba.win,SpaceGroupID)// change main copy of xtal
		LatticePanelParamProc(sva)
		StructPut/S xtal, crystalStructStr				// update the local copy too
		OverOccupyList(xtal,printIt=1)					// print notice if some sites have occ>1
	elseif (stringmatch(ctrlName,"buttonFindSpaceGroup"))
		String id = SelectNewSG("")
		if (strlen(id)<1)
			return 0
		endif
		SpaceGroupID = id
		dirty = 1												// source file no longer valid
		UpdatePanelLatticeConstControls(ba.win,SpaceGroupID)
		LatticePanelParamProc(sva)
	elseif (stringmatch(ctrlName,"buttonEditAtomPositions"))
		// there is no longer a "buttonEditAtomPositions", we get here from LatticeEditPopMenuProc()
		if (EditAtomPositions(xtal)>=0)					// xtal has been MODIFIED
			if (ForceLatticeToStructure(xtal))			// sets everything including remaking atom0,atom1,...
				return NaN
			endif
			xtal.hashID = xtalHashID(xtal)
			StructPut/S xtal, crystalStructStr			// update the local copy
			dirty = 2											// source file still mostly OK
			UpdatePanelLatticeConstControls(ba.win,SpaceGroupID)
			LatticePanelParamProc(sva)
		endif
	elseif (stringmatch(ctrlName,"buttonLatticeFromFile"))
		if (readCrystalStructure(xtal,"",printIt=1))
			DoAlert 0,"nothing read"
		endif
		a = xtal.a  ;  b = xtal.b  ;  c = xtal.c
		alpha = xtal.alpha  ;  bet = xtal.beta  ;  gam = xtal.gam
		SpaceGroupID = xtal.SpaceGroupID
		desc = xtal.desc
		alphaT = xtal.alphaT
		T_C = xtal.Temperature
		T_C = (xtal.haveDebyeT && numtype(T_C)) ? 22.5 : T_C
		dirty = 2														// yes, it is dirty, but xtal.sourceFile is valid
		StructPut/S xtal, crystalStructStr
		UpdatePanelLatticeConstControls(ba.win,SpaceGroupID)
		LatticePanelParamProc(sva)
	elseif (stringmatch(ctrlName,"buttonPrintLattice"))	// print shown lattice parameters to the history
		xtal.a = a  ;  xtal.b = b  ;  xtal.c = c
		xtal.alpha = alpha  ;  xtal.beta = bet  ;  xtal.gam = gam
		xtal.SpaceGroupID = SpaceGroupID
		xtal.desc = desc[0,99]
		xtal.alphaT = alphaT
		xtal.Temperature = T_C
		ForceLatticeToStructure(xtal)
		print " "
		print_crystalStructure(xtal)
	elseif (stringmatch(ctrlName,"buttonWriteLattice"))	// write current lattice parameters to an xml file
		writeCrystalStructure2xmlFile("","")
	elseif (stringmatch(ctrlName,"buttonAtomView"))	// add support for AtomView
		String cmd
		sprintf cmd,"LatticeSym#UpdatePanelLatticeConstControls(\"%s\",\"%s\")",ba.win,SpaceGroupID
		Execute/P "INSERTINCLUDE \"AtomView\", version>=0.43"
		Execute/P "COMPILEPROCEDURES "
		Execute/P "Init_AtomViewLattice()"
		Execute/P cmd
	elseif (stringmatch(ctrlName,"buttonPowderPattern"))	// add support for PowderPatterns
		sprintf cmd,"LatticeSym#UpdatePanelLatticeConstControls(\"%s\",\"%s\")",ba.win,SpaceGroupID
		Execute/P "INSERTINCLUDE \"PowderPatterns\", version>=0.24"
		Execute/P "COMPILEPROCEDURES "
		Execute/P "Init_PowderPatternLattice()"
		Execute/P cmd
	endif
	return 0
End


Static Function LatticeEditPopMenuProc(pa) : PopupMenuControl		// used in the LatticeSym Panel
	STRUCT WMPopupAction &pa
	if (pa.eventCode != 2)
		return 0
	endif

	if (!StringMatch(pa.ctrlName,"popupLatticeEdit"))
		return 0
	elseif (strsearch(pa.popStr,"Edit Atom Positions",0,2)>=0)
		STRUCT WMButtonAction ba
		ba.eventCode = 2
		ba.win = pa.win
		ba.ctrlName = "buttonEditAtomPositions"
		LatticePanelButtonProc(ba)
	elseif (strsearch(pa.popStr,"Change Current Setting",0,2)>=0)
		ChangeSettingCurrentXtal("",printIt=1)
	elseif (strsearch(pa.popStr,"Find Closest hkl",0,2)>=0)
		findClosestHKL(NaN)									// help={"Knowing either the d-spacing or the Q, find closest hkl's"}
	elseif (strsearch(pa.popStr,"Space Group number <--> symmetry",0,2)>=0)
		symmtry2SG("")											// help={"find the Space Group number from symmetry string,  e.g. Pmma, or sym from number"}
	elseif (strsearch(pa.popStr,"Describe the Symmetry Operations",0,2)>=0)
		DescribeSymOps($"", printIt=1)
	elseif (strsearch(pa.popStr,"angle between two hkl's",0,2)>=0)
		angleBetweenHKLs(NaN,NaN,NaN,  NaN,NaN,NaN, printIt=1)
	endif
	return 0
End


Static Function/T LatticeSetPanelList()			// returns list of win#subWin for all of the "LatticePanel"'s
	String panelList=WinList("*",";","WIN:64")	// list of ALL Panels
	String win, list=""
	Variable i
	for (i=0;i<ItemsInList(panelList);i+=1)
		win = GetUserData(StringFromList(i,panelList), "", "LatticePanelName")
		if (strlen(win))
			list += win+";"
		endif
	endfor
	return list
End

//	End of lattice set panel
// =========================================================================
// =========================================================================


// =========================================================================
// =========================================================================
//	Start of reading/writing lattice paremeters to a file

Function LoadCrystal(fname)
	String fname
	STRUCT crystalStructure xtal					// this sruct is filled  by this routine
	return readCrystalStructure(xtal,fname,printIt=1)
End


//Function testReadCrystal()
//	STRUCT crystalStructure xtal					// this sruct is filled  by this routine
//	readCrystalStructure(xtal,"")
//	print_crystalStructure(xtal)
//End
Static Function readCrystalStructure_xtl(xtal,fname)
	STRUCT crystalStructure &xtal					// this sruct is filled  by this routine
	String fname

	String list = keyStrFromFile(fname,"CrystalStructure","materialsPath")
	if (strlen(list)<1)								// try old style
		return read_cri_fileOLD(xtal,fname)
	endif

	Variable a,b,c,alpha,bet,gam, SpaceGroup, alphaT
	sscanf StringByKey("latticeParameters",list,"="), "{ %g, %g, %g, %g, %g, %g }"  , a,b,c,alpha,bet,gam
	if (V_flag!=6)
		return 1
	endif
	SpaceGroup = str2num(StringByKey("SpaceGroup",list,"="))
	if (numtype(SpaceGroup))
		return 1
	endif
	String unit = StringByKey("lengthUnit",list,"=")
	if (stringmatch(unit,"Ang*") || stringmatch(unit, ARING+"*"))
		a*= 10  ;  b*= 10  ;  c*= 10
	endif
	alphaT = str2num(StringByKey("latticeAlphaT",list,"="))

	String fullFile = StringByKey("keyStrFileName",list,"=")
	Variable i0,i1
	i1 = strlen(fullFile)-1						// possibly trim file length to fit, keep last part of fullFile
	i0 = max(0,i1-MAX_FILE_LEN)
	fullFile = fullFile[i0,i1]
	xtal.sourceFile = fullFile

	xtal.hashID = ""
	xtal.a = a  ;  xtal.b = b  ;  xtal.c = c
	xtal.alpha = alpha  ;  xtal.beta = bet  ;  xtal.gam = gam
	xtal.SpaceGroup = SpaceGroup
	String id = FindDefaultIDforSG(xtal.SpaceGroup)
	xtal.SpaceGroupID = id
	xtal.SpaceGroupIDnum = SpaceGroupID2num(id)		// change id to id number in [1,530]
	xtal.alphaT = !(alphaT>0) ? 0 : alphaT
	ForceLatticeToStructure(xtal)

	String str = StringByKey("structureDesc",list,"=")
	if (strlen(str)==0)
		str = StringByKey("latticeDesc",list,"=")
	endif
	xtal.desc = str[0,99]
	Variable N=0
	Variable i,xx,yy,zz, occ,DebyeT=0
	String item
	do
		item = StringByKey("AtomDesctiption"+num2istr(N+1),list,"=")
		sscanf item, "{%s %g %g %g %g}",str,xx,yy,zz,occ
		if (V_flag!=5)
			sscanf item, "{%s %g %g %g}",str,xx,yy,zz
			if (V_flag!=4)
				break
			endif
			occ = 1
		endif
		xtal.atom[N].name = str[0,59]
		xtal.atom[N].Zatom = ZfromLabel(str)
		xtal.atom[N].x = xx
		xtal.atom[N].y = yy
		xtal.atom[N].z = zz
		xtal.atom[N].occ = occ
		xtal.atom[N].DebyeT = DebyeT
		N += 1
	while(strlen(str) && N<STRUCTURE_ATOMS_MAX-1)
	xtal.N = N

	xtal.Unconventional00=NaN;  xtal.Unconventional01=NaN;  xtal.Unconventional02=NaN	// transform matrix for an unconventional unit cel
	xtal.Unconventional10=NaN;  xtal.Unconventional11=NaN;  xtal.Unconventional12=NaN	// default to a conventional cell
	xtal.Unconventional20=NaN;  xtal.Unconventional21=NaN;  xtal.Unconventional22=NaN
	str = StringByKey("Unconventional",list,"=")
	if (strlen(str))										// found an $Unconventional tag
		Variable u00,u10,u20, u01,u11,u21, u02,u12,u22
		sscanf str,"{ {%g,%g,%g}, {%g,%g,%g}, {%g,%g,%g} }",u00,u10,u20, u01,u11,u21, u02,u12,u22
		if (V_flag==9 && numtype(u00+u10+u20+u01+u11+u21+u02+u12+u22)==0)
			xtal.Unconventional00 = u00;	xtal.Unconventional01 = u01;	xtal.Unconventional02 = u02
			xtal.Unconventional10 = u10;	xtal.Unconventional11 = u11;	xtal.Unconventional12 = u12
			xtal.Unconventional20 = u20;	xtal.Unconventional21 = u21;	xtal.Unconventional22 = u22
 		endif
	endif

	if (N<1)
		printf "Failed to Load Crystal Structure from '%s'\r",StringByKey("keyStrFileName",list,"=")
		return 1
	else
		CleanOutCrystalStructure(xtal)
		UpdateCrystalStructureDefaults(xtal)
		printf "Loaded Crystal Structure from '%s'     using units of '%s'\r",StringByKey("keyStrFileName",list,"="),unit
	endif
	return 0
End
//
Function read_cri_fileOLD(xtal,fname)
	STRUCT crystalStructure &xtal					// this sruct is filled  by this routine
	String fname

	Variable f
	Open/R/M="the .cri file"/T=".cri" f as fname
	if (f==0)
		return 1
	endif

	Variable SpaceGroup, a,b,c,alpha,bet,gam
	String line
	xtal.hashID = ""

	FReadLine f, line
	line = ReplaceString("\r",line,"")
	xtal.desc = line[0,99]
	FReadLine f, line
	SpaceGroup = str2num(line)
	if (!isValidSpaceGroup(SpaceGroup))
		Close f
		return 1
	endif
	xtal.SpaceGroup = SpaceGroup

	FReadLine f, line
	sscanf line, "%g %g %g  %g %g %g",a,b,c,alpha,bet,gam
	if (V_flag!=6)
		Close f
		return 1
	endif
	xtal.a = a
	xtal.b = b
	xtal.c = c
	xtal.alpha = alpha
	xtal.beta = bet
	xtal.gam = gam

	FReadLine f, line
	Variable N = str2num(line)
	if (!(N>=1 && N<=STRUCTURE_ATOMS_MAX))
		Close f
		return 1
	endif
	xtal.N = N
	Variable i,xx,yy,zz, occ
	String str
	for (i=0;i<N;i+=1)
		FReadLine f, line

		sscanf line, "%s %g %g %g %g",str,xx,yy,zz,occ
		if (V_flag!=5)
			sscanf line, "%s %g %g %g",str,xx,yy,zz
			if (V_flag!=4)
				Close f
				return 1
			endif
			occ = 1
		endif
		xtal.atom[i].name = str[0,59]
		xtal.atom[i].Zatom = ZfromLabel(str)
		xtal.atom[i].x = xx
		xtal.atom[i].y = yy
		xtal.atom[i].z = zz
		xtal.atom[i].occ = occ
	endfor
	Close f
	if (ItemsInList(GetRTStackInfo(0))<2)
		printf "read structure from file  '%s'\r",S_fileName
	endif
	xtal.Unconventional00=NaN;  xtal.Unconventional01=NaN;  xtal.Unconventional02=NaN	// transform matrix for an unconventional unit cel
	xtal.Unconventional10=NaN;  xtal.Unconventional11=NaN;  xtal.Unconventional12=NaN	// default to a conventional cell
	xtal.Unconventional20=NaN;  xtal.Unconventional21=NaN;  xtal.Unconventional22=NaN
	ForceLatticeToStructure(xtal)					// set recip and direct and makes lattice constants valid
	UpdateCrystalStructureDefaults(xtal)
	return 0
End
//
ThreadSafe Static Function ZfromLabel(symb)	// returns Z for an atomic symbol (NOT case sensitive), used for the atomic strucure factor
	String symb					// atomic symbol
	symb = symb[0,1]
	symb[0,0] = UpperStr(symb[0,0])		// ensure first char is upper case
	symb[1,1] = LowerStr(symb[1,1])		// and second character is lower
	Variable c=char2num(symb[1,1])
	if (c<97 || c>122)
		symb = symb[0,0]
	endif
	Variable iz = WhichListItem(symb,ELEMENT_Symbols)+1
	return ((iz>0) ? iz : NaN)
End
//
ThreadSafe Static Function/T Z2symbol(Z)	// returns chemical symbol from atomic number Z.
	Variable Z					// atomic number
	if (!(Z>=1 && Z<=ELEMENT_Zmax))
		return ""
	endif
	return StringFromLIst(Z-1,ELEMENT_Symbols)
End


Static Function/T MinimalChemFormula(xtal,[maximal])
	STRUCT crystalStructure &xtal
	Variable maximal
	maximal = ParamIsDefault(maximal) ? 0 : maximal
	maximal = numtype(maximal) ? 0 : !(!maximal)
	Variable N = xtal.N
	if (N < 1)
		return ""
	endif

	Make/N=(N)/T/FREE symbs="", symbs0=""
	Make/N=(N)/D/FREE nums=NaN, Zs=NaN, nums0=NaN

	reMakeAtomXYZs(xtal)
	String symb
	Variable i, occ, Zatom
	for (i=0;i<N;i+=1)						// for each atom type
		occ = xtal.atom[i].occ
		occ = occ>0 ? occ : 1
		occ *= DimSize($("root:Packages:Lattices:atom"+num2istr(i)),0)
		Zatom = xtal.atom[i].Zatom
		symb = Z2symbol(Zatom)
		if (!(occ > 0) || !(Zatom>0) || strlen(symb)<1)
			continue
		endif
		symbs[i] = symb
		nums[i] = occ
		Zs[i] = Zatom
	endfor

	Variable j, m=0
	for (i=0;i<N;i+=1)						// combine all occurances of each element
		occ = 0
		Zatom = Zs[i]
		if (!(Zatom>0))
			continue
		endif
		Zs[i] = NaN
		occ = nums[i]
		for (j=i+1;j<N;j+=1)
			if (Zs[j]==Zatom)
				occ += nums[j]
				Zs[j] = NaN
			endif
		endfor
		if (occ>0)
			symbs0[m] = symbs[i]
			nums0[m] = occ
			m += 1
		endif
	endfor
	N = m
	Redimension/N=(N) symbs0, nums0

	if (!maximal)						// if !maximal, then remove common factors
		WaveStats/Q/M=1 nums0
		Variable small=V_min
		if (V_min==V_max)				// all nums0 are the same, so use 1
			nums0 = 1
		elseif (small >= 2)			// try to remove common factors
			for (j=small;j>=2;j-=1)
				MatrixOP/FREE/O remain = nums0 / j
				remain = mod(remain,1)
				if (sum(remain)==0)
					nums0 /= j
				endif
			endfor
		endif
	endif

	String formula="", str
	for (m=0;m<N;m+=1)			// create the formula string
		if (nums0[m]==1)
			str = symbs0[m]
		else
			sprintf str, "%s%g",symbs0[m],nums0[m]
		endif
		formula += str
	endfor
	return formula
End
//
//Function testRead()
//	STRUCT crystalStructure xtal					// this sruct is filled  by this routine
//	read_cri_fileOLD(xtal,"")
//	print_crystalStructure(xtal)
//End
//	Al2O3 crystal (hexagonal axis)
//	   167
//	   0.47580   0.47580   1.2991  90.00000  90.00000 120.00000
//	     2
//	Al001    0.00000   0.00000   0.35230   1.00000
//	O0001    0.30640   0.00000   0.25000   1.00000
//
//Function ReadLatticeFromKeyFile(fileName,path,xtal)	// read in the geometry structure from a keyword file
//	String fileName										// full path name to the file
//	String path											// name of an Igor path to use
//	STRUCT crystalStructure &xtal						// the structure to fill from the file
//
//	Variable tempPath = 0
//	if (strlen(path)<1)
//		path = "userMaterials"
//		PathInfo userMaterials
//		if (V_flag)
//		else
//			PathInfo Igor
//			NewPath/Z/Q userMaterials,S_path+"User Procedures:materials:"
//			tempPath = !V_flag
//			if (V_flag)
//				path = "home"							// no materials folder, just look in home
//			endif
//		endif
//	endif
//	String list = keyStrFromFile(fileName,"LatticeFile",path)// read in all of the tagged values into a keyword list
//	if (strlen(list)<1)
//		return 1
//	endif
//	if (tempPath)
//		KillPath/Z $path
//	endif
//
//	Variable a,b,c,alpha,beta,gam, SpaceGroup, alphaT
//	sscanf StringByKey("latticeParameters",list,"="), "{ %g, %g, %g, %g, %g, %g }"  , a,b,c,alpha,beta,gam
//	if (V_flag!=6)
//		return 1
//	endif
//	SpaceGroup = str2num(StringByKey("SpaceGroup",list,"="))
//	if (numtype(SpaceGroup))
//		return 1
//	endif
//	String unit = StringByKey("lengthUnit",list,"=")
//	if (stringmatch(unit,"Ang*"))
//		a*= 10  ;  b*= 10  ;  c*= 10
//	endif
//	alphaT = str2num(StringByKey("latticeAlphaT",list,"="))
//
//	xtal.a = a  ;  xtal.b = b  ;  xtal.c = c
//	xtal.alpha = alpha  ;  xtal.beta = beta  ;  xtal.gam = gam
//	xtal.SpaceGroup = SpaceGroup
//	xtal.desc = StringByKey("latticeDesc",list,"=")
////	xtal.alphaT = !(alphaT>0) ? 0 : alphaT
//	xtal.alphaT = !(abs(alphaT)<0.1) ? 0 : alphaT	// alphaT can be negative, but it must be small
//	ForceLatticeToStructure(xtal)
//	printf "Loaded lattice information from '%s'     using units of '%s'\r",StringByKey("keyStrFileName",list,"="),unit
//	return 0
//End
Static Function/S keyStrFromFile(fname,ftype,path)	// read in all of the tagged values from a file into a keyword list
	String fname									// full path name to file with tagged geometry values
	String ftype									// the required file identifier, included as a tag (ftype is optional)
	String path									// name of Igor path

	Variable refNum
	String fullName=""

	PathInfo materialsPath
	path = SelectString(V_flag, "", path)		// set path="" if it does not exist

	GetFileFolderInfo/P=$path/Q/Z fname
	fullName = SelectString(V_flag,S_Path,"")
	if (!V_isFile)
		Open/D/M="file with a lattice keyword list"/P=$path/R/T=".xtl" refNum as fname
		fullName = SelectString(V_flag,S_fileName,"")
	endif
	if (strlen(fullName)<1)
		return ""
	else
		Open/P=$path/R/Z=1 refNum as fullName
		fullName = SelectString(V_flag,S_fileName,"")
	endif
	if (strlen(fullName)<1 || !refNum)
		return ""
	endif
// the above 12 lines should not be necessary, the /Z=2 switch does not work right

	String tagStr,valStr,line
	if (strlen(ftype))							// an ftype is present, so check the file
		FReadLine refNum, line
		extractTagValueFromLine(line,tagStr,valStr)	// check both ways
		if (!stringmatch(tagStr,ftype) && !(stringmatch(tagStr,"filetype") && stringmatch(valStr,ftype)))
			printf "The file must start with '%s', but the first line is '%s'\r",ftype,line
			Close refNum
			return ""
		endif
	endif

	String list = ReplaceStringByKey("keyStrFileName","",fullName,"=")
	Variable i,dollar = char2num("$")
	do
		FReadLine refNum, line
		extractTagValueFromLine(line,tagStr,valStr)
		if (strlen(tagStr))
			// check if tag is already in list, only take the first instance of a tag, not the last
			if (keyInList(tagStr,list,"=",""))	// check if this key is already in list
				continue
			endif
			list = ReplaceStringByKey(tagStr,list,valStr,"=")
		endif
	while (strlen(line))							// go until EOF
	Close refNum
	return list
End
//
Static Function extractTagValueFromLine(line,tagStr,valStr)
	String line										// input line, presumably starting with a %tag
	String &tagStr								// returned
	String &valStr								// returned
	tagStr = ""									// init both to empty
	valStr = ""

	if (char2num(line)!=char2num("$"))		// if it does not starts with a $, do not process
		return 0
	endif

	Variable i = strsearch(line,"//",0)			// strip off comments
	if (i>=0)
		line = line[0,i-1]
	endif
	for (i=0;char2num(line[i+1])>32;i+=1)	// find end of tag, it ends with a space or lower
	endfor
	tagStr = ReplaceString(";",line[1,i],":")
	tagStr = ReplaceString("=",tagStr,"_")
	if (strlen(tagStr)<1)
		return 1
	endif

	for (i=i+1;char2num(line[i])<=32;i+=1)	// find first non-white space
	endfor
	valStr = line[i,Inf]							// value associated with tagStr

	for (i=strlen(valStr)-1;char2num(valStr[i])<=32 && i>0;i-=1)	// strip off trailing whitespace
	endfor
	valStr = ReplaceString(";",valStr[0,i],":")
	valStr = ReplaceString("=",valStr,"_")
	return 1
End

//	End of of reading/writing lattice paremeters to a file
// =========================================================================
// =========================================================================



//  ========================================================================= //
//  ================== Start of reading/writing xml files ==================  //

Function writeCrystalStructure2xmlFile(path,fname)	// writes the current xtal to an xml file
	String path									// path to write to
	String fname								// name of file to write

	STRUCT crystalStructure xtal				// this sructure is written in this routine
	FillCrystalStructDefault(xtal)
	ForceLatticeToStructure(xtal)
	String cif = crystalStructure2xml(xtal,NEW_LINE)	// convert xtal to cif

	Variable f
	Open/C="R*ch"/F="XML Files (*.xml):.xml;"/M="file to write"/P=$path/Z=2 f as fname
	if (V_flag==0)
		fprintf f,  "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>"+NEW_LINE+NEW_LINE
		FBinWrite f, cif
		Close f
	endif
	return V_flag
End


Function readCrystalStructure(xtal,fname,[printIt])
	STRUCT crystalStructure &xtal					// this sruct is filled  by this routine
	String fname
	Variable printIt
	printIt = ParamIsDefault(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt

	fname = FindMaterialsFile(fname)				// find full path to fname, and optionally set materials path

	PathInfo materialsPath
	if (V_flag==0)											// make it if it does not exist
		NewPath/Z materialsPath, ParseFilePath(1,FunctionPath("CrystalsAreHere"),":",1,0)
	endif
	PathInfo materialsPath
	if (V_flag==0)											// make it if it does not exist
		NewPath/Z materialsPath, SpecialDirPath("Documents",0,0,0)+"materials"
		if (V_flag)											// path still does not exist, use Documents folder
			NewPath/Z materialsPath, SpecialDirPath("Documents",0,0,0)
		endif
	endif

	// find the file type, this first looks at the top of the file, if that does not work, it looks at the file name extension
	String fileType = CrystalFileType(fname)	// "xml", "xtl", "cif", "", or the actual file extension (without the '.')
	Variable err=1
	if (stringmatch(fileType,"xml"))
		err = readCrystalStructureXML(xtal,fname)
	elseif (stringmatch(fileType,"xtl"))
		err = readCrystalStructure_xtl(xtal,fname)
	elseif (stringmatch(fileType,"cif"))
		err = readCrystalStructureCIF(xtal,fname)
	endif
	if (printIt)
		printf "\r%s xtal structure from '%s'\r\r",SelectString(err,"Loaded","ERROR -- Loading"),fname
		if (!err)
			OverOccupyList(xtal,printIt=1)						// print notice if some sites have occ>1
		endif
	endif
	return err
End
//
Static Function/S FindMaterialsFile(fname)		// returns full path to a materials file
	String fname
	String fileFilters = "XML Files (*.xml):.xml;XTL Files (*.xtl):.xtl;CIF Files (*.cif):.cif;old cri Files (*.cri):.cri;All Files:.*;"

	GetFileFolderInfo/Q/Z fname						// full path name passed, that was east
	if (V_isfile)
		return fname
	endif

	String dirString, dirList=""						// a list of possible materials directories
	PathInfo materialsPath
	String materialsDir = SelectString(V_flag,"__empty__",S_path)			// materialsPath is already set, the 1st choice

	dirString = SpecialDirPath("Igor Pro User Files",0,0,0)+"materials:"	// user's local copy, 2nd choice
	GetFileFolderInfo/Q/Z dirString
	String usersPath = SelectString(V_isFolder && !V_flag,"",dirString)

	dirString = SpecialDirPath("Documents",0,0,0)+"materials:"				// in the "Documents" folder, 3rd choice
	GetFileFolderInfo/Q/Z dirString
	String docPath = SelectString(V_isFolder && !V_flag,"",dirString)

	dirString=ParseFilePath(1,FunctionPath("CrystalsAreHere"),":",1,0)	// the standard distribution, 4th choice
	String stdPath = SelectString(strlen(dirString),"",dirString)

	dirString = SpecialDirPath("Igor Application",0,0,0)+"materials:"		// system local copy, 5th choice
	GetFileFolderInfo/Q/Z dirString
	String appPath = SelectString(V_isFolder && !V_flag,"",dirString)

	String dirKeyList=""
	if (StringMatch(materialsDir,usersPath))
		dirList = "User files;"
		dirKeyList = ReplaceStringByKey("User files",dirKeyList,materialsDir)
	elseif (StringMatch(materialsDir,docPath))
		dirList = "Documents;"
		dirKeyList = ReplaceStringByKey("Documents",dirKeyList,materialsDir)
	elseif (StringMatch(materialsDir,stdPath))
		dirList = "Standard Distribution;"
		dirKeyList = ReplaceStringByKey("Standard Distribution",dirKeyList,materialsDir)
	elseif (StringMatch(materialsDir,appPath))
		dirList = "Igor application folder;"
		dirKeyList = ReplaceStringByKey("Igor application folder",dirKeyList,materialsDir)
	else
		dirList = SelectString(StringMatch(materialsDir,"__empty__"), materialsDir+";","")
	endif

	if (WhichListItem("User files",dirList)<0 && strlen(usersPath))
		dirList += "User files;"
		dirKeyList = ReplaceStringByKey("User files",dirKeyList,usersPath)
	endif
	if (WhichListItem("Documents",dirList)<0 && strlen(docPath))
		dirList += "Documents;"
		dirKeyList = ReplaceStringByKey("Documents",dirKeyList,docPath)
	endif
	if (WhichListItem("Standard Distribution",dirList)<0 && strlen(stdPath))
		dirList += "Standard Distribution;"
		dirKeyList = ReplaceStringByKey("Standard Distribution",dirKeyList,stdPath)
	endif
	if (WhichListItem("Igor application folder",dirList)<0 && strlen(appPath))
		dirList += "Igor application folder;"
		dirKeyList = ReplaceStringByKey("Igor application folder",dirKeyList,appPath)
	endif

	fname = ParseFilePath(0,fname,":",1,0)			// check for fname in each of dirList
	if (strlen(fname))
		Variable i
		for (i=0;i<ItemsInList(dirList);i+=1)
			dirString = StringFromList(i,dirList)
			GetFileFolderInfo/Q/Z dirString+fname
			if (V_isFile)									// found fname in dirString
				UpdateMaterialsPath(dirString)
				return dirString+fname 
			endif
		endfor
	endif

	// could not find the materials file, time to start asking user for help
	dirString = StringFromList(0,dirList)
	if (itemsInList(dirList)>1)
		Prompt dirString,"Folder with materials",popup,dirList
		DoPrompt "Directory",dirString
		if (V_flag)
			return ""
		endif
	endif

	UpdateMaterialsPath(StringByKey(dirString,dirKeyList))

	PathInfo materialsPath
	String path=SelectString(V_flag, "", "materialsPath")
	Variable f
	Open/D=2/R/F=fileFilters/M="File with Crystal Information"/P=$path f as fname
	fname = S_fileName
	return fname
End
//
Static Function UpdateMaterialsPath(dirString)
	String dirString
	GetFileFolderInfo/Q/Z dirString
	if (!V_isFolder)
		return 1
	endif
	PathInfo materialsPath
	S_path = SelectString(V_flag, "", S_path)
	if (V_flag && stringmatch(S_path,dirString))		// nothing to change, all is OK
		return 0
	endif
	NewPath/Z/O materialsPath, dirString
	return 0
End
//
Static Function/T CrystalFileType(fname,[path])
	// Return the file type, of a crystal file.
	//	First look at top of the file, if that doesn't work, look at the file name extension, 
	//	Returned values are one of: "xml", "xtl", "cif", "", or the actual file extension (without the '.')
	String fname
	String path				// name of path, defaults to materialsPath
	path = SelectString(ParamIsDefault(path), path, "materialsPath")

	PathInfo $path
	path = SelectString(V_flag, "", path)		// set path="" if it does not exist

	String extension="", line=""
	Variable f
	Open/R/P=$path/R/Z=1 f as fname
	if (V_flag==0)
		FReadLine/N=100 f, line
		Close f
	endif

	// first try to determine the extension by looking at the first line of the file
	if (strsearch(line,"<?xml",0,2)>=0)
		extension = "xml"
	elseif (strsearch(line,"CrystalStructure",0,2)>=0)
		extension = "xtl"
	elseif (strsearch(line,"data_",0,2)>=0)
		extension = "cif"
	endif

	if (strlen(extension) < 1)	// need to actually look at the file name extension
		extension = ParseFilePath(4,fname,":",0,0)
	endif
	return extension
End

Static Function readCrystalStructureXML(xtal,fileName,[path])
	STRUCT crystalStructure &xtal
	String fileName
	String path
	if (ParamIsDefault(path))
		path = "home"
	endif

	Variable f
	String fileFilters = "XML Files (*.xml):.xml;TEXT Files (*.txt):.txt;All Files:.*;"
	Open/R/F=fileFilters/M="Select xml file with Crystal Description"/P=$path/R/Z=2 f as fileName
	if (strlen(S_fileName)<1 || V_flag)
		return 1
	endif
	String fullFile = S_fileName

	FStatus f
	String buf=PadString("",V_logEOF,0x20)
	FBinRead f, buf
	Close f
	buf = XMLremoveComments(buf)

	String cif = XMLtagContents("cif",buf)
	if (strlen(cif)<10)
		return 1
	endif

	Variable i0,i1
	i1 = strlen(fullFile)-1						// possibly trim file length to fit, keep last part of fullFile
	i0 = max(0,i1-MAX_FILE_LEN)
	fullFile = fullFile[i0,i1]
	xtal.sourceFile = fullFile

	xtal.hashID = ""
	String str = XMLtagContents("chemical_name_common",cif)
	xtal.desc = str[0,99]

	// 1st try for the id directly
	Variable SG=-1
	xtal.SpaceGroupID = TrimFrontBackWhiteSpace(XMLtagContents("space_group_id",cif))
	// 2nd try using the H-M symbol to set the SpaceGroupID & SpaceGroupIDnum
	if (strlen(xtal.SpaceGroupID)<1)
		String HMsym = XMLtagContents("H-M",cif)			//	Hermann-Manguin symbol, e.g. 'I 1 2/a 1'
		String nlist = SymString2SGtype(HMsym,2,0)			// checks in getHMsym2
		if (strlen(nlist)<1)
			nlist = SymString2SGtype(HMsym,2,1)				// checks in getHMsym2 again, but ignoring minus signs
		endif
		if (ItemsInList(nlist)==1)								// just one result, use it
			Variable idNum = str2num(StringFromList(0,nlist))
			if (isValidSpaceGroupIDnum(idNum))					// found valid SpaceGroupIDnum
				String allIDs=MakeAllIDs()
				xtal.SpaceGroupID = StringFromList(idNum-1, allIDs)
				printf "Setting Space Group from H-M = \"%s\"\r", HMsym
			endif
		endif
	endif
	// 3rd try, look at the symmetry operations
	if (strlen(xtal.SpaceGroupID)<1)
		String symLine = GetSymLinesFromXMLbuffer(cif)	// get the symmetry operations
		xtal.SpaceGroupID = FindIDfromSymOps(symLine)
	endif
	// 4th try, look for the Space Group number [1,230], and use the default id for that Space Group
	if (strlen(xtal.SpaceGroupID)<1)
		SG = str2num(XMLtagContents("space_group_IT_number",cif))
		if (SG>0)
			xtal.SpaceGroupID = FindDefaultIDforSG(SG)		// try to get id from SG [1,230]
		endif
	endif
	if (!isValidSpaceGroupID(xtal.SpaceGroupID))			// give up
		Abort "cannot find the Space Group from CIF file"
	endif
	xtal.SpaceGroupIDnum = SpaceGroupID2num(xtal.SpaceGroupID)		// change id to id number in [1,530]
	SG = str2num(xtal.SpaceGroupID)
	xtal.SpaceGroup = SG
	String SpaceGroupID = xtal.SpaceGroupID

	String cell = XMLtagContents("cell",cif)			// cell group
	String unit
	xtal.a = str2num(XMLtagContents("a",cell))
	unit = StringByKey("unit", XMLattibutes2KeyList("a",cell),"=")
	xtal.a *= ConvertUnits2meters(unit,defaultLen=1e-10)*1e9	// want length in nm
	xtal.b = str2num(XMLtagContents("b",cell))
	unit = StringByKey("unit", XMLattibutes2KeyList("b",cell),"=")
	xtal.b *= ConvertUnits2meters(unit,defaultLen=1e-10)*1e9
	xtal.c = str2num(XMLtagContents("c",cell))
	unit = StringByKey("unit", XMLattibutes2KeyList("c",cell),"=")
	xtal.c *= ConvertUnits2meters(unit,defaultLen=1e-10)*1e9
	xtal.alpha = str2num(XMLtagContents("alpha",cell))
	xtal.beta = str2num(XMLtagContents("beta",cell))
	xtal.gam = str2num(XMLtagContents("gamma",cell))
	xtal.alphaT = str2num(XMLtagContents("alphaT",cell))
	Variable Temperature = str2num(XMLtagContents("temperature",cell))
	unit = StringByKey("unit", XMLattibutes2KeyList("temperature",cell),"=")
	unit = SelectString(strlen(unit),"C",unit)						// default Temperature units are C
	xtal.Temperature = ConvertTemperatureUnits(Temperature,unit,"C")
	xtal.Unconventional00=NaN;  xtal.Unconventional01=NaN;  xtal.Unconventional02=NaN	// transform matrix for an unconventional unit cel
	xtal.Unconventional10=NaN;  xtal.Unconventional11=NaN;  xtal.Unconventional12=NaN	// default to a conventional cell
	xtal.Unconventional20=NaN;  xtal.Unconventional21=NaN;  xtal.Unconventional22=NaN
	String uncoStr = XMLtagContents("Unconventional",cell)
	if (strlen(uncoStr))								// found an $Unconventional tag
		Variable u00,u10,u20, u01,u11,u21, u02,u12,u22
		sscanf uncoStr,"{ {%g,%g,%g}, {%g,%g,%g}, {%g,%g,%g} }",u00,u10,u20, u01,u11,u21, u02,u12,u22
		if (V_flag==9 && numtype(u00+u10+u20+u01+u11+u21+u02+u12+u22)==0)
			xtal.Unconventional00 = u00;	xtal.Unconventional01 = u01;	xtal.Unconventional02 = u02
			xtal.Unconventional10 = u10;	xtal.Unconventional11 = u11;	xtal.Unconventional12 = u12
			xtal.Unconventional20 = u20;	xtal.Unconventional21 = u21;	xtal.Unconventional22 = u22
 		endif
	endif

	// collect the atom sites values
	String atomSite, atomLabel, symbol, WyckoffSymbol
	Variable Zatom, fracX,fracY,fracZ,occupy,valence,DebyeT
	Variable Biso, Uiso, aU11,aU22,aU33, aU12,aU13,aU23, mult
	Variable start=0, N=0
	do
		atomSite = XMLtagContents("atom_site",cif, start=start)	// next atom site
		if (strlen(atomSite)<10)										// could not find another atom_site
			break
		endif
		atomLabel = XMLtagContents("label",atomSite)			// label for this atom type
		symbol = XMLtagContents("symbol",atomSite)				// atomic symbol for this atom type
		atomLabel = SelectString(strlen(atomLabel), symbol, atomLabel)	// if no label given, use the symbol (may not be unique)
		Zatom = ZfromLabel(SelectString(strlen(symbol),atomLabel,symbol))	// get Z
		WyckoffSymbol = XMLtagContents("WyckoffSymbol",atomSite)
		Wave xyz = str2vec(XMLtagContents2List("fract_xyz",atomSite))
		if (WaveExists(xyz) && numpnts(xyz)==3)
			fracX = xyz[0]
			fracY = xyz[1]
			fracZ = xyz[2]
		else
			fracX = str2num(XMLtagContents("fract_x",atomSite))
			fracY = str2num(XMLtagContents("fract_y",atomSite))
			fracZ = str2num(XMLtagContents("fract_z",atomSite))
		endif
		if (numtype(fracX+fracY+fracZ) && strlen(WyckoffSymbol))		// try to get x,y,z from Wyckoff symbol
			if (ForceXYZtoWyckoff(SpaceGroupID,WyckoffSymbol,fracX,fracY,fracZ))
				fracX = NaN													// will cause a break in next if()
			endif
		endif
		if (numtype(fracX+fracY+fracZ))									// give up, cannot determine coordinates
			break
		endif
		mult = 0
		if (strlen(WyckoffSymbol)==0)									// try to set Wyckoff symbol from coordinates
			WyckoffSymbol = FindWyckoffSymbol(SpaceGroupID,fracX,fracY,fracZ,mult)
		endif

		valence = str2num(XMLtagContents("valence",atomSite))
		valence = (numtype(valence)==0 && abs(valence)<10) ? round(valence) : 0
		occupy = str2num(XMLtagContents("occupancy",atomSite))
		occupy = (occupy>=0 && occupy<=1) ? occupy : 1
		DebyeT = str2num(XMLtagContents("DebyeTemperature",atomSite))
		DebyeT = DebyeT>0 ? DebyeT : NaN
		unit = StringByKey("unit", XMLattibutes2KeyList("DebyeTemperature",atomSite),"=")
		unit = SelectString(strlen(unit),"K",unit)				// default Debye Temperature units are K
		DebyeT = ConvertTemperatureUnits(DebyeT,unit,"K")	// DebyeT is always stored as 

		Biso = str2num(XMLtagContents("B_iso",atomSite))
		Biso = Biso <=0 || numtype(Biso) ?  NaN : Biso
		unit = StringByKey("unit", XMLattibutes2KeyList("B_iso",atomSite),"=")
		Biso *= ConvertUnits2meters(unit,defaultLen=1e-20)*1e18	// want length in nm^2

		Uiso = str2num(XMLtagContents("U_iso",atomSite))
		Uiso = Uiso <=0 || numtype(Uiso) ?  NaN : Uiso
		unit = StringByKey("unit", XMLattibutes2KeyList("U_iso",atomSite),"=")
		Uiso *= ConvertUnits2meters(unit,defaultLen=1e-20)*1e18	// want length in nm^2

		aU11 = str2num(XMLtagContents("aniso_U_11",atomSite))
		unit = StringByKey("unit", XMLattibutes2KeyList("aniso_U_11",atomSite),"=")
		aU11 *= ConvertUnits2meters(unit,defaultLen=1e-20)*1e18	// want length in nm^2
		aU22 = str2num(XMLtagContents("aniso_U_22",atomSite))
		unit = StringByKey("unit", XMLattibutes2KeyList("aniso_U_22",atomSite),"=")
		aU22 *= ConvertUnits2meters(unit,defaultLen=1e-20)*1e18
		aU33 = str2num(XMLtagContents("aniso_U_33",atomSite))
		unit = StringByKey("unit", XMLattibutes2KeyList("aniso_U_33",atomSite),"=")
		aU33 *= ConvertUnits2meters(unit,defaultLen=1e-20)*1e18
		aU12 = str2num(XMLtagContents("aniso_U_12",atomSite))
		unit = StringByKey("unit", XMLattibutes2KeyList("aniso_U_12",atomSite),"=")
		aU12 *= ConvertUnits2meters(unit,defaultLen=1e-20)*1e18
		aU13 = str2num(XMLtagContents("aniso_U_13",atomSite))
		unit = StringByKey("unit", XMLattibutes2KeyList("aniso_U_13",atomSite),"=")
		aU13 *= ConvertUnits2meters(unit,defaultLen=1e-20)*1e18
		aU23 = str2num(XMLtagContents("aniso_U_23",atomSite))
		unit = StringByKey("unit", XMLattibutes2KeyList("aniso_U_23",atomSite),"=")
		aU23 *= ConvertUnits2meters(unit,defaultLen=1e-20)*1e18

		xtal.atom[N].name = atomLabel[0,59]
		xtal.atom[N].Zatom = Zatom
		xtal.atom[N].x = fracX
		xtal.atom[N].y = fracY
		xtal.atom[N].z = fracZ
		xtal.atom[N].valence = valence
		xtal.atom[N].occ = occupy
		xtal.atom[N].WyckoffSymbol = WyckoffSymbol[0]
		xtal.atom[N].mult = mult
		xtal.atom[N].DebyeT = DebyeT
		xtal.atom[N].Biso = Biso
		xtal.atom[N].Uiso = Uiso
		xtal.atom[N].U11 = numtype(aU11) ?  NaN : aU11
		xtal.atom[N].U22 = numtype(aU22) ?  NaN : aU22
		xtal.atom[N].U33 = numtype(aU33) ?  NaN : aU33
		xtal.atom[N].U12 = numtype(aU12) ?  NaN : aU12
		xtal.atom[N].U13 = numtype(aU13) ?  NaN : aU13
		xtal.atom[N].U23 = numtype(aU23) ?  NaN : aU23
		N += 1
	while(N<STRUCTURE_ATOMS_MAX)
	xtal.N = N												// number of atoms described here
	ForceXtalAtomNamesUnique(xtal)					// forces all of the xtal atom names to be unique

	Variable i, unitsConvert, Nlen, start0, Nbond=0
	String bondKeys, label0,label1,list, BondBody
	start = 0
	do
		start0 = start
		BondBody = XMLtagContents("bond_chemical",cif, start=start0)	// next chemical bond site
		if (strlen(BondBody)<1)
			break
		endif
		Wave blen = str2vec(BondBody)				// bond length(s)
		bondKeys = XMLattibutes2KeyList("bond_chemical",cif, start=start)	// next bond
		label0 = StringByKey("n0",bondKeys,"=")
		label1 = StringByKey("n1",bondKeys,"=")
		Nlen = min(numpnts(blen),5)					// number of lengths for this bond
		if (strlen(label0) && strlen(label1) && numtype(sum(blen))==0 && Nlen>0)
			xtal.bond[Nbond].label0 = label0[0,59]
			xtal.bond[Nbond].label1 = label1[0,59]
			xtal.bond[Nbond].N = Nlen
			unit = StringByKey("unit",bondKeys,"=")
			unitsConvert = ConvertUnits2meters(unit,defaultLen=1e-10)*1e9	// want length in nm
			for (i=0;i<Nlen;i+=1)
				xtal.bond[Nbond].len[i] = blen[i]*unitsConvert
			endfor
			Nbond += 1
		endif
	while(Nbond<(2*STRUCTURE_ATOMS_MAX))
	xtal.Nbonds = Nbond
	ForceLatticeToStructure(xtal)
	UpdateCrystalStructureDefaults(xtal)
	return 0
End
//
Static Function/T GetSymLinesFromXMLbuffer(buf)		// get the sym ops from an xml buffer
	String buf			// content of an xml file

	String symOps = XMLtagContents2List("symmetry_equiv_pos_as_xyz",buf)
	if (strlen(symOps)<1)
		symOps = XMLtagContents2List("space_group_symop_operation_xyz",buf)
	endif
	symOps = ReplaceString(" ",symOps,"")					// no spaces in the symmetry operations
	return symOps
End
//
Static Function ForceXtalAtomNamesUnique(xtal)		// forces all of the xtal atom names to be unique
	STRUCT crystalStructure &xtal

	Variable N=xtal.N
	Make/T/N=(N)/FREE all=xtal.atom[p].name	// holds all the names

	String namej, base, nameTest
	Variable i,j,num
	for (j=0;j<(N-1);j+=1)
		namej = all[j]								// check this name against others
		if (strlen(namej)<1)						// skip empty names
			continue
		elseif (countDuplicateNames(all,namej)>1)		// will found duplicates (always find 1)
			splitLabel(namej,base,num)
			if (numtype(num))						// try to change "Cu" -> "Cu1"
				num = 1
				nameTest = AddNum2Base(base,num)
				if (countDuplicateNames(all,nameTest)<1)	// base+"1" does not exist
					all[j] = nameTest
				endif
			endif

			for (i=j+1;i<N;i+=1)					// look for matches to namej
				if (StringMatch(all[i],namej))	// need to change all[i]
					do
						num += 1
						nameTest = AddNum2Base(base,num)
					while (countDuplicateNames(all,nameTest))
					all[i] = nameTest
				endif
			endfor
		endif
	endfor

	for (j=0;j<N;j+=1)
		xtal.atom[j].name = all[j]
	endfor
End
//
Static Function splitLabel(name,base,num)
	String name		// for name of form "Ab0001"
	String &base		// returns "Ab"
	Variable &num		// returns 1

	Variable i, N=strlen(name)
	for (i=N-1;i>=0;i-=1)			// find number at end
		if (!isdigit(name[i]))
			break
		endif
	endfor
	num = (i+1>=N) ? NaN : str2num(name[i+1,Inf])
	base = ""
	if (i>=0)
		base = name[0,i]
	endif
End
//
Static Function countDuplicateNames(ws,name)
	// counts how many times name appers in ws
	Wave/T ws						// a string wave
	String name					// a string to search for
	Variable N=DimSize(ws,0)
	Make/N=(N)/FREE dup=0
	dup = StringMatch(ws[p],name)
	return sum(dup)
End
//
Static Function/T AddNum2Base(base,num)
	String base
	Variable num
	if (numtype(num) || num<0)
		return base
	endif
	String name = base + num2istr(num)
	if (strlen(name)>59)			// too long, trim base
		String nStr = num2istr(num)
		Variable i = 59 - strlen(nStr)
		name = base[0,i]+nStr
	endif
	return name
End


Static Function/T crystalStructure2xml(xtal,NL)	// convert contents of xtal structure to xml sting (suitable for a file)
	STRUCT crystalStructure &xtal				// this sruct is printed in this routine
	String NL											// new line string (probably "\r" or "\n")

	String cif="<cif>"+NL
	String str, unit=" unit=\"nm\""

	if (strlen(xtal.desc))
		cif += "\t<chemical_name_common>"+xtal.desc+"</chemical_name_common>"+NL
	endif
	if (numtype(xtal.SpaceGroup)==0)
		cif += "\t<space_group_IT_number>"+num2istr(xtal.SpaceGroup)+"</space_group_IT_number>"+NL
	endif
	if (strlen(xtal.SpaceGroupID)>0)
		cif += "\t<space_group_id>"+xtal.SpaceGroupID+"</space_group_id>"+NL
	endif
	Variable alphaT = xtal.alphaT
	alphaT = alphaT>0 ? alphaT : NaN

	cif += "\t<cell>"+NL
	cif += "\t\t<a"+unit+">"+num2StrFull(xtal.a)+"</a>"+NL
	cif += "\t\t<b"+unit+">"+num2StrFull(xtal.b)+"</b>"+NL
	cif += "\t\t<c"+unit+">"+num2StrFull(xtal.c)+"</c>"+NL
	cif += "\t\t<alpha>"+num2StrFull(xtal.alpha)+"</alpha>"+NL
	cif += "\t\t<beta>"+num2StrFull(xtal.beta)+"</beta>"+NL
	cif += "\t\t<gamma>"+num2StrFull(xtal.gam)+"</gamma>"+NL
	if (xtal.Vc > 0)
		cif += "\t\t<volume>"+num2StrFull(xtal.Vc)+"</volume>"+NL
	endif
	if (numtype(xtal.Temperature)==0)
		Variable Temperature = xtal.Temperature
		String Tunits = "C"
		if (Temperature<-100)						// for low Temperatures, use Kelvin
			Tunits = "K"
			Temperature = ConvertTemperatureUnits(Temperature,"C","K")
		endif
		sprintf str, "<temperature unit=\"%s\">%s</temperature>",Tunits,num2StrFull(Temperature)
		cif += "\t\t"+str+NL
	endif
	if (numtype(alphaT)==0)
		cif += "\t\t<alphaT>"+num2StrFull(alphaT)+"</alphaT>"
		cif += "\t\t	<!-- a = ao*(1+alphaT*(TempC-22.5)) -->"+NL
	endif
	if (numtype(xtal.Unconventional00)==0)
		str="{"
		str += "{"+num2StrFull(xtal.Unconventional00)+","+num2StrFull(xtal.Unconventional10)+","+num2StrFull(xtal.Unconventional20)+","+"}, "
		str += "{"+num2StrFull(xtal.Unconventional01)+","+num2StrFull(xtal.Unconventional11)+","+num2StrFull(xtal.Unconventional21)+","+"}, "
		str += "{"+num2StrFull(xtal.Unconventional02)+","+num2StrFull(xtal.Unconventional12)+","+num2StrFull(xtal.Unconventional22)+","+"} }"
		cif += "\t\t<Unconventional>"+str+"</Unconventional>"+NL
		//	<Unconventional>{ {0,1,0}, {0,0,1}, {1,0,0} }</Unconventional>
	endif
	cif += "\t</cell>"+NL

	Variable i, itemp
	for(i=0; i<xtal.N; i+=1)
		cif += "\t<atom_site>"+NL
		cif += "\t\t<label>"+xtal.atom[i].name+"</label>"+NL
		cif += "\t\t<symbol>"+Z2symbol(xtal.atom[i].Zatom)+"</symbol>"+NL
		cif += "\t\t<fract_xyz>"+num2StrFull(xtal.atom[i].x)+" "+num2StrFull(xtal.atom[i].y)+" "+num2StrFull(xtal.atom[i].z)+"</fract_xyz>"+NL
		if (numtype(xtal.atom[i].valence)==0 && abs(xtal.atom[i].valence)<10 && abs(xtal.atom[i].valence)>0)
			cif += "\t\t<valence>"+num2StrFull(xtal.atom[i].valence)+"</valence>"+NL
		endif
		if (xtal.atom[i].occ<1 && numtype(xtal.atom[i].occ)==0)
			cif += "\t\t<occupancy>"+num2StrFull(xtal.atom[i].occ)+"</occupancy>"+NL
		endif
		if (strlen(xtal.atom[i].WyckoffSymbol)>0)
			cif += "\t\t<WyckoffSymbol>"+xtal.atom[i].WyckoffSymbol+"</WyckoffSymbol>"+NL
		endif

		itemp = atomThermalInfo(xtal.atom[i])
		if (itemp==1)								// write ONE kind of thermal vibartion info
			cif += "\t\t<DebyeTemperature unit=\"K\">"+num2StrFull(xtal.atom[i].DebyeT)+"</DebyeTemperature>"+NL
		elseif (itemp==2)
			cif += "\t\t<B_iso unit=\""+ARING+"^2\">"+num2StrFull(xtal.atom[i].Biso * 100)+"</B_iso>"+NL
		elseif (itemp==3)
			cif += "\t\t<U_iso unit=\""+ARING+"^2\">"+num2StrFull(xtal.atom[i].Uiso * 100)+"</U_iso>"+NL
		elseif (itemp>=4)
			cif += "\t\t<aniso_U_11 unit=\"nm^2\">"+num2StrFull(xtal.atom[i].U11)+"</aniso_U_11>"+NL
			cif += "\t\t<aniso_U_22 unit=\"nm^2\">"+num2StrFull(xtal.atom[i].U22)+"</aniso_U_22>"+NL
			cif += "\t\t<aniso_U_33 unit=\"nm^2\">"+num2StrFull(xtal.atom[i].U33)+"</aniso_U_33>"+NL
			if (itemp==5)
				cif += "\t\t<aniso_U_12 unit=\"nm^2\">"+num2StrFull(xtal.atom[i].U12)+"</aniso_U_12>"+NL
				cif += "\t\t<aniso_U_13 unit=\"nm^2\">"+num2StrFull(xtal.atom[i].U13)+"</aniso_U_13>"+NL
				cif += "\t\t<aniso_U_23 unit=\"nm^2\">"+num2StrFull(xtal.atom[i].U23)+"</aniso_U_23>"+NL
			endif
		endif
		cif += "\t</atom_site>"+NL
	endfor

	String wStr
	for(i=0; i<xtal.Nbonds; i+=1)
		Make/N=(xtal.bond[i].N)/FREE/D lens
		lens = xtal.bond[i].len[p]
		wStr = vec2str(lens,bare=1,sep=" ")
		sprintf str, "\t<bond_chemical%s n0=\"%s\" n1=\"%s\">%s</bond_chemical>",unit,xtal.bond[i].label0,xtal.bond[i].label1,wStr
		cif += str+NL
	endfor

	cif += "</cif>"+NL
	return cif
End

//  =================== End of reading/writing xml files ===================  //
//  ========================================================================= //



// =========================================================================
// =========================================================================
//	Start of xtl --> xml conversion

// convert one xtl file to a newer xml file, also moves *.xtl file to folder "materials/old_xtl_files"
Function ConverXTLfile2XMLfile(xtlName)
	String xtlName

	PathInfo materialsPath
	if (V_flag==0)										// make it if it does not exist
		NewPath/Z materialsPath, ParseFilePath(1,FunctionPath("CrystalsAreHere"),":",1,0)
	endif
	PathInfo materialsPath
	if (V_flag==0)										// make it if it does not exist
		NewPath/Z materialsPath, SpecialDirPath("Documents",0,0,0)+"materials"
		if (V_flag)										// path still does not exist, use Documents folder
			NewPath/Z materialsPath, SpecialDirPath("Documents",0,0,0)
		endif
	endif
	String fileFilters = "XTL Files (*.xtl):.xtl;old cri Files (*.cri):.cri;All Files:.*;"
	Variable f
	Open/D=2/R/F=fileFilters/M="File with Crystal Information"/P=materialsPath f as xtlName
	xtlName = S_fileName
	if (strlen(xtlName)==0)
		return 1
	endif

	// check file extension
	String extension = ParseFilePath(4,xtlName,":",0,0)
	if (!stringmatch(extension,"xtl"))
		DoAlert 1,"This does not look like an old '*.xtl' file, stop now?"
		if (V_flag!=2)
			return 1
		endif
	endif
	printf "xtlName = '%s'\r",xtlName

	String cif = convertOneXTL2XML(xtlName,NEW_LINE)
	if (strlen(cif)<1)
		DoAlert 0,"Unable to properly interpret file '"+xtlName+"'"
		return 1
	endif

	// wite the xml file
	String xmlName = ParseFilePath(1,xtlName,":",1,0)+ParseFilePath(3,xtlName,":",0,0)+".xml"
	Open/C="R*ch"/T="TEXT"/Z f as xmlName
	if (V_flag==0)
		fprintf f,  "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>"+NEW_LINE+NEW_LINE
		FBinWrite f, cif
		Close f
		Variable err = 0
	else
		DoAlert 0, "UNable to write to file '"+xmlName+"'"
		err = 1
	endif

	// move xtl file to folder "old_xtl_files",  make sure the folder "old_xtl_files" exists, then move the *.xtl file
	String oldXTLpath = ParseFilePath(1,xtlName,":",1,0)+"old_xtl_files:"
	NewPath/C/O/Q/Z old_xtl_files, oldXTLpath		// just use NewPath to create the folder
	KillPath/Z old_xtl_files
	NewPath/C/O/Q/Z original_xtl_path, ParseFilePath(1,xtlName,":",1,0)	// original path to the xtl file
	String name=ParseFilePath(0,xtlName,":",1,0)
	MoveFile/P=original_xtl_path name  as ":old_xtl_files:"+name			// actually move the file to the "old_xtl_files"
	KillPath/Z original_xtl_path
	return 0
End
//
Static Function/T convertOneXTL2XML(xtl,NL)	// read in the xtl file, returns string with contents for xml file
	String xtl
	String NL											// new line string (probably "\r" or "\n"

	String keyVals = keyStrFromFile(xtl,"CrystalStructure","home")
	if (strlen(keyVals)<1)
		return ""
	endif

	String sval, str, cif="<cif>"+NL

	sval = StringByKey("structureDesc",keyVals,"=")
	if (strlen(sval))
		cif += "\t<chemical_name_common>"+sval+"</chemical_name_common>"+NL
	endif
	sval = StringByKey("SpaceGroup",keyVals,"=")
	if (numtype(str2num(sval))==0)
		cif += "\t<space_group_IT_number>"+sval+"</space_group_IT_number>"+NL
	endif

	String unit= StringByKey("lengthUnit",keyVals,"=")	// length unit for lattice constants a,b,c
	if (strlen(unit))
		unit = " unit=\""+unit+"\""
	endif

	Variable alphaT = NumberByKey("latticeAlphaT",keyVals,"=")
	Variable temperature = numtype(alphaT) ? NaN : 22.5	// temperature (C)
	if (numtype(NumberByKey("temperatureK",keyVals,"="))==0)
		temperature = NumberByKey("temperatureK",keyVals,"=") - 273.15
	endif

	str = StringByKey("latticeParameters",keyVals,"=")
	Variable a,b,c,alpha,bet,gam
	sscanf str,"{%g, %g, %g, %g, %g, %g}",a,b,c,alpha,bet,gam
	if (V_flag!=6)
		return ""
	endif
	cif += "\t<cell>"+NL
	cif += "\t\t<a"+unit+">"+num2StrFull(a)+"</a>"+NL
	cif += "\t\t<b"+unit+">"+num2StrFull(b)+"</b>"+NL
	cif += "\t\t<c"+unit+">"+num2StrFull(c)+"</c>"+NL
	cif += "\t\t<alpha>"+num2StrFull(alpha)+"</alpha>"+NL
	cif += "\t\t<beta>"+num2StrFull(bet)+"</beta>"+NL
	cif += "\t\t<gamma>"+num2StrFull(gam)+"</gamma>"+NL
	// cif += "\t\t<volume>"+num2StrFull(volume)+"</volume>"
	if (numtype(temperature)==0)
		cif += "\t\t<temperature unit=\"C\">"+num2StrFull(temperature)+"</temperature>"+NL
	endif
	if (numtype(alphaT)==0)
		cif += "\t\t<alphaT>"+num2StrFull(alphaT)+"</alphaT>"
		cif += "\t\t	<!-- a = ao*(1+alphaT*(TempC-"+SelectString(numtype(temperature),num2StrFull(temperature),"22.5")+")) -->"+NL
	endif
	str = StringByKey("Unconventional",keyVals,"=")
	if (strlen(str)>0)
		cif += "\t\t<Unconventional>"+str+"</Unconventional>"+NL
	endif
	cif += "\t</cell>"+NL

	Variable DebyeTemperature = NumberByKey("DebyeTemp",keyVals,"=")

	String atomLabel, symbol
	Variable xFraction,yFraction,zFraction,occupancy, i=1
	do
		str = StringByKey("AtomDesctiption"+num2istr(i),keyVals,"=")
		if (strlen(str)<1)
			break
		endif
		str = SingleWhiteSpace(str)

		sscanf str, "{%s %g %g %g %g}", atomLabel,xFraction,yFraction,zFraction,occupancy
		if (V_flag<4)
			break
		endif
		occupancy = V_flag<5 ? 1 : occupancy

		cif += "\t<atom_site>"+NL
		symbol = SelectString(isdigit(atomLabel[1]),atomLabel[0,1],atomLabel[0])
		cif += "\t\t<label>"+atomLabel+"</label>"+NL
		cif += "\t\t<symbol>"+symbol+"</symbol>"+NL
//		cif += "\t\t<fract_x>"+num2StrFull(xFraction)+"</fract_x>"+NL
//		cif += "\t\t<fract_y>"+num2StrFull(yFraction)+"</fract_y>"+NL
//		cif += "\t\t<fract_z>"+num2StrFull(zFraction)+"</fract_z>"+NL
		cif += "\t\t<fract_xyz>"+num2StrFull(xFraction)+" "+num2StrFull(yFraction)+" "+num2StrFull(zFraction)+"</fract_xyz>"+NL
		if (occupancy<1)
			cif += "\t\t<occupancy>"+num2StrFull(occupancy)+"</occupancy>"+NL
		endif
		if (DebyeTemperature>0)
			cif += "\t\t<DebyeTemperature unit=\"K\">"+num2StrFull(DebyeTemperature)+"</DebyeTemperature>"+NL
		endif

		cif += "\t</atom_site>"+NL
		i += 1
	while (1)

	String citation=""
	str = StringByKey("citation",keyVals,"=")
	if (strlen(str)>0)
		citation += str + ";"
	endif
	str = StringByKey("citation1",keyVals,"=")
	if (strlen(str)>0)
		citation += str + ";"
	endif
	str = StringByKey("citation2",keyVals,"=")
	if (strlen(str)>0)
		citation += str + ";"
	endif
	str = StringByKey("citation3",keyVals,"=")
	if (strlen(str)>0)
		citation += str + ";"
	endif
	if (strlen(citation)>0)
		cif += "\t<citation>"+citation+"</citation>"+NL
	endif

	// treat $remarks as comments
	str = StringByKey("remarks",keyVals,"=")
	if (strlen(str)>0)
		cif += "\t<!-- "+str+" -->"+NL
	endif

	// transfer any old comment lines
	GetFileFolderInfo/P=home/Q/Z xtl
	Variable f
	if (V_isFile)
		Open/R/Z=2/P=home f as S_Path
		if (V_flag==0)
			String line=" "
			do
				FReadLine f, line
				line = TrimFrontBackWhiteSpace(line)
				if (strsearch(line,"//",0)==0 && strlen(line)>2)
					cif += "\t<!-- "+line[2,Inf]+" -->"+NL
				endif
			while (strlen(line))
			Close f
		endif
	endif
	cif += "</cif>"+NL
	return cif
End
//
ThreadSafe Static Function/T num2StrFull(val)
	Variable val
	Variable i = placesOfPrecision(val)
	Variable absVal = abs(val)
	i = (absVal>=10 && absVal<1e6) ? max(i,1+floor(log(absVal))) : i
	String str, fmt
	sprintf fmt, "%%.%dg",i
	sprintf str,fmt,val
	return str
End
//
ThreadSafe Static Function/T SingleWhiteSpace(str,[sep,chars,trim])
	String str			// string to change
	String sep			// separator to use
	String chars		// an optional list of characters to change to white space, default is space & tab
	Variable trim		// if true, trim leading or trailing sep

	if (ParamIsDefault(sep))					// default separator is a single space
		sep = " "
	endif
	if (ParamIsDefault(chars))					// default is to treat both sapce and tab as white space
		chars = "\t"
	endif
	trim = ParamIsDefault(trim) ? 1 : trim
	trim = numtype(trim) ? 1 : trim
	sep = sep[0]								// sep can only contain 1 character
	chars = ReplaceString(sep,chars,"")		// remove sep, since space is the default separator

	Variable i,N=strlen(chars)
	for (i=0;i<N;i+=1)							// replace all occurances of characters in chars with sep
		str = ReplaceString(chars[i],str,sep)
	endfor

	do
		str = ReplaceString(sep+sep,str,sep)	// replace all runs of sep with a single sep char
	while(strsearch(str,sep+sep,0)>=0)

	if (trim)									// trim any leading or trailing sep characters
		if (stringmatch(str[0],sep))
			str = str[1,Inf]						// trim a leading sep character
		endif
		i = strlen(str)-1
		if (stringmatch(str[i],sep))
			str = str[0,i-1]					// trim a trailing sep character
		endif
	endif

	return str
End

//	End of xtl --> xml conversion
// =========================================================================
// =========================================================================



// =========================================================================
// =========================================================================
//	Start of CIF read

Static Function readCrystalStructureCIF(xtal,fileName,[path])
	STRUCT crystalStructure &xtal						// this sruct is printed in this routine
	String fileName
	String path
	if (ParamIsDefault(path))
		path = "home"
	endif

	Variable f
	String fileFilters = "CIF Files (*.cif):.cif;TEXT Files (*.txt):.txt;All Files:.*;"
	Open/R/F=fileFilters/M="Select CIF file with Crystal Description"/P=$path/R/Z=2 f as fileName
	if (strlen(S_fileName)<1 || V_flag)
		return 1
	endif
	String fullFile = S_fileName

	FStatus f
	String buf=PadString("",V_logEOF,0x20)
	FBinRead f, buf
	Close f
	buf = ReplaceString("\r\n",buf,"\n")
	buf = ReplaceString("\n\r",buf,"\n")
	buf = ReplaceString("\r",buf,"\n")

	Variable i = strsearch(buf,"data_",0)		// find the _data section
	if (i>0)
		i = strsearch(buf,"\ndata_",0)			// _data is not first, it must start a line
	endif
	if (i<0)
		DoAlert 0,"Cannot find 'data_*' to start"
		return 1
	endif
	buf = buf[i,Inf]
	i = strsearch(buf,"\n",0)
	String desc = buf[5,i-1]

	if (StringMatch(buf[0,11], "data_general"))	// skip the "data_general"
		i = strsearch(buf,"\ndata_",0)		// goto next one
		buf = buf[i+1,Inf]
		i = strsearch(buf,"\n",0)
		desc = buf[5,i-1]
	endif
	buf = buf[i,Inf]

	i = strsearch(buf,"\ndata_",0)			// check for second data_ section
	if (i>1)
		buf = buf[0,i-1]					// trim off any following "data_" sections
	endif

	if (strlen(buf)<1)
		return 1
	endif

	Variable err = CIF_interpret(xtal,buf,desc=desc)
	if (!err)
		Variable i0,i1
		i1 = strlen(fullFile)-1				// possibly trim file length to fit, keep last part of fullFile
		i0 = max(0,i1-MAX_FILE_LEN)
		fullFile = fullFile[i0,i1]
		xtal.sourceFile = fullFile
		UpdateCrystalStructureDefaults(xtal)
	endif
	return err
End
//
Static Function CIF_interpret(xtal,buf,[desc])
	STRUCT crystalStructure &xtal
	String buf
	String desc
	desc = SelectString(ParamIsDefault(desc),desc,"")

	buf += "\n"
	String name=""
	name = CIF_readString("_chemical_formula_structural",buf)
	if (strlen(name)<1)
		name = CIF_readString("_chemical_name_systematic",buf)
	endif
	if (strlen(name)<1)
		name = CIF_readString("_chemical_name_mineral",buf)
	endif
	String str = SelectString(strlen(name),desc,name)
	xtal.desc = str[0,99]

	// find lattice constants
	xtal.a = CIF_readNumber("_cell_length_a",buf)/10		// want nm, cif is in Angstroms
	xtal.b = CIF_readNumber("_cell_length_b",buf)/10
	xtal.c = CIF_readNumber("_cell_length_c",buf)/10
	xtal.alpha = CIF_readNumber("_cell_angle_alpha",buf)
	xtal.beta = CIF_readNumber("_cell_angle_beta",buf)
	xtal.gam = CIF_readNumber("_cell_angle_gamma",buf)
	if (numtype(xtal.a + xtal.b + xtal.c + xtal.alpha + xtal.beta + xtal.gam))
		return 1
	endif

	// 1st try for the id directly
	xtal.SpaceGroupID = CIF_readString("_space_group_id",buf)
	// 2nd try using the H-M symbol to set the SpaceGroupID & SpaceGroupIDnum
	if (strlen(xtal.SpaceGroupID)<1)
		String HMsym = CIF_readString("_symmetry_space_group_name_H-M",buf)	//	_symmetry_space_group_name_H-M 'I 1 2/a 1'
		String nlist = SymString2SGtype(HMsym,2,0)			// checks in getHMsym2
		if (strlen(nlist)<1)
			nlist = SymString2SGtype(HMsym,2,1)				// checks in getHMsym2 again, but ignoring minus signs
		endif
		if (ItemsInList(nlist)==1)								// just one result, use it
			Variable idNum = str2num(StringFromList(0,nlist))
			if (isValidSpaceGroupIDnum(idNum))					// found valid SpaceGroupIDnum
				String allIDs=MakeAllIDs()
				xtal.SpaceGroupID = StringFromList(idNum-1, allIDs)
				printf "Setting Space Group from H-M = \"%s\"\r", HMsym
			endif
		endif
	endif
	// 3rd try, look at the symmetry operations
	if (strlen(xtal.SpaceGroupID)<1)
		String symLine = GetSymLinesFromCIFbuffer(buf)	// get the sym ops from a CIF file buffer
		xtal.SpaceGroupID = FindIDfromSymOps(symLine)
	endif
	// 4th try, look for the Space Group number [1,230], and use the default id for that Space Group
	if (strlen(xtal.SpaceGroupID)<1)
		Variable SG=CIF_readNumber("_symmetry_Int_Tables_number",buf)
		if (SG>0)
			xtal.SpaceGroupID = FindDefaultIDforSG(SG)		// try to get id from SG [1,230]
		endif
	endif
	if (!isValidSpaceGroupID(xtal.SpaceGroupID))			// give up
		Abort "cannot find the Space Group from CIF file"
	endif
	xtal.SpaceGroup = str2num(xtal.SpaceGroupID)
	xtal.SpaceGroupIDnum = SpaceGroupID2num(xtal.SpaceGroupID)

	xtal.Temperature = CIF_readNumber("_cell_measurement_temperature",buf)-273.15 	// temperature (K) for cell parameters

	xtal.N = 0
	// find the atoms
	Variable i, cc
	String loop_, line, list=""							// find loop_ with the fractional atom positions
	for (i=strsearch(buf,"\nloop_",0); i>=0; i=strsearch(buf,"\nloop_",i+5))
		list = CIF_loop_Labels(buf[i,Inf])
		if (WhichListItem("_atom_site_fract_x",list)>=0)
			break
		endif
	endfor
	if (strlen(list))				// found atom positions, interpret the list
		loop_ = buf[i+1,Inf]								// loop_ now starts after the "\n" preceeding loop_
		i = 0
		line = NextLineInBuf(loop_,i)					// skip the line with "loop_"
		do
			line = TrimBoth(NextLineInBuf(loop_,i))	// trim off any white space and the terminating "\n"
		while (char2num(line[0])==95)

		for(; strlen(line)<1 || char2num(line[0])==35;)	// skip any blank lines or lines starting with "#"
			line = TrimBoth(NextLineInBuf(loop_,i))	
		endfor

		String symb											// Element symbol, e.g. "Fe", "Co", ...
		Variable Z, occ, Biso, Uiso, Uij, mult, N
		for (N=0; N<STRUCTURE_ATOMS_MAX && strlen(line)>0; N+=1)	// loop until you get an empty line
			line = ChangeCIFline2List(line)			// make line semi-colon separated list
			xtal.atom[N].DebyeT = NaN
			xtal.atom[N].Uiso = NaN	;	xtal.atom[N].Biso = NaN
			xtal.atom[N].U11 = NaN		;	xtal.atom[N].U22 = NaN	;		xtal.atom[N].U33 = NaN
			xtal.atom[N].U12 = NaN		;	xtal.atom[N].U13 = NaN	;		xtal.atom[N].U23 = NaN
			xtal.atom[N].valence = 0

			name = StringFromList(WhichListItem("_atom_site_label",list),line)
			xtal.atom[N].name = name[0,59]
			symb = StringFromList(WhichListItem("_atom_site_type_symbol",list),line)
			symb = SelectString(strlen(symb),name,symb)	// if symbol not included, try to get Z from label
			Z = ZfromLabel(symb)
			xtal.atom[N].Zatom = Z>0 ? Z : 1
			xtal.atom[N].x = str2num(StringFromList(WhichListItem("_atom_site_fract_x",list),line))
			xtal.atom[N].y = str2num(StringFromList(WhichListItem("_atom_site_fract_y",list),line))
			xtal.atom[N].z = str2num(StringFromList(WhichListItem("_atom_site_fract_z",list),line))
			if (numtype(xtal.atom[N].x + xtal.atom[N].y + xtal.atom[N].z) || strlen(name)<1)
				break
			endif
			str = StringFromList(WhichListItem("_atom_site_Wyckoff_symbol",list),line)
			xtal.atom[N].WyckoffSymbol = str[0]
			occ = str2num(StringFromList(WhichListItem("_atom_site_occupancy",list),line))
			xtal.atom[N].occ = occ >=0 ? limit(occ,0,1) : 1
			mult = str2num(StringFromList(WhichListItem("_atom_site_symmetry_multiplicity",list),line))
			xtal.atom[N].mult = mult>0 ? mult : 1
			Biso = str2num(StringFromList(WhichListItem("_atom_site_B_iso_or_equiv",list),line))/100	// assume value in Angstrom^2
			xtal.atom[N].Biso = Biso>0 ? Biso : NaN

			Uiso = str2num(StringFromList(WhichListItem("_atom_site_U_iso_or_equiv",list),line))/100	// assume value in Angstrom^2
			xtal.atom[N].Uiso = Uiso>0 ? Uiso : NaN

			xtal.atom[N].U11 = str2num(StringFromList(WhichListItem("_atom_site_aniso_U_11",list),line))/100	// assume value in Angstrom^2
			xtal.atom[N].U22 = str2num(StringFromList(WhichListItem("_atom_site_aniso_U_22",list),line))/100
			xtal.atom[N].U33 = str2num(StringFromList(WhichListItem("_atom_site_aniso_U_33",list),line))/100
			xtal.atom[N].U12 = str2num(StringFromList(WhichListItem("_atom_site_aniso_U_12",list),line))/100
			xtal.atom[N].U13 = str2num(StringFromList(WhichListItem("_atom_site_aniso_U_13",list),line))/100
			xtal.atom[N].U23 = str2num(StringFromList(WhichListItem("_atom_site_aniso_U_23",list),line))/100

			line = TrimBoth(NextLineInBuf(loop_,i))		// get the next line
		endfor
		xtal.N = N
	endif

	xtal.Nbonds	 = 0
	xtal.Unconventional00 = NaN
	xtal.sourceFile = ""
	ForceLatticeToStructure(xtal)
	return 0
End

Function/T GetSymLinesFromCIFbuffer(buf)		// get the sym ops from a CIF file buffer
	String buf			// content of a CIF file
	// find loop_ with symmetry ops

	String list="", symTag="", loop_, line, str
	Variable i, N
	for (i=strsearch(buf,"\nloop_",0); i>=0; i=strsearch(buf,"\nloop_",i+5))
		list = CIF_loop_Labels(buf[i,Inf])
		if (WhichListItem("_symmetry_equiv_pos_as_xyz",list)>=0)
			symTag = "_symmetry_equiv_pos_as_xyz"
			break
		elseif (WhichListItem("_space_group_symop_operation_xyz",list)>=0)
			symTag = "_space_group_symop_operation_xyz"
			break
		endif
	endfor
	Variable isym = strlen(list) ? WhichListItem(symTag,list) : -1
	if (isym>=0)												// found symmetry operations, get them
		loop_ = buf[i+1,Inf]								// loop_ now starts after the "\n" preceeding loop_
		i = 0
		line = NextLineInBuf(loop_,i)					// skip the line with "loop_"
		do
			line = TrimBoth(NextLineInBuf(loop_,i))	// trim off any white space and the terminating "\n"
		while (char2num(line[0])==95)

		for(; strlen(line)<1 || char2num(line[0])==35;)	// skip any blank lines or lines starting with "#"
			line = TrimBoth(NextLineInBuf(loop_,i))	
		endfor

		String symLine=""
		for (N=0; N<192 && !CIFloopEnd(line); N+=1)// loop until you get an empty line, over all symmetry op
			line = ChangeCIFline2List(line)				// make line semi-colon separated list
			str = StringFromList(isym,line)
			str = ReplaceString(" ",str,"")
			if (strlen(str))
				symLine += str+";"
			endif
			line = TrimBoth(NextLineInBuf(loop_,i))	// get the next line
		endfor
		// printf "found %g sym ops: %s",ItemsInList(symLine), symLine[0,130]
		// printf "%s\r", SelectString(strlen(symLine)>131, "", "...")
	endif

	return symLine
End
//
Static Function/T FindIDfromSymOps(symList)
	String symList
	if (strlen(symList)<1)
		return ""
	endif

	Variable idNum, N=ItemsInList(symList)
	String str, allIDs=MakeAllIDs(), id
	for (idNum=1;idNum<=530;idNum+=1)
		id = StringFromList(idNum-1, allIDs)
		if (SymOpsMatchesID(id,symList))
			return id					// found a match
		endif
	endfor
	return ""							// give up, nothing matches
End
//
Static Function SymOpsMatchesID(id,symList,[printIt])
	// check that the given list of sym ops match my list using id
	// returns 1 if a match, 0 if not match
	String id				// SpaceGroupID
	String symList
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt

	String internal = setSymLineID(id)
	Variable i,N=ItemsInList(internal)
	if (ItemsInList(symList) != N)
		return 0											// number of operations differ
	endif 

	Make/N=(N)/D/FREE internalV=NaN, symV=NaN
	internalV = symOpList2number(StringFromList(p,internal))
	symV = symOpList2number(StringFromList(p,symList))
	Sort internalV,internalV
	Sort symV,symV
	Variable match = EqualWaves(symV,internalV,1)
	if (!match && printIt)
		print "Mismatch in sym ops"
		printf "  given ops:  %s\r",symList
		printf "  but %s has:  %s\r",id,internal
	endif
	return match
End
//
Static Function symOpList2number(symOps)	// takes a list of symOps, e.g. "x,y,z;-x,y,-z;x+1/2,y+1/2,z;-x+1/2,y+1/2,-z" and returns a unique number
	String symOps
	String op
	Variable crc,out,i,N=ItemsInList(symOps)
	for (i=0,out=0; i<N; i+=1)
		op = StringFromList(i,symOps)							// each sym op in symOps
		crc = expr2number(StringFromList(0,op,","),0)		// x term
		crc = expr2number(StringFromList(1,op,","),crc)	// y term
		crc = expr2number(StringFromList(2,op,","),crc)	// z term
		out += crc														// add because the order of the ops does not matter
	endfor
	return out
End
//
Static Function expr2number(expression,crc)
	String expression
	Variable crc
	Variable mx,my,mz,b
	ParseOneSymEquation(expression, mx,my,mz,b)	// sets mx,my,mz,b
	String str
	sprintf str, "%d%d%d%.4f", mx,my,mz,b
	return StringCRC(crc,str)
End
//
Static Function CIFloopEnd(line)
	String line
	if (strlen(line)<1)
		return 1
	elseif (char2num(line[0])==95)
		return 1
	elseif (strsearch(line,"loop_",0)==0)
		return 1
	endif
	return 0
End
//
Static Function/S ChangeCIFline2List(line)	// change a CIF data line to a semi-colon separated list
	String line
	line = TrimFront(line)							// remove any indents
	line = ReplaceString(";",line,":")			// semi-colon is special
	line = ReplaceString("\t",line," ")

	Variable i0=0,isp,iqu
	do
		isp = strsearch(line," ",i0)
		iqu = strsearch(line,"'",i0)
		isp = isp<0 ? Inf : isp					// turn not found (=-1) into Inf
		iqu = iqu<0 ? Inf : iqu
		if (isp<iqu)									// found both, but space is first
			line[isp,isp] = ";"
			i0 = isp + 1
		elseif (iqu<isp)								// found a quote separator
			if (iqu==0)
				line = line[1,Inf]
			else			
				line[iqu,iqu] = ";"
			endif
			iqu = strsearch(line,"'",iqu+1)		// find closing quote
			line[iqu,iqu] = ";"
			i0 = iqu + 1
		else
			break											// did not find space or quote, done
		endif
	while(1)

	do
		line = ReplaceString(";;",line,";")	// change all multiple ";" to single ";"
	while (strsearch(line,";;",0)>=0)
	return line
End
//Function/S ChangeCIFline2List(line)				// change a CIF data line to a semi-colon separated list
//	String line
//	line = TrimBoth(line)	
//	line = ReplaceString("\t",line,";")
//	line = ReplaceString(" ",line,";")
//	do
//		line = ReplaceString(";;",line,";")		// change all multiple ";" to single ";"
//	while (strsearch(line,";;",0)>=0)
//	return line
//End
//
Static Function/T CIF_loop_Labels(buf)
	String buf

	Variable i1, i0 = strsearch(buf,"\nloop_\n",0)
	i0 += 7
	String name, list=""
	for(;char2num(buf[i0])==95 || char2num(buf[i0])<=32;)
		if (char2num(buf[i0])<=32)
			i0 += 1							// skip leading white space
			continue
		endif
			i1 = strsearch(buf,"\n",i0+1)-1
		name = buf[i0,i1]
		list += name+";"
		i0 = i1 + 2
	endfor
	return list
End
//
Static Function CIF_readNumber(key,buf)
	String key
	String buf

	Variable len=strlen(key)
	Variable i0 = strsearch(buf,"\n"+key,0)
	if (i0<0)
		return NaN
	endif
	Variable i1 = strsearch(buf,"\n",i0+1)
	i1 = i1<i0 ? Inf : i1
	return str2num(buf[i0+len+1,i1-1])
End
//
Static Function/T CIF_readString(key,buf)
	String key
	String buf

	Variable i0 = strsearch(buf,"\n"+key,0), i1,i2
	if (i0<0)
		return ""
	endif

	i1 = strsearch(buf,"\"",i0+1)
	i2 = strsearch(buf,"'",i0+1)
	i1 = i1<0 ? Inf : i1
	i2 = i2<0 ? Inf : i2
	i0 = min(i1,i2)+1
	if (i1>i2)								// quoted with a single quote
		i1 = strsearch(buf,"'",i0+1)-1
	else
		i1 = strsearch(buf,"\"",i0+1)-1
	endif
	return TrimFrontBackWhiteSpace(buf[i0,i1])
End

//	End of CIF read
// =========================================================================
// =========================================================================



// =========================================================================
// =========================================================================
//	Start of crystal symmetry stuff

Static Constant TRICLINIC=0,MONOCLINIC=1,ORTHORHOMBIC=2,TETRAGONAL=3,TRIGONAL=4,HEXAGONAL=5,CUBIC=6
Static Constant P_CENTER=0,F_CENTER=1,B_CENTER=2,RHOMBOHEDRAL=3,C_CENTER=4,A_CENTER=5
Static Constant FCC=225,BCC=229,DIA=227,SIMPLE_CUBIC=195,SAPPHIRE=167	// generic Space Group numbers
Strconstant LatticeSystemNames="Triclinic;Monoclinic;Orthorhombic;Tetragonal;Trigonal;Hexagonal;Cubic"


Function/C getFstruct(h,k,l,[keV,T,printIt])	// user interface to getting F for current crystal structure
	Variable h,k,l
	Variable keV
	Variable T
	Variable printIt
	T = ParamIsDefault(T) ? NaN : T
	keV = ParamIsDefault(keV) ? NaN : keV
	printIt = ParamIsDefault(printIt) ? NaN : printIt
	printIt = numtype(printIt) ? (strlen(GetRTStackInfo(2))==0) : !(!printIt)

	STRUCT crystalStructure xtal						// this sruct is filled  by this routine
	FillCrystalStructDefault(xtal)
	Variable haveCromer = (exists("Get_f")==6)		// calculate fatom using Cromer-Liberman
	Variable haveDebyeT =  xtalhasDebye(xtal)
	if (numtype(h+k+l))
		h = numtype(h) ? 0 : h
		k = numtype(k) ? 0 : k
		l = numtype(l) ? 4 : l
		keV = keV>0 ? keV : NaN
		keV = haveCromer ? keV : NaN
		T = T>-273 ? T : xtal.Temperature
		T = haveDebyeT ? T : NaN
		Prompt h,"H"
		Prompt k,"K"
		Prompt l,"L"
		Prompt keV,"Energy (keV)"
		Prompt T,"Temperature (C)"
		if (haveCromer && haveDebyeT)
			DoPrompt "(hkl)",h,k,l,keV,T
		elseif (haveDebyeT)
			DoPrompt "(hkl)",h,k,l,T
		elseif (haveCromer)
			DoPrompt "(hkl)",h,k,l,keV,T
		else
			DoPrompt "(hkl)",h,k,l
		endif
		if (V_flag)
			return cmplx(NaN,NaN)
		endif
		if (numtype(h+k+l))
			return cmplx(NaN,NaN)
		endif
		keV = keV>0 ? keV : NaN
		T = T>-273.15 ? T : NaN

		printf "%sgetFstruct(%g,%g,%g",BULLET,h,k,l
		if (keV>0)
			printf ", keV=%g",keV
		endif
		if (T>0)
			printf ", T=%g",T
		endif
		printf ")\r"
		printIt = 1
	endif

	Variable T_K=ConvertTemperatureUnits(T,"C","K")	// Need T in K, but it was passed as C
	Variable/C Fc = Fstruct(xtal,h,k,l,keV=keV,T_K=T_K)
	Variable Freal=real(Fc),Fimag=imag(Fc)
	if (printIt)
		printf "F('%s', %d %d %d",xtal.desc,h,k,l
		if (keV>0)
			printf ", %gkeV",keV
		endif
		if (T>0)
			printf ", %g%sC",T,DEGREESIGN
		endif
		printf ") = %g %s i%g,      |F|^2 = %g\r",real(Fc),SelectString(imag(Fc)<0,"+","-"),abs(imag(Fc)),magsqr(Fc)
	endif
	return Fc
End
//Function/C getFstruct(h,k,l)				// user interface to getting F for current crystal structure
//	Variable h,k,l
//	if (numtype(h+k+l))
//		h = numtype(h) ? 0 : h
//		k = numtype(k) ? 0 : k
//		l = numtype(l) ? 2 : l
//		Prompt h,"H"
//		Prompt k,"K"
//		Prompt l,"L"
//		DoPrompt "(hkl)",h,k,l
//		if (V_flag)
//			return cmplx(NaN,NaN)
//		endif
//		if (numtype(h+k+l))
//			return cmplx(NaN,NaN)
//		endif
//	endif
//	STRUCT crystalStructure xtal						// this sruct is filled  by this routine
//	FillCrystalStructDefault(xtal)
//
//	Variable/C Fc = Fstruct(xtal,h,k,l)
//	Variable Freal=real(Fc),Fimag=imag(Fc)
//	if (strlen(GetRTStackInfo(2))==0)
//		printf "F('%s', %d %d %d) = %g %s i%g,      |F|^2 = %g\r",xtal.desc,h,k,l,real(Fc),SelectString(imag(Fc)<0,"+","-"),abs(imag(Fc)),magsqr(Fc)
//	endif
//	return Fc
//End


Function/C Fstruct(xtal,h,k,l,[keV,T_K])
	STRUCT crystalStructure &xtal				// sruct defining the crystal, not changed here
	Variable h,k,l
	Variable keV
	Variable T_K										// Temperature (K), used to calculate Debye-Waller factor
	keV = ParamIsDefault(keV) || keV<0 ? NaN : keV
	T_K = ParamIsDefault(T_K) ? NaN : max(0,T_K)

	String SpaceGroupID = xtal.SpaceGroupID
	if (!isValidSpaceGroupID(SpaceGroupID) || numtype(h+k+l))
		return cmplx(NaN,NaN)	// bad inputs
	endif
	Variable/C zero=cmplx(0,0)
	Variable usingHexAxes = (abs(90-xtal.alpha)+abs(90-xtal.beta)+abs(120-xtal.gam))<1e-6
	Variable system = latticeSystem(SpaceGroupID)

	if (!(mod(h,1) || mod(k,1) || mod(l,1)))	// non-integral, always allowed
		String sym = getHMsym2(xtal.SpaceGroupIDnum)	// get symmetry symbol
		strswitch (sym[0,0])
			case "F":
				if (!ALLOW_FC(h,k,l))
					return zero
				endif
				break
			case "I":
				if (!ALLOW_BC(h,k,l))
					return zero
				endif
				break
			case "C":
				if (!ALLOW_CC(h,k,l))
					return zero
				endif
				break
			case "A":
				if (!ALLOW_AC(h,k,l))
					return zero
				endif
				break
			case "R":
				if (usingHexAxes)
					if (!ALLOW_RHOM_HEX(h,k,l))		// rhombohedral cell with hexagonal axes
						return zero
					endif
				endif
				break											// using rhombohedral axes, so all are allowed
		endswitch
		if (system==HEXAGONAL)
			if (!ALLOW_HEXAGONAL(h,k,l)) 
				return zero
			endif
		endif
	endif

	Make/N=3/D/FREE hkl={h,k,l}
	Variable Qmag = 2*PI/dSpacing(xtal,h,k,l)	// |Q| vector (nm)
	if (numtype(xtal.atom[0].U11)==0)				// need Q-vector
		Wave recip = recipFrom_xtal(xtal)			// get reicprocal lattice from xtal
		MatrixOP/FREE qvec = recip x hkl
	endif
	keV = keV>0 ? keV : NumVarOrDefault("root:Packages:Lattices:keV",NaN)

	Variable/C fatomC
	Variable fatomMag, fatomArg=0
	String name
	Variable valence
	Variable m
	Variable/C c2PI=cmplx(0,2*PI), ifatomArg
	Variable/C Fc=cmplx(0,0)							// the result, complex structure factor

	reMakeAtomXYZs(xtal)
	if (!(xtal.N>=1))
		return cmplx(1,0)
	endif
	FUNCREF Get_f_proto fa = $"Get_f"				// if Get_f() does not exist, then Get_f_proto() will be used

	for (m=0;m<xtal.N;m+=1)							// loop over the defined atoms
		name="root:Packages:Lattices:atom"+num2istr(m)
		valence = numtype(xtal.atom[m].valence) ? 0 : xtal.atom[m].valence
		Wave ww = $name
		if (valence==0 )
			fatomC = fa(Z2symbol(xtal.atom[m].Zatom),Qmag/10, keV)
		else
			fatomC = fa(Z2symbol(xtal.atom[m].Zatom),Qmag/10, keV,valence=valence)
		endif
		fatomC = r2polar(fatomC)
		fatomMag = real(fatomC)
		fatomArg = imag(fatomC)
		fatomMag *= xtal.atom[m].occ

		if (xtal.Vibrate)
			Variable amu, DW, thetaM								// thetaM = Debye Temperature (K)
			Variable Biso, Uiso, itemp
			//	M = B*(sin^2(theta)/lam^2) --> M = B/(16 PI^2) * Q^2
			//	exp(-B * sin^2(theta)/lam^2)		B = 8 * PI^2 * <u^2> =  8 * PI^2 * Uiso
			itemp = atomThermalInfo(xtal.atom[m],T_K=T_K)	// Get kind of thermal info present (0 means none)
			if (itemp==1 && T_K>1e-9)								// has a DebyeWaller factor and non-zero Temperature
				thetaM = xtal.atom[m].DebyeT
				amu = Element_amu(xtal.atom[m].Zatom)			// mass of atom (amu)
				DW = exp(-DW_factor_M(T_K,thetaM,Qmag,amu))// calculates the M in exp(-M), no I/O
				fatomMag *= DW>0 ? DW : 1
			elseif (itemp==2)
				Biso = xtal.atom[m].Biso							// B-factor (nm2)
				fatomMag *= exp(-Biso * (Qmag/(4*PI))^2)	// exp[-B * (sin(theta)/lam)^2 ]
			elseif (itemp==3)
				Uiso = xtal.atom[m].Uiso
				fatomMag *= exp(-Qmag*Qmag*Uiso/2)				// B = (8*PI^2 U)
			elseif (itemp>=4)
				// qUq = q^t x U x q
				Variable qUq = (xtal.atom[m].U11)*qvec[0]^2 + (xtal.atom[m].U22)*qvec[1]^2 + (xtal.atom[m].U33)*qvec[2]^2
				if (itemp==5)
					qUq += 2*(xtal.atom[m].U23)*qvec[1]*qvec[2] + 2*(xtal.atom[m].U13)*qvec[2]*qvec[0] + 2*(xtal.atom[m].U12)*qvec[0]*qvec[1]
				endif
				DW = numtype(qUq) && qUq>0 ? 1 : exp(-qUq/2)
				fatomMag *= DW>0 ? DW : 1
			endif
		endif
		ifatomArg = cmplx(0,fatomArg)
		MatrixOP/FREE Fcm = fatomMag * sum(exp(c2PI*(ww x hkl) + ifatomArg))
		Fc += Fcm[0]										// accumulate for this atom
	endfor

	Variable Fr, Fi
#ifdef DO_HEXAGONAL_EXTRA
	if (system==HEXAGONAL || (usingHexAxes && system==TRIGONAL))
		Variable arg = 2*PI*((h+2*k)/3 + l*0.5)	// hexagonal has atoms at (0,0,0) and (1/3, 2/3, 1/2)
		Fr = 1. + cos(arg)
		Fi = sin(arg)
		Variable rr=real(Fc), ii=imag(Fc)
		Fc = cmplx(rr*Fr - ii*Fi, rr*Fi + ii*Fr)
		//  for hexagonal:
		//	h+2k=3n,		l=even;		F = 4*f			1
		//	h+2k=3n±1,	l=odd;		F = sqrt(3)*f   sqrt(3)/4
		//	h+2k=3n±1,	l=even;		F = f			1/4
		//	h+2k=3n,		l=odd; 		F = 0			0
	endif
#endif

	Fr = abs(real(Fc))<1e-8 ? 0 : real(Fc)		// set tiny values to zero
	Fi = abs(imag(Fc))<1e-8 ? 0 : imag(Fc)
	return cmplx(Fr,Fi)
End
//
Function/C Get_f_proto(AtomType,QAngstrom, keV, [valence])	// simple-minded fatom, just (Z-valence)
	string AtomType
	variable keV										// energy is ignored in this simple calculation
	Variable QAngstrom								// |Q| in (1/Angstrom), == 2*PI/d[Angstrom]
	variable valence									// optional integer for valence
	valence = ParamIsDefault(valence) || numtype(valence) ? 0 : valence

	Variable iz=LatticeSym#ZfromLabel(AtomType), Bval=0.6
	if (iz<1 || numtype(iz) || iz>102)
		return NaN
	endif

	Make/N=(102)/FREE Asum, Cvalues
	Asum[0]= {0.99735,1.98741,2.94052,5.76343,5.30249,5.75354,18.3887,7.69533,8.51545,9.89498,10.7273,11.7134,12.7675,14.6203,19.8767,25.3137,34.2962,17.6483,18.8237,19.9577,21.077,22.1475,23.0399,23.9394,24.5825}
	Asum[25]= {25.0759,25.4912,25.8733,26.3919,26.7177,27.4294,28.3022,29.2955,30.4209,31.889,31.5458,32.6372,33.2276,34.2758,35.3453,36.4462,37.5814,38.7151,39.8867,41.0858,42.315,43.5839,44.8458,46.1801,47.6055}
	Asum[50]= {49.2118,51.1507,53.7329,55.6953,62.8485,69.7236,66.8272,60.125,58.3443,56.8636,55.7896,55.1253,54.8513,54.7338,54.9527,55.3,55.7754,56.3494,57,57.6277,58.2686,58.9901,59.7891,60.6514,61.5686}
	Asum[75]= {62.5366,63.5314,64.5881,65.6985,66.8067,67.969,69.1621,70.3962,71.6794,73.0198,74.3032,75.8277,77.2375,78.732,79.5262,80.5945,81.6164,82.4483,83.4147,84.4334,85.2935,86.0964,86.9115,87.6827,88.4067}
	Asum[100]= {89.0798,89.5422}

	Cvalues[0]= {0.00125,0.01249,0.05988,-1.76339,-0.30409,0.24637,-11.3926,0.30413,0.48398,1.09912,1.26883,1.28489,1.23139,0.3787,-3.87732,-8.3143,-16.2972,1.35009,1.1743,1.03849,0.91762,0.84574,0.95226,1.05195}
	Cvalues[24]= {1.40818,1.91452,2.49899,3.11686,3.59838,4.26842,4.56068,4.69149,4.70026,4.57602,4.10864,5.43921,5.34841,5.76121,5.71886,5.65073,5.55047,5.41556,5.28192,5.11007,4.91427,4.68114,4.41158,4.14542}
	Cvalues[48]= {3.81227,3.38881,2.78462,1.84739,0.26635,-0.70709,-6.85854,-12.7404,-8.8456,-1.1393,1.64038,4.12018,6.19355,7.85731,9.12488,10.242,11.029,11.6817,12.2062,12.6322,12.9818,13.3489,13.7049,13.9824}
	Cvalues[72]= {14.1842,14.3239,14.4097,14.4449,14.4526,14.3992,14.2911,14.18,14.0203,13.8301,13.5991,13.3183,12.9796,12.6868,12.1642,11.7484,11.2497,11.4561,11.3864,11.3632,11.5358,11.5685,11.5431,11.6823}
	Cvalues[96]= {11.8853,12.0698,12.2983,12.5741,12.9008,13.4313}

	Variable Svector = QAngstrom/(4*PI)
	Variable f0 = Asum[iz-1] * exp(-Svector*Svector*Bval)
	f0 = numtype(f0) ? 1 : max(f0,0)			// always a valid number > 0
	if (valence)
		f0 += (f0/iz)^2 * valence
	endif
	return cmplx(f0,0)
End
//Function/C Get_f_proto(AtomType,QAngstrom, keV, [valence])	// simple-minded fatom, just (Z-valence)
//	string AtomType
//	variable keV										// energy is ignored in this simple calculation
//	Variable QAngstrom								// |Q| in (1/Angstrom), == 2*PI/d[Angstrom]
//	variable valence									// optional integer for valence
//	valence = ParamIsDefault(valence) || numtype(valence) ? 0 : valence
//
//	Variable iz=ZfromLabel(AtomType), Bval
//	if (iz<1 || numtype(iz))
//		return NaN
//	elseif (iz<=10)
//		Make/FREE BwTemp={58.3331,10.9071,4.33979,42.9165,23.0888,12.7188,0.02064,13.8964,11.2651, 9}	// use 9 for Ne
//		Bval = BwTemp[iz-1]
//	elseif (iz<=18)
//		Make/FREE coef={6.8736,0.011759,-0.025672}
//		Bval = poly(coef,iz-1)
//	elseif (iz<=50)
//		Make/FREE coef={51.647,-3.5557,0.085621,-0.00070461}
//		Bval = poly(coef,iz-1)
//	elseif (iz<=58)
//		Make/FREE BwTemp={5.24328,4.74225,4.27091,0.26422,0.23092,0.15152,0.1104,0.12335}
//		Bval = BwTemp[iz-51]
//	elseif (iz<=70)
//		Make/FREE coef={30.119,-0.85371,0.006597}
//		Bval = poly(coef,iz-1)
//	else
//		Make/FREE coef={-57.258,2.1161,-0.025226,9.8321e-05}
//		Bval = poly(coef,iz-1)
//	endif
//
//	Variable Svector = QAngstrom/(4*PI)
//	Variable f0 = iz * exp(-Svector*Svector*(Bval)) - valence
//	f0 = numtype(f0) ? 1 : max(f0,0)			// always a valid number > 0
//	return cmplx(f0,0)
//End
//Function/C Get_f_proto(AtomType,Qmag, keV, [valence])	// simple-minded fatom, just (Z-valence)
//	string AtomType
//	variable keV,Qmag
//	variable valence									// optional integer for valence
//	valence = ParamIsDefault(valence) ? 0 : valence
//	Variable freal = ZfromLabel(AtomType) - valence
//	freal = numtype(freal) ? 1 : max(freal,0)
//	return cmplx(freal,0)
//End
//
Static Function DW_factor_M(T,thetaM,Qmag,amu)		// calculates the M in exp(-M), no I/O
	Variable T			// Temperature (K)
	Variable thetaM	// Debye Temperature (K)
	Variable Qmag		// length of q vector (1/nm)
	Variable amu		// mass of atom (amu)

	if (T<0 || thetaM<0 || Qmag<0 || amu<0 || numtype(T+thetaM+Qmag+amu))
		return NaN
	elseif (T<=0)
		return 0
	endif

	Variable xx = thetaM / T
	if (xx>50000)
		return 0
	endif
	Variable Phi = Integrate1D(LatticeSym#PhiIntegrand,0,xx)/xx
	Variable B = (3*hbar^2 * T * c^2)/(2*(amu*amu_eV)*kB*thetaM^2)* (Phi + xx/4)
	Variable M = B *(Qmag*1.e9)^2		// Q is in 1/nm, but we need it in 1/m
	return M
End
//
Static Function PhiIntegrand(xx)		// this function should never be called with xx<0
	Variable xx
	return xx>0 ? xx/(exp(xx)-1) : 1
End


Function/C FstructMax(xtal,Qmag,[keV])
	STRUCT crystalStructure &xtal				// sruct defining the crystal, not changed here
	Variable Qmag										// Q (1/nm)
	Variable keV
	keV = ParamIsDefault(keV) || keV<=0 ? NumVarOrDefault("root:Packages:Lattices:keV",NaN) : keV

	String name
	Variable valence, m
	Variable/C fatomC, Fc=cmplx(0,0)				// the result, complex structure factor

	reMakeAtomXYZs(xtal)
	if (!(xtal.N>=1))
		return cmplx(1,0)
	endif
	FUNCREF Get_f_proto fa = $"Get_f"				// if Get_f() does not exist, then Get_f_proto() will be used

	for (m=0;m<xtal.N;m+=1)							// loop over the defined atoms
		name="root:Packages:Lattices:atom"+num2istr(m)
		valence = numtype(xtal.atom[m].valence) ? 0 : xtal.atom[m].valence
		Wave ww = $name
		if (valence==0 )
			fatomC = fa(Z2symbol(xtal.atom[m].Zatom),Qmag/10, keV)
		else
			fatomC = fa(Z2symbol(xtal.atom[m].Zatom),Qmag/10, keV,valence=valence)
		endif
		Fc += fatomC * cmplx(DimSize(ww,0) * xtal.atom[m].occ, 0)
	endfor
	return Fc
End


// to calculate mu, you need the CromerLiberman file
Function get_muOfXtal(keV, [printIt])					// returns mu of xtal (1/micron)
	Variable keV
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt
	if (strlen(WinList("CromerLiberman.ipf","","WIN:129"))<1)
		print "You must include the CromerLiberman.ipf to calculate mu"
		return NaN
	endif

	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))				//fill the lattice structure with test values
		DoAlert 0, "ERROR get_muOfXtal()\rno lattice structure found"
		return NaN
	endif

	if (numtype(keV) || keV<=0)
		keV = 10
		Prompt keV, "x-ray energy (keV)"
		DoPrompt "keV", keV
		if (V_flag)
			return NaN
		endif
		printIt = 1
	endif
	if (printIt)
		printf "get_muOfXtal(%g)\r", keV
	endif

	Variable mu = muOfXtal(xtal, keV)		// returns mu of xtal (1/micron)
	if (printIt)
		printf  "   mu['%s', %g keV] = %g (1/µm)  -->  absorption length %g (µm)\r",xtal.desc,keV,mu,1/mu
	endif
	return mu
End
//
Static Function muOfXtal(xtal, keV)	// returns mu of xtal (1/micron)
	STRUCT crystalStructure &xtal
	Variable keV

	FUNCREF Get_f_proto fa = $"Get_f"	// if Get_f() does not exist, then cannot calc mu
	if (numtype(keV) || keV<=0 || NumberByKey("ISPROTO",FuncRefInfo(fa)))
		return NaN								// return NaN if no Cromer or invalid energy
	endif
	// f" = sigma / (2 * re * lambda)
	// sigma = f" * (2 * re * lambda)
	// sigma = mu/n = mu*Vc		// 1/Vc = n, number density
	// mu = (2 * re * lambda * f") / Vc

	Variable/C fatomC
	Variable m, Fpp, mult, occ
	for (m=0,Fpp=0; m<xtal.N; m+=1)		// loop over the defined atoms, accumulate f"
		fatomC = fa(Z2symbol(xtal.atom[m].Zatom),0, keV,valence=(xtal.atom[m].valence))
		occ = xtal.atom[m].occ
		occ = occ <= 0 ? 1 : occ
		mult = xtal.atom[m].mult
		Wave watom = $("root:Packages:Lattices:atom"+num2istr(m))
		if (mult<1 && WaveExists(watom))
			mult = DimSize(watom,0)
		endif
		mult = max(mult,1)
		Fpp += imag(fatomC) * mult * occ
	endfor
	if (!(Fpp>=0))								// if Fpp is negative or NaN, then invalid
		return NaN
	endif
	Variable mu = 2*re_nm*(hc_keVnm/keV)*Fpp / (xtal.Vc)
	return (mu * 1000)						// convert mu from 1/nm --> 1/µm
End


Function setCromerEnergy(keV)
	Variable keV

	Prompt keV, "Energy (keV) for Cromer-Liberman"
	DoPrompt "Energy (keV)", keV
	if (V_flag)
		return keV
	endif
	Variable/G root:Packages:Lattices:keV = keV		// only used when calculating Cromer-Liberman values
	return keV
End


ThreadSafe Static Function Element_amu(Z)
	Variable Z
	String amuList
	amuList   = "1.00794;4.002602;6.941;9.012182;10.811;12.0107;14.0067;15.9994;18.9984032;20.1797;22.98977;24.305;"
	amuList += "26.981538;28.0855;30.973761;32.065;35.453;39.948;39.0983;40.078;44.95591;47.867;50.9415;51.9961;"
	amuList += "54.938049;55.845;58.9332;58.6934;63.546;65.409;69.723;72.64;74.9216;78.96;79.904;83.798;85.4678;87.62;"
	amuList += "88.90585;91.224;92.90638;95.94;98;101.07;102.9055;106.42;107.8682;112.411;114.818;118.71;121.76;"
	amuList += "127.6;126.90447;131.293;132.90545;137.327;138.9055;140.116;140.90765;144.24;145;150.36;151.964;157.25;"
	amuList += "158.92534;162.5;164.93032;167.259;168.93421;173.04;174.967;178.49;180.9479;183.84;186.207;190.23;"
	amuList += "192.217;195.078;196.96655;200.59;204.3833;207.2;208.98038;209;210;222;223;226;227;232.0381;231.03588;"
	amuList += "238.02891;237;244;243;247;247;251;252;257;258;259;262;261;262;266;264;277;268" 
	return (str2num(StringFromList(Z-1,amuList)))
End

Function reMakeAtomXYZs(xtal)
	STRUCT crystalStructure &xtal				// this sruct is filled  by this routine
	if (Exists("root:Packages:Lattices:atom0")==1)
		Wave ww = root:Packages:Lattices:atom0
		if (WaveExists(ww) && strlen(xtal.hashID))	// if first atom has correct hash, assume the rest are OK too
			if (stringmatch(StringByKey("ID",note(ww),"="),xtal.hashID))
				return 0									// waves have correct hash, so return
			endif
		endif
	endif
	if (strlen(xtal.hashID)<1)
		xtal.hashID = xtalHashID(xtal)			// re-set hash function to identify associated waves
	endif
	String wnote=ReplaceStringByKey("ID","",xtal.hashID,"="), name
	Variable m
	for (m=0;m<xtal.N;m+=1)						// loop over each atom type
		name = "root:Packages:Lattices:atom"+num2istr(m)
		Make/N=3/O/D $name
		Wave ww = $name
		positionsOfOneAtomType(xtal,xtal.atom[m].x,xtal.atom[m].y,xtal.atom[m].z,ww)
		wnote = ReplaceStringByKey("atomType",wnote,xtal.atom[m].name,"=")
		wnote = ReplaceNumberByKey("Zatom",wnote,xtal.atom[m].Zatom,"=")
		wnote = ReplaceNumberByKey("occupy",wnote,xtal.atom[m].occ,"=")
		wnote = ReplaceNumberByKey("valence",wnote,xtal.atom[m].valence,"=")
		Note/K ww, wnote
	endfor
	return 0
End
//Function testFstruct(h,k,l)
//	Variable h,k,l
//	STRUCT crystalStructure xtal				// this sruct is filled  by this routine
////	readCrystalStructure(xtal,"Macintosh HD:Users:tischler:Desktop:new sym calc:al2o3.cri")
////	readCrystalStructure(xtal,"Macintosh HD:Users:tischler:Desktop:new sym calc:si.cri")
////	readCrystalStructure(xtal,"")
////	readCrystalStructure(xtal,"silicon.xtl")
//	FillCrystalStructDefault(xtal)
//
//	Variable/C Fc = Fstruct(xtal,h,k,l)
//	printf "F(%d %d %d) = %g %s i%g,      |F| = %g\r",h,k,l,real(Fc),SelectString(imag(Fc)<0,"+","-"),abs(imag(Fc)),sqrt(magsqr(Fc))
//End
//
Static Function positionsOfOneAtomType(xtal,xx,yy,zz,xyzIN)
	STRUCT crystalStructure &xtal	// provides SpaceGroup, a,b,c
	Variable xx,yy,zz		// fractional coords of this kind of atom
	Wave xyzIN				// result, list of all equiv positions for this atom in fractional coords

	String SpaceGroupID=xtal.SpaceGroupID
	SetSymOpsForSpaceGroup(SpaceGroupID)		// ensure existance of symmetry op mats and vecs
	Wave mats = $("root:Packages:Lattices:SymOps:"+CleanupName("equivXYZM"+xtal.SpaceGroupID,0))
	Wave bvecs = $("root:Packages:Lattices:SymOps:"+CleanupName("equivXYZB"+xtal.SpaceGroupID,0))
	if (!WaveExists(mats) || !WaveExists(bvecs))
		Abort"Unable to get symmetry operations in positionsOfOneAtomType()"
	endif
	Wave direct = directFrom_xtal(xtal)		// get direct lattice from xtal
	Variable minDist2 = minPossibleBondLength^2

	Make/N=(3,3)/D/FREE mat
	Make/N=3/D/FREE bv, in={xx,yy,zz}, vec
	Variable m,Neq=NumberByKey("numSymOps", note(mats),"=")

	Make/N=(Neq,3)/D/FREE xyz=NaN				// internal copy of fractional coords
	Make/N=(Neq,3)/D/FREE xyznm=NaN				// real positions (in nm), NOT fractional coords (in sync with xyz[][])
	Variable N
	// printf "atom at  %.5f, %.5f, %.5f\r",in[0],in[1],in[2]
	for (m=0,N=0;m<Neq;m+=1)
		mat = mats[m][p][q]
		bv = bvecs[m][p]
		MatrixOp/FREE rr = mat x in + bv		// rr is relative coord of (xx,yy,zz) after operation
		MatrixOp/FREE rr = rr - floor(rr)		// reduce to [0,1), the first unit cell
		rr = mod(rr,1)
		rr = rr<1e-12 ? 0 : rr						// a fractional coord < 1e-12 is 0

		MatrixOP/FREE vec = direct x rr			// real space vector for rr
		if (Neq<2)
			MatrixOP/FREE isDup = maxVal( greater(minDist2, sumRows(magSqr(xyznm - vec^t))) )
		else
			MatrixOP/FREE isDup = maxVal( greater(minDist2, sumRows(magSqr(xyznm - rowRepeat(vec,Neq)))) )
		endif
		if (isDup[0]<1)								// not a duplicate, so add to the list of positions
			xyz[N][] = rr[q]							// rr is not an equivalent atom, save it to xyz[N]
			xyznm[N][] = vec[q]						// keep xyznm in sync with xyz
			N += 1
		endif
	endfor

	Wave Unconventional=root:Packages:Lattices:Unconventional
	if (WaveExists(Unconventional))					// Unconventional exists, transform all the fractional coords
		MatrixOP/FREE xyz = ( Unconventional x (xyz^t) )^t
	endif

	Redimension/N=(N,3) xyzIN						// remove extra space (it is all filled with NaN's)
	xyzIN = xyz[p][q]									// and, update xyzIN with correct values
	return N
End


//Function allowedHKL(h,k,l,xtal)					// NOT threadsafe, but this could use Cromer
//	Variable h,k,l
//	STRUCT crystalStructure &xtal
//	Variable/C Fc = Fstruct(xtal,h,k,l)
//	return (magsqr(Fc)/(xtal.N)^2 > 0.0001)	// allowed means more than 0.01 electron/atom
//End
//
// Example:
//	Make/N=(xtal.N)/WAVE/FREE atomWaves=$("root:Packages:Lattices:atom"+num2istr(p))	// cannot reach these waves from within a separate thread
//	Make/N=(Nz)/FREE/WAVE testAllowed
//	MultiThread testAllowed = allowedHKL(h[p],k[p],l[p],xtal, atomWaves=atomWaves)
//
ThreadSafe Function allowedHKL(h,k,l,xtal,[atomWaves])		// does NOT use Cromer, but can be multi-threaded
	Variable h,k,l											// the hkl may be non-integers
	STRUCT crystalStructure &xtal
	Wave/WAVE atomWaves									// ONLY needed when calling this with MultiThread

	if (!isValidSpaceGroupID(xtal.SpaceGroupID) || numtype(h+k+l))
		return 0												// bad inputs (not allowed)
	endif
	Variable usingHexAxes = (abs(90-xtal.alpha)+abs(90-xtal.beta)+abs(120-xtal.gam))<1e-6
	Variable system = latticeSystem(xtal.SpaceGroupID)

	if (!(mod(h,1) || mod(k,1) || mod(l,1)))	// non-integral, always allowed
//		String sym = getHMsym(xtal.SpaceGroupIDnum)	// get symmetry symbol
		String sym = getHMsym2(xtal.SpaceGroupIDnum)	// get symmetry symbol
		strswitch (sym[0,0])
			case "F":
				if (mod(h+k,2) || mod(k+l,2) )		// face-centered, hkl must be all even or all odd
					return 0
				endif
				break
			case "I":
				if (mod(round(h+k+l),2))				// body-centered, !mod(round(h+k+l),2), sum must be even
					return 0
				endif
				break
			case "C":
				if (mod(round(h+k),2))					// C-centered, !mod(round(h+k),2)
					return 0
				endif
				break
			case "A":
				if (mod(round(k+l),2))					// A-centered, !mod(round(k+l),2)
					return 0
				endif
				break
			case "R":
				if (usingHexAxes)							// rhombohedral cell with hexagonal axes
					if (mod(-h+k+l,3) && mod(h-k+l,3))	//		allowed are -H+K+L=3n or H-K+L=3n
						return 0
					endif
				endif
				break											// using rhombohedral axes, so all are allowed
		endswitch
		if (system==HEXAGONAL)
			if (!mod(h+2*k,3) && mod(l,2)) 			// hexagonal, forbidden are: H+2K=3N with L odd
				return 0
			endif
		endif
	endif

//	reMakeAtomXYZs(xtal)								// NOT ThreadSafe
	if (xtal.N<1)
		return 1												// No atom defined, but passed simple tests, it is allowed
	endif

	Variable fatomMag, m, FMax=0
	Variable/C c2PI=cmplx(0,2*PI)
	Variable/C Fc=cmplx(0,0)							// the result, complex structure factor
	Make/N=3/D/FREE hkl={h,k,l}
	for (m=0;m<xtal.N;m+=1)							// loop over the defined atoms
		if (WaveExists(atomWaves))
			Wave ww = atomWaves[m]						// needed when running in a separate thread
		else
			Wave ww = $("root:Packages:Lattices:atom"+num2istr(m))	// not available from a separate thread
		endif
		fatomMag = xtal.atom[m].Zatom * xtal.atom[m].occ	// just use Z for the f_atom (this is always REAL)
		if (WaveExists(ww))
			MatrixOP/O/FREE Fcm = fatomMag * sum(exp(c2PI*(ww x hkl)))
			Fc += Fcm[0]									// accumulate for this atom
			FMax += max(fatomMag,0.01) * DimSize(ww,0)
		else
			Fc += cmplx(fatomMag,0)					// no atom position, just make it in phase
			FMax += max(fatomMag,0.01)				// always at least 0.01 electrons/atom
		endif
	endfor

#ifdef DO_HEXAGONAL_EXTRA
	if (system==HEXAGONAL || (usingHexAxes && system==TRIGONAL))
		Variable arg, Fr, Fi
		arg = 2*PI*((h+2*k)/3 + l*0.5)				// hexagonal has atoms at (0,0,0) and (1/3, 2/3, 1/2)
		Fr = 1. + cos(arg)
		Fi = sin(arg)
		Variable rr=real(Fc), ii=imag(Fc)
		Fc = cmplx(rr*Fr - ii*Fi, rr*Fi + ii*Fr)
		//  for hexagonal:
		//	h+2k=3n,		l=even;		F = 4*f			1
		//	h+2k=3n±1,	l=odd;		F = sqrt(3)*f   sqrt(3)/4
		//	h+2k=3n±1,	l=even;		F = f			1/4
		//	h+2k=3n,		l=odd; 		F = 0			0
	endif
#endif

	return (magsqr(Fc)/(xtal.N)^2 > (FMax/50)^2)			// allowed means more than 0.01 electron/atom
//	return (magsqr(Fc)/(xtal.N)^2 > 0.0001)			// allowed means more than 0.01 electron/atom
End


ThreadSafe Function/WAVE equivalentHKLs(xtal,hkl0,[noNeg])
	// returns a wave of all the hkls that are symmetry equivalent to hkl0
	// e.g. hl0=(100) returns { (100), (010), (001), (-100), (0 -1 00, (0 0 -1) }, 6 hkl's
	// and with noNeg=1, returns { (100), (010), (0 0 -1) }, 3 hkl's

	STRUCT crystalStructure &xtal			// convert the hkls to q vectors
	Wave hkl0										// given hkl
	Variable noNeg									// flag, if true, then do not include (-1,-1,-1) for (111)
	noNeg = ParamIsDefault(noNeg) || numtype(noNeg) ? 0 : !(!noNeg)

	if (numtype(sum(hkl0))!=0)				// bad values
		return $""
	elseif (norm(hkl0)<=0)						// hkl0 == 0
		Make/N=3/D/FREE hkls={0,0,0}
		return hkls
	endif
	Variable thresh = norm(hkl0)/1e4		// error that triggers a difference

	Wave SymmetryOp = $"root:Packages:Lattices:SymOps:SymmetryOps"+num2istr(xtal.SpaceGroup)
	if (!WaveExists(SymmetryOp))
		Wave SymmetryOp = $MakeSymmetryOps(xtal)	// wave with all the symmetry operation
	endif

	Variable Nproper=NumberByKey("Nproper",note(SymmetryOp),"=")
	if (Nproper<1)
		return $""									// there should always be at least 1 proper rotation (identity)
	endif
	Make/N=(3,3)/D/FREE rot
	Make/N=(Nproper,3)/D/FREE hklEquiv=NaN	// holds all the symmetry equivalent hkl

	if (Nproper==1)
		hklEquiv[0][] = hkl0[q]				// only 1 proper rotation, must be the identity
		return hklEquiv							//   this also avoids error at RowRepeat(hkl,1)
	endif

	Variable isym,NsymOps, diffMin
	for (isym=0, NsymOps=0; isym<Nproper; isym+=1)
		rot = SymmetryOp[isym][p][q]
		MatrixOp/FREE/O hkl = rot x hkl0	// possible hkl to add to hklEquiv
		MatrixOp/FREE/O diffs = sumRows(magsqr(hklEquiv - RowRepeat(hkl,Nproper)))
		diffMin = WaveMin(diffs)
		if (noNeg)									// try again with -hkl
			MatrixOp/FREE/O diffs = sumRows(magsqr(hklEquiv - RowRepeat(-hkl,Nproper)))
			diffMin = min(diffMin,WaveMin(diffs))
		endif
		if (diffMin>thresh || NsymOps<1)	// no match in hklEquiv, add this hkl to list
			hklEquiv[NsymOps][] = hkl[q]
			NsymOps += 1
		endif
	endfor

	Redimension/N=(NsymOps,3) hklEquiv
	return hklEquiv
End
//	Function test_equivalentHKLs(hklStr)
//		String hklStr
//		Wave hkl0 = str2vec(hklStr)
//	
//		STRUCT crystalStructure xtal
//		if (FillCrystalStructDefault(xtal))				//fill the lattice structure with test values
//			DoAlert 0, "ERROR RadialCorrelationFrom3D()\rno lattice structure found"
//			return 1
//		endif
//	//	normalize(hkl0)
//		Wave qhats = equivalentHKLs(xtal,hkl0,noNeg=1)// list of desired directions dims of qhats=(N,3)
//		Variable Nqs = DimSize(qhats,0)
//		print "from hkl0 = ",vec2str(hkl0),"   Nqs = ",Nqs
//		printWave(qhats,name="",brief=1)
//	End


ThreadSafe Function/S MakeSymmetryOps(xtal)			// make a wave with the symmetry operation
	STRUCT crystalStructure &xtal

	String SG = num2istr(xtal.SpaceGroup), wName
	wName = ReplaceString(":",xtal.SpaceGroupID,"_")	// cannot use CleanUpName() here, not ThreadSafe
	wName = "root:Packages:Lattices:SymOps:equivXYZM"+ReplaceString("-",wName,"_")
	Wave equivXYZM = $wName
	if (!WaveExists(equivXYZM))
		return ""
	endif
	Variable Nequiv=DimSize(equivXYZM,0)					// number of all symmetry operations for this space group
	String/G root:Packages:Lattices:SymOps:SymmetryOpsPath="root:Packages:Lattices:SymOps:SymmetryOps"+SG

	wName = "root:Packages:Lattices:SymOps:SymmetryOps"+SG
	Make/N=(2*Nequiv,3,3)/D/O $wName=NaN
	Wave ops = $wName

	Make/N=(3,3)/D/FREE direct,mat
	direct = { {xtal.a0, xtal.a1, xtal.a2}, {xtal.b0, xtal.b1, xtal.b2}, {xtal.c0, xtal.c1, xtal.c2} }
	MatrixOp/FREE/O directI = Inv(direct)

	Variable i, N=0, Nproper=0
	for (i=0;i<Nequiv;i+=1)									// loop through all Nequiv, rejecting duplicates, and only accepting proper rotations
		mat = equivXYZM[i][p][q]
		MatrixOp/O/FREE mat = direct x mat x directI	// convert to cartesian, similarity transform
		if (isMatInMats(mat,ops) ||MatrixDet(mat)<0)	// skip duplicates and improper rotations
			continue
		endif
		ops[N][][] = mat[q][r]									// save this mat
		Nproper += MatrixDet(mat)>0 ? 1 : 0
		N += 1
	endfor

	// go through list again, this time only taking unique IMproper rotations
	for (i=0;i<Nequiv;i+=1)									// again, loop through all Nequiv, rejecting duplicates
		mat = equivXYZM[i][p][q]
		MatrixOp/O/FREE mat = direct x mat x directI	// convert to cartesian, similarity transform
		if (isMatInMats(mat,ops))								// skip duplicates
			continue
		endif
		ops[N][][] = mat[q][r]									// save this mat
		N += 1
	endfor
	Redimension/N=(N,-1,-1) ops

	ops = abs(ops[p][q][r])<1e-13 ? 0 : ops[p][q][r]		// remove the almost zeros
	ops = abs(1-ops[p][q][r])<1e-13 ? 1 : ops[p][q][r]	//  and make almost 1's equal to 1
	ops = abs(1+ops[p][q][r])<1e-13 ? -1 : ops[p][q][r]

	String wnote="waveClass=SymmetryOperations;"
	wnote = ReplaceNumberByKey("Nproper",wnote,Nproper,"=")
	wnote = ReplaceNumberByKey("SpaceGroup",wnote,xtal.SpaceGroup,"=")
	wnote = ReplaceStringByKey("SpaceGroupID",wnote,xtal.SpaceGroupID,"=")
	Note/K ops, wnote
	return GetWavesDataFolder(ops,2)
End
//
ThreadSafe Static Function isMatInMats(mat,ops)
	Wave mat, ops

	Make/N=(3,3)/FREE/D delta
	Variable i,N=DimSize(ops,0)							// Make/N=(2*N,3,3)/D/O ops=NaN
	for (i=0;i<N;i+=1)
		delta = mat[p][q] - ops[i][p][q]
		if (numtype(sum(delta)))							// extra mats in ops are filled with NaN
			break
		endif
		delta = abs(delta)
		if (sum(delta)<1e-13)
			return 1
		endif
	endfor
	return 0
End
//
//Function testAllSymmetryOps()
//	Make/N=(230,3)/O NtestAll=NaN
//	SetScale/P x 1,1,"", NtestAll
//
//	Variable SpaceGroup, err
//	for (SpaceGroup=1; SpaceGroup<=230 && !err; SpaceGroup+=1)
//		LatticeSym#SetSymOpsForSpaceGroup(SpaceGroup)
//		err = testOneSymmetryOp(SpaceGroup)
//		Wave Mwave = $("root:Packages:Lattices:SymOps:equivXYZM"+num2istr(SpaceGroup))
//		NtestAll[SpaceGroup-1][0] = DimSize(Mwave,0)
//		Wave Swave = $("root:Packages:Lattices:SymOps:SymmetryOps"+num2istr(SpaceGroup))
//		NtestAll[SpaceGroup-1][1] = DimSize(Swave,0)
//		NtestAll[SpaceGroup-1][21] = NumberByKey("Nproper",note(Swave),"=")
//	endfor
//	String str=SelectString(err,"checked all 230 Space Groups with no error","ERROR on Space Group "+num2istr(SpaceGroup-1))
//	print str
//	DoAlert 0, str
//	if (!err)
//		KillWaves/Z detTest, NtestAll
//	endif
//End
//Static Function testOneSymmetryOp(SpaceGroup)
//	Variable SpaceGroup
//	if (!isValidSpaceGroup(SpaceGroup))
//		printf "SpaceGroup = %g, not in [1,230]\r",SpaceGroup
//		return 1
//	endif
//
//	STRUCT crystalStructure xtal
//	FillCrystalStructDefault(xtal)
//	xtal.SpaceGroup = SpaceGroup
//	LatticeSym#setDirectRecip(xtal)
//	initSymmetryOperations()						// initialize all symmetry operations
//	indexLots#MakeSymmetryOps(xtal)				// make a wave with the symmetry operation
//	SVAR symName=root:Packages:Lattices:SymOps:SymmetryOpsPath
//	Wave sym=$symName
//
//	Variable Nproper = NumberByKey("Nproper",note(sym),"="), N=DimSize(sym,0), err=0
//	if (N!=2*Nproper && 0)
//		printf "SpaceGroup = %d,  N=%g, and Nproper = %g,   Nproper is not N/2\r",SpaceGroup,N,Nproper
//		err = 1
//	endif
//	Make/N=(N)/O detTest=NaN
//	Make/N=(3,3)/D/FREE mat
//	Variable i,det
//	for (i=0;i<N;i+=1)
//		mat = sym[i][p][q]
//		det = MatrixDet(mat)
//		detTest[i] = det
//		if (i<Nproper && abs(det-1)>1e-8)
//			printf "SpaceGroup = %d,  for i=%d, should be proper rotation, but det = %g\r",SpaceGroup,i,det
//			err = 1
//		elseif (i>=Nproper && abs(det+1)>1e-8)
//			printf "SpaceGroup = %d,  for i=%d, should be IMproper rotation, but det = %g\r",SpaceGroup,i,det
//			err = 1
//		elseif (abs(abs(det)-1)>1e-8)
//			printf "SpaceGroup = %d,  for i=%d, strange error, det = %g\r",SpaceGroup,i,det
//			err = 1
//		endif
//	endfor
//	return err
//End


Function/T DescribeSymOps(SymOps,[printIt])		// prints description of symmetry operations, returns result as a list too
	Wave SymOps
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : !(!printIt)

	if (!WaveExists(SymOps))
		String SymmetryOpsPath=StrVarOrDefault("root:Packages:Lattices:SymOps:SymmetryOpsPath",""), id=""
		Wave SymOps = $StrVarOrDefault("root:Packages:Lattices:SymOps:SymmetryOpsPath","")
		Variable SGw = 0
		if (WaveExists(SymOps))
			id = StringByKey("SpaceGroupID",note(SymOps),"=")
			SGw = str2num(id)
			SGw = numtype(SGw) ? NumberByKey("SpaceGroup",note(SymOps),"=") : SGw
		endif
		SGw = numtype(SGw) || SGw<1 ? 0 : SGw
		STRUCT crystalStructure xtal
		if (FillCrystalStructDefault(xtal) && !WaveExists(SymOps))
			return ""
		endif

		Variable useWave=0, useXtal=0
		if (isValidSpaceGroupID(id) && isValidSpaceGroupID(xtal.SpaceGroupID) && StringMatch(xtal.SpaceGroupID,id))
			useWave = 1
		elseif (isValidSpaceGroupID(id) && !isValidSpaceGroupID(xtal.SpaceGroupID))		// use wave
			useWave = 1
		elseif (!isValidSpaceGroupID(id) && isValidSpaceGroupID(xtal.SpaceGroupID))		// use xtal
			useXtal = 1
		elseif (isValidSpaceGroupID(id) && isValidSpaceGroupID(xtal.SpaceGroupID))		// ask
			DoAlert 2, "SymmetryOpsPath and Lattice Panel Disagree,\r  Remake SymmetryOps"+num2istr(SGw)+" to match Panel?"
			if (V_flag>2)
				return ""
			endif
			useWave = (V_flag==2)
			useXtal = (V_flag==1)
		endif
		if (useXtal)
			String wname = MakeSymmetryOps(xtal)
			String/G root:Packages:Lattices:SymOps:SymmetryOpsPath = wname
			Wave SymOps = $wname
		endif
	endif
	if (!WaveExists(SymOps))
		return ""
	elseif (DimSize(SymOps,1)!=3 || DimSize(SymOps,2)!=3)
		return ""
	endif
	Variable Nproper = NumberByKey("Nproper",note(SymOps),"=")
	String SpaceGroupID = StringByKey("SpaceGroupID",note(SymOps),"=")
	if (!(Nproper>0))
		return ""
	endif
	if (printIt)
		String system = StringFromList(latticeSystem(SpaceGroupID),LatticeSystemNames)
		Variable idNum = SpaceGroupID2num(SpaceGroupID)
		printf "For Space Group %s  (%s) %s,  %g proper rotations\r",SpaceGroupID,system,getHMsym2(idNum),Nproper
	endif

	Make/N=3/D/FREE axis, v3
	Make/N=(3,3)/D/FREE sym
	String out = ReplaceNumberByKey("Nproper","",Nproper,"=")
	out = ReplaceNumberByKey("SpaceGroup",out,str2num(SpaceGroupID),"=")
	out = ReplaceStringByKey("SpaceGroupID",out,SpaceGroupID,"=")
	String str, name								// name of axis
	Variable i, angle, div
	for (i=0;i<Nproper;i+=1)
		sym = SymOps[i][p][q]
		angle = axisOfMatrix(sym,axis)		// returns normalized axis
		angle = abs(angle)<0.1 ? 0 : angle

		if (abs(axis[0]-1)<0.02)
			name = "X-axis"
		elseif (abs(axis[1]-1)<0.02)
			name = "Y-axis"
		elseif (abs(axis[2]-1)<0.02)
			name = "Z-axis"
		elseif (abs(axis[0]+1)<0.02)
			angle = -angle
			name = "X-axis"
		elseif (abs(axis[1]+1)<0.02)
			angle = -angle
			name = "Y-axis"
		elseif (abs(axis[2]+1)<0.02)
			angle = -angle
			name = "Z-axis"
		elseif (abs(angle)>0)
			div = smallestNonZeroValue(axis)
			axis /= div
			name = vec2str(axis)+" axis"
		endif
		if (abs(angle)<0.1)
			str = "Identity (no rotation)"
		else
			sprintf str, "% 4.0f%s rotation about the %s",angle,DEGREESIGN,name
		endif
		out += str+";"
		if (printIt)
			print str
		endif
	endfor
	return out
End


ThreadSafe Static Function isValidSpaceGroup(SG)			// returns TRUE if SG is an int in range [1,230]
	Variable SG
	return ( SG == limit(round(SG), 1, 230) )
End


ThreadSafe Static Function isValidSpaceGroupID(id)		// returns TRUE if id is valid
	String id															// a space group id, e.g. "15" or "15:-b2"
	String allIDs=MakeAllIDs()
	return WhichListItem(id,allIDs,";",0,0)>=0
End


ThreadSafe Static Function isValidSpaceGroupIDnum(idNum)	// returns TRUE if SG is an int in range [1,530]
	Variable idNum
	return ( idNum == limit(round(idNum), 1, 530) )
End


Function SpaceGroupID2num(id)
	String id									// a space group id, e.g. "15" or "15:-b2"

	String allIDs=MakeAllIDs()
	Variable idNum = 1+WhichListItem(id,allIDs,";",0,0)
	idNum = isValidSpaceGroupIDnum(idNum) ? idNum : NaN
	if (numtype(idNum))
		idNum = FindDefaultIDnumForSG(round(str2num(id)))
	endif
	idNum = isValidSpaceGroupIDnum(idNum) ? idNum : NaN
	return idNum
End


ThreadSafe Static Function/T FindDefaultIDforSG(SG)
	Variable SG					// space group number [1,230]
	if (!isValidSpaceGroup(SG))		// invalid
		return ""
	endif

	// find first space group starting with "id:"
	string str=num2istr(SG)+":*", str2=num2istr(SG)
	String allIDs=MakeAllIDs(), id
	Variable i
	for (i=0;i<ItemsInList(allIDs);i+=1)
		id = StringFromList(i,allIDs)
		if (StringMatch(StringFromList(i,allIDs),str) || !cmpstr(id,str2))
			return StringFromList(i,allIDs)
		endif	
	endfor
	return ""
End


ThreadSafe Static Function FindDefaultIDnumForSG(SG)
	// look for exatly "SG", or the first occurance of "SG:"
	Variable SG					// space group number [1,230]
	if (!isValidSpaceGroup(SG))		// invalid
		return NaN
	endif

	// find first space group starting with "id:"
	string str=num2istr(SG)+":*", str2=num2istr(SG)
	String allIDs=MakeAllIDs(), id
	Variable i
	for (i=0;i<ItemsInList(allIDs);i+=1)
		id = StringFromList(i,allIDs)
		if (StringMatch(id,str) || !cmpstr(id,str2))
			return i+1
		endif	
	endfor
	return NaN
End




ThreadSafe Static Function/T MakeAllIDs()
	// Returns the list with all of the 530 Space Group types.
	//	In the 230 SpaceGroups, there are:
	//	  140 Space Groups of   1 types
	//	   30 Space Groups of   2 types
	//	   26 Space Groups of   3 types
	//	   25 Space Groups of   6 types
	//	    6 Space Groups of   9 types
	//	    1 Space Groups of  12 types
	//	    2 Space Groups of  18 types
	// for the full list, use  NumbersOfTypes(), which is shown below.

	String allIDs = "1;2;3:b;3:c;3:a;4:b;4:c;4:a;5:b1;5:b2;5:b3;5:c1;5:c2;5:c3;5:a1;5:a2;"
	allIDs += "5:a3;6:b;6:c;6:a;7:b1;7:b2;7:b3;7:c1;7:c2;7:c3;7:a1;7:a2;7:a3;8:b1;8:b2;8:b3;"
	allIDs += "8:c1;8:c2;8:c3;8:a1;8:a2;8:a3;9:b1;9:b2;9:b3;9:-b1;9:-b2;9:-b3;9:c1;9:c2;"
	allIDs += "9:c3;9:-c1;9:-c2;9:-c3;9:a1;9:a2;9:a3;9:-a1;9:-a2;9:-a3;10:b;10:c;10:a;11:b;"
	allIDs += "11:c;11:a;12:b1;12:b2;12:b3;12:c1;12:c2;12:c3;12:a1;12:a2;12:a3;13:b1;13:b2;"
	allIDs += "13:b3;13:c1;13:c2;13:c3;13:a1;13:a2;13:a3;14:b1;14:b2;14:b3;14:c1;14:c2;14:c3;"
	allIDs += "14:a1;14:a2;14:a3;15:b1;15:b2;15:b3;15:-b1;15:-b2;15:-b3;15:c1;15:c2;15:c3;"
	allIDs += "15:-c1;15:-c2;15:-c3;15:a1;15:a2;15:a3;15:-a1;15:-a2;15:-a3;16;17;17:cab;17:bca;"
	allIDs += "18;18:cab;18:bca;19;20;20:cab;20:bca;21;21:cab;21:bca;22;23;24;25;25:cab;25:bca;"
	allIDs += "26;26:ba-c;26:cab;26:-cba;26:bca;26:a-cb;27;27:cab;27:bca;28;28:ba-c;28:cab;"
	allIDs += "28:-cba;28:bca;28:a-cb;29;29:ba-c;29:cab;29:-cba;29:bca;29:a-cb;30;30:ba-c;"
	allIDs += "30:cab;30:-cba;30:bca;30:a-cb;31;31:ba-c;31:cab;31:-cba;31:bca;31:a-cb;32;"
	allIDs += "32:cab;32:bca;33;33:ba-c;33:cab;33:-cba;33:bca;33:a-cb;34;34:cab;34:bca;35;"
	allIDs += "35:cab;35:bca;36;36:ba-c;36:cab;36:-cba;36:bca;36:a-cb;37;37:cab;37:bca;38;"
	allIDs += "38:ba-c;38:cab;38:-cba;38:bca;38:a-cb;39;39:ba-c;39:cab;39:-cba;39:bca;39:a-cb;"
	allIDs += "40;40:ba-c;40:cab;40:-cba;40:bca;40:a-cb;41;41:ba-c;41:cab;41:-cba;41:bca;"
	allIDs += "41:a-cb;42;42:cab;42:bca;43;43:cab;43:bca;44;44:cab;44:bca;45;45:cab;45:bca;"
	allIDs += "46;46:ba-c;46:cab;46:-cba;46:bca;46:a-cb;47;48:1;48:2;49;49:cab;49:bca;50:1;"
	allIDs += "50:2;50:1cab;50:2cab;50:1bca;50:2bca;51;51:ba-c;51:cab;51:-cba;51:bca;51:a-cb;"
	allIDs += "52;52:ba-c;52:cab;52:-cba;52:bca;52:a-cb;53;53:ba-c;53:cab;53:-cba;53:bca;"
	allIDs += "53:a-cb;54;54:ba-c;54:cab;54:-cba;54:bca;54:a-cb;55;55:cab;55:bca;56;56:cab;"
	allIDs += "56:bca;57;57:ba-c;57:cab;57:-cba;57:bca;57:a-cb;58;58:cab;58:bca;59:1;59:2;"
	allIDs += "59:1cab;59:2cab;59:1bca;59:2bca;60;60:ba-c;60:cab;60:-cba;60:bca;60:a-cb;61;"
	allIDs += "61:ba-c;62;62:ba-c;62:cab;62:-cba;62:bca;62:a-cb;63;63:ba-c;63:cab;63:-cba;"
	allIDs += "63:bca;63:a-cb;64;64:ba-c;64:cab;64:-cba;64:bca;64:a-cb;65;65:cab;65:bca;66;"
	allIDs += "66:cab;66:bca;67;67:ba-c;67:cab;67:-cba;67:bca;67:a-cb;68:1;68:2;68:1ba-c;"
	allIDs += "68:2ba-c;68:1cab;68:2cab;68:1-cba;68:2-cba;68:1bca;68:2bca;68:1a-cb;68:2a-cb;"
	allIDs += "69;70:1;70:2;71;72;72:cab;72:bca;73;73:ba-c;74;74:ba-c;74:cab;74:-cba;74:bca;"
	allIDs += "74:a-cb;75;76;77;78;79;80;81;82;83;84;85:1;85:2;86:1;86:2;87;88:1;88:2;89;90;"
	allIDs += "91;92;93;94;95;96;97;98;99;100;101;102;103;104;105;106;107;108;109;110;111;112;"
	allIDs += "113;114;115;116;117;118;119;120;121;122;123;124;125:1;125:2;126:1;126:2;127;"
	allIDs += "128;129:1;129:2;130:1;130:2;131;132;133:1;133:2;134:1;134:2;135;136;137:1;137:2;"
	allIDs += "138:1;138:2;139;140;141:1;141:2;142:1;142:2;143;144;145;146:H;146:R;147;148:H;"
	allIDs += "148:R;149;150;151;152;153;154;155:H;155:R;156;157;158;159;160:H;160:R;161:H;"
	allIDs += "161:R;162;163;164;165;166:H;166:R;167:H;167:R;168;169;170;171;172;173;174;175;"
	allIDs += "176;177;178;179;180;181;182;183;184;185;186;187;188;189;190;191;192;193;194;"
	allIDs += "195;196;197;198;199;200;201:1;201:2;202;203:1;203:2;204;205;206;207;208;209;"
	allIDs += "210;211;212;213;214;215;216;217;218;219;220;221;222:1;222:2;223;224:1;224:2;"
	allIDs += "225;226;227:1;227:2;228:1;228:2;229;230"
	return allIDs
End
//
//	Function NumbersOfTypes()
//		String allIDs = LatticeSym#MakeAllIDs()
//		Make/N=230/I/O types=0
//		SetScale/P x 1,1,"SpaceGroup", types
//		Variable i,m
//		for (i=0;i<530;i+=1)
//			m = str2num(StringFromList(i,allIDs))
//			types[m-1] +=1
//		endfor
//	
//		Variable tMax=WaveMax(types)
//		Duplicate/FREE types cntTypes
//		Make/N=(tMax)/I/FREE sums
//		print "In the 230 SpaceGroups, there are:"
//		for (i=0;i<tMax;i+=1)
//			cntTypes = (types[p]==(i+1))
//			if (sum(cntTypes))
//				printf "  %3d Space Groups with %3d types\r",sum(cntTypes), i+1
//			endif
//		endfor
//
//		Display /W=(298,120,1440,504) types
//		ModifyGraph mode=4, marker=19, lStyle=1, tick=2, mirror=1, minor=1, grid=1
//		Label left "number of types"
//		SetAxis/A/N=1/E=1 left
//	End



// returns info about the symmetry of a structure
// sym holds the symmetry info on return, and xyz holds the atom positions.  It returns the number of atom positions put in xyz
// which is always at least 1.  If you call with a bad wave ref for xyz, then it only returns the sym string, and the returned value is 0
//
//	Extensions
//	----------
//		Monoclinic			unique axis b		unique axis c		unique axis a
//					 				abc  c-ba		abc   ba-c			abc	-acb
//					 			------------		------------ 		------------
//		cell choice 1		  :b1 	:-b1		:c1 	:-c1			:a1	:-a1
//				      2		  :b2 	:-b2		:c2  	:-c2			:a2	:-a2
//				      3		  :b3 	:-b3		:c3  	:-c3			:a3	:-a3
//
//    Orthorhombic	:ba-c		change of basis abc -> ba-c
//							:1			origin choice 1
//							:2ba-c	origin choice 2, change of basis abc -> ba-c
//
//    Tetragonal		:1	origin choice 1
//          Cubic		:2	origin choice 2
//
//    Trigonal		:H	hexagonal    axes
//						:R	rhombohedral axes
//
//
//							SG				  idNum			#IDs		#SpaceGroups
//	Cubic				[195,230]		[489,530]		 42			36
//	Hexagonal		[168,194]		[462,488]		 27			27
//	Trigonal			[143,167]		[430,461]		 32			25
//	Tetragonal		[75,142]			[349,429]		 81			68
//	Orthorhombic	[16,74]			[108,348]		241			59
//	Monoclinic		[3,15]			[3,107]			105			13
//	Triclinic		[1,2]				[1,2]				  2			 2
//
//											SG										idNum
//	Rhombohedral	[146,148,155,160,161,166,167]	[434,437,445,451,453,459,461]


ThreadSafe Function/S getHMboth(SpaceGroupIDnum)	// returns short and (full) Hermann-Mauguin symbol
	Variable SpaceGroupIDnum								//Space Group number, from International Tables

	String short = getHMsym(SpaceGroupIDnum)
	String full = getHMsym2(SpaceGroupIDnum)
	if (StringMatch(short,full))
		return short
	else
		return short + "  ("+full+")"
	endif
End


ThreadSafe Function/T getHMsym(idNum)		// returns short Hermann-Mauguin symbol
	Variable idNum									// index into the SpaceGroup IDs [1,530]
	if (!isValidSpaceGroupIDnum(idNum))
		return ""									// invalid SpaceGroup ID number
	endif

	String HM1=""									// there are 530 items in this list
	HM1  = "P1;P-1;P2:b;P2:c;P2:a;P21:b;P21:c;P21:a;C2:b1;C2:b2;C2:b3;C2:c1;C2:c2;C2:c3;C2:a1;C2:a2;C2:a3;Pm:b;Pm:c;Pm:a;Pc:b1;Pc:b2;Pc:b3;"
	HM1 += "Pc:c1;Pc:c2;Pc:c3;Pc:a1;Pc:a2;Pc:a3;Cm:b1;Cm:b2;Cm:b3;Cm:c1;Cm:c2;Cm:c3;Cm:a1;Cm:a2;Cm:a3;Cc:b1;Cc:b2;Cc:b3;Cc:-b1;Cc:-b2;Cc:-b3;"
	HM1 += "Cc:c1;Cc:c2;Cc:c3;Cc:-c1;Cc:-c2;Cc:-c3;Cc:a1;Cc:a2;Cc:a3;Cc:-a1;Cc:-a2;Cc:-a3;P2/m:b;P2/m:c;P2/m:a;P21/m:b;P21/m:c;P21/m:a;C2/m:b1;"
	HM1 += "C2/m:b2;C2/m:b3;C2/m:c1;C2/m:c2;C2/m:c3;C2/m:a1;C2/m:a2;C2/m:a3;P2/c:b1;P2/c:b2;P2/c:b3;P2/c:c1;P2/c:c2;P2/c:c3;P2/c:a1;P2/c:a2;"
	HM1 += "P2/c:a3;P21/c:b1;P21/c:b2;P21/c:b3;P21/c:c1;P21/c:c2;P21/c:c3;P21/c:a1;P21/c:a2;P21/c:a3;C2/c:b1;C2/c:b2;C2/c:b3;C2/c:-b1;C2/c:-b2;"
	HM1 += "C2/c:-b3;C2/c:c1;C2/c:c2;C2/c:c3;C2/c:-c1;C2/c:-c2;C2/c:-c3;C2/c:a1;C2/c:a2;C2/c:a3;C2/c:-a1;C2/c:-a2;C2/c:-a3;P222;P2221;P2122;"
	HM1 += "P2212;P21212;P22121;P21221;P212121;C2221;A2122;B2212;C222;A222;B222;F222;I222;I212121;Pmm2;P2mm;Pm2m;Pmc21;Pcm21;P21ma;P21am;Pb21m;"
	HM1 += "Pm21b;Pcc2;P2aa;Pb2b;Pma2;Pbm2;P2mb;P2cm;Pc2m;Pm2a;Pca21;Pbc21;P21ab;P21ca;Pc21b;Pb21a;Pnc2;Pcn2;P2na;P2an;Pb2n;Pn2b;Pmn21;Pnm21;"
	HM1 += "P21mn;P21nm;Pn21m;Pm21n;Pba2;P2cb;Pc2a;Pna21;Pbn21;P21nb;P21cn;Pc21n;Pn21a;Pnn2;P2nn;Pn2n;Cmm2;A2mm;Bm2m;Cmc21;Ccm21;A21ma;A21am;"
	HM1 += "Bb21m;Bm21b;Ccc2;A2aa;Bb2b;Amm2;Bmm2;B2mm;C2mm;Cm2m;Am2m;Abm2;Bma2;B2cm;C2mb;Cm2a;Ac2m;Ama2;Bbm2;B2mb;C2cm;Cc2m;Am2a;Aba2;Bba2;"
	HM1 += "B2cb;C2cb;Cc2a;Ac2a;Fmm2;F2mm;Fm2m;Fdd2;F2dd;Fd2d;Imm2;I2mm;Im2m;Iba2;I2cb;Ic2a;Ima2;Ibm2;I2mb;I2cm;Ic2m;Im2a;Pmmm;Pnnn:1;Pnnn:2;"
	HM1 += "Pccm;Pmaa;Pbmb;Pban:1;Pban:2;Pncb:1;Pncb:2;Pcna:1;Pcna:2;Pmma;Pmmb;Pbmm;Pcmm;Pmcm;Pmam;Pnna;Pnnb;Pbnn;Pcnn;Pncn;Pnan;Pmna;Pnmb;"
	HM1 += "Pbmn;Pcnm;Pncm;Pman;Pcca;Pccb;Pbaa;Pcaa;Pbcb;Pbab;Pbam;Pmcb;Pcma;Pccn;Pnaa;Pbnb;Pbcm;Pcam;Pmca;Pmab;Pbma;Pcmb;Pnnm;Pmnn;Pnmn;"
	HM1 += "Pmmn:1;Pmmn:2;Pnmm:1;Pnmm:2;Pmnm:1;Pmnm:2;Pbcn;Pcan;Pnca;Pnab;Pbna;Pcnb;Pbca;Pcab;Pnma;Pmnb;Pbnm;Pcmn;Pmcn;Pnam;Cmcm;Ccmm;Amma;"
	HM1 += "Amam;Bbmm;Bmmb;Cmca;Ccmb;Abma;Acam;Bbcm;Bmab;Cmmm;Ammm;Bmmm;Cccm;Amaa;Bbmb;Cmma;Cmmb;Abmm;Acmm;Bmcm;Bmam;Ccca:1;Ccca:2;Cccb:1;"
	HM1 += "Cccb:2;Abaa:1;Abaa:2;Acaa:1;Acaa:2;Bbcb:1;Bbcb:2;Bbab:1;Bbab:2;Fmmm;Fddd:1;Fddd:2;Immm;Ibam;Imcb;Icma;Ibca;Icab;Imma;Immb;Ibmm;"
	HM1 += "Icmm;Imcm;Imam;P4;P41;P42;P43;I4;I41;P-4;I-4;P4/m;P42/m;P4/n:1;P4/n:2;P42/n:1;P42/n:2;I4/m;I41/a:1;I41/a:2;P422;P4212;P4122;P41212;"
	HM1 += "P4222;P42212;P4322;P43212;I422;I4122;P4mm;P4bm;P42cm;P42nm;P4cc;P4nc;P42mc;P42bc;I4mm;I4cm;I41md;I41cd;P-42m;P-42c;P-421m;P-421c;"
	HM1 += "P-4m2;P-4c2;P-4b2;P-4n2;I-4m2;I-4c2;I-42m;I-42d;P4/mmm;P4/mcc;P4/nbm:1;P4/nbm:2;P4/nnc:1;P4/nnc:2;P4/mbm;P4/mnc;P4/nmm:1;P4/nmm:2;"
	HM1 += "P4/ncc:1;P4/ncc:2;P42/mmc;P42/mcm;P42/nbc:1;P42/nbc:2;P42/nnm:1;P42/nnm:2;P42/mbc;P42/mnm;P42/nmc:1;P42/nmc:2;P42/ncm:1;P42/ncm:2;"
	HM1 += "I4/mmm;I4/mcm;I41/amd:1;I41/amd:2;I41/acd:1;I41/acd:2;P3;P31;P32;R3:H;R3:R;P-3;R-3:H;R-3:R;P312;P321;P3112;P3121;P3212;P3221;R32:H;"
	HM1 += "R32:R;P3m1;P31m;P3c1;P31c;R3m:H;R3m:R;R3c:H;R3c:R;P-31m;P-31c;P-3m1;P-3c1;R-3m:H;R-3m:R;R-3c:H;R-3c:R;P6;P61;P65;P62;P64;P63;P-6;"
	HM1 += "P6/m;P63/m;P622;P6122;P6522;P6222;P6422;P6322;P6mm;P6cc;P63cm;P63mc;P-6m2;P-6c2;P-62m;P-62c;P6/mmm;P6/mcc;P63/mcm;P63/mmc;P23;F23;"
	HM1 += "I23;P213;I213;Pm-3;Pn-3:1;Pn-3:2;Fm-3;Fd-3:1;Fd-3:2;Im-3;Pa-3;Ia-3;P432;P4232;F432;F4132;I432;P4332;P4132;I4132;P-43m;F-43m;I-43m;"
	HM1 += "P-43n;F-43c;I-43d;Pm-3m;Pn-3n:1;Pn-3n:2;Pm-3n;Pn-3m:1;Pn-3m:2;Fm-3m;Fm-3c;Fd-3m:1;Fd-3m:2;Fd-3c:1;Fd-3c:2;Im-3m;Ia-3d;"
	return StringFromList(idNum-1,HM1)
End

ThreadSafe Function/T getHMsym2(idNum)	// returns longer Hermann-Mauguin symbol
	Variable idNum									// index into the SpaceGroup IDs [1,530]
	if (!isValidSpaceGroupIDnum(idNum))
		return ""									// invalid SpaceGroup ID number
	endif

	String HM2=""									// there are 530 items in this list
	HM2  = "P1;P-1;P121;P112;P211;P1211;P1121;P2111;C121;A121;I121;A112;B112;I112;B211;C211;I211;P1m1;P11m;Pm11;P1c1;P1n1;P1a1;P11a;P11n;P11b;"
	HM2 += "Pb11;Pn11;Pc11;C1m1;A1m1;I1m1;A11m;B11m;I11m;Bm11;Cm11;Im11;C1c1;A1n1;I1a1;A1a1;C1n1;I1c1;A11a;B11n;I11b;B11b;A11n;I11a;Bb11;Cn11;"
	HM2 += "Ic11;Cc11;Bn11;Ib11;P12/m1;P112/m;P2/m11;P121/m1;P1121/m;P21/m11;C12/m1;A12/m1;I12/m1;A112/m;B112/m;I112/m;B2/m11;C2/m11;I2/m11;"
	HM2 += "P12/c1;P12/n1;P12/a1;P112/a;P112/n;P112/b;P2/b11;P2/n11;P2/c11;P121/c1;P121/n1;P121/a1;P1121/a;P1121/n;P1121/b;P21/b11;P21/n11;"
	HM2 += "P21/c11;C12/c1;A12/n1;I12/a1;A12/a1;C12/n1;I12/c1;A112/a;B112/n;I112/b;B112/b;A112/n;I112/a;B2/b11;C2/n11;I2/c11;C2/c11;B2/n11;"
	HM2 += "I2/b11;P222;P2221;P2122;P2212;P21212;P22121;P21221;P212121;C2221;A2122;B2212;C222;A222;B222;F222;I222;I212121;Pmm2;P2mm;Pm2m;Pmc21;"
	HM2 += "Pcm21;P21ma;P21am;Pb21m;Pm21b;Pcc2;P2aa;Pb2b;Pma2;Pbm2;P2mb;P2cm;Pc2m;Pm2a;Pca21;Pbc21;P21ab;P21ca;Pc21b;Pb21a;Pnc2;Pcn2;P2na;P2an;"
	HM2 += "Pb2n;Pn2b;Pmn21;Pnm21;P21mn;P21nm;Pn21m;Pm21n;Pba2;P2cb;Pc2a;Pna21;Pbn21;P21nb;P21cn;Pc21n;Pn21a;Pnn2;P2nn;Pn2n;Cmm2;A2mm;Bm2m;"
	HM2 += "Cmc21;Ccm21;A21ma;A21am;Bb21m;Bm21b;Ccc2;A2aa;Bb2b;Amm2;Bmm2;B2mm;C2mm;Cm2m;Am2m;Abm2;Bma2;B2cm;C2mb;Cm2a;Ac2m;Ama2;Bbm2;B2mb;C2cm;"
	HM2 += "Cc2m;Am2a;Aba2;Bba2;B2cb;C2cb;Cc2a;Ac2a;Fmm2;F2mm;Fm2m;Fdd2;F2dd;Fd2d;Imm2;I2mm;Im2m;Iba2;I2cb;Ic2a;Ima2;Ibm2;I2mb;I2cm;Ic2m;Im2a;"
	HM2 += "Pmmm;Pnnn:1;Pnnn:2;Pccm;Pmaa;Pbmb;Pban:1;Pban:2;Pncb:1;Pncb:2;Pcna:1;Pcna:2;Pmma;Pmmb;Pbmm;Pcmm;Pmcm;Pmam;Pnna;Pnnb;Pbnn;Pcnn;Pncn;"
	HM2 += "Pnan;Pmna;Pnmb;Pbmn;Pcnm;Pncm;Pman;Pcca;Pccb;Pbaa;Pcaa;Pbcb;Pbab;Pbam;Pmcb;Pcma;Pccn;Pnaa;Pbnb;Pbcm;Pcam;Pmca;Pmab;Pbma;Pcmb;Pnnm;"
	HM2 += "Pmnn;Pnmn;Pmmn:1;Pmmn:2;Pnmm:1;Pnmm:2;Pmnm:1;Pmnm:2;Pbcn;Pcan;Pnca;Pnab;Pbna;Pcnb;Pbca;Pcab;Pnma;Pmnb;Pbnm;Pcmn;Pmcn;Pnam;Cmcm;"
	HM2 += "Ccmm;Amma;Amam;Bbmm;Bmmb;Cmca;Ccmb;Abma;Acam;Bbcm;Bmab;Cmmm;Ammm;Bmmm;Cccm;Amaa;Bbmb;Cmma;Cmmb;Abmm;Acmm;Bmcm;Bmam;Ccca:1;Ccca:2;"
	HM2 += "Cccb:1;Cccb:2;Abaa:1;Abaa:2;Acaa:1;Acaa:2;Bbcb:1;Bbcb:2;Bbab:1;Bbab:2;Fmmm;Fddd:1;Fddd:2;Immm;Ibam;Imcb;Icma;Ibca;Icab;Imma;Immb;"
	HM2 += "Ibmm;Icmm;Imcm;Imam;P4;P41;P42;P43;I4;I41;P-4;I-4;P4/m;P42/m;P4/n:1;P4/n:2;P42/n:1;P42/n:2;I4/m;I41/a:1;I41/a:2;P422;P4212;P4122;"
	HM2 += "P41212;P4222;P42212;P4322;P43212;I422;I4122;P4mm;P4bm;P42cm;P42nm;P4cc;P4nc;P42mc;P42bc;I4mm;I4cm;I41md;I41cd;P-42m;P-42c;P-421m;"
	HM2 += "P-421c;P-4m2;P-4c2;P-4b2;P-4n2;I-4m2;I-4c2;I-42m;I-42d;P4/mmm;P4/mcc;P4/nbm:1;P4/nbm:2;P4/nnc:1;P4/nnc:2;P4/mbm;P4/mnc;P4/nmm:1;"
	HM2 += "P4/nmm:2;P4/ncc:1;P4/ncc:2;P42/mmc;P42/mcm;P42/nbc:1;P42/nbc:2;P42/nnm:1;P42/nnm:2;P42/mbc;P42/mnm;P42/nmc:1;P42/nmc:2;P42/ncm:1;"
	HM2 += "P42/ncm:2;I4/mmm;I4/mcm;I41/amd:1;I41/amd:2;I41/acd:1;I41/acd:2;P3;P31;P32;R3:H;R3:R;P-3;R-3:H;R-3:R;P312;P321;P3112;P3121;P3212;"
	HM2 += "P3221;R32:H;R32:R;P3m1;P31m;P3c1;P31c;R3m:H;R3m:R;R3c:H;R3c:R;P-31m;P-31c;P-3m1;P-3c1;R-3m:H;R-3m:R;R-3c:H;R-3c:R;P6;P61;P65;P62;"
	HM2 += "P64;P63;P-6;P6/m;P63/m;P622;P6122;P6522;P6222;P6422;P6322;P6mm;P6cc;P63cm;P63mc;P-6m2;P-6c2;P-62m;P-62c;P6/mmm;P6/mcc;P63/mcm;"
	HM2 += "P63/mmc;P23;F23;I23;P213;I213;Pm-3;Pn-3:1;Pn-3:2;Fm-3;Fd-3:1;Fd-3:2;Im-3;Pa-3;Ia-3;P432;P4232;F432;F4132;I432;P4332;P4132;I4132;"
	HM2 += "P-43m;F-43m;I-43m;P-43n;F-43c;I-43d;Pm-3m;Pn-3n:1;Pn-3n:2;Pm-3n;Pn-3m:1;Pn-3m:2;Fm-3m;Fm-3c;Fd-3m:1;Fd-3m:2;Fd-3c:1;Fd-3c:2;Im-3m;"
	HM2 += "Ia-3d;"
	return StringFromList(idNum-1,HM2)
End


//// Full H-M symbols only differ from the regular ones (in getHMsym) for Space Groups: 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
//ThreadSafe Function/S getFullHMSym(SpaceGroup)	// returns full Hermann-Mauguin symbol
//	Variable SpaceGroup					//Space Group number, from International Tables
//	if (!isValidSpaceGroup(SpaceGroup))
//		return ""							// invalid SpaceGroup number
//	endif
//
//	// there are 230 items in this list
//	String fullSymms="P1;P-1;P121:b;P1211:b;C121:b1;P1m1:b;"
//	fullSymms += "P1c1:b1;C1m1:b1;C1c1:b1;P12/m1:b;P121/m1:b;"
//	fullSymms += "C12/m1:b1;P12/c1:b1;P121/c1:b1;C12/c1:b1;P222;"
//	fullSymms += "P2221;P21212;P212121;C2221;C222;F222;I222;I212121;"
//	fullSymms += "Pmm2;Pmc21;Pcc2;Pma2;Pca21;Pnc2;Pmn21;Pba2;Pna21;"
//	fullSymms += "Pnn2;Cmm2;Cmc21;Ccc2;Amm2;Abm2;Ama2;Aba2;Fmm2;"
//	fullSymms += "Fdd2;Imm2;Iba2;Ima2;Pmmm;Pnnn:1;Pccm;Pban:1;Pmma;"
//	fullSymms += "Pnna;Pmna;Pcca;Pbam;Pccn;Pbcm;Pnnm;Pmmn:1;Pbcn;"
//	fullSymms += "Pbca;Pnma;Cmcm;Cmca;Cmmm;Cccm;Cmma;Ccca:1;Fmmm;"
//	fullSymms += "Fddd:1;Immm;Ibam;Ibca;Imma;P4;P41;P42;P43;I4;I41"
//	fullSymms += ";P-4;I-4;P4/m;P42/m;P4/n:1;P42/n:1;I4/m;I41/a:1;"
//	fullSymms += "P422;P4212;P4122;P41212;P4222;P42212;P4322;P43212;"
//	fullSymms += "I422;I4122;P4mm;P4bm;P42cm;P42nm;P4cc;P4nc;P42mc;"
//	fullSymms += "P42bc;I4mm;I4cm;I41md;I41cd;P-42m;P-42c;P-421m;"
//	fullSymms += "P-421c;P-4m2;P-4c2;P-4b2;P-4n2;I-4m2;I-4c2;I-42m;"
//	fullSymms += "I-42d;P4/mmm;P4/mcc;P4/nbm:1;P4/nnc:1;P4/mbm;P4/mnc;"
//	fullSymms += "P4/nmm:1;P4/ncc:1;P42/mmc;P42/mcm;P42/nbc:1;P42/nnm:1;"
//	fullSymms += "P42/mbc;P42/mnm;P42/nmc:1;P42/ncm:1;I4/mmm;I4/mcm;"
//	fullSymms += "I41/amd:1;I41/acd:1;P3;P31;P32;R3:H;P-3;R-3:H;P312;"
//	fullSymms += "P321;P3112;P3121;P3212;P3221;R32:H;P3m1;P31m;P3c1;"
//	fullSymms += "P31c;R3m:H;R3c:H;P-31m;P-31c;P-3m1;P-3c1;R-3m:H;"
//	fullSymms += "R-3c:H;P6;P61;P65;P62;P64;P63;P-6;P6/m;P63/m;P622;"
//	fullSymms += "P6122;P6522;P6222;P6422;P6322;P6mm;P6cc;P63cm;P63mc;"
//	fullSymms += "P-6m2;P-6c2;P-62m;P-62c;P6/mmm;P6/mcc;P63/mcm;"
//	fullSymms += "P63/mmc;P23;F23;I23;P213;I213;Pm-3;Pn-3:1;Fm-3;"
//	fullSymms += "Fd-3:1;Im-3;Pa-3;Ia-3;P432;P4232;F432;F4132;I432;"
//	fullSymms += "P4332;P4132;I4132;P-43m;F-43m;I-43m;P-43n;F-43c;"
//	fullSymms += "I-43d;Pm-3m;Pn-3n:1;Pm-3n;Pn-3m:1;Fm-3m;Fm-3c;"
//	fullSymms += "Fd-3m:1;Fd-3c:1;Im-3m;Ia-3d"
//	return StringFromList(SpaceGroup-1,fullSymms)	// set the symmetry symbol
//End


ThreadSafe Function/S getHallSymbol(idNum)
	Variable idNum									// index into the SpaceGroup IDs [1,530]
	if (!isValidSpaceGroupIDnum(idNum))
		return ""									// invalid SpaceGroup ID number
	endif

	String Hall=""									// there are 530 items in this list
	Hall  = "P 1;-P 1;P 2y;P 2;P 2x;P 2yb;P 2c;P 2xa;C 2y;A 2y;I 2y;A 2;B 2;I 2;B 2x;C 2x;I 2x;P -2y;P -2;P -2x;P -2yc;P -2yac;P -2ya;P -2a;"
	Hall += "P -2ab;P -2b;P -2xb;P -2xbc;P -2xc;C -2y;A -2y;I -2y;A -2;B -2;I -2;B -2x;C -2x;I -2x;C -2yc;A -2yac;I -2ya;A -2ya;C -2ybc;I -2yc;"
	Hall += "A -2a;B -2bc;I -2b;B -2b;A -2ac;I -2a;B -2xb;C -2xbc;I -2xc;C -2xc;B -2xbc;I -2xb;-P 2y;-P 2;-P 2x;-P 2yb;-P 2c;-P 2xa;-C 2y;"
	Hall += "-A 2y;-I 2y;-A 2;-B 2;-I 2;-B 2x;-C 2x;-I 2x;-P 2yc;-P 2yac;-P 2ya;-P 2a;-P 2ab;-P 2b;-P 2xb;-P 2xbc;-P 2xc;-P 2ybc;-P 2yn;"
	Hall += "-P 2yab;-P 2ac;-P 2n;-P 2bc;-P 2xab;-P 2xn;-P 2xac;-C 2yc;-A 2yac;-I 2ya;-A 2ya;-C 2ybc;-I 2yc;-A 2a;-B 2bc;-I 2b;-B 2b;-A 2ac;"
	Hall += "-I 2a;-B 2xb;-C 2xbc;-I 2xc;-C 2xc;-B 2xbc;-I 2xb;P 2 2;P 2c 2;P 2a 2a;P 2 2b;P 2 2ab;P 2bc 2;P 2ac 2ac;P 2ac 2ab;C 2c 2;A 2a 2a;"
	Hall += "B 2 2b;C 2 2;A 2 2;B 2 2;F 2 2;I 2 2;I 2b 2c;P 2 -2;P -2 2;P -2 -2;P 2c -2;P 2c -2c;P -2a 2a;P -2 2a;P -2 -2b;P -2b -2;P 2 -2c;"
	Hall += "P -2a 2;P -2b -2b;P 2 -2a;P 2 -2b;P -2b 2;P -2c 2;P -2c -2c;P -2a -2a;P 2c -2ac;P 2c -2b;P -2b 2a;P -2ac 2a;P -2bc -2c;P -2a -2ab;"
	Hall += "P 2 -2bc;P 2 -2ac;P -2ac 2;P -2ab 2;P -2ab -2ab;P -2bc -2bc;P 2ac -2;P 2bc -2bc;P -2ab 2ab;P -2 2ac;P -2 -2bc;P -2ab -2;P 2 -2ab;"
	Hall += "P -2bc 2;P -2ac -2ac;P 2c -2n;P 2c -2ab;P -2bc 2a;P -2n 2a;P -2n -2ac;P -2ac -2n;P 2 -2n;P -2n 2;P -2n -2n;C 2 -2;A -2 2;B -2 -2;"
	Hall += "C 2c -2;C 2c -2c;A -2a 2a;A -2 2a;B -2 -2b;B -2b -2;C 2 -2c;A -2a 2;B -2b -2b;A 2 -2;B 2 -2;B -2 2;C -2 2;C -2 -2;A -2 -2;A 2 -2c;"
	Hall += "B 2 -2c;B -2c 2;C -2b 2;C -2b -2b;A -2c -2c;A 2 -2a;B 2 -2b;B -2b 2;C -2c 2;C -2c -2c;A -2a -2a;A 2 -2ac;B 2 -2bc;B -2bc 2;"
	Hall += "C -2bc 2;C -2bc -2bc;A -2ac -2ac;F 2 -2;F -2 2;F -2 -2;F 2 -2d;F -2d 2;F -2d -2d;I 2 -2;I -2 2;I -2 -2;I 2 -2c;I -2a 2;I -2b -2b;"
	Hall += "I 2 -2a;I 2 -2b;I -2b 2;I -2c 2;I -2c -2c;I -2a -2a;-P 2 2;P 2 2 -1n;-P 2ab 2bc;-P 2 2c;-P 2a 2;-P 2b 2b;P 2 2 -1ab;-P 2ab 2b;"
	Hall += "P 2 2 -1bc;-P 2b 2bc;P 2 2 -1ac;-P 2a 2c;-P 2a 2a;-P 2b 2;-P 2 2b;-P 2c 2c;-P 2c 2;-P 2 2a;-P 2a 2bc;-P 2b 2n;-P 2n 2b;-P 2ab 2c;"
	Hall += "-P 2ab 2n;-P 2n 2bc;-P 2ac 2;-P 2bc 2bc;-P 2ab 2ab;-P 2 2ac;-P 2 2bc;-P 2ab 2;-P 2a 2ac;-P 2b 2c;-P 2a 2b;-P 2ac 2c;-P 2bc 2b;"
	Hall += "-P 2b 2ab;-P 2 2ab;-P 2bc 2;-P 2ac 2ac;-P 2ab 2ac;-P 2ac 2bc;-P 2bc 2ab;-P 2c 2b;-P 2c 2ac;-P 2ac 2a;-P 2b 2a;-P 2a 2ab;-P 2bc 2c;"
	Hall += "-P 2 2n;-P 2n 2;-P 2n 2n;P 2 2ab -1ab;-P 2ab 2a;P 2bc 2 -1bc;-P 2c 2bc;P 2ac 2ac -1ac;-P 2c 2a;-P 2n 2ab;-P 2n 2c;-P 2a 2n;"
	Hall += "-P 2bc 2n;-P 2ac 2b;-P 2b 2ac;-P 2ac 2ab;-P 2bc 2ac;-P 2ac 2n;-P 2bc 2a;-P 2c 2ab;-P 2n 2ac;-P 2n 2a;-P 2c 2n;-C 2c 2;-C 2c 2c;"
	Hall += "-A 2a 2a;-A 2 2a;-B 2 2b;-B 2b 2;-C 2bc 2;-C 2bc 2bc;-A 2ac 2ac;-A 2 2ac;-B 2 2bc;-B 2bc 2;-C 2 2;-A 2 2;-B 2 2;-C 2 2c;-A 2a 2;"
	Hall += "-B 2b 2b;-C 2b 2;-C 2b 2b;-A 2c 2c;-A 2 2c;-B 2 2c;-B 2c 2;C 2 2 -1bc;-C 2b 2bc;C 2 2 -1bc;-C 2b 2c;A 2 2 -1ac;-A 2a 2c;"
	Hall += "A 2 2 -1ac;-A 2ac 2c;B 2 2 -1bc;-B 2bc 2b;B 2 2 -1bc;-B 2b 2bc;-F 2 2;F 2 2 -1d;-F 2uv 2vw;-I 2 2;-I 2 2c;-I 2a 2;-I 2b 2b;"
	Hall += "-I 2b 2c;-I 2a 2b;-I 2b 2;-I 2a 2a;-I 2c 2c;-I 2 2b;-I 2 2a;-I 2c 2;P 4;P 4w;P 4c;P 4cw;I 4;I 4bw;P -4;I -4;-P 4;-P 4c;P 4ab -1ab;"
	Hall += "-P 4a;P 4n -1n;-P 4bc;-I 4;I 4bw -1bw;-I 4ad;P 4 2;P 4ab 2ab;P 4w 2c;P 4abw 2nw;P 4c 2;P 4n 2n;P 4cw 2c;P 4nw 2abw;I 4 2;"
	Hall += "I 4bw 2bw;P 4 -2;P 4 -2ab;P 4c -2c;P 4n -2n;P 4 -2c;P 4 -2n;P 4c -2;P 4c -2ab;I 4 -2;I 4 -2c;I 4bw -2;I 4bw -2c;P -4 2;P -4 2c;"
	Hall += "P -4 2ab;P -4 2n;P -4 -2;P -4 -2c;P -4 -2ab;P -4 -2n;I -4 -2;I -4 -2c;I -4 2;I -4 2bw;-P 4 2;-P 4 2c;P 4 2 -1ab;-P 4a 2b;"
	Hall += "P 4 2 -1n;-P 4a 2bc;-P 4 2ab;-P 4 2n;P 4ab 2ab -1ab;-P 4a 2a;P 4ab 2n -1ab;-P 4a 2ac;-P 4c 2;-P 4c 2c;P 4n 2c -1n;-P 4ac 2b;"
	Hall += "P 4n 2 -1n;-P 4ac 2bc;-P 4c 2ab;-P 4n 2n;P 4n 2n -1n;-P 4ac 2a;P 4n 2ab -1n;-P 4ac 2ac;-I 4 2;-I 4 2c;I 4bw 2bw -1bw;-I 4bd 2;"
	Hall += "I 4bw 2aw -1bw;-I 4bd 2c;P 3;P 31;P 32;R 3;P 3*;-P 3;-R 3;-P 3*;P 3 2;P 3 2\";P 31 2c (0 0 1);P 31 2\";P 32 2c (0 0 -1);P 32 2\";"
	Hall += "R 3 2\";P 3* 2;P 3 -2\";P 3 -2;P 3 -2\"c;P 3 -2c;R 3 -2\";P 3* -2;R 3 -2\"c;P 3* -2n;-P 3 2;-P 3 2c;-P 3 2\";-P 3 2\"c;-R 3 2\";"
	Hall += "-P 3* 2;-R 3 2\"c;-P 3* 2n;P 6;P 61;P 65;P 62;P 64;P 6c;P -6;-P 6;-P 6c;P 6 2;P 61 2 (0 0 -1);P 65 2 (0 0 1);P 62 2c (0 0 1);"
	Hall += "P 64 2c (0 0 -1);P 6c 2c;P 6 -2;P 6 -2c;P 6c -2;P 6c -2c;P -6 2;P -6c 2;P -6 -2;P -6c -2c;-P 6 2;-P 6 2c;-P 6c 2;-P 6c 2c;P 2 2 3;"
	Hall += "F 2 2 3;I 2 2 3;P 2ac 2ab 3;I 2b 2c 3;-P 2 2 3;P 2 2 3 -1n;-P 2ab 2bc 3;-F 2 2 3;F 2 2 3 -1d;-F 2uv 2vw 3;-I 2 2 3;-P 2ac 2ab 3;"
	Hall += "-I 2b 2c 3;P 4 2 3;P 4n 2 3;F 4 2 3;F 4d 2 3;I 4 2 3;P 4acd 2ab 3;P 4bd 2ab 3;I 4bd 2c 3;P -4 2 3;F -4 2 3;I -4 2 3;P -4n 2 3;"
	Hall += "F -4c 2 3;I -4bd 2c 3;-P 4 2 3;P 4 2 3 -1n;-P 4a 2bc 3;-P 4n 2 3;P 4n 2 3 -1n;-P 4bc 2bc 3;-F 4 2 3;-F 4c 2 3;F 4d 2 3 -1d;"
	Hall += "-F 4vw 2vw 3;F 4d 2 3 -1cd;-F 4cvw 2vw 3;-I 4 2 3;-I 4bd 2c 3;"
	return StringFromList(idNum-1,Hall)
End


//ThreadSafe Function/S getHallSymbol(SpaceGroup)
//	Variable SpaceGroup					//Space Group number, from International Tables
//	if (!isValidSpaceGroup(SpaceGroup))
//		return ""							// invalid SpaceGroup number
//	endif
//
//	// there are 230 items in this list
//	String Hall="P 1;-P 1;P 2y:b;P 2yb:b;C 2y:b1;P -2y:b;P -2yc:b1;"
//	Hall += "C -2y:b1;C -2yc:b1;-P 2y:b;-P 2yb:b;-C 2y:b1;-P 2yc:b1;"
//	Hall += "-P 2ybc:b1;-C 2yc:b1;P 2 2;P 2c 2;P 2 2ab;P 2ac 2ab;"
//	Hall += "C 2c 2;C 2 2;F 2 2;I 2 2;I 2b 2c;P 2 -2;P 2c -2;"
//	Hall += "P 2 -2c;P 2 -2a;P 2c -2ac;P 2 -2bc;P 2ac -2;P 2 -2ab;"
//	Hall += "P 2c -2n;P 2 -2n;C 2 -2;C 2c -2;C 2 -2c;A 2 -2;A 2 -2c;"
//	Hall += "A 2 -2a;A 2 -2ac;F 2 -2;F 2 -2d;I 2 -2;I 2 -2c;I 2 -2a;"
//	Hall += "-P 2 2;P 2 2 -1n:1;-P 2 2c;P 2 2 -1ab:1;-P 2a 2a;"
//	Hall += "-P 2a 2bc;-P 2ac 2;-P 2a 2ac;-P 2 2ab;-P 2ab 2ac;"
//	Hall += "-P 2c 2b;-P 2 2n;P 2 2ab -1ab:1;-P 2n 2ab;-P 2ac 2ab;"
//	Hall += "-P 2ac 2n;-C 2c 2;-C 2bc 2;-C 2 2;-C 2 2c;-C 2b 2;"
//	Hall += "C 2 2 -1bc:1;-F 2 2;F 2 2 -1d:1;-I 2 2;-I 2 2c;-I 2b 2c;"
//	Hall += "-I 2b 2;P 4;P 4w;P 4c;P 4cw;I 4;I 4bw;P -4;I -4;-P 4;"
//	Hall += "-P 4c;P 4ab -1ab:1;P 4n -1n:1;-I 4;I 4bw -1bw:1;P 4 2;"
//	Hall += "P 4ab 2ab;P 4w 2c;P 4abw 2nw;P 4c 2;P 4n 2n;P 4cw 2c;"
//	Hall += "P 4nw 2abw;I 4 2;I 4bw 2bw;P 4 -2;P 4 -2ab;P 4c -2c;"
//	Hall += "P 4n -2n;P 4 -2c;P 4 -2n;P 4c -2;P 4c -2ab;I 4 -2;"
//	Hall += "I 4 -2c;I 4bw -2;I 4bw -2c;P -4 2;P -4 2c;P -4 2ab;"
//	Hall += "P -4 2n;P -4 -2;P -4 -2c;P -4 -2ab;P -4 -2n;I -4 -2;"
//	Hall += "I -4 -2c;I -4 2;I -4 2bw;-P 4 2;-P 4 2c;P 4 2 -1ab:1;"
//	Hall += "P 4 2 -1n:1;-P 4 2ab;-P 4 2n;P 4ab 2ab -1ab:1;"
//	Hall += "P 4ab 2n -1ab:1;-P 4c 2;-P 4c 2c;P 4n 2c -1n:1;"
//	Hall += "P 4n 2 -1n:1;-P 4c 2ab;-P 4n 2n;P 4n 2n -1n:1;"
//	Hall += "P 4n 2ab -1n:1;-I 4 2;-I 4 2c;I 4bw 2bw -1bw:1;"
//	Hall += "I 4bw 2aw -1bw:1;P 3;P 31;P 32;R 3:H;-P 3;-R 3:H;"
//	Hall += "P 3 2;P 3 2'';P 31 2c (0 0 1);P 31 2'';P 32 2c (0 0 -1);"
//	Hall += "P 32 2'';R 3 2'':H;P 3 -2'';P 3 -2;P 3 -2''c;P 3 -2c;R 3 -2'':H;"
//	Hall += "R 3 -2''c:H;-P 3 2;-P 3 2c;-P 3 2'';-P 3 2''c;-R 3 2'':H;"
//	Hall += "-R 3 2''c:H;P 6;P 61;P 65;P 62;P 64;P 6c;P -6;-P 6;-P 6c;"
//	Hall += "P 6 2;P 61 2 (0 0 -1);P 65 2 (0 0 1);P 62 2c (0 0 1);"
//	Hall += "P 64 2c (0 0 -1);P 6c 2c;P 6 -2;P 6 -2c;P 6c -2;P 6c -2c;"
//	Hall += "P -6 2;P -6c 2;P -6 -2;P -6c -2c;-P 6 2;-P 6 2c;-P 6c 2;"
//	Hall += "-P 6c 2c;P 2 2 3;F 2 2 3;I 2 2 3;P 2ac 2ab 3;I 2b 2c 3;"
//	Hall += "-P 2 2 3;P 2 2 3 -1n:1;-F 2 2 3;F 2 2 3 -1d:1;-I 2 2 3;"
//	Hall += "-P 2ac 2ab 3;-I 2b 2c 3;P 4 2 3;P 4n 2 3;F 4 2 3;F 4d 2 3;"
//	Hall += "I 4 2 3;P 4acd 2ab 3;P 4bd 2ab 3;I 4bd 2c 3;P -4 2 3;"
//	Hall += "F -4 2 3;I -4 2 3;P -4n 2 3;F -4c 2 3;I -4bd 2c 3;-P 4 2 3;"
//	Hall += "P 4 2 3 -1n:1;-P 4n 2 3;P 4n 2 3 -1n:1;-F 4 2 3;-F 4c 2 3;"
//	Hall += "F 4d 2 3 -1d:1;F 4d 2 3 -1cd:1;-I 4 2 3;-I 4bd 2c 3;"
//	return StringFromList(SpaceGroup-1,Hall)	// set the Hall symbol
//End


ThreadSafe Static Function latticeSystem(SpaceGroupID)
	String SpaceGroupID

	Variable SG=str2num(SpaceGroupID)	//Space Group number, from International Tables [1,230]
	if (SG>230)
		return -1					  	// invalid
	elseif (SG>=195)
		return CUBIC
	elseif (SG>=168)
		return HEXAGONAL
	elseif (SG>=143)
		return TRIGONAL				// Trigonal, (using the hexagonal cell axes)
	elseif (SG>=75)
		return TETRAGONAL
	elseif (SG>=16)
		return ORTHORHOMBIC
	elseif (SG>=3)
		return MONOCLINIC
	elseif  (SG>0)
		return TRICLINIC
	endif
	return -1							// invalid
End


Function/S symmtry2SG(strIN,[types,printIt])	// find the Space Group number from the symmetry string
	String strIN
	Variable types						// -1=Check All, 1=Hermann-Mauguin, 2=Full Hermann-Mauguin, 4=Hall, 8=Lattice System, 16=SpaceGroupID, 32=Ignore Minuses
	Variable printIt
	types = ParamIsDefault(types) ? -1 : round(types)
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt

	if (strlen(strIN)<1 || numtype(types))
		Prompt strIN, "Symmetry Symbol or Space Group Number, (e.g. \"Pmm*\"), wild cards allowed"
		Variable t1=!(!(types&1))+1, t2=!(!(types&2))+1, t4=!(!(types&4))+1, t8=!(!(types&8))+1, t16=!(!(types&16))+1, t32=!(!(types&32))+1
		Prompt t1, "Hermann-Mauguin", popup "-;Hermann-Mauguin"
		Prompt t2, "Full Hermann-Mauguin", popup "-;Full Hermann-Mauguin"
		Prompt t4, "Hall", popup "-;Hall"
		Prompt t8, "Lattice System", popup "-;Lattice System"
		Prompt t16, "SpaceGroup ID", popup "-;SpaceGroup ID"
		Prompt t32, "Ignore All Minus Signs", popup "-;Ignore Minus Signs"
		DoPrompt "Symmetry Symbol",strIN,t1,t2,t4,t8,t16,t32
		if (V_flag)
			return ""
		endif
		types = 0
		types += t1==2 ? 1 : 0
		types += t2==2 ? 2 : 0
		types += t4==2 ? 4 : 0
		types += t8==2 ? 8 : 0
		types += t16==2 ? 16 : 0
		types += t32==2 ? 32 : 0
		printIt = 1
	endif
	if (strlen(strIN)<1)
		return ""
	endif
	Variable ignoreMinus = !(!(types & 32))					// 1 means ignore minus signs in matching

	String list="", nameList=""
	Variable idNum
	if (types & 1)
		list += SymString2SGtype(strIN,1,ignoreMinus)		// 1 = Hermann-Mauguin
		nameList += "Hermann-Mauguin, "
	endif
	if (types & 2)
		list += SymString2SGtype(strIN,2,ignoreMinus)		// 2 = Full Hermann-Mauguin, HM2
		nameList += "FULL Hermann-Mauguin, "
	endif
	if (types & 4)
		list += SymString2SGtype(strIN,4,ignoreMinus)		// 4 = Hall
		nameList += "Hall, "
	endif
	if (types & 8)
		list += SymString2SGtype(strIN,8,ignoreMinus)		// 8 = Lattice System
		nameList += "Lattice System, "
	endif
	if (types & 16)
		list += SymString2SGtype(strIN,16,ignoreMinus)		// 16 = space group ID, e.g. "15:b3"
		nameList += "Space Group ID, "
	endif
	nameList = TrimBoth(nameList,chars=", ")
	list = RemoveDuplicatesFromList(list)
	list = SortList(list,";",1)

	if (printIt)	// print out information about each SpaceGroupIDnum in list
		Variable i, Nlist=ItemsInList(list), showDefault
		if (Nlist<1)
			printf "No matches of \"%s\" to a symbol in:{%s}\r",strIN,nameList
			return ""
		elseif (Nlist>1)
			printf "There are %g possible matches of  \"%s\"  to a symbol in: {%s}\r",Nlist,strIN,nameList
		endif
		String allIDs=MakeAllIDs()
		printf "\t\tSG id\t\t\t\tSystem\t\t\t\tH-M\t\t\tHall\t\t\tfull H-M\r"
		String id, tab,fullHM,HM, system, systemNames="Triclinic\t;Monoclinic\t;Orthorhombic;Tetragonal\t;Trigonal\t;Hexagonal\t;Cubic\t\t"
		for (i=0; i<Nlist; i+=1)
			idNum = str2num(StringFromList(i,list))
			if (isValidSpaceGroupIDnum(idNum))
				fullHM = getHMsym2(idNum)			// usually fullHM is the same as HM
				HM = getHMSym(idNum)
				fullHM = SelectString(StringMatch(fullHM,HM),"\t\t\t"+fullHM,"")
				tab = SelectString(strlen(getHMsym2(idNum))>5,"\t","")

				id = StringFromList(idNum-1,allIDs)
				system = StringFromList(latticeSystem(id),systemNames)
				showDefault = strsearch(id,":",0)>=0 && StringMatch(id,FindDefaultIDforSG(str2num(id)))
				printf "%8s\t-->\t\t%s\t\t%s\t\t%s%s%s%s\r", id,system,HM,tab,getHallSymbol(idNum),fullHM,SelectString(showDefault,"","\t--DEFAULT--")
			endif
		endfor
	endif
	return list
End


Static Function/S SymString2SGtype(symIN,type,ignoreMinus)
	// finds space group of a Hermann-Mauguin or Hall symbol, wild cards allowed, return list of idNums [1,530]
	String symIN						// requested symbol, if empty, then a dialog will come up
	Variable type						// 1=Hermann-Mauguin, 2=Full Hermann-Mauguin, 4=Hall, 8=Lattice System, 16=space group ID
	Variable ignoreMinus			// if true, then ignore any minus signs when matching

	String find = ReplaceString(" ",symIN,"")	// do not include spaces in search
	if (ignoreMinus)
		find = ReplaceString("-",find,"")			// ignore all minus signs too
	endif

	type = round(type)
	if (!(type & 31))					// return "" if type not 1, 2, 4, 8, or 16
		return ""
	endif

	String sym,list=""
	Variable idNum
	if (type==8)						// searching for a Lattice System
		if (StringMatch("Triclinic",find))
			list += expandRange("1-2",";")+";"
		endif
		if (StringMatch("Monoclinic",find))
			list += expandRange("3-107",";")+";"
		endif
		if (StringMatch("Orthorhombic",find))
			list += expandRange("108-348",";")+";"
		endif
		if (StringMatch("Tetragonal",find))
			list += expandRange("349-429",";")+";"
		endif
		if (StringMatch("Trigonal",find))
			list += expandRange("430-461",";")+";"
		elseif (StringMatch("Rhombohedral",find))
			list += expandRange("434,437,445,451,453,459,461",";")+";"
		endif
		if (StringMatch("Hexagonal",find))
			list += expandRange("462-488",";")+";"
		endif
		if (StringMatch("Cubic",find))
			list += expandRange("489-530",";")+";"
		endif
		return list
	endif

	if (type==16)											// searching for a Space Group ID, e.g. "15:b3"
		String allIDs=MakeAllIDs()
		for (idNum=1; idNum<=530; idNum+=1)
			sym = StringFromList(idNum-1,allIDs)
			if (ignoreMinus)
				sym = ReplaceString("-",sym,"")		// do not include minus signs in match
			endif
			if (StringMatch(sym, find))
				list += num2istr(idNum)+";"			// found a match, save it
			endif
		endfor
		return list
	endif

	if (type==1)						// searching for a Hermann-Mauguin
		FUNCREF getHMsym symbolFunc = getHMsym
	elseif (type==2)					// searching for a Full Hermann-Mauguin
		FUNCREF getHMsym symbolFunc = getHMsym2
	elseif (type==4)					// searching for a Hall symbol
		FUNCREF getHMsym symbolFunc = getHallSymbol
	endif

	// symbolFunc has been set, check all idNums in symbolFunc(idNum)
	for (idNum=1;idNum<=530;idNum+=1)			// check all 530 space group types using symbolFunc
		sym = ReplaceString(" ",symbolFunc(idNum),"")	// ignore spaces
		if (ignoreMinus)
			sym = ReplaceString("-",sym,"")		// optionally, do not include minus signs
		endif

		if (StringMatch(sym,find))				// look for find in sym
			list += num2istr(idNum)+";"			// found a match, save it
		endif
	endfor
	return list
End


//	DEPRECATED	DEPRECATED	DEPRECATED	DEPRECATED
// This function is DEPRECATED, it is just left here for old stuff, use getHMsym() instead
//	ThreadSafe Function/S getSymString(SpaceGroup)	// returns short Hermann-Mauguin symbol
//		Variable SpaceGroup							//Space Group number, from International Tables
//		return getHMsym(SpaceGroup)
//	End



//		=========================================================================
//		===================== Start Hex to Rhom Conversions =====================

// Convert direct lattices back and forth between Hexagonal & Rhombohedral
// see:  https://quantumwise.com/support/tutorials/item/510-rhombohedral-and-hexagonal-settings-of-trigonal-crystals
ThreadSafe Function/WAVE RhomLatticeFromHex(directH)			// returns a Rhombohedral direct lattice
	Wave directH												// Hexagonal direct lattice
	// should only be used for Space Groups, {146, 148, 155,160, 161, 166, 167}
	Make/N=(3,3)/D/FREE HexRhomMat = {{2,1,1}, {-1,1,1}, {-1,-2,1}}
	HexRhomMat /= 3
	MatrixOP/FREE directR = directH x HexRhomMat	// Obverse Rhombohedral direct lattice from Hexagonal
	return directR
End


ThreadSafe Function/WAVE HexLatticeFromRhom(directR)			// returns a Hexagonal direct lattice
	Wave directR												// Rhombohedral direct lattice
	// should only be used for Space Groups, {146, 148, 155,160, 161, 166, 167}
	Make/N=(3,3)/D/FREE RhomHexMat = { {1,-1,0}, {0,1,-1}, {1,1,1} }
	MatrixOP/FREE directH = directR x RhomHexMat	// Hexgonal direct lattice from Obverse Rhombohedral
	return directH
End


ThreadSafe Function/C Hex2Rhom(aH,cH)			// convert lattice constants
	Variable aH,cH					// Hexagonal lattice constants
	Variable aR,alpha				// Rhombohedral lattice constants
	aR = sqrt(3*aH*aH + cH*cH) / 3
	alpha = 2 * asin( 1.5 / sqrt(3+(cH/aH)^2) )
	return cmplx(aR, alpha * 180/PI)
End


ThreadSafe Function/C Rhom2Hex(aR,alpha)		// convert lattice constants
	Variable aR, alpha			// Rhombohedral lattice constants
	Variable aH,cH					// Hexagonal lattice constants
	alpha *= PI / 180.0			// convert degree --> radian
	aH = 2 * aR * sin(alpha/2)
	cH = aR * sqrt( 3.0 + 6.0 * cos(alpha) )
	return cmplx(aH,cH)
End


// This is just a demonstration, change it to make it useful
Static Function Rhom2HexFractonal(xtal)					// converts hexagonal --> rhombohedral, fractional coordinates
	STRUCT crystalStructure &xtal

	Variable/C a_c = Rhom2Hex(xtal.a,xtal.alpha)
	print "   Starting from Rhombohedral"
	printf "aHex = %g nm,   cHex = %g nm\r",real(a_c), imag(a_c)

	Wave directR = directFrom_xtal(xtal)			// Obverse Rhombohedral direct lattice
	Wave directH = HexLatticeFromRhom(directR)	// returns a Hexagonal direct lattice

	Make/N=3/D/FREE xyzR								// one atom's fractional rhombohedral coordinates
	Variable i
	for (i=0; i < xtal.N; i+=1)
		xyzR = {xtal.atom[i].x, xtal.atom[i].y, xtal.atom[i].z}
		MatrixOP/FREE xyzH = Inv(directH) x directR x xyzR
		MatrixOp/FREE xyzH = xyzH - floor(xyzH)	// atom's fractional hexagonal coordinates
		xyzH = mod(xyzH,1)
		xyzH = xyzH<1e-12 ? 0 : mod(xyzH,1)		// a fractional coord < 1e-12 is 0
		printf "fractional: Rhom=%s  -->  Hex=%s\r",vec2str(xyzR,zeroThresh=1e-12),vec2str(xyzH,zeroThresh=1e-12)
	endfor
End
//
// This is just a demonstration, change it to make it useful
Static Function Hex2RhomFractonal(xtal)			// converts hexagonal --> rhombohedral, fractional coordinates
	STRUCT crystalStructure &xtal

	Variable/C a_alpha = Hex2Rhom(xtal.a, xtal.c)
	print "   Starting from Hexagonal"
	printf "aRhom = %g nm,   alphaRhom = %g¡\r",real(a_alpha), imag(a_alpha)

	Wave directH = directFrom_xtal(xtal)			// Hexagonal direct lattice
	Wave directR = RhomLatticeFromHex(directH)	// Obverse Rhombohedral direct lattice from Hexagonal

	Make/N=3/D/FREE xyzH								// one atom's fractional hexagonal coordinates
	Variable i
	for (i=0; i < xtal.N; i+=1)
		xyzH = {xtal.atom[i].x, xtal.atom[i].y, xtal.atom[i].z}
		MatrixOP/FREE xyzR = Inv(directR) x directH x xyzH
		MatrixOp/FREE xyzR = xyzR - floor(xyzR)	// atom's fractional rhombohedral coordinates
		xyzR = xyzR<1e-12 ? 0 : mod(xyzR,1)		// a fractional coord < 1e-12 is 0
		printf "fractional: Hex=%s  -->  Rhom=%s\r",vec2str(xyzH,zeroThresh=1e-12),vec2str(xyzR,zeroThresh=1e-12)
	endfor
End
//
//	Static Function test_Hex_Rhom()
//		STRUCT crystalStructure xtal
//	
//		//	FillCrystalStructDefault(xtal)
//	 	// ************************ set xtal to Hexagonal values *************************
//		xtal.desc = "Al2O3 Sapphire (hexagonal)"
//		xtal.a = 0.4758	;	xtal.b = 0.4758	;	xtal.c = 1.2991
//		xtal.alpha = 90	;	xtal.beta = 90		;	xtal.gam = 120
//		xtal.SpaceGroup = 167
//		xtal.a0 = 0.4758	;	xtal.b0 = -0.2379					;	xtal.c0 = 0
//		xtal.a1 = 0			;	xtal.b1 = 0.41205488712064	;	xtal.c1 = 0
//		xtal.a2 = 0			;	xtal.b2 = 0							;	xtal.c2 = 1.2991
//		xtal.as0 = 13.205517669566	;	xtal.bs0 = 0		;	xtal.cs0 = 0
//		xtal.as1 = 7.6242091813124	;	xtal.bs1 = 15.248418362625	;	xtal.cs1 = -0
//		xtal.as2 = 0						;	xtal.bs2 = 0		;	xtal.cs2 = 4.8365678601952
//		xtal.Vc = 0.25469597973584	;	xtal.density = 3.9885352190818
//		xtal.Temperature = 22.5		;	xtal.alphaT = nan
//	
//		xtal.N = 2
//		xtal.atom[0].name = "Al"
//		xtal.atom[0].Zatom = 13
//		xtal.atom[0].x = 0				;	xtal.atom[0].y = 0	;	xtal.atom[0].z = 0.3523
//		xtal.atom[0].mult = 12			;	xtal.atom[0].occ = 1;	xtal.atom[0].WyckoffSymbol = "c"
//		xtal.atom[0].valence = 0		;	xtal.atom[0].DebyeT = 1047
//	
//		xtal.atom[1].name = "O"
//		xtal.atom[1].Zatom = 8
//		xtal.atom[1].x = 0.3064		;	xtal.atom[1].y = 0	;	xtal.atom[1].z = 0.25
//		xtal.atom[1].mult = 18			;	xtal.atom[1].occ = 1;	xtal.atom[1].WyckoffSymbol = "e"
//		xtal.atom[1].valence = 0		;	xtal.atom[1].DebyeT = 1047
//	
//		Hex2RhomFractonal(xtal)				// *********************************
//	
//	
//	 	// *********************** set xtal to Rhombohedral values ***********************
//		xtal.a = 0.512815510469192		;	xtal.b = 0.512815510469192	;	xtal.c = 0.512815510469192
//		xtal.alpha = 55.2793424027957	;	xtal.beta = 55.2793424027957	;	xtal.gam = 55.2793424027957
//	
//		Make/N=(3,3)/D/FREE DR
//		DR[0][0]= {0.2379,0.137351629040212,0.433033333333333}
//		DR[0][1]= {-0.2379,0.137351629040212,0.433033333333333}
//		DR[0][2]= {-5.6751082567348e-17,-0.274703258080424,0.433033333333333}
//		xtal.a0 = DR[0][0]	;	xtal.b0 = DR[0][1]	;		xtal.c0 = DR[0][2]
//		xtal.a1 = DR[1][0]	;	xtal.b1 = DR[1][1]	;		xtal.c1 = DR[1][2]
//		xtal.a2 = DR[2][0]	;	xtal.b2 = DR[2][1]	;		xtal.c2 = DR[2][2]
//	
//		xtal.atom[0].x = 0.3523	;	xtal.atom[0].y = 0.3523	;	xtal.atom[0].z = 0.3523
//		xtal.atom[1].x = 0.5564	;	xtal.atom[1].y = 0.9436	;	xtal.atom[1].z = 0.25
//		xtal.atom[1].y -= 1
//		print " "
//	
//		Rhom2HexFractonal(xtal)				// *********************************
//	End

//		====================== End Hex to Rhom Conversions ======================
//		=========================================================================



//	#define ALLOW_FC(H,K,L) (!(((H)+(K))%2) && !(((K)+(L))%2))		// H,K,L all even or all odd
//	#define ALLOW_BC(H,K,L) (!(((H)+(K)+(L))%2))					// !mod(round(h+k+l),2)
//	#define ALLOW_CC(H,K,L) (!((H)+(K))%2)							// !mod(round(h+k),2)
//	#define ALLOW_AC(H,K,L) (!((K)+(L))%2)							// !mod(round(k+l),2)
//	#define ALLOW_RHOM_HEX(H,K,L) (((-(H)+(K)+(L))%3)==0 || (((H)-(K)+(L))%3)==0)   // allowed are -H+K+L=3n or H-K+L=3n
//	#define ALLOW_HEXAGONAL(H,K,L) ((((H)+2*(K))%3) || !((L)%2))   // forbidden are: H+2K=3N with L odd
//
ThreadSafe Static Function ALLOW_FC(h,k,l)			// face-centered, hkl must be all even or all odd
	Variable h,k,l
	return !mod(h+k,2) && !mod(k+l,2)
End
//
ThreadSafe Static Function ALLOW_BC(h,k,l)			// body-centered, !mod(round(h+k+l),2), sum must be even
	Variable h,k,l
	return !mod(round(h+k+l),2)
End
//
ThreadSafe Static Function ALLOW_CC(h,k,l)			// C-centered, !mod(round(h+k),2)
	Variable h,k,l
	return !mod(round(h+k),2)
End
//
ThreadSafe Static Function ALLOW_AC(h,k,l)			// A-centered, !mod(round(k+l),2)
	Variable h,k,l
	return !mod(round(k+l),2)
End
//
ThreadSafe Static Function ALLOW_RHOM_HEX(h,k,l)	// rhombohedral hexagonal, allowed are -H+K+L=3n or H-K+L=3n
	Variable h,k,l
	return !mod(-h+k+l,3) || !mod(h-k+l,3)
End
//
ThreadSafe Static Function ALLOW_HEXAGONAL(h,k,l)	// hexagonal, // forbidden are: H+2K=3N with L odd
	Variable h,k,l
	return mod(h+2*k,3) || !mod(l,2)
End



ThreadSafe Function isRhombohedralXtal(xtal)
	STRUCT crystalStructure &xtal
	return (isRhombohedral(xtal.SpaceGroup) && ( abs(xtal.alpha-xtal.beta) + abs(xtal.alpha-xtal.gam) ) < 1e-5)
End

ThreadSafe Static Function isRhombohedral(SpaceGroup)
	Variable SpaceGroup
	switch(SpaceGroup)		// there are only 7 space groups that are Rhombohedral, sym for these all start with "R"
		case 146:
		case 148:
		case 155:
		case 160:
		case 161:
		case 166:
		case 167:
			return 1
	endswitch
	return 0
End



ThreadSafe Function PrimitiveCellFactor(xtal)		// number of primitive unit cells in conventional cell, or number of atoms in conventional cell
	STRUCT crystalStructure &xtal

	String sym = getHMsym(xtal.SpaceGroupIDnum)	// symmetry string
	strswitch(sym[0,0])
		case "F":										// Face Centered
			return 4
		case "I":										// Body Centered
		case "C":										// Side Centered
		case "A":										// Side Centered
		case "B":										// Side Centered
			return 2
		case "R":										// Rhombohedral
			if (abs(xtal.alpha-90)+abs(xtal.beta-90)+abs(xtal.gam-120) < 0.01)		// using Hexagonal axes
				return 3									// there are 3 rhombohedral cells / hexagonal cell
			else
				return 1
			endif
	endswitch

//	Variable SG = xtal.SpaceGroup
//	if (SG<=194 && SG>=168)						// Hexagonal
//		return 1											// Hexagonal is a primitive cell
//	elseif (SG<168 && SG>=143)					// Trigonal (Rhombohedral)
//		if (abs(xtal.alpha-90)+abs(xtal.beta-90)+abs(xtal.gam-120) < 0.01)		// using Hexagonal axes
//			return 3										// there are 3 rhombohedral cells / hexagonal cell
//		endif
//	endif
	return 1
End


ThreadSafe Function isValidLatticeConstants(xtal)	// returns 1 if lattice constants have valid symmetry for the SpaceGroup
	STRUCT crystalStructure &xtal							// this sruct is set in this routine
	if (!isValidSpaceGroupID(xtal.SpaceGroupID))		// check if Space Group ID is valid
		return 0
	endif

	Variable a=xtal.a, b=xtal.b, c=xtal.c
	Variable alpha=xtal.alpha, bet=xtal.beta, gam=xtal.gam
	Variable SG = str2num(xtal.SpaceGroupID)

	switch(latticeSystem(xtal.SpaceGroupID))
		case CUBIC:									// Cubic space groups [195,230]
			return isCUBIC_LC(a,b,c,alpha,bet,gam)

		case HEXAGONAL:							// Hexagonal space groups [168,194]
			return isHEXAGONAL_LC(a,b,c,alpha,bet,gam)

		case TRIGONAL:								// Trigonal space groups [143,167]
			// generally use hexagonal cell, for rhomohedral use rhomohedral cell, unless obviously the hexagonal
			if (isHEXAGONAL_LC(a,b,c,alpha,bet,gam))	// hexagonal is always permitted
				return 1
			elseif (isRhombohedral(SG))		// a rhomohedral space group, may have rhom axes
				return (isRHOMBOHEDRAL_LC(a,b,c,alpha,bet,gam))
			endif
			return 0

		case TETRAGONAL	:						// Tetragonal space groups [75,142]
			return isTETRAGONAL_LC(a,b,c,alpha,bet,gam)

		case ORTHORHOMBIC:						// Orthorhombic space groups [16,74]
			return isORTHORHOMBIC_LC(a,b,c,alpha,bet,gam)

		case MONOCLINIC:							// Monoclinic space groups [3,15]
			return isMONOCLINIC_LC(a,b,c,alpha,bet,gam)

		case TRICLINIC:							// Triclinic space groups 1,2]
			return isTRICLINIC_LC(a,b,c,alpha,bet,gam)

		default:
			return 0									// system is invalid, NOT valid
	endswitch
End
//
ThreadSafe Static Function isCUBIC_LC(a,b,c,alpha,bet,gam)	// returns 1 if LC are Cubic
	// requires alpha=beta=gamma=90, and a=b=c
	Variable a,b,c,alpha,bet,gam
	Variable lenTol = max(max(a,b),c) * 1e-4, angleTol=0.00573		// 1e-4 degree

	Variable basic = numtype(a+b+c+alpha+bet+gam)==0 && a>0
	Variable anglesOK = (abs(alpha-90)<angleTol) && (abs(bet-90)<angleTol) && (abs(gam-90)<angleTol)
	return basic & anglesOK && (abs(a-b)<lenTol) && (abs(a-c)<lenTol)
End
//
ThreadSafe Static Function isHEXAGONAL_LC(a,b,c,alpha,bet,gam)	// returns 1 if LC are Hexagonal
	// requires a=b, alpha=beta=90, gamma=120
	Variable a,b,c,alpha,bet,gam
	Variable lenTol = max(max(a,b),c) * 1e-4, angleTol=0.00573		// 1e-4 degree

	Variable basic = numtype(a+b+c+alpha+bet+gam)==0 && a>0 && c>0
	Variable anglesOK = abs(alpha-90)<angletol && abs(bet-90)<angletol && abs(gam-120)<angleTol
	return anglesOK &&  (abs(a-b)<lenTol)
End
//
ThreadSafe Static Function isRHOMBOHEDRAL_LC(a,b,c,alpha,bet,gam)	// returns 1 if LC are Rhombohedral
	// requires a=b=c, alpha=beta=gamma
	Variable a,b,c,alpha,bet,gam
	Variable lenTol = max(max(a,b),c) * 1e-4, angleTol=0.00573		// 1e-4 degree

	Variable basic = numtype(a+b+c+alpha+bet+gam)==0 && a>0 && alpha>0
	Variable anglesOK = abs(alpha-bet)<angleTol && (alpha-gam)<angleTol
	return anglesOK && abs(a-b)<lenTol && abs(a-c)<lenTol
End
//
ThreadSafe Static Function isTETRAGONAL_LC(a,b,c,alpha,bet,gam)	// returns 1 if LC are Tetragonal
	// requires alpha=beta=gamma=90, and either a=b or b=c or a=c
	Variable a,b,c,alpha,bet,gam
	Variable lenTol = max(max(a,b),c) * 1e-4, angleTol=0.00573		// 1e-4 degree

	Variable basic = numtype(a+b+c+alpha+bet+gam)==0 && a>0 && b>0 && c>0
	Variable anglesOK = (abs(alpha-90)<angleTol) && (abs(bet-90)<angleTol) && (abs(gam-90)<angleTol)
	return anglesOK && ( abs(a-b) || abs(a-c) || abs(b-c) )
End
//
ThreadSafe Static Function isORTHORHOMBIC_LC(a,b,c,alpha,bet,gam)	// returns 1 if LC are Orthorhombic
	// requires alpha=beta=gamma=90, no conditions on a,b,c
	Variable a,b,c,alpha,bet,gam
	Variable angleTol=0.00573		// 1e-4 degree

	Variable basic = numtype(a+b+c+alpha+bet+gam)==0 && a>0 && b>0 && c>0
	return (abs(alpha-90)<angleTol) && (abs(bet-90)<angleTol) && (abs(gam-90)<angleTol)
End
//
ThreadSafe Static Function isMONOCLINIC_LC(a,b,c,alpha,bet,gam)	// returns 1 if LC are Monoclinic
	// requires two of alpha,beta,gamma to be 90, no conditions on a,b,c
	Variable a,b,c,alpha,bet,gam
	Variable angleTol=0.00573		// 1e-4 degree

	Variable basic = numtype(a+b+c+alpha+bet+gam)==0 && a>0 && b>0 && c>0 && alpha>0 && bet>0 && gam>0
	Variable num90s = (abs(alpha-90)<angleTol) + (abs(bet-90)<angleTol) + (abs(gam-90)<angleTol)
	return num90s >= 2
End
//
ThreadSafe Static Function isTRICLINIC_LC(a,b,c,alpha,bet,gam)	// returns 1 if LC are Triclinic
	// no conditions, always 1 if no numbers are NaN, Inf, or negative
	Variable a,b,c,alpha,bet,gam
	return numtype(a+b+c+alpha+bet+gam)==0 && a>0 && b>0 && c>0 && alpha>0 && bet>0 && gam>0
End
//Function TestisValidLatticeConstants()	// returns 1 if lattice constants have valid symmetry for the SpaceGroup
//	STRUCT crystalStructure xtal
//	FillCrystalStructDefault(xtal)
//	print isValidLatticeConstants(xtal)
//End


Function CheckDirectRecip(direct,recip)	// returns True if recip & direct go together
	Wave direct,recip
	Variable tol = WaveType(direct) & 0x04 ? 1e-13 : 1e-6
	Variable tol2 = WaveType(recip) & 0x04 ? 1e-13 : 1e-6
	tol = max(tol,tol2)
	MatrixOP/FREE delta = (2*PI * (Inv(direct))^t) - recip
	MatrixOP/FREE err = sum(abs(delta))
	return err[0] < tol
End

//	End of crystal symmetry stuff
// =========================================================================
// =========================================================================

// =========================================================================
// =========================================================================
//	Start of some utility routines

// changes hkl[3] to the lowest order hkl, ignores whether a reflection is allowed, just removes common factors
ThreadSafe Function lowestOrderHKL(h,k,l)
	Variable &h,&k,&l									// these hkl are returned with all common factors removed
	Variable f, i, maxDiv = gcdZero(h,k)		// max possible divisor
	maxDiv = min(maxDiv, gcdZero(h,l))
	maxDiv = min(maxDiv, gcdZero(k,l))
	maxDiv = numtype(maxDiv) ? 1 : maxDiv		// this happens if h,k,l all zero
	for (i=maxDiv,f=1; i>=2; i-=1)				// check all divisors in range [maxDiv, 2]
		if (mod(h,i) || mod(k,i) || mod(l,i))	// i is not a factor of h, k, and l
			continue
		endif
		h /= i
		k /= i
		l /= i
		f *= i
	endfor
	h = h==0 ? 0 : h									// rmove "-0"
	k = k==0 ? 0 : k
	l = l==0 ? 0 : l
	return f
End
//
ThreadSafe Function gcdZero(a,b)		// this is needed because gcd returns NaN when a or b is 0
	Variable a,b
	Variable answer
	if (a==0 && b==0)
		answer = Inf
	elseif (a==0)
		answer = b
	elseif (b==0)
		answer = a
	else
		answer = gcd(a,b)						// note, gcd ignores the sign
	endif
	return abs(answer)
End
//ThreadSafe Function lowestOrderHKL(h,k,l)
//	Variable &h,&k,&l							// these hkl are returned with all common factors removed
//	Variable maxDiv							// max possible divisor
//	Variable i
//	maxDiv = max(abs(h),abs(k))			// the maximum divisor cannot be bigger than the smallest of hkl
//	maxDiv = max(maxDiv,abs(l))
//	for (i=maxDiv;i>=2;i-=1)				// check all divisorts in range [2, maxDiv]
//		if (mod(h,i) || mod(k,i) || mod(l,i))	// i is not a factor of h, k, and l
//			continue
//		endif
//		h /= i
//		k /= i
//		l /= i
//	endfor
//End


// changes hkl[3] to the lowest order allowed hkl (ie for FCC, 0,0,12 -> 002 not 001
Function lowestAllowedHKL(h,k,l)
	Variable &h,&k,&l

	STRUCT crystalStructure xtal		// temporary crystal structure
	FillCrystalStructDefault(xtal)	//fill the lattice structure with default values
//	ForceLatticeToStructure(xtal)

	Variable i
	Variable hh=h, kk=k, ll=l
	lowestOrderHKL(hh,kk,ll)			// remove all common factors

	for (i=1;i<16;i+=1)					// never need more than 16 to reach an allowed reflection
		h = i*hh								// try each of the multiples to reach an allowed reflection
		k = i*kk
		l = i*ll
		if (allowedHKL(h,k,l,xtal))
			return 0
		endif
	endfor
	return 0
End


Function NearestAllowedHKL(xtal,hkl,[dhklMax])
	// find the hkl of the allowed reflection that is closest to the given hkl
	STRUCT crystalStructure &xtal
	Wave hkl				// wave with input hkl, returned as nearest allowed hkl
	Variable dhklMax	// distance in hkl to search around hkl input
	dhklMax = ParamIsDefault(dhklMax) || numtype(dhklMax) || dhklMax<1 ? 1 : dhklMax

	Variable N=(2*dhklMax+1)^3			// number of hkl's to check
	Variable h0=round(hkl[0]), k0=round(hkl[1]), l0=round(hkl[2])
	if (N<=1)
		hkl={h0,k0,l0}
	endif
	Make/N=(N,3)/FREE/I hklTest
	Variable i, m, dh,dk,dl
	for (dl=0; dl<=dhklMax; dl=incrementIndex(dl))
		for (dk=0; dk<=dhklMax; dk=incrementIndex(dk))
			for (dh=0; dh<=dhklMax; dh=incrementIndex(dh), m+=1)
				hklTest[m][0] = dh+h0		// list of possibel nearest allowed hkl's
				hklTest[m][1] = dk+k0
				hklTest[m][2] = dl+l0
			endfor
		endfor
	endfor
	Wave recip = recipFrom_xtal(xtal)	// reciprocal lattice vectors
	MatrixOp/FREE qIn = recip x hkl		// q of input hkl
	MatrixOp/FREE dqs = sumRows(magSqr((recip x (hklTest)^t)^t - RowRepeat(qIn,N)))
	Make/N=(N)/FREE/I indexW
	MakeIndex dqs, indexW					// indexW contains hkl's ordered by closeness to hkl

	// find first allowed reflection checking in the order given by indexW
	for (i=0;i<N;i+=1)
		m = indexW[i]
		if (allowedHKL(hklTest[m][0],hklTest[m][1],hklTest[m][2],xtal))
			hkl = hklTest[m][p]				// found an allowed reflection, the answer
			return 0
		endif
	endfor
	dhklMax += 2								// increase search range and try again
	return NearestAllowedHKL(xtal,hkl,dhklMax=dhklMax)
End
//	Function test()
//		STRUCT crystalStructure xtal
//		FillCrystalStructDefault(xtal)
//	
//		Make/N=3/D/FREE hkl={0,0,1.45}
//		printf "starting at hkl = %s\r",vec2str(hkl)
//		Variable err = NearestAllowedHKL(xtal,hkl)
//		if (err)
//			print "ERROR"
//		else
//			printf "nearest allowed hkl = %s\r",vec2str(hkl)
//		endif
//		return 0
//	End


Function/WAVE recipFrom_xtal(xtal)					// returns a FREE wave with reciprocal lattice
	STRUCT crystalStructure &xtal
	Make/N=(3,3)/D/FREE RL								// the reciprocal lattice
	RL = { {xtal.as0,xtal.as1,xtal.as2}, {xtal.bs0,xtal.bs1,xtal.bs2}, {xtal.cs0,xtal.cs1,xtal.cs2} }
	if (numtype(sum(RL)) || WaveMax(RL)==0)		// bad numbers in RL
		setDirectRecip(xtal)							// re-make the as0, as1, ...
		RL = { {xtal.as0,xtal.as1,xtal.as2}, {xtal.bs0,xtal.bs1,xtal.bs2}, {xtal.cs0,xtal.cs1,xtal.cs2} }
	endif
	String wnote="waveClass=directLattice;"
	wnote = ReplaceNumberByKey("SpaceGroup",wnote,xtal.SpaceGroup,"=")
	wnote = ReplaceStringByKey("SpaceGroupID",wnote,xtal.SpaceGroupID,"=")
	String str
	sprintf str, "{%.7g,%.7g,%.7g,%.7g,%.7g,%.7g}", xtal.a,xtal.b,xtal.c,xtal.alpha,xtal.beta,xtal.gam
	wnote = ReplaceStringByKey("latticeParameters",wnote,str,"=")
	Note/K RL, wnote
	return RL
End


Function/WAVE directFrom_xtal(xtal)				// returns a FREE wave with real lattice
	STRUCT crystalStructure &xtal
	Make/N=(3,3)/D/FREE DL								// the direct lattice
	DL = { {xtal.a0,xtal.a1,xtal.a2}, {xtal.b0,xtal.b1,xtal.b2}, {xtal.c0,xtal.c1,xtal.c2} }
	if (numtype(sum(DL)) || WaveMax(DL)==0)		// bad numbers in DL
		setDirectRecip(xtal)							// re-make the a0, a1, ...
		DL = { {xtal.a0,xtal.a1,xtal.a2}, {xtal.b0,xtal.b1,xtal.b2}, {xtal.c0,xtal.c1,xtal.c2} }
	endif
	String wnote="waveClass=directLattice;"
	wnote = ReplaceNumberByKey("SpaceGroup",wnote,xtal.SpaceGroup,"=")
	wnote = ReplaceStringByKey("SpaceGroupID",wnote,xtal.SpaceGroupID,"=")
	String str
	sprintf str, "{%.7g,%.7g,%.7g,%.7g,%.7g,%.7g}", xtal.a,xtal.b,xtal.c,xtal.alpha,xtal.beta,xtal.gam
	wnote = ReplaceStringByKey("latticeParameters",wnote,str,"=")
	Note/K DL, wnote
	return DL
End


ThreadSafe Function/WAVE str2recip(str)		// returns a FREE wave with reciprocal lattice
	String str
	Variable as0,as1,as2,bs0,bs1,bs2,cs0,cs1,cs2
	str = ReplaceString("},{",str,"}{")		// sometimes string is like: "{{1,2,3},{4,5,6},{7,8,9}}"
	sscanf str, "{{%g,%g,%g}{%g,%g,%g}{%g,%g,%g}}",as0,as1,as2,bs0,bs1,bs2,cs0,cs1,cs2
	if (V_flag==9)
		Make/N=(3,3)/D/FREE RL						// the reciprocal lattice
		RL = { {as0,as1,as2}, {bs0,bs1,bs2}, {cs0,cs1,cs2} }
		return RL
	else
		return $""
	endif
End


ThreadSafe Function/S hkl2str(h,k,l, [bar])	// format h,k,l into a string of acceptable minimal length
	Variable h,k,l
	Variable bar														// OPTIONAL, return unicocde string with bars instead of negatives
	bar = ParamIsDefault(bar) || numtype(bar) ? 0 : bar	// default is negative signs
	bar = IgorVersion()<7 ? 0 : bar							// bar not available for Igor 6 (requires unicode support)
	h = abs(h)<1e-14 ? 0 : h
	k = abs(k)<1e-14 ? 0 : k
	l = abs(l)<1e-14 ? 0 : l
	String hklStr
	if (numtype(h+k+l))
		hklStr = "nan,nan,nan"
	elseif (abs(mod(h,1))+abs(mod(k,1))+abs(mod(l,1)) > 1e-6)
		sprintf hklStr,"%g, %g, %g",h,k,l					// hkl are non-integers
	elseif (bar)
		String sp = SelectString(abs(h)<10 && abs(k)<10 && abs(l)<10, " ", "")
		hklStr = minus2bar(num2istr(h)) + sp + minus2bar(num2istr(k)) + sp + minus2bar(num2istr(l))
	elseif (k<0 || l<0)
		sprintf hklStr,"%.0f, %.0f, %.0f",h,k,l
	elseif (abs(h)<10 && k<10 && l<10)
		sprintf hklStr,"%.0f%.0f%.0f",h,k,l
	else
		sprintf hklStr,"%.0f %.0f %.0f",h,k,l
	endif
#if (IgorVersion()>7)
	if (bar)														// only Arial and Tahoma work properly with OVERLINE
		hklStr = "\\[0\\F'"+BAR_FONT_ALWAYS+"'" + hklStr + "\\F]0"
	endif
#endif
	return hklStr
End
//ThreadSafe Function/S hkl2str(h,k,l)	// format h,k,l into a string of acceptable minimal length
//	Variable h,k,l
//	if (numtype(h+k+l))
//		return "nan,nan,nan"
//	endif
//	h = abs(h)<1e-14 ? 0 : h
//	k = abs(k)<1e-14 ? 0 : k
//	l = abs(l)<1e-14 ? 0 : l
//	String hkl
//	if (abs(mod(h,1))+abs(mod(k,1))+abs(mod(l,1)) > 1e-6)	// hkl are non-integers
//		sprintf hkl,"%g, %g, %g",h,k,l
//	elseif (k<0 || l<0)
//		sprintf hkl,"%.0f, %.0f, %.0f",h,k,l
//	elseif (abs(h)<10 && k<10 && l<10)
//		sprintf hkl,"%.0f%.0f%.0f",h,k,l
//	else
//		sprintf hkl,"%.0f %.0f %.0f",h,k,l
//	endif
//	return hkl
//End


ThreadSafe Function str2hkl(hklStr,h,k,l)
	// returns the hkl values from a string, pretty forgiving about format in string
	// moved to here from Dynamical.ipf versions <=1.14
	String hklStr
	Variable &h,&k,&l

	h = NaN ;	k = NaN ;	l = NaN
	hklStr = ReplaceString("+",hklStr," ")
	hklStr = ReplaceString("-",hklStr," -")
	hklStr = TrimFrontBackWhiteSpace(hklStr)

	Wave w3 = str2vec(hklStr)
	Variable N=numpnts(w3)
	if (N==1)
		hklStr[1] = ";"
		hklStr[3] = ";"
		Wave w3 = str2vec(hklStr)
		N=numpnts(w3)
	elseif (N==2)
		String str=hklStr
		str[1] = " "							// split at the start
		Wave w3 = str2vec(str)
		if (numpnts(w3)==2)
			str=hklStr
			str[strlen(hklStr)-1] = " "	// split at the end
			Wave w3 = str2vec(str)
		endif
		N=numpnts(w3)
	endif

	if (N>2 && numtype(sum(w3,0,2))==0)
		h = w3[0]
		k = w3[1]
		l = w3[2]
		return 0
	endif
	return 1
End
//	Function test(str)
//		String str
//		Variable h,k,l
//		Variable err = str2hkl(str,h,k,l)
//		printf "'%s' --> (%g, %g, %g)%s\r",str,h,k,l,SelectString(err,"","   ****** ERROR *****")
//	End


#if (IgorVersion() > 7)
ThreadSafe Function/T hkl2IgorBarStr(h,k,l)	// changes negatives to a bar over the number, only for displaying, not printing
	Variable h,k,l				// hkl value
	if (numtype(h+k+l))
		return num2str(h)+","+num2str(k)+","+num2str(l)
	endif
	//	String sp=SelectString(abs(h)>9 || abs(k)>9 || abs(l)>9 || h<0 || k<0 || l<0,""," ")
	String sp=SelectString(abs(h)>9 || abs(k)>9 || abs(l)>9,""," ")
	//	print h,k,l,"     ",minus2bar(num2istr(h)) + sp + minus2bar(num2istr(k)) + sp + minus2bar(num2istr(l))
	return minus2bar(num2istr(h)) + sp + minus2bar(num2istr(k)) + sp + minus2bar(num2istr(l))
End
//
ThreadSafe Static Function/T minus2bar(str,[single])		// change an Igor string that has minuses to one using a bar over the following character
	String str
	Variable single								// minus sign only applies to next character, otherwise assume integers
	single = ParamIsDefault(single) || numtype(single) ? 0 : single

	String c, out=""
	Variable i
	if (single)
		for (i=0;i<strlen(str);i+=1)
			if (CmpStr(str[i],"-")==0)		// for each "-" put OVERLINE after the following character
				i += 1
				out += str[i]+OVERLINE
			else
				out += str[i]
			endif
		endfor
	else
		Variable ic, neg=0
		for (i=0,neg=0; i<strlen(str); i+=1)
			c = str[i]
			ic = char2num(str[i])
			if (ic==45)							// a negative sign
				neg = 1
			elseif (48<=ic && ic<=57)		// c is a digit {0,1,2,3,4,5,6,7,8,9}
				out += c+SelectString(neg, "", OVERLINE)
			else
				out += c+SelectString(neg && ic>32, "", OVERLINE)
				neg = 0							// not a digit, reset neg and copy char
			endif
		endfor
	endif
	return out
End
#else
ThreadSafe Function/T hkl2IgorBarStr(h,k,l)	// changes negatives to a bar over the number, only for displaying, not printing
	Variable h,k,l				// hkl value

	if (numtype(h+k+l))
		return num2str(h)+","+num2str(k)+","+num2str(l)
	endif
	//	String extra=SelectString(abs(h)>9 || abs(k)>9 || abs(l)>9 || h<0 || k<0 || l<0,""," ")
	String extra=SelectString(abs(h)>9 || abs(k)>9 || abs(l)>9,""," ")
	String str=""
	str += minus2bar(num2istr(h),spaces=floor(log(abs(h)))) + extra
	str += minus2bar(num2istr(k),spaces=floor(log(abs(k)))) + extra
	str += minus2bar(num2istr(l),spaces=floor(log(abs(l))))
	return str
End
//
ThreadSafe Function/T minus2bar(str,[spaces,single])	// change an Igor string that has minuses to one using a bar over the following character
	String str
	Variable spaces
	Variable single								// minus sign only applies to next character, otherwise assume integers
	spaces = ParamIsDefault(spaces) ? 0 : spaces
	spaces = round(spaces)==limit(spaces,1,5) ? spaces : 0
	single = 0										// single is only here for compatibility with the Igor 7 version

	String sspaces = PadString("",spaces,0x20), bs
	sprintf bs, "\\[9\\S\\f01%s\\]9\\M\\X9",sspaces+BCHAR
	str = ReplaceString("-",str,bs)
	return str
End
#endif

//	Function copy_xtal(target,source)						// copy a crystalStructure source to target
//		STRUCT crystalStructure &source
//		STRUCT crystalStructure &target
//	
//		target.desc = source.desc
//	
//		target.a = source.a
//		target.b = source.b
//		target.c = source.c
//		target.alpha = source.alpha
//		target.beta = source.beta
//		target.gam = source.gam
//	
//		target.SpaceGroup = source.SpaceGroup
//		target.SpaceGroupID = source.SpaceGroupID
//		target.SpaceGroupIDnum = source.SpaceGroupIDnum
//		target.Vc = source.Vc
//		target.density = source.density
//		target.alphaT = source.alphaT
//	
//		target.Temperature = source.Temperature
//		target.Vibrate = source.Vibrate
//		target.haveDebyeT = source.haveDebyeT
//		target.hashID = source.hashID
//	
//		target.N = source.N
//		Variable i, N=source.N
//		for (i=0;i<N;i+=1)
//			copy_atomType(target.atom[i],source.atom[i])
//		endfor
//	
//		target.Nbonds = source.Nbonds
//		N = target.Nbonds
//		for (i=0;i<N;i+=1)
//			copy_bondType(target.bond[i],source.bond[i])
//		endfor
//	
//		target.a0 = source.a0
//		target.b0 = source.b0
//		target.c0 = source.c0
//		target.a1 = source.a1
//		target.b1 = source.b1
//		target.c1 = source.c1
//		target.a2 = source.a2
//		target.b2 = source.b2
//		target.c2 = source.c2
//	
//		target.as0 = source.as0
//		target.bs0 = source.bs0
//		target.cs0 = source.cs0
//		target.as1 = source.as1
//		target.bs1 = source.bs1
//		target.cs1 = source.cs1
//		target.as2 = source.as2
//		target.bs2 = source.bs2
//		target.cs2 = source.cs2
//	
//		target.Unconventional00 = source.Unconventional00
//		target.Unconventional01 = source.Unconventional01
//		target.Unconventional02 = source.Unconventional02
//		target.Unconventional10 = source.Unconventional10
//		target.Unconventional11 = source.Unconventional11
//		target.Unconventional12 = source.Unconventional12
//		target.Unconventional20 = source.Unconventional20
//		target.Unconventional21 = source.Unconventional21
//		target.Unconventional22 = source.Unconventional22
//	
//		String fullFile = source.sourceFile
//		fullFile = fullFile[0,MAX_FILE_LEN-1]
//		target.sourceFile = fullFile
//	End
//
//	Static Function copy_atomType(target,source)		// copy a atomTypeStructure source to target
//		STRUCT atomTypeStructure &source
//		STRUCT atomTypeStructure &target
//	
//		target.name = source.name
//		target.Zatom = source.Zatom
//		target.valence = source.valence
//		target.x = source.x
//		target.y = source.y
//		target.z = source.z
//		target.occ = source.occ
//		target.WyckoffSymbol = source.WyckoffSymbol
//		target.DebyeT = source.DebyeT
//		target.Biso = source.Biso
//		target.Uiso = source.Uiso
//		target.U11 = source.U11
//		target.U22 = source.U22
//		target.U33 = source.U33
//		target.U12 = source.U12
//		target.U13 = source.U13
//		target.U23 = source.U23
//	End
//
//	Static Function copy_bondType(target,source)		// copy a bondTypeStructure source to target
//		STRUCT bondTypeStructure &source
//		STRUCT bondTypeStructure &target
//	
//		target.label0 = source.label0
//		target.label1 = source.label1
//		target.N = source.N
//		Variable i,N = source.N
//		for (i=0;i<N;i+=1)
//			target.len[i] = source.len[i]
//		endfor
//	End


Function StructGet_xtal(strStruct,xtal6)		// take value of strStruct, and fill xtal6, whether strStruct was v5 or v6
	// this accepts strStruct for either the old v5 or new v6 crystalStructure
	String strStruct										// pre vers 6.0 len=25672, after 6.0 it is 25686
	STRUCT crystalStructure &xtal6					// structure that will be filled

	if (strlen(strStruct)>=25686)					// string is for a version 6 xtal
		StructGet/S/B=2 xtal6, strStruct
	else														// shorter strings assumed to be older version 5 xtal
		STRUCT crystalStructure5 xtal5
		StructGet/S/B=2 xtal5, strStruct
		copy_xtal56(xtal6,xtal5)						// copy xtal5 --> xtal
	endif
End
//
Static Function copy_xtal56(xtal6,xtal5)					// copy a crystalStructure xtal5 --> xtal6
	STRUCT crystalStructure &xtal6
	STRUCT crystalStructure5 &xtal5

	xtal6.desc = xtal5.desc
	xtal6.a = xtal5.a				;	xtal6.b = xtal5.b				;	xtal6.c = xtal5.c
	xtal6.alpha = xtal5.alpha	;	xtal6.beta = xtal5.beta	;	xtal6.gam = xtal5.gam
	xtal6.SpaceGroup = xtal5.SpaceGroup						// in range [1,230]

	String id = FindDefaultIDforSG(xtal5.SpaceGroup)
	xtal6.SpaceGroupID = id										// change SG number to id string
	xtal6.SpaceGroupIDnum = SpaceGroupID2num(id)			// change id to id number in [1,530]

	xtal6.Vc = xtal5.Vc
	xtal6.density = xtal5.density
	xtal6.alphaT = xtal5.alphaT

	xtal6.Temperature = xtal5.Temperature
	xtal6.Vibrate = xtal5.Vibrate
	xtal6.haveDebyeT = xtal5.haveDebyeT
	xtal6.hashID = xtal5.hashID

	xtal6.N = xtal5.N
	Variable i, N=xtal5.N
	for (i=0;i<N;i+=1)
		xtal6.atom[i] = xtal5.atom[i]							// was:  copy_atomType(xtal6.atom[i],xtal5.atom[i])
	endfor

	xtal6.Nbonds = xtal5.Nbonds
	N = xtal6.Nbonds
	for (i=0;i<N;i+=1)
		xtal6.bond[i] = xtal5.bond[i]							// was:  copy_bondType(xtal6.bond[i],xtal5.bond[i])
	endfor

	xtal6.a0 = xtal5.a0	;	xtal6.b0 = xtal5.b0	;	xtal6.c0 = xtal5.c0
	xtal6.a1 = xtal5.a1	;	xtal6.b1 = xtal5.b1	;	xtal6.c1 = xtal5.c1
	xtal6.a2 = xtal5.a2	;	xtal6.b2 = xtal5.b2	;	xtal6.c2 = xtal5.c2

	xtal6.as0 = xtal5.as0	;	xtal6.bs0 = xtal5.bs0	;	xtal6.cs0 = xtal5.cs0
	xtal6.as1 = xtal5.as1	;	xtal6.bs1 = xtal5.bs1	;	xtal6.cs1 = xtal5.cs1
	xtal6.as2 = xtal5.as2	;	xtal6.bs2 = xtal5.bs2	;	xtal6.cs2 = xtal5.cs2

	xtal6.Unconventional00 = xtal5.Unconventional00
	xtal6.Unconventional01 = xtal5.Unconventional01
	xtal6.Unconventional02 = xtal5.Unconventional02
	xtal6.Unconventional10 = xtal5.Unconventional10
	xtal6.Unconventional11 = xtal5.Unconventional11
	xtal6.Unconventional12 = xtal5.Unconventional12
	xtal6.Unconventional20 = xtal5.Unconventional20
	xtal6.Unconventional21 = xtal5.Unconventional21
	xtal6.Unconventional22 = xtal5.Unconventional22

	String fullFile = xtal5.sourceFile
	fullFile = fullFile[0,MAX_FILE_LEN-1]
	xtal6.sourceFile = fullFile
End

//	End of utility
// =========================================================================
// =========================================================================


//
//	================================================================================
//
//			This section creates the wave definitions used in    InitLatticeSymPackage()
//



Function InitLatticeSymPackage([showPanel])			// used to initialize this package
	Variable showPanel											// optionally show the Panel too
	showPanel = ParamIsDefault(showPanel) || numtype(showPanel) ? 0 : showPanel
	NewDataFolder/O root:Packages
	NewDataFolder/O root:Packages:Lattices
	NewDataFolder/O root:Packages:Lattices:SymOps
	if (!exists("root:Packages:Lattices:keV"))
		Variable/G root:Packages:Lattices:keV = 10		// only used when calculating Cromer-Liberman values
	endif
	if (showPanel)
		MakeLatticeParametersPanel("")
	endif
End



//
//	================================================================================
//
//			This section is for making the matricies that changes the settings for a Space Group
//

Static Function/WAVE GetSettingTransForm(id)
	String id					// SpaceGroup ID name, e.g. "146:H"

	Variable i = SpaceGroupID2num(id)		// in range [1,530]
	if (numtype(i))								// check for valid id
		return $""
	endif

	Make/N=(530)/T/FREE CBMs=""				// the transform from id to the "standard" setting
	CBMs[0]   = {"x,y,z","x,y,z","x,y,z","y,z,x","z,x,y","x,y,z","y,z,x","z,x,y","x,y,z","-x+z,y,-x","-z,y,x-z","y,z,x","x-y,z,-y","-x,z,-x+y","z,x,y"}
	CBMs[15]  = {"y-z,x,-z","-y,x,-y+z","x,y,z","y,z,x","z,x,y","x,y,z","-x+z,y,-x","-z,y,x-z","y,z,x","x-y,z,-y","-x,z,-x+y","z,x,y","y-z,x,-z"}
	CBMs[28]  = {"-y,x,-y+z","x,y,z","-x+z,y,-x","-z,y,x-z","y,z,x","x-y,z,-y","-x,z,-x+y","z,x,y","y-z,x,-z","-y,x,-y+z","x,y,z","-x+z,y,-x"}
	CBMs[40]  = {"-z,y,x-z","-x+z,y+1/4,-x","x,y-1/4,z","-z,y+1/4,x-z","y,z,x","x-y,z,-y","-x,z,-x+y","x-y,z+1/4,-y","y,z-1/4,x","-x,z+1/4,-x+y"}
	CBMs[50]  = {"z,x,y","y-z,x,-z","-y,x,-y+z","y-z,x+1/4,-z","z,x+1/4,y","-y,x+1/4,-y+z","x,y,z","y,z,x","z,x,y","x,y,z","y,z,x","z,x,y","x,y,z"}
	CBMs[63]  = {"-x+z,y,-x","-z,y,x-z","y,z,x","x-y,z,-y","-x,z,-x+y","z,x,y","y-z,x,-z","-y,x,-y+z","x,y,z","-x+z,y,-x","-z,y,x-z","y,z,x"}
	CBMs[75]  = {"x-y,z,-y","-x,z,-x+y","z,x,y","y-z,x,-z","-y,x,-y+z","x,y,z","-x+z,y,-x","-z,y,x-z","y,z,x","x-y,z,-y","-x,z,-x+y","z,x,y"}
	CBMs[87]  = {"y-z,x,-z","-y,x,-y+z","x,y,z","-x+z,y,-x","-z,y,x-z","-x+z-1/4,y+1/4,-x","x+1/4,y+1/4,z","-z-1/4,y+1/4,x-z","y,z,x","x-y,z,-y"}
	CBMs[97]  = {"-x,z,-x+y","x-y-1/4,z+1/4,-y","y+1/4,z+1/4,x","-x-1/4,z+1/4,-x+y","z,x,y","y-z,x,-z","-y,x,-y+z","y-z-1/4,x+1/4,-z","z-1/4,x+1/4,y"}
	CBMs[106] = {"-y-1/4,x+1/4,-y+z","x,y,z","x,y,z","y,z,x","z,x,y","x,y,z","y,z,x","z,x,y","x,y,z","x,y,z","y,z,x","z,x,y","x,y,z","y,z,x","z,x,y"}
	CBMs[121] = {"x,y,z","x,y,z","x,y,z","x,y,z","y,z,x","z,x,y","x,y,z","y,x,-z","y,z,x","z,y,-x","z,x,y","x,z,-y","x,y,z","y,z,x","z,x,y","x,y,z"}
	CBMs[137] = {"y,x,-z","y,z,x","z,y,-x","z,x,y","x,z,-y","x,y,z","y,x,-z","y,z,x","z,y,-x","z,x,y","x,z,-y","x,y,z","y,x,-z","y,z,x","z,y,-x"}
	CBMs[152] = {"z,x,y","x,z,-y","x,y,z","y,x,-z","y,z,x","z,y,-x","z,x,y","x,z,-y","x,y,z","y,z,x","z,x,y","x,y,z","y,x,-z","y,z,x","z,y,-x"}
	CBMs[167] = {"z,x,y","x,z,-y","x,y,z","y,z,x","z,x,y","x,y,z","y,z,x","z,x,y","x,y,z","y,x,-z","y,z,x","z,y,-x","z,x,y","x,z,-y","x,y,z","y,z,x"}
	CBMs[183] = {"z,x,y","x,y,z","y,x,-z","y,z,x","z,y,-x","z,x,y","x,z,-y","x,y,z","y,x,-z","y,z,x","z,y,-x","z,x,y","x,z,-y","x,y,z","y,x,-z"}
	CBMs[198] = {"y,z,x","z,y,-x","z,x,y","x,z,-y","x,y,z","y,x,-z","y,z,x","z,y,-x","z,x,y","x,z,-y","x,y,z","y,z,x","z,x,y","x,y,z","y,z,x"}
	CBMs[213] = {"z,x,y","x,y,z","y,z,x","z,x,y","x,y,z","y,z,x","z,x,y","x,y,z","y,x,-z","y,z,x","z,y,-x","z,x,y","x,z,-y","x,y,z","x,y,z"}
	CBMs[228] = {"x-1/4,y-1/4,z-1/4","x,y,z","y,z,x","z,x,y","x,y,z","x-1/4,y-1/4,z","y,z,x","y-1/4,z-1/4,x","z,x,y","z-1/4,x-1/4,y","x,y,z"}
	CBMs[239] = {"y,x,-z","y,z,x","z,y,-x","z,x,y","x,z,-y","x,y,z","y,x,-z","y,z,x","z,y,-x","z,x,y","x,z,-y","x,y,z","y,x,-z","y,z,x","z,y,-x"}
	CBMs[254] = {"z,x,y","x,z,-y","x,y,z","y,x,-z","y,z,x","z,y,-x","z,x,y","x,z,-y","x,y,z","y,z,x","z,x,y","x,y,z","y,z,x","z,x,y","x,y,z"}
	CBMs[269] = {"y,x,-z","y,z,x","z,y,-x","z,x,y","x,z,-y","x,y,z","y,z,x","z,x,y","x,y,z","x-1/4,y-1/4,z","y,z,x","y-1/4,z-1/4,x","z,x,y"}
	CBMs[282] = {"z-1/4,x-1/4,y","x,y,z","y,x,-z","y,z,x","z,y,-x","z,x,y","x,z,-y","x,y,z","y,x,-z","x,y,z","y,x,-z","y,z,x","z,y,-x","z,x,y"}
	CBMs[296] = {"x,z,-y","x,y,z","y,x,-z","y,z,x","z,y,-x","z,x,y","x,z,-y","x,y,z","y,x,-z","y,z,x","z,y,-x","z,x,y","x,z,-y","x,y,z","y,z,x"}
	CBMs[311] = {"z,x,y","x,y,z","y,z,x","z,x,y","x,y,z","x+1/4,y+1/4,z","y,z,x","y+1/4,z+1/4,x","z,x,y","z-1/4,x+1/4,y","x,y,z","x,y-1/4,z-1/4"}
	CBMs[323] = {"x,y,z","x+1/4,y,z-1/4","y,z,x","y,z-1/4,x-1/4","y,z,x","y+1/4,z,x-1/4","z,x,y","z,x+1/4,y-1/4","z,x,y","z-1/4,x,y-1/4","x,y,z"}
	CBMs[334] = {"x,y,z","x+1/8,y+1/8,z+1/8","x,y,z","x,y,z","y,z,x","z,x,y","x,y,z","x-1/4,y+1/4,z+1/4","x,y,z","x-1/4,y+1/4,z+1/4","y,z,x"}
	CBMs[345] = {"y-1/4,z+1/4,x+1/4","z,x,y","z-1/4,x+1/4,y+1/4","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z"}
	CBMs[358] = {"x,y,z","x-1/4,y+1/4,z","x,y,z","x+1/4,y+1/4,z+1/4","x,y,z","x,y,z","x+1/2,y-1/4,z+1/8","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z"}
	CBMs[370] = {"x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z"}
	CBMs[386] = {"x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z"}
	CBMs[402] = {"x-1/4,y-1/4,z","x,y,z","x-1/4,y-1/4,z-1/4","x,y,z","x,y,z","x,y,z","x-1/4,y+1/4,z","x,y,z","x-1/4,y+1/4,z","x,y,z","x,y,z","x,y,z"}
	CBMs[414] = {"x-1/4,y+1/4,z+1/4","x,y,z","x-1/4,y+1/4,z-1/4","x,y,z","x,y,z","x,y,z","x-1/4,y+1/4,z+1/4","x,y,z","x-1/4,y+1/4,z-1/4","x,y,z"}
	CBMs[424] = {"x,y,z","x,y,z","x,y-1/4,z+1/8","x,y,z","x,y-1/4,z+1/8","x,y,z","x,y,z","x,y,z","x,y,z"}
	CBMs[433] = {"2/3*x-1/3*y-1/3*z,1/3*x+1/3*y-2/3*z,1/3*x+1/3*y+1/3*z","x,y,z","x,y,z","2/3*x-1/3*y-1/3*z,1/3*x+1/3*y-2/3*z,1/3*x+1/3*y+1/3*z"}
	CBMs[437] = {"x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","2/3*x-1/3*y-1/3*z,1/3*x+1/3*y-2/3*z,1/3*x+1/3*y+1/3*z","x,y,z","x,y,z"}
	CBMs[447] = {"x,y,z","x,y,z","x,y,z","2/3*x-1/3*y-1/3*z,1/3*x+1/3*y-2/3*z,1/3*x+1/3*y+1/3*z","x,y,z"}
	CBMs[452] = {"2/3*x-1/3*y-1/3*z,1/3*x+1/3*y-2/3*z,1/3*x+1/3*y+1/3*z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z"}
	CBMs[458] = {"2/3*x-1/3*y-1/3*z,1/3*x+1/3*y-2/3*z,1/3*x+1/3*y+1/3*z","x,y,z","2/3*x-1/3*y-1/3*z,1/3*x+1/3*y-2/3*z,1/3*x+1/3*y+1/3*z","x,y,z"}
	CBMs[462] = {"x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z"}
	CBMs[478] = {"x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z"}
	CBMs[494] = {"x,y,z","x-1/4,y-1/4,z-1/4","x,y,z","x,y,z","x+1/8,y+1/8,z+1/8","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z"}
	CBMs[507] = {"x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x,y,z","x-1/4,y-1/4,z-1/4","x,y,z","x,y,z"}
	CBMs[521] = {"x+1/4,y+1/4,z+1/4","x,y,z","x,y,z","x,y,z","x+1/8,y+1/8,z+1/8","x,y,z","x+3/8,y+3/8,z+3/8","x,y,z","x,y,z"}

	String CBMx = CBMs[i-1]
	Wave mat = MatrixFromSymLine(CBMx,3,zeroBad=1)		// (3,3) matrix, NONE of the CBM have a constant part, so a (4,4) mat is pointless
	String wnote="waveClass=ChangeBasisMat"					// fill the wave note
	wnote = ReplaceStringByKey("SGid",wnote,id,"=")
	wnote = ReplaceStringByKey("CBMx",wnote,CBMx,"=")
	Note/K mat, wnote
	return mat
End
//
//	Function test_GetSettingTransForm()
//		String id, idList="146:R;146:H;346;-1"
//		Variable i, N=ItemsInList(idList)
//		for (i=0;i<N;i+=1)
//			id = StringFromList(i,idList)
//			Wave mat = GetSettingTransForm(id)
//			print " "
//			if (WaveExists(mat))
//				print note(mat)
//				printWave(mat,name="mat",brief=1)
//			else
//				printf "no settings transform matrix with Space Group ID = \"%s\"\r",id
//			endif
//		endfor
//	End
//
// makes one matrix either (3,3) or (4,4) from the symmetry op string
Static Function/WAVE MatrixFromSymLine(symOp,dim,[zeroBad])	// returns either 3x3 or 4x4 matrix
	String symOp													// something like "-x+1/2,-y+1/2,z" or "2/3*x-1/3*y-1/3*z,1/3*x+1/3*y-2/3*z,1/3*x+1/3*y+1/3*z"
	Variable dim													// either 3 or 4
	Variable zeroBad												// if 1, do not allow |mat3| near zero
	zeroBad = ParamIsDefault(zeroBad) || numtype(zeroBad) ? 0 : zeroBad

	if (dim!=3 && dim!=4)
		return $""
	endif
	Make/N=(dim,dim)/D/FREE mat=0							// either (3,3) or (4,4)

	String expression
	Variable i, err
	for (i=0,err=0; i<3; i+=1)
		expression = StringFromList(i,symOp,",")

		Wave vec = ParseOneSymExpression(expression)	// parse one expression of form "-x+y"  or "-x", or "-x+y", etc.
		if (!WaveExists(vec))
			err = 1
			break
		endif
		mat[i][] = vec[q]											// vec is always [4], but this still works for dim=3
	endfor

	if (!err && zeroBad)										// skip this if !zeroBad
		Make/N=(3,3)/D/FREE mat3 = mat[p][q]				// mat3 is ALWAYS (3,3)
		err = err || (MatrixDet(mat3)<0.001)
	endif
	if (err)
		WaveClear mat
	endif
	return mat
End
//
Static Function/WAVE ParseOneSymExpression(expression)	// parse one expression of form "-x+y"  or "-x", or "-x+y, etc.
	String expression

	expression = LowerStr(expression)		// only lower case
	expression = ReplaceString(" ",expression,"")	// no spaces
	expression = ReplaceString("-",expression,"+-")	// so all terms are joined by a "+"
	do
		expression = ReplaceString("++",expression,"+")	// change all multiple "+" to a single "+"
	while (strsearch(expression,"++",0)>=0)
	expression = TrimFront(expression,chars="+")	// remove any leading "+"

	Make/N=4/D/FREE vec=0
	String term
	Variable i=0, N=ItemsInList(expression,"+")
	do
		term = StringFromList(i,expression,"+")
		if (strlen(term)<1)
			break
		endif
		vec[0] += EvaluateOneVariable(term,"x")
		vec[1] += EvaluateOneVariable(term,"y")
		vec[2] += EvaluateOneVariable(term,"z")
		vec[3] += EvaluateOneVariable(term,"")
		i += 1
	while (1)
	return vec
End
//
Static Function EvaluateOneVariable(term,varName)
	String term		// an algebraic term, does not contain "+", and the "-" can ONLY be at the front
	String varName	// name of variable, e.g. "x", this NOT case sensitive

	Variable val=0
	term = ReplaceString(" ",term,"")
	if (strlen(term)==0)						// empty
		val = 0
	elseif (strlen(varName)==0)				// no variable, the constant term
		val = arithmetic(term,def=0)
	elseif (strsearch(term,varName,0,2)>=0)
		term = ReplaceString(varName,term,"*1*")
		do												// change all double "**" to single "*"
			term = ReplaceString("**",term,"*")
		while(strsearch(term,"**",0)>=0)
		term = TrimBoth(term, chars="*")	// remove any leading or trailing "*"
		term = ReplaceString("-*",term,"-")
		val = arithmetic(term,def=1)
	endif
	return val
End
//
//	Function test_ParseOneSymExpression()
//		String list = " z*2; z2; 2*z; 2z; 2/3*x-1/3*y-1/3*z + -1/7; 2/3*x-1/3*y-1/3*z -1/7; 2/3*x-1/3*y-1/3*z +1/7; 2/3*x-1/3*y; 2/3*x-1/3*y-1/3*z + -1/7 + 1/7"
//		Make/N=(4,9)/D/FREE answers
//		answers[0][0] = {0,0,2,0}	;	answers[0][1] = {0,0,2,0}	;	answers[0][2] = {0,0,2,0}
//		answers[0][3] = {0,0,2,0}	;					answers[0][4] = {2/3,-1/3,-1/3,-1/7}
//		answers[0][5] = {2/3,-1/3,-1/3,-1/7}	;	answers[0][6] = {2/3,-1/3,-1/3,1/7}
//		answers[0][7] = {2/3,-1/3,0,0}			;	answers[0][8] = {2/3,-1/3,-1/3,0}
//		Make/N=4/D/FREE answer
//		String expression
//		Variable i, bad
//		for (i=0; i<ItemsInList(list); i+=1)
//			expression = StringFromList(i,list)
//			Wave vec = LatticeSym#ParseOneSymExpression(expression)
//			answer = answers[p][i]
//			MatrixOp/FREE diff = Sum(Abs(vec-answer))
//			bad = abs(diff[0]>1e-8)
//			printf "%s\t%d\t%s  -->  %s",SelectString(bad,"\t\t","ERROR"),i,expression,vec2str(vec)
//			if (bad)
//				printf "\t\tanswer = %s",vec2str(answer)
//			endif
//			printf "\r"
//		endfor
//	End









//
//	================================================================================
//
//			This section is for making the matricies for the symmetry operations for Space Groups
//

Static Function SetSymOpsForSpaceGroup(SpaceGroupID)	// make the symmetry operations mats and vecs (if needed), returns number of operations
	String SpaceGroupID
	Variable SG=str2num(SpaceGroupID)
	Wave mats = $("root:Packages:Lattices:SymOps:"+CleanupName("equivXYZM"+SpaceGroupID,0))
	Wave bvecs = $("root:Packages:Lattices:SymOps:"+CleanupName("equivXYZB"+SpaceGroupID,0))
	Variable numSymOps
	if (WaveExists(mats) && WaveExists(bvecs))				// check if they exist
		if (StringMatch(StringByKey("SpaceGroupID",note(mats),"="), SpaceGroupID))
			numSymOps = NumberByKey("numSymOps",note(mats),"=")
			return numtype(numSymOps) ? 0 : numSymOps		// do not re-make, just return number of operations
		endif
	endif
	if (!isValidSpaceGroup(SG))								// Space Group must be in range [1, 230]
		DoAlert 0, "Bad Space Group = "+SpaceGroupID+", in SetSymOpsForSpaceGroup"
		return 1
	endif
	if (!DataFolderExists("root:Packages:Lattices:SymOps:"))
		DoAlert 0, "Cannot make symmetry matricies, the target data folder does not exist"
		print "Cannot make symmetry matricies, the target data folder does not exist,  'root:Packages:Lattices:SymOps:'\r"
		return 0
	endif

	String symOperations=setSymLineID(SpaceGroupID)	// a string like "x,y,z;-x,-y,z;-x,y,-z;x,-y,-z;x+1/2,y+1/2,z;-x+1/2,-y+1/2,z;-x+1/2,y+1/2,-z;x+1/2,-y+1/2,-z"

	Variable i,N=ItemsInList(symOperations)
	String wName = "root:Packages:Lattices:SymOps:"+CleanupName("equivXYZM"+SpaceGroupID,0)
	Make/N=(N,3,3)/O/B $wName									// this only holds 0 or ±1
	Wave equivM = $wName
	wName = "root:Packages:Lattices:SymOps:"+CleanupName("equivXYZB"+SpaceGroupID,0)
	Make/N=(N,3)/O/D $wName
	Wave equivB = $wName

	Make/N=(3,3)/O/D mat_SymItem
	Make/N=3/O/D vec_SymItem
	Wave mat=mat_SymItem, vec=vec_SymItem
	Variable err = 0
	for (i=0;i<N;i+=1)
		err = err || make1MatrixAndVecFromSymLine(StringFromList(i,symOperations))
		equivM[i][][] = mat[q][r]
		equivB[i][] = vec[q]
	endfor
	KillWaves/Z mat_SymItem,vec_SymItem
	if (err)
		Abort "error making symmetry matricies in SetSymOpsForSpaceGroup()"
	endif

	String wnote="waveClass=SymmetryOperations;"
	wnote = ReplaceNumberByKey("numSymOps",wnote,N,"=")
	wnote = ReplaceNumberByKey("SpaceGroup",wnote,SG,"=")
	wnote = ReplaceStringByKey("SpaceGroupID",wnote,SpaceGroupID,"=")
	Note/K equivM, wnote
	Note/K equivB, wnote
	return N
End
//
// makes one matrix[3][3] and one vector[3] from the expression like "x+1/3,y+2/3,z+2/3"
Static Function make1MatrixAndVecFromSymLine(symItem)		// returns result in mat_SymItem and vec_SymItem, return value is error flag
	String symItem													// something like "-x+1/2,-y+1/2,z"
	Wave mat=mat_SymItem, vec=vec_SymItem
	if (!WaveExists(mat) || !WaveExists(vec))
		Make/N=(3,3)/O/D mat_SymItem
		Make/N=3/O/D vec_SymItem
		Wave mat=mat_SymItem, vec=vec_SymItem
	endif

	Variable m0,m1,m2, b
	Variable i, err=0
	for (i=0;i<3;i+=1)
		err = err || ParseOneSymEquation(StringFromList(i,symItem,","),m0,m1,m2,b)
		mat[i][0] = m0
		mat[i][1] = m1
		mat[i][2] = m2
		vec[i] = b
	endfor
	return err
End
//
Static Function ParseOneSymEquation(expression, mx,my,mz,b)		// parse one expression of form "-x+y"  or "-x", or "-x+y, etc.
	String expression
	Variable &mx, &my, &mz, &b

	if (char2num(expression)>45)				// expression does not start with a '+' or '-', (43 or 45)
		expression = "+" + expression		//   so add a leading '+'
	endif
	expression = LowerStr(expression)		// only lower case
	expression = ReplaceString(" ",expression,"")	// no spaces

	mx=0  ;  my=0  ;  mz=0  ;  b=0
	// only certain fractions are allowed as constants
	String fractions="+1/2;-1/2;+1/3;-1/3;+2/3;-2/3;+1/4;-1/4;+3/4;-3/4;+1/6;-1/6;+5/6;-5/6;"
	Make/D/FREE fractionsV = {1/2,-1/2, 1/3,-1/3, 2/3,-2/3, 1/4,-1/4, 3/4,-3/4, 1/6,-1/6, 5/6,-5/6}
	Variable i,N=ItemsInList(fractions)
	for (i=0; i<N; i+=1)						// find the constant part, b
		if (strsearch(expression,StringFromList(i,fractions),0)>=0)
			b = fractionsV[i]
			break
		endif
	endfor

	if (strsearch(expression,"+x",0)>=0)	// find the x part mx
		mx = 1
	elseif (strsearch(expression,"-x",0)>=0)
		mx = -1
	endif

	if (strsearch(expression,"+y",0)>=0)	// find the y part my
		my = 1
	elseif (strsearch(expression,"-y",0)>=0)
		my = -1
	endif

	if (strsearch(expression,"+z",0)>=0)	// find the z part mz
		mz = 1
	elseif (strsearch(expression,"-z",0)>=0)
		mz = -1
	endif
	return 0
End
//Function testParseOne()
//	String expression = "-x+y-2/3"
//	Variable m0,m1,m2, b
//	LatticeSym#ParseOneSymEquation(expression,m0,m1,m2,b)
//	print expression,"    ",m0,"  ",m1,"  ",m2,"  ",b,"   ",expressionStr(m0,m1,m2,b)
//	expression = "x-y"
//	LatticeSym#ParseOneSymEquation(expression,m0,m1,m2,b)
//	print expression,"    ",m0,"  ",m1,"  ",m2,"  ",b,"   ",expressionStr(m0,m1,m2,b)
//	expression = "-z-1/4"
//	LatticeSym#ParseOneSymEquation(expression,m0,m1,m2,b)
//	print expression,"    ",m0,"  ",m1,"  ",m2,"  ",b,"   ",expressionStr(m0,m1,m2,b)
//
//	expression = "-1/4-z"
//	LatticeSym#ParseOneSymEquation(expression,m0,m1,m2,b)
//	print expression,"    ",m0,"  ",m1,"  ",m2,"  ",b,"   ",expressionStr(m0,m1,m2,b)
//
//	expression = "1/4-z"
//	LatticeSym#ParseOneSymEquation(expression,m0,m1,m2,b)
//	print expression,"    ",m0,"  ",m1,"  ",m2,"  ",b,"   ",expressionStr(m0,m1,m2,b)
//End
//
//
//
//Window Table_symLines() : Table
//	PauseUpdate; Silent 1		// building window...
//	Edit/W=(5,44,1379,833) numOps.xy,symLines.y
//	ModifyTable format(Point)=1,width(Point)=36,alignment(symLines.y)=0,width(symLines.y)=1248
//	ModifyTable width(numOps.x)=48,alignment(numOps.d)=1,width(numOps.d)=54
//EndMacro
//
Static Function/T expressionStr(mx,my,mz,b)	// turn a set of coefficients back into a string, opposite of ParseOneSymEquation()
	Variable mx,my,mz,b

	String str, out=""
	if (mx==1)
		out += "+x"
	elseif (mx==-1)
		out += "-x"
	elseif(mx)
		sprintf str, "%+.0fx", mx
		out += str
	endif

	if (my==1)
		out += "+y"
	elseif (my==-1)
		out += "-y"
	elseif(my)
		sprintf str, "%+.0fy", my
		out += str
	endif

	if (mz==1)
		out += "+z"
	elseif (mz==-1)
		out += "-z"
	elseif(mz)
		sprintf str, "%+.0fz", mz
		out += str
	endif
	out = RemoveLeadingString(out,"+",1)

	if (b)
		out += num2fraction(b,0.1,addPlus=1)	// turn b into a fraction string
	endif
	return out
End



Static Function/T setSymLineID(id)
	String id									// a space group id, e.g. "15" or "15:-b2"

	if (!isValidSpaceGroupID(id))		// perhaps only a number was passed
		Variable SG
		SG = str2num(id)
		SG = strsearch(id,":",0)>0 ? NaN : SG
		id = FindDefaultIDforSG(SG)		// find first space group starting with "id:"
	endif
	if (!isValidSpaceGroupID(id))
		return ""								// invalid
	endif

	return setSymLineIDnum(SpaceGroupID2num(id))
End
//
Static Function/T setSymLineIDnum(idNum)
	Variable idNum								// Space Group ID number [1,530]
	if (!isValidSpaceGroupIDnum(idNum))
		return ""								// invalid
	endif

	// These are the duplicate symmetry lines, I don't understand these duplicates????:
	//	324 & 322  ==  68:1ba-c & 68:1
	//	328 & 326  ==  68:1-cba & 68:1cab
	// 332 & 330  ==  68:1a-cb & 68:1bca
	String symLine322, symLine326, symLine330
	symLine322  = "x,y,z;-x,-y+1/2,-z+1/2;-x,-y,z;x,-y,-z;-x,y,-z;x,y+1/2,-z+1/2;-x,y+1/2,z+1/2;x,-y+1/2,z+1/2;x+1/2,y+1/2,z;"
	symLine322 += "-x+1/2,-y,-z+1/2;-x+1/2,-y+1/2,z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;x+1/2,y,-z+1/2;-x+1/2,y,z+1/2;x+1/2,-y,z+1/2"
	symLine326  = "x,y,z;-x+1/2,-y,-z+1/2;-x,-y,z;x,-y,-z;-x,y,-z;x+1/2,y,-z+1/2;-x+1/2,y,z+1/2;x+1/2,-y,z+1/2;x,y+1/2,z+1/2;"
	symLine326 += "-x+1/2,-y+1/2,-z;-x,-y+1/2,z+1/2;x,-y+1/2,-z+1/2;-x,y+1/2,-z+1/2;x+1/2,y+1/2,-z;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLine330  = "x,y,z;-x,-y+1/2,-z+1/2;-x,-y,z;x,-y,-z;-x,y,-z;x,y+1/2,-z+1/2;-x,y+1/2,z+1/2;x,-y+1/2,z+1/2;x+1/2,y,z+1/2;"
	symLine330 += "-x+1/2,-y+1/2,-z;-x+1/2,-y,z+1/2;x+1/2,-y,-z+1/2;-x+1/2,y,-z+1/2;x+1/2,y+1/2,-z;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z"

	Make/N=531/T/FREE symLines=""
	symLines[0] = ""
	// Triclinic SG[1,2]  SG_idNum [1-2]   (2 idNums)
	symLines[1]  = "x,y,z"
	symLines[2]  = "x,y,z;-x,-y,-z"
	// Monoclinic SG[3,15]  SG_idNum [3,107]   (105 idNums)
	symLines[3]  = "x,y,z;-x,y,-z"
	symLines[4]  = "x,y,z;-x,-y,z"
	symLines[5]  = "x,y,z;x,-y,-z"
	symLines[6]  = "x,y,z;-x,y+1/2,-z"
	symLines[7]  = "x,y,z;-x,-y,z+1/2"
	symLines[8]  = "x,y,z;x+1/2,-y,-z"
	symLines[9]  = "x,y,z;-x,y,-z;x+1/2,y+1/2,z;-x+1/2,y+1/2,-z"
	symLines[10]  = "x,y,z;-x,y,-z;x,y+1/2,z+1/2;-x,y+1/2,-z+1/2"
	symLines[11]  = "x,y,z;-x,y,-z;x+1/2,y+1/2,z+1/2;-x+1/2,y+1/2,-z+1/2"
	symLines[12]  = "x,y,z;-x,-y,z;x,y+1/2,z+1/2;-x,-y+1/2,z+1/2"
	symLines[13]  = "x,y,z;-x,-y,z;x+1/2,y,z+1/2;-x+1/2,-y,z+1/2"
	symLines[14]  = "x,y,z;-x,-y,z;x+1/2,y+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2"
	symLines[15]  = "x,y,z;x,-y,-z;x+1/2,y,z+1/2;x+1/2,-y,-z+1/2"
	symLines[16]  = "x,y,z;x,-y,-z;x+1/2,y+1/2,z;x+1/2,-y+1/2,-z"
	symLines[17]  = "x,y,z;x,-y,-z;x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,-z+1/2"
	symLines[18]  = "x,y,z;x,-y,z"
	symLines[19]  = "x,y,z;x,y,-z"
	symLines[20]  = "x,y,z;-x,y,z"
	symLines[21]  = "x,y,z;x,-y,z+1/2"
	symLines[22]  = "x,y,z;x+1/2,-y,z+1/2"
	symLines[23]  = "x,y,z;x+1/2,-y,z"
	symLines[24]  = "x,y,z;x+1/2,y,-z"
	symLines[25]  = "x,y,z;x+1/2,y+1/2,-z"
	symLines[26]  = "x,y,z;x,y+1/2,-z"
	symLines[27]  = "x,y,z;-x,y+1/2,z"
	symLines[28]  = "x,y,z;-x,y+1/2,z+1/2"
	symLines[29]  = "x,y,z;-x,y,z+1/2"
	symLines[30]  = "x,y,z;x,-y,z;x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[31]  = "x,y,z;x,-y,z;x,y+1/2,z+1/2;x,-y+1/2,z+1/2"
	symLines[32]  = "x,y,z;x,-y,z;x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[33]  = "x,y,z;x,y,-z;x,y+1/2,z+1/2;x,y+1/2,-z+1/2"
	symLines[34]  = "x,y,z;x,y,-z;x+1/2,y,z+1/2;x+1/2,y,-z+1/2"
	symLines[35]  = "x,y,z;x,y,-z;x+1/2,y+1/2,z+1/2;x+1/2,y+1/2,-z+1/2"
	symLines[36]  = "x,y,z;-x,y,z;x+1/2,y,z+1/2;-x+1/2,y,z+1/2"
	symLines[37]  = "x,y,z;-x,y,z;x+1/2,y+1/2,z;-x+1/2,y+1/2,z"
	symLines[38]  = "x,y,z;-x,y,z;x+1/2,y+1/2,z+1/2;-x+1/2,y+1/2,z+1/2"
	symLines[39]  = "x,y,z;x,-y,z+1/2;x+1/2,y+1/2,z;x+1/2,-y+1/2,z+1/2"
	symLines[40]  = "x,y,z;x+1/2,-y,z+1/2;x,y+1/2,z+1/2;x+1/2,-y+1/2,z"
	symLines[41]  = "x,y,z;x+1/2,-y,z;x+1/2,y+1/2,z+1/2;x,-y+1/2,z+1/2"
	symLines[42]  = "x,y,z;x+1/2,-y,z;x,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[43]  = "x,y,z;x,-y+1/2,z+1/2;x+1/2,y+1/2,z;x+1/2,-y,z+1/2"
	symLines[44]  = "x,y,z;x,-y,z+1/2;x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z"
	symLines[45]  = "x,y,z;x+1/2,y,-z;x,y+1/2,z+1/2;x+1/2,y+1/2,-z+1/2"
	symLines[46]  = "x,y,z;x,y+1/2,-z+1/2;x+1/2,y,z+1/2;x+1/2,y+1/2,-z"
	symLines[47]  = "x,y,z;x,y+1/2,-z;x+1/2,y+1/2,z+1/2;x+1/2,y,-z+1/2"
	symLines[48]  = "x,y,z;x,y+1/2,-z;x+1/2,y,z+1/2;x+1/2,y+1/2,-z+1/2"
	symLines[49]  = "x,y,z;x+1/2,y,-z+1/2;x,y+1/2,z+1/2;x+1/2,y+1/2,-z"
	symLines[50]  = "x,y,z;x+1/2,y,-z;x+1/2,y+1/2,z+1/2;x,y+1/2,-z+1/2"
	symLines[51]  = "x,y,z;-x,y+1/2,z;x+1/2,y,z+1/2;-x+1/2,y+1/2,z+1/2"
	symLines[52]  = "x,y,z;-x,y+1/2,z+1/2;x+1/2,y+1/2,z;-x+1/2,y,z+1/2"
	symLines[53]  = "x,y,z;-x,y,z+1/2;x+1/2,y+1/2,z+1/2;-x+1/2,y+1/2,z"
	symLines[54]  = "x,y,z;-x,y,z+1/2;x+1/2,y+1/2,z;-x+1/2,y+1/2,z+1/2"
	symLines[55]  = "x,y,z;-x,y+1/2,z+1/2;x+1/2,y,z+1/2;-x+1/2,y+1/2,z"
	symLines[56]  = "x,y,z;-x,y+1/2,z;x+1/2,y+1/2,z+1/2;-x+1/2,y,z+1/2"
	symLines[57]  = "x,y,z;-x,y,-z;-x,-y,-z;x,-y,z"
	symLines[58]  = "x,y,z;-x,-y,z;-x,-y,-z;x,y,-z"
	symLines[59]  = "x,y,z;x,-y,-z;-x,-y,-z;-x,y,z"
	symLines[60]  = "x,y,z;-x,y+1/2,-z;-x,-y,-z;x,-y+1/2,z"
	symLines[61]  = "x,y,z;-x,-y,z+1/2;-x,-y,-z;x,y,-z+1/2"
	symLines[62]  = "x,y,z;x+1/2,-y,-z;-x,-y,-z;-x+1/2,y,z"
	symLines[63]  = "x,y,z;-x,y,-z;-x,-y,-z;x,-y,z;x+1/2,y+1/2,z;-x+1/2,y+1/2,-z;-x+1/2,-y+1/2,-z;x+1/2,-y+1/2,z"
	symLines[64]  = "x,y,z;-x,y,-z;-x,-y,-z;x,-y,z;x,y+1/2,z+1/2;-x,y+1/2,-z+1/2;-x,-y+1/2,-z+1/2;x,-y+1/2,z+1/2"
	symLines[65]  = "x,y,z;-x,y,-z;-x,-y,-z;x,-y,z;x+1/2,y+1/2,z+1/2;-x+1/2,y+1/2,-z+1/2;-x+1/2,-y+1/2,-z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[66]  = "x,y,z;-x,-y,z;-x,-y,-z;x,y,-z;x,y+1/2,z+1/2;-x,-y+1/2,z+1/2;-x,-y+1/2,-z+1/2;x,y+1/2,-z+1/2"
	symLines[67]  = "x,y,z;-x,-y,z;-x,-y,-z;x,y,-z;x+1/2,y,z+1/2;-x+1/2,-y,z+1/2;-x+1/2,-y,-z+1/2;x+1/2,y,-z+1/2"
	symLines[68]  = "x,y,z;-x,-y,z;-x,-y,-z;x,y,-z;x+1/2,y+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;-x+1/2,-y+1/2,-z+1/2;x+1/2,y+1/2,-z+1/2"
	symLines[69]  = "x,y,z;x,-y,-z;-x,-y,-z;-x,y,z;x+1/2,y,z+1/2;x+1/2,-y,-z+1/2;-x+1/2,-y,-z+1/2;-x+1/2,y,z+1/2"
	symLines[70]  = "x,y,z;x,-y,-z;-x,-y,-z;-x,y,z;x+1/2,y+1/2,z;x+1/2,-y+1/2,-z;-x+1/2,-y+1/2,-z;-x+1/2,y+1/2,z"
	symLines[71]  = "x,y,z;x,-y,-z;-x,-y,-z;-x,y,z;x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,-z+1/2;-x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2"
	symLines[72]  = "x,y,z;-x,y,-z+1/2;-x,-y,-z;x,-y,z+1/2"
	symLines[73]  = "x,y,z;-x+1/2,y,-z+1/2;-x,-y,-z;x+1/2,-y,z+1/2"
	symLines[74]  = "x,y,z;-x+1/2,y,-z;-x,-y,-z;x+1/2,-y,z"
	symLines[75]  = "x,y,z;-x+1/2,-y,z;-x,-y,-z;x+1/2,y,-z"
	symLines[76]  = "x,y,z;-x+1/2,-y+1/2,z;-x,-y,-z;x+1/2,y+1/2,-z"
	symLines[77]  = "x,y,z;-x,-y+1/2,z;-x,-y,-z;x,y+1/2,-z"
	symLines[78]  = "x,y,z;x,-y+1/2,-z;-x,-y,-z;-x,y+1/2,z"
	symLines[79]  = "x,y,z;x,-y+1/2,-z+1/2;-x,-y,-z;-x,y+1/2,z+1/2"
	symLines[80]  = "x,y,z;x,-y,-z+1/2;-x,-y,-z;-x,y,z+1/2"
	symLines[81]  = "x,y,z;-x,y+1/2,-z+1/2;-x,-y,-z;x,-y+1/2,z+1/2"
	symLines[82]  = "x,y,z;-x+1/2,y+1/2,-z+1/2;-x,-y,-z;x+1/2,-y+1/2,z+1/2"
	symLines[83]  = "x,y,z;-x+1/2,y+1/2,-z;-x,-y,-z;x+1/2,-y+1/2,z"
	symLines[84]  = "x,y,z;-x+1/2,-y,z+1/2;-x,-y,-z;x+1/2,y,-z+1/2"
	symLines[85]  = "x,y,z;-x+1/2,-y+1/2,z+1/2;-x,-y,-z;x+1/2,y+1/2,-z+1/2"
	symLines[86]  = "x,y,z;-x,-y+1/2,z+1/2;-x,-y,-z;x,y+1/2,-z+1/2"
	symLines[87]  = "x,y,z;x+1/2,-y+1/2,-z;-x,-y,-z;-x+1/2,y+1/2,z"
	symLines[88]  = "x,y,z;x+1/2,-y+1/2,-z+1/2;-x,-y,-z;-x+1/2,y+1/2,z+1/2"
	symLines[89]  = "x,y,z;x+1/2,-y,-z+1/2;-x,-y,-z;-x+1/2,y,z+1/2"
	symLines[90]  = "x,y,z;-x,y,-z+1/2;-x,-y,-z;x,-y,z+1/2;x+1/2,y+1/2,z;-x+1/2,y+1/2,-z+1/2;-x+1/2,-y+1/2,-z;x+1/2,-y+1/2,z+1/2"
	symLines[91]  = "x,y,z;-x+1/2,y,-z+1/2;-x,-y,-z;x+1/2,-y,z+1/2;x,y+1/2,z+1/2;-x+1/2,y+1/2,-z;-x,-y+1/2,-z+1/2;x+1/2,-y+1/2,z"
	symLines[92]  = "x,y,z;-x+1/2,y,-z;-x,-y,-z;x+1/2,-y,z;x+1/2,y+1/2,z+1/2;-x,y+1/2,-z+1/2;-x+1/2,-y+1/2,-z+1/2;x,-y+1/2,z+1/2"
	symLines[93]  = "x,y,z;-x+1/2,y,-z;-x,-y,-z;x+1/2,-y,z;x,y+1/2,z+1/2;-x+1/2,y+1/2,-z+1/2;-x,-y+1/2,-z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[94]  = "x,y,z;-x,y+1/2,-z+1/2;-x,-y,-z;x,-y+1/2,z+1/2;x+1/2,y+1/2,z;-x+1/2,y,-z+1/2;-x+1/2,-y+1/2,-z;x+1/2,-y,z+1/2"
	symLines[95]  = "x,y,z;-x,y,-z+1/2;-x,-y,-z;x,-y,z+1/2;x+1/2,y+1/2,z+1/2;-x+1/2,y+1/2,-z;-x+1/2,-y+1/2,-z+1/2;x+1/2,-y+1/2,z"
	symLines[96]  = "x,y,z;-x+1/2,-y,z;-x,-y,-z;x+1/2,y,-z;x,y+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;-x,-y+1/2,-z+1/2;x+1/2,y+1/2,-z+1/2"
	symLines[97]  = "x,y,z;-x,-y+1/2,z+1/2;-x,-y,-z;x,y+1/2,-z+1/2;x+1/2,y,z+1/2;-x+1/2,-y+1/2,z;-x+1/2,-y,-z+1/2;x+1/2,y+1/2,-z"
	symLines[98]  = "x,y,z;-x,-y+1/2,z;-x,-y,-z;x,y+1/2,-z;x+1/2,y+1/2,z+1/2;-x+1/2,-y,z+1/2;-x+1/2,-y+1/2,-z+1/2;x+1/2,y,-z+1/2"
	symLines[99]  = "x,y,z;-x,-y+1/2,z;-x,-y,-z;x,y+1/2,-z;x+1/2,y,z+1/2;-x+1/2,-y+1/2,z+1/2;-x+1/2,-y,-z+1/2;x+1/2,y+1/2,-z+1/2"
	symLines[100]  = "x,y,z;-x+1/2,-y,z+1/2;-x,-y,-z;x+1/2,y,-z+1/2;x,y+1/2,z+1/2;-x+1/2,-y+1/2,z;-x,-y+1/2,-z+1/2;x+1/2,y+1/2,-z"
	symLines[101]  = "x,y,z;-x+1/2,-y,z;-x,-y,-z;x+1/2,y,-z;x+1/2,y+1/2,z+1/2;-x,-y+1/2,z+1/2;-x+1/2,-y+1/2,-z+1/2;x,y+1/2,-z+1/2"
	symLines[102]  = "x,y,z;x,-y+1/2,-z;-x,-y,-z;-x,y+1/2,z;x+1/2,y,z+1/2;x+1/2,-y+1/2,-z+1/2;-x+1/2,-y,-z+1/2;-x+1/2,y+1/2,z+1/2"
	symLines[103]  = "x,y,z;x,-y+1/2,-z+1/2;-x,-y,-z;-x,y+1/2,z+1/2;x+1/2,y+1/2,z;x+1/2,-y,-z+1/2;-x+1/2,-y+1/2,-z;-x+1/2,y,z+1/2"
	symLines[104]  = "x,y,z;x,-y,-z+1/2;-x,-y,-z;-x,y,z+1/2;x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,-z;-x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,z"
	symLines[105]  = "x,y,z;x,-y,-z+1/2;-x,-y,-z;-x,y,z+1/2;x+1/2,y+1/2,z;x+1/2,-y+1/2,-z+1/2;-x+1/2,-y+1/2,-z;-x+1/2,y+1/2,z+1/2"
	symLines[106]  = "x,y,z;x,-y+1/2,-z+1/2;-x,-y,-z;-x,y+1/2,z+1/2;x+1/2,y,z+1/2;x+1/2,-y+1/2,-z;-x+1/2,-y,-z+1/2;-x+1/2,y+1/2,z"
	symLines[107]  = "x,y,z;x,-y+1/2,-z;-x,-y,-z;-x,y+1/2,z;x+1/2,y+1/2,z+1/2;x+1/2,-y,-z+1/2;-x+1/2,-y+1/2,-z+1/2;-x+1/2,y,z+1/2"
	// Orthorhombic SG[16,74]  SG_idNum [108,348]   (241 idNums)
	symLines[108]  = "x,y,z;-x,-y,z;x,-y,-z;-x,y,-z"
	symLines[109]  = "x,y,z;-x,-y,z+1/2;x,-y,-z;-x,y,-z+1/2"
	symLines[110]  = "x,y,z;-x+1/2,-y,z;x+1/2,-y,-z;-x,y,-z"
	symLines[111]  = "x,y,z;-x,-y,z;x,-y+1/2,-z;-x,y+1/2,-z"
	symLines[112]  = "x,y,z;-x,-y,z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z"
	symLines[113]  = "x,y,z;-x,-y+1/2,z+1/2;x,-y,-z;-x,y+1/2,-z+1/2"
	symLines[114]  = "x,y,z;-x+1/2,-y,z+1/2;x+1/2,-y,-z+1/2;-x,y,-z"
	symLines[115]  = "x,y,z;-x+1/2,-y,z+1/2;x+1/2,-y+1/2,-z;-x,y+1/2,-z+1/2"
	symLines[116]  = "x,y,z;-x,-y,z+1/2;x,-y,-z;-x,y,-z+1/2;x+1/2,y+1/2,z;-x+1/2,-y+1/2,z+1/2;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z+1/2"
	symLines[117]  = "x,y,z;-x+1/2,-y,z;x+1/2,-y,-z;-x,y,-z;x,y+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;x+1/2,-y+1/2,-z+1/2;-x,y+1/2,-z+1/2"
	symLines[118]  = "x,y,z;-x,-y,z;x,-y+1/2,-z;-x,y+1/2,-z;x+1/2,y,z+1/2;-x+1/2,-y,z+1/2;x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2"
	symLines[119]  = "x,y,z;-x,-y,z;x,-y,-z;-x,y,-z;x+1/2,y+1/2,z;-x+1/2,-y+1/2,z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z"
	symLines[120]  = "x,y,z;-x,-y,z;x,-y,-z;-x,y,-z;x,y+1/2,z+1/2;-x,-y+1/2,z+1/2;x,-y+1/2,-z+1/2;-x,y+1/2,-z+1/2"
	symLines[121]  = "x,y,z;-x,-y,z;x,-y,-z;-x,y,-z;x+1/2,y,z+1/2;-x+1/2,-y,z+1/2;x+1/2,-y,-z+1/2;-x+1/2,y,-z+1/2"
	symLines[122]  = "x,y,z;-x,-y,z;x,-y,-z;-x,y,-z;x,y+1/2,z+1/2;-x,-y+1/2,z+1/2;x,-y+1/2,-z+1/2;-x,y+1/2,-z+1/2;x+1/2,y,z+1/2;"
	symLines[122] += "-x+1/2,-y,z+1/2;x+1/2,-y,-z+1/2;-x+1/2,y,-z+1/2;x+1/2,y+1/2,z;-x+1/2,-y+1/2,z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z"
	symLines[123]  = "x,y,z;-x,-y,z;x,-y,-z;-x,y,-z;x+1/2,y+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2"
	symLines[124]  = "x,y,z;-x,-y+1/2,z;x,-y,-z+1/2;-x+1/2,y,-z;x+1/2,y+1/2,z+1/2;-x+1/2,-y,z+1/2;x+1/2,-y+1/2,-z;-x,y+1/2,-z+1/2"
	symLines[125]  = "x,y,z;-x,-y,z;-x,y,z;x,-y,z"
	symLines[126]  = "x,y,z;x,-y,-z;x,y,-z;x,-y,z"
	symLines[127]  = "x,y,z;-x,y,-z;x,y,-z;-x,y,z"
	symLines[128]  = "x,y,z;-x,-y,z+1/2;-x,y,z;x,-y,z+1/2"
	symLines[129]  = "x,y,z;-x,-y,z+1/2;-x,y,z+1/2;x,-y,z"
	symLines[130]  = "x,y,z;x+1/2,-y,-z;x+1/2,y,-z;x,-y,z"
	symLines[131]  = "x,y,z;x+1/2,-y,-z;x,y,-z;x+1/2,-y,z"
	symLines[132]  = "x,y,z;-x,y+1/2,-z;x,y,-z;-x,y+1/2,z"
	symLines[133]  = "x,y,z;-x,y+1/2,-z;x,y+1/2,-z;-x,y,z"
	symLines[134]  = "x,y,z;-x,-y,z;-x,y,z+1/2;x,-y,z+1/2"
	symLines[135]  = "x,y,z;x,-y,-z;x+1/2,y,-z;x+1/2,-y,z"
	symLines[136]  = "x,y,z;-x,y,-z;x,y+1/2,-z;-x,y+1/2,z"
	symLines[137]  = "x,y,z;-x,-y,z;-x+1/2,y,z;x+1/2,-y,z"
	symLines[138]  = "x,y,z;-x,-y,z;-x,y+1/2,z;x,-y+1/2,z"
	symLines[139]  = "x,y,z;x,-y,-z;x,y+1/2,-z;x,-y+1/2,z"
	symLines[140]  = "x,y,z;x,-y,-z;x,y,-z+1/2;x,-y,z+1/2"
	symLines[141]  = "x,y,z;-x,y,-z;x,y,-z+1/2;-x,y,z+1/2"
	symLines[142]  = "x,y,z;-x,y,-z;x+1/2,y,-z;-x+1/2,y,z"
	symLines[143]  = "x,y,z;-x,-y,z+1/2;-x+1/2,y,z+1/2;x+1/2,-y,z"
	symLines[144]  = "x,y,z;-x,-y,z+1/2;-x,y+1/2,z;x,-y+1/2,z+1/2"
	symLines[145]  = "x,y,z;x+1/2,-y,-z;x,y+1/2,-z;x+1/2,-y+1/2,z"
	symLines[146]  = "x,y,z;x+1/2,-y,-z;x+1/2,y,-z+1/2;x,-y,z+1/2"
	symLines[147]  = "x,y,z;-x,y+1/2,-z;x,y+1/2,-z+1/2;-x,y,z+1/2"
	symLines[148]  = "x,y,z;-x,y+1/2,-z;x+1/2,y,-z;-x+1/2,y+1/2,z"
	symLines[149]  = "x,y,z;-x,-y,z;-x,y+1/2,z+1/2;x,-y+1/2,z+1/2"
	symLines[150]  = "x,y,z;-x,-y,z;-x+1/2,y,z+1/2;x+1/2,-y,z+1/2"
	symLines[151]  = "x,y,z;x,-y,-z;x+1/2,y,-z+1/2;x+1/2,-y,z+1/2"
	symLines[152]  = "x,y,z;x,-y,-z;x+1/2,y+1/2,-z;x+1/2,-y+1/2,z"
	symLines[153]  = "x,y,z;-x,y,-z;x+1/2,y+1/2,-z;-x+1/2,y+1/2,z"
	symLines[154]  = "x,y,z;-x,y,-z;x,y+1/2,-z+1/2;-x,y+1/2,z+1/2"
	symLines[155]  = "x,y,z;-x+1/2,-y,z+1/2;-x,y,z;x+1/2,-y,z+1/2"
	symLines[156]  = "x,y,z;-x,-y+1/2,z+1/2;-x,y+1/2,z+1/2;x,-y,z"
	symLines[157]  = "x,y,z;x+1/2,-y+1/2,-z;x+1/2,y+1/2,-z;x,-y,z"
	symLines[158]  = "x,y,z;x+1/2,-y,-z+1/2;x,y,-z;x+1/2,-y,z+1/2"
	symLines[159]  = "x,y,z;-x,y+1/2,-z+1/2;x,y,-z;-x,y+1/2,z+1/2"
	symLines[160]  = "x,y,z;-x+1/2,y+1/2,-z;x+1/2,y+1/2,-z;-x,y,z"
	symLines[161]  = "x,y,z;-x,-y,z;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[162]  = "x,y,z;x,-y,-z;x,y+1/2,-z+1/2;x,-y+1/2,z+1/2"
	symLines[163]  = "x,y,z;-x,y,-z;x+1/2,y,-z+1/2;-x+1/2,y,z+1/2"
	symLines[164]  = "x,y,z;-x,-y,z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z"
	symLines[165]  = "x,y,z;-x,-y,z+1/2;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z+1/2"
	symLines[166]  = "x,y,z;x+1/2,-y,-z;x,y+1/2,-z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[167]  = "x,y,z;x+1/2,-y,-z;x+1/2,y+1/2,-z+1/2;x,-y+1/2,z+1/2"
	symLines[168]  = "x,y,z;-x,y+1/2,-z;x+1/2,y+1/2,-z+1/2;-x+1/2,y,z+1/2"
	symLines[169]  = "x,y,z;-x,y+1/2,-z;x+1/2,y,-z+1/2;-x+1/2,y+1/2,z+1/2"
	symLines[170]  = "x,y,z;-x,-y,z;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[171]  = "x,y,z;x,-y,-z;x+1/2,y+1/2,-z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[172]  = "x,y,z;-x,y,-z;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2"
	symLines[173]  = "x,y,z;-x,-y,z;-x,y,z;x,-y,z;x+1/2,y+1/2,z;-x+1/2,-y+1/2,z;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[174]  = "x,y,z;x,-y,-z;x,y,-z;x,-y,z;x,y+1/2,z+1/2;x,-y+1/2,-z+1/2;x,y+1/2,-z+1/2;x,-y+1/2,z+1/2"
	symLines[175]  = "x,y,z;-x,y,-z;x,y,-z;-x,y,z;x+1/2,y,z+1/2;-x+1/2,y,-z+1/2;x+1/2,y,-z+1/2;-x+1/2,y,z+1/2"
	symLines[176]  = "x,y,z;-x,-y,z+1/2;-x,y,z;x,-y,z+1/2;x+1/2,y+1/2,z;-x+1/2,-y+1/2,z+1/2;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z+1/2"
	symLines[177]  = "x,y,z;-x,-y,z+1/2;-x,y,z+1/2;x,-y,z;x+1/2,y+1/2,z;-x+1/2,-y+1/2,z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z"
	symLines[178]  = "x,y,z;x+1/2,-y,-z;x+1/2,y,-z;x,-y,z;x,y+1/2,z+1/2;x+1/2,-y+1/2,-z+1/2;x+1/2,y+1/2,-z+1/2;x,-y+1/2,z+1/2"
	symLines[179]  = "x,y,z;x+1/2,-y,-z;x,y,-z;x+1/2,-y,z;x,y+1/2,z+1/2;x+1/2,-y+1/2,-z+1/2;x,y+1/2,-z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[180]  = "x,y,z;-x,y+1/2,-z;x,y,-z;-x,y+1/2,z;x+1/2,y,z+1/2;-x+1/2,y+1/2,-z+1/2;x+1/2,y,-z+1/2;-x+1/2,y+1/2,z+1/2"
	symLines[181]  = "x,y,z;-x,y+1/2,-z;x,y+1/2,-z;-x,y,z;x+1/2,y,z+1/2;-x+1/2,y+1/2,-z+1/2;x+1/2,y+1/2,-z+1/2;-x+1/2,y,z+1/2"
	symLines[182]  = "x,y,z;-x,-y,z;-x,y,z+1/2;x,-y,z+1/2;x+1/2,y+1/2,z;-x+1/2,-y+1/2,z;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[183]  = "x,y,z;x,-y,-z;x+1/2,y,-z;x+1/2,-y,z;x,y+1/2,z+1/2;x,-y+1/2,-z+1/2;x+1/2,y+1/2,-z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[184]  = "x,y,z;-x,y,-z;x,y+1/2,-z;-x,y+1/2,z;x+1/2,y,z+1/2;-x+1/2,y,-z+1/2;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2"
	symLines[185]  = "x,y,z;-x,-y,z;-x,y,z;x,-y,z;x,y+1/2,z+1/2;-x,-y+1/2,z+1/2;-x,y+1/2,z+1/2;x,-y+1/2,z+1/2"
	symLines[186]  = "x,y,z;-x,-y,z;-x,y,z;x,-y,z;x+1/2,y,z+1/2;-x+1/2,-y,z+1/2;-x+1/2,y,z+1/2;x+1/2,-y,z+1/2"
	symLines[187]  = "x,y,z;x,-y,-z;x,y,-z;x,-y,z;x+1/2,y,z+1/2;x+1/2,-y,-z+1/2;x+1/2,y,-z+1/2;x+1/2,-y,z+1/2"
	symLines[188]  = "x,y,z;x,-y,-z;x,y,-z;x,-y,z;x+1/2,y+1/2,z;x+1/2,-y+1/2,-z;x+1/2,y+1/2,-z;x+1/2,-y+1/2,z"
	symLines[189]  = "x,y,z;-x,y,-z;x,y,-z;-x,y,z;x+1/2,y+1/2,z;-x+1/2,y+1/2,-z;x+1/2,y+1/2,-z;-x+1/2,y+1/2,z"
	symLines[190]  = "x,y,z;-x,y,-z;x,y,-z;-x,y,z;x,y+1/2,z+1/2;-x,y+1/2,-z+1/2;x,y+1/2,-z+1/2;-x,y+1/2,z+1/2"
	symLines[191]  = "x,y,z;-x,-y,z;-x,y,z+1/2;x,-y,z+1/2;x,y+1/2,z+1/2;-x,-y+1/2,z+1/2;-x,y+1/2,z;x,-y+1/2,z"
	symLines[192]  = "x,y,z;-x,-y,z;-x,y,z+1/2;x,-y,z+1/2;x+1/2,y,z+1/2;-x+1/2,-y,z+1/2;-x+1/2,y,z;x+1/2,-y,z"
	symLines[193]  = "x,y,z;x,-y,-z;x,y,-z+1/2;x,-y,z+1/2;x+1/2,y,z+1/2;x+1/2,-y,-z+1/2;x+1/2,y,-z;x+1/2,-y,z"
	symLines[194]  = "x,y,z;x,-y,-z;x,y+1/2,-z;x,-y+1/2,z;x+1/2,y+1/2,z;x+1/2,-y+1/2,-z;x+1/2,y,-z;x+1/2,-y,z"
	symLines[195]  = "x,y,z;-x,y,-z;x,y+1/2,-z;-x,y+1/2,z;x+1/2,y+1/2,z;-x+1/2,y+1/2,-z;x+1/2,y,-z;-x+1/2,y,z"
	symLines[196]  = "x,y,z;-x,y,-z;x,y,-z+1/2;-x,y,z+1/2;x,y+1/2,z+1/2;-x,y+1/2,-z+1/2;x,y+1/2,-z;-x,y+1/2,z"
	symLines[197]  = "x,y,z;-x,-y,z;-x+1/2,y,z;x+1/2,-y,z;x,y+1/2,z+1/2;-x,-y+1/2,z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[198]  = "x,y,z;-x,-y,z;-x,y+1/2,z;x,-y+1/2,z;x+1/2,y,z+1/2;-x+1/2,-y,z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[199]  = "x,y,z;x,-y,-z;x,y+1/2,-z;x,-y+1/2,z;x+1/2,y,z+1/2;x+1/2,-y,-z+1/2;x+1/2,y+1/2,-z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[200]  = "x,y,z;x,-y,-z;x,y,-z+1/2;x,-y,z+1/2;x+1/2,y+1/2,z;x+1/2,-y+1/2,-z;x+1/2,y+1/2,-z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[201]  = "x,y,z;-x,y,-z;x,y,-z+1/2;-x,y,z+1/2;x+1/2,y+1/2,z;-x+1/2,y+1/2,-z;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2"
	symLines[202]  = "x,y,z;-x,y,-z;x+1/2,y,-z;-x+1/2,y,z;x,y+1/2,z+1/2;-x,y+1/2,-z+1/2;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2"
	symLines[203]  = "x,y,z;-x,-y,z;-x+1/2,y,z+1/2;x+1/2,-y,z+1/2;x,y+1/2,z+1/2;-x,-y+1/2,z+1/2;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[204]  = "x,y,z;-x,-y,z;-x,y+1/2,z+1/2;x,-y+1/2,z+1/2;x+1/2,y,z+1/2;-x+1/2,-y,z+1/2;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[205]  = "x,y,z;x,-y,-z;x,y+1/2,-z+1/2;x,-y+1/2,z+1/2;x+1/2,y,z+1/2;x+1/2,-y,-z+1/2;x+1/2,y+1/2,-z;x+1/2,-y+1/2,z"
	symLines[206]  = "x,y,z;x,-y,-z;x,y+1/2,-z+1/2;x,-y+1/2,z+1/2;x+1/2,y+1/2,z;x+1/2,-y+1/2,-z;x+1/2,y,-z+1/2;x+1/2,-y,z+1/2"
	symLines[207]  = "x,y,z;-x,y,-z;x,y+1/2,-z+1/2;-x,y+1/2,z+1/2;x+1/2,y+1/2,z;-x+1/2,y+1/2,-z;x+1/2,y,-z+1/2;-x+1/2,y,z+1/2"
	symLines[208]  = "x,y,z;-x,y,-z;x+1/2,y,-z+1/2;-x+1/2,y,z+1/2;x,y+1/2,z+1/2;-x,y+1/2,-z+1/2;x+1/2,y+1/2,-z;-x+1/2,y+1/2,z"
	symLines[209]  = "x,y,z;-x,-y,z;-x,y,z;x,-y,z;x,y+1/2,z+1/2;-x,-y+1/2,z+1/2;-x,y+1/2,z+1/2;x,-y+1/2,z+1/2;x+1/2,y,z+1/2;-x+1/2,-y,z+1/2;"
	symLines[209] += "-x+1/2,y,z+1/2;x+1/2,-y,z+1/2;x+1/2,y+1/2,z;-x+1/2,-y+1/2,z;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[210]  = "x,y,z;x,-y,-z;x,y,-z;x,-y,z;x,y+1/2,z+1/2;x,-y+1/2,-z+1/2;x,y+1/2,-z+1/2;x,-y+1/2,z+1/2;x+1/2,y,z+1/2;x+1/2,-y,-z+1/2;"
	symLines[210] += "x+1/2,y,-z+1/2;x+1/2,-y,z+1/2;x+1/2,y+1/2,z;x+1/2,-y+1/2,-z;x+1/2,y+1/2,-z;x+1/2,-y+1/2,z"
	symLines[211]  = "x,y,z;-x,y,-z;x,y,-z;-x,y,z;x,y+1/2,z+1/2;-x,y+1/2,-z+1/2;x,y+1/2,-z+1/2;-x,y+1/2,z+1/2;x+1/2,y,z+1/2;-x+1/2,y,-z+1/2;"
	symLines[211] += "x+1/2,y,-z+1/2;-x+1/2,y,z+1/2;x+1/2,y+1/2,z;-x+1/2,y+1/2,-z;x+1/2,y+1/2,-z;-x+1/2,y+1/2,z"
	symLines[212]  = "x,y,z;-x,-y,z;-x+1/4,y+1/4,z+1/4;x+1/4,-y+1/4,z+1/4;x,y+1/2,z+1/2;-x,-y+1/2,z+1/2;-x+1/4,y-1/4,z-1/4;x+1/4,-y-1/4,z-1/4;"
	symLines[212] += "x+1/2,y,z+1/2;-x+1/2,-y,z+1/2;-x-1/4,y+1/4,z-1/4;x-1/4,-y+1/4,z-1/4;x+1/2,y+1/2,z;-x+1/2,-y+1/2,z;-x-1/4,y-1/4,z+1/4;x-1/4,-y-1/4,z+1/4"
	symLines[213]  = "x,y,z;x,-y,-z;x+1/4,y+1/4,-z+1/4;x+1/4,-y+1/4,z+1/4;x,y+1/2,z+1/2;x,-y+1/2,-z+1/2;x+1/4,y-1/4,-z-1/4;x+1/4,-y-1/4,z-1/4;"
	symLines[213] += "x+1/2,y,z+1/2;x+1/2,-y,-z+1/2;x-1/4,y+1/4,-z-1/4;x-1/4,-y+1/4,z-1/4;x+1/2,y+1/2,z;x+1/2,-y+1/2,-z;x-1/4,y-1/4,-z+1/4;x-1/4,-y-1/4,z+1/4"
	symLines[214]  = "x,y,z;-x,y,-z;x+1/4,y+1/4,-z+1/4;-x+1/4,y+1/4,z+1/4;x,y+1/2,z+1/2;-x,y+1/2,-z+1/2;x+1/4,y-1/4,-z-1/4;-x+1/4,y-1/4,z-1/4;"
	symLines[214] += "x+1/2,y,z+1/2;-x+1/2,y,-z+1/2;x-1/4,y+1/4,-z-1/4;-x-1/4,y+1/4,z-1/4;x+1/2,y+1/2,z;-x+1/2,y+1/2,-z;x-1/4,y-1/4,-z+1/4;-x-1/4,y-1/4,z+1/4"
	symLines[215]  = "x,y,z;-x,-y,z;-x,y,z;x,-y,z;x+1/2,y+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[216]  = "x,y,z;x,-y,-z;x,y,-z;x,-y,z;x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,-z+1/2;x+1/2,y+1/2,-z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[217]  = "x,y,z;-x,y,-z;x,y,-z;-x,y,z;x+1/2,y+1/2,z+1/2;-x+1/2,y+1/2,-z+1/2;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2"
	symLines[218]  = "x,y,z;-x,-y,z;-x,y,z+1/2;x,-y,z+1/2;x+1/2,y+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[219]  = "x,y,z;x,-y,-z;x+1/2,y,-z;x+1/2,-y,z;x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,-z+1/2;x,y+1/2,-z+1/2;x,-y+1/2,z+1/2"
	symLines[220]  = "x,y,z;-x,y,-z;x,y+1/2,-z;-x,y+1/2,z;x+1/2,y+1/2,z+1/2;-x+1/2,y+1/2,-z+1/2;x+1/2,y,-z+1/2;-x+1/2,y,z+1/2"
	symLines[221]  = "x,y,z;-x,-y,z;-x+1/2,y,z;x+1/2,-y,z;x+1/2,y+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;-x,y+1/2,z+1/2;x,-y+1/2,z+1/2"
	symLines[222]  = "x,y,z;-x,-y,z;-x,y+1/2,z;x,-y+1/2,z;x+1/2,y+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;-x+1/2,y,z+1/2;x+1/2,-y,z+1/2"
	symLines[223]  = "x,y,z;x,-y,-z;x,y+1/2,-z;x,-y+1/2,z;x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,-z+1/2;x+1/2,y,-z+1/2;x+1/2,-y,z+1/2"
	symLines[224]  = "x,y,z;x,-y,-z;x,y,-z+1/2;x,-y,z+1/2;x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,-z+1/2;x+1/2,y+1/2,-z;x+1/2,-y+1/2,z"
	symLines[225]  = "x,y,z;-x,y,-z;x,y,-z+1/2;-x,y,z+1/2;x+1/2,y+1/2,z+1/2;-x+1/2,y+1/2,-z+1/2;x+1/2,y+1/2,-z;-x+1/2,y+1/2,z"
	symLines[226]  = "x,y,z;-x,y,-z;x+1/2,y,-z;-x+1/2,y,z;x+1/2,y+1/2,z+1/2;-x+1/2,y+1/2,-z+1/2;x,y+1/2,-z+1/2;-x,y+1/2,z+1/2"
	symLines[227]  = "x,y,z;-x,-y,z;x,-y,-z;-x,y,-z;-x,-y,-z;x,y,-z;-x,y,z;x,-y,z"
	symLines[228]  = "x,y,z;-x+1/2,-y+1/2,-z+1/2;-x,-y,z;x,-y,-z;-x,y,-z;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[229]  = "x,y,z;-x+1/2,-y+1/2,z;x,-y+1/2,-z+1/2;-x+1/2,y,-z+1/2;-x,-y,-z;x+1/2,y+1/2,-z;-x,y+1/2,z+1/2;x+1/2,-y,z+1/2"
	symLines[230]  = "x,y,z;-x,-y,z;x,-y,-z+1/2;-x,y,-z+1/2;-x,-y,-z;x,y,-z;-x,y,z+1/2;x,-y,z+1/2"
	symLines[231]  = "x,y,z;-x+1/2,-y,z;x,-y,-z;-x+1/2,y,-z;-x,-y,-z;x+1/2,y,-z;-x,y,z;x+1/2,-y,z"
	symLines[232]  = "x,y,z;-x,-y+1/2,z;x,-y+1/2,-z;-x,y,-z;-x,-y,-z;x,y+1/2,-z;-x,y+1/2,z;x,-y,z"
	symLines[233]  = "x,y,z;-x+1/2,-y+1/2,-z;-x,-y,z;x,-y,-z;-x,y,-z;x+1/2,y+1/2,-z;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[234]  = "x,y,z;-x+1/2,-y+1/2,z;x,-y+1/2,-z;-x+1/2,y,-z;-x,-y,-z;x+1/2,y+1/2,-z;-x,y+1/2,z;x+1/2,-y,z"
	symLines[235]  = "x,y,z;-x,-y+1/2,-z+1/2;-x,-y,z;x,-y,-z;-x,y,-z;x,y+1/2,-z+1/2;-x,y+1/2,z+1/2;x,-y+1/2,z+1/2"
	symLines[236]  = "x,y,z;-x,-y+1/2,z;x,-y+1/2,-z+1/2;-x,y,-z+1/2;-x,-y,-z;x,y+1/2,-z;-x,y+1/2,z+1/2;x,-y,z+1/2"
	symLines[237]  = "x,y,z;-x+1/2,-y,-z+1/2;-x,-y,z;x,-y,-z;-x,y,-z;x+1/2,y,-z+1/2;-x+1/2,y,z+1/2;x+1/2,-y,z+1/2"
	symLines[238]  = "x,y,z;-x+1/2,-y,z;x,-y,-z+1/2;-x+1/2,y,-z+1/2;-x,-y,-z;x+1/2,y,-z;-x,y,z+1/2;x+1/2,-y,z+1/2"
	symLines[239]  = "x,y,z;-x+1/2,-y,z;x+1/2,-y,-z;-x,y,-z;-x,-y,-z;x+1/2,y,-z;-x+1/2,y,z;x,-y,z"
	symLines[240]  = "x,y,z;-x,-y+1/2,z;x,-y,-z;-x,y+1/2,-z;-x,-y,-z;x,y+1/2,-z;-x,y,z;x,-y+1/2,z"
	symLines[241]  = "x,y,z;-x,-y,z;x,-y+1/2,-z;-x,y+1/2,-z;-x,-y,-z;x,y,-z;-x,y+1/2,z;x,-y+1/2,z"
	symLines[242]  = "x,y,z;-x,-y,z+1/2;x,-y,-z+1/2;-x,y,-z;-x,-y,-z;x,y,-z+1/2;-x,y,z+1/2;x,-y,z"
	symLines[243]  = "x,y,z;-x,-y,z+1/2;x,-y,-z;-x,y,-z+1/2;-x,-y,-z;x,y,-z+1/2;-x,y,z;x,-y,z+1/2"
	symLines[244]  = "x,y,z;-x,-y,z;x+1/2,-y,-z;-x+1/2,y,-z;-x,-y,-z;x,y,-z;-x+1/2,y,z;x+1/2,-y,z"
	symLines[245]  = "x,y,z;-x+1/2,-y,z;x,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2;-x,-y,-z;x+1/2,y,-z;-x,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[246]  = "x,y,z;-x,-y+1/2,z;x+1/2,-y+1/2,-z+1/2;-x+1/2,y,-z+1/2;-x,-y,-z;x,y+1/2,-z;-x+1/2,y+1/2,z+1/2;x+1/2,-y,z+1/2"
	symLines[247]  = "x,y,z;-x+1/2,-y+1/2,z+1/2;x,-y+1/2,-z;-x+1/2,y,-z+1/2;-x,-y,-z;x+1/2,y+1/2,-z+1/2;-x,y+1/2,z;x+1/2,-y,z+1/2"
	symLines[248]  = "x,y,z;-x+1/2,-y+1/2,z;x,-y,-z+1/2;-x+1/2,y+1/2,-z+1/2;-x,-y,-z;x+1/2,y+1/2,-z;-x,y,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[249]  = "x,y,z;-x+1/2,-y+1/2,z;x+1/2,-y+1/2,-z+1/2;-x,y,-z+1/2;-x,-y,-z;x+1/2,y+1/2,-z;-x+1/2,y+1/2,z+1/2;x,-y,z+1/2"
	symLines[250]  = "x,y,z;-x+1/2,-y+1/2,z+1/2;x,-y+1/2,-z+1/2;-x+1/2,y,-z;-x,-y,-z;x+1/2,y+1/2,-z+1/2;-x,y+1/2,z+1/2;x+1/2,-y,z"
	symLines[251]  = "x,y,z;-x+1/2,-y,z+1/2;x,-y,-z;-x+1/2,y,-z+1/2;-x,-y,-z;x+1/2,y,-z+1/2;-x,y,z;x+1/2,-y,z+1/2"
	symLines[252]  = "x,y,z;-x,-y+1/2,z+1/2;x,-y+1/2,-z+1/2;-x,y,-z;-x,-y,-z;x,y+1/2,-z+1/2;-x,y+1/2,z+1/2;x,-y,z"
	symLines[253]  = "x,y,z;-x+1/2,-y+1/2,z;x+1/2,-y+1/2,-z;-x,y,-z;-x,-y,-z;x+1/2,y+1/2,-z;-x+1/2,y+1/2,z;x,-y,z"
	symLines[254]  = "x,y,z;-x,-y,z;x+1/2,-y,-z+1/2;-x+1/2,y,-z+1/2;-x,-y,-z;x,y,-z;-x+1/2,y,z+1/2;x+1/2,-y,z+1/2"
	symLines[255]  = "x,y,z;-x,-y,z;x,-y+1/2,-z+1/2;-x,y+1/2,-z+1/2;-x,-y,-z;x,y,-z;-x,y+1/2,z+1/2;x,-y+1/2,z+1/2"
	symLines[256]  = "x,y,z;-x+1/2,-y+1/2,z;x,-y,-z;-x+1/2,y+1/2,-z;-x,-y,-z;x+1/2,y+1/2,-z;-x,y,z;x+1/2,-y+1/2,z"
	symLines[257]  = "x,y,z;-x+1/2,-y,z;x+1/2,-y,-z+1/2;-x,y,-z+1/2;-x,-y,-z;x+1/2,y,-z;-x+1/2,y,z+1/2;x,-y,z+1/2"
	symLines[258]  = "x,y,z;-x,-y+1/2,z;x,-y,-z+1/2;-x,y+1/2,-z+1/2;-x,-y,-z;x,y+1/2,-z;-x,y,z+1/2;x,-y+1/2,z+1/2"
	symLines[259]  = "x,y,z;-x+1/2,-y,z;x,-y+1/2,-z;-x+1/2,y+1/2,-z;-x,-y,-z;x+1/2,y,-z;-x,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[260]  = "x,y,z;-x+1/2,-y,z+1/2;x,-y,-z+1/2;-x+1/2,y,-z;-x,-y,-z;x+1/2,y,-z+1/2;-x,y,z+1/2;x+1/2,-y,z"
	symLines[261]  = "x,y,z;-x,-y+1/2,z+1/2;x,-y+1/2,-z;-x,y,-z+1/2;-x,-y,-z;x,y+1/2,-z+1/2;-x,y+1/2,z;x,-y,z+1/2"
	symLines[262]  = "x,y,z;-x,-y+1/2,z;x+1/2,-y+1/2,-z;-x+1/2,y,-z;-x,-y,-z;x,y+1/2,-z;-x+1/2,y+1/2,z;x+1/2,-y,z"
	symLines[263]  = "x,y,z;-x,-y,z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;-x,-y,-z;x,y,-z;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[264]  = "x,y,z;-x,-y+1/2,z+1/2;x,-y,-z;-x,y+1/2,-z+1/2;-x,-y,-z;x,y+1/2,-z+1/2;-x,y,z;x,-y+1/2,z+1/2"
	symLines[265]  = "x,y,z;-x+1/2,-y,z+1/2;x+1/2,-y,-z+1/2;-x,y,-z;-x,-y,-z;x+1/2,y,-z+1/2;-x+1/2,y,z+1/2;x,-y,z"
	symLines[266]  = "x,y,z;-x+1/2,-y+1/2,z;x+1/2,-y,-z+1/2;-x,y+1/2,-z+1/2;-x,-y,-z;x+1/2,y+1/2,-z;-x+1/2,y,z+1/2;x,-y+1/2,z+1/2"
	symLines[267]  = "x,y,z;-x+1/2,-y,z+1/2;x,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z;-x,-y,-z;x+1/2,y,-z+1/2;-x,y+1/2,z+1/2;x+1/2,-y+1/2,z"
	symLines[268]  = "x,y,z;-x,-y+1/2,z+1/2;x+1/2,-y+1/2,-z;-x+1/2,y,-z+1/2;-x,-y,-z;x,y+1/2,-z+1/2;-x+1/2,y+1/2,z;x+1/2,-y,z+1/2"
	symLines[269]  = "x,y,z;-x,-y,z+1/2;x,-y+1/2,-z;-x,y+1/2,-z+1/2;-x,-y,-z;x,y,-z+1/2;-x,y+1/2,z;x,-y+1/2,z+1/2"
	symLines[270]  = "x,y,z;-x,-y,z+1/2;x+1/2,-y,-z+1/2;-x+1/2,y,-z;-x,-y,-z;x,y,-z+1/2;-x+1/2,y,z+1/2;x+1/2,-y,z"
	symLines[271]  = "x,y,z;-x+1/2,-y,z+1/2;x+1/2,-y,-z;-x,y,-z+1/2;-x,-y,-z;x+1/2,y,-z+1/2;-x+1/2,y,z;x,-y,z+1/2"
	symLines[272]  = "x,y,z;-x,-y+1/2,z;x+1/2,-y,-z;-x+1/2,y+1/2,-z;-x,-y,-z;x,y+1/2,-z;-x+1/2,y,z;x+1/2,-y+1/2,z"
	symLines[273]  = "x,y,z;-x+1/2,-y,z;x+1/2,-y+1/2,-z;-x,y+1/2,-z;-x,-y,-z;x+1/2,y,-z;-x+1/2,y+1/2,z;x,-y+1/2,z"
	symLines[274]  = "x,y,z;-x,-y+1/2,z+1/2;x,-y,-z+1/2;-x,y+1/2,-z;-x,-y,-z;x,y+1/2,-z+1/2;-x,y,z+1/2;x,-y+1/2,z"
	symLines[275]  = "x,y,z;-x,-y,z;x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2;-x,-y,-z;x,y,-z;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[276]  = "x,y,z;-x+1/2,-y+1/2,z+1/2;x,-y,-z;-x+1/2,y+1/2,-z+1/2;-x,-y,-z;x+1/2,y+1/2,-z+1/2;-x,y,z;x+1/2,-y+1/2,z+1/2"
	symLines[277]  = "x,y,z;-x+1/2,-y+1/2,z+1/2;x+1/2,-y+1/2,-z+1/2;-x,y,-z;-x,-y,-z;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2;x,-y,z"
	symLines[278]  = "x,y,z;-x+1/2,-y+1/2,-z;-x,-y,z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;x+1/2,y+1/2,-z;-x,y,z;x,-y,z"
	symLines[279]  = "x,y,z;-x+1/2,-y+1/2,z;x+1/2,-y,-z;-x,y+1/2,-z;-x,-y,-z;x+1/2,y+1/2,-z;-x+1/2,y,z;x,-y+1/2,z"
	symLines[280]  = "x,y,z;-x,-y+1/2,-z+1/2;-x,-y+1/2,z+1/2;x,-y,-z;-x,y+1/2,-z+1/2;x,y,-z;-x,y+1/2,z+1/2;x,-y,z"
	symLines[281]  = "x,y,z;-x,-y,z+1/2;x,-y+1/2,-z+1/2;-x,y+1/2,-z;-x,-y,-z;x,y,-z+1/2;-x,y+1/2,z+1/2;x,-y+1/2,z"
	symLines[282]  = "x,y,z;-x+1/2,-y,-z+1/2;-x+1/2,-y,z+1/2;x+1/2,-y,-z+1/2;-x,y,-z;x,y,-z;-x,y,z;x+1/2,-y,z+1/2"
	symLines[283]  = "x,y,z;-x,-y,z+1/2;x+1/2,-y,-z;-x+1/2,y,-z+1/2;-x,-y,-z;x,y,-z+1/2;-x+1/2,y,z;x+1/2,-y,z+1/2"
	symLines[284]  = "x,y,z;-x+1/2,-y+1/2,z+1/2;x+1/2,-y+1/2,-z;-x,y,-z+1/2;-x,-y,-z;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z;x,-y,z+1/2"
	symLines[285]  = "x,y,z;-x+1/2,-y+1/2,z+1/2;x,-y,-z+1/2;-x+1/2,y+1/2,-z;-x,-y,-z;x+1/2,y+1/2,-z+1/2;-x,y,z+1/2;x+1/2,-y+1/2,z"
	symLines[286]  = "x,y,z;-x+1/2,-y,z;x+1/2,-y+1/2,-z+1/2;-x,y+1/2,-z+1/2;-x,-y,-z;x+1/2,y,-z;-x+1/2,y+1/2,z+1/2;x,-y+1/2,z+1/2"
	symLines[287]  = "x,y,z;-x,-y+1/2,z+1/2;x+1/2,-y+1/2,-z+1/2;-x+1/2,y,-z;-x,-y,-z;x,y+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y,z"
	symLines[288]  = "x,y,z;-x+1/2,-y,z+1/2;x,-y+1/2,-z;-x+1/2,y+1/2,-z+1/2;-x,-y,-z;x+1/2,y,-z+1/2;-x,y+1/2,z;x+1/2,-y+1/2,z+1/2"
	symLines[289]  = "x,y,z;-x,-y+1/2,z;x+1/2,-y,-z+1/2;-x+1/2,y+1/2,-z+1/2;-x,-y,-z;x,y+1/2,-z;-x+1/2,y,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[290]  = "x,y,z;-x+1/2,-y,z+1/2;x+1/2,-y+1/2,-z;-x,y+1/2,-z+1/2;-x,-y,-z;x+1/2,y,-z+1/2;-x+1/2,y+1/2,z;x,-y+1/2,z+1/2"
	symLines[291]  = "x,y,z;-x,-y+1/2,z+1/2;x+1/2,-y,-z+1/2;-x+1/2,y+1/2,-z;-x,-y,-z;x,y+1/2,-z+1/2;-x+1/2,y,z+1/2;x+1/2,-y+1/2,z"
	symLines[292]  = "x,y,z;-x+1/2,-y,z+1/2;x+1/2,-y+1/2,-z+1/2;-x,y+1/2,-z;-x,-y,-z;x+1/2,y,-z+1/2;-x+1/2,y+1/2,z+1/2;x,-y+1/2,z"
	symLines[293]  = "x,y,z;-x,-y+1/2,z+1/2;x+1/2,-y,-z;-x+1/2,y+1/2,-z+1/2;-x,-y,-z;x,y+1/2,-z+1/2;-x+1/2,y,z;x+1/2,-y+1/2,z+1/2"
	symLines[294]  = "x,y,z;-x,-y,z+1/2;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z+1/2;-x,-y,-z;x,y,-z+1/2;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z+1/2"
	symLines[295]  = "x,y,z;-x+1/2,-y+1/2,z+1/2;x+1/2,-y,-z+1/2;-x,y+1/2,-z;-x,-y,-z;x+1/2,y+1/2,-z+1/2;-x+1/2,y,z+1/2;x,-y+1/2,z"
	symLines[296]  = "x,y,z;-x+1/2,-y+1/2,z+1/2;x+1/2,-y,-z;-x,y+1/2,-z+1/2;-x,-y,-z;x+1/2,y+1/2,-z+1/2;-x+1/2,y,z;x,-y+1/2,z+1/2"
	symLines[297]  = "x,y,z;-x,-y,z+1/2;x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z;-x,-y,-z;x,y,-z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z"
	symLines[298]  = "x,y,z;-x,-y,z+1/2;x,-y,-z;-x,y,-z+1/2;-x,-y,-z;x,y,-z+1/2;-x,y,z;x,-y,z+1/2;x+1/2,y+1/2,z;-x+1/2,-y+1/2,z+1/2;"
	symLines[298] += "x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z+1/2;-x+1/2,-y+1/2,-z;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z+1/2"
	symLines[299]  = "x,y,z;-x,-y,z+1/2;x,-y,-z+1/2;-x,y,-z;-x,-y,-z;x,y,-z+1/2;-x,y,z+1/2;x,-y,z;x+1/2,y+1/2,z;-x+1/2,-y+1/2,z+1/2;"
	symLines[299] += "x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z;-x+1/2,-y+1/2,-z;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z"
	symLines[300]  = "x,y,z;-x+1/2,-y,z;x+1/2,-y,-z;-x,y,-z;-x,-y,-z;x+1/2,y,-z;-x+1/2,y,z;x,-y,z;x,y+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;"
	symLines[300] += "x+1/2,-y+1/2,-z+1/2;-x,y+1/2,-z+1/2;-x,-y+1/2,-z+1/2;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2;x,-y+1/2,z+1/2"
	symLines[301]  = "x,y,z;-x,-y,z;x+1/2,-y,-z;-x+1/2,y,-z;-x,-y,-z;x,y,-z;-x+1/2,y,z;x+1/2,-y,z;x,y+1/2,z+1/2;-x,-y+1/2,z+1/2;"
	symLines[301] += "x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2;-x,-y+1/2,-z+1/2;x,y+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[302]  = "x,y,z;-x,-y,z;x,-y+1/2,-z;-x,y+1/2,-z;-x,-y,-z;x,y,-z;-x,y+1/2,z;x,-y+1/2,z;x+1/2,y,z+1/2;-x+1/2,-y,z+1/2;"
	symLines[302] += "x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2;-x+1/2,-y,-z+1/2;x+1/2,y,-z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[303]  = "x,y,z;-x,-y+1/2,z;x,-y,-z;-x,y+1/2,-z;-x,-y,-z;x,y+1/2,-z;-x,y,z;x,-y+1/2,z;x+1/2,y,z+1/2;-x+1/2,-y+1/2,z+1/2;"
	symLines[303] += "x+1/2,-y,-z+1/2;-x+1/2,y+1/2,-z+1/2;-x+1/2,-y,-z+1/2;x+1/2,y+1/2,-z+1/2;-x+1/2,y,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[304]  = "x,y,z;-x,-y+1/2,z+1/2;x,-y,-z;-x,y+1/2,-z+1/2;-x,-y,-z;x,y+1/2,-z+1/2;-x,y,z;x,-y+1/2,z+1/2;x+1/2,y+1/2,z;"
	symLines[304] += "-x+1/2,-y,z+1/2;x+1/2,-y+1/2,-z;-x+1/2,y,-z+1/2;-x+1/2,-y+1/2,-z;x+1/2,y,-z+1/2;-x+1/2,y+1/2,z;x+1/2,-y,z+1/2"
	symLines[305]  = "x,y,z;-x,-y+1/2,z+1/2;x,-y+1/2,-z+1/2;-x,y,-z;-x,-y,-z;x,y+1/2,-z+1/2;-x,y+1/2,z+1/2;x,-y,z;x+1/2,y+1/2,z;"
	symLines[305] += "-x+1/2,-y,z+1/2;x+1/2,-y,-z+1/2;-x+1/2,y+1/2,-z;-x+1/2,-y+1/2,-z;x+1/2,y,-z+1/2;-x+1/2,y,z+1/2;x+1/2,-y+1/2,z"
	symLines[306]  = "x,y,z;-x+1/2,-y,z+1/2;x+1/2,-y,-z+1/2;-x,y,-z;-x,-y,-z;x+1/2,y,-z+1/2;-x+1/2,y,z+1/2;x,-y,z;x,y+1/2,z+1/2;"
	symLines[306] += "-x+1/2,-y+1/2,z;x+1/2,-y+1/2,-z;-x,y+1/2,-z+1/2;-x,-y+1/2,-z+1/2;x+1/2,y+1/2,-z;-x+1/2,y+1/2,z;x,-y+1/2,z+1/2"
	symLines[307]  = "x,y,z;-x,-y,z;x+1/2,-y,-z+1/2;-x+1/2,y,-z+1/2;-x,-y,-z;x,y,-z;-x+1/2,y,z+1/2;x+1/2,-y,z+1/2;x,y+1/2,z+1/2;"
	symLines[307] += "-x,-y+1/2,z+1/2;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;-x,-y+1/2,-z+1/2;x,y+1/2,-z+1/2;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[308]  = "x,y,z;-x,-y,z;x,-y+1/2,-z+1/2;-x,y+1/2,-z+1/2;-x,-y,-z;x,y,-z;-x,y+1/2,z+1/2;x,-y+1/2,z+1/2;x+1/2,y,z+1/2;"
	symLines[308] += "-x+1/2,-y,z+1/2;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;-x+1/2,-y,-z+1/2;x+1/2,y,-z+1/2;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[309]  = "x,y,z;-x,-y+1/2,z+1/2;x,-y,-z;-x,y+1/2,-z+1/2;-x,-y,-z;x,y+1/2,-z+1/2;-x,y,z;x,-y+1/2,z+1/2;x+1/2,y,z+1/2;"
	symLines[309] += "-x+1/2,-y+1/2,z;x+1/2,-y,-z+1/2;-x+1/2,y+1/2,-z;-x+1/2,-y,-z+1/2;x+1/2,y+1/2,-z;-x+1/2,y,z+1/2;x+1/2,-y+1/2,z"
	symLines[310]  = "x,y,z;-x,-y,z;x,-y,-z;-x,y,-z;-x,-y,-z;x,y,-z;-x,y,z;x,-y,z;x+1/2,y+1/2,z;-x+1/2,-y+1/2,z;x+1/2,-y+1/2,-z;"
	symLines[310] += "-x+1/2,y+1/2,-z;-x+1/2,-y+1/2,-z;x+1/2,y+1/2,-z;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[311]  = "x,y,z;-x,-y,z;x,-y,-z;-x,y,-z;-x,-y,-z;x,y,-z;-x,y,z;x,-y,z;x,y+1/2,z+1/2;-x,-y+1/2,z+1/2;x,-y+1/2,-z+1/2;"
	symLines[311] += "-x,y+1/2,-z+1/2;-x,-y+1/2,-z+1/2;x,y+1/2,-z+1/2;-x,y+1/2,z+1/2;x,-y+1/2,z+1/2"
	symLines[312]  = "x,y,z;-x,-y,z;x,-y,-z;-x,y,-z;-x,-y,-z;x,y,-z;-x,y,z;x,-y,z;x+1/2,y,z+1/2;-x+1/2,-y,z+1/2;x+1/2,-y,-z+1/2;"
	symLines[312] += "-x+1/2,y,-z+1/2;-x+1/2,-y,-z+1/2;x+1/2,y,-z+1/2;-x+1/2,y,z+1/2;x+1/2,-y,z+1/2"
	symLines[313]  = "x,y,z;-x,-y,z;x,-y,-z+1/2;-x,y,-z+1/2;-x,-y,-z;x,y,-z;-x,y,z+1/2;x,-y,z+1/2;x+1/2,y+1/2,z;-x+1/2,-y+1/2,z;"
	symLines[313] += "x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2;-x+1/2,-y+1/2,-z;x+1/2,y+1/2,-z;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[314]  = "x,y,z;-x+1/2,-y,z;x,-y,-z;-x+1/2,y,-z;-x,-y,-z;x+1/2,y,-z;-x,y,z;x+1/2,-y,z;x,y+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;"
	symLines[314] += "x,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2;-x,-y+1/2,-z+1/2;x+1/2,y+1/2,-z+1/2;-x,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[315]  = "x,y,z;-x,-y+1/2,z;x,-y+1/2,-z;-x,y,-z;-x,-y,-z;x,y+1/2,-z;-x,y+1/2,z;x,-y,z;x+1/2,y,z+1/2;-x+1/2,-y+1/2,z+1/2;"
	symLines[315] += "x+1/2,-y+1/2,-z+1/2;-x+1/2,y,-z+1/2;-x+1/2,-y,-z+1/2;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y,z+1/2"
	symLines[316]  = "x,y,z;-x,-y+1/2,z;x,-y,-z;-x,y+1/2,-z;-x,-y,-z;x,y+1/2,-z;-x,y,z;x,-y+1/2,z;x+1/2,y+1/2,z;-x+1/2,-y,z;x+1/2,-y+1/2,-z;"
	symLines[316] += "-x+1/2,y,-z;-x+1/2,-y+1/2,-z;x+1/2,y,-z;-x+1/2,y+1/2,z;x+1/2,-y,z"
	symLines[317]  = "x,y,z;-x,-y+1/2,z;x,-y+1/2,-z;-x,y,-z;-x,-y,-z;x,y+1/2,-z;-x,y+1/2,z;x,-y,z;x+1/2,y+1/2,z;-x+1/2,-y,z;x+1/2,-y,-z;"
	symLines[317] += "-x+1/2,y+1/2,-z;-x+1/2,-y+1/2,-z;x+1/2,y,-z;-x+1/2,y,z;x+1/2,-y+1/2,z"
	symLines[318]  = "x,y,z;-x,-y,z+1/2;x,-y,-z+1/2;-x,y,-z;-x,-y,-z;x,y,-z+1/2;-x,y,z+1/2;x,-y,z;x,y+1/2,z+1/2;-x,-y+1/2,z;x,-y+1/2,-z;"
	symLines[318] += "-x,y+1/2,-z+1/2;-x,-y+1/2,-z+1/2;x,y+1/2,-z;-x,y+1/2,z;x,-y+1/2,z+1/2"
	symLines[319]  = "x,y,z;-x,-y,z;x,-y,-z+1/2;-x,y,-z+1/2;-x,-y,-z;x,y,-z;-x,y,z+1/2;x,-y,z+1/2;x,y+1/2,z+1/2;-x,-y+1/2,z+1/2;x,-y+1/2,-z;"
	symLines[319] += "-x,y+1/2,-z;-x,-y+1/2,-z+1/2;x,y+1/2,-z+1/2;-x,y+1/2,z;x,-y+1/2,z"
	symLines[320]  = "x,y,z;-x,-y,z;x,-y,-z+1/2;-x,y,-z+1/2;-x,-y,-z;x,y,-z;-x,y,z+1/2;x,-y,z+1/2;x+1/2,y,z+1/2;-x+1/2,-y,z+1/2;x+1/2,-y,-z;"
	symLines[320] += "-x+1/2,y,-z;-x+1/2,-y,-z+1/2;x+1/2,y,-z+1/2;-x+1/2,y,z;x+1/2,-y,z"
	symLines[321]  = "x,y,z;-x,-y,z+1/2;x,-y,-z;-x,y,-z+1/2;-x,-y,-z;x,y,-z+1/2;-x,y,z;x,-y,z+1/2;x+1/2,y,z+1/2;-x+1/2,-y,z;x+1/2,-y,-z+1/2;"
	symLines[321] += "-x+1/2,y,-z;-x+1/2,-y,-z+1/2;x+1/2,y,-z;-x+1/2,y,z+1/2;x+1/2,-y,z"
	symLines[322]  = symLine322
	symLines[323]  = "x,y,z;-x,-y+1/2,z;x,-y+1/2,-z+1/2;-x,y,-z+1/2;-x,-y,-z;x,y+1/2,-z;-x,y+1/2,z+1/2;x,-y,z+1/2;x+1/2,y+1/2,z;-x+1/2,-y,z;"
	symLines[323] += "x+1/2,-y,-z+1/2;-x+1/2,y+1/2,-z+1/2;-x+1/2,-y+1/2,-z;x+1/2,y,-z;-x+1/2,y,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[324]  = symLine322
	symLines[325]  = "x,y,z;-x,-y+1/2,z;x,-y,-z+1/2;-x,y+1/2,-z+1/2;-x,-y,-z;x,y+1/2,-z;-x,y,z+1/2;x,-y+1/2,z+1/2;x+1/2,y+1/2,z;-x+1/2,-y,z;"
	symLines[325] += "x+1/2,-y+1/2,-z+1/2;-x+1/2,y,-z+1/2;-x+1/2,-y+1/2,-z;x+1/2,y,-z;-x+1/2,y+1/2,z+1/2;x+1/2,-y,z+1/2"
	symLines[326]  = symLine326
	symLines[327]  = "x,y,z;-x+1/2,-y,z;x,-y,-z+1/2;-x+1/2,y,-z+1/2;-x,-y,-z;x+1/2,y,-z;-x,y,z+1/2;x+1/2,-y,z+1/2;x,y+1/2,z+1/2;"
	symLines[327] += "-x+1/2,-y+1/2,z+1/2;x,-y+1/2,-z;-x+1/2,y+1/2,-z;-x,-y+1/2,-z+1/2;x+1/2,y+1/2,-z+1/2;-x,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[328]  = symLine326
	symLines[329]  = "x,y,z;-x+1/2,-y,z+1/2;x,-y,-z+1/2;-x+1/2,y,-z;-x,-y,-z;x+1/2,y,-z+1/2;-x,y,z+1/2;x+1/2,-y,z;x,y+1/2,z+1/2;"
	symLines[329] += "-x+1/2,-y+1/2,z;x,-y+1/2,-z;-x+1/2,y+1/2,-z+1/2;-x,-y+1/2,-z+1/2;x+1/2,y+1/2,-z;-x,y+1/2,z;x+1/2,-y+1/2,z+1/2"
	symLines[330]  = symLine330
	symLines[331]  = "x,y,z;-x,-y+1/2,z+1/2;x,-y+1/2,-z;-x,y,-z+1/2;-x,-y,-z;x,y+1/2,-z+1/2;-x,y+1/2,z;x,-y,z+1/2;x+1/2,y,z+1/2;"
	symLines[331] += "-x+1/2,-y+1/2,z;x+1/2,-y+1/2,-z+1/2;-x+1/2,y,-z;-x+1/2,-y,-z+1/2;x+1/2,y+1/2,-z;-x+1/2,y+1/2,z+1/2;x+1/2,-y,z"
	symLines[332]  = symLine330
	symLines[333]  = "x,y,z;-x,-y+1/2,z;x,-y+1/2,-z+1/2;-x,y,-z+1/2;-x,-y,-z;x,y+1/2,-z;-x,y+1/2,z+1/2;x,-y,z+1/2;x+1/2,y,z+1/2;"
	symLines[333] += "-x+1/2,-y+1/2,z+1/2;x+1/2,-y+1/2,-z;-x+1/2,y,-z;-x+1/2,-y,-z+1/2;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z;x+1/2,-y,z"
	symLines[334]  = "x,y,z;-x,-y,z;x,-y,-z;-x,y,-z;-x,-y,-z;x,y,-z;-x,y,z;x,-y,z;x,y+1/2,z+1/2;-x,-y+1/2,z+1/2;x,-y+1/2,-z+1/2;"
	symLines[334] += "-x,y+1/2,-z+1/2;-x,-y+1/2,-z+1/2;x,y+1/2,-z+1/2;-x,y+1/2,z+1/2;x,-y+1/2,z+1/2;x+1/2,y,z+1/2;-x+1/2,-y,z+1/2;"
	symLines[334] += "x+1/2,-y,-z+1/2;-x+1/2,y,-z+1/2;-x+1/2,-y,-z+1/2;x+1/2,y,-z+1/2;-x+1/2,y,z+1/2;x+1/2,-y,z+1/2;x+1/2,y+1/2,z;"
	symLines[334] += "-x+1/2,-y+1/2,z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;-x+1/2,-y+1/2,-z;x+1/2,y+1/2,-z;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[335]  = "x,y,z;-x+1/4,-y+1/4,-z+1/4;-x,-y,z;x,-y,-z;-x,y,-z;x+1/4,y+1/4,-z+1/4;-x+1/4,y+1/4,z+1/4;x+1/4,-y+1/4,z+1/4;"
	symLines[335] += "x,y+1/2,z+1/2;-x+1/4,-y-1/4,-z-1/4;-x,-y+1/2,z+1/2;x,-y+1/2,-z+1/2;-x,y+1/2,-z+1/2;x+1/4,y-1/4,-z-1/4;-x+1/4,y-1/4,z-1/4;"
	symLines[335] += "x+1/4,-y-1/4,z-1/4;x+1/2,y,z+1/2;-x-1/4,-y+1/4,-z-1/4;-x+1/2,-y,z+1/2;x+1/2,-y,-z+1/2;-x+1/2,y,-z+1/2;x-1/4,y+1/4,-z-1/4;"
	symLines[335] += "-x-1/4,y+1/4,z-1/4;x-1/4,-y+1/4,z-1/4;x+1/2,y+1/2,z;-x-1/4,-y-1/4,-z+1/4;-x+1/2,-y+1/2,z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;"
	symLines[335] += "x-1/4,y-1/4,-z+1/4;-x-1/4,y-1/4,z+1/4;x-1/4,-y-1/4,z+1/4"
	symLines[336]  = "x,y,z;-x+1/4,-y+1/4,z;x,-y+1/4,-z+1/4;-x+1/4,y,-z+1/4;-x,-y,-z;x-1/4,y-1/4,-z;-x,y-1/4,z-1/4;x-1/4,-y,z-1/4;"
	symLines[336] += "x,y+1/2,z+1/2;-x+1/4,-y-1/4,z+1/2;x,-y-1/4,-z-1/4;-x+1/4,y+1/2,-z-1/4;-x,-y+1/2,-z+1/2;x-1/4,y+1/4,-z+1/2;-x,y+1/4,z+1/4;"
	symLines[336] += "x-1/4,-y+1/2,z+1/4;x+1/2,y,z+1/2;-x-1/4,-y+1/4,z+1/2;x+1/2,-y+1/4,-z-1/4;-x-1/4,y,-z-1/4;-x+1/2,-y,-z+1/2;"
	symLines[336] += "x+1/4,y-1/4,-z+1/2;-x+1/2,y-1/4,z+1/4;x+1/4,-y,z+1/4;x+1/2,y+1/2,z;-x-1/4,-y-1/4,z;x+1/2,-y-1/4,-z+1/4;"
	symLines[336] += "-x-1/4,y+1/2,-z+1/4;-x+1/2,-y+1/2,-z;x+1/4,y+1/4,-z;-x+1/2,y+1/4,z-1/4;x+1/4,-y+1/2,z-1/4"
	symLines[337]  = "x,y,z;-x,-y,z;x,-y,-z;-x,y,-z;-x,-y,-z;x,y,-z;-x,y,z;x,-y,z;x+1/2,y+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;x+1/2,-y+1/2,-z+1/2;"
	symLines[337] += "-x+1/2,y+1/2,-z+1/2;-x+1/2,-y+1/2,-z+1/2;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[338]  = "x,y,z;-x,-y,z;x,-y,-z+1/2;-x,y,-z+1/2;-x,-y,-z;x,y,-z;-x,y,z+1/2;x,-y,z+1/2;x+1/2,y+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;"
	symLines[338] += "x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;-x+1/2,-y+1/2,-z+1/2;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[339]  = "x,y,z;-x+1/2,-y,z;x,-y,-z;-x+1/2,y,-z;-x,-y,-z;x+1/2,y,-z;-x,y,z;x+1/2,-y,z;x+1/2,y+1/2,z+1/2;-x,-y+1/2,z+1/2;"
	symLines[339] += "x+1/2,-y+1/2,-z+1/2;-x,y+1/2,-z+1/2;-x+1/2,-y+1/2,-z+1/2;x,y+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2;x,-y+1/2,z+1/2"
	symLines[340]  = "x,y,z;-x,-y+1/2,z;x,-y+1/2,-z;-x,y,-z;-x,-y,-z;x,y+1/2,-z;-x,y+1/2,z;x,-y,z;x+1/2,y+1/2,z+1/2;-x+1/2,-y,z+1/2;"
	symLines[340] += "x+1/2,-y,-z+1/2;-x+1/2,y+1/2,-z+1/2;-x+1/2,-y+1/2,-z+1/2;x+1/2,y,-z+1/2;-x+1/2,y,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[341]  = "x,y,z;-x,-y+1/2,z;x,-y,-z+1/2;-x+1/2,y,-z;-x,-y,-z;x,y+1/2,-z;-x,y,z+1/2;x+1/2,-y,z;x+1/2,y+1/2,z+1/2;-x+1/2,-y,z+1/2;"
	symLines[341] += "x+1/2,-y+1/2,-z;-x,y+1/2,-z+1/2;-x+1/2,-y+1/2,-z+1/2;x+1/2,y,-z+1/2;-x+1/2,y+1/2,z;x,-y+1/2,z+1/2"
	symLines[342]  = "x,y,z;-x+1/2,-y,z;x,-y+1/2,-z;-x,y,-z+1/2;-x,-y,-z;x+1/2,y,-z;-x,y+1/2,z;x,-y,z+1/2;x+1/2,y+1/2,z+1/2;-x,-y+1/2,z+1/2;"
	symLines[342] += "x+1/2,-y,-z+1/2;-x+1/2,y+1/2,-z;-x+1/2,-y+1/2,-z+1/2;x,y+1/2,-z+1/2;-x+1/2,y,z+1/2;x+1/2,-y+1/2,z"
	symLines[343]  = "x,y,z;-x,-y+1/2,z;x,-y,-z;-x,y+1/2,-z;-x,-y,-z;x,y+1/2,-z;-x,y,z;x,-y+1/2,z;x+1/2,y+1/2,z+1/2;-x+1/2,-y,z+1/2;"
	symLines[343] += "x+1/2,-y+1/2,-z+1/2;-x+1/2,y,-z+1/2;-x+1/2,-y+1/2,-z+1/2;x+1/2,y,-z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y,z+1/2"
	symLines[344]  = "x,y,z;-x+1/2,-y,z;x+1/2,-y,-z;-x,y,-z;-x,-y,-z;x+1/2,y,-z;-x+1/2,y,z;x,-y,z;x+1/2,y+1/2,z+1/2;-x,-y+1/2,z+1/2;"
	symLines[344] += "x,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2;-x+1/2,-y+1/2,-z+1/2;x,y+1/2,-z+1/2;-x,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[345]  = "x,y,z;-x,-y,z+1/2;x,-y,-z+1/2;-x,y,-z;-x,-y,-z;x,y,-z+1/2;-x,y,z+1/2;x,-y,z;x+1/2,y+1/2,z+1/2;-x+1/2,-y+1/2,z;"
	symLines[345] += "x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z+1/2;-x+1/2,-y+1/2,-z+1/2;x+1/2,y+1/2,-z;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z+1/2"
	symLines[346]  = "x,y,z;-x,-y,z;x,-y+1/2,-z;-x,y+1/2,-z;-x,-y,-z;x,y,-z;-x,y+1/2,z;x,-y+1/2,z;x+1/2,y+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;"
	symLines[346] += "x+1/2,-y,-z+1/2;-x+1/2,y,-z+1/2;-x+1/2,-y+1/2,-z+1/2;x+1/2,y+1/2,-z+1/2;-x+1/2,y,z+1/2;x+1/2,-y,z+1/2"
	symLines[347]  = "x,y,z;-x,-y,z;x+1/2,-y,-z;-x+1/2,y,-z;-x,-y,-z;x,y,-z;-x+1/2,y,z;x+1/2,-y,z;x+1/2,y+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;"
	symLines[347] += "x,-y+1/2,-z+1/2;-x,y+1/2,-z+1/2;-x+1/2,-y+1/2,-z+1/2;x+1/2,y+1/2,-z+1/2;-x,y+1/2,z+1/2;x,-y+1/2,z+1/2"
	symLines[348]  = "x,y,z;-x,-y,z+1/2;x,-y,-z;-x,y,-z+1/2;-x,-y,-z;x,y,-z+1/2;-x,y,z;x,-y,z+1/2;x+1/2,y+1/2,z+1/2;-x+1/2,-y+1/2,z;"
	symLines[348] += "x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z;-x+1/2,-y+1/2,-z+1/2;x+1/2,y+1/2,-z;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z"
	// Tetragonal SG[75,142]  SG_idNum [349,429]   (81 idNums)
	symLines[349]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z"
	symLines[350]  = "x,y,z;-y,x,z+1/4;-x,-y,z+1/2;y,-x,z-1/4"
	symLines[351]  = "x,y,z;-y,x,z+1/2;-x,-y,z;y,-x,z+1/2"
	symLines[352]  = "x,y,z;-y,x,z-1/4;-x,-y,z+1/2;y,-x,z+1/4"
	symLines[353]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x+1/2,y+1/2,z+1/2;-y+1/2,x+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;y+1/2,-x+1/2,z+1/2"
	symLines[354]  = "x,y,z;-y,x+1/2,z+1/4;-x,-y,z;y,-x+1/2,z+1/4;x+1/2,y+1/2,z+1/2;-y+1/2,x,z-1/4;-x+1/2,-y+1/2,z+1/2;y+1/2,-x,z-1/4"
	symLines[355]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z"
	symLines[356]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;x+1/2,y+1/2,z+1/2;y+1/2,-x+1/2,-z+1/2;-x+1/2,-y+1/2,z+1/2;-y+1/2,x+1/2,-z+1/2"
	symLines[357]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;-x,-y,-z;y,-x,-z;x,y,-z;-y,x,-z"
	symLines[358]  = "x,y,z;-y,x,z+1/2;-x,-y,z;y,-x,z+1/2;-x,-y,-z;y,-x,-z+1/2;x,y,-z;-y,x,-z+1/2"
	symLines[359]  = "x,y,z;-x+1/2,-y+1/2,-z;-y+1/2,x+1/2,z;-x,-y,z;y+1/2,-x+1/2,z;y,-x,-z;-y,x,-z;x+1/2,y+1/2,-z"
	symLines[360]  = "x,y,z;-y+1/2,x,z;-x+1/2,-y+1/2,z;y,-x+1/2,z;-x,-y,-z;y+1/2,-x,-z;x+1/2,y+1/2,-z;-y,x+1/2,-z"
	symLines[361]  = "x,y,z;-x+1/2,-y+1/2,-z+1/2;-y+1/2,x+1/2,z+1/2;-x,-y,z;y+1/2,-x+1/2,z+1/2;y,-x,-z;-y,x,-z;x+1/2,y+1/2,-z+1/2"
	symLines[362]  = "x,y,z;-y,x+1/2,z+1/2;-x+1/2,-y+1/2,z;y+1/2,-x,z+1/2;-x,-y,-z;y,-x+1/2,-z+1/2;x+1/2,y+1/2,-z;-y+1/2,x,-z+1/2"
	symLines[363]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;-x,-y,-z;y,-x,-z;x,y,-z;-y,x,-z;x+1/2,y+1/2,z+1/2;-y+1/2,x+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;"
	symLines[363] += "y+1/2,-x+1/2,z+1/2;-x+1/2,-y+1/2,-z+1/2;y+1/2,-x+1/2,-z+1/2;x+1/2,y+1/2,-z+1/2;-y+1/2,x+1/2,-z+1/2"
	symLines[364]  = "x,y,z;-x,-y+1/2,-z+1/4;-y,x+1/2,z+1/4;-x,-y,z;y,-x+1/2,z+1/4;y,-x,-z;-y,x,-z;x,y+1/2,-z+1/4;x+1/2,y+1/2,z+1/2;"
	symLines[364] += "-x+1/2,-y,-z-1/4;-y+1/2,x,z-1/4;-x+1/2,-y+1/2,z+1/2;y+1/2,-x,z-1/4;y+1/2,-x+1/2,-z+1/2;-y+1/2,x+1/2,-z+1/2;x+1/2,y,-z-1/4"
	symLines[365]  = "x,y,z;-y-1/4,x+1/4,z+1/4;-x,-y+1/2,z;y+1/4,-x+1/4,z+1/4;-x,-y,-z;y+1/4,-x-1/4,-z-1/4;x,y+1/2,-z;-y-1/4,x-1/4,-z-1/4;"
	symLines[365] += "x+1/2,y+1/2,z+1/2;-y+1/4,x-1/4,z-1/4;-x+1/2,-y,z+1/2;y-1/4,-x-1/4,z-1/4;-x+1/2,-y+1/2,-z+1/2;y-1/4,-x+1/4,-z+1/4;"
	symLines[365] += "x+1/2,y,-z+1/2;-y+1/4,x+1/4,-z+1/4"
	symLines[366]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x,-y,-z;-x,y,-z;y,x,-z;-y,-x,-z"
	symLines[367]  = "x,y,z;-y+1/2,x+1/2,z;-x,-y,z;y+1/2,-x+1/2,z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;y,x,-z;-y,-x,-z"
	symLines[368]  = "x,y,z;-y,x,z+1/4;-x,-y,z+1/2;y,-x,z-1/4;x,-y,-z+1/2;-x,y,-z;y,x,-z-1/4;-y,-x,-z+1/4"
	symLines[369]  = "x,y,z;-y+1/2,x+1/2,z+1/4;-x,-y,z+1/2;y+1/2,-x+1/2,z-1/4;x+1/2,-y+1/2,-z-1/4;-x+1/2,y+1/2,-z+1/4;y,x,-z;-y,-x,-z+1/2"
	symLines[370]  = "x,y,z;-y,x,z+1/2;-x,-y,z;y,-x,z+1/2;x,-y,-z;-x,y,-z;y,x,-z+1/2;-y,-x,-z+1/2"
	symLines[371]  = "x,y,z;-y+1/2,x+1/2,z+1/2;-x,-y,z;y+1/2,-x+1/2,z+1/2;x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2;y,x,-z;-y,-x,-z"
	symLines[372]  = "x,y,z;-y,x,z-1/4;-x,-y,z+1/2;y,-x,z+1/4;x,-y,-z+1/2;-x,y,-z;y,x,-z+1/4;-y,-x,-z-1/4"
	symLines[373]  = "x,y,z;-y+1/2,x+1/2,z-1/4;-x,-y,z+1/2;y+1/2,-x+1/2,z+1/4;x+1/2,-y+1/2,-z+1/4;-x+1/2,y+1/2,-z-1/4;y,x,-z;-y,-x,-z+1/2"
	symLines[374]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x,-y,-z;-x,y,-z;y,x,-z;-y,-x,-z;x+1/2,y+1/2,z+1/2;-y+1/2,x+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;"
	symLines[374] += "y+1/2,-x+1/2,z+1/2;x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2;y+1/2,x+1/2,-z+1/2;-y+1/2,-x+1/2,-z+1/2"
	symLines[375]  = "x,y,z;-y,x+1/2,z+1/4;-x,-y,z;y,-x+1/2,z+1/4;x,-y+1/2,-z+1/4;-x,y+1/2,-z+1/4;y,x,-z;-y,-x,-z;x+1/2,y+1/2,z+1/2;"
	symLines[375] += "-y+1/2,x,z-1/4;-x+1/2,-y+1/2,z+1/2;y+1/2,-x,z-1/4;x+1/2,-y,-z-1/4;-x+1/2,y,-z-1/4;y+1/2,x+1/2,-z+1/2;-y+1/2,-x+1/2,-z+1/2"
	symLines[376]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;-x,y,z;x,-y,z;-y,-x,z;y,x,z"
	symLines[377]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z;-y+1/2,-x+1/2,z;y+1/2,x+1/2,z"
	symLines[378]  = "x,y,z;-y,x,z+1/2;-x,-y,z;y,-x,z+1/2;-x,y,z+1/2;x,-y,z+1/2;-y,-x,z;y,x,z"
	symLines[379]  = "x,y,z;-y+1/2,x+1/2,z+1/2;-x,-y,z;y+1/2,-x+1/2,z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2;-y,-x,z;y,x,z"
	symLines[380]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;-x,y,z+1/2;x,-y,z+1/2;-y,-x,z+1/2;y,x,z+1/2"
	symLines[381]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2;-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2"
	symLines[382]  = "x,y,z;-y,x,z+1/2;-x,-y,z;y,-x,z+1/2;-x,y,z;x,-y,z;-y,-x,z+1/2;y,x,z+1/2"
	symLines[383]  = "x,y,z;-y,x,z+1/2;-x,-y,z;y,-x,z+1/2;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z;-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2"
	symLines[384]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;-x,y,z;x,-y,z;-y,-x,z;y,x,z;x+1/2,y+1/2,z+1/2;-y+1/2,x+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;"
	symLines[384] += "y+1/2,-x+1/2,z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2;-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2"
	symLines[385]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;-x,y,z+1/2;x,-y,z+1/2;-y,-x,z+1/2;y,x,z+1/2;x+1/2,y+1/2,z+1/2;-y+1/2,x+1/2,z+1/2;"
	symLines[385] += "-x+1/2,-y+1/2,z+1/2;y+1/2,-x+1/2,z+1/2;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z;-y+1/2,-x+1/2,z;y+1/2,x+1/2,z"
	symLines[386]  = "x,y,z;-y,x+1/2,z+1/4;-x,-y,z;y,-x+1/2,z+1/4;-x,y,z;x,-y,z;-y,-x+1/2,z+1/4;y,x+1/2,z+1/4;x+1/2,y+1/2,z+1/2;-y+1/2,x,z-1/4;"
	symLines[386] += "-x+1/2,-y+1/2,z+1/2;y+1/2,-x,z-1/4;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2;-y+1/2,-x,z-1/4;y+1/2,x,z-1/4"
	symLines[387]  = "x,y,z;-y,x+1/2,z+1/4;-x,-y,z;y,-x+1/2,z+1/4;-x,y,z+1/2;x,-y,z+1/2;-y+1/2,-x,z+1/4;y+1/2,x,z+1/4;x+1/2,y+1/2,z+1/2;"
	symLines[387] += "-y+1/2,x,z-1/4;-x+1/2,-y+1/2,z+1/2;y+1/2,-x,z-1/4;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z;-y,-x+1/2,z-1/4;y,x+1/2,z-1/4"
	symLines[388]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;x,-y,-z;-x,y,-z;-y,-x,z;y,x,z"
	symLines[389]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;x,-y,-z+1/2;-x,y,-z+1/2;-y,-x,z+1/2;y,x,z+1/2"
	symLines[390]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;-y+1/2,-x+1/2,z;y+1/2,x+1/2,z"
	symLines[391]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2;-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2"
	symLines[392]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;y,x,-z;-y,-x,-z;-x,y,z;x,-y,z"
	symLines[393]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;y,x,-z+1/2;-y,-x,-z+1/2;-x,y,z+1/2;x,-y,z+1/2"
	symLines[394]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;y+1/2,x+1/2,-z;-y+1/2,-x+1/2,-z;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[395]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;y+1/2,x+1/2,-z+1/2;-y+1/2,-x+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[396]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;y,x,-z;-y,-x,-z;-x,y,z;x,-y,z;x+1/2,y+1/2,z+1/2;y+1/2,-x+1/2,-z+1/2;-x+1/2,-y+1/2,z+1/2;"
	symLines[396] += "-y+1/2,x+1/2,-z+1/2;y+1/2,x+1/2,-z+1/2;-y+1/2,-x+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[397]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;y,x,-z+1/2;-y,-x,-z+1/2;-x,y,z+1/2;x,-y,z+1/2;x+1/2,y+1/2,z+1/2;y+1/2,-x+1/2,-z+1/2;"
	symLines[397] += "-x+1/2,-y+1/2,z+1/2;-y+1/2,x+1/2,-z+1/2;y+1/2,x+1/2,-z;-y+1/2,-x+1/2,-z;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[398]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;x,-y,-z;-x,y,-z;-y,-x,z;y,x,z;x+1/2,y+1/2,z+1/2;y+1/2,-x+1/2,-z+1/2;-x+1/2,-y+1/2,z+1/2;"
	symLines[398] += "-y+1/2,x+1/2,-z+1/2;x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2;-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2"
	symLines[399]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;x,-y+1/2,-z+1/4;-x,y+1/2,-z+1/4;-y,-x+1/2,z+1/4;y,x+1/2,z+1/4;x+1/2,y+1/2,z+1/2;"
	symLines[399] += "y+1/2,-x+1/2,-z+1/2;-x+1/2,-y+1/2,z+1/2;-y+1/2,x+1/2,-z+1/2;x+1/2,-y,-z-1/4;-x+1/2,y,-z-1/4;-y+1/2,-x,z-1/4;y+1/2,x,z-1/4"
	symLines[400]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x,-y,-z;-x,y,-z;y,x,-z;-y,-x,-z;-x,-y,-z;y,-x,-z;x,y,-z;-y,x,-z;-x,y,z;x,-y,z;-y,-x,z;y,x,z"
	symLines[401]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x,-y,-z+1/2;-x,y,-z+1/2;y,x,-z+1/2;-y,-x,-z+1/2;-x,-y,-z;y,-x,-z;x,y,-z;-y,x,-z;-x,y,z+1/2;"
	symLines[401] += "x,-y,z+1/2;-y,-x,z+1/2;y,x,z+1/2"
	symLines[402]  = "x,y,z;-x+1/2,-y+1/2,-z;-y,x,z;-x,-y,z;y,-x,z;y+1/2,-x+1/2,-z;-y+1/2,x+1/2,-z;x,-y,-z;-x,y,-z;y,x,-z;-y,-x,-z;"
	symLines[402] += "x+1/2,y+1/2,-z;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z;-y+1/2,-x+1/2,z;y+1/2,x+1/2,z"
	symLines[403]  = "x,y,z;-y+1/2,x,z;-x+1/2,-y+1/2,z;y,-x+1/2,z;x,-y+1/2,-z;-x+1/2,y,-z;y,x,-z;-y+1/2,-x+1/2,-z;-x,-y,-z;y+1/2,-x,-z;"
	symLines[403] += "x+1/2,y+1/2,-z;-y,x+1/2,-z;-x,y+1/2,z;x+1/2,-y,z;-y,-x,z;y+1/2,x+1/2,z"
	symLines[404]  = "x,y,z;-x+1/2,-y+1/2,-z+1/2;-y,x,z;-x,-y,z;y,-x,z;y+1/2,-x+1/2,-z+1/2;-y+1/2,x+1/2,-z+1/2;x,-y,-z;-x,y,-z;y,x,-z;-y,-x,-z;"
	symLines[404] += "x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2;-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2"
	symLines[405]  = "x,y,z;-y+1/2,x,z;-x+1/2,-y+1/2,z;y,-x+1/2,z;x,-y+1/2,-z+1/2;-x+1/2,y,-z+1/2;y,x,-z+1/2;-y+1/2,-x+1/2,-z+1/2;-x,-y,-z;"
	symLines[405] += "y+1/2,-x,-z;x+1/2,y+1/2,-z;-y,x+1/2,-z;-x,y+1/2,z+1/2;x+1/2,-y,z+1/2;-y,-x,z+1/2;y+1/2,x+1/2,z+1/2"
	symLines[406]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;y+1/2,x+1/2,-z;-y+1/2,-x+1/2,-z;-x,-y,-z;y,-x,-z;x,y,-z;"
	symLines[406] += "-y,x,-z;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z;-y+1/2,-x+1/2,z;y+1/2,x+1/2,z"
	symLines[407]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2;y+1/2,x+1/2,-z+1/2;-y+1/2,-x+1/2,-z+1/2;-x,-y,-z;"
	symLines[407] += "y,-x,-z;x,y,-z;-y,x,-z;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2;-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2"
	symLines[408]  = "x,y,z;-x+1/2,-y+1/2,-z;-y+1/2,x+1/2,z;-x,-y,z;y+1/2,-x+1/2,z;y,-x,-z;-y,x,-z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;y,x,-z;"
	symLines[408] += "-y,-x,-z;x+1/2,y+1/2,-z;-x,y,z;x,-y,z;-y+1/2,-x+1/2,z;y+1/2,x+1/2,z"
	symLines[409]  = "x,y,z;-y+1/2,x,z;-x+1/2,-y+1/2,z;y,-x+1/2,z;x+1/2,-y,-z;-x,y+1/2,-z;y+1/2,x+1/2,-z;-y,-x,-z;-x,-y,-z;y+1/2,-x,-z;"
	symLines[409] += "x+1/2,y+1/2,-z;-y,x+1/2,-z;-x+1/2,y,z;x,-y+1/2,z;-y+1/2,-x+1/2,z;y,x,z"
	symLines[410]  = "x,y,z;-x+1/2,-y+1/2,-z;-y+1/2,x+1/2,z;-x,-y,z;y+1/2,-x+1/2,z;y,-x,-z;-y,x,-z;x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2;"
	symLines[410] += "y,x,-z+1/2;-y,-x,-z+1/2;x+1/2,y+1/2,-z;-x,y,z+1/2;x,-y,z+1/2;-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2"
	symLines[411]  = "x,y,z;-y+1/2,x,z;-x+1/2,-y+1/2,z;y,-x+1/2,z;x+1/2,-y,-z+1/2;-x,y+1/2,-z+1/2;y+1/2,x+1/2,-z+1/2;-y,-x,-z+1/2;-x,-y,-z;"
	symLines[411] += "y+1/2,-x,-z;x+1/2,y+1/2,-z;-y,x+1/2,-z;-x+1/2,y,z+1/2;x,-y+1/2,z+1/2;-y+1/2,-x+1/2,z+1/2;y,x,z+1/2"
	symLines[412]  = "x,y,z;-y,x,z+1/2;-x,-y,z;y,-x,z+1/2;x,-y,-z;-x,y,-z;y,x,-z+1/2;-y,-x,-z+1/2;-x,-y,-z;y,-x,-z+1/2;x,y,-z;-y,x,-z+1/2;"
	symLines[412] += "-x,y,z;x,-y,z;-y,-x,z+1/2;y,x,z+1/2"
	symLines[413]  = "x,y,z;-y,x,z+1/2;-x,-y,z;y,-x,z+1/2;x,-y,-z+1/2;-x,y,-z+1/2;y,x,-z;-y,-x,-z;-x,-y,-z;y,-x,-z+1/2;x,y,-z;-y,x,-z+1/2;"
	symLines[413] += "-x,y,z+1/2;x,-y,z+1/2;-y,-x,z;y,x,z"
	symLines[414]  = "x,y,z;-x+1/2,-y+1/2,-z+1/2;-y+1/2,x+1/2,z+1/2;-x,-y,z;y+1/2,-x+1/2,z+1/2;y,-x,-z;-y,x,-z;x,-y,-z+1/2;-x,y,-z+1/2;"
	symLines[414] += "y+1/2,x+1/2,-z;-y+1/2,-x+1/2,-z;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z;-y,-x,z+1/2;y,x,z+1/2"
	symLines[415]  = "x,y,z;-y+1/2,x,z+1/2;-x+1/2,-y+1/2,z;y,-x+1/2,z+1/2;x,-y+1/2,-z;-x+1/2,y,-z;y,x,-z+1/2;-y+1/2,-x+1/2,-z+1/2;-x,-y,-z;"
	symLines[415] += "y+1/2,-x,-z+1/2;x+1/2,y+1/2,-z;-y,x+1/2,-z+1/2;-x,y+1/2,z;x+1/2,-y,z;-y,-x,z+1/2;y+1/2,x+1/2,z+1/2"
	symLines[416]  = "x,y,z;-x+1/2,-y+1/2,-z+1/2;-y+1/2,x+1/2,z+1/2;-x,-y,z;y+1/2,-x+1/2,z+1/2;y,-x,-z;-y,x,-z;x,-y,-z;-x,y,-z;"
	symLines[416] += "y+1/2,x+1/2,-z+1/2;-y+1/2,-x+1/2,-z+1/2;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2;-y,-x,z;y,x,z"
	symLines[417]  = "x,y,z;-y+1/2,x,z+1/2;-x+1/2,-y+1/2,z;y,-x+1/2,z+1/2;x,-y+1/2,-z+1/2;-x+1/2,y,-z+1/2;y,x,-z;-y+1/2,-x+1/2,-z;-x,-y,-z;"
	symLines[417] += "y+1/2,-x,-z+1/2;x+1/2,y+1/2,-z;-y,x+1/2,-z+1/2;-x,y+1/2,z+1/2;x+1/2,-y,z+1/2;-y,-x,z;y+1/2,x+1/2,z"
	symLines[418]  = "x,y,z;-y,x,z+1/2;-x,-y,z;y,-x,z+1/2;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;y+1/2,x+1/2,-z+1/2;-y+1/2,-x+1/2,-z+1/2;-x,-y,-z;"
	symLines[418] += "y,-x,-z+1/2;x,y,-z;-y,x,-z+1/2;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z;-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2"
	symLines[419]  = "x,y,z;-y+1/2,x+1/2,z+1/2;-x,-y,z;y+1/2,-x+1/2,z+1/2;x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2;y,x,-z;-y,-x,-z;-x,-y,-z;"
	symLines[419] += "y+1/2,-x+1/2,-z+1/2;x,y,-z;-y+1/2,x+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2;-y,-x,z;y,x,z"
	symLines[420]  = "x,y,z;-x+1/2,-y+1/2,-z+1/2;-y+1/2,x+1/2,z+1/2;-x,-y,z;y+1/2,-x+1/2,z+1/2;y,-x,-z;-y,x,-z;x+1/2,-y+1/2,-z+1/2;"
	symLines[420] += "-x+1/2,y+1/2,-z+1/2;y,x,-z;-y,-x,-z;x+1/2,y+1/2,-z+1/2;-x,y,z;x,-y,z;-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2"
	symLines[421]  = "x,y,z;-y+1/2,x,z+1/2;-x+1/2,-y+1/2,z;y,-x+1/2,z+1/2;x+1/2,-y,-z;-x,y+1/2,-z;y+1/2,x+1/2,-z+1/2;-y,-x,-z+1/2;-x,-y,-z;"
	symLines[421] += "y+1/2,-x,-z+1/2;x+1/2,y+1/2,-z;-y,x+1/2,-z+1/2;-x+1/2,y,z;x,-y+1/2,z;-y+1/2,-x+1/2,z+1/2;y,x,z+1/2"
	symLines[422]  = "x,y,z;-x+1/2,-y+1/2,-z+1/2;-y+1/2,x+1/2,z+1/2;-x,-y,z;y+1/2,-x+1/2,z+1/2;y,-x,-z;-y,x,-z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;"
	symLines[422] += "y,x,-z+1/2;-y,-x,-z+1/2;x+1/2,y+1/2,-z+1/2;-x,y,z+1/2;x,-y,z+1/2;-y+1/2,-x+1/2,z;y+1/2,x+1/2,z"
	symLines[423]  = "x,y,z;-y+1/2,x,z+1/2;-x+1/2,-y+1/2,z;y,-x+1/2,z+1/2;x+1/2,-y,-z+1/2;-x,y+1/2,-z+1/2;y+1/2,x+1/2,-z;-y,-x,-z;-x,-y,-z;"
	symLines[423] += "y+1/2,-x,-z+1/2;x+1/2,y+1/2,-z;-y,x+1/2,-z+1/2;-x+1/2,y,z+1/2;x,-y+1/2,z+1/2;-y+1/2,-x+1/2,z;y,x,z"
	symLines[424]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x,-y,-z;-x,y,-z;y,x,-z;-y,-x,-z;-x,-y,-z;y,-x,-z;x,y,-z;-y,x,-z;-x,y,z;x,-y,z;-y,-x,z;y,x,z;"
	symLines[424] += "x+1/2,y+1/2,z+1/2;-y+1/2,x+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;y+1/2,-x+1/2,z+1/2;x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2;"
	symLines[424] += "y+1/2,x+1/2,-z+1/2;-y+1/2,-x+1/2,-z+1/2;-x+1/2,-y+1/2,-z+1/2;y+1/2,-x+1/2,-z+1/2;x+1/2,y+1/2,-z+1/2;-y+1/2,x+1/2,-z+1/2;"
	symLines[424] += "-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2;-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2"
	symLines[425]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x,-y,-z+1/2;-x,y,-z+1/2;y,x,-z+1/2;-y,-x,-z+1/2;-x,-y,-z;y,-x,-z;x,y,-z;-y,x,-z;-x,y,z+1/2;"
	symLines[425] += "x,-y,z+1/2;-y,-x,z+1/2;y,x,z+1/2;x+1/2,y+1/2,z+1/2;-y+1/2,x+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;y+1/2,-x+1/2,z+1/2;"
	symLines[425] += "x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;y+1/2,x+1/2,-z;-y+1/2,-x+1/2,-z;-x+1/2,-y+1/2,-z+1/2;y+1/2,-x+1/2,-z+1/2;"
	symLines[425] += "x+1/2,y+1/2,-z+1/2;-y+1/2,x+1/2,-z+1/2;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z;-y+1/2,-x+1/2,z;y+1/2,x+1/2,z"
	symLines[426]  = "x,y,z;-x,-y+1/2,-z+1/4;-y,x+1/2,z+1/4;-x,-y,z;y,-x+1/2,z+1/4;y,-x,-z;-y,x,-z;x,-y+1/2,-z+1/4;-x,y+1/2,-z+1/4;y,x,-z;"
	symLines[426] += "-y,-x,-z;x,y+1/2,-z+1/4;-x,y,z;x,-y,z;-y,-x+1/2,z+1/4;y,x+1/2,z+1/4;x+1/2,y+1/2,z+1/2;-x+1/2,-y,-z-1/4;-y+1/2,x,z-1/4;"
	symLines[426] += "-x+1/2,-y+1/2,z+1/2;y+1/2,-x,z-1/4;y+1/2,-x+1/2,-z+1/2;-y+1/2,x+1/2,-z+1/2;x+1/2,-y,-z-1/4;-x+1/2,y,-z-1/4;"
	symLines[426] += "y+1/2,x+1/2,-z+1/2;-y+1/2,-x+1/2,-z+1/2;x+1/2,y,-z-1/4;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2;-y+1/2,-x,z-1/4;y+1/2,x,z-1/4"
	symLines[427]  = "x,y,z;-y+1/4,x-1/4,z+1/4;-x,-y+1/2,z;y+1/4,-x+1/4,z-1/4;x,-y,-z;-x,y+1/2,-z;y+1/4,x-1/4,-z+1/4;-y+1/4,-x+1/4,-z-1/4;"
	symLines[427] += "-x,-y,-z;y-1/4,-x+1/4,-z-1/4;x,y+1/2,-z;-y-1/4,x-1/4,-z+1/4;-x,y,z;x,-y+1/2,z;-y-1/4,-x+1/4,z-1/4;y-1/4,x-1/4,z+1/4;"
	symLines[427] += "x+1/2,y+1/2,z+1/2;-y-1/4,x+1/4,z-1/4;-x+1/2,-y,z+1/2;y-1/4,-x-1/4,z+1/4;x+1/2,-y+1/2,-z+1/2;-x+1/2,y,-z+1/2;"
	symLines[427] += "y-1/4,x+1/4,-z-1/4;-y-1/4,-x-1/4,-z+1/4;-x+1/2,-y+1/2,-z+1/2;y+1/4,-x-1/4,-z+1/4;x+1/2,y,-z+1/2;-y+1/4,x+1/4,-z-1/4;"
	symLines[427] += "-x+1/2,y+1/2,z+1/2;x+1/2,-y,z+1/2;-y+1/4,-x-1/4,z+1/4;y+1/4,x+1/4,z-1/4"
	symLines[428]  = "x,y,z;-x,-y+1/2,-z+1/4;-y,x+1/2,z+1/4;-x,-y,z;y,-x+1/2,z+1/4;y,-x,-z;-y,x,-z;x+1/2,-y,-z+1/4;-x+1/2,y,-z+1/4;y,x,-z+1/2;"
	symLines[428] += "-y,-x,-z+1/2;x,y+1/2,-z+1/4;-x,y,z+1/2;x,-y,z+1/2;-y+1/2,-x,z+1/4;y+1/2,x,z+1/4;x+1/2,y+1/2,z+1/2;-x+1/2,-y,-z-1/4;"
	symLines[428] += "-y+1/2,x,z-1/4;-x+1/2,-y+1/2,z+1/2;y+1/2,-x,z-1/4;y+1/2,-x+1/2,-z+1/2;-y+1/2,x+1/2,-z+1/2;x,-y+1/2,-z-1/4;"
	symLines[428] += "-x,y+1/2,-z-1/4;y+1/2,x+1/2,-z;-y+1/2,-x+1/2,-z;x+1/2,y,-z-1/4;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z;-y,-x+1/2,z-1/4;y,x+1/2,z-1/4"
	symLines[429]  = "x,y,z;-y+1/4,x-1/4,z+1/4;-x,-y+1/2,z;y+1/4,-x+1/4,z-1/4;x,-y,-z+1/2;-x+1/2,y,-z;y-1/4,x+1/4,-z+1/4;-y+1/4,-x+1/4,-z+1/4;"
	symLines[429] += "-x,-y,-z;y-1/4,-x+1/4,-z-1/4;x,y+1/2,-z;-y-1/4,x-1/4,-z+1/4;-x,y,z+1/2;x+1/2,-y,z;-y+1/4,-x-1/4,z-1/4;y-1/4,x-1/4,z-1/4;"
	symLines[429] += "x+1/2,y+1/2,z+1/2;-y-1/4,x+1/4,z-1/4;-x+1/2,-y,z+1/2;y-1/4,-x-1/4,z+1/4;x+1/2,-y+1/2,-z;-x,y+1/2,-z+1/2;"
	symLines[429] += "y+1/4,x-1/4,-z-1/4;-y-1/4,-x-1/4,-z-1/4;-x+1/2,-y+1/2,-z+1/2;y+1/4,-x-1/4,-z+1/4;x+1/2,y,-z+1/2;-y+1/4,x+1/4,-z-1/4;"
	symLines[429] += "-x+1/2,y+1/2,z;x,-y+1/2,z+1/2;-y-1/4,-x+1/4,z+1/4;y+1/4,x+1/4,z+1/4"
	// Trigonal SG[143,167]  SG_idNum [430,461]   (32 idNums)
	symLines[430]  = "x,y,z;-y,x-y,z;-x+y,-x,z"
	symLines[431]  = "x,y,z;-y,x-y,z+1/3;-x+y,-x,z-1/3"
	symLines[432]  = "x,y,z;-y,x-y,z-1/3;-x+y,-x,z+1/3"
	symLines[433]  = "x,y,z;-y,x-y,z;-x+y,-x,z;x-1/3,y+1/3,z+1/3;-y-1/3,x-y+1/3,z+1/3;-x+y-1/3,-x+1/3,z+1/3;x+1/3,y-1/3,z-1/3;"
	symLines[433] += "-y+1/3,x-y-1/3,z-1/3;-x+y+1/3,-x-1/3,z-1/3"
	symLines[434]  = "x,y,z;z,x,y;y,z,x"
	symLines[435]  = "x,y,z;-y,x-y,z;-x+y,-x,z;-x,-y,-z;y,-x+y,-z;x-y,x,-z"
	symLines[436]  = "x,y,z;-y,x-y,z;-x+y,-x,z;-x,-y,-z;y,-x+y,-z;x-y,x,-z;x-1/3,y+1/3,z+1/3;-y-1/3,x-y+1/3,z+1/3;-x+y-1/3,-x+1/3,z+1/3;"
	symLines[436] += "-x-1/3,-y+1/3,-z+1/3;y-1/3,-x+y+1/3,-z+1/3;x-y-1/3,x+1/3,-z+1/3;x+1/3,y-1/3,z-1/3;-y+1/3,x-y-1/3,z-1/3;"
	symLines[436] += "-x+y+1/3,-x-1/3,z-1/3;-x+1/3,-y-1/3,-z-1/3;y+1/3,-x+y-1/3,-z-1/3;x-y+1/3,x-1/3,-z-1/3"
	symLines[437]  = "x,y,z;z,x,y;y,z,x;-x,-y,-z;-z,-x,-y;-y,-z,-x"
	symLines[438]  = "x,y,z;-y,x-y,z;-x+y,-x,z;-y,-x,-z;-x+y,y,-z;x,x-y,-z"
	symLines[439]  = "x,y,z;-y,x-y,z;-x+y,-x,z;x-y,-y,-z;-x,-x+y,-z;y,x,-z"
	symLines[440]  = "x,y,z;-y,x-y,z+1/3;-x+y,-x,z-1/3;-y,-x,-z-1/3;-x+y,y,-z+1/3;x,x-y,-z"
	symLines[441]  = "x,y,z;-y,x-y,z+1/3;-x+y,-x,z-1/3;x-y,-y,-z-1/3;-x,-x+y,-z+1/3;y,x,-z"
	symLines[442]  = "x,y,z;-y,x-y,z-1/3;-x+y,-x,z+1/3;-y,-x,-z+1/3;-x+y,y,-z-1/3;x,x-y,-z"
	symLines[443]  = "x,y,z;-y,x-y,z-1/3;-x+y,-x,z+1/3;x-y,-y,-z+1/3;-x,-x+y,-z-1/3;y,x,-z"
	symLines[444]  = "x,y,z;-y,x-y,z;-x+y,-x,z;x-y,-y,-z;-x,-x+y,-z;y,x,-z;x-1/3,y+1/3,z+1/3;-y-1/3,x-y+1/3,z+1/3;-x+y-1/3,-x+1/3,z+1/3;"
	symLines[444] += "x-y-1/3,-y+1/3,-z+1/3;-x-1/3,-x+y+1/3,-z+1/3;y-1/3,x+1/3,-z+1/3;x+1/3,y-1/3,z-1/3;-y+1/3,x-y-1/3,z-1/3;"
	symLines[444] += "-x+y+1/3,-x-1/3,z-1/3;x-y+1/3,-y-1/3,-z-1/3;-x+1/3,-x+y-1/3,-z-1/3;y+1/3,x-1/3,-z-1/3"
	symLines[445]  = "x,y,z;z,x,y;y,z,x;-y,-x,-z;-x,-z,-y;-z,-y,-x"
	symLines[446]  = "x,y,z;-y,x-y,z;-x+y,-x,z;-x+y,y,z;x,x-y,z;-y,-x,z"
	symLines[447]  = "x,y,z;-y,x-y,z;-x+y,-x,z;y,x,z;x-y,-y,z;-x,-x+y,z"
	symLines[448]  = "x,y,z;-y,x-y,z;-x+y,-x,z;-x+y,y,z+1/2;x,x-y,z+1/2;-y,-x,z+1/2"
	symLines[449]  = "x,y,z;-y,x-y,z;-x+y,-x,z;y,x,z+1/2;x-y,-y,z+1/2;-x,-x+y,z+1/2"
	symLines[450]  = "x,y,z;-y,x-y,z;-x+y,-x,z;-x+y,y,z;x,x-y,z;-y,-x,z;x-1/3,y+1/3,z+1/3;-y-1/3,x-y+1/3,z+1/3;-x+y-1/3,-x+1/3,z+1/3;"
	symLines[450] += "-x+y-1/3,y+1/3,z+1/3;x-1/3,x-y+1/3,z+1/3;-y-1/3,-x+1/3,z+1/3;x+1/3,y-1/3,z-1/3;-y+1/3,x-y-1/3,z-1/3;"
	symLines[450] += "-x+y+1/3,-x-1/3,z-1/3;-x+y+1/3,y-1/3,z-1/3;x+1/3,x-y-1/3,z-1/3;-y+1/3,-x-1/3,z-1/3"
	symLines[451]  = "x,y,z;z,x,y;y,z,x;y,x,z;x,z,y;z,y,x"
	symLines[452]  = "x,y,z;-y,x-y,z;-x+y,-x,z;-x+y,y,z+1/2;x,x-y,z+1/2;-y,-x,z+1/2;x-1/3,y+1/3,z+1/3;-y-1/3,x-y+1/3,z+1/3;"
	symLines[452] += "-x+y-1/3,-x+1/3,z+1/3;-x+y-1/3,y+1/3,z-1/6;x-1/3,x-y+1/3,z-1/6;-y-1/3,-x+1/3,z-1/6;x+1/3,y-1/3,z-1/3;"
	symLines[452] += "-y+1/3,x-y-1/3,z-1/3;-x+y+1/3,-x-1/3,z-1/3;-x+y+1/3,y-1/3,z+1/6;x+1/3,x-y-1/3,z+1/6;-y+1/3,-x-1/3,z+1/6"
	symLines[453]  = "x,y,z;z,x,y;y,z,x;y+1/2,x+1/2,z+1/2;x+1/2,z+1/2,y+1/2;z+1/2,y+1/2,x+1/2"
	symLines[454]  = "x,y,z;-y,x-y,z;-x+y,-x,z;-y,-x,-z;-x+y,y,-z;x,x-y,-z;-x,-y,-z;y,-x+y,-z;x-y,x,-z;y,x,z;x-y,-y,z;-x,-x+y,z"
	symLines[455]  = "x,y,z;-y,x-y,z;-x+y,-x,z;-y,-x,-z+1/2;-x+y,y,-z+1/2;x,x-y,-z+1/2;-x,-y,-z;y,-x+y,-z;x-y,x,-z;y,x,z+1/2;x-y,-y,z+1/2;-x,-x+y,z+1/2"
	symLines[456]  = "x,y,z;-y,x-y,z;-x+y,-x,z;x-y,-y,-z;-x,-x+y,-z;y,x,-z;-x,-y,-z;y,-x+y,-z;x-y,x,-z;-x+y,y,z;x,x-y,z;-y,-x,z"
	symLines[457]  = "x,y,z;-y,x-y,z;-x+y,-x,z;x-y,-y,-z+1/2;-x,-x+y,-z+1/2;y,x,-z+1/2;-x,-y,-z;y,-x+y,-z;x-y,x,-z;-x+y,y,z+1/2;x,x-y,z+1/2;-y,-x,z+1/2"
	symLines[458]  = "x,y,z;-y,x-y,z;-x+y,-x,z;x-y,-y,-z;-x,-x+y,-z;y,x,-z;-x,-y,-z;y,-x+y,-z;x-y,x,-z;-x+y,y,z;x,x-y,z;-y,-x,z;"
	symLines[458] += "x-1/3,y+1/3,z+1/3;-y-1/3,x-y+1/3,z+1/3;-x+y-1/3,-x+1/3,z+1/3;x-y-1/3,-y+1/3,-z+1/3;-x-1/3,-x+y+1/3,-z+1/3;"
	symLines[458] += "y-1/3,x+1/3,-z+1/3;-x-1/3,-y+1/3,-z+1/3;y-1/3,-x+y+1/3,-z+1/3;x-y-1/3,x+1/3,-z+1/3;-x+y-1/3,y+1/3,z+1/3;"
	symLines[458] += "x-1/3,x-y+1/3,z+1/3;-y-1/3,-x+1/3,z+1/3;x+1/3,y-1/3,z-1/3;-y+1/3,x-y-1/3,z-1/3;-x+y+1/3,-x-1/3,z-1/3;"
	symLines[458] += "x-y+1/3,-y-1/3,-z-1/3;-x+1/3,-x+y-1/3,-z-1/3;y+1/3,x-1/3,-z-1/3;-x+1/3,-y-1/3,-z-1/3;y+1/3,-x+y-1/3,-z-1/3;"
	symLines[458] += "x-y+1/3,x-1/3,-z-1/3;-x+y+1/3,y-1/3,z-1/3;x+1/3,x-y-1/3,z-1/3;-y+1/3,-x-1/3,z-1/3"
	symLines[459]  = "x,y,z;z,x,y;y,z,x;-y,-x,-z;-x,-z,-y;-z,-y,-x;-x,-y,-z;-z,-x,-y;-y,-z,-x;y,x,z;x,z,y;z,y,x"
	symLines[460]  = "x,y,z;-y,x-y,z;-x+y,-x,z;x-y,-y,-z+1/2;-x,-x+y,-z+1/2;y,x,-z+1/2;-x,-y,-z;y,-x+y,-z;x-y,x,-z;-x+y,y,z+1/2;x,x-y,z+1/2;"
	symLines[460] += "-y,-x,z+1/2;x-1/3,y+1/3,z+1/3;-y-1/3,x-y+1/3,z+1/3;-x+y-1/3,-x+1/3,z+1/3;x-y-1/3,-y+1/3,-z-1/6;-x-1/3,-x+y+1/3,-z-1/6;"
	symLines[460] += "y-1/3,x+1/3,-z-1/6;-x-1/3,-y+1/3,-z+1/3;y-1/3,-x+y+1/3,-z+1/3;x-y-1/3,x+1/3,-z+1/3;-x+y-1/3,y+1/3,z-1/6;"
	symLines[460] += "x-1/3,x-y+1/3,z-1/6;-y-1/3,-x+1/3,z-1/6;x+1/3,y-1/3,z-1/3;-y+1/3,x-y-1/3,z-1/3;-x+y+1/3,-x-1/3,z-1/3;"
	symLines[460] += "x-y+1/3,-y-1/3,-z+1/6;-x+1/3,-x+y-1/3,-z+1/6;y+1/3,x-1/3,-z+1/6;-x+1/3,-y-1/3,-z-1/3;y+1/3,-x+y-1/3,-z-1/3;"
	symLines[460] += "x-y+1/3,x-1/3,-z-1/3;-x+y+1/3,y-1/3,z+1/6;x+1/3,x-y-1/3,z+1/6;-y+1/3,-x-1/3,z+1/6"
	symLines[461]  = "x,y,z;z,x,y;y,z,x;-y+1/2,-x+1/2,-z+1/2;-x+1/2,-z+1/2,-y+1/2;-z+1/2,-y+1/2,-x+1/2;-x,-y,-z;-z,-x,-y;-y,-z,-x;"
	symLines[461] += "y+1/2,x+1/2,z+1/2;x+1/2,z+1/2,y+1/2;z+1/2,y+1/2,x+1/2"
	// Hexagonal SG[168,194]  SG_idNum [462,488]   (25 idNums)
	symLines[462]  = "x,y,z;x-y,x,z;-y,x-y,z;-x,-y,z;-x+y,-x,z;y,-x+y,z"
	symLines[463]  = "x,y,z;x-y,x,z+1/6;-y,x-y,z+1/3;-x,-y,z+1/2;-x+y,-x,z-1/3;y,-x+y,z-1/6"
	symLines[464]  = "x,y,z;x-y,x,z-1/6;-y,x-y,z-1/3;-x,-y,z+1/2;-x+y,-x,z+1/3;y,-x+y,z+1/6"
	symLines[465]  = "x,y,z;x-y,x,z+1/3;-y,x-y,z-1/3;-x,-y,z;-x+y,-x,z+1/3;y,-x+y,z-1/3"
	symLines[466]  = "x,y,z;x-y,x,z-1/3;-y,x-y,z+1/3;-x,-y,z;-x+y,-x,z-1/3;y,-x+y,z+1/3"
	symLines[467]  = "x,y,z;x-y,x,z+1/2;-y,x-y,z;-x,-y,z+1/2;-x+y,-x,z;y,-x+y,z+1/2"
	symLines[468]  = "x,y,z;-x+y,-x,-z;-y,x-y,z;x,y,-z;-x+y,-x,z;-y,x-y,-z"
	symLines[469]  = "x,y,z;x-y,x,z;-y,x-y,z;-x,-y,z;-x+y,-x,z;y,-x+y,z;-x,-y,-z;-x+y,-x,-z;y,-x+y,-z;x,y,-z;x-y,x,-z;-y,x-y,-z"
	symLines[470]  = "x,y,z;x-y,x,z+1/2;-y,x-y,z;-x,-y,z+1/2;-x+y,-x,z;y,-x+y,z+1/2;-x,-y,-z;-x+y,-x,-z+1/2;y,-x+y,-z;x,y,-z+1/2;x-y,x,-z;-y,x-y,-z+1/2"
	symLines[471]  = "x,y,z;x-y,x,z;-y,x-y,z;-x,-y,z;-x+y,-x,z;y,-x+y,z;x-y,-y,-z;-x,-x+y,-z;y,x,-z;-y,-x,-z;-x+y,y,-z;x,x-y,-z"
	symLines[472]  = "x,y,z;x-y,x,z+1/6;-y,x-y,z+1/3;-x,-y,z+1/2;-x+y,-x,z-1/3;y,-x+y,z-1/6;x-y,-y,-z;-x,-x+y,-z-1/3;y,x,-z+1/3;-y,-x,-z-1/6;"
	symLines[472] += "-x+y,y,-z+1/2;x,x-y,-z+1/6"
	symLines[473]  = "x,y,z;x-y,x,z-1/6;-y,x-y,z-1/3;-x,-y,z+1/2;-x+y,-x,z+1/3;y,-x+y,z+1/6;x-y,-y,-z;-x,-x+y,-z+1/3;y,x,-z-1/3;-y,-x,-z+1/6;"
	symLines[473] += "-x+y,y,-z+1/2;x,x-y,-z-1/6"
	symLines[474]  = "x,y,z;x-y,x,z+1/3;-y,x-y,z-1/3;-x,-y,z;-x+y,-x,z+1/3;y,-x+y,z-1/3;x-y,-y,-z;-x,-x+y,-z+1/3;y,x,-z-1/3;-y,-x,-z-1/3;"
	symLines[474] += "-x+y,y,-z;x,x-y,-z+1/3"
	symLines[475]  = "x,y,z;x-y,x,z-1/3;-y,x-y,z+1/3;-x,-y,z;-x+y,-x,z-1/3;y,-x+y,z+1/3;x-y,-y,-z;-x,-x+y,-z-1/3;y,x,-z+1/3;-y,-x,-z+1/3;"
	symLines[475] += "-x+y,y,-z;x,x-y,-z-1/3"
	symLines[476]  = "x,y,z;x-y,x,z+1/2;-y,x-y,z;-x,-y,z+1/2;-x+y,-x,z;y,-x+y,z+1/2;x-y,-y,-z;-x,-x+y,-z;y,x,-z;-y,-x,-z+1/2;-x+y,y,-z+1/2;x,x-y,-z+1/2"
	symLines[477]  = "x,y,z;x-y,x,z;-y,x-y,z;-x,-y,z;-x+y,-x,z;y,-x+y,z;-x+y,y,z;x,x-y,z;-y,-x,z;y,x,z;x-y,-y,z;-x,-x+y,z"
	symLines[478]  = "x,y,z;x-y,x,z;-y,x-y,z;-x,-y,z;-x+y,-x,z;y,-x+y,z;-x+y,y,z+1/2;x,x-y,z+1/2;-y,-x,z+1/2;y,x,z+1/2;x-y,-y,z+1/2;-x,-x+y,z+1/2"
	symLines[479]  = "x,y,z;x-y,x,z+1/2;-y,x-y,z;-x,-y,z+1/2;-x+y,-x,z;y,-x+y,z+1/2;-x+y,y,z+1/2;x,x-y,z+1/2;-y,-x,z+1/2;y,x,z;x-y,-y,z;-x,-x+y,z"
	symLines[480]  = "x,y,z;x-y,x,z+1/2;-y,x-y,z;-x,-y,z+1/2;-x+y,-x,z;y,-x+y,z+1/2;-x+y,y,z;x,x-y,z;-y,-x,z;y,x,z+1/2;x-y,-y,z+1/2;-x,-x+y,z+1/2"
	symLines[481]  = "x,y,z;-x+y,-x,-z;-y,x-y,z;x,y,-z;-x+y,-x,z;-y,x-y,-z;-y,-x,-z;-x+y,y,-z;x,x-y,-z;-x+y,y,z;x,x-y,z;-y,-x,z"
	symLines[482]  = "x,y,z;-x+y,-x,-z+1/2;-y,x-y,z;x,y,-z+1/2;-x+y,-x,z;-y,x-y,-z+1/2;-y,-x,-z;-x+y,y,-z;x,x-y,-z;-x+y,y,z+1/2;x,x-y,z+1/2;-y,-x,z+1/2"
	symLines[483]  = "x,y,z;-x+y,-x,-z;-y,x-y,z;x,y,-z;-x+y,-x,z;-y,x-y,-z;x-y,-y,-z;-x,-x+y,-z;y,x,-z;y,x,z;x-y,-y,z;-x,-x+y,z"
	symLines[484]  = "x,y,z;-x+y,-x,-z+1/2;-y,x-y,z;x,y,-z+1/2;-x+y,-x,z;-y,x-y,-z+1/2;x-y,-y,-z;-x,-x+y,-z;y,x,-z;y,x,z+1/2;x-y,-y,z+1/2;-x,-x+y,z+1/2"
	symLines[485]  = "x,y,z;x-y,x,z;-y,x-y,z;-x,-y,z;-x+y,-x,z;y,-x+y,z;x-y,-y,-z;-x,-x+y,-z;y,x,-z;-y,-x,-z;-x+y,y,-z;x,x-y,-z;-x,-y,-z;"
	symLines[485] += "-x+y,-x,-z;y,-x+y,-z;x,y,-z;x-y,x,-z;-y,x-y,-z;-x+y,y,z;x,x-y,z;-y,-x,z;y,x,z;x-y,-y,z;-x,-x+y,z"
	symLines[486]  = "x,y,z;x-y,x,z;-y,x-y,z;-x,-y,z;-x+y,-x,z;y,-x+y,z;x-y,-y,-z+1/2;-x,-x+y,-z+1/2;y,x,-z+1/2;-y,-x,-z+1/2;-x+y,y,-z+1/2;"
	symLines[486] += "x,x-y,-z+1/2;-x,-y,-z;-x+y,-x,-z;y,-x+y,-z;x,y,-z;x-y,x,-z;-y,x-y,-z;-x+y,y,z+1/2;x,x-y,z+1/2;-y,-x,z+1/2;y,x,z+1/2;"
	symLines[486] += "x-y,-y,z+1/2;-x,-x+y,z+1/2"
	symLines[487]  = "x,y,z;x-y,x,z+1/2;-y,x-y,z;-x,-y,z+1/2;-x+y,-x,z;y,-x+y,z+1/2;x-y,-y,-z+1/2;-x,-x+y,-z+1/2;y,x,-z+1/2;-y,-x,-z;-x+y,y,-z;"
	symLines[487] += "x,x-y,-z;-x,-y,-z;-x+y,-x,-z+1/2;y,-x+y,-z;x,y,-z+1/2;x-y,x,-z;-y,x-y,-z+1/2;-x+y,y,z+1/2;x,x-y,z+1/2;-y,-x,z+1/2;y,x,z;"
	symLines[487] += "x-y,-y,z;-x,-x+y,z"
	symLines[488]  = "x,y,z;x-y,x,z+1/2;-y,x-y,z;-x,-y,z+1/2;-x+y,-x,z;y,-x+y,z+1/2;x-y,-y,-z;-x,-x+y,-z;y,x,-z;-y,-x,-z+1/2;-x+y,y,-z+1/2;"
	symLines[488] += "x,x-y,-z+1/2;-x,-y,-z;-x+y,-x,-z+1/2;y,-x+y,-z;x,y,-z+1/2;x-y,x,-z;-y,x-y,-z+1/2;-x+y,y,z;x,x-y,z;-y,-x,z;y,x,z+1/2;"
	symLines[488] += "x-y,-y,z+1/2;-x,-x+y,z+1/2"
	// Cubic SG[195,230]  SG_idNum [489,530]   (42 idNums)
	symLines[489]  = "x,y,z;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;-x,-y,z;x,-y,-z;-x,y,-z"
	symLines[490]  = "x,y,z;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;-x,-y,z;x,-y,-z;-x,y,-z;x,y+1/2,z+1/2;z,x+1/2,y+1/2;"
	symLines[490] += "y,z+1/2,x+1/2;-y,-z+1/2,x+1/2;z,-x+1/2,-y+1/2;-y,z+1/2,-x+1/2;-z,-x+1/2,y+1/2;-z,x+1/2,-y+1/2;y,-z+1/2,-x+1/2;"
	symLines[490] += "-x,-y+1/2,z+1/2;x,-y+1/2,-z+1/2;-x,y+1/2,-z+1/2;x+1/2,y,z+1/2;z+1/2,x,y+1/2;y+1/2,z,x+1/2;-y+1/2,-z,x+1/2;"
	symLines[490] += "z+1/2,-x,-y+1/2;-y+1/2,z,-x+1/2;-z+1/2,-x,y+1/2;-z+1/2,x,-y+1/2;y+1/2,-z,-x+1/2;-x+1/2,-y,z+1/2;x+1/2,-y,-z+1/2;"
	symLines[490] += "-x+1/2,y,-z+1/2;x+1/2,y+1/2,z;z+1/2,x+1/2,y;y+1/2,z+1/2,x;-y+1/2,-z+1/2,x;z+1/2,-x+1/2,-y;-y+1/2,z+1/2,-x;"
	symLines[490] += "-z+1/2,-x+1/2,y;-z+1/2,x+1/2,-y;y+1/2,-z+1/2,-x;-x+1/2,-y+1/2,z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z"
	symLines[491]  = "x,y,z;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;-x,-y,z;x,-y,-z;-x,y,-z;x+1/2,y+1/2,z+1/2;"
	symLines[491] += "z+1/2,x+1/2,y+1/2;y+1/2,z+1/2,x+1/2;-y+1/2,-z+1/2,x+1/2;z+1/2,-x+1/2,-y+1/2;-y+1/2,z+1/2,-x+1/2;-z+1/2,-x+1/2,y+1/2;"
	symLines[491] += "-z+1/2,x+1/2,-y+1/2;y+1/2,-z+1/2,-x+1/2;-x+1/2,-y+1/2,z+1/2;x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2"
	symLines[492]  = "x,y,z;z,x,y;y,z,x;-y+1/2,-z,x+1/2;z+1/2,-x+1/2,-y;-y,z+1/2,-x+1/2;-z+1/2,-x,y+1/2;-z,x+1/2,-y+1/2;y+1/2,-z+1/2,-x;"
	symLines[492] += "-x+1/2,-y,z+1/2;x+1/2,-y+1/2,-z;-x,y+1/2,-z+1/2"
	symLines[493]  = "x,y,z;z,x,y;y,z,x;-y,-z+1/2,x;z,-x,-y+1/2;-y+1/2,z,-x;-z,-x+1/2,y;-z+1/2,x,-y;y,-z,-x+1/2;-x,-y+1/2,z;x,-y,-z+1/2;"
	symLines[493] += "-x+1/2,y,-z;x+1/2,y+1/2,z+1/2;z+1/2,x+1/2,y+1/2;y+1/2,z+1/2,x+1/2;-y+1/2,-z,x+1/2;z+1/2,-x+1/2,-y;-y,z+1/2,-x+1/2;"
	symLines[493] += "-z+1/2,-x,y+1/2;-z,x+1/2,-y+1/2;y+1/2,-z+1/2,-x;-x+1/2,-y,z+1/2;x+1/2,-y+1/2,-z;-x,y+1/2,-z+1/2"
	symLines[494]  = "x,y,z;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;-x,-y,z;x,-y,-z;-x,y,-z;-x,-y,-z;-z,-x,-y;-y,-z,-x;"
	symLines[494] += "y,z,-x;-z,x,y;y,-z,x;z,x,-y;z,-x,y;-y,z,x;x,y,-z;-x,y,z;x,-y,z"
	symLines[495]  = "x,y,z;-x+1/2,-y+1/2,-z+1/2;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;-z+1/2,-x+1/2,-y+1/2;"
	symLines[495] += "-y+1/2,-z+1/2,-x+1/2;y+1/2,z+1/2,-x+1/2;-z+1/2,x+1/2,y+1/2;y+1/2,-z+1/2,x+1/2;z+1/2,x+1/2,-y+1/2;z+1/2,-x+1/2,y+1/2;"
	symLines[495] += "-y+1/2,z+1/2,x+1/2;-x,-y,z;x,-y,-z;-x,y,-z;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[496]  = "x,y,z;z,x,y;y,z,x;-y+1/2,-z+1/2,x;z,-x+1/2,-y+1/2;-y+1/2,z,-x+1/2;-z+1/2,-x+1/2,y;-z+1/2,x,-y+1/2;y,-z+1/2,-x+1/2;"
	symLines[496] += "-x+1/2,-y+1/2,z;x,-y+1/2,-z+1/2;-x+1/2,y,-z+1/2;-x,-y,-z;-z,-x,-y;-y,-z,-x;y+1/2,z+1/2,-x;-z,x+1/2,y+1/2;y+1/2,-z,x+1/2;"
	symLines[496] += "z+1/2,x+1/2,-y;z+1/2,-x,y+1/2;-y,z+1/2,x+1/2;x+1/2,y+1/2,-z;-x,y+1/2,z+1/2;x+1/2,-y,z+1/2"
	symLines[497]  = "x,y,z;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;-x,-y,z;x,-y,-z;-x,y,-z;-x,-y,-z;-z,-x,-y;-y,-z,-x;"
	symLines[497] += "y,z,-x;-z,x,y;y,-z,x;z,x,-y;z,-x,y;-y,z,x;x,y,-z;-x,y,z;x,-y,z;x,y+1/2,z+1/2;z,x+1/2,y+1/2;y,z+1/2,x+1/2;-y,-z+1/2,x+1/2;"
	symLines[497] += "z,-x+1/2,-y+1/2;-y,z+1/2,-x+1/2;-z,-x+1/2,y+1/2;-z,x+1/2,-y+1/2;y,-z+1/2,-x+1/2;-x,-y+1/2,z+1/2;x,-y+1/2,-z+1/2;"
	symLines[497] += "-x,y+1/2,-z+1/2;-x,-y+1/2,-z+1/2;-z,-x+1/2,-y+1/2;-y,-z+1/2,-x+1/2;y,z+1/2,-x+1/2;-z,x+1/2,y+1/2;y,-z+1/2,x+1/2;"
	symLines[497] += "z,x+1/2,-y+1/2;z,-x+1/2,y+1/2;-y,z+1/2,x+1/2;x,y+1/2,-z+1/2;-x,y+1/2,z+1/2;x,-y+1/2,z+1/2;x+1/2,y,z+1/2;z+1/2,x,y+1/2;"
	symLines[497] += "y+1/2,z,x+1/2;-y+1/2,-z,x+1/2;z+1/2,-x,-y+1/2;-y+1/2,z,-x+1/2;-z+1/2,-x,y+1/2;-z+1/2,x,-y+1/2;y+1/2,-z,-x+1/2;"
	symLines[497] += "-x+1/2,-y,z+1/2;x+1/2,-y,-z+1/2;-x+1/2,y,-z+1/2;-x+1/2,-y,-z+1/2;-z+1/2,-x,-y+1/2;-y+1/2,-z,-x+1/2;y+1/2,z,-x+1/2;"
	symLines[497] += "-z+1/2,x,y+1/2;y+1/2,-z,x+1/2;z+1/2,x,-y+1/2;z+1/2,-x,y+1/2;-y+1/2,z,x+1/2;x+1/2,y,-z+1/2;-x+1/2,y,z+1/2;x+1/2,-y,z+1/2;"
	symLines[497] += "x+1/2,y+1/2,z;z+1/2,x+1/2,y;y+1/2,z+1/2,x;-y+1/2,-z+1/2,x;z+1/2,-x+1/2,-y;-y+1/2,z+1/2,-x;-z+1/2,-x+1/2,y;"
	symLines[497] += "-z+1/2,x+1/2,-y;y+1/2,-z+1/2,-x;-x+1/2,-y+1/2,z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;-x+1/2,-y+1/2,-z;-z+1/2,-x+1/2,-y;"
	symLines[497] += "-y+1/2,-z+1/2,-x;y+1/2,z+1/2,-x;-z+1/2,x+1/2,y;y+1/2,-z+1/2,x;z+1/2,x+1/2,-y;z+1/2,-x+1/2,y;-y+1/2,z+1/2,x;"
	symLines[497] += "x+1/2,y+1/2,-z;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[498]  = "x,y,z;-x+1/4,-y+1/4,-z+1/4;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;-z+1/4,-x+1/4,-y+1/4;"
	symLines[498] += "-y+1/4,-z+1/4,-x+1/4;y+1/4,z+1/4,-x+1/4;-z+1/4,x+1/4,y+1/4;y+1/4,-z+1/4,x+1/4;z+1/4,x+1/4,-y+1/4;z+1/4,-x+1/4,y+1/4;"
	symLines[498] += "-y+1/4,z+1/4,x+1/4;-x,-y,z;x,-y,-z;-x,y,-z;x+1/4,y+1/4,-z+1/4;-x+1/4,y+1/4,z+1/4;x+1/4,-y+1/4,z+1/4;x,y+1/2,z+1/2;"
	symLines[498] += "-x+1/4,-y-1/4,-z-1/4;z,x+1/2,y+1/2;y,z+1/2,x+1/2;-y,-z+1/2,x+1/2;z,-x+1/2,-y+1/2;-y,z+1/2,-x+1/2;-z,-x+1/2,y+1/2;"
	symLines[498] += "-z,x+1/2,-y+1/2;y,-z+1/2,-x+1/2;-z+1/4,-x-1/4,-y-1/4;-y+1/4,-z-1/4,-x-1/4;y+1/4,z-1/4,-x-1/4;-z+1/4,x-1/4,y-1/4;"
	symLines[498] += "y+1/4,-z-1/4,x-1/4;z+1/4,x-1/4,-y-1/4;z+1/4,-x-1/4,y-1/4;-y+1/4,z-1/4,x-1/4;-x,-y+1/2,z+1/2;x,-y+1/2,-z+1/2;"
	symLines[498] += "-x,y+1/2,-z+1/2;x+1/4,y-1/4,-z-1/4;-x+1/4,y-1/4,z-1/4;x+1/4,-y-1/4,z-1/4;x+1/2,y,z+1/2;-x-1/4,-y+1/4,-z-1/4;"
	symLines[498] += "z+1/2,x,y+1/2;y+1/2,z,x+1/2;-y+1/2,-z,x+1/2;z+1/2,-x,-y+1/2;-y+1/2,z,-x+1/2;-z+1/2,-x,y+1/2;-z+1/2,x,-y+1/2;"
	symLines[498] += "y+1/2,-z,-x+1/2;-z-1/4,-x+1/4,-y-1/4;-y-1/4,-z+1/4,-x-1/4;y-1/4,z+1/4,-x-1/4;-z-1/4,x+1/4,y-1/4;y-1/4,-z+1/4,x-1/4;"
	symLines[498] += "z-1/4,x+1/4,-y-1/4;z-1/4,-x+1/4,y-1/4;-y-1/4,z+1/4,x-1/4;-x+1/2,-y,z+1/2;x+1/2,-y,-z+1/2;-x+1/2,y,-z+1/2;"
	symLines[498] += "x-1/4,y+1/4,-z-1/4;-x-1/4,y+1/4,z-1/4;x-1/4,-y+1/4,z-1/4;x+1/2,y+1/2,z;-x-1/4,-y-1/4,-z+1/4;z+1/2,x+1/2,y;y+1/2,z+1/2,x;"
	symLines[498] += "-y+1/2,-z+1/2,x;z+1/2,-x+1/2,-y;-y+1/2,z+1/2,-x;-z+1/2,-x+1/2,y;-z+1/2,x+1/2,-y;y+1/2,-z+1/2,-x;-z-1/4,-x-1/4,-y+1/4;"
	symLines[498] += "-y-1/4,-z-1/4,-x+1/4;y-1/4,z-1/4,-x+1/4;-z-1/4,x-1/4,y+1/4;y-1/4,-z-1/4,x+1/4;z-1/4,x-1/4,-y+1/4;z-1/4,-x-1/4,y+1/4;"
	symLines[498] += "-y-1/4,z-1/4,x+1/4;-x+1/2,-y+1/2,z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;x-1/4,y-1/4,-z+1/4;-x-1/4,y-1/4,z+1/4;x-1/4,-y-1/4,z+1/4"
	symLines[499]  = "x,y,z;z,x,y;y,z,x;-y+1/4,-z+1/4,x;z,-x+1/4,-y+1/4;-y+1/4,z,-x+1/4;-z+1/4,-x+1/4,y;-z+1/4,x,-y+1/4;y,-z+1/4,-x+1/4;"
	symLines[499] += "-x+1/4,-y+1/4,z;x,-y+1/4,-z+1/4;-x+1/4,y,-z+1/4;-x,-y,-z;-z,-x,-y;-y,-z,-x;y-1/4,z-1/4,-x;-z,x-1/4,y-1/4;y-1/4,-z,x-1/4;"
	symLines[499] += "z-1/4,x-1/4,-y;z-1/4,-x,y-1/4;-y,z-1/4,x-1/4;x-1/4,y-1/4,-z;-x,y-1/4,z-1/4;x-1/4,-y,z-1/4;x,y+1/2,z+1/2;z,x+1/2,y+1/2;"
	symLines[499] += "y,z+1/2,x+1/2;-y+1/4,-z-1/4,x+1/2;z,-x-1/4,-y-1/4;-y+1/4,z+1/2,-x-1/4;-z+1/4,-x-1/4,y+1/2;-z+1/4,x+1/2,-y-1/4;"
	symLines[499] += "y,-z-1/4,-x-1/4;-x+1/4,-y-1/4,z+1/2;x,-y-1/4,-z-1/4;-x+1/4,y+1/2,-z-1/4;-x,-y+1/2,-z+1/2;-z,-x+1/2,-y+1/2;"
	symLines[499] += "-y,-z+1/2,-x+1/2;y-1/4,z+1/4,-x+1/2;-z,x+1/4,y+1/4;y-1/4,-z+1/2,x+1/4;z-1/4,x+1/4,-y+1/2;z-1/4,-x+1/2,y+1/4;"
	symLines[499] += "-y,z+1/4,x+1/4;x-1/4,y+1/4,-z+1/2;-x,y+1/4,z+1/4;x-1/4,-y+1/2,z+1/4;x+1/2,y,z+1/2;z+1/2,x,y+1/2;y+1/2,z,x+1/2;"
	symLines[499] += "-y-1/4,-z+1/4,x+1/2;z+1/2,-x+1/4,-y-1/4;-y-1/4,z,-x-1/4;-z-1/4,-x+1/4,y+1/2;-z-1/4,x,-y-1/4;y+1/2,-z+1/4,-x-1/4;"
	symLines[499] += "-x-1/4,-y+1/4,z+1/2;x+1/2,-y+1/4,-z-1/4;-x-1/4,y,-z-1/4;-x+1/2,-y,-z+1/2;-z+1/2,-x,-y+1/2;-y+1/2,-z,-x+1/2;"
	symLines[499] += "y+1/4,z-1/4,-x+1/2;-z+1/2,x-1/4,y+1/4;y+1/4,-z,x+1/4;z+1/4,x-1/4,-y+1/2;z+1/4,-x,y+1/4;-y+1/2,z-1/4,x+1/4;"
	symLines[499] += "x+1/4,y-1/4,-z+1/2;-x+1/2,y-1/4,z+1/4;x+1/4,-y,z+1/4;x+1/2,y+1/2,z;z+1/2,x+1/2,y;y+1/2,z+1/2,x;-y-1/4,-z-1/4,x;"
	symLines[499] += "z+1/2,-x-1/4,-y+1/4;-y-1/4,z+1/2,-x+1/4;-z-1/4,-x-1/4,y;-z-1/4,x+1/2,-y+1/4;y+1/2,-z-1/4,-x+1/4;-x-1/4,-y-1/4,z;"
	symLines[499] += "x+1/2,-y-1/4,-z+1/4;-x-1/4,y+1/2,-z+1/4;-x+1/2,-y+1/2,-z;-z+1/2,-x+1/2,-y;-y+1/2,-z+1/2,-x;y+1/4,z+1/4,-x;"
	symLines[499] += "-z+1/2,x+1/4,y-1/4;y+1/4,-z+1/2,x-1/4;z+1/4,x+1/4,-y;z+1/4,-x+1/2,y-1/4;-y+1/2,z+1/4,x-1/4;x+1/4,y+1/4,-z;"
	symLines[499] += "-x+1/2,y+1/4,z-1/4;x+1/4,-y+1/2,z-1/4"
	symLines[500]  = "x,y,z;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;-x,-y,z;x,-y,-z;-x,y,-z;-x,-y,-z;-z,-x,-y;-y,-z,-x;"
	symLines[500] += "y,z,-x;-z,x,y;y,-z,x;z,x,-y;z,-x,y;-y,z,x;x,y,-z;-x,y,z;x,-y,z;x+1/2,y+1/2,z+1/2;z+1/2,x+1/2,y+1/2;y+1/2,z+1/2,x+1/2;"
	symLines[500] += "-y+1/2,-z+1/2,x+1/2;z+1/2,-x+1/2,-y+1/2;-y+1/2,z+1/2,-x+1/2;-z+1/2,-x+1/2,y+1/2;-z+1/2,x+1/2,-y+1/2;y+1/2,-z+1/2,-x+1/2;"
	symLines[500] += "-x+1/2,-y+1/2,z+1/2;x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2;-x+1/2,-y+1/2,-z+1/2;-z+1/2,-x+1/2,-y+1/2;"
	symLines[500] += "-y+1/2,-z+1/2,-x+1/2;y+1/2,z+1/2,-x+1/2;-z+1/2,x+1/2,y+1/2;y+1/2,-z+1/2,x+1/2;z+1/2,x+1/2,-y+1/2;z+1/2,-x+1/2,y+1/2;"
	symLines[500] += "-y+1/2,z+1/2,x+1/2;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[501]  = "x,y,z;z,x,y;y,z,x;-y+1/2,-z,x+1/2;z+1/2,-x+1/2,-y;-y,z+1/2,-x+1/2;-z+1/2,-x,y+1/2;-z,x+1/2,-y+1/2;y+1/2,-z+1/2,-x;"
	symLines[501] += "-x+1/2,-y,z+1/2;x+1/2,-y+1/2,-z;-x,y+1/2,-z+1/2;-x,-y,-z;-z,-x,-y;-y,-z,-x;y+1/2,z,-x+1/2;-z+1/2,x+1/2,y;y,-z+1/2,x+1/2;"
	symLines[501] += "z+1/2,x,-y+1/2;z,-x+1/2,y+1/2;-y+1/2,z+1/2,x;x+1/2,y,-z+1/2;-x+1/2,y+1/2,z;x,-y+1/2,z+1/2"
	symLines[502]  = "x,y,z;z,x,y;y,z,x;-y,-z+1/2,x;z,-x,-y+1/2;-y+1/2,z,-x;-z,-x+1/2,y;-z+1/2,x,-y;y,-z,-x+1/2;-x,-y+1/2,z;x,-y,-z+1/2;"
	symLines[502] += "-x+1/2,y,-z;-x,-y,-z;-z,-x,-y;-y,-z,-x;y,z+1/2,-x;-z,x,y+1/2;y+1/2,-z,x;z,x+1/2,-y;z+1/2,-x,y;-y,z,x+1/2;x,y+1/2,-z;"
	symLines[502] += "-x,y,z+1/2;x+1/2,-y,z;x+1/2,y+1/2,z+1/2;z+1/2,x+1/2,y+1/2;y+1/2,z+1/2,x+1/2;-y+1/2,-z,x+1/2;z+1/2,-x+1/2,-y;"
	symLines[502] += "-y,z+1/2,-x+1/2;-z+1/2,-x,y+1/2;-z,x+1/2,-y+1/2;y+1/2,-z+1/2,-x;-x+1/2,-y,z+1/2;x+1/2,-y+1/2,-z;-x,y+1/2,-z+1/2;"
	symLines[502] += "-x+1/2,-y+1/2,-z+1/2;-z+1/2,-x+1/2,-y+1/2;-y+1/2,-z+1/2,-x+1/2;y+1/2,z,-x+1/2;-z+1/2,x+1/2,y;y,-z+1/2,x+1/2;"
	symLines[502] += "z+1/2,x,-y+1/2;z,-x+1/2,y+1/2;-y+1/2,z+1/2,x;x+1/2,y,-z+1/2;-x+1/2,y+1/2,z;x,-y+1/2,z+1/2"
	symLines[503]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x,-z,y;x,-y,-z;x,z,-y;z,y,-x;-x,y,-z;-z,y,x;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;"
	symLines[503] += "-z,x,-y;y,-z,-x;y,x,-z;-y,-x,-z;-x,z,y;-x,-z,-y;z,-y,x;-z,-y,-x"
	symLines[504]  = "x,y,z;-y+1/2,x+1/2,z+1/2;-x,-y,z;y+1/2,-x+1/2,z+1/2;x+1/2,-z+1/2,y+1/2;x,-y,-z;x+1/2,z+1/2,-y+1/2;z+1/2,y+1/2,-x+1/2;"
	symLines[504] += "-x,y,-z;-z+1/2,y+1/2,x+1/2;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;y+1/2,x+1/2,-z+1/2;"
	symLines[504] += "-y+1/2,-x+1/2,-z+1/2;-x+1/2,z+1/2,y+1/2;-x+1/2,-z+1/2,-y+1/2;z+1/2,-y+1/2,x+1/2;-z+1/2,-y+1/2,-x+1/2"
	symLines[505]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x,-z,y;x,-y,-z;x,z,-y;z,y,-x;-x,y,-z;-z,y,x;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;"
	symLines[505] += "-z,x,-y;y,-z,-x;y,x,-z;-y,-x,-z;-x,z,y;-x,-z,-y;z,-y,x;-z,-y,-x;x,y+1/2,z+1/2;-y,x+1/2,z+1/2;-x,-y+1/2,z+1/2;"
	symLines[505] += "y,-x+1/2,z+1/2;x,-z+1/2,y+1/2;x,-y+1/2,-z+1/2;x,z+1/2,-y+1/2;z,y+1/2,-x+1/2;-x,y+1/2,-z+1/2;-z,y+1/2,x+1/2;z,x+1/2,y+1/2;"
	symLines[505] += "y,z+1/2,x+1/2;-y,-z+1/2,x+1/2;z,-x+1/2,-y+1/2;-y,z+1/2,-x+1/2;-z,-x+1/2,y+1/2;-z,x+1/2,-y+1/2;y,-z+1/2,-x+1/2;"
	symLines[505] += "y,x+1/2,-z+1/2;-y,-x+1/2,-z+1/2;-x,z+1/2,y+1/2;-x,-z+1/2,-y+1/2;z,-y+1/2,x+1/2;-z,-y+1/2,-x+1/2;x+1/2,y,z+1/2;"
	symLines[505] += "-y+1/2,x,z+1/2;-x+1/2,-y,z+1/2;y+1/2,-x,z+1/2;x+1/2,-z,y+1/2;x+1/2,-y,-z+1/2;x+1/2,z,-y+1/2;z+1/2,y,-x+1/2;"
	symLines[505] += "-x+1/2,y,-z+1/2;-z+1/2,y,x+1/2;z+1/2,x,y+1/2;y+1/2,z,x+1/2;-y+1/2,-z,x+1/2;z+1/2,-x,-y+1/2;-y+1/2,z,-x+1/2;"
	symLines[505] += "-z+1/2,-x,y+1/2;-z+1/2,x,-y+1/2;y+1/2,-z,-x+1/2;y+1/2,x,-z+1/2;-y+1/2,-x,-z+1/2;-x+1/2,z,y+1/2;-x+1/2,-z,-y+1/2;"
	symLines[505] += "z+1/2,-y,x+1/2;-z+1/2,-y,-x+1/2;x+1/2,y+1/2,z;-y+1/2,x+1/2,z;-x+1/2,-y+1/2,z;y+1/2,-x+1/2,z;x+1/2,-z+1/2,y;"
	symLines[505] += "x+1/2,-y+1/2,-z;x+1/2,z+1/2,-y;z+1/2,y+1/2,-x;-x+1/2,y+1/2,-z;-z+1/2,y+1/2,x;z+1/2,x+1/2,y;y+1/2,z+1/2,x;-y+1/2,-z+1/2,x;"
	symLines[505] += "z+1/2,-x+1/2,-y;-y+1/2,z+1/2,-x;-z+1/2,-x+1/2,y;-z+1/2,x+1/2,-y;y+1/2,-z+1/2,-x;y+1/2,x+1/2,-z;-y+1/2,-x+1/2,-z;"
	symLines[505] += "-x+1/2,z+1/2,y;-x+1/2,-z+1/2,-y;z+1/2,-y+1/2,x;-z+1/2,-y+1/2,-x"
	symLines[506]  = "x,y,z;-y+1/4,x+1/4,z+1/4;-x,-y,z;y+1/4,-x+1/4,z+1/4;x+1/4,-z+1/4,y+1/4;x,-y,-z;x+1/4,z+1/4,-y+1/4;z+1/4,y+1/4,-x+1/4;"
	symLines[506] += "-x,y,-z;-z+1/4,y+1/4,x+1/4;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;y+1/4,x+1/4,-z+1/4;"
	symLines[506] += "-y+1/4,-x+1/4,-z+1/4;-x+1/4,z+1/4,y+1/4;-x+1/4,-z+1/4,-y+1/4;z+1/4,-y+1/4,x+1/4;-z+1/4,-y+1/4,-x+1/4;x,y+1/2,z+1/2;"
	symLines[506] += "-y+1/4,x-1/4,z-1/4;-x,-y+1/2,z+1/2;y+1/4,-x-1/4,z-1/4;x+1/4,-z-1/4,y-1/4;x,-y+1/2,-z+1/2;x+1/4,z-1/4,-y-1/4;"
	symLines[506] += "z+1/4,y-1/4,-x-1/4;-x,y+1/2,-z+1/2;-z+1/4,y-1/4,x-1/4;z,x+1/2,y+1/2;y,z+1/2,x+1/2;-y,-z+1/2,x+1/2;z,-x+1/2,-y+1/2;"
	symLines[506] += "-y,z+1/2,-x+1/2;-z,-x+1/2,y+1/2;-z,x+1/2,-y+1/2;y,-z+1/2,-x+1/2;y+1/4,x-1/4,-z-1/4;-y+1/4,-x-1/4,-z-1/4;"
	symLines[506] += "-x+1/4,z-1/4,y-1/4;-x+1/4,-z-1/4,-y-1/4;z+1/4,-y-1/4,x-1/4;-z+1/4,-y-1/4,-x-1/4;x+1/2,y,z+1/2;-y-1/4,x+1/4,z-1/4;"
	symLines[506] += "-x+1/2,-y,z+1/2;y-1/4,-x+1/4,z-1/4;x-1/4,-z+1/4,y-1/4;x+1/2,-y,-z+1/2;x-1/4,z+1/4,-y-1/4;z-1/4,y+1/4,-x-1/4;"
	symLines[506] += "-x+1/2,y,-z+1/2;-z-1/4,y+1/4,x-1/4;z+1/2,x,y+1/2;y+1/2,z,x+1/2;-y+1/2,-z,x+1/2;z+1/2,-x,-y+1/2;-y+1/2,z,-x+1/2;"
	symLines[506] += "-z+1/2,-x,y+1/2;-z+1/2,x,-y+1/2;y+1/2,-z,-x+1/2;y-1/4,x+1/4,-z-1/4;-y-1/4,-x+1/4,-z-1/4;-x-1/4,z+1/4,y-1/4;"
	symLines[506] += "-x-1/4,-z+1/4,-y-1/4;z-1/4,-y+1/4,x-1/4;-z-1/4,-y+1/4,-x-1/4;x+1/2,y+1/2,z;-y-1/4,x-1/4,z+1/4;-x+1/2,-y+1/2,z;"
	symLines[506] += "y-1/4,-x-1/4,z+1/4;x-1/4,-z-1/4,y+1/4;x+1/2,-y+1/2,-z;x-1/4,z-1/4,-y+1/4;z-1/4,y-1/4,-x+1/4;-x+1/2,y+1/2,-z;"
	symLines[506] += "-z-1/4,y-1/4,x+1/4;z+1/2,x+1/2,y;y+1/2,z+1/2,x;-y+1/2,-z+1/2,x;z+1/2,-x+1/2,-y;-y+1/2,z+1/2,-x;-z+1/2,-x+1/2,y;"
	symLines[506] += "-z+1/2,x+1/2,-y;y+1/2,-z+1/2,-x;y-1/4,x-1/4,-z+1/4;-y-1/4,-x-1/4,-z+1/4;-x-1/4,z-1/4,y+1/4;-x-1/4,-z-1/4,-y+1/4;"
	symLines[506] += "z-1/4,-y-1/4,x+1/4;-z-1/4,-y-1/4,-x+1/4"
	symLines[507]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x,-z,y;x,-y,-z;x,z,-y;z,y,-x;-x,y,-z;-z,y,x;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;"
	symLines[507] += "-z,x,-y;y,-z,-x;y,x,-z;-y,-x,-z;-x,z,y;-x,-z,-y;z,-y,x;-z,-y,-x;x+1/2,y+1/2,z+1/2;-y+1/2,x+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;"
	symLines[507] += "y+1/2,-x+1/2,z+1/2;x+1/2,-z+1/2,y+1/2;x+1/2,-y+1/2,-z+1/2;x+1/2,z+1/2,-y+1/2;z+1/2,y+1/2,-x+1/2;-x+1/2,y+1/2,-z+1/2;"
	symLines[507] += "-z+1/2,y+1/2,x+1/2;z+1/2,x+1/2,y+1/2;y+1/2,z+1/2,x+1/2;-y+1/2,-z+1/2,x+1/2;z+1/2,-x+1/2,-y+1/2;-y+1/2,z+1/2,-x+1/2;"
	symLines[507] += "-z+1/2,-x+1/2,y+1/2;-z+1/2,x+1/2,-y+1/2;y+1/2,-z+1/2,-x+1/2;y+1/2,x+1/2,-z+1/2;-y+1/2,-x+1/2,-z+1/2;-x+1/2,z+1/2,y+1/2;"
	symLines[507] += "-x+1/2,-z+1/2,-y+1/2;z+1/2,-y+1/2,x+1/2;-z+1/2,-y+1/2,-x+1/2"
	symLines[508]  = "x,y,z;-y-1/4,x+1/4,z-1/4;-x+1/2,-y,z+1/2;y-1/4,-x-1/4,z+1/4;x-1/4,-z-1/4,y+1/4;x+1/2,-y+1/2,-z;x+1/4,z-1/4,-y-1/4;"
	symLines[508] += "z+1/4,y-1/4,-x-1/4;-x,y+1/2,-z+1/2;-z-1/4,y+1/4,x-1/4;z,x,y;y,z,x;-y+1/2,-z,x+1/2;z+1/2,-x+1/2,-y;-y,z+1/2,-x+1/2;"
	symLines[508] += "-z+1/2,-x,y+1/2;-z,x+1/2,-y+1/2;y+1/2,-z+1/2,-x;y+1/4,x-1/4,-z-1/4;-y+1/4,-x+1/4,-z+1/4;-x-1/4,z+1/4,y-1/4;"
	symLines[508] += "-x+1/4,-z+1/4,-y+1/4;z-1/4,-y-1/4,x+1/4;-z+1/4,-y+1/4,-x+1/4"
	symLines[509]  = "x,y,z;-y+1/4,x-1/4,z+1/4;-x+1/2,-y,z+1/2;y+1/4,-x+1/4,z-1/4;x+1/4,-z+1/4,y-1/4;x+1/2,-y+1/2,-z;x-1/4,z+1/4,-y+1/4;"
	symLines[509] += "z-1/4,y+1/4,-x+1/4;-x,y+1/2,-z+1/2;-z+1/4,y-1/4,x+1/4;z,x,y;y,z,x;-y+1/2,-z,x+1/2;z+1/2,-x+1/2,-y;-y,z+1/2,-x+1/2;"
	symLines[509] += "-z+1/2,-x,y+1/2;-z,x+1/2,-y+1/2;y+1/2,-z+1/2,-x;y-1/4,x+1/4,-z+1/4;-y-1/4,-x-1/4,-z-1/4;-x+1/4,z-1/4,y+1/4;"
	symLines[509] += "-x-1/4,-z-1/4,-y-1/4;z+1/4,-y+1/4,x-1/4;-z-1/4,-y-1/4,-x-1/4"
	symLines[510]  = "x,y,z;-y+1/4,x-1/4,z+1/4;-x,-y+1/2,z;y+1/4,-x+1/4,z-1/4;x+1/4,-z+1/4,y-1/4;x,-y,-z+1/2;x-1/4,z+1/4,-y+1/4;"
	symLines[510] += "z-1/4,y+1/4,-x+1/4;-x+1/2,y,-z;-z+1/4,y-1/4,x+1/4;z,x,y;y,z,x;-y,-z+1/2,x;z,-x,-y+1/2;-y+1/2,z,-x;-z,-x+1/2,y;"
	symLines[510] += "-z+1/2,x,-y;y,-z,-x+1/2;y-1/4,x+1/4,-z+1/4;-y+1/4,-x+1/4,-z+1/4;-x+1/4,z-1/4,y+1/4;-x+1/4,-z+1/4,-y+1/4;"
	symLines[510] += "z+1/4,-y+1/4,x-1/4;-z+1/4,-y+1/4,-x+1/4;x+1/2,y+1/2,z+1/2;-y-1/4,x+1/4,z-1/4;-x+1/2,-y,z+1/2;y-1/4,-x-1/4,z+1/4;"
	symLines[510] += "x-1/4,-z-1/4,y+1/4;x+1/2,-y+1/2,-z;x+1/4,z-1/4,-y-1/4;z+1/4,y-1/4,-x-1/4;-x,y+1/2,-z+1/2;-z-1/4,y+1/4,x-1/4;"
	symLines[510] += "z+1/2,x+1/2,y+1/2;y+1/2,z+1/2,x+1/2;-y+1/2,-z,x+1/2;z+1/2,-x+1/2,-y;-y,z+1/2,-x+1/2;-z+1/2,-x,y+1/2;-z,x+1/2,-y+1/2;"
	symLines[510] += "y+1/2,-z+1/2,-x;y+1/4,x-1/4,-z-1/4;-y-1/4,-x-1/4,-z-1/4;-x-1/4,z+1/4,y-1/4;-x-1/4,-z-1/4,-y-1/4;z-1/4,-y-1/4,x+1/4;-z-1/4,-y-1/4,-x-1/4"
	symLines[511]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;-x,z,-y;x,-y,-z;-x,-z,y;-z,-y,x;-x,y,-z;z,-y,-x;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;"
	symLines[511] += "-z,-x,y;-z,x,-y;y,-z,-x;-y,-x,z;y,x,z;x,-z,-y;x,z,y;-z,y,-x;z,y,x"
	symLines[512]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;-x,z,-y;x,-y,-z;-x,-z,y;-z,-y,x;-x,y,-z;z,-y,-x;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;"
	symLines[512] += "-z,-x,y;-z,x,-y;y,-z,-x;-y,-x,z;y,x,z;x,-z,-y;x,z,y;-z,y,-x;z,y,x;x,y+1/2,z+1/2;y,-x+1/2,-z+1/2;-x,-y+1/2,z+1/2;"
	symLines[512] += "-y,x+1/2,-z+1/2;-x,z+1/2,-y+1/2;x,-y+1/2,-z+1/2;-x,-z+1/2,y+1/2;-z,-y+1/2,x+1/2;-x,y+1/2,-z+1/2;z,-y+1/2,-x+1/2;"
	symLines[512] += "z,x+1/2,y+1/2;y,z+1/2,x+1/2;-y,-z+1/2,x+1/2;z,-x+1/2,-y+1/2;-y,z+1/2,-x+1/2;-z,-x+1/2,y+1/2;-z,x+1/2,-y+1/2;"
	symLines[512] += "y,-z+1/2,-x+1/2;-y,-x+1/2,z+1/2;y,x+1/2,z+1/2;x,-z+1/2,-y+1/2;x,z+1/2,y+1/2;-z,y+1/2,-x+1/2;z,y+1/2,x+1/2;x+1/2,y,z+1/2;"
	symLines[512] += "y+1/2,-x,-z+1/2;-x+1/2,-y,z+1/2;-y+1/2,x,-z+1/2;-x+1/2,z,-y+1/2;x+1/2,-y,-z+1/2;-x+1/2,-z,y+1/2;-z+1/2,-y,x+1/2;"
	symLines[512] += "-x+1/2,y,-z+1/2;z+1/2,-y,-x+1/2;z+1/2,x,y+1/2;y+1/2,z,x+1/2;-y+1/2,-z,x+1/2;z+1/2,-x,-y+1/2;-y+1/2,z,-x+1/2;"
	symLines[512] += "-z+1/2,-x,y+1/2;-z+1/2,x,-y+1/2;y+1/2,-z,-x+1/2;-y+1/2,-x,z+1/2;y+1/2,x,z+1/2;x+1/2,-z,-y+1/2;x+1/2,z,y+1/2;"
	symLines[512] += "-z+1/2,y,-x+1/2;z+1/2,y,x+1/2;x+1/2,y+1/2,z;y+1/2,-x+1/2,-z;-x+1/2,-y+1/2,z;-y+1/2,x+1/2,-z;-x+1/2,z+1/2,-y;"
	symLines[512] += "x+1/2,-y+1/2,-z;-x+1/2,-z+1/2,y;-z+1/2,-y+1/2,x;-x+1/2,y+1/2,-z;z+1/2,-y+1/2,-x;z+1/2,x+1/2,y;y+1/2,z+1/2,x;"
	symLines[512] += "-y+1/2,-z+1/2,x;z+1/2,-x+1/2,-y;-y+1/2,z+1/2,-x;-z+1/2,-x+1/2,y;-z+1/2,x+1/2,-y;y+1/2,-z+1/2,-x;-y+1/2,-x+1/2,z;"
	symLines[512] += "y+1/2,x+1/2,z;x+1/2,-z+1/2,-y;x+1/2,z+1/2,y;-z+1/2,y+1/2,-x;z+1/2,y+1/2,x"
	symLines[513]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;-x,z,-y;x,-y,-z;-x,-z,y;-z,-y,x;-x,y,-z;z,-y,-x;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;"
	symLines[513] += "-z,-x,y;-z,x,-y;y,-z,-x;-y,-x,z;y,x,z;x,-z,-y;x,z,y;-z,y,-x;z,y,x;x+1/2,y+1/2,z+1/2;y+1/2,-x+1/2,-z+1/2;"
	symLines[513] += "-x+1/2,-y+1/2,z+1/2;-y+1/2,x+1/2,-z+1/2;-x+1/2,z+1/2,-y+1/2;x+1/2,-y+1/2,-z+1/2;-x+1/2,-z+1/2,y+1/2;-z+1/2,-y+1/2,x+1/2;"
	symLines[513] += "-x+1/2,y+1/2,-z+1/2;z+1/2,-y+1/2,-x+1/2;z+1/2,x+1/2,y+1/2;y+1/2,z+1/2,x+1/2;-y+1/2,-z+1/2,x+1/2;z+1/2,-x+1/2,-y+1/2;"
	symLines[513] += "-y+1/2,z+1/2,-x+1/2;-z+1/2,-x+1/2,y+1/2;-z+1/2,x+1/2,-y+1/2;y+1/2,-z+1/2,-x+1/2;-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2;"
	symLines[513] += "x+1/2,-z+1/2,-y+1/2;x+1/2,z+1/2,y+1/2;-z+1/2,y+1/2,-x+1/2;z+1/2,y+1/2,x+1/2"
	symLines[514]  = "x,y,z;y+1/2,-x+1/2,-z+1/2;-x,-y,z;-y+1/2,x+1/2,-z+1/2;-x+1/2,z+1/2,-y+1/2;x,-y,-z;-x+1/2,-z+1/2,y+1/2;"
	symLines[514] += "-z+1/2,-y+1/2,x+1/2;-x,y,-z;z+1/2,-y+1/2,-x+1/2;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;"
	symLines[514] += "-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2;x+1/2,-z+1/2,-y+1/2;x+1/2,z+1/2,y+1/2;-z+1/2,y+1/2,-x+1/2;z+1/2,y+1/2,x+1/2"
	symLines[515]  = "x,y,z;y,-x,-z+1/2;-x,-y,z;-y,x,-z+1/2;-x,z,-y+1/2;x,-y,-z;-x,-z,y+1/2;-z,-y,x+1/2;-x,y,-z;z,-y,-x+1/2;z,x,y;y,z,x;"
	symLines[515] += "-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;-y,-x,z+1/2;y,x,z+1/2;x,-z,-y+1/2;x,z,y+1/2;-z,y,-x+1/2;z,y,x+1/2;"
	symLines[515] += "x,y+1/2,z+1/2;y,-x+1/2,-z;-x,-y+1/2,z+1/2;-y,x+1/2,-z;-x,z+1/2,-y;x,-y+1/2,-z+1/2;-x,-z+1/2,y;-z,-y+1/2,x;"
	symLines[515] += "-x,y+1/2,-z+1/2;z,-y+1/2,-x;z,x+1/2,y+1/2;y,z+1/2,x+1/2;-y,-z+1/2,x+1/2;z,-x+1/2,-y+1/2;-y,z+1/2,-x+1/2;-z,-x+1/2,y+1/2;"
	symLines[515] += "-z,x+1/2,-y+1/2;y,-z+1/2,-x+1/2;-y,-x+1/2,z;y,x+1/2,z;x,-z+1/2,-y;x,z+1/2,y;-z,y+1/2,-x;z,y+1/2,x;x+1/2,y,z+1/2;"
	symLines[515] += "y+1/2,-x,-z;-x+1/2,-y,z+1/2;-y+1/2,x,-z;-x+1/2,z,-y;x+1/2,-y,-z+1/2;-x+1/2,-z,y;-z+1/2,-y,x;-x+1/2,y,-z+1/2;z+1/2,-y,-x;"
	symLines[515] += "z+1/2,x,y+1/2;y+1/2,z,x+1/2;-y+1/2,-z,x+1/2;z+1/2,-x,-y+1/2;-y+1/2,z,-x+1/2;-z+1/2,-x,y+1/2;-z+1/2,x,-y+1/2;"
	symLines[515] += "y+1/2,-z,-x+1/2;-y+1/2,-x,z;y+1/2,x,z;x+1/2,-z,-y;x+1/2,z,y;-z+1/2,y,-x;z+1/2,y,x;x+1/2,y+1/2,z;y+1/2,-x+1/2,-z+1/2;"
	symLines[515] += "-x+1/2,-y+1/2,z;-y+1/2,x+1/2,-z+1/2;-x+1/2,z+1/2,-y+1/2;x+1/2,-y+1/2,-z;-x+1/2,-z+1/2,y+1/2;-z+1/2,-y+1/2,x+1/2;"
	symLines[515] += "-x+1/2,y+1/2,-z;z+1/2,-y+1/2,-x+1/2;z+1/2,x+1/2,y;y+1/2,z+1/2,x;-y+1/2,-z+1/2,x;z+1/2,-x+1/2,-y;-y+1/2,z+1/2,-x;"
	symLines[515] += "-z+1/2,-x+1/2,y;-z+1/2,x+1/2,-y;y+1/2,-z+1/2,-x;-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2;x+1/2,-z+1/2,-y+1/2;"
	symLines[515] += "x+1/2,z+1/2,y+1/2;-z+1/2,y+1/2,-x+1/2;z+1/2,y+1/2,x+1/2"
	symLines[516]  = "x,y,z;y+1/4,-x-1/4,-z+1/4;-x,-y+1/2,z;-y+1/4,x+1/4,-z-1/4;-x+1/4,z+1/4,-y-1/4;x,-y,-z+1/2;-x-1/4,-z+1/4,y+1/4;"
	symLines[516] += "-z-1/4,-y+1/4,x+1/4;-x+1/2,y,-z;z+1/4,-y-1/4,-x+1/4;z,x,y;y,z,x;-y,-z+1/2,x;z,-x,-y+1/2;-y+1/2,z,-x;-z,-x+1/2,y;"
	symLines[516] += "-z+1/2,x,-y;y,-z,-x+1/2;-y-1/4,-x+1/4,z+1/4;y+1/4,x+1/4,z+1/4;x+1/4,-z-1/4,-y+1/4;x+1/4,z+1/4,y+1/4;-z+1/4,y+1/4,-x-1/4;"
	symLines[516] += "z+1/4,y+1/4,x+1/4;x+1/2,y+1/2,z+1/2;y-1/4,-x+1/4,-z-1/4;-x+1/2,-y,z+1/2;-y-1/4,x-1/4,-z+1/4;-x-1/4,z-1/4,-y+1/4;"
	symLines[516] += "x+1/2,-y+1/2,-z;-x+1/4,-z-1/4,y-1/4;-z+1/4,-y-1/4,x-1/4;-x,y+1/2,-z+1/2;z-1/4,-y+1/4,-x-1/4;z+1/2,x+1/2,y+1/2;"
	symLines[516] += "y+1/2,z+1/2,x+1/2;-y+1/2,-z,x+1/2;z+1/2,-x+1/2,-y;-y,z+1/2,-x+1/2;-z+1/2,-x,y+1/2;-z,x+1/2,-y+1/2;y+1/2,-z+1/2,-x;"
	symLines[516] += "-y+1/4,-x-1/4,z-1/4;y-1/4,x-1/4,z-1/4;x-1/4,-z+1/4,-y-1/4;x-1/4,z-1/4,y-1/4;-z-1/4,y-1/4,-x+1/4;z-1/4,y-1/4,x-1/4"
	symLines[517]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x,-z,y;x,-y,-z;x,z,-y;z,y,-x;-x,y,-z;-z,y,x;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;"
	symLines[517] += "-z,x,-y;y,-z,-x;y,x,-z;-y,-x,-z;-x,z,y;-x,-z,-y;z,-y,x;-z,-y,-x;-x,-y,-z;y,-x,-z;x,y,-z;-y,x,-z;-x,z,-y;-x,y,z;-x,-z,y;"
	symLines[517] += "-z,-y,x;x,-y,z;z,-y,-x;-z,-x,-y;-y,-z,-x;y,z,-x;-z,x,y;y,-z,x;z,x,-y;z,-x,y;-y,z,x;-y,-x,z;y,x,z;x,-z,-y;x,z,y;-z,y,-x;z,y,x"
	symLines[518]  = "x,y,z;-x+1/2,-y+1/2,-z+1/2;-y,x,z;-x,-y,z;y,-x,z;x,-z,y;x,-y,-z;x,z,-y;z,y,-x;-x,y,-z;-z,y,x;y+1/2,-x+1/2,-z+1/2;"
	symLines[518] += "-y+1/2,x+1/2,-z+1/2;-x+1/2,z+1/2,-y+1/2;-x+1/2,-z+1/2,y+1/2;-z+1/2,-y+1/2,x+1/2;z+1/2,-y+1/2,-x+1/2;z,x,y;y,z,x;-y,-z,x;"
	symLines[518] += "z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;-z+1/2,-x+1/2,-y+1/2;-y+1/2,-z+1/2,-x+1/2;y+1/2,z+1/2,-x+1/2;-z+1/2,x+1/2,y+1/2;"
	symLines[518] += "y+1/2,-z+1/2,x+1/2;z+1/2,x+1/2,-y+1/2;z+1/2,-x+1/2,y+1/2;-y+1/2,z+1/2,x+1/2;y,x,-z;-y,-x,-z;-x,z,y;-x,-z,-y;z,-y,x;"
	symLines[518] += "-z,-y,-x;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2;-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2;"
	symLines[518] += "x+1/2,-z+1/2,-y+1/2;x+1/2,z+1/2,y+1/2;-z+1/2,y+1/2,-x+1/2;z+1/2,y+1/2,x+1/2"
	symLines[519]  = "x,y,z;-y+1/2,x,z;-x+1/2,-y+1/2,z;y,-x+1/2,z;x,-z+1/2,y;x,-y+1/2,-z+1/2;x,z,-y+1/2;z,y,-x+1/2;-x+1/2,y,-z+1/2;-z+1/2,y,x;"
	symLines[519] += "z,x,y;y,z,x;-y+1/2,-z+1/2,x;z,-x+1/2,-y+1/2;-y+1/2,z,-x+1/2;-z+1/2,-x+1/2,y;-z+1/2,x,-y+1/2;y,-z+1/2,-x+1/2;y,x,-z+1/2;"
	symLines[519] += "-y+1/2,-x+1/2,-z+1/2;-x+1/2,z,y;-x+1/2,-z+1/2,-y+1/2;z,-y+1/2,x;-z+1/2,-y+1/2,-x+1/2;-x,-y,-z;y+1/2,-x,-z;x+1/2,y+1/2,-z;"
	symLines[519] += "-y,x+1/2,-z;-x,z+1/2,-y;-x,y+1/2,z+1/2;-x,-z,y+1/2;-z,-y,x+1/2;x+1/2,-y,z+1/2;z+1/2,-y,-x;-z,-x,-y;-y,-z,-x;"
	symLines[519] += "y+1/2,z+1/2,-x;-z,x+1/2,y+1/2;y+1/2,-z,x+1/2;z+1/2,x+1/2,-y;z+1/2,-x,y+1/2;-y,z+1/2,x+1/2;-y,-x,z+1/2;y+1/2,x+1/2,z+1/2;"
	symLines[519] += "x+1/2,-z,-y;x+1/2,z+1/2,y+1/2;-z,y+1/2,-x;z+1/2,y+1/2,x+1/2"
	symLines[520]  = "x,y,z;-y+1/2,x+1/2,z+1/2;-x,-y,z;y+1/2,-x+1/2,z+1/2;x+1/2,-z+1/2,y+1/2;x,-y,-z;x+1/2,z+1/2,-y+1/2;z+1/2,y+1/2,-x+1/2;"
	symLines[520] += "-x,y,-z;-z+1/2,y+1/2,x+1/2;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;y+1/2,x+1/2,-z+1/2;"
	symLines[520] += "-y+1/2,-x+1/2,-z+1/2;-x+1/2,z+1/2,y+1/2;-x+1/2,-z+1/2,-y+1/2;z+1/2,-y+1/2,x+1/2;-z+1/2,-y+1/2,-x+1/2;-x,-y,-z;"
	symLines[520] += "y+1/2,-x+1/2,-z+1/2;x,y,-z;-y+1/2,x+1/2,-z+1/2;-x+1/2,z+1/2,-y+1/2;-x,y,z;-x+1/2,-z+1/2,y+1/2;-z+1/2,-y+1/2,x+1/2;x,-y,z;"
	symLines[520] += "z+1/2,-y+1/2,-x+1/2;-z,-x,-y;-y,-z,-x;y,z,-x;-z,x,y;y,-z,x;z,x,-y;z,-x,y;-y,z,x;-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2;"
	symLines[520] += "x+1/2,-z+1/2,-y+1/2;x+1/2,z+1/2,y+1/2;-z+1/2,y+1/2,-x+1/2;z+1/2,y+1/2,x+1/2"
	symLines[521]  = "x,y,z;-x+1/2,-y+1/2,-z+1/2;-y+1/2,x+1/2,z+1/2;-x,-y,z;y+1/2,-x+1/2,z+1/2;x+1/2,-z+1/2,y+1/2;x,-y,-z;x+1/2,z+1/2,-y+1/2;"
	symLines[521] += "z+1/2,y+1/2,-x+1/2;-x,y,-z;-z+1/2,y+1/2,x+1/2;y,-x,-z;-y,x,-z;-x,z,-y;-x,-z,y;-z,-y,x;z,-y,-x;z,x,y;y,z,x;-y,-z,x;"
	symLines[521] += "z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;-z+1/2,-x+1/2,-y+1/2;-y+1/2,-z+1/2,-x+1/2;y+1/2,z+1/2,-x+1/2;-z+1/2,x+1/2,y+1/2;"
	symLines[521] += "y+1/2,-z+1/2,x+1/2;z+1/2,x+1/2,-y+1/2;z+1/2,-x+1/2,y+1/2;-y+1/2,z+1/2,x+1/2;y+1/2,x+1/2,-z+1/2;-y+1/2,-x+1/2,-z+1/2;"
	symLines[521] += "-x+1/2,z+1/2,y+1/2;-x+1/2,-z+1/2,-y+1/2;z+1/2,-y+1/2,x+1/2;-z+1/2,-y+1/2,-x+1/2;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2;"
	symLines[521] += "x+1/2,-y+1/2,z+1/2;-y,-x,z;y,x,z;x,-z,-y;x,z,y;-z,y,-x;z,y,x"
	symLines[522]  = "x,y,z;-y,x+1/2,z+1/2;-x+1/2,-y+1/2,z;y+1/2,-x,z+1/2;x+1/2,-z,y+1/2;x,-y+1/2,-z+1/2;x+1/2,z+1/2,-y;z+1/2,y+1/2,-x;"
	symLines[522] += "-x+1/2,y,-z+1/2;-z,y+1/2,x+1/2;z,x,y;y,z,x;-y+1/2,-z+1/2,x;z,-x+1/2,-y+1/2;-y+1/2,z,-x+1/2;-z+1/2,-x+1/2,y;"
	symLines[522] += "-z+1/2,x,-y+1/2;y,-z+1/2,-x+1/2;y+1/2,x+1/2,-z;-y,-x,-z;-x,z+1/2,y+1/2;-x,-z,-y;z+1/2,-y,x+1/2;-z,-y,-x;-x,-y,-z;"
	symLines[522] += "y,-x+1/2,-z+1/2;x+1/2,y+1/2,-z;-y+1/2,x,-z+1/2;-x+1/2,z,-y+1/2;-x,y+1/2,z+1/2;-x+1/2,-z+1/2,y;-z+1/2,-y+1/2,x;"
	symLines[522] += "x+1/2,-y,z+1/2;z,-y+1/2,-x+1/2;-z,-x,-y;-y,-z,-x;y+1/2,z+1/2,-x;-z,x+1/2,y+1/2;y+1/2,-z,x+1/2;z+1/2,x+1/2,-y;"
	symLines[522] += "z+1/2,-x,y+1/2;-y,z+1/2,x+1/2;-y+1/2,-x+1/2,z;y,x,z;x,-z+1/2,-y+1/2;x,z,y;-z+1/2,y,-x+1/2;z,y,x"
	symLines[523]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x,-z,y;x,-y,-z;x,z,-y;z,y,-x;-x,y,-z;-z,y,x;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;"
	symLines[523] += "-z,x,-y;y,-z,-x;y,x,-z;-y,-x,-z;-x,z,y;-x,-z,-y;z,-y,x;-z,-y,-x;-x,-y,-z;y,-x,-z;x,y,-z;-y,x,-z;-x,z,-y;-x,y,z;-x,-z,y;"
	symLines[523] += "-z,-y,x;x,-y,z;z,-y,-x;-z,-x,-y;-y,-z,-x;y,z,-x;-z,x,y;y,-z,x;z,x,-y;z,-x,y;-y,z,x;-y,-x,z;y,x,z;x,-z,-y;x,z,y;-z,y,-x;"
	symLines[523] += "z,y,x;x,y+1/2,z+1/2;-y,x+1/2,z+1/2;-x,-y+1/2,z+1/2;y,-x+1/2,z+1/2;x,-z+1/2,y+1/2;x,-y+1/2,-z+1/2;x,z+1/2,-y+1/2;"
	symLines[523] += "z,y+1/2,-x+1/2;-x,y+1/2,-z+1/2;-z,y+1/2,x+1/2;z,x+1/2,y+1/2;y,z+1/2,x+1/2;-y,-z+1/2,x+1/2;z,-x+1/2,-y+1/2;"
	symLines[523] += "-y,z+1/2,-x+1/2;-z,-x+1/2,y+1/2;-z,x+1/2,-y+1/2;y,-z+1/2,-x+1/2;y,x+1/2,-z+1/2;-y,-x+1/2,-z+1/2;-x,z+1/2,y+1/2;"
	symLines[523] += "-x,-z+1/2,-y+1/2;z,-y+1/2,x+1/2;-z,-y+1/2,-x+1/2;-x,-y+1/2,-z+1/2;y,-x+1/2,-z+1/2;x,y+1/2,-z+1/2;-y,x+1/2,-z+1/2;"
	symLines[523] += "-x,z+1/2,-y+1/2;-x,y+1/2,z+1/2;-x,-z+1/2,y+1/2;-z,-y+1/2,x+1/2;x,-y+1/2,z+1/2;z,-y+1/2,-x+1/2;-z,-x+1/2,-y+1/2;"
	symLines[523] += "-y,-z+1/2,-x+1/2;y,z+1/2,-x+1/2;-z,x+1/2,y+1/2;y,-z+1/2,x+1/2;z,x+1/2,-y+1/2;z,-x+1/2,y+1/2;-y,z+1/2,x+1/2;"
	symLines[523] += "-y,-x+1/2,z+1/2;y,x+1/2,z+1/2;x,-z+1/2,-y+1/2;x,z+1/2,y+1/2;-z,y+1/2,-x+1/2;z,y+1/2,x+1/2;x+1/2,y,z+1/2;-y+1/2,x,z+1/2;"
	symLines[523] += "-x+1/2,-y,z+1/2;y+1/2,-x,z+1/2;x+1/2,-z,y+1/2;x+1/2,-y,-z+1/2;x+1/2,z,-y+1/2;z+1/2,y,-x+1/2;-x+1/2,y,-z+1/2;"
	symLines[523] += "-z+1/2,y,x+1/2;z+1/2,x,y+1/2;y+1/2,z,x+1/2;-y+1/2,-z,x+1/2;z+1/2,-x,-y+1/2;-y+1/2,z,-x+1/2;-z+1/2,-x,y+1/2;"
	symLines[523] += "-z+1/2,x,-y+1/2;y+1/2,-z,-x+1/2;y+1/2,x,-z+1/2;-y+1/2,-x,-z+1/2;-x+1/2,z,y+1/2;-x+1/2,-z,-y+1/2;z+1/2,-y,x+1/2;"
	symLines[523] += "-z+1/2,-y,-x+1/2;-x+1/2,-y,-z+1/2;y+1/2,-x,-z+1/2;x+1/2,y,-z+1/2;-y+1/2,x,-z+1/2;-x+1/2,z,-y+1/2;-x+1/2,y,z+1/2;"
	symLines[523] += "-x+1/2,-z,y+1/2;-z+1/2,-y,x+1/2;x+1/2,-y,z+1/2;z+1/2,-y,-x+1/2;-z+1/2,-x,-y+1/2;-y+1/2,-z,-x+1/2;y+1/2,z,-x+1/2;"
	symLines[523] += "-z+1/2,x,y+1/2;y+1/2,-z,x+1/2;z+1/2,x,-y+1/2;z+1/2,-x,y+1/2;-y+1/2,z,x+1/2;-y+1/2,-x,z+1/2;y+1/2,x,z+1/2;x+1/2,-z,-y+1/2;"
	symLines[523] += "x+1/2,z,y+1/2;-z+1/2,y,-x+1/2;z+1/2,y,x+1/2;x+1/2,y+1/2,z;-y+1/2,x+1/2,z;-x+1/2,-y+1/2,z;y+1/2,-x+1/2,z;x+1/2,-z+1/2,y;"
	symLines[523] += "x+1/2,-y+1/2,-z;x+1/2,z+1/2,-y;z+1/2,y+1/2,-x;-x+1/2,y+1/2,-z;-z+1/2,y+1/2,x;z+1/2,x+1/2,y;y+1/2,z+1/2,x;-y+1/2,-z+1/2,x;"
	symLines[523] += "z+1/2,-x+1/2,-y;-y+1/2,z+1/2,-x;-z+1/2,-x+1/2,y;-z+1/2,x+1/2,-y;y+1/2,-z+1/2,-x;y+1/2,x+1/2,-z;-y+1/2,-x+1/2,-z;"
	symLines[523] += "-x+1/2,z+1/2,y;-x+1/2,-z+1/2,-y;z+1/2,-y+1/2,x;-z+1/2,-y+1/2,-x;-x+1/2,-y+1/2,-z;y+1/2,-x+1/2,-z;x+1/2,y+1/2,-z;"
	symLines[523] += "-y+1/2,x+1/2,-z;-x+1/2,z+1/2,-y;-x+1/2,y+1/2,z;-x+1/2,-z+1/2,y;-z+1/2,-y+1/2,x;x+1/2,-y+1/2,z;z+1/2,-y+1/2,-x;"
	symLines[523] += "-z+1/2,-x+1/2,-y;-y+1/2,-z+1/2,-x;y+1/2,z+1/2,-x;-z+1/2,x+1/2,y;y+1/2,-z+1/2,x;z+1/2,x+1/2,-y;z+1/2,-x+1/2,y;"
	symLines[523] += "-y+1/2,z+1/2,x;-y+1/2,-x+1/2,z;y+1/2,x+1/2,z;x+1/2,-z+1/2,-y;x+1/2,z+1/2,y;-z+1/2,y+1/2,-x;z+1/2,y+1/2,x"
	symLines[524]  = "x,y,z;-y,x,z+1/2;-x,-y,z;y,-x,z+1/2;x,-z,y+1/2;x,-y,-z;x,z,-y+1/2;z,y,-x+1/2;-x,y,-z;-z,y,x+1/2;z,x,y;y,z,x;-y,-z,x;"
	symLines[524] += "z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;y,x,-z+1/2;-y,-x,-z+1/2;-x,z,y+1/2;-x,-z,-y+1/2;z,-y,x+1/2;-z,-y,-x+1/2;-x,-y,-z;"
	symLines[524] += "y,-x,-z+1/2;x,y,-z;-y,x,-z+1/2;-x,z,-y+1/2;-x,y,z;-x,-z,y+1/2;-z,-y,x+1/2;x,-y,z;z,-y,-x+1/2;-z,-x,-y;-y,-z,-x;y,z,-x;"
	symLines[524] += "-z,x,y;y,-z,x;z,x,-y;z,-x,y;-y,z,x;-y,-x,z+1/2;y,x,z+1/2;x,-z,-y+1/2;x,z,y+1/2;-z,y,-x+1/2;z,y,x+1/2;x,y+1/2,z+1/2;"
	symLines[524] += "-y,x+1/2,z;-x,-y+1/2,z+1/2;y,-x+1/2,z;x,-z+1/2,y;x,-y+1/2,-z+1/2;x,z+1/2,-y;z,y+1/2,-x;-x,y+1/2,-z+1/2;-z,y+1/2,x;"
	symLines[524] += "z,x+1/2,y+1/2;y,z+1/2,x+1/2;-y,-z+1/2,x+1/2;z,-x+1/2,-y+1/2;-y,z+1/2,-x+1/2;-z,-x+1/2,y+1/2;-z,x+1/2,-y+1/2;"
	symLines[524] += "y,-z+1/2,-x+1/2;y,x+1/2,-z;-y,-x+1/2,-z;-x,z+1/2,y;-x,-z+1/2,-y;z,-y+1/2,x;-z,-y+1/2,-x;-x,-y+1/2,-z+1/2;y,-x+1/2,-z;"
	symLines[524] += "x,y+1/2,-z+1/2;-y,x+1/2,-z;-x,z+1/2,-y;-x,y+1/2,z+1/2;-x,-z+1/2,y;-z,-y+1/2,x;x,-y+1/2,z+1/2;z,-y+1/2,-x;"
	symLines[524] += "-z,-x+1/2,-y+1/2;-y,-z+1/2,-x+1/2;y,z+1/2,-x+1/2;-z,x+1/2,y+1/2;y,-z+1/2,x+1/2;z,x+1/2,-y+1/2;z,-x+1/2,y+1/2;"
	symLines[524] += "-y,z+1/2,x+1/2;-y,-x+1/2,z;y,x+1/2,z;x,-z+1/2,-y;x,z+1/2,y;-z,y+1/2,-x;z,y+1/2,x;x+1/2,y,z+1/2;-y+1/2,x,z;"
	symLines[524] += "-x+1/2,-y,z+1/2;y+1/2,-x,z;x+1/2,-z,y;x+1/2,-y,-z+1/2;x+1/2,z,-y;z+1/2,y,-x;-x+1/2,y,-z+1/2;-z+1/2,y,x;z+1/2,x,y+1/2;"
	symLines[524] += "y+1/2,z,x+1/2;-y+1/2,-z,x+1/2;z+1/2,-x,-y+1/2;-y+1/2,z,-x+1/2;-z+1/2,-x,y+1/2;-z+1/2,x,-y+1/2;y+1/2,-z,-x+1/2;y+1/2,x,-z;"
	symLines[524] += "-y+1/2,-x,-z;-x+1/2,z,y;-x+1/2,-z,-y;z+1/2,-y,x;-z+1/2,-y,-x;-x+1/2,-y,-z+1/2;y+1/2,-x,-z;x+1/2,y,-z+1/2;-y+1/2,x,-z;"
	symLines[524] += "-x+1/2,z,-y;-x+1/2,y,z+1/2;-x+1/2,-z,y;-z+1/2,-y,x;x+1/2,-y,z+1/2;z+1/2,-y,-x;-z+1/2,-x,-y+1/2;-y+1/2,-z,-x+1/2;"
	symLines[524] += "y+1/2,z,-x+1/2;-z+1/2,x,y+1/2;y+1/2,-z,x+1/2;z+1/2,x,-y+1/2;z+1/2,-x,y+1/2;-y+1/2,z,x+1/2;-y+1/2,-x,z;y+1/2,x,z;"
	symLines[524] += "x+1/2,-z,-y;x+1/2,z,y;-z+1/2,y,-x;z+1/2,y,x;x+1/2,y+1/2,z;-y+1/2,x+1/2,z+1/2;-x+1/2,-y+1/2,z;y+1/2,-x+1/2,z+1/2;"
	symLines[524] += "x+1/2,-z+1/2,y+1/2;x+1/2,-y+1/2,-z;x+1/2,z+1/2,-y+1/2;z+1/2,y+1/2,-x+1/2;-x+1/2,y+1/2,-z;-z+1/2,y+1/2,x+1/2;"
	symLines[524] += "z+1/2,x+1/2,y;y+1/2,z+1/2,x;-y+1/2,-z+1/2,x;z+1/2,-x+1/2,-y;-y+1/2,z+1/2,-x;-z+1/2,-x+1/2,y;-z+1/2,x+1/2,-y;"
	symLines[524] += "y+1/2,-z+1/2,-x;y+1/2,x+1/2,-z+1/2;-y+1/2,-x+1/2,-z+1/2;-x+1/2,z+1/2,y+1/2;-x+1/2,-z+1/2,-y+1/2;z+1/2,-y+1/2,x+1/2;"
	symLines[524] += "-z+1/2,-y+1/2,-x+1/2;-x+1/2,-y+1/2,-z;y+1/2,-x+1/2,-z+1/2;x+1/2,y+1/2,-z;-y+1/2,x+1/2,-z+1/2;-x+1/2,z+1/2,-y+1/2;"
	symLines[524] += "-x+1/2,y+1/2,z;-x+1/2,-z+1/2,y+1/2;-z+1/2,-y+1/2,x+1/2;x+1/2,-y+1/2,z;z+1/2,-y+1/2,-x+1/2;-z+1/2,-x+1/2,-y;"
	symLines[524] += "-y+1/2,-z+1/2,-x;y+1/2,z+1/2,-x;-z+1/2,x+1/2,y;y+1/2,-z+1/2,x;z+1/2,x+1/2,-y;z+1/2,-x+1/2,y;-y+1/2,z+1/2,x;"
	symLines[524] += "-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2;x+1/2,-z+1/2,-y+1/2;x+1/2,z+1/2,y+1/2;-z+1/2,y+1/2,-x+1/2;z+1/2,y+1/2,x+1/2"
	symLines[525]  = "x,y,z;-x+1/4,-y+1/4,-z+1/4;-y+1/4,x+1/4,z+1/4;-x,-y,z;y+1/4,-x+1/4,z+1/4;x+1/4,-z+1/4,y+1/4;x,-y,-z;x+1/4,z+1/4,-y+1/4;"
	symLines[525] += "z+1/4,y+1/4,-x+1/4;-x,y,-z;-z+1/4,y+1/4,x+1/4;y,-x,-z;-y,x,-z;-x,z,-y;-x,-z,y;-z,-y,x;z,-y,-x;z,x,y;y,z,x;-y,-z,x;"
	symLines[525] += "z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;-z+1/4,-x+1/4,-y+1/4;-y+1/4,-z+1/4,-x+1/4;y+1/4,z+1/4,-x+1/4;-z+1/4,x+1/4,y+1/4;"
	symLines[525] += "y+1/4,-z+1/4,x+1/4;z+1/4,x+1/4,-y+1/4;z+1/4,-x+1/4,y+1/4;-y+1/4,z+1/4,x+1/4;y+1/4,x+1/4,-z+1/4;-y+1/4,-x+1/4,-z+1/4;"
	symLines[525] += "-x+1/4,z+1/4,y+1/4;-x+1/4,-z+1/4,-y+1/4;z+1/4,-y+1/4,x+1/4;-z+1/4,-y+1/4,-x+1/4;x+1/4,y+1/4,-z+1/4;-x+1/4,y+1/4,z+1/4;"
	symLines[525] += "x+1/4,-y+1/4,z+1/4;-y,-x,z;y,x,z;x,-z,-y;x,z,y;-z,y,-x;z,y,x;x,y+1/2,z+1/2;-x+1/4,-y-1/4,-z-1/4;-y+1/4,x-1/4,z-1/4;"
	symLines[525] += "-x,-y+1/2,z+1/2;y+1/4,-x-1/4,z-1/4;x+1/4,-z-1/4,y-1/4;x,-y+1/2,-z+1/2;x+1/4,z-1/4,-y-1/4;z+1/4,y-1/4,-x-1/4;"
	symLines[525] += "-x,y+1/2,-z+1/2;-z+1/4,y-1/4,x-1/4;y,-x+1/2,-z+1/2;-y,x+1/2,-z+1/2;-x,z+1/2,-y+1/2;-x,-z+1/2,y+1/2;-z,-y+1/2,x+1/2;"
	symLines[525] += "z,-y+1/2,-x+1/2;z,x+1/2,y+1/2;y,z+1/2,x+1/2;-y,-z+1/2,x+1/2;z,-x+1/2,-y+1/2;-y,z+1/2,-x+1/2;-z,-x+1/2,y+1/2;"
	symLines[525] += "-z,x+1/2,-y+1/2;y,-z+1/2,-x+1/2;-z+1/4,-x-1/4,-y-1/4;-y+1/4,-z-1/4,-x-1/4;y+1/4,z-1/4,-x-1/4;-z+1/4,x-1/4,y-1/4;"
	symLines[525] += "y+1/4,-z-1/4,x-1/4;z+1/4,x-1/4,-y-1/4;z+1/4,-x-1/4,y-1/4;-y+1/4,z-1/4,x-1/4;y+1/4,x-1/4,-z-1/4;-y+1/4,-x-1/4,-z-1/4;"
	symLines[525] += "-x+1/4,z-1/4,y-1/4;-x+1/4,-z-1/4,-y-1/4;z+1/4,-y-1/4,x-1/4;-z+1/4,-y-1/4,-x-1/4;x+1/4,y-1/4,-z-1/4;-x+1/4,y-1/4,z-1/4;"
	symLines[525] += "x+1/4,-y-1/4,z-1/4;-y,-x+1/2,z+1/2;y,x+1/2,z+1/2;x,-z+1/2,-y+1/2;x,z+1/2,y+1/2;-z,y+1/2,-x+1/2;z,y+1/2,x+1/2;"
	symLines[525] += "x+1/2,y,z+1/2;-x-1/4,-y+1/4,-z-1/4;-y-1/4,x+1/4,z-1/4;-x+1/2,-y,z+1/2;y-1/4,-x+1/4,z-1/4;x-1/4,-z+1/4,y-1/4;"
	symLines[525] += "x+1/2,-y,-z+1/2;x-1/4,z+1/4,-y-1/4;z-1/4,y+1/4,-x-1/4;-x+1/2,y,-z+1/2;-z-1/4,y+1/4,x-1/4;y+1/2,-x,-z+1/2;-y+1/2,x,-z+1/2;"
	symLines[525] += "-x+1/2,z,-y+1/2;-x+1/2,-z,y+1/2;-z+1/2,-y,x+1/2;z+1/2,-y,-x+1/2;z+1/2,x,y+1/2;y+1/2,z,x+1/2;-y+1/2,-z,x+1/2;"
	symLines[525] += "z+1/2,-x,-y+1/2;-y+1/2,z,-x+1/2;-z+1/2,-x,y+1/2;-z+1/2,x,-y+1/2;y+1/2,-z,-x+1/2;-z-1/4,-x+1/4,-y-1/4;"
	symLines[525] += "-y-1/4,-z+1/4,-x-1/4;y-1/4,z+1/4,-x-1/4;-z-1/4,x+1/4,y-1/4;y-1/4,-z+1/4,x-1/4;z-1/4,x+1/4,-y-1/4;z-1/4,-x+1/4,y-1/4;"
	symLines[525] += "-y-1/4,z+1/4,x-1/4;y-1/4,x+1/4,-z-1/4;-y-1/4,-x+1/4,-z-1/4;-x-1/4,z+1/4,y-1/4;-x-1/4,-z+1/4,-y-1/4;z-1/4,-y+1/4,x-1/4;"
	symLines[525] += "-z-1/4,-y+1/4,-x-1/4;x-1/4,y+1/4,-z-1/4;-x-1/4,y+1/4,z-1/4;x-1/4,-y+1/4,z-1/4;-y+1/2,-x,z+1/2;y+1/2,x,z+1/2;"
	symLines[525] += "x+1/2,-z,-y+1/2;x+1/2,z,y+1/2;-z+1/2,y,-x+1/2;z+1/2,y,x+1/2;x+1/2,y+1/2,z;-x-1/4,-y-1/4,-z+1/4;-y-1/4,x-1/4,z+1/4;"
	symLines[525] += "-x+1/2,-y+1/2,z;y-1/4,-x-1/4,z+1/4;x-1/4,-z-1/4,y+1/4;x+1/2,-y+1/2,-z;x-1/4,z-1/4,-y+1/4;z-1/4,y-1/4,-x+1/4;"
	symLines[525] += "-x+1/2,y+1/2,-z;-z-1/4,y-1/4,x+1/4;y+1/2,-x+1/2,-z;-y+1/2,x+1/2,-z;-x+1/2,z+1/2,-y;-x+1/2,-z+1/2,y;-z+1/2,-y+1/2,x;"
	symLines[525] += "z+1/2,-y+1/2,-x;z+1/2,x+1/2,y;y+1/2,z+1/2,x;-y+1/2,-z+1/2,x;z+1/2,-x+1/2,-y;-y+1/2,z+1/2,-x;-z+1/2,-x+1/2,y;"
	symLines[525] += "-z+1/2,x+1/2,-y;y+1/2,-z+1/2,-x;-z-1/4,-x-1/4,-y+1/4;-y-1/4,-z-1/4,-x+1/4;y-1/4,z-1/4,-x+1/4;-z-1/4,x-1/4,y+1/4;"
	symLines[525] += "y-1/4,-z-1/4,x+1/4;z-1/4,x-1/4,-y+1/4;z-1/4,-x-1/4,y+1/4;-y-1/4,z-1/4,x+1/4;y-1/4,x-1/4,-z+1/4;-y-1/4,-x-1/4,-z+1/4;"
	symLines[525] += "-x-1/4,z-1/4,y+1/4;-x-1/4,-z-1/4,-y+1/4;z-1/4,-y-1/4,x+1/4;-z-1/4,-y-1/4,-x+1/4;x-1/4,y-1/4,-z+1/4;-x-1/4,y-1/4,z+1/4;"
	symLines[525] += "x-1/4,-y-1/4,z+1/4;-y+1/2,-x+1/2,z;y+1/2,x+1/2,z;x+1/2,-z+1/2,-y;x+1/2,z+1/2,y;-z+1/2,y+1/2,-x;z+1/2,y+1/2,x"
	symLines[526]  = "x,y,z;-y,x+1/4,z+1/4;-x+1/4,-y+1/4,z;y+1/4,-x,z+1/4;x+1/4,-z,y+1/4;x,-y+1/4,-z+1/4;x+1/4,z+1/4,-y;z+1/4,y+1/4,-x;"
	symLines[526] += "-x+1/4,y,-z+1/4;-z,y+1/4,x+1/4;z,x,y;y,z,x;-y+1/4,-z+1/4,x;z,-x+1/4,-y+1/4;-y+1/4,z,-x+1/4;-z+1/4,-x+1/4,y;"
	symLines[526] += "-z+1/4,x,-y+1/4;y,-z+1/4,-x+1/4;y+1/4,x+1/4,-z;-y,-x,-z;-x,z+1/4,y+1/4;-x,-z,-y;z+1/4,-y,x+1/4;-z,-y,-x;-x,-y,-z;"
	symLines[526] += "y,-x-1/4,-z-1/4;x-1/4,y-1/4,-z;-y-1/4,x,-z-1/4;-x-1/4,z,-y-1/4;-x,y-1/4,z-1/4;-x-1/4,-z-1/4,y;-z-1/4,-y-1/4,x;"
	symLines[526] += "x-1/4,-y,z-1/4;z,-y-1/4,-x-1/4;-z,-x,-y;-y,-z,-x;y-1/4,z-1/4,-x;-z,x-1/4,y-1/4;y-1/4,-z,x-1/4;z-1/4,x-1/4,-y;"
	symLines[526] += "z-1/4,-x,y-1/4;-y,z-1/4,x-1/4;-y-1/4,-x-1/4,z;y,x,z;x,-z-1/4,-y-1/4;x,z,y;-z-1/4,y,-x-1/4;z,y,x;x,y+1/2,z+1/2;"
	symLines[526] += "-y,x-1/4,z-1/4;-x+1/4,-y-1/4,z+1/2;y+1/4,-x+1/2,z-1/4;x+1/4,-z+1/2,y-1/4;x,-y-1/4,-z-1/4;x+1/4,z-1/4,-y+1/2;"
	symLines[526] += "z+1/4,y-1/4,-x+1/2;-x+1/4,y+1/2,-z-1/4;-z,y-1/4,x-1/4;z,x+1/2,y+1/2;y,z+1/2,x+1/2;-y+1/4,-z-1/4,x+1/2;z,-x-1/4,-y-1/4;"
	symLines[526] += "-y+1/4,z+1/2,-x-1/4;-z+1/4,-x-1/4,y+1/2;-z+1/4,x+1/2,-y-1/4;y,-z-1/4,-x-1/4;y+1/4,x-1/4,-z+1/2;-y,-x+1/2,-z+1/2;"
	symLines[526] += "-x,z-1/4,y-1/4;-x,-z+1/2,-y+1/2;z+1/4,-y+1/2,x-1/4;-z,-y+1/2,-x+1/2;-x,-y+1/2,-z+1/2;y,-x+1/4,-z+1/4;x-1/4,y+1/4,-z+1/2;"
	symLines[526] += "-y-1/4,x+1/2,-z+1/4;-x-1/4,z+1/2,-y+1/4;-x,y+1/4,z+1/4;-x-1/4,-z+1/4,y+1/2;-z-1/4,-y+1/4,x+1/2;x-1/4,-y+1/2,z+1/4;"
	symLines[526] += "z,-y+1/4,-x+1/4;-z,-x+1/2,-y+1/2;-y,-z+1/2,-x+1/2;y-1/4,z+1/4,-x+1/2;-z,x+1/4,y+1/4;y-1/4,-z+1/2,x+1/4;"
	symLines[526] += "z-1/4,x+1/4,-y+1/2;z-1/4,-x+1/2,y+1/4;-y,z+1/4,x+1/4;-y-1/4,-x+1/4,z+1/2;y,x+1/2,z+1/2;x,-z+1/4,-y+1/4;x,z+1/2,y+1/2;"
	symLines[526] += "-z-1/4,y+1/2,-x+1/4;z,y+1/2,x+1/2;x+1/2,y,z+1/2;-y+1/2,x+1/4,z-1/4;-x-1/4,-y+1/4,z+1/2;y-1/4,-x,z-1/4;x-1/4,-z,y-1/4;"
	symLines[526] += "x+1/2,-y+1/4,-z-1/4;x-1/4,z+1/4,-y+1/2;z-1/4,y+1/4,-x+1/2;-x-1/4,y,-z-1/4;-z+1/2,y+1/4,x-1/4;z+1/2,x,y+1/2;y+1/2,z,x+1/2;"
	symLines[526] += "-y-1/4,-z+1/4,x+1/2;z+1/2,-x+1/4,-y-1/4;-y-1/4,z,-x-1/4;-z-1/4,-x+1/4,y+1/2;-z-1/4,x,-y-1/4;y+1/2,-z+1/4,-x-1/4;"
	symLines[526] += "y-1/4,x+1/4,-z+1/2;-y+1/2,-x,-z+1/2;-x+1/2,z+1/4,y-1/4;-x+1/2,-z,-y+1/2;z-1/4,-y,x-1/4;-z+1/2,-y,-x+1/2;-x+1/2,-y,-z+1/2;"
	symLines[526] += "y+1/2,-x-1/4,-z+1/4;x+1/4,y-1/4,-z+1/2;-y+1/4,x,-z+1/4;-x+1/4,z,-y+1/4;-x+1/2,y-1/4,z+1/4;-x+1/4,-z-1/4,y+1/2;"
	symLines[526] += "-z+1/4,-y-1/4,x+1/2;x+1/4,-y,z+1/4;z+1/2,-y-1/4,-x+1/4;-z+1/2,-x,-y+1/2;-y+1/2,-z,-x+1/2;y+1/4,z-1/4,-x+1/2;"
	symLines[526] += "-z+1/2,x-1/4,y+1/4;y+1/4,-z,x+1/4;z+1/4,x-1/4,-y+1/2;z+1/4,-x,y+1/4;-y+1/2,z-1/4,x+1/4;-y+1/4,-x-1/4,z+1/2;y+1/2,x,z+1/2;"
	symLines[526] += "x+1/2,-z-1/4,-y+1/4;x+1/2,z,y+1/2;-z+1/4,y,-x+1/4;z+1/2,y,x+1/2;x+1/2,y+1/2,z;-y+1/2,x-1/4,z+1/4;-x-1/4,-y-1/4,z;"
	symLines[526] += "y-1/4,-x+1/2,z+1/4;x-1/4,-z+1/2,y+1/4;x+1/2,-y-1/4,-z+1/4;x-1/4,z-1/4,-y;z-1/4,y-1/4,-x;-x-1/4,y+1/2,-z+1/4;"
	symLines[526] += "-z+1/2,y-1/4,x+1/4;z+1/2,x+1/2,y;y+1/2,z+1/2,x;-y-1/4,-z-1/4,x;z+1/2,-x-1/4,-y+1/4;-y-1/4,z+1/2,-x+1/4;-z-1/4,-x-1/4,y;"
	symLines[526] += "-z-1/4,x+1/2,-y+1/4;y+1/2,-z-1/4,-x+1/4;y-1/4,x-1/4,-z;-y+1/2,-x+1/2,-z;-x+1/2,z-1/4,y+1/4;-x+1/2,-z+1/2,-y;"
	symLines[526] += "z-1/4,-y+1/2,x+1/4;-z+1/2,-y+1/2,-x;-x+1/2,-y+1/2,-z;y+1/2,-x+1/4,-z-1/4;x+1/4,y+1/4,-z;-y+1/4,x+1/2,-z-1/4;"
	symLines[526] += "-x+1/4,z+1/2,-y-1/4;-x+1/2,y+1/4,z-1/4;-x+1/4,-z+1/4,y;-z+1/4,-y+1/4,x;x+1/4,-y+1/2,z-1/4;z+1/2,-y+1/4,-x-1/4;"
	symLines[526] += "-z+1/2,-x+1/2,-y;-y+1/2,-z+1/2,-x;y+1/4,z+1/4,-x;-z+1/2,x+1/4,y-1/4;y+1/4,-z+1/2,x-1/4;z+1/4,x+1/4,-y;z+1/4,-x+1/2,y-1/4;"
	symLines[526] += "-y+1/2,z+1/4,x-1/4;-y+1/4,-x+1/4,z;y+1/2,x+1/2,z;x+1/2,-z+1/4,-y-1/4;x+1/2,z+1/2,y;-z+1/4,y+1/2,-x-1/4;z+1/2,y+1/2,x"
	symLines[527]  = "x,y,z;-x+1/4,-y+1/4,-z-1/4;-y+1/4,x+1/4,z+1/4;-x,-y,z;y+1/4,-x+1/4,z+1/4;x+1/4,-z+1/4,y+1/4;x,-y,-z;x+1/4,z+1/4,-y+1/4;"
	symLines[527] += "z+1/4,y+1/4,-x+1/4;-x,y,-z;-z+1/4,y+1/4,x+1/4;y,-x,-z+1/2;-y,x,-z+1/2;-x,z,-y+1/2;-x,-z,y+1/2;-z,-y,x+1/2;z,-y,-x+1/2;"
	symLines[527] += "z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;-z+1/4,-x+1/4,-y-1/4;-y+1/4,-z+1/4,-x-1/4;y+1/4,z+1/4,-x-1/4;"
	symLines[527] += "-z+1/4,x+1/4,y-1/4;y+1/4,-z+1/4,x-1/4;z+1/4,x+1/4,-y-1/4;z+1/4,-x+1/4,y-1/4;-y+1/4,z+1/4,x-1/4;y+1/4,x+1/4,-z+1/4;"
	symLines[527] += "-y+1/4,-x+1/4,-z+1/4;-x+1/4,z+1/4,y+1/4;-x+1/4,-z+1/4,-y+1/4;z+1/4,-y+1/4,x+1/4;-z+1/4,-y+1/4,-x+1/4;x+1/4,y+1/4,-z-1/4;"
	symLines[527] += "-x+1/4,y+1/4,z-1/4;x+1/4,-y+1/4,z-1/4;-y,-x,z+1/2;y,x,z+1/2;x,-z,-y+1/2;x,z,y+1/2;-z,y,-x+1/2;z,y,x+1/2;x,y+1/2,z+1/2;"
	symLines[527] += "-x+1/4,-y-1/4,-z+1/4;-y+1/4,x-1/4,z-1/4;-x,-y+1/2,z+1/2;y+1/4,-x-1/4,z-1/4;x+1/4,-z-1/4,y-1/4;x,-y+1/2,-z+1/2;"
	symLines[527] += "x+1/4,z-1/4,-y-1/4;z+1/4,y-1/4,-x-1/4;-x,y+1/2,-z+1/2;-z+1/4,y-1/4,x-1/4;y,-x+1/2,-z;-y,x+1/2,-z;-x,z+1/2,-y;-x,-z+1/2,y;"
	symLines[527] += "-z,-y+1/2,x;z,-y+1/2,-x;z,x+1/2,y+1/2;y,z+1/2,x+1/2;-y,-z+1/2,x+1/2;z,-x+1/2,-y+1/2;-y,z+1/2,-x+1/2;-z,-x+1/2,y+1/2;"
	symLines[527] += "-z,x+1/2,-y+1/2;y,-z+1/2,-x+1/2;-z+1/4,-x-1/4,-y+1/4;-y+1/4,-z-1/4,-x+1/4;y+1/4,z-1/4,-x+1/4;-z+1/4,x-1/4,y+1/4;"
	symLines[527] += "y+1/4,-z-1/4,x+1/4;z+1/4,x-1/4,-y+1/4;z+1/4,-x-1/4,y+1/4;-y+1/4,z-1/4,x+1/4;y+1/4,x-1/4,-z-1/4;-y+1/4,-x-1/4,-z-1/4;"
	symLines[527] += "-x+1/4,z-1/4,y-1/4;-x+1/4,-z-1/4,-y-1/4;z+1/4,-y-1/4,x-1/4;-z+1/4,-y-1/4,-x-1/4;x+1/4,y-1/4,-z+1/4;-x+1/4,y-1/4,z+1/4;"
	symLines[527] += "x+1/4,-y-1/4,z+1/4;-y,-x+1/2,z;y,x+1/2,z;x,-z+1/2,-y;x,z+1/2,y;-z,y+1/2,-x;z,y+1/2,x;x+1/2,y,z+1/2;-x-1/4,-y+1/4,-z+1/4;"
	symLines[527] += "-y-1/4,x+1/4,z-1/4;-x+1/2,-y,z+1/2;y-1/4,-x+1/4,z-1/4;x-1/4,-z+1/4,y-1/4;x+1/2,-y,-z+1/2;x-1/4,z+1/4,-y-1/4;"
	symLines[527] += "z-1/4,y+1/4,-x-1/4;-x+1/2,y,-z+1/2;-z-1/4,y+1/4,x-1/4;y+1/2,-x,-z;-y+1/2,x,-z;-x+1/2,z,-y;-x+1/2,-z,y;-z+1/2,-y,x;"
	symLines[527] += "z+1/2,-y,-x;z+1/2,x,y+1/2;y+1/2,z,x+1/2;-y+1/2,-z,x+1/2;z+1/2,-x,-y+1/2;-y+1/2,z,-x+1/2;-z+1/2,-x,y+1/2;-z+1/2,x,-y+1/2;"
	symLines[527] += "y+1/2,-z,-x+1/2;-z-1/4,-x+1/4,-y+1/4;-y-1/4,-z+1/4,-x+1/4;y-1/4,z+1/4,-x+1/4;-z-1/4,x+1/4,y+1/4;y-1/4,-z+1/4,x+1/4;"
	symLines[527] += "z-1/4,x+1/4,-y+1/4;z-1/4,-x+1/4,y+1/4;-y-1/4,z+1/4,x+1/4;y-1/4,x+1/4,-z-1/4;-y-1/4,-x+1/4,-z-1/4;-x-1/4,z+1/4,y-1/4;"
	symLines[527] += "-x-1/4,-z+1/4,-y-1/4;z-1/4,-y+1/4,x-1/4;-z-1/4,-y+1/4,-x-1/4;x-1/4,y+1/4,-z+1/4;-x-1/4,y+1/4,z+1/4;x-1/4,-y+1/4,z+1/4;"
	symLines[527] += "-y+1/2,-x,z;y+1/2,x,z;x+1/2,-z,-y;x+1/2,z,y;-z+1/2,y,-x;z+1/2,y,x;x+1/2,y+1/2,z;-x-1/4,-y-1/4,-z-1/4;-y-1/4,x-1/4,z+1/4;"
	symLines[527] += "-x+1/2,-y+1/2,z;y-1/4,-x-1/4,z+1/4;x-1/4,-z-1/4,y+1/4;x+1/2,-y+1/2,-z;x-1/4,z-1/4,-y+1/4;z-1/4,y-1/4,-x+1/4;"
	symLines[527] += "-x+1/2,y+1/2,-z;-z-1/4,y-1/4,x+1/4;y+1/2,-x+1/2,-z+1/2;-y+1/2,x+1/2,-z+1/2;-x+1/2,z+1/2,-y+1/2;-x+1/2,-z+1/2,y+1/2;"
	symLines[527] += "-z+1/2,-y+1/2,x+1/2;z+1/2,-y+1/2,-x+1/2;z+1/2,x+1/2,y;y+1/2,z+1/2,x;-y+1/2,-z+1/2,x;z+1/2,-x+1/2,-y;-y+1/2,z+1/2,-x;"
	symLines[527] += "-z+1/2,-x+1/2,y;-z+1/2,x+1/2,-y;y+1/2,-z+1/2,-x;-z-1/4,-x-1/4,-y-1/4;-y-1/4,-z-1/4,-x-1/4;y-1/4,z-1/4,-x-1/4;"
	symLines[527] += "-z-1/4,x-1/4,y-1/4;y-1/4,-z-1/4,x-1/4;z-1/4,x-1/4,-y-1/4;z-1/4,-x-1/4,y-1/4;-y-1/4,z-1/4,x-1/4;y-1/4,x-1/4,-z+1/4;"
	symLines[527] += "-y-1/4,-x-1/4,-z+1/4;-x-1/4,z-1/4,y+1/4;-x-1/4,-z-1/4,-y+1/4;z-1/4,-y-1/4,x+1/4;-z-1/4,-y-1/4,-x+1/4;x-1/4,y-1/4,-z-1/4;"
	symLines[527] += "-x-1/4,y-1/4,z-1/4;x-1/4,-y-1/4,z-1/4;-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2;x+1/2,-z+1/2,-y+1/2;x+1/2,z+1/2,y+1/2;"
	symLines[527] += "-z+1/2,y+1/2,-x+1/2;z+1/2,y+1/2,x+1/2"
	symLines[528]  = "x,y,z;-y,x+1/4,z-1/4;-x+1/4,-y+1/4,z;y+1/4,-x,z-1/4;x+1/4,-z,y-1/4;x,-y+1/4,-z+1/4;x+1/4,z-1/4,-y;z+1/4,y-1/4,-x;"
	symLines[528] += "-x+1/4,y,-z+1/4;-z,y+1/4,x-1/4;z,x,y;y,z,x;-y+1/4,-z+1/4,x;z,-x+1/4,-y+1/4;-y+1/4,z,-x+1/4;-z+1/4,-x+1/4,y;"
	symLines[528] += "-z+1/4,x,-y+1/4;y,-z+1/4,-x+1/4;y+1/4,x-1/4,-z;-y,-x,-z+1/2;-x,z+1/4,y-1/4;-x,-z,-y+1/2;z+1/4,-y,x-1/4;-z,-y,-x+1/2;"
	symLines[528] += "-x,-y,-z;y,-x-1/4,-z+1/4;x-1/4,y-1/4,-z;-y-1/4,x,-z+1/4;-x-1/4,z,-y+1/4;-x,y-1/4,z-1/4;-x-1/4,-z+1/4,y;-z-1/4,-y+1/4,x;"
	symLines[528] += "x-1/4,-y,z-1/4;z,-y-1/4,-x+1/4;-z,-x,-y;-y,-z,-x;y-1/4,z-1/4,-x;-z,x-1/4,y-1/4;y-1/4,-z,x-1/4;z-1/4,x-1/4,-y;"
	symLines[528] += "z-1/4,-x,y-1/4;-y,z-1/4,x-1/4;-y-1/4,-x+1/4,z;y,x,z+1/2;x,-z-1/4,-y+1/4;x,z,y+1/2;-z-1/4,y,-x+1/4;z,y,x+1/2;"
	symLines[528] += "x,y+1/2,z+1/2;-y,x-1/4,z+1/4;-x+1/4,-y-1/4,z+1/2;y+1/4,-x+1/2,z+1/4;x+1/4,-z+1/2,y+1/4;x,-y-1/4,-z-1/4;"
	symLines[528] += "x+1/4,z+1/4,-y+1/2;z+1/4,y+1/4,-x+1/2;-x+1/4,y+1/2,-z-1/4;-z,y-1/4,x+1/4;z,x+1/2,y+1/2;y,z+1/2,x+1/2;-y+1/4,-z-1/4,x+1/2;"
	symLines[528] += "z,-x-1/4,-y-1/4;-y+1/4,z+1/2,-x-1/4;-z+1/4,-x-1/4,y+1/2;-z+1/4,x+1/2,-y-1/4;y,-z-1/4,-x-1/4;y+1/4,x+1/4,-z+1/2;"
	symLines[528] += "-y,-x+1/2,-z;-x,z-1/4,y+1/4;-x,-z+1/2,-y;z+1/4,-y+1/2,x+1/4;-z,-y+1/2,-x;-x,-y+1/2,-z+1/2;y,-x+1/4,-z-1/4;"
	symLines[528] += "x-1/4,y+1/4,-z+1/2;-y-1/4,x+1/2,-z-1/4;-x-1/4,z+1/2,-y-1/4;-x,y+1/4,z+1/4;-x-1/4,-z-1/4,y+1/2;-z-1/4,-y-1/4,x+1/2;"
	symLines[528] += "x-1/4,-y+1/2,z+1/4;z,-y+1/4,-x-1/4;-z,-x+1/2,-y+1/2;-y,-z+1/2,-x+1/2;y-1/4,z+1/4,-x+1/2;-z,x+1/4,y+1/4;"
	symLines[528] += "y-1/4,-z+1/2,x+1/4;z-1/4,x+1/4,-y+1/2;z-1/4,-x+1/2,y+1/4;-y,z+1/4,x+1/4;-y-1/4,-x-1/4,z+1/2;y,x+1/2,z;x,-z+1/4,-y-1/4;"
	symLines[528] += "x,z+1/2,y;-z-1/4,y+1/2,-x-1/4;z,y+1/2,x;x+1/2,y,z+1/2;-y+1/2,x+1/4,z+1/4;-x-1/4,-y+1/4,z+1/2;y-1/4,-x,z+1/4;"
	symLines[528] += "x-1/4,-z,y+1/4;x+1/2,-y+1/4,-z-1/4;x-1/4,z-1/4,-y+1/2;z-1/4,y-1/4,-x+1/2;-x-1/4,y,-z-1/4;-z+1/2,y+1/4,x+1/4;"
	symLines[528] += "z+1/2,x,y+1/2;y+1/2,z,x+1/2;-y-1/4,-z+1/4,x+1/2;z+1/2,-x+1/4,-y-1/4;-y-1/4,z,-x-1/4;-z-1/4,-x+1/4,y+1/2;-z-1/4,x,-y-1/4;"
	symLines[528] += "y+1/2,-z+1/4,-x-1/4;y-1/4,x-1/4,-z+1/2;-y+1/2,-x,-z;-x+1/2,z+1/4,y+1/4;-x+1/2,-z,-y;z-1/4,-y,x+1/4;-z+1/2,-y,-x;"
	symLines[528] += "-x+1/2,-y,-z+1/2;y+1/2,-x-1/4,-z-1/4;x+1/4,y-1/4,-z+1/2;-y+1/4,x,-z-1/4;-x+1/4,z,-y-1/4;-x+1/2,y-1/4,z+1/4;"
	symLines[528] += "-x+1/4,-z+1/4,y+1/2;-z+1/4,-y+1/4,x+1/2;x+1/4,-y,z+1/4;z+1/2,-y-1/4,-x-1/4;-z+1/2,-x,-y+1/2;-y+1/2,-z,-x+1/2;"
	symLines[528] += "y+1/4,z-1/4,-x+1/2;-z+1/2,x-1/4,y+1/4;y+1/4,-z,x+1/4;z+1/4,x-1/4,-y+1/2;z+1/4,-x,y+1/4;-y+1/2,z-1/4,x+1/4;"
	symLines[528] += "-y+1/4,-x+1/4,z+1/2;y+1/2,x,z;x+1/2,-z-1/4,-y-1/4;x+1/2,z,y;-z+1/4,y,-x-1/4;z+1/2,y,x;x+1/2,y+1/2,z;-y+1/2,x-1/4,z-1/4;"
	symLines[528] += "-x-1/4,-y-1/4,z;y-1/4,-x+1/2,z-1/4;x-1/4,-z+1/2,y-1/4;x+1/2,-y-1/4,-z+1/4;x-1/4,z+1/4,-y;z-1/4,y+1/4,-x;"
	symLines[528] += "-x-1/4,y+1/2,-z+1/4;-z+1/2,y-1/4,x-1/4;z+1/2,x+1/2,y;y+1/2,z+1/2,x;-y-1/4,-z-1/4,x;z+1/2,-x-1/4,-y+1/4;"
	symLines[528] += "-y-1/4,z+1/2,-x+1/4;-z-1/4,-x-1/4,y;-z-1/4,x+1/2,-y+1/4;y+1/2,-z-1/4,-x+1/4;y-1/4,x+1/4,-z;-y+1/2,-x+1/2,-z+1/2;"
	symLines[528] += "-x+1/2,z-1/4,y-1/4;-x+1/2,-z+1/2,-y+1/2;z-1/4,-y+1/2,x-1/4;-z+1/2,-y+1/2,-x+1/2;-x+1/2,-y+1/2,-z;y+1/2,-x+1/4,-z+1/4;"
	symLines[528] += "x+1/4,y+1/4,-z;-y+1/4,x+1/2,-z+1/4;-x+1/4,z+1/2,-y+1/4;-x+1/2,y+1/4,z-1/4;-x+1/4,-z-1/4,y;-z+1/4,-y-1/4,x;"
	symLines[528] += "x+1/4,-y+1/2,z-1/4;z+1/2,-y+1/4,-x+1/4;-z+1/2,-x+1/2,-y;-y+1/2,-z+1/2,-x;y+1/4,z+1/4,-x;-z+1/2,x+1/4,y-1/4;"
	symLines[528] += "y+1/4,-z+1/2,x-1/4;z+1/4,x+1/4,-y;z+1/4,-x+1/2,y-1/4;-y+1/2,z+1/4,x-1/4;-y+1/4,-x-1/4,z;y+1/2,x+1/2,z+1/2;"
	symLines[528] += "x+1/2,-z+1/4,-y+1/4;x+1/2,z+1/2,y+1/2;-z+1/4,y+1/2,-x+1/4;z+1/2,y+1/2,x+1/2"
	symLines[529]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x,-z,y;x,-y,-z;x,z,-y;z,y,-x;-x,y,-z;-z,y,x;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;"
	symLines[529] += "-z,x,-y;y,-z,-x;y,x,-z;-y,-x,-z;-x,z,y;-x,-z,-y;z,-y,x;-z,-y,-x;-x,-y,-z;y,-x,-z;x,y,-z;-y,x,-z;-x,z,-y;-x,y,z;-x,-z,y;"
	symLines[529] += "-z,-y,x;x,-y,z;z,-y,-x;-z,-x,-y;-y,-z,-x;y,z,-x;-z,x,y;y,-z,x;z,x,-y;z,-x,y;-y,z,x;-y,-x,z;y,x,z;x,-z,-y;x,z,y;-z,y,-x;"
	symLines[529] += "z,y,x;x+1/2,y+1/2,z+1/2;-y+1/2,x+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;y+1/2,-x+1/2,z+1/2;x+1/2,-z+1/2,y+1/2;x+1/2,-y+1/2,-z+1/2;"
	symLines[529] += "x+1/2,z+1/2,-y+1/2;z+1/2,y+1/2,-x+1/2;-x+1/2,y+1/2,-z+1/2;-z+1/2,y+1/2,x+1/2;z+1/2,x+1/2,y+1/2;y+1/2,z+1/2,x+1/2;"
	symLines[529] += "-y+1/2,-z+1/2,x+1/2;z+1/2,-x+1/2,-y+1/2;-y+1/2,z+1/2,-x+1/2;-z+1/2,-x+1/2,y+1/2;-z+1/2,x+1/2,-y+1/2;y+1/2,-z+1/2,-x+1/2;"
	symLines[529] += "y+1/2,x+1/2,-z+1/2;-y+1/2,-x+1/2,-z+1/2;-x+1/2,z+1/2,y+1/2;-x+1/2,-z+1/2,-y+1/2;z+1/2,-y+1/2,x+1/2;-z+1/2,-y+1/2,-x+1/2;"
	symLines[529] += "-x+1/2,-y+1/2,-z+1/2;y+1/2,-x+1/2,-z+1/2;x+1/2,y+1/2,-z+1/2;-y+1/2,x+1/2,-z+1/2;-x+1/2,z+1/2,-y+1/2;-x+1/2,y+1/2,z+1/2;"
	symLines[529] += "-x+1/2,-z+1/2,y+1/2;-z+1/2,-y+1/2,x+1/2;x+1/2,-y+1/2,z+1/2;z+1/2,-y+1/2,-x+1/2;-z+1/2,-x+1/2,-y+1/2;-y+1/2,-z+1/2,-x+1/2;"
	symLines[529] += "y+1/2,z+1/2,-x+1/2;-z+1/2,x+1/2,y+1/2;y+1/2,-z+1/2,x+1/2;z+1/2,x+1/2,-y+1/2;z+1/2,-x+1/2,y+1/2;-y+1/2,z+1/2,x+1/2;"
	symLines[529] += "-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2;x+1/2,-z+1/2,-y+1/2;x+1/2,z+1/2,y+1/2;-z+1/2,y+1/2,-x+1/2;z+1/2,y+1/2,x+1/2"
	symLines[530]  = "x,y,z;-y+1/4,x-1/4,z+1/4;-x,-y+1/2,z;y+1/4,-x+1/4,z-1/4;x+1/4,-z+1/4,y-1/4;x,-y,-z+1/2;x-1/4,z+1/4,-y+1/4;"
	symLines[530] += "z-1/4,y+1/4,-x+1/4;-x+1/2,y,-z;-z+1/4,y-1/4,x+1/4;z,x,y;y,z,x;-y,-z+1/2,x;z,-x,-y+1/2;-y+1/2,z,-x;-z,-x+1/2,y;"
	symLines[530] += "-z+1/2,x,-y;y,-z,-x+1/2;y-1/4,x+1/4,-z+1/4;-y+1/4,-x+1/4,-z+1/4;-x+1/4,z-1/4,y+1/4;-x+1/4,-z+1/4,-y+1/4;"
	symLines[530] += "z+1/4,-y+1/4,x-1/4;-z+1/4,-y+1/4,-x+1/4;-x,-y,-z;y-1/4,-x+1/4,-z-1/4;x,y+1/2,-z;-y-1/4,x-1/4,-z+1/4;-x-1/4,z-1/4,-y+1/4;"
	symLines[530] += "-x,y,z+1/2;-x+1/4,-z-1/4,y-1/4;-z+1/4,-y-1/4,x-1/4;x+1/2,-y,z;z-1/4,-y+1/4,-x-1/4;-z,-x,-y;-y,-z,-x;y,z+1/2,-x;"
	symLines[530] += "-z,x,y+1/2;y+1/2,-z,x;z,x+1/2,-y;z+1/2,-x,y;-y,z,x+1/2;-y+1/4,-x-1/4,z-1/4;y-1/4,x-1/4,z-1/4;x-1/4,-z+1/4,-y-1/4;"
	symLines[530] += "x-1/4,z-1/4,y-1/4;-z-1/4,y-1/4,-x+1/4;z-1/4,y-1/4,x-1/4;x+1/2,y+1/2,z+1/2;-y-1/4,x+1/4,z-1/4;-x+1/2,-y,z+1/2;"
	symLines[530] += "y-1/4,-x-1/4,z+1/4;x-1/4,-z-1/4,y+1/4;x+1/2,-y+1/2,-z;x+1/4,z-1/4,-y-1/4;z+1/4,y-1/4,-x-1/4;-x,y+1/2,-z+1/2;"
	symLines[530] += "-z-1/4,y+1/4,x-1/4;z+1/2,x+1/2,y+1/2;y+1/2,z+1/2,x+1/2;-y+1/2,-z,x+1/2;z+1/2,-x+1/2,-y;-y,z+1/2,-x+1/2;-z+1/2,-x,y+1/2;"
	symLines[530] += "-z,x+1/2,-y+1/2;y+1/2,-z+1/2,-x;y+1/4,x-1/4,-z-1/4;-y-1/4,-x-1/4,-z-1/4;-x-1/4,z+1/4,y-1/4;-x-1/4,-z-1/4,-y-1/4;"
	symLines[530] += "z-1/4,-y-1/4,x+1/4;-z-1/4,-y-1/4,-x-1/4;-x+1/2,-y+1/2,-z+1/2;y+1/4,-x-1/4,-z+1/4;x+1/2,y,-z+1/2;-y+1/4,x+1/4,-z-1/4;"
	symLines[530] += "-x+1/4,z+1/4,-y-1/4;-x+1/2,y+1/2,z;-x-1/4,-z+1/4,y+1/4;-z-1/4,-y+1/4,x+1/4;x,-y+1/2,z+1/2;z+1/4,-y-1/4,-x+1/4;"
	symLines[530] += "-z+1/2,-x+1/2,-y+1/2;-y+1/2,-z+1/2,-x+1/2;y+1/2,z,-x+1/2;-z+1/2,x+1/2,y;y,-z+1/2,x+1/2;z+1/2,x,-y+1/2;z,-x+1/2,y+1/2;"
	symLines[530] += "-y+1/2,z+1/2,x;-y-1/4,-x+1/4,z+1/4;y+1/4,x+1/4,z+1/4;x+1/4,-z-1/4,-y+1/4;x+1/4,z+1/4,y+1/4;-z+1/4,y+1/4,-x-1/4;z+1/4,y+1/4,x+1/4"

	return symLines[idNum]
End






//
//	================================================================================
//
//			This section is for making the matricies for the symmetry operations for Wyckoff Symbols
//			Wyckoff info can be found at:	http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-wp-list
//
// 		the Funcitons associated with Wyckoff symbols are:
//			used by regular routines (above this):
//				WyckoffMultiplicity(), WyckoffMenuStr(), FindWyckoffSymbol(), ForceXYZtoWyckoff()
//
//			ONLY used by other Wyckoff routines (below this):
//				GetWyckoffSymStrings()

Static Function WyckoffMultiplicity(SpaceGroupID,letter)
	String SpaceGroupID
	String letter

	Wave/T Wlist=GetWyckoffSymStrings(SpaceGroupID)	// col0=letter, col1=symOp, col2=mult
	Variable i, N=DimSize(Wlist,0), mult=NaN
	for (i=0;i<N;i+=1)
		if (cmpstr(Wlist[i][0],letter)==0)
			mult = str2num(Wlist[i][2])
		endif
	endfor
	if (strsearch(SpaceGroupID,":R",0,2)>0)				// ONLY for Rhombohedral
		mult /= 3													//   divide multiplicity by 3
	endif
	return mult
End
//
Static Function/T WyckoffMenuStr(SpaceGroupID)
	String SpaceGroupID
	Wave/T Wlist=GetWyckoffSymStrings(SpaceGroupID)	// col0=letter, col1=symOp, col2=mult
	Variable i, N=DimSize(Wlist,0)
	String mStr=""
	for (i=0;i<N;i+=1)
		mStr += Wlist[i][0]+";"								// add each character to the menu str
	endfor
	return mStr
End
//
//Function testWyckoff(SG,letter)
//	Variable SG
//	String letter
//
//	String mStr = LatticeSym#WyckoffMenuStr(SG)
//	Variable mult = LatticeSym#WyckoffMultiplicity(SG,letter)
//	printf "mStr = '%s'\r",mStr
//	printf "for '%s',  mult = %g\r",letter,mult
//	return 0
//End










Static Function/T FindWyckoffSymbol(SpaceGroupID, x0,y0,z0, mult)
	String SpaceGroupID
	Variable x0,y0,z0
	Variable &mult

	Make/D/FREE vec4 = {x0,y0,z0,1}
	Make/D/FREE vec0 = {x0,y0,z0,0}
	Wave/T WyckList = GetWyckoffSymStrings(SpaceGroupID)
	String symbol=""
	mult = 0
	Variable m, N=DimSize(WyckList,0)
	for (m=0;m<N;m+=1)
		Wave mat = MatrixFromSymLine(WyckList[m][1], 4, zeroBad=0)
		MatrixOP/FREE diff = Sum(Abs( (mat x vec4) - vec0))
		if (diff[0]<1e-4)
			symbol = WyckList[m][0]
			mult = round(str2num(WyckList[m][2]))
			break
		endif
	endfor
	return symbol
End
//
//			for 47
//		test_FindWyckoffSymbol("47",  0, 0.5, 0.3765)
//
//			for 167:H
//		test_FindWyckoffSymbol("167:H", 0, 0, 0.25)		// "a", 6		0,0,1/4
//		test_FindWyckoffSymbol("167:H", 0,0,1/4)			// "a", 6		0,0,1/4
//		test_FindWyckoffSymbol("167:H", 0,0,0)				// "b", 6		0,0,0
//		test_FindWyckoffSymbol("167:H", 0,0,0.12)			// "c", 12		0,0,z
//		test_FindWyckoffSymbol("167:H", 1/2,0,0)			// "d", 18		1/2,0,0
//		test_FindWyckoffSymbol("167:H", 0.1,0,1/4)			// "e", 18		x,0,1/4
//		test_FindWyckoffSymbol("167:H", 0.3064, 0, 0.25)	// "e", 18		x,0,1/4
//		test_FindWyckoffSymbol("167:H", 0.1,0.1,1/4)		// "f", 36		x,y,z
//
//			for 167:R
//		test_FindWyckoffSymbol("167:R", 1/4,1/4,1/4)		// "a", 2		1/4,1/4,1/4
//		test_FindWyckoffSymbol("167:R", 0,0,0)				// "b", 2		0,0,0
//		test_FindWyckoffSymbol("167:R", 0.1,0.1,0.1)		// "c", 4		x,x,x
//		test_FindWyckoffSymbol("167:R", 1/2,0,1/2)			// "d", 6		1/2, 0, 1/2	
//		test_FindWyckoffSymbol("167:R", 1/2,0,0)			// "d", 6		1/2, 0, 0
//			xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//			xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//			xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//
//	Function test_FindWyckoffSymbol(SpaceGroupID, x0,y0,z0)
//		String SpaceGroupID
//		Variable x0,y0,z0
//		Variable mult
//		String letter = LatticeSym#FindWyckoffSymbol(SpaceGroupID, x0,y0,z0, mult)
//		printf "{%g, %g, %g} --> %s (%g),   for %s\r",x0,y0,z0, letter, mult, SpaceGroupID
//	End
//
//	for 		http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-wp-list
// for		http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-wp-list?gnum=167&grha=rhombohedral
// there appears to be a problem, also my values seem to be wrong.
//
//
//		167:R
//		a		1/2,1/4,0					 2		should be (1/4, 1/4, 1/4)
//		b		0,0,0							 2		OK
//		c		z*2,z,0						 4		should be (x,x,x)
//		d		1/2,0,1/2					 6		OK (would prefer (1/2,0,0)
//		e		x+1/2,-x*2+1/4,x			 6		should be (x,-x+1/2,1/4), (1/4,x,-x+1/2), (-x+1/2,1/4,x), (-x,x+1/2,3/4), (3/4,-x,x+1/2), or (x+1/2,3/4,-x)
//		f		x-y+z*2,-x*2+z,x-y*2	12		should be one of:
//														(x,y,z), (z,x,y), (y,z,x)
//														(-z+1/2,-y+1/2,-x+1/2), (-y+1/2,-x+1/2,-z+1/2), (-x+1/2,-z+1/2,-y+1/2)
//														(-x,-y,-z), (-z,-x,-y), (-y,-z,-x)
//														(z+1/2,y+1/2,x+1/2), (y+1/2,x+1/2,z+1/2), (x+1/2,z+1/2,y+1/2)





Static Function ForceXYZtoWyckoff(SpaceGroupID,symbol,x0,y0,z0)
	String SpaceGroupID
	String symbol
	Variable &x0,&y0,&z0
	Wave/T WyckList=GetWyckoffSymStrings(SpaceGroupID)

	String item, symOp=""
	Variable xop,yop,zop
	Variable m, N=DimSize(WyckList,0)
	for (m=0;m<N;m+=1)
		if (strsearch(WyckList[m][0], symbol,0)==0)
			symOp = WyckList[m][1]
		endif
	endfor
	if (strlen(symOp)<1)
		return 1
	endif
	symOp = ReplaceString("2x",symOp,"2*x")
	symOp = ReplaceString("2y",symOp,"2*y")
	symOp = ReplaceString("2z",symOp,"2*z")
	symOp = ReplaceString("x",symOp,num2str(x0))
	symOp = ReplaceString("y",symOp,num2str(y0))
	symOp = ReplaceString("z",symOp,num2str(z0))

	x0 = arithmetic(StringFromList(0,symOp,","))
	y0 = arithmetic(StringFromList(1,symOp,","))
	z0 = arithmetic(StringFromList(2,symOp,","))

	return 0
End









//  ======================================================================================  //
//  =========================== Start of Change Crystal Setting ==========================  //

Function ChangeSettingCurrentXtal(id, [printIt])	// change the setting of the CURRENT xtal
	String id										// the requested target id, if "", then prompt the user
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : !(!printIt)
	if (printIt && strlen(GetRTStackInfo(2)))
		printf "%sChangeSettingCurrentXtal(\"%s\")\r", BULLET,id
	endif

	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))	//fill the lattice structure with test values
		DoAlert 0, "no crystal structure found"
		return 1
	endif

	if (ChangeXtalSetting(xtal, id, printIt=printIt))	// failed to change the setting
		return 1
	endif

	// update everything, first the Panel values
	SVAR crystalStructStr = root:Packages:Lattices:PanelValues:crystalStructStr
	StructPut/S xtal, crystalStructStr

	NVAR dirty = root:Packages:Lattices:PanelValues:dirty
	SVAR desc = root:Packages:Lattices:PanelValues:desc
	SVAR SpaceGroupID = root:Packages:Lattices:PanelValues:SpaceGroupID
	NVAR a = root:Packages:Lattices:PanelValues:a
	NVAR b = root:Packages:Lattices:PanelValues:b
	NVAR c = root:Packages:Lattices:PanelValues:c
	NVAR alpha = root:Packages:Lattices:PanelValues:alpha
	NVAR bet = root:Packages:Lattices:PanelValues:bet
	NVAR gam = root:Packages:Lattices:PanelValues:gam
	if (NVAR_Exists(a) && NVAR_Exists(b) && NVAR_Exists(c))
		a = xtal.a	;	b = xtal.b	;	c = xtal.c
	endif
	if (NVAR_Exists(alpha) && NVAR_Exists(bet) && NVAR_Exists(gam))
		alpha = xtal.alpha	;	bet = xtal.beta	;	gam = xtal.gam
	endif
	if (NVAR_Exists(dirty))
		dirty = 1
	endif
	if (SVAR_Exists(SpaceGroupID))
		SpaceGroupID = xtal.SpaceGroupID
	endif
	if (SVAR_Exists(desc))
		desc = xtal.desc
	endif

	String wList=LatticeSetPanelList(), win
	if (ItemsInList(wList)>0)						// a #LatticePanel is displayed
		dirty = 1	;	SpaceGroupID = xtal.SpaceGroupID
		a = xtal.a	;	b = xtal.b	;	c = xtal.c
		alpha = xtal.alpha	;	bet = xtal.beta	;	gam = xtal.gam
		Variable i
		STRUCT WMSetVariableAction sva
		sva.eventCode = 2
		STRUCT WMButtonAction ba
		ba.eventCode = 2
		ba.ctrlName = "buttonLatticeSave"
		for (i=0;i<ItemsInList(wList);i+=1)	// update all occurances of the #LatticePanel
			win = StringFromList(i,wList)
			sva.win = win
			LatticePanelParamProc(sva)
			ba.win = win
			LatticePanelButtonProc(ba)
		endfor
	else													// NO #LatticePanel is displayed
		ForceLatticeToStructure(xtal)
		UpdateCrystalStructureDefaults(xtal)	// update the default xtal
	endif
	return 0
End
//
Static Function ChangeXtalSetting(xtal, id, [printIt])	// change the setting of given xtal structure 
	STRUCT crystalStructure &xtal			// this gets MODIFIED
	String id										// the requested target id, if not valid user is prompted
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : !(!printIt)

	Variable SG = xtal.SpaceGroup			// Space Group number from international tables, allowed range is [1, 230]
	String idSource = xtal.SpaceGroupID
	if (!isValidSpaceGroup(SG))
		printf "ERROR -- ChangeXtalSetting(\"%s\") xtal is invalid, nothing done.\r",id
		return 1
	elseif (numtype(str2num(id))==0 && SG != str2num(id))
		printf "ERROR -- ChangeXtalSetting(\"%s\"), Cannot convert setting \"%s\" --> \"%s\"\r",id, idSource,id
		return 1
	endif

	if (!(SG==str2num(id)))					// no id given, ask the user
		String defaultID = FindDefaultIDforSG(xtal.SpaceGroup)
		String allIDs=MakeAllIDs()
		String idList="", idNumList = symmtry2SG(num2str(SG)+"*", types=16)	// find the Space Group number from the symmetry string
		if (ItemsInList(idNumList)<1)
			if (printIt)
				print "Only one setting for this Space Group, nothing to do"
			endif
			return 0
		endif
		Variable i, idNum, def
		for (i=0;i<ItemsInList(idNumList);i+=1)
			idNum = str2num(StringFromList(i,idNumList))
			id = StringFromList(idNum-1,allIDs)
			def = cmpstr(id, defaultID) == 0
			if (cmpstr(id, idSource) == 0)
				continue								// skip the current setting
			endif
			id += ",  " + getHMsym2(idNum)
			if (def)
				id += "   --Default--"
			endif
			idList += id+";"
		endfor
		id = defaultID								// preset to default id
		String str
		sprintf str, "Target Space Group ID  (currently it is '%s')",idSource
		Prompt id, str, popup, idList
		DoPrompt "Target Space Group", id
		if (V_flag)
			return 1
		endif
		id = StringFromList(0,id,",")
		printIt = 1
	endif

	if (printIt)
		printf "ChangeXtalSetting(\"%s\")\r", id
	endif

	if (numtype(str2num(id))==0 && SG != str2num(id))
		printf "ERROR -- ChangeXtalSetting(\"%s\"), Cannot convert setting \"%s\" --> \"%s\"\r",id, idSource,id
		return 1
	endif

	if (printIt)
		print "Before:"
		print_crystalStructure(xtal, brief=1)
	endif

	ConvertSetting(xtal, id)					// change setting of xtal from: idSource --> id

	if (printIt)
		print "\r  After:"
		print_crystalStructure(xtal, brief=1)
	endif
	return 0
End
//
Static Function ConvertSetting(xtal, target)	// change the setting of given xtal structure, lattice & all atoms
	STRUCT crystalStructure &xtal		// this is CHANGED from the SpaceGroupID settings to the setting of "target"
	String target								// the new desired setting, if target="", then convert to default setting
	if (!isValidLatticeConstants(xtal))
		return 1											// given xtal is invalid	
	endif

	String source = xtal.SpaceGroupID			// current SpaceGroupID
	String defalt = FindDefaultIDforSG(xtal.SpaceGroup)	// default SpaceGroupID
	if (!isValidSpaceGroupID(source))			// if source is invalid, do nothing
		return 1
	endif
	if (!isValidSpaceGroupID(target))			// if target is invalid, set to the default setting for this space group
		target = defalt
	endif
	if (StringMatch(source,target))				// converting to itself, nothing to do
		return 0
	endif

	xtal.SpaceGroupID = target
	xtal.SpaceGroupIDnum = SpaceGroupID2num(target)
	Wave CBM = GetSettingTransForm(target)	// converts Defalt --> Target
	Wave CBM0 = GetSettingTransForm(source)	// converts Defalt --> Source
	Wave DL = directFrom_xtal(xtal)				// returns a FREE wave with real lattice
	MatrixOp/FREE GS = DL^t x DL					// source metrical matrix, see description below
	MatrixOp/FREE GD = Inv(CBM0^t) x GS x Inv(CBM0)		//	1st convert source --> default,  GD is default metrical matrix
	MatrixOp/FREE GT = CBM^t x GD x CBM		//	2nd convert default --> target,  target metrical matrix

	// calculate the new lattice constants from target metrical matrix
	Variable a=sqrt(GT[0][0]), b=sqrt(GT[1][1]), c=sqrt(GT[2][2])
	xtal.a = a	;	xtal.b = b	;	xtal.c = c
	xtal.alpha = acos( GT[1][2]/(b*c) ) * 180/PI
	xtal.beta  = acos( GT[0][2]/(a*c) ) * 180/PI
	xtal.gam   = acos( GT[0][1]/(a*b) ) * 180/PI
	setDirectRecip(xtal)							// update the other xtal values
	//
	//	The formula for the transformation of the metrical matrix G(A) of Setting A to the metrical matrix G(B) of Setting B is (see e.g. Boisen & Gibbs [2]):
	//			G(B) = transpose(InvCBMx) * G(A) * InvCBMx
	//
	//										a¥a		a¥b		a¥c
	//	metrical matrix = 			a¥b		b¥b		b¥c
	//										a¥c		b¥c		c¥c
	//
	//										a^2			a*b*cos(gamma)		a*c*cos(beta)
	//	metrical matrix = 	a*b*cos(gamma)				b^2			b*c*cos(alpha)
	//								a*c*cos(beta)		b*c*cos(alpha)			c^2
	//
	//		if DL = {a,b,c}, where a,b,c are all column vectors
	//		metrical matrix = G = DL^t x DL

	// transform the atomic coordinates:
	Make/N=3/D/FREE fracS
	Variable i, Natom=xtal.N
	for (i=0; i<Natom; i+=1)
		fracS = {xtal.atom[i].x, xtal.atom[i].y, xtal.atom[i].z}
		//	MatrixOp/FREE fracD = CBM0 x fracS				// converts source --> default
		//	MatrixOp/FREE fracT = Inv(CBM) x fracD		// converts default --> target
		MatrixOp/FREE fracT = Inv(CBM) x CBM0 x fracS	// convert: source --> default --> target
		MatrixOp/FREE fracT = fracT - floor(fracT)		// reduce to first cell
		fracT = fracT<1e-12 ? 0 : mod(fracT,1)			// a fractional coord < 1e-12 is 0
		xtal.atom[i].x = fracT[0]
		xtal.atom[i].y = fracT[1]
		xtal.atom[i].z = fracT[2]
	endfor

	// try to fix name when changing between Hexagonal and Rhombohedral
	if (strsearch(source,":H",0)>0 && strsearch(target,":R",0)>0)		// Hexagonal --> Rhombohedral
		xtal.desc = ReplaceString("Hexagonal", xtal.desc,"Rhombohedral")
		xtal.desc = ReplaceString("Hex", xtal.desc,"Rhom")
	elseif (strsearch(source,":R",0)>0 && strsearch(target,":H",0)>0)	// Rhombohedral --> Hexagonal
		xtal.desc = ReplaceString("Rhombohedral", xtal.desc,"Hexagonal")
		xtal.desc = ReplaceString("Rhom", xtal.desc,"Hex")
	endif

	return 0
End












Static Function/WAVE GetWyckoffSymStrings(SpaceGroupID)
	String SpaceGroupID
	if (!isValidSpaceGroupID(SpaceGroupID))
		return $""
	endif

	Variable SG_idNum = SpaceGroupID2num(SpaceGroupID)
	Make/N=(530)/T/FREE WyckoffSyms
	// Triclinic SG[1,2]  SG_idNum [1-2]   (2 idNums)
	WyckoffSyms[0]   = "a:x,y,z:1;"
	WyckoffSyms[1]   = "a:0,0,0:1;b:0,0,1/2:1;c:0,1/2,0:1;d:1/2,0,0:1;e:1/2,1/2,0:1;f:1/2,0,1/2:1;g:0,1/2,1/2:1;h:1/2,1/2,1/2:1;i:x,y,z:2;"
	// Monoclinic SG[3,15]  SG_idNum [3,107]   (105 idNums)
	WyckoffSyms[2]   = "a:0,y,0:1;b:0,y,1/2:1;c:1/2,y,0:1;d:1/2,y,1/2:1;e:x,y,z:2;"
	WyckoffSyms[3]   = "a:0,0,y:1;b:1/2,0,y:1;c:0,1/2,y:1;d:1/2,1/2,y:1;e:z,x,y:2;"
	WyckoffSyms[4]   = "a:y,0,0:1;b:y,1/2,0:1;c:y,0,1/2:1;d:y,1/2,1/2:1;e:y,z,x:2;"
	WyckoffSyms[5]   = "a:x,y,z:2;"
	WyckoffSyms[6]   = "a:z,x,y:2;"
	WyckoffSyms[7]   = "a:y,z,x:2;"
	WyckoffSyms[8]   = "a:0,y,0:2;b:0,y,1/2:2;c:x,y,z:4;"
	WyckoffSyms[9]   = "a:0,y,0:2;b:1/2,y,1/2:2;c:-z,y,x-z:4;"
	WyckoffSyms[10]  = "a:0,y,0:2;b:1/2,y,0:2;c:-x+z,y,-x:4;"
	WyckoffSyms[11]  = "a:0,0,y:2;b:1/2,0,y:2;c:z,x,y:4;"
	WyckoffSyms[12]  = "a:0,0,y:2;b:1/2,1/2,y:2;c:x-z,-z,y:4;"
	WyckoffSyms[13]  = "a:0,0,y:2;b:0,1/2,y:2;c:-x,-x+z,y:4;"
	WyckoffSyms[14]  = "a:y,0,0:2;b:y,1/2,0:2;c:y,z,x:4;"
	WyckoffSyms[15]  = "a:y,0,0:2;b:y,1/2,1/2:2;c:y,x-z,-z:4;"
	WyckoffSyms[16]  = "a:y,0,0:2;b:y,0,1/2:2;c:y,-x,-x+z:4;"
	WyckoffSyms[17]  = "a:x,0,z:1;b:x,1/2,z:1;c:x,y,z:2;"
	WyckoffSyms[18]  = "a:z,x,0:1;b:z,x,1/2:1;c:z,x,y:2;"
	WyckoffSyms[19]  = "a:0,z,x:1;b:1/2,z,x:1;c:y,z,x:2;"
	WyckoffSyms[20]  = "a:x,y,z:2;"
	WyckoffSyms[21]  = "a:-z,y,x-z:2;"
	WyckoffSyms[22]  = "a:-x+z,y,-x:2;"
	WyckoffSyms[23]  = "a:z,x,y:2;"
	WyckoffSyms[24]  = "a:x-z,-z,y:2;"
	WyckoffSyms[25]  = "a:-x,-x+z,y:2;"
	WyckoffSyms[26]  = "a:y,z,x:2;"
	WyckoffSyms[27]  = "a:y,x-z,-z:2;"
	WyckoffSyms[28]  = "a:y,-x,-x+z:2;"
	WyckoffSyms[29]  = "a:x,0,z:2;b:x,y,z:4;"
	WyckoffSyms[30]  = "a:-z,0,x-z:2;b:-z,y,x-z:4;"
	WyckoffSyms[31]  = "a:-x+z,0,-x:2;b:-x+z,y,-x:4;"
	WyckoffSyms[32]  = "a:z,x,0:2;b:z,x,y:4;"
	WyckoffSyms[33]  = "a:x-z,-z,0:2;b:x-z,-z,y:4;"
	WyckoffSyms[34]  = "a:-x,-x+z,0:2;b:-x,-x+z,y:4;"
	WyckoffSyms[35]  = "a:0,z,x:2;b:y,z,x:4;"
	WyckoffSyms[36]  = "a:0,x-z,-z:2;b:y,x-z,-z:4;"
	WyckoffSyms[37]  = "a:0,-x,-x+z:2;b:y,-x,-x+z:4;"
	WyckoffSyms[38]  = "a:x,y,z:4;"
	WyckoffSyms[39]  = "a:-z,y,x-z:4;"
	WyckoffSyms[40]  = "a:-x+z,y,-x:4;"
	WyckoffSyms[41]  = "a:-z,y,x-z:4;"
	WyckoffSyms[42]  = "a:x,y,z:4;"
	WyckoffSyms[43]  = "a:-x+z,y,-x:4;"
	WyckoffSyms[44]  = "a:z,x,y:4;"
	WyckoffSyms[45]  = "a:x-z,-z,y:4;"
	WyckoffSyms[46]  = "a:-x,-x+z,y:4;"
	WyckoffSyms[47]  = "a:x-z,-z,y:4;"
	WyckoffSyms[48]  = "a:z,x,y:4;"
	WyckoffSyms[49]  = "a:-x,-x+z,y:4;"
	WyckoffSyms[50]  = "a:y,z,x:4;"
	WyckoffSyms[51]  = "a:y,x-z,-z:4;"
	WyckoffSyms[52]  = "a:y,-x,-x+z:4;"
	WyckoffSyms[53]  = "a:y,x-z,-z:4;"
	WyckoffSyms[54]  = "a:y,z,x:4;"
	WyckoffSyms[55]  = "a:y,-x,-x+z:4;"
	WyckoffSyms[56]  = "a:0,0,0:1;b:0,1/2,0:1;c:0,0,1/2:1;d:1/2,0,0:1;e:1/2,1/2,0:1;f:0,1/2,1/2:1;g:1/2,0,1/2:1;h:1/2,1/2,1/2:1;i:0,y,0:2;j:1/2,y,0:2;"
	WyckoffSyms[56] += "k:0,y,1/2:2;l:1/2,y,1/2:2;m:x,0,z:2;n:x,1/2,z:2;o:x,y,z:4;"
	WyckoffSyms[57]  = "a:0,0,0:1;b:0,0,1/2:1;c:1/2,0,0:1;d:0,1/2,0:1;e:0,1/2,1/2:1;f:1/2,0,1/2:1;g:1/2,1/2,0:1;h:1/2,1/2,1/2:1;i:0,0,y:2;j:0,1/2,y:2;"
	WyckoffSyms[57] += "k:1/2,0,y:2;l:1/2,1/2,y:2;m:z,x,0:2;n:z,x,1/2:2;o:z,x,y:4;"
	WyckoffSyms[58]  = "a:0,0,0:1;b:1/2,0,0:1;c:0,1/2,0:1;d:0,0,1/2:1;e:1/2,0,1/2:1;f:1/2,1/2,0:1;g:0,1/2,1/2:1;h:1/2,1/2,1/2:1;i:y,0,0:2;j:y,0,1/2:2;"
	WyckoffSyms[58] += "k:y,1/2,0:2;l:y,1/2,1/2:2;m:0,z,x:2;n:1/2,z,x:2;o:y,z,x:4;"
	WyckoffSyms[59]  = "a:0,0,0:2;b:1/2,0,0:2;c:0,0,1/2:2;d:1/2,0,1/2:2;e:x,1/4,z:2;f:x,y,z:4;"
	WyckoffSyms[60]  = "a:0,0,0:2;b:0,1/2,0:2;c:1/2,0,0:2;d:1/2,1/2,0:2;e:z,x,1/4:2;f:z,x,y:4;"
	WyckoffSyms[61]  = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,0:2;d:0,1/2,1/2:2;e:1/4,z,x:2;f:y,z,x:4;"
	WyckoffSyms[62]  = "a:0,0,0:2;b:0,1/2,0:2;c:0,0,1/2:2;d:0,1/2,1/2:2;e:1/4,1/4,0:4;f:1/4,1/4,1/2:4;g:0,y,0:4;h:0,y,1/2:4;i:x,0,z:4;j:x,y,z:8;"
	WyckoffSyms[63]  = "a:0,0,0:2;b:0,1/2,0:2;c:1/2,0,1/2:2;d:1/2,1/2,1/2:2;e:0,1/4,1/4:4;f:1/2,1/4,1/4:4;g:0,y,0:4;h:1/2,y,1/2:4;i:-z,0,x-z:4;"
	WyckoffSyms[63] += "j:-z,y,x-z:8;"
	WyckoffSyms[64]  = "a:0,0,0:2;b:0,1/2,0:2;c:1/2,0,0:2;d:1/2,1/2,0:2;e:1/4,1/4,1/4:4;f:1/4,1/4,1/4:4;g:0,y,0:4;h:1/2,y,0:4;i:-x+z,0,-x:4;"
	WyckoffSyms[64] += "j:-x+z,y,-x:8;"
	WyckoffSyms[65]  = "a:0,0,0:2;b:0,0,1/2:2;c:1/2,0,0:2;d:1/2,0,1/2:2;e:0,1/4,1/4:4;f:1/2,1/4,1/4:4;g:0,0,y:4;h:1/2,0,y:4;i:z,x,0:4;j:z,x,y:8;"
	WyckoffSyms[66]  = "a:0,0,0:2;b:0,0,1/2:2;c:1/2,1/2,0:2;d:1/2,1/2,1/2:2;e:1/4,0,1/4:4;f:1/4,1/2,1/4:4;g:0,0,y:4;h:1/2,1/2,y:4;i:x-z,-z,0:4;"
	WyckoffSyms[66] += "j:x-z,-z,y:8;"
	WyckoffSyms[67]  = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,0:2;d:0,1/2,1/2:2;e:1/4,1/4,1/4:4;f:1/4,1/4,1/4:4;g:0,0,y:4;h:0,1/2,y:4;i:-x,-x+z,0:4;"
	WyckoffSyms[67] += "j:-x,-x+z,y:8;"
	WyckoffSyms[68]  = "a:0,0,0:2;b:1/2,0,0:2;c:0,1/2,0:2;d:1/2,1/2,0:2;e:1/4,0,1/4:4;f:1/4,1/2,1/4:4;g:y,0,0:4;h:y,1/2,0:4;i:0,z,x:4;j:y,z,x:8;"
	WyckoffSyms[69]  = "a:0,0,0:2;b:1/2,0,0:2;c:0,1/2,1/2:2;d:1/2,1/2,1/2:2;e:1/4,1/4,0:4;f:1/4,1/4,1/2:4;g:y,0,0:4;h:y,1/2,1/2:4;i:0,x-z,-z:4;"
	WyckoffSyms[69] += "j:y,x-z,-z:8;"
	WyckoffSyms[70]  = "a:0,0,0:2;b:1/2,0,0:2;c:0,0,1/2:2;d:1/2,0,1/2:2;e:1/4,1/4,1/4:4;f:1/4,1/4,1/4:4;g:y,0,0:4;h:y,0,1/2:4;i:0,-x,-x+z:4;"
	WyckoffSyms[70] += "j:y,-x,-x+z:8;"
	WyckoffSyms[71]  = "a:0,0,0:2;b:1/2,1/2,0:2;c:0,1/2,0:2;d:1/2,0,0:2;e:0,y,1/4:2;f:1/2,y,1/4:2;g:x,y,z:4;"
	WyckoffSyms[72]  = "a:0,0,0:2;b:0,1/2,1/2:2;c:0,1/2,0:2;d:0,0,1/2:2;e:1/4,y,1/4:2;f:1/4,y,1/4:2;g:-z,y,x-z:4;"
	WyckoffSyms[73]  = "a:0,0,0:2;b:1/2,1/2,1/2:2;c:0,1/2,0:2;d:1/2,0,1/2:2;e:1/4,y,0:2;f:1/4,y,1/2:2;g:-x+z,y,-x:4;"
	WyckoffSyms[74]  = "a:0,0,0:2;b:0,1/2,1/2:2;c:0,0,1/2:2;d:0,1/2,0:2;e:1/4,0,y:2;f:1/4,1/2,y:2;g:z,x,y:4;"
	WyckoffSyms[75]  = "a:0,0,0:2;b:1/2,0,1/2:2;c:0,0,1/2:2;d:1/2,0,0:2;e:1/4,1/4,y:2;f:1/4,1/4,y:2;g:x-z,-z,y:4;"
	WyckoffSyms[76]  = "a:0,0,0:2;b:1/2,1/2,1/2:2;c:0,0,1/2:2;d:1/2,1/2,0:2;e:0,1/4,y:2;f:1/2,1/4,y:2;g:-x,-x+z,y:4;"
	WyckoffSyms[77]  = "a:0,0,0:2;b:1/2,0,1/2:2;c:1/2,0,0:2;d:0,0,1/2:2;e:y,1/4,0:2;f:y,1/4,1/2:2;g:y,z,x:4;"
	WyckoffSyms[78]  = "a:0,0,0:2;b:1/2,1/2,0:2;c:1/2,0,0:2;d:0,1/2,0:2;e:y,1/4,1/4:2;f:y,1/4,1/4:2;g:y,x-z,-z:4;"
	WyckoffSyms[79]  = "a:0,0,0:2;b:1/2,1/2,1/2:2;c:1/2,0,0:2;d:0,1/2,1/2:2;e:y,0,1/4:2;f:y,1/2,1/4:2;g:y,-x,-x+z:4;"
	WyckoffSyms[80]  = "a:0,0,0:2;b:1/2,0,0:2;c:0,0,1/2:2;d:1/2,0,1/2:2;e:x,y,z:4;"
	WyckoffSyms[81]  = "a:0,0,0:2;b:0,0,1/2:2;c:1/2,0,1/2:2;d:1/2,0,0:2;e:-z,y,x-z:4;"
	WyckoffSyms[82]  = "a:0,0,0:2;b:1/2,0,1/2:2;c:1/2,0,0:2;d:0,0,1/2:2;e:-x+z,y,-x:4;"
	WyckoffSyms[83]  = "a:0,0,0:2;b:0,1/2,0:2;c:1/2,0,0:2;d:1/2,1/2,0:2;e:z,x,y:4;"
	WyckoffSyms[84]  = "a:0,0,0:2;b:1/2,0,0:2;c:1/2,1/2,0:2;d:0,1/2,0:2;e:x-z,-z,y:4;"
	WyckoffSyms[85]  = "a:0,0,0:2;b:1/2,1/2,0:2;c:0,1/2,0:2;d:1/2,0,0:2;e:-x,-x+z,y:4;"
	WyckoffSyms[86]  = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,0:2;d:0,1/2,1/2:2;e:y,z,x:4;"
	WyckoffSyms[87]  = "a:0,0,0:2;b:0,1/2,0:2;c:0,1/2,1/2:2;d:0,0,1/2:2;e:y,x-z,-z:4;"
	WyckoffSyms[88]  = "a:0,0,0:2;b:0,1/2,1/2:2;c:0,0,1/2:2;d:0,1/2,0:2;e:y,-x,-x+z:4;"
	WyckoffSyms[89]  = "a:0,0,0:4;b:0,1/2,0:4;c:1/4,1/4,0:4;d:1/4,1/4,1/2:4;e:0,y,1/4:4;f:x,y,z:8;"
	WyckoffSyms[90]  = "a:0,0,0:4;b:0,1/2,0:4;c:0,1/4,1/4:4;d:1/2,1/4,1/4:4;e:1/4,y,1/4:4;f:-z,y,x-z:8;"
	WyckoffSyms[91]  = "a:0,0,0:4;b:0,1/2,0:4;c:1/4,1/4,1/4:4;d:1/4,1/4,1/4:4;e:1/4,y,0:4;f:-x+z,y,-x:8;"
	WyckoffSyms[92]  = "a:0,0,0:4;b:0,1/2,0:4;c:0,1/4,1/4:4;d:1/2,1/4,1/4:4;e:1/4,y,1/4:4;f:-z,y,x-z:8;"
	WyckoffSyms[93]  = "a:0,0,0:4;b:0,1/2,0:4;c:1/4,1/4,0:4;d:1/4,1/4,1/2:4;e:0,y,1/4:4;f:x,y,z:8;"
	WyckoffSyms[94]  = "a:0,0,0:4;b:0,1/2,0:4;c:1/4,1/4,1/4:4;d:1/4,1/4,1/4:4;e:1/4,y,0:4;f:-x+z,y,-x:8;"
	WyckoffSyms[95]  = "a:0,0,0:4;b:0,0,1/2:4;c:0,1/4,1/4:4;d:1/2,1/4,1/4:4;e:1/4,0,y:4;f:z,x,y:8;"
	WyckoffSyms[96]  = "a:0,0,0:4;b:0,0,1/2:4;c:1/4,0,1/4:4;d:1/4,1/2,1/4:4;e:1/4,1/4,y:4;f:x-z,-z,y:8;"
	WyckoffSyms[97]  = "a:0,0,0:4;b:0,0,1/2:4;c:1/4,1/4,1/4:4;d:1/4,1/4,1/4:4;e:0,1/4,y:4;f:-x,-x+z,y:8;"
	WyckoffSyms[98]  = "a:0,0,0:4;b:0,0,1/2:4;c:1/4,0,1/4:4;d:1/4,1/2,1/4:4;e:1/4,1/4,y:4;f:x-z,-z,y:8;"
	WyckoffSyms[99]  = "a:0,0,0:4;b:0,0,1/2:4;c:0,1/4,1/4:4;d:1/2,1/4,1/4:4;e:1/4,0,y:4;f:z,x,y:8;"
	WyckoffSyms[100]  = "a:0,0,0:4;b:0,0,1/2:4;c:1/4,1/4,1/4:4;d:1/4,1/4,1/4:4;e:0,1/4,y:4;f:-x,-x+z,y:8;"
	WyckoffSyms[101]  = "a:0,0,0:4;b:1/2,0,0:4;c:1/4,0,1/4:4;d:1/4,1/2,1/4:4;e:y,1/4,0:4;f:y,z,x:8;"
	WyckoffSyms[102]  = "a:0,0,0:4;b:1/2,0,0:4;c:1/4,1/4,0:4;d:1/4,1/4,1/2:4;e:y,1/4,1/4:4;f:y,x-z,-z:8;"
	WyckoffSyms[103]  = "a:0,0,0:4;b:1/2,0,0:4;c:1/4,1/4,1/4:4;d:1/4,1/4,1/4:4;e:y,0,1/4:4;f:y,-x,-x+z:8;"
	WyckoffSyms[104]  = "a:0,0,0:4;b:1/2,0,0:4;c:1/4,1/4,0:4;d:1/4,1/4,1/2:4;e:y,1/4,1/4:4;f:y,x-z,-z:8;"
	WyckoffSyms[105]  = "a:0,0,0:4;b:1/2,0,0:4;c:1/4,0,1/4:4;d:1/4,1/2,1/4:4;e:y,1/4,0:4;f:y,z,x:8;"
	WyckoffSyms[106]  = "a:0,0,0:4;b:1/2,0,0:4;c:1/4,1/4,1/4:4;d:1/4,1/4,1/4:4;e:y,0,1/4:4;f:y,-x,-x+z:8;"
	WyckoffSyms[107]  = "a:0,0,0:1;b:1/2,0,0:1;c:0,1/2,0:1;d:0,0,1/2:1;e:1/2,1/2,0:1;f:1/2,0,1/2:1;g:0,1/2,1/2:1;h:1/2,1/2,1/2:1;i:x,0,0:2;j:x,0,1/2:2;"
	WyckoffSyms[107] += "k:x,1/2,0:2;l:x,1/2,1/2:2;m:0,y,0:2;n:0,y,1/2:2;o:1/2,y,0:2;p:1/2,y,1/2:2;q:0,0,z:2;r:1/2,0,z:2;s:0,1/2,z:2;t:1/2,1/2,z:2;"
	// Orthorhombic SG[16,74]  SG_idNum [108,348]   (241 idNums)
	WyckoffSyms[107] += "u:x,y,z:4;"
	WyckoffSyms[108]  = "a:x,0,0:2;b:x,1/2,0:2;c:0,y,1/4:2;d:1/2,y,1/4:2;e:x,y,z:4;"
	WyckoffSyms[109]  = "a:0,x,0:2;b:0,x,1/2:2;c:1/4,0,y:2;d:1/4,1/2,y:2;e:z,x,y:4;"
	WyckoffSyms[110]  = "a:0,0,x:2;b:1/2,0,x:2;c:y,1/4,0:2;d:y,1/4,1/2:2;e:y,z,x:4;"
	WyckoffSyms[111]  = "a:0,0,z:2;b:0,1/2,z:2;c:x,y,z:4;"
	WyckoffSyms[112]  = "a:z,0,0:2;b:z,0,1/2:2;c:z,x,y:4;"
	WyckoffSyms[113]  = "a:0,z,0:2;b:1/2,z,0:2;c:y,z,x:4;"
	WyckoffSyms[114]  = "a:x,y,z:4;"
	WyckoffSyms[115]  = "a:x,0,0:4;b:0,y,1/4:4;c:x,y,z:8;"
	WyckoffSyms[116]  = "a:0,x,0:4;b:1/4,0,y:4;c:z,x,y:8;"
	WyckoffSyms[117]  = "a:0,0,x:4;b:y,1/4,0:4;c:y,z,x:8;"
	WyckoffSyms[118]  = "a:0,0,0:2;b:0,1/2,0:2;c:1/2,0,1/2:2;d:0,0,1/2:2;e:x,0,0:4;f:x,0,1/2:4;g:0,y,0:4;h:0,y,1/2:4;i:0,0,z:4;j:0,1/2,z:4;k:1/4,1/4,z:4;"
	WyckoffSyms[118] += "l:x,y,z:8;"
	WyckoffSyms[119]  = "a:0,0,0:2;b:0,0,1/2:2;c:1/2,1/2,0:2;d:1/2,0,0:2;e:0,x,0:4;f:1/2,x,0:4;g:0,0,y:4;h:1/2,0,y:4;i:z,0,0:4;j:z,0,1/2:4;k:z,1/4,1/4:4;"
	WyckoffSyms[119] += "l:z,x,y:8;"
	WyckoffSyms[120]  = "a:0,0,0:2;b:1/2,0,0:2;c:0,1/2,1/2:2;d:0,1/2,0:2;e:0,0,x:4;f:0,1/2,x:4;g:y,0,0:4;h:y,1/2,0:4;i:0,z,0:4;j:1/2,z,0:4;k:1/4,z,1/4:4;"
	WyckoffSyms[120] += "l:y,z,x:8;"
	WyckoffSyms[121]  = "a:0,0,0:4;b:0,0,1/2:4;c:1/4,1/4,1/4:4;d:1/4,1/4,3/4:4;e:x,0,0:8;f:0,y,0:8;g:0,0,z:8;h:1/4,1/4,z:8;i:1/4,y,1/4:8;j:x,1/4,1/4:8;"
	WyckoffSyms[121] += "k:x,y,z:16;"
	WyckoffSyms[122]  = "a:0,0,0:2;b:1/2,0,0:2;c:0,0,1/2:2;d:0,1/2,0:2;e:x,0,0:4;f:x,0,1/2:4;g:0,y,0:4;h:1/2,y,0:4;i:0,0,z:4;j:0,1/2,z:4;k:x,y,z:8;"
	WyckoffSyms[123]  = "a:x,0,1/4:4;b:1/4,y,0:4;c:0,1/4,z:4;d:x,y,z:8;"
	WyckoffSyms[124]  = "a:0,0,z:1;b:0,1/2,z:1;c:1/2,0,z:1;d:1/2,1/2,z:1;e:x,0,z:2;f:x,1/2,z:2;g:0,y,z:2;h:1/2,y,z:2;i:x,y,z:4;"
	WyckoffSyms[125]  = "a:z,0,0:1;b:z,0,1/2:1;c:z,1/2,0:1;d:z,1/2,1/2:1;e:z,x,0:2;f:z,x,1/2:2;g:z,0,y:2;h:z,1/2,y:2;i:z,x,y:4;"
	WyckoffSyms[126]  = "a:0,z,0:1;b:1/2,z,0:1;c:0,z,1/2:1;d:1/2,z,1/2:1;e:0,z,x:2;f:1/2,z,x:2;g:y,z,0:2;h:y,z,1/2:2;i:y,z,x:4;"
	WyckoffSyms[127]  = "a:0,y,z:2;b:1/2,y,z:2;c:x,y,z:4;"
	WyckoffSyms[128]  = "a:y,0,-z:2;b:y,1/2,-z:2;c:y,x,-z:4;"
	WyckoffSyms[129]  = "a:z,0,y:2;b:z,1/2,y:2;c:z,x,y:4;"
	WyckoffSyms[130]  = "a:-z,y,0:2;b:-z,y,1/2:2;c:-z,y,x:4;"
	WyckoffSyms[131]  = "a:y,z,0:2;b:y,z,1/2:2;c:y,z,x:4;"
	WyckoffSyms[132]  = "a:0,-z,y:2;b:1/2,-z,y:2;c:x,-z,y:4;"
	WyckoffSyms[133]  = "a:0,0,z:2;b:0,1/2,z:2;c:1/2,0,z:2;d:1/2,1/2,z:2;e:x,y,z:4;"
	WyckoffSyms[134]  = "a:z,0,0:2;b:z,0,1/2:2;c:z,1/2,0:2;d:z,1/2,1/2:2;e:z,x,y:4;"
	WyckoffSyms[135]  = "a:0,z,0:2;b:1/2,z,0:2;c:0,z,1/2:2;d:1/2,z,1/2:2;e:y,z,x:4;"
	WyckoffSyms[136]  = "a:0,0,z:2;b:0,1/2,z:2;c:1/4,y,z:2;d:x,y,z:4;"
	WyckoffSyms[137]  = "a:0,0,-z:2;b:1/2,0,-z:2;c:y,1/4,-z:2;d:y,x,-z:4;"
	WyckoffSyms[138]  = "a:z,0,0:2;b:z,0,1/2:2;c:z,1/4,y:2;d:z,x,y:4;"
	WyckoffSyms[139]  = "a:-z,0,0:2;b:-z,1/2,0:2;c:-z,y,1/4:2;d:-z,y,x:4;"
	WyckoffSyms[140]  = "a:0,z,0:2;b:1/2,z,0:2;c:y,z,1/4:2;d:y,z,x:4;"
	WyckoffSyms[141]  = "a:0,-z,0:2;b:0,-z,1/2:2;c:1/4,-z,y:2;d:x,-z,y:4;"
	WyckoffSyms[142]  = "a:x,y,z:4;"
	WyckoffSyms[143]  = "a:y,x,-z:4;"
	WyckoffSyms[144]  = "a:z,x,y:4;"
	WyckoffSyms[145]  = "a:-z,y,x:4;"
	WyckoffSyms[146]  = "a:y,z,x:4;"
	WyckoffSyms[147]  = "a:x,-z,y:4;"
	WyckoffSyms[148]  = "a:0,0,z:2;b:1/2,0,z:2;c:x,y,z:4;"
	WyckoffSyms[149]  = "a:0,0,-z:2;b:0,1/2,-z:2;c:y,x,-z:4;"
	WyckoffSyms[150]  = "a:z,0,0:2;b:z,1/2,0:2;c:z,x,y:4;"
	WyckoffSyms[151]  = "a:-z,0,0:2;b:-z,0,1/2:2;c:-z,y,x:4;"
	WyckoffSyms[152]  = "a:0,z,0:2;b:0,z,1/2:2;c:y,z,x:4;"
	WyckoffSyms[153]  = "a:0,-z,0:2;b:1/2,-z,0:2;c:x,-z,y:4;"
	WyckoffSyms[154]  = "a:0,y,z:2;b:x,y,z:4;"
	WyckoffSyms[155]  = "a:y,0,-z:2;b:y,x,-z:4;"
	WyckoffSyms[156]  = "a:z,0,y:2;b:z,x,y:4;"
	WyckoffSyms[157]  = "a:-z,y,0:2;b:-z,y,x:4;"
	WyckoffSyms[158]  = "a:y,z,0:2;b:y,z,x:4;"
	WyckoffSyms[159]  = "a:0,-z,y:2;b:x,-z,y:4;"
	WyckoffSyms[160]  = "a:0,0,z:2;b:0,1/2,z:2;c:x,y,z:4;"
	WyckoffSyms[161]  = "a:z,0,0:2;b:z,0,1/2:2;c:z,x,y:4;"
	WyckoffSyms[162]  = "a:0,z,0:2;b:1/2,z,0:2;c:y,z,x:4;"
	WyckoffSyms[163]  = "a:x,y,z:4;"
	WyckoffSyms[164]  = "a:y,x,-z:4;"
	WyckoffSyms[165]  = "a:z,x,y:4;"
	WyckoffSyms[166]  = "a:-z,y,x:4;"
	WyckoffSyms[167]  = "a:y,z,x:4;"
	WyckoffSyms[168]  = "a:x,-z,y:4;"
	WyckoffSyms[169]  = "a:0,0,z:2;b:0,1/2,z:2;c:x,y,z:4;"
	WyckoffSyms[170]  = "a:z,0,0:2;b:z,0,1/2:2;c:z,x,y:4;"
	WyckoffSyms[171]  = "a:0,z,0:2;b:1/2,z,0:2;c:y,z,x:4;"
	WyckoffSyms[172]  = "a:0,0,z:2;b:0,1/2,z:2;c:1/4,1/4,z:4;d:x,0,z:4;e:0,y,z:4;f:x,y,z:8;"
	WyckoffSyms[173]  = "a:z,0,0:2;b:z,0,1/2:2;c:z,1/4,1/4:4;d:z,x,0:4;e:z,0,y:4;f:z,x,y:8;"
	WyckoffSyms[174]  = "a:0,z,0:2;b:1/2,z,0:2;c:1/4,z,1/4:4;d:0,z,x:4;e:y,z,0:4;f:y,z,x:8;"
	WyckoffSyms[175]  = "a:0,y,z:4;b:x,y,z:8;"
	WyckoffSyms[176]  = "a:y,0,-z:4;b:y,x,-z:8;"
	WyckoffSyms[177]  = "a:z,0,y:4;b:z,x,y:8;"
	WyckoffSyms[178]  = "a:-z,y,0:4;b:-z,y,x:8;"
	WyckoffSyms[179]  = "a:y,z,0:4;b:y,z,x:8;"
	WyckoffSyms[180]  = "a:0,-z,y:4;b:x,-z,y:8;"
	WyckoffSyms[181]  = "a:0,0,z:4;b:0,1/2,z:4;c:1/4,1/4,z:4;d:x,y,z:8;"
	WyckoffSyms[182]  = "a:z,0,0:4;b:z,0,1/2:4;c:z,1/4,1/4:4;d:z,x,y:8;"
	WyckoffSyms[183]  = "a:0,z,0:4;b:1/2,z,0:4;c:1/4,z,1/4:4;d:y,z,x:8;"
	WyckoffSyms[184]  = "a:0,0,z:2;b:1/2,0,z:2;c:x,0,z:4;d:0,y,z:4;e:1/2,y,z:4;f:x,y,z:8;"
	WyckoffSyms[185]  = "a:0,0,-z:2;b:0,1/2,-z:2;c:0,x,-z:4;d:y,0,-z:4;e:y,1/2,-z:4;f:y,x,-z:8;"
	WyckoffSyms[186]  = "a:z,0,0:2;b:z,1/2,0:2;c:z,x,0:4;d:z,0,y:4;e:z,1/2,y:4;f:z,x,y:8;"
	WyckoffSyms[187]  = "a:-z,0,0:2;b:-z,0,1/2:2;c:-z,0,x:4;d:-z,y,0:4;e:-z,y,1/2:4;f:-z,y,x:8;"
	WyckoffSyms[188]  = "a:0,z,0:2;b:0,z,1/2:2;c:0,z,x:4;d:y,z,0:4;e:y,z,1/2:4;f:y,z,x:8;"
	WyckoffSyms[189]  = "a:0,-z,0:2;b:1/2,-z,0:2;c:x,-z,0:4;d:0,-z,y:4;e:1/2,-z,y:4;f:x,-z,y:8;"
	WyckoffSyms[190]  = "a:0,0,z:4;b:1/2,0,z:4;c:x,1/4,z:4;d:x,y,z:8;"
	WyckoffSyms[191]  = "a:0,0,-z:4;b:0,1/2,-z:4;c:1/4,x,-z:4;d:y,x,-z:8;"
	WyckoffSyms[192]  = "a:z,0,0:4;b:z,1/2,0:4;c:z,x,1/4:4;d:z,x,y:8;"
	WyckoffSyms[193]  = "a:-z,0,0:4;b:-z,0,1/2:4;c:-z,1/4,x:4;d:-z,y,x:8;"
	WyckoffSyms[194]  = "a:0,z,0:4;b:0,z,1/2:4;c:1/4,z,x:4;d:y,z,x:8;"
	WyckoffSyms[195]  = "a:0,-z,0:4;b:1/2,-z,0:4;c:x,-z,1/4:4;d:x,-z,y:8;"
	WyckoffSyms[196]  = "a:0,0,z:4;b:1/4,y,z:4;c:x,y,z:8;"
	WyckoffSyms[197]  = "a:0,0,-z:4;b:y,1/4,-z:4;c:y,x,-z:8;"
	WyckoffSyms[198]  = "a:z,0,0:4;b:z,1/4,y:4;c:z,x,y:8;"
	WyckoffSyms[199]  = "a:-z,0,0:4;b:-z,y,1/4:4;c:-z,y,x:8;"
	WyckoffSyms[200]  = "a:0,z,0:4;b:y,z,1/4:4;c:y,z,x:8;"
	WyckoffSyms[201]  = "a:0,-z,0:4;b:1/4,-z,y:4;c:x,-z,y:8;"
	WyckoffSyms[202]  = "a:0,0,z:4;b:x,y,z:8;"
	WyckoffSyms[203]  = "a:0,0,-z:4;b:y,x,-z:8;"
	WyckoffSyms[204]  = "a:z,0,0:4;b:z,x,y:8;"
	WyckoffSyms[205]  = "a:-z,0,0:4;b:-z,y,x:8;"
	WyckoffSyms[206]  = "a:0,z,0:4;b:y,z,x:8;"
	WyckoffSyms[207]  = "a:0,-z,0:4;b:x,-z,y:8;"
	WyckoffSyms[208]  = "a:0,0,z:4;b:1/4,1/4,z:8;c:0,y,z:8;d:x,0,z:8;e:x,y,z:16;"
	WyckoffSyms[209]  = "a:z,0,0:4;b:z,1/4,1/4:8;c:z,0,y:8;d:z,x,0:8;e:z,x,y:16;"
	WyckoffSyms[210]  = "a:0,z,0:4;b:1/4,z,1/4:8;c:y,z,0:8;d:0,z,x:8;e:y,z,x:16;"
	WyckoffSyms[211]  = "a:0,0,z:8;b:x,y,z:16;"
	WyckoffSyms[212]  = "a:z,0,0:8;b:z,x,y:16;"
	WyckoffSyms[213]  = "a:0,z,0:8;b:y,z,x:16;"
	WyckoffSyms[214]  = "a:0,0,z:2;b:0,1/2,z:2;c:x,0,z:4;d:0,y,z:4;e:x,y,z:8;"
	WyckoffSyms[215]  = "a:z,0,0:2;b:z,0,1/2:2;c:z,x,0:4;d:z,0,y:4;e:z,x,y:8;"
	WyckoffSyms[216]  = "a:0,z,0:2;b:1/2,z,0:2;c:0,z,x:4;d:y,z,0:4;e:y,z,x:8;"
	WyckoffSyms[217]  = "a:0,0,z:4;b:0,1/2,z:4;c:x,y,z:8;"
	WyckoffSyms[218]  = "a:z,0,0:4;b:z,0,1/2:4;c:z,x,y:8;"
	WyckoffSyms[219]  = "a:0,z,0:4;b:1/2,z,0:4;c:y,z,x:8;"
	WyckoffSyms[220]  = "a:0,0,z:4;b:1/4,y,z:4;c:x,y,z:8;"
	WyckoffSyms[221]  = "a:0,0,-z:4;b:y,1/4,-z:4;c:y,x,-z:8;"
	WyckoffSyms[222]  = "a:z,0,0:4;b:z,1/4,y:4;c:z,x,y:8;"
	WyckoffSyms[223]  = "a:-z,0,0:4;b:-z,y,1/4:4;c:-z,y,x:8;"
	WyckoffSyms[224]  = "a:0,z,0:4;b:y,z,1/4:4;c:y,z,x:8;"
	WyckoffSyms[225]  = "a:0,-z,0:4;b:1/4,-z,y:4;c:x,-z,y:8;"
	WyckoffSyms[226]  = "a:0,0,0:1;b:1/2,0,0:1;c:0,0,1/2:1;d:1/2,0,1/2:1;e:0,1/2,0:1;f:1/2,1/2,0:1;g:0,1/2,1/2:1;h:1/2,1/2,1/2:1;i:x,0,0:2;j:x,0,1/2:2;"
	WyckoffSyms[226] += "k:x,1/2,0:2;l:x,1/2,1/2:2;m:0,y,0:2;n:0,y,1/2:2;o:1/2,y,0:2;p:1/2,y,1/2:2;q:0,0,z:2;r:0,1/2,z:2;s:1/2,0,z:2;t:1/2,1/2,z:2;"
	WyckoffSyms[226] += "u:0,y,z:4;v:1/2,y,z:4;w:x,0,z:4;x:x,1/2,z:4;y:x,y,0:4;z:x,y,1/2:4;A:x,y,z:8;"
	WyckoffSyms[227]  = "a:1/4,1/4,1/4:2;b:3/4,1/4,1/4:2;c:1/4,1/4,3/4:2;d:1/4,3/4,1/4:2;e:1/2,1/2,1/2:4;f:0,0,0:4;g:x,1/4,1/4:4;h:x,1/4,3/4:4;"
	WyckoffSyms[227] += "i:1/4,y,1/4:4;j:3/4,y,1/4:4;k:1/4,1/4,z:4;l:1/4,3/4,z:4;m:x,y,z:8;"
	WyckoffSyms[228]  = "a:1/4,1/4,1/4:2;b:3/4,1/4,1/4:2;c:1/4,1/4,3/4:2;d:1/4,3/4,1/4:2;e:1/2,1/2,1/2:4;f:0,0,0:4;g:x,1/4,1/4:4;h:x,1/4,3/4:4;"
	WyckoffSyms[228] += "i:1/4,y,1/4:4;j:3/4,y,1/4:4;k:1/4,1/4,z:4;l:1/4,3/4,z:4;m:x,y,z:8;"
	WyckoffSyms[229]  = "a:0,0,0:2;b:1/2,1/2,0:2;c:0,1/2,0:2;d:1/2,0,0:2;e:0,0,1/4:2;f:1/2,0,1/4:2;g:0,1/2,1/4:2;h:1/2,1/2,1/4:2;i:x,0,1/4:4;"
	WyckoffSyms[229] += "j:x,1/2,1/4:4;k:0,y,1/4:4;l:1/2,y,1/4:4;m:0,0,z:4;n:1/2,1/2,z:4;o:0,1/2,z:4;p:1/2,0,z:4;q:x,y,0:4;r:x,y,z:8;"
	WyckoffSyms[230]  = "a:0,0,0:2;b:0,1/2,1/2:2;c:0,0,1/2:2;d:0,1/2,0:2;e:1/4,0,0:2;f:1/4,1/2,0:2;g:1/4,0,1/2:2;h:1/4,1/2,1/2:2;i:1/4,x,0:4;"
	WyckoffSyms[230] += "j:1/4,x,1/2:4;k:1/4,0,y:4;l:1/4,1/2,y:4;m:z,0,0:4;n:z,1/2,1/2:4;o:z,0,1/2:4;p:z,1/2,0:4;q:0,x,y:4;r:z,x,y:8;"
	WyckoffSyms[231]  = "a:0,0,0:2;b:1/2,0,1/2:2;c:1/2,0,0:2;d:0,0,1/2:2;e:0,1/4,0:2;f:0,1/4,1/2:2;g:1/2,1/4,0:2;h:1/2,1/4,1/2:2;i:0,1/4,x:4;"
	WyckoffSyms[231] += "j:1/2,1/4,x:4;k:y,1/4,0:4;l:y,1/4,1/2:4;m:0,z,0:4;n:1/2,z,1/2:4;o:1/2,z,0:4;p:0,z,1/2:4;q:y,0,x:4;r:y,z,x:8;"
	WyckoffSyms[232]  = "a:1/4,1/4,0:2;b:3/4,1/4,0:2;c:3/4,1/4,1/2:2;d:1/4,1/4,1/2:2;e:0,0,0:4;f:0,0,1/2:4;g:x,1/4,0:4;h:x,1/4,1/2:4;i:1/4,y,0:4;"
	WyckoffSyms[232] += "j:1/4,y,1/2:4;k:1/4,1/4,z:4;l:1/4,3/4,z:4;m:x,y,z:8;"
	WyckoffSyms[233]  = "a:1/4,1/4,0:2;b:3/4,1/4,0:2;c:3/4,1/4,1/2:2;d:1/4,1/4,1/2:2;e:0,0,0:4;f:0,0,1/2:4;g:x,1/4,0:4;h:x,1/4,1/2:4;i:1/4,y,0:4;"
	WyckoffSyms[233] += "j:1/4,y,1/2:4;k:1/4,1/4,z:4;l:1/4,3/4,z:4;m:x,y,z:8;"
	WyckoffSyms[234]  = "a:0,1/4,1/4:2;b:0,3/4,1/4:2;c:1/2,3/4,1/4:2;d:1/2,1/4,1/4:2;e:0,0,0:4;f:1/2,0,0:4;g:0,x,1/4:4;h:1/2,x,1/4:4;i:0,1/4,y:4;"
	WyckoffSyms[234] += "j:1/2,1/4,y:4;k:z,1/4,1/4:4;l:z,1/4,3/4:4;m:z,x,y:8;"
	WyckoffSyms[235]  = "a:0,1/4,1/4:2;b:0,3/4,1/4:2;c:1/2,3/4,1/4:2;d:1/2,1/4,1/4:2;e:0,0,0:4;f:1/2,0,0:4;g:0,x,1/4:4;h:1/2,x,1/4:4;i:0,1/4,y:4;"
	WyckoffSyms[235] += "j:1/2,1/4,y:4;k:z,1/4,1/4:4;l:z,1/4,3/4:4;m:z,x,y:8;"
	WyckoffSyms[236]  = "a:1/4,0,1/4:2;b:1/4,0,3/4:2;c:1/4,1/2,3/4:2;d:1/4,1/2,1/4:2;e:0,0,0:4;f:0,1/2,0:4;g:1/4,0,x:4;h:1/4,1/2,x:4;i:y,0,1/4:4;"
	WyckoffSyms[236] += "j:y,1/2,1/4:4;k:1/4,z,1/4:4;l:3/4,z,1/4:4;m:y,z,x:8;"
	WyckoffSyms[237]  = "a:1/4,0,1/4:2;b:1/4,0,3/4:2;c:1/4,1/2,3/4:2;d:1/4,1/2,1/4:2;e:0,0,0:4;f:0,1/2,0:4;g:1/4,0,x:4;h:1/4,1/2,x:4;i:y,0,1/4:4;"
	WyckoffSyms[237] += "j:y,1/2,1/4:4;k:1/4,z,1/4:4;l:3/4,z,1/4:4;m:y,z,x:8;"
	WyckoffSyms[238]  = "a:0,0,0:2;b:0,1/2,0:2;c:0,0,1/2:2;d:0,1/2,1/2:2;e:1/4,0,z:2;f:1/4,1/2,z:2;g:0,y,0:4;h:0,y,1/2:4;i:x,0,z:4;j:x,1/2,z:4;"
	WyckoffSyms[238] += "k:1/4,y,z:4;l:x,y,z:8;"
	WyckoffSyms[239]  = "a:0,0,0:2;b:1/2,0,0:2;c:0,0,1/2:2;d:1/2,0,1/2:2;e:0,1/4,-z:2;f:1/2,1/4,-z:2;g:y,0,0:4;h:y,0,1/2:4;i:0,x,-z:4;j:1/2,x,-z:4;"
	WyckoffSyms[239] += "k:y,1/4,-z:4;l:y,x,-z:8;"
	WyckoffSyms[240]  = "a:0,0,0:2;b:0,0,1/2:2;c:1/2,0,0:2;d:1/2,0,1/2:2;e:z,1/4,0:2;f:z,1/4,1/2:2;g:0,0,y:4;h:1/2,0,y:4;i:z,x,0:4;j:z,x,1/2:4;"
	WyckoffSyms[240] += "k:z,1/4,y:4;l:z,x,y:8;"
	WyckoffSyms[241]  = "a:0,0,0:2;b:0,1/2,0:2;c:1/2,0,0:2;d:1/2,1/2,0:2;e:-z,0,1/4:2;f:-z,1/2,1/4:2;g:0,y,0:4;h:1/2,y,0:4;i:-z,0,x:4;j:-z,1/2,x:4;"
	WyckoffSyms[241] += "k:-z,y,1/4:4;l:-z,y,x:8;"
	WyckoffSyms[242]  = "a:0,0,0:2;b:1/2,0,0:2;c:0,1/2,0:2;d:1/2,1/2,0:2;e:0,z,1/4:2;f:1/2,z,1/4:2;g:y,0,0:4;h:y,1/2,0:4;i:0,z,x:4;j:1/2,z,x:4;"
	WyckoffSyms[242] += "k:y,z,1/4:4;l:y,z,x:8;"
	WyckoffSyms[243]  = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,0:2;d:0,1/2,1/2:2;e:1/4,-z,0:2;f:1/4,-z,1/2:2;g:0,0,y:4;h:0,1/2,y:4;i:x,-z,0:4;j:x,-z,1/2:4;"
	WyckoffSyms[243] += "k:1/4,-z,y:4;l:x,-z,y:8;"
	WyckoffSyms[244]  = "a:0,0,0:4;b:0,0,1/2:4;c:1/4,0,z:4;d:x,1/4,1/4:4;e:x,y,z:8;"
	WyckoffSyms[245]  = "a:0,0,0:4;b:0,0,1/2:4;c:0,1/4,-z:4;d:1/4,x,1/4:4;e:y,x,-z:8;"
	WyckoffSyms[246]  = "a:0,0,0:4;b:1/2,0,0:4;c:z,1/4,0:4;d:1/4,x,1/4:4;e:z,x,y:8;"
	WyckoffSyms[247]  = "a:0,0,0:4;b:1/2,0,0:4;c:-z,0,1/4:4;d:1/4,1/4,x:4;e:-z,y,x:8;"
	WyckoffSyms[248]  = "a:0,0,0:4;b:0,1/2,0:4;c:0,z,1/4:4;d:1/4,1/4,x:4;e:y,z,x:8;"
	WyckoffSyms[249]  = "a:0,0,0:4;b:0,1/2,0:4;c:1/4,-z,0:4;d:x,1/4,1/4:4;e:x,-z,y:8;"
	WyckoffSyms[250]  = "a:0,0,0:2;b:1/2,0,0:2;c:1/2,1/2,0:2;d:0,1/2,0:2;e:x,0,0:4;f:x,1/2,0:4;g:1/4,y,1/4:4;h:0,y,z:4;i:x,y,z:8;"
	WyckoffSyms[251]  = "a:0,0,0:2;b:0,1/2,0:2;c:1/2,1/2,0:2;d:1/2,0,0:2;e:0,x,0:4;f:1/2,x,0:4;g:y,1/4,1/4:4;h:y,0,-z:4;i:y,x,-z:8;"
	WyckoffSyms[252]  = "a:0,0,0:2;b:0,1/2,0:2;c:0,1/2,1/2:2;d:0,0,1/2:2;e:0,x,0:4;f:0,x,1/2:4;g:1/4,1/4,y:4;h:z,0,y:4;i:z,x,y:8;"
	WyckoffSyms[253]  = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,1/2:2;d:0,1/2,0:2;e:0,0,x:4;f:0,1/2,x:4;g:1/4,y,1/4:4;h:-z,y,0:4;i:-z,y,x:8;"
	WyckoffSyms[254]  = "a:0,0,0:2;b:0,0,1/2:2;c:1/2,0,1/2:2;d:1/2,0,0:2;e:0,0,x:4;f:1/2,0,x:4;g:y,1/4,1/4:4;h:y,z,0:4;i:y,z,x:8;"
	WyckoffSyms[255]  = "a:0,0,0:2;b:1/2,0,0:2;c:1/2,0,1/2:2;d:0,0,1/2:2;e:x,0,0:4;f:x,0,1/2:4;g:1/4,1/4,y:4;h:0,-z,y:4;i:x,-z,y:8;"
	WyckoffSyms[256]  = "a:0,0,0:4;b:0,1/2,0:4;c:0,y,1/4:4;d:1/4,0,z:4;e:1/4,1/2,z:4;f:x,y,z:8;"
	WyckoffSyms[257]  = "a:0,0,0:4;b:1/2,0,0:4;c:y,0,1/4:4;d:0,1/4,-z:4;e:1/2,1/4,-z:4;f:y,x,-z:8;"
	WyckoffSyms[258]  = "a:0,0,0:4;b:0,0,1/2:4;c:1/4,0,y:4;d:z,1/4,0:4;e:z,1/4,1/2:4;f:z,x,y:8;"
	WyckoffSyms[259]  = "a:0,0,0:4;b:0,1/2,0:4;c:1/4,y,0:4;d:-z,0,1/4:4;e:-z,1/2,1/4:4;f:-z,y,x:8;"
	WyckoffSyms[260]  = "a:0,0,0:4;b:1/2,0,0:4;c:y,1/4,0:4;d:0,z,1/4:4;e:1/2,z,1/4:4;f:y,z,x:8;"
	WyckoffSyms[261]  = "a:0,0,0:4;b:0,0,1/2:4;c:0,1/4,y:4;d:1/4,-z,0:4;e:1/4,-z,1/2:4;f:x,-z,y:8;"
	WyckoffSyms[262]  = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,0:2;d:0,1/2,1/2:2;e:0,0,z:4;f:0,1/2,z:4;g:x,y,0:4;h:x,y,1/2:4;i:x,y,z:8;"
	WyckoffSyms[263]  = "a:0,0,0:2;b:1/2,0,0:2;c:0,0,1/2:2;d:1/2,0,1/2:2;e:z,0,0:4;f:z,0,1/2:4;g:0,x,y:4;h:1/2,x,y:4;i:z,x,y:8;"
	WyckoffSyms[264]  = "a:0,0,0:2;b:0,1/2,0:2;c:1/2,0,0:2;d:1/2,1/2,0:2;e:0,z,0:4;f:1/2,z,0:4;g:y,0,x:4;h:y,1/2,x:4;i:y,z,x:8;"
	WyckoffSyms[265]  = "a:0,0,0:4;b:0,0,1/2:4;c:1/4,1/4,z:4;d:1/4,3/4,z:4;e:x,y,z:8;"
	WyckoffSyms[266]  = "a:0,0,0:4;b:1/2,0,0:4;c:z,1/4,1/4:4;d:z,1/4,3/4:4;e:z,x,y:8;"
	WyckoffSyms[267]  = "a:0,0,0:4;b:0,1/2,0:4;c:1/4,z,1/4:4;d:3/4,z,1/4:4;e:y,z,x:8;"
	WyckoffSyms[268]  = "a:0,0,0:4;b:1/2,0,0:4;c:x,1/4,0:4;d:x,y,1/4:4;e:x,y,z:8;"
	WyckoffSyms[269]  = "a:0,0,0:4;b:0,1/2,0:4;c:1/4,x,0:4;d:y,x,1/4:4;e:y,x,-z:8;"
	WyckoffSyms[270]  = "a:0,0,0:4;b:0,1/2,0:4;c:0,x,1/4:4;d:1/4,x,y:4;e:z,x,y:8;"
	WyckoffSyms[271]  = "a:0,0,0:4;b:0,0,1/2:4;c:0,1/4,x:4;d:1/4,y,x:4;e:-z,y,x:8;"
	WyckoffSyms[272]  = "a:0,0,0:4;b:0,0,1/2:4;c:1/4,0,x:4;d:y,1/4,x:4;e:y,z,x:8;"
	WyckoffSyms[273]  = "a:0,0,0:4;b:1/2,0,0:4;c:x,0,1/4:4;d:x,1/4,y:4;e:x,-z,y:8;"
	WyckoffSyms[274]  = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,0:2;d:0,1/2,1/2:2;e:0,0,z:4;f:0,1/2,z:4;g:x,y,0:4;h:x,y,z:8;"
	WyckoffSyms[275]  = "a:0,0,0:2;b:1/2,0,0:2;c:0,0,1/2:2;d:1/2,0,1/2:2;e:z,0,0:4;f:z,0,1/2:4;g:0,x,y:4;h:z,x,y:8;"
	WyckoffSyms[276]  = "a:0,0,0:2;b:0,1/2,0:2;c:1/2,0,0:2;d:1/2,1/2,0:2;e:0,z,0:4;f:1/2,z,0:4;g:y,0,x:4;h:y,z,x:8;"
	WyckoffSyms[277]  = "a:1/4,1/4,z:2;b:1/4,3/4,z:2;c:0,0,0:4;d:0,0,1/2:4;e:1/4,y,z:4;f:x,1/4,z:4;g:x,y,z:8;"
	WyckoffSyms[278]  = "a:1/4,1/4,z:2;b:1/4,3/4,z:2;c:0,0,0:4;d:0,0,1/2:4;e:1/4,y,z:4;f:x,1/4,z:4;g:x,y,z:8;"
	WyckoffSyms[279]  = "a:z,1/4,1/4:2;b:z,1/4,3/4:2;c:0,0,0:4;d:1/2,0,0:4;e:z,1/4,y:4;f:z,x,1/4:4;g:z,x,y:8;"
	WyckoffSyms[280]  = "a:z,1/4,1/4:2;b:z,1/4,3/4:2;c:0,0,0:4;d:1/2,0,0:4;e:z,1/4,y:4;f:z,x,1/4:4;g:z,x,y:8;"
	WyckoffSyms[281]  = "a:1/4,z,1/4:2;b:3/4,z,1/4:2;c:0,0,0:4;d:0,1/2,0:4;e:y,z,1/4:4;f:1/4,z,x:4;g:y,z,x:8;"
	WyckoffSyms[282]  = "a:1/4,z,1/4:2;b:3/4,z,1/4:2;c:0,0,0:4;d:0,1/2,0:4;e:y,z,1/4:4;f:1/4,z,x:4;g:y,z,x:8;"
	WyckoffSyms[283]  = "a:0,0,0:4;b:0,1/2,0:4;c:0,y,1/4:4;d:x,y,z:8;"
	WyckoffSyms[284]  = "a:0,0,0:4;b:1/2,0,0:4;c:y,0,1/4:4;d:y,x,-z:8;"
	WyckoffSyms[285]  = "a:0,0,0:4;b:0,0,1/2:4;c:1/4,0,y:4;d:z,x,y:8;"
	WyckoffSyms[286]  = "a:0,0,0:4;b:0,1/2,0:4;c:1/4,y,0:4;d:-z,y,x:8;"
	WyckoffSyms[287]  = "a:0,0,0:4;b:1/2,0,0:4;c:y,1/4,0:4;d:y,z,x:8;"
	WyckoffSyms[288]  = "a:0,0,0:4;b:0,0,1/2:4;c:0,1/4,y:4;d:x,-z,y:8;"
	WyckoffSyms[289]  = "a:0,0,0:4;b:0,0,1/2:4;c:x,y,z:8;"
	WyckoffSyms[290]  = "a:0,0,0:4;b:0,0,1/2:4;c:y,x,-z:8;"
	WyckoffSyms[291]  = "a:0,0,0:4;b:0,0,1/2:4;c:x,1/4,z:4;d:x,y,z:8;"
	WyckoffSyms[292]  = "a:0,0,0:4;b:0,0,1/2:4;c:1/4,x,-z:4;d:y,x,-z:8;"
	WyckoffSyms[293]  = "a:0,0,0:4;b:1/2,0,0:4;c:z,x,1/4:4;d:z,x,y:8;"
	WyckoffSyms[294]  = "a:0,0,0:4;b:1/2,0,0:4;c:-z,1/4,x:4;d:-z,y,x:8;"
	WyckoffSyms[295]  = "a:0,0,0:4;b:0,1/2,0:4;c:1/4,z,x:4;d:y,z,x:8;"
	WyckoffSyms[296]  = "a:0,0,0:4;b:0,1/2,0:4;c:x,-z,1/4:4;d:x,-z,y:8;"
	WyckoffSyms[297]  = "a:0,0,0:4;b:0,1/2,0:4;c:0,y,1/4:4;d:1/4,1/4,0:8;e:x,0,0:8;f:0,y,z:8;g:x,y,1/4:8;h:x,y,z:16;"
	WyckoffSyms[298]  = "a:0,0,0:4;b:1/2,0,0:4;c:y,0,1/4:4;d:1/4,1/4,0:8;e:0,x,0:8;f:y,0,-z:8;g:y,x,1/4:8;h:y,x,-z:16;"
	WyckoffSyms[299]  = "a:0,0,0:4;b:0,0,1/2:4;c:1/4,0,y:4;d:0,1/4,1/4:8;e:0,x,0:8;f:z,0,y:8;g:1/4,x,y:8;h:z,x,y:16;"
	WyckoffSyms[300]  = "a:0,0,0:4;b:0,1/2,0:4;c:1/4,y,0:4;d:0,1/4,1/4:8;e:0,0,x:8;f:-z,y,0:8;g:1/4,y,x:8;h:-z,y,x:16;"
	WyckoffSyms[301]  = "a:0,0,0:4;b:1/2,0,0:4;c:y,1/4,0:4;d:1/4,0,1/4:8;e:0,0,x:8;f:y,z,0:8;g:y,1/4,x:8;h:y,z,x:16;"
	WyckoffSyms[302]  = "a:0,0,0:4;b:0,0,1/2:4;c:0,1/4,y:4;d:1/4,0,1/4:8;e:x,0,0:8;f:0,-z,y:8;g:x,1/4,y:8;h:x,-z,y:16;"
	WyckoffSyms[303]  = "a:0,0,0:4;b:1/2,0,0:4;c:1/4,1/4,0:8;d:x,0,0:8;e:1/4,y,1/4:8;f:0,y,z:8;g:x,y,z:16;"
	WyckoffSyms[304]  = "a:0,0,0:4;b:0,1/2,0:4;c:1/4,1/4,0:8;d:0,x,0:8;e:y,1/4,1/4:8;f:y,0,-z:8;g:y,x,-z:16;"
	WyckoffSyms[305]  = "a:0,0,0:4;b:0,1/2,0:4;c:0,1/4,1/4:8;d:0,x,0:8;e:1/4,1/4,y:8;f:z,0,y:8;g:z,x,y:16;"
	WyckoffSyms[306]  = "a:0,0,0:4;b:0,0,1/2:4;c:0,1/4,1/4:8;d:0,0,x:8;e:1/4,y,1/4:8;f:-z,y,0:8;g:-z,y,x:16;"
	WyckoffSyms[307]  = "a:0,0,0:4;b:0,0,1/2:4;c:1/4,0,1/4:8;d:0,0,x:8;e:y,1/4,1/4:8;f:y,z,0:8;g:y,z,x:16;"
	WyckoffSyms[308]  = "a:0,0,0:4;b:1/2,0,0:4;c:1/4,0,1/4:8;d:x,0,0:8;e:1/4,1/4,y:8;f:0,-z,y:8;g:x,-z,y:16;"
	WyckoffSyms[309]  = "a:0,0,0:2;b:1/2,0,0:2;c:1/2,0,1/2:2;d:0,0,1/2:2;e:1/4,1/4,0:4;f:1/4,1/4,1/2:4;g:x,0,0:4;h:x,0,1/2:4;i:0,y,0:4;j:0,y,1/2:4;"
	WyckoffSyms[309] += "k:0,0,z:4;l:0,1/2,z:4;m:1/4,1/4,z:8;n:0,y,z:8;o:x,0,z:8;p:x,y,0:8;q:x,y,1/2:8;r:x,y,z:16;"
	WyckoffSyms[310]  = "a:0,0,0:2;b:0,1/2,0:2;c:1/2,1/2,0:2;d:1/2,0,0:2;e:0,1/4,1/4:4;f:1/2,1/4,1/4:4;g:0,x,0:4;h:1/2,x,0:4;i:0,0,y:4;j:1/2,0,y:4;"
	WyckoffSyms[310] += "k:z,0,0:4;l:z,0,1/2:4;m:z,1/4,1/4:8;n:z,0,y:8;o:z,x,0:8;p:0,x,y:8;q:1/2,x,y:8;r:z,x,y:16;"
	WyckoffSyms[311]  = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,1/2:2;d:0,1/2,0:2;e:1/4,0,1/4:4;f:1/4,1/2,1/4:4;g:0,0,x:4;h:0,1/2,x:4;i:y,0,0:4;j:y,1/2,0:4;"
	WyckoffSyms[311] += "k:0,z,0:4;l:1/2,z,0:4;m:1/4,z,1/4:8;n:y,z,0:8;o:0,z,x:8;p:y,0,x:8;q:y,1/2,x:8;r:y,z,x:16;"
	WyckoffSyms[312]  = "a:0,0,1/4:4;b:0,1/2,1/4:4;c:0,0,0:4;d:0,1/2,0:4;e:1/4,1/4,0:4;f:1/4,3/4,0:4;g:x,0,1/4:8;h:0,y,1/4:8;i:0,0,z:8;j:0,1/2,z:8;"
	WyckoffSyms[312] += "k:1/4,1/4,z:8;l:x,y,0:8;m:x,y,z:16;"
	WyckoffSyms[313]  = "a:1/4,0,0:4;b:1/4,0,1/2:4;c:0,0,0:4;d:0,0,1/2:4;e:0,1/4,1/4:4;f:0,1/4,3/4:4;g:1/4,x,0:8;h:1/4,0,y:8;i:z,0,0:8;j:z,0,1/2:8;"
	WyckoffSyms[313] += "k:z,1/4,1/4:8;l:0,x,y:8;m:z,x,y:16;"
	WyckoffSyms[314]  = "a:0,1/4,0:4;b:1/2,1/4,0:4;c:0,0,0:4;d:1/2,0,0:4;e:1/4,0,1/4:4;f:3/4,0,1/4:4;g:0,1/4,x:8;h:y,1/4,0:8;i:0,z,0:8;j:1/2,z,0:8;"
	WyckoffSyms[314] += "k:1/4,z,1/4:8;l:y,0,x:8;m:y,z,x:16;"
	WyckoffSyms[315]  = "a:1/4,0,0:4;b:1/4,0,1/2:4;c:0,0,0:4;d:0,0,1/2:4;e:1/4,1/4,0:4;f:1/4,1/4,1/2:4;g:0,1/4,z:4;h:x,0,0:8;i:x,0,1/2:8;j:1/4,y,0:8;"
	WyckoffSyms[315] += "k:1/4,y,1/2:8;l:1/4,0,z:8;m:0,y,z:8;n:x,1/4,z:8;o:x,y,z:16;"
	WyckoffSyms[316]  = "a:1/4,0,0:4;b:1/4,0,1/2:4;c:0,0,0:4;d:0,0,1/2:4;e:1/4,1/4,0:4;f:1/4,1/4,1/2:4;g:0,1/4,z:4;h:x,0,0:8;i:x,0,1/2:8;j:1/4,y,0:8;"
	WyckoffSyms[316] += "k:1/4,y,1/2:8;l:1/4,0,z:8;m:0,y,z:8;n:x,1/4,z:8;o:x,y,z:16;"
	WyckoffSyms[317]  = "a:0,1/4,0:4;b:1/2,1/4,0:4;c:0,0,0:4;d:1/2,0,0:4;e:0,1/4,1/4:4;f:1/2,1/4,1/4:4;g:z,0,1/4:4;h:0,x,0:8;i:1/2,x,0:8;j:0,1/4,y:8;"
	WyckoffSyms[317] += "k:1/2,1/4,y:8;l:z,1/4,0:8;m:z,0,y:8;n:z,x,1/4:8;o:z,x,y:16;"
	WyckoffSyms[318]  = "a:0,1/4,0:4;b:1/2,1/4,0:4;c:0,0,0:4;d:1/2,0,0:4;e:0,1/4,1/4:4;f:1/2,1/4,1/4:4;g:z,0,1/4:4;h:0,x,0:8;i:1/2,x,0:8;j:0,1/4,y:8;"
	WyckoffSyms[318] += "k:1/2,1/4,y:8;l:z,1/4,0:8;m:z,0,y:8;n:z,x,1/4:8;o:z,x,y:16;"
	WyckoffSyms[319]  = "a:0,0,1/4:4;b:0,1/2,1/4:4;c:0,0,0:4;d:0,1/2,0:4;e:1/4,0,1/4:4;f:1/4,1/2,1/4:4;g:1/4,z,0:4;h:0,0,x:8;i:0,1/2,x:8;j:y,0,1/4:8;"
	WyckoffSyms[319] += "k:y,1/2,1/4:8;l:0,z,1/4:8;m:y,z,0:8;n:1/4,z,x:8;o:y,z,x:16;"
	WyckoffSyms[320]  = "a:0,0,1/4:4;b:0,1/2,1/4:4;c:0,0,0:4;d:0,1/2,0:4;e:1/4,0,1/4:4;f:1/4,1/2,1/4:4;g:1/4,z,0:4;h:0,0,x:8;i:0,1/2,x:8;j:y,0,1/4:8;"
	WyckoffSyms[320] += "k:y,1/2,1/4:8;l:0,z,1/4:8;m:y,z,0:8;n:1/4,z,x:8;o:y,z,x:16;"
	WyckoffSyms[321]  = "a:0,1/4,1/4:4;b:0,1/4,3/4:4;c:1/4,3/4,0:8;d:0,0,0:8;e:x,1/4,1/4:8;f:0,y,1/4:8;g:0,1/4,z:8;h:1/4,0,z:8;i:x,y,z:16;"
	WyckoffSyms[322]  = "a:0,1/4,1/4:4;b:0,1/4,3/4:4;c:1/4,3/4,0:8;d:0,0,0:8;e:x,1/4,1/4:8;f:0,y,1/4:8;g:0,1/4,z:8;h:1/4,0,z:8;i:x,y,z:16;"
	WyckoffSyms[323]  = "a:0,1/4,1/4:4;b:0,1/4,3/4:4;c:1/4,3/4,0:8;d:0,0,0:8;e:x,1/4,1/4:8;f:0,y,1/4:8;g:0,1/4,z:8;h:1/4,0,z:8;i:x,y,z:16;"
	WyckoffSyms[324]  = "a:0,1/4,1/4:4;b:0,1/4,3/4:4;c:1/4,3/4,0:8;d:0,0,0:8;e:x,1/4,1/4:8;f:0,y,1/4:8;g:0,1/4,z:8;h:1/4,0,z:8;i:x,y,z:16;"
	WyckoffSyms[325]  = "a:1/4,0,1/4:4;b:3/4,0,1/4:4;c:0,1/4,3/4:8;d:0,0,0:8;e:1/4,x,1/4:8;f:1/4,0,y:8;g:z,0,1/4:8;h:z,1/4,0:8;i:z,x,y:16;"
	WyckoffSyms[326]  = "a:1/4,0,1/4:4;b:3/4,0,1/4:4;c:0,1/4,3/4:8;d:0,0,0:8;e:1/4,x,1/4:8;f:1/4,0,y:8;g:z,0,1/4:8;h:z,1/4,0:8;i:z,x,y:16;"
	WyckoffSyms[327]  = "a:1/4,0,1/4:4;b:3/4,0,1/4:4;c:0,1/4,3/4:8;d:0,0,0:8;e:1/4,x,1/4:8;f:1/4,0,y:8;g:z,0,1/4:8;h:z,1/4,0:8;i:z,x,y:16;"
	WyckoffSyms[328]  = "a:1/4,0,1/4:4;b:3/4,0,1/4:4;c:0,1/4,3/4:8;d:0,0,0:8;e:1/4,x,1/4:8;f:1/4,0,y:8;g:z,0,1/4:8;h:z,1/4,0:8;i:z,x,y:16;"
	WyckoffSyms[329]  = "a:1/4,1/4,0:4;b:1/4,3/4,0:4;c:3/4,0,1/4:8;d:0,0,0:8;e:1/4,1/4,x:8;f:y,1/4,0:8;g:1/4,z,0:8;h:0,z,1/4:8;i:y,z,x:16;"
	WyckoffSyms[330]  = "a:1/4,1/4,0:4;b:1/4,3/4,0:4;c:3/4,0,1/4:8;d:0,0,0:8;e:1/4,1/4,x:8;f:y,1/4,0:8;g:1/4,z,0:8;h:0,z,1/4:8;i:y,z,x:16;"
	WyckoffSyms[331]  = "a:1/4,1/4,0:4;b:1/4,3/4,0:4;c:3/4,0,1/4:8;d:0,0,0:8;e:1/4,1/4,x:8;f:y,1/4,0:8;g:1/4,z,0:8;h:0,z,1/4:8;i:y,z,x:16;"
	WyckoffSyms[332]  = "a:1/4,1/4,0:4;b:1/4,3/4,0:4;c:3/4,0,1/4:8;d:0,0,0:8;e:1/4,1/4,x:8;f:y,1/4,0:8;g:1/4,z,0:8;h:0,z,1/4:8;i:y,z,x:16;"
	WyckoffSyms[333]  = "a:0,0,0:4;b:0,0,1/2:4;c:0,1/4,1/4:8;d:1/4,0,1/4:8;e:1/4,1/4,0:8;f:1/4,1/4,1/4:8;g:x,0,0:8;h:0,y,0:8;i:0,0,z:8;j:1/4,1/4,z:16;"
	WyckoffSyms[333] += "k:1/4,y,1/4:16;l:x,1/4,1/4:16;m:0,y,z:16;n:x,0,z:16;o:x,y,0:16;p:x,y,z:32;"
	WyckoffSyms[334]  = "a:1/8,1/8,1/8:8;b:1/8,1/8,2/3:8;c:0,0,0:16;d:1/2,1/2,1/2:16;e:x,1/8,1/8:16;f:1/8,y,1/8:16;g:1/8,1/8,z:16;h:x,y,z:32;"
	WyckoffSyms[335]  = "a:1/8,1/8,1/8:8;b:1/8,1/8,2/3:8;c:0,0,0:16;d:1/2,1/2,1/2:16;e:x,1/8,1/8:16;f:1/8,y,1/8:16;g:1/8,1/8,z:16;h:x,y,z:32;"
	WyckoffSyms[336]  = "a:0,0,0:2;b:0,1/2,1/2:2;c:1/2,1/2,0:2;d:1/2,0,1/2:2;e:x,0,0:4;f:x,1/2,0:4;g:0,y,0:4;h:0,y,1/2:4;i:0,0,z:4;j:1/2,0,z:4;"
	WyckoffSyms[336] += "k:1/4,1/4,1/4:8;l:0,y,z:8;m:x,0,z:8;n:x,y,0:8;o:x,y,z:16;"
	WyckoffSyms[337]  = "a:0,0,1/4:4;b:1/2,0,1/4:4;c:0,0,0:4;d:1/2,0,0:4;e:1/4,1/4,1/4:8;f:x,0,1/4:8;g:0,y,1/4:8;h:0,0,z:8;i:0,1/2,z:8;j:x,y,0:8;"
	WyckoffSyms[337] += "k:x,y,z:16;"
	WyckoffSyms[338]  = "a:1/4,0,0:4;b:1/4,1/2,0:4;c:0,0,0:4;d:0,1/2,0:4;e:1/4,1/4,1/4:8;f:1/4,x,0:8;g:1/4,0,y:8;h:z,0,0:8;i:z,0,1/2:8;j:0,x,y:8;"
	WyckoffSyms[338] += "k:z,x,y:16;"
	WyckoffSyms[339]  = "a:0,1/4,0:4;b:0,1/4,1/2:4;c:0,0,0:4;d:0,0,1/2:4;e:1/4,1/4,1/4:8;f:0,1/4,x:8;g:y,1/4,0:8;h:0,z,0:8;i:1/2,z,0:8;j:y,0,x:8;"
	WyckoffSyms[339] += "k:y,z,x:16;"
	WyckoffSyms[340]  = "a:0,0,0:8;b:1/4,1/4,1/4:8;c:x,0,1/4:8;d:1/4,y,0:8;e:0,1/4,z:8;f:x,y,z:16;"
	WyckoffSyms[341]  = "a:0,0,0:8;b:1/4,1/4,1/4:8;c:x,0,1/4:8;d:1/4,y,0:8;e:0,1/4,z:8;f:x,y,z:16;"
	WyckoffSyms[342]  = "a:0,0,0:4;b:0,0,1/2:4;c:1/4,1/4,1/4:4;d:1/4,1/4,3/4:4;e:0,1/4,z:4;f:x,0,0:8;g:1/4,y,1/4:8;h:0,y,z:8;i:x,1/4,z:8;j:x,y,z:16;"
	WyckoffSyms[343]  = "a:0,0,0:4;b:0,0,1/2:4;c:1/4,1/4,1/4:4;d:1/4,1/4,3/4:4;e:0,1/4,z:4;f:x,0,0:8;g:1/4,y,1/4:8;h:0,y,z:8;i:x,1/4,z:8;j:x,y,z:16;"
	WyckoffSyms[344]  = "a:0,0,0:4;b:1/2,0,0:4;c:1/4,1/4,1/4:4;d:3/4,1/4,1/4:4;e:z,0,1/4:4;f:0,x,0:8;g:1/4,1/4,y:8;h:z,0,y:8;i:z,x,1/4:8;j:z,x,y:16;"
	WyckoffSyms[345]  = "a:0,0,0:4;b:1/2,0,0:4;c:1/4,1/4,1/4:4;d:3/4,1/4,1/4:4;e:z,0,1/4:4;f:0,x,0:8;g:1/4,1/4,y:8;h:z,0,y:8;i:z,x,1/4:8;j:z,x,y:16;"
	WyckoffSyms[346]  = "a:0,0,0:4;b:0,1/2,0:4;c:1/4,1/4,1/4:4;d:1/4,3/4,1/4:4;e:1/4,z,0:4;f:0,0,x:8;g:y,1/4,1/4:8;h:y,z,0:8;i:1/4,z,x:8;j:y,z,x:16;"
	WyckoffSyms[347]  = "a:0,0,0:4;b:0,1/2,0:4;c:1/4,1/4,1/4:4;d:1/4,3/4,1/4:4;e:1/4,z,0:4;f:0,0,x:8;g:y,1/4,1/4:8;h:y,z,0:8;i:1/4,z,x:8;j:y,z,x:16;"
	// Tetragonal SG[75,142]  SG_idNum [349,429]   (81 idNums)
	WyckoffSyms[348]  = "a:0,0,z:1;b:1/2,1/2,z:1;c:0,1/2,z:2;d:x,y,z:4;"
	WyckoffSyms[349]  = "a:x,y,z:4;"
	WyckoffSyms[350]  = "a:0,0,z:2;b:1/2,1/2,z:2;c:0,1/2,z:2;d:x,y,z:4;"
	WyckoffSyms[351]  = "a:x,y,z:4;"
	WyckoffSyms[352]  = "a:0,0,z:2;b:0,1/2,z:4;c:x,y,z:8;"
	WyckoffSyms[353]  = "a:0,0,z:4;b:x,y,z:8;"
	WyckoffSyms[354]  = "a:0,0,0:1;b:0,0,1/2:1;c:1/2,1/2,0:1;d:1/2,1/2,1/2:1;e:0,0,z:2;f:1/2,1/2,z:2;g:0,1/2,z:2;h:x,y,z:4;"
	WyckoffSyms[355]  = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,1/4:2;d:0,1/2,3/4:2;e:0,0,z:4;f:0,1/2,z:4;g:x,y,z:8;"
	WyckoffSyms[356]  = "a:0,0,0:1;b:0,0,1/2:1;c:1/2,1/2,0:1;d:1/2,1/2,1/2:1;e:0,1/2,0:2;f:0,1/2,1/2:2;g:0,0,z:2;h:1/2,1/2,z:2;i:0,1/2,z:4;j:x,y,0:4;"
	WyckoffSyms[356] += "k:x,y,1/2:4;l:x,y,z:8;"
	WyckoffSyms[357]  = "a:0,0,0:2;b:1/2,1/2,0:2;c:0,1/2,0:2;d:0,1/2,1/2:2;e:0,0,1/4:2;f:1/2,1/2,1/4:2;g:0,0,z:4;h:1/2,1/2,z:4;i:0,1/2,z:4;j:x,y,0:4;"
	WyckoffSyms[357] += "k:x,y,z:8;"
	WyckoffSyms[358]  = "a:1/4,3/4,0:2;b:1/4,3/4,1/2:2;c:1/4,1/4,z:2;d:0,0,0:4;e:0,0,1/2:4;f:1/4,3/4,z:4;g:x,y,z:8;"
	WyckoffSyms[359]  = "a:1/4,3/4,0:2;b:1/4,3/4,1/2:2;c:1/4,1/4,z:2;d:0,0,0:4;e:0,0,1/2:4;f:1/4,3/4,z:4;g:x,y,z:8;"
	WyckoffSyms[360]  = "a:1/4,1/4,1/4:2;b:1/4,1/4,3/4:2;c:0,0,0:4;d:0,0,1/2:4;e:3/4,1/4,z:4;f:1/4,1/4,z:4;g:x,y,z:8;"
	WyckoffSyms[361]  = "a:1/4,1/4,1/4:2;b:1/4,1/4,3/4:2;c:0,0,0:4;d:0,0,1/2:4;e:3/4,1/4,z:4;f:1/4,1/4,z:4;g:x,y,z:8;"
	WyckoffSyms[362]  = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,0:4;d:0,1/2,1/4:4;e:0,0,z:4;f:1/4,1/4,1/4:8;g:0,1/2,z:8;h:x,y,0:8;i:x,y,z:16;"
	WyckoffSyms[363]  = "a:0,1/4,1/8:4;b:0,1/4,2/3:4;c:0,0,0:8;d:0,0,1/2:8;e:0,1/4,z:8;f:x,y,z:16;"
	WyckoffSyms[364]  = "a:0,1/4,1/8:4;b:0,1/4,2/3:4;c:0,0,0:8;d:0,0,1/2:8;e:0,1/4,z:8;f:x,y,z:16;"
	WyckoffSyms[365]  = "a:0,0,0:1;b:0,0,1/2:1;c:1/2,1/2,0:1;d:1/2,1/2,1/2:1;e:1/2,0,0:2;f:1/2,0,1/2:2;g:0,0,z:2;h:1/2,1/2,z:2;i:0,1/2,z:4;j:x,x,0:4;"
	WyckoffSyms[365] += "k:x,x,1/2:4;l:x,0,0:4;m:x,1/2,1/2:4;n:x,0,1/2:4;o:x,1/2,0:4;p:x,y,z:8;"
	WyckoffSyms[366]  = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,z:2;d:0,0,z:4;e:x,x,0:4;f:x,x,1/2:4;g:x,y,z:8;"
	WyckoffSyms[367]  = "a:0,y,0:4;b:1/2,y,0:4;c:x,x,1/3:4;d:x,y,z:8;"
	WyckoffSyms[368]  = "a:x,x,0:4;b:x,y,z:8;"
	WyckoffSyms[369]  = "a:0,0,0:2;b:1/2,1/2,0:2;c:0,1/2,0:2;d:0,1/2,1/2:2;e:0,0,1/4:2;f:1/2,1/2,1/4:2;g:0,0,z:4;h:1/2,1/2,z:4;i:0,1/2,z:4;j:x,0,0:4;"
	WyckoffSyms[369] += "k:x,1/2,1/2:4;l:x,0,1/2:4;m:x,1/2,0:4;n:x,x,1/4:4;o:x,x,3/4:4;p:x,y,z:8;"
	WyckoffSyms[370]  = "a:0,0,0:2;b:0,0,1/2:2;c:0,0,z:4;d:0,1/2,z:4;e:x,x,0:4;f:x,x,1/2:4;g:x,y,z:8;"
	WyckoffSyms[371]  = "a:0,y,0:4;b:1/2,y,0:4;c:x,x,2/3:4;d:x,y,z:8;"
	WyckoffSyms[372]  = "a:x,x,0:4;b:x,y,z:8;"
	WyckoffSyms[373]  = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,0:4;d:0,1/2,1/4:4;e:0,0,z:4;f:0,1/2,z:8;g:x,x,0:8;h:x,0,0:8;i:x,0,1/2:8;j:x,x+1/2,1/4:8;"
	WyckoffSyms[373] += "k:x,y,z:16;"
	WyckoffSyms[374]  = "a:0,0,0:4;b:0,0,1/2:4;c:0,0,z:8;d:x,x,0:8;e:-x,x,0:8;f:x,1/4,1/8:8;g:x,y,z:16;"
	WyckoffSyms[375]  = "a:0,0,z:1;b:1/2,1/2,z:1;c:1/2,0,z:2;d:x,x,z:4;e:x,0,z:4;f:x,1/2,z:4;g:x,y,z:8;"
	WyckoffSyms[376]  = "a:0,0,z:2;b:1/2,0,z:2;c:x,x+1/2,z:4;d:x,y,z:8;"
	WyckoffSyms[377]  = "a:0,0,z:2;b:1/2,1/2,z:2;c:0,1/2,z:4;d:x,x,z:4;e:x,y,z:8;"
	WyckoffSyms[378]  = "a:0,0,z:2;b:0,1/2,z:4;c:x,x,z:4;d:x,y,z:8;"
	WyckoffSyms[379]  = "a:0,0,z:2;b:1/2,1/2,z:2;c:0,1/2,z:4;d:x,y,z:8;"
	WyckoffSyms[380]  = "a:0,0,z:2;b:0,1/2,z:4;c:x,y,z:8;"
	WyckoffSyms[381]  = "a:0,0,z:2;b:1/2,1/2,z:2;c:0,1/2,z:2;d:x,0,z:4;e:x,1/2,z:4;f:x,y,z:8;"
	WyckoffSyms[382]  = "a:0,0,z:4;b:0,1/2,z:4;c:x,y,z:8;"
	WyckoffSyms[383]  = "a:0,0,z:2;b:0,1/2,z:4;c:x,x,z:8;d:x,0,z:8;e:x,y,z:16;"
	WyckoffSyms[384]  = "a:0,0,z:4;b:1/2,0,z:4;c:x,x+1/2,z:8;d:x,y,z:16;"
	WyckoffSyms[385]  = "a:0,0,z:4;b:0,y,z:8;c:x,y,z:16;"
	WyckoffSyms[386]  = "a:0,0,z:8;b:x,y,z:16;"
	WyckoffSyms[387]  = "a:0,0,0:1;b:1/2,1/2,1/2:1;c:0,0,1/2:1;d:1/2,1/2,0:1;e:1/2,0,0:2;f:1/2,0,1/2:2;g:0,0,z:2;h:1/2,1/2,z:2;i:x,0,0:4;j:x,1/2,1/2:4;"
	WyckoffSyms[387] += "k:x,0,1/2:4;l:x,1/2,0:4;m:0,1/2,z:4;n:x,x,z:4;o:x,y,z:8;"
	WyckoffSyms[388]  = "a:0,0,1/4:2;b:1/2,0,1/4:2;c:1/2,1/2,1/4:2;d:0,1/2,1/4:2;e:0,0,0:2;f:1/2,1/2,0:2;g:x,0,1/4:4;h:1/2,y,1/4:4;i:x,1/2,1/4:4;"
	WyckoffSyms[388] += "j:0,y,1/4:4;k:0,0,z:4;l:1/2,1/2,z:4;m:0,1/2,z:4;n:x,y,z:8;"
	WyckoffSyms[389]  = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,z:2;d:0,0,z:4;e:x,x+1/2,z:4;f:x,y,z:8;"
	WyckoffSyms[390]  = "a:0,0,0:2;b:0,0,1/2:2;c:0,0,z:4;d:0,1/2,z:4;e:x,y,z:8;"
	WyckoffSyms[391]  = "a:0,0,0:1;b:1/2,1/2,0:1;c:1/2,1/2,1/2:1;d:0,0,1/2:1;e:0,0,z:2;f:1/2,1/2,z:2;g:0,1/2,z:2;h:x,x,0:4;i:x,x,1/2:4;j:x,0,z:4;"
	WyckoffSyms[391] += "k:x,1/2,z:4;l:x,y,z:8;"
	WyckoffSyms[392]  = "a:0,0,1/4:2;b:1/2,1/2,1/4:2;c:0,0,0:2;d:1/2,1/2,0:2;e:x,x,1/4:4;f:x,x,3/4:4;g:0,0,z:4;h:1/2,1/2,z:4;i:0,1/2,z:4;j:x,y,z:8;"
	WyckoffSyms[393]  = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,0:2;d:0,1/2,1/2:2;e:0,0,z:4;f:0,1/2,z:4;g:x,x+1/2,0:4;h:x,x+1/2,1/2:4;i:x,y,z:8;"
	WyckoffSyms[394]  = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,1/4:2;d:0,1/2,3/4:2;e:0,0,z:4;f:x,-x+1/2,1/4:4;g:x,x+1/2,1/4:4;h:0,1/2,z:4;i:x,y,z:8;"
	WyckoffSyms[395]  = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,1/4:2;d:0,1/2,3/4:2;e:0,0,z:4;f:0,1/2,z:4;g:x,x,0:8;h:x,x+1/2,1/4:8;i:x,0,z:8;j:x,y,z:16;"
	WyckoffSyms[396]  = "a:0,0,1/4:4;b:0,0,0:4;c:0,1/2,1/4:4;d:0,1/2,0:4;e:x,x,1/4:8;f:0,0,z:8;g:0,1/2,z:8;h:x,x+1/2,0:8;i:x,y,z:16;"
	WyckoffSyms[397]  = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,0:4;d:0,1/2,1/4:4;e:0,0,z:4;f:x,0,0:8;g:x,0,1/2:8;h:0,1/2,z:8;i:x,x,z:8;j:x,y,z:16;"
	WyckoffSyms[398]  = "a:0,0,0:4;b:0,0,1/2:4;c:0,0,z:8;d:x,1/4,1/8:8;e:x,y,z:16;"
	WyckoffSyms[399]  = "a:0,0,0:1;b:0,0,1/2:1;c:1/2,1/2,0:1;d:1/2,1/2,1/2:1;e:0,1/2,1/2:2;f:0,1/2,0:2;g:0,0,z:2;h:1/2,1/2,z:2;i:0,1/2,z:4;j:x,x,0:4;"
	WyckoffSyms[399] += "k:x,x,1/2:4;l:x,0,0:4;m:x,0,1/2:4;n:x,1/2,0:4;o:x,1/2,1/2:4;p:x,y,0:8;q:x,y,1/2:8;r:x,x,z:8;s:x,0,z:8;t:x,1/2,z:8;u:x,y,z:16;"
	WyckoffSyms[400]  = "a:0,0,1/4:2;b:0,0,0:2;c:1/2,1/2,1/4:2;d:1/2,1/2,0:2;e:0,1/2,0:4;f:0,1/2,1/4:4;g:0,0,z:4;h:1/2,1/2,z:4;i:0,1/2,z:8;j:x,x,1/4:8;"
	WyckoffSyms[400] += "k:x,0,1/4:8;l:x,1/2,1/4:8;m:x,y,0:8;n:x,y,z:16;"
	WyckoffSyms[401]  = "a:1/4,1/4,0:2;b:1/4,1/4,1/2:2;c:3/4,1/4,0:2;d:3/4,1/4,1/2:2;e:0,0,0:4;f:0,0,1/2:4;g:1/4,1/4,z:4;h:3/4,1/4,z:4;i:x,x,0:8;"
	WyckoffSyms[401] += "j:x,x,1/2:8;k:x,1/4,0:8;l:x,1/4,1/2:8;m:x,-x,z:8;n:x,y,z:16;"
	WyckoffSyms[402]  = "a:1/4,1/4,0:2;b:1/4,1/4,1/2:2;c:3/4,1/4,0:2;d:3/4,1/4,1/2:2;e:0,0,0:4;f:0,0,1/2:4;g:1/4,1/4,z:4;h:3/4,1/4,z:4;i:x,x,0:8;"
	WyckoffSyms[402] += "j:x,x,1/2:8;k:x,1/4,0:8;l:x,1/4,1/2:8;m:x,-x,z:8;n:x,y,z:16;"
	WyckoffSyms[403]  = "a:1/4,1/4,1/4:2;b:1/4,1/4,3/4:2;c:1/4,3/4,3/4:4;d:1/4,3/4,0:4;e:1/4,1/4,z:4;f:0,0,0:8;g:1/4,3/4,z:8;h:x,x,1/4:8;i:x,1/4,1/4:8;"
	WyckoffSyms[403] += "j:x,3/4,1/4:8;k:x,y,z:16;"
	WyckoffSyms[404]  = "a:1/4,1/4,1/4:2;b:1/4,1/4,3/4:2;c:1/4,3/4,3/4:4;d:1/4,3/4,0:4;e:1/4,1/4,z:4;f:0,0,0:8;g:1/4,3/4,z:8;h:x,x,1/4:8;i:x,1/4,1/4:8;"
	WyckoffSyms[404] += "j:x,3/4,1/4:8;k:x,y,z:16;"
	WyckoffSyms[405]  = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,1/2:2;d:0,1/2,0:2;e:0,0,z:4;f:0,1/2,z:4;g:x,x+1/2,0:4;h:x,x+1/2,1/2:4;i:x,y,0:8;j:x,y,1/2:8;"
	WyckoffSyms[405] += "k:x,x+1/2,z:8;l:x,y,z:16;"
	WyckoffSyms[406]  = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,0:4;d:0,1/2,1/4:4;e:0,0,z:4;f:0,1/2,z:8;g:x,x+1/2,1/4:8;h:x,y,0:8;i:x,y,z:16;"
	WyckoffSyms[407]  = "a:3/4,1/4,0:2;b:3/4,1/4,1/2:2;c:1/4,1/4,z:2;d:0,0,0:4;e:0,0,1/2:4;f:3/4,1/4,z:4;g:x,-x,0:8;h:x,-x,1/2:8;i:1/4,y,z:8;j:x,x,z:8;"
	WyckoffSyms[407] += "k:x,y,z:16;"
	WyckoffSyms[408]  = "a:3/4,1/4,0:2;b:3/4,1/4,1/2:2;c:1/4,1/4,z:2;d:0,0,0:4;e:0,0,1/2:4;f:3/4,1/4,z:4;g:x,-x,0:8;h:x,-x,1/2:8;i:1/4,y,z:8;j:x,x,z:8;"
	WyckoffSyms[408] += "k:x,y,z:16;"
	WyckoffSyms[409]  = "a:3/4,1/4,1/4:4;b:3/4,1/4,0:4;c:1/4,1/4,z:4;d:0,0,0:8;e:3/4,1/4,z:8;f:x,-x,1/4:8;g:x,y,z:16;"
	WyckoffSyms[410]  = "a:3/4,1/4,1/4:4;b:3/4,1/4,0:4;c:1/4,1/4,z:4;d:0,0,0:8;e:3/4,1/4,z:8;f:x,-x,1/4:8;g:x,y,z:16;"
	WyckoffSyms[411]  = "a:0,0,0:2;b:1/2,1/2,0:2;c:0,1/2,0:2;d:0,1/2,1/2:2;e:0,0,1/4:2;f:1/2,1/2,1/4:2;g:0,0,z:4;h:1/2,1/2,z:4;i:0,1/2,z:4;j:x,0,0:4;"
	WyckoffSyms[411] += "k:x,1/2,1/2:4;l:x,0,1/2:4;m:x,1/2,0:4;n:x,x,1/4:8;o:0,y,z:8;p:1/2,y,z:8;q:x,y,0:8;r:x,y,z:16;"
	WyckoffSyms[412]  = "a:0,0,0:2;b:0,0,1/4:2;c:1/2,1/2,0:2;d:1/2,1/2,1/4:2;e:0,1/2,1/4:4;f:0,1/2,0:4;g:0,0,z:4;h:1/2,1/2,z:4;i:x,x,0:4;j:x,x,1/2:4;"
	WyckoffSyms[412] += "k:0,1/2,z:8;l:x,0,1/4:8;m:x,1/2,1/4:8;n:x,y,0:8;o:x,x,z:8;p:x,y,z:16;"
	WyckoffSyms[413]  = "a:1/4,1/4,0:4;b:3/4,1/4,0:4;c:1/4,1/4,1/4:4;d:3/4,1/4,3/4:4;e:0,0,0:8;f:1/4,1/4,z:8;g:3/4,1/4,z:8;h:x,1/4,0:8;i:x,1/4,1/2:8;"
	WyckoffSyms[413] += "j:x,x,1/4:8;k:x,y,z:16;"
	WyckoffSyms[414]  = "a:1/4,1/4,0:4;b:3/4,1/4,0:4;c:1/4,1/4,1/4:4;d:3/4,1/4,3/4:4;e:0,0,0:8;f:1/4,1/4,z:8;g:3/4,1/4,z:8;h:x,1/4,0:8;i:x,1/4,1/2:8;"
	WyckoffSyms[414] += "j:x,x,1/4:8;k:x,y,z:16;"
	WyckoffSyms[415]  = "a:1/4,3/4,1/4:2;b:3/4,1/4,1/4:2;c:1/4,1/4,1/4:4;d:1/4,1/4,0:4;e:0,0,1/2:4;f:0,0,0:4;g:3/4,1/4,z:4;h:1/4,1/4,z:8;i:x,1/4,3/4:8;"
	WyckoffSyms[415] += "j:x,1/4,1/4:8;k:x,x,0:8;l:x,x,1/2:8;m:x,-x,z:8;n:x,y,z:16;"
	WyckoffSyms[416]  = "a:1/4,3/4,1/4:2;b:3/4,1/4,1/4:2;c:1/4,1/4,1/4:4;d:1/4,1/4,0:4;e:0,0,1/2:4;f:0,0,0:4;g:3/4,1/4,z:4;h:1/4,1/4,z:8;i:x,1/4,3/4:8;"
	WyckoffSyms[416] += "j:x,1/4,1/4:8;k:x,x,0:8;l:x,x,1/2:8;m:x,-x,z:8;n:x,y,z:16;"
	WyckoffSyms[417]  = "a:0,0,0:4;b:0,0,1/4:4;c:0,1/2,0:4;d:0,1/2,1/4:4;e:0,0,z:8;f:0,1/2,z:8;g:x,x+1/2,1/4:8;h:x,y,0:8;i:x,y,z:16;"
	WyckoffSyms[418]  = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,0:4;d:0,1/2,1/4:4;e:0,0,z:4;f:x,x,0:4;g:x,-x,0:4;h:0,1/2,z:8;i:x,y,0:8;j:x,x,z:8;k:x,y,z:16;"
	WyckoffSyms[419]  = "a:3/4,1/4,3/4:2;b:3/4,1/4,1/4:2;c:3/4,1/4,z:4;d:1/4,1/4,z:4;e:0,0,0:8;f:x,-x,1/4:8;g:1/4,y,z:8;h:x,y,z:16;"
	WyckoffSyms[420]  = "a:3/4,1/4,3/4:2;b:3/4,1/4,1/4:2;c:3/4,1/4,z:4;d:1/4,1/4,z:4;e:0,0,0:8;f:x,-x,1/4:8;g:1/4,y,z:8;h:x,y,z:16;"
	WyckoffSyms[421]  = "a:3/4,1/4,0:4;b:3/4,1/4,3/4:4;c:0,0,1/2:4;d:0,0,0:4;e:1/4,1/4,z:4;f:3/4,1/4,z:8;g:x,-x,1/2:8;h:x,-x,0:8;i:x,x,z:8;j:x,y,z:16;"
	WyckoffSyms[422]  = "a:3/4,1/4,0:4;b:3/4,1/4,3/4:4;c:0,0,1/2:4;d:0,0,0:4;e:1/4,1/4,z:4;f:3/4,1/4,z:8;g:x,-x,1/2:8;h:x,-x,0:8;i:x,x,z:8;j:x,y,z:16;"
	WyckoffSyms[423]  = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,0:4;d:0,1/2,1/4:4;e:0,0,z:4;f:1/4,1/4,1/4:8;g:0,1/2,z:8;h:x,x,0:8;i:x,0,0:8;j:x,1/2,0:8;"
	WyckoffSyms[423] += "k:x,x+1/2,1/4:16;l:x,y,0:16;m:x,x,z:16;n:0,y,z:16;o:x,y,z:32;"
	WyckoffSyms[424]  = "a:0,0,1/4:4;b:0,1/2,1/4:4;c:0,0,0:4;d:0,1/2,0:4;e:1/4,1/4,1/4:8;f:0,0,z:8;g:0,1/2,z:8;h:x,x+1/2,0:8;i:x,x,1/4:16;j:x,0,1/4:16;"
	WyckoffSyms[424] += "k:x,y,0:16;l:x,x+1/2,z:16;m:x,y,z:32;"
	WyckoffSyms[425]  = "a:0,3/4,1/8:4;b:0,1/4,1/3:4;c:0,0,0:8;d:0,0,1/2:8;e:0,1/4,z:8;f:x,0,0:16;g:x,x+1/4,7/8:16;h:0,y,z:16;i:x,y,z:32;"
	WyckoffSyms[426]  = "a:0,3/4,1/8:4;b:0,1/4,1/3:4;c:0,0,0:8;d:0,0,1/2:8;e:0,1/4,z:8;f:x,0,0:16;g:x,x+1/4,7/8:16;h:0,y,z:16;i:x,y,z:32;"
	WyckoffSyms[427]  = "a:0,1/4,1/3:8;b:0,1/4,1/8:8;c:0,0,0:16;d:0,1/4,z:16;e:x,0,1/4:16;f:x,x+1/4,1/8:16;g:x,y,z:32;"
	WyckoffSyms[428]  = "a:0,1/4,1/3:8;b:0,1/4,1/8:8;c:0,0,0:16;d:0,1/4,z:16;e:x,0,1/4:16;f:x,x+1/4,1/8:16;g:x,y,z:32;"
	// Trigonal SG[143,167]  SG_idNum [430,461]   (32 idNums)
	WyckoffSyms[429]  = "a:0,0,z:1;b:1/3,2/3,z:1;c:2/3,1/3,z:1;d:x,y,z:3;"
	WyckoffSyms[430]  = "a:x,y,z:3;"
	WyckoffSyms[431]  = "a:x,y,z:3;"
	WyckoffSyms[432]  = "a:0,0,z:3;b:x,y,z:9;"
	WyckoffSyms[433]  = "a:z,z,z:1;b:x+z,-x+y+z,-y+z:3;"
	WyckoffSyms[434]  = "a:0,0,0:1;b:0,0,1/2:1;c:0,0,z:2;d:1/3,2/3,z:2;e:1/2,0,0:3;f:1/2,0,1/2:3;g:x,y,z:6;"
	WyckoffSyms[435]  = "a:0,0,0:3;b:0,0,1/2:3;c:0,0,z:6;d:1/2,0,1/2:9;e:1/2,0,0:9;f:x,y,z:18;"
	WyckoffSyms[436]  = "a:0,0,0:1;b:1/2,1/2,1/2:1;c:z,z,z:2;d:1,0,1/2:3;e:1/2,1/2,0:3;f:x+z,-x+y+z,-y+z:6;"
	WyckoffSyms[437]  = "a:0,0,0:1;b:0,0,1/2:1;c:1/3,2/3,0:1;d:1/3,2/3,1/2:1;e:2/3,1/3,0:1;f:2/3,1/3,1/2:1;g:0,0,z:2;h:1/3,2/3,z:2;i:2/3,1/3,z:2;"
	WyckoffSyms[437] += "j:x,-x,0:3;k:x,-x,1/2:3;l:x,y,z:6;"
	WyckoffSyms[438]  = "a:0,0,0:1;b:0,0,1/2:1;c:0,0,z:2;d:1/3,2/3,z:2;e:x,0,0:3;f:x,0,1/2:3;g:x,y,z:6;"
	WyckoffSyms[439]  = "a:x,-x,1/3:3;b:x,-x,5/6:3;c:x,y,z:6;"
	WyckoffSyms[440]  = "a:x,0,1/3:3;b:x,0,5/6:3;c:x,y,z:6;"
	WyckoffSyms[441]  = "a:x,-x,2/3:3;b:x,-x,1/6:3;c:x,y,z:6;"
	WyckoffSyms[442]  = "a:x,0,2/3:3;b:x,0,1/6:3;c:x,y,z:6;"
	WyckoffSyms[443]  = "a:0,0,0:3;b:0,0,1/2:3;c:0,0,z:6;d:x,0,0:9;e:x,0,1/2:9;f:x,y,z:18;"
	WyckoffSyms[444]  = "a:0,0,0:1;b:1/2,1/2,1/2:1;c:z,z,z:2;d:x,-x,0:3;e:x+1/2,-x+1/2,1/2:3;f:x+z,-x+y+z,-y+z:6;"
	WyckoffSyms[445]  = "a:0,0,z:1;b:1/3,2/3,z:1;c:2/3,1/3,z:1;d:x,-x,z:3;e:x,y,z:6;"
	WyckoffSyms[446]  = "a:0,0,z:1;b:1/3,2/3,z:2;c:x,0,z:3;d:x,y,z:6;"
	WyckoffSyms[447]  = "a:0,0,z:2;b:1/3,2/3,z:2;c:2/3,1/3,z:2;d:x,y,z:6;"
	WyckoffSyms[448]  = "a:0,0,z:2;b:1/3,2/3,z:2;c:x,y,z:6;"
	WyckoffSyms[449]  = "a:0,0,z:3;b:x,-x,z:9;c:x,y,z:18;"
	WyckoffSyms[450]  = "a:z,z,z:1;b:x+z,-x*2+z,x+z:3;c:x+z,-x+y+z,-y+z:6;"
	WyckoffSyms[451]  = "a:0,0,z:6;b:x,y,z:18;"
	WyckoffSyms[452]  = "a:z,z,z:2;b:x+z,-x+y+z,-y+z:6;"
	WyckoffSyms[453]  = "a:0,0,0:1;b:0,0,1/2:1;c:1/3,2/3,0:2;d:1/3,2/3,1/2:2;e:0,0,z:2;f:1/2,0,0:3;g:1/2,0,1/2:3;h:1/3,2/3,z:4;i:x,-x,0:6;j:x,-x,1/2:6;"
	WyckoffSyms[453] += "k:x,0,z:6;l:x,y,z:12;"
	WyckoffSyms[454]  = "a:0,0,1/4:2;b:0,0,0:2;c:1/3,2/3,1/4:2;d:2/3,1/3,1/4:2;e:0,0,z:4;f:1/3,2/3,z:4;g:1/2,0,0:6;h:x,-x,1/4:6;i:x,y,z:12;"
	WyckoffSyms[455]  = "a:0,0,0:1;b:0,0,1/2:1;c:0,0,z:2;d:1/3,2/3,z:2;e:1/2,0,0:3;f:1/2,0,1/2:3;g:x,0,0:6;h:x,0,1/2:6;i:x,-x,z:6;j:x,y,z:12;"
	WyckoffSyms[456]  = "a:0,0,1/4:2;b:0,0,0:2;c:0,0,z:4;d:1/3,2/3,z:4;e:1/2,0,0:6;f:x,0,1/4:6;g:x,y,z:12;"
	WyckoffSyms[457]  = "a:0,0,0:3;b:0,0,1/2:3;c:0,0,z:6;d:1/2,0,1/2:9;e:1/2,0,0:9;f:x,0,0:18;g:x,0,1/2:18;h:x,-x,z:18;i:x,y,z:36;"
	WyckoffSyms[458]  = "a:0,0,0:1;b:1/2,1/2,1/2:1;c:z,z,z:2;d:1,0,1/2:3;e:1/2,1/2,0:3;f:x,-x,0:6;g:x+1/2,-x+1/2,1/2:6;h:x+z,-x*2+z,x+z:6;"
	WyckoffSyms[458] += "i:x+z,-x+y+z,-y+z:12;"
	WyckoffSyms[459]  = "a:0,0,1/4:6;b:0,0,0:6;c:0,0,z:12;d:1/2,0,0:18;e:x,0,1/4:18;f:x,y,z:36;"
	WyckoffSyms[460]  = "a:1/4,1/4,1/4:2;b:0,0,0:2;c:z,z,z:4;d:1/2,1/2,0:6;e:x+1/4,-x+1/4,1/4:6;f:x+z,-x+y+z,-y+z:12;"
	// Hexagonal SG[168,194]  SG_idNum [462,488]   (25 idNums)
	WyckoffSyms[461]  = "a:0,0,z:1;b:1/3,2/3,z:2;c:1/2,0,z:3;d:x,y,z:6;"
	WyckoffSyms[462]  = "a:x,y,z:6;"
	WyckoffSyms[463]  = "a:x,y,z:6;"
	WyckoffSyms[464]  = "a:0,0,z:3;b:1/2,1/2,z:3;c:x,y,z:6;"
	WyckoffSyms[465]  = "a:0,0,z:3;b:1/2,1/2,z:3;c:x,y,z:6;"
	WyckoffSyms[466]  = "a:0,0,z:2;b:1/3,2/3,z:2;c:x,y,z:6;"
	WyckoffSyms[467]  = "a:0,0,0:1;b:0,0,1/2:1;c:1/3,2/3,0:1;d:1/3,2/3,1/2:1;e:2/3,1/3,0:1;f:2/3,1/3,1/2:1;g:0,0,z:2;h:1/3,2/3,z:2;i:2/3,1/3,z:2;"
	WyckoffSyms[467] += "j:x,y,0:3;k:x,y,1/2:3;l:x,y,z:6;"
	WyckoffSyms[468]  = "a:0,0,0:1;b:0,0,1/2:1;c:1/3,2/3,0:2;d:1/3,2/3,1/2:2;e:0,0,z:2;f:1/2,0,0:3;g:1/2,0,1/2:3;h:1/3,2/3,z:4;i:1/2,0,z:6;j:x,y,0:6;"
	WyckoffSyms[468] += "k:x,y,1/2:6;l:x,y,z:12;"
	WyckoffSyms[469]  = "a:0,0,1/4:2;b:0,0,0:2;c:1/3,2/3,1/4:2;d:2/3,1/3,1/4:2;e:0,0,z:4;f:1/3,2/3,z:4;g:1/2,0,0:6;h:x,y,1/4:6;i:x,y,z:12;"
	WyckoffSyms[470]  = "a:0,0,0:1;b:0,0,1/2:1;c:1/3,2/3,0:2;d:1/3,2/3,1/2:2;e:0,0,z:2;f:1/2,0,0:3;g:1/2,0,1/2:3;h:1/3,2/3,z:4;i:1/2,0,z:6;j:x,0,0:6;"
	WyckoffSyms[470] += "k:x,0,1/2:6;l:x,-x,0:6;m:x,-x,1/2:6;n:x,y,z:12;"
	WyckoffSyms[471]  = "a:x,0,0:6;b:x,x,1/4:6;c:x,y,z:12;"
	WyckoffSyms[472]  = "a:x,0,0:6;b:x,x,3/4:6;c:x,y,z:12;"
	WyckoffSyms[473]  = "a:0,0,0:3;b:0,0,1/2:3;c:1/2,0,0:3;d:1/2,0,1/2:3;e:0,0,z:6;f:1/2,0,z:6;g:x,0,0:6;h:x,0,1/2:6;i:x,x,0:6;j:x,x,1/2:6;k:x,y,z:12;"
	WyckoffSyms[474]  = "a:0,0,0:3;b:0,0,1/2:3;c:1/2,0,0:3;d:1/2,0,1/2:3;e:0,0,z:6;f:1/2,0,z:6;g:x,0,0:6;h:x,0,1/2:6;i:x,x,0:6;j:x,x,1/2:6;k:x,y,z:12;"
	WyckoffSyms[475]  = "a:0,0,0:2;b:0,0,1/4:2;c:1/3,2/3,1/4:2;d:1/3,2/3,3/4:2;e:0,0,z:4;f:1/3,2/3,z:4;g:x,0,0:6;h:x,x,1/4:6;i:x,y,z:12;"
	WyckoffSyms[476]  = "a:0,0,z:1;b:1/3,2/3,z:2;c:1/2,0,z:3;d:x,0,z:6;e:x,-x,z:6;f:x,y,z:12;"
	WyckoffSyms[477]  = "a:0,0,z:2;b:1/3,2/3,z:4;c:1/2,0,z:6;d:x,y,z:12;"
	WyckoffSyms[478]  = "a:0,0,z:2;b:1/3,2/3,z:4;c:x,0,z:6;d:x,y,z:12;"
	WyckoffSyms[479]  = "a:0,0,z:2;b:1/3,2/3,z:2;c:x,-x,z:6;d:x,y,z:12;"
	WyckoffSyms[480]  = "a:0,0,0:1;b:0,0,1/2:1;c:1/3,2/3,0:1;d:1/3,2/3,1/2:1;e:2/3,1/3,0:1;f:2/3,1/3,1/2:1;g:0,0,z:2;h:1/3,2/3,z:2;i:2/3,1/3,z:2;"
	WyckoffSyms[480] += "j:x,-x,0:3;k:x,-x,1/2:3;l:x,y,0:6;m:x,y,1/2:6;n:x,-x,z:6;o:x,y,z:12;"
	WyckoffSyms[481]  = "a:0,0,0:2;b:0,0,1/4:2;c:1/3,2/3,0:2;d:1/3,2/3,1/4:2;e:2/3,1/3,0:2;f:2/3,1/3,1/4:2;g:0,0,z:4;h:1/3,2/3,z:4;i:2/3,1/3,z:4;"
	WyckoffSyms[481] += "j:x,-x,0:6;k:x,y,1/4:6;l:x,y,z:12;"
	WyckoffSyms[482]  = "a:0,0,0:1;b:0,0,1/2:1;c:1/3,2/3,0:2;d:1/3,2/3,1/2:2;e:0,0,z:2;f:x,0,0:3;g:x,0,1/2:3;h:1/3,2/3,z:4;i:x,0,z:6;j:x,y,0:6;"
	WyckoffSyms[482] += "k:x,y,1/2:6;l:x,y,z:12;"
	WyckoffSyms[483]  = "a:0,0,0:2;b:0,0,1/4:2;c:1/3,2/3,1/4:2;d:2/3,1/3,1/4:2;e:0,0,z:4;f:1/3,2/3,z:4;g:x,0,0:6;h:x,y,1/4:6;i:x,y,z:12;"
	WyckoffSyms[484]  = "a:0,0,0:1;b:0,0,1/2:1;c:1/3,2/3,0:2;d:1/3,2/3,1/2:2;e:0,0,z:2;f:1/2,0,0:3;g:1/2,0,1/2:3;h:1/3,2/3,z:4;i:1/2,0,z:6;j:x,0,0:6;"
	WyckoffSyms[484] += "k:x,0,1/2:6;l:x,x,0:6;m:x,x,1/2:6;n:x,0,z:12;o:x,x,z:12;p:x,y,0:12;q:x,y,1/2:12;r:x,y,z:24;"
	WyckoffSyms[485]  = "a:0,0,1/4:2;b:0,0,0:2;c:1/3,2/3,1/4:4;d:1/3,2/3,0:4;e:0,0,z:4;f:1/2,0,1/4:6;g:1/2,0,0:6;h:1/3,2/3,z:8;i:1/2,0,z:12;j:x,0,1/4:12;"
	WyckoffSyms[485] += "k:x,x,1/4:12;l:x,y,0:12;m:x,y,z:24;"
	WyckoffSyms[486]  = "a:0,0,1/4:2;b:0,0,0:2;c:1/3,2/3,1/4:4;d:1/3,2/3,0:4;e:0,0,z:4;f:1/2,0,0:6;g:x,0,1/4:6;h:1/3,2/3,z:8;i:x,x,0:12;j:x,y,1/4:12;"
	WyckoffSyms[486] += "k:x,0,z:12;l:x,y,z:24;"
	WyckoffSyms[487]  = "a:0,0,0:2;b:0,0,1/4:2;c:1/3,2/3,1/4:2;d:1/3,2/3,3/4:2;e:0,0,z:4;f:1/3,2/3,z:4;g:1/2,0,0:6;h:x,x,1/4:6;i:x,0,0:12;j:x,y,1/4:12;"
	WyckoffSyms[487] += "k:x,x,z:12;l:x,y,z:24;"
	// Cubic SG[195,230]  SG_idNum [489,530]   (42 idNums)
	WyckoffSyms[488]  = "a:0,0,0:1;b:1/2,1/2,1/2:1;c:0,1/2,1/2:3;d:1/2,0,0:3;e:x,x,x:4;f:x,0,0:6;g:x,0,1/2:6;h:x,1/2,0:6;i:x,1/2,1/2:6;j:x,y,z:12;"
	WyckoffSyms[489]  = "a:0,0,0:4;b:1/2,1/2,1/2:4;c:1/4,1/4,1/4:4;d:3/4,3/4,3/4:4;e:x,x,x:16;f:x,0,0:24;g:x,1/4,1/4:24;h:x,y,z:48;"
	WyckoffSyms[490]  = "a:0,0,0:2;b:0,1/2,1/2:6;c:x,x,x:8;d:x,0,0:12;e:x,1/2,0:12;f:x,y,z:24;"
	WyckoffSyms[491]  = "a:x,x,x:4;b:x,y,z:12;"
	WyckoffSyms[492]  = "a:x,x,x:8;b:x,0,1/4:12;c:x,y,z:24;"
	WyckoffSyms[493]  = "a:0,0,0:1;b:1/2,1/2,1/2:1;c:0,1/2,1/2:3;d:1/2,0,0:3;e:x,0,0:6;f:x,0,1/2:6;g:x,1/2,0:6;h:x,1/2,1/2:6;i:x,x,x:8;j:0,y,z:12;"
	WyckoffSyms[493] += "k:1/2,y,z:12;l:x,y,z:24;"
	WyckoffSyms[494]  = "a:1/4,1/4,1/4:2;b:0,0,0:4;c:1/2,1/2,1/2:4;d:1/4,3/4,3/4:6;e:x,x,x:8;f:x,1/4,1/4:12;g:x,3/4,1/4:12;h:x,y,z:24;"
	WyckoffSyms[495]  = "a:1/4,1/4,1/4:2;b:0,0,0:4;c:1/2,1/2,1/2:4;d:1/4,3/4,3/4:6;e:x,x,x:8;f:x,1/4,1/4:12;g:x,3/4,1/4:12;h:x,y,z:24;"
	WyckoffSyms[496]  = "a:0,0,0:4;b:1/2,1/2,1/2:4;c:1/4,1/4,1/4:8;d:0,1/4,1/4:24;e:x,0,0:24;f:x,x,x:32;g:x,1/4,1/4:48;h:0,y,z:48;i:x,y,z:96;"
	WyckoffSyms[497]  = "a:1/8,1/8,1/8:8;b:2/3,2/3,2/3:8;c:0,0,0:16;d:1/2,1/2,1/2:16;e:x,x,x:32;f:x,1/8,1/8:48;g:x,y,z:96;"
	WyckoffSyms[498]  = "a:1/8,1/8,1/8:8;b:2/3,2/3,2/3:8;c:0,0,0:16;d:1/2,1/2,1/2:16;e:x,x,x:32;f:x,1/8,1/8:48;g:x,y,z:96;"
	WyckoffSyms[499]  = "a:0,0,0:2;b:0,1/2,1/2:6;c:1/4,1/4,1/4:8;d:x,0,0:12;e:x,0,1/2:12;f:x,x,x:16;g:0,y,z:24;h:x,y,z:48;"
	WyckoffSyms[500]  = "a:0,0,0:4;b:1/2,1/2,1/2:4;c:x,x,x:8;d:x,y,z:24;"
	WyckoffSyms[501]  = "a:0,0,0:8;b:1/4,1/4,1/4:8;c:x,x,x:16;d:x,0,1/4:24;e:x,y,z:48;"
	WyckoffSyms[502]  = "a:0,0,0:1;b:1/2,1/2,1/2:1;c:0,1/2,1/2:3;d:1/2,0,0:3;e:x,0,0:6;f:x,1/2,1/2:6;g:x,x,x:8;h:x,1/2,0:12;i:0,y,y:12;j:1/2,y,y:12;"
	WyckoffSyms[502] += "k:x,y,z:24;"
	WyckoffSyms[503]  = "a:0,0,0:2;b:1/4,1/4,1/4:4;c:3/4,3/4,3/4:4;d:0,1/2,1/2:6;e:1/4,0,1/2:6;f:1/4,1/2,0:6;g:x,x,x:8;h:x,0,0:12;i:x,0,1/2:12;"
	WyckoffSyms[503] += "j:x,1/2,0:12;k:1/4,y,-y+1/2:12;l:1/4,y,y+1/2:12;m:x,y,z:24;"
	WyckoffSyms[504]  = "a:0,0,0:4;b:1/2,1/2,1/2:4;c:1/4,1/4,1/4:8;d:0,1/4,1/4:24;e:x,0,0:24;f:x,x,x:32;g:0,y,y:48;h:1/2,y,y:48;i:x,1/4,1/4:48;"
	WyckoffSyms[504] += "j:x,y,z:96;"
	WyckoffSyms[505]  = "a:0,0,0:8;b:1/2,1/2,1/2:8;c:1/8,1/8,1/8:16;d:2/3,2/3,2/3:16;e:x,x,x:32;f:x,0,0:48;g:1/8,y,-y+1/4:48;h:x,y,z:96;"
	WyckoffSyms[506]  = "a:0,0,0:2;b:0,1/2,1/2:6;c:1/4,1/4,1/4:8;d:1/4,1/2,0:12;e:x,0,0:12;f:x,x,x:16;g:x,1/2,0:24;h:0,y,y:24;i:1/4,y,-y+1/2:24;"
	WyckoffSyms[506] += "j:x,y,z:48;"
	WyckoffSyms[507]  = "a:1/8,1/8,1/8:4;b:2/3,2/3,2/3:4;c:x,x,x:8;d:1/8,y,-y+1/4:12;e:x,y,z:24;"
	WyckoffSyms[508]  = "a:1/3,1/3,1/3:4;b:7/8,7/8,7/8:4;c:x,x,x:8;d:1/8,y,y+1/4:12;e:x,y,z:24;"
	WyckoffSyms[509]  = "a:1/8,1/8,1/8:8;b:7/8,7/8,7/8:8;c:1/8,0,1/4:12;d:2/3,0,1/4:12;e:x,x,x:16;f:x,0,1/4:24;g:1/8,y,y+1/4:24;h:1/8,y,-y+1/4:24;"
	WyckoffSyms[509] += "i:x,y,z:48;"
	WyckoffSyms[510]  = "a:0,0,0:1;b:1/2,1/2,1/2:1;c:0,1/2,1/2:3;d:1/2,0,0:3;e:x,x,x:4;f:x,0,0:6;g:x,1/2,1/2:6;h:x,1/2,0:12;i:x,x,z:12;j:x,y,z:24;"
	WyckoffSyms[511]  = "a:0,0,0:4;b:1/2,1/2,1/2:4;c:1/4,1/4,1/4:4;d:3/4,3/4,3/4:4;e:x,x,x:16;f:x,0,0:24;g:x,1/4,1/4:24;h:x,x,z:48;i:x,y,z:96;"
	WyckoffSyms[512]  = "a:0,0,0:2;b:0,1/2,1/2:6;c:x,x,x:8;d:1/4,1/2,0:12;e:x,0,0:12;f:x,1/2,0:24;g:x,x,z:24;h:x,y,z:48;"
	WyckoffSyms[513]  = "a:0,0,0:2;b:0,1/2,1/2:6;c:1/4,1/2,0:6;d:1/4,0,1/2:6;e:x,x,x:8;f:x,0,0:12;g:x,1/2,0:12;h:x,0,1/2:12;i:x,y,z:24;"
	WyckoffSyms[514]  = "a:0,0,0:8;b:1/4,1/4,1/4:8;c:0,1/4,1/4:24;d:1/4,0,0:24;e:x,x,x:32;f:x,0,0:48;g:x,1/4,1/4:48;h:x,y,z:96;"
	WyckoffSyms[515]  = "a:1/3,0,1/4:12;b:7/8,0,1/4:12;c:x,x,x:16;d:x,0,1/4:24;e:x,y,z:48;"
	WyckoffSyms[516]  = "a:0,0,0:1;b:1/2,1/2,1/2:1;c:0,1/2,1/2:3;d:1/2,0,0:3;e:x,0,0:6;f:x,1/2,1/2:6;g:x,x,x:8;h:x,1/2,0:12;i:0,y,y:12;j:1/2,y,y:12;"
	WyckoffSyms[516] += "k:0,y,z:24;l:1/2,y,z:24;m:x,x,z:24;n:x,y,z:48;"
	WyckoffSyms[517]  = "a:1/4,1/4,1/4:2;b:3/4,1/4,1/4:6;c:0,0,0:8;d:0,3/4,1/4:12;e:x,1/4,1/4:12;f:x,x,x:16;g:x,3/4,1/4:24;h:1/4,y,y:24;i:x,y,z:48;"
	WyckoffSyms[518]  = "a:1/4,1/4,1/4:2;b:3/4,1/4,1/4:6;c:0,0,0:8;d:0,3/4,1/4:12;e:x,1/4,1/4:12;f:x,x,x:16;g:x,3/4,1/4:24;h:1/4,y,y:24;i:x,y,z:48;"
	WyckoffSyms[519]  = "a:0,0,0:2;b:0,1/2,1/2:6;c:1/4,0,1/2:6;d:1/4,1/2,0:6;e:1/4,1/4,1/4:8;f:x,0,0:12;g:x,0,1/2:12;h:x,1/2,0:12;i:x,x,x:16;"
	WyckoffSyms[519] += "j:1/4,y,y+1/2:24;k:0,y,z:24;l:x,y,z:48;"
	WyckoffSyms[520]  = "a:1/4,1/4,1/4:2;b:0,0,0:4;c:1/2,1/2,1/2:4;d:1/4,3/4,3/4:6;e:x,x,x:8;f:1/2,1/4,3/4:12;g:x,1/4,1/4:12;h:x,1/4,3/4:24;"
	WyckoffSyms[520] += "i:1/2,y,y+1/2:24;j:1/2,y,-y:24;k:x,x,z:24;l:x,y,z:48;"
	WyckoffSyms[521]  = "a:1/4,1/4,1/4:2;b:0,0,0:4;c:1/2,1/2,1/2:4;d:1/4,3/4,3/4:6;e:x,x,x:8;f:1/2,1/4,3/4:12;g:x,1/4,1/4:12;h:x,1/4,3/4:24;"
	WyckoffSyms[521] += "i:1/2,y,y+1/2:24;j:1/2,y,-y:24;k:x,x,z:24;l:x,y,z:48;"
	WyckoffSyms[522]  = "a:0,0,0:4;b:1/2,1/2,1/2:4;c:1/4,1/4,1/4:8;d:0,1/4,1/4:24;e:x,0,0:24;f:x,x,x:32;g:x,1/4,1/4:48;h:0,y,y:48;i:1/2,y,y:48;"
	WyckoffSyms[522] += "j:0,y,z:96;k:x,x,z:96;l:x,y,z:192;"
	WyckoffSyms[523]  = "a:1/4,1/4,1/4:8;b:0,0,0:8;c:1/4,0,0:24;d:0,1/4,1/4:24;e:x,0,0:48;f:x,1/4,1/4:48;g:x,x,x:64;h:1/4,y,y:96;i:0,y,z:96;j:x,y,z:192;"
	WyckoffSyms[524]  = "a:1/8,1/8,1/8:8;b:1/3,1/3,1/3:8;c:0,0,0:16;d:1/2,1/2,1/2:16;e:x,x,x:32;f:x,1/8,1/8:48;g:x,x,z:96;h:0,y,-y:96;i:x,y,z:192;"
	WyckoffSyms[525]  = "a:1/8,1/8,1/8:8;b:1/3,1/3,1/3:8;c:0,0,0:16;d:1/2,1/2,1/2:16;e:x,x,x:32;f:x,1/8,1/8:48;g:x,x,z:96;h:0,y,-y:96;i:x,y,z:192;"
	WyckoffSyms[526]  = "a:1/8,1/8,1/8:16;b:1/4,1/4,1/4:32;c:0,0,0:32;d:7/8,1/8,1/8:48;e:x,x,x:64;f:x,1/8,1/8:96;g:1/4,y,-y:96;h:x,y,z:192;"
	WyckoffSyms[527]  = "a:1/8,1/8,1/8:16;b:1/4,1/4,1/4:32;c:0,0,0:32;d:7/8,1/8,1/8:48;e:x,x,x:64;f:x,1/8,1/8:96;g:1/4,y,-y:96;h:x,y,z:192;"
	WyckoffSyms[528]  = "a:0,0,0:2;b:0,1/2,1/2:6;c:1/4,1/4,1/4:8;d:1/4,0,1/2:12;e:x,0,0:12;f:x,x,x:16;g:x,0,1/2:24;h:0,y,y:24;i:1/4,y,-y+1/2:48;"
	WyckoffSyms[528] += "j:0,y,z:48;k:x,x,z:48;l:x,y,z:96;"
	WyckoffSyms[529]  = "a:0,0,0:16;b:1/8,1/8,1/8:16;c:1/8,0,1/4:24;d:1/3,0,1/4:24;e:x,x,x:32;f:x,0,1/4:48;g:1/8,y,-y+1/4:48;h:x,y,z:96;"

	String list=WyckoffSyms[SG_idNum-1]
	Make/N=(ItemsInList(list),3)/FREE/T WyckList
	WyckList[][0] = StringFromList(0,StringFromList(p,list),":")
	WyckList[][1] = StringFromList(1,StringFromList(p,list),":")
	WyckList[][2] = StringFromList(2,StringFromList(p,list),":")
	SetDimLabel 1,0,letter,WyckList
	SetDimLabel 1,1,symOp,WyckList
	SetDimLabel 1,2,mult,WyckList

	return WyckList
End
//	Function test_GetWyckoffSymStrings(SpaceGroupID)
//		String SpaceGroupID
//		Wave WyckList=LatticeSym#GetWyckoffSymStrings230(SpaceGroupID)
//		// print WyckList
//		Duplicate/O WyckList, WyckListView
//		DisplayTableOfWave(WyckListView)
//	End
