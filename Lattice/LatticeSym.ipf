#pragma TextEncoding = "UTF-8"		// For details execute DisplayHelpTopic "The TextEncoding Pragma"
#pragma ModuleName=LatticeSym
#pragma version = 5.18
#include "Utility_JZT" version>=3.78
#include "xtl_Locate"										// used to find the path to the materials files (only contains CrystalsAreHere() )

// #define 	OLD_LATTICE_ORIENTATION					// used to get old direct lattice orientation (pre version 5.00)

Static strConstant NEW_LINE="\n"						//	was NL="\r"
Static Constant minPossibleBondLength = 0.050		// 0.050 nm = 50 pm, minimum possible distance between atoms (smallest known bond is 74 pm)
Static Constant ELEMENT_Zmax = 116

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
// with version 4.15, more changes along with vers 4.14, changed: setLattice(), FillLatticeParametersPanel(), readFileXML(),
//		LatticePanelButtonProc(), readCrystalStructure_xtl(), readCrystalStructureXML(), CIF_interpret(), read_cri_fileOLD()
// with version 4.16, added OverOccupyList() and OverOccupyList2Str(), and changed other routines to print a warning when occupy>1
//		changed: print_crystalStructure(), readCrystalStructure90, LatticePanelButtonProc()
// with version 4.17, changed MinimalChemFormula(), it better removes common factors
// with version 4.18, fixed an ERROR in Fstruct, it now properly checks for valid xtal.SpaceGroup
//		also added recipFrom_xtal(), a little utility that make a free recip matrix from xtal
// with version 4.19, small change in reMakeAtomXYZs(), avoids error when debug on "NVAR SVAR WAVE Checking" is on
// with version 4.20, the xtl file is now carried along in the crystalStructure structure
// with version 4.21, added directFrom_xtal(), which is a lot like recipFrom_xtal()
// with version 4.23, added str2hkl(hklStr,h,k,l)
// with version 4.24, moved ConvertUnits2meters(), ConvertTemperatureUnits() , and placesOfPrecision() to Utility_JZT.ipf
//		also added NearestAllowedHKL(xtal,hkl)
// with version 4.25, added equivalentHKLs(xtal,hkl0,[noNeg]), get list of hkl's symmetry equivalent to hkl0
// with version 4.26, added bond info to the AtomView menu
// with version 4.27, when writing files, use "\n" instead of "\r" for line termination. 
//		also changed crystalStructure2xml() and convertOneXTL2XML() to take a new line argument.
// with version 4.28, in positionsOfOneAtomType() duplicate atoms are now within 50pm (in x,y,&z) not just 1e-3
//		also positionsOfOneAtomType() uses free waves rather than real temp waves
// with version 4.30, in str2recip(str), now handles both "}{" and "},{" type strings
//		also added some helpful comments about convert recip <--> direct using MatrixOP
// with version 4.31, added wave note info in directFrom_xtal(xtal) and in recipFrom_xtal(xtal)
// with version 4.32, added direct2LatticeConstants(direct)
// with version 4.33, added isValidLatticeConstants(direct), and imporved formatting in print_crystalStructure()
// with version 4.34, added DescribeSymOps(direct), CheckDirectRecip()
// with version 4.35, MatrixOP in Fstruct(), change variable Q --> Qmag, and calc of F uses complex math now.
//		also in positionsOfOneAtomType, make condition based on true atom distances, & use MatrixOP too.
// with version 4.37, changed allowedHKL(), now a ThreadSafe version of allowedHKL()
//		allowedHKL() now cannot use Cromer, but it is ThreadSafe (done for indexing routines).
//		to use allowedHKL() with MultiThread you MUST pass the wave refs to the atomN waves
//		also added isValidSpaceGroup()
// with version 4.38, added ELEMENT_Symbols and changed Get_f_proto() to be a better.
// with version 4.39, Fixed ERROR in Fstruct()
// with version 4.40, Get_f_proto() was ThreadSafe, it must NOT be ThreadSafe (since Get_F is not)
// with version 4.41, removed ELEMENT_Symbols (it is now in Utility_JZT.ipf)
// with version 4.42, fixed problem with Rhombohedal lattices & GetLatticeConstants()
// with version 4.43, added showPanel to InitLatticeSymPackage()
// with version 4.44, definition of MenuItemIfWindowAbsent() was changed
// with version 4.45, added SetToDummyXTAL() for use by FillCrystalStructDefault() for startup bug found by Jan
// with version 4.46, added SetToDummyATOM() and cleaned up v4.45
// with version 4.48, fixed a bug in positionsOfOneAtomType() & equivalentHKLs() occured when rowRepeat(vec,1)
//
// with version 5.00, changed definition of direct lattice to a||x and c*|| z (was b||y), added #define OLD_LATTICE_ORIENTATION to get old way
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

// Rhombohedral Transformation:
//
// for Rhombohedral (hkl), and Hexagonal (HKL)
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
// Although not used here, note that the following also works:
//		NOTE Inv(A^t) = Inv(A)^t, so the order of Inv() and ^t do not matter
//		MatrixOP recipLatice   = 2*PI * (Inv(directLattice))^t
//		MatrixOP directLattice = 2*PI * Inv(recipLatice^t)
//		Vc      = MatrixDet(directLattice)		// Volume of unit cell
//		VcRecip = MatrixDet(recipLatice)		// Volume of reciprocal lattice cell
//
//		kf^ = ki^ - 2*(ki^ . q^)*q^		note: (ki^ . q^) must be NEGATIVE (both Bragg & Laue)


Menu "Analysis"
	SubMenu "Lattice"
		"<BLattice Set Panel...",MakeLatticeParametersPanel("")
		help={"Define the crystal structure and lattice to be used"}
		"Show Crystal Structure",showCrystalStructure()
		help={"Shows the crystal structure and lattice that is currently defined"}
		"  Edit the Atom Positions...",EditAtomPositionsMenu()
		help={"Manually set/change the atom positions"}
		"Load a new Crystal Structure...",LoadCrystal("")
		help={"load a rystal structure from a fie"}
		"Write Current Crystal to XML File...",writeCrystalStructure2xmlFile("","")
		help={"takes the current crystal structure and writes it to an xml file"}
		"d[hkl]",/Q,get_dhkl(NaN,NaN,NaN,T=NaN)
		help={"show d-spacing for the hkl reflection"}
		"Fstructure [hkl]", /Q ,getFstruct(NaN,NaN,NaN)
		help={"Crude Structure Factor using current lattice"}
		"Find Closest hkl from d-spacing or Q",findClosestHKL(NaN)
		help={"Knowing either the d-spacing or the Q, find closest hkl's"}
		"\\M0Space Group number <--> symmetry",symmtry2SG("")
		help={"find the Space Group number from symmetry string,  e.g. Pmma, or sym from number"}
		"Describe the Symmetry Operations", DescribeSymOps($"")
		"angle between two hkl's",angleBetweenHKLs(NaN,NaN,NaN,  NaN,NaN,NaN)
		"  Convert old xtl files to new xml files",ConverXTLfile2XMLfile("")
		"-"
		MenuItemIfWindowAbsent("Include Powder Patterns Support","PowderPatterns.ipf","WIN:128"), Execute/P "INSERTINCLUDE  \"PowderPatterns\", version>=0.10";Execute/P "COMPILEPROCEDURES ";Execute/P "Init_PowderPatternLattice()"
		help = {"Load Function used to compute Powder Patterns from Loaded Lattice"}
		MenuItemIfWindowAbsent("Include Atom View Support","AtomView.ipf","WIN:128"), Execute/P "INSERTINCLUDE  \"AtomView\", version>=0.17";Execute/P "COMPILEPROCEDURES ";Execute/P "Init_AtomViewLattice()"
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

#if (IgorVersion()<7)
	strConstant DEGREESIGN = "\241"			// option-shift-8
	strConstant BULLET = "\245"				// option-8
	strConstant ARING = "\201"				// Angstrom sign, option-shift-A
#if StringMatch(IgorInfo(2),"Windows")
	strConstant BCHAR = "\257"
#else
	strConstant BCHAR = "\321"
#endif
#else
	strConstant DEGREESIGN = "\xC2\xB0"	// UTF8, DEGREE SIGN
	strConstant BULLET = "\xE2\x80\xA2"
	strConstant ARING = "\xC3\x85"			// Aring, Angstrom sign
	strConstant BCHAR = "\xE2\x80\x94"		// EM DASH
#endif

// =========================================================================
// =========================================================================
//	Start of Structure definitions

Structure crystalStructure	// structure definition for a crystal lattice
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
	NewPanel /K=1 /W=(675,60,675+220,60+508)
	DoWindow/C LatticeSet
	FillLatticeParametersPanel(strStruct,"LatticeSet",0,0)
End


Function showCrystalStructure()						// prints the structure that is currently being used
	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))					//fill the lattice structure with test values
		DoAlert 0, "no crystal structure found"
		return 1
	endif
	String str, sym = getHMboth(xtal.SpaceGroup)
	sprintf str, "'%s'  %d atoms,  #%d   %s\r%.9g, %.9g, %.9gnm,\r%g%s, %g%s, %g%s",xtal.desc,xtal.N,xtal.SpaceGroup,sym,xtal.a,xtal.b,xtal.c,xtal.alpha,DEGREESIGN,xtal.beta,DEGREESIGN,xtal.gam,DEGREESIGN
	Variable netCharge = NetChargeCell(xtal)
	if (netCharge)
		str += "\r  *** Charge imbalance in cell = "+num2str(netCharge)+" ***"
	endif
	DoAlert 1,str+"\r\tprint all details to history too?"
	if (V_flag==1)
		print_crystalStructure(xtal)					// prints out the value in a crystalStructure structure
	else
		printf "currently using  '%s'  lattice is  #%d   %s     %.9gnm, %.9gnm, %.9gnm,   %g%s, %g%s, %g%s",xtal.desc,xtal.SpaceGroup,sym,xtal.a,xtal.b,xtal.c,xtal.alpha,DEGREESIGN,xtal.beta,DEGREESIGN,xtal.gam,DEGREESIGN
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
	copy_xtal(xtal,xtal_IN)						// copy xtal_IN to a working copy

	SetSymOpsForSpaceGroup(xtal.SpaceGroup)		// ensure that symmetry ops are right
	Variable Natoms = round(xtal.N)				// number of predefined atoms
	Natoms = numtype(Natoms) ? 0 : limit(Natoms,0,STRUCTURE_ATOMS_MAX)
	String wyckMenuStr = "\" ;"+WyckoffMenuStr(xtal.SpaceGroup)+"\""

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
	PauseForUser $win										// WAIT here for Panel to close then proceed
	KillWaves/Z atomNameList

	// user is done entering atom type information, save it
	StructGet/S/B=2 xtal, strStruct						// retrieve xtal from the global string
	KillStrings/Z strStruct									// done with this global string
	if (LatticeBad(xtal,atomsToo=1))
		DoAlert 0, "ERROR -- Invalid Lattice"
		return NaN
	endif

	DoAlert 1,"Entered Values are Valid,\rRepalce Curret Values with These?"
	if (V_flag!=1)
		return NaN
	endif
	xtal.Vibrate = xtalVibrates(xtal)						// True if some Thermal vibration info present in xtal
	xtal.haveDebyeT = xtalHasDebye(xtal)					// True if some one of the atoms has a Debye Temperature
	xtal.hashID = xtalHashID(xtal)
	copy_xtal(xtal_IN,xtal)								// copy working copy to the passed copy
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
	StructGet/S/B=2 xtal, strStruct
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
		StructGet/S/B=2 xtal, strStruct
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
	StructGet/S/B=2 xtal, strStruct

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
		if (!ForceXYZtoWyckoff(xtal.SpaceGroup,symbol,xx,yy,zz))
			SetVariable setAtomRelX,win=$win,value= _NUM:xx
			SetVariable setAtomRelY,win=$win,value= _NUM:yy
			SetVariable setAtomRelZ,win=$win,value= _NUM:zz
		endif
	elseif (StringMatch(ctrlName,"setAtomRel*") && strlen(symbol)==0)	// changed xyz, find Wyckoff symbol
		Variable mm
		symbol = FindWyckoffSymbol(xtal.SpaceGroup,xx,yy,zz,mult)
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
		mult = WyckoffMultiplicity(xtal.SpaceGroup,symbol)
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
	StructGet/S/B=2 xtal, strStruct
	print_crystalStructure(xtal)
End


Function print_crystalStructure(xtal)			// prints out the value in a crystalStructure
	STRUCT crystalStructure &xtal				// this sruct is printed in this routine
	Variable i,N=xtal.N
	if (WhichListItem("LatticePanelButtonProc",GetRTStackInfo(0))>=0)
		printf "%sshowCrystalStructure()\r",BULLET
	endif
	if (strlen(xtal.desc))
		printf "'%s'   \t",xtal.desc
	else
		printf "\t\t\t"
	endif
	Variable SG = xtal.SpaceGroup
	printf "Space Group=%d   %s   %s       Vc = %g (nm^3)",SG,StringFromList(latticeSystem(SG),LatticeSystemNames), getHMboth(SG),xtal.Vc
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
	copy_xtal(xtal,xtalIN)					// not changing xtalIN

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
	T_K = ParamIsDefault(T_K) ? 0 : T_K			// only consider T_K if it is passed (0 is valid)
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
	if (thetaM>0 && T_K>=0)		// have a valid and Debye Temperature and Temperature
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

	Variable SpaceGroup=xtal.SpaceGroup			// local value for convienence
	if (!isValidSpaceGroup(SpaceGroup))			// invalid SpaceGroup, it must be in range [1,230]
		DoAlert 0, "invalid Space Group number "+num2str(SpaceGroup)
		return 1
	endif
	//	Cubic				[195,230]		//	a
	//	Hexagonal		[168,194]			//	a,c
	//	Trigonal			[143,167]		//	a,alpha
	//	Tetragonal		[75,142]		//	a,c
	//	Orthorhombic	[16,74]			//	a,b,c
	//	Monoclinic		[3,15]		//	a,b,c,gamma
	//	Triclinic		[1,2]				//	a,b,c,alpha,beta,gamma

	if (SpaceGroup>=195)				// Cubic
		xtal.b = xtal.a
		xtal.c = xtal.a
		xtal.alpha=90;  xtal.beta=90;  xtal.gam=90
	elseif(SpaceGroup>=168)			// Hexagonal
		xtal.b = xtal.a
		xtal.alpha=90;  xtal.beta=90;  xtal.gam=120

	elseif(SpaceGroup>=143)			// Trigonal (generally hexagonal cell), for rhomohedral use rhomohedral cell, unless obviously the hexagonal
		if (isRhombohedral(SpaceGroup))	// rhombohedral structure
			if ((abs(90-xtal.alpha)+abs(90-xtal.beta)+abs(120-xtal.gam))<1e-3)	// obviously hexagonal
				xtal.b = xtal.a
				xtal.alpha=90  ;  xtal.beta=90  ;  xtal.gam=120
			else								// rhombohedral with rhombohedral cell
				xtal.b = xtal.a
				xtal.c = xtal.a
				xtal.beta = xtal.alpha
				xtal.gam = xtal.alpha
			endif
		else									// "P" trigonal, use hexagonal cell
			xtal.b = xtal.a
			xtal.alpha=90  ;  xtal.beta=90  ;  xtal.gam=120
		endif

	elseif(SpaceGroup>=75)				// Tetragonal
		xtal.b = xtal.a
		xtal.alpha=90;  xtal.beta=90;  xtal.gam=90
	elseif(SpaceGroup>=16)				// Orthorhombic
		xtal.alpha=90;  xtal.beta=90;  xtal.gam=90
	elseif(SpaceGroup>=3)				// Monoclinic
//		xtal.alpha=90;  xtal.beta=90
		xtal.alpha=90;  xtal.gam=90
//	else										// Triclinic
	endif

	String str
	if (xtal.a<=0 || xtal.b<=0 || xtal.c<=0 || numtype(xtal.a+xtal.b+xtal.c))	// check for valid a,b,c
		sprintf str,"invalid, (a,b,c) = (%g,%g,%g)",xtal.a,xtal.b,xtal.c
		DoAlert 0, str
		return 1
	endif
	Variable alpha=xtal.alpha,bet=xtal.beta,gam=xtal.gam
	if (alpha<=0 || bet<=0 || gam<=0 || alpha>=180 || bet>=180 || gam>=180 || numtype(alpha+bet+gam))	// check for valid angles
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
	STRUCT crystalStructure &xtal						// this sruct is filled  by this routine
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

Static Function NetChargeCell(xtal)						// find the net charge in a cell (from valences), should be zero
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
	xtal.SpaceGroup=1
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

Function angleBetweenHKLs(h1,k1,l1,  h2,k2,l2)		// find angle between (h1,k1,l1) and (h2,k2,l2)
	Variable h1,k1,l1,h2,k2,l2
	if (numtype(h1+k1+l1+h2+k2+l2) && ItemsInList(GetRTStackInfo(0))<2 )
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
	if (ItemsInList(GetRTStackInfo(0))<2)
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
		StructGet/S/B=2 xtal, strStruct
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
	xtal.SpaceGroup=195									// simple cubic
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
	copy_xtal(xtal,xtalIN)
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
	Variable new = !(NumVarOrDefault("root:Packages:Lattices:PanelValues:SpaceGroup",-1)>0)
	String/G root:Packages:Lattices:PanelValues:desc		// create, but don't fill in the values for the panel
	Variable/G root:Packages:Lattices:PanelValues:SpaceGroup
	Variable/G root:Packages:Lattices:PanelValues:a
	Variable/G root:Packages:Lattices:PanelValues:b
	Variable/G root:Packages:Lattices:PanelValues:c
	Variable/G root:Packages:Lattices:PanelValues:alpha
	Variable/G root:Packages:Lattices:PanelValues:bet
	Variable/G root:Packages:Lattices:PanelValues:gam
	Variable/G root:Packages:Lattices:PanelValues:alphaT
	Variable/G root:Packages:Lattices:PanelValues:dirty		// =1 (xtal.sourceFile bad too),  =2 (values changed, but xtal.sourceFile is OK)
	Variable/G root:Packages:Lattices:PanelValues:T_C
	NVAR SpaceGroup=root:Packages:Lattices:PanelValues:SpaceGroup
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
		StructGet/S/B=2 xtal, strStruct				// found passed structure information, load into xtal
		CleanOutCrystalStructure(xtal)
		dirty = 2											// xtal.sourceFile should be OK
		crystalStructStr = strStruct
	elseif(new)												// no old values present, use usual defaults
		FillCrystalStructDefault(xtal)
		a=xtal.a  ;  b=xtal.b  ;  c=xtal.c
		alpha=xtal.alpha  ;  bet=xtal.beta  ;  gam=xtal.gam
		SpaceGroup=xtal.SpaceGroup
		alphaT=xtal.alphaT
		desc=xtal.desc
		T_C = xtal.Temperature
		T_C = numtype(T_C) || T_C<-273.14 ? 22.5 : T_C
		T_C = (xtal.haveDebyeT) ? T_C : NaN		// if no Debye Temperatuers, no temperature needed
		dirty = 0
//		StructPut/S xtal, strStruct
		StructPut/S xtal, crystalStructStr
	endif

	SetWindow kwTopWin,userdata(LatticePanelName)=hostWin+"#LatticePanel"
	NewPanel/K=1/W=(left,top,left+221,top+508)/HOST=$hostWin
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
	SetVariable set_gamma, fSize=12, title="\xCE\xB3"+DEGREESIGN
#endif

	Button buttonLatticeSave,pos={35,233},size={150,20},proc=LatticePanelButtonProc,title="Save"
	Button buttonLatticeSave,help={"Save these values as the current values"}
	Button buttonLatticeRevert,pos={35,258},size={150,20},proc=LatticePanelButtonProc,title="Revert"
	Button buttonLatticeRevert,help={"revert to the current default"}
	Button buttonEditAtomPositions,pos={35,283},size={150,20},proc=LatticePanelButtonProc,title="Edit Atom Positions"
	Button buttonEditAtomPositions,help={"Manually edit the atom positions"}
	Button buttonLatticeFromFile,pos={35,308},size={150,20},proc=LatticePanelButtonProc,title="Load from a file"
	Button buttonLatticeFromFile,help={"Fill the values in this panel from a file"}
	Button buttonPrintLattice,pos={35,333},size={150,20},proc=LatticePanelButtonProc,title="Print Lattice to History"
	Button buttonPrintLattice,help={"print lattice values to the history"}
	Button buttonWriteLattice,pos={35,358},size={150,20},proc=LatticePanelButtonProc,title="Write Lattice to File"
	Button buttonWriteLattice,help={"write current lattice values to an xml file"}

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
	SetVariable T_LatticeVar,limits={-273.151,inf,0},value=root:Packages:Lattices:PanelValues:T_C
	SetVariable T_LatticeVar,help={"Temperature (C) used when Thermal factors are given"}
	ValDisplay Fr_3atticeDisp,pos={20,483},size={80,17},title="F =",value= #"root:Packages:Lattices:PanelValues:Fr"
	ValDisplay Fr_3atticeDisp,font="Lucida Grande",fSize=12,format="%.3f",limits={0,0,0},barmisc={0,1000}
	ValDisplay Fr_3atticeDisp,help={"real part of F(hkl) calulated from the lattice"},frame=0
	ValDisplay Fi_3atticeDisp,pos={104,483},size={80,17},title="+i",value= #"root:Packages:Lattices:PanelValues:Fi"
	ValDisplay Fi_3atticeDisp,font="Lucida Grande",fSize=12,format="%.3f",limits={0,0,0},barmisc={0,1000}
	ValDisplay Fi_3atticeDisp,help={"imag part of F(hkl) calulated from the lattice"},frame=0

	String subWin = GetUserData("","","LatticePanelName")
	UpdatePanelLatticeConstControls(subWin,SpaceGroup)

	STRUCT WMSetVariableAction sva
	sva.win = subWin
	sva.eventCode = 2
	LatticePanelParamProc(sva)
	return "#LatticePanel"
End
//
Static Function UpdatePanelLatticeConstControls(subWin,SpaceGroup)
	// update a,b,c, alpha,beta,gamma in LatticeSet Panel, changes who is enabled, not the values
	String subWin
	Variable SpaceGroup

	SpaceGroup = round(SpaceGroup)
	if (!isValidSpaceGroup(SpaceGroup))
		return 1
	endif

	NVAR a=root:Packages:Lattices:PanelValues:a
	NVAR b=root:Packages:Lattices:PanelValues:b
	NVAR c=root:Packages:Lattices:PanelValues:c
	NVAR alpha=root:Packages:Lattices:PanelValues:alpha
	NVAR bet=root:Packages:Lattices:PanelValues:bet
	NVAR gam=root:Packages:Lattices:PanelValues:gam
	NVAR T_C=root:Packages:Lattices:PanelValues:T_C

	String titleStr="\\JC#"+num2istr(SpaceGroup)+" "
	if (SpaceGroup>=195)												// Cubic, a
		SetVariable set_a_nm,noedit=0,frame=1,win=$subWin	// enable
		SetVariable set_b_nm,noedit=1,frame=0,win=$subWin	// disable
		SetVariable set_c_nm,noedit=1,frame=0,win=$subWin
		SetVariable set_alpha,noedit=1,frame=0,win=$subWin
		SetVariable set_beta,noedit=1,frame=0,win=$subWin
		SetVariable set_gamma,noedit=1,frame=0,win=$subWin
		b=a  ;  c=a  ;  alpha=90 ; bet=90 ; gam=90
		titleStr += "Cubic"
	elseif (SpaceGroup>=168)									// Hexagonal, a, c
		SetVariable set_a_nm,noedit=0,frame=1,win=$subWin	// enable
		SetVariable set_c_nm,noedit=0,frame=1,win=$subWin
		SetVariable set_b_nm,noedit=1,frame=0,win=$subWin	// disable
		SetVariable set_alpha,noedit=1,frame=0,win=$subWin
		SetVariable set_beta,noedit=1,frame=0,win=$subWin
		SetVariable set_gamma,noedit=1,frame=0,win=$subWin
		b=a ; alpha=90 ; bet=90 ; gam=120
		titleStr += "Hexagonal"
	elseif (isRhombohedral(SpaceGroup) && !((abs(90-alpha)+abs(90-bet)+abs(120-gam))<1e-6))	// Rhombohedral, with rhombohedral cell
		SetVariable set_a_nm,noedit=0,frame=1,win=$subWin	// enable
		SetVariable set_alpha,noedit=0,frame=1,win=$subWin
		SetVariable set_b_nm,noedit=1,frame=0,win=$subWin	// disable
		SetVariable set_c_nm,noedit=1,frame=0,win=$subWin
		SetVariable set_beta,noedit=1,frame=0,win=$subWin
		SetVariable set_gamma,noedit=1,frame=0,win=$subWin
		b=a  ;  c=a
		bet=alpha  ;  gam=alpha
		titleStr += "Rhombohedral"
	elseif (SpaceGroup>=143)									// Trigonal, with hexagonal cell
		SetVariable set_a_nm,noedit=0,frame=1,win=$subWin	// enable
		SetVariable set_c_nm,noedit=0,frame=1,win=$subWin
		SetVariable set_b_nm,noedit=1,frame=0,win=$subWin	// disable
		SetVariable set_alpha,noedit=1,frame=0,win=$subWin
		SetVariable set_beta,noedit=1,frame=0,win=$subWin
		SetVariable set_gamma,noedit=1,frame=0,win=$subWin
		b=a
		alpha=90  ;  bet=90  ;  gam=120
		titleStr += "Trigonal"
	elseif (SpaceGroup>=75)									// Tetragonal, a, c
		SetVariable set_a_nm,noedit=0,frame=1,win=$subWin	// enable
		SetVariable set_c_nm,noedit=0,frame=1,win=$subWin
		SetVariable set_b_nm,noedit=1,frame=0,win=$subWin	// disable
		SetVariable set_alpha,noedit=1,frame=0,win=$subWin
		SetVariable set_beta,noedit=1,frame=0,win=$subWin
		SetVariable set_gamma,noedit=1,frame=0,win=$subWin
		b=a ; alpha=90 ; bet=90 ; gam=90
		titleStr += "Tetragonal"
	elseif (SpaceGroup>=16)									// Orthorhombic, a, b, c
		SetVariable set_a_nm,noedit=0,frame=1,win=$subWin	// enable
		SetVariable set_b_nm,noedit=0,frame=1,win=$subWin
		SetVariable set_c_nm,noedit=0,frame=1,win=$subWin
		SetVariable set_alpha,noedit=1,frame=0,win=$subWin	// disable
		SetVariable set_beta,noedit=1,frame=0,win=$subWin
		SetVariable set_gamma,noedit=1,frame=0,win=$subWin
		alpha=90 ; bet=90 ; gam=90
		titleStr += "Orthorhombic"
	elseif (SpaceGroup>=3)										// Monoclinic, a, b, c, gamma
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
	titleStr += "   \\F'Courier'"+getHMboth(SpaceGroup)
	titleStr = minus2bar(titleStr)								// change all minuses to a bar over following character
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
Static Function SelectNewSG(find)
	String find
	if (strlen(find)<1)
		find = StrVarOrDefault("root:Packages:Lattices:PanelValues:SpaceGroupSearch","")
		Prompt find, "Space Group Search, use * for wild card"
		DoPrompt "Search String", find
		if (V_flag)
			return NaN
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

	String list = symmtry2SG(find,type=0,printIt=0)
	String system, systemNames="Triclinic\t;Monoclinic\t;Orthorhombic;Tetragonal\t;Trigonal\t;Hexagonal\t;Cubic\t\t"
	Variable SG, Nlist=ItemsInList(list)
	for (i=0; i<Nlist; i+=1)
		SG = str2num(StringFromList(i,list))
		if (isValidSpaceGroup(SG))
			system = StringFromList(latticeSystem(SG),systemNames)
			//	sprintf str,"%d  %s  %s;", SG,system,getFullHMSym(SG)
			sprintf str, "%d  %s  %s  [%s];", SG,system,getFullHMSym(SG),getHallSymbol(SG)
			symList += str
		endif
	endfor
	Variable N=ItemsInList(symList)
	if (N<1)
		return NaN
	elseif (N==1)
		return str2num(symList)
	endif

	String sym
	Prompt sym,"Space Group",popup,symList
	DoPrompt "Space Group",sym
	if (V_flag)
		return NaN
	endif
	return str2num(sym)
End


Function LatticePanelButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	if (ba.eventCode != 2)			// not mouse up
		return 0
	endif

	String ctrlName=ba.ctrlName
	SVAR desc=root:Packages:Lattices:PanelValues:desc
	NVAR SpaceGroup=root:Packages:Lattices:PanelValues:SpaceGroup
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
	StructGet/S/B=2 xtal, crystalStructStr				// pre-load with current information
	CleanOutCrystalStructure(xtal)
	STRUCT WMSetVariableAction sva
	sva.win = ba.win
	sva.eventCode = 2

	if (stringmatch(ctrlName,"buttonLatticeRevert"))	// close window, do not save numbers
		FillCrystalStructDefault(xtal)					//fill the lattice structure with default values
		a=xtal.a  ;  b=xtal.b  ;  c=xtal.c
		alpha=xtal.alpha  ;  bet=xtal.beta  ;  gam=xtal.gam
		SpaceGroup=xtal.SpaceGroup
		alphaT = xtal.alphaT
		desc = xtal.desc
		T_C = xtal.Temperature
		dirty = 0
//		SetVariable set_SpaceGroupSearch, userdata(oldValue)=num2istr(SpaceGroup),win=$ GetUserData("","","LatticePanelName")
		UpdatePanelLatticeConstControls(ba.win,SpaceGroup)
		LatticePanelParamProc(sva)
		StructPut/S xtal, crystalStructStr
	elseif (stringmatch(ctrlName,"buttonLatticeSave") && dirty)
		xtal.a = a  ;  xtal.b = b  ;  xtal.c = c
		xtal.alpha = alpha  ;  xtal.beta = bet  ;  xtal.gam = gam
		xtal.SpaceGroup = SpaceGroup
		xtal.alphaT = alphaT
		xtal.Temperature = T_C
		xtal.desc = desc[0,99]								// desc is limited to 100 chars
		if (dirty==1)											// dirty==1, means xtal.sourceFile is bad
			xtal.sourceFile = ""
		endif
		ForceLatticeToStructure(xtal)
		UpdateCrystalStructureDefaults(xtal)
		dirty = 0
		UpdatePanelLatticeConstControls(ba.win,SpaceGroup)// change main copy of xtal
		LatticePanelParamProc(sva)
		StructPut/S xtal, crystalStructStr				// update the local copy too
		OverOccupyList(xtal,printIt=1)					// print notice if some sites have occ>1
	elseif (stringmatch(ctrlName,"buttonFindSpaceGroup"))
		Variable SG = SelectNewSG("")
		if (numtype(SG))
			return 0
		endif
		SpaceGroup = SG
		dirty = 1												// source file no longer valid
		UpdatePanelLatticeConstControls(ba.win,SpaceGroup)
		LatticePanelParamProc(sva)
	elseif (stringmatch(ctrlName,"buttonEditAtomPositions"))
		if (EditAtomPositions(xtal)>=0)					// xtal has been MODIFIED
			if (ForceLatticeToStructure(xtal))			// sets everything including remaking atom0,atom1,...
				return NaN
			endif
			xtal.hashID = xtalHashID(xtal)
			StructPut/S xtal, crystalStructStr			// update the local copy
			dirty = 2											// source file still mostly OK
			UpdatePanelLatticeConstControls(ba.win,SpaceGroup)
			LatticePanelParamProc(sva)
		endif
	elseif (stringmatch(ctrlName,"buttonLatticeFromFile"))
		if (readCrystalStructure(xtal,"",printIt=1))
			DoAlert 0,"nothing read"
		endif
		a = xtal.a  ;  b = xtal.b  ;  c = xtal.c
		alpha = xtal.alpha  ;  bet = xtal.beta  ;  gam = xtal.gam
		SpaceGroup = xtal.SpaceGroup
		desc = xtal.desc
		alphaT = xtal.alphaT
		T_C = xtal.Temperature
		T_C = (xtal.haveDebyeT && numtype(T_C)) ? 22.5 : T_C
		dirty = 2														// yes, it is dirty, but xtal.sourceFile is valid
		StructPut/S xtal, crystalStructStr
		UpdatePanelLatticeConstControls(ba.win,SpaceGroup)
		LatticePanelParamProc(sva)
	elseif (stringmatch(ctrlName,"buttonPrintLattice"))	// print shown lattice parameters to the history
		xtal.a = a  ;  xtal.b = b  ;  xtal.c = c
		xtal.alpha = alpha  ;  xtal.beta = bet  ;  xtal.gam = gam
		xtal.SpaceGroup = SpaceGroup
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
		sprintf cmd,"LatticeSym#UpdatePanelLatticeConstControls(\"%s\",%g)",ba.win,SpaceGroup
		Execute/P "INSERTINCLUDE \"AtomView\", version>=0.17"
		Execute/P "COMPILEPROCEDURES "
		Execute/P "Init_AtomViewLattice()"
		Execute/P cmd
	elseif (stringmatch(ctrlName,"buttonPowderPattern"))	// add support for PowderPatterns
		sprintf cmd,"LatticeSym#UpdatePanelLatticeConstControls(\"%s\",%g)",ba.win,SpaceGroup
		Execute/P "INSERTINCLUDE \"PowderPatterns\", version>=0.10"
		Execute/P "COMPILEPROCEDURES "
		Execute/P "Init_PowderPatternLattice()"
		Execute/P cmd
	endif
	return 0
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

	PathInfo materials
	if (strlen(S_path)<1)							// make it if it does not exist
		NewPath/Z materials, ParseFilePath(1,FunctionPath("CrystalsAreHere"),":",1,0)
	endif
	PathInfo materials
	if (strlen(S_path)<1)							// make it if it does not exist
		PathInfo Igor
		NewPath/Z materials, S_path+"User Procedures:materials"
	endif
	String list = keyStrFromFile(fname,"CrystalStructure","materials")
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

	PathInfo/S materials
	GetFileFolderInfo/P=$path/Q/Z fname
	fullName = S_Path
	if (!V_isFile)
		Open/D/M="file with a lattice keyword list"/P=$path/R/T=".xtl" refNum as fname
		fullName = S_fileName
	endif
	if (strlen(fullName)<1)
		return ""
	else
		Open/P=$path/R/Z=1 refNum as fullName
		fullName = S_fileName
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
	printIt = ParamIsDefault(printIt) ? 0 : printIt

	fname = FindMaterialsFile(fname)				// find full path to fname, and optionally set materials path
//	PathInfo materials
//	if (strlen(S_path)<1)							// make it if it does not exist
//		NewPath/Z materials, ParseFilePath(1,FunctionPath("CrystalsAreHere"),":",1,0)
//	endif
//	PathInfo materials
//	if (strlen(S_path)<1)							// make it if it does not exist
//		PathInfo Igor
//		NewPath/Z materials, S_path+"User Procedures:materials"
//	endif
//	String fileFilters = "XML Files (*.xml):.xml;XTL Files (*.xtl):.xtl;old cri Files (*.cri):.cri;All Files:.*;"
//	Variable f
//	Open/D=2/R/F=fileFilters/M="File with Crystal Information"/P=materials f as fname
//	fname = S_fileName
//	if (strlen(fname)==0)
//		return 1
//	endif

	// get extension, and extension to choose the appropriate read crystal function
	String extension = ParseFilePath(4,fname,":",0,0)
	Variable err=1
	if (stringmatch(extension,"xml"))
		err = readCrystalStructureXML(xtal,fname)
	elseif (stringmatch(extension,"xtl"))
		err = readCrystalStructure_xtl(xtal,fname)
	elseif (stringmatch(extension,"cif"))
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
Static Function/S FindMaterialsFile(fname)					// returns full path to a materials file
	String fname
	String fileFilters = "XML Files (*.xml):.xml;XTL Files (*.xtl):.xtl;CIF Files (*.cif):.cif;old cri Files (*.cri):.cri;All Files:.*;"

	GetFileFolderInfo/Q/Z fname						// full path name passed, that was east
	if (V_isfile)
		return fname
	endif

	String dirString, dirList=""							// a list of possible materials directories
	PathInfo materials
	String materialsPath = SelectString(V_flag,"__empty__",S_path)			// materials path is already set, the 1st choice

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
	if (StringMatch(materialsPath,usersPath))
		dirList = "User files;"
		dirKeyList = ReplaceStringByKey("User files",dirKeyList,materialsPath)
	elseif (StringMatch(materialsPath,docPath))
		dirList = "Documents;"
		dirKeyList = ReplaceStringByKey("Documents",dirKeyList,materialsPath)
	elseif (StringMatch(materialsPath,stdPath))
		dirList = "Standard Distribution;"
		dirKeyList = ReplaceStringByKey("Standard Distribution",dirKeyList,materialsPath)
	elseif (StringMatch(materialsPath,appPath))
		dirList = "Igor application folder;"
		dirKeyList = ReplaceStringByKey("Igor application folder",dirKeyList,materialsPath)
	else
		dirList = SelectString(StringMatch(materialsPath,"__empty__"), materialsPath+";","")
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
				LatticeSym#UpdateMaterialsPath(dirString)
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

	LatticeSym#UpdateMaterialsPath(StringByKey(dirString,dirKeyList))

	Variable f
	Open/D=2/R/F=fileFilters/M="File with Crystal Information"/P=materials f as fname
	fname = S_fileName
	return fname
End
//Static Function/S FindMaterialsFile(fname)					// returns full path to a materials file
//	String fname
////	String fileFilters = "XML Files (*.xml):.xml;XTL Files (*.xtl):.xtl;old cri Files (*.cri):.cri;All Files:.*;"
//	String fileFilters = "XML Files (*.xml):.xml;XTL Files (*.xtl):.xtl;CIF Files (*.cif):.cif;old cri Files (*.cri):.cri;All Files:.*;"
//	String dirString, dirList=""							// a list of possible materials directories
//
//	PathInfo materials
//	dirList += SelectString(V_flag,"",S_path+";")							// materials path is already set, the 1st choice
//
//	dirString = SpecialDirPath("Igor Pro User Files",0,0,0)+"materials:"	// user's local copy, 2nd choice
//	GetFileFolderInfo/Q/Z dirString
//	dirList += SelectString(WhichListItem(dirString,dirList)<0 && V_isFolder,"",dirString+";")
//
//	PathInfo Igor
//	dirString = SelectString(V_flag,"",S_path+"materials:;")				// system local copy, 3rd choice
//	GetFileFolderInfo/Q/Z dirString
//	dirList += SelectString(WhichListItem(dirString,dirList)<0 && V_isFolder,"",dirString+";")
//
//	dirString=ParseFilePath(1,FunctionPath("CrystalsAreHere"),":",1,0)	// the default distribution, 4th choice
//	dirList += SelectString(WhichListItem(dirString,dirList)<0 && strlen(dirString),"",dirString+";")
//	// printf "dirList =  '%s'\r",dirList
//
//	GetFileFolderInfo/Q/Z fname						// full path name passed, that was east
//	if (V_isfile)
//		return fname
//	endif
//
//	fname = ParseFilePath(0,fname,":",1,0)			// check for fname in each of dirList
//	Variable i
//	for (i=0;i<ItemsInList(dirList);i+=1)
//		dirString = StringFromList(i,dirList)
//		GetFileFolderInfo/Q/Z dirString+fname
//		if (V_isFile)									// found fname in dirString
//			UpdateMaterialsPath(dirString)
//			return dirString+fname 
//		endif
//	endfor
//
//	// could not find the materials file, time to start asking user for help
//	dirString = StringFromList(0,dirList)
//	if (itemsInList(dirList)>1)
//		Prompt dirString,"Folder with materials",popup,dirList
//		DoPrompt "Directory",dirString
//		if (V_flag)
//			return ""
//		endif
//	endif
//	UpdateMaterialsPath(dirString)
//
//	Variable f
//	Open/D=2/R/F=fileFilters/M="File with Crystal Information"/P=materials f as fname
//
//	fname = S_fileName
//
//	// printf "dirList =  '%s'\r",dirList
//	// printf "fname =  '%s'\r",fname
//	return fname
//End
//
Static Function UpdateMaterialsPath(dirString)
	String dirString
	GetFileFolderInfo/Q/Z dirString
	if (!V_isFolder)
		return 1
	endif
	PathInfo materials
	if (V_flag && stringmatch(S_path,dirString))		// nothing to change, all is OK
		return 0
	endif
//	NewPath/Z materials, dirString
	NewPath/Z/O materials, dirString
	return 0
End

Static Function readCrystalStructureXML(xtal,fname)
	STRUCT crystalStructure &xtal						// this sruct is filled  by this routine
	String fname

	PathInfo materialsXML
	if (strlen(S_path)<1)								// make it if it does not exist
//		NewPath/Z materialsXML, ParseFilePath(1,FunctionPath("CrystalsAreHere"),":",1,0)+"xml"
		NewPath/Z materialsXML, ParseFilePath(1,FunctionPath("CrystalsAreHere"),":",1,0)
	endif
	PathInfo materialsXML
	if (strlen(S_path)<1)								// make it if it does not exist
		PathInfo Igor
		NewPath/Z materialsXML, S_path+"User Procedures:materials:xml"
	endif

	Variable err = readFileXML(xtal,fname,path="materialsXML")
	if (!err)
		UpdateCrystalStructureDefaults(xtal)
	endif
	return err
End
//
Static Function readFileXML(xtal,fileName,[path])
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

	Variable i0,i1
	i0 = strsearch(buf,"<cif>",0,2)
	if (i0<0)
		return 1
	endif
	i0 += 5												// start just after <cif>
	i1 = strsearch(buf,"</cif>",0,2)
	if (i1<0)
		return 1
	endif
	i1 -= 1												// ends just before </cif>
	String cif = buf[i0,i1]

	i1 = strlen(fullFile)-1						// possibly trim file length to fit, keep last part of fullFile
	i0 = max(0,i1-MAX_FILE_LEN)
	fullFile = fullFile[i0,i1]
	xtal.sourceFile = fullFile

	xtal.hashID = ""
	String str = XMLtagContents("chemical_name_common",cif)
	xtal.desc = str[0,99]
	Variable SG = str2num(XMLtagContents("space_group_IT_number",cif))
	xtal.SpaceGroup = isValidSpaceGroup(SG) ? SG : 0

	String cell = XMLtagContents("cell",cif)				// cell group
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

	// collect the atom sites
	String atomSite, atomLabel, symbol, WyckoffSymbol
	Variable Zatom, fracX,fracY,fracZ,occupy,valence,DebyeT,N=0
	Variable Biso, Uiso, aU11,aU22,aU33, aU12,aU13,aU23, mult
	i0 = strsearch(cif,"<atom_site",i0,2)
	do
		atomSite = XMLtagContents("atom_site",cif[i0,Inf])	// one atom site
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
			if (ForceXYZtoWyckoff(SG,WyckoffSymbol,fracX,fracY,fracZ))
				fracX = NaN													// will cause a break in next if()
			endif
		endif
		if (numtype(fracX+fracY+fracZ))									// give up, cannot determine coordinates
			break
		endif
		mult = 0
		if (strlen(WyckoffSymbol)==0)									// try to set Wyckoff symbol from coordinates
			WyckoffSymbol = FindWyckoffSymbol(SG,fracX,fracY,fracZ,mult)
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
		i0 = strsearch(cif,"<atom_site",i0+10,2)
		N += 1
	while(N<STRUCTURE_ATOMS_MAX && i0>0)
	xtal.N = N											// number of atoms described here
	ForceXtalAtomNamesUnique(xtal)					// forces all of the xtal atom names to be unique

	Variable i, unitsConvert, Nbond=0, Nlen
	String bondKeys, label0,label1,list
	i0 = strsearch(cif,"<bond_chemical",i0,2)
	do
		Wave blen = str2vec(XMLtagContents("bond_chemical",cif[i0,Inf]))	// bond length(s)
		bondKeys = XMLattibutes2KeyList("bond_chemical",cif[i0,Inf])
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
		i0 = strsearch(cif,"<bond_chemical",i0+10,2)
	while(Nbond<(2*STRUCTURE_ATOMS_MAX) && i0>0)
	xtal.Nbonds = Nbond
	ForceLatticeToStructure(xtal)
	return 0
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
			cif += "\t\t<WyckoffSymbol>"+xtal.atom[i].WyckoffSymbol+"</WyckoffSymb>"+NL
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

	PathInfo materials
	if (strlen(S_path)<1)							// make it if it does not exist
		NewPath/Z materials, ParseFilePath(1,FunctionPath("CrystalsAreHere"),":",1,0)
	endif
	PathInfo materials
	if (strlen(S_path)<1)							// make it if it does not exist
		PathInfo Igor
		NewPath/Z materials, S_path+"User Procedures:materials"
	endif
	String fileFilters = "XTL Files (*.xtl):.xtl;old cri Files (*.cri):.cri;All Files:.*;"
	Variable f
	Open/D=2/R/F=fileFilters/M="File with Crystal Information"/P=materials f as xtlName
	xtlName = S_fileName
	if (strlen(xtlName)==0)
		return 1
	endif

	// check extension
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

Static Function readCrystalStructureCIF(xtal,fname)
	STRUCT crystalStructure &xtal						// this sruct is filled  by this routine
	String fname

	PathInfo materialsXML
	if (strlen(S_path)<1)								// make it if it does not exist
		NewPath/Z materialsXML, ParseFilePath(1,FunctionPath("CrystalsAreHere"),":",1,0)
	endif
	PathInfo materialsXML
	if (strlen(S_path)<1)								// make it if it does not exist
		PathInfo Igor
		NewPath/Z materialsXML, S_path+"User Procedures:materials:xml"
	endif

	Variable err = readFileCIF(xtal,fname,path="materialsXML")
	if (!err)
		UpdateCrystalStructureDefaults(xtal)
	endif
	return err
End
//
Static Function readFileCIF(xtal,fileName,[path])
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
//	i = strsearch(buf,"\n",i+1)+1
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
	endif
	return err
End

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

	Variable SG=CIF_readNumber("_symmetry_Int_Tables_number",buf)
	if (numtype(SG))
		String SGstr=CIF_readString("_symmetry_space_group_name_H-M",buf)
		SG = str2num(SymString2SG(SGstr,1))
	endif
	xtal.SpaceGroup = SG > 0 ? SG : 1
	xtal.Temperature = CIF_readNumber("_cell_measurement_temperature",buf)-273.15 	// temperature (K) for cell parameters

	xtal.N = 0
	// find the atoms
	Variable i, cc
	String list=""					// find loop_ with the fractional atom positions
	for (i=strsearch(buf,"\nloop_",0); i>=0; i=strsearch(buf,"\nloop_",i+5))
		list = CIF_loop_Labels(buf[i,Inf])
		if (WhichListItem("_atom_site_fract_x",list)>=0)
			break
		endif
	endfor

	if (strlen(list))				// found atom positions, interpret the list
		buf = buf[i,Inf]
		String line, last=StringFromList(ItemsInList(list)-1,list)
		i = strsearch(buf,"\n"+last+"\n",0)+strlen(last)+2		// i is now start of table values
		buf = buf[i,Inf]
		for (i=strsearch(buf,"\n",0); i>=0; i=strsearch(buf,"\n",i+1))	// find end of list
			cc = char2num(buf[i+1])				// values end with a line having an "_", an empty line, or line starting with "#"
			if (cc<=10 || cc==35 || cc==95)
				break
			endif
		endfor
		if (i>0)
			buf = buf[0,i]
		endif

		buf = buf[0,i]								// buf now contains ONLY the list values
		buf = ReplaceString(" ",buf,"\t")
		do
			buf = ReplaceString("\t\t",buf,"\t")	// change all multiple tabs to single tabs
		while (strsearch(buf,"\t\t",0)>=0)
		buf = ReplaceString("\n\t",buf,"\n")
		buf = ReplaceString("\t\n",buf,"\n")
		if (char2num(buf)==9)
			buf = ReplaceString("\t",buf,"",0,1)		// remove the leading tab
		endif
		buf = ReplaceString("\t",buf,";")			// turns white space separated table into a set of lists

		String symb
		Variable Z, occ, Biso, Uiso, Uij

		Variable N, i0,i1
		for (N=0,i0=0; N<STRUCTURE_ATOMS_MAX && i0<strlen(buf); N+=1)
			i1 = strsearch(buf,"\n",i0)
			line = buf[i0,i1-1]								// un-terminated line
			xtal.atom[N].DebyeT = NaN
			xtal.atom[N].Uiso = NaN	;	xtal.atom[N].Biso = NaN
			xtal.atom[N].U11 = NaN		;	xtal.atom[N].U22 = NaN	;		xtal.atom[N].U33 = NaN
			xtal.atom[N].U12 = NaN		;	xtal.atom[N].U13 = NaN	;		xtal.atom[N].U23 = NaN

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
			i0 = i1+1										// start of next line
		endfor
		xtal.N = N
	endif
	xtal.Nbonds	 = 0
	xtal.Unconventional00 = NaN
	xtal.sourceFile = ""
	ForceLatticeToStructure(xtal)
//	xtal.hashID = xtalHashID(xtal)
	return 0
End
//
Static Function/T CIF_loop_Labels(buf)
	String buf

	Variable i1, i0 = strsearch(buf,"\nloop_\n",0)
	i0 += 7
	String name, list=""
	for(;char2num(buf[i0])==95;)
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

	Variable len=strlen(key)
	Variable i0 = strsearch(buf,"\n"+key,0)

	if (i0<0)
		return ""
	endif

	i0 = strsearch(buf,"'",i0+1)+1
	Variable i1 = strsearch(buf,"'",i0+1)-1
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
	STRUCT crystalStructure &xtal				// this sruct is filled  by this routine
	Variable h,k,l
	Variable keV
	Variable T_K										// Temperature (K), used to calculate Debye-Waller factor
	keV = ParamIsDefault(keV) ? NaN : keV
	T_K = ParamIsDefault(T_K) ? NaN : max(0,T_K)
	T_K = T_K>0 ? T_K : NaN

	Variable SpaceGroup=xtal.SpaceGroup
	if (!isValidSpaceGroup(SpaceGroup) || numtype(h+k+l))
		return cmplx(NaN,NaN)	// bad inputs
	endif
	Variable/C zero=cmplx(0,0)
	Variable usingHexAxes = (abs(90-xtal.alpha)+abs(90-xtal.beta)+abs(120-xtal.gam))<1e-6
	Variable system = latticeSystem(SpaceGroup)

	if (!(mod(h,1) || mod(k,1) || mod(l,1)))	// non-integral, always allowed
		String sym = getHMsym(SpaceGroup)			// get symmetry symbol
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
	keV = keV>0 ? keV : NumVarOrDefault("root:Packages:Lattices:keV",10)

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

	for (m=0;m<xtal.N;m+=1)							// loop over the defined atoms
		name="root:Packages:Lattices:atom"+num2istr(m)
		valence = numtype(xtal.atom[m].valence) ? 0 : xtal.atom[m].valence
		Wave ww = $name
		FUNCREF Get_f_proto fa= $"Get_f"			// if Get_f() does not exist, then Get_f_proto() will be used
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
			if (itemp==1)
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

	Variable Fr, Fi, arg
	if (system==HEXAGONAL || (usingHexAxes && system==TRIGONAL))
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
	valence = ParamIsDefault(valence) ? 0 : valence

	Variable iz=ZfromLabel(AtomType), Bval
	if (iz<1 || numtype(iz))
		return NaN
	elseif (iz<=10)
		Make/FREE BwTemp={58.3331,10.9071,4.33979,42.9165,23.0888,12.7188,0.02064,13.8964,11.2651, 9}	// use 9 for Ne
		Bval = BwTemp[iz-1]
	elseif (iz<=18)
		Make/FREE coef={6.8736,0.011759,-0.025672}
		Bval = poly(coef,iz-1)
	elseif (iz<=50)
		Make/FREE coef={51.647,-3.5557,0.085621,-0.00070461}
		Bval = poly(coef,iz-1)
	elseif (iz<=58)
		Make/FREE BwTemp={5.24328,4.74225,4.27091,0.26422,0.23092,0.15152,0.1104,0.12335}
		Bval = BwTemp[iz-51]
	elseif (iz<=70)
		Make/FREE coef={30.119,-0.85371,0.006597}
		Bval = poly(coef,iz-1)
	else
		Make/FREE coef={-57.258,2.1161,-0.025226,9.8321e-05}
		Bval = poly(coef,iz-1)
	endif

	Variable Svector = QAngstrom/(4*PI)/10	// must be in (1/Angstrom)
	Variable f0 = iz * exp(-Svector*Svector*(Bval)) - valence
	f0 = numtype(f0) ? 1 : max(f0,0)			// always a valid number > 0
	return cmplx(f0,0)
End
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
//
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
		name="root:Packages:Lattices:atom"+num2istr(m)
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

	Variable SpaceGroup=xtal.SpaceGroup
	SetSymOpsForSpaceGroup(SpaceGroup)			// ensure existance of symmetry op mats and vecs
	Wave mats = $("root:Packages:Lattices:SymOps:equivXYZM"+num2istr(SpaceGroup))
	Wave bvecs = $("root:Packages:Lattices:SymOps:equivXYZB"+num2istr(SpaceGroup))
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
		rr += abs(floor(rr))						// translate to a unit cell with only positive values
		rr = mod(rr,1)									//   and restrict values to values to [0,1), the first unit cell

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

	if (!isValidSpaceGroup(xtal.SpaceGroup) || numtype(h+k+l))
		return 0												// bad inputs (not allowed)
	endif
	Variable usingHexAxes = (abs(90-xtal.alpha)+abs(90-xtal.beta)+abs(120-xtal.gam))<1e-6
	Variable system = LatticeSym#latticeSystem(xtal.SpaceGroup)

	if (!(mod(h,1) || mod(k,1) || mod(l,1)))	// non-integral, always allowed
		String sym = getHMsym(xtal.SpaceGroup)	// get symmetry symbol
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

	Variable fatomMag, m, Fr, Fi
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
		else
			Fc += cmplx(fatomMag,0)					// no atom position, just make it in phase
		endif
	endfor

	if (system==HEXAGONAL || (usingHexAxes && system==TRIGONAL))
		Variable arg
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

	return (magsqr(Fc)/(xtal.N)^2 > 0.0001)			// allowed means more than 0.01 electron/atom
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
		Wave SymmetryOp = $LatticeSym#MakeSymmetryOps(xtal)	// wave with all the symmetry operation
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


ThreadSafe Function/S MakeSymmetryOps(xtal)				// make a wave with the symmetry operation
	STRUCT crystalStructure &xtal

	String SG = num2istr(xtal.SpaceGroup)
	Wave equivXYZM = $("root:Packages:Lattices:SymOps:equivXYZM"+SG)
	if (!WaveExists(equivXYZM))
		return ""
	endif
	Variable Nequiv=DimSize(equivXYZM,0)				// number of all symmetry operations for this space group
	String/G root:Packages:Lattices:SymOps:SymmetryOpsPath="root:Packages:Lattices:SymOps:SymmetryOps"+SG

	String wName = "root:Packages:Lattices:SymOps:SymmetryOps"+SG
	Make/N=(2*Nequiv,3,3)/D/O $wName=NaN
	Wave ops = $wName

	Make/N=(3,3)/D/FREE direct,mat
	direct = { {xtal.a0, xtal.a1, xtal.a2}, {xtal.b0, xtal.b1, xtal.b2}, {xtal.c0, xtal.c1, xtal.c2} }
	MatrixOp/FREE/O directI = Inv(direct)

	Variable i, N=0, Nproper=0
	for (i=0;i<Nequiv;i+=1)									// loop through all Nequiv, rejecting duplicates, and only accepting proper rotations
		mat = equivXYZM[i][p][q]
		MatrixOp/O/FREE mat = direct x mat x directI		// convert to cartesian, similarity transform
		if (isMatInMats(mat,ops) ||MatrixDet(mat)<0)	// skip duplicates and improper rotations
			continue
		endif
		ops[N][][] = mat[q][r]								// save this mat
		Nproper += MatrixDet(mat)>0 ? 1 : 0
		N += 1
	endfor

	// go through list again, this time only taking unique IMproper rotations
	for (i=0;i<Nequiv;i+=1)									// again, loop through all Nequiv, rejecting duplicates
		mat = equivXYZM[i][p][q]
		MatrixOp/O/FREE mat = direct x mat x directI		// convert to cartesian, similarity transform
		if (isMatInMats(mat,ops))							// skip duplicates
			continue
		endif
		ops[N][][] = mat[q][r]								// save this mat
		N += 1
	endfor
	Redimension/N=(N,-1,-1) ops

	ops = abs(ops[p][q][r])<1e-13 ? 0 : ops[p][q][r]		// remove the almost zeros
	ops = abs(1-ops[p][q][r])<1e-13 ? 1 : ops[p][q][r]	//  and make almost 1's equal to 1
	ops = abs(1+ops[p][q][r])<1e-13 ? -1 : ops[p][q][r]

	String wnote="waveClass=SymmetryOperations;"
	wnote = ReplaceNumberByKey("Nproper",wnote,Nproper,"=")
	wnote = ReplaceNumberByKey("SpaceGroup",wnote,xtal.SpaceGroup,"=")
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
		String SymmetryOpsPath=StrVarOrDefault("root:Packages:Lattices:SymOps:SymmetryOpsPath","")
		Wave SymOps = $StrVarOrDefault("root:Packages:Lattices:SymOps:SymmetryOpsPath","")
		Variable SGw = WaveExists(SymOps) ? NumberByKey("SpaceGroup",note(SymOps),"=") : 0
		STRUCT crystalStructure xtal
		if (FillCrystalStructDefault(xtal) && !WaveExists(SymOps))
			return ""
		endif

		Variable useWave=0, useXtal=0
		if (isValidSpaceGroup(SGw) && isValidSpaceGroup(xtal.SpaceGroup) && xtal.SpaceGroup == SGw)
			useWave = 1
		elseif (isValidSpaceGroup(SGw) && !isValidSpaceGroup(xtal.SpaceGroup))		// use wave
			useWave = 1
		elseif (!isValidSpaceGroup(SGw) && isValidSpaceGroup(xtal.SpaceGroup))		// use xtal
			useXtal = 1
		elseif (isValidSpaceGroup(SGw) && isValidSpaceGroup(xtal.SpaceGroup))		// ask
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
	Variable SpaceGroup = NumberByKey("SpaceGroup",note(SymOps),"=")
	if (!(Nproper>0))
		return ""
	endif
	if (printIt)
		String system = StringFromList(latticeSystem(SpaceGroup),LatticeSystemNames)
		printf "For Space Group %g  (%s) %s,  %g proper rotations\r",SpaceGroup,system,getFullHMSym(SpaceGroup),Nproper
	endif

	Make/N=3/D/FREE axis, v3
	Make/N=(3,3)/D/FREE sym
	String out = ReplaceNumberByKey("Nproper","",Nproper,"=")
	out = ReplaceNumberByKey("SpaceGroup",out,SpaceGroup,"=")
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

ThreadSafe Function/S getHMboth(SpaceGroup)	// returns short and (full) Hermann-Mauguin symbol
	Variable SpaceGroup						//Space Group number, from International Tables

	String short = getHMsym(SpaceGroup)
	String full = getFullHMSym(SpaceGroup)
	if (StringMatch(short,full))
		return short
	else
		return short + "  ("+full+")"
	endif
End


ThreadSafe Function/S getHMsym(SpaceGroup)	// returns short Hermann-Mauguin symbol
	Variable SpaceGroup						//Space Group number, from International Tables
	if (!isValidSpaceGroup(SpaceGroup))
		return ""								// invalid SpaceGroup number
	endif

	// there are 230 items in this list
	String symms="P1;P-1;P2:b;P21:b;C2:b1;Pm:b;Pc:b1;"
	symms += "Cm:b1;Cc:b1;P2/m:b;P21/m:b;C2/m:b1;P2/c:b1;"
	symms += "P21/c:b1;C2/c:b1;P222;P2221;P21212;P212121;"
	symms += "C2221;C222;F222;I222;I212121;Pmm2;Pmc21;"
	symms += "Pcc2;Pma2;Pca21;Pnc2;Pmn21;Pba2;Pna21;"
	symms += "Pnn2;Cmm2;Cmc21;Ccc2;Amm2;Abm2;Ama2;Aba2;"
	symms += "Fmm2;Fdd2;Imm2;Iba2;Ima2;Pmmm;Pnnn:1;Pccm;"
	symms += "Pban:1;Pmma;Pnna;Pmna;Pcca;Pbam;Pccn;Pbcm;"
	symms += "Pnnm;Pmmn:1;Pbcn;Pbca;Pnma;Cmcm;Cmca;Cmmm;"
	symms += "Cccm;Cmma;Ccca:1;Fmmm;Fddd:1;Immm;Ibam;Ibca;"
	symms += "Imma;P4;P41;P42;P43;I4;I41;P-4;I-4;"
	symms += "P4/m;P42/m;P4/n:1;P42/n:1;I4/m;I41/a:1;P422;"
	symms += "P4212;P4122;P41212;P4222;P42212;P4322;P43212;"
	symms += "I422;I4122;P4mm;P4bm;P42cm;P42nm;P4cc;"
	symms += "P4nc;P42mc;P42bc;I4mm;I4cm;I41md;I41cd;"
	symms += "P-42m;P-42c;P-421m;P-421c;P-4m2;P-4c2;P-4b2;"
	symms += "P-4n2;I-4m2;I-4c2;I-42m;I-42d;P4/mmm;P4/mcc;"
	symms += "P4/nbm:1;P4/nnc:1;P4/mbm;P4/mnc;P4/nmm:1;P4/ncc:1;"
	symms += "P42/mmc;P42/mcm;P42/nbc:1;P42/nnm:1;P42/mbc;P42/mnm;"
	symms += "P42/nmc:1;P42/ncm:1;I4/mmm;I4/mcm;I41/amd:1;"
	symms += "I41/acd:1;P3;P31;P32;R3:H;P-3;R-3:H;P312;"
	symms += "P321;P3112;P3121;P3212;P3221;R32:H;P3m1;P31m;"
	symms += "P3c1;P31c;R3m:H;R3c:H;P-31m;P-31c;P-3m1;P-3c1;"
	symms += "R-3m:H;R-3c:H;P6;P61;P65;P62;P64;P63;P-6;"
	symms += "P6/m;P63/m;P622;P6122;P6522;P6222;P6422;P6322;"
	symms += "P6mm;P6cc;P63cm;P63mc;P-6m2;P-6c2;P-62m;P-62c;"
	symms += "P6/mmm;P6/mcc;P63/mcm;P63/mmc;P23;F23;I23;P213;"
	symms += "I213;Pm-3;Pn-3:1;Fm-3;Fd-3:1;Im-3;Pa-3;Ia-3;"
	symms += "P432;P4232;F432;F4132;I432;P4332;P4132;I4132;"
	symms += "P-43m;F-43m;I-43m;P-43n;F-43c;I-43d;Pm-3m;"
	symms += "Pn-3n:1;Pm-3n;Pn-3m:1;Fm-3m;Fm-3c;Fd-3m:1;Fd-3c:1;"
	symms += "Im-3m;Ia-3d;"
	return StringFromList(SpaceGroup-1,symms)	// set the symmetry symbol
End


// Full H-M symbols only differ from the regular ones (in getHMsym) for Space Groups: 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
ThreadSafe Function/S getFullHMSym(SpaceGroup)	// returns full Hermann-Mauguin symbol
	Variable SpaceGroup					//Space Group number, from International Tables
	if (!isValidSpaceGroup(SpaceGroup))
		return ""							// invalid SpaceGroup number
	endif

	// there are 230 items in this list
	String fullSymms="P1;P-1;P121:b;P1211:b;C121:b1;P1m1:b;"
	fullSymms += "P1c1:b1;C1m1:b1;C1c1:b1;P12/m1:b;P121/m1:b;"
	fullSymms += "C12/m1:b1;P12/c1:b1;P121/c1:b1;C12/c1:b1;P222;"
	fullSymms += "P2221;P21212;P212121;C2221;C222;F222;I222;I212121;"
	fullSymms += "Pmm2;Pmc21;Pcc2;Pma2;Pca21;Pnc2;Pmn21;Pba2;Pna21;"
	fullSymms += "Pnn2;Cmm2;Cmc21;Ccc2;Amm2;Abm2;Ama2;Aba2;Fmm2;"
	fullSymms += "Fdd2;Imm2;Iba2;Ima2;Pmmm;Pnnn:1;Pccm;Pban:1;Pmma;"
	fullSymms += "Pnna;Pmna;Pcca;Pbam;Pccn;Pbcm;Pnnm;Pmmn:1;Pbcn;"
	fullSymms += "Pbca;Pnma;Cmcm;Cmca;Cmmm;Cccm;Cmma;Ccca:1;Fmmm;"
	fullSymms += "Fddd:1;Immm;Ibam;Ibca;Imma;P4;P41;P42;P43;I4;I41"
	fullSymms += ";P-4;I-4;P4/m;P42/m;P4/n:1;P42/n:1;I4/m;I41/a:1;"
	fullSymms += "P422;P4212;P4122;P41212;P4222;P42212;P4322;P43212;"
	fullSymms += "I422;I4122;P4mm;P4bm;P42cm;P42nm;P4cc;P4nc;P42mc;"
	fullSymms += "P42bc;I4mm;I4cm;I41md;I41cd;P-42m;P-42c;P-421m;"
	fullSymms += "P-421c;P-4m2;P-4c2;P-4b2;P-4n2;I-4m2;I-4c2;I-42m;"
	fullSymms += "I-42d;P4/mmm;P4/mcc;P4/nbm:1;P4/nnc:1;P4/mbm;P4/mnc;"
	fullSymms += "P4/nmm:1;P4/ncc:1;P42/mmc;P42/mcm;P42/nbc:1;P42/nnm:1;"
	fullSymms += "P42/mbc;P42/mnm;P42/nmc:1;P42/ncm:1;I4/mmm;I4/mcm;"
	fullSymms += "I41/amd:1;I41/acd:1;P3;P31;P32;R3:H;P-3;R-3:H;P312;"
	fullSymms += "P321;P3112;P3121;P3212;P3221;R32:H;P3m1;P31m;P3c1;"
	fullSymms += "P31c;R3m:H;R3c:H;P-31m;P-31c;P-3m1;P-3c1;R-3m:H;"
	fullSymms += "R-3c:H;P6;P61;P65;P62;P64;P63;P-6;P6/m;P63/m;P622;"
	fullSymms += "P6122;P6522;P6222;P6422;P6322;P6mm;P6cc;P63cm;P63mc;"
	fullSymms += "P-6m2;P-6c2;P-62m;P-62c;P6/mmm;P6/mcc;P63/mcm;"
	fullSymms += "P63/mmc;P23;F23;I23;P213;I213;Pm-3;Pn-3:1;Fm-3;"
	fullSymms += "Fd-3:1;Im-3;Pa-3;Ia-3;P432;P4232;F432;F4132;I432;"
	fullSymms += "P4332;P4132;I4132;P-43m;F-43m;I-43m;P-43n;F-43c;"
	fullSymms += "I-43d;Pm-3m;Pn-3n:1;Pm-3n;Pn-3m:1;Fm-3m;Fm-3c;"
	fullSymms += "Fd-3m:1;Fd-3c:1;Im-3m;Ia-3d"
	return StringFromList(SpaceGroup-1,fullSymms)	// set the symmetry symbol
End


ThreadSafe Function/S getHallSymbol(SpaceGroup)
	Variable SpaceGroup					//Space Group number, from International Tables
	if (!isValidSpaceGroup(SpaceGroup))
		return ""							// invalid SpaceGroup number
	endif

	// there are 230 items in this list
	String Hall="P 1;-P 1;P 2y:b;P 2yb:b;C 2y:b1;P -2y:b;P -2yc:b1;"
	Hall += "C -2y:b1;C -2yc:b1;-P 2y:b;-P 2yb:b;-C 2y:b1;-P 2yc:b1;"
	Hall += "-P 2ybc:b1;-C 2yc:b1;P 2 2;P 2c 2;P 2 2ab;P 2ac 2ab;"
	Hall += "C 2c 2;C 2 2;F 2 2;I 2 2;I 2b 2c;P 2 -2;P 2c -2;"
	Hall += "P 2 -2c;P 2 -2a;P 2c -2ac;P 2 -2bc;P 2ac -2;P 2 -2ab;"
	Hall += "P 2c -2n;P 2 -2n;C 2 -2;C 2c -2;C 2 -2c;A 2 -2;A 2 -2c;"
	Hall += "A 2 -2a;A 2 -2ac;F 2 -2;F 2 -2d;I 2 -2;I 2 -2c;I 2 -2a;"
	Hall += "-P 2 2;P 2 2 -1n:1;-P 2 2c;P 2 2 -1ab:1;-P 2a 2a;"
	Hall += "-P 2a 2bc;-P 2ac 2;-P 2a 2ac;-P 2 2ab;-P 2ab 2ac;"
	Hall += "-P 2c 2b;-P 2 2n;P 2 2ab -1ab:1;-P 2n 2ab;-P 2ac 2ab;"
	Hall += "-P 2ac 2n;-C 2c 2;-C 2bc 2;-C 2 2;-C 2 2c;-C 2b 2;"
	Hall += "C 2 2 -1bc:1;-F 2 2;F 2 2 -1d:1;-I 2 2;-I 2 2c;-I 2b 2c;"
	Hall += "-I 2b 2;P 4;P 4w;P 4c;P 4cw;I 4;I 4bw;P -4;I -4;-P 4;"
	Hall += "-P 4c;P 4ab -1ab:1;P 4n -1n:1;-I 4;I 4bw -1bw:1;P 4 2;"
	Hall += "P 4ab 2ab;P 4w 2c;P 4abw 2nw;P 4c 2;P 4n 2n;P 4cw 2c;"
	Hall += "P 4nw 2abw;I 4 2;I 4bw 2bw;P 4 -2;P 4 -2ab;P 4c -2c;"
	Hall += "P 4n -2n;P 4 -2c;P 4 -2n;P 4c -2;P 4c -2ab;I 4 -2;"
	Hall += "I 4 -2c;I 4bw -2;I 4bw -2c;P -4 2;P -4 2c;P -4 2ab;"
	Hall += "P -4 2n;P -4 -2;P -4 -2c;P -4 -2ab;P -4 -2n;I -4 -2;"
	Hall += "I -4 -2c;I -4 2;I -4 2bw;-P 4 2;-P 4 2c;P 4 2 -1ab:1;"
	Hall += "P 4 2 -1n:1;-P 4 2ab;-P 4 2n;P 4ab 2ab -1ab:1;"
	Hall += "P 4ab 2n -1ab:1;-P 4c 2;-P 4c 2c;P 4n 2c -1n:1;"
	Hall += "P 4n 2 -1n:1;-P 4c 2ab;-P 4n 2n;P 4n 2n -1n:1;"
	Hall += "P 4n 2ab -1n:1;-I 4 2;-I 4 2c;I 4bw 2bw -1bw:1;"
	Hall += "I 4bw 2aw -1bw:1;P 3;P 31;P 32;R 3:H;-P 3;-R 3:H;"
	Hall += "P 3 2;P 3 2'';P 31 2c (0 0 1);P 31 2'';P 32 2c (0 0 -1);"
	Hall += "P 32 2'';R 3 2'':H;P 3 -2'';P 3 -2;P 3 -2''c;P 3 -2c;R 3 -2'':H;"
	Hall += "R 3 -2''c:H;-P 3 2;-P 3 2c;-P 3 2'';-P 3 2''c;-R 3 2'':H;"
	Hall += "-R 3 2''c:H;P 6;P 61;P 65;P 62;P 64;P 6c;P -6;-P 6;-P 6c;"
	Hall += "P 6 2;P 61 2 (0 0 -1);P 65 2 (0 0 1);P 62 2c (0 0 1);"
	Hall += "P 64 2c (0 0 -1);P 6c 2c;P 6 -2;P 6 -2c;P 6c -2;P 6c -2c;"
	Hall += "P -6 2;P -6c 2;P -6 -2;P -6c -2c;-P 6 2;-P 6 2c;-P 6c 2;"
	Hall += "-P 6c 2c;P 2 2 3;F 2 2 3;I 2 2 3;P 2ac 2ab 3;I 2b 2c 3;"
	Hall += "-P 2 2 3;P 2 2 3 -1n:1;-F 2 2 3;F 2 2 3 -1d:1;-I 2 2 3;"
	Hall += "-P 2ac 2ab 3;-I 2b 2c 3;P 4 2 3;P 4n 2 3;F 4 2 3;F 4d 2 3;"
	Hall += "I 4 2 3;P 4acd 2ab 3;P 4bd 2ab 3;I 4bd 2c 3;P -4 2 3;"
	Hall += "F -4 2 3;I -4 2 3;P -4n 2 3;F -4c 2 3;I -4bd 2c 3;-P 4 2 3;"
	Hall += "P 4 2 3 -1n:1;-P 4n 2 3;P 4n 2 3 -1n:1;-F 4 2 3;-F 4c 2 3;"
	Hall += "F 4d 2 3 -1d:1;F 4d 2 3 -1cd:1;-I 4 2 3;-I 4bd 2c 3;"
	return StringFromList(SpaceGroup-1,Hall)	// set the Hall symbol
End


ThreadSafe Static Function latticeSystem(SpaceGroup)
	Variable SpaceGroup				//Space Group number, from International Tables
	if (SpaceGroup>230)
		return -1					  	 // invalid
	elseif (SpaceGroup>=195)
		return CUBIC
	elseif (SpaceGroup>=168)
		return HEXAGONAL
	elseif (SpaceGroup>=143)
		return TRIGONAL				// Trigonal, (using the hexagonal cell axes)
	elseif (SpaceGroup>=75)
		return TETRAGONAL
	elseif (SpaceGroup>=16)
		return ORTHORHOMBIC
	elseif (SpaceGroup>=3)
		return MONOCLINIC
	elseif  (SpaceGroup>0)
		return TRICLINIC
	endif
	return -1							// invalid
End


Function/S symmtry2SG(strIN,[type,printIt])	// find the Space Group number from the symmetry string
	String strIN
	Variable type						// 0=Check All, 1=Hermann-Mauguin, 2=Full Hermann-Mauguin, 3=Hall, 4=Lattice System, 5=SpaceGroup Number, 
	Variable printIt
	type = ParamIsDefault(type) ? 0 : round(type)
	type = type<0 || type>5 ? NaN : type
	printIt = ParamIsDefault(printIt) ? NaN : printIt
	printIt = numtype(printIt) ?  strlen(GetRTStackInfo(2))<1 : !(!printIt)

	if (strlen(strIN)<1 || numtype(type))
		Prompt strIN, "Symmetry Symbol or Space Group Number, (e.g. \"Pmm*\"), wild cards allowed"
		Prompt type,"Symbol Type",popup,"Check All Symbol Types;Hermann-Mauguin;Full Hermann-Mauguin;Hall;Lattice System;SpaceGroup Number"
		type += 1
		DoPrompt "Symmetry Symbol",strIN,type
		if (V_flag)
			return ""
		endif
		type -= 1
		printIt = 1
	endif
	if (strlen(strIN)<1)
		return ""
	endif
	type = round(type)

	String list="", name=StringFromList(type,"any type of;Hermann-Mauguin;FULL Hermann-Mauguin;Hall;Lattice System;")
	Variable SG
	if (type==5 || type==0)
		SG = round(str2num(strIN))
		SG = (SG>=1 && SG<=230) ? SG : NaN
		if (numtype(SG)==0)
			sprintf list,"%d;",SG
		endif
	endif
	if (type==1 || type==0)
		list += SymString2SG(strIn,1)
	endif
	if (type==2 || type==0)
		list += SymString2SG(strIn,2)
	endif
	if (type==3 || type==0)
		list += SymString2SG(strIn,3)
	endif
	if (type==4 || type==0)
		list += SymString2SG(strIn,4)
	endif

	Variable i,Nlist=ItemsInList(list), last, now
	String str = SortList(list,";",2)
	list = ""
	for (last=-Inf,i=0; i<Nlist; i+=1)		// remove duplicates
		now = str2num(StringFromList(i,str))
		now = numtype(now) ? NaN : now
		now = (now>=1 && now<=230) ? now : NaN
		if (now>last)
			list += num2istr(now)+";"
			last = now
		endif
	endfor
	Nlist = ItemsInList(list)

	// finally print out information about each SG in list
	if (printIt)
		if (Nlist<1)
			printf "No matches of \"%s\" to a %s symbol\r",strIN,name
			return ""
		elseif (Nlist>1)
			printf "There are %g possible matches of  \"%s\"  to %s symbol\r",Nlist,strIN,name
		endif
		printf "SG\t\t\t\tSystem\t\t\t\tH-M\t\t\tHall\r"
		String tab,fullHM,HM, system, systemNames="Triclinic\t;Monoclinic\t;Orthorhombic;Tetragonal\t;Trigonal\t;Hexagonal\t;Cubic\t\t"
		for (i=0; i<Nlist; i+=1)
			SG = str2num(StringFromList(i,list))
			if (isValidSpaceGroup(SG))
				fullHM = getFullHMSym(SG)			// usually fullHM is the same as HM
				HM = getHMSym(SG)
				fullHM = SelectString(StringMatch(fullHM,HM),"\t\tfull H-M = ["+fullHM+"]","")
				tab = SelectString(strlen(getFullHMSym(SG))>5,"\t","")
				system = StringFromList(latticeSystem(SG),systemNames)
				printf "%d   \t-->\t\t%s\t\t%s\t\t%s%s%s\r", SG,system,HM,tab,getHallSymbol(SG),fullHM
			endif
		endfor
	endif
	return list
End


Static Function/S SymString2SG(symIN,type)	// finds space group of a Hermann-Mauguin or Hall symbol, wild cards allowed
	String symIN						// requested symbol, if empty, then a dialog will come up
	Variable type						// 1=Hermann-Mauguin, 2=Full Hermann-Mauguin, 3=Hall, 4=Lattice System

	String find = ReplaceString(" ",symIN,"")	// do not include spaces in search
	String list=""

	type = round(type)
	if (type==1)
		FUNCREF getHMsym symbolFunc = getHMsym
	elseif (type==2)
		FUNCREF getHMsym symbolFunc = getFullHMSym
	elseif (type==3)
		FUNCREF getHMsym symbolFunc = getHallSymbol
	elseif (type==4)
		if (StringMatch("Triclinic",find))
			list += expandRange("1-2",";")+";"
		endif
		if (StringMatch("Monoclinic",find))
			list += expandRange("3-15",";")+";"
		endif
		if (StringMatch("Orthorhombic",find))
			list += expandRange("16-74",";")+";"
		endif
		if (StringMatch("Tetragonal",find))
			list += expandRange("75-142",";")+";"
		endif
		if (StringMatch("Trigonal",find))
			list += expandRange("143-167",";")+";"
		elseif (StringMatch("Rhombohedral",find))
			list += expandRange("143-167",";")+";"
		endif
		if (StringMatch("Hexagonal",find))
			list += expandRange("168-193",";")+";"
		endif
		if (StringMatch("Cubic",find))
			list += expandRange("195-230",";")+";"
		endif
	else
		return ""
	endif

	String sym
	Variable SG
	for (SG=1;SG<=230;SG+=1)					// check all 230 space goups
		sym = symbolFunc(SG)
		if (StringMatch(ReplaceString(" ",sym,"")	,find))// ignore spaces
			list += num2istr(SG)+";"				// found a match, save it
		endif
	endfor
	return list
End


//	DEPRECATED	DEPRECATED	DEPRECATED	DEPRECATED
// This function is DEPRECATED, it is just left here for old stuff, use getHMsym() instead
ThreadSafe Function/S getSymString(SpaceGroup)	// returns short Hermann-Mauguin symbol
	Variable SpaceGroup							//Space Group number, from International Tables
	return getHMsym(SpaceGroup)
End


ThreadSafe Function/C Hex2Rhom(aH,cH)			// convert lattice constants
	Variable aH,cH					// Hexagonal lattice constants
	Variable aR,alpha				// Rhombohedral lattice constants
	aR = (1/3) * sqrt(3*aH^2 + cH^2)
	alpha = asin( 3/2 / sqrt(3+(cH/aH)^2) ) * 2 * 180/PI
	return cmplx(aR,alpha)
End


ThreadSafe Function/C Rhom2Hex(aR,alpha)		// convert lattice constants
	Variable aR, alpha				// Rhombohedral lattice constants
	Variable aH,cH					// Hexagonal lattice constants
	Variable  ca2					// (cH/aH)^2
	ca2 = (3/2 / sin(alpha*PI/180/2))^2 - 3
	aH  = sqrt( (3*aR)^2/(ca2 + 3) )
	cH = aH*sqrt(ca2)
	return cmplx(aH,cH)
End



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

ThreadSafe Function isRhombohedral(SpaceGroup)
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

	String sym = getHMsym(xtal.SpaceGroup)		// symmetry string
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
	STRUCT crystalStructure &xtal					// this sruct is set in this routine

	Variable SpaceGroup=xtal.SpaceGroup			// local value for convienence
	if (!isValidSpaceGroup(SpaceGroup))			// invalid SpaceGroup, it must be in range [1,230]
		return 0
	elseif (numtype(xtal.a + xtal.b + xtal.c + xtal.alpha + xtal.beta + xtal.gam))
		return 0								// no Inf or NaN
	elseif (xtal.a <= 0 || xtal.b <= 0 || xtal.c <= 0 || xtal.alpha <= 0 || xtal.beta <= 0 || xtal.gam <= 0)
		return 0								// must be > 0
	endif

	Variable lenTol=max(max(xtal.a,xtal.b),xtal.c) * 1e-4, angleTol=1e-4*180/PI
	Variable num90s = (abs(xtal.alpha - 90)<angleTol) + (abs(xtal.beta - 90)<angleTol) + (abs(xtal.gam - 90)<angleTol)

	if (SpaceGroup>=195)					// Cubic [195,230]
		if (num90s==3 && (abs(xtal.a - xtal.b)<lenTol) && (abs(xtal.a - xtal.c)<lenTol) )
			return 1							// require alpha=beta=gamma=90, and a=b=c
		endif
	elseif(SpaceGroup>=168)				// Hexagonal [168,194]
		if ( num90s==2 && abs(xtal.gam - 120) > angleTol && abs(xtal.a - xtal.b)<lenTol )
			return 1							// require a==b, alpha=beta=90, gamma=190
		endif
	elseif(SpaceGroup>=143)				// Trigonal [143,167] (generally hexagonal cell), for rhomohedral use rhomohedral cell, unless obviously the hexagonal
		if ( num90s==2 && abs(xtal.gam - 120)<angleTol && abs(xtal.a - xtal.b)<lenTol )
			return 1							// require a==b, alpha=beta=90, gamma=190
		elseif ( abs(xtal.alpha - xtal.beta)<angleTol && abs(xtal.alpha - xtal.gam)<angleTol )
			if ( abs(xtal.a - xtal.b)<lenTol && abs(xtal.a - xtal.c)<lenTol )	
				return 1						// Rhombohedral axes, require a=b=c
			endif
		endif
	elseif(SpaceGroup>=75)				// Tetragonal [75,142]
		if (num90s==3 && abs(xtal.a - xtal.b)<lenTol)
			return 1							// require a==b, alpha=beta=gamma=90
		endif
	elseif(SpaceGroup>=16)				// Orthorhombic [16,74]
		if (num90s==3)
			return 1							// require alpha=beta=gamma=90
		endif
	elseif(SpaceGroup>=3)					// Monoclinic [3,15]
		if ( abs(xtal.alpha - 90)<angleTol && abs(xtal.gam - 90)<angleTol )
			return 1
		endif
	else										// Triclinic [1,2]
		return 1								// no requirements on lattice constants
	endif

	return 0									// NOT valid
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
	Variable &h,&k,&l							// these hkl are returned with all common factors removed
	Variable maxDiv								// max possible divisor
	Variable i
	maxDiv = max(abs(h),abs(k))				// the maximum divisor cannot be bigger than the smallest of hkl
	maxDiv = max(maxDiv,abs(l))
	for (i=maxDiv;i>=2;i-=1)					// check all divisorts in range [2, maxDiv]
		if (mod(h,i) || mod(k,i) || mod(l,i))	// i is not a factor of h, k, and l
			continue
		endif
		h /= i
		k /= i
		l /= i
	endfor
End


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


ThreadSafe Function/S hkl2str(h,k,l)	// format h,k,l into a string of acceptable minimal length
	Variable h,k,l
	if (numtype(h+k+l))
		return "nan,nan,nan"
	endif
	h = abs(h)<1e-14 ? 0 : h
	k = abs(k)<1e-14 ? 0 : k
	l = abs(l)<1e-14 ? 0 : l
	String hkl
	if (abs(mod(h,1))+abs(mod(k,1))+abs(mod(l,1)) > 1e-6)	// hkl are non-integers
		sprintf hkl,"%g, %g, %g",h,k,l
	elseif (k<0 || l<0)
		sprintf hkl,"%.0f, %.0f, %.0f",h,k,l
	elseif (abs(h)<10 && k<10 && l<10)
		sprintf hkl,"%.0f%.0f%.0f",h,k,l
	else
		sprintf hkl,"%.0f %.0f %.0f",h,k,l
	endif
	return hkl
End


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


ThreadSafe Function/T hkl2IgorBarStr(h,k,l)	// changes negatives to a bar over the number, only for displaying, not printing
	Variable h,k,l				// hkl value

	if (numtype(h+k+l))
		return num2str(h)+","+num2str(k)+","+num2str(l)
	endif
	String extra=SelectString(abs(h)>9 || abs(k)>9 || abs(l)>9 || h<0 || k<0 || l<0,""," ")
	String str=""
	str += minus2bar(num2istr(h),spaces=floor(log(abs(h)))) + extra
	str += minus2bar(num2istr(k),spaces=floor(log(abs(k)))) + extra
	str += minus2bar(num2istr(l),spaces=floor(log(abs(l))))
	return str
End


ThreadSafe Function/T minus2bar(str,[spaces])		// change an Igor string that has minuses to one using a bar over the following character
	String str
	Variable spaces
	spaces = ParamIsDefault(spaces) ? 0 : spaces
	spaces = round(spaces)==limit(spaces,1,5) ? spaces : 0

	String sspaces = PadString("",spaces,0x20), bs
	for (; strsearch(str,"-",0)>=0;)
		sprintf bs, "\\[9\\S\\f01%s\\]9\\M\\X9",sspaces+BCHAR
		str = ReplaceString("-",str,bs,0,1)
	endfor
	return str
End


Function copy_xtal(target,source)						// copy a crystalStructure source to target
	STRUCT crystalStructure &source
	STRUCT crystalStructure &target

	target.desc = source.desc

	target.a = source.a
	target.b = source.b
	target.c = source.c
	target.alpha = source.alpha
	target.beta = source.beta
	target.gam = source.gam

	target.SpaceGroup = source.SpaceGroup
	target.Vc = source.Vc
	target.density = source.density
	target.alphaT = source.alphaT

	target.Temperature = source.Temperature
	target.Vibrate = source.Vibrate
	target.haveDebyeT = source.haveDebyeT
	target.hashID = source.hashID

	target.N = source.N
	Variable i, N=source.N
	for (i=0;i<N;i+=1)
		copy_atomType(target.atom[i],source.atom[i])
	endfor

	target.Nbonds = source.Nbonds
	N = target.Nbonds
	for (i=0;i<N;i+=1)
		copy_bondType(target.bond[i],source.bond[i])
	endfor

	target.a0 = source.a0
	target.b0 = source.b0
	target.c0 = source.c0
	target.a1 = source.a1
	target.b1 = source.b1
	target.c1 = source.c1
	target.a2 = source.a2
	target.b2 = source.b2
	target.c2 = source.c2

	target.as0 = source.as0
	target.bs0 = source.bs0
	target.cs0 = source.cs0
	target.as1 = source.as1
	target.bs1 = source.bs1
	target.cs1 = source.cs1
	target.as2 = source.as2
	target.bs2 = source.bs2
	target.cs2 = source.cs2

	target.Unconventional00 = source.Unconventional00
	target.Unconventional01 = source.Unconventional01
	target.Unconventional02 = source.Unconventional02
	target.Unconventional10 = source.Unconventional10
	target.Unconventional11 = source.Unconventional11
	target.Unconventional12 = source.Unconventional12
	target.Unconventional20 = source.Unconventional20
	target.Unconventional21 = source.Unconventional21
	target.Unconventional22 = source.Unconventional22

	String fullFile = source.sourceFile
	fullFile = fullFile[0,MAX_FILE_LEN-1]
	target.sourceFile = fullFile
End
//
Static Function copy_atomType(target,source)		// copy a atomTypeStructure source to target
	STRUCT atomTypeStructure &source
	STRUCT atomTypeStructure &target

	target.name = source.name
	target.Zatom = source.Zatom
	target.valence = source.valence
	target.x = source.x
	target.y = source.y
	target.z = source.z
	target.occ = source.occ
	target.WyckoffSymbol = source.WyckoffSymbol
	target.DebyeT = source.DebyeT
	target.Biso = source.Biso
	target.Uiso = source.Uiso
	target.U11 = source.U11
	target.U22 = source.U22
	target.U33 = source.U33
	target.U12 = source.U12
	target.U13 = source.U13
	target.U23 = source.U23
End
//
Static Function copy_bondType(target,source)		// copy a bondTypeStructure source to target
	STRUCT bondTypeStructure &source
	STRUCT bondTypeStructure &target

	target.label0 = source.label0
	target.label1 = source.label1
	target.N = source.N
	Variable i,N = source.N
	for (i=0;i<N;i+=1)
		target.len[i] = source.len[i]
	endfor
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
	if (showPanel)
		MakeLatticeParametersPanel("")
	endif
End



//
//	================================================================================
//
//			This section is for making the matricies for the symmetry operations for Space Groups
//



Static Function SetSymOpsForSpaceGroup(SpaceGroup)		// make the symmetry operations mats and vecs (if needed), returns number of operations
	Variable SpaceGroup
	Wave mats = $("root:Packages:Lattices:SymOps:equivXYZM"+num2istr(SpaceGroup))
	Wave bvecs = $("root:Packages:Lattices:SymOps:equivXYZB"+num2istr(SpaceGroup))
	Variable numSymOps
	if (WaveExists(mats) && WaveExists(bvecs))				// check if they exist
		numSymOps = NumberByKey("numSymOps",note(mats),"=")
		return numtype(numSymOps) ? 0 : numSymOps		// do not re-make, just return number of operations
	endif
	if (!isValidSpaceGroup(SpaceGroup))						// Space Group must be in range [1, 230]
		DoAlert 0, "Bad Space Group = "+num2str(SpaceGroup)+", in SetSymOpsForSpaceGroup"
		return 1
	endif
	if (!DataFolderExists("root:Packages:Lattices:SymOps:"))
		DoAlert 0, "Cannot make symmetry matricies, the target data folder does not exist"
		print "Cannot make symmetry matricies, the target data folder does not exist,  'root:Packages:Lattices:SymOps:'\r"
		return 0
	endif

	String symOperations=setSymLine(SpaceGroup)			// a string like "x,y,z;-x,-y,z;-x,y,-z;x,-y,-z;x+1/2,y+1/2,z;-x+1/2,-y+1/2,z;-x+1/2,y+1/2,-z;x+1/2,-y+1/2,-z"
	Variable i,N=ItemsInList(symOperations)
	String wName = "root:Packages:Lattices:SymOps:equivXYZM"+num2istr(SpaceGroup)
	Make/N=(N,3,3)/O/B $wName										// this only holds 0 or ±1
	Wave equivM = $wName
	wName = "root:Packages:Lattices:SymOps:equivXYZB"+num2istr(SpaceGroup)
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
	Note/K equivM, ReplaceNumberByKey("numSymOps","",N,"=")
	Note/K equivB, ReplaceNumberByKey("numSymOps","",N,"=")
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
Static Function ParseOneSymEquation(expression,m0,m1,m2,b)		// parse one expression of form "-x+y"  or "-x", or "-x+y, etc.
	String expression
	Variable &m0,&m1,&m2, &b
	m0=0  ;  m1=0  ;  m2=0 ; b=0	// init, some may not be set in the strswitch
	Variable num=NaN, denom=NaN

	Variable i,is, N=strlen(expression), op
	for (i=0, op=1,is=1; i<N && op; i+=1)
		strswitch(expression[i,i])
			case "-":
				is = -1
				break
			case "+":
				is = 1
				break
			case "x":
				m0 = is
				is = 1
				break
			case "y":
				m1 = is
				is = 1
				break
			case "z":
				m2 = is
				is = 1
				break
			default:					//should be a digit
				op=0					// halt loop, and set pointer back
				num=str2num(expression[i,Inf])
		endswitch
	endfor
	if (op)								// ran out of things to do, no b
		return 0
	endif
	if (char2num("/")!=char2num(expression[i]))
		DoAlert 0, "ERROR"
		print "error on expression = ",expression
		return 1
	endif
	b = is * num/str2num(expression[i+1]) // eval the constant
	if (numtype(b))
		DoAlert 0, "ERROR"
		print "error on expression = ",expression
		return 1
	endif
	return 0							// return no-error
End
//Function testParseOne()
//	String expression = "-x+y-2/3"
//	Variable m0,m1,m2, b
//	ParseOneSymEquation(expression,m0,m1,m2,b)
//	print expression,"    ",m0,"  ",m1,"  ",m2,"  ",b
//	expression = "x-y"
//	ParseOneSymEquation(expression,m0,m1,m2,b)
//	print expression,"    ",m0,"  ",m1,"  ",m2,"  ",b
//	expression = "-z-1/4"
//	ParseOneSymEquation(expression,m0,m1,m2,b)
//	print expression,"    ",m0,"  ",m1,"  ",m2,"  ",b
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
Static Function/T setSymLine(SpaceGroup)
	Variable SpaceGroup								// Space Group number [1,230]
	if (!isValidSpaceGroup(SpaceGroup))		// invalid
		return ""
	endif
	Make/N=230/O/T symLines_temp__
	Wave/T symLines = symLines_temp__
	// Triclinic [1,2]  lines 0-1
	symLines[0]  = "x,y,z"
	symLines[1]  = "x,y,z;-x,-y,-z"
	// Monoclinic [3,15]  lines 2-14
	symLines[2]  = "x,y,z;-x,y,-z"
	symLines[3]  = "x,y,z;-x,y+1/2,-z"
	symLines[4]  = "x,y,z;-x,y,-z;x+1/2,y+1/2,z;-x+1/2,y+1/2,-z"
	symLines[5]  = "x,y,z;x,-y,z"
	symLines[6]  = "x,y,z;x,-y,z+1/2"
	symLines[7]  = "x,y,z;x,-y,z;x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[8]  = "x,y,z;x,-y,z+1/2;x+1/2,y+1/2,z;x+1/2,-y+1/2,z+1/2"
	symLines[9]  = "x,y,z;-x,y,-z;-x,-y,-z;x,-y,z"
	symLines[10]  = "x,y,z;-x,y+1/2,-z;-x,-y,-z;x,-y+1/2,z"
	symLines[11]  = "x,y,z;-x,y,-z;-x,-y,-z;x,-y,z;x+1/2,y+1/2,z;-x+1/2,y+1/2,-z;-x+1/2,-y+1/2,-z;x+1/2,-y+1/2,z"
	symLines[12]  = "x,y,z;-x,y,-z+1/2;-x,-y,-z;x,-y,z+1/2"
	symLines[13]  = "x,y,z;-x,y+1/2,-z+1/2;-x,-y,-z;x,-y+1/2,z+1/2"
	symLines[14]  = "x,y,z;-x,y,-z+1/2;-x,-y,-z;x,-y,z+1/2;x+1/2,y+1/2,z;-x+1/2,y+1/2,-z+1/2;-x+1/2,-y+1/2,-z;x+1/2,-y+1/2,z+1/2"
	// Orthorhombic [16,74]  lines 15-73
	symLines[15]  = "x,y,z;-x,-y,z;x,-y,-z;-x,y,-z"
	symLines[16]  = "x,y,z;-x,-y,z+1/2;x,-y,-z;-x,y,-z+1/2"
	symLines[17]  = "x,y,z;-x,-y,z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z"
	symLines[18]  = "x,y,z;-x+1/2,-y,z+1/2;x+1/2,-y+1/2,-z;-x,y+1/2,-z+1/2"
	symLines[19]  = "x,y,z;-x,-y,z+1/2;x,-y,-z;-x,y,-z+1/2;x+1/2,y+1/2,z;-x+1/2,-y+1/2,z+1/2;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z+1/2"
	symLines[20]  = "x,y,z;-x,-y,z;x,-y,-z;-x,y,-z;x+1/2,y+1/2,z;-x+1/2,-y+1/2,z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z"
	symLines[21]  = "x,y,z;-x,-y,z;x,-y,-z;-x,y,-z;x,y+1/2,z+1/2;-x,-y+1/2,z+1/2;x,-y+1/2,-z+1/2;-x,y+1/2,-z+1/2;x+1/2,y,z+1/2;"
	symLines[21] += "-x+1/2,-y,z+1/2;x+1/2,-y,-z+1/2;-x+1/2,y,-z+1/2;x+1/2,y+1/2,z;-x+1/2,-y+1/2,z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z"
	symLines[22]  = "x,y,z;-x,-y,z;x,-y,-z;-x,y,-z;x+1/2,y+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2"
	symLines[23]  = "x,y,z;-x,-y+1/2,z;x,-y,-z+1/2;-x+1/2,y,-z;x+1/2,y+1/2,z+1/2;-x+1/2,-y,z+1/2;x+1/2,-y+1/2,-z;-x,y+1/2,-z+1/2"
	symLines[24]  = "x,y,z;-x,-y,z;-x,y,z;x,-y,z"
	symLines[25]  = "x,y,z;-x,-y,z+1/2;-x,y,z;x,-y,z+1/2"
	symLines[26]  = "x,y,z;-x,-y,z;-x,y,z+1/2;x,-y,z+1/2"
	symLines[27]  = "x,y,z;-x,-y,z;-x+1/2,y,z;x+1/2,-y,z"
	symLines[28]  = "x,y,z;-x,-y,z+1/2;-x+1/2,y,z+1/2;x+1/2,-y,z"
	symLines[29]  = "x,y,z;-x,-y,z;-x,y+1/2,z+1/2;x,-y+1/2,z+1/2"
	symLines[30]  = "x,y,z;-x+1/2,-y,z+1/2;-x,y,z;x+1/2,-y,z+1/2"
	symLines[31]  = "x,y,z;-x,-y,z;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[32]  = "x,y,z;-x,-y,z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z"
	symLines[33]  = "x,y,z;-x,-y,z;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[34]  = "x,y,z;-x,-y,z;-x,y,z;x,-y,z;x+1/2,y+1/2,z;-x+1/2,-y+1/2,z;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[35]  = "x,y,z;-x,-y,z+1/2;-x,y,z;x,-y,z+1/2;x+1/2,y+1/2,z;-x+1/2,-y+1/2,z+1/2;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z+1/2"
	symLines[36]  = "x,y,z;-x,-y,z;-x,y,z+1/2;x,-y,z+1/2;x+1/2,y+1/2,z;-x+1/2,-y+1/2,z;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[37]  = "x,y,z;-x,-y,z;-x,y,z;x,-y,z;x,y+1/2,z+1/2;-x,-y+1/2,z+1/2;-x,y+1/2,z+1/2;x,-y+1/2,z+1/2"
	symLines[38]  = "x,y,z;-x,-y,z;-x,y,z+1/2;x,-y,z+1/2;x,y+1/2,z+1/2;-x,-y+1/2,z+1/2;-x,y+1/2,z;x,-y+1/2,z"
	symLines[39]  = "x,y,z;-x,-y,z;-x+1/2,y,z;x+1/2,-y,z;x,y+1/2,z+1/2;-x,-y+1/2,z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[40]  = "x,y,z;-x,-y,z;-x+1/2,y,z+1/2;x+1/2,-y,z+1/2;x,y+1/2,z+1/2;-x,-y+1/2,z+1/2;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[41]  = "x,y,z;-x,-y,z;-x,y,z;x,-y,z;x,y+1/2,z+1/2;-x,-y+1/2,z+1/2;-x,y+1/2,z+1/2;x,-y+1/2,z+1/2;x+1/2,y,z+1/2;"
	symLines[41] += "-x+1/2,-y,z+1/2;-x+1/2,y,z+1/2;x+1/2,-y,z+1/2;x+1/2,y+1/2,z;-x+1/2,-y+1/2,z;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[42]  = "x,y,z;-x,-y,z;-x+1/4,y+1/4,z+1/4;x+1/4,-y+1/4,z+1/4;x,y+1/2,z+1/2;-x,-y+1/2,z+1/2;-x+1/4,y+3/4,z+3/4;"
	symLines[42] += "x+1/4,-y+3/4,z+3/4;x+1/2,y,z+1/2;-x+1/2,-y,z+1/2;-x+3/4,y+1/4,z+3/4;x+3/4,-y+1/4,z+3/4;x+1/2,y+1/2,z;"
	symLines[42] += "-x+1/2,-y+1/2,z;-x+3/4,y+3/4,z+1/4;x+3/4,-y+3/4,z+1/4"
	symLines[43]  = "x,y,z;-x,-y,z;-x,y,z;x,-y,z;x+1/2,y+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[44]  = "x,y,z;-x,-y,z;-x,y,z+1/2;x,-y,z+1/2;x+1/2,y+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[45]  = "x,y,z;-x,-y,z;-x+1/2,y,z;x+1/2,-y,z;x+1/2,y+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;-x,y+1/2,z+1/2;x,-y+1/2,z+1/2"
	symLines[46]  = "x,y,z;-x,-y,z;x,-y,-z;-x,y,-z;-x,-y,-z;x,y,-z;-x,y,z;x,-y,z"
	symLines[47]  = "x,y,z;-x+1/2,-y+1/2,-z+1/2;-x,-y,z;x,-y,-z;-x,y,-z;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[48]  = "x,y,z;-x,-y,z;x,-y,-z+1/2;-x,y,-z+1/2;-x,-y,-z;x,y,-z;-x,y,z+1/2;x,-y,z+1/2"
	symLines[49]  = "x,y,z;-x+1/2,-y+1/2,-z;-x,-y,z;x,-y,-z;-x,y,-z;x+1/2,y+1/2,-z;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[50]  = "x,y,z;-x+1/2,-y,z;x+1/2,-y,-z;-x,y,-z;-x,-y,-z;x+1/2,y,-z;-x+1/2,y,z;x,-y,z"
	symLines[51]  = "x,y,z;-x+1/2,-y,z;x,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2;-x,-y,-z;x+1/2,y,-z;-x,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[52]  = "x,y,z;-x+1/2,-y,z+1/2;x,-y,-z;-x+1/2,y,-z+1/2;-x,-y,-z;x+1/2,y,-z+1/2;-x,y,z;x+1/2,-y,z+1/2"
	symLines[53]  = "x,y,z;-x+1/2,-y,z;x+1/2,-y,-z+1/2;-x,y,-z+1/2;-x,-y,-z;x+1/2,y,-z;-x+1/2,y,z+1/2;x,-y,z+1/2"
	symLines[54]  = "x,y,z;-x,-y,z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;-x,-y,-z;x,y,-z;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[55]  = "x,y,z;-x+1/2,-y+1/2,z;x+1/2,-y,-z+1/2;-x,y+1/2,-z+1/2;-x,-y,-z;x+1/2,y+1/2,-z;-x+1/2,y,z+1/2;x,-y+1/2,z+1/2"
	symLines[56]  = "x,y,z;-x,-y,z+1/2;x,-y+1/2,-z;-x,y+1/2,-z+1/2;-x,-y,-z;x,y,-z+1/2;-x,y+1/2,z;x,-y+1/2,z+1/2"
	symLines[57]  = "x,y,z;-x,-y,z;x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2;-x,-y,-z;x,y,-z;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[58]  = "x,y,z;-x+1/2,-y+1/2,-z;-x,-y,z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;x+1/2,y+1/2,-z;-x,y,z;x,-y,z"
	symLines[59]  = "x,y,z;-x+1/2,-y+1/2,z+1/2;x+1/2,-y+1/2,-z;-x,y,-z+1/2;-x,-y,-z;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z;x,-y,z+1/2"
	symLines[60]  = "x,y,z;-x+1/2,-y,z+1/2;x+1/2,-y+1/2,-z;-x,y+1/2,-z+1/2;-x,-y,-z;x+1/2,y,-z+1/2;-x+1/2,y+1/2,z;x,-y+1/2,z+1/2"
	symLines[61]  = "x,y,z;-x+1/2,-y,z+1/2;x+1/2,-y+1/2,-z+1/2;-x,y+1/2,-z;-x,-y,-z;x+1/2,y,-z+1/2;-x+1/2,y+1/2,z+1/2;x,-y+1/2,z"
	symLines[62]  = "x,y,z;-x,-y,z+1/2;x,-y,-z;-x,y,-z+1/2;-x,-y,-z;x,y,-z+1/2;-x,y,z;x,-y,z+1/2;x+1/2,y+1/2,z;-x+1/2,-y+1/2,z+1/2;"
	symLines[62] += "x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z+1/2;-x+1/2,-y+1/2,-z;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z+1/2"
	symLines[63]  = "x,y,z;-x,-y+1/2,z+1/2;x,-y,-z;-x,y+1/2,-z+1/2;-x,-y,-z;x,y+1/2,-z+1/2;-x,y,z;x,-y+1/2,z+1/2;x+1/2,y+1/2,z;"
	symLines[63] += "-x+1/2,-y,z+1/2;x+1/2,-y+1/2,-z;-x+1/2,y,-z+1/2;-x+1/2,-y+1/2,-z;x+1/2,y,-z+1/2;-x+1/2,y+1/2,z;x+1/2,-y,z+1/2"
	symLines[64]  = "x,y,z;-x,-y,z;x,-y,-z;-x,y,-z;-x,-y,-z;x,y,-z;-x,y,z;x,-y,z;x+1/2,y+1/2,z;-x+1/2,-y+1/2,z;x+1/2,-y+1/2,-z;"
	symLines[64] += "-x+1/2,y+1/2,-z;-x+1/2,-y+1/2,-z;x+1/2,y+1/2,-z;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[65]  = "x,y,z;-x,-y,z;x,-y,-z+1/2;-x,y,-z+1/2;-x,-y,-z;x,y,-z;-x,y,z+1/2;x,-y,z+1/2;x+1/2,y+1/2,z;-x+1/2,-y+1/2,z;"
	symLines[65] += "x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2;-x+1/2,-y+1/2,-z;x+1/2,y+1/2,-z;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[66]  = "x,y,z;-x,-y+1/2,z;x,-y,-z;-x,y+1/2,-z;-x,-y,-z;x,y+1/2,-z;-x,y,z;x,-y+1/2,z;x+1/2,y+1/2,z;-x+1/2,-y,z;"
	symLines[66] += "x+1/2,-y+1/2,-z;-x+1/2,y,-z;-x+1/2,-y+1/2,-z;x+1/2,y,-z;-x+1/2,y+1/2,z;x+1/2,-y,z"
	symLines[67]  = "x,y,z;-x,-y+1/2,-z+1/2;-x,-y,z;x,-y,-z;-x,y,-z;x,y+1/2,-z+1/2;-x,y+1/2,z+1/2;x,-y+1/2,z+1/2;x+1/2,y+1/2,z;"
	symLines[67] += "-x+1/2,-y,-z+1/2;-x+1/2,-y+1/2,z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;x+1/2,y,-z+1/2;-x+1/2,y,z+1/2;x+1/2,-y,z+1/2"
	symLines[68]  = "x,y,z;-x,-y,z;x,-y,-z;-x,y,-z;-x,-y,-z;x,y,-z;-x,y,z;x,-y,z;x,y+1/2,z+1/2;-x,-y+1/2,z+1/2;x,-y+1/2,-z+1/2;"
	symLines[68] += "-x,y+1/2,-z+1/2;-x,-y+1/2,-z+1/2;x,y+1/2,-z+1/2;-x,y+1/2,z+1/2;x,-y+1/2,z+1/2;x+1/2,y,z+1/2;-x+1/2,-y,z+1/2;"
	symLines[68] += "x+1/2,-y,-z+1/2;-x+1/2,y,-z+1/2;-x+1/2,-y,-z+1/2;x+1/2,y,-z+1/2;-x+1/2,y,z+1/2;x+1/2,-y,z+1/2;x+1/2,y+1/2,z;"
	symLines[68] += "-x+1/2,-y+1/2,z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;-x+1/2,-y+1/2,-z;x+1/2,y+1/2,-z;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[69]  = "x,y,z;-x+1/4,-y+1/4,-z+1/4;-x,-y,z;x,-y,-z;-x,y,-z;x+1/4,y+1/4,-z+1/4;-x+1/4,y+1/4,z+1/4;x+1/4,-y+1/4,z+1/4;"
	symLines[69] += "x,y+1/2,z+1/2;-x+1/4,-y+3/4,-z+3/4;-x,-y+1/2,z+1/2;x,-y+1/2,-z+1/2;-x,y+1/2,-z+1/2;x+1/4,y+3/4,-z+3/4;"
	symLines[69] += "-x+1/4,y+3/4,z+3/4;x+1/4,-y+3/4,z+3/4;x+1/2,y,z+1/2;-x+3/4,-y+1/4,-z+3/4;-x+1/2,-y,z+1/2;x+1/2,-y,-z+1/2;"
	symLines[69] += "-x+1/2,y,-z+1/2;x+3/4,y+1/4,-z+3/4;-x+3/4,y+1/4,z+3/4;x+3/4,-y+1/4,z+3/4;x+1/2,y+1/2,z;-x+3/4,-y+3/4,-z+1/4;"
	symLines[69] += "-x+1/2,-y+1/2,z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;x+3/4,y+3/4,-z+1/4;-x+3/4,y+3/4,z+1/4;x+3/4,-y+3/4,z+1/4"
	symLines[70]  = "x,y,z;-x,-y,z;x,-y,-z;-x,y,-z;-x,-y,-z;x,y,-z;-x,y,z;x,-y,z;x+1/2,y+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;"
	symLines[70] += "x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2;-x+1/2,-y+1/2,-z+1/2;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2;"
	symLines[70] += "x+1/2,-y+1/2,z+1/2"
	symLines[71]  = "x,y,z;-x,-y,z;x,-y,-z+1/2;-x,y,-z+1/2;-x,-y,-z;x,y,-z;-x,y,z+1/2;x,-y,z+1/2;x+1/2,y+1/2,z+1/2;"
	symLines[71] += "-x+1/2,-y+1/2,z+1/2;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;-x+1/2,-y+1/2,-z+1/2;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z;"
	symLines[71] += "x+1/2,-y+1/2,z"
	symLines[72]  = "x,y,z;-x,-y+1/2,z;x,-y,-z+1/2;-x+1/2,y,-z;-x,-y,-z;x,y+1/2,-z;-x,y,z+1/2;x+1/2,-y,z;x+1/2,y+1/2,z+1/2;"
	symLines[72] += "-x+1/2,-y,z+1/2;x+1/2,-y+1/2,-z;-x,y+1/2,-z+1/2;-x+1/2,-y+1/2,-z+1/2;x+1/2,y,-z+1/2;-x+1/2,y+1/2,z;"
	symLines[72] += "x,-y+1/2,z+1/2"
	symLines[73]  = "x,y,z;-x,-y+1/2,z;x,-y,-z;-x,y+1/2,-z;-x,-y,-z;x,y+1/2,-z;-x,y,z;x,-y+1/2,z;x+1/2,y+1/2,z+1/2;-x+1/2,-y,z+1/2;"
	// Tetragonal [75,142]  lines 74-141
	symLines[73] += "x+1/2,-y+1/2,-z+1/2;-x+1/2,y,-z+1/2;-x+1/2,-y+1/2,-z+1/2;x+1/2,y,-z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y,z+1/2"
	symLines[74]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z"
	symLines[75]  = "x,y,z;-y,x,z+1/4;-x,-y,z+1/2;y,-x,z+3/4"
	symLines[76]  = "x,y,z;-y,x,z+1/2;-x,-y,z;y,-x,z+1/2"
	symLines[77]  = "x,y,z;-y,x,z+3/4;-x,-y,z+1/2;y,-x,z+1/4"
	symLines[78]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x+1/2,y+1/2,z+1/2;-y+1/2,x+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;y+1/2,-x+1/2,z+1/2"
	symLines[79]  = "x,y,z;-y,x+1/2,z+1/4;-x,-y,z;y,-x+1/2,z+1/4;x+1/2,y+1/2,z+1/2;-y+1/2,x,z+3/4;-x+1/2,-y+1/2,z+1/2;y+1/2,-x,z+3/4"
	symLines[80]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z"
	symLines[81]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;x+1/2,y+1/2,z+1/2;y+1/2,-x+1/2,-z+1/2;-x+1/2,-y+1/2,z+1/2;-y+1/2,x+1/2,-z+1/2"
	symLines[82]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;-x,-y,-z;y,-x,-z;x,y,-z;-y,x,-z"
	symLines[83]  = "x,y,z;-y,x,z+1/2;-x,-y,z;y,-x,z+1/2;-x,-y,-z;y,-x,-z+1/2;x,y,-z;-y,x,-z+1/2"
	symLines[84]  = "x,y,z;-x+1/2,-y+1/2,-z;-y+1/2,x+1/2,z;-x,-y,z;y+1/2,-x+1/2,z;y,-x,-z;-y,x,-z;x+1/2,y+1/2,-z"
	symLines[85]  = "x,y,z;-x+1/2,-y+1/2,-z+1/2;-y+1/2,x+1/2,z+1/2;-x,-y,z;y+1/2,-x+1/2,z+1/2;y,-x,-z;-y,x,-z;x+1/2,y+1/2,-z+1/2"
	symLines[86]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;-x,-y,-z;y,-x,-z;x,y,-z;-y,x,-z;x+1/2,y+1/2,z+1/2;-y+1/2,x+1/2,z+1/2;"
	symLines[86] += "-x+1/2,-y+1/2,z+1/2;y+1/2,-x+1/2,z+1/2;-x+1/2,-y+1/2,-z+1/2;y+1/2,-x+1/2,-z+1/2;x+1/2,y+1/2,-z+1/2;"
	symLines[86] += "-y+1/2,x+1/2,-z+1/2"
	symLines[87]  = "x,y,z;-x,-y+1/2,-z+1/4;-y,x+1/2,z+1/4;-x,-y,z;y,-x+1/2,z+1/4;y,-x,-z;-y,x,-z;x,y+1/2,-z+1/4;x+1/2,y+1/2,z+1/2;"
	symLines[87] += "-x+1/2,-y,-z+3/4;-y+1/2,x,z+3/4;-x+1/2,-y+1/2,z+1/2;y+1/2,-x,z+3/4;y+1/2,-x+1/2,-z+1/2;-y+1/2,x+1/2,-z+1/2;"
	symLines[87] += "x+1/2,y,-z+3/4"
	symLines[88]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x,-y,-z;-x,y,-z;y,x,-z;-y,-x,-z"
	symLines[89]  = "x,y,z;-y+1/2,x+1/2,z;-x,-y,z;y+1/2,-x+1/2,z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;y,x,-z;-y,-x,-z"
	symLines[90]  = "x,y,z;-y,x,z+1/4;-x,-y,z+1/2;y,-x,z+3/4;x,-y,-z+1/2;-x,y,-z;y,x,-z+3/4;-y,-x,-z+1/4"
	symLines[91]  = "x,y,z;-y+1/2,x+1/2,z+1/4;-x,-y,z+1/2;y+1/2,-x+1/2,z+3/4;x+1/2,-y+1/2,-z+3/4;-x+1/2,y+1/2,-z+1/4;y,x,-z;-y,-x,-z+1/2"
	symLines[92]  = "x,y,z;-y,x,z+1/2;-x,-y,z;y,-x,z+1/2;x,-y,-z;-x,y,-z;y,x,-z+1/2;-y,-x,-z+1/2"
	symLines[93]  = "x,y,z;-y+1/2,x+1/2,z+1/2;-x,-y,z;y+1/2,-x+1/2,z+1/2;x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2;y,x,-z;-y,-x,-z"
	symLines[94]  = "x,y,z;-y,x,z+3/4;-x,-y,z+1/2;y,-x,z+1/4;x,-y,-z+1/2;-x,y,-z;y,x,-z+1/4;-y,-x,-z+3/4"
	symLines[95]  = "x,y,z;-y+1/2,x+1/2,z+3/4;-x,-y,z+1/2;y+1/2,-x+1/2,z+1/4;x+1/2,-y+1/2,-z+1/4;-x+1/2,y+1/2,-z+3/4;y,x,-z;-y,-x,-z+1/2"
	symLines[96]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x,-y,-z;-x,y,-z;y,x,-z;-y,-x,-z;x+1/2,y+1/2,z+1/2;-y+1/2,x+1/2,z+1/2;"
	symLines[96] += "-x+1/2,-y+1/2,z+1/2;y+1/2,-x+1/2,z+1/2;x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2;y+1/2,x+1/2,-z+1/2;"
	symLines[96] += "-y+1/2,-x+1/2,-z+1/2"
	symLines[97]  = "x,y,z;-y,x+1/2,z+1/4;-x,-y,z;y,-x+1/2,z+1/4;x,-y+1/2,-z+1/4;-x,y+1/2,-z+1/4;y,x,-z;-y,-x,-z;x+1/2,y+1/2,z+1/2;"
	symLines[97] += "-y+1/2,x,z+3/4;-x+1/2,-y+1/2,z+1/2;y+1/2,-x,z+3/4;x+1/2,-y,-z+3/4;-x+1/2,y,-z+3/4;y+1/2,x+1/2,-z+1/2;"
	symLines[97] += "-y+1/2,-x+1/2,-z+1/2"
	symLines[98]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;-x,y,z;x,-y,z;-y,-x,z;y,x,z"
	symLines[99]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z;-y+1/2,-x+1/2,z;y+1/2,x+1/2,z"
	symLines[100]  = "x,y,z;-y,x,z+1/2;-x,-y,z;y,-x,z+1/2;-x,y,z+1/2;x,-y,z+1/2;-y,-x,z;y,x,z"
	symLines[101]  = "x,y,z;-y+1/2,x+1/2,z+1/2;-x,-y,z;y+1/2,-x+1/2,z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2;-y,-x,z;y,x,z"
	symLines[102]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;-x,y,z+1/2;x,-y,z+1/2;-y,-x,z+1/2;y,x,z+1/2"
	symLines[103]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2;-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2"
	symLines[104]  = "x,y,z;-y,x,z+1/2;-x,-y,z;y,-x,z+1/2;-x,y,z;x,-y,z;-y,-x,z+1/2;y,x,z+1/2"
	symLines[105]  = "x,y,z;-y,x,z+1/2;-x,-y,z;y,-x,z+1/2;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z;-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2"
	symLines[106]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;-x,y,z;x,-y,z;-y,-x,z;y,x,z;x+1/2,y+1/2,z+1/2;-y+1/2,x+1/2,z+1/2;"
	symLines[106] += "-x+1/2,-y+1/2,z+1/2;y+1/2,-x+1/2,z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2;-y+1/2,-x+1/2,z+1/2;"
	symLines[106] += "y+1/2,x+1/2,z+1/2"
	symLines[107]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;-x,y,z+1/2;x,-y,z+1/2;-y,-x,z+1/2;y,x,z+1/2;x+1/2,y+1/2,z+1/2;-y+1/2,x+1/2,z+1/2;"
	symLines[107] += "-x+1/2,-y+1/2,z+1/2;y+1/2,-x+1/2,z+1/2;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z;-y+1/2,-x+1/2,z;y+1/2,x+1/2,z"
	symLines[108]  = "x,y,z;-y,x+1/2,z+1/4;-x,-y,z;y,-x+1/2,z+1/4;-x,y,z;x,-y,z;-y,-x+1/2,z+1/4;y,x+1/2,z+1/4;x+1/2,y+1/2,z+1/2;"
	symLines[108] += "-y+1/2,x,z+3/4;-x+1/2,-y+1/2,z+1/2;y+1/2,-x,z+3/4;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2;-y+1/2,-x,z+3/4;"
	symLines[108] += "y+1/2,x,z+3/4"
	symLines[109]  = "x,y,z;-y,x+1/2,z+1/4;-x,-y,z;y,-x+1/2,z+1/4;-x,y,z+1/2;x,-y,z+1/2;-y+1/2,-x,z+1/4;y+1/2,x,z+1/4;"
	symLines[109] += "x+1/2,y+1/2,z+1/2;-y+1/2,x,z+3/4;-x+1/2,-y+1/2,z+1/2;y+1/2,-x,z+3/4;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z;"
	symLines[109] += "-y,-x+1/2,z+3/4;y,x+1/2,z+3/4"
	symLines[110]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;x,-y,-z;-x,y,-z;-y,-x,z;y,x,z"
	symLines[111]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;x,-y,-z+1/2;-x,y,-z+1/2;-y,-x,z+1/2;y,x,z+1/2"
	symLines[112]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;-y+1/2,-x+1/2,z;y+1/2,x+1/2,z"
	symLines[113]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2;-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2"
	symLines[114]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;y,x,-z;-y,-x,-z;-x,y,z;x,-y,z"
	symLines[115]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;y,x,-z+1/2;-y,-x,-z+1/2;-x,y,z+1/2;x,-y,z+1/2"
	symLines[116]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;y+1/2,x+1/2,-z;-y+1/2,-x+1/2,-z;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[117]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;y+1/2,x+1/2,-z+1/2;-y+1/2,-x+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[118]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;y,x,-z;-y,-x,-z;-x,y,z;x,-y,z;x+1/2,y+1/2,z+1/2;y+1/2,-x+1/2,-z+1/2;"
	symLines[118] += "-x+1/2,-y+1/2,z+1/2;-y+1/2,x+1/2,-z+1/2;y+1/2,x+1/2,-z+1/2;-y+1/2,-x+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2;"
	symLines[118] += "x+1/2,-y+1/2,z+1/2"
	symLines[119]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;y,x,-z+1/2;-y,-x,-z+1/2;-x,y,z+1/2;x,-y,z+1/2;x+1/2,y+1/2,z+1/2;"
	symLines[119] += "y+1/2,-x+1/2,-z+1/2;-x+1/2,-y+1/2,z+1/2;-y+1/2,x+1/2,-z+1/2;y+1/2,x+1/2,-z;-y+1/2,-x+1/2,-z;-x+1/2,y+1/2,z;"
	symLines[119] += "x+1/2,-y+1/2,z"
	symLines[120]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;x,-y,-z;-x,y,-z;-y,-x,z;y,x,z;x+1/2,y+1/2,z+1/2;y+1/2,-x+1/2,-z+1/2;"
	symLines[120] += "-x+1/2,-y+1/2,z+1/2;-y+1/2,x+1/2,-z+1/2;x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2;-y+1/2,-x+1/2,z+1/2;"
	symLines[120] += "y+1/2,x+1/2,z+1/2"
	symLines[121]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;x,-y+1/2,-z+1/4;-x,y+1/2,-z+1/4;-y,-x+1/2,z+1/4;y,x+1/2,z+1/4;x+1/2,y+1/2,z+1/2;"
	symLines[121] += "y+1/2,-x+1/2,-z+1/2;-x+1/2,-y+1/2,z+1/2;-y+1/2,x+1/2,-z+1/2;x+1/2,-y,-z+3/4;-x+1/2,y,-z+3/4;-y+1/2,-x,z+3/4;"
	symLines[121] += "y+1/2,x,z+3/4"
	symLines[122]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x,-y,-z;-x,y,-z;y,x,-z;-y,-x,-z;-x,-y,-z;y,-x,-z;x,y,-z;-y,x,-z;-x,y,z;x,-y,z;-y,-x,z;y,x,z"
	symLines[123]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x,-y,-z+1/2;-x,y,-z+1/2;y,x,-z+1/2;-y,-x,-z+1/2;-x,-y,-z;y,-x,-z;x,y,-z;-y,x,-z;"
	symLines[123] += "-x,y,z+1/2;x,-y,z+1/2;-y,-x,z+1/2;y,x,z+1/2"
	symLines[124]  = "x,y,z;-x+1/2,-y+1/2,-z;-y,x,z;-x,-y,z;y,-x,z;y+1/2,-x+1/2,-z;-y+1/2,x+1/2,-z;x,-y,-z;-x,y,-z;y,x,-z;-y,-x,-z;"
	symLines[124] += "x+1/2,y+1/2,-z;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z;-y+1/2,-x+1/2,z;y+1/2,x+1/2,z"
	symLines[125]  = "x,y,z;-x+1/2,-y+1/2,-z+1/2;-y,x,z;-x,-y,z;y,-x,z;y+1/2,-x+1/2,-z+1/2;-y+1/2,x+1/2,-z+1/2;x,-y,-z;-x,y,-z;"
	symLines[125] += "y,x,-z;-y,-x,-z;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2;-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2"
	symLines[126]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;y+1/2,x+1/2,-z;-y+1/2,-x+1/2,-z;-x,-y,-z;y,-x,-z;"
	symLines[126] += "x,y,-z;-y,x,-z;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z;-y+1/2,-x+1/2,z;y+1/2,x+1/2,z"
	symLines[127]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2;y+1/2,x+1/2,-z+1/2;-y+1/2,-x+1/2,-z+1/2;"
	symLines[127] += "-x,-y,-z;y,-x,-z;x,y,-z;-y,x,-z;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2;-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2"
	symLines[128]  = "x,y,z;-x+1/2,-y+1/2,-z;-y+1/2,x+1/2,z;-x,-y,z;y+1/2,-x+1/2,z;y,-x,-z;-y,x,-z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;"
	symLines[128] += "y,x,-z;-y,-x,-z;x+1/2,y+1/2,-z;-x,y,z;x,-y,z;-y+1/2,-x+1/2,z;y+1/2,x+1/2,z"
	symLines[129]  = "x,y,z;-x+1/2,-y+1/2,-z;-y+1/2,x+1/2,z;-x,-y,z;y+1/2,-x+1/2,z;y,-x,-z;-y,x,-z;x+1/2,-y+1/2,-z+1/2;"
	symLines[129] += "-x+1/2,y+1/2,-z+1/2;y,x,-z+1/2;-y,-x,-z+1/2;x+1/2,y+1/2,-z;-x,y,z+1/2;x,-y,z+1/2;-y+1/2,-x+1/2,z+1/2;"
	symLines[129] += "y+1/2,x+1/2,z+1/2"
	symLines[130]  = "x,y,z;-y,x,z+1/2;-x,-y,z;y,-x,z+1/2;x,-y,-z;-x,y,-z;y,x,-z+1/2;-y,-x,-z+1/2;-x,-y,-z;y,-x,-z+1/2;x,y,-z;"
	symLines[130] += "-y,x,-z+1/2;-x,y,z;x,-y,z;-y,-x,z+1/2;y,x,z+1/2"
	symLines[131]  = "x,y,z;-y,x,z+1/2;-x,-y,z;y,-x,z+1/2;x,-y,-z+1/2;-x,y,-z+1/2;y,x,-z;-y,-x,-z;-x,-y,-z;y,-x,-z+1/2;x,y,-z;"
	symLines[131] += "-y,x,-z+1/2;-x,y,z+1/2;x,-y,z+1/2;-y,-x,z;y,x,z"
	symLines[132]  = "x,y,z;-x+1/2,-y+1/2,-z+1/2;-y+1/2,x+1/2,z+1/2;-x,-y,z;y+1/2,-x+1/2,z+1/2;y,-x,-z;-y,x,-z;x,-y,-z+1/2;"
	symLines[132] += "-x,y,-z+1/2;y+1/2,x+1/2,-z;-y+1/2,-x+1/2,-z;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z;-y,-x,z+1/2;"
	symLines[132] += "y,x,z+1/2"
	symLines[133]  = "x,y,z;-x+1/2,-y+1/2,-z+1/2;-y+1/2,x+1/2,z+1/2;-x,-y,z;y+1/2,-x+1/2,z+1/2;y,-x,-z;-y,x,-z;x,-y,-z;-x,y,-z;"
	symLines[133] += "y+1/2,x+1/2,-z+1/2;-y+1/2,-x+1/2,-z+1/2;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2;-y,-x,z;y,x,z"
	symLines[134]  = "x,y,z;-y,x,z+1/2;-x,-y,z;y,-x,z+1/2;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;y+1/2,x+1/2,-z+1/2;-y+1/2,-x+1/2,-z+1/2;"
	symLines[134] += "-x,-y,-z;y,-x,-z+1/2;x,y,-z;-y,x,-z+1/2;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z;-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2"
	symLines[135]  = "x,y,z;-y+1/2,x+1/2,z+1/2;-x,-y,z;y+1/2,-x+1/2,z+1/2;x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2;y,x,-z;-y,-x,-z;"
	symLines[135] += "-x,-y,-z;y+1/2,-x+1/2,-z+1/2;x,y,-z;-y+1/2,x+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2;-y,-x,z;y,x,z"
	symLines[136]  = "x,y,z;-x+1/2,-y+1/2,-z+1/2;-y+1/2,x+1/2,z+1/2;-x,-y,z;y+1/2,-x+1/2,z+1/2;y,-x,-z;-y,x,-z;x+1/2,-y+1/2,-z+1/2;"
	symLines[136] += "-x+1/2,y+1/2,-z+1/2;y,x,-z;-y,-x,-z;x+1/2,y+1/2,-z+1/2;-x,y,z;x,-y,z;-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2"
	symLines[137]  = "x,y,z;-x+1/2,-y+1/2,-z+1/2;-y+1/2,x+1/2,z+1/2;-x,-y,z;y+1/2,-x+1/2,z+1/2;y,-x,-z;-y,x,-z;x+1/2,-y+1/2,-z;"
	symLines[137] += "-x+1/2,y+1/2,-z;y,x,-z+1/2;-y,-x,-z+1/2;x+1/2,y+1/2,-z+1/2;-x,y,z+1/2;x,-y,z+1/2;-y+1/2,-x+1/2,z;y+1/2,x+1/2,z"
	symLines[138]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x,-y,-z;-x,y,-z;y,x,-z;-y,-x,-z;-x,-y,-z;y,-x,-z;x,y,-z;-y,x,-z;-x,y,z;x,-y,z;"
	symLines[138] += "-y,-x,z;y,x,z;x+1/2,y+1/2,z+1/2;-y+1/2,x+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;y+1/2,-x+1/2,z+1/2;x+1/2,-y+1/2,-z+1/2;"
	symLines[138] += "-x+1/2,y+1/2,-z+1/2;y+1/2,x+1/2,-z+1/2;-y+1/2,-x+1/2,-z+1/2;-x+1/2,-y+1/2,-z+1/2;y+1/2,-x+1/2,-z+1/2;"
	symLines[138] += "x+1/2,y+1/2,-z+1/2;-y+1/2,x+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2;-y+1/2,-x+1/2,z+1/2;"
	symLines[138] += "y+1/2,x+1/2,z+1/2"
	symLines[139]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x,-y,-z+1/2;-x,y,-z+1/2;y,x,-z+1/2;-y,-x,-z+1/2;-x,-y,-z;y,-x,-z;x,y,-z;-y,x,-z;"
	symLines[139] += "-x,y,z+1/2;x,-y,z+1/2;-y,-x,z+1/2;y,x,z+1/2;x+1/2,y+1/2,z+1/2;-y+1/2,x+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;"
	symLines[139] += "y+1/2,-x+1/2,z+1/2;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;y+1/2,x+1/2,-z;-y+1/2,-x+1/2,-z;-x+1/2,-y+1/2,-z+1/2;"
	symLines[139] += "y+1/2,-x+1/2,-z+1/2;x+1/2,y+1/2,-z+1/2;-y+1/2,x+1/2,-z+1/2;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z;-y+1/2,-x+1/2,z;"
	symLines[139] += "y+1/2,x+1/2,z"
	symLines[140]  = "x,y,z;-x,-y+1/2,-z+1/4;-y,x+1/2,z+1/4;-x,-y,z;y,-x+1/2,z+1/4;y,-x,-z;-y,x,-z;x,-y+1/2,-z+1/4;-x,y+1/2,-z+1/4;"
	symLines[140] += "y,x,-z;-y,-x,-z;x,y+1/2,-z+1/4;-x,y,z;x,-y,z;-y,-x+1/2,z+1/4;y,x+1/2,z+1/4;x+1/2,y+1/2,z+1/2;-x+1/2,-y,-z+3/4;"
	symLines[140] += "-y+1/2,x,z+3/4;-x+1/2,-y+1/2,z+1/2;y+1/2,-x,z+3/4;y+1/2,-x+1/2,-z+1/2;-y+1/2,x+1/2,-z+1/2;x+1/2,-y,-z+3/4;"
	symLines[140] += "-x+1/2,y,-z+3/4;y+1/2,x+1/2,-z+1/2;-y+1/2,-x+1/2,-z+1/2;x+1/2,y,-z+3/4;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2;"
	symLines[140] += "-y+1/2,-x,z+3/4;y+1/2,x,z+3/4"
	symLines[141]  = "x,y,z;-x,-y+1/2,-z+1/4;-y,x+1/2,z+1/4;-x,-y,z;y,-x+1/2,z+1/4;y,-x,-z;-y,x,-z;x+1/2,-y,-z+1/4;-x+1/2,y,-z+1/4;"
	symLines[141] += "y,x,-z+1/2;-y,-x,-z+1/2;x,y+1/2,-z+1/4;-x,y,z+1/2;x,-y,z+1/2;-y+1/2,-x,z+1/4;y+1/2,x,z+1/4;x+1/2,y+1/2,z+1/2;"
	symLines[141] += "-x+1/2,-y,-z+3/4;-y+1/2,x,z+3/4;-x+1/2,-y+1/2,z+1/2;y+1/2,-x,z+3/4;y+1/2,-x+1/2,-z+1/2;-y+1/2,x+1/2,-z+1/2;"
	symLines[141] += "x,-y+1/2,-z+3/4;-x,y+1/2,-z+3/4;y+1/2,x+1/2,-z;-y+1/2,-x+1/2,-z;x+1/2,y,-z+3/4;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z;"
	symLines[141] += "-y,-x+1/2,z+3/4;y,x+1/2,z+3/4"
	// Trigonal [143,167]  lines 142-166
	symLines[142]  = "x,y,z;-y,x-y,z;-x+y,-x,z"
	symLines[143]  = "x,y,z;-y,x-y,z+1/3;-x+y,-x,z+2/3"
	symLines[144]  = "x,y,z;-y,x-y,z+2/3;-x+y,-x,z+1/3"
	symLines[145]  = "x,y,z;-y,x-y,z;-x+y,-x,z;x+2/3,y+1/3,z+1/3;-y+2/3,x-y+1/3,z+1/3;-x+y+2/3,-x+1/3,z+1/3;x+1/3,y+2/3,z+2/3;"
	symLines[145] += "-y+1/3,x-y+2/3,z+2/3;-x+y+1/3,-x+2/3,z+2/3"
	symLines[146]  = "x,y,z;-y,x-y,z;-x+y,-x,z;-x,-y,-z;y,-x+y,-z;x-y,x,-z"
	symLines[147]  = "x,y,z;-y,x-y,z;-x+y,-x,z;-x,-y,-z;y,-x+y,-z;x-y,x,-z;x+2/3,y+1/3,z+1/3;-y+2/3,x-y+1/3,z+1/3;"
	symLines[147] += "-x+y+2/3,-x+1/3,z+1/3;-x+2/3,-y+1/3,-z+1/3;y+2/3,-x+y+1/3,-z+1/3;x-y+2/3,x+1/3,-z+1/3;x+1/3,y+2/3,z+2/3;"
	symLines[147] += "-y+1/3,x-y+2/3,z+2/3;-x+y+1/3,-x+2/3,z+2/3;-x+1/3,-y+2/3,-z+2/3;y+1/3,-x+y+2/3,-z+2/3;x-y+1/3,x+2/3,-z+2/3"
	symLines[148]  = "x,y,z;-y,x-y,z;-x+y,-x,z;-y,-x,-z;-x+y,y,-z;x,x-y,-z"
	symLines[149]  = "x,y,z;-y,x-y,z;-x+y,-x,z;x-y,-y,-z;-x,-x+y,-z;y,x,-z"
	symLines[150]  = "x,y,z;-y,x-y,z+1/3;-x+y,-x,z+2/3;-y,-x,-z+2/3;-x+y,y,-z+1/3;x,x-y,-z"
	symLines[151]  = "x,y,z;-y,x-y,z+1/3;-x+y,-x,z+2/3;x-y,-y,-z+2/3;-x,-x+y,-z+1/3;y,x,-z"
	symLines[152]  = "x,y,z;-y,x-y,z+2/3;-x+y,-x,z+1/3;-y,-x,-z+1/3;-x+y,y,-z+2/3;x,x-y,-z"
	symLines[153]  = "x,y,z;-y,x-y,z+2/3;-x+y,-x,z+1/3;x-y,-y,-z+1/3;-x,-x+y,-z+2/3;y,x,-z"
	symLines[154]  = "x,y,z;-y,x-y,z;-x+y,-x,z;x-y,-y,-z;-x,-x+y,-z;y,x,-z;x+2/3,y+1/3,z+1/3;-y+2/3,x-y+1/3,z+1/3;"
	symLines[154] += "-x+y+2/3,-x+1/3,z+1/3;x-y+2/3,-y+1/3,-z+1/3;-x+2/3,-x+y+1/3,-z+1/3;y+2/3,x+1/3,-z+1/3;x+1/3,y+2/3,z+2/3;"
	symLines[154] += "-y+1/3,x-y+2/3,z+2/3;-x+y+1/3,-x+2/3,z+2/3;x-y+1/3,-y+2/3,-z+2/3;-x+1/3,-x+y+2/3,-z+2/3;y+1/3,x+2/3,-z+2/3"
	symLines[155]  = "x,y,z;-y,x-y,z;-x+y,-x,z;-x+y,y,z;x,x-y,z;-y,-x,z"
	symLines[156]  = "x,y,z;-y,x-y,z;-x+y,-x,z;y,x,z;x-y,-y,z;-x,-x+y,z"
	symLines[157]  = "x,y,z;-y,x-y,z;-x+y,-x,z;-x+y,y,z+1/2;x,x-y,z+1/2;-y,-x,z+1/2"
	symLines[158]  = "x,y,z;-y,x-y,z;-x+y,-x,z;y,x,z+1/2;x-y,-y,z+1/2;-x,-x+y,z+1/2"
	symLines[159]  = "x,y,z;-y,x-y,z;-x+y,-x,z;-x+y,y,z;x,x-y,z;-y,-x,z;x+2/3,y+1/3,z+1/3;-y+2/3,x-y+1/3,z+1/3;-x+y+2/3,-x+1/3,z+1/3;"
	symLines[159] += "-x+y+2/3,y+1/3,z+1/3;x+2/3,x-y+1/3,z+1/3;-y+2/3,-x+1/3,z+1/3;x+1/3,y+2/3,z+2/3;-y+1/3,x-y+2/3,z+2/3;"
	symLines[159] += "-x+y+1/3,-x+2/3,z+2/3;-x+y+1/3,y+2/3,z+2/3;x+1/3,x-y+2/3,z+2/3;-y+1/3,-x+2/3,z+2/3"
	symLines[160]  = "x,y,z;-y,x-y,z;-x+y,-x,z;-x+y,y,z+1/2;x,x-y,z+1/2;-y,-x,z+1/2;x+2/3,y+1/3,z+1/3;-y+2/3,x-y+1/3,z+1/3;"
	symLines[160] += "-x+y+2/3,-x+1/3,z+1/3;-x+y+2/3,y+1/3,z+5/6;x+2/3,x-y+1/3,z+5/6;-y+2/3,-x+1/3,z+5/6;x+1/3,y+2/3,z+2/3;"
	symLines[160] += "-y+1/3,x-y+2/3,z+2/3;-x+y+1/3,-x+2/3,z+2/3;-x+y+1/3,y+2/3,z+1/6;x+1/3,x-y+2/3,z+1/6;-y+1/3,-x+2/3,z+1/6"
	symLines[161]  = "x,y,z;-y,x-y,z;-x+y,-x,z;-y,-x,-z;-x+y,y,-z;x,x-y,-z;-x,-y,-z;y,-x+y,-z;x-y,x,-z;y,x,z;x-y,-y,z;-x,-x+y,z"
	symLines[162]  = "x,y,z;-y,x-y,z;-x+y,-x,z;-y,-x,-z+1/2;-x+y,y,-z+1/2;x,x-y,-z+1/2;-x,-y,-z;y,-x+y,-z;x-y,x,-z;y,x,z+1/2;"
	symLines[162] += "x-y,-y,z+1/2;-x,-x+y,z+1/2"
	symLines[163]  = "x,y,z;-y,x-y,z;-x+y,-x,z;x-y,-y,-z;-x,-x+y,-z;y,x,-z;-x,-y,-z;y,-x+y,-z;x-y,x,-z;-x+y,y,z;x,x-y,z;-y,-x,z"
	symLines[164]  = "x,y,z;-y,x-y,z;-x+y,-x,z;x-y,-y,-z+1/2;-x,-x+y,-z+1/2;y,x,-z+1/2;-x,-y,-z;y,-x+y,-z;x-y,x,-z;-x+y,y,z+1/2;"
	symLines[164] += "x,x-y,z+1/2;-y,-x,z+1/2"
	symLines[165]  = "x,y,z;-y,x-y,z;-x+y,-x,z;x-y,-y,-z;-x,-x+y,-z;y,x,-z;-x,-y,-z;y,-x+y,-z;x-y,x,-z;-x+y,y,z;x,x-y,z;-y,-x,z;"
	symLines[165] += "x+2/3,y+1/3,z+1/3;-y+2/3,x-y+1/3,z+1/3;-x+y+2/3,-x+1/3,z+1/3;x-y+2/3,-y+1/3,-z+1/3;-x+2/3,-x+y+1/3,-z+1/3;"
	symLines[165] += "y+2/3,x+1/3,-z+1/3;-x+2/3,-y+1/3,-z+1/3;y+2/3,-x+y+1/3,-z+1/3;x-y+2/3,x+1/3,-z+1/3;-x+y+2/3,y+1/3,z+1/3;"
	symLines[165] += "x+2/3,x-y+1/3,z+1/3;-y+2/3,-x+1/3,z+1/3;x+1/3,y+2/3,z+2/3;-y+1/3,x-y+2/3,z+2/3;-x+y+1/3,-x+2/3,z+2/3;"
	symLines[165] += "x-y+1/3,-y+2/3,-z+2/3;-x+1/3,-x+y+2/3,-z+2/3;y+1/3,x+2/3,-z+2/3;-x+1/3,-y+2/3,-z+2/3;y+1/3,-x+y+2/3,-z+2/3;"
	symLines[165] += "x-y+1/3,x+2/3,-z+2/3;-x+y+1/3,y+2/3,z+2/3;x+1/3,x-y+2/3,z+2/3;-y+1/3,-x+2/3,z+2/3"
	symLines[166]  = "x,y,z;-y,x-y,z;-x+y,-x,z;x-y,-y,-z+1/2;-x,-x+y,-z+1/2;y,x,-z+1/2;-x,-y,-z;y,-x+y,-z;x-y,x,-z;-x+y,y,z+1/2;"
	symLines[166] += "x,x-y,z+1/2;-y,-x,z+1/2;x+2/3,y+1/3,z+1/3;-y+2/3,x-y+1/3,z+1/3;-x+y+2/3,-x+1/3,z+1/3;x-y+2/3,-y+1/3,-z+5/6;"
	symLines[166] += "-x+2/3,-x+y+1/3,-z+5/6;y+2/3,x+1/3,-z+5/6;-x+2/3,-y+1/3,-z+1/3;y+2/3,-x+y+1/3,-z+1/3;x-y+2/3,x+1/3,-z+1/3;"
	symLines[166] += "-x+y+2/3,y+1/3,z+5/6;x+2/3,x-y+1/3,z+5/6;-y+2/3,-x+1/3,z+5/6;x+1/3,y+2/3,z+2/3;-y+1/3,x-y+2/3,z+2/3;"
	symLines[166] += "-x+y+1/3,-x+2/3,z+2/3;x-y+1/3,-y+2/3,-z+1/6;-x+1/3,-x+y+2/3,-z+1/6;y+1/3,x+2/3,-z+1/6;-x+1/3,-y+2/3,-z+2/3;"
	symLines[166] += "y+1/3,-x+y+2/3,-z+2/3;x-y+1/3,x+2/3,-z+2/3;-x+y+1/3,y+2/3,z+1/6;x+1/3,x-y+2/3,z+1/6;-y+1/3,-x+2/3,z+1/6"
	// Hexagonal [168,194]  lines 167-193
	symLines[167]  = "x,y,z;x-y,x,z;-y,x-y,z;-x,-y,z;-x+y,-x,z;y,-x+y,z"
	symLines[168]  = "x,y,z;x-y,x,z+1/6;-y,x-y,z+1/3;-x,-y,z+1/2;-x+y,-x,z+2/3;y,-x+y,z+5/6"
	symLines[169]  = "x,y,z;x-y,x,z+5/6;-y,x-y,z+2/3;-x,-y,z+1/2;-x+y,-x,z+1/3;y,-x+y,z+1/6"
	symLines[170]  = "x,y,z;x-y,x,z+1/3;-y,x-y,z+2/3;-x,-y,z;-x+y,-x,z+1/3;y,-x+y,z+2/3"
	symLines[171]  = "x,y,z;x-y,x,z+2/3;-y,x-y,z+1/3;-x,-y,z;-x+y,-x,z+2/3;y,-x+y,z+1/3"
	symLines[172]  = "x,y,z;x-y,x,z+1/2;-y,x-y,z;-x,-y,z+1/2;-x+y,-x,z;y,-x+y,z+1/2"
	symLines[173]  = "x,y,z;-x+y,-x,-z;-y,x-y,z;x,y,-z;-x+y,-x,z;-y,x-y,-z"
	symLines[174]  = "x,y,z;x-y,x,z;-y,x-y,z;-x,-y,z;-x+y,-x,z;y,-x+y,z;-x,-y,-z;-x+y,-x,-z;y,-x+y,-z;x,y,-z;x-y,x,-z;-y,x-y,-z"
	symLines[175]  = "x,y,z;x-y,x,z+1/2;-y,x-y,z;-x,-y,z+1/2;-x+y,-x,z;y,-x+y,z+1/2;-x,-y,-z;-x+y,-x,-z+1/2;y,-x+y,-z;x,y,-z+1/2;"
	symLines[175] += "x-y,x,-z;-y,x-y,-z+1/2"
	symLines[176]  = "x,y,z;x-y,x,z;-y,x-y,z;-x,-y,z;-x+y,-x,z;y,-x+y,z;x-y,-y,-z;-x,-x+y,-z;y,x,-z;-y,-x,-z;-x+y,y,-z;x,x-y,-z"
	symLines[177]  = "x,y,z;x-y,x,z+1/6;-y,x-y,z+1/3;-x,-y,z+1/2;-x+y,-x,z+2/3;y,-x+y,z+5/6;x-y,-y,-z;-x,-x+y,-z+2/3;y,x,-z+1/3;"
	symLines[177] += "-y,-x,-z+5/6;-x+y,y,-z+1/2;x,x-y,-z+1/6"
	symLines[178]  = "x,y,z;x-y,x,z+5/6;-y,x-y,z+2/3;-x,-y,z+1/2;-x+y,-x,z+1/3;y,-x+y,z+1/6;x-y,-y,-z;-x,-x+y,-z+1/3;y,x,-z+2/3;"
	symLines[178] += "-y,-x,-z+1/6;-x+y,y,-z+1/2;x,x-y,-z+5/6"
	symLines[179]  = "x,y,z;x-y,x,z+1/3;-y,x-y,z+2/3;-x,-y,z;-x+y,-x,z+1/3;y,-x+y,z+2/3;x-y,-y,-z;-x,-x+y,-z+1/3;y,x,-z+2/3;"
	symLines[179] += "-y,-x,-z+2/3;-x+y,y,-z;x,x-y,-z+1/3"
	symLines[180]  = "x,y,z;x-y,x,z+2/3;-y,x-y,z+1/3;-x,-y,z;-x+y,-x,z+2/3;y,-x+y,z+1/3;x-y,-y,-z;-x,-x+y,-z+2/3;y,x,-z+1/3;"
	symLines[180] += "-y,-x,-z+1/3;-x+y,y,-z;x,x-y,-z+2/3"
	symLines[181]  = "x,y,z;x-y,x,z+1/2;-y,x-y,z;-x,-y,z+1/2;-x+y,-x,z;y,-x+y,z+1/2;x-y,-y,-z;-x,-x+y,-z;y,x,-z;-y,-x,-z+1/2;"
	symLines[181] += "-x+y,y,-z+1/2;x,x-y,-z+1/2"
	symLines[182]  = "x,y,z;x-y,x,z;-y,x-y,z;-x,-y,z;-x+y,-x,z;y,-x+y,z;-x+y,y,z;x,x-y,z;-y,-x,z;y,x,z;x-y,-y,z;-x,-x+y,z"
	symLines[183]  = "x,y,z;x-y,x,z;-y,x-y,z;-x,-y,z;-x+y,-x,z;y,-x+y,z;-x+y,y,z+1/2;x,x-y,z+1/2;-y,-x,z+1/2;y,x,z+1/2;x-y,-y,z+1/2;-x,-x+y,z+1/2"
	symLines[184]  = "x,y,z;x-y,x,z+1/2;-y,x-y,z;-x,-y,z+1/2;-x+y,-x,z;y,-x+y,z+1/2;-x+y,y,z+1/2;x,x-y,z+1/2;-y,-x,z+1/2;y,x,z;x-y,-y,z;-x,-x+y,z"
	symLines[185]  = "x,y,z;x-y,x,z+1/2;-y,x-y,z;-x,-y,z+1/2;-x+y,-x,z;y,-x+y,z+1/2;-x+y,y,z;x,x-y,z;-y,-x,z;y,x,z+1/2;x-y,-y,z+1/2;-x,-x+y,z+1/2"
	symLines[186]  = "x,y,z;-x+y,-x,-z;-y,x-y,z;x,y,-z;-x+y,-x,z;-y,x-y,-z;-y,-x,-z;-x+y,y,-z;x,x-y,-z;-x+y,y,z;x,x-y,z;-y,-x,z"
	symLines[187]  = "x,y,z;-x+y,-x,-z+1/2;-y,x-y,z;x,y,-z+1/2;-x+y,-x,z;-y,x-y,-z+1/2;-y,-x,-z;-x+y,y,-z;x,x-y,-z;-x+y,y,z+1/2;"
	symLines[187] += "x,x-y,z+1/2;-y,-x,z+1/2"
	symLines[188]  = "x,y,z;-x+y,-x,-z;-y,x-y,z;x,y,-z;-x+y,-x,z;-y,x-y,-z;x-y,-y,-z;-x,-x+y,-z;y,x,-z;y,x,z;x-y,-y,z;-x,-x+y,z"
	symLines[189]  = "x,y,z;-x+y,-x,-z+1/2;-y,x-y,z;x,y,-z+1/2;-x+y,-x,z;-y,x-y,-z+1/2;x-y,-y,-z;-x,-x+y,-z;y,x,-z;y,x,z+1/2;"
	symLines[189] += "x-y,-y,z+1/2;-x,-x+y,z+1/2"
	symLines[190]  = "x,y,z;x-y,x,z;-y,x-y,z;-x,-y,z;-x+y,-x,z;y,-x+y,z;x-y,-y,-z;-x,-x+y,-z;y,x,-z;-y,-x,-z;-x+y,y,-z;x,x-y,-z;"
	symLines[190] += "-x,-y,-z;-x+y,-x,-z;y,-x+y,-z;x,y,-z;x-y,x,-z;-y,x-y,-z;-x+y,y,z;x,x-y,z;-y,-x,z;y,x,z;x-y,-y,z;-x,-x+y,z"
	symLines[191]  = "x,y,z;x-y,x,z;-y,x-y,z;-x,-y,z;-x+y,-x,z;y,-x+y,z;x-y,-y,-z+1/2;-x,-x+y,-z+1/2;y,x,-z+1/2;-y,-x,-z+1/2;"
	symLines[191] += "-x+y,y,-z+1/2;x,x-y,-z+1/2;-x,-y,-z;-x+y,-x,-z;y,-x+y,-z;x,y,-z;x-y,x,-z;-y,x-y,-z;-x+y,y,z+1/2;x,x-y,z+1/2;"
	symLines[191] += "-y,-x,z+1/2;y,x,z+1/2;x-y,-y,z+1/2;-x,-x+y,z+1/2"
	symLines[192]  = "x,y,z;x-y,x,z+1/2;-y,x-y,z;-x,-y,z+1/2;-x+y,-x,z;y,-x+y,z+1/2;x-y,-y,-z+1/2;-x,-x+y,-z+1/2;y,x,-z+1/2;-y,-x,-z;"
	symLines[192] += "-x+y,y,-z;x,x-y,-z;-x,-y,-z;-x+y,-x,-z+1/2;y,-x+y,-z;x,y,-z+1/2;x-y,x,-z;-y,x-y,-z+1/2;-x+y,y,z+1/2;"
	symLines[192] += "x,x-y,z+1/2;-y,-x,z+1/2;y,x,z;x-y,-y,z;-x,-x+y,z"
	symLines[193]  = "x,y,z;x-y,x,z+1/2;-y,x-y,z;-x,-y,z+1/2;-x+y,-x,z;y,-x+y,z+1/2;x-y,-y,-z;-x,-x+y,-z;y,x,-z;-y,-x,-z+1/2;"
	symLines[193] += "-x+y,y,-z+1/2;x,x-y,-z+1/2;-x,-y,-z;-x+y,-x,-z+1/2;y,-x+y,-z;x,y,-z+1/2;x-y,x,-z;-y,x-y,-z+1/2;-x+y,y,z;"
	symLines[193] += "x,x-y,z;-y,-x,z;y,x,z+1/2;x-y,-y,z+1/2;-x,-x+y,z+1/2"
	// Cubic [195,230]  lines 194-229
	symLines[194]  = "x,y,z;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;-x,-y,z;x,-y,-z;-x,y,-z"
	symLines[195]  = "x,y,z;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;-x,-y,z;x,-y,-z;-x,y,-z;x,y+1/2,z+1/2;"
	symLines[195] += "z,x+1/2,y+1/2;y,z+1/2,x+1/2;-y,-z+1/2,x+1/2;z,-x+1/2,-y+1/2;-y,z+1/2,-x+1/2;-z,-x+1/2,y+1/2;-z,x+1/2,-y+1/2;"
	symLines[195] += "y,-z+1/2,-x+1/2;-x,-y+1/2,z+1/2;x,-y+1/2,-z+1/2;-x,y+1/2,-z+1/2;x+1/2,y,z+1/2;z+1/2,x,y+1/2;y+1/2,z,x+1/2;"
	symLines[195] += "-y+1/2,-z,x+1/2;z+1/2,-x,-y+1/2;-y+1/2,z,-x+1/2;-z+1/2,-x,y+1/2;-z+1/2,x,-y+1/2;y+1/2,-z,-x+1/2;"
	symLines[195] += "-x+1/2,-y,z+1/2;x+1/2,-y,-z+1/2;-x+1/2,y,-z+1/2;x+1/2,y+1/2,z;z+1/2,x+1/2,y;y+1/2,z+1/2,x;-y+1/2,-z+1/2,x;"
	symLines[195] += "z+1/2,-x+1/2,-y;-y+1/2,z+1/2,-x;-z+1/2,-x+1/2,y;-z+1/2,x+1/2,-y;y+1/2,-z+1/2,-x;-x+1/2,-y+1/2,z;"
	symLines[195] += "x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z"
	symLines[196]  = "x,y,z;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;-x,-y,z;x,-y,-z;-x,y,-z;x+1/2,y+1/2,z+1/2;"
	symLines[196] += "z+1/2,x+1/2,y+1/2;y+1/2,z+1/2,x+1/2;-y+1/2,-z+1/2,x+1/2;z+1/2,-x+1/2,-y+1/2;-y+1/2,z+1/2,-x+1/2;"
	symLines[196] += "-z+1/2,-x+1/2,y+1/2;-z+1/2,x+1/2,-y+1/2;y+1/2,-z+1/2,-x+1/2;-x+1/2,-y+1/2,z+1/2;x+1/2,-y+1/2,-z+1/2;"
	symLines[196] += "-x+1/2,y+1/2,-z+1/2"
	symLines[197]  = "x,y,z;z,x,y;y,z,x;-y+1/2,-z,x+1/2;z+1/2,-x+1/2,-y;-y,z+1/2,-x+1/2;-z+1/2,-x,y+1/2;-z,x+1/2,-y+1/2;"
	symLines[197] += "y+1/2,-z+1/2,-x;-x+1/2,-y,z+1/2;x+1/2,-y+1/2,-z;-x,y+1/2,-z+1/2"
	symLines[198]  = "x,y,z;z,x,y;y,z,x;-y,-z+1/2,x;z,-x,-y+1/2;-y+1/2,z,-x;-z,-x+1/2,y;-z+1/2,x,-y;y,-z,-x+1/2;-x,-y+1/2,z;"
	symLines[198] += "x,-y,-z+1/2;-x+1/2,y,-z;x+1/2,y+1/2,z+1/2;z+1/2,x+1/2,y+1/2;y+1/2,z+1/2,x+1/2;-y+1/2,-z,x+1/2;z+1/2,-x+1/2,-y;"
	symLines[198] += "-y,z+1/2,-x+1/2;-z+1/2,-x,y+1/2;-z,x+1/2,-y+1/2;y+1/2,-z+1/2,-x;-x+1/2,-y,z+1/2;x+1/2,-y+1/2,-z;-x,y+1/2,-z+1/2"
	symLines[199]  = "x,y,z;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;-x,-y,z;x,-y,-z;-x,y,-z;-x,-y,-z;-z,-x,-y;"
	symLines[199] += "-y,-z,-x;y,z,-x;-z,x,y;y,-z,x;z,x,-y;z,-x,y;-y,z,x;x,y,-z;-x,y,z;x,-y,z"
	symLines[200]  = "x,y,z;-x+1/2,-y+1/2,-z+1/2;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;-z+1/2,-x+1/2,-y+1/2;"
	symLines[200] += "-y+1/2,-z+1/2,-x+1/2;y+1/2,z+1/2,-x+1/2;-z+1/2,x+1/2,y+1/2;y+1/2,-z+1/2,x+1/2;z+1/2,x+1/2,-y+1/2;"
	symLines[200] += "z+1/2,-x+1/2,y+1/2;-y+1/2,z+1/2,x+1/2;-x,-y,z;x,-y,-z;-x,y,-z;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2;"
	symLines[200] += "x+1/2,-y+1/2,z+1/2"
	symLines[201]  = "x,y,z;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;-x,-y,z;x,-y,-z;-x,y,-z;-x,-y,-z;-z,-x,-y;"
	symLines[201] += "-y,-z,-x;y,z,-x;-z,x,y;y,-z,x;z,x,-y;z,-x,y;-y,z,x;x,y,-z;-x,y,z;x,-y,z;x,y+1/2,z+1/2;z,x+1/2,y+1/2;"
	symLines[201] += "y,z+1/2,x+1/2;-y,-z+1/2,x+1/2;z,-x+1/2,-y+1/2;-y,z+1/2,-x+1/2;-z,-x+1/2,y+1/2;-z,x+1/2,-y+1/2;y,-z+1/2,-x+1/2;"
	symLines[201] += "-x,-y+1/2,z+1/2;x,-y+1/2,-z+1/2;-x,y+1/2,-z+1/2;-x,-y+1/2,-z+1/2;-z,-x+1/2,-y+1/2;-y,-z+1/2,-x+1/2;"
	symLines[201] += "y,z+1/2,-x+1/2;-z,x+1/2,y+1/2;y,-z+1/2,x+1/2;z,x+1/2,-y+1/2;z,-x+1/2,y+1/2;-y,z+1/2,x+1/2;x,y+1/2,-z+1/2;"
	symLines[201] += "-x,y+1/2,z+1/2;x,-y+1/2,z+1/2;x+1/2,y,z+1/2;z+1/2,x,y+1/2;y+1/2,z,x+1/2;-y+1/2,-z,x+1/2;z+1/2,-x,-y+1/2;"
	symLines[201] += "-y+1/2,z,-x+1/2;-z+1/2,-x,y+1/2;-z+1/2,x,-y+1/2;y+1/2,-z,-x+1/2;-x+1/2,-y,z+1/2;x+1/2,-y,-z+1/2;"
	symLines[201] += "-x+1/2,y,-z+1/2;-x+1/2,-y,-z+1/2;-z+1/2,-x,-y+1/2;-y+1/2,-z,-x+1/2;y+1/2,z,-x+1/2;-z+1/2,x,y+1/2;"
	symLines[201] += "y+1/2,-z,x+1/2;z+1/2,x,-y+1/2;z+1/2,-x,y+1/2;-y+1/2,z,x+1/2;x+1/2,y,-z+1/2;-x+1/2,y,z+1/2;x+1/2,-y,z+1/2;"
	symLines[201] += "x+1/2,y+1/2,z;z+1/2,x+1/2,y;y+1/2,z+1/2,x;-y+1/2,-z+1/2,x;z+1/2,-x+1/2,-y;-y+1/2,z+1/2,-x;-z+1/2,-x+1/2,y;"
	symLines[201] += "-z+1/2,x+1/2,-y;y+1/2,-z+1/2,-x;-x+1/2,-y+1/2,z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;-x+1/2,-y+1/2,-z;"
	symLines[201] += "-z+1/2,-x+1/2,-y;-y+1/2,-z+1/2,-x;y+1/2,z+1/2,-x;-z+1/2,x+1/2,y;y+1/2,-z+1/2,x;z+1/2,x+1/2,-y;z+1/2,-x+1/2,y;"
	symLines[201] += "-y+1/2,z+1/2,x;x+1/2,y+1/2,-z;-x+1/2,y+1/2,z;x+1/2,-y+1/2,z"
	symLines[202]  = "x,y,z;-x+1/4,-y+1/4,-z+1/4;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;-z+1/4,-x+1/4,-y+1/4;"
	symLines[202] += "-y+1/4,-z+1/4,-x+1/4;y+1/4,z+1/4,-x+1/4;-z+1/4,x+1/4,y+1/4;y+1/4,-z+1/4,x+1/4;z+1/4,x+1/4,-y+1/4;"
	symLines[202] += "z+1/4,-x+1/4,y+1/4;-y+1/4,z+1/4,x+1/4;-x,-y,z;x,-y,-z;-x,y,-z;x+1/4,y+1/4,-z+1/4;-x+1/4,y+1/4,z+1/4;"
	symLines[202] += "x+1/4,-y+1/4,z+1/4;x,y+1/2,z+1/2;-x+1/4,-y+3/4,-z+3/4;z,x+1/2,y+1/2;y,z+1/2,x+1/2;-y,-z+1/2,x+1/2;"
	symLines[202] += "z,-x+1/2,-y+1/2;-y,z+1/2,-x+1/2;-z,-x+1/2,y+1/2;-z,x+1/2,-y+1/2;y,-z+1/2,-x+1/2;-z+1/4,-x+3/4,-y+3/4;"
	symLines[202] += "-y+1/4,-z+3/4,-x+3/4;y+1/4,z+3/4,-x+3/4;-z+1/4,x+3/4,y+3/4;y+1/4,-z+3/4,x+3/4;z+1/4,x+3/4,-y+3/4;"
	symLines[202] += "z+1/4,-x+3/4,y+3/4;-y+1/4,z+3/4,x+3/4;-x,-y+1/2,z+1/2;x,-y+1/2,-z+1/2;-x,y+1/2,-z+1/2;x+1/4,y+3/4,-z+3/4;"
	symLines[202] += "-x+1/4,y+3/4,z+3/4;x+1/4,-y+3/4,z+3/4;x+1/2,y,z+1/2;-x+3/4,-y+1/4,-z+3/4;z+1/2,x,y+1/2;y+1/2,z,x+1/2;"
	symLines[202] += "-y+1/2,-z,x+1/2;z+1/2,-x,-y+1/2;-y+1/2,z,-x+1/2;-z+1/2,-x,y+1/2;-z+1/2,x,-y+1/2;y+1/2,-z,-x+1/2;"
	symLines[202] += "-z+3/4,-x+1/4,-y+3/4;-y+3/4,-z+1/4,-x+3/4;y+3/4,z+1/4,-x+3/4;-z+3/4,x+1/4,y+3/4;y+3/4,-z+1/4,x+3/4;"
	symLines[202] += "z+3/4,x+1/4,-y+3/4;z+3/4,-x+1/4,y+3/4;-y+3/4,z+1/4,x+3/4;-x+1/2,-y,z+1/2;x+1/2,-y,-z+1/2;-x+1/2,y,-z+1/2;"
	symLines[202] += "x+3/4,y+1/4,-z+3/4;-x+3/4,y+1/4,z+3/4;x+3/4,-y+1/4,z+3/4;x+1/2,y+1/2,z;-x+3/4,-y+3/4,-z+1/4;z+1/2,x+1/2,y;"
	symLines[202] += "y+1/2,z+1/2,x;-y+1/2,-z+1/2,x;z+1/2,-x+1/2,-y;-y+1/2,z+1/2,-x;-z+1/2,-x+1/2,y;-z+1/2,x+1/2,-y;y+1/2,-z+1/2,-x;"
	symLines[202] += "-z+3/4,-x+3/4,-y+1/4;-y+3/4,-z+3/4,-x+1/4;y+3/4,z+3/4,-x+1/4;-z+3/4,x+3/4,y+1/4;y+3/4,-z+3/4,x+1/4;"
	symLines[202] += "z+3/4,x+3/4,-y+1/4;z+3/4,-x+3/4,y+1/4;-y+3/4,z+3/4,x+1/4;-x+1/2,-y+1/2,z;x+1/2,-y+1/2,-z;-x+1/2,y+1/2,-z;"
	symLines[202] += "x+3/4,y+3/4,-z+1/4;-x+3/4,y+3/4,z+1/4;x+3/4,-y+3/4,z+1/4"
	symLines[203]  = "x,y,z;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;-x,-y,z;x,-y,-z;-x,y,-z;-x,-y,-z;-z,-x,-y;"
	symLines[203] += "-y,-z,-x;y,z,-x;-z,x,y;y,-z,x;z,x,-y;z,-x,y;-y,z,x;x,y,-z;-x,y,z;x,-y,z;x+1/2,y+1/2,z+1/2;z+1/2,x+1/2,y+1/2;"
	symLines[203] += "y+1/2,z+1/2,x+1/2;-y+1/2,-z+1/2,x+1/2;z+1/2,-x+1/2,-y+1/2;-y+1/2,z+1/2,-x+1/2;-z+1/2,-x+1/2,y+1/2;"
	symLines[203] += "-z+1/2,x+1/2,-y+1/2;y+1/2,-z+1/2,-x+1/2;-x+1/2,-y+1/2,z+1/2;x+1/2,-y+1/2,-z+1/2;-x+1/2,y+1/2,-z+1/2;"
	symLines[203] += "-x+1/2,-y+1/2,-z+1/2;-z+1/2,-x+1/2,-y+1/2;-y+1/2,-z+1/2,-x+1/2;y+1/2,z+1/2,-x+1/2;-z+1/2,x+1/2,y+1/2;"
	symLines[203] += "y+1/2,-z+1/2,x+1/2;z+1/2,x+1/2,-y+1/2;z+1/2,-x+1/2,y+1/2;-y+1/2,z+1/2,x+1/2;x+1/2,y+1/2,-z+1/2;"
	symLines[203] += "-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2"
	symLines[204]  = "x,y,z;z,x,y;y,z,x;-y+1/2,-z,x+1/2;z+1/2,-x+1/2,-y;-y,z+1/2,-x+1/2;-z+1/2,-x,y+1/2;-z,x+1/2,-y+1/2;"
	symLines[204] += "y+1/2,-z+1/2,-x;-x+1/2,-y,z+1/2;x+1/2,-y+1/2,-z;-x,y+1/2,-z+1/2;-x,-y,-z;-z,-x,-y;-y,-z,-x;y+1/2,z,-x+1/2;"
	symLines[204] += "-z+1/2,x+1/2,y;y,-z+1/2,x+1/2;z+1/2,x,-y+1/2;z,-x+1/2,y+1/2;-y+1/2,z+1/2,x;x+1/2,y,-z+1/2;-x+1/2,y+1/2,z;"
	symLines[204] += "x,-y+1/2,z+1/2"
	symLines[205]  = "x,y,z;z,x,y;y,z,x;-y,-z+1/2,x;z,-x,-y+1/2;-y+1/2,z,-x;-z,-x+1/2,y;-z+1/2,x,-y;y,-z,-x+1/2;-x,-y+1/2,z;"
	symLines[205] += "x,-y,-z+1/2;-x+1/2,y,-z;-x,-y,-z;-z,-x,-y;-y,-z,-x;y,z+1/2,-x;-z,x,y+1/2;y+1/2,-z,x;z,x+1/2,-y;z+1/2,-x,y;"
	symLines[205] += "-y,z,x+1/2;x,y+1/2,-z;-x,y,z+1/2;x+1/2,-y,z;x+1/2,y+1/2,z+1/2;z+1/2,x+1/2,y+1/2;y+1/2,z+1/2,x+1/2;"
	symLines[205] += "-y+1/2,-z,x+1/2;z+1/2,-x+1/2,-y;-y,z+1/2,-x+1/2;-z+1/2,-x,y+1/2;-z,x+1/2,-y+1/2;y+1/2,-z+1/2,-x;"
	symLines[205] += "-x+1/2,-y,z+1/2;x+1/2,-y+1/2,-z;-x,y+1/2,-z+1/2;-x+1/2,-y+1/2,-z+1/2;-z+1/2,-x+1/2,-y+1/2;-y+1/2,-z+1/2,-x+1/2;"
	symLines[205] += "y+1/2,z,-x+1/2;-z+1/2,x+1/2,y;y,-z+1/2,x+1/2;z+1/2,x,-y+1/2;z,-x+1/2,y+1/2;-y+1/2,z+1/2,x;x+1/2,y,-z+1/2;"
	symLines[205] += "-x+1/2,y+1/2,z;x,-y+1/2,z+1/2"
	symLines[206]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x,-z,y;x,-y,-z;x,z,-y;z,y,-x;-x,y,-z;-z,y,x;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;"
	symLines[206] += "-z,-x,y;-z,x,-y;y,-z,-x;y,x,-z;-y,-x,-z;-x,z,y;-x,-z,-y;z,-y,x;-z,-y,-x"
	symLines[207]  = "x,y,z;-y+1/2,x+1/2,z+1/2;-x,-y,z;y+1/2,-x+1/2,z+1/2;x+1/2,-z+1/2,y+1/2;x,-y,-z;x+1/2,z+1/2,-y+1/2;"
	symLines[207] += "z+1/2,y+1/2,-x+1/2;-x,y,-z;-z+1/2,y+1/2,x+1/2;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;"
	symLines[207] += "y+1/2,x+1/2,-z+1/2;-y+1/2,-x+1/2,-z+1/2;-x+1/2,z+1/2,y+1/2;-x+1/2,-z+1/2,-y+1/2;z+1/2,-y+1/2,x+1/2;"
	symLines[207] += "-z+1/2,-y+1/2,-x+1/2"
	symLines[208]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x,-z,y;x,-y,-z;x,z,-y;z,y,-x;-x,y,-z;-z,y,x;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;"
	symLines[208] += "-z,-x,y;-z,x,-y;y,-z,-x;y,x,-z;-y,-x,-z;-x,z,y;-x,-z,-y;z,-y,x;-z,-y,-x;x,y+1/2,z+1/2;-y,x+1/2,z+1/2;"
	symLines[208] += "-x,-y+1/2,z+1/2;y,-x+1/2,z+1/2;x,-z+1/2,y+1/2;x,-y+1/2,-z+1/2;x,z+1/2,-y+1/2;z,y+1/2,-x+1/2;-x,y+1/2,-z+1/2;"
	symLines[208] += "-z,y+1/2,x+1/2;z,x+1/2,y+1/2;y,z+1/2,x+1/2;-y,-z+1/2,x+1/2;z,-x+1/2,-y+1/2;-y,z+1/2,-x+1/2;-z,-x+1/2,y+1/2;"
	symLines[208] += "-z,x+1/2,-y+1/2;y,-z+1/2,-x+1/2;y,x+1/2,-z+1/2;-y,-x+1/2,-z+1/2;-x,z+1/2,y+1/2;-x,-z+1/2,-y+1/2;z,-y+1/2,x+1/2;"
	symLines[208] += "-z,-y+1/2,-x+1/2;x+1/2,y,z+1/2;-y+1/2,x,z+1/2;-x+1/2,-y,z+1/2;y+1/2,-x,z+1/2;x+1/2,-z,y+1/2;x+1/2,-y,-z+1/2;"
	symLines[208] += "x+1/2,z,-y+1/2;z+1/2,y,-x+1/2;-x+1/2,y,-z+1/2;-z+1/2,y,x+1/2;z+1/2,x,y+1/2;y+1/2,z,x+1/2;-y+1/2,-z,x+1/2;"
	symLines[208] += "z+1/2,-x,-y+1/2;-y+1/2,z,-x+1/2;-z+1/2,-x,y+1/2;-z+1/2,x,-y+1/2;y+1/2,-z,-x+1/2;y+1/2,x,-z+1/2;"
	symLines[208] += "-y+1/2,-x,-z+1/2;-x+1/2,z,y+1/2;-x+1/2,-z,-y+1/2;z+1/2,-y,x+1/2;-z+1/2,-y,-x+1/2;x+1/2,y+1/2,z;-y+1/2,x+1/2,z;"
	symLines[208] += "-x+1/2,-y+1/2,z;y+1/2,-x+1/2,z;x+1/2,-z+1/2,y;x+1/2,-y+1/2,-z;x+1/2,z+1/2,-y;z+1/2,y+1/2,-x;-x+1/2,y+1/2,-z;"
	symLines[208] += "-z+1/2,y+1/2,x;z+1/2,x+1/2,y;y+1/2,z+1/2,x;-y+1/2,-z+1/2,x;z+1/2,-x+1/2,-y;-y+1/2,z+1/2,-x;-z+1/2,-x+1/2,y;"
	symLines[208] += "-z+1/2,x+1/2,-y;y+1/2,-z+1/2,-x;y+1/2,x+1/2,-z;-y+1/2,-x+1/2,-z;-x+1/2,z+1/2,y;-x+1/2,-z+1/2,-y;z+1/2,-y+1/2,x;"
	symLines[208] += "-z+1/2,-y+1/2,-x"
	symLines[209]  = "x,y,z;-y+1/4,x+1/4,z+1/4;-x,-y,z;y+1/4,-x+1/4,z+1/4;x+1/4,-z+1/4,y+1/4;x,-y,-z;x+1/4,z+1/4,-y+1/4;"
	symLines[209] += "z+1/4,y+1/4,-x+1/4;-x,y,-z;-z+1/4,y+1/4,x+1/4;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;"
	symLines[209] += "y+1/4,x+1/4,-z+1/4;-y+1/4,-x+1/4,-z+1/4;-x+1/4,z+1/4,y+1/4;-x+1/4,-z+1/4,-y+1/4;z+1/4,-y+1/4,x+1/4;"
	symLines[209] += "-z+1/4,-y+1/4,-x+1/4;x,y+1/2,z+1/2;-y+1/4,x+3/4,z+3/4;-x,-y+1/2,z+1/2;y+1/4,-x+3/4,z+3/4;x+1/4,-z+3/4,y+3/4;"
	symLines[209] += "x,-y+1/2,-z+1/2;x+1/4,z+3/4,-y+3/4;z+1/4,y+3/4,-x+3/4;-x,y+1/2,-z+1/2;-z+1/4,y+3/4,x+3/4;z,x+1/2,y+1/2;"
	symLines[209] += "y,z+1/2,x+1/2;-y,-z+1/2,x+1/2;z,-x+1/2,-y+1/2;-y,z+1/2,-x+1/2;-z,-x+1/2,y+1/2;-z,x+1/2,-y+1/2;y,-z+1/2,-x+1/2;"
	symLines[209] += "y+1/4,x+3/4,-z+3/4;-y+1/4,-x+3/4,-z+3/4;-x+1/4,z+3/4,y+3/4;-x+1/4,-z+3/4,-y+3/4;z+1/4,-y+3/4,x+3/4;"
	symLines[209] += "-z+1/4,-y+3/4,-x+3/4;x+1/2,y,z+1/2;-y+3/4,x+1/4,z+3/4;-x+1/2,-y,z+1/2;y+3/4,-x+1/4,z+3/4;x+3/4,-z+1/4,y+3/4;"
	symLines[209] += "x+1/2,-y,-z+1/2;x+3/4,z+1/4,-y+3/4;z+3/4,y+1/4,-x+3/4;-x+1/2,y,-z+1/2;-z+3/4,y+1/4,x+3/4;z+1/2,x,y+1/2;"
	symLines[209] += "y+1/2,z,x+1/2;-y+1/2,-z,x+1/2;z+1/2,-x,-y+1/2;-y+1/2,z,-x+1/2;-z+1/2,-x,y+1/2;-z+1/2,x,-y+1/2;y+1/2,-z,-x+1/2;"
	symLines[209] += "y+3/4,x+1/4,-z+3/4;-y+3/4,-x+1/4,-z+3/4;-x+3/4,z+1/4,y+3/4;-x+3/4,-z+1/4,-y+3/4;z+3/4,-y+1/4,x+3/4;"
	symLines[209] += "-z+3/4,-y+1/4,-x+3/4;x+1/2,y+1/2,z;-y+3/4,x+3/4,z+1/4;-x+1/2,-y+1/2,z;y+3/4,-x+3/4,z+1/4;x+3/4,-z+3/4,y+1/4;"
	symLines[209] += "x+1/2,-y+1/2,-z;x+3/4,z+3/4,-y+1/4;z+3/4,y+3/4,-x+1/4;-x+1/2,y+1/2,-z;-z+3/4,y+3/4,x+1/4;z+1/2,x+1/2,y;"
	symLines[209] += "y+1/2,z+1/2,x;-y+1/2,-z+1/2,x;z+1/2,-x+1/2,-y;-y+1/2,z+1/2,-x;-z+1/2,-x+1/2,y;-z+1/2,x+1/2,-y;y+1/2,-z+1/2,-x;"
	symLines[209] += "y+3/4,x+3/4,-z+1/4;-y+3/4,-x+3/4,-z+1/4;-x+3/4,z+3/4,y+1/4;-x+3/4,-z+3/4,-y+1/4;z+3/4,-y+3/4,x+1/4;"
	symLines[209] += "-z+3/4,-y+3/4,-x+1/4"
	symLines[210]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x,-z,y;x,-y,-z;x,z,-y;z,y,-x;-x,y,-z;-z,y,x;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;"
	symLines[210] += "-z,-x,y;-z,x,-y;y,-z,-x;y,x,-z;-y,-x,-z;-x,z,y;-x,-z,-y;z,-y,x;-z,-y,-x;x+1/2,y+1/2,z+1/2;-y+1/2,x+1/2,z+1/2;"
	symLines[210] += "-x+1/2,-y+1/2,z+1/2;y+1/2,-x+1/2,z+1/2;x+1/2,-z+1/2,y+1/2;x+1/2,-y+1/2,-z+1/2;x+1/2,z+1/2,-y+1/2;"
	symLines[210] += "z+1/2,y+1/2,-x+1/2;-x+1/2,y+1/2,-z+1/2;-z+1/2,y+1/2,x+1/2;z+1/2,x+1/2,y+1/2;y+1/2,z+1/2,x+1/2;"
	symLines[210] += "-y+1/2,-z+1/2,x+1/2;z+1/2,-x+1/2,-y+1/2;-y+1/2,z+1/2,-x+1/2;-z+1/2,-x+1/2,y+1/2;-z+1/2,x+1/2,-y+1/2;"
	symLines[210] += "y+1/2,-z+1/2,-x+1/2;y+1/2,x+1/2,-z+1/2;-y+1/2,-x+1/2,-z+1/2;-x+1/2,z+1/2,y+1/2;-x+1/2,-z+1/2,-y+1/2;"
	symLines[210] += "z+1/2,-y+1/2,x+1/2;-z+1/2,-y+1/2,-x+1/2"
	symLines[211]  = "x,y,z;-y+3/4,x+1/4,z+3/4;-x+1/2,-y,z+1/2;y+3/4,-x+3/4,z+1/4;x+3/4,-z+3/4,y+1/4;x+1/2,-y+1/2,-z;"
	symLines[211] += "x+1/4,z+3/4,-y+3/4;z+1/4,y+3/4,-x+3/4;-x,y+1/2,-z+1/2;-z+3/4,y+1/4,x+3/4;z,x,y;y,z,x;-y+1/2,-z,x+1/2;"
	symLines[211] += "z+1/2,-x+1/2,-y;-y,z+1/2,-x+1/2;-z+1/2,-x,y+1/2;-z,x+1/2,-y+1/2;y+1/2,-z+1/2,-x;y+1/4,x+3/4,-z+3/4;"
	symLines[211] += "-y+1/4,-x+1/4,-z+1/4;-x+3/4,z+1/4,y+3/4;-x+1/4,-z+1/4,-y+1/4;z+3/4,-y+3/4,x+1/4;-z+1/4,-y+1/4,-x+1/4"
	symLines[212]  = "x,y,z;-y+1/4,x+3/4,z+1/4;-x+1/2,-y,z+1/2;y+1/4,-x+1/4,z+3/4;x+1/4,-z+1/4,y+3/4;x+1/2,-y+1/2,-z;"
	symLines[212] += "x+3/4,z+1/4,-y+1/4;z+3/4,y+1/4,-x+1/4;-x,y+1/2,-z+1/2;-z+1/4,y+3/4,x+1/4;z,x,y;y,z,x;-y+1/2,-z,x+1/2;"
	symLines[212] += "z+1/2,-x+1/2,-y;-y,z+1/2,-x+1/2;-z+1/2,-x,y+1/2;-z,x+1/2,-y+1/2;y+1/2,-z+1/2,-x;y+3/4,x+1/4,-z+1/4;"
	symLines[212] += "-y+3/4,-x+3/4,-z+3/4;-x+1/4,z+3/4,y+1/4;-x+3/4,-z+3/4,-y+3/4;z+1/4,-y+1/4,x+3/4;-z+3/4,-y+3/4,-x+3/4"
	symLines[213]  = "x,y,z;-y+1/4,x+3/4,z+1/4;-x,-y+1/2,z;y+1/4,-x+1/4,z+3/4;x+1/4,-z+1/4,y+3/4;x,-y,-z+1/2;x+3/4,z+1/4,-y+1/4;"
	symLines[213] += "z+3/4,y+1/4,-x+1/4;-x+1/2,y,-z;-z+1/4,y+3/4,x+1/4;z,x,y;y,z,x;-y,-z+1/2,x;z,-x,-y+1/2;-y+1/2,z,-x;-z,-x+1/2,y;"
	symLines[213] += "-z+1/2,x,-y;y,-z,-x+1/2;y+3/4,x+1/4,-z+1/4;-y+1/4,-x+1/4,-z+1/4;-x+1/4,z+3/4,y+1/4;-x+1/4,-z+1/4,-y+1/4;"
	symLines[213] += "z+1/4,-y+1/4,x+3/4;-z+1/4,-y+1/4,-x+1/4;x+1/2,y+1/2,z+1/2;-y+3/4,x+1/4,z+3/4;-x+1/2,-y,z+1/2;"
	symLines[213] += "y+3/4,-x+3/4,z+1/4;x+3/4,-z+3/4,y+1/4;x+1/2,-y+1/2,-z;x+1/4,z+3/4,-y+3/4;z+1/4,y+3/4,-x+3/4;-x,y+1/2,-z+1/2;"
	symLines[213] += "-z+3/4,y+1/4,x+3/4;z+1/2,x+1/2,y+1/2;y+1/2,z+1/2,x+1/2;-y+1/2,-z,x+1/2;z+1/2,-x+1/2,-y;-y,z+1/2,-x+1/2;"
	symLines[213] += "-z+1/2,-x,y+1/2;-z,x+1/2,-y+1/2;y+1/2,-z+1/2,-x;y+1/4,x+3/4,-z+3/4;-y+3/4,-x+3/4,-z+3/4;-x+3/4,z+1/4,y+3/4;"
	symLines[213] += "-x+3/4,-z+3/4,-y+3/4;z+3/4,-y+3/4,x+1/4;-z+3/4,-y+3/4,-x+3/4"
	symLines[214]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;-x,z,-y;x,-y,-z;-x,-z,y;-z,-y,x;-x,y,-z;z,-y,-x;z,x,y;y,z,x;-y,-z,x;z,-x,-y;"
	symLines[214] += "-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;-y,-x,z;y,x,z;x,-z,-y;x,z,y;-z,y,-x;z,y,x"
	symLines[215]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;-x,z,-y;x,-y,-z;-x,-z,y;-z,-y,x;-x,y,-z;z,-y,-x;z,x,y;y,z,x;-y,-z,x;z,-x,-y;"
	symLines[215] += "-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;-y,-x,z;y,x,z;x,-z,-y;x,z,y;-z,y,-x;z,y,x;x,y+1/2,z+1/2;y,-x+1/2,-z+1/2;"
	symLines[215] += "-x,-y+1/2,z+1/2;-y,x+1/2,-z+1/2;-x,z+1/2,-y+1/2;x,-y+1/2,-z+1/2;-x,-z+1/2,y+1/2;-z,-y+1/2,x+1/2;"
	symLines[215] += "-x,y+1/2,-z+1/2;z,-y+1/2,-x+1/2;z,x+1/2,y+1/2;y,z+1/2,x+1/2;-y,-z+1/2,x+1/2;z,-x+1/2,-y+1/2;-y,z+1/2,-x+1/2;"
	symLines[215] += "-z,-x+1/2,y+1/2;-z,x+1/2,-y+1/2;y,-z+1/2,-x+1/2;-y,-x+1/2,z+1/2;y,x+1/2,z+1/2;x,-z+1/2,-y+1/2;x,z+1/2,y+1/2;"
	symLines[215] += "-z,y+1/2,-x+1/2;z,y+1/2,x+1/2;x+1/2,y,z+1/2;y+1/2,-x,-z+1/2;-x+1/2,-y,z+1/2;-y+1/2,x,-z+1/2;-x+1/2,z,-y+1/2;"
	symLines[215] += "x+1/2,-y,-z+1/2;-x+1/2,-z,y+1/2;-z+1/2,-y,x+1/2;-x+1/2,y,-z+1/2;z+1/2,-y,-x+1/2;z+1/2,x,y+1/2;y+1/2,z,x+1/2;"
	symLines[215] += "-y+1/2,-z,x+1/2;z+1/2,-x,-y+1/2;-y+1/2,z,-x+1/2;-z+1/2,-x,y+1/2;-z+1/2,x,-y+1/2;y+1/2,-z,-x+1/2;"
	symLines[215] += "-y+1/2,-x,z+1/2;y+1/2,x,z+1/2;x+1/2,-z,-y+1/2;x+1/2,z,y+1/2;-z+1/2,y,-x+1/2;z+1/2,y,x+1/2;x+1/2,y+1/2,z;"
	symLines[215] += "y+1/2,-x+1/2,-z;-x+1/2,-y+1/2,z;-y+1/2,x+1/2,-z;-x+1/2,z+1/2,-y;x+1/2,-y+1/2,-z;-x+1/2,-z+1/2,y;"
	symLines[215] += "-z+1/2,-y+1/2,x;-x+1/2,y+1/2,-z;z+1/2,-y+1/2,-x;z+1/2,x+1/2,y;y+1/2,z+1/2,x;-y+1/2,-z+1/2,x;z+1/2,-x+1/2,-y;"
	symLines[215] += "-y+1/2,z+1/2,-x;-z+1/2,-x+1/2,y;-z+1/2,x+1/2,-y;y+1/2,-z+1/2,-x;-y+1/2,-x+1/2,z;y+1/2,x+1/2,z;x+1/2,-z+1/2,-y;"
	symLines[215] += "x+1/2,z+1/2,y;-z+1/2,y+1/2,-x;z+1/2,y+1/2,x"
	symLines[216]  = "x,y,z;y,-x,-z;-x,-y,z;-y,x,-z;-x,z,-y;x,-y,-z;-x,-z,y;-z,-y,x;-x,y,-z;z,-y,-x;z,x,y;y,z,x;-y,-z,x;z,-x,-y;"
	symLines[216] += "-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;-y,-x,z;y,x,z;x,-z,-y;x,z,y;-z,y,-x;z,y,x;x+1/2,y+1/2,z+1/2;"
	symLines[216] += "y+1/2,-x+1/2,-z+1/2;-x+1/2,-y+1/2,z+1/2;-y+1/2,x+1/2,-z+1/2;-x+1/2,z+1/2,-y+1/2;x+1/2,-y+1/2,-z+1/2;"
	symLines[216] += "-x+1/2,-z+1/2,y+1/2;-z+1/2,-y+1/2,x+1/2;-x+1/2,y+1/2,-z+1/2;z+1/2,-y+1/2,-x+1/2;z+1/2,x+1/2,y+1/2;"
	symLines[216] += "y+1/2,z+1/2,x+1/2;-y+1/2,-z+1/2,x+1/2;z+1/2,-x+1/2,-y+1/2;-y+1/2,z+1/2,-x+1/2;-z+1/2,-x+1/2,y+1/2;"
	symLines[216] += "-z+1/2,x+1/2,-y+1/2;y+1/2,-z+1/2,-x+1/2;-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2;x+1/2,-z+1/2,-y+1/2;"
	symLines[216] += "x+1/2,z+1/2,y+1/2;-z+1/2,y+1/2,-x+1/2;z+1/2,y+1/2,x+1/2"
	symLines[217]  = "x,y,z;y+1/2,-x+1/2,-z+1/2;-x,-y,z;-y+1/2,x+1/2,-z+1/2;-x+1/2,z+1/2,-y+1/2;x,-y,-z;-x+1/2,-z+1/2,y+1/2;"
	symLines[217] += "-z+1/2,-y+1/2,x+1/2;-x,y,-z;z+1/2,-y+1/2,-x+1/2;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;"
	symLines[217] += "-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2;x+1/2,-z+1/2,-y+1/2;x+1/2,z+1/2,y+1/2;-z+1/2,y+1/2,-x+1/2;"
	symLines[217] += "z+1/2,y+1/2,x+1/2"
	symLines[218]  = "x,y,z;y,-x,-z+1/2;-x,-y,z;-y,x,-z+1/2;-x,z,-y+1/2;x,-y,-z;-x,-z,y+1/2;-z,-y,x+1/2;-x,y,-z;z,-y,-x+1/2;z,x,y;"
	symLines[218] += "y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;-y,-x,z+1/2;y,x,z+1/2;x,-z,-y+1/2;x,z,y+1/2;-z,y,-x+1/2;"
	symLines[218] += "z,y,x+1/2;x,y+1/2,z+1/2;y,-x+1/2,-z;-x,-y+1/2,z+1/2;-y,x+1/2,-z;-x,z+1/2,-y;x,-y+1/2,-z+1/2;-x,-z+1/2,y;"
	symLines[218] += "-z,-y+1/2,x;-x,y+1/2,-z+1/2;z,-y+1/2,-x;z,x+1/2,y+1/2;y,z+1/2,x+1/2;-y,-z+1/2,x+1/2;z,-x+1/2,-y+1/2;"
	symLines[218] += "-y,z+1/2,-x+1/2;-z,-x+1/2,y+1/2;-z,x+1/2,-y+1/2;y,-z+1/2,-x+1/2;-y,-x+1/2,z;y,x+1/2,z;x,-z+1/2,-y;x,z+1/2,y;"
	symLines[218] += "-z,y+1/2,-x;z,y+1/2,x;x+1/2,y,z+1/2;y+1/2,-x,-z;-x+1/2,-y,z+1/2;-y+1/2,x,-z;-x+1/2,z,-y;x+1/2,-y,-z+1/2;"
	symLines[218] += "-x+1/2,-z,y;-z+1/2,-y,x;-x+1/2,y,-z+1/2;z+1/2,-y,-x;z+1/2,x,y+1/2;y+1/2,z,x+1/2;-y+1/2,-z,x+1/2;"
	symLines[218] += "z+1/2,-x,-y+1/2;-y+1/2,z,-x+1/2;-z+1/2,-x,y+1/2;-z+1/2,x,-y+1/2;y+1/2,-z,-x+1/2;-y+1/2,-x,z;y+1/2,x,z;"
	symLines[218] += "x+1/2,-z,-y;x+1/2,z,y;-z+1/2,y,-x;z+1/2,y,x;x+1/2,y+1/2,z;y+1/2,-x+1/2,-z+1/2;-x+1/2,-y+1/2,z;"
	symLines[218] += "-y+1/2,x+1/2,-z+1/2;-x+1/2,z+1/2,-y+1/2;x+1/2,-y+1/2,-z;-x+1/2,-z+1/2,y+1/2;-z+1/2,-y+1/2,x+1/2;"
	symLines[218] += "-x+1/2,y+1/2,-z;z+1/2,-y+1/2,-x+1/2;z+1/2,x+1/2,y;y+1/2,z+1/2,x;-y+1/2,-z+1/2,x;z+1/2,-x+1/2,-y;"
	symLines[218] += "-y+1/2,z+1/2,-x;-z+1/2,-x+1/2,y;-z+1/2,x+1/2,-y;y+1/2,-z+1/2,-x;-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2;"
	symLines[218] += "x+1/2,-z+1/2,-y+1/2;x+1/2,z+1/2,y+1/2;-z+1/2,y+1/2,-x+1/2;z+1/2,y+1/2,x+1/2"
	symLines[219]  = "x,y,z;y+1/4,-x+3/4,-z+1/4;-x,-y+1/2,z;-y+1/4,x+1/4,-z+3/4;-x+1/4,z+1/4,-y+3/4;x,-y,-z+1/2;-x+3/4,-z+1/4,y+1/4;"
	symLines[219] += "-z+3/4,-y+1/4,x+1/4;-x+1/2,y,-z;z+1/4,-y+3/4,-x+1/4;z,x,y;y,z,x;-y,-z+1/2,x;z,-x,-y+1/2;-y+1/2,z,-x;"
	symLines[219] += "-z,-x+1/2,y;-z+1/2,x,-y;y,-z,-x+1/2;-y+3/4,-x+1/4,z+1/4;y+1/4,x+1/4,z+1/4;x+1/4,-z+3/4,-y+1/4;"
	symLines[219] += "x+1/4,z+1/4,y+1/4;-z+1/4,y+1/4,-x+3/4;z+1/4,y+1/4,x+1/4;x+1/2,y+1/2,z+1/2;y+3/4,-x+1/4,-z+3/4;-x+1/2,-y,z+1/2;"
	symLines[219] += "-y+3/4,x+3/4,-z+1/4;-x+3/4,z+3/4,-y+1/4;x+1/2,-y+1/2,-z;-x+1/4,-z+3/4,y+3/4;-z+1/4,-y+3/4,x+3/4;"
	symLines[219] += "-x,y+1/2,-z+1/2;z+3/4,-y+1/4,-x+3/4;z+1/2,x+1/2,y+1/2;y+1/2,z+1/2,x+1/2;-y+1/2,-z,x+1/2;z+1/2,-x+1/2,-y;"
	symLines[219] += "-y,z+1/2,-x+1/2;-z+1/2,-x,y+1/2;-z,x+1/2,-y+1/2;y+1/2,-z+1/2,-x;-y+1/4,-x+3/4,z+3/4;y+3/4,x+3/4,z+3/4;"
	symLines[219] += "x+3/4,-z+1/4,-y+3/4;x+3/4,z+3/4,y+3/4;-z+3/4,y+3/4,-x+1/4;z+3/4,y+3/4,x+3/4"
	symLines[220]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x,-z,y;x,-y,-z;x,z,-y;z,y,-x;-x,y,-z;-z,y,x;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;"
	symLines[220] += "-z,-x,y;-z,x,-y;y,-z,-x;y,x,-z;-y,-x,-z;-x,z,y;-x,-z,-y;z,-y,x;-z,-y,-x;-x,-y,-z;y,-x,-z;x,y,-z;-y,x,-z;"
	symLines[220] += "-x,z,-y;-x,y,z;-x,-z,y;-z,-y,x;x,-y,z;z,-y,-x;-z,-x,-y;-y,-z,-x;y,z,-x;-z,x,y;y,-z,x;z,x,-y;z,-x,y;-y,z,x;"
	symLines[220] += "-y,-x,z;y,x,z;x,-z,-y;x,z,y;-z,y,-x;z,y,x"
	symLines[221]  = "x,y,z;-x+1/2,-y+1/2,-z+1/2;-y,x,z;-x,-y,z;y,-x,z;x,-z,y;x,-y,-z;x,z,-y;z,y,-x;-x,y,-z;-z,y,x;"
	symLines[221] += "y+1/2,-x+1/2,-z+1/2;-y+1/2,x+1/2,-z+1/2;-x+1/2,z+1/2,-y+1/2;-x+1/2,-z+1/2,y+1/2;-z+1/2,-y+1/2,x+1/2;"
	symLines[221] += "z+1/2,-y+1/2,-x+1/2;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;-z+1/2,-x+1/2,-y+1/2;"
	symLines[221] += "-y+1/2,-z+1/2,-x+1/2;y+1/2,z+1/2,-x+1/2;-z+1/2,x+1/2,y+1/2;y+1/2,-z+1/2,x+1/2;z+1/2,x+1/2,-y+1/2;"
	symLines[221] += "z+1/2,-x+1/2,y+1/2;-y+1/2,z+1/2,x+1/2;y,x,-z;-y,-x,-z;-x,z,y;-x,-z,-y;z,-y,x;-z,-y,-x;x+1/2,y+1/2,-z+1/2;"
	symLines[221] += "-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2;-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2;x+1/2,-z+1/2,-y+1/2;"
	symLines[221] += "x+1/2,z+1/2,y+1/2;-z+1/2,y+1/2,-x+1/2;z+1/2,y+1/2,x+1/2"
	symLines[222]  = "x,y,z;-y+1/2,x+1/2,z+1/2;-x,-y,z;y+1/2,-x+1/2,z+1/2;x+1/2,-z+1/2,y+1/2;x,-y,-z;x+1/2,z+1/2,-y+1/2;"
	symLines[222] += "z+1/2,y+1/2,-x+1/2;-x,y,-z;-z+1/2,y+1/2,x+1/2;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;"
	symLines[222] += "y+1/2,x+1/2,-z+1/2;-y+1/2,-x+1/2,-z+1/2;-x+1/2,z+1/2,y+1/2;-x+1/2,-z+1/2,-y+1/2;z+1/2,-y+1/2,x+1/2;"
	symLines[222] += "-z+1/2,-y+1/2,-x+1/2;-x,-y,-z;y+1/2,-x+1/2,-z+1/2;x,y,-z;-y+1/2,x+1/2,-z+1/2;-x+1/2,z+1/2,-y+1/2;-x,y,z;"
	symLines[222] += "-x+1/2,-z+1/2,y+1/2;-z+1/2,-y+1/2,x+1/2;x,-y,z;z+1/2,-y+1/2,-x+1/2;-z,-x,-y;-y,-z,-x;y,z,-x;-z,x,y;y,-z,x;"
	symLines[222] += "z,x,-y;z,-x,y;-y,z,x;-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2;x+1/2,-z+1/2,-y+1/2;x+1/2,z+1/2,y+1/2;"
	symLines[222] += "-z+1/2,y+1/2,-x+1/2;z+1/2,y+1/2,x+1/2"
	symLines[223]  = "x,y,z;-x+1/2,-y+1/2,-z+1/2;-y+1/2,x+1/2,z+1/2;-x,-y,z;y+1/2,-x+1/2,z+1/2;x+1/2,-z+1/2,y+1/2;x,-y,-z;"
	symLines[223] += "x+1/2,z+1/2,-y+1/2;z+1/2,y+1/2,-x+1/2;-x,y,-z;-z+1/2,y+1/2,x+1/2;y,-x,-z;-y,x,-z;-x,z,-y;-x,-z,y;-z,-y,x;"
	symLines[223] += "z,-y,-x;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;-z+1/2,-x+1/2,-y+1/2;-y+1/2,-z+1/2,-x+1/2;"
	symLines[223] += "y+1/2,z+1/2,-x+1/2;-z+1/2,x+1/2,y+1/2;y+1/2,-z+1/2,x+1/2;z+1/2,x+1/2,-y+1/2;z+1/2,-x+1/2,y+1/2;"
	symLines[223] += "-y+1/2,z+1/2,x+1/2;y+1/2,x+1/2,-z+1/2;-y+1/2,-x+1/2,-z+1/2;-x+1/2,z+1/2,y+1/2;-x+1/2,-z+1/2,-y+1/2;"
	symLines[223] += "z+1/2,-y+1/2,x+1/2;-z+1/2,-y+1/2,-x+1/2;x+1/2,y+1/2,-z+1/2;-x+1/2,y+1/2,z+1/2;x+1/2,-y+1/2,z+1/2;-y,-x,z;y,x,z;"
	symLines[223] += "x,-z,-y;x,z,y;-z,y,-x;z,y,x"
	symLines[224]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x,-z,y;x,-y,-z;x,z,-y;z,y,-x;-x,y,-z;-z,y,x;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;"
	symLines[224] += "-z,-x,y;-z,x,-y;y,-z,-x;y,x,-z;-y,-x,-z;-x,z,y;-x,-z,-y;z,-y,x;-z,-y,-x;-x,-y,-z;y,-x,-z;x,y,-z;-y,x,-z;"
	symLines[224] += "-x,z,-y;-x,y,z;-x,-z,y;-z,-y,x;x,-y,z;z,-y,-x;-z,-x,-y;-y,-z,-x;y,z,-x;-z,x,y;y,-z,x;z,x,-y;z,-x,y;-y,z,x;"
	symLines[224] += "-y,-x,z;y,x,z;x,-z,-y;x,z,y;-z,y,-x;z,y,x;x,y+1/2,z+1/2;-y,x+1/2,z+1/2;-x,-y+1/2,z+1/2;y,-x+1/2,z+1/2;"
	symLines[224] += "x,-z+1/2,y+1/2;x,-y+1/2,-z+1/2;x,z+1/2,-y+1/2;z,y+1/2,-x+1/2;-x,y+1/2,-z+1/2;-z,y+1/2,x+1/2;z,x+1/2,y+1/2;"
	symLines[224] += "y,z+1/2,x+1/2;-y,-z+1/2,x+1/2;z,-x+1/2,-y+1/2;-y,z+1/2,-x+1/2;-z,-x+1/2,y+1/2;-z,x+1/2,-y+1/2;y,-z+1/2,-x+1/2;"
	symLines[224] += "y,x+1/2,-z+1/2;-y,-x+1/2,-z+1/2;-x,z+1/2,y+1/2;-x,-z+1/2,-y+1/2;z,-y+1/2,x+1/2;-z,-y+1/2,-x+1/2;"
	symLines[224] += "-x,-y+1/2,-z+1/2;y,-x+1/2,-z+1/2;x,y+1/2,-z+1/2;-y,x+1/2,-z+1/2;-x,z+1/2,-y+1/2;-x,y+1/2,z+1/2;-x,-z+1/2,y+1/2;"
	symLines[224] += "-z,-y+1/2,x+1/2;x,-y+1/2,z+1/2;z,-y+1/2,-x+1/2;-z,-x+1/2,-y+1/2;-y,-z+1/2,-x+1/2;y,z+1/2,-x+1/2;-z,x+1/2,y+1/2;"
	symLines[224] += "y,-z+1/2,x+1/2;z,x+1/2,-y+1/2;z,-x+1/2,y+1/2;-y,z+1/2,x+1/2;-y,-x+1/2,z+1/2;y,x+1/2,z+1/2;x,-z+1/2,-y+1/2;"
	symLines[224] += "x,z+1/2,y+1/2;-z,y+1/2,-x+1/2;z,y+1/2,x+1/2;x+1/2,y,z+1/2;-y+1/2,x,z+1/2;-x+1/2,-y,z+1/2;y+1/2,-x,z+1/2;"
	symLines[224] += "x+1/2,-z,y+1/2;x+1/2,-y,-z+1/2;x+1/2,z,-y+1/2;z+1/2,y,-x+1/2;-x+1/2,y,-z+1/2;-z+1/2,y,x+1/2;z+1/2,x,y+1/2;"
	symLines[224] += "y+1/2,z,x+1/2;-y+1/2,-z,x+1/2;z+1/2,-x,-y+1/2;-y+1/2,z,-x+1/2;-z+1/2,-x,y+1/2;-z+1/2,x,-y+1/2;y+1/2,-z,-x+1/2;"
	symLines[224] += "y+1/2,x,-z+1/2;-y+1/2,-x,-z+1/2;-x+1/2,z,y+1/2;-x+1/2,-z,-y+1/2;z+1/2,-y,x+1/2;-z+1/2,-y,-x+1/2;"
	symLines[224] += "-x+1/2,-y,-z+1/2;y+1/2,-x,-z+1/2;x+1/2,y,-z+1/2;-y+1/2,x,-z+1/2;-x+1/2,z,-y+1/2;-x+1/2,y,z+1/2;-x+1/2,-z,y+1/2;"
	symLines[224] += "-z+1/2,-y,x+1/2;x+1/2,-y,z+1/2;z+1/2,-y,-x+1/2;-z+1/2,-x,-y+1/2;-y+1/2,-z,-x+1/2;y+1/2,z,-x+1/2;-z+1/2,x,y+1/2;"
	symLines[224] += "y+1/2,-z,x+1/2;z+1/2,x,-y+1/2;z+1/2,-x,y+1/2;-y+1/2,z,x+1/2;-y+1/2,-x,z+1/2;y+1/2,x,z+1/2;x+1/2,-z,-y+1/2;"
	symLines[224] += "x+1/2,z,y+1/2;-z+1/2,y,-x+1/2;z+1/2,y,x+1/2;x+1/2,y+1/2,z;-y+1/2,x+1/2,z;-x+1/2,-y+1/2,z;y+1/2,-x+1/2,z;"
	symLines[224] += "x+1/2,-z+1/2,y;x+1/2,-y+1/2,-z;x+1/2,z+1/2,-y;z+1/2,y+1/2,-x;-x+1/2,y+1/2,-z;-z+1/2,y+1/2,x;z+1/2,x+1/2,y;"
	symLines[224] += "y+1/2,z+1/2,x;-y+1/2,-z+1/2,x;z+1/2,-x+1/2,-y;-y+1/2,z+1/2,-x;-z+1/2,-x+1/2,y;-z+1/2,x+1/2,-y;y+1/2,-z+1/2,-x;"
	symLines[224] += "y+1/2,x+1/2,-z;-y+1/2,-x+1/2,-z;-x+1/2,z+1/2,y;-x+1/2,-z+1/2,-y;z+1/2,-y+1/2,x;-z+1/2,-y+1/2,-x;"
	symLines[224] += "-x+1/2,-y+1/2,-z;y+1/2,-x+1/2,-z;x+1/2,y+1/2,-z;-y+1/2,x+1/2,-z;-x+1/2,z+1/2,-y;-x+1/2,y+1/2,z;-x+1/2,-z+1/2,y;"
	symLines[224] += "-z+1/2,-y+1/2,x;x+1/2,-y+1/2,z;z+1/2,-y+1/2,-x;-z+1/2,-x+1/2,-y;-y+1/2,-z+1/2,-x;y+1/2,z+1/2,-x;-z+1/2,x+1/2,y;"
	symLines[224] += "y+1/2,-z+1/2,x;z+1/2,x+1/2,-y;z+1/2,-x+1/2,y;-y+1/2,z+1/2,x;-y+1/2,-x+1/2,z;y+1/2,x+1/2,z;x+1/2,-z+1/2,-y;"
	symLines[224] += "x+1/2,z+1/2,y;-z+1/2,y+1/2,-x;z+1/2,y+1/2,x"
	symLines[225]  = "x,y,z;-y,x,z+1/2;-x,-y,z;y,-x,z+1/2;x,-z,y+1/2;x,-y,-z;x,z,-y+1/2;z,y,-x+1/2;-x,y,-z;-z,y,x+1/2;z,x,y;y,z,x;"
	symLines[225] += "-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;y,x,-z+1/2;-y,-x,-z+1/2;-x,z,y+1/2;-x,-z,-y+1/2;z,-y,x+1/2;"
	symLines[225] += "-z,-y,-x+1/2;-x,-y,-z;y,-x,-z+1/2;x,y,-z;-y,x,-z+1/2;-x,z,-y+1/2;-x,y,z;-x,-z,y+1/2;-z,-y,x+1/2;x,-y,z;"
	symLines[225] += "z,-y,-x+1/2;-z,-x,-y;-y,-z,-x;y,z,-x;-z,x,y;y,-z,x;z,x,-y;z,-x,y;-y,z,x;-y,-x,z+1/2;y,x,z+1/2;x,-z,-y+1/2;"
	symLines[225] += "x,z,y+1/2;-z,y,-x+1/2;z,y,x+1/2;x,y+1/2,z+1/2;-y,x+1/2,z;-x,-y+1/2,z+1/2;y,-x+1/2,z;x,-z+1/2,y;x,-y+1/2,-z+1/2;"
	symLines[225] += "x,z+1/2,-y;z,y+1/2,-x;-x,y+1/2,-z+1/2;-z,y+1/2,x;z,x+1/2,y+1/2;y,z+1/2,x+1/2;-y,-z+1/2,x+1/2;z,-x+1/2,-y+1/2;"
	symLines[225] += "-y,z+1/2,-x+1/2;-z,-x+1/2,y+1/2;-z,x+1/2,-y+1/2;y,-z+1/2,-x+1/2;y,x+1/2,-z;-y,-x+1/2,-z;-x,z+1/2,y;"
	symLines[225] += "-x,-z+1/2,-y;z,-y+1/2,x;-z,-y+1/2,-x;-x,-y+1/2,-z+1/2;y,-x+1/2,-z;x,y+1/2,-z+1/2;-y,x+1/2,-z;-x,z+1/2,-y;"
	symLines[225] += "-x,y+1/2,z+1/2;-x,-z+1/2,y;-z,-y+1/2,x;x,-y+1/2,z+1/2;z,-y+1/2,-x;-z,-x+1/2,-y+1/2;-y,-z+1/2,-x+1/2;"
	symLines[225] += "y,z+1/2,-x+1/2;-z,x+1/2,y+1/2;y,-z+1/2,x+1/2;z,x+1/2,-y+1/2;z,-x+1/2,y+1/2;-y,z+1/2,x+1/2;-y,-x+1/2,z;"
	symLines[225] += "y,x+1/2,z;x,-z+1/2,-y;x,z+1/2,y;-z,y+1/2,-x;z,y+1/2,x;x+1/2,y,z+1/2;-y+1/2,x,z;-x+1/2,-y,z+1/2;y+1/2,-x,z;"
	symLines[225] += "x+1/2,-z,y;x+1/2,-y,-z+1/2;x+1/2,z,-y;z+1/2,y,-x;-x+1/2,y,-z+1/2;-z+1/2,y,x;z+1/2,x,y+1/2;y+1/2,z,x+1/2;"
	symLines[225] += "-y+1/2,-z,x+1/2;z+1/2,-x,-y+1/2;-y+1/2,z,-x+1/2;-z+1/2,-x,y+1/2;-z+1/2,x,-y+1/2;y+1/2,-z,-x+1/2;y+1/2,x,-z;"
	symLines[225] += "-y+1/2,-x,-z;-x+1/2,z,y;-x+1/2,-z,-y;z+1/2,-y,x;-z+1/2,-y,-x;-x+1/2,-y,-z+1/2;y+1/2,-x,-z;x+1/2,y,-z+1/2;"
	symLines[225] += "-y+1/2,x,-z;-x+1/2,z,-y;-x+1/2,y,z+1/2;-x+1/2,-z,y;-z+1/2,-y,x;x+1/2,-y,z+1/2;z+1/2,-y,-x;-z+1/2,-x,-y+1/2;"
	symLines[225] += "-y+1/2,-z,-x+1/2;y+1/2,z,-x+1/2;-z+1/2,x,y+1/2;y+1/2,-z,x+1/2;z+1/2,x,-y+1/2;z+1/2,-x,y+1/2;-y+1/2,z,x+1/2;"
	symLines[225] += "-y+1/2,-x,z;y+1/2,x,z;x+1/2,-z,-y;x+1/2,z,y;-z+1/2,y,-x;z+1/2,y,x;x+1/2,y+1/2,z;-y+1/2,x+1/2,z+1/2;"
	symLines[225] += "-x+1/2,-y+1/2,z;y+1/2,-x+1/2,z+1/2;x+1/2,-z+1/2,y+1/2;x+1/2,-y+1/2,-z;x+1/2,z+1/2,-y+1/2;z+1/2,y+1/2,-x+1/2;"
	symLines[225] += "-x+1/2,y+1/2,-z;-z+1/2,y+1/2,x+1/2;z+1/2,x+1/2,y;y+1/2,z+1/2,x;-y+1/2,-z+1/2,x;z+1/2,-x+1/2,-y;-y+1/2,z+1/2,-x;"
	symLines[225] += "-z+1/2,-x+1/2,y;-z+1/2,x+1/2,-y;y+1/2,-z+1/2,-x;y+1/2,x+1/2,-z+1/2;-y+1/2,-x+1/2,-z+1/2;-x+1/2,z+1/2,y+1/2;"
	symLines[225] += "-x+1/2,-z+1/2,-y+1/2;z+1/2,-y+1/2,x+1/2;-z+1/2,-y+1/2,-x+1/2;-x+1/2,-y+1/2,-z;y+1/2,-x+1/2,-z+1/2;"
	symLines[225] += "x+1/2,y+1/2,-z;-y+1/2,x+1/2,-z+1/2;-x+1/2,z+1/2,-y+1/2;-x+1/2,y+1/2,z;-x+1/2,-z+1/2,y+1/2;-z+1/2,-y+1/2,x+1/2;"
	symLines[225] += "x+1/2,-y+1/2,z;z+1/2,-y+1/2,-x+1/2;-z+1/2,-x+1/2,-y;-y+1/2,-z+1/2,-x;y+1/2,z+1/2,-x;-z+1/2,x+1/2,y;"
	symLines[225] += "y+1/2,-z+1/2,x;z+1/2,x+1/2,-y;z+1/2,-x+1/2,y;-y+1/2,z+1/2,x;-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2;"
	symLines[225] += "x+1/2,-z+1/2,-y+1/2;x+1/2,z+1/2,y+1/2;-z+1/2,y+1/2,-x+1/2;z+1/2,y+1/2,x+1/2"
	symLines[226]  = "x,y,z;-x+1/4,-y+1/4,-z+1/4;-y+1/4,x+1/4,z+1/4;-x,-y,z;y+1/4,-x+1/4,z+1/4;x+1/4,-z+1/4,y+1/4;x,-y,-z;"
	symLines[226] += "x+1/4,z+1/4,-y+1/4;z+1/4,y+1/4,-x+1/4;-x,y,-z;-z+1/4,y+1/4,x+1/4;y,-x,-z;-y,x,-z;-x,z,-y;-x,-z,y;-z,-y,x;"
	symLines[226] += "z,-y,-x;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;-z+1/4,-x+1/4,-y+1/4;-y+1/4,-z+1/4,-x+1/4;"
	symLines[226] += "y+1/4,z+1/4,-x+1/4;-z+1/4,x+1/4,y+1/4;y+1/4,-z+1/4,x+1/4;z+1/4,x+1/4,-y+1/4;z+1/4,-x+1/4,y+1/4;"
	symLines[226] += "-y+1/4,z+1/4,x+1/4;y+1/4,x+1/4,-z+1/4;-y+1/4,-x+1/4,-z+1/4;-x+1/4,z+1/4,y+1/4;-x+1/4,-z+1/4,-y+1/4;"
	symLines[226] += "z+1/4,-y+1/4,x+1/4;-z+1/4,-y+1/4,-x+1/4;x+1/4,y+1/4,-z+1/4;-x+1/4,y+1/4,z+1/4;x+1/4,-y+1/4,z+1/4;-y,-x,z;y,x,z;"
	symLines[226] += "x,-z,-y;x,z,y;-z,y,-x;z,y,x;x,y+1/2,z+1/2;-x+1/4,-y+3/4,-z+3/4;-y+1/4,x+3/4,z+3/4;-x,-y+1/2,z+1/2;"
	symLines[226] += "y+1/4,-x+3/4,z+3/4;x+1/4,-z+3/4,y+3/4;x,-y+1/2,-z+1/2;x+1/4,z+3/4,-y+3/4;z+1/4,y+3/4,-x+3/4;-x,y+1/2,-z+1/2;"
	symLines[226] += "-z+1/4,y+3/4,x+3/4;y,-x+1/2,-z+1/2;-y,x+1/2,-z+1/2;-x,z+1/2,-y+1/2;-x,-z+1/2,y+1/2;-z,-y+1/2,x+1/2;"
	symLines[226] += "z,-y+1/2,-x+1/2;z,x+1/2,y+1/2;y,z+1/2,x+1/2;-y,-z+1/2,x+1/2;z,-x+1/2,-y+1/2;-y,z+1/2,-x+1/2;-z,-x+1/2,y+1/2;"
	symLines[226] += "-z,x+1/2,-y+1/2;y,-z+1/2,-x+1/2;-z+1/4,-x+3/4,-y+3/4;-y+1/4,-z+3/4,-x+3/4;y+1/4,z+3/4,-x+3/4;"
	symLines[226] += "-z+1/4,x+3/4,y+3/4;y+1/4,-z+3/4,x+3/4;z+1/4,x+3/4,-y+3/4;z+1/4,-x+3/4,y+3/4;-y+1/4,z+3/4,x+3/4;"
	symLines[226] += "y+1/4,x+3/4,-z+3/4;-y+1/4,-x+3/4,-z+3/4;-x+1/4,z+3/4,y+3/4;-x+1/4,-z+3/4,-y+3/4;z+1/4,-y+3/4,x+3/4;"
	symLines[226] += "-z+1/4,-y+3/4,-x+3/4;x+1/4,y+3/4,-z+3/4;-x+1/4,y+3/4,z+3/4;x+1/4,-y+3/4,z+3/4;-y,-x+1/2,z+1/2;y,x+1/2,z+1/2;"
	symLines[226] += "x,-z+1/2,-y+1/2;x,z+1/2,y+1/2;-z,y+1/2,-x+1/2;z,y+1/2,x+1/2;x+1/2,y,z+1/2;-x+3/4,-y+1/4,-z+3/4;"
	symLines[226] += "-y+3/4,x+1/4,z+3/4;-x+1/2,-y,z+1/2;y+3/4,-x+1/4,z+3/4;x+3/4,-z+1/4,y+3/4;x+1/2,-y,-z+1/2;x+3/4,z+1/4,-y+3/4;"
	symLines[226] += "z+3/4,y+1/4,-x+3/4;-x+1/2,y,-z+1/2;-z+3/4,y+1/4,x+3/4;y+1/2,-x,-z+1/2;-y+1/2,x,-z+1/2;-x+1/2,z,-y+1/2;"
	symLines[226] += "-x+1/2,-z,y+1/2;-z+1/2,-y,x+1/2;z+1/2,-y,-x+1/2;z+1/2,x,y+1/2;y+1/2,z,x+1/2;-y+1/2,-z,x+1/2;z+1/2,-x,-y+1/2;"
	symLines[226] += "-y+1/2,z,-x+1/2;-z+1/2,-x,y+1/2;-z+1/2,x,-y+1/2;y+1/2,-z,-x+1/2;-z+3/4,-x+1/4,-y+3/4;-y+3/4,-z+1/4,-x+3/4;"
	symLines[226] += "y+3/4,z+1/4,-x+3/4;-z+3/4,x+1/4,y+3/4;y+3/4,-z+1/4,x+3/4;z+3/4,x+1/4,-y+3/4;z+3/4,-x+1/4,y+3/4;"
	symLines[226] += "-y+3/4,z+1/4,x+3/4;y+3/4,x+1/4,-z+3/4;-y+3/4,-x+1/4,-z+3/4;-x+3/4,z+1/4,y+3/4;-x+3/4,-z+1/4,-y+3/4;"
	symLines[226] += "z+3/4,-y+1/4,x+3/4;-z+3/4,-y+1/4,-x+3/4;x+3/4,y+1/4,-z+3/4;-x+3/4,y+1/4,z+3/4;x+3/4,-y+1/4,z+3/4;"
	symLines[226] += "-y+1/2,-x,z+1/2;y+1/2,x,z+1/2;x+1/2,-z,-y+1/2;x+1/2,z,y+1/2;-z+1/2,y,-x+1/2;z+1/2,y,x+1/2;x+1/2,y+1/2,z;"
	symLines[226] += "-x+3/4,-y+3/4,-z+1/4;-y+3/4,x+3/4,z+1/4;-x+1/2,-y+1/2,z;y+3/4,-x+3/4,z+1/4;x+3/4,-z+3/4,y+1/4;x+1/2,-y+1/2,-z;"
	symLines[226] += "x+3/4,z+3/4,-y+1/4;z+3/4,y+3/4,-x+1/4;-x+1/2,y+1/2,-z;-z+3/4,y+3/4,x+1/4;y+1/2,-x+1/2,-z;-y+1/2,x+1/2,-z;"
	symLines[226] += "-x+1/2,z+1/2,-y;-x+1/2,-z+1/2,y;-z+1/2,-y+1/2,x;z+1/2,-y+1/2,-x;z+1/2,x+1/2,y;y+1/2,z+1/2,x;-y+1/2,-z+1/2,x;"
	symLines[226] += "z+1/2,-x+1/2,-y;-y+1/2,z+1/2,-x;-z+1/2,-x+1/2,y;-z+1/2,x+1/2,-y;y+1/2,-z+1/2,-x;-z+3/4,-x+3/4,-y+1/4;"
	symLines[226] += "-y+3/4,-z+3/4,-x+1/4;y+3/4,z+3/4,-x+1/4;-z+3/4,x+3/4,y+1/4;y+3/4,-z+3/4,x+1/4;z+3/4,x+3/4,-y+1/4;"
	symLines[226] += "z+3/4,-x+3/4,y+1/4;-y+3/4,z+3/4,x+1/4;y+3/4,x+3/4,-z+1/4;-y+3/4,-x+3/4,-z+1/4;-x+3/4,z+3/4,y+1/4;"
	symLines[226] += "-x+3/4,-z+3/4,-y+1/4;z+3/4,-y+3/4,x+1/4;-z+3/4,-y+3/4,-x+1/4;x+3/4,y+3/4,-z+1/4;-x+3/4,y+3/4,z+1/4;"
	symLines[226] += "x+3/4,-y+3/4,z+1/4;-y+1/2,-x+1/2,z;y+1/2,x+1/2,z;x+1/2,-z+1/2,-y;x+1/2,z+1/2,y;-z+1/2,y+1/2,-x;z+1/2,y+1/2,x"
	symLines[227]  = "x,y,z;-x+1/4,-y+1/4,-z+3/4;-y+1/4,x+1/4,z+1/4;-x,-y,z;y+1/4,-x+1/4,z+1/4;x+1/4,-z+1/4,y+1/4;x,-y,-z;"
	symLines[227] += "x+1/4,z+1/4,-y+1/4;z+1/4,y+1/4,-x+1/4;-x,y,-z;-z+1/4,y+1/4,x+1/4;y,-x,-z+1/2;-y,x,-z+1/2;-x,z,-y+1/2;"
	symLines[227] += "-x,-z,y+1/2;-z,-y,x+1/2;z,-y,-x+1/2;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;-z,-x,y;-z,x,-y;y,-z,-x;"
	symLines[227] += "-z+1/4,-x+1/4,-y+3/4;-y+1/4,-z+1/4,-x+3/4;y+1/4,z+1/4,-x+3/4;-z+1/4,x+1/4,y+3/4;y+1/4,-z+1/4,x+3/4;"
	symLines[227] += "z+1/4,x+1/4,-y+3/4;z+1/4,-x+1/4,y+3/4;-y+1/4,z+1/4,x+3/4;y+1/4,x+1/4,-z+1/4;-y+1/4,-x+1/4,-z+1/4;"
	symLines[227] += "-x+1/4,z+1/4,y+1/4;-x+1/4,-z+1/4,-y+1/4;z+1/4,-y+1/4,x+1/4;-z+1/4,-y+1/4,-x+1/4;x+1/4,y+1/4,-z+3/4;"
	symLines[227] += "-x+1/4,y+1/4,z+3/4;x+1/4,-y+1/4,z+3/4;-y,-x,z+1/2;y,x,z+1/2;x,-z,-y+1/2;x,z,y+1/2;-z,y,-x+1/2;z,y,x+1/2;"
	symLines[227] += "x,y+1/2,z+1/2;-x+1/4,-y+3/4,-z+1/4;-y+1/4,x+3/4,z+3/4;-x,-y+1/2,z+1/2;y+1/4,-x+3/4,z+3/4;x+1/4,-z+3/4,y+3/4;"
	symLines[227] += "x,-y+1/2,-z+1/2;x+1/4,z+3/4,-y+3/4;z+1/4,y+3/4,-x+3/4;-x,y+1/2,-z+1/2;-z+1/4,y+3/4,x+3/4;y,-x+1/2,-z;"
	symLines[227] += "-y,x+1/2,-z;-x,z+1/2,-y;-x,-z+1/2,y;-z,-y+1/2,x;z,-y+1/2,-x;z,x+1/2,y+1/2;y,z+1/2,x+1/2;-y,-z+1/2,x+1/2;"
	symLines[227] += "z,-x+1/2,-y+1/2;-y,z+1/2,-x+1/2;-z,-x+1/2,y+1/2;-z,x+1/2,-y+1/2;y,-z+1/2,-x+1/2;-z+1/4,-x+3/4,-y+1/4;"
	symLines[227] += "-y+1/4,-z+3/4,-x+1/4;y+1/4,z+3/4,-x+1/4;-z+1/4,x+3/4,y+1/4;y+1/4,-z+3/4,x+1/4;z+1/4,x+3/4,-y+1/4;"
	symLines[227] += "z+1/4,-x+3/4,y+1/4;-y+1/4,z+3/4,x+1/4;y+1/4,x+3/4,-z+3/4;-y+1/4,-x+3/4,-z+3/4;-x+1/4,z+3/4,y+3/4;"
	symLines[227] += "-x+1/4,-z+3/4,-y+3/4;z+1/4,-y+3/4,x+3/4;-z+1/4,-y+3/4,-x+3/4;x+1/4,y+3/4,-z+1/4;-x+1/4,y+3/4,z+1/4;"
	symLines[227] += "x+1/4,-y+3/4,z+1/4;-y,-x+1/2,z;y,x+1/2,z;x,-z+1/2,-y;x,z+1/2,y;-z,y+1/2,-x;z,y+1/2,x;x+1/2,y,z+1/2;"
	symLines[227] += "-x+3/4,-y+1/4,-z+1/4;-y+3/4,x+1/4,z+3/4;-x+1/2,-y,z+1/2;y+3/4,-x+1/4,z+3/4;x+3/4,-z+1/4,y+3/4;x+1/2,-y,-z+1/2;"
	symLines[227] += "x+3/4,z+1/4,-y+3/4;z+3/4,y+1/4,-x+3/4;-x+1/2,y,-z+1/2;-z+3/4,y+1/4,x+3/4;y+1/2,-x,-z;-y+1/2,x,-z;-x+1/2,z,-y;"
	symLines[227] += "-x+1/2,-z,y;-z+1/2,-y,x;z+1/2,-y,-x;z+1/2,x,y+1/2;y+1/2,z,x+1/2;-y+1/2,-z,x+1/2;z+1/2,-x,-y+1/2;"
	symLines[227] += "-y+1/2,z,-x+1/2;-z+1/2,-x,y+1/2;-z+1/2,x,-y+1/2;y+1/2,-z,-x+1/2;-z+3/4,-x+1/4,-y+1/4;-y+3/4,-z+1/4,-x+1/4;"
	symLines[227] += "y+3/4,z+1/4,-x+1/4;-z+3/4,x+1/4,y+1/4;y+3/4,-z+1/4,x+1/4;z+3/4,x+1/4,-y+1/4;z+3/4,-x+1/4,y+1/4;"
	symLines[227] += "-y+3/4,z+1/4,x+1/4;y+3/4,x+1/4,-z+3/4;-y+3/4,-x+1/4,-z+3/4;-x+3/4,z+1/4,y+3/4;-x+3/4,-z+1/4,-y+3/4;"
	symLines[227] += "z+3/4,-y+1/4,x+3/4;-z+3/4,-y+1/4,-x+3/4;x+3/4,y+1/4,-z+1/4;-x+3/4,y+1/4,z+1/4;x+3/4,-y+1/4,z+1/4;-y+1/2,-x,z;"
	symLines[227] += "y+1/2,x,z;x+1/2,-z,-y;x+1/2,z,y;-z+1/2,y,-x;z+1/2,y,x;x+1/2,y+1/2,z;-x+3/4,-y+3/4,-z+3/4;-y+3/4,x+3/4,z+1/4;"
	symLines[227] += "-x+1/2,-y+1/2,z;y+3/4,-x+3/4,z+1/4;x+3/4,-z+3/4,y+1/4;x+1/2,-y+1/2,-z;x+3/4,z+3/4,-y+1/4;z+3/4,y+3/4,-x+1/4;"
	symLines[227] += "-x+1/2,y+1/2,-z;-z+3/4,y+3/4,x+1/4;y+1/2,-x+1/2,-z+1/2;-y+1/2,x+1/2,-z+1/2;-x+1/2,z+1/2,-y+1/2;"
	symLines[227] += "-x+1/2,-z+1/2,y+1/2;-z+1/2,-y+1/2,x+1/2;z+1/2,-y+1/2,-x+1/2;z+1/2,x+1/2,y;y+1/2,z+1/2,x;-y+1/2,-z+1/2,x;"
	symLines[227] += "z+1/2,-x+1/2,-y;-y+1/2,z+1/2,-x;-z+1/2,-x+1/2,y;-z+1/2,x+1/2,-y;y+1/2,-z+1/2,-x;-z+3/4,-x+3/4,-y+3/4;"
	symLines[227] += "-y+3/4,-z+3/4,-x+3/4;y+3/4,z+3/4,-x+3/4;-z+3/4,x+3/4,y+3/4;y+3/4,-z+3/4,x+3/4;z+3/4,x+3/4,-y+3/4;"
	symLines[227] += "z+3/4,-x+3/4,y+3/4;-y+3/4,z+3/4,x+3/4;y+3/4,x+3/4,-z+1/4;-y+3/4,-x+3/4,-z+1/4;-x+3/4,z+3/4,y+1/4;"
	symLines[227] += "-x+3/4,-z+3/4,-y+1/4;z+3/4,-y+3/4,x+1/4;-z+3/4,-y+3/4,-x+1/4;x+3/4,y+3/4,-z+3/4;-x+3/4,y+3/4,z+3/4;"
	symLines[227] += "x+3/4,-y+3/4,z+3/4;-y+1/2,-x+1/2,z+1/2;y+1/2,x+1/2,z+1/2;x+1/2,-z+1/2,-y+1/2;x+1/2,z+1/2,y+1/2;"
	symLines[227] += "-z+1/2,y+1/2,-x+1/2;z+1/2,y+1/2,x+1/2"
	symLines[228]  = "x,y,z;-y,x,z;-x,-y,z;y,-x,z;x,-z,y;x,-y,-z;x,z,-y;z,y,-x;-x,y,-z;-z,y,x;z,x,y;y,z,x;-y,-z,x;z,-x,-y;-y,z,-x;"
	symLines[228] += "-z,-x,y;-z,x,-y;y,-z,-x;y,x,-z;-y,-x,-z;-x,z,y;-x,-z,-y;z,-y,x;-z,-y,-x;-x,-y,-z;y,-x,-z;x,y,-z;-y,x,-z;"
	symLines[228] += "-x,z,-y;-x,y,z;-x,-z,y;-z,-y,x;x,-y,z;z,-y,-x;-z,-x,-y;-y,-z,-x;y,z,-x;-z,x,y;y,-z,x;z,x,-y;z,-x,y;-y,z,x;"
	symLines[228] += "-y,-x,z;y,x,z;x,-z,-y;x,z,y;-z,y,-x;z,y,x;x+1/2,y+1/2,z+1/2;-y+1/2,x+1/2,z+1/2;-x+1/2,-y+1/2,z+1/2;"
	symLines[228] += "y+1/2,-x+1/2,z+1/2;x+1/2,-z+1/2,y+1/2;x+1/2,-y+1/2,-z+1/2;x+1/2,z+1/2,-y+1/2;z+1/2,y+1/2,-x+1/2;"
	symLines[228] += "-x+1/2,y+1/2,-z+1/2;-z+1/2,y+1/2,x+1/2;z+1/2,x+1/2,y+1/2;y+1/2,z+1/2,x+1/2;-y+1/2,-z+1/2,x+1/2;"
	symLines[228] += "z+1/2,-x+1/2,-y+1/2;-y+1/2,z+1/2,-x+1/2;-z+1/2,-x+1/2,y+1/2;-z+1/2,x+1/2,-y+1/2;y+1/2,-z+1/2,-x+1/2;"
	symLines[228] += "y+1/2,x+1/2,-z+1/2;-y+1/2,-x+1/2,-z+1/2;-x+1/2,z+1/2,y+1/2;-x+1/2,-z+1/2,-y+1/2;z+1/2,-y+1/2,x+1/2;"
	symLines[228] += "-z+1/2,-y+1/2,-x+1/2;-x+1/2,-y+1/2,-z+1/2;y+1/2,-x+1/2,-z+1/2;x+1/2,y+1/2,-z+1/2;-y+1/2,x+1/2,-z+1/2;"
	symLines[228] += "-x+1/2,z+1/2,-y+1/2;-x+1/2,y+1/2,z+1/2;-x+1/2,-z+1/2,y+1/2;-z+1/2,-y+1/2,x+1/2;x+1/2,-y+1/2,z+1/2;"
	symLines[228] += "z+1/2,-y+1/2,-x+1/2;-z+1/2,-x+1/2,-y+1/2;-y+1/2,-z+1/2,-x+1/2;y+1/2,z+1/2,-x+1/2;-z+1/2,x+1/2,y+1/2;"
	symLines[228] += "y+1/2,-z+1/2,x+1/2;z+1/2,x+1/2,-y+1/2;z+1/2,-x+1/2,y+1/2;-y+1/2,z+1/2,x+1/2;-y+1/2,-x+1/2,z+1/2;"
	symLines[228] += "y+1/2,x+1/2,z+1/2;x+1/2,-z+1/2,-y+1/2;x+1/2,z+1/2,y+1/2;-z+1/2,y+1/2,-x+1/2;z+1/2,y+1/2,x+1/2"
	symLines[229]  = "x,y,z;-y+1/4,x+3/4,z+1/4;-x,-y+1/2,z;y+1/4,-x+1/4,z+3/4;x+1/4,-z+1/4,y+3/4;x,-y,-z+1/2;x+3/4,z+1/4,-y+1/4;"
	symLines[229] += "z+3/4,y+1/4,-x+1/4;-x+1/2,y,-z;-z+1/4,y+3/4,x+1/4;z,x,y;y,z,x;-y,-z+1/2,x;z,-x,-y+1/2;-y+1/2,z,-x;-z,-x+1/2,y;"
	symLines[229] += "-z+1/2,x,-y;y,-z,-x+1/2;y+3/4,x+1/4,-z+1/4;-y+1/4,-x+1/4,-z+1/4;-x+1/4,z+3/4,y+1/4;-x+1/4,-z+1/4,-y+1/4;"
	symLines[229] += "z+1/4,-y+1/4,x+3/4;-z+1/4,-y+1/4,-x+1/4;-x,-y,-z;y+3/4,-x+1/4,-z+3/4;x,y+1/2,-z;-y+3/4,x+3/4,-z+1/4;"
	symLines[229] += "-x+3/4,z+3/4,-y+1/4;-x,y,z+1/2;-x+1/4,-z+3/4,y+3/4;-z+1/4,-y+3/4,x+3/4;x+1/2,-y,z;z+3/4,-y+1/4,-x+3/4;-z,-x,-y;"
	symLines[229] += "-y,-z,-x;y,z+1/2,-x;-z,x,y+1/2;y+1/2,-z,x;z,x+1/2,-y;z+1/2,-x,y;-y,z,x+1/2;-y+1/4,-x+3/4,z+3/4;"
	symLines[229] += "y+3/4,x+3/4,z+3/4;x+3/4,-z+1/4,-y+3/4;x+3/4,z+3/4,y+3/4;-z+3/4,y+3/4,-x+1/4;z+3/4,y+3/4,x+3/4;"
	symLines[229] += "x+1/2,y+1/2,z+1/2;-y+3/4,x+1/4,z+3/4;-x+1/2,-y,z+1/2;y+3/4,-x+3/4,z+1/4;x+3/4,-z+3/4,y+1/4;x+1/2,-y+1/2,-z;"
	symLines[229] += "x+1/4,z+3/4,-y+3/4;z+1/4,y+3/4,-x+3/4;-x,y+1/2,-z+1/2;-z+3/4,y+1/4,x+3/4;z+1/2,x+1/2,y+1/2;y+1/2,z+1/2,x+1/2;"
	symLines[229] += "-y+1/2,-z,x+1/2;z+1/2,-x+1/2,-y;-y,z+1/2,-x+1/2;-z+1/2,-x,y+1/2;-z,x+1/2,-y+1/2;y+1/2,-z+1/2,-x;"
	symLines[229] += "y+1/4,x+3/4,-z+3/4;-y+3/4,-x+3/4,-z+3/4;-x+3/4,z+1/4,y+3/4;-x+3/4,-z+3/4,-y+3/4;z+3/4,-y+3/4,x+1/4;"
	symLines[229] += "-z+3/4,-y+3/4,-x+3/4;-x+1/2,-y+1/2,-z+1/2;y+1/4,-x+3/4,-z+1/4;x+1/2,y,-z+1/2;-y+1/4,x+1/4,-z+3/4;"
	symLines[229] += "-x+1/4,z+1/4,-y+3/4;-x+1/2,y+1/2,z;-x+3/4,-z+1/4,y+1/4;-z+3/4,-y+1/4,x+1/4;x,-y+1/2,z+1/2;z+1/4,-y+3/4,-x+1/4;"
	symLines[229] += "-z+1/2,-x+1/2,-y+1/2;-y+1/2,-z+1/2,-x+1/2;y+1/2,z,-x+1/2;-z+1/2,x+1/2,y;y,-z+1/2,x+1/2;z+1/2,x,-y+1/2;"
	symLines[229] += "z,-x+1/2,y+1/2;-y+1/2,z+1/2,x;-y+3/4,-x+1/4,z+1/4;y+1/4,x+1/4,z+1/4;x+1/4,-z+3/4,-y+1/4;x+1/4,z+1/4,y+1/4;"
	symLines[229] += "-z+1/4,y+1/4,-x+3/4;z+1/4,y+1/4,x+1/4"
	String str = symLines[SpaceGroup-1]
	KillWaves/Z symLines_temp__
	return str
End



Static Function WyckoffMultiplicity(SG,letter)
	Variable SG
	String letter
	String list = WyckoffStrFromSG(SG)
	if (strlen(list)<1)
		return 1
	endif
	Variable mult = NumberByKey(letter,list,":",";",1)
	return mult
End
//
Static Function/T WyckoffMenuStr(SG)
	Variable SG
	String list = WyckoffStrFromSG(SG)
	if (strlen(list)<1)
		return ""
	endif

	Variable i, N=ItemsInList(list)
	String mStr=""
	for (i=1;i<=N;i+=1)
		mStr += num2char(i<27 ? i+96 : i+38)+";"
	endfor
	return mStr
End
//
Static Function/T WyckoffStrFromSG(SG)
	Variable SG
	if (!isValidSpaceGroup(SG))
		return ""
	endif
	SG = round(SG)

	Make/N=230/FREE/T WyckoffWave=""
	WyckoffWave[0] = {"a:1;","a:1;b:1;c:1;d:1;e:1;f:1;g:1;h:1;i:2;","a:1;b:1;c:1;d:1;e:2;","a:2;","a:2;b:2;c:4;","a:1;b:1;c:2;"}
	WyckoffWave[6] = {"a:2;","a:2;b:4;","a:4;",":1;b:1;c:1;d:1;e:1;f:1;g:1;h:1;i:2;j:2;k:2;l:2;m:2;n:2;o:4;","a:2;b:2;c:2;d:2;e:2;f:4;"}
	WyckoffWave[11] = {"a:2;b:2;c:2;d:2;e:4;f:4;g:4;h:4;i:4;j:8;","a:2;b:2;c:2;d:2;e:2;f:2;g:4;","a:2;b:2;c:2;d:2;e:4;"}
	WyckoffWave[14] = {"a:4;b:4;c:4;d:4;e:4;f:8;","a:1;b:1;c:1;d:1;e:1;f:1;g:1;h:1;i:2;j:2;k:2;l:2;m:2;n:2;o:2;p:2;q:2;r:2;s:2;t:2;u:4;"}
	WyckoffWave[16] = {"a:2;b:2;c:2;d:2;e:4;","a:2;b:2;c:4;","a:4;","a:4;b:4;c:8;","a:2;b:2;c:2;d:2;e:4;f:4;g:4;h:4;i:4;j:4;k:4;l:8;"}
	WyckoffWave[21] = {"a:4;b:4;c:4;d:4;e:8;j:8;f:8;i:8;g:8;h:8;k:16;","a:2;b:2;c:2;d:2;e:4;f:4;g:4;h:4;i:4;j:4;k:8;","a:4;b:4;c:4;d:8;"}
	WyckoffWave[24] = {"a:1;b:1;c:1;d:1;e:2;f:2;g:2;h:2;i:4;","a:2;b:2;c:4;","a:2;b:2;c:2;d:2;e:4;","a:2;b:2;c:2;d:4;","a:4;","a:2;b:2;c:4;"}
	WyckoffWave[30] = {"a:2;b:4;","a:2;b:2;c:4;","a:4;","a:2;b:2;c:4;","a:2;b:2;c:4;d:4;e:4;f:8;","a:4;b:8;","a:4;b:4;c:4;d:8;"}
	WyckoffWave[37] = {"a:2;b:2;c:4;d:4;e:4;f:8;","a:4;b:4;c:4;d:8;","a:4;b:4;c:8;"}
	WyckoffWave[40] = {"a:4;b:8;","a:4;b:8;c:8;d:8;e:16;","a:8;b:16;","a:2;b:2;c:4;d:4;e:8;","a:4;b:4;c:8;","a:4;b:4;c:8;"}
	WyckoffWave[46] = {"a:1;b:1;c:1;d:1;e:1;f:1;g:1;h:1;i:2;j:2;k:2;l:2;m:2;n:2;o:2;p:2;q:2;r:2;s:2;t:2;u:4;v:4;w:4;x:4;y:4;z:4;A:8;"}
	WyckoffWave[47] = {"a:2;b:2;c:2;d:2;e:4;f:4;g:4;h:4;i:4;j:4;k:4;l:4;m:8;","a:2;b:2;c:2;d:2;e:2;f:2;g:2;h:2;i:4;j:4;k:4;l:4;m:4;n:4;o:4;p:4;q:4;r:8;"}
	WyckoffWave[49] = {"a:2;b:2;c:2;d:2;e:4;f:4;g:4;h:4;i:4;j:4;k:4;l:4;m:8;","a:2;b:2;c:2;d:2;e:2;f:2;g:4;h:4;i:4;j:4;k:4;l:8;","a:4;b:4;c:4;d:4;e:8;"}
	WyckoffWave[52] = {"a:2;b:2;c:2;d:2;e:4;f:4;g:4;h:4;i:8;","a:4;b:4;c:4;d:4;e:4;f:8;","a:2;b:2;c:2;d:2;e:4;f:4;g:4;h:4;i:8;","a:4;b:4;c:4;d:4;e:8;"}
	WyckoffWave[56] = {"a:4;b:4;c:4;d:4;e:8;","a:2;b:2;c:2;d:2;e:4;f:4;g:4;h:8;","a:2;b:2;c:4;d:4;e:4;f:4;g:8;","a:4;b:4;c:4;d:8;"}
	WyckoffWave[60] = {"a:4;b:4;c:8;","a:4;b:4;c:4;d:8;","a:4;b:4;c:4;d:8;e:8;f:8;g:8;h:16;","a:4;b:4;c:8;d:8;e:8;f:8;g:16;"}
	WyckoffWave[64] = {"a:2;b:2;c:2;d:2;e:4;f:4;g:4;h:4;i:4;j:4;k:4;l:4;m:8;n:8;o:8;p:8;q:8;r:16;","a:4;b:4;c:4;d:4;e:4;f:4;g:8;h:8;i:8;j:8;k:8;l:8;m:16;"}
	WyckoffWave[66] = {"a:4;b:4;c:4;d:4;e:4;f:4;g:4;h:8;i:8;j:8;k:8;l:8;m:8;n:8;o:16;","a:4;b:4;c:8;d:8;e:8;f:8;g:8;h:8;i:16;"}
	WyckoffWave[68] = {"a:4;b:4;c:8;d:8;e:8;f:8;g:8;h:8;i:8;j:16;k:16;l:16;m:16;n:16;o:16;p:32;","a:8;b:8;c:16;d:16;e:16;f:16;g:16;h:32;"}
	WyckoffWave[70] = {"a:2;b:2;c:2;d:2;e:4;f:4;g:4;h:4;i:4;j:4;k:8;l:8;m:8;n:8;o:16;","a:4;b:4;c:4;d:4;e:8;f:8;g:8;h:8;i:8;j:8;k:16;"}
	WyckoffWave[72] = {"a:8;b:8;c:8;d:8;e:8;f:16;","a:4;b:4;c:4;d:4;e:4;f:8;g:8;h:8;i:8;j:16;","a:1;b:1;c:2;d:4;","a:4;","a:2;b:2;c:2;d:4;","a:4;"}
	WyckoffWave[78] = {"a:2;b:4;c:8;","a:4;b:8;","a:1;b:1;c:1;d:1;e:2;f:2;g:2;h:4;","a:2;b:2;c:2;d:2;e:4;f:4;g:8;"}
	WyckoffWave[82] = {"a:1;b:1;c:1;d:1;e:2;f:2;g:2;h:2;i:4;j:4;k:4;l:8;","a:2;b:2;c:2;d:2;e:2;f:2;g:4;h:4;i:4;j:4;k:8;"}
	WyckoffWave[84] = {"a:2;b:2;c:2;d:4;e:4;f:4;g:8;","a:2;b:2;c:4;d:4;e:4;f:4;g:8;","a:2;b:2;c:4;d:4;e:4;f:8;g:8;h:8;i:16;"}
	WyckoffWave[87] = {"a:4;b:4;c:8;d:8;e:8;f:16;","a:1;b:1;c:1;d:1;e:2;f:2;g:2;h:2;i:4;j:4;k:4;l:4;m:4;n:4;o:4;p:8;","a:2;b:2;c:2;d:4;e:4;f:4;g:8;"}
	WyckoffWave[90] = {"a:4;b:4;c:4;d:8;","a:4;b:8;","a:2;b:2;c:2;d:2;e:2;f:2;g:4;h:4;i:4;j:4;k:4;l:4;m:4;n:4;o:4;p:8;","a:2;b:2;c:4;d:4;e:4;f:4;g:8;"}
	WyckoffWave[94] = {"a:4;b:4;c:4;d:8;","a:4;b:8;","a:2;b:2;c:4;d:4;e:4;f:8;g:8;h:8;i:8;j:8;k:16;","a:4;b:4;c:8;d:8;e:8;f:8;g:16;"}
	WyckoffWave[98] = {"a:1;b:1;c:2;d:4;e:4;f:4;g:8;","a:2;b:2;c:4;d:8;","a:2;b:2;c:4;d:4;e:8;","a:2;b:4;c:4;d:8;","a:2;b:2;c:4;d:8;"}
	WyckoffWave[103] = {"a:2;b:4;c:8;","a:2;b:2;c:2;d:4;e:4;f:8;","a:4;b:4;c:8;","a:2;b:4;c:8;d:8;e:16;","a:4;b:4;c:8;d:16;"}
	WyckoffWave[108] = {"a:4;b:8;c:16;","a:8;b:16;","a:1;b:1;c:1;d:1;e:2;f:2;g:2;h:2;i:4;j:4;k:4;l:4;m:4;n:4;o:8;"}
	WyckoffWave[111] = {"a:2;c:2;b:2;d:2;e:2;f:2;g:4;h:4;i:4;j:4;k:4;l:4;m:4;n:8;","a:2;b:2;c:2;d:4;e:4;f:8;","a:2;b:2;c:4;d:4;e:8;"}
	WyckoffWave[114] = {"a:1;b:1;c:1;d:1;e:2;f:2;g:2;h:4;i:4;j:4;k:4;l:8;","a:2;b:2;c:2;d:2;e:4;f:4;g:4;h:4;i:4;j:8;"}
	WyckoffWave[116] = {"a:2;b:2;c:2;d:2;e:4;f:4;g:4;h:4;i:8;","a:2;b:2;c:2;d:2;e:4;f:4;g:4;h:4;i:8;"}
	WyckoffWave[118] = {"a:2;b:2;c:2;d:2;e:4;f:4;g:8;h:8;i:8;j:16;","a:4;d:4;b:4;c:4;e:8;h:8;f:8;g:8;i:16;"}
	WyckoffWave[120] = {"a:2;b:2;c:4;d:4;e:4;f:8;g:8;h:8;i:8;j:16;","a:4;b:4;c:8;d:8;e:16;"}
	WyckoffWave[122] = {"a:1;b:1;c:1;d:1;e:2;f:2;g:2;h:2;i:4;j:4;k:4;l:4;m:4;n:4;o:4;p:8;q:8;r:8;s:8;t:8;u:16;"}
	WyckoffWave[123] = {"a:2;c:2;b:2;d:2;e:4;f:4;g:4;h:4;i:8;j:8;k:8;l:8;m:8;n:16;"}
	WyckoffWave[124] = {"a:2;b:2;c:2;d:2;e:4;f:4;g:4;h:4;i:8;j:8;k:8;l:8;m:8;n:16;","a:2;b:2;c:4;d:4;e:4;f:8;g:8;h:8;i:8;j:8;k:16;"}
	WyckoffWave[126] = {"a:2;b:2;c:2;d:2;e:4;f:4;g:4;h:4;i:8;j:8;k:8;l:16;","a:2;b:2;c:4;d:4;e:4;f:8;g:8;h:8;i:16;"}
	WyckoffWave[128] = {"a:2;b:2;c:2;d:4;e:4;f:4;g:8;h:8;i:8;j:8;k:16;","a:4;b:4;c:4;d:8;e:8;f:8;g:16;"}
	WyckoffWave[130] = {"a:2;b:2;c:2;d:2;e:2;f:2;g:4;h:4;i:4;j:4;k:4;l:4;m:4;n:8;o:8;p:8;q:8;r:16;","a:2;c:2;b:2;d:2;e:4;f:4;g:4;h:4;i:4;j:4;k:8;l:8;m:8;n:8;o:8;p:16;"}
	WyckoffWave[132] = {"a:4;b:4;c:4;d:4;e:8;f:8;g:8;h:8;i:8;j:8;k:16;","a:2;b:2;c:4;d:4;e:4;f:4;g:4;h:8;i:8;j:8;k:8;l:8;m:8;n:16;"}
	WyckoffWave[134] = {"a:4;b:4;c:4;d:4;e:8;f:8;g:8;h:8;i:16;","a:2;b:2;c:4;d:4;e:4;f:4;g:4;h:8;i:8;j:8;k:16;"}
	WyckoffWave[136] = {"a:2;b:2;c:4;d:4;e:8;f:8;g:8;h:16;","a:4;b:4;c:4;d:4;e:4;f:8;g:8;h:8;i:8;j:16;"}
	WyckoffWave[138] = {"a:2;b:2;c:4;d:4;e:4;f:8;g:8;h:8;i:8;j:8;k:16;l:16;m:16;n:16;o:32;","a:4;b:4;c:4;d:4;e:8;f:8;g:8;h:8;i:16;j:16;k:16;l:16;m:32;"}
	WyckoffWave[140] = {"a:4;b:4;c:8;d:8;e:8;f:16;g:16;h:16;i:32;","a:8;b:8;c:16;d:16;e:16;f:16;g:32;","a:1;b:1;c:1;d:3;","a:3;"}
	WyckoffWave[144] = {"a:3;","a:3;b:9;","a:1;b:1;c:2;d:2;e:3;f:3;g:6;","a:3;b:3;c:6;d:9;e:9;f:18;"}
	WyckoffWave[148] = {"a:1;b:1;c:1;d:1;e:1;f:1;g:2;h:2;i:2;j:3;k:3;l:6;","a:1;b:1;c:2;d:2;e:3;f:3;g:6;","a:3;b:3;c:6;","a:3;b:3;c:6;"}
	WyckoffWave[152] = {"a:3;b:3;c:6;","a:3;b:3;c:6;","a:3;b:3;c:6;d:9;e:9;f:18;","a:1;b:1;c:1;d:3;e:6;","a:1;b:2;c:3;d:6;","a:2;b:2;c:2;d:6;"}
	WyckoffWave[158] = {"a:2;b:2;c:6;","a:3;b:9;c:18;","a:6;b:18;","a:1;b:1;c:2;d:2;e:2;f:3;g:3;h:4;i:6;j:6;k:6;l:12;"}
	WyckoffWave[162] = {"a:2;b:2;c:2;d:2;e:4;f:4;g:6;h:6;i:12;","a:1;b:1;c:2;d:2;e:3;f:3;g:6;h:6;i:6;j:12;","a:2;b:2;c:4;d:4;e:6;f:6;g:12;"}
	WyckoffWave[165] = {"a:3;b:3;c:6;e:9;d:9;f:18;g:18;h:18;i:36;","a:6;b:6;c:12;d:18;e:18;f:36;","a:1;b:2;c:3;d:6;","a:6;","a:6;"}
	WyckoffWave[170] = {"a:3;b:3;c:6;","a:3;b:3;c:6;","a:2;b:2;c:6;","a:1;b:1;c:1;d:1;e:1;f:1;g:2;h:2;i:2;j:3;k:3;l:6;"}
	WyckoffWave[174] = {"a:1;b:1;c:2;d:2;e:2;f:3;g:3;h:4;i:6;j:6;k:6;l:12;","a:2;b:2;c:2;d:2;e:4;f:4;g:6;h:6;i:12;"}
	WyckoffWave[176] = {"a:1;b:1;c:2;d:2;e:2;f:3;g:3;h:4;i:6;j:6;k:6;l:6;m:6;n:12;","a:6;b:6;c:12;","a:6;b:6;c:12;"}
	WyckoffWave[179] = {"a:3;b:3;c:3;d:3;e:6;f:6;g:6;h:6;i:6;j:6;k:12;","a:3;b:3;c:3;d:3;e:6;f:6;g:6;h:6;i:6;j:6;k:12;"}
	WyckoffWave[181] = {"a:2;b:2;c:2;d:2;e:4;f:4;g:6;h:6;i:12;","a:1;b:2;c:3;d:6;e:6;f:12;","a:2;b:4;c:6;d:12;","a:2;b:4;c:6;d:12;"}
	WyckoffWave[185] = {"a:2;b:2;c:6;d:12;","a:1;b:1;c:1;d:1;e:1;f:1;g:2;h:2;i:2;j:3;k:3;l:6;m:6;n:6;o:12;"}
	WyckoffWave[187] = {"a:2;c:2;e:2;b:2;d:2;f:2;g:4;h:4;i:4;j:6;k:6;l:12;","a:1;b:1;c:2;d:2;e:2;f:3;g:3;h:4;i:6;j:6;k:6;l:12;"}
	WyckoffWave[189] = {"a:2;b:2;c:2;d:2;e:4;f:4;g:6;h:6;i:12;","a:1;b:1;c:2;d:2;e:2;f:3;g:3;h:4;i:6;j:6;k:6;l:6;m:6;n:12;o:12;p:12;q:12;r:24;"}
	WyckoffWave[191] = {"a:2;b:2;c:4;d:4;e:4;f:6;g:6;h:8;i:12;j:12;k:12;l:12;m:24;","a:2;b:2;c:4;d:4;e:4;f:6;g:6;h:8;i:12;j:12;k:12;l:24;"}
	WyckoffWave[193] = {"a:2;b:2;c:2;d:2;e:4;f:4;g:6;h:6;i:12;j:12;k:12;l:24;","a:1;b:1;c:3;d:3;e:4;f:6;i:6;g:6;h:6;j:12;"}
	WyckoffWave[195] = {"a:4;b:4;c:4;d:4;e:16;f:24;g:24;h:48;","a:2;b:6;c:8;d:12;e:12;f:24;","a:4;b:12;"}
	WyckoffWave[198] = {"a:8;b:12;c:24;","a:1;b:1;c:3;d:3;e:6;h:6;f:6;g:6;i:8;j:12;k:12;l:24;"}
	WyckoffWave[200] = {"a:2;b:4;c:4;d:6;e:8;f:12;g:12;h:24;","a:4;b:4;c:8;d:24;e:24;f:32;g:48;h:48;i:96;"}
	WyckoffWave[202] = {"a:8;b:8;c:16;d:16;e:32;f:48;g:96;","a:2;b:6;c:8;d:12;e:12;f:16;g:24;h:48;"}
	WyckoffWave[204] = {"a:4;b:4;c:8;d:24;","a:8;b:8;c:16;d:24;e:48;"}
	WyckoffWave[206] = {"a:1;b:1;c:3;d:3;e:6;f:6;g:8;h:12;i:12;j:12;k:24;","a:2;b:4;c:4;d:6;e:6;f:6;g:8;h:12;i:12;j:12;k:12;l:12;m:24;"}
	WyckoffWave[208] = {"a:4;b:4;c:8;d:24;e:24;f:32;g:48;h:48;i:48;j:96;","a:8;b:8;c:16;d:16;e:32;f:48;g:48;h:96;"}
	WyckoffWave[210] = {"a:2;b:6;c:8;d:12;e:12;f:16;g:24;h:24;i:24;j:48;","a:4;b:4;c:8;d:12;e:24;"}
	WyckoffWave[212] = {"a:4;b:4;c:8;d:12;e:24;","a:8;b:8;c:12;d:12;e:16;f:24;h:24;g:24;i:48;"}
	WyckoffWave[214] = {"a:1;b:1;c:3;d:3;e:4;f:6;g:6;h:12;i:12;j:24;","a:4;b:4;c:4;d:4;e:16;f:24;g:24;h:48;i:96;"}
	WyckoffWave[216] = {"a:2;b:6;c:8;d:12;e:12;f:24;g:24;h:48;","a:2;b:6;c:6;d:6;e:8;f:12;g:12;h:12;i:24;"}
	WyckoffWave[218] = {"a:8;b:8;c:24;d:24;e:32;f:48;g:48;h:96;","a:12;b:12;c:16;d:24;e:48;"}
	WyckoffWave[220] = {"a:1;b:1;c:3;d:3;e:6;f:6;g:8;h:12;i:12;j:12;k:24;l:24;m:24;n:48;","a:2;b:6;c:8;d:12;e:12;f:16;g:24;h:24;i:48;"}
	WyckoffWave[222] = {"a:2;b:6;c:6;d:6;e:8;f:12;g:12;h:12;i:16;j:24;k:24;l:48;","a:2;b:4;c:4;d:6;e:8;f:12;g:12;h:24;i:24;j:24;k:24;l:48;"}
	WyckoffWave[224] = {"a:4;b:4;c:8;d:24;e:24;f:32;g:48;h:48;i:48;j:96;k:96;l:192;","a:8;b:8;c:24;d:24;e:48;f:48;g:64;h:96;i:96;j:192;"}
	WyckoffWave[226] = {"a:8;b:8;c:16;d:16;e:32;f:48;g:96;h:96;i:192;","a:16;b:32;c:32;d:48;e:64;f:96;g:96;h:192;"}
	WyckoffWave[228] = {"a:2;b:6;c:8;d:12;e:12;f:16;g:24;h:24;i:48;j:48;k:48;l:96;","a:16;b:16;c:24;d:24;e:32;f:48;g:48;h:96;"}
	return WyckoffWave[SG-1]
End
//Function testWyckoff(SG,letter)
//	Variable SG
//	String letter
//
//	String mStr = LatticeSym#WyckoffMenuStr(SG)
//	Variable mult = LatticeSym#WyckoffMultiplicity(SG,letter)
//	print "list = ",LatticeSym#WyckoffStrFromSG(SG)
//	printf "mStr = '%s'\r",mStr
//	printf "for '%s',  mult = %g\r",letter,mult
//	return 0
//End

























// this is not neeed, already done by WyckoffMultiplicity(47,"A")

//Function MultiplicityFromWyckoff(SG,symbol)
//	// returns multiplicity for the Wyckoff symbol
//	Variable SG
//	String symbol
//
//	Wave/T WyckList=GetWyckoffSymStrings(SG)
//	Variable i
//	for (i=0;i<DimSize(WyckList,0);i+=1)
//		if (strsearch(WyckList[i][0],symbol,0)==0)
//			return str2num(WyckList[i][2])
//		endif
//	endfor
//	return 0
//End

Static Function/T FindWyckoffSymbol(SG,x0,y0,z0,mult)
	Variable SG
	Variable x0,y0,z0
	Variable &mult
	Wave/T WyckList=GetWyckoffSymStrings(SG)

	String sx=num2str(x0), sy=num2str(y0), sz=num2str(z0)
	String item, symOp, symbol=""

	Variable xop,yop,zop
	Variable m, N=DimSize(WyckList,0)
	mult = 0
	for (m=0;m<N;m+=1)
		symOp = WyckList[m][1]
		symOp = ReplaceString("2x",symOp,"2*x")
		symOp = ReplaceString("2y",symOp,"2*y")
		symOp = ReplaceString("2z",symOp,"2*z")
		symOp = ReplaceString("x",symOp,sx)
		symOp = ReplaceString("y",symOp,sy)
		symOp = ReplaceString("z",symOp,sz)

		Execute "Variable/G root:temp_temp_ = "+StringFromList(0,symOp,",")
		xop = NumVarOrDefault("root:temp_temp_",NaN)
		Execute "Variable/G root:temp_temp_ = "+StringFromList(1,symOp,",")
		yop = NumVarOrDefault("root:temp_temp_",NaN)
		Execute "Variable/G root:temp_temp_ = "+StringFromList(2,symOp,",")
		zop = NumVarOrDefault("root:temp_temp_",NaN)

		if ((abs(xop-x0) + abs(yop-y0) + abs(zop-z0)) < 1e-4)
			symbol = WyckList[m][0]
			mult = round(str2num(WyckList[m][2]))
			break
		endif
	endfor
	KillVariables/Z root:temp_temp_
	return symbol
End
//		//	test_FindWyckoffSymbol(47,  0, 0.5, 0.3765)
//	Function test_FindWyckoffSymbol(SG,x0,y0,z0)
//		Variable SG
//		Variable x0,y0,z0
//	
//		Variable mult
//		String symbol = FindWyckoffSymbol(SG,x0,y0,z0,mult)
//		printf "Wyckoff = \"%s\",  mult = %g\r",symbol,mult
//	End


Static Function ForceXYZtoWyckoff(SG,symbol,x0,y0,z0)
	Variable SG
	String symbol
	Variable &x0,&y0,&z0
	Wave/T WyckList=GetWyckoffSymStrings(SG)

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

	Execute "Variable/G root:temp_temp_ = "+StringFromList(0,symOp,",")
	x0 = NumVarOrDefault("root:temp_temp_",x0)
	Execute "Variable/G root:temp_temp_ = "+StringFromList(1,symOp,",")
	y0 = NumVarOrDefault("root:temp_temp_",y0)
	Execute "Variable/G root:temp_temp_ = "+StringFromList(2,symOp,",")
	z0 = NumVarOrDefault("root:temp_temp_",z0)
	KillVariables/Z root:temp_temp_
	return 0
End
//	Function test_ForceXYZtoWyckoff(SG,symbol,x0,y0,z0)
//		Variable SG
//		String symbol
//		Variable x0,y0,z0
//	
//		String out, str
//		sprintf out,"#%g, Wyckoff=\"%s\",   {%g, %g, %g} --> ",SG,symbol,x0,y0,z0
//		Variable err=ForceXYZtoWyckoff(SG,symbol,x0,y0,z0)
//		sprintf str,"{%g, %g, %g}",x0,y0,z0
//		out += str
//		out += SelectString(err,"","   ****** ERROR ******")
//		print out
//	End
//

Static Function/WAVE GetWyckoffSymStrings(SG)
	Variable SG
	if (!isValidSpaceGroup(SG))
		return $""
	endif

	Make/N=(230)/T/FREE WyckoffSyms
	// Triclinic [1,2]  lines 0-1
	WyckoffSyms[0] = "a:x,y,z:1;"
	WyckoffSyms[1] = "a:0,0,0:1;b:0,0,1/2:1;c:0,1/2,0:1;d:1/2,0,0:1;e:1/2,1/2,0:1;f:1/2,0,1/2:1;g:0,1/2,1/2:1;h:1/2,1/2,1/2:1;i:x,y,z:2;"
	// Monoclinic [3,15]  lines 2-14
	WyckoffSyms[2] = "a:0,y,0:1;b:0,y,1/2:1;c:1/2,y,0:1;d:1/2,y,1/2:1;e:x,y,z:2;"
	WyckoffSyms[3] = "a:x,y,z:2;"
	WyckoffSyms[4] = "a:0,y,0:2;b:0,y,1/2:2;c:x,y,z:4;"
	WyckoffSyms[5] = "a:x,0,z:1;b:x,1/2,z:1;c:x,y,z:2;"
	WyckoffSyms[6] = "a:x,y,z:2;"
	WyckoffSyms[7] = "a:x,0,z:2;b:x,y,z:4;"
	WyckoffSyms[8] = "a:x,y,z:4;"
	WyckoffSyms[9] = "a:0,0,0:1;b:0,1/2,0:1;c:0,0,1/2:1;d:1/2,0,0:1;e:1/2,1/2,0:1;f:0,1/2,1/2:1;g:1/2,0,1/2:1;h:1/2,1/2,1/2:1;i:0,y,0:2;j:1/2,y,0:2;k:0,y,1/2:2;l:1/2,y,1/2:2;m:x,0,z:2;n:x,1/2,z:2;o:x,y,z:4;"
	WyckoffSyms[10] = "a:0,0,0:2;b:1/2,0,0:2;c:0,0,1/2:2;d:1/2,0,1/2:2;e:x,1/4,z:2;f:x,y,z:4;"
	WyckoffSyms[11] = "a:0,0,0:2;b:0,1/2,0:2;c:0,0,1/2:2;d:0,1/2,1/2:2;e:1/4,1/4,0:4;f:1/4,1/4,1/2:4;g:0,y,0:4;h:0,y,1/2:4;i:x,0,z:4;j:x,y,z:8;"
	WyckoffSyms[12] = "a:0,0,0:2;b:1/2,1/2,0:2;c:0,1/2,0:2;d:1/2,0,0:2;e:0,y,1/4:2;f:1/2,y,1/4:2;g:x,y,z:4;"
	WyckoffSyms[13] = "a:0,0,0:2;b:1/2,0,0:2;c:0,0,1/2:2;d:1/2,0,1/2:2;e:x,y,z:4;"
	WyckoffSyms[14] = "a:0,0,0:4;b:0,1/2,0:4;c:1/4,1/4,0:4;d:1/4,1/4,1/2:4;e:0,y,1/4:4;f:x,y,z:8;"
	// Orthorhombic [16,74]  lines 15-73
	WyckoffSyms[15] = "a:0,0,0:1;b:1/2,0,0:1;c:0,1/2,0:1;d:0,0,1/2:1;e:1/2,1/2,0:1;f:1/2,0,1/2:1;g:0,1/2,1/2:1;h:1/2,1/2,1/2:1;i:x,0,0:2;j:x,0,1/2:2;k:x,1/2,0:2;l:x,1/2,1/2:2;m:0,y,0:2;n:0,y,1/2:2;o:1/2,y,0:2;p:1/2,y,1/2:2;q:0,0,z:2;r:1/2,0,z:2;s:0,1/2,z:2;t:1/2,1/2,z:2;u:x,y,z:4;"
	WyckoffSyms[16] = "a:x,0,0:2;b:x,1/2,0:2;c:0,y,1/4:2;d:1/2,y,1/4:2;e:x,y,z:4;"
	WyckoffSyms[17] = "a:0,0,z:2;b:0,1/2,z:2;c:x,y,z:4;"
	WyckoffSyms[18] = "a:x,y,z:4;"
	WyckoffSyms[19] = "a:x,0,0:4;b:0,y,1/4:4;c:x,y,z:8;"
	WyckoffSyms[20] = "a:0,0,0:2;b:0,1/2,0:2;c:1/2,0,1/2:2;d:0,0,1/2:2;e:x,0,0:4;f:x,0,1/2:4;g:0,y,0:4;h:0,y,1/2:4;i:0,0,z:4;j:0,1/2,z:4;k:1/4,1/4,z:4;l:x,y,z:8;"
	WyckoffSyms[21] = "a:0,0,0:4;b:0,0,1/2:4;c:1/4,1/4,1/4:4;d:1/4,1/4,3/4:4;e:x,0,0:8;f:0,y,0:8;g:0,0,z:8;h:1/4,1/4,z:8;i:1/4,y,1/4:8;j:x,1/4,1/4:8;k:x,y,z:16;"
	WyckoffSyms[22] = "a:0,0,0:2;b:1/2,0,0:2;c:0,0,1/2:2;d:0,1/2,0:2;e:x,0,0:4;f:x,0,1/2:4;g:0,y,0:4;h:1/2,y,0:4;i:0,0,z:4;j:0,1/2,z:4;k:x,y,z:8;"
	WyckoffSyms[23] = "a:x,0,1/4:4;b:1/4,y,0:4;c:0,1/4,z:4;d:x,y,z:8;"
	WyckoffSyms[24] = "a:0,0,z:1;b:0,1/2,z:1;c:1/2,0,z:1;d:1/2,1/2,z:1;e:x,0,z:2;f:x,1/2,z:2;g:0,y,z:2;h:1/2,y,z:2;i:x,y,z:4;"
	WyckoffSyms[25] = "a:0,y,z:2;b:1/2,y,z:2;c:x,y,z:4;"
	WyckoffSyms[26] = "a:0,0,z:2;b:0,1/2,z:2;c:1/2,0,z:2;d:1/2,1/2,z:2;e:x,y,z:4;"
	WyckoffSyms[27] = "a:0,0,z:2;b:0,1/2,z:2;c:1/4,y,z:2;d:x,y,z:4;"
	WyckoffSyms[28] = "a:x,y,z:4;"
	WyckoffSyms[29] = "a:0,0,z:2;b:1/2,0,z:2;c:x,y,z:4;"
	WyckoffSyms[30] = "a:0,y,z:2;b:x,y,z:4;"
	WyckoffSyms[31] = "a:0,0,z:2;b:0,1/2,z:2;c:x,y,z:4;"
	WyckoffSyms[32] = "a:x,y,z:4;"
	WyckoffSyms[33] = "a:0,0,z:2;b:0,1/2,z:2;c:x,y,z:4;"
	WyckoffSyms[34] = "a:0,0,z:2;b:0,1/2,z:2;c:1/4,1/4,z:4;d:x,0,z:4;e:0,y,z:4;f:x,y,z:8;"
	WyckoffSyms[35] = "a:0,y,z:4;b:x,y,z:8;"
	WyckoffSyms[36] = "a:0,0,z:4;b:0,1/2,z:4;c:1/4,1/4,z:4;d:x,y,z:8;"
	WyckoffSyms[37] = "a:0,0,z:2;b:1/2,0,z:2;c:x,0,z:4;d:0,y,z:4;e:1/2,y,z:4;f:x,y,z:8;"
	WyckoffSyms[38] = "a:0,0,z:4;b:1/2,0,z:4;c:x,1/4,z:4;d:x,y,z:8;"
	WyckoffSyms[39] = "a:0,0,z:4;b:1/4,y,z:4;c:x,y,z:8;"
	WyckoffSyms[40] = "a:0,0,z:4;b:x,y,z:8;"
	WyckoffSyms[41] = "a:0,0,z:4;b:1/4,1/4,z:8;c:0,y,z:8;d:x,0,z:8;e:x,y,z:16;"
	WyckoffSyms[42] = "a:0,0,z:8;b:x,y,z:16;"
	WyckoffSyms[43] = "a:0,0,z:2;b:0,1/2,z:2;c:x,0,z:4;d:0,y,z:4;e:x,y,z:8;"
	WyckoffSyms[44] = "a:0,0,z:4;b:0,1/2,z:4;c:x,y,z:8;"
	WyckoffSyms[45] = "a:0,0,z:4;b:1/4,y,z:4;c:x,y,z:8;"
	WyckoffSyms[46] = "a:0,0,0:1;b:1/2,0,0:1;c:0,0,1/2:1;d:1/2,0,1/2:1;e:0,1/2,0:1;f:1/2,1/2,0:1;g:0,1/2,1/2:1;h:1/2,1/2,1/2:1;i:x,0,0:2;j:x,0,1/2:2;k:x,1/2,0:2;l:x,1/2,1/2:2;m:0,y,0:2;n:0,y,1/2:2;o:1/2,y,0:2;p:1/2,y,1/2:2;q:0,0,z:2;r:0,1/2,z:2;s:1/2,0,z:2;t:1/2,1/2,z:2;u:0,y,z:4;v:1/2,y,z:4;w:x,0,z:4;x:x,1/2,z:4;y:x,y,0:4;z:x,y,1/2:4;A:x,y,z:8;"
	WyckoffSyms[47] = "a:1/4,1/4,1/4:2;b:3/4,1/4,1/4:2;c:1/4,1/4,3/4:2;d:1/4,3/4,1/4:2;e:1/2,1/2,1/2:4;f:0,0,0:4;g:x,1/4,1/4:4;h:x,1/4,3/4:4;i:1/4,y,1/4:4;j:3/4,y,1/4:4;k:1/4,1/4,z:4;l:1/4,3/4,z:4;m:x,y,z:8;"
	WyckoffSyms[48] = "a:0,0,0:2;b:1/2,1/2,0:2;c:0,1/2,0:2;d:1/2,0,0:2;e:0,0,1/4:2;f:1/2,0,1/4:2;g:0,1/2,1/4:2;h:1/2,1/2,1/4:2;i:x,0,1/4:4;j:x,1/2,1/4:4;k:0,y,1/4:4;l:1/2,y,1/4:4;m:0,0,z:4;n:1/2,1/2,z:4;o:0,1/2,z:4;p:1/2,0,z:4;q:x,y,0:4;r:x,y,z:8;"
	WyckoffSyms[49] = "a:1/4,1/4,0:2;b:3/4,1/4,0:2;c:3/4,1/4,1/2:2;d:1/4,1/4,1/2:2;e:0,0,0:4;f:0,0,1/2:4;g:x,1/4,0:4;h:x,1/4,1/2:4;i:1/4,y,0:4;j:1/4,y,1/2:4;k:1/4,1/4,z:4;l:1/4,3/4,z:4;m:x,y,z:8;"
	WyckoffSyms[50] = "a:0,0,0:2;b:0,1/2,0:2;c:0,0,1/2:2;d:0,1/2,1/2:2;e:1/4,0,z:2;f:1/4,1/2,z:2;g:0,y,0:4;h:0,y,1/2:4;i:x,0,z:4;j:x,1/2,z:4;k:1/4,y,z:4;l:x,y,z:8;"
	WyckoffSyms[51] = "a:0,0,0:4;b:0,0,1/2:4;c:1/4,0,z:4;d:x,1/4,1/4:4;e:x,y,z:8;"
	WyckoffSyms[52] = "a:0,0,0:2;b:1/2,0,0:2;c:1/2,1/2,0:2;d:0,1/2,0:2;e:x,0,0:4;f:x,1/2,0:4;g:1/4,y,1/4:4;h:0,y,z:4;i:x,y,z:8;"
	WyckoffSyms[53] = "a:0,0,0:4;b:0,1/2,0:4;c:0,y,1/4:4;d:1/4,0,z:4;e:1/4,1/2,z:4;f:x,y,z:8;"
	WyckoffSyms[54] = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,0:2;d:0,1/2,1/2:2;e:0,0,z:4;f:0,1/2,z:4;g:x,y,0:4;h:x,y,1/2:4;i:x,y,z:8;"
	WyckoffSyms[55] = "a:0,0,0:4;b:0,0,1/2:4;c:1/4,1/4,z:4;d:1/4,3/4,z:4;e:x,y,z:8;"
	WyckoffSyms[56] = "a:0,0,0:4;b:1/2,0,0:4;c:x,1/4,0:4;d:x,y,1/4:4;e:x,y,z:8;"
	WyckoffSyms[57] = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,0:2;d:0,1/2,1/2:2;e:0,0,z:4;f:0,1/2,z:4;g:x,y,0:4;h:x,y,z:8;"
	WyckoffSyms[58] = "a:1/4,1/4,z:2;b:1/4,3/4,z:2;c:0,0,0:4;d:0,0,1/2:4;e:1/4,y,z:4;f:x,1/4,z:4;g:x,y,z:8;"
	WyckoffSyms[59] = "a:0,0,0:4;b:0,1/2,0:4;c:0,y,1/4:4;d:x,y,z:8;"
	WyckoffSyms[60] = "a:0,0,0:4;b:0,0,1/2:4;c:x,y,z:8;"
	WyckoffSyms[61] = "a:0,0,0:4;b:0,0,1/2:4;c:x,1/4,z:4;d:x,y,z:8;"
	WyckoffSyms[62] = "a:0,0,0:4;b:0,1/2,0:4;c:0,y,1/4:4;d:1/4,1/4,0:8;e:x,0,0:8;f:0,y,z:8;g:x,y,1/4:8;h:x,y,z:16;"
	WyckoffSyms[63] = "a:0,0,0:4;b:1/2,0,0:4;c:1/4,1/4,0:8;d:x,0,0:8;e:1/4,y,1/4:8;f:0,y,z:8;g:x,y,z:16;"
	WyckoffSyms[64] = "a:0,0,0:2;b:1/2,0,0:2;c:1/2,0,1/2:2;d:0,0,1/2:2;e:1/4,1/4,0:4;f:1/4,1/4,1/2:4;g:x,0,0:4;h:x,0,1/2:4;i:0,y,0:4;j:0,y,1/2:4;k:0,0,z:4;l:0,1/2,z:4;m:1/4,1/4,z:8;n:0,y,z:8;o:x,0,z:8;p:x,y,0:8;q:x,y,1/2:8;r:x,y,z:16;"
	WyckoffSyms[65] = "a:0,0,1/4:4;b:0,1/2,1/4:4;c:0,0,0:4;d:0,1/2,0:4;e:1/4,1/4,0:4;f:1/4,3/4,0:4;g:x,0,1/4:8;h:0,y,1/4:8;i:0,0,z:8;j:0,1/2,z:8;k:1/4,1/4,z:8;l:x,y,0:8;m:x,y,z:16;"
	WyckoffSyms[66] = "a:1/4,0,0:4;b:1/4,0,1/2:4;c:0,0,0:4;d:0,0,1/2:4;e:1/4,1/4,0:4;f:1/4,1/4,1/2:4;g:0,1/4,z:4;h:x,0,0:8;i:x,0,1/2:8;j:1/4,y,0:8;k:1/4,y,1/2:8;l:1/4,0,z:8;m:0,y,z:8;n:x,1/4,z:8;o:x,y,z:16;"
	WyckoffSyms[67] = "a:0,1/4,1/4:4;b:0,1/4,3/4:4;c:1/4,3/4,0:8;d:0,0,0:8;e:x,1/4,1/4:8;f:0,y,1/4:8;g:0,1/4,z:8;h:1/4,0,z:8;i:x,y,z:16;"
	WyckoffSyms[68] = "a:0,0,0:4;b:0,0,1/2:4;c:0,1/4,1/4:8;d:1/4,0,1/4:8;e:1/4,1/4,0:8;f:1/4,1/4,1/4:8;g:x,0,0:8;h:0,y,0:8;i:0,0,z:8;j:1/4,1/4,z:16;k:1/4,y,1/4:16;l:x,1/4,1/4:16;m:0,y,z:16;n:x,0,z:16;o:x,y,0:16;p:x,y,z:32;"
	WyckoffSyms[69] = "a:1/8,1/8,1/8:8;b:1/8,1/8,5/8:8;c:0,0,0:16;d:1/2,1/2,1/2:16;e:x,1/8,1/8:16;f:1/8,y,1/8:16;g:1/8,1/8,z:16;h:x,y,z:32;"
	WyckoffSyms[70] = "a:0,0,0:2;b:0,1/2,1/2:2;c:1/2,1/2,0:2;d:1/2,0,1/2:2;e:x,0,0:4;f:x,1/2,0:4;g:0,y,0:4;h:0,y,1/2:4;i:0,0,z:4;j:1/2,0,z:4;k:1/4,1/4,1/4:8;l:0,y,z:8;m:x,0,z:8;n:x,y,0:8;o:x,y,z:16;"
	WyckoffSyms[71] = "a:0,0,1/4:4;b:1/2,0,1/4:4;c:0,0,0:4;d:1/2,0,0:4;e:1/4,1/4,1/4:8;f:x,0,1/4:8;g:0,y,1/4:8;h:0,0,z:8;i:0,1/2,z:8;j:x,y,0:8;k:x,y,z:16;"
	WyckoffSyms[72] = "a:0,0,0:8;b:1/4,1/4,1/4:8;c:x,0,1/4:8;d:1/4,y,0:8;e:0,1/4,z:8;f:x,y,z:16;"
	WyckoffSyms[73] = "a:0,0,0:4;b:0,0,1/2:4;c:1/4,1/4,1/4:4;d:1/4,1/4,3/4:4;e:0,1/4,z:4;f:x,0,0:8;g:1/4,y,1/4:8;h:0,y,z:8;i:x,1/4,z:8;j:x,y,z:16;"
	// Tetragonal [75,142]  lines 74-141
	WyckoffSyms[74] = "a:0,0,z:1;b:1/2,1/2,z:1;c:0,1/2,z:2;d:x,y,z:4;"
	WyckoffSyms[75] = "a:x,y,z:4;"
	WyckoffSyms[76] = "a:0,0,z:2;b:1/2,1/2,z:2;c:0,1/2,z:2;d:x,y,z:4;"
	WyckoffSyms[77] = "a:x,y,z:4;"
	WyckoffSyms[78] = "a:0,0,z:2;b:0,1/2,z:4;c:x,y,z:8;"
	WyckoffSyms[79] = "a:0,0,z:4;b:x,y,z:8;"
	WyckoffSyms[80] = "a:0,0,0:1;b:0,0,1/2:1;c:1/2,1/2,0:1;d:1/2,1/2,1/2:1;e:0,0,z:2;f:1/2,1/2,z:2;g:0,1/2,z:2;h:x,y,z:4;"
	WyckoffSyms[81] = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,1/4:2;d:0,1/2,3/4:2;e:0,0,z:4;f:0,1/2,z:4;g:x,y,z:8;"
	WyckoffSyms[82] = "a:0,0,0:1;b:0,0,1/2:1;c:1/2,1/2,0:1;d:1/2,1/2,1/2:1;e:0,1/2,0:2;f:0,1/2,1/2:2;g:0,0,z:2;h:1/2,1/2,z:2;i:0,1/2,z:4;j:x,y,0:4;k:x,y,1/2:4;l:x,y,z:8;"
	WyckoffSyms[83] = "a:0,0,0:2;b:1/2,1/2,0:2;c:0,1/2,0:2;d:0,1/2,1/2:2;e:0,0,1/4:2;f:1/2,1/2,1/4:2;g:0,0,z:4;h:1/2,1/2,z:4;i:0,1/2,z:4;j:x,y,0:4;k:x,y,z:8;"
	WyckoffSyms[84] = "a:1/4,3/4,0:2;b:1/4,3/4,1/2:2;c:1/4,1/4,z:2;d:0,0,0:4;e:0,0,1/2:4;f:1/4,3/4,z:4;g:x,y,z:8;"
	WyckoffSyms[85] = "a:1/4,1/4,1/4:2;b:1/4,1/4,3/4:2;c:0,0,0:4;d:0,0,1/2:4;e:3/4,1/4,z:4;f:1/4,1/4,z:4;g:x,y,z:8;"
	WyckoffSyms[86] = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,0:4;d:0,1/2,1/4:4;e:0,0,z:4;f:1/4,1/4,1/4:8;g:0,1/2,z:8;h:x,y,0:8;i:x,y,z:16;"
	WyckoffSyms[87] = "a:0,1/4,1/8:4;b:0,1/4,5/8:4;c:0,0,0:8;d:0,0,1/2:8;e:0,1/4,z:8;f:x,y,z:16;"
	WyckoffSyms[88] = "a:0,0,0:1;b:0,0,1/2:1;c:1/2,1/2,0:1;d:1/2,1/2,1/2:1;e:1/2,0,0:2;f:1/2,0,1/2:2;g:0,0,z:2;h:1/2,1/2,z:2;i:0,1/2,z:4;j:x,x,0:4;k:x,x,1/2:4;l:x,0,0:4;m:x,1/2,1/2:4;n:x,0,1/2:4;o:x,1/2,0:4;p:x,y,z:8;"
	WyckoffSyms[89] = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,z:2;d:0,0,z:4;e:x,x,0:4;f:x,x,1/2:4;g:x,y,z:8;"
	WyckoffSyms[90] = "a:0,y,0:4;b:1/2,y,0:4;c:x,x,3/8:4;d:x,y,z:8;"
	WyckoffSyms[91] = "a:x,x,0:4;b:x,y,z:8;"
	WyckoffSyms[92] = "a:0,0,0:2;b:1/2,1/2,0:2;c:0,1/2,0:2;d:0,1/2,1/2:2;e:0,0,1/4:2;f:1/2,1/2,1/4:2;g:0,0,z:4;h:1/2,1/2,z:4;i:0,1/2,z:4;j:x,0,0:4;k:x,1/2,1/2:4;l:x,0,1/2:4;m:x,1/2,0:4;n:x,x,1/4:4;o:x,x,3/4:4;p:x,y,z:8;"
	WyckoffSyms[93] = "a:0,0,0:2;b:0,0,1/2:2;c:0,0,z:4;d:0,1/2,z:4;e:x,x,0:4;f:x,x,1/2:4;g:x,y,z:8;"
	WyckoffSyms[94] = "a:0,y,0:4;b:1/2,y,0:4;c:x,x,5/8:4;d:x,y,z:8;"
	WyckoffSyms[95] = "a:x,x,0:4;b:x,y,z:8;"
	WyckoffSyms[96] = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,0:4;d:0,1/2,1/4:4;e:0,0,z:4;f:0,1/2,z:8;g:x,x,0:8;h:x,0,0:8;i:x,0,1/2:8;j:x,x+1/2,1/4:8;k:x,y,z:16;"
	WyckoffSyms[97] = "a:0,0,0:4;b:0,0,1/2:4;c:0,0,z:8;d:x,x,0:8;e:-x,x,0:8;f:x,1/4,1/8:8;g:x,y,z:16;"
	WyckoffSyms[98] = "a:0,0,z:1;b:1/2,1/2,z:1;c:1/2,0,z:2;d:x,x,z:4;e:x,0,z:4;f:x,1/2,z:4;g:x,y,z:8;"
	WyckoffSyms[99] = "a:0,0,z:2;b:1/2,0,z:2;c:x,x+1/2,z:4;d:x,y,z:8;"
	WyckoffSyms[100] = "a:0,0,z:2;b:1/2,1/2,z:2;c:0,1/2,z:4;d:x,x,z:4;e:x,y,z:8;"
	WyckoffSyms[101] = "a:0,0,z:2;b:0,1/2,z:4;c:x,x,z:4;d:x,y,z:8;"
	WyckoffSyms[102] = "a:0,0,z:2;b:1/2,1/2,z:2;c:0,1/2,z:4;d:x,y,z:8;"
	WyckoffSyms[103] = "a:0,0,z:2;b:0,1/2,z:4;c:x,y,z:8;"
	WyckoffSyms[104] = "a:0,0,z:2;b:1/2,1/2,z:2;c:0,1/2,z:2;d:x,0,z:4;e:x,1/2,z:4;f:x,y,z:8;"
	WyckoffSyms[105] = "a:0,0,z:4;b:0,1/2,z:4;c:x,y,z:8;"
	WyckoffSyms[106] = "a:0,0,z:2;b:0,1/2,z:4;c:x,x,z:8;d:x,0,z:8;e:x,y,z:16;"
	WyckoffSyms[107] = "a:0,0,z:4;b:1/2,0,z:4;c:x,x+1/2,z:8;d:x,y,z:16;"
	WyckoffSyms[108] = "a:0,0,z:4;b:0,y,z:8;c:x,y,z:16;"
	WyckoffSyms[109] = "a:0,0,z:8;b:x,y,z:16;"
	WyckoffSyms[110] = "a:0,0,0:1;b:1/2,1/2,1/2:1;c:0,0,1/2:1;d:1/2,1/2,0:1;e:1/2,0,0:2;f:1/2,0,1/2:2;g:0,0,z:2;h:1/2,1/2,z:2;i:x,0,0:4;j:x,1/2,1/2:4;k:x,0,1/2:4;l:x,1/2,0:4;m:0,1/2,z:4;n:x,x,z:4;o:x,y,z:8;"
	WyckoffSyms[111] = "a:0,0,1/4:2;b:1/2,0,1/4:2;c:1/2,1/2,1/4:2;d:0,1/2,1/4:2;e:0,0,0:2;f:1/2,1/2,0:2;g:x,0,1/4:4;h:1/2,y,1/4:4;i:x,1/2,1/4:4;j:0,y,1/4:4;k:0,0,z:4;l:1/2,1/2,z:4;m:0,1/2,z:4;n:x,y,z:8;"
	WyckoffSyms[112] = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,z:2;d:0,0,z:4;e:x,x+1/2,z:4;f:x,y,z:8;"
	WyckoffSyms[113] = "a:0,0,0:2;b:0,0,1/2:2;c:0,0,z:4;d:0,1/2,z:4;e:x,y,z:8;"
	WyckoffSyms[114] = "a:0,0,0:1;b:1/2,1/2,0:1;c:1/2,1/2,1/2:1;d:0,0,1/2:1;e:0,0,z:2;f:1/2,1/2,z:2;g:0,1/2,z:2;h:x,x,0:4;i:x,x,1/2:4;j:x,0,z:4;k:x,1/2,z:4;l:x,y,z:8;"
	WyckoffSyms[115] = "a:0,0,1/4:2;b:1/2,1/2,1/4:2;c:0,0,0:2;d:1/2,1/2,0:2;e:x,x,1/4:4;f:x,x,3/4:4;g:0,0,z:4;h:1/2,1/2,z:4;i:0,1/2,z:4;j:x,y,z:8;"
	WyckoffSyms[116] = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,0:2;d:0,1/2,1/2:2;e:0,0,z:4;f:0,1/2,z:4;g:x,x+1/2,0:4;h:x,x+1/2,1/2:4;i:x,y,z:8;"
	WyckoffSyms[117] = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,1/4:2;d:0,1/2,3/4:2;e:0,0,z:4;f:x,-x+1/2,1/4:4;g:x,x+1/2,1/4:4;h:0,1/2,z:4;i:x,y,z:8;"
	WyckoffSyms[118] = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,1/4:2;d:0,1/2,3/4:2;e:0,0,z:4;f:0,1/2,z:4;g:x,x,0:8;h:x,x+1/2,1/4:8;i:x,0,z:8;j:x,y,z:16;"
	WyckoffSyms[119] = "a:0,0,1/4:4;b:0,0,0:4;c:0,1/2,1/4:4;d:0,1/2,0:4;e:x,x,1/4:8;f:0,0,z:8;g:0,1/2,z:8;h:x,x+1/2,0:8;i:x,y,z:16;"
	WyckoffSyms[120] = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,0:4;d:0,1/2,1/4:4;e:0,0,z:4;f:x,0,0:8;g:x,0,1/2:8;h:0,1/2,z:8;i:x,x,z:8;j:x,y,z:16;"
	WyckoffSyms[121] = "a:0,0,0:4;b:0,0,1/2:4;c:0,0,z:8;d:x,1/4,1/8:8;e:x,y,z:16;"
	WyckoffSyms[122] = "a:0,0,0:1;b:0,0,1/2:1;c:1/2,1/2,0:1;d:1/2,1/2,1/2:1;e:0,1/2,1/2:2;f:0,1/2,0:2;g:0,0,z:2;h:1/2,1/2,z:2;i:0,1/2,z:4;j:x,x,0:4;k:x,x,1/2:4;l:x,0,0:4;m:x,0,1/2:4;n:x,1/2,0:4;o:x,1/2,1/2:4;p:x,y,0:8;q:x,y,1/2:8;r:x,x,z:8;s:x,0,z:8;t:x,1/2,z:8;u:x,y,z:16;"
	WyckoffSyms[123] = "a:0,0,1/4:2;b:0,0,0:2;c:1/2,1/2,1/4:2;d:1/2,1/2,0:2;e:0,1/2,0:4;f:0,1/2,1/4:4;g:0,0,z:4;h:1/2,1/2,z:4;i:0,1/2,z:8;j:x,x,1/4:8;k:x,0,1/4:8;l:x,1/2,1/4:8;m:x,y,0:8;n:x,y,z:16;"
	WyckoffSyms[124] = "a:1/4,1/4,0:2;b:1/4,1/4,1/2:2;c:3/4,1/4,0:2;d:3/4,1/4,1/2:2;e:0,0,0:4;f:0,0,1/2:4;g:1/4,1/4,z:4;h:3/4,1/4,z:4;i:x,x,0:8;j:x,x,1/2:8;k:x,1/4,0:8;l:x,1/4,1/2:8;m:x,-x,z:8;n:x,y,z:16;"
	WyckoffSyms[125] = "a:1/4,1/4,1/4:2;b:1/4,1/4,3/4:2;c:1/4,3/4,3/4:4;d:1/4,3/4,0:4;e:1/4,1/4,z:4;f:0,0,0:8;g:1/4,3/4,z:8;h:x,x,1/4:8;i:x,1/4,1/4:8;j:x,3/4,1/4:8;k:x,y,z:16;"
	WyckoffSyms[126] = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,1/2:2;d:0,1/2,0:2;e:0,0,z:4;f:0,1/2,z:4;g:x,x+1/2,0:4;h:x,x+1/2,1/2:4;i:x,y,0:8;j:x,y,1/2:8;k:x,x+1/2,z:8;l:x,y,z:16;"
	WyckoffSyms[127] = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,0:4;d:0,1/2,1/4:4;e:0,0,z:4;f:0,1/2,z:8;g:x,x+1/2,1/4:8;h:x,y,0:8;i:x,y,z:16;"
	WyckoffSyms[128] = "a:3/4,1/4,0:2;b:3/4,1/4,1/2:2;c:1/4,1/4,z:2;d:0,0,0:4;e:0,0,1/2:4;f:3/4,1/4,z:4;g:x,-x,0:8;h:x,-x,1/2:8;i:1/4,y,z:8;j:x,x,z:8;k:x,y,z:16;"
	WyckoffSyms[129] = "a:3/4,1/4,1/4:4;b:3/4,1/4,0:4;c:1/4,1/4,z:4;d:0,0,0:8;e:3/4,1/4,z:8;f:x,-x,1/4:8;g:x,y,z:16;"
	WyckoffSyms[130] = "a:0,0,0:2;b:1/2,1/2,0:2;c:0,1/2,0:2;d:0,1/2,1/2:2;e:0,0,1/4:2;f:1/2,1/2,1/4:2;g:0,0,z:4;h:1/2,1/2,z:4;i:0,1/2,z:4;j:x,0,0:4;k:x,1/2,1/2:4;l:x,0,1/2:4;m:x,1/2,0:4;n:x,x,1/4:8;o:0,y,z:8;p:1/2,y,z:8;q:x,y,0:8;r:x,y,z:16;"
	WyckoffSyms[131] = "a:0,0,0:2;b:0,0,1/4:2;c:1/2,1/2,0:2;d:1/2,1/2,1/4:2;e:0,1/2,1/4:4;f:0,1/2,0:4;g:0,0,z:4;h:1/2,1/2,z:4;i:x,x,0:4;j:x,x,1/2:4;k:0,1/2,z:8;l:x,0,1/4:8;m:x,1/2,1/4:8;n:x,y,0:8;o:x,x,z:8;p:x,y,z:16;"
	WyckoffSyms[132] = "a:1/4,1/4,0:4;b:3/4,1/4,0:4;c:1/4,1/4,1/4:4;d:3/4,1/4,3/4:4;e:0,0,0:8;f:1/4,1/4,z:8;g:3/4,1/4,z:8;h:x,1/4,0:8;i:x,1/4,1/2:8;j:x,x,1/4:8;k:x,y,z:16;"
	WyckoffSyms[133] = "a:1/4,3/4,1/4:2;b:3/4,1/4,1/4:2;c:1/4,1/4,1/4:4;d:1/4,1/4,0:4;e:0,0,1/2:4;f:0,0,0:4;g:3/4,1/4,z:4;h:1/4,1/4,z:8;i:x,1/4,3/4:8;j:x,1/4,1/4:8;k:x,x,0:8;l:x,x,1/2:8;m:x,-x,z:8;n:x,y,z:16;"
	WyckoffSyms[134] = "a:0,0,0:4;b:0,0,1/4:4;c:0,1/2,0:4;d:0,1/2,1/4:4;e:0,0,z:8;f:0,1/2,z:8;g:x,x+1/2,1/4:8;h:x,y,0:8;i:x,y,z:16;"
	WyckoffSyms[135] = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,0:4;d:0,1/2,1/4:4;e:0,0,z:4;f:x,x,0:4;g:x,-x,0:4;h:0,1/2,z:8;i:x,y,0:8;j:x,x,z:8;k:x,y,z:16;"
	WyckoffSyms[136] = "a:3/4,1/4,3/4:2;b:3/4,1/4,1/4:2;c:3/4,1/4,z:4;d:1/4,1/4,z:4;e:0,0,0:8;f:x,-x,1/4:8;g:1/4,y,z:8;h:x,y,z:16;"
	WyckoffSyms[137] = "a:3/4,1/4,0:4;b:3/4,1/4,3/4:4;c:0,0,1/2:4;d:0,0,0:4;e:1/4,1/4,z:4;f:3/4,1/4,z:8;g:x,-x,1/2:8;h:x,-x,0:8;i:x,x,z:8;j:x,y,z:16;"
	WyckoffSyms[138] = "a:0,0,0:2;b:0,0,1/2:2;c:0,1/2,0:4;d:0,1/2,1/4:4;e:0,0,z:4;f:1/4,1/4,1/4:8;g:0,1/2,z:8;h:x,x,0:8;i:x,0,0:8;j:x,1/2,0:8;k:x,x+1/2,1/4:16;l:x,y,0:16;m:x,x,z:16;n:0,y,z:16;o:x,y,z:32;"
	WyckoffSyms[139] = "a:0,0,1/4:4;b:0,1/2,1/4:4;c:0,0,0:4;d:0,1/2,0:4;e:1/4,1/4,1/4:8;f:0,0,z:8;g:0,1/2,z:8;h:x,x+1/2,0:8;i:x,x,1/4:16;j:x,0,1/4:16;k:x,y,0:16;l:x,x+1/2,z:16;m:x,y,z:32;"
	WyckoffSyms[140] = "a:0,3/4,1/8:4;b:0,1/4,3/8:4;c:0,0,0:8;d:0,0,1/2:8;e:0,1/4,z:8;f:x,0,0:16;g:x,x+1/4,7/8:16;h:0,y,z:16;i:x,y,z:32;"
	WyckoffSyms[141] = "a:0,1/4,3/8:8;b:0,1/4,1/8:8;c:0,0,0:16;d:0,1/4,z:16;e:x,0,1/4:16;f:x,x+1/4,1/8:16;g:x,y,z:32;"
	// Trigonal [143,167]  lines 142-166
	WyckoffSyms[142] = "a:0,0,z:1;b:1/3,2/3,z:1;c:2/3,1/3,z:1;d:x,y,z:3;"
	WyckoffSyms[143] = "a:x,y,z:3;"
	WyckoffSyms[144] = "a:x,y,z:3;"
	WyckoffSyms[145] = "a:0,0,z:3;b:x,y,z:9;"
	WyckoffSyms[146] = "a:0,0,0:1;b:0,0,1/2:1;c:0,0,z:2;d:1/3,2/3,z:2;e:1/2,0,0:3;f:1/2,0,1/2:3;g:x,y,z:6;"
	WyckoffSyms[147] = "a:0,0,0:3;b:0,0,1/2:3;c:0,0,z:6;d:1/2,0,1/2:9;e:1/2,0,0:9;f:x,y,z:18;"
	WyckoffSyms[148] = "a:0,0,0:1;b:0,0,1/2:1;c:1/3,2/3,0:1;d:1/3,2/3,1/2:1;e:2/3,1/3,0:1;f:2/3,1/3,1/2:1;g:0,0,z:2;h:1/3,2/3,z:2;i:2/3,1/3,z:2;j:x,-x,0:3;k:x,-x,1/2:3;l:x,y,z:6;"
	WyckoffSyms[149] = "a:0,0,0:1;b:0,0,1/2:1;c:0,0,z:2;d:1/3,2/3,z:2;e:x,0,0:3;f:x,0,1/2:3;g:x,y,z:6;"
	WyckoffSyms[150] = "a:x,-x,1/3:3;b:x,-x,5/6:3;c:x,y,z:6;"
	WyckoffSyms[151] = "a:x,0,1/3:3;b:x,0,5/6:3;c:x,y,z:6;"
	WyckoffSyms[152] = "a:x,-x,2/3:3;b:x,-x,1/6:3;c:x,y,z:6;"
	WyckoffSyms[153] = "a:x,0,2/3:3;b:x,0,1/6:3;c:x,y,z:6;"
	WyckoffSyms[154] = "a:0,0,0:3;b:0,0,1/2:3;c:0,0,z:6;d:x,0,0:9;e:x,0,1/2:9;f:x,y,z:18;"
	WyckoffSyms[155] = "a:0,0,z:1;b:1/3,2/3,z:1;c:2/3,1/3,z:1;d:x,-x,z:3;e:x,y,z:6;"
	WyckoffSyms[156] = "a:0,0,z:1;b:1/3,2/3,z:2;c:x,0,z:3;d:x,y,z:6;"
	WyckoffSyms[157] = "a:0,0,z:2;b:1/3,2/3,z:2;c:2/3,1/3,z:2;d:x,y,z:6;"
	WyckoffSyms[158] = "a:0,0,z:2;b:1/3,2/3,z:2;c:x,y,z:6;"
	WyckoffSyms[159] = "a:0,0,z:3;b:x,-x,z:9;c:x,y,z:18;"
	WyckoffSyms[160] = "a:0,0,z:6;b:x,y,z:18;"
	WyckoffSyms[161] = "a:0,0,0:1;b:0,0,1/2:1;c:1/3,2/3,0:2;d:1/3,2/3,1/2:2;e:0,0,z:2;f:1/2,0,0:3;g:1/2,0,1/2:3;h:1/3,2/3,z:4;i:x,-x,0:6;j:x,-x,1/2:6;k:x,0,z:6;l:x,y,z:12;"
	WyckoffSyms[162] = "a:0,0,1/4:2;b:0,0,0:2;c:1/3,2/3,1/4:2;d:2/3,1/3,1/4:2;e:0,0,z:4;f:1/3,2/3,z:4;g:1/2,0,0:6;h:x,-x,1/4:6;i:x,y,z:12;"
	WyckoffSyms[163] = "a:0,0,0:1;b:0,0,1/2:1;c:0,0,z:2;d:1/3,2/3,z:2;e:1/2,0,0:3;f:1/2,0,1/2:3;g:x,0,0:6;h:x,0,1/2:6;i:x,-x,z:6;j:x,y,z:12;"
	WyckoffSyms[164] = "a:0,0,1/4:2;b:0,0,0:2;c:0,0,z:4;d:1/3,2/3,z:4;e:1/2,0,0:6;f:x,0,1/4:6;g:x,y,z:12;"
	WyckoffSyms[165] = "a:0,0,0:3;b:0,0,1/2:3;c:0,0,z:6;d:1/2,0,1/2:9;e:1/2,0,0:9;f:x,0,0:18;g:x,0,1/2:18;h:x,-x,z:18;i:x,y,z:36;"
	WyckoffSyms[166] = "a:0,0,1/4:6;b:0,0,0:6;c:0,0,z:12;d:1/2,0,0:18;e:x,0,1/4:18;f:x,y,z:36;"
	// Hexagonal [168,194]  lines 167-193
	WyckoffSyms[167] = "a:0,0,z:1;b:1/3,2/3,z:2;c:1/2,0,z:3;d:x,y,z:6;"
	WyckoffSyms[168] = "a:x,y,z:6;"
	WyckoffSyms[169] = "a:x,y,z:6;"
	WyckoffSyms[170] = "a:0,0,z:3;b:1/2,1/2,z:3;c:x,y,z:6;"
	WyckoffSyms[171] = "a:0,0,z:3;b:1/2,1/2,z:3;c:x,y,z:6;"
	WyckoffSyms[172] = "a:0,0,z:2;b:1/3,2/3,z:2;c:x,y,z:6;"
	WyckoffSyms[173] = "a:0,0,0:1;b:0,0,1/2:1;c:1/3,2/3,0:1;d:1/3,2/3,1/2:1;e:2/3,1/3,0:1;f:2/3,1/3,1/2:1;g:0,0,z:2;h:1/3,2/3,z:2;i:2/3,1/3,z:2;j:x,y,0:3;k:x,y,1/2:3;l:x,y,z:6;"
	WyckoffSyms[174] = "a:0,0,0:1;b:0,0,1/2:1;c:1/3,2/3,0:2;d:1/3,2/3,1/2:2;e:0,0,z:2;f:1/2,0,0:3;g:1/2,0,1/2:3;h:1/3,2/3,z:4;i:1/2,0,z:6;j:x,y,0:6;k:x,y,1/2:6;l:x,y,z:12;"
	WyckoffSyms[175] = "a:0,0,1/4:2;b:0,0,0:2;c:1/3,2/3,1/4:2;d:2/3,1/3,1/4:2;e:0,0,z:4;f:1/3,2/3,z:4;g:1/2,0,0:6;h:x,y,1/4:6;i:x,y,z:12;"
	WyckoffSyms[176] = "a:0,0,0:1;b:0,0,1/2:1;c:1/3,2/3,0:2;d:1/3,2/3,1/2:2;e:0,0,z:2;f:1/2,0,0:3;g:1/2,0,1/2:3;h:1/3,2/3,z:4;i:1/2,0,z:6;j:x,0,0:6;k:x,0,1/2:6;l:x,-x,0:6;m:x,-x,1/2:6;n:x,y,z:12;"
	WyckoffSyms[177] = "a:x,0,0:6;b:x,2x,1/4:6;c:x,y,z:12;"
	WyckoffSyms[178] = "a:x,0,0:6;b:x,2x,3/4:6;c:x,y,z:12;"
	WyckoffSyms[179] = "a:0,0,0:3;b:0,0,1/2:3;c:1/2,0,0:3;d:1/2,0,1/2:3;e:0,0,z:6;f:1/2,0,z:6;g:x,0,0:6;h:x,0,1/2:6;i:x,2x,0:6;j:x,2x,1/2:6;k:x,y,z:12;"
	WyckoffSyms[180] = "a:0,0,0:3;b:0,0,1/2:3;c:1/2,0,0:3;d:1/2,0,1/2:3;e:0,0,z:6;f:1/2,0,z:6;g:x,0,0:6;h:x,0,1/2:6;i:x,2x,0:6;j:x,2x,1/2:6;k:x,y,z:12;"
	WyckoffSyms[181] = "a:0,0,0:2;b:0,0,1/4:2;c:1/3,2/3,1/4:2;d:1/3,2/3,3/4:2;e:0,0,z:4;f:1/3,2/3,z:4;g:x,0,0:6;h:x,2x,1/4:6;i:x,y,z:12;"
	WyckoffSyms[182] = "a:0,0,z:1;b:1/3,2/3,z:2;c:1/2,0,z:3;d:x,0,z:6;e:x,-x,z:6;f:x,y,z:12;"
	WyckoffSyms[183] = "a:0,0,z:2;b:1/3,2/3,z:4;c:1/2,0,z:6;d:x,y,z:12;"
	WyckoffSyms[184] = "a:0,0,z:2;b:1/3,2/3,z:4;c:x,0,z:6;d:x,y,z:12;"
	WyckoffSyms[185] = "a:0,0,z:2;b:1/3,2/3,z:2;c:x,-x,z:6;d:x,y,z:12;"
	WyckoffSyms[186] = "a:0,0,0:1;b:0,0,1/2:1;c:1/3,2/3,0:1;d:1/3,2/3,1/2:1;e:2/3,1/3,0:1;f:2/3,1/3,1/2:1;g:0,0,z:2;h:1/3,2/3,z:2;i:2/3,1/3,z:2;j:x,-x,0:3;k:x,-x,1/2:3;l:x,y,0:6;m:x,y,1/2:6;n:x,-x,z:6;o:x,y,z:12;"
	WyckoffSyms[187] = "a:0,0,0:2;b:0,0,1/4:2;c:1/3,2/3,0:2;d:1/3,2/3,1/4:2;e:2/3,1/3,0:2;f:2/3,1/3,1/4:2;g:0,0,z:4;h:1/3,2/3,z:4;i:2/3,1/3,z:4;j:x,-x,0:6;k:x,y,1/4:6;l:x,y,z:12;"
	WyckoffSyms[188] = "a:0,0,0:1;b:0,0,1/2:1;c:1/3,2/3,0:2;d:1/3,2/3,1/2:2;e:0,0,z:2;f:x,0,0:3;g:x,0,1/2:3;h:1/3,2/3,z:4;i:x,0,z:6;j:x,y,0:6;k:x,y,1/2:6;l:x,y,z:12;"
	WyckoffSyms[189] = "a:0,0,0:2;b:0,0,1/4:2;c:1/3,2/3,1/4:2;d:2/3,1/3,1/4:2;e:0,0,z:4;f:1/3,2/3,z:4;g:x,0,0:6;h:x,y,1/4:6;i:x,y,z:12;"
	WyckoffSyms[190] = "a:0,0,0:1;b:0,0,1/2:1;c:1/3,2/3,0:2;d:1/3,2/3,1/2:2;e:0,0,z:2;f:1/2,0,0:3;g:1/2,0,1/2:3;h:1/3,2/3,z:4;i:1/2,0,z:6;j:x,0,0:6;k:x,0,1/2:6;l:x,2x,0:6;m:x,2x,1/2:6;n:x,0,z:12;o:x,2x,z:12;p:x,y,0:12;q:x,y,1/2:12;r:x,y,z:24;"
	WyckoffSyms[191] = "a:0,0,1/4:2;b:0,0,0:2;c:1/3,2/3,1/4:4;d:1/3,2/3,0:4;e:0,0,z:4;f:1/2,0,1/4:6;g:1/2,0,0:6;h:1/3,2/3,z:8;i:1/2,0,z:12;j:x,0,1/4:12;k:x,2x,1/4:12;l:x,y,0:12;m:x,y,z:24;"
	WyckoffSyms[192] = "a:0,0,1/4:2;b:0,0,0:2;c:1/3,2/3,1/4:4;d:1/3,2/3,0:4;e:0,0,z:4;f:1/2,0,0:6;g:x,0,1/4:6;h:1/3,2/3,z:8;i:x,2x,0:12;j:x,y,1/4:12;k:x,0,z:12;l:x,y,z:24;"
	WyckoffSyms[193] = "a:0,0,0:2;b:0,0,1/4:2;c:1/3,2/3,1/4:2;d:1/3,2/3,3/4:2;e:0,0,z:4;f:1/3,2/3,z:4;g:1/2,0,0:6;h:x,2x,1/4:6;i:x,0,0:12;j:x,y,1/4:12;k:x,2x,z:12;l:x,y,z:24;"
	// Cubic [195,230]  lines 194-229
	WyckoffSyms[194] = "a:0,0,0:1;b:1/2,1/2,1/2:1;c:0,1/2,1/2:3;d:1/2,0,0:3;e:x,x,x:4;f:x,0,0:6;g:x,0,1/2:6;h:x,1/2,0:6;i:x,1/2,1/2:6;j:x,y,z:12;"
	WyckoffSyms[195] = "a:0,0,0:4;b:1/2,1/2,1/2:4;c:1/4,1/4,1/4:4;d:3/4,3/4,3/4:4;e:x,x,x:16;f:x,0,0:24;g:x,1/4,1/4:24;h:x,y,z:48;"
	WyckoffSyms[196] = "a:0,0,0:2;b:0,1/2,1/2:6;c:x,x,x:8;d:x,0,0:12;e:x,1/2,0:12;f:x,y,z:24;"
	WyckoffSyms[197] = "a:x,x,x:4;b:x,y,z:12;"
	WyckoffSyms[198] = "a:x,x,x:8;b:x,0,1/4:12;c:x,y,z:24;"
	WyckoffSyms[199] = "a:0,0,0:1;b:1/2,1/2,1/2:1;c:0,1/2,1/2:3;d:1/2,0,0:3;e:x,0,0:6;f:x,0,1/2:6;g:x,1/2,0:6;h:x,1/2,1/2:6;i:x,x,x:8;j:0,y,z:12;k:1/2,y,z:12;l:x,y,z:24;"
	WyckoffSyms[200] = "a:1/4,1/4,1/4:2;b:0,0,0:4;c:1/2,1/2,1/2:4;d:1/4,3/4,3/4:6;e:x,x,x:8;f:x,1/4,1/4:12;g:x,3/4,1/4:12;h:x,y,z:24;"
	WyckoffSyms[201] = "a:0,0,0:4;b:1/2,1/2,1/2:4;c:1/4,1/4,1/4:8;d:0,1/4,1/4:24;e:x,0,0:24;f:x,x,x:32;g:x,1/4,1/4:48;h:0,y,z:48;i:x,y,z:96;"
	WyckoffSyms[202] = "a:1/8,1/8,1/8:8;b:5/8,5/8,5/8:8;c:0,0,0:16;d:1/2,1/2,1/2:16;e:x,x,x:32;f:x,1/8,1/8:48;g:x,y,z:96;"
	WyckoffSyms[203] = "a:0,0,0:2;b:0,1/2,1/2:6;c:1/4,1/4,1/4:8;d:x,0,0:12;e:x,0,1/2:12;f:x,x,x:16;g:0,y,z:24;h:x,y,z:48;"
	WyckoffSyms[204] = "a:0,0,0:4;b:1/2,1/2,1/2:4;c:x,x,x:8;d:x,y,z:24;"
	WyckoffSyms[205] = "a:0,0,0:8;b:1/4,1/4,1/4:8;c:x,x,x:16;d:x,0,1/4:24;e:x,y,z:48;"
	WyckoffSyms[206] = "a:0,0,0:1;b:1/2,1/2,1/2:1;c:0,1/2,1/2:3;d:1/2,0,0:3;e:x,0,0:6;f:x,1/2,1/2:6;g:x,x,x:8;h:x,1/2,0:12;i:0,y,y:12;j:1/2,y,y:12;k:x,y,z:24;"
	WyckoffSyms[207] = "a:0,0,0:2;b:1/4,1/4,1/4:4;c:3/4,3/4,3/4:4;d:0,1/2,1/2:6;e:1/4,0,1/2:6;f:1/4,1/2,0:6;g:x,x,x:8;h:x,0,0:12;i:x,0,1/2:12;j:x,1/2,0:12;k:1/4,y,-y+1/2:12;l:1/4,y,y+1/2:12;m:x,y,z:24;"
	WyckoffSyms[208] = "a:0,0,0:4;b:1/2,1/2,1/2:4;c:1/4,1/4,1/4:8;d:0,1/4,1/4:24;e:x,0,0:24;f:x,x,x:32;g:0,y,y:48;h:1/2,y,y:48;i:x,1/4,1/4:48;j:x,y,z:96;"
	WyckoffSyms[209] = "a:0,0,0:8;b:1/2,1/2,1/2:8;c:1/8,1/8,1/8:16;d:5/8,5/8,5/8:16;e:x,x,x:32;f:x,0,0:48;g:1/8,y,-y+1/4:48;h:x,y,z:96;"
	WyckoffSyms[210] = "a:0,0,0:2;b:0,1/2,1/2:6;c:1/4,1/4,1/4:8;d:1/4,1/2,0:12;e:x,0,0:12;f:x,x,x:16;g:x,1/2,0:24;h:0,y,y:24;i:1/4,y,-y+1/2:24;j:x,y,z:48;"
	WyckoffSyms[211] = "a:1/8,1/8,1/8:4;b:5/8,5/8,5/8:4;c:x,x,x:8;d:1/8,y,-y+1/4:12;e:x,y,z:24;"
	WyckoffSyms[212] = "a:3/8,3/8,3/8:4;b:7/8,7/8,7/8:4;c:x,x,x:8;d:1/8,y,y+1/4:12;e:x,y,z:24;"
	WyckoffSyms[213] = "a:1/8,1/8,1/8:8;b:7/8,7/8,7/8:8;c:1/8,0,1/4:12;d:5/8,0,1/4:12;e:x,x,x:16;f:x,0,1/4:24;g:1/8,y,y+1/4:24;h:1/8,y,-y+1/4:24;i:x,y,z:48;"
	WyckoffSyms[214] = "a:0,0,0:1;b:1/2,1/2,1/2:1;c:0,1/2,1/2:3;d:1/2,0,0:3;e:x,x,x:4;f:x,0,0:6;g:x,1/2,1/2:6;h:x,1/2,0:12;i:x,x,z:12;j:x,y,z:24;"
	WyckoffSyms[215] = "a:0,0,0:4;b:1/2,1/2,1/2:4;c:1/4,1/4,1/4:4;d:3/4,3/4,3/4:4;e:x,x,x:16;f:x,0,0:24;g:x,1/4,1/4:24;h:x,x,z:48;i:x,y,z:96;"
	WyckoffSyms[216] = "a:0,0,0:2;b:0,1/2,1/2:6;c:x,x,x:8;d:1/4,1/2,0:12;e:x,0,0:12;f:x,1/2,0:24;g:x,x,z:24;h:x,y,z:48;"
	WyckoffSyms[217] = "a:0,0,0:2;b:0,1/2,1/2:6;c:1/4,1/2,0:6;d:1/4,0,1/2:6;e:x,x,x:8;f:x,0,0:12;g:x,1/2,0:12;h:x,0,1/2:12;i:x,y,z:24;"
	WyckoffSyms[218] = "a:0,0,0:8;b:1/4,1/4,1/4:8;c:0,1/4,1/4:24;d:1/4,0,0:24;e:x,x,x:32;f:x,0,0:48;g:x,1/4,1/4:48;h:x,y,z:96;"
	WyckoffSyms[219] = "a:3/8,0,1/4:12;b:7/8,0,1/4:12;c:x,x,x:16;d:x,0,1/4:24;e:x,y,z:48;"
	WyckoffSyms[220] = "a:0,0,0:1;b:1/2,1/2,1/2:1;c:0,1/2,1/2:3;d:1/2,0,0:3;e:x,0,0:6;f:x,1/2,1/2:6;g:x,x,x:8;h:x,1/2,0:12;i:0,y,y:12;j:1/2,y,y:12;k:0,y,z:24;l:1/2,y,z:24;m:x,x,z:24;n:x,y,z:48;"
	WyckoffSyms[221] = "a:1/4,1/4,1/4:2;b:3/4,1/4,1/4:6;c:0,0,0:8;d:0,3/4,1/4:12;e:x,1/4,1/4:12;f:x,x,x:16;g:x,3/4,1/4:24;h:1/4,y,y:24;i:x,y,z:48;"
	WyckoffSyms[222] = "a:0,0,0:2;b:0,1/2,1/2:6;c:1/4,0,1/2:6;d:1/4,1/2,0:6;e:1/4,1/4,1/4:8;f:x,0,0:12;g:x,0,1/2:12;h:x,1/2,0:12;i:x,x,x:16;j:1/4,y,y+1/2:24;k:0,y,z:24;l:x,y,z:48;"
	WyckoffSyms[223] = "a:1/4,1/4,1/4:2;b:0,0,0:4;c:1/2,1/2,1/2:4;d:1/4,3/4,3/4:6;e:x,x,x:8;f:1/2,1/4,3/4:12;g:x,1/4,1/4:12;h:x,1/4,3/4:24;i:1/2,y,y+1/2:24;j:1/2,y,-y:24;k:x,x,z:24;l:x,y,z:48;"
	WyckoffSyms[224] = "a:0,0,0:4;b:1/2,1/2,1/2:4;c:1/4,1/4,1/4:8;d:0,1/4,1/4:24;e:x,0,0:24;f:x,x,x:32;g:x,1/4,1/4:48;h:0,y,y:48;i:1/2,y,y:48;j:0,y,z:96;k:x,x,z:96;l:x,y,z:192;"
	WyckoffSyms[225] = "a:1/4,1/4,1/4:8;b:0,0,0:8;c:1/4,0,0:24;d:0,1/4,1/4:24;e:x,0,0:48;f:x,1/4,1/4:48;g:x,x,x:64;h:1/4,y,y:96;i:0,y,z:96;j:x,y,z:192;"
	WyckoffSyms[226] = "a:1/8,1/8,1/8:8;b:3/8,3/8,3/8:8;c:0,0,0:16;d:1/2,1/2,1/2:16;e:x,x,x:32;f:x,1/8,1/8:48;g:x,x,z:96;h:0,y,-y:96;i:x,y,z:192;"
	WyckoffSyms[227] = "a:1/8,1/8,1/8:16;b:1/4,1/4,1/4:32;c:0,0,0:32;d:7/8,1/8,1/8:48;e:x,x,x:64;f:x,1/8,1/8:96;g:1/4,y,-y:96;h:x,y,z:192;"
	WyckoffSyms[228] = "a:0,0,0:2;b:0,1/2,1/2:6;c:1/4,1/4,1/4:8;d:1/4,0,1/2:12;e:x,0,0:12;f:x,x,x:16;g:x,0,1/2:24;h:0,y,y:24;i:1/4,y,-y+1/2:48;j:0,y,z:48;k:x,x,z:48;l:x,y,z:96;"
	WyckoffSyms[229] = "a:0,0,0:16;b:1/8,1/8,1/8:16;c:1/8,0,1/4:24;d:3/8,0,1/4:24;e:x,x,x:32;f:x,0,1/4:48;g:1/8,y,-y+1/4:48;h:x,y,z:96;"

	String list=WyckoffSyms[SG-1]
	Make/N=(ItemsInList(list),3)/FREE/T WyckList
	WyckList[][0] = StringFromList(0,StringFromList(p,list),":")
	WyckList[][1] = StringFromList(1,StringFromList(p,list),":")
	WyckList[][2] = StringFromList(2,StringFromList(p,list),":")
	return WyckList
End
//	Function test_GetWyckoffSymStrings(SG)
//		Variable SG
//		Wave WyckList=GetWyckoffSymStrings(SG)
//		print WyckList
//	End
