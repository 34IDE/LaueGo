#pragma rtGlobals=1		// Use modern global access method.
#pragma IgorVersion = 4.0
#pragma version = 2.24
#include "Elements", version>=2.00
#include "CromerLiberman", version>=1.9

//	Sept 14, 2005
//		also changed Get_MuFormula() and MuOverRhoFormula() so they will not prompt if called by another function
//	Apr 30, 2014
//		added the USE_OBSOLETE, added some comments, changed line termination to CR -> LF
//	June 4, 2014
//		updated a few things
//	May 3, 2016
//		added the calls to ElementDataInitPackage()


//	#define USE_OBSOLETE		// use this define to use the obsolete functions at the end of this procedure file


Menu "Analysis"
	Submenu "X-ray"
		"mu of material",Get_MuFormula("",NaN)
		"mu \  rho", MuOverRhoFormula("",NaN)
	End
End



Function Get_MuFormula(formula,keV)			// get mu (1/µm)
	// includes the data base of densities
	String formula		// chemical formula, Elements start with a capital letter.
	Variable keV		// photon energy (keV)
	Variable topLevel = ItemsInList(GetRTStackInfo(0))<2
	if (strlen(formula)<1 || numtype(keV) || keV<1e-3)
		if (!topLevel)
			return NaN
		endif
		Prompt formula, "chemical formula"
		Prompt keV, "energy (keV)"
		DoPrompt "formula and energy",formula,keV
		if (V_flag)
			return NaN
		endif
	endif
	String formula_In = formula
	if (strlen(formula)<1 || numtype(keV) || keV<1e-3)
		return NaN
	endif
	if (Exists("root:Packages:Elements:name")!=1)
		ElementDataInitPackage()
	endif

	Variable density
	String list = ChemicalFormulaToList(formula)	// parse the chemical formula into a list
	if (ItemsInList(list)==1)				// a single element, look it up
		density = Element_density(element2Z(StringFromList(0,list)))
	else										// check in the density list
		String densityList=StrVarOrDefault("root:Packages:Elements:densityList","TiO2:4.26;NiO:6.67;CoO:6.45;Fe2O3:5.24;")
		density = str2num(StringByKey(formula,densityList))
	endif
	if (numtype(density))
		list = LookUpMaterial(formula)
		density = NumberByKey("density",list,"=")
		formula = StringByKey("formula",list,"=")
	endif

	Variable mu = MuOverRhoFormula(formula,keV)*density/1e4
	if (topLevel)
		if (numtype(mu))
			printf  "   unable to compute mu of '%s' at %g keV\r",formula_In,keV
		else
			printf  "   mu of '%s' at %g keV  is %g (1/micron)  -->  absorption length %g (micron)\r",formula_In,keV,mu,1/mu
		endif
	endif
	return mu
End


Function MuOverRhoFormula(formula,keV)	// returns mu/rho  (cm^2/g) the given material
	String formula						// chemical formula, Elements start with a capital letter.
	Variable keV						// photon energy (keV)
	Variable topLevel = ItemsInList(GetRTStackInfo(0))<2
	if (strlen(formula)<1 || numtype(keV) || keV<1e-3)
		if (!topLevel)
			return NaN
		endif
		Prompt formula, "chemical formula"
		Prompt keV, "energy (keV)"
		DoPrompt "formula and energy",formula,keV
		if (V_flag)
			return NaN
		endif
	endif
	String formula_In = formula
	if (strlen(formula)<1 || numtype(keV) || keV<1e-3)
		return NaN
	endif
	if (Exists("root:Packages:Elements:name")!=1)
		ElementDataInitPackage()
	endif

	String symb						// element symbol
	String item							// one element with its number of atoms
	Variable molWt						// molecular weight
	Variable i=0
	formula = ChemicalFormulaToList(formula)	// parse the chemical formula into a list
	molWt = Formula2amu(formula)	// molecular weight of formula (amu)
	Variable mu_rho=0					// mu/rho for the molecule  (cm^2/g)
	for (i=0;strlen(StringFromList(i,formula));i+=1)
		item = StringFromList(i,formula)
		symb = StringFromList(0,item,",")
		mu_rho += Get_MuOverRho(symb,keV)*Formula2amu(item)		// accumulate amu * (mu/rho)
	endfor
	if (topLevel)
		printf  "   mu/rho of '%s' at %g keV  is %g (cm^2/g)\r",formula_In,keV,mu_rho/molWt
	endif
	return mu_rho/molWt
End


#ifdef USE_OBSOLETE

// ===========================================================================================
// ===========================================================================================
//		Everything below here is DEPRECATED and OBSOLETE, do NOT use the following functions
// ===========================================================================================
// ===========================================================================================


//  Obsolete function only for compatibility
Function GetfoXray(symb,K)	// compute fo(Q)
	String symb				// name of atom desired
	Variable K					// 2*pi/d	 (1/Angstroms)

	Variable fo = Get_f0(symb,K)
	if (numtype(fo))
		printf "  fo(%s,  K=%g) = Unknown symbol\r",symb,K
	else
		printf "  fo(%s,  K=%g) = %g  (electrons)\r",symb,K, fo
	endif
	return fo
End


//  Obsolete function only for compatibility
Function foXray(symb,K)		// compute fo(Q)
	String symb				// name of atom desired
	Variable K					// 2*pi/d	 (1/Angstroms)

	return Get_f0(symb,K)
End


//	DEPRECATED	DEPRECATED	DEPRECATED	DEPRECATED
// THIS ROUTINE IS DEPRECATED, use the set lattice panel instead.  This is ONLY called by setLattice()
//
// gets lattice constants for known structures via dialogs,  it uses structure numbers from the International Tables
// implemented structures are:
//	FCC				225
//	BCC				229
//	Diamond		227
//	Simple Cubic	221
//	Hexagonal		194
//	Sapphire		167
//	Triclinic		1
Function getLatticeConstants(structureNum,structureName,aa,bb,cc,aalpha,bbeta,ggam)
	Variable &structureNum							// same as for Internationl tables
	String &structureName								// see structures below for list of valid names
	Variable &aa,&bb,&cc,&aalpha,&bbeta,&ggam		// lattice constants, Angstroms and angles in degrees
	aa = numtype(aa)||aa<=0 ? 4.05 : aa				// default is aluminum
	bb = numtype(bb)||bb<=0 ? aa : bb
	cc = numtype(cc)||cc<=0 ? aa : cc
	aalpha = numtype(aalpha)||aalpha<=0 ? 90 : aalpha	// for invalid angles, use 90¡
	bbeta = numtype(bbeta)||bbeta<=0 ? 90 : bbeta
	ggam = numtype(ggam)||ggam<=0 ? 90 : ggam
	//	Cubic			[195,230]		//	a
	//	Hexagonal		[168,194]		//	a,c
	//	Trigonal		[143,167]		//	a,alpha
	//	Tetragonal		[75,142]		//	a,c
	//	Orthorhombic	[16,74]		//	a,b,c
	//	Monoclinic		[3,15]			//	a,b,c,gamma
	//	Triclinic		[1,2]			//	a,b,c,alpha,beta,gamma
	String item="",structures="FCC:225;BCC:229;Diamond:227;Simple Cubic:221;Hexagonal:194;Sapphire:167;Triclinic:1"
	Variable i,N=ItemsInList(structures)
	if (structureNum>=1 && structureNum<=230)	// valid structureNum
		structureName = ""
		for(i=0;i<N;i+=1)
			item = StringFromList(i,structures)
			if (structureNum==str2num(StringFromList(1,item,":")))
				structureName = StringFromList(0,item,":")
			endif
		endfor
		if (strlen(structureName)<1)
			DoAlert 0, "unknown structure number "+num2str(structureNum)
		endif
	endif
	if (strlen(structureName)>1)					// a structureName passed
		structureNum = -1
		for(i=0;i<N;i+=1)
			item = StringFromList(i,structures)
			if (stringmatch(structureName,StringFromList(0,item,":")))
				structureNum = str2num(StringFromList(1,item,":"))
			endif
		endfor
		if (structureNum<1)
			DoAlert 0, "unknown structure name '"+structureName+"'"
		endif
	endif
	item = structureName+":"+num2istr(structureNum)
	Prompt item, "lattice structure", popup, structures
	DoPrompt "lattice structure",item
	if (V_flag)
		return 1
	endif
	structureNum=str2num(StringFromList(1,item,":"))
	structureNum = max(min(structureNum,230),1)				// force into range [1,230]
	structureName = StringFromList(0,item,":")

	Variable a=aa,b=bb,c=cc,alpha=aalpha,beta=bbeta,gam=ggam	// local values of lattice constants
	Prompt a, "a ()"
	Prompt b, "b ()"
	Prompt c, "c ()"
	Prompt alpha, "alpha (¡)"
	Prompt beta, "beta (¡)"
	Prompt gam, "gamma (¡)"

	if (structureNum>=195)		// Cubic
		DoPrompt "Cubic lattice constants",a
		b=a;  c=a
		alpha=90;  beta=90;  gam=90
	elseif(structureNum>=168)	// Hexagonal
		c = (a==c) ? 3*a : c
		DoPrompt "Hexagonal lattice constants",a,c
		b = a
		alpha=90;  beta=90;  gam=120
	elseif(structureNum>=143)	// Trigonal, use the hexagonal cell
		c = (a==c) ? 3*a : c
		DoPrompt "Trigonal lattice constants",a,c
		b = a
		alpha=90;  beta=90;  gam=120
	elseif(structureNum>=75)		// Tetragonal
		DoPrompt "Tetragonal lattice constants",a,c
		b = a
		alpha=90;  beta=90;  gam=90
	elseif(structureNum>=16)		// Orthorhombic
		DoPrompt "Orthorhombic lattice constants",a,b,c
		alpha=90;  beta=90;  gam=90
	elseif(structureNum>=3)		// Monoclinic
		DoPrompt "Monoclinic lattice constants",a,b,c,gam
		alpha=90;  beta=90
	else								// Triclinic
		DoPrompt "Triclinic lattice constants",a,b,c,alpha,beta,gam
	endif
	ForceLatticeToStructure(structureNum,a,b,c,alpha,beta,gam)
	aa = a;			bb = b;			cc = c
	aalpha = alpha;	bbeta = beta;	ggam = gam
	return 0
End

// forces lattice constants to match the structure number (e.g. for cubic, forces b and c to be a, and all angles 90)
Static Function ForceLatticeToStructure(structureNum,a,b,c,alpha,beta,gam)
	Variable structureNum								// same as for Internationl tables
	Variable a,&b,&c,&alpha,&beta,&gam				// lattice constants (note, a is never changed), angles in degrees

	if (numtype(structureNum) || structureNum<0 || structureNum>230)	// invalid structureNum, it must be in [1,230]
		DoAlert 0, "invalid structure number "+num2str(structureNum)
		return 1
	endif
	//	Cubic			[195,230]		//	a
	//	Hexagonal		[168,194]		//	a,c
	//	Trigonal		[143,167]		//	a,alpha
	//	Tetragonal		[75,142]		//	a,c
	//	Orthorhombic	[16,74]		//	a,b,c
	//	Monoclinic		[3,15]			//	a,b,c,gamma
	//	Triclinic		[1,2]			//	a,b,c,alpha,beta,gamma

	if (structureNum>=195)		// Cubic
		b = a
		c = a
		alpha=90;  beta=90;  gam=90
	elseif(structureNum>=168)	// Hexagonal
		b = a
		alpha=90;  beta=90;  gam=120
	elseif(structureNum>=143)	// Trigonal, use the hexagonal cell
		b = a
		alpha=90;  beta=90;  gam=120
	elseif(structureNum>=75)		// Tetragonal
		b = a
		alpha=90;  beta=90;  gam=90
	elseif(structureNum>=16)		// Orthorhombic
		alpha=90;  beta=90;  gam=90
	elseif(structureNum>=3)		// Monoclinic
		alpha=90;  beta=90
//	else								// Triclinic
	endif

	String str
	if (a<=0 || b<=0 || c<=0 || numtype(a+b+c))			// check for valid a,b,c
		sprintf str,"invalid, (a,b,c) = (%g,%g,%g)",a,b,c
		DoAlert 0, str
		return 1
	endif
	if (alpha<=0 || beta<=0 || gam<=0 || alpha>=180 || beta>=180 || gam>=180 || numtype(alpha+beta+gam))	// check for valid angles
		sprintf str,"invalid, (alpha,beta,gam) = (%g,%g,%g)",alpha,beta,gam
		DoAlert 0, str
		return 1
	endif
	return 0
End


#endif
