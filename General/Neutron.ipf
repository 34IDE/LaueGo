#pragma rtGlobals=1		// Use modern global access method.
#pragma IgorVersion = 4.0
#pragma version = 1.3
#include  "Elements", version>=1.0
#include "PhysicalConstants", version>=2.02


Menu "Analysis"
	Submenu "Neutron"
		"Neutron Absorption", muNeutron("",NaN,NaN)
		"Neutron Data", NeutronData("","")
		"-"
		"velocity -> eV", velocity2eV(NaN)
		"velocity -> lambda", velocity2lambda(NaN)
		"eV -> velocity", eV2velocity(NaN)
		"eV -> lambda", eV2lambda(NaN)
		"wavelength -> eV", lambda2eV(NaN)
		"wavelength -> velocity", lambda2velocity(NaN)
	End
End

Static Constant NA			= 6.02214076e+23		// Avagadro's number
Static Constant mn_kg		= 1.67492749804e-27	// mass of neutron (kg)
Static Constant e_C		= 1.602176634e-19		// elementary Charge (C)
Static Constant c_ms		= 299792458				// speed of light (m/s) (exact)
Static Constant h_Js		= 6.62607015e-34		// Planck constant (J s)



//				Unit		Quantity
//Isotope					Isotope
//Zatom					Z
//mass						atomic mass number for individual isotopes
//conc						Natural abundance
//Coh_b		fm			bound coherent scattering length
//Inc_b		fm			bound incoherent scattering length
//Coh_xs		barn		bound coherent scattering cross section
//Inc_xs		barn		bound incoherent scattering cross section
//Scatt_xs		barn		total bound scattering cross section
//Abs_xs		barn		absorption cross section for 2200 m/s neutrons


Function muNeutron(symb,lambda,rho)		// returns mu (1/mm)
	String symb
	Variable lambda							// wavelength (�)
	Variable rho							// density (g/cc), if NaN or <0, use value for element

	Variable iprint=0
	Variable jN = Isotope2index(symb)
	Variable lambda0 = velocity2lambda(2200)
	if (numtype(jN) || numtype(rho) || numtype(lambda) || lambda<0)
		WAVE/T isotopes=root:Packages:Neutron:Isotope
		String isoList=""
		Variable i
		for (i=0;i<numpnts(isotopes);i+=1)
			isoList = AddListItem(isotopes[i],isoList,";",inf)
		endfor
		rho = (numtype(rho)|| rho<0) ? Element_density(Isotope2Z(symb)) : rho
		lambda = (numtype(lambda) || lambda<0) ? lambda0 : lambda
		Prompt rho, "density of atom (g/cc)"
		Prompt symb, "element or isotope", popup, isoList
		Prompt lambda, "neutron wavelength (�)"
		DoPrompt "select", symb,lambda,rho
		if (V_flag)
			return NaN
		endif
		iprint = 1
		jN = Isotope2index(symb)
	endif
	Variable Z=Isotope2Z(symb)
	if (numtype(Z))
		return NaN
	endif

	rho = (numtype(rho)|| rho<0) ? Element_density(Z) : rho
	Variable amu = Element_amu(Z)
//	Variable NA=6.02214199e23			// Avagadro's number
	Variable mu								// mu in (1/cm)
	Wave Abs_xs = root:Packages:Neutron:Abs_xs
	Variable sigma = Abs_xs[jN]*1e-24

	mu = rho*sigma*NA/amu				//	mu = rho*sigma*NA/amu  (1/cm)
	mu *= lambda/lambda0
	if(iprint)
		printf "for '%s' at %g (g/cc),  %g (�),   �=%g (1/cm),   absorption length=%g (mm)\r",symb,rho,lambda,mu,10/mu
	endif
	return mu/10							// absorption factor  (1/mm)
End


Function/C NeutronData(symb,type)		// returns neutron data
	String symb
	String type
	Variable rho							// density (g/cc), if NaN or <0, use value for element

	String promptList, typeList = "Coh_b;Inc_b;Coh_xs;Inc_xs;Scatt_xs;Abs_xs;conc;amu"
	Variable itype=WhichListItem(LowerStr(type),LowerStr(typeList))

	Variable iprint=0
	Variable jN = Isotope2index(symb)
	if (numtype(jN) || numtype(jN) || numtype(itype) || itype<0 || itype>7)
		WAVE/T isotopes=root:Packages:Neutron:Isotope
		String isoList=""
		Variable i
		for (i=0;i<numpnts(isotopes);i+=1)
			isoList = AddListItem(isotopes[i],isoList,";",inf)
		endfor
		promptList = "bound coherent scattering length;bound incoherent scattering length;bound coherent scattering cross section;"
		promptList += "bound incoherent scattering cross section;total bound scattering cross section;absorption cross section;concentration;atomic mass"
		Prompt symb, "element or isotope", popup, isoList
		Prompt itype,"information",popup,promptList
		itype += 1
		DoPrompt "select", symb,itype
		if (V_flag)
			return NaN
		endif
		itype -= 1
		iprint = 1
		jN = Isotope2index(symb)
	endif

	if (numtype(jN) || numtype(itype) || itype<0 || itype>7)
		return NaN
	endif

	String line=""
	Variable/C cval
	switch(itype)
		case 0:						// Coh_b
			Wave/C Coh_b=root:Packages:Neutron:Coh_b
			cval = Coh_b[jN]
			sprintf line, "bound coherent scattering length = (%g,%g) (fm)",real(cval),imag(cval)
			break
		case 1:						// Inc_b
			Wave/C Inc_b=root:Packages:Neutron:Inc_b
			cval = Inc_b[jN]
			sprintf line, "bound incoherent scattering length = (%g,%g) (fm)",real(cval),imag(cval)
			break
		case 2:						// Coh_xs
			Wave Coh_xs=root:Packages:Neutron:Coh_xs
			cval = cmplx(Coh_xs[jN],0)
			sprintf line, "bound coherent scattering cross section = %g (barn)",real(cval)
			break
		case 3:						// Inc_xs
			Wave Inc_xs=root:Packages:Neutron:Inc_xs
			cval = cmplx(Inc_xs[jN],0)
			sprintf line, "bound incoherent scattering cross section = %g (barn)",real(cval)
			break
		case 4:						// Scatt_xs
			Wave Scatt_xs=root:Packages:Neutron:Scatt_xs
			cval = cmplx(Scatt_xs[jN],0)
			sprintf line, "total bound scattering cross section = %g (barn)",real(cval)
			break
		case 5:						// Abs_xs
			Wave Abs_xs=root:Packages:Neutron:Abs_xs
			cval = cmplx(Abs_xs[jN],0)
			sprintf line, "absorption cross section (for 2200 m/s neutrons) = %g (barn)",real(cval)
			break
		case 6:						// concentration 
			Wave conc=root:Packages:Neutron:conc
			cval = cmplx(conc[jN],0)
			sprintf line, "Natural abundance = %g %%",real(cval)
			cval /= 100
			break
		case 7:						// atomic mass (amu)
			Wave amuIsotope=root:Packages:Neutron:amuIsotope
			cval = cmplx(amuIsotope[jN],0)
			sprintf line, "atomic mass = %g (amu)",real(cval)
			break
		default:
			sprintf line, "Illegal values given"
			cval = cmplx(NaN,NaN)
	endswitch

	if (iprint)
		print line
	endif
	return cval
End



Function Isotope2index(iso)
	String iso

	Wave/T isotope=root:Packages:Neutron:Isotope
	Variable N=numpnts(isotope),i
	for(i=0;i<N;i+=1)
		if (!cmpstr(iso,isotope[i]))
			break
		endif
	endfor
	return (i<N) ? i : NaN
End



Function Isotope2Z(iso)			// given an isotope (like 29Si) return the Z (=14 for Si)
	String iso

	String ele
	Variable m,z,j
	m = str2num(iso)
	if (numtype(m))			// no mass number
		ele = Iso
	else
		j = floor(log(m))+1
		ele = iso[j,10]
	endif
	z = element2Z(ele)
	return Z
End



Function NeutronDataInitPackage()
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:Neutron

	ElementDataInitPackage()
//	Execute "doIncludePhysicalConstants()"

	String/G root:Packages:Neutron:symbols
	SVAR symbols=root:Packages:Neutron:symbols
	symbols   = "H;He;Li;Be;B;C;N;O;F;Ne;Na;Mg;Al;Si;P;S;Cl;Ar;"
	symbols += "K;Ca;Sc;Ti;V;Cr;Mn;Fe;Co;Ni;Cu;Zn;Ga;Ge;As;Se;Br;Kr;"
	symbols += "Rb;Sr;Y;Zr;Nb;Mo;Tc;Ru;Rh;Pd;Ag;Cd;In;Sn;Sb;Te;I;Xe;"
	symbols += "Cs;Ba;La;Ce;Pr;Nd;Pm;Sm;Eu;Gd;Tb;Dy;Ho;Er;Tm;Yb;Lu;"
	symbols += "Hf;Ta;W;Re;Os;Ir;Pt;Au;Hg;Tl;Pb;Bi;Po;At;Rn;Fr;Ra;Ac;Th;Pa;U;"
	symbols += "Np;Pu;Am;Cm;Bk;Cf;Es;Fm;Md;No;Lr;Rf;Db;Sg;Bh;Hs;Mt"

	MakeNeutronWaves()
End
//Proc doIncludePhysicalConstants()
//	Execute/P "INSERTINCLUDE  \"PhysicalConstants\""
//	Execute/P "COMPILEPROCEDURES "
//	Execute/P "PhysicalConstantsInitPackage()"
//	Execute/P "DELETEINCLUDE  \"PhysicalConstants\""
//	Execute/P "COMPILEPROCEDURES "
//EndMacro



static Function element2Z(symb)			// returns Z for an atomic symbol (NOT case sensitive)
	String symb							// atomic symbol
	SVAR symbols=root:Packages:Neutron:symbols
	symb[0,0] = UpperStr(symb[0,0])		// ensure first char is upper case
	symb[1,1] = LowerStr(symb[1,1])		// and second character is lower
	Variable iz = WhichListItem(symb,symbols)+1
	return ((iz>0) ? iz : NaN)
End





// units of interest are:  velocity (m/s), energy (eV), wavelength (�)

Function velocity2eV(v)				// kinetic energy (eV) for a neutron
	Variable v							// neutron velocity (m/s)
	Variable iprint=0
	if (numtype(v) || v<0)
		v = 2200						// default value
		Prompt v, "neutron velocity (m/s)"
		DoPrompt "velocity (m/s)", v
		if (V_flag)
			return NaN
		endif
		iprint=1
	endif

	Variable m = mn_kg
	Variable eC = e_C					// electron charge (C)
	Variable c = c_ms
	Variable K							// kinetic energy J initiall, later changed to eV
	Variable p = m*v/sqrt(1-(v/c)^2)	// no Newtonian calc, but done for no roundoff at 0
	Variable x2 = (p/(m*c))^2
	K = p^2/m / (1+sqrt(x2+1))
	K /= eC								// change from J to eV
	if (iprint)
		printf "for a neutron going %g (m/s), the kinetic energy is %g (eV)\r", v,K
	endif
	return K							// kinetic energy in eV
End

Function velocity2lambda(v)			// wavelength for a neutron (Angstron)
	Variable v							// neutron velocity (m/s)
	Variable iprint=0
	if (numtype(v) || v<0)
		v = 2200						// default value
		Prompt v, "neutron velocity (m/s)"
		DoPrompt "velocity (m/s)", v
		if (V_flag)
			return NaN
		endif
		iprint=1
	endif
	// p = h/lambda   -->  lambda = h/p
	Variable m = mn_kg
	Variable c = c_ms
	Variable h = h_Js
	Variable p							// momentum
	p = m*v/sqrt(1-(v/c)^2)			// relativity works OK, even for small v
	Variable lambda = h/p * 1e10		// wavelength (Angstrom)
	if (iprint)
		printf "for a neutron going %g (m/s), the wavelength is %g (�)\r", v,lambda
	endif
	return lambda
End


Function eV2velocity(eV)
	Variable eV							// neutron energy (eV)
	Variable iprint=0
	if (numtype(eV) || eV<0)
		eV = .0253						// default value
		Prompt eV, "neutron kinetic energy (eV)"
		DoPrompt "energy (eV)", eV
		if (V_flag)
			return NaN
		endif
		iprint=1
	endif
	Variable m = mn_kg
	Variable c = c_ms
	Variable eC = e_C					// electron charge (C)
	Variable K = eV*eC					// kinetic energy (J)
	Variable kmc2 = K/(m*c^2)
	// K = mc^2 * (gamma - 1)
	Variable v = sqrt(K/m*(2+kmc2)/(1+kmc2)^2)
	if (iprint)
		printf "for a neutron with kinetic energy %g (eV), the velocity is %g (m/s)\r", eV,v
	endif
	return v
End


Function eV2lambda(eV)
	Variable eV							// neutron energy (eV)
	Variable iprint=0
	if (numtype(eV) || eV<0)
		eV = .0253						// default value
		Prompt eV, "neutron kinetic energy (eV)"
		DoPrompt "energy (eV)", eV
		if (V_flag)
			return NaN
		endif
		iprint=1
	endif
	Variable c = c_ms
	Variable m = mn_kg
	Variable eC = e_C					// electron charge (C)
	Variable h = h_Js
	Variable K = eV*eC					// kinetic energy (J)
	Variable p							// momentum

	// (mc^2 + K)^2 = (pc)^2 + (mc2)^2
	// p = sqrt[2mK + (K/c)^2],   this works for small a large p
	p = sqrt(2*m*K + (K/c)^2)
	Variable lambda = h/p * 1e10		// wavelength (�)
	if (iprint)
		printf "for a neutron with kinetic energy %g (eV), the wavelength is %g (�)\r", eV,lambda
	endif
	return lambda
End


Function lambda2eV(lambda)			// neutron energy (eV) from the wavelength
	Variable lambda						// wavelength � on entry, change to (m)
	Variable iprint=0
	if (numtype(lambda) || lambda<0)
		lambda = 1.8					// default value
		Prompt lambda, "neutron wavelength (�)"
		DoPrompt "wavelength (�)", lambda
		if (V_flag)
			return NaN
		endif
		iprint=1
	endif
	lambda *= 1e-10
	Variable c = c_ms
	Variable h = h_Js
	Variable m = mn_kg
	Variable eC = e_C					// electron charge (C)
	Variable p = h/lambda
	Variable K							// kinetic energy
	Variable x2 = (p/m/c)^2
	K = p^2/m / (1+sqrt(x2+1))
	Variable eV = K/eC
	if (iprint)
		printf "for a neutron with wavelength %g (�), the kinetic energy is %g (eV)\r", lambda*1e10,eV
	endif
	return eV
End


Function lambda2velocity(lambda)		// neutron velocity from the wavelength
	Variable lambda						// wavelength � on entry, change to (m)
	Variable iprint=0
	if (numtype(lambda) || lambda<0)
		lambda = 1.8					// default value
		Prompt lambda, "neutron wavelength (�)"
		DoPrompt "wavelength (�)", lambda
		if (V_flag)
			return NaN
		endif
		iprint=1
	endif
	lambda *= 1e-10
	Variable c = c_ms
	Variable h = h_Js
	Variable m = mn_kg
	Variable p = h/lambda
	Variable x2 = (p/m/c)^2
	Variable v = p/m / sqrt(1+x2)		// neutron velocity (m/s)
	if (iprint)
		printf "for a neutron with wavelength %g (�), the velocity is %g (m/s)\r", lambda*1e10,v
	endif
	return v
End






Function MakeNeutronWaves()
	Make/O/N=371/T root:Packages:Neutron:Isotope
	Wave/T Isotope=root:Packages:Neutron:Isotope
	Isotope[0,15] = {"H","1H","2H","3H","He","3He","4He","Li","6Li","7Li","Be","B","10B","11B","C","12C"}
	Isotope[16,30] = {"13C","N","14N","15N","O","16O","17O","18O","F","Ne","20Ne","21Ne","22Ne","Na","Mg"}
	Isotope[31,44] = {"24Mg","25Mg","26Mg","Al","Si","28Si","29Si","30Si","P","S","32S","33S","34S","36S"}
	Isotope[45,58] = {"Cl","35Cl","37Cl","Ar","36Ar","38Ar","40Ar","K","39K","40K","41K","Ca","40Ca","42Ca"}
	Isotope[59,71] = {"43Ca","44Ca","46Ca","48Ca","Sc","Ti","46Ti","47Ti","48Ti","49Ti","50Ti","V","50V"}
	Isotope[72,84] = {"51V","Cr","50Cr","52Cr","53Cr","54Cr","Mn","Fe","54Fe","56Fe","57Fe","58Fe","Co"}
	Isotope[85,97] = {"Ni","58Ni","60Ni","61Ni","62Ni","64Ni","Cu","63Cu","65Cu","Zn","64Zn","66Zn","67Zn"}
	Isotope[98,110] = {"68Zn","70Zn","Ga","69Ga","71Ga","Ge","70Ge","72Ge","73Ge","74Ge","76Ge","As","Se"}
	Isotope[111,123] = {"74Se","76Se","77Se","78Se","80Se","82Se","Br","79Br","81Br","Kr","78Kr","80Kr","82Kr"}
	Isotope[124,136] = {"83Kr","84Kr","86Kr","Rb","85Rb","87Rb","Sr","84Sr","86Sr","87Sr","88Sr","Y","Zr"}
	Isotope[137,149] = {"90Zr","91Zr","92Zr","94Zr","96Zr","Nb","Mo","92Mo","94Mo","95Mo","96Mo","97Mo","98Mo"}
	Isotope[150,162] = {"100Mo","Tc","Ru","96Ru","98Ru","99Ru","100Ru","101Ru","102Ru","104Ru","Rh","Pd","102Pd"}
	Isotope[163,173] = {"104Pd","105Pd","106Pd","108Pd","110Pd","Ag","107Ag","109Ag","Cd","106Cd","108Cd"}
	Isotope[174,184] = {"110Cd","111Cd","112Cd","113Cd","114Cd","116Cd","In","113In","115In","Sn","112Sn"}
	Isotope[185,195] = {"114Sn","115Sn","116Sn","117Sn","118Sn","119Sn","120Sn","122Sn","124Sn","Sb","121Sb"}
	Isotope[196,207] = {"123Sb","Te","120Te","122Te","123Te","124Te","125Te","126Te","128Te","130Te","I","Xe"}
	Isotope[208,218] = {"124Xe","126Xe","128Xe","129Xe","130Xe","131Xe","132Xe","134Xe","136Xe","Cs","Ba"}
	Isotope[219,229] = {"130Ba","132Ba","134Ba","135Ba","136Ba","137Ba","138Ba","La","138La","139La","Ce"}
	Isotope[230,240] = {"136Ce","138Ce","140Ce","142Ce","Pr","Nd","142Nd","143Nd","144Nd","145Nd","146Nd"}
	Isotope[241,251] = {"148Nd","150Nd","Pm","Sm","144Sm","147Sm","148Sm","149Sm","150Sm","152Sm","154Sm"}
	Isotope[252,262] = {"Eu","151Eu","153Eu","Gd","152Gd","154Gd","155Gd","156Gd","157Gd","158Gd","160Gd"}
	Isotope[263,274] = {"Tb","Dy","156Dy","158Dy","160Dy","161Dy","162Dy","163Dy","164Dy","Ho","Er","162Er"}
	Isotope[275,285] = {"164Er","166Er","167Er","168Er","170Er","Tm","Yb","168Yb","170Yb","171Yb","172Yb"}
	Isotope[286,296] = {"173Yb","174Yb","176Yb","Lu","175Lu","176Lu","Hf","174Hf","176Hf","177Hf","178Hf"}
	Isotope[297,309] = {"179Hf","180Hf","Ta","180Ta","181Ta","W","180W","182W","183W","184W","186W","Re","185Re"}
	Isotope[310,320] = {"187Re","Os","184Os","186Os","187Os","188Os","189Os","190Os","192Os","Ir","191Ir"}
	Isotope[321,332] = {"193Ir","Pt","190Pt","192Pt","194Pt","195Pt","196Pt","198Pt","Au","Hg","196Hg","198Hg"}
	Isotope[333,343] = {"199Hg","200Hg","201Hg","202Hg","204Hg","Tl","203Tl","205Tl","Pb","204Pb","206Pb"}
	Isotope[344,358] = {"207Pb","208Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","233U","234U","235U"}
	Isotope[359,370] = {"238U","Np","Pu","238Pu","239Pu","240Pu","242Pu","Am","Cm","244Cm","246Cm","248Cm"}

	Make/O/N=371/C root:Packages:Neutron:Coh_b
	Wave/C Coh_b=root:Packages:Neutron:Coh_b
	Coh_b[0,5] = {cmplx(-3.739,0),cmplx(-3.7406,0),cmplx(6.671,0),cmplx(4.792,0),cmplx(3.26,0),cmplx(5.74,-1.483)}
	Coh_b[6,11] = {cmplx(3.26,0),cmplx(-1.9,0),cmplx(2,-0.261),cmplx(-2.22,0),cmplx(7.79,0),cmplx(5.3,-0.213)}
	Coh_b[12,17] = {cmplx(-0.1,-1.066),cmplx(6.65,0),cmplx(6.646,0),cmplx(6.6511,0),cmplx(6.19,0),cmplx(9.36,0)}
	Coh_b[18,23] = {cmplx(9.37,0),cmplx(6.44,0),cmplx(5.803,0),cmplx(5.803,0),cmplx(5.78,0),cmplx(5.84,0)}
	Coh_b[24,29] = {cmplx(5.654,0),cmplx(4.566,0),cmplx(4.631,0),cmplx(6.66,0),cmplx(3.87,0),cmplx(3.63,0)}
	Coh_b[30,35] = {cmplx(5.375,0),cmplx(5.66,0),cmplx(3.62,0),cmplx(4.89,0),cmplx(3.449,0),cmplx(4.1491,0)}
	Coh_b[36,41] = {cmplx(4.107,0),cmplx(4.7,0),cmplx(4.58,0),cmplx(5.13,0),cmplx(2.847,0),cmplx(2.804,0)}
	Coh_b[42,47] = {cmplx(4.74,0),cmplx(3.48,0),cmplx(3,0),cmplx(9.577,0),cmplx(11.65,0),cmplx(3.08,0)}
	Coh_b[48,53] = {cmplx(1.909,0),cmplx(24.9,0),cmplx(3.5,0),cmplx(1.83,0),cmplx(3.67,0),cmplx(3.74,0)}
	Coh_b[54,60] = {cmplx(3,0),cmplx(2.69,0),cmplx(4.7,0),cmplx(4.8,0),cmplx(3.36,0),cmplx(-1.56,0),cmplx(1.42,0)}
	Coh_b[61,66] = {cmplx(3.6,0),cmplx(0.39,0),cmplx(12.29,0),cmplx(-3.438,0),cmplx(4.93,0),cmplx(3.63,0)}
	Coh_b[67,72] = {cmplx(-6.08,0),cmplx(1.04,0),cmplx(6.18,0),cmplx(-0.3824,0),cmplx(7.6,0),cmplx(-0.402,0)}
	Coh_b[73,78] = {cmplx(3.635,0),cmplx(-4.5,0),cmplx(4.92,0),cmplx(-4.2,0),cmplx(4.55,0),cmplx(-3.73,0)}
	Coh_b[79,85] = {cmplx(9.45,0),cmplx(4.2,0),cmplx(9.94,0),cmplx(2.3,0),cmplx(15,0),cmplx(2.49,0),cmplx(10.3,0)}
	Coh_b[86,91] = {cmplx(14.4,0),cmplx(2.8,0),cmplx(7.6,0),cmplx(-8.7,0),cmplx(-0.37,0),cmplx(7.718,0)}
	Coh_b[92,97] = {cmplx(6.43,0),cmplx(10.61,0),cmplx(5.68,0),cmplx(5.22,0),cmplx(5.97,0),cmplx(7.56,0)}
	Coh_b[98,103] = {cmplx(6.03,0),cmplx(6,0),cmplx(7.288,0),cmplx(7.88,0),cmplx(6.4,0),cmplx(8.185,0)}
	Coh_b[104,110] = {cmplx(10,0),cmplx(8.51,0),cmplx(5.02,0),cmplx(7.58,0),cmplx(8.2,0),cmplx(6.58,0),cmplx(7.97,0)}
	Coh_b[111,116] = {cmplx(0.8,0),cmplx(12.2,0),cmplx(8.25,0),cmplx(8.24,0),cmplx(7.48,0),cmplx(6.34,0)}
	Coh_b[117,122] = {cmplx(6.795,0),cmplx(6.8,0),cmplx(6.79,0),cmplx(7.81,0),cmplx(Nan,Nan),cmplx(Nan,Nan)}
	Coh_b[123,128] = {cmplx(Nan,Nan),cmplx(Nan,Nan),cmplx(Nan,Nan),cmplx(8.1,0),cmplx(7.09,0),cmplx(7.03,0)}
	Coh_b[129,135] = {cmplx(7.23,0),cmplx(7.02,0),cmplx(7,0),cmplx(5.67,0),cmplx(7.4,0),cmplx(7.15,0),cmplx(7.75,0)}
	Coh_b[136,142] = {cmplx(7.16,0),cmplx(6.4,0),cmplx(8.7,0),cmplx(7.4,0),cmplx(8.2,0),cmplx(5.5,0),cmplx(7.054,0)}
	Coh_b[143,148] = {cmplx(6.715,0),cmplx(6.91,0),cmplx(6.8,0),cmplx(6.91,0),cmplx(6.2,0),cmplx(7.24,0)}
	Coh_b[149,154] = {cmplx(6.58,0),cmplx(6.73,0),cmplx(6.8,0),cmplx(7.03,0),cmplx(Nan,Nan),cmplx(Nan,Nan)}
	Coh_b[155,160] = {cmplx(Nan,Nan),cmplx(Nan,Nan),cmplx(Nan,Nan),cmplx(Nan,Nan),cmplx(Nan,Nan),cmplx(5.88,0)}
	Coh_b[161,167] = {cmplx(5.91,0),cmplx(7.7,0),cmplx(7.7,0),cmplx(5.5,0),cmplx(6.4,0),cmplx(4.1,0),cmplx(7.7,0)}
	Coh_b[168,173] = {cmplx(5.922,0),cmplx(7.555,0),cmplx(4.165,0),cmplx(4.87,-0.7),cmplx(5,0),cmplx(5.4,0)}
	Coh_b[174,180] = {cmplx(5.9,0),cmplx(6.5,0),cmplx(6.4,0),cmplx(-8,-5.73),cmplx(7.5,0),cmplx(6.3,0),cmplx(4.065,-0.0539)}
	Coh_b[181,186] = {cmplx(5.39,0),cmplx(4.01,-0.0562),cmplx(6.225,0),cmplx(6,0),cmplx(6.2,0),cmplx(6,0)}
	Coh_b[187,192] = {cmplx(5.93,0),cmplx(6.48,0),cmplx(6.07,0),cmplx(6.12,0),cmplx(6.49,0),cmplx(5.74,0)}
	Coh_b[193,198] = {cmplx(5.97,0),cmplx(5.57,0),cmplx(5.71,0),cmplx(5.38,0),cmplx(5.8,0),cmplx(5.3,0)}
	Coh_b[199,204] = {cmplx(3.8,0),cmplx(-0.05,-0.116),cmplx(7.96,0),cmplx(5.02,0),cmplx(5.56,0),cmplx(5.89,0)}
	Coh_b[205,210] = {cmplx(6.02,0),cmplx(5.28,0),cmplx(4.92,0),cmplx(Nan,Nan),cmplx(Nan,Nan),cmplx(Nan,Nan)}
	Coh_b[211,216] = {cmplx(Nan,Nan),cmplx(Nan,Nan),cmplx(Nan,Nan),cmplx(Nan,Nan),cmplx(Nan,Nan),cmplx(Nan,Nan)}
	Coh_b[217,222] = {cmplx(5.42,0),cmplx(5.07,0),cmplx(-3.6,0),cmplx(7.8,0),cmplx(5.7,0),cmplx(4.67,0)}
	Coh_b[223,229] = {cmplx(4.91,0),cmplx(6.83,0),cmplx(4.84,0),cmplx(8.24,0),cmplx(8,0),cmplx(8.24,0),cmplx(4.84,0)}
	Coh_b[230,235] = {cmplx(5.8,0),cmplx(6.7,0),cmplx(4.84,0),cmplx(4.75,0),cmplx(4.58,0),cmplx(7.69,0)}
	Coh_b[236,242] = {cmplx(7.7,0),cmplx(14,0),cmplx(2.8,0),cmplx(14,0),cmplx(8.7,0),cmplx(5.7,0),cmplx(5.3,0)}
	Coh_b[243,248] = {cmplx(12.6,0),cmplx(0.8,-1.65),cmplx(-3,0),cmplx(14,0),cmplx(-3,0),cmplx(-19.2,-11.7)}
	Coh_b[249,254] = {cmplx(14,0),cmplx(-5,0),cmplx(9.3,0),cmplx(7.22,-1.26),cmplx(6.13,-2.53),cmplx(8.22,0)}
	Coh_b[255,260] = {cmplx(6.5,-13.82),cmplx(10,0),cmplx(10,0),cmplx(6,-17),cmplx(6.3,0),cmplx(-1.14,-71.9)}
	Coh_b[261,266] = {cmplx(9,0),cmplx(9.15,0),cmplx(7.38,0),cmplx(16.9,-0.276),cmplx(6.1,0),cmplx(6,0)}
	Coh_b[267,272] = {cmplx(6.7,0),cmplx(10.3,0),cmplx(-1.4,0),cmplx(5,0),cmplx(49.4,-0.79),cmplx(8.01,0)}
	Coh_b[273,279] = {cmplx(7.79,0),cmplx(8.8,0),cmplx(8.2,0),cmplx(10.6,0),cmplx(3,0),cmplx(7.4,0),cmplx(9.6,0)}
	Coh_b[280,285] = {cmplx(7.07,0),cmplx(12.43,0),cmplx(-4.07,-0.62),cmplx(6.77,0),cmplx(9.66,0),cmplx(9.43,0)}
	Coh_b[286,291] = {cmplx(9.56,0),cmplx(19.3,0),cmplx(8.72,0),cmplx(7.21,0),cmplx(7.24,0),cmplx(6.1,-0.57)}
	Coh_b[292,298] = {cmplx(7.7,0),cmplx(10.9,0),cmplx(6.61,0),cmplx(0.8,0),cmplx(5.9,0),cmplx(7.46,0),cmplx(13.2,0)}
	Coh_b[299,305] = {cmplx(6.91,0),cmplx(7,0),cmplx(6.91,0),cmplx(4.86,0),cmplx(5,0),cmplx(6.97,0),cmplx(6.53,0)}
	Coh_b[306,312] = {cmplx(7.48,0),cmplx(-0.72,0),cmplx(9.2,0),cmplx(9,0),cmplx(9.3,0),cmplx(10.7,0),cmplx(10,0)}
	Coh_b[313,319] = {cmplx(11.6,0),cmplx(10,0),cmplx(7.6,0),cmplx(10.7,0),cmplx(11,0),cmplx(11.5,0),cmplx(10.6,0)}
	Coh_b[320,325] = {cmplx(Nan,Nan),cmplx(Nan,Nan),cmplx(9.6,0),cmplx(9,0),cmplx(9.9,0),cmplx(10.55,0)}
	Coh_b[326,331] = {cmplx(8.83,0),cmplx(9.89,0),cmplx(7.8,0),cmplx(7.63,0),cmplx(12.692,0),cmplx(30.3,0)}
	Coh_b[332,337] = {cmplx(Nan,Nan),cmplx(16.9,0),cmplx(Nan,Nan),cmplx(Nan,Nan),cmplx(Nan,Nan),cmplx(Nan,Nan)}
	Coh_b[338,343] = {cmplx(8.776,0),cmplx(6.99,0),cmplx(9.52,0),cmplx(9.405,0),cmplx(9.9,0),cmplx(9.22,0)}
	Coh_b[344,349] = {cmplx(9.28,0),cmplx(9.5,0),cmplx(8.532,0),cmplx(Nan,Nan),cmplx(Nan,Nan),cmplx(Nan,Nan)}
	Coh_b[350,355] = {cmplx(Nan,Nan),cmplx(10,0),cmplx(Nan,Nan),cmplx(10.31,0),cmplx(9.1,0),cmplx(8.417,0)}
	Coh_b[356,361] = {cmplx(10.1,0),cmplx(12.4,0),cmplx(10.47,0),cmplx(8.402,0),cmplx(10.55,0),cmplx(Nan,Nan)}
	Coh_b[362,368] = {cmplx(14.1,0),cmplx(7.7,0),cmplx(3.5,0),cmplx(8.1,0),cmplx(8.3,0),cmplx(Nan,Nan),cmplx(9.5,0)}
	Coh_b[369,370] = {cmplx(9.3,0),cmplx(7.7,0)}

	Make/O/N=371/C root:Packages:Neutron:Inc_b
	Wave/C Inc_b=root:Packages:Neutron:Inc_b
	Inc_b[0,5] = {cmplx(Nan,Nan),cmplx(25.274,0),cmplx(4.04,0),cmplx(-1.04,0),cmplx(Nan,Nan),cmplx(-2.5,2.568)}
	Inc_b[6,11] = {cmplx(0,0),cmplx(Nan,Nan),cmplx(-1.89,0.26),cmplx(-2.49,0),cmplx(0.12,0),cmplx(Nan,Nan)}
	Inc_b[12,17] = {cmplx(-4.7,1.231),cmplx(-1.3,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(-0.52,0),cmplx(Nan,Nan)}
	Inc_b[18,24] = {cmplx(2,0),cmplx(-0.02,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(0.18,0),cmplx(0,0),cmplx(-0.082,0)}
	Inc_b[25,31] = {cmplx(Nan,Nan),cmplx(0,0),cmplx(0.6,0),cmplx(0,0),cmplx(3.59,0),cmplx(Nan,Nan),cmplx(0,0)}
	Inc_b[32,38] = {cmplx(1.48,0),cmplx(0,0),cmplx(0.256,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(0.09,0),cmplx(0,0)}
	Inc_b[39,45] = {cmplx(0.2,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(1.5,0),cmplx(0,0),cmplx(0,0),cmplx(Nan,Nan)}
	Inc_b[46,52] = {cmplx(6.1,0),cmplx(0.1,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(0,0),cmplx(0,0),cmplx(Nan,Nan)}
	Inc_b[53,59] = {cmplx(1.4,0),cmplx(Nan,Nan),cmplx(1.5,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(0,0),cmplx(Nan,Nan)}
	Inc_b[60,66] = {cmplx(0,0),cmplx(0,0),cmplx(0,0),cmplx(-6,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(-3.5,0)}
	Inc_b[67,73] = {cmplx(0,0),cmplx(5.1,0),cmplx(0,0),cmplx(Nan,Nan),cmplx(Nan,Nan),cmplx(6.35,0),cmplx(Nan,Nan)}
	Inc_b[74,80] = {cmplx(0,0),cmplx(0,0),cmplx(6.87,0),cmplx(0,0),cmplx(1.79,0),cmplx(Nan,Nan),cmplx(0,0)}
	Inc_b[81,87] = {cmplx(0,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(-6.2,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(0,0)}
	Inc_b[88,94] = {cmplx(3.9,0),cmplx(0,0),cmplx(0,0),cmplx(Nan,Nan),cmplx(0.22,0),cmplx(1.79,0),cmplx(Nan,Nan)}
	Inc_b[95,101] = {cmplx(0,0),cmplx(0,0),cmplx(-1.5,0),cmplx(0,0),cmplx(0,0),cmplx(Nan,Nan),cmplx(-0.85,0)}
	Inc_b[102,108] = {cmplx(-0.82,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(0,0),cmplx(3.4,0),cmplx(0,0),cmplx(0,0)}
	Inc_b[109,115] = {cmplx(-0.69,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(0,0),cmplx(0.6,0),cmplx(0,0),cmplx(0,0)}
	Inc_b[116,122] = {cmplx(0,0),cmplx(Nan,Nan),cmplx(-1.1,0),cmplx(0.6,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(0,0)}
	Inc_b[123,129] = {cmplx(0,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(0,0),cmplx(Nan,Nan),cmplx(Nan,Nan),cmplx(Nan,Nan)}
	Inc_b[130,136] = {cmplx(Nan,Nan),cmplx(0,0),cmplx(0,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(1.1,0),cmplx(Nan,Nan)}
	Inc_b[137,143] = {cmplx(0,0),cmplx(-1.08,0),cmplx(0,0),cmplx(0,0),cmplx(0,0),cmplx(-0.139,0),cmplx(Nan,Nan)}
	Inc_b[144,150] = {cmplx(0,0),cmplx(0,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(0,0)}
	Inc_b[151,157] = {cmplx(Nan,Nan),cmplx(Nan,Nan),cmplx(0,0),cmplx(0,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(Nan,Nan)}
	Inc_b[158,164] = {cmplx(0,0),cmplx(0,0),cmplx(Nan,Nan),cmplx(Nan,Nan),cmplx(0,0),cmplx(0,0),cmplx(-2.6,0)}
	Inc_b[165,171] = {cmplx(0,0),cmplx(0,0),cmplx(0,0),cmplx(Nan,Nan),cmplx(1,0),cmplx(-1.6,0),cmplx(Nan,Nan)}
	Inc_b[172,178] = {cmplx(0,0),cmplx(0,0),cmplx(0,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(Nan,Nan),cmplx(0,0)}
	Inc_b[179,185] = {cmplx(0,0),cmplx(Nan,Nan),cmplx(0.017,0),cmplx(-2.1,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(0,0)}
	Inc_b[186,192] = {cmplx(Nan,Nan),cmplx(0,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(0,0)}
	Inc_b[193,199] = {cmplx(0,0),cmplx(Nan,Nan),cmplx(-0.05,0),cmplx(-0.1,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(0,0)}
	Inc_b[200,206] = {cmplx(-2.04,0),cmplx(0,0),cmplx(-0.26,0),cmplx(0,0),cmplx(0,0),cmplx(0,0),cmplx(1.58,0)}
	Inc_b[207,213] = {cmplx(3.04,0),cmplx(0,0),cmplx(0,0),cmplx(0,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(Nan,Nan)}
	Inc_b[214,220] = {cmplx(0,0),cmplx(0,0),cmplx(0,0),cmplx(1.29,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(0,0)}
	Inc_b[221,227] = {cmplx(0,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(Nan,Nan),cmplx(Nan,Nan)}
	Inc_b[228,234] = {cmplx(3,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(0,0),cmplx(0,0),cmplx(0,0),cmplx(-0.35,0)}
	Inc_b[235,241] = {cmplx(Nan,Nan),cmplx(0,0),cmplx(21,0),cmplx(0,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(0,0)}
	Inc_b[242,248] = {cmplx(0,0),cmplx(3.2,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(11,0),cmplx(0,0),cmplx(31.4,-10.3)}
	Inc_b[249,255] = {cmplx(0,0),cmplx(0,0),cmplx(0,0),cmplx(Nan,Nan),cmplx(4.5,-2.14),cmplx(3.2,0),cmplx(Nan,Nan)}
	Inc_b[256,262] = {cmplx(0,0),cmplx(0,0),cmplx(5,-13.16),cmplx(0,0),cmplx(5,-55.8),cmplx(0,0),cmplx(0,0)}
	Inc_b[263,269] = {cmplx(-0.17,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(0,0),cmplx(0,0),cmplx(4.9,0),cmplx(0,0)}
	Inc_b[270,276] = {cmplx(1.3,0),cmplx(0,0),cmplx(-1.7,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(0,0),cmplx(0,0)}
	Inc_b[277,283] = {cmplx(1,0),cmplx(0,0),cmplx(0,0),cmplx(0.9,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(0,0)}
	Inc_b[284,290] = {cmplx(-5.59,0),cmplx(0,0),cmplx(-5.3,0),cmplx(0,0),cmplx(0,0),cmplx(Nan,Nan),cmplx(2.2,0)}
	Inc_b[291,297] = {cmplx(3,0.61),cmplx(Nan,Nan),cmplx(0,0),cmplx(0,0),cmplx(0.9,0),cmplx(0,0),cmplx(1.06,0)}
	Inc_b[298,303] = {cmplx(0,0),cmplx(Nan,Nan),cmplx(Nan,Nan),cmplx(-0.29,0),cmplx(Nan,Nan),cmplx(0,0)}
	Inc_b[304,310] = {cmplx(0,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(0,0),cmplx(Nan,Nan),cmplx(2,0),cmplx(2.8,0)}
	Inc_b[311,317] = {cmplx(Nan,Nan),cmplx(0,0),cmplx(0,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(Nan,Nan),cmplx(0,0)}
	Inc_b[318,323] = {cmplx(0,0),cmplx(Nan,Nan),cmplx(Nan,Nan),cmplx(Nan,Nan),cmplx(Nan,Nan),cmplx(0,0)}
	Inc_b[324,330] = {cmplx(0,0),cmplx(0,0),cmplx(-1,0),cmplx(0,0),cmplx(0,0),cmplx(-1.84,0),cmplx(Nan,Nan)}
	Inc_b[331,337] = {cmplx(0,0),cmplx(0,0),cmplx(15.5,0),cmplx(0,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(0,0)}
	Inc_b[338,343] = {cmplx(Nan,Nan),cmplx(1.06,0),cmplx(-0.242,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(0,0)}
	Inc_b[344,349] = {cmplx(0.14,0),cmplx(0,0),cmplx(Nan,Nan),cmplx(0.259,0),cmplx(Nan,Nan),cmplx(Nan,Nan)}
	Inc_b[350,355] = {cmplx(Nan,Nan),cmplx(0,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(Nan,Nan),cmplx(Nan,Nan)}
	Inc_b[356,362] = {cmplx(1,0),cmplx(0,0),cmplx(1.3,0),cmplx(0,0),cmplx(Nan,Nan),cmplx(Nan,Nan),cmplx(0,0)}
	Inc_b[363,369] = {cmplx(1.3,0),cmplx(0,0),cmplx(0,0),cmplx(2,0),cmplx(Nan,Nan),cmplx(0,0),cmplx(0,0)}
	Inc_b[370,370] = {cmplx(0,0)}

	Make/O/N=371 root:Packages:Neutron:conc
	Wave conc=root:Packages:Neutron:conc
	conc[0,18] = {Nan,99.985,0.015,0,Nan,0.00014,100,Nan,7.5,92.5,100,Nan,20,80,Nan,98.9,1.1,Nan,99.63}
	conc[19,36] = {0.37,Nan,99.762,0.038,0.2,100,Nan,90.51,0.27,9.22,100,Nan,78.99,10,11.01,100,Nan,92.23}
	conc[37,53] = {4.67,3.1,100,Nan,95.02,0.75,4.21,0.02,Nan,75.77,24.23,Nan,0.337,0.063,99.6,Nan,93.258}
	conc[54,70] = {0.012,6.73,Nan,96.941,0.647,0.135,2.086,0.004,0.187,100,Nan,8.2,7.4,73.8,5.4,5.2,Nan}
	conc[71,88] = {0.25,99.75,Nan,4.35,83.79,9.5,2.36,100,Nan,5.8,91.7,2.2,0.3,100,Nan,68.27,26.1,1.13}
	conc[89,106] = {3.59,0.91,Nan,69.17,30.83,Nan,48.6,27.9,4.1,18.8,0.6,Nan,60.1,39.9,Nan,20.5,27.4,7.8}
	conc[107,125] = {36.5,7.8,100,Nan,0.9,9,7.6,23.5,49.6,9.4,Nan,50.69,49.31,Nan,0.35,2.25,11.6,11.5,57}
	conc[126,142] = {17.3,Nan,72.17,27.83,Nan,0.56,9.86,7,82.58,100,Nan,51.45,11.32,17.19,17.28,2.76,100}
	conc[143,160] = {Nan,14.84,9.25,15.92,16.68,9.55,24.13,9.63,0,Nan,5.5,1.9,12.7,12.6,17,31.6,18.7,100}
	conc[161,176] = {Nan,1.02,11.14,22.33,27.33,26.46,11.72,Nan,51.83,48.17,Nan,1.25,0.89,12.51,12.81,24.13}
	conc[177,195] = {12.22,28.72,7.47,Nan,4.3,95.7,Nan,1,0.7,0.4,14.7,7.7,24.3,8.6,32.4,4.6,5.6,Nan,57.3}
	conc[196,212] = {42.7,Nan,0.096,2.6,0.908,4.816,7.14,18.95,31.69,33.8,100,Nan,0.1,0.09,1.91,26.4,4.1}
	conc[213,230] = {21.2,26.9,10.4,8.9,100,Nan,0.11,0.1,2.42,6.59,7.85,11.23,71.7,Nan,0.09,99.91,Nan,0.19}
	conc[231,247] = {0.25,88.48,11.08,100,Nan,27.16,12.18,23.8,8.29,17.19,5.75,5.63,0,Nan,3.1,15.1,11.3}
	conc[248,265] = {13.9,7.4,26.6,22.6,Nan,47.8,52.2,Nan,0.2,2.1,14.8,20.6,15.7,24.8,21.8,100,Nan,0.06}
	conc[266,283] = {0.1,2.34,19,25.5,24.9,28.1,100,Nan,0.14,1.56,33.4,22.9,27.1,14.9,100,Nan,0.14,3.06}
	conc[284,300] = {14.3,21.9,16.1,31.8,12.7,Nan,97.39,2.61,Nan,0.2,5.2,18.6,27.1,13.7,35.2,Nan,0.012}
	conc[301,317] = {99.988,Nan,0.1,26.3,14.3,30.7,28.6,Nan,37.4,62.6,Nan,0.02,1.58,1.6,13.3,16.1,26.4}
	conc[318,336] = {41,Nan,37.3,62.7,Nan,0.01,0.79,32.9,33.8,25.3,7.2,100,Nan,0.2,10.1,17,23.1,13.2,29.6}
	conc[337,356] = {6.8,Nan,29.524,70.476,Nan,1.4,24.1,22.1,52.4,100,Nan,Nan,Nan,Nan,0,Nan,100,0,Nan,0}
	conc[357,370] = {0.005,0.72,99.275,0,Nan,0,0,0,0,0,Nan,0,0,0}

	Make/O/N=371 root:Packages:Neutron:Zatom
	Wave Zatom=root:Packages:Neutron:Zatom
	Zatom[0,35] = {1,1,1,1,2,2,2,3,3,3,4,5,5,5,6,6,6,7,7,7,8,8,8,8,9,10,10,10,10,11,12,12,12,12,13,14}
	Zatom[36,63] = {14,14,14,15,16,16,16,16,16,17,17,17,18,18,18,18,19,19,19,19,20,20,20,20,20,20,20,21}
	Zatom[64,91] = {22,22,22,22,22,22,23,23,23,24,24,24,24,24,25,26,26,26,26,26,27,28,28,28,28,28,28,29}
	Zatom[92,119] = {29,29,30,30,30,30,30,30,31,31,31,32,32,32,32,32,32,33,34,34,34,34,34,34,34,35,35,35}
	Zatom[120,147] = {36,36,36,36,36,36,36,37,37,37,38,38,38,38,38,39,40,40,40,40,40,40,41,42,42,42,42,42}
	Zatom[148,175] = {42,42,42,43,44,44,44,44,44,44,44,44,45,46,46,46,46,46,46,46,47,47,47,48,48,48,48,48}
	Zatom[176,203] = {48,48,48,48,49,49,49,50,50,50,50,50,50,50,50,50,50,50,51,51,51,52,52,52,52,52,52,52}
	Zatom[204,231] = {52,52,53,54,54,54,54,54,54,54,54,54,54,55,56,56,56,56,56,56,56,56,57,57,57,58,58,58}
	Zatom[232,259] = {58,58,59,60,60,60,60,60,60,60,60,61,62,62,62,62,62,62,62,62,63,63,63,64,64,64,64,64}
	Zatom[260,287] = {64,64,64,65,66,66,66,66,66,66,66,66,67,68,68,68,68,68,68,68,69,70,70,70,70,70,70,70}
	Zatom[288,315] = {70,71,71,71,72,72,72,72,72,72,72,73,73,73,74,74,74,74,74,74,75,75,75,76,76,76,76,76}
	Zatom[316,343] = {76,76,76,77,77,77,78,78,78,78,78,78,78,79,80,80,80,80,80,80,80,80,81,81,81,82,82,82}
	Zatom[344,370] = {82,82,83,84,85,86,87,88,89,90,91,92,92,92,92,92,93,94,94,94,94,94,95,96,96,96,96}

	Make/O/N=371 root:Packages:Neutron:mass
	Wave mass=root:Packages:Neutron:mass
	mass[0,26] = {Nan,1,2,3,Nan,3,4,Nan,6,7,Nan,Nan,10,11,Nan,12,13,Nan,14,15,Nan,16,17,18,Nan,Nan,20}
	mass[27,51] = {21,22,Nan,Nan,24,25,26,Nan,Nan,28,29,30,Nan,Nan,32,33,34,36,Nan,35,37,Nan,36,38,40}
	mass[52,77] = {Nan,39,40,41,Nan,40,42,43,44,46,48,Nan,Nan,46,47,48,49,50,Nan,50,51,Nan,50,52,53,54}
	mass[78,102] = {Nan,Nan,54,56,57,58,Nan,Nan,58,60,61,62,64,Nan,63,65,Nan,64,66,67,68,70,Nan,69,71}
	mass[103,128] = {Nan,70,72,73,74,76,Nan,Nan,74,76,77,78,80,82,Nan,79,81,Nan,78,80,82,83,84,86,Nan,85}
	mass[129,153] = {87,Nan,84,86,87,88,Nan,Nan,90,91,92,94,96,Nan,Nan,92,94,95,96,97,98,100,Nan,Nan,96}
	mass[154,174] = {98,99,100,101,102,104,Nan,Nan,102,104,105,106,108,110,Nan,107,109,Nan,106,108,110}
	mass[175,195] = {111,112,113,114,116,Nan,113,115,Nan,112,114,115,116,117,118,119,120,122,124,Nan,121}
	mass[196,216] = {123,Nan,120,122,123,124,125,126,128,130,Nan,Nan,124,126,128,129,130,131,132,134,136}
	mass[217,237] = {Nan,Nan,130,132,134,135,136,137,138,Nan,138,139,Nan,136,138,140,142,Nan,Nan,142,143}
	mass[238,258] = {144,145,146,148,150,Nan,Nan,144,147,148,149,150,152,154,Nan,151,153,Nan,152,154,155}
	mass[259,279] = {156,157,158,160,Nan,Nan,156,158,160,161,162,163,164,Nan,Nan,162,164,166,167,168,170}
	mass[280,300] = {Nan,Nan,168,170,171,172,173,174,176,Nan,175,176,Nan,174,176,177,178,179,180,Nan,180}
	mass[301,321] = {181,Nan,180,182,183,184,186,Nan,185,187,Nan,184,186,187,188,189,190,192,Nan,191,193}
	mass[322,342] = {Nan,190,192,194,195,196,198,Nan,Nan,196,198,199,200,201,202,204,Nan,203,205,Nan,204}
	mass[343,363] = {206,207,208,Nan,Nan,Nan,Nan,Nan,Nan,Nan,Nan,Nan,Nan,233,234,235,238,Nan,Nan,238,239}
	mass[364,370] = {240,242,Nan,Nan,244,246,248}

	Make/O/N=371 root:Packages:Neutron:amuIsotope
	Wave amuIsotope=root:Packages:Neutron:amuIsotope
	amuIsotope[0,11] = {1.00797,1.00783,2.0141,3.01605,4.0026,3.01603,4.0026,6.941,6.01512,7.016,9.01218,10.81}
	amuIsotope[12,22] = {10.0129,11.0093,12.011,12.0989,13.0034,14.0067,14.0031,15.0001,15.9994,15.9949,16.9991}
	amuIsotope[23,33] = {17.9992,18.9984,20.179,19.9924,20.9938,21.9914,22.9898,24.305,23.985,24.9858,25.9826}
	amuIsotope[34,44] = {26.9815,28.0855,27.9769,28.9765,29.9738,30.9738,32.06,31.9721,32.9715,33.9679,35.9671}
	amuIsotope[45,55] = {35.453,34.9689,36.9659,39.948,35.9675,37.9627,39.9624,39.0983,38.9637,39.964,40.9618}
	amuIsotope[56,66] = {40.08,39.9626,41.9586,42.9588,43.9555,45.9537,47.9525,44.9559,47.88,45.9526,46.9518}
	amuIsotope[67,77] = {47.9479,48.9479,49.9448,50.9415,49.9472,50.944,51.996,49.946,51.9405,52.9407,53.9389}
	amuIsotope[78,88] = {54.938,55.847,53.9396,55.9349,56.9354,57.9333,58.9332,58.69,57.9353,59.9308,60.9311}
	amuIsotope[89,99] = {61.9283,63.928,63.546,62.9296,64.9278,65.38,63.9291,65.926,66.9271,67.9249,69.9253}
	amuIsotope[100,110] = {69.72,68.9256,70.9247,72.59,69.9242,71.9221,72.9235,73.9212,75.9214,74.9216,78.96}
	amuIsotope[111,121] = {73.9225,75.9192,76.9199,77.9173,79.9165,81.9167,79.904,78.9183,80.9163, 83.8,77.9204}
	amuIsotope[122,132] = {79.9164,81.9135,82.9141,83.9115,85.9106,85.4678,84.9118,86.9092,87.62,83.9134,85.9093}
	amuIsotope[133,143] = {86.9089,87.9056,88.9059,91.22,89.9047,90.9056,91.905,93.9063,95.9083,92.9064,95.94}
	amuIsotope[144,154] = {91.9068,93.9051,94.9058,95.9047,96.906,97.9054,99.9075,   98,101.07,95.9076,97.9053}
	amuIsotope[155,165] = {98.9059,99.9042,100.906,101.904,103.905,102.906,106.42,101.906,103.904,104.905,105.903}
	amuIsotope[166,176] = {107.904,109.905,107.868,106.905,108.905,112.41,105.906,107.904,109.903,110.904,111.903}
	amuIsotope[177,187] = {112.904,113.903,115.905,114.82,112.904,114.904,118.69,111.905,113.903,114.903,115.902}
	amuIsotope[188,198] = {116.903,117.902,118.903,119.902,121.903,123.905,121.75,120.904,122.904,127.6,119.904}
	amuIsotope[199,209] = {121.903,122.904,123.903,124.904,125.903,127.904,129.906,126.905,131.29,123.906,125.904}
	amuIsotope[210,220] = {127.904,128.905,129.904,130.905,131.904,133.905,135.907,132.905,137.33,129.906,131.905}
	amuIsotope[221,231] = {133.905,134.906,135.905,136.906,137.905,138.906,137.907,138.906,140.12,135.907,137.906}
	amuIsotope[232,242] = {139.905,141.909,140.908,144.24,141.908,142.91,143.91,144.913,145.913,147.917,149.921}
	amuIsotope[243,253] = {  145,150.36,143.912,146.915,147.915,148.917,149.917,151.92,153.922,151.96,150.92}
	amuIsotope[254,264] = {152.921,157.25,151.92,153.921,154.923,155.922,156.924,157.924,159.927,158.925,162.5}
	amuIsotope[265,275] = {155.924,157.924,159.925,160.927,161.927,162.929,163.929,164.93,167.26,161.929,163.929}
	amuIsotope[276,286] = {165.93,166.932,167.932,169.935,168.934,173.04,167.934,169.935,170.936,171.936,172.938}
	amuIsotope[287,297] = {173.939,175.943,174.967,174.941,175.943,178.49,173.94,175.941,176.943,177.944,178.946}
	amuIsotope[298,308] = {179.947,180.948,179.947,180.948,183.85,179.947,181.948,182.95,183.951,185.954,186.207}
	amuIsotope[309,319] = {184.953,186.956,190.2,183.952,185.954,186.956,187.956,188.958,189.958,191.961,192.22}
	amuIsotope[320,330] = {190.961,192.963,195.08,189.96,191.961,193.963,194.965,195.965,197.968,196.967,200.59}
	amuIsotope[331,341] = {195.966,197.967,198.968,199.968,200.97,201.971,203.973,204.383,202.972,204.974,207.2}
	amuIsotope[342,353] = {203.973,205.974,206.976,207.977,208.98,  209,  210,  222,  223,226.025,227.028,232.038}
	amuIsotope[354,364] = {231.036,238.029,233.04,234.041,235.044,238.051,  237,  242,238.05,239.052,240.054}
	amuIsotope[365,370] = {242.059,  243,  247,244.063,246.067,248.072}

	Make/O/N=371 root:Packages:Neutron:Coh_xs
	Wave Coh_xs=root:Packages:Neutron:Coh_xs
	Coh_xs[0,14] = {1.7568,1.7583,5.592,2.89,1.34,4.42,1.34,0.454,0.51,0.619,7.63,3.54,0.144,5.56,5.551}
	Coh_xs[15,30] = {5.559,4.81,11.01,11.03,5.21,4.232,4.232,4.2,4.29,4.017,2.62,2.695,5.6,1.88,1.66,3.631}
	Coh_xs[31,46] = {4.03,1.65,3,1.495,2.163,2.12,2.78,2.64,3.307,1.0186,0.988,2.8,1.52,1.1,11.526,17.06}
	Coh_xs[47,63] = {1.19,0.458,77.9,1.5,0.421,1.69,1.76,1.1,0.91,2.78,2.9,1.42,0.31,0.25,1.6,0.019,19}
	Coh_xs[64,79] = {1.485,3.05,1.66,4.65,0.14,4.8,0.0184,7.3,0.0203,1.66,2.54,3.042,2.22,2.6,1.75,11.22}
	Coh_xs[80,96] = {2.2,12.42,0.66,28,0.779,13.3,26.1,0.99,7.26,9.5,0.017,7.485,5.2,14.1,4.054,3.42,4.48}
	Coh_xs[97,115] = {7.18,4.57,4.5,6.675,7.8,5.15,8.42,12.6,9.1,3.17,7.2,8,5.44,7.98,0.1,18.7,8.6,8.5,7.03}
	Coh_xs[116,134] = {5.05,5.8,5.81,5.79,7.67,Nan,Nan,Nan,Nan,Nan,8.2,6.32,6.2,6.6,6.19,6,4.04,6.88,6.42}
	Coh_xs[135,153] = {7.55,6.44,5.1,9.5,6.9,8.4,3.8,6.253,5.67,6,5.81,6,4.83,6.59,5.44,5.69,5.8,6.21,Nan}
	Coh_xs[154,172] = {Nan,Nan,Nan,Nan,Nan,Nan,4.34,4.39,7.5,7.5,3.8,5.1,2.1,7.5,4.407,7.17,2.18,3.04,3.1}
	Coh_xs[173,191] = {3.7,4.4,5.3,5.1,12.1,7.1,5,2.08,3.65,2.02,4.871,4.5,4.8,4.5,4.42,5.28,4.63,4.71,5.29}
	Coh_xs[192,210] = {4.14,4.48,3.9,4.1,3.64,4.23,3.5,1.8,0.002,8,3.17,3.88,4.36,4.55,3.5,2.96,Nan,Nan,Nan}
	Coh_xs[211,229] = {Nan,Nan,Nan,Nan,Nan,Nan,3.69,3.23,1.6,7.6,4.08,2.74,3.03,5.86,2.94,8.53,8,8.53,2.94}
	Coh_xs[230,250] = {4.23,5.64,2.94,2.84,2.64,7.43,7.5,25,1,25,9.5,4.1,3.5,20,0.422,1,25,1,63.5,25,3.1}
	Coh_xs[251,271] = {11,6.57,5.5,8.5,29.3,13,13,40.8,5,650,10,10.52,6.84,35.9,4.7,5,5.6,13.3,0.25,3.1,307}
	Coh_xs[272,289] = {8.06,7.63,9.7,8.4,14.1,1.1,6.9,11.6,6.28,19.42,2.13,5.8,11.7,11.2,11.5,46.8,9.6,6.53}
	Coh_xs[290,310] = {6.59,4.7,7.6,15,5.5,0.1,4.4,7,21.9,6,6.2,6,2.97,3,6.1,5.36,7.03,0.065,10.6,10.2,10.9}
	Coh_xs[311,330] = {14.4,13,17,13,7.3,14.4,15.2,16.6,14.1,Nan,Nan,11.58,10,12.3,14,9.8,12.3,7.6,7.32,20.24}
	Coh_xs[331,347] = {115,Nan,36,Nan,Nan,Nan,Nan,9.678,6.14,11.39,11.115,12.3,10.68,10.82,11.34,9.148,0}
	Coh_xs[348,368] = {0,0,0,13,0,13.36,10.4,8.903,12.8,19.3,13.78,8.871,14,Nan,25,7.5,1.54,8.2,8.7,0,11.3}
	Coh_xs[369,370] = {10.9,7.5}

	Make/O/N=371 root:Packages:Neutron:Inc_xs
	Wave Inc_xs=root:Packages:Neutron:Inc_xs
	Inc_xs[0,18] = {80.26,80.27,2.05,0.14,0,1.6,0,0.92,0.46,0.78,0.0018,1.7,3,0.21,0.001,0,0.034,0.5,0.5}
	Inc_xs[19,37] = {5e-05,0.0008,0,0.004,0,0.0008,0.008,0,0.05,0,1.62,0.08,0,0.28,0,0.0082,0.004,0,0.001}
	Inc_xs[38,60] = {0,0.005,0.007,0,0.3,0,0,5.3,4.7,0.001,0.225,0,0,0,0.27,0.25,0.5,0.3,0.05,0,0,0.5,0}
	Inc_xs[61,85] = {0,0,4.5,2.87,0,1.5,0,3.3,0,5.08,0.5,5.07,1.83,0,0,5.93,0,0.4,0.4,0,0,0.3,0,4.8,5.2}
	Inc_xs[86,109] = {0,0,1.9,0,0,0.55,0.006,0.4,0.077,0,0,0.28,0,0,0.16,0.091,0.084,0.18,0,0,1.5,0,0,0.06}
	Inc_xs[110,135] = {0.32,0,0,0.05,0,0,0,0.1,0.15,0.05,0.01,0,0,0,Nan,0,0,0.5,0.5,0.5,0.06,0,0,0.5,0,0.15}
	Inc_xs[136,161] = {0.02,0,0.15,0,0,0,0.0024,0.04,0,0,0.5,0,0.5,0,0,0.5,0.4,0,0,Nan,0,Nan,0,0,0.3,0.093}
	Inc_xs[162,185] = {0,0,0.8,0,0,0,0.58,0.13,0.32,3.46,0,0,0,0.3,0,0.3,0,0,0.54,3.7e-05,0.55,0.022,0,0}
	Inc_xs[186,210] = {0.3,0,0.3,0,0.3,0,0,0,0.007,0.0003,0.001,0.09,0,0,0.52,0,0.008,0,0,0,0.31,0,0,0,0}
	Inc_xs[211,235] = {Nan,0,Nan,0,0,0,0.21,0.15,0,0,0,0.5,0,0.5,0,1.13,0.5,1.13,0.001,0,0,0,0,0.015,9.2}
	Inc_xs[236,264] = {0,55,0,5,0,0,0,1.3,39,0,143,0,137,0,0,0,2.5,3.1,1.3,151,0,0,25,0,394,0,0,0.004,54.4}
	Inc_xs[265,293] = {0,0,0,3,0,0.21,0,0.36,1.1,0,0,0,0.13,0,0,0.1,4,0,0,3.9,0,3.5,0,0,0.7,0.6,1.2,2.6,0}
	Inc_xs[294,320] = {0,0.1,0,0.14,0,0.01,0.5,0.011,1.63,0,0,0.3,0,0,0.9,0.5,1,0.3,0,0,0.3,0,0.5,0,0,0,Nan}
	Inc_xs[321,344] = {Nan,0.13,0,0,0,0.13,0,0,0.43,6.6,0,0,30,0,Nan,0,0,0.21,0.14,0.007,0.003,0,0,0.002}
	Inc_xs[345,368] = {0,0.0084,Nan,Nan,Nan,Nan,0,Nan,0,0.1,0.005,0.1,0,0.2,0,0.5,Nan,0,0.2,0,0,0.3,Nan,0}
	Inc_xs[369,370] = {0,0}

	Make/O/N=371 root:Packages:Neutron:Scatt_xs
	Wave Scatt_xs=root:Packages:Neutron:Scatt_xs
	Scatt_xs[0,16] = {82.02,82.03,7.64,3.03,1.34,6,1.34,1.37,0.97,1.4,7.63,5.24,3.1,5.77,5.551,5.559,4.84}
	Scatt_xs[17,32] = {11.51,11.53,5.21,4.232,4.232,4.2,4.29,4.018,2.628,2.695,5.7,1.88,3.28,3.71,4.03,1.93}
	Scatt_xs[33,49] = {3,1.503,2.167,2.12,2.78,2.64,3.312,1.026,0.988,3.1,1.52,1.1,16.8,21.8,1.19,0.683,77.9}
	Scatt_xs[50,67] = {1.5,0.421,1.96,2.01,1.6,1.2,2.83,2.9,1.42,0.8,0.25,1.6,0.019,23.5,4.35,3.05,3.2,4.65}
	Scatt_xs[68,86] = {3.4,4.8,5.1,7.8,5.09,3.49,2.54,3.042,8.15,2.6,2.15,11.62,2.2,12.42,1,28,5.6,18.5,26.1}
	Scatt_xs[87,103] = {0.99,9.2,9.5,0.017,8.03,5.2,14.5,4.131,3.42,4.48,7.46,4.57,4.5,6.83,7.89,5.23,8.6}
	Scatt_xs[104,122] = {12.6,9.1,4.7,7.2,8,5.5,8.3,0.1,18.7,8.65,8.5,7.03,5.05,5.9,5.96,5.84,7.68,Nan,Nan}
	Scatt_xs[123,142] = {Nan,Nan,6.6,8.2,6.8,6.7,7.1,6.25,6,4.04,7.4,6.42,7.7,6.46,5.1,9.7,6.9,8.4,3.8,6.255}
	Scatt_xs[143,161] = {5.71,6,5.81,6.5,4.83,7.1,5.44,5.69,6.3,6.6,Nan,Nan,Nan,Nan,Nan,144.8,4.483,4.6,4.48}
	Scatt_xs[162,181] = {7.5,7.5,4.6,5.1,2.1,7.5,4.99,7.3,2.5,6.5,3.1,3.7,4.4,5.6,5.1,12.4,7.1,5,2.62,3.65}
	Scatt_xs[182,200] = {2.57,4.892,4.5,4.8,4.8,4.42,5.6,4.63,5,5.29,4.14,4.48,3.9,4.1,3.64,4.32,3.5,1.8,0.52}
	Scatt_xs[201,220] = {8,3.18,3.88,4.36,4.55,3.81,Nan,Nan,Nan,Nan,Nan,Nan,Nan,Nan,Nan,Nan,3.9,3.38,1.6,7.6}
	Scatt_xs[221,239] = {4.08,3.2,3.03,6.4,2.94,9.66,8.5,9.66,2.94,4.23,5.64,2.94,2.84,2.66,16.6,7.5,80,1,30}
	Scatt_xs[240,262] = {9.5,4.1,3.5,21.3,39,1,39,1,200,25,3.1,11,9.2,8.6,9.8,180,13,13,66,5,1044,10,10.52}
	Scatt_xs[263,282] = {6.84,90.3,4.7,5,5.6,16,0.25,3.3,307,8.42,8.7,9.7,8.4,14.1,1.2,6.9,11.6,6.38,23.4,2.13}
	Scatt_xs[283,302] = {5.8,15.6,11.2,15,46.8,9.6,7.2,7.2,5.9,10.2,15,5.5,0.2,4.4,7.1,21.9,6.01,7,6.01,4.6}
	Scatt_xs[303,322] = {3,6.1,5.7,7.03,0.065,11.5,10.7,11.9,14.7,13,17,13,7.3,14.9,15.2,16.6,14,Nan,Nan,11.71}
	Scatt_xs[323,341] = {10,12.3,14,9.9,12.3,7.6,7.75,26.8,115,Nan,66,Nan,Nan,9.828,Nan,9.89,6.28,11.4,11.118}
	Scatt_xs[342,358] = {12.3,10.68,10.82,11.34,9.156,Nan,Nan,12.6,Nan,13,Nan,13.36,10.5,8.908,12.9,19.3,14}
	Scatt_xs[359,370] = {8.871,14.5,Nan,25,7.7,1.54,8.2,9,Nan,11.3,10.9,7.5}

	Make/O/N=371 root:Packages:Neutron:Abs_xs
	Wave Abs_xs=root:Packages:Neutron:Abs_xs
	Abs_xs[0,14] = {0.3326,0.3326,0.000519,0,0.00747,5333,0,70.5,940,0.0454,0.0076,767,3835,0.0055,0.0035}
	Abs_xs[15,27] = {0.00353,0.00137,1.9,1.91,2.4e-05,0.00019,1e-04,0.236,0.00016,0.0096,0.039,0.036,0.67}
	Abs_xs[28,42] = {0.046,0.53,0.063,0.05,0.19,0.0382,0.231,0.171,0.177,0.101,0.107,0.172,0.53,0.54,0.54}
	Abs_xs[43,60] = {0.227,0.15,33.5,44.1,0.433,0.675,5.2,0.8,0.66,2.1,2.1,35,1.46,0.43,0.41,0.68,6.2,0.88}
	Abs_xs[61,78] = {0.74,1.09,27.5,6.09,0.59,1.7,7.84,2.2,0.179,5.08,60,4.9,3.05,15.8,0.76,18.1,0.36,13.3}
	Abs_xs[79,95] = {2.56,2.25,2.59,2.48,1.28,37.18,4.49,4.6,2.9,2.5,14.5,1.52,3.78,4.5,2.17,1.11,0.93}
	Abs_xs[96,114] = {0.62,6.8,1.1,0.092,2.75,2.18,3.61,2.2,3,0.8,15.1,0.4,0.16,4.5,11.7,51.8,85,42,0.43}
	Abs_xs[115,132] = {0.61,0.044,6.9,11,2.7,25,6.4,11.8,29,185,0.113,0.003,0.38,0.48,0.12,1.28,0.87,1.04}
	Abs_xs[133,148] = {16,0.058,1.28,0.185,0.011,1.17,0.22,0.0499,0.0229,1.15,2.48,0.019,0.015,13.1,0.5,2.5}
	Abs_xs[149,167] = {0.127,0.4,20,2.56,0.28,Nan,6.9,4.8,3.3,1.17,0.31,144.8,6.9,3.4,0.6,20,0.304,8.55,0.226}
	Abs_xs[168,187] = {63.3,37.6,91,2520,1,1.1,11,24,2.2,20600,0.34,0.075,193.8,12,202,0.626,1,0.114,30,0.14}
	Abs_xs[188,205] = {2.3,0.22,2.2,0.14,0.18,0.133,4.91,5.75,3.8,4.7,2.3,3.4,418,6.8,1.55,1.04,0.215,0.29}
	Abs_xs[206,226] = {6.15,23.9,165,3.5,Nan,21,Nan,85,0.45,0.265,0.26,29,1.1,30,7,2,5.8,0.68,3.6,0.27,8.97}
	Abs_xs[227,245] = {57,8.93,0.63,7.3,1.1,0.57,0.95,11.5,50.5,18.7,337,3.6,42,1.4,2.5,1.2,168.4,5922,0.7}
	Abs_xs[246,263] = {57,2.4,42080,104,206,8.4,4530,9100,312,49700,735,85,61100,1.5,2.59e+05,2.2,0.77,23.4}
	Abs_xs[264,283] = {994,33,43,56,600,194,124,2840,64.7,159,19,13,19.6,659,2.74,5.8,100,34.8,2230,11.4}
	Abs_xs[284,302] = {48.6,0.8,17.1,69.4,2.85,74,21,2065,104.1,561,23.5,373,84,41,13.04,20.6,563,20.5,18.3}
	Abs_xs[303,322] = {30,20.7,10.1,1.7,37.9,89.7,112,76.4,16,3000,80,320,4.7,25,13.1,2,425,954,111,10.3}
	Abs_xs[323,340] = {152,10,1.44,27.5,0.72,3.66,98.65,372.3,3080,2,2150,Nan,7.8,4.89,0.43,3.43,11.4,0.104}
	Abs_xs[341,356] = {0.171,0.65,0.03,0.699,0.00048,0.0338,Nan,Nan,Nan,Nan,12.8,Nan,7.37,200.6,7.57,574.7}
	Abs_xs[357,370] = {100.1,680.9,2.68,175.9,Nan,558,1017.3,289.6,18.5,75.3,Nan,16.2,1.36,3}
End


