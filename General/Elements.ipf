#pragma rtGlobals=2		// Use modern global access method.
#pragma IgorVersion = 4.0
#pragma version = 1.80
#pragma ModuleName=elements
#include "MaterialsLocate"								// used to find the path to the materials files
Constant MIN_LINE_SEPARATION_FRACTION = 0.15	// you can over ride this in your main procedure window.
//StrConstant ELEMENT_Symbols = "H;He;Li;Be;B;C;N;O;F;Ne;Na;Mg;Al;Si;P;S;Cl;Ar;K;Ca;Sc;Ti;V;Cr;Mn;Fe;Co;Ni;Cu;Zn;Ga;Ge;As;Se;Br;Kr;Rb;Sr;Y;Zr;Nb;Mo;Tc;Ru;Rh;Pd;Ag;Cd;In;Sn;Sb;Te;I;Xe;Cs;Ba;La;Ce;Pr;Nd;Pm;Sm;Eu;Gd;Tb;Dy;Ho;Er;Tm;Yb;Lu;Hf;Ta;W;Re;Os;Ir;Pt;Au;Hg;Tl;Pb;Bi;Po;At;Rn;Fr;Ra;Ac;Th;Pa;U;Np;Pu;Am;Cm;Bk;Cf;Es;Fm;Md;No;Lr;Rf;Db;Sg;Bh;Hs;Mt;Ds;Rg;Cn;;Fl;;Lv"
Constant ELEMENT_Zmax = 116

//	Sept 14, 2005
//		updated amu using NIST values
//		also changed Formula2amu() and BindingEnergy() so they will not prompt if called by another function
//
//	Oct 29, 2005
//		changed so that Function/T ChemicalFormulaToList("") returns an empty string
//
//	Feb 7, 2007
//		changed how it looks for the 'materials' folder
//
//	Mar 18, 2007
//		made the 'materials' folder work properly
//
//	Jan 3, 2007
//		suppress some of the alerts in ChemicalFormulaToList()
//
//	Mar 7, 2011
//		modified ChemicalFormulaToList(), it now accepts parenthesis sub-formulas
//
//	Mar 8, 2011		1.64
//		fixed error in  ChemicalFormulaToList(), it now accepts fractional atoms correctly
//		fixed LookUpMaterial() too, similar problem
//
//	Jul 20, 2011		1.65
//		fixed error in  ChemicalFormulaToList(), it now deal with compound formulats, e.g. "TiO2 2(P 3(CO))"
//
//	Oct 16, 2012		1.66
//		Added emission lines, MakeFullEmissionInfo(), uses EmissionLineStruct (MakeAllElementLines() is deprecated)
//
//	Jun 12, 2013		1.71
//		Added Debye Temperatures (K)
//
//	Jun 13, 2013		1.70
//		Added elements up to 116, and made the string constant ELEMENT_Symbols and ELEMENT_Zmax, and removed references to "root:Packages:Elements:symbols"
//
//	Jun 19, 2013		1.73
//		fixed up LookUpMaterial(), added niceFloatStr(), suppressed use of global symbols everywhere, only use ELEMENT_Symbols
//
//	Jun 19, 2013		1.74
//		changed LookUpMaterial() to read a mtl file (materials) that is in an xml format
//
//	Jun 25, 2013		1.75
//		added FindBestElementFromEmissionLine(), to identify the element from a line
//
//	Jun 28, 2014		1.77
//		modified FindBestElementFromEmissionLine(), added an excitation energy
//
//	Sep 30, 2015		1.78
//		moved the constant string "ELEMENT_Symbols" from here to Utility_JZT.ipf
//
//	Nov 14, 2015		1.79
//		moved isotope data to data/isotopes.xml
//
//	Nov 16, 2015		1.80
//		moved all the other element data to isotope data to data/elementData.xml

Menu "Analysis"
      Submenu "Element"
		"Element Data", ElementData("","")
		"Atomic mass of molecule",Formula2amu("")
	End

        Submenu "X-ray"
		"Binding Energies",BindingEnergy("","")
		"Emission Energies",EmissionEnergies("","")
        End
End


Function ElementData(symb,property)
	String symb				// atomic symbol
	String property				// item to return

	Variable topLevel = ItemsInList(GetRTStackInfo(0))<2
	Variable val=NaN
	String keys = "Z;amu;density;valence;state;firstIon;heatFusion;heatVapor;thermConduc;specificHeat;"
	keys += "atomRadius;meltPoint;boilPoint;covRadius;elecConduc;electroneg;atomVol;name"
	keys = LowerStr(keys)
	property = LowerStr(property)
	Variable ip = WhichListItem(property,keys)		// index to property
	ip = (ip>=0) ? ip+1 : NaN
	Variable Z = element2Z(symb)
	String propertys = "atomic number;atomic mass;density;valence;state;first ionization energy;heat of fusion;"
	propertys += "heat of vaporization;thermal conductivity;specific heat;atomic radius;melting point;"
	propertys += "boiling point;covalent radius;electrical conductivity;electronegativity;atomic volume;DebyeT (K);name"
	if ((numtype(Z) || numtype(ip)) && topLevel)
		Prompt ip, "property to return",popup, propertys
		Prompt symb, "atomic symbol", popup, ELEMENT_Symbols
		DoPrompt "select", symb,ip
		if (V_flag)
			return NaN
		endif
		Z = element2Z(symb)
	endif
	if (numtype(Z) || numtype(ip))
		return NaN
	endif
	String out=""									// optional string out
	if (ip==1)
		val = Z
	elseif (ip==2)
		val = Element_amu(Z)
	elseif (ip==3)
		val = Element_density(Z)
	elseif (ip==4)
		val = Element_valence(Z)
	elseif (ip==5)
		val = WhichListItem(Element_state(Z),"s;l;g")
		out = Element_state(Z)
	elseif (ip==6)
		val = Element_firstIon(Z)
	elseif (ip==7)
		val = Element_heatFusion(Z)
	elseif (ip==8)
		val = Element_heatVapor(Z)
	elseif (ip==9)
		val = Element_thermConduc(Z)
	elseif (ip==10)
		val = Element_specificHeat(Z)
	elseif (ip==11)
		val = Element_atomRadius(Z)
	elseif (ip==12)
		val = Element_meltPoint(Z)
	elseif (ip==13)
		val = Element_boilPoint(Z)
	elseif (ip==14)
		val = Element_covRadius(Z)
	elseif (ip==15)
		val = Element_elecConduc(Z)
	elseif (ip==16)
		val = Element_electroneg(Z)
	elseif (ip==17)
		val = Element_atomVol(Z)
	elseif (ip==18)
		val = Element_DebyeT(Z)
	elseif (ip==19)
		out = Element_names(Z)
	else
		return NaN
	endif
	if (topLevel)
		printf "for %s, the %s is %s\r",symb,StringFromList(ip-1,propertys),SelectString(strlen(out),num2str(val),"'"+out+"'")
	endif
	return val
End


Function element2Z(symb)		// returns Z for an atomic symbol (NOT case sensitive)
	String symb				// atomic symbol
	symb[0,0] = UpperStr(symb[0,0])		// ensure first char is upper case
	symb[1,1] = LowerStr(symb[1,1])		// and second character is lower
	Variable iz = WhichListItem(symb,ELEMENT_Symbols)+1
	if (!(iz>0))
		iz = str2num(symb)
	endif
	return ((iz>0) ? iz : NaN)
End


Function Formula2amu(formula)		// returns mol. wt of a chemical formula
	String formula		// chemical formula, Elements start with a capital letter.
	Variable topLevel = ItemsInList(GetRTStackInfo(0))<2
	if (strlen(formula)<1)
		if (!topLevel)
			return NaN
		endif
		Prompt formula, "chemical formula"
		DoPrompt "chemical formula",formula
		if (V_flag)
			return NaN
		endif
	endif
	String formula_In = formula
	if (strsearch(formula,";",0)<0 && strsearch(formula,",",0)<0)
		formula = ChemicalFormulaToList(formula)
	endif
	String symb		// element symbol
	String item			// one element with its number of atoms
	Variable i,Natoms, molWt=0
	for (i=0;strlen(StringFromList(i,formula));i+=1)
		item = StringFromList(i,formula)
		symb = StringFromList(0,item,",")
		Natoms = str2num(StringFromList(1,item,","))
		Natoms = numtype(Natoms) ? 1 : Natoms
		molWt += Natoms * Element_amu(element2Z(symb))
	endfor
	if (topLevel)
		printf  "   atomic mass of '%s' is %g (amu)\r",formula_In,molWt
	endif
	return molWt
End


Function/T ChemicalFormulaToList(formula,[force1s])	// parse a chemical formula into a list
	// understands formulas of type:    "TiO2 2(P 3(CO))".  Note atomic symbols must have proper Upper/Lower case
	String formula		// chemical formula, Elements start with a capital letter.
	Variable force1s	// forces a 1 for single elements in output list
	force1s = ParamIsDefault(force1s) ? 0 : force1s

	formula = ReplaceString("\t",formula," ")
	formula = ReplaceString(",",formula," ")
	formula = ReplaceString("  ",formula," ")
	Variable N=strlen(formula)
	Variable ascii					// value of one ascii character
	Variable i=0					// index into formula
	String symb					// atomic symbol
	Variable num					// number of atoms for symb (may be a fraction)
	Variable multiplier			// multiplier for sub formula
	String list=""					// parsed formula as a list of the form "TiO2" -> "Ti;O,2"
	i = 0
	do
		ascii = char2num(formula[i])		// all done, no more ascii
		if (numtype(ascii))
			break
		elseif (ascii==32)					// skip spaces
			i += 1
			continue
		endif

		multiplier = str2num(formula[i,N])
		if (numtype(multiplier)==0)		// a multiplier with a sub formula
			Variable i0,i1
			i0 = strsearch(formula,"(",i)
			i1 = FindMatchingBracket(formula,i0)
			if (i1<=i0)
				if (ItemsInList(GetRTStackInfo(0))<2)
					DoAlert 0, "ChemicalFormulaToList(), cannot find matching parenthsis in  '"+formula+"'"
				endif
				return ""
			endif

			// printf "sub formula:      %g(%s)\r",multiplier,formula[i0+1,i1-1]
			String subList = ChemicalFormulaToList(formula[i0+1,i1-1],force1s=1)
			subLIst = multiplyElemetnList(subList,multiplier)	// increment subList by multiplier
			list += subList
			i = i1+1
			continue
		endif

		if (!(ascii>=char2num("A") && char2num("Z")>=ascii))	// formula is empty, all done (all elements start with a capital letter)
			if (i>0 && ItemsInList(GetRTStackInfo(0))<2)
				DoAlert 0, "ChemicalFormulaToList(), atomic symbols must start with a Capital letter, could not parse '"+formula+"'"
			endif
			return ""
		endif
		symb = formula[i]
		i += 1									// move pointer ahead one
		ascii = char2num(formula[i])
		if (char2num("a")<=ascii && ascii<=char2num("z"))	// next char is a lower case letter, make it part of symb
			symb += formula[i]				// optional second character of atomic symbol
			i += 1								// move pointer ahead one
		endif

		ascii = char2num(formula[i])		// after symbol will either be a number or another symbol, or end of string
		if (ascii>=char2num("A") && ascii<=char2num("Z"))	// found a capital letter
			num = 1
		elseif (i>=N)							// found end of string
			num = 1
		elseif (ascii==32 || ascii==9)		// found whitespace
			num = 1
			i += 1
		else
			num = str2num(formula[i,N])	// assume that it must be a number
			if  (numtype(num))
				if (ItemsInList(GetRTStackInfo(0))<2)
					DoAlert 0, "ChemicalFormulaToList(), found invalid number, could not parse '"+formula+"'"
				endif
				return ""
			endif
			i0 = i-1							// got the number after a symbol, find char after the number
			do
				i += 1
			while (strsearch("+-.0123456789", formula[i],0)>=0)
		endif
		if (num==1 && !force1s)			// add the number of this element to the list
			list += symb+";"
		else
			list += symb+","+num2str(num)+";"
		endif
	while(i<N)
	return list
End
//
Static Function/T multiplyElemetnList(in,mult)
	String in
	Variable mult
	if (strlen(in)<1 || !(mult>0) && numtype(mult)!=0)
		return ""
	endif

	String item,symb,out=""
	Variable i,num
	for (i=0;i<ItemsInLIst(in);i+=1)
		item = StringFromList(i,in)
		symb = StringFromList(0,item,",")
		num = str2num(StringFromList(1,item,","))*mult
		out += symb+","+num2str(num)+";"
	endfor
	return out
End
//
Static Function FindMatchingBracket(str,istart)
	String str
	Variable istart			// location of b0 in str

	String b0=str[istart]		// starting bracket, can be "(" "[", or "{"
	String bOpen="([{", bClose=")]}"
	Variable ib=strsearch(bOpen,b0,0), iclose
	if (ib<0)
		return -1
	endif
	String b1=bClose[ib]
	do
		iclose = strsearch(str,b1,istart)
		if (strsearch(str[istart+1,iclose],b0,0)<1)
			break
		endif
		istart = iclose+1
	while(1)
	return iclose
End


Function/T LookUpMaterial(material)
	String material				// name of the material, look in file materials/material.mtl

	Variable refNum
	String list					// results

	PathInfo materials
	if (strlen(S_path)<1)							// make it if it does not exist
		NewPath/Z materials, ParseFilePath(1,FunctionPath("MaterialsAreHere"),":",1,0)
	endif
	PathInfo materials
	if (strlen(S_path)<1)							// make it if it does not exist
		PathInfo Igor
		NewPath/Z materials, S_path+"User Procedures:materials"
	endif
	PathInfo materials
	if (strsearch(material,".mtl",0,2)<0)
		material += ".mtl"
	endif
	Open/P=igor/R/T="TEXT"/Z refNum S_path+material
	if (V_flag)					// file not found
		// if not found in igor, try in path=home
		Open/P=home/R/T="TEXT"/Z refNum material
		if (V_flag)				// file not found
		return ""				// give up
		endif
	endif
	FStatus refNum
	String buf=PadString("",V_logEOF,0x20)
	FBinRead refNum, buf
	Close refNum

	if (strsearch(buf, "<?xml",0)>=0)
		list = ProcessMTLfileContentsXML(buf)
	else
		list = ProcessOldMTLfileContents(buf)
	endif
	return list
End
//
Static Function/T ProcessMTLfileContentsXML(buf)
	String buf						// contents of a new mtl file in xml format

	buf = XMLremoveComments(buf)
	String mix = XMLtagContents("chemical_mixture",buf)

	String name = StringByKey("name", XMLattibutes2KeyList("chemical_mixture",buf),"="), list=""
	if (strlen(name))
		list = ReplaceStringByKey("name",list,name,"=")
	endif

	Variable density=str2num(XMLtagContents("density",mix))
	String unit = StringByKey("unit", XMLattibutes2KeyList("density",mix),"=")
	unit = SelectString(strlen(unit),"g/cm^3",unit)	// no unit specified, set to default
	if (!StringMatch(unit,"g/cm^3"))					// only understand density in g/cm^3
		printf "ERROR -- Invalid units for density = '%s'\r",unit
		density = NaN
	endif
	if (numtype(density)==0 && density>0)
		list = ReplaceNumberByKey("density",list,density,"=")
	endif

	String part, symbol, formula=""
	Variable i=0, fraction, Z
	do
		part = XMLattibutes2KeyList("part",buf,occurance=i)
		symbol = StringByKey("symbol",part,"=")
		if (strlen(symbol)<1)
			Z = NumberByKey("Z",part,"=")
			symbol = StringFromList(Z-1,ELEMENT_Symbols)
		endif
		fraction = NumberByKey("fraction",part,"=")
		fraction = numtype(fraction) ? 1 : fraction
		if (strlen(symbol))
			formula += symbol
			if (fraction!=1)
				formula += elements#niceFloatStr(fraction)
			endif
		endif
		i += 1
	while (strlen(symbol))
	if (strlen(formula))
		list += "formula="+formula+";"
	endif
	return list
End
//
Static Function/T ProcessOldMTLfileContents(buf)
	String buf						// contents of an old mtl file
		//	CdWO4
		//	density, ncomp, zed, fract, zed, fract
		//	7.9,3,48,1.0,74,1.0,8,4.0
	buf = ReplaceString("\r\n", buf,"\n")
	buf = ReplaceString("\n\r", buf,"\n")
	buf = ReplaceString("\r", buf,"\n")

	Variable i1, i0 = strsearch(buf,"\n",0)+1	// position at start of second line
	i0 = strsearch(buf,"\n",i0)+1			// position at start of third line
	i1 = strsearch(buf,"\n",i0)-1			// position before terminator at end of third line
	String line = buf[i0,i1]

	line = ReplaceString("\t",line,",")	// replace all tabs with commas
	line = ZapControlCodes(line)
	line = ReplaceString(" ",line,",")	// replace all spaces with commas
	line = ReplaceString(";",line,",")	// replace all semi-colons with commas
	do
		line = ReplaceString(",,",line,",")	// change all multiple commas to one commas
	while (strsearch(line, ",,",0)>=0)
	if (strsearch(line, ",",0)==0)			// remove a possibly leading comma
		line = line[1,Inf]
	endif
	// line should now have each item sepearated with only one comma, and no white space

	Variable density,N,i,Z,f
	sscanf line, "%g,%g",density,N
	sprintf list, "density=%g",density

	String formula="", list=""
	for (i=0;i<2*N;i+=2)
		Z = str2num(StringFromList(i+2,line,","))
		f = str2num(StringFromList(i+3,line,","))
		if (numtype(Z) || numtype(f))
			return ""
		endif
		formula += StringFromList(Z-1,ELEMENT_Symbols)
		if (f!=1)
			formula += niceFloatStr(f)
		endif
	endfor
	if (strlen(formula))
		list += ";formula="+formula+";"
	endif
	return list
End
//
Static Function/T niceFloatStr(v)	// format float string to minimal decimal size (with NO E-5 type notation)
	Variable v

	if (numtype(v))
		return num2str(v)
	elseif (v==0)
		return "0"
	endif
	Variable av = abs(v)			// av is now >0, not negative, not 0

	String fmt, str=""
	Variable after = max(placesOfPrecision(av)-floor(log(av))-1, 0)
	sprintf fmt,"%%.%df",after
	sprintf str,fmt,v
	return str
End


// ========================================================================
// ========================================================================
// ========================================================================


Function Element_amu(Z)
	Variable Z
	String amuList
	amuList   = "1.00794;4.002602;6.941;9.012182;10.811;12.0107;14.0067;15.9994;18.9984032;20.1797;22.98977;24.305;"
	amuList += "26.981538;28.0855;30.973761;32.065;35.453;39.948;39.0983;40.078;44.95591;47.867;50.9415;51.9961;"
	amuList += "54.938049;55.845;58.9332;58.6934;63.546;65.409;69.723;72.64;74.9216;78.96;79.904;83.798;85.4678;87.62;"
	amuList += "88.90585;91.224;92.90638;95.94;98;101.07;102.9055;106.42;107.8682;112.411;114.818;118.71;121.76;"
	amuList += "127.6;126.90447;131.293;132.90545;137.327;138.9055;140.116;140.90765;144.24;145;150.36;151.964;157.25;"
	amuList += "158.92534;162.5;164.93032;167.259;168.93421;173.04;174.967;178.49;180.9479;183.84;186.207;190.23;"
	amuList += "192.217;195.078;196.96655;200.59;204.3833;207.2;208.98038;209;210;222;223;226;227;232.0381;231.03588;"
	amuList += "238.02891;237;244;243;247;247;251;252;257;258;259;262;261;262;266;264;277;268;281;280;285;284;289;288;293;" 
	return (str2num(StringFromList(Z-1,amuList)))
End

Function Element_density(Z)
	Variable Z
	String densityList
	densityList   = "0.071;0.126;0.533;1.845;2.34;2.26;0.001165;0.001331;0.001696;0.0008391;0.969;1.735;2.6941;2.32;"
	densityList += "1.82;2.07;0.003214;0.00166;0.86;1.55;2.98;4.53;6.1;7.18;7.43;7.86;8.9;8.876;8.94;7.112;5.877;"
	densityList += "5.307;5.72;4.78;3.11;0.003484;1.529;2.54;4.456;6.494;8.55;10.2;11.48;12.39;12.39;12;10.48;"
	densityList += "8.63;7.3;7.3;6.679;6.23;4.92;0.005458;1.87;3.5;6.127;6.637;6.761;6.994;7.2;7.51;5.228;7.877;"
	densityList += "8.214;8.525;8.769;9.039;9.294;6.953;9.811;13.29;16.624;19.3;20.98;22.53;22.39;21.41;18.85;"
	densityList += "13.522;11.83;11.33;9.73;9.3;1;0.00923;1;5;10.05;11.7;15.34;18.92;"
	densityList += "20.25;19.84;11.7;13.51;14.;"
	return (str2num(StringFromList(Z-1,densityList)))
End

Function Element_valence(Z)
	Variable Z
	String valence
	valence   = "1;0;1;2;3;4;3;2;1;0;1;2;3;4;5;6;1;0;1;2;3;4;5;3;2;"
	valence += "3;2;2;2;2;3;4;3;4;1;0;1;2;3;4;5;6;7;3;3;2;1;2;3;4;"
	valence += "3;4;1;0;1;2;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;4;5;6;7;"
	valence += "4;4;4;3;2;1;2;3;2;1;0;1;2;3;4;5;6;5;4;3;"
	return (str2num(StringFromList(Z-1,valence)))
End

Function/T Element_state(Z)
	Variable Z
	String state
	state   = "g;g;s;s;s;s;g;g;g;g;s;s;s;s;s;s;g;g;s;s;s;s;s;s;s;"
	state += "s;s;s;s;s;l;s;s;s;l;g;s;s;s;s;s;s;s;s;s;s;s;s;s;s;"
	state += "s;s;s;g;l;s;s;s;s;s;s;s;s;s;s;s;s;s;s;s;s;s;s;s;s;"
	state += "s;s;s;s;l;s;s;s;s;s;g;l;s;s;s;s;s;s;s;s;"
	return (StringFromList(Z-1,state))
End

Function Element_firstIon(Z)
	Variable Z
	String firstIon
	firstIon   = "13.598;24.587;5.392;9.322;8.298;11.26;14.534;13.618;17.422;21.564;"
	firstIon += "5.139;7.646;5.986;8.151;10.486;10.36;12.967;15.759;4.341;6.113;"
	firstIon += "6.54;6.82;6.74;6.766;7.435;7.87;7.86;7.635;7.726;9.394;"
	firstIon += "5.999;7.899;9.81;9.752;11.814;13.999;4.177;5.695;6.38;6.84;"
	firstIon += "6.88;7.099;7.28;7.37;7.46;8.34;7.576;8.993;5.786;7.344;"
	firstIon += "8.641;9.009;10.451;12.13;3.894;5.212;5.58;5.47;5.42;5.49;"
	firstIon += "5.55;5.63;5.67;6.15;5.86;5.93;6.02;6.101;6.184;6.254;"
	firstIon += "5.43;6.65;7.89;7.98;7.88;8.7;9.1;9;9.225;10.437;"
	firstIon += "6.108;7.416;7.289;8.42;;10.748;0;5.279;5.17;6.08;"
	firstIon += "5.88;6.05;6.19;6.06;6;6.02;6.23;6.3;6.42;6.5;"
	return (str2num(StringFromList(Z-1,firstIon)))
End

Function Element_heatFusion(Z)
	Variable Z
	String heatFusion
	heatFusion   = "0.0585;0.021;3;11.71;22.6;;0.36;0.222;0.26;0.34;"
	heatFusion += "2.601;8.95;10.7;50.2;0.63;1.73;3.21;1.188;2.33;8.53;"
	heatFusion += "16.11;18.6;20.8;20;14.64;13.8;16.19;17.2;13.14;7.38;"
	heatFusion += "5.59;31.8;27.7;5.54;5.286;1.638;2.34;8.2;17.5;21;"
	heatFusion += "26.9;36;23;25.52;21.76;16.74;11.3;6.07;3.26;7.2;"
	heatFusion += "19.83;17.49;7.76;2.3;2.092;8.01;11.3;9.2;10.04;10.88;"
	heatFusion += " ;11.09;9.21;10.46;15.48;;11.06;17.15;16.8;7.7;"
	heatFusion += "18.6;21.76;36;35.4;33.05;29.29;26.36;19.66;12.36;2.292;"
	heatFusion += "4.27;4.77;11;13;12;2.9;2.1;8.37;;15.65;"
	heatFusion += " ;15.48;"
	return (str2num(StringFromList(Z-1,heatFusion)))
End

Function Element_heatVapor(Z)
	Variable Z
	String heatVapor
	heatVapor   = "0.4581;0.084;147.1;297;507.8;715;2.7928;3.4109;3.2698;1.77;"
	heatVapor += "98.01;127.6;290.8;359;12.4;10;10.2;6.506;76.9;154.67;"
	heatVapor += "304.8;425.2;446.7;339.5;219.74;349.5;373.3;377.5;300.5;115.3;"
	heatVapor += "256.06;334.3;32.4;26.32;14.725;9.029;69.2;136.9;363.3;590.5;"
	heatVapor += "690.1;590.4;502;567.77;495.39;393.3;250.63;99.87;226.35;290.37;"
	heatVapor += "67.97;50.63;20.9;12.64;67.74;140.2;399.57;313.8;332.63;283.68;"
	heatVapor += " ;191.63;175.73;311.71;;230;251.04;292.88;191;128;"
	heatVapor += "355;661.07;737;422.58;707.1;627.6;563.58;510.45;324.43;59.3;"
	heatVapor += "162.09;177.9;179;120;30;16.4;64;136.82;;543.92;"
	heatVapor += " ;422.58;"
	return (str2num(StringFromList(Z-1,heatVapor)))
End

Function Element_thermConduc(Z)
	Variable Z
	String thermConduc
	thermConduc   = "0.1815;0.152;84.7;200;27;155;0.02598;0.2674;0.0279;0.0493;"
	thermConduc += "141;156;237;148;0.235;0.269;0.0089;0.0177;102.5;200;"
	thermConduc += "15.8;21.9;30.7;93.7;7.82;80.2;100;90.7;401;116;"
	thermConduc += "40.6;59.9;50;2.04;0.122;0.00949;58.2;35.3;17.2;22.7;"
	thermConduc += "53.7;138;50.6;117;150;71.8;429;96.8;81.6;66.6;"
	thermConduc += "24.3;2.35;0.449;0.00569;35.9;18.4;13.5;11.4;12.5;16.5;"
	thermConduc += "17.9;13.3;13.9;10.6;11.1;10.7;16.2;14.3;16.8;34.9;"
	thermConduc += "16.4;23;57.5;174;47.9;87.6;147;71.6;317;8.34;"
	thermConduc += "46.1;35.3;7.87;20;1.7;0.00364;15;18.6;12;54;"
	thermConduc += "47;27.6;6.3;6.74;10;10;10;10;10;10;"
	return (str2num(StringFromList(Z-1,thermConduc)))
End

Function Element_specificHeat(Z)
	Variable Z
	String specificHeat
	specificHeat   = "14.304;5.193;3.582;1.825;1.026;0.709;1.042;0.92;0.824;1.03;"
	specificHeat += "1.23;1.02;0.9;0.70;0.769;0.71;0.48;0.52;0.757;0.647;"
	specificHeat += "0.568;0.523;0.489;0.449;0.48;0.449;0.421;0.444;0.385;0.388;"
	specificHeat += "0.371;0.32;0.33;0.32;0.226;0.248;0.363;0.3;0.3;0.278;"
	specificHeat += "0.265;0.25;0.24;0.238;0.242;0.244;0.232;0.233;0.233;0.228;"
	specificHeat += "0.207;0.202;0.145;0.158;0.24;0.204;0.19;0.19;0.193;0.19;"
	specificHeat += " ;0.197;0.182;0.236;0.18;0.173;0.165;0.168;0.16;0.155;"
	specificHeat += "0.15;0.14;0.14;0.13;0.137;0.13;0.13;0.13;0.128;0.140;"
	specificHeat += "0.129;0.129;0.122;;;0.094;;0.094;0.12;0.113;"
	specificHeat += " ;0.12;;0.13;"
	return (str2num(StringFromList(Z-1,specificHeat)))
End

Function Element_atomRadius(Z)
	Variable Z
	String atomRadius
	atomRadius   = "2.08;;1.55;1.12;0.98;0.91;0.92;0.65;0.57;0.51;1.9;1.6;1.43;1.32;1.28;"
	atomRadius += "1.27;0.97;0.88;2.35;1.97;1.62;1.45;1.34;1.3;1.35;1.26;1.25;1.24;1.28;1.38;"
	atomRadius += "1.41;1.37;1.39;1.4;1.12;1.03;2.48;2.15;1.78;1.6;1.46;1.39;1.36;1.34;1.34;"
	atomRadius += "1.37;1.44;1.71;1.66;1.62;1.59;1.42;1.32;1.24;2.67;2.22;1.38;1.81;1.82;1.82;"
	atomRadius += " ;1.81;1.99;1.8;1.8;1.8;1.79;1.78;1.77;1.94;1.75;1.67;1.49;1.41;1.37;"
	atomRadius += "1.35;1.36;1.39;1.46;1.6;1.71;1.75;1.7;1.67;1.45;1.34;2.7;2.33;1.88;1.8;1.63;1.56;"
	return (str2num(StringFromList(Z-1,atomRadius)))
End

Function Element_meltPoint(Z)
	Variable Z
	String meltPoint
	meltPoint   = "13.81;0.95;453.7;1560;2365;3825;63.15;54.8;53.55;24.55;"
	meltPoint += "371;922;933.5;1683;317.3;392.2;172.17;83.95;336.8;1112;"
	meltPoint += "1814;1945;2163;2130;1518;1808;1768;1726;1356.6;692.73;"
	meltPoint += "302.92;1211.5;1090;494;265.95;116;312.63;1042;1795;2128;"
	meltPoint += "2742;2896;2477;2610;2236;1825;1235.08;594.26;429.78;505.12;"
	meltPoint += "903.91;722.72;386.7;161.39;301.54;1002;1191;1071;1204;1294;"
	meltPoint += "1315;1347;1095;1585;1629;1685;1747;1802;1818;1092;"
	meltPoint += "1936;2504;3293;3695;3455;3300;2720;2042.1;1337.58;234.31;"
	meltPoint += "577;600.65;544.59;527;575;202;300;973;1324;2028;"
	meltPoint += "1845;1408;912;913;1449;1620;;1170;1130;1800;"
	return (str2num(StringFromList(Z-1,meltPoint)))
End

Function Element_boilPoint(Z)
	Variable Z
	String boilPoint
	boilPoint   = "20.28;4.216;1615;3243;4275;5100;77.344;90.188;85;27.1;"
	boilPoint += "1156;1380;2740;2630;553;717.82;239.18;87.45;1033;1757;"
	boilPoint += "3109;3560;3650;2945;2335;3023;3143;3005;2840;1180;"
	boilPoint += "2478;3107;876;958;331.85;120.85;961;1655;3611;4682;"
	boilPoint += "5015;4912;4538;4425;3970;3240;2436;1040;2350;2876;"
	boilPoint += "1860;1261;457.5;165.1;944;2078;3737;3715;3785;3347;"
	boilPoint += "3273;2067;1800;3545;3500;2840;2968;3140;2223;1469;"
	boilPoint += "3668;4875;5730;5825;5870;5300;4700;4100;3130;629.88;"
	boilPoint += "1746;2023;1837;;610;211.4;950;1413;3470;5060;"
	boilPoint += "4300;4407;4175;3505;2880;"
	return (str2num(StringFromList(Z-1,boilPoint)))
End

Function Element_covRadius(Z)
	Variable Z
	String covRadius
	covRadius   = "0.32;0.93;1.23;0.9;0.82;0.77;0.75;0.73;0.72;0.71;1.54;1.36;1.18;1.11;1.06;"
	covRadius += "1.02;0.99;0.98;2.03;1.74;1.44;1.32;1.22;1.18;1.17;1.17;1.16;1.15;1.17;1.25;"
	covRadius += "1.26;1.22;1.2;1.16;1.14;1.89;2.16;1.91;1.62;1.45;1.34;1.3;1.27;1.25;1.25;"
	covRadius += "1.28;1.34;1.41;1.44;1.41;1.4;1.36;1.33;1.31;2.35;1.98;1.25;1.65;1.65;1.64;"
	covRadius += "1.63;1.62;1.85;1.61;1.59;1.59;1.58;1.57;1.56;1.7;1.56;1.44;1.34;1.3;1.28;"
	covRadius += "1.26;1.27;1.3;1.34;1.49;1.48;1.47;1.46;1.53;1.47;;;;;1.65;2.0;1.96;"
	return (str2num(StringFromList(Z-1,covRadius)))
End

Function Element_elecConduc(Z)
	Variable Z
	String elecConduc
	elecConduc   = " ;;11.7;25;5.E-12;0.07;;;;;20.1;22.4;37.7;0.0004;1.E-16;"
	elecConduc += "5.E-16;;;16.4;31.3;1.5;2.6;4;7.9;0.5;11.2;17.9;14.6;60.7;16.9;"
	elecConduc += "1.8;3.E-6;3.8;8;1.E-16;;47.8;5;1.8;2.3;6.6;17.3;0.001;14.9;23;"
	elecConduc += "10;62.9;14.7;3.4;8.7;2.6;0.0002;1.E-11;;5.3;2.8;1.9;1.4;1.5;1.6;"
	elecConduc += "2;1.1;1.1;0.8;0.9;1.1;1.1;1.2;1.3;3.7;1.5;3.4;8.1;18.2;5.8;"
	elecConduc += "12.3;21.3;9.4;48.8;1;5.6;4.8;0.9;0.7;;;;1;;7.1;"
	return (str2num(StringFromList(Z-1,elecConduc)))
End

Function Element_electroneg(Z)
	Variable Z
	String electroneg
	electroneg   = "2.1;0;0.98;1.57;2.04;2.55;3.04;3.44;3.98;0;0.93;1.31;1.61;1.9;2.19;"
	electroneg += "2.58;3.16;0;0.82;1;1.36;1.54;1.63;1.66;1.55;1.83;1.88;1.91;1.9;1.65;"
	electroneg += "1.81;2.01;2.18;2.55;2.96;0;0.82;0.95;1.22;1.33;1.6;2.16;1.9;2.2;2.28;"
	electroneg += "2.2;1.93;1.69;1.78;1.96;2.05;2.1;2.66;2.6;0.79;0.89;1.1;1.12;1.13;1.14;"
	electroneg += "1.13;1.17;1.2;1.2;1.1;1.22;1.23;1.24;1.25;1.1;1.27;1.3;1.5;2.36;1.9;"
	electroneg += "2.2;2.2;2.28;2.54;2;2.04;2.33;2.02;2;2.2;0;0.7;0.89;1.1;1.3;"
	return (str2num(StringFromList(Z-1,electroneg)))
End

Function Element_atomVol(Z)
	Variable Z
	String atomVol
	atomVol   = "14.1;31.8;13.1;5;4.6;5.3;17.3;14;17.1;16.9;23.7;14;10;12.1;17;"
	atomVol += "15.5;18.7;24.2;45.3;29.9;15;10.6;8.35;7.23;7.39;7.1;6.7;6.6;7.1;9.2;"
	atomVol += "11.8;13.6;13.1;16.5;23.5;32.2;55.9;33.7;19.8;14.1;10.8;9.4;8.5;8.3;8.3;"
	atomVol += "8.9;10.3;13.1;15.7;16.3;18.4;20.5;25.7;42.9;70;39;22.5;21;20.8;20.6;"
	atomVol += "22.4;19.9;28.9;19.9;19.2;19;18.7;18.4;18.1;24.8;17.8;13.6;10.9;9.53;8.85;"
	atomVol += "8.43;8.54;9.1;10.2;14.8;17.2;18.3;21.3;22.7;;50.5;;45.2;22.5;19.9;"
	return (str2num(StringFromList(Z-1,atomVol)))
End

Function Element_DebyeT(Z)	// Debye Temperature (K)
	// room temperature values of DebyeT from: http://www.knowledgedoor.com/2/elements_handbook/debye_temperature.html
	// Al value from Acta Cryst. (1972). A28, 22-27    [ doi:10.1107/S056773947200004X ]
	Variable Z
	String DebyeT
	DebyeT = "122;NaN;448;1031;1362;1550;NaN;NaN;NaN;74.6;155;330;393;692;576;527;"
	DebyeT += "NaN;92;100;230;476;380;390;424;363;373;386;345;310;237;240;403;275;"
	DebyeT += "152.5;NaN;71.9;59;148;214;250;260;377;422;415;350;275;221;221;129;"
	DebyeT += "254;200;152;109;64;43;116;135;138;138;148;NaN;184;118;155;158;158;"
	DebyeT += "161;163;167;118;116;213;225;312;275;400;228;225;178;92;96;87;116;"
	DebyeT += "NaN;NaN;NaN;39;89;NaN;100;262;300;"
	return (str2num(StringFromList(Z-1,DebyeT)))
End

Function/T Element_names(Z)
	Variable Z
	String names
	names   = "Hydrogen;Helium;Lithium;Beryllium;Boron;Carbon;Nitrogen;Oxygen;"
	names += "Fluorine;Neon;Sodium;Magnesium;Aluminum;Silicon;Phosphorus;Sulfur;"
	names += "Chlorine;Argon;Potassium;Calcium;Scandium;Titanium;Vanadium;"
	names += "Chromium;Manganese;Iron;Cobalt;Nickel;Copper;Zinc;Gallium;"
	names += "Germanium;Arsenic;Selenium;Bromine;Krypton;Rubidium;Strontium;"
	names += "Yttrium;Zirconium;Niobium;Molybdenum;Technetium;Ruthenium;"
	names += "Rhodium;Palladium;Silver;Cadmium;Indium;Tin;Antimony;Tellurium;"
	names += "Iodine;Xenon;Cesium;Barium;Lanthanum;Cerium;Praseodymium;"
	names += "Neodymium;Promethium;Samarium;Europium;Gadolinium;Terbium;"
	names += "Dysprosium;Holmium;Erbium;Thulium;Ytterbium;Lutetium;"
	names += "Hafnium;Tantalum;Tungsten;Rhenium;Osmium;Iridium;Platinum;Gold;"
	names += "Mercury;Thallium;Lead;Bismuth;Polonium;Astatine;Radon;Francium;"
	names += "Radium;Actinium;Thorium;Protactinium;Uranium;Neptunium;Plutonium;"
	names += "Americium;Curium;Berkelium;Californium;Einsteinium;Fermium;"
	names += "Mendelevium;Nobelium;Lawrencium;Rutherfordium;Dubnium;Seaborgium;"
	names += "Bohrium;Hassium;Meitnerium;Darmstadtium;Roentgenium;Copernicium;"
	names += "Ununtrium;Flerovium;Ununpentium;Livermorium;"
	return (StringFromList(Z-1,names))
End

Function BindingEnergy(symb,edgeType)
	String symb				// atomic symbol
	String edgeType				// item to return

	Variable topLevel = ItemsInList(GetRTStackInfo(0))<2
	if ((strlen(symb)<1 || strlen(edgeType)<1) && topLevel)
		String edgeTypes="K;L1;L2;L3;M1;M2;M3;M4;M5;N1;N2;N3;N4;N5;N6;N7;O1;O2;O3;O4;O5;P1;P2;P3"
		Prompt edgeType, "edge type",popup, edgeTypes
		Prompt symb, "atomic symbol", popup, ELEMENT_Symbols
		DoPrompt "select", symb,edgeType
		if (V_flag)
			return NaN
		endif
	endif
	if (strlen(symb)<1 || strlen(edgeType)<1)
		return NaN
	endif

	Variable Z=element2Z(symb)
	Wave/T edgeStrings=root:Packages:Elements:edgeStrings
	Variable eV = NumberByKey(edgeType, edgeStrings[Z])
	eV = numtype(Z) ? NaN : eV
	if (topLevel)
		printf "the %s(%s) edge is at %g eV\r",symb,edgeType,eV
	endif
	return eV
End

Function Make_edgeStrings()
	Make/N=93/O/T root:Packages:Elements:edgeStrings
	Wave/T edgeStrings=root:Packages:Elements:edgeStrings
	edgeStrings[1] = "K:13.6;"
	edgeStrings[2] = "K:24.6;"
	edgeStrings[3] = "K:54.7;"
	edgeStrings[4] = "K:111.5;"
	edgeStrings[5] = "K:188;"
	edgeStrings[6] = "K:284.2;"
	edgeStrings[7] = "K:409.9;L1:37.3;"
	edgeStrings[8] = "K:543.1;L1:41.6;"
	edgeStrings[9] = "K:696.7;"
	edgeStrings[10] = "K:870.2;L1:48.5;L2:21.7;L3:21.6;"
	edgeStrings[11] = "K:1070.8;L1:63.5;L2:30.65;L3:30.81;"
	edgeStrings[12] = "K:1303;L1:88.7;L2:49.78;L3:49.5;"
	edgeStrings[13] = "K:1559.6;L1:117.8;L2:72.95;L3:72.55;"
	edgeStrings[14] = "K:1839;L1:149.7;L2:99.82;L3:99.42;"
	edgeStrings[15] = "K:2145.5;L1:189;L2:136;L3:135;"
	edgeStrings[16] = "K:2472;L1:230.9;L2:163.6;L3:162.5;"
	edgeStrings[17] = "K:2822.4;L1:270;L2:202;L3:200;"
	edgeStrings[18] = "K:3205.9;L1:326.3;L2:250.6;L3:248.4;M1:29.3;M2:15.9;M3:15.7;"
	edgeStrings[19] = "K:3608.4;L1:378.6;L2:297.3;L3:294.6;M1:34.8;M2:18.3;M3:18.3;"
	edgeStrings[20] = "K:4038.5;L1:438.4;L2:349.7;L3:346.2;M1:44.3;M2:25.4;M3:25.4;"
	edgeStrings[21] = "K:4492;L1:498;L2:403.6;L3:398.7;M1:51.1;M2:28.3;M3:28.3;"
	edgeStrings[22] = "K:4966;L1:560.9;L2:460.2;L3:453.8;M1:58.7;M2:32.6;M3:32.6;"
	edgeStrings[23] = "K:5465;L1:626.7;L2:519.8;L3:512.1;M1:66.3;M2:37.2;M3:37.2;"
	edgeStrings[24] = "K:5989;L1:696;L2:583.8;L3:574.1;M1:74.1;M2:42.2;M3:42.2;"
	edgeStrings[25] = "K:6539;L1:769.1;L2:649.9;L3:638.7;M1:82.3;M2:47.2;M3:47.2;"
	edgeStrings[26] = "K:7112;L1:844.6;L2:719.9;L3:706.8;M1:91.3;M2:52.7;M3:52.7;"
	edgeStrings[27] = "K:7709;L1:925.1;L2:793.2;L3:778.1;M1:101;M2:58.9;M3:59.9;"
	edgeStrings[28] = "K:8333;L1:1008.6;L2:870;L3:852.7;M1:110.8;M2:68;M3:66.2;"
	edgeStrings[29] = "K:8979;L1:1096.7;L2:952.3;L3:932.7;M1:122.5;M2:77.3;M3:75.1;"
	edgeStrings[30] = "K:9659;L1:1196.2;L2:1044.9;L3:1021.8;M1:139.8;M2:91.4;M3:88.6;M4:10.2;M5:10.1;"
	edgeStrings[31] = "K:10367;L1:1299;L2:1143.2;L3:1116.4;M1:159.5;M2:103.5;M3:100;M4:18.7;M5:18.7;"
	edgeStrings[32] = "K:11103;L1:1414.6;L2:1248.1;L3:1217;M1:180.1;M2:124.9;M3:120.8;M4:29.8;M5:29.2;"
	edgeStrings[33] = "K:11867;L1:1527;L2:1359.1;L3:1323.6;M1:204.7;M2:146.2;M3:141.2;M4:41.7;M5:41.7;"
	edgeStrings[34] = "K:12658;L1:1652;L2:1474.3;L3:1433.9;M1:229.6;M2:166.5;M3:160.7;M4:55.5;M5:54.6;"
	edgeStrings[35] = "K:13474;L1:1782;L2:1596;L3:1550;M1:257;M2:189;M3:182;M4:70;M5:69;"
	edgeStrings[36] = "K:14326;L1:1921;L2:1730.9;L3:1678.4;M1:292.8;M2:222.2;M3:214.4;M4:95;M5:93.8;N1:27.5;N2:14.1;N3:14.1;"
	edgeStrings[37] = "K:15200;L1:2065;L2:1864;L3:1804;M1:326.7;M2:248.7;M3:239.1;M4:113;M5:112;N1:30.5;N2:16.3;N3:15.3;"
	edgeStrings[38] = "K:16105;L1:2216;L2:2007;L3:1940;M1:358.7;M2:280.3;M3:270;M4:136;M5:134.2;N1:38.9;N2:21.3;N3:20.1;"
	edgeStrings[39] = "K:17038;L1:2373;L2:2156;L3:2080;M1:392;M2:310.6;M3:298.8;M4:157.7;M5:155.8;N1:43.8;N2:24.4;N3:23.1;"
	edgeStrings[40] = "K:17998;L1:2532;L2:2307;L3:2223;M1:430.3;M2:343.5;M3:329.8;M4:181.1;M5:178.8;N1:50.6;N2:28.5;N3:27.1;"
	edgeStrings[41] = "K:18986;L1:2698;L2:2465;L3:2371;M1:466.6;M2:376.1;M3:360.6;M4:205;M5:202.3;N1:56.4;N2:32.6;N3:30.8;"
	edgeStrings[42] = "K:20000;L1:2866;L2:2625;L3:2520;M1:506.3;M2:411.6;M3:394;M4:231.1;M5:227.9;N1:63.2;N2:37.6;N3:35.5;"
	edgeStrings[43] = "K:21044;L1:3043;L2:2793;L3:2677;M1:544;M2:447.6;M3:417.7;M4:257.6;M5:253.9;N1:69.5;N2:42.3;N3:39.9;"
	edgeStrings[44] = "K:22117;L1:3224;L2:2967;L3:2838;M1:586.1;M2:483.5;M3:461.4;M4:284.2;M5:280;N1:75;N2:46.3;N3:43.2;"
	edgeStrings[45] = "K:23220;L1:3412;L2:3146;L3:3004;M1:628.1;M2:521.3;M3:496.5;M4:311.9;M5:307.2;N1:81.4;N2:50.5;N3:47.3;"
	edgeStrings[46] = "K:24350;L1:3604;L2:3330;L3:3173;M1:671.6;M2:559.9;M3:532.3;M4:340.5;M5:335.2;N1:87.1;N2:55.7;N3:50.9;"
	edgeStrings[47] = "K:25514;L1:3806;L2:3524;L3:3351;M1:719;M2:603.8;M3:573;M4:374;M5:368.3;N1:97;N2:63.7;N3:58.3;"
	edgeStrings[48] = "K:26711;L1:4018;L2:3727;L3:3538;M1:772;M2:652.6;M3:618.4;M4:411.9;M5:405.2;N1:109.8;N2:63.9;N3:63.9;N4:11.7;N5:10.7;"
	edgeStrings[49] = "K:27940;L1:4238;L2:3938;L3:3730;M1:827.2;M2:703.2;M3:665.3;M4:451.4;M5:443.9;N1:122.9;N2:73.5;N3:73.5;N4:17.7;N5:16.9;"
	edgeStrings[50] = "K:29200;L1:4465;L2:4156;L3:3929;M1:884.7;M2:756.5;M3:714.6;M4:493.2;M5:484.9;N1:137.1;N2:83.6;N3:83.6;N4:24.9;N5:23.9;"
	edgeStrings[51] = "K:30491;L1:4698;L2:4380;L3:4132;M1:946;M2:812.7;M3:766.4;M4:537.5;M5:528.2;N1:153.2;N2:95.6;N3:95.6;N4:33.3;N5:32.1;"
	edgeStrings[52] = "K:31814;L1:4939;L2:4612;L3:4341;M1:1006;M2:870.8;M3:820;M4:583.4;M5:573;N1:169.4;N2:103.3;N3:103.3;N4:41.9;N5:40.4;"
	edgeStrings[53] = "K:33169;L1:5188;L2:4852;L3:4557;M1:1072;M2:931;M3:875;M4:630.8;M5:619.3;N1:186;N2:123;N3:123;N4:50.6;N5:48.9;"
	edgeStrings[54] = "K:34561;L1:5453;L2:5107;L3:4786;M1:1148.7;M2:1002.1;M3:940.6;M4:689;M5:676.4;N1:213.2;N2:146.7;N3:145.5;N4:69.5;N5:67.5;O1:23.3;O2:13.4;O3:12.1;"
	edgeStrings[55] = "K:35985;L1:5714;L2:5359;L3:5012;M1:1211;M2:1071;M3:1003;M4:740.5;M5:726.6;N1:232.3;N2:172.4;N3:161.3;N4:79.8;N5:77.5;O1:22.7;O2:14.2;O3:12.1;"
	edgeStrings[56] = "K:37441;L1:5989;L2:5624;L3:5247;M1:1293;M2:1137;M3:1063;M4:795.7;M5:780.5;N1:253.5;N2:192;N3:178.6;N4:92.6;N5:89.9;O1:30.3;O2:17;O3:14.8;"
	edgeStrings[57] = "K:38925;L1:6266;L2:5891;L3:5483;M1:1362;M2:1209;M3:1128;M4:853;M5:836;N1:274.7;N2:205.8;N3:196;N4:105.3;N5:102.5;O1:34.3;O2:19.3;O3:16.8;"
	edgeStrings[58] = "K:40443;L1:6549;L2:6164;L3:5723;M1:1436;M2:1274;M3:1187;M4:902.4;M5:883.8;N1:291;N2:223.2;N3:206.5;N4:109;N6:0.1;N7:0.1;O1:37.8;O2:19.8;O3:17;"
	edgeStrings[59] = "K:41991;L1:6835;L2:6440;L3:5964;M1:1511;M2:1337;M3:1242;M4:948.3;M5:928.8;N1:304.5;N2:236.3;N3:217.6;N4:115.1;N5:115.1;N6:2;N7:2;O1:37.4;O2:22.3;O3:22.3;"
	edgeStrings[60] = "K:43569;L1:7126;L2:6722;L3:6208;M1:1575;M2:1403;M3:1297;M4:1003.3;M5:980.4;N1:319.2;N2:243.3;N3:224.6;N4:120.5;N5:120.5;N6:1.5;N7:1.5;O1:37.5;O2:21.1;O3:21.1;"
	edgeStrings[61] = "K:45184;L1:7428;L2:7013;L3:6459;M2:1471;M3:1357;M4:1052;M5:1027;N2:242;N3:242;N4:120;N5:120;"
	edgeStrings[62] = "K:46834;L1:7737;L2:7312;L3:6716;M1:1723;M2:1541;M3:1420;M4:1110.9;M5:1083.4;N1:347.2;N2:265.6;N3:247.4;N4:129;N5:129;N6:5.2;N7:5.2;O1:37.4;O2:21.3;O3:21.3;"
	edgeStrings[63] = "K:48519;L1:8052;L2:7617;L3:6977;M1:1800;M2:1614;M3:1481;M4:1158.6;M5:1127.5;N1:360;N2:284;N3:257;N4:133;N5:127.7;N6:0;N7:0;O1:32;O2:22;O3:22;"
	edgeStrings[64] = "K:50239;L1:8376;L2:7930;L3:7243;M1:1881;M2:1688;M3:1544;M4:1221.9;M5:1189.6;N1:378.6;N2:286;N3:271;N5:142.6;N6:8.6;N7:8.6;O1:36;O2:28;O3:21;"
	edgeStrings[65] = "K:51996;L1:8708;L2:8252;L3:7514;M1:1968;M2:1768;M3:1611;M4:1276.9;M5:1241.1;N1:396;N2:322.4;N3:284.1;N4:150.5;N5:150.5;N6:7.7;N7:2.4;O1:45.6;O2:28.7;O3:22.6;"
	edgeStrings[66] = "K:53789;L1:9046;L2:8581;L3:7790;M1:2047;M2:1842;M3:1676;M4:1333;M5:1292.6;N1:414.2;N2:333.5;N3:293.2;N4:153.6;N5:153.6;N6:8;N7:4.3;O1:49.9;O2:26.3;O3:26.3;"
	edgeStrings[67] = "K:55618;L1:9394;L2:8918;L3:8071;M1:2128;M2:1923;M3:1741;M4:1392;M5:1351;N1:432.4;N2:343.5;N3:308.2;N4:160;N5:160;N6:8.6;N7:5.2;O1:49.3;O2:30.8;O3:24.1;"
	edgeStrings[68] = "K:57486;L1:9751;L2:9264;L3:8358;M1:2207;M2:2006;M3:1812;M4:1453;M5:1409;N1:449.8;N2:366.2;N3:320.2;N4:167.6;N5:167.6;N7:4.7;O1:50.6;O2:31.4;O3:24.7;"
	edgeStrings[69] = "K:59390;L1:10116;L2:9617;L3:8648;M1:2307;M2:2090;M3:1885;M4:1515;M5:1468;N1:470.9;N2:385.9;N3:332.6;N4:175.5;N5:175.5;N7:4.6;O1:54.7;O2:31.8;O3:25;"
	edgeStrings[70] = "K:61332;L1:10486;L2:9978;L3:8944;M1:2398;M2:2173;M3:1950;M4:1576;M5:1528;N1:480.5;N2:388.7;N3:339.7;N4:191.2;N5:182.4;N6:2.5;N7:1.3;O1:52;O2:30.3;O3:24.1;"
	edgeStrings[71] = "K:63314;L1:10870;L2:10349;L3:9244;M1:2491;M2:2264;M3:2024;M4:1639;M5:1589;N1:506.8;N2:412.4;N3:359.2;N4:206.1;N5:196.3;N6:8.9;N7:7.5;O1:57.3;O2:33.6;O3:26.7;"
	edgeStrings[72] = "K:65351;L1:11271;L2:10739;L3:9561;M1:2601;M2:2365;M3:2108;M4:1716;M5:1662;N1:538;N2:438.2;N3:380.7;N4:220;N5:211.5;N6:15.9;N7:14.2;O1:64.2;O2:38;O3:29.9;"
	edgeStrings[73] = "K:67416;L1:11682;L2:11136;L3:9881;M1:2708;M2:2469;M3:2194;M4:1793;M5:1735;N1:563.4;N2:463.4;N3:400.9;N4:237.9;N5:226.4;N6:23.5;N7:21.6;O1:69.7;O2:42.2;O3:32.7;"
	edgeStrings[74] = "K:69525;L1:12100;L2:11544;L3:10207;M1:2820;M2:2575;M3:2281;M4:1872;M5:1809;N1:594.1;N2:490.4;N3:423.6;N4:255.9;N5:243.5;N6:33.6;N7:31.4;O1:75.6;O2:45.3;O3:36.8;"
	edgeStrings[75] = "K:71676;L1:12527;L2:11959;L3:10535;M1:2932;M2:2682;M3:2367;M4:1949;M5:1883;N1:625.4;N2:518.7;N3:446.8;N4:273.9;N5:260.5;N6:42.9;N7:40.5;O1:83;O2:45.6;O3:34.6;"
	edgeStrings[76] = "K:73871;L1:12968;L2:12385;L3:10871;M1:3049;M2:2792;M3:2457;M4:2031;M5:1960;N1:658.2;N2:549.1;N3:470.7;N4:293.1;N5:278.5;N6:53.4;N7:50.7;O1:84;O2:58;O3:44.5;"
	edgeStrings[77] = "K:76111;L1:13419;L2:12824;L3:11215;M1:3174;M2:2909;M3:2551;M4:2116;M5:2040;N1:691.1;N2:577.8;N3:495.8;N4:311.9;N5:296.3;N6:63.8;N7:60.8;O1:95.2;O2:63;O3:48;"
	edgeStrings[78] = "K:78395;L1:13880;L2:13273;L3:11564;M1:3296;M2:3027;M3:2645;M4:2202;M5:2122;N1:725.4;N2:609.1;N3:519.4;N4:331.6;N5:314.6;N6:74.5;N7:71.2;O1:101.7;O2:65.3;O3:51.7;"
	edgeStrings[79] = "K:80725;L1:14353;L2:13734;L3:11919;M1:3425;M2:3148;M3:2743;M4:2291;M5:2206;N1:762.1;N2:642.7;N3:546.3;N4:353.2;N5:335.1;N6:87.6;N7:84;O1:107.2;O2:74.2;O3:57.2;"
	edgeStrings[80] = "K:83102;L1:14839;L2:14209;L3:12284;M1:3562;M2:3279;M3:2847;M4:2385;M5:2295;N1:802.2;N2:680.2;N3:576.6;N4:378.2;N5:358.8;N6:104;N7:99.9;O1:127;O2:83.1;O3:64.5;O4:9.6;O5:7.8;"
	edgeStrings[81] = "K:85530;L1:15347;L2:14698;L3:12658;M1:3704;M2:3416;M3:2957;M4:2485;M5:2389;N1:846.2;N2:720.5;N3:609.5;N4:405.7;N5:385;N6:122.2;N7:117.8;O1:136;O2:94.6;O3:73.5;O4:14.7;O5:12.5;"
	edgeStrings[82] = "K:88005;L1:15861;L2:15200;L3:13035;M1:3851;M2:3554;M3:3066;M4:2586;M5:2484;N1:891.8;N2:761.9;N3:643.5;N4:434.3;N5:412.2;N6:141.7;N7:136.9;O1:147;O2:106.4;O3:83.3;O4:20.7;O5:18.1;"
	edgeStrings[83] = "K:90526;L1:16388;L2:15711;L3:13419;M1:3999;M2:3696;M3:3177;M4:2688;M5:2580;N1:939;N2:805.2;N3:678.8;N4:464;N5:440.1;N6:162.3;N7:157;O1:159.3;O2:119;O3:92.6;O4:26.9;O5:23.8;"
	edgeStrings[84] = "K:93105;L1:16939;L2:16244;L3:13814;M1:4149;M2:3854;M3:3302;M4:2798;M5:2683;N1:995;N2:851;N3:705;N4:500;N5:473;N6:184;N7:184;O1:177;O2:132;O3:104;O4:31;O5:31;"
	edgeStrings[85] = "K:95730;L1:17493;L2:16785;L3:14214;M1:4317;M2:4008;M3:3426;M4:2909;M5:2787;N1:1042;N2:886;N3:740;N4:533;N5:507;N6:210;N7:210;O1:195;O2:148;O3:115;O4:40;O5:40;"
	edgeStrings[86] = "K:98404;L1:18049;L2:17337;L3:14619;M1:4482;M2:4159;M3:3538;M4:3022;M5:2892;N1:1097;N2:929;N3:768;N4:567;N5:541;N6:238;N7:238;O1:214;O2:164;O3:127;O4:48;O5:48;P1:26;"
	edgeStrings[87] = "K:101137;L1:18639;L2:17907;L3:15031;M1:4652;M2:4327;M3:3663;M4:3136;M5:3000;N1:1153;N2:980;N3:810;N4:603;N5:577;N6:268;N7:268;O1:234;O2:182;O3:140;O4:58;O5:58;P1:34;P2:15;P3:15;"
	edgeStrings[88] = "K:103922;L1:19237;L2:18484;L3:15444;M1:4822;M2:4490;M3:3792;M4:3248;M5:3105;N1:1208;N2:1058;N3:879;N4:636;N5:603;N6:299;N7:299;O1:254;O2:200;O3:153;O4:68;O5:68;P1:44;P2:19;P3:19;"
	edgeStrings[89] = "K:106755;L1:19840;L2:19083;L3:15871;M1:5002;M2:4656;M3:3909;M4:3370;M5:3219;N1:1269;N2:1080;N3:890;N4:675;N5:639;N6:319;N7:319;O1:272;O2:215;O3:167;O4:80;O5:80;"
	edgeStrings[90] = "K:109651;L1:20472;L2:19693;L3:16300;M1:5182;M2:4830;M3:4046;M4:3491;M5:3332;N1:1330;N2:1168;N3:966.4;N4:712.1;N5:675.2;N6:342.4;N7:333.1;O1:290;O2:229;O3:182;O4:92.5;O5:85.4;P1:41.4;P2:24.5;P3:16.6;"
	edgeStrings[91] = "K:112601;L1:21105;L2:20314;L3:16733;M1:5367;M2:5001;M3:4174;M4:3611;M5:3442;N1:1387;N2:1224;N3:1007;N4:743;N5:708;N6:371;N7:360;O1:310;O2:232;O3:232;O4:94;O5:94;"
	edgeStrings[92] = "K:115606;L1:21757;L2:20948;L3:17166;M1:5548;M2:5182;M3:4303;M4:3728;M5:3552;N1:1439;N2:1271;N3:1043;N4:778.3;N5:736.2;N6:388.2;N7:377.4;O1:321;O2:257;O3:192;O4:102.8;O5:94.2;P1:43.9;P2:26.8;P3:16.8;"
End



//  ============================================================================  //
//  =========================== Start of emission line functions ==========================  //

Function EmissionEnergies(symb,edgeType)		// display or return a single emission line energy
	String symb				// atomic symbol
	String edgeType			// item to return

	Wave/T FullEmissionLineInfo = root:Packages:Elements:FullEmissionLineInfo
	if (!WaveExists(FullEmissionLineInfo))
		ElementDataInitPackage()
	endif

	String edgeTypes="All;Ka1;Ka2;Ka1,2;Kb1;Kb2;Kb3;L1;La1;La2;La1,2;Lb1;Lb2;Lg1;Lb2,15;Ma1;"
	Variable ask=0, topLevel=strlen(GetRTStackInfo(2))<1
	if (topLevel)
		ask = WhichListItem(symb, ELEMENT_Symbols )<0
		ask = WhichListItem(edgeType, edgeTypes )<0
	endif
	if (ask)
		Prompt edgeType, "edge type",popup, edgeTypes
		Prompt symb, "atomic symbol", popup, ELEMENT_Symbols
		DoPrompt "select", symb,edgeType
		if (V_flag)
			return NaN
		endif
	endif
	if (WhichListItem(symb, ELEMENT_Symbols )<0 || WhichListItem(edgeType, edgeTypes )<0)
		return NaN
	endif

	Variable Z=element2Z(symb)
	STRUCT EmissionLineStruct em
	StructGet/S em, FullEmissionLineInfo[Z]

	Variable eV = NaN
	if (stringmatch(edgeType,"All"))
		printEmissionLineStruct(em)
		if (em.N >0)
			eV = em.line[em.N - 1].eV	// return the highest energy
		endif
	else
		Variable i, rel=NaN
		for (i=0;i<em.N;i+=1)
			if (stringmatch(edgeType, em.line[i].name))
				eV = em.line[i].eV
				rel = em.line[i].rel
				break
			endif
		endfor
		if ((ask || topLevel) && numtype(eV)==0)
			printf "  %s(%s) is at %g eV  (relative strength = %g)\r",symb,edgeType,eV,rel
		endif
	endif
	return eV
End


// returns the edge name for a given emission line
Function/T emissionEdgeName(emissionLine)
	String emissionLine
	//	emission					edge
	//	Kxxx						K		Ka1;Ka2;Ka1,2;Kb1;Kb2;Kb3
	//	Lb4, Lb3, Lg2,Lg3		L1
	//	Lb1, Lg1 				L2
	//	L1, La1, La2, Lb2		L3
	//	Max						M5
	//	Mbx						M4
	String Lkeys="Lb4:L1;Lb3:L1;Lg2:L1;Lg3:L1;Lb1:L2;Lg1:L2;L1:L3;La1:L3;La2:L3;Lb2:L3"

	emissionLine = emissionLine[0,2]			// only need first 3 characters
	if (StringMatch(emissionLine,"K*"))
		return "K"
	elseif (StringMatch(emissionLine,"Ma*"))
		return "M5"
	elseif (StringMatch(emissionLine,"Mb*"))
		return "M4"
	endif
	return StringByKey(emissionLine,Lkeys)
End

Static Function/T emissionLineList(Z)			// get list of lines to fit
	Variable Z
	// retruns a ";" separated list, each element is of form "eV:strength:name"
	// where the strength is relative to 100, and name is the name of the emission line

	Variable minSep = NumVarOrDefault("root:Packages:Elements:defaultPeakFWHM",0)*MIN_LINE_SEPARATION_FRACTION
	STRUCT EmissionLineStruct em
	Wave/T FullEmissionLineInfo = root:Packages:Elements:FullEmissionLineInfo
	if (!WaveExists(FullEmissionLineInfo))
		return ""
	endif
	StructGet/S em, FullEmissionLineInfo[Z]
	Variable N=em.N
	if (N<1)
		return ""
	endif

	Make/N=(N)/FREE eV0, rel0, rels=NaN,eVs=NaN
	Make/N=(N)/T/FREE name0
	eV0 = em.line[p].eV
	rel0 = em.line[p].rel
	name0 = em.line[p].name
	Sort eV0, eV0,rel0,name0
	String name,names=""
	Variable i, eV, eVi,deV,relsum=0, eVsum=0, m=0, Nsum=0
	for (i=N-1; i>=0; i-=1)		// combine nearby lines into one line
		eVi = eV0[i]
		if (Nsum==0)					// a new point
			eV = eVi
			eVsum = eVi
			relsum = rel0[i]
			name = name0[i]
			Nsum = 1
		elseif ((eV-eVi)<minSep)	// too close, accumulate
			name += "+"+name0[i]
			eVsum += eVi
			relsum += rel0[i]
			Nsum += 1
		else								// a big gap, save last point & reset next one			
			eVs[m] = eVsum / Nsum	// average energy
			rels[m] = relsum
			names += name+";"
			m += 1

			eV = eVi						// reset for next point
			eVsum = eVi
			relsum = rel0[i]
			name = name0[i]
			Nsum = 1
		endif
	endfor
	eVs[m] = eVsum / Nsum			// save the last point
	rels[m] = relsum
	names += name+";"
	N = m+1
	Redimension/N=(N) eVs,rels
	Reverse eVs,rels
	names = reverseList(names)
	String out=""						// form the string
	for (i=0;i<N;i+=1)
		out += num2str(eVs[i])+":"+num2str(rels[i])+":"+StringFromList(i,names)+";"
	endfor
	return out
End


Static Function MakeFullEmissionInfo()
	Make/T/N=93/O root:Packages:Elements:FullEmissionLineInfo=""
	Wave/T FullEmissionLineInfo = root:Packages:Elements:FullEmissionLineInfo

	String eStr = ""						// This contains the input data
	eStr += "H:;"
	eStr += "He:;"
	eStr += "Li:54.3,Ka12,150:;"
	eStr += "Be:108.5,Ka12,150:;"
	eStr += "B:183.3,Ka12,151:;"
	eStr += "C:277,Ka12,147:;"
	eStr += "N:392.4,Ka12,150:;"
	eStr += "O:524.9,Ka12,151:;"
	eStr += "F:676.8,Ka12,148:;"
	eStr += "Ne:848.6,Ka12,150:;"
	eStr += "Na:1041,Ka12,150:;"
	eStr += "Mg:1253.6,Ka12,150:;"
	eStr += "Al:1486.3,Ka2,50:1486.7,Ka1,100:1557.4,Kb1,1:;"
	eStr += "Si:1739.4,Ka2,50:1740,Ka1,100:1835.9,Kb1,2:;"
	eStr += "P:2012.7,Ka2,50:2013.7,Ka1,100:2139.1,Kb1,3:;"
	eStr += "S:2306.6,Ka2,50:2307.8,Ka1,100:2464,Kb1,5:;"
	eStr += "Cl:2620.8,Ka2,50:2622.4,Ka1,100:2815.6,Kb1,6:;"
	eStr += "Ar:2955.6,Ka2,50:2957.7,Ka1,100:3190.5,Kb13,10:;"
	eStr += "K:3311.1,Ka2,50:3313.8,Ka1,100:3589.6,Kb13,11:;"
	eStr += "Ca:3688.1,Ka2,50:3691.7,Ka1,100:4012.7,Kb13,13:;"
	eStr += "Sc:348.3,L1,21:395.4,La12,111:399.6,Lb1,77:4086.1,Ka2,50:4090.6,Ka1,100:4460.5,Kb13,15:;"
	eStr += "Ti:395.3,L1,46:452.2,La12,111:458.4,Lb1,79:4504.9,Ka2,50:4510.8,Ka1,100:4931.8,Kb13,15:;"
	eStr += "V:446.5,L1,28:511.3,La12,111:519.2,Lb1,80:4944.6,Ka2,50:4952.2,Ka1,100:5427.3,Kb13,15:;"
	eStr += "Cr:500.3,L1,17:572.8,La12,111:582.8,Lb1,79:5405.5,Ka2,50:5414.7,Ka1,100:5946.7,Kb13,15:;"
	eStr += "Mn:556.3,L1,15:637.4,La12,111:648.8,Lb1,77:5887.6,Ka2,50:5898.8,Ka1,100:6490.4,Kb13,17:;"
	eStr += "Fe:615.2,L1,10:705,La12,111:718.5,Lb1,66:6390.8,Ka2,50:6403.8,Ka1,100:7058,Kb13,17:;"
	eStr += "Co:677.8,L1,10:776.2,La12,111:791.4,Lb1,76:6915.3,Ka2,51:6930.3,Ka1,100:7649.4,Kb13,17:;"
	eStr += "Ni:742.7,L1,9:851.5,La12,111:868.8,Lb1,68:7460.9,Ka2,51:7478.2,Ka1,100:8264.7,Kb13,17:;"
	eStr += "Cu:811.1,L1,8:929.7,La12,111:949.8,Lb1,65:8027.8,Ka2,51:8047.8,Ka1,100:8905.3,Kb13,17:;"
	eStr += "Zn:884,L1,7:1011.7,La12,111:1034.7,Lb1,65:8615.8,Ka2,51:8638.9,Ka1,100:9572,Kb13,17:;"
	eStr += "Ga:957.2,L1,7:1097.9,La12,111:1124.8,Lb1,66:9224.8,Ka2,51:9251.7,Ka1,100:10260.3,Kb3,5:10264.2,Kb1,66:;"
	eStr += "Ge:1036.2,L1,6:1188,La12,111:1218.5,Lb1,60:9855.3,Ka2,51:9886.4,Ka1,100:10978,Kb3,6:10982.1,Kb1,60:;"
	eStr += "As:1120,L1,6:1282,La12,111:1317,Lb1,60:10508,Ka2,51:10543.7,Ka1,100:11720.3,Kb3,6:11726.2,Kb1,13:11864,Kb2,1:;"
	eStr += "Se:1204.4,L1,6:1379.1,La12,111:1419.2,Lb1,59:11181.4,Ka2,52:11222.4,Ka1,100:12489.6,Kb3,6:12495.9,Kb1,13:12652,Kb2,1:;"
	eStr += "Br:1293.5,L1,5:1480.4,La12,111:1525.9,Lb1,59:11877.6,Ka2,52:11924.2,Ka1,100:13284.5,Kb3,7:13291.4,Kb1,14:13469.5,Kb2,1:;"
	eStr += "Kr:1386,L1,5:1586,La12,111:1636.6,Lb1,57:2598,Ka2,52:12649,Ka1,100:14104,Kb3,7:14112,Kb1,14:14315,Kb2,2:;"
	eStr += "Rb:1482.4,L1,5:1692.6,La2,11:1694.1,La1,100:1752.2,Lb1,58:13335.8,Ka2,52:13395.3,Ka1,100:14951.7,Kb3,7:14961.3,Kb1,14:15185,Kb2,2:;"
	eStr += "Sr:1582.2,L1,5:1804.7,La2,11:1806.6,La1,100:1871.7,Lb1,58:14097.9,Ka2,52:14165,Ka1,100:15824.9,Kb3,7:15835.7,Kb1,14:16084.6,Kb2,3:;"
	eStr += "Y:1685.4,L1,5:1920.5,La2,11:1922.6,La1,100:1995.8,Lb1,57:14882.9,Ka2,52:14958.4,Ka1,100:16725.8,Kb3,8:16737.8,Kb1,15:17015.4,Kb2,3:;"
	eStr += "Zr:1792,L1,5:2039.9,La2,11:2042.4,La1,100:2124.4,Lb1,54:2219.4,Lb215,1:2302.7,Lg1,2:15690.9,Ka2,52:15775.1,Ka1,100:17654,Kb3,8:17667.8,Kb1,15:17970,Kb2,3:;"
	eStr += "Nb:1902.2,L1,5:2163,La2,11:2165.9,La1,100:2257.4,Lb1,52:2367,Lb215,3:2461.8,Lg1,2:16521,Ka2,52:16615.1,Ka1,100:18606.3,Kb3,8:18622.5,Kb1,15:18953,Kb2,3:;"
	eStr += "Mo:2015.7,L1,5:2289.8,La2,11:2293.2,La1,100:2394.8,Lb1,53:2518.3,Lb215,5:2623.5,Lg1,3:17374.3,Ka2,52:17479.3,Ka1,100:19590.3,Kb3,8:19608.3,Kb1,15:19965.2,Kb2,3:;"
	eStr += "Tc:2122,L1,5:2420,La2,11:2424,La1,100:2538,Lb1,54:2674,Lb215,7:2792,Lg1,3:18250.8,Ka2,53:18367.1,Ka1,100:20599,Kb3,8:20619,Kb1,16:21005,Kb2,4:;"
	eStr += "Ru:2252.8,L1,4:2554.3,La2,11:2558.6,La1,100:2683.2,Lb1,54:2836,Lb215,10:2964.5,Lg1,4:19150.4,Ka2,53:19279.2,Ka1,100:21634.6,Kb3,8:21656.8,Kb1,16:22074,Kb2,4:;"
	eStr += "Rh:2376.5,L1,4:2692,La2,11:2696.7,La1,100:2834.4,Lb1,52:3001.3,Lb215,10:3143.8,Lg1,5:20073.7,Ka2,53:20216.1,Ka1,100:22698.9,Kb3,8:22723.6,Kb1,16:23172.8,Kb2,4:;"
	eStr += "Pd:2503.4,L1,4:2833.3,La2,11:2838.6,La1,100:2990.2,Lb1,53:3171.8,Lb215,12:3328.7,Lg1,6:21020.1,Ka2,53:21177.1,Ka1,100:23791.1,Kb3,8:23818.7,Kb1,16:24299.1,Kb2,4:;"
	eStr += "Ag:2633.7,L1,4:2978.2,La2,11:2984.3,La1,100:3150.9,Lb1,56:3347.8,Lb215,13:3519.6,Lg1,6:21990.3,Ka2,53:22162.9,Ka1,100:24911.5,Kb3,9:24942.4,Kb1,16:25456.4,Kb2,4:;"
	eStr += "Cd:2767.4,L1,4:3126.9,La2,11:3133.7,La1,100:3316.6,Lb1,58:3528.1,Lb215,15:3716.9,Lg1,6:22984.1,Ka2,53:23173.6,Ka1,100:26061.2,Kb3,9:26095.5,Kb1,17:26643.8,Kb2,4:;"
	eStr += "In:2904.4,L1,4:3279.3,La2,11:3286.9,La1,100:3487.2,Lb1,58:3713.8,Lb215,15:3920.8,Lg1,6:24002,Ka2,53:24209.7,Ka1,100:27237.7,Kb3,9:27275.9,Kb1,17:27860.8,Kb2,5:;"
	eStr += "Sn:3045,L1,4:3435.4,La2,11:3444,La1,100:3662.8,Lb1,60:3904.9,Lb215,16:4131.1,Lg1,7:25044,Ka2,53:25271.3,Ka1,100:28444,Kb3,9:28486,Kb1,17:29109.3,Kb2,5:;"
	eStr += "Sb:3188.6,L1,4:3595.3,La2,11:3604.7,La1,100:3843.6,Lb1,61:4100.8,Lb215,17:4347.8,Lg1,8:26110.8,Ka2,54:26359.1,Ka1,100:29679.2,Kb3,9:29725.6,Kb1,18:30389.5,Kb2,5:;"
	eStr += "Te:3335.6,L1,4:3758.8,La2,11:3769.3,La1,100:4029.6,Lb1,61:4301.7,Lb215,18:4570.9,Lg1,8:27201.7,Ka2,54:27472.3,Ka1,100:30944.3,Kb3,9:30995.7,Kb1,18:31700.4,Kb2,5:;"
	eStr += "I:3485,L1,4:3926,La2,11:3937.6,La1,100:4220.7,Lb1,61:4507.5,Lb215,19:4800.9,Lg1,8:28317.2,Ka2,54:28612,Ka1,100:32239.4,Kb3,9:32294.7,Kb1,18:33042,Kb2,5:;"
	eStr += "Xe:3636,L1,4:4093,La2,11:4109.9,La1,100:4414,Lb1,60:4714,Lb215,20:5034,Lg1,8:29458,Ka2,54:29779,Ka1,100:33562,Kb3,9:33624,Kb1,18:34415,Kb2,5:;"
	eStr += "Cs:3795,L1,4:4272.2,La2,11:4286.5,La1,100:4619.8,Lb1,61:4935.9,Lb215,20:5280.4,Lg1,8:30625.1,Ka2,54:30972.8,Ka1,100:34919.4,Kb3,9:34986.9,Kb1,18:35822,Kb2,6:;"
	eStr += "Ba:3954.1,L1,4:4450.9,La2,11:4466.3,La1,100:4827.5,Lb1,60:5156.5,Lb215,20:5531.1,Lg1,9:31817.1,Ka2,54:32193.6,Ka1,100:36304,Kb3,10:36378.2,Kb1,18:37257,Kb2,6:;"
	eStr += "La:833,Ma1,100:4124,L1,4:4634.2,La2,11:4651,La1,100:5042.1,Lb1,60:5383.5,Lb215,21:5788.5,Lg1,9:33034.1,Ka2,54:33441.8,Ka1,100:37720.2,Kb3,10:37801,Kb1,19:38729.9,Kb2,6:;"
	eStr += "Ce:883,Ma1,100:4287.5,L1,4:4823,La2,11:4840.2,La1,100:5262.2,Lb1,61:5613.4,Lb215,21:6052,Lg1,9:34278.9,Ka2,55:34719.7,Ka1,100:39170.1,Kb3,10:39257.3,Kb1,19:40233,Kb2,6:;"
	eStr += "Pr:929.2,Ma1,100:4453.2,L1,4:5013.5,La2,11:5033.7,La1,100:5488.9,Lb1,61:5850,Lb215,21:6322.1,Lg1,9:35550.2,Ka2,55:36026.3,Ka1,100:40652.9,Kb3,10:40748.2,Kb1,19:41773,Kb2,6:;"
	eStr += "Nd:978,Ma1,100:4633,L1,4:5207.7,La2,11:5230.4,La1,100:5721.6,Lb1,60:6089.4,Lb215,21:6602.1,Lg1,10:36847.4,Ka2,55:37361,Ka1,100:42166.5,Kb3,10:42271.3,Kb1,19:43335,Kb2,6:;"
	eStr += "Pm:4809,L1,4:5408,La2,11:5432,La1,100:5961,Lb1,61:6339,Lb2,21:6892,Lg1,10:38171.2,Ka2,55:38724.7,Ka1,100:43713,Kb3,10:43826,Kb1,19:44942,Kb2,6:;"
	eStr += "Sm:1081,Ma1,100:4994.5,L1,4:5609,La2,11:5636.1,La1,100:6205.1,Lb1,61:6587,Lb215,21:7178,Lg1,10:39522.4,Ka2,55:40118.1,Ka1,100:45289,Kb3,10:45413,Kb1,19:46578,Kb2,6:;"
	eStr += "Eu:1131,Ma1,100:5177.2,L1,4:5816.6,La2,11:5845.7,La1,100:6456.4,Lb1,62:6843.2,Lb215,21:7480.3,Lg1,10:40901.9,Ka2,56:41542.2,Ka1,100:46903.6,Kb3,10:47037.9,Kb1,19:48256,Kb2,6:;"
	eStr += "Gd:1185,Ma1,100:5362.1,L1,4:6025,La2,11:6057.2,La1,100:6713.2,Lb1,62:7102.8,Lb215,21:7785.8,Lg1,11:42308.9,Ka2,56:42996.2,Ka1,100:48555,Kb3,10:48697,Kb1,20:49959,Kb2,7:;"
	eStr += "Tb:1240,Ma1,100:5546.7,L1,4:6238,La2,11:6272.8,La1,100:6978,Lb1,61:7366.7,Lb215,21:8102,Lg1,11:43744.1,Ka2,56:44481.6,Ka1,100:50229,Kb3,10:50382,Kb1,20:51698,Kb2,7:;"
	eStr += "Dy:1293,Ma1,100:5743.1,L1,4:6457.7,La2,11:6495.2,La1,100:7247.7,Lb1,62:7635.7,Lb2,20:8418.8,Lg1,11:45207.8,Ka2,56:45998.4,Ka1,100:51957,Kb3,10:52119,Kb1,20:53476,Kb2,7:;"
	eStr += "Ho:1348,Ma1,100:5943.4,L1,4:6679.5,La2,11:6719.8,La1,100:7525.3,Lb1,64:7911,Lb215,20:8747,Lg1,11:46699.7,Ka2,56:47546.7,Ka1,100:53711,Kb3,11:53877,Kb1,20:55293,Kb2,7:;"
	eStr += "Er:1406,Ma1,100:6152,L1,4:6905,La2,11:6948.7,La1,100:7810.9,Lb1,64:8189,Lb215,20:9089,Lg1,11:48221.1,Ka2,56:49127.7,Ka1,100:55494,Kb3,11:55681,Kb1,21:57210,Kb2,7:;"
	eStr += "Tm:1462,Ma1,100:6341.9,L1,4:7133.1,La2,11:7179.9,La1,100:8101,Lb1,64:8468,Lb215,20:9426,Lg1,12:49772.6,Ka2,57:50741.6,Ka1,100:57304,Kb3,11:57517,Kb1,21:59090,Kb2,7:;"
	eStr += "Yb:1521.4,Ma1,100:6545.5,L1,4:7367.3,La2,11:7415.6,La1,100:8401.8,Lb1,65:8758.8,Lb215,20:9780.1,Lg1,12:51354,Ka2,57:52388.9,Ka1,100:59140,Kb3,11:59370,Kb1,21:60980,Kb2,7:;"
	eStr += "Lu:1581.3,Ma1,100:6752.8,L1,4:7604.9,La2,11:7655.5,La1,100:8709,Lb1,66:9048.9,Lb2,19:10143.4,Lg1,12:52965,Ka2,57:54069.8,Ka1,100:61050,Kb3,11:61283,Kb1,21:62970,Kb2,7:;"
	eStr += "Hf:1644.6,Ma1,100:6959.6,L1,5:7844.6,La2,11:7899,La1,100:9022.7,Lb1,67:9347.3,Lb2,20:10515.8,Lg1,12:54611.4,Ka2,57:55790.2,Ka1,100:62980,Kb3,11:63234,Kb1,22:64980,Kb2,7:;"
	eStr += "Ta:1709.6,Ma1,100:7173.1,L1,5:8087.9,La2,11:8146.1,La1,100:9343.1,Lb1,67:9651.8,Lb2,20:10895.2,Lg1,12:56277,Ka2,57:57532,Ka1,100:64948.8,Kb3,11:65223,Kb1,22:66990,Kb2,7:;"
	eStr += "W:1775.4,Ma1,100:7387.8,L1,5:8335.2,La2,11:8397.6,La1,100:9672.4,Lb1,67:9961.5,Lb2,21:11285.9,Lg1,13:57981.7,Ka2,58:59318.2,Ka1,100:66951.4,Kb3,11:67244.3,Kb1,22:69067,Kb2,8:;"
	eStr += "Re:1842.5,ma1,100:7603.6,L1,5:8586.2,La2,11:8652.5,La1,100:10010,Lb1,66:10275.2,Lb2,22:11685.4,Lg1,13:59717.9,Ka2,58:61140.3,Ka1,100:68994,Kb3,12:69310,Kb1,22:71232,Kb2,8:;"
	eStr += "Os:1910.2,Ma1,100:7822.2,L1,5:8841,La2,11:8911.7,La1,100:10355.3,Lb1,67:10598.5,Lb2,22:12095.3,Lg1,13:61486.7,Ka2,58:63000.5,Ka1,100:71077,Kb3,12:71413,Kb1,23:73363,Kb2,8:;"
	eStr += "Ir:1979.9,Ma1,100:8045.8,L1,5:9099.5,La2,11:9175.1,La1,100:10708.3,Lb1,66:10920.3,Lb2,22:12512.6,Lg1,13:63286.7,Ka2,58:64895.6,Ka1,100:73202.7,Kb3,12:73560.8,Kb1,23:75575,Kb2,8:;"
	eStr += "Pt:2050.5,Ma1,100:8268,L1,5:9361.8,La2,11:9442.3,La1,100:11070.7,Lb1,67:11250.5,Lb2,23:12942,Lg1,13:65112,Ka2,58:66832,Ka1,100:75368,Kb3,12:75748,Kb1,23:77850,Kb2,8:;"
	eStr += "Au:2122.9,Ma1,100:8493.9,L1,5:9628,La2,11:9713.3,La1,100:11442.3,Lb1,67:11584.7,Lb2,23:13381.7,Lg1,13:66989.5,Ka2,59:68803.7,Ka1,100:77580,Kb3,12:77984,Kb1,23:80150,Kb2,8:;"
	eStr += "Hg:2195.3,Ma1,100:8721,L1,5:9897.6,La2,11:9988.8,La1,100:11822.6,Lb1,67:11924.1,Lb2,24:13830.1,Lg1,14:68895,Ka2,59:70819,Ka1,100:79822,Kb3,12:80253,Kb1,23:82515,Kb2,8:;"
	eStr += "Tl:2270.6,Ma1,100:8953.2,L1,6:10172.8,La2,11:10268.5,La1,100:12213.3,Lb1,67:12271.5,Lb2,25:14291.5,Lg1,14:70831.9,Ka2,60:72871.5,Ka1,100:82118,Kb3,12:82576,Kb1,23:84910,Kb2,8:;"
	eStr += "Pb:2345.5,Ma1,100:9184.5,L1,6:10449.5,La2,11:10551.5,La1,100:12613.7,Lb1,66:12622.6,Lb2,25:14764.4,Lg1,14:72804.2,Ka2,60:74969.4,Ka1,100:84450,Kb3,12:84936,Kb1,23:87320,Kb2,8:;"
	eStr += "Bi:2422.6,Ma1,100:9420.4,L1,6:10730.9,La2,11:10838.8,La1,100:12979.9,Lb2,25:13023.5,Lb1,67:15247.7,Lg1,14:74814.8,Ka2,60:77107.9,Ka1,100:86834,Kb3,12:87343,Kb1,23:89830,Kb2,9:;"
	eStr += "Po:;"
	eStr += "At:;"
	eStr += "Rn:;"
	eStr += "Fr:;"
	eStr += "Ra:;"
	eStr += "Ac:;"
	eStr += "Th:2996.1,Ma1,100:11118.6,L1,6:12809.6,La2,11:12968.7,La1,100:15623.7,Lb2,26:16202.2,Lb1,69:18982.5,Lg1,16:89953,Ka2,62:93350,Ka1,100:104831,Kb3,12:105609,Kb1,24:108640,Kb2,9:;"
	eStr += "Pa:;"
	eStr += "U:3170.8,Ma1,100:11618.3,L1,7:13438.8,La2,11:13614.7,La1,100:16428.3,Lb2,26:17220,Lb1,61:20167.1,Lg1,15:94665,Ka2,62:98439,Ka1,100:110406,Kb3,13:111300,Kb1,24:114530,Kb2,9:;"
 
	String substitutions="Kb13:Kb1,3;La12:La1,2"		// used to fix line names

	STRUCT EmissionLineStruct em
	String symb,list="", item, name, newName, strStruct
	Variable Z,i
	for (Z=1;Z<=92;Z+=1)
		symb = StringFromList(Z-1,ELEMENT_Symbols)
		list = StringByKey(symb,eStr)
		em.symb = symb
		em.Z = Z
		em.N = ItemsInList(list,":")
		for (i=0;i<ItemsInList(list,":");i+=1)			// loop over each line for this element
			item = StringFromList(i,list,":")
			name = StringFromList(1,item,",")
			newName = StringByKey(name, substitutions)
			newName = SelectString(strlen(newName),name,newName)	
			em.line[i].eV = str2num(StringFromList(0,item,","))
			em.line[i].name = newName
			em.line[i].rel = str2num(StringFromList(2,item,","))
		endfor
		StructPut/S em, strStruct
		FullEmissionLineInfo[Z] = strStruct
	endfor
End


Structure EmissionLineStruct
	char	symb[2]
	int16	Z
	int16	N			// number of emission lines stored
	STRUCT oneEmissionLineStruct line[20]
EndStructure
//
Static Structure oneEmissionLineStruct
	char	name[6]	// name of line, e.g. "Ka1"
	double	eV			// energy of line
	double	rel			// relative strength
EndStructure
//
Function printEmissionLineStruct(em)
	STRUCT EmissionLineStruct &em
	printf "for '%s' (%g):\r",em.symb,em.Z
	if (em.N<1)
		print "\tempty\r"
	else
		Variable i
		for (i=0;i<em.N;i+=1)
			printf "\t%s \t%g eV \t%g\r",em.line[i].name,em.line[i].eV,em.line[i].rel
		endfor
	endif
End



//  ============================================================================  //
//  ========================== Start of Find Element From Line ===========================  //

Function/T FindBestElementFromEmissionLine(eV,acceptableLines,[dE,excitation])	// searches for the right element
	Variable eV
	String acceptableLines	// acceptable lines, check all if "",  examples: "Ka*;Kb*", or "K*", "", "Lb1",  or "*L*"
	Variable dE					// resolution, find strongest line within eV±dE, for dE=0, just use closest line
	Variable excitation		// excitation energy (eV), defaults to Inf
	dE = ParamIsDefault(dE) || numtype(dE) || dE<0? 0 : dE
	excitation = ParamIsDefault(excitation) || !(excitation>0) ? Inf : excitation

	String list, line, item, lineMin=""
	Variable dEmin=Inf,eMin,Zmin, strengthBest=-Inf
	Variable iline, Z, eni, strength, deltaE, better
	for (Z=1;Z<=92;Z+=1)											// for each element
		list = emissionLineList(Z)
		for (iline=0;iline<ItemsInList(list);iline+=1)	// for each line in an element
			item = StringFromList(iline,list)					// a single line for element Z
			eni = str2num(StringFromList(0,item,":"))		// energy of this line
			line = StringFromList(2,item,":")					// name of this line
			strength = str2num(StringFromList(1,item,":"))
			if (!acceptableLine(acceptableLines,line))		// if line is not of desired type, continue
				continue
			endif

			deltaE = abs(eni-eV)
			better = (deltaE<dE && strength>strengthBest)	// eni is close enough, look for best line, dE>0
			better = better || (dE==0 && deltaE<dEmin)		// find closest line, no width specified
			if (numtype(excitation)==0)							// check that we can excite this line
				better = BindingEnergy(num2str(Z),emissionEdgeName(line)) >= excitation ? better : 0
			endif
			if (better)
				dEmin = deltaE
				eMin = eni
				lineMin = line
				Zmin = Z
				strengthBest = strength
			endif
		endfor
	endfor

	String str
	sprintf str, "%s;%d;%s;%g:%g", StringFromList(Zmin-1,ELEMENT_Symbols),Zmin,lineMin,eMin,strengthBest
	return str
End
//
Static Function acceptableLine(acceptableLines,line)	// returns 1 if line is in acceptableLines
	String acceptableLines, line
	if (strlen(acceptableLines)<1)					// if no lines specified, accept all
		return 1
	endif
	Variable i
	for (i=0;i<ItemsInList(acceptableLines);i+=1)
		if (StringMatch(line,StringFromList(i,acceptableLines)))
			return 1
		endif
	endfor
	return 0
End

//  =========================== End of Find Element From Line ===========================  //
//  ============================================================================  //






//		10543.7	33	As	Ka1		100
//		10508.0	33	As	Ka2		51
//		11720.3	33	As	Kb3	6
//		11726.2	33	As	Kb1	13
//		11864		33	As	Kb2	1
//		1282.0		33	As	La1,2	111
//		1317.0		33	As	Lb1		60
//		1120		33	As	L1		6
//
//	elementLinesZ[0] = Kalpha1(Z)
//	elementLinesZ[1] = Kalpha2(Z)
//	elementLinesZ[2] = Kbeta3(Z)
//	elementLinesZ[3] = Kbeta1(Z)
//	elementLinesZ[4] = Kbeta2(Z)
//	elementLinesZ[5] = Lalpha1(Z)
//	elementLinesZ[6] = Lalpha2(Z)
//	elementLinesZ[7] = Lbeta1(Z)
//	elementLinesZ[8] = Lbeta2(Z)	





//	DEPRECATED     DEPRECATED     DEPRECATED     DEPRECATED     DEPRECATED     DEPRECATED

//	DEPRECATED     DEPRECATED     DEPRECATED     DEPRECATED     DEPRECATED     DEPRECATED

//	DEPRECATED     DEPRECATED     DEPRECATED     DEPRECATED     DEPRECATED     DEPRECATED

//	DEPRECATED     DEPRECATED     DEPRECATED     DEPRECATED     DEPRECATED     DEPRECATED



// The following set of functions is not as good as the previous ones using FullEmissionLineInfo[] (Deprecated)
//
//Static Function MakeAllElementLines()
//	String/G root:Packages:Elements:symbShow="Fe"
//	Make/N=(93,9)/O root:Packages:Elements:allEmissionLInes
//	Wave allEmissionLInes=root:Packages:Elements:allEmissionLInes
//	Variable Z
//	for (Z=1;Z<=92;Z+=1)
//		Wave ww = MakeElementLines(Z)
//		allEmissionLInes[Z][] = ww[q]
//		KIllWaves/Z ww
//	endfor
//End
//
Static Function/WAVE MakeElementLines(Z)
	Variable Z
	if (!(Z>=1 && Z<=92))
		return $""
	endif
	Make/N=9/FREE elementLinesZ=NaN
	elementLinesZ[0] = Kalpha1(Z)
	elementLinesZ[1] = Kalpha2(Z)
	elementLinesZ[2] = Kbeta3(Z)
	elementLinesZ[3] = Kbeta1(Z)
	elementLinesZ[4] = Kbeta2(Z)
	elementLinesZ[5] = Lalpha1(Z)
	elementLinesZ[6] = Lalpha2(Z)
	elementLinesZ[7] = Lbeta1(Z)
	elementLinesZ[8] = Lbeta2(Z)	
	return elementLinesZ
End
//
Static Function Kalpha1(Z)			// K1 - L2
	Variable Z
	WAVE/T edgeString = root:Packages:Elements:edgeStrings
	String list = edgeString[Z]
	Variable eV = NumberByKey("K",list) - NumberByKey("L2",list)
	return eV
End
//
Static Function Kalpha2(Z)			// K1 - L3
	Variable Z
	WAVE/T edgeString = root:Packages:Elements:edgeStrings
	String list = edgeString[Z]
	Variable eV = NumberByKey("K",list) - NumberByKey("L3",list)
	return eV
End
//
Static Function Kbeta3(Z)			// K1 - M2
	Variable Z
	WAVE/T edgeString = root:Packages:Elements:edgeStrings
	String list = edgeString[Z]
	Variable eV = NumberByKey("K",list) - NumberByKey("M2",list)
	return eV
End
//
Static Function Kbeta1(Z)			// K1 - M3
	Variable Z
	WAVE/T edgeString = root:Packages:Elements:edgeStrings
	String list = edgeString[Z]
	Variable eV = NumberByKey("K",list) - NumberByKey("M3",list)
	return eV
End
//
Static Function Kbeta2(Z)			// K1 - (N2+N3)
	Variable Z
	WAVE/T edgeString = root:Packages:Elements:edgeStrings
	String list = edgeString[Z]
	Variable eV = NumberByKey("K",list) - NumberByKey("N2",list)
	eV += NumberByKey("K",list) - NumberByKey("N3",list)
	eV /= 2
	return eV
End
//
Static Function Lalpha1(Z)			// L3 - M5
	Variable Z
	WAVE/T edgeString = root:Packages:Elements:edgeStrings
	String list = edgeString[Z]
	Variable eV = NumberByKey("L3",list) - NumberByKey("M5",list)
	return eV
End
//
Static Function Lalpha2(Z)			// L3 - M4
	Variable Z
	WAVE/T edgeString = root:Packages:Elements:edgeStrings
	String list = edgeString[Z]
	Variable eV = NumberByKey("L3",list) - NumberByKey("M4",list)
	return eV
End
//
Static Function Lbeta1(Z)			// L2 - M4
	Variable Z
	WAVE/T edgeString = root:Packages:Elements:edgeStrings
	String list = edgeString[Z]
	Variable eV = NumberByKey("L2",list) - NumberByKey("M4",list)
	return eV
End
//
Static Function Lbeta2(Z)			// L3 - N5
	Variable Z
	WAVE/T edgeString = root:Packages:Elements:edgeStrings
	String list = edgeString[Z]
	Variable eV = NumberByKey("L3",list) - NumberByKey("N5",list)
	return eV
End

//  =========================== End of emission line functions ===========================  //
//  ============================================================================  //


//  ============================================================================  //
//  ============================= Start of isotope functions =============================  //

Function/WAVE Make_IsotopesList()
	// Use:
	//	print StringByKey("4", IsotopeLists[2])
	// to return the  list pair   "amu,naturalFraction"

	Variable f
	String fileName = ParseFilePath(1,FunctionPath("Make_IsotopesList"),":",1,0)+"data:isotopes.xml"
	Open/R/M="isotope info file"/T=".xml"/Z f as fileName
	if (V_flag)
		DoAlert 0, "ERROR -- could not find 'LaueGo:General:data:isotopes.xml'"
		return $""		
	endif
	FStatus f
	String buf = PadString("",V_logEOF,0x20)
	FBinRead f, buf
	Close f

	String isotopes = XMLtagContents("isotopes",buf)
	if (strlen(isotopes)<100)
		return $""
	endif
	Variable Zmax = ItemsInList(ELEMENT_Symbols)
	Make/N=(Zmax+1)/T/O root:Packages:Elements:IsotopeLists/WAVE=IsotopeLists = ""
	Note/K IsotopeLists, "waveClass=IsotopesList;"
	IsotopeLists[0] = "MassNumber:AtomicMass,NaturalFraction;"

	String oneZ, item, line
	Variable iZ, m, massNumber,atomicMass,Fraction
	for (iZ=1;iZ<=Zmax;iZ+=1)
		oneZ = XMLtagContents("Z"+num2istr(iZ),isotopes)
		line = ""
		item = "X"
		for (m=0; ;m+=1)
			item = XMLtagContents("is",oneZ,occurance=m)
			if (strlen(item)<1)
				break
			endif
			item = ReplaceString(" ",item,":",0,1)
			item = ReplaceString(" ",item,",")
			line += item+";"
		endfor
		IsotopeLists[iZ] = line
	endfor
	return IsotopeLists
End

//  ============================= End of isotope functions ==============================  //
//  ============================================================================  //



Static Function/T ZapControlCodes(str)			// remove parts of string with ASCII code < 32
	String str
	Variable i = 0
	do
		if (char2num(str[i,i])<32)
			str[i,i+1] = str[i+1,i+1]
		endif
		i += 1
	while(i<strlen(str))
	return str
End



Function ElementDataInitPackage([fwhm])
	Variable fwhm
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:Elements
//	String/G root:Packages:Elements:symbols = ELEMENT_Symbols	// use of this global is DEPRECATED
//	SVAR symbols=root:Packages:Elements:symbols
//	symbols = ELEMENT_Symbols
//	symbols   = "H;He;Li;Be;B;C;N;O;F;Ne;Na;Mg;Al;Si;P;S;Cl;Ar;"
//	symbols += "K;Ca;Sc;Ti;V;Cr;Mn;Fe;Co;Ni;Cu;Zn;Ga;Ge;As;Se;Br;Kr;"
//	symbols += "Rb;Sr;Y;Zr;Nb;Mo;Tc;Ru;Rh;Pd;Ag;Cd;In;Sn;Sb;Te;I;Xe;"
//	symbols += "Cs;Ba;La;Ce;Pr;Nd;Pm;Sm;Eu;Gd;Tb;Dy;Ho;Er;Tm;Yb;Lu;"
//	symbols += "Hf;Ta;W;Re;Os;Ir;Pt;Au;Hg;Tl;Pb;Bi;Po;At;Rn;Fr;Ra;Ac;Th;Pa;U;"
//	symbols += "Np;Pu;Am;Cm;Bk;Cf;Es;Fm;Md;No;Lr;Rf;Db;Sg;Bh;Hs;Mt;Ds;Rg;Cn;;Fl;Lv;"

	Make_edgeStrings()
	if (Exists("root:Packages:Elements:defaultPeakFWHM")!=2 || !ParamIsDefault(fwhm))
		fwhm = (numtype(fwhm) || fwhm<0 || ParamIsDefault(fwhm)) ? 0 : fwhm
		Variable /G root:Packages:Elements:defaultPeakFWHM=fwhm
	endif
	MakeFullEmissionInfo()
//	MakeAllElementLines()
	Make_IsotopesList()
End
