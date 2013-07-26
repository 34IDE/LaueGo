#pragma rtGlobals=2		// Use modern global access method.
#pragma IgorVersion = 4.0
#pragma version = 1.75
#pragma ModuleName=elements
#include "MaterialsLocate"								// used to find the path to the materials files
Constant MIN_LINE_SEPARATION_FRACTION = 0.15	// you can over ride this in your main procedure window.
StrConstant ELEMENT_Symbols = "H;He;Li;Be;B;C;N;O;F;Ne;Na;Mg;Al;Si;P;S;Cl;Ar;K;Ca;Sc;Ti;V;Cr;Mn;Fe;Co;Ni;Cu;Zn;Ga;Ge;As;Se;Br;Kr;Rb;Sr;Y;Zr;Nb;Mo;Tc;Ru;Rh;Pd;Ag;Cd;In;Sn;Sb;Te;I;Xe;Cs;Ba;La;Ce;Pr;Nd;Pm;Sm;Eu;Gd;Tb;Dy;Ho;Er;Tm;Yb;Lu;Hf;Ta;W;Re;Os;Ir;Pt;Au;Hg;Tl;Pb;Bi;Po;At;Rn;Fr;Ra;Ac;Th;Pa;U;Np;Pu;Am;Cm;Bk;Cf;Es;Fm;Md;No;Lr;Rf;Db;Sg;Bh;Hs;Mt;Ds;Rg;Cn;;Fl;;Lv"
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
//Function Element_amu(Z)
//	Variable Z
//	String amuList
//	amuList   = "1.00797;4.0026;6.941;9.01218;10.81;12.011;14.0067;15.9994;18.9984;20.179;22.9898;24.305;26.9815;28.0855;"
//	amuList += "30.9738;32.06;35.453;39.948;39.0983;40.08;44.9559;47.88;50.9415;51.996;54.938;55.847;58.9332;58.69;"
//	amuList += "63.546;65.38;69.72;72.59;74.9216;78.96;79.904;83.8;85.4678;87.62;88.9059;91.22;92.9064;95.94;"
//	amuList += "98;101.07;102.906;106.42;107.868;112.41;114.82;118.69;121.75;127.6;126.905;131.29;132.905;137.33;"
//	amuList += "138.906;140.12;140.908;144.24;145;150.36;151.96;157.25;158.925;162.5;164.93;167.26;168.934;173.04;"
//	amuList += "174.967;178.49;180.948;183.85;186.207;190.2;192.22;195.08;196.967;200.59;204.383;207.2;208.98;209;"
//	amuList += "210;222;223;226.025;227.028;232.038;231.036;238.029;"
//	amuList += "237;242;243;247;247;249;254;253;256;254;257"
//	return (str2num(StringFromList(Z-1,amuList)))
//End

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
	atomRadius += "1.35;1.36;1.39;1.46;1.6;1.71;1.75;1.7;1.67;1.45;1.34;2.7;2.33;1.88;1.8;"
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
	covRadius += "1.26;1.27;1.3;1.34;1.49;1.48;1.47;1.46;1.53;1.47;;;;;1.65;"
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

	Variable z=element2Z(symb)
	Wave/T edgeStrings=root:Packages:Elements:edgeStrings
	Variable eV = NumberByKey(edgeType, edgeStrings[z])
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


Static Function/T emissionLineList(Z)			// get list of lines to fit
	Variable Z

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
	for (i=N-1; i>=0; i-=1)				// combine nearby lines into one line
		eVi = eV0[i]
		if (Nsum==0)						// a new point
			eV = eVi
			eVsum = eVi
			relsum = rel0[i]
			name = name0[i]
			Nsum = 1
		elseif ((eV-eVi)<minSep)		// too close, accumulate
			name += "+"+name0[i]
			eVsum += eVi
			relsum += rel0[i]
			Nsum += 1
		else								// a big gap, save last point & reset next one			
			eVs[m] = eVsum / Nsum		// average energy
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
	eVs[m] = eVsum / Nsum				// save the last point
	rels[m] = relsum
	names += name+";"
	N = m+1
	Redimension/N=(N) eVs,rels
	Reverse eVs,rels
	names = reverseList(names)
	String out=""							// form the string
	for (i=0;i<N;i+=1)
		out += num2str(eVs[i])+":"+num2str(rels[i])+":"+StringFromList(i,names)+";"
	endfor
	return out
End


Static Function MakeFullEmissionInfo()
	Make/T/N=92/O root:Packages:Elements:FullEmissionLineInfo=""
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

Function/T FindBestElementFromEmissionLine(eV,acceptableLines,[dE])		// searches for the right element
	Variable eV
	String acceptableLines		// acceptable lines, check all if "",  examples: "Ka*;Kb*", or "K*", "", "Lb1",  or "*L*"
	Variable dE					// resolution, find strongest line within eV±dE, for dE=0, just use closest line
	dE = ParamIsDefault(dE) || numtype(dE) || dE<0? 0 : dE

	String list, line, item, lineMin=""
	Variable dEmin=Inf,eMin,Zmin, strengthBest=-Inf
	Variable iline, Z, eni, strength, deltaE, better
	for (Z=1;Z<=92;Z+=1)									// for each element
		list = emissionLineList(Z)
		for (iline=0;iline<ItemsInList(list);iline+=1)		// for each line in an element
			item = StringFromList(iline,list)				// a single line for element Z
			eni = str2num(StringFromList(0,item,":"))	// energy of this line
			line = StringFromList(2,item,":")				// name of this line
			strength = str2num(StringFromList(1,item,":"))
			if (!acceptableLine(acceptableLines,line))		// if line is not of desired type, continue
				continue
			endif

			deltaE = abs(eni-eV)
			better = (deltaE<dE && strength>strengthBest)	// eni is close enough, look for best line, dE>0
			better = better || (dE==0 && deltaE<dEmin)		// find closest line, no width specified
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

Function Make_IsotopesList()
	// Use:
	//	print StringByKey("4", IsotopeLists[2])
	// to return the  list pair   "amu,naturalFraction"
	Make/N=117/O/T root:Packages:Elements:IsotopeLists
	Wave/T IsotopeLists=root:Packages:Elements:IsotopeLists
	IsotopeLists[0]="MassNumber:AtomicMass,NaturalFraction;"
	IsotopeLists[1]="1:1.0078250321,0.999885;2:2.014101778,0.000115;3:3.0160492675,0;4:4.02783,0;5:5.03954,0;6:6.04494,0;"
	IsotopeLists[2]="3:3.0160293097,1.37e-06;4:4.0026032497,0.99999863;5:5.01222,0;6:6.0188881,0;7:7.02803,0;8:8.033922,0;9:9.04382,0;10:10.0524,0;"
	IsotopeLists[3]="4:4.02718,0;5:5.01254,0;6:6.0151223,0.0759;7:7.016004,0.9241;8:8.0224867,0;9:9.0267891,0;10:10.035481,0;11:11.043796,0;12:12.05378,0;"
	IsotopeLists[4]="5:5.04079,0;6:6.019726,0;7:7.0169292,0;8:8.00530509,0;9:9.0121821,1;10:10.0135337,0;11:11.021658,0;12:12.026921,0;13:13.03613,0;14:14.04282,0;"
	IsotopeLists[5]="7:7.02992,0;8:8.0246067,0;9:9.0133288,0;10:10.012937,0.199;11:11.0093055,0.801;12:12.0143521,0;13:13.0177803,0;14:14.025404,0;15:15.031097,0;16:16.0398"
	IsotopeLists[5]+="1,0;17:17.04693,0;18:18.05617,0;19:19.06373,0;"
	IsotopeLists[6]="8:8.037675,0;9:9.0310401,0;10:10.0168531,0;11:11.0114338,0;12:12,0.9893;13:13.0033548378,0.0107;14:14.003241988,0;15:15.0105993,0;16:16.014701,0;17:17."
	IsotopeLists[6]+="022584,0;18:18.02676,0;19:19.03525,0;20:20.04032,0;21:21.04934,0;22:22.05645,0;"
	IsotopeLists[7]="10:10.04262,0;11:11.0268,0;12:12.0186132,0;13:13.00573858,0;14:14.0030740052,0.99632;15:15.0001088984,0.00368;16:16.0061014,0;17:17.00845,0;18:18.01408"
	IsotopeLists[7]+="2,0;19:19.017027,0;20:20.02337,0;21:21.02709,0;22:22.03444,0;23:23.04051,0;24:24.0505,0;"
	IsotopeLists[8]="12:12.034405,0;13:13.02481,0;14:14.00859529,0;15:15.0030654,0;16:15.9949146221,0.99757;17:16.9991315,0.00038;18:17.9991604,0.00205;19:19.003579,0;20:20"
	IsotopeLists[8]+=".0040762,0;21:21.008655,0;22:22.00997,0;23:23.01569,0;24:24.02037,0;25:25.02914,0;26:26.03775,0;"
	IsotopeLists[9]="14:14.03608,0;15:15.01801,0;16:16.011466,0;17:17.00209524,0;18:18.0009377,0;19:18.9984032,1;20:19.99998132,0;21:20.9999489,0;22:22.002999,0;23:23.00357"
	IsotopeLists[9]+=",0;24:24.0081,0;25:25.01209,0;26:26.01963,0;27:27.02689,0;28:28.03567,0;29:29.04326,0;"
	IsotopeLists[10]="16:16.025757,0;17:17.0177,0;18:18.0056971,0;19:19.0018798,0;20:19.9924401759,0.9048;21:20.99384674,0.0027;22:21.99138551,0.0925;23:22.99446734,0;24:23."
	IsotopeLists[10]+="993615,0;25:24.99779,0;26:26.00046,0;27:27.00762,0;28:28.01211,0;29:29.01935,0;30:30.02387,0;31:31.03311,0;32:32.03991,0;"
	IsotopeLists[11]="18:18.02718,0;19:19.013879,0;20:20.007348,0;21:20.9976551,0;22:21.9944368,0;23:22.98976967,1;24:23.99096333,0;25:24.9899544,0;26:25.99259,0;27:26.99401"
	IsotopeLists[11]+=",0;28:27.99889,0;29:29.00281,0;30:30.00923,0;31:31.0136,0;32:32.01965,0;33:33.02739,0;34:34.0349,0;35:35.04418,0;"
	IsotopeLists[12]="20:20.018863,0;21:21.011714,0;22:21.9995741,0;23:22.9941249,0;24:23.9850419,0.7899;25:24.98583702,0.1;26:25.98259304,0.1101;27:26.98434074,0;28:27.9838"
	IsotopeLists[12]+="767,0;29:28.98855,0;30:29.99046,0;31:30.99655,0;32:31.99915,0;33:33.00559,0;34:34.00907,0;35:35.01749,0;36:36.02245,0;37:37.03124,0;"
	IsotopeLists[13]="21:21.02804,0;22:22.01952,0;23:23.007265,0;24:23.999941,0;25:24.9904286,0;26:25.98689166,0;27:26.98153844,1;28:27.98191018,0;29:28.9804448,0;30:29.9829"
	IsotopeLists[13]+="6,0;31:30.983946,0;32:31.98812,0;33:32.99087,0;34:33.99693,0;35:34.99994,0;36:36.00635,0;37:37.01031,0;38:38.0169,0;39:39.0219,0;"
	IsotopeLists[14]="22:22.03453,0;23:23.02552,0;24:24.011546,0;25:25.004107,0;26:25.99233,0;27:26.98670476,0;28:27.9769265327,0.922297;29:28.97649472,0.046832;30:29.973770"
	IsotopeLists[14]+="22,0.030872;31:30.97536327,0;32:31.9741481,0;33:32.978001,0;34:33.978576,0;35:34.98458,0;36:35.98669,0;37:36.993,0;38:37.99598,0;39:39.0023,0;40:40.005"
	IsotopeLists[14]+="8,0;41:41.0127,0;42:42.0161,0;"
	IsotopeLists[15]="24:24.03435,0;25:25.02026,0;26:26.01178,0;27:26.99919,0;28:27.992312,0;29:28.9818014,0;30:29.9783138,0;31:30.97376151,1;32:31.97390716,0;33:32.9717253,"
	IsotopeLists[15]+="0;34:33.973636,0;35:34.9733142,0;36:35.97826,0;37:36.97961,0;38:37.98447,0;39:38.98642,0;40:39.99105,0;41:40.9948,0;42:42.00009,0;43:43.00331,0;44:44.0"
	IsotopeLists[15]+="0988,0;45:45.01514,0;46:46.02383,0;"
	IsotopeLists[16]="26:26.02788,0;27:27.0188,0;28:28.00437,0;29:28.99661,0;30:29.984903,0;31:30.9795544,0;32:31.97207069,0.9493;33:32.9714585,0.0076;34:33.96786683,0.0429;"
	IsotopeLists[16]+="35:34.96903214,0;36:35.96708088,0.0002;37:36.97112572,0;38:37.971163,0;39:38.97514,0;40:39.97547,0;41:40.98003,0;42:41.98149,0;43:42.9866,0;44:43.98832"
	IsotopeLists[16]+=",0;45:44.99482,0;46:45.99957,0;47:47.00762,0;48:48.01299,0;49:49.02201,0;"
	IsotopeLists[17]="28:28.02851,0;29:29.01411,0;30:30.00477,0;31:30.99242,0;32:31.985689,0;33:32.9774518,0;34:33.97376197,0;35:34.96885271,0.7578;36:35.96830695,0;37:36.96"
	IsotopeLists[17]+="59026,0.2422;38:37.96801055,0;39:38.9680077,0;40:39.97042,0;41:40.97065,0;42:41.97317,0;43:42.9742,0;44:43.97854,0;45:44.9797,0;46:45.98412,0;47:46.987"
	IsotopeLists[17]+="95,0;48:47.99485,0;49:48.99989,0;50:50.00773,0;51:51.01353,0;"
	IsotopeLists[18]="30:30.02156,0;31:31.01213,0;32:31.99766,0;33:32.98993,0;34:33.98027,0;35:34.9752567,0;36:35.96754628,0.003365;37:36.9667759,0;38:37.9627322,0.000632;39"
	IsotopeLists[18]+=":38.964313,0;40:39.962383123,0.996003;41:40.9645008,0;42:41.96305,0;43:42.96567,0;44:43.965365,0;45:44.96809,0;46:45.96809,0;47:46.97219,0;48:47.97507,"
	IsotopeLists[18]+="0;49:48.98218,0;50:49.98594,0;51:50.99324,0;52:51.99817,0;53:53.00623,0;"
	IsotopeLists[19]="32:32.02192,0;33:33.00726,0;34:33.99841,0;35:34.988012,0;36:35.981293,0;37:36.97337691,0;38:37.9690801,0;39:38.9637069,0.932581;40:39.96399867,0.000117"
	IsotopeLists[19]+=";41:40.96182597,0.067302;42:41.9624031,0;43:42.960716,0;44:43.96156,0;45:44.9607,0;46:45.961976,0;47:46.961678,0;48:47.965513,0;49:48.96745,0;50:49.972"
	IsotopeLists[19]+="78,0;51:50.97638,0;52:51.98261,0;53:52.98712,0;54:53.99399,0;55:54.99939,0;"
	IsotopeLists[20]="34:34.01412,0;35:35.00477,0;36:35.99309,0;37:36.985872,0;38:37.976319,0;39:38.9707177,0;40:39.9625912,0.96941;41:40.9622783,0;42:41.9586183,0.00647;43:"
	IsotopeLists[20]+="42.9587668,0.00135;44:43.9554811,0.02086;45:44.9561859,0;46:45.9536928,4e-05;47:46.9545465,0;48:47.952534,0.00187;49:48.955673,0;50:49.957518,0;51:50.9"
	IsotopeLists[20]+="6147,0;52:51.9651,0;53:52.97005,0;54:53.97468,0;55:54.98055,0;56:55.98579,0;57:56.99236,0;"
	IsotopeLists[21]="36:36.01492,0;37:37.00305,0;38:37.9947,0;39:38.98479,0;40:39.977964,0;41:40.9692513,0;42:41.9655168,0;43:42.961151,0;44:43.959403,0;45:44.9559102,1;46:"
	IsotopeLists[21]+="45.9551703,0;47:46.952408,0;48:47.952235,0;49:48.950024,0;50:49.952187,0;51:50.953603,0;52:51.95665,0;53:52.95924,0;54:53.963,0;55:54.96743,0;56:55.972"
	IsotopeLists[21]+="66,0;57:56.97704,0;58:57.98307,0;59:58.98804,0;"
	IsotopeLists[22]="38:38.00977,0;39:39.00132,0;40:39.9905,0;41:40.98313,0;42:41.973032,0;43:42.968523,0;44:43.9596902,0;45:44.9581243,0;46:45.9526295,0.0825;47:46.9517638"
	IsotopeLists[22]+=",0.0744;48:47.9479471,0.7372;49:48.9478708,0.0541;50:49.9447921,0.0518;51:50.946616,0;52:51.946898,0;53:52.94973,0;54:53.95087,0;55:54.95512,0;56:55.95"
	IsotopeLists[22]+="799,0;57:56.9629,0;58:57.96611,0;59:58.97196,0;60:59.97564,0;61:60.98202,0;"
	IsotopeLists[23]="40:40.01109,0;41:40.99974,0;42:41.99123,0;43:42.98065,0;44:43.9744,0;45:44.965782,0;46:45.9601995,0;47:46.9549069,0;48:47.9522545,0;49:48.9485169,0;50:"
	IsotopeLists[23]+="49.9471628,0.0025;51:50.9439637,0.9975;52:51.9447797,0;53:52.944343,0;54:53.946444,0;55:54.94724,0;56:55.95036,0;57:56.95236,0;58:57.95665,0;59:58.9593"
	IsotopeLists[23]+=",0;60:59.9645,0;61:60.96741,0;62:61.97314,0;63:62.97675,0;"
	IsotopeLists[24]="42:42.00643,0;43:42.99771,0;44:43.98547,0;45:44.97916,0;46:45.968362,0;47:46.962907,0;48:47.954036,0;49:48.9513411,0;50:49.9460496,0.04345;51:50.944771"
	IsotopeLists[24]+="8,0;52:51.9405119,0.83789;53:52.9406538,0.09501;54:53.9388849,0.02365;55:54.9408442,0;56:55.940645,0;57:56.94375,0;58:57.94425,0;59:58.94863,0;60:59.94"
	IsotopeLists[24]+="973,0;61:60.95409,0;62:61.9558,0;63:62.96186,0;64:63.9642,0;65:64.97037,0;"
	IsotopeLists[25]="44:44.00687,0;45:44.99451,0;46:45.98672,0;47:46.9761,0;48:47.96887,0;49:48.959623,0;50:49.954244,0;51:50.9482155,0;52:51.9455701,0;53:52.9412947,0;54:5"
	IsotopeLists[25]+="3.9403632,0;55:54.9380496,1;56:55.9389094,0;57:56.938287,0;58:57.93999,0;59:58.94045,0;60:59.94319,0;61:60.94446,0;62:61.94797,0;63:62.94981,0;64:63.95"
	IsotopeLists[25]+="373,0;65:64.9561,0;66:65.96082,0;67:66.96382,0;"
	IsotopeLists[26]="45:45.01456,0;46:46.00081,0;47:46.99289,0;48:47.98056,0;49:48.97361,0;50:49.96299,0;51:50.956825,0;52:51.948117,0;53:52.9453123,0;54:53.9396148,0.05845"
	IsotopeLists[26]+=";55:54.938298,0;56:55.9349421,0.91754;57:56.9353987,0.02119;58:57.9332805,0.00282;59:58.9348805,0;60:59.934077,0;61:60.936749,0;62:61.93677,0;63:62.940"
	IsotopeLists[26]+="12,0;64:63.94087,0;65:64.94494,0;66:65.94598,0;67:66.95,0;68:67.95251,0;69:68.9577,0;"
	IsotopeLists[27]="48:48.00176,0;49:48.98972,0;50:49.98154,0;51:50.97072,0;52:51.96359,0;53:52.954225,0;54:53.9484641,0;55:54.9420031,0;56:55.9398439,0;57:56.9362962,0;58"
	IsotopeLists[27]+=":57.9357576,0;59:58.9332002,1;60:59.9338222,0;61:60.9324794,0;62:61.934054,0;63:62.933615,0;64:63.935814,0;65:64.936485,0;66:65.93983,0;67:66.94061,0;6"
	IsotopeLists[27]+="8:67.94436,0;69:68.9452,0;70:69.94981,0;71:70.95173,0;72:71.95641,0;"
	IsotopeLists[28]="50:49.99593,0;51:50.98772,0;52:51.97568,0;53:52.96846,0;54:53.95791,0;55:54.951336,0;56:55.942136,0;57:56.9398,0;58:57.9353479,0.680769;59:58.9343516,0"
	IsotopeLists[28]+=";60:59.9307906,0.262231;61:60.9310604,0.011399;62:61.9283488,0.036345;63:62.9296729,0;64:63.9279696,0.009256;65:64.930088,0;66:65.929115,0;67:66.93157,"
	IsotopeLists[28]+="0;68:67.931845,0;69:68.93518,0;70:69.93614,0;71:70.94,0;72:71.9413,0;73:72.94608,0;74:73.94791,0;75:74.95297,0;76:75.95533,0;77:76.96083,0;78:77.9638,0;"
	IsotopeLists[29]="52:51.99718,0;53:52.98555,0;54:53.97671,0;55:54.96605,0;56:55.95856,0;57:56.949216,0;58:57.9445407,0;59:58.9395041,0;60:59.9373681,0;61:60.9334622,0;62"
	IsotopeLists[29]+=":61.932587,0;63:62.9296011,0.6917;64:63.9297679,0;65:64.9277937,0.3083;66:65.928873,0;67:66.92775,0;68:67.92964,0;69:68.929425,0;70:69.932409,0;71:70.9"
	IsotopeLists[29]+="3262,0;72:71.93552,0;73:72.93649,0;74:73.9402,0;75:74.9417,0;76:75.94599,0;77:76.94795,0;78:77.95281,0;79:78.95528,0;80:79.96189,0;"
	IsotopeLists[30]="54:53.99295,0;55:54.98398,0;56:55.97238,0;57:56.96491,0;58:57.9546,0;59:58.94927,0;60:59.941832,0;61:60.939514,0;62:61.934334,0;63:62.9332156,0;64:63.9"
	IsotopeLists[30]+="291466,0.4863;65:64.9292451,0;66:65.9260368,0.279;67:66.9271309,0.041;68:67.9248476,0.1875;69:68.9265535,0;70:69.925325,0.0062;71:70.927727,0;72:71.926"
	IsotopeLists[30]+="861,0;73:72.92978,0;74:73.92946,0;75:74.93294,0;76:75.93339,0;77:76.93709,0;78:77.93857,0;79:78.94268,0;80:79.94441,0;81:80.95048,0;82:81.95484,0;"
	IsotopeLists[31]="56:55.99491,0;57:56.98293,0;58:57.97425,0;59:58.96337,0;60:59.95706,0;61:60.94917,0;62:61.94418,0;63:62.93914,0;64:63.936838,0;65:64.9327393,0;66:65.93"
	IsotopeLists[31]+="1592,0;67:66.9282049,0;68:67.9279835,0;69:68.925581,0.60108;70:69.926028,0;71:70.924705,0.39892;72:71.9263694,0;73:72.92517,0;74:73.92694,0;75:74.92650"
	IsotopeLists[31]+="1,0;76:75.92893,0;77:76.92928,0;78:77.93166,0;79:78.93292,0;80:79.93659,0;81:80.93775,0;82:81.94316,0;83:82.94687,0;84:83.95234,0;"
	IsotopeLists[32]="58:57.99101,0;59:58.98175,0;60:59.97019,0;61:60.96379,0;62:61.95465,0;63:62.94964,0;64:63.94157,0;65:64.93944,0;66:65.93385,0;67:66.932738,0;68:67.9280"
	IsotopeLists[32]+="97,0;69:68.927972,0;70:69.9242504,0.2084;71:70.924954,0;72:71.9220762,0.2754;73:72.9234594,0.0773;74:73.9211782,0.3628;75:74.9228595,0;76:75.9214027,0."
	IsotopeLists[32]+="0761;77:76.9235485,0;78:77.922853,0;79:78.9254,0;80:79.925445,0;81:80.92882,0;82:81.92955,0;83:82.93451,0;84:83.93731,0;85:84.94269,0;86:85.94627,0;"
	IsotopeLists[33]="60:59.99313,0;61:60.98062,0;62:61.9732,0;63:62.96369,0;64:63.95757,0;65:64.94948,0;66:65.94437,0;67:66.93919,0;68:67.93679,0;69:68.93228,0;70:69.93093,"
	IsotopeLists[33]+="0;71:70.927115,0;72:71.926753,0;73:72.923825,0;74:73.9239291,0;75:74.9215964,1;76:75.9223939,0;77:76.9206477,0;78:77.921829,0;79:78.920948,0;80:79.9225"
	IsotopeLists[33]+="78,0;81:80.922133,0;82:81.9245,0;83:82.92498,0;84:83.92906,0;85:84.93181,0;86:85.93623,0;87:86.93958,0;88:87.94456,0;89:88.94923,0;"
	IsotopeLists[34]="65:64.96466,0;66:65.95521,0;67:66.95009,0;68:67.94187,0;69:68.93956,0;70:69.9335,0;71:70.93227,0;72:71.927112,0;73:72.926767,0;74:73.9224766,0.0089;75:"
	IsotopeLists[34]+="74.9225236,0;76:75.9192141,0.0937;77:76.9199146,0.0763;78:77.9173095,0.2377;79:78.9184998,0;80:79.9165218,0.4961;81:80.9179929,0;82:81.9167,0.0873;83:8"
	IsotopeLists[34]+="2.919119,0;84:83.918465,0;85:84.92224,0;86:85.924271,0;87:86.92852,0;88:87.93142,0;89:88.93602,0;90:89.93942,0;91:90.94537,0;92:91.94933,0;"
	IsotopeLists[35]="67:66.96479,0;68:67.95825,0;69:68.95018,0;70:69.94462,0;71:70.93925,0;72:71.9365,0;73:72.93179,0;74:73.929891,0;75:74.925776,0;76:75.924542,0;77:76.921"
	IsotopeLists[35]+="38,0;78:77.921146,0;79:78.9183376,0.5069;80:79.91853,0;81:80.916291,0.4931;82:81.916805,0;83:82.91518,0;84:83.916504,0;85:84.915608,0;86:85.918797,0;87"
	IsotopeLists[35]+=":86.920711,0;88:87.92407,0;89:88.92639,0;90:89.93063,0;91:90.93397,0;92:91.93926,0;93:92.9431,0;94:93.94868,0;"
	IsotopeLists[36]="69:68.96532,0;70:69.95601,0;71:70.95051,0;72:71.94191,0;73:72.93893,0;74:73.93326,0;75:74.931034,0;76:75.925948,0;77:76.924668,0;78:77.920386,0.0035;79"
	IsotopeLists[36]+=":78.920083,0;80:79.916378,0.0228;81:80.916592,0;82:81.9134846,0.1158;83:82.914136,0.1149;84:83.911507,0.57;85:84.912527,0;86:85.9106103,0.173;87:86.913"
	IsotopeLists[36]+="3543,0;88:87.914447,0;89:88.91763,0;90:89.919524,0;91:90.92344,0;92:91.926153,0;93:92.93127,0;94:93.93436,0;95:94.93984,0;96:95.94307,0;97:96.94856,0;"
	IsotopeLists[37]="71:70.96532,0;72:71.95908,0;73:72.95037,0;74:73.94447,0;75:74.938569,0;76:75.935071,0;77:76.930407,0;78:77.928141,0;79:78.923997,0;80:79.922519,0;81:80"
	IsotopeLists[37]+=".918994,0;82:81.918208,0;83:82.915112,0;84:83.914385,0;85:84.9117893,0.7217;86:85.9111671,0;87:86.9091835,0.2783;88:87.911319,0;89:88.91228,0;90:89.914"
	IsotopeLists[37]+="809,0;91:90.916534,0;92:91.919725,0;93:92.922033,0;94:93.926407,0;95:94.929319,0;96:95.934284,0;97:96.93734,0;98:97.9417,0;99:98.94542,0;100:99.94987,0"
	IsotopeLists[37]+=";101:100.9532,0;102:101.95921,0;"
	IsotopeLists[38]="73:72.96597,0;74:73.95631,0;75:74.94992,0;76:75.94161,0;77:76.93776,0;78:77.932179,0;79:78.929707,0;80:79.924525,0;81:80.923213,0;82:81.918401,0;83:82."
	IsotopeLists[38]+="917555,0;84:83.913425,0.0056;85:84.912933,0;86:85.9092624,0.0986;87:86.9088793,0.07;88:87.9056143,0.8258;89:88.9074529,0;90:89.9077376,0;91:90.91021,0;"
	IsotopeLists[38]+="92:91.91103,0;93:92.914022,0;94:93.91536,0;95:94.919358,0;96:95.92168,0;97:96.926149,0;98:97.928471,0;99:98.93332,0;100:99.93535,0;101:100.94052,0;102:"
	IsotopeLists[38]+="101.94302,0;103:102.94895,0;104:103.95233,0;"
	IsotopeLists[39]="77:76.94962,0;78:77.9435,0;79:78.93735,0;80:79.93434,0;81:80.92913,0;82:81.92679,0;83:82.92235,0;84:83.92039,0;85:84.916427,0;86:85.914888,0;87:86.9108"
	IsotopeLists[39]+="778,0;88:87.9095034,0;89:88.9058479,1;90:89.9071514,0;91:90.907303,0;92:91.908947,0;93:92.909582,0;94:93.911594,0;95:94.912824,0;96:95.915898,0;97:96.9"
	IsotopeLists[39]+="18131,0;98:97.92222,0;99:98.924635,0;100:99.92776,0;101:100.93031,0;102:101.93356,0;103:102.93694,0;104:103.94145,0;105:104.94509,0;106:105.95022,0;"
	IsotopeLists[40]="79:78.94916,0;80:79.94055,0;81:80.93682,0;82:81.93109,0;83:82.92865,0;84:83.92325,0;85:84.92147,0;86:85.91647,0;87:86.914817,0;88:87.910226,0;89:88.908"
	IsotopeLists[40]+="889,0;90:89.9047037,0.5145;91:90.905645,0.1122;92:91.9050401,0.1715;93:92.9064756,0;94:93.9063158,0.1738;95:94.9080427,0;96:95.908276,0.028;97:96.91095"
	IsotopeLists[40]+="1,0;98:97.912746,0;99:98.916511,0;100:99.91776,0;101:100.92114,0;102:101.92298,0;103:102.9266,0;104:103.92878,0;105:104.93305,0;106:105.93591,0;107:106"
	IsotopeLists[40]+=".94086,0;108:107.94428,0;"
	IsotopeLists[41]="81:80.94905,0;82:81.94313,0;83:82.9367,0;84:83.93357,0;85:84.92791,0;86:85.92504,0;87:86.92036,0;88:87.91796,0;89:88.9135,0;90:89.911264,0;91:90.906991"
	IsotopeLists[41]+=",0;92:91.9071932,0;93:92.9063775,1;94:93.9072835,0;95:94.9068352,0;96:95.9081,0;97:96.9080971,0;98:97.910331,0;99:98.911618,0;100:99.914181,0;101:100.9"
	IsotopeLists[41]+="15252,0;102:101.91804,0;103:102.91914,0;104:103.92246,0;105:104.92393,0;106:105.92819,0;107:106.93031,0;108:107.93501,0;109:108.93763,0;110:109.94268,0;"
	IsotopeLists[42]="83:82.94874,0;84:83.94009,0;85:84.93659,0;86:85.9307,0;87:86.92733,0;88:87.921953,0;89:88.919481,0;90:89.913936,0;91:90.911751,0;92:91.90681,0.1484;93:"
	IsotopeLists[42]+="92.906812,0;94:93.9050876,0.0925;95:94.9058415,0.1592;96:95.9046789,0.1668;97:96.906021,0.0955;98:97.9054078,0.2413;99:98.9077116,0;100:99.907477,0.096"
	IsotopeLists[42]+="3;101:100.910347,0;102:101.910297,0;103:102.9132,0;104:103.91376,0;105:104.91697,0;106:105.918134,0;107:106.92169,0;108:107.92358,0;109:108.92781,0;110"
	IsotopeLists[42]+=":109.92973,0;111:110.93451,0;112:111.93684,0;113:112.94203,0;"
	IsotopeLists[43]="85:84.94894,0;86:85.94288,0;87:86.93653,0;88:87.93283,0;89:88.92754,0;90:89.92356,0;91:90.91843,0;92:91.91526,0;93:92.910248,0;94:93.909656,0;95:94.907"
	IsotopeLists[43]+="656,0;96:95.907871,0;97:96.906365,0;98:97.907216,0;99:98.9062546,0;100:99.9076576,0;101:100.907314,0;102:101.909213,0;103:102.909179,0;104:103.91144,0;"
	IsotopeLists[43]+="105:104.91166,0;106:105.914355,0;107:106.91508,0;108:107.91848,0;109:108.91963,0;110:109.92339,0;111:110.92505,0;112:111.92924,0;113:112.93133,0;114:11"
	IsotopeLists[43]+="3.93588,0;115:114.93828,0;"
	IsotopeLists[44]="87:86.94918,0;88:87.94042,0;89:88.93611,0;90:89.92978,0;91:90.92638,0;92:91.92012,0;93:92.91705,0;94:93.91136,0;95:94.910413,0;96:95.907598,0.0554;97:9"
	IsotopeLists[44]+="6.907555,0;98:97.905287,0.0187;99:98.9059393,0.1276;100:99.9042197,0.126;101:100.9055822,0.1706;102:101.9043495,0.3155;103:102.9063237,0;104:103.90543,"
	IsotopeLists[44]+="0.1862;105:104.90775,0;106:105.907327,0;107:106.90991,0;108:107.91019,0;109:108.9132,0;110:109.91397,0;111:110.91756,0;112:111.91855,0;113:112.92254,0;"
	IsotopeLists[44]+="114:113.924,0;115:114.92831,0;116:115.93016,0;117:116.93479,0;118:117.93703,0;"
	IsotopeLists[45]="89:88.94938,0;90:89.94287,0;91:90.93655,0;92:91.93198,0;93:92.92574,0;94:93.9217,0;95:94.9159,0;96:95.914518,0;97:96.91134,0;98:97.910716,0;99:98.90813"
	IsotopeLists[45]+="2,0;100:99.908117,0;101:100.906164,0;102:101.906843,0;103:102.905504,1;104:103.906655,0;105:104.905692,0;106:105.907285,0;107:106.906751,0;108:107.9087"
	IsotopeLists[45]+="3,0;109:108.908736,0;110:109.91095,0;111:110.91166,0;112:111.91461,0;113:112.91542,0;114:113.91885,0;115:114.92012,0;116:115.92371,0;117:116.92535,0;11"
	IsotopeLists[45]+="8:117.92943,0;119:118.93136,0;120:119.93578,0;121:120.93808,0;"
	IsotopeLists[46]="91:90.94948,0;92:91.94042,0;93:92.93591,0;94:93.92877,0;95:94.92469,0;96:95.91822,0;97:96.91648,0;98:97.912721,0;99:98.911768,0;100:99.908505,0;101:100"
	IsotopeLists[46]+=".908289,0;102:101.905608,0.0102;103:102.906087,0;104:103.904035,0.1114;105:104.905084,0.2233;106:105.903483,0.2733;107:106.905128,0;108:107.903894,0.26"
	IsotopeLists[46]+="46;109:108.905954,0;110:109.905152,0.1172;111:110.90764,0;112:111.907313,0;113:112.91015,0;114:113.910365,0;115:114.91368,0;116:115.91416,0;117:116.917"
	IsotopeLists[46]+="84,0;118:117.91898,0;119:118.92268,0;120:119.92403,0;121:120.92818,0;122:121.9298,0;123:122.93426,0;"
	IsotopeLists[47]="94:93.94278,0;95:94.93548,0;96:95.93068,0;97:96.924,0;98:97.92176,0;99:98.9176,0;100:99.91607,0;101:100.9128,0;102:101.912,0;103:102.908972,0;104:103.9"
	IsotopeLists[47]+="08628,0;105:104.906528,0;106:105.906666,0;107:106.905093,0.51839;108:107.905954,0;109:108.904756,0.48161;110:109.90611,0;111:110.905295,0;112:111.90700"
	IsotopeLists[47]+="4,0;113:112.906566,0;114:113.908808,0;115:114.90876,0;116:115.91136,0;117:116.91168,0;118:117.91458,0;119:118.91567,0;120:119.91879,0;121:120.91985,0;1"
	IsotopeLists[47]+="22:121.92332,0;123:122.9249,0;124:123.92853,0;125:124.93054,0;126:125.9345,0;127:126.93688,0;"
	IsotopeLists[48]="96:95.93977,0;97:96.93494,0;98:97.92758,0;99:98.92501,0;100:99.92023,0;101:100.91868,0;102:101.91478,0;103:102.913419,0;104:103.909848,0;105:104.909468"
	IsotopeLists[48]+=",0;106:105.906458,0.0125;107:106.906614,0;108:107.904183,0.0089;109:108.904986,0;110:109.903006,0.1249;111:110.904182,0.128;112:111.9027572,0.2413;113:"
	IsotopeLists[48]+="112.9044009,0.1222;114:113.9033581,0.2873;115:114.905431,0;116:115.904755,0.0749;117:116.907218,0;118:117.906914,0;119:118.90992,0;120:119.909851,0;121"
	IsotopeLists[48]+=":120.91298,0;122:121.9135,0;123:122.917,0;124:123.91765,0;125:124.92125,0;126:125.92235,0;127:126.92643,0;128:127.92776,0;129:128.93226,0;130:129.93398,0;"
	IsotopeLists[49]="98:97.94224,0;99:98.93461,0;100:99.93115,0;101:100.92656,0;102:101.92471,0;103:102.919914,0;104:103.91834,0;105:104.914673,0;106:105.913461,0;107:106.9"
	IsotopeLists[49]+="10292,0;108:107.90972,0;109:108.907154,0;110:109.907169,0;111:110.905111,0;112:111.905533,0;113:112.904061,0.0429;114:113.904917,0;115:114.903878,0.957"
	IsotopeLists[49]+="1;116:115.90526,0;117:116.904516,0;118:117.906355,0;119:118.905846,0;120:119.90796,0;121:120.907849,0;122:121.91028,0;123:122.910439,0;124:123.91318,0;"
	IsotopeLists[49]+="125:124.9136,0;126:125.91646,0;127:126.91734,0;128:127.92017,0;129:128.92166,0;130:129.92485,0;131:130.92677,0;132:131.93292,0;133:132.93834,0;134:133."
	IsotopeLists[49]+="94466,0;"
	IsotopeLists[50]="100:99.93895,0;101:100.93606,0;102:101.93049,0;103:102.92813,0;104:103.92319,0;105:104.92139,0;106:105.91688,0;107:106.91567,0;108:107.91197,0;109:108."
	IsotopeLists[50]+="911287,0;110:109.907853,0;111:110.907735,0;112:111.904821,0.0097;113:112.905173,0;114:113.902782,0.0066;115:114.903346,0.0034;116:115.901744,0.1454;117"
	IsotopeLists[50]+=":116.902954,0.0768;118:117.901606,0.2422;119:118.903309,0.0859;120:119.9021966,0.3258;121:120.9042369,0;122:121.9034401,0.0463;123:122.9057219,0;124:12"
	IsotopeLists[50]+="3.9052746,0.0579;125:124.9077849,0;126:125.907654,0;127:126.910351,0;128:127.910535,0;129:128.91344,0;130:129.91385,0;131:130.91692,0;132:131.917744,0;"
	IsotopeLists[50]+="133:132.92381,0;134:133.92846,0;135:134.93473,0;136:135.93934,0;137:136.94579,0;"
	IsotopeLists[51]="103:102.94012,0;104:103.93629,0;105:104.93153,0;106:105.92876,0;107:106.92415,0;108:107.92216,0;109:108.918136,0;110:109.91676,0;111:110.91321,0;112:11"
	IsotopeLists[51]+="1.912395,0;113:112.909378,0;114:113.9091,0;115:114.906599,0;116:115.906797,0;117:116.90484,0;118:117.905532,0;119:118.903946,0;120:119.905074,0;121:120"
	IsotopeLists[51]+=".903818,0.5721;122:121.9051754,0;123:122.9042157,0.4279;124:123.9059375,0;125:124.905248,0;126:125.90725,0;127:126.906915,0;128:127.909167,0;129:128.90"
	IsotopeLists[51]+="915,0;130:129.911546,0;131:130.91195,0;132:131.914413,0;133:132.91524,0;134:133.92055,0;135:134.92517,0;136:135.93066,0;137:136.93531,0;138:137.94096,0"
	IsotopeLists[51]+=";139:138.94571,0;"
	IsotopeLists[52]="106:105.9377,0;107:106.93504,0;108:107.92949,0;109:108.92746,0;110:109.92241,0;111:110.92112,0;112:111.91706,0;113:112.91593,0;114:113.91206,0;115:114."
	IsotopeLists[52]+="91158,0;116:115.90842,0;117:116.908634,0;118:117.905825,0;119:118.906408,0;120:119.90402,0.0009;121:120.90493,0;122:121.9030471,0.0255;123:122.904273,0"
	IsotopeLists[52]+=".0089;124:123.9028195,0.0474;125:124.9044247,0.0707;126:125.9033055,0.1884;127:126.905217,0;128:127.9044614,0.3174;129:128.906596,0;130:129.9062228,0.3"
	IsotopeLists[52]+="408;131:130.9085219,0;132:131.908524,0;133:132.91094,0;134:133.91154,0;135:134.91645,0;136:135.9201,0;137:136.92532,0;138:137.92922,0;139:138.93473,0;1"
	IsotopeLists[52]+="40:139.9387,0;141:140.94439,0;142:141.9485,0;"
	IsotopeLists[53]="108:107.94329,0;109:108.93819,0;110:109.93521,0;111:110.93028,0;112:111.92797,0;113:112.92364,0;114:113.92185,0;115:114.91792,0;116:115.91674,0;117:116"
	IsotopeLists[53]+=".91365,0;118:117.91338,0;119:118.91018,0;120:119.910048,0;121:120.907366,0;122:121.907592,0;123:122.905598,0;124:123.9062114,0;125:124.9046241,0;126:12"
	IsotopeLists[53]+="5.905619,0;127:126.904468,1;128:127.905805,0;129:128.904987,0;130:129.906674,0;131:130.9061242,0;132:131.907995,0;133:132.907806,0;134:133.909877,0;135"
	IsotopeLists[53]+=":134.91005,0;136:135.91466,0;137:136.917873,0;138:137.92238,0;139:138.92609,0;140:139.93121,0;141:140.93483,0;142:141.94018,0;143:142.94407,0;144:143.9"
	IsotopeLists[53]+="4961,0;"
	IsotopeLists[54]="110:109.94448,0;111:110.94163,0;112:111.93567,0;113:112.93338,0;114:113.92815,0;115:114.92654,0;116:115.92174,0;117:116.92056,0;118:117.91657,0;119:118"
	IsotopeLists[54]+=".91555,0;120:119.91215,0;121:120.911386,0;122:121.90855,0;123:122.908471,0;124:123.9058958,0.0009;125:124.9063982,0;126:125.904269,0.0009;127:126.90518"
	IsotopeLists[54]+=",0;128:127.9035304,0.0192;129:128.9047795,0.2644;130:129.9035079,0.0408;131:130.9050819,0.2118;132:131.9041545,0.2689;133:132.905906,0;134:133.9053945,"
	IsotopeLists[54]+="0.1044;135:134.907207,0;136:135.90722,0.0887;137:136.911563,0;138:137.91399,0;139:138.918787,0;140:139.92164,0;141:140.92665,0;142:141.9297,0;143:142.9"
	IsotopeLists[54]+="3489,0;144:143.93823,0;145:144.94367,0;146:145.9473,0;147:146.95301,0;"
	IsotopeLists[55]="112:111.95033,0;113:112.94454,0;114:113.94142,0;115:114.93594,0;116:115.93291,0;117:116.92864,0;118:117.926555,0;119:118.922371,0;120:119.920678,0;121:"
	IsotopeLists[55]+="120.917184,0;122:121.916122,0;123:122.91299,0;124:123.912246,0;125:124.909725,0;126:125.909448,0;127:126.907418,0;128:127.907748,0;129:128.906063,0;130"
	IsotopeLists[55]+=":129.906706,0;131:130.90546,0;132:131.90643,0;133:132.905447,1;134:133.906713,0;135:134.905972,0;136:135.907306,0;137:136.907084,0;138:137.911011,0;139"
	IsotopeLists[55]+=":138.913358,0;140:139.917277,0;141:140.920044,0;142:141.924292,0;143:142.92733,0;144:143.93203,0;145:144.93539,0;146:145.94016,0;147:146.94386,0;148:14"
	IsotopeLists[55]+="7.9489,0;149:148.95272,0;150:149.95797,0;151:150.962,0;"
	IsotopeLists[56]="114:113.95094,0;115:114.94771,0;116:115.94168,0;117:116.93886,0;118:117.93344,0;119:118.93105,0;120:119.92605,0;121:120.92449,0;122:121.92026,0;123:122"
	IsotopeLists[56]+=".91885,0;124:123.915088,0;125:124.91462,0;126:125.911244,0;127:126.91112,0;128:127.908309,0;129:128.908674,0;130:129.90631,0.00106;131:130.906931,0;132"
	IsotopeLists[56]+=":131.905056,0.00101;133:132.906002,0;134:133.904503,0.02417;135:134.905683,0.06592;136:135.90457,0.07854;137:136.905821,0.11232;138:137.905241,0.71698;"
	IsotopeLists[56]+="139:138.908835,0;140:139.910599,0;141:140.914406,0;142:141.916448,0;143:142.920617,0;144:143.92294,0;145:144.92692,0;146:145.93011,0;147:146.93399,0;14"
	IsotopeLists[56]+="8:147.93768,0;149:148.94246,0;150:149.94562,0;151:150.9507,0;152:151.95416,0;153:152.95961,0;"
	IsotopeLists[57]="117:116.95001,0;118:117.94657,0;119:118.94099,0;120:119.93807,0;121:120.93301,0;122:121.93071,0;123:122.92624,0;124:123.92453,0;125:124.92067,0;126:125"
	IsotopeLists[57]+=".91937,0;127:126.91616,0;128:127.91545,0;129:128.91267,0;130:129.91232,0;131:130.91011,0;132:131.91011,0;133:132.9084,0;134:133.90849,0;135:134.906971,"
	IsotopeLists[57]+="0;136:135.90765,0;137:136.90647,0;138:137.907107,0.0009;139:138.906348,0.9991;140:139.909473,0;141:140.910957,0;142:141.914074,0;143:142.916059,0;144:1"
	IsotopeLists[57]+="43.91959,0;145:144.92164,0;146:145.9257,0;147:146.92782,0;148:147.93219,0;149:148.93437,0;150:149.93857,0;151:150.94156,0;152:151.94611,0;153:152.94945"
	IsotopeLists[57]+=",0;154:153.9544,0;155:154.95813,0;"
	IsotopeLists[58]="119:118.95276,0;120:119.94664,0;121:120.94367,0;122:121.93801,0;123:122.93551,0;124:123.93052,0;125:124.92854,0;126:125.9241,0;127:126.92275,0;128:127."
	IsotopeLists[58]+="91887,0;129:128.91809,0;130:129.91469,0;131:130.91442,0;132:131.91149,0;133:132.91155,0;134:133.90903,0;135:134.909146,0;136:135.90714,0.00185;137:136."
	IsotopeLists[58]+="90778,0;138:137.905986,0.00251;139:138.906647,0;140:139.905434,0.8845;141:140.908271,0;142:141.90924,0.11114;143:142.912381,0;144:143.913643,0;145:144."
	IsotopeLists[58]+="91723,0;146:145.91869,0;147:146.92251,0;148:147.92439,0;149:148.92829,0;150:149.93023,0;151:150.93404,0;152:151.93638,0;153:152.94058,0;154:153.94332,0"
	IsotopeLists[58]+=";155:154.94804,0;156:155.95126,0;157:156.95634,0;"
	IsotopeLists[59]="121:120.95536,0;122:121.95165,0;123:122.94596,0;124:123.94296,0;125:124.93783,0;126:125.93531,0;127:126.93083,0;128:127.9288,0;129:128.92486,0;130:129."
	IsotopeLists[59]+="92338,0;131:130.92006,0;132:131.91912,0;133:132.9162,0;134:133.91567,0;135:134.91314,0;136:135.91265,0;137:136.91068,0;138:137.910749,0;139:138.908932,"
	IsotopeLists[59]+="0;140:139.909071,0;141:140.907648,1;142:141.91004,0;143:142.910812,0;144:143.913301,0;145:144.914507,0;146:145.91759,0;147:146.91898,0;148:147.92218,0;"
	IsotopeLists[59]+="149:148.923791,0;150:149.927,0;151:150.92823,0;152:151.9316,0;153:152.93365,0;154:153.93739,0;155:154.93999,0;156:155.94412,0;157:156.94717,0;158:157.9"
	IsotopeLists[59]+="5178,0;159:158.95523,0;"
	IsotopeLists[60]="126:125.94307,0;127:126.9405,0;128:127.93539,0;129:128.93325,0;130:129.92878,0;131:130.9271,0;132:131.92312,0;133:132.92221,0;134:133.91865,0;135:134.9"
	IsotopeLists[60]+="1824,0;136:135.91502,0;137:136.91464,0;138:137.91193,0;139:138.91192,0;140:139.90931,0;141:140.909605,0;142:141.907719,0.272;143:142.90981,0.122;144:14"
	IsotopeLists[60]+="3.910083,0.238;145:144.912569,0.083;146:145.913112,0.172;147:146.916096,0;148:147.916889,0.057;149:148.920144,0;150:149.920887,0.056;151:150.923825,0;1"
	IsotopeLists[60]+="52:151.92468,0;153:152.927695,0;154:153.92948,0;155:154.93263,0;156:155.9352,0;157:156.93927,0;158:157.94187,0;159:158.94639,0;160:159.94939,0;161:160."
	IsotopeLists[60]+="95433,0;"
	IsotopeLists[61]="128:127.94826,0;129:128.94316,0;130:129.94045,0;131:130.9358,0;132:131.93375,0;133:132.92972,0;134:133.92849,0;135:134.92462,0;136:135.92345,0;137:136."
	IsotopeLists[61]+="92071,0;138:137.91945,0;139:138.91676,0;140:139.9158,0;141:140.913607,0;142:141.91295,0;143:142.910928,0;144:143.912586,0;145:144.912744,0;146:145.9146"
	IsotopeLists[61]+="92,0;147:146.915134,0;148:147.917468,0;149:148.918329,0;150:149.920979,0;151:150.921203,0;152:151.92349,0;153:152.924113,0;154:153.92655,0;155:154.9281"
	IsotopeLists[61]+=",0;156:155.93106,0;157:156.9332,0;158:157.93669,0;159:158.93913,0;160:159.94299,0;161:160.94586,0;162:161.95029,0;163:162.95352,0;"
	IsotopeLists[62]="130:129.94863,0;131:130.94589,0;132:131.94082,0;133:132.93873,0;134:133.93402,0;135:134.93235,0;136:135.9283,0;137:136.92705,0;138:137.92354,0;139:138."
	IsotopeLists[62]+="922302,0;140:139.918991,0;141:140.918469,0;142:141.915193,0;143:142.914624,0;144:143.911995,0.0307;145:144.913406,0;146:145.913037,0;147:146.914893,0.1"
	IsotopeLists[62]+="499;148:147.914818,0.1124;149:148.91718,0.1382;150:149.917271,0.0738;151:150.919928,0;152:151.919728,0.2675;153:152.922094,0;154:153.922205,0.2275;155:"
	IsotopeLists[62]+="154.924636,0;156:155.925526,0;157:156.92835,0;158:157.92999,0;159:158.9332,0;160:159.93514,0;161:160.93883,0;162:161.94122,0;163:162.94536,0;164:163.94"
	IsotopeLists[62]+="828,0;165:164.95298,0;"
	IsotopeLists[63]="132:131.95416,0;133:132.9489,0;134:133.94632,0;135:134.94172,0;136:135.9395,0;137:136.93521,0;138:137.93345,0;139:138.92984,0;140:139.92808,0;141:140.9"
	IsotopeLists[63]+="2489,0;142:141.9234,0;143:142.920287,0;144:143.918774,0;145:144.916261,0;146:145.9172,0;147:146.916741,0;148:147.918154,0;149:148.917926,0;150:149.9196"
	IsotopeLists[63]+="98,0;151:150.919846,0.4781;152:151.92174,0;153:152.921226,0.5219;154:153.922975,0;155:154.922889,0;156:155.924751,0;157:156.925419,0;158:157.92784,0;15"
	IsotopeLists[63]+="9:158.929084,0;160:159.93197,0;161:160.93368,0;162:161.93704,0;163:162.93921,0;164:163.94299,0;165:164.94572,0;166:165.94997,0;167:166.95305,0;"
	IsotopeLists[64]="136:135.94707,0;137:136.94465,0;138:137.93997,0;139:138.93808,0;140:139.93395,0;141:140.93221,0;142:141.92823,0;143:142.92674,0;144:143.92279,0;145:144"
	IsotopeLists[64]+=".92169,0;146:145.918305,0;147:146.919089,0;148:147.91811,0;149:148.919336,0;150:149.918655,0;151:150.920344,0;152:151.919788,0.002;153:152.921746,0;154"
	IsotopeLists[64]+=":153.920862,0.0218;155:154.922619,0.148;156:155.92212,0.2047;157:156.923957,0.1565;158:157.924101,0.2484;159:158.926385,0;160:159.927051,0.2186;161:160"
	IsotopeLists[64]+=".929666,0;162:161.930981,0;163:162.93399,0;164:163.93586,0;165:164.93938,0;166:165.9416,0;167:166.94557,0;168:167.94836,0;169:168.95287,0;"
	IsotopeLists[65]="138:137.95287,0;139:138.94803,0;140:139.94554,0;141:140.94116,0;142:141.93886,0;143:142.93475,0;144:143.93253,0;145:144.92888,0;146:145.92718,0;147:146"
	IsotopeLists[65]+=".924037,0;148:147.9243,0;149:148.923242,0;150:149.923654,0;151:150.923098,0;152:151.92407,0;153:152.923431,0;154:153.92469,0;155:154.9235,0;156:155.924"
	IsotopeLists[65]+="744,0;157:156.924021,0;158:157.92541,0;159:158.925343,1;160:159.927164,0;161:160.927566,0;162:161.92948,0;163:162.930644,0;164:163.93335,0;165:164.9348"
	IsotopeLists[65]+="8,0;166:165.93805,0;167:166.94005,0;168:167.94364,0;169:168.94622,0;170:169.95025,0;171:170.9533,0;"
	IsotopeLists[66]="140:139.95379,0;141:140.95119,0;142:141.94627,0;143:142.94383,0;144:143.93907,0;145:144.93695,0;146:145.93272,0;147:146.93088,0;148:147.92718,0;149:148"
	IsotopeLists[66]+=".927334,0;150:149.92558,0;151:150.92618,0;152:151.924714,0;153:152.925761,0;154:153.924422,0;155:154.925749,0;156:155.924278,0.0006;157:156.925461,0;15"
	IsotopeLists[66]+="8:157.924405,0.001;159:158.925736,0;160:159.925194,0.0234;161:160.92693,0.1891;162:161.926795,0.2551;163:162.928728,0.249;164:163.929171,0.2818;165:164"
	IsotopeLists[66]+=".9317,0;166:165.932803,0;167:166.93565,0;168:167.93723,0;169:168.9403,0;170:169.94267,0;171:170.94648,0;172:171.94911,0;173:172.95344,0;"
	IsotopeLists[67]="142:141.95986,0;143:142.95469,0;144:143.95164,0;145:144.94688,0;146:145.9441,0;147:146.93984,0;148:147.93727,0;149:148.93379,0;150:149.93335,0;151:150."
	IsotopeLists[67]+="931681,0;152:151.93174,0;153:152.930195,0;154:153.930596,0;155:154.929079,0;156:155.92971,0;157:156.92819,0;158:157.92895,0;159:158.927709,0;160:159.92"
	IsotopeLists[67]+="8726,0;161:160.927852,0;162:161.929092,0;163:162.92873,0;164:163.930231,0;165:164.930319,1;166:165.932281,0;167:166.933126,0;168:167.9355,0;169:168.936"
	IsotopeLists[67]+="868,0;170:169.93961,0;171:170.94146,0;172:171.94482,0;173:172.94729,0;174:173.95115,0;175:174.95405,0;"
	IsotopeLists[68]="144:143.96059,0;145:144.95746,0;146:145.95212,0;147:146.94931,0;148:147.94444,0;149:148.94217,0;150:149.93776,0;151:150.93746,0;152:151.93508,0;153:152"
	IsotopeLists[68]+=".935093,0;154:153.932777,0;155:154.9332,0;156:155.93102,0;157:156.93195,0;158:157.92991,0;159:158.930681,0;160:159.92908,0;161:160.930001,0;162:161.928"
	IsotopeLists[68]+="775,0.0014;163:162.930029,0;164:163.929197,0.0161;165:164.930723,0;166:165.93029,0.3361;167:166.932045,0.2293;168:167.932368,0.2678;169:168.934588,0;17"
	IsotopeLists[68]+="0:169.93546,0.1493;171:170.938026,0;172:171.939352,0;173:172.9424,0;174:173.94434,0;175:174.94793,0;176:175.95029,0;177:176.95437,0;"
	IsotopeLists[69]="146:145.9665,0;147:146.96108,0;148:147.95755,0;149:148.95265,0;150:149.94967,0;151:150.94543,0;152:151.9443,0;153:152.942028,0;154:153.94142,0;155:154."
	IsotopeLists[69]+="939192,0;156:155.93901,0;157:156.93676,0;158:157.937,0;159:158.93481,0;160:159.93509,0;161:160.9334,0;162:161.93397,0;163:162.932648,0;164:163.933451,0"
	IsotopeLists[69]+=";165:164.932432,0;166:165.933553,0;167:166.932849,0;168:167.93417,0;169:168.934211,1;170:169.935798,0;171:170.936426,0;172:171.938396,0;173:172.9396,0;"
	IsotopeLists[69]+="174:173.94216,0;175:174.94383,0;176:175.94699,0;177:176.94904,0;178:177.95264,0;179:178.95534,0;"
	IsotopeLists[70]="148:147.96676,0;149:148.96348,0;150:149.95799,0;151:150.95525,0;152:151.95017,0;153:152.94921,0;154:153.94624,0;155:154.94579,0;156:155.94285,0;157:156"
	IsotopeLists[70]+=".94266,0;158:157.939858,0;159:158.94015,0;160:159.93756,0;161:160.93785,0;162:161.93575,0;163:162.93627,0;164:163.93452,0;165:164.935398,0;166:165.9338"
	IsotopeLists[70]+="8,0;167:166.934947,0;168:167.933894,0.0013;169:168.935187,0;170:169.934759,0.0304;171:170.936322,0.1428;172:171.9363777,0.2183;173:172.9382068,0.1613;1"
	IsotopeLists[70]+="74:173.9388581,0.3183;175:174.9412725,0;176:175.942568,0.1276;177:176.945257,0;178:177.946643,0;179:178.95017,0;180:179.95233,0;181:180.95615,0;"
	IsotopeLists[71]="150:149.97267,0;151:150.96715,0;152:151.96361,0;153:152.95869,0;154:153.9571,0;155:154.95423,0;156:155.95291,0;157:156.950102,0;158:157.94917,0;159:158"
	IsotopeLists[71]+=".94662,0;160:159.94602,0;161:160.94354,0;162:161.94322,0;163:162.9412,0;164:163.94122,0;165:164.93961,0;166:165.93976,0;167:166.93831,0;168:167.9387,0;"
	IsotopeLists[71]+="169:168.937649,0;170:169.938472,0;171:170.93791,0;172:171.939082,0;173:172.938927,0;174:173.9403335,0;175:174.9407679,0.9741;176:175.9426824,0.0259;177"
	IsotopeLists[71]+=":176.943755,0;178:177.945951,0;179:178.947324,0;180:179.94988,0;181:180.95197,0;182:181.95521,0;183:182.95757,0;184:183.96117,0;"
	IsotopeLists[72]="154:153.96425,0;155:154.96276,0;156:155.95925,0;157:156.95813,0;158:157.95465,0;159:158.954,0;160:159.95071,0;161:160.95033,0;162:161.947203,0;163:162."
	IsotopeLists[72]+="94706,0;164:163.94442,0;165:164.94454,0;166:165.94225,0;167:166.9426,0;168:167.94063,0;169:168.94116,0;170:169.93965,0;171:170.94049,0;172:171.93946,0;"
	IsotopeLists[72]+="173:172.94065,0;174:173.94004,0.0016;175:174.941503,0;176:175.9414018,0.0526;177:176.94322,0.186;178:177.9436977,0.2728;179:178.9458151,0.1362;180:179."
	IsotopeLists[72]+="9465488,0.3508;181:180.9490991,0;182:181.950553,0;183:182.95353,0;184:183.95545,0;185:184.95878,0;186:185.96092,0;"
	IsotopeLists[73]="156:155.97169,0;157:156.96815,0;158:157.96637,0;159:158.96291,0;160:159.96136,0;161:160.95837,0;162:161.95715,0;163:162.95432,0;164:163.95357,0;165:164"
	IsotopeLists[73]+=".95082,0;166:165.95047,0;167:166.94797,0;168:167.94779,0;169:168.94592,0;170:169.94609,0;171:170.94446,0;172:171.94474,0;173:172.94354,0;174:173.94417,"
	IsotopeLists[73]+="0;175:174.94365,0;176:175.94474,0;177:176.944472,0;178:177.94575,0;179:178.945934,0;180:179.947466,0.00012;181:180.947996,0.99988;182:181.950152,0;183:"
	IsotopeLists[73]+="182.951373,0;184:183.954009,0;185:184.955559,0;186:185.95855,0;187:186.96041,0;188:187.96371,0;"
	IsotopeLists[74]="158:157.97394,0;159:158.97228,0;160:159.96837,0;161:160.96709,0;162:161.96334,0;163:162.96253,0;164:163.95898,0;165:164.95834,0;166:165.95502,0;167:166"
	IsotopeLists[74]+=".95467,0;168:167.95186,0;169:168.95176,0;170:169.94929,0;171:170.94946,0;172:171.94742,0;173:172.94783,0;174:173.94616,0;175:174.94677,0;176:175.94559,"
	IsotopeLists[74]+="0;177:176.94662,0;178:177.94585,0;179:178.947072,0;180:179.946706,0.0012;181:180.948198,0;182:181.948206,0.265;183:182.9502245,0.1431;184:183.9509326,0"
	IsotopeLists[74]+=".3064;185:184.9534206,0;186:185.954362,0.2843;187:186.957158,0;188:187.958487,0;189:188.96191,0;190:189.96318,0;"
	IsotopeLists[75]="160:159.98149,0;161:160.97766,0;162:161.97571,0;163:162.97197,0;164:163.97032,0;165:164.96705,0;166:165.9658,0;167:166.96256,0;168:167.96161,0;169:168."
	IsotopeLists[75]+="95883,0;170:169.95816,0;171:170.95555,0;172:171.95529,0;173:172.95306,0;174:173.95311,0;175:174.95139,0;176:175.95157,0;177:176.95027,0;178:177.95085,0"
	IsotopeLists[75]+=";179:178.94998,0;180:179.95079,0;181:180.950065,0;182:181.95121,0;183:182.950821,0;184:183.952524,0;185:184.9529557,0.374;186:185.954987,0;187:186.9557"
	IsotopeLists[75]+="508,0.626;188:187.9581123,0;189:188.959228,0;190:189.96182,0;191:190.963124,0;192:191.96596,0;"
	IsotopeLists[76]="162:161.98382,0;163:162.98205,0;164:163.97793,0;165:164.97648,0;166:165.97253,0;167:166.97155,0;168:167.96783,0;169:168.96708,0;170:169.96357,0;171:170"
	IsotopeLists[76]+=".96304,0;172:171.96008,0;173:172.95979,0;174:173.95712,0;175:174.95708,0;176:175.95495,0;177:176.95505,0;178:177.95335,0;179:178.95395,0;180:179.95235,"
	IsotopeLists[76]+="0;181:180.95327,0;182:181.952186,0;183:182.95311,0;184:183.952491,0.0002;185:184.954043,0;186:185.953838,0.0159;187:186.9557479,0.0196;188:187.955836,0"
	IsotopeLists[76]+=".1324;189:188.9581449,0.1615;190:189.958445,0.2626;191:190.960928,0;192:191.961479,0.4078;193:192.964148,0;194:193.965179,0;195:194.96812,0;196:195.96962,0;"
	IsotopeLists[77]="165:164.98758,0;166:165.98551,0;167:166.98154,0;168:167.97997,0;169:168.97639,0;170:169.97503,0;171:170.97178,0;172:171.97064,0;173:172.96771,0;174:173"
	IsotopeLists[77]+=".9668,0;175:174.96428,0;176:175.96351,0;177:176.96117,0;178:177.96108,0;179:178.95915,0;180:179.95925,0;181:180.95764,0;182:181.95813,0;183:182.95681,0"
	IsotopeLists[77]+=";184:183.95739,0;185:184.95659,0;186:185.957951,0;187:186.957361,0;188:187.958852,0;189:188.958716,0;190:189.96059,0;191:190.960591,0.373;192:191.96260"
	IsotopeLists[77]+="2,0;193:192.962924,0.627;194:193.965076,0;195:194.965977,0;196:195.96838,0;197:196.969636,0;198:197.97228,0;199:198.97379,0;"
	IsotopeLists[78]="168:167.98804,0;169:168.98642,0;170:169.98233,0;171:170.98125,0;172:171.97738,0;173:172.9765,0;174:173.972811,0;175:174.97228,0;176:175.969,0;177:176.9"
	IsotopeLists[78]+="6845,0;178:177.96571,0;179:178.96548,0;180:179.96322,0;181:180.96318,0;182:181.96127,0;183:182.96173,0;184:183.9599,0;185:184.96075,0;186:185.95943,0;1"
	IsotopeLists[78]+="87:186.96056,0;188:187.959396,0;189:188.960832,0;190:189.95993,0.00014;191:190.961685,0;192:191.961035,0.00782;193:192.962985,0;194:193.962664,0.32967;"
	IsotopeLists[78]+="195:194.964774,0.33832;196:195.964935,0.25242;197:196.967323,0;198:197.967876,0.07163;199:198.970576,0;200:199.971424,0;201:200.9745,0;202:201.97574,0;"
	IsotopeLists[79]="171:170.99177,0;172:171.99011,0;173:172.9864,0;174:173.98492,0;175:174.98155,0;176:175.98027,0;177:176.97722,0;178:177.97598,0;179:178.97341,0;180:179."
	IsotopeLists[79]+="9724,0;181:180.96995,0;182:181.96962,0;183:182.96762,0;184:183.96747,0;185:184.96581,0;186:185.966,0;187:186.96456,0;188:187.96509,0;189:188.96389,0;19"
	IsotopeLists[79]+="0:189.964699,0;191:190.96365,0;192:191.96481,0;193:192.964132,0;194:193.965339,0;195:194.965018,0;196:195.966551,0;197:196.966552,1;198:197.968225,0;19"
	IsotopeLists[79]+="9:198.968748,0;200:199.97072,0;201:200.971641,0;202:201.97379,0;203:202.975137,0;204:203.97771,0;205:204.97961,0;"
	IsotopeLists[80]="175:174.99141,0;176:175.98741,0;177:176.98634,0;178:177.982476,0;179:178.98178,0;180:179.97832,0;181:180.97781,0;182:181.97475,0;183:182.97456,0;184:18"
	IsotopeLists[80]+="3.9719,0;185:184.97198,0;186:185.96946,0;187:186.96979,0;188:187.96756,0;189:188.96813,0;190:189.96628,0;191:190.96706,0;192:191.96557,0;193:192.966644"
	IsotopeLists[80]+=",0;194:193.965382,0;195:194.96664,0;196:195.965815,0.0015;197:196.967195,0;198:197.966752,0.0997;199:198.968262,0.1687;200:199.968309,0.231;201:200.970"
	IsotopeLists[80]+="285,0.1318;202:201.970626,0.2986;203:202.972857,0;204:203.973476,0.0687;205:204.976056,0;206:205.977499,0;207:206.98258,0;208:207.98594,0;"
	IsotopeLists[81]="177:176.99688,0;178:177.99523,0;179:178.99147,0;180:179.99019,0;181:180.9869,0;182:181.98561,0;183:182.9827,0;184:183.98176,0;185:184.9791,0;186:185.97"
	IsotopeLists[81]+="855,0;187:186.97617,0;188:187.97592,0;189:188.97369,0;190:189.97379,0;191:190.97189,0;192:191.97214,0;193:192.97055,0;194:193.97105,0;195:194.96965,0;1"
	IsotopeLists[81]+="96:195.97052,0;197:196.96954,0;198:197.97047,0;199:198.96981,0;200:199.970945,0;201:200.970804,0;202:201.972091,0;203:202.972329,0.29524;204:203.973849"
	IsotopeLists[81]+=",0;205:204.974412,0.70476;206:205.976095,0;207:206.977408,0;208:207.982005,0;209:208.985349,0;210:209.990066,0;"
	IsotopeLists[82]="181:180.99671,0;182:181.992676,0;183:182.99193,0;184:183.9882,0;185:184.98758,0;186:185.9843,0;187:186.98403,0;188:187.98106,0;189:188.98088,0;190:189."
	IsotopeLists[82]+="97818,0;191:190.9782,0;192:191.97576,0;193:192.97608,0;194:193.97397,0;195:194.97447,0;196:195.97271,0;197:196.97338,0;198:197.97198,0;199:198.97291,0;"
	IsotopeLists[82]+="200:199.971816,0;201:200.97285,0;202:201.972144,0;203:202.973375,0;204:203.973029,0.014;205:204.974467,0;206:205.974449,0.241;207:206.975881,0.221;208:"
	IsotopeLists[82]+="207.976636,0.524;209:208.981075,0;210:209.984173,0;211:210.988731,0;212:211.9918875,0;213:212.9965,0;214:213.9997981,0;"
	IsotopeLists[83]="185:184.99771,0;186:185.99648,0;187:186.99346,0;188:187.99217,0;189:188.98951,0;190:189.98852,0;191:190.98605,0;192:191.98537,0;193:192.98306,0;194:193"
	IsotopeLists[83]+=".98275,0;195:194.98075,0;196:195.98061,0;197:196.97893,0;198:197.97902,0;199:198.97758,0;200:199.97814,0;201:200.97697,0;202:201.97767,0;203:202.976868"
	IsotopeLists[83]+=",0;204:203.977805,0;205:204.977375,0;206:205.978483,0;207:206.978455,0;208:207.979727,0;209:208.980383,1;210:209.984105,0;211:210.987258,0;212:211.9912"
	IsotopeLists[83]+="72,0;213:212.994375,0;214:213.998699,0;215:215.00183,0;216:216.0062,0;"
	IsotopeLists[84]="190:189.99511,0;191:190.99465,0;192:191.99152,0;193:192.9911,0;194:193.98828,0;195:194.98805,0;196:195.98551,0;197:196.98557,0;198:197.98334,0;199:198."
	IsotopeLists[84]+="9836,0;200:199.98174,0;201:200.98221,0;202:201.9807,0;203:202.98141,0;204:203.980307,0;205:204.98117,0;206:205.980465,0;207:206.981578,0;208:207.981231"
	IsotopeLists[84]+=",0;209:208.982416,0;210:209.982857,0;211:210.986637,0;212:211.988852,0;213:212.992843,0;214:213.995186,0;215:214.999415,0;216:216.0019052,0;217:217.006"
	IsotopeLists[84]+="25,0;218:218.0089658,0;"
	IsotopeLists[85]="193:193.00019,0;194:193.99897,0;195:194.99655,0;196:195.9957,0;197:196.99329,0;198:197.99275,0;199:198.99063,0;200:199.99029,0;201:200.98849,0;202:201."
	IsotopeLists[85]+="98845,0;203:202.98685,0;204:203.98726,0;205:204.98604,0;206:205.9866,0;207:206.985776,0;208:207.986583,0;209:208.986159,0;210:209.987131,0;211:210.9874"
	IsotopeLists[85]+="81,0;212:211.990735,0;213:212.992921,0;214:213.996356,0;215:214.998641,0;216:216.002409,0;217:217.00471,0;218:218.008681,0;219:219.0113,0;220:220.0153,"
	IsotopeLists[85]+="0;221:221.01814,0;222:222.02233,0;223:223.02534,0;"
	IsotopeLists[86]="196:196.00231,0;197:197.00166,0;198:197.99878,0;199:198.99831,0;200:199.99568,0;201:200.99554,0;202:201.99322,0;203:202.99332,0;204:203.99137,0;205:204"
	IsotopeLists[86]+=".99167,0;206:205.99016,0;207:206.99073,0;208:207.989631,0;209:208.99038,0;210:209.98968,0;211:210.990585,0;212:211.990689,0;213:212.993868,0;214:213.99"
	IsotopeLists[86]+="5346,0;215:214.998729,0;216:216.000258,0;217:217.003915,0;218:218.005586,0;219:219.009475,0;220:220.0113841,0;221:221.01546,0;222:222.0175705,0;223:223"
	IsotopeLists[86]+=".02179,0;224:224.02409,0;225:225.02844,0;226:226.03089,0;227:227.03541,0;228:228.03808,0;"
	IsotopeLists[87]="200:200.0065,0;201:201.00399,0;202:202.00329,0;203:203.00105,0;204:204.00059,0;205:204.99866,0;206:205.99849,0;207:206.99686,0;208:207.99713,0;209:208."
	IsotopeLists[87]+="99592,0;210:209.996398,0;211:210.995529,0;212:211.996195,0;213:212.996175,0;214:213.998955,0;215:215.000326,0;216:216.003188,0;217:217.004616,0;218:218"
	IsotopeLists[87]+=".007563,0;219:219.009241,0;220:220.012313,0;221:221.014246,0;222:222.017544,0;223:223.0197307,0;224:224.02324,0;225:225.025607,0;226:226.02934,0;227:22"
	IsotopeLists[87]+="7.03183,0;228:228.03572,0;229:229.03843,0;230:230.04251,0;231:231.04541,0;232:232.04965,0;"
	IsotopeLists[88]="203:203.00921,0;204:204.00648,0;205:205.00619,0;206:206.00378,0;207:207.00373,0;208:208.00178,0;209:209.00194,0;210:210.00045,0;211:211.00089,0;212:211"
	IsotopeLists[88]+=".999783,0;213:213.00035,0;214:214.000091,0;215:215.002704,0;216:216.003518,0;217:217.006306,0;218:218.007124,0;219:219.010069,0;220:220.011015,0;221:22"
	IsotopeLists[88]+="1.013908,0;222:222.015362,0;223:223.018497,0;224:224.020202,0;225:225.023604,0;226:226.0254026,0;227:227.0291707,0;228:228.0310641,0;229:229.03482,0;23"
	IsotopeLists[88]+="0:230.03708,0;231:231.04122,0;232:232.04369,0;233:233.048,0;234:234.05055,0;"
	IsotopeLists[89]="207:207.01209,0;208:208.01149,0;209:209.00957,0;210:210.00926,0;211:211.00765,0;212:212.00781,0;213:213.00657,0;214:214.00689,0;215:215.00645,0;216:216"
	IsotopeLists[89]+=".008721,0;217:217.009333,0;218:218.01163,0;219:219.0124,0;220:220.01475,0;221:221.01558,0;222:222.017829,0;223:223.019126,0;224:224.021708,0;225:225.02"
	IsotopeLists[89]+="3221,0;226:226.02609,0;227:227.027747,0;228:228.0310148,0;229:229.03293,0;230:230.03603,0;231:231.03855,0;232:232.04202,0;233:233.04455,0;234:234.04842"
	IsotopeLists[89]+=",0;235:235.0511,0;236:236.05518,0;"
	IsotopeLists[90]="210:210.01503,0;211:211.01486,0;212:212.01292,0;213:213.01296,0;214:214.01145,0;215:215.01173,0;216:216.011051,0;217:217.01307,0;218:218.013268,0;219:2"
	IsotopeLists[90]+="19.01552,0;220:220.015733,0;221:221.018171,0;222:222.018454,0;223:223.020795,0;224:224.021459,0;225:225.023941,0;226:226.024891,0;227:227.027699,0;228:"
	IsotopeLists[90]+="228.0287313,0;229:229.031755,0;230:230.0331266,0;231:231.0362971,0;232:232.0380504,1;233:233.0415769,0;234:234.043595,0;235:235.0475,0;236:236.04971,0;"
	IsotopeLists[90]+="237:237.05389,0;238:238.05624,0;"
	IsotopeLists[91]="213:213.02118,0;214:214.02074,0;215:215.0191,0;216:216.01911,0;217:217.01829,0;218:218.02001,0;219:219.01988,0;220:220.02188,0;221:221.02186,0;222:222."
	IsotopeLists[91]+="02373,0;223:223.02396,0;224:224.02561,0;225:225.02612,0;226:226.027933,0;227:227.028793,0;228:228.031037,0;229:229.032089,0;230:230.034533,0;231:231.03"
	IsotopeLists[91]+="58789,1;232:232.038582,0;233:233.0402402,0;234:234.043302,0;235:235.04544,0;236:236.04868,0;237:237.05114,0;238:238.0545,0;239:239.05713,0;240:240.06098,0;"
	IsotopeLists[92]="218:218.02349,0;219:219.02492,0;220:220.02471,0;221:221.02635,0;222:222.02607,0;223:223.02772,0;224:224.02759,0;225:225.02938,0;226:226.02934,0;227:227"
	IsotopeLists[92]+=".03114,0;228:228.031366,0;229:229.033496,0;230:230.033927,0;231:231.036289,0;232:232.0371463,0;233:233.039628,0;234:234.0409456,5.5e-05;235:235.0439231"
	IsotopeLists[92]+=",0.0072;236:236.0455619,0;237:237.048724,0;238:238.0507826,0.992745;239:239.0542878,0;240:240.056586,0;241:241.06033,0;242:242.06293,0;"
	IsotopeLists[93]="225:225.0339,0;226:226.03513,0;227:227.03496,0;228:228.03618,0;229:229.03625,0;230:230.03781,0;231:231.03823,0;232:232.0401,0;233:233.04073,0;234:234.0"
	IsotopeLists[93]+="42889,0;235:235.0440559,0;236:236.04656,0;237:237.0481673,0;238:238.0509405,0;239:239.0529314,0;240:240.056169,0;241:241.05825,0;242:242.06164,0;243:24"
	IsotopeLists[93]+="3.06427,0;244:244.06785,0;"
	IsotopeLists[94]="228:228.03873,0;229:229.04014,0;230:230.039646,0;231:231.04126,0;232:232.041179,0;233:233.04299,0;234:234.043305,0;235:235.045282,0;236:236.0460481,0;2"
	IsotopeLists[94]+="37:237.0484038,0;238:238.0495534,0;239:239.0521565,0;240:240.0538075,0;241:241.0568453,0;242:242.0587368,0;243:243.061997,0;244:244.064198,0;245:245.06"
	IsotopeLists[94]+="7739,0;246:246.070198,0;247:247.07407,0;"
	IsotopeLists[95]="231:231.04556,0;232:232.04659,0;233:233.04647,0;234:234.04779,0;235:235.04803,0;236:236.04957,0;237:237.04997,0;238:238.05198,0;239:239.053018,0;240:24"
	IsotopeLists[95]+="0.055288,0;241:241.0568229,0;242:242.059543,0;243:243.0613727,0;244:244.0642794,0;245:245.066445,0;246:246.069768,0;247:247.07209,0;248:248.07575,0;249"
	IsotopeLists[95]+=":249.07848,0;"
	IsotopeLists[96]="233:233.0508,0;234:234.05024,0;235:235.05159,0;236:236.05141,0;237:237.05289,0;238:238.05302,0;239:239.05495,0;240:240.055519,0;241:241.0576467,0;242:2"
	IsotopeLists[96]+="42.0588293,0;243:243.0613822,0;244:244.0627463,0;245:245.0654856,0;246:246.0672176,0;247:247.070347,0;248:248.072342,0;249:249.075947,0;250:250.078351,"
	IsotopeLists[96]+="0;251:251.082278,0;252:252.08487,0;"
	IsotopeLists[97]="235:235.05658,0;236:236.05733,0;237:237.05713,0;238:238.05827,0;239:239.05836,0;240:240.05975,0;241:241.06022,0;242:242.06205,0;243:243.063002,0;244:24"
	IsotopeLists[97]+="4.065168,0;245:245.0663554,0;246:246.06867,0;247:247.070299,0;248:248.07308,0;249:249.07498,0;250:250.078311,0;251:251.080753,0;252:252.0843,0;253:253."
	IsotopeLists[97]+="08688,0;254:254.0906,0;"
	IsotopeLists[98]="237:237.06207,0;238:238.06141,0;239:239.06258,0;240:240.0623,0;241:241.06372,0;242:242.06369,0;243:243.06542,0;244:244.06599,0;245:245.06804,0;246:246."
	IsotopeLists[98]+="0687988,0;247:247.070992,0;248:248.072178,0;249:249.074847,0;250:250.0764,0;251:251.07958,0;252:252.08162,0;253:253.085127,0;254:254.087316,0;255:255.0"
	IsotopeLists[98]+="9104,0;256:256.09344,0;"
	IsotopeLists[99]="240:240.06892,0;241:241.06866,0;242:242.0697,0;243:243.06963,0;244:244.07097,0;245:245.07132,0;246:246.07297,0;247:247.07365,0;248:248.07546,0;249:249."
	IsotopeLists[99]+="07641,0;250:250.07865,0;251:251.079984,0;252:252.08297,0;253:253.084818,0;254:254.088016,0;255:255.090266,0;256:256.09359,0;257:257.09598,0;"
	IsotopeLists[100]="242:242.07343,0;243:243.07451,0;244:244.07408,0;245:245.07538,0;246:246.07528,0;247:247.07682,0;248:248.077184,0;249:249.07902,0;250:250.079515,0;251:2"
	IsotopeLists[100]+="51.081566,0;252:252.08246,0;253:253.085176,0;254:254.086848,0;255:255.089955,0;256:256.091767,0;257:257.095099,0;258:258.09707,0;259:259.10059,0;"
	IsotopeLists[101]="245:245.08102,0;246:246.08193,0;247:247.0818,0;248:248.08291,0;249:249.083,0;250:250.08449,0;251:251.08492,0;252:252.08663,0;253:253.08728,0;254:254.08"
	IsotopeLists[101]+="973,0;255:255.091075,0;256:256.09405,0;257:257.095535,0;258:258.098425,0;259:259.1005,0;260:260.10365,0;"
	IsotopeLists[102]="249:249.08782,0;250:250.08749,0;251:251.08896,0;252:252.088966,0;253:253.09065,0;254:254.090949,0;255:255.093232,0;256:256.094276,0;257:257.09685,0;258"
	IsotopeLists[102]+=":258.0982,0;259:259.10102,0;260:260.10264,0;261:261.10574,0;262:262.10752,0;"
	IsotopeLists[103]="251:251.09436,0;252:252.09533,0;253:253.09526,0;254:254.09659,0;255:255.09677,0;256:256.09876,0;257:257.09961,0;258:258.10188,0;259:259.10299,0;260:260"
	IsotopeLists[103]+=".10557,0;261:261.10694,0;262:262.10969,0;263:263.11139,0;"
	IsotopeLists[104]="253:253.10068,0;254:254.10017,0;255:255.10149,0;256:256.10118,0;257:257.10307,0;258:258.10357,0;259:259.10563,0;260:260.10643,0;261:261.10875,0;262:262"
	IsotopeLists[104]+=".10992,0;263:263.11254,0;264:264.11398,0;"
	IsotopeLists[105]="255:255.1074,0;256:256.10811,0;257:257.10786,0;258:258.10944,0;259:259.10972,0;260:260.11143,0;261:261.11211,0;262:262.11415,0;263:263.11508,0;264:264."
	IsotopeLists[105]+="11747,0;265:265.11866,0;"
	IsotopeLists[106]="258:258.11315,0;259:259.11465,0;260:260.11444,0;261:261.1162,0;262:262.11648,0;263:263.11831,0;264:264.11892,0;265:265.12107,0;266:266.12193,0;"
	IsotopeLists[107]="260:260.1218,0;261:261.1218,0;262:262.12301,0;263:263.12315,0;264:264.12473,0;265:265.1252,0;266:266.12701,0;267:267.12774,0;"
	IsotopeLists[108]="263:263.12871,0;264:264.12841,0;265:265.13,0;266:266.13004,0;267:267.13177,0;268:268.13216,0;269:269.13411,0;277:nan,0;"
	IsotopeLists[109]="265:265.13657,0;266:266.13794,0;267:267.13753,0;268:268.13882,0;269:269.13911,0;270:270.14072,0;271:271.14123,0;"
	IsotopeLists[110]="267:267.14396,0;268:268.14353,0;269:269.14514,0;270:270.14463,0;271:271.14608,0;272:272.14631,0;273:273.14925,0;281:nan,0;"
	IsotopeLists[111]="272:272.15348,0;"
	IsotopeLists[112]="285:nan,0;"
	IsotopeLists[113]=""
	IsotopeLists[114]="289:nan,0;"
	IsotopeLists[115]=""
	IsotopeLists[116]="292:nan,0;"
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
	String/G root:Packages:Elements:symbols = ELEMENT_Symbols	// use of this global is DEPRECATED
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
