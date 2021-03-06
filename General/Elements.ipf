#pragma rtGlobals=2		// Use modern global access method.
#pragma IgorVersion = 4.0
#pragma version = 2.14
#pragma ModuleName=elements
#if strlen(WinList("LaueGoFirst.ipf",";","INDEPENDENTMODULE:1"))
#include "MaterialsLocate"						// used to find the path to the materials files, moved to ElementDataInitPackage()
#endif
#initFunctionName "ElementDataInitPackage()"
Constant MIN_LINE_SEPARATION_FRACTION = 0.15	// you can over ride this in your main procedure window.
Constant ELEMENT_Zmax = 118
Constant ELEMENT_MAX_N_EMISSION = 20

Static strConstant edgeTypes = "K;L1;L2;L3;M1;M2;M3;M4;M5;N1;N2;N3;N4;N5;N6;N7;O1;O2;O3;O4;O5;P1;P2;P3"
Static strConstant electronStates = "1s;2s;2p 1/2;2p 3/2;3s;3p 1/2;3p 3/2;3d 3/2;3d 5/2;4s;4p 1/2;4p 3/2;4d 3/2;4d 5/2;4f 5/2;4f 7/2;5s;5p 1/2;5p 3/2;5d 3/2;5d 5/2;6s;6p 1/2;6p 3/2"
Static strConstant emissionTypes = "Ka1;Ka2;Ka1,2;Kb1;Kb2;Kb3;L1;La1;La2;La1,2;Lb1;Lb2;Lg1;Lb2,15;Ma1;"

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
//
//	Nov 17, 2015		2.00
//		ALL other element data now comes from data/elementData.xml
//
//	May 3, 2015			2.01
//		moved #include "MaterialsLocate" to the ElementDataInitPackage() funciton
//
//	Jun 2, 2017			2.04
//		modified EmissionEnergies(), It can now return the average K-emission, or L, or M, or even <Ka>
//
//	Dec 1, 2017			2.06
//		added ChemFormula2IgorStr(), returns a chemical formula suitable for a graph annotation, versions for Igor 6 & 7 
//
//	Mar 12, 2018		2.07
//		added line:   #initFunctionName "ElementDataInitPackage()"
//
//	Mar 19, 2018		2.08
//		now also reads in valenceList from the xml file
//
//	Jul 2, 2018			2.09
//		fixed ProcessMTLfileContentsXML()
//
//	Jul 3, 2018			2.10
//		changed ProcessMTLfileContentsXML(), it now can use fraction or massfraction (fraction is atomic)
//
//	Aug 11, 2018		2.11
//		added electronStates (that match up with edgeTypes)
//
//	Aug 12, 2018		2.12
//		added Function MakePeriodicTablePanel(list)
//
//	Aug 18, 2018		2.13
//		changed Make_ElementDataWaves() and  Make_IsotopesList(), changed the call to XMLtagContents() to speed things up.
//
//	Aug 12, 2018		2.14
//		changed ProcessMTLfileContentsXML() it can now use <chemical_formula>, also fraction is now either mass or atomic.

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


Function ElementData(symb,property,[printIt])
	String symb					// atomic symbol
	String property			// item to return
	Variable printIt
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? strlen(GetRTStackInfo(2))==0 : printIt

	String keys = "Z;amu;density;valence;valenceList;state;firstIon;heatFusion;heatVaporization;thermalConductivity;specificHeat;"
	keys += "atomRadius;meltPoint;boilPoint;covRadius;elecConduc;electroneg;atomVol;DebyeT;name"
	keys = LowerStr(keys)
	property = LowerStr(property)
	Variable ip = WhichListItem(property,keys)		// index to property
	ip = (ip>=0) ? ip+1 : NaN
	Variable Z = element2Z(symb)
	String propertys = "atomic number;atomic mass;density;valence;valenceList;state;first ionization energy;heat of fusion;"
	propertys += "heat of vaporization;thermal conductivity;specific heat;atomic radius;melting point;"
	propertys += "boiling point;covalent radius;electrical conductivity;electronegativity;atomic volume;DebyeT (K);name"
	if ((numtype(Z) || numtype(ip)) && printIt)
		Prompt ip, "property to return",popup, propertys
		Prompt symb, "atomic symbol", popup, ELEMENT_Symbols
		DoPrompt "select", symb,ip
		if (V_flag)
			return NaN
		endif
		property = StringFromList(ip-1,propertys)
		Z = element2Z(symb)
		printIt = 1
	endif
	if (numtype(Z) || numtype(ip))
		return NaN
	endif
	if (printIt)
		printf "ElementData(\"%s\",\"%s\")\r",symb,property
	endif

	String unit="", out="", key=StringFromList(ip-1,keys)
	Variable val=NaN
	if (ip==1)
		val = Z
		Wave wFloat = $""
	elseif (ip==5 || ip==19)
		Wave/T wStr = $("root:Packages:Elements:"+key)
		out = wStr[Z]
		if (ip==5)
			val = WhichListItem(out,"s;l;g")
		endif
		unit = WaveUnits(wStr,-1)
	else
		Wave wFloat = $("root:Packages:Elements:"+key)
		val = wFloat[Z]
		unit = WaveUnits(wFloat,-1)
	endif

	if (printIt)
		unit = SelectString(strlen(unit),"", " ("+unit+")")
		printf "for %s, the %s is %s%s\r",symb,StringFromList(ip-1,propertys),SelectString(strlen(out),num2str(val),"'"+out+"'"),unit
	endif
	return val
End



//  ======================================================================================  //
//  ============================= Start of chemical formuals =============================  //

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
	String symb			// element symbol
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
	String symb						// atomic symbol
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
		elseif (ascii==32 || ascii==9)	// found whitespace
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
			i0 = i-1								// got the number after a symbol, find char after the number
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

	String b0=str[istart]	// starting bracket, can be "(" "[", or "{"
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


#if (IgorVersion()>7)
Function/T ChemFormula2IgorStr(formula)	// returns a chemical formula suitable for a graph annotation
	String formula				// somthing like Ca10(PO4)6(OH)2

	formula = ReplaceString("(", formula, " (")
	formula = TrimBoth(formula)
	String c
	Variable i
	do
		c = formula[i]
		if (isdigit(c))
			formula[i,i] = num2char(8320+str2num(c))	// change digit to sub-script digit
		endif
		i += 1
	while (i<strlen(formula))
	do																// change all multiple spaces to a single space
		formula = ReplaceString("  ", formula, " ")
	while (strsearch(formula,"  ",0) >= 0)
	return formula
End
#else
Function/T ChemFormula2IgorStr(formula)	// returns a chemical formula suitable for a graph annotation
	String formula				// somthing like Ca10(PO4)6(OH)2

	formula = ReplaceString("(", formula, " (")
	String c
	Variable i, onDigit=0
	do
		c = formula[i]
		if (isdigit(c) || CmpStr(c, ".")==0)
			if (!onDigit)
				onDigit = 1
				formula[i] = "\\B"		// this inserts the string, does not replace
				i += 1
			endif	
		else
			if (onDigit)
				onDigit = 0
				formula[i] = "\\M"		// this inserts the string, does not replace
				i += 1
			endif	
		endif
		i += 1
	while (i<strlen(formula))
	if (onDigit)
		formula += "\\M"
	endif	
	do										// change all multiple spaces to a single space
		formula = ReplaceString("  ", formula, " ")
	while (strsearch(formula,"  ",0) >= 0)
	formula = TrimBoth(formula)
	return formula
End
#endif


Function/T LookUpMaterial(material)
	String material				// name of the material, look in file materials/material.mtl

	Variable refNum
	String list						// results

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
	String name = StringByKey("name", XMLattibutes2KeyList("chemical_mixture",buf),"="), list=""
	if (strlen(name))
		list = ReplaceStringByKey("name",list,name,"=")
	endif

	String mix = XMLtagContents("chemical_mixture",buf)
	Variable density=str2num(XMLtagContents("density",mix))
	String unit = StringByKey("unit", XMLattibutes2KeyList("density",mix),"=")
	unit = SelectString(strlen(unit),"g/cm^3",unit)	// no unit specified, set to default
	if (!StringMatch(unit,"g/cm^3"))						// only understand density in g/cm^3
		printf "ERROR -- Invalid units for density = '%s'\r",unit
		density = NaN
	endif
	if (numtype(density)==0 && density>0)
		list = ReplaceNumberByKey("density",list,density,"=")
	endif

	String formula = XMLtagContents("chemical_formula",mix)
	if (strlen(formula))
		list += "formula="+formula+";"
		return list													// have the formula, no need to look further
	endif

	String item="x", nodeList=XMLNodeList(mix)
	Variable m=0, Npart=0
	for (m=0;strlen(item);m+=1)								// count the number of <part> tags
		item = StringFromList(m,nodeList)
		Npart += cmpstr(item,"part")==0
	endfor

	Make/N=(Npart)/D/FREE fractions=0
	Make/N=(Npart)/T/FREE symbols=""
	String part, symbol
	Variable i=0, fraction, Z, massFlag=0, atomicFlag=0
	do
		part = XMLattibutes2KeyList("part",mix,occurance=i)
		if (strlen(part)<1)
			break
		endif
		symbol = StringByKey("symbol",part,"=")
		if (strlen(symbol)<1)
			Z = NumberByKey("Z",part,"=")
			symbol = StringFromList(Z-1,ELEMENT_Symbols)
		endif
		if (strlen(StringByKey("atomic",part,"=")))
			fraction = NumberByKey("atomic",part,"=")	// atomic fraction
			atomicFlag = 1
		else
			fraction = NumberByKey("mass",part,"=")		// mass fraction
			fraction /= Element_amu( WhichListItem(symbol,ELEMENT_Symbols)+1 )	// convert mass fraction --> atomic
			massFlag = 1
		endif
		if (numtype(fraction)==0)
			fractions[i] = fraction
			symbols[i] = symbol
		endif
		i += 1
	while (strlen(symbol) && numtype(fraction)==0)
	if (massFlag && atomicFlag)								// cannot mix atomic and mass fractions
		return ""
	endif

	Variable total = sum(fractions)							// normalize so the total fraction is 1
	fractions /= total
	formula = ""
	for (i=0;i<Npart;i+=1)
		formula += symbols[i]
		fraction = fractions[i]
		if (fraction!=1)
			formula += elements#niceFloatStr(fraction)
		endif
	endfor

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
	line = ReplaceString(" ",line,",")		// replace all spaces with commas
	line = ReplaceString(";",line,",")		// replace all semi-colons with commas
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
	Variable av = abs(v)				// av is now >0, not negative, not 0

	String fmt, str=""
	Variable after = max(placesOfPrecision(av)-floor(log(av))-1, 0)
	sprintf fmt,"%%.%df",after
	sprintf str,fmt,v
	return str
End

//  ============================== End of chemical formuals ==============================  //
//  ======================================================================================  //



//  ======================================================================================  //
//  ======================= Start of Convienent element functions ========================  //


Function element2Z(symb)		// returns Z for an atomic symbol (NOT case sensitive)
	String symb						// atomic symbol
	symb[0,0] = UpperStr(symb[0,0])		// ensure first char is upper case
	symb[1,1] = LowerStr(symb[1,1])		// and second character is lower
	Variable iz = WhichListItem(symb,ELEMENT_Symbols)+1
	if (!(iz>0))
		iz = str2num(symb)
	endif
	return ((iz>0) ? iz : NaN)
End


Function Element_amu(Z)
	Variable Z
	Wave amu = root:Packages:Elements:amu
	return amu[Z]
End

Function Element_density(Z)
	Variable Z
	Wave density = root:Packages:Elements:density
	return density[Z]
End

Function Element_valence(Z)
	Variable Z
	Wave valence = root:Packages:Elements:valence
	return valence[Z]
End

Function/T Element_valenceList(Z)
	Variable Z
	Wave/T valenceList = root:Packages:Elements:valenceList
	return valenceList[Z]
End

Function/T Element_state(Z)
	Variable Z
	Wave/T state = root:Packages:Elements:state
	return state[Z]
End

Function Element_firstIon(Z)
	Variable Z
	Wave firstIon = root:Packages:Elements:firstIon
	return firstIon[Z]
End

Function Element_heatFusion(Z)
	Variable Z
	Wave heatFusion = root:Packages:Elements:heatFusion
	return heatFusion[Z]
End

Function Element_heatVapor(Z)
	Variable Z
	Wave heatVaporization = root:Packages:Elements:heatVaporization
	return heatVaporization[Z]
End

Function Element_thermConduc(Z)
	Variable Z
	Wave thermalConductivity = root:Packages:Elements:thermalConductivity
	return thermalConductivity[Z]
End

Function Element_specificHeat(Z)
	Variable Z
	Wave specificHeat = root:Packages:Elements:specificHeat
	return specificHeat[Z]
End

Function Element_atomRadius(Z)
	Variable Z
	Wave atomRadius = root:Packages:Elements:atomRadius
	return atomRadius[Z]
End

Function Element_meltPoint(Z)
	Variable Z
	Wave meltPoint = root:Packages:Elements:meltPoint
	return meltPoint[Z]
End

Function Element_boilPoint(Z)
	Variable Z
	Wave boilPoint = root:Packages:Elements:boilPoint
	return boilPoint[Z]
End

Function Element_covRadius(Z)
	Variable Z
	Wave covRadius = root:Packages:Elements:covRadius
	return covRadius[Z]
End

Function Element_elecConduc(Z)
	Variable Z
	Wave elecConduc = root:Packages:Elements:elecConduc
	return elecConduc[Z]
End

Function Element_electroneg(Z)
	Variable Z
	Wave electroneg = root:Packages:Elements:electroneg
	return electroneg[Z]
End

Function Element_atomVol(Z)
	Variable Z
	Wave atomVol = root:Packages:Elements:atomVol
	return atomVol[Z]
End

Function Element_DebyeT(Z)	// Debye Temperature (K)
	// room temperature values of DebyeT from: http://www.knowledgedoor.com/2/elements_handbook/debye_temperature.html
	// Al value from Acta Cryst. (1972). A28, 22-27����[ doi:10.1107/S056773947200004X ]
	Variable Z
	Wave DebyeT = root:Packages:Elements:DebyeT
	return DebyeT[Z]
End

Function/T Element_names(Z)
	Variable Z
	Wave/T name = root:Packages:Elements:name
	return name[Z]
End

Function BindingEnergy(symb,edgeType)
	String symb				// atomic symbol
	String edgeType		// item to return

	Variable topLevel = ItemsInList(GetRTStackInfo(0))<2
	if ((strlen(symb)<1 || strlen(edgeType)<1) && topLevel)
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


Function/T edgeType2electronState(type)	// takes "L2" returns 2p1/2
	String type
	Variable i = WhichListItem(type, edgeTypes,";",0,0)
	if (i<0)
		return ""
	endif
	return StringFromList(i,electronStates)
End

//  ======================== End of Convienent element functions =========================  //
//  ======================================================================================  //



//  ======================================================================================  //
//  ========================== Start of emission line functions ==========================  //

Function EmissionEnergies(symb,linetype,[printIt])		// display or return am average or a single emission line energy
	String symb				// atomic symbol
	String linetype		// item to return
	Variable printIt
	Variable topLevel=strlen(GetRTStackInfo(2))==0
	printIt = ParamIsDefault(printIt) || numtype(printIt) ? topLevel : printIt
	linetype = ReplaceString("<",linetype,"")
	linetype = ReplaceString(">",linetype,"")

	Wave/T FullEmissionLineInfo = root:Packages:Elements:FullEmissionLineInfo
	if (!WaveExists(FullEmissionLineInfo))
		ElementDataInitPackage()
	endif

	String emissionTypesMenu="All;<K>;<Ka>;<Kb>;<Kg>;Ka1;Ka2;Ka1,2;Kb1;Kb2;Kb3;<L>;<La>;<Lb>;<Lg>;L1;La1;La2;La1,2;Lb1;Lb2;Lg1;Lb2,15;Ma1;"

	Variable Z=element2Z(symb)
	if (topLevel && (numtype(Z) || !edgeTypeInLines(linetype,emissionTypes)))
		Prompt linetype, "emission line type",popup, emissionTypesMenu
		Prompt symb, "atomic symbol", popup, ELEMENT_Symbols
		DoPrompt "select", symb,linetype
		if (V_flag)
			return NaN
		endif
		linetype = ReplaceString("<",linetype,"")
		linetype = ReplaceString(">",linetype,"")
		printIt = 1
	endif
	Z = element2Z(symb)
	if (numtype(Z) || !edgeTypeInLines(linetype,emissionTypes))
		return NaN
	endif

	STRUCT EmissionLineStruct em
	StructGet/S em, FullEmissionLineInfo[Z]

	Variable i, eV, rel
	if (stringmatch(linetype,"All"))		// print out all the emission lines for symb
		printEmissionLineStruct(em)
		if (em.N >0)
			eV = em.line[em.N - 1].eV			// return the highest energy
		endif

	else												// looking for a single or an average
		Variable relSum, N
		for (i=0,eV=0,relSum=0,N=0; i<em.N; i+=1)
			if (strsearch(em.line[i].name, linetype,0,2)==0)	// starts with linetype
				rel = em.line[i].rel
				eV += em.line[i].eV * rel		// accumulate for average
				relSum += rel
				N += 1
			endif
		endfor
		eV = N ? eV/relSum : NaN
		if (printIt && numtype(eV)==0)
			String str=""
			if (N>1)
				sprintf str, ",  the weighted average of %d lines",N
			endif	
			printf "  %s(%s) is at %g eV  (relative strength = %g%s)\r",symb,linetype,eV,relSum,str
		endif
	endif
	return eV
End
//
Static Function edgeTypeInLines(edgeType,lines)
	String edgeType
	String lines
	if (cmpstr(edgeType,"All")==0)
		return 1
	endif
	Variable i, N=ItemsInList(lines)
	for (i=0;i<N;i+=1)
		if (strsearch(StringFromList(i,lines),edgeType,0,2)==0)
			return 1
		endif
	endfor
	return 0
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

Static Function/T emissionLineList(Z,[minSep])			// get list of lines to fit
	Variable Z
	Variable minSep						// minimum separation between lines for grouping
	minSep = ParamIsDefault(minSep) || numtype(minSep) || minSep<0 ? NumVarOrDefault("root:Packages:Elements:defaultPeakFWHM",0)*MIN_LINE_SEPARATION_FRACTION : minSep
	// retruns a ";" separated list, each element is of form "eV:strength:name"
	// where the strength is relative to 100, and name is the name of the emission line

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



Structure EmissionLineStruct
	char	symb[2]
	int16	Z
	int16	N			// number of emission lines stored
	STRUCT oneEmissionLineStruct line[ELEMENT_MAX_N_EMISSION]
EndStructure
//
Static Structure oneEmissionLineStruct
	char	name[6]	// name of line, e.g. "Ka1"
	double	eV		// energy of line
	double	rel	// relative strength
EndStructure
//
Function printEmissionLineStruct(em)
	STRUCT EmissionLineStruct &em
	printf "for '%s' (%g):\r",em.symb,em.Z
	if (em.N<1)
		print "\tempty\r"
	else
		Variable i, N=min(em.N,ELEMENT_MAX_N_EMISSION)
		for (i=0;i<em.N;i+=1)
			printf "\t%s \t%g eV \t%g\r",em.line[i].name,em.line[i].eV,em.line[i].rel
		endfor
	endif
End
//
Static Function initEmissionLineStruct(em)
	STRUCT EmissionLineStruct &em
	em.Z = 0
	em.N = 0
	Variable i
	for (i=0;i<ELEMENT_MAX_N_EMISSION;i+=1)
		em.line[i].name = ""
		em.line[i].eV = NaN
		em.line[i].rel = NaN
	endfor
End

//  =========================== End of emission line functions ===========================  //
//  ======================================================================================  //



//  ======================================================================================  //
//  ========================== Start of Find Element From Line ===========================  //

Function/T FindBestElementFromEmissionLine(eV,acceptableLines,[dE,excitation])	// searches for the right element
	Variable eV
	String acceptableLines	// acceptable lines, check all if "",  examples: "Ka*;Kb*", or "K*", "", "Lb1",  or "*L*"
	Variable dE					// HW resolution, find strongest line within eV�dE, for dE=0, just use closest line
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
			if (!acceptableType(acceptableLines,line))		// if line is not of desired type, continue
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
Static Function acceptableType(acceptableTypes,type)	// returns 1 if type is in acceptableTypes
	String acceptableTypes, type
	if (strlen(acceptableTypes)<1)								// if no types specified, accept all
		return 1
	endif
	Variable i
	for (i=0;i<ItemsInList(acceptableTypes);i+=1)
		if (StringMatch(type,StringFromList(i,acceptableTypes)))
			return 1
		endif
	endfor
	return 0
End
//Static Function acceptableLine(acceptableLines,line)	// returns 1 if line is in acceptableLines
//	String acceptableLines, line
//	if (strlen(acceptableLines)<1)								// if no lines specified, accept all
//		return 1
//	endif
//	Variable i
//	for (i=0;i<ItemsInList(acceptableLines);i+=1)
//		if (StringMatch(line,StringFromList(i,acceptableLines)))
//			return 1
//		endif
//	endfor
//	return 0
//End


Function/T FindBestElementFromEdge(eV,acceptableEdges,[symbols,dE,all])	// searches for the right element
	// returns symbol:Z:edge:eV   e.g.   "Ti:22:L2:460.2"
	Variable eV
	String acceptableEdges	// acceptable edges, check all if "",  examples: "Ka*;Kb*", or "K*", "", "Lb1",  or "*L*"
	String symbols				// list of allowed element symbols, "" means all
	Variable dE					// HW resolution, find strongest edge within eV�dE, for dE=0, just use closest edge
	Variable all				// return list of all elements within �dE (dE is required for all)
	symbols = SelectString(ParamIsDefault(symbols),symbols,"")
	dE = ParamIsDefault(dE) || numtype(dE) || dE<0? 0 : dE
	all = ParamIsDefault(all) || numtype(all) || dE<=0 ? 0 : all
	acceptableEdges = ReplaceString(",", acceptableEdges, ";")
	acceptableEdges = ReplaceString(" ", acceptableEdges, ";")
	symbols = ReplaceString(",", symbols, ";")
	symbols = ReplaceString(" ", symbols, ";")

	Wave/T edgeStrings = root:Packages:Elements:edgeStrings
	Make/N=1000/T/FREE edges=""
	Make/N=1000/D/FREE deltas=Inf
	String list, edge, item, edgeBest="", symb, str
	Variable dEmin=Inf,eBest,Zbest=0, Nedges=0
	Variable iedge, Z, eni, deltaE, better
	for (Z=1;Z<=92;Z+=1)											// for each element
		list = edgeStrings(Z)
		if (strlen(symbols)>0)
			symb = StringFromList(Z-1, ELEMENT_Symbols)
			if (WhichListItem(symb, symbols)<0)
				continue
			endif
		endif
		for (iedge=0;iedge<ItemsInList(list);iedge+=1)	// for each edge in an element
			item = StringFromList(iedge,list)					// a single edge for element Z
			edge = StringFromList(0,item,":")					// name of this edge
			eni = str2num(StringFromList(1,item,":"))		// energy of this edge (eV)
			if (!elements#acceptableType(acceptableEdges,edge))		// if edge is not of desired type, continue
				continue
			endif

			deltaE = abs(eni-eV)
			if (all && deltaE<dE)
				sprintf str, "%s:%d:%s:%g", StringFromList(Z-1,ELEMENT_Symbols),Z,edge,eni
				edges[Nedges] = str
				deltas[Nedges] = eV-eni
				Nedges += 1
				if (Nedges>=1000)
					break
				endif
			endif
			better = (deltaE<dEmin)								// eni is close enough, look for best edge, dE>0
//			better = better || (dE==0 && deltaE<dEmin)		// find closest edge, no width specified
			if (better)
				dEmin = deltaE
				eBest = eni
				edgeBest = edge
				Zbest = Z
			endif
		endfor
	endfor

	String out=""
	if (!all && Zbest>0)
		sprintf out, "%s:%d:%s:%g:%g", StringFromList(Zbest-1,ELEMENT_Symbols),Zbest,edgeBest,eBest,abs(eBest-eV)
	elseif (all)
		Redimension/N=(Nedges) edges, deltas
		Duplicate/FREE deltas, ABSdeltas
		ABSdeltas = abs(deltas)
		Sort ABSdeltas, edges ,deltas
		Variable i
		for (i=0;i<Nedges;i+=1)
			sprintf str, "%s:%g;", edges[i],deltas[i]
			out += str
		endfor
	endif
	return out
End
//  =========================== End of Find Element From Line ============================  //
//  ======================================================================================  //






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
//  ======================================================================================  //



//  ======================================================================================  //
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
	Variable iZ, m, massNumber,atomicMass,Fraction, start=0
	for (iZ=1;iZ<=Zmax;iZ+=1)
		oneZ = XMLtagContents("Z"+num2istr(iZ),isotopes, start=start)
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

//  ============================== End of isotope functions ==============================  //
//  ======================================================================================  //


//  ======================================================================================  //
//  ======================== Start of Make Element Table Selector ========================  //
//
//		To use MakePeriodicTablePanel(), try:
//
//		String list = "Sr;Ti;O"			// preselected starting values (semi-colon separated)
//		list = MakePeriodicTablePanel(list, wait=1)		// list will contain the user selected values
//
// 	or do not use wait=1, and use:  GetUserData("ElementsPanel","","listON"))
//
//		or make your own:  Function ElementPanelFunc(symb)
//
Function/T MakePeriodicTablePanel(list, [wait, extra])
	String list									// list of elements to start with
	Variable wait
	String extra								// text to store in the extra userdata, e.g. "key=value;"
	wait = ParamIsDefault(wait) || numtype(wait) ? 0 : wait
	extra = SelectString(ParamIsDefault(extra),extra,"")

	Variable N=ItemsInList(ELEMENT_Symbols)
	Make/N=(N,2)/FREE/I rc=-1

	rc[0][0] = {0,0}							// first row, [H,He]
	rc[0][1] = {0,17}

	rc[2][0] = {1,1,1,1,1,1,1,1}			// second row, [Li,Ne]
	rc[2][1] = {0,1,12,13,14,15,16,17}

	rc[10][0] = {2,2,2,2,2,2,2,2}		// third row, [Na,Ar]
	rc[10][1] = {0,1,12,13,14,15,16,17}

	rc[18][0] = {3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3}	// fourth row, [K,Kr]
	rc[18][1] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17}

	rc[36][0] = {4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4}		// fifth row, [Rb,Xe]
	rc[36][1] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17}

	rc[54][0] = {5,5,5}						// sixth row, [Cs,Rn]
	rc[54][1] = {0,1,2}
	rc[57][0] = {8,8,8,8,8,8,8,8,8,8,8,8,8,8}
	rc[57][1] = {2,3,4,5,6,7,8,9,10,11,12,13,14,15}
	rc[71][0] = {5,5,5,5,5,5,5,5,5,5,5,5,5,5,5}
	rc[71][1] = {3,4,5,6,7,8,9,10,11,12,13,14,15,16,17}

	rc[86][0] = {6,6,6}						// seventh row, [Fr,Og]
	rc[86][1] = {0,1,2}
	rc[89][0] = {9,9,9,9,9,9,9,9,9,9,9,9,9,9}
	rc[89][1] = {2,3,4,5,6,7,8,9,10,11,12,13,14,15}
	rc[103][0] = {6,6,6,6,6,6,6,6,6,6,6,6,6,6,6}
	rc[103][1] = {3,4,5,6,7,8,9,10,11,12,13,14,15,16,17}

	String symb, name
	Variable i,top,left
	NewPanel /W=(342,65,893,349)/N=ElementsPanel/K=1
	for (i=0;i<N;i+=1)
		left = 5 + 30*rc[i][1]
		top = 5 + 29*rc[i][0]
		top -= top >230 ? 17 : 0			// shift for the Lanthanides
		symb = StringFromList(i,ELEMENT_Symbols)
		name = "button_"+symb
		Button $name,pos={left,top},size={30,30},title=symb, proc=elements#ElementsPanelButtonProc
		if (WhichListItem(symb,list) >= 0)
			Button $name, userdata(ON_OFF)="1", fColor=(65535,32768,32768)
		endif
	endfor
	Button buttonClearAll,pos={201,64},size={60,20},fColor=(65535,65535,65535),proc=elements#ElementsPanelButtonProc,title="Clear All"
	SetWindow ElementsPanel userdata(listON)=list
	SetWindow ElementsPanel userdata(extra)=extra					// store extra information, as text

	if (wait)
		DoWindow /T ElementsPanel, "Close When Done Seclecting Elements..."
		String/G ElementsPanelList_JZT = list
		PauseForUser ElementsPanel
		list = ElementsPanelList_JZT
	endif
	return list
End
//
Static Function ElementsPanelButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	if (ba.eventCode == 2)
		String extra = GetUserData(ba.win,"","extra")
		FUNCREF ElementPanelProtoFunc extraFunc = $StringByKey("func",extra,"=")// assign extra function
		if (cmpstr(ba.ctrlName,"buttonClearAll",0)==0 && cmpstr(ba.win,"ElementsPanel")==0)
			SVAR gList = ElementsPanelList_JZT
			String list, name
			if (SVAR_Exists(gList))
				gList = ""
			endif
			SetWindow $(ba.win) userdata(listON)=""

			list = ControlNameList(ba.win)
			Variable i, N=ItemsInList(list)
			for (i=0;i<N;i+=1)
				name = StringFromList(i,list)
				if (strsearch(name,"button_",0)==0)
				Button $name, win=$(ba.win), userdata(ON_OFF)="0", fColor=(65535,65535,65535)
				endif
			endfor
			extraFunc(ba.win,"")						// call the extra function, defaults to ElementPanelProtoFunc()
			return 0

		elseif (strsearch(ba.ctrlName,"button_",0)!=0 || cmpstr(ba.win,"ElementsPanel")!=0)
			return 0
		endif

		String symb = ReplaceString("button_",ba.ctrlName,"")
		list = GetUserData(ba.win,"","listON" )
		Variable on = str2num(GetUserData(ba.win,ba.ctrlName,"ON_OFF"))
		on = numtype(on) ? 0 : on
		if (on)
			Button $(ba.ctrlName), win=$(ba.win), userdata(ON_OFF)="0", fColor=(65535,65535,65535)
			list = RemoveFromList(symb,list)
		else
			Button $(ba.ctrlName), win=$(ba.win), userdata(ON_OFF)="1", fColor=(65535,32768,32768)
			//	list += symb+";"
			list = AddListItem(symb,list,";",Inf)
		endif
		SetWindow $(ba.win) userdata(listON)=list
		SVAR gList = ElementsPanelList_JZT
		if (SVAR_Exists(gList))
			gList = list
		endif
		extraFunc(ba.win,symb)						// call the extra function, defaults to ElementPanelProtoFunc()
	endif
	return 0
End
//
Function ElementPanelProtoFunc(win,symb)
	String win
	String symb
	// print "in ElementPanelProtoFunc with symb= ",symb
End

//  ========================= End of Make Element Table Selector =========================  //
//  ======================================================================================  //



//  ======================================================================================  //
//  ======================== Start of Load/Initializing Wave data ========================  //

Static Function Make_ElementDataWaves()
	// loads waves of element data from data file "data:elementData.xml"
	Variable f
	String fileName = ParseFilePath(1,FunctionPath("Make_ElementDataWaves"),":",1,0)+"data:elementData.xml"
	Open/R/M="Element info file"/T=".xml"/Z f as fileName
	if (V_flag)
		DoAlert 0, "ERROR -- could not find 'LaueGo:General:data:elementData.xml'"
		return 1
	endif
	FStatus f
	String buf = PadString("",V_logEOF,0x20)
	FBinRead f, buf
	Close f

	String elementData = XMLtagContents("element_xray",buf)
	if (strlen(elementData)<100)
		return 1
	endif
	elementData = XMLremoveComments(elementData)

	Variable Zmax = ItemsInList(ELEMENT_Symbols)
	Make/N=(Zmax+1)/T/O root:Packages:Elements:name/WAVE=name = ""
	Make/N=(Zmax+1)/D/O root:Packages:Elements:amu/WAVE=amu = NaN
	Make/N=(Zmax+1)/D/O root:Packages:Elements:density/WAVE=density = NaN
	Make/N=(Zmax+1)/D/O root:Packages:Elements:valence/WAVE=valence = NaN
	Make/N=(Zmax+1)/T/O root:Packages:Elements:valenceList/WAVE=valenceList = ""
	Make/N=(Zmax+1)/D/O root:Packages:Elements:firstIon/WAVE=firstIon = NaN
	Make/N=(Zmax+1)/D/O root:Packages:Elements:heatFusion/WAVE=heatFusion = NaN
	Make/N=(Zmax+1)/D/O root:Packages:Elements:heatVaporization/WAVE=heatVaporization = NaN
	Make/N=(Zmax+1)/D/O root:Packages:Elements:thermalConductivity/WAVE=thermalConductivity = NaN
	Make/N=(Zmax+1)/D/O root:Packages:Elements:specificHeat/WAVE=specificHeat = NaN
	Make/N=(Zmax+1)/D/O root:Packages:Elements:meltPoint/WAVE=meltPoint = NaN
	Make/N=(Zmax+1)/D/O root:Packages:Elements:boilPoint/WAVE=boilPoint = NaN
	Make/N=(Zmax+1)/D/O root:Packages:Elements:atomRadius/WAVE=atomRadius = NaN
	Make/N=(Zmax+1)/D/O root:Packages:Elements:covRadius/WAVE=covRadius = NaN
	Make/N=(Zmax+1)/D/O root:Packages:Elements:electroneg/WAVE=electroneg = NaN
	Make/N=(Zmax+1)/D/O root:Packages:Elements:atomVol/WAVE=atomVol = NaN
	Make/N=(Zmax+1)/D/O root:Packages:Elements:elecConduc/WAVE=elecConduc = NaN
	Make/N=(Zmax+1)/D/O root:Packages:Elements:DebyeT/WAVE=DebyeT = NaN
	Make/N=(Zmax+1)/T/O root:Packages:Elements:state/WAVE=state = ""
	Make/N=(Zmax+1)/T/O root:Packages:Elements:edgeStrings/WAVE=edgeStrings = ""
	Make/N=(Zmax+1)/T/O root:Packages:Elements:FullEmissionLineInfo/WAVE=FullEmissionLineInfo = ""

	SetScale d 0,0,"amu", amu
	SetScale d 0,0,"g cm^-3", density
	SetScale d 0,0,"V", firstIon
	SetScale d 0,0,"kJ mole^-1", heatFusion
	SetScale d 0,0,"kJ mole^-1", heatVaporization
	SetScale d 0,0,"W m^-1 K^-1", thermalConductivity
	SetScale d 0,0,"J g^-1 K^-1", specificHeat
	SetScale d 0,0,"K", meltPoint
	SetScale d 0,0,"K", boilPoint
	SetScale d 0,0,"�", atomRadius
	SetScale d 0,0,"�", covRadius
	SetScale d 0,0,"cm^3 mole^-1", atomVol
	SetScale d 0,0,"10^6 Ohm^-1 m^-1", elecConduc
	SetScale d 0,0,"K", DebyeT

	Note/K name, "waveClass=ElementInfo_name;"
	Note/K amu, "waveClass=ElementInfo_amu;"
	Note/K density, "waveClass=ElementInfo_density;"
	Note/K valence, "waveClass=ElementInfo_valence;"
	Note/K valenceList, "waveClass=ElementInfo_valenceList;"
	Note/K firstIon, "waveClass=ElementInfo_firstIon;"
	Note/K heatFusion, "waveClass=ElementInfo_heatFusion;"
	Note/K heatVaporization, "waveClass=ElementInfo_heatVaporization;"
	Note/K thermalConductivity, "waveClass=ElementInfo_thermalConductivity;"
	Note/K specificHeat, "waveClass=ElementInfo_specificHeat;"
	Note/K meltPoint, "waveClass=ElementInfo_meltPoint;"
	Note/K boilPoint, "waveClass=ElementInfo_boilPoint;"
	Note/K atomRadius, "waveClass=ElementInfo_atomRadius;"
	Note/K covRadius, "waveClass=ElementInfo_covRadius;"
	Note/K electroneg, "waveClass=ElementInfo_electroneg;"
	Note/K atomVol, "waveClass=ElementInfo_atomVol;"
	Note/K elecConduc, "waveClass=ElementInfo_elecConductivity;"
	Note/K DebyeT, "waveClass=ElementInfo_DebyeT;"
	Note/K state, "waveClass=ElementInfo_state;"
	Note/K edgeStrings, "waveClass=ElementInfo_xrayEdges;"
	Note/K FullEmissionLineInfo, "waveClass=ElementInfo_xrayEmission;"

	STRUCT EmissionLineStruct em
	String oneZ, edges, emission, ename, str, strStruct
	Variable iZ, i, N, eV,rel, start=0
	for (iZ=1;iZ<=Zmax;iZ+=1)
		oneZ = XMLtagContents("Z"+num2istr(iZ),elementData, start=start)

		name[iZ] = XMLtagContents("name",oneZ)
		amu[iZ] = str2num(XMLtagContents("amu",oneZ))
		density[iZ] = str2num(XMLtagContents("density",oneZ))
		valenceList[iZ] = XMLtagContents("valence",oneZ)
		valenceList[iZ] = ReplaceString(" ", valenceList[iZ], ";")
		valence[iZ] = str2num(valenceList[iZ])
//		valence[iZ] = str2num(XMLtagContents("valence",oneZ))
		firstIon[iZ] = str2num(XMLtagContents("firstIon",oneZ))
		heatFusion[iZ] = str2num(XMLtagContents("fusion",oneZ))
		heatVaporization[iZ] = str2num(XMLtagContents("vapor",oneZ))
		thermalConductivity[iZ] = str2num(XMLtagContents("thermalConduc",oneZ))
		specificHeat[iZ] = str2num(XMLtagContents("specificHeat",oneZ))
		meltPoint[iZ] = str2num(XMLtagContents("melt",oneZ))
		boilPoint[iZ] = str2num(XMLtagContents("boil",oneZ))
		atomRadius[iZ] = str2num(XMLtagContents("radius",oneZ))
		covRadius[iZ] = str2num(XMLtagContents("covRadius",oneZ))
		electroneg[iZ] = str2num(XMLtagContents("electroneg",oneZ))
		atomVol[iZ] = str2num(XMLtagContents("atomVol",oneZ))
		DebyeT[iZ] = str2num(XMLtagContents("DebyeT",oneZ))
		elecConduc[iZ] = str2num(XMLtagContents("elecConduc",oneZ))
		state[iZ] = XMLtagContents("state",oneZ)

		edges = XMLtagContents("edges",oneZ)
		N = str2num(XMLtagContents("N",edges))
		N = N > 0 ? N : 0
		for (i=0;i<N;i+=1)
			str = XMLtagContents("edge",edges,occurance=i)
			edgeStrings[iZ] += ReplaceString(" ",str,":")+";"
		endfor

		emission = XMLtagContents("emissionLines",oneZ)
		N = str2num(XMLtagContents("N",emission))
		N = N > 0 ? N : 0
		initEmissionLineStruct(em)
		em.symb = StringFromList(iZ-1,ELEMENT_Symbols)
		em.Z = iZ
		em.N = N
		for (i=0;i<N;i+=1)
			str = XMLtagContents("line",emission,occurance=i)
			sscanf str, "%s %g %g",ename,eV,rel
			if (V_flag==3)
				em.line[i].name = ename
				em.line[i].eV = eV
				em.line[i].rel = rel
			endif
		endfor
		StructPut/S em, strStruct
		FullEmissionLineInfo[iZ] = strStruct
	endfor
	return 0
End

//  ========================= End of Load/Initializing Wave data =========================  //
//  ======================================================================================  //





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

	Make_ElementDataWaves()
	if (Exists("root:Packages:Elements:defaultPeakFWHM")!=2 || !ParamIsDefault(fwhm))
		fwhm = (numtype(fwhm) || fwhm<0 || ParamIsDefault(fwhm)) ? 0 : fwhm
		Variable /G root:Packages:Elements:defaultPeakFWHM=fwhm
	endif
	Make_IsotopesList()

	if (!exists("MaterialsAreHere"))
		// this provides the MaterialsAreHere() that locates the default location for materials
		//	it is done this way incase MaterialsLocate.ipf file does not exist.
		//	I would rather just put the include at the top of this file.
		Execute/P/Q/Z "INSERTINCLUDE \"MaterialsLocate\""
		Execute/P/Q/Z "COMPILEPROCEDURES "
	endif
End
