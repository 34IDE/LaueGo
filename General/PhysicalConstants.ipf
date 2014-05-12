#pragma rtGlobals=1							// turns�on�runtime�lookup�of�globals
#pragma ModuleName=PhysicalConstants
#pragma version = 2.04
#pragma IgorVersion = 6.11

//	By Jon Tischler (ORNL)  Aug 12, 2010

Static StrConstant PhysicalConstantServer="http://physics.nist.gov/cuu/Constants/Table/allascii.txt"

// Note, this list of constants is note really needed, it is just a convienence for copying and pasting into other experiments.
//
// Set to 2006 CODATA values, Aug 11, 2010
//  from	http://physics.nist.gov/cuu/Constants
Constant c_ms		= 299792458				// speed of light (m/s) (exact)
Constant eo_Fm		= 8.854187817e-12		// permitivity of vacuum (F/m)
Constant G_MKS		= 6.67428E-11			// Gravitational constant (m^3 1/kg 1/s)
Constant h_Js		= 6.62606896E-34		// Plank constant (J s)
Constant h_eVs		= 4.13566733e-15		// Plank constant (eV s)
Constant e_C		= 1.602176487e-19		// Charge on electron (C)
Constant alpha_fine = 7.2973525376e-3	// fine structure constant
Constant aB_m		= 0.52917720859e-10	// Bohr radius (m)
Constant Rydberg	= 13.60569193			// Rydberg (eV) ( = alpha^2*me*c/2/h * hc)
Constant me_kg		= 9.10938215e-31		// mass of electron (kg)
Constant me_keV	= 510.998910			// mass of electron (keV)
Constant me_mp		= 5.4461702177e-4		// mass ratio m(electron)/m(proton)
Constant mp_kg		= 1.672621637e-27		// mass of proton (kg)
Constant mn_kg		= 1.674927211e-27		// mass of neutron (kg)
Constant mn_eV		= 9.39565346e8			// mass of neutron (eV)
Constant NA			= 6.02214179e23		// Avagadro's number
Constant kB_JK		= 1.3806504e-23		// Boltzman constant (J/K)
Constant kB_eVK	= 8.617343e-5			// Boltzman constant (eV/K)
Constant J_eV		= 1.602176487e-19		// Joules/eV
Constant atm_Pa	= 101325					// Number of Pa in one atm (exact)
Constant hc_keVA	= 12.39841930			// h*c (keV-�)
Constant hc_keVnm = 1.239841930			// h*c (keV-nm)
Constant re_m		= 2.8179403267e-15	// Thompson radius (m)
Constant re_A		= 2.8179403267e-5		// Thompson radius (�)
Constant CuKa1		= 1.540593226			// wavelength of Cu Kalpha(1) (�)
Constant MoKa1		= 0.709317155			// wavelength of Mo Kalpha(1) (�)
Constant aSi_A		= 5.43102064			// lattice constant of Si at 22.5� (�), 2006 CODATA,  alpha=2.56E-6 (deg/K)


Menu "Analysis"
	SubMenu "Physical Constants"
		"New Static Constant, Physical Constant...",PhysicalConstant_InsertStatic("*",web=0)
		"<BGet a Physical Constant...",LookUpPhysicalConstant("*")
		"<I  update your local copy [Rarely needed]",UpdateLocalCopyOfConstants()
	End
End



Static Constant PhysicalConstantMaxStr=60

Static Structure PhysicalConstantStructure
	int16	valid
	char name[PhysicalConstantMaxStr+1]
	double value
	double err
	char unit[PhysicalConstantMaxStr+1]
	int16 exact
EndStructure
//
Function copyPhysicalConstantStructure(f,i)
	STRUCT PhysicalConstantStructure &f, &i
	f.valid	= i.valid
	f.name	= i.name
	f.value	= i.value
	f.err		= i.err
	f.unit	= i.unit
	f.exact	= i.exact
End



Function PhysicalConstant_InsertStatic(name,[web,printIt])
	String name
	Variable printIt
	Variable web							// 1 try to use web first, 0 only use local
	if (ParamIsDefault(printIt))
		printIt = strlen(GetRTStackInfo(2)) ? printIt : 1
	endif
	web = ParamIsDefault(web) ? 0 : web
	if (strlen(name)<1)
		return NaN
	endif
	STRUCT PhysicalConstantStructure c
	LookUpPhysicalConstant(name,c=c,web=web,printIt=0)
	if (!(c.valid))
		return NaN
	endif
	name = c.name
	String unit = c.unit
	Variable value = c.value

	String conversions="inverse meter-electron volt relationship=hc_keVnm,1e6,h*c (keV-nm):hc_keVA,1e7,h*c (keV-�);"
	conversions += "classical electron radius=re_m,1,Thompson radius (m):re_nm,1e9,Thompson radius (nm):re_A,1e10,Thompson radius (�);"
	conversions += "speed of light in vacuum=c_ms,1:c_nms,1e9;"
	conversions += "{220} lattice spacing of silicon=aoSi022_m,1:aoSi022_nm,1e9:aoSi022_A,1e10;"
	conversions += "electron mass energy equivalent in MeV=me_eV,1e6:me_keV,1e3:me_keV,1;"
	conversions += "Avogadro constant=NA,1;"

	String lists=StringByKey(name,conversions,"="), outName, comment=""
	Variable N=ItemsInList(lists,":"),i
	if (N<1)
		outName = CleanupName(name,0)
		unit = CleanupName(c.unit,0)
		i = strlen(unit)
		outName = outName[0,30-i-1]
		for (i=strlen(outName)-1; stringmatch(outName[i],"_") && i>=0; i-=1)	// find last non-underscore
		endfor
		outName = outName[0,i]+"_"+unit
	elseif (N==1)
		outName = StringFromList(0,lists,",")
		value *= str2num(StringFromList(1,lists,","))
		comment = StringFromList(2,lists,",")
	else
		String choice
		lists = ReplaceString(":",lists,";")
		Prompt choice,"desired output",popup, lists 
		DoPrompt "output",choice
		if (V_flag)
			return 1
		endif
		outName = StringFromList(0,choice,",")
		value *= str2num(StringFromList(1,choice,","))
		comment = StringFromList(2,choice,",")
	endif
	outName = CleanupName(outName,0)
	outName = ReplaceString("__",outName,"_")
	comment = SelectString(strlen(comment),"", "\t\t\t// "+comment)

	if (1)
		printf "from \"%s\"  =  %.13g  (%s)\r",c.name,c.value,c.unit
		printf "insert:\t\tStatic Constant %s = %.13g%s\r",outName,value,comment
		return value
	else
		String cmd
		sprintf cmd, "INSERTINCLUDE \"Static Constant %s = %.15g%s\"", outName,value,comment
		DoAlert/T="Insert Constant\r" 1, "Insert line\r"+cmd
		V_flag = 1
		if (V_flag==1)
			printf "from \"%s\"  =  %.13g  (%s)\r",c.name,c.value,c.unit
			printf "insert:\t\tStatic Constant %s = %.13g%s\r",outName,value,comment
			Execute/P cmd
		else
			print "nothing done"
		endif
	endif
	return 0
End



Function SiLatticeConst(TempC)				// computes temperature dependent Si Lattice constant
	Variable TempC							// temperature in degrees C
	Variable alpha=2.56E-6				// coefficient of expansion (1/�K)
	Variable dT = TempC - 22.5
	return aSi_A*(1+alpha*dT)
End

Function LookUpPhysicalConstant(name,[c,web,printIt])
	String name
	STRUCT PhysicalConstantStructure &c
	Variable printIt
	Variable web							// 1 try to use web first, 0 only use local
	c.valid = 0
	if (ParamIsDefault(printIt))
		printIt = strlen(GetRTStackInfo(2)) ? printIt : 1
	endif

	web = ParamIsDefault(web) ? 1 : web
	if (strlen(name)<1)
		return NaN
	endif

	String buf=""
	if (web)
		String sValue = FetchURL(PhysicalConstantServer)
		String errMsg = GetRTErrMessage()
		if (GetRTError(1)==0)
			Variable i1 = char2num(sValue[0])==char2num("\"") ? 1 : 0
			Variable i2=strlen(sValue)-1
			i2 = char2num(sValue[i2])==char2num("\"") ? i2-1 : i2
			web = strsearch(sValue,"Fundamental Physical Constants",0)>=0
			buf = SelectString(web,"",sValue[i1,i2])		// web version
		else
			printf "Could not get information from web, '%s'\r",errMsg
			web = 0
			buf = ""
		endif
	endif
	if (!web)												// not using web, use local version
		Variable f=0
		Open/R/Z=1 f as ParseFilePath(1,FunctionPath("LookUpPhysicalConstant"),":",1,0)+"Physical Constants.txt"
		if (V_flag)
			DoAlert 0,"Could not find file 'Physical Constants.txt'"
			return NaN
		endif
		FStatus f
		buf=PadString("",V_logEOF,0)
		FBinRead f, buf
		Close f
	endif
	if (strsearch(buf,"Fundamental Physical Constants",0)<0)
		return NaN
	endif

	buf = ReplaceString("\r\n",buf,"\n")
	buf = ReplaceString("\n\r",buf,"\n")
	buf = ReplaceString("\r",buf,"\n")
	Variable i=strsearch(buf,"----------------------------------",0)
	if (i<0)
		return NaN
	endif
	i = strsearch(buf,"\n",i+1)
	if (i<0)
		return NaN
	endif
	buf = buf[i+1,Inf]

	name = chooseConstant(buf,name)
	if (strlen(name)<1)
		return NaN
	endif

	STRUCT PhysicalConstantStructure clocal
	if (PhysicalConstantFromBuf(buf,name,clocal))
		return NaN
	endif

	if (!ParamIsDefault(c))
		copyPhysicalConstantStructure(c,clocal)
	endif
//	Variable/G PhysicalConstant_value=clocal.value, PhysicalConstant_err=clocal.err
//	String/G PhysicalConstant_unit=clocal.unit, PhysicalConstant_name=clocal.name
	if (printIt)
		print formatPhysicalConstantStructure(clocal)+SelectString(web,"   from local file","   from web server")
	endif
	return clocal.value
End

Static Function/T chooseConstant(buf,name)
	String buf
	String name

	String line,allNamesList="",choice=""
	Variable i1=0,i2
	i2 = strsearch(buf,"\n",i1)
	do
		line = TrimFrontBackWhiteSpace(buf[i1,i1+54])
		if (stringmatch(line,name))
			if (strlen(line))
				allNamesList += line+";"
			endif
		endif
		i1 = i2+1
		i2 = strsearch(buf,"\n",i1)
	while(i2>0)

	Variable N=ItemsInList(allNamesList)
	if (N<1)
		choice = ""
	elseif (N==1)
		choice = StringFromList(0,allNamesList)
	else
		Prompt choice,"Physical Constant Name",popup,"_search_;"+allNamesList
		DoPrompt "Physical Constant",choice
		if (V_flag)
			choice = ""
		endif
	endif
	if (stringmatch(choice,"_search_"))
		Prompt name,"enter search string, can use * and !"
		DoPrompt "Search",name
		if (V_flag)
			return ""
		endif
		choice = chooseConstant(buf,name)
	endif
	return choice
End
//
Static Function PhysicalConstantFromBuf(buf,name,c)
	// This routine is really stupid, but that is because the text file that I download from NIST is really stupid.  It has no rules,
	//	and there is little about the file that is standard.  If you can find either a more standard text file or an xml file that would be better.
	String buf
	String name
	STRUCT PhysicalConstantStructure &c

	c.valid = 0
	String strVal, strErr, unit=""
	Variable i1=0,i2
	i2 = strsearch(buf,"\n",i1)
	do
		if (stringmatch(TrimFrontBackWhiteSpace(buf[i1,i1+54]),name))
			String ss = buf[i1+55,i2-1]						// from start of value to end of line
			Variable isexp = stringmatch(ss[0,44],"*e*")		// the number has an 'e'. check all of numbers but none of units
			if (isexp)
				Variable j
				j = strsearch(ss,"e",0)
				j = strsearch(ss," ",j)
				strVal = ReplaceString(" ",ss[0,j],"")
				strVal = ReplaceString("...",strVal,"")			// sometimes used with "exact" constants
				c.value = str2num(strVal)
				strErr = ss[j,Inf]
				j = strsearch(strErr,"e",0)
				j = strsearch(strErr," ",j)
				unit = strErr[j+1,Inf]
				strErr = ReplaceString(" ",strErr[0,j],"")
				unit = TrimFrontBackWhiteSpace(unit)
				c.unit =  unit[0,PhysicalConstantMaxStr]
			else
				strVal = ReplaceString(" ",ss[0,22],"")
				strVal = ReplaceString("...",strVal,"")			// sometimes used with "exact" constants
				c.value = str2num(strVal)
				strErr = ReplaceString(" ",ss[23,43],"")
				c.err = str2num(strErr)
				c.exact = stringmatch(strErr,"*exact*")
				unit = TrimFrontBackWhiteSpace(ss[44,Inf])
				c.unit =  unit[0,PhysicalConstantMaxStr]
			endif
			name = TrimFrontBackWhiteSpace(buf[i1,i1+54])
			c.name = name[0,PhysicalConstantMaxStr]
			c.valid = 1
			break
		endif
		i1 = i2+1
		i2 = strsearch(buf,"\n",i1)
	while(i2>0)
	return !(c.valid)
End

Static Function/T formatPhysicalConstantStructure(c)
	STRUCT PhysicalConstantStructure &c

	if (!(c.valid))
		return ""
	endif
	String full, str, fmt
	fmt = "%."+num2str(placesOfPrecision(c.value))+"g"
	sprintf str,fmt,c.value
	full = c.name+" = "+str
	fmt = "%."+num2str(placesOfPrecision(c.err))+"g"
	sprintf str,fmt,c.err
	String pm = SelectString(stringmatch(IgorInfo(2),"Macintosh"),"+/-","�")
	full += SelectString(numtype(c.err)," � "+str,"")
	full += SelectString(strlen(c.unit),""," ("+c.unit+")")
	full += SelectString(c.exact,""," (exact)")
	return full
End



// creates (or overwrites) a text file containing all of the constants from PhysicalConstantServer
// The text file is named "Physical Constants.txt" and will be in the same folder as this ipf file
Function UpdateLocalCopyOfConstants()
	String sValue = FetchURL(PhysicalConstantServer)
	String errMsg = GetRTErrMessage()
	if (GetRTError(1))
		printf "Could not get information from web, '%s'\r",errMsg
		return 1
	endif

	Variable i1 = char2num(sValue[0])==char2num("\"") ? 1 : 0
	Variable i2=strlen(sValue)-1
	i2 = char2num(sValue[i2])==char2num("\"") ? i2-1 : i2
	if (i2<1000)					// too short to be right
		return 1
	endif
	String buf=sValue[i1,i2]		// buf now has the results to write
	DoAlert 1,"Update your local file of Physical Constants with values from Web?"
	if (V_flag!=1)
		return 1
	endif

	String fileName = ParseFilePath(1,FunctionPath("UpdateLocalCopyOfConstants"),":",1,0)+"Physical Constants.txt"
	Variable f=0
	Open/Z=1 f as fileName
	if (V_flag)
		DoAlert 0,"Could not find file 'Physical Constants.txt'"
		return NaN
	endif
	FBinWrite f, buf
	Close f
	printf "Your local copy of 'Physical Constants.txt' has been updated from '%s'\r",PhysicalConstantServer
End

Static Function/T TrimFrontBackWhiteSpace(str)
	String str
	str = TrimLeadingWhiteSpace(str)
	str = TrimTrailingWhiteSpace(str)
	return str
End
//
Static Function/T TrimLeadingWhiteSpace(str)
	String str
	Variable i, N=strlen(str)
	for (i=0;char2num(str[i])<=32 && i<N;i+=1)	// find first non-white space
	endfor
	return str[i,Inf]
End
//
Static Function/T 	TrimTrailingWhiteSpace(str)
	String str
	Variable i
	for (i=strlen(str)-1; char2num(str[i])<=32 && i>=0; i-=1)	// find last non-white space
	endfor
	return str[0,i]
End
//
Static Function placesOfPrecision(a)		// number of significant figures in a number (at most 16)
	Variable a
	a = roundSignificant(abs(a),17)
	Variable i
	for (i=1;i<18;i+=1)
		if (abs(a-roundSignificant(a,i))/a<1e-15)
			break
		endif
	endfor
	return i
End
//
// This routine is much faster than going through an [sprintf str,"%g",val] conversion
Static Function roundSignificant(val,N)	// round val to N significant figures
	Variable val			// input value to round
	Variable N			// number of significant figures

	if (val==0 || numtype(val))
		return val
	endif
	Variable is,tens
	is = sign(val) 
	val = abs(val)
	tens = 10^(N-floor(log(val))-1)
	return is*round(val*tens)/tens
End
