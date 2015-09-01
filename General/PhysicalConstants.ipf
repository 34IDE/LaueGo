#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName=PhysicalConstants
#pragma version = 2.11
#pragma IgorVersion = 6.3
#include "Utility_JZT", version>=3.51		// supplies:  TrimFrontBackWhiteSpace(str), TrimLeadingWhiteSpace(str), TrimTrailingWhiteSpace(str), placesOfPrecision(a), roundSignificant(val,N)
Static StrConstant PhysicalConstantServerURL="http://physics.nist.gov/cuu/Constants/Table/allascii.txt"

//	By Jon Tischler (ORNL)  Aug 12, 2010
//
// functions of interest to users are:
//		GetPhysicalConstant(name,[c,printIt])
//		PhysicalConstant_InsertStatic(name,[printIt])
//
//			rarely someone will want to use:
//		DateOfLocalPhysicalConstants([printIt])
//		UpdateLocalCopyOfConstants()
//		SiLatticeConst(TempC)
//
//			very rarely someone will want to use:
//		formatPhysicalConstantStructure(c,[always])


//	This list of constants is not really needed, it is just a convienence for copying and pasting into other experiments.
//	Set to 2010 CODATA values, Jan 7, 2015,   from  http://physics.nist.gov/cuu/Constants
//
//	Constant c_ms		= 299792458				// speed of light (m/s) (exact)
//	Constant eo_Fm		= 8.85418781762039e-12		// permitivity of vacuum (F/m),  = 1/(4*PI*1e-7)/(299792458)^2 (exact)
//	Constant G_MKS		= 6.67384e-11			// Gravitational constant (m^3 1/kg 1/s)
//	Constant h_Js		= 6.62606957e-34		// Planck constant (J s)
//	Constant h_eVs		= 4.135667516e-15		// Plank constant (eV s)
//	Constant e_C		= 1.602176565e-19		// Charge on electron (C)
//	Constant alpha_fine = 7.2973525698e-3	// fine structure constant
//	Constant aB_m		= 0.52917721092e-10	// Bohr radius (m)
//	Constant Rydberg	= 13.60569253			// Rydberg (eV) ( = alpha^2*me*c/2/h * hc)
//	Constant me_kg		= 9.10938291e-31		// mass of electron (kg)
//	Constant me_keV	= 510.998928			// mass of electron (keV)
//	Constant me_mp		= 5.4461702178e-4		// mass ratio m(electron)/m(proton)
//	Constant mp_kg		= 1.672621777e-27		// mass of proton (kg)
//	Constant mn_kg		= 1.674927351e-27		// mass of neutron (kg)
//	Constant mn_eV		= 9.39565379e8			// mass of neutron (eV)
//	Constant NA			= 6.02214129e+23		// Avagadro's number
//	Constant kB_JK		= 1.3806488e-23		// Boltzman constant (J/K)
//	Constant kB_eVK	= 8.6173324e-5			// Boltzman constant (eV/K)
//	Constant J_eV		= 1.602176565e-19		// Joules/eV
//	Constant atm_Pa	= 101325					// Number of Pa in one atm (exact)
//	Constant hc_keVA	= 12.3984193			// h*c (keV-Å)
//	Constant hc_keVnm = 1.23984193			// h*c (keV-nm)
//	Constant re_m		= 2.8179403267e-15	// Thompson radius (m)
//	Constant re_A		= 2.8179403267e-5		// Thompson radius (Å)
//	Constant CuKa1		= 1.540593226			// wavelength of Cu Kalpha(1) (Å)
//	Constant MoKa1		= 0.709317155			// wavelength of Mo Kalpha(1) (Å)
//	Constant aSi_A		= 5.4310205052			// lattice constant of Si at 22.5° (Å), 2010 CODATA,  alpha=2.56E-6 (deg/K)


Menu "Analysis"
	SubMenu "Physical Constants"
		"New Static Constant, Physical Constant string",PhysicalConstant_InsertStatic("*")
		"<BGet a Physical Constant...",GetPhysicalConstant("*")
		"  Browse NIST web site...", PhysicalConstants#BrowseConstantWebSite()
		"<I  update your local copy [Rarely needed]",PhysicalConstants#UpdateLocalCopyOfConstants()
		"  (constants updated:  "+PhysicalConstants#DateOfLastUpdate()
		//	"  date of your local copy",PhysicalConstants#DateOfLocalPhysicalConstants()
	End
End
Static Function/T DateOfLastUpdate()	// provides dynamic text for menu item
	Variable epoch = DateOfLocalPhysicalConstants()
	return Secs2Date(epoch,1)+",  "+Secs2Time(epoch,0)
End


// Browse the NIST cPhysical Constants web site"
Static Function BrowseConstantWebSite()
	String url = ReplaceString("Table/allascii.txt",PhysicalConstantServerURL,"")
	BrowseURL url
End


// writes to history a string that you can copy/paste into your procedure.
Function PhysicalConstant_InsertStatic(name,[printIt])
	String name
	Variable printIt
	if (ParamIsDefault(printIt) || numtype(printIt))
		printIt = strlen(GetRTStackInfo(2))==0
	endif
	printIt = !(!printIt)
	if (strlen(name)<1)
		return NaN
	endif
	STRUCT PhysicalConstantStructure c
	GetPhysicalConstant(name,c=c,printIt=0)
	if (!(c.valid))
		return NaN
	endif
	name = c.name
	Variable value = c.value

	// other names (and scalings) for commonly used constants
	String conversions="inverse meter-electron volt relationship=hc_keVnm,1e6,h*c (keV-nm):hc_keVA,1e7,h*c (keV-Å);"
	conversions += "classical electron radius=re_m,1,Thompson radius (m):re_nm,1e9,Thompson radius (nm):re_A,1e10,Thompson radius (Å);"
	conversions += "speed of light in vacuum=c_ms,1,speed of light in vacuum (m s^-1):c_nms,1e9,speed of light in vacuum (nm s^-1);"
	conversions += "{220} lattice spacing of silicon=aoSi022_m,1,Si {220} (m):aoSi022_nm,1e9,Si {220} (nm):aoSi022_A,1e10,Si {220} (Å);"
	conversions += "Avogadro constant=NA,1,Avogadro constant (1 mole);"
	conversions += "electron mass energy equivalent in MeV="
		conversions += "me_eV,1e6,electron mass energy equivalent (eV):"
		conversions += "me_keV,1e3,electron mass energy equivalent (keV):"
		conversions += "me_MeV,1,electron mass energy equivalent (MeV);"

	String lists=StringByKey(name,conversions,"="), outName, comment
	Variable N=ItemsInList(lists,":"),i
	if (N<1)
		outName = CleanupName(name,0)
		comment = name
		if (strlen(c.unit))
			String sunit = CleanupName(c.unit,0)
			i = strlen(sunit)
			outName = outName[0,30-i-1]
			for (i=strlen(outName)-1; stringmatch(outName[i],"_") && i>=0; i-=1)	// find last non-underscore
			endfor
			outName = outName[0,i]+"_"+sunit
			comment += " ("+c.unit+")"
		endif
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
		if (strlen(c.unit))
			printf "from \"%s\"  =  %.15g  (%s)\r",c.name,c.value,c.unit
		else
			printf "from \"%s\"  =  %.15g\r",c.name,c.value
		endif
		printf "insert:\t\tStatic Constant %s = %.15g%s\r",outName,value,comment
		return value
//	else
//		String cmd
//		sprintf cmd, "INSERTINCLUDE \"Static Constant %s = %.15g%s\"", outName,value,comment
//		DoAlert/T="Insert Constant\r" 1, "Insert line\r"+cmd
//		V_flag = 1
//		if (V_flag==1)
//			printf "from \"%s\"  =  %.13g  (%s)\r",c.name,c.value,c.unit
//			printf "insert:\t\tStatic Constant %s = %.13g%s\r",outName,value,comment
//			Execute/P cmd
//		else
//			print "nothing done"
//		endif
	endif
	return 0
End



//Function test()
//	GetPhysicalConstant("speed of light in vacuum",printIt=1)
//	GetPhysicalConstant("speed of light*",printIt=1)
//	GetPhysicalConstant("mag. constant",printIt=1)
//	GetPhysicalConstant("molar mass constant",printIt=1)
//	GetPhysicalConstant("electron-muon mass ratio",printIt=1)
//	GetPhysicalConstant("{220}*",printIt=1)
//	GetPhysicalConstant("*220*",printIt=1)
//	GetPhysicalConstant("Wien wavelength*",printIt=1)
//	GetPhysicalConstant("Wien wavelength displacement law constant",printIt=1)
//	GetPhysicalConstant("",printIt=1)
//	GetPhysicalConstant("error",printIt=1)
//	GetPhysicalConstant("\n",printIt=1)
//	GetPhysicalConstant("hc",printIt=1)
//End
//
Function GetPhysicalConstant(name,[c,printIt])	// returns value of constant
	String name
	STRUCT PhysicalConstantStructure &c	// structure that is filled in this routine
	Variable printIt
	if (!ParamIsDefault(c))
		initPhysicalConstantStructure(c)	// mainly set c.valid=0
		c.name = name[0,PhysicalConstantMaxStrLen]
	endif
	if (ParamIsDefault(printIt) || numtype(printIt))
		printIt = strlen(GetRTStackInfo(2))==0
	endif

	if (strlen(name)<1)
		if (printIt)
			print "*** no name given"
		endif
		return NaN
	endif

	STRUCT PhysicalConstantStructureAll cAll
	LoadPackagePreferences/MIS=1 "PhysicalConstantsJZT" , "PhysicalConstantsPrefs", 0, cAll
	if (V_bytesRead != V_structSize || V_flag)
		DoAlert/T="Physical Constants" 1, "Update Physical Constants from NIST web site?"
		if (V_flag==1)									// could not load from Preferences, goto web
			UpdateLocalCopyOfConstants()			// get constants from NIST web site, and try again
			LoadPackagePreferences "PhysicalConstantsJZT" , "PhysicalConstantsPrefs", 0, cAll
		else
			return NaN
		endif
	endif
	Variable ic = chooseConstant(cAll,name)	// ask user which constant?
	STRUCT PhysicalConstantStructure clocal	// store result into clocal
	getStruct_i(cAll,ic,clocal)
	name = clocal.name

	if (strlen(name)<1 || !(clocal.valid))
		if (printIt)
			print "*** name = '"+name+"',  is INVALID"
		endif
		return NaN
	endif

	if (!ParamIsDefault(c))
		copyPhysicalConstantStructure(c,clocal)	// struct c was supplied, so fill it
	endif
	if (printIt)
		print formatPhysicalConstantStructure(clocal,always=1)
	endif
	return clocal.value
End
//
Static Function chooseConstant(cAll,name)
	STRUCT PhysicalConstantStructureAll &cAll
	String name

	if (strsearch(name,"hc",0)==0 || strsearch(name,"h*c",0)==0)
		name = "inverse meter-electron volt relationship"		// this is the "official" name of "hc"
	elseif (StringMatch(name,"c"))
		name = "speed of light in vacuum"								// this is the "official" name of "c"
	endif
	String line, allNamesList="",choice=""

	Variable N = cAll.N1 + cAll.N2 + cAll.N3 + cAll.N4
	STRUCT PhysicalConstantStructure clocal
	Variable i, ic=-1
	for (i=0;i<N;i+=1)
		getStruct_i(cAll,i,clocal)
		if (stringmatch(clocal.name,name) && strlen(clocal.name))
			if (clocal.valid)
				allNamesList += clocal.name+";"
			endif
		endif
	endfor
	Variable Nlist = ItemsInList(allNamesList)
	if (Nlist<1)
		choice = ""
	elseif (Nlist==1)
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
			return -1
		endif
		ic = chooseConstant(cAll,name)
	endif

	for (i=0;i<N;i+=1)
		getStruct_i(cAll,i,clocal)
		if (stringmatch(clocal.name,choice) && strlen(clocal.name))
			ic = i
			break
		endif
	endfor
	return ic
End
//
// returns a nice string showing contents of c
ThreadSafe Static Function/T formatPhysicalConstantStructure(c,[always])
	STRUCT PhysicalConstantStructure &c
	Variable always
	always = ParamIsDefault(always) || numtype(always) ? 0 : !(!always)

	String full="", str, fmt
	if (!(c.valid))
		if (always)
			sprintf full, "name = '%s',  is INVALID\r",c.name
		endif
		return full
	endif
	fmt = "%."+num2str(placesOfPrecision(c.value))+"g"
	sprintf str,fmt,c.value
	full = c.name+" = "+str
	fmt = "%."+num2str(placesOfPrecision(c.err))+"g"
	sprintf str,fmt,c.err
	String pm = SelectString(stringmatch(IgorInfo(2),"Macintosh"),"+/-","±")
	full += SelectString(numtype(c.err)," ± "+str,"")
	full += SelectString(strlen(c.unit),""," ("+c.unit+")")
	full += SelectString(c.exact,""," (exact)")
	return full
End


// returns epoch of current copy of PhysicalConstants, and optionally prints too
Static Function DateOfLocalPhysicalConstants([printIt])
	Variable printIt
	printIt = (ParamIsDefault(printIt) || numtype(printIt)) ? strlen(GetRTStackInfo(2))==0 : printIt

	STRUCT PhysicalConstantStructureAll cAll
	LoadPackagePreferences/MIS=1 "PhysicalConstantsJZT" , "PhysicalConstantsPrefs", 0, cAll
	if (V_bytesRead != V_structSize || V_flag)
		print "ERROR -- no current copy of Physical Constants found"
		return NaN
	endif
	if (printIt)
		printf "Local copy of constants last updated from NIST web server on:  %s,  %s\r",Secs2Date(cAll.updateEpoch,1), Secs2Time(cAll.updateEpoch,1)
	endif
	return cAll.updateEpoch
End


// Creates (or overwrites) the Igor Package Prefs file containing all of the constants from PhysicalConstantServerURL
Static Function UpdateLocalCopyOfConstants()
	String buf = getFullASCIIfromWeb()		// ascii buffer with all of the constants info, must have a terminating "\n"
	if (strlen(buf)<200)
		print buf
		return 1
	endif

	STRUCT PhysicalConstantStructureAll cAll
	Variable N = FillConstantStucturesFromBuf(buf, cAll)
	SavePackagePreferences/FLSH=1 "PhysicalConstantsJZT" , "PhysicalConstantsPrefs", 0 , cAll
	printf "Your local copy of 'Physical Constants' has been updated from '%s'\r",PhysicalConstantServerURL
//	String fileName = ParseFilePath(1,FunctionPath("UpdateLocalCopyOfConstants"),":",1,0)+"Physical Constants.txt"
End
//
// Return the ascii buffer with all of the constants info from web (the NIST web server)
Static Function/T getFullASCIIfromWeb()
	String buf=""
	String sValue = FetchURL(PhysicalConstantServerURL)
	String errMsg = GetRTErrMessage()
	if (GetRTError(1))
		printf "ERROR -- Could not get information from web, '%s'\r",errMsg
		return ""
	endif

	Variable i1 = char2num(sValue[0])==char2num("\"") ? 1 : 0
	Variable i2=strlen(sValue)-1
	i2 = char2num(sValue[i2])==char2num("\"") ? i2-1 : i2
	if (strsearch(sValue,"Fundamental Physical Constants",0)<0)
		printf "Could not get information from web, '%s'\r",errMsg
		return ""
	endif
	buf = sValue[i1,i2]
	if (strsearch(buf,"Fundamental Physical Constants",0)<0)
		return "ERROR -- Downloaded 'allascii.txt' file is INVALID"
	endif

	buf = ReplaceString("\r\n",buf,"\n")
	buf = ReplaceString("\n\r",buf,"\n")
	buf = ReplaceString("\r",buf,"\n")
	Variable i=strsearch(buf,"----------------------------------",0)
	if (i<0)
		return "ERROR -- 'allascii.txt' file is INVALID"
	endif
	i = strsearch(buf,"\n",i+1)
	if (i<0)
		return "ERROR -- 'allascii.txt' file is INVALID"
	endif
	buf = TrimFrontBackWhiteSpace(buf[i+1,Inf])
	buf += "\n"				// ensure terminating <NL>
	return buf
End
//
// take ASCII result from web server and fill cAll
Static Function FillConstantStucturesFromBuf(buf,cAll)
	String buf
	STRUCT PhysicalConstantStructureAll &cAll

	Variable val0=NaN, err0=NaN, unit0=NaN
	// This routine is really stupid, but that is because the text file that I download from NIST is really stupid.  It has no rules,
	//	and there is little about the file that is standard.  If you can find either a more standard text file or an xml file that would be better.
	Variable i = strsearch(buf,"\nmolar mass constant ",0)
	if (i<0)
		return 0
	endif
	String line=buf[i+1,i+300]
	i = strsearch(line,"\n",0)-1
	if (i<0)
		return 0
	endif
	line = line[0,i]

	// assuming fixed format lengths for name, value, error, & units
	val0 = nextNonSpace(line,strlen("molar mass constant "))
	err0 = nextNonSpace(line,val0+6)
	unit0 = nextNonSpace(line,err0+8)
	if (numtype(val0+err0+unit0))
		return 0
	endif

	STRUCT PhysicalConstantStructure clocal
	String strVal, strErr, unit, name
	Variable N=ItemsInList(buf,"\n")
	Variable nConstants=0, j,m
	cAll.N1 = 0
	cAll.N2 = 0
	cAll.N3 = 0
	cAll.N4 = 0
	for (i=0;i<100;i+=1)
		initPhysicalConstantStructure(call.c1[i])
		initPhysicalConstantStructure(call.c2[i])
		initPhysicalConstantStructure(call.c3[i])
		initPhysicalConstantStructure(call.c4[i])
	endfor

	for (i=0,line="xxx"; i<N && strlen(line); i+=1)
		line = TrimFrontBackWhiteSpace(StringFromList(i,buf,"\n"))
		if (strlen(line)<2)
			continue
		endif

		strVal = line[val0,err0-1]
		strVal = ReplaceString("...",strVal,"")	// sometimes used with "exact" constants
		strVal = ReplaceString(" ",strVal,"")		// no spaces allowed in value
		clocal.value = str2num(strVal)

		strErr = line[err0,unit0-1]
		clocal.exact = stringmatch(strErr,"*(exact)*")
		strErr = ReplaceString(" ",strErr,"")
		clocal.err = str2num(strErr)

		unit = TrimFrontBackWhiteSpace(line[unit0,Inf])
		clocal.unit = unit[0,PhysicalConstantMaxStrLen]

		name = TrimFrontBackWhiteSpace(line[0,val0-1])
		name = ReplaceString("mom.",name,"moment")	// remove abbreviations
		name = ReplaceString("mag.",name,"magnetic")
		clocal.name = name[0,PhysicalConstantMaxStrLen]
		clocal.valid = 1									// set valid to true

		j = mod(nConstants,100)
		m = floor(nConstants/100)
		nConstants += 1
		if (m==0)
			cAll.N1 += 1
			copyPhysicalConstantStructure(cAll.c1[j],clocal)
		elseif (m==1)
			cAll.N2 += 1
			copyPhysicalConstantStructure(cAll.c2[j],clocal)
		elseif (m==2)
			cAll.N3 += 1
			copyPhysicalConstantStructure(cAll.c3[j],clocal)
		elseif (m==3)
			cAll.N4 += 1
			copyPhysicalConstantStructure(cAll.c4[j],clocal)
		else
			break
		endif
	endfor
	cAll.updateEpoch = DateTime
	return nConstants
End
//
Static Function nextNonSpace(str,start)
	String str
	Variable start
	Variable i,N=strlen(str)
	for (i=start;i<N;i+=1)
		if (char2num(str[i])>32)
			return i
		endif
	endfor
	return i
End


// Returns Si Lattice constant (Å) at a particular temperature
ThreadSafe Function SiLatticeConst(TempC)
	Variable TempC							// temperature in degrees C
	Variable alpha=2.56E-6				// coefficient of expansion (1/°K)
	Variable dT = TempC - 22.5
	return 5.4310205052*(1+alpha*dT)	// 5.4310205052 = aoSi at 22.5° (Å), 2010 CODATA
End



// PhysicalConstant structures
Static Constant PhysicalConstantMaxStrLen=60
Static Structure PhysicalConstantStructure
	int16	valid
	char name[PhysicalConstantMaxStrLen+1]
	double value
	double err
	char unit[PhysicalConstantMaxStrLen+1]
	int16 exact
EndStructure
//
Static Structure PhysicalConstantStructureAll
	double updateEpoch				// Igor epoch when last updated
	int16	N1
	STRUCT PhysicalConstantStructure c1[100]
	int16	N2
	STRUCT PhysicalConstantStructure c2[100]
	int16	N3
	STRUCT PhysicalConstantStructure c3[100]
	int16	N4
	STRUCT PhysicalConstantStructure c4[100]			// c4 only needs 35 of these 100
EndStructure
//
ThreadSafe Static Function copyPhysicalConstantStructure(f,i)
	STRUCT PhysicalConstantStructure &f, &i
	f.valid	= i.valid
	f.name	= i.name
	f.value	= i.value
	f.err		= i.err
	f.unit	= i.unit
	f.exact	= i.exact
End
//
ThreadSafe Static Function initPhysicalConstantStructure(c)
	STRUCT PhysicalConstantStructure &c
	c.valid	= 0			// init to NOT valid
	c.name	= ""
	c.value	= NaN
	c.err		= NaN
	c.unit	= ""
	c.exact	= 0
End
//
// Copy constant structure i from cAll into ci
ThreadSafe Static Function getStruct_i(cAll,i,ci)
	STRUCT PhysicalConstantStructureAll &cAll
	Variable i
	STRUCT PhysicalConstantStructure &ci

	Variable j = mod(i,100)		// index into c1, c2, c3, or c4
	Variable m = floor(i/100)		// which group of 100 to choose
	if (m==0)
		copyPhysicalConstantStructure(ci,cAll.c1[j])
	elseif (m==1)
		copyPhysicalConstantStructure(ci,cAll.c2[j])
	elseif (m==2)
		copyPhysicalConstantStructure(ci,cAll.c3[j])
	elseif (m==3)
		copyPhysicalConstantStructure(ci,cAll.c4[j])
	else
		initPhysicalConstantStructure(ci)
	endif
End
