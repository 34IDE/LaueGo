#pragma rtGlobals= 2
#pragma version = 3.03
#pragma ModuleName = JZTgeneral
#pragma hide = 1
//	#pragma IndependentModule=JZTgeneral
//	DefaultFont "Helvetica"		// This is in "JonFirst.ipf", that is enough

Menu "Analysis"
	Submenu "Packages"
		"-"
		"Physical Constants",Execute/P "INSERTINCLUDE  \"PhysicalConstants\"";Execute/P "COMPILEPROCEDURES "
		help = {"Load Package of Physical Constants."}
		"FWHM",Execute/P "INSERTINCLUDE  \"FWHM\"";Execute/P "COMPILEPROCEDURES "
		help = {"Load commands for getting FWHM."}
		Submenu "spec"
			"spec with images",Execute/P "INSERTINCLUDE  \"specImages\", version>=0.44";Execute/P "COMPILEPROCEDURES ";Execute/P "init_specImage(\"\")"
			help = {"Support for spec files plus support for spec with image files (like a Pilatus)."}
			"spec files",Execute/P "INSERTINCLUDE  \"spec\", version>=2.27";Execute/P "COMPILEPROCEDURES ";Execute/P "specInitPackage()"
			help = {"Load everything for reading spec files."}
		End
		"MDA files",Execute/P "INSERTINCLUDE  \"mdaFiles\"";Execute/P "COMPILEPROCEDURES "
		help = {"Load mda files (from APS)."}
		"BURT Files",Execute/P "INSERTINCLUDE  \"BurtFiles\"";Execute/P "COMPILEPROCEDURES "
		help = {"Load Burt files (from APS)."}
		"Gizmo Markers & Zoom", Execute/P "INSERTINCLUDE \"GizmoZoomTranslate\", version>=1.36" ; 	Execute/P "COMPILEPROCEDURES ";Execute/P "InitGizmoZoomTranslate()"
		help = {"Provide support of Zooming Gizmos and putting a 'Cursor' on a Gizmo"}
//		"Color Names", Execute/P "INSERTINCLUDE \"ColorNames\"" ; 	Execute/P "COMPILEPROCEDURES "
//		help = {"Load procedures that will provide names for RGBA colors"}
		"Image Display Range", Execute/P "INSERTINCLUDE \"ImageDisplayScaling\"" ; 	Execute/P "COMPILEPROCEDURES "
		help = {"Load procedures that set range of color table based on a Marquee region"}

		"WinView Reader", Execute/P "INSERTINCLUDE \"WinView\", version>=2.01"
		help = {"Load procedures reading and looking at WinView images"}
		SubMenu "X-ray"
			"X-ray Data",Execute/P "INSERTINCLUDE  \"Xray\", version>=2.21";Execute/P "COMPILEPROCEDURES ";Execute/P "ElementDataInitPackage()"
			help = {"Load procedures for providing X-ray data"}
			"  Elements Data",Execute/P "INSERTINCLUDE  \"Elements\", version>=1.69";Execute/P "COMPILEPROCEDURES ";Execute/P "ElementDataInitPackage()"
			help = {"Load procedures for providing data on the elements"}
			"  Cromer-Liberman",Execute/P "INSERTINCLUDE  \"CromerLiberman\", version>=1.7";Execute/P "COMPILEPROCEDURES "
			help = {"Load procedures for providing Cromer-Liberman"}
			"Ion Chamber",Execute/P "INSERTINCLUDE  \"IonChamber\", version>=3.2";Execute/P "COMPILEPROCEDURES ";Execute/P "ionChamberInitPackage()"
			help = {"Load procedures for evaluating ion chamber output"}
			"Lattices",Execute/P "INSERTINCLUDE  \"LatticeSym\", version>=3.77";Execute/P "COMPILEPROCEDURES ";Execute/P "InitLatticeSymPackage()"
			help = {"Load lattice symmetry procedures"}
		End
		"Neutron Scattering Data",Execute/P "INSERTINCLUDE  \"Neutron\", version>=1.1";Execute/P "COMPILEPROCEDURES ";Execute/P "NeutronDataInitPackage()"
		help = {"Load procedures for providing Neutron Scattering Data"}
		SubMenu "Add Local User Packages"
			"(Moved to:  File -> Add Local User Packages"
		End
		"-"
	End
End
//
Menu "File"
	"-"
	SubMenu "Add Local User Package -->"
		JZTgeneral_CheckForUserPackages(""), /Q,JZTgeneral#JZTgeneral_StartUpLocalPackage()
	End
End
//
Menu "Data"
	SubMenu "Packages"
		"MDA files from APS",Execute/P "INSERTINCLUDE  \"mdaFiles\"";Execute/P "COMPILEPROCEDURES "
		help = {"Load mda files (from APS)."}
		"BURT files from APS",Execute/P "INSERTINCLUDE  \"BurtFiles\"";Execute/P "COMPILEPROCEDURES "
		help = {"Load Burt files (from APS)."}
	End
End
Menu "Load Waves"
	SubMenu "Packages"
		"MDA files from APS",Execute/P "INSERTINCLUDE  \"mdaFiles\"";Execute/P "COMPILEPROCEDURES "
		help = {"Load mda files (from APS)."}
		"BURT files from APS",Execute/P "INSERTINCLUDE  \"BurtFiles\"";Execute/P "COMPILEPROCEDURES "
		help = {"Load Burt files (from APS)."}
	End
End
//
Menu "Misc"
	"StopAllTimers",StopAllTimers()
End




//	This provides for an easy way to load any package in the folder "Documents:WaveMetrics:Igor Pro 6 User Files:User Procedures/Local Packages"
//	just put the ipf file in there and it will be available to this function for loading.
//	if the first 50 lines of the file contain a line such as   "#requiredPackages "HDF5images;microGeometryN;", then it will only load that file if the 
//	ipf files in the list are already loaded.
//	if the package contains the line "#excludeFromPackageMenu", then it will not show up in the list of packages to load (useful for subroutines)
//	Note if the ipf file contains both #requiredPackages and #excludeFromPackageMenu, then the behavior depends upon order, but that situation is stupid.
Function/T JZTgeneral_CheckForUserPackages(dirPath)
	String dirPath				// full path to a directory

	String rootPath = SpecialDirPath("Documents", 0, 0, 0 )+"WaveMetrics:Igor Pro 6 User Files:User Procedures:Local Packages:"
	dirPath = SelectString(strlen(dirPath),rootPath,dirPath)
	String pre = dirPath[strlen(rootPath),Inf]
	pre += SelectString(strlen(pre),"",":")
	String pathName = UniqueName("lpath",12,0)
	NewPath/O/Q/Z $pathName ,  dirPath
	if (V_flag)
		return ""
	endif
	String fname,menuList="",required
	Variable i=0
	do
		fname = IndexedFile($pathName,i,".ipf")
		if (strlen(fname)==0)
			break
		endif
		if (allProcsArePresent(ReplaceString(".ipf",fname,"")))
			i += 1													// skip if fname is already included
			continue
		endif
		required = getListOfRequiredProcs(pathName,fname)
		if (allProcsArePresent(required))							// check if required procs are already included
			menuList = AddListItem(pre+fname,menuList,";",Inf)
		endif
		i += 1
	while (1)

	//	check subdirectories, but exclude those labeled as 'always', 'old', or 'subroutine*'
	String dirName="", dirList="", addMenus
	for (i=0,dirName=IndexedDir($pathName,0,1); strlen(dirName); i+=1, dirName=IndexedDir($pathName,i,1))
		if (strsearch(dirName,"subroutine",0,2)>=0)
			continue
		elseif (strsearch(dirName,":old",0,2)>=0)
			continue
		elseif (strsearch(dirName," old",0,2)>=0)
			continue
		elseif (strsearch(dirName,"(old",0,2)>=0)
			continue
		elseif (StringMatch(dirName,"old"))
			continue
		elseif (StringMatch(dirName,"always"))
			continue
		endif
		addMenus = JZTgeneral_CheckForUserPackages(dirName)	// check inside this sub-directory recursively
		menuList += SelectString(strlen(addMenus),"","-;") + addMenus
	endfor
	KillPath/Z $pathName
	return menuList
End


Static Function/T getListOfRequiredProcs(path,ipf)
	String path, ipf
	String line, list=""
	Variable f, i,i0,i1
	Open/R/P=$path/T=".ipf"/Z f as ipf
	for (i=0;i<50;i+=1)							// check only first 50 lines
		FReadLine/N=500 f, line
		if (strlen(line)==0)						// end of file
			Close f
			return list
		endif
		if (strsearch(line,"#requiredPackages",0)==0)
			i0 = strsearch(line,"\"",0)
			i1 = strsearch(line,"\"",i0+1)
			if (i0==-1 || i1==-1 || i1<=i0)
				break
			endif
			list = line[i0+1,i1-1]
			break
		endif
		if (strsearch(line,"#excludeFromPackageMenu",0)==0)
			list = Hash("NON-EXISTANT",1)		// this "package" will not exist
			break
		endif
	endfor
	Close f
	list = ReplaceString(",",list,";")				// change commas to semi-colons
	return list
End
//
Static Function JZTgeneral_StartUpLocalPackage()
	GetLastUserMenuInfo							// sets S_value, V_value, etc.
	String pkg=S_value
	if (V_flag || !stringmatch(pkg,"*.ipf" ))
		return 1
	endif
	String initFuncs=getInitFunctionsName(pkg), cmd
	pkg = ReplaceString(".ipf",pkg,"")
	Variable i=strsearch(pkg,":",inf,1)+1				// start of name part, strip off any leading path part
	if (i>=strlen(pkg))
		DoAlert 0, "ERROR:\rCould not include package '"+pkg+"'"
		return 0
	endif
	pkg = pkg[i,inf]									// trim off any path part
	print "\r  adding local package  ",pkg
	sprintf cmd, "Execute/P \"INSERTINCLUDE  \\\"%s\\\"\";Execute/P \"COMPILEPROCEDURES \"", pkg
	if (strlen(initFuncs))
		cmd += ";Execute/P \""+initFuncs+" \""
	endif
	//	print "cmd = ",cmd
	Execute cmd
	return 0
End
//
Static Function allProcsArePresent(list)
	String list
	String ipf
	Variable i,N=ItemsInList(list)
	for (i=0;i<N;i+=1)
		ipf = StringFromList(i,list)+".ipf"
		if (ItemsInList(WinList(ipf,";","WIN:128"))<1)
			return 0
		endif
	endfor
	return 1
End
//
Static Function/T getInitFunctionsName(ipf)
	String ipf										// name of name of ipf file
	String pathName = UniqueName("lpath",12,0)
	NewPath/O/Q/Z $pathName ,  SpecialDirPath("Documents", 0, 0, 0 )+"WaveMetrics:Igor Pro 6 User Files:User Procedures:Local Packages:"
	GetFileFolderInfo/P=$pathName/Q/Z=1 ipf
	if (!V_isFile && strsearch(ipf,".ipf",0,2)<0)
		ipf += ".ipf"
	endif

	String line, initFunc=""
	Variable f, i,i0,i1
	Open/R/P=$pathName/T=".ipf"/Z f as ":"+ipf
	if (V_flag)
		DoAlert 0, "ERROR:\rCould not re-open file '"+ipf+"'"
		KillPath/Z $pathName
		return ""
	endif
	for (i=0;i<50;i+=1)
		FReadLine/N=500 f, line
		if (strlen(line)==0)						// end of file
			Close f
			return initFunc
		endif
		if (strsearch(line,"#initFunctionName",0,2)==0)
			i0 = strsearch(line,"\"",0)
			i1 = strsearch(line,"\"",i0+1)
			if (i0==-1 || i1==-1 || i1<=i0)
				break
			endif
			initFunc = line[i0+1,i1-1]
			break
		endif
	endfor
	Close f
	KillPath/Z $pathName
	return initFunc
End
