#pragma rtGlobals=1		// Use modern global access method.
#pragma version = 0.36
#pragma ModuleName=HDF5images

// Dec 12, 2009, version 0.200		Added support for multiple images in one HDF5 file
// Dec 12, 2009, version 0.208		Added support for wire being in its own folder
// Oct 8, 2012, version 0.213		Changed the way slices were dealt with
// Nov 15, 2012, version 0.23		Changed load commands so that dark and slice are passed as "slice:1;dark:fullNameofWave;" in extras 
// Nov 29, 2012, version 0.24		added HDFSaveImage()
//	Jun 19, 2014, version 0.30		added optional string extras to argument of ReadHDF5header()
//	Mar 19, 2015, version 0.31		added support in GuessDetectorID() for new Perkin-Elmer side detectors, but this routine is DEPRECATED
//												also added NewImageGraphHDF5proto(), in cases where ImageDisplayScaling.ipf has not been included
//	Apr  5, 2015, version 0.32		NewImageGraph() now has withButtons as an optional parameter
//	Mar  5, 2016, version 0.33		LoadHDF5imageFile(...), supports read of an roi, also, ReadHDF5header(...) & HDF5ReadROI(...), will work too
//	Apr 11, 2016, version 0.34		renamed HDF5DataInfo to HDF5DataInfoLaueGo
//	Jan 12, 2017, version 0.35		also from /entry1/user, read:  {email, proposal, phone, affiliation, badge}
//	Apr 26, 2018, version 0.36		fixed multibyte writing problem
//	Jun  4, 2021, version 0.36		with Igor 9, change HDF5DataInfoLaueGo --> HDF5DataInfo (now in built in to Igor)

Static Constant SKIP_FIRST_N = 2		// skip the first SKIP_FIRST_N points in a vector

strConstant HDFfileFilters = "HDF Files (*.h5, *.hdf5, *.hdf):.h5,.hdf5,.hdf,;All Files:.*;"
strConstant HDFfileFiltersStrict = "HDF Files (*.h5, *.hdf5, *.hdf):.h5,.hdf5,.hdf,;"


Menu "Load Waves"
	"   Load Jon's HDF Image...",LoadHDF5imageFile("")
End

// units for the metadata that usually accompanies an image
Strconstant HDF5_MetaUnits ="X1:µm;Y1:µm;Z1:µm;H1:µm;F1:µm;X2:µm;Y2:µm;Z2:µm;H2:µm;F2:µm;wirebaseX:µm;wirebaseY:µm;wirebaseZ:µm;AerotechH:µm;depth:µm;keV:keV;undulatorGap:mm;undulatorTaper:mm;exposure:s;ScalerClockFreq:Hz;HutchTemperature:¡C;ringCurrent:mA;sampleDistance:µm;ScalerCountTime:s;"


Static Constant H5S_MAX_RANK = 32
Static Constant kHDF5DataInfoVersion = 1000		// 1000 means 1.000.
//Static Structure HDF5DataInfo					// Use with HDF5DatasetInfo and HDF5AttributeInfo functions
#if IgorVersion() < 9
Structure HDF5DataInfoLaueGo			//rxadd
	// Input fields (inputs to HDF5 XOP)
	uint32 version							// Must be set to kHDF5DataInfoVersion
	char structName[16]					// Must be "HDF5DataInfo".

	// Output fields (outputs from HDF5 XOP)
	double datatype_class;				// e.g., H5T_INTEGER, H5T_FLOAT.
	char datatype_class_str[32];		// String with class spelled out. e.g., "H5T_INTEGER", "H5T_FLOAT".
	double datatype_size;				// Size in bytes of one element.
	double datatype_sign;				// H5T_SGN_NONE (unsigned), H5T_SGN_2 (signed), H5T_SGN_ERROR (this type does not have a sign, i.e., it is not an integer type).
	double datatype_order;				// H5T_ORDER_LE, H5T_ORDER_BE, H5T_ORDER_VAX
	char datatype_str[64];				// Human-readable string, e.g., "16-bit unsigned integer"
	double dataspace_type;				// H5S_NO_CLASS (-1), H5S_SCALAR (0), H5S_SIMPLE (1), H5S_COMPLEX (2).
	double ndims;							// Zero for H5S_SCALAR. Number of dimensions in the dataset for H5S_SIMPLE.
	double dims[H5S_MAX_RANK];		// Size of each dimension.
	double maxdims[H5S_MAX_RANK];	// Maximum size of each dimension.
EndStructure
//
Static Function InitHDF5DataInfo(di)		// Sets input fields.
	STRUCT HDF5DataInfoLaueGo &di
	// HDF5XOP uses these fields to make sure the structure passed in to it is compatible.
	di.version = kHDF5DataInfoVersion
	di.structName = "HDF5DataInfo"
End
#endif

Function/S LoadHDF5imageFile(fName,[extras])
	String fName													// fully qualified name of file to open
	String extras
	extras = SelectString(ParamIsDefault(extras),extras,"")
	Wave dark = $StringByKey("dark",extras)				// optional darkImage to subtract (should match full image, not roi)
	Variable slice=NumberByKey("slice",extras)			// slice (only used for mult-image, i.e. 3D arrays)
	Variable noHeader=NumberByKey("noHeader",extras)
	noHeader = numtype(noHeader) ? 0 : !(!noHeader)	// do not read the header, just the image
	Wave roi=str2vec(StringByKey("roi",extras,"="))	// a 4-vector with {iLo,iHi,jLo,jHi}, or a NULL wave

	Variable f, askSlice=(numtype(slice)>0)
	if (strlen(ParseFilePath(3,fName,":",0,0))<1)		// call dialog if no file name passed
		Open /D/M="HDF5 file"/R/F=HDFfileFilters f		// use /D to get full path name
		fName = S_filename
		askSlice = 1
		// if had to ask for image, also ask about dark
		if (strlen(WaveListClass("darkImage","*","DIMS:2"))>1 && !WaveExists(dark))
			String darkName="_none_"
			Prompt darkName,"name of wave with dark image",popup,"_none_;"+WaveListClass("darkImage","*","DIMS:2")
			DoPrompt "Dark Image",darkName
			if (V_flag)
				return ""
			endif
			Wave dark = $ReplaceString("_none_",darkName,"")
		endif
	endif
	if (strlen(fName)<1)			// no file name, quit
		return ""
	endif

#if IgorVersion() < 9
	STRUCT HDF5DataInfoLaueGo di
#else
	STRUCT HDF5DataInfo di			// Defined in HDF5 Browser.ipf.
#endif
	InitHDF5DataInfo(di)			// Initialize structure.
	HDF5OpenFile/R f as fName
	Variable dummy = HDF5DatasetInfo(f,"/entry1/data/data",0,di)
	HDF5CloseFile f
	Variable nx,ny,Nslices, ndims=di.ndims
	if (ndims!=3 && ndims!=2)		// must be dim 2 or 3
		return ""
	endif
	if (ndims==2)
		nx = di.dims[1]
		ny = di.dims[0]
		Nslices = 1
		slice = 0						// only choice for 2d data file
		askSlice = 0
	else
		nx = di.dims[2]
		ny = di.dims[1]
		Nslices = di.dims[0]
		if (slice >= Nslices)
			return ""
		endif
	endif

	if (Ndims==3 && askSlice)
		slice = numtype(slice) ? 0 : slice
		String str
		sprintf str,"image index, allowed range = [0,%d]",Nslices-1
		Prompt slice,str
		DoPrompt "Slice in 3D array",slice
		if (V_flag)
			return ""
		endif
		slice = round(slice)
		if (slice<0 || slice>=Nslices)
			sprintf str,"slice=%g, allowed range = [0,%d]",slice,Nslices-1
			DoAlert 0,str
			return ""
		endif
	endif

	String wnote=""
	if (!noHeader)
		wnote = ReadHDF5header(fName,extras=extras)
	endif
	HDF5OpenFile/R/Z f as fName
	if (V_Flag)
		return ""
	endif

	Make/FREE/N=(3,4) slab=0
	slab[0][0]= {slice,0,0}		// start
	slab[][1]= 1						// stride
	slab[][2]= 1						// count
	if (ndims==2)
		slab[0][3]= {ny,nx,0}		// block
		slab[2][] = 0
	else
		slab[0][3]= {1,ny,nx}		// block
	endif

	String wName = ParseFilePath(3, fName,":",0,0)
	wName = CleanupName(wName,0)						// wave name based on the file name
	wName += SelectString(ndims==3,"","_S"+num2istr(slice))
//	if (exists(wName))										// if wave already exists, create unique name
//		wName = UniqueName( wName+"_", 1, 1)			// by commenting this out we get overwriting
//	endif


//print slab
//printf "nx=%d, ny=%d, Nslices=%d\r",nx,ny,Nslices
//print " di.dims = ",di.dims[0],di.dims[1],di.dims[2],di.dims[3]
	if (ndims==3)
		HDF5LoadData/IGOR=-1/N=$wName/O/SLAB=slab/TYPE=2/Q/Z f,"/entry1/data/data"  
	else
		HDF5LoadData/IGOR=-1/N=$wName/O/TYPE=2/Q/Z f,"/entry1/data/data"  
	endif
	if (V_flag)
		print "V_Flag=",V_Flag
		return ""
	endif
	HDF5CloseFile f

	Wave image = $StringFromLIst(0,S_waveNames)
	Redimension/N=(nx*ny) image
	Redimension/N=(ny,nx) image
	MatrixOp/O image = image^t		// WHY !!!!!!
	if (WaveExists(dark))
		Redimension/I image								// change from 16bit unsigned to 32 bit signed
		image -= dark
		if (!noHeader)
			wnote = ReplaceStringByKey("bkg",wnote,NameOfWave(dark),"=")
		endif
	endif

	Variable i0=0,j0=0										// start of the  sub-region (in binned pixels)
	if (WaveExists(roi))
		i0 = roi[0]												// start of sub-region
		j0 = roi[2]	
		Variable i1=roi[1], j1=roi[3]					// end of sub-region
		nx = DimSize(image,0)
		ny = DimSize(image,1)
		i1 = (i1==-1) ? nx-1 : i1							// "-1" means the last point
		j1 = (j1==-1) ? ny-1 : j1
		if (i1>=nx || j1>=ny || i1<i0 || j1<j0 || i0<0 || j0<0 || numtype(i0+j0+i1+j1))
			printf "ERROR -- LoadHDF5imageFile(), the given roi = [%g,%g][%g,%g] does not fit within the image\r",i0,i1,j0,j1
			return ""
		elseif (i0>0 || j0>0 || i1<(nx-1) || j1<(ny-1))
			Duplicate/FREE image, imageOrig				// store the original image (temporarily)
			Redimension/N=(i1-i0+1, j1-j0+1) image	// re-dimension original image to the roi
			image = imageOrig[i0+p][j0+q]				// copy original data back into smaller image
			WaveClear imageOrig
		endif
	endif

	if (!noHeader && Ndims==3 && slice>=0)			// extract position of an individual slice (if we are extracting one slice)
		wnote = ReplaceNumberByKey("slice",wnote,slice,"=")
		str = StringFromList(slice,StringByKey("X2", wnote,"="),",")
		if (strlen(str))
			wnote = ReplaceStringByKey("X2",wnote,str,"=")
		endif
		str = StringFromList(slice,StringByKey("Y2", wnote,"="),",")
		if (strlen(str))
			wnote = ReplaceStringByKey("Y2",wnote,str,"=")
		endif
		str = StringFromList(slice,StringByKey("Z2", wnote,"="),",")
		if (strlen(str))
			wnote = ReplaceStringByKey("Z2",wnote,str,"=")
		endif
		str = StringFromList(slice,StringByKey("AerotechH", wnote,"="),",")
		if (strlen(str))
			wnote = ReplaceStringByKey("AerotechH",wnote,str,"=")
		endif
		str = StringFromList(slice,StringByKey("H2", wnote,"="),",")
		if (strlen(str))
			wnote = ReplaceStringByKey("H2",wnote,str,"=")
		endif
		str = StringFromList(slice,StringByKey("F2", wnote,"="),",")
		if (strlen(str))
			wnote = ReplaceStringByKey("F2",wnote,str,"=")
		endif
		str = StringFromList(slice,StringByKey("keV", wnote,"="),",")
		if (strlen(str))
			wnote = ReplaceStringByKey("keV",wnote,str,"=")
		endif
	endif

//		Not needed, the header is alredy adjusted in 
//	if (!noHeader && (i0!=0 || j0!=0))					// sub-region of image does not start at (0,0)
//		Variable xdim=DimSize(image,0), ydim=DimSize(image,1)
//		Variable startx, endx, groupx, starty, endy, groupy
//		startx	= NumberByKey("startx",wnote,"=")
//		endx	= NumberByKey("endx",wnote,"=")
//		groupx	= NumberByKey("groupx",wnote,"=")
//		starty	= NumberByKey("starty",wnote,"=")
//		endy	= NumberByKey("endy",wnote,"=")
//		groupy	= NumberByKey("groupy",wnote,"=")
//		startx += i0*groupx									// re-set to match the requested sub-region read
//		starty += j0*groupy
//		endx = xdim*groupx+startx-1
//		endy = ydim*groupy+starty-1
//		wnote = ReplaceNumberByKey("xdim",wnote,xdim,"=")
//		wnote = ReplaceNumberByKey("ydim",wnote,ydim,"=")
//		wnote = ReplaceNumberByKey("startx",wnote,startx,"=")
//		wnote = ReplaceNumberByKey("endx",wnote,endx,"=")
//		wnote = ReplaceNumberByKey("starty",wnote,starty,"=")
//		wnote = ReplaceNumberByKey("endy",wnote,endy,"=")
//		wnote = ReplaceNumberByKey("startxRead",wnote,i0,"=")
//		wnote = ReplaceNumberByKey("startyRead",wnote,j0,"=")
//	endif
	SetScale/P x 0,1,"pixel", image
	SetScale/P y 0,1,"pixel", image
	Note/K image,wnote

	//	String wName = WinViewReadROI(fName,0,-1,0,-1)	// load file into wName
	if (ItemsInList(GetRTStackInfo(2))==0)
		printf "total length = %d x %d  = %d points\r", DimSize(image,0),DimSize(image,1),numpnts(image)
		print "number type is  '"+IgorFileTypeString(NumberByKey("NUMTYPE",WaveInfo(image,0)))+"'"
		if (WaveExists(dark))
			printf "subtracted the dark image  '%s'\r",NameOfWave(dark)
		endif

		print "Created a 2-d wave    '"+NameOfWave(image)+"'"
		DoAlert 1, "Display this image"
		if (V_Flag==1)
			FUNCREF NewImageGraphHDF5proto fnew= $"NewImageGraph"
			fnew(image)
		endif
	endif
	return GetWavesDataFolder(image,2)
End
//
Function/S NewImageGraphHDF5proto(image,[withButtons])
	Wave image
	Variable withButtons
	return ""
End



//strConstant testfileNameStr = "Macintosh HD:Users:Shared:temp_tischler:dev:reconstructC_Big:Ge-fly_4.h5"
//Function test_LoadHDF5imageFile(fname,slice)
//	String fname
//	Variable slice
//	Variable f
//	Open/R/D=2/F=HDFfileFilters/M="pick multi-image HDF5 file"/P=home f as fname
//	fname = S_fileName
//	String imageName =LoadHDF5imageFile(fName,slice=slice)
//	printf "from file '%s',   loaded image '%s'\r",fname,imageName
//End
//
//
//Function/S LoadHDF5imageFile(fName,[dark])
//	String fName												// fully qualified name of file to open
//	Wave dark													// optional darkImage to subtract
//	if (ParamIsDefault(dark))
//		Wave dark = $""
//	endif
//
//	Variable f
//	if (strlen(ParseFilePath(3,fName,":",0,0))<1)		// call dialog if no file name passed
////		Open /D/M="HDF5 file"/R/T=".h5 " f				// use /D to get full path name
//		Open /D/M="HDF5 file"/R/F=HDFfileFilters f		// use /D to get full path name
//		fName = S_filename
//		// if had to ask for image, also ask about dark
//		if (strlen(WaveListClass("darkImage","*","DIMS:2"))>1 && !WaveExists(dark))
//			String darkName="_none_"
//			Prompt darkName,"name of wave with dark image",popup,"_none_;"+WaveListClass("darkImage","*","DIMS:2")
//			DoPrompt "Dark Image",darkName
//			if (V_flag)
//				return ""
//			endif
//			Wave dark = $ReplaceString("_none_",darkName,"")
//		endif
//	endif
//	if (strlen(fName)<1)										// no file name, quit
//		return ""
//	endif
//
//	String wnote = ReadHDF5header(fName)
//	HDF5OpenFile/R/Z f as fName
//	if (V_Flag)
//		return ""
//	endif
//
////	HDF5ListGroup/F/TYPE=15 f, "/entry1/data"
////print "V_Flag=",V_Flag
////print "S_HDF5ListGroup = ",S_HDF5ListGroup
//
//	HDF5ListGroup/F/TYPE=15 f, "/entry1/data"
////print "V_Flag=",V_Flag
////print "S_HDF5ListGroup = ",S_HDF5ListGroup
//
//	String dataName = StringFromList(0,S_HDF5ListGroup)
////	print "dataName = ",dataName
////	HDF5LoadGroup/IGOR=-1/L=7/T=HDF5data  : , f , "/entry1/data" 
////	print "V_Flag =",V_Flag
////	print "S_groupPaths=",S_groupPaths
////	print "S_dataFolderPaths=",S_dataFolderPaths
////	print "S_objectPaths=",S_objectPaths
////	HDF5OpenGroup [ /Z] f , "entry" , groupID
//
//
////	String wName = ParseFilePath(3, fName,":",0,0)+"_h5"
//	String wName = ParseFilePath(3, fName,":",0,0)
//	wName = CleanupName(wName,0)						// wave name based on the file name
//	if (exists(wName))										// if wave already exists, create unique name
//		wName = UniqueName( wName+"_", 1, 1)
//	endif
//	HDF5LoadData/IGOR=-1/N=$wName/O/TYPE=2/Q/Z f,"/entry1/data/data"  
//	if (V_flag)
//		printf "in LoadHDF5imageFile(), looking for '/entry1/data/data',    got V_Flag=%g\r",V_Flag
////		print "V_Flag=",V_Flag
//		return ""
//	endif
//	Wave image = $StringFromLIst(0,S_waveNames)
//	MatrixOp/O image = image^t		// WHY !!!!!!
//	HDF5CloseFile f
//
//	Variable i0=0,j0=0										// start of sub-region (binned pixels)
//	if (i0!=0 || j0!=0)										// sub-region of image does not start at (0,0)
//		Variable xdim=DimSize(image,0), ydim=DimSize(image,1)
//		Variable startx, endx, groupx, starty, endy, groupy
//		startx	= NumberByKey("startx",wnote,"=")
//		endx	= NumberByKey("endx",wnote,"=")
//		groupx	= NumberByKey("groupx",wnote,"=")
//		starty	= NumberByKey("starty",wnote,"=")
//		endy	= NumberByKey("endy",wnote,"=")
//		groupy	= NumberByKey("groupy",wnote,"=")
//		startx += i0*groupx									// re-set to match the sub-region read
//		starty += j0*groupy
//		endx = xdim*groupx+startx-1
//		endy = ydim*groupy+starty-1
//		wnote = ReplaceNumberByKey("xdim",wnote,xdim,"=")
//		wnote = ReplaceNumberByKey("ydim",wnote,ydim,"=")
//		wnote = ReplaceNumberByKey("startx",wnote,startx,"=")
//		wnote = ReplaceNumberByKey("endx",wnote,endx,"=")
//		wnote = ReplaceNumberByKey("starty",wnote,starty,"=")
//		wnote = ReplaceNumberByKey("endy",wnote,endy,"=")
//		wnote = ReplaceNumberByKey("startxRead",wnote,i0,"=")
//		wnote = ReplaceNumberByKey("startyRead",wnote,j0,"=")
//	endif
//	SetScale/P x 0,1,"pixel", image
//	SetScale/P y 0,1,"pixel", image
//	if (WaveExists(dark))
//		Redimension/I image
//		image -= dark
//		wnote = ReplaceStringByKey("bkg",wnote,NameOfWave(dark),"=")
//	endif
//	Note/K image,wnote
//
////	HDF5OpenGroup [ /Z] locationID , nameStr , groupID
//
////	String wName = WinViewReadROI(fName,0,-1,0,-1)	// load file into wName
//	if (ItemsInList(GetRTStackInfo(0))<=1)
//		printf "total length = %d x %d  = %d points\r", DimSize(image,0),DimSize(image,1),numpnts(image)
//		print "number type is  '"+IgorFileTypeString(NumberByKey("NUMTYPE",WaveInfo(image,0)))+"'"
//		if (WaveExists(dark))
//			printf "subtracted the dark image  '%s'\r",NameOfWave(dark)
//		endif
//
//		print "Created a 2-d wave    '"+NameOfWave(image)+"'"
//		DoAlert 1, "Display this image"
//		if (V_Flag==1)
//			NewImageGraph(image)
//		endif
//	endif
//	return GetWavesDataFolder(image,2)
//End
//
// NOTE, Igor does not support this, it is just here for completeness
// This will return an roi, but it reads the whole image and then takes a piece of it
Function/T HDF5ReadROI(fileName,i0,i1,j0,j1,[extras])
	String fileName					// fully qualified name of file to open (will not prompt)
	Variable i0,i1,j0,j1				// pixel range of ROI (if i1 or j1<0 then use whole image)
	String extras
	extras = SelectString(ParamIsDefault(extras),extras,"")
	String str
	sprintf str, "%g,%g,%g,%g",i0,i1,j0,j1
	extras = ReplaceStringByKey("roi",extras,str,"=")
	if (strlen(extras))
		return LoadHDF5imageFile(fileName,extras=extras)
	else
		return LoadHDF5imageFile(fileName)
	endif
End



Function/T ReadHDF5header(fName,[extras])
	String fName					// fully qualified name of file to open (will not prompt)
	String extras					// optional switches (only supports EscanOnly in this routine)
	extras = SelectString(ParamIsDefault(extras),extras,"")

	if (strlen(fName)<1)
		return ""
	endif

	Variable EscanOnly=NumberByKey("EscanOnly",extras)
	EscanOnly = numtype(EscanOnly) ? 0 : !(!EscanOnly)	// only read "X1;Y1;Z1;depth;keV;startx;endx;groupx;starty;endy;groupy;I0;I0gain;ScalerCountTime;exposure;

	String wnote="", str, model=""
	Variable f, value

	wnote = ReplaceStringByKey("waveClass",wnote,"rawImage","=")
	wnote = ReplaceStringByKey("imageFileName",wnote,ParseFilePath(0,fName,":",1,0),"=")
	wnote = ReplaceStringByKey("imageFilePath",wnote,ParseFilePath(1,fName,":",1,0),"=")

	HDF5OpenFile/R/Z f as fName
	if (V_Flag)
		return ""
	endif

#if IgorVersion() < 9
	STRUCT HDF5DataInfoLaueGo di	// Defined in HDF5 Browser.ipf.
#else
	STRUCT HDF5DataInfo di	// Defined in HDF5 Browser.ipf.
#endif
	InitHDF5DataInfo(di)			// Initialize structure.
	Variable dummy = HDF5DatasetInfo(f,"/entry1/data/data",0,di)
	Variable Nslices = (di.ndims == 3) ? di.dims[0] : 0	// 0 means 2-d data
	if (Nslices)
		wnote= ReplaceNumberByKey("Nslices",wnote,round(Nslices),"=")
	endif

	str = getStrVecHDF5dataNum(f,"entry1/depth",places=2,Nww=Nslices)
	if (strlen(str))
		wnote= ReplaceStringByKey("depth",wnote,str,"=")
	endif

	Variable reconstructed = Nslices==1 && numtype(str2num(str))==0		// a single reconstructed image from flyScan, 1 image + depth
	if (!reconstructed && !EscanOnly)
		// wire positions
		HDF5ListGroup/Z f , "/entry1/wire"
		String wireFolder = SelectString(V_flag,"entry1/wire/","entry1/")
		str = getStrVecHDF5dataNum(f,wireFolder+"wireX",places=2,Nww=Nslices)
		if (strlen(str))
			wnote= ReplaceStringByKey("X2",wnote,str,"=")
		endif
		str = getStrVecHDF5dataNum(f,wireFolder+"wireY",places=2,Nww=Nslices)
		if (strlen(str))
			wnote= ReplaceStringByKey("Y2",wnote,str,"=")
			Wave y2wave = $getWaveHDF5dataNum(f,wireFolder+"wireY",Nww=Nslices)
		endif
		str = getStrVecHDF5dataNum(f,wireFolder+"wireZ",places=2,Nww=Nslices)
		if (strlen(str))
			wnote= ReplaceStringByKey("Z2",wnote,str,"=")
			Wave z2wave = $getWaveHDF5dataNum(f,wireFolder+"wireZ",Nww=Nslices)
		endif
		value = get1HDF5dataNum(f,wireFolder+"wirebaseX")
		if (numtype(value)==0)
			value = round(value*1000)/1000			// only three places for wirebase
			wnote= ReplaceNumberByKey("wirebaseX",wnote,value,"=")
		endif
		value = get1HDF5dataNum(f,wireFolder+"wirebaseY")
		if (numtype(value)==0)
			value = round(value*1000)/1000
			wnote= ReplaceNumberByKey("wirebaseY",wnote,value,"=")
		endif
		value = get1HDF5dataNum(f,wireFolder+"wirebaseZ")
		if (numtype(value)==0)
			value = round(value*1000)/1000
			wnote= ReplaceNumberByKey("wirebaseZ",wnote,value,"=")
		endif
		str = getStrVecHDF5dataNum(f,wireFolder+"AerotechH",places=3,Nww=Nslices)
		if (strlen(str))
			wnote= ReplaceStringByKey("AerotechH",wnote,str,"=")
			// Wave Aerotechwave = $getWaveHDF5dataNum(f,wireFolder+"wireY")
		endif
	endif

	// sample positions
	str = getStrVecHDF5dataNum(f,"entry1/sample/sampleX",places=2,Nww=Nslices)
	if (strlen(str))
		wnote= ReplaceStringByKey("X1",wnote,str,"=")
	endif
	str = getStrVecHDF5dataNum(f,"entry1/sample/sampleY",places=2,Nww=Nslices)
	if (strlen(str))
		wnote= ReplaceStringByKey("Y1",wnote,str,"=")
		Wave y1wave = $getWaveHDF5dataNum(f,"entry1/sample/sampleY",Nww=Nslices)
	endif
	str = getStrVecHDF5dataNum(f,"entry1/sample/sampleZ",places=2,Nww=Nslices)
	if (strlen(str))
		wnote= ReplaceStringByKey("Z1",wnote,str,"=")
		Wave z1wave = $getWaveHDF5dataNum(f,"entry1/sample/sampleZ",Nww=Nslices,Nww=Nslices)
	endif

	// detector pixels & ROI
	Variable startx,endx,groupx,xdim
	startx = get1HDF5dataNum(f,"entry1/detector/startx")
	if (numtype(startx)==0)
		wnote= ReplaceNumberByKey("startx",wnote,round(startx),"=")
	endif
	endx = get1HDF5dataNum(f,"entry1/detector/endx")
	if (numtype(endx)==0)
		wnote= ReplaceNumberByKey("endx",wnote,round(endx),"=")
	endif
	groupx = get1HDF5dataNum(f,"entry1/detector/binx")
	if (numtype(groupx)==0)
		wnote= ReplaceNumberByKey("groupx",wnote,round(groupx),"=")
	endif
	xdim = (endx-startx+1)/groupx
	if (numtype(xdim)==0)
		wnote= ReplaceNumberByKey("xdim",wnote,round(xdim),"=")
	endif
	if (!EscanOnly)
		value = get1HDF5dataNum(f,"entry1/detector/Nx")
		if (numtype(value)==0)
			wnote= ReplaceNumberByKey("xDimDet",wnote,round(value),"=")
		endif
	endif

	Variable starty,endy,groupy, ydim
	starty = get1HDF5dataNum(f,"entry1/detector/starty")
	if (numtype(starty)==0)
		wnote= ReplaceNumberByKey("starty",wnote,round(starty),"=")
	endif
	endy = get1HDF5dataNum(f,"entry1/detector/endy")
	if (numtype(endy)==0)
		wnote= ReplaceNumberByKey("endy",wnote,round(endy),"=")
	endif
	groupy = get1HDF5dataNum(f,"entry1/detector/biny")
	if (numtype(groupy)==0)
		wnote= ReplaceNumberByKey("groupy",wnote,round(groupy),"=")
	endif
	ydim = (endy-starty+1)/groupy
	if (numtype(ydim)==0)
		wnote= ReplaceNumberByKey("ydim",wnote,round(ydim),"=")
	endif
	if (!EscanOnly)
		value = get1HDF5dataNum(f,"entry1/detector/Ny")
		if (numtype(value)==0)
			wnote= ReplaceNumberByKey("yDimDet",wnote,round(value),"=")
		endif
	endif

	wnote = ResetHeaderToROI(wnote,extras)		// this corrects for when extras contains an roi

	if (!EscanOnly)
		// strings: detectorID, Model, beamLine, title, userName, sampleName, file_name, file_time
		str = get1HDF5dataStr(f,"entry1/detector/ID")
		if (strlen(str))
			wnote= ReplaceStringByKey("detectorID",wnote,str,"=")
		endif
		str = get1HDF5dataStr(f,"entry1/detector/model")
		if (strlen(str))
			model = str
			wnote= ReplaceStringByKey("detectorModel",wnote,model,"=")
		endif
		str = get1HDF5dataStr(f,"Facility/facility_beamline")
		if (strlen(str))
			wnote= ReplaceStringByKey("beamLine",wnote,str,"=")
		endif
		str = get1HDF5dataStr(f,"entry1/title")
		if (strlen(str))
			wnote= ReplaceStringByKey("title",wnote,str,"=")
		endif

		str = get1HDF5dataStr(f,"entry1/user/name")			// read the "/entry1/user" group
		if (strlen(str))
			wnote= ReplaceStringByKey("userName",wnote,str,"=")
		endif
		str = get1HDF5dataStr(f,"entry1/user/email")
		if (strlen(str))
			wnote= ReplaceStringByKey("email",wnote,str,"=")
		endif
		str = get1HDF5dataStr(f,"entry1/user/proposal")
		if (strlen(str))
			wnote= ReplaceStringByKey("proposal",wnote,str,"=")
		endif
		str = get1HDF5dataStr(f,"entry1/user/telephone_number")
		if (strlen(str))
			wnote= ReplaceStringByKey("phone",wnote,str,"=")
		endif
		str = get1HDF5dataStr(f,"entry1/user/affiliation")
		if (strlen(str))
			wnote= ReplaceStringByKey("affiliation",wnote,str,"=")
		endif
		value = get1HDF5dataNum(f,"entry1/user/facility_user_id")
		if (numtype(value)==0)
			wnote= ReplaceNumberByKey("badge",wnote,round(value),"=")
		endif

		str = get1HDF5dataStr(f,"entry1/sample/name")
		if (strlen(str))
			wnote= ReplaceStringByKey("sampleName",wnote,str,"=")
		endif
		str = get1HDF5AttrStr(f,"/","file_name")
		if (strlen(str))
			wnote= ReplaceStringByKey("file_name",wnote,str,"=")
		endif
		str = get1HDF5AttrStr(f,"/","file_time")
		if (strlen(str))
			wnote= ReplaceStringByKey("file_time",wnote,str,"=")
		endif
	endif

	// monitor  I0, I_start, I_final, I0_calc, I_final_calc, I_start_calc, ScalerClockFreq, ScalerClock_calc, ScalerCountTime
	str = getStrVecHDF5dataNum(f,"entry1/monitor/I0",Nww=Nslices,places=0)
	if (strlen(str))
		wnote= ReplaceStringByKey("I0",wnote,str,"=")
	endif
	str = getStrVecHDF5dataNum(f,"entry1/monitor/ScalerCountTime",Nww=Nslices,places=5)
	if (strlen(str))
		wnote= ReplaceStringByKey("ScalerCountTime",wnote,str,"=")
	endif

	if (!EscanOnly)
		str = getStrVecHDF5dataNum(f,"entry1/monitor/I_start",Nww=Nslices,places=0)
		if (strlen(str))
			wnote= ReplaceStringByKey("Istart",wnote,str,"=")
		endif
		str = getStrVecHDF5dataNum(f,"entry1/monitor/I_final",Nww=Nslices,places=0)
		if (strlen(str))
			wnote= ReplaceStringByKey("Ifinal",wnote,str,"=")
		endif
		str = getStrVecHDF5dataNum(f,"entry1/monitor/I0_calc",Nww=Nslices,places=5)
		if (strlen(str))
			wnote= ReplaceStringByKey("I0_calc",wnote,str,"=")
		endif
		str = getStrVecHDF5dataNum(f,"entry1/monitor/I_final_calc",Nww=Nslices,places=5)
		if (strlen(str))
			wnote= ReplaceStringByKey("Ifinal_calc",wnote,str,"=")
		endif
		str = getStrVecHDF5dataNum(f,"entry1/monitor/I_start_calc",Nww=Nslices,places=0)
		if (strlen(str))
			wnote= ReplaceStringByKey("Istart_calc",wnote,str,"=")
		endif
		str = getStrVecHDF5dataNum(f,"entry1/monitor/ScalerClockFreq",Nww=Nslices,places=0)
		if (strlen(str))
			wnote= ReplaceStringByKey("ScalerClockFreq",wnote,str,"=")
		endif
		str = getStrVecHDF5dataNum(f,"entry1/monitor/ScalerClock_calc",Nww=Nslices,places=5)
		if (strlen(str))
			wnote= ReplaceStringByKey("ScalerClock_calc",wnote,str,"=")
		endif
	endif

	// energy, exposure, scanNum, sampleDistance, detectorGain
	str = getStrVecHDF5dataNum(f,"entry1/sample/incident_energy",places=4,Nww=Nslices)
	if (strlen(str))
		wnote= ReplaceStringByKey("keV",wnote,str,"=")
	endif

	value = get1HDF5dataNum(f,"entry1/detector/exposure")
	if (numtype(value)==0)
//		// This is for the older software, the new ones write the actual exposure time
//		if (value == limit(round(value),0,7))		// special for the silly Perkin Elmer detectors
//			if (strsearch(model,"XRD1621",0)<=0)
//				value = str2num(StringFromList(value,"66.5;82.7;99.9;124.9;166.6;249.8;499.8;999.9;"))*1e-3
//			elseif (strsearch(model,"XRD0820",0)<=0)
//				value = str2num(StringFromList(value,"66.5;79.9;99.8;133.2;199.9;400.0;999.8;1999.8;"))*1e-3
//			endif
//		endif
		wnote= ReplaceNumberByKey("exposure",wnote,value,"=")
	endif

	if (!EscanOnly)
		value = get1HDF5dataNum(f,"entry1/detector/gain")
		if (numtype(value)==0)
			str = num2str(value)
			if (value == limit(round(value),0,5))		// special for the silly Perkin Elmer detectors
				str = StringFromList(value,"0.25;0.5;1;2;4;8;")+"pF"
			endif
			wnote= ReplaceStringByKey("detectorGain",wnote,str,"=")
		endif
		value = get1HDF5dataNum(f,"entry1/scanNum")
		if (numtype(value)==0)
			wnote= ReplaceNumberByKey("scanNum",wnote,value,"=")
		endif
		value = get1HDF5dataNum(f,"entry1/sample/distance")
		if (numtype(value)==0)
			value = 1e-2 *round(value*1e5)	// convert from mm to µm
			wnote= ReplaceNumberByKey("sampleDistance",wnote,value,"=")
		endif
	endif

	if (!EscanOnly)
		// BeamBad, CCDshutter, HutchTemperature, LightOn, MonoMode
		value = get1HDF5dataNum(f,"entry1/microDiffraction/BeamBad")
		if (numtype(value)==0)
			wnote= ReplaceNumberByKey("BeamBad",wnote,value,"=")
		endif

		value = get1HDF5dataNum(f,"entry1/microDiffraction/CCDshutter")
		if (numtype(value)==0)
			str = num2str(value)							// should be 0 or 1
			str = SelectString(value==0,str,"in")		// 0 means in
			str = SelectString(value==1,str,"out")	// 1 means out
			wnote= ReplaceStringByKey("CCDshutter",wnote,str,"=")
		endif

		value = get1HDF5dataNum(f,"entry1/microDiffraction/HutchTemperature")
		if (numtype(value)==0)
			value = round(value*100)/100
			wnote= ReplaceNumberByKey("HutchTemperature",wnote,value,"=")
		endif

		value = get1HDF5dataNum(f,"entry1/microDiffraction/LightOn")
		if (numtype(value)==0)
			wnote= ReplaceNumberByKey("LightOn",wnote,value,"=")
		endif
		str = get1HDF5dataStr(f,"entry1/microDiffraction/MonoMode")
		if (strlen(str))
			wnote= ReplaceStringByKey("MonoMode",wnote,str,"=")
		endif

		value = get1HDF5dataNum(f,"entry1/microDiffraction/source/current")
		if (numtype(value)==0)
			wnote= ReplaceNumberByKey("ringCurrent",wnote,value,"=")
		endif
		str = get1HDF5dataStr(f,"entry1/microDiffraction/source/gap")
		if (strlen(str))
			wnote= ReplaceStringByKey("undulatorGap",wnote,str,"=")
		endif
		str = get1HDF5dataStr(f,"entry1/microDiffraction/source/taper")
		if (strlen(str))
			wnote= ReplaceStringByKey("undulatorTaper",wnote,str,"=")
		endif
		str = get1HDF5dataStr(f,"entry1/microDiffraction/source/top_up")
		if (strlen(str))
			wnote= ReplaceStringByKey("topUp",wnote,str,"=")
		endif
	endif

	if (!EscanOnly)
		if (exists("ReadHDF5headerMore")==6)
			FUNCREF ReadHDF5headerMoreProto  func=$"ReadHDF5headerMore"
			wnote = func(f,wnote)
		endif
	endif

	HDF5CloseFile f

#if (Exists("YZ2H")==6)
	// optionally add H&F from Y&Z
	//	if (numtype(Y1+Z1)==0)
	//		wnote = ReplaceNumberByKey("H1",wnote,1e-3*round(YZ2H(Y1,Z1)*1e3),"=")
	//		wnote = ReplaceNumberByKey("F1",wnote,1e-3*round(YZ2F(Y1,Z1)*1e3),"=")
	//	elseif (WaveExists(y1wave) && WaveExists(z1wave))
	if (WaveExists(y1wave) && WaveExists(z1wave))			// add sample positions H1 & F1 if Y1 & Z1 exist
		Duplicate/FREE y1wave,f1wave,h1wave
		h1wave = YZ2H(y1wave[p],z1wave[p])
		f1wave = YZ2F(y1wave[p],z1wave[p])
		wnote = ReplaceStringByKey("H1",wnote,vec2str(h1wave,places=3),"=")
		wnote = ReplaceStringByKey("F1",wnote,vec2str(f1wave,places=3),"=")
	endif
	//	if (numtype(Y2+Z2)==0)
	//		wnote = ReplaceNumberByKey("H2",wnote,1e-3*round(YZ2H(Y2,Z2)*1e3),"=")
	//		wnote = ReplaceNumberByKey("F2",wnote,1e-3*round(YZ2F(Y2,Z2)*1e3),"=")
	//	elseif (WaveExists(y2wave) && WaveExists(z2wave))
	if (WaveExists(y2wave) && WaveExists(z2wave))			// add wire positions H2 & F2 if Y2 & Z2 exist
		Duplicate/FREE y2wave,f2wave,h2wave
		h2wave = YZ2H(y2wave[p],z2wave[p])
		f2wave = YZ2F(y2wave[p],z2wave[p])
		wnote = ReplaceStringByKey("H2",wnote,vec2str(h2wave,places=3),"=")
		wnote = ReplaceStringByKey("F2",wnote,vec2str(f2wave,places=3),"=")
	endif
#endif
	KillWaves/Z y1wave, z1wave,y2wave, z2wave
	return wnote
End
//
Static Function/T ResetHeaderToROI(list,extras)
	String list
	String extras

	Wave roi=str2vec(StringByKey("roi",extras,"="))	// a 4-vector with {iLo,iHi,jLo,jHi}, or a NULL wave
	if (!WaveExists(roi))
		return list
	elseif (DimSize(roi,0)<4)				// not enough data, give up
		return list
	endif
	Variable i0=roi[0], i1=roi[1], j0=roi[2], j1=roi[3]

	Variable xdim=NumberByKey("xdim",list,"="), ydim=NumberByKey("ydim",list,"=")		// size of image (NOT detector)
	Variable Nx=NumberByKey("xDimDet",list,"="), Ny=NumberByKey("yDimDet",list,"=")	// full size of detector (NOT image)
	Variable startx=NumberByKey("startx",list,"="), starty=NumberByKey("starty",list,"=")
	Variable endx=NumberByKey("endx",list,"="), endy=NumberByKey("endy",list,"=")
	Variable groupx=NumberByKey("groupx",list,"="), groupy=NumberByKey("groupy",list,"=")
	groupx = numtype(groupx) ? 1 : groupx
	groupy = numtype(groupy) ? 1 : groupy
	xdim = numtype(xdim) ? floor((endx-startx+1)/groupx) : xdim
	ydim = numtype(ydim) ? floor((endy-starty+1)/groupy) : ydim
	xdim = i1 - i0 + 1						// re-set to match the requested sub-region
	ydim = j1 - j0 + 1
	startx += i0*groupx
	starty += j0*groupy
	endx = startx + xdim*groupx - 1
	endy = starty + ydim*groupy - 1
	if (numtype(startx+endx + starty+endy + xdim+ydim) || startx<0 || starty<0 || endx>Nx || endy>Ny)
		return list
	endif

	list = ReplaceNumberByKey("startx",list,startx,"=")
	list = ReplaceNumberByKey("starty",list,starty,"=")
	list = ReplaceNumberByKey("endx",list,endx,"=")
	list = ReplaceNumberByKey("endx",list,endx,"=")
	list = ReplaceNumberByKey("xdim",list,xdim,"=")
	list = ReplaceNumberByKey("ydim",list,ydim,"=")
	return list
End



//Function test_ReadHDF5header(fname)
//	String fname
//	Open/R/D=2/F=HDFfileFilters/M="pick multi-image HDF5 file"/P=home f as fname
//	fname = S_fileName
//	String str = ReadHDF5header(fName)
//	print fname
//	print str
//End
//
//
//		// VO2,   Current, Epoch, Temperature, Volt, Resistance
//		Function/T ReadHDF5headerMore(f,wnote)
//			Variable f
//			String wnote
//		
//			if (f<=0)
//				return wnote
//			endif
//			Variable Voltage, Current, value
//			value = HDF5images#get1HDF5dataNum(f,"entry1/sample/VO2_Current")
//			Current = value
//			if (numtype(value)==0)
//				wnote= ReplaceNumberByKey("VO2_Current",wnote,value,"=")
//			endif
//			value = HDF5images#get1HDF5dataNum(f,"entry1/sample/VO2_Epoch")
//			if (numtype(value)==0)
//				wnote= ReplaceNumberByKey("VO2_Epoch",wnote,value,"=")
//			endif
//			value = HDF5images#get1HDF5dataNum(f,"entry1/sample/VO2_Temp")
//			if (numtype(value)==0)
//				wnote= ReplaceNumberByKey("VO2_Temp",wnote,value,"=")
//			endif
//			value = HDF5images#get1HDF5dataNum(f,"entry1/sample/VO2_Volt")
//			Voltage = value
//			if (numtype(value)==0)
//				wnote= ReplaceNumberByKey("VO2_Volt",wnote,value,"=")
//			endif
//			value = 	Voltage/Current
//			if (numtype(value)==0)
//				wnote= ReplaceNumberByKey("VO2_Reistance",wnote,value,"=")
//			endif
//			return wnote
//		End
//
Function/T ReadHDF5headerMoreProto(f,wnote)
		Variable f
		String wnote
	return wnote
End


//
Static Function get1HDF5dataNum(f,dataName)
	Variable f
	String dataName

	String wName = UniqueName("temp",1,0)
	HDF5LoadData/IGOR=-1/N=$wName/O/TYPE=2/Q/Z f, dataName
	if (V_Flag)
		return NaN
	endif
	Wave ww = $wName
	Variable value = ww[0]
	KillWaves $wName
	return value
End
//
Static Function/T get1HDF5dataStr(f,dataName)
	Variable f
	String dataName

	String wName = UniqueName("temp",1,0), value
	HDF5LoadData/IGOR=-1/N=$wName/O/TYPE=2/Q/Z f, dataName
	if (V_Flag)
		return ""
	endif

	if (WaveType($wName)==0)		// string
		Wave/T wstring = $wName
		value = wstring[0]
	else									// number
		Wave wnum = $wName
		value = num2str(wnum[0])
	endif
	KillWaves $wName
	return value
End
//
// read a vector from HDF file and return it as a comma separated list (not a wave)
Static Function/T getStrVecHDF5dataNum(f,dataName,[places,Nww])
	Variable f
	String dataName
	Variable places
	Variable Nww				// expected length of ww, used with SKIP_FIRST_N

	Wave ww = $getWaveHDF5dataNum(f,dataName,Nww=Nww)
	String str=""
	if (ParamIsDefault(places))
		str = vec2str(ww)
	else
		str = vec2str(ww,places=places)
	endif
	KillWaves ww
	Variable i = strlen(str)-1
	str = SelectString(stringmatch(str[i,i],","),str,str,str[0,i-1])	// remove a trailing comma
	return str
End
//
//Static Function/T getStrVecHDF5dataNum(f,dataName,[places])
//	Variable f
//	String dataName
//	Variable places
//
//	String wName = UniqueName("temp",1,0)
//	HDF5LoadData/IGOR=-1/N=$wName/O/TYPE=2/Q/Z f,dataName
//	if (V_flag)
//		print "V_Flag=",V_Flag
//		return ""
//	endif
//
//	Wave ww = $wName
//	String str=""
//	if (ParamIsDefault(places))
//		str = vec2str(vec)
//	else
//		str = vec2str(vec,places=places)
//	endif
////	String str=""
////	Variable i,N=numpnts(ww)
////	for (i=0;i<N;i+=1)
////		str += num2str(ww[i])+","
////	endfor
//	KillWaves $wName
//	str = SelectString(strlen(str),"",str[0,strlen(str)-2])	// remove trailing comma
//	return str
//End
//
//Static Function/T getStrVecHDF5dataNum(f,dataName,[places])
//	Variable f
//	String dataName
//	Variable places
//
//	Wave ww = $getWaveHDF5dataNum(f,dataName)
//	String str=""
//	if (ParamIsDefault(places))
//		str = vec2str(ww)
//	else
//		str = vec2str(ww,places=places)
//	endif
//	KillWaves ww
//	Variable i = strlen(str)-1
//	str = SelectString(stringmatch(str[i,i],","),str,str,str[0,i-1])	// remove a trailing comma
//	return str
//End
//
Static Function/T vec2str(vec,[places])
	Wave vec
	Variable places
	if (!WaveExists(vec))
		return ""
	endif

	places = ParamIsDefault(places) ? -1 : limit(round(places),0,30)
	String fmt = "%."+num2istr(places)+"f"
	String str,list=""
	Variable i,N=numpnts(vec)
	for (i=0;i<N;i+=1)
		sprintf str,fmt,vec[i]
		list += str+","
	endfor
	list = SelectString(strlen(list),"",list[0,strlen(list)-2])	// remove trailing comma
	return list
End
//Static Function/T vec2str(vec,[places])
//	Wave vec
//	Variable places
//	if (!WaveExists(vec))
//		return ""
//	endif
//
//	places = ParamIsDefault(places) ? -1 : limit(round(places),0,30)
//	Variable factor = places>0 ? 10^places : 1
//	String str=""
//	Variable i,N=numpnts(vec)
//	for (i=0;i<N;i+=1)
//		str += num2str(round(vec[i]*factor)/factor)+","
//	endfor
//	str = SelectString(strlen(str),"",str[0,strlen(str)-2])	// remove trailing comma
//	return str
//End
//
// read a vector from HDF file and return name of the wave containing the result
Static Function/T getWaveHDF5dataNum(f,dataName, [wName,Nww])
	Variable f
	String dataName
	String wName
	Variable Nww						// optional desired number of points in vector (maybe skip first SKIP_FIRST_N)
	if (ParamIsDefault(wName))
		wName = UniqueName("temp",1,0)
	endif
	Nww = ParamIsDefault(Nww) ? NaN : Nww
	Nww = Nww>SKIP_FIRST_N ? Nww : NaN

	HDF5LoadData/IGOR=-1/N=$wName/O/TYPE=2/Q/Z f,dataName
	if (V_flag)
		// printf "looking for '%s',    got V_Flag=%g\r",dataName,V_Flag
		return ""
	endif

	Wave ww = $wName
	Variable skip = (numpnts(ww)==(Nww+SKIP_FIRST_N) && Nww>0) ? SKIP_FIRST_N : 0
	if (skip>0)
		ww = ww[p+SKIP_FIRST_N]
		Redimension/N=(Nww) ww
	endif

	return GetWavesDataFolder($wName,2)
End
//
//Static Function/T getWaveHDF5dataNum(f,dataName, [wName])
//	Variable f
//	String dataName
//	String wName
//	if (ParamIsDefault(wName))
//		wName = UniqueName("temp",1,0)
//	endif
//
//	HDF5LoadData/IGOR=-1/N=$wName/O/TYPE=2/Q/Z f,dataName
//	if (V_flag)
//		// printf "looking for '%s',    got V_Flag=%g\r",dataName,V_Flag
//		return ""
//	endif
//
//	return GetWavesDataFolder($wName,2)
//End
//
Static Function/T get1HDF5AttrStr(f,dataName,attrName)		// get an attribute string
	Variable f
	String dataName, attrName

	String wName=UniqueName("temp",1,0)
	HDF5LoadData/A=attrName/N=$wName/O/Q/TYPE=1/Z f, dataName
	if (V_Flag)
		return ""
	endif
	Wave/T wstring = $wName
	String str = wstring[0]
	KillWaves $wName
	return str
End



// This routine is DEPRECATED.  It is OLD, and makes bad assumptions.
Static Function/T GuessDetectorID(wnote,[fast])
	String wnote
	Variable fast
	if (strlen(wnote)<1)
		return ""
	endif
	fast = ParamIsDefault(fast) ? 0 : fast

	String detectorID=""
	Variable endxy = max(NumberByKey("endx",wnote,"="),NumberByKey("endy",wnote,"="))
	String fileName = StringByKey("fileName",wnote,"=")
	if (stringmatch(fileName,"Orange") || endxy>2046)
		detectorID = "PE1621, 723-3335"				// Orange detector
	elseif (stringmatch(fileName,"*Yellow*") && endxy>1022)
		detectorID = "PE0820, 763-1807"				// Yellow detector
	elseif (stringmatch(fileName,"*Purple*") && endxy>1022)
		detectorID = "PE0820, 763-1850"				// Purple detector
	elseif (endxy>1022 && fast)
		detectorID = "PE0820, 763-1807"				// Yellow detector, fast assumes 1K x 1K is always Yellow
	elseif (endxy>1022)
		detectorID = "_none_"
		String knowDetectorIDs = "PE1621 723-3335 (Orange);PE0820 763-1807 (Yellow);PE0820 763-1850 (Purple);"
		knowDetectorIDs += "PE0822 883-4841 (Yellow);PE0822 883-4843 (Purple);"
		Prompt detectorID,"detector ID",popup, "_none_;"+knowDetectorIDs
		DoPrompt "detector ID",detectorID
		if (V_flag)
			detectorID = ""
		endif
		detectorID = ReplaceString("_none_",detectorID,"")
		Variable i = strsearch(detectorID," (",0)
		if (i>0)
			detectorID = detectorID[0,i-1]
		endif
	else
		DoAlert 0, "Unable to set detectorID for image "+fileName
	endif
	return detectorID
End
//
//Static Function/T GuessDetectorID(image,[fast])
//	Wave image
//	Variable fast
//	if (!WaveExists(Image))
//		return ""
//	endif
//	fast = ParamIsDefault(fast) ? 0 : fast
//
//	String wnote=note(image), detectorID=""
//	Variable endxy = max(NumberByKey("endx",wnote,"="),NumberByKey("endy",wnote,"="))
//	String fileName = StringByKey("fileName",wnote,"=")
//	if (stringmatch(fileName,"Orange") || endxy>2046)
//		detectorID = "PE1621, 723-3335"				// Orange detector
//	elseif (stringmatch(fileName,"*Yellow*") && endxy>1022)
//		detectorID = "PE0820, 763-1807"				// Yellow detector
//	elseif (stringmatch(fileName,"*Purple*") && endxy>1022)
//		detectorID = "PE0820, 763-1850"				// Purple detector
//	elseif (endxy>1022 && fast)
//		detectorID = "PE0820, 763-1807"				// Yellow detector, fast assumes 1K x 1K is always Yellow
//	elseif (endxy>1022)
//		detectorID = "_none_"
//		Prompt detectorID,"detector ID",popup, "_none_;PE1621, 723-3335 (Orange);PE0820, 763-1807 (Yellow);PE0820, 763-1850 (Purple)"
//		DoPrompt "detector ID",detectorID
//		if (V_flag)
//			detectorID = ""
//		endif
//		detectorID = ReplaceString("_none_",detectorID,"")
//		Variable i = strsearch(detectorID," (",0)
//		if (i>0)
//			detectorID = detectorID[0,i-1]
//		endif
//	else
//		DoAlert 0, "Unable to set detectorID for image "+NameOfWave(image)
//	endif
//	return detectorID
//End



Function/T HDFSaveImage(fName,image)	// saves an image to an HDF5 file
	String fName							// fully qualified name of file to open (will not prompt)
	Wave image
	if (!WaveExists(image))
		print "ERROR -- HDFSaveImage(), the image wave does not exist"
		return ""
	elseif (WaveDims(image)!=2)
		printf "ERROR -- HDFSaveImage(), only writes 2-d waves, '%s' is %d-d\r",WaveDims(image),NameOfWave(image)
		return ""
	elseif (strlen(fName)<1)
		printf "ERROR -- HDFSaveImage(), No file name passed\r"
		return ""
	endif

	Variable f
	PathInfo home
	if (V_flag)
		HDF5CreateFile/P=home/O/Z f as fName
	else
		HDF5CreateFile/O/Z f as fName		//rxadd
	endif
	if (V_flag != 0)
		printf "ERROR -- HDF5CreateFile failed saving image '%s' to file '%s'\r",NameOfWave(image),fName
		return ""
	endif
	Variable entry1, entry1_data, entry1_wire
	HDF5CreateGroup/Z f ,"entry1",entry1
	if (V_Flag)
		print "ERROR -- HDF5CreateGroup failed to make /entry1"
		HDF5CloseFile/Z f
		return ""
	endif
	HDF5CreateGroup/Z entry1 ,"data",entry1_data
	if (V_Flag)
		print "ERROR -- HDF5CreateGroup failed to make /entry1/data"
		HDF5CloseFile/Z f
		return ""
	endif
	MatrixOP/FREE it = image^t
	HDF5SaveData/IGOR=0/O/Z it, entry1_data, "data"
	if (V_flag)
		print "ERROR -- HDF5SaveData failed to save image '%s' into file '%s'\r",NameOfWave(image),fName
		HDF5CloseFile/Z f
		return ""
	endif

	// other things to store in the file
	Make/N=1/FREE/D array
	Make/N=1/FREE/I iarray
	Make/N=1/T/FREE sarray
	String wnote = note(image)

	// top level attributes
	sarray = fName
	HDF5SaveData/A="file_name"/O/Z sarray, f, "/"
	if (V_Flag)
		print "ERROR -- HDF5SaveData failed to make attribute 'file_name'"
		HDF5CloseFile/Z f
		return ""
	endif

	Variable epoch=DateTime
	sarray = Secs2Date(epoch,-2)+"T"+Secs2Time(epoch,3)+Secs2Time(date2secs(-1,-1,-1),4)
	HDF5SaveData/A="file_time"/O/Z sarray, f, "/"
	if (V_Flag)
		print "ERROR -- HDF5SaveData failed to make attribute 'file_time'"
		HDF5CloseFile/Z f
		return ""
	endif

	sarray = "Igor, HDF5images.ipf, v0.36"
	HDF5SaveData/A="creator"/O/Z sarray, f, "/"
	if (V_Flag)
		print "ERROR -- HDF5SaveData failed to make attribute 'creator'"
		HDF5CloseFile/Z f
		return ""
	endif

	// group entry1
	String title = StringByKey("title",wnote,"=")
	Variable scanNum = NumberByKey("scanNum",wnote,"=")
	Variable X2=NumberByKey("X2",wnote,"="), Y2=NumberByKey("Y2",wnote,"="), Z2=NumberByKey("Z2",wnote,"=")
	Variable wirebaseX=NumberByKey("wirebaseX",wnote,"="), wirebaseY=NumberByKey("wirebaseY",wnote,"="), wirebaseZ=NumberByKey("wirebaseZ",wnote,"=")
	Variable AerotechH=NumberByKey("AerotechH",wnote,"=")
	if (numtype(scanNum)==0)
		iarray = scanNum
		HDF5SaveData/IGOR=0/O/Z iarray, entry1, "scanNum"
	endif
	if (numtype(X2+Y2+Z2)==0)
		HDF5CreateGroup/Z entry1 ,"wire",entry1_wire
		if (V_Flag)
			print "ERROR -- HDF5CreateGroup failed to make /entry1/wire"
			HDF5CloseFile/Z f
			return ""
		endif
		array = X2
		HDF5SaveData/IGOR=0/O/Z array, entry1_wire, "wireX"
		array = Y2
		HDF5SaveData/IGOR=0/O/Z array, entry1_wire, "wireY"
		array = Z2
		HDF5SaveData/IGOR=0/O/Z array, entry1_wire, "wireZ"
		if (numtype(wirebaseX+wirebaseY+wirebaseZ)==0)
			array = wirebaseX
			HDF5SaveData/IGOR=0/O/Z array, entry1_wire, "wirebaseX"
			array = wirebaseY
			HDF5SaveData/IGOR=0/O/Z array, entry1_wire, "wirebaseY"
			array = wirebaseZ
			HDF5SaveData/IGOR=0/O/Z array, entry1_wire, "wirebaseZ"
		endif
		if (numtype(AerotechH)==0)
			array = AerotechH
			HDF5SaveData/IGOR=0/O/Z array, entry1_wire, "AerotechH"
		endif
	endif
	if (strlen(title))
		sarray = title
		HDF5SaveData/IGOR=0/O/Z sarray, entry1, "title"
	endif

	// group entry1/user
	String userName = StringByKey("userName",wnote,"=")
	if (strlen(userName))
		Variable entry1_user
		HDF5CreateGroup/Z entry1 ,"user",entry1_user
		if (V_Flag)
			print "ERROR -- HDF5CreateGroup failed to make /entry1/user"
			HDF5CloseFile/Z f
			return fName
		endif
		sarray = userName
		HDF5SaveData/IGOR=0/O/Z sarray, entry1_user, "name"
	endif

	// group entry1/microDiffraction
	Variable BeamBad = NumberByKey("BeamBad",wnote,"=")
	Variable CCDshutter = NumberByKey("CCDshutter",wnote,"=")
	Variable HutchTemperature = NumberByKey("HutchTemperature",wnote,"=")
	Variable LightOn = NumberByKey("LightOn",wnote,"=")
	Variable makeMicroDiff = numtype(BeamBad)==0 || numtype(CCDshutter)==0 || numtype(HutchTemperature)==0 || numtype(LightOn)==0

	Variable ringCurrent = NumberByKey("ringCurrent",wnote,"=")
	Variable gap = NumberByKey("undulatorGap",wnote,"=")
	Variable taper = NumberByKey("undulatorTaper",wnote,"=")
	Variable topUp = NumberByKey("topUp",wnote,"=")
	Variable makeSource = numtype(ringCurrent)==0 || numtype(gap)==0 || numtype(taper)==0

	if (makeMicroDiff || makeSource)
		Variable entry1_microDiffraction
		HDF5CreateGroup/Z entry1 ,"microDiffraction",entry1_microDiffraction
		if (V_Flag)
			print "ERROR -- HDF5CreateGroup failed to make /entry1/microDiffraction"
			HDF5CloseFile/Z f
			return fName
		endif
		if (numtype(BeamBad)==0)
			iarray = BeamBad
			HDF5SaveData/IGOR=0/O/Z iarray, entry1_microDiffraction, "BeamBad"
		endif
		if (numtype(CCDshutter)==0)
			iarray = CCDshutter
			HDF5SaveData/IGOR=0/O/Z iarray, entry1_microDiffraction, "CCDshutter"
		endif
		if (numtype(HutchTemperature)==0)
			array = HutchTemperature
			HDF5SaveData/IGOR=0/O/Z array, entry1_microDiffraction, "HutchTemperature"
		endif
		if (numtype(LightOn)==0)
			iarray = LightOn
			HDF5SaveData/IGOR=0/O/Z iarray, entry1_microDiffraction, "LightOn"
		endif

		// group entry1/microDiffraction/source
		if (makeSource)
			Variable entry1_microDiffraction_source
			HDF5CreateGroup/Z entry1_microDiffraction ,"source",entry1_microDiffraction_source
			if (V_Flag)
				print "ERROR -- HDF5CreateGroup failed to make /entry1/microDiffraction/source"
				HDF5CloseFile/Z f
				return fName
			endif
			if (numtype(ringCurrent)==0)
				array = ringCurrent
				HDF5SaveData/IGOR=0/O/Z array, entry1_microDiffraction_source, "current"
			endif
			if (numtype(gap)==0)
				array = gap
				HDF5SaveData/IGOR=0/O/Z array, entry1_microDiffraction_source, "gap"
			endif
			if (numtype(taper)==0)
				array = taper
				HDF5SaveData/IGOR=0/O/Z array, entry1_microDiffraction_source, "taper"
			endif
			if (numtype(topUp)==0)
				iarray = topUp
				HDF5SaveData/IGOR=0/O/Z iarray, entry1_microDiffraction_source, "topUp"
			endif
		endif

	endif

	// group entry1/monitor
	Variable I0 = NumberByKey("I0",wnote,"=")
	Variable Istart = NumberByKey("Istart",wnote,"=")
	Variable Ifinal = NumberByKey("Ifinal",wnote,"=")
	Variable I0_calc = NumberByKey("I0_calc",wnote,"=")
	Variable Ifinal_calc = NumberByKey("Ifinal_calc",wnote,"=")
	Variable Istart_calc = NumberByKey("Istart_calc",wnote,"=")
	Variable ScalerClockFreq = NumberByKey("ScalerClockFreq",wnote,"=")
	Variable ScalerClock_calc = NumberByKey("ScalerClock_calc",wnote,"=")
	Variable ScalerCountTime = NumberByKey("ScalerCountTime",wnote,"=")
	if (numtype(I0)==0 || numtype(Istart)==0 || numtype(Ifinal)==0)		// if any are OK, then write
		Variable entry1_monitor
		HDF5CreateGroup/Z entry1 ,"monitor",entry1_monitor
		if (V_Flag)
			print "ERROR -- HDF5CreateGroup failed to make /entry1/monitor"
			HDF5CloseFile/Z f
			return fName
		endif
		if (numtype(I0)==0)
			iarray = I0
			HDF5SaveData/IGOR=0/O/Z iarray, entry1_monitor, "I0"
		endif
		if (numtype(Istart)==0)
			iarray = Istart
			HDF5SaveData/IGOR=0/O/Z iarray, entry1_monitor, "I_start"
		endif
		if (numtype(Ifinal)==0)
			iarray = Ifinal
			HDF5SaveData/IGOR=0/O/Z iarray, entry1_monitor, "I_final"
		endif
		if (numtype(I0_calc)==0)
			array = I0_calc
			HDF5SaveData/IGOR=0/O/Z array, entry1_monitor, "I0_calc"
		endif
		if (numtype(Ifinal_calc)==0)
			array = Ifinal_calc
			HDF5SaveData/IGOR=0/O/Z array, entry1_monitor, "I_final_calc"
		endif
		if (numtype(Istart_calc)==0)
			array = Istart_calc
			HDF5SaveData/IGOR=0/O/Z array, entry1_monitor, "I_start_calc"
		endif
		if (numtype(ScalerClockFreq)==0)
			iarray = ScalerClockFreq
			HDF5SaveData/IGOR=0/O/Z iarray, entry1_monitor, "ScalerClockFreq"
		endif
		if (numtype(ScalerClock_calc)==0)
			iarray = ScalerClock_calc
			HDF5SaveData/IGOR=0/O/Z iarray, entry1_monitor, "ScalerClock_calc"
		endif
		if (numtype(ScalerCountTime)==0)
			iarray = ScalerCountTime
			HDF5SaveData/IGOR=0/O/Z iarray, entry1_monitor, "ScalerCountTime"
		endif
	endif

	// group entry1/sample
	Variable keV = NumberByKey("keV",wnote,"=")
	Variable X1 = NumberByKey("X1",wnote,"=")
	Variable Y1 = NumberByKey("Y1",wnote,"=")
	Variable Z1 = NumberByKey("Z1",wnote,"=")
	Variable sampleDistance = NumberByKey("sampleDistance",wnote,"=")
	String sampleName = StringByKey("sampleName",wnote,"=")
	Variable doit = strlen(sampleName)>0
	doit += numtype(X1+Y1+Z1)==0 + numtype(keV)==0 + numtype(sampleDistance)==0
	if (doit)
		Variable entry1_sample
		HDF5CreateGroup/Z entry1 ,"sample",entry1_sample
		if (V_Flag)
			print "ERROR -- HDF5CreateGroup failed to make /entry1/sample"
			HDF5CloseFile/Z f
			return fName
		endif
		if (numtype(keV)==0)
			array = keV
			HDF5SaveData/IGOR=0/O/Z array, entry1_sample, "incident_energy"
		endif
		if (numtype(X1+Y1+Z1)==0)
			array = X1
			HDF5SaveData/IGOR=0/O/Z array, entry1_sample, "sampleX"
			array = Y1
			HDF5SaveData/IGOR=0/O/Z array, entry1_sample, "sampleY"
			array = Z1
			HDF5SaveData/IGOR=0/O/Z array, entry1_sample, "sampleZ"
		endif
		if (numtype(sampleDistance)==0)
			array = sampleDistance
			HDF5SaveData/IGOR=0/O/Z array, entry1_sample, "distance"
		endif
		if (strlen(sampleName))
			sarray = sampleName
			HDF5SaveData/IGOR=0/O/Z sarray, entry1_sample, "name"
		endif
	endif

	// group entry1/detector
	Variable startx = NumberByKey("startx",wnote,"=")
	Variable endx = NumberByKey("endx",wnote,"=")
	Variable groupx = NumberByKey("groupx",wnote,"=")
	Variable starty = NumberByKey("starty",wnote,"=")
	Variable endy = NumberByKey("endy",wnote,"=")
	Variable groupy = NumberByKey("groupy",wnote,"=")
	Variable Nx = NumberByKey("xDimDet",wnote,"=")
	Variable Ny = NumberByKey("yDimDet",wnote,"=")
	Variable exposure = NumberByKey("exposure",wnote,"=")
	Variable gain = str2num(StringByKey("detectorGain",wnote,"="))
	Make/FREE gains = {0.25,0.5,1,2,4,8}
	gain = BinarySearch(gains,gain)
	gain = gain<0 ? NaN : gain
	String ID = StringByKey("detectorID",wnote,"=")

	if (numtype(Nx+Ny))			// for all images, really want to save:   Nx,Ny,  startx,endx,groupsx,  starty,endy,groupsy,
		Nx = DimSize(image,0)
		Ny = DimSize(image,1)
	endif
	if (numtype(startx+endx+groupx))
		startx = 0
		endx = DimSize(image,0)-1
		groupx = 1
	endif
	if (numtype(starty+endy+groupy))
		starty = 0
		endy = DimSize(image,1)-1
		groupy = 1
	endif
	doit = numtype(gain)==0 || strlen(ID)>0
	doit += numtype(startx+endx+groupx+starty+endy+groupy)==0
	doit += numtype(Nx+Ny)==0 + numtype(exposure)==0
	if (doit)
		Variable entry1_detector
		HDF5CreateGroup/Z entry1 ,"detector",entry1_detector
		if (V_Flag)
			print "ERROR -- HDF5CreateGroup failed to make /entry1/detector"
			HDF5CloseFile/Z f
			return fName
		endif
		if (numtype(Nx+Ny)==0)
			iarray = Nx
			HDF5SaveData/IGOR=0/O/Z iarray, entry1_detector, "Nx"
			iarray = Ny
			HDF5SaveData/IGOR=0/O/Z iarray, entry1_detector, "Ny"
		endif
		if (numtype(startx+endx+groupx+starty+endy+groupy)==0)
			iarray = startx
			HDF5SaveData/IGOR=0/O/Z iarray, entry1_detector, "startx"
			iarray = endx
			HDF5SaveData/IGOR=0/O/Z iarray, entry1_detector, "endx"
			iarray = groupx
			HDF5SaveData/IGOR=0/O/Z iarray, entry1_detector, "binx"
			iarray = starty
			HDF5SaveData/IGOR=0/O/Z iarray, entry1_detector, "starty"
			iarray = endy
			HDF5SaveData/IGOR=0/O/Z iarray, entry1_detector, "endy"
			iarray = groupy
			HDF5SaveData/IGOR=0/O/Z iarray, entry1_detector, "biny"
		endif
		if (numtype(exposure)==0)
			array = exposure
			HDF5SaveData/IGOR=0/O/Z array, entry1_detector, "exposure"
		endif
		if (strlen(ID))
			sarray = ID
			HDF5SaveData/IGOR=0/O/Z sarray, entry1_detector, "ID"
		endif
		if (numtype(gain)==0)
			array = gain
			HDF5SaveData/IGOR=0/O/Z array, entry1_detector, "gain"
		endif
	endif

	HDF5CloseFile/Z f
	if (V_flag)
		print "ERROR -- HDF5CloseFile failed to close image '%s' in file '%s', file may be good(?)\r",NameOfWave(image),fName
	endif
	return fName
End


// Do not use Write1HDF5imageFile(), instead use GenericSaveImage() in ImageDisplayScaling.ipf
//	DEPRECATED        DEPRECATED        DEPRECATED        DEPRECATED        DEPRECATED        DEPRECATED        DEPRECATED
// This is an old style that is not used anywhere
Function/S Write1HDF5imageFile(image,fName)
	Wave image
	String fName												// fully qualified name of file to open
	if (!WaveExists(image))
		String imageName=""
		Prompt imageName,"image to save",popup,WaveList("*",";","DIMS:2")
		DoPrompt "image",imageName
		if (V_flag)
			return ""
		endif
		Wave image = $imageName
	endif

	Variable f
	if (strlen(ParseFilePath(3,fName,":",0,0))<1)				// call dialog if no file name passed
		Open /D/M="HDF5 file"/T=".h5 " f						// use /D to get full path name
		fName = S_filename
	endif
	if (strlen(fName)<1)										// no file name, quit
		return ""
	endif
	if (strsearch(fName, ".",0)<0)
		fName += ".h5"
	endif
	fName = TrimTrailingWhiteSpace(fName)					// remove trailing spaces
	// printf "Saving image '%s' [%d x %d] to file '%s'...\r",NameOfWave(image),DimSize(image,0),DimSize(image,1),fName

	return HDFSaveImage(fName,image)						// saves an image to an HDF5 file
End
//
//Function/S Write1HDF5imageFile(image,fName)
//	Wave image
//	String fName												// fully qualified name of file to open
//	if (!WaveExists(image))
//		String imageName=""
//		Prompt imageName,"image to save",popup,WaveList("*",";","DIMS:2")
//		DoPrompt "image",imageName
//		if (V_flag)
//			return ""
//		endif
//		Wave image = $imageName
//	endif
//
//	Variable f
//	if (strlen(ParseFilePath(3,fName,":",0,0))<1)				// call dialog if no file name passed
//		Open /D/M="HDF5 file"/T=".h5 " f						// use /D to get full path name
//		fName = S_filename
//	endif
//	if (strlen(fName)<1)										// no file name, quit
//		return ""
//	endif
//	if (strsearch(fName, ".",0)<0)
//		fName += ".h5"
//	endif
//	fName = TrimTrailingWhiteSpace(fName)					// remove trailing spaces
//	// printf "Saving image '%s' [%d x %d] to file '%s'...\r",NameOfWave(image),DimSize(image,0),DimSize(image,1),fName
//
//	PathInfo home
//	if (V_flag)
//		HDF5CreateFile/P=home/O/Z f as fname
//	else
//		HDF5CreateFile/O/Z f as fname		//rxadd
//	endif
//	if (V_flag != 0)
//		printf "ERROR -- HDF5CreateFile failed saving image '%s' to file '%s'\r",NameOfWave(image),fName
//		return ""
//	endif
//	Variable entry1, entry1_data, entry1_wire
//	HDF5CreateGroup/Z f ,"entry1",entry1
//	if (V_Flag)
//		print "ERROR -- HDF5CreateGroup failed to make /entry1"
//		HDF5CloseFile/Z f
//		return ""
//	endif
//	HDF5CreateGroup/Z entry1 ,"data",entry1_data
//	if (V_Flag)
//		print "ERROR -- HDF5CreateGroup failed to make /entry1/data"
//		HDF5CloseFile/Z f
//		return ""
//	endif
//	MatrixOP/FREE it = image^t
//	HDF5SaveData/IGOR=0/O/Z it, entry1_data, "data"
//	if (V_flag)
//		print "ERROR -- HDF5SaveData failed to save image '%s' into file '%s'\r",NameOfWave(image),fName
//		HDF5CloseFile/Z f
//		return ""
//	endif
//
//	// other things to store in the file
//	Make/N=1/FREE/D array
//	Make/N=1/FREE/I iarray
//	Make/N=1/T/FREE sarray
//	String wnote = note(image)
//
//	// top level attributes
//	sarray = fName
//	HDF5SaveData/A="file_name"/O/Z sarray, f, "/"
//	if (V_Flag)
//		print "ERROR -- HDF5SaveData failed to make attribute 'file_name'"
//		HDF5CloseFile/Z f
//		return ""
//	endif
//
//	Variable epoch=DateTime
//	sarray = Secs2Date(epoch,-2)+"T"+Secs2Time(epoch,3)+Secs2Time(date2secs(-1,-1,-1),4)
//	HDF5SaveData/A="file_time"/O/Z sarray, f, "/"
//	if (V_Flag)
//		print "ERROR -- HDF5SaveData failed to make attribute 'file_time'"
//		HDF5CloseFile/Z f
//		return ""
//	endif
//
//	sarray = "Igor, HDF5images.ipf, v0.209"
//	HDF5SaveData/A="creator"/O/Z sarray, f, "/"
//	if (V_Flag)
//		print "ERROR -- HDF5SaveData failed to make attribute 'creator'"
//		HDF5CloseFile/Z f
//		return ""
//	endif
//
//	// group entry1
//	String title = StringByKey("title",wnote,"=")
//	Variable scanNum = NumberByKey("scanNum",wnote,"=")
//	Variable X2=NumberByKey("X2",wnote,"="), Y2=NumberByKey("Y2",wnote,"="), Z2=NumberByKey("Z2",wnote,"=")
//	Variable wirebaseX=NumberByKey("wirebaseX",wnote,"="), wirebaseY=NumberByKey("wirebaseY",wnote,"="), wirebaseZ=NumberByKey("wirebaseZ",wnote,"=")
//	Variable AerotechH=NumberByKey("AerotechH",wnote,"=")
//	if (numtype(scanNum)==0)
//		iarray = scanNum
//		HDF5SaveData/IGOR=0/O/Z iarray, entry1, "scanNum"
//	endif
//	if (numtype(X2+Y2+Z2)==0)
//		HDF5CreateGroup/Z entry1 ,"wire",entry1_wire
//		if (V_Flag)
//			print "ERROR -- HDF5CreateGroup failed to make /entry1/wire"
//			HDF5CloseFile/Z f
//			return ""
//		endif
//		array = X2
//		HDF5SaveData/IGOR=0/O/Z array, entry1_wire, "wireX"
//		array = Y2
//		HDF5SaveData/IGOR=0/O/Z array, entry1_wire, "wireY"
//		array = Z2
//		HDF5SaveData/IGOR=0/O/Z array, entry1_wire, "wireZ"
//		if (numtype(wirebaseX+wirebaseY+wirebaseZ)==0)
//			array = wirebaseX
//			HDF5SaveData/IGOR=0/O/Z array, entry1_wire, "wirebaseX"
//			array = wirebaseY
//			HDF5SaveData/IGOR=0/O/Z array, entry1_wire, "wirebaseY"
//			array = wirebaseZ
//			HDF5SaveData/IGOR=0/O/Z array, entry1_wire, "wirebaseZ"
//		endif
//		if (numtype(AerotechH)==0)
//			array = AerotechH
//			HDF5SaveData/IGOR=0/O/Z array, entry1_wire, "AerotechH"
//		endif
//	endif
//	if (strlen(title))
//		sarray = title
//		HDF5SaveData/IGOR=0/O/Z sarray, entry1, "title"
//	endif
//
//	// group entry1/user
//	String userName = StringByKey("userName",wnote,"=")
//	if (strlen(userName))
//		Variable entry1_user
//		HDF5CreateGroup/Z entry1 ,"user",entry1_user
//		if (V_Flag)
//			print "ERROR -- HDF5CreateGroup failed to make /entry1/user"
//			HDF5CloseFile/Z f
//			return fName
//		endif
//		sarray = userName
//		HDF5SaveData/IGOR=0/O/Z sarray, entry1_user, "name"
//	endif
//
//	// group entry1/microDiffraction
//	Variable BeamBad = NumberByKey("BeamBad",wnote,"=")
//	Variable CCDshutter = NumberByKey("CCDshutter",wnote,"=")
//	Variable HutchTemperature = NumberByKey("HutchTemperature",wnote,"=")
//	Variable LightOn = NumberByKey("LightOn",wnote,"=")
//	Variable makeMicroDiff = numtype(BeamBad)==0 || numtype(CCDshutter)==0 || numtype(HutchTemperature)==0 || numtype(LightOn)==0
//
//	Variable ringCurrent = NumberByKey("ringCurrent",wnote,"=")
//	Variable gap = NumberByKey("undulatorGap",wnote,"=")
//	Variable taper = NumberByKey("undulatorTaper",wnote,"=")
//	Variable topUp = NumberByKey("topUp",wnote,"=")
//	Variable makeSource = numtype(ringCurrent)==0 || numtype(gap)==0 || numtype(taper)==0
//
//	if (makeMicroDiff || makeSource)
//		Variable entry1_microDiffraction
//		HDF5CreateGroup/Z entry1 ,"microDiffraction",entry1_microDiffraction
//		if (V_Flag)
//			print "ERROR -- HDF5CreateGroup failed to make /entry1/microDiffraction"
//			HDF5CloseFile/Z f
//			return fName
//		endif
//		if (numtype(BeamBad)==0)
//			iarray = BeamBad
//			HDF5SaveData/IGOR=0/O/Z iarray, entry1_microDiffraction, "BeamBad"
//		endif
//		if (numtype(CCDshutter)==0)
//			iarray = CCDshutter
//			HDF5SaveData/IGOR=0/O/Z iarray, entry1_microDiffraction, "CCDshutter"
//		endif
//		if (numtype(HutchTemperature)==0)
//			array = HutchTemperature
//			HDF5SaveData/IGOR=0/O/Z array, entry1_microDiffraction, "HutchTemperature"
//		endif
//		if (numtype(LightOn)==0)
//			iarray = LightOn
//			HDF5SaveData/IGOR=0/O/Z iarray, entry1_microDiffraction, "LightOn"
//		endif
//
//		// group entry1/microDiffraction/source
//		if (makeSource)
//			Variable entry1_microDiffraction_source
//			HDF5CreateGroup/Z entry1_microDiffraction ,"source",entry1_microDiffraction_source
//			if (V_Flag)
//				print "ERROR -- HDF5CreateGroup failed to make /entry1/microDiffraction/source"
//				HDF5CloseFile/Z f
//				return fName
//			endif
//			if (numtype(ringCurrent)==0)
//				array = ringCurrent
//				HDF5SaveData/IGOR=0/O/Z array, entry1_microDiffraction_source, "current"
//			endif
//			if (numtype(gap)==0)
//				array = gap
//				HDF5SaveData/IGOR=0/O/Z array, entry1_microDiffraction_source, "gap"
//			endif
//			if (numtype(taper)==0)
//				array = taper
//				HDF5SaveData/IGOR=0/O/Z array, entry1_microDiffraction_source, "taper"
//			endif
//			if (numtype(topUp)==0)
//				iarray = topUp
//				HDF5SaveData/IGOR=0/O/Z iarray, entry1_microDiffraction_source, "topUp"
//			endif
//		endif
//
//	endif
//
//	// group entry1/monitor
//	Variable I0 = NumberByKey("I0",wnote,"=")
//	Variable Istart = NumberByKey("Istart",wnote,"=")
//	Variable Ifinal = NumberByKey("Ifinal",wnote,"=")
//	Variable I0_calc = NumberByKey("I0_calc",wnote,"=")
//	Variable Ifinal_calc = NumberByKey("Ifinal_calc",wnote,"=")
//	Variable Istart_calc = NumberByKey("Istart_calc",wnote,"=")
//	Variable ScalerClockFreq = NumberByKey("ScalerClockFreq",wnote,"=")
//	Variable ScalerClock_calc = NumberByKey("ScalerClock_calc",wnote,"=")
//	Variable ScalerCountTime = NumberByKey("ScalerCountTime",wnote,"=")
//	if (numtype(I0)==0 || numtype(Istart)==0 || numtype(Ifinal)==0)		// if any are OK, then write
//		Variable entry1_monitor
//		HDF5CreateGroup/Z entry1 ,"monitor",entry1_monitor
//		if (V_Flag)
//			print "ERROR -- HDF5CreateGroup failed to make /entry1/monitor"
//			HDF5CloseFile/Z f
//			return fName
//		endif
//		if (numtype(I0)==0)
//			iarray = I0
//			HDF5SaveData/IGOR=0/O/Z iarray, entry1_monitor, "I0"
//		endif
//		if (numtype(Istart)==0)
//			iarray = Istart
//			HDF5SaveData/IGOR=0/O/Z iarray, entry1_monitor, "I_start"
//		endif
//		if (numtype(Ifinal)==0)
//			iarray = Ifinal
//			HDF5SaveData/IGOR=0/O/Z iarray, entry1_monitor, "I_final"
//		endif
//		if (numtype(I0_calc)==0)
//			array = I0_calc
//			HDF5SaveData/IGOR=0/O/Z array, entry1_monitor, "I0_calc"
//		endif
//		if (numtype(Ifinal_calc)==0)
//			array = Ifinal_calc
//			HDF5SaveData/IGOR=0/O/Z array, entry1_monitor, "I_final_calc"
//		endif
//		if (numtype(Istart_calc)==0)
//			array = Istart_calc
//			HDF5SaveData/IGOR=0/O/Z array, entry1_monitor, "I_start_calc"
//		endif
//		if (numtype(ScalerClockFreq)==0)
//			iarray = ScalerClockFreq
//			HDF5SaveData/IGOR=0/O/Z iarray, entry1_monitor, "ScalerClockFreq"
//		endif
//		if (numtype(ScalerClock_calc)==0)
//			iarray = ScalerClock_calc
//			HDF5SaveData/IGOR=0/O/Z iarray, entry1_monitor, "ScalerClock_calc"
//		endif
//		if (numtype(ScalerCountTime)==0)
//			iarray = ScalerCountTime
//			HDF5SaveData/IGOR=0/O/Z iarray, entry1_monitor, "ScalerCountTime"
//		endif
//	endif
//
//	// group entry1/sample
//	Variable keV = NumberByKey("keV",wnote,"=")
//	Variable X1 = NumberByKey("X1",wnote,"=")
//	Variable Y1 = NumberByKey("Y1",wnote,"=")
//	Variable Z1 = NumberByKey("Z1",wnote,"=")
//	Variable sampleDistance = NumberByKey("sampleDistance",wnote,"=")
//	String sampleName = StringByKey("sampleName",wnote,"=")
//	Variable doit = strlen(sampleName)>0
//	doit += numtype(X1+Y1+Z1)==0 + numtype(keV)==0 + numtype(sampleDistance)==0
//	if (doit)
//		Variable entry1_sample
//		HDF5CreateGroup/Z entry1 ,"sample",entry1_sample
//		if (V_Flag)
//			print "ERROR -- HDF5CreateGroup failed to make /entry1/sample"
//			HDF5CloseFile/Z f
//			return fName
//		endif
//		if (numtype(keV)==0)
//			array = keV
//			HDF5SaveData/IGOR=0/O/Z array, entry1_sample, "incident_energy"
//		endif
//		if (numtype(X1+Y1+Z1)==0)
//			array = X1
//			HDF5SaveData/IGOR=0/O/Z array, entry1_sample, "sampleX"
//			array = Y1
//			HDF5SaveData/IGOR=0/O/Z array, entry1_sample, "sampleY"
//			array = Z1
//			HDF5SaveData/IGOR=0/O/Z array, entry1_sample, "sampleZ"
//		endif
//		if (numtype(sampleDistance)==0)
//			array = sampleDistance
//			HDF5SaveData/IGOR=0/O/Z array, entry1_sample, "distance"
//		endif
//		if (strlen(sampleName))
//			sarray = sampleName
//			HDF5SaveData/IGOR=0/O/Z sarray, entry1_sample, "name"
//		endif
//	endif
//
//	// group entry1/detector
//	Variable startx = NumberByKey("startx",wnote,"=")
//	Variable endx = NumberByKey("endx",wnote,"=")
//	Variable groupx = NumberByKey("groupx",wnote,"=")
//	Variable starty = NumberByKey("starty",wnote,"=")
//	Variable endy = NumberByKey("endy",wnote,"=")
//	Variable groupy = NumberByKey("groupy",wnote,"=")
//	Variable Nx = NumberByKey("xDimDet",wnote,"=")
//	Variable Ny = NumberByKey("yDimDet",wnote,"=")
//	Variable exposure = NumberByKey("exposure",wnote,"=")
//	Variable gain = str2num(StringByKey("detectorGain",wnote,"="))
//	Make/FREE gains = {0.25,0.5,1,2,4,8}
//	gain = BinarySearch(gains,gain)
//	gain = gain<0 ? NaN : gain
//	String ID = StringByKey("detectorID",wnote,"=")
//
//	if (numtype(Nx+Ny))			// for all images, really want to save:   Nx,Ny,  startx,endx,groupsx,  starty,endy,groupsy,
//		Nx = DimSize(image,0)
//		Ny = DimSize(image,1)
//	endif
//	if (numtype(startx+endx+groupx))
//		startx = 0
//		endx = DimSize(image,0)-1
//		groupx = 1
//	endif
//	if (numtype(starty+endy+groupy))
//		starty = 0
//		endy = DimSize(image,1)-1
//		groupy = 1
//	endif
//	doit = numtype(gain)==0 || strlen(ID)>0
//	doit += numtype(startx+endx+groupx+starty+endy+groupy)==0
//	doit += numtype(Nx+Ny)==0 + numtype(exposure)==0
//	if (doit)
//		Variable entry1_detector
//		HDF5CreateGroup/Z entry1 ,"detector",entry1_detector
//		if (V_Flag)
//			print "ERROR -- HDF5CreateGroup failed to make /entry1/detector"
//			HDF5CloseFile/Z f
//			return fName
//		endif
//		if (numtype(Nx+Ny)==0)
//			iarray = Nx
//			HDF5SaveData/IGOR=0/O/Z iarray, entry1_detector, "Nx"
//			iarray = Ny
//			HDF5SaveData/IGOR=0/O/Z iarray, entry1_detector, "Ny"
//		endif
//		if (numtype(startx+endx+groupx+starty+endy+groupy)==0)
//			iarray = startx
//			HDF5SaveData/IGOR=0/O/Z iarray, entry1_detector, "startx"
//			iarray = endx
//			HDF5SaveData/IGOR=0/O/Z iarray, entry1_detector, "endx"
//			iarray = groupx
//			HDF5SaveData/IGOR=0/O/Z iarray, entry1_detector, "binx"
//			iarray = starty
//			HDF5SaveData/IGOR=0/O/Z iarray, entry1_detector, "starty"
//			iarray = endy
//			HDF5SaveData/IGOR=0/O/Z iarray, entry1_detector, "endy"
//			iarray = groupy
//			HDF5SaveData/IGOR=0/O/Z iarray, entry1_detector, "biny"
//		endif
//		if (numtype(exposure)==0)
//			array = exposure
//			HDF5SaveData/IGOR=0/O/Z array, entry1_detector, "exposure"
//		endif
//		if (strlen(ID))
//			sarray = ID
//			HDF5SaveData/IGOR=0/O/Z sarray, entry1_detector, "ID"
//		endif
//		if (numtype(gain)==0)
//			array = gain
//			HDF5SaveData/IGOR=0/O/Z array, entry1_detector, "gain"
//		endif
//	endif
//
//	HDF5CloseFile/Z f
//	if (V_flag)
//		print "ERROR -- HDF5CloseFile failed to close image '%s' in file '%s', file may be good(?)\r",NameOfWave(image),fName
//	endif
//	return fName
//End
//				End of
//	DEPRECATED        DEPRECATED        DEPRECATED        DEPRECATED        DEPRECATED        DEPRECATED        DEPRECATED




 Function fileTime2Epoch(fileTime,[UTC])		// find epoch from "file_time" in HDF5 image
	String fileTime								// string with the local time
	Variable UTC								// set to TRUE if you want seconds from 1/1/1904 UTC, otherwise just ignores time zones
	UTC = ParamIsDefault(UTC) ? 0 : UTC
	UTC = numtype(UTC) ? 0 : UTC
	if (strlen(fileTime)<1)
		return NaN
	endif

	Variable year,month,day,hr,mi,se,tzh=0,tzm=0
	Variable num=UTC ? 8 : 6					// number of numbers to interpret
	if (UTC)
		sscanf fileTime, "%4d-%02d-%02dT%02d:%02d:%02d%d:%d",year,month,day,hr,mi,se,tzh,tzm
	else
		sscanf fileTime, "%4d-%02d-%02dT%02d:%02d:%02d",year,month,day,hr,mi,se
	endif
	if (V_flag!=num)
		return NaN
	endif

	Variable seconds = date2secs(year,month,day)
	seconds += hr*3600 + mi*60 + se - (tzh*3600 + tzm*60)
	return seconds
End
