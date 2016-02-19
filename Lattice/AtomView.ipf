#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma version = 0.38
#pragma IgorVersion = 6.3
#pragma ModuleName=AtomView
#include "Elements", version>=1.77
//	#include "GizmoZoomTranslate", version>=1.37
#include "GizmoZoomTranslate", version>=2.00
#include "GizmoClip", version>=2.00
#include "GizmoMarkers", version>=2.00
#include "LatticeSym", version>=5.10
#requiredPackages "LatticeSym;"
#initFunctionName "Init_AtomViewLattice()"

// These Constant values can be OverRidden by adding the line in your Main Procedure window.  Don't change this file.
Constant AtomView_GrayBkg = 0.75		// Can OverRide with :OverRide Constant AtomView_GrayBkg=0.95
Constant AtomView_BondLineWidth = 5	// Can OverRide with :OverRide Constant AtomView_BondLineWidth=3
StrConstant AtomView_CellOutLineColor = "0.733333,0.733333,0.733333,1"
Constant AtomView_CellOutLineColorR = 0.733333
Constant AtomView_CellOutLineColorG = 0.733333
Constant AtomView_CellOutLineColorB = 0.733333
Constant AtomView_CellOutLineColorA = 1


#if (IgorVersion()<7)
StrConstant AtomView_BondColor = "0.4,0.4,0.4,1"
#else
Constant AtomView_BondColorR = 0.2
Constant AtomView_BondColorG = 0.2
Constant AtomView_BondColorB = 0.2
Constant AtomView_BondColorA = 0.8
Constant AtomView_BondDia = 0.008
#endif
Constant AtomView_UseCovalent = 1		// Can OverRide with :OverRide Constant AtomView_UseCovalent=1
												// turn this flag on to use covalent radius instead of atomic radius.
Constant AtomView_SphereQuality = 50	// Can OverRide with :OverRide Constant AtomView_SphereQuality=100
Constant AtomView_zero= 1e-10			// distances less than this (=1e-19 m) are considered zero
Constant AtomView_minBondLen = 0.050		// 0.050 nm = 50 pm, minimum possible distance between atoms (smallest known bond is 74 pm)
Constant AtomView_maxBondLen = 0.310		// 0.310 nm = 310 pm, maximum possible distance between atoms
Constant AtomView_UseBlend = -1		// Can OverRide with :OverRide Constant AtomView_UseBlend=1, 0=NoBlend, 1=Blend, -1=Auto

//Static Constant GizmoScaleSize_BASE=7.5
Static Constant GizmoScaleSize_BASE=3.75

Menu "Analysis"
	SubMenu "Lattice"
		"-"
		"Make Cells of Atoms...",MakeCellsOfLattice(NaN,NaN,NaN)
		MenuItemFolderWithClassExists("Print Bond Info","atomViewBonds","DIMS:2,MINCOLS:3,MAXCOLS:3"),AllUniqueBonds("")
		MenuItemFolderWithClassExists("Gizmo of Atoms","atomViewXYZ","DIMS:2,MINCOLS:3,MAXCOLS:3"),MakeAtomViewGizmo($"")
		MenuItemFolderWithClassExists("Atom Type at Cursor","atomViewBonds","DIMS:2,MINCOLS:3,MAXCOLS:3"),print AtomView#ShowAtomViewInfo(1)
	End
End


Static Function AfterFileOpenHook(refNum,file,pathName,type,creator,kind)
	Variable refNum, kind
	String file,pathName,type,creator
	if ((kind==1) || (kind==2))		// an experiment (packed or unpacked)
		Init_AtomViewLattice()
	endif
	return 0
End
Static Function IgorStartOrNewHook(IgorApplicationNameStr)
	String IgorApplicationNameStr
	Init_AtomViewLattice()
	return 0
End



Static Function AtomViewPopMenuProc(pa) : PopupMenuControl		// used in the LatticeSym Panel
	STRUCT WMPopupAction &pa
	if (pa.eventCode != 2)
		return 0
	endif

	if (strsearch(pa.popStr,"Make Cells of Atoms",0,2)>=0)
		//	printf "%sMakeCellsOfLattice(NaN,NaN,NaN)\r", BULLET
		MakeCellsOfLattice(NaN,NaN,NaN)
	elseif (strsearch(pa.popStr,"Bond Info",0,2)>=0)
		printf "%sMakeAtomViewGizmo($\"\")\r", BULLET
		AllUniqueBonds("", printIt=1)
	elseif (strsearch(pa.popStr,"Gizmo of Atoms",0,2)>=0)
		printf "%sMakeAtomViewGizmo($\"\")\r", BULLET
		MakeAtomViewGizmo($"")
	elseif (strsearch(pa.popStr,"Atom Type at Cursor",0,2)>=0)
		printf "%sShowAtomViewInfo(1)\r", BULLET
		print ShowAtomViewInfo(1)
	endif

	return 0
End



Static Function/T ShowAtomViewInfo(alert)
	Variable alert
	alert = !(!alert)

	String panelWin="GizmoScatterMarkerPanel"

	if (strlen(WinList(panelWin,"", "WIN:64 " ))<1)
		return ""
	endif

	ControlInfo/W=$panelWin waveSelectPopup
	if (abs(V_flag)!=3)
		return ""
	endif
	Wave xyz=$S_Value
	if (!WaveExists(xyz))
		return ""
	endif

	ControlInfo/W=$panelWin alongWave
	if (abs(V_flag)!=5)
		return ""
	endif
	Variable pnt = round(V_Value)
	if (numtype(pnt) || pnt<0 || pnt>=DimSize(xyz,0))
		return ""
	endif

	String wNote=note(xyz)
	Wave Zwave = $StringByKey("ZWave",wNote,"=")
	Wave/T AtomTypewave = $StringByKey("atomAtomTypeWave",wNote,"=")
	if (!WaveExists(Zwave) || !WaveExists(AtomTypewave))
		return ""
	endif
	Variable Z=round(Zwave[pnt])
	String type=AtomTypewave[pnt]
	Variable bondLenMax = NumberByKey("bondLenMax",wNote,"=")

	String str
	if (bondLenMax>0)
		sprintf str, "atom#%d,   %s [Z=%g]\rbond max = %g nm",pnt,type,Z,bondLenMax
	else
		sprintf str, "atom#%d,   %s [Z=%g]",pnt,type,Z
	endif

	if (alert)
		DoAlert 0,str
	endif
	return str
End



//  ============================================================================  //
//  ==================== Start Make Cell(s) info to Display ====================  //

Function/WAVE MakeCellsOfLattice(Na,Nb,Nc,[blen,GizmoScaleSize])
	Variable Na,Nb,Nc				// number of cells in each direction, default is 1
	Variable bLen					// maximum distance that gets a bond
	Variable GizmoScaleSize	// 1 gives nice pictures
	bLen = ParamIsDefault(bLen) ? NaN : bLen
	bLen = bLen>0 ? bLen : NaN
	GizmoScaleSize = ParamIsDefault(GizmoScaleSize) ? 1 : GizmoScaleSize
	GizmoScaleSize = numtype(GizmoScaleSize) || GizmoScaleSize<=0 ? 1 : GizmoScaleSize

	STRUCT crystalStructure xtal
	if (FillCrystalStructDefault(xtal))					//fill the lattice structure with test values
		DoAlert 0, "no crystal structure found"
		return $""
	endif

	if (numtype(Na+Nb+Nc) || Na<0.1 || Nb<0.1 || Nc<0.1)
		Na = numtype(Na) || Na<0.1 ? 1 : Na
		Nb = numtype(Nb) || Nb<0.1 ? 1 : Nb
		Nc = numtype(Nc) || Nc<0.1 ? 1 : Nc
		Prompt Na,"Number of Cells along a"
		Prompt Nb,"Number of Cells along b"
		Prompt Nc,"Number of Cells along c"
		Prompt blen,"Maximum bond length (nm)"
		Prompt GizmoScaleSize,"Gizmo Atom Scale Size (>0)"
		if (xtal.Nbonds >= 1)
			DoPrompt "Cells to Make",Na,Nb,Nc,GizmoScaleSize
		else
			DoPrompt "Cells to Make",Na,Nb,Nc,blen,GizmoScaleSize
		endif
		if (V_flag)
			return $""
		endif
		printf "%sMakeCellsOfLattice(%g,%g,%g", BULLET, Na,Nb,Nc
		if (numtype(blen)==0 || !ParamIsDefault(blen))
			printf ", blen=%g",blen
		endif
		if (!ParamIsDefault(GizmoScaleSize) || GizmoScaleSize!=1)
			printf ", GizmoScaleSize=%g",GizmoScaleSize
		endif
		printf ")\r"
	endif
	GizmoScaleSize = numtype(GizmoScaleSize) || GizmoScaleSize<=0 ? 1 : GizmoScaleSize
	if (!(Na>0 && Nb>0 && Nc>0 && numtype(Na+Nb+Nc)==0))
		return $""
	endif

	Wave xyz=MakeOneCellsAtoms(xtal,Na,Nb,Nc,blen=blen,GizmoScaleSize=GizmoScaleSize)

	if (strlen(GetRTStackInfo(2))<1)
		String wNote=note(xyz)
		String prefix=StringByKey("prefix",wNote,"=")
		String sourceFldr=StringByKey("sourceFldr",wNote,"=")
		String desc=StringByKey("desc",wNote,"=")
		printf "Created atoms for '%s'  in '%s'\r",desc,sourceFldr
		Variable Nbonds=NumberByKey("Nbonds",note(xyz),"=")
		Variable bondLenMax=NumberByKey("bondLenMax",wNote,"=")
		if (!(Nbonds>0))
			printf "  prefix = '%s',  NO bonds were created\r",prefix
		elseif (bondLenMax > 0)
			printf "  prefix = '%s',  using max bond length = %g nm,  made %g bonds\r",prefix,bondLenMax,Nbonds
		else
			printf "  prefix = '%s',  %g bonds (using given bond lengths)\r",prefix,Nbonds
		endif
	endif
	MakeAtomViewGizmo(xyz)
	return xyz
End


Function/WAVE MakeOneCellsAtoms(xtal,Na,Nb,Nc,[blen,GizmoScaleSize])
	STRUCT crystalStructure &xtal
	Variable Na,Nb,Nc				// number of cells in each direction, default is 1
	Variable bLen					// maximum distance that gets a bond
	Variable GizmoScaleSize	// 1 gives nice pictures
	Na = numtype(Na) || Na<=0 ? 1 : Na
	Nb = numtype(Nb) || Nb<=0 ? 1 : Nb
	Nc = numtype(Nc) || Nc<=0 ? 1 : Nc
	bLen = ParamIsDefault(bLen) ? NaN : bLen
	bLen = bLen>0 ? bLen : NaN
	GizmoScaleSize = ParamIsDefault(GizmoScaleSize) ? 1 : GizmoScaleSize
	GizmoScaleSize = numtype(GizmoScaleSize) || GizmoScaleSize<=0 ? 1 : GizmoScaleSize

	Variable xtalN = xtal.N
	if (xtalN<1)
		return $""
	endif

	Make/N=(3,3)/D/FREE direct
	direct[0][0] = xtal.a0;		direct[0][1] = xtal.b0;		direct[0][2] = xtal.c0
	direct[1][0] = xtal.a1;		direct[1][1] = xtal.b1;		direct[1][2] = xtal.c1
	direct[2][0] = xtal.a2;		direct[2][1] = xtal.b2;		direct[2][2] = xtal.c2
	Make/N=3/D/FREE Ns={ceil(Na),ceil(Nb),ceil(Nc)}
	MatrixOp/FREE/O maxSize0 = maxVal(abs((direct x Ns)))		// find biggest size
	Variable maxSize=maxSize0[0]

	String fldrSav= GetDataFolder(1)
	String atomFldr = "root:Packages:Lattices:"	// location of generic atom positions
	String desc = xtal.desc
	String prefix=CleanupName(desc[0,16],0)			// only keep the first 17 characters
	for (;strsearch(prefix,"__",0)>=0;)
		prefix = ReplaceString("__",prefix,"_")		// change multiple underscores to singles
	endfor
	Variable i=strsearch(prefix,"_",Inf,1)			// also, trim off a trailing underscore
	i = i>0 && i==strlen(prefix)-1 ? i-1 : strlen(prefix)-1
	prefix = prefix[0,i]

	//	Function/T AddEndingToWaveName(wName,waveNameEnd)
	String fldr = GetDataFolder(1)+prefix+"_AtomView"
	if (!stringmatch(GetDataFolder(0),prefix+"_AtomView"))		// not already in fldr
		NewDataFolder/O/S $fldr
	endif
	reMakeAtomXYZs(xtal)									// ensure that atom positions are current

	Variable Nmax=200
	Make/N=(Nmax,3)/O $(prefix+"_XYZ")/WAVE=xyz = NaN
	Make/N=(Nmax,3)/O $(prefix+"_Size")/WAVE=size = 1
	Make/N=(Nmax,4)/O $(prefix+"_RGBA")/WAVE=rgba = 0
	Make/N=(Nmax)/O $(prefix+"_Z")/WAVE=zw = NaN
	Make/N=(Nmax)/O $(prefix+"_Occupy")/WAVE=occw = NaN
	Make/N=(Nmax)/O/T $(prefix+"_Type")/WAVE=types = ""

	Make/N=4/FREE rgbai
	String wNote="", atomType
	Variable rad, m, Zatom, occupy, Natom, Ni, covRad, atomRad
	for (Natom=0,m=0; m<xtalN; m+=1)					// loop over each atom type
		Wave atomi=$(atomFldr + "atom"+num2istr(m))
		if (!WaveExists(atomi))
			break
		endif
		wNote = note(atomi)
		Zatom = NumberByKey("Zatom",wNote,"=")
		Zatom = Zatom>0 ? Zatom : 1
		covRad = Element_covRadius(Zatom)/10			// covalent radius in nm (not Angstrom)
		atomRad = Element_atomRadius(Zatom)/10
		rad = AtomView_UseCovalent && numtype(covRad)==0 ? covRad : atomRad
		rad = rad>0 ? rad : 2.0							// default radius, probably for trans-uranics
		occupy = NumberByKey("occupy",wNote,"=")
		occupy = occupy>=0 ? limit(occupy,0,1) : 1
		atomType = StringByKey("atomType",wNote,"=")
		Wave rgb = atomRGB(Zatom,1)						// returns the Jmol color for an atom
		rgbai[0,2] = rgb[p]
		rgbai[3] = occupy*0.5

		Wave xyzi = RepeatOneAtomSet(atomi,Na,Nb,Nc)
		if (!WaveExists(xyzi))
			continue
		endif

		MatrixOp/FREE/O xyzi = (direct x xyzi^t)^t	// convert fractional to real distances
		Ni = DimSize(xyzi,0)
		if (Natom+Ni>=Nmax)
			Nmax += 300
			Redimension/N=(Nmax,-1) xyz, size, rgba
			Redimension/N=(Nmax) zw, types, occw
		endif
		xyz[Natom,Natom+Ni-1][] = xyzi[p-Natom][q]
		types[Natom,Natom+Ni-1] = atomType
		size[Natom,Natom+Ni-1] = 2*rad/maxSize * (GizmoScaleSize*GizmoScaleSize_BASE)
		zw[Natom,Natom+Ni-1] = Zatom
		occw[Natom,Natom+Ni-1] = occupy
		rgba[Natom,Natom+Ni-1][] = rgbai[q]
		Natom += Ni
	endfor
	Redimension/N=(Natom,-1) xyz, size, rgba
	Redimension/N=(Natom) zw, types, occw
	xyz = abs(xyz)<AtomView_zero ? 0 : xyz

	if (!(blen>0))
		Wave bonds = MakeBondList_Given(prefix,xtal,xyz)
	endif
	if (!WaveExists(bonds))
//		if (!(blen>0))				// find max length to use for bonds
//			// blen = FindMinSeparation(xyz)*1.01
//			blen = FindMinSeparation(xyz)*1.05
//		endif
		// blen is set in MakeBondList_blen()
		Wave bonds = MakeBondList_blen(prefix,xyz,blen=blen)
	endif

	Wave cell = MakeCellOutline(prefix,xtal,Na=Na,Nb=Nb,Nc=Nc)
	Wave cell0 = MakeCellOutline(prefix,xtal,name=NameOfWave(cell)+"0")
//	if ((Na>1 || Nb >1 || Nc>1) && WaveExists(cell))
//		Wave cell0 = MakeCellOutline(prefix,xtal,name=NameOfWave(cell)+"0")
//	endif
	Wave corners = MakeGizmocubeCorners(xyz)

	SetDataFolder fldrSav

	wNote = "waveClass=atomViewXYZ;"
	wNote = ReplaceStringByKey("sourceFldr",wNote,GetWavesDataFolder(xyz,1),"=")
	wNote = ReplaceStringByKey("desc",wNote,desc,"=")
	wNote = ReplaceStringByKey("prefix",wNote,prefix,"=")
	wNote = ReplaceNumberByKey("Natom",wNote,Natom,"=")
	String str
	sprintf str,"%g %g %g",Na,Nb,Nc
	wNote = ReplaceStringByKey("Nabc",wNote,str,"=")

	Make/N=3/D/FREE a={xtal.a0,xtal.a1,xtal.a2}, b={xtal.b0,xtal.b1,xtal.b2}, c={xtal.c0,xtal.c1,xtal.c2}
	a = abs(a)<1e-12 ? 0 : a
	b = abs(b)<1e-12 ? 0 : b
	c = abs(c)<1e-12 ? 0 : c
	wNote = ReplaceStringByKey("aVec",wNote,vec2str(a,bare=1,sep=","),"=")
	wNote = ReplaceStringByKey("bVec",wNote,vec2str(b,bare=1,sep=","),"=")
	wNote = ReplaceStringByKey("cVec",wNote,vec2str(c,bare=1,sep=","),"=")

	if (blen>0)
		wNote = ReplaceNumberByKey("bondLenMax",wNote,blen,"=")
	endif
	wNote = ReplaceStringByKey("sizeWave",wNote,GetWavesDataFolder(size,2),"=")
	wNote = ReplaceStringByKey("rgbaWave",wNote,GetWavesDataFolder(rgba,2),"=")
	wNote = ReplaceStringByKey("ZWave",wNote,GetWavesDataFolder(zw,2),"=")
	wNote = ReplaceStringByKey("OccupyWave",wNote,GetWavesDataFolder(occw,2),"=")
	wNote = ReplaceStringByKey("atomAtomTypeWave",wNote,GetWavesDataFolder(types,2),"=")
	if (WaveExists(bonds))
		wNote = ReplaceStringByKey("bondsWave",wNote,GetWavesDataFolder(bonds,2),"=")
		Variable Nbonds=NumberByKey("Nbonds",note(bonds),"=")
		if (Nbonds>0)
			wNote = ReplaceNumberByKey("Nbonds",wNote,Nbonds,"=")
		endif
	endif
	if (WaveExists(corners))
		wNote = ReplaceStringByKey("cornersWave",wNote,GetWavesDataFolder(corners,2),"=")
	endif
	if (WaveExists(cell))
		wNote = ReplaceStringByKey("cellOutlineWave",wNote,GetWavesDataFolder(cell,2),"=")
	endif
	if (WaveExists(cell0))
		wNote = ReplaceStringByKey("cellOutlineWave0",wNote,GetWavesDataFolder(cell0,2),"=")
	endif
	Note/K xyz, wNote

	wNote = ReplaceStringByKey("source","waveClass=;",GetWavesDataFolder(xyz,2),"=")
	Note/K size, ReplaceStringByKey("waveClass",wNote,"atomViewSize","=")
	Note/K types, ReplaceStringByKey("waveClass",wNote,"atomViewAtomType","=")
	Note/K zw, ReplaceStringByKey("waveClass",wNote,"atomViewZ","=")
	Note/K occw, ReplaceStringByKey("waveClass",wNote,"atomViewOccupy","=")
	Note/K rgba, ReplaceStringByKey("waveClass",wNote,"atomViewRGBA","=")
	return xyz
End


Static Function/WAVE RepeatOneAtomSet(atomi,Na,Nb,Nc)
	Wave atomi						// these are RELATIVE coordinates, assumed in range [0,1)
	Variable Na,Nb,Nc				// number of cells in each direction, default is 1
	Na = numtype(Na) || Na<0.1 ? 1 : Na
	Nb = numtype(Nb) || Nb<0.1 ? 1 : Nb
	Nc = numtype(Nc) || Nc<0.1 ? 1 : Nc

	Variable Ni=DimSize(atomi,0)
	if (Ni<1)
		return $""
	endif

	Variable Nmax=300, N
	Make/N=(Nmax,3)/FREE/D vecs

	Make/N=3/D/FREE ahat={1,0,0}, bhat={0,1,0}, chat={0,0,1}, v
	Variable ia,ib,ic, i
	for (N=0,ic=0; ic<=ceil(Nc); ic+=1)
		for (ib=0; ib<=ceil(Nb); ib+=1)
			for (ia=0; ia<=ceil(Na); ia+=1)
				for (i=0;i<Ni;i+=1)	// loop over each atom in atomi
					v = atomi[i][p] + ia*ahat[p] + ib*bhat[p] + ic*chat[p]
					if (v[0]<=Na && v[1]<=Nb && v[2]<=Nc)	// keep this atom
						if (N>=Nmax)
							Nmax += 300
							Redimension/N=(Nmax,-1) vecs
						endif
						vecs[N][] = v[q]
						N += 1
					endif
				endfor
			endfor
		endfor
	endfor
	Redimension/N=(N,-1) vecs
	return vecs
End


Static Function/WAVE MakeCellOutline(prefix,xtal,[Na,Nb,Nc,name])	// makes a gizmo path wave that outlines conventional cell
	String prefix
	STRUCT crystalStructure &xtal
	Variable Na,Nb,Nc					// number of cells in each direction
	String name							// name of wave to make
	Na = ParamIsDefault(Na) ? 1 : Na
	Nb = ParamIsDefault(Nb) ? 1 : Nb
	Nc = ParamIsDefault(Nc) ? 1 : Nc
	Na = numtype(Na) || Na<1 ? 1 : ceil(Na)
	Nb = numtype(Nb) || Nb<1 ? 1 : ceil(Nb)
	Nc = numtype(Nc) || Nc<1 ? 1 : ceil(Nc)
	name = SelectString(ParamIsDefault(name),name,"")

	Make/N=3/FREE/D a,b,c			// the lattice vectors
	a = {xtal.a0, xtal.a1, xtal.a2}
	b = {xtal.b0, xtal.b1, xtal.b2}
	c = {xtal.c0, xtal.c1, xtal.c2}

	name = SelectString(strlen(name),prefix+"_CellOutline",name)
	Make/N=(23,3)/D/O $name/WAVE=cell = NaN

	cell[0][] = 0							// form the a-b base
	cell[1][] = Na*a[q]
	cell[2][] = Na*a[q] + Nb*b[q]
	cell[3][] = Nb*b[q]
	cell[4][] = 0
	cell[5][] = NaN

	cell[6][] = 0							// the 4 connectors between c=0 and c=1 planes
	cell[7][] = Nc*c[q]					// from 000
	cell[8][] = NaN

	cell[9][] = Na*a[q]					// from a
	cell[10][] = Na*a[q] + Nc*c[q]
	cell[11][] = NaN

	cell[12][] = Nb*b[q]				// from b
	cell[13][] = Nb*b[q] + Nc*c[q]
	cell[14][] = NaN

	cell[15][] = Na*a[q] +Nb*b[q]	// from a-b corner
	cell[16][] = Na*a[q] + Nb*b[q] + Nc*c[q]
	cell[17][] = NaN

	cell[18][] = Nc*c[q]				// form the a-b at c
	cell[19][] = Na*a[q] + Nc*c[q]
	cell[20][] = Na*a[q] + Nb*b[q] + Nc*c[q]
	cell[21][] = Nb*b[q] + Nc*c[q]
	cell[22][] = Nc*c[q]

	cell = abs(cell)<AtomView_zero ? 0 : cell

	String wNote = ReplaceStringByKey("prefix","waveClass=atomViewCellOutline;",prefix,"=")
	Note/K cell, wNote

	return cell
End


Static Function/WAVE MakeBondList_Given(prefix,xtal,xyz)	// This makes the bond list from given bond lengths
	String prefix
	STRUCT crystalStructure &xtal
	Wave xyz				// list of atom xyz positions

	if (!WaveExists(xyz))
		return $""
	elseif (DimSize(xyz,0)<2 || DimSize(xyz,1)<3)
		return $""
	elseif (xtal.Nbonds < 1)
		return $""
	endif
	Variable N=DimSize(xyz,0), Ntypes=(xtal.Nbonds)

	String name = GetWavesDataFolder(xyz,1)+prefix+"_Type"
	Wave/T types = $name
	if (!WaveExists(types))
		return $""
	endif

	name = GetWavesDataFolder(xyz,1)+prefix+"_Bonds"
	Variable Nmax=300, Nbonds=0
	Make/N=(Nmax,3)/O $name/WAVE=bonds = NaN
	name += "_Source"
	Make/N=(Nmax)/O $name/WAVE=bsource = NaN

	String wNote="waveClass=atomViewBonds;"
	wNote = ReplaceStringByKey("source",wNote,GetWavesDataFolder(xyz,2),"=")
	wNote = ReplaceStringByKey("prefix",wNote,prefix,"=")

	String type, n0,n1
	Make/N=3/D/FREE xyz0, dxyz

	Variable i,j, m, mswap, blenMax, Nb, len
	for (Nb=0,m=0; m<Ntypes; m+=1)			// loop over each defined bond type
		n0 = xtal.bond[m].label0
		n1 = xtal.bond[m].label1
		Make/N=(xtal.bond[m].N)/FREE/D blens
		blens = xtal.bond[m].len[p]
		blens = numtype(blens) ? -Inf : blens[p]
		blenMax = WaveMax(blens)*1.05		// the largest valid distance, use 5% over
		for (mswap=0;mswap<2;mswap+=1)
			for (j=0;j<(N-1);j+=1)
				xyz0 = xyz[j][p]
				type = types[j]
				if (!stringmatch(n0,type))
					continue							// first type does not match, skip
				endif
				for (i=j+1;i<N;i+=1)
					if (!stringmatch(n1,types[i]))
						continue						// second type does not match, skip
					endif
					dxyz = xyz0[p] - xyz[i][p]
					len = norm(dxyz)
					if (AtomView_zero<len && len<=blenMax)	// found a bond
						if ((Nb+3)>=Nmax)			// need more room
							Nmax += 300	
							Redimension/N=(Nmax,-1) bonds, bsource
						endif
						bonds[Nb+0][] = xyz0[q]
						bonds[Nb+1][] = xyz[i][q]
						bonds[Nb+2][] = NaN
						bsource[Nb+0] = j
						bsource[Nb+1] = i
						bsource[Nb+2] = NaN
						Nb += 3
						Nbonds += 1
					endif
				endfor
			endfor

			if (stringmatch(n0,n1))			// same atoms, don't repeat
				mswap +=1
			else
				n1 = xtal.bond[m].label0		// swap and try again
				n0 = xtal.bond[m].label1
			endif
		endfor
	endfor

	if (Nb<1)
		Redimension/N=(0,-1) bonds, bsource
		KillWaves/Z bonds, bsource
		return $""
	endif

	Redimension/N=(Nb-1,-1) bonds, bsource
	bonds = abs(bonds)<AtomView_zero ? 0 : bonds
	wNote = ReplaceNumberByKey("Nbonds",wNote,Nbonds,"=")
	Note/K bonds, wNote
	Note/K bsource, ReplaceStringByKey("waveClass",wNote,"atomViewBonds_Source","=")
	return bonds
End


Static Function/WAVE MakeBondList_blen(prefix,xyz,[bLen])	// This is a guess when no bonds are given
	String prefix
	Wave xyz				// list of atom xyz positions
	Variable bLen		// maximum distance that gets a bond
	bLen = ParamIsDefault(bLen) ? NaN : bLen
	bLen = bLen>0 ? bLen : NaN

	if (!WaveExists(xyz))
		return $""
	elseif (DimSize(xyz,0)<2 || DimSize(xyz,1)<3)
		return $""
	endif
	Variable N=DimSize(xyz,0)
	String name = GetWavesDataFolder(xyz,2)
	name = ReplaceString("_XYZ",name,"_Type")
	Wave/T types = $name

	Make/N=(N)/FREE eNeg=Element_electroneg(LatticeSym#ZfromLabel(types[p]))
	WaveStats/M=1/Q eNeg
	Variable isMetalic = (V_max-V_min) < 0.4

	if (!(blen>0))				// find max length to use for bonds
		if (isMetalic)
			blen = AtomView#FindMinSeparation(xyz)*1.05
		else
			blen = AtomView_maxBondLen
		else
		endif
	endif
	//	print "blen =",blen

	name = GetWavesDataFolder(xyz,1)+prefix+"_Bonds"
	Variable Nmax=300, Nbonds=0
	Make/N=(Nmax,3)/O $name/WAVE=bonds = NaN
	name += "_Source"
	Make/N=(Nmax)/O $name/WAVE=bsource = NaN
	String buniqueName = GetWavesDataFolder(xyz,1)+prefix+"_Bonds_Unique"

	String wNote="waveClass=atomViewBonds;"
	wNote = ReplaceStringByKey("source",wNote,GetWavesDataFolder(xyz,2),"=")
	wNote = ReplaceStringByKey("prefix",wNote,prefix,"=")
	wNote = ReplaceNumberByKey("bondLenMax",wNote,blen,"=")

	Make/N=3/D/FREE xyz0, dxyz
	Variable i,j, Nb, len, eNegj, eNegi
	for (Nb=0,j=0; j<(N-1); j+=1)
		xyz0 = xyz[j][p]
		eNegj = eNeg[j]
		for (i=j+1;i<N;i+=1)
			eNegi = eNeg[i]
			if (!isMetalic && abs(eNegi-eNegj)<0.5)	// a covalent or ionic bond needs deltaElector > 0.5
				continue
			endif
			dxyz = xyz0[p] - xyz[i][p]
			len = norm(dxyz)
			if (AtomView_minBondLen<len && len<=blen)	// found a bond
				if ((Nb+3)>=Nmax)				// need more room
					Nmax += 300	
					Redimension/N=(Nmax,-1) bonds, bsource
				endif
				bonds[Nb+0][] = xyz0[q]
				bonds[Nb+1][] = xyz[i][q]
				bonds[Nb+2][] = NaN
				bsource[Nb+0] = j
				bsource[Nb+1] = i
				bsource[Nb+2] = NaN
				Nb += 3
				Nbonds += 1
			endif
		endfor
	endfor

	if (Nb<1)
		Redimension/N=(0,-1) bonds, bsource
		KillWaves/Z bonds, bsource
		return $""
	endif
	// print "Nbonds =",Nbonds, "   ",SelectString(isMetalic,"Ionic","Metalic")

	Redimension/N=(Nb-1,-1) bonds, bsource
	bonds = abs(bonds)<AtomView_minBondLen ? 0 : bonds
	wNote = ReplaceNumberByKey("Nbonds",wNote,Nbonds,"=")
	Note/K bonds, wNote
	Note/K bsource, ReplaceStringByKey("waveClass",wNote,"atomViewBonds_Source","=")

	// Store the unique bonds found here
	Make/N=(Nb,3)/O/T $buniqueName/WAVE=bunique = ""
	Make/N=3/D/FREE bvec
	Variable dup,m,Nu=0		// Nu is number of unique bonds
	String b0,b1
	for (i=0,Nu=0; i<Nb; i+=3)
		b0 = types[bsource[i]]
		b1 = types[bsource[i+1]]

		for (m=0,dup=0; m<Nu && !dup; m+=1)	// search if b0,b1 already in bunique
			dup += stringmatch(bunique[m][0],b0) && stringmatch(bunique[m][1],b1)
			dup += stringmatch(bunique[m][0],b1) && stringmatch(bunique[m][1],b0)
		endfor
		if (!dup)										// no match, add this bond
			bunique[Nu][0] = b0
			bunique[Nu][1] = b1
			bvec = bonds[i][p] - bonds[i+1][p]
			bunique[Nu][2] = num2str(norm(bvec))
			Nu += 1
		endif
	endfor
	Redimension/N=(Nu,-1) bunique
	wNote = RemoveByKey("Nbonds",wNote,"=")
	wNote = ReplaceNumberByKey("NbondsUnique",wNote,Nu,"=")
	wNote = ReplaceStringByKey("waveClass",wNote,"atomViewBonds_Unique","=")
	Note/K bunique, wNote

	return bonds
End
//Static Function/WAVE MakeBondList_blen(prefix,xyz,[bLen])	// This is a guess when no bonds are given
//	String prefix
//	Wave xyz				// list of atom xyz positions
//	Variable bLen		// maximum distance that gets a bond
//	bLen = ParamIsDefault(bLen) ? NaN : bLen
//	bLen = bLen>0 ? bLen : NaN
//
//	if (!WaveExists(xyz))
//		return $""
//	elseif (DimSize(xyz,0)<2 || DimSize(xyz,1)<3)
//		return $""
//	endif
//	if (!(blen>0))				// find max length to use for bonds
//		blen = FindMinSeparation(xyz)*1.05
//	endif
//	Variable N=DimSize(xyz,0)
//
//	String name = GetWavesDataFolder(xyz,2)
//	name = ReplaceString("_XYZ",name,"_Type")
//	Wave/T types = $name
//
//	name = GetWavesDataFolder(xyz,1)+prefix+"_Bonds"
//	Variable Nmax=300, Nbonds=0
//	Make/N=(Nmax,3)/O $name/WAVE=bonds = NaN
//	name += "_Source"
//	Make/N=(Nmax)/O $name/WAVE=bsource = NaN
//	String buniqueName = GetWavesDataFolder(xyz,1)+prefix+"_Bonds_Unique"
//
//	String wNote="waveClass=atomViewBonds;"
//	wNote = ReplaceStringByKey("source",wNote,GetWavesDataFolder(xyz,2),"=")
//	wNote = ReplaceStringByKey("prefix",wNote,prefix,"=")
//	wNote = ReplaceNumberByKey("bondLenMax",wNote,blen,"=")
//
//	Make/N=3/D/FREE xyz0, dxyz
//	Variable i,j, Nb, len
//	for (Nb=0,j=0; j<(N-1); j+=1)
//		xyz0 = xyz[j][p]
//
//		for (i=j+1;i<N;i+=1)
//			dxyz = xyz0[p] - xyz[i][p]
//			len = norm(dxyz)
//			if (AtomView_minBondLen<len && len<=blen)	// found a bond
//				if ((Nb+3)>=Nmax)				// need more room
//					Nmax += 300	
//					Redimension/N=(Nmax,-1) bonds, bsource
//				endif
//				bonds[Nb+0][] = xyz0[q]
//				bonds[Nb+1][] = xyz[i][q]
//				bonds[Nb+2][] = NaN
//				bsource[Nb+0] = j
//				bsource[Nb+1] = i
//				bsource[Nb+2] = NaN
//				Nb += 3
//				Nbonds += 1
//			endif
//		endfor
//	endfor
//
//	if (Nb<1)
//		Redimension/N=(0,-1) bonds, bsource
//		KillWaves/Z bonds, bsource
//		return $""
//	endif
//
//	Redimension/N=(Nb-1,-1) bonds, bsource
//	bonds = abs(bonds)<AtomView_minBondLen ? 0 : bonds
//	wNote = ReplaceNumberByKey("Nbonds",wNote,Nbonds,"=")
//	Note/K bonds, wNote
//	Note/K bsource, ReplaceStringByKey("waveClass",wNote,"atomViewBonds_Source","=")
//
//	// Store the unique bonds found here
//	Make/N=(Nb,3)/O/T $buniqueName/WAVE=bunique = ""
//	Make/N=3/D/FREE bvec
//	Variable dup,m,Nu=0		// Nu is number of unique bonds
//	String b0,b1
//	for (i=0,Nu=0; i<Nb; i+=3)
//		b0 = types[bsource[i]]
//		b1 = types[bsource[i+1]]
//
//		for (m=0,dup=0; m<Nu && !dup; m+=1)	// search if b0,b1 already in bunique
//			dup += stringmatch(bunique[m][0],b0) && stringmatch(bunique[m][1],b1)
//			dup += stringmatch(bunique[m][0],b1) && stringmatch(bunique[m][1],b0)
//		endfor
//		if (!dup)										// no match, add this bond
//			bunique[Nu][0] = b0
//			bunique[Nu][1] = b1
//			bvec = bonds[i][p] - bonds[i+1][p]
//			bunique[Nu][2] = num2str(norm(bvec))
//			Nu += 1
//		endif
//	endfor
//	Redimension/N=(Nu,-1) bunique
//	wNote = RemoveByKey("Nbonds",wNote,"=")
//	wNote = ReplaceNumberByKey("NbondsUnique",wNote,Nu,"=")
//	wNote = ReplaceStringByKey("waveClass",wNote,"atomViewBonds_Unique","=")
//	Note/K bunique, wNote
//
//	return bonds
//End
//	Function test_MakeBondList_blen()
//		Wave xyz=Si_XYZ
//		Wave wb = MakeBondList_blen("Si",xyz)
//	
//		print "wave name =",NameOfWave(wb)
//	End

Function PrintComputedBondListXML(fldr)		// print the computed bond list from MakeBondList_blen()
	String fldr
	fldr = SelectString(strlen(fldr),"",":")+fldr

	String bList = WaveListClass("atomViewBonds_Unique","*","MINCOLS:3,MAXCOLS:3,TEXT:1",fldr=fldr)
	if (ItemsInList(bList)==1)
		Wave/T bunique = $StringFromList(0,bList)
	endif
	if (!WaveExists(bunique))
		print "Could not find *_Bonds_Unique wave"
		return 1
	endif
	Variable i,Nb=DimSize(bunique,0), blen
	for (i=0;i<Nb;i+=1)
		blen = str2num(bunique[i][2])
		printf"\t<bond_chemical unit=\"nm\" n0=\"%s\" n1=\"%s\">%.4g</bond_chemical>\r",bunique[i][0],bunique[i][1],blen
	endfor
	return 0
End




Static Function/WAVE atomRGB(Z,maxVal)		// returns the Jmol color for an atom
	Variable Z					// atomic number
	Variable maxVal			// max value for a color, the saturated value
	if (!(Z>=1 && Z<=109))
		return $""
	endif
	maxVal = numtype(maxVal)||maxVal<=0 ? 1 : maxVal
	String rgbs = "255,255,255;217,255,255;204,128,255;194,255,0;255,181,181;144,144,144;48,80,248;255,13,13;144,224,80;179,227,245;"
	rgbs += "171,92,242;138,255,0;191,166,166;240,200,160;255,128,0;255,255,48;31,240,31;128,209,227;143,64,212;61,255,0;230,230,230;"
	rgbs += "191,194,199;166,166,171;138,153,199;156,122,199;224,102,51;240,144,160;80,208,80;200,128,51;125,128,176;194,143,143;"
	rgbs += "102,143,143;189,128,227;255,161,0;166,41,41;92,184,209;112,46,176;0,255,0;148,255,255;148,224,224;115,194,201;84,181,181;"
	rgbs += "59,158,158;36,143,143;10,125,140;0,105,133;192,192,192;255,217,143;166,117,115;102,128,128;158,99,181;212,122,0;148,0,148;"
	rgbs += "66,158,176;87,23,143;0,201,000C;112,212,255;255,255,199;217,255,199;199,255,199;163,255,199;143,255,199;97,255,199;"
	rgbs += "69,255,199;48,255,199;31,255,199;0,255,156;0,230,117;0,212,82;0,191,56;0,171,36;77,194,255;77,166,255;33,148,214;38,125,171;"
	rgbs += "38,102,150;23,84,135;208,208,224;255,209,35;184,184,208;166,84,77;87,89,97;158,79,181;171,92,0;117,79,69;66,130,150;"
	rgbs += "66,0,102;0,125,0007;112,171,250;0,186,255;0,161,255;0,143,255;0,128,255;0,107,255;84,92,242;120,92,227;138,79,227;"
	rgbs += "161,54,212;179,31,212;179,31,186;179,13,166;189,13,135;199,0,102;204,0,89;209,0,79;217,0,69;224,0,56;230,0,46;235,0,38;"

	String rgb = StringFromList(Z-1,rgbs)
	Variable r = str2num(StringFromList(0,rgb,",")) / 255 * maxVal
	Variable g = str2num(StringFromList(1,rgb,",")) / 255 * maxVal
	Variable b = str2num(StringFromList(2,rgb,",")) / 255 * maxVal
	if (numtype(r+g+b))
		return $""
	endif
	Make/N=3/FREE w={r,g,b}
	return w
End


Static Function FindMinSeparation(xyz)	// find the closest distance between two atoms (but>0)
	Wave xyz				// list of atom xyz positions

	Variable N=DimSize(xyz,0)
	Make/N=3/D/FREE xyz0								// the test position
	Variable i,iN, dmin=Inf
	for (i=0;i<(N-1);i+=1)
		xyz0 = xyz[i][p]
		iN = N-i-1
		if (iN>1)											// many atom pairs to check, use MatrixOP
			Make/N=(iN,3)/FREE/D xyzi
			xyzi = xyz[p+i+1][q]
			MatrixOP/FREE/O dxyz = sqrt(sumRows(magSqr(xyzi - rowRepeat(xyz0,iN))))
			dxyz = dxyz<AtomView_minBondLen ? Inf : dxyz	// don't permit zero distances
			dmin = min(dmin,WaveMin(dxyz))
		else
			xyz0 -= xyz[i+1][p]							// only 1 atom pair to check
			Variable dlast = norm(xyz0)
			dmin = dlast > AtomView_minBondLen && dlast<dmin ? dlast : dmin
		endif
	endfor
	return dmin
End
//Function test_FindMinSeparation()
//	Make/N=(18,3)/FREE/D xyz
//	xyz[0][0]= {0,0.25,0,0.25,0.5,0.75,0.5,0.75,1,0,1,0,1,0,1,1,0.5,0.5}
//	xyz[0][1]= {0,0.25,0.5,0.75,0,0.25,0.5,0.75,0,1,1,0,0,1,1,0.5,1,0.5}
//	xyz[0][2]= {0,0.25,0.5,0.75,0.5,0.75,0,0.25,0,0,0,1,1,1,1,0.5,0.5,1}
//	xyz *= 0.54310206
//
//	Variable dmin = FindMinSeparation(xyz)
//	printf "dmin = %g\r",dmin
//End

//  ===================== End Make Cell(s) info to Display =====================  //
//  ============================================================================  //



//  ============================================================================  //
//  ======================== Start Get Info About Bonds ========================  //

Function/T AllUniqueBonds(fldrName,[printIt])
	// returns all the unique bonds found in this AtomView, 
	// structure of result is:
	//	Nbonds,bmin,bmax,;len1,id1,id1;len2,id2,id2;...
	String fldrName
	Variable printIt
	printIt = ParamIsDefault(printIt) ? NaN : printIt
	printIt = numtype(printIt) ? (strlen(GetRTStackInfo(2))==0) : !(!printIt)

	String dflist=FoldersWithWaveClass(fldrName,"atomViewBonds","*","DIMS:2,MINCOLS:3,MAXCOLS:3")
	if (ItemsInList(dflist))
		if (ItemsInList(dflist)<1)
			return ""
		elseif (ItemsInList(dflist)==1)
			fldrName = StringFromList(0,dflist)
			printIt = 1
		else
			Prompt fldrName,"Atom View Folder with Waves", popup, reverseList(dflist)
			DoPrompt "Atom View?", fldrName
			if (V_flag)
				return ""
			endif
			printIt =1
		endif
	endif
	fldrName += SelectString(StringMatch(fldrName,"*:"),":","")	// ensure fldrName ends with ":"
	fldrName = SelectString(strsearch(fldrName,"root:",0) && strsearch(fldrName,":",0),"",":")+fldrName	// add leading ":" unless "root:..."
	if (printIt)
		printf "AllUniqueBonds(\"%s\")\r",fldrName
	endif
	if (ItemsInList(FoldersWithWaveClass(fldrName,"atomViewBonds","*","")))
		return ""
	endif

	String key=ParseFilePath(0,fldrName,":",1,0)
	Variable i=strsearch(key,"_AtomView",-Inf)
	key = SelectString(i>1,key,key[0,i-1])
	String prefix = fldrName+key

//	Wave bonds = $(prefix+"_Bonds")
//	Wave bondSource = $(prefix+"_Bonds_Source")
	Wave bonds = $AddEndingToWaveName(prefix,"_Bonds")
	Wave bondSource = $AddEndingToWaveName(prefix,"_Bonds_Source")
	Wave/T type = $(prefix+"_Type")
	if (!WaveExists(bonds) || !WaveExists(bondSource) || !WaveExists(type))
		return ""
	endif

	Make/N=3/D/FREE diff
	Variable N=DimSize(bonds,0), m, blen
	Variable Nbonds=floor(N/3) + 1
	if (Nbonds<0 || Nbonds>1e9 || numtype(Nbonds))
		return ""
	endif
	Make/N=(Nbonds)/D/FREE bondLens=NaN
	Make/N=(Nbonds)/T/FREE bondTypes0="", bondTypes1=""
	for (i=0;i<N;i+=3)
		diff = bonds[i+1][p] - bonds[i][p]
		blen = norm(diff)
		m = i/3
		bondLens[m] = blen
		bondTypes0[m] = type[bondSource[i]]
		bondTypes1[m] = type[bondSource[i+1]]
	endfor

	Variable Ndups=0
	for (m=0;m<Nbonds;m+=1)					// remove duplicates
		for (i=m+1;i<Nbonds;i+=1)
			if (abs(bondLens[m]-bondLens[i])>1e-6)
				continue
			elseif (!StringMatch(bondTypes0[i],bondTypes0[m]))
				continue
			elseif (!StringMatch(bondTypes1[i],bondTypes1[m]))
				continue
			endif
			Ndups += numtype(bondLens[i]) ? 0 : 1
			bondLens[i] = NaN
			bondTypes0[i] = ""
			bondTypes1[i] = ""
		endfor
	endfor
	Make/N=(Nbonds)/FREE index
	MakeIndex bondLens, index
	IndexSort index, bondLens, bondTypes0, bondTypes1

	WaveStats/Q/M=1 bondLens
	Nbonds = V_npnts
	Variable bmin=V_min, bmax=V_max
	if (Nbonds<0 || Nbonds>1e9 || numtype(Nbonds))
		return ""
	endif
	Redimension/N=(Nbonds), bondLens, bondTypes0, bondTypes1
	if (printIt)
		printf "in \"%s\"\r",fldrName
		for (m=0;m<Nbonds;m+=1)
			printf "  bond [%s, %s] = %.4f nm\r",bondTypes0[m],bondTypes1[m],bondLens[m]
		endfor
		printf "range of bond lengths = [%g, %g]nm\r",bmin,bmax
	endif

	String str, out=""
	sprintf out, "%d,%g,%g;", Nbonds,bmin,bmax
	for (m=0;m<Nbonds;m+=1)
		sprintf str, "%g,%s,%s;", bondLens[m], bondTypes0[m],bondTypes1[m]
		out += str
	endfor

	return out
End

//  ========================= End Get Info About Bonds =========================  //
//  ============================================================================  //



//  ============================= Start Make Gizmo =============================  //
//  ============================================================================  //

Function/T MakeAtomViewGizmo(xyz,[showNames,scaleFactor,useBlend])	// returns name of Gizmo
	Wave xyz
	Variable showNames					// if true, show a,b,c labels on lattice vectors
	Variable scaleFactor				// scale up model in Gizmo Window
	Variable useBlend					// 0=no blend, 1=blend, -1=auto
	scaleFactor = ParamIsDefault(scaleFactor) ? 1.25 : scaleFactor
	scaleFactor = numtype(scaleFactor) || scaleFactor<=0 ? 1.25 : scaleFactor
	useBlend = ParamIsDefault(useBlend) || numtype(useBlend) ? AtomView_UseBlend : useBlend
	if(exists("NewGizmo")!=4)			// Do nothing if the Gizmo XOP is not available.
		DoAlert 0, "Gizmo XOP must be installed"
		return ""
	endif

	Variable ilist,i
	String name=""
	if (!WaveExists(xyz))
		String list = WaveListClass("atomViewXYZ","*","DIMS:2,MINCOLS:3,MAXCOLS:3")
		ilist = ItemsInList(list)
		if (ilist==1)
			name = StringFromList(0,list)
		elseif (ilist>1)
			Prompt name,"Wave with Atom Positions",popup,list
			DoPrompt "Atom Positions",name
			if (V_flag)
				return ""
			endif
		else									// maybe check a folder?
			String fldrName
			list=FoldersWithWaveClass("","atomViewXYZ","*","DIMS:2,MINCOLS:3,MAXCOLS:3")
			if (ItemsInList(list)<1)
				return ""
			elseif (ItemsInList(list)==1)
				fldrName = StringFromList(0,list)
			else
				Prompt fldrName,"Atom View Folder with Waves", popup, reverseList(list)
				DoPrompt "Atom View?", fldrName
				if (V_flag)
					return ""
				endif
			endif
			fldrName += SelectString(StringMatch(fldrName,"*:"),":","")	// ensure fldrName ends with ":"
			fldrName = SelectString(strsearch(fldrName,"root:",0) && strsearch(fldrName,":",0),"",":")+fldrName	// add leading ":" unless "root:..."
			String key=ParseFilePath(0,fldrName,":",1,0)
			i = strsearch(key,"_AtomView",-Inf)
			key = SelectString(i>1,key,key[0,i-1])		// remove trailing "_AtomView"
			name = fldrName+key+"_XYZ"						// add the "_XYZ"
		endif
		Wave xyz = $name
	endif
	if (!WaveExists(xyz))
		return ""
	endif

	String gizName = StringFromlist(0,FindGizmosWithWave(xyz))		// find the Gizmo which contains the specified wave
	if (strlen(gizName))
		DoWindow/F $gizName
		return gizName
	endif

	String wNote=note(xyz)
	String prefix=StringByKey("prefix",wNote,"=")
	if (strlen(prefix)<1)
		return ""
	endif
	gizName = CleanupName("Gizmo_"+prefix,0)

	Wave size = $StringByKey("sizeWave",wNote,"=")
	Wave rgba = $StringByKey("rgbaWave",wNote,"=")
	Wave Zwave = $StringByKey("ZWave",wNote,"=")
	Wave/T AtomTypewave = $StringByKey("atomAtomTypeWave",wNote,"=")
	Wave bonds = $StringByKey("bondsWave",wNote,"=")
	if (DimSize(bonds,0)<1)
		Wave bonds = $""
	endif
	Wave corners = $StringByKey("cornersWave",wNote,"=")
	Wave cell = $StringByKey("cellOutlineWave",wNote,"=")
	Wave cell0 = $StringByKey("cellOutlineWave0",wNote,"=")
	Variable bondLenMax = NumberByKey("bondLenMax",wNote,"=")
	String sourceFldr=StringByKey("sourceFldr",wNote,"=")
	String desc=StringByKey("desc",wNote,"=")
	Variable Na,Nb,Nc
	String title2="", title3=""
	sscanf StringByKey("Nabc",wNote,"="),"%g %g %g", Na,Nb,Nc
	if ((Na>0.1 || Nb>0.1 || Nc>0.1) && V_flag==3)
		sprintf title2,"%g x %g x %g cells",Na,Nb,Nc
	endif
	if (bondLenMax>0)
		sprintf title3,"max bond length = %g nm",bondLenMax
	endif
	Wave a = str2vec(StringByKey("aVec",wNote,"="),sep=",")
	Wave b = str2vec(StringByKey("bVec",wNote,"="),sep=",")
	Wave c = str2vec(StringByKey("cVec",wNote,"="),sep=",")

#if (IgorVersion()>=7)
	// for Igor7 prefer blend, but turn off blend when long labels are used
	useBlend = useBlend<0 && !StringMatch(AtomTypewave[0],"*001") ? 1 : useBlend
#endif
	if (useBlend < 0)			// auto was chosen, decide on blending
		Variable m, N=DimSize(xyz,0)
		Make/N=3/D/FREE xyz0
		// when useBlend == -1, then turn on blending when two atoms are at the same location
		for (m=0;m<N && useBlend<0;m+=1)				// search for two atoms at same position if useBlend is auto (-1)
			xyz0 = xyz[m][p]
			MatrixOP/FREE/O dxyz = greater(AtomView_zero,sumRows(magSqr(xyz-rowRepeat(xyz0,N))))
			useBlend = sum(dxyz)>=2 ? 1 : useBlend
		endfor
		useBlend = useBlend<0 ? 0 : !(!useBlend)	// if no duplicate atoms found, useBlend=0 (no blending)
	endif

	String str, objectList="", attributeList="", scaleBarGroup=""
#if (IgorVersion()<7)
	Execute "NewGizmo/N="+gizName+"/T=\""+gizName+"\" /W=(234,45,992,803)"
	Execute "ModifyGizmo startRecMacro"
#else
	NewGizmo/N=$gizName/W=(234,45,992,803)/T=gizName
	ModifyGizmo startRecMacro
#endif
	String titleGroup = AddGizmoTitleGroup("",desc,title2=title2,title3=title3,title4=sourceFldr)

	MatrixOP/FREE maxX = maxVal(col(xyz,0))
	MatrixOP/FREE maxY = maxVal(col(xyz,1))
	MatrixOP/FREE maxZ = maxVal(col(xyz,2))
	Variable maxLength = max(max(maxX[0],maxY[0]),maxZ[0])
	scaleBarGroup = AddScaleBarGroup("",maxLength,"nm",scaleFactor=scaleFactor)

	Variable fixedLight = 0
#if (IgorVersion()<7)
	if (strlen(titleGroup+scaleBarGroup) && fixedLight)
		Execute "AppendToGizmo light=Directional,name=lightFixed"				// this light is fixed (just like title)
		Execute "ModifyGizmo light=lightFixed property={ position,0,0,-1,0}"
		Execute "ModifyGizmo light=lightFixed property={ direction,0,0,-1}"
		Execute "ModifyGizmo light=lightFixed property={ ambient,0.45,0.45,0.45,1}"
		Execute "ModifyGizmo light=lightFixed property={ specular,1,1,1,1.0}"
		Execute "ModifyGizmo light=lightFixed property={ diffuse,0.45,0.45,0.45,1.0}"
	endif
#else
	if (strlen(scaleBarGroup) && fixedLight)
		AppendToGizmo light=Directional,name=lightFixed				// this light is fixed (just like title)
		ModifyGizmo modifyObject=lightFixed objectType=light, property={ position,0,0,-1,0}
		ModifyGizmo modifyObject=lightFixed objectType=light, property={ direction,0,0,-1}
		ModifyGizmo modifyObject=lightFixed objectType=light, property={ ambient,0.45,0.45,0.45,1}
		ModifyGizmo modifyObject=lightFixed objectType=light, property={ specular,1,1,1,1.0}
		ModifyGizmo modifyObject=lightFixed objectType=light, property={ diffuse,0.45,0.45,0.45,1.0}
	endif
#endif

	objectList += "lightMoving;"
#if (IgorVersion()<7)
	Execute "AppendToGizmo light=Directional,name=lightMoving"						// this light is spins with the atoms
	Execute "ModifyGizmo light=lightMoving property={ position,0,0,-1,0}"
	Execute "ModifyGizmo light=lightMoving property={ direction,0,0,-1}"
	Execute "ModifyGizmo light=lightMoving property={ specular,1,1,1,1}"
//	Execute "ModifyGizmo light=lightMoving property={ diffuse,0,0,0,1}"
	Execute "ModifyGizmo light=lightMoving property={ ambient,0.866667,0.866667,0.866667,1}"
	Execute "ModifyGizmo light=lightMoving property={ diffuse,0.866667,0.866667,0.866667,1}"
#else
	AppendToGizmo light=Directional,name=lightMoving							// this light is spins with the atoms
	ModifyGizmo modifyObject=lightMoving objectType=light, property={ position,0,0,-1,0}
	ModifyGizmo modifyObject=lightMoving objectType=light, property={ direction,0,0,-1}
	ModifyGizmo modifyObject=lightMoving objectType=light, property={ specular,1,1,1,1}
//	ModifyGizmo modifyObject=lightMoving objectType=light, property={ diffuse,0,0,0,1}
	ModifyGizmo modifyObject=lightMoving objectType=light, property={ ambient,0.866667,0.866667,0.866667,1}
	ModifyGizmo modifyObject=lightMoving objectType=light, property={ diffuse,0.866667,0.866667,0.866667,1}
#endif

	if (numpnts(a)==3 && numpnts(b)==3 && numpnts(c)==3)
		if (ParamIsDefault(showNames))
			objectList += AddRealLatticeAxesGroup("",a,b,c)+";"
		else
			objectList += AddRealLatticeAxesGroup("",a,b,c,showNames=showNames)+";"
		endif
	else
		// ************************* Group Object Start *******************
		objectList += "groupAxisCue;"
#if (IgorVersion()<7)
		Execute "AppendToGizmo group,name=groupAxisCue"
		Execute "ModifyGizmo currentGroupObject=\"groupAxisCue\""
		Execute "AppendToGizmo freeAxesCue={0,0,0,1},name=freeAxesCue0"
		Execute "AppendToGizmo attribute lineWidth=2, name=lineWidthCue"
		Execute "AppendToGizmo attribute specular={0,0,0,1,1032},name=specularOFF"
		Execute "AppendToGizmo attribute shininess={0,1032},name=shininessOFF"
		Execute "ModifyGizmo setDisplayList=0, attribute=lineWidthCue"
		Execute "ModifyGizmo setDisplayList=1, attribute=specularOFF"
		Execute "ModifyGizmo setDisplayList=1, attribute=shininessOFF"
		Execute "ModifyGizmo setDisplayList=2, object=freeAxesCue0"
		Execute "ModifyGizmo currentGroupObject=\"::\""
#else
		AppendToGizmo group,name=groupAxisCue
		ModifyGizmo currentGroupObject="groupAxisCue"
		AppendToGizmo freeAxesCue={0,0,0,1},name=freeAxesCue0
		AppendToGizmo attribute lineWidth=2, name=lineWidthCue
		AppendToGizmo attribute specular={0,0,0,1,1032},name=specularOFF
		AppendToGizmo attribute shininess={0,1032},name=shininessOFF
		ModifyGizmo setDisplayList=0, attribute=lineWidthCue
		ModifyGizmo setDisplayList=1, attribute=specularOFF
		ModifyGizmo setDisplayList=1, attribute=shininessOFF
		ModifyGizmo setDisplayList=2, object=freeAxesCue0
		ModifyGizmo currentGroupObject="::"
#endif
		// ************************* Group Object End *******************
	endif

#if (IgorVersion()<7)
	sprintf str,"AppendToGizmo sphere={0.1,%d,%d},name=generalAtom",AtomView_SphereQuality,AtomView_SphereQuality
	Execute str									// This is a generic atom, with shininess & some specular, no color or size
	Execute "AppendToGizmo attribute shininess={100,1032},name=shininessAtom0"
	Execute "AppendToGizmo attribute specular={1.0,1.0,0.8,1,1032},name=specularAtom0"
	Execute "ModifyGizmo setObjectAttribute={generalAtom,shininessAtom0}"
	Execute "ModifyGizmo setObjectAttribute={generalAtom,specularAtom0}"
	Execute "ModifyGizmo modifyObject=generalAtom,property={colorType,0}"
#else
	// This is a generic atom, with shininess & some specular, no color or size
	AppendToGizmo sphere={0.1,AtomView_SphereQuality,AtomView_SphereQuality},name=generalAtom
	AppendToGizmo attribute shininess={100,1032},name=shininessAtom0
	AppendToGizmo attribute specular={1.0,1.0,0.8,1,1032},name=specularAtom0
	ModifyGizmo setObjectAttribute={generalAtom,shininessAtom0}
	ModifyGizmo setObjectAttribute={generalAtom,specularAtom0}
	ModifyGizmo modifyObject=generalAtom, objectType=sphere, property={colorType,0}
#endif

#if (IgorVersion()<7)
	Execute "AppendToGizmo Scatter="+GetWavesDataFolder(xyz,2)+",name=atomViewAtoms"
	Execute "ModifyGizmo ModifyObject=atomViewAtoms property={ scatterColorType,1}"
	Execute "ModifyGizmo ModifyObject=atomViewAtoms property={ markerType,0}"
	Execute "ModifyGizmo ModifyObject=atomViewAtoms property={ sizeType,1}"
	Execute "ModifyGizmo ModifyObject=atomViewAtoms property={ rotationType,0}"
	Execute "ModifyGizmo ModifyObject=atomViewAtoms property={ size,0.2}"
	Execute "ModifyGizmo ModifyObject=atomViewAtoms property={ colorWave,"+GetWavesDataFolder(rgba,2)+"}"
	Execute "ModifyGizmo ModifyObject=atomViewAtoms property={ sizeWave,"+GetWavesDataFolder(size,2)+"}"
	Execute "ModifyGizmo ModifyObject=atomViewAtoms property={ Shape,7}"			// 7 means an object, set in next line
	Execute "ModifyGizmo ModifyObject=atomViewAtoms property={ objectName,generalAtom}"
	objectList += "atomViewAtoms;"
#else
//	if (useBlend)				// when using blend, you can also show the atom labels
	if (useBlend && Na*Nb*Nc <= 2)		// when using blend, you can also show the atom labels
		AppendToGizmo Scatter=$GetWavesDataFolder(xyz,2),name=atomViewAtomsLabels
		ModifyGizmo ModifyObject=atomViewAtomsLabels,objectType=scatter,property={ scatterColorType,0}
		ModifyGizmo ModifyObject=atomViewAtomsLabels,objectType=scatter,property={ markerType,0}
		ModifyGizmo ModifyObject=atomViewAtomsLabels,objectType=scatter,property={ sizeType,0}
		ModifyGizmo ModifyObject=atomViewAtomsLabels,objectType=scatter,property={ rotationType,0}
		ModifyGizmo ModifyObject=atomViewAtomsLabels,objectType=scatter,property={ Shape,8}
		ModifyGizmo ModifyObject=atomViewAtomsLabels,objectType=scatter,property={ size,1.5}
		ModifyGizmo ModifyObject=atomViewAtomsLabels,objectType=scatter,property={ color,0,0,0,1}
		ModifyGizmo ModifyObject=atomViewAtomsLabels,objectType=scatter,property={ TextWave,$GetWavesDataFolder(AtomTypewave,2)}
		objectList += "atomViewAtomsLabels;"
	endif
	AppendToGizmo Scatter=$GetWavesDataFolder(xyz,2),name=atomViewAtoms
	ModifyGizmo ModifyObject=atomViewAtoms objectType=scatter, property={ scatterColorType,1}
	ModifyGizmo ModifyObject=atomViewAtoms objectType=scatter, property={ markerType,0}
	ModifyGizmo ModifyObject=atomViewAtoms objectType=scatter, property={ sizeType,1}
	ModifyGizmo ModifyObject=atomViewAtoms objectType=scatter, property={ rotationType,0}
	ModifyGizmo ModifyObject=atomViewAtoms objectType=scatter, property={ size,0.2}
	ModifyGizmo ModifyObject=atomViewAtoms objectType=scatter, property={ colorWave, $GetWavesDataFolder(rgba,2)}
	ModifyGizmo ModifyObject=atomViewAtoms objectType=scatter, property={ sizeWave, $GetWavesDataFolder(size,2)}
	ModifyGizmo ModifyObject=atomViewAtoms objectType=scatter, property={ Shape,7}			// 7 means an object, set in next line
	ModifyGizmo ModifyObject=atomViewAtoms objectType=scatter, property={ objectName,generalAtom}
	objectList += "atomViewAtoms;"
#endif

	if (WaveExists(bonds))
		objectList += "atomViewBonds;"
		Variable Nbonds = NumberByKey("Nbonds",note(bonds),"=")
		Variable lineWidth = limit(2+3.9*exp(-0.0046*Nbonds),1,8)
		lineWidth = numtype(lineWidth) ? 3 : round(lineWidth*5)/5
#if (IgorVersion()<7)
		Execute "AppendToGizmo Path="+GetWavesDataFolder(bonds,2)+",name=atomViewBonds"
		Execute "ModifyGizmo ModifyObject=atomViewBonds property={ pathColorType,1}"
		Execute "ModifyGizmo ModifyObject=atomViewBonds property={ lineWidthType,1}"
		Execute "ModifyGizmo ModifyObject=atomViewBonds property={ lineWidth,"+num2str(lineWidth)+"}"
		Execute "ModifyGizmo ModifyObject=atomViewBonds property={ pathColor,"+AtomView_BondColor+"}"
		Execute "ModifyGizmo setObjectAttribute={atomViewBonds,specularBond0}"
		Execute "AppendToGizmo attribute specular={0.1,0.1,0.1,1,1032},name=specularBond0"
#else
		AppendToGizmo Path=$GetWavesDataFolder(bonds,2),name=atomViewBonds
		ModifyGizmo ModifyObject=atomViewBonds objectType=path, property={ pathColorType,1}
		ModifyGizmo ModifyObject=atomViewBonds objectType=path, property={ lineWidthType,1}
		ModifyGizmo ModifyObject=atomViewBonds objectType=path, property={ lineWidth,lineWidth}
		ModifyGizmo ModifyObject=atomViewBonds objectType=path, property={ pathColor,AtomView_BondColorR,AtomView_BondColorG,AtomView_BondColorB,AtomView_BondColorA}
		ModifyGizmo ModifyObject=atomViewBonds,objectType=path,property={ drawTube,1}
		ModifyGizmo ModifyObject=atomViewBonds,objectType=path,property={ fixedRadius,AtomView_BondDia}
			ModifyGizmo setObjectAttribute={atomViewBonds,specularBond0}
		AppendToGizmo attribute specular={0.1,0.1,0.1,1,1032},name=specularBond0
#endif
	endif

	if (WaveExists(cell))
		objectList += "cellOutline;"
#if (IgorVersion()<7)
		Execute "AppendToGizmo Path="+GetWavesDataFolder(cell,2)+",name=cellOutline"
		Execute "ModifyGizmo ModifyObject=cellOutline property={ pathColorType,1}"
		Execute "ModifyGizmo ModifyObject=cellOutline property={ lineWidthType,1}"
		Execute "ModifyGizmo ModifyObject=cellOutline property={ lineWidth,0.5}"
		Execute "ModifyGizmo ModifyObject=cellOutline property={ pathColor,"+AtomView_CellOutLineColor+"}"
#else
		AppendToGizmo Path=$GetWavesDataFolder(cell,2),name=cellOutline
		ModifyGizmo ModifyObject=cellOutline objectType=path, property={ pathColorType,1}
		ModifyGizmo ModifyObject=cellOutline objectType=path, property={ lineWidthType,1}
		ModifyGizmo ModifyObject=cellOutline objectType=path, property={ lineWidth,0.5}
		ModifyGizmo ModifyObject=cellOutline objectType=path, property={ pathColor, AtomView_CellOutLineColorR,AtomView_CellOutLineColorG,AtomView_CellOutLineColorB,AtomView_CellOutLineColorA}
#endif
	endif
	if (WaveExists(cell0))
		objectList += "cellOutline0;"
#if (IgorVersion()<7)
		Execute "AppendToGizmo Path="+GetWavesDataFolder(cell0,2)+",name=cellOutline0"
		Execute "ModifyGizmo ModifyObject=cellOutline0 property={ pathColorType,1}"
		Execute "ModifyGizmo ModifyObject=cellOutline0 property={ lineWidthType,1}"
		Execute "ModifyGizmo ModifyObject=cellOutline0 property={ lineWidth,"+num2str(AtomView_BondLineWidth)+"}"
		Execute "ModifyGizmo ModifyObject=cellOutline0 property={ pathColor,"+AtomView_CellOutLineColor+"}"
#else
		AppendToGizmo Path=$GetWavesDataFolder(cell0,2),name=cellOutline0
		ModifyGizmo ModifyObject=cellOutline0 objectType=path, property={ pathColorType,1}
		ModifyGizmo ModifyObject=cellOutline0 objectType=path, property={ lineWidthType,1}
		ModifyGizmo ModifyObject=cellOutline0 objectType=path, property={ lineWidth,AtomView_BondLineWidth}
		ModifyGizmo ModifyObject=cellOutline0 objectType=path, property={ pathColor,AtomView_CellOutLineColorR,AtomView_CellOutLineColorG,AtomView_CellOutLineColorB,AtomView_CellOutLineColorA}
#endif
	endif

	if (WaveExists(corners))
		objectList += "AtomViewCubeCorners;"
#if (IgorVersion()<7)
		Execute "AppendToGizmo Scatter="+GetWavesDataFolder(corners,2)+",name=AtomViewCubeCorners"
		Execute "ModifyGizmo ModifyObject=AtomViewCubeCorners property={ scatterColorType,0}"
		Execute "ModifyGizmo ModifyObject=AtomViewCubeCorners property={ markerType,0}"
		Execute "ModifyGizmo ModifyObject=AtomViewCubeCorners property={ sizeType,0}"
		Execute "ModifyGizmo ModifyObject=AtomViewCubeCorners property={ rotationType,0}"
		Execute "ModifyGizmo ModifyObject=AtomViewCubeCorners property={ Shape,1}"
		Execute "ModifyGizmo ModifyObject=AtomViewCubeCorners property={ size,1}"
		Execute "ModifyGizmo ModifyObject=AtomViewCubeCorners property={ color,0,0,0,1}"
		Execute "ModifyGizmo userString={CubeCorners,\"AtomViewCubeCorners\"}"	// save name of cube corner object
#else
		AppendToGizmo Scatter=$GetWavesDataFolder(corners,2),name=AtomViewCubeCorners
		ModifyGizmo ModifyObject=AtomViewCubeCorners,objectType=scatter,property={ Shape,1}
		ModifyGizmo ModifyObject=AtomViewCubeCorners,objectType=scatter,property={ size,0}
		SetWindow $gizName userdata(CubeCorners)="AtomViewCubeCorners"
#endif
	endif

#if (IgorVersion()<7)
	Execute "AppendToGizmo attribute blendFunc={770,771},name=blendingFunction"
	sprintf str,"{%g,%g,%g,1}",AtomView_GrayBkg,AtomView_GrayBkg,AtomView_GrayBkg
	Execute "ModifyGizmo setDisplayList=-1, opName=clearColor0, operation=clearColor, data="+str
#else
	AppendToGizmo attribute blendFunc={770,771},name=blendingFunction
	ModifyGizmo setDisplayList=-1, opName=clearColor0, operation=clearColor, data={AtomView_GrayBkg,AtomView_GrayBkg,AtomView_GrayBkg,1}
#endif

#if (IgorVersion()<7)
	if (strlen(titleGroup+scaleBarGroup))
		if (fixedLight)
			Execute "ModifyGizmo setDisplayList=-1, object=lightFixed"
		endif
		if (strlen(titleGroup))
			Execute "ModifyGizmo setDisplayList=-1, object="+titleGroup
		endif
		if (strlen(scaleBarGroup))
			Execute "ModifyGizmo setDisplayList=-1, object="+scaleBarGroup
		endif
		Execute "ModifyGizmo setDisplayList=-1, opName=MainTransform, operation=mainTransform"
	endif
	Execute "ModifyGizmo setDisplayList=-1, opName=orthoBase, operation=ortho, data={-2,2,-2,2,-3,3}"
	sprintf str, "ModifyGizmo setDisplayList=-1, opName=scaleBase, operation=scale, data={%g,%g,%g}",scaleFactor,scaleFactor,scaleFactor
	Execute str
	if (useBlend)
		Execute "ModifyGizmo setDisplayList=-1, attribute=blendingFunction"
		Execute "ModifyGizmo setDisplayList=-1, opName=enableBlend, operation=enable, data=3042"
	endif
#else
	if (strlen(scaleBarGroup))
		if (fixedLight)
			ModifyGizmo setDisplayList=-1, object=lightFixed
		endif
		if (strlen(scaleBarGroup))
			ModifyGizmo setDisplayList=-1, object=$scaleBarGroup
		endif
		ModifyGizmo setDisplayList=-1, opName=MainTransform, operation=mainTransform
	endif
	ModifyGizmo setDisplayList=-1, opName=orthoBase, operation=ortho, data={-2,2,-2,2,-3,3}
	ModifyGizmo setDisplayList=-1, opName=scaleBase, operation=scale, data={scaleFactor,scaleFactor,scaleFactor}
	if (useBlend)
		ModifyGizmo setDisplayList=-1, attribute=blendingFunction
		ModifyGizmo setDisplayList=-1, opName=enableBlend, operation=enable, data=3042
	endif
#endif

#if (IgorVersion()<7)
	for (i=0;i<ItemsInList(attributeList);i+=1)
		Execute "ModifyGizmo setDisplayList=-1, attribute="+StringFromList(i,attributeList)
	endfor
	for (i=0;i<ItemsInList(objectList);i+=1)
		Execute "ModifyGizmo setDisplayList=-1, object="+StringFromList(i,objectList)
	endfor
#else
	for (i=0;i<ItemsInList(attributeList);i+=1)
		ModifyGizmo setDisplayList=-1, attribute=$StringFromList(i,attributeList)
	endfor
	for (i=0;i<ItemsInList(objectList);i+=1)
		ModifyGizmo setDisplayList=-1, object=$StringFromList(i,objectList)
	endfor
#endif

#if (IgorVersion()<7)
	Execute "ModifyGizmo SETQUATERNION={0.398347,0.525683,0.596496,0.457349}"
	Execute "ModifyGizmo autoscaling=1"
	Execute "ModifyGizmo currentGroupObject=\"\""
	if (strlen(scaleBarGroup))
		Execute "ModifyGizmo namedHookStr={ScaleBarHook,\"GizmoUtil#GzimoReSetScaleBarHookProc\"}"
	endif
	if (strlen(title2))
		Execute "ModifyGizmo namedHookStr={TitleHook,\"AtomView#AtomViewGizmoFixHookProc\"}"
	endif
	Execute "ModifyGizmo compile"
	//	Execute "ModifyGizmo showInfo"
	Execute "ModifyGizmo bringToFront"
	Execute "ModifyGizmo endRecMacro"
#else
	ModifyGizmo SETQUATERNION={0.398347,0.525683,0.596496,0.457349}
	ModifyGizmo autoscaling=1
	ModifyGizmo currentGroupObject=""
	if (strlen(scaleBarGroup))
		ModifyGizmo namedHookStr={ScaleBarHook,"GizmoUtil#GzimoReSetScaleBarHookProc"}
	endif
	if (strlen(title2))
		ModifyGizmo namedHookStr={TitleHook,"AtomView#AtomViewGizmoFixHookProc"}
	endif
	ModifyGizmo compile
	//	ModifyGizmo showInfo
	ModifyGizmo bringToFront
	ModifyGizmo endRecMacro
#endif
	return gizName
End



Function/T AddRealLatticeAxesGroup(groupName,ain,bin,cin,[font,showNames])
	String groupName			// probably "RealLatticeAxesGroup0", If this is empty, then a unique name will be assigned
	Wave ain,bin,cin
	String font
	Variable showNames		// if true, show a,b,c labels on lattice vectors
	font = SelectString(ParamIsDefault(font),font,"Geneva")
	font = SelectString(strlen(font),"Geneva",font)
	showNames = ParamIsDefault(showNames) ? 0 : showNames
	showNames = numtype(showNames) ? 0 : !(!showNames)

	if (numpnts(ain)!=3 || numpnts(bin)!=3 || numpnts(cin)!=3)
		return ""
	endif
	Duplicate/FREE ain,a
	Duplicate/FREE bin,b
	Duplicate/FREE cin,c
	Variable len=max(norm(a),norm(b))
	len = max(len,norm(c))
	a /= len
	b /= len
	c /= len
	Make/N=3/D/FREE aM=-a, bM=-b, cM=-c

	Variable angle
	Make/N=4/D/FREE rotX, rotY, rotZ
	Make/N=3/D/FREE zhat={0,0,1}, hat
	Cross zhat,a
	Wave W_Cross=W_Cross
	W_Cross /= norm(a)
	angle = asin(normalize(W_Cross))*180/PI		// angle (degree)
	rotX[1,3] = W_Cross[p-1]
	rotX[0] = angle

	Cross zhat,b
	W_Cross /= norm(b)
	angle = asin(normalize(W_Cross))*180/PI		// angle (degree)
	rotY[1,3] = W_Cross[p-1]
	rotY[0] = angle

	Cross zhat,c
	W_Cross /= norm(c)
	angle = asin(normalize(W_Cross))*180/PI		// angle (degree)
	rotZ[1,3] = W_Cross[p-1]
	rotZ[0] = angle
	KillWaves/Z W_Cross

	if (CheckName(groupName,5))						// invalid groupName passed, create one
		if (strlen(groupName)<1)						// create a unique group name
			NewDataFolder/O root:Packages				// ensure Packages exists
			NewDataFolder/O root:Packages:JZT_GizmoUtility	// ensure geometry exists
			if (exists("root:Packages:JZT_GizmoUtility:RealLatticeNumber")!=2)
				Variable/G root:Packages:JZT_GizmoUtility:RealLatticeNumber=-1
			endif
			NVAR RealLatticeNumber = root:Packages:JZT_GizmoUtility:RealLatticeNumber
			RealLatticeNumber = numtype(RealLatticeNumber) ? -1 : limit(round(RealLatticeNumber),-1,Inf)
			RealLatticeNumber += 1
			groupName = "RealLatticeAxesGroup"+num2istr(RealLatticeNumber)
		endif
		groupName = CleanupName(groupName,0)
	endif
	if (CheckName(groupName,5))						// invalid groupName passed, give up
		return ""
	endif

	// ************************* Group Object Start *******************
#if (IgorVersion()<7)
	Execute "AppendToGizmo group,name="+groupName
	Execute "ModifyGizmo currentGroupObject=\""+groupName+"\""
	Execute "AppendToGizmo cylinder={0.02,0,0.1,25,25},name=cylinderArrow"
	Execute "AppendToGizmo line={0,0,0,"+vec2str(a,bare=1)+"}, name=lineX"
	Execute "AppendToGizmo line={0,0,0,"+vec2str(b,bare=1)+"}, name=lineY"
	Execute "AppendToGizmo line={0,0,0,"+vec2str(c,bare=1)+"}, name=lineZ"
	if (showNames)
		Execute "AppendToGizmo string=\"a\",strFont=\""+font+"\",name=string_a"
		Execute "ModifyGizmo modifyObject=string_a property={Clipped,0}"
		Execute "AppendToGizmo string=\"b\",strFont=\""+font+"\",name=string_b"
		Execute "ModifyGizmo modifyObject=string_b property={Clipped,0}"
		Execute "AppendToGizmo string=\"c\",strFont=\""+font+"\",name=string_c"
		Execute "ModifyGizmo modifyObject=string_c property={Clipped,0}"
	endif
	Execute "AppendToGizmo attribute lineWidth=2, name=lineWidthArrow"
	Execute "AppendToGizmo attribute color={1,0,0,1},name=colorRed"
	Execute "AppendToGizmo attribute color={0,1,0,1},name=colorGreen"
	Execute "AppendToGizmo attribute color={0,0,1,1},name=colorBlue"
	Execute "AppendToGizmo attribute specular={0,0,0,1,1032},name=specularOFF"
	Execute "AppendToGizmo attribute shininess={0,1032},name=shininessOFF"

	Execute "ModifyGizmo setDisplayList=0, attribute=shininessOFF"
	Execute "ModifyGizmo setDisplayList=1, attribute=lineWidthArrow"


	Execute "ModifyGizmo setDisplayList=-1, attribute=colorRed"
	Execute "ModifyGizmo setDisplayList=-1, object=lineX"
	Execute "ModifyGizmo setDisplayList=-1, opName=translateX, operation=translate, data={"+vec2str(a,bare=1)+"}"
	Execute "ModifyGizmo setDisplayList=-1, opName=rotateX, operation=rotate, data={"+vec2str(rotX,bare=1)+"}"
	Execute "ModifyGizmo setDisplayList=-1, object=cylinderArrow"
	if (showNames)
		Execute "ModifyGizmo setDisplayList=-1, opName=scale_a, operation=scale, data={0.1,0.1,0.1}"
		Execute "ModifyGizmo setDisplayList=-1, object=string_a"
		Execute "ModifyGizmo setDisplayList=-1, opName=scaleM, operation=scale, data={10,10,10}"
	endif
	rotX[0] = -rotX[0]
	Execute "ModifyGizmo setDisplayList=-1, opName=rotateMX, operation=rotate, data={"+vec2str(rotX,bare=1)+"}"
	Execute "ModifyGizmo setDisplayList=-1, opName=translateMX, operation=translate, data={"+vec2str(aM,bare=1)+"}"

	Execute "ModifyGizmo setDisplayList=-1, attribute=colorGreen"
	Execute "ModifyGizmo setDisplayList=-1, object=lineY"
	Execute "ModifyGizmo setDisplayList=-1, opName=translateY, operation=translate, data={"+vec2str(b,bare=1)+"}"
	Execute "ModifyGizmo setDisplayList=-1, opName=rotateY, operation=rotate, data={-90,1,0,0}"
	Execute "ModifyGizmo setDisplayList=-1, object=cylinderArrow"
	if (showNames)
		Execute "ModifyGizmo setDisplayList=-1, opName=scale_b, operation=scale, data={0.1,0.1,0.1}"
		Execute "ModifyGizmo setDisplayList=-1, object=string_b"
		Execute "ModifyGizmo setDisplayList=-1, opName=scaleMb, operation=scale, data={10,10,10}"
	endif
	Execute "ModifyGizmo setDisplayList=-1, opName=rotateMY, operation=rotate, data={90,1,0,0}"
	Execute "ModifyGizmo setDisplayList=-1, opName=translateMY, operation=translate, data={"+vec2str(bM,bare=1)+"}"

	Execute "ModifyGizmo setDisplayList=-1, attribute=colorBlue"
	Execute "ModifyGizmo setDisplayList=-1, object=lineZ"
	Execute "ModifyGizmo setDisplayList=-1, opName=translateZ, operation=translate, data={"+vec2str(c,bare=1)+"}"
	Execute "ModifyGizmo setDisplayList=-1, object=cylinderArrow"
	if (showNames)
		Execute "ModifyGizmo setDisplayList=-1, opName=scale_c, operation=scale, data={0.1,0.1,0.1}"
		Execute "ModifyGizmo setDisplayList=-1, object=string_c"
		Execute "ModifyGizmo setDisplayList=-1, opName=scaleMc, operation=scale, data={10,10,10}"
	endif
	Execute "ModifyGizmo setDisplayList=-1, opName=translateMZ, operation=translate, data={"+vec2str(b,bare=1)+"}"

	Execute "ModifyGizmo currentGroupObject=\"::\""
#else
	AppendToGizmo group,name=$groupName
	ModifyGizmo currentGroupObject=groupName

	AppendToGizmo line={0,0,0,a[0],a[1],a[2]}, name=lineX
	Modifygizmo modifyobject=lineX,objectType=line,property={colorType,2}
	Modifygizmo modifyobject=lineX,objectType=line,property={colorValue,0,1,0,0,1}
	Modifygizmo modifyobject=lineX,objectType=line,property={colorValue,1,1,0,0,1}
	Modifygizmo modifyobject=lineX,objectType=line,property={arrowMode,2}
	Modifygizmo modifyobject=lineX,objectType=line,property={endArrowHeight,0.1}
	Modifygizmo modifyobject=lineX,objectType=line,property={endArrowBase,0.03}

	AppendToGizmo line={0,0,0,b[0],b[1],b[2]}, name=lineY
	Modifygizmo modifyobject=lineY,objectType=line,property={colorType,2}
	Modifygizmo modifyobject=lineY,objectType=line,property={colorValue,0,0,1,0,1}
	Modifygizmo modifyobject=lineY,objectType=line,property={colorValue,1,0,1,0,1}
	Modifygizmo modifyobject=lineY,objectType=line,property={arrowMode,2}
	Modifygizmo modifyobject=lineY,objectType=line,property={endArrowHeight,0.1}
	Modifygizmo modifyobject=lineY,objectType=line,property={endArrowBase,0.03}

	AppendToGizmo line={0,0,0,c[0],c[1],c[2]}, name=lineZ
	Modifygizmo modifyobject=lineZ,objectType=line,property={colorType,2}
	Modifygizmo modifyobject=lineZ,objectType=line,property={colorValue,0,0,0,1,1}
	Modifygizmo modifyobject=lineZ,objectType=line,property={colorValue,1,0,0,1,1}
	Modifygizmo modifyobject=lineZ,objectType=line,property={arrowMode,2}
	Modifygizmo modifyobject=lineZ,objectType=line,property={endArrowHeight,0.1}
	Modifygizmo modifyobject=lineZ,objectType=line,property={endArrowBase,0.03}

	if (showNames)
		AppendToGizmo string="a",strFont=font,name=string_a
		AppendToGizmo string="b",strFont=font,name=string_b
		AppendToGizmo string="c",strFont=font,name=string_c
	endif
	AppendToGizmo attribute lineWidth=2, name=lineWidthArrow
	AppendToGizmo attribute specular={0,0,0,1,1032},name=specularOFF
	AppendToGizmo attribute shininess={0,1032},name=shininessOFF

	ModifyGizmo setDisplayList=0, attribute=shininessOFF
	ModifyGizmo setDisplayList=1, attribute=lineWidthArrow

	ModifyGizmo setDisplayList=-1, object=lineX
	ModifyGizmo setDisplayList=-1, opName=translateX, operation=translate, data={a[0],a[1],a[2]}
	ModifyGizmo setDisplayList=-1, opName=rotateX, operation=rotate, data={rotX[0],rotX[1],rotX[2],rotX[3]}
	if (showNames)
		ModifyGizmo setDisplayList=-1, opName=scale_a, operation=scale, data={0.1,0.1,0.1}
		ModifyGizmo setDisplayList=-1, object=string_a
		ModifyGizmo setDisplayList=-1, opName=scaleM, operation=scale, data={10,10,10}
	endif
	rotX[0] = -rotX[0]
	ModifyGizmo setDisplayList=-1, opName=rotateMX, operation=rotate, data={rotX[0],rotX[1],rotX[2],rotX[3]}
	ModifyGizmo setDisplayList=-1, opName=translateMX, operation=translate, data={-a[0],-a[1],-a[2]}
	ModifyGizmo setDisplayList=-1, object=lineY
	ModifyGizmo setDisplayList=-1, opName=translateY, operation=translate, data={b[0],b[1],b[2]}
	ModifyGizmo setDisplayList=-1, opName=rotateY, operation=rotate, data={-90,1,0,0}
	if (showNames)
		ModifyGizmo setDisplayList=-1, opName=scale_b, operation=scale, data={0.1,0.1,0.1}
		ModifyGizmo setDisplayList=-1, object=string_b
		ModifyGizmo setDisplayList=-1, opName=scaleMb, operation=scale, data={10,10,10}
	endif
	ModifyGizmo setDisplayList=-1, opName=rotateMY, operation=rotate, data={90,1,0,0}
	ModifyGizmo setDisplayList=-1, opName=translateMY, operation=translate, data={-b[0],-b[1],-b[2]}

	ModifyGizmo setDisplayList=-1, object=lineZ
	ModifyGizmo setDisplayList=-1, opName=translateZ, operation=translate, data={c[0],c[1],c[2]}
	if (showNames)
		ModifyGizmo setDisplayList=-1, opName=scale_c, operation=scale, data={0.1,0.1,0.1}
		ModifyGizmo setDisplayList=-1, object=string_c
		ModifyGizmo setDisplayList=-1, opName=scaleMc, operation=scale, data={10,10,10}
	endif
	ModifyGizmo setDisplayList=-1, opName=translateMZ, operation=translate, data={b[0],b[1],b[2]}

	ModifyGizmo currentGroupObject="::"
#endif
	// ************************* Group Object End *******************

	return groupName
End



Static Function AtomViewGizmoFixHookProc(s)
	STRUCT WMGizmoHookStruct &s
	if (!StringMatch(s.eventName,"scale"))
		return 0
	endif
	String win=s.winName

	Wave xyz = $GetItemFromGizmoObject(win,"atomViewAtoms","scatter")
	if (!WaveExists(xyz))
		return 0
	endif

	String wnote=note(xyz), title2="", str
	str = StringByKey("Nabc",wnote,"=")
	if (strlen(str)>1)
		title2 = ReplaceString(" ",str," x ")+" cells"	// new title2 string
	endif
	if (strlen(title2)<1)
		return 0
	endif

#if (IgorVersion()<7)
	Execute "GetGizmo/N="+win+"/Z objectNameList"
	String ObjectNames = StrVarOrDefault("S_ObjectNames","")
	KillStrings/Z S_ObjectNames

	String titleGroupName=""									// find the name of the title group
	Variable i
	for (i=0;i<ItemsInList(ObjectNames);i+=1)
		if (stringmatch(StringFromList(i,ObjectNames),"gizmoStringGroup*"))
		titleGroupName = StringFromList(i,ObjectNames)
		endif
	endfor
	if (strlen(titleGroupName)<1)
		return 0
	endif

	Execute "ModifyGizmo/N="+win+"/Z startRecMacro"
	Execute "ModifyGizmo/N="+win+"/Z currentGroupObject=\""+titleGroupName+"\""
	Execute "ModifyGizmo/N="+win+"/Z modifyObject=Title2, property={string,\""+title2+"\"}"
	Execute "ModifyGizmo/N="+win+"/Z currentGroupObject=\"::\""
	Execute "ModifyGizmo/N="+win+"/Z endRecMacro"
#else
	String desc = StringByKey("desc",wnote,"="), title3=""
	Variable bondLenMax=NumberByKey("bondLenMax",wnote,"=")
	if (bondLenMax>0)
		sprintf title3,"max bond length = %g nm",bondLenMax
	endif
	String sourceFldr = StringByKey("sourceFldr",wnote,"=")
	String titleGroup = AddGizmoTitleGroup("textTitle",desc,title2=title2,title3=title3,title4=sourceFldr)
#endif
	return 0	 
End

//  ============================================================================  //
//  ============================== End Make Gizmo ==============================  //



//  ============================================================================  //
//  =========================== Start Initialization  =========================== //

Function Init_AtomViewLattice()
	InitLatticeSymPackage()
	ElementDataInitPackage()
	GMarkers#InitGizmoMarkers()
	GZoomTrans#InitGizmoZoomTranslate()
	GClipPlanes#InitGizmoClipPlanes()
	return 0
End

//  ============================ End Initialization  ============================ //
//  ============================================================================  //

