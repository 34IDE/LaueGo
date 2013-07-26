#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma version = 2.0

Static Function IgorStartOrNewHook(IgorApplicationNameStr)
	String IgorApplicationNameStr

	Variable g = NumVarOrDefault("root:Packages:MICRO_GEOMETRY_VERSION",NaN)
	if (numtype(g) && !DataFolderExists("root:Packages:"))
		Execute/P "JonFirst#setMICRO_GEOMETRY_VERSION_PATH()"
//		Execute/P "JonFirst#setMICRO_GEOMETRY_VERSION_PATH(preset=1)"
		Execute/P "INSERTINCLUDE  \"mdaFiles\", version>=1.08";Execute/P "COMPILEPROCEDURES "
	endif
End
