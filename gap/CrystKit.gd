
#############################################################################
##
#F  CaratNormalizedInputFile( filename ) . . . . . returns a new normalized 'resfile' in Carat tmp directory
##
DeclareGlobalFunction( "CaratNormalizedInputFile" );

#############################################################################
##
#F  CaratName( S ) . . . . . Call Carat program `Name'
##
DeclareGlobalFunction( "CaratName" );

#############################################################################
##
#F  CaratReverse_name( name ) . . . . . Call Carat program `Reverse_name'
##
DeclareGlobalFunction( "CaratReverse_name" );
DeclareGlobalFunction( "CaratQ_catalog" );

DeclareGlobalFunction( "IdentifyGroupGenerators3d" );
DeclareGlobalFunction( "IdentifyGroupGenerators4d" );
DeclareGlobalFunction( "AugmentedMatrixOnLeft" );

DeclareGlobalFunction( "PhipnInverse" );
DeclareGlobalFunction( "NewFiniteOrdersOfGLNZ" );
DeclareGlobalFunction( "OrbitOfStandardAffineCrystGroupOnLeftByNormalizerPointGroup" );

DeclareGlobalFunction( "AffineIsomorphismSpaceGroups" );
DeclareGlobalFunction( "IdentifySpaceGroup" );
DeclareGlobalFunction( "MinimalGeneratingSetSolvableSpaceGroupByPcpGroup" );


DeclareGlobalVariable( "cryst2names" );
DeclareGlobalVariable( "cryst3names" );
DeclareGlobalVariable( "aff3names" );
DeclareGlobalVariable( "aff4names" );
