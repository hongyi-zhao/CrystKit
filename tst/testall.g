#
# CrystKit: The crystallographic groups kit based on GAP related packages and interfaces to other 3rd tools.
#
# This file runs package tests. It is also referenced in the package
# metadata in PackageInfo.g.
#
LoadPackage( "CrystKit" );

TestDirectory(DirectoriesPackageLibrary( "CrystKit", "tst" ),
  rec(exitGAP := true));

FORCE_QUIT_GAP(1); # if we ever get here, there was an error
