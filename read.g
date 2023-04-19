#
# CrystKit: The crystallographic groups kit based on GAP related packages and interfaces to other 3rd tools.
#

# read some help datasets to work with crystallographic groups
ReadPackage( "CrystKit", "gap/names/cryst2names.gi" );
ReadPackage( "CrystKit", "gap/names/cryst3names.gi" );
ReadPackage( "CrystKit", "gap/names/aff3names.gi" );
ReadPackage( "CrystKit", "gap/names/aff4names.gi" );

#[GAP Forum] Needed vs. suggested packages
#https://mail.google.com/mail/u/0/?ogbl#inbox/FMfcgzGsmDrQkbKDFWndHDTZkxXRwMxB
# Read Browse applications only if the Browse package will be loaded.
#if IsPackageMarkedForLoading( "Browse", ">= 1.8.3" ) then
#  ReadPackage( "atlasrep", "gap/brmindeg.g" );
#  if IsPackageMarkedForLoading( "CTblLib", "" ) then
#    ReadPackage( "atlasrep", "gap/brspor.g" );
#  fi;
#fi;

# Reading the implementation part of the package.
if IsPackageMarkedForLoading( "RepnDecomp", ">= 1.3.0" ) then
  ReadPackage( "CrystKit", "gap/CrystKit.gi");
fi;

