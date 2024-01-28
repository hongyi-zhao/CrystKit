##The GAP package CrystKit

#Description

CrystKit is a GAP package designed for the manipulation and analysis of crystallographic groups. It offers a variety of functionalities, including space group simplification on any basis, minimal generating set computation for solvable space groups, isomorphism detection between space groups, enantiomorphic pair existence checking, efficient orbit computation for the space group vector system (SNoT), and rapid identification of space groups (up to six-dimension) based on CARAT. It also encompasses the calculation of allowed orders of point groups in any dimensional space groups (an extension of the crystallographic restriction theorem) and conjugacy determination between two finite matrix groups.

#Main Features

Simplification of Space Groups Represented on Any Basis: Transparent simplification to more manageable generator forms.

```gap
SGGenSet227me:=[ [ [ 0, -1, 0, 1/2 ], [ 0, 0, -1, 1/2 ], [ -1, 0, 0, 1/2 ], [ 0, 0, 0, 1 ] ], [ [ -15/4, 29/4, -15/4, -15/16 ], [ -33/8, 55/8, -25/8, -25/32 ], 
      [ -25/8, 55/8, -33/8, -41/32 ], [ 0, 0, 0, 1 ] ] ];
S:=AffineCrystGroupOnLeft(SGGenSet227me);
conj1:=ConjugatorSpaceGroupSimplification(S);
S1:=S^(conj1^-1);
GeneratorsOfGroup(S1);
```

Minimal Generating Set Calculation for Solvable Space Groups: Identifying the minimal number of generators for 3D space group 227.

```gap
S2:=SpaceGroupOnLeftIT(3,227);
Length(GeneratorsOfGroup(S2));
minSgen:=MinimalGeneratingSetAffineCrystGroup(S2);
Length(minSgen.mgs);
```

Isomorphism Determination Between Space Groups Represented on Any Basis: Checking isomorphism and directly comparing two space groups.

```gap
AffineIsomorphismSpaceGroups(S,S2);
S^(last^-1)=S2;
```

Existence Check for Enantiomorphic Pairs in Any Basis Represented Space Groups: Finding enantiomorphic partners for given space groups.

```gap
S4:=SpaceGroupOnLeftIT(3,212);
conj2:=ConjugatorSpaceGroupEnantiomorphicPartner(S4);
S5:=S4^(conj2^-1);

```

Efficient Calculation of Orbits of the Space Group Vector System (SNoT): Orbit calculation using the normalizer of the point group.

```gap
P4:=PointGroup(S4);
d:= DimensionOfMatrixGroup(S4) - 1;
norm := GeneratorsOfGroup(Normalizer(GL(d,Integers), P4)); 
orbnpg := OrbitSpaceGroupStdByNormalizerPointGroup(S4, GeneratorsOfGroup(P4), norm);
```

Rapid Identification of Space Groups (Up to Six Dimensions) Based on Any Basis (using CARAT): Quick identification methods for space groups.

```gap
IdentifySpaceGroup(S4);
IdentifySpaceGroup(S5);
```

Allowed Orders Calculation for Point Groups in Any Dimensional Space Groups: Extends the crystallographic restriction theorem by calculating new finite orders of GLNZ.

```gap
NewFiniteOrdersOfGLNZ(2);
NewFiniteOrdersOfGLNZ(4);
```

Conjugacy Determination Between Two Finite Matrix Groups (using RepnDecomp): Determines if two point groups are conjugate.

```gap
P:=PointGroup(S);;
S6:=SpaceGroupOnLeftIT(3,222);;
P6:=PointGroup(S6);;
ConjugatorMatrixGroups(P, P6); 
```

#Install

CrystKit is distributed with GAP, and does not require any installation. It is loaded with the GAP command

```
gap> LoadPackage( "CrystKit" ); 
```

#Contact

For bug reports, suggestions and other comments please use the issue tracker on the GitHub page of the package:

https://github.com/hongyi-zhao/CrystKit/issues

For further information, contact us at: <hongshengzhao@xpc.edu.cn> or <hongyi.zhao@gmail.com>.

