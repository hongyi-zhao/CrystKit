# https://www.cryst.ehu.es/cgi-bin/cryst/programs/checkgr.pl?tipog=gesp
# https://iso.byu.edu/iso/findssghelp.php
# The output notation for superspace-group operators will match that of the input: There are three choices: (x,y,z,t,u,v), (x1,x2,x3,x4,x5,x6), and (xs1,xs2,xs3,xs4,xs5,xs6). See the ISO(3+d)D help page for more information about these notations. 

# Define variable-length `Indeterminate's in function.
# https://mail.google.com/mail/u/0/?ogbl#search/in%3Asent+gap/QgrcJHsBqxpXbHCSnPqBGGBZRDqphmxvpZb
InstallGlobalFunction( IdentifyGroupGenerators, function( S )
  
  local gens, vecname, d, vec, g, i;
  
  # BCS IDENTIFY GROUP 和 findssg 都是工作于左作用下面的标准表示的。
  if IsAffineCrystGroupOnRight(S) then
    S:=TransposedMatrixGroup(S);
  fi;

  if not IsStandardSpaceGroup(S) then
    S:=StandardAffineCrystGroup(S);
  fi;

  gens:= GeneratorsOfGroup(S);
  vecname:=["x","y","z","t","u","v"];
  d := First(Set(DimensionsMat(gens[1]))) - 1;
  
  if d > Size(vecname) then
    Error("Not yet supported");
  fi;

  vec := List([1..d],x -> X(Rationals,x)); 
  gens:=List( gens , g -> g{[1..d]}{[1..d]} * vec +g{[1..d]}[d+1] );
 
  for i in [1..d] do
    SetName(vec[i],vecname[i]);
  od;
  
  return gens;

end );


InstallGlobalFunction( AugmentedMatrixOnLeft, function(m, b)
  
  local d, g, i;

  if Size( Set(DimensionsMat(m)) ) <> 1 then
    Error( "linear part is not a square matrix" );    
  fi;
  d := First(Unique(DimensionsMat(m)));
  g := List([1..d], i -> Concatenation(m[i], [b[i]]));
  Add(g, Concatenation(Zero([1..d]), [1]));
  return g;

end );


#############################################################################
##
#F  CaratNormalizedInputFile( filename ) . . . . . returns a new normalized
#  'resfile' in Carat tmp directory
##
InstallGlobalFunction( CaratNormalizedInputFile, function( filename )

  local resfile, input, str, x, y;

  resfile := CaratTmpFile( "res" );  
  input := InputTextFile( filename );
  # RewindStream(input);
  str := ReadAll( input );
  CloseStream( input );

  str := Filtered(SplitString(str, "\n"), x -> not IsEmpty(x));
  str := Filtered(str, x -> not ForAll(SplitString(x, " "), y -> IsEmpty(y)));    
  str:=Chomp(JoinStringsWithSeparator(str,"\n")); 
  # CaratShowFile( resfile );
  PrintTo(resfile, str);
  
  return resfile;

end );

# Do the following given matrices generate a space group? If so, find its name and normal representation. 
# https://lbfm-rwth.github.io/carat/doc/examples/Ex5.html
#############################################################################
##
#F  CaratName( S ) . . . . . Call Carat program `Name'
##
InstallGlobalFunction( CaratName, function( S )
  
  local Str, Sgen, out, c, C,
        dir, shell, d, cmd, program, cs, args, 
        _in, _out, res, err, x; 

  # https://docs.gap-system.org/pkg/caratinterface/htm/CHAP003.htm#SECT001
  # In crystallography, the convention usually is that matrix groups act from the left on column vectors. This convention is adopted also in CARAT.
  if IsAffineCrystGroupOnRight( S ) then
    Str := TransposedMatrixGroup( S );
    return CaratName( Str );
  fi;

  # get temporary file names
  Sgen := CaratTmpFile( "Sgen" );  
  out := CaratTmpFile( "out"  );

  d := DimensionOfMatrixGroup( S ) - 1;

  # `Name' is based on the CARAT database of all Q-classes of finite unimodular groups of degree up to 6. 
  if not IsSpaceGroup(S) or d > 6 then
    Error( "only works for space groups of degree up to 6" );
  fi;


  # compatiable with the result returned by `Name' in this case:
  if IsTrivial( PointGroup( S ) ) then
    return fail;
  fi;

  # `Name' is not able to perform calculations directly using proper cyclotomics. Therefore, a conversion process is carried out beforehand:

  # 用下面的解决方法来首先彻底简化已给空间群的表示：
  C:=ConjugatorReducedSpaceGroup(S); 
  S:=S^(C^-1);

  CaratWriteMatrixFile(Sgen, GeneratorsOfGroup( S ));
  # CaratShowFile(Sgen);

  dir := DirectoryCurrent();
  shell := Filename(DirectoriesSystemPrograms(), "sh");
  _in:=InputTextNone();
  _out:=OutputTextFile(out, false);

  cmd := "Name";
  program := Filename(CARAT_BIN_DIR, cmd);
  cs := "-c";
  args :=Concatenation(program, " ", Sgen, " ", "-M -c 2>/dev/null");
  err := Process(dir, shell, _in, _out, [cs, args]);
  
  CloseStream( _in );
  CloseStream( _out );

  # S is not a space group: 
  if err <> 0 then
    return fail;
  fi;

  # read back the result
  _out := InputTextFile( out );
  # RewindStream( _out );
  res:=Chomp( ReadAll( _out ) );
  CloseStream( _out );

  res:=SplitString(res, "\n");
  if First(res) = Last(res) then
    res:=[First(res)];
  else
    # For d = 2 or 3
    res:=[First(res),Last(res)]; 
  fi;
    
  res[1] := SplitString(res[1], "-");
  res[1,2] := JoinStringsWithSeparator(SplitString(res[1,2],"."), " ");

  # remove temporary files
  # borrowed from GAP package GUAVA
  #    This only works on unix
  #    Exec(Concatenation("rm -f ",f));
  #    This is better:
  List([Sgen, out], RemoveFile);
 
  return res;

end );


# Ref. BravaisGroupsCrystalFamily
#############################################################################
##
#F  CaratReverse_name( name ) . . . . . Call Carat program `Reverse_name'
##
InstallGlobalFunction( CaratReverse_name, function( name )

    local out, _in, _out,
          dir, shell, cmd, program, cs, args,
          err, gens, x; 
  
    # get temporary file name 
    out := CaratTmpFile( "out" ); 
    
    # execute in the current directory
    dir := DirectoryCurrent();
 
    # select the shell, bourne shell is the default: sh -c cmd
    shell := Filename( DirectoriesSystemPrograms(), "sh" );

    _in := InputTextString( Concatenation(List( name, x -> Concatenation(x, "\n") ))  ); 

    #  OutputTextFile(filename, append)
    _out := OutputTextFile( out, false ); 

    # find executable 
    cmd := "Reverse_name"; 
    # CaratHelp(command);
    program := Filename( CARAT_BIN_DIR, cmd ); 
 
    cs := "-c";
    args := Concatenation(program, " ", "2>/dev/null");
    # Run the following command to facilitate repeated debugging
    RewindStream( _in );
    # execute the command
    err := Process( dir, shell, _in, _out, [ cs, args ] );
    # CaratShowFile( out );
    CloseStream( _in );
    CloseStream( _out );

    # did it work? 
    if err <> 0  then 
      Error( Concatenation( "Carat program ", cmd, 
                            " failed with error code ", String(err) ) ); 
    fi; 

    # read Carat result from file, and remove temporary file 
    gens := CaratReadMatrixFile( out ); 
    RemoveFile( out ); 
    
    return gens; 

end );

# Find all Z-classes, affine classes and torsion free space groups in the Q-class of min.108. 
# https://lbfm-rwth.github.io/carat/doc/examples/Ex12.html
InstallGlobalFunction( CaratQ_catalog, function( str, out )

  local dir, shell, instr, _in, _out,
        cmd, program, cs, args,
        err, x; 
       
    # execute in the current directory
    dir := DirectoryCurrent();
 
    # select the shell, bourne shell is the default: sh -c cmd
    shell := Filename( DirectoriesSystemPrograms(), "sh" );
    
    # _in := InputTextString( str );
    # 该命令的 str 末尾要求必须以空格隔开，后接 out：
    str := Filtered(SplitString(str, "\n"), x -> not IsEmpty(x));
    str := JoinStringsWithSeparator(str, "\n");
    instr := Concatenation( str, " ", out, "\nq"  );
    _in := InputTextString( instr );

    # get temporary file name 
    # out := CaratTmpFile( "out" ); 
    #  OutputTextFile(filename, append)
    _out := OutputTextNone();
    # _out := OutputTextFile( out, false ); 

    # find executable 
    cmd := "Q_catalog"; 
    program := Filename( CARAT_BIN_DIR, cmd ); 
 
    cs := "-c";
    args := Concatenation(program, " ", ">/dev/null");
    # Run the following command to facilitate repeated debugging
    RewindStream( _in );
    # execute the command
    err := Process( dir, shell, _in, _out, [ cs, args ] );
    # CaratShowFile( out );
    CloseStream( _in );
    CloseStream( _out );
    
    if err = 0 then
      return out;
    else
      return fail;
    fi;

end );


InstallGlobalFunction( 
IdentifySpaceGroup, function( S )
  
  local d, res, name, nr, nrs, Sref, c, C, i;

  if IsAffineCrystGroupOnRight( S ) then
    S := TransposedMatrixGroup( S );
    # 此时的进一步递归处理：
    res := IdentifySpaceGroup( S );
    if res <> fail then
      res[2] := TransposedMat(res[2]) ^-1;
    fi;
    return res;
  fi;

  d := DimensionOfMatrixGroup( S ) - 1;
  # CrystCatRecord(TransposedMatrixGroup(S)).parameters;
  
  if d > 6 or not IsSpaceGroup( S ) then
    Error("only applicable to space groups in dimensions up to 6");
  fi;

  name:=CaratName( S );
  
  # S is not a space group:
  if name = fail then
    return fail;
  fi;

  if d = 2 then 
    nr:=Position(cryst2names, name);
    Sref:=SpaceGroupOnLeftIT(d,nr);
    C := AffineIsomorphismSpaceGroups(S, Sref);
  elif d = 3 then
      
    # cryst3names 为最精细分类，故aff3names没有必要使用：
    # Positions(aff3names, name );
    nrs:=Positions(cryst3names, name );
    
    # determine from the enantiomorphic_pairs
    if Size(nrs) = 1 or 0 > DeterminantMat(AffineIsomorphismSpaceGroups(S,SpaceGroupOnLeftIT(d,nrs[2]))) then
      nr := nrs[1];
    else
      nr := nrs[2];
    fi;
    
    Sref:=SpaceGroupOnLeftIT(d,nr);
    C := AffineIsomorphismSpaceGroups(S, Sref);

  elif d = 4 then
    
    nr:=Position(aff4names, name);
    Sref:=TransposedMatrixGroup(SpaceGroup(d,nr));
    C := AffineIsomorphismSpaceGroups(S, Sref);

  else
    nr := fail;
    c := TransposedMat(InternalBasis(S));
    C := AugmentedMatrixOnLeft( c, 0 *[1..d] );
  fi;
 
  if nr <> fail then
    res := [Concatenation([nr], name), C];
    return res;
  else
    res := [name[1], C];
    return res;
  fi;

end );


# Create the parity transformation matrices. 
# https://community.wolfram.com/groups/-/m/t/2829458?p_p_auth=mg22BhW6
# Create the odd parity transformation matrices.
# https://mail.google.com/mail/u/0/?ogbl#sent/QgrcJHsbjCZSdRXDtbTvBvxRmDjWSjVzqXL


# https://www.britannica.com/science/parity-particle-physics
# parity, in physics, property important in the quantum-mechanical description of a physical system. In most cases it relates to the symmetry of the wave function representing a system of fundamental particles. A parity transformation replaces such a system with a type of mirror image. Stated mathematically, the spatial coordinates describing the system are inverted through the point at the origin; that is, the coordinates x, y, and z are replaced with −x, −y, and −z. In general, if a system is identical to the original system after a parity transformation, the system is said to have even parity. If the final formulation is the negative of the original, its parity is odd. For either parity the physical observables, which depend on the square of the wave function, are unchanged. A complex system has an overall parity that is the product of the parities of its components.


# https://en.wikipedia.org/wiki/Parity_(physics)
# In physics, a parity transformation (also called parity inversion) is the flip in the sign of one spatial coordinate. In three dimensions, it can also refer to the simultaneous flip in the sign of all three spatial coordinates (a point reflection)
# It can also be thought of as a test for chirality of a physical phenomenon, in that a parity inversion transforms a phenomenon into its mirror image.

# Construct all the possiable parity transformation corresponding to the space group dimension under consideration:
# InstallGlobalFunction( AffineTransformationListOfParityOperator, function( d ) 
  
#   local par, ppar, apar, i, x;
#   # d:=3;
#   # parity transformation (also called parity inversion)
#   # par:=List(Filtered([1..d], IsOddInt), x -> Concatenation(-List([1..x], i -> i^0), List([1..d-x], i -> i^0)));
#   # or
#   par:=List(Filtered([1..d], IsOddInt), x -> Concatenation(-List([1..x], i -> 1), List([1..d-x], i -> 1)));
#   ppar:=Union(List(par, x -> PermutationsList(x)));
#   # conjugator of parity 
#   apar:=List(ppar, x -> DiagonalMat(Concatenation(x,[1])));

#   return apar;
# end );

# Find the highest finite order group element of an infinite group.
# https://mail.google.com/mail/u/0/?ogbl#sent/KtbxLwHLtgJlPSGRwNfgCnBVmRxJDfvqdV
#  The highest finite order of an element of your group is 6, cf. https://en.wikipedia.org/wiki/Crystallographic_restriction_theorem .

# Hope this helps,

#     Stefan


#  Thank you very much for your explanations!

# Well -- I don't know why a couple of numbers are missing in the entry for n = 24
# in the paper you refer to. -- Most likely it was just a clerical error
# (I don't see how an algorithm may plausibly miss exactly the last 9 numbers,
# but forgetting to paste a line of text may happen relatively easily).  

# By the way -- I had a brief look what your function OrderByLatticeDimensionCrt needs to do,
# and wrote an ad-hoc function NewFiniteOrdersOfGLNZ with the same functionality:

# Thanks to Stefan Kohl <sk239@st-andrews.ac.uk> for providing the following function:
InstallGlobalFunction( PhipnInverse,  function ( n )

  local  qs, p, k;

  qs := [];
  k  := 0;
  repeat
    k := k + 1;
    p := RootInt(n,k) + 1;
    if IsPrime(p) then
      if n = (p-1) * p^(k-1) then
        Add(qs,p^k);
      elif p = 2 then
        break;
      fi;
    fi;
  until false;
  if SmallestRootInt(n) = 2 then
    Add(qs,2*n);
  fi;
  return Set(qs);
end );

# Thanks to Stefan Kohl <sk239@st-andrews.ac.uk> for providing the following function:
InstallGlobalFunction( NewFiniteOrdersOfGLNZ, function ( n )

  local  orders, build, P, P2, choices, invs;

  build := function ( choices, pos, ord )

    local  q;

    for q in choices[pos] do
      if Gcd(ord,q) = 1 then
        if pos = Length(choices) then
          Add(orders,q*ord);
        else
          build(choices,pos+1,q*ord);
        fi;
      fi;
    od;
  end;

  if not IsPosInt(n) or n mod 2 = 1 then return []; fi;
  invs := List([1..n],PhipnInverse);
  orders := [];
  for P in Partitions(n/2) do
    P2 := 2 * P;
    choices := invs{P2};
    build(choices,1,1);
  od;
  orders := Union(orders,2*Filtered(orders,n->n mod 2 = 1));
  return orders;

end );

# For n = 100, I observe similar timings as those you report --
# and the above code is not yet optimized in any way.

# NewFiniteOrdersOfGLNZ(200); time;

# Best wishes,

#     Stefan


# 3.4 Methods provided by CARAT
# https://docs.gap-system.org/pkg/caratinterface/htm/CHAP003.htm
# CARAT implements methods for the following functions and operations declared in the GAP library. For a detailed description of these functions, please consult the GAP manual (section Matrix Groups in Characteristic 0). 
# https://docs.gap-system.org/doc/ref/chap44.html#X7FB0138F79E8C5E7


# https://en.wikipedia.org/wiki/Space_group#Bieberbach's_theorems
# https://en.wikipedia.org/wiki/Hilbert%27s_eighteenth_problem
# A GEOMETRIC PROOF OF BIEBERBACH'S THEOREMS
# ON CRYSTALLOGRAPHIC GROUPS
# by Peter Buser
# Pour Ariane et Georges
# https://www.e-periodica.ch/cntmng?pid=ens-001:1985:31::54
# In 1910 Bieberbach proved two celebrated theorems in response to
# Hilbert's 18th problem.

# Groups which satisfy the hypothesis of Theorem I are called n-dimensional
# crystallographic groups.

# Theorem I. Every discrete group of isometries acting on the n-dimensional
# euclidean space R" with compact fundamental domain contains n linearly
# independent translations.

# Theorem II. For each fixed n there are only finitely many isomorphism
# classes of n-dimensional crystallographic groups.

InstallGlobalFunction( OrbitSpaceGroupStdByNormalizerPointGroup, function( S, Pgen,norm )

  local P, d, Sgen, M, hom,
        orb, rep, t, tau, newtau, newrep,
        Snewgen, len, neworb, n, o, pos,
        res, g, x, y,
        # For debug:
        reppos, tgen, Simg, Snew, lst, i, ll;

  if not IsStandardAffineCrystGroup( S ) then
    Error("only work with StandardAffineCrystGroup");
  fi;

  if IsAffineCrystGroupOnRight( S ) then
    S := TransposedMatrixGroup( S );
    res := OrbitSpaceGroupStdByNormalizerPointGroup( S );
    res.M := TransposedMat( res.M );
    res.rep := List(res.rep, TransposedMat);
    return res;
  fi;

  hom := PointHomomorphism(S);
  d := DimensionOfMatrixGroup(S) - 1;
  P := PointGroup(S);

  if IsTrivial(P) then
    # P:=Group(IdentityMat(d));
    P:=Group( One(GL(d, Integers)) );
  fi;


  # Pgen := GeneratorsOfGroup(P);
  Sgen := List(Pgen, x -> PreImagesRepresentative(hom, x));
  Sgen := List(Sgen, x -> AugmentedMatrixOnLeft(x{[1..d]}{[1..d]}, List(x{[1..d]}[d+1], FractionModOne)));
  

  # norm := GeneratorsOfGroup(Normalizer(GL(d, Integers), P));
  norm := List(Filtered(norm, x -> not x in P), y -> AugmentedMatrixOnLeft(y, 0*[1..d]));


  # Definition of symmorphic space groups:
  # Algorithms for Crystallographic Groups
  # BETTINA EICK,1 BERND SOUVIGNIER2
  # page 318
  # 7. Definition.
  # 2. Space groups containing a subgroup isomorphic to their full point group are called symmorphic space groups. This is the case if and
  # only if the image of \tau lies in Z^n.
  t := List(Concatenation(List( Sgen, x -> x{[1..d]}[d+1] )), FractionModOne); 
  
  orb := [ Sgen ];
  tau := [ t ];
  rep := [ One(S) ];
  M := Concatenation( List( Pgen, g -> g - IdentityMat(d) ) );

  # For debug:
  reppos := [ [One(S), 1] ];
  tgen := List(IdentityMat(d), x -> AugmentedMatrixOnLeft(IdentityMat(d), x));

  # catch the trivial cases and 
  # SymmorphicSpaceGroup，直接返回结果即可。
  if not ( 
          IsTrivial(P) or 
          IsEmpty(norm) or 
          # SymmorphicSpaceGroup
          ForAll( tau, IsZero ) 
          ) then

    # 因为 norm 是作用在点群上的，所以必须首先枚举完
    # norm 作用下的 Pgen 的所有可能变化，
    # 然后，再进一步处理。
    repeat
      len := Size(orb);
      for n in norm do
        for o in orb do
          pos := Position(orb, o);
          neworb := OnTuples(o, n);
          neworb := List(neworb, x -> AugmentedMatrixOnLeft(x{[1..d]}{[1..d]}, List(x{[1..d]}[d+1], FractionModOne)));

          if not neworb in orb then
            newrep := rep[pos] * n;
               
            # Simg:=AffineCrystGroupOnLeft(Concatenation(neworb, tgen));
            # 对点群部分生成元用 newrep{[1..d]}{[1..d]}^-1 作用，
            # 此时，必然点群不变，然后再经过 hom 在 S 中找出对应的 Sgen，
            # 这样实际上，必然仍生成 S，
            # 由此，进一步得到维持和原Sgen的点群部分对应的 Snewgen 的方法，
            # 同时使得 Simg = Snew。
            # 这样，就实现了保持点群部分全同的轨道表示。
            # Snewgen:=List(Pgen, x -> PreImagesRepresentative(hom, x ^ (newrep{[1..d]}{[1..d]}^-1) ) ^ newrep );
            # or
            Snewgen:=List( OnTuples(Pgen, newrep{[1..d]}{[1..d]}^-1), x -> PreImagesRepresentative( hom, x ) ^ newrep );
            
            # Snew:=AffineCrystGroupOnLeft( Concatenation( Snewgen, tgen) );
            # Print( S^(newrep^-1) = Simg, " ", Simg = Snew, "\n" );
                
            newtau:=List(Concatenation(List( Snewgen, x -> x{[1..d]}[d+1] )), FractionModOne);
            
            # tau 中的两项所对应的 SG 之间可以通过纯平移共轭同构，其中包括了它们相等的情况（零解）。
            # 因此单独使用 ForAll 也是可以的，但是基于第一个条件可以提高效率，避免不必要的计算：
            if not newtau in tau and ForAll(List(tau, x -> SolveInhomEquationsModZ( M, newtau - x, false)[1] ), IsEmpty) then
            
              Add(orb, neworb);
              Add(tau, newtau);
              Add(rep, newrep);
              
              # For debug
              Add(reppos, [n, pos]);

            fi;
          fi;
        od;
      od;

    until len = Size(orb);
  
  fi;

  # The following relationships hold:
  # 1. Simg = S ^ rep[i];

  # for o in orb do
  #   pos := Position(orb, o);
  #   Simg := AffineCrystGroupOnLeft( Concatenation(o, tgen) );

  #   Print( S ^( rep[pos]^-1 )  = Simg, "\n");
  # od;

  # 2. The normalizer elements in the rep list
  #  correspond to the product of a specific ordered subset of the 
  #  first elements of each item in the reppos list.
  #
  # lst:=List(reppos, Last);
  # for i in [1..Size(lst)] do
  #   pos := i;
  #   ll := [pos];
  #   while lst[pos] <> 1 do
  #     Add(ll, lst[pos]);
  #     pos:= lst[pos];
  #   od;
  #   ll := Reversed(ll);
  #   Print(i, " ",rep[i] = Product(List(reppos, First){ll}),"\n");
  # od;

  res := rec( 
              tau  := tau,
              rep  := rep,
              M     := M 
            );

  return res;

end );

# 和 OrbitSpaceGroupStdByNormalizerPointGroup 的结果进行比较：
InstallGlobalFunction(
OrbitSpaceGroupStdByCollectEquivExtensions, function( S )
  local P, d, norm, Pgen, I, 
        N, F, rels, mat, ext, orbcee, orb,
        Sgen, t, M, pos, x;

  # S:= SpaceGroup(4, 834);
  # S:= SpaceGroupOnRightIT(3,74);

  if not IsStandardAffineCrystGroup( S ) then
    Error("only work with StandardAffineCrystGroup");
  fi;

  # By default, the related functions called here 
  # implemented in Cryst package are designed for matrices 
  # acting on the right 
  if IsAffineCrystGroupOnLeft( S ) then
    S := TransposedMatrixGroup(S);
  fi;

  P := PointGroup(S);
  d := DimensionOfMatrixGroup(S)-1;
  I := IdentityMat(d);
  norm := Filtered(GeneratorsOfGroup( Normalizer( GL(d, Integers), P ) ), x -> not x in P);
 
  if IsTrivial( P ) then
    # d:=3;
    # S := MakeSpaceGroup( d, [], [], false );
    # S := MakeSpaceGroup( d, [], [], transpose );
    # t := Concatenation(List( GeneratorsOfGroup(S), x -> x[1+d]{[1..d]} ));

    # In this algorithm, Sgen is obtianed as follows: 
    # PreImagesRepresentative( PointHomomorphism(S) , IdentityMat(d));
    # 
    Sgen := [ IdentityMat(d + 1) ];
    t := List( Concatenation(List(Sgen, x ->x{[1..d]}[d+1])), FractionModOne );
    orb := [t];
    return orb;
  fi;

  # first get group relators for grp
  N := NiceObject( P );
  F := Image( IsomorphismFpGroupByGenerators( N, GeneratorsOfGroup( N ) ) );
  rels := List( RelatorsOfFpGroup( F ), ExtRepOfObj );

  Pgen := GeneratorsOfGroup( P );

  # construct equations which determine the non-primitive translations
  # an alternative would be
  #  mat := MatJacobianMatrix( F, Pgen );
  mat := GroupExtEquations( d, Pgen, rels );


  # now solve them modulo integers
  ext := SolveHomEquationsModZ( mat ); 
  # 很多情况，轨道尺寸过长，比如 SpaceGroup(4, 834) 的结果如下：
  # 从而造成过滤出结果的相关工作非常耗时，使得这个方法没有实用价值了。
  # [ 4096, 2 ]
  # List(ext, Size); 

  #F  . . . . collect extensions equivalent by conjugation with elems from norm
  # collect group extensions which are equivalent as space groups
  # 下面的结果就是在norm作用下按 orbits 归类的结果：

  #  该函数中的下面的逻辑，会清除 ext[1] 变量的内容，
  #  故该函数无法第二次调用：
  #  while ll<>[] do
  #   SubtractSet( ll, orb );
  #  od; 
  orbcee := CollectEquivExtensions( ext[1], ext[2], norm, P );
  # For debug
  # orbcee := dev1CollectEquivExtensions( ext[1], ext[2], norm, PZ );
  

  # orbnpg:=OrbitSpaceGroupStdByNormalizerPointGroup(S);
  # tau := orbnpg.tau;
  # rep := orbnpg.rep;
  # M := orbnpg.M;

  # # 不完全相同，但是它们之间以纯平移共轭一一对应：
  # for i in Difference( orb, orbcee[2] ) do
  #   for j in Difference( orbcee[2], orb ) do
  #     sol := SolveInhomEquationsModZ(M, i - j, false)[1];
  #     if not IsEmpty(sol) then
  #       Print(i, " ",j,"\n");
  #     fi; 
  #   od;
  # od;

  # 如下，即可从 orbcee 中提取出当前SG在 std 表示下的被norm 共轭的所有可能的 orbits： 
  # 注意和当前的惯例 matrices acting on the left or right 的情况对应：
  # For matrices acting on the right
  # 这里的用法和我的写法是一致的：
  # https://www.math.colostate.edu/~hulpke/GAPQA/qa7.html
  M:=TransposedMat(Concatenation(List( Pgen, x -> TransposedMat(x) - IdentityMat(d) )));

  # # 和 ConjugatorSpaceGroupsStdSamePG 中的算法比较：
  # M1 := List( [1..d], i->[] ); i := 0;

  # for g in Pgen do
  #     g := g - I;
  #     M1{[1..d]}{[1..d]+i*d} := g{[1..d]}{[1..d]};
  #     i := i+1;
  # od;
  # # 相等：
  # M = M1;

  # 另一种方法：
  # List( Pgen, x -> x - IdentityMat(d) );
  # List([1..d], i -> Flat(List( last, x ->  x[i])))=M;

  Sgen:=List( Pgen, x -> PreImagesRepresentative(PointHomomorphism( S ), x ) );  
  t:= List(Concatenation(List(Sgen, x ->x[d+1]{[1..d]})), FractionModOne);
  pos:=First([1..Size(orbcee)], i -> First( orbcee[i], x -> not IsEmpty( SolveInhomEquationsModZ(M, t -x, true)[1] ) ) <> fail );

  orb := orbcee[pos];
  return orb;

end );  


# 这个是进一步深入系统研究crystallography space groups的非常好的工具集：
# 基于carat的进一步研究：
# https://lbfm-rwth.github.io/carat/doc/introduction.html#examples
InstallGlobalFunction( AffineIsomorphismSpaceGroups, function( S1, S2 )

    local d, P1, P2, S1tr, S2tr,
          S1s, S2s, S3s, 
          P1s, P2s, S3sgen, t3s, t, sol, pos,
          c1, C1, c2, C2, c3, C3, C4, C,
          Pgen, norm, orbnpg, tau, rep, M;

    # Affine crystallographic groups vs space groups.
    # https://github.com/gap-packages/cryst/issues/36#issuecomment-1472348928
    # Space groups have a translation subgroup of full rank d, whereas general affine crystallographic groups may have a translation subgroup of smaller rank.
    if not IsSpaceGroup( S1 ) or not IsSpaceGroup( S2 ) then
      Error("S1 and S2 must be space groups");
    fi;

    # We work with SpaceGroupOnLeft
    if IsAffineCrystGroupOnLeft( S1 ) <> IsAffineCrystGroupOnLeft( S2 ) then
      return fail;
    else
      if IsAffineCrystGroupOnRight( S1 ) then
        S1tr := TransposedMatrixGroup( S1 );
        S2tr := TransposedMatrixGroup( S2 );
        C    := AffineIsomorphismSpaceGroups( S1tr, S2tr );
        if C = fail then
            return fail;
        else
            return TransposedMat( C )^-1;
        fi;
      fi;
    fi;

    d := DimensionOfMatrixGroup( S1 ) - 1; 
    P1 := PointGroup( S1 );
    P2 := PointGroup( S2 );

    # If the point groups of S1 and S2 are not isomorphic, then S1 and S2 belong to different Q-classes.
    if Size( P1 ) <> Size( P2 ) or fail = IsomorphismGroups(P1,P2) then
      return fail;
    fi;
 
    # go to standard representation
    # 此处采用我的记号，
    # S1^C1 = S1s
    # For matrices acting on the left, 对应于 cryst 的如下记号：
    # S1^(C1^-1) = S1s
    C1    := ConjugatorReducedSpaceGroup( S1 );
    S1s := S1^(C1^-1);
    P1s := PointGroup( S1s );

    # S2^C2 = S2s
    C2  := ConjugatorReducedSpaceGroup( S2 );
    S2s := S2^(C2^-1);
    P2s := PointGroup( S2s );

    # P1s^c3 = P2s; 
    # Check whether S1 and S2 belong to the same Z-class:
    c3  := RepresentativeAction( GL(d,Integers), P1s, P2s );
    if c3 = fail then
        return fail;
    fi;
    C3 := AugmentedMatrixOnLeft( c3, 0*[1..d] );

    # 按照我的记法：
    # S1s ^ C3 = S3s 
    S3s := S1s ^ (C3^-1);
    
    # 求 S2s 和 S3s 之间的同构关系：
    Pgen := GeneratorsOfGroup(P2s);

    S3sgen:= List(Pgen, x -> PreImagesRepresentative(PointHomomorphism(S3s), x ));
    t3s:= List(Concatenation(List(S3sgen, x ->x{[1..d]}[d+1])), FractionModOne);

    # 汇总所有共轭关系如下：
    # S1^C1 = S1s
    # S2^C2 = S2s
    # S1s ^ C3 = S3s 
    # S2s ^ C4 = S3s
    # S1s ^ C3 = S2s ^ C4 
    # S1s ^ C3 * (C4 ^ -1) = S2s = S2^C2
    # S1 ^ (C1 * C3 * C4 ^ -1 * C2 ^ -1) = S2
    norm := GeneratorsOfGroup(Normalizer(GL(d,Integers), P2s)); 
    orbnpg := OrbitSpaceGroupStdByNormalizerPointGroup(S2s, Pgen, norm);
    tau := orbnpg.tau;
    rep := orbnpg.rep;
    M := orbnpg.M;

    # 验证结果的正确性：
    # S1^(C^-1) = S2;
    # AffineCrystGroupOnLeft(OnTuples( GeneratorsOfGroup(S1), C ))=S2;

    # Check whether S1 and S2 belong to the same space group type:
    for t in tau do
      #  Using fr package:
      # sol := SolutionMatMod1( TransposedMat(M), t - t1 );
      # if sol <> fail then

      # S2s ^ C4 = S3s
      # ( S2s ^ rep[pos]) ^ {E|sol}) = S3s
      sol := SolveInhomEquationsModZ( M, t3s - t, false )[1];
      if not IsEmpty(sol) then
        pos := Position(tau, t);
        # so we have
        # C4:= rep[pos] * {E|sol} 
        # = [[ rep[pos]{[1..d]}{[1..d]}, (0 * [1..d]) ^ T ], [ 0* [1..d], 1 ]] *  [[IdentityMat(d), sol ^ T], [0 * [1..d],1]]
        # = [[ rep[pos]{[1..d]}{[1..d]}, rep[pos]{[1..d]}{[1..d]} * sol ^ T], [ 0 * [1..d],1]]
        # so, we have the following:
        # C4 := rep[pos] * AugmentedMatrixOnLeft( IdentityMat(d), sol[1] );
        # or
        C4 := AugmentedMatrixOnLeft( rep[pos]{[1..d]}{[1..d]}, rep[pos]{[1..d]}{[1..d]} * sol[1] );
        C := C1 * C3 * C4 ^ -1 * C2 ^ -1;
        # Catch the case of a trivial point group in ConjugatorSpaceGroups 
        # https://github.com/gap-packages/cryst/commit/51f53da7de4f1e697d5dc7fa1fc687c1e2e43b23
        # return [C, C1, C2, C3, C4];
        # Print(C); 
        return C;
      fi;
    od;

    # If we have arrived here, it means that S1 and S2 must belong to different space group types.
    # https://github.com/gap-packages/cryst/issues/38#issuecomment-1498458435
    return fail;

end );


# About my research results on the minimal generating set of space groups and my thanks to you.
# https://mail.google.com/mail/u/0/?ogbl#drafts/KtbxLwHDlCCBldrFJCXLMDsMVGZBLcklzL

# Ask for your comments and suggestions about my implementation of the MinimalGeneratingSetAffineCrystGroup.
# https://mail.google.com/mail/u/0/?ogbl#inbox/QgrcJHsHqgRmwLtpXtLbsrWKXXWRXfgncNL
#  Dear Zhao,

# well -- the minimal generating set problem for infinite non-nilpotent groups
# is not only very hard, but in its generality it is algorithmically undecidable.
# As to the minimal generating set problem for pcp groups (or in fact, as to
# pcp groups in general), I think Max is a much more suitable person to talk to
# than me.

# Best wishes,

#     Stefan

# P.S.: As to the finite orders of elements of GL(n,Z) -- of course no algorithm
# can be faster than linear in the output length.

# 下面的算法在寻找 tmgs 方法也存在进一步改进和精细搜索的可能，比如：
# 在找出tmgs后，如果当前的 长度仍然大于 b，则可以在 bound 到 Size( tmgs) - 1 
# 的范围内进行tupgens的更彻底的搜索。
# 但是，考虑到目前的结果合理性，Stefan在上面的评注以及进一步细化实现的不易。
# 似乎并没有进一步处理的必要。
InstallGlobalFunction( MinimalGeneratingSetAffineCrystGroup, function( S )
  local d, ntgens,  iso, epi, H, Hgens, sgs, bound, cmgs, sch, cmb,
        tupgens, tmgs, mgs, res, i, k, s, t, x;

  # S := SpaceGroup(4,4565);
  # S := SpaceGroup(4,4);
  # S:= AffineCrystGroupOnLeft(SGGenSetBC[229]);

  if not IsSolvable(S) then
    Error( "only work with solvable AffineCrystGroup" );
  fi;
  
  # do a few basic checks
  if IsAffineCrystGroupOnRight( S ) then
    S := TransposedMatrixGroup( S );
    res := MinimalGeneratingSetAffineCrystGroup( S );
    res.mgs := List(res.mgs, TransposedMat);
    return res;
  fi;
  
  d := DimensionOfMatrixGroup( S ) - 1;

  iso:=IsomorphismPcpGroup( S );
  H:=Image(iso);
  Hgens:=GeneratorsOfGroup(H);
  sgs:=SmallGeneratingSet(H);

  # 修正 AbelianInvariants 作为 the highest theoretical lower bound 的相关计算：
  # lower bound
  # Try to obtain the MinimalGeneratingSet of an FP group with the clues of the relationship analysis among the relators based on graph theory.
  # https://github.com/gap-packages/grape/issues/45#issuecomment-1345264858
  epi:=MaximalAbelianQuotient(H);
  bound:=Maximum(List(List([ Image(epi), PointGroup(S) ], MinimalGeneratingSet), Size));

  # As long as it is not a cyclic group, at least two generators are required.
  if bound=1 and not IsAbelian(H) then
    bound:=2;
  fi;
  
  sch := []; # The already searched set so far.

  # The key is to choose the appropriate method to construct the generating sets for enumeration. Here's what we do:
  if Size(sgs) > bound then 
    # In this case, Set is much more efficient than Unique:
    cmgs:=Set(Filtered(List( Combinations(Hgens, Size(Hgens) -1 ), x -> Product(x) ), x -> not IsOne(x) ));
     
    for i in [bound..Size(sgs)-1] do
      cmb:=Combinations(cmgs,i);
      Add(sch, [i, cmb]);
      mgs:=First(cmb, x -> Index(H, Subgroup(H, x))=1);

      if mgs <> fail then
        break;
      fi;
    od;
   
    if mgs = fail then
      mgs:=sgs;
    fi;

    # Perform further search to minimize the size of the generating set as much as possible:
    if Size(mgs) > bound then
      tupgens:=Set(Filtered(List( Tuples(Hgens, Size(Hgens) -1 ), x -> Product(x) ), x -> not IsOne(x) ));

      for i in [bound..Size(mgs)-1] do
        # This is a key trick: Construction of an appropriate number of generating sets.

        k:=i;
        s := Last(First(sch, x -> x[1] = i));
        repeat
          k:=k+1;
          # Until one of the following conditions is met:
          # 1. k is already equal to Size(tupgens), in this cae, tupgens will be used entirely to construct the generating sets of length i;
          # 2. The number of generating sets of length i constructed with a subset of k elements selected in tupgens - Size(s)  >= Size(tupgens)，
          # That is, the number of generating sets that have not been enumerated >= Size(tupgens).
        until k = Size(tupgens) or Size( Combinations( [1..k], i ) ) - Size(s) >= Size(sgs) * Size(tupgens); 
        t := Combinations( Shuffle(tupgens){[1 ..k ]}, i );
        # Avoid re-searching generating sets that have already been searched:
        SubtractSet(t, s); 

        tmgs:=First(t, x -> Index(H, Subgroup(H, x))=1);
        
        if tmgs <> fail then 
          mgs:=tmgs;
          break;
        fi;
      od; 
    fi;

  else
    mgs:=sgs;
  fi;

  mgs:=List( mgs, x -> PreImagesRepresentative(iso, x) );

  res := rec( mgs  := mgs,
              bound  := bound );
  return res;

end );


InstallGlobalFunction( 
ConjugatorMatrixGroups, function( G1, G2 )

  local cc1, cc2, ccr2, gen1, gen2, lst, hom1, hom2, m, r, x, c, g;

  if not IsMatrixGroup( G1 ) or 
     not IsMatrixGroup( G2 ) or
     not IsFinite( G1 ) or 
     not IsFinite( G2 ) then
    Error("only work with finite matrix groups" );
  fi;

  # Check if G1 and G2 are isomorphic groups of the same size:
  if Size( G1 ) <> Size( G2 ) or fail = IsomorphismGroups(G1,G2) then
    return fail;
  fi;
  
  # Using the minimal generating set to improve algorithm efficiency.
  gen1 := MinimalGeneratingSet(G1);
  G1:= GroupWithGenerators(gen1);
  # hom1:=GroupHomomorphismByImagesNC(G1, G1);
  # or simply
  hom1:=IdentityMapping(G1);

  # Define polynomial variable
  x := Indeterminate(Rationals, "x");

  # Calculate elementary divisors and conjugacy class sizes for each generator in G1's minimal generating set
  cc1 := ConjugacyClasses(G1);
  lst := List(gen1, g -> [
    ElementaryDivisorsMat(PolynomialRing(Rationals, 1), x * g^0 - g * One(x)),
    Size(First(cc1, c -> g in c))
  ]);

  # Calculate conjugacy classes representatives for G2
  cc2 := ConjugacyClasses(G2);

  ccr2 := Cartesian(List(lst, l -> Filtered(cc2, c -> 
    ElementaryDivisorsMat(PolynomialRing(Rationals, 1), x * Representative(c)^0 - Representative(c) * One(x)) = l[1] and Size(c) = l[2]
  )));

  # Find a set of generators of G2 that ensures LinearRepresentationIsomorphism between G1 and G2
  for r in ccr2 do
    for gen2 in Cartesian(List(r, List)) do
      if Subgroup(G2, gen2) = G2 then
        G2:=GroupWithGenerators(gen2);
        hom2 := GroupHomomorphismByImagesNC(G1, G2);
        # compute m, such that `G1^(m^-1)=G2`
        m:=LinearRepresentationIsomorphism(hom1, hom2 : use_kronecker);
        if fail <> m then
          # Normalize `m`
          m := m/Gcd(Flat(m));
          # Return the inverse of `m`, 
          return m^-1;
        fi;  
      fi;
    od;
  od;

  # If such a isomorphism cannot be found, return fail
  return fail;

end );


InstallGlobalFunction( 
EnantiomorphicPairOfSpaceGroup, function( S )
  local d, trT, P, Pgen, N, norm, diag, C, res,
        A, B, z, Ugen, orbnpg, orbnpg_posi, i, x;

  if IsAffineCrystGroupOnRight( S ) then
    S := TransposedMatrixGroup( S );
    res := EnantiomorphicPairOfSpaceGroup( S );

    if res <> fail then
      res := TransposedMat( res );
    fi;

    return res;
  fi;

  if not IsStandardAffineCrystGroup( S ) then
    S := StandardAffineCrystGroup(S);
  fi;

  d:= DimensionOfMatrixGroup(S) - 1;
  trT := TransposedMat(TranslationBasis( S ));
  P:=PointGroup(S);
  Pgen:=GeneratorsOfGroup(P);
  
  N:=Normalizer(GL(d,Integers), P);
  norm :=Set(Filtered(GeneratorsOfGroup(N), x -> not x in P));
  
  # 便于后续构造 Schreier generators:
  # https://en.wikipedia.org/wiki/Schreier%27s_lemma
  # https://github.com/gap-packages/cryst/issues/23#issuecomment-844364463
  AddSet(norm, IdentityMat(d));

  # ref AffineNormalizer
  # For right action:
  # we work in a standard representation
  # if not IsStandardAffineCrystGroup( S ) then
  #     invT := T^-1;
  #     gens := List( GeneratorsOfGroup( N ), x -> T * x * invT );
  #     Pgens := List( Pgens, x -> T * x * invT );
  #     Sgens := List( Sgens, x -> S!.lconj * x * S!.rconj );
  # else
  #     gens := GeneratorsOfGroup( N );
  # fi;
 
  diag:=List( [1..d+1], i -> 1 );
  diag[1]:=-1;
  
  C:=fail;
 
  if ForAll(Pgen, x -> DeterminantMat(x)=1) then
    if ForAll(norm, x -> DeterminantMat(x)=1) then
      # Print("case 1: ", i, "\n");
      # The conjugator to get the enantiomorphic partner
      C:=DiagonalMat(diag);
    else
  
      A:=Filtered(norm, x -> DeterminantMat(x)=1);
      B:=Filtered(norm, x -> DeterminantMat(x)=-1);
      z:=First(B);
      Ugen:=Union(A, z * B, B/z, z * A/z );
     
      orbnpg:=OrbitSpaceGroupStdByNormalizerPointGroup(S, Pgen, norm);
      orbnpg_posi:=OrbitSpaceGroupStdByNormalizerPointGroup(S, Pgen, Ugen);
      
      # The following is enough for quick check:
      if Size(orbnpg.tau)/2=Size(orbnpg_posi.tau) then
        # Print("case 2: ", i, "\n");
        # The conjugator to get the enantiomorphic partner
        C:=AugmentedMatrixOnLeft(z, 0 * [1..d]);
      fi;

      # 进一步彻底验证轨道分裂关系：
      # if Size( orbnpg.rep) >=2 and ForAny(orbnpg.rep, x -> DeterminantMat(x) = -1) then

      #   pos:=First([1..Size(orbnpg.rep)], i -> DeterminantMat(orbnpg.rep[i]) = -1);
      #   tau:=orbnpg.tau[pos];
   
      #   P1gen:=List( [1..Size(Pgen)], i -> AugmentedMatrixOnLeft(Pgen[i],tau{[1 + d*(i-1)..i* d]}) );
              
      #   S_nega:=AffineCrystGroupOnLeft(Concatenation(P1gen, tgen));

      #   orbnpg_nega:=OrbitSpaceGroupStdByNormalizerPointGroup(S_nega, Pgen, Ugen);
            
      #   M:=orbnpg.M;
      #   orb:=Difference(Set(Concatenation(orbnpg_posi.tau, orbnpg_nega.tau)), orbnpg.tau);
      #   # if Size(orbnpg.tau)/2=Size(orbnpg_posi.tau) and ForAll(orb, x -> fail <> First(orbnpg.tau, y -> not IsEmpty(SolveInhomEquationsModZ(M, x -y,false)[1]) ) ) then
      #     Print("case 2: ", i, "\n");
      #     Add(lst3,i);
      #   fi;

      # fi;
    fi; 
  fi;

  if C = fail then 
    return C;
  else
    return C^(AugmentedMatrixOnLeft(trT, 0*[1..d])^-1);
  fi;

end );


# Add the function DirectSumDecompositionMatrix.
# https://github.com/gap-packages/utils/issues/64#issuecomment-1591422710
# Direct sum decomposition of a matrix.
# https://mail.google.com/mail/u/0/?ogbl#sent/KtbxLwHHpRfqHzrslftpfCbXhvQJngJvxV
# Dear forum,

# This is a possible approach, inspired in
# https://github.com/gap-packages/numericalsgps/blob/fcde379b01bd44b1fa80cd69d7ddd6a8acdcfe2f/gap/catenary-tame.gi#LL803C1-L831C4

InstallGlobalFunction( 
DirectSumDecompositionMatrix, function(l)
   
  local nr,nc,i,j,nzr,nzb,rest,bls,nbls;

	if not(IsMatrix(l))  then
		Error("The argument must be a matrix.");
	fi;
  rest:=StructuralCopy(l);
  nr:=Length(rest);
  nc:=Length(rest[1]);
  bls:=[];
  i:=1;
  j:=1;
  while nr>0 and nc>0 do
      nzr:=true;
      nzb:=true;
      while nzr or nzb do
          nzr:=ForAny([1..i],i1->ForAny([j+1..nc],j1-> not IsZero(rest[i1][j1])));
          nzb:=ForAny([i+1..nr],i1->ForAny([1..j],j1-> not IsZero(rest[i1][j1])));
          if nzr then
              j:=j+1;
          fi;
          if nzb then
              i:=i+1;
          fi;
      od;
      Add(bls,List([1..i],i1->rest[i1]{[1..j]}));
      rest:=List([i+1..nr],i1->rest[i1]{[j+1..nc]});
      nr:=nr-i;
      nc:=nc-j;
      i:=1; j:=1;
  od;

  nbls:=Length(bls);

  # check if we have filled all columns and rows
  if nc>0 and nbls>0 then # add zeroes at the end of the last block
      bls[nbls]:=List(bls[nbls], l->Concatenation(l,ListWithIdenticalEntries(nc,0)));
  fi;

  if nr>0 and nbls>0 then # add zero rows at the end of the last block
      bls[nbls]:=Concatenation(bls[nbls], ListWithIdenticalEntries(nr,ListWithIdenticalEntries(Length(bls[nbls]), 0)
      ));
  fi;

  return bls;

end );


#############################################################################
##
#F  LLLTranslationBasis( S ) . . . . . determine basis of translation lattice 
# using LLLReducedBasis 
##

# 1. 似乎应该首先转到标准表示，再进行这里的简化处理。
# 2. 采用LLLReducedBasis的方法不太好处理 lllrb.transformation 这部分。
# 故改为基于 基于LLLReducedGramMat的方法
# InstallGlobalFunction( 
# LLLTranslationBasis, function ( S )

#     local d, P, Sgens, Pgens, trans, g, m, F, Fgens, rel, new,
#           lllrb, N, T1, T2, T;

#     if IsAffineCrystGroupOnRight( S ) then
#       T := LLLTranslationBasis( TransposedMatrixGroup( S ) );
#       return T;
#     fi;

#     d := DimensionOfMatrixGroup( S ) - 1;
#     P := PointGroup( S );
#     Pgens := [];
#     Sgens := [];
#     trans := [];

#     # first the obvious translations
#     for g in GeneratorsOfGroup( S ) do
#         m := g{[1..d]}{[1..d]};
#         if IsOne( m ) then
#             Add( trans, g{[1..d]}[d+1] );
#         else
#             Add( Sgens, g );
#             Add( Pgens, m );
#         fi;
#     od;

#     # then the hidden translations
#     if not IsTrivial( P ) then
#         F := Image( IsomorphismFpGroupByGenerators( P, Pgens ) );
#         Fgens := GeneratorsOfGroup( FreeGroupOfFpGroup( F ) );
#         for rel in RelatorsOfFpGroup( F ) do
#             new := MappedWord( rel, Fgens, Sgens );
#             Add( trans, new{[1..d]}[d+1] );
#         od;
#     fi;
    
#     # 似乎不用下面的处理，就可以保证最后的T能转到标准表示：
#     # make translations invariant under point group
#     # trans := Set( Union( Orbits( TransposedMatrixGroup(P), trans ) ) );

#     # return ReducedLatticeBasis( trans );

#     # 按列矢量形式的对应的格基变换：
#     # basis_lattice * TransposedMat(lllrb.transformation ) = TransposedMat(lllrb.basis);
#     lllrb:=LLLReducedBasis( trans, "linearcomb" );
#     T1:=TransposedMat(Filtered(TransposedMat(lllrb.transformation), x -> not IsZero(x)));
    
#     # 即使不采用orbit 计算的处理，得到的T1 仍旧可能大于d维。
#     # 但是，下面的处理的到底好不好，也不知道。
#     # if Last(DimensionsMat(T1))>d then
#     #   N := LinearIndependentColumns( T1 );
#     #   T1 :=TransposedMat(TransposedMat(T1){N});
#     # fi;

#     # TransposedMat(lllrb.transformation)*trans=lllrb.basis;
#     T2:=lllrb.basis;
#     T:=TransposedMat(T1^-1*T2);

#     # 确保 LLLTranslationBasis 可以转到标准表示。
#     return T * TransposedMat(TranslationBasis(S^(AugmentedMatrixOnLeft(T, 0 *[1..d])^-1)));

# end );


InstallGlobalFunction( 
LLLTranslationBasis, function ( S )

    local d, P, F, llg, T1, T2, T;

    if IsAffineCrystGroupOnRight( S ) then
      T := LLLTranslationBasis( TransposedMatrixGroup( S ) );
      return T;
    fi;
    
    d := DimensionOfMatrixGroup( S ) - 1;
    # note down the original base.
    T1:= TransposedMat(InternalBasis(S));
    # Then switch to the standard representation.
    S := StandardAffineCrystGroup(S);
    P := PointGroup( S );

    F:=Sum(List(P, g-> TransposedMat(g) * g ));
    llg:=LLLReducedGramMat(F);
    T2:=TransposedMat(llg.transformation);
    
    T:=T1*T2;

    # Ensure the transformation matrix obtained can be used to converted to standard representation.
    # return T * TransposedMat(TranslationBasis(S^(AugmentedMatrixOnLeft(T, 0 *[1..d])^-1)));
    return T;

end );



# About the three classes translations related to a specific space group.
# https://mail.google.com/mail/u/0/?ogbl#search/branton%40byu.edu+origin+shift+/QgrcJHsbjCgGxkTcpwdpcRTMdWjmWTPHncg

# Some further discussions based on your lecture note "Group theory applied to crystallography" and some of your papers.
# https://mail.google.com/mail/u/0/?ogbl#sent/KtbxLvHgMkZlGBgrVFcHgjClrgRRwWDdtg

# Bernd Souvignier的 lecture note，page 27， Theorem 43:
# https://www.math.ru.nl/~souvi/krist_09/cryst.pdf
# The following theorem (which is not hard to prove) states that by an appropriate shift of the
# origin, the coordinates of a SNoT become rational numbers with denominators at most the order
# |P| of the point group.


# 由此，可以分如下步骤，简化空间群的表示:
# 1. 基于 LLLTranslationBasis 整数化点群（线性）部分，
# 2. 基于 Theorem 43，有理化矢量系统，

# 用下面的解决方法来首先彻底简化已给空间群的表示：
InstallGlobalFunction( 
ConjugatorReducedSpaceGroup, function( S )
  local d, C, P, hom, trans, v, x;
  
  # 这样处理是否合适？
  if IsAffineCrystGroupOnRight( S ) then
    C := ConjugatorReducedSpaceGroup( TransposedMatrixGroup( S ) );
    return TransposedMat(C) ^ -1;
  fi;

  d:=DimensionOfMatrixGroup(S) - 1;
  C:=AugmentedMatrixOnLeft(LLLTranslationBasis(S), 0*[1..d]);

  S:=S^(C^-1);
  P:=PointGroup(S);

  if not ForAll(Flat(List(GeneratorsOfGroup(S), x ->x{[1..d]}[d+1]) ), IsRat) then
    hom:=PointHomomorphism(S);
    trans:=List(P, x -> PreImagesRepresentative(hom, x) );
    v:=Sum(List(trans, x -> x{[1..d]}[d+1]))/Size(trans);

    C:=C * AugmentedMatrixOnLeft(IdentityMat(d), v);

  fi;

  return C;

end );

