# https://www.cryst.ehu.es/cgi-bin/cryst/programs/checkgr.pl?tipog=gesp
# https://iso.byu.edu/iso/findssghelp.php
# The output notation for superspace-group operators will match that of the input: There are three choices: (x,y,z,t,u,v), (x1,x2,x3,x4,x5,x6), and (xs1,xs2,xs3,xs4,xs5,xs6). See the ISO(3+d)D help page for more information about these notations. 

# Define variable-length `Indeterminate's in function.
# https://mail.google.com/mail/u/0/?ogbl#search/in%3Asent+gap/QgrcJHsBqxpXbHCSnPqBGGBZRDqphmxvpZb
InstallGlobalFunction( IdentifyGroupGenerators, function( gens )
  
  local d, vecname, vec, g;

  if not ForAll(gens, IsAffineMatrixOnLeft) then
    Error("AffineMatrixOnLeft test failed");
  fi;

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
  
  local Str, Sgen, out,
        dir, shell, d, cmd, program, cs, args, 
        _in, _out, res, err, x; 

  # get temporary file names
  Sgen := CaratTmpFile( "Sgen" );  
  out := CaratTmpFile( "out"  );
  
  d := DimensionOfMatrixGroup( S ) - 1;

  # `Name' is based on the CARAT database of all Q-classes of finite unimodular groups of degree up to 6. 
  if not IsSpaceGroup(S) or d > 6 then
    Error( "only works for space groups of degree up to 6" );
  fi;

  # https://docs.gap-system.org/pkg/caratinterface/htm/CHAP003.htm#SECT001
  # In crystallography, the convention usually is that matrix groups act from the left on column vectors. This convention is adopted also in CARAT.
  if IsAffineCrystGroupOnRight( S ) then
    Str := TransposedMatrixGroup( S );
    return CaratName( Str );
  fi;

  # compatiable with the result returned by `Name' in this case:
  if IsTrivial( PointGroup( S ) ) then
    return fail;
  fi;

  # `Name' is not able to perform calculations directly using proper cyclotomics. Therefore, a conversion process is carried out beforehand:
  if not IsStandardSpaceGroup( S ) then
    S := StandardAffineCrystGroup( S ); 
  fi;
  CaratWriteMatrixFile(Sgen, GeneratorsOfGroup( S ));

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

  d := DimensionOfMatrixGroup( S ) - 1;
  # CrystCatRecord(TransposedMatrixGroup(S)).parameters;
  
  if d > 6 or not IsSpaceGroup( S ) then
    Error("only applicable to space groups in dimensions up to 6");
  fi;

  if IsAffineCrystGroupOnRight( S ) then
    S := TransposedMatrixGroup( S );
    # 此时的进一步递归处理：
    res := IdentifySpaceGroup( S );
    res[2] := TransposedMat(res[2]) ^-1;
    return res;
  fi;

  name:=CaratName( S );

  if d = 2 then 
    nr:=Position(cryst2names, name);
    Sref:=SpaceGroupOnLeftIT(d,nr);
    C := AffineIsomorphismSpaceGroups(S, Sref);
  elif d = 3 then
      
    # cryst3names 为最精细分类，故aff3names没有必要使用：
    # Positions(aff3names, name );
    nrs:=Positions(cryst3names, name );
    
    # determine from the enantiomorphic_pairs
    if Size(nrs) = 1 or 1 <> DeterminantMat(AffineIsomorphismSpaceGroups(S,SpaceGroupOnLeftIT(d,nrs[2]))) then
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

# 基于 CRT 的算法实现：
# Orbits: These SGs are affine isomorphic (but not equivalent/identical) by a conjugator whose linear part is a element of NormalizerPointGroupInGLnZ
# and the transation part is zero.


# 轨道长度必须小于CRT约束的点群元的 order：
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


# 因为 cryst 的相关bug已经修复，故可以不再使用fr包：
# https://github.com/gap-packages/cryst/pull/33#issuecomment-1427933014
# https://github.com/gap-packages/cryst/commit/e761436db386107f3c0b60a927bbc352b3220d58
# LoadPackage("fr");
InstallGlobalFunction( OrbitCrystStdByNormalizerPointGroup, function( S )
  local P, d, Pgen, Sgen, NPgen, NP, I, M, hom, t, t1, t2, sol,
        orb, orbs, norms, nelm, cnt, sch, conv, threshold, maxord, 
        ord, n, N, NSgen, iso, F, Fgen, res, g, i, j;

  # Call the PackageName variable from with the corresponding package.        
  # https://mail.google.com/mail/u/0/?ogbl#sent/KtbxLxGcCbLxdbmDRqglCFpbRFhGQNpKLq
  if not IsStandardAffineCrystGroup( S ) then
    Error("only work with StandardAffineCrystGroup");
  fi;

  if IsAffineCrystGroupOnRight( S ) then
    S := TransposedMatrixGroup( S );
    res := OrbitCrystStdByNormalizerPointGroup( S );
    res.M := TransposedMat( res.M );
    res.norms := List(res.norms, TransposedMat);
    return res;          
  fi;

  d := DimensionOfMatrixGroup(S) - 1;
  I := IdentityMat(d);

  # 基于 cryst 包提供的一些函数：
  hom := PointHomomorphism(S);
  P:= PointGroup( S );

  # Why can't there be something like `TrivialGroup( IsMatrixGroup )` in GAP?
  # https://mail.google.com/mail/u/0/?ogbl#sent/KtbxLrjNdcnZznrRkHcNjCnCsKqBGBMczg
  # TrivialGroup( IsMatrixGroup ); 

  # 这里的处理，也解决了下面的问题：
  # Catch the case of a trivial point group in ConjugatorSpaceGroups 
  # https://github.com/gap-packages/cryst/commit/51f53da7de4f1e697d5dc7fa1fc687c1e2e43b23
  # 如果两个标准表示的SG 的点群都为 IsTrivial, 则因为格基都为： IdentityMat(d)，
  # 它们之间的变换必为 IdentityMat(d + 1)。
  # About the result of IsSpaceGroup on some edge cases.
  # https://github.com/gap-packages/cryst/issues/35
  if IsTrivial(P) then
    # P:=Group(IdentityMat(d));
    P:=Group( One(GL(d, Integers)) );
  fi;

  Pgen := GeneratorsOfGroup( P );
  Sgen := List( Pgen, x -> PreImagesRepresentative( hom, x ) );
  t1:= List(Concatenation(List(Sgen, x ->x{[1..d]}[d+1])), FractionModOne);
  
  # In trivial case, sol is always 0 * [1..d]:
  # b:=RandomInvertibleMat(Size(t1))[1];
  # # or
  # b:=RandomUnimodularMat(Size(t1):domain:=[-1000..1000])[1];
  # # M * sol = b  (mod Z)
  # sol:=SolveInhomEquationsModZ(M, b, false)[1];

  
  # 基于 MappedWord 方法(效率也差不多)：
  # Sgen := Filtered( Sgen, x -> not IsOne( x{[1..d]}{[1..d]} ) );
  # Pgen := List( Sgen, x -> x{[1..d]}{[1..d]} );
  # t1:= List(Concatenation(List(Sgen, x ->x{[1..d]}[d+1])), FractionModOne);

  # P:= Group(Pgen, IdentityMat(d) );
  # iso:=IsomorphismFpGroupByGeneratorsNC(P, Pgen); 
  # F := Image(iso);
  # Fgen := FreeGeneratorsOfFpGroup(F);
  # # or
  # # Fgen := GeneratorsOfGroup( FreeGroupOfFpGroup( F ) );
  # # 用法：
  # # List(List( Pgen, g -> g^-1 ), x -> UnderlyingElement(x^iso));
  # # List(last, x -> MappedWord(x, Fgen, Sgen)); 

  M:= Concatenation( List( Pgen, g -> g - I ) );
  
  NPgen:=Filtered(GeneratorsOfGroup( Normalizer( GL(d, Integers), P ) ), x -> not x in P);

  orbs:=[t1];
  norms:=[I]; 

  # catch the trivial cases
  if IsTrivial(P) or IsEmpty(NPgen) then
    res := rec( 
                orbs  := orbs,
                norms  := norms,
                M     := M 
              );
    return res;
  fi; 
   
  NP:=Group(NPgen);

  # 因为要保序，以便for枚举，故不使用Set进行处理：
  # PG Normalizer elements used for enumeration so far:
  nelm:=[I];
  cnt := 0;
  sch := []; # The valid searched subset of nelm.

  if IsOddInt(d) then 
    maxord := Maximum(NewFiniteOrdersOfGLNZ(d+1));
  else 
    maxord := Maximum(NewFiniteOrdersOfGLNZ(d));
  fi; 

  for i in nelm do

    # 1. cnt 的数值可以用来适当控制外层循环的次数：根据当前观察到的 conv 在内层循环
    # 后不再变化的次数的阈值，来决定是否退出外层循环。
    # 2. 若整个群已经枚举完成，则满足 Size(nelm) = Order(NP)。
    # 3. 需要注意的是：这里的 orbs 可能是由不同的 norm 共轭SG 得到的，但是 CRT 中的 order
    # 则对应于用此处的同一个 norm / conjuator 去连续变换表示（change-of-basis，变换格基）情况下的共轭结果平移部分的变化周期的长度的约束。 
    # 因此， 在不是由同一个 norm / conjuator 去连续变换表示而产生全部 orbs 的情况下，
    # Size(orbs) 和 CRT 约束的最高有限阶之间并没有固定的关系。
    # 使用本算法对 4 维空间群进行的测试，印证了上述结论。

    # convergence criteria and thresholds
    conv := [Size(orbs), Size(sch), Size(nelm)];
    threshold := maxord * Size(NPgen);

    # Print( cnt," ", threshold, " ", maxord," ", Size(NPgen), " ", Size(nelm) ," ",Size(orbs)," ", Size(sch),"\n");

    # Order(NP) = Size(nelm) is equivalent to:
    # Order(NP) <> infinity and Order(NP) = Size(nelm)
    # The following condition is more expensive:
    # if Order(NP) = Size(nelm) or (Order(NP) = infinity and cnt >= threshold)  then 
    # The following condition is sufficient for obtaining a complete orbit
    if cnt >= threshold or Order(NP) = Size(nelm)  then 
      res := rec( 
                  orbs  := orbs,
                  norms  := norms,
                  M     := M 
                );
      return res;
      # break;
    fi;

    for j in NPgen do
      g := i * j;

      if not IsOne(g) and not g in nelm then
        Add(nelm,g);
        
        # t2:=List(Concatenation( List( List( OnTuples(Pgen, g^-1), x -> PreImagesRepresentative( hom, x )^AugmentedMatrixOnLeft(g, 0*[1..d]) ), x -> x{[1..d]}[d+1] ) ), FractionModOne); 
        
        # Use MappedWord method:
        # t1:=List(Concatenation( List(
        # List( List( OnTuples(Pgen, g^-1), x -> UnderlyingElement( x ^ iso) ), x ->  MappedWord(x, Fgen, Sgen)^AugmentedMatrixOnLeft(g, 0*[1..d]) ), x -> x{[1..d]}[d+1] ) ), FractionModOne);

        # 我所采用的形式对应如下关系： S1 ^ conj = S2  
        # sol:= SolveInhomEquationsModZ( M, t2-t1, false)[1];

        # 共轭作用对应的表示惯例： 
        # {R, t1}^{E,x} = {R, t2}
        # (R-I) x + t1 = t2 
        # M x + t1 =t2 
        # M x = t2 -t1 

        ord:=1;
        orb:=[];
        while true do 
      
          n := g^ord; 

          # In case g has finite order
          if IsOne(n) then break; fi;

          # 另：t2的求法，参考下面相关部分的注释，
          # 在写文章时，可以作为相关理论的详细描述。
          # 道理上都是相同或类似的。只是在描述和步骤上有些差别而已。
          # devCollectEquivExtensions
          # dev1CollectEquivExtensions
           
          # In case g has infinite order, based on the following theorems:
          # https://en.wikipedia.org/wiki/Crystallographic_restriction_theorem
          # Bieberbach's Theorem II
          # https://en.wikipedia.org/wiki/Space_group#Bieberbach's_theorems

          # The following algorithm is derived:
          ord := ord + 1;  
          N := AugmentedMatrixOnLeft(n, 0*[1..d]);
          NSgen:= List( OnTuples(Pgen, n^-1), x -> PreImagesRepresentative( hom, x )^N );
          
          # Use MappedWord method:
          # NSgen:= List( List( OnTuples(Pgen, n^-1), x -> UnderlyingElement( x ^ iso) ), x ->  MappedWord(x, Fgen, Sgen)^N );
          t2:=List(Concatenation(List(NSgen, x ->x{[1..d]}[d+1])), FractionModOne);
          if not IsEmpty(orb) and t2 = First(orb) then break; fi;
          
          Add(orb,t2);
          
          # 一个非纯平移 conjugator 的幂可能等价于一个纯平移 conjugator；
          # 在 update orbs 列表前，进一步排除此可能。
          # 有点类似于 Gram–Schmidt process 过程中的处理。

          # t1 和 t2 所对应的SG之间可以通过纯平移共轭同构，其中包括了它们相等的情况（有平凡（零）解）。
          # 因此单独使用 ForAll 也是可以的，但是基于第一个条件可以提高效率，避免不必要的计算：
          if not t2 in orbs and ForAll(List(orbs, t -> SolveInhomEquationsModZ( M, t2-t, false)[1] ), IsEmpty) then
          # or using fr package
          # if not t2 in orbs and ForAll(List(orbs, t -> SolutionMatMod1(TransposedMat(M), t2-t) ), x -> x = fail) then

            Add(norms, n);
            Add(orbs, t2);
            AddSet(sch, g);
          fi;

        od; 
      fi;
    od;

    if conv{[1,2]} = [Size(orbs), Size(sch)] and conv[3] < Size(nelm) then
      cnt := cnt + 1;
    fi;

  od;

end );

# 和 OrbitCrystStdByNormalizerPointGroup 的结果进行比较：
InstallGlobalFunction(
OrbitCrystStdByCollectEquivExtensions, function( S )
  local P, d, norm, Pgen, I, 
        N, F, rels, mat, ext, oscee, orbs,
        Sgen, t, M, pos, x;

  # S:= SpaceGroup(4, 834);
  # S:= SpaceGroupOnRightIT(3,74);

  if not IsStandardAffineCrystGroup( S ) then
    Error("only work with StandardAffineCrystGroup");
  fi;

  # By default, the related functions called here implemented in
  # Cryst package are designed for matrices acting on the right 
  # 
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
    orbs := [t];
    return orbs;
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
  oscee := CollectEquivExtensions( ext[1], ext[2], norm, P );
  # For debug
  # oscee := dev1CollectEquivExtensions( ext[1], ext[2], norm, PZ );
  

  # osnpg:=OrbitCrystStdByNormalizerPointGroup(S);
  # M := osnpg.M;
  # orb := osnpg.orbs;
  # norm := osnpg.norms;

  # # 不完全相同，但是它们之间以纯平移共轭一一对应：
  # for i in Difference( orb, oscee[2] ) do
  #   for j in Difference( oscee[2], orb ) do
  #     sol := SolveInhomEquationsModZ(M, i - j, false)[1];
  #     if not IsEmpty(sol) then
  #       Print(i, " ",j,"\n");
  #     fi; 
  #   od;
  # od;

  # 如下，即可从 oscee 中提取出当前SG在 std 表示下的被norm 共轭的所有可能的 orbits： 
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
  pos:=First([1..Size(oscee)], i -> First( oscee[i], x -> not IsEmpty( SolveInhomEquationsModZ(M, t -x, true)[1] ) ) <> fail );

  orbs := oscee[pos];
  return orbs;

end );  


# 这个是进一步深入系统研究crystallography space groups的非常好的工具集：
# 基于carat的进一步研究：
# https://lbfm-rwth.github.io/carat/doc/introduction.html#examples
InstallGlobalFunction( AffineIsomorphismSpaceGroups, function( S1, S2 )

    local d, P1, P2, ls1, ls2, S1tr, S2tr,
          S1s, S2s, S3s, 
          P1s, P2s, S3sgen, t3s, t, sol, pos,
          c1, C1, c2, C2,  c3, C3, C4, C,
          osnpg, M, orb, norm;

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

    # some short cuts
    if Size( P1 ) <> Size( P2 ) then
        return fail;
    fi;
    ls1 := AsSortedList( List( ConjugacyClasses( P1 ),
             x -> [ Size(x), TraceMat( Representative(x) ),
                    DeterminantMat( Representative(x) ) ] ) );
    ls2 := AsSortedList( List( ConjugacyClasses( P2 ),
             x -> [ Size(x), TraceMat( Representative(x) ),
                    DeterminantMat( Representative(x) ) ] ) );
    if ls1 <> ls2 then
        return fail;
    fi;

    # go to standard representation
    # 此处采用我的记号，
    # S1^C1 = S1s
    # For matrices acting on the left, 对应于 cryst 的如下记号：
    # S1^(C1^-1) = S1s
    S1s := StandardAffineCrystGroup(S1);
    P1s := PointGroup( S1s );
    c1 := TransposedMat(InternalBasis( S1 ));
    C1    := AugmentedMatrixOnLeft( c1, 0*[1..d] );

    # S2^C2 = S2s
    S2s := StandardAffineCrystGroup(S2);
    P2s := PointGroup( S2s );
    c2 := TransposedMat(InternalBasis( S2 ));
    C2    := AugmentedMatrixOnLeft( c2, 0*[1..d] );

    # P1s^c3 = P2s; 
    c3  := RepresentativeAction( GL(d,Integers), P1s, P2s );
    if c3 = fail then
        return fail;
    fi;
    C3 := AugmentedMatrixOnLeft( c3, 0*[1..d] );

    # 按照我的记法：
    # S1s ^ C3 = S3s 
    S3s := S1s ^ (C3^-1);
    
    # 求 S2s 和 S3s 之间的同构关系：
    S3sgen:= List( GeneratorsOfGroup(P2s), x -> PreImagesRepresentative(PointHomomorphism(S3s), x ) );
    t3s:= List(Concatenation(List(S3sgen, x ->x{[1..d]}[d+1])), FractionModOne);

    # 汇总所有共轭关系如下：
    # S1^C1 = S1s
    # S2^C2 = S2s
    # S1s ^ C3 = S3s 
    # S2s ^ C4 = S3s
    # S1s ^ C3 = S2s ^ C4 
    # S1s ^ C3 * (C4 ^ -1) = S2s = S2^C2
    # S1 ^ (C1 * C3 * C4 ^ -1 * C2 ^ -1) = S2
    osnpg := OrbitCrystStdByNormalizerPointGroup( S2s );
    M := osnpg.M;
    orb := osnpg.orbs;
    norm := osnpg.norms;

    # 验证结果的正确性：
    # S1^(C^-1) = S2;
    # AffineCrystGroupOnLeft(OnTuples( GeneratorsOfGroup(S1), C ))=S2;
    for t in orb do
      #  Using fr package:
      # sol := SolutionMatMod1( TransposedMat(M), t - t1 );
      # if sol <> fail then

      sol := SolveInhomEquationsModZ( M, t3s - t, false )[1];
      if not IsEmpty(sol) then
        pos := Position(orb, t);
        # S2s ^ C4 = S3s
        # The following relationship holds for matrices acting on the left
        C4 := AugmentedMatrixOnLeft( norm[pos], norm[pos] * sol[1] );
        C := C1 * C3 * C4 ^ -1 * C2 ^ -1;
        # Catch the case of a trivial point group in ConjugatorSpaceGroups 
        # https://github.com/gap-packages/cryst/commit/51f53da7de4f1e697d5dc7fa1fc687c1e2e43b23
        # return [C, C1, C2, C3, C4];
        return C;
      fi;
    od;

    # The incorrect logic implemented in the ConjugatorSpaceGroupsStdSamePG and ConjugatorSpaceGroups.
    # https://github.com/gap-packages/cryst/issues/38
    # return fail;

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
# 但是，考虑到目前的结果合理性， Stefan在上面的评注以及进一步细化实现的不易。似乎并没有进一步处理的必要。
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