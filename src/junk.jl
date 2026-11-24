# not fully translated from GAP3
"""
RelsViaCosetTable(G,cosets)  . . . . . . . construct relators for the
RelsViaCosetTable(G,cosets,ggens)  . . . . . . . . . given concrete
RelsViaCosetTable(G,cosets,F,words,H,F1)  . . . . . . . group

`RelsViaCosetTable`  constructs a defining set of relators  for the given
concrete group using John Cannon's relations finding algrotithm.

It is a  subroutine  of function  `PresentationViaCosetTable`.  Hence its
input  and  output  are   specifically  designed  for  this  purpose.  In
particular, it does not check the arguments.

G,            # given group
cosets,       # right cosets of G with respect to H
"""
function RelsViaCosetTable(G,cosets,arg...)
  #  F,        # given free group
  #  words,    # given words for the generators of H
  #  H,        # subgroup, if specified
  #  F1,       # f.p. group isomorphic to H
  #  F2,       # f.p. group isomorphic to G
  #  ng1,      # position number of identity element in G
  #  nh1,      # position number of identity element in H
  #  perms,    # permutations induced by the gens on the cosets
  #  stage,    # 1 or 2
  #  table,    # columns in the table for gens
  #  rels,     # representatives of the relators
  #  relsGen,  # relators sorted by start generator
  #  subgroup, # rows for the subgroup gens
  #  rels1,    # list of relators
  #  app,      # arguments list for `MakeConsequences`
  #  index,    # index of the table
  #  col,      # generator col in auxiliary table
  #  perm,     # permutations induced by a generator on the cosets
  #  gens,     # abstract gens in which the relators are written
  #  gens2,    # the above abstract gens and their inverses
  #  ggens,    # concrete generators of G
  #  ngens,    # number of generators of G
  #  ngens2,   # twice the above number
  #  order,    # order of a generator
  #  actcos,   # part 1 of Schreier vector of G by H
  #  actgen,   # part 2 of Schreier vector of G by H
  #  tab0,     # auxiliary table in parallel to table <table>
  #  cosRange, # range from 1 to index (= number of cosets)
  #  genRange, # range of the odd integers from 1 to 2*ngens-1
  #  geners,   # order in which the table cols are worked off
  #  next,     # local coset number
  #  left1,    # part 1 of Schreier vector of H by trivial group
  #  right1,   # part 2 of Schreier vector of H by trivial group
  #  n,        # number of subgroup element
  #  words2,   # words for the generators of H and their inverses
  #  h         # subgroup element
  if length(arg)==1 ggens=arg[1]
  else ggens=G[:generators]
  end
  ngens=length(ggens)
  ngens2=ngens * 2
  if cosets[1] in G ng1=PositionSorted(cosets, cosets[1] ^ 0)
  else ng1=1
  end
  index=length(cosets)
  tab0=[]
  table=[]
  subgroup=[]
  cosRange=1:index
  genRange=map((i->begin 2i-1 end), 1:ngens)
  if length(arg)<2
      stage=1
      F2=FreeGroup(ngens)
      rels=[]
  else
      stage=2
      F=arg[1]
      words=arg[2]
      F2=Group(F[:generators], IdWord)
      if haskey(F, :namesGenerators)
          F2[:namesGenerators]=F[:namesGenerators]
      end
      H=arg[3]
      nh1=PositionSorted(H[:elements], (H[:elements])[1] ^ 0)
      F1=arg[4]
      rels=map(rel->MappedWord(rel, F1[:generators], words), F1[:relators])
      left1=F1[:actcos]
      right1=F1[:actgen]
      words2=[]
      for i in 1:length(F1[:generators])
          push!(words2, words[i])
          push!(words2, words[i] ^ -1)
      end
  end
  gens=F2[:generators]
  gens2=[]
  perms=map((gen->begin Permutation(gen, cosets, OnRight) end), ggens)
  for i in 1:ngens
      push!(gens2, gens[i])
      push!(gens2, gens[i] ^ -1)
      perm=perms[i]
      col=OnTuples(cosRange, perm)
      gen=0*cosRange
      push!(tab0, col)
      push!(table, gen)
      order=Order(G, ggens[i])
      if order==2 push!(rels, gens[i] ^ 2)
      else col=OnTuples(cosRange, perm ^ -1)
          gen=0*cosRange
      end
      push!(tab0, col)
      push!(table, gen)
  end
  cosets=0*cosRange
  actcos=0*cosRange
  actgen=0*cosRange
  cosets[1]=ng1
  actcos[ng1]=ng1
  j=1
  i=0
  while i<index
      i+=1
      c=cosets[i]
      g=0
      while g<ngens2
          g+=1
          next=(tab0[g])[c]
          if next>0 && actcos[next]==0
              g1=(g+2 * mod(g, 2))-1
              table[g][c]=next
              table[g1][next]=c
              tab0[g][c]=0
              tab0[g1][next]=0
              actcos[next]=c
              actgen[next]=g
              j=j+1
              cosets[j]=next
              if j==index
                  g=ngens2
                  i=index
              end
          end
      end
  end
  rels=RelatorRepresentatives(rels)
  app=fill(0,11)
  app[1]=table
  app[5]=subgroup
  if stage==2
      relsGen=RelsSortedByStartGen(gens, rels, table)
      app[4]=relsGen
      for g in genRange
          gen0=tab0[g]
          for c=cosRange
              if gen0[c]==0
                  app[10]=g
                  app[11]=c
                  n=MakeConsequences(app)
              end
          end
      end
  end
  geners=1:ngens2
  for i in cosets
      for j=geners
          if table[j][i]==0
              g=(j+2 * mod(j, 2))-1
              c=tab0[j][i]
              table[j][i]=c
              table[g][c]=i
              tab0[j][i]=0
              tab0[g][c]=0
              rel=IdWord
              while c != ng1
                  g=actgen[c]
                  rel=rel * gens2[g] ^ -1
                  c=actcos[c]
              end
              rel=rel ^ -1 * gens2[j] ^ -1
              c=i
              while c != ng1
                  g=actgen[c]
                  rel=rel * gens2[g] ^ -1
                  c=actcos[c]
              end
              if stage==2
                  h=MappedWord(rel, gens, ggens)
                  n=PositionSorted(H[:elements], h)
                  while n != nh1
                      g=right1[n]
                      rel=rel * words2[g] ^ -1
                      n=left1[n]
                  end
              end
              rels1=RelatorRepresentatives([rel])
              if length(rels1)>0
                  rel=rels1[1]
                  if !rel in rels
                      push!(rels, rel)
                  end
              end
              relsGen=RelsSortedByStartGen(gens, rels, table)
              app[4]=relsGen
              for g in genRange
                  gen=table[g]
                  gen0=tab0[g]
                  inv0=tab0[g+1]
                  for c=cosRange
                      if gen[c]>0
                          gen0[c]=0
                          inv0[gen[c]]=0
                      end
                  end
              end
              for g in genRange
                  gen0=tab0[g]
                  for c=cosRange
                      if gen0[c]==0
                          app[10]=g
                          app[11]=c
                          n=MakeConsequences(app)
                      end
                  end
              end
          end
      end
  end
  F2=F2 // rels
  F2[:actcos]=actcos
  F2[:actgen]=actgen
  return F2
end

function Ignore end
if !@isdefined InfoFpGroup1
  InfoFpGroup1=Ignore
end
if !@isdefined InfoFpGroup2
  InfoFpGroup2=Ignore
end

"""
`PresentationViaCosetTable( G[, F, words] )`

`PresentationViaCosetTable`  constructs a presentation  for the given group
`G`.  The method  being used  is John  Cannon's relations finding algorithm
which has been described in cite{Can73} or in cite{Neu82}.

In  its first form,  if only the  group `G` has  been specified, it applies
Cannon's  single stage  algorithm which,  by plain  element multiplication,
computes a coset table of `G` with respect to its trivial subgroup and then
uses coset enumeration methods to find a defining set of relators for `G`.

|    gap> G := GeneralLinearGroup( 2, 7 );
    GL(2,7)
    gap> G.generators;
    [ [ [ Z(7), 0*Z(7) ], [ 0*Z(7), Z(7)^0 ] ],
      [ [ Z(7)^3, Z(7)^0 ], [ Z(7)^3, 0*Z(7) ] ] ]
    gap> Size( G );
    2016
    gap> P := PresentationViaCosetTable( G );
    << presentation with 2 gens and 5 rels of total length 46 >>
    gap> relators( P );
    &I  1. f.2^3
    &I  2. f.1^6
    &I  3. f.1*f.2*f.1*f.2*f.1*f.2*f.1*f.2*f.1*f.2*f.1*f.2
    &I  4. f.1*f.2*f.1^-1*f.2*f.1*f.2^-1*f.1^-1*f.2*f.1*f.2*f.1^-1*f.2^-1
    &I  5. f.1^2*f.2*f.1*f.2*f.1*f.2^-1*f.1^-1*f.2^-1*f.1^3*f.2^-1 |
    gap> DisplayPresentation(P);
    1: b³=1
    2: a³=a⁻³
    3: bababa=a⁻¹b⁻¹a⁻¹b⁻¹a⁻¹b⁻¹
    4: bab⁻¹a⁻¹ba=ab⁻¹a⁻¹bab⁻¹
    5: a²babab⁻¹a⁻¹b⁻¹a³b⁻¹=1 |
The  second form  allows to  call Cannon's  two stage algorithm which first
applies  the single stage  algorithm to an  appropriate subgroup `H` of `G`
and  then uses the resulting relators of `H`  and a coset table of `G` with
respect  to `H` to find relators of  `G`. In this case the second argument,
`F`,  is assumed to be  a free group with  the same number of generators as
`G`, and `words` is expected to be a list of words in the generators of `F`
which, when being evaluated in the corresponding generators of `G`, provide
subgroup generators for `H`.

|    gap> M12 := MathieuGroup( 12 );;
    gap> M12.generators;
    [ ( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11), ( 3, 7,11, 8)( 4,10, 5, 6),
      ( 1,12)( 2,11)( 3, 6)( 4, 8)( 5, 9)( 7,10) ]
    gap> F := FreeGroup( "a", "b", "c" );
    Group( a, b, c )
    gap> words := [ F.1, F.2 ];
    [ a, b ]
    gap> P := PresentationViaCosetTable( M12, F, words );
    << presentation with 3 gens and 10 rels of total length 97 >>
    gap> DisplayPresentation(P);
    1: c=c⁻¹
    2: b²=b⁻²
    3: cac=a⁻¹c⁻¹a⁻¹
    4: ab⁻²ab⁻²ab⁻²=1
    5: a¹¹=1
    6: a²b⁻¹a²b=ba⁻¹b²a²
    7: ba⁻¹ba⁻¹b⁻¹a=a⁻¹bab⁻¹ab⁻¹
    8: a²ba²b⁻²a⁻¹ba⁻¹b⁻¹a⁻¹b⁻¹=1
    9: a²b⁻¹a⁻¹b⁻¹acbcabab=1
    10: a³ba²ba⁻²caba⁻¹ca=1
    gap> G := FpGroup(P);
    Group( a, b, c )
    gap> G.relators;
    [ c^2, b^4, a*c*a*c*a*c, a*b^-2*a*b^-2*a*b^-2, a^11,
      a^2*b*a^-2*b^-2*a*b^-1*a^2*b^-1,
      a*b*a^-1*b*a^-1*b^-1*a*b*a^-1*b*a^-1*b^-1,
      a^2*b*a^2*b^-2*a^-1*b*a^-1*b^-1*a^-1*b^-1,
      a^2*b^-1*a^-1*b^-1*a*c*b*c*a*b*a*b, a^3*b*a^2*b*a^-2*c*a*b*a^-1*c*a
     ] |

Before  it is returned,  the resulting presentation  is being simplified by
appropriate  calls  of  the  function  `simplify`  (see "Tietze
Transformations"),  but without  allowing it  to eliminate  any generators.
This  restriction guarantees  that we  get a  bijection between the list of
generators of `G` and the list of generators in the presentation. Hence, if
the  generators  of  `G`  are  redundant  and  if  you  don't  care for the
bijection, it may be convenient to apply the function
`simplify` again.

julia> P=Presentation("
1: a³=1
2: a⁻¹f⁻¹b⁻¹=1
3: b³=1
4: c³=1
5: c⁻¹e⁻¹d=1
6: d³=1
7: af=ba
8: ac=ca
9: bc=cb
10: ad=da
11: bd=db
12: dc=ce⁻¹
")
Presentation: 6 generators, 12 relators of total length 42
    gap> simplify( P );
    &I  there are 4 generators and 10 relators of total length 36 |

PresentationViaCosetTable(G[,F,words]) . . . . .
construct a presentation for the given group

`PresentationViaCosetTable`   constructs   a  presentation  for  a  given
concrete  group.  It applies  John  Cannon's  relations finding algorithm
which has been described in

  John J. Cannon:  Construction of  defining relators  for finite groups.
  Discrete Math. 5 (1973), 105-129, and in

  Joachim Neubueser: An elementary introduction to coset table methods in
  computational  group theory.  Groups-St. Andrews 1981,  edited by C. M.
  Campbell and E. F. Robertson, pp. 1-45.  London Math. Soc. Lecture Note
  Series no. 71, Cambridge Univ. Press, Cambridge, 1982.

If only a group  `G`  has been  specified,  the single stage algorithm is
applied.

If the  two stage algorithm  is to  be used,  `PresentationViaCosetTable`
expects a subgroup `H` of `G` to be described by two additional arguments
`F`  and  `words`,  where  `F`  is a  free group  with the same number of
generators as  `G`,  and  `words` is a list of words in the generators of
`F`  which supply  a list of generators of  `H`  if they are evaluated as
words in the corresponding generators of `G`.
"""
function PresentationViaCosetTable(G,arg...)
  #  F,          # given free group
  #  fgens,      # generators of F
  #  fwords,     # given words in the generators of F
  #  words,      # tidied up words for the subgroup generators
  #  H,          # subgroup
  #  elts,       # elements of G or H
  #  cosets,     # right cosets of G with respect to H
  #  F1,         # f.p. group isomorphic to H
  #  F2,         # f.p. group isomorphic to G
  #  P,          # resulting presentation
  #  ggens,      # concrete generators of G
  #  hgens,      # concrete generators of H
  #  ngens,      # number of generators of G
  if !(IsRec(G) && (haskey(G, :isGroup) && G[:isGroup]))
      error("first argument must be a group")
  end
  ggens=G[:generators]
  ngens=length(ggens)
  if length(arg)==0
    InfoFpGroup1("#I  calling the single stage relations finding algorithm\n")
    elts=elements(G)
    F2=RelsViaCosetTable(G, elts)
  else
    F=arg[1]
    if !(IsRec(F) && (haskey(F, :isFpGroup) && F[:isFpGroup]))
        error("second argument must be a free group")
    end
    fgens=F[:generators]
    if length(fgens) != ngens
        error("the given groups have different number of generators")
    end
    fwords=arg[2]
    hgens=map(w->MappedWord(w, fgens, ggens), fwords)
    H=Subgroup(G, hgens)
    if Size(H)==1 || Size(H)==Size(G)
      InfoFpGroup1("#I  calling the single stage relations finding algorithm\n")
      elts=elements(G)
      F2=RelsViaCosetTable(G, elts)
    else
      InfoFpGroup1("#I  ","calling the two stage relations finding algorithm\n",
                   "#I  using a subgroup of size ",Size(H),"  &&  index ",
                   Size(G) // Size(H), "\n")
      words=[]
      hgens=[]
      for w=fwords
          h=MappedWord(w, fgens, ggens)
          if h != G.identity && (!h in hgens && !(h ^ -1) in hgens)
              push!(words, w)
              push!(hgens, h)
          end
      end
      hgens=map(w->MappedWord(w, fgens, ggens), words)
      elts=elements(H)
      F1=RelsViaCosetTable(H, elts, hgens)
      cosets=RightCosets(G, H)
      F2=RelsViaCosetTable(G, cosets, F, words, H, F1)
    end
  end
  P=Presentation(F2, 0)
  P.generatorsLimit=ngens
  GoGo(P)
  P.generatorsLimit=0
  P.debug=1
  return P
end
