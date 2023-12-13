
<a id='GroupPresentations'></a>

<a id='GroupPresentations-1'></a>

# GroupPresentations

- [`GroupPresentations`](index.md#GroupPresentations)
- [`GroupPresentations.AbsWord`](index.md#GroupPresentations.AbsWord)
- [`GroupPresentations.FpGroup`](index.md#GroupPresentations.FpGroup)
- [`GroupPresentations.Presentation`](index.md#GroupPresentations.Presentation-Tuple{FpGroup})
- [`GroupPresentations.conjugate`](index.md#GroupPresentations.conjugate)
- [`GroupPresentations.relators`](index.md#GroupPresentations.relators)
- [`GroupPresentations.showgens`](index.md#GroupPresentations.showgens)
- [`GroupPresentations.simplify`](index.md#GroupPresentations.simplify)
- [`GroupPresentations.tryconjugate`](index.md#GroupPresentations.tryconjugate)
- [`GroupPresentations.@AbsWord`](index.md#GroupPresentations.@AbsWord)

<a id='GroupPresentations' href='#GroupPresentations'>#</a>
**`GroupPresentations`** &mdash; *Module*.



This  is a  port of  some GAP3/VkCurve  functionality on *presentations* of *finitely presented groups*.

We  have defined just enough functionality  on finitely presented groups so that  presentations  can  be  translated  to  finitely presented groups and vice-versa. The focus is on presentations, the goal being to simplify them.

The  elements of finitely presented groups are `AbsWord` or abstract words, representing elements of a free group. 

```julia-repl
julia> @AbsWord a,b,c,d,e,f # same as a=AbsWord(:a);b=AbsWord(:b)...

julia> F=FpGroup([a,b,c,d,e,f])
FreeGroup(a,b,c,d,e,f)

julia> G=F/[a^2,b^2,d*f^-1,e^2,f^2,a*b^-1*c,a*e*c^-1,b*d^-1*c,c*d*e^-1,a*f*c^-2,c^4]
FreeGroup(a,b,c,d,e,f)/[a²,b²,df⁻¹,e²,f²,ab⁻¹c,aec⁻¹,bd⁻¹c,cde⁻¹,afc⁻²,c⁴]

julia> simplify(F) # the main function of this package
Presentation: 2 generators, 4 relators, total length 16
Presentation: 2 generators, 3 relators, total length 10
FreeGroup(a,c)/[a²,ac⁻¹ac⁻¹,c⁴]
```

The simplification is done by the following process:

```julia-rep1
julia> P=Presentation(G);simplify(P);G=FpGroup(P)
```

The  functions `Presentation`  and `FpGroup`  create a  presentation from a finitely presented group and vice versa.

In order to speed up the algorithms, the relators in a presentation are not represented  internally by `AbsWord`s, but by lists of positive or negative generator  numbers which  we call  *Tietze words*.  Here is another example with a few functions to explore presentations.

```julia-repl
julia> @AbsWord a,b 

julia> F=FpGroup([a,b])
FreeGroup(a,b)

julia> G=F/[a^2,b^7,comm(a,a^b),comm(a,a^(b^2))*inv(b^a)]
FreeGroup(a,b)/[a²,b⁷,a⁻¹b⁻¹a⁻¹bab⁻¹ab,a⁻¹b⁻²a⁻¹b²ab⁻²ab²a⁻¹b⁻¹a]

julia> P=Presentation(G) # by default give a summary
Presentation: 2 generators, 4 relators, total length 30

julia> relators(P)
4-element Vector{AbsWord}:
 a²
 b⁷
 ab⁻¹abab⁻¹ab
 b⁻²ab²ab⁻²ab²ab⁻¹
```

```julia-rep1
julia> showgens(P)
1. a 10 occurrences involution
2. b 20 occurrences

julia> dump(P) # here in relators inverses are represented by capitalizing
# F relator
1:3 aa
2:0 bbbbbbb
3:0 aBabaBab
4:0 abbaBBabbaBBB
gens=AbsWord[a, b] involutions:AbsWord[a] modified=false numredunds=0

julia> display_balanced(P)
1: a=A
2: bbbbbbb=1
3: aBab=BAbA
4: BBabbaBBabbaB=1

julia> P=tryconjugate(P) # try to conjugate the generators
Presentation: 2 generators, 4 relators, total length 30
Bab=> Presentation: 2 generators, 3 relators, total length 28
# Bab gives Presentation: 2 generators, 3 relators, total length 28
Presentation: 2 generators, 3 relators, total length 28

julia> FpGroup(P) # slightly simplified group
FreeGroup(a,b)/[b⁷,bab⁻¹abab⁻¹a,b⁻¹ab²ab⁻²ab²ab⁻²]
```

for  more  information  look  at  the  help  strings  of `AbsWord, FpGroup, Presentation,     relators,    display_balanced,    simplify,    conjugate, tryconjugate`. 

A  minimal thing to add to this package so it would be a reasonable package for finitely preented groups is the Coxeter-Todd algorithm.


<a target='_blank' href='https://github.com/jmichel7/GroupPresentations.jl/blob/585b41e9aa1ce7d458ce958ba86918348a8f1d5f/src/GroupPresentations.jl#L1-L96' class='documenter-source'>source</a><br>

<a id='GroupPresentations.AbsWord' href='#GroupPresentations.AbsWord'>#</a>
**`GroupPresentations.AbsWord`** &mdash; *Type*.



An  `AbsWord` represents an  element of the  free group on some generators. The  generators  are  indexed  by  `Symbols`.  For  example  the  `Absword` representing `a³b⁻²a` is represented internally as `[:a => 3, :b => -2, :a => 1]`. The mulitiplcation follows the group rule:

```julia-repl
julia> w=AbsWord([:a => 3, :b => -2, :a => 1])
a³b⁻²a

julia> w*AbsWord([:a=>-1,:b=>1])
a³b⁻¹

```

A positive `AbsWord` may be obtained by giving `Symbols` as arguments

```julia-repl
julia> AbsWord(:b,:a,:a,:b)
ba²b
```


<a target='_blank' href='https://github.com/jmichel7/GroupPresentations.jl/blob/585b41e9aa1ce7d458ce958ba86918348a8f1d5f/src/GroupPresentations.jl#L183-L201' class='documenter-source'>source</a><br>

<a id='GroupPresentations.@AbsWord' href='#GroupPresentations.@AbsWord'>#</a>
**`GroupPresentations.@AbsWord`** &mdash; *Macro*.



`@AbsWord x,y` is the same as `x=AbsWord(:x);y=AbsWord(y)`


<a target='_blank' href='https://github.com/jmichel7/GroupPresentations.jl/blob/585b41e9aa1ce7d458ce958ba86918348a8f1d5f/src/GroupPresentations.jl#L225' class='documenter-source'>source</a><br>

<a id='GroupPresentations.FpGroup' href='#GroupPresentations.FpGroup'>#</a>
**`GroupPresentations.FpGroup`** &mdash; *Type*.



`FpGroup(P::Presentation)`

returns the finitely presented group defined  by the presentation `P`.


<a target='_blank' href='https://github.com/jmichel7/GroupPresentations.jl/blob/585b41e9aa1ce7d458ce958ba86918348a8f1d5f/src/GroupPresentations.jl#L381-L385' class='documenter-source'>source</a><br>

<a id='GroupPresentations.Presentation-Tuple{FpGroup}' href='#GroupPresentations.Presentation-Tuple{FpGroup}'>#</a>
**`GroupPresentations.Presentation`** &mdash; *Method*.



`Presentation( G::FpGroup[, debug=1])`

returns  the  presentation  corresponding  to  the given finitely presented group `G`.

The  optional `debug` parameter  can be used  to restrict or  to extend the amount  of output  provided by  Tietze transformation  functions when being applied  to the created  presentation. The default  value 1 is designed for interactive  use and implies  explicit messages to  be displayed by most of these functions. A `debug` value of 0 will suppress these messages, whereas a `debug` value of 2 will enforce some additional output.


<a target='_blank' href='https://github.com/jmichel7/GroupPresentations.jl/blob/585b41e9aa1ce7d458ce958ba86918348a8f1d5f/src/GroupPresentations.jl#L431-L443' class='documenter-source'>source</a><br>

<a id='GroupPresentations.relators' href='#GroupPresentations.relators'>#</a>
**`GroupPresentations.relators`** &mdash; *Function*.



`relators(P::Presentation)` relators of `P` as `AbsWord`s. 


<a target='_blank' href='https://github.com/jmichel7/GroupPresentations.jl/blob/585b41e9aa1ce7d458ce958ba86918348a8f1d5f/src/GroupPresentations.jl#L378' class='documenter-source'>source</a><br>

<a id='GroupPresentations.showgens' href='#GroupPresentations.showgens'>#</a>
**`GroupPresentations.showgens`** &mdash; *Function*.



`showgens(P,list=eachindex(P.generators))`

prints  the generators of `P` with the total number of their occurrences in the  relators, and notes involutions. A  second `list` argument prints only those generators.


<a target='_blank' href='https://github.com/jmichel7/GroupPresentations.jl/blob/585b41e9aa1ce7d458ce958ba86918348a8f1d5f/src/GroupPresentations.jl#L566-L572' class='documenter-source'>source</a><br>

<a id='GroupPresentations.simplify' href='#GroupPresentations.simplify'>#</a>
**`GroupPresentations.simplify`** &mdash; *Function*.



`simplify(G)`

`simplify`  applies  Tietze  transformations  to a  copy of  the presentation of the given finitely presented group `G` in order to reduce it with respect to  the number of generators, the number of relators, and the relator lengths.

`simplify` returns the resulting finitely presented group (which is isomorphic to `G`).

```julia-repl
julia> @AbsWord a,b,c,d,e,f

julia> F=FpGroup([a,b,c,d,e,f])
FreeGroup(a,b,c,d,e,f)

julia> G=F/[a^2,b^2,d*f^-1,e^2,f^2,a*b^-1*c,a*e*c^-1,b*d^-1*c,c*d*e^-1,a*f*c^-2,c^4]
FreeGroup(a,b,c,d,e,f)/[a²,b²,df⁻¹,e²,f²,ab⁻¹c,aec⁻¹,bd⁻¹c,cde⁻¹,afc⁻²,c⁴]

julia> simplify(G)
FreeGroup(a,c)/[a²,ac⁻¹ac⁻¹,c⁴]
```

In fact, the call

```julia-rep1
julia> simplify(G)
```

is an abbreviation of the call sequence

```julia-rep1
julia> P=Presentation(G,0);simplify(P);FpGroup(P)
```

which applies  a rather simple-minded strategy of  Tietze transformations to the intermediate presentation `P`. If for  some  concrete group the resulting presentation  is unsatisfying, then  you  should  try  a  more  sophisticated,  interactive  use of  the available Tietze transformation functions  (see "Tietze Transformations").


<a target='_blank' href='https://github.com/jmichel7/GroupPresentations.jl/blob/585b41e9aa1ce7d458ce958ba86918348a8f1d5f/src/GroupPresentations.jl#L1238-L1279' class='documenter-source'>source</a><br>


`simplify(p [,tries])`

simplify  the  presentation  `p`.  We  have  found heuristics which make it somewhat  efficient, but the algorithm depends  on random numbers so is not reproducible.  The main  idea is  to rotate  relators between  calls to the basic  `GroupPresentations.Go` function. By default 100 such rotations are tried (unless  the  presentation  is  so  small  that  less rotations exhaust all possible  ones), but the actual number tried  can be controlled by giving a second  parameter `tries` to the function. Another useful tool to deal with presentations is `tryconjugate`.

```julia-rep1
julia> display_balanced(p)
1: ab=ba
2: dbd=bdb
3: bcb=cbc
4: cac=aca
5: adca=cadc
6: dcdc=cdcd
7: adad=dada
8: Dbdcbd=cDbdcb
9: adcDad=dcDadc
10: dcdadc=adcdad
11: dcabdcbda=adbcbadcb
12: caCbdcbad=bdcbadBcb
13: cbDadcbad=bDadcbadc
14: cdAbCadBc=bdcAbCdBa
15: cdCbdcabdc=bdcbadcdaD
16: DDBcccbdcAb=cAbCdcBCddc
17: CdaBdbAdcbCad=abdcAbDadBCbb
18: bdbcabdcAADAdBDa=cbadcbDadcBDABDb
19: CbdbadcDbbdCbDDadcBCDAdBCDbdaDCDbdcbadcBCDAdBCDBBdacDbdccb=abdbcabdcAdcbCDDBCDABDABDbbdcbDadcbCDAdBCabDACbdBadcaDbAdd

julia> simplify(p)
Presentation: 4 generators, 18 relators, total length 304
Presentation: 4 generators, 18 relators, total length 284
Presentation: 4 generators, 17 relators, total length 264
Presentation: 4 generators, 16 relators, total length 256
Presentation: 4 generators, 15 relators, total length 244
Presentation: 4 generators, 15 relators, total length 240
Presentation: 4 generators, 15 relators, total length 226
Presentation: 4 generators, 15 relators, total length 196
Presentation: 4 generators, 15 relators, total length 178
Presentation: 4 generators, 15 relators, total length 172
Presentation: 4 generators, 14 relators, total length 158

julia> display_balanced(p)
1: ab=ba
2: dbd=bdb
3: bcb=cbc
4: cac=aca
5: adAc=cadA
6: dcdc=cdcd
7: adad=dada
8: CdBcbd=bCdBcb
9: adcDad=dcDadc
10: dcdadc=adcdad
11: cbdcbdc=dcbdcbd
12: dcbadcbda=adcbcadcb
13: cbCDadcab=DadcbadcD
14: caDCbdBcADbda=bDBaDbADcbadc
```


<a target='_blank' href='https://github.com/jmichel7/GroupPresentations.jl/blob/585b41e9aa1ce7d458ce958ba86918348a8f1d5f/src/GroupPresentations.jl#L3146-L3209' class='documenter-source'>source</a><br>

<a id='GroupPresentations.conjugate' href='#GroupPresentations.conjugate'>#</a>
**`GroupPresentations.conjugate`** &mdash; *Function*.



`conjugate(p,conjugation)`

This program modifies a presentation by conjugating a generator by another. The  conjugation to  apply is  described by  a length-3  string of the same style  as  the  result  of  `display_balanced`,  that is `"abA"` means replace  the second generator by its  conjugate by the first, and  `"Aba"` means replace it by its conjugate by the inverse of the first.

```julia-rep1
julia> display_balanced(P)
1: dabcd=abcda
2: dabcdb=cabcda
3: bcdabcd=dabcdbc

julia> display_balanced(conjugate(P,"Cdc"))
<< presentation with 4 generators, 3 relators of total length 36>>
1: dcabdc=cabdca
2: abdcab=cabdca
3: bdcabd=cabdca
```


<a target='_blank' href='https://github.com/jmichel7/GroupPresentations.jl/blob/585b41e9aa1ce7d458ce958ba86918348a8f1d5f/src/GroupPresentations.jl#L3248-L3269' class='documenter-source'>source</a><br>

<a id='GroupPresentations.tryconjugate' href='#GroupPresentations.tryconjugate'>#</a>
**`GroupPresentations.tryconjugate`** &mdash; *Function*.



`tryconjugate(p[,goal[,printlevel]])`

This program tries to simplify group presentations by applying conjugations to  the  generators.  The  algorithm  depends  on  random  numbers,  and on tree-searching,  so is  not reproducible.  By default  the program stops as soon  as a shorter presentation is found.  Sometimes this does not give the desired  presentation.  One  can  give  a  second argument `goal`, then the program  will only stop when  a presentation of length  less than `goal` is found.  Finally, a third  argument can be  given and then all presentations the  programs runs  over which  are of  length less  than or  equal to this argument are displayed. Due to the non-deterministic nature of the program, it  may be useful to  run it several times  on the same input. Upon failure (to improve the presentation), the program returns `p`.

```julia-rep1
julia> display_balanced(p)
1: ba=ab
2: dbd=bdb
3: cac=aca
4: bcb=cbc
5: dAca=Acad
6: dcdc=cdcd
7: adad=dada
8: dcDbdc=bdcbdB
9: dcdadc=adcdad
10: adcDad=dcDadc
11: BcccbdcAb=dcbACdddc
julia> p=tryconjugate(p)
Presentation: 4 generators, 11 relators, total length 100
dcD=> Presentation: 4 generators, 10 relators, total length 90
# dcD gives Presentation: 4 generators, 10 relators, total length 90
Presentation: 4 generators, 10 relators, total length 90

julia> p=tryconjugate(p)
Dcd=> Presentation: 4 generators, 10 relators, total length 88
# Dcd gives Presentation: 4 generators, 10 relators, total length 88
Presentation: 4 generators, 10 relators, total length 88

julia> p=tryconjugate(p)
dcD=> Presentation: 4 generators, 10 relators, total length 90
Dbd=> Presentation: 4 generators, 10 relators, total length 96
Aca=> Presentation: 4 generators, 9 relators, total length 84
Presentation: 4 generators, 8 relators, total length 76
# Aca gives Presentation: 4 generators, 8 relators, total length 76
Presentation: 4 generators, 8 relators, total length 76

julia> p=tryconjugate(p)
Bcb=> Presentation: 4 generators, 8 relators, total length 70
# Bcb gives Presentation: 4 generators, 8 relators, total length 70
Presentation: 4 generators, 8 relators, total length 70

julia> p=tryconjugate(p)
Cac=> Presentation: 4 generators, 8 relators, total length 64
# Cac gives Presentation: 4 generators, 8 relators, total length 64
Presentation: 4 generators, 8 relators, total length 64

julia> p=tryconjugate(p)
caC=> Presentation: 4 generators, 8 relators, total length 58
# caC gives Presentation: 4 generators, 8 relators, total length 58
Presentation: 4 generators, 8 relators, total length 58

julia> p=tryconjugate(p)
Cac=> Presentation: 4 generators, 8 relators, total length 64
Cbc=> Presentation: 4 generators, 7 relators, total length 50
# Cbc gives Presentation: 4 generators, 7 relators, total length 50
Presentation: 4 generators, 7 relators, total length 50

julia> p=tryconjugate(p)
cdC=> Presentation: 4 generators, 7 relators, total length 56
Dcd=> Presentation: 4 generators, 7 relators, total length 54
Cac=> Presentation: 4 generators, 7 relators, total length 48
# Cac gives Presentation: 4 generators, 7 relators, total length 48
Presentation: 4 generators, 7 relators, total length 48

julia> p=tryconjugate(p)
caC=> Presentation: 4 generators, 7 relators, total length 50
Cdc=> Presentation: 4 generators, 7 relators, total length 50
Dbd=> Dcd=> Presentation: 4 generators, 7 relators, total length 60
Bab=> Aba=> Aca=> Presentation: 4 generators, 7 relators, total length 46
# Aca gives Presentation: 4 generators, 7 relators, total length 46
Presentation: 4 generators, 7 relators, total length 46

julia> display_balanced(p)
1: db=bd
2: ba=ab
3: cac=aca
4: ada=dad
5: bcb=cbc
6: cdcd=dcdc
7: AdCacd=cAdCac
```


<a target='_blank' href='https://github.com/jmichel7/GroupPresentations.jl/blob/585b41e9aa1ce7d458ce958ba86918348a8f1d5f/src/GroupPresentations.jl#L3289-L3381' class='documenter-source'>source</a><br>

