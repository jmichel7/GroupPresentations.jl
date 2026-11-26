"""
This  is a  port of  some GAP3/VkCurve  functionality on *presentations* of
*finitely presented groups*.

We  have defined just enough functionality  on finitely presented groups so
that  presentations  can  be  translated  to  finitely presented groups and
vice-versa. The focus is on presentations, the aim is to simplify them.

This package depends on `PermGroups`.

The elements of `FpGroup`s (an abbreviation for "finitely presented groups"
are  [`AbsWord`](@ref), an abbreviation  for "abstract words", representing
elements of a free group.

```julia-repl
julia> @AbsWord a,b,c,d,e,f # same as a=AbsWord(:a);b=AbsWord(:b)...

julia> F=FpGroup([a,b,c,d,e,f])
FreeGroup(a,b,c,d,e,f)

julia> G=F/[a^2,b^2,d*f^-1,e^2,f^2,a*b^-1*c,a*e*c^-1,b*d^-1*c,c*d*e^-1,a*f*c^-2,c^4]
FreeGroup(a,b,c,d,e,f)/[a²,b²,df⁻,e²,f²,ab⁻c,aec⁻,bd⁻c,cde⁻,afc⁻²,c⁴]
```

```julia-rep1
julia> simplify(G) # the main function of this package
Presentation: 2 generators, 4 relators of total length 16
Presentation: 2 generators, 3 relators of total length 10
FreeGroup(a,c)/[a²,ac⁻ac⁻,c⁴]
```
Note that when displaying an `AbsWord`, an exponent ⁻¹ is printed ⁻ to save
space.

The simplification is done by the following process:

```julia-rep1
julia> P=Presentation(G);simplify(P);G=FpGroup(P)
```
The   functions  [`Presentation`](@ref)  and   [`FpGroup`](@ref)  create  a
presentation from a finitely presented group and vice versa.

In order to speed up the algorithms, the relators in a presentation are not
represented  internally by `AbsWord`s, but by lists of positive or negative
generator  numbers which  we call  *Tietze words*.  Here is another example
with a few functions to explore presentations.

```julia-repl
julia> @AbsWord a,b

julia> F=FpGroup([a,b])
FreeGroup(a,b)

julia> G=F/[a^2,b^7,comm(a,a^b),comm(a,a^(b^2))*inv(b^a)]
FreeGroup(a,b)/[a²,b⁷,a⁻b⁻a⁻bab⁻ab,a⁻b⁻²a⁻b²ab⁻²ab²a⁻b⁻a]

julia> P=Presentation(G) # default printing gives a summary
Presentation: 2 generators, 4 relators of total length 30

julia> relators(P)
4-element Vector{AbsWord}:
 a²
 b⁷
 ab⁻abab⁻ab
 b⁻²ab²ab⁻²ab²ab⁻
```

```julia-rep1
julia> showgens(P)
1. a 10 occurrences involution
2. b 20 occurrences

julia> display_balanced(P)
1: a=a⁻
2: b⁷=1
3: ab⁻ab=b⁻a⁻ba⁻
4: b⁻²ab²ab⁻²ab²ab⁻=1

julia> dump(P)
# F relator
1:3 a²
2:1 b⁷
3:1 ab⁻abab⁻ab
4:1 b⁻²ab²ab⁻²ab²ab⁻
gens=a,b involutions:a modified=false numredunds=0
```
for  more  information  look  at the docstrings of
[`AbsWord`](@ref),
[`FpGroup`](@ref),
[`Presentation`](@ref),
[`relators`](@ref),
[`display_balanced`](@ref),
[`simplify`](@ref),
[`conjugate`](@ref),
[`tryconjugate`](@ref).

Finally, the Todd-Coxeter algorithm computes, if finite, a `PermGroup` from
an  `FpGroup`. The permutations  represent the multiplication  table by the
generators so they are on a large set of integers.
```julia-repl
julia> @AbsWord a,b

julia> G=FpGroup([a,b])/[a^2,b^2,(a*b)^3]
FreeGroup(a,b)/[a²,b²,ababab]

julia> PermGroup(G)
Group((1,5)(2,4)(3,6),(1,2)(3,5)(4,6))
```
"""
module GroupPresentations
## Changing Presentations
#
#The  functions `add_generator`, `add_relator`, `remove_relator` can be used to
#change a presentation. In general, they will change the isomorphism type of
#the  group defined  by the  presentation, hence,  though they are sometimes
#used  as subroutines by Tietze transformations functions like `Substitute`,
#they do *not* perform Tietze transformations themselves.
#
## Tietze Transformations
#
#The  functions described in this section can be used to modify a
#group presentation by Tietze transformations.
#
#In  general, the aim of such modifications  will be to *simplify* the given
#presentation,  i.e., to reduce  the number of  generators and the number of
#relators  without increasing too much the  sum of all relator lengths which
#we  will  call  the  *total  length*  of the presentation. Depending on the
#concrete presentation under investigation one may end up with a nice, short
#presentation   or  with  a   very  huge  one.
#
#There  is no way  to find the  shortest presentation which  can be obtained
#from  a given  one. Therefore,  what we  offer are  some lower-level Tietze
#transformation   functions  and,  in  addition,  a  heuristic  higher-level
#function   (which  of  course   cannot  be  the   optimal  choice  for  all
#presentations).
#
#The design of these functions follows closely the concept of the ANU Tietze
#transformation  program designed by George Havas cite{Hav69} which has been
#available  from Canberra since 1977 in a stand-alone version implemented by
#Peter   Kenne  and  James  Richardson  and   later  on  revised  by  Edmund
#F.~Robertson (see cite{HKRR84}, cite{Rob88}).
#
#The  higher-level  function  is  `simplify`.  The lower-level functions are
#`Eliminate`, `Search`, `SearchEqual`, and `FindCyclicJoins`.
#
#Some  of  these  functions  may  eliminate  generators,  but  they  do *not*
#introduce  new generators. However,  sometimes you will  need to substitute
#certain  words as  new generators  in order  to improve  your presentation.
#Therefore  there are the  functions `Substitute` and `SubstituteCyclicJoins`
#which introduce new generators.
#
#Finally  the functions `tracing`  and `images` can
#be used to determine and to display the images or preimages of the involved
#generators under the isomorphism which is defined by the sequence of Tietze
#transformations which are applied to a presentation.
#
#The  functions, `showpairs`, and `PrintOptions`,  can be useful. There are
#also  the  *Tietze  options*:  parameters  which  essentially influence the
#performance  of the  functions mentioned  above; they  are not specified as
#arguments of function calls. Instead, they are stored in the presentation.
rio(io::IO=stdout;p...)=IOContext(io,:limit=>true,p...)
function stringind(io::IO,n::Integer)
  if get(io,:TeX,false)
    n in 0:9 ? "_"*string(n) : "_{"*string(n)*"}"
  elseif get(io,:limit,false)
    if n<0 res=['₋']; n=-n else res=Char[] end
    for i in reverse(digits(n)) push!(res,Char(0x2080+i)) end
    String(res)
  else "_"*string(n)
  end
end

const supvec=collect("⁰¹²³⁴⁵⁶⁷⁸⁹⁻")

function stringexp(io::IO,n::Integer)
  if isone(n) ""
  elseif get(io,:TeX,false)
    n in 0:9 ? "^"*string(n) : "^{"*string(n)*"}"
  elseif get(io,:limit,false)
    if n==-1 return "⁻" end
    if n<0 res=['⁻']; n=-n else res=Char[] end
    for i in reverse(digits(n)) push!(res,supvec[i+1]) end
    String(res)
  else "^"*string(n)
  end
end

using PermGroups
using Combinat: tally
export AbsWord, @AbsWord, Presentation, FpGroup, Go, GoGo, conjugate,
       tryconjugate, simplify, relators, display_balanced, tracing, images,
       showgens, showpairs

plural(n,w)=string(n)*" "*w*(n==1 ? "" : "s")

#------------------ Abstract Words ----------------------------------
"""
An  `AbsWord` represents an  element of the  free group on some generators.
The  generators  are  indexed  by  `Symbols`.  For  example  the  `Absword`
representing `a³b⁻²a` is represented internally as
`[:a => 3, :b => -2, :a => 1]`. The mulitiplcation follows the group rule:
```julia-repl
julia> w=AbsWord([:a => 3, :b => -2, :a => 1])
a³b⁻²a

julia> w[4]
:b => -1

julia> w[3:4]
ab⁻

julia> length(w)
6

julia> w[begin:end]
a³b⁻²a

julia> w*AbsWord([:a=>-1,:b=>1])
a³b⁻

```
A positive `AbsWord` may be obtained by giving `Symbols` as arguments
```julia-repl
julia> AbsWord(:b,:a,:a,:b)
ba²b
```
"""
struct AbsWord
  d::Vector{Pair{Symbol,Int}}
  function AbsWord(v::Vector{Pair{Symbol,Int}};check=true)
    if check && length(v)>0
      ri=1
      for i in 2:length(v)
        if ri==0 || first(v[i])!=first(v[ri])
          ri+=1
          v[ri]=v[i]
        else c=last(v[ri])+last(v[i])
          if iszero(c) && ri>0 ri-=1
          else v[ri]=first(v[ri])=>c
          end
        end
      end
      resize!(v,ri)
    end
    new(v)
  end
end

AbsWord(x::Symbol...)=isempty(x) ? AbsWord("") : AbsWord([s=>1 for s in x])

"`@AbsWord x,y` define `AbsWord`s. The same as `x=AbsWord(:x);y=AbsWord(y)`"
macro AbsWord(t)
  if t isa Expr
    for v in t.args
      Base.eval(Main,:($v=AbsWord($(Core.QuoteNode(Symbol(v))))))
    end
  elseif t isa Symbol
    Base.eval(Main,:($t=AbsWord($(Core.QuoteNode(t)))))
  end
end

"""
`AbsWord(s::String)` computes an `AbsWord` from a line of
`display_balanced`  or from the display of an `AbsWord`. It is assumed that
all the symbols in the `AbsWord` print as one (non-superscript) character.
```julia-repl
julia> AbsWord("aec⁻")
aec⁻

julia> AbsWord("c²=c⁻²")
c⁴
```
"""
function AbsWord(s::AbstractString)
  ss=split(s,"=")
  if ss[end]=="1" pop!(ss) end
  if length(ss)==2 return AbsWord(ss[1])*inv(AbsWord(ss[2])) end
  s=ss[1]
  res=Pair{Symbol,Int}[]
  prevl=nothing
  exp=Int[]
  function num(exp::Vector{Int})
    sgn=1;n=0
    for i in exp
      if i==-1 sgn=i
      else n=10*n+i
      end
    end
    sgn*(n==0 ? 1 : n)
  end
  for c in s
    i=findfirst(==(c),supvec)
    if isnothing(i)
      if !isnothing(prevl) push!(res,Symbol(prevl)=>num(exp)) end
      prevl=c
      exp=Int[]
    else push!(exp,i==11 ? -1 : i-1)
    end
  end
  if !isnothing(prevl) push!(res,Symbol(prevl)=>num(exp)) end
  AbsWord(res)
end

"returns symbol for an AbsWord of length 1"
function mon(w::AbsWord)
  if length(w.d)!=1 || last(w.d[1])!=1 error("not a generator") end
  first(w.d[1])
end

"""
`AbsWord(word::Vector{Int}, gens::Vector{AbsWord})`

`AbsWord(word::Vector{Int}, gens::Vector{Symbol})`

`AbsWord(word::Vector{Int}, P::Presentation)`

`AbsWord(word::Vector{Int})`

Transforms Tietze word `word` to an absword, using given generators `gens`.
In the third case, use `gens(P)`.
In the fourth case, the list `Symbol.('a':'z')` is used.

```julia-repl
julia> AbsWord([-1,-2,1,-3,-3],AbsWord.([:a,:b,:c]))
a⁻b⁻ac⁻²

julia> AbsWord([-1,-2,1,-3,-3],[:a,:b,:c])
a⁻b⁻ac⁻²

julia> AbsWord([-1,-2,1,-3,-3])
a⁻b⁻ac⁻²
```
"""
AbsWord(tz::Vector{Int},gens::Vector{AbsWord})=AbsWord(tz,mon.(gens))
AbsWord(tz::Vector{Int},ss::Vector{Symbol})=AbsWord(map(i->ss[abs(i)]=>sign(i),tz))
AbsWord(tz::Vector{Int})=AbsWord(tz,Symbol.('a':'z'))

function Base.show(io::IO,a::AbsWord)
  if isone(a) print(io,".") end
  for (s,c) in a.d
    print(io,string(s))
    if c!=1 print(io,stringexp(io,c)) end
  end
end

Base.one(::Type{AbsWord})=AbsWord(Pair{Symbol,Int}[])
Base.isone(a::AbsWord)=isempty(a.d)
Base.:*(a::AbsWord,b::AbsWord)=AbsWord(vcat(a.d,b.d))
Base.inv(a::AbsWord)=AbsWord([k=>-v for (k,v) in reverse(a.d)];check=false)
Base.:^(a::AbsWord, n::Integer)=n>=0 ? Base.power_by_squaring(a,n) :
                                       Base.power_by_squaring(inv(a),-n)
Base.:^(a::AbsWord, b::AbsWord)=inv(b)*a*b
Base.:/(a::AbsWord, b::AbsWord)=a*inv(b)
Base.:\(a::AbsWord, b::AbsWord)=inv(a)*b
Base.length(a::AbsWord)=sum(x->abs(last(x)),a.d;init=0)
Base.copy(a::AbsWord)=AbsWord(copy(a.d))
Base.:(==)(a::AbsWord,b::AbsWord)=a.d==b.d
Base.hash(a::AbsWord,h::UInt)=hash(a.d,h)

function Base.getindex(w::AbsWord,i::Integer)
  if i<=0 error("index out of bounds") end
  for (s,n) in w.d
    if i>abs(n) i-=abs(n)
    elseif n<0 return s=>-1
    else return s=>1
    end
  end
  error("index out of bounds")
end

Base.getindex(w::AbsWord,i::AbstractVector{<:Integer})=AbsWord(getindex.(Ref(w),i))

Base.firstindex(w::AbsWord)=1
Base.lastindex(w::AbsWord)=length(w)
#-------------------------- FpGroups -----------------------------
@GapObj struct FpGroup <: Group{AbsWord}
  gens::Vector{AbsWord}
  rels::Vector{AbsWord}
end

FpGroup(gens::Vector{AbsWord})=FpGroup(gens,AbsWord[],Dict{Symbol,Any}())

FpGroup(s::Symbol...)=FpGroup(AbsWord.(collect(s)),AbsWord[],Dict{Symbol,Any}())

function Base.show(io::IO,G::FpGroup)
  print(io,"FreeGroup(",join(gens(G),","),")")
  if length(G.rels)>0 print(io,"/[");join(io,G.rels,",");print(io,"]") end
end

function Base.:/(G::FpGroup,rel::Vector{AbsWord})
  FpGroup(G.gens,vcat(G.rels,rel),Dict{Symbol,Any}())
end

# Tietze word from AbsWord
function Groups.word(G::FpGroup,w::AbsWord)
  res=Vector{Int}(undef,length(w))
  x=1
  for (s,m) in w.d
    p=findfirst(x->x.d[1][1]==s,gens(G))
    if isnothing(p) error(w," is not a word for ",G) end
    for i in 1:abs(m)
      res[x]=m<0 ? -p : p
      x+=1
    end
  end
  res
end

function Base.length(F::FpGroup)
  get!(F,:length)do
    if isempty(F.rels) error("infinite") end
    m=Base.get_extension(GroupPresentations,:Gap4GP)
    if isnothing(m)
      error("'length' needs GAP; do 'using GAP' to make extension available")
    end
    m.length(F)
  end
end

function Groups.elements(F::FpGroup)
  get!(g,:elements)do
    m=Base.get_extension(GroupPresentations,:Gap4GP)
    if isnothing(m)
      error("'elements' needs GAP; do 'using GAP' to make extension available")
    end
    m.elements(F)
  end
end

function Groups.elements(F::FpGroup,i)
  if isempty(F.rels)
    if i==0 return AbsWord() end
    l1=map(x->x.d[1],gens(F))
    l=map(x->x[1]=>-x[2],l1)
    l=vcat(l,l1)
    res=AbsWord.(collect.(Iterators.product(fill(l,i)...)))
    filter(x->length(x)==i,res)
  else
    filter(x->length(x)==i,elements(F))
  end
end

#--------------------------  Tietze words and structs -------------
"""
`TietzeWord(word::AbsWord, generators::Vector{AbsWord})`

Let  `generators` be a  list of abstract  generators and `word` an abstract
word   in  these   generators.  The   function  `TietzeWord`   returns  the
corresponding (reduced) Tietze word.

```julia-repl
julia> F=FpGroup(:a,:b,:c)
FreeGroup(a,b,c)

julia> GroupPresentations.TietzeWord(comm(F(1),F(2))*inv(F(3)^2*F(2)),gens(F))
5-element Vector{Int64}:
 -1
 -2
  1
 -3
 -3
```
"""
function TietzeWord(w::AbsWord,gens::Vector{AbsWord})
  ss=mon.(gens)
  res=Int[]
  for (s,c) in w.d
    p=findfirst(==(s),ss)
    if p===nothing error(s," is not in ",gens) end
    if c>0 append!(res,fill(p,c))
    else append!(res,fill(-p,-c))
    end
  end
  res
end

@GapObj mutable struct Presentation
  generators::Vector{AbsWord} # copy of initial gens
  inverses::Vector{Int}
  relators::Vector{Vector{Int}}
  flags::Vector{Int}
  modified::Bool
  numredunds::Int
end
# possible .prop fields: tietze, eliminationsLimit,
# expandLimit, generatorsLimit, lengthLimit, loopLimit, debug,
# saveLimit, searchSimultaneous, protected, status

AbsWord(tz::Vector{Int},P::Presentation)=AbsWord(tz,gens(P))

Base.copy(P::Presentation)=Presentation(copy(P.generators),
 copy(P.inverses),copy(copy.(P.relators)),copy(P.flags),P.modified,
 P.numredunds,copy(P.prop))

totallength(P::Presentation)=sum(length,P.relators;init=0)
Groups.gens(P::Presentation)=P.generators
Groups.ngens(P::Presentation)=length(P.generators)
(P::Presentation)(l::Int...)=isempty(l) ? AbsWord() : AbsWord(collect(l),P)

" `relators(P::Presentation)` return the relators of `P` as `AbsWord`s. "
relators(P::Presentation)=AbsWord.(P.relators,Ref(gens(P)))

"""
`FpGroup(P::Presentation)`

returns the finitely presented group defined  by the presentation `P`.
"""
function FpGroup(P::Presentation)
  debug(P,3,"# converting presentation to FpGroup")
  if P.numredunds>0 RemoveGenerators(P) end
  sort!(P)
  G=FpGroup(copy(gens(P)),copy(relators(P)),Dict{Symbol,Any}())
  if haskey(P, :imagesOldGens)
    G.imagesOldGens=deepcopy(P.imagesOldGens)
  end
  if haskey(P, :preImagesNewGens)
    G.preImagesNewGens=deepcopy(P.preImagesNewGens)
  end
  G
end

"""
`Presentation(G::FpGroup; debug=1)`

returns  the  presentation  corresponding  to  the given finitely presented
group `G`.

The  `debug` keyword controls the amount of informative printing
done  by Tietze transformations when operating on the created presentation.
The  default value  1 implies  messages to  be displayed  by most  of these
functions.  A `debug`  value of  0 will  suppress these messages, whereas a
`debug` value of 2 will cause some additional output.
"""
Presentation(G::FpGroup;debug::Integer=1)=Presentation(G.gens,G.rels;debug)

"""
`Presentation(gens::Vector{AbsWord},rels::Vector{AbsWord};debug::Int=1)`

constructs  a `Presentation`  from a  list `gens`  of generators (which can
also be given as a list of `Symbol`s) and a list `rels` of relators.
"""
function Presentation(gens::Vector{AbsWord},grels::Vector{AbsWord};debug::Int=1)
  numgens=length(gens)
  rels=map(w->reduceword(TietzeWord(w,gens);cyclically=true),grels)
  P=Presentation(deepcopy(gens),numgens:-1:-numgens,rels,fill(1,length(rels)),
                 false, 0, Dict{Symbol,Any}())
  P.eliminationsLimit=100
  P.expandLimit=150
  P.generatorsLimit=0
  P.lengthLimit=0 # infinity
  P.loopLimit=0 # infinity
  P.debug=debug
  P.saveLimit=10
  P.searchSimultaneous=20
  PrintStatus(P,2)
  P.protected=ngens(P);HandleLength1Or2Relators(P);P.protected=0
  sort!(P)
  PrintStatus(P,2)
  P
end

Presentation(s::Vector{Symbol},r::Vector{AbsWord};debug=1)=
  Presentation(AbsWord.(s),r;debug)

"""
`Presentation(s::String; debug=1)`

computes  a  presentation  from  the  output  of  `display_balanced`  or of
`relators`.  It is assumed  that writing any  generator of the presentation
takes one (non-superscript) character.
"""
function Presentation(s::String;debug=1)
  s=replace(s,r"\n\s*="=>"=")
  l=split(s,"\n")
  l=map(s->replace(s,r"^ *[0-9]*: *"=>""),l)
  l=filter(x->match(r"^ *$",x)===nothing,l)
  l=map(x->replace(x,r"\s*$"=>""),l)
  rels=AbsWord.(l)
  l=map(x->first.(x.d),rels)
  if isempty(l) return Presentation(AbsWord[],AbsWord[];debug) end
  atoms=length(l)==1 ? unique(l[1]) : union(map(x->first.(x.d),rels)...)
  sort!(atoms)
  Presentation(AbsWord.(atoms),rels;debug)
end

"""
`Presentation(v::Vector{Vector{Int}};debug=1)`

Constructs a presentation from `v`, a list of Tietze words, using for generators
the alphabet `a:z`, then `x_i` for integers `i`.
"""
function Presentation(v::Vector{Vector{Int}};debug=1)
  gens=sort(unique(abs(x) for r in v for x in r))
  if isempty(gens) return Presentation(AbsWord[],AbsWord[];debug) end
  g=setdiff(1:maximum(gens),gens)
  if !isempty(g) error("generators $g missing") end
  vars=length(gens)<=26 ? Symbol.('a':'a'+length(gens)-1) :
      map(i->Symbol("x",stringind(rio(),i)),1:length(gens))
  Presentation(AbsWord.(vars),map(r->AbsWord(r,vars),v);debug)
end

Base.getindex(T::Presentation,i)=T.inverses[ngens(T)+1-i]
Base.setindex!(T::Presentation,j,i)=T.inverses[ngens(T)+1-i]=j

# reduce cyclically a Tietze word
function reduceword(w::Vector{Int};cyclically=false)
  i=1;
  res=Int[]
  for j in eachindex(w)
    if isempty(res) || res[end]!=-w[j] if w[j]!=0 push!(res,w[j]) end
    else pop!(res)
    end
  end
  b=1;e=length(res)
  if cyclically while b<e && res[b]==-res[e] b+=1;e-=1 end end
  res[b:e]
end

# reduce cyclically a Tietze word using the table of inverses
function reduceword(w::Vector{Int},T::Presentation)
  i=1;
  res=Int[]
  for j in eachindex(w)
    if T[w[j]]==0 continue end
    if isempty(res) || T[res[end]]!=T[-w[j]] push!(res,T[w[j]])
    else pop!(res)
    end
  end
  b=1;e=length(res)
  while b<e && T[res[b]]==T[-res[e]] b+=1;e-=1 end
  res[b:e]
end

# sorts by length relators and remove duplicates and empty ones
function Base.sort!(T::Presentation)
  if isempty(T.relators) return end
  p=sortperm(eachindex(T.relators),by=i->(length(T.relators[i]),T.relators[i],-T.flags[i]))
  p=unique(i->T.relators[i],p)
  while !isempty(p) && isempty(T.relators[p[1]]) p=p[2:end] end
  if p==eachindex(T.relators) return end
  if T.debug>=3 println("#sort! the relators") end
  numrels=length(p)
  T.relators[1:numrels].=T.relators[p];resize!(T.relators,numrels)
  T.flags[1:numrels].=T.flags[p];resize!(T.flags,numrels)
end

function debug(P,i,arg...)
  if P.debug>=i println("#",arg...) end
end

include("toddcoxeter.jl")
"""
`PermGroups.PermGroup(F::FpGroup;maxcosets=1000)`
`PermGroups.PermGroup(P::Presentation;maxcosets=1000)`

Computes  a `PermGroup` for `F` or `P`. The algorithm may need to enumerate
a number of cosets equal to several times the size of the result group.
"""
PermGroups.PermGroup(p::Presentation;opt...)=toddcoxeter(p.relators;opt...)
PermGroups.PermGroup(F::FpGroup;opt...)=PermGroup(Presentation(F);opt...)

"""
`PrintStatus(P::Presentation)`

updates  and prints  the current  status of  a presentation  `P`, i.e., the
number  of generators, the number of relators,  and the total length of all
relators.

The  printing is suppressed if  none of the three  values has changed since
the last call.
"""
function PrintStatus(P::Presentation,debug=0)
  if P.debug<debug return end
  status=[ngens(P)-P.numredunds, length(P.relators),totallength(P)]
  if !haskey(P,:status) || status!=P.status
    P.status=status
    show(stdout,P)
    println()
  end
end

function Base.show(io::IO, p::Presentation)
  print(io,"Presentation")
  if haskey(p,:name) print(io," ",p.name) end
  print(io,": ",plural(ngens(p)-p.numredunds,"generator"),
   ", ",plural(length(p.relators),"relator")," of total length ",totallength(p))
end

"""
A `Presentation` internally has:

  - `generators`: the generators, as `AbsWord`s of length 1;

  - `relators`: the relators as Tietze words, i.e. as lists of positive or
    negative  generator numbers; a  generator known to  be an involution is
    always stored as a positive number.

  - `inverses`:  has length `2n+1`, where `n=ngens(P)`, and is read through
    `P[i]=inverses[n+1-i]`;  thus  `P[0]`  is  the  middle  entry, which is
    always `0`. For a generator `i` there are 3 cases:

  - `P[i]==i` and `P[-i]==-i`: normal case
  - `P[i]==P[-i]==i`: `i` is known to be an involution;
  - `P[i]==j` and `P[-i]==-j` with `j!=i`: `i` has been identified with `j`;
    this happens only between a transformation and `RemoveGenerators`. If
    `j==0` this means `i` has been eliminated.

    So `inverses` is a projection onto the surviving generators: we have
    `P[P[i]]==P[i]`, and `P[i]`, `P[-i]` are equal or opposite.

The *use* of `inverses` is in `reduceword(w,T)`, which replaces each `i` by
`P[i]`  and cancels `i` against `P[-i]`; the entry `P[-i]==P[i]==i` is what
makes  it use `i²=1`. The invariant is `P[-i]==P[i]==i` if and only if some
relator is `[i,i]` and is flagged `3`. Everything which rewrites, merges or
drops  a relator has to maintain it, which  is why `sort!` keeps the flag 3
of duplicate relators.

  - `flags`: is parallel to `relators` (has the same length)
    and records the state of each relator while simplifying:

- `0`: not changed since an earlier round of `Search`.
- `1`: new, or changed since it was last looked at: eligible for `Search`.
- `2`: changed during the current round of `Search`, so not to be used
    again in that round;
- `3`: this relator is `[i,i]` and is the witness of `P[-i]==P[i]==i`.

The  values `0`, `1`, `2` are only  a scheduling heuristic in `Search`. The
value  `3`  matters:  a  relator  carrying  it  must  not  be reduced away,
rewritten  or dropped.  `HandleLength1Or2Relators` sets  it when it detects
`i²=1`,  `sort!` keeps it when it  merges duplicate relators, and `SearchC`
leaves such relators alone.

  - `modified`: whether the last transformation changed anything;

  - `numredunds`: how many generators have been eliminated but not yet
     renumbered out (see `RemoveGenerators`).
"""
function Base.dump(T::Presentation)
  l=filter(i->T[-i]!=-T[i],eachindex(gens(T)))
  l1=filter(i->T[-i]==T[i],l)
  println("# F relator")
  for i in eachindex(T.relators)
    println(rio(),"$i:",T.flags[i]," ",T(T.relators[i]...))
  end
  print("gens=");join(rio(),gens(T),",")
  if length(l1)>0 print(" involutions:");join(rio(),gens(T)[l1],",") end
  println(" modified=",T.modified," numredunds=",T.numredunds)
  if l1!=l display(stack(map(i->[gens(T)[i],T[i],T[-i]],l))) end
  println("inverses:",T.inverses)
end

"""
`display_balanced(p)`

Tries  to display the  presentation `p` as  equalities between words of the
same length.

```julia-rep1
julia> relators(p)
7-element Vector{AbsWord}:
 dbd⁻b⁻
 bab⁻a⁻
 caca⁻c⁻a⁻
 adad⁻a⁻d⁻
 bcbc⁻b⁻c⁻
 cdcdc⁻d⁻c⁻d⁻
 a⁻dc⁻acdc⁻a⁻cd⁻ac⁻

julia> display_balanced(p)
1: db=bd
2: ba=ab
3: cac=aca
4: ada=dad
5: bcb=cbc
6: cdcd=dcdc
7: a⁻dc⁻acd=ca⁻dc⁻ac
```
"""
function display_balanced(T;dumb=false)
  f(i,w)=vcat(w,w)[i:i+div(length(w),2)-1]
  used=Set{Int}();for x in T.relators, y in x push!(used,abs(y)) end
  if ngens(T)>length(used)
    println(plural(ngens(T)-length(used)," free generator"))
  end
  for (i,w) in enumerate(T.relators)
    lw=length(w)
    if isodd(lw) || dumb println(rio(),i, ": ", T(w...), "=1")
    elseif lw>0
      m=argmax(map(i->count(>(0), f(i,w)), 1:lw))
      println(rio(),"$i: ",T(f(m,w)...),"=",T(-reverse(f(m+div(lw,2),w))...))
    end
  end
end

"""
`showgens(P,list=eachindex(gens(P)))`

prints  the generators of `P` with the total number of their occurrences in
the  relators, and notes involutions. A  second `list` argument prints only
those generators.
"""
function showgens(T::Presentation,list=eachindex(gens(T)))
  if ngens(T)==0 println("#I  there are no generators");return end
  occur=Occurrences(T)
  if list isa Integer list=[list] end
  for i in list
    print(i, ". ", T(i)," ",plural(occur[1][i],"occurrence"))
    print(" (",plural(occur[3][i],"time")," in relator ",occur[2][i],")")
    if T[-i]>0 print(" involution") end
    println()
  end
end

"""
`showpairs(P;n=10)`

shows  the  `n`  most  frequently  occurring squarefree relator subwords of
length 2 with their number of occurrences.`n=0` is interpreted as infinity.

This is useful information in the context of the `Substitute` command.
"""
function showpairs(P::Presentation;n::Integer=10)
  for (m,(num,i,j,k)) in enumerate(MostFrequentPairs(P,n))
    if k>1 i=-i end
    if isodd(k) j=-j end
    println(rio(),"#I $m. ",plural(num,"occurrence")," of ",P(i,j))
  end
end

"""
`add_generator(P[, generator])`

`add_generator` adds a new generator to the list of generators.

If  you don't specify a second  argument, then `add_generator` will define a
new  abstract generator which is the first unused letter in `a:z` if it exists,
otherwise the first unused `x_i` for `i` starting with 1.

If  you specify  a second  argument, then  `generator` must  be an abstract
generator which does not yet occur in the presentation. `add_generator` will
add it to the presentation.
"""
function add_generator(P::Presentation,gen=nothing)
  Check(P)
  if gen===nothing
    gen=newgen(P)
    debug(P,1,"add_generator  new generator is ", gen)
  else
    if !(gen isa AbsWord && length(gen)==1)
      error("second argument must be an abstract generator")
    end
    if gen in P.generators || inv(gen) in P.generators
      println("#I generator ",gen," is already in the presentation")
      return
    end
    push!(P.generators,gen)
    P.inverses=vcat([length(P.generators)],P.inverses, [-length(P.generators)])
  end
  tracing(P,false)
  P.modified=true
end

"""
newgen(P::Presentation)

let  'x' be the first unused letter in `a:z` for `gens(P)` if it exists (if
tracing  the first letter never used), otherwise the first unused `x_i` for
`i`   starting  with  1.   A  new  abstract   generator  `x`  is  added  to
`P.generators`.
"""
function newgen(P::Presentation)
  if haskey(P, :imagesOldGens) || haskey(P, :preImagesNewGens)
    i=maximum(findall(j->AbsWord(Symbol('a'+j)) in gens(P),0:25))
    if i==26 i=nothing end
  else
    i=findfirst(j->!(AbsWord(Symbol('a'+j)) in gens(P)),0:25)-1
  end
  if isnothing(i)
    i=1
    while AbsWord(Symbol("x_",i)) in gens(P) i+=1 end
    gen=AbsWord(Symbol("x_",i))
  else
    gen=AbsWord(Symbol('a'+i))
  end
  push!(P.generators,gen)
  P.inverses=vcat(length(P.generators),P.inverses,-length(P.generators))
  gen
end

"""
`add_relator(P::Presentation, word::AbsWord)`

adds  the word `word` to the list of  relators.  `word` must
be a word in the generators of the given presentation.
"""
function add_relator(P::Presentation, word::AbsWord;stop=true)
  Check(P)
  debug(P,3,"add_relator adding ", word)
  rel=reduceword(TietzeWord(word,gens(P));cyclically=true)
  if !isempty(rel)
    push!(P.relators,rel)
    push!(P.flags,1)
    P.modified=true
  end
  if stop tracing(P,false) end
end

"""
HandleLength1Or2Relators(presentation) . . . handle short relators

`HandleLength1Or2Relators`  searches  for  relators  of  length  1 or 2 and
performs suitable Tietze transformations for each of them:

Generators occurring in relators of length 1 are eliminated.

Generators  occurring  in  square  relators  of  length  2 are marked to be
involutions.

If  a  relator  of  length  2  involves  two different generators, then the
generator  with the larger  number is substituted  by the other  one in all
relators and finally eliminated from the set of generators.
"""
function HandleLength1Or2Relators(T::Presentation)
# dump(T)
  debug(T,3,"Handle length 1 or 2 relators")
  gens=T.generators
  rels=T.relators
  done=false
  while !done
    done=true
    i=0
    while i<length(T.relators)
      i+=1
      lg=length(rels[i])
      if !(lg in (1,2)) || T.flags[i]>2 continue end
      rep1=T[rels[i][1]]
      if lg==1
        rep1=abs(rep1)
        if rep1>T.protected
          if T[rep1]==rep1 T.numredunds+=1 end
          T[rep1]=T[-rep1]=0
          debug(T,2,"Handle12 eliminating ",gens[rep1]," redund=",T.numredunds)
          UpdateGeneratorImages(T, rep1, Int[])
          done=false
        end
      else
        rep2=T[rels[i][2]]
        if abs(rep2)<abs(rep1) rep1,rep2=rep2,rep1 end
        if rep1<0 rep1,rep2=(-rep1,-rep2) end
        if rep1==0 # already eliminated
          rep2=abs(rep2)
          if rep2>T.protected
            if T[rep2]==rep2 T.numredunds+=1 end
            T[rep2]=T[-rep2]=0
            debug(T,2,"Handle12 eliminating ",gens[rep2]," redund=",T.numredunds)
            UpdateGeneratorImages(T, rep2, Int[])
            done=false
          end
        elseif rep1!=-rep2  # otherwise not reduced
          if rep1!=rep2 # not an involution
            if T[rep2]==T[-rep2] && T[-rep1]<0 # rep2^2=1 => rep1^2=1
              push!(rels,[rep1, rep1])
              push!(T.flags,3)
              T[-rep1]=rep1
            end
            if abs(rep2)>T.protected
              if T[rep2]==rep2 T.numredunds+=1 end
              T[rep2]=T[-rep1];T[-rep2]=rep1
              if haskey(T,:imagesOldGens) || T.debug>=2
                if rep2>0 rep1=T[-rep1] end
                debug(T,2,"Handle12 eliminating ",gens[abs(rep2)],"=",
                    AbsWord([rep1],T)," redund=",T.numredunds)
                UpdateGeneratorImages(T,abs(rep2),[rep1])
              end
              done=false
            end
          elseif T[-rep1]<0  # an involution not yet detected
            rels[i]=[rep1,rep1]
            T.flags[i]=3
            T[-rep1]=rep1
            done=false
          end
        end
      end
    end
    if !done
      for i in eachindex(T.generators)
        if T[i]!=i T[i],T[-i]=T[T[i]],T[-T[i]] end
      end
# the next loop is FunTzReplaceGens
      for i in eachindex(T.relators)
        if length(T.relators[i])==2 && T.flags[i]==3 && abs(T[T.relators[i][1]])== abs(T.relators[i][1])continue end#don't remove involutions
        T.relators[i]=reduceword(T.relators[i],T)
        T.flags[i]=1
      end
    end
  end
  if T.numredunds>0 RemoveGenerators(T) end
# dump(T)
end

"""
`remove_relator(P,n)`

removes the `n`th relator of `P`.
"""
function remove_relator(T::Presentation, n)
  if !(n in eachindex(T.relators)) error("relator number out of range") end
  rels=T.relators
  debug(T,3,"remove_relator removing $n-th relator")
  leng=length(T.relators[n])
  if leng==2 && rels[n][1]==rels[n][2] # forget_involution
    num=abs(rels[n][1])
    if T[-num]==num T[-num]=-num end
  end
  rels[n]=Int[]
  sort!(T)
  PrintStatus(T,2)
  tracing(T,false)
end

"""
`simplify(G)`

`simplify`  applies  Tietze  transformations  to a  copy of  the
presentation of the given finitely presented group `G` in order to reduce
it with respect to  the number of generators, the number of relators, and
the relator lengths.

`simplify` returns the resulting finitely presented group (which
is isomorphic to `G`).

```julia-repl
julia> @AbsWord a,b,c,d,e,f

julia> F=FpGroup([a,b,c,d,e,f])
FreeGroup(a,b,c,d,e,f)

julia> G=F/[a^2,b^2,d*f^-1,e^2,f^2,a*b^-1*c,a*e*c^-1,b*d^-1*c,c*d*e^-1,a*f*c^-2,c^4]
FreeGroup(a,b,c,d,e,f)/[a²,b²,df⁻,e²,f²,ab⁻c,aec⁻,bd⁻c,cde⁻,afc⁻²,c⁴]

julia> simplify(G)
FreeGroup(a,c)/[a²,ac⁻ac⁻,c⁴]
```

In fact, the call

```julia-rep1
julia> simplify(G)
```

is an abbreviation of the call sequence

```julia-rep1
julia> P=Presentation(G;debug=0);simplify(P);FpGroup(P)
```
which applies  a rather simple-minded strategy of  Tietze transformations
to the intermediate presentation `P`.
If for  some  concrete group the resulting presentation  is unsatisfying,
then  you  should  try  a  more  sophisticated,  interactive  use of  the
available Tietze transformation functions  (see "Tietze Transformations").
"""
function simplify(G::FpGroup)
  T=Presentation(G;debug=0)
  GoGo(T)
  FpGroup(T)
end

"""
`Check(presentation)`

checks some fields of the  given presentation to be consistent. 
"""
function Check(P::Presentation)
  n=length(P.generators)
  length(P.inverses)==2n+1 || error("inconsistent inverses")
  length(P.flags)==length(P.relators) || error("inconsistent relators")
  P.inverses[n+1]==0 || error("P[0] should be 0")
  for i in 1:n
    a=P[i];b=P[-i]
    (0<=a<=n && -n<=b<=n) || error("inconsistent image for generator $i")
    (a==0)==(b==0) || error("generator $i is eliminated on one side only")
    (b==a || b==-a) || error("inconsistent images $a and $b for generator $i")
    a>0 && P[a]!=a && error("P[$i]=$a is not a fixed point of the table")
  end
  for r in P.relators, x in r
    P[x]!=0 && continue
    error("relator $r mentions eliminated generator $x")
  end
  for k in eachindex(P.flags)
    P.flags[k] in 0:3 || error("flag $(P.flags[k]) on relator $k")
    if P.flags[k]==3
      r=P.relators[k]
      (length(r)==2 && r[1]==r[2] && abs(P[r[1]])==abs(r[1])) ||
        error("relator $k = $r has flag 3 but is not an involution witness")
    end
  end
  for i in 1:n
    P[-i]==P[i]==i || continue
    any(k->begin r=P.relators[k]
            P.flags[k]==3 && length(r)==2 && r[1]==r[2] && abs(P[r[1]])==i end,
        eachindex(P.relators)) ||
      error("generator $i is recorded as an involution, but no relator [i,i] with flag 3 witnesses it")
  end
end

"""
`Eliminate(P,gen::AbsWord)`
`Eliminate(P,n::Int=1)`

`Eliminate`  tries to eliminate one or  more generators from a presentation
`P` via Tietze transformations.

Any  relator  which  contains  some  generator  just  once  can  be used to
substitute  that generator by  a word in  the remaining generators. If such
generators  and relators  exist, then  `Eliminate` chooses  a generator for
which  the  product  of  its  number  of  occurrences and the length of the
substituting  word is minimal,  and then it  eliminates this generator from
the  presentation, provided that the resulting total length of the relators
does  not exceed  the Tietze  option parameter  `P.spaceLimit`. The default
value  of `P.spaceLimit` is `infinity`, but  you may alter it appropriately
(see Tietze options below).

If you specify a generator `gen` as  second argument, then  `Eliminate`
only tries to eliminate that generator.

If you specify an integer `n` as second argument, then `Eliminate` tries to
eliminate   up   to   `n`   generators.   The   calls   `Eliminate(P)`  and
`Eliminate(P,1)` are equivalent.
"""
function Eliminate(P::Presentation,n=1)
  debug(P,3,"Eliminate n=$n")
  Check(P)
  gennum=0
  if !(n isa Integer)
    gennum=findfirst(==(n),gens(P))
    if gennum===nothing
      error("usage: Eliminate(presentation[,nbgens or gen])")
    end
    n=1
  end
  if n==1
    if gennum==0 EliminateGen1(P)
    else EliminateGen(P, gennum)
    end
    PrintStatus(P)
    if P.numredunds>0 RemoveGenerators(P) end
    HandleLength1Or2Relators(P)
    sort!(P)
    PrintStatus(P)
  else
    saveelim=P.eliminationsLimit;P.eliminationsLimit=n
    savegenlim=P.generatorsLimit;P.generatorsLimit=length(P.generators)-n
    EliminateGens(P)
    P.eliminationsLimit=saveelim
    P.generatorsLimit=savegenlim
  end
end

# substitute gen by w
function SubstituteGen(T::Presentation, gen, w)
  if gen<0 w=-reverse(w); gen=-gen end
  if T[gen]==gen && T[-gen]==gen T[-gen]=-gen end # forget involution
  iw=-reverse(w)
  for i in eachindex(T.relators)
    res=Int[]
    occ=0
    for j in T.relators[i]
      if j==gen append!(res,w);occ+=1
      elseif j==-gen append!(res,iw);occ+=1
      else push!(res,j)
      end
    end
    if occ>0 T.flags[i]=1 end
    if length(res)==2 && allequal(res) T.relators[i]=res
    else T.relators[i]=reduceword(res,T)
    end
  end
  UpdateGeneratorImages(T, abs(gen), gen<0 ? iw : w)
end

"""
EliminateGen( presentation,n) . . . eliminates the nth generator

`EliminateGen` eliminates the Tietze generator tietze.generators[n]
if possible, i. e. if that generator can be isolated  in some appropriate
Tietze relator.  However,  the elimination  will not be  performed if the
resulting total length of the relators cannot be guaranteed to not exceed
the parameter P.lengthLimit (if !=0).
"""
function EliminateGen(P::Presentation, num::Int)
 spacelimit=totallength(P)+P.lengthLimit
  P.modified=false
  if !(num in eachindex(P.generators))
    error("EliminateGen: second argument is not a valid generator number")
  end
  occTotal,occRelNum,occMult=Occurrences(P, num)
  if occMult==1
    lg=length(P.relators[occRelNum])
    space=(occTotal-1)*(lg-1)-lg
    if P.lengthLimit!=0 && totallength(P)+space>spacelimit return end
    gen=num
    rel=P.relators[occRelNum]
    pos=findfirst(==(gen),rel)
    if pos===nothing
      gen=-gen
      pos=findfirst(==(gen),rel)
    end
    word=vcat(rel[pos+1:end],rel[1:pos-1])
    debug(P,2,"EliminateGen eliminating ", P.generators[num], "=",
         let a=gen>0 ? -reverse(word) : word;  AbsWord(word,P) end)
    SubstituteGen(P,-gen, word)
    P[num]=P[-num]=0
    P.numredunds+=1
    P.modified=true
  end
end

# returns 3 vectors:
# 1: nb. occurences of each generator
# 2: n0 of a relator of minimal length among those where occurs with minimal>0 multiplicity
# 3: the corresponding multiplicity
function Occurrences(T::Presentation)
  occ=fill(0,length(T.generators))
  min=fill(0,length(T.generators))
  mult=fill(0,length(T.generators))
  for j in eachindex(T.relators)
    lococc=fill(0,length(T.generators))
    for i in T.relators[j] lococc[abs(i)]+=1 end
    occ+=lococc
    for i in eachindex(T.generators)
      if lococc[i]==0 continue end
      if mult[i]==0 || lococc[i]<mult[i] ||
        (lococc[i]==mult[i] && length(T.relators[j])<length(T.relators[min[i]]))
          min[i]=j
          mult[i]=lococc[i]
      end
    end
  end
  [occ,min,mult]
end

Occurrences(T::Presentation,gen)=getindex.(Occurrences(T),gen)

"""
EliminateGen1(presentation)  . . . . . . .  eliminates a generator

`EliminateGen1`  tries to  eliminate a  Tietze generator:  If there are
Tietze generators which occur just once in certain Tietze relators,  then
one of them is chosen  for which the product of the length of its minimal
defining word  and the  number of its  occurrences  is minimal.  However,
the elimination  will not be performed  if the resulting  total length of
the  relators   cannot  be  guaranteed   to  not  exceed   the  parameter
P.lengthLimit (if !=0).
"""
function EliminateGen1(T::Presentation)
  occTotals,occRelNo,occMultiplicities=Occurrences(T)
  modified=false
  num=0
  space=0
  for i in T.protected+1:length(T.generators)
    if occMultiplicities[i]==1
      len=length(T.relators[occRelNo[i]])
      ispace=(occTotals[i]-1)*(len-1)-len
      if num==0 || ispace<=space
        num=i
        space=ispace
      end
    end
  end
  if num>0 && (T.lengthLimit==0 || space<=T.lengthLimit)
    gen=num
    occRelNo=occRelNo[num]
    rel=T.relators[occRelNo]
    len=length(T.relators[occRelNo])
    pos=findfirst(==(gen),rel)
    if pos===nothing
      gen=-gen
      pos=findfirst(==(gen),rel)
    end
#   display_balanced(T)
#   @show rel
    word=vcat(rel[pos+1:end],rel[1:pos-1])
    debug(T,2,"EliminateGen1 eliminating ",T.generators[num],"=",
      let a=gen>0 ? -reverse(word) : word;  AbsWord(word,T) end)
    SubstituteGen(T, -gen, word)
    T[num]=T[-num]=0
    T.numredunds+=1
    modified=true
  end
  T.modified=modified
end

"""
EliminateGens(P::presentation)

`EliminateGens`  repeatedly eliminates generators  from the presentation of
the given group until one of the following conditions is violated:

(1) The  current  number of  generators  is greater than P.generatorsLimit.
(2) The number of generators eliminated so far is less than P.eliminationsLimit.
(3) The  total length of the relators  has not yet grown  to a percentage
    greater than the parameter P.expandLimit.
(4) The  next  elimination  will  not  extend the total length to a value
    greater than the parameter P.lengthLimit.

the function will not eliminate any protected generators.
"""
function EliminateGens(P::Presentation)
  debug(P,3,"EliminateGens eliminating generators")
  redundantsLimit=5
  bound=totallength(P)*P.expandLimit//100
  modified=false
  P.modified=true
  num=0
  while P.modified && num<P.eliminationsLimit && totallength(P)<=bound &&
   length(P.generators)-P.numredunds>P.generatorsLimit
    EliminateGen1(P)
#   print("a ");PrintStatus(P);@show P.flags
    if P.numredunds==redundantsLimit
      RemoveGenerators(P)
#     print("b ");PrintStatus(P);@show P.flags
    end
    modified=modified || P.modified
    num+=1
  end
  P.modified=modified
  if P.numredunds>0
    RemoveGenerators(P)
#   print("c ");PrintStatus(P);@show P.flags
  end
  if modified
    HandleLength1Or2Relators(P)
#   print("d ");PrintStatus(P);@show P.flags
    sort!(P)
    PrintStatus(P,2)
  end
end

"""
`FindCyclicJoins(P)`

`FindCyclicJoins` performs Tietze transformations on a presentation `P`. It
searches  for pairs of  generators which generate  the same cyclic subgroup
and eliminates one of the two generators of each such pair it finds.

More  precisely: `FindCyclicJoins` searches for pairs of generators `a` and
`b`  such that (possibly after inverting  or conjugating some relators) the
set  of  relators  contains  the  commutator  `[a,b]`, a power `a^n`, and a
product  of the form `a^s  b^t` with `s` prime  to `n`. For each such pair,
`FindCyclicJoins` uses the Euclidian algorithm to express `a` as a power of
`b`, and then it eliminates `a`.
"""
function FindCyclicJoins(P::Presentation)
  debug(P,3,"FindCyclicJoins searching for cyclic joins")
  T=P
  T.modified=false
  newstart=true
  while newstart
    exponents=find_exponents(P)
    if iszero(exponents) return end
    newstart=false
    gens=T.generators
    numgens=length(T.generators)
    rels=T.relators
    numrels=length(rels)
    flags=T.flags
    i=0
    while i<numrels
      i+=1
      rel=rels[i]
      if length(T.relators[i])!=4 || rel[1]!=T[-rel[3]] || rel[2]!=T[-rel[4]] continue end
      num=abs.(rel[1:2])
      exp=exponents[num]
      if iszero(exp) continue end
      fac=[0, 0]
      e=[0, 0]
      j=0
      while j<numrels
        j+=1
        if length(T.relators[j])==0 || j==i continue end
        if !all(x->abs(x) in num,rels[j]) continue end
        rel=rels[j]
        e.=0
        for x in rel, i in 1:2
          if x==num[i] e[i]+=1 elseif x==-num[i] e[i]-=1 end
        end
        powers=count(!iszero,e)
        if powers<2 continue end
        for k in 1:2
          fac[k]=num[k]
          if exp[k]>0
            e[k]=mod(e[k], exp[k])
            if 2e[k]>exp[k]
              e[k]=exp[k]-e[k]
              fac[k]=-fac[k]
            end
          elseif e[k]<0
            e[k]=-e[k]
            fac[k]=-fac[k]
          end
          if fac[k]<0 fac[k]=T[fac[k]] end
        end
        for k in 1:2
          if e[k]>0 && e[3-k]==0
            exp[k]=gcd(e[k], exp[k])
            if exp[k] != exponents[num[k]]
              exponents[num[k]]=exp[k]
              e[k]=exp[k]
            end
          end
        end
        rel=vcat(fill(fac[1],e[1]), fill(fac[2],e[2]))
        if rels[j]!=rel
          rels[j]=rel
          flags[j]=1
          T.modified=true
          debug(P,3,"FindCyclicJoins rels[$j] reduced to ",rels[j])
        end
        if e[1]==1 n=num[1]
        elseif e[2]==1 n=num[2]
        else
          n=0
          for k in 1:2
            if n==0 && e[k]>1 && gcd(e[k], exp[k])==1
              gen=gens[num]
              if fac[1]<0 gen[1]=inv(gen[1]) end
              if fac[2]<0 gen[2]=inv(gen[2]) end
              word=gen[k]*gen[3-k]^(e[3-k]*gcdx(e[k], exp[k])[2])
              add_relator(P, word)
              numrels=length(T.relators)
              n=num[k]
            end
          end
        end
        if n!=0 && P.generatorsLimit<numgens
          Eliminate(P)
          T.modified=true
          j=numrels
          i=numrels
          if 1<numgens newstart=true end
        end
      end
    end
  end
  if T.modified
    HandleLength1Or2Relators(P)
    sort!(P)
    PrintStatus(P,1)
  end
end

"""
`Go(P::Presentation[,silent])`

`Go`  performs Tietze transformations on a  presentation `P`. It is perhaps
the  most convenient of the  interactive Tietze transformation commands. It
offers  a  kind  of  default  strategy  which,  in  general, saves you from
explicitly calling the lower-level commands it involves.

Roughly  speaking, `Go` consists of a  loop over a procedure which involves
two  phases:  In  the  *search  phase*  it calls `Search` and `SearchEqual`
described  below which  try to  reduce the  relator lengths by substituting
common  subwords  of  relators,  in  the  *elimination  phase* it calls the
command  `Eliminate` described below  (or, more precisely,  a subroutine of
`Eliminate`  in order to save some  administrative overhead) which tries to
eliminate  generators  that  can  be  expressed  as  words in the remaining
generators.

If  `Go`  succeeds  in  reducing  the  number  of generators, the number of
relators,  or the total  length of all  relators, then it  displays the new
status  before returning (provided that you did  not set debug to
zero).  However, it does not  provide any output if  all these three values
have  remained unchanged,  even if  the `SearchEqual`  command involved has
changed  the  presentation  such  that  another  call of `Go` might provide
further  progress. Hence, in such a case  it makes sense to repeat the call
of  the command for  several times (or  to call instead  the `GoGo` command
which we will describe next).

As an example  we  compute  a presentation of a  subgroup of index 408 in
`PSL(2,17)`.

```
|    gap> F2 := FreeGroup( "a", "b" );;
    gap> G := F2 / [ F2.1^9, F2.2^2, (F2.1*F2.2)^4, (F2.1^2*F2.2)^3 ];;
    gap> a := G.1;;  b := G.2;;
    gap> H := Subgroup( G, [ (a*b)^2, (a^-1*b)^2 ] );;
    gap> Index( G, H );
    408
    gap> P := PresentationSubgroup( G, H );
    << presentation with 8 gens and 36 rels of total length 111 >>|
    julia> P=Presentation("
    1: a=a⁻
    2: g=g⁻
    3: c=c⁻
    4: f=f⁻
    5: h=h⁻
    6: d=d⁻
    7: fbc=1
    8: db⁻a=1
    9: gdb=1
    10: ceb=1
    11: eab=1
    12: cab⁻=1
    13: b⁻³=1
    14: cb⁻d=1
    15: acb=1
    16: afb⁻=1
    17: b³=1
    18: bac=1
    19: cba=1
    20: ab⁻c=1
    21: hab⁻=1
    22: dab=1
    23: b⁻ec=1
    24: b⁻ad=1
    25: b⁻ca=1
    26: agb=1
    27: ab=b⁻a⁻
    28: ab=bc⁻
    29: ac=c⁻e
    30: ac=b²
    31: cb⁻=bc⁻
    32: aca=c⁻a⁻c⁻")
    Presentation: 8 generators, 32 relators of total length 98
    julia> P.protected=2;
    julia> P.debug=2;
    julia> simplify(P)
#Handle12 eliminating f=E redund=1
#Handle12 eliminating e=d redund=2
#Handle12 eliminating e=d redund=2
#Handle12 eliminating h=c redund=3
#Handle12 eliminating g=c redund=4
Presentation: 4 generators, 8 relators of total length 21
Presentation: 4 generators, 7 relators of total length 18
#EliminateGen1 eliminating d=bc
#EliminateGen1 eliminating c=aB
Presentation: 2 generators, 4 relators of total length 14
Presentation: 2 generators, 3 relators of total length 9
    julia> relators(P)
3-element Vector{AbsWord}:
 a²
 b³
 baba
```
Note that the number  of loops over the two phases as well  as the number
of  subword   searches  or  generator  eliminations  in  each  phase  are
determined by a set of option parameters which may  heavily influence the
resulting presentation and the computing time (see Tietze options below).

If "silent" is specified as true, then the printing of the status line by
`Go` in case of T.debug=1 is suppressed.
"""
function Go(P::Presentation,silent=false)
  printstatus=P.debug==1 && !silent
  Check(P)
  Search(P)
  count=0
  while (P.loopLimit==0 || count<P.loopLimit) && !all(isempty,P.relators)
#   println("Go loop")
    SearchEqual(P)
#   @show P.flags
    if P.modified Search(P) end
#   @show P.flags
    EliminateGens(P)
#   print("before search");dump(P)
    if P.modified
      Search(P)
      count+=1
    end
#   print("after search");dump(P)
    if printstatus print("#");PrintStatus(P) end
    if !P.modified break end
  end
  if !all(isempty,P.relators)
    FindCyclicJoins(P)
    if P.modified Search(P) end
    if printstatus print("#");PrintStatus(P) end
  end
  P
end

"""
`GoGo(P)`

`GoGo` performs Tietze transformations on a presentation `P`. It repeatedly
calls  the  `Go`  command  until  neither  the number of generators nor the
number of relators nor the total length of all relators have changed during
five consecutive calls of `Go`.
"""
function GoGo(T::Presentation)
  numgens=length(T.generators)
  numrels=length(T.relators)
  total=totallength(T)
  count=0
  while count<5
    Go(T, T.debug==1)
    count+=1
    if T.debug==1 && (length(T.generators)<numgens || length(T.relators)<numrels)
        PrintStatus(T)
    end
    if length(T.generators)<numgens || length(T.relators)<numrels ||
      totallength(T)<total
      numgens=length(T.generators)
      numrels=length(T.relators)
      total=totallength(T)
      count=0
    end
  end
  if T.debug==1 PrintStatus(T) end
  T
end

"""
OccurrencesPairs(T::Presentation,gen)
counts occurrences of pairs gen,geni or their inverses (where 0<gen<geni)
circularly in the relators or their inverses
returns a 4×numgens matrix; each line is zero until gen+1
the first line  counts the occurrences of (gen,geni)
the second             the occurrences of (gen,-geni)
the third              the occurrences of (-gen,geni)
the last               the occurrences of (-gen,-geni)
"""
function OccurrencesPairs(T::Presentation,gen)
  res=fill(0,4,length(T.generators))
  for rel in T.relators, i in eachindex(rel)
    if rel[i]==gen add=0
    elseif rel[i]==-gen add=2
    else continue
    end
    geni=rel[i==end ? 1 : i+1]
    if geni>gen res[1+add,geni]+=1
    elseif geni<-gen res[2+add,-geni]+=1
    end
    geni=rel[i==1 ? end : i-1]
    if geni>gen res[4-add,geni]+=1
    elseif geni<-gen res[3-add,-geni]+=1
    end
  end
  res
end

"""
MostFrequentPairs(presentation, n) ....  occurrences of pairs

`MostFrequentPairs`  returns a list  describing the  n  most frequently
occurring relator subwords of the form  g1*g2,  where  g1  and  g2 are
different generators or their inverses. n=0 interpreted as infinity
"""
function MostFrequentPairs(T::Presentation, nmax::Int)
  if nmax<0 error("second argument must be≥0") end
  gens=T.generators
  numgens=length(T.generators)
  pairs=Vector{Int}[]
  n=0
  if nmax==1
    max=0
    for i in 1:numgens-1
      occlist=OccurrencesPairs(T, i)
      for j in i+1:numgens, k in 4:-1:1
        if occlist[k,j]>=max
          max=occlist[k,j]
          if isempty(pairs) push!(pairs,[max, i, j, k-1])
          else pairs[1]=[max, i, j, k-1]
          end
          n=1
        end
      end
    end
  else
    for i in 1:numgens-1
      occlist=OccurrencesPairs(T, i)
      for j in i+1:numgens, k in 1:4
        if occlist[k,j]>0
          n+=1
          if length(pairs)==n-1 push!(pairs,[occlist[k,j], i, j, k-1])
          else pairs[n]=[occlist[k,j], i, j, k-1]
          end
        end
      end
      if nmax!=0 && n>nmax
        sort!(pairs,rev=true)
        resize!(pairs,nmax)
        n=nmax
      end
    end
    sort!(pairs,rev=true)
  end
  pairs
end

"""
`PrintOptions(P)`

Several   of  the  Tietze  transformation   commands  described  above  are
controlled  by certain parameters, the *Tietze options*, which often have a
tremendous  influence on  their performance  and results.  However, in each
application  of  the  commands,  an  appropriate  choice  of  these  option
parameters  will depend  on the  concrete presentation under investigation.
Therefore  we have implemented the  Tietze options in such  a way that they
are  associated to the presentation: Each presentation keeps its own set of
Tietze option parameters in the form of ordinary fields. In particular, you
may  alter the value of any of these Tietze options by just assigning a new
value to the respective component.

`PrintOptions`  prints the Tietze  option fields of  the  specified
presentation `P`.

The Tietze options have the following meaning.

`protected`:
  The  first `P.protected` generators  in a presentation  `P` are protected
  from  being eliminated by the Tietze transformations functions. There are
  only two exceptions: The option `P.protected` is ignored by the functions
  `Eliminate(P,gen)` and `Substitute(P,n,eliminate)`    because   they
  explicitly  specify the generator to be  eliminated. The default value of
  `protected` is 0.

`eliminationsLimit`:
  Whenever  the  elimination  phase  of  the  `Go` command is entered for a
  presentation  `P`, then  it will  eliminate at most `P.eliminationsLimit`
  generators (except for further ones which have turned out to be trivial).
  Hence  you may use the `eliminationsLimit` parameter as a break criterion
  for  the  `Go`  command.  Note,  however,  that  it  is  ignored  by  the
  `Eliminate` command. The default value of `eliminationsLimit` is 100.

`expandLimit`:
  Whenever the routine for eliminating more than 1 generators is called for
  a presentation `P` by the `Eliminate` command or the elimination phase of
  the  `Go` command, then it saves the  given total length of the relators,
  and  subsequently it  checks the  current total  length against its value
  before  each elimination. If the total  length has increased to more than
  `P.expandLimit`  per cent of its original value, then the routine returns
  instead   of  eliminating  another  generator.  Hence  you  may  use  the
  `expandLimit`  parameter as a  break criterion for  the `Go` command. The
  default value of `expandLimit` is 150.

`generatorsLimit`:
  Whenever  the  elimination  phase  of  the  `Go` command is entered for a
  presentation  `P` with  `n` generators,  then it  will eliminate  at most
  `n-P.generatorsLimit` generators (except for generators which turn out to
  be trivial). Hence you may use the `generatorsLimit` parameter as a break
  criterion for the `Go` command. The default `generatorsLimit` is 0.

`lengthLimit`:
  The  Tietze transformation commands will never eliminate a generator of a
  presentation  `P`,  if  they  cannot  exclude  the  possibility  that the
  resulting   total   length   of   the   relators  exceeds  the  value  of
  `P.lengthLimit`. The default value of `lengthLimit` is `infinity`.

`loopLimit`:
  Whenever  the `Go` command is called for a presentation `P`, then it will
  loop over at most `P.loopLimit` of its basic steps. Hence you may use the
  `loopLimit`  parameter as  a break  criterion for  the `Go`  command. The
  default value of `loopLimit` is `infinity`.

`debug`:
  Whenever Tietze transformation commands are called for a presentation `P`
  with  `P.debug` `= 0`, they will not  provide any output except for error
  messages. If `P.debug` `= 1`, they will display some reasonable amount of
  output  which allows you to watch the  progress of the computation and to
  decide  about your next commands.  In the case `P.debug`  `= 2`, you will
  get  a much more generous amount of  output. Finally, if `P.debug` `= 3`,
  various  messages on internal details will be added. The default value of
  `debug` is 1.

`saveLimit`:
  Whenever  the  `Search`  command  has  finished  its  main  loop over all
  relators  of a presentation `P`, then  it checks whether during this loop
  the   total  length  of  the  relators  has  been  reduced  by  at  least
  `P.saveLimit`  per cent. If  this is the  case, then `Search` repeats its
  procedure  instead  of  returning.  Hence  you  may  use  the `saveLimit`
  parameter  as  a  break  criterion  for  the  `Search`  command  and,  in
  particular,  for the search phase of  the `Go` command. The default value
  of `saveLimit` is 10.

`searchSimultaneous`:
  Whenever  the  `Search`  or  the  `SearchEqual`  command  is called for a
  presentation `P`, then it is allowed to handle up to `searchSimultaneous`
  short  relators simultaneously (see  for the description  of the `Search`
  command  for  more  details).  The  choice  of this parameter may heavily
  influence  the performance as well as the  result of the `Search` and the
  `SearchEqual`  commands and  hence also  of the  search phase of the `Go`
  command. The default value of `searchSimultaneous` is 20.

As soon as a presentation has been defined, you may alter any of its Tietze
option  parameters  at  any  time  by  just  assigning  a  new value to the
respective component.

To  demonstrate the  effect of  the `eliminationsLimit`  parameter, we will
give  an example in which we  handle a subgroup of index  240 in a group of
order  40320  given  by  a  presentation  due  to  B.~H.  Neumann. First we
construct  a presentation  of the  subgroup, and  then we  apply to  it the
`GoGo`  command for  different values  of the `eliminationsLimit` parameter
(including  the  default  value  100).  In  fact, we also alter the `debug`
parameter,  but this is only done in  order to suppress most of the output.
In  all cases  the resulting  presentations cannot  be improved any more by
applying the `GoGo` command again, i.e., they are the best results which we
can get without substituting new generators.
```julia-rep1
    gap> F3 := FreeGroup( "a", "b", "c" );;
    gap> G := F3 / [ F3.1^3, F3.2^3, F3.3^3, (F3.1*F3.2)^5,
            (F3.1^-1*F3.2)^5, (F3.1*F3.3)^4, (F3.1*F3.3^-1)^4,
            F3.1*F3.2^-1*F3.1*F3.2*F3.3^-1*F3.1*F3.3*F3.1*F3.3^-1,
            (F3.2*F3.3)^3, (F3.2^-1*F3.3)^4 ];;
    gap> H := Subgroup( G, [ G.1, G.3 ] );;
    gap> P := PresentationSubgroup( G, H );
    << presentation with 224 gens and 593 rels of total length 2769 >>
    julia> for i in [ 28, 29, 30, 94, 100 ]
    Pi=deepcopy(P)
    Pi.eliminationsLimit=i
    println("# eliminationsLimit set to ",i)
      Pi.debug=0
      GroupPresentations.GoGo(Pi)
      print(Pi)
    end
    #  eliminationsLimit set to 28
    #  there are 2 generators and 95 relators of total length 10817
    #  eliminationsLimit set to 29
    #  there are 2 generators and 5 relators of total length 35
    #  eliminationsLimit set to 30
    #  there are 3 generators and 98 relators of total length 2928
    #  eliminationsLimit set to 94
    #  there are 4 generators and 78 relators of total length 1667
    #  eliminationsLimit set to 100
    #  there are 3 generators and 90 relators of total length 3289 |
```

Similarly,  we demonstrate  the influence  of the  `saveLimit` parameter by
just  continuing the  preceding example  for some  different values  of the
`saveLimit`  parameter  (including  its  default  value  10),  but  without
changing  the `eliminationsLimit`  parameter which  keeps its default value
100.

```julia-rep1
    julia> for i in [ 9, 10, 11, 12, 15 ]
       Pi=deepcopy(P)
       Pi.saveLimit=i
       println("# saveLimit set to ",i)
       Pi.debug=0
       GoGo(Pi)
       print(Pi)
    end
    &I  saveLimit set to 9
    &I  there are 3 generators and 97 relators of total length 5545
    &I  saveLimit set to 10
    &I  there are 3 generators and 90 relators of total length 3289
    &I  saveLimit set to 11
    &I  there are 3 generators and 103 relators of total length 3936
    &I  saveLimit set to 12
    &I  there are 2 generators and 4 relators of total length 21
    &I  saveLimit set to 15
    &I  there are 3 generators and 143 relators of total length 18326 |
```
"""
function PrintOptions(P::Presentation)
  names=[:eliminationsLimit,:expandLimit,:generatorsLimit,:lengthLimit,
         :loopLimit,:debug,:saveLimit,:searchSimultaneous]
  lst=sort(collect(intersect(keys(P.prop),names)))
  len=maximum(length.(string.(lst)))
  foreach(lst)do nam
    println(rpad(nam,len),"= ",getproperty(P,nam))
  end
end

"""
`RemoveGenerators(P::presentation)`

deletes  the redundant  Tietze generators  and renumbers  the non-redundant
ones  accordingly. The  redundant generators  are assumed  to be  marked in
`P.inverses` by an entry `P[i]!=i`.
"""
function RemoveGenerators(P::Presentation)
  debug(P,3,"RemoveGenerators renumbering the Tietze generators")
  if P.numredunds==0 return end
  ngens=length(P.generators)
  j=0
# @show j,P.inverses, P.numredunds
  for i in 1:ngens
    if P[i]==i
      j+=1
      if j<i
        P[i]=j
        if P[-i]>0 P[-i]=j else P[-i]=-j end
      end
    else
      P[i]=P[-i]=0
    end
  end
# @show j,P.inverses
  if j!=ngens-P.numredunds
    error("j=$j ngens=$ngens redunds=$(P.numredunds). You should never get here.\n")
  end
# @show j, P.inverses, ngens
  for rel in P.relators, i in eachindex(rel) rel[i]=P[rel[i]] end
  if haskey(P, :imagesOldGens)
    for i in 1:length(P.imagesOldGens)
      image=P.imagesOldGens[i]
      newim=Int[]
      for j=1:length(image) push!(newim, P[image[j]]) end
      P.imagesOldGens[i]=reduceword(newim)
    end
  end
  for i in 1:ngens
    j=P[i]
    if j<i && j>0
      P.generators[j]=P.generators[i]
      P[j]=j
      P[-j]=P[-i]
      if haskey(P, :preImagesNewGens)
          P.preImagesNewGens[j]=P.preImagesNewGens[i]
      end
    end
  end
  j=ngens
  ngens1=ngens+1
  ngens-=P.numredunds
  P.inverses=P.inverses[ngens1.+(-ngens:ngens)]
  resize!(P.generators,ngens)
  if haskey(P, :preImagesNewGens) resize!(P.preImagesNewGens,ngens) end
  P.numredunds=0
end

# hash of cyclic subword of relator w starting at i and of length l
# i in 1:length(w) or i negative: backward subword
function hashword(w,i,l,T)
  res=UInt(0)
  if i<0
    for j in -i:-1:max(1,-i+1-l) res=hash(T[-w[j]],res) end
    for j in length(w):-1:length(w)-l-i+1 res=hash(T[-w[j]],res) end
  else
    for j in i:min(i+l-1,length(w)) res=hash(T[w[j]],res) end
    for j in 1:i+l-1-length(w) res=hash(T[w[j]],res) end
  end
  res
end

# finds intersection of 2 sorted lists; returns list of pairs (i1,i2)
# such that by(a[i1])==by(b[i2])
function commonsorted(a,b;by=x->x)
  res=Tuple{Int,Int}[]
  ai=bi=1
  ri=0
@inbounds while ai<=length(a) && bi<=length(b)
    c=cmp(by(a[ai]),by(b[bi]))
    if     c>0 bi+=1
    elseif c<0 ai+=1
    else push!(res,(ai,bi)); ai+=1; bi+=1
    end
  end
  res
end

# perhaps w1[abs(i1)]*sign(i1)=w2[abs(i2)]*sign(i2), extend the match
# as much as possible; returns (new i1, new i2, length of found match)
function maxmatch(w1,i1,w2,i2,T)
  l=0;lw=length(w1)
  prev(i,w)=i>0 ? (i==1 ? length(w) : i-1) : (i==-length(w) ? -1 : i-1)
  next(i,w)=i>0 ? (i==length(w) ? 1 : i+1) : (i==-1 ? -length(w) : i+1)
  eq(i1,i2)=T[w1[abs(i1)]*sign(i1)]==T[w2[abs(i2)]*sign(i2)]
  n1=i1;n2=i2
  while l<lw
    if !eq(n1,n2) break end
    j1=n1;j2=n2;l+=1
    n1=next(j1,w1);n2=next(j2,w2)
  end
  if l>0
  while l<lw
    n1=prev(i1,w1);n2=prev(i2,w2)
    if !eq(n1,n2) break end
    i1=n1;i2=n2;l+=1
  end
  end
  (i1,i2,l)
end

# cyclic subword of relator w starting at i and of length l; takes i mod
# length(w) but remembers sign (positive: forward word, negative: -backward)
function subword(w,i,l)
  neg=sign(i)
  i=mod1(abs(i),length(w))
  if neg<0 vcat(-w[i:-1:max(1,i+1-l)],-w[length(w):-1:length(w)-l+i+1])
  else vcat(w[i:min(i+l-1,length(w))],w[1:i+l-1-length(w)])
  end
end

# doreplace(w1,w2,pos1,pos2,l) replace in relator w2 substring of length l
# starting at pos2 with complement in w1 of identical substring of length l
# starting at pos1 in w1.
# return new w2
function doreplace(w1,w2,pos1,pos2,l)
  l1=length(w1)
  compl=subword(w1,pos1>0 ? mod(1-pos1,l1)-l1 : mod1(-pos1+1,l1),l1-l)
  vcat(compl,subword(w2,pos2+l,length(w2)-l))
# circshift(vcat(compl,subword(w2,pos2+l,length(w2)-l)),pos2-1)
end

# search from i-th relator, hashing subwords until j-th
function SearchC(T::Presentation,i,j,equal=false)
# @show i,j,equal
  rels=T.relators
# println("SearchC(",join(AbsWord.(rels[i:j],Ref(T)),","),",",equal,")")
  len=length(rels[i])
  if !equal l=div(len,2)+1;lmin=len-len&1;lmax=lmin+1
  elseif isodd(len) error("searchequal and odd length")
  else l=div(len,2);lmin=lmax=len
  end
# @show l,lmin,lmax
  hh=NamedTuple{(:hash, :pos, :relno), Tuple{UInt, Int, Int}}[]
# if T.inverses!=T.numgens:-1:-T.numgens error("numgens=",T.numgens," invs=",
#                    T.inverses) end
  altered=Int[]
  for k in i:length(rels)
    if T.flags[k]==3 continue end # keep involution witness
    rel=rels[k]
    if length(rel)<lmin continue end # sort should be called before SearchC
    if length(rel)!=2 || rel[1]!=rel[2] rel=reduceword(rel,T) end
    while true
      hk=map(u->(hash=hashword(rel,u,l,T),pos=u),eachindex(rel))
      sort!(hk,by=x->x.hash)
      pp=commonsorted(hh,hk;by=x->x.hash)
      if isempty(pp) break end
      hk=map(pp)do (p1,p2)
        (pos1,pos2,w)=(hh[p1].pos,hk[p2].pos,hh[p1].relno)
        (match=maxmatch(rels[w],pos1,rel,pos2,T),relno1=w,relno2=k)
      end
#     @show hk
      (pos1,pos2,l1),w,k=hk[argmax(map(x->(x[1][3],-x[1][2]),hk))]
#     if l1<l println("!!!!hash",hk) end
      if l1<l error("!")
      else
#      print(rio(),"$k:",AbsWord(rel,T),"[$pos2:$(pos2+l1-1)]->$w:",AbsWord(rels[w],T),"\\[$pos1:$(pos1+l1-1)]")
        rel=doreplace(rels[w],rel,pos1,pos2,l1)
#       println(rio(),"=>",AbsWord(rel,T))
      end
      rel=reduceword(rel,T)
#     println("=>",AbsWord(rel,T))
#     if T.flags[k]==1 T.flags[k]=2 end
#     if equal T.flags[k]=1 end
      if equal || length(rel)<l break end
    end
    if k<=j && length(rel) in lmin:lmax
      for u in 1:length(rel)
        push!(hh,(hash=hashword(rel,u,l,T),pos=u,relno=k))
        push!(hh,(hash=hashword(rel,-u,l,T),pos=-u,relno=k))
      end
      sort!(hh,by=x->x.hash)
    end
    if rel!=rels[k]
      push!(altered,k)
      rels[k]=rel
    end
  end
  altered
end

"""
`Search(P)`

`Search`  performs Tietze transformations on a presentation `P`. It tries
to  reduce the relator lengths by  substituting common subwords of relators
by shorter words.

The idea is to find pairs of relators `r₁` length `l₁` and `r₂` length `l₂`
with  `l₁≤l₂` such that `r₁` and `r₂` coincide (possibly after inverting or
conjugating  one of them)  on some maximal  subword `w`, of length `>l₁/2`,
and  then to substitute each copy of  `w` in `r₂` by the inverse complement
of `w` in `r₁`.

Two  of the Tietze  option parameters which  are listed at  the end of this
section  may  strongly  influence  the  performance  and the results of the
`Search`   command.   These   are   the   parameters   `P.saveLimit`  and
`P.searchSimultaneous`. The first of them has the following effect.

When  Search  has  finished  its  main  loop  over all relators, then, in
general,  there are relators which have changed and hence should be handled
again in another run through the whole procedure. However, experience shows
that  it really does  not pay to  continue this way  until no more relators
change.  Therefore,  `Search`  starts  a  new  loop only if the loop just
finished  has  reduced  the  total  length  of  the  relators  by  at least
`P.saveLimit` per cent.

The default value of `P.saveLimit` is 10.

To  understand the effect of  the parameter `P.searchSimultaneous`, we have
to look in more detail at how `Search` proceeds.

First,  it  sorts  the  list  of  relators  by  increasing lengths. Then it
performs  a loop  over this  list. In  each step  of this loop, the current
relator  is treated  as *short  relator* `r₁`,  and a  subroutine is called
which  loops over the succeeding relators, treating them as *long relators*
`r₂` and performing the respective comparisons and substitutions.

As  this  subroutine  performs  a  very  expensive  process,  it  has  been
implemented as SearchC. For the given relator `r₁` of length `l₁`, it first
determines  the *minimal  match length*  `l=⌊l₁/2⌋+1`. Then  it builds up a
hash  list for all  subwords of length  `l` occurring in  the conjugates of
`r₁`  or  `r₁⁻¹`,  and  finally  it  loops  over all long relators `r₂` and
compares the hash values of their subwords of length `l` against this list.
A  true comparison is  only done if a hash match has been found.

To improve the efficiency of this process we allow the subroutine to handle
several  short  relators  simultaneously  provided  that they have the same
minimal  match  length.  If,  for  example,  it  handles `n` short relators
simultaneously,  then you save `n-1` loops over the long relators `r₂`, but
you pay for it by additional fruitless subword comparisons. In general, you
will  not get the best performance  by always choosing the maximal possible
number of short relators to be handled simultaneously. In fact, the optimal
choice  of  the  number  will  depend  on  the  concrete presentation under
investigation.   You  can  use   the  parameter  `P.searchSimultaneous`  to
prescribe  an upper bound  for the number  of short relators  to be handled
simultaneously.

The default value of `P.searchSimultaneous` is 20.
"""
function Search(P::Presentation)
  debug(P,3,"Search subwords")
  Check(P)
  T=P
  rels=T.relators
  T.modified=false
  save=P.saveLimit//100
  while totallength(T)>0
    sort!(P)
    modified=false
    oldtotal=totallength(T)
    flag=0
    for i in eachindex(T.relators)
      if T.flags[i]>1 || isempty(rels[i]) continue end
      len=length(rels[i])
      lmax=iseven(len) ? len+1 : len
      if flag<T.flags[i] flag=T.flags[i] end
      simultan=1
      j=i
      lastj=0
      k=i+1
      while k<=length(T.relators) && length(rels[k])<=lmax && simultan<P.searchSimultaneous
        if T.flags[k]<=1
          lastj=j
          j=k
          simultan+=1
        end
        k+=1
      end
      while k<=length(T.relators) && (T.flags[k]>1 || (flag==0 && T.flags[k]==0))
        k+=1
      end
      if k>length(T.relators) j=lastj end
      if i<=j
#       @show T.flags,i,j
        altered=SearchC(T,i,j)
        if !isempty(altered)
          debug(P,3,"SearchC($i,$j): ",altered," altered")
#         dump(P)
          modified=true
          T.flags[altered]=map(x->x==3 ? x : 2,T.flags[altered])
 #        @show T.flags,altered
        end
        i=j
      end
    end
    for i in eachindex(T.flags) if T.flags[i] in 1:2 T.flags[i]-=1 end end
    if modified
     if totallength(T)<oldtotal
        T.modified=true
#       print("in search avant handle");display_balanced(P)
        HandleLength1Or2Relators(P)
#       print("in search apres handle");display_balanced(P)
        sort!(P)
        PrintStatus(P,2)
      end
    end
    if !modified break end
  end
end

"""
`SearchEqual( P )`

`SearchEqual`  performs Tietze transformations on a  presentation  `P`.
It tries to alter relators by substituting common subwords of relators by
subwords of equal length.

The  idea is  to find  pairs of  relators `r₁`  length `l₁` and `r₂` length
`l₂`, such that `l₁` is even, `l₁≤l₂`, and `r₁` and `r₂` coincide (possibly
after  inverting or conjugating one of them) in some maximal subword `w` of
length  `≥l₁/2`. Let `l` be the length  of `w`. Then, if `l>l₁/2`, the pair
is  handled  as  in  `Search`.  Otherwise,  if `l=l₁/2`, then `SearchEqual`
substitutes  each copy of `w`  in `r₂` by the  inverse complement of `w` in
`r₁`.

The  Tietze  option   parameter   `P.searchSimultaneous`  is   used  by
`SearchEqual` in the same way as described for `Search`.

However, `SearchEqual` does  not use the  parameter `P.saveLimit`:
The loop over the relators is executed exactly once.
"""
function SearchEqual(P::Presentation)
  debug(P,3,"SearchEqual length subwords")
  Check(P)
  T=P
  simultanlimit=P.searchSimultaneous
  sort!(P)
  rels=T.relators
  modified=false
  oldtotal=totallength(T)
  i=1
  while i<length(T.relators)
    leng=length(rels[i])
    if leng>3 && iseven(leng)
      simultan=1
      j=i
      lastj=0
      k=i+1
      while k<=length(T.relators) && (length(rels[k])<=leng && simultan<simultanlimit)
        if length(rels[k])==leng
          lastj=j
          j=k
          simultan+=1
        end
        k+=1
      end
      while k<=length(T.relators) && length(rels[k])<leng k+=1 end
      if k>length(T.relators) j=lastj end
      if i<=j
#       display_balanced(P,true)
        altered=SearchC(T, i, j, true)
        if !isempty(altered)
          debug(P,3,"SearchCeq($i,$j) altered=",altered)
          T.flags[altered].=1
          modified=true
        end
        i=j
      end
    end
    i+=1
  end
  if modified
   if totallength(T)<oldtotal
      T.modified=true
      HandleLength1Or2Relators(P)
      sort!(P)
      PrintStatus(P,2)
    end
  end
end

"""
`Substitute(P::Presentation,n=1,eliminate=0)`

Basically,  it substitutes a squarefree word of length 2 as a new generator
and  then eliminates a generator from  the extended generator list. We will
describe this process in more detail.

`n` is expected to be a positive integer, and `eliminate` is expected to be
0, 1, or 2.

`Substitute`  first determines the `n` most frequently occurring squarefree
relator  subwords  of  length  2  and  sorts  them by decreasing numbers of
occurrences.  Let `ab` be the  `n`th word in that  list, and let `x` be the
first  available name for a new generator (see `add_generator` for details)
Then  `Substitute` defines a new generator `x`, adds it to the presentation
together  with a new relator `x\\ab`,  and replaces all occurrences of `ab`
in the given relators by `x`.

Finally,  it eliminates some generator  from  the  extended presentation.
The  choice  of that  generator  depends  on  the  actual  value  of  the
`eliminate` parameter:

If `eliminate` is zero, then the generator to be eliminated is  chosen as
by the  `Eliminate` command.  This means  that in this case it may well
happen that  it is the generator `x` just introduced which  is  now
deleted  again  so  that  you do  not  get  any  remarkable  progress  in
transforming  your  presentation.   On  the other  hand,  this  procedure
guaranties that the total length of the relators will not be increased by
a call of `Substitute` with `eliminate = 0`.

Otherwise,  if  `eliminate`  is  1  or  2, then `Substitute` eliminates the
respective factor of the substituted word `ab`, i.e., `a` for `eliminate=1`
or  `b` for `eliminate=2`. In this case,  it may well happen that the total
length  of  the  relators  increases,  but  sometimes  such an intermediate
extension is the only way to finally reduce a given presentation.

In  order to decide which arguments might  be appropriate for the next call
of  `Substitute`, it is helpful to print  out a list of the most frequently
occurring  squarefree relator subwords of length 2. The `showpairs` command
does this.

As an example we handle a subgroup of index 266 in the Janko group `J₁`.

|    gap> F2 := FreeGroup( "a", "b" );;
    gap> J1 := F2 / [ F2.1^2, F2.2^3, (F2.1*F2.2)^7,
         Comm(F2.1,F2.2)^10, Comm(F2.1,F2.2^-1*(F2.1*F2.2)^2)^6 ];;
    gap> a := J1.1;;  b := J1.2;;
    gap> H := Subgroup ( J1, [ a, b^(a*b*(a*b^-1)^2) ] );;
    gap> P := PresentationSubgroup( J1, H );
```
julia> P=Presentation("
1: a=a⁻
2: r=r⁻
3: b=b⁻
4: g=g⁻
5: p=p⁻
6: q=q⁻
7: l=l⁻
8: j=j⁻
9: ied=1
10: c³=1
11: dca=1
12: i³=1
13: e³=1
14: tao=1
15: fac=1
16: w³=1
17: vnl=1
18: mh⁻j=1
19: utf=1
20: wup=1
21: m⁻kg=1
22: snq=1
23: ts⁻f=1
24: v³=1
25: gc=mn⁻
26: rg=ks⁻
27: lo=op⁻
28: gb=q⁻g⁻
29: h⁻f=fi
30: bh=do⁻
31: gj=l⁻g⁻
32: lb=v⁻f
33: gr=r⁻g⁻
34: gc⁻beb=1
35: we²du⁻=1
36: rhn²g=1
37: u⁵=1
38: kop=ha⁻b⁻
39: gal=feb⁻
40: gsf=sfb⁻
41: agf=l⁻b⁻e⁻
42: ahn=ib⁻g⁻
43: hab=kop⁻
44: wu⁻d⁻bgbf=1
45: jc⁻abed⁻b=1
46: wuo⁻nfd⁻u⁻=1
47: kh⁻bdi²b=1
48: jbde⁻bac=1
49: gbca=h⁻b⁻a⁻c
50: qlaf=fi⁻u⁻p⁻
51: babc=e⁻cb⁻a⁻
52: abac=cb⁻a⁻b⁻
53: k⁻gh⁻m=m⁻hg⁻k
54: rbed=fdf⁻t⁻
55: km⁻ha=r⁻b⁻df⁻
56: f⁻afp=p⁻f⁻a⁻f
57: rlns=hk⁻r⁻h
58: grgs=sl⁻q⁻l⁻
59: f⁻h⁻be⁻df⁻acb=1
60: sfo⁻fi⁻f⁻ah⁻g=1
61: vlgbac²sl=1
62: vbi⁻e²brgs=1
63: bdi⁻f⁻a=kfw⁻uf⁻
64: bgbdi=di²a⁻i⁻
65: pu⁻th⁻m=u⁻th⁻mg⁻
66: ababab=b⁻a⁻b⁻a⁻b⁻a⁻
67: lgrkgl=v⁻g⁻r⁻g⁻sn
68: at⁻f⁻aft=t⁻f⁻a⁻fta⁻
69: lnqn⁻ln=nq⁻n⁻l⁻nq⁻
70: gqgbdp=r⁻g⁻q⁻g⁻b⁻d
71: jrglql=mg⁻m⁻j⁻r⁻g⁻
72: bhgm⁻jmg=df⁻l⁻fd⁻b⁻h
73: km⁻jmk⁻rkm⁻j=r⁻km⁻j⁻mk⁻r⁻km⁻
74: jmgm⁻jmgm⁻jm=mg⁻m⁻j⁻mg⁻m⁻j⁻mg⁻")
Presentation: 23 generators, 74 relators of total length 461
|    gap> GoGo( P );
    &I  there are 3 generators and 47 relators of total length 1368
    &I  there are 2 generators and 46 relators of total length 3773
    &I  there are 2 generators and 46 relators of total length 2570
    gap> GoGo( P );
    &I  there are 2 generators and 46 relators of total length 2568
    gap> GoGo( P );
    gap> & We do not get any more progress without substituting a new
    gap> & generator
    gap> Substitute( P );
    &I  substituting new generator _x28 defined by _x6*_x23^-1
    &I  eliminating _x28 = _x6*_x23^-1
    gap> & GAP cannot substitute a new generator without extending the
    gap> & total length, so we have to explicitly ask for it
    gap> showpairs( P );
    &I  1.  504  occurrences of  _x6 * _x23^-1
    &I  2.  504  occurrences of  _x6^-1 * _x23
    &I  3.  448  occurrences of  _x6 * _x23
    &I  4.  448  occurrences of  _x6^-1 * _x23^-1
    gap> Substitute( P, 2, 1 );
    &I  substituting new generator _x29 defined by _x6^-1*_x23
    &I  eliminating _x6 = _x23*_x29^-1
    &I  there are 2 generators and 46 relators of total length 2867
    gap> GoGo( P );
    &I  there are 2 generators and 45 relators of total length 2417
    &I  there are 2 generators and 45 relators of total length 2122
    gap> Substitute( P, 1, 2 );
    &I  substituting new generator _x30 defined by _x23*_x29^-1
    &I  eliminating _x29 = _x30^-1*_x23
    &I  there are 2 generators and 45 relators of total length 2192
    gap> GoGo( P );
    &I  there are 2 generators and 42 relators of total length 1637
    &I  there are 2 generators and 40 relators of total length 1286
    &I  there are 2 generators and 36 relators of total length 807
    &I  there are 2 generators and 32 relators of total length 625
    &I  there are 2 generators and 22 relators of total length 369
    &I  there are 2 generators and 18 relators of total length 213
    &I  there are 2 generators and 13 relators of total length 141
    &I  there are 2 generators and 12 relators of total length 121
    &I  there are 2 generators and 10 relators of total length 101
    gap> showpairs( P );
    &I  1.  19  occurrences of  _x23 * _x30^-1
    &I  2.  19  occurrences of  _x23^-1 * _x30
    &I  3.  14  occurrences of  _x23 * _x30
    &I  4.  14  occurrences of  _x23^-1 * _x30^-1
    gap> & If we save a copy of the current presentation, then later we
    gap> & will be able to restart the computation from the current state
    gap> P1 := copy( P );;
    gap> & Just for demonstration, let's make an inconvenient choice
    gap> Substitute( P, 3, 1 );
    &I  substituting new generator _x31 defined by _x23*_x30
    &I  eliminating _x23 = _x31*_x30^-1
    &I  there are 2 generators and 10 relators of total length 122
    gap> GoGo( P );
    &I  there are 2 generators and 9 relators of total length 105
    gap> & The presentation is worse than the one we have saved, so let's
    gap> & restart from that one again
    gap> P := copy( P1 );
    << presentation with 2 gens and 10 rels of total length 101 >>
    gap> Substitute( P, 2, 1);
    &I  substituting new generator _x31 defined by _x23^-1*_x30
    &I  eliminating _x23 = _x30*_x31^-1
    &I  there are 2 generators and 10 relators of total length 107
    gap> GoGo( P );
    &I  there are 2 generators and 9 relators of total length 84
    &I  there are 2 generators and 8 relators of total length 75
    gap> Substitute( P, 2, 1);
    &I  substituting new generator _x32 defined by _x30^-1*_x31
    &I  eliminating _x30 = _x31*_x32^-1
    &I  there are 2 generators and 8 relators of total length 71
    gap> GoGo( P );
    &I  there are 2 generators and 7 relators of total length 56
    &I  there are 2 generators and 5 relators of total length 36
    gap> relators(P);
    &I  1. _x32^5
    &I  2. _x31^5
    &I  3. _x31^-1*_x32^-1*_x31^-1*_x32^-1*_x31^-1*_x32^-1
    &I  4. _x31*_x32*_x31^-1*_x32*_x31^-1*_x32*_x31*_x32^-2
    &I  5. _x31^-1*_x32^2*_x31*_x32^-1*_x31^2*_x32^-1*_x31*_x32^2 |
```

As shown in the preceding example, you can use the `copy` command to save
a copy of a presentation and to restart from it  again if you want
to try an alternative strategy.
"""
function Substitute(P::Presentation,n::Int=1,elim::Int=0)
  if !(0<=elim<=2) error("last argument must be 0, 1, or 2") end
  T=P
  if length(T.generators)<2 return end
  pairs=MostFrequentPairs(P, n)
  if length(pairs)<n
    debug(P,1,"Substitute: there are less than $n pairs available")
    return
  end
  m,i,j,k=pairs[n]
  if k>1 i=T[-i] end
  if isodd(k) j=T[-j] end
  pair=[i, j]
  gen=newgen(P)
  word=T(pair...)
  debug(P,1,"Substitute new generator ",gen," defined by ", word)
  add_relator(P,gen\word;stop=false)
  UpdateGeneratorImages(P, pair)
  sdebug=P.debug;P.debug=0;Search(P);
  P.debug=2
  if elim>0 EliminateGen(P, abs(pair[elim])) else EliminateGen1(P) end
  P.debug=sdebug
  if T.modified Search(P) end
  if T.numredunds>0 RemoveGenerators(P) end
  PrintStatus(P,1)
end

"""
`find_exponents(presentation)`  tries  to  find  exponents  for  the Tietze
generators  to  use  in  SubstituteCyclicJoins  and  returns them in a list
parallel to the list of the generators.
"""
function find_exponents(T::Presentation)
  debug(T,3,"find generator exponents")
  rels=T.relators
  exponents=fill(0,length(T.generators))
  for i in eachindex(rels)
    if isempty(rels[i]) || !allequal(rels[i]) continue end
    if rels[i][1]<0 rels[i]=-rels[i] end
    num=rels[i][1]
    exponents[num]=exp=gcd(exponents[num], length(rels[i]))
    if exp<length(rels[i])
      rels[i]=fill(num,exp)
      T.flags[i]=1
      T.modified=true
    end
  end
  exponents
end

"""
`SubstituteCyclicJoins(presentation)`

`SubstituteCyclicJoins`    performs   Tietze   transformations   on   a
presentation `P`.  It tries to find pairs of generators `a` and `b`, say,
for which among  the  relators (possibly after  inverting  or conjugating
some of them) there are the commutator `[a,b]` and powers `a^m` and `b^n`
with mutually  prime exponents  `m` and  `n`.   For each  such  pair,  it
substitutes the product `ab`  as  a new generator, and then it eliminates
the generators `a` and `b`.
"""
function SubstituteCyclicJoins(P::Presentation)
  debug(P,3,"SubstituteCyclicJoins")
  exponents=find_exponents(P)
  P.modified=false
  if iszero(exponents) return end
  sort!(P)
  gens=P.generators
  rels=P.relators
  i=1
  while i<=length(P.relators) && length(P.relators[i])<=4
    rel=rels[i]
    if length(P.relators[i])==4 && rel[1]==P[-rel[3]] && rel[2]==P[-rel[4]]
      num1=abs(rel[1])
      exp1=exponents[num1]
      num2=abs(rel[2])
      exp2=exponents[num2]
      if exp1>0 && exp2>0 && gcd(exp1, exp2)==1
        gen=newgen(P)
        gen2=gens[num2]
        debug(P,1,"I substituting new generator ",gen," defined by ",
                gens[num1]*gen2)
        add_relator(P, gens[num1]/gen^exp2;stop=false)
        add_relator(P, gen2/gen^exp1;stop=false)
        UpdateGeneratorImages(P,[num1, num2])
        debug=P.debug;P.debug=2;
        Eliminate(P,gens[num1]);Eliminate(P,gen2);
        P.debug=debug
        P.modified=true
        i=0
      end
    end
    i+=1
  end
  if P.modified
    HandleLength1Or2Relators(P)
    sort!(P)
    PrintStatus(P,1)
  end
end

"""
`Substitute(P::Presentation,word[,string::String])`

`word`  is be either an abstract word or a Tietze word in the generators of
`P`.  It substitutes the given word as a new generator of `P`. This is done
as follows.

First,  `Substitute` creates a new abstract generator, `g` say, and adds it
to the presentation `P`, then it adds a new relator `g\\word` to `P`.

If  a `string` has been specified as  third argument, the new generator `g`
will  be named by `string`,  otherwise it will get  a default name as
described in the function `add_generator`.

More precisely: If, for instance, `word` is an abstract word, a call

`Substitute( P, word )`

is more or less equivalent to

`add_generator(P);add_relator(P,P.generators[end]\\word)|`

whereas a call `Substitute(P,word,string)`
is more or less equivalent to

`g=AbsWord(string);add_generator(P,g);add_relator(P,g\\word )`

The  essential difference is, that `Substitute`, as a Tietze transformation
of  `P`, saves and updates  the lists of generator  images and preimages if
they  are being traced under the Tietze transformations applied to `P` (see
the   function   `tracing`   below),   whereas   a  call  of  the  function
`add_generator` (which does not perform Tietze transformations) will delete
these lists and hence terminate the tracing.

Example. (Presentation(PerfectGroup(960,1)) in GAP3).
```julia-rep1
P=Presentation("
1: a=a⁻
2: e=e⁻
3: f=f⁻
4: c=c⁻
5: d=d⁻
6: b³=1
7: fc=c⁻f⁻
8: dc=c⁻d⁻
9: fe=e⁻f⁻
10: ec=c⁻e⁻
11: da=a⁻f⁻
12: ed=d⁻e⁻
13: fd=d⁻f⁻
14: ca=a⁻e⁻
15: ea=a⁻c⁻
16: fa=a⁻d⁻
17: fb=be⁻
18: b⁻cbfd=1
19: b⁻ebfe=1
20: b⁻dbfedc=1
21: babab=a⁻b⁻a⁻b⁻a⁻
")
Presentation: 6 generators, 21 relators of total length 84

julia> gens(P)
6-element Vector{AbsWord}:
 a
 b
 c
 d
 e
 f
julia> simplify(P)
Presentation: 3 generators, 10 relators of total length 81
Presentation: 3 generators, 10 relators of total length 80

julia> showgens(P)
1. a 31 occurrences (2 times in relator 1) involution
2. b 26 occurrences (3 times in relator 3)
3. d 23 occurrences (2 times in relator 2) involution

julia> a,b=P.generators[1:2]
2-element Vector{AbsWord}:
 a
 b

julia> GroupPresentations.Substitute( P, a*b, :c )
#Substitute new generator c defined by ab
Presentation: 4 generators, 11 relators of total length 83

julia> simplify(P)
Presentation: 3 generators, 10 relators of total length 74
Presentation: 3 generators, 9 relators of total length 58

julia> showgens(P)
1. a 17 occurrences (2 times in relator 1) involution
2. d 19 occurrences (2 times in relator 2) involution
3. c 22 occurrences (3 times in relator 4)
```
"""
function Substitute(P::Presentation,word::Union{Vector{Int},AbsWord},s::Symbol=
   Symbol(""))
  Check(P)
  if word isa Vector
    tzword=reduceword(word)
    word=AbsWord(tzword,P)
  else
    tzword=TietzeWord(word,P.generators)
  end
  images=0
  if haskey(P,:imagesOldGens)
    images=P.imagesOldGens
    delete!(P,:imagesOldGens)
  end
  if s==Symbol("")
    add_generator(P)
    gen=gens(P)[end]
  else
    gen=AbsWord(s)
    add_generator(P, gen)
  end
  debug(P,1,"Substitute new generator ", gen," defined by ", word)
  add_relator(P,gen\word)
  if images isa Vector
    P.imagesOldGens=images
    UpdateGeneratorImages(P, tzword)
  end
  PrintStatus(P,1)
end

"""
`tracing(P::Presentation,nostop=true)`

a  sequence of  Tietze transformations  applied to  a presentation `P₁` and
ending  up with a presentation `P₂`, defines an isomorphism `φ` between the
groups defined by `P₁` and `P₂`, respectively. To know `φ` (resp. `φ⁻¹`) we
need  to know `φ(old generators)` (resp. `φ⁻¹(new generators)`). The Tietze
transformations  functions can  trace these.  This is  not done  by default
since  the involved words may grow to tremendous length; it will be done on
request  by  calling  `tracing(P)`.  A  call  `tracing(P,false)`  will stop
tracing.

`tracing` initializes three fields of `P`:

`P.oldGenerators`:  is  initialized  to  a  copy of `P.generators`.

`P.imagesOldGens`:  is the list `φ(old generators)`  as Tietze words in the
new generators. Its `i`-th entry is initialized to the Tietze word `[i]`.

`P.preImagesNewGens`: The list `φ⁻¹(new generators)` as Tietze words in the
old generators. Its `i`-th entry is initialized by the Tietze word `[i]`.

The  existence of  these fields  will cause  the Tietze  transformations to
update the lists `imagesOldGens` and `preImagesNewGens`.

There are a few restrictions concerning the tracing of generator images:

The  functions `add_generator`, `add_relator`, and `remove_relator` may change
the  isomorphism type  of the  presentation. Therefore,  if any  of them is
called, it will call `tracing(P,false)`.

You  can reinitialize  tracing at  any later  state by  calling `tracing(P)`
again: if the above fields do already exist when `tracing` is being called,
they will be initialized again.

See an example with `images`.
"""
function tracing(P::Presentation,nostop=true)
  if !nostop
    if !haskey(P, :imagesOldGens) return end
    delete!(P, :imagesOldGens)
    delete!(P, :preImagesNewGens)
    delete!(P, :oldGenerators)
    debug(P,1,"terminated tracing generator images")
  else
    P.oldGenerators=deepcopy(P.generators)
    P.imagesOldGens=map(i->[i],eachindex(P.generators))
    P.preImagesNewGens=map(i->[i],eachindex(P.generators))
    debug(P,1,"started tracing generator images")
  end
end

"""
`UpdateGeneratorImages(T::Presentation, word )`

assumes  that a new  generator defined by  the Tietze word  `word` has just
been  added to the presentation.  It converts `word` from  a Tietze word in
the  new generators to  a Tietze word  in the old  generators and adds that
word to the list of preimages.
"""
function UpdateGeneratorImages(T::Presentation, word)
  if !haskey(T,:preImagesNewGens) return end
  preim=T.preImagesNewGens
  newim=Int[]
  for num in word append!(newim, num>0 ? preim[num] : -reverse(preim[-num])) end
  push!(preim, reduceword(newim))
end

"""
`UpdateGeneratorImages(T::Presentation, n::Integer, word )`

assumes  that  the  `n`-th  generator  has  just  been replaced by the word
`word`.  It  updates  the  images  of  the old generators by replacing each
occurrence of the `n`-th generator by the given Tietze word `word`. This is
all if `word` does not involve `n` itself. If `word` involves `n` once then
the  preimages are updated  in turn. If  `word` involves `n`  more than one
time then the preimages are no more maintained.
"""
function UpdateGeneratorImages(T::Presentation, n::Integer, word)
  if n<0 error("i did not know it was possible") end
  invert(w)=-reverse(w)
  if !haskey(T,:imagesOldGens) return end
  invword=invert(word)
  for (i,image) in enumerate(T.imagesOldGens)
    newim=Int[]
    for (j,x) in enumerate(image)
      if x==n append!(newim, word)
      elseif x==-n append!(newim, invword)
      else push!(newim, x)
      end
    end
    T.imagesOldGens[i]=reduceword(newim)
  end
  if !haskey(T,:preImagesNewGens) return end
  p=findall(x->x==n || x==-n, word)
  if length(p)!=1
   if length(p)!=0 delete!(T,preImagesNewGens) end
    return
  end
  p=only(p)
  n1=word[p]
  word=vcat(invert(word[1:p-1]),n,invert(word[p+1:end]))
  word=vcat(map(x->x>0 ? T.preImagesNewGens[x] : invert(T.preImagesNewGens[-x]),
                word)...)
  T.preImagesNewGens[abs(n1)]=reduceword(n1>0 ? word : invert(word))
# images(T)
end

"""
`images(P::Presentation)`

If  `P` is a presentation in which generator images and preimages are being
traced  through all Tietze transformations  applied to `P`, `images` prints
the  preimages of the current generators as words in the old generators and
the images of the old generators as words in the current generators.

```julia-repl
julia> P=Presentation("
1: a=a⁻
2: e=e⁻
3: f=f⁻
4: c=c⁻
5: d=d⁻
6: b³=1
7: fc=c⁻f⁻
8: dc=c⁻d⁻
9: fe=e⁻f⁻
10: ec=c⁻e⁻
11: da=a⁻f⁻
12: ed=d⁻e⁻
13: fd=d⁻f⁻
14: ca=a⁻e⁻
15: ea=a⁻c⁻
16: fa=a⁻d⁻
17: fb=be⁻
18: b⁻cbfd=1
19: b⁻ebfe=1
20: b⁻dbfedc=1
21: babab=a⁻b⁻a⁻b⁻a⁻
")
Presentation: 6 generators, 21 relators of total length 84
julia> tracing(P)
julia> GroupPresentations.Go(P)
Presentation: 3 generators, 10 relators of total length 81

julia> P.imagesOldGens
6-element Vector{Vector{Int64}}:
 [1]
 [2]
 [1, -2, 1, 3, 1, 2, 1]
 [3]
 [-2, 1, 3, 1, 2]
 [1, 3, 1]

julia> P.preImagesNewGens
3-element Vector{Vector{Int64}}:
 [1]
 [2]
 [4]
```

```julia-rep1
julia> GroupPresentations.images(P)
current generators in terms of the old ones:
  a=a
  b=b
  d=d
old generators in terms of the current ones:
  a=a
  b=b
  c=ab⁻¹adaba
  d=d
  e=b⁻¹adab
  f=ada
```
"""
function images(P::Presentation)
  Check(P)
  if !haskey(P, :imagesOldGens)
    println("#I presentation is not tracing generator images"); return
  end
  if haskey(P,:preImagesNewGens)
    println("current generators in terms of the old ones:")
    for i in eachindex(P.generators) println(rio(),"  ",P.generators[i],"=",
           AbsWord(P.preImagesNewGens[i],P.oldGenerators)) end
  end
  println("old generators in terms of the current ones:")
  for i in eachindex(P.oldGenerators) println(rio(),"  ",P.oldGenerators[i],"=",
         AbsWord(P.imagesOldGens[i],P))
  end
end

"""
RenumberGenerators( presentation, sort list )

renumbers the generators of the given presentation
according to the given sort list.
"""
function RenumberGenerators(P::Presentation, L::AbstractVector)
  debug(P,3,"RenumberGenerators")
  RemoveGenerators(P)
  numgens=length(P.generators)
  if sort(L)!=1:numgens
    error("<L> must determine a permutation of the generator numbers")
  end
  invsort=invperm(L)
  gens=P.generators
  rels=P.relators
  oldgens=deepcopy(gens)
  oldinvs=copy(P.inverses)
  new=copy(P.inverses)
  println("gens==",gens,"\ninvs==",P.inverses,"\nrels==",rels)
  for i in 1:numgens
    j=L[i]
    gens[i]=oldgens[j]
    if oldinvs[numgens+1+j]==j P[-i]=i else P[-i]=-i end
    new[numgens+1+i]=invsort[i]
    new[numgens+1-i]=-invsort[i]
  end
  println("gens==",gens,"\ninvs==",P.inverses,"\nnew==",new)
  for rel in P.relators
    for i in eachindex(rel) rel[i]=new[numgens+1+rel[i]] end
  end
  println("rels==", rels)
  if haskey(P, :imagesOldGens)
    images(P)
    for i in 1:length(P.imagesOldGens)
      P.imagesOldGens[i]=map(j->new[numgens+1+j],P.imagesOldGens[i])
    end
    permute!(P.preImagesNewGens, L)
    println("renumbered to")
    images(P)
  end
end

"""
`simplify(p [,tries])`

simplify  the  presentation  `p`.  We  have  found heuristics which make it
somewhat  efficient, but the algorithm depends  on random numbers so is not
reproducible.  The main  idea is  to rotate  relators between  calls to the
basic  `GroupPresentations.Go` function. By default 100 such rotations are tried
(unless  the  presentation  is  so  small  that  less rotations exhaust all
possible  ones), but the actual number tried  can be controlled by giving a
second  parameter `tries` to the function. Another useful tool to deal with
presentations is `tryconjugate`.

```julia-rep1
julia> P=Presentation("
1: ab=ba
2: dbd=bdb
3: bcb=cbc
4: cac=aca
5: adca=cadc
6: dcdc=cdcd
7: adad=dada
8: d⁻bdcbd=cd⁻bdcb
9: adcd⁻ad=dcd⁻adc
10: dcdadc=adcdad
11: dcabdcbda=adbcbadcb
12: cac⁻bdcbad=bdcbadb⁻cb
13: cbd⁻adcbad=bd⁻adcbadc
14: cda⁻bc⁻adb⁻c=bdca⁻bc⁻db⁻a
15: cdc⁻bdcabdc=bdcbadcdad⁻
16: d⁻²b⁻c³bdca⁻b=ca⁻bc⁻dcb⁻c⁻d²c
17: c⁻dab⁻dba⁻dcbc⁻ad=abdca⁻bd⁻adb⁻c⁻b²
18: bdbcabdca⁻²d⁻a⁻db⁻d⁻a=cbadcbd⁻adcb⁻d⁻a⁻b⁻d⁻b
19: c⁻bdbadcd⁻b²dc⁻bd⁻²adcb⁻c⁻d⁻a⁻db⁻c⁻d⁻bdad⁻c⁻d⁻bdcbadcb⁻c⁻d⁻a⁻db⁻c⁻d⁻b⁻²dacd⁻bdc²b=abdbcabdca⁻dcbc⁻d⁻²b⁻c⁻d⁻a⁻b⁻d⁻a⁻b⁻d⁻b²dcbd⁻adcbc⁻d⁻a⁻db⁻c⁻abd⁻a⁻c⁻bdb⁻adcad⁻ba⁻d²
")
Presentation: 4 generators, 19 relators of total length 370
julia> simplify(p)
Presentation: 4 generators, 18 relators of total length 304
Presentation: 4 generators, 18 relators of total length 284
Presentation: 4 generators, 17 relators of total length 264
Presentation: 4 generators, 16 relators of total length 256
Presentation: 4 generators, 15 relators of total length 244
Presentation: 4 generators, 15 relators of total length 240
Presentation: 4 generators, 15 relators of total length 226
Presentation: 4 generators, 15 relators of total length 196
Presentation: 4 generators, 15 relators of total length 178
Presentation: 4 generators, 15 relators of total length 172
Presentation: 4 generators, 14 relators of total length 158

julia> display_balanced(p)
1: ab=ba
2: dbd=bdb
3: bcb=cbc
4: cac=aca
5: ada⁻c=cada⁻
6: dcdc=cdcd
7: adad=dada
8: c⁻db⁻cbd=bc⁻db⁻cb
9: adcd⁻ad=dcd⁻adc
10: dcdadc=adcdad
11: cbdcbdc=dcbdcbd
12: dcbadcbda=adcbcadcb
13: cbc⁻d⁻adcab=d⁻adcbadcd⁻
14: cad⁻c⁻bdb⁻ca⁻d⁻bda=bd⁻b⁻ad⁻ba⁻d⁻cbadc
```
"""
function simplify(P::Presentation,lim=100)
  if isempty(P.relators) return end
  rot(tt,i)=tt[i]=circshift(tt[i],-1)
  function test()
    if prod(filter(!iszero,big.(length.(P.relators))))<lim
      tt=P.relators
      v=fill(0,length(tt))
      if length(v)==0 return false end
      while true
        before=[length(tt),sum(length,tt)]
        j=length(v)
        while v[j]==length(tt[j])-1
          rot(tt,j)
          v[j]=0
          j-=1
          if j==0 return false end
        end
        rot(tt,j)
        v[j]+=1
        GoGo(P)
        tt=P.relators
        if [length(tt),sum(length,tt)]<before return true end
      end
    else
      for k in 1:lim
        tt=P.relators
        before=[length(tt),sum(length,tt)]
        rot(tt,rand(eachindex(tt)))
        GoGo(P)
        tt=P.relators
        if [length(tt),sum(length,tt;init=0)]<before return true end
      end
    end
    return false
  end
  while test() end
end

"""
`conjugate(p,conjugation)`

This program modifies a presentation by conjugating a generator by another.
The  conjugation to  apply is  described by  a length-3  string of the same
style  as the result of `display_balanced`,  that is `"aba⁻"` means replace
the  second generator  by its  conjugate by  the first, and  `"a⁻ba"` means
replace it by its conjugate by the inverse of the first.

```julia-rep1
julia> P=Presentation("
1: dabcd=abcda
2: dabcdb=cabcda
3: bcdabcd=dabcdbc
")
Presentation: 4 generators, 3 relators of total length 36

julia> display_balanced(conjugate(P,"c⁻dc"))
<< presentation with 4 generators, 3 relators of total length 36>>
1: dcabdc=cabdca
2: abdcab=cabdca
3: bdcabd=cabdca
```
"""
function conjugate(p::Presentation, s::Vector{<:Integer})
  if length(s)!=3 || s[1]!=-s[3] error("should be conjugate like aba⁻") end
  p=deepcopy(p);SubstituteGen(p,s[2],s)
  simplify(p, 100)
  p
end

conjugate(p::Presentation, s::String)=conjugate(p,TietzeWord(AbsWord(s),gens(p)))

"""
`tryconjugate(p;goal=[0,0],info=[0,0])`

This program tries to simplify group presentations by applying conjugations
to  the  generators.  The  algorithm  depends  on  random  numbers,  and on
tree-searching,  so is  not reproducible.  By default  the program stops as
soon  as a shorter presentation is found.  Sometimes this does not give the
desired  presentation.  One  can  give  a  second argument `goal`, then the
program  will only stop when  a presentation of length  less than `goal` is
found.  Finally, a third  argument can be  given and then all presentations
the  programs runs  over which  are of  length less  than or  equal to this
argument are displayed. Due to the non-deterministic nature of the program,
it  may be useful to  run it several times  on the same input. Upon failure
(to improve the presentation), the program returns `p`.

```julia-rep1
julia> P=Presentation("
1: ba=ab
2: dbd=bdb
3: cac=aca
4: bcb=cbc
5: da⁻ca=a⁻cad
6: dcdc=cdcd
7: adad=dada
8: dcd⁻bdc=bdcbdb⁻
9: dcdadc=adcdad
10: adcd⁻ad=dcd⁻adc
11: b⁻c³bdca⁻b=dcba⁻c⁻d³c
")
Presentation: 4 generators, 11 relators of total length 100

julia> p=tryconjugate(p)
Presentation: 4 generators, 11 relators of total length 100
dcd⁻=> Presentation: 4 generators, 10 relators of total length 90
# dcd⁻ gives Presentation: 4 generators, 10 relators of total length 90
Presentation: 4 generators, 10 relators of total length 90

julia> p=tryconjugate(p)
d⁻cd=> Presentation: 4 generators, 10 relators of total length 88
# d⁻cd gives Presentation: 4 generators, 10 relators of total length 88
Presentation: 4 generators, 10 relators of total length 88

julia> p=tryconjugate(p)
dcd⁻=> Presentation: 4 generators, 10 relators of total length 90
d⁻bd=> Presentation: 4 generators, 10 relators of total length 96
a⁻ca=> Presentation: 4 generators, 9 relators of total length 84
Presentation: 4 generators, 8 relators of total length 76
# a⁻ca gives Presentation: 4 generators, 8 relators of total length 76
Presentation: 4 generators, 8 relators of total length 76

julia> p=tryconjugate(p)
b⁻cb=> Presentation: 4 generators, 8 relators of total length 70
# b⁻cb gives Presentation: 4 generators, 8 relators of total length 70
Presentation: 4 generators, 8 relators of total length 70

julia> p=tryconjugate(p)
c⁻ac=> Presentation: 4 generators, 8 relators of total length 64
# c⁻ac gives Presentation: 4 generators, 8 relators of total length 64
Presentation: 4 generators, 8 relators of total length 64

julia> p=tryconjugate(p)
cac⁻=> Presentation: 4 generators, 8 relators of total length 58
# cac⁻ gives Presentation: 4 generators, 8 relators of total length 58
Presentation: 4 generators, 8 relators of total length 58

julia> p=tryconjugate(p)
c⁻ac=> Presentation: 4 generators, 8 relators of total length 64
c⁻bc=> Presentation: 4 generators, 7 relators of total length 50
# c⁻bc gives Presentation: 4 generators, 7 relators of total length 50
Presentation: 4 generators, 7 relators of total length 50

julia> p=tryconjugate(p)
cdc⁻=> Presentation: 4 generators, 7 relators of total length 56
d⁻cd=> Presentation: 4 generators, 7 relators of total length 54
c⁻ac=> Presentation: 4 generators, 7 relators of total length 48
# c⁻ac gives Presentation: 4 generators, 7 relators of total length 48
Presentation: 4 generators, 7 relators of total length 48

julia> p=tryconjugate(p)
cac⁻=> Presentation: 4 generators, 7 relators of total length 50
c⁻dc=> Presentation: 4 generators, 7 relators of total length 50
d⁻bd=> d⁻cd=> Presentation: 4 generators, 7 relators of total length 60
b⁻ab=> a⁻ba=> a⁻ca=> Presentation: 4 generators, 7 relators of total length 46
# a⁻ca gives Presentation: 4 generators, 7 relators of total length 46
Presentation: 4 generators, 7 relators of total length 46

julia> display_balanced(p)
1: da=ad
2: ca=ac
3: aba=bab
4: bcb=cbc
5: cdc=dcd
6: dbdb=bdbd
7: bdcbdc=dcbdcb
```
"""
function tryconjugate(p::Presentation;goal=[0,0],info=[0,0])
  perf(p)=[length(p.relators), totallength(p)]
  function triples(p)
    res=Vector{Int}[]
    for v in p.relators, i in 1:length(v)-2
      if v[i]==-v[i+2]
        if v[i+1]>0 push!(res,v[i:i+2])
        else        push!(res,-reverse(v[i:i+2]))
        end
      end
    end
    first.(sort(tally(res),by=x->-x[2]))
  end
  function expand(p)
    res=map(x->[x, conjugate(p, x)], triples(p))
    sort!(res,by=x->perf(x[2]))
  end
  p=deepcopy(p)
  GoGo(p)
  if iszero(goal) goal=perf(p) end
  p1=[]
  for c in triples(p)
    print(rio(),AbsWord(c,p),"=> ")
    n=deepcopy(p);SubstituteGen(n,c[2],[-c[1],c[2],-c[3]]);GoGo(n)
    if perf(n)<=info
      println(rio(),"# ",AbsWord(c,p)," gives rels/len ", perf(n))
      display_balanced(n)
    end
    if perf(n)<goal
      println(rio(),"# ",AbsWord(c,p)," gives ",n)
      return n
    end
    push!(p1,(c,n))
  end
  sort!(p1,by=x->perf(x[2]))
  for p2 in p1, c in triples(p2[2])
    n=deepcopy(p2[2]);SubstituteGen(n,c[2],[-c[1],c[2],-c[3]]);GoGo(n)
    print(rio(),AbsWord(p2[1],p),"+",AbsWord(c,p),"=>")
    if perf(n)<=info
      println(rio(),"# ",AbsWord(p2[1],p),"->",AbsWord(c,p)," gives rels/len ", perf(n))
      display_balanced(n)
    end
    if perf(n)<goal
      println(rio(),"# ",AbsWord(p2[1],p),"->",AbsWord(c,p)," gives ",n)
      return n
    end
  end
  println("\n# could not shrink ",p)
  p
end

p1=Presentation(
[[1, 1], [7, 7], [3, 3], [6, 6], [8, 8], [4, 4], [6, 2, 3], [4, -2, 1], [7, 4,
 2], [3, 5, 2], [5, 1, 2], [3, 1, -2], [-2, -2, -2], [3, -2, 4], [1, 3, 2],
 [1, 6, -2], [2, 2, 2], [2, 1, 3], [3, 2, 1], [1, -2, 3], [8, 1, -2], [4, 1,
 2], [-2, 5, 3], [-2, 1, 4], [-2, 3, 1], [1, 7, 2], [1, 2, 1, 2], [1, 2, 3,
 -2], [1, 3, -5, 3], [1, 3, -2, -2], [3, -2, 3, -2], [1, 3, 1, 3, 1, 3]])
p2=Presentation(
[[1, 1], [5, 5], [6, 6], [3, 3], [4, 4], [2, 2, 2], [6, 3, 6, 3], [4, 3, 4,
 3], [6, 5, 6, 5], [5, 3, 5, 3], [4, 1, 6, 1], [5, 4, 5, 4], [6, 4, 6, 4], [3,
 1, 5, 1], [5, 1, 3, 1], [6, 1, 4, 1], [6, 2, 5, -2], [-2, 3, 2, 6, 4], [-2,
 5, 2, 6, 5], [-2, 4, 2, 6, 5, 4, 3], [2, 1, 2, 1, 2, 1, 2, 1, 2, 1]])
p3=Presentation(
[[1, 1], [18, 18], [2, 2], [7, 7], [16, 16], [17, 17], [12, 12], [10, 10], [9,
 5, 4], [3, 3, 3], [4, 3, 1], [9, 9, 9], [5, 5, 5], [20, 1, 15], [6, 1, 3],
 [23, 23, 23], [22, 14, 12], [13, -8, 10], [21, 20, 6], [23, 21, 16], [-13,
 11, 7], [19, 14, 17], [20, -19, 6], [22, 22, 22], [7, 3, 14, -13], [18, 7,
 19, -11], [12, 15, 16, -15], [7, 2, 7, 17], [-8, 6, -9, -6], [2, 8, 15, -4],
 [7, 10, 7, 12], [12, 2, -6, 22], [7, 18, 7, 18], [7, -3, 2, 5, 2], [23, 5, 5,
 4, -21], [18, 8, 14, 14, 7], [21, 21, 21, 21, 21], [11, 15, 16, 2, 1, -8],
 [7, 1, 12, 2, -5, -6], [7, 19, 6, 2, -6, -19], [1, 7, 6, 5, 2, 12], [1, 8,
 14, 7, 2, -9], [8, 1, 2, 16, -15, -11], [23, -21, -4, 2, 7, 2, 6], [10, -3,
 1, 2, 5, -4, 2], [23, 21, -15, 14, 6, -4, -21], [11, -8, 2, 4, 9, 9, 2], [10,
 2, 4, -5, 2, 1, 3], [7, 2, 3, 1, -3, 1, 2, 8], [17, 12, 1, 6, 16, 21, 9, -6],
 [2, 1, 2, 3, 1, 2, -3, 5], [1, 2, 1, 3, 2, 1, 2, -3], [-11, 7, -8, 13, -11,
 7, -8, 13], [18, 2, 5, 4, 20, 6, -4, -6], [11, -13, 8, 1, 6, -4, 2, 18], [-6,
 1, 6, 16, -6, 1, 6, 16], [18, 12, 14, 19, -8, 18, 11, -8], [7, 18, 7, 19, 12,
 17, 12, -19], [-6, -8, 2, -5, 4, -6, 1, 3, 2], [19, 6, -15, 6, -9, -6, 1, -8,
 7], [22, 12, 7, 2, 1, 3, 3, 19, 12], [22, 2, -9, 5, 5, 2, 18, 7, 19], [2, 4,
 -9, -6, 1, 6, -21, 23, -6, -11], [2, 7, 2, 4, 9, 9, 1, -9, -9, -4], [16, -21,
 20, -8, 13, 7, -13, 8, -20, 21], [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2], [12,
 7, 18, 11, 7, 12, -14, -19, 7, 18, 7, 22], [1, -20, -6, 1, 6, 20, 1, -20, -6,
 1, 6, 20], [12, 14, 17, -14, 12, 14, 17, -14, 12, 14, 17, -14], [7, 17, 7, 2,
 4, 16, -4, 2, 7, 17, 7, 18], [10, 18, 7, 12, 17, 12, 7, 18, 10, 13, 7, -13],
 [2, 8, 7, -13, 10, 13, 7, -8, 2, 4, -6, 12, 6, -4], [11, -13, 10, 13, -11,
 18, 11, -13, 10, 13, -11, 18, 11, -13, 10, 13, -11, 18], [10, 13, 7, -13, 10,
 13, 7, -13, 10, 13, 7, -13, 10, 13, 7, -13, 10, 13, 7, -13]])
p3simp=Presentation(
[[2, 2, 2], [2, -1, 2, -1, 2, -1], [1, 1, -2, -1, -2, 1, 1, -2, -1, -2], [1,
 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], [1, 1, -2, -1, -1, -2, 1, 1, 2, -1, -1, -1,
 2], [1, 1, -2, -1, -1, -2, -1, 2, 1, -2, -1, 2, 1]])
p4=Presentation(
[[2, 3, -2, -3], [5, 3, 5, -3, -5, -3], [3, 4, 3, -4, -3, -4], [4, 2, 4, -2,
 -4, -2], [2, 5, 4, 2, -4, -5, -2, -4], [5, 4, 5, 4, -5, -4, -5, -4], [2, 5,
 2, 5, -2, -5, -2, -5], [5, 4, 5, 2, 5, 4, 1, -5, -2, -5, -4, -5, -2], [-5, 3,
 5, 4, 3, 5, 1, 1, -3, -4, -5, -3, 5, -4], [2, 5, 4, -5, 2, 5, 1, 1, -4, -5,
 -2, 5, -4, -5], [5, 4, 2, 3, 5, 4, 3, 5, 2, -3, -4, -5, -2, -3, -4, -3, -5,
 -2], [4, 2, -4, 3, 5, 4, 3, 2, 5, -3, -4, 3, -5, -2, -3, -4, -5, -3], [4, 3,
 -5, 2, 5, 4, 3, 2, 5, -4, -5, -2, -3, -4, -5, -2, 5, -3], [4, 5, -2, 3, -4,
 2, 5, -3, 4, -2, 3, -5, 4, -3, 2, -4, -5, -3], [4, 5, -4, 3, 5, 4, 2, 3, 5,
 4, 5, -2, -5, -4, -5, -2, -3, -4, -5, -3], [-5, -5, -3, 4, 4, 4, 3, 5, 4, -2,
 3, -4, -5, -5, 4, 3, -4, -5, 4, -3, 2, -4], [-4, 5, 2, -3, 5, 3, -2, 5, 4, 3,
 -4, 2, 5, -3, -3, 4, 3, -5, -2, 5, -3, 2, -4, -5, -3, -2], [3, 5, 3, 4, 2, 3,
 5, 4, -2, -2, -5, -2, 5, -3, -5, 2, -3, 5, 3, 2, 5, 3, -4, -5, -2, 5, -3, -4,
 -5, -2, -3, -4], [-4, 3, 5, 3, 2, 5, 4, -5, 3, 3, 5, -4, 3, -5, -5, 2, 5, 4,
 -3, -4, -5, -2, 5, -3, -4, -5, 3, 5, 2, -5, -4, -5, 3, 5, 4, 3, 2, 5, 4, -3,
 -4, -5, -2, 5, -3, -4, -5, -3, -3, 5, 2, 4, -5, 3, 5, 4, 4, 3, -5, -5, 2, -3,
 5, -2, -4, -5, -2, 3, -5, -3, 4, 2, 5, -3, -2, 4, 3, -5, 2, 5, 4, -3, -4, -5,
 -2, 5, -3, -4, -5, -3, -3, 5, 3, 2, 5, 3, 2, 5, 4, 3, 5, 5, 4, -3, -4, -5, 2,
 -4, -5, -3, -2, -4, -3, -5, -3, -2]])
p5=Presentation(
[[3, -1], [2, -3], [2, 3, -1, -3], [1, 4, -1, -4], [4, 1, -4, -1], [1, 3, -2,
 -3], [1, 2, -3, -2], [4, 3, -4, -3], [3, 4, -3, -4], [2, 3, 4, 1], [3, 2, 3,
 -2, -3, -2], [-1, 4, 1, 2, -1, -4, 1, 4, -2, -4], [-1, 4, 1, 3, -1, -4, 1, 4,
 -3, -4]])
p6=Presentation(
[[1, 2, -1, -2], [4, 2, 4, -2, -4, -2], [2, 3, 2, -3, -2, -3], [3, 1, 3, -1,
 -3, -1], [1, 4, 3, 1, -3, -4, -1, -3], [4, 3, 4, 3, -4, -3, -4, -3], [1, 4,
 1, 4, -1, -4, -1, -4], [-4, 2, 4, 3, 2, 4, -2, -3, -4, -2, 4, -3], [1, 4, 3,
 -4, 1, 4, -3, -4, -1, 4, -3, -4], [4, 3, 4, 1, 4, 3, -4, -1, -4, -3, -4, -1],
 [4, 3, 1, 2, 4, 3, 2, 4, 1, -2, -3, -4, -1, -2, -3, -2, -4, -1], [3, 1, -3,
 2, 4, 3, 2, 1, 4, -2, -3, 2, -4, -1, -2, -3, -4, -2], [3, 2, -4, 1, 4, 3, 2,
 1, 4, -3, -4, -1, -2, -3, -4, -1, 4, -2], [3, 4, -1, 2, -3, 1, 4, -2, 3, -1,
 2, -4, 3, -2, 1, -3, -4, -2], [3, 4, -3, 2, 4, 3, 1, 2, 4, 3, 4, -1, -4, -3,
 -4, -1, -2, -3, -4, -2], [-4, -4, -2, 3, 3, 3, 2, 4, 3, -1, 2, -3, -4, -4, 3,
 2, -3, -4, 3, -2, 1, -3], [-3, 4, 1, -2, 4, 2, -1, 4, 3, 2, -3, 1, 4, -2, -2,
 3, 2, -4, -1, 4, -2, 1, -3, -4, -2, -1], [2, 4, 2, 3, 1, 2, 4, 3, -1, -1, -4,
 -1, 4, -2, -4, 1, -2, 4, 2, 1, 4, 2, -3, -4, -1, 4, -2, -3, -4, -1, -2, -3],
 [-3, 2, 4, 2, 1, 4, 3, -4, 2, 2, 4, -3, 2, -4, -4, 1, 4, 3, -2, -3, -4, -1,
 4, -2, -3, -4, 2, 4, 1, -4, -3, -4, 2, 4, 3, 2, 1, 4, 3, -2, -3, -4, -1, 4,
 -2, -3, -4, -2, -2, 4, 1, 3, -4, 2, 4, 3, 3, 2, -4, -4, 1, -2, 4, -1, -3, -4,
 -1, 2, -4, -2, 3, 1, 4, -2, -1, 3, 2, -4, 1, 4, 3, -2, -3, -4, -1, 4, -2, -3,
 -4, -2, -2, 4, 2, 1, 4, 2, 1, 4, 3, 2, 4, 4, 3, -2, -3, -4, 1, -3, -4, -2,
 -1, -3, -2, -4, -2, -1]])
end
