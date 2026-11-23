# auto-generated tests from julia-repl docstrings
using Test, GroupPresentations, PermGroups
function mytest(file::String,cmd::String,man::String)
  println(file," ",cmd)
  exec=repr(MIME("text/plain"),eval(Meta.parse(cmd)),context=:limit=>true)
  if endswith(cmd,";") return true end
  exec=replace(exec,r"\s*$"m=>"",r"\s*$"s=>"",r"^\s*"=>"")
  if exec==man return true end
  inds=collect(eachindex(exec))
  i=findfirst(i->i<=lastindex(man) && exec[i]!=man[i],inds)
  if i==nothing i=ncodeunits(man)+1 else i=inds[i] end
  print("exec=$(repr(exec[i:end]))\nmanl=$(repr(man[i:end]))\n")
  false
end
@testset verbose=true "GroupPresentations.jl" begin
@test mytest("GroupPresentations.jl","@AbsWord a,b,c,d,e,f","nothing")
@test mytest("GroupPresentations.jl","F=FpGroup([a,b,c,d,e,f])","FreeGroup(a,b,c,d,e,f)")
@test mytest("GroupPresentations.jl","G=F/[a^2,b^2,d*f^-1,e^2,f^2,a*b^-1*c,a*e*c^-1,b*d^-1*c,c*d*e^-1,a*f*c^-2,c^4]","FreeGroup(a,b,c,d,e,f)/[a²,b²,df⁻,e²,f²,ab⁻c,aec⁻,bd⁻c,cde⁻,afc⁻²,c⁴]")
@test mytest("GroupPresentations.jl","@AbsWord a,b","nothing")
@test mytest("GroupPresentations.jl","F=FpGroup([a,b])","FreeGroup(a,b)")
@test mytest("GroupPresentations.jl","G=F/[a^2,b^7,comm(a,a^b),comm(a,a^(b^2))*inv(b^a)]","FreeGroup(a,b)/[a²,b⁷,a⁻b⁻a⁻bab⁻ab,a⁻b⁻²a⁻b²ab⁻²ab²a⁻b⁻a]")
@test mytest("GroupPresentations.jl","P=Presentation(G)","Presentation: 2 generators, 4 relators of total length 30")
@test mytest("GroupPresentations.jl","relators(P)","4-element Vector{AbsWord}:\n a²\n b⁷\n ab⁻abab⁻ab\n b⁻²ab²ab⁻²ab²ab⁻")
@test mytest("GroupPresentations.jl","w=AbsWord([:a => 3, :b => -2, :a => 1])","a³b⁻²a")
@test mytest("GroupPresentations.jl","w[4]",":b => -1")
@test mytest("GroupPresentations.jl","w[3:4]","ab⁻")
@test mytest("GroupPresentations.jl","length(w)","6")
@test mytest("GroupPresentations.jl","w[begin:end]","a³b⁻²a")
@test mytest("GroupPresentations.jl","w*AbsWord([:a=>-1,:b=>1])","a³b⁻")
@test mytest("GroupPresentations.jl","AbsWord(:b,:a,:a,:b)","ba²b")
@test mytest("GroupPresentations.jl","AbsWord(\"aec⁻\")","aec⁻")
@test mytest("GroupPresentations.jl","AbsWord(\"c²=c⁻²\")","c⁴")
@test mytest("GroupPresentations.jl","AbsWord([-1,-2,1,-3,-3],AbsWord.([:a,:b,:c]))","a⁻b⁻ac⁻²")
@test mytest("GroupPresentations.jl","AbsWord([-1,-2,1,-3,-3],[:a,:b,:c])","a⁻b⁻ac⁻²")
@test mytest("GroupPresentations.jl","F=FpGroup(:a,:b,:c)","FreeGroup(a,b,c)")
@test mytest("GroupPresentations.jl","GroupPresentations.TietzeWord(comm(F(1),F(2))*inv(F(3)^2*F(2)),gens(F))","5-element Vector{Int64}:\n -1\n -2\n  1\n -3\n -3")
@test mytest("GroupPresentations.jl","@AbsWord a,b,c,d,e,f","nothing")
@test mytest("GroupPresentations.jl","F=FpGroup([a,b,c,d,e,f])","FreeGroup(a,b,c,d,e,f)")
@test mytest("GroupPresentations.jl","G=F/[a^2,b^2,d*f^-1,e^2,f^2,a*b^-1*c,a*e*c^-1,b*d^-1*c,c*d*e^-1,a*f*c^-2,c^4]","FreeGroup(a,b,c,d,e,f)/[a²,b²,df⁻,e²,f²,ab⁻c,aec⁻,bd⁻c,cde⁻,afc⁻²,c⁴]")
@test mytest("GroupPresentations.jl","simplify(G)","FreeGroup(a,c)/[a²,ac⁻ac⁻,c⁴]")
@test mytest("GroupPresentations.jl","P=Presentation(\"\n1: a=a⁻\n2: e=e⁻\n3: f=f⁻\n4: c=c⁻\n5: d=d⁻\n6: b³=1\n7: fc=c⁻f⁻\n8: dc=c⁻d⁻\n9: fe=e⁻f⁻\n10: ec=c⁻e⁻\n11: da=a⁻f⁻\n12: ed=d⁻e⁻\n13: fd=d⁻f⁻\n14: ca=a⁻e⁻\n15: ea=a⁻c⁻\n16: fa=a⁻d⁻\n17: fb=be⁻\n18: b⁻cbfd=1\n19: b⁻ebfe=1\n20: b⁻dbfedc=1\n21: babab=a⁻b⁻a⁻b⁻a⁻\n\")","Presentation: 6 generators, 21 relators of total length 84")
@test mytest("GroupPresentations.jl","tracing(P)","nothing")
@test mytest("GroupPresentations.jl","GroupPresentations.Go(P)","Presentation: 3 generators, 10 relators of total length 81")
@test mytest("GroupPresentations.jl","P.imagesOldGens","6-element Vector{Vector{Int64}}:\n [1]\n [2]\n [1, -2, 1, 3, 1, 2, 1]\n [3]\n [-2, 1, 3, 1, 2]\n [1, 3, 1]")
@test mytest("GroupPresentations.jl","P.preImagesNewGens","3-element Vector{Vector{Int64}}:\n [1]\n [2]\n [4]")
end
