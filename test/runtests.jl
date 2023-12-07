# auto-generated tests from julia-repl docstrings
using Test, GroupPresentations, PermGroups
function mytest(file::String,cmd::String,man::String)
  println(file," ",cmd)
  exec=repr(MIME("text/plain"),eval(Meta.parse(cmd)),context=:limit=>true)
  if endswith(cmd,";") exec="nothing" 
  else exec=replace(exec,r"\s*$"m=>"")
       exec=replace(exec,r"\s*$"s=>"")
  end
  if exec!=man 
    i=1
    while i<=lastindex(exec) && i<=lastindex(man) && exec[i]==man[i]
      i=nextind(exec,i)
    end
    print("exec=$(repr(exec[i:end]))\nmanl=$(repr(man[i:end]))\n")
  end
  exec==man
end
@testset verbose=true "GroupPresentations.jl" begin
@test mytest("GroupPresentations.jl","@AbsWord a,b","nothing")
@test mytest("GroupPresentations.jl","F=FpGroup([a,b])","FreeGroup(a,b)")
@test mytest("GroupPresentations.jl","G=F/[a^2,b^7,comm(a,a^b),comm(a,a^(b^2))*inv(b^a)]","FreeGroup(a,b)/[a²,b⁷,a⁻¹b⁻¹a⁻¹bab⁻¹ab,a⁻¹b⁻²a⁻¹b²ab⁻²ab²a⁻¹b⁻¹a]")
@test mytest("GroupPresentations.jl","P=Presentation(G)","Presentation: 2 generators, 4 relators, total length 30")
@test mytest("GroupPresentations.jl","relators(P)","4-element Vector{AbsWord}:\n a²\n b⁷\n ab⁻¹abab⁻¹ab\n b⁻²ab²ab⁻²ab²ab⁻¹")
@test mytest("GroupPresentations.jl","w=AbsWord([:a => 3, :b => -2, :a => 1])","a³b⁻²a")
@test mytest("GroupPresentations.jl","w*AbsWord([:a=>-1,:b=>1])","a³b⁻¹")
@test mytest("GroupPresentations.jl","AbsWord(:b,:a,:a,:b)","ba²b")
@test mytest("GroupPresentations.jl","F=FpGroup(:a,:b,:c)","FreeGroup(a,b,c)")
@test mytest("GroupPresentations.jl","GroupPresentations.TietzeWord(comm(F(1),F(2))*inv(F(3)^2*F(2)),gens(F))","5-element Vector{Int64}:\n -1\n -2\n  1\n -3\n -3")
@test mytest("GroupPresentations.jl","AbsWord([-1,-2,1,-3,-3],AbsWord.([:a,:b,:c]))","a⁻¹b⁻¹ac⁻²")
@test mytest("GroupPresentations.jl","@AbsWord a,b,c,d,e,f","nothing")
@test mytest("GroupPresentations.jl","F=FpGroup([a,b,c,d,e,f])","FreeGroup(a,b,c,d,e,f)")
@test mytest("GroupPresentations.jl","G=F/[a^2,b^2,d*f^-1,e^2,f^2,a*b^-1*c,a*e*c^-1,b*d^-1*c,c*d*e^-1,a*f*c^-2,c^4]","FreeGroup(a,b,c,d,e,f)/[a²,b²,df⁻¹,e²,f²,ab⁻¹c,aec⁻¹,bd⁻¹c,cde⁻¹,afc⁻²,c⁴]")
@test mytest("GroupPresentations.jl","simplify(G)","FreeGroup(a,c)/[a²,ac⁻¹ac⁻¹,c⁴]")
@test mytest("GroupPresentations.jl","P=Presentation(\"\n1: a=A\n2: e=E\n3: f=F\n4: c=C\n5: d=D\n6: bbb=1\n7: fc=CF\n8: dc=CD\n9: fe=EF\n10: ec=CE\n11: da=AF\n12: ed=DE\n13: fd=DF\n14: ca=AE\n15: ea=AC\n16: fa=AD\n17: fb=bE\n18: Bcbfd=1\n19: Bebfe=1\n20: Bdbfedc=1\n21: babab=ABABA\n\")","Presentation: 6 generators, 21 relators, total length 84")
@test mytest("GroupPresentations.jl","tracing(P)","nothing")
@test mytest("GroupPresentations.jl","GroupPresentations.Go(P)","Presentation: 3 generators, 10 relators, total length 81")
@test mytest("GroupPresentations.jl","P.imagesOldGens","6-element Vector{Vector{Int64}}:\n [1]\n [2]\n [1, -2, 1, 3, 1, 2, 1]\n [3]\n [-2, 1, 3, 1, 2]\n [1, 3, 1]")
@test mytest("GroupPresentations.jl","P.preImagesNewGens","3-element Vector{Vector{Int64}}:\n [1]\n [2]\n [4]")
end
