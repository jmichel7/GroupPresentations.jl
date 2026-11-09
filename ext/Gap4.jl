module Gap4

using GAP, GroupPresentations, PermGroups

function GAP.Obj(g::FpGroup)
  get!(g,:GAP)do
    v=Symbol.(g.gens)
    gapg=GAP.Globals.FreeGroup(GAP.Obj(v,recursive=true))
    gapgens=getproperty.(Ref(gapg),eachindex(v))
    gaprels=map(g.rels)do rel
      res=gapgens[1]^0
      for p in rel.d
        res*=gapgens[findfirst(==(p[1]),v)]^p[2]
      end
      res
    end
    gapg/GAP.Obj(gaprels)
  end
end

function Base.length(g::FpGroup)
  get!(g,:length)do
    GAP.Globals.Size(GAP.Obj(g))
  end
end

function Groups.elements(g::FpGroup)
  get!(g,:elements)do
    l=collect(GAP.Globals.Elements(GAP.Obj(g)))
    l=collect.(GAP.Globals.ExtRepOfObj.(l))
    map(w->prod(x->g.gens[x[1]]^x[2],Iterators.partition(w,2);init=one(g)),l)
  end
end

end
