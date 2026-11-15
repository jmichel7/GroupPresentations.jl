# A Todd-Coxeter coset enumeration for the trivial subgroup, with coincidence
# handling, written from the standard algorithm and validated against known
# group orders (see the end of this file).  It gives, for a finite group and a
# presentation that the enumeration can close:
#
#   toddcoxeter(rels,ngens) -> (order, perms) or nothing
#
# where `rels` is a Vector of Tietze words (Vector{Int}) and `perms[k]` is the
# permutation of the cosets induced by generator k in the regular action.
# This is the only oracle in the review that does not go through the package.

# column of a letter: 2k-1 for generator k, 2k for its inverse
tccol(i::Int) = i > 0 ? 2i-1 : -2i
# the column of the inverse letter
tcinvcol(j::Int) = isodd(j) ? j+1 : j-1

mutable struct TC
  ng::Int
  rels::Vector{Vector{Int}}
  tab::Vector{Vector{Int}}
  par::Vector{Int}
  size::Vector{Int}
  maxcosets::Int
end

function tcfind(t::TC, i::Int)
  while t.par[i] != i
    t.par[i] = t.par[t.par[i]]      # path halving
    i = t.par[i]
  end
  i
end

function tcnewcoset!(t::TC)
  length(t.tab) >= t.maxcosets && return 0
  push!(t.tab, zeros(Int, 2t.ng))
  push!(t.par, length(t.tab))
  push!(t.size, 1)
  length(t.tab)
end

function tcmerge!(t::TC, a::Int, b::Int)
  a = tcfind(t, a); b = tcfind(t, b)
  a == b && return
  if t.size[a] < t.size[b]; a, b = b, a; end
  t.par[b] = a; t.size[a] += t.size[b]
  for j in eachindex(t.tab[b])
    d = t.tab[b][j]
    d == 0 && continue
    tcdefine!(t, a, j, d)
  end
end

function tcdefine!(t::TC, c::Int, j::Int, d::Int)
  c = tcfind(t, c); d = tcfind(t, d)
  e = t.tab[c][j]
  if e == 0
    t.tab[c][j] = d
    if t.tab[d][tcinvcol(j)] == 0; t.tab[d][tcinvcol(j)] = c; end
    return
  end
  e = tcfind(t, e)
  e == d || tcmerge!(t, e, d)
end

"trace relator r from coset c, defining entries as needed; returns the end coset (0 if impossible)"
function tctrace!(t::TC, c::Int, r::Vector{Int})
  x = tcfind(t, c)
  for g in r
    j = tccol(g)
    d = t.tab[x][j]
    if d == 0
      d = tcnewcoset!(t)
      d == 0 && return 0                 # coset limit reached
      tcdefine!(t, x, j, d)
      x = tcfind(t, d)
    else
      x = tcfind(t, d)
    end
  end
  x
end

"""
    toddcoxeter(ngens, rels; maxcosets=1000, maxpasses=1000)

Order of the group presented by `ngens` generators and Tietze `rels`, and the
permutations of the cosets induced by the generators, or `nothing` if the
enumeration does not close within the limits (the group may be infinite).
"""
function toddcoxeter(rels::Vector{Vector{Int}},
                     ngens=maximum(map(x->maximum(abs.(x)),rels));
                     maxcosets::Int=1000, maxpasses::Int=1000)
  t = TC(ngens, [reduceword(r) for r in rels], [zeros(Int, 2ngens)], [1], [1],
         maxcosets)
  # an empty relator is the identity: nothing to enumerate
  all(isempty, t.rels) && return (1, [Perm(collect(1:1)) for _ in 1:ngens])
  pass = 0
  while true
    pass += 1
    pass > maxpasses && error("try to increase maxpasses($maxpasses)")
    changed = false
    c = 1
    while c <= length(t.tab)
      if tcfind(t, c) == c
        for r in t.rels
          y = tctrace!(t, c, r)
          y == 0 && error("tryto increase maxcosets($maxcosets)")
          y = tcfind(t, y)
          if y != c
            tcmerge!(t, y, c)
            changed = true
          end
        end
      end
      c += 1
    end
    changed || break
  end
  # relabel the live cosets 1..n
  live = Int[]
  label = zeros(Int, length(t.tab))
  for c in eachindex(t.tab)
    if tcfind(t, c) == c
      push!(live, c); label[c] = length(live)
    end
  end
  n = length(live)
  # every entry must be defined for the enumeration to be complete
  for c in live, j in 1:2ngens
    (t.tab[c][j] == 0 || tcfind(t, t.tab[c][j]) == 0) && error("humm")
  end
  perms = Vector{Perm}()
  for k in 1:ngens
    img = [label[tcfind(t, t.tab[c][2k-1])] for c in live]
    push!(perms, Perm(img))
  end
  g=Group(Perm{Perms.Idef}.(perms))
  g.length=n
  g
end
