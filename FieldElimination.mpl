FieldCoefficientRestriction := proc(ideal_gens, coeff_vars, vars_to_elim)
    poly_vars := {op(indets(ideal_gens))} minus {op(coeff_vars)}:
    vars_to_keep := {op(coeff_vars)} minus {op(vars_to_elim)}:
    gens := map(poly -> poly * denom(poly), ideal_gens):
    primes := [PolynomialIdeals[PrimeDecomposition](PolynomialIdeals[PolynomialIdeal](op(gens)))]:
    primes_to_keep := []:
    for P in primes do
        if Groebner[Basis](P, tdeg(op(poly_vars))) <> [1] then
            primes_to_keep := [op(primes_to_keep), P]
        end: 
    end do:
    cleaned_ideal := PolynomialIdeals[Intersect](op(primes_to_keep)):
    gb := Groebner[Basis](cleaned_ideal, plex(op(vars_to_elim), op(poly_vars))):
    result := []:
    for poly in gb do
        if {op(indets(poly))} subset (poly_vars union vars_to_keep) then
            result := [op(result), poly]:
        end if:
    end do:
    return result:
end proc:

#------------------------------------------------------------------------------

FieldToIdeal := proc(gens)
    all_vars := indets(gens):
    subs_dupl := map(v -> v = cat(v, _aux), all_vars):
    all_dupl := map(v -> subs(subs_dupl, v), all_vars):
    common_denom := 1:
    polys := []:
    for f in gens do
        common_denom := lcm(common_denom, denom(f)):
        polys := [op(polys), numer(f) * subs(subs_dupl, denom(f)) - subs(subs_dupl, numer(f)) * denom(f)]:
    end do:
    gb := Groebner[Basis]([op(polys), common_denom * t - 1], plex(t, op(all_dupl))):
    return select(p -> not (t in indets(p)), gb):
end proc:

#------------------------------------------------------------------------------

FieldElimination := proc(gens, vars_to_elim)
    all_vars := indets(gens):
    vars_to_keep := {op(all_vars)} minus {op(vars_to_elim)}:
    subs_dupl := map(v -> v = cat(v, _aux), all_vars):
    elim_dupl := subs(subs_dupl, vars_to_elim):
    keep_dupl := subs(subs_dupl, vars_to_keep):

    gb := Groebner[Basis](FieldToIdeal(gens), plex(op(elim_dupl), op(keep_dupl))):
    J_gens := select(p -> {op(indets(p))} subset ({op(all_vars)} union keep_dupl), gb):
    J_restriction := FieldCoefficientRestriction(J_gens, all_vars, vars_to_elim):
    result := {}:
    for p in J_restriction do
        cf := [coeffs(p, keep_dupl)]:
        lc := cf[1]:
        cf := map(c -> c / lc, cf):
        result := {op(result), op(cf[2..nops(cf)])}:
    end do:
    return result:
end proc:

#------------------------------------------------------------------------------

# for testing
CompareFields := proc(gens_l, gens_r)
    return evalb(FieldToIdeal(gens_l) = FieldToIdeal(gens_r)) 
end proc:
