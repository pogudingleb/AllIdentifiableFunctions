###############################################################################
# Part 1: Algorithms for computation with subfields of rational functions
###############################################################################

with(PolynomialIdeals):

#------------------------------------------------------------------------------

IdealsEq := proc(A, B)
    # Checks whether polynomial ideals are equal
    return evalb((A subset B) and (B subset A)):
end:

#------------------------------------------------------------------------------

FieldToIdeal := proc(gens)
    # Input: generators of a subfield of the field of rational functions
    # Computes the MSQ ideal of the field with the new variables of the form x_aux
    # See: https://mediatum.ub.tum.de/doc/685465/document.pdf Definition 2.16
    #      https://doi.org/10.1016/j.jsc.2005.09.010          Lemma 2.1
    local all_vars, subs_dupl, all_dupl, common_denom, polys, f, gb:
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
    gb := select(p -> not (t in indets(p)), gb):
    return PolynomialIdeal(gb, variables=all_dupl):
end proc:

#------------------------------------------------------------------------------

FieldCoefficientRestriction := proc(J, msq_for_subfield)
    # Input: J - a polynomial ideals over a field of rational functions
    #        msq_for_subfield - the MSQ ideal for a subfield E of coefficients (see FieldToIdeal)
    # Computes the radical of the restriction of the ideal to the subfield E 
    # in the sense of https://doi.org/10.1016/j.jsc.2005.09.010 (MSQ-paper in what follows)
    #
    # NOTE: unlike the algorithm in the MSQ-paper, we compute the radical, not the restriction itself
    # one can obtain the algorithm MSQ-paper by replacing PrimeDecomposition with PrimaryDecomposition 
    # in the code below
    local poly_vars, coeff_vars, subs_aux, coeff_aux, gens, subs_aux_msq, gens_msq, msq_ideal_aux, 
    msq_components, J_ext, components, primes_to_keep, P, elim_P, comp, cleaned_ideal:

    poly_vars := IdealInfo[Variables](J):
    coeff_vars := IdealInfo[Parameters](J) union IdealInfo[Parameters](msq_for_subfield):
    subs_aux := map(v -> v = parse(cat(v, _aux_aux)), coeff_vars):
    coeff_aux := subs(subs_aux, coeff_vars):

    # list F of polynomials in the notation of the MSQ-paper, page 375
    gens := map(p -> subs(subs_aux, p * denom(p)), IdealInfo[Generators](J)):

    # Substitution to avoid names clashing between the aux variables in the msq ideal and the variables in J
    subs_aux_msq := map(v -> v = parse(cat(v, _aux)), IdealInfo[Variables](msq_for_subfield)):
    gens_msq := map(p -> subs(subs_aux_msq, p), IdealInfo[Generators](msq_for_subfield)):
    msq_ideal_aux := PolynomialIdeal(gens_msq, variables=map(s -> rhs(s), subs_aux_msq)):
    msq_components := [PrimeDecomposition(msq_ideal_aux)]:

    J_ext := PolynomialIdeal([op(gens), op(gens_msq)], variables=poly_vars union coeff_aux): 
    components := [PrimeDecomposition(J_ext)]:
    
    # Selecting prime components as in Remark on page 377 in MSQ-paper
    primes_to_keep := []:
    for P in components do
        # Going through the prime decomposition of the MSQ-deal mimics the
        # working over the restricted field (which is what has been done in the MSQ-paper)
        elim_P := EliminationIdeal(P, coeff_aux):
        for comp in msq_components do
            if elim_P subset comp then
                primes_to_keep := [op(primes_to_keep), P]:
            end if:
        end do: 
    end do:
    if nops(primes_to_keep) > 0 then
        cleaned_ideal := Intersect(op(primes_to_keep)):
    else
        cleaned_ideal := PolynomialIdeal([0], variables = poly_vars):
    end if:

    # Applying Lemma 2.2 from the MSQ-paper
    return EliminationIdeal(cleaned_ideal, poly_vars):
end proc:

#------------------------------------------------------------------------------

FilterGenerators := proc(ideal)
    # Input: ideal over a rational function field
    # Computes a simplified set of generators of the field of definition
    local gb, gens, p, cf, lc, gsorted, result, big_ideal, cur_ideal, g, new_ideal:
    gb := Groebner[Basis](ideal, tdeg(op(IdealInfo[Variables](ideal)))):
    gens := {}:
    for p in gb do
        cf := sort([coeffs(p, IdealInfo[Variables](ideal))], (a, b) -> length(convert(a, string)) < length(convert(b, string))):
        lc := cf[1]:
        cf := map(c -> c / lc, cf):
        gens := {op(gens), op(cf[2..nops(cf)])}:
    end do:
    gsorted := sort([op(gens)], (a, b) -> length(convert(a, string)) < length(convert(b, string)));
    result := {}:
    big_ideal := FieldToIdeal(gens):
    cur_ideal := FieldToIdeal(result):
    for g in gsorted do
        if big_ideal = cur_ideal then
            return result:
        end if:
        new_ideal := FieldToIdeal({op(result), g}):
        if new_ideal <> cur_ideal then
            result := {op(result), g}:
            cur_ideal := new_ideal:
        end if:
    end do:
    return result:
end proc:

#------------------------------------------------------------------------------

FieldIntersection := proc(gens_left, gens_right)
    # Input: gens_left and gens_right - generators of a subfields E and F of a field of rational functions
    # If terminates, resturns an ideal with field of definition being the intersection of generators of E and F
    # Is guaranteed to terminate if at least one of E and F is algebraically closed in rational functions (see REF)
    # Implementation below is a version of Algorithm 2.38 from https://mediatum.ub.tum.de/doc/685465/document.pdf
    local msq_left, msq_right, poly_vars, coeff_vars, Ii, Ji, gb, result, p, cf, lc;

    msq_left := FieldToIdeal(gens_left):
    msq_right := FieldToIdeal(gens_right):
    poly_vars := IdealInfo[Variables](msq_left) union IdealInfo[Variables](msq_right):
    coeff_vars := IdealInfo[Parameters](msq_left) union IdealInfo[Parameters](msq_right):

    Ii := PolynomialIdeal([1], variables=poly_vars):
    Ji := msq_left:

    while not IdealsEq(Ii, Ji) do
        Ii := FieldCoefficientRestriction(Ji, msq_right):
        Ji := FieldCoefficientRestriction(Ii, msq_left):
    end do:

    return Ji:
end proc:

#------------------------------------------------------------------------------

CompareFields := proc(gens_l, gens_r)
    # Checks whether gens_l and gens_r generate the same subfield of rational functions
    return IdealsEq(FieldToIdeal(gens_l), FieldToIdeal(gens_r)):
end proc:

#------------------------------------------------------------------------------

###############################################################################
# Part 2: Algorithm for computing  identifiable functions
###############################################################################

with(DifferentialAlgebra):

#------------------------------------------------------------------------------

ExtractDenominator := proc(model)
    # Input: model a list of rational functions
    # returns the function multiplied by their denominators and 
    # an inequality corresponding to the LCM of the denominators
    local common_denom, r;
    common_denom := 1:
    for r in model do
        common_denom := lcm(common_denom, denom(r)):
    end do:
    return [op(map(p -> denom(p) * p, model)), common_denom <> 0]:
end proc:

#------------------------------------------------------------------------------

CheckReducibilitySet := proc(polys, charset)
    # Input: polys - list of differential polynomials
    #        charset - a characteristic set
    # Returns true iff all the polynomials are reduced to zero wrt to the charset
    local result:
    result := true:
    for e in polys do
        if NormalForm(e, charset) <> 0 then
            result := false:
            break:
        end if:
    end do:
    return result:
end proc:

#------------------------------------------------------------------------------

GetIOEquations := proc(model, states, ios, params, infolevel)
    # Input: model - list of differential polynomials defining the model
    #        states, ios, params - list of names of state variables, input-output variables
    #                              and parameter, respectively
    # Computes a list of input-output equations of the model
    local Relim, Rorig, charsets, chset_orig, general_comps, general, c, e, gen_comp, io_eqs:

    Relim := DifferentialRing(blocks = [[op(states)], op(ios)], derivations = [t], parameters = params):
    Rorig := DifferentialRing(blocks = [op(ios), op(states)], derivations = [t], parameters = params):
    chset_orig := RosenfeldGroebner(model, Rorig)[1]:
 
    if infolevel > 0 then
        printf("    Computing the characteristic set (singsol = none)\n"):
    end if:
    charsets := RosenfeldGroebner(model, Relim, singsol = none):
    
    if CheckReducibilitySet(Equations(charsets[1]), chset_orig) then
        gen_comp := charsets[1]:
    else
        if infolevel > 0 then
            printf("    Did not pick the right component. Using singsol = all\n"):
        end if:
        charsets := RosenfeldGroebner(model, Relim):
    
        if infolevel > 0 then
            printf("    Selecting the general component\n"):
        end if:
        general_comps := []:
        for c in charsets do
            general := true:
            for e in Equations(c) do
                if NormalForm(e, chset_orig) <> 0 then
                    general := false:
                    break:
                end if:
            end do:
            if general then
                general_comps := [op(general_comps), c]:
            end if:
        end do:
        if nops(general_comps) > 1 then
            print("More than one component is considered general!", general_comps):
        end if:
        gen_comp := general_comps[1]:
    end if:
    io_eqs := Equations(gen_comp, leader < parse(cat(states[-1], "(t)"))):
    return io_eqs:
end proc:

#------------------------------------------------------------------------------

# Adapted from 
# https://www.mapleprimes.com/questions/36772-Extract-Specific-Coefficients-Of-A-Multivariate
# by Kitonum 15364
coefff:=proc(P, t)
    local L, H, i, k:
    L:=[coeffs(P, indets(P), 'h')]: H:=[h]: k:=0:
    for i from 1 to nops(H) do
        if H[i]=t then k:=L[i] fi:
    end do:
    return k;
end proc:

#------------------------------------------------------------------------------

DecomposePolynomial := proc(p, vars_main, vars_coef, infolevel)
    # Input: p - polynomial in two groups of variables: vars_main and vars_coef
    # Computes a decomposition of minimal length of p as a linear combination 
    # of products A * B, where A is a polynomial in vars_main and B 
    # is a polynomial in vars_coef return two lists: list of A's and list of B's
    local cf, monoms, result_cf, result_monom, i, c, m, j, lc, lm, coeff_in_c:
    cf := [coeffs(collect(p, vars_main, 'distributed'), vars_main, 'monoms')]:
    monoms := [monoms]:
    result_cf := []:
    result_monom := Array([]):
    if infolevel > 0 then
        printf("        Number of monomials: %a\n", nops(cf)):
    end:
    for i from 1 to nops(cf) do
        c := cf[i]:
        m := monoms[i]:
        for j from 1 to nops(result_cf) do
            lc, lm := Groebner[LeadingTerm](result_cf[j], plex(op(vars_coef))):
            coeff_in_c := coefff(c, lm):
            c := c - coeff_in_c / lc * result_cf[j]:
            result_monom[j] := result_monom[j] + coeff_in_c / lc * m:
        end do:
        if c <> 0 then
            result_cf := [op(result_cf), c]:
            ArrayTools[Append](result_monom, m):
        end if:
    end do:
    if infolevel > 0 then
        printf("        Reduced to: %a\n", nops(result_cf)):
    end:
    return [result_cf, convert(result_monom, list)]:
end proc:

#------------------------------------------------------------------------------

ConstructWronskian := proc(io_eq, model, states, ios, params, subs_param, infolevel)
    # Input - the same as for GetIOEquations + one IO-equation + flag subs_param
    # Computes the Wronskian for this equation using the representation
    # given by DecomposePolynomial. Return a pair of the Wronskian
    # reduced modulo the original system and a list of coefficients
    # of the compressed io_eq
    local diff_to_ord, v, vt, h, v_ord, vd, p, decomp, diff_polys, Rorig, chset_orig,
    M, yus, yus_reduced, M_sub:

    diff_to_ord := {}:
    for v in ios do
        vt := parse(cat(v, "(t)")):
        diff_to_ord := {op(diff_to_ord), vt = v}:
        for h from 1 to nops(states) + 1 do
            v_ord := parse(cat(v, "_", h)):
            vd := diff(vt, t$h):
            diff_to_ord := {op(diff_to_ord), vd = v_ord}:
        end do:
    end do:
    p := subs(diff_to_ord, io_eq):
    if infolevel > 0 then
        printf("    Combining monomials to reduce the dimension\n"):
    end if:
    decomp := DecomposePolynomial(p, map(e -> rhs(e), diff_to_ord), params, infolevel):
    diff_polys := map(p -> subs(map(e -> rhs(e) = lhs(e), diff_to_ord), p), decomp[2]):
    Rorig := DifferentialRing(blocks = [op(ios), op(states)], derivations = [t], parameters = params):
    chset_orig := RosenfeldGroebner(model, Rorig)[1]:
    
    if infolevel > 0 then
        printf("    Computing the Wronskian\n"):
    end if:
    M := VectorCalculus[Wronskian](diff_polys, t):
    yus := indets(M) minus {t}:

    if infolevel > 0 then
        printf("    Reducing the Wronskian\n"):
    end if:
    if subs_param then
        roll := rand(1..15):
        params_sub := map(v -> v = roll(), params):
        Rorig := DifferentialRing(blocks = [op(ios), op(states)], derivations = [t]):
        chset_orig := RosenfeldGroebner(subs(params_sub, model), Rorig)[1]:
        M := subs(params_sub, M):
    else
        Rorig := DifferentialRing(blocks = [op(ios), op(states)], derivations = [t], parameters = params):
        chset_orig := RosenfeldGroebner(model, Rorig)[1]:
    end if:
    yus_reduced := map(p -> p = NormalForm(p, chset_orig), yus):
    M_sub := subs(yus_reduced, M):
    M_sub := subs(map(x -> parse(cat(x, "(t)")) = x, states), M_sub):
    return [M_sub, decomp[1]]:
end proc:

#------------------------------------------------------------------------------

SingleExperimentIdentifiableFunctions := proc(model, {infolevel := 0})
    # Input: model - ODE model represented as a list of differential polynomials
    # Computes generators of the field of single-identifiable functions
    local states, ios, params, io_eqs, si_gens, eq, wrnsk, echelon_form: 

    states, ios, params := op(ParseInput(model)):

    # Step 1
    if infolevel > 0 then
        printf("Step 1: Computing input-output equations\n"):
    end if:
    model_denomfree := ExtractDenominator(model):
    io_eqs := GetIOEquations(model_denomfree, states, ios, params, infolevel):
    if infolevel > 0 then
        printf("Total number of io-equations: %a\n", nops(io_eqs)):
    end if:
 
    si_gens := {}:
    for eq in io_eqs do
        # Step 2
        if infolevel > 0 then
            printf("Step 2: Constructing the Wronskian\n"):
        end if:
        wrnsk := ConstructWronskian(eq, model_denomfree, states, ios, params, false, infolevel)[1]:
        # Step 3
        if infolevel > 0 then
            printf("Step 3: Computing the reduced row echelon form of the Wronskian\n"):
        end if:
        echelon_form := LinearAlgebra[ReducedRowEchelonForm](wrnsk):
        si_gens := {op(si_gens), op(select(x -> not type(x, numeric), convert(echelon_form, list)))}:
    end do:

    # Step 4
    if infolevel > 0 then
        printf("Step 4: restricting to the field of parameters"):
    end if:
    return FilterGenerators(FieldIntersection(si_gens, params)):
end proc:

#------------------------------------------------------------------------------

# Adapted from https://github.com/pogudingleb/SIAN
FunctionToVariable := proc(f):
    convert(convert(f, string)[1..-4], symbol):
end proc:

ParseInput := proc(model)
   local all_symbols, x_functions, io_functions, params, states, ios:
   all_symbols := foldl(`union`, op( map(e -> indets(e), model) )) minus {t}:
   x_functions := map(f -> int(f, t), select( f -> type(int(f, t), function(name)), all_symbols )):
   io_functions := select( f -> not type(int(f, t), function(name)) and type(f, function(name)) and not f in x_functions, all_symbols ):
   params := [op(select(f -> not type(f, function(name)) and not type(int(f, t), function(name)), all_symbols))]:
   states := [op(map(FunctionToVariable, x_functions))]:
   ios := [op(map(FunctionToVariable, io_functions))]:
   return [states, ios, params]:
end proc:

#------------------------------------------------------------------------------

MultiExperimentIdentifiableFunctions := proc(model, {infolevel := 0, simplified_generators := false, no_bound := false})
    # Input: model - ODE model defined as a list of differential polynomials
    # Computes the bound from Theorem REF. Returns a list containing
    # 1. the bound
    # 2. the list of lists of coefficients of the IO-equations (after compression)
    # 3. (optional) simplified set of generators of the field of multi-experiment identifiable functions

    local states, ios, params, io_eqs, io_coeffs, io_coef, wrnsk, s, roll, wrnsk_sub, r, bound, 
    generators, i, eq, result:
    states, ios, params := op(ParseInput(model)):

    if infolevel > 0 then
        printf("Computing input-output equations\n"):
    end if:
    model_denomfree := ExtractDenominator(model):
    io_eqs := GetIOEquations(model_denomfree, states, ios, params, infolevel):
    if infolevel > 0 then
        printf("Total number of io-equations: %a\n", nops(io_eqs)):
    end if:
 
    io_coeffs := []:
    bound := 0:
    for eq in io_eqs do
        if no_bound then
            io_coef := DecomposePolynomial(eq, indets(eq) minus {op(params)}, params, infolevel)[1]:
        else
            if infolevel > 0 then
                printf("Constructing the Wronskian\n"):
            end if:
            wrnsk, io_coef := op(ConstructWronskian(eq, model_denomfree, states, ios, params, true, infolevel)):
    
            # in the notation of the theorem
            s := nops(io_coef) - 1:
            # substitution does not increase the rank, so the resulting bound will be correct
            roll := rand(1..15):
            wrnsk_sub := map(v -> v = roll(), indets(wrnsk)):
            r := LinearAlgebra[Rank](subs(wrnsk_sub, wrnsk)):
            bound := max(bound, s - r + 1)
        end if:
        io_coeffs := [op(io_coeffs), io_coef]:
    end do:

    result := [bound, io_coeffs]:

    if simplified_generators then
        generators := {}:
        for io_coef in io_coeffs do
             io_coef := sort(io_coef, (a, b) -> length(convert(a, string)) < length(convert(b, string)));
             for i from 2 to nops(io_coef) do
                 generators := {op(generators), io_coef[i] / io_coef[1]}:
             end do:
        end do:
        result := [op(result), FilterGenerators(FieldToIdeal(generators))]:
    end if:

    return result:
end proc:

#------------------------------------------------------------------------------
