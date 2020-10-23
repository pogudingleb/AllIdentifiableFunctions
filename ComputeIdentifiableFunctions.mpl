with(DifferentialAlgebra):

read cat(currentdir(), "/FieldElimination.mpl"):

#------------------------------------------------------------------------------

GetIOEquations := proc(model, states, ios, params)
    Relim := DifferentialRing(blocks = [[op(states)], op(ios)], derivations = [t], parameters = params):
    Rorig := DifferentialRing(blocks = [op(ios), op(states)], derivations = [t], parameters = params):
    printf("    Computing the characteristic set\n"):
    charsets := RosenfeldGroebner(model, Relim):
    chset_orig := RosenfeldGroebner(model, Rorig)[1]:
    
    # picking the general component
    printf("     Selecting the general component\n"):
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

DecomposePolynomial := proc(p, vars_main, vars_coef)
    cf := [coeffs(collect(p, vars_main, 'distributed'), vars_main, 'monoms')]:
    monoms := [monoms]:
    result_cf := []:
    result_monom := []:
    printf("        Number of monomials: %a\n", nops(cf)):
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
            result_monom := [op(result_monom), m]:
        end if:
    end do:
    printf("        Reduced to: %a\n", nops(result_cf)):
    return [result_cf, result_monom]:
end proc:

#------------------------------------------------------------------------------

ConstructWronskian := proc(io_eq, model, states, ios, params)
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
    printf("    Combining monomials to reduce the dimension\n"):
    decomp := DecomposePolynomial(p, map(e -> rhs(e), diff_to_ord), params):
    diff_polys := map(p -> subs(map(e -> rhs(e) = lhs(e), diff_to_ord), p), decomp[2]):
    print(diff_polys):
    Rorig := DifferentialRing(blocks = [op(ios), op(states)], derivations = [t], parameters = params):
    chset_orig := RosenfeldGroebner(model, Rorig)[1]:
    
    printf("    Computing the Wronskian\n"):
    M := VectorCalculus[Wronskian](diff_polys, t):
    yus := indets(M) minus {t}:
    printf("    Reducing the Wronskian\n"):
    yus_reduced := map(p -> p = NormalForm(p, chset_orig), yus):
    M_sub := subs(yus_reduced, M):
    M_sub := subs(map(x -> parse(cat(x, "(t)")) = x, states), M_sub):
    return [M_sub, decomp[1]]:
end proc:

#------------------------------------------------------------------------------

SingleExperimentIdentifiableFunction := proc(model)
    model_data := ParseInput(model):
    states := model_data[1]:
    ios := model_data[2]:
    params := model_data[3]:

    # Step 1
    printf("Step 1: Computing input-output equations\n"):
    io_eqs := GetIOEquations(model, states, ios, params):
    printf("Total number of io-equations: %a\n", nops(io_eqs)):
 
    si_gens := {}:
    for eq in io_eqs do
        # Step 2
        printf("Step 2: Constructing the Wronskian\n"):
        wrnsk := ConstructWronskian(eq, model, states, ios, params)[1]:
        # Step 3
        printf("Step 3: Computing the reduced row echelon form of the Wronskian\n"):
        echelon_from := LinearAlgebra[ReducedRowEchelonForm](wrnsk):
        si_gens := {op(si_gens), op(select(x -> not type(x, numeric), convert(echelon_from, list)))}:
    end do:

    return FieldElimination(si_gens, params):
end proc:

#------------------------------------------------------------------------------

# Adapted from https://github.com/pogudingleb/SIAN
FunctionToVariable := proc(f):
    convert(convert(f, string)[1..-4], symbol):
end proc:

ParseInput := proc(model)
   all_symbols := foldl(`union`, op( map(e -> indets(e), model) )) minus {t}: 
   x_functions := map(f -> int(f, t), select( f -> type(int(f, t), function(name)), all_symbols )):
   io_functions := select( f -> not type(int(f, t), function(name)) and type(f, function(name)) and not f in x_functions, all_symbols ):
   params := [op(select(f -> not type(f, function(name)) and not type(int(f, t), function(name)), all_symbols))]:
   states := [op(map(FunctionToVariable, x_functions))]:
   ios := [op(map(FunctionToVariable, io_functions))]:
   return [states, ios, params]:
end proc:

#------------------------------------------------------------------------------
