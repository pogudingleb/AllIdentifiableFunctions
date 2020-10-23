with(DifferentialAlgebra):

#------------------------------------------------------------------------------

GetIOEquations := proc(model, states, outputs, inputs, params)
    Relim := DifferentialRing(blocks = [[op(states)], op(outputs), op(inputs)], derivations = [t], parameters = params):
    Rorig := DifferentialRing(blocks = [op(inputs), op(outputs), op(states)], derivations = [t], parameters = params):
    charsets := RosenfeldGroebner(model, Relim):
    chset_orig := RosenfeldGroebner(model, Rorig)[1]:
    
    # picking the general component
    general_comps := []:
    for c in charsets do
        general := true:
        for e in Equations(c) do
            if NormalForm(e, chset_orig) <> 0 then
                print(e, NormalForm(e, chset_orig)):
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

DecomposePolynomial := proc(p, vars_main, vars_coef)
    cf := [coeffs(collect(p, vars_main, 'distributed'), vars_main, 'monoms')]:
    result_cf := []:
    result_monom := []:
    for i from 1 to nops(cf) do
        c := cf[i]:
        m := monoms[i]:
        for j from 1 to nops(result_cf) do
            lc, lm := Groebner[LeadingTerm](result_cf[j], plex(op(vars_coef))):
            coeff_in_c := coeff(c, lm):
            c := c - coeff_in_c / lc * result_cf[j]:
            result_monom[j] := result_monom[j] + coeff_in_c / lc * m:
        end do:
        if c <> 0 then
            result_cf := [op(result_cf), c]:
            result_monom := [op(result_monom), m]:
        end if:
    end do:
    return [result_cf, result_monom]:
end proc:

#------------------------------------------------------------------------------

ConstructWronskian := proc(io_eq, model, states, outputs, inputs, params)
    diff_to_ord := {}:
    for v in [op(inputs), op(outputs)] do
        vt := parse(cat(v, "(t)")):
        diff_to_ord := {op(diff_to_ord), vt = v}:
        for h from 1 to nops(states) + 1 do
            v_ord := parse(cat(v, "_", h)):
            vd := diff(vt, t$h):
            diff_to_ord := {op(diff_to_ord), vd = v_ord}:
        end do:
    end do:
    p := subs(diff_to_ord, io_eq):
    decomp := DecomposePolynomial(p, map(e -> rhs(e), diff_to_ord), params):
    diff_polys := map(p -> subs(map(e -> rhs(e) = lhs(e), diff_to_ord), p), decomp[2]):
    Rorig := DifferentialRing(blocks = [op(inputs), op(outputs), op(states)], derivations = [t], parameters = params):
    chset_orig := RosenfeldGroebner(model, Rorig)[1]:
    
    M := VectorCalculus[Wronskian](diff_polys, t):
    yus := indets(M) minus {t}:
    yus_reduced := map(p -> p = NormalForm(p, chset_orig), yus):
    M_sub := subs(yus_reduced, M):
    return [M_sub, decomp[1]]:
end proc:

#------------------------------------------------------------------------------
