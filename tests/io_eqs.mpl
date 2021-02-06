read "../ComputeIdentifiableFunctions.mpl";

cases := [
    [
        [[
            diff(x1(t), t) - (a * x1(t) - b * x1(t) * x2(t)),
            diff(x2(t), t) - (-c * x2(t) + d * x1(t) * x2(t) + u(t)),
            y(t) = x1(t)
        ],
        [x1, x2], [y], [u], [a, b, c, d]],
        [y(t) * diff(diff(y(t), t), t) - diff(y(t), t)^2 - y(t)^2 * diff(y(t), t) * d + y(t) * diff(y(t), t) * c + a * d * y(t)^3 + u(t) * y(t)^2 * b - a * c * y(t)^2]
    ],
    [
        [[
            diff(xa(t), t) + k1 * xa(t),
            diff(xb(t), t) - k1 * xa(t) + k2 * xb(t),
            diff(xc(t), t) - k2 * xb(t),
            diff(eA(t), t),
            diff(eC(t), t),
            y1(t) - xc(t),
            y2(t) - eA(t) * xa(t) - eB * xb(t) - eC(t) * xc(t),
            y3(t) - eA(t), 
            y4(t) - eC(t)
        ],
        [xa, xb, xc, eA, eC], [y2, y1, y3, y4], [], [k1, k2, eB]],
        [
            k1 * k2 * (y2(t) - y1(t) * y4(t)) - eB * k1 * diff(y1(t), t) - k2 * y3(t) * diff(y1(t), t) - y3(t) * diff(diff(y1(t), t), t),
            diff(diff(diff(y1(t), t), t), t) + (k1 + k2) * diff(diff(y1(t), t), t) + k1 * k2 * diff(y1(t), t),
            diff(y3(t), t),
            diff(y4(t), t)
        ]
    ] 
]:

num_passed := 0:
num_failed := 0:

for case in cases do
    input := case[1]:
    print(input);
    correct := case[2]:
    result := GetIOEquations(op(input), 0):
    passed := true:
    if nops(result) <> nops(correct) then
        passed := false:
    else
        for i from 1 to nops(result) do
            if simplify(result[i] - correct[i]) <> 0 then
                passed := false: 
            end if:
        end do:
    end if:
    if passed = true then
        printf("PASSED");
        num_passed := num_passed + 1:
    else
        printf("FAILED");
        print(result);
        print(correct);
        num_failed := num_failed + 1:
    end if:
end do:

printf("Passed: %a, failed %a \n", num_passed, num_failed);
