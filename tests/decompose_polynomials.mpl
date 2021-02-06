read "../ComputeIdentifiableFunctions.mpl":

cases := [
    [
        [a * x + a * y, [x, y], [a]], 1
    ],
    [
        [(a + b) * (x + y) - (b + c) * (y + z) + 3 * (a + c) * (x + z), [x, y, z], [a, b, c]], 3
    ]
];

num_passed := 0:
num_failed := 0:

for case in cases do
    input := case[1]:
    correct := case[2]:
    passed := true:
    result := DecomposePolynomial(op(input), 0):
    if nops(result[1]) <> correct then
        passed := false:
    end if:
    orig_poly := 0:
    for i from 1 to nops(result[1]) do
        orig_poly := orig_poly + result[1][i] * result[2][i]:
    end do:
    if simplify(orig_poly - input[1]) <> 0 then
        passed := false:
    end if:
    if passed = true then
        printf("PASSED");
        num_passed := num_passed + 1:
    else
        printf("FAILED");
        num_failed := num_failed + 1:
        print(result);
    end if:
end do:

printf("Passed: %a, failed %a \n", num_passed, num_failed);
