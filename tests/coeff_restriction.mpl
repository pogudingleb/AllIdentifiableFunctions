read "../ComputeIdentifiableFunctions.mpl":

cases := [
    [
        [[a * x + b * y + c * y^2, d * y], [a, b, c, d], [c, d]],
        [x, y]
    ],
    [
        [[a * x + b * z + c * y^2, d * y], [a, b, c, d], [c, d]],
        [a * x + b * z, y]
    ],
    # MQS, Example page 7, radical version
    [
        [[x1 * z1  - z2, x2 * z2 - z3, z3^3], [x1, x2], []],
        [z1, z2, z3]
    ]
];

num_passed := 0:
num_failed := 0:

for case in cases do
    input := case[1]:
    correct := case[2]:
    poly_vars := {op(indets(input[1]))} minus {op(input[2])}:
    if Groebner[Basis](FieldCoefficientRestriction(op(input)), tdeg(op(poly_vars))) = Groebner[Basis](correct, tdeg(op(poly_vars))) then
        printf("PASSED");
        num_passed := num_passed + 1:
    else
        printf("FAILED");
        num_failed := num_failed + 1:
        print("Expected: ", Groebner[Basis](FieldCoefficientRestriction(op(input)), tdeg(op(poly_vars))));
        print("Got: ", Groebner[Basis](correct, tdeg(op(poly_vars))));
    end if:
end do:

printf("Passed: %a, failed %a \n", num_passed, num_failed);
