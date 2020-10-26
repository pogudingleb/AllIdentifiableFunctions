read "../ComputeIdentifiableFunctions.mpl":

cases := [
    [
        [PolynomialIdeal([a * x + b * y + c * y^2, d * y], variables={x, y}), PolynomialIdeal([a - a_aux, b - b_aux], variables={a_aux, b_aux})],
        PolynomialIdeal([x, y])
    ],
    [
        [PolynomialIdeal([a * x + b * z + c * y^2, d * y], variables={x, y, z}), PolynomialIdeal([a - a_aux, b - b_aux], variables={a_aux, b_aux})],
        PolynomialIdeal([a * x + b * z, y], variables={x, y, z})
    ],
    [
        [PolynomialIdeal([a - a_aux], variables={a_aux}), PolynomialIdeal([a - a_aux, b - b_aux], variables={a_aux, b_aux})],
        PolynomialIdeal([a - a_aux], variables={a_aux})
    ],
    [
        [PolynomialIdeal([-a - b + a_aux + b_aux, -a * a_aux + a * b + a_aux^2 - a_aux * b], variables={a_aux, b_aux}), PolynomialIdeal([-a + a_aux, -b + b_aux], variables={a_aux, b_aux})],
        PolynomialIdeal([-a - b + a_aux + b_aux, -a * a_aux + a * b + a_aux^2 - a_aux * b], variables={a_aux, b_aux})
    ],
    [
        [
            PolynomialIdeal([-a - b + a_aux + b_aux, -a * a_aux + a * b + a_aux^2 - a_aux * b], variables={a_aux, b_aux}), 
            PolynomialIdeal([-a - b + a_aux + b_aux, -a * a_aux + a * b + a_aux^2 - a_aux * b, c - c_aux], variables={a_aux, b_aux, c_aux})
        ],
        PolynomialIdeal([-a - b + a_aux + b_aux, -a * a_aux + a * b + a_aux^2 - a_aux * b], variables={a_aux, b_aux})
    ]
];

num_passed := 0:
num_failed := 0:

for case in cases do
    input := case[1]:
    correct := case[2]:
    if IdealsEq(FieldCoefficientRestriction(op(input)), correct) then
        printf("PASSED\n");
        num_passed := num_passed + 1:
    else
        printf("FAILED\n");
        num_failed := num_failed + 1:
        print("Expected: ", correct);
        print("Got: ", FieldCoefficientRestriction(op(input)));
    end if:
end do:

printf("Passed: %a, failed %a \n", num_passed, num_failed);
