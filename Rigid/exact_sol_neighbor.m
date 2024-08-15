function exact_sol_neighbor()
    % plots an "exact" 2x2 configuration solution (ie. k = 0) using the
    % omega function and the calculate_neighbor function.

    % "exact" meaning that the rotation maticies are calculated using
    % explicit formulas.
    k = 0;
    h = 2^(-k);
    N = 2^(k+1);
    L = N^2;
    t11 = [0.87; -0.87; 1];
    t11 = 0.5 * h * t11 / norm(t11, 2);

    % major angle (approx.)
    blue_green = 36.9811 / 360 * (2 * pi);
    
    % minor angle (approx.)
    blue_pink = 31.3686 / 360 * (2 * pi);

    % verify: minor angle calculated via the omega function
    omega(blue_green, 1, 0.87/2) - blue_pink

    U0 = init_configuration(k);

    p1 = U0(1:3);
    Q1 = reshape(U0(4:12), 3, 3);

    % calculate x_12 to be east of x_11 with major angle
    [p2_new, Q2_new] = calculate_neighbor(h, t11, 1, "e", p1, Q1, blue_green);
    
    % calculate x_21 to be south of x_22 with minor angle
    [p3_new, Q3_new] = calculate_neighbor(h, t11, 1, "s", p1, Q1, blue_pink);
    
    % calculate x_22 to be south of x_12 with negative minor angle (see
    % "ambiguity measuring angles" for an explanation for the sign)
    [p4_new, Q4_new] = calculate_neighbor(h, t11, 2, "s", p2_new, Q2_new, -blue_pink);

    U = vertcat(p1, ...
                reshape(Q1, [], 1), ...
                ...
                p2_new, ...
                reshape(Q2_new, [], 1), ...
                ...
                p3_new, ...
                reshape(Q3_new, [], 1), ...
                ...
                p4_new, ...
                reshape(Q4_new, [], 1));

    plot_skeleton_square(U, k, t11, 1, 1);
    
    A = assemble_A(k, t11);
    % hinge constraint violation
    norm(A * U, 2)
end