function test_stiffness()
    h = 1/2;
    h_T = 1/4;
    alpha = 20.9;
    wall_incline = tan(alpha / 360 * 2 * pi);
    t11 = [wall_incline; -wall_incline; 1];

    z1 = [h_T/2, t11(2)/t11(3) * h_T/2];
    z4 = -z1;
    z3 = [h_T/2, h/2];
    z6 = [-h_T/2, h/2];
    z2 = 0.5 * (z3 + z1);
    z5 = 0.5 * (z4 + z6);

    z = [z1', z2', z3', z4', z5', z6'];

    c4n = z';
    n4e = [4, 1, 2; 5, 4, 2; 5, 2, 3; 6, 5, 3];

    [c4n, n4e, ~, ~] = red_refine(c4n, n4e, [], []);
    [c4n, n4e, ~, ~] = red_refine(c4n, n4e, [], []);
    %[c4n, n4e, ~, ~] = red_refine(c4n, n4e, [], []);
    %[c4n, n4e, ~, ~] = red_refine(c4n, n4e, [], []);
    %[c4n, n4e, ~, ~] = red_refine(c4n, n4e, [], []);

    nC = size(c4n, 1);
    x = sym('x', 'real');
    y = sym('y', 'real');
    t = (y + h_T/2) / h_T;

    c = t11(2)/t11(3) * h_T / (h/2 - 2 * (t - 0.5) * t11(2)/t11(3) * h_T / 2);

    g_plus = (1 - t) * (h/2 - x) * c;
    g_minus = -t * (h/2 - x) * c;

    L1 = t*(1 - (1-t) * c)*(x - z2(2)) / (z1(2) - z2(2)) * (x - z3(2)) / (z1(2) - z3(2));
    L2 = t*(1 - (1-t) * c)*(x - z1(2)) / (z2(2) - z1(2)) * (x - z3(2)) / (z2(2) - z3(2));
    L3 = t*(1 - (1-t) * c)*(x - z1(2)) / (z3(2) - z1(2)) * (x - z2(2)) / (z3(2) - z2(2));

    L4 = (1 - t)*(1 + t * c)*(x - z5(2)) / (z4(2) - z5(2)) * (x - z6(2)) / (z4(2) - z6(2));
    L5 = (1 - t)*(1 + t * c)*(x - z4(2)) / (z5(2) - z4(2)) * (x - z6(2)) / (z5(2) - z6(2));
    L6 = (1 - t)*(1 + t * c)*(x - z4(2)) / (z6(2) - z4(2)) * (x - z5(2)) / (z6(2) - z5(2));

    L1 = subs(L1, x, x + g_plus);
    L2 = subs(L2, x, x + g_plus);
    L3 = subs(L3, x, x + g_plus);

    L4 = subs(L4, x, x + g_minus);
    L5 = subs(L5, x, x + g_minus);
    L6 = subs(L6, x, x + g_minus);

    DL1 = [diff(L1, x); diff(L1, y)];
    DL2 = [diff(L2, x); diff(L2, y)];
    DL3 = [diff(L3, x); diff(L3, y)];

    DL4 = [diff(L4, x); diff(L4, y)];
    DL5 = [diff(L5, x); diff(L5, y)];
    DL6 = [diff(L6, x); diff(L6, y)];

    DL = [DL1, DL2, DL3, DL4, DL5, DL6];
    
    submatrix = 1:6;
    sz = size(submatrix, 2)
    index = 1:sz;
    S = zeros(sz, sz);

    for i=index
        for j=index
            disp("Calculating Stiffness Matrix...")
            disp([num2str(round((j + sz * (i - 1))/(sz^2)*100)), "%"])
            if i <= j
                product = DL(:,submatrix(i))' * DL(:,submatrix(j));
                U = arrayfun(@(k) subs(subs(product, y, c4n(k,1)), x, c4n(k,2)), 1:nC);

                [area, integral] = integrate_fun_trimesh(c4n', n4e', U);
                S(i,j) = integral;
                if i < j
                    S(j,i) = S(i,j);
                end
            end
        end
    end

    for q=index
        det(S(1:q,1:q))
    end


    [area, integral] = integrate_fun_trimesh(c4n', n4e', LL1)

    u = @(y,x) 0.05 * y - 3 * y * x ^ 3;

    U = arrayfun(@(i) u(c4n(i,1), c4n(i,2)), 1:nC);

    figure(1);
    trisurf(n4e,c4n(:,1),c4n(:,2), U);
end