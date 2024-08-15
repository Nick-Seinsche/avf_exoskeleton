function [S, B, corrector] = stiffness_and_load(h, h_T, alpha, orientation)
    % calculates the Stiffness matrix and Load vector.
    % Takes s input:
    % size parameters h, h_T > 0
    % wall inclination alpha
    % and orientation (either 1 or 2)
    % Returns:
    % - Stiffness Matrix S for the canonical basis on G_D_h
    % - Load vector B where B_i is the integral of the i-th canonical basis
    % vector on the plate.
    % - the corrector which is a 2d transformation that transformes
    % the vector starting at z_5 pointing in the direction where
    % any function of W_h,h_T is affine; to the first unit vector.
    % (see thesis)
    wall_incline = tan(alpha / 360 * 2 * pi);
    t11 = [-(-1)^(orientation)*wall_incline; (-1)^(orientation)*wall_incline; 1];
    t11_v = rot_mat(atan(t11(2)),0,0) * t11;

    z1 = [h_T/2, -t11_v(1)/t11_v(3) * h_T/2];
    z4 = -z1;
    z3 = [h_T/2, h/2];
    z6 = [-h_T/2, h/2];
    z2 = 0.5 * (z3 + z1);
    z5 = 0.5 * (z4 + z6);

    % create simple triangulation of plate
    z = [z1', z2', z3', z4', z5', z6'];
    c4n = z';
    n4e = [4, 1, 2; 5, 4, 2; 5, 2, 3; 6, 5, 3];

    [c4n, n4e, ~, ~] = red_refine(c4n, n4e, [], []);
    [c4n, n4e, ~, ~] = red_refine(c4n, n4e, [], []);

    nC = size(c4n, 1);
    x = sym('x', 'real');
    y = sym('y', 'real');
    t = (y + h_T/2) / h_T;

    max_shift = -t11_v(1)/t11_v(3) * h_T;
    c = max_shift / (h/2 - 2 * (t - 0.5) * max_shift / 2);

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
    
    % assemble basis vectors and their gradients in a matrix
    DL = [DL1, DL2, DL3, DL4, DL5, DL6];
    LL = [L1; L2; L3; L4; L5; L6];
    
    submatrix = [2,3,5,6];
    sz = size(submatrix, 2);
    index = 1:sz;
    S = zeros(sz, sz);  % stiffness matrix
    B = zeros(sz, 1);  % load vector

    formatSpec = "Calculating Stiffness Matrix and Load Vector ... (%d percent)\n";

    for i=index
        for j=index
            fprintf(formatSpec, round((j + sz * (i - 1))/(sz^2)*100));
            if i <= j
                product = DL(:,submatrix(i))' * DL(:,submatrix(j));
                U = arrayfun(@(k) subs(subs(product, y, c4n(k,1)), x, c4n(k,2)), 1:nC);

                [~, integral] = integrate_fun_trimesh(c4n', n4e', U);
                S(i,j) = integral;
                if i < j
                    S(j,i) = S(i,j);
                end
            end
        end

        U = arrayfun(@(k) eval(subs(subs(LL(submatrix(i)), y, c4n(k,1)), x, c4n(k,2))), 1:nC);
        [~, integral] = integrate_fun_trimesh(c4n', n4e', U);
        B(i) = integral;
    end
    
    % calculate the corrector matrix
    corrector = [1, 0 ; eval(subs(subs(g_plus, y, -h_T/2), x, z5(2))) , h_T];
end

function Q = rot_mat(a,b,c)
    Q = [1, 0, 0; 0, cos(a), -sin(a); 0, sin(a), cos(a)];
    Q = Q * [cos(b), 0, sin(b); 0, 1, 0; -sin(b), 0, cos(b)];
    Q = Q * [cos(c), -sin(c), 0; sin(c), cos(c), 0; 0,0,1];
end