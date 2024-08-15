function [LAMBDA, toplength, bottomlength] = discrete_nabla_matrix(h, h_T, alpha, orientation)
    % calculates the stiffness matrix LAMBDA as well as the length of the
    % top and bottom edge of the plate.
    wall_incline = tan(alpha / 360 * 2 * pi);
    t11 = [-(-1)^(orientation)*wall_incline; (-1)^(orientation)*wall_incline; 1];
    t11_v = rot_mat(atan(t11(2)),0,0) * t11;

    z1 = [h_T/2, -t11_v(1)/t11_v(3) * h_T/2];
    z4 = -z1;
    z3 = [h_T/2, h/2];
    z6 = [-h_T/2, h/2];
    z2 = 0.5 * (z3 + z1);
    z5 = 0.5 * (z4 + z6);

    z = [z1', z2', z3', z4', z5', z6'];

    toplength = h/2 + t11_v(1)/t11_v(3) * h_T/2;
    bottomlength = h/2 - t11_v(1)/t11_v(3) * h_T/2; 

    x = sym('x', 'real');

    % Q32 Basis
    temp_Lhat = (x - z1(2)) * (x - z2(2)) * (x - z3(2));
    factorhat = subs(diff(temp_Lhat), z1(2));
    Lhat = temp_Lhat / factorhat;

    temp_Lcheck = (x - z4(2)) * (x - z5(2)) * (x - z6(2));
    factorcheck = subs(diff(temp_Lcheck), z4(2));
    Lcheck = temp_Lcheck / factorcheck;

    L1 = (x - z2(2)) / (z1(2) - z2(2)) * (x - z3(2)) / (z1(2) - z3(2));
    L2 = (x - z1(2)) / (z2(2) - z1(2)) * (x - z3(2)) / (z2(2) - z3(2));
    L3 = (x - z1(2)) / (z3(2) - z1(2)) * (x - z2(2)) / (z3(2) - z2(2));

    L4 = (x - z5(2)) / (z4(2) - z5(2)) * (x - z6(2)) / (z4(2) - z6(2));
    L5 = (x - z4(2)) / (z5(2) - z4(2)) * (x - z6(2)) / (z5(2) - z6(2));
    L6 = (x - z4(2)) / (z6(2) - z4(2)) * (x - z5(2)) / (z6(2) - z5(2));

    factorL1 = subs(diff(L1), z1(2));
    factorL2 = subs(diff(L2), z1(2));
    factorL3 = subs(diff(L3), z1(2));

    factorL4 = subs(diff(L4), z4(2));
    factorL5 = subs(diff(L5), z4(2));
    factorL6 = subs(diff(L6), z4(2));
    
    B1 = L1 - factorL1 * Lhat;
    B2 = L2 - factorL2 * Lhat;
    B3 = L3 - factorL3 * Lhat;

    B4 = L4 - factorL4 * Lcheck;
    B5 = L5 - factorL5 * Lcheck;
    B6 = L6 - factorL6 * Lcheck;

    X = linspace(z1(2), z3(2), 100);
    Y = arrayfun(@(t) subs(Lhat, t), X);

    %figure(1);
    %plot(X, Y); hold on
    %plot(z1(2), 0, 'ro');
    %plot(z2(2), 0, 'ro');
    %plot(z3(2), 0, 'ro');
    %plot(z1(2), 1, 'ro');
    %plot(z2(2), 1, 'ro');
    %plot(z3(2), 1, 'ro');
    %plot(X, arrayfun(@(t) subs(B1, t), X));
    %plot(X, arrayfun(@(t) subs(B2, t), X));
    %plot(X, arrayfun(@(t) subs(B3, t), X)); hold off

    %figure(2);
    %X = linspace(z4(2), z6(2), 100);

    %plot(z4(2), 0, 'bo'); hold on
    %plot(z5(2), 0, 'bo');
    %plot(z6(2), 0, 'bo');
    %plot(z4(2), 1, 'bo');
    %plot(z5(2), 1, 'bo');
    %plot(z6(2), 1, 'bo');
    %plot(X, arrayfun(@(t) subs(Lcheck, t), X));
    %plot(X, arrayfun(@(t) subs(B4, t), X));
    %plot(X, arrayfun(@(t) subs(B5, t), X));
    %plot(X, arrayfun(@(t) subs(B6, t), X)); hold off

    B1d = diff(B1);
    B2d = diff(B2);
    B3d = diff(B3);
    Lhatd = diff(Lhat);

    B4d = diff(B4);
    B5d = diff(B5);
    B6d = diff(B6);
    Lcheckd = diff(Lcheck);

    Lambda_1 = zeros(3, 4);
    Lambda_2 = zeros(3, 4);

    Lambda_1(1,:) = [0, 0, 0, 1];

    Lambda_1(2, 1) = subs(B1d, z2(2));
    Lambda_1(2, 2) = subs(B2d, z2(2));
    Lambda_1(2, 3) = subs(B3d, z2(2));
    Lambda_1(2, 4) = subs(Lhatd, z2(2));

    Lambda_1(3, 1) = subs(B1d, z3(2));
    Lambda_1(3, 2) = subs(B2d, z3(2));
    Lambda_1(3, 3) = subs(B3d, z3(2));
    Lambda_1(3, 4) = subs(Lhatd, z3(2));

    Lambda_2(1,:) = [0, 0, 0, 1];

    Lambda_2(2, 1) = subs(B4d, z5(2));
    Lambda_2(2, 2) = subs(B5d, z5(2));
    Lambda_2(2, 3) = subs(B6d, z5(2));
    Lambda_2(2, 4) = subs(Lcheckd, z5(2));

    Lambda_2(3, 1) = subs(B4d, z6(2));
    Lambda_2(3, 2) = subs(B5d, z6(2));
    Lambda_2(3, 3) = subs(B6d, z6(2));
    Lambda_2(3, 4) = subs(Lcheckd, z6(2));

    LAMBDA = [Lambda_1, zeros(3, 4); zeros(3,4), Lambda_2];    
end

function Q = rot_mat(a,b,c)
    Q = [1, 0, 0; 0, cos(a), -sin(a); 0, sin(a), cos(a)];
    Q = Q * [cos(b), 0, sin(b); 0, 1, 0; -sin(b), 0, cos(b)];
    Q = Q * [cos(c), -sin(c), 0; sin(c), cos(c), 0; 0,0,1];
end