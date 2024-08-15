function test_G_interpolate()
    k = 0;
    orientation = 2;

    h = 2^(-k); N = 2^(k+1); %L = N^2;
    h_T = h/2;
    alpha = 20.9;
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

    c4n = z';
    n4e = [4, 1, 2; 5, 4, 2; 5, 2, 3; 6, 5, 3];

    [c4n, n4e, ~, ~] = red_refine(c4n, n4e, [], []);
    [c4n, n4e, ~, ~] = red_refine(c4n, n4e, [], []);
    [c4n, n4e, ~, ~] = red_refine(c4n, n4e, [], []);
    %[c4n, n4e, ~, ~] = red_refine(c4n, n4e, [], []);
    %[c4n, n4e, ~, ~] = red_refine(c4n, n4e, [], []);
    nC = size(c4n, 1);

    u = G_interpolate_new([0,0,0.25,0,0,0], h, h_T, alpha, 2);

    U = arrayfun(@(i) u(c4n(i,1), c4n(i,2)), 1:nC);

    figure(1);
    trisurf(n4e,c4n(:,1),c4n(:,2),U);
    zlim([-1 1]);
end

function Q = rot_mat(a,b,c)
    Q = [1, 0, 0; 0, cos(a), -sin(a); 0, sin(a), cos(a)];
    Q = Q * [cos(b), 0, sin(b); 0, 1, 0; -sin(b), 0, cos(b)];
    Q = Q * [cos(c), -sin(c), 0; sin(c), cos(c), 0; 0,0,1];
end