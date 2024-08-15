function test_crease
    k = 2; % determines the size of the configuration
    h = 2^(-k); % size parameter
    N = 2 / h;
    h_T = h/4; % height parameter for a cross-like structure
    alpha = 20.6;
    alpha_rad = alpha / 360 * 2 * pi;
    wall_incline = tan(alpha_rad);
    i = 1;
    j = 2;
    t11 = [-(-1)^(i+j)*wall_incline; (-1)^(i+j)*wall_incline; 1]
    t11_v = rot_mat(atan(t11(2)),0,0) * t11

    z1 = [h_T/2, -t11_v(1)/t11_v(3) * h_T/2];
    z4 = -z1;
    z3 = [h_T/2, h/2];
    z6 = [-h_T/2, h/2];
    z2 = 0.5 * (z3 + z1);
    z5 = 0.5 * (z4 + z6);

    z = [z1', z2', z3', z4', z5', z6'];
end

function Q = rot_mat(a,b,c)
    Q = [1, 0, 0; 0, cos(a), -sin(a); 0, sin(a), cos(a)];
    Q = Q * [cos(b), 0, sin(b); 0, 1, 0; -sin(b), 0, cos(b)];
    Q = Q * [cos(c), -sin(c), 0; sin(c), cos(c), 0; 0,0,1];
end