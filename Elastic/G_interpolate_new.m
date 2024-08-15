function [u, g_plus, g_minus] = G_interpolate_new(V, h, h_T, alpha, orientation)
    % Takes:
    % - V in R^6
    % - size parameters h, h_T > 0
    % - wall inclination alpha
    % - an orientation (either 1 or 2) which distinguishes northern/western
    % blade from southern/eastern blade

    % Returns:
    % - A function u in G_h,h_T with the values V_i at z_i on the reference
    % blade.
    % - The function g_plus and g_minus (see thesis)
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

    x = sym('x', 'real');
    y = sym('y', 'real');
    t = (y + h_T/2) / h_T;

    L1 = (x - z2(2)) / (z1(2) - z2(2)) * (x - z3(2)) / (z1(2) - z3(2));
    L2 = (x - z1(2)) / (z2(2) - z1(2)) * (x - z3(2)) / (z2(2) - z3(2));
    L3 = (x - z1(2)) / (z3(2) - z1(2)) * (x - z2(2)) / (z3(2) - z2(2));

    L4 = (x - z5(2)) / (z4(2) - z5(2)) * (x - z6(2)) / (z4(2) - z6(2));
    L5 = (x - z4(2)) / (z5(2) - z4(2)) * (x - z6(2)) / (z5(2) - z6(2));
    L6 = (x - z4(2)) / (z6(2) - z4(2)) * (x - z5(2)) / (z6(2) - z5(2));

    max_shift = -t11_v(1)/t11_v(3) * h_T;
    c = max_shift / (h/2 - 2 * (t - 0.5) * max_shift / 2);

    ptop = V(1) * L1 + V(2) * L2 + V(3) * L3; 
    pbot = V(4) * L4 + V(5) * L5 + V(6) * L6; 

    g_plus = (1 - t) * (h/2 - x) * c;
    g_minus = -t * (h/2 - x) * c;

   p = subs(ptop, x, x + g_plus) * t * (1 - (1 - t) * c) + subs(pbot, x,x + g_minus) * (1-t) * (1 + t * c);

   u = @(yy,xx) eval(subs(subs(p, x, xx), y, yy));

   g_plus = @(yy,xx) eval(subs(subs(g_plus, x, xx), y, yy));
   g_minus = @(yy,xx) eval(subs(subs(g_minus, x, xx), y, yy));
end

function Q = rot_mat(a,b,c)
    Q = [1, 0, 0; 0, cos(a), -sin(a); 0, sin(a), cos(a)];
    Q = Q * [cos(b), 0, sin(b); 0, 1, 0; -sin(b), 0, cos(b)];
    Q = Q * [cos(c), -sin(c), 0; sin(c), cos(c), 0; 0,0,1];
end