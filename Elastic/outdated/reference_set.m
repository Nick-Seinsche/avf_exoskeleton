function reference_set()
    h = 1/2;
    h_T = 1/4;
    orientation = 2;
    alpha = 21.6;
    wall_incline = tan(alpha / 360 * 2 * pi)
    t11 = [-(-1)^(orientation)*wall_incline; (-1)^(orientation)*wall_incline; 1]
    t11_v = rot_mat(atan(t11(2)),0,0) * t11

    z1 = [h_T/2, -t11_v(1)/t11_v(3) * h_T/2];
    z4 = -z1;
    z3 = [h_T/2, h/2];
    z6 = [-h_T/2, h/2];
    z2 = 0.5 * (z3 + z1);
    z5 = 0.5 * (z4 + z6);

    z = [z1', z2', z3', z4', z5', z6'];

    plot(z(2,:), z(1,:), "bo"); hold on
    line([0, t11(2)], [0, t11(3)], "color", "r")
    line(z(2,[1, 2, 3, 6, 5, 4, 1]), z(1,[1, 2, 3, 6, 5, 4, 1])); hold off
    xlim([-h h]);
    ylim([-h h]);
    set(gca, 'XDir','reverse');
end

function Q = rot_mat(a,b,c)
    Q = [1, 0, 0; 0, cos(a), -sin(a); 0, sin(a), cos(a)];
    Q = Q * [cos(b), 0, sin(b); 0, 1, 0; -sin(b), 0, cos(b)];
    Q = Q * [cos(c), -sin(c), 0; sin(c), cos(c), 0; 0,0,1];
end