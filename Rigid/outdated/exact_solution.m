function exact_solution()
    k = 0;
    h = 2^(-k);
    N = 2^(k+1);
    L = N^2;
    t11 = [-0.87; 0.87; 1];
    t11 = 0.5 * h * t11 / norm(t11, 2);

    blue_green = 15.9811 / 360 * (2 * pi);
    blue_pink = omega(blue_green, 1, 0);

    U0 = init_configuration(k);

    p1 = U0(1:3);
    Q1 = reshape(U0(4:12), 3, 3);
    p2 = U0(13:15);
    Q2 = reshape(U0(16:24), 3, 3);
    p3 = U0(25:27);
    Q3 = reshape(U0(28:36), 3, 3);
    p4 = U0(37:39);
    Q4 = reshape(U0(40:48), 3, 3);
    
    mp12 = 0.5 * (p1 + p2);
    mp13 = 0.5 * (p1 + p3);

    Q = rot_mat(0,atan(t11(1)/t11(3)),0);
    QQ = rot_mat(atan(-t11(2)/t11(3)),0,0);

    p2_new = Q' * rot_mat(0,0,blue_green) * Q * (p2 - mp12) + mp12;
    Q2_new = Q' * rot_mat(0,0,blue_green) * Q * Q2;

    p3_new = QQ' * rot_mat(0,0,blue_pink) * QQ * (p3 - mp13) + mp13;
    Q3_new = QQ' * rot_mat(0,0,blue_pink)* QQ * Q3;

    p4_new = QQ' * rot_mat(0,0,blue_pink) * QQ * (p4 - mp13) + mp13;
    Q4_new = QQ' * rot_mat(0,0,blue_pink) * QQ * Q4;

    %p4_new_alt = Q' * rot_mat(0,0,blue_green) * Q * (p4 - mp12) + mp12;
    %Q4_new_alt = Q' * rot_mat(0,0,blue_green) * Q * Q4;

    qhinge34 = QQ' * rot_mat(0,0,blue_pink)* QQ * [t11(1); 0; t11(3)];
    mp34 = 0.5 * (p3_new + p4_new);
    %mp24 = 0.5 * (p2_new + p4_new_alt);

    qhinge24 = Q' * rot_mat(0,0,blue_green) * Q * [0; t11(2); t11(3)];
   
    qqhinge34 = rot_mat(0,-atan(qhinge34(1)/qhinge34(3)),0) * qhinge34;
    %qqqhinge34 = rot_mat(atan(qqhinge34(2)/qqhinge34(3)),0,0) * qqhinge34;

    %qqhinge24 = rot_mat(0,-atan(qhinge24(1)/qhinge24(3)),0) * qhinge24;

    % rotates the hinge vector parallel to the x3 axis
    QQQ = rot_mat(atan(qqhinge34(2)/qqhinge34(3)),0,0) * rot_mat(0,-atan(qhinge34(1)/qhinge34(3)),0);

    %QQQ_alt = rot_mat(atan(qqhinge24(2)/qqhinge24(3)),0,0) * rot_mat(0,-atan(qhinge24(1)/qhinge24(3)),0);

    p4_new_new = QQQ' * rot_mat(0,0,-blue_green) * QQQ * (p4_new - mp34) + mp34;
    Q4_new_new = QQQ' * rot_mat(0,0,-blue_green) * QQQ * Q4_new;

    %p4_new_new_alt = inv(QQQ_alt) * rot_mat(0,0,-blue_pink) * QQQ_alt * (p4_new_alt - mp24) + mp24;
    %Q4_new_new_alt = inv(QQQ_alt) * rot_mat(0,0,-blue_pink) * QQQ_alt * Q4_new_alt;

    U = vertcat(p1, ...
                reshape(Q1, [], 1), ...
                ...
                p2_new, ...
                reshape(Q2_new, [], 1), ...
                ...
                p3_new, ...
                reshape(Q3_new, [], 1), ...
                ...
                p4_new_new, ...
                reshape(Q4_new_new, [], 1));


    plot_skeleton_square(U, k, t11, 1, 1);
    draw_vector(p1, p2_new - p1)
    mp12n = 0.5 * (p1 + p2_new);
    aa = p1
    p1
    p1
end

%draw_vector(p1, Q1 * [0, h/2, 0]');
    %draw_vector(p1 + Q1 * [0, h/2, 0]', Q2_new * [0, h/2, 0]');
    %draw_vector(p1 + Q1 * [0, h/2, 0]' + Q2_new * [0, h/2, 0]', Q2_new * [h/2, 0, 0]');
    %draw_vector(p1 + Q1 * [0, h/2, 0]' + Q2_new * [0, h/2, 0]' + Q2_new * [h/2, 0, 0]', Q4_new_new * [h/2, 0, 0]');
    %draw_vector(p1 + Q1 * [0, h/2, 0]' + Q2_new * [0, h/2, 0]' + Q2_new * [h/2, 0, 0]' + Q4_new_new * [h/2, 0, 0]', Q4_new_new * [0, -h/2, 0]');
    %draw_vector(p1 + Q1 * [0, h/2, 0]' + Q2_new * [0, h/2, 0]' + Q2_new * [h/2, 0, 0]' + Q4_new_new * [h/2, 0, 0]' + Q4_new_new * [0, -h/2, 0]', Q3_new * [0, -h/2, 0]');
    %draw_vector(p1 + Q1 * [0, h/2, 0]' + Q2_new * [0, h/2, 0]' + Q2_new * [h/2, 0, 0]' + Q4_new_new * [h/2, 0, 0]' + Q4_new_new * [0, -h/2, 0]' + Q3_new * [0, -h/2, 0]', Q3_new * [-h/2,0, 0]');
    %draw_vector(p1 + Q1 * [0, h/2, 0]' + Q2_new * [0, h/2, 0]' + Q2_new * [h/2, 0, 0]' + Q4_new_new * [h/2, 0, 0]' + Q4_new_new * [0, -h/2, 0]' + Q3_new * [0, -h/2, 0]' + Q3_new * [-h/2,0, 0]', Q1 * [-h/2, 0, 0]');

    %draw_vector(mp24, qhinge24);
    %draw_vector(mp34, qhinge34);
    %draw_vector(mp34, mp24 - mp34);

    %Q1 * [0, h/2, 0]' + Q2_new * [0, h/2, 0]' + Q2_new * [h/2, 0, 0]' + Q4_new_new * [h/2, 0, 0]' + Q4_new_new * [0, -h/2, 0]' + Q3_new * [0, -h/2, 0]' + Q3_new * [-h/2,0, 0]' + Q1 * [-h/2, 0, 0]';

    %A = assemble_A(N, L, h, t11);
    %norm(A * U, 2);

function draw_vector(start, direction)
    line([start(2), start(2) + direction(2)], [start(1), start(1) + direction(1)], ...
        [start(3), start(3) + direction(3)], "color", "red");
end

function draw_vector2(start, direction)
    line([start(2), start(2) + direction(2)], [start(1), start(1) + direction(1)], ...
        [start(3), start(3) + direction(3)], "color", "green");
end

function Q = rot_mat(a,b,c)
    Q = [1, 0, 0; 0, cos(a), -sin(a); 0, sin(a), cos(a)];
    Q = Q * [cos(b), 0, sin(b); 0, 1, 0; -sin(b), 0, cos(b)];
    Q = Q * [cos(c), -sin(c), 0; sin(c), cos(c), 0; 0,0,1];
end