function exact_solution_old()
    k = 0;
    h = 2^(-k);
    N = 2^(k+1);
    L = N^2;
    %t11 = [0.00; -0.00; 1];
    t11 = [0.25; -0.25; 1];
    t11 = 0.5 * h * t11 / norm(t11, 2);
    
    p1 = [-0.5; -0.5; 0];
    Q1 = eye(3);
    p2 = [-0.5; 0.5; 0];
    Q2 = eye(3);
    p3 = [0.5; -0.5; 0];
    Q3 = eye(3);
    p4 = [0.5; 0.5; 0];
    Q4 = eye(3);

    v = 0.5 * (t(1,1,t11) + t(2,1,t11));

    % rotates the hinge v up
    QQ = rot_mat(0,atan(-v(1)/v(3)),0);

    theta_blue_green = 0.7;
    p2_new = inv(QQ) * rot_mat(0,0,theta_blue_green) * (p2 - [-0.5;0;0]) + [-0.5;0;0];
    Q2_new = inv(QQ) * rot_mat(0,0,theta_blue_green) * QQ * Q2;

    v = 0.5 * (t(1,1,t11) + t(1,2,t11));
    QQQ = rot_mat(atan(v(2)/v(3)),0,0);
    
    %theta_blue_pink = 0.4;
    theta_blue_pink = 0.6860148;

    p3_new = inv(QQQ) * rot_mat(0,0,theta_blue_pink) * (p3 - [0;-0.5;0]) + [0;-0.5;0];
    Q3_new = inv(QQQ) * rot_mat(0,0,theta_blue_pink) * QQQ * Q3;

    theta_pink_cyan = 0.7;

    p4_new = inv(QQQ) * rot_mat(0,0,theta_blue_pink) * QQQ * (p4 - [0;-0.5;0]) + [0;-0.5;0];
    Q4_new = inv(QQQ) * rot_mat(0,0,theta_blue_pink) * QQQ * Q4;

    qhinge = Q4_new * 0.5 * (t(2,2,t11) + t(1,2,t11));

    SS = rot_mat(0,atan(-qhinge(1)/qhinge(3)),0);
    qhingeS = SS * qhinge;
    SSS = rot_mat(atan(qhingeS(2)/qhingeS(3)),0,0);
    qhingeSS = SSS * qhingeS;
    S = SSS * SS;

    mp_new = 0.5 * (p3_new + p4_new);

    p4_new_new = inv(S) * rot_mat(0,0,-theta_pink_cyan) * S * (p4_new - mp_new) + mp_new;
    Q4_new_new = inv(S) * rot_mat(0,0,-theta_pink_cyan) * S * Q4_new;

    U = vertcat(p1, reshape(Q1, [], 1), p2_new, reshape(Q2_new, [], 1), ...
        p3_new, reshape(Q3_new, [], 1), p4_new_new, reshape(Q4_new_new, [], 1));

    plot_box_skeleton(U, k, t11);
    %line([p4_new_new(1), p4_new_new(1)+1], [p4_new_new(2), p4_new_new(2)+1], ...
    %    [p4_new_new(3), p4_new_new(3)+1]);

    A = assemble_A(N, L, h, t11);
    norm(A * U, 2)
end

function tij = t(i,j, t11)
    tij = [t11(1) * (-1)^(j-1); t11(2)*(-1)^(i-1); t11(3)];
end