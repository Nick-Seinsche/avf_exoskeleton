function test_transform3()
    h = 1/2;
    h_T = 1/4;
    alpha = 20.9;
    alpha_rad = alpha / 360 * 2 * pi;
    wall_incline = tan(alpha_rad);
    %t11 = [wall_incline; -wall_incline; 1];
    t11 = [-tan(alpha_rad); tan(alpha_rad); 1];

    t11_v = rot_mat(atan(t11(2)),0,0) * t11;

    z1 = [h_T/2, -t11_v(1)/t11_v(3) * h_T/2];
    z4 = -z1;
    z3 = [h_T/2, h/2];
    z6 = [-h_T/2, h/2];
    z2 = 0.5 * (z3 + z1);
    z5 = 0.5 * (z4 + z6);

    z = [z1', z2', z3', z4', z5', z6'];

    x = zeros(3, 1);
    Q = eye(3);

    QdirN = [0, -1; 0, 0; 1, 0];
    QdirS = [0, 1; 0, 0; 1, 0];
    QdirW = [0, 0; 0, -1; -1, 0];
    QdirE = [0, 0; 0, 1; -1, 0];

    plateN = rot_mat(atan(t11(1)),0,0) * QdirN * z;
    
    plateS = rot_mat(atan(t11(2)),0,0) * QdirS * z;
    
    plateW = rot_mat(0,atan(t11(1)),0) * QdirW * z;
    
    plateE = rot_mat(0,atan(t11(2)),0) * QdirE * z;
    
    for i=1:6
        plot3(plateN(2,i), plateN(1,i), plateN(3,i), "bo"); hold on
    end
    fc =fill3(plateN(2,[1,2,3,6,5,4]), plateN(1,[1,2,3,6,5,4]), plateN(3,[1,2,3,6,5,4]), "b");
    fc(1).FaceAlpha = 0.25;

    for i=1:6
        plot3(plateS(2,i), plateS(1,i), plateS(3,i), "ro"); hold on
    end
    fc =fill3(plateS(2,[1,2,3,6,5,4]), plateS(1,[1,2,3,6,5,4]), plateS(3,[1,2,3,6,5,4]), "r");
    fc(1).FaceAlpha = 0.25;

    for i=1:6
        plot3(plateW(2,i), plateW(1,i), plateW(3,i), "go"); hold on
    end
    fc =fill3(plateW(2,[1,2,3,6,5,4]), plateW(1,[1,2,3,6,5,4]), plateW(3,[1,2,3,6,5,4]), "g");
    fc(1).FaceAlpha = 0.25;

    for i=1:6
        plot3(plateE(2,i), plateE(1,i), plateE(3,i), "mo"); hold on
    end
    fc =fill3(plateE(2,[1,2,3,6,5,4]), plateE(1,[1,2,3,6,5,4]), plateE(3,[1,2,3,6,5,4]), "m");
    fc(1).FaceAlpha = 0.25;

    draw_vector([0,0,0], 0.3 * t11);

    %draw_vector([0,0,0], 0.3 * rot_mat(-alpha_rad,0,0) * t11);
    %draw_vector([0,0,0], 0.3 * [0, t11(2), t11(3)]);
    %draw_vector([0,0,0], 0.3 * [t11(1), 0, t11(3)]);

    %[t11(1), 0, t11(3)];
    %rot_mat(0,alpha_rad,0) * [0, t11(2), t11(3)]';


    t22 = rot_mat(-atan(t11(1)/t11(3)),0,0) * t11;
    rot_mat(0, -atan(t22(1)/t22(3)),0) * t22;

    rot_mat(-atan(t11(1)/t11(3)),0,0) * rot_mat(0, -atan(t22(1)/t22(3)),0);
    rot_mat(0, -atan(t22(1)/t22(3)),0) * rot_mat(-atan(t11(1)/t11(3)),0,0);

    xlabel('x-axis');
    ylabel('y-axis');
    zlabel('z-axis');
end

function Q = rot_mat(a,b,c)
    Q = [1, 0, 0; 0, cos(a), -sin(a); 0, sin(a), cos(a)];
    Q = Q * [cos(b), 0, sin(b); 0, 1, 0; -sin(b), 0, cos(b)];
    Q = Q * [cos(c), -sin(c), 0; sin(c), cos(c), 0; 0,0,1];
end

function draw_vector(start, direction)
    line([start(2), start(2) + direction(2)], [start(1), start(1) + direction(1)], ...
        [start(3), start(3) + direction(3)], "color", "green");
end