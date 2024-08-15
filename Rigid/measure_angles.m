function measure_angles()
    k = 0;
    h_T = 1;
    dx = 0.3;
    alpha = 21.8; % 21.8;

    U = bidirectional_bending_nonlinear(k, h_T, dx, alpha, 1); %k, delta_x, crease_vector

    % unpack the solution
    p1 = U(1:3);
    Q1 = reshape(U(4:12), 3, 3);
    p2 = U(13:15);
    Q2 = reshape(U(16:24), 3, 3);
    p3 = U(25:27);
    Q3 = reshape(U(28:36), 3, 3);
    p4 = U(37:39);
    Q4 = reshape(U(40:48), 3, 3);

    % calculate the hinge angles
    %angle_blue_green = angle((Q1 * [0;1;0]), -(Q2 * [0;-1;0])) / (2 * pi) * 360
    %angle_blue_pink = angle((Q1 * [1;0;0]), -(Q3 * [-1;0;0])) / (2 * pi) * 360
    %angle_green_cyan = angle((Q2 * [1;0;0]), -(Q4 * [-1;0;0])) / (2 * pi) * 360
    %angle_cyan_pink = angle((Q3 * [0;1;0]), -(Q4 * [0;-1;0])) / (2 * pi) * 360

    angle_blue_green = angle((Q1 * [0;1;0]), -(Q2 * [0;-1;0]))
    angle_blue_pink = angle((Q1 * [1;0;0]), -(Q3 * [-1;0;0]))
    angle_green_cyan = angle((Q2 * [1;0;0]), -(Q4 * [-1;0;0]))
    angle_cyan_pink = angle((Q3 * [0;1;0]), -(Q4 * [0;-1;0]))

    omega(angle_blue_green, 1, 0.2)

    %h = 2^(-0);
    %Q1 * [0, h/2, 0]' + Q2 * [0, h/2, 0]' + Q2 * [h/2, 0, 0]' + Q4 * [h/2, 0, 0]' + Q4 * [0, -h/2, 0]' + Q3 * [0, -h/2, 0]' + Q3 * [-h/2,0, 0]' + Q1 * [-h/2, 0, 0]';
    %Q1(:,2) + Q2(:,2) + Q2(:,1) + Q4(:,1) - Q4(:,2) - Q3(:,2) - Q3(:,1) - Q1(:,1);
end

function draw_vector(start, direction)
    line([start(1), start(1) + direction(1)], [start(2), start(2) + direction(2)], ...
        [start(3), start(3) + direction(3)], "color", "red");
end

function alpha = angle(v,w)
    alpha = acos(dot(v,w) / (norm(v, 2) * norm(w, 2)));
end