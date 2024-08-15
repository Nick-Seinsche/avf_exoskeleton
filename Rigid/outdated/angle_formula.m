function angle_formula()
    k = 0;
    h = 2^(-k);
    N = 2^(k+1);
    L = N^2;
    t11 = [0.87; -0.87; 1];
    t11 = 0.5 * h * t11 / norm(t11, 2);

    a = 36.9811 / 360 * (2 * pi);
    b = 31.3686 / 360 * (2 * pi);

    tau = 1 / sqrt(1 + (t11(1) / t11(3))^2);
    tau_tilde = t11(1) / t11(3) * tau;

    j = 1 / sqrt(1 + (t11(2) / t11(3))^2);
    j_tilde = -t11(2) / t11(3) * j;

    Q2 = [tau^2 * cos(a) + tau_tilde^2, -tau*sin(a), tau * tau_tilde * (cos(a) - 1);
        tau * sin(a), cos(a), tau_tilde * sin(a);
        tau * tau_tilde * (cos(a) - 1), - tau_tilde * sin(a), tau_tilde^2 * cos(a) + tau^2];

    Q3 = [cos(b), -j*sin(b), j_tilde*sin(b);
        j*sin(b), j^2*cos(b) + j_tilde^2, j*j_tilde * (1-cos(b));
        -j_tilde*sin(b), j*j_tilde*(1-cos(b)), j_tilde^2*cos(b)+j^2];

    qhinge34 = [t11(1) * cos(b) + j_tilde * t11(3) * sin(b); j * t11(1) * sin(b) + j * j_tilde * t11(3) * (1 - cos(b)); -j_tilde * t11(1) * sin(b) + j_tilde^2 * t11(3) * cos(b) + j^2 * t11(3)];
    
    s = sqrt(1 + ((  t11(1) * cos(b) + j_tilde * t11(3) * sin(b))/(-j_tilde * t11(1) * sin(b) + j_tilde^2 * t11(3) * cos(b) + j^2 * t11(3)))^2);
    qqhinge34 = [0; j * t11(1) * sin(b) + j * j_tilde * t11(3) * (1 - cos(b)); (-j_tilde * t11(1) * sin(b) + j_tilde^2 * t11(3) * cos(b) + j^2 * t11(3) + (t11(1) * cos(b) + j_tilde * t11(3) * sin(b))^2/(-j_tilde * t11(1) * sin(b) + j_tilde^2 * t11(3) * cos(b) + j^2 * t11(3)))/s];


    w = Q3 * [t11(1), 0, t11(3)]';
    i = 1 / sqrt(1 + (w(1) / w(3))^2);
    i_tilde = -(w(1) / w(3)) * i;
    G = [i, 0, i_tilde; 0, 1, 0; -i_tilde, 0, i];
    v = G * w;
    k = 1 / sqrt(1 + (v(2) / v(3))^2);
    k_tilde = (v(2) / v(3)) * k;
    H = [1,0,0;0,k,-k_tilde;0,k_tilde,k];
    D = H * G;

    U = init_configuration(0);

    Q4 = inv(D) * rot_mat(0,0,-a) * D * Q3;

    t1 = t11(1);
    t2 = t11(2);
    t3 = t11(3);
    
    q = t1^2/t3^2;

    %for ell=0:5:45

    %a = 36.9811 / 360 * (2 * pi);
    %a = 19.5215 / 360 * (2 * pi);
    a = 15 / 360 * (2 * pi);
    
    f = @(b) abs(t1/(2*t3*(t1^2/t3^2 + 1)) - t2/(2*t3*(t2^2/t3^2 + 1)) - (t1*cos(a))/(2*t3*(t1^2/t3^2 + 1)) + (t2*cos(b))/(2*t3*(t2^2/t3^2 + 1)) + (t1*sin(a))/(2*t3*(t1^2/t3^2 + 1)^(1/2)) + (t2*sin(b))/(2*t3*(t2^2/t3^2 + 1)^(1/2)))^2 + abs(cos(b)/2 - cos(a)/(2*(t1^2/t3^2 + 1)) + sin(a)/(2*(t1^2/t3^2 + 1)^(1/2)) - sin(b)/(2*(t2^2/t3^2 + 1)^(1/2)) - t1^2/(2*t3^2*(t1^2/t3^2 + 1)) + 1/2)^2 + abs(cos(a)/2 - cos(b)/(2*(t2^2/t3^2 + 1)) + sin(a)/(2*(t1^2/t3^2 + 1)^(1/2)) - sin(b)/(2*(t2^2/t3^2 + 1)^(1/2)) - t2^2/(2*t3^2*(t2^2/t3^2 + 1)) + 1/2)^2 - 1/2;
    g = @(b) q * ((2 - cos(b) - cos(a))/(q+1) + (sin(a) - sin(b))/(sqrt(q + 1)))^2 + (cos(b) + (1 - cos(a))/(q + 1)+(sin(a) - sin(b))/(sqrt(1 + q)))^2 + (cos(a) + (1 - cos(b))/(q + 1)+(sin(a) - sin(b))/(sqrt(1 + q)))^2 - 2;


    omega([0.3, 0.4], 0, 0)
    omega([-0.4, -0.3], 0, 0)
    [X,Y] = meshgrid(-pi/2:0.1:pi/2, -pi/2:0.1:pi/2);
    Z = zeros(size(X));
    for i = 1:size(X,1)
        for j = 1:size(X, 2)
            Z(i,j) = omega([X(i,j), Y(i,j)], 0, 0);
        end
    end

    surfc(X,Y,Z);
    xlabel("xaxis")
    ylabel("yaxis")
end

function draw_vector(start, direction)
    line([start(2), start(2) + direction(2)], [start(1), start(1) + direction(1)], ...
        [start(3), start(3) + direction(3)], "color", "red");
end

function Q = rot_mat(a,b,c)
    Q = [1, 0, 0; 0, cos(a), -sin(a); 0, sin(a), cos(a)];
    Q = Q * [cos(b), 0, sin(b); 0, 1, 0; -sin(b), 0, cos(b)];
    Q = Q * [cos(c), -sin(c), 0; sin(c), cos(c), 0; 0,0,1];
end