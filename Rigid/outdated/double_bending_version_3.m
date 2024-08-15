function double_bending()
    k = 2;
    h = 2^(-k);
    N = 2^(k+1);
    L = N^2;
    t11 = [0.25; -0.25; 1];
    t11 = 0.5 * h * t11 / norm(t11, 2);

    A = assemble_A(N,L,h,t11);
    %figure(1);
    %plot([0,0], [0,0]);
    spy(A);

    top_left = vertcat([-0.5; -0.5; 0], reshape(eye(3), [], 1));
    top_right = vertcat([-0.5; 0.5; 0], reshape(eye(3), [], 1));
    bottom_left = vertcat([0.5; -0.5; 0], reshape(eye(3), [], 1));
    bottom_right = vertcat([0.5; 0.5; 0], reshape(eye(3), [], 1));

    %top_left2 = [0.3000, -0.6611, 0.0710, 0.9363, 0.0705, -0.3441, ...
    %             -0.0653, 0.9975, 0.0267, 0.3451, -0.0025, 0.9386]';
    %top_right2 = [-0.1914, 0.0563, 0.1237, 0.3952, 0.8730, -0.2859, ...
    %              -0.9052, 0.4230, 0.0404, 0.1562, 0.2429, 0.9574]';
    %bottom_left2 = [0.9585, -0.1739, -0.0701, 0.3831, 0.9182, 0.1004, ...
    %                -0.9149, 0.3621, 0.1787, 0.1277, -0.1603, 0.9788]';
    %bottom_right2 = [0.5, 0.5, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1]';

    U0 = init_configuration(k);

    %dp = vertcat(1, (37:48)');
    dirichlet_point = indices(3,3,N);
    dirichlet_point = dirichlet_point(1);
    
    dp = vertcat(indices(2,2,N)', dirichlet_point);
    fp = setdiff(1:(12 * L), dp);
    nfp = size(fp, 1);
    %dpv = vertcat(0.3, bottom_right);
    dpv = vertcat(U0(indices(2,2,N)), U0(dirichlet_point)+0.2);

    A_lumped = A(:,fp);
    b = -A(:,dp) * dpv;

    %U0 = vertcat(top_left, top_right, bottom_left, bottom_right);

    %Y = init_configuration(2)
    %size(Y)
    %plot_box_skeleton(Y, 2, t11);

    % import casadi
    addpath('casadi')
    import casadi.*

    UU = SX.sym('x', 12 * L);

    alpha = 0;
    c = 1;
    for i = 1:N
        for j=1:N
                alpha = alpha + (UU((4:6) + 12 * (j-1) + 12 * N * (i-1))' * UU((4:6) + 12 * (j-1) + 12 * N * (i-1)) - 1)^2;
                alpha = alpha + c * (UU((4:6) + 12 * (j-1) + 12 * N * (i-1))' * UU((7:9) + 12 * (j-1) + 12 * N * (i-1)) - 0)^2;
                alpha = alpha + c * (UU((4:6) + 12 * (j-1) + 12 * N * (i-1))' * UU((10:12) + 12 * (j-1) + 12 * N * (i-1)) - 0)^2;
                alpha = alpha + (UU((7:9) + 12 * (j-1) + 12 * N * (i-1))' * UU((7:9) + 12 * (j-1) + 12 * N * (i-1)) - 1)^2;
                alpha = alpha + c * (UU((7:9) + 12 * (j-1) + 12 * N * (i-1))' * UU((10:12) + 12 * (j-1) + 12 * N * (i-1)) - 0)^2;
                alpha = alpha + (UU((10:12) + 12 * (j-1) + 12 * N * (i-1))' * UU((10:12) + 12 * (j-1) + 12 * N * (i-1)) - 1)^2;
        end
    end

    beta = SX.sym('y', L);
    for l = 1:L
        beta(l) = det33(reshape(UU((4:12) + 12 * (l-1)), 3, 3));
    end

    lbx = -ones(12 * L, 1) * inf;
    lbx(dp) = dpv; 
    ubx = ones(12 * L, 1) * inf;
    ubx(dp) = dpv; 

    % norm(1./ linspace(1,nfp,nfp)'
    nlp = struct('x', UU, 'f', 1 * norm(A_lumped * UU(fp) - b, 2)^2, 'g', vertcat(alpha, beta));
    S = nlpsol('S', 'ipopt', nlp);

    sol = S('x0', U0, 'lbg', 0, 'ubg', vertcat(0, inf * ones(L, 1)), 'lbx', lbx, 'ubx', ubx);

    U = sol.x.full();

    plot_box_skeleton(U, k, t11);
end


function A = assemble_A(N, L, h, t11)
    ctr = 1; 
    A = sparse(12 * N * (N-1), 12 * L);
    for i = 1:N
        for j = 1:N
            if i < N
                for k = 0:2
                    A(ctr+k, 1+12*(j-1)+(i-1) * 12 * N+k) = 1;
                    A(ctr+k, 1+12*(j-1)+(i-1) * 12 * N+3+k) = 0.5 * h;
                    A(ctr+k, 1+12*(j-1)+i * 12 * N+k) = -1;
                    A(ctr+k, 1+12*(j-1)+i * 12 * N+3+k) = 0.5 * h;
                end
                ctr = ctr + 3;
            end
            if j < N
                for k = 0:2
                    A(ctr+k, 1+12*(j-1)+(i-1) * 12 * N+k) = 1;
                    A(ctr+k, 1+12*(j-1)+(i-1) * 12 * N+6+k) = 0.5 * h;
                    A(ctr+k, 1+12*(j-1)+(i-1) * 12 * N+k+12) = -1;
                    A(ctr+k, 1+12*(j-1)+(i-1) * 12 * N+k+12+6) = 0.5 * h;
                end
                ctr = ctr + 3;
            end
            if i < N
                v = 0.5 * (t(j,i,t11) + t(j,j+1,t11))';
                for k = 0:2
                    A(ctr+k, 1+12*(j-1)+(i-1)* 12 * N+k+3) = v * [1;0;0];
                    A(ctr+k, 1+12*(j-1)+(i-1)* 12 * N+k+6) = v * [0;1;0];
                    A(ctr+k, 1+12*(j-1)+(i-1) * 12 * N+k+9) = v * [0;0;1];
                    A(ctr+k, 1+12*(j-1)+i*12*N+k+3) = -v * [1; 0; 0];
                    A(ctr+k, 1+12*(j-1)+i*12*N+k+6) = -v * [0; 1; 0];
                    A(ctr+k, 1+12*(j-1)+i*12*N+k+9) = -v * [0; 0; 1];
                end
                ctr = ctr + 3;
            end
            if j < N
                v = 0.5 * (t(j,i,t11) + t(j+1,i,t11))';
                for k = 0:2
                    A(ctr+k, 1+12*(j-1)+(i-1)*12*N+k+3) = v * [1;0;0];
                    A(ctr+k, 1+12*(j-1)+(i-1)*12*N+k+6) = v * [0;1;0];
                    A(ctr+k, 1+12*(j-1)+(i-1)*12*N+k+9) = v * [0;0;1];
                    A(ctr+k, 1+12*(j-1)+(i-1)*12*N+k+3+12) = -v * [1;0;0];
                    A(ctr+k, 1+12*(j-1)+(i-1)*12*N+k+6+12) = -v * [0;1;0];
                    A(ctr+k, 1+12*(j-1)+(i-1)*12*N+k+9+12) = -v * [0;0;1];
                end
                ctr = ctr + 3;
            end
        end
    end
end


function plot_line_skeleton(U, N, h, t11)
    for i=1:N
        for j=1:N
            plot3(U(1 + 12 * (j-1) + 12 * N * (i-1)), ...
                U(2 + 12 * (j-1) + 12 * N * (i-1)), ...
                U(3 + 12 * (j-1) + 12 * N * (i-1)), "bo"); hold on
            Q = reshape(U((4:12) + 12 * (j-1) + 12 * N * (i-1)), [3, 3]);
            x0 = U((1:3) + 12 * (j-1) + 12 * N * (i-1));
            north = Q * [-h/2; 0; 0];
            west = Q * [0; -h/2; 0];
            east = Q * [0; h/2; 0];
            south = Q * [h/2; 0; 0];
            plot3([x0(1), x0(1) + north(1)], [x0(2), x0(2) + north(2)], ...
                [x0(3), x0(3) + north(3)]);
            plot3([x0(1), x0(1) + south(1)], [x0(2), x0(2) + south(2)], ...
                [x0(3), x0(3) + south(3)]);
            plot3([x0(1), x0(1) + west(1)], [x0(2), x0(2) + west(2)], ...
                [x0(3), x0(3) + west(3)]);
            plot3([x0(1), x0(1) + east(1)], [x0(2), x0(2) + east(2)], ...
                [x0(3), x0(3) + east(3)]);
            tij = t(i,j,t11);
            plot3([x0(1), x0(1) + tij(1)], [x0(2), x0(2) + tij(2)], ...
                [x0(3), x0(3) + tij(3)]);     
        end
    end
    zlim([-1,1]);
end

function plot_box_skeleton(U, k, t11)
    N = 2^(k+1);
    L = N^2;
    h = 2^(-k);
    colors = ['r', 'b', 'g', 'm', 'c'];
    for i=1:N
        for j=1:N
            plot3(U(1 + 12 * (j-1) + 12 * N * (i-1)), ...
                U(2 + 12 * (j-1) + 12 * N * (i-1)), ...
                U(3 + 12 * (j-1) + 12 * N * (i-1)), "bo"); hold on
            Q = reshape(U((4:12) + 12 * (j-1) + 12 * N * (i-1)), [3, 3]);
            x0 = U((1:3) + 12 * (j-1) + 12 * N * (i-1));
            north = Q * [-h/2; 0; 0];
            west = Q * [0; -h/2; 0];
            east = Q * [0; h/2; 0];
            south = Q * [h/2; 0; 0];
            dirCross= [north, east, south, west];
            tij = t(j,i,t11);
            hNQ = Q * 0.5 * (t(j,i+1,t11) + tij);
            hSQ = Q * 0.5 * (t(j,i-1,t11) + tij);
            hWQ = Q * 0.5 * (t(j+1,i,t11) + tij);
            hEQ = Q * 0.5 * (t(j-1,i,t11) + tij);
            dirHinge = [hNQ, hEQ, hSQ, hWQ];
            tijQ = Q * tij;
            
            for l = 1:4
                plot3([x0(1) - tijQ(1), x0(1) + tijQ(1), ...
                x0(1) + dirCross(1,l) + dirHinge(1,l), ...
                x0(1) + dirCross(1,l) - dirHinge(1,l), ...
                x0(1) - tijQ(1)], ...
                [x0(2) - tijQ(2), x0(2) + tijQ(2), ...
                x0(2) + dirCross(2,l) + dirHinge(2,l), ...
                x0(2) + dirCross(2,l) - dirHinge(2,l), ...
                x0(2) - tijQ(2)], ...
                [x0(3) - tijQ(3), x0(3) + tijQ(3), ...
                x0(3) + dirCross(3,l) + dirHinge(3,l), ...
                x0(3) + dirCross(3,l) - dirHinge(3,l), ...
                x0(3) - tijQ(3)], colors(1 + mod(j + (i-1) * N, 5)));
            end
        end
    end
    set(gca, 'YDir','reverse')
    grid on
    zlim([-1,1]);

    %{
    for i=1:N
        for j=1:N
            plot3(U(1 + 12 * (j-1) + 12 * N * (i-1)), ...
                U(2 + 12 * (j-1) + 12 * N * (i-1)), ...
                U(3 + 12 * (j-1) + 12 * N * (i-1)), "bo"); hold on
            Q = reshape(U((4:12) + 12 * (j-1) + 12 * N * (i-1)), [3, 3]);
            x0 = U((1:3) + 12 * (j-1) + 12 * N * (i-1));
            north = Q * [-h/2; 0; 0];
            west = Q * [0; -h/2; 0];
            east = Q * [0; h/2; 0];
            south = Q * [h/2; 0; 0];
            tij = t(j,i,t11);
            hNQ = Q * 0.5 * (t(j,i+1,t11) + tij);
            hSQ = Q * 0.5 * (t(j,i-1,t11) + tij);
            hWQ = Q * 0.5 * (t(j+1,i,t11) + tij);
            hEQ = Q * 0.5 * (t(j-1,i,t11) + tij);
            tijQ = Q * tij;

            plot3([x0(1) - tijQ(1), x0(1) + tijQ(1), ...
                x0(1) + north(1) + hNQ(1), x0(1) + north(1) - hNQ(1), ...
                x0(1) - tijQ(1)], ...
                [x0(2) - tijQ(2), x0(2) + tijQ(2), ...
                x0(2) + north(2) + hNQ(2), x0(2) + north(2) - hNQ(2), ...
                x0(2) - tijQ(2)], ...
                [x0(3) - tijQ(3), x0(3) + tijQ(3), ...
                x0(3) + north(3) + hNQ(3), x0(3) + north(3) - hNQ(3), ...
                x0(3) - tijQ(3)])

            plot3([x0(1) - tijQ(1), x0(1) + tijQ(1), x0(1) + south(1) + hSQ(1), x0(1) + south(1) - hSQ(1), x0(1) - tijQ(1)], ...
                [x0(2) - tijQ(2), x0(2) + tijQ(2), x0(2) + south(2) + hSQ(2), x0(2) + south(2) - hSQ(2), x0(2) - tijQ(2)], ...
                [x0(3) - tijQ(3), x0(3) + tijQ(3), x0(3) + south(3) + hSQ(3), x0(3) + south(3) - hSQ(3), x0(3) - tijQ(3)])

            plot3([x0(1) - tijQ(1), x0(1) + tijQ(1), x0(1) + west(1) + hWQ(1), x0(1) + west(1) - hWQ(1), x0(1) - tijQ(1)], ...
                [x0(2) - tijQ(2), x0(2) + tijQ(2), x0(2) + west(2) + hWQ(2), x0(2) + west(2) - hWQ(2), x0(2) - tijQ(2)], ...
                [x0(3) - tijQ(3), x0(3) + tijQ(3), x0(3) + west(3) + hWQ(3), x0(3) + west(3) - hWQ(3), x0(3) - tijQ(3)])

            plot3([x0(1) - tijQ(1), x0(1) + tijQ(1), x0(1) + east(1) + hEQ(1), x0(1) + east(1) - hEQ(1), x0(1) - tijQ(1)], ...
                [x0(2) - tijQ(2), x0(2) + tijQ(2), x0(2) + east(2) + hEQ(2), x0(2) + east(2) - hEQ(2), x0(2) - tijQ(2)], ...
                [x0(3) - tijQ(3), x0(3) + tijQ(3), x0(3) + east(3) + hEQ(3), x0(3) + east(3) - hEQ(3), x0(3) - tijQ(3)])
        end
    end
    set(gca, 'YDir','reverse')
    grid on
    zlim([-1,1]);
    %}
end


function crease_vectors()
    k = 3;
    h = 2^(-k);
    N = 2^(k+1);
    L = N^2;

    %U = zeros(12 * L, 1);
    t11 = [0.5; 0.5; 1];

    x = linspace(1, N, N);

    [x,y] = meshgrid(x,x);

    x = reshape(x, [], 1);
    y = reshape(y, [], 1);

    UU = zeros(N, N);
    VV = zeros(N, N);
    WW = zeros(N, N);

    for i=1:size(x, 1)
            v = t(x(i), y(i), t11);
            UU(x(i), y(i)) = v(1);
            VV(x(i), y(i)) = v(2);
            WW(x(i), y(i)) = 0;
    end

    xx = linspace(1-h, -1, N)+h/2;

    [XX, YY] = meshgrid(xx, xx);
    ZZ = zeros(size(XX));

    quiver(XX, YY, UU, VV, 0.25);
    xlim([-1 1]);
    ylim([-1 1]);
end

function U0 = init_configuration(k)
    N = 2^(k+1);
    h = 2^(-k);
    L = N^2;
    U0 = zeros(12 * L, 1);
    for i = 1:N
        for j = 1:N
            U0((1:2) + 12 * (j-1) + 12 * N * (i-1)) = [-1 + (i-1)*h + h/2, -1 + (j-1)*h + h/2];
            U0((4:12) + 12 * (j-1) + 12 * N * (i-1)) = reshape(eye(3), 3, 3);
        end
    end
end

function tij = t(i,j, t11)
    tij = [t11(1) * (-1)^(j-1); t11(2)*(-1)^(i-1); t11(3)];
end

function v = det33(M)
    v = M(1,1) * M(2,2) * M(3,3) + M(1, 2) * M(2, 3) * M(3, 1) + ...
        M(1,3) * M(2,1) * M(3,2) - M(1,3) * M(2,2) * M(3,1) - ...
        M(1,2) * M(2,1) * M(3,3) - M(1,1) * M(2,3) * M(3,2);
end

function Q = rot_mat(a,b,c)
    Q = [1, 0, 0; 0, cos(a), -sin(a); 0, sin(a), cos(a)];
    Q = Q * [cos(b), 0, sin(b); 0, 1, 0; -sin(b), 0, cos(b)];
    Q = Q * [cos(c), -sin(c), 0; sin(c), cos(c), 0; 0,0,1];
end

function v = indices(i, j, N)
    v = ((1:12) + 12 * (j - 1) + 12 * N * (i - 1));
end