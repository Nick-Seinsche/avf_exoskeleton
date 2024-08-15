function double_bending()
    k = 0;
    h = 2^(-k);
    N = 2^(k+1);
    L = N^2;
    t11 = 0.4 * [-0.25; -0.25; 1];

    A = assemble_A(N,L,h,t11);
    %figure(1);
    %plot([0,0], [0,0]);
    %spy(A);
    
    dp = vertcat(1, (37:48)');
    fp = setdiff(1:48, dp);

    top_left = vertcat([-0.5; -0.5; 0], reshape(eye(3), [], 1));
    top_right = vertcat([-0.5; 0.5; 0], reshape(eye(3), [], 1));
    bottom_left = vertcat([0.5; -0.5; 0], reshape(eye(3), [], 1));
    bottom_right = vertcat([0.5; 0.5; 0], reshape(eye(3), [], 1));

    dpv = vertcat(-0.4, bottom_right);

    %U = vertcat(top_left, top_right, bottom_left, bottom_right);

    A_lumped = A(:,fp);

    b = -A(:,dp) * dpv;

    %U_fp = A_lumped \ b;

    %U = zeros(12 * L, 1);
    %U(fp) = U_fp;
    %U(dp) = vertcat(top_left, bottom_right);
    
    %A_lumped * vertcat(top_right, bottom_left)

    %-A(:,dp) * vertcat(top_left, bottom_right)

    U0 = vertcat(top_left, top_right, bottom_left, bottom_right);

    % import casadi
    addpath('casadi linux')
    import casadi.*

    UU = SX.sym('x', 12 * L);

    %UU(fp) = UU_fp;

    %lbg = ones(12 * L) * inf;

    %UU = zeros(12 * L, 1);
    %UU(fp) = UU_fp;
    %UU(dp) = dpv;
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

    lbx = -ones(12 * L, 1) * inf;
    lbx(dp) = dpv; 
    ubx = ones(12 * L, 1) * inf;
    ubx(dp) = dpv; 

    nlp = struct('x', UU, 'f', 10 * alpha, 'g', A_lumped * UU(fp) - b);
    S = nlpsol('S', 'ipopt', nlp);

    sol = S('x0', U0, 'lbg', 0, 'ubg', 0, 'lbx', lbx, 'ubx', ubx);

    U = sol.x.full();

    A * U

    %U = U0

    %disp(U);

    %figure(2);
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

            %plot3([x0(1), x0(1) + north(1)], [x0(2), x0(2) + north(2)], ...
            %    [x0(3), x0(3) + north(3)]);
            plot3([x0(1) - tijQ(1), x0(1) + tijQ(1), x0(1) + north(1) + hNQ(1), x0(1) + north(1) - hNQ(1), x0(1) - tijQ(1)], ...
                [x0(2) - tijQ(2), x0(2) + tijQ(2), x0(2) + north(2) + hNQ(2), x0(2) + north(2) - hNQ(2), x0(2) - tijQ(2)], ...
                [x0(3) - tijQ(3), x0(3) + tijQ(3), x0(3) + north(3) + hNQ(3), x0(3) + north(3) - hNQ(3), x0(3) - tijQ(3)])

            plot3([x0(1) - tijQ(1), x0(1) + tijQ(1), x0(1) + south(1) + hSQ(1), x0(1) + south(1) - hSQ(1), x0(1) - tijQ(1)], ...
                [x0(2) - tijQ(2), x0(2) + tijQ(2), x0(2) + south(2) + hSQ(2), x0(2) + south(2) - hSQ(2), x0(2) - tijQ(2)], ...
                [x0(3) - tijQ(3), x0(3) + tijQ(3), x0(3) + south(3) + hSQ(3), x0(3) + south(3) - hSQ(3), x0(3) - tijQ(3)])

            plot3([x0(1) - tijQ(1), x0(1) + tijQ(1), x0(1) + west(1) + hWQ(1), x0(1) + west(1) - hWQ(1), x0(1) - tijQ(1)], ...
                [x0(2) - tijQ(2), x0(2) + tijQ(2), x0(2) + west(2) + hWQ(2), x0(2) + west(2) - hWQ(2), x0(2) - tijQ(2)], ...
                [x0(3) - tijQ(3), x0(3) + tijQ(3), x0(3) + west(3) + hWQ(3), x0(3) + west(3) - hWQ(3), x0(3) - tijQ(3)])

            plot3([x0(1) - tijQ(1), x0(1) + tijQ(1), x0(1) + east(1) + hEQ(1), x0(1) + east(1) - hEQ(1), x0(1) - tijQ(1)], ...
                [x0(2) - tijQ(2), x0(2) + tijQ(2), x0(2) + east(2) + hEQ(2), x0(2) + east(2) - hEQ(2), x0(2) - tijQ(2)], ...
                [x0(3) - tijQ(3), x0(3) + tijQ(3), x0(3) + east(3) + hEQ(3), x0(3) + east(3) - hEQ(3), x0(3) - tijQ(3)])

            %plot3([x0(1), x0(1) + south(1)], [x0(2), x0(2) + south(2)], ...
            %    [x0(3), x0(3) + south(3)]);
            %plot3([x0(1), x0(1) + west(1)], [x0(2), x0(2) + west(2)], ...
            %    [x0(3), x0(3) + west(3)]);
            %plot3([x0(1), x0(1) + east(1)], [x0(2), x0(2) + east(2)], ...
            %    [x0(3), x0(3) + east(3)]);
            %plot3([x0(1), x0(1) + tijQ(1)], [x0(2), x0(2) + tijQ(2)], ...
            %    [x0(3), x0(3) + tijQ(3)]);     

            %[x0 - tij, x0+tij, x0 + north +tij, x0 + north - tij]
        end
    end
    set(gca, 'YDir','reverse')
    zlim([-1,1]);
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


function plot_skeleton(U, N, h, t11)
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


function tij = t(i,j, t11)
    tij = [t11(1) * (-1)^(j-1); t11(2)*(-1)^(i-1); t11(3)];
end