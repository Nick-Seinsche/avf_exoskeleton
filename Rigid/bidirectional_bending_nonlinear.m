function [U, objective] = bidirectional_bending_nonlinear(k, h_T, dx, alpha, do_plot)
    % k - determines the size of the configuration
    % h_T - determines the thickness of the configuration
    % dx - determines the perturbation for the boundary condition
    % alpha - determines the wall inclination
    % do_plot - true/false
    h = 2^(-k);
    N = 2^(k+1);
    L = N^2;
    wall_incline = h_T * tan(alpha / 360 * 2 * pi) / 2;
    t11 = [wall_incline; -wall_incline; h_T / 2];

    A = assemble_A(k,t11);
    figure(1);
    spy(A);

    U0 = init_configuration(k);

    fixed_point = indices(2,2,N);
    fixed_point = fixed_point(1);
    
    dp = vertcat(indices(1,1,N)', fixed_point);
    fp = setdiff(1:(12 * L), dp);
    nfp = size(fp, 1);
    dpv = vertcat(U0(indices(1,1,N)), U0(fixed_point)-dx);

    A_lumped = A(:,fp);
    b = -A(:,dp) * dpv;

    % import casadi
    addpath('../casadi_windows')
    import casadi.*

    UU = SX.sym('x', 12 * L);

    alpha = SX.sym('alpha', 6 * L);
    for i = 1:N
        for j=1:N
            alpha(1+6*(j-1)+6*N*(i-1)) = (UU((4:6) + 12 * (j-1) + 12 * N * (i-1))' * UU((4:6) + 12 * (j-1) + 12 * N * (i-1)) - 1);
            alpha(2+6*(j-1)+6*N*(i-1)) = (UU((4:6) + 12 * (j-1) + 12 * N * (i-1))' * UU((7:9) + 12 * (j-1) + 12 * N * (i-1)) - 0);
            alpha(3+6*(j-1)+6*N*(i-1)) = (UU((4:6) + 12 * (j-1) + 12 * N * (i-1))' * UU((10:12) + 12 * (j-1) + 12 * N * (i-1)) - 0);
            alpha(4+6*(j-1)+6*N*(i-1)) = (UU((7:9) + 12 * (j-1) + 12 * N * (i-1))' * UU((7:9) + 12 * (j-1) + 12 * N * (i-1)) - 1);
            alpha(5+6*(j-1)+6*N*(i-1)) = (UU((7:9) + 12 * (j-1) + 12 * N * (i-1))' * UU((10:12) + 12 * (j-1) + 12 * N * (i-1)) - 0);
            alpha(6+6*(j-1)+6*N*(i-1)) = (UU((10:12) + 12 * (j-1) + 12 * N * (i-1))' * UU((10:12) + 12 * (j-1) + 12 * N * (i-1)) - 1);
        end
    end

    beta = SX.sym('beta', L);
    for l = 1:L
        beta(l) = det33(reshape(UU((4:12) + 12 * (l-1)), 3, 3));
    end
    
    lbx = -ones(12 * L, 1) * inf;
    lbx(dp) = dpv; 
    ubx = ones(12 * L, 1) * inf;
    ubx(dp) = dpv; 
    
    nlp = struct('x', UU, 'f', 1 * norm((A_lumped * UU(fp) - b), 2)^2, 'g', vertcat(alpha, beta));
    S = nlpsol('S', 'ipopt', nlp);

    sol = S('x0', U0, 'lbg', vertcat(zeros(6 * L, 1), 0 * ones(L, 1)), ...
                      'ubg', vertcat(zeros(6 * L, 1), inf * ones(L, 1)), ...
                      'lbx', lbx, 'ubx', ubx);

    U = sol.x.full();
    objective = sol.f.full();

    if do_plot
        figure(2);
        plot_skeleton_square(U, k, t11, 1, 1);
    end
end

function v = det33(M)
    v = M(1,1) * M(2,2) * M(3,3) + M(1, 2) * M(2, 3) * M(3, 1) + ...
        M(1,3) * M(2,1) * M(3,2) - M(1,3) * M(2,2) * M(3,1) - ...
        M(1,2) * M(2,1) * M(3,3) - M(1,1) * M(2,3) * M(3,2);
end

function v = indices(i, j, N)
    v = ((1:12) + 12 * (j - 1) + 12 * N * (i - 1));
end