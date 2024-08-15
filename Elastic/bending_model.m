function [energy, violation, affine_transform, U] = bending_model(k, dx, alpha, mu, do_plot)
    % k - determines the size of the configuration
    % dx - perturbation for the boundary condition
    % alpha - parameter for the wall inclination
    h = 2^(-k); % size parameter
    N = 2^(k+1); % size of the configuration
    L = N^2; % number of cross-like-structures
    h_T = h/4; % height parameter for a cross-like structure
    d = (3 + 9 + 4 * (2 * 2 + 4 * 1)); % DOF's per cross

    % IMPORTANT: There are different versions of Casadi for Windows
    % and Linux
    addpath('../casadi_windows') % addpath('../casadi_linux')
    import casadi.*
    
    U0 = init_configuration(k);

    % (x11, Q11) indices
    dirichlet_cross = indices(1,1,N);
    dirichlet_cross = dirichlet_cross(1:12);
    
    % z-component of xNN
    dirichlet_point = indices(N,N,N);
    dirichlet_point = dirichlet_point(3);

    % set indices above as fixed
    dp = vertcat(dirichlet_cross', dirichlet_point);
    fp = setdiff(1:(d * L), dp);
    nfp = size(fp, 1);
    % add a perturbation dx to the fixation
    dpv = vertcat(U0(dirichlet_cross), U0(dirichlet_point)+dx);
    
    % retrieve hinge condition sparse matrix
    A3 = assemble_A3(k);
    if do_plot
        figure(1);
        spy(A3);
    end
    
    % BUILD THE NLP:

    UU = SX.sym('x', d * L);

    % Rotation matrix determinant condition
    beta = SX.sym('beta', L);
    for l = 1:L
        beta(l) = det33(reshape(UU((4:12) + d * (l-1)), 3, 3));
    end

    % Rotation matrix in O(3) - condition
    theta = SX.sym('alpha', 6 * L);

    % condition: hinge conditions are satisfied
    UU_hinge = SX.zeros(24 * N ^ 2);

    % in the following the (relevant) DOFS from UU are transformed from 
    % their reference plate to the plate in the configuration and saved
    % in UU3. Finally A3 * UU3 = 0 guarantees that all neighboring
    % elements share a hinge:

    for i = 1:N
        for j=1:N
            theta(1+6*(j-1)+6*N*(i-1)) = (UU((4:6) + d * (j-1) + d * N * (i-1))' * UU((4:6) + d * (j-1) + d * N * (i-1)) - 1);
            theta(2+6*(j-1)+6*N*(i-1)) = (UU((4:6) + d * (j-1) + d * N * (i-1))' * UU((7:9) + d * (j-1) + d * N * (i-1)) - 0);
            theta(3+6*(j-1)+6*N*(i-1)) = (UU((4:6) + d * (j-1) + d * N * (i-1))' * UU((10:12) + d * (j-1) + d * N * (i-1)) - 0);
            theta(4+6*(j-1)+6*N*(i-1)) = (UU((7:9) + d * (j-1) + d * N * (i-1))' * UU((7:9) + d * (j-1) + d * N * (i-1)) - 1);
            theta(5+6*(j-1)+6*N*(i-1)) = (UU((7:9) + d * (j-1) + d * N * (i-1))' * UU((10:12) + d * (j-1) + d * N * (i-1)) - 0);
            theta(6+6*(j-1)+6*N*(i-1)) = (UU((10:12) + d * (j-1) + d * N * (i-1))' * UU((10:12) + d * (j-1) + d * N * (i-1)) - 1);

            % - Hinge conditions are satisfied ----
            alpha_rad = alpha / 360 * 2 * pi;
            wall_incline = tan(alpha_rad);
            t11 = [-(-1)^(i+j)*wall_incline; (-1)^(i+j)*wall_incline; 1];
            t11_v = rot_mat(atan(t11(2)),0,0) * t11;
        
            z1 = [h_T/2, -t11_v(1)/t11_v(3) * h_T/2];
            z4 = -z1;
            z3 = [h_T/2, h/2];
            z6 = [-h_T/2, h/2];
            z2 = 0.5 * (z3 + z1);
            z5 = 0.5 * (z4 + z6);
        
            z = [z1', z2', z3', z4', z5', z6'];
            z3D = [z; zeros(1, 6)];
            
            % unpack the entries for the ij-th cross
            xij = UU((1:3) + d * (j-1) + d * N * (i-1));
            Qij = reshape (UU((4:12) + d * (j-1) + d * N * (i-1)), 3, 3);
            
            sN1 = UU((13:14) + d * (j-1) + d * N * (i-1));
            sN2 = UU((15:16) + d * (j-1) + d * N * (i-1));
            uN = UU((17:20) + d * (j-1) + d * N * (i-1));
            
            sE1 = UU((21:22) + d * (j-1) + d * N * (i-1));
            sE2 = UU((23:24) + d * (j-1) + d * N * (i-1));
            uE = UU((25:28) + d * (j-1) + d * N * (i-1));

            sS1 = UU((29:30) + d * (j-1) + d * N * (i-1));
            sS2 = UU((31:32) + d * (j-1) + d * N * (i-1));
            uS = UU((33:36) + d * (j-1) + d * N * (i-1));

            sW1 = UU((37:38) + d * (j-1) + d * N * (i-1));
            sW2 = UU((39:40) + d * (j-1) + d * N * (i-1));
            uW = UU((41:44) + d * (j-1) + d * N * (i-1));

            % plate in reference configuration
            % it has the following form
            % z = [z1', z2', z3', z4', z5', z6'] in 2 times 6
            % there the zi are in 2d space, the z3D are embedded into 3d
            % space

            % QdirN transforms the 2d reference plate embeded into 3d space
            % to the northern blade without wall inclination
            % tdirN then adds the wall inclination
            [tdirN, QdirN] = transform3("N", t11);
            [tdirE, QdirE] = transform3("E", t11);
            [tdirS, QdirS] = transform3("S", t11);
            [tdirW, QdirW] = transform3("W", t11);

            % the four blades in reference configuration
            plateN = repmat(xij, 1,6) + Qij * tdirN * QdirN * z3D;
            plateS = repmat(xij, 1,6) + Qij * tdirS * QdirS * z3D;
            plateW = repmat(xij, 1,6) + Qij * tdirW * QdirW * z3D;
            plateE = repmat(xij, 1,6) + Qij * tdirE * QdirE * z3D;

            % Adds the in-plane and out-of-plane deflection to the 3d
            % embeded reference blade. Then transforms this deviation
            % into the configuration space
            deviationN = Qij * tdirN * QdirN * [zeros(3, 1), [0.5 * sN1; uN(1)], [sN1; uN(2)], zeros(3, 1), [0.5 * sN2; uN(3)], [sN2; uN(4)]];
            deviationE = Qij * tdirE * QdirE * [zeros(3, 1), [0.5 * sE1; uE(1)], [sE1; uE(2)], zeros(3, 1), [0.5 * sE2; uE(3)], [sE2; uE(4)]];
            deviationW = Qij * tdirW * QdirW * [zeros(3, 1), [0.5 * sW1; uW(1)], [sW1; uW(2)], zeros(3, 1), [0.5 * sW2; uW(3)], [sW2; uW(4)]];
            deviationS = Qij * tdirS * QdirS * [zeros(3, 1), [0.5 * sS1; uS(1)], [sS1; uS(2)], zeros(3, 1), [0.5 * sS2; uS(3)], [sS2; uS(4)]];

            % Since the in plane deformation is P1 discretized, the
            % inplane deformation at 0.5 * (p2 + p1) is 0.5 * sNi

            % the deflection boundary condition gives us a 0 constraint on
            % the derivative on z1 and z4 which will not be plotted  for
            % now

            % adding the deformation
            plateN = plateN + deviationN;
            plateE = plateE + deviationE;
            plateS = plateS + deviationS;
            plateW = plateW + deviationW;

            UU_hinge((1:24) + 24 * (j-1) + 24 * N * (i-1)) = reshape([plateN(:, [3, 6]), plateE(:, [3, 6]), plateS(:, [3, 6]), plateW(:, [3, 6])], [24, 1]);
        end
    end
    gamma = A3 * UU_hinge;
    
    affine_transform = Function('at', {UU}, {UU_hinge}, {'UU'}, {'UUh'});

    constraints = vertcat(theta, beta, gamma);

    lbx = -ones(d * L, 1) * inf;
    lbx(dp) = dpv; 
    ubx = ones(d * L, 1) * inf;
    ubx(dp) = dpv; 

    lbg = vertcat(zeros(6 * L, 1), 0 * ones(L, 1), zeros(12 * N * (N-1), 1));
    ubg = vertcat(zeros(6 * L, 1), inf * ones(L, 1), zeros(12 * N * (N-1), 1));

    %[Splus, Bplus, correctorplus] = stiffness_and_load(1, 0.25, alpha, 2);
    %[Sminus, Bminus, correctorminus] = stiffness_and_load(1, 0.25, alpha, 1);
    
    [LAMBDA_plus, toplength_plus, bottomlength_plus] = discrete_nabla_matrix(h, h_T, alpha, 2);
    [LAMBDA_minus, toplength_minus, bottomlength_minus] = discrete_nabla_matrix(h, h_T, alpha, 1);
    
    %writematrix(Splus, "cache/S_plus.txt")
    %writematrix(Bplus, "cache/B_plus.txt")
    %writematrix(correctorplus, "cache/corrector_plus.txt")

    %writematrix(Sminus, "cache/S_minus.txt")
    %writematrix(Bminus, "cache/B_minus.txt")
    %writematrix(correctorminus, "cache/corrector_minus.txt")


    % Reads the Stiffness, load and corrector matrix/vector from cache for
    % h = 1, h_T = 1/4. For a recalculation, uncomment the lines
    % above.
    Splus = readmatrix("cache/S_plus.txt");
    Bplus = readmatrix("cache/B_plus.txt");
    
    Bplus = h^2 * Bplus; %Apply the scaling
    
    correctorplus = readmatrix("cache/corrector_plus.txt");
    
    correctorplus = [1, 0; 0, h] * correctorplus; % Apply the scaling

    Sminus = readmatrix("cache/S_minus.txt");
    Bminus = readmatrix("cache/B_minus.txt");

    Bminus = h^2 * Bminus; % Apply the scaling

    correctorminus = readmatrix("cache/corrector_minus.txt");

    correctorminus = [1, 0; 0, h] * correctorminus; % Apply the scaling

    % indices to extract non-boundary condition points from W and G space
    % i.e Lambda_D is exdtracted from Lambda via the following indices:
    submatrixW = [2, 3, 6, 7];
    submatrixG = [2, 3, 5, 6];

    LAMBDA_plus = LAMBDA_plus(submatrixG, submatrixW);
    LAMBDA_minus = LAMBDA_minus(submatrixG, submatrixW);

    energyBend = SX.zeros(N,N);
    energyStretch = SX.zeros(N,N);
    e_gamma = 1;
    mu = 1;
    
    % The energy is determined in two steps:
    % First the matricies are chosen with the correct orientation
    % then loops over north west east south blade and sums up
    % the energy for one blade
    for i=1:N
        for j=1:N
            if mod(i+j,2) == 0
                S = Splus;
                B = Bplus;
                LAMBDA = LAMBDA_plus;
                corrector_inv = inv(correctorplus);
                toplength = toplength_plus;
                bottomlength = bottomlength_plus;
            elseif mod(i+j,2) == 1
                S = Sminus;
                B = Bminus;
                LAMBDA = LAMBDA_minus;
                corrector_inv = inv(correctorminus);
                toplength = toplength_minus;
                bottomlength = bottomlength_minus;
            end
            Uij = UU(indices(i,j,N));
            
            for p = [12, 20, 28, 36]
                % pure bending term ----------------------------------
                nabla_h_U = LAMBDA * Uij(p+5:p+8);
                energyBend(i,j) = energyBend(i,j) + e_gamma / 2 * nabla_h_U' * S * nabla_h_U;

                % stretch term -----------------------------------
                % coefficient vector for the stretch term
                ZD = SX.zeros(4, 1);

                % ZD(1) is Z_2 in the thesis
                
                % the following is the gradient of p at z2.
                % By the construction of the discrete, nabla,
                % the derivative in x direction at the top and
                % bottom edge of the plate is exact. Therefore
                % we may use nabla_h_U(1) for nabla_x w(z_{2,2})
                nabla_w = corrector_inv * [nabla_h_U(1), (Uij(p+5) - Uij(p+7))]';
                
                % the following is the gradient of s at z2
                % Uij(p+1)/2/(toplength/2) is the derivative of s1 at
                % z2 in x direction, where we use linearity to calculate
                % the constant slope.
                % Uij(p+1)/2 - Uij(p+3)/2 is just s(z_{2}) - w(z_{5})
                % in the thesis
                nabla_v = corrector_inv * [Uij(p+1)/2/(toplength/2), Uij(p+2)/2/(toplength/2); Uij(p+1)/2 - Uij(p+3)/2, Uij(p+2)/2 - Uij(p+4)/2]';
                
                % calculate the stretch term from the nablas
                stretch = nabla_w * nabla_w' + (nabla_v + nabla_v');

                % calculate Z_2 as in the thesis
                ZD(1) = stretch(1)^2 + stretch(2)^2 + stretch(3)^2 + stretch(4)^2;

                % ZD(2) is Z(3) in the thesis
                % the following is the gradient of p at z3
                % nabla_h_U(2) is nabla_x w(z_{3,2})
                % since the corrector matrix is diagonal, we just need
                % to divide the second entry by h_T
                nabla_w = [nabla_h_U(2), (Uij(p+6) - Uij(p+8))/h_T];

                % the following is the gradient of s at z_3
                % again the corrector matrix is diagonal, we just need
                % to divide by h_T
                nabla_v = [Uij(p+1)/toplength, Uij(p+2)/toplength; (Uij(p+1) - Uij(p+3))/h_T, (Uij(p+2) - Uij(p+4))/h_T]';
    
                stretch = nabla_w * nabla_w' + (nabla_v + nabla_v');
                ZD(2) = stretch(1)^2 + stretch(2)^2 + stretch(3)^2 + stretch(4)^2;

                % ZD(3) is Z(5) in the thesis
                % the following is the gradient of p at z5
                nabla_w = corrector_inv * [nabla_h_U(3), (Uij(p+5) - Uij(p+7))]';

                nabla_v = corrector_inv * [Uij(p+3)/2/(bottomlength/2), Uij(p+4)/2/(bottomlength/2); Uij(p+1)/2 - Uij(p+3)/2, Uij(p+2)/2 - Uij(p+4)/2]';
    
                stretch = nabla_w * nabla_w' + (nabla_v + nabla_v');
                ZD(3) = stretch(1)^2 + stretch(2)^2 + stretch(3)^2 + stretch(4)^2;
    
                % ZD(4) is Z(6)
                % the followng is the gradient of p at z6
                nabla_w = [nabla_h_U(4), (Uij(p+6) - Uij(p+8))/h_T];

                nabla_v = [Uij(p+3)/bottomlength, Uij(p+4)/bottomlength; (Uij(p+1) - Uij(p+3))/h_T, (Uij(p+2) - Uij(p+4))/h_T]';
    
                stretch = nabla_w * nabla_w' + (nabla_v + nabla_v');
                ZD(4) = stretch(1)^2 + stretch(2)^2 + stretch(3)^2 + stretch(4)^2;

                energyStretch(i,j) = energyStretch(i,j) + 0.5 * B' * ZD;

            end
        end
    end

    energy = sum(sum(energyBend)) + mu * sum(sum(energyStretch));

    nlp = struct('x', UU, 'f', energy, 'g', constraints);
    S = nlpsol('S', 'ipopt', nlp);

    sol = S('x0', U0, 'lbg', lbg, 'ubg', ubg, 'lbx', lbx, 'ubx', ubx);

    U = sol.x.full();

    measureEnergy = Function('m', {UU}, {energyBend + energyStretch}, {'U'}, {'e'});
    measureEnergyBend = Function('m', {UU}, {energyBend}, {'U'}, {'e'});
    measureEnergyStretch = Function('m', {UU}, {energyStretch}, {'U'}, {'e'});

    % measure raw constraint violation
    beta_prime = SX.sym('beta_p', L);
    for l = 1:L
        beta_prime(l) = det33(reshape(UU((4:12) + d * (l-1)), 3, 3)) - 1;
    end
    measureConstraintViolation = Function('c', {UU}, {norm(vertcat(theta, beta_prime, gamma), 2)}, 'U', 'v');
    %writematrix(U(indices(16,8,N)), "cache/solution_cross.txt")

    if do_plot
        figure(2);
        plot_exo_skeleton(U, alpha, h, h_T);
        xlim([-1.2, 1.2]);
        ylim([-1.2, 1.2]);
        zlim([-1.2, 1.2]);
    
        EnergyHM = full(measureEnergy(U));
        figure(3);
        heatmap(EnergyHM, 'ColorLimits', [min(min(EnergyHM)), max(max(EnergyHM))], 'Title', "Energy");
    
        EnergyBendHM = full(measureEnergyBend(U));
        figure(4);
        heatmap(EnergyBendHM, 'ColorLimits', [min(min(EnergyBendHM)), max(max(EnergyBendHM))], 'Title', 'Energy: Bend');
    
        EnergyStretchHM = full(measureEnergyStretch(U));
        figure(5);
        heatmap(EnergyStretchHM, 'ColorLimits', [min(min(EnergyStretchHM)), max(max(EnergyStretchHM))], 'Title', "Energy: Stretch");
    end

    energy = sum(sum(full(measureEnergy(U))));
    violation = full(measureConstraintViolation(U));
end

function v = indices(i,j, N)
    d = (3 + 9 + 4 * (2 * 2 + 4 * 1));
    v = ((1:d) + d * (j - 1) + d * N * (i - 1));
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
