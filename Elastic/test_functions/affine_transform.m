function affine_transform = phi(alpha, k, h, h_T)
    % in the following the (relevant) DOFS from UU are transformed from 
    % their reference plate to the plate in the configuration and saved
    % in UU3. Finally A3 * UU3 = 0 guarantees that all neighboring
    % elements share a hinge:
    addpath('../casadi_windows') % addpath('../casadi_linux')
    import casadi.*

    N = 2^(k+1); % size of the configuration
    L = N^2; % number of cross-like-structures
    d = (3 + 9 + 4 * (2 * 2 + 4 * 1)); % DOF's per cross

    UU = SX.sym('x', d * L);
    UU_hinge = SX.zeros(24 * N ^ 2);

    % in the following the (relevant) DOFS from UU are transformed from 
    % their reference plate to the plate in the configuration and saved
    % in UU3. Finally A3 * UU3 = 0 guarantees that all neighboring
    % elements share a hinge:

    for i = 1:N
        for j=1:N
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
    affine_transform = Function('at', {UU}, {UU_hinge}, {'UU'}, {'UUh'});
end

function Q = rot_mat(a,b,c)
    Q = [1, 0, 0; 0, cos(a), -sin(a); 0, sin(a), cos(a)];
    Q = Q * [cos(b), 0, sin(b); 0, 1, 0; -sin(b), 0, cos(b)];
    Q = Q * [cos(c), -sin(c), 0; sin(c), cos(c), 0; 0,0,1];
end