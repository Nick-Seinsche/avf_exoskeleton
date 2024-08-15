function plot_exo_skeleton(U, alpha, h, h_T)
    % Given a configuration U, a wall inclination alpha and problem size
    % parameters h and h_T, creates a 3D plot of the configuration
    % The plates are displayed using linear patches
    N = 2 / h;
    d = (3 + 9 + 4 * (2 * 2 + 4 * 1));
    alpha_rad = alpha / 360 * 2 * pi;
    wall_incline = tan(alpha_rad);
   
    for i=1:N
        for j=1:N
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
            xij = U((1:3) + d * (j-1) + d * N * (i-1));
            Qij = reshape (U((4:12) + d * (j-1) + d * N * (i-1)), 3, 3);
            
            sN1 = U((13:14) + d * (j-1) + d * N * (i-1));
            sN2 = U((15:16) + d * (j-1) + d * N * (i-1));
            uN = U((17:20) + d * (j-1) + d * N * (i-1));
            
            sE1 = U((21:22) + d * (j-1) + d * N * (i-1));
            sE2 = U((23:24) + d * (j-1) + d * N * (i-1));
            uE = U((25:28) + d * (j-1) + d * N * (i-1));

            sS1 = U((29:30) + d * (j-1) + d * N * (i-1));
            sS2 = U((31:32) + d * (j-1) + d * N * (i-1));
            uS = U((33:36) + d * (j-1) + d * N * (i-1));

            sW1 = U((37:38) + d * (j-1) + d * N * (i-1));
            sW2 = U((39:40) + d * (j-1) + d * N * (i-1));
            uW = U((41:44) + d * (j-1) + d * N * (i-1));

            % plate in reference configuration
            % it has the following form
            % z = [z1', z2', z3', z4', z5', z6'] is a 2x6 matrix
            % z3D has zeros filled for the z-axis ie. a 3x6 matrix

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
                
            % inplane deformation at 0.5 * (p2 + p1) is 0.5 * sNi

            % the deflection boundary condition gives us a 0 constraint on
            % the derivative at z1 and z4 which is neglected due to using
            % patches. Plotting the actual deformation of the plate
            % would be to computationally intensive.

            % adding the deformation
            plateN = plateN + deviationN;
            plateE = plateE + deviationE;
            plateS = plateS + deviationS;
            plateW = plateW + deviationW;

            for k=1:6
                plot3(plateN(2,k), plateN(1,k), plateN(3,k), "bo"); hold on
                plot3(plateE(2,k), plateE(1,k), plateE(3,k), "ro");
                plot3(plateS(2,k), plateS(1,k), plateS(3,k), "go");
                plot3(plateW(2,k), plateW(1,k), plateW(3,k), "mo");
            end
                % plot x-y-z axis
                fc =fill3(plateN(2,[1,2,3,6,5,4]), plateN(1,[1,2,3,6,5,4]), plateN(3,[1,2,3,6,5,4]), "b");
                fc(1).FaceAlpha = 0.25;

                fc =fill3(plateE(2,[1,2,3,6,5,4]), plateE(1,[1,2,3,6,5,4]), plateE(3,[1,2,3,6,5,4]), "r");
                fc(1).FaceAlpha = 0.25;

                fc =fill3(plateS(2,[1,2,3,6,5,4]), plateS(1,[1,2,3,6,5,4]), plateS(3,[1,2,3,6,5,4]), "g");
                fc(1).FaceAlpha = 0.25;

                fc =fill3(plateW(2,[1,2,3,6,5,4]), plateW(1,[1,2,3,6,5,4]), plateW(3,[1,2,3,6,5,4]), "m");
                fc(1).FaceAlpha = 0.25;

            % draw vector perpendicular to the bidirectional bending plane
            draw_vector(xij, Qij * [0;0;1] / 5);
        end
    end
    set(gca, 'YDir','reverse'); grid on
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
        [start(3), start(3) + direction(3)], "color", "red");
end