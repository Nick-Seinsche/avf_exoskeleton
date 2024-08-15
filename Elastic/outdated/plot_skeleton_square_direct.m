function plot_skeleton_square_direct(U, alpha, h, h_T)
    N = 2 / h;
    d = (12 + 4 * (2 * 2 + 4 * 1));
    alpha_rad = alpha / 360 * 2 * pi;
    wall_incline = tan(alpha_rad);
    t11 = [-tan(alpha_rad); tan(alpha_rad); 1];
    t11_v = rot_mat(atan(t11(2)),0,0) * t11;

    z1 = [h_T/2, -t11_v(1)/t11_v(3) * h_T/2];
    z4 = -z1;
    z3 = [h_T/2, h/2];
    z6 = [-h_T/2, h/2];
    z2 = 0.5 * (z3 + z1);
    z5 = 0.5 * (z4 + z6);

    z = [z1', z2', z3', z4', z5', z6'];

    QdirN = [0, -1; 0, 0; 1, 0];
    QdirS = [0, 1; 0, 0; 1, 0];
    QdirW = [0, 0; 0, -1; -1, 0];
    QdirE = [0, 0; 0, 1; -1, 0];
    for i=1:N
        for j=1:N
            %(-1)^(i+j)
            t11 = [-(-1)^(i+j)*tan(alpha_rad); (-1)^(i+j)*tan(alpha_rad); 1];
            t11_v = rot_mat(atan(t11(2)),0,0) * t11;
        
            z1 = [h_T/2, -t11_v(1)/t11_v(3) * h_T/2];
            z4 = -z1;
            z3 = [h_T/2, h/2];
            z6 = [-h_T/2, h/2];
            z2 = 0.5 * (z3 + z1);
            z5 = 0.5 * (z4 + z6);
        
            z = [z1', z2', z3', z4', z5', z6'];
            
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
            % z = [z1', z2', z3', z4', z5', z6'] in 3 times 6
            % there the zi are in 3d space

            plateN = repmat(xij, 1,6) + Qij * rot_mat(atan(t11(1)),0,0) * QdirN * z;
            plateS = repmat(xij, 1,6) + Qij * rot_mat(atan(t11(2)),0,0) * QdirS * z;
            plateW = repmat(xij, 1,6) + Qij * rot_mat(0,atan(t11(1)),0) * QdirW * z;
            plateE = repmat(xij, 1,6) + Qij * rot_mat(0,atan(t11(2)),0) * QdirE * z;

            % Since the in plane deformation is P1 discretized, the
            % inplane deformation at 0.5 * (p2 + p1) is 0.5 * sNi

            % the deflection boundary condition gives us a 0 constraint on
            % the derivative on z1 and z4 which will not be plotted  for
            % now

            % adding the deformation - id
            plateN = plateN + [zeros(3, 1), [0.5 * sN1; uN(1)], [sN1; uN(2)], zeros(3, 1), [0.5 * sN2; uN(3)], [sN2; uN(4)]];
            plateE = plateE + [zeros(3, 1), [0.5 * sE1; uE(1)], [sE1; uE(2)], zeros(3, 1), [0.5 * sE2; uE(3)], [sE2; uE(4)]];
            plateS = plateS + [zeros(3, 1), [0.5 * sS1; uS(1)], [sS1; uS(2)], zeros(3, 1), [0.5 * sS2; uS(3)], [sS2; uS(4)]];
            plateW = plateW + [zeros(3, 1), [0.5 * sW1; uW(1)], [sW1; uW(2)], zeros(3, 1), [0.5 * sW2; uW(3)], [sW2; uW(4)]];

            for k=1:6
                plot3(plateN(2,k), plateN(1,k), plateN(3,k), "bo"); hold on
                plot3(plateE(2,k), plateE(1,k), plateE(3,k), "ro");
                plot3(plateS(2,k), plateS(1,k), plateS(3,k), "go");
                plot3(plateW(2,k), plateW(1,k), plateW(3,k), "mo");

                % for now just fill w, later add the last two basis vectors
                % (which are zero), actually only one is needed
                % since one is positional and zero and the other
                % controls the derivative
                fc =fill3(plateN(2,[1,2,3,6,5,4]), plateN(1,[1,2,3,6,5,4]), plateN(3,[1,2,3,6,5,4]), "b");
                fc(1).FaceAlpha = 0.25;

                fc =fill3(plateE(2,[1,2,3,6,5,4]), plateE(1,[1,2,3,6,5,4]), plateE(3,[1,2,3,6,5,4]), "r");
                fc(1).FaceAlpha = 0.25;

                fc =fill3(plateS(2,[1,2,3,6,5,4]), plateS(1,[1,2,3,6,5,4]), plateS(3,[1,2,3,6,5,4]), "g");
                fc(1).FaceAlpha = 0.25;

                fc =fill3(plateW(2,[1,2,3,6,5,4]), plateW(1,[1,2,3,6,5,4]), plateW(3,[1,2,3,6,5,4]), "m");
                fc(1).FaceAlpha = 0.25;
            end

            draw_vector(xij, t11)
        end
    end

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