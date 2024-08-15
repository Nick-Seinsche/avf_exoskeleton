function plot_exo_skeleton_animation(UT, alpha, h, h_T, M, obj)
    % Given a configuration U, a wall inclination alpha and problem size
    % parameters h and h_T, creates a 3D plot of the configuration
    % The plates are displayed using linear patches
    N = 2 / h;
    d = (3 + 9 + 4 * (2 * 2 + 4 * 1));
    alpha_rad = alpha / 360 * 2 * pi;
    wall_incline = tan(alpha_rad);

    plateNT = zeros(3, 6 * N^2 * M);
    plateET = zeros(3, 6 * N^2 * M);
    plateST = zeros(3, 6 * N^2 * M);
    plateWT = zeros(3, 6 * N^2 * M);
   
    for tt=1:M
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
                xij = UT((1:3) + d * (j-1) + d * N * (i-1), tt);
                Qij = reshape (UT((4:12) + d * (j-1) + d * N * (i-1), tt), 3, 3);
                
                sN1 = UT((13:14) + d * (j-1) + d * N * (i-1), tt);
                sN2 = UT((15:16) + d * (j-1) + d * N * (i-1), tt);
                uN = UT((17:20) + d * (j-1) + d * N * (i-1), tt);
                
                sE1 = UT((21:22) + d * (j-1) + d * N * (i-1), tt);
                sE2 = UT((23:24) + d * (j-1) + d * N * (i-1), tt);
                uE = UT((25:28) + d * (j-1) + d * N * (i-1), tt);
    
                sS1 = UT((29:30) + d * (j-1) + d * N * (i-1), tt);
                sS2 = UT((31:32) + d * (j-1) + d * N * (i-1), tt);
                uS = UT((33:36) + d * (j-1) + d * N * (i-1), tt);
    
                sW1 = UT((37:38) + d * (j-1) + d * N * (i-1), tt);
                sW2 = UT((39:40) + d * (j-1) + d * N * (i-1), tt);
                uW = UT((41:44) + d * (j-1) + d * N * (i-1), tt);
    
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
        
                % Adds the in-plane and out-of-plane deflection to the 3d
                % embeded reference blade. Then transforms this deviation
                % into the configuration space
                plateNT(:, (1:6) + 6 * j + 6 * N * (i - 1) + (tt - 1) * 6 * N^2) = repmat(xij, 1,6) + Qij * tdirN * QdirN * z3D + Qij * tdirN * QdirN * [zeros(3, 1), [0.5 * sN1; uN(1)], [sN1; uN(2)], zeros(3, 1), [0.5 * sN2; uN(3)], [sN2; uN(4)]];
                plateET(:, (1:6) + 6 * j + 6 * N * (i - 1) + (tt - 1) * 6 * N^2) = repmat(xij, 1,6) + Qij * tdirE * QdirE * z3D + Qij * tdirE * QdirE * [zeros(3, 1), [0.5 * sE1; uE(1)], [sE1; uE(2)], zeros(3, 1), [0.5 * sE2; uE(3)], [sE2; uE(4)]];
                plateST(:, (1:6) + 6 * j + 6 * N * (i - 1) + (tt - 1) * 6 * N^2) = repmat(xij, 1,6) + Qij * tdirS * QdirS * z3D + Qij * tdirS * QdirS * [zeros(3, 1), [0.5 * sS1; uS(1)], [sS1; uS(2)], zeros(3, 1), [0.5 * sS2; uS(3)], [sS2; uS(4)]];
                plateWT(:, (1:6) + 6 * j + 6 * N * (i - 1) + (tt - 1) * 6 * N^2) = repmat(xij, 1,6) + Qij * tdirW * QdirW * z3D + Qij * tdirW * QdirW * [zeros(3, 1), [0.5 * sW1; uW(1)], [sW1; uW(2)], zeros(3, 1), [0.5 * sW2; uW(3)], [sW2; uW(4)]];
            end
        end
    end
    for tt=1:M
        for i=1:N
            for j=1:N               
                % plot x-y-z axis
                fc =fill3(plateNT(2,[1,2,3,6,5,4] + 6 * j + 6 * N * (i - 1) + (tt - 1) * 6 * N^2), plateNT(1,[1,2,3,6,5,4] + 6 * j + 6 * N * (i - 1) + (tt - 1) * 6 * N^2), plateNT(3,[1,2,3,6,5,4] + 6 * j + 6 * N * (i - 1) + (tt - 1) * 6 * N^2), "b");  hold on
                fc(1).FaceAlpha = 0.25;

                set(gca, 'YDir','reverse'); grid on
                xlabel('x-axis');
                ylabel('y-axis');
                zlabel('z-axis');
                %view(-35,1);
                %view(-45,15);
                view(-45,85);
                    
                fc =fill3(plateET(2,[1,2,3,6,5,4] + 6 * j + 6 * N * (i - 1) + (tt - 1) * 6 * N^2), plateET(1,[1,2,3,6,5,4] + 6 * j + 6 * N * (i - 1) + (tt - 1) * 6 * N^2), plateET(3,[1,2,3,6,5,4] + 6 * j + 6 * N * (i - 1) + (tt - 1) * 6 * N^2), "r");
                fc(1).FaceAlpha = 0.25;

                fc =fill3(plateST(2,[1,2,3,6,5,4] + 6 * j + 6 * N * (i - 1) + (tt - 1) * 6 * N^2), plateST(1,[1,2,3,6,5,4] + 6 * j + 6 * N * (i - 1) + (tt - 1) * 6 * N^2), plateST(3,[1,2,3,6,5,4] + 6 * j + 6 * N * (i - 1) + (tt - 1) * 6 * N^2), "g");
                fc(1).FaceAlpha = 0.25;

                fc =fill3(plateWT(2,[1,2,3,6,5,4] + 6 * j + 6 * N * (i - 1) + (tt - 1) * 6 * N^2), plateWT(1,[1,2,3,6,5,4] + 6 * j + 6 * N * (i - 1) + (tt - 1) * 6 * N^2), plateWT(3,[1,2,3,6,5,4] + 6 * j + 6 * N * (i - 1) + (tt - 1) * 6 * N^2), "m");
                fc(1).FaceAlpha = 0.25;

                xlim([-1.4, 1.4]);
                ylim([-1.4, 1.4]);
                zlim([-1.4, 1.4]);
            end
        end
        frame = getframe(gcf);
        writeVideo(obj, frame);
        pause(0.005);
        if tt < M
            clf;
        end
    end
    hold off
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