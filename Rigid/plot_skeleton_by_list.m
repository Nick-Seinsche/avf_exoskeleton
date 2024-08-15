function plot_skeleton_by_list(U, h, orientation, t11, visual_height, plot_full_skeleton)
    L = size(U, 1) / 12;
    colors = ['r', 'b', 'g', 'm', 'c'];
    for ell=1:L
        plot3(U(2 + 12 * (ell-1)), ...
            U(1 + 12 * (ell-1)), ...
            U(3 + 12 * (ell-1)), "bo"); hold on
        
        Q = reshape(U((4:12) + 12 * (ell-1)), [3, 3]);
        x0 = U((1:3) + 12 * (ell-1));
        
        north = Q * [-h/2; 0; 0];
        west = Q * [0; -h/2; 0];
        east = Q * [0; h/2; 0];
        south = Q * [h/2; 0; 0];
        dirCross= [north, west, south, east];
        tij = vertcat((-1)^(orientation(ell)) * t11(1:2), t11(3));
        hNQ = Q * [0; tij(2); tij(3)];
        hWQ = Q * [tij(1); 0; tij(3)];
        hSQ = Q * [0; -tij(2); tij(3)];
        hEQ = Q * [-tij(1); 0; tij(3)];
        dirHinge = visual_height * [hNQ, hWQ, hSQ, hEQ];
        tijQ = visual_height * Q * tij;
        tij_minus_Q = visual_height * Q * [-tij(1); -tij(2); tij(3)];
        
        for l = 1:4
            if l == 3
                tijQ = tij_minus_Q;
            end

            if plot_full_skeleton
                plot3([x0(2) - tijQ(2), x0(2) + tijQ(2), ...
                x0(2) + dirCross(2,l) + dirHinge(2,l), ...
                x0(2) + dirCross(2,l) - dirHinge(2,l), ...
                x0(2) - tijQ(2)], ...
                [x0(1) - tijQ(1), x0(1) + tijQ(1), ...
                x0(1) + dirCross(1,l) + dirHinge(1,l), ...
                x0(1) + dirCross(1,l) - dirHinge(1,l), ...
                x0(1) - tijQ(1)], ...
                [x0(3) - tijQ(3), x0(3) + tijQ(3), ...
                x0(3) + dirCross(3,l) + dirHinge(3,l), ...
                x0(3) + dirCross(3,l) - dirHinge(3,l), ...
                x0(3) - tijQ(3)], colors(1 + mod(ell, 5)));
            end
            
            line([x0(2), x0(2)+dirCross(2,l)], ...
                [x0(1), x0(1)+dirCross(1,l)], ...
                [x0(3), x0(3)+dirCross(3,l)], ...
                    'Color', colors(1 + mod(ell, 5)));
        end
    end
    set(gca, 'YDir','reverse');
    grid on
    %xlim([-1.5, 1.5]);
    %ylim([-1.5, 1.5]);
    zlim([-1.5, 1.5]);
    xlabel('x-axis');
    ylabel('y-axis');
    zlabel('z-axis');
end

function draw_vector(start, direction)
    line([start(2), start(2) + direction(2)], [start(1), start(1) + direction(1)], ...
        [start(3), start(3) + direction(3)], "color", "red");
end