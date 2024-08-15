function plot_skeleton_old(U, k, t11)
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
                
                line([x0(1), x0(1)+dirCross(1,l)], ...
                    [x0(2), x0(2)+dirCross(2,l)], ...
                    [x0(3), x0(3)+dirCross(3,l)], ...
                        'Color', colors(1 + mod(j + (i-1) * N, 5)));
            end
        end
    end
    set(gca, 'YDir','reverse')
    grid on
    xlim([-1.5, 1.5]);
    ylim([-1.5, 1.5]);
    zlim([-1.5,1.5]);
end