function cross_like_structure()
    % Plots a cross-like-structure for visualization
    k = -1;
    h = 2^(-k);
    h_T = 0.5 * h;
    N = 2^(k+1);
    L = N^2;
    alpha_rad = 21.6 / 360 * 2 * pi;
    w = h_T * tan(alpha_rad) / 2;
    t11 = [-w; w; h_T/2];
    %t11 = h_T * t11;
    
    plot_skeleton_square(init_configuration(k), k, t11, 1, 1);

    draw_vector([-1,-1,-h_T/2], [2,0,0]);
    draw_vector([1,-1,-h_T/2], [0,2,0]);
    draw_vector([1,1,-h_T/2], [-2,0,0]);
    draw_vector([-1,1,-h_T/2], [0,-2,0]);

    draw_vector([-1,-1,h_T/2], [2,0,0]);
    draw_vector([1,-1,h_T/2], [0,2,0]);
    draw_vector([1,1,h_T/2], [-2,0,0]);
    draw_vector([-1,1,h_T/2], [0,-2,0]);

    draw_vector([-1,-1,-h_T/2], [0,0,h_T]);
    draw_vector([-1,1,-h_T/2], [0,0,h_T]);
    draw_vector([1,-1,-h_T/2], [0,0,h_T]);
    draw_vector([1,1,-h_T/2], [0,0,h_T]);

    draw_vector2([0,0,0], t11);
end

function draw_vector(start, direction)
    line([start(2), start(2) + direction(2)], [start(1), start(1) + direction(1)], ...
        [start(3), start(3) + direction(3)], "color", "black");
end

function draw_vector2(start, direction)
    line([start(2), start(2) + direction(2)], [start(1), start(1) + direction(1)], ...
        [start(3), start(3) + direction(3)], "color", "green", 'LineWidth',3);
end
