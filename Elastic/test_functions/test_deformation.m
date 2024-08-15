function test_deformation()
    h = 2;
    h_T = h/4; % height parameter for a cross-like structure
    alpha = 20.9;
    U = readmatrix("cache/solution_cross.txt");

    U(1:3) = 0;
    U(4:12) = reshape(eye(3), 9, []);

    U(13:end) = 50000 * U(13:end) ;

    plot_exo_skeleton(U, alpha, h, h_T);
end