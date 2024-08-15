function animation()
    k = 1;
    N = 2^(k+1);
    L = N^2; 
    d = (3 + 9 + 4 * (2 * 2 + 4 * 1));
    alpha = 20.9;
    h = 2^(-k); % size parameter
    h_T = h/4; % height parameter for a cross-like structure

    M = 80;

    T = linspace(-0.9, 0.9, M);

    %UT = zeros(d * L, M);

    %i = 1;
    
    %for dx = T
    %    [~, ~, ~, U] = bending_model(1, dx, 20.9, 1, 0);
    %    UT(:,i) = U;
    %    i = i+1;
    %end

    %writematrix(UT, 'cache/aniU.txt');

    UT = readmatrix('cache/aniU.txt');

    newVid = VideoWriter('video3', 'MPEG-4'); % New
    newVid.FrameRate = 30;
    newVid.Quality = 100;
    open(newVid);

    plot_exo_skeleton_animation(UT, alpha, h, h_T, M, newVid);
    pause(0.5);
    plot_exo_skeleton_animation(UT(:,M:-1:1), alpha, h, h_T, M, newVid);
    pause(0.5);
    plot_exo_skeleton_animation(UT, alpha, h, h_T, M, newVid);

    close(newVid);
end