function test_skeleton_square()
    %k = -1;
    %h = 2^(-k);
    %h_T = 0.5 * h;
    %alpha = 20.9;
    %U0 = init_configuration(k);
   
    %xij = zeros(3, 1);
    %Qij = reshape (eye(3,3), 9, 1);
    
    %sN1 = zeros(2,1);
    %sN2 = zeros(2,1);
    %uN = [0;0.1;0;0.1];
    
    %sE1 = zeros(2,1);
    %sE2 = zeros(2,1);
    %uE = [0;0.2;0;0.2];
    
    %sS1 = zeros(2,1);
    %sS2 = zeros(2,1);
    %uS = [0;0.1;0;0.1];

    %sW1 = zeros(2,1);
    %sW2 = zeros(2,1);
    %uW = [0;0.1;0;0.1];

    %U = [xij; Qij; sN1; sN2; uN; sE1; sE2; uE; sS1; sS2; uS; sW1; sW2; uW];

    %plot_skeleton_square_old(U0, alpha, h, h_T);
    %plot_skeleton_square_by_reference(U0, alpha, h, h_T);
    
    %% Test Case: Reference Configuration
    
    k = 0;
    U00 = init_configuration(k);

    h = 2^(-k);
    h_T = 0.5 * h;
    alpha = 20.9;

    plot_exo_skeleton(U00, alpha, h, h_T);
    
    %% Test Case: Non Reference Configuration 11East-12West

    k = 0;
    h = 2^(-k);
    N = 2^(k+1);
    L = N^2;
    U00 = init_configuration(k);
    d = (3 + 9 + 4 * (2 * 2 + 4 * 1));
    
    % u11,E
    U00((21:22) + d * (1-1) + d * N * (1-1)) = [0.03,0.05];
    U00((23:24) + d * (1-1) + d * N * (1-1)) = -[-0.25,-0.1];
    %U00((25:28) + d * (1-1) + d * N * (1-1)) = [0,0.1,0,0.1];

    %u12,W
    U00((37:38) + d * (2-1) + d * N * (1-1)) = [0.03,-0.05];
    U00((39:40) + d * (2-1) + d * N * (1-1)) = -[-0.25,0.1];
    %U00((41:44) + d * (2-1) + d * N * (1-1)) = [0,-0.1,0,-0.1];

    h = 2^(-k);
    h_T = 0.5 * h;
    alpha = 20.9;
    plot_exo_skeleton(U00, alpha, h, h_T);

    %% Test Case: Non Reference Configuration 21East-22West

    k = 0;
    h = 2^(-k);
    N = 2^(k+1);
    L = N^2;
    U00 = init_configuration(k);
    d = (3 + 9 + 4 * (2 * 2 + 4 * 1));
    
    % u11,E
    U00((21:22) + d * (1-1) + d * N * (2-1)) = [0.03,0.05];
    U00((23:24) + d * (1-1) + d * N * (2-1)) = [-0.25,-0.1];
    U00((25:28) + d * (1-1) + d * N * (2-1)) = [0,0.1,0,0.1];

    %u12,W
    U00((37:38) + d * (2-1) + d * N * (2-1)) = [0.03,-0.05];
    U00((39:40) + d * (2-1) + d * N * (2-1)) = [-0.25,0.1];
    U00((41:44) + d * (2-1) + d * N * (2-1)) = [0,-0.1,0,-0.1];

    h = 2^(-k);
    h_T = 0.5 * h;
    alpha = 20.9;
    plot_exo_skeleton(U00, alpha, h, h_T);

    %% Test Case: Non Reference Configuration 11South-21North

    k = 0;
    h = 2^(-k);
    N = 2^(k+1);
    L = N^2;
    U00 = init_configuration(k);
    d = (3 + 9 + 4 * (2 * 2 + 4 * 1));
    
    % u11,S
    U00((29:30) + d * (1-1) + d * N * (1-1)) = [-0.01, 0.2];
    U00((31:32) + d * (1-1) + d * N * (1-1)) = [0.01, -0.2];
    U00((33:36) + d * (1-1) + d * N * (1-1)) = [0,0.15,0,0.01];

    %u21,N
    U00((13:14) + d * (1-1) + d * N * (2-1)) = [-0.01, -0.2];
    U00((15:16) + d * (1-1) + d * N * (2-1)) = [0.01, 0.2];
    U00((17:20) + d * (1-1) + d * N * (2-1)) = [0,-0.15,0,-0.01]';

    h = 2^(-k);
    h_T = 0.5 * h;
    alpha = 20.9;
    plot_exo_skeleton(U00, alpha, h, h_T);

    %% Test Case: Non Reference Configuration 12South-22North

    k = 0;
    h = 2^(-k);
    N = 2^(k+1);
    L = N^2;
    U00 = init_configuration(k);
    d = (3 + 9 + 4 * (2 * 2 + 4 * 1));
    
    % u12,S
    U00((29:30) + d * (2-1) + d * N * (1-1)) = [-0.01, 0.2];
    U00((31:32) + d * (2-1) + d * N * (1-1)) = [0.01, -0.2];
    U00((33:36) + d * (2-1) + d * N * (1-1)) = [0,0.15,0,0.01];

    %u22,N
    U00((13:14) + d * (2-1) + d * N * (2-1)) = [-0.01, -0.2];
    U00((15:16) + d * (2-1) + d * N * (2-1)) = [0.01, 0.2];
    U00((17:20) + d * (2-1) + d * N * (2-1)) = [0,-0.15,0,-0.01]';

    h = 2^(-k);
    h_T = 0.5 * h;
    alpha = 20.9;
    plot_exo_skeleton(U00, alpha, h, h_T);
end