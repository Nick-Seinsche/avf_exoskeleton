function test_A3D()
    %% Visualize the sparse matrix
    k = 0;
    A = assemble_A3D(k);
    spy(A);

    %% Check if AX = 0 for init config
    k = 1;
    h = 2^(-k); % size parameter
    N = 2^(k+1); % size of the configuration
    L = N^2; % number of cross-like-structures
    h_T = h/4; % height parameter for a cross-like structure
    
    addpath('../casadi_windows') % addpath('../casadi_linux')
    import casadi.*
    at = affine_transform(alpha, k, h, h_T);

    U0 = init_configuration(k);
    A = assemble_A3(k);

    any(A * at(U0))

    %% Check if AX = 0 for 11East-12West non init config
    k = 1;
    h = 2^(-k); % size parameter
    N = 2^(k+1); % size of the configuration
    L = N^2; % number of cross-like-structures
    h_T = h/4; % height parameter for a cross-like structure
    d = (3 + 9 + 4 * (2 * 2 + 4 * 1)); % DOF's per cross

    U00 = init_configuration(k);
    
    % u11,E
    U00((21:22) + d * (1-1) + d * N * (1-1)) = [0.03,0.05];
    U00((23:24) + d * (1-1) + d * N * (1-1)) = [-0.25,-0.1];
    U00((25:28) + d * (1-1) + d * N * (1-1)) = [0,0.1,0,0.1];

    %u12,W
    U00((37:38) + d * (2-1) + d * N * (1-1)) = [0.03,-0.05];
    U00((39:40) + d * (2-1) + d * N * (1-1)) = [-0.25,0.1];
    U00((41:44) + d * (2-1) + d * N * (1-1)) = [0,-0.1,0,-0.1];

    addpath('../casadi_windows') % addpath('../casadi_linux')
    import casadi.*
    affine_transform = phi(alpha, k, h, h_T);

    A = assemble_A3(k);
    any(A * affine_transform(U00))
    
    %% Check if AX = 0 for 11South-21North non init config
    k = 1;
    h = 2^(-k); % size parameter
    N = 2^(k+1); % size of the configuration
    L = N^2; % number of cross-like-structures
    h_T = h/4; % height parameter for a cross-like structure
    d = (3 + 9 + 4 * (2 * 2 + 4 * 1)); % DOF's per cross

    U00 = init_configuration(k);
    
    % u11,S
    U00((29:30) + d * (1-1) + d * N * (1-1)) = [-0.01, 0.2];
    U00((31:32) + d * (1-1) + d * N * (1-1)) = [0.01, -0.2];
    U00((33:36) + d * (1-1) + d * N * (1-1)) = [0,0.15,0,0.01];

    %u21,N
    U00((13:14) + d * (1-1) + d * N * (2-1)) = [-0.01, -0.2];
    U00((15:16) + d * (1-1) + d * N * (2-1)) = [0.01, 0.2];
    U00((17:20) + d * (1-1) + d * N * (2-1)) = [0,-0.15,0,-0.01]';

    addpath('../casadi_windows') % addpath('../casadi_linux')
    import casadi.*
    affine_transform = phi(alpha, k, h, h_T);

    A = assemble_A3(k);
    any(A * affine_transform(U00))
    
end