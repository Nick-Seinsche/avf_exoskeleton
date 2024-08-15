function U0 = init_configuration(k)
    % given a problem size k, generates the planar / trivial configuration
    % for a bidirectional bending layer of that size
    h = 2^(-k);
    N = 2^(k+1);
    L = N^2;
    
    % The deflection is discretized using G_h (6 pts) and W_h (8pts)
    % The in-plane deflection is discretized using P1 (4 pts x 2)

    % Applying the boundary conditions, we have W_h (4 pts), G_h (4pts),
    % and P1 (2 pts x 2)

    % we have N^2 crosses. For each cross:
    % INDEX 1-3: POSITION xij
    % INDEX 4-12: ROTATION QIJ
    % IN PLANE DEFORMATION NORTH at p2 (2d vector)
    % IN PLANE DEFORMATION NORTH at p4 (2d vector)
    % DEFLECTION NORTH at z2, z3, z5, z6 (4d vector)
    % SAME FOR EAST, SOUTH, WEST

    % Once assembled, the vector U could be dissasembled as follows:
    % FOR i,j = 1...N
        %xij = U((1:3) + d * (j-1) + d * N * (i-1));
        %Qij = reshape (U((4:12) + d * (j-1) + d * N * (i-1)), 3, 3);
        
        %sN1 = U((13:14) + d * (j-1) + d * N * (i-1));
        %sN2 = U((15:16) + d * (j-1) + d * N * (i-1));
        %uN = U((17:20) + d * (j-1) + d * N * (i-1));
        
        %sE1 = U((21:22) + d * (j-1) + d * N * (i-1));
        %sE2 = U((23:24) + d * (j-1) + d * N * (i-1));
        %uE = U((25:28) + d * (j-1) + d * N * (i-1));
    
        %sS1 = U((29:30) + d * (j-1) + d * N * (i-1));
        %sS2 = U((31:32) + d * (j-1) + d * N * (i-1));
        %uS = U((33:36) + d * (j-1) + d * N * (i-1));
    
        %sW1 = U((37:38) + d * (j-1) + d * N * (i-1));
        %sW2 = U((39:40) + d * (j-1) + d * N * (i-1));
        %uW = U((41:44) + d * (j-1) + d * N * (i-1));
    % END

    d = (3 + 9 + 4 * (2 * 2 + 4 * 1));
    U0 = zeros(d * L, 1);

    for i = 1:N
        for j = 1:N
            U0((1:2) + d * (j-1) + d * N * (i-1)) = [-1 + (i-1) * h + h/2, -1 + (j-1) * h + h/2];
            % U0(3) = 0;
            U0((4:12) + d * (j-1) + d * N * (i-1)) = reshape(eye(3), 9, 1);
            U0((13:44) + d * (j-1) + d * N * (i-1)) = zeros(32, 1);
        end
    end
end