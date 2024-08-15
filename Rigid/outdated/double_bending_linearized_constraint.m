function double_bending_linearized_constraint()
    k = 0;
    h = 2^(-k);
    N = 2^(k+1);
    L = N^2;
    
    t11 = [0.15; -0.15; 1];
    t11 = 0.5 * h * t11 / norm(t11, 2);
    
    [Alin, b0] = assemble_A_linearized(N, L, h, t11);
    A = assemble_A(N, L, h, t11);

    U0_lin = init_configuration_lin(k);

    dirichlet_point_lin = indices_lin(2,2,N);
    dirichlet_point_lin = dirichlet_point_lin(1);
    
    dp_lin = vertcat(indices_lin(1,1,N)', dirichlet_point_lin); 
    fp_lin = setdiff(1:(6 * L), dp_lin);
    %nfp_lin = size(fp_lin, 1);
    dpv_lin = vertcat(U0_lin(indices_lin(1,1,N)), U0_lin(dirichlet_point_lin)-0.2);
    
    Alin_lumped = Alin(:,fp_lin);
    b_lin = -Alin(:,dp_lin) * dpv_lin + b0;

    Usollinlumped = (Alin_lumped' * Alin_lumped) \ (Alin_lumped' * b_lin);
    
    Usollin = zeros(6 * L , 1);
    Usollin(dp_lin) = dpv_lin;
    Usollin(fp_lin) = Usollinlumped;

    Usol = convert_lin_to_full(Usollin, k);

    norm(A * Usol, 2)

    plot_skeleton(Usol, k, t11, 1)
end

function U0 = init_configuration_lin(k)
    N = 2^(k+1);
    h = 2^(-k);
    L = N^2;
    U0 = zeros(6 * L, 1);
    for i = 1:N
        for j = 1:N
            U0((1:2) + 6 * (j-1) + 6 * N * (i-1)) = [-1 + (i-1)*h + h/2, -1 + (j-1)*h + h/2];
            U0((4:6) + 6 * (j-1) + 6 * N * (i-1)) = [0, 0, 0];
        end
    end
end

function UU = convert_lin_to_full(U, k)
    N = 2^(k+1);
    h = 2^(-k);
    L = N^2;
    UU = zeros(12 * L, 1);
    for i = 1:N
        for j = 1:N
            UU((1:3) + 12 * (j-1) + 12 * N * (i-1)) = U((1:3) + 6 * (j-1) + 6 * N * (i-1));

            UU((4:12) + 12 * (j-1) + 12 * N * (i-1)) = reshape(eye(3) + [0, U(4 + 6 * (j-1) + 6 * N * (i-1)), ...
                U(5 + 6 * (j-1) + 6 * N * (i-1)); ...
                -U(4 + 6 * (j-1) + 6 * N * (i-1)), 0, U(6 + 6 * (j-1) + 6 * N * (i-1)); ...
                -U(5 + 6 * (j-1) + 6 * N * (i-1)), -U(6 + 6 * (j-1) + 6 * N * (i-1)), 0], 9, 1);
        end
    end
end

function v = indices_lin(i, j, N)
    v = ((1:6) + 6 * (j - 1) + 6 * N * (i - 1));
end