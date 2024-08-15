function A = assemble_A(k, t11)
    % Assembles the hinge condition matrix.
    % k - size of the configuration
    % t11 - crease vector
    h = 2^(-k); N = 2^(k+1); L = N^2;
    eq_ctr = 1; 
    
    I = zeros(20 * L); J = zeros(20 * L); V = zeros(20 * L);
    idx_ctr = 1;

    for i = 1:N
        for j = 1:N
            if i < N
                for k = 0:2
                    I(idx_ctr) = eq_ctr+k;
                    J(idx_ctr) = 1+12*(j-1)+(i-1) * 12 * N+k;
                    V(idx_ctr) = 1;
                    idx_ctr = idx_ctr + 1;
                    %A(eq_ctr+k, 1+12*(j-1)+(i-1) * 12 * N+k) = 1;

                    I(idx_ctr) = eq_ctr+k;
                    J(idx_ctr) = 1+12*(j-1)+(i-1) * 12 * N+3+k;
                    V(idx_ctr) = 0.5 * h;
                    idx_ctr = idx_ctr + 1;
                    %A(eq_ctr+k, 1+12*(j-1)+(i-1) * 12 * N+3+k) = 0.5 * h;
                    
                    I(idx_ctr) = eq_ctr+k;
                    J(idx_ctr) = 1+12*(j-1)+i * 12 * N+k;
                    V(idx_ctr) = -1;
                    idx_ctr = idx_ctr + 1;
                    %A(eq_ctr+k, 1+12*(j-1)+i * 12 * N+k) = -1;

                    I(idx_ctr) = eq_ctr+k;
                    J(idx_ctr) = 1+12*(j-1)+i * 12 * N+3+k;
                    V(idx_ctr) = 0.5 * h;
                    idx_ctr = idx_ctr + 1;
                    %A(eq_ctr+k, 1+12*(j-1)+i * 12 * N+3+k) = 0.5 * h;
                end
                eq_ctr = eq_ctr + 3;
            end
            if j < N
                for k = 0:2
                    I(idx_ctr) = eq_ctr+k;
                    J(idx_ctr) = 1+12*(j-1)+(i-1) * 12 * N+k;
                    V(idx_ctr) = 1;
                    idx_ctr = idx_ctr + 1;
                    %A(eq_ctr+k, 1+12*(j-1)+(i-1) * 12 * N+k) = 1;

                    I(idx_ctr) = eq_ctr+k;
                    J(idx_ctr) = 1+12*(j-1)+(i-1) * 12 * N+6+k;
                    V(idx_ctr) = 0.5 * h;
                    idx_ctr = idx_ctr + 1;
                    %A(eq_ctr+k, 1+12*(j-1)+(i-1) * 12 * N+6+k) = 0.5 * h;
                    
                    I(idx_ctr) = eq_ctr+k;
                    J(idx_ctr) = 1+12*(j-1)+(i-1) * 12 * N+k+12;
                    V(idx_ctr) = -1;
                    idx_ctr = idx_ctr + 1;
                    %A(eq_ctr+k, 1+12*(j-1)+(i-1) * 12 * N+k+12) = -1;

                    I(idx_ctr) = eq_ctr+k;
                    J(idx_ctr) = 1+12*(j-1)+(i-1) * 12 * N+k+12+6;
                    V(idx_ctr) = 0.5 * h;
                    idx_ctr = idx_ctr + 1;
                    %A(eq_ctr+k, 1+12*(j-1)+(i-1) * 12 * N+k+12+6) = 0.5 * h;
                end
                eq_ctr = eq_ctr + 3;
            end
            if i < N
                tij = vertcat((-1)^(i+j) * t11(1:2), t11(3));
                v = [0, -tij(2), tij(3)];

                for k = 0:2
                    I(idx_ctr) = eq_ctr+k;
                    J(idx_ctr) = 1+12*(j-1)+(i-1)* 12 * N+k+3;
                    V(idx_ctr) = v * [1;0;0];
                    idx_ctr = idx_ctr + 1;
                    %A(eq_ctr+k, 1+12*(j-1)+(i-1)* 12 * N+k+3) = v * [1;0;0];
                    
                    I(idx_ctr) = eq_ctr+k;
                    J(idx_ctr) = 1+12*(j-1)+(i-1)* 12 * N+k+6;
                    V(idx_ctr) = v * [0;1;0];
                    idx_ctr = idx_ctr + 1;
                    %A(eq_ctr+k, 1+12*(j-1)+(i-1)* 12 * N+k+6) =  v * [0;1;0];
                    
                    I(idx_ctr) = eq_ctr+k;
                    J(idx_ctr) = 1+12*(j-1)+(i-1) * 12 * N+k+9;
                    V(idx_ctr) = v * [0;0;1];
                    idx_ctr = idx_ctr + 1;
                    %A(eq_ctr+k, 1+12*(j-1)+(i-1) * 12 * N+k+9) = v * [0;0;1];
                    
                    I(idx_ctr) = eq_ctr+k;
                    J(idx_ctr) = 1+12*(j-1)+i*12*N+k+3;
                    V(idx_ctr) = -v * [1; 0; 0];
                    idx_ctr = idx_ctr + 1;
                    %A(eq_ctr+k, 1+12*(j-1)+i*12*N+k+3) = -v * [1; 0; 0];
                    
                    I(idx_ctr) = eq_ctr+k;
                    J(idx_ctr) = 1+12*(j-1)+i*12*N+k+6;
                    V(idx_ctr) = -v * [0; 1; 0];
                    idx_ctr = idx_ctr + 1;
                    %A(eq_ctr+k, 1+12*(j-1)+i*12*N+k+6) = -v * [0; 1; 0];
                    
                    I(idx_ctr) = eq_ctr+k;
                    J(idx_ctr) = 1+12*(j-1)+i*12*N+k+9;
                    V(idx_ctr) = -v * [0; 0; 1];
                    idx_ctr = idx_ctr + 1;
                    %A(eq_ctr+k, 1+12*(j-1)+i*12*N+k+9) = -v * [0; 0; 1];
                end
                eq_ctr = eq_ctr + 3;
            end
            if j < N
                tij = vertcat((-1)^(i+j) * t11(1:2), t11(3));
                v = [-tij(1), 0, tij(3)];

                for k = 0:2
                    I(idx_ctr) = eq_ctr+k;
                    J(idx_ctr) = 1+12*(j-1)+(i-1)*12*N+k+3;
                    V(idx_ctr) = v * [1;0;0];
                    idx_ctr = idx_ctr + 1;
                    %A(eq_ctr+k, 1+12*(j-1)+(i-1)*12*N+k+3) = v * [1;0;0];
                    
                    I(idx_ctr) = eq_ctr+k;
                    J(idx_ctr) = 1+12*(j-1)+(i-1)*12*N+k+6;
                    V(idx_ctr) = v * [0;1;0];
                    idx_ctr = idx_ctr + 1;
                    %A(eq_ctr+k, 1+12*(j-1)+(i-1)*12*N+k+6) = v * [0;1;0];
                    
                    I(idx_ctr) = eq_ctr+k;
                    J(idx_ctr) = 1+12*(j-1)+(i-1)*12*N+k+9;
                    V(idx_ctr) = v * [0;0;1];
                    idx_ctr = idx_ctr + 1;
                    %A(eq_ctr+k, 1+12*(j-1)+(i-1)*12*N+k+9) = v * [0;0;1];
                    
                    I(idx_ctr) = eq_ctr+k;
                    J(idx_ctr) = 1+12*(j-1)+(i-1)*12*N+k+3+12;
                    V(idx_ctr) = -v * [1;0;0];
                    idx_ctr = idx_ctr + 1;
                    %A(eq_ctr+k, 1+12*(j-1)+(i-1)*12*N+k+3+12) = -v * [1;0;0];
                    
                    I(idx_ctr) = eq_ctr+k;
                    J(idx_ctr) = 1+12*(j-1)+(i-1)*12*N+k+6+12;
                    V(idx_ctr) = -v * [0;1;0];
                    idx_ctr = idx_ctr + 1;
                    %A(eq_ctr+k, 1+12*(j-1)+(i-1)*12*N+k+6+12) = -v * [0;1;0];
                    
                    I(idx_ctr) = eq_ctr+k;
                    J(idx_ctr) = 1+12*(j-1)+(i-1)*12*N+k+9+12;
                    V(idx_ctr) = -v * [0;0;1];
                    idx_ctr = idx_ctr + 1;
                    %A(eq_ctr+k, 1+12*(j-1)+(i-1)*12*N+k+9+12) = -v * [0;0;1];
                end
                eq_ctr = eq_ctr + 3;
            end
        end
    end
    I = I(1:idx_ctr-1);
    J = J(1:idx_ctr-1);
    V = V(1:idx_ctr-1);
    A = sparse(I, J, V, 12 * N * (N-1), 12 * L);
end