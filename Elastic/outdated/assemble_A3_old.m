function A = assemble_A3_old(k)
    % Assembles the sparse matrix for the hinge conditions
    % k is the size of the double bending configuration
    % A is the sparse matrix
    N = 2^(k+1);
    L = N^2;
    d = (12 + 4 * (2 * 2 + 4 * 1));

    A = sparse(12 * N * (N - 1), 24 * L);
    equation_counter = 1;

    I = zeros(24 * N ^ 2); J = zeros(24 * N ^ 2); V = zeros(24 * N ^ 2);
    index_counter = 1;
    
    % Indices of positional coordinates inside U

    % N3 | 1 2 3    
    % N6 | 4 5 6    
    %
    % E3 | 7 8 9    
    % E6 | 10 11 12
    %
    % S3 | 13 14 15    
    % S6 | 16 17 18
    %
    % W3 | 19 20 21    
    % W6 | 22 23 24    

    % i_j_S3 means:
    % cross (i,j), southern blade, position at z3
    % i_j+1_N6 means:
    % cross (i,j+1), northern blade, position at z6

    for i = 1:N
        for j = 1:N
            if i < N
                %i_j_S3 = (i+1)_j_N3
                I(index_counter) = equation_counter;
                J(index_counter) = 13 + 24 * (j-1) + 24 * N * (i-1);
                V(index_counter) = 1;
                index_counter = index_counter + 1;

                I(index_counter) = equation_counter;
                J(index_counter) = 1 + 24 * (j-1) + 24 * N * (i-0);
                V(index_counter) = -1;
                index_counter = index_counter + 1;

                A(equation_counter, 13 + 24 * (j-1) + 24 * N * (i-1)) = 1;
                A(equation_counter, 1 + 24 * (j-1) + 24 * N * (i-0)) = -1;
                equation_counter = equation_counter + 1;

                % ----

                I(index_counter) = equation_counter;
                J(index_counter) = 14 + 24 * (j-1) + 24 * N * (i-1);
                V(index_counter) = 1;
                index_counter = index_counter + 1;

                I(index_counter) = equation_counter;
                J(index_counter) = 2 + 24 * (j-1) + 24 * N * (i-0);
                V(index_counter) = -1;
                index_counter = index_counter + 1;

                A(equation_counter, 14 + 24 * (j-1) + 24 * N * (i-1)) = 1;
                A(equation_counter, 2 + 24 * (j-1) + 24 * N * (i-0)) = -1;
                equation_counter = equation_counter + 1;

                % ---

                I(index_counter) = equation_counter;
                J(index_counter) = 15 + 24 * (j-1) + 24 * N * (i-1);
                V(index_counter) = 1;
                index_counter = index_counter + 1;

                I(index_counter) = equation_counter;
                J(index_counter) = 3 + 24 * (j-1) + 24 * N * (i-0);
                V(index_counter) = -1;
                index_counter = index_counter + 1;

                A(equation_counter, 15 + 24 * (j-1) + 24 * N * (i-1)) = 1;
                A(equation_counter, 3 + 24 * (j-1) + 24 * N * (i-0)) = -1;
                equation_counter = equation_counter + 1;

                % ---

                %i_j_S6 = (i+1)_j_N6

                I(index_counter) = equation_counter;
                J(index_counter) = 16 + 24 * (j-1) + 24 * N * (i-1);
                V(index_counter) = 1;
                index_counter = index_counter + 1;

                I(index_counter) = equation_counter;
                J(index_counter) = 4 + 24 * (j-1) + 24 * N * (i-0);
                V(index_counter) = -1;
                index_counter = index_counter + 1;

                A(equation_counter, 16 + 24 * (j-1) + 24 * N * (i-1)) = 1;
                A(equation_counter, 4 + 24 * (j-1) + 24 * N * (i-0)) = -1;
                equation_counter = equation_counter + 1;

                % ---

                I(index_counter) = equation_counter;
                J(index_counter) = 17 + 24 * (j-1) + 24 * N * (i-1);
                V(index_counter) = 1;
                index_counter = index_counter + 1;

                I(index_counter) = equation_counter;
                J(index_counter) = 5 + 24 * (j-1) + 24 * N * (i-0);
                V(index_counter) = -1;
                index_counter = index_counter + 1;

                A(equation_counter, 17 + 24 * (j-1) + 24 * N * (i-1)) = 1;
                A(equation_counter, 5 + 24 * (j-1) + 24 * N * (i-0)) = -1;
                equation_counter = equation_counter + 1;

                % ---

                I(index_counter) = equation_counter;
                J(index_counter) = 18 + 24 * (j-1) + 24 * N * (i-1);
                V(index_counter) = 1;
                index_counter = index_counter + 1;

                I(index_counter) = equation_counter;
                J(index_counter) = 6 + 24 * (j-1) + 24 * N * (i-0);
                V(index_counter) = -1;
                index_counter = index_counter + 1;

                A(equation_counter, 18 + 24 * (j-1) + 24 * N * (i-1)) = 1;
                A(equation_counter, 6 + 24 * (j-1) + 24 * N * (i-0)) = -1;
                equation_counter = equation_counter + 1;
            end
        
            if j < N
                %i_j_E3 = i_(j+1)_W3

                I(index_counter) = equation_counter;
                J(index_counter) = 7 + 24 * (j-1) + 24 * N * (i-1);
                V(index_counter) = 1;
                index_counter = index_counter + 1;

                I(index_counter) = equation_counter;
                J(index_counter) = 19 + 24 * (j-0) + 24 * N * (i-1);
                V(index_counter) = -1;
                index_counter = index_counter + 1;

                A(equation_counter, 7 + 24 * (j-1) + 24 * N * (i-1)) = 1;
                A(equation_counter, 19 + 24 * (j-0) + 24 * N * (i-1)) = -1;
                equation_counter = equation_counter + 1;

                % ---

                I(index_counter) = equation_counter;
                J(index_counter) = 8 + 24 * (j-1) + 24 * N * (i-1);
                V(index_counter) = 1;
                index_counter = index_counter + 1;

                I(index_counter) = equation_counter;
                J(index_counter) = 20 + 24 * (j-0) + 24 * N * (i-1);
                V(index_counter) = -1;
                index_counter = index_counter + 1;

                A(equation_counter, 8 + 24 * (j-1) + 24 * N * (i-1)) = 1;
                A(equation_counter, 20 + 24 * (j-0) + 24 * N * (i-1)) = -1;
                equation_counter = equation_counter + 1;

                % ---

                I(index_counter) = equation_counter;
                J(index_counter) = 9 + 24 * (j-1) + 24 * N * (i-1);
                V(index_counter) = 1;
                index_counter = index_counter + 1;

                I(index_counter) = equation_counter;
                J(index_counter) = 21 + 24 * (j-0) + 24 * N * (i-1);
                V(index_counter) = -1;
                index_counter = index_counter + 1;

                A(equation_counter, 9 + 24 * (j-1) + 24 * N * (i-1)) = 1;
                A(equation_counter, 21 + 24 * (j-0) + 24 * N * (i-1)) = -1;
                equation_counter = equation_counter + 1;

                %i_j_E6 = i_(j+1)_W6
                I(index_counter) = equation_counter;
                J(index_counter) = 10 + 24 * (j-1) + 24 * N * (i-1);
                V(index_counter) = 1;
                index_counter = index_counter + 1;

                I(index_counter) = equation_counter;
                J(index_counter) = 22 + 24 * (j-0) + 24 * N * (i-1);
                V(index_counter) = -1;
                index_counter = index_counter + 1;

                A(equation_counter, 10 + 24 * (j-1) + 24 * N * (i-1)) = 1;
                A(equation_counter, 22 + 24 * (j-0) + 24 * N * (i-1)) = -1;
                equation_counter = equation_counter + 1;

                % ---

                I(index_counter) = equation_counter;
                J(index_counter) = 11 + 24 * (j-1) + 24 * N * (i-1);
                V(index_counter) = 1;
                index_counter = index_counter + 1;

                I(index_counter) = equation_counter;
                J(index_counter) = 23 + 24 * (j-0) + 24 * N * (i-1);
                V(index_counter) = -1;
                index_counter = index_counter + 1;

                A(equation_counter, 11 + 24 * (j-1) + 24 * N * (i-1)) = 1;
                A(equation_counter, 23 + 24 * (j-0) + 24 * N * (i-1)) = -1;
                equation_counter = equation_counter + 1;

                % ---

                I(index_counter) = equation_counter;
                J(index_counter) = 12 + 24 * (j-1) + 24 * N * (i-1);
                V(index_counter) = 1;
                index_counter = index_counter + 1;

                I(index_counter) = equation_counter;
                J(index_counter) = 24 + 24 * (j-0) + 24 * N * (i-1);
                V(index_counter) = -1;
                index_counter = index_counter + 1;

                A(equation_counter, 12 + 24 * (j-1) + 24 * N * (i-1)) = 1;
                A(equation_counter, 24 + 24 * (j-0) + 24 * N * (i-1)) = -1;
                equation_counter = equation_counter + 1;
            end
        end
    end

    I = I(1:index_counter-1);
    J = J(1:index_counter-1);
    V = V(1:index_counter-1);
    A_sparse = sparse(I, J, V);
end