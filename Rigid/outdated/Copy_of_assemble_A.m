function A = Copy_of_assemble_A(k, t11)
    % Assembles the hinge condition matrix.
    % k - size of the configuration
    %
    h = 2^(-k);
    N = 2^(k+1);
    L = N^2;
    ctr = 1; 
    
    

    A = sparse(12 * N * (N-1), 12 * L);
    for i = 1:N
        for j = 1:N
            if i < N
                for k = 0:2
                    A(ctr+k, 1+12*(j-1)+(i-1) * 12 * N+k) = 1;
                    A(ctr+k, 1+12*(j-1)+(i-1) * 12 * N+3+k) = 0.5 * h;
                    A(ctr+k, 1+12*(j-1)+i * 12 * N+k) = -1;
                    A(ctr+k, 1+12*(j-1)+i * 12 * N+3+k) = 0.5 * h;
                end
                ctr = ctr + 3;
            end
            if j < N
                for k = 0:2
                    A(ctr+k, 1+12*(j-1)+(i-1) * 12 * N+k) = 1;
                    A(ctr+k, 1+12*(j-1)+(i-1) * 12 * N+6+k) = 0.5 * h;
                    A(ctr+k, 1+12*(j-1)+(i-1) * 12 * N+k+12) = -1;
                    A(ctr+k, 1+12*(j-1)+(i-1) * 12 * N+k+12+6) = 0.5 * h;
                end
                ctr = ctr + 3;
            end
            if i < N
                tij = vertcat((-1)^(i+j) * t11(1:2), t11(3));
                v = [0, -tij(2), tij(3)];

                for k = 0:2
                    A(ctr+k, 1+12*(j-1)+(i-1)* 12 * N+k+3) = v * [1;0;0];
                    A(ctr+k, 1+12*(j-1)+(i-1)* 12 * N+k+6) =  v * [0;1;0];
                    A(ctr+k, 1+12*(j-1)+(i-1) * 12 * N+k+9) = v * [0;0;1];
                    A(ctr+k, 1+12*(j-1)+i*12*N+k+3) = -v * [1; 0; 0];
                    A(ctr+k, 1+12*(j-1)+i*12*N+k+6) = -v * [0; 1; 0];
                    A(ctr+k, 1+12*(j-1)+i*12*N+k+9) = -v * [0; 0; 1];
                end
                ctr = ctr + 3;
            end
            if j < N
                tij = vertcat((-1)^(i+j) * t11(1:2), t11(3));
                v = [-tij(1), 0, tij(3)];

                for k = 0:2
                    A(ctr+k, 1+12*(j-1)+(i-1)*12*N+k+3) = v * [1;0;0];
                    A(ctr+k, 1+12*(j-1)+(i-1)*12*N+k+6) = v * [0;1;0];
                    A(ctr+k, 1+12*(j-1)+(i-1)*12*N+k+9) = v * [0;0;1];
                    A(ctr+k, 1+12*(j-1)+(i-1)*12*N+k+3+12) = -v * [1;0;0];
                    A(ctr+k, 1+12*(j-1)+(i-1)*12*N+k+6+12) = -v * [0;1;0];
                    A(ctr+k, 1+12*(j-1)+(i-1)*12*N+k+9+12) = -v * [0;0;1];
                end
                ctr = ctr + 3;
            end
        end
    end
end