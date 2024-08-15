 function A = assemble_A3D(k)
    N = 2^(k+1);
    L = N^2;
    d = (12 + 4 * (2 * 2 + 4 * 1));
    A = sparse(12 * N * (N - 1), d * L);
    ctr = 1;

    for i = 1:N
        for j = 1:N
            if i < N
                % sS1,i,j,1 = sN1,i+1,j,1
                % sS1,i,j,2 = -sN1,i+1,j,2
                
                % sS2,i,j,1 = sN2,i+1,j,1
                % sS2,i,j,2 = -sN2,i+1,j,2

                % uS,i,j,2 = -uN,i+1,j,2
                % uS,i,j,4 = -uN,i+1,j,4

                % sS1,i,j = U((29:30) + d * (j-1) + d * N * (i-1))
                % sN1,i+1,j = U((13:14) + d * (j-1) + d * N * (i-0))
                
                % sS2,i,j = U((31:32) + d * (j-1) + d * N * (i-1))
                % sN2,i+1,j = U((15:16) + d * (j-1) + d * N * (i-0))
                
                % uS,i,j = U((33:36) + d * (j-1) + d * N * (i-1))
                % uN,i+1,j = U((17:20) + d * (j-1) + d * N * (i-0))

                A(ctr, 29 + d * (j-1) + d * N * (i-1)) = 1;
                A(ctr, 13 + d * (j-1) + d * N * (i-0)) = -1;
                ctr = ctr + 1;

                A(ctr, 30 + d * (j-1) + d * N * (i-1)) = 1;
                A(ctr, 14 + d * (j-1) + d * N * (i-0)) = +1;
                ctr = ctr + 1;

                A(ctr, 31 + d * (j-1) + d * N * (i-1)) = 1;
                A(ctr, 15 + d * (j-1) + d * N * (i-0)) = -1;
                ctr = ctr + 1;

                A(ctr, 32 + d * (j-1) + d * N * (i-1)) = 1;
                A(ctr, 16 + d * (j-1) + d * N * (i-0)) = +1;
                ctr = ctr + 1;
                
                A(ctr, 34 + d * (j-1) + d * N * (i-1)) = 1;
                A(ctr, 18 + d * (j-1) + d * N * (i-0)) = +1;
                ctr = ctr + 1;

                A(ctr, 36 + d * (j-1) + d * N * (i-1)) = 1;
                A(ctr, 20 + d * (j-1) + d * N * (i-0)) = +1;
                ctr = ctr + 1;
            end

            if j < N
                % sE1,i,j,1 = sW1,i,j+1,1
                % sE1,i,j,2 = -sW1,i,j+1,2
                
                % sE2,i,j,1 = sW2,i,j+1,1
                % sE2,i,j,2 = -sW2,i,j+1,2

                % uE,i,j,2 = -uW,i,j+1,2
                % uE,i,j,4 = -uW,i,j+1,4
                
                % sE1,i,j = U((21:22) + d * (j-1) + d * N * (i-1))
                % sW1,i,j+1 = U((37:38) + d * (j-0) + d * N * (i-1))

                % sE2,i,j = U((23:24) + d * (j-1) + d * N * (i-1))
                % sW2,i,j+1 = U((39:40) + d * (j-0) + d * N * (i-1))

                % uE,i,j = U((25:28) + d * (j-1) + d * N * (i-1))
                % uW,i,j+1 = U((41:44) + d * (j-0) + d * N * (i-1))

                A(ctr, 21 + d * (j-1) + d * N * (i-1)) = 1;
                A(ctr, 37 + d * (j-0) + d * N * (i-1)) = -1;
                ctr = ctr + 1;

                A(ctr, 22 + d * (j-1) + d * N * (i-1)) = 1;
                A(ctr, 38 + d * (j-0) + d * N * (i-1)) = +1;
                ctr = ctr + 1;

                A(ctr, 23 + d * (j-1) + d * N * (i-1)) = 1;
                A(ctr, 39 + d * (j-0) + d * N * (i-1)) = -1;
                ctr = ctr + 1;

                A(ctr, 24 + d * (j-1) + d * N * (i-1)) = 1;
                A(ctr, 40 + d * (j-0) + d * N * (i-1)) = +1;
                ctr = ctr + 1;


                A(ctr, 26 + d * (j-1) + d * N * (i-1)) = 1;
                A(ctr, 42 + d * (j-0) + d * N * (i-1)) = +1;
                ctr = ctr + 1;

                A(ctr, 28 + d * (j-1) + d * N * (i-1)) = 1;
                A(ctr, 44 + d * (j-0) + d * N * (i-1)) = +1;
                ctr = ctr + 1;
            end
        end
    end

    ctr
    size(A)
end