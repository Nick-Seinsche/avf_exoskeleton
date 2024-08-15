function U0 = init_configuration(k)
    N = 2^(k+1);
    h = 2^(-k);
    L = N^2;
    U0 = zeros(12 * L, 1);
    for i = 1:N
        for j = 1:N
            U0((1:2) + 12 * (j-1) + 12 * N * (i-1)) = [-1 + (i-1)*h + h/2, -1 + (j-1)*h + h/2];
            U0((4:12) + 12 * (j-1) + 12 * N * (i-1)) = reshape(eye(3), 3, 3);
        end
    end
end