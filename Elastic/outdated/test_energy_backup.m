function test_energy()
    k = -1;
    %UU = init_configuration(k);

    xij = zeros(3, 1);
    Qij = reshape (eye(3,3), 9, 1);
    
    sN1 = [-1.1;0];
    sN2 = [-1.1;0];
    uN = [0;0;0;0];
    
    sE1 = zeros(2,1);
    sE2 = zeros(2,1);
    uE = [0;0;0;0];
    
    sS1 = zeros(2,1);
    sS2 = zeros(2,1);
    uS = [0;0;0;0];

    sW1 = zeros(2,1);
    sW2 = zeros(2,1);
    uW = [0;0;0;0];

    UU = [xij; Qij; sN1; sN2; uN; sE1; sE2; uE; sS1; sS2; uS; sW1; sW2; uW];

    h = 2^(-k);
    h_T = 0.5 * h;
    N = 2 ^ (k + 1);
    alpha = 20.9;
    plot_skeleton_square_by_reference(UU, alpha, h, h_T);

    %[Splus, Bplus, correctorplus] = stiffness_and_load(h, h_T, alpha, 2);
    %[Sminus, Bminus, correctorminus] = stiffness_and_load(h, h_T, alpha, 1);
    
    [LAMBDA_plus, toplength_plus, bottomlength_plus] = discrete_nabla_matrix(h, h_T, alpha, 2);
    [LAMBDA_minus, toplength_minus, bottomlength_minus] = discrete_nabla_matrix(h, h_T, alpha, 1);

    %writematrix(Splus, "cache/S_plus.txt")
    %writematrix(Bplus, "cache/B_plus.txt")
    %writematrix(correctorplus, "cache/corrector_plus.txt")

    %writematrix(Sminus, "cache/S_minus.txt")
    %writematrix(Bminus, "cache/B_minus.txt")
    %writematrix(correctorminus, "cache/corrector_minus.txt")

    Splus = readmatrix("cache/S_plus.txt");
    Bplus = readmatrix("cache/B_plus.txt");
    correctorplus = readmatrix("cache/corrector_plus.txt");

    Sminus = readmatrix("cache/S_minus.txt");
    Bminus = readmatrix("cache/B_minus.txt");
    correctorminus = readmatrix("cache/corrector_minus.txt");

    submatrixW = [2,3,6,7];
    submatrixG = [2, 3, 5, 6];

    LAMBDA_plus = LAMBDA_plus(submatrixG, submatrixW);
    LAMBDA_minus = LAMBDA_minus(submatrixG, submatrixW);
    
    energy = 0;
    
    for i=1:N
        for j=1:N
            if mod(i+j,2) == 0
                S = Splus;
                B = Bplus;
                LAMBDA = LAMBDA_plus;
                corrector_inv = inv(correctorplus);
                toplength = toplength_plus;
                bottomlength = bottomlength_plus;
            elseif mod(i+j,2) == 1
                S = Sminus;
                B = Bminus;
                LAMBDA = LAMBDA_minus;
                corrector_inv = inv(correctorminus);
                toplength = toplength_minus;
                bottomlength = bottomlength_minus;
            end
            
            Uij = UU(indices(i,j,N));
            
            for p = [12, 20, 28, 36]
                wD = zeros(4, 1);

                % bending term ----------------------------------
                nabla_h_U = LAMBDA * Uij(p+5:p+8);
                energy = energy + nabla_h_U' * S * nabla_h_U;
                
                % sheer term -----------------------------------
                
                % wD(1) is w(2) in the thesis
                % the following is the gradient of p at z2
                nabla_p = corrector_inv * [nabla_h_U(1), (Uij(p+5) - Uij(p+7))]';
                nabla_s = corrector_inv * [Uij(p+1)/2/(toplength/2), Uij(p+1)/2 - Uij(p+2)/2; Uij(p+3)/2/(bottomlength/2), Uij(p+3)/2 - Uij(p+4)/2]';
    
                wD(1) = norm(nabla_p * nabla_p' + (nabla_s + nabla_s'), "fro");
    
                % wD(2) is w(3) in the thesis
                % the following is the gradient of p at z3
                nabla_p = [nabla_h_U(2), (Uij(p+6) - Uij(p+8))/h_T];
                nabla_s = [Uij(p+1)/toplength, Uij(p+1) - Uij(p+2); Uij(p+3)/bottomlength, Uij(p+3) - Uij(p+4)]';
    
                wD(2) = norm(nabla_p * nabla_p' + (nabla_s + nabla_s'), "fro");
    
                % wD(3) is w(5) in the thesis
                % the following is the gradient of p at z5
                nabla_p = corrector_inv * [nabla_h_U(3), (Uij(p+5) - Uij(p+7))]';
                nabla_s = corrector_inv * [Uij(p+2)/2/(toplength/2), Uij(p+1)/2 - Uij(p+2)/2; Uij(p+4)/2/(bottomlength/2), Uij(p+3)/2 - Uij(p+4)/2]';
    
                wD(3) = norm(nabla_p * nabla_p' + (nabla_s + nabla_s'), "fro");
    
                % wD(4) is w(6)
                % the followng is the gradient of p at z6
                nabla_p = [nabla_h_U(4), (Uij(p+6) - Uij(p+8))/h_T];
                nabla_s = corrector_inv * [Uij(p+2)/toplength, Uij(p+1) - Uij(p+2); Uij(p+4)/bottomlength, Uij(p+3) - Uij(p+4)]';
    
                wD(4) = norm(nabla_p * nabla_p' + (nabla_s + nabla_s'), "fro");
    
                energy = energy + B' * wD;

            end
        end
    end

    energy

end

function v = indices(i,j, N)
    d = (3 + 9 + 4 * (2 * 2 + 4 * 1));
    v = ((1:d) + d * (j - 1) + d * N * (i - 1));
end
