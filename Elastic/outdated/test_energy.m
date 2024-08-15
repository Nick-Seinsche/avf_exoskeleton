function test_energy()
    k = -1;
    L = 1;
    d = (3 + 9 + 4 * (2 * 2 + 4 * 1)); % DOF's per cross
    %UU = init_configuration(k);

    addpath('../casadi')
    import casadi.*

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

    U = [xij; Qij; sN1; sN2; uN; sE1; sE2; uE; sS1; sS2; uS; sW1; sW2; uW];

    h = 2^(-k);
    h_T = 0.5 * h;
    N = 2 ^ (k + 1);
    alpha = 20.9;
    plot_skeleton_square_by_reference(U, alpha, h, h_T);

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

    UU = SX.sym('x', d * L);
    
    energyBend = SX.zeros(N,N);
    energySheer = SX.zeros(N,N);
    e_gamma = 1;
    
    % The energy is determined in two steps:
    % First the matricies are chosen with the correct orientation
    % then loops over north west east south blade and sums up
    % the energy for one blade
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
                % pure bending term ----------------------------------
                nabla_h_U = LAMBDA * Uij(p+5:p+8);
                energyBend(i,j) = energyBend(i,j) + e_gamma / 2 * nabla_h_U' * S * nabla_h_U;

                % sheer term -----------------------------------
                % coefficient vector for the sheer term
                ZD = SX.zeros(4, 1);

                % ZD(1) is Z_2 in the thesis
                
                % the following is the gradient of p at z2.
                % By the construction of the discrete, nabla,
                % the derivative in x direction at the top and
                % bottom edge of the plate is exact. Therefore
                % we may use nabla_h_U(1) for nabla_x w(z_{2,2})
                nabla_w = corrector_inv * [nabla_h_U(1), (Uij(p+5) - Uij(p+7))]';
                
                % the following is the gradient of s at z2
                % Uij(p+1)/2/(toplength/2) is the derivative of s1 at
                % z2 in x direction, where we use linearity to calculate
                % the constant slope.
                % Uij(p+1)/2 - Uij(p+3)/2 is just s(z_{2}) - w(z_{5})
                % in the thesis
                nabla_v = corrector_inv * [Uij(p+1)/2/(toplength/2), Uij(p+2)/2/(toplength/2); Uij(p+1)/2 - Uij(p+3)/2, Uij(p+2)/2 - Uij(p+4)/2]';
                
                % calculate the sheer term from the nablas
                sheer = nabla_w * nabla_w' + (nabla_v + nabla_v');

                % calculate Z_2 as in the thesis
                ZD(1) = sheer(1)^2 + sheer(2)^2 + sheer(3)^2 + sheer(4)^2;

                % wD(2) is w(3) in the thesis
                % the following is the gradient of p at z3
                % nabla_h_U(2) is nabla_x w(z_{3,2})
                % since the corrector matrix is diagonal, we just need
                % to divide the second entry by h_T
                nabla_w = [nabla_h_U(2), (Uij(p+6) - Uij(p+8))/h_T];

                % the following is the gradient of s at z_3
                % again the corrector matrix is diagonal, we just need
                % to divide by h_T
                nabla_v = [Uij(p+1)/toplength, Uij(p+2)/toplength; (Uij(p+1) - Uij(p+3))/h_T, (Uij(p+2) - Uij(p+4))/h_T]';
    
                sheer = nabla_w * nabla_w' + (nabla_v + nabla_v');
                ZD(2) = sheer(1)^2 + sheer(2)^2 + sheer(3)^2 + sheer(4)^2;

                % wD(3) is w(5) in the thesis
                % the following is the gradient of p at z5
                nabla_w = corrector_inv * [nabla_h_U(3), (Uij(p+5) - Uij(p+7))]';

                nabla_v = corrector_inv * [Uij(p+3)/2/(bottomlength/2), Uij(p+4)/2/(bottomlength/2); Uij(p+1)/2 - Uij(p+3)/2, Uij(p+2)/2 - Uij(p+4)/2]';
    
                sheer = nabla_w * nabla_w' + (nabla_v + nabla_v');
                ZD(3) = sheer(1)^2 + sheer(2)^2 + sheer(3)^2 + sheer(4)^2;
    
                % wD(4) is w(6)
                % the followng is the gradient of p at z6
                nabla_w = [nabla_h_U(4), (Uij(p+6) - Uij(p+8))/h_T];

                nabla_v = [Uij(p+3)/bottomlength, Uij(p+4)/bottomlength; (Uij(p+1) - Uij(p+3))/h_T, (Uij(p+2) - Uij(p+4))/h_T]';
    
                sheer = nabla_w * nabla_w' + (nabla_v + nabla_v');
                ZD(4) = sheer(1)^2 + sheer(2)^2 + sheer(3)^2 + sheer(4)^2;

                energySheer(i,j) = energySheer(i,j) + 1 * 0.5 * B' * ZD;

            end
        end
    end

    energy = sum(sum(energyBend)) + sum(sum(energySheer));
    
    measureEnergy = Function('m', {UU}, {energyBend + energySheer}, {'U'}, {'e'});
    measureEnergyBend = Function('m', {UU}, {energyBend}, {'U'}, {'e'});
    measureEnergySheer = Function('m', {UU}, {energySheer}, {'U'}, {'e'});

    measureEnergy(U)
    measureEnergyBend(U)
    measureEnergySheer(U)

end

function v = indices(i,j, N)
    d = (3 + 9 + 4 * (2 * 2 + 4 * 1));
    v = ((1:d) + d * (j - 1) + d * N * (i - 1));
end
