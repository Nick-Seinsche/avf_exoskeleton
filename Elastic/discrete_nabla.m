function [n1,n2, LAMBDA, toplength, bottomlength] = discrete_nabla(W, h, h_T, alpha, orientation, submatrixG, submatrixW)
    [LAMBDA, toplength, bottomlength] = discrete_nabla_matrix(h, h_T, alpha, orientation);
    LAMBDA = LAMBDA(submatrixG, submatrixW);
    n1 = LAMBDA * W;
    n2 = (W(6) - W(3)) / h_T;
end