function test_double_bending_nonlinear()    
    obj = zeros(4,1);
    for k = 0:3
        h = 2^(-k);
        h_T = h/2;
        dx = 0.3 / (2^k);
        alpha = 21.8;
        figure(k+1);
        [~, obj0] = bidirectional_bending_nonlinear(k, h_T, dx, alpha, 1); %k, delta_x, crease_vector
        obj(k + 1) = obj0;
    end
   figure(5);
   plot(obj)
end
