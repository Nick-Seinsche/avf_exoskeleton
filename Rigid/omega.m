function x1 = omega(input, index, wi)
    if wi == 0
        wi = 0.87;
    end
    t11 = [wi; -wi; 1/2];

    t1 = t11(1);
    t2 = t11(2);
    t3 = t11(3);
    
    tau = t1^2/t3^2;

    if index == 1
        a = input;
        f = @(b) tau*((cos(a) + cos(b) - 2)/(tau + 1) - (sin(a) - sin(b))/(tau + 1)^(1/2))^2 + (cos(a) - (cos(b) - 1)/(tau + 1) + (sin(a) - sin(b))/(tau + 1)^(1/2))^2 + (cos(b) - (cos(a) - 1)/(tau + 1) + (sin(a) - sin(b))/(tau + 1)^(1/2))^2 - 2;
        df = @(b) 2*tau*(cos(b)/(tau + 1)^(1/2) - sin(b)/(tau + 1))*((cos(a) + cos(b) - 2)/(tau + 1) - (sin(a) - sin(b))/(tau + 1)^(1/2)) - 2*(sin(b) + cos(b)/(tau + 1)^(1/2))*(cos(b) - (cos(a) - 1)/(tau + 1) + (sin(a) - sin(b))/(tau + 1)^(1/2)) - 2*(cos(b)/(tau + 1)^(1/2) - sin(b)/(tau + 1))*(cos(a) - (cos(b) - 1)/(tau + 1) + (sin(a) - sin(b))/(tau + 1)^(1/2));

    elseif index == 2
        b = input;
        f = @(a) tau*((cos(a) + cos(b) - 2)/(tau + 1) - (sin(a) - sin(b))/(tau + 1)^(1/2))^2 + (cos(a) - (cos(b) - 1)/(tau + 1) + (sin(a) - sin(b))/(tau + 1)^(1/2))^2 + (cos(b) - (cos(a) - 1)/(tau + 1) + (sin(a) - sin(b))/(tau + 1)^(1/2))^2 - 2;
        df = @(a) 2*(cos(a)/(tau + 1)^(1/2) + sin(a)/(tau + 1))*(cos(b) - (cos(a) - 1)/(tau + 1) + (sin(a) - sin(b))/(tau + 1)^(1/2)) - 2*(sin(a) - cos(a)/(tau + 1)^(1/2))*(cos(a) - (cos(b) - 1)/(tau + 1) + (sin(a) - sin(b))/(tau + 1)^(1/2)) - 2*tau*(cos(a)/(tau + 1)^(1/2) + sin(a)/(tau + 1))*((cos(a) + cos(b) - 2)/(tau + 1) - (sin(a) - sin(b))/(tau + 1)^(1/2));
    else
        a = input(1);
        b = input(2);
        x1 = tau*((cos(a) + cos(b) - 2)/(tau + 1) - (sin(a) - sin(b))/(tau + 1)^(1/2))^2 + (cos(a) - (cos(b) - 1)/(tau + 1) + (sin(a) - sin(b))/(tau + 1)^(1/2))^2 + (cos(b) - (cos(a) - 1)/(tau + 1) + (sin(a) - sin(b))/(tau + 1)^(1/2))^2 - 2;
        return;
    end

    x1 = input * 0.8;
    while abs(f(x1)) > 1e-7
        x1 = x1 - f(x1) / df(x1);
    end
end