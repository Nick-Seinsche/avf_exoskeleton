function omega_plot()
    % creates a surface plot of omega and a plot of different w^(1) for
    % different wall inclinations
    X = -0.8:0.007:0.8;
    wi = 0.2;
    Y1 = arrayfun(@(x) omega(x, 1, wi), X);
    Y2 = arrayfun(@(x) omega(x, 1, 0.43), X);
    %Y2 = arrayfun(@(x) -omega(-x, 2, wi), X);
    
    figure(1);
    plot(X,X); hold on
    plot(X, Y1);
    plot(X, Y2);
    %plot(X, -2 * abs(X));
    %plot(X, 0.5 * abs(X));

    plot(30.5123 / 360 * 2 * pi, omega(30.5123 / 360 * 2 * pi, 1, wi), "bo");

    legend("identity", "\omega^{(1)}(\alpha), \tau = 1,16",  "\omega^{(1)}(\alpha), \tau = 1,7396", "(\alpha, \omega^{(1)}(\alpha)), \tau = 1,16", "location", "northwest")

    figure(2);
    X = -pi/2:0.1:pi/2;
    [XX, YY] = meshgrid(X);
    OMEGA = zeros(size(XX));

    for i=1:size(XX, 1)
        for j=1:size(XX, 2)
            OMEGA(i,j) = omega([XX(i,j), YY(i,j)], 0, wi);
        end
    end
    
    surfc(XX, YY, OMEGA);


    %Z = p(1) * X.^2 + p(2) * X + p(3);
    %plot(X, 0.25 * X.^2); hold off
    % plot(X(2:end), dY);
end