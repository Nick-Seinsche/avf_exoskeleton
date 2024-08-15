function energy_plot()
    X = linspace(0, 0.7, 20);
    
    e1 = zeros(size(X,2), 1);
    v1 = zeros(size(X,2), 1);

    e2 = zeros(size(X,2), 1);
    v2 = zeros(size(X,2), 1);

    e3 = zeros(size(X,2), 1);
    v3 = zeros(size(X,2), 1);

    e4 = zeros(size(X,2), 1);
    v4 = zeros(size(X,2), 1);
%{
    for i = 1:size(X, 2)
        [e1(i), v1(i)] = bending_model(1, X(i), 20.9, 0);
        [e2(i), v2(i)] = bending_model(1, X(i), 15.8, 0);
        [e3(i), v3(i)] = bending_model(1, X(i), 12.2, 0);
        [e4(i), v4(i)] = bending_model(1, X(i), 11.5, 0);
    end

    writematrix(e1, "cache/e1.txt");
    writematrix(e2, "cache/e2.txt");
    writematrix(e3, "cache/e3.txt");
    writematrix(e4, "cache/e4.txt");

    writematrix(v1, "cache/v1.txt");
    writematrix(v2, "cache/v2.txt");
    writematrix(v3, "cache/v3.txt");
    writematrix(v4, "cache/v4.txt");
%}
    e1 = readmatrix("cache/e1.txt");
    e2 = readmatrix("cache/e2.txt");
    e3 = readmatrix("cache/e3.txt");
    e4 = readmatrix("cache/e4.txt");

    v1 = readmatrix("cache/v1.txt");
    v2 = readmatrix("cache/v2.txt");
    v3 = readmatrix("cache/v3.txt");
    v4 = readmatrix("cache/v4.txt");

    figure(2);
    %yyaxis left
    plot(X, e1, "b"); hold on
    plot(X, e2, "b--");
    plot(X, e3, "g");
    plot(X, e4, "g--");
    ylabel("von Kárman Energy")
    %yyaxis right
    %semilogy(X, v1);
    %semilogy(X, v2);
    %semilogy(X, v3); 
    %semilogy(X, v4); hold off
    xlabel("\delta_z");
    %ylabel("Raw Constraint Violation");
    legend(["\alpha = 20.9", "\alpha = 15.8", "\alpha = 12.2", "\alpha = 11.5"], "Location", "northwest");

    %ax = axes; % or ax = gca();
    %ax.YAxis(1).Color = [0 0 0];
    %ax.YAxis(2).Color = [0 0 0];

end