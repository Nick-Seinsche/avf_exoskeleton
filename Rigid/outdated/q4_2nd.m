function q4_2nd
    a = sym('a', 'real');
    %a = 0;
    b = sym('b', 'real');
    t1 = sym('t1', 'real');
    t2 = -t1;%sym('t2', 'real');
    t3 = sym('t3', 'real');
    t = [t1, t2, t3];
    q = sym('tau', 'real');
    
    tau = cos(atan(t(1) / t(3)));
    tau_t = sin(atan(t(1) / t(3)));

    j = cos(atan(-t(2)/t(3)));
    j_t = -(t(2)/t(3)) * j;

    B = [tau, 0, tau_t; 0, 1, 0; -tau_t, 0, tau];
    C = [1, 0, 0; 0, j, -j_t; 0, j_t, j];

    Q1 = eye(3);
    Q2 = B' * [cos(a), -sin(a), 0; sin(a), cos(a), 0; 0, 0, 1] * B;
    Q3 = C' * [cos(b), -sin(b), 0; sin(b), cos(b), 0; 0, 0, 1] * C;

    w = Q3 * [t(1); 0; t(3)];
    w_alt = Q2 * [0; t(2); t(3)];

    i = cos(atan(-w(1) / w(3)));
    i_t = sin(atan(-w(1) / w(3)));

    i_alt = cos(atan(-w_alt(1) / w_alt(3)));
    i_t_alt = sin(atan(-w_alt(1) / w_alt(3)));

    G = [i, 0, i_t; 0, 1, 0; -i_t, 0, i];

    G_alt = [i_alt, 0, i_t_alt; 0, 1, 0; -i_t_alt, 0, i_alt];
    
    v = G * w;

    v_alt = G_alt * w_alt;

    k = cos(atan(v(2)/v(3)));
    k_t = sin(atan(v(2)/v(3)));

    k_alt = cos(atan(v_alt(2)/v_alt(3)));
    k_t_alt = sin(atan(v_alt(2)/v_alt(3)));
    
    H = [1, 0, 0; 0, k, -k_t; 0, k_t, k];
    
    H_alt = [1, 0, 0; 0, k_alt, -k_t_alt; 0, k_t_alt, k_alt];
    
    D = H * G;

    D_alt = H_alt * G_alt;

    Q4 = D' * [cos(a), sin(a), 0; -sin(a), cos(a), 0; 0,0,1] * D * Q3;
    
    Q4_alt = D_alt' * [cos(b), sin(b), 0; -sin(b), cos(b), 0; 0,0,1] * D_alt * Q2;

    %f0 = subs(Q4(2,1) - Q4_alt(2,1), a, 0)
    %df0 = diff(Q4(2,1) - Q4_alt(2,1), a);
    %ddf0 = diff(df0, a);
    %subs(df0, a, 0)
    %subs(ddf0, a, 0)

    U0 = init_configuration(0);

    p1 = U0(1:3);
    %Q1 = reshape(U0(4:12), 3, 3);
    p2 = U0(13:15);
    %Q2 = reshape(U0(16:24), 3, 3);
    p3 = U0(25:27);
    %Q3 = reshape(U0(28:36), 3, 3);
    p4 = U0(37:39);
    %Q4 = reshape(U0(40:48), 3, 3);

    mp12 = 0.5 * (p1 + p2);
    mp13 = 0.5 * (p1 + p3);

    p2_new = Q2 * (p2 - mp12) + mp12;
    p3_new = Q3 * (p3 - mp13) + mp13;
    p4_new = Q3 * (p4 - mp13) + mp13;
    p4_new_alt = Q2 * (p4 - mp12) + mp12;

    mp34 = 0.5 * (p3_new + p4_new);
    mp24 = 0.5 * (p2_new + p4_new_alt);
    
    %f = norm(mp24 - mp34, 2)^2 - 0.5;
    f = (mp24(1) - mp34(1))^2 + (mp24(2) - mp34(2))^2 + (mp24(3) - mp34(3))^2 -0.5;
    dfb = diff(f, b);
    dfa = diff(f, a);
    latex(f);

    %latex(df)

    express = q * ((2 - cos(b) - cos(a))/(q+1) + (sin(a) - sin(b))/(sqrt(q + 1)))^2 + (cos(b) + (1 - cos(a))/(q + 1)+(sin(a) - sin(b))/(sqrt(1 + q)))^2 + (cos(a) + (1 - cos(b))/(q + 1)+(sin(a) - sin(b))/(sqrt(1 + q)))^2 - 2;
    dfb = diff(express, b);
    dfa = diff(express, a);
    dfba = diff(dfb, a);
    dfaa = diff(dfa, a);
    dfbb = diff(dfb, b);

    subs(subs(dfa, a, 0), b, 0)
    subs(subs(dfb, a, 0), b, 0)
    subs(subs(dfaa, a, 0), b, 0)
    subs(subs(dfba, a, 0), b, 0)
    subs(subs(dfbb, a, 0), b, 0)


end