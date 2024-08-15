function [plateN, plateE, plateS, plateW] = reference_cross(xij, Qij, alpha)
    alpha = 20.9;
    alpha_rad = alpha / 360 * 2 * pi;
    t11 = [-tan(alpha_rad); tan(alpha_rad); 1];
    t11_v = rot_mat(atan(t11(2)),0,0) * t11;

    z1 = [h_T/2, -t11_v(1)/t11_v(3) * h_T/2];
    z4 = -z1;
    z3 = [h_T/2, h/2];
    z6 = [-h_T/2, h/2];
    z2 = 0.5 * (z3 + z1);
    z5 = 0.5 * (z4 + z6);

    z = [z1', z2', z3', z4', z5', z6'];

    QdirN = [0, -1; 0, 0; 1, 0];
    QdirS = [0, 1; 0, 0; 1, 0];
    QdirW = [0, 0; 0, -1; -1, 0];
    QdirE = [0, 0; 0, 1; -1, 0];

    plateN = repmat(xij, 1,6) + Qij * rot_mat(atan(t11(1)),0,0) * QdirN * z;
    plateS = repmat(xij, 1,6) + Qij * rot_mat(atan(t11(2)),0,0) * QdirS * z;
    plateW = repmat(xij, 1,6) + Qij * rot_mat(0,atan(t11(1)),0) * QdirW * z;
    plateE = repmat(xij, 1,6) + Qij * rot_mat(0,atan(t11(2)),0) * QdirE * z;
end