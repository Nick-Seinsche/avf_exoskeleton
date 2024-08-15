function [p, Q] = calculate_neighbor(h, t11, orientation, dir, p0, Q0, alpha)
    % Given
    % Problem size - h,
    % crease vector - t11
    % orientation of the cross - orientation
    % direction (n,s,w,e) - dir
    % position - p0
    % rotation matrix - Q0
    % wall inclination - alpha

    % Calculates for a given cross (p0, Q0) the neighbor (p,Q) in direction dir
    % such that the hinge condition is satisfied and the hinge angle is
    % alpha
    tij = vertcat((-1)^(orientation) * t11(1:2), t11(3));
    if dir == "s"
        Qr = rot_mat(atan(-tij(2)/tij(3)),0,0);
        p_hat = Q0 * [h/2; 0; 0] + (p0 + Q0 * [h/2; 0; 0]);
        %Q_hat = Q0;

        p_hat_hat = p0 + Q0' * (p_hat - p0);
        %Q_hat_hat = Qr';

        p_hat_hat_hat = Qr * rot_mat(0,0,alpha) * Qr' * [h/2; 0; 0] + (p_hat_hat - Qr' * [h/2; 0; 0]);
        Q_hat_hat_hat = Qr * rot_mat(0,0,alpha) * Qr';

        p = p0 + Q0 * (p_hat_hat_hat - p0);
        Q = Q0 * Q_hat_hat_hat;
    elseif dir == "n"
        Qr = rot_mat(atan(tij(2)/tij(3)),0,0);
        p_hat = Q0 * [-h/2; 0; 0] + (p0 + Q0 * [-h/2; 0; 0]);
        %Q_hat = Q0;

        p_hat_hat = p0 + Q0' * (p_hat - p0);
        %Q_hat_hat = Qr';

        p_hat_hat_hat = Qr * rot_mat(0,0,alpha) * Qr' * [-h/2; 0; 0] + (p_hat_hat - Qr' * [-h/2; 0; 0]);
        Q_hat_hat_hat = Qr * rot_mat(0,0,alpha) * Qr';

        p = p0 + Q0 * (p_hat_hat_hat - p0);
        Q = Q0 * Q_hat_hat_hat;
    elseif dir == "w"
        Qr = rot_mat(0,atan(-tij(1)/tij(3)),0);
        p_hat = Q0 * [0; -h/2; 0] + (p0 + Q0 * [0; -h/2; 0]);
        %Q_hat = Q0;

        p_hat_hat = p0 + Q0' * (p_hat - p0);
        %Q_hat_hat = Qr';

        p_hat_hat_hat = Qr * rot_mat(0,0,alpha) * Qr' * [0; -h/2; 0] + (p_hat_hat - Qr' * [0; -h/2; 0]);
        Q_hat_hat_hat = Qr * rot_mat(0,0,alpha) * Qr';

        p = p0 + Q0 * (p_hat_hat_hat - p0);
        Q = Q0 * Q_hat_hat_hat;
    elseif dir == "e"
        Qr = rot_mat(0,atan(tij(1)/tij(3)),0); %-
        p_hat = Q0 * [0; h/2; 0] + (p0 + Q0 * [0; h/2; 0]);
        %Q_hat = Q0;

        p_hat_hat = p0 + Q0' * (p_hat - p0);
        %Q_hat_hat = Qr';

        p_hat_hat_hat = Qr * rot_mat(0,0,alpha) * Qr' * [0; h/2; 0] + (p_hat_hat - Qr' * [0; h/2; 0]);
        Q_hat_hat_hat = Qr * rot_mat(0,0,alpha) * Qr';

        p = p0 + Q0 * (p_hat_hat_hat - p0);
        Q = Q0 * Q_hat_hat_hat;
    else
        p = 0;
        Q = 0;
    end
end

function Q = rot_mat(a,b,c)
    Q = [1, 0, 0; 0, cos(a), -sin(a); 0, sin(a), cos(a)];
    Q = Q * [cos(b), 0, sin(b); 0, 1, 0; -sin(b), 0, cos(b)];
    Q = Q * [cos(c), -sin(c), 0; sin(c), cos(c), 0; 0,0,1];
end