function [tdir, Qdir] = transform3(dir, t11)
    % creates the transformation matricies for the affine functions
    % phi_ij_hat in the thesis, that transform the reference plate
    % to the plate in the configuration.

    % Qdir is the transformation does the transformation without the
    % wall inclination. applyinjg tdir afterwards applied the wall 
    % inclination to the plate
    
    % dir is the direction of the plate i.e N, S, W, E
    % t11 is the crease vector (NOT the modified crease vector; important!)
    
    % the transformation from the reference plate given by z (3x6)
    % to the plate in the configuration is:

    %plateN = repmat(xij, 1,6) + Qij * rot_mat(atan(t11(1)),0,0) * QdirN * z;
    %plateS = repmat(xij, 1,6) + Qij * rot_mat(atan(t11(2)),0,0) * QdirS * z;
    %plateW = repmat(xij, 1,6) + Qij * rot_mat(0,atan(t11(1)),0) * QdirW * z;
    %plateE = repmat(xij, 1,6) + Qij * rot_mat(0,atan(t11(2)),0) * QdirE * z;

    if dir == "N"
        Qdir = [0, -1, 0; 
                0, 0, -1; 
                1, 0, 0];
        tdir = rot_mat(atan(t11(1)),0,0);
    elseif dir == "S"
        Qdir = [0, 1, 0; 
                0, 0, 1; 
                1, 0, 0];
        tdir = rot_mat(atan(t11(2)),0,0);
    elseif dir == "W"
        Qdir = [0, 0, 1; 
                0, -1, 0; 
                -1, 0, 0];
        tdir = rot_mat(0,atan(t11(1)),0);
    elseif dir == "E"
        Qdir = [0, 0, -1; 
                0, 1, 0; 
                -1, 0, 0];
        tdir = rot_mat(0,atan(t11(2)),0);
    else
        Qdir = [];
        tdir = [];
    end
end

function Q = rot_mat(a,b,c)
    Q = [1, 0, 0; 0, cos(a), -sin(a); 0, sin(a), cos(a)];
    Q = Q * [cos(b), 0, sin(b); 0, 1, 0; -sin(b), 0, cos(b)];
    Q = Q * [cos(c), -sin(c), 0; sin(c), cos(c), 0; 0,0,1];
end
