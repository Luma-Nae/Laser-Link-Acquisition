%% quaternion to matrix function

function Q = quat2matrix(q)

    Q = zeros(4,4); 
    
    qs = q(1);
    qv = q(2:4);
    
    Q(1,2:4) = -qv;
    Q(2:4, 1) = qv;
    
    S = diag([qs, qs, qs]);
    S(1,2) = q(4);
    S(1,3) = -q(3);
    S(2,1) = -q(4);
    S(3,1) = q(3);
    S(3,2) = -q(2);
    S(2,3) = q(2);
    Q(1,1) = qs;
    Q(2:4, 2:4) = S;
end