function R = q2R(v)

    I = eye(3);
    
    qs = v(1);
    qv = v(2:4);
    
    A = diag([0,0,0]);
    A(2,1) = -qv(3);
    A(3,1)  = qv(2);
    A(1,2) = qv(3);
    A(1,3) = -qv(2);
    A(3,2) = -qv(1);
    A(2,3) = qv(1);

    R = I + 2*qs*A + 2*A^2;

end
