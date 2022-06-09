function T = RFrame(R)
   
    T = zeros(4);
    T(1:3,1:3) = R;
    T(4,4) = 1;
end

