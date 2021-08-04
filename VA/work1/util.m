function [x] = eta(j)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    global L
    if j==1
        x = 1.875/L;
    end
    if j>1
        x = (j-0.5)/L;
    end

end

function [x] = D(j)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    global L
    x = cos(eta(j)*L);

end
