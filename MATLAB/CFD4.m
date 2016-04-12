function [yx] = CFD4(y, hx)
    % use 4th order CFD
    yx = 1/(12*hx)*(-y(:,5:end) + 8*y(:,4:end-1) - 8*y(:,2:end-3) +  y(:,1:end-4));
end
