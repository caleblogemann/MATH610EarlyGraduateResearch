function [yx] = FFD2(y, hx)
    % use 2nd order FFD
    yx = 1/(2*hx)*(-y(:,3:end) + 4*y(:,2:end-1) - 3*y(:,1:end-2));
end
