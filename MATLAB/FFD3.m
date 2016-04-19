function [yx] = FFD3(y, hx)
    % use 2nd order FFD
    yx = 1/(2*hx)*(-11*y(:,1:end-3) + 18*y(:,2:end-2) - 9*y(:,3:end-1) + 2*y(:,4:end));
end
