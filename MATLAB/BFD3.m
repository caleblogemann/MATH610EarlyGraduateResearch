function [yx] = BFD3(y, hx)
    % use 2nd order BFD
    yx = 1/(6*hx)*(-2*y(:,1:end-3) + 9*y(:,2:end-2) - 18*y(:,3:end-1) + 11*y(:,4:end));
end
