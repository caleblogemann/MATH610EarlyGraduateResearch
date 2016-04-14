function [yx] = BFD2(y, hx)
    % use 2nd order BFD
    yx = 1/(2*hx)*(y(:,1:end-2) - 4*y(:,2:end-1) + 3*y(:,3:end));
end
