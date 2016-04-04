function [yx] = CFD2(y, hx)
    % use 2nd order CFD
    yx = 1/(2*hx)*(y(:,3:end) - y(:,1:end-2));
    % use 2nd order FFD and BFD for endpoints
    yx = [-3/2*y(:,1) + 2*y(:,2) - 1/2*y(:,3), yx, 3/2*y(:,end) - 2*y(:,end-1) + 1/2*y(:,end-2)];
end
