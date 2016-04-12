function [yx] = CFD2(y, hx)
    % use 2nd order CFD
    yx = 1/(2*hx)*(y(:,3:end) - y(:,1:end-2));
    % Assume y is constant at endpoints, that is y(:,0) = y(:,1) and y(:,end+1) = y(:,end)
    yx = [1/(2*hx)*(y(:,2) - y(:,1)), yx, 1/(2*hx)*(y(:,end) - y(:,end-1))];
end
