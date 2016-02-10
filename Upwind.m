a = 0;
b = 1;
N = 1000;
deltaX = (b-a)/N;
deltaT = .005;
lambda = .01;
nTimeSteps = round(.1/deltaT);
disp(lambda*deltaT/deltaX);

H = @(p) (1/2)*p.^2;

ul = 2;
ur = 1;

x = ((1:N)-(1/2))*deltaX;

% set up initial conditions
p = zeros(N, nTimeSteps+1);
w = zeros(N, nTimeSteps+1);
for i = 1:N
    if(x(i) < .5)
        p(i,1) = ul;
    else
        p(i,1) = ur;
    end
    w(i,1) = H(p(i));
end

pn = @(i1, i2, i3, m, p, w) 1/(2*deltaX)*(lambda*deltaT*p(i1,m) + 2*(-lambda*deltaT + deltaX)*p(i2,m) + deltaT*(lambda*p(i3,m) + w(i1,m)- w(i3,m)));
%wn = @(i1, i2, i3, m, p, w) (1/(2*deltaX))*(lambda^2*deltaT*p(i1,m)-lambda^2*deltaT*p(i3,m) + lambda*deltaT*w(i1,m) - 2*lambda*deltaT*w(i2,m) + 2*deltaX*w(i2,m) + lambda*deltaT*w(i3,m));
% time stepping
for n = 1:nTimeSteps
    p(1,n+1)=pn(1,1,2,n, p, w);
    %w(1,n+1)=wn(1,1,2,n, p, w);
    for i = 2:N-1
        p(i,n+1)=pn(i-1,i,i+1,n, p, w);
        %w(i,n+1)=wn(i-1,i,i+1,n, p, w);
    end
   p(N,n+1)=pn(N-1,N, N, n, p, w);
   %w(N,n+1)=wn(N-1,N, N, n, p, w);
   
   % instant relaxation
   %for j = 1:N
   %    w(j,n+1) = H(p(j,n+1));
   %end
   w(:,n+1) = H(p(:,n+1));
end

% plot each time step
for n=1:nTimeSteps+1
    plot(x,p(:,n));
    pause
end
