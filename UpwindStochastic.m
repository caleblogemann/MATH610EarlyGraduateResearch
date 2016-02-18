H = @(p) (1/2)*p.^2;
ul = 1;
ur = -1;
sigma = .1;
z = 0;

a = 0;
b = 1;
Nx = 200;
deltaX = (b-a)/Nx;
x = ((1:Nx)-(1/2))*deltaX;
%u0 = @(x, z) (x <= .5)*((ul + sigma*z)*x) + (x > .5)*((ur + sigma*z)*x + .5*(ul - ur));

% exact solution
pExact = @(x, z, t) (x <= .5 + sigma*z*t)*(ul + sigma*z) + ...
    (x > .5 + sigma*z*t)*(ur + sigma*z);
uExact = @(x, z, t) (x <= .5 + sigma*z*t)*((ul + sigma*z)*x ...
    -(ul + sigma*z)^2/2*t) + (x > .5 + sigma*z*t)*((ur + sigma*z)*x + ...
    .5*(ul - ur) - (ur + sigma*z)^2/2*t);

% two conditions for stability
% sqrt(alpha) > |H_p|
alpha = max(ul, ur)^2+.1;

% sqrt(alpha)*deltaT/deltaX < 1
% nu = sqrt(alpha)*deltaT/deltaX = .9
nu = .9;
deltaT = (nu*deltaX)/sqrt(alpha);
T = .3;
nTimeSteps = round(T/deltaT);

N = 7;
% set up initial conditions
q = zeros(2*(N+1), Nx, nTimeSteps+1);
q(1,:,1) = pExact(x, z, 0);
q(2,:,1) = H(q(1,:,1));

A = [zeros(N+1), eye(N+1); alpha*eye(N+1), zeros(N+1)];
delta = deltaT/deltaX;

%pn = @(i1, i2, i3, m, p, w) 1/(2*deltaX)*(sqrt(alpha)*deltaT*p(i1,m) + 2*(-sqrt(alpha)*deltaT + deltaX)*p(i2,m) + deltaT*(sqrt(alpha)*p(i3,m) + w(i1,m)- w(i3,m)));
%wn = @(i1, i2, i3, m, p, w) (1/(2*deltaX))*(sqrt(alpha)^2*deltaT*p(i1,m)-sqrt(alpha)^2*deltaT*p(i3,m) + sqrt(alpha)*deltaT*w(i1,m) - 2*sqrt(alpha)*deltaT*w(i2,m) + 2*deltaX*w(i2,m) + sqrt(alpha)*deltaT*w(i3,m));
% time stepping
for n = 1:nTimeSteps
    q(:,1,n+1)=q(:,1,n) - delta/2*(A*q(:,2,n) - A*q(:,1,n)) + delta^2/2*(A^2*(q(:,2,n) - q(:,1,n)));
    for i = 2:Nx-1
        q(:,i,n+1) = q(:,i,n) - delta/2*(A*q(:,i+1,n) - A*q(:,i-1,n)) + delta^2/2*(A^2*(q(:,i+1,n) - 2*q(:,i,n) + q(:,i-1,n)));
    end
    q(:,Nx,n+1)=q(:,Nx,n) - delta/2*(A*q(:,Nx,n) - A*q(:,Nx-1,n)) + delta^2/2*(A^2*(-q(:,Nx,n) + q(:,Nx-1,n)));
   
    % instant relaxation
    % w = H(p)
    q(N+2:end,:,n+1) = H(q(1:N+1,:,n+1));
end
% compute u at time T
u = cumsum(q(1,:,end)*deltaX);
figure(1);
plot(x,u);

% plot at time T
figure(2);
plot(x',q(1,:,end),'ro', 'Linewidth', 2);
hold on
plot(x,pExact(x,z,T),'Linewidth', 2);
%xlim([.3,.7]);

p = @(n, i, z) sum(q(1:N+1,i,n).*repmat(legendreP(0:N,z)',1,length(i)));