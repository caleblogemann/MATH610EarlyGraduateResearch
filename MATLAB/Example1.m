% specify problem
ul = 1;
ur = -1;
sigma = .1;
if ~exist('S')
    SImport
end

% set up parameters for method
alpha = 1.5;
N = 7;
Nx = 200;
CFL = .6;
epsilon = 1e-6;

% spacial discritization
a = 0;
b = 1;
deltaX = (b-a)/Nx;
x = ((1:Nx)-(1/2))*deltaX;

% set up exact solutions
pExact = @(x, z, t) (x <= .5 + sigma*z*t)*(ul + sigma*z) + ...
    (x > .5 + sigma*z*t)*(ur + sigma*z);
uExact = @(x, z, t) (x <= .5 + sigma*z*t)*((ul + sigma*z)*x ...
    -(ul + sigma*z)^2/2*t) + (x > .5 + sigma*z*t)*((ur + sigma*z)*x + ...
    .5*(ul - ur) - (ur + sigma*z)^2/2*t);

% get initial values of p
p0 = zeros(N+1,Nx);
p0(1,:) = pExact(x, 0, 0);
if N > 0
    p0(2,:) = sigma;
end

% set ending time
T = .001;
q = HamiltonJacobiInstantRelaxation(alpha, @(p)HExample1(p,S), N, Nx, deltaX, p0, T, CFL, epsilon, @imex3Matrix);

figure
subplot(2,2,1)
hold on
plot(x,deltaX*cumsum(q(1,:,end)), 'r')
xlim([.3,.7]);
title('Eu');
subplot(2,2,2)
plot(x,deltaX*cumsum(q(2,:,end)),'r')
title('Su');
xlim([.3,.7]);
hold on
subplot(2,2,3)
plot(x,q(1,:,end),'r')
title('Ep');
xlim([.3,.7]);
hold on
subplot(2,2,4)
plot(x,q(2,:,end),'r')
title('Sp')
xlim([.3,.7])
hold off
