% specify problem
ul = 1;
ur = -1;
sigma = .1;
%H = @(p) arrayfun(1/2, 1:length(p));

% spacial discritization
a = 0;
b = 1;
Nx = 200;
deltaX = (b-a)/Nx;
x = ((1:Nx)-(1/2))*deltaX;

% set up exact solutions
pExact = @(x, z, t) (x <= .5 + sigma*z*t)*(ul + sigma*z) + ...
    (x > .5 + sigma*z*t)*(ur + sigma*z);
uExact = @(x, z, t) (x <= .5 + sigma*z*t)*((ul + sigma*z)*x ...
    -(ul + sigma*z)^2/2*t) + (x > .5 + sigma*z*t)*((ur + sigma*z)*x + ...
    .5*(ul - ur) - (ur + sigma*z)^2/2*t);

% set up parameters for method
alpha = 1.5;
N = 7;

% get initial values of p
p0 = zeros(N+1,Nx);
p0(1,:) = pExact(x, 0, 0);
p0(2,:) = sigma;

% set ending time
T = .3;
q = HamiltonJacobiInstantRelaxation(alpha, @(p)H(p,S), N, Nx, deltaX, p0, T);
