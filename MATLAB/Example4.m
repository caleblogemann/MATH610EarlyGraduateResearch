% specify problem
sigma = .2;

% set up parameters for method
alpha = 1.5;
N = 7;
Nx = 200;
CFL = .2;
epsilon = 1e-6;

% spacial discritization
a = 0;
b = 2*pi;
deltaX = (b-a)/Nx;
x = ((1:Nx)-(1/2))*deltaX;

% get initial values of p
p0 = zeros(N+1,Nx);
p0(1,:) = cos(x);

% set up Hamiltonian/Numerical integration for Hamiltonian
nGaussLegendrePoints = 20;
[z, w] = lgwt(nGaussLegendrePoints, -1, 1);
phi = @(z) sqrt((2*(0:N)+1)/2)'.*(legendreP(0:N, z))';
phiz = cell2mat(arrayfun(phi, z', 'UniformOutput', false));

% set ending time
T = 1;
q = HamiltonJacobiInstantRelaxation(alpha, @(p)HExample4(p,sigma, z, w, phiz), N, Nx, deltaX, p0, T, CFL, epsilon, @imex1);

figure
subplot(2,2,1)
hold on
plot(x,deltaX*cumsum(q(1,:,end)), 'r')
%xlim([.3,.7]);
title('Eu');
subplot(2,2,2)
plot(x,deltaX*cumsum(q(2,:,end)),'r')
title('Su');
%xlim([.3,.7]);
hold on
subplot(2,2,3)
plot(x,q(1,:,end),'r')
title('Ep');
%xlim([.3,.7]);
hold on
subplot(2,2,4)
plot(x,q(2,:,end),'r')
title('Sp')
%xlim([.3,.7])
hold off
