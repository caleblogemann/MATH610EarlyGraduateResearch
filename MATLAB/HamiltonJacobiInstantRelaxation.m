function [q] = HamiltonJacobiInstantRelaxation(alpha, H, N, Nx, hx, p0, T, CFL, epsilon, numericalMethod)
    % Input checking/parsing requires my personal utility functions
    %p = inputParser();
    %p.addRequired('alpha', @Utils.isPositiveNumber);
    %p.addRequired('H', @Utils.isFunctionHandle);
    %p.addRequired('N', @Utils.isPositiveInteger);
    %p.addRequired('Nx', @Utils.isPositiveInteger);
    %p.addRequired('hx', @Utils.isPositiveNumber);
    %p.addRequired('p0', @Utils.isMatrix);
    %p.addRequired('T', @Utils.isPositiveNumber);
    %p.parse(alpha, H, N, Nx, hx, p0, T);

    % make sure p0 is correct size
    if(~isequal(size(p0),[N+1,Nx]))
        error('p0 is incorrect size');
    end

    % ht is the size of the time step
    % Nt is the number of time steps ht to get to end time T
    % determine ht in so make method stable
    % CFL = sqrt(alpha)*ht/hx < 1
    % ht must be smaller than htMax in order for nu condition to be met
    htMax = (CFL*hx)/sqrt(alpha);
    % round Nt up
    Nt = floor(T/htMax + 1);
    % find exact deltaT to get to time=T in nTimesSteps time steps
    ht = T/Nt;

    % set up matrix to store all values of p and w
    q = zeros(2*(N+1), Nx, Nt + 1);
    % input initial conditions
    q(1:N+1,:,1) = p0;
    % find initial values of w based on H(p)
    q(N+2:end,:,1) = cell2mat(arrayfun(@(i) H(q(1:N+1,i,1)), 1:Nx, 'UniformOutput', false));

    % set up matrix to determine system
    % system is q_t + A*q_x = 0
    A = [zeros(N+1), eye(N+1); alpha*eye(N+1), zeros(N+1)];

    q = numericalMethod(q, A, H, N, Nt, Nx, ht, hx, epsilon);
    % time stepping
    %for n = 1:Nt
    %    q(:,1,n+1)=q(:,1,n) - ht/(2*hx)*(A*q(:,2,n) - A*q(:,1,n)) + (ht/hx)^2/2*(A^2*(q(:,2,n) - q(:,1,n)));
    %    for i = 2:Nx-1
    %        q(:,i,n+1) = q(:,i,n) - ht/(2*hx)*(A*q(:,i+1,n) - A*q(:,i-1,n)) + (ht/hx)^2/2*(A^2*(q(:,i+1,n) - 2*q(:,i,n) + q(:,i-1,n)));
    %    end
    %    q(:,Nx,n+1)=q(:,Nx,n) - ht/(2*hx)*(A*q(:,Nx,n) - A*q(:,Nx-1,n)) + (ht/hx)^2/2*(A^2*(-q(:,Nx,n) + q(:,Nx-1,n)));

    %    % instant relaxation
    %    % w = H(p)
    %    q(N+2:end,:,n+1) = cell2mat(arrayfun(@(i) H(q(1:N+1,i,n)), 1:Nx, 'UniformOutput', false));
    %end
end
