function [q] = HJIMEX(alpha, H, N, Nx, hx, p0, T, epsilon)

    % make sure p0 is correct size
    if(~isequal(size(p0),[N+1,Nx]))
        error('p0 is incorrect size');
    end

    % ht is the size of the time step
    % Nt is the number of time steps ht to get to end time T
    % determine ht in so make method stable
    % sqrt(alpha)*ht/hx < 1
    % nu = sqrt(alpha)*ht/hx = .9
    nu = .2;
    % ht must be smaller than htMax in order for nu < .9
    htMax = (nu*hx)/sqrt(alpha);
    % round Nt up
    Nt = floor(T/htMax + 1);
    % find exact deltaT to get to time=T in nTimesSteps time steps
    ht = T/Nt;

    % set up matrix to store all values of p and w
    p = zeros(N+1, Nx, Nt + 1);
    w = zeros(N+1,Nx,Nt + 1);
    % input initial conditions
    p(:,:,1) = p0;
    % find initial values of w based on H(p)
    w(:,:,1) = cell2mat(arrayfun(@(i) H(p(:,i,1)), 1:Nx, 'UniformOutput', false));

    % set up matrix to determine system
    % system is q_t + A*q_x = 0
    %A = [zeros(N+1), eye(N+1); alpha*eye(N+1), zeros(N+1)];

%     aTilde = zeros(3);
%     aTilde(3,2) = 1;
%     a = 1/2*eye(3);
%     a(2,1) = -1/2;
%     a(3,2) = 1/2;
%     wTilde = [0, 1/2, 1/2];
%     w = [0, 1/2, 1/2];
% 
%     C = eye((N+1)^2); 
%     C(1:N+1, 1:N+1) = 0;

    % time stepping
    for n = 1:Nt
        Hpn = cell2mat(arrayfun(@(i) H(p(:,i,n)), 1:Nx, 'UniformOutput', false));
        
        w1 = (w(:,:,n) + 1/epsilon * ht * Hpn)/(1 + ht/epsilon);
        
        p(:,:,n+1) = p(:,:,n) + ht*(1/hx)*[w1(:,2:Nx) - w1(:,1:Nx-1), w1(:,Nx) - w1(:,Nx-1)];
        w(:,:,n+1) = w(:,:,n) + ht*alpha*(1/hx)*[p(:,2:Nx,n) - p(:,1:Nx-1,n), p(:,Nx,n) - p(:,Nx-1,n)] + ...
            1/epsilon * ht * (Hpn - w1);
    end
    
    q = [p; w];
end
