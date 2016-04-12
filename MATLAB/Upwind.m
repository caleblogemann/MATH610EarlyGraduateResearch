function [q] = Upwind(q, A, H, N, Nt, Nx, ht, hx, epsilon)
    [V, D] = eig(A);
    Dplus = max(zeros(size(D)), D);
    Dminus = min(zeros(size(D)), D);
    Aplus = V*Dplus*inv(V);
    Aminus = V*Dminus*inv(V);
    for n = 1:Nt
        q(:,1,n+1) = q(:,1,n) - ht/hx*(Aminus*(q(:,2,n) - q(:,1,n)));
        for i = 2:Nx-1
            q(:,i,n+1) = q(:,i,n) - ht/hx*(Aplus*(q(:,i,n) - q(:, i-1,n)) + Aminus*(q(:,i+1,n) - q(:,i,n)));
        end
        q(:,Nx,n+1) = q(:,Nx,n) - ht/hx*Aplus*(q(:,Nx,n) - q(:, Nx-1,n));

        % instant relaxation
        % w = H(p)
        q(N+2:end,:,n+1) = cell2mat(arrayfun(@(i) H(q(1:N+1,i,n)), 1:Nx, 'UniformOutput', false));
    end
end
