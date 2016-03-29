function [q] = LaxWendroff(q, A, H, N, Nt, Nx, ht, hx, epsilon)
    for n = 1:Nt
        q(:,1,n+1)=q(:,1,n) - ht/(2*hx)*(A*q(:,2,n) - A*q(:,1,n)) + (ht/hx)^2/2*(A^2*(q(:,2,n) - q(:,1,n)));
        for i = 2:Nx-1
            q(:,i,n+1) = q(:,i,n) - ht/(2*hx)*(A*q(:,i+1,n) - A*q(:,i-1,n)) + (ht/hx)^2/2*(A^2*(q(:,i+1,n) - 2*q(:,i,n) + q(:,i-1,n)));
        end
        q(:,Nx,n+1)=q(:,Nx,n) - ht/(2*hx)*(A*q(:,Nx,n) - A*q(:,Nx-1,n)) + (ht/hx)^2/2*(A^2*(-q(:,Nx,n) + q(:,Nx-1,n)));

        % instant relaxation
        % w = H(p)
        q(N+2:end,:,n+1) = cell2mat(arrayfun(@(i) H(q(1:N+1,i,n)), 1:Nx, 'UniformOutput', false));
    end
end
