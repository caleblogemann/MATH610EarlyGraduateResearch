function [q] = imex2Vector(q, A, H, N, Nt, Nx, ht, hx, epsilon)
    p = q(1:N+1,:,:);
    w = q(N+2:end,:,:);

    alpha = A(N+2,1);

    for n = 1:Nt
        Hpn = cell2mat(arrayfun(@(i) H(p(:,i,n)), 1:Nx, 'UniformOutput', false));

        w1 = (w(:,:,n) + 1/epsilon * ht * Hpn)/(1 + ht/epsilon);

        p(:,:,n+1) = p(:,:,n) + ht*(1/hx)*[w1(:,2:Nx) - w1(:,1:Nx-1), w1(:,Nx) - w1(:,Nx-1)];
        w(:,:,n+1) = w(:,:,n) + ht*alpha*(1/hx)*[p(:,2:Nx,n) - p(:,1:Nx-1,n), p(:,Nx,n) - p(:,Nx-1,n)] + ...
            1/epsilon * ht * (Hpn - w1);
    end

    q = [p; w];
end
