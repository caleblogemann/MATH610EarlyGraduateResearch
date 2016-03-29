function [q] = imex1Matrix(q, A, H, N, Nt, Nx, ht, hx, epsilon)
    for n = 1:Nt
        Hpn = cell2mat(arrayfun(@(i) H(q(1:N+1,i,n)), 1:Nx, 'UniformOutput', false));

        q1 = [q(1:N+1,:,n); (q(N+2:end,:,n) + ht/epsilon * Hpn)/(1 + ht/epsilon)];

        q1x = (1/hx) * diff(q1, 1, 2);
        q1x = [q1x, q1x(:,end)];
        q(:,:,n+1) = q(:, :, n) + ht * A * q1x + ht/epsilon * [zeros(N+1, Nx); Hpn - q1(N+2:end,:)];
    end
end
