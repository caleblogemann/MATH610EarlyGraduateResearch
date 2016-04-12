function [q] = imex3(q, A, H, N, Nt, Nx, ht, hx, epsilon)
    p = q(1:N+1,:,:);
    w = q(N+2:end,:,:);

    alpha = A(N+2,1);

    a = 0.24169426078821;
    b = 0.06042356519705;
    eta = 0.1291528696059;

    c = epsilon/(epsilon - a*ht);
    d = a*ht/(epsilon - a*ht);
    for n = 1:Nt
        Hpn = cell2mat(arrayfun(@(i) H(p(:,i,n)), 1:Nx, 'UniformOutput', false));

        % first stage
        % p1 = p(:, :, n) = pn
        w1 = c*w(:,:,n) + d*Hpn;

        % second stage
        % p2 = p(:, :, n) = pn
        w2 = c*w(:,:,n) + d*w1;

        % third stage
        w2x = CFD4([w2(:,1), w2(:,1), w2, w2(:,end), w2(:, end)], hx);
        p3 = p(:,:,n) - ht*w2x;

        Hp3 = cell2mat(arrayfun(@(i) H(p3(:,i)), 1:Nx, 'UniformOutput', false));
        pnx = CFD4([p(:,1,n), p(:,1,n), p(:,:,n), p(:,end,n), p(:,end,n)], hx);
        w3 = c*(w(:,:,n) - ht*alpha*pnx) + d*((1/a - 1)*(Hpn - w2) + Hp3);

        % fourth stage
        w3x = CFD4([w3(:,1), w3(:,1), w3, w3(:,end), w3(:, end)], hx);
        p4 = p(:,:,n) - ht/4 *(w2x + w3x);

        Hp4 = cell2mat(arrayfun(@(i) H(p4(:,i)), 1:Nx, 'UniformOutput', false));
        p3x = CFD4([p3(:,1), p3(:,1), p3, p3(:,end), p3(:, end)], hx);
        w4 = c*(w(:,:,n) - alpha*ht/4*(pnx + p3x)) + d*(b/a*(Hpn - w1) + ...
            eta/a*(Hpn - w2) + (1/2 - b - eta - a)/a*(Hp3 - w3) + Hp4);

        % update
        w4x = CFD4([w4(:,1), w4(:,1), w4, w4(:,end), w4(:, end)], hx);
        p(:,:,n+1) = p(:,:,n) - ht/6*(w2x + w3x + 4*w4x);

        p4x = CFD4([p4(:,1), p4(:,1), p4, p4(:,end), p4(:, end)], hx);
        w(:,:,n+1) = w(:,:,n) - ht/6*alpha*(pnx + p3x + 4*p4x) + ...
            ht/(6*epsilon)*(Hpn + Hp3 - (w2 + w3) + 4*(Hp4 - w4));
    end

    q = [p; w];
end
