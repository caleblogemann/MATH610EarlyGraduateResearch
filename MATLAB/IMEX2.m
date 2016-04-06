function [q] = imex2(q, A, H, N, Nt, Nx, ht, hx, epsilon)
    p = q(1:N+1,:,:);
    w = q(N+2:end,:,:);

    alpha = A(N+2,1);

    a = 2*epsilon/(2*epsilon + ht);
    b = ht/(2*epsilon + ht);
    for n = 1:Nt
        Hpn = cell2mat(arrayfun(@(i) H(p(:,i,n)), 1:Nx, 'UniformOutput', false));

        % first stage
        % p1 = p(:, :, n) = pn
        w1 = a*w(:,:,n) + b*Hpn;

        % second stage
        % p2 = p(:, :, n) = pn
        w2 = a*w(:,:,n) + b*w1;

        % third stage
        % 2nd order CFD
        w2x = CFD2(w2, hx);
        p3 = p(:,:,n) - ht*w2x;

        Hp3 = cell2mat(arrayfun(@(i) H(p3(:,i)), 1:Nx, 'UniformOutput', false));
        pnx = CFD2(p(:,:,n), hx);
        w3 = a*(w(:,:,n) - ht*alpha*pnx) + b*(Hpn + Hp3 - w2);

        % update
        w3x = CFD2(w3, hx);
        p3x = CFD2(p3,hx);
        p(:,:,n+1) = p(:,:,n) - ht/2*(w2x + w3x);
        w(:,:,n+1) = w(:,:,n) - ht/2*alpha*(pnx + p3x) + ht/(2*epsilon)*(Hpn + Hp3 - (w2 + w3));
    end

    q = [p; w];
end
