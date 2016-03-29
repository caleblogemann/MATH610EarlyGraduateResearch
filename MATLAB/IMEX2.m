function [q] = imex2(q, A, H, N, Nt, Nx, ht, hx, epsilon)
    p = q(1:N+1,:,:);
    w = q(N+2:end,:,:);

    alpha = A(N+2,1);

    for n = 1:Nt
        Hpn = cell2mat(arrayfun(@(i) H(p(:,i,n)), 1:Nx, 'UniformOutput', false));

        w1 = (w(:,:,n) + ht/(2*epsilon) * Hpn)/(1 + ht/(2*epsilon));
        w2 = (w(:,:,n) + ht/(2*epsilon) * w1)/(1 + ht/(2*epsilon));

        % 2nd order CFD
        w2x = CFD2(w2, hx);

        p3 = p(:,:,n) + ht* w2x;

        Hp3 = cell2mat(arrayfun(@(i) H(p3(:,i)), 1:Nx, 'UniformOutput', false));

        p3x = CFD2(p3,hx);

        w3 = (w(:,:,n) + ht*alpha*p3x + ht/(2*epsilon)*(Hpn + Hp3 - w2))/(1 + ht/(2*epsilon));

        w3x = CFD2(w3, hx);
        p2x = CFD2(p(:,:,n), hx);

        p(:,:,n+1) = p(:,:,n) + ht/2 *(w2x + w3x);
        w(:,:,n+1) = w(:,:,n) + ht/2 *alpha*(p2x + p3x) + ht/(2*epsilon) * (Hpn + Hp3 - (w2 + w3));
    end

    q = [p; w];
end
