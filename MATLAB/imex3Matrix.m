function [q] = imex3Matrix(q, A, H, N, Nt, Nx, ht, hx, epsilon)
    [V, D] = eig(A);
    Dplus = max(zeros(size(D)), D);
    Dminus = min(zeros(size(D)), D);
    Aplus = V*Dplus*inv(V);
    Aminus = V*Dminus*inv(V);

    a = 0.24169426078821;
    b = 0.06042356519705;
    eta = 0.1291528696059;

    c = epsilon/(epsilon - a*ht);
    d = a*ht/(epsilon - a*ht);

    B = [eye(N+1), zeros(N+1); zeros(N+1), c*eye(N+1)];
    for n = 1:Nt
        Hpn = cell2mat(arrayfun(@(i) H(q(1:N+1,i,n)), 1:Nx, 'UniformOutput', false));

        % first stage
        q1 = B*q(:,:,n) + [zeros(N+1, Nx); d*Hpn];

        % second stage
        q2 = B*q(:,:,n) + [zeros(N+1, Nx); d*q1(N+2:end, :)];

        % third stage
        % solve so p3 is known
        q2xPlus = BFD3([q2(:,1), q2(:,1), q2(:,1), q2], hx);
        q2xMinus = FFD3([q2, q2(:,end), q2(:,end), q2(:,end)], hx);
        q3 = B*(q(:,:,n) - ht*(Aplus*q2xPlus + Aminus*q2xMinus));
        % add Hp3 in to w3
        Hp3 = cell2mat(arrayfun(@(i) H(q3(1:N+1,i)), 1:Nx, 'UniformOutput', false));
        q3(N+2:end,:) = q3(N+2:end,:) + d*((1/a - 1)*(Hpn - q2(N+2:end, :)) + Hp3);

        % fourth stage
        q3xPlus = BFD3([q3(:,1), q3(:,1), q3(:,1), q3], hx);
        q3xMinus = FFD3([q3, q3(:,end), q3(:,end), q3(:,end)], hx);
        q4 = B*(q(:,:,n) - ht/4*(Aplus*(q2xPlus + q3xPlus) + Aminus*(q2xMinus + q2xPlus)));
        % add Hp4 in to w4
        Hp4 = cell2mat(arrayfun(@(i) H(q4(1:N+1,i)), 1:Nx, 'UniformOutput', false));
        q4(N+2:end,:) = q4(N+2:end,:) + d*(b/a*(Hpn - q1(N+2:end,:)) + ...
            eta/a*(Hpn - q2(N+2:end,:)) + (1/2 - b - eta - a)/a*(Hp3 - q3(N+2:end,:)) +...
            Hp4);

        % update
        q4xPlus = BFD3([q4(:,1), q4(:,1), q4(:,1), q4], hx);
        q4xMinus = FFD3([q4, q4(:,end), q4(:,end), q4(:,end)], hx);
        q(:,:,n+1) = q(:,:,n) - ht/6*(Aplus*(q2xPlus + q3xPlus + 4*q4xPlus) ...
            + Aminus*(q2xMinus + q3xMinus + 4*q4xMinus)) ...
            + ht/(6*epsilon)*[zeros(N+1, Nx); Hpn + Hp3 - (q2(N+2:end,:) ...
            + q3(N+2:end,:)) + 4*(Hp4 - q4(N+2:end,:))];
    end
end
