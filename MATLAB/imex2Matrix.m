function [q] = imex2Matrix(q, A, H, N, Nt, Nx, ht, hx, epsilon)
    [V, D] = eig(A);
    Dplus = max(zeros(size(D)), D);
    Dminus = min(zeros(size(D)), D);
    Aplus = V*Dplus*inv(V);
    Aminus = V*Dminus*inv(V);

    a = 2*epsilon/(2*epsilon + ht);
    b = ht/(2*epsilon + ht);

    B = [eye(N+1), zeros(N+1); zeros(N+1), a*eye(N+1)];
    for n = 1:Nt
        Hpn = cell2mat(arrayfun(@(i) H(q(1:N+1,i,n)), 1:Nx, 'UniformOutput', false));

        % first stage
        q1 = B*q(:,:,n) + [zeros(N+1, Nx); b*Hpn];

        % second stage
        q2 = B*q(:,:,n) + [zeros(N+1, Nx); b*q1(N+2:end, :)];

        % third stage
        % solve so p3 is known
        q2xPlus = BFD2([q2(:,1), q2(:,1), q2], hx);
        q2xMinus = FFD2([q2, q2(:,end), q2(:,end)], hx);
        q3 = B*(q(:,:,n) - ht*(Aplus*q2xPlus + Aminus*q2xMinus));
        % add Hp3 in to w3
        Hp3 = cell2mat(arrayfun(@(i) H(q3(1:N+1,i)), 1:Nx, 'UniformOutput', false));
        q3 = q3 + [zeros(N+1, Nx); b*(Hpn - q2(N+2:end, :) + Hp3)];

        % update
        q3xPlus = BFD2([q3(:,1), q3(:,1), q3], hx);
        q3xMinus = FFD2([q3, q3(:,end), q3(:,end)], hx);
        q(:,:,n+1) = q(:,:,n) - ht/2*(Aplus*(q2xPlus + q3xPlus) + Aminus*(q2xMinus + q3xMinus)) + ht/(2*epsilon)*[zeros(N+1, Nx); Hpn + Hp3 - (q2(N+2:end, :) + q3(N+2:end, :))];
    end
end
