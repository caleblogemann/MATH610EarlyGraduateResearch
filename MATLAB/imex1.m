function [q] = imex1(q, A, H, N, Nt, Nx, ht, hx, epsilon)
    [V, D] = eig(A);
    % normalize V
    %for i = 1:2*N+2
    %    V(i,:) = V(i,:)/norm(V(i,:));
    %end
    Dplus = max(zeros(size(D)), D);
    Dminus = min(zeros(size(D)), D);
    Aplus = V*Dplus*inv(V);
    Aminus = V*Dminus*inv(V);
    for n = 1:Nt
        Hpn = cell2mat(arrayfun(@(i) H(q(1:N+1,i,n)), 1:Nx, 'UniformOutput', false));

        q1 = [q(1:N+1,:,n); (q(N+2:end,:,n) + ht/epsilon * Hpn)/(1 + ht/epsilon)];
        q1x = [q1(:,2:end), q1(:,end)] - q1(:,1:end);
        q1xPlus = q1(:,1:end) - [q1(:,1), q1(:,1:end-1)];
        q1xMinus = [q1(:,2:end), q1(:,end)] - q1(:,1:end);
        q(:,:,n+1) = q(:, :, n) + (ht/hx) * (Aplus * q1xPlus + Aminus * q1xMinus) + ht/epsilon * [zeros(N+1, Nx); Hpn - q1(N+2:end,:)];
        %q(:,:,n+1) = q(:, :, n) + (ht/hx) * A * q1x + ht/epsilon * [zeros(N+1, Nx); Hpn - q1(N+2:end,:)];
    end
end
