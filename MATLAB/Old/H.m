function [Hp] = H_Example1(p,S)
    N = length(p);
    Hp = zeros(N,1);
    for k = 1:N
        for i = 1:N
            for j = 1:N
                Hp(k) = Hp(k) + p(i)*p(j)*S(i,j,k);
            end
        end
    end
    Hp = 1/2*Hp;
end

