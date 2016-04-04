function [ Hvalue ] = HExample4(p, sigma, z, w, phiz)
H = @(p, i) (1 + sigma*z(i))*abs(p);

integrand = @(p, i) H(p'*phiz(:,i), i)*phiz(:,i);

Hvalue = cell2mat(arrayfun(@(i)integrand(p, i), 1:length(z), 'UniformOutput', false))*w;

end

