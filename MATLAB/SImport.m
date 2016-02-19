S = zeros(8,8,8);
for i=1:8
S(:,:,i) = eval(strcat('Expression',num2str(i)));
end