

clear all;close all;clc
n = 3;
p_all = perms(1:n);

for i = 1:size(p_all,1)
    mat{i} = zeros(n,n);
    for j = 1:size(p_all,2)
        mat{i}(j,p_all(i,j)) = 1;
    end
end
for i = 1:length(mat)
    result(i) = all(eig(mat{i})>=0);
end
p_all((find(result)),:)
mat{find(result)}



