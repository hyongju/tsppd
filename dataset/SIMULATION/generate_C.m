function [ c ] = generate_C( vert,v )
% this function helps to generate the cost matrix based on the given map
c = zeros(v,v);
for i = 1:v
    for j = 1:v
        c(i,j) = norm(vert(i,:)-vert(j,:));
    end
end

% make c PSD
for i = 1:v
    sum_c(i) = 0;
    for j = 1:v
        sum_c(i) = sum_c(i) + abs(c(i,j));
    end
    c(i,i) = c(i,i) + sum_c(i);
    cM(i) = c(i,i);
end
for i = 1:c
    c(i,i) = max(cM)*0.5;
end

end

