function [ F,phi ] = GEM( data,label,k,theta,gamma)
%GEM Summary of this function goes here
%k is the number of clusters we want to get
[dimension,number] = size(data);
%F = zeros(dimension,k);
%F = zeros(0,784);
C = zeros(dimension,dimension,k);
c = zeros(dimension,dimension);
sumup = zeros(dimension);
sumdown = zeros(dimension);
for m = 1:k
    for i = 1:number
        if(label(i) == m)
            sumup = sumup + data(:,i)*data(:,i)';
            sumdown = sumdown+ones(dimension);
        end
    end
    C(:,:,m) = sumup./sumdown;
    sumup = zeros(dimension);
    sumdown = zeros(dimension);
end
f = zeros(0,dimension);
F = zeros(0,dimension);
for i = 1:k
    for j = 1:k
        if j ~= i
            [V,D] = eig(C(:,:,i),C(:,:,j)+(gamma./dimension).*trace(C(:,:,j)).*eye(dimension));
            [rows,cols] = find(D >= theta);
            f = union(F,V(:,rows)','rows');
            F = f;
        end
        clear f
    end
end
delta = 1;
alpha =1;

phi = abs(delta.*F*data);
phi = phi.^(alpha/2);
%omiga = mnrfit(phi',label);


end
