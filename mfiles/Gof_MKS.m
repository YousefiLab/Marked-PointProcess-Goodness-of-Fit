function MKS = Gof_MKS(X,C_a)

% size of samples
n = size(X,1);
d = size(X,2);
% create perms
comb = perms(1:d);
% generate KS
KS = zeros(size(comb,1),n);
FS = zeros(1,n);
% loop
for i=1:n
    bound = X(i,:);
    for j=1:size(comb,1)
        temp = bound(comb(j,:));
        mat = zeros(n,d);
        for k=1:d
            mat(:,k) = X(:,k)<=temp(k);
        end
        I = ismember(mat,ones(1,d),'rows');
        KS(j,i) = length(find(I));
        KS(j,i) = KS(j,i)/n;
    end
    FS(i) = prod(bound);  
end
[FS,~]=sort(FS);
for i=1:size(comb,1)
    D_KS(i) = max(abs(FS-KS(i,:)));
end
D_ks = max(D_KS);
MKS = D_ks<=C_a;






