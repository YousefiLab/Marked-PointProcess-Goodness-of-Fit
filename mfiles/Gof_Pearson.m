%% Uniformity tests
% Pearson X^2 -Statistic Test
function chi_squ = Gof_Pearson(X,DF,P)

n = size(X,1); %% The number of samples
d = size(X,2);
ind1   = 1;
sub_r = 0.2;
for i = sub_r:sub_r:1
    mat = zeros(n,d);
    for j=1:d
        mat(:,j) = (X(:,d)>=i-sub_r & X(:,d)<i);
    end
    I = ismember(mat,ones(1,d),'rows');
    numT = length(find(I));
    n_r_T(ind1) = numT;
    ind1   = ind1+1;
end
P_E_T   = (n/length(n_r_T))*ones(1,length(n_r_T));
chi_s   = sum( ((n_r_T-P_E_T).^2)./P_E_T );
C_a = chi2inv(1-P,DF);
chi_squ = chi_s<=C_a;
