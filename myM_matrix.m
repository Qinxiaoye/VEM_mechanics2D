function M = myM_matrix(node,k,hK,c)
% calculate scaled monomials(matrix)
% c: center

x = node(:,1);
y = node(:,2);
xK = c(1);
yK = c(2);

xi = (x-xK)./hK; % Eq.(23)
eta = (y-yK)./hK;% Eq.(23)


if k == 1
    M = [[1.0+0*x;0*x],[0*x;1.0+0*x],[-eta;xi],[eta;xi],[xi;0*x],[0*x;eta]]; % Eq.(22)
else
    disp('We are implementing case k > 1');
end

