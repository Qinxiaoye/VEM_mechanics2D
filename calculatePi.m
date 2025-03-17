function [dM,Pis,Pi,area] = calculatePi(node)

Nv = size(node,1);

hK = max(pdist(node));   

verts = node; verts1 = circshift(verts,-1); % verts([2:end,1],:);
area_components = verts(:,1).*verts1(:,2)-verts1(:,1).*verts(:,2);
area = 0.5*abs(sum(area_components));
centroid = sum((verts+verts1).*repmat(area_components,1,2))/(6*area);



x = node(: ,1); y = node(: ,2);


D = myM([x,y],1,hK,centroid,2);
v1 = 1:Nv; v2 = [2:Nv,1]; % loop index for vertices or edges
Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % 边长*外法线
% calculate B
Gradm = [0 0; 1./hK*[1, 0]; 1./hK*[0, 1]]; % k = 1
rotid1 = [Nv,1:Nv-1]; rotid2 = [2:Nv,1]; % ending and starting indices
normVec = 0.5*[y(rotid2)-y(rotid1), x(rotid1)-x(rotid2)]'; % a rotation of edge vector
B = Gradm*normVec; % B

% constraint
Bs = B;  Bs(1,:) = 1/Nv;
% consistency relation
G = B*D;  Gs = Bs*D;

Pis = Gs\Bs;
Pi = D*Pis;


dM = myGradmc(centroid,1,hK,centroid,2);


