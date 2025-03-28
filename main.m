% Copyright:
% Bingbing Xu
% Contact:
% xubingbingdut@foxmail.com, bingbing.xu@ikm.uni-hannover.de
% Institute of Continuum Mechanics, Leibniz Universität Hannover

% reference: https://www.sciencedirect.com/science/article/pii/S0045782519301215
%            https://github.com/Terenceyuyue/mVEM

clear;
load("cook.mat"); % you can generate mesh by your self

k = 1; % first order VEM for elasticity

E = 200000; nu = 0.3; % material parameters
mat.E = E; mat.nu = nu;

sumNode = size(node,1);

p = 100; % press

% boundary condition
nodeL = find(node(:,1)<0.001);
sL = size(nodeL,1);
nodeR  = find(node(:,1)>48-0.001);
% find face by node
pface = findFace(node,elem,nodeR);
press = [pface,ones(length(pface),1)*p];
fixNode = [nodeL,ones(sL,1),zeros(sL,1);nodeL,2*ones(sL,1),zeros(sL,1)];
fixMes = [fixNode(:,1)+(fixNode(:,2)-1)*sumNode,fixNode(:,3)];

nodeForce = getForce(node,elem,press,'y');

% global stiffness matrix
GK = globalK(node,elem,mat);

[GK,F] = boundaryCondition(GK,nodeForce,fixMes(:,1),fixMes(:,2),1);

uh = GK\F;

uh = full(uh);

ux = uh(1:sumNode);
uy = uh(sumNode+1:end);


% stress
[stress,mises] = calculateStress(node,elem,uh,mat);

figure;
showsolution(node+[ux,uy],elem,mises);
figure;
showsolution(node+[ux,uy],elem,uy);



