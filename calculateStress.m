function [stress,mises] = calculateStress(node,elem,uh,mat)

sumElem = size(elem,1); % the number of element
sumNode = size(node,1);

ux = uh(1:sumNode);
uy = uh(sumNode+1:end);

stress = zeros(sumNode,3);

nodeUsed = zeros(sumNode,1);

for n = 1:sumElem
    index = elem{n};

    coor = node(index,:);
    [dM,Pis] = calculatePi(coor);

    A1 = [1,0;0,0;0,1];
    A2 = [0,0;0,1;1,0];
    A = [A1,A2];

    Pis = blkdiag(Pis,Pis);

    E = mat.E; nu =mat.nu;
    % D = E*(1-nu)/((1+nu)*(1-2*nu))*[1,nu/(1-nu),0;nu/(1-nu),1,0;0,0,(1-2*nu)/(2*(1-nu))]; % plane strain
    C = E/(1-nu^2)*[1,nu,0;nu,1,0;0,0,(1-nu)/2];

    dM1 = blkdiag(dM',dM');
    eS = C*A*dM1*Pis*[ux(index);uy(index)];


    stress(index,:) = stress(index,:)+eS';
    nodeUsed(index) = nodeUsed(index)+1;
end


stress = stress./nodeUsed;

smax = (stress(:,1)+stress(:,2))./2+sqrt(((stress(:,1)-stress(:,2))./2).^2+stress(:,3).^2);
smin = (stress(:,1)+stress(:,2))./2-sqrt(((stress(:,1)-stress(:,2))./2).^2+stress(:,3).^2);


mises = sqrt((smax.^2+smin.^2+(smax-smin).^2)./2);
