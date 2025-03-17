function ff = getForce(node,elem,press,direction)
% trans press to load ff
% you can define different direction of the press, see direction

k = 1; % first order

sumNode = size(node,1);

ff = zeros(sumNode*2,1);
sumP = size(press,1);

nip = k+2;
[x,w] = gaussInt(nip);

for n = 1:sumP
    p = zeros(k+1,1);
    elemID = press(n,1);
    faceID = press(n,2);
    value = press(n,3);
    index = elem{elemID};
    Nv = length(index);
    v1 = 1:Nv; v2 = [2:Nv,1]; % loop index for vertices or edges
    
    elem1 = [v1(:), v2(:)];
    faceNode = elem1(faceID,:);
    faceNodeID = index(faceNode);
    edgeNodeCoor = node(faceNodeID,:);
    
    L = node(faceNodeID(1:2),:);
    L = (L(1,:)-L(2,:));
    Normal = [L(2),-L(1)]/norm(L);
    for m = 1:nip
        N = fun(x(m),0,0,1,k+1);
        dN = dfunc(x(m),0,0,1,k+1);
        J = sqrt((dN'*edgeNodeCoor(:,1))^2+(dN'*edgeNodeCoor(:,2))^2);
        p = p+w(m)*N*value*J; % see Eq.(12)
    end
    if strcmp(direction,'normal')
        ff(faceNodeID) = ff(faceNodeID)+p*Normal(1);
        ff(faceNodeID+sumNode) = ff(faceNodeID+sumNode)+p*Normal(2);
    elseif strcmp(direction,'x')
        ff(faceNodeID) = ff(faceNodeID)+p;
    elseif strcmp(direction,'y')
        ff(faceNodeID+sumNode) = ff(faceNodeID+sumNode)+p;
    end
end
ff = sparse(ff);