function K = elemK(elemID,centroid,diameter,node,elem,mat)
% calculate elem stiffness matrix
k=1; % order

index = elem{elemID};
Nv = length(index);
hK = diameter(elemID);

x = node(index ,1); y = node(index ,2);

Ndof = Nv; Nm = 3;

%% calculate D
D = myM_matrix([x,y],k,hK,centroid(elemID,:)); % Eq.(31)
v1 = 1:Nv; v2 = [2:Nv,1]; % loop index for vertices or edges
Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % normal*length

%% calculate B
E = mat.E; nu =mat.nu;
% plane strain
% C = E*(1-nu)/((1+nu)*(1-2*nu))*[1,nu/(1-nu),0;nu/(1-nu),1,0;0,0,(1-2*nu)/(2*(1-nu))]; % Eq.(18)
% plane stress
C = E/(1-nu^2)*[1,nu,0;nu,1,0;0,0,(1-nu)/2]; % Eq.(17)

elem1 = [v1(:), v2(:)];
e11 = 1/hK*[0,0,0,0,1,0]';
e12 = 1/hK*[0,0,0,1,0,0]';
e21 = 1/hK*[0,0,0,1,0,0]';
e22 = 1/hK*[0,0,0,0,0,1]';

B = zeros(2*Nm,2*Ndof);
for i = 1:2*Nm
    for n = 1:Nv
        phi = zeros(Nv,1);
        phi(elem1(n,1)) = 1;
        phi(elem1(n,2)) = 1;
        phi = blkdiag(phi,phi);
        NE = [Ne(n,1),0;0,Ne(n,2);Ne(n,2),Ne(n,1)];% Eq.(33), normal*length
        B(i,:) = B(i,:)+1/2*[e11(i),e22(i),e12(i)+e21(i)]*C*NE*phi';% Eq.(32)
    end
end

%% constraint, Eq.(27)
B01 = zeros(1,2*Nv);
for n = 1:Nv
    le = [-Ne(n,2),Ne(n,1)];
    Phi = zeros(2,Nv*2);
    Phi(1,elem1(n,:)) = 1; Phi(2,elem1(n,:)+Nv) = 1;
    B01 = B01+1/2*le*Phi; % Eqs.(35) or (37)
end
B02 = zeros(2,2*Nv);
for n = 1:Nv
    Phi = zeros(2,Nv*2);
    Phi(1,elem1(n,:)) = 1; Phi(2,elem1(n,:)+Nv) = 1;
    B02 = B02+1/2*Phi; % Eq.(38)
end
B0 = [B01;B02];
% Bs
Bs = B; Bs([1,2,3],:) = B0;
 
% consistency relation
G = B*D;     Gs = Bs*D;

%% --------- local stiffness matrix ---------
Pis = Gs\Bs;       % Eq.(40)
Pi  = D*Pis;       % Eq.(40)
I = eye(size(Pi));


Ke_c = Pis'*G*Pis; % Eq.(41)

% stabilization term
% To ensure accuracy, you have other options
alpha = 1/Nv*trace(Ke_c); % Eq.(42)
Ke_s = alpha * ((I - Pi)' * (I - Pi)); % Eq.(42)

K  = Ke_c + Ke_s; % Eq.(42)


