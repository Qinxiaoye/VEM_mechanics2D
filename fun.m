function N = fun(xi,eta,zeta,ndim,mnode)

N = zeros(mnode,1);

if ndim == 1
    if mnode==2
        N(1) = (1-xi)/2;
        N(2) = (1+xi)/2;
    elseif mnode==3
        N(1) = xi*(xi-1)/2;
        N(2) = -(xi-1)*(xi+1);
        N(3) = xi*(xi+1)/2;
    end
elseif ndim == 2
    if mnode == 4
        N(1) = (1-xi)*(1-eta)/4;
        N(2) = (1+xi)*(1-eta)/4;
        N(3) = (1+xi)*(1+eta)/4;
        N(4) = (1-xi)*(1+eta)/4;
    end
end