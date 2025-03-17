function dN = dfunc(xi,eta,zeta,ndim,mnode)

dN = zeros(mnode,ndim);

if ndim == 1
    if mnode==2
        dN(1) = -1/2;
        dN(2) = 1/2;
    elseif mnode==3
        dN(1) = 1/2*(2*xi-1);
        dN(2) = -2*xi;
        dN(3) = 1/2*(2*xi+1);
    end
end