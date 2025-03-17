function[KNew,FNew] = boundaryCondition(K,F,fixDof,gievDisp,m)


KNew = K;
FNew = F;

if m == 1

    FNew(fixDof) = gievDisp;

    KNew(fixDof,:) = 0;
    KNew = KNew+sparse(fixDof,fixDof,ones(size(fixDof,1),1),size(K,1),size(K,1));

else % 乘大数法
    alpha=max(K(:))*1e5;
    for n = 1:length(fixDof)
        dof = fixDof(n);
        KNew(dof,dof) = alpha*K(dof,dof);
        FNew(dof) = alpha*K(dof,dof)*gievDisp(n);
    end
end