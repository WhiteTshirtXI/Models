function [Min,Max] = SE_MinMax_Eig(Shat,lam,T)
for k = 1:6
    M(:,:,:,k) = real(ifftn(Shat(:,:,:,k)));
end
for a = 1:128
    for b = 1:128
        for c = 1:128
            C = (lam/T)*[M(a,b,c,1), M(a,b,c,2),M(a,b,c,3); M(a,b,c,2), M(a,b,c,4),M(a,b,c,5);...
                         M(a,b,c,3), M(a,b,c,5),M(a,b,c,6)] + eye(3,3);
            E(a,b,c,:) = eig(C);
        end
    end
end 
Min = min(E(:));
ax = max(E(:));
end

