function [n] = stress_2norm(Shat)
for i = 1:6
S(:,:,:,i) = real(ifftn(Shat(:,:,:,i)));
end
trs = S(:,:,:,1) + S(:,:,:,4) + S(:,:,:,6);
n = sqrt(sum(abs(trs(:))*(2/128)^3));
end

