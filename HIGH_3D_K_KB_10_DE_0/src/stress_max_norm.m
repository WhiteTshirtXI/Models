function [n] = stress_max_norm(Shat)
for i = 1:6
S(:,:,:,i) = real(ifftn(Shat(:,:,:,i)));
end
trs = S(:,:,:,1) + S(:,:,:,4) + S(:,:,:,6);
n = max(trs(:));
end

