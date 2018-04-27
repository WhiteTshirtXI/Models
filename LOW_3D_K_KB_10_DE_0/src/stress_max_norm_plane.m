function [n] = stress_max_norm_plane(Shat)
for i = 1:6
S(:,:,:,i) = real(ifftn(Shat(:,:,:,i)));
end
trs = S(:,:,:,1) + S(:,:,:,4) + S(:,:,:,6);
n = max(max(trs(:,:,65)));
end



