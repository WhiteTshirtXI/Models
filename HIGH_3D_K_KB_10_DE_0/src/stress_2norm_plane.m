function [n] = stress_2norm_plane(Shat)
for i = 1:6
S(:,:,:,i) = real(ifftn(Shat(:,:,:,i)));
end
trs = S(:,:,:,1) + S(:,:,:,4) + S(:,:,:,6);
B = trs(:,:,65);
n = sqrt(sum(abs(B(:)))*(2/128)^2);
end

