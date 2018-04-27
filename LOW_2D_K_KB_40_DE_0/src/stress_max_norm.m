function [n] = stress_max_norm(Shat)
S = real(ifft2(Shat));
trs = S(:,:,1) + S(:,:,3);
n = max(trs(:));
end

