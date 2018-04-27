function [n] = stress_2norm(Shat)
S = real(ifft2(Shat));
trs = S(:,:,1) + S(:,:,3);
n = sqrt(sum(abs(trs(:))*(2/128)^3));
end

