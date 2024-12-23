N = 64;  % number of antennas
L = 1;
f = [0.13, 0.41, 0.32];  % f = cos(theta)
A = exp(1j * 2 * pi * (0 : N-1)' * f);  % A = [a(f1), ..., a(fK)]
s = ones(3, 1);  % 3 sources
% z = A * s;
snr = 10;
y = awgn(A * s,snr,'measured');
lambda = sqrt(N * (L + log(N) + sqrt( 2 * L * log(N)))*10^(-snr/10));

cvx_begin sdp quiet 
cvx_solver
variable T(N, N) hermitian toeplitz
variable x 
variable z(N) complex
minimize 0.9*(0.5 * x + 0.5 * T(1,1))+0.5*sum_square_abs(y-z)
[ x z'; z T] >= 0;
cvx_end

[Phi,P] = rootmusic(T, 3, 'corr'); %%等效于对T矩阵进行范德蒙德分解
Phi / 2 / pi