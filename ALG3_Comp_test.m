clear;close;
a = @(theta,K) exp(1j*2*pi*(0:K-1)'*theta);
L = 100;  % symbols
theta = [15,20];
f = sin(theta/180*pi)/2;  % f = cos(theta)
K = 32;
scale = 256;
SNR = 1e-1;
num_s = 2;
i=1;
A = a(f,K);
OMP_scaler = 4;% Control the scaler coefficient of OMP algorithm

Data = randi([0,3],num_s,L);
s = 1/sqrt(2)*qammod(Data,4,'gray');  % 2 sources
n = sqrt(SNR/2)*(randn(K,L)+1j*randn(K,L));
y = A*s+n;
Estm = estimator(y,s,A*s,theta,num_s);
% [t_dft,est_dft,error_dft] = Estm.DFT(scale);
% [t_music,est_music,error_music] = Estm.MUSIC(scale);
% [t_rmusic,est_rmusic,error_rmusic] = Estm.RMUSIC();
[t_es,est_es,error_es] = Estm.ES();
% [t_capon,est_capon,error_capon] = Estm.CAPON(scale);
% [t_anm,est_anm,error_anm] = Estm.ANM();
% [t_ml,est_ml,error_ml] = Estm.ML(scale);
% [t_omp,est_omp,error_omp] = Estm.OMP(OMP_scaler*scale);
