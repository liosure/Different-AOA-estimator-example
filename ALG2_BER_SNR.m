clear;
N = 32;  % number of antennas
L = 50;
f = [0.13, 0.41];  % f = cos(theta)
A = exp(1j * 2 * pi * (0 : N-1)' * f);  % A = [a(f1), ..., a(fK)]
scale = 128;
LOOP_NUM = 100; % ���ؿ���ʵ�����
SNR_range = -15:3:9; % ����ȷ�Χ
true_param1 = 5; % ʾ������1
true_param2 = 10; % ʾ������2

% Ԥ����������
num_SNR = length(SNR_range);
anm_parameter_errors = zeros(num_SNR, 1);
dft_parameter_errors = zeros(num_SNR, 1);
mus_parameter_errors = zeros(num_SNR, 1);
cap_parameter_errors = zeros(num_SNR, 1);
anm_time = zeros(num_SNR, 1);
dft_time = zeros(num_SNR, 1);
mus_time = zeros(num_SNR, 1);
cap_time = zeros(num_SNR, 1);

% ���ؿ���ʵ��
for idx = 1:num_SNR
    SNR_dB = SNR_range(idx); % ��ǰ�����
    SNR_linear = 10^(SNR_dB / 10); % ���������
    errors_anm = zeros(LOOP_NUM, 1); % �洢ÿ��ʵ��Ĳ����������
    errors_dft = zeros(LOOP_NUM, 1);
    errors_music = zeros(LOOP_NUM, 1);
    errors_capon = zeros(LOOP_NUM, 1);
    time_anm = zeros(LOOP_NUM, 1); 
    time_dft = zeros(LOOP_NUM, 1); 
    time_music = zeros(LOOP_NUM, 1); 
    time_capon = zeros(LOOP_NUM, 1); 
    parfor loop = 1:LOOP_NUM
        % ����ģ�����ݺ�����
        Data = randi([0,3],2,L);
        s = 1/sqrt(2)*qammod(Data,4,'gray');  % 3 sources
        y = awgn(A * s,SNR_dB,'measured');
        [estimated_params,t_anm] = anm(N,L,SNR_dB,y); 
        [estimated_music,t_music] = music(y); 
        [estimated_dft,t_dft] = dft(y,N,scale); 
        [estimated_capon,t_capon] = capon(y,N,scale); 
        % ��������������
        errors_anm(loop) = norm(estimated_params - f'); % 
        errors_music(loop) = norm(estimated_music - f');
        errors_dft(loop) = norm(estimated_dft - f'); % 
        errors_capon(loop) = norm(estimated_capon - f');
        time_anm(loop) = t_anm;
        time_dft(loop) = t_dft;
        time_music(loop) = t_music;
        time_capon(loop) = t_capon;
    end
    
    % ����ÿ��SNR�����µ�ƽ���������
    anm_parameter_errors(idx) = asin(sqrt(mean(abs(errors_anm).^2)));
    dft_parameter_errors(idx) = asin(sqrt(mean(abs(errors_dft).^2)));
    mus_parameter_errors(idx) = asin(sqrt(mean(abs(errors_music).^2)));
    cap_parameter_errors(idx) = asin(sqrt(mean(abs(errors_capon).^2)));
    anm_time(idx) = mean(time_anm);
    dft_time(idx) = mean(time_dft);
    mus_time(idx) = mean(time_music+time_capon);
    cap_time(idx) = mean(time_capon);
end
semilogy(SNR_range,anm_parameter_errors,'b--x')
hold on
semilogy(SNR_range,dft_parameter_errors,'r--x')
semilogy(SNR_range,mus_parameter_errors,'m--x')
semilogy(SNR_range,cap_parameter_errors,'k--x')
xlabel("SNR/dB")
ylabel("MMSE/��")
figure()
semilogy(SNR_range,anm_time,'b--x')
hold on
semilogy(SNR_range,dft_time,'r--x')
semilogy(SNR_range,mus_time,'m--x')
semilogy(SNR_range,cap_time,'k--x')
xlabel("SNR/dB")
ylabel("time")
% ������
save('snr-mse1.mat')
% �Զ����������ɺ���

% �Զ���������ƺ���
function [est_anm,t1] = anm(N,L,snr,y)
    tic
%     lambda = sqrt(N * (L + log(N) + sqrt( 2 * L * log(N)))*10^(-snr/10)/2);
    cvx_begin sdp quiet 
    cvx_solver
    variable T(N, N) complex hermitian toeplitz
    variable x(L, L) complex hermitian
    variable z(N, L) complex
    minimize 0.5*(trace(x) + trace(T))+0.5*sum_square_abs(vec(y-z))
    [ x z'; z T] >= 0;
    cvx_end
    [Phi,~] = rootmusic(T, 2, 'corr'); %%��Ч�ڶ�T������з����ɵ·ֽ�
    est_anm = sort([Phi / 2 / pi]);
    t1 = toc;
end
function [est,t1] = music(y)
    tic
    [Phi,~] = rootmusic(y*y', 2, 'corr');
    est = sort([Phi / 2 / pi]);
    t1 = toc;
end
function [est,t1] = dft(y,M,scale)
    tic
    doa_spectrum=sum(abs(fft(y,M*scale)),2); 
    t1 = toc;
    [~, peak_idx] = findpeaks(doa_spectrum);
    [~, sorted_idx] = sort(doa_spectrum(peak_idx), 'descend');
    top_idx = peak_idx(sorted_idx);
    est = sort((top_idx(1:2))/(M*scale));
end
function [est,t1] = capon(y,M,scale)
    tic
    doa_spectrum = abs(1./diag(ifft(fft(inv(y*y'),M*scale),M*scale,2))); 
    t1 = toc;
    [~, peak_idx] = findpeaks(doa_spectrum);
    [~, sorted_idx] = sort(doa_spectrum(peak_idx), 'descend');
    top_idx = peak_idx(sorted_idx);
    est = sort((top_idx(1:2))/(M*scale));
end