clear;close;
a = @(theta,K) exp(1j*2*pi*(0:K-1)'*theta);
L = 100;  % symbols
scale = 1024;
SNR = 1e-2;
K = 32; % number of antennas
Data = randi([0,3],2,L);
s = 1/sqrt(2)*qammod(Data,4,'gray');  % 3 sources
n = sqrt(SNR/2)*(randn(K,L)+1j*randn(K,L));
theta = [-5,10];
f = sin(theta/180*pi)/2;  % f = cos(theta)
A = a(f,K);  % A = [a(f1), ..., a(fK)]
y = A*s+n;
p_dft(:) = sum(abs(fft(y,scale)),2);
p_dft(:) = p_dft(:)./max(p_dft(:));
f = figure(1);
axes1 = axes('Parent',f);
hold(axes1,'on')
view(axes1,[-30 19.2]);
xlabel('入射角/°')
ylabel('入射角度差/°')
zlabel('归一化功率/°')
plot3(axes1,asin(2*(-0.5:1/1024:0.5-1/1024))/pi*180,15*ones(scale,1),...
    [p_dft(scale/2:end),p_dft(1:scale/2-1)],'b-');
[L,~] = eig(y*y');
p_music = 1./diag(abs(ifft(fft(L(:,1:end-2)*L(:,1:end-2)',scale),scale,2)));
p_music = p_music(:)./max(p_music(:));
plot3(axes1,asin(2*(-0.5:1/1024:0.5-1/1024))/pi*180,15*ones(scale,1),...
    [p_music(scale/2:end);p_music(1:scale/2-1)],'r-');
legend('DFT','MUSIC','AutoUpdate','off','Location','northeast')

theta = [0,10];
f = sin(theta/180*pi)/2;  % f = cos(theta)
A = a(f,K);  % A = [a(f1), ..., a(fK)]
y = A*s+n;
p_dft(:) = sum(abs(fft(y,scale)),2);
p_dft(:) = p_dft(:)./max(p_dft(:));
[L,~] = eig(y*y');
p_music = 1./diag(abs(ifft(fft(L(:,1:end-2)*L(:,1:end-2)',scale),scale,2)));
p_music = p_music(:)./max(p_music(:));
plot3(axes1,asin(2*(-0.5:1/1024:0.5-1/1024))/pi*180,10*ones(scale,1),...
    [p_music(scale/2:end);p_music(1:scale/2-1)],'r-');
plot3(axes1,asin(2*(-0.5:1/1024:0.5-1/1024))/pi*180,10*ones(scale,1),...
    [p_dft(scale/2:end),p_dft(1:scale/2-1)],'b-');



theta = [5,10];
f = sin(theta/180*pi)/2;  % f = cos(theta)
A = a(f,K);  % A = [a(f1), ..., a(fK)]
y = A*s+n;
p_dft(:) = sum(abs(fft(y,scale)),2);
p_dft(:) = p_dft(:)./max(p_dft(:));
[L,~] = eig(y*y');
p_music = 1./diag(abs(ifft(fft(L(:,1:end-2)*L(:,1:end-2)',scale),scale,2)));
p_music = p_music(:)./max(p_music(:));
plot3(axes1,asin(2*(-0.5:1/1024:0.5-1/1024))/pi*180,5*ones(scale,1),...
    [p_music(scale/2:end);p_music(1:scale/2-1)],'r-');
plot3(axes1,asin(2*(-0.5:1/1024:0.5-1/1024))/pi*180,5*ones(scale,1),...
    [p_dft(scale/2:end),p_dft(1:scale/2-1)],'b-');


at = 2;
theta = [10-at,10];
f = sin(theta/180*pi)/2;  % f = cos(theta)
A = a(f,K);  % A = [a(f1), ..., a(fK)]
y = A*s+n;
p_dft(:) = sum(abs(fft(y,scale)),2);
p_dft(:) = p_dft(:)./max(p_dft(:));
[L,~] = eig(y*y');
p_music = 1./diag(abs(ifft(fft(L(:,1:end-2)*L(:,1:end-2)',scale),scale,2)));
p_music = p_music(:)./max(p_music(:));
plot3(axes1,asin(2*(-0.5:1/1024:0.5-1/1024))/pi*180,at*ones(scale,1),...
    [p_music(scale/2:end);p_music(1:scale/2-1)],'r-');
plot3(axes1,asin(2*(-0.5:1/1024:0.5-1/1024))/pi*180,at*ones(scale,1),...
    [p_dft(scale/2:end),p_dft(1:scale/2-1)],'b-');