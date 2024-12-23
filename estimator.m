classdef estimator
    %ESTIMATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        rec_sig
        trans_data
        trans_sig
        sig_size
        real_value
        a = @(theta,K) exp(1j*2*pi*(0:K-1)'*theta)
        source_num
    end
    
    methods
        function obj = estimator(r,x,y,v,num_s)
            %ESTIMATOR Construct an instance of this class
            %   Detailed explanation goes here
            obj.rec_sig = r;
            obj.trans_data = x;
            obj.trans_sig = y;
            obj.source_num = num_s;
            obj.sig_size = size(r);
            obj.real_value = sort(v);
        end
        
        function [time,est,error] = DFT(obj,scale)
            tic
            doa_spectrum = sum(abs(fft(obj.rec_sig,scale)),2);
            [~, peak_idx] = findpeaks(doa_spectrum);
            [~, sorted_idx] = sort(doa_spectrum(peak_idx), 'descend');
            top_idx = peak_idx(sorted_idx);
            est = real(asin(2*sort((top_idx(1:obj.source_num)-1)/(scale)))/pi)*180;
            time = toc;
            error = sort(est) - obj.real_value';
        end
        
        function [time,est,error] = MUSIC(obj,scale)
            tic
            R = obj.rec_sig*obj.rec_sig';
            [L,~] = eig(R);
            doa_spectrum = 1./diag(abs(ifft(fft(L(:,1:obj.sig_size(1)-obj.source_num)*L(:,1:obj.sig_size(1)-obj.source_num)',scale),scale,2)));
            [~, peak_idx] = findpeaks(doa_spectrum);
            [~, sorted_idx] = sort(doa_spectrum(peak_idx), 'descend');
            top_idx = peak_idx(sorted_idx);
            est = real(asin(2*sort((top_idx(1:obj.source_num)-1)/(scale)))/pi)*180;
            time = toc;
            error = sort(est) - obj.real_value';
        end
        
        function [time,est,error] = ES(obj)
            tic
            R = obj.rec_sig*obj.rec_sig';
            [L,~] = eig(R);
            U_f = L(1:obj.sig_size(1)-1,obj.sig_size(1)-obj.source_num+1:end);
            U_b = L(2:end,obj.sig_size(1)-obj.source_num+1:end);
            Psi = (U_b'*U_b)\U_b'*U_f;
            [~,lambda] = eig(Psi);
            est = sort(asin(angle(diag(lambda')) / pi)/pi)*180;
            time = toc;
            error = sort(est) - obj.real_value';
        end
        
        function [time,est,error] = RMUSIC(obj)
            tic
            R = obj.rec_sig*obj.rec_sig';
            [L,~] = eig(R);
            En = L(:,1:obj.sig_size(1)-obj.source_num)*L(:,1:obj.sig_size(1)-obj.source_num)';
            % Coefficient of the polynomial
            M = obj.sig_size(1);
            P = zeros(2 * M - 1, 1);
            for i = 1:2 * M - 1
                    P(i) = sum(diag(En,i-M));
            end
            % ROOT
            roots_P = roots(P);
            % ROOT on UNIT CIRCLE
            [~, sorted_idx] = sort(abs(abs(roots_P) - 1));
            unit_roots = roots_P(sorted_idx((1:obj.source_num)*2-1));
            est = sort(asin(angle(unit_roots') / pi)/pi)*180;
            time = toc;
            error = est - obj.real_value';
        end
        
        function [time,est,error] = CAPON(obj,scale)
            tic
            R = obj.rec_sig*obj.rec_sig';
            doa_spectrum = 1./diag(abs(ifft(fft(inv(R),scale),scale,2)));
            [~, peak_idx] = findpeaks(doa_spectrum);
            [~, sorted_idx] = sort(doa_spectrum(peak_idx), 'descend');
            top_idx = peak_idx(sorted_idx);
            est = real(asin(2*sort((top_idx(1:obj.source_num)-1)/(scale)))/pi)*180;
            time = toc;
            error = sort(est) - obj.real_value';
        end
        
        function [time,est,error] = OMP(obj,scale,varargin)
            tic
            est = zeros(obj.source_num,1);
            y = obj.rec_sig;
            for i = 1:obj.source_num
                % estimate the point with the largest energy
                doa_spectrum = sum(abs(fft(y,scale)),2);
                [~,pos] = max(doa_spectrum);
                % reconstruct signal
                sig_recst = obj.a((pos-1)/scale,obj.sig_size(1));
                y = y-sig_recst*((sig_recst'*y)/(sig_recst'*sig_recst));
                est(i) = real(asin(2*(pos-1)/scale))/pi*180;
            end
            time = toc;
            est = sort(est);
            if numel(varargin)~=0
                disp(cell2mat(varargin(1)))
            end
            error = sort(est) - obj.real_value';
        end
        
        function [time,est,error] = ML(obj,scale)
            tic
            DICT_ML = obj.a((-scale/2:scale/2-1)/scale,obj.sig_size(1));
            A_inv = DICT_ML'*DICT_ML;
            eval(['P_ML = zeros(',char(kron(ones(1,obj.source_num-1),'scale,')),'scale);'])
            for i = 1:scale
                for j = 1:scale
                    P_ML(i,j) = sum(sum(obj.rec_sig'*DICT_ML(:,[i,j])/(A_inv([i,j],[i,j]))*DICT_ML(:,[i,j])'*obj.rec_sig));
                end
            end
            [~,pos] = max(reshape(abs(P_ML),[scale^2,1]));
            est(1) = ((ceil(pos/scale)-1)/scale-1/2);
            est(2) = ((mod(pos,scale)-1)/scale-1/2);
            time = toc;
            est = sort(real(asin(2*est))/pi*180);
            error = est - obj.real_value';
        end
        
        function [time,est,error] = ANM(obj)
            y = obj.rec_sig;
            N = obj.sig_size(1);
            L = obj.sig_size(2);
            tic
            cvx_begin sdp quiet 
            cvx_solver 
            variable T(N, N) complex hermitian toeplitz
            variable x(L, L) complex hermitian
            variable z(N, L) complex
            minimize 0.9*(trace(x) + trace(T))+0.1*sum_square_abs(vec(y-z))
            [ x z'; z T] >= 0;
            cvx_end
            [Phi,~] = rootmusic(T, obj.source_num , 'corr'); %%等效于对T矩阵进行范德蒙德分解
            est = real(asin(Phi/ pi))/pi*180;
            time = toc;
            error = sort(est) - obj.real_value';
        end
    end
end

