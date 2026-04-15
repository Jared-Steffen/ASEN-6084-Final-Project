function [Imb,Bhat,sigma,sigmaB] = BackgroundSubtraction(ImgSeq,sigma_psf)
%{
This function implements Levesque's background estimation algoroithm for
each frame and subtracts it from each frame

%}

% Genreate window sizes, make widths odd
WMS = ceil(5*6*sigma_psf)+1;
WML = ceil(10*6*sigma_psf)+1;    
Wpsf = ceil(6*sigma_psf)+1;


% Initialize
h = size(ImgSeq,1);
w = size(ImgSeq,2);
N = size(ImgSeq,3);
Imb = zeros(h,w,N);
Bhat = zeros(h,w,N);
sigma = zeros(1,N);
sigmaB = zeros(1,N);
for i = 1:N

    % Current image
    Im_i = ImgSeq(:,:,i);
    k = 0;

    % Levesque background estimation
    while k < 5

        Bhatk_small = medfilt2(Im_i, [WMS WMS]);
        Bhatk_large = medfilt2(Im_i, [WML WML]);
        Bhatk = min(Bhatk_small,Bhatk_large);
        
        Ibk = Im_i - Bhatk;

        sigma_k = 1.4826 * mad(Ibk(:),1);        
        vark = sigma_k^2;

        Mbar = imboxfilt(Ibk, [Wpsf Wpsf]);
        E = abs(Mbar) > sigma_k;  % exclude bright and dark structure
        Nmask = ~E;
        
        M = double(Nmask);
        num = conv2(M .* Ibk, ones(Wpsf), 'same');
        den = conv2(M, ones(Wpsf), 'same');
        mu_r = num ./ max(den,1);
        
        varBk = var(mu_r(:),0);
        sigmaBk = sqrt(varBk);
    
        k = k + 1;
    end

    % Save results
    Imb(:,:,i) = Im_i - Bhatk;
    Bhat(:,:,i) = Bhatk;
    sigma(:,i) = sigma_k;
    sigmaB(:,i) = sigmaBk;
end

end

