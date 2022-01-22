N=1024;
Ncp=16;
Fs=1e3;
Ts=1/Fs;
Tsample=Ts;
Fd=0.1;
Np=32;
M=4;
k=log2(M);
Nframes=1e3;
pathDelays = 0;
avgPathGains = 0;


awgnchan = comm.AWGNChannel(...
'NoiseMethod','Signal to noise ratio (SNR)');

rayChan = comm.RayleighChannel('SampleRate',Fs, ...
    'PathDelays',pathDelays, ...
    'AveragePathGains',avgPathGains, ...
    'MaximumDopplerShift',Fd);

ricianChan = comm.RicianChannel('SampleRate',Fs, ...
    'PathDelays',pathDelays, ...
    'AveragePathGains',avgPathGains, ...
    'KFactor',3, ...
    'MaximumDopplerShift',Fd);


mimo2v2Chan = comm.MIMOChannel(...
    'SampleRate',Fs, ...
    'PathDelays',0, ...
    'AveragePathGains',0, ...
    'MaximumDopplerShift',Fd, ...
    'SpatialCorrelationSpecification','None', ...
    'PathGainsOutputPort',true,...    
    'NumTransmitAntennas',2, ...
    'NumReceiveAntennas',2);

mimo2v4Chan = comm.MIMOChannel(...
    'SampleRate',Fs, ...
    'PathDelays',pathDelays, ...
    'AveragePathGains',avgPathGains, ...
    'MaximumDopplerShift',Fd, ...
    'SpatialCorrelationSpecification','None', ...
    'PathGainsOutputPort',true,...
    'NumTransmitAntennas',2, ...
    'NumReceiveAntennas',4);


OSTBCEnc = comm.OSTBCEncoder('NumTransmitAntennas',2);

OSTBCComb2 = comm.OSTBCCombiner(...
    'NumTransmitAntennas',2,...
    'NumReceiveAntennas',2);

OSTBCComb4 = comm.OSTBCCombiner(...
    'NumTransmitAntennas',2,...
    'NumReceiveAntennas',4);

Data= randi([0,M-1],N,Nframes);
Data(1:Np:N,:)=0;
Data_mod=qammod(Data,M);
Data_mod(1:Np:N,:)=1;

IFFT_Data = (N/sqrt(N-N/Np))*ifft(Data_mod,N);
% Cyclic prefix
TxCy = [IFFT_Data((N-Ncp+1):N,:); IFFT_Data];       % Cyclic prefix
x=TxCy(:);

SNR = 0:2:16;
EsNo = SNR + 10*log10(Ts/Tsample);
EbNo = EsNo - 10*log10(k);

OFDM_BER_awgn=zeros(size(SNR));
OFDM_BER_ray=zeros(size(SNR));
OFDM_BER_rician=zeros(size(SNR));
OFDM_BER_mimo2v2=zeros(size(SNR));
OFDM_BER_mimo2v4=zeros(size(SNR));


% Rayleigh
ray_h=rayChan(ones((N+Ncp)*Nframes,1));
y_ray=x.*ray_h;

% Rician
rician_h=ricianChan(ones((N+Ncp)*Nframes,1));
y_rician=x.*rician_h;

%MIMO
x_OSTBC=OSTBCEnc(x);
[y_mimo2v2,pathGains2v2]=mimo2v2Chan(x_OSTBC);
[y_mimo2v4,pathGains2v4]=mimo2v4Chan(x_OSTBC);

for i = 1:length(SNR)
    
    awgnchan.SNR=SNR(i);
    % AWGN
    y_awgn = awgnchan(x);
    y_ray_ = awgnchan(y_ray);
    y_rician_ = awgnchan(y_rician);
    y_mimo2v2_ = awgnchan(y_mimo2v2);
    y_mimo2v4_ = awgnchan(y_mimo2v4);

    % equalization
    z_ray=y_ray_./ray_h;
    z_rician=y_rician_./rician_h;

    chEst2v2 = squeeze(sum(pathGains2v2,2));
    chEst2v4 = squeeze(sum(pathGains2v4,2));
    z_mimo2v2 = OSTBCComb2(y_mimo2v2_,chEst2v2);
    z_mimo2v4 = OSTBCComb4(y_mimo2v4_,chEst2v4);

    y_awgn=reshape(y_awgn,(N+Ncp),Nframes);
    z_ray=reshape(z_ray,(N+Ncp),Nframes);
    z_rician=reshape(z_rician,(N+Ncp),Nframes);
    z_mimo2v2=reshape(z_mimo2v2,(N+Ncp),Nframes);
    z_mimo2v4=reshape(z_mimo2v4,(N+Ncp),Nframes);

    % Remove CP
    y_awgn = y_awgn(Ncp+1:end,:);
    z_ray = z_ray(Ncp+1:end,:);
    z_rician = z_rician(Ncp+1:end,:);
    z_mimo2v2 = z_mimo2v2(Ncp+1:end,:);
    z_mimo2v4 = z_mimo2v4(Ncp+1:end,:);

    % Receiver
    R_Data_mod_awgn=((N-N/Np)/sqrt(N))*fft(y_awgn,N);
    R_Data_awgn=qamdemod(R_Data_mod_awgn,M);

    R_Data_mod_ray=((N-N/Np)/sqrt(N))*fft(z_ray,N);
    R_Data_ray=qamdemod(R_Data_mod_ray,M);
    R_Data_mod_rician=((N-N/Np)/sqrt(N))*fft(z_rician,N);
    R_Data_rician=qamdemod(R_Data_mod_rician,M);
    R_Data_mod_mimo2v2=((N-N/Np)/sqrt(N))*fft(z_mimo2v2,N);
    R_Data_mimo2v2=qamdemod(R_Data_mod_mimo2v2,M);
    R_Data_mod_mimo2v4=((N-N/Np)/sqrt(N))*fft(z_mimo2v4,N);
    R_Data_mimo2v4=qamdemod(R_Data_mod_mimo2v4,M);

    %end
    R_Data_awgn(1:Np:N,:)=0;
    R_Data_ray(1:Np:N,:)=0;
    R_Data_rician(1:Np:N,:)=0;
    R_Data_mimo2v2(1:Np:N,:)=0;
    R_Data_mimo2v4(1:Np:N,:)=0;
    
    Total_bits=(N-N/Np)*Nframes*k;
    OFDM_BER_awgn(i)=biterr(Data,R_Data_awgn,log2(M))./Total_bits;
    OFDM_BER_ray(i)=biterr(Data,R_Data_ray,log2(M))./Total_bits;
    OFDM_BER_rician(i)=biterr(Data,R_Data_rician,log2(M))./Total_bits;
    OFDM_BER_mimo2v2(i)=biterr(Data,R_Data_mimo2v2,log2(M))./Total_bits;
    OFDM_BER_mimo2v4(i)=biterr(Data,R_Data_mimo2v4,log2(M))./Total_bits;
end
th_ber=berawgn(EbNo,'qam',M);

figure
semilogy(SNR,th_ber,'--ob','linewidth',2);
hold on
grid on
semilogy(SNR,OFDM_BER_awgn,'--or','linewidth',2)
semilogy(SNR,OFDM_BER_ray,'--og','linewidth',2)
semilogy(SNR,OFDM_BER_rician,'--ok','linewidth',2)
semilogy(SNR,OFDM_BER_mimo2v2,':or','linewidth',2)
semilogy(SNR,OFDM_BER_mimo2v4,':ob','linewidth',2)
legend('4-QAM','OFDM(4-QAM) AWGN','OFDM(4-QAM) Rayleigh','OFDM(4-QAM) Rician','OFDM(4-QAM) Rayleigh MIMO2x2','OFDM(4-QAM) Rayleigh MIMO2x4')
title('OFDM BER vs SNR');
xlabel('SNR');
ylabel('BER');
axis([0 max(SNR) 1e-6 1])
