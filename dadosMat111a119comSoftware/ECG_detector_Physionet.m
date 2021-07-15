%% ECG_detector_Physionet.m
% Disciplina: PTC3456 - Processamento de Sinais Biomédicos
% Grupo 7
close all;
clear;
clc;


%% Parâmetros
% Sinta-se livre para editar!
low = true; % habilita filtro passa-baixa
high = true; % habilita filtro passa-alta
classic = false; % true : habilita o uso de filtros clássico desenvolvidos por 
                % Pan-Tompkins para Fs = 200 Hz 
                % false : habilita o uso de filtros customizados para Fs =
                % 360 Hz
graphics = true; % habilita exibição de gráficos
N_int = 30*360/200; % número de elementos usados na integração
signal = 1; % sinal estudade (1 ou 2)
save = true; % salva tábela


%% Design de Filtros
Fs = 360;
f = logspace(0,3,100);
% Passa-baixa
fc_low = 20; % 20 Hz
if classic
    num_low = [1, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0, 1];
    den_low = 36 .* [1, -2, 1];
else
    num_low = ones(1,8);
    den_low = 8;
end
if graphics
    figure;
    h = freqz(num_low, den_low, f, Fs);
    h_abs = abs(h);
    h_phase = angle(h);
    subplot(2,1,1)
    semilogx(f, 20*log10(h_abs), 'DisplayName', 'magnitude')
    hold on
    semilogx([1, 1000], [-3, -3], 'k:', 'DisplayName', '-3 dB')
    hold on
    semilogx([fc_low, fc_low], [-120, 60], 'k--', 'DisplayName', 'frequência de corte')
    hold on
    semilogx([Fs/2, Fs/2], [-120, 60], 'm--', 'DisplayName', 'frequência de Nyquist')
    hold on
    semilogx([Fs, Fs], [-120, 60], 'r--', 'DisplayName', 'frequência de amostragem')
    hold off
    grid on
    xticks(cat(2, cat(2, 1:9, cat(2, 10:10:90, 100:100:Fs)), 400:600:1000))
    yticks(-2000:20:2000)
    ylim([-100, 20])
    ylabel('Magnitude (dB)')
    xlabel('Frequência (Hz)')
    title('Diagrama de Bode de filtro passa-baixa')
    legend('Location', 'southeast')
    
    subplot(2,1,2)
    semilogx(f, (180/pi) .* h_phase, 'DisplayName', 'fase')
    hold on
    semilogx([fc_low, fc_low], [-200, 200], 'k--', 'DisplayName', 'frequência de corte')
    hold on
    semilogx([Fs/2, Fs/2], [-200, 200], 'm--', 'DisplayName', 'frequência de Nyquist')
    hold on
    semilogx([Fs, Fs], [-200, 200], 'r--', 'DisplayName', 'frequência de amostragem')
    hold off
    grid on
    xticks(cat(2, cat(2, 1:9, cat(2, 10:10:90, 100:100:Fs)), 400:600:1000))
    yticks(-360:90:360)
    ylim([-200, 200])
    ylabel('Fase (degraus)')
    xlabel('Frequência (Hz)')
    legend('Location', 'southeast')
end

% Passa-alta
fc_high = 5; % 5 Hz
if classic
    num_high = zeros(1,33);
    num_high(1,1) = -1;
    num_high(1,17) = 32;
    num_high(1,18) = -32;
    num_high(1,33) = 1;
    den_high = 32 .* [1, -1];
else
    order = 4;
    filter_type = 'high';
    [num_high, den_high] = butter(order, 2*fc_high/Fs, filter_type);
end
if graphics
    figure;
    h = freqz(num_high, den_high, f, Fs);
    h_abs = abs(h);
    h_phase = angle(h);
    subplot(2,1,1)
    semilogx(f, 20*log10(h_abs), 'DisplayName', 'magnitude')
    hold on
    semilogx([1, 1000], [-3, -3], 'k:', 'DisplayName', '-3 dB')
    hold on
    semilogx([fc_high, fc_high], [-120, 60], 'k--', 'DisplayName', 'frequência de corte')
    hold on
    semilogx([Fs/2, Fs/2], [-120, 60], 'm--', 'DisplayName', 'frequência de Nyquist')
    hold on
    semilogx([Fs, Fs], [-120, 60], 'r--', 'DisplayName', 'frequência de amostragem')
    hold off
    grid on
    xticks(cat(2, cat(2, 1:9, cat(2, 10:10:90, 100:100:Fs)), 400:600:1000))
    yticks(-2000:20:2000)
    ylim([-100, 20])
    ylabel('Magnitude (dB)')
    xlabel('Frequência (Hz)')
    title('Diagrama de Bode de filtro passa-alta')
    legend('Location', 'southeast')
    
    subplot(2,1,2)
    semilogx(f, (180/pi) .* h_phase, 'DisplayName', 'fase')
    hold on
    semilogx([fc_high, fc_high], [-200, 200], 'k--', 'DisplayName', 'frequência de corte')
    hold on
    semilogx([Fs/2, Fs/2], [-200, 200], 'm--', 'DisplayName', 'frequência de Nyquist')
    hold on
    semilogx([Fs, Fs], [-200, 200], 'r--', 'DisplayName', 'frequência de amostragem')
    hold off
    grid on
    xticks(cat(2, cat(2, 1:9, cat(2, 10:10:90, 100:100:Fs)), 400:600:1000))
    yticks(-360:90:360)
    ylim([-200, 200])
    ylabel('Fase (degraus)')
    xlabel('Frequência (Hz)')
    legend('Location', 'southeast')
    
    num_res = conv(num_low, num_high);
    den_res = conv(den_low, den_high);
    figure;
    h = freqz(num_res, den_res, f, Fs);
    h_abs = abs(h);
    h_phase = angle(h);
    subplot(2,1,1)
    semilogx(f, 20*log10(h_abs), 'DisplayName', 'magnitude')
    hold on
    semilogx([1, 1000], [-3, -3], 'k:', 'DisplayName', '-3 dB')
    hold on
    semilogx([fc_low, fc_low], [-120, 60], 'k--', 'DisplayName', 'frequência de corte do passa-baixa')
    hold on
    semilogx([fc_high, fc_high], [-120, 60], 'k--', 'DisplayName', 'frequência de corte do passa-alta')
    hold on
    semilogx([Fs/2, Fs/2], [-120, 60], 'm--', 'DisplayName', 'frequência de Nyquist')
    hold on
    semilogx([Fs, Fs], [-120, 60], 'r--', 'DisplayName', 'frequência de amostragem')
    hold off
    grid on
    xticks(cat(2, cat(2, 1:9, cat(2, 10:10:90, 100:100:Fs)), 400:600:1000))
    yticks(-2000:20:2000)
    ylim([-100, 20])
    ylabel('Magnitude (dB)')
    xlabel('Frequência (Hz)')
    title('Diagrama de Bode de filtro passa-banda')
    legend('Location', 'southeast')
    
    subplot(2,1,2)
    semilogx(f, (180/pi) .* h_phase, 'DisplayName', 'fase')
    hold on
    semilogx([fc_low, fc_low], [-200, 200], 'k--', 'DisplayName', 'frequência de corte do passa-baixa')
    hold on
    semilogx([fc_high, fc_high], [-200, 200], 'k--', 'DisplayName', 'frequência de corte do passa-alta')
    hold on
    semilogx([Fs/2, Fs/2], [-200, 200], 'm--', 'DisplayName', 'frequência de Nyquist')
    hold on
    semilogx([Fs, Fs], [-200, 200], 'r--', 'DisplayName', 'frequência de amostragem')
    hold off
    grid on
    xticks(cat(2, cat(2, 1:9, cat(2, 10:10:90, 100:100:Fs)), 400:600:1000))
    yticks(-360:90:360)
    ylim([-200, 200])
    ylabel('Fase (degraus)')
    xlabel('Frequência (Hz)')
    legend('Location', 'southeast')
end

% Derivação
num_d_dt = [2, 1, 0, -1, -2];
den_d_dt = 8;

% Integração
num_int = ones(1,N_int);
den_int = N_int;


ann_total = 0;
pred_total = [0; 0];
avg_ecg = zeros(9,1);
std_ecg = zeros(9,1);
avg_ecg_norm = zeros(9,1);
std_ecg_norm = zeros(9,1);
avg_IRR = zeros(9,1);
std_IRR = zeros(9,1);
TP = zeros(9,1);
FP = zeros(9,1);
FN = zeros(9,1);
for n =111:119
    %% Importação de dados
    arqnum =sprintf('%3d.mat', n);
    load(arqnum);      
    % save(arq_mat,'ecgs','ts','Fs','ann','type');  Ver documentacao abaixo sobre as variáveis
    % ecgs(Nt,2):[double] matriz com Nt linhas por 2 colunas. Cada coluna contém sinal com Nt amostras já calibradas em mV
    % ts(Nt,1): [double] vetor com Nt instantes de tempo em signal
    % Fs: frequencia de amostragem em Hz
    % ann(Na,1):[double] vetor contendo a posicao (em amostra) da anotacao. Tem 'Na' anotacoes. Ex.: suponha que a terceira 
    %       anotacao ocorreu na amostra 270, entao p=ann(3) retornará p=270.
    % type(Na,1):[char] vetor com a anotacao para cada uma das 'Na' anotacoes.
    N = size(ecgs,1);
    fprintf('%s sinal com %d amostras\n', arqnum, N);
    Ntypes = numel(type);
    delay = 0;
    
    
    %% Média e desvio-padrão do sinal original
    avg_ecg(n-110,1) = mean(ecgs(:,signal));
    std_ecg(n-110,1) = std(ecgs(:,signal));
    
    
    %% Filtragem
    % Passa-baixa
    if low
        ecgs_filt = filter(num_low, den_low, ecgs);
    else
        ecgs_filt = ecgs;
    end
    % Passa-alta
    if high
        ecgs_filt = filter(num_high, den_high, ecgs_filt);
    end
    
    
    %% Normalização
    ecgs_norm = ecgs_filt ./ abs(max(ecgs_filt));
    
    
    %% Média e desvio-padrão do sinal normalizado
    avg_ecg_norm(n-110,1) = mean(ecgs_norm(:,signal));
    std_ecg_norm(n-110,1) = std(ecgs_norm(:,signal));


    %% Derivada
    decgs_dt = filter(num_d_dt, den_d_dt, ecgs_norm);
    
    
    %% Quadrado
    decgs_dt_2 = decgs_dt .^ 2;
    
    
    %% Integração
    ecgs_int = filter(num_int, den_int, decgs_dt_2);
    
    
    %% Limiar Adaptativo
    pred = zeros(N,1);
    pred_size = 0;
    IRR = zeros(N,1);
    % Inicialização
    SPKI = prctile(ecgs_int(1:N,signal), 95); % estimação dos valores máximos de R
    NPKI = 0; % estimação do valor máximo do ruído
    THRESHOLD_I1 = SPKI/2; % limiar de detecção de onda R
    THRESHOLD_I2 = THRESHOLD_I1 / 2; % limiar retroativo de detecção de onda R
    RR_AVERAGE_1 = Fs;
    RR_AVERAGE_2 = Fs;
    RR_LOW_LIMIT = Fs/4;
    RR_HIGH_LIMIT = 2*Fs;
    RR_MISSED_LIMIT = 2*Fs;
    % Loop
    for i=3:N
        if (ecgs_int(i,signal) - ecgs_int(i-1,signal))*(ecgs_int(i-1,signal) - ecgs_int(i-2,signal)) < 0 % mudança de sinal
            if pred_size > 1
                if (ecgs_int(i,signal) > THRESHOLD_I1) && (i - pred(pred_size,1) > RR_LOW_LIMIT)% QRS encontrado
                    SPKI = 0.125*ecgs_int(i,signal) + 0.875*SPKI;
                    pred(pred_size+1,1) = i;
                    pred_size = pred_size + 1;
                    IRR(pred_size-1,1) = pred(pred_size,1) - pred(pred_size-1,1);
                else % ruído encontrado
                    NPKI = 0.125*ecgs_int(i,signal) + 0.875*NPKI;
                end
            else
                if ecgs_int(i,signal) > THRESHOLD_I1 % QRS encontrado
                    SPKI = 0.125*ecgs_int(i,signal) + 0.875*SPKI;
                    pred(pred_size+1,1) = i;
                    pred_size = pred_size + 1;
                else % ruído encontrado
                    NPKI = 0.125*ecgs_int(i,signal) + 0.875*NPKI;
                end
            end
            THRESHOLD_I1 = NPKI + 0.25*(SPKI - NPKI);
            THRESHOLD_I2 = 0.5*THRESHOLD_I1;
        end
        if pred_size > 8
            RR_AVERAGE_1 = sum(IRR(pred_size-8:pred_size-1,1)) / 8;
            IRR_ok = IRR((IRR(1:pred_size-1,1) >= RR_LOW_LIMIT) & (IRR(1:pred_size-1,1) <= RR_HIGH_LIMIT));
            if length(IRR_ok) > 8
                RR_AVERAGE_2 = sum(IRR_ok(length(IRR_ok)-7:length(IRR_ok),1)) / 8;
                RR_LOW_LIMIT = 0.92*RR_AVERAGE_2;
                RR_HIGH_LIMIT = 1.16*RR_AVERAGE_2;
                RR_MISSED_LIMIT = 1.66*RR_AVERAGE_2;
                if i - pred(pred_size) > RR_MISSED_LIMIT
                    % Implementar detecção retroativa
                    [peak, peak_i] = max(ecgs_int(pred(pred_size)+round(RR_LOW_LIMIT):i,signal));
                    peak_i = peak_i + pred(pred_size) + round(RR_LOW_LIMIT) - 1;
                    SPKI = 0.125*peak + 0.875*SPKI;
                    pred(pred_size+1,1) = peak_i;
                    pred_size = pred_size + 1;
                    IRR(pred_size-1,1) = pred(pred_size,1) - pred(pred_size-1,1);
                    THRESHOLD_I1 = NPKI + 0.25*(SPKI - NPKI);
                    THRESHOLD_I2 = 0.5*THRESHOLD_I1;
                end
            end
        end
    end
    
    
    %% Avaliação
    pred = pred(pred(:,1) > 0,1);
    IRR = IRR(IRR(:,1) > 0, 1);
    avg_IRR(n-110,1) = mean(IRR);
    std_IRR(n-110,1) = std(IRR);
    ann_total = ann_total + length(ann); 
    pred_total(1) = pred_total(1) + length(pred); 
    j = 1;
    fn = true;
    if pred(j,1) - delay <= (ann(1) + ann(2))/2
        TP(n-110,1) = TP(n-110,1) + 1;
        j = j + 1;
        fn = false; 
    end
    for i=2:length(ann)-1
        if fn
            FN(n-110,1) = FN(n-110,1) + 1;
        end
        fn = true;
        while (j <= pred_size) && (pred(min(j,pred_size),1) - delay <= (ann(i) + ann(i+1))/2) && fn
            if pred(j,1) - delay >= (ann(i) + ann(i-1))/2
                TP(n-110,1) = TP(n-110,1) + 1;
                fn = false;
            else
                FP(n-110,1) = FP(n-110,1) + 1;
            end
            j = j + 1;
        end
    end
    i = length(ann);
    if (j > pred_size)
        FN(n-110,1) = FN(n-110,1) + 1;
    else
        while (j <= pred_size) && (pred(j,1) - delay < (ann(length(ann)) + ann(length(ann)-1))/2)
            FN(n-110,1) = FN(n-110,1) + 1;
            j = j + 1;
        end
        if j > pred_size
            FN(n-110,1) = FN(n-110,1) + 1;
        else
            TP(n-110,1) = TP(n-110,1) + 1;
            FP(n-110,1) = FP(n-110,1) + pred_size - j;
        end
    end
    
    if length(ann) ~= TP(n-110,1) + FN(n-110,1)
        disp('AVISO: Há erros no sistema de avaliação !!!')
        disp(['length(ann) == TP + FN : ', num2str(length(ann) == TP(n-110,1) + FN(n-110,1)), ' : ', num2str(length(ann)), ' == ', num2str(TP(n-110,1) + FN(n-110,1))])
    end
    if pred_size ~= TP(n-110,1) + FP(n-110,1)
        disp('AVISO: Há erros no sistema de avaliação !!!')
        disp(['pred_size == TP + FP : ', num2str(pred_size == TP(n-110,1) + FP(n-110,1)), ' : ', num2str(pred_size), ' == ', num2str(TP(n-110,1) + FP(n-110,1))])
    end
    
    
    %% Visualização
    if graphics
        figure;
        max_ecg = max(ecgs(:,signal))*ones(Ntypes,1);
        subplot(2,1,1);
        plot(ts(1:N), ecgs(1:N,signal));
        hold on;
        plot(ts(ann(ann<N)+1), ecgs(ann(ann<N)+1,signal), 'ro');
        hold on;
        plot(ts(pred(pred<N)-delay), ecgs(pred(pred<N, 1)-delay,signal), 'mp');
        hold on;
        text(ts(ann(ann<N)+1), max_ecg, type);                    
        hold off;
        ylabel('mV');
        xlabel('signal');
        grid on;
        title([arqnum, ' sinal ', num2str(signal)]);

        subplot(2,1,2);
        plot(ts(1:N), ecgs_int(1:N,signal));
        hold on;
        plot(ts(pred(pred < N)), ecgs_int(pred(pred < N, 1),signal), 'ro');
        hold off;
        xlabel('s');
        grid on;
        title([arqnum, ' sinal ', num2str(signal), ' Pan-Thompkins']);
        drawnow;
    end
end


%% Estatísticas
% Tabela
casos = {'111.mat'; '112.mat'; '113.mat'; '114.mat'; '115.mat'; '116.mat'; '117.mat'; '118.mat'; '119.mat'; 'media'; 'dp'};
media_original = zeros(11,1);
media_original(1:9,1) = avg_ecg(:,1);
media_original(10,1) = mean(media_original(1:9,1));
media_original(11,1) = std(media_original(1:9,1));
dp_original = zeros(11,1);
dp_original(1:9,1) = std_ecg(:,1);
dp_original(10,1) = mean(dp_original(1:9,1));
dp_original(11,1) = std(dp_original(1:9,1));
media_normalizado = zeros(11,1);
media_normalizado(1:9,1) = avg_ecg_norm(:,1);
media_normalizado(10,1) = mean(media_normalizado(1:9,1));
media_normalizado(11,1) = std(media_normalizado(1:9,1));
dp_normalizado = zeros(11,1);
dp_normalizado(1:9,1) = std_ecg_norm(:,1);
dp_normalizado(10,1) = mean(dp_normalizado(1:9,1));
dp_normalizado(11,1) = std(dp_normalizado(1:9,1));
media_IRR_s = zeros(11,1);
media_IRR_s(1:9,1) = avg_IRR(:,1)./Fs;
media_IRR_s(10,1) = mean(media_IRR_s(1:9,1));
media_IRR_s(11,1) = std(media_IRR_s(1:9,1));
dp_IRR_s = zeros(11,1);
dp_IRR_s(1:9,1) = std_IRR(:,1)./Fs;
dp_IRR_s(10,1) = mean(dp_IRR_s(1:9,1));
dp_IRR_s(11,1) = std(dp_IRR_s(1:9,1));
bpm = media_IRR_s*60;
TP_aux = TP(:,1);
TP = zeros(11,1);
TP(1:9,1) = TP_aux(:,1);
TP(10,1) = mean(TP(1:9,1));
TP(11,1) = std(TP(1:9,1));
FP_aux = FP(:,1);
FP = zeros(11,1);
FP(1:9,1) = FP_aux(:,1);
FP(10,1) = mean(FP(1:9,1));
FP(11,1) = std(FP(1:9,1));
FN_aux = FN(:,1);
FN = zeros(11,1);
FN(1:9,1) = FN_aux(:,1);
FN(10,1) = mean(FN(1:9,1));
FN(11,1) = std(FN(1:9,1));
Nv = zeros(11,1); % número de QRSs verdadeiros
Nv(1:9,1) = TP(1:9,1) + FN(1:9,1);
Nv(10,1) = mean(Nv(1:9,1));
Nv(11,1) = std(Nv(1:9,1));
Nd = zeros(11,1); % número de QRSs detectados
Nd(1:9,1) = TP(1:9,1) + FP(1:9,1);
Nd(10,1) = mean(Nd(1:9,1));
Nd(11,1) = std(Nd(1:9,1));
prct_FP = zeros(11,1); % porcentagem de falso-positivos
prct_FP(1:9,1) = 100 .* FP(1:9,1) ./ Nv(1:9,1);
prct_FP(10,1) = mean(prct_FP(1:9,1));
prct_FP(11,1) = std(prct_FP(1:9,1));
prct_FN = zeros(11,1); % porcentagem de falso-negativos
prct_FN(1:9,1) = 100 .* FN(1:9,1) ./ Nv(1:9,1);
prct_FN(10,1) = mean(prct_FN(1:9,1));
prct_FN(11,1) = std(prct_FN(1:9,1));

res_table = table(casos, media_original, dp_original, ...
    media_normalizado, dp_normalizado, Nv, ...
    Nd, FP, prct_FP, FN, prct_FN, media_IRR_s, dp_IRR_s, bpm);

if save
    writetable(res_table, ['table_signal_', num2str(signal),'.csv']);
end

fprintf('\n Fim \n');

