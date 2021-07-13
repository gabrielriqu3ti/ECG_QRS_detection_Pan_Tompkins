%% ECG_detector_Physionet.m
% Disciplina: PTC3456 - Processamento de Sinais Biomédicos
% Grupo 7
close all;
clear;

%% Parâmetros
low = true;
high = true;


%% Design de Filtros
% Passa-baixa
num_low = [1, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0, 1];
den_low = 36 .* [1, -2, 1];

% Passa-alta
num_high = zeros(1,33);
num_high(1,1) = -1;
num_high(1,17) = 32;
num_high(1,18) = -32;
num_high(1,33) = 1;
den_high = 32 .* [1, -1];


avg_ecg = zeros(9,2);
std_ecg = zeros(9,2);
avg_ecg_norm = zeros(9,2);
std_ecg_norm = zeros(9,2);
n = 111;
for n =111:119
% if true
    %% Importação de dados
    arqnum =sprintf('%3d.mat', n);
    load(arqnum);      
    % save(arq_mat,'ecgs','ts','Fs','ann','type');  Ver documentacao abaixo sobre as variáveis
    % ecgs(Nt,2):[double] matriz com Nt linhas por 2 colunas. Cada coluna contém sinal com Nt amostras já calibradas em mV
    % ts(Nt,1): [double] vetor com Nt instantes de tempo em s
    % Fs: frequencia de amostragem em Hz
    % ann(Na,1):[double] vetor contendo a posicao (em amostra) da anotacao. Tem 'Na' anotacoes. Ex.: suponha que a terceira 
    %       anotacao ocorreu na amostra 270, entao p=ann(3) retornará p=270.
    % type(Na,1):[char] vetor com a anotacao para cada uma das 'Na' anotacoes.
    N = size(ecgs,1);
    fprintf('\n %s com %d amostras', arqnum, N);
    Ntypes = numel(type);
    
    
    %% Média e desvio-padrão do sinal original
    avg_ecg(n-110,:) = mean(ecgs);
    std_ecg(n-110,:) = std(ecgs);
    
    
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
    avg_ecg_norm(n-110,:) = mean(ecgs_norm);
    std_ecg_norm(n-110,:) = std(ecgs_norm);

    %% Predição
    
    
    %% Visualização
    figure;
    max_ecg = max(ecgs(:,1))*ones(Ntypes,1);
    subplot(4,1,1);
    plot(ts(1:N), ecgs(1:N,1));
    ylabel('mV');
    xlabel('s');
    grid on;
    title([arqnum, ' sinal 1']);
    hold on;
    plot(ts(ann(ann<N)+1), ecgs(ann(ann<N)+1), 'ro');
    text(ts(ann(ann<N)+1), max_ecg, type);                    
    
    subplot(4,1,2);
    max_ecg_norm = max(ecgs_norm(:,1))*ones(Ntypes,1);
    plot(ts(1:N), ecgs_norm(1:N,1));
    ylabel('mV');
    xlabel('s');
    grid on;
    title([arqnum, ' sinal 1 Pan-Thompkins']);
    hold on;
    plot(ts(ann(ann<N)+1+16+5), ecgs_norm(ann(ann<N)+1+16+5), 'ro');
    text(ts(ann(ann<N)+1+16+5), max_ecg_norm, type);                    
    
    subplot(4,1,3);
    plot(ts(1:N), ecgs(1:N,2));
    ylabel('mV');
    xlabel('s');
    grid on;
    title([arqnum, ' sinal 2']);
    
    subplot(4,1,4);
    plot(ts(1:N), ecgs_norm(1:N,2));
    ylabel('mV');
    xlabel('s');
    grid on;
    title([arqnum, ' sinal 2 Pan-Thompkins']);
    drawnow;
end
fprintf('\n Fim \n');

