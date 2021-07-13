% sf_ler_mat.m  Programa para ler .mat e display 
% EPUSP/PTC 5/06/2021 SF
%
%echo on
close all;
display('Exemplo de leitura dos sinais de ECG no formato .mat')
display('-Cada arquivo .mat(paciente) tem 2 ECGs com 2 min de duracao')
display('-Os sinais já estão calibrados em mV e o tempo em s')
display('-A anotacao de cada evento foi tambem colocada nas figuras')
Fa          =360;  %Hz
for n =111:119
    arqnum =sprintf('%3d.mat',n);
    load (arqnum);      %    save(arq_mat,'ecgs','ts','Fs','ann','type');  Ver documentacao abaixo sobre as variáveis
       % ecgs(Nt,2):[double] matriz com Nt linhas por 2 colunas. Cada coluna contém sinal com Nt amostras já calibradas em mV
       % ts(Nt,1): [double] vetor com Nt instantes de tempo em s
       % Fs: frequencia de amostragem em Hz
       % ann(Na,1):[double] vetor contendo a posicao (em amostra) da anotacao. Tem 'Na' anotacoes. Ex.: suponha que a terceira 
       %       anotacao ocorreu na amostra 270, entao p=ann(3) retornará p=270.
       % type(Na,1):[char] vetor com a anotacao para cada uma das 'Na' anotacoes.      
    N   =size(ecgs,1);
    fprintf('\n %s com %d amostras',arqnum,N);
    Ntypes=numel(type);
    %Plot 2D version of signal and labels
    figure; 
    maxEcg1 =max(ecgs(:,1))*ones(Ntypes,1);
    subplot(2,1,1); plot(ts(1:N),ecgs(1:N,1));ylabel('mV'); xlabel('s'); grid on;title([arqnum ' sinal 1']); hold on;
                    plot(ts(ann(ann<N)+1),ecgs(ann(ann<N)+1),'ro');
                    text(ts(ann(ann<N)+1),maxEcg1,type);                    
    subplot(2,1,2); plot(ts(1:N),ecgs(1:N,2));ylabel('mV'); xlabel('s'); grid on;title([arqnum ' sinal 2' ]);
    drawnow;
end
fprintf('\n Fim \n');
