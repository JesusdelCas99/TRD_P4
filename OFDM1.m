clc
clear variables
close all

% Empleamos latex como interprete de texto en las figuras
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%Cargamos el fichero OFDM1
fileID = fopen('OFDM1','r');
newrow=1;
datareal=fread(fileID,'double');
Ndatareal=length(datareal);

%Agrupamos los valores devuelto por el vector real 'Ndatareal' en parejas
%de valores complejos
for row=1:2:Ndatareal
    data(newrow)=complex(datareal(row),datareal(row+1));
    newrow=newrow+1;
end

%Cerramos el fichero fileID
fclose(fileID);

%% Cálculo del MAXIMUM CORRELATION METRIC METHOD

%Número de portadoras teóricas del modelo que queremos chequar
NFFT_K_=[2048 8192];

PNfft=zeros(1,2);
ANfft=zeros(1,2);

for m=1:length(NFFT_K_)

    NFFT_K=NFFT_K_(m);
    
    % Tomamos ND muestras de señal OFDM que comprendan un número suficiente de
    % simbolos OFDM
    Nd=3*NFFT_K;

    %El tamaño de la ventana de correlacion Nw se elige de forma que sea igual
    %o menor que el menor intervalo de guardo previsto en los simbolos OFDM
    Nw=NFFT_K*(1/32); %Valor mas pequeño del intervalo de guarda en DVB-T

    %Calculamos 'Correlation Metric Method'
    rho=Max_Corr_Method(Nw,data,Nd,NFFT_K);

    if m==1
        rho1=rho;
    else
        rho2=rho;
    end

    %Calculo de la potencia y la media
    PNfft1=(1/Nd)*sum(abs(rho).*abs(rho));
    ANfft1=(1/Nd)*sum(abs(rho));
    
    PNfft(m)=PNfft1;
    ANfft(m)=ANfft1;

end

%Graficamos la funcion rho 
figure(1)
plot(abs(rho1))
hold on
plot(abs(rho2))
xlabel('muestra (n)')
title('Funcion metrica de correlacion, Senal OFDM1')
legend('$NFFt=2048$','$NFFt=8192$')
grid on;

%Calculamos la variace-to-average-ratio metric
NFFT_1=(PNfft(1)-ANfft(1)^(2))/ANfft(1);
NFFT_2=(PNfft(2)-ANfft(2)^(2))/ANfft(2);

max_OFDM1=max(NFFT_1,NFFT_2); %Vemos que el modo de trabajo es 2K

if NFFT_1>NFFT_2
    fprintf('\nEl modo óptimo es 2K')    
else
    fprintf('\nEl modo óptimo es 8K')
end

%Para el fichero OFDM1 hemos obtenido el que el modo es 2k


%% Estimación de la duración del intervalo de guarda en DVB-T

%Numero de portadoras en 2k
NFFT_2K=2048;

%ND muestras de señal OFDM que comprendan un número suficiente de simbolos OFDM
Nd=3*NFFT_2K;

%Duracion estimada en numero de muestras del intervalo de guarda Ng
[Ng_est,thr,rho_max] = guard_interval(Nd,rho1,NFFT_2K);

%Graficamos las funciones
figure(2)
plot(abs(rho1)/rho_max)
hold on
plot(thr)
hold off
grid on;
xlabel('muestra (n)')
legend('$|\rho_{\hat{N}_{FFT,N_{W}}}(n)|/\rho^{\prime}_{\hat{N}_{FFT,N_{W}}}$','$thr(|\rho_{\hat{N}_{FFT,N_{W}}}(n)|,\acute{\rho}_{\hat{N}_{FFT,N_{W}}}/3)$')
title('Senal OFDM1','interpreter','latex')

%% Sincronización temporal en DVB-T

len_OFDM=NFFT_2K+Ng_est; %Longitud del simbolo OFDM
len_total=length(data);  %Longitud total del vector de simbolos PSK

rho=zeros(1,len_OFDM);

for i=0:len_OFDM-1
    suma=0;
    for j=0:Ng_est-1
        suma=conj(data(j+i+1))*data(j+i+NFFT_2K+1)+suma;
    end
    rho(i+1)=(suma)/Nw;

end

rho_abs=abs(rho);

%Buscamos el valor maximo así como su indice
max_rho=find(rho_abs==max(rho_abs));
maximo=max(rho_abs);

%Generamos el eje x para nuestra grafica, que seran los indices 'n'
n=0:1:length(rho_abs)-1;

%Graficamos los resultados obtenidos
figure(87)
hold on
stem(n(max_rho),rho_abs(max_rho),'o')
plot(n,rho_abs)
hold off
grid on;
title('Maximum correlation metric method')
xlabel('Muestra (n)')
ylabel('$|\rho_{\hat{N}_{FFT,N_{W}}}(n)|$')
legend('Muestra inicial del prefijo ciclico')
ini_n=max_rho-1+Ng_est;
sprintf("El primer simbolo OFDM comienza en la muestra "+ini_n)

%% Segmentamos la señal en simbolos OFDM
ini_n=ini_n+1;%Devolvemos a la muestra inicial del bucle
n_sim=1;
l_archivo=length(data);

while((ini_n+(n_sim-1)*(NFFT_2K+Ng_est)<=l_archivo) && (ini_n+(n_sim-1)*(NFFT_2K+Ng_est)+NFFT_2K-1<=l_archivo))
    if n_sim==1
        sim_vec=data(ini_n:ini_n+NFFT_2K-1);
        ini_pos=ini_n;
        end_pos=ini_n+NFFT_2K-1;
    else
        sim_vec=vertcat(sim_vec,data(ini_n+(n_sim-1)*(NFFT_2K+Ng_est):ini_n+(n_sim-1)*(NFFT_2K+Ng_est)+NFFT_2K-1));
        ini_pos=horzcat(ini_pos,ini_n+(n_sim-1)*(NFFT_2K+Ng_est));
        end_pos=horzcat(end_pos,ini_n+(n_sim-1)*(NFFT_2K+Ng_est)+NFFT_2K-1);
    end
    n_sim=n_sim+1;
end

ini_pos=ini_pos-1;
end_pos=end_pos-1;

sprintf("Numero completo de simbolos OFDM: "+(n_sim-1))