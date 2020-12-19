function [Ng_est,thr,rho_max] = guard_interval(Nd,rho,NFFT_2K)

    rho_max=max(abs(rho)); %Cálculo de la relación eq(4)
    
    suma=0; %Inicializamos la sumatorio del la eq(5) a cero
    
    %Funcion metrica de correlacion con el valor NFFt estimado, eq(5)
    thr=zeros(1,Nd);
    
    for i=1:Nd
        if abs(rho(i))>(rho_max/3)
            thr(i)=1;
        else
            thr(i)=0;
        end
        suma=sum(thr);
    end
    
    sigma=(1/Nd)*suma;
    
    %Calculo de la duracion estimada en número de muestras del intervalo de
    %guarda
    
    Ng=NFFT_2K*[1/32 1/16 1/8 1/4]; %Posibles duraciones del intervalo de 
                                    %guarda
    
    Ng_prima=zeros(1,length(Ng));   %Vector que alojara las estimaciones de
                                    %la duracion del intervalo de guarda
    
    for i=1:length(Ng)
        Ng_prima(i)=abs(sigma-((Ng(i))/(Ng(i)+NFFT_2K)));
    end
    
    %Estimación del argumento que minimiza la ecuación, eq(7)
    [~,pos]=min(Ng_prima);
    Ng_est=Ng(pos);
    
end
