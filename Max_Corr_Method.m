function rho=Max_Corr_Method(Nw,r,Nd,Nfft)

rho=zeros(1,Nd);

for i=0:Nd-1
    suma=0;
    for j=0:Nw-1
        suma=conj(r(j+i+1))*r(j+i+Nfft+1)+suma;
    end
    rho(i+1)=(suma)/Nw;

end
end