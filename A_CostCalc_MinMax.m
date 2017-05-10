%HEAT PUMP
capacity= [1500,25000]
PHP =14901.0.*(capacity./95).^0.9
%BOILER
duty = [1000,25000]
%k1 = 2.291;
%k2 = 0.738;
PB =10.^(k1 + k2*log10(duty))
%SOFC
Cap = [1000,10000]
PSOFC = ((1.9*0.95)^.95)*1.02*1000*Cap
%ICE
%k1 = 2.7635;
%k2 = 0.8574;
%k3 = -0.0098;
CAPICE = [100,15000]
PICE = 10.^(k1 + k2*log10(CAPICE) + k3.*(log10(CAPICE)).^2)
PICE2 = -0.0266.*(3*CAPICE/4).^2 + 578.84*(3*CAPICE/4) + 208174
%Turbine
CAPGT = [100,15000]
k1=3.4171;
k2=0.6112;
PGT = 10.^(k1 + k2.*log10(CAPGT))