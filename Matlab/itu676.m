function [y, y_wv, y_ox]=itu676(f,p,T,rho)
%%
% Atenuación específica causada por gases
% De la recomendación 676-11
% p: Dry air pressure hPa
% t: Temperatura ºC
% e: Presión parcial de vapor de agua hPa
% f: frecuencia GHz

% Cargamos parámetros espectroscópicos
load('oxy676.mat');
load('wat676.mat');

theta=300./T;
e=rho.*T./216.7;
%p=P-e;


%
% Line strength Si
%
Sio=a1.*(1e-7).*p.*(theta.^3.0).*exp(a2.*(1-theta));
Siw=b1.*(1e-1).*e.*(theta.^3.5).*exp(b2.*(1-theta));
%
% Line shape factor
%
deltao=(a5+a6.*theta).*(1e-4).*(p+e).*theta.^(0.8);         % water-vapor partial pressure, e, is added
deltaw=0;
deltafo_sinDop=a3.*1e-4.*(p.*theta.^(0.8-a4)+1.1.*e.*theta);    % width of the oxygen line without Doppler broadening
deltafw_sinDop=b3.*1e-4.*(p.*theta.^b4+b5.*e.*theta.^b6);       % width of the water vapor line without Doppler broadening

deltafo=sqrt(deltafo_sinDop.^2+2.25e-6);                                                                  % Doppler broadening correction
deltafw=0.535.*deltafw_sinDop + sqrt(0.217.*deltafw_sinDop.^2+((2.1316e-12.*fw.^2)./theta));              % Doppler broadening correction

%
part1Fio=(deltafo-deltao.*(fo-f))./((fo-f).^2+deltafo.^2);
part2Fio=(deltafo-deltao.*(fo+f))./((fo+f).^2+deltafo.^2);
%
%
Fio=(f./fo).*(part1Fio+part2Fio);
%
part1Fiw=(deltafw-deltaw.*(fw-f))./((fw-f).^2+deltafw.^2);
part2Fiw=(deltafw-deltaw.*(fw+f))./((fw+f).^2+deltafw.^2);
%
Fiw=(f./fw).*(part1Fiw+part2Fiw);

% Dry air continuum N1Df ( N'D(f) )%%%%%
d=5.6e-4.*(p+e).*theta.^(0.8);                 % water-vapor partial pressure, e, is not used
N2Df=f.*p.*theta.^2.*(((6.14e-5)./(d.*(1+(f./d).^2)))+...
        ((1.4e-12.*p.*theta.^1.5)./(1+1.9e-5.*f.^1.5)));

%
%
N2f=sum(Fio.*Sio)+sum(Fiw.*Siw)+N2Df;       
%
y=0.1820.*f.*N2f;

% Se incluyen los nuevos parámetros de salida
y_ox = 0.1820 .* f.* (sum(Fio.*Sio) + N2Df);
y_wv = 0.1820 .* f.* (sum(Fiw.*Siw));
%y_wc = 0.1820 .* f.* (Fiw(35).*Siw(35));     % water vapor continuum contribution