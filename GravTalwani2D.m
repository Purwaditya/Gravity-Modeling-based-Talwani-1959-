clc;close all;clear

G=6.673*10^-11;
z0=0; %topografi
dx=50; % besar grid arah x
dz=50; % besar grid arah z
x=0:50:2000; % titik observasi
sudut=4;

%% membuat model
blokX=0:dx:2000; %banyak blok arah X
blokZ=0.1:dz:1000; %banyak blok arah Z
density=zeros(length(blokZ)-1,length(blokX)-1);
density(8:13,10:18)=1000; %model dalam bentuk matriks 19x40
density(8:13,19:25)=1500;
densv1=reshape(density,760,[]); %model dalam bentuk vektor 1D (m)
%% menghitung gobs
gobs=zeros(1,length(x));
kernel=zeros(length(x),(length(blokZ)-1)*(length(blokX)-1));
for p=1:length(x);
    total=0;
    r=1;
    for i=1:length(blokX)-1
        for j=1:length(blokZ)-1
            xT=[blokX(i);blokX(i+1);blokX(i+1);blokX(i)];
            zT=[blokZ(j);blokZ(j);blokZ(j+1);blokZ(j+1)];
            
            kernel(p,r)=2*G*Talwani(x(p),z0,xT,zT,sudut); %matrik kernel (G)
            r=r+1;
        end
    end
    anomali=kernel*densv1; %d=Gm
end
%% Figure
figure(1)
subplot(2,1,1)
plot(x,anomali,'xb')
xlabel('Jarak(m)','FontWeight', 'bold', 'FontSize', 12)
ylabel('Gravitasi Observasi(mGal)','FontWeight', 'bold', 'FontSize', 12)

subplot(2,1,2)
for i=1:length(blokX)-1
        for j=1:length(blokZ)-1
            xf=[blokX(i);blokX(i+1);blokX(i+1);blokX(i)];
            zf=[blokZ(j);blokZ(j);blokZ(j+1);blokZ(j+1)];
            C=density(j,i);
            patch(xf,zf,C)
        end
end
xlim([0 max(x)])
ylim([0 max(zf)])
set(gca,'ydir','reverse')
xlabel('Jarak(m)','FontWeight', 'bold', 'FontSize', 12)
ylabel('Kedalaman(m)','FontWeight', 'bold', 'FontSize', 12)
colormap(jet)

% scale bar
cb = colorbar('southoutside');
cb.Label.String = 'Densitas (kg/m^3)';
cb.Label.FontWeight = 'bold';
cb.Label.FontSize = 12;