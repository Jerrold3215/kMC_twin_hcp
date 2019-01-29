close all;
clear all;
clc;
% for T=100:100:1000;
%     T_hcp(T/100,1)=T;
%     T_twin(T/100,1)=T;
% for Stress=0.1:0.1:2.0;
%     S_hcp(round(Stress/0.1),1)=Stress;
%     S_twin(round(Stress/0.1),1)=Stress;
%     Temp(round(Stress/0.1),1)=T;
Stress=5;
% S=1*10^12;
S=Stress*10^12;
T=300;
K=1;
N=30;
time=100;
L=zeros(3*N,6);
C_hcp=zeros(50,1);
C_twin=zeros(time,K);
G_eV_list=zeros(10,1);
for i=1:1:N; 
    L(3*i,1)=3;
    L(3*i-1,1)=2;
    L(3*i-2,1)=1;
end
for i=1:1:time;
    C_hcp(i,K)=0;
    C_twin(i,K)=0;
end
a=3.57;
Tm=1650;
H=1*10^-6;
R=1*10^-9;
alpha=0.123;
miu=80;
R=10^-6;
burger=sqrt(1/6)*a*10^-10;
chi=alpha*(miu*10^12)*(burger)^2;
Eb1=350; Eb2=310; Eb3=339; Eb4=345; Eb5=339;
Er1=385; Er2=428; Er3=286; Er4=282; Er5=477;
de1=-35; de2=-183; de3=18; de4=-120;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rc1=sqrt(-chi/(pi*S-pi*350/burger));
% rc2=sqrt(-chi/(pi*S-pi*275/burger));
% rc3=sqrt(-chi/(pi*S-pi*304/burger));
% rc4=sqrt(-chi/(pi*S-pi*162/burger));
% rc5=sqrt(-chi/(pi*S-pi*357/burger));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rc1=chi/(-350+S*burger);
rc2=chi/(-310+S*burger);
rc3=chi/(-339+S*burger);
rc4=chi/(-345+S*burger);
rc5=chi/(-339+S*burger);
rc1r=chi/(-385+S*burger);
rc2r=chi/(-458+S*burger);
rc3r=chi/(-286+S*burger);
rc4r=chi/(-282+S*burger);
rc5r=chi/(-477+S*burger);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rc1=(Eb1 + (Eb1^2 + 6*S*chi)^(1/2))/(3*S);
% rc2=(Eb2 + (Eb2^2 + 6*S*chi)^(1/2))/(3*S);
% rc3=(Eb3 + (Eb3^2 + 6*S*chi)^(1/2))/(3*S);
% rc4=(Eb4 + (Eb4^2 + 6*S*chi)^(1/2))/(3*S);
% rc5=(Eb5 + (Eb5^2 + 6*S*chi)^(1/2))/(3*S);
% rc1r=(Er1 + (Er1^2 + 6*S*chi)^(1/2))/(3*S);
% rc2r=(Er2 + (Er2^2 + 6*S*chi)^(1/2))/(3*S);
% rc3r=(Er3 + (Er3^2 + 6*S*chi)^(1/2))/(3*S);
% rc4r=(Er4 + (Er4^2 + 6*S*chi)^(1/2))/(3*S);
% rc5r=(Er5 + (Er5^2 + 6*S*chi)^(1/2))/(3*S);
E1=pi*R*R*(de1*(1-(T/Tm))-burger*S)*6.242e+15;
E2=pi*R*R*(de2*(1-(T/Tm))-burger*S)*6.242e+15;
E3=pi*R*R*(de3*(1-(T/Tm))-burger*S)*6.242e+15;
E4=pi*R*R*(de4*(1-(T/Tm))-burger*S)*6.242e+15;
E_eV_list=[E1 E2 E3 E4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G1=pi*rc1*rc1*(Eb1*(1-(T/Tm))-burger*S)+2*pi*rc1*chi;
G2=pi*rc2*rc2*(Eb2*(1-(T/Tm))-burger*S)+2*pi*rc2*chi;
G3=pi*rc3*rc3*(Eb3*(1-(T/Tm))-burger*S)+2*pi*rc3*chi;
G4=pi*rc4*rc4*(Eb4*(1-(T/Tm))-burger*S)+2*pi*rc4*chi;
G5=pi*rc5*rc5*(Eb5*(1-(T/Tm))-burger*S)+2*pi*rc5*chi;
% G1r=pi*rc1r*rc1r*(Er1*(1-(T/Tm))-burger*S)+2*pi*rc1r*chi;
% G2r=pi*rc2r*rc2r*(Er2*(1-(T/Tm))-burger*S)+2*pi*rc2r*chi;
% G3r=pi*rc3r*rc3r*(Er3*(1-(T/Tm))-burger*S)+2*pi*rc3r*chi;
% G4r=pi*rc4r*rc4r*(Er4*(1-(T/Tm))-burger*S)+2*pi*rc4r*chi;
% G5r=pi*rc5r*rc5r*(Er5*(1-(T/Tm))-burger*S)+2*pi*rc5r*chi;
G1r=pi*rc1*rc1*(Eb1*(1-(T/Tm))-burger*S)+2*pi*rc1*chi-E1;
G2r=pi*rc2*rc2*(Eb2*(1-(T/Tm))-burger*S)+2*pi*rc2*chi-E2;
G3r=pi*rc3*rc3*(Eb3*(1-(T/Tm))-burger*S)+2*pi*rc3*chi-E3;
G4r=pi*rc4*rc4*(Eb4*(1-(T/Tm))-burger*S)+2*pi*rc4*chi-E4;
G5r=pi*rc5*rc5*(Eb5*(1-(T/Tm))-burger*S)+2*pi*rc5*chi-E4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G_eV_list(1)=G1*6.242e+15;
G_eV_list(2)=G2*6.242e+15;
G_eV_list(3)=G3*6.242e+15;
G_eV_list(4)=G4*6.242e+15;
G_eV_list(5)=G5*6.242e+15;
G_eV_list(6)=G1r*6.242e+15;
G_eV_list(7)=G2r*6.242e+15;
G_eV_list(8)=G3r*6.242e+15;
G_eV_list(9)=G4r*6.242e+15;
G_eV_list(10)=G5r*6.242e+15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% E0=pi*R*R*(0*(1-(T/Tm))-S*burger);
% E1=pi*R*R*(-35*(1-(T/Tm))-S*burger);
% E2=pi*R*R*(-183*(1-(T/Tm))-S*burger);
% E3=pi*R*R*(18*(1-(T/Tm))-S*burger);
% E4=pi*R*R*(-120*(1-(T/Tm))-S*burger);
% E5=pi*R*R*(-120*(1-(T/Tm))-S*burger);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=10;
% stress=S*10^-12;
b=a*sqrt(6)/6;
k=8.617*10^-5;
m=6.242*10^15;
l=20;
A=(1:1:3*l);
B=(1.333:1:3*l+0.333);
C=(1.666:1:3*l+0.666);
L(:,2)=0;
dE1=G1;
dE2=G2;
dE3=G3;
dE4=G4;
dE5=G5;
dE1r=G1r;
dE2r=G2r;
dE3r=G3r;
dE4r=G4r;
dE5r=G5r;
v1=(10^-9)*exp(-dE1*m/(T*k));
v2=(10^-9)*exp(-dE2*m/(T*k));
v3=(10^-9)*exp(-dE3*m/(T*k));
v4=(10^-9)*exp(-dE4*m/(T*k));
v5=(10^-9)*exp(-dE5*m/(T*k));
v1r=(10^-9)*exp(-dE1r*m/(T*k));
v2r=(10^-9)*exp(-dE2r*m/(T*k));
v3r=(10^-9)*exp(-dE3r*m/(T*k));
v4r=(10^-9)*exp(-dE4r*m/(T*k));
v5r=(10^-9)*exp(-dE5r*m/(T*k));
%%%%%%%%%%%%%%%%%%%%%%%%% Structure checking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:1:time;
    close;
for i=5:1:(3*N-5);
    if [L(i+2,1) L(i+1,1)+1 L(i,1)+2 L(i-1,1) L(i-2,1)+1 L(i-3,1)+2]==[3 3 3 3 3 3];
        L(i,3)=v1;
    end
    if [L(i+2,1)+2 L(i+1,1) L(i,1)+1 L(i-1,1)+2 L(i-2,1) L(i-3,1)+1]==[3 3 3 3 3 3];
        L(i,3)=v1;
    end
    if [L(i+2,1)+1 L(i+1,1)+2 L(i,1) L(i-1,1)+1 L(i-2,1)+2 L(i-3,1)]==[3 3 3 3 3 3];
        L(i,3)=v1;
    end
    if [L(i+2,1)+2 L(i+1,1) L(i,1)+1 L(i-1,1) L(i-2,1)+1 L(i-3,1)+2]==[3 3 3 3 3 3];
        L(i,3)=v1r;
    end
    if [L(i+2,1)+1 L(i+1,1)+2 L(i,1) L(i-1,1)+2 L(i-2,1) L(i-3,1)+1]==[3 3 3 3 3 3];
        L(i,3)=v1r;
    end
    if [L(i+2,1) L(i+1,1)+1 L(i,1)+2 L(i-1,1)+1 L(i-2,1)+2 L(i-3,1)]==[3 3 3 3 3 3];
        L(i,3)=v1r;
    end
    if [L(i+2,1) L(i+1,1)+1 L(i,1)+2 L(i-1,1) L(i-2,1)+1 L(i-3,1)]==[3 3 3 3 3 3];
        L(i,3)=v2;
    end
    if [L(i+2,1)+2 L(i+1,1) L(i,1)+1 L(i-1,1)+2 L(i-2,1) L(i-3,1)+2]==[3 3 3 3 3 3];
        L(i,3)=v2;
    end
    if [L(i+2,1)+1 L(i+1,1)+2 L(i,1) L(i-1,1)+1 L(i-2,1)+2 L(i-3,1)+1]==[3 3 3 3 3 3];
        L(i,3)=v2;
    end
    if [L(i+2,1)+2 L(i+1,1) L(i,1)+1 L(i-1,1) L(i-2,1)+1 L(i-3,1)]==[3 3 3 3 3 3];
        L(i,3)=v2r;
    end
    if [L(i+2,1)+1 L(i+1,1)+2 L(i,1) L(i-1,1)+2 L(i-2,1) L(i-3,1)+2]==[3 3 3 3 3 3];
        L(i,3)=v2r;
    end
    if [L(i+2,1) L(i+1,1)+1 L(i,1)+2 L(i-1,1)+1 L(i-2,1)+2 L(i-3,1)+1]==[3 3 3 3 3 3];
        L(i,3)=v2r;
    end
    if [L(i+2,1)+1 L(i+1,1)+2 L(i,1) L(i-1,1)+1 L(i-2,1) L(i-3,1)+1]==[3 3 3 3 3 3];
        L(i,3)=v3;
    end
    if [L(i+2,1) L(i+1,1)+1 L(i,1)+2 L(i-1,1) L(i-2,1)+2 L(i-3,1)]==[3 3 3 3 3 3];
        L(i,3)=v3;
    end
    if [L(i+2,1)+2 L(i+1,1) L(i,1)+1 L(i-1,1)+2 L(i-2,1)+1 L(i-3,1)+2]==[3 3 3 3 3 3];
        L(i,3)=v3;
    end
    if [L(i+2,1) L(i+1,1)+1 L(i,1)+2 L(i-1,1)+1 L(i-2,1) L(i-3,1)+1]==[3 3 3 3 3 3];
        L(i,3)=v3r;
    end
    if [L(i+2,1)+2 L(i+1,1) L(i,1)+1 L(i-1,1) L(i-2,1)+2 L(i-3,1)]==[3 3 3 3 3 3];
        L(i,3)=v3r;
    end
    if [L(i+2,1)+1 L(i+1,1)+2 L(i,1) L(i-1,1)+2 L(i-2,1)+1 L(i-3,1)+2]==[3 3 3 3 3 3]
        L(i,3)=v3r;
    end
    if [L(i+2,1) L(i+1,1)+1 L(i,1) L(i-1,1)+1 L(i-2,1) L(i-3,1)+1]==[3 3 3 3 3 3];
        L(i,3)=v4;
    end
    if [L(i+2,1)+2 L(i+1,1) L(i,1)+2 L(i-1,1) L(i-2,1)+2 L(i-3,1)]==[3 3 3 3 3 3];
        L(i,3)=v4;
    end
    if [L(i+2,1)+1 L(i+1,1)+2 L(i,1)+1 L(i-1,1)+2 L(i-2,1)+1 L(i-3,1)+2]==[3 3 3 3 3 3];
        L(i,3)=v4;
    end
    if [L(i+2,1)+2 L(i+1,1) L(i,1)+2 L(i-1,1)+1 L(i-2,1) L(i-3,1)+1]==[3 3 3 3 3 3];
        L(i,3)=v4r;
    end
    if [L(i+2,1)+1 L(i+1,1)+2 L(i,1)+1 L(i-1,1) L(i-2,1)+2 L(i-3,1)]==[3 3 3 3 3 3];
        L(i,3)=v4r;
    end
    if [L(i+2,1) L(i+1,1)+1 L(i,1) L(i-1,1)+2 L(i-2,1)+1 L(i-3,1)+2]==[3 3 3 3 3 3];
        L(i,3)=v4r;
    end
    if [L(i+2,1)+2 L(i+1,1) L(i,1)+1 L(i-1,1)+2 L(i-2,1)+1 L(i-3,1)]==[3 3 3 3 3 3];
        L(i,3)=v5;
    end
    if [L(i+2,1)+1 L(i+1,1)+2 L(i,1) L(i-1,1)+1 L(i-2,1) L(i-3,1)+2]==[3 3 3 3 3 3];
        L(i,3)=v5;
    end
    if [L(i+2,1) L(i+1,1)+1 L(i,1)+2 L(i-1,1) L(i-2,1)+2 L(i-3,1)+1]==[3 3 3 3 3 3];
        L(i,3)=v5;
    end
    if [L(i+2,1)+1 L(i+1,1)+2 L(i,1) L(i-1,1)+2 L(i-2,1)+1 L(i-3,1)]==[3 3 3 3 3 3];
        L(i,3)=v5r;
    end
    if [L(i+2,1) L(i+1,1)+1 L(i,1)+2 L(i-1,1)+1 L(i-2,1) L(i-3,1)+2]==[3 3 3 3 3 3];
        L(i,3)=v5r;
    end
    if [L(i+2,1)+2 L(i+1,1) L(i,1)+1 L(i-1,1) L(i-2,1)+2 L(i-3,1)+1]==[3 3 3 3 3 3];
        L(i,3)=v5r;
    end
end
p1=v1./sum(L(5:(3*N-5),3));
p1r=v1r./sum(L(5:(3*N-5),3));
p2=v2./sum(L(5:(3*N-5),3));
p2r=v2r./sum(L(5:(3*N-5),3));
p3=v3./sum(L(5:(3*N-5),3));
p3r=v3r./sum(L(5:(3*N-5),3));
p4=v4./sum(L(5:(3*N-5),3));
p4r=v4r./sum(L(5:(3*N-5),3));
p5=v5./sum(L(5:(3*N-5),3));
p5r=v5r./sum(L(5:(3*N-5),3));
L(5:(3*N-5),4)=L(5:(3*N-5),3)./sum(L(5:(3*N-5),3));
for i=5:1:(3*N-5);
    L(i,5)=L(i,4)+sum(L(1:(i-1),4));
end
seed=rand;
% t_seed=rand
for i=5:1:(3*N-5);
    if L(i,5)>seed && L(i-1,5)<seed;
        L(i,6)=1;
        if L(i-1,1)==L(i,1)+1;
            L(i:length(L),1)=L(i:length(L),1)-1;
            continue
        end
        if L(i-1,1)==L(i,1)-1;
            L(i:length(L),1)=L(i:length(L),1)+1;
            continue
        end
        if L(i-1,1)==L(i,1)-2;
            L(i:length(L),1)=L(i:length(L),1)-1;
            continue
        end
        if L(i-1,1)==L(i,1)+2;
            L(i:length(L),1)=L(i:length(L),1)+1;
            continue
        end
    end
end
for i=1:1:3*N;
            if L(i,1)<1;
                L(i,1)=L(i,1)+3;
            end
            if L(i,1)>3;
                L(i,1)=L(i,1)-3;
            end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=5:1:(3*N-5);
%     if [L(i+1,1) L(i,1) L(i-1,1)]==[1 2 3];
%         plot(i,B,'ok','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',8);
%         axis([0 3*N 0 3*l+1]);
%         hold on;
%     end
%     if [L(i+1,1) L(i,1) L(i-1,1)]==[2 3 1];
%         plot(i,C,'ok','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',8);
%         axis([0 3*N 0 3*l+1]);
%         hold on;
%     end
%     if [L(i+1,1) L(i,1) L(i-1,1)]==[3 1 2];
%         plot(i,A,'ok','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',8);
%         axis([0 3*N 0 3*l+1]);
%         hold on;
%     end
%     if [L(i+1,1) L(i,1) L(i-1,1)]==[3 2 3];
%         plot(i,B,'ok','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',8);
%         axis([0 3*N 0 3*l+1]);
%         hold on;
%     end
%     if [L(i+1,1) L(i,1) L(i-1,1)]==[1 3 1];
%         plot(i,C,'ok','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',8);
%         axis([0 3*N 0 3*l+1]);
%         hold on;
%     end
%     if [L(i+1,1) L(i,1) L(i-1,1)]==[2 1 2];
%         plot(i,A,'ok','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',8);
%         axis([0 3*N 0 3*l+1]);
%         hold on;
%     end
%     if [L(i+1,1) L(i,1) L(i-1,1)]==[2 3 2];
%         plot(i,C,'ok','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',8);
%         axis([0 3*N 0 3*l+1]);
%         hold on;
%     end
%     if [L(i+1,1) L(i,1) L(i-1,1)]==[3 1 3];
%         plot(i,A,'ok','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',8);
%         axis([0 3*N 0 3*l+1]);
%         hold on;
%     end
%     if [L(i+1,1) L(i,1) L(i-1,1)]==[1 2 1];
%         plot(i,B,'ok','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',8);
%         axis([0 3*N 0 3*l+1]);
%         hold on;
%     end
%     if [L(i+1,1) L(i,1) L(i-1,1)]==[3 2 1];
%         plot(i,B,'ok','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',8);
%         axis([0 3*N 0 3*l+1]);
%         hold on;
%     end
%     if [L(i+1,1) L(i,1) L(i-1,1)]==[2 1 3];
%         plot(i,A,'ok','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',8);
%         axis([0 3*N 0 3*l+1]);
%         hold on;
%     end
%     if [L(i+1,1) L(i,1) L(i-1,1)]==[1 3 2];
%         plot(i,C,'ok','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',8);
%         axis([0 3*N 0 3*l+1]);
%         hold on;
%     end
% end
% set(gcf,'unit','normalized','position',[0.0,0.3,1.0,0.5]);
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 8]);
% saveas(gcf,['KMC_300K.',num2str(t),'.png']);
end
for i=5:1:(3*N-5);
    if [L(i+1,1) L(i,1) L(i-1,1)]==[1 2 3];
        plot(i,B,'ok','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',8);
        axis([0 3*N 0 3*l+1]);
        hold on;
    end
    if [L(i+1,1) L(i,1) L(i-1,1)]==[2 3 1];
        plot(i,C,'ok','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',8);
        axis([0 3*N 0 3*l+1]);
        hold on;
    end
    if [L(i+1,1) L(i,1) L(i-1,1)]==[3 1 2];
        plot(i,A,'ok','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',8);
        axis([0 3*N 0 3*l+1]);
        hold on;
    end
    if [L(i+1,1) L(i,1) L(i-1,1)]==[3 2 3];
        plot(i,B,'ok','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',8);
        axis([0 3*N 0 3*l+1]);
        hold on;
    end
    if [L(i+1,1) L(i,1) L(i-1,1)]==[1 3 1];
        plot(i,C,'ok','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',8);
        axis([0 3*N 0 3*l+1]);
        hold on;
    end
    if [L(i+1,1) L(i,1) L(i-1,1)]==[2 1 2];
        plot(i,A,'ok','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',8);
        axis([0 3*N 0 3*l+1]);
        hold on;
    end
    if [L(i+1,1) L(i,1) L(i-1,1)]==[2 3 2];
        plot(i,C,'ok','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',8);
        axis([0 3*N 0 3*l+1]);
        hold on;
    end
    if [L(i+1,1) L(i,1) L(i-1,1)]==[3 1 3];
        plot(i,A,'ok','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',8);
        axis([0 3*N 0 3*l+1]);
        hold on;
    end
    if [L(i+1,1) L(i,1) L(i-1,1)]==[1 2 1];
        plot(i,B,'ok','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',8);
        axis([0 3*N 0 3*l+1]);
        hold on;
    end
    if [L(i+1,1) L(i,1) L(i-1,1)]==[3 2 1];
        plot(i,B,'ok','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',8);
        axis([0 3*N 0 3*l+1]);
        hold on;
    end
    if [L(i+1,1) L(i,1) L(i-1,1)]==[2 1 3];
        plot(i,A,'ok','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',8);
        axis([0 3*N 0 3*l+1]);
        hold on;
    end
    if [L(i+1,1) L(i,1) L(i-1,1)]==[1 3 2];
        plot(i,C,'ok','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',8);
        axis([0 3*N 0 3*l+1]);
        hold on;
    end
end
    set(gcf,'unit','normalized','position',[0.0,0.3,1.0,0.5]);
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 8]);
    saveas(gcf,['KMC_',num2str(T),'K.',num2str(t),'.png']);
end
