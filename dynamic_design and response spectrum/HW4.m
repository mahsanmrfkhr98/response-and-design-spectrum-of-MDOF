clc
clear all
close all
beta=1/6;
gama=1/2;
delta_t=0.02;
Tn=0.001:0.01:5;
Wn=zeros(1,length(Tn));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kisay=0.02;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elcentro=xlsread('elcentro.xlsx');
t=elcentro(:,1);
p=-elcentro(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u=zeros(length(Tn),length(t));
udot=zeros(length(Tn),length(t));
uddot=zeros(length(Tn),length(t));
p_hat=zeros(length(Tn),length(t));
uddot_total=zeros(length(Tn),length(t));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  miu=1

for i=1:1:length(Tn)
    q=u(1,1);
    b=udot(1,2);
    a=uddot_total(1,2);
    Wn(i)=(2*pi)/Tn(i);
    c(i)=2*kisay*Wn(i);
    if i==1
        Wn=0
        c=0 
    end
    a1=(1/(beta*(delta_t^2)))+(gama/(beta*delta_t))*c(i);
    a2=(1/(beta*delta_t))+((gama/(beta))-1)*c(i);
    a3=((1/(2*beta))-1)*1+delta_t*((gama/(2*beta))-1)*c(i);
    k_hat(i)=(Wn(i)^2)+a1;
      for j=1:1:(length(t)-1)
           P=p';
           F(i,j+1)=P(1,j+1);
           p_hat(1,j+1)=P(1,j+1)+a1.*u(i,j)+a2.*udot(i,j)+a3.*uddot(i,j);
           u(i,j+1)=p_hat(1,j+1)/k_hat(i);
           udot(i,j+1)=(gama/(beta*delta_t)).*(u(i,j+1)-u(i,j))+(1-(gama/beta)).*udot(i,j)+delta_t.*(1-(gama/(2*beta))).*uddot(i,j);
           uddot(i,j+1)=(1/(beta*(delta_t^2))).*(u(i,j+1)-u(i,j))-(1/(beta*delta_t)).*udot(i,j)-(((1/(2*beta))-1).*uddot(i,j));
           uddot_total(i,j+1)=uddot(i,j+1)-F(i,j+1);
           
            
          if u(i,j+1)>q
              q=u(i,j+1);
          end
                if udot(i,j+1)>b
                   b=udot(i,j+1);
                end
           if uddot_total(i,j+1)>a
             a=uddot_total(i,j+1);
           end
      end
      
       umax2darsad(i)=q;
       udotmax2darsad(i)=b;
       uddotmax2darsad(i)=a;
       umax2darsad(1)=0;
       udotmax2darsad(1)=0;
       uddotmax2darsad(1)=0;
       f0(i)=(Wn(i)^2)*umax2darsad(i);
   
end 
%figure(1)
%plot(Tn,umax2darsad*981,'b');grid on;hold on
%legend('umax2darsad')
%xlabel('period(sec)','FontSize',12);
%ylabel(' max Displacement','FontSize',12);

%figure(2)
%plot(Tn,udotmax2darsad*981,'b');grid on;hold on
%legend('udotmax2darsad')
%xlabel('period(sec)','FontSize',12);
%ylabel(' max velocity','FontSize',12);
%figure(3)
%plot(Tn,uddotmax2darsad,'b');grid on;hold on
%legend('uddotmax2darsad')
%xlabel('period(sec)','FontSize',12);
%ylabel(' max Acceleration','FontSize',12);


%% mu>1
%% step1
m=1;
mu22=2;
mu88=8;
Tn=0.001:0.01:5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u=zeros(length(t),length(Tn));
ud=zeros(length(t),length(Tn));
udd=zeros(length(t),length(Tn));
fs=zeros(length(t),length(Tn));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for j=1:1:length(Tn)
     Wn(j)=(2*pi)/Tn(j);
a1=(m/(beta*delta_t))+gama*2*m*Wn(j)*kisay/beta;
a2=(gama*2*Wn(j)*m*kisay*delta_t/(2*beta))+(0.5*m/beta)-2*Wn(j)*m*kisay*delta_t;
    for l=1:10000
        fbary(l,j)=0.0001*(l);
        fy(l,j)=fbary(l,j)*f0(1,j);
         for i=1:1:(length(t)-1)
              if abs(fs(i,j))<fy(l,j)
                k(i,j)=m*(Wn(j)^2);
             else
               k(i,j)=0.8 ;
             end
             khat(i,j)=k(i,j)+m/(beta*delta_t^2)+(2*1*m*Wn(j)*kisay*gama/(beta*delta_t));
             deltaU(i,j)=(p(i+1)-p(i)+a1*ud(i,j)+a2*udd(i,j))/khat(i,j);
             deltaUd(i,j)=(1-0.5*gama/beta)*delta_t*udd(i,j)+(gama/(beta*delta_t))*deltaU(i,j)-gama*ud(i,j)/beta;
             deltaUdd(i,j)=(1/(beta*(delta_t^2)))*deltaU(i,j)-(1/(beta*delta_t))*ud(i,j)-0.5*udd(i,j)/beta;
             u(i+1,j)=u(i,j)+deltaU(i,j);
             ud(i+1,j)=deltaUd(i,j)+ud(i,j);
             fs(i+1,j)=fs(i,j)+k(i,j)*deltaU(i,j);
             if fs(i+1,j)>fy(l,j)
                 fs(i+1,j)=fy(l,j);
             end
             if fs(i+1,j)<-fy(l,j)
                fs(i+1,j)=-fy(l,j);
             end
             udd(i+1,j)=(p(i+1)-2*m*Wn(j)*kisay*ud(i+1,j)-fs(i+1,j))/m;
         end
            mu2(l,j)=(max(abs(u(:,j))/umax2darsad(j)))/fbary(l,j);
      end
end
fbary2=zeros(1,length(Tn));
fbary8=zeros(1,length(Tn));
for j=1:1:length(Tn)
    for l=1:1:9999
        if mu2(l,j)>mu22 && mu2(l+1,j)<mu22
            fbary2(j)=(0.0001/(mu2(l,j)-mu2(l+1,j)))*(mu22-mu2(l+1,j))+0.0001*l;
        end
         if mu2(l,j)>mu88 && mu2(l+1,j)<mu88
            fbary8(j)=(0.0001/(mu2(l,j)-mu2(l+1,j)))*(mu88-mu2(l+1,j))+0.0001*l;
        end
    end
end
  %  for j=2:1:length(Tn)
  %      figure(j-1);
 % plot(mu2(:,j),1./fbary(:,j),'b','LineWidth',2);hold on;
%  plot([mu22,mu22],[0,20],'r','LineWidth',2);hold on;
%  plot([mu88,mu88],[0,20],'g','LineWidth',2);
%  xlabel('\mu','FontSize',20);
 % ylabel(' R','FontSize',20);
%  title('\mu-R relation','FontSize',20);
%  legend('\mu-R relation','\mu=2','\mu=8');
%     end
    figure(length(Tn));
     plot(Tn(1,2:length(Tn)),fbary2(1,2:length(Tn)),'b','LineWidth',2);hold on;
     plot(Tn(1,2:length(Tn)),fbary8(1,2:length(Tn)),'r','LineWidth',2);
     xlabel('T_n','FontSize',20);
     ylabel(' normalized strength(1/R)','FontSize',20);
     title(' normalized strength for each period','FontSize',20);
     legend('\mu=2','\mu=8');
     figure(length(Tn)+1);
     plot(Tn,1./fbary2,'b','LineWidth',2);hold on;
     plot(Tn,1./fbary8,'r','LineWidth',2);
     xlabel('T_n','FontSize',20);
     ylabel(' R','FontSize',20);
     title(' reduction factor for each period','FontSize',20);
     legend('\mu=2','\mu=8');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5 
    %% step2
    Tn=0.001:0.01:5;
u2=zeros(length(t),length(Tn));
ug=zeros(length(t),length(Tn));
ud2=zeros(length(t),length(Tn));
udd2=zeros(length(t),length(Tn));
deltaU2=zeros(length(t)-1,length(Tn));
deltaUd2=zeros(length(t)-1,length(Tn));
deltaUdd2=zeros(length(t)-1,length(Tn));
fs2=zeros(length(t),length(Tn));
k2=zeros(length(t),length(Tn));
khat2=zeros(length(t),length(Tn));
fy2=zeros(1,length(Tn));
u8=zeros(length(t),length(Tn));
ug8=zeros(length(t),length(Tn));
ud8=zeros(length(t),length(Tn));
udd8=zeros(length(t),length(Tn));
deltaU8=zeros(length(t)-1,length(Tn));
deltaUd8=zeros(length(t)-1,length(Tn));
deltaUdd8=zeros(length(t)-1,length(Tn));
fs8=zeros(length(t),length(Tn));
k8=zeros(length(t),length(Tn));
khat8=zeros(length(t),length(Tn));
fy8=zeros(1,length(Tn));
U2=zeros(1,length(Tn));
Ud2=zeros(1,length(Tn));
Udd2=zeros(1,length(Tn));
U8=zeros(1,length(Tn));
Ud8=zeros(1,length(Tn));
Udd8=zeros(1,length(Tn));
U2s=zeros(1,length(Tn));
Ud2s=zeros(1,length(Tn));
Udd2s=zeros(1,length(Tn));
U8s=zeros(1,length(Tn));
Ud8s=zeros(1,length(Tn));
Udd8s=zeros(1,length(Tn));
uddg=zeros(length(t),length(Tn));
 
for j=1:1:length(Tn)
    Wn(j)=(2*pi)/Tn(j);
   a1=(m/(beta*delta_t))+gama*2*m*Wn(j)*kisay/beta;
   a2=(gama*2*Wn(j)*m*kisay*delta_t/(2*beta))+(0.5*m/beta)-2*Wn(j)*m*kisay*delta_t;
   fy2(j)=fbary2(j)*umax2darsad(j)*m*(Wn(j)^2);
   fy8(j)=fbary8(j)*umax2darsad(j)*m*(Wn(j)^2);
for i=1:1:(length(t)-1)
     if abs(fs2(i,j))<fy2(j)
        k2(i,j)=m*(Wn(j)^2);
    else
        k2(i,j)=0.08 ;
     end  
    if abs(fs8(i,j))<fy8(j)
        k8(i,j)=m*(Wn(j)^2);
    else
        k8(i,j)=0.08 ;
    end
    khat2(i,j)=k2(i,j)+m/(beta*delta_t^2)+(2*1*m*Wn(j)*kisay*gama/(beta*delta_t));
    deltaU2(i,j)=(p(i+1)-p(i)+a1*ud2(i,j)+a2*udd2(i,j))/khat2(i,j);
    deltaUd2(i,j)=(1-0.5*gama/beta)*delta_t*udd2(i,j)+(gama/(beta*delta_t))*deltaU2(i,j)-gama*ud2(i,j)/beta;
    deltaUdd2(i,j)=(1/(beta*(delta_t^2)))*deltaU2(i,j)-(1/(beta*delta_t))*ud2(i,j)-0.5*udd2(i,j)/beta;
    u2(i+1,j)=u2(i,j)+deltaU2(i,j);
    ud2(i+1,j)=deltaUd2(i,j)+ud2(i,j);
    fs2(i+1,j)=fs2(i,j)+k2(i,j)*deltaU2(i,j);
    khat8(i,j)=k8(i,j)+m/(beta*delta_t^2)+(2*1*m*Wn(j)*kisay*gama/(beta*delta_t));
    deltaU8(i,j)=(p(i+1)-p(i)+a1*ud8(i,j)+a2*udd8(i,j))/khat8(i,j);
    deltaUd8(i,j)=(1-0.5*gama/beta)*delta_t*udd8(i,j)+(gama/(beta*delta_t))*deltaU8(i,j)-gama*ud8(i,j)/beta;
    deltaUdd8(i,j)=(1/(beta*(delta_t^2)))*deltaU8(i,j)-(1/(beta*delta_t))*ud8(i,j)-0.5*udd8(i,j)/beta;
    u8(i+1,j)=u8(i,j)+deltaU8(i,j);
    ud8(i+1,j)=deltaUd8(i,j)+ud8(i,j);
    fs8(i+1,j)=fs8(i,j)+k8(i,j)*deltaU8(i,j);    
   if fs2(i+1,j)>fy2(j)
       fs2(i+1,j)=fy2(j);
   end
   if fs2(i+1,j)<-fy2(j)
       fs2(i+1,j)=-fy2(j);
   end
   if fs8(i+1,j)>fy8(j)
       fs8(i+1,j)=fy8(j);
   end
   if fs8(i+1,j)<-fy8(j)
       fs8(i+1,j)=-fy8(j);
   end    
       udd2(i+1,j)=(p(i+1)-2*m*Wn(j)*kisay*ud2(i+1,j)-fs2(i+1,j))/m;
       udd8(i+1,j)=(p(i+1)-2*m*Wn(j)*kisay*ud8(i+1,j)-fs8(i+1,j))/m;
       uddg(i+1,j)=-p(i+1)/m;
end
end
Uddt2=uddg+udd2;
Uddt8=uddg+udd8;
for j=1:1:length(Tn)
U2(j)=max(abs(u2(:,j)));
Ud2(j)=max(abs(ud2(:,j)));
Udd2(j)=max(abs(Uddt2(:,j)));
U8(j)=max(abs(u8(:,j)));
Ud8(j)=max(abs(ud8(:,j)));
Udd8(j)=max(abs(Uddt8(:,j)));
U2s(j)=U2(j)/mu22;
Ud2s(j)=Wn(j)*U2s(j);
Udd2s(j)=(Wn(j)^2)*U2s(j);
U8s(j)=U8(j)/mu88;
Ud8s(j)=Wn(j)*U8s(j);
Udd8s(j)=(Wn(j)^2)*U8s(j);
end

bilinear=xlsread('bilinear.xlsx');
Tnseismo=0:0.01:4.99
U2seismo=bilinear(:,3);
Ud2seismo=bilinear(:,8);
Udd2seismo=bilinear(:,13);
U8seismo=bilinear(:,4);
Ud8seismo=bilinear(:,9);
Udd8seismo=bilinear(:,14);
figure(length(Tn)+2);
     plot(Tn,U2*981,'b','LineWidth',2);hold on;
     plot(Tnseismo,U2seismo,'g--','LineWidth',2);hold on;
     plot(Tn,U8*981,'r','LineWidth',2);hold on;
     plot(Tnseismo,U8seismo,'k--','LineWidth',2);
     title(' response spectrum for displasement','FontSize',20);
     xlabel('Tn(sec)','FontSize',20);
     ylabel('displacement(m)','FontSize',20);
     legend('Numerical,\mu=2','SeismoSignal,\mu=2','Numerical,\mu=8','SeismoSignal,\mu=8'); 
figure(length(Tn)+3);
     plot(Tn,Ud2*981,'b','LineWidth',2);hold on;
     plot(Tnseismo,Ud2seismo,'b--','LineWidth',2);hold on;
     plot(Tn,Ud8*981,'r','LineWidth',2);hold on;
     plot(Tnseismo,Ud8seismo,'r--','LineWidth',2);
     title(' response spectrum for velocity','FontSize',20);
     xlabel('Tn(sec)','FontSize',20);
     ylabel('velocity(m/sec)','FontSize',20);
     legend('Numerical,\mu=2','SeismoSignal,\mu=2','Numerical,\mu=8','SeismoSignal,\mu=8');
figure(length(Tn)+4);
     plot(Tn,Udd2,'b','LineWidth',2);hold on;
     plot(Tnseismo,Udd2seismo,'b--','LineWidth',2);hold on;
     plot(Tn,Udd8,'r','LineWidth',2);hold on;
     plot(Tnseismo,Udd8seismo,'r--','LineWidth',2);
     title(' response spectrum for acceleration','FontSize',20);
     xlabel('Tn(sec)','FontSize',20);
     ylabel('Total acceleration(g)','FontSize',20);
     legend('Numerical,\mu=2','SeismoSignal,\mu=2','Numerical,\mu=8','SeismoSignal,\mu=8');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step3
ug0=0.281;
udg0=0.363485;
uddg0=0.31882;
alphad=2.42;
alphav=2.92;
alphaa=3.66;
g=9.81;
Tn2=[0.02,1/33,1/8,udg0*2*pi*alphav*((2*mu22-1)^0.5)/(mu22*alphaa*uddg0*g),2*pi*ug0*alphad/(udg0*alphav),10,33,50];
Adesign2=[uddg0,uddg0,uddg0*alphaa/((2*mu22-1)^0.5),uddg0*alphaa/((2*mu22-1)^0.5),alphav*udg0*2*pi/(g*mu22*Tn2(5)),((2*pi/Tn2(6))^2)*alphad*ug0/(g*mu22),((2*pi/Tn2(7))^2)*ug0/(g*mu22),((2*pi/Tn2(8))^2)*ug0/(g*mu22)];
Vdesign2=(1/(2*pi)).*Tn2.*Adesign2;
Ddesign2=((1/(2*pi))^2).*(Tn2.^2).*Adesign2;
Tn8=[0.02,1/33,1/8,udg0*2*pi*alphav*((2*mu88-1)^0.5)/(mu88*alphaa*uddg0*g),2*pi*ug0*alphad/(udg0*alphav),10,33,50];
Adesign8=[uddg0,uddg0,uddg0*alphaa/((2*mu88-1)^0.5),uddg0*alphaa/((2*mu88-1)^0.5),alphav*udg0*2*pi/(g*mu88*Tn2(5)),((2*pi/Tn2(6))^2)*alphad*ug0/(g*mu88),((2*pi/Tn2(7))^2)*ug0/(g*mu88),((2*pi/Tn2(8))^2)*ug0/(g*mu88)];
Vdesign8=(1/(2*pi)).*Tn8.*Adesign8;
Ddesign8=((1/(2*pi))^2).*(Tn8.^2).*Adesign8;
figure(length(Tn)+5);
loglog(Tn2,Adesign2,'b','LineWidth',2);hold on;
loglog(Tn8,Adesign8,'r','LineWidth',2);
title(' Inelastic (pseudo-acceleration) design spectrum','FontSize',20);
xlabel('Tn(sec)','FontSize',20);
ylabel('A_y(g)','FontSize',20);
legend('\mu=2','\mu=8');
figure(length(Tn)+6);
loglog(Tn2,Vdesign2,'b','LineWidth',2);hold on;
loglog(Tn8,Vdesign8,'r','LineWidth',2);
title(' Inelastic (pseudo-velocity) design spectrum','FontSize',20);
xlabel('Tn(sec)','FontSize',20);
ylabel('V_y(g)','FontSize',20);
legend('\mu=2','\mu=8');
figure(length(Tn)+7);
loglog(Tn2,Ddesign2,'b','LineWidth',2);hold on;
loglog(Tn8,Ddesign8,'r','LineWidth',2);
title(' Inelastic (yield displacement) design spectrum','FontSize',20);
xlabel('Tn(sec)','FontSize',20);
ylabel('D_y(g)','FontSize',20);
legend('\mu=2','\mu=8');
Ay2=zeros(1,length(Tn));
Ay8=zeros(1,length(Tn));
Vy2=zeros(1,length(Tn));
Vy8=zeros(1,length(Tn));
Dy2=zeros(1,length(Tn));
Dy8=zeros(1,length(Tn));

for i=1:1:length(Tn)
    if Tn(i)<Tn2(2)
        Ay2(i)=uddg0;
    elseif Tn(i)>Tn2(2) && Tn(i)<Tn2(3)
        Ay2(i)=1.42*(Tn(i)^0.54);
    elseif Tn(i)>Tn2(3) && Tn(i)<Tn2(4) 
        Ay2(i)=0.672;
    elseif Tn(i)>Tn2(4) 
        Ay2(i)=0.355*(Tn(i)^(-1.005));
    end
    if Tn(i)<Tn8(2)
        Ay8(i)=uddg0;
    elseif Tn(i)>Tn8(2) && Tn(i)<Tn8(3)
          Ay8(i)=0.279*Tn(i)^(-0.035);
    elseif Tn(i)>Tn8(3) && Tn(i)<Tn8(4)
          Ay8(i)=0.3;
    elseif Tn(i)>Tn8(4) 
          Ay8(i)=0.0886*Tn(i)^(-1.00215);
    end
    Vy2(i)=Ay2(i)/Wn(i);
    Dy2(i)=Ay2(i)/(Wn(i)^2);
    Vy8(i)=Ay8(i)/Wn(i);
    Dy8(i)=Ay8(i)/(Wn(i)^2); 
   
end
design=xlsread('design.xlsx');
Tnseismo=0:0.01:4.99
Udd2s=design(6:505,5)
Udd8s=design(6:505,8)
figure(length(Tn)+8);
     plot(Tn,Ay2,'b','LineWidth',2);hold on;
     plot(Tn,Udd2s,'g--','LineWidth',2);hold on;
     plot(Tn,Ay8,'r','LineWidth',2);hold on;
    plot(Tn,Udd8s,'k--','LineWidth',2);
     title('Design & response spectrum for pseudo acceleration','FontSize',20);
     xlabel('Tn(sec)','FontSize',20);
     ylabel('pseudo acceleration(g)','FontSize',20);
     legend('Design spectrum,\mu=2','response spectrum,\mu=2','Design spectrum,\mu=8','response spectrum,\mu=8');
% figure(length(Tn)+7);
%      plot(Tn,Vy2,'b','LineWidth',2);hold on;
%      plot(Tn,Ud2s,'b--','LineWidth',2);hold on;
%      plot(Tn,Vy8,'r','LineWidth',2);hold on;
%      plot(Tn,Ud8s,'r--','LineWidth',2);
%      title('Design & response spectrum for pseudo velocity','FontSize',20);
%      xlabel('Tn(sec)','FontSize',20);
%      ylabel('pseudo velocity(m/sec)','FontSize',20);
%      legend('Design spectrum,\mu=2','response spectrum,\mu=2','Design spectrum,\mu=8','response spectrum,\mu=8');    
% figure(length(Tn)+8);
%      plot(Tn,Dy2,'b','LineWidth',2);hold on;
%      plot(Tn,U2s,'b--','LineWidth',2);hold on;
%      plot(Tn,Dy8,'r','LineWidth',2);hold on;
%      plot(Tn,U8s,'r--','LineWidth',2);
     
     
