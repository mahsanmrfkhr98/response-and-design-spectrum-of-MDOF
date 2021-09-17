clc
clear all
close all
beta=1/6;
gama=1/2;
delta_t=0.02;
Tn=0.001:0.1:15;
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
       f0(i)=(Wn(i)^2)*umax2darsad(i)
   
end 
figure(1)
plot(Tn,umax2darsad*981,'b');grid on;hold on
legend('umax2darsad')
xlabel('period(sec)','FontSize',12);
ylabel(' max Displacement','FontSize',12);

figure(2)
plot(Tn,udotmax2darsad*981,'b');grid on;hold on
legend('udotmax2darsad')
xlabel('period(sec)','FontSize',12);
ylabel(' max velocity','FontSize',12);
figure(3)
plot(Tn,uddotmax2darsad,'b');grid on;hold on
legend('uddotmax2darsad')
xlabel('period(sec)','FontSize',12);
ylabel(' max Acceleration','FontSize',12);
figure(4)
plot(umax2darsad*981,f0*981)
%% miu>1
R=1:0.1:10;
t_total=2;
delta_t=0.02;

a1=1/(beta*delta_t^2)+(c*gama/(beta*delta_t));
a2=(1/(beta*delta_t))+gama*c/beta;
a3=(gama*c*delta_t/(2*beta))+(0.5*1/beta)-c*delta_t;
n=2/1
 R=zeros(length(Tn),n);
 ubi=zeros(length(Tn),length(t),n);
 ubidot=zeros(length(Tn),length(t),n);
 ubiddot=zeros(length(Tn),length(t),n);
 fs=zeros(length(Tn),length(t),n);
 delta_U=zeros(length(Tn),length(t)+2,n);
 delta_Udot=zeros(length(Tn),length(t)+2,n);
 delta_Uddot=zeros(length(Tn),length(t)+2,n);
 

 for j=1:1:length(t) 
    for l=1:1:n
         for i=1:150;
             Pl(i,j,l)=P(1,j)*9.81;
         end
    end
 end
 


for l=1:1:n 
         for i=1:1:length(Tn)
              ubiddot(i,1,l)=Pl(1,1,1);
             fy(i,l)=(f0(i)*9.81)/R(1,l);
              
             if i==1
                 fy(i,l)=f0(i+1)*9.81/R(1,l);
                 k(i,j,l)=(Wn(i)^2)+a1(1,i);
                 uy(i,l)=fy(i,l)/k(i,j,l);
             end
             fbary(i,l)=1/R(1,l);
             k0=(Wn(1,i)^2)+a1(1,i);
             for j=1:1:(length(t)-1);
                
                 if abs(fs(i,j,l))<=fy(i,l);
                    k(i,j,l)=(Wn(i)^2)+a1(1,i);
                    
                 else
                    k(i,j,l)=0.0001 ;
                    uy(i,l)=fy(i,l)/k(i,j-1,l);
                 end
                    k_hat(i,j,l)=k(i,j,l)+a1(1,i);
                    delta_U(i,j,l)=(Pl(i,j+1,l)-Pl(i,j,l)+a2(1,i).*ubidot(i,j,l)+a3(1,i).*ubiddot(i,j,l))/k_hat(i,j,l);
                    delta_Udot(i,j,l)=(1-0.5*gama/beta)*delta_t*ubiddot(i,j,l)+(gama/(beta*delta_t))*delta_U(i,j,l)-gama*ubidot(i,j,l)/beta;
                    delta_Uddot(i,j,l)=(1/(beta*(delta_t^2)))*delta_U(i,j,l)-(1/(beta*delta_t))*ubidot(i,j,l)-0.5*ubiddot(i,j,l)/beta;
                    ubi(i,j+1,l)=ubi(i,j,l)+delta_U(i,j,l);
                    ubidot(i,j+1,l)=delta_Udot(i,j,l)+ubidot(i,j,l);
                    fs(i,j+1,l)=fs(i,j,l)+k(i,j,l)*delta_U(i,j,l);
                    
                    if fs(i,j+1,l)>fy(i,l)
                    fs(i,j+1,l)=fy(i,l);
                    elseif fs(i,j+1,l)<-fy(i,l)
                    fs(i,j+1,l)=-fy(i,l);
                    end
                    ubiddot(i,j+1,l)=(Pl(i,j+1,l)-c(i)*ubidot(i,j+1,l)-fs(i,j+1,l));
                    %s=uy(i,l);
                   %if ubi(i,j,l)>s
                      % s=ubi(i,j,l);
                   %end
             end
             %um(i,l)=s;
            % miu(i,l)= abs(um(i,l))/uy(i,l);
             %if miu(i,l)==2
              %   um2(l)=um(i,l);
              %   uy2(l)=uy(i,l);
            % elseif miu==8
               %  um8(l)=um(i,l);
              %   uy8(l)=uy(i,l);
             %end  
      end
end

figure(5)
%plot(Tn,um(:,:,2),'--');hold on
%plot(Tn,uy(:,:,2),"b")
