
%in the name of god
%project1-simulation
%Dynamicmachine
%Sajjad Sarrafzadeh
%ST code: 810698267


clear 
clc
close all
tic
syms re4 te4 rd4 thd4 rdd4 thdd4 re9 te11 rd9 rdd9 thd11 thdd11
%analyzing motion in n frames
n=200;
% t=linspace(0,2*pi,n);
time=linspace(0,20,n);
O2B=250;
O4C=800;
CD=350;
O2O4=500;
%define vectors and magnitudes
r10=O4C*ones(1,n);
r9=zeros(1,n);
theta9=zeros(1,n);
theta8=3*pi/2*ones(1,n);
theta11=zeros(1,n);
r2=O2B*ones(1,n);
r6=O2O4*ones(1,n);
r11=CD*ones(1,n);
r8=900*ones(1,n);
theta2=linspace(deg2rad(225),20+deg2rad(225),n);
theta6=3*pi/2*ones(1,n);
r4=zeros(1,n);
theta4=zeros(1,n);

thd2=ones(1,n);
rdot4=zeros(1,n);
thetadot4=zeros(1,n);
thetaddot4=zeros(1,n);
rddot4=zeros(1,n);
theta10=zeros(1,n);
thetadot10=zeros(1,n);
thetaddot10=zeros(1,n);
rdot9=zeros(1,n);
rddot9=zeros(1,n);
thetadot11=zeros(1,n);
thetaddot11=zeros(1,n);
%forces
f12y=zeros(1,n);
f12x=zeros(1,n);
f32x=zeros(1,n);
f32y=zeros(1,n);
f14x=zeros(1,n);
f14y=zeros(1,n);
f54x=zeros(1,n);
f54y=zeros(1,n);
f34x=zeros(1,n);
f34y=zeros(1,n);
f65x=zeros(1,n);
f65y=zeros(1,n);
f16x=zeros(1,n);
f16y=zeros(1,n);
t2=zeros(1,n);
t3=zeros(1,n);
t6=zeros(1,n);









for i=1:n
    
    
    eq1=r2(1,i).*cos(theta2(1,i))+re4.*cos(te4)+r6(1,i).*cos(theta6(1,i));
    eq2=r2(1,i).*sin(theta2(1,i))+re4.*sin(te4)+r6(1,i).*sin(theta6(1,i));
    rr4=zeros(1,11);
    thett4=zeros(1,11);
    rr9=zeros(1,11);
    thett11=zeros(1,11);
    if i==1
        rr9(1)=124.3653;
        thett11(1)=2.7734;
    else
        rr9(1)=r9(i-1);
        thett11(1)=theta11(i-1);
    end
    if i==1   
        rr4=699.4832;
        thett4=1.3153;
    else
        rr4=r4(i-1);
        thett4=theta4(i-1);
    end
    j1=[diff(eq1,re4),diff(eq1,te4);diff(eq2,re4),diff(eq2,te4)];
    f1=[eq1;eq2];
    v1=j1\f1;
    

    %newton-raphson code for solving r4 and theta4 in each frame (loop1)
    m=1;
    while m<4
        
        A=[rr4;thett4]-subs(v1,[re4,te4],[rr4,thett4]);
        rr4=A(1);
        thett4=A(2);
        
        

       
        
        m=m+1;
        
        
        
        
        
        
        
        
    end
    r4(i)=rr4;
    theta4(i)=thett4;

   %analyzing velocity and calculating rdot4 and thetadot4 in each frame(loop1)
   eq3=rd4.*cos(theta4(i))-r2(i).*thd2(i).*sin(theta2(i))-r4(i).*thd4.*sin(theta4(i))==0;
   eq4=rd4.*sin(theta4(i))+r2(i).*thd2(i).*cos(theta2(i))+r4(i).*thd4.*cos(theta4(i))==0;
   s2=vpasolve([eq3,eq4],[rd4,thd4]);
   thetadot4(i)=s2.thd4;
   rdot4(i)=s2.rd4;
   %analyzing acceleration and calculating rdoubledot4 and thetadoubledot4 in each frame(loop1)
   eq5=-r2(i).*thd2(i)^2.*cos(theta2(i))+rdd4.*cos(theta4(i))-2*thetadot4(i).*rdot4(i).*sin(theta4(i))-r4(i)*thdd4.*sin(theta4(i))-r4(i).*thetadot4(i)^2.*cos(theta4(i))==0;
   eq6=-r2(i).*thd2(i)^2.*sin(theta2(i))+rdd4.*sin(theta4(i))+2*thetadot4(i).*rdot4(i).*cos(theta4(i))+r4(i)*thdd4.*cos(theta4(i))-r4(i).*thetadot4(i)^2.*sin(theta4(i))==0;
   s3=vpasolve([eq5,eq6],[rdd4,thdd4]);
   thetaddot4(i)=s3.thdd4;
   rddot4(i)=s3.rdd4;
   theta10(i)=theta4(i);
   thetadot10(i)=thetadot4(i);
   thetaddot10(i)=thetaddot4(i);
   eq7=r8(1,i).*cos(theta8(1,i))+re9.*cos(theta9(1,i))+r10(1,i).*cos(theta10(1,i))+r11(1,i).*cos(te11);
   eq8=r8(1,i).*sin(theta8(1,i))+re9.*sin(theta9(1,i))+r10(1,i).*sin(theta10(1,i))+r11(1,i).*sin(te11);
   if i==1
        rr9=124.3653;
        thett11=2.7734;
    else
        rr9=r9(i-1);
        thett11=theta11(i-1);
   end
  
    j2=[diff(eq7,re9),diff(eq7,te11);diff(eq8,re9),diff(eq8,te11)];
    f2=[eq7;eq8];
    v2=j2\f2;
    k=1;
    while k<4
        
        B=[rr9;thett11]-subs(v2,[re9,te11],[rr9,thett11]);
        rr9=B(1);
        thett11=B(2);
        k=k+1;
    
    
    
    end
    r9(i)=rr9;
    theta11(i)=thett11;
    eq9=rd9.*cos(theta9(i))-r10(i).*thetadot10(i).*sin(theta10(i))-r11(i).*thd11.*sin(theta11(i))==0;
    eq10=rd9.*sin(theta9(i))+r10(i).*thetadot10(i).*cos(theta10(i))+r11(i).*thd11.*cos(theta11(i))==0;
    s4=vpasolve([eq9,eq10],[rd9,thd11]);
    rdot9(i)=s4.rd9;
    thetadot11(i)=s4.thd11;
    eq11=rdd9.*cos(theta9(i))-thetaddot10(i).*r10(i).*sin(theta10(i))-thetadot10(i)^2.*r10(i).*cos(theta10(i))-thdd11.*r11(i).*sin(theta11(i))-thetadot11(i)^2.*r11(i).*cos(theta11(i))==0;
    eq12=rdd9.*sin(theta9(i))+thetaddot10(i).*r10(i).*cos(theta10(i))-thetadot10(i)^2.*r10(i).*sin(theta10(i))+thdd11.*r11(i).*cos(theta11(i))-thetadot11(i)^2.*r11(i).*sin(theta11(i))==0;
    s5=vpasolve([eq11,eq12],[rdd9,thdd11]);
    rddot9(i)=s5.rdd9;
    thetaddot11(i)=s5.thdd11;
    
end


x2=r2.*cos(theta2);
y2=r2.*sin(theta2);
x4=r4.*cos(theta4);
y4=r4.*sin(theta4)-500;
y10=r10.*sin(theta4)-500;
x10=r10.*cos(theta4);
x11=r11.*cos(theta11);
y11=r11.*sin(theta11);

for l=1:n
    
    
   X2=[0,-x2(l)];
   Y2=[0,-y2(l)];
   X4=[0,x4(l)];
   Y4=[-500,y4(l)];
   X10=[0,x10(l)];
   Y10=[-500,y10(l)];
   Y11=[y10(l),400];
   X11=[x10(l),-r9(l)];
   plot(X2,Y2,'r','linewidth',10)
   hold on
 
   plot(X10,Y10,'g','linewidth',10)
   plot(X4,Y4,'blue','linewidth',10)
   plot(X11,Y11,'linewidth',10)
   hold off
 
    
    xlim([-700 500])
    ylim([-500,500])

    
    grid on
    
    
    
     pause(0.1)
    
    
end

G2O2x=abs(O2B/2.*cos(theta2));
G2Bx=G2O2x;
G2O2y=O2B/2.*sin(theta2);
G2By=G2O2y;
G4O4x=abs(O4C/2.*cos(theta4));
G4Cx=G4O4x;
G4O4y=abs(O4C/2.*sin(theta4));
G4Cy=G4O4y;
G5Dx=abs(CD/2.*cos(theta11));
G5Cx=G5Dx;
G5Dy=abs(CD/2.*sin(theta11));
G5Cy=G5Dy;
m4=O4C/1000;
m5=CD/1000;
m2=O2B/1000;
m3=0.3;
m6=0.3;
mo=0.25;
syms f12xx f12yy f32xx f32yy f14xx f14yy f54xx f54yy f65xx f65yy T2 f34 f34xx f34yy T3 fs6 mo6 aG3x f16xx f16yy T6
% eqf1=f12xx+ff32xx==m2
aG2x=r2/2.*(cos(theta2));
% plot(time,aG2x)
aG2y=r2/2.*(sin(theta2));
% plot(time,aG2y)
aG4x=-O4C/2.*(thetaddot4(1,:).*sin(theta4(1,:))+abs(thetadot4(1,:)).^2.*cos(theta4(1,:)));
% plot(time,aG4x)
aG4y=O4C/2.*(thetaddot4(1,:).*cos(theta4(1,:))-abs(thetadot4(1,:)).^2.*sin(theta4(1,:)));
% plot(time,aG4y)
aG5x=O4C.*(thetaddot4(1,:).*sin(theta4(1,:))+abs(thetadot4(1,:)).^2.*(cos(theta4(1,:))))+CD/2.*(thetaddot11(1,:).*(sin(theta11(1,:)))+thetadot11(1,:).^2.*(cos(theta11(1,:))));
aG5x=-aG5x;
% plot(time,aG5x)
aG5y=-O4C.*(thetaddot4(1,:)).*(cos(theta4(1,:)))+(thetadot4(1,:)).^2.*((cos(theta4(1,:))))+CD/2.*(thetaddot11(1,:).*(cos(theta11(1,:)))+thetadot11(1,:).^2.*(sin(theta11(1,:))));
% plot(time,aG5y)
aBx=O2B.*(cos(theta2(1,:)));
aG3x=aBx;
% plot(time,aG3x)
aBy=O2B.*(sin(theta2(1,:)));
aG3y=aBy;
% plot(time,aG3y)
aG6=-rddot9;
% plot(time,aG6);
% plot(time,rdot9)
rdot9=-rdot9;


IG2=1/3*m2*O2B^2;
IG4=1/3*m4*O4C;
IG5=1/12*m5*CD^2;
IG3=1/12*m3*(60^2+100^2);
IG6=1/12*m6*(200^2+60^2);
p=2000;



for j=1:n
    
eqf1=f12xx+f32xx==m2*aG2x(j);
eqf2=f12yy+f32yy==m2*aG2y(j);
eqf3=T2+O2B*(cos(theta2(j)))*f32yy-O2B*(sin(theta2(j)))*f32xx==0;



eqf5=f14xx+f54xx-f34xx==m4*aG4x(j);
eqf6=f14yy+f34yy+f54yy==m4*aG4y(j);
eqf7=-T3+f34xx*r4(j)*sin(theta4(j))+f34yy*r4(j)*cos(theta4(j))+(O4C*(cos(theta4(j)))*f54yy-O4C*sin(theta4(j))*f54xx)==IG4*thetaddot4(j);


eqf8=f65xx-f54xx==m5*aG5x(j);
eqf9=f65yy-f54yy==m5*aG5y(j);
eqf10=(-G5Dx(j)*f65yy-G5Dy(j)*f65xx)+(-G5Cx(j)*f54yy-G5Cy(j)*f54xx)==IG5*thetaddot11(j);


eqf11=-f32xx+f34xx==m3*aG3x(j);
eqf12=-f32yy-f34yy==m3*aG3y(j);
eqf13=T3==IG3*thetaddot4(j);

eqf14=(f16xx)==abs(mo*f16yy);
eqf15=-f16yy-f65yy==0;
if rdot9(j)<0
    eqf16=f16xx+p-f65xx==m6*rddot9(j);
else
     eqf16=-f16xx+p-f65xx==m6*rddot9(j);
end

eqf17=T6+p*130==0;
% eqf18=f34xx/f34yy==abs(tan(theta4(j)));
    
eqss=[eqf1,eqf2,eqf3,eqf5,eqf6,eqf7,eqf8,eqf9,eqf10,eqf11,eqf11,eqf12,eqf13,eqf14,eqf15,eqf16,eqf17];
vars=[f14xx,f14yy,f12xx,f12yy,f34xx,f34yy,f32xx,f32yy,f54xx,f54yy,f65xx,f65yy,f16xx,f16yy,T6,T3,T2];
sff=solve(eqss,vars);

f14x(j)=sff.f14xx;
f14y(j)=sff.f14yy;
f12x(j)=sff.f12xx;
f12y(j)=sff.f12yy;
f34x(j)=sff.f34xx;
f34y(j)=sff.f34yy;
f32x(j)=sff.f32xx;
f32y(j)=sff.f32yy;
f54x(j)=sff.f54xx;
f54y(j)=sff.f54yy;
f65x(j)=sff.f65xx;
f65y(j)=sff.f65yy;
f16x(j)=sff.f16xx;
f16y(j)=sff.f16yy;
t6(j)=sff.T6;
t3(j)=sff.T3;
t2(j)=sff.T2;









end


figure;

subplot(3,3,1)
plot(time,-rddot9,'black')
grid on
xlabel('time')
ylabel('acceleration(mm/s^2)')
title('acceleration of part 6')


subplot(3,3,2)
plot(time,-rdot9,'black')
grid on
xlabel('time')
ylabel('velocity (mm/s)')
title('velocity of part 6')

subplot(3,3,3)
plot(time,-r9,'black')
grid on
xlabel('time')
ylabel('position (mm)')
title('position of part 6')

subplot(3,3,4)
plot(time,theta11,'black')
grid on
xlabel('time')
ylabel('theta (radian)')
title('angle of part 5')

subplot(3,3,5)
plot(time,thetadot11,'black')
grid on
xlabel('time')
ylabel('thetadot (radian/s)')
title('angular velocity of part 5')

subplot(3,3,6)
plot(time,thetaddot11,'black')
grid on
xlabel('time')
ylabel('thetaddot (radian/s)')
title('angular acceleration of part 5')

subplot(3,3,7)
plot(time,theta10,'black')
grid on
xlabel('time')
ylabel('theta (radian)')
title('angle of part 4')


subplot(3,3,8)
plot(time,thetadot10,'black')
grid on
xlabel('time')
ylabel('thetadot (radian/s)')
title('angular velocity of part 4')


subplot(3,3,9)
plot(time,thetaddot10,'black')
grid on
xlabel('time')
ylabel('thetaddot (radian/s^2)')
title('angular acceleration of part 4')

figure;
plot(time,f12x)
grid on
title('x-force of O2')
xlabel('time')
ylabel('newton')

figure;
plot(time,f12y)
title('y-force of O2')
xlabel('time')
ylabel('newton')


figure;
plot(time,f32y)
title('y-force of B')
xlabel('time')
ylabel('newton')


figure;
plot(time,f32x)
title('x-force of B')
xlabel('time')
ylabel('newton')


figure;
plot(time,f14x)
title('x-force of O4')
xlabel('time')
ylabel('newton')

figure;
plot(time,f14y)
title('y-force of O4')
xlabel('time')
ylabel('newton')

figure;
plot(time,f54x)
title('x-force of c')
xlabel('time')
ylabel('newton')



figure;
plot(time,f54y)
title('y-force of c')
xlabel('time')
ylabel('newton')


figure;
plot(time,f16y)
title('y-force of d')
xlabel('time')
ylabel('newton')

figure;
plot(time,f16x)
title('x-force of d')
xlabel('time')
ylabel('newton')

toc





