%% 数值解
clc;clear;close all

% 输入初始条件
% 预设条件
L=1;  % 棒子长度
Nx=180;  % 水平节点数
x=linspace(0,1,Nx);
dx=L/(Nx-1);  % 格点间距
K=1.e-2;  % 热传导系数
A=5;  % 振幅
dt=10;  % 时间步长
mu=(K*dt)/(dx^2);

% 动图参数
Fname='D:\研究生作业\数模作业1\';
GifName='FTCS.gif';
delay=0.1;  % 延迟时间

% 初始时刻，t=0时温度分布
% T0=A*cos((pi/L).*x);  % 按照初始条件得出T0温度分布
% T0=0.5*A*(1+cos(1*x/L*pi));  % PPT中T0温度分布
T0=A*exp( - (x-0.5*L).^4*1.e5 );

% 绘制并写入第一张图
figure
area(x,T0)
grid minor
drawnow
set(gca,'xlim',[0 L],'ylim',[0 A])
title('第1步')
xlabel('x')
ylabel('T')
print([Fname,'1','-dpng']);

V=imread([Fname,num2str(count),'.png']);
[X,map]=rgb2ind(V,256);
imwrite(X,map,[Fname,GifName],'gif','LoopCount',inf,'DelayTime',delay)

% 初始参数
count=0;

T1=zeros(size(T0));  % 储存t=1时刻温度计算值
T2=zeros(size(T0));  % 储存t=2时刻温度计算值


% 开始计算
% 使用FTCS(Forward Time Central Space)
xi=2:Nx-1;  % 头尾格点由初始条件和边界条件给出，xi为需计算格点
while count<5000
    count=count+1;

    T1(xi)=(1-2*mu)*T0(xi) + mu*( T0(xi+1)+T0(xi-1) );
    T1(1)=T1(2);
    T1(end)=T1(end-1);
    
    T0=T1;

    if mod(count,50)==0
        figure(1)
        grid minor
        hold on
        title(['第',num2str(count),'步'])
        xlabel('x')
        ylabel('T')

        area(x,T0)
        hold on
        drawnow
        set(gca,'xlim',[0 L],'ylim',[0 A])
        print([Fname,num2str(count),'-dpng']);

        V=imread([Fname,num2str(count),'.png']);
        [X,map]=rgb2ind(V,256);
        imwrite(X,map,[Fname,GifName],'gif','WriteMode','append','DelayTime',delay)
        
    end

end


%% 解析解
clc;clear;close all

% 输入初始条件
% 预设条件
L=1;  % 棒子长度
Nx=180;  % 水平节点数
x=linspace(0,1,Nx);
dx=L/(Nx-1);  % 格点间距
K=1.e-2;  % 热传导系数
A=5;  % 振幅
dt=10;  % 时间步长
mu=(K*dt)/(dx^2);

% 动图参数
Fname='D:\研究生作业\数模作业1\';
GifName='JieXi.gif';
delay=0.1;  % 延迟时间

% 初始时刻，t=0时温度分布
t=0;
T0=A*exp(-1*K*(pi/L).^2*t)*cos(pi/L*x);

% 绘制并写入第一张图
figure
area(x,T0)
grid minor
drawnow
set(gca,'xlim',[0 L],'ylim',[0 A])
title('第1步')
xlabel('x')
ylabel('T')
print([Fname,'1','-dpng']);

V=imread([Fname,num2str(count),'.png']);
[X,map]=rgb2ind(V,256);
imwrite(X,map,[Fname,GifName],'gif','LoopCount',inf,'DelayTime',delay)

% 开始计算
% 使用数值解
xi=2:Nx-1;  % 头尾格点由初始条件和边界条件给出，xi为需计算格点
while t<50000
    t=t+dt;

    T=A*exp(-1*K*(pi/L).^2*t)*cos(pi/L*x);

    if mod(t,500)==0
        figure(1)
        grid minor
        hold on
        title(['第',num2str(t/10),'步'])
        xlabel('x')
        ylabel('T')

        area(x,T)
        hold on
        drawnow
        set(gca,'xlim',[0 L],'ylim',[0 A])
        print([Fname,num2str(count),'-dpng']);

        V=imread([Fname,num2str(count),'.png']);
        [X,map]=rgb2ind(V,256);
        imwrite(X,map,[Fname,GifName],'gif','WriteMode','append','DelayTime',delay)
        
    end

end




















