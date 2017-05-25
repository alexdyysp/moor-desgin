%算yita
g=9.8;                  
Lc=22.05; lc=105e-3;  
mc=7*105e-3;            
nc=ceil(Lc/lc);
ro=1.025e3;            
v=36;                
M=3000;               
DH=0;
u=1.5;
H=20;
ee=zeros(1,10001);
zz=zeros(1,10001);

for yita=0:0.0001:1
	z=round((yita-0)*10000+1);
	h=[2,1,1,1,1,1, lc*ones(1,nc)];	%各物体长
	m=[1000,10,10,10,10,100, mc*ones(1,nc)];	%节点质量
	d=[2,5e-2,5e-2,5e-2,5e-2,0.3,zeros(1,nc)];	%节点半径
	N=length(h);	%共有多少节点
	[FT,FW,theta]=deal(zeros(1,N));	%初始化各变量
	FB=pi*(d/2).^2.*h*ro*g - m*g; 	%净浮力
	FB(1)=pi*(d(1)/2).^2.*h(1)*ro*g*yita-m(1)*g;	%浮标的净浮力	
	S=h.*d;	%受水阻力或风力的面积
	FW(1)=0.625*(1-yita)*S(1)*v^2+374*S(1)*u^2*(1-yita);	%第一个横向力，浮标风力
	FT(1)=0;	%第一个相互作用力，浮标不受另外物体作用
	theta(1)=0;	%第一个物体的倾角，浮标无倾角
	FT(2)=sqrt(FB(1)^2+FW(1)^2);	%第二个相互作用力，浮标对第一个钢管的弹力
	theta(2)=acos(FB(1)/FT(2));	%第二个物体倾角，第一个钢管的倾角
		for i=2:N-1
			fx=FT(i)*sin(theta(i))+374*S(i)*u^2*(1-yita)*sin(theta(i));	%每个相互作用力的X向分量
			fy=FB(i)+FT(i)*cos(theta(i));	%每个相互作用力的Y向分量
				if i==6
					fy=fy-M*g;		%第6个节点是有重球悬挂，需要加上
				end
  			FT(i+1)=sqrt(fx^2+fy^2);	%锚链上，i+1个相互作用力由i个X与Y向力勾股定理求得
			theta(i+1)=acos(fy/FT(i+1));	%i+1个倾角由i个的Y向力与i+1个相互作用力求得
			if theta(i+1)>pi/2
				theta(i+1)=pi/2; 	%可能有锚链会卧底与海床上，倾角为90度
			end
		end
	x=h.*sin(theta);
	x=cumsum([0 fliplr(x)]);	%x的坐标及反向累加运算
	y=h.*cos(theta);
	y=cumsum([0 fliplr(y)]);	%y的坐标及反向累加运算
	dh(1)=h(1)*yita;
	DH=y(end-1)+dh(1);
	ee(z)=abs(DH-H);
	zz(z)=yita;
end

for t=1:1:10001
	if ee(t)==min(ee)
		yita0=zz(t);ee0=ee(t);	
		fprintf('%f,%f,%d\n',min(ee),zz(t),t)	%输出最小误差的yita
	end
end

%主
g=9.8;                  % 重力加速度
Lc=22.05; lc=105e-3;  % 锚链长度和每环长度
mc=7*105e-3;            % 线质量
nc=ceil(Lc/lc);	  % 共有多少环
ro=1.025e3;            % 海水密度
v=36;                  % 24 m/s
M=3000;                 % 重物球
u=1.5;
H=20;	%水深
yita=yita0;

h=[2,1,1,1,1,1, lc*ones(1,nc)];	%各物体长
m=[1000,10,10,10,10,100, mc*ones(1,nc)];	%节点质量
d=[2,5e-2,5e-2,5e-2,5e-2,0.3,zeros(1,nc)];	%节点半径

N=length(h);	%共有多少节点
[FT,FW,theta]=deal(zeros(1,N));	%初始化各变量

FB=pi*(d/2).^2.*h*ro*g-m*g; 	%净浮力
FB(1)=pi*(d(1)/2).^2.*h(1)*ro*g*yita-m(1)*g;	%浮标的净浮力	

S=h.*d;	%受水阻力或风力的法相面积
FW(1)=0.625*(1-yita)*S(1)*v^2+374*S(1)*u^2*(1-yita);	%第一个横向力，浮标风力
FT(1)=0;	%第一个相互作用力，浮标不受另外物体作用
theta(1)=0;	%第一个物体的倾角，浮标无倾角

FT(2)=sqrt(FB(1)^2+FW(1)^2);	%第二个相互作用力，浮标对第一个钢管的弹力
theta(2)=acos(FB(1)/FT(2));	%第二个物体倾角，第一个钢管的倾角

for i=2:N-1
	fx=FT(i)*sin(theta(i))+374*S(i)*u^2*(1-yita)*sin(theta(i));	%每个相互作用力的X向分量
	fy=FB(i)+FT(i)*cos(theta(i));	%每个相互作用力的Y向分量

	if i==6
		fy=fy-M*g;		%第6个节点是有重球悬挂，需要加上
	end
   
	FT(i+1)=sqrt(fx^2+fy^2);	%锚链上，i+1个相互作用力由i个X与Y向力勾股定理求得
	theta(i+1)=acos(fy/FT(i+1));	%i+1个倾角由i个的Y向力与i+1个相互作用力求得
	if theta(i+1)>pi/2
		theta(i+1)=pi/2; 	%可能有锚链会卧底与海床上，倾角为90度
	end
end

y=h.*cos(theta); 
y=cumsum([0 fliplr(y)]);	%y的坐标及反向累加运算
x=h.*sin(theta);
x=cumsum([0 fliplr(x)]);	%x的坐标及反向累加运算
THETA=theta*180/3.14;	%换算成角度
%作图，先画浮标
x1=x(end)-1;x2=x1+2;x3=x(end-1)+1;x4=x3-2;	%画浮标，求4个顶点横坐标
y1=y(end);y2=y1;y3=y(end-1);y4=y3;	%求4个顶点纵坐标
fill([x1 x2 x3 x4], [y1 y2 y3 y4], 'g')	;		%画出封闭曲线的浮标
hold on
plot([0,max(x)+5],[H,H],'--');	%海平面线
plot(x(end-5:end-1),y(end-5:end-1),'o-','linewidth',2, 'markersize',3);	%画出4根钢管和钢桶
plot(x(end-6:end-5),y(end-6:end-5),'r-','linewidth',3, 'markersize',4);	%画出4根钢管和钢桶
plot(x(1:end-6),y(1:end-6),'-k','linewidth',2);		%画出锚链曲线
plot([x(end-6) x(end-6)],[y(end-6) y(end-6)-2],'-b','linewidth',1);	%
plot(x(end-6),y(end-6)-2,'.r','linewidth',1,'markersize',50);
ylim([0,H+h(1)*(1-yita)]); xlabel('x'); ylabel('y');	%确定图像范围上下限
box on	%加边框
axis image	%使画出的图紧紧围绕周边


	
