%算yita
g=9.8;                  
Lc=22.05; lc=105e-3;  
mc=7*105e-3;            
nc=ceil(Lc/lc);
ro=1.025e3;            
v=36;
H=18;                          
DH=0;	%初始化参数
ee=zeros(1,10001);
zz=zeros(1,10001);

for M=1770:1:1780
	tt=round((M-1770)/1+1);
for yita=0:0.0001:1
	z=round((yita-0)*10000+1);
	h=[2,1,1,1,1,1, lc*ones(1,nc)];	%各物体长
	m=[1000,10,10,10,10,100, mc*ones(1,nc)];	%节点质量
	d=[2,5e-2,5e-2,5e-2,5e-2,0.3,zeros(1,nc)];	%节点半径
	N=length(h);	%共有多少节点
	[FT,FW,theta]=deal(zeros(1,N));	%初始化各变量
	FB=pi*(d/2).^2.*h*ro*g - m*g; 	%净浮力
	FB(1)=pi*(d(1)/2).^2.*h(1)*ro*g*yita-m(1)*g;	%浮标的净浮力	
	S=h.*d;	%受水阻力或风力的法相面积
	FW(1)=0.625*(1-yita)*S(1)*v^2;	%第一个横向力，浮标风力
	FT(1)=0;	%第一个相互作用力，浮标不受另外物体作用
	theta(1)=0;	%第一个物体的倾角，浮标无倾角
	FT(2)=sqrt(FB(1)^2+FW(1)^2);	%第二个相互作用力，浮标对第一个钢管的弹力
	theta(2)=acos(FB(1)/FT(2));	%第二个物体倾角，第一个钢管的倾角
		for i=2:N-1
			fx=FT(i)*sin(theta(i));	%每个相互作用力的X向分量
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
	ee(z)=abs(DH-H);	%误差矩阵
	zz(z)=yita;		%相应Yita矩阵
end

%输出yita
for t=1:1:10001
	if ee(t)==min(ee)
		yita0(tt)=zz(t);ee0(tt)=ee(t);	
		fprintf('%f,%f,%d\n',min(ee),zz(t),t)	%输出最小误差的yita
	end
end

end

%主程序
g=9.8;                 % 重力加速度
Lc=22.05; lc=105e-3;   % 锚链长度和每环长度
mc=7*105e-3;           % 线质量
nc=ceil(Lc/lc);	       % 共有多少环
ro=1.025e3;            % 海水密度
H=18;
v=36;                  

for M=1770:1:1780  
	j=(M-1770)/1+1;
	yita=yita0(j);
	h=[2,1,1,1,1,1, lc*ones(1,nc)];	%各物体长
	m=[1000,10,10,10,10,100, mc*ones(1,nc)];	%节点质量
	d=[2, 5e-2, 5e-2, 5e-2, 5e-2, 0.3,   zeros(1,nc)];	%节点半径
	N=length(h);	%共有多少节点
	[FT, FW, theta]=deal(zeros(1,N));	%初始化各变量
	FB=pi*(d/2).^2.*h*ro*g-m*g; 	%净浮力
	FB(1)=pi*(d(1)/2).^2.*h(1)*ro*g*yita-m(1)*g;	%浮标的净浮力	
	S=h.*d;	%受水阻力或风力的法相面积
	FW(1)=0.625*(1-yita)*S(1)*v^2;	%第一个横向力，浮标风力
	FT(1)=0;	%第一个相互作用力，浮标不受另外物体作用
	theta(1)=0;	%第一个物体的倾角，浮标无倾角
	FT(2)=sqrt(FB(1)^2+FW(1)^2);	%第二个相互作用力，浮标对第一个钢管的弹力
	theta(2)=acos(FB(1)/FT(2));	%第二个物体倾角，第一个钢管的倾角
		for i=2:N-1
			fx=FT(i)*sin(theta(i));	%每个相互作用力的X向分量	
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
	THETA=theta*180/pi;	%换算成角度
	sita1(j)=THETA(6);
	sita2(j)=THETA(end);
	maxx(j)=max(x);
end

%输出图像
M=1770:1:1780 ;
plot(M,sita1,'.-')	%钢桶倾角
hold on
plot(M,sita2,'.-')	%锚链底部倾角
plot(M,maxx,'.-')	%浮动范围
plot(M,yita0,'.-')	%浸水百分比