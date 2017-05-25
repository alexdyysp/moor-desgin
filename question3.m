%��yita
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
	h=[2,1,1,1,1,1, lc*ones(1,nc)];	%�����峤
	m=[1000,10,10,10,10,100, mc*ones(1,nc)];	%�ڵ�����
	d=[2,5e-2,5e-2,5e-2,5e-2,0.3,zeros(1,nc)];	%�ڵ�뾶
	N=length(h);	%���ж��ٽڵ�
	[FT,FW,theta]=deal(zeros(1,N));	%��ʼ��������
	FB=pi*(d/2).^2.*h*ro*g - m*g; 	%������
	FB(1)=pi*(d(1)/2).^2.*h(1)*ro*g*yita-m(1)*g;	%����ľ�����	
	S=h.*d;	%��ˮ��������������
	FW(1)=0.625*(1-yita)*S(1)*v^2+374*S(1)*u^2*(1-yita);	%��һ�����������������
	FT(1)=0;	%��һ���໥�����������겻��������������
	theta(1)=0;	%��һ���������ǣ����������
	FT(2)=sqrt(FB(1)^2+FW(1)^2);	%�ڶ����໥������������Ե�һ���ֹܵĵ���
	theta(2)=acos(FB(1)/FT(2));	%�ڶ���������ǣ���һ���ֹܵ����
		for i=2:N-1
			fx=FT(i)*sin(theta(i))+374*S(i)*u^2*(1-yita)*sin(theta(i));	%ÿ���໥��������X�����
			fy=FB(i)+FT(i)*cos(theta(i));	%ÿ���໥��������Y�����
				if i==6
					fy=fy-M*g;		%��6���ڵ������������ң���Ҫ����
				end
  			FT(i+1)=sqrt(fx^2+fy^2);	%ê���ϣ�i+1���໥��������i��X��Y�������ɶ������
			theta(i+1)=acos(fy/FT(i+1));	%i+1�������i����Y������i+1���໥���������
			if theta(i+1)>pi/2
				theta(i+1)=pi/2; 	%������ê�����Ե��뺣���ϣ����Ϊ90��
			end
		end
	x=h.*sin(theta);
	x=cumsum([0 fliplr(x)]);	%x�����꼰�����ۼ�����
	y=h.*cos(theta);
	y=cumsum([0 fliplr(y)]);	%y�����꼰�����ۼ�����
	dh(1)=h(1)*yita;
	DH=y(end-1)+dh(1);
	ee(z)=abs(DH-H);
	zz(z)=yita;
end

for t=1:1:10001
	if ee(t)==min(ee)
		yita0=zz(t);ee0=ee(t);	
		fprintf('%f,%f,%d\n',min(ee),zz(t),t)	%�����С����yita
	end
end

%��
g=9.8;                  % �������ٶ�
Lc=22.05; lc=105e-3;  % ê�����Ⱥ�ÿ������
mc=7*105e-3;            % ������
nc=ceil(Lc/lc);	  % ���ж��ٻ�
ro=1.025e3;            % ��ˮ�ܶ�
v=36;                  % 24 m/s
M=3000;                 % ������
u=1.5;
H=20;	%ˮ��
yita=yita0;

h=[2,1,1,1,1,1, lc*ones(1,nc)];	%�����峤
m=[1000,10,10,10,10,100, mc*ones(1,nc)];	%�ڵ�����
d=[2,5e-2,5e-2,5e-2,5e-2,0.3,zeros(1,nc)];	%�ڵ�뾶

N=length(h);	%���ж��ٽڵ�
[FT,FW,theta]=deal(zeros(1,N));	%��ʼ��������

FB=pi*(d/2).^2.*h*ro*g-m*g; 	%������
FB(1)=pi*(d(1)/2).^2.*h(1)*ro*g*yita-m(1)*g;	%����ľ�����	

S=h.*d;	%��ˮ����������ķ������
FW(1)=0.625*(1-yita)*S(1)*v^2+374*S(1)*u^2*(1-yita);	%��һ�����������������
FT(1)=0;	%��һ���໥�����������겻��������������
theta(1)=0;	%��һ���������ǣ����������

FT(2)=sqrt(FB(1)^2+FW(1)^2);	%�ڶ����໥������������Ե�һ���ֹܵĵ���
theta(2)=acos(FB(1)/FT(2));	%�ڶ���������ǣ���һ���ֹܵ����

for i=2:N-1
	fx=FT(i)*sin(theta(i))+374*S(i)*u^2*(1-yita)*sin(theta(i));	%ÿ���໥��������X�����
	fy=FB(i)+FT(i)*cos(theta(i));	%ÿ���໥��������Y�����

	if i==6
		fy=fy-M*g;		%��6���ڵ������������ң���Ҫ����
	end
   
	FT(i+1)=sqrt(fx^2+fy^2);	%ê���ϣ�i+1���໥��������i��X��Y�������ɶ������
	theta(i+1)=acos(fy/FT(i+1));	%i+1�������i����Y������i+1���໥���������
	if theta(i+1)>pi/2
		theta(i+1)=pi/2; 	%������ê�����Ե��뺣���ϣ����Ϊ90��
	end
end

y=h.*cos(theta); 
y=cumsum([0 fliplr(y)]);	%y�����꼰�����ۼ�����
x=h.*sin(theta);
x=cumsum([0 fliplr(x)]);	%x�����꼰�����ۼ�����
THETA=theta*180/3.14;	%����ɽǶ�
%��ͼ���Ȼ�����
x1=x(end)-1;x2=x1+2;x3=x(end-1)+1;x4=x3-2;	%�����꣬��4�����������
y1=y(end);y2=y1;y3=y(end-1);y4=y3;	%��4������������
fill([x1 x2 x3 x4], [y1 y2 y3 y4], 'g')	;		%����������ߵĸ���
hold on
plot([0,max(x)+5],[H,H],'--');	%��ƽ����
plot(x(end-5:end-1),y(end-5:end-1),'o-','linewidth',2, 'markersize',3);	%����4���ֹܺ͸�Ͱ
plot(x(end-6:end-5),y(end-6:end-5),'r-','linewidth',3, 'markersize',4);	%����4���ֹܺ͸�Ͱ
plot(x(1:end-6),y(1:end-6),'-k','linewidth',2);		%����ê������
plot([x(end-6) x(end-6)],[y(end-6) y(end-6)-2],'-b','linewidth',1);	%
plot(x(end-6),y(end-6)-2,'.r','linewidth',1,'markersize',50);
ylim([0,H+h(1)*(1-yita)]); xlabel('x'); ylabel('y');	%ȷ��ͼ��Χ������
box on	%�ӱ߿�
axis image	%ʹ������ͼ����Χ���ܱ�


	
