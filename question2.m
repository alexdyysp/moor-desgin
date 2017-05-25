%��yita
g=9.8;                  
Lc=22.05; lc=105e-3;  
mc=7*105e-3;            
nc=ceil(Lc/lc);
ro=1.025e3;            
v=36;
H=18;                          
DH=0;	%��ʼ������
ee=zeros(1,10001);
zz=zeros(1,10001);

for M=1770:1:1780
	tt=round((M-1770)/1+1);
for yita=0:0.0001:1
	z=round((yita-0)*10000+1);
	h=[2,1,1,1,1,1, lc*ones(1,nc)];	%�����峤
	m=[1000,10,10,10,10,100, mc*ones(1,nc)];	%�ڵ�����
	d=[2,5e-2,5e-2,5e-2,5e-2,0.3,zeros(1,nc)];	%�ڵ�뾶
	N=length(h);	%���ж��ٽڵ�
	[FT,FW,theta]=deal(zeros(1,N));	%��ʼ��������
	FB=pi*(d/2).^2.*h*ro*g - m*g; 	%������
	FB(1)=pi*(d(1)/2).^2.*h(1)*ro*g*yita-m(1)*g;	%����ľ�����	
	S=h.*d;	%��ˮ����������ķ������
	FW(1)=0.625*(1-yita)*S(1)*v^2;	%��һ�����������������
	FT(1)=0;	%��һ���໥�����������겻��������������
	theta(1)=0;	%��һ���������ǣ����������
	FT(2)=sqrt(FB(1)^2+FW(1)^2);	%�ڶ����໥������������Ե�һ���ֹܵĵ���
	theta(2)=acos(FB(1)/FT(2));	%�ڶ���������ǣ���һ���ֹܵ����
		for i=2:N-1
			fx=FT(i)*sin(theta(i));	%ÿ���໥��������X�����
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
	ee(z)=abs(DH-H);	%������
	zz(z)=yita;		%��ӦYita����
end

%���yita
for t=1:1:10001
	if ee(t)==min(ee)
		yita0(tt)=zz(t);ee0(tt)=ee(t);	
		fprintf('%f,%f,%d\n',min(ee),zz(t),t)	%�����С����yita
	end
end

end

%������
g=9.8;                 % �������ٶ�
Lc=22.05; lc=105e-3;   % ê�����Ⱥ�ÿ������
mc=7*105e-3;           % ������
nc=ceil(Lc/lc);	       % ���ж��ٻ�
ro=1.025e3;            % ��ˮ�ܶ�
H=18;
v=36;                  

for M=1770:1:1780  
	j=(M-1770)/1+1;
	yita=yita0(j);
	h=[2,1,1,1,1,1, lc*ones(1,nc)];	%�����峤
	m=[1000,10,10,10,10,100, mc*ones(1,nc)];	%�ڵ�����
	d=[2, 5e-2, 5e-2, 5e-2, 5e-2, 0.3,   zeros(1,nc)];	%�ڵ�뾶
	N=length(h);	%���ж��ٽڵ�
	[FT, FW, theta]=deal(zeros(1,N));	%��ʼ��������
	FB=pi*(d/2).^2.*h*ro*g-m*g; 	%������
	FB(1)=pi*(d(1)/2).^2.*h(1)*ro*g*yita-m(1)*g;	%����ľ�����	
	S=h.*d;	%��ˮ����������ķ������
	FW(1)=0.625*(1-yita)*S(1)*v^2;	%��һ�����������������
	FT(1)=0;	%��һ���໥�����������겻��������������
	theta(1)=0;	%��һ���������ǣ����������
	FT(2)=sqrt(FB(1)^2+FW(1)^2);	%�ڶ����໥������������Ե�һ���ֹܵĵ���
	theta(2)=acos(FB(1)/FT(2));	%�ڶ���������ǣ���һ���ֹܵ����
		for i=2:N-1
			fx=FT(i)*sin(theta(i));	%ÿ���໥��������X�����	
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
	THETA=theta*180/pi;	%����ɽǶ�
	sita1(j)=THETA(6);
	sita2(j)=THETA(end);
	maxx(j)=max(x);
end

%���ͼ��
M=1770:1:1780 ;
plot(M,sita1,'.-')	%��Ͱ���
hold on
plot(M,sita2,'.-')	%ê���ײ����
plot(M,maxx,'.-')	%������Χ
plot(M,yita0,'.-')	%��ˮ�ٷֱ�