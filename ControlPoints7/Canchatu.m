function Canchatu
KZHD=load('���Ƶ�����.txt');
Forward=load('ǰ������.txt');
Bundle=load('������.txt');
set(0,'defaultfigurecolor','w');
%% ����в�
dX=KZHD(:,2)-Forward(:,2);   dX=dX*10;
dY=KZHD(:,3)-Forward(:,3);   dY=dY*10;
dZ=KZHD(:,4)-Forward(:,4);   dZ=dZ*40;
figure('Name','ǰ������');
hold on
xlabel('X'); ylabel('Y'); title('ǰ������');
for i=1:3%���ƿ��Ƶ�
    plot(KZHD(i,2),KZHD(i,3),'^');
end
for i=4:42%���Ƽ�˵�
    plot(KZHD(i,2),KZHD(i,3),'*');
end
for i=4:42
    line([KZHD(i,2),KZHD(i,2)+dX(i,1)],[KZHD(i,3),KZHD(i,3)+dY(i,1)]);
    line([KZHD(i,2),KZHD(i,2)],[KZHD(i,3),KZHD(i,3)+dZ(i,1)]);
end
box on
%% 
%% ����в�
dX=KZHD(:,2)-Bundle(:,2);   
dY=KZHD(:,3)-Bundle(:,3);   
dZ=KZHD(:,4)-Bundle(:,4);   
dX=dX*10;
dY=dY*10;
dZ=dZ*40;
figure('Name','������');
hold on
xlabel('X'); ylabel('Y'); title('������ƽ��');
for i=1:1%���ƿ��Ƶ�
    plot(KZHD(i,2),KZHD(i,3),'^');
end
for i=4:42%���Ƽ�˵�
    plot(KZHD(i,2),KZHD(i,3),'*');
end
for i=4:42
    line([KZHD(i,2),KZHD(i,2)+dX(i,1)],[KZHD(i,3),KZHD(i,3)+dY(i,1)]);
    line([KZHD(i,2),KZHD(i,2)],[KZHD(i,3),KZHD(i,3)+dZ(i,1)]);
end
box on

end

