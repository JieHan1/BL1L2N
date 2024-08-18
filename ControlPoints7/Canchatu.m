function Canchatu
KZHD=load('控制点坐标.txt');
Forward=load('前方交会.txt');
Bundle=load('光束法.txt');
set(0,'defaultfigurecolor','w');
%% 计算残差
dX=KZHD(:,2)-Forward(:,2);   dX=dX*10;
dY=KZHD(:,3)-Forward(:,3);   dY=dY*10;
dZ=KZHD(:,4)-Forward(:,4);   dZ=dZ*40;
figure('Name','前方交会');
hold on
xlabel('X'); ylabel('Y'); title('前方交会');
for i=1:3%绘制控制点
    plot(KZHD(i,2),KZHD(i,3),'^');
end
for i=4:42%绘制检核点
    plot(KZHD(i,2),KZHD(i,3),'*');
end
for i=4:42
    line([KZHD(i,2),KZHD(i,2)+dX(i,1)],[KZHD(i,3),KZHD(i,3)+dY(i,1)]);
    line([KZHD(i,2),KZHD(i,2)],[KZHD(i,3),KZHD(i,3)+dZ(i,1)]);
end
box on
%% 
%% 计算残差
dX=KZHD(:,2)-Bundle(:,2);   
dY=KZHD(:,3)-Bundle(:,3);   
dZ=KZHD(:,4)-Bundle(:,4);   
dX=dX*10;
dY=dY*10;
dZ=dZ*40;
figure('Name','光束法');
hold on
xlabel('X'); ylabel('Y'); title('光束法平差');
for i=1:1%绘制控制点
    plot(KZHD(i,2),KZHD(i,3),'^');
end
for i=4:42%绘制检核点
    plot(KZHD(i,2),KZHD(i,3),'*');
end
for i=4:42
    line([KZHD(i,2),KZHD(i,2)+dX(i,1)],[KZHD(i,3),KZHD(i,3)+dY(i,1)]);
    line([KZHD(i,2),KZHD(i,2)],[KZHD(i,3),KZHD(i,3)+dZ(i,1)]);
end
box on

end

