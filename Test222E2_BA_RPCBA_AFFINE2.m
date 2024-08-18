function BLH=Test222E2_BA_RPCBA_AFFINE2

tic
%%%%%%%%%%************************************************%%%%%%%%%%%

%基于基准影像的区域网平差：
%1.将第一景影像作为基准，也就是其前像点坐标作为真值，无像方改正；
%2.利用第一景原始RPC进行计算，获得区域网平差结果；
%3.利用第一景影像RPC经过像方放射变换的值进行计算，获得区域网平差结果
%4.精度评定采用像方残差的形式

%%%%%%%%%%************************************************%%%%%%%%%%%

set(0,'defaultfigurecolor','w');
%%%%------数据读取-------%%%%
% %读取左片RPC参数
close all;
% [filename,pathname]=uigetfile('.txt','选择相应的RPC文件');
% fnopen=fopen(strcat(pathname,filename),'rt');
% if(fnopen==-1)
%     msgbox('Input file or path is not correct',waring','warm');
% end
% RPCsLR(1,1) = readrpc([pathname,filename]);
RPCsLR(1,1) = readrpc('zy301a_fwd_381_6138_rpc.txt');

%读取右片RPC参数
% [filename,pathname]=uigetfile('.txt','选择相应的RPC文件');
% fnopen=fopen(strcat(pathname,filename),'rt');
% if(fnopen==-1)
%     msgbox('Input file or path is not correct',waring','warm');
% end
% RPCsLR(2,1) = readrpc([pathname,filename]);
RPCsLR(2,1) = readrpc('zy301a_bwd_381_6138_rpc.txt');

%读取左右像点坐标
% [filename,pathname]=uigetfile('.txt','选择同名点坐标');
% fnopen=fopen(strcat(pathname,filename),'rt');
% if(fnopen==-1)
%     msgbox('Input file or path is not correct',waring','warm');
% end

fnopen=fopen('平高7个均匀\多同名点.txt','rt');
nIP=fscanf(fnopen,'%d',1); %像点个数
for i=1:nIP
    fscanf(fnopen,'%d',1);
    Left_rc(i,1)=fscanf(fnopen,'%f',1);
    Left_rc(i,2)=fscanf(fnopen,'%f',1);
    Right_rc(i,1)=fscanf(fnopen,'%f',1);
    Right_rc(i,2)=fscanf(fnopen,'%f',1);
    ImgUV(i,1)=Left_rc(i,1);   ImgUV(i,2)=Left_rc(i,2);
    ImgUV(i,3)=Right_rc(i,1);  ImgUV(i,4)=Right_rc(i,2);
end

%%%%读取检核点
% [filename,pathname]=uigetfile('.txt','选择检核点坐标');
% fnopen=fopen(strcat(pathname,filename),'rt');
% if(fnopen==-1)
%     msgbox('Input file or path is not correct',waring','warm');
% end
fnopen=fopen('平高7个均匀\检核点.txt','rt');
nChk=fscanf(fnopen,'%d',1);
for i=1:nChk
    fscanf(fnopen,'%d',1);
    chkBLH(i,1)=fscanf(fnopen,'%f',1);
    chkBLH(i,2)=fscanf(fnopen,'%f',1);
    chkBLH(i,3)=fscanf(fnopen,'%f',1);
end
fclose('all');
%%%%------数据读取完毕-------%%%%

%%直接利用RPC计算地面坐标并与实际值比较:包括控制点+待定点%
for i=1:nIP
     BBB=RFM_ForwardResection(RPCsLR(1,1),RPCsLR(2,1),Left_rc(i,:),Right_rc(i,:));
     XYZ(i,1)=BBB(1,1);  XYZ(i,2)=BBB(1,2);   XYZ(i,3)=BBB(1,3);
end



fid=fopen('Result_BL1L2N.txt','w');       %结果分析文件




if(fid==-1)
    msgbox('Input file or path is not correct','waring','warm');
end
fprintf(fid,'******************************自由网_直接定向结果分析******************************\n');
fprintf(fid,'前方交会地面点坐标XYZ：\n');
fprintf(fid,'点号\t\t\tX\t\t\t\t\t\tY\t\t\t\t\t\tZ\n'); 
for i=1:nIP  %chkXYZ：检核点XYZ
    temp_XYZ=BLH2XYZ(XYZ(i,1),XYZ(i,2),XYZ(i,3));
    resectXYZ(i,1)=temp_XYZ(1,1);  % totalXYZ: 控制点XYZ+待定点XYZ 
    resectXYZ(i,2)=temp_XYZ(1,2);
    resectXYZ(i,3)=temp_XYZ(1,3);
    fprintf(fid,'%d   %20.8f   %20.8f   %20.8f\n',i,resectXYZ(i,1),resectXYZ(i,2),resectXYZ(i,3));
end

fprintf(fid,'地面点点坐标XYZ：\n');
fprintf(fid,'点号\t\t\tX\t\t\t\t\t\tY\t\t\t\t\t\tZ\n'); 

for i=1:nChk  %chkXYZ：检核点XYZ
    temp_XYZ=BLH2XYZ(chkBLH(i,1),chkBLH(i,2),chkBLH(i,3));
    totalXYZ(i,1)=temp_XYZ(1,1);  % totalXYZ: 控制点XYZ+待定点XYZ 
    totalXYZ(i,2)=temp_XYZ(1,2);
    totalXYZ(i,3)=temp_XYZ(1,3);
    fprintf(fid,'%d   %20.8f   %20.8f   %20.8f\n',i,totalXYZ(i,1),totalXYZ(i,2),totalXYZ(i,3));
end

fprintf(fid,'控制点+检核点坐标与前方交会结果比较：\n');
fprintf(fid,'点号\t\t\tdX\t\t\t\t\t\tdY\t\t\t\t\t\tdZ\n');  %检核点比较
for i=1:nIP
    resect_dXYZ(i,1)=totalXYZ(i,1)-resectXYZ(i,1);
    resect_dXYZ(i,2)=totalXYZ(i,2)-resectXYZ(i,2);
    resect_dXYZ(i,3)=totalXYZ(i,3)-resectXYZ(i,3);
    fprintf(fid,'%d   %20.8f   %20.8f   %20.8f\n',i,resect_dXYZ(i,1),resect_dXYZ(i,2),resect_dXYZ(i,3));
end

Xmax=max(abs(resect_dXYZ(:,1)));    Ymax=max(abs(resect_dXYZ(:,2)));
Zmax=max(abs(resect_dXYZ(:,3)));
Xmean=mean(resect_dXYZ(:,1));
Ymean=mean(resect_dXYZ(:,2));
Zmean=mean(resect_dXYZ(:,3));

fprintf(fid,'物方最大误差\n');
fprintf(fid,'%20.8f   %20.8f   %20.8f\n',Xmax,Ymax,Zmax);

XRMS=sqrt((sum(resect_dXYZ(:,1).*resect_dXYZ(:,1)))/nChk);    YRMS=sqrt((sum(resect_dXYZ(:,2).*resect_dXYZ(:,2)))/nChk);
ZRMS=sqrt((sum(resect_dXYZ(:,3).*resect_dXYZ(:,3)))/nChk);  
fprintf(fid,'物方中误差\n');
fprintf(fid,'%20.8f   %20.8f   %20.8f\n',XRMS,YRMS,ZRMS);

fprintf(fid,'物方误差均值\n');
fprintf(fid,'%20.8f   %20.8f   %20.8f\n',Xmean,Ymean,Zmean);
fprintf(fid,'******************************直接定向结果分析end******************************\n');
fprintf(fid,'\n\n');
%———————————————————直接定向结果分析完毕———————————————————————
RP=zeros(1,6);
re0=RP(1,1);re1=RP(1,2);re2=RP(1,3);rf0=RP(1,4);rf1=RP(1,5);rf2=RP(1,6);
%le0=RP(1,1);le1=RP(1,2);le2=RP(1,3);lf0=RP(1,4);lf1=RP(1,5);lf2=RP(1,6);
le0=6.57488519;le1=-0.00014762;le2=-0.00003674;lf0=22.95222291;lf1=0.00013474;lf2=0.00008032;  %基准影像以进行放射变换

%------计算仿射变换系数初值（6个）end-----%	

%------光束法平差主程序--------%
count=0;
MaxVale=1000;  LimitVal=0.0000001;
NumGcp=2;  %确定像方变换模型
figure; %L曲线
tic
while count<50&&MaxVale>LimitVal
     if NumGcp==0||NumGcp==1   %只有平移
        A=zeros(2*nIP,2); %每幅影像有两个平移参数
    elseif NumGcp==2  %平移和缩放
        A=zeros(2*nIP,4);%每幅影像有两个平移参数和两个缩放参数
    elseif NumGcp>=3
        A=zeros(2*nIP,6);  %仿射变换
     end
    B=zeros(4*nIP,3*nChk);
    L=zeros(4*nIP,1);
    Mx=zeros(6,1);      %系统误差补偿改正数
    Mt=zeros(3*nChk,1);    %待定点坐标改正数
    uv1=BackProj_RFM(RPCsLR(1,1),XYZ);  %用控制点反算像点坐标
    uv2=BackProj_RFM(RPCsLR(2,1),XYZ);  %用控制点反算像点坐标
    for i=1:nIP  %uvF：反投影仿射变换后坐标
        uvF(i,1)=uv1(i,1)+le0+le1*uv1(i,1)+le2*uv1(i,2);
        uvF(i,2)=uv1(i,2)+lf0+lf1*uv1(i,1)+lf2*uv1(i,2);
        uvF(i,3)=uv2(i,1)+re0+re1*uv2(i,1)+re2*uv2(i,2);
        uvF(i,4)=uv2(i,2)+rf0+rf1*uv2(i,1)+rf2*uv2(i,2);
        duv(i,1)=ImgUV(i,1)-uvF(i,1);
        duv(i,2)=ImgUV(i,2)-uvF(i,2);
        duv(i,3)=ImgUV(i,3)-uvF(i,3);
        duv(i,4)=ImgUV(i,4)-uvF(i,4);
    end
    for i=1:nIP
        j=4*i;
        %-------------左像-------------------%
        Xt=(XYZ(i,1)-RPCsLR(1,1).XOffset)/RPCsLR(1,1).XScale;
        Yt=(XYZ(i,2)-RPCsLR(1,1).YOffset)/RPCsLR(1,1).YScale;
        Zt=(XYZ(i,3)-RPCsLR(1,1).ZOffset)/RPCsLR(1,1).ZScale;
        for k=1:4
            PP(k)=NCHLCalP(RPCsLR(1,1).P{k},Xt,Yt,Zt);
            dPdX(k)=CaldPdX(RPCsLR(1,1).P{k},Xt,Yt,Zt);   %P1、P2、P3、P4对X求偏导
            dPdY(k)=CaldPdY(RPCsLR(1,1).P{k},Xt,Yt,Zt);   %P1、P2、P3、P4对Y求偏导
            dPdZ(k)=CaldPdZ(RPCsLR(1,1).P{k},Xt,Yt,Zt);	%P1、P2、P3、P4对Z求偏导
        end
        if NumGcp==0||NumGcp==1   %只有平移
            A(j-3,1)=0;
            A(j-2,2)=0;
        elseif NumGcp==2  %平移和缩放
            A(j-3,1)=0;
            A(j-3,2)=0;
            A(j-2,3)=0;
            A(j-2,4)=0;
        elseif NumGcp>=3
            A(j-3,1)=0;
            A(j-3,2)=0;
            A(j-3,3)=0;
            A(j-2,4)=0;
            A(j-2,5)=0;
            A(j-2,6)=0;
        end
        
        B(j-3,3*i-2)=(RPCsLR(1,1).xScale/RPCsLR(1,1).XScale)*((dPdX(1)*PP(2)-dPdX(2)*PP(1))/PP(2)/PP(2));
        B(j-3,3*i-1)=(RPCsLR(1,1).xScale/RPCsLR(1,1).YScale)*((dPdY(1)*PP(2)-dPdY(2)*PP(1))/PP(2)/PP(2));
        B(j-3,3*i)=(RPCsLR(1,1).xScale/RPCsLR(1,1).ZScale)*((dPdZ(1)*PP(2)-dPdZ(2)*PP(1))/PP(2)/PP(2));
        B(j-2,3*i-2)=(RPCsLR(1,1).yScale/RPCsLR(1,1).XScale)*((dPdX(3)*PP(4)-dPdX(4)*PP(3))/PP(4)/PP(4));
        B(j-2,3*i-1)=(RPCsLR(1,1).yScale/RPCsLR(1,1).YScale)*((dPdY(3)*PP(4)-dPdY(4)*PP(3))/PP(4)/PP(4));
        B(j-2,3*i)=(RPCsLR(1,1).yScale/RPCsLR(1,1).ZScale)*((dPdZ(3)*PP(4)-dPdZ(4)*PP(3))/PP(4)/PP(4));
        
        L(j-3,1)=duv(i,1);
        L(j-2,1)=duv(i,2);
        %-------------左像end-------------------%
        
        %-------------右像-------------------%
        Xt=(XYZ(i,1)-RPCsLR(2,1).XOffset)/RPCsLR(2,1).XScale;
        Yt=(XYZ(i,2)-RPCsLR(2,1).YOffset)/RPCsLR(2,1).YScale;
        Zt=(XYZ(i,3)-RPCsLR(2,1).ZOffset)/RPCsLR(2,1).ZScale;
        for k=1:4
            PP(k)=NCHLCalP(RPCsLR(2,1).P{k},Xt,Yt,Zt);
            dPdX(k)=CaldPdX(RPCsLR(2,1).P{k},Xt,Yt,Zt);   %P1、P2、P3、P4对X求偏导
            dPdY(k)=CaldPdY(RPCsLR(2,1).P{k},Xt,Yt,Zt);   %P1、P2、P3、P4对Y求偏导
            dPdZ(k)=CaldPdZ(RPCsLR(2,1).P{k},Xt,Yt,Zt);	%P1、P2、P3、P4对Z求偏导
        end
        
        if NumGcp==0||NumGcp==1   %只有平移
            A(j-1,1)=1;
            A(j,2)=1;
        elseif NumGcp==2  %平移和缩放
            A(j-1,1)=1;
            A(j-1,2)=uv2(i,1);
            A(j,3)=1;
            A(j,4)=uv2(i,1);
        elseif NumGcp>=3
            A(j-1,1)=1;
            A(j-1,2)=uv2(i,1);
            A(j-1,3)=uv2(i,2);
            A(j,4)=1;
            A(j,5)=uv2(i,1);
            A(j,6)=uv2(i,2);
        end
        
        B(j-1,3*i-2)=(RPCsLR(2,1).xScale/RPCsLR(2,1).XScale)*((dPdX(1)*PP(2)-dPdX(2)*PP(1))/PP(2)/PP(2));
        B(j-1,3*i-1)=(RPCsLR(2,1).xScale/RPCsLR(2,1).YScale)*((dPdY(1)*PP(2)-dPdY(2)*PP(1))/PP(2)/PP(2));
        B(j-1,3*i)=(RPCsLR(2,1).xScale/RPCsLR(2,1).ZScale)*((dPdZ(1)*PP(2)-dPdZ(2)*PP(1))/PP(2)/PP(2));
        B(j,3*i-2)=(RPCsLR(2,1).yScale/RPCsLR(2,1).XScale)*((dPdX(3)*PP(4)-dPdX(4)*PP(3))/PP(4)/PP(4));
        B(j,3*i-1)=(RPCsLR(2,1).yScale/RPCsLR(2,1).YScale)*((dPdY(3)*PP(4)-dPdY(4)*PP(3))/PP(4)/PP(4));
        B(j,3*i)=(RPCsLR(2,1).yScale/RPCsLR(2,1).ZScale)*((dPdZ(3)*PP(4)-dPdZ(4)*PP(3))/PP(4)/PP(4));
        L(j-1,1)=duv(i,3);
        L(j,1)=duv(i,4);
    end

    
    Nx=A'*A-A'*B*(inv(B'*B))*(B'*A);
    Rx=A'*L-A'*B*(inv(B'*B))*B'*L;
    Nt=B'*B-(B'*A)*(inv(A'*A))*(A'*B);
    Rt=B'*L-(B'*A)*(inv(A'*A))*A'*L;
    



%%%%% - Proposed BL1L2N Block
addpath('HJ\')
     N=[A B];
    TX = (N'*N+0.01*eye(size(N,2)))\N'*L;
    if count == 0
        D = Construction_Toeplitz([TX;TX], [1,-1]);
    end
    [TX] = balanced_L2_L1d2(N, L, 1e-3, TX);%Proposed BL1L2N 

    Mx=TX(1:2*NumGcp);
    Mt=TX(2*NumGcp+1:end);
    MaxC=max(abs(Mx));  %不满足条件就求解的参数就舍去
    if MaxVale>MaxC
        MaxVale=MaxC;
    else
        break;
    end
%%%%%%%% - 韩杰改
    
    if NumGcp==0||NumGcp==1   %只有平移
        re0=re0+Mx(1,1);
        rf0=rf0+Mx(2,1);
    elseif NumGcp==2  %平移和缩放
        re0=re0+Mx(1,1);	 re1=re1+Mx(2,1);
        rf0=rf0+Mx(3,1);	 rf1=rf1+Mx(4,1);
    elseif NumGcp>=3
        re0=re0+Mx(1,1);	 re1=re1+Mx(2,1);	 re2=re2+Mx(3,1);
        rf0=rf0+Mx(4,1);	 rf1=rf1+Mx(5,1);   rf2=rf2+Mx(6,1);
    end
    
    
    j=0;
    for i=1:nIP
        XYZ(i,1)=XYZ(i,1)+Mt(j+1,1);
        XYZ(i,2)=XYZ(i,2)+Mt(j+2,1);
        XYZ(i,3)=XYZ(i,3)+Mt(j+3,1);
        j=j+3;
    end
    count=count+1; 
end
toc
% % % %-----利用检核点结果分析------% % % %
fprintf(fid,'******************************光束法平差后结果分析******************************\n');
fprintf(fid,'反投影差值：\n');      %反投影差值
fprintf(fid,'点号\t\t\t左片du\t\t\t\t\t左片dv\t\t\t\t\t右片du\t\t\t\t\t右片du\n');
j=1;
for i=1:nIP
    V(j,1)=(XYZ(i,1)-chkBLH(i,1))*3600;
    V(j,2)=(XYZ(i,2)-chkBLH(i,2))*3600;
    V(j,3)=XYZ(i,3)-chkBLH(i,3);
    dduv(j,1)=duv(i,1);
    dduv(j,2)=duv(i,2);
    dduv(j,3)=duv(i,3);
    dduv(j,4)=duv(i,4);
    fprintf(fid,'%d   %20.8f   %20.8f   %20.8f   %20.8f\n',i,dduv(j,1),dduv(j,2),dduv(j,3),dduv(j,4));
    j=j+1;
end
LAbsRowmax=max(abs(dduv(:,1)));    LAbsColmax=max(abs(dduv(:,2)));
RAbsRowmax=max(abs(dduv(:,3)));    RAbsColmax=max(abs(dduv(:,4)));
fprintf(fid,'像方最大误差\n');
fprintf(fid,'%20.8f   %20.8f   %20.8f   %20.8f\n',LAbsRowmax,LAbsColmax,RAbsRowmax,RAbsColmax);

LRowRMS=sqrt((sum(dduv(:,1).*dduv(:,1)))/nIP);    LColRMS=sqrt((sum(dduv(:,2).*dduv(:,2)))/nIP);
RRowRMS=sqrt((sum(dduv(:,3).*dduv(:,3)))/nIP);   RColRMS=sqrt((sum(dduv(:,4).*dduv(:,4)))/nIP);
fprintf(fid,'像方中误差\n');
fprintf(fid,'%20.8f   %20.8f   %20.8f   %20.8f\n',LRowRMS,LColRMS,RRowRMS,RColRMS);
fprintf(fid,'----------------------------------------------------------------------------------\n');
fprintf(fid,'地面点坐标：\n');   %待定点坐标
fprintf(fid,'点号\t\t\tB\t\t\t\t\t\tL\t\t\t\t\t\tH\n');
for i=1:nIP
     fprintf(fid,'%d   %20.8f   %20.8f   %20.8f\n',i,XYZ(i,1),XYZ(i,2),XYZ(i,3));
end

fprintf(fid,'----------------------------------------------------------------------------------\n');
fprintf(fid,'点号\t\t\tdB\t\t\t\t\t\tdL\t\t\t\t\t\tdH\n');  %检核点比较
for i=1:length(chkBLH)
     fprintf(fid,'%d   %20.8f   %20.8f   %20.8f\n',i,V(i,1),V(i,2),V(i,3));
end

Bmax=max(abs(V(:,1)));    Lmax=max(abs(V(:,2)));
Hmax=max(abs(V(:,3)));
fprintf(fid,'物方最大误差\n');
fprintf(fid,'%20.8f   %20.8f   %20.8f\n',Bmax,Lmax,Hmax);

BRMS=sqrt((sum(V(:,1).*V(:,1)))/nIP);    LRMS=sqrt((sum(V(:,2).*V(:,2)))/nIP);
HRMS=sqrt((sum(V(:,3).*V(:,3)))/nIP);  
fprintf(fid,'物方中误差\n');
fprintf(fid,'%20.8f   %20.8f   %20.8f\n',BRMS,LRMS,HRMS);
fprintf(fid,'----------------------------------------------------------------------------------\n');

%XYZ相关信息
fprintf(fid,'\n');  fprintf(fid,'\n');
fprintf(fid,'地面点坐标XYZ：\n');
for i=1:nIP  %ctpXYZ:控制点XYZ
    temp_XYZ=BLH2XYZ(chkBLH(i,1),chkBLH(i,2),chkBLH(i,3));
    ctpXYZ(i,1)=temp_XYZ(1,1);
    ctpXYZ(i,2)=temp_XYZ(1,2);
    ctpXYZ(i,3)=temp_XYZ(1,3);
    totalXYZ(i,1)=ctpXYZ(i,1);  % totalXYZ: 控制点XYZ+待定点XYZ  
    totalXYZ(i,2)=ctpXYZ(i,2); 
    totalXYZ(i,3)=ctpXYZ(i,3); 
end

fprintf(fid,'点号\t\t\tdX\t\t\t\t\t\tdY\t\t\t\t\t\tdZ\n');  %检核点比较
for i=1:nIP
    temp_XYZ=BLH2XYZ(XYZ(i,1),XYZ(i,2),XYZ(i,3));
    calXYZ(i,1)=temp_XYZ(1,1);
    calXYZ(i,2)=temp_XYZ(1,2);
    calXYZ(i,3)=temp_XYZ(1,3);
    dXYZ(i,1)=totalXYZ(i,1)-calXYZ(i,1);
    dXYZ(i,2)=totalXYZ(i,2)-calXYZ(i,2);
    dXYZ(i,3)=totalXYZ(i,3)-calXYZ(i,3);
    fprintf(fid,'%d   %20.8f   %20.8f   %20.8f\n',i,dXYZ(i,1),dXYZ(i,2),dXYZ(i,3));
end

Xmax=max(abs(dXYZ(:,1)));    Ymax=max(abs(dXYZ(:,2)));
Zmax=max(abs(dXYZ(:,3)));

Xmean=mean(dXYZ(:,1));
Ymean=mean(dXYZ(:,2));
Zmean=mean(dXYZ(:,3));

fprintf(fid,'物方最大误差\n');
fprintf(fid,'%20.8f   %20.8f   %20.8f\n',Xmax,Ymax,Zmax);

XRMS=sqrt((sum(dXYZ(:,1).*dXYZ(:,1)))/nIP);    YRMS=sqrt((sum(dXYZ(:,2).*dXYZ(:,2)))/nIP);
ZRMS=sqrt((sum(dXYZ(:,3).*dXYZ(:,3)))/nIP);  
fprintf(fid,'物方中误差\n');
fprintf(fid,'%20.8f   %20.8f   %20.8f\n',XRMS,YRMS,ZRMS);

fprintf(fid,'物方误差均值\n');
fprintf(fid,'%20.8f   %20.8f   %20.8f\n',Xmean,Ymean,Zmean);

fprintf(fid,'******************************光束法平差后结果分析end******************************\n');
fclose(fid);
msgbox('结算完毕','提示');

%************************残差图显示***********************************%
figure('Name','前方交会');
hold on
xlabel('X'); ylabel('Y'); title('前方交会');

for i=1:nIP%绘制检核点
    plot(ctpXYZ(i,1),ctpXYZ(i,2),'*');
end
scalexy=50; scaleh=100;
quiver(ctpXYZ(:,1),ctpXYZ(:,2),resect_dXYZ(1:end,1)*scalexy,resect_dXYZ(1:end,2)*scalexy,0);
dX1=zeros(nIP,1);
quiver(ctpXYZ(:,1),ctpXYZ(:,2),dX1,resect_dXYZ(1:end,3)*scaleh,0);

scalexy=50; scaleh=100;
figure('Name','光束法平差');
hold on
xlabel('X'); ylabel('Y'); title('光束法平差');
for i=1:nIP%绘制控制点
    plot(ctpXYZ(i,1),ctpXYZ(i,2),'*');
end
quiver(ctpXYZ(:,1),ctpXYZ(:,2),dXYZ(1:end,1)*scalexy,dXYZ(1:end,2)*scalexy,0);
dX1=zeros(nIP,1);
quiver(ctpXYZ(:,1),ctpXYZ(:,2),dX1,dXYZ(1:end,3)*scaleh,0);


%%%%--------------------一块显示------------

%前方交会和区域网平差一块显示
axis ij
subplot(1,2,1)
hold on
box on
xlabel('X(m)'); ylabel('Y(m)'); title('{Before Adjustment(\times200)}');
for i=1:nIP%绘制检核点
    plot(ctpXYZ(i,1),ctpXYZ(i,2),'*k',LineWidth=1.5);
end
scalexy=200; scaleh=200;
quiver(ctpXYZ(:,1),ctpXYZ(:,2),resect_dXYZ(:,1)*scalexy,resect_dXYZ(:,2)*scalexy,0,LineWidth=1.5,Color='b');
dX1=zeros(nIP,1);
quiver(ctpXYZ(:,1),ctpXYZ(:,2),dX1,resect_dXYZ(:,3)*scaleh,0,LineWidth=1.5,Color='r');

subplot(1,2,2)
hold on
box on
xlabel('X(m)'); ylabel('Y(m)'); title('{After Adjustment(\times200)}');
for i=1:nIP%绘制检核点
    plot(ctpXYZ(i,1),ctpXYZ(i,2),'*k',LineWidth=1.5);
end
scalexy=200; scaleh=200;
quiver(ctpXYZ(:,1),ctpXYZ(:,2),dXYZ(:,1)*scalexy,dXYZ(:,2)*scalexy,0,LineWidth=1.5,Color='b');
dX1=zeros(nIP,1);
quiver(ctpXYZ(:,1),ctpXYZ(:,2),dX1,dXYZ(:,3)*scaleh,0,LineWidth=1.5,Color='r');
legend('XY-RMSE','Z-RMSE');


figure
hold on
box on
axis equal 
set(gca, 'FontName','times', 'FontSize',18, 'FontWeight','bold','LineWidth',1.5);
xlabel('{\it{u}}{(pixels\times 1000)}','fontname','times','fontsize',25,'FontWeight','bold'); 
ylabel('{\it{v}}{(pixels\times 1000)}','fontname','times','fontsize',25,'FontWeight','bold'); 

for i=1:nIP%绘制检核点
    plot(ctpXYZ(i,1),ctpXYZ(i,2),'*b','LineWidth',1);
end
scalexy=1000; scaleh=1000;
quiver(ctpXYZ(:,1),ctpXYZ(:,2),dduv(:,3)*scalexy,dduv(:,4)*scalexy,0,'b','LineWidth',1);

toc
end



%---------被调用函数及其定义函数定义---------------%

function uv=BackProj_RFM(img,BLH_Back)
%单景影像多个点反投影求解像点坐标
nRow=length(BLH_Back);
for i=1:nRow
   Xt=BLH_Back(i,1);
   Yt=BLH_Back(i,2);
   Zt=BLH_Back(i,3);
   X=(Xt-img.XOffset)/img.XScale;
   Y=(Yt-img.YOffset)/img.YScale;
   Z=(Zt-img.ZOffset)/img.ZScale;
   for j=1:4
       P(j)=NCHLCalP(img.P{j},X,Y,Z);
   end
  uv(i,1)=P(1)/P(2)*img.xScale+img.xOffset;
  uv(i,2)=P(3)/P(4)*img.yScale+img.yOffset;
end
end


function val=NCHLCalP(P,X,Y,Z)
%计算P(i)的值
val=P(1)+P(2)*Y+P(3)*X+P(4)*Z+P(5)*Y*X...       
			+P(6)*Y*Z+P(7)*X*Z+P(8)*Y*Y+P(9)*X*X...
			+P(10)*Z*Z+P(11)*Y*X*Z+P(12)*Y*Y*X+P(13)*Y*Y*Z...
			+P(14)*Y*X*X+P(15)*X*X*Z+P(16)*Y*Z*Z+P(17)*X*Z*Z...
			+P(18)*Y*Y*Y+P(19)*X*X*X+P(20)*Z*Z*Z;
end

function val=CaldPdX(P,X,Y,Z)
val=P(3)+P(5)*Y+P(7)*Z+2*P(9)*X+P(11)*Y*Z...
			+P(12)*Y*Y+2*P(14)*Y*X+2*P(15)*X*Z...
			+P(17)*Z*Z+3*P(19)*X*X;
end
function val=CaldPdY(P,X,Y,Z)
val=P(2)+P(5)*X+P(6)*Z+2*P(8)*Y+P(11)*X*Z...
			+2*P(12)*Y*X+2*P(13)*Y*Z+P(14)*X*X...
			+P(16)*Z*Z+3*P(18)*Y*Y;
end
function val=CaldPdZ(P,X,Y,Z)
val=P(4)+P(6)*Y+P(7)*X+2*P(10)*Z+P(11)*X*Y...
			+P(13)*Y*Y+P(15)*X*X+2*P(16)*Y*Z...
			+2*P(17)*X*Z+3*P(20)*Z*Z;
end

function XYZ=BLH2XYZ(B,L,H)
%大地坐标转换到空间直角坐标,参数B,L单位为度
a=6378137.0;  e2=0.0066943799803493275474544732430902;
B=d2rad(B);
L=d2rad(L);
sinB=sin(B);
cosB=cos(B);
sinL=sin(L);
cosL=cos(L);
N=a/sqrt(1-e2*sinB*sinB);
XYZ(1,1)=(N+H)*cosB*cosL;
XYZ(1,2)=(N+H)*cosB*sinL;
XYZ(1,3)=(N*(1-e2)+H)*sinB;
end
function radian=d2rad(du)
radian=du*pi/180;
end