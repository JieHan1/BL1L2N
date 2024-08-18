function BLH=Test222E2_BA_RPCBA_AFFINE2

tic
%%%%%%%%%%************************************************%%%%%%%%%%%

%���ڻ�׼Ӱ���������ƽ�
%1.����һ��Ӱ����Ϊ��׼��Ҳ������ǰ���������Ϊ��ֵ�����񷽸�����
%2.���õ�һ��ԭʼRPC���м��㣬���������ƽ������
%3.���õ�һ��Ӱ��RPC�����񷽷���任��ֵ���м��㣬���������ƽ����
%4.�������������񷽲в����ʽ

%%%%%%%%%%************************************************%%%%%%%%%%%

set(0,'defaultfigurecolor','w');
%%%%------���ݶ�ȡ-------%%%%
% %��ȡ��ƬRPC����
close all;
% [filename,pathname]=uigetfile('.txt','ѡ����Ӧ��RPC�ļ�');
% fnopen=fopen(strcat(pathname,filename),'rt');
% if(fnopen==-1)
%     msgbox('Input file or path is not correct',waring','warm');
% end
% RPCsLR(1,1) = readrpc([pathname,filename]);
RPCsLR(1,1) = readrpc('zy301a_fwd_381_6138_rpc.txt');

%��ȡ��ƬRPC����
% [filename,pathname]=uigetfile('.txt','ѡ����Ӧ��RPC�ļ�');
% fnopen=fopen(strcat(pathname,filename),'rt');
% if(fnopen==-1)
%     msgbox('Input file or path is not correct',waring','warm');
% end
% RPCsLR(2,1) = readrpc([pathname,filename]);
RPCsLR(2,1) = readrpc('zy301a_bwd_381_6138_rpc.txt');

%��ȡ�����������
% [filename,pathname]=uigetfile('.txt','ѡ��ͬ��������');
% fnopen=fopen(strcat(pathname,filename),'rt');
% if(fnopen==-1)
%     msgbox('Input file or path is not correct',waring','warm');
% end

fnopen=fopen('ƽ��7������\��ͬ����.txt','rt');
nIP=fscanf(fnopen,'%d',1); %������
for i=1:nIP
    fscanf(fnopen,'%d',1);
    Left_rc(i,1)=fscanf(fnopen,'%f',1);
    Left_rc(i,2)=fscanf(fnopen,'%f',1);
    Right_rc(i,1)=fscanf(fnopen,'%f',1);
    Right_rc(i,2)=fscanf(fnopen,'%f',1);
    ImgUV(i,1)=Left_rc(i,1);   ImgUV(i,2)=Left_rc(i,2);
    ImgUV(i,3)=Right_rc(i,1);  ImgUV(i,4)=Right_rc(i,2);
end

%%%%��ȡ��˵�
% [filename,pathname]=uigetfile('.txt','ѡ���˵�����');
% fnopen=fopen(strcat(pathname,filename),'rt');
% if(fnopen==-1)
%     msgbox('Input file or path is not correct',waring','warm');
% end
fnopen=fopen('ƽ��7������\��˵�.txt','rt');
nChk=fscanf(fnopen,'%d',1);
for i=1:nChk
    fscanf(fnopen,'%d',1);
    chkBLH(i,1)=fscanf(fnopen,'%f',1);
    chkBLH(i,2)=fscanf(fnopen,'%f',1);
    chkBLH(i,3)=fscanf(fnopen,'%f',1);
end
fclose('all');
%%%%------���ݶ�ȡ���-------%%%%

%%ֱ������RPC����������겢��ʵ��ֵ�Ƚ�:�������Ƶ�+������%
for i=1:nIP
     BBB=RFM_ForwardResection(RPCsLR(1,1),RPCsLR(2,1),Left_rc(i,:),Right_rc(i,:));
     XYZ(i,1)=BBB(1,1);  XYZ(i,2)=BBB(1,2);   XYZ(i,3)=BBB(1,3);
end



fid=fopen('Result_BL1L2N.txt','w');       %��������ļ�




if(fid==-1)
    msgbox('Input file or path is not correct','waring','warm');
end
fprintf(fid,'******************************������_ֱ�Ӷ���������******************************\n');
fprintf(fid,'ǰ��������������XYZ��\n');
fprintf(fid,'���\t\t\tX\t\t\t\t\t\tY\t\t\t\t\t\tZ\n'); 
for i=1:nIP  %chkXYZ����˵�XYZ
    temp_XYZ=BLH2XYZ(XYZ(i,1),XYZ(i,2),XYZ(i,3));
    resectXYZ(i,1)=temp_XYZ(1,1);  % totalXYZ: ���Ƶ�XYZ+������XYZ 
    resectXYZ(i,2)=temp_XYZ(1,2);
    resectXYZ(i,3)=temp_XYZ(1,3);
    fprintf(fid,'%d   %20.8f   %20.8f   %20.8f\n',i,resectXYZ(i,1),resectXYZ(i,2),resectXYZ(i,3));
end

fprintf(fid,'����������XYZ��\n');
fprintf(fid,'���\t\t\tX\t\t\t\t\t\tY\t\t\t\t\t\tZ\n'); 

for i=1:nChk  %chkXYZ����˵�XYZ
    temp_XYZ=BLH2XYZ(chkBLH(i,1),chkBLH(i,2),chkBLH(i,3));
    totalXYZ(i,1)=temp_XYZ(1,1);  % totalXYZ: ���Ƶ�XYZ+������XYZ 
    totalXYZ(i,2)=temp_XYZ(1,2);
    totalXYZ(i,3)=temp_XYZ(1,3);
    fprintf(fid,'%d   %20.8f   %20.8f   %20.8f\n',i,totalXYZ(i,1),totalXYZ(i,2),totalXYZ(i,3));
end

fprintf(fid,'���Ƶ�+��˵�������ǰ���������Ƚϣ�\n');
fprintf(fid,'���\t\t\tdX\t\t\t\t\t\tdY\t\t\t\t\t\tdZ\n');  %��˵�Ƚ�
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

fprintf(fid,'�﷽������\n');
fprintf(fid,'%20.8f   %20.8f   %20.8f\n',Xmax,Ymax,Zmax);

XRMS=sqrt((sum(resect_dXYZ(:,1).*resect_dXYZ(:,1)))/nChk);    YRMS=sqrt((sum(resect_dXYZ(:,2).*resect_dXYZ(:,2)))/nChk);
ZRMS=sqrt((sum(resect_dXYZ(:,3).*resect_dXYZ(:,3)))/nChk);  
fprintf(fid,'�﷽�����\n');
fprintf(fid,'%20.8f   %20.8f   %20.8f\n',XRMS,YRMS,ZRMS);

fprintf(fid,'�﷽����ֵ\n');
fprintf(fid,'%20.8f   %20.8f   %20.8f\n',Xmean,Ymean,Zmean);
fprintf(fid,'******************************ֱ�Ӷ���������end******************************\n');
fprintf(fid,'\n\n');
%��������������������������������������ֱ�Ӷ�����������ϡ���������������������������������������������
RP=zeros(1,6);
re0=RP(1,1);re1=RP(1,2);re2=RP(1,3);rf0=RP(1,4);rf1=RP(1,5);rf2=RP(1,6);
%le0=RP(1,1);le1=RP(1,2);le2=RP(1,3);lf0=RP(1,4);lf1=RP(1,5);lf2=RP(1,6);
le0=6.57488519;le1=-0.00014762;le2=-0.00003674;lf0=22.95222291;lf1=0.00013474;lf2=0.00008032;  %��׼Ӱ���Խ��з���任

%------�������任ϵ����ֵ��6����end-----%	

%------������ƽ��������--------%
count=0;
MaxVale=1000;  LimitVal=0.0000001;
NumGcp=2;  %ȷ���񷽱任ģ��
figure; %L����
tic
while count<50&&MaxVale>LimitVal
     if NumGcp==0||NumGcp==1   %ֻ��ƽ��
        A=zeros(2*nIP,2); %ÿ��Ӱ��������ƽ�Ʋ���
    elseif NumGcp==2  %ƽ�ƺ�����
        A=zeros(2*nIP,4);%ÿ��Ӱ��������ƽ�Ʋ������������Ų���
    elseif NumGcp>=3
        A=zeros(2*nIP,6);  %����任
     end
    B=zeros(4*nIP,3*nChk);
    L=zeros(4*nIP,1);
    Mx=zeros(6,1);      %ϵͳ����������
    Mt=zeros(3*nChk,1);    %���������������
    uv1=BackProj_RFM(RPCsLR(1,1),XYZ);  %�ÿ��Ƶ㷴���������
    uv2=BackProj_RFM(RPCsLR(2,1),XYZ);  %�ÿ��Ƶ㷴���������
    for i=1:nIP  %uvF����ͶӰ����任������
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
        %-------------����-------------------%
        Xt=(XYZ(i,1)-RPCsLR(1,1).XOffset)/RPCsLR(1,1).XScale;
        Yt=(XYZ(i,2)-RPCsLR(1,1).YOffset)/RPCsLR(1,1).YScale;
        Zt=(XYZ(i,3)-RPCsLR(1,1).ZOffset)/RPCsLR(1,1).ZScale;
        for k=1:4
            PP(k)=NCHLCalP(RPCsLR(1,1).P{k},Xt,Yt,Zt);
            dPdX(k)=CaldPdX(RPCsLR(1,1).P{k},Xt,Yt,Zt);   %P1��P2��P3��P4��X��ƫ��
            dPdY(k)=CaldPdY(RPCsLR(1,1).P{k},Xt,Yt,Zt);   %P1��P2��P3��P4��Y��ƫ��
            dPdZ(k)=CaldPdZ(RPCsLR(1,1).P{k},Xt,Yt,Zt);	%P1��P2��P3��P4��Z��ƫ��
        end
        if NumGcp==0||NumGcp==1   %ֻ��ƽ��
            A(j-3,1)=0;
            A(j-2,2)=0;
        elseif NumGcp==2  %ƽ�ƺ�����
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
        %-------------����end-------------------%
        
        %-------------����-------------------%
        Xt=(XYZ(i,1)-RPCsLR(2,1).XOffset)/RPCsLR(2,1).XScale;
        Yt=(XYZ(i,2)-RPCsLR(2,1).YOffset)/RPCsLR(2,1).YScale;
        Zt=(XYZ(i,3)-RPCsLR(2,1).ZOffset)/RPCsLR(2,1).ZScale;
        for k=1:4
            PP(k)=NCHLCalP(RPCsLR(2,1).P{k},Xt,Yt,Zt);
            dPdX(k)=CaldPdX(RPCsLR(2,1).P{k},Xt,Yt,Zt);   %P1��P2��P3��P4��X��ƫ��
            dPdY(k)=CaldPdY(RPCsLR(2,1).P{k},Xt,Yt,Zt);   %P1��P2��P3��P4��Y��ƫ��
            dPdZ(k)=CaldPdZ(RPCsLR(2,1).P{k},Xt,Yt,Zt);	%P1��P2��P3��P4��Z��ƫ��
        end
        
        if NumGcp==0||NumGcp==1   %ֻ��ƽ��
            A(j-1,1)=1;
            A(j,2)=1;
        elseif NumGcp==2  %ƽ�ƺ�����
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
    MaxC=max(abs(Mx));  %���������������Ĳ�������ȥ
    if MaxVale>MaxC
        MaxVale=MaxC;
    else
        break;
    end
%%%%%%%% - ���ܸ�
    
    if NumGcp==0||NumGcp==1   %ֻ��ƽ��
        re0=re0+Mx(1,1);
        rf0=rf0+Mx(2,1);
    elseif NumGcp==2  %ƽ�ƺ�����
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
% % % %-----���ü�˵�������------% % % %
fprintf(fid,'******************************������ƽ���������******************************\n');
fprintf(fid,'��ͶӰ��ֵ��\n');      %��ͶӰ��ֵ
fprintf(fid,'���\t\t\t��Ƭdu\t\t\t\t\t��Ƭdv\t\t\t\t\t��Ƭdu\t\t\t\t\t��Ƭdu\n');
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
fprintf(fid,'��������\n');
fprintf(fid,'%20.8f   %20.8f   %20.8f   %20.8f\n',LAbsRowmax,LAbsColmax,RAbsRowmax,RAbsColmax);

LRowRMS=sqrt((sum(dduv(:,1).*dduv(:,1)))/nIP);    LColRMS=sqrt((sum(dduv(:,2).*dduv(:,2)))/nIP);
RRowRMS=sqrt((sum(dduv(:,3).*dduv(:,3)))/nIP);   RColRMS=sqrt((sum(dduv(:,4).*dduv(:,4)))/nIP);
fprintf(fid,'�������\n');
fprintf(fid,'%20.8f   %20.8f   %20.8f   %20.8f\n',LRowRMS,LColRMS,RRowRMS,RColRMS);
fprintf(fid,'----------------------------------------------------------------------------------\n');
fprintf(fid,'��������꣺\n');   %����������
fprintf(fid,'���\t\t\tB\t\t\t\t\t\tL\t\t\t\t\t\tH\n');
for i=1:nIP
     fprintf(fid,'%d   %20.8f   %20.8f   %20.8f\n',i,XYZ(i,1),XYZ(i,2),XYZ(i,3));
end

fprintf(fid,'----------------------------------------------------------------------------------\n');
fprintf(fid,'���\t\t\tdB\t\t\t\t\t\tdL\t\t\t\t\t\tdH\n');  %��˵�Ƚ�
for i=1:length(chkBLH)
     fprintf(fid,'%d   %20.8f   %20.8f   %20.8f\n',i,V(i,1),V(i,2),V(i,3));
end

Bmax=max(abs(V(:,1)));    Lmax=max(abs(V(:,2)));
Hmax=max(abs(V(:,3)));
fprintf(fid,'�﷽������\n');
fprintf(fid,'%20.8f   %20.8f   %20.8f\n',Bmax,Lmax,Hmax);

BRMS=sqrt((sum(V(:,1).*V(:,1)))/nIP);    LRMS=sqrt((sum(V(:,2).*V(:,2)))/nIP);
HRMS=sqrt((sum(V(:,3).*V(:,3)))/nIP);  
fprintf(fid,'�﷽�����\n');
fprintf(fid,'%20.8f   %20.8f   %20.8f\n',BRMS,LRMS,HRMS);
fprintf(fid,'----------------------------------------------------------------------------------\n');

%XYZ�����Ϣ
fprintf(fid,'\n');  fprintf(fid,'\n');
fprintf(fid,'���������XYZ��\n');
for i=1:nIP  %ctpXYZ:���Ƶ�XYZ
    temp_XYZ=BLH2XYZ(chkBLH(i,1),chkBLH(i,2),chkBLH(i,3));
    ctpXYZ(i,1)=temp_XYZ(1,1);
    ctpXYZ(i,2)=temp_XYZ(1,2);
    ctpXYZ(i,3)=temp_XYZ(1,3);
    totalXYZ(i,1)=ctpXYZ(i,1);  % totalXYZ: ���Ƶ�XYZ+������XYZ  
    totalXYZ(i,2)=ctpXYZ(i,2); 
    totalXYZ(i,3)=ctpXYZ(i,3); 
end

fprintf(fid,'���\t\t\tdX\t\t\t\t\t\tdY\t\t\t\t\t\tdZ\n');  %��˵�Ƚ�
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

fprintf(fid,'�﷽������\n');
fprintf(fid,'%20.8f   %20.8f   %20.8f\n',Xmax,Ymax,Zmax);

XRMS=sqrt((sum(dXYZ(:,1).*dXYZ(:,1)))/nIP);    YRMS=sqrt((sum(dXYZ(:,2).*dXYZ(:,2)))/nIP);
ZRMS=sqrt((sum(dXYZ(:,3).*dXYZ(:,3)))/nIP);  
fprintf(fid,'�﷽�����\n');
fprintf(fid,'%20.8f   %20.8f   %20.8f\n',XRMS,YRMS,ZRMS);

fprintf(fid,'�﷽����ֵ\n');
fprintf(fid,'%20.8f   %20.8f   %20.8f\n',Xmean,Ymean,Zmean);

fprintf(fid,'******************************������ƽ���������end******************************\n');
fclose(fid);
msgbox('�������','��ʾ');

%************************�в�ͼ��ʾ***********************************%
figure('Name','ǰ������');
hold on
xlabel('X'); ylabel('Y'); title('ǰ������');

for i=1:nIP%���Ƽ�˵�
    plot(ctpXYZ(i,1),ctpXYZ(i,2),'*');
end
scalexy=50; scaleh=100;
quiver(ctpXYZ(:,1),ctpXYZ(:,2),resect_dXYZ(1:end,1)*scalexy,resect_dXYZ(1:end,2)*scalexy,0);
dX1=zeros(nIP,1);
quiver(ctpXYZ(:,1),ctpXYZ(:,2),dX1,resect_dXYZ(1:end,3)*scaleh,0);

scalexy=50; scaleh=100;
figure('Name','������ƽ��');
hold on
xlabel('X'); ylabel('Y'); title('������ƽ��');
for i=1:nIP%���ƿ��Ƶ�
    plot(ctpXYZ(i,1),ctpXYZ(i,2),'*');
end
quiver(ctpXYZ(:,1),ctpXYZ(:,2),dXYZ(1:end,1)*scalexy,dXYZ(1:end,2)*scalexy,0);
dX1=zeros(nIP,1);
quiver(ctpXYZ(:,1),ctpXYZ(:,2),dX1,dXYZ(1:end,3)*scaleh,0);


%%%%--------------------һ����ʾ------------

%ǰ�������������ƽ��һ����ʾ
axis ij
subplot(1,2,1)
hold on
box on
xlabel('X(m)'); ylabel('Y(m)'); title('{Before Adjustment(\times200)}');
for i=1:nIP%���Ƽ�˵�
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
for i=1:nIP%���Ƽ�˵�
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

for i=1:nIP%���Ƽ�˵�
    plot(ctpXYZ(i,1),ctpXYZ(i,2),'*b','LineWidth',1);
end
scalexy=1000; scaleh=1000;
quiver(ctpXYZ(:,1),ctpXYZ(:,2),dduv(:,3)*scalexy,dduv(:,4)*scalexy,0,'b','LineWidth',1);

toc
end



%---------�����ú������䶨�庯������---------------%

function uv=BackProj_RFM(img,BLH_Back)
%����Ӱ�����㷴ͶӰ����������
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
%����P(i)��ֵ
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
%�������ת�����ռ�ֱ������,����B,L��λΪ��
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