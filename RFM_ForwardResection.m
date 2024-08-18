function BLH=RFM_ForwardResection(img1,img2,Left_rc,Right_rc)
Xt=(img1.XOffset+img2.XOffset)/2;
Yt=(img1.YOffset+img2.YOffset)/2;
Zt=(img1.ZOffset+img2.ZOffset)/2;
u1=Left_rc(1,1);
v1=Left_rc(1,2);
u2=Right_rc(1,1);
v2=Right_rc(1,2);
count=0;
while 1
    %img1
    X=(Xt-img1.XOffset)/img1.XScale;
    Y=(Yt-img1.YOffset)/img1.YScale;
    Z=(Zt-img1.ZOffset)/img1.ZScale;
    for j=1:4
        P(j)=NCHLCalP(img1.P{j},X,Y,Z);
        dPdX(j)=CaldPdX(img1.P{j},X,Y,Z);
        dPdY(j)=CaldPdY(img1.P{j},X,Y,Z);
        dPdZ(j)=CaldPdZ(img1.P{j},X,Y,Z);
    end
    A(1,1)=(dPdX(1)*P(2)-dPdX(2)*P(1))/P(2)/P(2)/img1.XScale*img1.xScale;
    A(1,2)=(dPdY(1)*P(2)-dPdY(2)*P(1))/P(2)/P(2)/img1.YScale*img1.xScale;
    A(1,3)=(dPdZ(1)*P(2)-dPdZ(2)*P(1))/P(2)/P(2)/img1.ZScale*img1.xScale;
    A(2,1)=(dPdX(3)*P(4)-dPdX(4)*P(3))/P(4)/P(4)/img1.XScale*img1.yScale;
    A(2,2)=(dPdY(3)*P(4)-dPdY(4)*P(3))/P(4)/P(4)/img1.YScale*img1.yScale;
    A(2,3)=(dPdZ(3)*P(4)-dPdZ(4)*P(3))/P(4)/P(4)/img1.ZScale*img1.yScale;

    L(1,1)=u1-P(1)/P(2)*img1.xScale-img1.xOffset;
    L(2,1)=v1-P(3)/P(4)*img1.yScale-img1.yOffset;
		
    %img2
    X=(Xt-img2.XOffset)/img2.XScale;
    Y=(Yt-img2.YOffset)/img2.YScale;
    Z=(Zt-img2.ZOffset)/img2.ZScale;
    for j=1:4
        P(j)=NCHLCalP(img2.P{j},X,Y,Z);
        dPdX(j)=CaldPdX(img2.P{j},X,Y,Z);
        dPdY(j)=CaldPdY(img2.P{j},X,Y,Z);
        dPdZ(j)=CaldPdZ(img2.P{j},X,Y,Z);
    end
    A(3,1)=(dPdX(1)*P(2)-dPdX(2)*P(1))/P(2)/P(2)/img2.XScale*img2.xScale;
    A(3,2)=(dPdY(1)*P(2)-dPdY(2)*P(1))/P(2)/P(2)/img2.YScale*img2.xScale;
    A(3,3)=(dPdZ(1)*P(2)-dPdZ(2)*P(1))/P(2)/P(2)/img2.ZScale*img2.xScale;
    A(4,1)=(dPdX(3)*P(4)-dPdX(4)*P(3))/P(4)/P(4)/img2.XScale*img2.yScale;
    A(4,2)=(dPdY(3)*P(4)-dPdY(4)*P(3))/P(4)/P(4)/img2.YScale*img2.yScale;
    A(4,3)=(dPdZ(3)*P(4)-dPdZ(4)*P(3))/P(4)/P(4)/img2.ZScale*img2.yScale;

    L(3,1)=u2-P(1)/P(2)*img2.xScale-img2.xOffset;
    L(4,1)=v2-P(3)/P(4)*img2.yScale-img2.yOffset;
    
    re=inv(A'*A)*A'*L;
    Xt=Xt+re(1,1);
    Yt=Yt+re(2,1);
    Zt=Zt+re(3,1);
    count=count+1;
    if abs(re(1,1))<0.00000024&&abs(re(2,1))<0.00000024&&abs(re(3,1))<5||count>50
        break;
    end
end
BLH(1,1)=Xt;
BLH(1,2)=Yt;
BLH(1,3)=Zt;
end

function val=NCHLCalP(P,X,Y,Z)
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