function RPCs = readrpc(strfilepath)
% 2013.4.3
% Read ZY-3's rpc file(测绘卫星中心)
if strcmp(strfilepath(end-6:end),'rpc.txt')
    RPCs=[];
    fid=fopen(strfilepath);
    while(~feof(fid))
       strl=fgetl(fid);
       if(isempty(strl)),break;end
       CC=textscan(strl,'%s %s %s');
       if(isempty(CC{2})),continue;end
       term = strtrim(CC{1}{1});
       if size(term,2) > 15
           term = term(1:14);
       end
       val=str2double(CC{2}{1});
       switch term
           case('LINE_OFF:')
               RPCs.xOffset=val;
           case('SAMP_OFF:')
               RPCs.yOffset=val;
           case('LAT_OFF:')
               RPCs.XOffset=val;
           case('LONG_OFF:')
               RPCs.YOffset=val;
           case('HEIGHT_OFF:')
               RPCs.ZOffset=val;
           case('LINE_SCALE:')
               RPCs.xScale=val;
           case('SAMP_SCALE:')
               RPCs.yScale=val;
           case('LAT_SCALE:')
               RPCs.XScale=val;
           case('LONG_SCALE:')
               RPCs.YScale=val;
           case('HEIGHT_SCALE:')
               RPCs.ZScale=val;           
           case{'LINE_NUM_COEFF','LINE_DEN_COEFF','SAMP_NUM_COEFF','SAMP_DEN_COEFF'}
               P=zeros(1,20);
               P(1) = str2num(CC{2}{1});
               for i=2:20
                   strl=fgetl(fid);
                   CC=textscan(strl,'%s %s');
                   P(i)=str2double(CC{2}{1});
               end
               if(isequal(term,'LINE_NUM_COEFF'))
                   RPCs.P{1}=P;
               elseif(isequal(term,'LINE_DEN_COEFF'))
                   RPCs.P{2}=P;
               elseif(isequal(term,'SAMP_NUM_COEFF'))
                   RPCs.P{3}=P;
               elseif(isequal(term,'SAMP_DEN_COEFF'))
                   RPCs.P{4}=P;
               end
       end
    end
    fclose(fid);
end
% 资源卫星中心RPC
if strcmp(strfilepath(end-2:end),'rpb')
    RPCs=[];
    fid=fopen(strfilepath);
    for i=1:6
        strl=fgetl(fid);
    end
    while(~feof(fid))
        strl=fgetl(fid);
        if(isempty(strl)),break;end
        CC=textscan(strl,'%s %s %s');
        if(isempty(CC{3})),continue;end
        term = strtrim(CC{1}{1});
        val = CC{3}{1};
        val=str2double(val(1:end-1));
        switch term
           case('lineOffset')
               RPCs.xOffset=val;
           case('sampOffset')
               RPCs.yOffset=val;
           case('latOffset')
               RPCs.XOffset=val;
           case('longOffset')
               RPCs.YOffset=val;
           case('heightOffset')
               RPCs.ZOffset=val;
           case('lineScale')
               RPCs.xScale=val;
           case('sampScale')
               RPCs.yScale=val;
           case('latScale')
               RPCs.XScale=val;
           case('longScale')
               RPCs.YScale=val;
           case('heightScale')
               RPCs.ZScale=val;
           case{'lineNumCoef','lineDenCoef','sampNumCoef','sampDenCoef'}
               P=zeros(1,20);
               for i=1:20
                   strl=fgetl(fid);
                   CC=textscan(strl,'%s');
                   P(i)=str2double(CC{1}{1}(1:22));
               end
               if(isequal(term,'lineNumCoef'))
                   RPCs.P{1}=P;
               elseif(isequal(term,'lineDenCoef'))
                   RPCs.P{2}=P;
               elseif(isequal(term,'sampNumCoef'))
                   RPCs.P{3}=P;
               elseif(isequal(term,'sampDenCoef'))
                   RPCs.P{4}=P;
               end
        end
    end
    fclose(fid);
end


end