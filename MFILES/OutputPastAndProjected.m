function output=OutputPastAndProjectedSL(TGdata,samps,targyears,params)

% OutputPastAndProjectedSL(TGdata,samps,targyears,params)
%
% 
% TGdata: data structure including (at a minimum) meantime, Y, dY, Ycv
% and sitenames
%
% samps: cell array of matrices of projection samples
%
% targyears: vector of years corresponding to columns of samps elements
%
% params: optional structure of parameters including
%    - smoothwin: smoothing window to use for datum shift
%    - datumyr: year for baseline datum
%    - unitscale: factor to multiple mm by to produce desired units (default: .1)
%    - unitlabel: label for units (default: cm)
%    - plevs: probability levels for historical data plot (default: [.83 .95])
%    - timelims: time limits for plot (default: [1900 2100])
%    - ylims: ylimits for plot (default: [-50 200])
%    - baseyear: base year for projections prior to datum shift (default: 2000)
%    - doplot: produce plot (default: 1)
%    - optimizemode: optimization mode for GP fit for smoothing (default: 1.1)
%    - scens: scenario labels for samps (default: {'RCP8.5','RCP6.0','RCP4.5','RCP2.6'})
%    - doscens: scenarios to include on plot (default: [1 4])
%    - tablefn: if non-empty, base filename to use for table output
%    - tableqlevs: quantile levels to include in projections table (default: [.05 .17 .5 .83 .95 .99 .995 .999])
%
% Last updated by Bob Kopp, 2017-08-17 00:02:30 -0400

defval('params',[]);
defval('smoothwin',[]);
defval('datumyr',[]);
defval('unitscale',.1);
defval('unitlabel','cm')
defval('plevs',[.83 .95]);
defval('timelims',[1900 2100]);
defval('ylims',[-50 300]);
defval('baseyear',2000);
defval('doplot',1);
defval('optimizemode',1);
defval('scens',{'RCP8.5','RCP6.0','RCP4.5','RCP2.6'});
defval('doscens',[1 4]);
defval('showsmooth',0);
defval('scencolors',[1 0 0; 0 0 1]);

defval('tablefn',[]);
defval('tableqlevs',[.05 .17 .5 .83 .95 .99 .995 .999]);

if isstruct(params)
    parseFields(params);
end

output=[];

if ~isempty(smoothwin)
    subt=find(abs(TGdata.meantime-baseyear)<=smoothwin/2);
    if length(subt)>.5*smoothwin
        Mdiff=eye(length(TGdata.meantime));
        Mdiff(:,subt)=Mdiff(:,subt)-1/length(subt);
        TGdata.Y=Mdiff*TGdata.Y;
        TGdata.Ycv=Mdiff*TGdata.Ycv*Mdiff';
        TGdata.dY=sqrt(diag(TGdata.Ycv));
    end 
end

datumoffset=0; datumoffsetsd=0;
if ~isempty(smoothwin)
    TGdata2s = GPSmoothTideGauges(TGdata,smoothwin,optimizemode,[1 1]);
    if ~isempty(datumyr)
        subt=find(TGdata2s.meantime==datumyr);
        if length(subt)==1
            Mdiff=eye(length(TGdata2s.meantime));
            Mdiff(:,subt)=Mdiff(:,subt)-1;
            TGdata2s.Y=Mdiff*TGdata2s.Y;
            TGdata2s.Ycv=Mdiff*TGdata2s.Ycv*Mdiff';
            TGdata2s.dY=sqrt(diag(TGdata2s.Ycv));
            subt2=find(TGdata2s.meantime==baseyear);
            datumoffset=-TGdata2s.Y(subt2);
            datumoffsetsd=TGdata2s.dY(subt2);
        end
    end
end

if doplot
    clear datas;
    hp = subplot(2,1,1);    
    datas.x=TGdata.meantime;
    datas.y=(TGdata.Y-datumoffset)*unitscale;
    datas.dy=sqrt(TGdata.dY.^2+datumoffsetsd^2)*unitscale*norminv(plevs(1));
    datas2=datas;  datas2.dy= datas2.dy*norminv(plevs(2));
    PlotWithShadedErrors(datas2,[.5 .5 .5 ],.9,'none',':');
    hl=PlotWithShadedErrors(datas,[.5 .5 .5],.7,'-','--');

    targyears=[baseyear targyears];
    for rrr=1:size(samps,2)
        samps{rrr}=[zeros(size(samps{rrr},1),1) samps{rrr}];
        seeds=linspace(0,1,size(samps{rrr},1)+2); seeds=seeds(2:end-1);
        seeds=norminv(seeds(randperm(length(seeds))));
        datumshiftsamp=datumoffset+seeds*datumoffsetsd;
        samps{rrr}=bsxfun(@minus,samps{rrr},datumshiftsamp(:));
    end

    clear projs projs2;
    for rrr=1:length(doscens)
        projs{rrr}.x=targyears';
        projs{rrr}.y=quantile(samps{doscens(rrr)},.5,1)'*unitscale;
        projs{rrr}.dy=abs(bsxfun(@minus,quantile(samps{doscens(rrr)},[1-plevs(1) plevs(1)],1)'*unitscale,projs{rrr}.y));
        projs2{rrr}=projs{rrr};
        projs2{rrr}.dy=abs(bsxfun(@minus,quantile(samps{doscens(rrr)},[1-plevs(2) plevs(2)],1)'*unitscale,projs{rrr}.y));
    end

    hl=PlotWithShadedErrors(projs,scencolors,.9,'-','--',[baseyear timelims(2)]);
    PlotWithShadedErrors(projs2,scencolors,'none','none',':',[baseyear timelims(2)]);
    hl=PlotWithShadedErrors(projs,scencolors,'none','none','--',[baseyear timelims(2)]);
    hl=PlotWithShadedErrors(datas,[.5 .5 .5],'none','-','none');
    if showsmooth
        if exist('TGdata2s','var')
            plot(TGdata2s.meantime,(TGdata2s.Y)*unitscale,'k','linewidth',.5)
        end
    end

    set(hp,'xlim',timelims,'ylim',ylims);
    box on; longticks(gca,2);
    ylabel(['Sea level (' unitlabel ')']);
end

if ~isempty(tablefn)
    fid=fopen([tablefn '_hist.tsv'],'w');
    if (datumoffset==0)&&(datumoffsetsd==0)
        fprintf(fid,[TGdata.sitenames{1} ' (' unitlabel ')\n']);
    else
        fprintf(fid,[TGdata.sitenames{1} ' (' unitlabel ' above %0.0f-yr window centered at %0.0f)\n'],[smoothwin datumyr]);
    end
    fprintf(fid,'Historical observations\n\n');
    fprintf(fid,'Year\tAnnual\t1s\tSmoothed\t1s\n');
    for ttt=1:length(TGdata.meantime)
        fprintf(fid,'%0.0f',TGdata.meantime(ttt));
        fprintf(fid,'\t%0.1f',(TGdata.Y(ttt)-datumoffset)*unitscale);
        fprintf(fid,'\t%0.1f',unitscale*sqrt(TGdata.dY(ttt)^2+datumoffsetsd^2));
        if exist('TGdata2s','var')
        selyr=find(TGdata2s.meantime==TGdata.meantime(ttt));
        if (length(selyr)==1)
            fprintf(fid,'\t%0.1f',unitscale*[TGdata2s.Y(selyr) TGdata2s.dY(selyr)]);
        else
            fprintf(fid,'\tNaN\tNaN');
        end
        else
            fprintf(fid,'\tNaN\tNaN');
        end
        fprintf(fid,'\n');
    end
    fclose(fid);

    for rrr=1:length(scens)
        fid=fopen([tablefn '_' scens{rrr} '.tsv'],'w');
        if (datumoffset==0)&&(datumoffsetsd==0)
            fprintf(fid,[TGdata.sitenames{1} ' (' unitlabel 'above %0.0f)\n'],baseyear);
        else
            fprintf(fid,[TGdata.sitenames{1} ' (' unitlabel ' above %0.0f-yr window centered at %0.0f)\n'],[smoothwin datumyr]);
        end
        fprintf(fid,[scens{rrr} '\n\n']);
        fprintf(fid,'Year');
        fprintf(fid,'\t%0.1f',tableqlevs*100);
        fprintf(fid,'\n');
        for ttt=1:length(targyears)
            fprintf(fid,'%0.0f',targyears(ttt));
            fprintf(fid,'\t%0.0f',quantile(samps{rrr}(:,ttt),tableqlevs,1)'*unitscale);
            fprintf(fid,'\n');
        end
        fclose(fid);
    end

end

function parseFields(params)

    flds=fieldnames(params);
    for qqq=1:length(flds)
        if flds{qqq} 
            assignin('caller',flds{qqq},params.(flds{qqq}));
        end
    end