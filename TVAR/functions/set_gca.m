function [ax1 ax2 myytick] = set_gca(dates,C,col,ystep,freq,varargin)
%UNTITLED5 sets the axes so that they will print the the appropriate labels
%and limits
%   dates is a column of date numbers used to plot the graphs
%   C is the array containing all of the data that is plotted
%   col is entered as an array of the columns used in plotting the graph,
%   for example if columns 4,5, and 8 are plotted, then the input would be
%   {4,5,8}
%   ystep is the distance between labels on the y-axis. the limits are done
%   automatically, this is just the space to put between ticks
%   freq is the freq of data. 1 - for graphs printed with months on the
%   x-axis, 4 - for graphs printed with quarterly months (every 4th) on the
%   x-axis, 12 - for graphs printed with years on the x-axis
%   there are two optional arguments:
        % 'off' - turns off the right-hand side y-axis (eg if there are
        % multiple graphs on the page)
        % 2 - will print every other year instead of every year if space is
        % contrained

tmp=[kron((1801:2200)',ones(12,1)) kron(ones(400,1),(1:12)')];
tmp(:,3)=eomday(tmp(:,1),tmp(:,2));
tmpmi=datenum(tmp)+0.5;
tmpqi=datenum(tmp(3:3:end,:))+0.5;
tmpyi=datenum(tmp(12:12:end,:))+0.5;
tmpmm=datenum(tmp)+365/24+1.5;
tmpqm=datenum(tmp(3:3:end,:))+365/4;
tmpym=datenum(tmp(12:12:end,:))+365/2;
xtickmi=tmpmi((tmpmi>=min(dates))&(tmpmi<=max(dates)),:);
xtickqi=tmpqi((tmpqi>=min(dates))&(tmpqi<=max(dates)),:);
xtickyi=tmpyi((tmpyi>=min(dates))&(tmpyi<=max(dates)),:);
xtickmm=tmpmm((tmpmm>=min(dates))&(tmpmm<=max(dates)),:);
xtickqm=tmpqm((tmpqm>=min(dates))&(tmpqm<=max(dates)),:);
xtickym=tmpym((tmpym>=min(dates))&(tmpym<=max(dates)),:);

box on;
grid on;

xmin = dates(1);
xmax = dates(end);
ymax = 0;
ymin = 10000;


if isempty(varargin);
    right_axis = 'on';
    year_freq = 1;

elseif length(varargin)==1
    if strcmp(varargin{1},'off')
        right_axis = 'off';
        year_freq = 1;

    elseif varargin{1}==2
        year_freq = 2;
        right_axis = 'on';

    elseif strcmp(varargin{1},'3month')
        right_axis = 'on';
        year_freq = 1;

        
    end
else
    right_axis = varargin{1};
    year_freq = varargin{2};

end

for i = 1:length(col);
    max_1 = max(C{col{i}});
    if max_1>ymax
        ymax = max_1;
    end
end
for i = 1:length(col);
    min_1 = min(C{col{i}});
    if min_1<ymin
        ymin = min_1;
    end
end    



    ymin = floor(ymin/ystep)*ystep;
    ymax = ceil(ymax/ystep)*ystep;
    myytick = (ymin:ystep:ymax);

%
ax1=gca;

ax2=axes('Position',get(ax1,'Position'),...
   'YAxisLocation','right',...
   'Color','none'); 
set([ax1 ax2],'Xtick',[]);

if strcmp(right_axis,'on');
    set([ax1 ax2],'Ytick',myytick);
elseif strcmp(right_axis,'off')
    set(ax2,'Ytick',[]);
    set(ax1,'Ytick',myytick);
end

set([ax1 ax2],'xlim',[xmin xmax]);
set([ax1 ax2],'ylim',[ymin ymax]);
set([ax1 ax2],'TickLength',[0 0]);
myylim = get(gca,'yLim');
myxtick=myylim;
myxtick(2)=myxtick(1)+(myxtick(2)-myxtick(1))*0.01;

switch freq
   case 1
      set(ax2,'XTick',xtickmm,'TickLength',[0 0],'XtickLabel',[]);
      xlab = datestr(xtickmm,'mmm');
      for I=1:size(xlab,1)
          if isequal(xlab(I,:),'Jan')
             xlab(I,:)='Gen';
          elseif isequal(xlab(I,:),'May')
             xlab(I,:)='Mag';
          elseif isequal(xlab(I,:),'Jun')
             xlab(I,:)='Giu';
          elseif isequal(xlab(I,:),'Jul')
             xlab(I,:)='Lug';
          elseif isequal(xlab(I,:),'Aug')
             xlab(I,:)='Ago';
          elseif isequal(xlab(I,:),'Sep')
             xlab(I,:)='Set';
          elseif isequal(xlab(I,:),'Oct')
             xlab(I,:)='Ott';
          elseif isequal(xlab(I,:),'Dec')
             xlab(I,:)='Dic';
          end   
      end
      set(ax2,'XTickLabel',xlab);
      %datetick('x','mmm','keepticks','keeplimits');
      tmp1=union(xmin,union(xmax,xtickyi));
      tmp2=tmp1(2:end)-tmp1(1:end-1);
      tmp3=(tmp1(1:end-1)+tmp1(2:end))/2;
      tmp4=get(gca,'Position');
      tmp5=-0.04+tmp4(1)+tmp4(3)*(tmp3-xmin)./(xmax-xmin);
      tmp6=datevec(tmp3);
      tmp6=tmp6(:,1);

        for I=1:length(tmp3)
            if tmp2(I)>180
                annotation('textbox', [tmp5(I) tmp4(2)-0.05 0.075 0.05],...
                'String',num2str(tmp6(I)),...
                'LineStyle','none',...
                'HorizontalAlignment','center',...
                'FontSize',8,...
                'FontName','Times New Roman');
            end
        end
        line(ones(2,1)*xtickmi',myylim'*ones(1,length(xtickmi)),'Color','k','LineStyle',':');
%       set(ax2,'XTick',xtickym);
%       datetick('x','yyyy','keepticks','keeplimits');
%       line(ones(2,1)*xtickyi',myylim'*ones(1,length(xtickyi)),'Color','k');
%       line(ones(2,1)*xtickqi',myxtick'*ones(1,length(xtickqi)),'Color','k');
   case 4
      set(ax2,'XTick',xtickqm,'TickLength',[0 0],'XtickLabel',[]);
      %datetick('x','mmm','keepticks','keeplimits');
      xlab = datestr(xtickqm,'mmm');
      for I=1:size(xlab,1)
          if isequal(xlab(I,:),'Jan')
             xlab(I,:)='Gen';
          elseif isequal(xlab(I,:),'May')
             xlab(I,:)='Mag';
          elseif isequal(xlab(I,:),'Jun')
             xlab(I,:)='Giu';
          elseif isequal(xlab(I,:),'Jul')
             xlab(I,:)='Lug';
          elseif isequal(xlab(I,:),'Aug')
             xlab(I,:)='Ago';
          elseif isequal(xlab(I,:),'Sep')
             xlab(I,:)='Set';
          elseif isequal(xlab(I,:),'Oct')
             xlab(I,:)='Ott';
          elseif isequal(xlab(I,:),'Dec')
             xlab(I,:)='Dic';
          end   
      end
      set(ax2,'XTickLabel',xlab);
      %tmp1=union(xmin,union(xmax,xtickyi));
      %tmp2=tmp1(2:end)-tmp1(1:end-1);
      %tmp3=(tmp1(1:end-1)+tmp1(2:end))/2;
      %tmp4=get(gca,'Position');
      %tmp5=-0.04+tmp4(1)+tmp4(3)*(tmp3-xmin)./(xmax-xmin);
      %tmp6=datevec(tmp3);
      %tmp6=tmp6(:,1);

      %  for I=1:4:length(tmp3)
      %      if tmp2(I)>180
      %          annotation('textbox', [tmp5(I) tmp4(2)-0.05 0.075 0.05],...
      %          'String',num2str(tmp6(I)),...
      %          'LineStyle','none',...
      %          'HorizontalAlignment','center',...
      %          'FontSize',8,...
      %          'FontName','Times New Roman');
      %      end
      %  end
        line(ones(2,1)*xtickqi',myylim'*ones(1,length(xtickqi)),'Color','k','LineStyle',':');
        
        tmp1=union(xmin,union(xmax,xtickyi));
        tmp2=tmp1(2:end)-tmp1(1:end-1);
        tmp3=(tmp1(1:end-1)+tmp1(2:end))/2;
        tmp4=get(gca,'Position');
        w = tmp4(3)*tmp2./(xmax-xmin);
        tmp5=tmp4(1)+tmp4(3)*(tmp3-xmin)./(xmax-xmin)-w./2;    %%%%tmp5=-0.04+tmp4(1)+tmp4(3)*(tmp3-xmin)./(xmax-xmin);
        tmp6=datevec(tmp3);
        tmp6=tmp6(:,1);
        for I=1:length(tmp3)
            if tmp2(I)>180
            annotation('textbox',[tmp5(I) tmp4(2)-0.06 w(I) 0.05],...   %%%%annotation('textbox', [tmp5(I) tmp4(2)-0.06 0.075 0.05],...
            'String',num2str(tmp6(I)),...
            'LineStyle','none',...
            'HorizontalAlignment','center',...
            'FontSize',8,...
            'FontName','Times New Roman');
            end
        end
        

   case 12
      set(ax2,'XTick',xtickqm,'TickLength',[0.005 0],'XtickLabel',[]);
%       set(ax2,'XTick',xtickmm);
      tmp1=union(xmin,union(xmax,xtickyi));
      tmp2=tmp1(2:end)-tmp1(1:end-1);
      tmp3=(tmp1(1:end-1)+tmp1(2:end))/2;
      tmp4=get(gca,'Position');
      tmp5=-0.04+tmp4(1)+tmp4(3)*(tmp3-xmin)./(xmax-xmin);
      tmp6=datevec(tmp3);
      tmp6=tmp6(:,1);
      if year_freq ==1
        for I=1:length(tmp3)
            if tmp2(I)>180
                annotation('textbox', [tmp5(I) tmp4(2)-0.05 0.075 0.05],...
                'String',num2str(tmp6(I)),...
                'LineStyle','none',...
                'HorizontalAlignment','center',...
                'FontSize',8,...
                'FontName','Times New Roman');
            end
        end
      elseif year_freq==2
          for I=1:2:length(tmp3)
              if tmp2(I)>180
                annotation('textbox', [tmp5(I) tmp4(2)-0.05 0.075 0.05],...
                'String',num2str(tmp6(I)),...
                'LineStyle','none',...
                'HorizontalAlignment','center',...
                'FontSize',8,...
                'FontName','Times New Roman');
              end
          end
      end
      
      line(ones(2,1)*xtickyi',myylim'*ones(1,length(xtickyi)),'Color','k','LineStyle',':');
     % line(ones(2,1)*xtickmi',myylim'*ones(1,length(xtickmi)),'Color','k','LineStyle',':');
end

  
    
    
    


    

% 
set([ax1 ax2],...
     'FontName','Times New Roman',...
     'FontSize',8);
 end







