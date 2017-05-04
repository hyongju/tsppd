function []= draw_tour( tour,n,vert,name)
    if nargin==3,
        name = ' ';
    end
    figure,
    plot(vert(:,1),vert(:,2),'ok','MarkerSize',10,'LineWidth',2); hold on;

    for i = 1:length(tour)-1
        dif = vert(tour(i+1),:)-vert(tour(i),:);
        quiver(vert(tour(i),1),vert(tour(i),2),0.1*dif(1)/norm(dif),0.1 *dif(2)/norm(dif),0, 'MaxHeadSize', 1/norm(dif),'LineWidth',2);hold on;
        line([vert(tour(i),1) vert(tour(i+1),1)],[vert(tour(i),2) vert(tour(i+1),2)],'Color','black');hold on;
    end

    for i = 1:size(vert,1)
        if (mod(i-1,n)) == 0 && i ~= 1
            prtVal = n;
        else
            prtVal = mod(i-1,n);
        end
        if i == 1
            str = 'I';
            col = 'black';
        elseif i <= n+1
            str = 'P';
            col = 'red';
        else
            str = 'D';
            col = 'blue';
        end
        text(vert(i,1)+0.02,vert(i,2)+0.02, sprintf('%.0f %s',prtVal,str),'Color',sprintf('%s',col),'FontSize',16);
    end

    axis('equal');
    axis([0 1 0 1]);
    set(gca,'FontSize',16);
    title(name)


end

