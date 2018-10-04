function makePDF(file)
    input = strcat( file, '.fig');
    h = open(input);
    ax = gca;
    ax.FontWeight = 'normal';
    ax.FontSize = 12;
%     ylabel('y position');
%     ylabel('Jain''s index');
    title('');
%     ylabel('Transmission Rate(b/s/HZ)');
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    output = strcat(file,'.pdf');
    print(h,output,'-dpdf','-r0');
    close;
end