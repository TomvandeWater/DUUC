%plotEdited
function out = plotEdited(x,y,color,linewidth,xlab,ylab,heading,gridset)

plot(x,y,color,'LineWidth',linewidth)
xlabel(xlab);ylabel(ylab);title(heading);
if strcmp(gridset,'on')
grid on;
else
    grid off;
end
out = 1;



