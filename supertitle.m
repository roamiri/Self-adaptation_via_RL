function supertitle(mytitle,varargin)
% When using subtitle('MY TITLE','PorpertyName','PropertyValue'...), or
% subtitle('MY TITLE') after a group of subplots, then it provides a title
% MY TITLE with any property used that is defined in the original title
% function in Matlab, but without affecting the titles, xlables and ylabels
% of any of the subplots. 
%
% Make sure use the function after the group of subplots.
%
% Example
% x = 0:0.01:6;
% subplot(221), plot(x,sin(x)), xlabel('x'), ylabel('sin(x)'), title('sin(x)')
% subplot(222), plot(x,cos(x)), xlabel('x'), ylabel('cos(x)'), title('cos(x)')
% subplot(223), plot(x,sin(2*x)), xlabel('x'), ylabel('sin(2x)'), title('sin(2x)')
% subplot(224), plot(x,cos(2*x)), xlabel('x'), ylabel('cos(2x)'), title('cos(2x)')
% subtitle('Single title on top','FontSize',12,'Color','r')
%
% Copyright @ Md Shoaibur Rahman (shaoibur@bcm.edu)

axes('Units','Normal');
h1 = title(mytitle,varargin{1:length(varargin)});
h2 = xlabel('FBS Numbers','FontSize',12);
h3 = ylabel('Mean transmission rate (b/s/Hz)','FontSize',12);
set(gca,'visible','off');
set(h1,'visible','on');
set(h2,'visible','on');
set(h3,'visible','on');

