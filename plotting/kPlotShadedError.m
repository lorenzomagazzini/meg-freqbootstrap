function kPlotShadedError(t,m1,s1,cs1)

% t: An array of values for the x axis
%
% m1: The y values for the line e.g. the mean
% s1: The height above and below the line to shade e.g. the standard error
% cs1: The RGB colour for the shading for e.g. [0.8 0.8 1.0] is light blue
%  
%
% example usage:
%---------------------------------------------------------------
%	figure;
%	hold on;
%
%	kPlotShadedError(FreqArray,mean1,error1,[0.8 0.8 1.0]); %plot the error for mean 1 in light blue
%	kPlotShadedError(FreqArray,mean2,error2,[1.0 0.8 0.8]); %plot the error for mean2 in pink
%
%	plot(FreqArray,mean1,'Color',[0 0 1]); %plot the solid line for mean 1 in solid blue
%	plot(FreqArray,mean2,'Color',[1 0 0]); %plot the solid line for mean 2 in solid red
%
%	title('Power versus Frequency');
%	xlim([0 max(FreqArray)];
%	ylim([0 100]);
%	hold off;
%---------------------------------------------------------------

	deltat=t(2)-t(1);
	lt=length(t);	

	for i=1:lt,
		revi=(lt+1)-i;
		x(i)=t(i);
		x(i+lt)=t(revi);
		y(i)=m1(i)+s1(i);
		y(i+lt)=m1(revi)-s1(revi);
	end
	x(isnan(x)) = 0;
	y(isnan(y)) = 0;
	fill(x,y,cs1,'EdgeColor',cs1);

end
