function h = hverrorbar(xMean,yMean,hError,vError,color)


herrorbar(xMean,yMean,hError,hError,color);
r = errorbar(xMean,yMean,vError,color,'LineWidth',2);
h = r(1);

end