function [waist_x, waist_y]=guppyAnalysis_rev6(pic)

image_matrix=double(pic);
image_matrix=sum(image_matrix,3);

data_vertical=transpose(sum(image_matrix,2));
data_horizontal=sum(image_matrix,1);

pixelsize=8.3E-3; %pixel size in mm for Guppy

[outparamsgaussx, ~]=getFitParams(data_horizontal);
[outparamsgaussy, ~]=getFitParams(data_vertical);

%Compute the beam parameters.

waist_x=abs(outparamsgaussx(1));%outparamsgauss currently gives the sigma which is 1/4 of the 1/e^2 diamter
waist_y=abs(outparamsgaussy(1));%outparamsgauss currently gives the sigma which is 1/4 of the 1/e^2 diamter

    function [outparamsgauss, fitfunc]=getFitParams(summedOD)
        mythisx=(1:length(summedOD))*pixelsize;
        sigmaguess=max(mythisx)/10;
        ampguess=max(summedOD(~isinf(summedOD)));
        centerguess=find(summedOD==max(summedOD));
        centerguess=centerguess(1);
        centerguess=mythisx(centerguess);
        
        offsetguess=0;
        guessparams=[sigmaguess ampguess centerguess offsetguess];
        
        [outparamsgauss,fitfunc] = fitgaussian(mythisx,summedOD,guessparams);
    end
end

function [outparamsgauss, fitfunc] = fitgaussian(xvec,yvec,guessparams)
% this function fits a gaussian to a vector
% params is as follows: [gausswidth gaussheight centerx yoffset]
outparamsgauss=fminsearch(@(params) lsqmingetgauss(params,xvec,yvec),guessparams);
fitfunc=makegauss(outparamsgauss,xvec);
end

function [rmserror]=lsqmingetgauss(params,xvec,yvec)
fitfunc=makegauss(params,xvec);
yvec(isinf(yvec))=fitfunc(isinf(yvec));
yvec(isnan(yvec))=fitfunc(isnan(yvec));
error1=fitfunc-yvec;
% rmserror=sqrt(sum(error1.^2));
rmserror=(sum(error1.^4))^(0.25);
end

function [outgauss]=makegauss(params,xvec)
waist=params(1);
gaussheight=params(2);
centerx=params(3);
yoffset=params(4);

x=xvec-centerx;
outgauss = yoffset + gaussheight * exp(-(2*x.^2)./waist^2);
end