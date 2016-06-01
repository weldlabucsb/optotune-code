function [waist_x, waist_y]=guppyAnalysis_plots(data)
%Made 02/25/15
%Modified to rev3 to ask for ROI 5/20/2015, also flipped
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Select a file
    %close All;
    %[FileName,PathName] = uigetfile('*.tif','Select your Image');
    %file=strcat(PathName, FileName);
    %pic=imread(file);

    pic=data;
    
    image_matrix=double(pic);  
    image_matrix=sum(image_matrix,3);
    
    data_vertical=transpose(sum(image_matrix,2));    
    data_horizontal=sum(image_matrix,1);
    
    pixelsize=8.3E-3; %pixel size in mm for Guppy
    
    [outparamsgaussx, fitfunc_x]=getFitParams(data_horizontal);
    [outparamsgaussy, fitfunc_y]=getFitParams(data_vertical); 
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Compute the beam parameters.  
    waist_x=outparamsgaussx(1);%outparamsgauss currently gives the sigma which is 1/4 of the 1/e^2 diamter
    height_x=outparamsgaussx(2);
    center_x=outparamsgaussx(3);
    bkgd_x=outparamsgaussx(4);
    
    waist_y=outparamsgaussy(1);%outparamsgauss currently gives the sigma which is 1/4 of the 1/e^2 diamter
    height_y=outparamsgaussy(2);
    center_y=outparamsgaussy(3);   
    bkgd_y=outparamsgaussx(4);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Make Plots in X and Y Direction

    delete(figure(1));
    figure(1);
    set(gcf, 'Visible', 'Off'); 
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    set(gcf,'numbertitle','off','name','Beam Analysis'); 
    
    %Make the horizontal (summed X) plot.
    subplot(2,2,1);
    thisx=(1:length(data_horizontal/height_x))*pixelsize;
    set(gca,'FontSize',12);
    p=plot(thisx,data_horizontal/height_x,'k');
    %axis equal;
    hold on;
    hx=plot(thisx,fitfunc_x/height_x,'r'); set(hx,'LineWidth',2);
    %axis equal;
    xlabel('Position (mm)');
    ylabel('Normalized Counts');
    title(['Summed Intensity Profile X']);
    
    
    str1='Pixel Size = 8.3\mum';% for Guppy
    text(.05, .95, str1, 'Units', 'Normalized');
    str2=['Waist X = ' num2str(waist_x) ' mm'];
    text(.05, .9, str2, 'Units', 'Normalized');
    str3=['Bkgd = ' num2str(bkgd_x/height_x)];
    text(.05, .85, str3, 'Units', 'Normalized');

    %Make the vertical (summed y) plot.
    subplot(2,2,4);
    data_vertflip=fliplr(data_vertical);
    thisy=(1:length(data_vertflip/height_y))*pixelsize;
    set(gca,'FontSize',12);
    p=plot(data_vertflip/height_y,thisy,'k');
    %axis equal;
    hold on;
    hx=plot(fliplr(fitfunc_y/height_y),thisy,'r'); set(hx,'LineWidth',2);
    %axis equal;
    ylabel('Position (mm)');
    xlabel('Normalized Counts');
    title(['Summed Intensity Profile Y']);
    
    str1='Pixel Size = 8.3\mum'; %for Guppy
    text(.7, .95, str1, 'Units', 'Normalized');
    str2=['Waist Y = ' num2str(waist_y) ' mm'];
    text(.7, .9, str2, 'Units', 'Normalized');
    str3=['Bkgd = ' num2str(bkgd_y/height_y)];
    text(.7, .85, str3, 'Units', 'Normalized');

    %Show the image and draw a circle demonstrating the beam waist.
    subplot(2, 2, 3) 
    imagesc(pic);  
    caxis([0 330]);
    axis equal tight;
    xlabel('x px');
    ylabel('y px');
    colorbar
    
    set(gcf,'PaperPositionMode','auto');     
    
    keyboard
    
    function [outparamsgauss, fitfunc]=getFitParams(summedOD) 
    mythisx=(1:length(summedOD))*pixelsize;
    sigmaguess=max(mythisx)/10;
    ampguess=max(summedOD(~isinf(summedOD)));
    xguess=find(summedOD==max(summedOD)); xguess=xguess(1); xguess=mythisx(xguess);

    offsetguess=0; 
    guessparams=[sigmaguess ampguess xguess, offsetguess]; % 1/23/15 DMW removing offset from params    

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
    rmserror=sqrt(sum(error1.^2));
end
    
function [outgauss]=makegauss(params,xvec)
    waist=params(1);
    gaussheight=params(2);
    centerx=params(3);
    yoffset=params(4);

    x=xvec-centerx;
    outgauss = yoffset + gaussheight * exp(-(2*x.^2)./waist^2);
end



