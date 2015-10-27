function [] = gPSF(Data,NoR)    
    %Function to display the Spot Diagram of the last surface
    
    if nargin < 2
        NoR=25;             %Number of rays by default
    end
    
    ia=Data{end}(1);   %Rotational angle
    spc=Data{end}(2);   %Rotational angle   
    
    Ints=Data{1,end-2};     %Intersection with the last surface
    D=domain(Ints);         %Domain
    
    figure,hold on
    for k=0:size(Ints,2)-1; %for each meridian
    S=Ints(:,k+1);          %Data at meridian k          
    th=k*spc;               %Angle
    X=chebfun(@(x) real(S(x))*cos(th*pi/180)-imag(S(x))*sin(th*pi/180),D); 
    Y=chebfun(@(x) real(S(x))*sin(th*pi/180)+imag(S(x))*cos(th*pi/180),D);
    plot(X(linspace(D(2),D(1),NoR)),Y(linspace(D(2),D(1),NoR)),'.'),xlabel('mm'),ylabel('mm')   
    end
   
    gcf,title('Spot Diagram')
end
