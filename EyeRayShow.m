function [] = EyeRayShow(Data)
    %Function to plot chebfun-based ray tracing tridimensionally
     
     figure,view(3),hold on
     
     
     %% Representation of the surfaces on a circular domain
     Rpupil=1; %Pupil's radius
     x=linspace(-Rpupil,Rpupil,100);
     [X2,Y2]=meshgrid(x);
     ind=(X2.^2+Y2.^2)>=Rpupil;
     
     R1=7.77;  %Radius
     Q1=-0.18; %Asphericity
     Z1=(2*R1 + sqrt(4*R1^2 - 4*(1+Q1).*((X2).^2 + (Y2).^2)))/(2*(1+Q1));
     Z1=Z1-max(max(Z1));
     Z1(ind)=NaN;
     
     R2=6.4;  %Radius
     Q2=-0.6; %Asphericity
     D2=-0.5;
     Z2=(2*R2 + sqrt(4*R2^2 - 4*(1+Q2).*((X2).^2 + (Y2).^2)))/(2*(1+Q2));
     Z2=Z2-max(max(Z2))+D2;
     Z2(ind)=NaN;
          
     R3=12.4;  %Radius
     Q3=-0.94; %Asphericity
     D3=D2-3.16;
     Z3=(2*R3 + sqrt(4*R3^2 - 4*(1+Q3).*((X2).^2 + (Y2).^2)))/(2*(1+Q3));
     Z3=Z3-max(max(Z3))+D3;
     Z3(ind)=NaN;
     
     R4=9.04;  %Radius
     Q4=-0.94; %Asphericity
     D4=D3-1.2;
     Z4=(2*R4 + sqrt(4*R4^2 - 4*(1+Q4).*((X2).^2 + (Y2).^2)))/(2*(1+Q4));
     Z4=Z4-max(max(Z4))+D4;
     Z4(ind)=NaN;
     
     R5=-8.10;  %Radius
     Q5=0.96;   %Asphericity
     D5=D4-1.02;
     Z5=(2*R5 + sqrt(4*R5^2 - 4*(1+Q5).*((X2).^2 + (Y2).^2)))/(2*(1+Q5));
     Z5=Z5-max(max(Z5))+D5;
     Z5(ind)=NaN;
     
     R6=-8.1;   %Radius 
     Q6=0.96;   %Asphericity
     D6=D5-1.8;
     Z6=(2*R6 + sqrt(4*R6^2 - 4*(1+Q6).*((X2).^2 + (Y2).^2)))/(2*(1+Q6));
     Z6=Z6-max(max(Z6))+D6;
     Z6(ind)=NaN;
     
     R7=-12.82;  %Radius
     Q7=0.26;    %Asphericity
     D7=-25.32;
     Z7=(2*R7 + sqrt(4*R7^2 - 4*(1+Q7).*((X2).^2 + (Y2).^2)))/(2*(1+Q7));
     Z7=Z7-max(max(Z7))+D7;
     Z7(ind)=NaN;
      
     surf(X2,Y2,Z1),shading interp,alpha(.5) 
     surf(X2,Y2,Z2),shading interp,alpha(.5)
     surf(X2,Y2,Z3),shading interp,alpha(.5) 
     surf(X2,Y2,Z4),shading interp,alpha(.5)
     surf(X2,Y2,Z5),shading interp,alpha(.5) 
     surf(X2,Y2,Z6),shading interp,alpha(.5) 
     surf(X2,Y2,Z7),shading interp,alpha(.5) 

%     Representation of the chebfun surfaces on a rectangle domain. In 
%     order to run it you should comment the previous code and uncomment
%     the following code:

%     for k=1:length(Data{end-1})
%         surface=Data{end-1}{k};
%         plot(surface)
%     end

%%  Representation of some rays and intersection points 
    
    spc=Data{end}(2); %Rotational angle
    ia=Data{end}(1);
    D=domain(Data{1,1}(:,1));  %Data domain
    X0=0;                      %Source point coordinates
    Y0=0;
    Z0=1e7;
    
    %Source point
    plot3(X0,Y0,Z0,'.k','markersize',40),xlabel('mm'),ylabel('mm'),zlabel('mm')
    zlim([-26 2])
    
    
    for k=1:length(Data)-2 %for each surface
        Lens=Data{end-1}{k};
        
        for m=1:size(Data{1,1},2) %for each meridian
            for h=linspace(D(1),D(2),8)
                
                %Intersection points
                S=Data{1,k}(:,m);
                th=ia + (m-1)*spc;
                X=chebfun(@(x) real(S(x))*cos(th*pi/180)+imag(S(x))*sin(th*pi/180),D);  %Domain's r %%%%%%%%
                Y=chebfun(@(x) -real(S(x))*sin(th*pi/180)+imag(S(x))*cos(th*pi/180),D);
                plot3(X,Y,Lens(X,Y),'.k','markersize',20)
                
                if k==1 %Initial Beam
                    Z=Lens(X,Y);
                    plot3([X0 X(h)],[Y0 Y(h)],[Z0 Z(h)],'b','linewidth',2)
                else %Difracted Rays
                    Z=Lens(X,Y);
                    So=Data{1,k-1}(:,m);
                    X0=chebfun(@(x) real(So(x))*cos(th*pi/180)+imag(So(x))*sin(th*pi/180),D);  %Domain's r %%%%%%%%
                    Y0=chebfun(@(x) -real(So(x))*sin(th*pi/180)+imag(So(x))*cos(th*pi/180),D);
                    Z0=Lens0(X0,Y0);
                    plot3([X0(h) X(h)],[Y0(h) Y(h)],[Z0(h) Z(h)],'b','linewidth',2)
                end
            end
        end
        Lens0=Data{end-1}{k};
    end
end
