function [Data]=EyeRay(Z,NoM)
%In this function we propose an extension of Gagnon et al. method (Appl. 
%Optics, 2014). The presented method is non-symmetrical to optical axis 
%and all the optical surfaces are represented two-dimensionally by 
%Chebyshev polynomials as well as the rays of light that interact with the
%optical system.

%Algorithm:
% 1. Initial parameteres are defined;
% 2. Ray Tracing is performed for each meridian
%   2.1 Rotation of the surfaces;
%   2.2 Intersection (as a function of beta)
%   2.3 Refraction (for all surfaces except the last one)
%
%Input:
% Z: (1e7 by default)
% NoM: Number of Meridians (5 by default)
%
%Output:
%  Data - Intersections and surfaces
%
%===========
% Examples:
%===========
% [Data]=EyeRay(0,0,1e7,5); %Optical system centered at (0,0) for a source
% point at 1e7 over 5 meridians
% [Data]=EyeRay(0.02,0.01,1e7,10); %Optical system centered at (0.02,0.01)
% for a source point at 1e7 over 10 meridians
% 
% Note 1: For an off-axis analysis you can, for example, reduce the angle 
% of aperture to dombeta = [maxbeta-maxbeta/2, maxbeta] and to limit
% the range of the meridional angle (ia=80, fa=110)
%
% Note 2: The focal point changes when you move the optical system. Therefore
% if you want to see the spot diagram at the focal point you should
% calculate it and ajust retina position. The new focal point can easily be
% calculated through the circle of confusion.
%
% Note 3: The lateral movements of the optical system are restricted by
% surface's equations, namely by the zero at square root. For example, if you
% define EyeRay(0.2,0.3,1e7,20) you will get the following error message,
% "A change of sign/zero has been detected, unable to represent the
% result." In order to get free of it you should change the radius/asphericity
% 
%
% Graphical representation
% gPSF(Data,NoR) - Spot Diagram where NoR is the number of rays per meridian
% EyeRayShow(Data) - Schematic representation

% v1.01
% in: AGEYE, 608049.


%% Initial parameters
if nargin==0     
    Z=1e7;  %Source distance
    NoM=5;  %Number of Meridians
end

%Chebfun construction-time preferences
chebfunpref.setDefaults('factory');
chebfunpref.setDefaults('splitting', true)


dom=[-1 1 -1 1]; %Surface's domain
app=1.7; %Aperture diameter
r=app/2; %Radius of the optical aperture
maxbeta = atan(r/Z); %Angle between optical axis and periphery rays
dombeta = [-maxbeta, maxbeta]; % Angle of Aperture (Domain of beta) 
beta=chebfun(@(beta) beta,dombeta); %All possible angles defined as chebfun
p=0*beta + 1i*Z; %Light's initial position
I=exp(1i*(beta-pi/2)); %Light's initial direction
x=chebfun2(@(x,y) x,dom);  %Function necessary for the plane
NoS=7; %Number of surfaces

iI=I; %Ray's initial direction
iP=p; %Source's inital position

%% Ray tracing

%Allocates the variable Data in memory
for u=1:NoS
    Data{u}=chebfun(@(x) x,dombeta);
end

k=0; %Meridian

%Range of the Meridional angle
ia=0; %Inital angle
fa=180; %final angle
stp=(fa-ia)/NoM; %Interval between meridians

for th =ia:stp:fa-stp
k=k+1; %Meridian 

%Rotational matrix - there's a rotation of the 2D surface for each meridian
XX=chebfun2(@(x,y) x*cos(th*pi/180)+y*sin(th*pi/180),dom);
YY=chebfun2(@(x,y) -x*sin(th*pi/180)+y*cos(th*pi/180),dom);

%%  Model structurally similar to the human eye

%Anterior Cornea
R_ac=7.77;  %Anterior Radius (mm)
Q_ac=-0.18; %Asphericity
a_cornea = (2*R_ac + sqrt(4*R_ac^2 - 4*(1+Q_ac).*(XX.^2+YY.^2)))/(2*(1+Q_ac));
a_cornea=a_cornea-a_cornea(0,0);

% Real data acquired with Medmont 
% Taylor Coefficients for a Normal Cornea - 007N_sv_01
% Tcf=[0.0622092743495633,4.28916520878144e-06,-2.25938086697261e-07,-0.0646219022874312,5.29463643720144e-05,-0.0650068713872062,-0.000121476561591488,0.000387685723654671,0.000132182854115000,0.000236902022731613,0.000225589661496479,0.000718211202071288,3.20153620613121e-05,0.000339153780550980,-4.57729809123390e-05];

% Taylor Coefficients for a Astigmatic Cornea - 022N_KJ_01
% Tcf=[0.0619761031959612,-4.06224442136272e-05,-7.41344909029362e-06,-0.0644138873352629,-0.000605304329412603,-0.0663317170698532,0.000178259922118490,-4.26558757658122e-05,0.000158430177065628,-1.35003548524190e-05,-6.74371640006900e-05,0.000866473526160976,0.000154675889712618,-3.19585866015711e-05,-0.000201039194372111];

% %Taylor Coefficients for a Keratoconic Cornea - 006K_JA_OD_01
% Tcf=[0.0886236990114500,2.99804260817018e-05,-8.56434104791753e-05,-0.0791639864318572,0.0115421157238940,-0.0899294044157041,0.00241134995001902,0.00147824417505293,0.00331223834655794,0.00496445356471593,0.000614892657591541,0.00161976537063895,0.00281984225723166,-0.000422803683318245,0.00157177776565605];
%  
% % Taylor approximation using a polynomial up to 4th order
% A_Cornea = Tcf(1) + Tcf(2)*XX + Tcf(3)*YY + Tcf(4)*(XX.^2) + Tcf(5)*XX.*YY + Tcf(6)*(YY.^2) + Tcf(7)*(XX.^3) + Tcf(8)*YY.*(XX.^2) + Tcf(9)*XX.*(YY.^2) + Tcf(10)*(YY.^3) + Tcf(11)*(XX.^4) +Tcf(12)*(XX.^3).*YY +Tcf(13)*(XX.^2).*(YY.^2) +Tcf(14)*XX.*(YY.^3) +Tcf(15)*(YY.^4); 
% a_cornea=A_Cornea-A_Cornea(0,0);

%Posterior Cornea
R_pc=6.4;  %Radius(mm)
Q_pc=-0.6; %Asphericity
D_pc=-0.5; %Distance(mm)
p_cornea = (2*R_pc + sqrt(4*R_pc^2 - 4*(1+Q_pc).*(XX.^2 + YY.^2)))/(2*(1+Q_pc));
p_cornea = p_cornea-p_cornea(0,0) + D_pc;

%Anterior lens (cortex)
R_acl=12.4;  %Radius(mm)
Q_acl=-0.94; %Asphericity
D_acl=-3.66; %Distance(mm)
ac_lens=(2*R_acl + sqrt(4*R_acl^2 - 4*(1+Q_acl).*(XX.^2 + YY.^2)))/(2*(1+Q_acl));
ac_lens=ac_lens-ac_lens(0,0) + D_acl;

%Anterior lens (nucleus)
R_anl=3.04;  %Radius(mm)
Q_anl=-0.94; %Asphericity
D_anl=-4.86; %Distance(mm)
an_lens=(2*R_anl + sqrt(4*R_anl^2 - 4*(1+Q_anl).*(XX.^2 + YY.^2)))/(2*(1+Q_anl));
an_lens=an_lens-an_lens(0,0) + D_anl;
 
%Posterior lens (nucleus)
R_pnl=-2.10; %Radius(mm)
Q_pnl=0.96;  %Asphericity
D_pnl=-5.88; %Distance(mm)
pn_lens=(2*R_pnl - sqrt(4*R_pnl^2 - 4*(1+Q_pnl).*(XX.^2 + YY.^2)))/(2*(1+Q_pnl));
pn_lens=pn_lens-pn_lens(0,0) + D_pnl;

%Posterior lens (cortex)
R_pcl=-8.1;  %Radius(mm) 
Q_pcl=0.96;  %Asphericity
D_pcl=-7.68; %Distance(mm)
pc_lens=(2*R_pcl - sqrt(4*R_pcl^2 - 4*(1+Q_pcl).*(XX.^2 + YY.^2)))/(2*(1+Q_pcl));
pc_lens=pc_lens-pc_lens(0,0) + D_pcl;

%Retina
R_r=-12.82; %Radius(mm)
Q_r=0.26;   %Asphericity
D_r=-25.32; %Distance(mm)
retina=(2*R_r - sqrt(4*R_r^2 - 4*(1+Q_r).*(XX.^2 + YY.^2)))/(2*(1+Q_r));
retina=retina-retina(0,0) + D_r;

%Optical system
opt_sys={a_cornea,p_cornea,ac_lens,an_lens,pn_lens,pc_lens,retina};

%Refractive indices
ri=chebfun({@(x) 1, @(x) 1.3777, @(x) 1.3371,@(x) 1.3976,...
    @(x) 1.4033,@(x) 1.3976,@(x) 1.3377},-1*[1 a_cornea(0,0)...
     p_cornea(0,0) ac_lens(0,0) an_lens(0,0) pn_lens(0,0)...
      pc_lens(0,0) retina(0,0)],'splitting','on'); 
  
    for t = 1:NoS    
        
        Lens=opt_sys{t};
        p = chebfun(@(beta) intersection(Lens, p(beta), I(beta), x),dombeta,10,'vectorize');
        Data{t}(:,k)=p; % Intersection(X,Y)
        p=real(p)+1i*Lens(p); %Intersection(X,Z) 
        
        %Refraction 
        if t < NoS %Don't refract through the last surface  
            n1=(ri(-1*(Lens(0,0)+.0001))); %Refractive index 1 (before the Lens)
            n2=(ri(-1*(Lens(0,0)-.0001))); %Refractive index 2 (after the Lens)
            dz=diff(Lens,1,2); %Derivative on x
            N=atan(dz)+pi/2; %Normal of the surface
            Ni = N(Data{t}(:,k)); %Angle of the normal and x axis
            alpha1  = Ni-angle(I); %Angle between the normal and the ray before refraction
            alpha2 = asin(n1*sin(alpha1)/n2); %Refracted angle;
            theta = (Ni+alpha2)-pi; %Angle of the refracted ray;
            I = exp(1i*theta); %Light's new direction; 
        end       
    end
    
    I=iI; %Direction reset
    p=iP; %Position reset
end

%Add sufaces to the Output cell
Data{u+1}=opt_sys;
Data{u+2}=[ia,stp];

%Intersection function
function [p]=intersection(Lens, pb, Ib, x)
    a = imag(Ib)/real(Ib); % the slope of the ray's function
    f = imag(pb)+a*(x-real(pb)); % the ray's function    
    x2 = roots(Lens-f);  % the x value where the interface and ray intersect  
    p=x2(0);
end

end
