close all

%For a 1/2 scale, the spar weight using caps is about 0.6% for a given mass
%if we assume it's a tapered spar, otherwise it's about 1.5% of the mass

n = 8; %safety factor

%wing properties
crootrange = 1./ [2 2 2];
spanrange = 12 ./ [2 2 2];
massrange = [10 15 25];

taper = .75; 
t_c = .13; %thickness to chord ratio


%material properties
syield = 1200*10^6; %Kevlar 1280Mpa  , Carbon Fiber ~1200MPa
rho = 1600;  %Kevlar 1380 kg/m^3 , Carbon Fiber 1600 kg/m^3
min_t = 0.8e-3; %minimum thickness of material
max_t = 10e-3;

%spar section properties
spartype = 3;% 1 for tubular spar, 2 for rectangular spar, 3 for sparcaps
sparwidth = .10; %width of spar caps as percentage of the root thickness
non_tapered_spar = 1; %0 for tapering spar, 1 for fixed width



npanels = 50; %number of panels in y direction


syms y; syms z;
syms t; syms b; syms W;



for l = 1:length(crootrange)
croot = crootrange(l);    
span = spanrange(l);
mass = massrange(l);

%Lift distribution
Ly = (4*W/b/pi) * sqrt(1-(2*y/b)^2);
%bending moment distribution
My = (4*W/b/pi) * b/2 * int((b/2*z-y)*sqrt(1-z^2),z,2*y/b,1);

%Available thickness distribution
Dy = .9 * (t_c) * croot * (1 - (1-taper)*2*y/b);
Dy_fun = inline(char(.98 * (t_c) * croot * (1 - (1-taper)*2*y/span)));

if(spartype == 1)
    Ixx = pi/64*(Dy^4 - (Dy-2*t)^4); %moment of inertia for hollow cylinder
elseif(spartype == 2)
    Ixx = 1/12*t*Dy^3; %moment of inertia for rectangle
elseif(spartype == 3) 
    %moment of inertia for two sparcaps
    Ixx = 1/12*(sparwidth*.9*t_c*croot)*(Dy^3-(Dy-2*t)^3); 
end

sy = My * Dy/2 / Ixx; %maximum stress 
sy_fun =  (inline (char (subs(sy-syield/n,{W,b},{mass*9.81,span})),'t','y'));

dy = (span/2)/npanels;
yrange = 0:dy:(span/2);
trange = zeros(1,length(yrange)+1);

tguess = (t_c)*croot*.5;

nrange = length(yrange);
if(non_tapered_spar)
    nrange = 1;
end
for k = 1:nrange
    
    try %Solve for required thickness to satisfy yield
        trange(k)  = ...
	fminbnd(@(t) abs(sy_fun(t,yrange(k))),0,tguess,OPTIMSET('TolX',1e-8));
    catch ex
        fprintf('error in solving for thickness ...\n');
    end
    if(isnan(trange(k)))
        fprintf('error in solving for thickness ...\n');
        trange(k:end)=trange(k-1);
        break;
    end
    %imposing a minimum gage thickness
    trange(k) = max(trange(k),min_t);

end
yrange = yrange(1:npanels);
trange = trange(1:npanels);
if(non_tapered_spar) %if non-tapered use the root size all across the span
    trange = trange(1)*ones(size(trange));
end
%Computing the mass. Factor of 2 to account for two 'half' wings
if(spartype == 1)
        sparmass = 2*Dy_fun(yrange)*trange'*pi*dy*rho;
elseif(spartype == 2)
        sparmass = 2*Dy_fun(yrange)*trange'*dy*rho;
elseif(spartype == 3)        %spar caps are at the top and bottom, so an extra *2
        sparmass = 4*sparwidth*Dy_fun(yrange)*trange'*dy*rho;
end

if(spartype == 3)
    fprintf('sparwidth = %f mm , ',sparwidth*0.9*t_c*croot*1e3);
end
fprintf('sparmass = %f (%f%% of mass)\n',sparmass, sparmass/mass*100);


Dyrange = subs(Dy,{b,y},{span,yrange});
Myrange = subs(My,{b,W,y},{span,mass*9.81,yrange});Myrange = Myrange/Myrange(1);
Lyrange = subs(Ly,{b,W,y},{span,mass*9.81,yrange});Lyrange = Lyrange/Lyrange(1);
syrange = subs(sy/syield,{b,W,y,t},{span,mass*9.81,yrange,trange});

figure;
subplot(2,1,1);
hold on;
plot(yrange,Lyrange,'b');
plot(yrange,Myrange,'k'); 
plot(yrange,syrange,'r');
plot([0 yrange(end)],[1/n 1/n],'--');
legend('Lift distribution','Bending Moment','1/(safety factor)');

subplot(2,1,2);
plot(yrange,trange*1e3,'o'); hold on;
xlabel('y'); ylabel('thickness (mm)');

end
