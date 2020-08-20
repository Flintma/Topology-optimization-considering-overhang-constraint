function [OverhangC,doc_dx1] = OC(x,P1)
%%rotation
rotate=0;
%bottomtotop rotate=0
%lefttoright  rotate=1
%toptobottom  rotate=2
%righttoleft rotate=3
x=rot90(x,rotate);
%%Preparation&Base density calculation
R=0.97; %constraint critical rate
si=size(x);%compute size of the design domain [row column]
%%Prepare xi domain for analysis
xi=zeros(si(1)+1, si(2)+4);
xi(end,:)=ones(1,si(2)+4);
for i = 1:si(1)
for j=3:si(2)+2
xi(i,j)=x(i,j-2);
end
end
%%Calculate maximum density of base elements
N=700;
for i=1:si(1)
for j=3:si(2)+2
m=i+1;
for  n = j-2:j+2
yy(n)=xi(m,n);%record base element density 
end
yg(i,j-2)=(yy(j-1)*exp(N*yy(j-1))+yy(j)*exp(N*yy(j))+yy(j+1)*exp(N*yy(j+1)))/(exp(N*yy(j-1))+exp(N*yy(j))+exp(N*yy(j+1)));
end
end
%% Calculate OC
Phi=zeros(si(1),si(2));Sup=0;N0_sup=0;
for i=1:si(1)
for j=1:si(2)
Phi(i,j)=x(i,j)-yg(i,j);
if Phi(i,j)>0
Sup=Sup+Phi(i,j)*(yg(i,j)+x(i,j))*P1;
else
N0_sup=N0_sup-Phi(i,j)*(yg(i,j)+x(i,j)); 
end
end
end
OverhangC=R-N0_sup/(Sup+N0_sup);
%% Sensitivity analysis
dNdx=zeros(si(1),si(2));dSdx=zeros(si(1),si(2));
for o=1:si(1)
for p=3:si(2)+2
if Phi(o,p-2)<-10^(-4) || Phi(o,p-2)>10^(-4)
%dS/Ndx (in every element)
%Analytical solution
if o==1
if Phi(o,p-2)>0
dxyg(o,p-2)=P1*2*x(o,p-2);
dNdx(o,p-2)=dNdx(o,p-2)+dxyg(o,p-2);
else
dxyg(o,p-2)=-1*(2*x(o,p-2));
dSdx(o,p-2)=dxyg(o,p-2);
end
else
for i=p-2:p+2
dy(i)=xi(o,i);%record base element density
end
yg1(o,p-2)=(dy(p-2)*exp(N*dy(p-2))+dy(p-1)*exp(N*dy(p-1))+dy(p)*exp(N*dy(p)))/(exp(N*dy(p-2))+exp(N*dy(p-1))+exp(N*dy(p)));
dyg1(o,p-2)=(exp(N*dy(p))*(1+N*(dy(p)-yg1(o,p-2))))/(exp(N*dy(p-2))+exp(N*dy(p-1))+exp(N*dy(p)));  
yg2(o,p-2)=(dy(p-1)*exp(N*dy(p-1))+dy(p)*exp(N*dy(p))+dy(p+1)*exp(N*dy(p+1)))/(exp(N*dy(p-1))+exp(N*dy(p))+exp(N*dy(p+1)));
dyg2(o,p-2)=(exp(N*dy(p))*(1+N*(dy(p)-yg2(o,p-2))))/(exp(N*dy(p-1))+exp(N*dy(p))+exp(N*dy(p+1)));  
yg3(o,p-2)=(dy(p)*exp(N*dy(p))+dy(p+1)*exp(N*dy(p+1))+dy(p+2)*exp(N*dy(p+2)))/(exp(N*dy(p))+exp(N*dy(p+1))+exp(N*dy(p+2)));
dyg3(o,p-2)=(exp(N*dy(p))*(1+N*(dy(p)-yg3(o,p-2))))/(exp(N*dy(p))+exp(N*dy(p+1))+exp(N*dy(p+2)));
if p==3
dyg1(o,p-2)=0;
if Phi(o-1,p-2)>0
dyg2(o,p-2)=-P1*2*yg2(o,p-2)*dyg2(o,p-2);
dNdx(o,p-2)=dNdx(o,p-2)+dyg2(o,p-2);
else
dyg2(o,p-2)=-(-2*yg2(o,p-2)*dyg2(o,p-2));
dSdx(o,p-2)=dSdx(o,p-2)+dyg2(o,p-2);
end
if Phi(o-1,p-1)>0
dyg3(o,p-2)=-P1*2*yg2(o,p-2)*(dyg3(o,p-2));
dNdx(o,p-2)=dNdx(o,p-2)+dyg3(o,p-2);
else
dyg3(o,p-2)=-(-2*yg3(o,p-2)*dyg3(o,p-2));
dSdx(o,p-2)=dSdx(o,p-2)+dyg3(o,p-2);
end
elseif p==si(2)+2
dyg3(o,p-2)=0;
if Phi(o-1,p-2)>0
dyg2(o,p-2)=-P1*2*yg2(o,p-2)*(dyg2(o,p-2));
dNdx(o,p-2)=dNdx(o,p-2)+dyg2(o,p-2);
else
dyg2(o,p-2)=-(-2*yg2(o,p-2)*dyg2(o,p-2));
dSdx(o,p-2)=dSdx(o,p-2)+dyg2(o,p-2);
end
if Phi(o-1,p-3)>0
dyg1(o,p-2)=-P1*2*yg2(o,p-2)*(dyg1(o,p-2));
dNdx(o,p-2)=dNdx(o,p-2)+dyg1(o,p-2);
else
dyg1(o,p-2)=-(-2*yg1(o,p-2)*dyg1(o,p-2));
dSdx(o,p-2)=dSdx(o,p-2)+dyg1(o,p-2);
end
else
if Phi(o-1,p-2)>0
dyg2(o,p-2)=-P1*2*yg2(o,p-2)*(dyg2(o,p-2));
dNdx(o,p-2)=dNdx(o,p-2)+dyg2(o,p-2);
else
dyg2(o,p-2)=-(-2*yg2(o,p-2)*dyg2(o,p-2));
dSdx(o,p-2)=dSdx(o,p-2)+dyg2(o,p-2);
end
if Phi(o-1,p-1)>0
dyg3(o,p-2)=-P1*2*yg2(o,p-2)*(dyg3(o,p-2));
dNdx(o,p-2)=dNdx(o,p-2)+dyg3(o,p-2);
else
dyg3(o,p-2)=-(-2*yg3(o,p-2)*dyg3(o,p-2));
dSdx(o,p-2)=dSdx(o,p-2)+dyg3(o,p-2);
end
if Phi(o-1,p-3)>0
dyg1(o,p-2)=-P1*2*yg2(o,p-2)*(dyg1(o,p-2));
dNdx(o,p-2)=dNdx(o,p-2)+dyg1(o,p-2);
else
dyg1(o,p-2)=-(-2*yg1(o,p-2)*dyg1(o,p-2));
dSdx(o,p-2)=dSdx(o,p-2)+dyg1(o,p-2);
end
end
if Phi(o,p-2)>0
dNdx(o,p-2)=dNdx(o,p-2)+P1*2*x(o,p-2);
else
dSdx(o,p-2)=dSdx(o,p-2)+(-1*(2*x(o,p-2)));
end
end         
doc_dx1(o,p-2)=dNdx(o,p-2)/(N0_sup+Sup)-(dNdx(o,p-2)+dSdx(o,p-2))*Sup/(N0_sup+Sup)^2;    
else
%Numerical method
            x(o,p-2)=x(o,p-2)+10^(-12);
%%Prepare xi domain for analysis
xi=zeros(si(1)+1, si(2)+4);
xi(end,:)=ones(1,si(2)+4);
for i = 1:si(1)
    for j=3:si(2)+2
        xi(i,j)=x(i,j-2);
    end
end
%%Calculate maximum density of base elements
N=700;
for i=1:si(1)
    for j=3:si(2)+2
        m=i+1;
        for  n = j-2:j+2
                yy(n)=xi(m,n);%record base element density 
        end
         yg(i,j-2)=(yy(j-1)*exp(N*yy(j-1))+yy(j)*exp(N*yy(j))+yy(j+1)*exp(N*yy(j+1)))/(exp(N*yy(j-1))+exp(N*yy(j))+exp(N*yy(j+1)));
    end
end
%% Calculate OC
Phi2=zeros(si(1),si(2));Sup=0;N0_sup=0;
for i=1:si(1)
    for j=1:si(2)
        
           Phi2(i,j)=x(i,j)-yg(i,j);
        if Phi2(i,j)>0
            Sup=Sup+Phi2(i,j)*(yg(i,j)+x(i,j))*P1;
        else
            N0_sup=N0_sup-Phi2(i,j)*(yg(i,j)+x(i,j));  
        end
        
    end
end

         OverhangCC=R-N0_sup/(Sup+N0_sup);
          dxyg(o,p-2)=(OverhangCC-OverhangC)/10^(-12);    
          x(o,p-2)=x(o,p-2)-10^(-12);
          doc_dx1(o,p-2)=dxyg(o,p-2);
           
        end

    end
end
doc_dx1=rot90(doc_dx1,-rotate); 

