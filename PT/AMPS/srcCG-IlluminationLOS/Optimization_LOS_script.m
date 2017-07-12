clc;clear all; close all; 

%%%%%%%%%%%%%%%%%%%%%%Definitions%%%%%%%%%%%%%%%%%%%%%%%%%%
H2O=0;
CO2=1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%Modes%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modeSpecies=CO2; %H2O, CO2

nVar=5; %5, 3
nRows=129; %129 %104; %133;
nColumns=256;

nPixelsX=256;
nPixelsY=256;

rowStart=62; %62 %87; %58;
nDev=25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Cmatrix.dat')

medfilt2(Cmatrix);

for i=1:nPixelsX
    for j=1:nPixelsY
        C(i,j) = Cmatrix(i+(j-1)*nPixelsY,1);
    end
end

%fileID = fopen('I1_0037442903_H2O_CO2.dat'); %2015-04-12T07:14:00 - 133 58
%fileID = fopen('I1_00387481450.CAL_H2O_CO2_f3.dat'); %2015-04-12T17:57:00 - 104 87
%fileID = fopen('I1_00387554894.CAL_H2O_CO2_f3.dat'); %2015-04-13T14:21:00 - 133 58
fileID = fopen('I1_00387379010.CAL_H2O_CO2_f3.dat'); %2015-04-11T13:40:00 - 129 62
A=fread(fileID,nVar*nRows*nColumns,'float');

array=zeros(nVar,nRows,nColumns);

for i=1:nVar
    for j=1:nRows
        for k=1:nColumns
            array(i,j,k)=A((i-1)*nColumns*nRows+(j-1)*nColumns+k);
        end
    end
end

if (modeSpecies==H2O)
    species(:,:)=array(2,:,:)*1.0e20;
elseif (modeSpecies==CO2)
    species(:,:)=array(5,:,:)*1.0e19; %CO2
end
    

species = flipud(species);

colatitude=(0:5*pi/180:pi);
longitude=(0:5*2*pi/360:2*pi);

ct=1;

for i=1:length(colatitude)
    for j=1:length(longitude)
        sphericalHarmonic(ct,1)=1; %Y00
        sphericalHarmonic(ct,2)=sin(colatitude(i))*sin(longitude(j));%Y1-1
        sphericalHarmonic(ct,3)=cos(colatitude(i)); %Y10
        sphericalHarmonic(ct,4)=sin(colatitude(i))*cos(longitude(j)); %Y11
        sphericalHarmonic(ct,5)=sin(colatitude(i))^2.0*sin(2*longitude(j));%Y2-2
        sphericalHarmonic(ct,6)=cos(colatitude(i))*sin(colatitude(i))*sin(longitude(j));%Y2-1
        sphericalHarmonic(ct,7)=3*cos(colatitude(i))^2.0-1.0;%Y20
        sphericalHarmonic(ct,8)=cos(colatitude(i))*sin(colatitude(i))*cos(longitude(j));%Y21
        sphericalHarmonic(ct,9)=sin(colatitude(i))^2.0*cos(2*longitude(j));%Y22
        sphericalHarmonic(ct,10)=sin(colatitude(i))^3.0*sin(3*longitude(j));%Y3-3
        sphericalHarmonic(ct,11)=sin(colatitude(i))^2.0*cos(colatitude(i))*sin(2*longitude(j));%Y3-2
        sphericalHarmonic(ct,12)=(5*cos(colatitude(i))^2.0-1)*sin(colatitude(i))*sin(longitude(j));%Y3-1
        sphericalHarmonic(ct,13)=5*cos(colatitude(i))^3.0-3*cos(colatitude(i));%Y30
        sphericalHarmonic(ct,14)=(5*cos(colatitude(i))^2.0-1)*sin(colatitude(i))*cos(longitude(j));%Y31
        sphericalHarmonic(ct,15)=sin(colatitude(i))^2.0*cos(colatitude(i))*cos(2*longitude(j));%Y32
        sphericalHarmonic(ct,16)=sin(colatitude(i))^3.0*cos(3*longitude(j));%Y33
        sphericalHarmonic(ct,17)=sin(colatitude(i))^4.0*sin(4*longitude(j));%Y4-4
        sphericalHarmonic(ct,18)=sin(colatitude(i))^3.0*cos(colatitude(i))*sin(3*longitude(j));%Y4-3
        sphericalHarmonic(ct,19)=sin(colatitude(i))^2.0*(7*cos(colatitude(i))^2.0-1)*sin(2*longitude(j));%Y4-2
        sphericalHarmonic(ct,20)=sin(colatitude(i))*(7*cos(colatitude(i))^3.0-3*cos(colatitude(i)))*sin(longitude(j));%Y4-1
        sphericalHarmonic(ct,21)=35*cos(colatitude(i))^4.0-30*cos(colatitude(i))^2.0+3;%Y40
        sphericalHarmonic(ct,22)=sin(colatitude(i))*(7*cos(colatitude(i))^3.0-3*cos(colatitude(i)))*cos(longitude(j));%Y41
        sphericalHarmonic(ct,23)=sin(colatitude(i))^2.0*(7*cos(colatitude(i))^2.0-1)*cos(2*longitude(j));%Y42
        sphericalHarmonic(ct,24)=sin(colatitude(i))^3.0*cos(colatitude(i))*cos(3*longitude(j));%Y43
        sphericalHarmonic(ct,25)=sin(colatitude(i))^4.0*cos(4*longitude(j));%Y44
        ct=ct+1;
    end
end

for i=1:nRows
    for j=1:nColumns
        if(species(i,j)~=species(i,j) || species(i,j)<=0) 
            data(i+(j-1)*nRows) = 0.0;
        else
            data(i+(j-1)*nRows) = species(i,j);
        end
    end
end

for i=1:nRows
    for j=1:nColumns
        for k=1:nDev
            if(species(i,j)~=species(i,j) || species(i,j)<=0)
                Coptim(i+(j-1)*nRows,k) = 0.0;
            else
                Coptim(i+(j-1)*nRows,k) = Cmatrix(rowStart+i+(j-1)*nPixelsY,k);%Cmatrix(nPixelsY*(rowStart+i-1)+j-1,k);
            end
        end
    end
end
 

[m, n]=size(Coptim);

A=-0.0*ones(length(sphericalHarmonic(:,1)),1);

[x,resnorm,residual,exitflag,output,lambda]=lsqlin(Coptim,data,-sphericalHarmonic,A,[],[],[],[],[],optimoptions(@lsqlin,'Algorithm','interior-point','MaxIter',500)); %%TolFun default 2.2204e-14

Cfinal=zeros(nRows,nColumns);
for i=1:nRows
    for j=1:nColumns
        for k=1:nDev
            Cfinal(i,j) = Cfinal(i,j) + x(k)*Cmatrix(rowStart+i+(j-1)*nPixelsY,k);
        end
    end 
end

exit;
