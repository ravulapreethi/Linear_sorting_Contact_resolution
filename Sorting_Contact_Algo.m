% Code for implementing linear search algorithm with DEM
% In this version (code author: Preethi Ravula, 2017) no need to use sorting
% Matlab parallel toobox parfor loop - ContactalgoV6_Preethi, 
% Cuda threads ContactalgoV7_Preethi

close all
clear all
clc
%% limits for computational grid that the particles might be in over the time
xmin = 1; ymin = 1; 
xmax = 1.2; ymax = 1.2; Rmax = 0.0025;

%% read particle positions and radius
R = xlsread('Pos_array');
%R = xlsread('Hexsample');
R= R(R(:,1)<=1.194,:); %% 1.1744
R= R(R(:,2)<=1.2,:);
% R= R(R(:,1)>=1.09,:); %% 1.1744
% R= R(R(:,2)>=1.09,:);
X = R(:,1); Y = R(:,2);
pos = [X(:) Y(:)]; Np = length(X(:));

%% Initalize forces on particles
force = zeros(Np,1);
force2 = zeros(Np,1); % this is to compare with N^2 algo
%%
Lx=xmax-xmin; Ly=ymax-ymin;                 % total lenght of the grid                                    
Nx=floor((xmax-xmin)/0.0035);               % No of boxes in x -direction
Ny=floor((ymax-ymin)/0.0035);               % No of boxes in y -direction

%% find coordinates (Ix,Iy) of particles in the grid domain (e.g., (1,1), (1,2)....(Nx,Ny))
Xp=(Nx*(pos(:,1)-xmin)/(xmax-xmin))+1;     
Yp=1+Ny*(pos(:,2)-ymin)/(ymax-ymin);
Ix=floor(Xp); Iy=floor(Yp); %nz = ones(Nb,1);  
Ib=Nx*Iy+Ix-Nx; Ibmax = Nx*Ny; Nmax = max(Ib); % Ib is the box number of the particles

%%
%pos2 = [pos Ib]; [Ibs,pindx]=sort(Ib);         
%[Ibs2,pindx2]=unique(Ibs,'first');
%S_box = [Ibs pindx];   
% Ibs is the sorted box numbers, pindx is the indx of the box number in the
% original Ib matrix, e.g., if Ib is [2;9;4;6],Ibs=[2;4;6;9], pindx=[1;3;4;2]
% In this version, no need to sort

%%
count = 1; count2 = 1; % for counting and comparing number of contacts determined from both algo's

%figure(1)

tic

for ind1 = 1:Np
  i2 = Ib(ind1); % ind4 = 1;
    %scatter(pos(ind1,1),pos(ind1,2));
     %  hold on
  Ib_nei = ([i2-1,i2+1,i2,i2+Nx,i2+Nx-1,i2+Nx+1,i2-Nx-1,i2-Nx+1,i2-Nx]); % Id's of 9 neighbours
    for ind2 = 1:9
        Nb2 = find(Ib==Ib_nei(ind2)); % Nb2(:) is indx of neighbors in original pos matrix
        N_Nb2 = numel(Nb2);           % N_Nb2 are total number of contacts for ind1 particle
      if N_Nb2>0
          for ind3 = 1:N_Nb2
              if Nb2(ind3)>ind1  % this is to avoid duplicate searching
                %Nb2 = pindx(Ib_nei(ind2));  % indx of actual postions for neighbours
                p1 = pos(ind1,:);        p2 = pos(Nb2(ind3),:);
                r1 = 0.002;  r2 = 0.002;  %rest length
                dx   = p2(1)-p1(1);        dy   = p2(2)-p1(2);
                distSq = dx.*dx+dy.*dy;             %distance^2
                rSq    = (r1+r2);               %interaction distance squared
                rSq    = rSq.*rSq;
                    if distSq<rSq
                        count = count+1;              % for comparing no of contacts
%                        contact(ind1,ind4)=Nb2(ind3); % determine contact pairs for comparing
 %                       ind4 = ind4+1;
%                         plot(pos(ind1,1),pos(ind1,1));
%                         hold on
                      %  scatter(pos(Nb2(ind3),1),pos(Nb2(ind3),2));
                    %hold on
                    f1 = 1; force(ind1,1) = force(ind1,1) + f1;     
                    force(Nb2(ind3),1) = force(Nb2(ind3),1) - f1; % complementary force
                    end
               end
            end
         end 
         %end
    end
end
toc
count

%%
% Con_nei = [1:1:Nx ; j*Nx,j*Nx+1]; %% for walls grid 
% if ind1<=Nx || mod(Ib,Nx)<=1  % Implement search for contacts only if particles are along boundary of wall
    %    for ind1 = 1:Np 

%% below is original algorithm by searching the paricle against all other particles
% cost - O(N^2))
%figure(2)
tic
for ind1 = 1:Np %loop over all free particles
 %   ind5 = 1;  % 
    %scatter(pos(ind1,1),pos(ind1,2));
    %hold on
    for ind2 = ind1+1:Np %loop over 2nd half of free particles
        
        p1 = pos(ind1,:); %position
        p2 = pos(ind2,:);
        r1 = 0.002;
        r2 = 0.002;  %rest length
        dx   = p2(1)-p1(1);
        dy   = p2(2)-p1(2);
        distSq = dx.*dx+dy.*dy;             %distance^2
        rSq    = (r1+r2);               %interaction distance squared
        rSq    = rSq.*rSq;
        
        if distSq<rSq            
          count2 = count2+1;  
  %         contact2(ind1,ind5)=ind2;
      %  scatter(pos(ind2,1),pos(ind2,2));
       % hold on   
       f2 = 1; force2(ind1,1) = force2(ind1,1) + f2;
       force2(ind2,1) = force2(ind2,1) - f2;
   %     ind5 = ind5+1;
        end
        
    end


end
toc
count2
%%
% 
% %         if Ib_nei(ind2)==0 || Ib_nei(ind2)<0 ||Ib_nei(ind2)>Ibmax
% %             Ib_nei(ind2)=1;
% %             Nb2(ind2) = pindx(Ib_nei); % indx of actual postions for neighbouring elements
% %         else
% %             Nb2(ind2) = pindx(Ib_nei); % indx of actual postions for neighbouring elements
% %         end