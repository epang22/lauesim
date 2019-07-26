%SIMULATEPATTERN_PSUEDOSYMCOMPARE2
% Compare patterns for 3 pseudosymmetry variants on same plot
% Simulate back-reflection Laue pattern
% Reads in list of hkl and |F| from .csv file
% Using coordinate system of NorthStar:
    % incident beam propagates in -z direction
    % +y direction is up in lab coordinates
    % +x is to the left (when viewing the detector from the sample)
    % (0,0,0) is where beam hits sample
    % (0,0,L) is center of detector
    % For U=identity, unit cell alignment to lab coordinates: +a=+x, +b=+y, +c=+z
    % Patterns are displayed as viewed from the sample

clear
    
%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify B
a = 3.61;
c = 5.2;
B = [1/a 0 0;0 1/a 0;0 0 1/c];

% Define orientation (choose 1 option only, comment out other blocks)
%     %%% Option 1: Input U %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     U = eye(3);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Option 2: Input hkl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hkl1 = [1;2;0];     % hkl you want pattern centered around
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     %%% Option 3: Input UB directly %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     UB = [	0.180496,	-0.203529,	0.036278;
%             0.097182,	0.019776,	-0.179560;
%             0.186307,	0.186865,	0.058516];
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specify pseudosymmetry rotation
axis1 = [1 1 0];        % Pseudocubic a-axis (note that tetragonal unit cell here is primitive)
axis2 = [1 -1 0];       % Pseudocubic a-axis (note that tetragonal unit cell here is primitive)
pseudoangle = 90;      % in degrees

% some diffraction parameters
L = 200;    % Detector distance (mm)
sizex = 300;    % Detector width (mm)
sizey = 300;    % Detector height (mm)
Emax = 20;      % Maximum energy (keV)
Rcenter = 10;   % radius (mm) of center hole

% some plotting options
minFplot = 15;   % only plot spot if |F| greater than this value
Fexp = 1;     % Raise |F| to this power for marker size and color

labelhkl = 1;           % 1 if you want hkl labeled on pattern, 0 otherwise
labeloffsety = -3;      % vertical position of hkl labels relative to spots
maxhkl = 1000;           % maximum h^2+k^2+l^2 to label
minFlabel = 15;          % label hkl and |F| only if |F| greater than this value

labelF = 0;     % 1 if you want to label value of |F|, 0 otherwise


% .csv file containing list of hkl (from WinLaue)
filename = 'ICSD41579_ZrO2-12CeO2_tetragonal_Yashima1995_reflectionlist.csv';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% DO NOT EDIT BELOW THIS LINE (unless you are doing code development) %


% figure out which option specified and calculate UB
if exist('U','var')     %%% Option 1
    UB = U*B;
elseif exist('hkl1','var')  %%% Option 2
    if hkl1(1)==0 && hkl1(2)==0   % if hkl is 001,00-1,002, etc., method doesn't work. Instead, manually specify U.
        if hkl1(3)>0    % hkl=001, no rotation needed
            U = eye(3);
        elseif hkl1(3)<0    % hkl=00-1, 180deg rotation needed
            U = [1 0 0;0 -1 0;0 0 -1];     % Note that I arbitarily chose to rotate 180deg about x instead of y
        else
            error('hkl cannot be (000).');
        end
    else    % if hkl is not 001
        U = eye(3); UB = U*B;
        n0 = UB*[0;0;1]; n0 = n0/norm(n0);      % calc unit vector of plane normal
        n1 = UB*hkl1; n1 = n1/norm(n1);

        axisdiffrot = cross(n1,n0);       % axis of rotation
        if norm(axisdiffrot)>1
            if norm(axisdiffrot)-1 < 1e-8
                angle = asin(1);    % numerical errors sometimes make norm>1, which gives imaginary result
            end
        else
            angle = asin(norm(axisdiffrot));   % angle of rotation (radians)
        end
        U = axang2rotm([axisdiffrot(1) axisdiffrot(2) axisdiffrot(3) angle]);    % rotation matrix for specified hkl
    end
    UB = U*B;
elseif exist('UB','var')    %%% Option 3
    % Don't need to do anything, UB already supplied
else
    error('Invalid orientation input.');
end

% Calculate UB for pseudosymmetry variants
Rrot1 = axang2rotm([axis1 pseudoangle*pi/180]);  % Matrix for specified rotation 1
U = UB/B;
U1=U*Rrot1;
UB1 = U1*B;

Rrot2 = axang2rotm([axis2 pseudoangle*pi/180]);  % Matrix for specified rotation 1
U2=U*Rrot2;
UB2 = U2*B;


% Read in hkl data
hkl = dlmread(filename,',',1,0);


% calculations
Lhat = [0; 0; 1];  % normal vector from sample to detector (-incident beam)
detectorside = 1;       % 1=pattern viewed from sample side, 0=other side
lambdamin = 1e10*4.135668e-15*2.998e8/(1000*Emax);  % minimum wavelength radiation (A)

% initialize to values that are outside of detector
spots = max(sizex,sizey)*ones(size(hkl,1),2);     % coordinates of diffracted spots on detector (original variant)
spots1 = spots;     % coordinates of diffracted spots on detector (pseudosym variant 1)
spots2 = spots;     % coordinates of diffracted spots on detector (pseudosym variant 1)
for ii=1:size(hkl,1)
    %%% Original variant
    % find normal vector for hkl plane
    n = UB*[hkl(ii,1); hkl(ii,2); hkl(ii,3)];   % normal vector
    n = n/norm(n);      % normalize to unit vector
    
    % calculate diffracted beam
    axisdiffrot = cross(Lhat,n);       % axis of rotation
    
    if norm(axisdiffrot)>1
        if norm(axisdiffrot)-1 < 1e-8
            angle = asin(1);    % numerical errors sometimes make norm>1, which gives imaginary result
        end
    else
        angle = asin(norm(axisdiffrot));   % angle of rotation (radians)
    end
    
    Udiff = axang2rotm([axisdiffrot(1) axisdiffrot(2) axisdiffrot(3) 2*angle]);    % rotation matrix to go from -incident beam to diffracted beam
    q = Udiff*Lhat;             % diffracted beam (unit vector)
    
    if q(3) > 0     % if diffracted beam is pointing towards detector
        k = L/q(3); q = k*q;       % scale diffracted beam vector length so that it reaches detector
        
        % determine whether or not this spot is diffracted by incident wavelengths
        if dot(Lhat,n)>1
            if dot(Lhat,n)-1 < 1e-8
                theta = 0;  % numerical errors sometimes make dp>1, which gives imaginary result
            end
        else
            theta = pi/2-acos(dot(Lhat,n));
        end
        lambda = 2*hkl(ii,4)*sin(theta);    % required wavelength for this spot to be diffracted (A)
        
        % save coordinates of diffracted spot on detector if diffracted
        if lambda>lambdamin
            if detectorside==0
                spots(ii,1) = q(1);    % flip pattern horizontally
            else
                spots(ii,1) = -q(1);
            end
            spots(ii,2) = q(2);
        end
    end
    
    %%% Pseudosymmetry variant 1
    % find normal vector for hkl plane
    n = UB1*[hkl(ii,1); hkl(ii,2); hkl(ii,3)];   % normal vector
    n = n/norm(n);      % normalize to unit vector
    
    % calculate diffracted beam
    axisdiffrot = cross(Lhat,n);       % axis of rotation
    
    if norm(axisdiffrot)>1
        if norm(axisdiffrot)-1 < 1e-8
            angle = asin(1);    % numerical errors sometimes make norm>1, which gives imaginary result
        end
    else
        angle = asin(norm(axisdiffrot));   % angle of rotation (radians)
    end
    
    Udiff = axang2rotm([axisdiffrot(1) axisdiffrot(2) axisdiffrot(3) 2*angle]);    % rotation matrix to go from -incident beam to diffracted beam
    q = Udiff*Lhat;             % diffracted beam (unit vector)
    
    if q(3) > 0     % if diffracted beam is pointing towards detector
        k = L/q(3); q = k*q;       % scale diffracted beam vector length so that it reaches detector
        
        % determine whether or not this spot is diffracted by incident wavelengths
        if dot(Lhat,n)>1
            if dot(Lhat,n)-1 < 1e-8
                theta = 0;  % numerical errors sometimes make dp>1, which gives imaginary result
            end
        else
            theta = pi/2-acos(dot(Lhat,n));
        end
        lambda = 2*hkl(ii,4)*sin(theta);    % required wavelength for this spot to be diffracted (A)
        
        % save coordinates of diffracted spot on detector if diffracted
        if lambda>lambdamin
            if detectorside==0
                spots1(ii,1) = q(1);    % flip pattern horizontally
            else
                spots1(ii,1) = -q(1);
            end
            spots1(ii,2) = q(2);
        end
    end
    
    %%% Pseudosymmetry variant 2
    % find normal vector for hkl plane
    n = UB2*[hkl(ii,1); hkl(ii,2); hkl(ii,3)];   % normal vector
    n = n/norm(n);      % normalize to unit vector
    
    % calculate diffracted beam
    axisdiffrot = cross(Lhat,n);       % axis of rotation
    
    if norm(axisdiffrot)>1
        if norm(axisdiffrot)-1 < 1e-8
            angle = asin(1);    % numerical errors sometimes make norm>1, which gives imaginary result
        end
    else
        angle = asin(norm(axisdiffrot));   % angle of rotation (radians)
    end
    
    Udiff = axang2rotm([axisdiffrot(1) axisdiffrot(2) axisdiffrot(3) 2*angle]);    % rotation matrix to go from -incident beam to diffracted beam
    q = Udiff*Lhat;             % diffracted beam (unit vector)
    
    if q(3) > 0     % if diffracted beam is pointing towards detector
        k = L/q(3); q = k*q;       % scale diffracted beam vector length so that it reaches detector
        
        % determine whether or not this spot is diffracted by incident wavelengths
        if dot(Lhat,n)>1
            if dot(Lhat,n)-1 < 1e-8
                theta = 0;  % numerical errors sometimes make dp>1, which gives imaginary result
            end
        else
            theta = pi/2-acos(dot(Lhat,n));
        end
        lambda = 2*hkl(ii,4)*sin(theta);    % required wavelength for this spot to be diffracted (A)
        
        % save coordinates of diffracted spot on detector if diffracted
        if lambda>lambdamin
            if detectorside==0
                spots2(ii,1) = q(1);    % flip pattern horizontally
            else
                spots2(ii,1) = -q(1);
            end
            spots2(ii,2) = q(2);
        end
    end
end


% plot with intensity and labels
    
    % calculate structure factor
    F = hkl(:,5);
    Fmax = max(F);

figure('Position', [200 200 550 550]);
for ii=1:size(hkl,1)
    % Plot
    %%% Original variant
    if spots(ii,1)>(-sizex/2) && spots(ii,1)<(sizex/2) && spots(ii,2)>(-sizey/2) && spots(ii,2)<(sizey/2)   % if spot falls within detector
        spotcolor = [1-(F(ii)/Fmax)^Fexp 1-(F(ii)/Fmax)^Fexp 1-(F(ii)/Fmax)^Fexp];  % figure out spot color (proportional to F^2)
        spotsize = 20*(F(ii)/Fmax)^Fexp;   % figure out spot size (proportional to |F|^Fexp)
        
        % plot
        if spots(ii,1)^2+spots(ii,2)^2 > Rcenter^2      % only plot if outside of center hole
            if hkl(ii,5) > minFplot    % only plot if |F| larger than a certain value
                plot(spots(ii,1),spots(ii,2),'o','MarkerSize',spotsize,'MarkerEdgeColor',0.9*spotcolor,'MarkerFaceColor',spotcolor); hold on;
            end
        end
    end
    
    %%% psuedosym variant 1
    if spots1(ii,1)>(-sizex/2) && spots1(ii,1)<(sizex/2) && spots1(ii,2)>(-sizey/2) && spots1(ii,2)<(sizey/2)   % if spot falls within detector
        spotcolor = [1 1-(F(ii)/Fmax)^Fexp 1-(F(ii)/Fmax)^Fexp];    % figure out spot color (proportional to F^2)
        spotsize = 20*(F(ii)/Fmax)^Fexp;   % figure out spot size (proportional to |F|^Fexp)
        
        % plot
        if spots1(ii,1)^2+spots1(ii,2)^2 > Rcenter^2      % only plot if outside of center hole
            if hkl(ii,5) > minFplot    % only plot if |F| larger than a certain value
                plot(spots1(ii,1),spots1(ii,2),'o','MarkerSize',spotsize,'MarkerEdgeColor',0.9*spotcolor,'MarkerFaceColor',spotcolor); hold on;
            end
        end
    end
    
    %%% psuedosym variant 2
    if spots2(ii,1)>(-sizex/2) && spots2(ii,1)<(sizex/2) && spots2(ii,2)>(-sizey/2) && spots2(ii,2)<(sizey/2)   % if spot falls within detector
        spotcolor = [1-(F(ii)/Fmax)^Fexp 1-(F(ii)/Fmax)^Fexp 1];    % figure out spot color (proportional to F^2)
        spotsize = 20*(F(ii)/Fmax)^Fexp;   % figure out spot size (proportional to |F|^Fexp)
        
        % plot
        if spots2(ii,1)^2+spots2(ii,2)^2 > Rcenter^2      % only plot if outside of center hole
            if hkl(ii,5) > minFplot    % only plot if |F| larger than a certain value
                plot(spots2(ii,1),spots2(ii,2),'o','MarkerSize',spotsize,'MarkerEdgeColor',0.9*spotcolor,'MarkerFaceColor',spotcolor); hold on;
            end
        end
    end
    
    
    % Label
    %%% Original variant
    if spots(ii,1)>(-sizex/2) && spots(ii,1)<(sizex/2) && spots(ii,2)>(-sizey/2) && spots(ii,2)<(sizey/2)   % if spot falls within detector
        % label hkl
        if labelhkl==1      % label if option specified
            if hkl(ii,1)^2+hkl(ii,2)^2+hkl(ii,3)^2 < maxhkl    % only label hkl if h^2+k^2+l^2 below specified
                if hkl(ii,5) > minFlabel    % only label hkl if |F| larger than a certain value
                    if spots(ii,1)^2+spots(ii,2)^2 > Rcenter^2      % only label hkl if outside of center hole
                        if hkl(ii,1)<0 && hkl(ii,2)<0 && hkl(ii,3)<0
                            s1 = sprintf('$\\bar{%s}\\bar{%s}\\bar{%s}$',num2str(-hkl(ii,1)),num2str(-hkl(ii,2)),num2str(-hkl(ii,3)));
                            text(spots(ii,1),spots(ii,2)+labeloffsety-spotsize/2,s1, 'Interpreter', 'LaTeX','FontSize',8,'HorizontalAlignment','center');
                        elseif hkl(ii,1)>=0 && hkl(ii,2)<0 && hkl(ii,3)<0
                            s1 = sprintf('$%s\\bar{%s}\\bar{%s}$',num2str(hkl(ii,1)),num2str(-hkl(ii,2)),num2str(-hkl(ii,3)));
                            text(spots(ii,1),spots(ii,2)+labeloffsety-spotsize/2,s1, 'Interpreter', 'LaTeX','FontSize',8,'HorizontalAlignment','center');
                        elseif hkl(ii,1)<0 && hkl(ii,2)>=0 && hkl(ii,3)<0
                            s1 = sprintf('$\\bar{%s}%s\\bar{%s}$',num2str(-hkl(ii,1)),num2str(hkl(ii,2)),num2str(-hkl(ii,3)));
                            text(spots(ii,1),spots(ii,2)+labeloffsety-spotsize/2,s1, 'Interpreter', 'LaTeX','FontSize',8,'HorizontalAlignment','center');
                        elseif hkl(ii,1)<0 && hkl(ii,2)<0 && hkl(ii,3)>=0
                            s1 = sprintf('$\\bar{%s}\\bar{%s}%s$',num2str(-hkl(ii,1)),num2str(-hkl(ii,2)),num2str(hkl(ii,3)));
                            text(spots(ii,1),spots(ii,2)+labeloffsety-spotsize/2,s1, 'Interpreter', 'LaTeX','FontSize',8,'HorizontalAlignment','center');
                        elseif hkl(ii,1)>=0 && hkl(ii,2)>=0 && hkl(ii,3)<0
                            s1 = sprintf('$%s%s\\bar{%s}$',num2str(hkl(ii,1)),num2str(hkl(ii,2)),num2str(-hkl(ii,3)));
                            text(spots(ii,1),spots(ii,2)+labeloffsety-spotsize/2,s1, 'Interpreter', 'LaTeX','FontSize',8,'HorizontalAlignment','center');
                        elseif hkl(ii,1)>=0 && hkl(ii,2)<0 && hkl(ii,3)>=0
                            s1 = sprintf('$%s\\bar{%s}%s$',num2str(hkl(ii,1)),num2str(-hkl(ii,2)),num2str(hkl(ii,3)));
                            text(spots(ii,1),spots(ii,2)+labeloffsety-spotsize/2,s1, 'Interpreter', 'LaTeX','FontSize',8,'HorizontalAlignment','center');
                        elseif hkl(ii,1)<0 && hkl(ii,2)>=0 && hkl(ii,3)>=0
                            s1 = sprintf('$\\bar{%s}%s%s$',num2str(-hkl(ii,1)),num2str(hkl(ii,2)),num2str(hkl(ii,3)));
                            text(spots(ii,1),spots(ii,2)+labeloffsety-spotsize/2,s1, 'Interpreter', 'LaTeX','FontSize',8,'HorizontalAlignment','center');
                        else
                            s1 = [num2str(hkl(ii,1)),num2str(hkl(ii,2)),num2str(hkl(ii,3))];
                            text(spots(ii,1),spots(ii,2)+labeloffsety-spotsize/2,s1, 'Interpreter', 'LaTeX','FontSize',8,'HorizontalAlignment','center');
                        end
                        
                        if labelF == 1  % if specified, label value of structure factor |F|
                            s2 = sprintf('%.1f',F(ii));
                            text(spots(ii,1),spots(ii,2)-labeloffsety+spotsize/2,s2,'FontSize',8,'HorizontalAlignment','center');
                        end
                    end
                end
            end
        end
    end
    
    %%% psuedosym variant 1
    if spots1(ii,1)>(-sizex/2) && spots1(ii,1)<(sizex/2) && spots1(ii,2)>(-sizey/2) && spots1(ii,2)<(sizey/2)   % if spot falls within detector
        % label hkl
        if labelhkl==1      % label if option specified
            if hkl(ii,1)^2+hkl(ii,2)^2+hkl(ii,3)^2 < maxhkl    % only label hkl if h^2+k^2+l^2 below specified
                if hkl(ii,5) > minFlabel    % only label hkl if |F| larger than a certain value
                    if spots1(ii,1)^2+spots1(ii,2)^2 > Rcenter^2      % only label hkl if outside of center hole
                        if hkl(ii,1)<0 && hkl(ii,2)<0 && hkl(ii,3)<0
                            s1 = sprintf('$\\bar{%s}\\bar{%s}\\bar{%s}$',num2str(-hkl(ii,1)),num2str(-hkl(ii,2)),num2str(-hkl(ii,3)));
                            text(spots1(ii,1),spots1(ii,2)+labeloffsety-spotsize/2,s1, 'Interpreter', 'LaTeX','FontSize',8,'HorizontalAlignment','center','Color','red');
                        elseif hkl(ii,1)>=0 && hkl(ii,2)<0 && hkl(ii,3)<0
                            s1 = sprintf('$%s\\bar{%s}\\bar{%s}$',num2str(hkl(ii,1)),num2str(-hkl(ii,2)),num2str(-hkl(ii,3)));
                            text(spots1(ii,1),spots1(ii,2)+labeloffsety-spotsize/2,s1, 'Interpreter', 'LaTeX','FontSize',8,'HorizontalAlignment','center','Color','red');
                        elseif hkl(ii,1)<0 && hkl(ii,2)>=0 && hkl(ii,3)<0
                            s1 = sprintf('$\\bar{%s}%s\\bar{%s}$',num2str(-hkl(ii,1)),num2str(hkl(ii,2)),num2str(-hkl(ii,3)));
                            text(spots1(ii,1),spots1(ii,2)+labeloffsety-spotsize/2,s1, 'Interpreter', 'LaTeX','FontSize',8,'HorizontalAlignment','center','Color','red');
                        elseif hkl(ii,1)<0 && hkl(ii,2)<0 && hkl(ii,3)>=0
                            s1 = sprintf('$\\bar{%s}\\bar{%s}%s$',num2str(-hkl(ii,1)),num2str(-hkl(ii,2)),num2str(hkl(ii,3)));
                            text(spots1(ii,1),spots1(ii,2)+labeloffsety-spotsize/2,s1, 'Interpreter', 'LaTeX','FontSize',8,'HorizontalAlignment','center','Color','red');
                        elseif hkl(ii,1)>=0 && hkl(ii,2)>=0 && hkl(ii,3)<0
                            s1 = sprintf('$%s%s\\bar{%s}$',num2str(hkl(ii,1)),num2str(hkl(ii,2)),num2str(-hkl(ii,3)));
                            text(spots1(ii,1),spots1(ii,2)+labeloffsety-spotsize/2,s1, 'Interpreter', 'LaTeX','FontSize',8,'HorizontalAlignment','center','Color','red');
                        elseif hkl(ii,1)>=0 && hkl(ii,2)<0 && hkl(ii,3)>=0
                            s1 = sprintf('$%s\\bar{%s}%s$',num2str(hkl(ii,1)),num2str(-hkl(ii,2)),num2str(hkl(ii,3)));
                            text(spots1(ii,1),spots1(ii,2)+labeloffsety-spotsize/2,s1, 'Interpreter', 'LaTeX','FontSize',8,'HorizontalAlignment','center','Color','red');
                        elseif hkl(ii,1)<0 && hkl(ii,2)>=0 && hkl(ii,3)>=0
                            s1 = sprintf('$\\bar{%s}%s%s$',num2str(-hkl(ii,1)),num2str(hkl(ii,2)),num2str(hkl(ii,3)));
                            text(spots1(ii,1),spots1(ii,2)+labeloffsety-spotsize/2,s1, 'Interpreter', 'LaTeX','FontSize',8,'HorizontalAlignment','center','Color','red');
                        else
                            s1 = [num2str(hkl(ii,1)),num2str(hkl(ii,2)),num2str(hkl(ii,3))];
                            text(spots1(ii,1),spots1(ii,2)+labeloffsety-spotsize/2,s1, 'Interpreter', 'LaTeX','FontSize',8,'HorizontalAlignment','center','Color','red');
                        end
                        
                        if labelF == 1  % if specified, label value of structure factor |F|
                            s2 = sprintf('%.1f',F(ii));
                            text(spots1(ii,1),spots1(ii,2)-labeloffsety+spotsize/2,s2,'FontSize',8,'HorizontalAlignment','center','Color','red');
                        end
                    end
                end
            end
        end
    end
    
    %%% psuedosym variant 2
    if spots2(ii,1)>(-sizex/2) && spots2(ii,1)<(sizex/2) && spots2(ii,2)>(-sizey/2) && spots2(ii,2)<(sizey/2)   % if spot falls within detector
        % label hkl
        if labelhkl==1      % label if option specified
            if hkl(ii,1)^2+hkl(ii,2)^2+hkl(ii,3)^2 < maxhkl    % only label hkl if h^2+k^2+l^2 below specified
                if hkl(ii,5) > minFlabel    % only label hkl if |F| larger than a certain value
                    if spots2(ii,1)^2+spots2(ii,2)^2 > Rcenter^2      % only label hkl if outside of center hole
                        if hkl(ii,1)<0 && hkl(ii,2)<0 && hkl(ii,3)<0
                            s1 = sprintf('$\\bar{%s}\\bar{%s}\\bar{%s}$',num2str(-hkl(ii,1)),num2str(-hkl(ii,2)),num2str(-hkl(ii,3)));
                            text(spots2(ii,1),spots2(ii,2)+labeloffsety-spotsize/2,s1, 'Interpreter', 'LaTeX','FontSize',8,'HorizontalAlignment','center','Color','blue');
                        elseif hkl(ii,1)>=0 && hkl(ii,2)<0 && hkl(ii,3)<0
                            s1 = sprintf('$%s\\bar{%s}\\bar{%s}$',num2str(hkl(ii,1)),num2str(-hkl(ii,2)),num2str(-hkl(ii,3)));
                            text(spots2(ii,1),spots2(ii,2)+labeloffsety-spotsize/2,s1, 'Interpreter', 'LaTeX','FontSize',8,'HorizontalAlignment','center','Color','blue');
                        elseif hkl(ii,1)<0 && hkl(ii,2)>=0 && hkl(ii,3)<0
                            s1 = sprintf('$\\bar{%s}%s\\bar{%s}$',num2str(-hkl(ii,1)),num2str(hkl(ii,2)),num2str(-hkl(ii,3)));
                            text(spots2(ii,1),spots2(ii,2)+labeloffsety-spotsize/2,s1, 'Interpreter', 'LaTeX','FontSize',8,'HorizontalAlignment','center','Color','blue');
                        elseif hkl(ii,1)<0 && hkl(ii,2)<0 && hkl(ii,3)>=0
                            s1 = sprintf('$\\bar{%s}\\bar{%s}%s$',num2str(-hkl(ii,1)),num2str(-hkl(ii,2)),num2str(hkl(ii,3)));
                            text(spots2(ii,1),spots2(ii,2)+labeloffsety-spotsize/2,s1, 'Interpreter', 'LaTeX','FontSize',8,'HorizontalAlignment','center','Color','blue');
                        elseif hkl(ii,1)>=0 && hkl(ii,2)>=0 && hkl(ii,3)<0
                            s1 = sprintf('$%s%s\\bar{%s}$',num2str(hkl(ii,1)),num2str(hkl(ii,2)),num2str(-hkl(ii,3)));
                            text(spots2(ii,1),spots2(ii,2)+labeloffsety-spotsize/2,s1, 'Interpreter', 'LaTeX','FontSize',8,'HorizontalAlignment','center','Color','blue');
                        elseif hkl(ii,1)>=0 && hkl(ii,2)<0 && hkl(ii,3)>=0
                            s1 = sprintf('$%s\\bar{%s}%s$',num2str(hkl(ii,1)),num2str(-hkl(ii,2)),num2str(hkl(ii,3)));
                            text(spots2(ii,1),spots2(ii,2)+labeloffsety-spotsize/2,s1, 'Interpreter', 'LaTeX','FontSize',8,'HorizontalAlignment','center','Color','blue');
                        elseif hkl(ii,1)<0 && hkl(ii,2)>=0 && hkl(ii,3)>=0
                            s1 = sprintf('$\\bar{%s}%s%s$',num2str(-hkl(ii,1)),num2str(hkl(ii,2)),num2str(hkl(ii,3)));
                            text(spots2(ii,1),spots2(ii,2)+labeloffsety-spotsize/2,s1, 'Interpreter', 'LaTeX','FontSize',8,'HorizontalAlignment','center','Color','blue');
                        else
                            s1 = [num2str(hkl(ii,1)),num2str(hkl(ii,2)),num2str(hkl(ii,3))];
                            text(spots2(ii,1),spots2(ii,2)+labeloffsety-spotsize/2,s1, 'Interpreter', 'LaTeX','FontSize',8,'HorizontalAlignment','center','Color','blue');
                        end
                        
                        if labelF == 1  % if specified, label value of structure factor |F|
                            s2 = sprintf('%.1f',F(ii));
                            text(spots2(ii,1),spots2(ii,2)-labeloffsety+spotsize/2,s2,'FontSize',8,'HorizontalAlignment','center','Color','blue');
                        end
                    end
                end
            end
        end
    end
end

% plot center circle
pos = [-Rcenter -Rcenter 2*Rcenter 2*Rcenter]; 
rectangle('Position',pos,'Curvature',[1 1],'FaceColor','k','EdgeColor','k');

% axis options
axis square;
axis([-sizex/2 sizex/2 -sizey/2 sizey/2]);
h = gca; h.XAxis.Visible = 'off'; h.YAxis.Visible = 'off';
set(gca,'position',[0 0 1 1],'units','normalized');
