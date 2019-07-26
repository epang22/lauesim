%MAKECSV
% Make .csv input file for laue_simulatepattern3.m from VESTA reflections .txt file
% only works for tetragonal symmetry
% top of VESTA .txt file should look like this:
%    h    k    l      d (Å)      F(real)      F(imag)          |F|         2?          I    M ID(?) Phase
%    1    0    1   2.980385    66.512451     5.531018       66.742   11.55414  100.00000    8     1     1

clear

%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputfile = 'ICSD41579_ZrO2-12CeO2_tetragonal_Yashima1995_reflectionlist.txt';  % can also specify full path
outputfile = 'ICSD41579_ZrO2-12CeO2_tetragonal_Yashima1995_reflectionlist.csv';
minF = 0.01;    % Minimum structure factor |F| to be outputted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% DO NOT EDIT BELOW THIS LINE (unless you are doing code development) %


% read in data
fileidin = fopen(inputfile);
cw = [4,5,5,11,13,13,13,11,11,5,6,6];     % specify field widths
fmt = sprintf('%%%df',cw);  % specify floating point fields
data = cell2mat(textscan(fileidin,fmt,'headerLines', 1));
fclose(fileidin);


% Loop through hkl, calculate multiplicity, output to file
fileidout = fopen(outputfile,'w');     % Open file
totalhkl = 0;   % initialize running count of total hkl outputted
for ii=1:size(data,1)
    counter = 0;
    if data(ii,7)>minF
        h = data(ii,1);
        k = data(ii,2);
        l = data(ii,3);
        d = data(ii,4);
        F = data(ii,7);
        M = data(ii,10);
        
        % figure out multiplicity and print h,k,l,d,|F|
        if h==0 && k==0 && l~=0     % (0,0,l)
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',0,0,l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',0,0,-l,d,F);
            if M~=2
                fprintf('ii=%.0f, (0,0,l) loop run, M actual = %.0f\n',ii,M);
            end
            counter = 1;    % change to 1 when loop runs (debug purposes only)
            totalhkl = totalhkl + 2;  % running count of total hkl outputted
        elseif h==0 && k~=0 && l==0     % (0,k,0)
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',0,k,0,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',0,-k,0,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',k,0,0,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',-k,0,0,d,F);
            if M~=4
                fprintf('ii=%.0f, (0,k,0) loop run, M actual = %.0f\n',ii,M);
            end
            counter = 1;    % change to 1 when loop runs (debug purposes only)
            totalhkl = totalhkl + 4;  % running count of total hkl outputted
        elseif h~=0 && k==0 && l==0     % (h,0,0)
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',h,0,0,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',-h,0,0,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',0,h,0,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',0,-h,0,d,F);
            if M~=4
                fprintf('ii=%.0f, (h,0,0) loop run, M actual = %.0f\n',ii,M);
            end
            counter = 1;    % change to 1 when loop runs (debug purposes only)
            totalhkl = totalhkl + 4;  % running count of total hkl outputted
        elseif h==0 && k~=0 && l~=0     % (0,k,l)
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',0,k,l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',0,-k,l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',k,0,l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',-k,0,l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',0,k,-l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',0,-k,-l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',k,0,-l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',-k,0,-l,d,F);
            if M~=8
                fprintf('ii=%.0f, (0,k,l) loop run, M actual = %.0f\n',ii,M);
            end
            counter = 1;    % change to 1 when loop runs (debug purposes only)
            totalhkl = totalhkl + 8;  % running count of total hkl outputted
        elseif h~=0 && k==0 && l~=0     % (h,0,l)
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',h,0,l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',-h,0,l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',0,h,l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',0,-h,l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',h,0,-l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',-h,0,-l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',0,h,-l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',0,-h,-l,d,F);
            if M~=8
                fprintf('ii=%.0f, (h,0,l) loop run, M actual = %.0f\n',ii,M);
            end
            counter = 1;    % change to 1 when loop runs (debug purposes only)
            totalhkl = totalhkl + 8;  % running count of total hkl outputted
        elseif h~=0 && k~=0 && l==0 && h~=k     % (h,k,0)
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',h,k,0,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',-h,k,0,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',h,-k,0,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',-h,-k,0,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',k,h,0,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',-k,h,0,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',k,-h,0,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',-k,-h,0,d,F);
            if M~=8
                fprintf('ii=%.0f, (h,k,0) loop run, M actual = %.0f\n',ii,M);
            end
            counter = 1;    % change to 1 when loop runs (debug purposes only)
            totalhkl = totalhkl + 8;  % running count of total hkl outputted
        elseif h~=0 && k~=0 && l==0 && h==k     % (h,h,0)
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',h,h,0,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',-h,h,0,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',h,-h,0,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',-h,-h,0,d,F);
            if M~=4
                fprintf('ii=%.0f, (h,h,0) loop run, M actual = %.0f\n',ii,M);
            end
            counter = 1;    % change to 1 when loop runs (debug purposes only)
            totalhkl = totalhkl + 4;  % running count of total hkl outputted
        elseif h~=0 && k~=0 && l~=0 && h~=k     % (h,k,l)
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',h,k,l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',-h,k,l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',h,-k,l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',-h,-k,l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',k,h,l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',-k,h,l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',k,-h,l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',-k,-h,l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',h,k,-l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',-h,k,-l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',h,-k,-l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',-h,-k,-l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',k,h,-l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',-k,h,-l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',k,-h,-l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',-k,-h,-l,d,F);
            if M~=16
                fprintf('ii=%.0f, (h,k,l) loop run, M actual = %.0f\n',ii,M);
            end
            counter = 1;    % change to 1 when loop runs (debug purposes only)
            totalhkl = totalhkl + 16;  % running count of total hkl outputted
        elseif h~=0 && k~=0 && l~=0 && h==k     % (h,h,l)
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',h,h,l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',-h,h,l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',h,-h,l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',-h,-h,l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',h,h,-l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',-h,h,-l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',h,-h,-l,d,F);
            fprintf(fileidout,'%.0f,%.0f,%.0f,%.6f,%.6f\n',-h,-h,-l,d,F);
            if M~=8
                fprintf('ii=%.0f, (h,h,l) loop run, M actual = %.0f\n',ii,M);
            end
            counter = 1;    % change to 1 when loop runs (debug purposes only)
            totalhkl = totalhkl + 8;  % running count of total hkl outputted
        end

        % print if no loop run (debug purposes only)
        if counter==0
            fprintf('ii=%.0f, no loop run\n',ii);
        end
    end
end

% count how many hkl should have been outputted
totalhklM = 0;   % initialize running count of total hkl based on multiplicity M in .txt input file
for ii=1:size(data,1)
    if data(ii,7)>minF
        totalhklM = totalhklM + data(ii,10);
    end
end

% print hkl stats
fprintf('# hkl outputted: %.0f\n',totalhkl);
fprintf('# hkl based on multiplicity: %.0f\n',totalhklM);
fprintf('^These numbers should be equal. If not, the calculation may be incorrect.\n');

fclose('all');


