function results5(dat, datfile, roi, oename, sinfo, filename)
% results: Matlab function saves results computed by sfrmat5.m
%  Data stored is in an Excel format file.
%
%  Usage: status =  results(dat, datfile, roi, memo, filename)
%  dat =      spatial frequency, (srf array) in either a 2 column or
%             4 column (frq, r, g, b) form
%  datfile =  image file used for input data
%  roi =      pixel coordinates that define the Region of Interest in
%             datfile
%  oename =    name of OECF file applied (or 'none').
%  sinfo =     structure with various fields for sfrmat4 results
%  filename = optional name of file where results are to be saved. If
%             this is not supplied, the user is prompted by a dialog
%             window.
% 19 July 2019
%
% Copyright (c) Peter D. Burns 2019
%

if nargin<5
    filename = '';
end

if nargin<3
    disp('* Error in results function, at least 3 arguments needed *');
    return;
end

pdatfile = datfile;

disp(['* Writing results to file: ',filename])

[~, cols] = size(dat);

edgedat = sinfo.edgedat;
misreg = edgedat(:,end)';
samp = sinfo.samp;
nbin = samp(4);
se = sinfo.se;
sfr1050 = sinfo.sfr1050;
npol = sinfo.npol;
sunit = sinfo.sunit;
funit = sinfo.funit;

% Slope in degrees
dslope = edgedat(:,2);

line1 = 'Output from Matlab function sfrmat5.m';
line3 = 'Analysis:  Spatial Frequency Response';
line4 = datetime('today', 'Format', 'dd-MM-yyyy');
line5 = {'Image/data evaluated', pdatfile};
line6 = {'This output file', filename};
line7 = {'Selected region' ,['(', num2str(roi(1)), ', ',num2str(roi(2)),...
    '), to (',num2str(roi(3)),', ',num2str(roi(4)),')']};
line8 = {'OECF applied', oename};

if strcmp(sunit, 'mm')
    line8a = {['Sampling (Image), ',sunit], num2str(samp(1),3),'PPI',num2str(round(25.4/samp(1)),3) };
else
    line8a = {['Sampling (Image), ',sunit], num2str(samp(1),3)};
end
line8b = {'Edge fit order', num2str(npol,2)};
line8c = {'Binning', num2str(nbin, 2)};
line8d = {['ESF Sampling, ',sunit], num2str(samp(3),3)};
if cols>2
    line8e = {'  ','Red', 'Green','Blue', 'Lum'};
    line8f = {'Slope, degrees',num2str(dslope(1),3),num2str(dslope(2),3),...
        num2str(dslope(3),3),num2str(dslope(4),3)};
    line9  = {'Color Misreg, pixels',num2str(misreg(1),2),...
        num2str(misreg(2),2),num2str(misreg(3),2),num2str(misreg(4),2)};
    line10 = {'Sampling Efficiency',num2str(se(1,1),3),...
        num2str(se(1,2),3),num2str(se(1,3),3),num2str(se(1,4),3)};
    line11 = {['SFR50, ',funit],num2str(sfr1050(2,1),3),...
        num2str(sfr1050(2,2),3),num2str(sfr1050(2,3),3),num2str(sfr1050(2,4),3)};

    line12 = {['SFR10, ',funit],num2str(sfr1050(1,1),3),...
        num2str(sfr1050(1,2),3),num2str(sfr1050(1,3),3),num2str(sfr1050(1,4),3)};
else
    line8e = ' ';
    line8f = {'Slope, degrees',num2str(dslope(1),3)};
    line9 = ' ';
    line10 = {'Sampling Efficiency',num2str(se(1),3)};
    line11 = {['SFR50, ',funit],num2str(sfr1050(2),3)};
    line12 = {['SFR10, ',funit],num2str(sfr1050(1),3)};
end
if cols<4
    line13 = {['Frequency, ',funit],  'SFR'};
else
    line13 = {['Frequency, ',funit],'SFR-r','SFR-g','SFR-b','Lum'};
end

if exist(filename, 'file')==2
    disp(['Overwriting: ', filename])
    delete(filename);
end

writematrix(line1, filename, 'Sheet', 1, 'Range', 'B1');
writematrix(line3, filename, 'Sheet', 1, 'Range', 'B2');
writematrix(line4, filename, 'Sheet', 1, 'Range', 'B4');
writecell(line5, filename, 'Sheet', 1, 'Range', 'B5');
writecell(line6, filename, 'Sheet', 1, 'Range', 'B6');
writecell(line7, filename, 'Sheet', 1, 'Range', 'B7');
writecell(line8, filename, 'Sheet', 1, 'Range', 'B8');
writecell(line8a, filename, 'Sheet', 1, 'Range', 'B9');
writecell(line8b, filename, 'Sheet', 1, 'Range', 'B10');
writecell(line8c, filename, 'Sheet', 1, 'Range', 'B11');
writecell(line8d, filename, 'Sheet', 1, 'Range', 'B12');
writecell(line8e, filename, 'Sheet', 1, 'Range', 'B13');
writecell(line8f, filename, 'Sheet', 1, 'Range', 'B14');
writecell(line9, filename, 'Sheet', 1, 'Range', 'B15');
writecell(line10, filename, 'Sheet', 1, 'Range', 'B16');
writecell(line11, filename, 'Sheet', 1, 'Range', 'B17');
writecell(line12, filename, 'Sheet', 1, 'Range', 'B18');
writecell(line13, filename, 'Sheet', 1, 'Range', 'B20');
writematrix(round(dat,3), filename, 'Sheet', 1, 'Range', 'B21');
pause(2);
