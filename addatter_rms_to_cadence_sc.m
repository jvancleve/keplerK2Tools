function addatter_rms_to_cadence_sc(campaign, tweakBad, threshold, persistence, singleThreshold)

%NAME:
%  addatter_rms_to_cadence_sc.m
%PURPOSE:
%  Calculate COARSE_POINT cadences for Data Anomaly Flags with
%    persistent low threshold and non-persistent high-threshold exclusion
%USAGE:
%  Run in directory containing data.  
%CALLED BY:
%  None
%CALLS
%  parse_tcad_time.m
%  read_adatterr.m to generate inputs from TCAD text files
%INPUTS
%  Arguments
%  1.  campaign as a string
%  2.  vector tweakBad containing cadences before attitude tweak or otherwise
%    known to be bad evenif ADATTERR is good.  Recommend tweakBad if tweak
%    > 4" = 1.0 science pixels
%  3.  threshold, typically 1.5 pixels
%  4.  persistence above threshold, typically 4 cadences
%  5.  singleThreshold, mark cadence bad even if only one cadence.
%  Typically 2.5 pixels
%  Files
%  1.  c15_adatterr.mat containg adatterr[x,y,z] and mjd for ADCS TLM
%  2.  Equation of Time
%OUTPUTS
%  pseudo-xml file listing excluded cadences and text to paste into DAF XML
%    file
%NOTES
%TO-DO list
%REVISION HISTORY:
%Engineer          Org     Date        Description
%J. Van Cleve     SETI     01/02/2018  Created
%                          01/17/2018  Commented and functionalized 
%                                      tweakBad array for attitude tweaks
%                                      campaign is a string
%                                      protect against first/last cadence
%                                        bad
%J. Van Cleve     SETI      01/30/2018 Identify cadences w/o S/C TLM
%                           07/31/2018 fix nCadences bug (no effect on xml)
%                           10/31/2018 LONG --> SHORT in output 
close all
load(['c' campaign '_adatterr.mat'])
%convert to arcsec
x = 206000*adatterrmx;
y = 206000*adatterrmy;
z = 206000*adatterrmz;
fovRadius = 8/57;
%Ball RMS definition
r = sqrt(x.^2 + y.^2 + z.^2);
%science pixel motion
rho = sqrt(fovRadius^2*x.^2 + y.^2 + z.^2)/3.98;
baseMjd = floor(mjds(1)/100)*100;
figure('Position',[1 1 1200 500])
plot(mjds - baseMjd, r)
title(['C' campaign ' RMS Attitude Error from S/C ADDATTERR TLM'],'FontSize',14)
xlabel(['MJD - ' num2str(baseMjd)],'FontSize',12)
ylabel('Arcsec','FontSize',12)
grid
saveas(gcf,['c' campaign '_addatterr_mjd_sc.fig'])
saveas(gcf,['c' campaign '_addatterr_mjd_sc.png'])
%%
figure('Position',[1 1 1200 500])
plot(mjds - baseMjd, rho)
title(['C' campaign ' RMS Pixel Motion Error from S/C ADDATTERR TLM'],'FontSize',14)
xlabel(['MJD - ' num2str(baseMjd)],'FontSize',12)
ylabel('Science PIxels','FontSize',12)
grid on
saveas(gcf,['c' campaign '_adatterr2pix_mjd_sc.fig'])
saveas(gcf,['c' campaign '_adatterr2pix_mjd_sc.png'])
%%
% convert MJD to cadence using Equation of Time
load(['c' campaign '_eot_sc.mat'])
excludeBufferSize = 1;
cadenceFromEngReal = polyval(p_mjd_to_cad,mjds);
cadenceFromEngInteger = round(cadenceFromEngReal);
listOfCadences = unique(cadenceFromEngInteger);
minCadence = min(listOfCadences);
maxCadence = max(listOfCadences);
%fig bug here
nCadences = maxCadence - minCadence + 1;
cadencesWithNoEngData = setdiff(minCadence:maxCadence, listOfCadences);

%%
daFlags = zeros(length(listOfCadences),1);
maxErr = zeros(length(listOfCadences),1);
for i = 1:length(listOfCadences)
    if ~mod(i,1000)
        display(num2str(i))
    end
%absolute value sign    
    maxErr(i) = max(abs(rho(cadenceFromEngInteger == listOfCadences(i))));
end
save(['c' campaign '_max_rms_sc.mat'], 'cadencesWithNoEngData', 'listOfCadences','maxErr','campaign','excludeBufferSize')
%%
baseCadence = floor(listOfCadences(1)/1000)*1000;
momentumDumps = dlmread(['c' campaign '_isMmntmDmp_SC.txt']);
adcsNotFinePoint = dlmread(['c' campaign '_isNotFinePoint_SC.txt']);
adcsNfpNotDump = setdiff(adcsNotFinePoint,momentumDumps);
maxErr(ismember(listOfCadences, momentumDumps)) = 0;
figure('Position',[1 1 1200 500])
plot(listOfCadences - baseCadence, maxErr)
title(['C' campaign ' Max RMS Error over Cadence from S/C TLM excluding resats, threshold: ' num2str(threshold)  ' persistence: ' ...
    num2str(persistence)],'FontSize',14)
xlabel(['SC Cadence - ' num2str(baseCadence)],'FontSize',12)
ylabel('Science Pixels','FontSize',12)
grid on
hold on
plot(momentumDumps - baseCadence, ones(length(momentumDumps),1),'m+')
plot(adcsNfpNotDump - baseCadence, 1.2*ones(length(adcsNfpNotDump),1),'r+')
s = axis;
axis([s(1) s(2) 0 3])
%legend('Science Pixels at FOV edge', 'momentumDump','ADCS notFinePoint notDump')
%%
tweakIndices = ismember(listOfCadences, tweakBad);
gapIndicators = (maxErr > threshold) | tweakIndices;
firstDifference = diff(gapIndicators);
startGaps = find( firstDifference > 0 ) + 1;
%first/last cadence bad case
if gapIndicators(1)
    startGaps = [1; startGaps];
end
lastInGaps = find( firstDifference < 0 );
if gapIndicators(end)
    lastInGaps = [lastInGaps; nCadences];
end
gapLength = lastInGaps - startGaps + 1;
bigGaps = find(gapLength >= persistence);
listOfExcludedCadences = [];
display('start list for xml file')
fid = fopen(['c' campaign '_coarse_point_DAF_SC_ADATTER_' num2str(threshold) '_' ...
    num2str(persistence) '.xml'],'w');
fprintf(fid,'ADATTERR-derived SC COARSE_POINT, program addatter_rms_to_cadence_sc.m\n');
fprintf(fid,['Engineer:  J. Van Cleve  Date: ' datestr(now)  '\n']);
fprintf(fid,'start Cadence %7.0f\n', listOfCadences(1));
fprintf(fid,'end Cadence %7.0f\n', listOfCadences(end));
fprintf(fid,'number of pre-tweak bad cadences %7.0f\n', length(tweakBad))
fprintf(fid,['Persistent Exceedances: threshold ' num2str(threshold) ' persistence ' num2str(persistence) '\n']);
fprintf(fid,'Relative Cadences \n');
for j = 1:length(bigGaps)
%max/min cover first/last cadence bad case    
   fprintf(fid,'%5.0f %5.0f\n', max(1, listOfCadences(startGaps(bigGaps(j))) - listOfCadences(1) + 1 - excludeBufferSize), ...
       min(nCadences, listOfCadences(lastInGaps(bigGaps(j))) - listOfCadences(1) + 1 + excludeBufferSize));
   plotGapStart = max(minCadence, listOfCadences(startGaps(bigGaps(j))) - excludeBufferSize);
   plotGapEnd = min(maxCadence, listOfCadences(lastInGaps(bigGaps(j))) + excludeBufferSize);
   plot([plotGapStart plotGapEnd] - baseCadence,threshold*ones(2,1),'g','LineWidth',10)
   listOfExcludedCadences = [listOfExcludedCadences plotGapStart:plotGapEnd];
   if (j == 1)
       legend('Science Pixels at FOV edge', 'momentumDump','ADCS notFinePoint notDump','SO COARSE POINT','Location','Best','AutoUpdate','off')
   end
end
fprintf(fid,'Excluded unique persistent cadences  %3.0f = %3.1f percent \n', length(unique(listOfExcludedCadences)), (100*length(unique(listOfExcludedCadences))/length(listOfCadences)));

fprintf(fid,'copy and paste text below into xml file\n');
for i = 1:length(bigGaps)
sprintf(' <dataAnomaly type="COARSE_POINT" cadenceType="SHORT" startCadence="%5.0f" endCadence="%5.0f"/>\n',...
     max(minCadence, listOfCadences(startGaps(bigGaps(i))) - excludeBufferSize), ...
     min(maxCadence, listOfCadences(lastInGaps(bigGaps(i))) + excludeBufferSize))
fprintf(fid, ' <dataAnomaly type="COARSE_POINT" cadenceType="SHORT" startCadence="%5.0f" endCadence="%5.0f"/>\n',...
     max(minCadence, listOfCadences(startGaps(bigGaps(i))) - excludeBufferSize), ...
     min(maxCadence, listOfCadences(lastInGaps(bigGaps(i))) + excludeBufferSize));
end

%NO cadences allowed above this threshold
highBadCadences = setdiff(listOfCadences(find(maxErr > singleThreshold)), listOfExcludedCadences);
persistentExcludedCadences = listOfExcludedCadences;
listOfExcludedCadences = union(listOfExcludedCadences, highBadCadences);
fprintf(fid, '\n');
fprintf(fid, ['high single bad cadences, threshold = ' num2str(singleThreshold) '\n']);
for j = 1:length(highBadCadences)
    fprintf(fid,  ' <dataAnomaly type="COARSE_POINT" cadenceType="SHORT" startCadence="%5.0f" endCadence="%5.0f"/>\n',...
        highBadCadences(j), highBadCadences(j));
    plot(highBadCadences(j) - baseCadence, singleThreshold,'g+')
end

fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid,['Total number of unique excluded cadences %3.0f = %3.1f percent \n'], length(unique(listOfExcludedCadences)), (100*length(unique(listOfExcludedCadences))/length(listOfCadences)));
saveas(gcf,['c' campaign '_addatterr2pix_cadence_nodump_' num2str(threshold) '_' ...
    num2str(persistence) '_sc.fig'])
saveas(gcf,['c' campaign '_addatterr2pix_cadence_nodump_' num2str(threshold) '_' ...
    num2str(persistence) '_sc.png'])

save(['c' campaign '_sc_excluded_cadences.mat'],'listOfCadences','listOfExcludedCadences')
fclose(fid);
end

