function analyzeResults

% Initialize parameters
resultFolderName = 'userOpinionScores\';

% Analyze the results for the first session
resultStr_S1 = analyzeSession(resultFolderName, 'session01');

% Analyze the results for the second session
resultStr_S2 = analyzeSession(resultFolderName, 'session02');

% Merge two sessions using the common set
[resultStr, indsOrig] = mergeSessions(resultStr_S1, resultStr_S2);

% Plot results
plotResults(resultStr, indsOrig);

% Save results
% save('vsenseVVDB2_subjectiveScores.mat', 'resultStr', 'indsOrig');

end

%--------------------------------
% Analyze session with sessionLabel
%--------------------------------
function resultStr = analyzeSession(folderName, sessionLabel)    
    % --- Read the CSV and find the votes
    resultStr = findVotes(folderName, sessionLabel);
    resultStr = treatVotes(resultStr);
end

%--------------------------------
% Find votes per session
%--------------------------------
function resultStr = findVotes(folderName, sessionLabel)
    % Find first session CSV files
    dirS = dir([folderName  '*_' sessionLabel '_*.csv']);
    
    % --- Read the CSV and process
    % For each CSV folder
    for ind = 1:length(dirS)
        % Read the table
        Tab = readtable([pwd filesep folderName dirS(ind).name]);

        % Find sorted names and sortInd
        stNames = Tab.StimuliName(7:end);
        [stNamesS, sortInd] = sort(stNames);

        % Get the votes 
        votesSubj = Tab.Vote(7:end);
        votesSubjSorted_Raw = votesSubj(sortInd);

        % Store the values
        sVotes(:,ind) = votesSubjSorted_Raw;
        sNames(:,ind) = stNamesS;
    end
    
    % Return the results
    resultStr.sVotes = sVotes;
    resultStr.sNames = sNames;
end

%--------------------------------
% Treat votes with resultStr
%--------------------------------
function resultStr = treatVotes(resultStr)
    % --- Process and find z scores for all the stimuli
    s_ij = resultStr.sVotes;
    [MOS, STD, CI] = treatVotesForSij(s_ij);
    resultStr.MOS = MOS;
    resultStr.STD = STD;
    resultStr.CI  = CI;
end

%--------------------------------
% Treat votes for give s_ij
%--------------------------------
function [MOS, STD, CI] = treatVotesForSij(s_ij)
    % Compute MOS and CI
    [MOS_j, STD_j, CI_j] = findMosCi(s_ij, size(s_ij,2));
    mu_rMos  = mean(MOS_j);
    std_rMos = std(MOS_j);

    % -- try without subtracting the references value
    mu_i  = mean(s_ij,1);
    std_i = std(s_ij, [], 1);
    z_ij  = (s_ij - repmat(mu_i,size(s_ij,1),1))./...
             repmat(std_i,size(s_ij,1),1);
    
    % Remove outliers
    [zOR_ij, ~] = outlierRemoval(z_ij);
    userCount = size(zOR_ij,2);
    [MOS_zj, STD_zj, CI_zj] = findMosCi(zOR_ij, userCount);
    mu_zMos  = mean(MOS_zj);
    std_zMos = std(MOS_zj);

    % Compute final MOS?=
    zOR_ij_scaled = std_rMos*( (zOR_ij - mu_zMos)/std_zMos ) + mu_rMos;
    [MOS, STD, CI] = findMosCi(zOR_ij_scaled, size(zOR_ij_scaled,2));
    display("MOS")
    display(MOS)
    display("std:")
    display(STD)
    display("end std")
end

%--------------------------------
% Merge two sessions using the common set
%--------------------------------
function [resultStr, indsOrig] = mergeSessions(resStr_S1, resStr_S2)

    % Find indices first for the common set
    indStr = findIndices(resStr_S1, resStr_S2);

    % find aVal and bVal from y = ax + b
    session1votes_subset = resStr_S1.MOS(indStr.srcSubsSs1);
    y = session1votes_subset;
    session2votes_subset = resStr_S2.MOS(indStr.srcSubsSs2);
    x = session2votes_subset;
    aVal = corr(x,y)*std(y)/std(x);
    bVal = mean(y) - aVal*mean(x);
    
    % --- Merge all scores on session1 scale
    % Map the session2 votes over session1 scale!
    mappedMOS_S2 = aVal*resStr_S2.MOS + bVal;
    mappedSTD_S2 = aVal*resStr_S2.STD;
    mappedCI_S2  = aVal*resStr_S2.CI;
    
    
    % Init mappedMOS and mappedCI
    mappedMOS = NaN(length(indStr.sortedNamesAll),1);
    mappedCI  = NaN(length(indStr.sortedNamesAll),1);
    mappedSTD = NaN(length(indStr.sortedNamesAll),1);%zelf toegevoegd

    % Take values directly from session1
    mappedMOS(indStr.tarDifsSs1) = resStr_S1.MOS(indStr.srcDifsSs1);
    mappedCI(indStr.tarDifsSs1)  = resStr_S1.CI(indStr.srcDifsSs1);
    mappedSTD(indStr.tarDifsSs1) = resStr_S1.STD(indStr.srcDifsSs1);

    % Take mapped values from session 2 directly as well
    mappedMOS(indStr.tarDifsSs2) = mappedMOS_S2(indStr.srcDifsSs2);
    mappedCI(indStr.tarDifsSs2)  = mappedCI_S2(indStr.srcDifsSs2);
    mappedSTD(indStr.tarDifsSs2) = mappedSTD_S2(indStr.srcDifsSs2);

    % Take average for the common subset
    mappedMOS(indStr.tarSubsAll) = (resStr_S1.MOS(indStr.srcSubsSs1) +...
                                     mappedMOS_S2(indStr.srcSubsSs2))./2;
    mappedCI(indStr.tarSubsAll)  = (resStr_S1.CI(indStr.srcSubsSs1) +...
                                     mappedCI_S2(indStr.srcSubsSs2))./2;
    mappedSTD(indStr.tarSubsAll) = (resStr_S1.STD(indStr.srcSubsSs1)+...
                                     mappedSTD_S2(indStr.srcSubsSs2))./2;

    display(mappedMOS)
    display(mappedSTD)
    filename = "mos_std.csv"
    csvwrite(filename, [mappedMOS, mappedSTD])
	% Return
    resultStr.wRef_names = indStr.sortedNamesAll;
    resultStr.wRef_MOS   = mappedMOS;
    resultStr.wRef_CI    = mappedCI;
    indsOrig = indStr.indsOrig;
    
    % Compute values for the distortion cases as well
    distortedNames = indStr.sortedNamesAll;
    distortedNames(indStr.indsOrig) = [];
    mappedMOS_Dist = mappedMOS;
    mappedMOS_Dist(indStr.indsOrig) = [];
    mappedCI_Dist  = mappedCI;
    mappedCI_Dist(indStr.indsOrig) = [];
    
    % Find content, representation, compression names as a list
    namesStr =  uniqueNames(distortedNames);
    
    % Return only distortion as well
    resultStr.Dist_names     = distortedNames;
    resultStr.Dist_compTypes = namesStr.compTypNames;
    resultStr.Dist_compTypeQ = namesStr.compTypNamesQ;
    resultStr.Dist_contentQ  = namesStr.contentNamesQ;
    resultStr.Dist_MOS       = mappedMOS_Dist;
    resultStr.Dist_CI        = mappedCI_Dist;
    
end

%--------------------------------
% Find indices for both sides of the equation for merging operation
%--------------------------------
function indStr = findIndices(resStr_S1, resStr_S2)
    % ====== Merge all MOS
    sortedNamesAll = unique(cat(1,resStr_S1.sNames,...
                                  resStr_S2.sNames));
    sortedNamesInt = intersect(resStr_S1.sNames,...
                               resStr_S2.sNames);
    
    % Find source indices
    [~, srcSubsSs1, ~] = intersect(resStr_S1.sNames,...
                                   sortedNamesInt);
    [~, srcSubsSs2, ~] = intersect(resStr_S2.sNames,...
                                   sortedNamesInt);
    srcDifsSs1 = setdiff(1:length(resStr_S1.sNames), srcSubsSs1);
    srcDifsSs2 = setdiff(1:length(resStr_S2.sNames), srcSubsSs2);
    
    % Find target indices
    [~, tarSubsAll] = intersect(sortedNamesAll, sortedNamesInt);
    [~, tarDifsSs1] = intersect(sortedNamesAll,...
                                resStr_S1.sNames(srcDifsSs1));
    [~, tarDifsSs2] = intersect(sortedNamesAll,...
                                resStr_S2.sNames(srcDifsSs2));
                            
    % Find the indices for the uncompressed (i.e., hidden references)
    indsOrig = [];
    for ind = 1:length(sortedNamesAll)
        if strfind(sortedNamesAll{ind}, '_orig')
            indsOrig = cat(1, indsOrig, ind);
        end
    end
                            
	% Return
    indStr.sortedNamesAll = sortedNamesAll;
    indStr.srcSubsSs1 = srcSubsSs1;
    indStr.srcSubsSs2 = srcSubsSs2;
    indStr.srcDifsSs1 = srcDifsSs1;
    indStr.srcDifsSs2 = srcDifsSs2;
    indStr.tarSubsAll = tarSubsAll;
    indStr.tarDifsSs1 = tarDifsSs1;
    indStr.tarDifsSs2 = tarDifsSs2;
    indStr.indsOrig   = indsOrig;
end

%--------------------------------
% Find unique names
%--------------------------------
function namesStr =  uniqueNames(allNames)
    % Find content, representation, compression names as a list
    contentNames = {};
    represnNames = {};
    compTypNames = {};
    for ind = 1:length(allNames)
        splStr = strsplit(allNames{ind},{'_','.'});

        contentNames = cat(1, contentNames, splStr{1});
        represnNames = cat(1, represnNames, splStr{2});
        compTypNames = cat(1, compTypNames, splStr{3});
    end

    % Find unique content, representation, compression names
    contentNamesQ = unique(contentNames);
    represnNamesQ = unique(represnNames);
    compTypNamesQ = unique(compTypNames);
    
    % Return
    namesStr.contentNames  = contentNames;
    namesStr.represnNames  = represnNames;
    namesStr.compTypNames  = compTypNames;
    namesStr.contentNamesQ = contentNamesQ;
    namesStr.represnNamesQ = represnNamesQ;
    namesStr.compTypNamesQ = compTypNamesQ;
end

%--------------------------------
% Find MOS, standard deviation, and confidence interval
%--------------------------------
function [MOS, Std, CI] = findMosCi(data, userCount)
    % Check data dimensions
    if size(data,2) ~= userCount
        data = data';
        assert(size(data,2) == userCount);
    end
    
    % Find MOS
    MOS = mean(data,2);
    
    % Find Std
    Std = std(data,0,2);
    
    % Find Confidence Interval
    % CI as described in ITU-T P.1401 Appendix III
    % ---
    % confLevel = 0.95;
    % critProb = 1 - ( 1 - confLevel )/2;
    subjNo = size(data, 2);
    degOfFrd = subjNo - 1;
    % - These tScores are ONLY for 0.95 confLevel ---!!!!
    tScores = [12.71; 4.303; 3.182; 2.776; 2.571; 2.447; 2.365; 2.306; 
               2.262; 2.228; 2.201; 2.179; 2.160; 2.145; 2.131; 2.120; 
               2.110; 2.101; 2.093; 2.086; 2.080; 2.074; 2.069; 2.064; 
               2.060; 2.056; 2.052; 2.048; 2.045; 2.042; 2.021; 2.009; 
               2.000; 1.990; 1.984; 1.980; 1.960];
    if degOfFrd > length(tScores)
        tScore = 1.96;
    else
        tScore = tScores(degOfFrd);
        disp("tscore:")
        disp(tScore)
    end
 
    % Compute confidence interval
    standError = Std/sqrt(subjNo);
    CI = standError*tScore;
    
end

%--------------------------------
% Outlier removal
%--------------------------------
function [resArrOut, outliers] = outlierRemoval(resArr)
    % Check ITU-R BT.500-13 Section 2.3.1 of Annex 2

    % Find Mean and Std before the procedure
    [resMean, resStd] = findMosCi(resArr, size(resArr,2));

    % Find Beta
    betaVal = moment(double(resArr),4,2)./...
             (moment(double(resArr),2,2).^2); 

    % Find indices where 2<B<4 and otherwise
    betaBtw24 = find( (betaVal < 4) & (betaVal > 2) );
    betaOthw = find( (betaVal >= 4) | (betaVal <= 2) );

    % If 2<B<4
    P_jkr1 = bsxfun(@lt, resMean(betaBtw24) + 2*resStd(betaBtw24),...
                         resArr(betaBtw24,:));
    Q_jkr1 = bsxfun(@gt, resMean(betaBtw24) - 2*resStd(betaBtw24),...
                         resArr(betaBtw24,:));

    % Otherwise
    P_jkr2 = bsxfun(@lt, resMean(betaOthw) + sqrt(20)*resStd(betaOthw),...
                         resArr(betaOthw,:));
    Q_jkr2 = bsxfun(@gt, resMean(betaOthw) - sqrt(20)*resStd(betaOthw),...
                         resArr(betaOthw,:));

    % Append two cases
    P_jkr = cat(1, P_jkr1, P_jkr2);
    Q_jkr = cat(1, Q_jkr1, Q_jkr2);

    % Find the sum
    P_i = sum(P_jkr,1);
    Q_i = sum(Q_jkr,1);

    % Find ratios
    firstRatio = (P_i + Q_i)/size(resArr,1);
    secondRatio = abs( (P_i - Q_i)./(P_i + Q_i) );

    % Find outliers
    outliers = find( (firstRatio > 0.05) & (secondRatio < 0.3) );
    if isempty(outliers)
        numOut = 0; 
    else
        numOut = length(outliers); 
    end

    % Prepare output
    resArrOut = resArr;
    resArrOut(:, outliers) = [];

    disp([' -- ' num2str(numOut) ' outliers detected and removed! --']);
    if ~isempty(outliers)
        disp(outliers)
    end

end

%--------------------------------
% Find Matlab version
%--------------------------------
function outFlag = newerThan2016b
    verString = version;
    verNum = regexp(verString, 'R([0-9])+([a-b])', 'match');
    if and(str2double(verNum(2:end-1)) >= 2016, strcmp(verNum(end),'b'))
        % Newer than 2016b
        outFlag = 1;
    else
        outFlag = 0;
    end
end

%--------------------------------
% Plot results
%--------------------------------
function plotResults(resultStr, indsOrig)

    % ---- Plot all together
    figure, errorbar(resultStr.wRef_MOS, resultStr.wRef_CI)
    if newerThan2016b
        xticks(1:size(resultStr.wRef_MOS)), xtickangle(90); 
        xticklabels(resultStr.wRef_names);
    else
        set(gca, 'xtick', 1:size(resultStr.wRef_MOS),...
                 'xticklabelrotation', 90,...
                 'xticklabel', resultStr.wRef_names);
    end
    title ('MOS for and CI for the all votes');
    
    % Highlight references
    vecY = NaN(size(resultStr.wRef_MOS)); 
    vecY(indsOrig) = resultStr.wRef_MOS(indsOrig);
    hold on; plot(vecY, 'ro')
    axis([1 length(resultStr.wRef_names) 0 90])
    
    % ---- Plot individual graphs
    includeOBJ = true;
    plotIndividualResults(resultStr, includeOBJ);
    
    includeOBJ = false;
    plotIndividualResults(resultStr, includeOBJ);
end

%--------------------------------
% Plot for each content
%--------------------------------
function plotIndividualResults(resultStr, includeOBJ)

% Load bitrates
load('vsenseVVDB2_bitrates.mat', 'namesVecDist',...
                                 'fileSizeVecDist',...
                                 'rateKbpsVecDist',...
                                 'rateMbpsVecDist');

% Prepare the lines, markers, and legend entries
colsL = {'ro-','bs:','md--','kv-.', 'rs--','bs--','ms--','ks--'};
for indL = 1:length(resultStr.Dist_compTypeQ)
    if ~isempty(strfind(resultStr.Dist_compTypeQ{indL},'Draco'))
        legendList{indL} = 'Draco+JPEG';
    elseif ~isempty(strfind(resultStr.Dist_compTypeQ{indL},'TMC1'))
        legendList{indL} = 'G-PCC (RAHT)';
    elseif ~isempty(strfind(resultStr.Dist_compTypeQ{indL},'TMC2-allIn'))
        legendList{indL} = 'V-PCC (AI)';
    elseif ~isempty(strfind(resultStr.Dist_compTypeQ{indL},'TMC2-rand'))
        legendList{indL} = 'V-PCC (RA)';
    else
        error('Something''s wrong');
    end
end

% Flag to either include or exclude OBJ results
if includeOBJ == true
    incOBJ = 1; % 1 to include OBJ, 2 to exclude
    includeOBJstr = 'compare all';
else
    incOBJ = 2; % 1 to include OBJ, 2 to exclude
    includeOBJstr = 'only point clouds';
end

% For each content
for ind = 1:length(resultStr.Dist_contentQ)
    % Prepare figure with name
    eval(['h_' num2str(ind) ' = figure(''name'',''' ...
         resultStr.Dist_contentQ{ind} ' - MOS vs Bitrate - ' ...
         includeOBJstr ''');']);

    % find the content indice
    indCont = find(~cellfun('isempty', strfind(resultStr.Dist_names, ...
                                          resultStr.Dist_contentQ{ind})));

    % Get different types (-orig)
    for indT = incOBJ:length(resultStr.Dist_compTypeQ)
        indType = find(~cellfun('isempty', ...
                                  strfind(resultStr.Dist_compTypes,...
                                  resultStr.Dist_compTypeQ {indT})));
        indIntr = intersect(indCont, indType);

        % Call relevant figure
        eval(['figure(h_' num2str(ind) '),']);
        errorbar(rateMbpsVecDist(indIntr), resultStr.Dist_MOS(indIntr), ...
                      resultStr.Dist_CI(indIntr), colsL{indT}); hold on;
        grid on; xlabel('Bitrate (Mbps)'), ylabel('MOS');
        title(resultStr.Dist_contentQ{ind})
        if ind == 1, 
            legend(legendList(incOBJ:end), 'location','southeast'); 
        end
    end
end

end