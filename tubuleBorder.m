function  [approxLines lineMap flag]= tubuleBorder(edg,m2o,LGM)
% function  [approxLines lineMap flag]= tubuleBorder(edg,m2o,LGM)
% function to model border of tubule as described in section 4.3 "Detect the
% tubules' border" which approximate edges responses in the border of the
% tubules with line and then connect them. 
% If the function called with no input argument, it will run the algorithm
% for all image in the default folder.
%  Input, 
%   edg, the binary image of edges, can be computed using PreProcess.
%   m2o, an image that all regions in edg are labeled with the second
%   dimesion of Podczeck value.
%   LGM, a binary image of lumen and germ-mass.
%  Output:
%   approxLines, a data structure that hold infomation about lines that are
%   fitted on the edges.
%   lineMap, an image with the same size as edg. the approximated lines are
%   represented by value one and connecting lines that are proposed by the
%   algorithm to fill the gap in the border of the tubules are represented
%   by two. This is called boundary Map in the report.
%   flag, boolean shows wheather the boundary map is acceptable or not.
% 
% See also PreProcess, measure, combineGMLumen
% 
% Author: A.Rahim Kadkhodamohammadi (r.k.mohammadi@gmail.com)
% 01 June 2012 CBA, Uppsala University
%--------------------------------------------------------------

if nargin < 3
    verbose =1;
    folder = '../images/GATA-4/';
    files = dir([ folder '*.TIF']);
    for ii = 5 : numel(files)% : -1 : 2
        en_im = imread([folder files(ii).name]);
        inFile =['../../tmpResult/mat/' files(ii).name(1:end-4) 'edg.mat'];
        load(inFile);
        outFile = ['../../tmpResult/' files(ii).name(1:end-4) 'lines.tif'];
        approxLines = findBoarder(edg,m2o, verbose, LGM, en_im, outFile, inFile);
    end
else
    [approxLines lineMap flag] = findBoarder(edg,m2o, 0, LGM);
end
end
    



function [approxLines lineMap2 flag] = findBoarder(edg,m2o, verbose,LGM, en_im, outFile, inFile)
    
    %internal parameter
    maxDist = 100; % maximum neighborhood distance
    angDist = 20; % maximum angular distance (ang+20>a>ang-20 is accepted)
    
    sz = size(edg);
    lines= im2mat(and(m2o> 0, m2o<0.08));
%     dipshow(2,lines)
    skLine = bwmorph(lines, 'skel', Inf); % all the lines be only one pixel thick
%     dipshow(skLine)
    %branchP = bwmorph(skLine, 'branchpoints');
    fk = ones(3,3);
    fk(2,2) = 8;
    filterLine = imfilter(double(skLine), fk, 'same');
    branchP = filterLine > 10;
    skLine(branchP) = 0;
%     dipshow(skLine);
    
    skLbl = label(skLine, 2, 10);
    podoz = measure(skLbl,  [], 'PodczeckShapes', [],2);
    m2skLbl = msr2obj(skLbl, podoz, 'PodczeckShapes', 2);
    finalLines = and(m2skLbl>0, m2skLbl<0.08);
%     dipshow(finalLines)
    
    % labe lnew lines and normalize them by PCA
    bwL = bwlabel(im2mat(finalLines));
    if verbose 
        h = figure; imshow(en_im); % en_im
    end
    stats = regionprops(bwL, 'PixelIdxList');
    approxLines = repmat(struct, [numel(stats),1]);
    errorThreshold = 50;
    for li =1 : numel(stats)
        tl = stats(li).PixelIdxList;
        % convert 1d inedx to 2d index
        xy = fix(tl/sz(1))+1;
        xy(:,2) = mod(tl,sz(1));
        
        % convert to zero mean for PCA
        meanxy = mean(xy);
        xy = xy - repmat(meanxy, size(xy,1),1);
        [coef, score, latent] = princomp(xy);
        pr = minmax(score(:,1)');
        p1 = coef(:,1)'*pr(1) + meanxy;
        p2 = coef(:,1)'*pr(2) + meanxy;
        angd = atand((p1(2)-p2(2))/(p1(1)-p2(1)));
        PSlope = coef(2,2)/coef(1,2);
        if latent(2) < errorThreshold && verbose
            line([p1(1) p2(1)], [p1(2) p2(2)], 'color','r', 'LineWidth',2);
        end
%         if latent(2) >= errorThreshold
%             line([p1(1) p2(1)], [p1(2) p2(2)], 'color','w', 'LineWidth',2);
%         end
        approxLines(li).p = [p1; p2];
%         approxLines(li).p2 = p2;
        approxLines(li).angd = [angd, angd];
        approxLines(li).prependicularSlope = [PSlope PSlope];
        approxLines(li).latent = latent;
        approxLines(li).coef = coef;
        approxLines(li).meanxy = meanxy;
        approxLines(li).score = score;
        
        % ignore lines that their length is more or less qeual with the 
        % smallest latent
%         if 0.8 * sqrt(sum((p1-p2).^2)) < latent(2)
%              approxLines(li).ignore = true;
%         else
            approxLines(li).ignore = false;
%         end
    end
    approxLines = splitErrorProneLines(approxLines, errorThreshold ,verbose);
    
%     maxDist = 100;
%     angDist = 20;
    [neighborLines numConnectingLines]= ...
        getHeadandTailN(approxLines, maxDist, angDist);
    lineMap = drawLines(approxLines, neighborLines,sz, zeros(numel(approxLines),1));
    numConnections= numConnectedLines(neighborLines);
    filteredLines = eliminateLines(numConnections,approxLines, LGM);
    lineMap2 = drawLines(approxLines, neighborLines,sz, filteredLines);
    %dipshow(cat(3,lineMap,lineMap2),'labels');
    SE = strel('disk',7,8);
    
%    check wheather the boundary map is accpetable or not by :
%     look at the total number of edges and connection between edges.   
    if numel(approxLines) < 200 || (.5 *numel(approxLines))>numConnectingLines
        flag = 0;
    else
        flag =1;
    end
    lineMapClosed= imclose(lineMap2>0,SE);
    if verbose
        if ~flag 
            imwrite(lineMapClosed, [outFile(1:end-4) 'reject.tif']);
        else
            imwrite(lineMapClosed, outFile);
        end
        save(inFile, 'lineMap', 'lineMap2', 'approxLines','neighborLines','-append');
    end
%     saveas(h, outFile);
    close all
end


function lineMap = drawLines(lines, connectingLines,sizeIm, filtLine)
lineMap = zeros(sizeIm(1:2));
for li = 1 :  numel(lines)
    if ~lines(li).ignore && ~filtLine(li)
        % draw the line
        p1 = lines(li).p(1,:);
        p2 = lines(li).p(2,:);
        if lines(li).angd(1) == lines(li).angd(2)
            [x, y, ~] = improfile(lineMap, [p1(1), p2(1)], [p1(2), p2(2)]);
            x = max(round(x),1);
            y = max(round(y),1);
            for p = 1 : length(x)
                %             try
                lineMap(y(p),x(p)) = 1;
                %             catch err
                %                 disp(['stopped', err]);
                %             end
            end
        else % the is split into two parts
            p11 = lines(li).slines(1,:);
            p12 = lines(li).slines(2,:);
            p21 = lines(li).slines(3,:);
            p22 = lines(li).slines(4,:);
            [x, y, ~] = improfile(lineMap, [p11(1), p12(1)], [p11(2), p12(2)]);
            x = max(round(x),1);
            y = max(round(y),1);
            for p = 1 : length(x)
               lineMap(y(p),x(p)) = 1;
            end
            [x, y, ~] = improfile(lineMap, [p21(1), p22(1)], [p21(2), p22(2)]);
            x = max(round(x),1);
            y = max(round(y),1);
            for p = 1 : length(x)
               lineMap(y(p),x(p)) = 1;
            end
        end
        % draw conecting lines
        connecting2P1= connectingLines(li).p1;
        for ii = 1: numel(connecting2P1)
            if filtLine(connecting2P1(ii).lineInd) == 0
                tp1 = connecting2P1(ii).nearstP;
                [x, y, ~] = improfile(lineMap, [p1(1), tp1(1)], [p1(2), tp1(2)]);
                x = max(round(x),1);
                y = max(round(y),1);
                for p = 1 : length(x)
                    lineMap(y(p),x(p)) = 2;
                end
            end
        end
        
        connecting2P2= connectingLines(li).p2;
        for ii = 1: numel(connecting2P2)
            if filtLine(connecting2P2(ii).lineInd) == 0
                tp2 = connecting2P2(ii).nearstP;
                [x, y, ~] = improfile(lineMap, [p2(1), tp2(1)], [p2(2), tp2(2)]);
                x = max(round(x),1);
                y = max(round(y),1);
                for p = 1 : length(x)
                    lineMap(y(p),x(p)) = 2;
                end
            end
        end
    end
end
lineMap = lineMap(1:sizeIm(1), 1:sizeIm(2));
end

function lines = splitErrorProneLines(lines,errorThreshold,verbose)

% find error prone lines
latent = cat(2,lines.latent);
ind = find(latent(2,:) > errorThreshold);

% ignore lines that their approximated line shorter than 1.5 of its pixels
% it means that the line somewhat similar to circle
for i =1 : numel(ind)
    p = lines(ind(i)).p;
    len = sqrt((p(1,1)-p(2,1))^2+ (p(1,2)-p(2,2))^2);
    if 1.5*len < numel(lines(ind(i)).score(:,1));
        lines(ind(i)).ignore = true;
    end
end
% split lines into two line to reduce error
for i = 1 : numel(ind)
    ci = ind(i);
    [~ ,sId] = sort(lines(ci).score(:,1));
    score = lines(ci).score(sId,:);
    ls = size(score,1);
%     xy1 = (lines(ci).coef'*score(1:fix(ls/2),:)')' + ...
%         repmat(lines(ci).meanxy, fix(ls/2),1);
%     xy2 = (lines(ci).coef'*score(fix(ls/2):end, :)' + ...
%         repmat(lines(ci).meanxy, ls-fix(ls/2)+1,1)')';
    xy1 = (score(1:fix(ls/2),:)*lines(ci).coef') + ...
        repmat(lines(ci).meanxy, fix(ls/2),1);
    xy2 = (score(fix(ls/2):end,:)*lines(ci).coef') + ...
        repmat(lines(ci).meanxy, ls-fix(ls/2)+1,1);
    meanxy1 = mean(xy1);
    xy1 = xy1 - repmat(meanxy1, size(xy1,1),1);
    [coef, score, latent1] = princomp(xy1);
    pr = minmax(score(:,1)');
    p11 = coef(:,1)'*pr(1) + meanxy1;
    p12 = coef(:,1)'*pr(2) + meanxy1;
    angd1 = atand((p11(2)-p12(2))/(p11(1)-p12(1)));
    PS1 = coef(2,2)/coef(1,2);
    
    meanxy2 = mean(xy2);
    xy2 = xy2 - repmat(meanxy2, size(xy2,1),1);
    [coef, score, latent2] = princomp(xy2);
    pr = minmax(score(:,1)');
    p21 = coef(:,1)'*pr(1) + meanxy2;
    p22 = coef(:,1)'*pr(2) + meanxy2;
    angd2 = atand((p21(2)-p22(2))/(p21(1)-p22(1)));
    PS2 = coef(2,2)/coef(1,2);
    
    if verbose
       line([p11(1) p12(1)], [p11(2) p12(2)], 'color','g', 'LineWidth',2);
       line([p21(1) p22(1)], [p21(2) p22(2)], 'color','g', 'LineWidth',2);
    end
    
    %store only the beging of first line and end of second line
    lines(ci).p = [p11; p22];
%     lines(ci).p2 = p22;
    lines(ci).angd = [angd1, angd2];
    lines(ci).prependicularSlope = [PS1, PS2];
    lines(ci).latent = [latent1,latent2];
    lines(ci).coef = coef;
    lines(ci).slines = [p11; p12; p21;p22];
end

end

function [neighborLines count] = getHeadandTailN(lines, maxDist, angDist)
numLine = length(lines);
neighborLines = struct;
count = 0;
for i = 1 : numLine
    p1 = lines(i).p(1,:);
    p2 = lines(i).p(2,:);
    m = lines(i).prependicularSlope(1);
    angd = lines(i).angd(1);
    neighborLines(i).p1 = getNeighboringLines(lines,i,numLine, p1,p2,m,angd,maxDist, angDist);
    count = count + numel(neighborLines(i).p1);
    m = lines(i).prependicularSlope(2);
    angd = lines(i).angd(2);
    neighborLines(i).p2 = getNeighboringLines(lines,i,numLine, p2,p1,m,angd,maxDist, angDist);
    count = count + numel(neighborLines(i).p2);
end

end


function TN = getNeighboringLines(lines,cLine,numLine,p1,p2,m,angd, maxDist, angDist)
TN = [];
lAng = angd-angDist;
upAng  = angd + angDist;
if lAng>=-90 && upAng<=90
    angCondition = @(x)(x>lAng && x< upAng);
elseif lAng < -90
    angCondition = @(x)(x>(180 +lAng) || x< upAng);
else
    angCondition = @(x)(x>lAng || x< (upAng-180));
end
b = p1(2) - m*p1(1);
acceptSign = -1 * sign(m*p2(1) + b -p2(2));
for j = cLine+1 : numLine
    if ~lines(j).ignore
        p = lines(j).p;
        d1 = sqrt((p(1,1)-p1(1))^2 + (p(1,2)-p1(2))^2);
        d2 = sqrt((p(2,1)-p1(1))^2 + (p(2,2)-p1(2))^2);
        minInd = 1;
        if d1>d2; d1=d2; minInd = 2; end
        
        % check whether line is in acceptable distnace and complitely fell
        % in one side of the current line
        line2Ang = lines(j).angd(minInd);
        l2lAng = atand((p1(2)-p(minInd,2))/(p1(1)-p(minInd,1)));
        if (d1<maxDist) && angCondition(line2Ang) && angCondition(l2lAng) ...
                && sign(m*p(1,1) + b -p(1,2))==acceptSign ...
                && sign(m*p(2,1) + b -p(2,2))==acceptSign
            t.lineInd= j;
            t.pointInd = minInd;
            t.dist = d1;
            t.nearstP = p(minInd, :);
            TN= [TN t];
            line([p1(1) p(minInd,1)], [p1(2) p(minInd,2)], 'color','b', 'LineWidth',2);
        end
    end
end
end

function numConnections= numConnectedLines(neighborLines)
num = numel(neighborLines);
numConnections = zeros(num,2);

for i =1 : num
    tlen = numel(neighborLines(i).p1);
    numConnections(i,1) = numConnections(i,1) + tlen;
    for ii =1 : tlen
        pointInd = neighborLines(i).p1(ii).pointInd;
        lineInd = neighborLines(i).p1(ii).lineInd;
        numConnections(lineInd,pointInd) =numConnections(lineInd,pointInd)+1; 
    end
    tlen = numel(neighborLines(i).p2);
    numConnections(i,2) = numConnections(i,2) + tlen;
    for ii =1 : tlen
        pointInd = neighborLines(i).p2(ii).pointInd;
        lineInd = neighborLines(i).p2(ii).lineInd;
        numConnections(lineInd,pointInd) =numConnections(lineInd,pointInd)+1; 
    end
end
end

function filteredLines = eliminateLines(numConns,Lines, LGM)
thre = 3;

% check Lines that are disconnected from one side and they have only few(2)
% connection on the other side

tmp1 = numConns(:,1)*1000+numConns(:,2);
tmp1 = and(tmp1>0,tmp1<thre);
tmp2 = numConns(:,1)+numConns(:,2)*1000;
tmp2 = and(tmp2>0, tmp2<thre);

tmp = or(tmp1, tmp2);
Ind1 = find(tmp);

filteredLines = checkwithGM(Ind1, Lines, LGM);
% check lines that have zero connection.
Ind2 = find(numConns(:,1)+numConns(:,2) == 0);
for i = 1 : numel(Ind2)
    if Lines(Ind2(i)).angd(1) ==  Lines(Ind2(i)).angd(1) % it is not divied line
        pp = Lines(Ind2(i)).p;
        lenL = sqrt(sum((pp(1,:)-pp(2,:)).^2));
        if lenL < 40
            filteredLines(Ind2(i)) = 1;
            Ind2(i) = 0;
        end
    end
end
remainInd = find(Ind2);
rInd2 = Ind2(remainInd);
filteredLines = checkwithGM(rInd2, Lines, LGM, filteredLines);

end

function filteredLines = checkwithGM(Ind, Lines, LGM, filteredLines)
MaxGMDist = 160; % 160/2
LGM=bwareaopen(LGM, 20000); %remove small reions which should be done before
boundary = bwboundaries(LGM, 'noholes');
numGM = numel(boundary);
if nargin < 4
    filteredLines = zeros(numel(Lines),1);
end
for i = 1 : numel(Ind)
    ci = Ind(i);
    if Lines(ci).angd(1) ==  Lines(ci).angd(1) % it is not divied line
        p1 = Lines(ci).p(1,:);
        p2 = Lines(ci).p(2,:);
        d1 = zeros(numGM,1);
        d2 = zeros(numGM,1);
        mp1 = zeros(numGM,1);
        mp2 = zeros(numGM,1);
        for ii =1 : numGM
           tb = boundary{ii};
           td = sqrt((tb(:,2)-p1(1)).^2 + (tb(:,1)-p1(2)).^2);
           [d1(ii) mp1(ii)]= min(td);
           [d2(ii) mp2(ii)]= min(sqrt((tb(:,2)-p2(1)).^2 + (tb(:,1)-p2(2)).^2));
        end
        if (min(d1) + min(d2)) > MaxGMDist
            filteredLines(ci) = 1;
        else
           [am ai] = min(d1);
           [bm bi] = min(d2);
           if am<bm
               gmId = ai;
           else
               gmId = bi;
           end
           b1 = boundary{gmId}(mp1(gmId),:);
           b2 = boundary{gmId}(mp2(gmId),:);
           db = sqrt((b1-b2).^2);
           dp = sqrt((p1-p2).^2);
           
           if db < (0.8*dp)
               filteredLines(ci) = 1;
           end
               
        end
            
    end
end
end