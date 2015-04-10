% delivers errors on whole dataset reading the output files from a folder1
% and folder2 -- see run_pwrs_red parameter: storeFolder
function fullResult= evalSceneFlowKittiFolder ( )

compareGT = 0;
numList = 0:193;

folder2 = '../test1/';
folder1 = '../test2/';

doDifference = 2;

% disable certain images:
%  [a,b] = find(numList == 0);numList(b) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% at 3 pixel occ/nocev
allScores=[];

epes   = [0,0,0,0];
epes1  = [0,0,0,0];
epeOcc = [0,0];
epeNoc = [0,0];
epeOHit = [0,0];
epeNHit = [0,0];
allEpe1=[];

bestCombi = zeros(4);
listPic = [];
elms = 0;
sum = zeros(4);
elms2 = 0;
sum2 = zeros(4);
failImage = 0;
failImages = [];
for i=1:numel(numList)
  innerF = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  scoreFile  = dir(fullfile(folder1, innerF, sprintf('/RESULTS_K%03d_10*.txt', numList(i)) ));
  fName = fullfile(folder1, innerF, scoreFile.name);
  fid = fopen ( fName );
  if fid<0
    scoreFile  = dir(fullfile(folder1, innerF, sprintf('/RESULTS_%03d*.txt', numList(i)) ));
    fName = fullfile(folder1, innerF, scoreFile.name);
    fid = fopen ( fName );
  end 
  
  scoreFile2 = dir(fullfile(folder2, innerF, sprintf('/RESULTS_K%03d_10*.txt', numList(i)) ));
  fName2 = fullfile(folder2, innerF, scoreFile2.name);
  fid2 = fopen ( fName2 );
  if fid2<0
    scoreFile2  = dir(fullfile(folder2, innerF, sprintf('/RESULTS_%03d*.txt', numList(i)) ));
    fName2 = fullfile(folder2, innerF, scoreFile2.name);
    fid2 = fopen ( fName2 );
  end
  
  if fid > -1 && fid2 > -1
    try
      res = reshape( fscanf(fid, '%*s 2/3/4/5 %f & %f & %f & %f'), [4,4]);
    catch
      breakHere = 1;
      fclose( fid );fclose( fid2 );
      continue;
    end

    frewind(fid);
    for iiii=1:6 fgets(fid);end
    res4 = fscanf(fid, 'DispEPE %f & %f \nFlowEPE %f & %f');
    if ~isempty(res4)
      if any(isnan(res4))
        stophere = 1;
      else
        epes1 = epes1 + res4';
      end
    end
    try
      res2 = reshape( fscanf(fid2, '%*s 2/3/4/5 %f & %f & %f & %f'), [4,4]);
    catch
      breakHere = 1;
      fclose( fid );fclose( fid2 );
    end
    res = res(:,[2,4,1,3]);
    sum = sum + res;
    elms = elms +1;
    
    res2 = res2(:,[2,4,1,3]);
    sum2 = sum2 + res2;
    elms2 = elms2 +1;

    frewind(fid2);
    for iiii=1:6 fgets(fid2);end
    res3 = fscanf(fid2, 'DispEPE %f & %f \nFlowEPE %f & %f');

    if ~isempty(res3)
      if any(isnan(res3))
        stophere = 1;
      else
        epes = epes + res3';
      end
    end
    
    fclose( fid );fclose( fid2 );
    
    if res(2,2)< res2(2,2)
      bestCombi = bestCombi + res;
    else
      bestCombi = bestCombi + res2;
    end
    %%%%%%%%% differnece of methods :
    if doDifference == 1
      res = res2 - res;
      res4 = res3-res4;      
    end
    if doDifference == 2
      res = res - res2;
      res4 = res4-res3;      
    end    
    %%%%%%%%%%%%%%
    %%%%%%%%%%%%%
    listPic(end+1)     =  numList(i);
    allScores(1,end+1) = res(1,1);
    allScores(2,end)   = res(2,1);
    allScores(3,end)   = res(3,1);
    allScores(4,end)   = res(4,1);
    allEpe1(1:4,end+1) = res4';    
    pos=4;
    for ii=2:4
      for iii = 1:4
        pos = pos + 1;
        allScores(pos,end) = res(iii,ii);
      end
    end
    %%%%%%%%%%%%%
    
  else
    failImage = failImage+1;
    failImages(end+1) = numList(i);
    if fid > -1
      if compareGT
        res = reshape( fscanf(fid, '%*s 2/3/4/5 %f & %f & %f & %f'), [4,4]);
        res2 = res;
      end
      fclose( fid );
    end
    if fid2 > -1
      if compareGT
        res2 = reshape( fscanf(fid2, '%*s 2/3/4/5 %f & %f & %f & %f'), [4,4]);
        res = res2;
      end
      fclose( fid2 );
    end
    if compareGT
      sum2 = sum2 + res2;
      elms2 = elms2 +1;
      sum = sum + res;
      elms = elms +1;
    end
  end
end
fullResult1 = sum ./ elms;
fullResult1.*100

fullResult2 = sum2 ./ elms2;
fullResult2.*100
fprintf('Elements%03d\n', elms2);

fprintf('%.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f\n', 100*fullResult1(:));
fprintf('%.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f\n', 100*fullResult2(:));
% 
bestCombi = bestCombi ./ elms;
fprintf('%.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f\n', 100*bestCombi(:));

endNr = min(7, numel(numList));
[Worst5 Worst5Order] = sort(allEpe1(4,:));
Worst5Order = listPic(Worst5Order);
fprintf('Worst 5 flowEPE occ:\n %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f\n', cat( 1, Worst5Order(end:-1:end-endNr), Worst5(end:-1:end-endNr) ));
[Worst5 Worst5Order] = sort(allEpe1(3,:));
Worst5Order = listPic(Worst5Order);
fprintf('Worst 5 flowEPE noc:\n %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f\n', cat( 1, Worst5Order(end:-1:end-endNr), Worst5(end:-1:end-endNr) ));
[Worst5 Worst5Order] = sort(allEpe1(2,:));
Worst5Order = listPic(Worst5Order);
fprintf('Worst 5 dispEPE occ:\n %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f\n', cat( 1, Worst5Order(end:-1:end-endNr), Worst5(end:-1:end-endNr) ));
[Worst5 Worst5Order] = sort(allEpe1(1,:));
Worst5Order = listPic(Worst5Order);
fprintf('Worst 5 dispEPE noc:\n %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f\n', cat( 1, Worst5Order(end:-1:end-endNr), Worst5(end:-1:end-endNr) ));
%%%
endNr = min(7, numel(numList));
[Worst5 Worst5Order] = sort(allScores(2,:));
Worst5Order = listPic(Worst5Order);
fprintf('Worst 5 flow occ:\n %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f\n', cat( 1, Worst5Order(end:-1:end-endNr), 100*Worst5(end:-1:end-endNr) ));
[Worst5 Worst5Order] = sort(allScores(6,:));
Worst5Order = listPic(Worst5Order);
fprintf('Worst 5 flow noc:\n %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f\n', cat( 1, Worst5Order(end:-1:end-endNr), 100*Worst5(end:-1:end-endNr) ));
[Worst5 Worst5Order] = sort(allScores(10,:));
Worst5Order = listPic(Worst5Order);
fprintf('Worst 5 disp occ:\n %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f\n', cat( 1, Worst5Order(end:-1:end-endNr), 100*Worst5(end:-1:end-endNr) ));
[Worst5 Worst5Order] = sort(allScores(14,:));
Worst5Order = listPic(Worst5Order);
fprintf('Worst 5 disp noc:\n %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f\n', cat( 1, Worst5Order(end:-1:end-endNr), 100*Worst5(end:-1:end-endNr) ));

epes  = epes./ elms;
epes1 = epes1./ elms;
fprintf('EPED: %.3f, %.3f EPEF: %.3f %.3f\nEPED: %.3f, %.3f EPEF: %.3f %.3f\n', epes1, epes);

return
fullResult = cat( 3, fullResult, fullResult1(:,:,1) );