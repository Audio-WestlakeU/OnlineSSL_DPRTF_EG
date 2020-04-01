function TDOA = compute_TDOA(micPos,canPos,MP)

micNum = size(micPos,1);
canNum = size(canPos,1);
mpNum = size(MP,2);

TDOA = zeros(mpNum,canNum);

dis = sqrt(sum((repmat(reshape(micPos,[micNum,1,3]),[1,canNum,1])-repmat(reshape(canPos,[1,canNum,3]),[micNum,1,1])).^2,3));

for mp = 1:mpNum
    TDOA(mp,:) = (dis(MP(2,mp),:)-dis(MP(1,mp),:))/344;
end

% figure;imagesc(log(abs(TDOA)))