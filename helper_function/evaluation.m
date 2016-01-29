function [Acc, Sp, Se, MCC] = evaluation(vmask, rmask, cmask)
% Evaluation of vessel detection results
% Inputs:
% vmask: ground truth vessel mask
% rmask: result mask
% cmask: area of consideration
%
% Outputs:
% Acc: Accuracy
% Sp: Specificity
% Se: Sensitivity
% MCC: Mathews correlation coefficient
%% Obtain TP, TN, FP, FN
vmask  = (vmask >0); 
TP_mask = (vmask==1).*(rmask==1).*cmask;
TP = sum(TP_mask(:));
TN_mask = (vmask==0).*(rmask==0).*cmask;
TN = sum(TN_mask(:));
FP_mask = (vmask==0).*(rmask==1).*cmask;
FP = sum(FP_mask(:));
FN_mask = (vmask==1).*(rmask==0).*cmask;
FN = sum(FN_mask(:));

%% Display intermediate results for checking
% figure; subplot(2,2,1); imagesc(vmask); colormap gray; axis image; axis off; title('ground truth');
% subplot(2,2,2); imagesc(rmask); colormap gray; axis image; axis off; title('result mask');
% subplot(2,2,3); imagesc(vmask - rmask); colormap gray; axis image; axis off; title('difference');
% subplot(2,2,4); imagesc(FN_mask); colormap gray; axis image; axis off; title('FN');

%% calculation of metrics
N = TP+TN+FP+FN;
Acc = (TP + TN)/N;
Se = TP/(TP + FN);
Sp = TN/(TN + FP);

S = (TP + FN)/N;
P = (TP + FP)/N;
MCC = (TP/N - S*P)/sqrt(P*S*(1-S)*(1-P));

end

