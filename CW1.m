% C O D E   T A K E S  A R O U N D  1 5  M I N S

%Load the EMNIST letters databse
load('emnist-letters.mat')

%extract features from train data
images=dataset.train.images;

%extract labels from train data
labels=dataset.train.labels;

%------------------------- K N N ------------------------------------------
%train KNN model 
tic
knnmodel= fitcknn(images,labels);
toc

%use trained model to predict
tic
predicted = predict(knnmodel, dataset.test.images(:,:));
toc

%put predicted label and actual label side by side
evaluate=[predicted dataset.test.labels(:,:)];

%compare predicted and actual
evaluate=[evaluate (evaluate(:,1)==evaluate(:,2))];

%find accuracy of knn
knn_accuracy=sum(evaluate(:,3))/length(evaluate(:,3))

%------------------------- Naive Bayes ------------------------------------
%train N.B model 
tic
NBmodel = fitcnb(double(images), labels, "DistributionNames","mn");
toc

%use trained model to predict
tic
predicted2 = predict(NBmodel, double(dataset.test.images(:,:)));
toc

%put predicted label and actual label side by side
evaluate2=[predicted2 dataset.test.labels(:,:)];

%compare predicted and actual
evaluate2=[evaluate2 (evaluate2(:,1)==evaluate2(:,2))];

%find accuracy of knn
NB_accuracy=sum(evaluate2(:,3))/length(evaluate2(:,3))

%------------------------- Decision Tree ------------------------------------
%train D.T model 
tic
DTmodel = fitctree(double(images), labels);
toc

%use trained model to predict
tic
predicted3 = predict(DTmodel, double(dataset.test.images(:,:)));
toc

%put predicted label and actual label side by side
evaluate3=[predicted3 dataset.test.labels(:,:)];

%compare predicted and actual
evaluate3=[evaluate3 (evaluate3(:,1)==evaluate3(:,2))];

%find accuracy of knn
DT_accuracy=sum(evaluate3(:,3))/length(evaluate3(:,3))


%--------------------ERROR ANALYSIS----------------------------------------

%Confusion Matrix for DT
DTodelCM = confusionchart(dataset.test.labels(:,:),predicted3)









