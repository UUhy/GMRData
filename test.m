%This is a test script for GMRData

dat = GMRData(); 
dat = dat.setFull(true);

Chip1 = dat.loadDataSet('data/C1');
figure('name','Chip 1');
Chip1.plotAll;