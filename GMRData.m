%GMRData Class
%
%The GMRData class is designed to facilitate analysis of GMR data taken
%from our GMR tester.  
%
%The data file format has the following criteria:
%   1)  ASCII file
%   2)  Tab delimited floating point numbers
%   3)  Each pair of columns contain f and resistance data
%   4)  The first pair of columns is the average f and resistance data
%
%The functions of this class are:
%   loadData        =   Loads a single gmr data file
%   loadDataSet     =   Loads a set of gmr data file with same base name
%                       and 1 unique numeric identifier
%   loadDataSetRep  =   Loads a set of gmr data file with same base name
%                       and 2 unique numeric identifier
%   setFull         =   Specify full loop or virgin loop
%   plotAll         =   Plot all data
%   plotOne         =   Plot one data
%   plotOne3D       =   Plot the data in 3D view
%   plotCompare     =   Plot all sensors on same figure
%   plotDifference  =   Plot difference of data
%   plotRBase       =   Plot base resistance with time
%
%Long Chang, UH, 3/2/2016

classdef GMRData
    
    properties(GetAccess=public, SetAccess=protected, Hidden=false)
        
        %Type of data
        dataType;
        
        %Virgin loop or Full loop
        full;
        
        %Number of datasets
        nData;
        
        %Field data
        %Data:          [N R]
        %DataSet:       [N S R]
        %DataSetRep:    [N S R]
        %   N = data points
        %   R = repeat
        %   S = sensor number
        f;
        
        %Resistance data
        %Data:          [N R]
        %DataSet:       [N S R]
        %DataSetRep:    [N S R]
        %   N = data points
        %   R = repeat
        %   S = sensor number
        r;
        
        %Base Resistance
        %Data:          [1 R]
        %DataSet:       [S R]
        %DataSetRep:    [S R]
        %   N = data points
        %   R = repeat
        %   S = sensor number
        rbase;
        
        %Field data average
        %Data:          [N 1]
        %DataSet:       [N S]
        %DataSetRep:    [N S]
        %   N = data points
        %   S = sensor number
        fA;
        
        %Resistance data average
        %Data:          [N 1]
        %DataSet:       [N S]
        %DataSetRep:    [N S]
        %   N = data points
        %   S = sensor number
        rA;
        
        %Base resistance average
        %Data:          [1 1]
        %DataSet:       [1 S]
        %DataSetRep:    [1 S]
        %   S = sensor number
        rbaseA;
        
        %Color scheme
        c;
    end
    
    properties(GetAccess=public, SetAccess=protected, Hidden=true)
        %No Private Variables
    end
    
    properties(GetAccess=private, Constant=true, Hidden=true)
        %No Constants
    end
    
    methods (Access=public, Sealed=true)
        
        function obj = GMRData()
        %obj = GMRData()
        %
        %   Output:
        %       obj     	:   GMRData instance
        %
        %   Input:
        %       filename    :   string  =   name of file with GMR data
        %
        %   Description:
        %   Constructor for GMRData class.  A constructor initializes all
        %   the variables in the class.
        %
        %   See also loadData, loadDataSet, loadDataSetRep, setFull  
        %   plotAll, plotOne, plotOne3D, plotCompare, plotDifference
        %   plotRBase
     
            obj.dataType = 0;
            obj.full = true;
            obj.nData = [];
            obj.f = [];
            obj.r = [];
            obj.rbase = [];
            obj.fA = [];
            obj.rA = [];
            obj.rbaseA = [];
            obj.c = {'b' 'g' 'r' 'm' 'c' 'k' 'b--' 'g--' 'r--' 'm--' 'c--' 'k--'};
        end
        
        function obj = loadData(obj,filename)
        %obj = loadData(obj, filename)
        %
        %   Output:
        %       obj         :   GMRData instance
        %
        %   Input:
        %       obj         :   GMRData instance
        %       filename    :   string  =   name of file with GMR data
        %
        %   Description:
        %   Load the GMR measurement data from the specified file.  
        %
        %   f is an [N R] matrix
        %   r is an [N R] matrix
        %   rbase is an [1 R] matrix
        %   fA is and [N 1] matrix
        %   rA is an [N 1] matrix
        %   rbaseA contains a single value
        %       N is the number of data points
        %       R is the number of repeated measurements
        %
        %   Example:
        %   tmp = GMRData();
        %   dat = tmp.loadData('GMR_S1.txt');
        %
        %   See also loadData, loadDataSet, loadDataSetRep, setFull  
        %   plotAll, plotOne, plotOne3D, plotCompare, plotDifference
            
            obj = obj.setDataType(0);
            dat = load(filename);

            if obj.full
                %Determine the indices of a full loop
                tmp = smooth(dat(:,1));
                len = length(tmp);
                ind1 = 1:floor(len/3);
                ind3 = floor(len/3*2):len;
                t1 = abs(tmp(ind1));
                i1 = find(t1 == max(t1));
                t3 = abs(tmp(ind3));
                i3 = find(t3 == max(t3))+floor(len/3*2);
                index = i1(1):i3(1);
                
                %Store data into class variables
                obj.f = dat(index,3:2:end);
                obj.r = dat(index,4:2:end);
                obj.fA = dat(index,1);
                obj.rA = dat(index,2);
                
                %Calculates base resistances
                obj.nData = size(obj.f,2);
                obj.rbase = zeros([1 obj.nData]);
                obj.rbaseA = 0;
                obj = obj.calcRBase();
            else
                %Store data into class variables
                obj.f = dat(:,3:2:end);
                obj.r = dat(:,4:2:end);
                obj.fA = dat(:,1);
                obj.rA = dat(:,2);
                
                %Calculates base resistances
                obj.nData = size(obj.f,2);
                obj.rbase = zeros([1 obj.nData]);
                obj.rbaseA = 0;
                obj = obj.calcRBase();
            end
        end
        
        function obj = loadDataSet(obj, basename)
        %obj = loadDataSet(obj, basename)
        %
        %   Output:
        %       obj         :   GMRData instance
        %
        %   Input:
        %       obj         :   GMRData instance
        %       basename    :   string  =   base name of the set of files
        %                                   containing GMR data
        %
        %   Description:
        %   The objective of this function is to load a single set of GMR
        %   measurements for more than one sensors.
        %
        %   Load a set of GMR measurement data from the directory.
        %   It is assumed that the name of the file contains a numeric
        %   identifier to distinguish the file.  For example a basename of
        %   'Blah_Blah_' is suitable to load the following files:
        %       Blah_Blah_2.txt
        %       Blah_Blah_5.txt
        %       Blah_Blah_12.txt
        %   For example, a basename of 'subfolder/Blah_' is suitable
        %   to load the following files
        %       subfolder/Blah_2.txt
        %       subfolder/Blah_5.txt
        %       subfolder/Blah_12.txt
        %
        %   f is an [N S R] matrix
        %   r is an [N S R] matrix
        %   rbase is an [S R] matrix
        %   fA is and [N S] matrix
        %   rA is an [N S] matrix
        %   rbaseA is an [S 1] matrix
        %       N is the number of data points
        %       S is the number of sensors
        %       R is the number of repeated measurements
        %
        %   Example:
        %   tmp = GMRData();
        %   dat = tmp.loadDataSet('Blah_Blah_');
        %
        %   See also loadData, loadDataSet, loadDataSetRep, setFull  
        %   plotAll, plotOne, plotOne3D, plotCompare, plotDifference
        %   plotRBase
            
            obj = obj.setDataType(1);
            [filename, sensorNumber] = obj.getDirListing(basename);
            dat = load(filename{1});
            nRep = size(dat,2)/2-1;
            obj.nData = [max(sensorNumber) nRep];
            
            if obj.full
                %Determine the indices of a full loop
                
                tmp = smooth(dat(:,1));
                len = length(tmp);
                ind1 = 1:floor(len/3);
                ind3 = floor(len/3*2):len;
                t1 = abs(tmp(ind1));
                i1 = find(t1 == max(t1));
                t3 = abs(tmp(ind3));
                i3 = find(t3 == max(t3))+floor(len/3*2);
                index = i1(1):i3(1);
                
                %Store data
                obj.f = zeros([length(index) obj.nData(1) obj.nData(2)]);
                obj.r = zeros([length(index) obj.nData(1) obj.nData(2)]);
                
                for i = 1:length(sensorNumber)
                    dat = load(filename{i});
                    obj.f(:,sensorNumber(i),:) = dat(index,3:2:end);
                    obj.r(:,sensorNumber(i),:) = dat(index,4:2:end);
                    obj.fA(:,sensorNumber(i)) = dat(index,1);
                    obj.rA(:,sensorNumber(i)) = dat(index,2);
                end
                
                %Calculates base resistance
                obj.rbase = zeros([obj.nData(1) obj.nData(2)]);
                obj.rbaseA = zeros([1 obj.nData(1)]);
                obj = obj.calcRBase();
            else
                %Store data
                %Store data
                obj.f = zeros([length(dat(:,1)) obj.nData(1) obj.nData(2)]);
                obj.r = zeros([length(dat(:,1)) obj.nData(1) obj.nData(2)]);
                
                for i = 1:length(sensorNumber)
                    dat = load(filename{i});
                    obj.f(:,sensorNumber(i),:) = dat(:,3:2:end);
                    obj.r(:,sensorNumber(i),:) = dat(:,4:2:end);
                    obj.fA(:,sensorNumber(i)) = dat(:,1);
                    obj.rA(:,sensorNumber(i)) = dat(:,2);
                end
                
                %Calculates base resistance
                obj.rbase = zeros([obj.nData(1) obj.nData(2)]);
                obj.rbaseA = zeros([1 obj.nData(1)]);
                obj = obj.calcRBase();
            end
        end
        
        function obj = loadDataSetRep(obj, basename)
        %obj = loadDataSetRep(obj, basename)
        %
        %   Output:
        %       obj         :   GMRData instance
        %
        %   Input:
        %       obj         :   GMRData instance
        %       basename    :   string  =   base name of the set of files
        %                                   containing GMR data
        %
        %   Description:
        %   The objective of this function is to load a series of sets of
        %   GMR measurements.  An example is a dataset generated from
        %   measuring the response of 12 sensors every 15 minutes for 2
        %   hours.
        %
        %   Load a set of GMR measurement data from the directory.
        %   It is assumed that the name of the file contains a numeric
        %   identifier to distinguish the file.  For example a basename of
        %   'Blah_Blah_' is suitable to load the following files:
        %       Blah_Blah_1_S2.txt
        %       Blah_Blah_2_S2.txt
        %       Blah_Blah_3_S2.txt
        %       Blah_Blah_1_S5.txt
        %       Blah_Blah_2_S5.txt
        %       Blah_Blah_3_S5.txt
        %
        %   f is an [N S R] matrix
        %   r is an [N S R] matrix
        %   rbase is an [S R] matrix
        %   fA is and [N S] matrix
        %   rA is an [N S] matrix
        %   rbaseA is an [S 1] matrix
        %       N is the number of data points
        %       S is the number of sensors
        %       R is the number of repeated measurements
        %
        %   Example:
        %   tmp = GMRData();
        %   data1 = tmp.loadDataSetRep('Blah_Blah_');
        %
        %   See also loadData, loadDataSet, loadDataSetRep, setFull  
        %   plotAll, plotOne, plotOne3D, plotCompare, plotDifference
        %   plotRBase
            
            obj = obj.setDataType(2);
            [filename, sensorNumber, repNumber] = obj.getDirListing(basename);
            nRep = length(unique(repNumber));
            obj.nData = [max(sensorNumber) nRep];
                
            if obj.full
                %Determine the indices of a full loop
                dat = load(filename{1});
                tmp = smooth(dat(:,1));
                len = length(tmp);
                ind1 = 1:floor(len/3);
                ind3 = floor(len/3*2):len;
                t1 = abs(tmp(ind1));
                i1 = find(t1 == max(t1));
                t3 = abs(tmp(ind3));
                i3 = find(t3 == max(t3))+floor(len/3*2);
                index = i1(1):i3(1);
                
                %Store data

                obj.f = zeros([length(index) obj.nData(1) obj.nData(2)]);
                obj.r = zeros([length(index) obj.nData(1) obj.nData(2)]);
                
                for i = 1:length(sensorNumber)
                    dat = load(filename{i});
                    obj.f(:,sensorNumber(i),repNumber(i)) = dat(index,1);
                    obj.r(:,sensorNumber(i),repNumber(i)) = dat(index,2);
                end
                obj.fA = mean(obj.f,3);
                obj.rA = mean(obj.r,3);
                
                %Calculates base resistance
                obj.rbase = zeros([obj.nData(1) obj.nData(2)]);
                obj.rbaseA = zeros([1 obj.nData(1)]);
                obj = obj.calcRBase();
            else
                %Store data
                for i = 1:length(sensorNumber)
                    dat = load(filename{i});
                    
                    obj.f(:,sensorNumber(i),repNumber(i)) = dat(:,1);
                    obj.r(:,sensorNumber(i),repNumber(i)) = dat(:,2);
                end
                obj.fA = mean(obj.f,3);
                obj.rA = mean(obj.r,3);

                %Calculates base resistances
                obj.rbase = zeros([obj.nData(1) obj.nData(2)]);
                obj.rbaseA = zeros([1 obj.nData(1)]);
                obj = obj.calcRBase();
            end
        end
        
        function obj = setFull(obj, full)
        %obj = setFull(obj, full)
        %
        %   Output:
        %       obj         :   GMRData instance
        %
        %   Input:
        %       obj         :   GMRData instance
        %       full        :   boolean
        %                       false   =   virgin loop
        %                       true    =   full loop
        %
        %   Description:
        %   Set the loop property between virgin loop or full loop.  A
        %   virgin loop is the raw measurement which measures from zero to
        %   +F to -F to +F to 0.  A full loop is a section from the virgin
        %   loop which contains data from +F to -F to +F.
        %
        %   Example:
        %   tmp = GMRData();
        %   tmp = tmp.setFull(false)
        %   dat1 = tmp.loadData('GMR_S1.txt');
        %   dat2 = tmp.loadData('GMR_S2.txt');
        %
        %   See also loadData, loadDataSet, loadDataSetRep, setFull  
        %   plotAll, plotOne, plotOne3D, plotCompare, plotDifference
        %   plotRBase
        
            if ~isa(full,'logical')
                error('GMRData:setFull','The input full must be a boolean.');
            end
        
            obj.full = full;  
        end
        
        function obj = plotAll(obj)
        %obj = plotAll(obj)
        %
        %   Output:
        %       obj         :   GMRData instance
        %
        %   Input:
        %       obj         :   GMRData instance
        %
        %   Description:
        %   Plots all the data for review.  The plot generated depends on
        %   the dataset:
        %   loadData:       Plot average data on top of all data.
        %   loadDataSet:    Plot all average data with different line
        %                   color and type.
        %   loadDataSetRep: Plot all repitition of each sensor in a 3x4
        %                   subplot
        %
        %   Example:
        %   tmp = GMRData();
        %   data1 = tmp.loadData('GMR_S1.txt');
        %   data1.plotAll();
        %
        %   See also loadData, loadDataSet, loadDataSetRep, setFull  
        %   plotAll, plotOne, plotOne3D, plotCompare, plotDifference
        %   plotRBase
        
            if obj.dataType == 0
                hold on
                for i = 1:obj.nData
                    plot(obj.f(:,i),obj.r(:,i)-obj.rbase(i),'b');
                end
                plot(obj.fA,obj.rA-obj.rbaseA,'r');
                hold off
                xlabel('f [Oe]','fontsize',14);
                ylabel('R [Ohm]','fontsize',14);
                rMean = obj.rbaseA;
                rSpread = std(obj.rbase)/rMean*100;
                title(['Rm(' num2str(rMean,3) ')  Spread(' num2str(rSpread,1) '%)']);
            elseif obj.dataType == 1
                for i = 1:obj.nData(1)
                    subplot(3,4,i)
                    hold on
                    for j = 1:obj.nData(2)
                        plot(obj.f(:,i,j),obj.r(:,i,j)-obj.rbase(i,j),'b');
                    end
                    plot(obj.fA(:,i),obj.rA(:,i)-obj.rbaseA(i),'r');
                    hold off
                    xlabel('f [Oe]','fontsize',14);
                    ylabel('R [Ohm]','fontsize',14);
                    rMean = obj.rbaseA(i);
                    rSpread = std(obj.rbase(i,:))/rMean*100;
                    title(['S#' num2str(i) ': Rm(' num2str(rMean,3) ')  Spr(' num2str(rSpread,1) '%)']);
                end
            elseif obj.dataType == 2
                for i = 1:obj.nData(1)
                    subplot(3,4,i)
                    hold on
                    for j = 1:obj.nData(2)
                        plot(obj.f(:,i,j),obj.r(:,i,j)-obj.rbase(i,j),'b');
                    end
                    plot(obj.fA(:,i),obj.rA(:,i)-obj.rbaseA(i),'r');
                    hold off
                    xlabel('f [Oe]','fontsize',14);
                    ylabel('R [Ohm]','fontsize',14);
                    rMean = obj.rbaseA(i);
                    rSpread = std(obj.rbase(i,:))/rMean*100;
                    title(['S#' num2str(i) ': Rm(' num2str(rMean,3) ')  Spr(' num2str(rSpread,1) '%)']);
                end
            end
        end
        
        function obj = plotOne(obj, varargin)
        %obj = plotOne(obj, sensorNumber*, avgOnly*)
        %
        %   Output:
        %       obj         	:   GMRData instance
        %
        %   Input:
        %       obj         	:   GMRData instance
        %       sensorNumber*   :   integer     =   sensor number to plot
        %       avgOnly         :   boolen      =   show average only
        %                                       =   false (default)
        %
        %   Description:
        %   Plot one dataset for review.  The plot generated depends on the
        %   dataset:
        %   loadData:       Plot average data
        %   loadDataSet:    Plot specified sensor, all data
        %   loadDataSetRep: Plot specified sensor, all data
        %
        %   Example:
        %   tmp = GMRData();
        %   dat = tmp.loadData('GMR_S1.txt');
        %   dat.plotAll();
        %
        %   See also loadData, loadDataSet, loadDataSetRep, setFull  
        %   plotAll, plotOne, plotOne3D, plotCompare, plotDifference
        %   plotRBase
        
            sensorNumber = 1;
            avgOnly = false;
            if nargin > 1
                sensorNumber = varargin{1};
            end
            if nargin > 2
                avgOnly = varargin{2};
            end
        
            if obj.dataType == 0
                plot(obj.fA,obj.rA-obj.rbaseA,'b');
                xlabel('f [Oe]','fontsize',14);
                ylabel('R [Ohm]','fontsize',14);
                rMean = obj.rbaseA;
                rSpread = std(obj.rbase)/rMean*100;
                title(['Rm(' num2str(rMean,3) ')  Spread(' num2str(rSpread,1) '%)']);
            elseif any(obj.dataType == [1,2])
                hold on
                if ~avgOnly
                    for j = 1:obj.nData(2)
                        plot(obj.f(:,sensorNumber,j),obj.r(:,sensorNumber,j)-obj.rbase(sensorNumber,j),'b');
                    end
                end
                plot(obj.fA(:,sensorNumber),obj.rA(:,sensorNumber)-obj.rbaseA(sensorNumber),'r');
                hold off
                xlabel('f [Oe]','fontsize',14);
                ylabel('R [Ohm]','fontsize',14);
                rMean = obj.rbaseA(sensorNumber);
                rSpread = std(obj.rbase(sensorNumber,:))/rMean*100;
                title(['S#' num2str(sensorNumber) ': Rm(' num2str(rMean,3) ')  Spr(' num2str(rSpread,1) '%)']);
            end
        end
        
        function obj = plotOne3D(obj, varargin)
        %obj = plotOne3D(obj, sensorNumber*)
        %
        %   Output:
        %       obj         	:   GMRData instance
        %
        %   Input:
        %       obj         	:   GMRData instance
        %       sensorNumber*   :   integer     =   sensor number to plot
        %
        %   Description:
        %   loadData:       Nothing
        %   loadDataSet:    Plot specified sensor in 3D
        %   loadDataSetRep: Plot specified sensor in 3D
        %
        %   Example:
        %   tmp = GMRData();
        %   tmp = tmp.setFull(true)
        %   data1 = tmp.loadData('GMR_S1.txt');
        %   data1.plotAll();
        %
        %   See also loadData, loadDataSet, loadDataSetRep, setFull  
        %   plotAll, plotOne, plotOne3D, plotCompare, plotDifference
        %   plotRBase
        
            sensorNumber = 1;
            if nargin > 1
                sensorNumber = varargin{1};
            end
        
            if obj.dataType == 0
                %Do Nothing
            elseif any(obj.dataType == [1, 2])
                z = ones(size(obj.f(:,1,1)));
                hold on
                for i = 1:size(obj.f,3)
                    plot3(z*i,obj.f(:,sensorNumber,i),obj.r(:,sensorNumber,i))
                end
                hold off
                view(3)
                xlabel('Trial [#]','fontsize',14);
                ylabel('Field [Oe]','fontsize',14);
                zlabel('R [Ohm]','fontsize',14);
            end
        end
        
        function obj = plotRBase(obj)
        %obj = plotRBase(obj)
        %
        %   Output:
        %       obj         :   GMRData instance
        %
        %   Input:
        %       obj         :   GMRData instance
        %
        %   Description:
        %   Plots the base resistance to track it's change over time
        %   loadData:       Plot base resistance of same sensor over time
        %   loadDataSet:    Plot base resistance of all sensors
        %   loadDataSetRep: Plot base resistance of all sensors over time
        %
        %   Example:
        %   tmp = GMRData();
        %   dat = tmp.loadData('GMR_S1.txt');
        %   dat.plotRBase();
        %
        %   See also loadData, loadDataSet, loadDataSetRep, setFull  
        %   plotAll, plotOne, plotOne3D, plotCompare, plotDifference
        %   plotRBase
        
            if obj.dataType == 0
                plot(obj.rbase,'b');
                xlabel('Data Number','fontsize',14);
                ylabel('R [Ohm]','fontsize',14);
                rMean = obj.rbaseA;
                rSpread = std(obj.rbase)/rMean*100;
                title(['Rm(' num2str(rMean,3) ')  Spread(' num2str(rSpread,1) '%)']);
            elseif any(obj.dataType == [1, 2])
                for i = 1:obj.nData(1)
                    subplot(3,4,i)
                    plot(obj.rbase(i,:),'b');
                    xlabel('Data Number','fontsize',14);
                    ylabel('R [Ohm]','fontsize',14);
                    rMean = obj.rbaseA(i);
                    rSpread = std(obj.rbase(i,:))/rMean*100;
                    title(['S#' num2str(i) ': Rm(' num2str(rMean,3) ')  Spr(' num2str(rSpread,1) '%)']);
                end
            end
        end
        
        function obj = plotDifference(obj)
        %obj = plotDifference(obj)
        %
        %   Output:
        %       obj         :   GMRData instance
        %
        %   Input:
        %       obj         :   GMRData instance
        %
        %   Description:
        %   Plots the difference between the average data and individual
        %   data.  Not sure how useful this is yet, but the idea is to
        %   quantify changes in the GMR loop.
        %   loadData:       Plot difference of (raw data) - (average data)
        %   loadDataSet:    Plot difference of data - average data for each
        %                   sensor
        %   loadDataSetRep: Plot difference of data - average data for each
        %                   sensor
        %
        %   Example:
        %   tmp = GMRData();
        %   dat = tmp.loadData('GMR_S1.txt');
        %   dat.plotDifference();
        %
        %   See also loadData, loadDataSet, loadDataSetRep, setFull  
        %   plotAll, plotOne, plotOne3D, plotCompare, plotDifference
        %   plotRBase
            
            if obj.dataType == 0
                resr = zeros([1 obj.nData]);
                for i = 1:obj.nData
                    resr(i) = sum(((obj.r(:,i)-obj.rbase(i))-(obj.rA-obj.rbaseA)).^2);
                end
                plot(resr,'b');
                xlabel('Trial [#]','fontsize',14);
                ylabel('Residual [a.u.]','fontsize',14);
            elseif any(obj.dataType == [1, 2])
                resr = zeros([1 obj.nData(2)]);
                for i = 1:obj.nData(1)
                    subplot(3,4,i)
                    hold on
                    for j = 1:obj.nData(2)
                        resr(j) = sum(((obj.r(:,i,j)-obj.rbase(i,j))-(obj.rA(:,i)-obj.rbaseA(i))).^2);
                    end
                    plot(resr,'b');
                    hold off
                    xlabel('Trial [#]','fontsize',14);
                    ylabel('Residual [a.u.]','fontsize',14);
                end
            end
        end
        
        function obj = plotCompare(obj)
        %obj = plotCompare(obj)
        %
        %   Output:
        %       obj         :   GMRData instance
        %
        %   Input:
        %       obj         :   GMRData instance
        %
        %   Description:
        %   Plots all the different sensor data on a single figure for
        %   comparison.
        %   loadData:       Plot average data
        %   loadDataSet:    Plot average data for all sensors
        %   loadDataSetRep: Plot average data for all sensors
        %
        %   Example:
        %   tmp = GMRData();
        %   dat = tmp.loadDataSet('GMR_S');
        %   dat.plotCompare();
        %
        %   See also loadData, loadDataSet, loadDataSetRep, setFull  
        %   plotAll, plotOne, plotOne3D, plotCompare, plotDifference
        %   plotRBase
        
            if obj.dataType == 0
                hold on
                plot(obj.fA,obj.rA,'r');
                hold off
                xlabel('f [Oe]','fontsize',14);
                ylabel('R [Ohm]','fontsize',14);
            elseif any(obj.dataType == [1, 2])
                leg = {};
                for i = 1:obj.nData(1)
                    hold on
                    plot(obj.fA(:,i),obj.rA(:,i)-obj.rbaseA(i),obj.c{i});
                    leg{end+1} = ['S' num2str(i)];
                    hold off
                    xlabel('f [Oe]','fontsize',14);
                    ylabel('R [Ohm]','fontsize',14);
                    rMean = obj.rbaseA(i);
                    rSpread = std(obj.rbase(i,:))/rMean*100;
                    title(['S#' num2str(i) ': Rm(' num2str(rMean,3) ')  Spr(' num2str(rSpread,1) '%)']);
                    legend(leg,'location','ne')
                end
            end
        end
    end
    
    methods (Access=protected, Sealed=true)
        
        function obj = setDataType(obj,type)
        %obj = setDataType(obj, type)
        %
        %   Output:
        %       obj         :   GMRData instance
        %
        %   Input:
        %       obj         :   GMRData instance
        %       type        :   integer     =   identifier for data type
        %                       0           =   Single GMR sensor
        %                       1           =   Set of GMR sensor
        %                       2           =   Repeated sets of GMR sensor
        %
        %   Description:
        %   Set the data type parameter.  The data type parameter is used
        %   to specify how other functions interact with the data.
            
            if any(type == [0, 1, 2])
                obj.dataType = type;
            else
                error('GMRData:setDataType','The type parameter must be in the set [0 1].');
            end
        end
        
        function obj = calcRBase(obj)
        %obj = calcRBase(obj)
        %
        %   Output:
        %       obj         =   GMRData instance
        %
        %   Input:
        %       obj         =   GMRData instance
        %
        %   Description:
        %   Determine the base resistance.  The base resistance is the
        %   lowest resistance in the measurement.
            
            if obj.dataType == 0
                for i = 1:obj.nData
                    obj.rbase(i) = min(smooth(obj.r(:,i)));
                end
                obj.rbaseA = min(obj.rA);
            elseif any(obj.dataType == [1, 2])
                for i = 1:obj.nData(1)
                    for j = 1:obj.nData(2)
                        obj.rbase(i,j) = min(smooth(obj.r(:,i,j)));
                    end
                    obj.rbaseA(i) = min(obj.rA(:,i));
                end
            end
        end
        
        function [filename, sensorNumber, repNumber] = getDirListing(~, basename)
        %obj = getDirListing(obj, basename)
        %
        %   Output:
        %       obj         =   GMRData instance
        %
        %   Input:
        %       obj         =   GMRData instance
        %       basename    :   string  =   base name of files
        %
        %   Description:
        %   Determines the filenames which contain the specified basename.
        %   The base name may also include a directory.  For example, the
        %   following basenames are valid:
        %       'Blah_blah'
        %       'AFolder/Blah_blah'
        
            directory = '';
            iSlash = find(basename == '/');
            if ~isempty(iSlash)
                directory = basename(1:iSlash(end));
                basename = basename(iSlash(end)+1:end);
            end
            
            listing = dir([directory basename '*.txt']);
            
            filename = cell([1 length(listing)]);
            sensorNumber = zeros([0 length(filename)]);
            repNumber = sensorNumber;
            for i = 1:length(listing)
                filename{i} = [directory listing(i).name];
                tmp = regexp(filename{i},'\d+','match');
                sensorNumber(i) = str2double(tmp(end));
                try
                    repNumber(i) = str2double(tmp(end-1));
                catch
                    continue;
                end
            end
            
        end
        
    end
    
end