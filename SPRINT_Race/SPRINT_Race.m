classdef SPRINT_Race < handle
% SPRINT_Race class

% COPYRIGHT
%  (C) 2015 Tiantian Zhang, Machine Learning Lab, University of Central Florida
% zhangtt@knights.ucf.edu
% Reference
% T. Zhang, M. Georgiopoulos, G. C. Anagnostopoulos, "Pareto Optimal Model
% Selection via SPRINT-Race", Under Review
    
    %% Public Properties
    properties
        M               % number of initial models
        D               % number of objectives
        toDel           % index of the models have been removed so far
        alpha           % overall Type I of SPRINT-Race
        beta            % overall Type II of SPRINT-Race
        delta           % the parameter of indifference zone
        models          % the final set of non-dominated models returned by SPRINT-Race
        no_samples      % number of samples used 
    end

    %% Protected Properties
    properties(GetAccess = protected, SetAccess = protected)
        stop            % the termination matrix, stop(i,j,1) indicates if the SPRT_1
                        % between the i-th model and j-th model is active
                        % or not; stop(i,j,2) indicates if the SPRT_2 is
                        % active or not. (0 - active; 1 - terminated)
        dom             % the dominance matrix, dom(i,j) is the number that the i-th
                        % model dominates the j-th model
        t               % the step index of SPRINT-Race
        data            % store the performance vectors needed at the current step 
        As              % the set of lower boundary values obtained via Sequential Holm's FWER control method
        Bs              % the set of upper boundary values obtained via Sequential Holm's FWER control method
        a               % the number of null hypotheses have been accepted
        r               % the number of null hypotheses have been rejected
    end
    
    %% Public Methods
    methods
        
        % Constructor
        function obj = SPRINT_Race(M, D, alpha, beta, delta)
            if nargin == 5
                obj.M = M;
                obj.D = D;
                obj.alpha = alpha;
                obj.beta = beta;
                obj.delta = delta;
                obj.models = 1:obj.M;
                obj.toDel = [];
                obj.stop = zeros(M, M, 2);
                obj.stop(:,:,1) = eye(M);
                obj.stop(:,:,2) = eye(M);
                obj.dom = zeros(M, M);
                obj.data = cell(1,M);
                obj.As = zeros(1,M * (M - 1));
                obj.Bs = zeros(1,M * (M - 1));
                obj.a = 0;
                obj.r = 0;
                obj.no_samples = 0;
            else
                %error('Too few/many arguments');
            end
        end % SPRINT_Race
        
        function dataAcquisition(obj)
            % obtain the performance vectors needed for the current step
            % this function is user-specified
            % in this illustration, we simply generate performance vectors
            % from a set of predefined multivariate Gaussian distributions
            obj.data = cell(1,obj.M);
            load('gaussians.mat'); % load the matrix containing the means and variances of M 
            for i = 1:obj.M
                for j = i + 1:obj.M
                    if(min(obj.stop(i,j,:)) == 0) % meaning that performance vectors of the i-th and j-th model are needed
                        if isempty(obj.data{i})
                            obj.data{i} = mvnrnd(gaussians{i,1},gaussians{i,2},1);
                        end
                        if isempty(obj.data{j})
                            obj.data{j} = mvnrnd(gaussians{j,1},gaussians{j,2},1);
                        end                        
                    end
                end
            end            
        end % dataAcquisition                   
        
        % SPRINT-Race
        function Racing(obj)
            compCritical(obj);
            while(min(min(min(obj.stop)) == 0)) % at least one test is active
                dataAcquisition(obj);
                for i = 1:obj.M
                    for j = i + 1:obj.M
                        if(min(obj.stop(i,j,:)) == 0) % start the test
                            obj.no_samples = obj.no_samples + 1;
                            % data acquisition                            
                            [w1, w2] = obj.dominates(obj.data{i}, obj.data{j});
                            obj.dom(i,j) = obj.dom(i,j) + w1;
                            obj.dom(j,i) = obj.dom(j,i) + w2;
                            w1 = obj.dom(i,j);
                            w2 = obj.dom(j,i);                            
                            % compute the test statistics
                            t1 = w1 * log(1/(1 - 2 * obj.delta)) - w2 * log(1 + 2 * obj.delta);
                            t2 = w1 * log(1 + 2 * obj.delta) - w2 * log(1/(1 - 2 * obj.delta));
                            if obj.stop(i,j,1) == 0 && t1 <= obj.As(obj.a + 1)
                                obj.stop(i,j,1) = 1;
                                obj.a = obj.a + 1;
                            end
                            if obj.stop(i,j,1) == 0 && t1 >= obj.Bs(obj.r + 1)
                                obj.stop(i,j,1) = 2;
                                obj.r = obj.r + 1;
                            end
                            if obj.stop(i,j,2) == 0 && t2 <= obj.As(obj.a + 1)
                                obj.stop(i,j,2) = 1;
                                obj.a = obj.a + 1;
                            end
                            if obj.stop(i,j,2) == 0 && t2 >= obj.Bs(obj.r + 1)
                                obj.stop(i,j,2) = 2;
                                obj.r = obj.r + 1;
                            end
                            % see if any model is eliminated 
                            if(obj.stop(i,j,1) == 1 && obj.stop(i,j,2) == 1) % i is removed from racing
                                obj.toDel = [obj.toDel,i];
                                obj.stop(i,:,1) = ones(1,obj.M); % stop all the tests involving model i
                                obj.stop(i,:,2) = ones(1,obj.M); 
                                obj.stop(:,i,1) = ones(obj.M,1);
                                obj.stop(:,i,2) = ones(obj.M,1);
                            elseif(obj.stop(i,j,1) == 2 && obj.stop(i,j,2) == 2) % j is removed from racing
                                obj.toDel = [obj.toDel,j];
                                obj.stop(j,:,1) = ones(1,obj.M); % stop all the tests involving model j
                                obj.stop(j,:,2) = ones(1,obj.M); 
                                obj.stop(:,j,1) = ones(obj.M,1);
                                obj.stop(:,j,2) = ones(obj.M,1);
                            end                           
                        end
                    end
                end                
            end
            % find out the models that are returned by sequential racing
            obj.toDel = unique(obj.toDel);
            obj.models(obj.toDel) = [];    
        end % Racing
                
    end
    
    
    methods (Access = private)
        function compCritical(obj)
            % Compute the critical values used in Sequential Holm's Step-Down
            % procedure. Please refers to paper "Sequential tests of multiple 
            % hypotheses controlling type I and II family-wise error rate (Jay 
            % Bartroff, 2014)" for detail
            % output arguments
            % @obj.As - a 1 x m vector containing the lower bounds
            % @obj.Bs - a 1 x m vector containing the upper bounds
            % input arguments
            % @obj.M - total number of models
            % @obj.alpha - overall Type I error probability
            % @obj.beta - overall Type II error probability
            m = obj.M * (obj.M - 1);
            for s = 1 : m
                tmp_a = (m - s + 1 - obj.beta) * obj.alpha/(m - s + 1)/(m - obj.beta);
                tmp_b = (m - s + 1 - obj.alpha) * obj.beta/(m - s + 1)/(m - obj.alpha);
                obj.As(s) = log(obj.beta / (1 - tmp_a) / (m - s + 1));
                obj.Bs(s) = log((1 - tmp_b) * (m - s + 1) / obj.alpha);
            end        
        end % compCritical
        
        % find the number of dominance 
        function [w1, w2] = dominates(~, data1, data2)
            w1 = 0;
            w2 = 0;
            tmp = data1 - data2;
            % for maximization problem
            if min(tmp) >=0 && max(tmp) > 0
                w1 = w1 + 1;
            elseif max(tmp) <=0 && min(tmp) < 0
                w2 = w2 + 1;
            end
        end % dominates        
       
    end
    
end