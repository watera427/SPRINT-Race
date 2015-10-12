SPRINT_Race is a simple implementation of a multi-objective racing 
procedure based on Sequential Probability Ratio Test with Indifference Zone.

For details, please refer to the following paper:

T. Zhang, M. Georgiopoulos, G. C. Anagnostopoulos, "Pareto Optimal Model 
Selection via SPRINT-Race", Under Review 

SPRINT_Race is available at https://github.com/watera427/SPRINT-Race_v2

Author contact: Tiantian Zhang 

Email: zhangtt@knights.ucf.edu

COPYRIGHT (C) 2015 Tiantian Zhang, Machine Learning Lab, University of Central Florida
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are 
met:

+ Redistributions of source code must retain the above copyright 
  notice, this list of conditions and the following disclaimer.
+ Redistributions in binary form must reproduce the above copyright 
  notice, this list of conditions and the following disclaimer in 
  the documentation and/or other materials provided with the distribution
      
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGE.

=======================================================================================



SPRINT_Race method:
================

+ obj = SPRINT_Race(M, D, alpha, beta, delta)

  Create a S_Race object before starting a racing

  ' M ' is the number of initial models
  ' D ' is the number of objectives
  ' alpha ' is the overall Type I of SPRINT-Race
  ' beta ' is the overall Type II of SPRINT-Race
  ' delta ' is the indifference zone parameter

+ dataAcquisition(obj)
  
  Generate performance vectors needed at each step of racing
  
  
+ Racing(obj)
  
  The SPRINT-Race procedure. Once started, it will continue racing until no
  more comparison and sampling is needed. The final set of non-dominated 
  models returned by SPRINT-Race will be stored in obj.models.

+ [w1, w2] = dominates(~,data1, data2)

  Count the number of dominance
  
  ' data1 ' is the performance vector of the first model
  ' data2 ' is the performance vector of the second model
  ' w1 ' is 1 if data1 dominates data2; 0 otherwise.
  ' w2 ' is 1 if data2 dominates data1; 0 otherwise.  

  
S_Race variables:
==================
        M               % number of initial models
        D               % number of objectives
        toDel           % index of the models have been removed so far
        alpha           % the overall Type I of SPRINT-Race
        beta            % the overall Type II of SPRINT-Race
        delta           % the parameter of indifference zone
        models          % the final set of non-dominated models returned by SPRINT-Race
        no_samples      % number of samples used 
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


Tips on Practical Use:
==================
+ Please modify function 'dataAcuquisition' based on your needs
+ SPRINT_Race is designed for multi-objective maximization problem only
+ Do not reinitialize SPRINT_Race object before racing is complete
+ Remember to add path in default search path ('STARTUP.m')
+ Author used MATLAB R2013b. 



Examples:
==================
Please refer to 'example.m'





