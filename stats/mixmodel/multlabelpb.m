classdef multlabelpb
    % The abstract base class to represent a multi-class label distribution 
    % inference problem
    %
    
    % Created by Dahua Lin, on April 21, 2010
    %
    
    methods
        Q = infer_q(obj, L)
        % Infer the probabilities of labels given the likelihoods
        %
        %   Q = obj.infer_q(L);
        %       infer the probabilities of labels given the likelihoods.
        %       Here, L is the likelihood matrix. Suppose there are 
        %       m models and n samples, then L should be a real matrix
        %       of size m x n, and L(i, j) is the likelihood of the 
        %       j-th sample with respect to the i-th model.
        %
        %       In the output, Q is also a m x n real matrix, where 
        %       Q(i, j) is the inferred posterior probability of 
        %       that the j-th sample is generated from the i-th model.
        %
        
        v = eval_qscore(obj, Q)
        % Evaluates the score of the labeling given by Q
        %
        %   v = obj.eval_qscore(Q);
        %       Evaluates the score of labeling given by Q. The score
        %       depends on the particular labeling model. Generally,
        %       a greater score indicates that q is more consistent
        %       with the prior labeling model.
        %
        
    end
    
end