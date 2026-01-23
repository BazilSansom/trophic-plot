function [W, meta] = tflExample(name)

switch lower(name)

    case 'dyad'
        % a -> b
        W = [0 1 ;
             0 0 ];

        meta.regime = 1;
        meta.description = 'Dyad';
        meta.q_expected = '0';
    
    case 'star'
        % a -> b, a -> c
        W = [0 1 1 ;
             0 0 0;
             0 0 0 ];

        meta.regime = 1;
        meta.description = 'minimal star topology';
        meta.q_expected = '0';

    case 'chain'
        % a -> b -> c
        W = [0 1 0;
             0 0 1;
             0 0 0];

        meta.regime = 1;
        meta.description = 'Simple chain';
        meta.q_expected = '0';

    case 'ffl'
        % a -> b -> c, a -> c
        W = [0 1 1;
             0 0 1;
             0 0 0];

        meta.regime = 2;
        meta.description = 'Simple feed-forward loop';
        meta.q_expected = '> 0';

    case 'fbl3_circulation'
        % a -> b -> c -> a 
        W = [0 1 0;
             0 0 1;
             1 0 0];

        meta.regime = 3;
        meta.description = 'Simple 3 node feedback loop';
        meta.q_expected = '1';

    
    case 'fbl3'
        % a -> b -> c -> a 
        W = [0 2 0;
             0 0 2;
             1 0 0];

        meta.regime = 3;
        meta.description = 'Simple 3 node feedback loop with some net flow';
        meta.q_expected = '>1';
    
    case 'coherent'
        
        E=[1,4; 1,5; 2,5; 3,5; 3,6; 4,7; 5,7; 5,8; 5,9; 6,9];
        W = edgelistToAdjacency(E);

        meta.regime = 1;
        meta.description = 'Small perfectly coherent network (3 layrers)';
        meta.q_expected = '0';

    case 'semi-coherent'
        
        E=[1,4; 1,5; 2,5; 3,5; 3,6; 4,7; 5,7; 5,8; 5,9; 6,9; 1,7; 2,8; 3,9; 1,8;2,7;2,9];
        W = edgelistToAdjacency(E);

        meta.regime = 2;
        meta.description = 'Small semi-coherent network (3 layrers)';
        meta.q_expected = '>0';

    case 'hierarchical_feedback_modules'
        % Two strongly connected components with unidirectional coupling
        E=[2,1; 3,1; 1,2; 3,2; 1,3; 2,3; 5,4; 6,4; 3,5; 4,5; 6,5;4,6;5,6];
        W = edgelistToAdjacency(E);
        
        meta.regime = 3;
        meta.description = 'Two strongly connected components with unidirectional coupling';
        meta.q_expected = '>0';

    case 'dense_scc_tail'
        % A large strongly connected component with a tail
        data = load(fullfile(fileparts(mfilename('fullpath')), 'badcase_denseSCC_tail.mat'));
        W = data.W_bad;
        
        meta.regime = 3;
        meta.description = 'A large strongly connected component with a tail';
        meta.q_expected = '>0';

    case

    otherwise
        error('Unknown example "%s".', name);
end
