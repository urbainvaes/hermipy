function [out,varargout] = cont(hFun,u0,lambda0,Lf,varargin)
% CONT Arclength continuation with Keller and Moore-Penrose algorithms
%   CONT(HFUN,U0,LAMBDA0,LF,VARARGIN) performs arclength continuation on
%   a system of algebraic equations given by HFUN(U,LAMBDA) where LAMBDA
%   is the continuation parameter. Here U0 and LAMBDA0 are the initial
%   lambda and a guess to initialize continuation. LF corresponds to the
%   final length of the curve. VARARGIN can be one of the following
%   optional arguments, paired with their corresponding parameters.
%
%       'MonitorFcn'    :   Monitor function (e.g. the bifurcation curve).
%       'Step'          :   Step along the arclength
%       'Algorithm'     :   Continuation algorithm. Two options are allowed
%                           'Keller' and 'Moore-Penrose'. {'Keller'}
%       'Jacobian'      :   If set to 'on' means that the objective
%                           function outputs  the analytic Jacobian.{'off'}
%       'MaxIter'       :   Maximum number of Newton iterations.
%       'TolFun'        :   Tolerance on the residual of the norm {1e-6}.
%       'direction'     :   Whether to move forward (+1) or backward(-1)
%                           on the continuation curve.
%       'MaxStep'       :   The maximum allowed step {inf}.
%       'MinStep'       :   The minumum allowed step {eps}.
%       'Adaptive'      :   Step adaptation. {[]}
%                           A 3 or 4-element vector. {Recommended values}
%                            1. Number of Iterations below which the step
%                               size is extended. {length of U0}
%                            2. Extension/reduction factor {Must be >1}
%                            3. Norm of |U_n - U_(n-1)| above which ds is
%                               reduced {1e-n, n positive}
%                            4. (Optional) Cosine of angle between tangent
%                               vectors dot(V_n,V_(n-1)).
%                               {Keller, 0; Moore-Penrose, 0.9+}
%       'Stop'          :   Function handle to terminate continuation when
%                           the value returned becomes 0. If set to 'on' it
%                           creates a stop button in case a monitor
%                           function has any graphical output.
%       'OutputSel'     :   A vector specifying which components of the
%                           solution are returned.
%
%   The outputs depend on whether a MonitorFcn is specified.
%
%   [LAMBDA,U,STAB] = CONT(...) outputs the parameter LAMBDA, the
%       corresponding solution, U, and the 'stability' of the solutions
%       given by STAB, depending on the real part of the eigenvalues of the
%       Jacobian:
%         +1 : All positive
%         -1 : All negative
%          0 : Otherwise
%
%   [LAMBDA,MON,U,STAB] = CONT(...) with the 'MonitorFcn' property set
%       to a function handle, returns the abobe, while also returning the
%       outputs of the monitor function, MON.
%
% (c)COPYRIGHT v.1.0 February 2011 Nikos Savva

    % ARGUMENT CHECKS
    error(nargchk(0,0,mod(nargin,2)));
    error(nargoutchk(2,4,nargout));
    lastwarn('');
    hStopBtn = [];

    epsilon = sqrt(eps);
    % DEFAULT CONTINUATION OPTIONS
    opts = struct('monitorfcn',[],'step',0.1,'adaptive',[],...
        'maxstep',inf,'minstep',eps,'jacobian',[],'algorithm','keller',...
        'maxiter',4*length(u0),'direction',1,'tolfun',1e-10,...
        'stop',[],'outputsel',1:length(u0));

    % READ CONTINUATION OPTIONS
    for i=1:2:nargin-4,
        opts.(lower(varargin{i}))=varargin{i+1};
    end

    ds = opts.step;
    hMonitor = opts.monitorfcn;
    TolFun = opts.tolfun;
    maxiter = opts.maxiter;
    OutSel = opts.outputsel;
    Nsteps = floor(Lf/ds) + 1;
    N = length(u0);
    iStep = 1;
    lambda = lambda0;

    % CHECK ADAPTIVITY
    if isempty(opts.adaptive)
        isAdaptive = 0;
    else
        isAdaptive = 1;
        Nthresh = opts.adaptive(1);
        f_ds = opts.adaptive(2);
        tolX = opts.adaptive(3);

        try
            tolTheta = opts.adaptive(4);
        catch %#ok
            switch lower(opts.algorithm)
                case 'keller'
                    tolTheta = 0;
                case 'moore-penrose'
                    tolTheta = 0.99;
            end
        end
    end

    % STOPPING CRITERION
    if ~isempty(opts.stop),
       isStop = 1;
       hStop = opts.stop;
    else
       isStop = 0;
    end

    % CHECK JACOBIAN INPUT
    if strcmp(opts.jacobian,'on')
        hFunJac = @FunJacAna;
    else
        hFunJac = @FunJac;
    end

    % CHOOSE OUTPUTS
    outLam = zeros(Nsteps,1);

    if nargout>=2
        if isempty(hMonitor),
            outU = zeros(length(OutSel),Nsteps);
            outMon = [];
        else
            outMon = zeros(Nsteps,1);
            outU = [];
        end
    end

    if (nargout==3 && isempty(hMonitor)) || nargout==4,
        outStab = zeros(Nsteps,1);
    else
        outStab = [];
        outU = zeros(length(OutSel),Nsteps);
    end

    if nargout==4,
       outU = zeros(length(OutSel),Nsteps);
    end

    % SOLVE FOR THE FIRST TIME & GET THE JACOBIAN
    [f,Jac] = hFunJac(u0,lambda0);
    isOK = 0;

    for i = 1:maxiter,
        % CHECK CONVERGENCE
        if(norm(f,2) <= TolFun), isOK = 1; break; end;

        % NEWTON STEP
        u0 = u0 - (Jac(:,1:N)\f);

        % NEW FUNCTION & JACOBIAN VALUES
        [f,Jac] = hFunJac(u0,lambda0);
    end

    if ~isOK, error('MATLAB:Cont','Bad initial condition'); end
    u0 = [u0;lambda0]; u = u0;
    AssignOutputs;

    % NEW DIRECTION
    dz_ds0 = null(Jac);
    dz_ds0 = (opts.direction*sign(dz_ds0(end)+eps))*dz_ds0;
    S = Nsteps; iStep = 1; L = 0; reject = 0;

    switch lower(opts.algorithm),
        case 'moore-penrose'
            while L<Lf,
                if ~reject, iStep = iStep + 1; end;
                % INITIALIZE
                isTerminal = 0; reject = 1; isOK = 0; lastwarn('');
                dz_ds = dz_ds0;

                % PREDICTOR STEP
                u = u0 + dz_ds0*ds;

                % NEWTON ITERATION
                for NewtonStep = 1:maxiter,

                    % EVALUATE F & JACOBIAN
                    [f,Jac] = hFunJac(u(1:end-1),u(end));

                    % CHECK CONVERGENCE
                    if(norm(f,2) <= TolFun), isOK = 1; break; end

                    % CORRECTOR
                    delta = -[Jac ; dz_ds']\[f Jac*dz_ds; 0 0];
                    u = u + delta(:,1);
                    dz_ds = dz_ds + delta(:,2);
                    dz_ds = dz_ds/norm(dz_ds,2);

                    % CHECK IF WARNING IS ISSUED
                    if ~isempty(lastwarn),
                        lastwarn('');
                        break;
                    end

                    % MAKE SURE WE MOVE IN THE CORRECT DIRECTION
                    if dot(dz_ds0,dz_ds)<0, dz_ds = - dz_ds; end
                end

                % POST PROCESSING STEP
                if PostProcessing, break; end;
            end
        case 'keller'
            % RHS OF DIRECTION ESTIMATION
            RHS = zeros(N+1,1); RHS(N+1) = 1;

            while L<Lf,
                if ~reject,
                    iStep = iStep + 1;
                else
                    [f,Jac] = hFunJac(u0(1:end-1),u0(end)); %#ok
                end

                % INITIALIZE
                isTerminal = 0; reject = 1; isOK = 0; lastwarn('');
                dz_ds = dz_ds0;

                % COMPUTE DIRECTION ALONG HYPERPLANE
                v = [Jac;dz_ds']\RHS;
                dz_ds = v/norm(v,2);

                % MAKE SURE WE MOVE IN THE CORRECT DIRECTION
                if dot(u - u0,dz_ds)<0,
                    dz_ds = - dz_ds;
                end

                % PREDICTOR STEP
                u = u0 + dz_ds*ds;

                % EVALUATE F & JACOBIAN
                [f,Jac] = hFunJac(u(1:end-1),u(end));

                % NEWTON ITERATION
                for NewtonStep = 1:maxiter,
                    % CHECK CONVERGENCE
                    if(norm(f,2) <= TolFun), isOK = 1; break; end

                    % INVERSION
                    delta = [Jac;dz_ds']\[-f; ds - dot(u-u0,dz_ds)];
                    u = u + delta;

                    % EVALUATE F & JACOBIAN
                    [f,Jac] = hFunJac(u(1:end-1),u(end));
                end

                % POST PROCESSING STEP
                if PostProcessing, break; end;
            end
        otherwise
            error('MATLAB:CONT',['Invalid algorithm.',...
                ' Only ''Moore-Penrose'' and ''Keller'' are acceptable']);
    end

    % PREPARE OUTPUTS TO EXIT
    out = outLam(1:iStep);

    if nargout>=2 && ~isempty(hMonitor),
        varargout{1} = outMon(1:iStep);
    else
        varargout{1} = outU(:,1:iStep);
    end

    if (nargout==3 && isempty(hMonitor))
        varargout{2} = outStab(1:iStep);
    else
        varargout{2} = outU(:,1:iStep);
    end

    if nargout==4,
        varargout{3} = outStab(1:iStep);
    end

    % FUNCTION ASSIGNMENT
    function AssignOutputs
       outLam(iStep) = u(end);
       if ~isempty(outU), outU(:,iStep) = u(OutSel); end
       if ~isempty(outStab), outStab(iStep) = isStable(Jac); end
       if ~isempty(outMon), outMon(iStep) = hMonitor(u(1:end-1),u(end));end
    end

    % JACOBIAN IS ESTIMATED
    function [FF,JJ] = FunJac(x,lambda)
       FF = hFun(x,lambda);

       % NUMERICAL JACOBIAN
       x = [x;lambda];
       JJ = zeros(N,N+1);
       for n=1:N+1,
           x_new = x;  x_new(n) = x(n) + epsilon;
           f_new = hFun(x_new(1:N),x_new(N+1));
           JJ(:,n) = (f_new-FF)/epsilon;
       end
    end

    % JACOBIAN IS SUPPLIED
    function [FF,JJ] = FunJacAna(x,lambda)
        [FF,JJ] = hFun(x,lambda);

        if size(JJ,2)==N,
           JJ = [JJ, (hFun(x,lambda + epsilon) - FF)/epsilon];
        end
    end

    % POST PROCESSING
    function out = PostProcessing
        out = 0;
        ds0 = ds;

        % ADAPTIVE STEP
        if isAdaptive,
            m = dot(dz_ds,dz_ds0);

            if ((m<tolTheta || norm(u-u0,2)>tolX) && isOK) || ~isOK,
                ds = ds/f_ds;
                if ds < opts.minstep,
                    isTerminal = 1;
                    warning('MATLAB:Continuation',...
                        'Smallest step size reached');
                end
            elseif isOK && NewtonStep<Nthresh,
                ds = min(ds*f_ds,opts.maxstep);
                reject = 0;
            else
                reject = 0;
            end
        elseif isOK,
            reject = 0;
        else
            isTerminal = 1;
            reject = 1;
        end

        % STORE IF NOT REJECTED
        if ~reject,
            AssignOutputs;
            u0 = u;
            dz_ds0 = dz_ds;
            L = L + ds0;
        end

        % TERMINATE ON STOP
        if isStop
            if isa(hStop,'function_handle') && hStop(u,lambda)
                isTerminal = 1;
            elseif ~isempty(hMonitor),
                hFig = findobj('type','figure');

                if ~isempty(hFig) && isempty(hStopBtn),
                    delete(findobj('tag','ContStop'));
                    pos = get(0,'DefaultUicontrolPosition');
                    pos(1) = pos(1) - 15;
                    pos(2) = pos(2) - 15;
                    hStopBtn = uicontrol('Parent',hFig,...
                          'Style','pushbutton', ...
                          'Callback',@StopBtnCallback,...
                          'String','Stop', ...
                          'Tag','ContStop',...
                          'Position',pos);
                    set(hFig,'UserData',0);
                elseif get(hStopBtn,'UserData'),
                    isTerminal = 1;
                end
            end
        end

        % TERMINATE
        if isTerminal,
            iStep = iStep-(~isOK);
            warning('MATLAB:Cont','Continuation terminated.');
            out = 1;
        end

        % RE-ALLOCATE OUTPUTS IF ARRAY GROWS
        if iStep>S,
            S = S+round(2*(Lf-L)/L*iStep);
            temp = outLam;
            outLam = zeros(S,1);
            outLam(1:iStep) = temp;

            if ~isempty(outU),
                temp = outU;
                outU = zeros(length(OutSel),S);
                outU(:,1:iStep) = temp;
            end

            if ~isempty(outStab),
                temp = outStab;
                outStab = zeros(S,1);
                outStab(1:iStep) = temp;
            end

            if ~isempty(outMon),
                temp = outMon;
                outMon = zeros(S,1);
                outMon(1:iStep) = temp;
            end
        end
    end
end

% STOP BUTTON
function StopBtnCallback(gco,varargin)
    set(gco,'UserData',1);
end

% DETERMINE STABILITY
function out = isStable(Jac)
    eigvals = real(eig(Jac));
    if all(eigvals>0), out = 1;
    elseif all(eigvals<0), out = -1;
    else out = 0;
    end
end
