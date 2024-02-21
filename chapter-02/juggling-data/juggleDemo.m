%% Analyses of juggling data including registration. 
%%
%  The juggling data are from ten sequences, records or trials of a juggler 
%  juggling three balls.  
%  The measured variable is the position in three-space of the tip
%  of the forefinger in his right hand, taken 200 times per second:  
%    The X-coordinate is the horizontal position in the frontal plane,
%    that is, across the front of the body, measured from left to right.
%    The Y-coordinate is the position in the sagittal plane, that is the
%    body viewed from the side, larger values corresponding to away from
%    the body.
%    The Z-cocordinate is vertical.
%  Each of the ten sequences consists of from 11 to 13 cycles, a cycle
%  consisting of the ball being thrown into the air, the hand moving in
%  an arc across the body and down in preparation for receiving the 
%  incoming ball from the left hand, and then the transporting of the
%  new ball into position for its throw.  The total number of these
%  cycles across the ten records is 123, and the average duration of the
%  cycles is 713 milliseconds.

%  This Matlab file is designed to display the juggling data in various 
%  ways, including displaying the velocity and the acceleration of the
%  forefinger in addition to its position.  

%  The duration of the juggling cycles varies from throw to throw as the
%  juggler adjusts to variations in the trajectories of the ball.  The
%  main objective of these analyses is to display this variation in
%  timing, or what we call phase variation.  This achieved by the 
%  registration process, involving a smooth nonlinear strictly increasing
%  transformation of time called time-warping.  The function used here
%  for this purpose is the register_fd function in the functional data
%  analysis function library. 

%% Loading the required fda functions

addpath ('../fdaM')

%%  Load the data from ascii files Xcoord.txt, Ycoord.txt and Zcoord.txt
%  There are data from 10 trials, and the durations, number of
%  observations, and the number of cycles vary over trials.
%  Load the X, Y and Z coordinates for each trial, each matrix is
%  1876 by 10, 1876 being the maximum number of observations in any
%  trial.  The records begin and end at points corresponding to
%  minima of the tangential acceleration vector, as discussed in
%  Ramsay & Gribble, 1999. 

load Xcoord.txt
load Ycoord.txt
load Zcoord.txt

coord = zeros(size(Xcoord,1),3,10);  %  Set up array of coordinate values
coord(:,1,:) = Xcoord;
coord(:,2,:) = Ycoord;
coord(:,3,:) = Zcoord;

%%  Load number of cycles and number of observation points per record
%  1 record has 11 cycles, 5 have 12 cycles, and 4 have 13 cycles

load seqn.txt     %  the number of cycles per trial
disp('Number of cycles per trial:')
disp(num2str(seqn))

% Number of cycles per trial:
% 12  13  12  12  12  13  11  13  13  12

load ni.txt   % the number of observations in each trial
disp('Number of observations per trial:')
disp(num2str(ni'))

% disp('Number of observations per trial:')
% 1743  1876  1718  1712  1697  1824  1592  1834  1833  1712

%% Observation times in seconds
%  The data were recorded at 200 herz, or every five milliseconds.
%  Set the time to be in milliseconds, this time being 
%  about the finest inter-knot interval we'd be needing to use.

timei = 5*(ni-1);  %  duration of entire record in milliseconds

timepercycle = timei./seqn';  %  mean duration per cycle
disp('Mean cycle duration for each record in milliseconds:')
disp(num2str(round(timepercycle)'))

% Mean cycle duration for each record in milliseconds:
% 726  721  715  713  707  701  723  705  705  713

meantimepercycle = mean(timepercycle);  % mean duration over all records
disp(['Mean cycle duration over all records: ', ...
      num2str(round(meantimepercycle)), ' milliseconds'])

% Mean cycle duration over all records: 713 milliseconds

%% Record number, number of cycles, number of observations, 
%     time per cycle, and duration per cycle (msec):
%            1          12        1743        8710         726
%            2          13        1876        9375         721
%            3          12        1718        8585         715
%            4          12        1712        8555         713
%            5          12        1697        8480         707
%            6          13        1824        9115         701
%            7          11        1592        7955         723
%            8          13        1834        9165         705
%            9          13        1833        9160         705
%           10          12        1712        8555         713

%%  Plot each coordinate separately for each trial

index = 1;

figure(1)
for i=index
    subplot(3,1,1)
    tveci = linspace(0,timei(i),ni(i))';
    plot(tveci, coord(1:ni(i),1,i), '.-',[0,9500],[0,0],'r:')
    xlabel('')
    ylabel('\fontsize{13} X')
    axis([0,9500,-.1,.1])
    subplot(3,1,2)
    plot(tveci, coord(1:ni(i),2,i), '.-',[0,9500],[0,0],'r:')
    xlabel('')
    ylabel('\fontsize{13} Y')
    axis([0,9500,-.1,.1])
    subplot(3,1,3)
    plot(tveci, coord(1:ni(i),3,i), '.-',[0,9500],[0,0],'r:')
    xlabel('\fontsize{13} millisec.')
    ylabel('\fontsize{13} Z')
    axis([0,9500,-.2,.2])
    title(['\fontsize{13} trial ',num2str(i)])
    if length(index) > 1
        pause
    end
end

%%  Load timings of features for tangential velocity and acceleration.
%  By "tangential" is meant the velocity or acceleration along the
%  trajectory of the forefinger.  This is the length of the velocity or
%  acceleration vector.
%  
%  The tangential velocity features for each cycle as follows:
%  1. a minor local maximum in tangential velocity positioned between 
%     the minimum and the primary maximum
%  2. a major local maximum in tangential velocity, and
%  3. a primary local minimum in tangential velocity
%  The tangential acceleration features for each cycle are as follows:
%  1. a major acceleration peak 
%  2. a minor peak
%  3. a minimum
%  4. a small maximum, but small only because the ball is also being
%     accelerated
%  5. a minimum

load Vfeaturetimes  % tangential velocity     feature times
load Afeaturetimes  % tangential acceleration feature times

%% Display feature times for velocity

% disp(['    launch', '    drop', '      handoff', ])
% Vfeaturemean(1:3,:)'
%
%     launch    drop      handoff
% 
%     0.0471    0.2662    0.4939
%     0.7702    0.9760    1.2151
%     1.4939    1.6901    1.9294
%     2.1880    2.3997    2.6406
%     2.9008    3.1169    3.3566
%     3.6366    3.8396    4.1350
%     4.3422    4.5422    4.7089
%     5.0537    5.2578    5.5072
%     5.7699    5.9808    6.2116
%     6.4810    6.6859    6.9274
%     7.2035    7.4097    7.6314
%     7.9129    8.1116    8.3482
%     8.5845    8.7817    9.0020

%% Display feature times for acceleration

% disp(['centripital', '  brake', '  handoff', '  launch'])
% Afeaturestd(1:4,:)'
% 
% centripital   brake     handoff   launch
% 
%     0.0278    0.0329    0.0335    0.0394
%     0.0315    0.0363    0.0326    0.0477
%     0.0380    0.0375    0.0452    0.0511
%     0.0382    0.0414    0.0440    0.0516
%     0.0428    0.0436    0.0497    0.0551
%     0.0464    0.0486    0.0416    0.0503
%     0.0508    0.0537    0.0833    0.0920
%     0.0624    0.0680    0.0736    0.0733
%     0.0707    0.0797    0.0826    0.0828
%     0.0773    0.0826    0.0761    0.0854
%     0.0877    0.0932    0.0856    0.1050
%     0.0900    0.0927    0.0932    0.1028 
%     0.1212    0.1230    0.1223    0.1151

%%  Cell array ContuousRegCell stores smooths of the the data by 
%  unpenalized least squares using a minimal basis.
%  The part of each record before the first velocity feature time for
%  the first cycle is removed.

%  The basis for this smooth is constructed by using nine knots
%  per cycle, and the smoothing is by unpenalized least squares.
%  The derivative estimates are unstable at the end points, but this
%  basis is much more suitable for registering the data
%  Because the records vary in length, this basis must be set up
%  separately for each record

load recordfdCell0

%%  Landmark registration of each curve using first velocity peak.
%  The first velocity peak is at or just before the throw, and used
%  to define the beginning and end of each cycle.  This registration
%  approximately normalizes the cycle durations.

%  load cellarray containing results

load landmarkregCell

%  Plot the unregistered and registered records along with the 
%  time warpings for each

index = 1:10;
for i=index
    tveci    = recordfdCell0{i,1};
    fdobji   = recordfdCell0{i,3};    %  unregistered function
    regfdi   = landmarkregCell{i,1};  %  registered function
    warpfdi  = landmarkregCell{i,2};  %  warping function
    Wfdi     = landmarkregCell{i,3};  %  W(t) = log-deriv of h(t)
    %  plot the unregistered and registered functions
    figure(1)
    ymati    = eval_fd(tveci, fdobji); 
    yregmati = eval_fd(tveci, regfdi);
    %  loop through coordinates
    for j=1:3
        subplot(3,1,j)
        plot(tveci, ymati(:,1,j), 'r-', ...
             tveci, yregmati(:,1,j), 'b-')
        if j==1
            title(['Record ',num2str(i)])
        end
    end
    %  plot the deformation functions, d(t) = h(t) - t
    warpmat = eval_fd(tveci, warpfdi);
    rngi    = [0,max(recordfdCell0{i,1})];
    figure(2)
    subplot(1,1,1)
    plot(tveci, warpmat-tveci, 'b-', ...
         [rngi(1), rngi(2)], [0,0], 'r:')
    %  pause if there are multiple records
    if length(index) > 1
        pause
    end
end

  
%%  A strictly periodic smooth of each record.
%  This part of this demo file is intended to show how we can use
%  registration to register each record to a strictly periodic image
%  of itself, where the period is set to the mean period across all
%  cycles, 713 milliseconds.

%  Re-smooth data in two ways:  first with the B-spline basis with
%  ten knots per cycle, and then with a fourier basis with a period
%  equal to the mean cycle duration over all cycles, 712.9 milliseconds.

ContinuousRegCell = cell(10,8);

nbasis0 = 3;  %  number of basis functions per cycle for warping

norder  = 6;
nknots0 = 9;
period  = meantimepercycle;
for i=1:10
    disp(['Record ',num2str(i)])
    %  set up interval for the cycle as the period times the
    %  number of cycles
    Ti       = period*seqn(i);  %  end of interval
    rngi     = [0,Ti];          %  range for registered functions
    startime = Vfeaturetimes(1,1,i);  % time of first velocity peak
    %  select observation times starting with startime
    tvec     = linspace(0,timei(i)/1000,ni(i))';
    indi     = find(tvec >= startime);
    tvec     = tvec(indi) - startime;
    ntvec    = length(tvec);
    tvec     = linspace(0,Ti,ntvec)';
    %  rescale selected times to span range
    sclfac   = Ti/(tvec(ntvec)-tvec(1));
    tvec     = tvec*sclfac;
    %  splinebasis and fdPar objects for this interval
    Nnbasis   = seqn(i)*nknots0 + 5;
    Nbasisobj = create_bspline_basis(rngi, Nnbasis, norder);
    NfdParobj = fdPar(Nbasisobj);
    %  non-periodic smooth of the data 
    data     = squeeze(coord(indi,:,i));
    Nfdobj   = smooth_basis(tvec,data,NfdParobj);
    %  coerce the coefficient matrix to array format
    %    so as to plot as multivariate functional data
    coefmat  = getcoef(Nfdobj);
    coefdim  = size(coefmat);
    Ncoefarry = reshape(coefmat,coefdim(1),1,coefdim(2));
    Nfdobj  = putcoef(Nfdobj, Ncoefarry);
    Nfdobj  = putnames(Nfdobj, jugglenames);
    %  set up a fourier basis with period = average cycle length
    Pnbasis = 10;
    Pbasisobj = create_fourier_basis([0,Ti], Pnbasis, period);
    PfdParobj = fdPar(Pbasisobj);
    %  periodic smooth of the data
    Pfdobj   = smooth_basis(tvec,data,PfdParobj);
    %  coerce the coefficient matrix to array format
    coefmat  = getcoef(Pfdobj);
    coefdim  = size(coefmat);
    Pcoefarry = reshape(coefmat,coefdim(1),1,coefdim(2));
    Pfdobj  = putcoef(Pfdobj, Pcoefarry);
    Pfdobj  = putnames(Pfdobj, jugglenames);
    %  set up a basis for the function W defining warping function h
    wbasisi = create_bspline_basis(rngi, seqn(i)*nbasis0+3);
    WfdPari = fdPar(wbasisi,2,lambda);
    %  continuously register the objects, with the periodic basis 
    %  smooth being the target
    [Rfdobj, warpfd, Wfd] = ...
              register_fd(Pfdobj, Nfdobj, WfdPari);
    %  record the results
    ContinuousRegCell{i,1} = tvec;
    ContinuousRegCell{i,2} = data;
    ContinuousRegCell{i,3} = Nfdobj;
    ContinuousRegCell{i,4} = Pfdobj;
    ContinuousRegCell{i,5} = Rfdobj;
    ContinuousRegCell{i,6} = warpfd;
    ContinuousRegCell{i,7} = Wfd;
    ContinuousRegCell{i,8} = sclfac;
end

%  plot the results

index = 1:10;
for i=index
    Ti       = period*seqn(i);  %  end of interval
    tveci   = ContinuousRegCell{i,1};    
    Nfdobji = ContinuousRegCell{i,3};  %  unregistered function
    Rfdobji = ContinuousRegCell{i,5};  %  registered function
    warpfdi = ContinuousRegCell{i,6};  %  warping function
    Wfdi    = ContinuousRegCell{i,7};  %  W(t) = log-deriv of h(t)
    sclfac  = ContinuousRegCell{i,8};  %  rescaling factor
    %  plot the unregistered and registered functions
    ymati    = eval_fd(tveci, Nfdobji); 
    yregmati = eval_fd(tveci, Rfdobji);
    figure(1)
    for j=1:3
        subplot(3,1,j)
        plot(tveci, ymati(:,1,j), 'r-', ...
             tveci, yregmati(:,1,j), 'b-', ...
             [0, Ti], [0,0], 'r:')
        if j==3
            xlabel('\fontsize{13} time (msec)')
        end
        ylabel('\fontsize{13} Meters')
        title(['\fontsize{13} Coordinate ', num2str(j)])
        axis([0,Ti,-0.2,0.2])
        if j==1
            legend('non-periodic', 'periodic')
        end
    end
    %  plot the deformation functions, d(t) = h(t) - t
    deformi = eval_fd(tveci, warpfdi) - tveci;
    figure(2)
    subplot(1,1,1)
    plot(tveci, deformi, 'b-', ...
        [0, Ti], [0,0], 'r:')
    hold on
    %  plot velocity features
    if ~isempty(Vfeaturetimes)
        Vnfeature = size(Vfeaturetimes,1);
        for j=1:Vnfeature
            Vtime = sclfac*Vfeaturetimes(j,1:seqn(i),i)*1000;
            Veval = interp1(tveci, deformi, Vtime);
            plot(Vtime, Veval, [colorvec(j),'o'])
        end
    end
    for k=1:seqn(i)
        Vtime1k = sclfac*Vfeaturetimes(1,k,i)*1000;
        plot([Vtime1k,Vtime1k],[-90,90],'b--')
    end
    if ~isempty(Afeaturetimes)
        %  plot acceleration features
        if ~isempty(Afeaturetimes)
            Anfeature = size(Afeaturetimes,1);
            for j=1:Anfeature
                Atime = sclfac*Afeaturetimes(j,1:seqn(i),i)*1000;
                Aeval = interp1(tveci, deformi, Atime);
                plot(Atime, Aeval, [colorvec(j),'*'])
            end
        end
    end
    hold off
    axis([0,Ti,-90,90])
    xlabel('\fontsize{13} time (msec)')
    ylabel('\fontsize{13} Deformation d(t) = h(t) - t (msec)')
    title(['\fontsize{13} Record ', num2str(i), ...
           '  Number of cycles = ', num2str(seqn(i))])
    %  pause if there are multiple records
    if length(index) > 1
        pause
    end
end

