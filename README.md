# Notes #

## Git Stuff ##
1. To see which m-files and simulink files are being tracked by git on branch `master` issue
    * `git ls-tree -r master | grep -E "(*\.m|*.slx)" `

## Comments ##
1. `minseg_project.m`
    I tried to keep everything in one file, `minseg_project.m`, but there are two different deliverables - 
        * The report
        * The controller
    So, it seems easier to have the report generator be `minseg_project.m`
        * we don't need to use everything, but it can generate pieces we want to insert into the report
    And the controller can be run by, `run_MinSeg_Controller.m`
2. `run_MinSeg_Controller.m`
    * This can load parameters and state-space matrices into the workspace before calling the simulink model.
    * The only caveat would be to remember, before the final report to place the final parameters, etc. back into the report
3. `MinSeg_Controller.slx`
    * I'm not sure if this is correct... are the other state variables not needed to balance? 


## Good References ##
[Inverted Pendulum: State-Space Methods for Controller Design](http://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=ControlStateSpace)

[Inverted Pendulum: System Modeling](http://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=SystemModeling)

[Linear-Quadratic-Gaussian Control](http://www.mathworks.com/help/control/multi-input-multi-output-control-design.html)

[Design an LQG Regulator](http://www.mathworks.com/help/control/getstart/functions-for-compensator-design.html#f2-1031766)

[State Feedback Control in Simulink](http://ctms.engin.umich.edu/CTMS/index.php?example=AircraftPitch&section=SimulinkModeling)

[Digital State Space Controller Design](http://ctms.engin.umich.edu/CTMS/index.php?example=AircraftPitch&section=ControlDigital#2)

## Emails ##
1. 
    From: Schulenberg, Christopher S 
    Sent: Thursday, March 05, 2015 9:09 PM
    To: Adams, Paul R
    Subject: Running the minseg untethered

    I was able to get MATLAB 2014b installed with the Arduino support package and run the model un-tethered.

    1.  Use the ‘deploy to target hardware’ button to deploy the model to Arduino:
     
    2.  Do not press the ‘play’ button
    3.  Untether the minseg
    4.  Power on the minseg

    I think when you press the play button it loads the ‘external mode’ code which relies on triggers from PC Simulink to run. If you deploy to target HW it is ready to go at the next power cycle.

    The code with no modifications gives me ~2 seconds before scooting off to oblivion. Tomorrow’s task for me is to see how modify things to 1) reflect the analytical portion of the project in order to 2) actually stabilize the minseg, and 3) clean up the god awful project directory that was unzipped from the download.

    I’d also like to include a few features like cutting the motor when the thing is laying on its side.

    Thanks,
    Chris

2. 
    At this stage, we run the risk of code conflicts, you are welcome to modify any of this of course. I’ve taken my stab at it. 

    Simulink working model
    1.  I copied, en masse, the model DualMinSegBalance_Code_V31.slx to MinSeg_Controller.slx 
    a.  This gives us the following
    i.  Angle > theta kill switch
    ii. Driver conversion factors
    iii.    Gyro calibration template
    iv. Other tidbits
    b.  We need to modify/add, at the least
    i.  An estimator loop
    ii. Gain matrix K, matched to our minseg/model/state-matrices (LQR probably)
    iii.    Possibly other gotchas
    Parameter file
    2.  I updated/renamed the file MinSeg_Parameters.m to be loaded by MinSeg_Controller.slx and modified the parameters for our minseg/model/state-matrices. This file also contains the conversion factors used for DualMinSegBalance_Code_V31.slx and contains the code to compute K (currently implemented with LQR). 

    Report file
    3.  I renamed minseg_project.m to MinSeg_Report.m This is the Matlab which follows the pdf document step-by-step, generates plots, figures, text in a published format. I recommend we use this medium for the report itself as opposed to Word. If you want you could do the text write-up in Word and I can embed it in Matlab markup (I’m not implying you do all the write-up yourself ;)) I feel comfortable quickly adding/embedding anything into that file, then publishing it produces a report-ready document. We can discuss Saturday if you have strong feelings either way. 

    Code Organization/Git tips
    4.  You’ll need to pull everything down. There might be merge conflicts since you were working last night. You could create a branch called “chris�? to reserve your changes that you want to manually integrate. 
    5.  I’ve adopted the convention of naming all of our primary source files with a MinSeg pre-appendage to make them easily recognizable. 
    6.  Folders
    a.  The Report-related files are in the root (minseg) directory (and in the html directory, but we can ignore that folder). 
    b.  The MinSeg related files are in the simulink  directory
    7.  There are many files (c source, Simulink-generated, etc.) that we don’t want to be passing around via git. 
    a.  To see all the files currently managed by git, issue
    i.  git ls-tree -r master
    b.  To see all the files managed by git that are source code, issue
    i.  git ls-tree -r master | grep -E "(*\.m|*\.slx)"
    c.  If there are other files, tell me and I’ll add them, or help you add them. The .gitignore is set up such that everything is ignored, except what we explicitly declare not to ignore. It’s more of a notgitignore…
    MinSeg balancing thoughts
    8.  I ran the same experiment as you last night, got it running and began trying to figure out how to make it balance. One key is finding the right k. We can wait til tomorrow to discuss in depth. 
    9.  I got it to the point where it tried to balance, but was responding wayyy to slowly. As mentioned above, we need to implement an observer. 
    10. Good tip about ‘deploy to hardware’ vs. ‘run in external mode’


