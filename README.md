# GAMGnew
## expanding the feather of GAMG solver in OpenFOAM with more solving strategies
add F ,W and V cycle according to https://en.wikipedia.org/wiki/Multigrid_method

how to use:
1. go into root directory and use the command: wmake libso (OpenFOAM v8 environment needed)

2. Add library "libGAMGnew.so" in ControlDict of wanted case as follow:

      libs 
      (
            "libOpenFOAM.so"
      )
      
3. change keyword solver in i.e. for pressure solver in fvSolution as follow, where cycle type is specified by keyword "cycleType" and relaxation factor for prolongation by keyword "relaxationFactor":

    p
    {
        solver          GAMGnew;//GAMG1;
        cycleType       VcycleNew;//FcycleNoRecursion;//Vcycle;//VcycleNew;//Wcycle;//FcycleNoRecursion;//FcycleNoRecursion;
        tolerance       1e-08;
        relTol          0.1;
        smoother        GaussSeidel;
        relaxationFactor  0.5;
        //minIter         0;
        //nPostSweeps     4;
        //maxPostSweeps   8;
        //FinestSweeps    1;
    }



4. run the case as normal
