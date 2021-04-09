# GAMGnew
## expanding the feather of GAMG solver in OpenFOAM with more solving strategies
add F and V cycle according to https://en.wikipedia.org/wiki/Multigrid_method

how to use:
1. go into root directory and use the command: wmake libso (OpenFOAM v8 environment needed)

2. Add library "libGAMGnew.so" in ControlDict of wanted case as follow:

      libs 
      (
            "libOpenFOAM.so"
      )
      
3. change keyword solver in i.e. for pressure solver in fvSolution as follow:

    p
    {
        solver          GAMGnew;
        cycleType       Wcycle;//Vcycle, VcycleNew
        tolerance       1e-08;
        relTol          0.1;
        smoother        GaussSeidel;
    }

4. run the case as normal
