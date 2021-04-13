void Foam::GAMGSolverNew::VcycleNew
(
    const PtrList<lduMatrix::smoother>& smoothers,

    scalarField& scratch1,
    scalarField& scratch2,

    PtrList<scalarField>& coarseCorrFields,
    PtrList<scalarField>& coarseSources,
    const direction cmpt,
    label leveli
) const
{
    static PtrList<scalarField> coarseEpsFields(coarseCorrFields);

    const label coarsestLevel = matrixLevels_.size() - 1;

    // Residual restriction (going to coarser levels)
    if (coarseSources.set(leveli + 1))
    {
        //Info << "current level in former stage: " << leveli << endl;

        // If the optional pre-smoothing sweeps are selected
        // smooth the coarse-grid field for the restricted source
        if (nPreSweeps_)
        {
            coarseCorrFields[leveli] = 0.0;

            smoothers[leveli + 1].smooth
            (
                coarseCorrFields[leveli],
                coarseSources[leveli],
                cmpt,
                min
                (
                    nPreSweeps_ +  preSweepsLevelMultiplier_*leveli,
                    maxPreSweeps_
                )
            );

            scalarField::subField ACf
            (
                scratch1,
                coarseCorrFields[leveli].size()
            );

            // Scale coarse-grid correction field
            // but not on the coarsest level because it evaluates to 1
            if (scaleCorrection_ && leveli < coarsestLevel - 1)
            {
                scale
                (
                    coarseCorrFields[leveli],
                    const_cast<scalarField&>
                    (
                        ACf.operator const scalarField&()
                    ),
                    matrixLevels_[leveli],
                    interfaceLevelsBouCoeffs_[leveli],
                    interfaceLevels_[leveli],
                    coarseSources[leveli],
                    cmpt
                );
            }

            // Correct the residual with the new solution
            matrixLevels_[leveli].Amul
            (
                const_cast<scalarField&>
                (
                    ACf.operator const scalarField&()
                ),
                coarseCorrFields[leveli],
                interfaceLevelsBouCoeffs_[leveli],
                interfaceLevels_[leveli],
                cmpt
            );

            coarseSources[leveli] -= ACf;
        }

        // Info << "restrictionV in level: " << leveli << endl;

        // Residual is equal to source
        agglomeration_.restrictField
        (
            coarseSources[leveli + 1],
            coarseSources[leveli],
            leveli + 1,
            true
        );

        coarseEpsFields[leveli + 1] = 0.0;
    }




    if (debug >= 2 && nPreSweeps_)
    {
        Pout<< endl;
    }

    if (leveli == matrixLevels_.size() - 2) {
        if (coarseCorrFields.set(coarsestLevel))
        {
            //Info << "do smoothing in " << leveli + 1 << "  " << coarsestLevel << endl;
            solveCoarsestLevel
            (
                coarseEpsFields[coarsestLevel],
                coarseSources[coarsestLevel]
            );
        }       
    }
    else {
        Foam::GAMGSolverNew::VcycleNew
        (
            smoothers,
            
            scratch1,
            scratch2,

            coarseEpsFields,
            coarseSources,
            cmpt,
            leveli + 1
        );
    }
    
    scalarField dummyField(0);

    if (coarseCorrFields.set(leveli))
    {
        //Info << "current level in latter stage: " << leveli << endl;

        // Create a field for the pre-smoothed correction field
        // as a sub-field of the finestCorrection which is not
        // currently being used
        scalarField::subField preSmoothedCoarseCorrField
        (
            scratch2,
            coarseCorrFields[leveli].size()
        );

        // Only store the preSmoothedCoarseCorrField if pre-smoothing is
        // used
        if (nPreSweeps_)
        {
            preSmoothedCoarseCorrField = coarseCorrFields[leveli];
        }

        coarseEpsFields[leveli + 1] *= relaxationFactor_; 

        // Info << "prolongV in level: " << leveli << endl;

        agglomeration_.prolongField
        (
            coarseCorrFields[leveli],
            (
                coarseEpsFields.set(leveli + 1)
                ? coarseEpsFields[leveli + 1]
                : dummyField              // dummy value
            ),
            leveli + 1,
            true
        );


        // Create A.psi for this coarse level as a sub-field of Apsi
        scalarField::subField ACf
        (
            scratch1,
            coarseCorrFields[leveli].size()
        );
        scalarField& ACfRef =
            const_cast<scalarField&>(ACf.operator const scalarField&());

        if (interpolateCorrection_) //&& leveli < coarsestLevel - 2)
        {
            if (coarseCorrFields.set(leveli+1))
            {
                interpolate
                (
                    coarseCorrFields[leveli],
                    ACfRef,
                    matrixLevels_[leveli],
                    interfaceLevelsBouCoeffs_[leveli],
                    interfaceLevels_[leveli],
                    agglomeration_.restrictAddressing(leveli + 1),
                    coarseCorrFields[leveli + 1],
                    cmpt
                );
            }
            else
            {
                interpolate
                (
                    coarseCorrFields[leveli],
                    ACfRef,
                    matrixLevels_[leveli],
                    interfaceLevelsBouCoeffs_[leveli],
                    interfaceLevels_[leveli],
                    cmpt
                );
            }
        }

        // Scale coarse-grid correction field
        // but not on the coarsest level because it evaluates to 1
        if
        (
            scaleCorrection_
            && (interpolateCorrection_ || leveli < coarsestLevel - 1)
        )
        {
            scale
            (
                coarseCorrFields[leveli],
                ACfRef,
                matrixLevels_[leveli],
                interfaceLevelsBouCoeffs_[leveli],
                interfaceLevels_[leveli],
                coarseSources[leveli],
                cmpt
            );
        }

        // Only add the preSmoothedCoarseCorrField if pre-smoothing is
        // used
        if (nPreSweeps_)
        {
            coarseCorrFields[leveli] += preSmoothedCoarseCorrField;
        }
        //Info << "number of smoothing: " <<min(nPostSweeps_ + postSweepsLevelMultiplier_*leveli, maxPostSweeps_)<< "in level " << leveli << endl;
        smoothers[leveli + 1].smooth
        (
            coarseCorrFields[leveli],
            coarseSources[leveli],
            cmpt,
            min
            (
                nPostSweeps_ + postSweepsLevelMultiplier_*leveli,
                maxPostSweeps_
            )
        );
    }
}