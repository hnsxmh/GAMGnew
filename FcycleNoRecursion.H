void Foam::GAMGSolverNew::FcycleNoRecursion
(
    const PtrList<lduMatrix::smoother>& smoothers,
    scalarField& psi,
    const scalarField& source,
    scalarField& Apsi,
    scalarField& finestCorrection,
    scalarField& finestResidual,

    scalarField& scratch1,
    scalarField& scratch2,

    PtrList<scalarField>& coarseCorrFields,
    PtrList<scalarField>& coarseSources,
    const direction cmpt
) const
{
    //Info << "W cycle begins!" << endl;
    // debug = 2;

    const label coarsestLevel = matrixLevels_.size() - 1;

    // Restrict finest grid residual for the next level up.
    agglomeration_.restrictField(coarseSources[0], finestResidual, 0, true); 

    coarseCorrFields[0] = 0.0;

    if (debug >= 2 && nPreSweeps_)
    {
        Pout<< "Pre-smoothing scaling factors: ";
    }

    // Residual restriction (going to coarser levels)
    for (label leveli = 0; leveli < coarsestLevel; leveli++)
    {
        if (coarseSources.set(leveli + 1))
        {
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

            // Info << "residual in  "<< leveli <<" : " << gSumMag(coarseSources[leveli], matrix().mesh().comm()) << endl;

            // Residual is equal to source
            agglomeration_.restrictField
            (
                coarseSources[leveli + 1],
                coarseSources[leveli],
                leveli + 1,
                true
            );

            coarseCorrFields[leveli + 1] = 0.0;
        }
    }

    if (debug >= 2 && nPreSweeps_)
    {
        Pout<< endl;
    }

    // Solve Coarsest level with either an iterative or direct solver
    
    if (coarseCorrFields.set(coarsestLevel))
    {
        solveCoarsestLevel
        (
            coarseCorrFields[coarsestLevel],//coarseCorrFields[coarsestLevel],
            coarseSources[coarsestLevel]
        );
    }

    if (debug >= 2)
    {
        Pout<< "Post-smoothing scaling factors: ";
    }

    // Smoothing and prolongation of the coarse correction fields
    // (going to finer levels)

    scalarField dummyField(0);

    for (label curMaxlevel = coarsestLevel - 1; curMaxlevel >= 0; curMaxlevel--) 
    {
        // for (label leveli = coarsestLevel - 1; leveli >= curMaxlevel; leveli--)
        // {
        //     if (coarseCorrFields.set(leveli))
        //     {
        //         //Info << "current level for prolognation: " << leveli << endl;
        //         // Create a field for the pre-smoothed correction field
        //         // as a sub-field of the finestCorrection which is not
        //         // currently being used
        //         scalarField::subField preSmoothedCoarseCorrField
        //         (
        //             scratch2,
        //             coarseCorrFields[leveli].size()
        //         );

        //         // Only store the preSmoothedCoarseCorrField if pre-smoothing is
        //         // used
        //         if (nPreSweeps_)
        //         {
        //             preSmoothedCoarseCorrField = coarseCorrFields[leveli];
        //         }

        //         //Info << "correction in level " << leveli + 1 << " is: " << gSumMag(coarseCorrFields[leveli + 1], matrix().mesh().comm()) << endl;

        //         agglomeration_.prolongField
        //         (
        //             coarseCorrFields[leveli],
        //             (
        //                 coarseCorrFields.set(leveli + 1)
        //             ? coarseCorrFields[leveli + 1]//coarseCorrFields[leveli + 1]
        //             : dummyField              // dummy value
        //             ),
        //             leveli + 1,
        //             true
        //         );

        //         // Create A.psi for this coarse level as a sub-field of Apsi
        //         scalarField::subField ACf
        //         (
        //             scratch1,
        //             coarseCorrFields[leveli].size()
        //         );
        //         scalarField& ACfRef =
        //             const_cast<scalarField&>(ACf.operator const scalarField&());

        //         if (interpolateCorrection_) //&& leveli < coarsestLevel - 2)
        //         {
        //             if (coarseCorrFields.set(leveli+1))
        //             {
        //                 interpolate
        //                 (
        //                     coarseCorrFields[leveli],
        //                     ACfRef,
        //                     matrixLevels_[leveli],
        //                     interfaceLevelsBouCoeffs_[leveli],
        //                     interfaceLevels_[leveli],
        //                     agglomeration_.restrictAddressing(leveli + 1),
        //                     coarseCorrFields[leveli + 1],
        //                     cmpt
        //                 );
        //             }
        //             else
        //             {
        //                 interpolate
        //                 (
        //                     coarseCorrFields[leveli],
        //                     ACfRef,
        //                     matrixLevels_[leveli],
        //                     interfaceLevelsBouCoeffs_[leveli],
        //                     interfaceLevels_[leveli],
        //                     cmpt
        //                 );
        //             }
        //         }

        //         // Scale coarse-grid correction field
        //         // but not on the coarsest level because it evaluates to 1
        //         if
        //         (
        //             scaleCorrection_
        //         && (interpolateCorrection_ || leveli < coarsestLevel - 1)
        //         )
        //         {
        //             scale
        //             (
        //                 coarseCorrFields[leveli],
        //                 ACfRef,
        //                 matrixLevels_[leveli],
        //                 interfaceLevelsBouCoeffs_[leveli],
        //                 interfaceLevels_[leveli],
        //                 coarseSources[leveli],
        //                 cmpt
        //             );
        //         }

        //         // Only add the preSmoothedCoarseCorrField if pre-smoothing is
        //         // used
        //         if (nPreSweeps_)
        //         {
        //             coarseCorrFields[leveli] += preSmoothedCoarseCorrField;
        //         }

        //         //scalar normFactor = this->normFactor(psi, source, Apsi, coarseSources[leveli]);

        //         //Info << "residual before smoothing: "<< gSumMag(coarseSources[leveli], matrix().mesh().comm())/normFactor; << endl;
        //         //Info << "residual before smoothing: "<< gSumMag(coarseSources[leveli], matrix().mesh().comm()) << endl;

        //         smoothers[leveli + 1].smooth
        //         (
        //             coarseCorrFields[leveli],
        //             coarseSources[leveli],
        //             cmpt,
        //             min
        //             (
        //                 nPostSweeps_ + postSweepsLevelMultiplier_*leveli,
        //                 maxPostSweeps_
        //             )
        //         );
        //         //Info << "correction in  "<< leveli <<" : " << gSumMag(coarseSources[leveli], matrix().mesh().comm()) << endl;
        //         /////////////////////////////add source correct
        //         //Info <<"doing smooth"<<min(nPostSweeps_ + postSweepsLevelMultiplier_*leveli,maxPostSweeps_)<< "times in current level: " << leveli + 1 << endl;
        //         if (leveli == curMaxlevel) 
        //         {
        //             matrixLevels_[leveli].Amul
        //             (
        //                 const_cast<scalarField&>
        //                 (
        //                     ACf.operator const scalarField&()
        //                 ),
        //                 coarseCorrFields[leveli],
        //                 interfaceLevelsBouCoeffs_[leveli],
        //                 interfaceLevels_[leveli],
        //                 cmpt
        //             );

        //             coarseSources[leveli] -= ACf;

        //             // for (label index = curMaxlevel + 1; index <= coarsestLevel; index++) 
        //             // {   
        //             //     // Info << "index: " << index << endl;
        //             //     coarseCorrFields[index] = 0.0;
        //             // }
        //         }
        //         else 
        //         {
        //             //Info << "index: " << leveli << endl;
        //             coarseCorrFields[leveli] = 0.0;
        //         }

        //         // coarseCorrFields[leveli] = 0.0;

        //         //Info << "residual after smoothing: "<< gSumMag(coarseSources[leveli], matrix().mesh().comm())/normFactor; << endl;
        //         //Info << "residual after smoothing: "<< gSumMag(coarseSources[leveli], matrix().mesh().comm()) << endl;
        //     }
        // }


        if (coarseCorrFields.set(curMaxlevel))
        {
            //Info << "current level for prolognation: " << leveli << endl;
            // Create a field for the pre-smoothed correction field
            // as a sub-field of the finestCorrection which is not
            // currently being used
            scalarField::subField preSmoothedCoarseCorrField
            (
                scratch2,
                coarseCorrFields[curMaxlevel].size()
            );

            // Only store the preSmoothedCoarseCorrField if pre-smoothing is
            // used
            if (nPreSweeps_)
            {
                preSmoothedCoarseCorrField = coarseCorrFields[curMaxlevel];
            }

            //Info << "correction in level " << leveli + 1 << " is: " << gSumMag(coarseCorrFields[leveli + 1], matrix().mesh().comm()) << endl;

            agglomeration_.prolongField
            (
                coarseCorrFields[curMaxlevel],
                (
                    coarseCorrFields.set(curMaxlevel + 1)
                ? coarseCorrFields[curMaxlevel + 1]//coarseCorrFields[leveli + 1]
                : dummyField              // dummy value
                ),
                curMaxlevel + 1,
                true
            );

            // Create A.psi for this coarse level as a sub-field of Apsi
            scalarField::subField ACf
            (
                scratch1,
                coarseCorrFields[curMaxlevel].size()
            );
            scalarField& ACfRef =
                const_cast<scalarField&>(ACf.operator const scalarField&());

            if (interpolateCorrection_) //&& leveli < coarsestLevel - 2)
            {
                if (coarseCorrFields.set(curMaxlevel+1))
                {
                    interpolate
                    (
                        coarseCorrFields[curMaxlevel],
                        ACfRef,
                        matrixLevels_[curMaxlevel],
                        interfaceLevelsBouCoeffs_[curMaxlevel],
                        interfaceLevels_[curMaxlevel],
                        agglomeration_.restrictAddressing(curMaxlevel + 1),
                        coarseCorrFields[curMaxlevel + 1],
                        cmpt
                    );
                }
                else
                {
                    interpolate
                    (
                        coarseCorrFields[curMaxlevel],
                        ACfRef,
                        matrixLevels_[curMaxlevel],
                        interfaceLevelsBouCoeffs_[curMaxlevel],
                        interfaceLevels_[curMaxlevel],
                        cmpt
                    );
                }
            }

            // Scale coarse-grid correction field
            // but not on the coarsest level because it evaluates to 1
            if
            (
                scaleCorrection_
            && (interpolateCorrection_ || curMaxlevel < coarsestLevel - 1)
            )
            {
                scale
                (
                    coarseCorrFields[curMaxlevel],
                    ACfRef,
                    matrixLevels_[curMaxlevel],
                    interfaceLevelsBouCoeffs_[curMaxlevel],
                    interfaceLevels_[curMaxlevel],
                    coarseSources[curMaxlevel],
                    cmpt
                );
            }

            // Only add the preSmoothedCoarseCorrField if pre-smoothing is
            // used
            if (nPreSweeps_)
            {
                coarseCorrFields[curMaxlevel] += preSmoothedCoarseCorrField;
            }

            //scalar normFactor = this->normFactor(psi, source, Apsi, coarseSources[leveli]);

            //Info << "residual before smoothing: "<< gSumMag(coarseSources[leveli], matrix().mesh().comm())/normFactor; << endl;
            //Info << "residual before smoothing: "<< gSumMag(coarseSources[leveli], matrix().mesh().comm()) << endl;

            smoothers[curMaxlevel + 1].smooth
            (
                coarseCorrFields[curMaxlevel],
                coarseSources[curMaxlevel],
                cmpt,
                min
                (
                    nPostSweeps_ + postSweepsLevelMultiplier_*curMaxlevel,
                    maxPostSweeps_
                )
            );
            //Info << "correction in  "<< leveli <<" : " << gSumMag(coarseSources[leveli], matrix().mesh().comm()) << endl;
            /////////////////////////////add source correct
            //Info <<"doing smooth"<<min(nPostSweeps_ + postSweepsLevelMultiplier_*leveli,maxPostSweeps_)<< "times in current level: " << leveli + 1 << endl;

            matrixLevels_[curMaxlevel].Amul
            (
                const_cast<scalarField&>
                (
                    ACf.operator const scalarField&()
                ),
                coarseCorrFields[curMaxlevel],
                interfaceLevelsBouCoeffs_[curMaxlevel],
                interfaceLevels_[curMaxlevel],
                cmpt
            );

            coarseSources[curMaxlevel] -= ACf;
        }





        for (label leveli = curMaxlevel; leveli < coarsestLevel; leveli++)
        {
            if (coarseSources.set(leveli + 1))
            {
                //Info << "current level for restriction: " << leveli << endl;
                // If the optional pre-smoothing sweeps are selected
                // smooth the coarse-grid field for the restricted source
                if (nPreSweeps_)
                {
                    // coarseCorrFields[leveli] = 0.0;

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

                // Residual is equal to source
                //Info << "correction in  "<< leveli <<" : " << gSumMag(coarseCorrFields[leveli], matrix().mesh().comm()) << endl;
                agglomeration_.restrictField
                (
                    coarseSources[leveli + 1],
                    coarseSources[leveli],
                    leveli + 1,
                    true
                );
            }
        }

        if (debug >= 2 && nPreSweeps_)
        {
            Pout<< endl;
        }

        // Info << "residual in  "<< coarsestLevel <<" : " << gSumMag(coarseSources[coarsestLevel], matrix().mesh().comm()) << endl;

        // Solve Coarsest level with either an iterative or direct solver
        if (coarseCorrFields.set(coarsestLevel))
        {
            //coarseCorrFields[coarsestLevel] = 0.0;

            solveCoarsestLevel
            (
                coarseCorrFields[coarsestLevel],
                coarseSources[coarsestLevel]
            );
        }

        if (debug >= 2)
        {
            Pout<< "Post-smoothing scaling factors: ";
        }      
    }


    for (label leveli = coarsestLevel - 1; leveli >= 0; leveli--)
    {
        //Info << "current level for prologation: " << leveli << endl;
        if (coarseCorrFields.set(leveli))
        {
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

            // Info << "correction in "<< leveli <<" : " << gSumMag(coarseSources[leveli], matrix().mesh().comm()) << endl;

            agglomeration_.prolongField
            (
                coarseCorrFields[leveli],
                (
                    coarseCorrFields.set(leveli + 1)
                  ? coarseCorrFields[leveli + 1]
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

    // Prolong the finest level correction
    agglomeration_.prolongField
    (
        finestCorrection,
        coarseCorrFields[0],
        0,
        true
    );

    if (interpolateCorrection_)
    {
        interpolate
        (
            finestCorrection,
            Apsi,
            matrix_,
            interfaceBouCoeffs_,
            interfaces_,
            agglomeration_.restrictAddressing(0),
            coarseCorrFields[0],
            cmpt
        );
    }

    if (scaleCorrection_)
    {
        // Scale the finest level correction
        scale
        (
            finestCorrection,
            Apsi,
            matrix_,
            interfaceBouCoeffs_,
            interfaces_,
            finestResidual,
            cmpt
        );
    }

    forAll(psi, i)
    {
        psi[i] += finestCorrection[i];
    }

    smoothers[0].smooth
    (
        psi,
        source,
        cmpt,
        nFinestSweeps_
    );

    // for (label leveli = 0; leveli < coarseCorrFields.size(); leveli++)
    // {
    //     coarseCorrFields[leveli] = 0.0;
    // }
}