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

//Info << "smoothing at finest level." << endl;

smoothers[0].smooth
(
    psi,
    source,
    cmpt,
    nFinestSweeps_
);
