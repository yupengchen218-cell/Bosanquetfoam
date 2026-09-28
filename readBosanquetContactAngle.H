/*---------------------------------------------------------------------------*\
  Read an area-weighted contact angle from one alpha.water patch.
\*---------------------------------------------------------------------------*/

scalar readBosanquetContactAngle
(
    const fvMesh& mesh,
    const volScalarField& alpha1,
    const volVectorField& U,
    const label patchi
)
{
    const volScalarField::Boundary& alphaBf = alpha1.boundaryField();

    if
    (
        !isA<alphaContactAngleTwoPhaseFvPatchScalarField>(alphaBf[patchi])
    )
    {
        FatalErrorInFunction
            << "Patch '" << mesh.boundary()[patchi].name() << "' in "
            << alpha1.name() << " is not an alpha contact-angle boundary "
            << "condition." << nl << exit(FatalError);
    }

    const alphaContactAngleTwoPhaseFvPatchScalarField& alphaCA =
        refCast<const alphaContactAngleTwoPhaseFvPatchScalarField>
        (
            alphaBf[patchi]
        );

    tmp<volVectorField> tGradAlpha = fvc::grad(alpha1, "nHat");
    surfaceVectorField gradAlphaf(fvc::interpolate(tGradAlpha()));

    const dimensionedScalar deltaN
    (
        "deltaN",
        1e-8/cbrt(average(mesh.V()))
    );

    surfaceVectorField nHat
    (
        IOobject
        (
            "nHatBosanquet",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        gradAlphaf/(mag(gradAlphaf) + deltaN)
    );

    tmp<scalarField> tTheta =
        alphaCA.theta
        (
            U.boundaryField()[patchi],
            nHat.boundaryField()[patchi]
        );

    const scalarField& theta = tTheta();
    const scalarField& magSf = mesh.magSf().boundaryField()[patchi];
    const scalar area = gSum(magSf);

    if (area <= SMALL)
    {
        FatalErrorInFunction
            << "Contact-angle patch '" << mesh.boundary()[patchi].name()
            << "' has zero area." << nl << exit(FatalError);
    }

    return gSum(theta*magSf)/area;
}

// ************************************************************************* //
