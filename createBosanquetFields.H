/*---------------------------------------------------------------------------*\
  Bosanquet model input, diagnostic fields and global connectivity update.
\*---------------------------------------------------------------------------*/

IOdictionary bosanquetProperties
(
    IOobject
    (
        "bosanquetProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const scalar connectivityAlpha
(
    bosanquetProperties.lookupOrDefault<scalar>("connectivityAlpha", 0.5)
);

if (connectivityAlpha <= 0 || connectivityAlpha >= 1)
{
    FatalErrorInFunction
        << "connectivityAlpha must lie strictly between 0 and 1; value = "
        << connectivityAlpha << nl << exit(FatalError);
}

const wordList seedPatches(bosanquetProperties.lookup("seedPatches"));

if (seedPatches.empty())
{
    FatalErrorInFunction
        << "seedPatches is empty in constant/bosanquetProperties. "
        << "At least one donor/reservoir patch is required." << nl
        << exit(FatalError);
}

const PtrList<dictionary> bosRegions
(
    bosanquetProperties.lookup("regions")
);

if (bosRegions.empty())
{
    FatalErrorInFunction
        << "No regions are defined in constant/bosanquetProperties."
        << nl << exit(FatalError);
}

autoPtr<surfaceTensionModel> bosanquetSigmaModel
(
    surfaceTensionModel::New(mixture, mesh)
);

volScalarField pCapBos
(
    IOobject
    (
        "pCapBos",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimPressure, 0)
);

volScalarField mainWaterMask
(
    IOobject
    (
        "mainWaterMask",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0)
);

volScalarField bosanquetRegionMask
(
    IOobject
    (
        "bosanquetRegionMask",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0)
);

surfaceScalarField mainMeniscusMask
(
    IOobject
    (
        "mainMeniscusMask",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0)
);

auto updateBosanquetFields = [&]()
{
    pCapBos.primitiveFieldRef() = scalar(0);
    bosanquetRegionMask.primitiveFieldRef() = scalar(0);

    tmp<volScalarField> tSigmaBos = bosanquetSigmaModel->sigma();
    const volScalarField& sigmaBos = tSigmaBos();

    forAll(bosRegions, regioni)
    {
        const dictionary& regionDict = bosRegions[regioni];
        const word regionName(regionDict.lookup("name"));
        const word wallPatchName(regionDict.lookup("wallPatch"));
        const word cellZoneName(regionDict.lookup("cellZone"));

        const dimensionedScalar radius("radius", dimLength, regionDict);

        if (radius.value() <= SMALL)
        {
            FatalErrorInFunction
                << "Region '" << regionName << "' has invalid radius "
                << radius << nl << exit(FatalError);
        }

        const label patchi = mesh.boundaryMesh().findPatchID(wallPatchName);
        const label zonei = mesh.cellZones().findZoneID(cellZoneName);

        if (patchi < 0)
        {
            FatalErrorInFunction
                << "Cannot find wallPatch '" << wallPatchName
                << "' for region '" << regionName << "'." << nl
                << exit(FatalError);
        }

        if (zonei < 0)
        {
            FatalErrorInFunction
                << "Cannot find cellZone '" << cellZoneName
                << "' for region '" << regionName << "'." << nl
                << exit(FatalError);
        }

        const scalar thetaDeg =
            readBosanquetContactAngle(mesh, alpha1, U, patchi);
        const scalar thetaRad =
            thetaDeg*constant::mathematical::pi/180.0;
        const cellZone& zone = mesh.cellZones()[zonei];

        forAll(zone, i)
        {
            const label celli = zone[i];
            pCapBos[celli] =
                2.0*sigmaBos[celli]*Foam::cos(thetaRad)/radius.value();
            bosanquetRegionMask[celli] = 1.0;
        }

        Info<< "Bosanquet region '" << regionName << "': wallPatch="
            << wallPatchName << ", cellZone=" << cellZoneName
            << ", cells=" << returnReduce(zone.size(), sumOp<label>())
            << ", R=" << radius.value() << " m, theta=" << thetaDeg
            << " deg" << nl;
    }

    pCapBos.correctBoundaryConditions();
    bosanquetRegionMask.correctBoundaryConditions();

    // Classify water cells for connectivity. This threshold defines numerical
    // connectivity only; it is not used to reconstruct the VOF interface.
    boolList waterCell(mesh.nCells(), false);
    labelList waterLabel(mesh.nCells(), 0);

    forAll(waterCell, celli)
    {
        waterCell[celli] = alpha1[celli] >= connectivityAlpha;
        waterLabel[celli] = waterCell[celli] ? 1 : 0;
    }

    // Block faces crossing the threshold. Water-water and air-air faces stay
    // open, so regionSplit produces distinct global water components.
    boolList blockedFace(mesh.nFaces(), false);
    const labelUList& owner = mesh.faceOwner();
    const labelUList& neighbour = mesh.faceNeighbour();

    forAll(neighbour, facei)
    {
        blockedFace[facei] =
            waterCell[owner[facei]] != waterCell[neighbour[facei]];
    }

    labelList neighbourWater;
    syncTools::swapBoundaryCellList(mesh, waterLabel, neighbourWater);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    forAll(patches, patchi)
    {
        const polyPatch& patch = patches[patchi];

        if (patch.coupled())
        {
            forAll(patch, patchFacei)
            {
                const label facei = patch.start() + patchFacei;
                const label boundaryFacei = facei - mesh.nInternalFaces();
                blockedFace[facei] =
                    waterLabel[owner[facei]]
                 != neighbourWater[boundaryFacei];
            }
        }
    }

    syncTools::syncFaceList(mesh, blockedFace, orEqOp<bool>());

    const regionSplit cellRegion(mesh, blockedFace);
    boolList seedRegion(cellRegion.nRegions(), false);

    forAll(seedPatches, seedi)
    {
        const label patchi = mesh.boundaryMesh().findPatchID(seedPatches[seedi]);

        if (patchi < 0)
        {
            FatalErrorInFunction
                << "Cannot find connectivity seed patch '"
                << seedPatches[seedi] << "'." << nl << exit(FatalError);
        }

        const labelUList& faceCells = mesh.boundary()[patchi].faceCells();
        forAll(faceCells, i)
        {
            const label celli = faceCells[i];
            if (waterCell[celli])
            {
                seedRegion[cellRegion[celli]] = true;
            }
        }
    }

    Pstream::listCombineReduce(seedRegion, orEqOp<bool>());

    label nSeedRegions = 0;
    forAll(seedRegion, regioni)
    {
        nSeedRegions += seedRegion[regioni] ? 1 : 0;
    }

    mainWaterMask.primitiveFieldRef() = scalar(0);
    label nLocalMainWaterCells = 0;

    forAll(waterCell, celli)
    {
        if (waterCell[celli] && seedRegion[cellRegion[celli]])
        {
            mainWaterMask[celli] = 1.0;
            ++nLocalMainWaterCells;
        }
    }

    mainWaterMask.correctBoundaryConditions();

    // A face is a Bosanquet face only if it touches the seed-connected water
    // component and lies in a configured Bosanquet cellZone. Multiplication
    // by snGrad(alpha1) in pEqn.H naturally restricts the actual force to the
    // diffuse water-air interface.
    mainMeniscusMask =
        pos(fvc::interpolate(mainWaterMask) - scalar(SMALL))
       *pos(fvc::interpolate(bosanquetRegionMask) - scalar(SMALL));
    mainMeniscusMask.correctBoundaryConditions();

    const label nMainWaterCells =
        returnReduce(nLocalMainWaterCells, sumOp<label>());

    if (!nSeedRegions || !nMainWaterCells)
    {
        WarningInFunction
            << "No water component connected to seedPatches at time "
            << runTime.timeName() << ". Bosanquet mask is zero; original CSF "
            << "is used everywhere." << nl;
        mainMeniscusMask = scalar(0);
    }

    Info<< "Connectivity: globalRegions=" << cellRegion.nRegions()
        << ", seededWaterRegions=" << nSeedRegions
        << ", mainWaterCells=" << nMainWaterCells
        << ", threshold=" << connectivityAlpha << nl << endl;
};

// ************************************************************************* //
