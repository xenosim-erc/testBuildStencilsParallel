/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    testBuildStencilsParallel

Description
    Simple utility for testing the construction of face stencils in parallel
    as a proof of concept for least squares higher-order finite volume solvers.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Pstream.H"
#include "PstreamBuffers.H"
#include "treeBoundBox.H"
#include "boundBox.H"
#include "DynamicList.H"
#include "HashSet.H"
#include "PtrList.H"
#include "Map.H"
#include "labelPair.H"
#include "processorPolyPatch.H"


// * * * * * * * * * * * * * * Helper Functions  * * * * * * * * * * * * * * //


// Local-only bound box of owned cell centres (NO parallel reduction)
static boundBox localCellCentreBounds(const fvMesh& mesh)
{
    boundBox bb
    (
        vector(GREAT, GREAT, GREAT),
        vector(-GREAT, -GREAT, -GREAT)
    );

    const vectorField& C = mesh.C();
    forAll(C, celli)
    {
        bb.min() = min(bb.min(), C[celli]);
        bb.max() = max(bb.max(), C[celli]);
    }

    // Handle empty (possible in pathological decompositions)
    if (bb.min().x() > bb.max().x())
    {
        const vector z(0, 0, 0);
        return boundBox(z, z);
    }

    return bb;
}

static treeBoundBox localExpandedProcBox
(
    const fvMesh& mesh, const scalar radius
)
{
    boundBox bb = localCellCentreBounds(mesh);

    treeBoundBox tbb(bb.min(), bb.max());
    tbb.min() -= vector(radius, radius, radius);
    tbb.max() += vector(radius, radius, radius);
    return tbb;
}

static List<treeBoundBox> allGatherExpandedProcBoxes
(
    const fvMesh& mesh, const scalar radius
)
{
    List<treeBoundBox> procBoxes(Pstream::nProcs());
    procBoxes[Pstream::myProcNo()] = localExpandedProcBox(mesh, radius);

    // ESI v2412: convenient allGather
    Pstream::allGatherList(procBoxes);

    return procBoxes;
}

// Collect owned face centres (internal faces + optional processor patch faces
// on this rank)
static void buildOwnedFaceCentres
(
    const fvMesh& mesh,
    const bool includeInternalFaces,
    const bool includeProcPatchFaces,
    DynamicList<vector>& faceCentres
)
{
    faceCentres.clear();

    const vectorField& Cf = mesh.Cf();
    const label nInternal = mesh.nInternalFaces();

    if (includeInternalFaces)
    {
        for (label f = 0; f < nInternal; ++f)
        {
            faceCentres.append(Cf[f]);
        }
    }

    if (includeProcPatchFaces)
    {
        forAll(mesh.boundaryMesh(), patchi)
        {
            const polyPatch& pp = mesh.boundaryMesh()[patchi];
            if (isA<processorPolyPatch>(pp))
            {
                const label start = pp.start();
                forAll(pp, i)
                {
                    faceCentres.append(Cf[start + i]);
                }
            }
        }
    }
}

// Union a set of face-centre "query cubes" into a single box
static treeBoundBox unionFaceQueryBox
(
    const UList<vector>& faceCentres,
    const scalar radius
)
{
    treeBoundBox bb
    (
        vector(GREAT, GREAT, GREAT),
        vector(-GREAT, -GREAT, -GREAT)
    );

    const vector d(radius, radius, radius);

    forAll(faceCentres, i)
    {
        const vector& x = faceCentres[i];
        bb.min() = min(bb.min(), x - d);
        bb.max() = max(bb.max(), x + d);
    }

    // Empty safety
    if (bb.min().x() > bb.max().x())
    {
        const vector z(0, 0, 0);
        return treeBoundBox(z, z);
    }

    return bb;
}

// Build collapsed queries: one union box per target proc that overlaps that
// proc's expanded box.
// This is a "collapsed" comms pattern (coarse), typically producing a superset.
static void buildCollapsedQueriesToProc
(
    const List<treeBoundBox>& procBoxesExpanded,
    const treeBoundBox& myUnionQueryBox,
    Map<treeBoundBox>& queriesToProc
)
{
    queriesToProc.clear();

    const label myProc = Pstream::myProcNo();
    const label nProcs = procBoxesExpanded.size();

    for (label p = 0; p < nProcs; ++p)
    {
        if (p == myProc) continue;

        if (procBoxesExpanded[p].overlaps(myUnionQueryBox))
        {
            queriesToProc.insert(p, myUnionQueryBox);
        }
    }
}

// Mark local cells inside a query box (local IDs)
static void markCellsInBox
(
    const fvMesh& mesh,
    const treeBoundBox& box,
    labelHashSet& marked
)
{
    const vectorField& C = mesh.C();

    forAll(C, celli)
    {
        if (box.contains(C[celli]))
        {
            marked.insert(celli);
        }
    }
}

// Convert hashset to labelList (sorted for reproducibility)
static labelList toSortedLabelList(const labelHashSet& s)
{
    labelList lst(s.toc());
    Foam::sort(lst);
    return lst;
}

// Two-phase exchange:
//   Phase 1: send collapsed query boxes (one per target proc); receive boxes
//            from any proc.
//   Phase 2: respond with local cell IDs-in-box; receive responses for the
//            procs we queried.
//
// Outputs:
//   - cellsUsedByProc[p]  : local cells on THIS rank that belong to stencils of
//                           proc p (inverse map)
//   - neededRemoteCells[p]: on THIS rank, the REMOTE cell IDs (owned by proc p)
//                           that THIS proc needs

static void exchangeQueriesAndBuildCaches
(
    const fvMesh& mesh,
    const Map<treeBoundBox>& queriesToProc,
    const treeBoundBox& myUnionQueryBox,
    List<labelHashSet>& cellsUsedByProc,
    Map<labelList>& neededRemoteCells
)
{
    const label myProc = Pstream::myProcNo();
    const label nProcs = Pstream::nProcs();

    cellsUsedByProc.setSize(nProcs);
    for (label p = 0; p < nProcs; ++p) cellsUsedByProc[p].clear();

    neededRemoteCells.clear();

    // Also mark local cells for my own stencils (for output field
    // stencilProc<myProc>)
    {
        labelHashSet& mySet = cellsUsedByProc[myProc];
        if (!myUnionQueryBox.empty())
        {
            markCellsInBox(mesh, myUnionQueryBox, mySet);
        }
    }

    // --------------------------
    // Phase 1: send/recv queries
    // --------------------------
    Map<treeBoundBox> incomingQueryFromProc;

    {
        PstreamBuffers qBufs(Pstream::commsTypes::nonBlocking);

        // Send my collapsed queries
        forAllConstIter(Map<treeBoundBox>, queriesToProc, it)
        {
            const label p = it.key();
            UOPstream toProc(p, qBufs);
            toProc << it();
        }

        qBufs.finishedSends();

        // Receive queries from any proc
        for (label p = 0; p < nProcs; ++p)
        {
            if (p == myProc) continue;
            if (!qBufs.recvDataCount(p)) continue;

            UIPstream fromProc(p, qBufs);
            treeBoundBox qb;
            fromProc >> qb;

            incomingQueryFromProc.insert(p, qb);
        }
    }

    // ----------------------------------------
    // Phase 2: respond (cell IDs) + recv replies
    // ----------------------------------------
    {
        PstreamBuffers rBufs(Pstream::commsTypes::nonBlocking);

        // Respond to incoming queries: compute local hits and send back
        forAllConstIter(Map<treeBoundBox>, incomingQueryFromProc, it)
        {
            const label sender = it.key();
            const treeBoundBox& qb = it();

            labelHashSet& usedBySender = cellsUsedByProc[sender];
            markCellsInBox(mesh, qb, usedBySender);

            // Send back the local cell IDs to sender (these are "remote IDs"
            // from sender's POV)
            labelList hitCells = toSortedLabelList(usedBySender);

            UOPstream back(sender, rBufs);
            back << hitCells;
        }

        rBufs.finishedSends();

        // Receive replies for the procs we queried; cache as "needed remote
        // cells"
        forAllConstIter(Map<treeBoundBox>, queriesToProc, it)
        {
            const label p = it.key();
            if (!rBufs.recvDataCount(p)) continue;

            UIPstream fromProc(p, rBufs);
            labelList remoteCellsOnP;
            fromProc >> remoteCellsOnP;

            neededRemoteCells.insert(p, remoteCellsOnP);
        }
    }
}

// Build (proc, remoteCellID) -> ghostSlot cache and total ghost count
static void buildRemoteToGhostCache
(
    const Map<labelList>& neededRemoteCells,
    HashTable<label, labelPair, labelPair::Hash<>>& remoteToGhost,
    label& nGhost
)
{
    remoteToGhost.clear();
    nGhost = 0;

    forAllConstIter(Map<labelList>, neededRemoteCells, it)
    {
        const label proc = it.key();
        const labelList& cells = it();

        forAll(cells, i)
        {
            const labelPair key(proc, cells[i]);
            if (!remoteToGhost.found(key))
            {
                remoteToGhost.insert(key, nGhost++);
            }
        }
    }
}

static PtrList<volScalarField> writeStencilFields
(
    const fvMesh& mesh,
    const List<labelHashSet>& cellsUsedByProc
)
{
    const label nProcs = Pstream::nProcs();
    PtrList<volScalarField> fields(nProcs);

    for (label p = 0; p < nProcs; ++p)
    {
        word fldName("stencilProc" + Foam::name(p));

        fields.set
        (
            p,
            new volScalarField
            (
                IOobject
                (
                    fldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimless, 0.0)
            )
        );

        volScalarField& f = fields[p];

        forAllConstIter(labelHashSet, cellsUsedByProc[p], iter)
        {
            f[iter.key()] = 1.0;
        }

        f.correctBoundaryConditions();
        f.write();
    }

    return fields;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Collapsed (box-based) stencil marking + remote cell ID cache builder "
        "(ESI v2412). "
        "Writes stencilProc0..stencilProcN-1 as volScalarFields (1=cell used "
        "by proc's stencils)."
    );

    argList::addOption
    (
        "radius", "scalar", "Stencil radius (mesh units). Required."
    );
    argList::addBoolOption
    (
        "procFaces", "Also include processor patch faces in the owned face set."
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    if (!args.found("radius"))
    {
        FatalErrorInFunction
            << "Missing required option: -radius <scalar>" << nl
            << exit(FatalError);
    }

    const scalar radius = readScalar(args.lookup("radius")());
    const bool includeProcFaces = args.found("procFaces");

    if (radius <= SMALL)
    {
        FatalErrorInFunction
            << "radius must be > 0" << nl
            << exit(FatalError);
    }

    Info<< "radius            = " << radius << nl
        << "includeProcFaces  = " << includeProcFaces << nl
        << "nProcs            = " << Pstream::nProcs() << nl
        << "myProc            = " << Pstream::myProcNo() << nl << endl;

    // 1) Gather expanded proc boxes (Approach C candidate pruning)
    List<treeBoundBox> procBoxesExpanded
    (
        allGatherExpandedProcBoxes(mesh, radius)
    );

    // 2) Build owned face centres + a single union query box for this proc
    DynamicList<vector> ownedFaceCentres;
    buildOwnedFaceCentres
    (
        mesh, /*internal*/true, /*procPatch*/includeProcFaces, ownedFaceCentres
    );

    Info<< "Owned faces on proc " << Pstream::myProcNo()
        << " = " << ownedFaceCentres.size() << nl << endl;

    treeBoundBox myUnionQueryBox = unionFaceQueryBox(ownedFaceCentres, radius);

    // 3) Collapsed queries (one box per target proc)
    Map<treeBoundBox> queriesToProc;
    if (!myUnionQueryBox.empty())
    {
        buildCollapsedQueriesToProc
        (
            procBoxesExpanded, myUnionQueryBox, queriesToProc
        );
    }

    Info<< "Target procs queried from proc " << Pstream::myProcNo()
        << " = " << queriesToProc.size() << nl << endl;

    // 4) Exchange queries and build:
    //    - cellsUsedByProc (inverse map for writing fields)
    //    - neededRemoteCells (forward map: what THIS proc needs from each
    //      remote proc)
    List<labelHashSet> cellsUsedByProc;
    Map<labelList> neededRemoteCells;

    exchangeQueriesAndBuildCaches
    (
        mesh,
        queriesToProc,
        myUnionQueryBox,
        cellsUsedByProc,
        neededRemoteCells
    );

    // 5) Build remote-to-ghost cache (solver-ready)
    HashTable<label, labelPair, labelPair::Hash<>> remoteToGhost;
    label nGhost = 0;

    buildRemoteToGhostCache(neededRemoteCells, remoteToGhost, nGhost);

    Info<< "neededRemoteCells procs = " << neededRemoteCells.size() << nl
        << "nGhost (unique remote cells) = " << nGhost << nl << endl;

    // 6) Write stencil fields stencilProc0..stencilProcN-1
    writeStencilFields(mesh, cellsUsedByProc);

    Info<< "Done." << endl;
    return 0;
}


// ************************************************************************* //
