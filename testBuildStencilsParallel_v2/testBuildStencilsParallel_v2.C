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

    Points to discuss with Philip:

    - Right now i build stencils for faces. Maybe it make sense to build
      stencils for cells and boundary faces. Internal faces stencil is then
      combination of face cells stencil. This may reduce load becouse
      mesh have much more faces than cells.

    - When using layer aproach to build stencil i loop over all remote cells to
      find which one are the candidates. Maybe it make sense to exchange also
      cellCells data to be able to loop in the similar manner.

    - When collecting remote cells coordinate I'm looping over all cells. Maybe
      it is faster to construct a tree and loop over tree (and maybe later to
      reuse that tree for stencil construction)

Author
    Ivan Batistic, UCD.
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Pstream.H"
#include "PstreamBuffers.H"
#include "treeBoundBox.H"
#include "boundBox.H"
#include "DynamicList.H"
#include "HashSet.H"
#include "HashTable.H"
#include "PtrList.H"
#include "Map.H"
#include "labelPair.H"
#include "processorPolyPatch.H"
#include "globalIndex.H"


// * * * * * * * * * * * * * * Helper Functions  * * * * * * * * * * * * * * //


static labelList collectCellsUsingIndexedOctree
(
    const label faceI,
    const fvMesh& mesh,
    const label N,
    const labelList& remoteGlobalCells,
    const vectorField& remoteCentres,
    const globalIndex& globalCells
)
{
    // To be added later for testing purpose. Maybe it can be faster than
    // layer approach.

    labelList faceStencil(N);

    return faceStencil;
}


static labelList collectCellsUsingLayers
(
    const label faceI,
    const fvMesh& mesh,
    const label N,
    const labelList& remoteGlobalCells,
    const vectorField& remoteCentres,
    const globalIndex& globalCells,
    const boolList& procPatchesCells
)
{
    // Number of cell layers fetched to get N closest cells
    // Number of layers is controled in terms of overSampleFactor
    label overSampleFactor = 4;

    const label nbOfCandidates = N * overSampleFactor;

    // When building stencils, stencil is larger than N if there are cells with
    // equall distance. We use relTol to detect such cases. This can help to
    // have symmetric stencils on structured meshes.
    const scalar relTol = 1e-10;

    // Get potential candidates using complete layers looping.
    // Cell ID is in local indexing format.
    labelHashSet localCandidates;

    localCandidates.clear();

    labelHashSet currentLayer;
    labelHashSet nextLayer;

    const label ownerCell = mesh.faceOwner()[faceI];
    currentLayer.insert(ownerCell);
    localCandidates.insert(ownerCell);

    if (faceI < mesh.nInternalFaces())
    {
        const label nei = mesh.faceNeighbour()[faceI];
        currentLayer.insert(nei);
        localCandidates.insert(nei);
    }

    const labelListList& cellCells = mesh.cellCells();

    while
    (
        !currentLayer.empty()
     && localCandidates.size() < nbOfCandidates
    )
    {
        nextLayer.clear();

        forAllConstIter(labelHashSet, currentLayer, it)
        {
            const label c = it.key();
            const labelList& neigbours = cellCells[c];

            forAll(neigbours, i)
            {
                const label cellID = neigbours[i];
                if (localCandidates.insert(cellID))
                {
                    nextLayer.insert(cellID);
                }
            }
        }

        currentLayer.transfer(nextLayer);
    }

    // Above algorithm works only on owned faces since for remote cells we did
    // not transfered cellCells structure.
    // Probably it make sense to transfer cellCells structure, to be done later.
    //
    // Algorithm: 1) Sort current candidates
    //            2) Loop over existing candidates and check if some cell is
    //               at processor boundary. Looping is done up to N position.
    //            3) Get radius of the cell that is at position 1.5*N
    //            4) Check if any remote cell fits into this radius, in case it
    //               is inside this radius add it to the candidates list.

    // Phase 1: Build distance list for local candidates
    //          Using squared distance for efficiency
    const vector faceCentre = mesh.Cf()[faceI];
    const vectorField& C = mesh.C();

    labelList localList(localCandidates.size());
    {
        label k = 0;
        forAllConstIter(labelHashSet, localCandidates, it)
        {
            localList[k++] = it.key();
        }
    }

    List<Tuple2<label, scalar>> localDist(localList.size());

    forAll(localList, i)
    {
        const label cI = localList[i];
        localDist[i] = Tuple2<label, scalar>(cI, magSqr(C[cI] - faceCentre));
    }

    Foam::stableSort
    (
        localDist,
	[](auto& A, auto& B)
	{
	    return A.second() < B.second();
	}
    );

    // Phase 2: Check first N local candidates for processor boundary touch
    //
    bool localStencil = true;
    const label nCheck = min(N, localDist.size());

    for (label pos = 0; pos < nCheck; ++pos)
    {
	const label& cellID = localDist[pos].first();
	if (procPatchesCells[cellID])
	{
	    localStencil = false;
	    break;
	}
    }

    // If stencil is local we can already return the stencil
    if (localStencil)
    {
	if ( localList.size() < N )
	{
	    FatalErrorInFunction
		<< "Local candidates for stencil have size of: "
		<< localList.size()
		<< " but required minimum  stencil size is larger: " << N << nl
		<< "Increase the mesh size!"
		<< exit(FatalError);
	}

	// Enlarge stencil to include cells with the same distance
        const scalar cut = localDist[N-1].second();
        const scalar cutTol = cut*(1.0 + relTol);

        label nPick = N;
        while (nPick < localDist.size() && localDist[nPick].second() <= cutTol)
        {
            ++nPick;
        }

        labelList stencil(nPick);
        for (label i = 0; i < nPick; ++i)
        {
            const label localCellID = localDist[i].first();
            stencil[i] = globalCells.toGlobal(localCellID);
        }

        return stencil;
    }

    // Phase 3: Radius from position posR = ceil(1.5*N)
    label posR =  N + N/2 + (N%2);
    posR = min(max(posR, 0), localDist.size()-1);

    const scalar R2 = localDist[posR].second();

    DynamicList<Tuple2<label, scalar>> allDist;
    allDist.reserve((posR + 1) + 256);

    // Local cells within extended stencil radius squared R2
    for (label i = 0; i <= posR; ++i)
    {
        const label cellID = localDist[i].first();
        const label globalCellID = globalCells.toGlobal(cellID);
        allDist.append(Tuple2<label, scalar>(globalCellID, localDist[i].second()));
    }

    // Add remote cells within radius squared R2
    forAll(remoteGlobalCells, i)
    {
        const scalar d2 = magSqr(remoteCentres[i] - faceCentre);
        if (d2 <= R2)
        {
            allDist.append(Tuple2<label, scalar>(remoteGlobalCells[i], d2));
        }
    }

    // Sort again local and remote cells together
    Foam::stableSort
    (
        allDist,
        [](const Tuple2<label, scalar>& A, const Tuple2<label, scalar>& B)
        {
            return A.second() < B.second();
        }
    );

    // Check minimum stencil size
    if ( allDist.size() < N )
    {
	FatalErrorInFunction
	    << "Candidates for stencil have size of: " << allDist.size()
	    << " but required minimum stencil size is larger: " << N << nl
	    << "Increase the mesh size!"
	    << exit(FatalError);
    }

    // Enlarge stencil to include cells with the same distance
    const scalar cut = allDist[N-1].second();
    const scalar cutTol = cut*(1.0 + relTol);

    label nPick = N;
    while (nPick < allDist.size() && allDist[nPick].second() <= cutTol)
    {
	++nPick;
    }

    labelList faceStencil(nPick);
    for (label i = 0; i < nPick; ++i)
    {
	faceStencil[i] = allDist[i].first();
    }

    return faceStencil;
}


static labelList buildFacesStencil
(
    const label faceI,
    const fvMesh& mesh,
    const globalIndex& globalCells,
    const labelList& remoteGlobalCells,
    const vectorField& remoteCentres,
    const label N
)
{
    labelList  faceStencil;

    if (true)
    {
	// Flag cells at processor boundary, this flag is used to check if
	// stencil is on owned cells or it is necesary to loop over remote cells
	boolList procPatchesCells(mesh.nCells(), false);

	forAll(mesh.boundaryMesh(), patchI)
	{
	    const polyPatch& pp = mesh.boundaryMesh()[patchI];
	    if (isA<processorPolyPatch>(pp))
	    {
		const labelUList& faceCells = pp.faceCells();
		forAll(faceCells, i)
		{
		    procPatchesCells[faceCells[i]] = true;
		}
	    }
	}

	faceStencil  =
	    collectCellsUsingLayers
	    (
	        faceI,
		mesh,
		N,
		remoteGlobalCells,
		remoteCentres,
		globalCells,
		procPatchesCells
	    );
    }
    else
    {
	faceStencil  =
	    collectCellsUsingIndexedOctree
	    (
	        faceI,
		mesh,
		N,
		remoteGlobalCells,
		remoteCentres,
		globalCells
	    );
    }

    return faceStencil;
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
        "Collapsed (box-based) stencil marking + remote cell ID cache builder"
        "Writes processor stencils stencilProc0..stencilProcN-1 as "
	"volScalarFields."
    );

    argList::addOption
    (
        "debug", "label", "Debug level"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    label debug = 0;
    if (args.found("debug"))
    {
	debug = readLabel(args["debug"]);
    }

    Info<< "Debug level = " << debug << nl << endl;
    Info<< "Testing on nProcs: " << Pstream::nProcs() << endl;

    // ---------------------------------------------------------------------- //
    // 0. step - Basic preliminaries
    // ---------------------------------------------------------------------- //

    const vectorField& C = mesh.C();
    const scalarField& V = mesh.V();
    const label& nLocalCells = mesh.nCells();

    if (debug)
    {
	Pout<< "Processor: " << Pstream::myProcNo() << " owns "
	    << nLocalCells << " cells" << endl;
    }

    // ---------------------------------------------------------------------- //
    // 1. step - calculate and exchange halo thickness. Total halo depth is
    //           average depth of first layer multiplied with 5
    //           Mulitplication with 5 avoid calculating 5 layer depth but
    //           insert assumption that cells are uniform near proc boundary
    // ---------------------------------------------------------------------- //


    // Halo cells that owns this processor, we will take average over multiple
    // processor boundaries
    labelHashSet haloCells;
    scalar procPatchArea = 0.0;
    scalar haloVol = 0.0;

    forAll(mesh.boundaryMesh(), patchI)
    {
	const polyPatch& pp = mesh.boundaryMesh()[patchI];
	if (isA<processorPolyPatch>(pp))
	{
	    const vectorField& Sf = pp.faceAreas();
	    const labelUList& faceCells = pp.faceCells();
	    forAll(faceCells, i)
	    {
		haloCells.insert(faceCells[i]);
		procPatchArea += mag(Sf[i]);
		haloVol += V[faceCells[i]];
	    }
	}
    }

    // Halo depth on this processor
    const scalar haloDepthLocal  =
    	(procPatchArea > 0 ? haloVol/procPatchArea : 0.0);

    if (debug)
    {
	Pout<< "Processor: " << Pstream::myProcNo() << " halo thicnkess depth "
	    << "on this processor: " << haloDepthLocal
	    << ", halo volume: " << haloVol << ", procPatch area: "
	    << procPatchArea << endl;
    }

    // We need to exchance halo depths, becouse this processor needs halo depth
    // on neigbouring processors
    labelHashSet nbrSet;
    forAll(mesh.boundaryMesh(), patchI)
    {
	const polyPatch& pp = mesh.boundaryMesh()[patchI];
	if (isA<processorPolyPatch>(pp))
	{
	    const processorPolyPatch& ppp =
		refCast<const processorPolyPatch>(pp);
	    nbrSet.insert(ppp.neighbProcNo());
	}
    }

    labelList nbrs(nbrSet.toc());

    // Point to point exchange of haloDepth
    List<scalar> neighHaloDepth(Pstream::nProcs(), 0.0);

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
    forAll(nbrs, i)
    {
	UOPstream os(nbrs[i], pBufs);
	os << haloDepthLocal;
    }
    pBufs.finishedSends();

    forAll(nbrs, i)
    {
	const label from = nbrs[i];
	if (pBufs.recvDataCount(from))
	{
	    UIPstream is(from, pBufs);
	    scalar lenN;
	    is >> lenN;
	    neighHaloDepth[from] = lenN;
	}
    }

    const scalar haloDepth = 5 * max(neighHaloDepth);

    if (debug)
    {
	Pout<< "Processor: " << Pstream::myProcNo() << " halo thicnkess depth "
	    << "on neigbours multiplied with scaling: " << haloDepth << endl;
    }

    // ---------------------------------------------------------------------- //
    // 2. step - build treeBoundBox for owned faces (inflated with halo depth)
    // ---------------------------------------------------------------------- //

    treeBoundBox ownedFacesBox;
    {
	vector minPt(GREAT, GREAT, GREAT);
	vector maxPt(-GREAT, -GREAT, -GREAT);

	forAll(mesh.faces(), faceI)
	{
	    const point& faceCentre = mesh.faceCentres()[faceI];
	    minPt = min(minPt, faceCentre);
	    maxPt = max(maxPt, faceCentre);
	}

	// Expand domain for halo thickness
	minPt -= vector(haloDepth, haloDepth, haloDepth);
	maxPt += vector(haloDepth, haloDepth, haloDepth);

	ownedFacesBox = treeBoundBox(minPt, maxPt);
    }

    // ---------------------------------------------------------------------- //
    // 3. step - build treeBoundBox for owned cells and build global list of
    //           treeBoundBox
    // ---------------------------------------------------------------------- //

    treeBoundBox ownedCellsBox;
    {
	vector minPt(GREAT, GREAT, GREAT);
	vector maxPt(-GREAT, -GREAT, -GREAT);

	forAll(C, celli)
	{
	    minPt = min(minPt, C[celli]);
	    maxPt = max(maxPt, C[celli]);
	}

        ownedCellsBox = treeBoundBox(minPt, maxPt);
    }

    List<treeBoundBox> allOwnedCellsBox(Pstream::nProcs());
    allOwnedCellsBox[Pstream::myProcNo()] = ownedCellsBox;

    Pstream::allGatherList(allOwnedCellsBox);

    if (debug)
    {
	Info<< "Gathered processor cell boxes from all ranks." << nl << endl;
    }

    // ---------------------------------------------------------------------- //
    // 4. step - decide which processor will comunicate with this one
    // ---------------------------------------------------------------------- //

    DynamicList<label> procToQuery;

    for (label proc = 0; proc < Pstream::nProcs(); ++proc)
    {
	if (proc == Pstream::myProcNo())
	{
	    continue;
	}
	if (allOwnedCellsBox[proc].overlaps(ownedFacesBox))
	{
	    procToQuery.append(proc);
	}
    }

    if (debug)
    {
	Pout << "Processors " << Pstream::myProcNo() << " will query: "
	     << procToQuery.size() << " processors: " << procToQuery << endl;
    }

    // ---------------------------------------------------------------------- //
    // 5. step - send ownedFacesBox to other processors,
    //           receive boxes from others and respond with cell lists
    // ---------------------------------------------------------------------- //

    // We will store remote cell ID with global indexing
    globalIndex globalCells(mesh.nCells());

    // Inverse map for writing fields (which processors are using this one data)
    List<labelHashSet> cellsUsedByProc;

    //  Needed remote cells (with global indexing)
    List<labelList> remoteGlobalCellsPerProc;

    cellsUsedByProc.setSize(Pstream::nProcs());
    remoteGlobalCellsPerProc.setSize(Pstream::nProcs());

    forAll(remoteGlobalCellsPerProc, procI)
    {
	cellsUsedByProc[procI].clear();
	remoteGlobalCellsPerProc[procI].clear();
    }

    // Phase 1: Exchange treeBoundBoxes between processors
    Map<treeBoundBox> incomingBoxesFromProc;
    {
	PstreamBuffers sBufs(Pstream::commsTypes::nonBlocking);

	forAll(procToQuery, i)
	{
	    const label toProc = procToQuery[i];
	    UOPstream os(toProc, sBufs);
	    os << ownedFacesBox;
	}

	sBufs.finishedSends();

	for (label p = 0; p < Pstream::nProcs(); ++p)
	{
	    if (!sBufs.recvDataCount(p) || p == Pstream::myProcNo())
	    {
		continue;
	    }
	    UIPstream is(p, sBufs);
	    treeBoundBox qb;
	    is >> qb;

	    incomingBoxesFromProc.insert(p, qb);
	}
    }

    // Phase 2: Mark cells in overlaping boxes and fill out remoteGlobalCells list
    {
	PstreamBuffers rBufs(Pstream::commsTypes::nonBlocking);

	forAllConstIter(Map<treeBoundBox>, incomingBoxesFromProc, it)
	{
	    const label sender = it.key();
	    const treeBoundBox& qb = it();

	    labelHashSet& usedBySender = cellsUsedByProc[sender];
	    forAll(C, cellI)
	    {
		if (qb.contains(C[cellI]))
		{
		    usedBySender.insert(cellI);
		}
	    }

	    // Send back as a compact list
	    labelList markedCells(usedBySender.size());
	    label i = 0;
	    forAllConstIter(labelHashSet, usedBySender, iter)
	    {
		const label localCell = iter.key();
		markedCells[i++] = globalCells.toGlobal(localCell);
	    }

	    //Foam::sort(markedCells);
	    UOPstream os(sender, rBufs);
	    os << markedCells;
	}

	rBufs.finishedSends();

	forAll(procToQuery, i)
	{
	    const label fromProc = procToQuery[i];
	    if (!rBufs.recvDataCount(fromProc))
	    {
		continue;
	    }
	    UIPstream is(fromProc, rBufs);
	    labelList lst;
	    is >> lst;

	    remoteGlobalCellsPerProc[fromProc] = lst;
	}
    }

    if (debug)
    {
        label totalNeeded = 0;
        forAll(procToQuery, i)
        {
            totalNeeded += remoteGlobalCellsPerProc[procToQuery[i]].size();
        }

        Pout<< "Proc " << Pstream::myProcNo()
            << " queried " << procToQuery.size()
            << " procs, received total remote cells: " << totalNeeded
            << nl;
    }

    // ---------------------------------------------------------------------- //
    // 6. step - write processors stencils for visualisation purposes
    // ---------------------------------------------------------------------- //

    writeStencilFields(mesh, cellsUsedByProc);

    // ---------------------------------------------------------------------- //
    // 7. step - Construct one list of all halo cells using global cell IDs
    // ---------------------------------------------------------------------- //

    label nRemote = 0;
    forAll(remoteGlobalCellsPerProc, p)
    {
	nRemote += remoteGlobalCellsPerProc[p].size();
    }
    labelList remoteGlobalCells(nRemote);

    label k = 0;

    forAll(remoteGlobalCellsPerProc, p)
    {
	const labelList& ID = remoteGlobalCellsPerProc[p];

	forAll(ID, i)
	{
	    remoteGlobalCells[k] = ID[i];
	    ++k;
	}
    }

    // ---------------------------------------------------------------------- //
    // 8. step - Get remote cells coordinates
    // ---------------------------------------------------------------------- //

    List<vectorField> remoteCellCentresPerProc;
    remoteCellCentresPerProc.setSize(Pstream::nProcs());

    {
	// Phase 1: send globalIds to each processor in contact
	PstreamBuffers reqBufs(Pstream::commsTypes::nonBlocking);

	for (label p = 0; p < Pstream::nProcs(); ++p)
	{
	    if (p == Pstream::myProcNo() || remoteGlobalCellsPerProc[p].empty())
	    {
		continue;
	    }
	    const labelList& procGlobaCellID = remoteGlobalCellsPerProc[p];

	    UOPstream os(p, reqBufs);
	    os << procGlobaCellID;
	}

	reqBufs.finishedSends();

	// Phase 2: respond with cell centre coordinates
	PstreamBuffers repBufs(Pstream::commsTypes::nonBlocking);

	for (label fromProc = 0; fromProc < Pstream::nProcs(); ++fromProc)
	{
	    if
	    (
	        fromProc == Pstream::myProcNo()
	     || !reqBufs.recvDataCount(fromProc)
	    )
	    {
		continue;
	    }

	    UIPstream is(fromProc, reqBufs);

	    labelList procGlobaCellID;
	    is >> procGlobaCellID;

	    vectorField centres(procGlobaCellID.size());


	     forAll(procGlobaCellID, i)
	     {
		 const label globaCellID = procGlobaCellID[i];

		 const label localCell = globalCells.toLocal(globaCellID);

		 if (localCell < 0 || localCell >= mesh.nCells())
		 {
		     FatalErrorInFunction
			 << "Invalid global->local mapping: globaCellID="
			 << globaCellID << " localCell=" << localCell
			 << " on proc " << Pstream::myProcNo() << nl
			 << exit(FatalError);
		 }

		 centres[i] = C[localCell];
	     }

	     UOPstream os(fromProc, repBufs);
	     os << centres;
	}

	repBufs.finishedSends();

	// Phase 3: Send back populated centres lists
	for (label p = 0; p < Pstream::nProcs(); ++p)
	{
	    if (p == Pstream::myProcNo() || remoteGlobalCellsPerProc[p].empty())
	    {
		continue;
	    }

	    if (!repBufs.recvDataCount(p))
	    {
		FatalErrorInFunction
		    << "Did not receive centres from proc " << p << nl
		    << exit(FatalError);
	    }

	    UIPstream is(p, repBufs);

	    vectorField centres;
	    is >> centres;

	    if (centres.size() != remoteGlobalCellsPerProc[p].size())
	    {
		FatalErrorInFunction
		    << "Centres reply size mismatch from proc " << p
		    << ": got " << centres.size()
		    << " expected " << remoteGlobalCellsPerProc[p].size() << nl
		    << exit(FatalError);
	    }

	    remoteCellCentresPerProc[p].transfer(centres);
	}

	if (debug)
	{
	    label nTot = 0;
	    for (label p = 0; p < Pstream::nProcs(); ++p)
	    {
		nTot += remoteCellCentresPerProc[p].size();
	    }

	    Pout<< "Proc " << Pstream::myProcNo()
		<< " received remote centres for " << nTot
		<< " cells (sum over procs)" << nl;
	}
    }

    // ---------------------------------------------------------------------- //
    // 9. step - Construct one list of all halo cells centres
    // ---------------------------------------------------------------------- //

    nRemote = 0;
    forAll(remoteCellCentresPerProc, p)
    {
	nRemote += remoteCellCentresPerProc[p].size();
    }
    vectorField remoteCellCentres(nRemote);

    k = 0;

    forAll(remoteCellCentresPerProc, p)
    {
	const vectorList& cellCentre = remoteCellCentresPerProc[p];

	forAll(cellCentre, i)
	{
	    remoteCellCentres[k] = cellCentre[i];
	    ++k;
	}
    }

    if (debug)
    {
	if (remoteCellCentresPerProc.size() == remoteGlobalCellsPerProc.size())
	{
	    FatalErrorInFunction
		<< "The size of cell centre list and cellGlobalID list are "
		<< "not the same. Something went wrong."
		<< exit(FatalError);
	}

	forAll(remoteGlobalCellsPerProc, p)
	{
	    if
	    (
	        remoteGlobalCellsPerProc[p].size()
	     != remoteCellCentresPerProc[p].size()
	    )
	    {
		FatalErrorInFunction
		    << "The size of cellCentre list and cellID list for proc: "
		    << p << " are not matching. Something went wrong"
		    << exit(FatalError);
	    }
	}
    }

    // ---------------------------------------------------------------------- //
    // 10. step - Construct stencils using sphere bounded layer technique
    //            or indexedOctree.
    // ---------------------------------------------------------------------- //

    // Number of cells in stencil
    const label N = 30;

    // Number of required layer to be fetched depends on:
    //    a) the number of cells (interpolation order)
    //    b) mesh type (2D/3D)
    //    c) cells type
    //    d) cells aspect ratio

    // Here I will use overSampleFactor = 4...
    // We should chat about this.

    List<labelList> faceStencil(mesh.nFaces());

    // Build stencil for internal faces
    for (label faceI = 0; faceI < mesh.nInternalFaces(); ++faceI)
    {
	faceStencil[faceI] =
	    buildFacesStencil
	    (
	        faceI,
		mesh,
		globalCells,
		remoteGlobalCells,
		remoteCellCentres,
		N
	    );
    }

    // Build stencil for boundary faces
    // Skip empty faces and non-owned processor faces to avoid double counting
    forAll(mesh.boundaryMesh(), patchI)
    {
	const polyPatch& pp = mesh.boundaryMesh()[patchI];

	if (isA<emptyPolyPatch>(pp))
	{
	    continue;
	}

	if (isA<processorPolyPatch>(pp))
	{
	    const processorPolyPatch& ppp =
		refCast<const processorPolyPatch>(pp);

	    if (!ppp.owner())
	    {
		continue;
	    }
	}

        forAll(pp, i)
        {
            const label faceI = pp.start() + i;

            faceStencil[faceI] =
                buildFacesStencil
                (
		    faceI,
		    mesh,
		    globalCells,
		    remoteGlobalCells,
		    remoteCellCentres,
		    N
                );
        }
    }


    Info<< "End." << endl;
    return 0;
}


// ************************************************************************* //
