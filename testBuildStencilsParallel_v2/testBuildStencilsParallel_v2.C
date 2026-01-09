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


static void collectCellsUsingIndexedOctree
(
    const label faceI,
    const fvMesh& mesh,
    const label nbOfCandidates,
    labelHashSet& candidates,
    const labelList& remoteGlobalCells,
    const vectorField& remoteCentres,
    const globalIndex& globalCells
)
{
    // To be added later for testing purpose. Maybe it can be faster than
    // layer approach.
}


static void collectCellsUsingLayers
(
    const label faceI,
    const fvMesh& mesh,
    const label nbOfCandidates,
    labelHashSet& candidates,
    const labelList& remoteGlobalCells,
    const vectorField& remoteCentres,
    const globalIndex& globalCells
)
{
    // Insert face owner
    const label ownerCell = mesh.faceOwner()[faceI];
    candidates.clear();
    candidates.insert(ownerCell);

    // Insert face neighbour
    if (faceI < mesh.nInternalFaces())
    {
	const label neiCell = mesh.faceNeighbour()[faceI];
	candidates.insert(neiCell)
    }

    const labelListList& cellCells = mesh.cellCells();

    labelHashSet prevLayer;

    // Store layer per layer until required number of candidates is reached
    while (candidates .size() < nbOfCandidates)
    {
       const labelList& currentLayer = cellCells[cellI];

       // Add first layer of cells
       forAll(curCellCells, cI)
       {
           stencilCells.insert(curCellCells[cI]);
           prevLayer.insert(curCellCells[cI]);
       }
    }


}


static labelList buildFacesStencil
(
    const label faceI,
    const fvMesh& mesh,
    const globalIndex& globalCells,
    const labelList& remoteGlobalCells,
    const vectorField& remoteCentres,
    const label N,
    const label overSampleFactor
)
{
    // Phase 1: Get potential candidates
    labelHashSet candidates;

    if (true)
    {
	const label nbOfCandidates = N * overSampleFactor;

	collectCellsUsingLayers
	(
	     faceI,
	     mesh,
	     nbOfCandidates,
	     candidates,
	     remoteGlobalCells,
	     remoteCentre,
	     globalCells
        );
    }
    else
    {
	collectCellsUsingIndexedOctree
	(
	     faceI,
	     mesh,
	     nbOfCandidates,
	     candidates,
	     remoteGlobalCells,
	     remoteCentre,
	     globalCells
        );
    }

    // Phase 2: Filter false candidates using distance
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
    const vectorField& Cf = mesh.Cf();
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

    if (ownedFaces.empty())
    {
	// Empty rank - something went wrong?
	ownedFacesBox = treeBoundBox(vector::zero, vector::zero);
    }
    else
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
    vectorList remoteCellCentres(nRemote);

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

    // Number of cell layers fetched to get N closest cells
    // Number of layers is controled in terms of overSampleFactor
    label overSampleFactor = 4;

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
		remoteCentres,
		N,
		overSampleFactor
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
		    remoteCentres,
		    N,
		    overSampleFactor
                );
        }
    }


    Info<< "End." << endl;
    return 0;
}


// ************************************************************************* //
