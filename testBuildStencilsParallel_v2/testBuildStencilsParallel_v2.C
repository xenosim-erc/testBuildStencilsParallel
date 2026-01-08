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
    // 0. step - Get faces for which we will build stencil. On processor
    // boundary, processor that is face owner holds the face to avoid double
    // counting.
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

    // I should use dynamic list only for processor faces? Why copy all faces?
    DynamicList<label> ownedFaces;
    {
	ownedFaces.reserve(mesh.nInternalFaces());

	// Add internal faces
	for (label faceI = 0; faceI < mesh.nInternalFaces(); ++faceI)
	{
	    ownedFaces.append(faceI);
	}

	// For debug purpose
	label nProcFacesOwnedLocal = 0;

	// Add processor boundary faces
	forAll(mesh.boundaryMesh(), patchI)
	{
	    const polyPatch& pp = mesh.boundaryMesh()[patchI];

	    if (isA<processorPolyPatch>(pp))
	    {
		const processorPolyPatch& ppp =
		    refCast<const processorPolyPatch>(pp);

		if (ppp.owner())
		{
		    nProcFacesOwnedLocal += pp.size();

		    forAll(pp, i)
		    {
		        ownedFaces.append(pp.start() + i);
		    }
		}
	    }
	}

	if (debug)
	{
	    const label nInternalLocal = mesh.nInternalFaces();
	    const label nOwnedLocal    = ownedFaces.size();

	    Pout<< "Processor: " << Pstream::myProcNo() << " owns "
		<< nOwnedLocal << " faces. Of which " << nInternalLocal
		<< " are internal faces." << endl;

	    const label nInternalGlobal =
		returnReduce(nInternalLocal, sumOp<label>());

	    const label nProcFacesOwnedGlobal =
		returnReduce(nProcFacesOwnedLocal, sumOp<label>());

	    const label nOwnedGlobal =
		returnReduce(nOwnedLocal, sumOp<label>());

	    Info<< nl << "Total number of owned faces: " << nOwnedGlobal << nl
		<< "Total number of owned processor faces: "
		<< nProcFacesOwnedGlobal
		<< nl << "Total numner of internal faces: " << nInternalGlobal
		<< endl;
	}
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

	forAll(ownedFaces, faceI)
	{
	    const vector& p = Cf[ownedFaces[faceI]];
	    minPt = min(minPt, p);
	    maxPt = max(maxPt, p);
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
    List<labelList> remoteCellsPerProc;

    cellsUsedByProc.setSize(Pstream::nProcs());
    remoteCellsPerProc.setSize(Pstream::nProcs());

    forAll(remoteCellsPerProc, procI)
    {
	cellsUsedByProc[procI].clear();
	remoteCellsPerProc[procI].clear();
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

    // Phase 2: Mark cells in overlaping boxes and fill out remoteCells list
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

	    remoteCellsPerProc[fromProc] = lst;
	}
    }

    if (debug)
    {
        label totalNeeded = 0;
        forAll(procToQuery, i)
        {
            totalNeeded += remoteCellsPerProc[procToQuery[i]].size();
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

    labelList remoteCells;

    forAll(remoteCells, p)
    {
	remoteCells.append(remoteCellsPerProc[p]);
    }
    //Foam::sort(remoteCellsGlobal);


    Info<< "End." << endl;
    return 0;
}


// ************************************************************************* //
