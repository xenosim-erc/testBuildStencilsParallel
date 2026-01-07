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
#include "PtrList.H"
#include "Map.H"
#include "labelPair.H"
#include "processorPolyPatch.H"


// * * * * * * * * * * * * * * Helper Functions  * * * * * * * * * * * * * * //


static void writeStencilPerProcFields
(
    const fvMesh& mesh,
    const Time& runTime,
    const List<labelList>& cellsUsedByProc,
    const label debug = 0
)
{
    for (label p = 0; p < Pstream::nProcs(); ++p)
    {
        if (cellsUsedByProc[p].empty())
	{
	    continue;
	}

        const word fieldName("stencilProc_" + Foam::name(p));

        volScalarField f
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0.0)
        );

        const labelList& lst = cellsUsedByProc[p];
        forAll(lst, i)
        {
	    const label c = lst[i];
	    if (c >= 0 && c < mesh.nCells())
	    {
		f[c] = 1.0;
	    }
        }

        f.correctBoundaryConditions();
        f.write();

        if (debug)
        {
            Info<< "Wrote field: " << fieldName << " size=" << lst.size() << nl;
        }
    }
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

    //
    // 0. step - Get faces for which we will build stencil. On processor
    // boundary, processor that is face owner holds the face to avoid double
    // counting.
    //
    const vectorField& Cf = mesh.Cf();
    const vectorField& C = mesh.C();
    const scalarField& V = mesh.V();
    const label& nLocalCells = mesh.nCells();

    if (debug)
    {
	Pout<< "Processor: " << Pstream::myProcNo() << " owns "
	    << nLocalCells << " cells" << endl;
    }

    DynamicList<label> ownedFaces;
    {
	ownedFaces.reserve(mesh.nInternalFaces());

	// Add internal faces
	for (label faceI = 0; faceI < mesh.nInternalFaces(); ++faceI)
	{
	    ownedFaces.append(faceI);
	}

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
		    forAll(pp, i)
		    {
		        ownedFaces.append(pp.start() + i);
		    }
		}
	    }
	}

	if (debug)
	{
	    Pout<< "Processor: " << Pstream::myProcNo() << " owns "
		<< ownedFaces.size() << " faces" << endl;
	}
    }

    //
    // 1. step - calculate and exchange halo thickness
    //

    // Halo cells that owns this processor, we will take average over multiple
    // processor boundaries
    labelHashSet haloCells;
    forAll(mesh.boundaryMesh(), patchi)
    {
	const polyPatch& pp = mesh.boundaryMesh()[patchi];
	if (!isA<processorPolyPatch>(pp))
	{
	    continue;
	}

	const labelUList& faceCells = pp.faceCells();
	forAll(faceCells, i)
	{
	    haloCells.insert(faceCells[i]);
	}
    }

    scalar sumV = 0.0;
    label  nV   = 0;

    forAllConstIter(labelHashSet, haloCells, iter)
    {
	const label c = iter.key();
	sumV += V[c];
	++nV;
    }

    const scalar haloAvgVLocal =
	(nV > 0 ? sumV / scalar(nV) : 0.0);
    const scalar haloDepthLocal  =
	(haloAvgVLocal > 0 ? Foam::cbrt(haloAvgVLocal) : 0.0);

    labelHashSet nbrSet;
    forAll(mesh.boundaryMesh(), patchi)
    {
	const polyPatch& pp = mesh.boundaryMesh()[patchi];
	if (!isA<processorPolyPatch>(pp))
	{
	    continue;
	}
	const processorPolyPatch& ppp = refCast<const processorPolyPatch>(pp);
	nbrSet.insert(ppp.neighbProcNo());
    }

    labelList nbrs(nbrSet.toc());

    // Exchange haloDepth with other processors
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

    const scalar maxHaloDepth = max(neighHaloDepth);
    const scalar haloRadius = 5 * maxHaloDepth;

    //
    // 2. step - build treeBoundBox for processor faces
    //
    treeBoundBox ownedFacesBox;

    if (ownedFaces.empty())
    {
	// Empty rank
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
	minPt -= vector(haloRadius, haloRadius, haloRadius);
	maxPt += vector(haloRadius, haloRadius, haloRadius);

	ownedFacesBox = treeBoundBox(minPt, maxPt);
    }

    //
    // 3. step - build treeBoundBox for processor cells
    //

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

    // TO DO: This should be replaced with neighbor-only layering.
    List<treeBoundBox> allOwnedCellsBox(Pstream::nProcs());
    allOwnedCellsBox[Pstream::myProcNo()] = ownedCellsBox;

    Pstream::allGatherList(allOwnedCellsBox);

    if (debug)
    {
	Info<< "Gathered processor cell boxes from all ranks." << nl << endl;
    }

    //
    // 4. step - decide which processor will comunicate with this one
    //
    DynamicList<label> procToQuery;

    for (label proc = 0; proc < Pstream::nProcs(); ++proc)
    {
	if (proc == Pstream::myProcNo()) continue;

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

    // Send ownedFacesBox to other processors, receive boxes from others and
    // respond with cell lists
    List<labelList> replyToProc(Pstream::nProcs());
    {
	// Send my processor faces box to neighbouring processors
	PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

	forAll(procToQuery, i)
	{
	    const label& toProc = procToQuery[i];
	    UOPstream os(toProc, pBufs);
	    os << ownedFacesBox;
	}

        pBufs.finishedSends();

	// Recieve processor faces boxes and answer with list of cells inside
	for (label fromProc = 0; fromProc < Pstream::nProcs(); ++fromProc)
	{
	    if (fromProc == Pstream::myProcNo())
	    {
		continue;
	    }

	    if (pBufs.recvDataCount(fromProc))
	    {
		UIPstream is(fromProc, pBufs);
		treeBoundBox qb;
		is >> qb;

		DynamicList<label> hits;

		// Guess: reserve size is 10% of cells at this proc
		const label reserveSize  =
		    min(max(label(0.1*nLocalCells), 1024), 200000);

		hits.reserve(reserveSize);

		// This is location where we can acheve speed-up becouse right
		// now, we check the whole domain.
		forAll(C, cellI)
		{
		    if (qb.contains(C[cellI]))
		    {
			hits.append(cellI);
		    }
		}

		replyToProc[fromProc].transfer(hits);
	    }
	}
    }

    // Send and receieve cell ID lists from neighbouring processors
    List<labelList> remoteCells(Pstream::nProcs());
    {
	PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

	for (label toProc = 0; toProc < Pstream::nProcs(); ++toProc)
	{
	    if (toProc == Pstream::myProcNo())
	    {
		continue;
	    }

	    if (replyToProc[toProc].size())
	    {
		UOPstream os(toProc, pBufs);
		os << replyToProc[toProc];
	    }
	}

	pBufs.finishedSends();

	forAll(procToQuery, i)
	{
	    const label fromProc = procToQuery[i];

	    if (pBufs.recvDataCount(fromProc))
	    {
		UIPstream is(fromProc, pBufs);
		labelList lst;
		is >> lst;
		remoteCells[fromProc] = lst;
	    }
	}
    }

    if (debug)
    {
	label totalRemote = 0;
	forAll(remoteCells, proc)
	{
	    totalRemote += remoteCells[proc].size();
	}

	Pout<< "Total remote candidate cells received: " << totalRemote
	    <<  endl;
    }


    if (debug)
    {
	Pout <<  "Processors " << Pstream::myProcNo() << " have remote cells on"
	     << remoteCells.size() <<" processors. " << endl;
    }
    // Write processor stencils for visualisation purpose
    //writeStencilPerProcFields(mesh, runTime, replyToProc, debug);

    Info<< "End." << endl;
    return 0;
}


// ************************************************************************* //
