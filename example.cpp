static const char help[] = "Simple Tests";

#include <petsc.h>
#include <petsc/private/dmpleximpl.h>
#include <egadsTypes.h>
#include <egads.h>
#include "petscdmplex.h"
#include <petscviewerhdf5.h>
#include <iostream>



int main(int argc, char **argv) {
    PetscCall(PetscInitialize(&argc, &argv, NULL, help));

    PetscViewer petscViewer = nullptr;
    PetscViewerHDF5Open(PETSC_COMM_SELF, "/Users/mcgurn/scratch/testinput/_meshMappingTest/domain.hdf5", FILE_MODE_READ, &petscViewer);


    DM dmLoad;
    DMCreate(PETSC_COMM_SELF, &dmLoad);
    DMSetType(dmLoad, DMPLEX);
    DMLoad(dmLoad, petscViewer);
    DMView(dmLoad, PETSC_VIEWER_STDOUT_WORLD);

    // determine the number of cells in the
    PetscInt dmCells;
    // extract the connected cells and store them
    PetscInt cStart, cEnd;
    DMPlexGetHeightStratum(dmLoad, 0, &cStart, &cEnd);
    PetscInt numberCells = cEnd-cStart;

    DM dm;
    DMClone(dmLoad, &dm);

    // create an empty vec
    Vec glob;
    VecCreate(PETSC_COMM_WORLD, &glob);
    PetscViewerHDF5PushTimestepping(petscViewer);
    PetscCall(PetscObjectSetName((PetscObject)glob, "/cell_fields/solution_fieldB"));
    PetscCall(VecLoad(glob, petscViewer));

    // sanity checks
    PetscInt blockSize;
    VecGetBlockSize(glob, &blockSize);
    printf("Block Size: %d\n", blockSize);
    PetscInt vecSize;
    VecGetSize(glob, &vecSize);
    printf("Vec Size: %d\n", vecSize);
    printf("Vec cells: %d\n", (vecSize/blockSize));
    printf("Number Cells: %d\n", numberCells);

    // Setup the DM for this field
    PetscInt dim;
    DMGetDimension(dm, &dim);
    PetscBool simplex;
    DMPlexIsSimplex(dm, &simplex);
    PetscFE fe;
    PetscFECreateLagrange(PETSC_COMM_SELF, dim, blockSize, simplex, 0, PETSC_DETERMINE, &fe);
    DMSetField(dm, 0, NULL, (PetscObject)fe);
    PetscFEDestroy(&fe);
    DMCreateDS(dm);



    DMInterpolationInfo interpolant;
    DMInterpolationCreate(PETSC_COMM_SELF, &interpolant);
    DMInterpolationSetDim(interpolant, 2);
    DMInterpolationSetDof(interpolant, blockSize);

    // Copy over the np of particles
    PetscReal pt[2] = {0.0189287, 0.00845254};
    DMInterpolationAddPoints(interpolant, 1, pt);

    /* Particles that lie outside the domain should be dropped,
    whereas particles that move to another partition should trigger a migration */
    DMInterpolationSetUp(interpolant, dm, PETSC_FALSE, PETSC_TRUE);

    // Create a vec to hold the information
    Vec eulerianFieldAtParticles;
    VecCreateSeq(PETSC_COMM_SELF, 1 * blockSize, &eulerianFieldAtParticles);

    // interpolate
    PetscCall(DMInterpolationEvaluate(interpolant, dm, glob, eulerianFieldAtParticles));

    // Now cleanup
    DMInterpolationDestroy(&interpolant);
    VecView(eulerianFieldAtParticles, PETSC_VIEWER_STDOUT_WORLD);








//    VecView(glob, PETSC_VIEWER_STDOUT_WORLD);

//    // try to create global mesh
//    Vec glob;
//    DMCreateGlobalVector(dm, &glob);
//    VecView(glob, PETSC_VIEWER_STDOUT_WORLD);
//    PetscInt size;
//    VecGetSize(glob, &size);
//    printf("Size: %d", size);

//    VecCreate(PETSC_COMM_WORLD, &glob);
////    PetscViewerHDF5PushGroup(petscViewer, "cell_fields");
//    PetscBool has;
//    PetscCall(PetscViewerHDF5HasGroup(petscViewer, "cell_fields", &has));
//
//    PetscCall(PetscObjectSetName((PetscObject)glob, "solution_euler"));
//    PetscViewerHDF5PopGroup(petscViewer);
//
//    PetscCall(VecLoad(glob, petscViewer));
//
//

//    PetscViewerDestroy(&petscViewer);

    PetscFinalize();
    return 0;
}