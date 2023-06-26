static char help[] = "Simplified example for conduction";

#include <petscdmplex.h>
#include <petscds.h>
#include <petscts.h>
#include <set>
#include <iostream>

/*
The heat transfer equation assuming constant density, specific heat and, conductivity.  There are no source/sink terms.

    \[\rho C \frac{\partial u}{\partial dt} - k \Delta u = 0\]

In weak form:

   \[ \int_{\Omega} \left[w \rho C \frac{\partial T}{\partial t} + k \frac{\partial w}{\partial x_i} \frac{\partial T}{\partial x_j}  \right] d\Omega + \int_{\Gamma} q w d\Gamma \]

Therefore the integrand for the test function w term is

   \[\rho C \frac{\partial T}{\partial t}\]

The integrand for the test function gradient term w is

   \[k \frac{\partial T}{\partial x_j} \]

The gradient \(\frac{\partial F}{\partial T} \)
for g3 - integrand for the test function gradient and basis function gradient term is \(k \) for the diagonal terms. All others are zero.

The gradient \(\frac{\partial F}{\partial \partial T/\partial t} \) for g0 - integrand for the test and basis function term is \(\rho C \).  This uses the u\_tShift scalar in petsc.  All others are zero.
*/


/*
 * TODO:
// * create a box mesh so I can the dimensions
 * move the flags into hard code
 * figure out how to change the boundary condition
 */
typedef enum {
    specificHeat, conductivity, density, total
} ConductionProperties;

std::vector<PetscInt> bcNodes;

PetscInt coupledWallId;


/**
 * Compute the test function integrated.  Note there is only a single field.
 */
static void
wIntegrandTestFunction(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[],
                       const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
                       const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
                       const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscInt numConstants,
                       const PetscScalar constants[], PetscScalar f0[]) {
    f0[0] = constants[density] * constants[specificHeat] * u_t[0];
}

/**
 * Compute the test function integrated.  Note there is only a single field.
 */
static void wIntegrandTestGradientFunction(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[],
                                           const PetscInt uOff_x[],
                                           const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                                           const PetscInt aOff[],
                                           const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
                                           const PetscScalar a_x[], PetscReal t, const PetscReal x[],
                                           PetscInt numConstants,
                                           const PetscScalar constants[], PetscScalar f1[]) {
    for (PetscInt d = 0; d < dim; ++d) {
        f1[d] = constants[conductivity] * u_x[d];
    }
}

/**
 * Compute the jacobian term g0 - integrand for the test and basis function term i
 */
static void jacobianG0Term(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[],
                           const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[],
                           const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
                           const PetscScalar a_x[],
                           PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants,
                           const PetscScalar constants[], PetscScalar g3[]) {
    for (PetscInt d = 0; d < dim; ++d) { g3[d * dim + d] = u_tShift * constants[density] * constants[specificHeat]; }
}


/**
 * Compute the jacobian term g3 - integrand for the test function gradient and basis function gradient term
 */
static void jacobianG3Term(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[],
                           const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[],
                           const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
                           const PetscScalar a_x[],
                           PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants,
                           const PetscScalar constants[], PetscScalar g3[]) {
    for (PetscInt d = 0; d < dim; ++d) { g3[d * dim + d] = constants[conductivity]; }
}

static PetscErrorCode
EssentialCoupledWallBC(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u,
                       void *ctx) {
    *u = 500.0;
    std::cout << "calling EssentialCoupledWallBC" << std::endl;
    return PETSC_SUCCESS;
}

static void
NaturalCoupledWallBC(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[],
                     const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
                     const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                     PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants,
                     const PetscScalar constants[], PetscScalar f0[]) {
    // The normal is facing out, so scale the heat flux by -1
    f0[0] = -10000.0;
    std::cout << "calling NaturalCoupledWallBC" << std::endl;

}

static PetscErrorCode
EssentialFarFieldWallBC(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u,
                        void *ctx) {
    *u = 300.0;
    return PETSC_SUCCESS;
}


PetscErrorCode UpdateBoundaryCondition(TS ts) {
    PetscFunctionBeginUser;
    DM dm;
    PetscCall(TSGetDM(ts, &dm));

    // Get the solution vector to see if it should be natural or essential
    DMLabel label;
    PetscCall(DMGetLabel(dm, "marker", &label));

    // Get the current global vector
    Vec currentGlobalVec;
    PetscCall(TSGetSolution(ts, &currentGlobalVec));

    // Get the current time
    PetscReal time;
    PetscCall(TSGetTime(ts, &time));

    // Get the local vector and fill in any boundary values
    Vec locVec;
    PetscCall(DMGetLocalVector(dm, &locVec));
    PetscCall(DMPlexInsertBoundaryValues(dm, PETSC_TRUE, locVec, time, nullptr, nullptr, nullptr));
    PetscCall(DMGlobalToLocal(dm, currentGlobalVec, INSERT_VALUES, locVec));

    // Get the array
    const PetscScalar *locVecArray;
    PetscCall(VecGetArrayRead(locVec, &locVecArray));

    // Get the array
    const PetscScalar *temperaturePointer;
    PetscCall(DMPlexPointLocalFieldRead(dm, bcNodes.front(), 0, locVecArray, &temperaturePointer));
    const PetscScalar temperature = temperaturePointer[0];

    // cleanup
    PetscCall(VecRestoreArrayRead(locVec, &locVecArray));

    VecView(locVec, PETSC_VIEWER_STDOUT_WORLD);

    // Now determine what kind of boundary we need
    DMBoundaryConditionType neededBcType = DM_BC_NATURAL;
    if (temperature > 1000) {
        neededBcType = DM_BC_ESSENTIAL;
    }

    // Assume that the bcId is 0 for now.
    PetscDS ds;
    PetscCall(DMGetDS(dm, &ds));
    DMBoundaryConditionType currentBcType;

    // Get the current bc
    PetscCall(
            PetscDSGetBoundary(ds, coupledWallId, NULL, &currentBcType, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                               NULL, NULL));
    // Change the boundary if needed
    std::cout << "temperature " << temperature << std::endl;

    if (currentBcType != neededBcType) {
        const PetscInt leftWallId = 1;
        switch (neededBcType) {
            case DM_BC_ESSENTIAL:
                std::cout << "switching bc to DM_BC_ESSENTIAL" << std::endl;
                PetscCall(
                        PetscDSUpdateBoundary(ds, coupledWallId, DM_BC_ESSENTIAL, "coupledWall", label, 1, &leftWallId, 0, 0, NULL,
                                              (void (*)(void)) EssentialCoupledWallBC, NULL,
                                              NULL));
                break;
            case DM_BC_NATURAL:
                std::cout << "switching bc to DM_BC_NATURAL" << std::endl;
                PetscCall(
                        PetscDSUpdateBoundary(ds, coupledWallId, DM_BC_NATURAL, "coupledWall", label, 1, &leftWallId, 0, 0, NULL,
                                              NULL, NULL,
                                              NULL));
                break;
            default:
                throw std::invalid_argument("Unable to handle BC type");
        }

        // Reset the TS
        PetscCall(TSReset(ts));

        // Create a new global vector
        Vec newGlobalVector;
        PetscCall(DMCreateGlobalVector(dm, &newGlobalVector));

        // Map from the local vector back to the global
        PetscCall(DMLocalToGlobal(dm, locVec, INSERT_VALUES, newGlobalVector));

        // Set in the TS
        PetscCall(TSSetSolution(ts, newGlobalVector));
    }


    // Cleanup
    PetscCall(DMRestoreLocalVector(dm, &locVec));

    PetscFunctionReturn(PETSC_SUCCESS);
}

//PetscErrorCode Ch(TS, PetscInt, PetscReal, Vec, void *)



static PetscErrorCode CreateMesh(MPI_Comm comm, PetscOptions options, DM *dm) {
    PetscFunctionBeginUser;

    PetscCall(DMCreate(comm, dm));
    PetscCall(DMSetType(*dm, DMPLEX));
    PetscCall(PetscObjectSetOptions((PetscObject) *dm, options));
    PetscCall(PetscObjectSetName((PetscObject) *dm, "oneDimMesh"));
    PetscCall(DMSetFromOptions(*dm));
    PetscCall(DMViewFromOptions(*dm, NULL, "-dm_view"));

    PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SetupProblem(DM dm) {
    PetscDS ds;
    DMLabel label;

    PetscFunctionBeginUser;
    PetscCall(DMGetLabel(dm, "marker", &label));
    PetscCall(DMPlexLabelComplete(dm, label));
    PetscCall(DMGetDS(dm, &ds));
    PetscCall(PetscDSSetJacobian(ds, 0, 0, jacobianG0Term, NULL, NULL, jacobianG3Term));
    PetscCall(PetscDSSetResidual(ds, 0, wIntegrandTestFunction, wIntegrandTestGradientFunction));

    // hard code a boundary condition for now
    const PetscInt leftWallId = 1;
//    PetscCall(DMAddBoundary(dm, DM_BC_NATURAL, "coupledWall", label, 1, &leftWallId, 0, 0, NULL,
//                            NULL, NULL,
//                            NULL, &coupledWallId));
//    PetscWeakForm wf;
//    PetscCall(PetscDSGetBoundary(ds, coupledWallId, &wf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
//                                 NULL));
//    PetscCall(PetscWeakFormSetIndexBdResidual(wf, label, leftWallId, 0, 0, 0, NaturalCoupledWallBC, 0, NULL));

    printf("CoupledBoundaryWall: %" PetscInt_FMT, coupledWallId);
    // hard code a boundary condition for now
    const PetscInt rightWallId = 2;
    PetscCall(PetscDSAddBoundary(ds, DM_BC_ESSENTIAL, "farFieldWall", label, 1, &rightWallId, 0, 0, NULL,
                                 (void (*)(void)) EssentialFarFieldWallBC, NULL,
                                 NULL, NULL));

    // Set the constant values for this problem
    PetscReal parameterArray[total] = {1000, 1, 1};
    PetscCall(PetscDSSetConstants(ds, total, parameterArray));


    // Start out getting all the faces
    PetscInt depth;
    PetscCall(DMPlexGetDepth(dm, &depth));

    // Get the label of points in this
    PetscInt pStart, pEnd;
    PetscCall(DMPlexGetChart(dm, &pStart, &pEnd));

    // get the global section
    PetscSection section;
    PetscCall(DMGetSection(dm, &section));


    // Determine the BC Node
    for (PetscInt p = pStart; p < pEnd; ++p) {
        // Get the dof here
        PetscInt dof;
        PetscCall(PetscSectionGetDof(section, p, &dof));
        // Get the label here
        PetscInt bcValue;
        PetscCall(DMLabelGetValue(label, p, &bcValue));

        if (dof && bcValue == 1) {
            bcNodes.push_back(p);
        }
    }

    if (bcNodes.size() != 1) {
        throw std::invalid_argument("There should be a single bc node");
    }


    PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SetupDiscretization(DM dm) {
    DM cdm = dm;
    PetscFE fe;
    PetscInt dim;

    PetscFunctionBeginUser;
    PetscCall(DMGetDimension(dm, &dim));
    /* Create finite element */
    PetscCall(
            PetscFECreateLagrange(PETSC_COMM_SELF, dim, 1, PETSC_TRUE, 1 /*degree  of space */ , PETSC_DETERMINE, &fe));
    PetscCall(PetscObjectSetName((PetscObject) fe, "temperature"));
    /* Set discretization and boundary conditions for each mesh */
    PetscCall(DMSetField(dm, 0, NULL, (PetscObject) fe));
    PetscCall(DMCreateDS(dm));

    // Setup the problem
    PetscCall(SetupProblem(dm));

    // copy over the discratation
    while (cdm) {
        PetscCall(DMCopyDisc(dm, cdm));
        PetscCall(DMGetCoarseDM(cdm, &cdm));
    }
    PetscCall(PetscFEDestroy(&fe));
    PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SetInitialConditions(TS ts, Vec u) {
    DM dm;
    PetscReal t;

    PetscFunctionBeginUser;
    PetscCall(TSGetDM(ts, &dm));
    PetscCall(TSGetTime(ts, &t));

    // Set the initial condition
    VecSet(u, 300.0);

    PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc, char **argv) {
    DM dm;
    TS ts;
    Vec u;
    PetscOptions options;
    PetscFunctionBeginUser;
    PetscCall(PetscInitialize(&argc, &argv, NULL, help));

    // Create a petsc options and set the required options
    PetscCall(PetscOptionsCreate(&options));

    // for now, hard code the petsc options
    PetscOptionsSetValue(options, "-ts_type", "beuler");
    PetscOptionsSetValue(options, "-ts_max_steps", "500");
    PetscOptionsSetValue(options, "-ts_dt", "0.1");
    PetscOptionsSetValue(options, "-snes_error_if_not_converged", NULL);
    PetscOptionsSetValue(options, "-dm_view", "hdf5:/Users/mcgurn/scratch/results/testOutput/sol.h5");
    PetscOptionsSetValue(options, "-ts_monitor_solution",
                         "hdf5:/Users/mcgurn/scratch/results/testOutput/sol.h5::append");
    PetscOptionsSetValue(options, "-pc_type", "lu");
    // Set the mesh parameters
    PetscOptionsSetValue(options, "-dm_plex_separate_marker", NULL);
    PetscOptionsSetValue(options, "-dm_plex_dim", "1");
    PetscOptionsSetValue(options, "-dm_plex_box_faces", "10");
    PetscOptionsSetValue(options, "-dm_plex_box_upper", "0.1");

    PetscCall(CreateMesh(PETSC_COMM_WORLD, options, &dm));
    PetscCall(SetupDiscretization(dm));

    PetscCall(TSCreate(PETSC_COMM_WORLD, &ts));
    PetscCall(PetscObjectSetOptions((PetscObject) ts, options));
    PetscCall(TSSetDM(ts, dm));
    PetscCall(DMTSSetBoundaryLocal(dm, DMPlexTSComputeBoundary, NULL));
    PetscCall(DMTSSetIFunctionLocal(dm, DMPlexTSComputeIFunctionFEM, NULL));
    PetscCall(DMTSSetIJacobianLocal(dm, DMPlexTSComputeIJacobianFEM, NULL));
    PetscCall(TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP));
    PetscCall(TSSetFromOptions(ts));
    PetscCall(TSSetComputeInitialCondition(ts, SetInitialConditions));
    PetscCall(TSSetPreStep(ts, UpdateBoundaryCondition));
    PetscCall(DMCreateGlobalVector(dm, &u));
    PetscCall(DMTSCheckFromOptions(ts, u));
    PetscCall(SetInitialConditions(ts, u));
    PetscCall(PetscObjectSetName((PetscObject) u, "temperature"));
    PetscCall(TSSolve(ts, u));

    PetscCall(VecDestroy(&u));
    PetscCall(TSDestroy(&ts));
    PetscCall(DMDestroy(&dm));

    PetscCall(PetscOptionsLeft(options));
    PetscCall(PetscOptionsDestroy(&options));

    PetscObjectsDump(PETSC_STDOUT, PETSC_TRUE);
    PetscCall(PetscFinalize());
    return 0;
}

